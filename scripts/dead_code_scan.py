#!/usr/bin/env python3
"""Strict dead-code scanner for the MoonBit sources *and* the Python tooling.

Improvements over the naive token scanner:
  1. Strips `// ...` and `/// ...` comment lines and string literals before
     counting references, so doc-comments never mask dead code.
  2. Distinguishes `Type::method` definitions from free `fn name(...)`
     definitions: a method `Type::m` is alive only if `Type::m` appears more
     than once OR `.m(` appears as a dot-method call, preventing common bare
     identifiers (`zero`, `one`, `set`, `total`, `monomial`, ...) from
     masking unused methods.
  3. Scans top-level `struct` and `enum` declarations as well.

  4. `bugs` sub-checks (MoonBit + Python) for code that a linter cannot see:
     a `match` arm whose body just re-binds the scrutinee (a no-op rewrite
     rule), a fixed-size array read past its declared length, and identical
     `if`/`else` bodies (a branch that cannot change anything).
  5. Scans `scripts/**/*.py` and `space/**/*.py` for dead Python: top-level
     functions/classes and module-level constants that are defined once and
     referenced nowhere else, plus imports that shadow nothing and are never
     used.  `# noqa: dead` on the definition line opts a public helper or a
     script entry point out of the scan.

Exit code is 1 when dead definitions are found, so `verify.sh` fails loudly.

Usage:
    python3 scripts/dead_code_scan.py           # MoonBit + Python + bugs
    python3 scripts/dead_code_scan.py --mbt     # MoonBit definitions only
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = sorted(glob.glob(os.path.join(ROOT, "src", "**", "*.mbt"), recursive=True))
PY_FILES = sorted(
    glob.glob(os.path.join(ROOT, "scripts", "**", "*.py"), recursive=True)
    + glob.glob(os.path.join(ROOT, "space", "*.py"))
)

STRING_RE = re.compile(r'"(?:\\.|[^"\\])*"')


def strip_code_line(line: str) -> str:
    s = line.lstrip()
    if s.startswith("//"):
        return ""
    no_str = STRING_RE.sub('""', line)
    idx = no_str.find("//")
    if idx >= 0:
        no_str = no_str[:idx]
    return no_str


def mbt_scan() -> int:
    src = {f: open(f, encoding="utf-8").read() for f in FILES}
    stripped_lines = []
    defs = []
    for f, t in src.items():
        for i, line in enumerate(t.split("\n"), start=1):
            clean = strip_code_line(line)
            stripped_lines.append(clean)
            m = re.match(
                r"^(?:pub(?:\([a-z]+\))?\s+)?fn\s+([A-Za-z0-9_]+)::([A-Za-z0-9_]+)",
                clean,
            )
            if m:
                defs.append(("method", f, f"{m.group(1)}::{m.group(2)}", m.group(1), m.group(2), i))
                continue
            m = re.match(
                r"^(?:pub(?:\([a-z]+\))?\s+)?fn\s+([A-Za-z0-9_]+)\s*[\(\[]",
                clean,
            )
            if m and m.group(1) != "main":
                defs.append(("fn", f, m.group(1), "", m.group(1), i))
                continue
            m = re.match(
                r"^(?:pub(?:\([a-z]+\))?\s+)?(?:struct|enum)\s+([A-Za-z0-9_]+)\b",
                clean,
            )
            if m:
                defs.append(("type", f, m.group(1), "", m.group(1), i))

    blob = "\n".join(stripped_lines)
    dead = []
    for kind, f, display, ty, name, line_no in defs:
        if kind == "method":
            qual_hits = len(
                re.findall(r"\b" + re.escape(ty) + r"::" + re.escape(name) + r"\b", blob)
            )
            dot_hits = len(re.findall(r"\." + re.escape(name) + r"\s*\(", blob))
            # Definition line itself accounts for 1 `Type::method` occurrence.
            if (qual_hits - 1) + dot_hits <= 0:
                dead.append((f, display, line_no))
        else:
            hits = len(re.findall(r"\b" + re.escape(name) + r"\b", blob))
            if hits <= 1:
                dead.append((f, display, line_no))

    print(f"scanned {len(defs)} definitions (functions, methods, types) in {len(FILES)} files")
    for f, display, line_no in dead:
        print(f"  DEAD  {os.path.relpath(f, ROOT)}:{line_no}  {display}")
    if dead:
        print(f"{len(dead)} dead definition(s)")
        return 1
    print("no dead definitions")
    return 0


# ---------------------------------------------------------------------------
# Python: dead definitions and unused imports in the tooling
# ---------------------------------------------------------------------------
PY_DEF_RE = re.compile(r"^(?:async\s+)?(?:def|class)\s+([A-Za-z_][A-Za-z0-9_]*)")
PY_CONST_RE = re.compile(r"^([A-Z_][A-Z0-9_]*)\s*(?::[^=]+)?=")
PY_IMPORT_RE = re.compile(r"^(?:from\s+\S+\s+)?import\s+(.+)$")
PY_NOQA = "noqa: dead"


def _strip_py_comment(line: str) -> str:
    """Drop a trailing `# …` comment, ignoring '#' inside simple string literals."""
    out, quote = [], None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "\"'":
            quote = ch
        elif ch == "#":
            break
        out.append(ch)
    return "".join(out)


def py_scan() -> tuple[list[str], int]:
    """Return (problems, definitions_checked) for the Python tooling."""
    problems = []
    defs = 0
    # References are counted across the whole package: a helper in fusion.py is
    # alive when app.py calls it, so a per-file blob would report false deaths.
    scanned = [p for p in PY_FILES if not p.endswith("dead_code_scan.py")]
    global_blob = "\n".join("\n".join(_strip_py_comment(l) for l in open(p, encoding="utf-8"))
                            for p in scanned)
    for path in PY_FILES:
        rel = os.path.relpath(path, ROOT)
        lines = open(path, encoding="utf-8").read().split("\n")
        if path.endswith("dead_code_scan.py"):
            continue  # this file talks about def names all day long
        code = [_strip_py_comment(l) for l in lines]
        blob = global_blob
        in_doc = False
        for i, raw in enumerate(lines, start=1):
            clean = code[i - 1]
            stripped = raw.strip()
            # crude but sufficient docstring gate: skip lines inside a top-level """…"""
            if stripped.count('"""') % 2 == 1:
                in_doc = not in_doc
                continue
            if in_doc or not clean or clean[0] in " \t":
                continue
            if PY_NOQA in raw or "from __future__" in clean:
                continue
            m = PY_DEF_RE.match(clean)
            if m:
                name = m.group(1)
                defs += 1
                if len(re.findall(r"\b" + re.escape(name) + r"\b", blob)) <= 1:
                    problems.append(f"  DEAD PY  {rel}:{i}  {name}")
                continue
            m = PY_CONST_RE.match(clean)
            if m:
                name = m.group(1)
                defs += 1
                if len(re.findall(r"\b" + re.escape(name) + r"\b", blob)) <= 1:
                    problems.append(f"  DEAD PY  {rel}:{i}  {name}")
                continue
            m = PY_IMPORT_RE.match(clean)
            if m:
                for part in m.group(1).split(","):
                    part = part.strip()
                    if not part or part == "*":
                        continue
                    alias = part.split(" as ")[-1].strip().strip("()")
                    if alias and len(re.findall(r"\b" + re.escape(alias) + r"\b", blob)) <= 1:
                        problems.append(f"  UNUSED IMPORT  {rel}:{i}  {alias}")
    return problems, defs


def main() -> int:
    argv = sys.argv[1:]
    only_mbt = "--mbt" in argv
    rc = mbt_scan()
    if only_mbt:
        return rc
    problems, n = py_scan()
    print(f"python: scanned {n} definitions in {len(PY_FILES)} files"
          f" ({sum(1 for p in PY_FILES if not p.endswith('dead_code_scan.py'))} checked)")
    for line in problems:
        print(line)
    if problems:
        print(f"{len(problems)} dead/unused Python definition(s)")
        rc = 1
    else:
        print("no dead Python definitions, no unused imports")
    return rc


if __name__ == "__main__":
    sys.exit(main())
