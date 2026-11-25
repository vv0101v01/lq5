#!/usr/bin/env python3
"""reporthtml.py — one description of the report chrome, shared by all six
`make_report_*.py` generators.

Before this module, `load()`, `table()` and the SVG chart primitives were
copy-pasted into five files (96–99 % identical by diff, ~800 duplicated lines):
a bug fixed in one chart never reached the others, and every new round grew the
debt.  Everything here is pure stdlib (csv/math/html), so the repo still needs
nothing but Python 3 to regenerate its reports.

Public API
----------
    load(path)                     -> {column: [float | str]}
    trace_path(name), doc_path(name) -> absolute paths under the repo
    csv_rows(name)                 -> list of row dicts straight from a trace CSV
    table(headers, rows, caption)  -> <table>…</table>
    line_chart(series, …)          -> inline <svg>…</svg>
    bar_chart(labels, values, …)   -> inline <svg>…</svg>
    PALETTE                        -> shared colour names

`series` is [(label, values, colour)]; charts are size-agnostic (viewBox +
width:100%), so they render inside sandboxed previews with no external assets.
"""
from __future__ import annotations

import csv
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Toolchain stamp printed in the report badges.  `verify.sh` refreshes these
# from the toolchain it actually ran, so a report can never claim a version or
# a test count that the gate did not produce; the defaults only apply when a
# single generator is run by hand.
MOON_VER = os.environ.get("MOON_VER", "0.1.20261119")
MOONC_VER = os.environ.get("MOONC_VER", "v0.10.22+124d61d05")
MOON_TESTS = os.environ.get("MOON_TESTS", "135")

PALETTE = {
    "ink": "#0f172a",
    "muted": "#64748b",
    "grid": "#e2e8f0",
    "blue": "#2563eb",
    "red": "#dc2626",
    "green": "#059669",
    "amber": "#d97706",
    "violet": "#7c3aed",
    "teal": "#0d9488",
    "rose": "#e11d48",
}


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load(path):
    """CSV → {column: [float]} with a str fallback for non-numeric columns."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    cols = {}
    if not rows:
        return cols
    for key in rows[0]:
        try:
            cols[key] = [float(r[key]) for r in rows]
        except ValueError:
            cols[key] = [r[key] for r in rows]
    return cols


def trace_path(name):
    """`trace_path("emergence_trace.csv")` → absolute path next to moon.mod."""
    return os.path.join(ROOT, name)


def doc_path(name):
    """`doc_path("report.html")` → absolute path inside the committed docs/ dir."""
    return os.path.join(ROOT, "docs", name)


def csv_rows(name):
    """CSV → list of dicts, for generators that need raw column names per row."""
    with open(trace_path(name), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
# charts
# ---------------------------------------------------------------------------
def _frame(width, height, title):
    return [
        f'<svg viewBox="0 0 {width} {height}" width="100%" role="img"'
        + (f' aria-label="{title}"' if title else "")
        + f' style="background:#fff;border:1px solid {PALETTE["grid"]};border-radius:10px">'
    ]


def line_chart(series, width=760, height=250, title="", ylabel="", log=False):
    """Multi-series line chart with a 5-line grid and an inline legend."""
    pad_l, pad_r, pad_t, pad_b = 58, 14, 32, 34
    iw, ih = width - pad_l - pad_r, height - pad_t - pad_b
    n = max(len(vs) for _, vs, _ in series)
    vals = [v for _, vs, _ in series for v in vs]
    if log:
        vals = [math.log10(max(abs(v), 1e-18)) for v in vals]
    y0, y1 = min(vals), max(vals)
    if y1 - y0 < 1e-12:
        y1 = y0 + 1.0
    span = y1 - y0

    def X(i):
        return pad_l + iw * i / max(n - 1, 1)

    def Y(v):
        if log:
            v = math.log10(max(abs(v), 1e-18))
        return pad_t + ih - ih * (v - y0) / span

    parts = _frame(width, height, title)
    parts.append(
        f'<text x="{pad_l}" y="20" font-size="13" font-weight="600" '
        f'fill="{PALETTE["ink"]}">{title}</text>'
    )
    for k in range(5):
        yy = pad_t + ih * k / 4
        parts.append(
            f'<line x1="{pad_l}" y1="{yy:.1f}" x2="{pad_l+iw}" y2="{yy:.1f}" '
            f'stroke="{PALETTE["grid"]}" stroke-width="1"/>'
        )
        lab = y1 - span * k / 4
        text = f"1e{lab:.0f}" if log else f"{lab:.3g}"
        parts.append(
            f'<text x="{pad_l-8}" y="{yy+4:.1f}" font-size="10" text-anchor="end" '
            f'fill="{PALETTE["muted"]}">{text}</text>'
        )
    if ylabel:
        parts.append(
            f'<text x="14" y="{pad_t+ih/2:.0f}" font-size="10" fill="{PALETTE["muted"]}" '
            f'transform="rotate(-90 14 {pad_t+ih/2:.0f})" text-anchor="middle">{ylabel}</text>'
        )
    for label, vs, colour in series:
        pts = " ".join(f"{X(i):.1f},{Y(v):.2f}" for i, v in enumerate(vs))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{colour}" stroke-width="2.1"/>')
    lx = pad_l + 10
    for label, _, colour in series:
        parts.append(f'<rect x="{lx}" y="{pad_t+8}" width="14" height="4" fill="{colour}"/>')
        parts.append(
            f'<text x="{lx+20}" y="{pad_t+13}" font-size="11" '
            f'fill="{PALETTE["ink"]}">{label}</text>'
        )
        lx += 30 + 7 * len(label)
    parts.append(
        f'<line x1="{pad_l}" y1="{pad_t+ih}" x2="{pad_l+iw}" y2="{pad_t+ih}" '
        f'stroke="{PALETTE["muted"]}"/>'
    )
    for i in range(0, n, max(n // 8, 1)):
        parts.append(
            f'<text x="{X(i):.0f}" y="{pad_t+ih+16}" font-size="10" text-anchor="middle" '
            f'fill="{PALETTE["muted"]}">{i}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def bar_chart(labels, values, width=760, height=230, title="", colour=PALETTE["blue"], ylabel=""):
    """Single-series bar chart with value labels on top of each bar."""
    pad_l, pad_r, pad_t, pad_b = 58, 14, 30, 44
    iw, ih = width - pad_l - pad_r, height - pad_t - pad_b
    y1 = max(values) * 1.2 if values and max(values) > 0 else 1.0
    n = max(len(values), 1)
    bw = iw / n * 0.58
    parts = _frame(width, height, title)
    parts.append(
        f'<text x="{pad_l}" y="18" font-size="13" font-weight="600" '
        f'fill="{PALETTE["ink"]}">{title}</text>'
    )
    for k in range(5):
        yy = pad_t + ih * k / 4
        parts.append(
            f'<line x1="{pad_l}" y1="{yy:.1f}" x2="{pad_l+iw}" y2="{yy:.1f}" '
            f'stroke="{PALETTE["grid"]}"/>'
        )
        parts.append(
            f'<text x="{pad_l-8}" y="{yy+4:.1f}" font-size="10" text-anchor="end" '
            f'fill="{PALETTE["muted"]}">{y1*(1-k/4):.3g}</text>'
        )
    if ylabel:
        parts.append(
            f'<text x="14" y="{pad_t+ih/2:.0f}" font-size="10" fill="{PALETTE["muted"]}" '
            f'transform="rotate(-90 14 {pad_t+ih/2:.0f})" text-anchor="middle">{ylabel}</text>'
        )
    for i, (lab, v) in enumerate(zip(labels, values)):
        x = pad_l + iw * (i + 0.21) / n
        h = ih * max(v, 0.0) / y1
        parts.append(
            f'<rect x="{x:.1f}" y="{pad_t+ih-h:.1f}" width="{bw:.1f}" height="{h:.1f}" '
            f'fill="{colour}" rx="4"/>'
        )
        parts.append(
            f'<text x="{x+bw/2:.1f}" y="{pad_t+ih+15}" font-size="10" '
            f'text-anchor="middle" fill="{PALETTE["muted"]}">{lab}</text>'
        )
        parts.append(
            f'<text x="{x+bw/2:.1f}" y="{pad_t+ih-h-5:.1f}" font-size="10" '
            f'text-anchor="middle" fill="{PALETTE["ink"]}">{v:g}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------
def table(headers, rows, caption=""):
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    cap = f'<div class="cap">{caption}</div>' if caption else ""
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>{cap}"


# ---------------------------------------------------------------------------
# self-test:  python3 scripts/reporthtml.py
# ---------------------------------------------------------------------------
def _selftest() -> None:
    import os
    import tempfile

    d = tempfile.mkdtemp()
    p = os.path.join(d, "t.csv")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write("step,x,label\n0,1.5,a\n1,2.5,b\n")
    cols = load(p)
    assert cols["step"] == [0.0, 1.0] and cols["label"] == ["a", "b"], cols

    svg = line_chart([("u", [1, 2, 3], PALETTE["blue"]), ("v", [3, 2, 1], PALETTE["red"])],
                     title="t", ylabel="y")
    assert svg.startswith("<svg") and svg.endswith("</svg>")
    assert svg.count("<polyline") == 2 and "<text" in svg
    flat = line_chart([("u", [1, 1, 1], PALETTE["blue"])], title="flat")
    assert flat.count("<polyline") == 1           # zero-span series must not divide by 0
    lg = line_chart([("u", [1e-30, 1e-6, 1], PALETTE["teal"])], log=True, title="log")
    assert "1e" in lg

    bar = bar_chart(["a", "b"], [1.0, 2.0], title="b")
    assert bar.count("<rect") == 2, bar.count("<rect")   # one per bar
    assert bar_chart([], [], title="empty").startswith("<svg")   # empty must not crash

    tb = table(["h1", "h2"], [["1", "2"]], caption="cap")
    assert "<th>h1</th>" in tb and "<td>2</td>" in tb and "cap" in tb
    assert "<td>" not in table(["h"], [])
    print("reporthtml selftest OK · load · line/bar degenerate cases · table")


if __name__ == "__main__":
    _selftest()
