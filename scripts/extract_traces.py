#!/usr/bin/env python3
"""Split `run_output.txt` into the CSV traces used by the reports.

Both pipelines print their data between markers:

    # BEGIN_CSV ... # END_CSV                (round 1: supergeometry)
    # BEGIN_CSV_LANG ... # END_CSV_LANG      (round 2: Langlands)
    # ... through # BEGIN_CSV_CORRESPONDENCE  (round 6: correspondence)
    # ... and   # BEGIN_CSV_OCTAVE            (round 7: ten channels)

This script turns those blocks into `emergence_trace.csv` and
`langlands_trace.csv`, so the CSVs are always a *derived* artefact of a real
run instead of a hand-edited file.

Usage: python3 scripts/extract_traces.py
"""
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "run_output.txt")

BLOCKS = [
    ("# BEGIN_CSV", "# END_CSV", "emergence_trace.csv"),
    ("# BEGIN_CSV_LANG", "# END_CSV_LANG", "langlands_trace.csv"),
    ("# BEGIN_CSV_FRONTIER", "# END_CSV_FRONTIER", "frontier_trace.csv"),
    ("# BEGIN_CSV_UNIVERSE", "# END_CSV_UNIVERSE", "universe_trace.csv"),
    ("# BEGIN_CSV_HEXAD", "# END_CSV_HEXAD", "hexad_trace.csv"),
    (
        "# BEGIN_CSV_CORRESPONDENCE",
        "# END_CSV_CORRESPONDENCE",
        "correspondence_trace.csv",
    ),
    ("# BEGIN_CSV_OCTAVE", "# END_CSV_OCTAVE", "octave_trace.csv"),
]


def main():
    text = open(OUT, encoding="utf-8").read()
    for begin, end, name in BLOCKS:
        m = re.search(re.escape(begin) + r"\n(.*?)\n" + re.escape(end), text, re.S)
        if not m:
            raise SystemExit(f"marker {begin} not found in {OUT}")
        rows = [ln for ln in m.group(1).strip().split("\n") if ln.strip()]
        path = os.path.join(ROOT, name)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(rows) + "\n")
        print(f"wrote {path} ({len(rows) - 1} data rows)")


if __name__ == "__main__":
    main()
