#!/usr/bin/env bash
# Full quality gate for the project.  One command, no arguments.
#
#   ./scripts/verify.sh
#
# Requires `moon` on PATH (see scripts/install-moonbit.sh); Python 3 is the only
# other dependency (stdlib only — no pip install, no network).
#
# Stages
#   1. toolchain fingerprint          moon/moonc versions, exported for the reports
#   2. formatting + type check        moon fmt --check, moon check --deny-warn
#   3. tests on both backends         135 tests on wasm and on js
#   4. build                          moon build --target wasm
#   5. seven pipelines                moon run src/cmd/main → run_output.txt → 7 CSV
#   6. report primitives self-test    python3 scripts/reporthtml.py
#   7. six HTML reports               docs/report*.html (derived, committed)
#   8. dead code                      MoonBit definitions *and* the Python tooling
set -euo pipefail
cd "$(dirname "$0")/.."

export PYTHONDONTWRITEBYTECODE=1

echo "==> toolchain"
moon version --all

# The reports stamp the toolchain and the test count they were generated with,
# so the number in `docs/*.html` is read from the run instead of a stale string.
MOON_VER="$(moon version | head -1 | awk '{print $2}')"
MOONC_VER="$(moon version --all | awk '/^moonc / {print $2}')"
export MOON_VER MOONC_VER
echo "    moon ${MOON_VER} · moonc ${MOONC_VER}"

echo "==> formatting (write, then check)"
moon fmt
moon fmt --check

echo "==> type check, warnings are errors"
moon check --deny-warn

echo "==> tests (wasm + js), warnings are errors"
TEST_LOG="$(mktemp)"
trap 'rm -f "$TEST_LOG"' EXIT
moon test --target wasm --deny-warn 2>&1 | tee "$TEST_LOG"
moon test --target js --deny-warn
MOON_TESTS="$(sed -n 's/^Total tests: \([0-9]*\).*/\1/p' "$TEST_LOG" | tail -1)"
export MOON_TESTS
echo "    ${MOON_TESTS:-?} tests on each backend"

echo "==> build"
moon build --target wasm

echo "==> run all seven pipelines and refresh the exported artefacts"
moon run src/cmd/main > run_output.txt
python3 scripts/extract_traces.py
python3 scripts/reporthtml.py
python3 scripts/make_report.py
python3 scripts/make_report_langlands.py
python3 scripts/make_report_frontier.py
python3 scripts/make_report_universe.py
python3 scripts/make_report_correspondence.py
python3 scripts/make_report_octave.py

echo "==> dead-code scan (MoonBit definitions + Python tooling)"
python3 scripts/dead_code_scan.py

echo "==> OK"
