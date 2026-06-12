#!/usr/bin/env bash
# Run the analysis smoke suite (setupfiles/analysis_smoketests/test_analyses.py)
# headless, from anywhere.
#
# This is the CI entry point: it cd's to the testsuite root, forces a non-
# interactive matplotlib backend, puts setupfiles/ on PYTHONPATH (so the
# IC_analysis_* modules import), and runs the suite with unittest. Exits non-zero
# if any test fails -- so a CI job is just
# `bash setupfiles/analysis_smoketests/run.sh`.
#
# Usage:
#   bash setupfiles/analysis_smoketests/run.sh   # whole suite
#   bash setupfiles/analysis_smoketests/run.sh test_analyses.TestSedov
#       # a single case (passed through; run dir = analysis_smoketests)
set -euo pipefail

# testsuite root = two levels up from this script's directory
# (testsuite/setupfiles/analysis_smoketests/run.sh).
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$(dirname "$here")")"
cd "$root"

export MPLBACKEND=Agg
export PYTHONPATH="$root/setupfiles${PYTHONPATH:+:$PYTHONPATH}"

if [ "$#" -gt 0 ]; then
    cd "$here"
    exec python -m unittest "$@"
fi
exec python -m unittest discover -s setupfiles/analysis_smoketests -p "test_*.py" -v
