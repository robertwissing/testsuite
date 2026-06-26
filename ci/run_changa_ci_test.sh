#!/usr/bin/env bash
# Build the ChaNGa fork in CI config, then run testsuite test(s) at +p4 -- either
# ALL tests in ci/matrix.txt (when no <test> is given) or ONE explicit test --
# and GATE each against its committed reference, or (--bless) generate those
# references. Blessing just drives `runtest.sh -r`, which writes the reference to
# test_cases/<test>/references/<runname>/ automatically (runtest.sh cd's into
# test_cases/<test>/ and -r saves relative to there). --bless-readydata does the
# same blessing but on an EXISTING run dir -- it skips the build and the sim and
# only runs the analysis --save-reference (use it after a run, or to re-bless
# from data you already have).
#
# LOCAL helper (the CI workflow is self-contained/inline; it does NOT call this).
# Use it to build the CI config and bless references before pushing, or to
# reproduce a gate locally. ci/matrix.txt lists the same tests as the workflow.
#
# Optional `-a..-h value` flags (after --bless, before <changa-src>) are the
# test's IC/analysis params; they are forwarded UNCHANGED to BOTH runtest.sh (to
# build the IC) and IC_analysis_<test>.py (same -a..-h mapping). E.g. sedov's -a
# selects the case (0 Pin/Rin, 1 point, 2 kernel) for the sim AND the analysis.
#
# Global options (before <changa-src>): --bless | --bless-readydata (mode);
#   --matrix <file> (alternative test list for matrix runs, e.g. a larger nightly
#   ci/matrix_full.txt; default ci/matrix.txt); --extra "<flags>" (single-test
#   runtest.sh [extra], e.g. --extra "-n 50" to cap steps -- matrix rows use the
#   per-row 'extra:' column instead).
#
# Usage:
#   run_changa_ci_test.sh                   <changa-src>                                            # gate ALL ci/matrix.txt tests
#   run_changa_ci_test.sh --bless           <changa-src>                                            # bless ALL ci/matrix.txt tests
#   run_changa_ci_test.sh --bless-readydata <changa-src>                                            # bless ALL from existing run dirs (no sim)
#   run_changa_ci_test.sh --matrix ci/matrix_full.txt <changa-src>                                  # gate ALL rows of an alternative list
#   run_changa_ci_test.sh         [-a..-h v] [--extra "-n 50"] <changa-src> <test> <Nsmooth> <Nres> <distrib> <reg-tol>   # gate one
#   run_changa_ci_test.sh --bless [-a..-h v] [--extra "-n 50"] <changa-src> <test> <Nsmooth> <Nres> <distrib>             # bless one
#   run_changa_ci_test.sh --bless-readydata  [-a..-h v]        <changa-src> <test> <Nsmooth> <Nres> <distrib>             # bless one from existing data (no sim)
# (The build always runs first and is idempotent; there is no separate build-only
#  mode -- a bare <changa-src> builds, then runs the matrix.)
#
# Env: PE=4, CHANGA_CI_EXT=~/.changa_ci_ext
set -uo pipefail

# Mode: gate (default, run sim + compare), bless (run sim + save reference), or
# blessready (data already present -> skip the sim, only analyze + save reference).
# --matrix <file> points the matrix run at an alternative test list (e.g. a larger
# nightly suite, ci/matrix_full.txt), defaulting to ci/matrix.txt.
MODE=gate
MATRIX_FILE=""
EXTRA_CLI=""        # single-test runtest.sh [extra] (e.g. --extra "-n 50"); matrix rows carry their own
while [ $# -gt 0 ]; do
    case "$1" in
        --bless)           MODE=bless;        shift ;;
        --bless-readydata) MODE=blessready;   shift ;;
        --matrix)          MATRIX_FILE="${2:?--matrix needs a file}"; shift 2 ;;
        --extra)           EXTRA_CLI="${2:?--extra needs a value}";   shift 2 ;;
        *) break ;;
    esac
done

# Collect leading test-param flags (-a..-h value) to mirror to runtest.sh + the
# analysis. getopts stops at the first non-option (the <changa-src> positional).
declare -a PARAM_FLAGS=()
while getopts ":a:b:c:d:e:f:g:h:" opt; do
    case "$opt" in
        [a-h]) PARAM_FLAGS+=("-$opt" "$OPTARG") ;;
        :)  echo "FATAL: option -$OPTARG requires a value"; exit 1 ;;
        \?) echo "FATAL: unknown option -$OPTARG"; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

CHANGA_SRC="${1:?usage: run_changa_ci_test.sh [--bless | --bless-readydata] [--matrix <file>] [--extra <flags>] [-a..-h v] <changa-src> [<test> <Nsmooth> <Nres> <distrib> [reg-tol]]}"
TEST="${2:-}"; NSM="${3:-}"; NRES="${4:-}"; DIST="${5:-}"; TOL="${6:-}"
TS_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PE="${PE:-2}"        # ++local +p 2 by default, matching the CI launcher
EXT="${CHANGA_CI_EXT:-$HOME/.changa_ci_ext}"
CHARM_TARGET="netlrts-linux-x86_64-smp"
# ChaNGa branch to build: the MHD development branch for now; switch to the
# public 'master' once changaMHDTP is upstreamed. Override with CHANGA_BRANCH=...
CHANGA_BRANCH="${CHANGA_BRANCH:-MHD_tp_replication}"

# Resolve relative paths against the invocation cwd, since we cd to TS_ROOT below.
[ "${CHANGA_SRC#/}" = "$CHANGA_SRC" ] && CHANGA_SRC="$PWD/$CHANGA_SRC"
[ "${EXT#/}" = "$EXT" ] && EXT="$PWD/$EXT"
[ -d "$CHANGA_SRC" ] || { echo "FATAL: changa-src not found: $CHANGA_SRC"; exit 1; }
# runtest.sh reads these from the environment (no source edits): point CHANGA_DIR
# at the build dir, and use the CI-consistent ++local +p N charmrun launcher so a
# locally blessed reference matches what CI produces.
export CHANGA_DIR="$CHANGA_SRC"
export CHANGA_CHARMRUN_OPTS="++local +p $PE"
# A user-given --matrix is relative to the invocation cwd; the default lives under
# TS_ROOT (read after the cd below).
if [ -n "$MATRIX_FILE" ]; then
    [ "${MATRIX_FILE#/}" = "$MATRIX_FILE" ] && MATRIX_FILE="$PWD/$MATRIX_FILE"
else
    MATRIX_FILE="ci/matrix.txt"
fi

cd "$TS_ROOT"

# --- BUILD (idempotent: skip once ChaNGa.smp exists; skipped entirely for
# --bless-readydata, which neither runs nor needs ChaNGa) --------------------
if [ "$MODE" != blessready ] && [ ! -x "$CHANGA_SRC/ChaNGa.smp" ]; then
    mkdir -p "$EXT"
    # Development charm++ (default branch HEAD) + the utility repo (STRUCT_DIR).
    [ -d "$EXT/charm" ]   || git clone https://github.com/charmplusplus/charm.git "$EXT/charm"
    [ -d "$EXT/utility" ] || git clone https://github.com/N-bodyshop/utility "$EXT/utility"
    ( cd "$EXT/charm" && \
      { [ -x "$CHARM_TARGET/bin/charmrun" ] || \
        ./build ChaNGa netlrts-linux-x86_64 smp --with-production --enable-error-checking -j"$(nproc)"; } ) \
      || { echo "FATAL: charm build failed"; exit 1; }
    export CHARM_DIR="$EXT/charm" STRUCT_DIR="$EXT/utility/structures"
    # Build the intended ChaNGa branch (the source tree may be left on another).
    ( cd "$CHANGA_SRC" && git checkout "$CHANGA_BRANCH" ) \
      || { echo "FATAL: cannot checkout ChaNGa branch '$CHANGA_BRANCH' in $CHANGA_SRC"; exit 1; }
    ( cd "$CHANGA_SRC" && ./configure --enable-mhd --enable-mhdcleaning --enable-isph && make -j"$(nproc)" ) \
      || { echo "FATAL: ChaNGa build failed"; exit 1; }
    ln -sf "$CHANGA_SRC/ChaNGa" "$CHANGA_SRC/ChaNGa.smp"
    ln -sf "$CHARM_DIR/$CHARM_TARGET/bin/charmrun" "$CHANGA_SRC/charmrun"
    # (CHANGA_DIR + charmrun opts are passed to runtest.sh via the environment
    # above -- no runtest.sh source patching needed.)
fi

# --- per-test helpers (one call below for a single test; the matrix loop reuses
# them per row). Args: <test> <Nsmooth> <Nres> <distrib> <reg-tol> [-a..-h v ...]
# -- the trailing flags (if any) are forwarded to runtest.sh + the analysis. ---

# BLESS: runtest.sh -r runs the sim and saves the reference (auto-lands in
# test_cases/<test>/references/<runname>/; commit that folder). reg-tol unused.
# $6 is `extra`: runtest.sh's [extra] positional (e.g. "-n 50" to cap steps),
# forwarded ONLY to the run command (NOT the analysis), passed as vm=0 + extra.
do_bless() {
    local t=$1 ns=$2 nr=$3 d=$4 extra=$6; shift 6   # drop reg-tol ($5) + extra ($6)
    ./runtest.sh -r "$@" "$t" "$ns" "$nr" "$d" 1 0 "$extra" || { echo "FAIL[bless]: $t"; return 1; }
    echo "BLESSED: test_cases/$t/references/ (commit it)"
}

# BLESS-READYDATA: the run data already exists -- skip the sim, only run the
# analysis --save-reference on the existing run dir (this is exactly the step
# runtest.sh -r does after a run, so the reference auto-lands in
# test_cases/<test>/references/<runname>/). reg-tol unused.
do_bless_readydata() {
    local t=$1 ns=$2 nr=$3 d=$4; shift 6            # drop reg-tol ($5) + extra ($6, unused: no run)
    local rundir
    rundir=$(ls -dt "test_cases/$t/${t}${nr}_N${ns}"*CHANGA 2>/dev/null | head -1)
    [ -n "$rundir" ] || { echo "FAIL[no-rundir]: $t (no existing run to bless; run the sim first)"; return 1; }
    python "setupfiles/IC_analysis_${t}.py" "$rundir" "$@" --save-reference --save "$rundir/ci" \
        || { echo "FAIL[bless-readydata]: $t"; return 1; }
    echo "BLESSED (ready data): test_cases/$t/references/ (commit it)"
}

# GATE: run, then compare against the committed reference.
do_gate() {
    local t=$1 ns=$2 nr=$3 d=$4 tol=$5 extra=$6; shift 6   # $6 extra -> run only
    [ -n "$tol" ] || { echo "FAIL[no-tol]: $t (reg-tol required for gate mode)"; return 1; }
    ./runtest.sh "$@" "$t" "$ns" "$nr" "$d" 1 0 "$extra" || { echo "FAIL[run]: $t"; return 1; }
    local rundir
    rundir=$(ls -dt "test_cases/$t/${t}${nr}_N${ns}"*CHANGA 2>/dev/null | head -1)
    [ -n "$rundir" ] || { echo "FAIL[no-rundir]: $t"; return 1; }
    # --reference (no path) auto-discovers <rundir-parent>/references/<runname>/.
    python "setupfiles/IC_analysis_${t}.py" "$rundir" "$@" --reference --reg-tol "$tol" --save "$rundir/ci" \
        || { echo "FAIL[gate]: $t (residual > $tol)"; return 1; }
    echo "PASS: $t"
}

# Run every row of $MATRIX_FILE (mode: gate|bless|blessready). Rows are
# "<test> <Nsmooth> <Nres> <distrib> <reg-tol> [-a..-h v ...] [extra: <run-flags>]".
# The optional -a..-h flags go to BOTH runtest.sh + the analysis (so one test may
# appear on several rows, e.g. sedov -a 0/1/2). Anything after an "extra:" marker
# is runtest.sh's [extra] positional (e.g. "extra: -n 50" to cap steps) and is
# forwarded ONLY to the run, not the analysis. Runs ALL rows (does not stop on the
# first failure) and exits nonzero iff any failed, mirroring the CI workflow.
run_matrix() {
    local mode=$1 fails=0 ran=0 t ns nr d tol rest flags extra
    [ -f "$MATRIX_FILE" ] || { echo "FATAL: matrix file not found: $MATRIX_FILE"; exit 1; }
    while read -r t ns nr d tol rest <&3; do
        { [ -z "$t" ] || [ "${t#\#}" != "$t" ]; } && continue   # skip blanks/comments
        # Split the trailing tokens at the 'extra:' marker: flags before (-> run +
        # analysis), extra after (-> run only). No marker -> all flags, no extra.
        flags="$rest"; extra=""
        if [ "${rest#*extra:}" != "$rest" ]; then
            flags="${rest%%extra:*}"
            read -r extra <<<"${rest#*extra:}"          # trims surrounding space
        fi
        ran=$((ran + 1))
        case "$mode" in
            bless)      do_bless           "$t" "$ns" "$nr" "$d" "$tol" "$extra" $flags || fails=$((fails + 1)) ;;
            blessready) do_bless_readydata "$t" "$ns" "$nr" "$d" "$tol" "$extra" $flags || fails=$((fails + 1)) ;;
            *)          do_gate            "$t" "$ns" "$nr" "$d" "$tol" "$extra" $flags || fails=$((fails + 1)) ;;
        esac
    done 3< "$MATRIX_FILE"
    echo "----"
    echo "matrix $mode: $((ran - fails))/$ran passed"
    [ "$fails" -eq 0 ] || exit 1
}

# --- dispatch: no <test> -> run the whole matrix; else the single test -------
# (the build above already ran; per-row -a..-h live in matrix.txt, so a global
# -a..-h would be ambiguous across tests -> reject it in matrix mode.)
if [ -z "$TEST" ]; then
    [ "${#PARAM_FLAGS[@]}" -eq 0 ] || { echo "FATAL: -a..-h flags require an explicit <test> (matrix rows carry their own)"; exit 1; }
    run_matrix "$MODE"
    exit 0
fi

case "$MODE" in
    bless)      do_bless           "$TEST" "$NSM" "$NRES" "$DIST" "$TOL" "$EXTRA_CLI" "${PARAM_FLAGS[@]+"${PARAM_FLAGS[@]}"}" || exit 1 ;;
    blessready) do_bless_readydata "$TEST" "$NSM" "$NRES" "$DIST" "$TOL" "$EXTRA_CLI" "${PARAM_FLAGS[@]+"${PARAM_FLAGS[@]}"}" || exit 1 ;;
    *)          do_gate            "$TEST" "$NSM" "$NRES" "$DIST" "$TOL" "$EXTRA_CLI" "${PARAM_FLAGS[@]+"${PARAM_FLAGS[@]}"}" || exit 1 ;;
esac
