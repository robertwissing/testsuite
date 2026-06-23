#!/usr/bin/env bash
# Showcase sph_interp across every target geometry on one snapshot.
#
#   ./interp_demo.sh [SNAPSHOT] [options]
#   ./interp_demo.sh --help
#   ./interp_demo.sh            # no args -> a default high-contrast snapshot
#
# Renders particles | UniformGrid | AMRGrid | VoronoiGrid | Projection2D side by
# side. Thin wrapper: puts setupfiles/ on PYTHONPATH and runs the Python demo.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$SCRIPT_DIR/setupfiles:${PYTHONPATH:-}"
exec python "$SCRIPT_DIR/setupfiles/IC_interp_demo.py" "$@"
