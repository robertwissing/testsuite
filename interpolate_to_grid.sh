#!/usr/bin/env bash
# Interpolate a tipsy SPH snapshot onto a grid via the sph_interp package.
#
#   ./interpolate_to_grid.sh INPUT OUTPUT [options]
#   ./interpolate_to_grid.sh --help
#   ./interpolate_to_grid.sh            # no args -> interactive prompts
#
# Thin wrapper: puts setupfiles/ on PYTHONPATH and runs the Python CLI.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$SCRIPT_DIR/setupfiles:${PYTHONPATH:-}"
exec python "$SCRIPT_DIR/setupfiles/IC_interpolate_to_grid.py" "$@"
