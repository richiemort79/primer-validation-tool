#!/usr/bin/env bash
# Create (or update) the Python environment for the primer validation tool.
#
# The environment lives outside the repo (default ~/.venvs/primer_validation_tool)
# because a venv synced by Dropbox/OneDrive breaks on other machines.
#
#   ./setup.sh          # create/update the environment
#   ./setup.sh --dev    # also install pytest for running the tests
#
# Override the location with PRIMER_TOOL_VENV=/path/to/venv.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
VENV="${PRIMER_TOOL_VENV:-$HOME/.venvs/primer_validation_tool}"

PYTHON=""
for candidate in python3 python; do
    if command -v "$candidate" >/dev/null 2>&1 &&
       "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 9))' 2>/dev/null; then
        PYTHON="$candidate"
        break
    fi
done
if [ -z "$PYTHON" ]; then
    echo "❌ Python 3.9 or newer is required but was not found on PATH." >&2
    exit 1
fi

if [ ! -x "$VENV/bin/python" ]; then
    echo "Creating environment in $VENV (using $("$PYTHON" --version))"
    mkdir -p "$(dirname "$VENV")"
    "$PYTHON" -m venv "$VENV"
fi

echo "Installing requirements..."
"$VENV/bin/python" -m pip install --quiet --upgrade pip
"$VENV/bin/python" -m pip install --quiet -r "$REPO_DIR/requirements.txt"
if [ "${1:-}" = "--dev" ]; then
    "$VENV/bin/python" -m pip install --quiet pytest
fi

# Records which requirements.txt was installed, so ./check_primers can tell
# when it needs to re-run setup.
cp "$REPO_DIR/requirements.txt" "$VENV/.primer_tool_requirements"

echo "✓ Environment ready: $VENV"
echo "  Run the tool with: $REPO_DIR/check_primers primers.xlsx"
