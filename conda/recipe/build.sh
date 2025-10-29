#!/usr/bin/env bash
# Build script for conda (unix)
set -euo pipefail

# Install the package into the build environment using pip.
"${PYTHON}" -m pip install . --no-deps -vv
