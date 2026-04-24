#!/usr/bin/env bash
set -euo pipefail

LIBSTDCXX="/usr/lib/x86_64-linux-gnu/libstdc++.so.6"
LIBPYTHON="/usr/lib/x86_64-linux-gnu/libpython3.10.so"
MATLAB_CMD="${MATLAB_CMD:-matlab}"

export LD_PRELOAD="${LIBSTDCXX} ${LIBPYTHON}"

exec "$MATLAB_CMD"
