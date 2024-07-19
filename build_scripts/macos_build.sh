#!/bin/sh

# Creat build venvs if necessary
if [ ! -d "build_envs" ]; then
    mkdir build_envs

    python3.9 -m venv build_envs/env39
    source build_envs/env39/bin/activate
    pip install -U pip build

    python3.10 -m venv build_envs/env310
    source build_envs/env310/bin/activate
    pip install -U pip build

    python3.11 -m venv build_envs/env311
    source build_envs/env311/bin/activate
    pip install -U pip build

    python3.12 -m venv build_envs/env312
    source build_envs/env312/bin/activate
    pip install -U pip build
fi

# Python 3.9
source build_envs/env39/bin/activate
python -m build --wheel
deactivate

# Python 3.10
source build_envs/env310/bin/activate
python -m build --wheel
deactivate

# Python 3.11
source build_envs/env311/bin/activate
python -m build --wheel
deactivate

# Python 3.12
source build_envs/env312/bin/activate
python -m build --wheel
deactivate

# Cleanup
# rm -rf build_envs
