#!/bin/sh

# Creat build venvs if necessary
if [ ! -d "build_envs" ]; then
    mkdir build_envs

    python3.12 -m venv build_envs/env312
    source build_envs/env312/bin/activate
    pip install -U pip build

    python3.13 -m venv build_envs/env313
    source build_envs/env313/bin/activate
    pip install -U pip build

    python3.14 -m venv build_envs/env314
    source build_envs/env314/bin/activate
    pip install -U pip build

    python3.15 -m venv build_envs/env315
    source build_envs/env315/bin/activate
    pip install -U pip build

fi

# Python 3.12
source build_envs/env312/bin/activate
python -m build --wheel
deactivate

# Python 3.13
source build_envs/env313/bin/activate
python -m build --wheel
deactivate

# Python 3.14
source build_envs/env314/bin/activate
python -m build --wheel
deactivate

# Python 3.15
source build_envs/env315/bin/activate
python -m build --wheel
deactivate

# Cleanup
# rm -rf build_envs
