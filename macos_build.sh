#!/bin/sh

cd src

make cleanall

# Generate Python3.9
. ~/sw/pyenv/39/bin/activate
make build_ext
python setup.py bdist_wheel
deactivate

# Generate Python3.8
. ~/sw/pyenv/38/bin/activate
make build_ext
python setup.py bdist_wheel
deactivate

# Generate Python3.7
. ~/sw/pyenv/37/bin/activate
make build_ext
python setup.py bdist_wheel
deactivate


