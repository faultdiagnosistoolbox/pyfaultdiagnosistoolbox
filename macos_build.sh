#!/bin/sh

cd src
make cleanall
cd ..

# Generate Python3.12
source fdt_312/bin/activate
python -m build --wheel
deactivate

# Generate Python3.11
source fdt_311/bin/activate
python -m build --wheel
deactivate

# Generate Python3.10
source fdt_310/bin/activate
python -m build --wheel
deactivate
