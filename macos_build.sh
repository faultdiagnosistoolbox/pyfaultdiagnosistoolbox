#!/bin/sh

cd src

make cleanall

# Generate Python3.11
source ../fdt_311/bin/activate
make build_ext
make wheel
deactivate

# Generate Python3.10
source ../fdt_310/bin/activate
make build_ext
make wheel
deactivate

# Generate Python3.9
source ../fdt_39/bin/activate
make build_ext
make wheel
deactivate

# Generate Python3.8
source ../fdt_38/bin/activate
make build_ext
make wheel
deactivate


