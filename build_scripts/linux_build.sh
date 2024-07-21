#!/bin/sh
# docker build --tag build_fdt_x86 --platform=linux/amd64 -f Dockerfile_linux_x86 .
# docker run -v `pwd`:/fdt_build buildpyfdt
# docker run -it -v `pwd`:/fdt_build buildpyfdt /bin/bash

# Generate Python3.9
source /py_env/env39/bin/activate
python -m build --wheel > log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Generate Python3.10
source /py_env/env310/bin/activate
python -m build --wheel > log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Generate Python3.11
source /py_env/env311/bin/activate
python -m build --wheel > log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Generate Python3.12
. /py_env/env312/bin/activate
python -m build --wheel > log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Cleanup
rm -rf dist/wheelhouse
