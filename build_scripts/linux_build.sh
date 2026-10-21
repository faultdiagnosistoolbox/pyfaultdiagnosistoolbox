#!/bin/sh
# docker build --tag build_fdt_arm -f Dockerfile_linux_aarch64 .
# docker build --tag build_fdt_x86 --platform=linux/amd64 -f Dockerfile_linux_x86 .
# docker run --rm -it -v `pwd`:/fdt_build build_fdt_arm /bin/bash


# Generate Python3.12
. /py_env/env312/bin/activate
python -m build --wheel | tee log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Generate Python3.13
. /py_env/env313/bin/activate
python -m build --wheel | tee log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Generate Python3.14
. /py_env/env314/bin/activate
python -m build --wheel | tee log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Generate Python3.14t
. /py_env/env314t/bin/activate
python -m build --wheel | tee log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Generate Python3.15
. /py_env/env315/bin/activate
python -m build --wheel | tee log.txt
wheel_name=$(grep "Successfully" log.txt | sed "s/Successfully built \(.*\.whl\)/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

# Generate Python3.15t
. /py_env/env315t/bin/activate
python -m build --wheel | tee log.txt
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
