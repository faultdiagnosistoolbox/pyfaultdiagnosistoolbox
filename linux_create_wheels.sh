#!/bin/sh

# docker run -v `pwd`:/fdt_build buildpyfdt
# docker run -it -v `pwd`:/fdt_build buildpyfdt /bin/bash

cd src

make cleanall

if [ ! "$AUDITWHEEL_PLAT" == "manylinux1_x86_64" ]; then
# Generate Python3.6
. /py_env/env36/bin/activate
make build_ext
python setup.py bdist_wheel > log.txt
wheel_name=$(grep "creating 'dist/" log.txt | sed "s/.*dist\/\(.*\)' and .*/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate


# Generate Python3.7
. /py_env/env37/bin/activate
python setup.py build_ext --inplace
python setup.py bdist_wheel > log.txt
wheel_name=$(grep "creating 'dist/" log.txt | sed "s/.*dist\/\(.*\)' and .*/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate
fi

# Generate Python3.8
. /py_env/env38/bin/activate
python setup.py build_ext --inplace
python setup.py bdist_wheel > log.txt
wheel_name=$(grep "creating 'dist/" log.txt | sed "s/.*dist\/\(.*\)' and .*/\1/g")
rm -f log.txt
cd dist
auditwheel repair "$wheel_name"
rm -f "$wheel_name"
mv -f wheelhouse/* .
cd ..
deactivate

rm -rf dist/wheelhouse

