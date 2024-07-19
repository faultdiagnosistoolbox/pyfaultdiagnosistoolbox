#!/bin/bash
# build suitable docker conrtainers for linux buils on x86 and arm
# docker build --tag build_fdt_x86 --platform linux/amd64 -f Dockerfile_linux_x86 .
# docker build --tag build_fdt_aarch64 -f Dockerfile_linux_aarch64 .

# docker run -it -v $(pwd):/fdt_build buildpyfdt_2014 bash

docker run --rm -v $(pwd):/fdt_build build_fdt_x86
docker run --rm -v $(pwd):/fdt_build build_fdt_aarch64

