#!/bin/bash
# docker build -t buildpyfdt_2010 --platform linux/arm64 -f 2010.Dockerfile .
# docker build -t buildpyfdt_2014 --platform linux/arm64 -f 2014.Dockerfile .

# docker run -it -v `pwd`:/fdt_build buildpyfdt_2014 /bin/bash
# docker run -it -v `pwd`:/fdt_build buildpyfdt_2010 /bin/bash

docker run -v `pwd`:/fdt_build buildpyfdt_2010
docker run -v `pwd`:/fdt_build buildpyfdt_2014

