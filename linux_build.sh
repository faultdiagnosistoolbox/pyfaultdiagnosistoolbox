#!/bin/bash

# docker run -it -v `pwd`:/fdt_build buildpyfdt /bin/bash
docker run -v `pwd`:/fdt_build buildpyfdt_1
docker run -v `pwd`:/fdt_build buildpyfdt_2010
docker run -v `pwd`:/fdt_build buildpyfdt_2014

