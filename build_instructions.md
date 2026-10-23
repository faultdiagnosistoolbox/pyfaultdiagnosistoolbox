# Build instructions

## Build wheels
Binaries for the 4 latest Python releases shall be compiled for platforms MacOS, Windows, Linux-arm and Linux x86.

### Linux ARM and x86
Built on a MacOS host using Docker.
1. Update `Dockerfile_linux_aarch64`, `Dockerfile_linux_x86`, and `build_scripts/linux_build.sh`.
1. Create docker images for ARM and x86 by running
    ```bash
    cd Docker
    docker build --tag build_fdt_arm -f Dockerfile_linux_aarch64 .
    docker build --tag build_fdt_x86 --platform=linux/amd64 -f Dockerfile_linux_x86 .
    ```
3. Run the build scripts for ARM and x86 by running
    ```bash
    docker run --rm -it -v `pwd`:/fdt_build build_fdt_arm /bin/bash
    docker run --rm -it -v `pwd`:/fdt_build build_fdt_x86 /bin/bash
    ```

### MacOS
1. Update `build_scripts/macos_build.sh`.
1. Run the build script by running
    ```bash
    ./build_scripts/macos_build.sh
    ```

### Windows


### Send to pypi


## Build documentation
The documentation is built using Sphinx using the command
```bash
sphinx-build -M html docs doc_build
```
and build dependencies can be found in `docs/requirements-doc.txt`.

## Build development environment
Create a venv and install dependencies defined in `requirements-develop.txt`. Install the toolbox as editable package by running
```bash
pip install -e .
```

## Run test suite
With a development environment active, run the test-suite using `pytest` by running
```bash
python -m pytest -v --cov=faultdiagnosistoolbox
```
or
```bash
python -m pytest -v --cov=faultdiagnosistoolbox --cov-report=html
````
