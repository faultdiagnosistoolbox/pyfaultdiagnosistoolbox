FROM quay.io/pypa/manylinux2014_x86_64:latest

RUN mkdir /py_env && mkdir /fdt_build

RUN /opt/python/cp38-cp38/bin/python3 -m venv /py_env/env38 && \
    . /py_env/env38/bin/activate && \
    pip install -U pip && \
    pip install setuptools wheel==0.34.2 auditwheel && \
    deactivate

RUN /opt/python/cp39-cp39/bin/python3 -m venv /py_env/env39 && \
    . /py_env/env39/bin/activate && \
    pip install -U pip && \
    pip install setuptools wheel==0.34.2 auditwheel && \
    deactivate

RUN /opt/python/cp310-cp310/bin/python3 -m venv /py_env/env310 && \
    . /py_env/env310/bin/activate && \
    pip install -U pip && \
    pip install setuptools wheel auditwheel && \
    deactivate

RUN /opt/python/cp311-cp311/bin/python3 -m venv /py_env/env311 && \
    . /py_env/env311/bin/activate && \
    pip install -U pip && \
    pip install setuptools wheel auditwheel && \
    deactivate

WORKDIR /fdt_build

CMD ./linux_create_wheels.sh
