FROM quay.io/pypa/manylinux2010_x86_64:latest

RUN mkdir /py_env && mkdir /fdt_build

RUN /opt/python/cp36-cp36m/bin/python3 -m venv /py_env/env36 && \
    . /py_env/env36/bin/activate && \
    pip install -U pip && \
    pip install setuptools wheel auditwheel && \
    deactivate

RUN /opt/python/cp37-cp37m/bin/python3 -m venv /py_env/env37 && \
    . /py_env/env37/bin/activate && \
    pip install -U pip && \
    pip install setuptools wheel auditwheel && \
    deactivate

RUN /opt/python/cp38-cp38/bin/python3 -m venv /py_env/env38 && \
    . /py_env/env38/bin/activate && \
    pip install -U pip && \
    pip install setuptools wheel auditwheel && \
    deactivate

WORKDIR /fdt_build

CMD ./linux_build.sh
