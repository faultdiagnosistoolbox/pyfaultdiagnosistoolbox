from setuptools import setup, Extension
# To use a consistent encoding
# from codecs import open
from os import path  # , system
import numpy as np

here = path.abspath(path.dirname(__file__))

incdir = np.get_include()

batchfaultisolation = Extension('batchfaultisolation',
                                sources=['BatchFaultIsolation.cc'],
                                include_dirs=[incdir],
                                extra_compile_args=['-Wno-unused-function', '-Wno-unknown-pragmas'])

setup(
    name='batchfaultisolation',
    version='0.1',
    description='',

    # Author details
    author='Erik Frisk',
    author_email='erik.frisk@liu.se',

    ext_modules=[batchfaultisolation],
)
