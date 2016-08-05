from distutils.core import setup, Extension
import numpy as np

incdir = np.get_include()

module1 = Extension('test',
                    include_dirs=[incdir],
                    sources = ['testmodule.c'],
                    extra_compile_args=['-Wno-unused-function'])

setup (name = 'test',
       version = '0.1',
       description = 'Simple module working both in Python 2 and 3',
       ext_modules = [module1]
)
