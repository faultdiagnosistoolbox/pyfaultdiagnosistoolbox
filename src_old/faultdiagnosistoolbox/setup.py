from distutils.core import setup, Extension
import numpy as np

incdir = np.get_include()

module1 = Extension('structuralanalysis',
                    sources = ['structuralanalysismodule.cc'],
                    include_dirs=[incdir,'/Users/frisk/sw/CSParse/Include'],
                    extra_compile_args=['-Wno-unused-function'],
                    library_dirs=['/Users/frisk/sw/CSParse/Lib'],
                    libraries=['csparse'], extra_objects=['SparseMatrix.o', 'StructuralAnalysisModel.o', 'MSOalg.o'])

setup (name = 'structuralanalysis',
       version = '0.1',
       description = 'Simple module for basic structural analysis of sparse matrices',
       ext_modules = [module1]
)
