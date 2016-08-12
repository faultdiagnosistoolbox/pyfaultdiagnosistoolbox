from distutils.core import setup, Extension
import numpy as np

incdir = np.get_include()

module1 = Extension('structuralanalysis',
                    sources = ['structuralanalysismodule.cc', 'SparseMatrix.cc', 'StructuralAnalysisModel.cc', 'MSOAlg.cc'],
                    include_dirs=[incdir,'../CSparse/Include'],
                    extra_compile_args=['-Wno-unused-function', '-Wno-unknown-pragmas'],
                    library_dirs=['../CSparse/Lib'],
                    libraries=['csparse'])
# , extra_objects=['SparseMatrix.o', 'StructuralAnalysisModel.o', 'MSOAlg.o']
setup (name = 'structuralanalysis',
       version = '0.1',
       description = 'Simple module for basic structural analysis of sparse matrices',
       ext_modules = [module1]
)
