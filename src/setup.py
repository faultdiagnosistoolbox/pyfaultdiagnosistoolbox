from setuptools import setup, find_packages, Extension
# To use a consistent encoding
from codecs import open
from os import path, system
import platform

here = path.abspath(path.dirname(__file__))

# incdir = np.get_include()

class get_numpy_include(object):
    """Defer numpy.get_include() until after numpy is installed."""

    def __str__(self):
        import numpy
        return numpy.get_include()
    

if platform.system() == "Windows":
    extra_compile = []
else:
    extra_compile = ['-Wno-unused-function', '-Wno-unknown-pragmas']

strucanalysis_ext = Extension('faultdiagnosistoolbox.structuralanalysis',
                    sources = ['faultdiagnosistoolbox/structuralanalysismodule.cc', 'faultdiagnosistoolbox/SparseMatrix.cc',
                               'faultdiagnosistoolbox/StructuralAnalysisModel.cc', 'faultdiagnosistoolbox/MSOAlg.cc'],
                    include_dirs=[get_numpy_include(),'CSparse/Include'],
                    extra_compile_args=extra_compile,
                    library_dirs=['CSparse/Lib'],
                    libraries=['csparse'])

if platform.system() == "Windows":
    if not path.isfile('CSparse/Lib/libcsparse.lib'):
        system('cd CSparse/Lib & nmake -f Makefile.win')
else:
    if not path.isfile('CSparse/Lib/libcsparse.a'):
        system('(cd CSparse; MACOSX_DEPLOYMENT_TARGET=10.6 make)')

readme = open('README.md', 'r')
README_TEXT = readme.read()
readme.close()

setup(
    name='faultdiagnosistoolbox',
    version='0.12',

    description='A Fault Diagnosis Toolbox',
    long_description=README_TEXT,

    url='https://faultdiagnosistoolbox.github.io',

    # Author details
    author='Erik Frisk',
    author_email='erik.frisk@liu.se',

    ext_modules = [strucanalysis_ext],

    include_package_data = True, 
    
    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
    ],

    # What does your project relate to?
    keywords='model based fault diagnosis, signal processing',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    setup_requires=['cython', 'numpy', 'wheel'],
    install_requires=['sympy', 'numpy', 'scipy', 'matplotlib', 'sklearn', 'wheel'],
    python_requires='>=3.6'
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
#    extras_require={
#        'dev': ['check-manifest'],
#        'test': ['coverage'],
#    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
#    package_data={
#        'sample': ['package_data.dat'],
#    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
 #   data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #},
)
