from setuptools import setup, find_packages, Extension
# To use a consistent encoding
from codecs import open
from os import path, system
import platform
import re

here = path.abspath(path.dirname(__file__))

VERSIONFILE = "faultdiagnosistoolbox/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    ver_str = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


class get_numpy_include(object):
    """Defer numpy.get_include() until after numpy is installed."""

    def __str__(self):
        import numpy
        return numpy.get_include()


if platform.system() == "Windows":
    extra_compile = []
else:
    extra_compile = ['-Wno-unused-function', '-Wno-unknown-pragmas', '--std=c++11']

strucanalysis_ext = Extension('faultdiagnosistoolbox.structuralanalysis',
                              sources=['faultdiagnosistoolbox/structuralanalysismodule.cc',
                                       'faultdiagnosistoolbox/SparseMatrix.cc',
                                       'faultdiagnosistoolbox/StructuralAnalysisModel.cc',
                                       'faultdiagnosistoolbox/MSOAlg.cc'],
                              include_dirs=[get_numpy_include(), 'CSparse/Include'],
                              extra_compile_args=extra_compile,
                              library_dirs=['CSparse/Lib'],
                              libraries=['csparse'])

if platform.system() == "Windows":
    if not path.isfile('CSparse/Lib/libcsparse.lib'):
        system('cd CSparse/Lib & nmake -f Makefile.win')
elif platform.system() == "Darwin":
    if not path.isfile('CSparse/Lib/libcsparse.a'):
        if platform.processor() == "arm":
            system('(cd CSparse/Lib; TARGET_ARCH="-arch arm64 -arch x86_64" CFLAGS=-mmacosx-version-min=11 make)')
        else:
            system('(cd CSparse/Lib; CFLAGS=-mmacosx-version-min=11 make)')
else:  # Linux
    if not path.isfile('CSparse/Lib/libcsparse.a'):
        system('(cd CSparse/Lib; make)')

readme = open('README.md', 'r')
README_TEXT = readme.read()
readme.close()

setup(
    name='faultdiagnosistoolbox',
    version=ver_str,

    description='A toolbox for Model Based Fault Diagnosis of dynamic systems based on structural analysis',
    long_description_content_type="text/x-rst",
    long_description=README_TEXT,

    url='https://faultdiagnosistoolbox.github.io',

    # Author details
    author='Erik Frisk',
    author_email='erik.frisk@liu.se',

    ext_modules=[strucanalysis_ext],

    include_package_data=True,

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
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
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
    setup_requires=['numpy', 'wheel'],
    install_requires=['sympy', 'numpy', 'scipy', 'matplotlib', 'scikit-learn', 'wheel'],
    python_requires='>=3.8'
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
    # entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    # },
)
