# Fault Diagnosis Toolbox
<p>
<img alt="PyPI - Python Version" src="https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Ffaultdiagnosistoolbox%2Fpyfaultdiagnosistoolbox%2Fmaster%2Fpyproject.toml" />
<img alt="PyPI - License" src="https://img.shields.io/pypi/l/faultdiagnosistoolbox" />
</p>

Fault Diagnosis Toolbox is a Python package for the analysis and design of fault diagnosis systems for dynamic systems, 
primarily described by differential algebraic equations. Key features of the toolbox are extensive support for 
structural analysis of large-scale dynamic models, fault isolability analysis, sensor placement analysis, 
and code generation in C/C++ and Python. The toolbox is an adaptation of our Matlab toolbox that can be downloaded from https://faultdiagnosistoolbox.github.io.

The toolbox is freely available under an MIT license. 

## Documentation
The documentation is available at [readthedocs](https://faultdiagnosistoolbox.readthedocs.io/).

For a quick introduction, see the [use case](https://faultdiagnosistoolbox.readthedocs.io/en/latest/usecase.html) where an industrial 
size example, an automotive engine, is analyzed, C-code for residual generators is generated, and the resulting 
diagnosis system is evaluated on test-cell measurements from our engine laboratory.

## Installation 
The toolbox is available on [pip](https://pypi.org/project/faultdiagnosistoolbox/) and can be installed as:
```
pip install faultdiagnosistoolbox
```

## Publications

If you use this toolbox in your research, please cite 

   Erik Frisk, Mattias Krysander, and Daniel Jung. "_A Toolbox for Analysis and Design of Model Based Diagnosis Systems for Large Scale Models_" (https://doi.org/10.1016/j.ifacol.2017.08.504), IFAC World Congress. Toulouse, France, 2017.

and any relevant papers of ours. See a list of [key references](https://faultdiagnosistoolbox.readthedocs.io/) for details.
