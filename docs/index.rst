.. title:: faultdiagnosistoolbox

Fault Diagnosis Toolbox
=======================

``Fault Diagnosis Toolbox`` is a Python package for analysis and design of fault diagnosis systems for dynamic systems,
primarily described by differential-algebraic equations. Key features of the toolbox are extensive support for
structural analysis of large-scale dynamic models, fault isolability analysis, sensor placement analysis,
and code generation in C/C++ and Python.

For a quick introduction, see the :doc:`use-case </usecase>` where an industrial size example, an automotive engine, is analyzed, C-code for residual generators is generated, and the resulting diagnosis system is evaluated on test-cell measurements from our engine laboratory.

If you use this toolbox in your research, please cite 

   Erik Frisk, Mattias Krysander, and Daniel Jung. "`A Toolbox for Analysis and Design of Model Based Diagnosis Systems for Large Scale Models <https://doi.org/10.1016/j.ifacol.2017.08.504>`_",
   IFAC World Congress. Toulouse, France, 2017.

The toolbox is freely available under a MIT license.

The toolbox is an adaptation of our Matlab toolbox that can be downloaded from https://faultdiagnosistoolbox.github.io where you
can also find `documentation <https://faultdiagnosistoolbox.github.io/_releases/user-manual_2018-12-09.pdf>`_. The latest version of the python package can always be pip-installed.

.. raw:: html

   <hr>

Table of Contents
-----------------

.. toctree::
   :maxdepth: 1

   Introduction <introduction>
   Basic Usage <basicusage>
   Use case <usecase>
   Fault Diagnosis Toolbox Reference <faultdiagnosistoolbox>
