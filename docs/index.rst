.. title:: faultdiagnosistoolbox

Fault Diagnosis Toolbox
=======================

``Fault Diagnosis Toolbox`` is a Python package for the analysis and design of fault diagnosis systems for dynamic systems, primarily described by differential-algebraic equations. Key features of the toolbox are extensive support for structural analysis of large-scale dynamic models, fault isolability analysis, sensor placement analysis, and code generation in C/C++ and Python.

The code is available at `GitHub <https://github.com/faultdiagnosistoolbox/pyfaultdiagnosistoolbox>`_ and can be installed via `pip <https://pypi.org/project/faultdiagnosistoolbox/>`_. The toolbox is freely available under an MIT license.

If you use this toolbox in your research, please cite 

   Erik Frisk, Mattias Krysander, and Daniel Jung. "`A Toolbox for Analysis and Design of Model Based Diagnosis Systems for Large Scale Models <https://doi.org/10.1016/j.ifacol.2017.08.504>`_",
   IFAC World Congress. Toulouse, France, 2017.

For a quick introduction, see the :doc:`use-case </usecase>` where an industrial size example, an automotive engine, is analyzed, C-code for residual generators is generated, and the resulting diagnosis system is evaluated on test-cell measurements from our engine laboratory.

The toolbox is an adaptation of our Matlab toolbox that can be downloaded from https://faultdiagnosistoolbox.github.io/ where you can also find additional Matlab documentation.

Key references
--------------

    E. Frisk, M. Krysander, and D. Jung. "*A Toolbox for Analysis and Design of Model Based Diagnosis Systems for Large Scale Models*" (https://doi.org/10.1016/j.ifacol.2017.08.504), IFAC World Congress. Toulouse, France, 2017.

    M. Krysander, Åslund, J., Nyberg, M. (2007). "*An efficient algorithm for finding minimal overconstrained subsystems for model-based diagnosis*" (https://doi.org/10.1109/TSMCA.2007.909555). IEEE Transactions on Systems, Man, and Cybernetics-Part A: Systems and Humans, 38(1), 197-206.

    M. Krysander, Åslund, J., Frisk, E. (2010). "*A structural algorithm for finding testable sub-models and multiple fault isolability analysis*". In 21st International Workshop on Principles of Diagnosis (DX-10), Portland, Oregon, USA (pp. 17-18).

    M. Krysander and E. Frisk. "*Sensor placement for fault diagnosis*" (https://doi.org/10.1109/TSMCA.2008.2003968). Systems, Man and Cybernetics, Part A: Systems and Humans, IEEE Transactions on, 38(6):1398-1410, 2008.

    E. Frisk, A. Bregon, J. Åslund, M. Krysander, B. Pulido, and G. Biswas. "*Diagnosability analysis considering causal interpretations for differential constraints*" (https://doi.org/10.1109/TSMCA.2012.2189877). Systems, Man and Cybernetics, Part A: Systems and Humans, IEEE Transactions on, 42(5):1216-1229, 2012.

See also additional references at https://faultdiagnosistoolbox.github.io/references/.


.. raw:: html

   <hr>

Table of Contents
-----------------

.. toctree::
   :maxdepth: 1

   Basic Usage <basicusage>
   Examples <examples/examples>
   Larger Use case <usecase>
   API <faultdiagnosistoolbox>
