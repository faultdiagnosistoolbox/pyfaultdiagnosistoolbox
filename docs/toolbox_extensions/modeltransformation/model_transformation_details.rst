Model Transformation Extension
==============================

Overview
--------
The model transformation extensions enables equation optimization for the
given model. Repeated subexpressions are identified and replaced with new 
equation variables, allowing for better compatibility with Fault Diagnosis 
Toolbox. 


For example, if several relations contain ``x1*x2``, the transformed model can
introduce a new variable such as ``x0`` and add the relation
``-x0 + x1*x2``. The original relations are then rewritten in terms of
``x0``. The resulting model is still a standard
``faultdiagnosistoolbox.DiagnosisModel`` input.

Subexpression identification
----------------------------

The extension utilizes the ``cse`` function from the SymPy Python library for
symbolic computations in the function ``analyze_modelDef`` to identify
possible subexpression reductions and optimizations. All replacement equations
of at least 2 variables are kept and returned from the function, as well as
the reduced equations.


Only replacements that depend on at least two original symbols are kept.
Simpler replacements are substituted back into the reduced relations. This
keeps the transformed model from being cluttered with auxiliary variables that
do not improve the model structure.

Variable replacements
---------------------

The main entry point is ``optimize``. It first creates a deep copy of the input
dictionary with ``create_model_copy`` so that the original model definition is
not modified. It then calls ``analyze_modelDef`` and rebuilds the relation list
in the copied model:

1. each kept replacement symbol is added to ``model_def["x"]``,
2. each replacement is added as a new zero-equality relation, and
3. the reduced original relations are appended after the replacement
   relations.

The replacement relation for a pair ``(x0, x1*x2)`` is written as
``-x0 + x1*x2`` because relations in a symbolic toolbox model are expected to
equal zero.

The return value from ``optimize`` is a dictionary with two keys:

``model_def``
   The transformed model definition.

``summary``
   A printable summary of the replacement equations, reduced relations, and
   basic statistics.

For workflows that only need the transformed dictionary,
``get_optimized_modelDef`` returns the ``model_def`` value directly.

Progress output
---------------

The transformation functions print progress while running. This is useful for
larger symbolic models where common-subexpression elimination can take
noticeable time.

Usage and limitations
---------------------

The extension expects a symbolic ``model_def`` dictionary with ``rels`` stored
as a list. Relations that support SymPy replacement, such as SymPy
expressions, are transformed. Other relation objects are left unchanged by the
replacement helper.

The transformation changes the equation set by adding auxiliary variables and
relations. It should therefore be applied before constructing the final
``DiagnosisModel`` used for structural analysis, residual generation, or other
toolbox workflows.
