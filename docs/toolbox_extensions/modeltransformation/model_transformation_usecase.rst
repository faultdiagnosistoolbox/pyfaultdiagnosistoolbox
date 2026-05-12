Model Transformation Use Cases
==============================

The examples below show how to apply the model transformation extension before
creating a ``DiagnosisModel``. The workflow is:

1. define a symbolic ``model_def``
2. call ``fdt.ModelTransforming.optimize``
3. create the diagnosis model from the transformed dictionary
4. run the normal toolbox analysis methods

Minimal repeated-subexpression example
--------------------------------------

This example contains the repeated expression ``x1*x2`` in several relations.
The transformation introduces an auxiliary variable for the repeated expression
and rewrites the affected equations.

.. code-block:: python

   import faultdiagnosistoolbox as fdt
   import sympy as sym


   model_def = {
       "type": "Symbolic",
       "x": ["x1", "x2", "x3"],
       "f": ["f1", "f2", "f3"],
       "z": ["y1", "y2", "y3"],
   }

   sym.var(model_def["x"])
   sym.var(model_def["f"])
   sym.var(model_def["z"])

   model_def["rels"] = [
       -y1 + x1 * x2 + x3 + f1,
       -y2 + x1 * x2 + x3 + f2,
       -y3 + x1 * x2 + f3,
   ]

   transformed = fdt.ModelTransforming.optimize(model_def)
   print(transformed["summary"])

   model = fdt.DiagnosisModel(
       transformed["model_def"],
       name="Transformed repeated-subexpression model",
   )
   model.Lint()

The summary lists the generated replacement equations and the rewritten model
relations. The original ``model_def`` is copied before transformation, so it can
still be used unchanged if needed.

Using only the transformed model definition
-------------------------------------------

When the summary is not needed, use ``get_optimized_modelDef``:

.. code-block:: python

   optimized_model_def = fdt.ModelTransforming.get_optimized_modelDef(model_def)
   model = fdt.DiagnosisModel(optimized_model_def)


Electric motor example
----------------------

The electric motor example in ``src/code_examples/ElectricMotor.py`` applies
the transformation before running structural analysis, isolability analysis,
and residual generator design.

.. code-block:: python

   import faultdiagnosistoolbox as fdt
   import sympy as sym


   model_def = {
       "type": "Symbolic",
       "x": ["dIc", "dw", "dth", "Ic", "w", "th", "alpha", "DT", "Tm", "Tl"],
       "f": ["fR", "fi", "fw", "fD"],
       "z": ["V", "yi", "yw", "yd"],
       "parameters": ["Ka", "b", "R", "J", "L"],
   }

   sym.var(model_def["x"])
   sym.var(model_def["f"])
   sym.var(model_def["z"])
   sym.var(model_def["parameters"])

   model_def["rels"] = [
       -V + Ic * (R + fR) + L * dIc + Ka * Ic * w,
       -Tm + Ka * Ic**2,
       -J * dw + DT - b * w,
       -DT + Tm - Tl,
       -dth + w,
       -dw + alpha,
       -yi + Ic + fi,
       -yw + w + fw,
       -yd + DT + fD,
       fdt.DiffConstraint("dIc", "Ic"),
       fdt.DiffConstraint("dw", "w"),
       fdt.DiffConstraint("dth", "th"),
   ]

   transformed = fdt.ModelTransforming.optimize(model_def)
   print(transformed["summary"])

   model = fdt.DiagnosisModel(
       transformed["model_def"],
       name="Electric Motor",
   )

   model.Lint()
   model.IsolabilityAnalysis(causality="der")

After this point the transformed model can be used with the same toolbox
methods as any other ``DiagnosisModel``.
