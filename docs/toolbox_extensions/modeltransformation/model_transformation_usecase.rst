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


Induction motor example
-----------------------

The induction motor example in ``src/code_examples/InductionMotor.py`` can be
used in the same way. The model contains products between electrical and
mechanical variables, auxiliary torque/flux terms, and differential
constraints. Applying the transformation before constructing the
``DiagnosisModel`` lets the toolbox introduce any additional repeated
subexpressions it finds.

.. code-block:: python

   import faultdiagnosistoolbox as fdt
   import sympy as sym


   model_def = {
       "type": "Symbolic",
       "x": [
           "i_a", "i_b", "lambda_a", "lambda_b", "w",
           "di_a", "di_b", "dlambda_a", "dlambda_b", "dw",
           "q_a", "q_b", "Tl",
       ],
       "f": ["f_a", "f_b"],
       "z": ["u_a", "u_b", "y1", "y2", "y3"],
       "parameters": ["a", "b", "c", "d", "L_M", "k", "c_f", "c_t"],
   }

   sym.var(model_def["x"])
   sym.var(model_def["f"])
   sym.var(model_def["z"])
   sym.var(model_def["parameters"])

   model_def["rels"] = [
       -q_a + w * lambda_a,
       -q_b + w * lambda_b,
       -di_a - a * i_a + b * c * lambda_a + b * q_b + d * u_a,
       -di_b - a * i_b + b * c * lambda_b + b * q_a + d * u_b,
       -dlambda_a + L_M * c * i_a - c * lambda_a - q_b,
       -dlambda_b + L_M * c * i_b - c * lambda_b - q_a,
       -dw - k * c_f * w + k * c_t * (i_a * lambda_b - i_b * lambda_a) - k * Tl,
       fdt.DiffConstraint("di_a", "i_a"),
       fdt.DiffConstraint("di_b", "i_b"),
       fdt.DiffConstraint("dlambda_a", "lambda_a"),
       fdt.DiffConstraint("dlambda_b", "lambda_b"),
       fdt.DiffConstraint("dw", "w"),
       -y1 + i_a + f_a,
       -y2 + i_b + f_b,
       -y3 + w,
   ]

   transformed = fdt.ModelTransforming.optimize(model_def)
   print(transformed["summary"])

   model = fdt.DiagnosisModel(
       transformed["model_def"],
       name="Transformed induction motor",
   )
   model.Lint()

Once the transformed model has been created, the remaining analysis is the same
as in the original induction motor example:

.. code-block:: python

   msos = model.MSO()
   mtes = model.MTES()
   print(f"Found {len(msos)} MSO sets and {len(mtes)} MTES sets.")

   oi_mso = [model.IsObservable(m_i) for m_i in msos]
   li_mso = [model.IsLowIndex(m_i) for m_i in msos]
   print(
       f"Out of {len(msos)} MSO sets, "
       f"{sum(oi_mso)} observable, "
       f"{sum(li_mso)} low (structural) differential index"
   )

   model.MSOCausalitySweep(mtes[0])
   red_eq = mtes[0][10]
   matching_equations = [e for e in mtes[0] if e != red_eq]

   Gamma = model.Matching(matching_equations)
   model.SeqResGen(Gamma, red_eq, "ResGen", batch=True, language="C")
