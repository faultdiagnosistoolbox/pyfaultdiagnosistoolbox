External Simulation Use Cases
=============================

The examples below show how to generate FMUs from three different models. The
examples follow the same workflow:

1. import a model,
2. select a residual equation from an MTES,
3. compute a matching for the remaining equations, and
4. call ``GenerateFMU``.

The ``model_path`` argument must point to the Python file that defines the model
object used by the FMU. For imported models, ``inspect.getfile`` keeps that path
aligned with the imported model module. For a standalone script, ``__file__`` can
be used as the model path.

Selecting Indices
-----------------

``model.MTES()`` returns a list of minimal test equation supports. The MTES index
is the zero-based position in that list. For example, ``mtes[1]`` selects the
second MTES returned by the model.

Each MTES is itself a list of equation indices from the model. One equation in
the selected MTES is used as the residual equation, ``res_eq``. The remaining
equations are passed to ``model.Matching`` to compute ``Gamma``. In the examples
below, ``selected_mtes[0]`` means "use the first equation in the selected MTES as
the residual equation".

Electric Motor
--------------

This example defines the complete electric motor model in the same script that
generates ``generated/electric_motor_example.fmu``. It uses MTES index ``1`` and
the first equation in that MTES as the residual equation.

.. code-block:: python

   from pathlib import Path

   import faultdiagnosistoolbox as fdt
   import sympy as sym


   model_def = {
       "type": "Symbolic",
       "x": ["dI", "dw", "I", "w", "T", "Tm", "Tl"],
       "f": ["fR", "fi", "fw", "fT"],
       "z": ["V", "yi", "yw", "yT"],
       "parameters": ["Ka", "b", "R", "J", "L"],
   }

   dI, dw, I, w, T, Tm, Tl = sym.symbols(model_def["x"])
   fR, fi, fw, fT = sym.symbols(model_def["f"])
   V, yi, yw, yT = sym.symbols(model_def["z"])
   Ka, b, R, J, L = sym.symbols(model_def["parameters"])

   model_def["rels"] = [
       -V + I * (R + fR) + L * dI + Ka * I * w,
       -Tm + Ka * I**2,
       -J * dw + T - b * w,
       -T + Tm - Tl,
       -yi + I + fi,
       -yw + w + fw,
       -yT + T + fT,
       fdt.DiffConstraint("dI", "I"),
       fdt.DiffConstraint("dw", "w"),
   ]

   model = fdt.DiagnosisModel(model_def, name="Electric motor")

   mtes = model.MTES()
   # Select the second MTES. Python list indices start at zero.
   selected_mtes = mtes[1]
   # Select the first equation in that MTES as the residual equation.
   res_eq = selected_mtes[0]
   matching_equations = [eq for eq in selected_mtes if eq != res_eq]
   Gamma = model.Matching(matching_equations)

   fmu_path = model.GenerateFMU(
       model_path=Path(__file__),
       Gamma=Gamma,
       res_eq=res_eq,
       fmu_name="electric_motor_example",
   )

   print(fmu_path)

Three-Tank System
-----------------

This example imports the three-tank model and generates
``generated/three_tank_example.fmu``. It uses MTES index ``2``.

.. code-block:: python

   import inspect

   from faultdiagnosistoolbox.models import _three_tank_model


   model = _three_tank_model.model
   model_path = inspect.getfile(_three_tank_model)

   mtes = model.MTES()
   # Select the third MTES. Python list indices start at zero.
   selected_mtes = mtes[2]
   # Select the first equation in that MTES as the residual equation.
   res_eq = selected_mtes[0]
   matching_equations = [eq for eq in selected_mtes if eq != res_eq]
   Gamma = model.Matching(matching_equations)

   fmu_path = model.GenerateFMU(
       model_path=model_path,
       Gamma=Gamma,
       res_eq=res_eq,
       fmu_name="three_tank_example",
   )

   print(fmu_path)

Induction Motor
---------------

This example imports the induction motor model and generates
``generated/induction_motor_example.fmu``. It uses MTES index ``1``.

.. code-block:: python

   import inspect
   from pathlib import Path

   from faultdiagnosistoolbox.models import _induction_motor_model


   model = _induction_motor_model.model
   model_path = inspect.getfile(_induction_motor_model)

   mtes = model.MTES()
   # Select the second MTES. Python list indices start at zero.
   selected_mtes = mtes[1]
   # Select the first equation in that MTES as the residual equation.
   res_eq = selected_mtes[0]
   matching_equations = [eq for eq in selected_mtes if eq != res_eq]
   Gamma = model.Matching(matching_equations)

   fmu_path = model.GenerateFMU(
       model_path=Path(model_path),
       Gamma=Gamma,
       res_eq=res_eq,
       fmu_name="induction_motor_example",
   )

   print(fmu_path)
