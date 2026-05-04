Toolbox Extension
=================
This section describes three user accessibility features of the toolbox.

Model transformation
--------------------


External simulation
--------------------


The external simulation subsystem was implemented to enable residual computations from the Fault Diagnosis Toolbox to be exported as an FMU. The functionality was integrated into the central ``DiagnosisModel`` class, allowing the user to generate an FMU through a single method call::

    model.GenerateFMU(model_path, Gamma, res_eq)

Read more:

.. toctree::
   :maxdepth: 1

   Use cases </toolbox_extensions/externalsimulation/external_simulation_usecase>
   FMU generation details </toolbox_extensions/externalsimulation/external_simulation_details>
   FMU import in Simulink </toolbox_extensions/externalsimulation/run_fmu_in_matlab>


Neural network integration
---------------------------
  
