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




Neural Network Integration
--------------------------


The Neural Network Integration was created to allow automated generation of configuration files for the `Neural Residual Toolbox <https://github.com/westny/neural-residual>`_. The functionality was integrated into the central ``DiagnosisModel`` class, allowing the user to generate configuration files through a single method call::
   
    model.generate_config_file(gamma, path, res_eq, params, file_name, type)


Read more:


.. toctree::
   :maxdepth: 1


   Neural Network Integration details <toolbox_extensions/neural_network_integration/neural_network_integration_details>
   Use cases <toolbox_extensions/neural_network_integration/neural_network_integration_usecases>
