Neural Network Integration Use Cases
====================================

This guide demonstrates how to generate YAML configurations for the Neural Residual Toolbox 
using the NeuralNetworkIntegration extension. By leveraging structural analysis, the 
extension automates the creation of physics-informed neural network architectures.

The workflow is:

1. **Define Model and Signal Mappings**
   Import or define a ``DiagnosisModel`` and create a ``params`` dictionary. This dictionary 
   is crucial as it maps the symbolic model variables (``x_var``) to the dataset column 
   names (``signals`` and ``y_var``).

2. **Select MSO and Perform Causality Analysis**
   Select a Minimal Structurally Overdetermined (MSO) set and perform a causality analysis 
   (e.g., via a causality sweep) to determine a valid matching (``gamma``) and a 
   suitable residual equation (``res_eq``).

3. **Choose Modeling Philosophy**
   Decide on the ``model_type``:
   
   - Use ``'Latent'`` for a black-box neural ODE approach.
   - Use ``'Greybox'`` to enforce the explicit structural dependencies identified in the matching.

4. **Generate the Configuration File**
   Call ``GenerateConfigFile`` with the target path, matching, residual equation, and 
   parameters. The resulting YAML file is natively compatible with the toolbox training pipeline.

Prerequisites
-------------

To use this extension effectively, ensure you have:

- A valid ``DiagnosisModel`` object.
- A selected Matching (``gamma``) and a residual equation index (``res_eq``).
- A defined ``params`` dictionary containing the signal-to-variable mappings.

Example Scenarios
=================

This section provides practical examples of how to use the integration extension to generate different types of neural residual configurations.

Example 1 – Latent Neural ODE
-----------------------------

.. code-block:: python

   import numpy as np
   import os
   from tests.test_neural_network.test_engine_model import LiU_ICE_model
   from faultdiagnosistoolbox.NeuralNetworkIntegration.helpers import GenConfigParams

   # 1. Setup Model and Parameters
   model = LiU_ICE_model
   mso = model.MSO()[584] # MSO 585
   res_eq = mso[3]

   signals = ["intercooler_pressure", "intercooler_temperature", "intake_manifold_pressure", 
               "air_mass_flow", "engine_speed", "throttle_position", "wastegate_position", 
               "injected_fuel_mass", "ambient_temperature", "ambient_pressure"]
   y_vars = ["y_p_ic", "y_T_ic", "y_p_im", "y_W_af", "y_omega_e", "y_alpha_th", "y_u_wg", "y_wfc", "y_T_amb", "y_p_amb"]
   x_vars = ["p_ic", "T_ic", "p_im", "W_af", "omega_e", "alpha_th", "u_wg", "wfc", "T_amb", "p_amb"]
   
   params = GenConfigParams(y_vars, x_vars, signals)

   # 2. Perform Structural Analysis
   gamma = model.Matching(np.setdiff1d(mso, res_eq))

   # 3. Generate Configuration
   output_path = "output_configs"
   if not os.path.exists(output_path): os.makedirs(output_path)

   model.GenerateConfigFile(
      output_path, gamma, res_eq, params, 
      file_name="example_latent", model_type="latent"
   )

**Resulting latent configuration file:**

.. code-block:: yaml

   # General Description
   description: 'Generated latent configuration file: example_latent for model: LiU_ICE model'

   dataset: LiU_ICE model

   signals: {
      intercooler_pressure: 'y_p_ic', 
      intercooler_temperature: 'y_T_ic', 
      intake_manifold_pressure: 'y_p_im', 
      air_mass_flow: 'y_W_af', 
      engine_speed: 'y_omega_e', 
      throttle_position: 'y_alpha_th', 
      wastegate_position: 'y_u_wg', 
      injected_fuel_mass: 'y_wfc', 
      ambient_temperature: 'y_T_amb', 
      ambient_pressure: 'y_p_amb'
   }

   # Signals to be zeroed
   zeroed_signals: null

   # Dynamic Equations
   dynamic: {
      latent: {
         states: [ latent ], 
         num_latents: 13, 
         inputs: [ y_T_amb, y_W_af, y_alpha_th, y_omega_e, y_p_amb, y_p_ic, y_p_im, y_u_wg, y_wfc ]
      }
   }

   # Predictors
   predictors: {
      y_T_ic: {
         states: [ latent ], 
         use_latent: True, 
         inputs: [ y_p_ic, y_p_im, y_alpha_th, y_T_amb ]
      }
   }


Example 2 - Greybox Neural ODE
-------------------------------

The Greybox approach enforces the specific structural dependencies identified in the matching. Unlike the Latent approach, this creates a physics-informed neural network where states and dependencies are explicitly defined based on the model equations.

.. code-block:: python

   import numpy as np
   import os
   from tests.test_neural_network.test_engine_model import LiU_ICE_model
   from faultdiagnosistoolbox.NeuralNetworkIntegration.helpers import GenConfigParams

   # Setup Model and Parameters
   model = LiU_ICE_model
   mso = model.MSO()[584] # MSO 585
   res_eq = mso[3] # Using an equation known to have a valid causality

   # Mapping symbolic variables to dataset signals
   signals = ["intercooler_pressure", "intercooler_temperature", "intake_manifold_pressure", 
              "air_mass_flow", "engine_speed", "throttle_position", "wastegate_position", 
              "injected_fuel_mass", "ambient_temperature", "ambient_pressure"]
   y_vars = ["y_p_ic", "y_T_ic", "y_p_im", "y_W_af", "y_omega_e", "y_alpha_th", "y_u_wg", "y_wfc", "y_T_amb", "y_p_amb"]
   x_vars = ["p_ic", "T_ic", "p_im", "W_af", "omega_e", "alpha_th", "u_wg", "wfc", "T_amb", "p_amb"]
    
   params = GenConfigParams(y_vars, x_vars, signals)

   # Structural Analysis
   # The matching determines the internal computational flow
   gamma = model.Matching(np.setdiff1d(mso, res_eq))

   # Generate Greybox Configuration
   output_path = "output_configs"
   if not os.path.exists(output_path): os.makedirs(output_path)

   model.GenerateConfigFile(
      path=output_path, 
      gamma=gamma, 
      res_eq=res_eq, 
      params=params, 
      file_name="example_greybox", 
      model_type="greybox"
   )

**Resulting greybox configuration file:**

.. code-block:: yaml

   # General Description
   description: 'Generated greybox configuration file: example_greybox for model: LiU_ICE model'

   dataset: LiU_ICE model

   signals: {
      intercooler_pressure: 'y_p_ic', 
      intercooler_temperature: 'y_T_ic', 
      intake_manifold_pressure: 'y_p_im', 
      air_mass_flow: 'y_W_af', 
      engine_speed: 'y_omega_e', 
      throttle_position: 'y_alpha_th', 
      wastegate_position: 'y_u_wg', 
      injected_fuel_mass: 'y_wfc', 
      ambient_temperature: 'y_T_amb', 
      ambient_pressure: 'y_p_amb'
   }

   # Signals to be zeroed
   zeroed_signals: null

   # Dynamic Equations
   dynamic: {
      m_t: {
         states: [ m_em, T_em, m_t, T_t, wg_pos ], 
         inputs: [ y_p_amb ]
      }, 
      T_t: {
         states: [ m_em, T_em, m_t, T_t, omega_tc, wg_pos ], 
         inputs: [ y_p_amb ]
      }, 
      m_c: {
         states: [ m_af, T_af, m_c, T_c, omega_tc ], 
         inputs: [ y_p_ic ]
      }, 
      T_c: {
         states: [ m_af, T_af, m_c, T_c, omega_tc ], 
         inputs: [ y_p_ic ]
      }, 
      m_ic: {
         states: [ m_c, T_c, T_ic ], 
         inputs: [ y_p_ic, y_p_im, y_alpha_th ]
      }, 
      T_ic: {
         states: [ m_c, T_c, m_ic, T_ic ], 
         inputs: [ y_p_ic, y_p_im, y_alpha_th, y_T_amb ]
      }, 
      m_im: {
         states: [ m_im, m_em, T_em, T_ic ], 
         inputs: [ y_p_ic, y_p_im, y_omega_e, y_alpha_th ]
      }, 
      m_em: {
         states: [ m_im, m_em, T_em, m_t, T_t, wg_pos ], 
         inputs: [ y_p_im, y_omega_e, y_wfc ]
      }, 
      T_em: {
         states: [ m_im, m_em, T_em, m_t, T_t, wg_pos ], 
         inputs: [ y_p_im, y_omega_e, y_wfc, y_T_amb ]
      }, 
      m_af: {
         states: [ m_af, T_af, m_c, T_c, omega_tc ], 
         inputs: [ y_W_af ]
      }, 
      T_af: {
         states: [ m_af, T_af, m_c, T_c, omega_tc ], 
         inputs: [ y_W_af, y_T_amb ]
      }, 
      omega_tc: {
         states: [ m_af, T_af, m_c, T_c, m_em, T_em, m_t, T_t, omega_tc ], 
         inputs: [ ]
      }, 
      wg_pos: {
         states: [ wg_pos ], 
         inputs: [ y_u_wg ]
      }
   }

   # Predictors
   predictors: {
      y_T_ic: {
         states: [ m_c, T_c, m_ic, T_ic ], 
         inputs: [ y_p_ic, y_p_im, y_alpha_th, y_T_amb ]
      }
   }

Example 3 – Multiple residuals from a single MSO
---------------------------------------------

Selecting different residual equations from the same MSO yields alternative architectures for the same fault. In MSO 1590, both equations 85 and 86 provide valid causal matchings.

.. code-block:: python

   import numpy as np
      import os
      from tests.test_neural_network.test_engine_model import LiU_ICE_model
      from faultdiagnosistoolbox.NeuralNetworkIntegration.helpers import GenConfigParams

      # Setup Model and select a flexible MSO
      model = LiU_ICE_model
      mso = model.MSO()[1589] # MSO 1590
      
      signals = ["intercooler_pressure", "intercooler_temperature", "intake_manifold_pressure", 
                  "air_mass_flow", "engine_speed", "throttle_position", "wastegate_position", 
                  "injected_fuel_mass", "ambient_temperature", "ambient_pressure"]
      y_vars = ["y_p_ic", "y_T_ic", "y_p_im", "y_W_af", "y_omega_e", "y_alpha_th", "y_u_wg", "y_wfc", "y_T_amb", "y_p_amb"]
      x_vars = ["p_ic", "T_ic", "p_im", "W_af", "omega_e", "alpha_th", "u_wg", "wfc", "T_amb", "p_amb"]
      params = GenConfigParams(y_vars, x_vars, signals)

      # Generate multiple configurations for different residual equations
      # In MSO 1590, equations 85 and 86 are both causally valid
      output_path = "output_configs/multi"
      if not os.path.exists(output_path): os.makedirs(output_path)

      for eq_num in [85, 86]:
         gamma = model.Matching(np.setdiff1d(mso, eq_num))
         model.GenerateConfigFile(
               output_path, gamma, eq_num, params, 
               file_name=f"mso1590_res_eq_{eq_num}", 
               model_type="latent"
         )
         print(f"Generated config for equation {eq_num}")

First configurationfile (res_eq 85):
------------------------------------

.. code-block:: yaml

   # General Description
   description: 'Generated latent configuration file: mso1590_res_eq_85 for model: LiU_ICE model'

   dataset: LiU_ICE model

   signals: {
      intercooler_pressure: 'y_p_ic', 
      intercooler_temperature: 'y_T_ic', 
      intake_manifold_pressure: 'y_p_im', 
      air_mass_flow: 'y_W_af', 
      engine_speed: 'y_omega_e', 
      throttle_position: 'y_alpha_th', 
      wastegate_position: 'y_u_wg', 
      injected_fuel_mass: 'y_wfc', 
      ambient_temperature: 'y_T_amb', 
      ambient_pressure: 'y_p_amb'
   }

   # Signals to be zeroed
   zeroed_signals: null

   # Dynamic Equations
   dynamic: {
      latent: {
         states: [ latent ], 
         num_latents: 11, 
         inputs: [ y_T_amb, y_W_af, y_alpha_th, y_omega_e, y_p_amb, y_p_ic, y_p_im, y_u_wg, y_wfc ]
      }
   }

   # Predictors
   predictors: {
      y_T_ic: {
         states: [ latent ], 
         use_latent: True, 
         inputs: [ y_p_ic ]
      }
   }

Second configuration file (res_eq 86):
------------------------------------

.. code-block:: yaml

   # General Description
   description: 'Generated latent configuration file: mso1590_res_eq_86 for model: LiU_ICE model'

   dataset: LiU_ICE model

   signals: {
      intercooler_pressure: 'y_p_ic', 
      intercooler_temperature: 'y_T_ic', 
      intake_manifold_pressure: 'y_p_im', 
      air_mass_flow: 'y_W_af', 
      engine_speed: 'y_omega_e', 
      throttle_position: 'y_alpha_th', 
      wastegate_position: 'y_u_wg', 
      injected_fuel_mass: 'y_wfc', 
      ambient_temperature: 'y_T_amb', 
      ambient_pressure: 'y_p_amb'
   }

   # Signals to be zeroed
   zeroed_signals: null

   # Dynamic Equations
   dynamic: {
      latent: {
         states: [ latent ], 
         num_latents: 11, 
         inputs: [ y_T_amb, y_T_ic, y_W_af, y_alpha_th, y_omega_e, y_p_amb, y_p_ic, y_u_wg, y_wfc ]
      }
   }

   # Predictors
   predictors: {
      y_p_im: {
         states: [ latent ], 
         use_latent: True, 
         inputs: [ y_p_ic, y_T_ic, y_alpha_th ]
      }
   }