Neural Network Integration Details
==================================

Overview
--------
This extension bridges the gap between structural analysis and data-driven residual generation by leveraging **Neural Ordinary Differential Equations (Neural ODEs)**. Instead of manual configuration, it transforms a selected structural matching (``gamma``) and residual equation (``res_eq``) into optimized neural model configurations. 

By utilizing the structural information from the model, the extension ensures that the generated neural residuals are grounded in the system's underlying physical dependencies, even when the exact dynamics are learned from data.

Structural Extraction and Transformation
----------------------------------------
When the ``GenerateConfigFile`` method is invoked, the extension transforms the results of a previously conducted structural analysis into a formal model specification. The process relies on the following inputs:

* **Pre-computed Matching**: The method takes a specific matching ``gamma`` and the designated residual equation ``res_eq``, which have been identified through prior causality analysis (via a causality sweep).

* **Dependency Mapping**: Using the structural information embedded in the model, the extension traces dependencies to categorize variables into:

    - **State variables**: Internal dynamics required for the model.
    - **Input/Sensor variables**: Drivers identified from the known variable set (``Z``).
    - **Prediction Targets**: The specific signals the neural network is tasked to reconstruct.

* **Parameter Integration**: Global settings and signal mappings are injected via the ``params`` dictionary to ensure the final YAML is self-contained. This dictionary acts as the bridge between the structural model and the data files by defining:

    - ``x_var``: The symbolic variable names used within the structural model.
    - ``y_var``: The corresponding public names used for external representation.
    - ``signals``: The definitive mapping that connects these model variables to the specific column headers in the CSV data files.

  The resulting YAML file includes this ``signals`` dictionary as a lookup table, ensuring the Neural Residual Toolbox correctly identifies which data column corresponds to which physical input, state, or measurement identified during the structural extraction.

Structural vs. Latent Modeling
------------------------------
The extension supports two different modeling philosophies, defined by the ``model_type`` parameter:

* **Latent Mode**: Generates a configuration utilizing a hidden state-space container (``num_latents``). In this mode, the structural analysis is used to identify the necessary inputs and outputs, while the internal physical state equations are replaced by a neural ODE. This is ideal for black-box or partially known systems.

* **Greybox Mode**: Preserves the explicit structural state definitions and relationships defined in the provided matching (``gamma``). This creates a neural network architecture constrained by the known physical structure of the model, ensuring that the learned residuals respect the system's structural constraints.

Integration with Training Pipelines
-----------------------------------
The generated configurations are natively compatible with the PyTorch Lightning-based training pipeline in the Neural Residual Toolbox. They define the dynamic equations, predictors, and dataset descriptions required to train robust, physics-informed neural residuals for fault diagnosis.