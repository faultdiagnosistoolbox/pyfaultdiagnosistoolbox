Run FMU in MATLAB
=================

This guide describes the current workflow for running a generated FMU in MATLAB/Simulink with the toolbox Python environment.

1. Generate the FMU
-------------------

Make sure the FMU has been generated and exists in ``generated/``.

The example is called ``test.fmu`` in ``generated``.

.. code-block:: bash

   generated/test.fmu

2. Start MATLAB with the platform wrapper
-----------------------------------------

From the repository root, run:

Linux:

.. code-block:: bash

   ./build_scripts/start_matlab_with_preload.sh

Windows Command Prompt:

.. code-block:: cmd

   build_scripts\start_matlab_with_preload.bat

The Linux wrapper starts MATLAB with:

.. code-block:: bash

   LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libpython3.10.so"

The Windows wrapper starts MATLAB with the repository ``.venv`` added to ``PATH`` so MATLAB can find the virtual environment's Python DLLs.

3. Point MATLAB to the toolbox Python environment
-------------------------------------------------

In the MATLAB Command Window, run:

.. code-block:: matlab

   addpath(fullfile(pwd, "build_scripts"))
   setup_matlab_toolbox_env

This script:

* points MATLAB to the repo-local Python interpreter in ``.venv/bin/python`` on Linux or ``.venv\Scripts\python.exe`` on Windows
* validates the Python interpreter
* validates ``numpy``
* validates ``scipy``
* starts Simulink

4. Confirm the interpreter in MATLAB
------------------------------------

You should see output similar to:

.. code-block:: matlab

   /path/to/repo/.venv/bin/python
   2.x.x
   1.x.x

If MATLAB is already using another Python interpreter, restart MATLAB and run ``setup_matlab_toolbox_env`` before any other Python-related command.

5. Import the FMU in Simulink
-----------------------------

Once Simulink has opened:

1. Create a new model or open your existing model.
2. Insert an FMU block.
3. Browse to the FMU file.

Example FMU path:

.. code-block:: text

   generated/test.fmu

Or use the absolute path if needed.

6. Enable debug logging if the FMU fails
----------------------------------------

If the FMU block fails during import or simulation:

1. Open the FMU block parameters.
2. Enable FMU debug logging.
3. Run again and inspect the detailed error.

7. Useful manual checks in MATLAB
---------------------------------

If you want to verify the Python environment manually:

.. code-block:: matlab

   pyrun("import sys; print(sys.executable)")
   pyrun("import numpy; print(numpy.__version__)")
   pyrun("import scipy; print(scipy.__version__)")

8. Important note
-----------------

MATLAB must use the same Python environment as the toolbox FMU workflow.

The intended interpreter is:

.. code-block:: text

   <repo>/.venv/bin/python

on Linux, or:

.. code-block:: text

   <repo>\.venv\Scripts\python.exe

on Windows.

Not:

.. code-block:: text

   /usr/bin/python3.10

unless that system interpreter has been set up with the exact same compatible package versions.
