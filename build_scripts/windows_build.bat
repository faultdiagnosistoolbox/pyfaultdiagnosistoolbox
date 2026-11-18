@echo off
setlocal

REM Python versions to build for (minor versions of Python 3)
set VERSIONS=12 13 14 15

REM Create build venvs if necessary
if not exist "build_envs" (
    mkdir build_envs

    for %%v in (%VERSIONS%) do (
        py -3.%%v -m venv build_envs\env3%%v
        call build_envs\env3%%v\Scripts\activate.bat
        python -m pip install -U pip build
        call deactivate
    )
)

REM Build wheels
for %%v in (%VERSIONS%) do (
    echo === Building wheel for Python 3.%%v ===
    call build_envs\env3%%v\Scripts\activate.bat
    python -m build --wheel
    call deactivate
)

REM Cleanup
REM rmdir /s /q build_envs

endlocal
