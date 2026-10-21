@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%.." >nul
set "REPO_ROOT=%CD%"
popd >nul

if not defined MATLAB_CMD set "MATLAB_CMD=matlab"

if exist "%REPO_ROOT%\.venv\Scripts" (
    set "PATH=%REPO_ROOT%\.venv\Scripts;%REPO_ROOT%\.venv;%PATH%"
)

call "%MATLAB_CMD%" %*
