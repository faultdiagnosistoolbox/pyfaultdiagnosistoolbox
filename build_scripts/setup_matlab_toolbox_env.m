function setup_matlab_toolbox_env()
% Point MATLAB to the repo-local Python environment and validate it.

script_dir = fileparts(mfilename("fullpath"));
repo_root = fileparts(script_dir);

if ispc
    venv_python = fullfile(repo_root, ".venv", "Scripts", "python.exe");
else
    venv_python = fullfile(repo_root, ".venv", "bin", "python");
end

if ~isfile(venv_python)
    error("Python interpreter not found: %s", venv_python);
end

py = pyenv;
if py.Status == "NotLoaded"
    pyenv("Version", venv_python);
elseif string(py.Executable) ~= string(venv_python)
    warning([ ...
        "MATLAB Python is already loaded as %s. Restart MATLAB and run " ...
        "setup_matlab_toolbox_env first to switch to %s."], ...
        string(py.Executable), string(venv_python));
end

pyrun("import sys; print(sys.executable)")
pyrun("import numpy; print(numpy.__version__)")
pyrun("import scipy; print(scipy.__version__)")
simulink
end
