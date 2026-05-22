# Building the Documentation locally

The documentation is built with Sphinx.

Run these commands from the repository root, not from the `docs/` directory.

You can build it from the same virtual environment used for toolbox development or create a dedicated one:

```bash
python -m venv venv
```

Activate it:
- **Linux / macOS**
  ```bash
  source venv/bin/activate
  ```

- **Windows (PowerShell)**
  ```powershell
  .\venv\Scripts\Activate.ps1
  ```

- **Windows (Command Prompt)**
  ```cmd
  venv\Scripts\activate.bat
  ```

Install the documentation dependencies:

```bash
pip install -r docs/requirements-doc.txt
```

Build the HTML documentation:

```bash
sphinx-build -b html docs/ docs/doc_build/html
```

Open the generated documentation in a browser (Linux):

```bash
xdg-open docs/doc_build/html/index.html
```

The `xdg-open` command is for Linux. On other systems, open `docs/doc_build/html/index.html` directly in your browser.
