.PHONY: test-cov

test:
	python -m pytest -v

test-cov:
	python -m pytest -v --cov=faultdiagnosistoolbox

test-cov-report:
	python -m pytest -v --cov=faultdiagnosistoolbox --cov-report=html

build:
	python -m build --wheel

doc_build:
	sphinx-build -M html docs doc_build

# usernmame __token__, password is the token
pypi_upload:
	twine upload src/dist/*