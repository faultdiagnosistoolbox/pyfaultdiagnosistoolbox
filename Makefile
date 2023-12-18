.PHONY: test-cov

test:
	python -m pytest -v

test-cov:
	python -m pytest -v --cov=faultdiagnosistoolbox

test-cov-report:
	python -m pytest -v --cov=faultdiagnosistoolbox --cov-report=html

# usernmame __token__, password is the token
pypi_upload:
	twine upload src/dist/*