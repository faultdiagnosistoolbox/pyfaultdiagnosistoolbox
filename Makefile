.PHONY: test-cov

test:
	python -m pytest -v

test-cov:
	python -m pytest -v --cov=faultdiagnosistoolbox

test-cov-report:
	python -m pytest -v --cov=faultdiagnosistoolbox --cov-report=html
