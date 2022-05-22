.PHONY: build build_ext wheel sdist wheel clean cleanall alldist

build_ext: # CSparse/Lib/libcsparse.a  # Built by setup script
	python setup.py build_ext --inplace

wheel: # CSparse/Lib/libcsparse.a build_ext  # Built by setup script
	python -m build --wheel
#	python setup.py bdist_wheel

clean:
	rm -rf build *~ faultdiagnosistoolbox/*.so

cleanall:
	make clean
	cd CSparse && make purge
	rm -rf faultdiagnosistoolbox.egg-info

test:
	python -m pytest -v

upload:
	python -m twine upload dist/*

# Below not really needed
alldist: build wheel sdist

build: CSparse/Lib/libcsparse.a build_ext
	python setup.py build

sdist:
	python setup.py sdist

CSparse/Lib/libcsparse.a:
	cd CSparse && make library

CSparse/Lib/csparse.lib:
	cd CSparse/Lib && nmake -f Makefile.win

