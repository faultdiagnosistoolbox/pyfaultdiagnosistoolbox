.PHONY: build build_ext wheel sdist wheel clean cleanall alldist

build_ext: CSparse/Lib/libcsparse.a
	python setup.py build_ext --inplace

alldist: build wheel sdist

build: CSparse/Lib/libcsparse.a
	python setup.py build

wheel: CSparse/Lib/libcsparse.a
	python setup.py bdist_wheel

sdist:
	python setup.py sdist

CSparse/Lib/libcsparse.a:
	cd CSparse && make library

CSparse/Lib/csparse.lib:
	cd CSparse/Lib && nmake -f Makefile.win

clean:
	rm -rf build *~ faultdiagnosistoolbox/*.so

cleanall:
	make clean
	cd CSparse && make purge
	rm -rf faultdiagnosistoolbox.egg-info
	rm -rf dist

test:
	python -m pytest -v
