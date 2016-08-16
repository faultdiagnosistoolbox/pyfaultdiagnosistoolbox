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
# Anaconda python sets MACOSX_DEPLOYMENT_TARGET to 10.6, so to avoid
# warnings, this is needed also for supporting libraries
# Should not affect Linux builds
	cd CSparse && MACOSX_DEPLOYMENT_TARGET=10.6 make 

clean:
	rm -rf build *~ faultdiagnosistoolbox/*.so

cleanall:
	make clean
	cd CSparse && make purge
	rm -rf faultdiagnosistoolbox.egg-info
	rm -rf dist

