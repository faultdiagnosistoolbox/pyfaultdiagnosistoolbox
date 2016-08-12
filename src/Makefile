.PHONY: package clean cleanall

package: CSparse/Lib/libcsparse.a
	python setup.py build

CSparse/Lib/libcsparse.a:
	cd CSparse && MACOSX_DEPLOYMENT_TARGET=10.6 make

clean:
	rm -rf build/* *~ 

cleanall:
	make clean
	cd CSparse && make purge


