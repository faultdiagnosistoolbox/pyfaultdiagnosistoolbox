TARGET_VERSION=`sw_vers -productVersion` # 11.10
.PHONY: package clean

package: 
	MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python setup.py build

clean:
	rm -rf build/* *~ 

