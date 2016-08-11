TARGET_VERSION=`sw_vers -productVersion` # 11.10
.PHONY: module clean

module: structuralanalysismodule.cc setup.py SparseMatrix.cc StructuralAnalysisModel.cc MSOAlg.cc
	MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python setup.py build

clean:
	rm -rf build/* *~ 

