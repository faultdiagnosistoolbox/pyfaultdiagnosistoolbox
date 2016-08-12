#TARGET_VERSION=`sw_vers -productVersion` # 11.10
.PHONY: module clean

module: structuralanalysismodule.cc setup.py SparseMatrix.cc StructuralAnalysisModel.cc MSOAlg.cc
	python setup.py build_ext --inplace

clean:
	rm -rf build/* *~ structuralanalysis*.so

