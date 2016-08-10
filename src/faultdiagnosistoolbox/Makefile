CSPARSE=/Users/frisk/sw/CSParse
OPT=-DNDEBUG -O3
CPPFLAGS=-Wall $(OPT) -I$(CSPARSE)/Include 
TARGET_VERSION=`sw_vers -productVersion` # 11.10
.PHONY: module clean

module: structuralanalysismodule.cc setup.py SparseMatrix.o StructuralAnalysisModel.o MSOAlg.o
	CC=g++ MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python setup.py build

SparseMatrix.o: SparseMatrix.cc SparseMatrix.h
	g++ -c $< -o $@ $(CPPFLAGS)

MSOAlg.o: MSOAlg.cc MSOAlg.h
	g++ -c $< -o $@ $(CPPFLAGS)

StructuralAnalysisModel.o: StructuralAnalysisModel.cc StructuralAnalysisModel.h
	g++ -c $< -o $@ $(CPPFLAGS)

clean:
	rm -rf build/* *~ SparseMatrix.o StructuralAnalysisModel.o MSOAlg.o

