TARGET=build/lib.macosx-10.6-x86_64-2.7/structuralanalysis3.so
CSPARSE=/Users/frisk/sw/CSParse
OPT=-DNDEBUG -O3
CPPFLAGS=-Wall $(OPT) -I$(CSPARSE)/Include 
TARGET_VERSION=`sw_vers -productVersion` # 11.10

$(TARGET): structuralanalysis3module.cc setup.py SparseMatrix.o StructuralAnalysisModel.o MSOalg.o
	CC=g++ MACOSX_DEPLOYMENT_TARGET=10.11.5 python setup3.py build

SparseMatrix.o: SparseMatrix.cc SparseMatrix.h
	g++ -c $< -o $@ $(CPPFLAGS)

MSOalg.o: MSOalg.cc MSOalg.h
	g++ -c $< -o $@ $(CPPFLAGS)

StructuralAnalysisModel.o: StructuralAnalysisModel.cc StructuralAnalysisModel.h
	g++ -c $< -o $@ $(CPPFLAGS)

clean:
	rm -rf build *~ SparseMatrix.o StructuralAnalysisModel.o MSOalg.o

