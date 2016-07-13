TARGET=build/lib.macosx-10.6-x86_64-2.7/mypyc.so
CSPARSE=/Users/frisk/sw/CSParse
OPT=-DNDEBUG -O3
CPPFLAGS=-Wall $(OPT) -I$(CSPARSE)/Include 

$(TARGET): structuralanalysismodule.c setup.py SparseMatrix.o
	MACOSX_DEPLOYMENT_TARGET=10.11 python setup.py build

SparseMatrix.o: SparseMatrix.cc SparseMatrix.h
	g++ -c $< -o $@ $(CPPFLAGS)

clean:
	rm -rf build *~ SparseMatrix.o

