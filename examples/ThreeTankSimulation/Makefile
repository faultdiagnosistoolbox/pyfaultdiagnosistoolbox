SOURCE = ResGen1.c ResGen2.c ResGen3.c
#TARGET_DIR = $(patsubst %.cc, $(BUILD_DIR)/%.so, $(SOURCE))

OPT=-DNDEBUG -O3
CPPFLAGS=-Wall $(OPT)
#TARGET_VERSION=`sw_vers -productVersion` # 11.10

all:
	python $(patsubst %.cc, %_setup.py, $<) build_ext --inplace
	python $(patsubst %.cc, %_setup.py, $<) build_ext --inplace
	python $(patsubst %.cc, %_setup.py, $<) build_ext --inplace

clean:
	rm -rf *~ *.so

cleanall:
	make clean
	rm -f $(SOURCE) $(patsubst %.cc, %_setup.py, $(SOURCE))



