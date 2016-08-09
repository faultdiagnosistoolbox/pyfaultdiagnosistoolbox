SOURCE = ResGen_1649_1.cc ResGen_4017_84.cc ResGen_4477_86.cc ResGen_4011_1.cc ResGen_4066_1.cc ResGen_4016_1.cc ResGen_4074_84.cc

BUILD_DIR = build/lib.macosx-10.6-x86_64-3.5

OPT=-DNDEBUG -O3
CPPFLAGS=-Wall $(OPT)
TARGET_VERSION=`sw_vers -productVersion` # 11.10

all:
	@echo make clean/cleanall to clear all generated files

#$(BUILD_DIR)/ResGen_1649_1c_b.so: ResGen_1649_1c_b.cc
#	CC=g++ MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python $(patsubst %.cc, %_setup.py, $<) build

clean:
	rm -rf $(BUILD_DIR)/*.so extmodelfuns.o *~

cleanall:
	make clean
	rm -rf $(patsubst %.cc, %_setup.py, $(SOURCE)) $(SOURCE)


