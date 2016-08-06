SOURCE = ResGen_1649_1c.cc ResGen_1649_1c_b.cc
BUILD_DIR = build/lib.macosx-10.6-x86_64-2.7

TARGET_DIR = $(patsubst %.cc, $(BUILD_DIR)/%.so, $(SOURCE))

OPT=-DNDEBUG -O3
CPPFLAGS=-Wall $(OPT)
TARGET_VERSION=`sw_vers -productVersion` # 11.10

all: $(TARGET_DIR)

$(BUILD_DIR)/ResGen_1649_1c.so: ResGen_1649_1c.cc
	CC=g++ MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python $(patsubst %.cc, %_setup.py, $<) build

$(BUILD_DIR)/ResGen_1649_1c_b.so: ResGen_1649_1c_b.cc
	CC=g++ MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python $(patsubst %.cc, %_setup.py, $<) build

clean:
	rm -rf $(TARGET_DIR) extmodelfuns.o *~

cleanall:
	make clean
	rm -f $(SOURCE) $(patsubst %.c, %_setup.py, $(SOURCE))



