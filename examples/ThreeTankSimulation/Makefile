SOURCE = ResGen1.c ResGen2.c ResGen3.c
BUILD_DIR = build/lib.macosx-10.6-x86_64-2.7

TARGET_DIR = $(patsubst %.c, $(BUILD_DIR)/%.so, $(SOURCE))

OPT=-DNDEBUG -O3
CPPFLAGS=-Wall $(OPT)
TARGET_VERSION=`sw_vers -productVersion` # 11.10

all: $(TARGET_DIR)

$(BUILD_DIR)/ResGen1.so: ResGen1.c
	MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python $(patsubst %.c, %_setup.py, $<) build

$(BUILD_DIR)/ResGen2.so: ResGen2.c
	MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python $(patsubst %.c, %_setup.py, $<) build

$(BUILD_DIR)/ResGen3.so: ResGen3.c
	MACOSX_DEPLOYMENT_TARGET=$(TARGET_VERSION) python $(patsubst %.c, %_setup.py, $<) build

clean:
	rm -rf build *~

cleanall:
	make clean
	rm $(SOURCE)


