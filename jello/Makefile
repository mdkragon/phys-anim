# Makefile used to build jello sim

TARGETS = all clean
.PHONY: $(TARGETS)

CXX=g++
CC=gcc
LD=ld
STRIP=strip
CXX_FLAGS= -O3 -fpic


PROJ_NAME = jello
OUT_DIR = ./output


SRC_FILES = $(wildcard *.cpp) $(wildcard tinyxml/*.cpp)
OBJ_FILES = $(patsubst %.cpp, %.o,$(SRC_FILES))

INC_DIRS = -I./tinyxml -I/usr/local/include

LIB_DIRS = -L/usr/local/lib -L/usr/lib
LD_FLAGS += -lGL -lGLU -lglut -lILU -lILUT

# for tinyxml stl support
CXX_FLAGS += -DTIXML_USE_STL



ifndef OSTYPE
  OSTYPE = $(shell uname -s|awk '{print tolower($$0)}')
endif

ifeq ($(OSTYPE),linux)
	CXX_FLAGS+= -DLINUX
  SHLIBEXT= so
  LIBOPTS= -shared -fpic
  LIBRT= -lrt
endif
ifeq ($(OSTYPE),darwin)
  SHLIBEXT= dylib
  LIBOPTS= -bundle -undefined dynamic_lookup
	LD_FLAGS= -framework Carbon -framework OpenGL -framework GLUT -lIL -lILU -lILUT 
ifeq ($(MODE),32)
  CC=gcc -arch i386
  CXX=g++ -arch i386
  LD=g++ -arch i386
endif
  CXX_FLAGS+= -O2
  LIBRT=
endif

all: jello run_jello

%.o: %.cc
	$(CXX) $(CXX_FLAGS) $(INC_DIRS) -o $@ -c $<
%.o: %.cpp
	$(CXX) $(CXX_FLAGS) $(INC_DIRS) -o $@ -c $<
%.o: %.c
	$(CC) $(CXX_FLAGS) $(INC_DIRS) -o $@ -c $<

jello: $(OBJ_FILES)
	@mkdir -p $(OUT_DIR)
	$(CXX) -o $(OUT_DIR)/$@ $^ $(CXX_FLAGS) $(INC_DIRS) $(LD_FLAGS) $(LIB_DIRS)

run_jello:
	$(OUT_DIR)/jello
	

clean:
	rm -f *.o tinyxml/*.o $(OUT_DIR)/jello

