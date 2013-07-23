MAKEDIR_P = mkdir -p

#LIBS = -lgsl -L/usr/global/intel/mkl/10.2.6.038/lib/em64t \
-lmkl_lapack -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
LIBS = -lgsl -L/usr/lib64/atlas/ -lcblas -lclapack -latlas
COMPILER = g++ -fopenmp
CFLAGS = -g -std=c++11 -pedantic -Wall -O3

# Some directory variables
DIR_SRC = src
DIR_RUN = data
DIR_TARGET = bin
DIR_OBJ = obj
DIR_OBJ_ALL = $(patsubst $(DIR_SRC)/%, $(DIR_OBJ)/%, $(shell find $(DIR_SRC) -type d | tail -n +2))

# Find all libs and Targets
LIB_OBJ    := $(patsubst $(DIR_SRC)/%.hh, $(DIR_OBJ)/%.o, $(shell find $(DIR_SRC) -name *.hh))
TARGET     := $(patsubst $(DIR_SRC)/%.cc, $(DIR_TARGET)/%, $(shell find $(DIR_SRC) -maxdepth 1 -name *.cc))
TARGET_OBJ := $(patsubst $(DIR_TARGET)/%, $(DIR_OBJ)/%.o, $(TARGET))

compile: $(LIB_OBJ) $(TARGET_OBJ) $(TARGET)
compile_test: compile test

all: clean compile_test
default: compile

.PHONY: clean
.PRECIOUS: $(LIB_OBJ) $(TARGET_OBJ)

# Build a directory tree for the obj directory
$(DIR_SRC) $(DIR_RUN) $(DIR_TARGET) $(DIR_OBJ) $(DIR_OBJ_ALL):
	$(MAKEDIR_P) $@

$(DIR_OBJ)/%.o: $(DIR_SRC)/%.cc | $(DIR_OBJ) $(DIR_OBJ_ALL) $(DIR_TARGET)
	@echo Building $@
	@$(COMPILER) -c $(CFLAGS) $^ -o $@

$(DIR_TARGET)/%: $(DIR_OBJ)/%.o $(LIB_OBJ) | $(DIR_RUN) $(TARGET_OBJ) $(LIB_OBJ)
	@echo Linking $@
	@$(COMPILER) -Wall $^ $(LIBS) -o $@

test: compile
	@bin/./test_runner

clean:
	-rm -rf $(DIR_OBJ) $(DIR_TARGET)
