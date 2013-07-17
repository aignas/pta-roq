MAKEDIR_P = mkdir -p

#LIBS = -lgsl -L/usr/global/intel/mkl/10.2.6.038/lib/em64t \
	   -lmkl_lapack -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
LIBS = -lgsl -L/usr/lib64/atlas/ -lcblas -lclapack -latlas
COMPILER = g++ -fopenmp
CFLAGS = -g -std=c++11 -pedantic -Wall -O3

# Some directory variables
DIR_SRC  = src
DIR_OBJ  = obj
DIR_RUN  = data
DIR_SRC_T= src_test
DIR_OBJ_T= obj_test
DIR_RUN_T= data_test

DIR_ALL = $(DIR_SRC) $(DIR_OBJ) $(DIR_RUN) $(DIR_SRC_T) $(DIR_OBJ_T) $(DIR_RUN_T)

.PHONY: compile debug_compile clean

# Manually write down the headers and the main program files
SRC_LIB_HH    := $(shell find $(DIR_SRC) -name *.hh)
SRC_LIB_CC    := $(patsubst $(DIR_SRC)/%.hh, %.cc, $(SRC_LIB_HH))

SRC_TARGET_CC := data_generator.cc roq_generator.cc
SRC_TARGET_OBJ:= $(patsubst %.cc, $(DIR_OBJ)/%.o, $(SRC_TARGET_CC))
SRC_TARGET    := $(patsubst %.cc, %, $(SRC_TARGET_CC))

SRC_OBJECTS   := $(patsubst %.cc, $(DIR_OBJ)/%.o, $(SRC_LIB_CC))

# For test files as well (include the main headers here as well
TEST_LIB_HH   := linalg.hh roq.hh
TEST_LIB_CC   := $(patsubst %.hh, %.cc, $(TEST_LIB_HH))

TEST_TARGET_CC:= test_runner.cc
TEST_TARGET_OBJ:=$(patsubst %.cc, $(DIR_OBJ_T)/%.o, $(TEST_TARGET_CC))
TEST_TARGET   := $(patsubst %.cc, %, $(TEST_TARGET_CC))

TEST_OBJECTS  := $(patsubst %.cc, $(DIR_OBJ)/%.o, $(SRC_LIB_CC)) \
                 $(patsubst %.cc, $(DIR_OBJ_T)/%.o, $(TEST_LIB_CC))

# Append directories to the headers
SRC_LIB_HH    := $(patsubst %.hh, $(DIR_SRC)/%.hh, $(SRC_LIB_HH))
TEST_LIB_HH   := $(patsubst %.hh, $(DIR_SRC_T)/%.hh, $(TEST_LIB_HH)) \
                 $(patsubst %.hh, $(DIR_SRC)/%.hh, $(SRC_LIB_HH))

# All CC files to compile
ALL_CC        := $(DIR_SRC)/%.cc $(DIR_SRC_T)/%.cc

# Require directory creation
$(DIR_ALL):
	$(MAKEDIR_P) $@

$(SRC_OBJECTS):     | $(DIR_OBJ)
$(SRC_TARGET_OBJ):  | $(DIR_OBJ)
$(TEST_OBJECTS):    | $(DIR_OBJ) $(DIR_OBJ_T) $(DIR_RUN)
$(TEST_TARGET_OBJ): | $(DIR_OBJ) $(DIR_OBJ_T) $(DIR_RUN_T)

$(DIR_OBJ)/%.o: $(DIR_SRC)/%.cc
	@echo Building $@
	@$(COMPILER) -c $(CFLAGS) $< -o $@

$(DIR_OBJ_T)/%.o: $(DIR_SRC_T)/%.cc
	@echo Building $@
	@$(COMPILER) -c $(CFLAGS) $< -o $@

$(SRC_TARGET): $(SRC_TARGET_OBJ)
	$(COMPILER) $(SRC_OBJECTS) $(DIR_OBJ)/$@.o -Wall $(LIBS) -o $@

$(TEST_TARGET): $(TEST_TARGET_OBJ)
	$(COMPILER) $(TEST_OBJECTS) $(DIR_OBJ_T)/$@.o -Wall $(LIBS) -o $@

.PRECIOUS: $(SRC_TARGET) $(SRC_OBJECTS) $(TEST_TARGET) $(TEST_OBJECTS)

compile: $(SRC_TARGET)
compile_test: clean $(TEST_TARGET)

default: compile
all: default

debug_compile: cleanall compile compile_test

test: compile_test
	@echo Executing the test suite
	./test_runner

roq: clean compile
	@echo Running the Reduced Basis generator
	./roq_generator

cleanall:
	-rm -rf $(DIR_OBJ) $(DIR_OBJ_T) $(SRC_TARGET) $(TEST_TARGET)

clean:
	-rm -rf $(SRC_TARGET) $(TEST_TARGET)
