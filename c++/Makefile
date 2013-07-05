MAKEDIR_P = mkdir -p

TARGET = data_generator
LIBS = -lgsl -L/usr/global/intel/mkl/10.2.6.038/lib/em64t \
	   -lmkl_lapack -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
COMPILER = g++
CFLAGS = -g -pedantic -Wall -O3
SRCDIR = src
OBJDIR = obj
TESTDIR = test

.PHONY: default all clean

default: $(OBJDIR) $(TARGET)
all: clean default
runall: default

# Manually write down the headers and the main program files
SCR_HH  := linalg.hh pulsar.hh signal-model.hh signal-vectors.hh signal-sampling.hh random-helper.hh
SCR_CC  := $(patsubst %.hh, %.cc, $(SCR_HH)) data_generator.cc
TEST_HH := test-linalg.hh test-pulsar.hh test-signal-model.hh test-signal-vectors.hh
TEST_CC := $(patsubst %.hh, %.cc, $(TEST_HH)) test_runner.cc

OBJECTS := $(patsubst %.cc, $(OBJDIR)/%.o, $(SCR_CC))

$(OBJECTS): | $(OBJDIR)

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS) | $(RUNDIR)
	$(COMPILER) $(OBJECTS) -Wall $(LIBS) -o $@

$(OBJDIR):
	$(MAKEDIR_P) $@

$(RUNDIR):
	$(MAKEDIR_P) $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(HEADERS)
	@echo Building $@
	@$(COMPILER) -c $(CFLAGS) $< -o $@

test: default
	./$(TARGET) test

clean:
	-rm -rf $(OBJDIR) $(TARGET)

compare: default
	time ./$(TARGET)
	time python2 ../py2.7/data-gen.py

# DO NOT DELETE


