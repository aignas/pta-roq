MAKEDIR_P = mkdir -p

TARGET = program
LIBS = -lgsl -lgslcblas
CC = g++
CFLAGS = -g -pedantic -Wall -O2
SRCDIR = src
OBJDIR = obj
RUNDIR = data_raw/exp1 data_raw/exp2 data_raw/exp3 data_raw/exp4

.PHONY: default all clean

default: $(OBJDIR) $(TARGET)
all: clean default
runall: default

# Automatically include all cc files in the src directory
HEADERS  := $(wildcard $(SRCDIR)/*.hh) $(wildcard $(TESTDIR)/*.hh)
CC_FILES := $(wildcard $(SRCDIR)/*.cc) $(wildcard $(TESTDIR)/*.cc)
OBJECTS  := $(patsubst $(SRCDIR)/%.cc, $(OBJDIR)/%.o, $(CC_FILES)) \

$(OBJECTS): | $(OBJDIR)

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS) | $(RUNDIR)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

$(OBJDIR):
	$(MAKEDIR_P) $@

$(RUNDIR):
	$(MAKEDIR_P) $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(HEADERS)
	@echo Building $@
	@$(CC) -c $(CFLAGS) $< -o $@

plotall:
	plot/./elliptic-orbits.py
	./plot/plot-exp-parabolic-test.py
	./plot/plot-exp-parabolic.py
	./plot/plot-exp-test.py
	./plot/plot-exp-testpart.py

runall:
	./$(TARGET) run 1
	./$(TARGET) run 2
	./$(TARGET) run 3
	./$(TARGET) run 4

test: default
	./$(TARGET) test

clean:
	-rm -rf $(OBJDIR) $(TARGET)
