MAKEDIR_P = mkdir -p

TARGET = data_generator
LIBS = -L /usr/lib64/atlas -lcblaskk
CC = g++
CFLAGS = -g -pedantic -Wall -O3
SRCDIR = src
OBJDIR = obj
TESTDIR = test
RUNDIR = 

.PHONY: default all clean

default: $(OBJDIR) $(TARGET)
all: clean default
runall: default

# Manually write down the headers and the main program files
SCR_HH  := linalg.hh pulsar.hh signal-model.hh signal-vectors.hh
SCR_CC  := $(patsubst %.hh, %.cc, $(SCR_HH))
TEST_HH := test-linalg.hh test-pulsar.hh test-signal-model.hh test-signal-vectors.hh
TEST_CC := $(patsubst %.hh, %.cc, $(TEST_HH)) test_runner.cc

OBJECTS := $(patsubst %.cc, $(OBJDIR)/%.o, $(SCR_CC)) \
           $(patsubst %.cc, $(OBJDIR)/%.o, $(TEST_CC))

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

test: default
	./$(TARGET) test

clean:
	-rm -rf $(OBJDIR) $(TARGET)
