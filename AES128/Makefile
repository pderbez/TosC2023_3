# Generic Makefile for compiling a simple executable.

#Polynomial used to define the finite field GF(2^n) -
# AES Field
GF_POLY := 0x11b
# F_2
#GF_POLY := 0x02
# Piccolo
#GF_POLY := 0x13

CC := g++
#CC := icpc
SRCDIR := src
BUILDDIR := build
USERDEFINES := -DGF_POLY=$(GF_POLY)

#CFLAGS := -g -mtune=native -march=native -Wall -Wno-unused-variable -std=c++11 -DNDEBUG
CFLAGS := -O3 -mtune=native -march=native -Wall -Wno-unused-variable -std=c++17 -DNDEBUG -fopenmp
#CFLAGS := -g -O2 -mtune=native -march=native -Wall -Wno-unused-variable -std=c++17 -DNDEBUG -fopenmp
#CFLAGS := -g -Wall -Wno-unused-variable -std=c++11

#LIBS := -ltbb -lm
LIBS :=

TARGET := aesCM

SOURCES := $(shell find $(SRCDIR) -type f -name *.cpp)
HEADERS := $(shell find $(SRCDIR) -type f -name *.hpp)



#OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.o))

#DEPS := $(OBJECTS:.o=.deps)
DEPS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.deps))


all: $(TARGET)

$(TARGET): $(OBJECTS)
	@echo " Linking GDSSearch..."; $(CC) $(USERDEFINES) $(CFLAGS) $^ $(LIBS) -o $(TARGET)



$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp  $(HEADERS)
	@mkdir -p $(BUILDDIR)
	@echo " CC $<"; $(CC) $(USERDEFINES) $(CFLAGS) -MD -MF $(@:.o=.deps) -c -o $@ $<

clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET) *~

-include $(DEPS)

$(BUILDDIR)/GFElement.o:   .FORCE

.FORCE:

.PHONY: clean .FORCE
