####################
## Makefile Setup ##
####################

# Get the host-name if empty
ifeq ($(host-name),)
	host-name := $(shell hostname)
endif
# Get the kernel-name if empty
ifeq ($(kernel-name),)
	kernel-name := $(shell uname -s)
endif
# Get the arch-name if empty
ifeq ($(arch-name),)
	arch-name := $(shell uname -p)
endif

# Define the C++ compiler to use
CXX := $(shell which g++)

# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPSFILE = $(DEPSDIR)/$(notdir $*.d)
DEPSFLAGS = -MD -MF $(DEPSFILE)

# Define any directories containing header files
#   To include directories use -Ipath/to/files
INCLUDES += -I/home/vincent/software/mtl4 -I/home/vincent/software/thrust

# Define cxx compile flags
CXXFLAGS  := -std=c++11 -fopenmp -g -O0 -W -Wall -Wextra
CXXFLAGS  += -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP

# Define any directories containing libraries
#   To include directories use -Lpath/to/files
LDFLAGS +=

# Define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX
ifeq ($(kernel-name), Linux)
  LDLIBS += -lsfml-graphics -lsfml-window -lsfml-system -lX11 -lGL -lpthread
endif
ifeq ($(kernel-name), Darwin)
  LDLIBS += -lsfml-graphics -lsfml-window -lsfml-system -framework OpenGL
endif

######################
## Makefile Options ##
######################


####################
## Makefile Rules ##
####################

# Suffix replacement rules
#   $^: the name of the prereqs of the rule
#   $<: the name of the first prereq of the rule
#   $@: the name of the target of the rule

EXEC = $(basename $(wildcard *.cpp))

# 'make' - default rule
all: $(EXEC)

# Default rule for creating an exec of $(EXEC) from a .o file
$(EXEC): % : %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Default rule for creating a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -o $@ -c $<

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	-$(RM) *.o $(EXEC)
	$(RM) -r $(DEPSDIR)

# Define rules that do not actually generate the corresponding file
.PHONY: clean all

# Include the dependency files
-include $(wildcard $(DEPSDIR)/*.d)
