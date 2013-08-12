#
# 'make'        build all executable files
# 'make XXX'    build XXX executable
# 'make clean'  removes all .o and executable files
#

# Dependency directory
DEPSDIR := $(shell mkdir -p .deps; echo .deps)

# Define the C compiler to use
CXX := mpic++
LINK := $(CXX)

# Define any compile-time flags
CFLAGS := -fopenmp -funroll-loops -O3 -W #-Wall #-Wextra
# 'make DEBUG=1' - compile with debug flags
ifeq ($(DEBUG),1)
CFLAGS += -g -fno-inline
endif
# 'make PROFILE=1' - compile with profiling flags
ifeq ($(PROFILE),1)
CFLAGS += -g -pg
endif
# Dependency flags
DEPCFLAGS = -MD -MF $(DEPSDIR)/$*.d -MP

# Other in-code flags
CFLAGS +=

# define any directories containing header files other than /usr/include
#   include directories like -Ipath/to/files
INCLUDES = -I.

# define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX options
LDFLAGS +=

##################
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
##################

EXEC += Broadcast
EXEC += Generate
EXEC += Scatter
EXEC += Serial
EXEC += TeamScatter

# 'make' - default rule
all: $(EXEC)

# Rules for executables

Broadcast: Broadcast.o
	$(LINK) $(CFLAGS) -o $@ $^ $(LDFLAGS)

Generate: Generate.o
	$(LINK) $(CFLAGS) -o $@ $^ $(LDFLAGS)

Scatter: Scatter.o
	$(LINK) $(CFLAGS) -o $@ $^ $(LDFLAGS)

Serial: Serial.o
	$(LINK) $(CFLAGS) -o $@ $^ $(LDFLAGS)

TeamScatter: TeamScatter.o
	$(LINK) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# suffix replacement rule for building .o's from .cpp's
#   $<: the name of the prereq of the rule (a .cpp file)
#   $@: the name of the target of the rule (a .o file)
.cpp.o:
	$(CXX) $(CFLAGS) $(DEPCFLAGS) $(INCLUDES) -c -o $@ $<

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	-$(RM) *.o *~ */*~
	-$(RM) $(EXEC)
	$(RM) -r $(DEPSDIR)

DEPFILES := $(wildcard $(DEPSDIR)/*.d) $(wildcard $(DEPSDIR)/*/*.d)
ifneq ($(DEPFILES),)
-include $(DEPFILES)
endif

# define rules that do not actually generate the corresponding file
.PHONY: clean all
