###########################################################################
# MLFit                                                                   #
# Nick Danielson / Amir Farbin                                            #
#   Build targets                                                         #
#     lib   - make libMLFit.a                                             #
#     shlib - make libMLFit.so                                            #
#     clean - delete all intermediate and final build objects             #
###########################################################################

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

# --- External configuration ----
CC           = g++
CCFLAGS      = -g -fPIC -O2 -Wno-deprecated -O -ansi -D_GNU_SOURCE 
MFLAGS       = -MM -Wno-deprecated
INCLUDES     = 
WORKDIR      = workdir
LIBDIR       = $(WORKDIR)
# -------------------------------------------------------------

#CCFLAGS  += $(ROOTCFLAGS)

# Internal configuration
PACKAGE=MLFit
LD_LIBRARY_PATH:=$(ROOTSYS)/lib
OBJDIR=$(WORKDIR)/objects/
INCDIR=$(WORKDIR)/$(PACKAGE)/
VPATH=$(INCDIR) $(OBJDIR)

INCLUDES += -I.. -I$(WORKDIR)/ -I$(ROOTSYS)/include
ROOTSYS  ?= ERROR_RootSysIsNotDefined
SKIPRDLIST = $(PACKAGE)_LinkDef.rdl RooNag.rdl
RDLLIST   = $(filter-out $(SKIPRDLIST),$(wildcard *.rdl))
SKIPLIST  = RooNag.cc
CINTFILE  = $(PACKAGE)Cint.cc
CINTOBJ   = $(PACKAGE)Cint.o
LIBFILE   = $(LIBDIR)/lib$(PACKAGE).a
SHLIBFILE = $(LIBDIR)/lib$(PACKAGE).so

default: shlib

# List of all includes files to copy from rdl
HHLIST=$(patsubst %.rdl,%.hh,$(RDLLIST))

# List of all source files to build
CCLIST=$(filter-out $(SKIPLIST),$(wildcard *.cc))

# List of all object files to build
OLIST=$(patsubst %.cc,%.o,$(CCLIST))

# List of all dependency file to make
DLIST=$(patsubst %.rdl,%.d,$(RDLLIST))

# Implicit rule making all dependency Makefiles included at the end of this makefile
$(INCDIR)/%.d: %.cc $(HHLIST)
	@echo "Making $@"
	@set -e; $(CC) $(MFLAGS) $(CCFLAGS) $(INCLUDES) $< \
	          | sed 's/\($(notdir $*)\)\.o[ :]*/\1.o $(notdir $@) : /g' > $@; \
	        [ -s $@ ] || rm -f $@

# Implicit rule copying all RDL to INCDIR/HH
%.hh: %.rdl 
	@mkdir -p $(INCDIR)
	@cp $< $(INCDIR)/$@

# Implicit rule to compile all classes
%.o : %.cc 
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $(OBJDIR)/$(notdir $@) $(INCLUDES)


# Rule to make ROOTCINT output file
$(OBJDIR)/$(CINTOBJ): $(RDLLIST) 
	@mkdir -p $(INCDIR)
	@mkdir -p $(OBJDIR)
	@echo "Running rootcint"
	@cd $(INCDIR) ; $(ROOTSYS)/bin/rootcint -f $(CINTFILE) -c $(INCLUDES) $(HHLIST)
	@echo "Compiling $(CINTFILE)"
	@$(CC) $(CCFLAGS) -c $(INCDIR)/$(CINTFILE) -o $(OBJDIR)/$(CINTOBJ) $(INCLUDES)

# Rule to combine objects into a library
$(LIBFILE): $(OLIST) $(INCDIR)/$(CINTFILE)) $(patsubst %.cc,%.o,$(OBJDIR)/$(CINTFILE))
	@echo "Making $(LIBFILE)"
	@rm -f $(LIBFILE)
	@ar q $(LIBFILE) $(addprefix $(OBJDIR)/,$(OLIST) $(patsubst %.cc,%.o,$(CINTFILE))) 
	@ranlib $(LIBFILE)

# Rule to combine objects into a shared library
$(SHLIBFILE): $(OLIST) $(patsubst %.cc,%.o,$(OBJDIR)/$(CINTFILE))
	@echo "Making $(SHLIBFILE)"
	@rm -f $(SHLIBFILE)
	@$(CC) $(addprefix $(OBJDIR)/,$(OLIST) $(patsubst %.cc,%.o,$(CINTFILE))) -bundle -undefined dynamic_lookup -o $(SHLIBFILE)

# Useful build targets
lib: $(LIBFILE) 
shlib: $(SHLIBFILE)
clean:
	rm -f $(OBJDIR)/*
	rm -f $(INCDIR)/*
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE)

.PHONY : shlib lib default clean

-include $(addprefix $(INCDIR)/,$(DLIST))
