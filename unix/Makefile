#----------------------------------------------------------------------
#   wan3, weq3, (c) Ilya Valuev 1997-2007
#----------------------------------------------------------------------
# This file jumps into the architecture-specific directory and 
# builds a program.
#
# This is taken from Paul Smith, GNU
#
SHELL = /bin/tcsh
# The variable DEFAULT_ARCH contains the default value.
# To choose a different value, users can edit this, 
# or invoke make using
#     make _ARCH = <user's architecture>
#
# To build parallel version use:
#     make _MPI=true _ARCH = <user's architecture>

DEFAULT_ARCH = gnu

ifeq (,$(filter Obj_%,$(notdir $(CURDIR))))
  include Makefile.target
else

#----- End Boilerplate

# Make the first (i.e. default) target all, 
# then use special target to declare it phony to 
# avoid searching for non-existent file
all::
.PHONY: all


# This section is for the architecture-specific target directory

# This sub-section is common to all architectures

all:: equil
        
TESTNAME  = wequil2

ROOTDIR   = ../../..
MAINDIR = $(ROOTDIR)/wequil2
PROPDIR = $(MAINDIR)/props
UTILDIR = $(ROOTDIR)/ivutils
INCDIR  = $(UTILDIR)/include
TSTDIR  = $(ROOTDIR)/applications/$(TESTNAME)

# VPATH tells make where to look for source files.

VPATH   = $(MAINDIR):$(UTILDIR)/src:$(UTILDIR)/src/graph:$(UTILDIR)/src/four:$(UTILDIR)/src/rparam:$(PROPDIR)


SRC_EQ = \
        cequilr.cpp  \
        analyse.cpp  \
        interact.cpp  \
        pcycle.cpp \
        plasma.cpp \
        plasmar.cpp

SRC_AN = \
        analyser.cpp \
        ancurcur.cpp \
        andstruc.cpp \
        antrjlog.cpp \
        anpair.cpp \
        antherm.cpp \
        anmfield.cpp \
        anvcorr.cpp \
        anvdistr.cpp \
        interact.cpp \
        nanalyse.cpp \
        pcycle.cpp \
        plasma.cpp \
        plasmar.cpp

# add more headers for packing  
INC_EQ = \
        $(PHOTDIR)/analyse.h \

SRC_UTILS = \
        vector_3.cpp  \
        rparam.c  \
        common.cpp  \
        four.c  \
        erf_namd.c  \
        statist.cpp \
        record.cpp \

INC_UTILS = \
        $(INCDIR)/common.h \
        $(INCDIR)/vector_3.h \
        $(INCDIR)/rparam.h \
        $(INCDIR)/plane_3.h \
        $(INCDIR)/refobj.h \
        $(INCDIR)/logger.h \
        $(INCDIR)/pencil.h \
        $(INCDIR)/contour.h \
        $(INCDIR)/region.h \
        $(INCDIR)/seqpack.h \
        $(INCDIR)/grid.h \

             
SRC_AEQ = $(SRC_EQ) $(SRC_UTILS)
OBJ_AEQ1 = $(SRC_AEQ:.cpp=.o)
OBJ_AEQ = $(OBJ_AEQ1:.c=.o)
      

SRC_AAN = $(SRC_AN) $(SRC_UTILS)
OBJ_AAN1 = $(SRC_AAN:.cpp=.o)
OBJ_AAN = $(OBJ_AAN1:.c=.o)

OBJ = $(OBJ_AAN) $(OBJ_AEQ)

_MPI = false
ifeq (true,$(_MPI))
MESSAGE = 'Using MPI'
MPI_PREPROC_FLAG = -DUSE_MPI
endif

  

#
# Note: here are the preprocessor flags
# These need to be applied manually, using e.g. 
#
# touch a.cpp; make
#
#alg_parser.o: alg_parser.cpp
#       $(CPP) $(CPPFLAGS) -O0 -c $<


# This subsection contains branches to specific architecture definitions
#
# Initialize with values that are guaranteed to fail
# This is used to catch invalid _ARCH values
  C = "Invalid_ARCH_value"
  CPP = "Invalid_ARCH_value"
  F77 = "Invalid_ARCH_value"
# Add null defaults for rarely used variables
  RUNDIR = .
  SLEEP = 
  POSTLINK = 

include $(SRCDIR)/Makefile.arch


# This sub-section is common to all architectures
.SUFFIXES : .o .cpp .c .F
# inference rule to get .o targets from .cpp prerequisites
.cpp.o :
	$(CPP) $(CPPFLAGS) -c $<
# inference rule to get .o targets from .c prerequisites
.c.o :
	$(C) $(CFLAGS) -c $<

# Target rules for executables specify prerequisite objects
# and linking command.

equil: $(OBJ_AEQ)
	$(CPP) $(OFLAG) $(OBJ_AEQ) \
        $(USRLIB) -o equil.exe
        $(POSTLINK);


analyse: $(OBJ_AAN)
	$(CPP) $(OFLAG) $(OBJ_AAN) \
        $(USRLIB) -o analyse.exe; \
        $(POSTLINK);

#lib: $(OBJ)
#       ar -r libgrasp.a $(OBJ) $(OBJ_REAX)

#  include $(SRCDIR)/Makefile.testing


zip      :
	tar -cPvzf builder.tgz $(PHOTDIR)/*.cpp $(PHOTDIR)/*.h $(INCDIR)/*.h $(UTILDIR)/src/*.cpp \
         $(UTILDIR)/doc/*.txt $(UTILDIR)/src/graph/*.cpp $(UTILDIR)/src/rparam/*.c \
         $(TSTDIR)/doc/*.* $(TSTDIR)/vs9/*.* $(TSTDIR)/*.gnu $(TSTDIR)/*.h $(TSTDIR)/*.cpp \
         $(TSTDIR)/*.ini $(TSTDIR)/unix/Makefile*


tar_full:
	mkdir fdtd_tmp; \
        cd fdtd_tmp; \
        cvs -d :ext:valuev@guru.hepti.kiae.ru:/home/dev/Photon \
                checkout newfdtd; \
        tar -zcvf newfdtd.tgz  \
        mv newfdtd.tgz ../; \
        cd .. ; \
        \rm -r fdtd_tmp;


# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
   -include $(DEPENDS)
   
#----- Begin Boilerplate
endif
