
EXTRA_FLAGS = -O2

PROFILEFLAG = 


if BIT64
MYFFLAGS = -c -static -fno-automatic -m64 -fsecond-underscore $(EXTRA_FLAGS) 
MYLDFLAGS = -m64 $(PROFILEFLAG)
MYCFLAGS = -m64 $(EXTRA_FLAGS) -fcommon
else
# doubt we'll ever need 32 bit, but just in case will likely need to hack this
MYFFLAGS = -c -static -fno-automatic -fsecond-underscore $(EXTRA_FLAGS)
MYCFLAGS = -static -m64 $(EXTRA_FLAGS) 
MYLDFLAGS = -m64 -lgfortran -lpthread $(PROFILEFLAG)
endif

if USE_ROOT
ROOTLIBS = -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -ldl
IROOT = -I/$(ROOT_DIR)/include 
else
ROOTLIBS =
IROOT = 
endif

AM_FFLAGS = $(MYFFLAGS)
AM_CFLAGS = $(MYCFLAGS)
AM_LDFLAGS = $(MYLDFLAGS)

# tack on flag to provide git version to C/fortran codes
AM_CFLAGS += -DGIT_SNANA_VERSION=\"$(GIT_VERSION)\"

sntable_dump_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm $(GSLCBLASLIB) -lpthread -lz -lstdc++ $(ROOTLIBS)

sntable_combine_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm $(GSLCBLASLIB) -lpthread -lz -lstdc++ $(ROOTLIBS)

# --------------------------------------------
bin_PROGRAMS = sntable_dump sntable_combine


SNTOOLS_OUTPUT = \
        sntools_output.cpp

sntable_dump_SOURCES = sntable_dump.c sntools.c sntools_output.cpp


sntable_combine_SOURCES = sntable_combine.c sntools.c sntools_output.cpp
