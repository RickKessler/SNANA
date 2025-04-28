# $Source: /cvs/sdss/snana/src/Makefile,v $
# $Revision: 1.500 $    Last Modified: $Date: 2019/05/02 20:24:47 $
#
#
#  R.Kessler 
#
#  HISTORY 
# -----------------------------------
#
# Feb 05 2021: rename sntools_fitsio.c[h] -> sntools_dataformat_fits.c[h]
#
# Mar 16 2021: 
#   + remove obsolete install feature
#   + finall add debug flag;  make debug
#
# Jun 24 2021: 
#   + include sntools_genGauss_asym dependence for SALT2mu
#
# Jul 17 2021: make sntools_gridmap.o
# Dec 10 2021: incorporate sntools_data.o
# Dec 22 2021: incorporate sntools_sim_readme.c[h]
# Apr 03 2022: link SIMSED_xxx programs to sntools_output because 
#              sntools.c util now calls AUTOSTORE util.
#
# Jun 10 2022: for deciding C-python interace, check $VERSION_LIBPYTHON
#    instead of $SNANA_PYTHON_DIR. 
#
# Sep 15 2022: define EXTRA_FLAGS_CPLUS = -std=c++17 (for ROOT v6 at NERSC)
# Nov 22 2022: add USE_BAYESN flag to use libyaml
# Nov 23 2022: add genmag_extinction
# Nov 29 2022: remove kcordump
# Dec 02 2022: append CFLAGS with snana git version
# May 25 2023: incluce yaml.o ingested from github 
#          (original need is genmag_BAYESN.c)
# May 26 2023: link sim to sntools_sim_atmosphere
# Sep 14 2023: add genmag_extrap
# Feb 22 2024: add sntools_wgtmap
# Feb 24 2024: switch -O1 to -O2 optimizer (20% faster photo-z fits)
# Mar 22 2024: fix kcor to link with genmag_SEDtools.o
# Apr 12 2024: few tweaks to handle sntools_wgtmap.c[h]
# Apr 13 2025: add cnpy for wfit to read cov matrix in optional npz format
# Apr 24 2025: fix make for wfit to work on both midway and perlmutter
#             
# ------------------------------------------------------------------------

SHELL = /bin/sh

BITNESS = $(shell uname -p)
OSNAME  = $(shell uname -s)

# Ubuntu returns unknown for uname -p, so fall back to uname -m
ifeq ($(BITNESS), unknown)
   BITNESS = $(shell uname -m)
endif


# ================================

FFC  = gfortran
#CC   = cc
CC   = gcc

 EXTRA_FLAGS =
#  EXTRA_FLAGS = -O1    # -ftree-vectorize -ftree-vectorizer-verbose=2
  PROFILEFLAG = 

EXTRA_FLAGS = -O2  # xxx RESTORE
#EXTRA_FLAGS = -g -fbounds-check  # valgrind --leak-check=full 

# check if gcc version is 10 or higher (9.2020)
GCCVERSION10 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 10)

ifeq ($(GCCVERSION10),0)
   EXTRA_FLAGS_FORTRAN  = $(EXTRA_FLAGS)
   EXTRA_FLAGS_C        = $(EXTRA_FLAGS)
   EXTRA_FLAGS_CPLUS    = 
else 
# for gcc v10, allow argument-mismatch errors (9.2020)
   EXTRA_FLAGS_FORTRAN = $(EXTRA_FLAGS) -fallow-argument-mismatch -fcommon
   EXTRA_FLAGS_C       = $(EXTRA_FLAGS) -fcommon
   EXTRA_FLAGS_CPLUS   = -std=c++17
endif

# after running code with -pg, "gprof $SNANA_DIR/bin/[code].exe"
#  PROFILEFLAG = -pg
# ------------------------
# compile default storage options for both ROOT and HBOOK ....
USE_ROOT   = 1
USE_HBOOK  = 1
USE_BAYESN = 1

# ... unless the product is not installed.
ifndef CERN_DIR
   USE_HBOOK = 0
endif

ifndef ROOT_DIR
   USE_ROOT = 0
endif


# ------------------------------------------
# define libraries

# default lib-subdir name is just 'lib'.
# To see where include files are pulled from, add -H to compiler flags

  LIB_SDIR  = lib

  LGSL    =  -L$(GSL_DIR)/$(LIB_SDIR) -lgsl  \
	     -L$(GSL_DIR)/$(LIB_SDIR) -lgslcblas

  IGSL    =  -I$(GSL_DIR)/include

  LCFITSIO =  -L$(CFITSIO_DIR)/$(LIB_SDIR) -lcfitsio
  ICFITSIO =  -I/$(CFITSIO_DIR)/include

# ---------------

ifeq ($(USE_HBOOK),0)
  LCERN  = 
else
  LCERN  = -L$(CERN_DIR)/$(LIB_SDIR) -lpacklib -lmathlib -lkernlib -lncurses
endif

# ----------------

# USE_ROOT flag affects the compilation of sntools_ouptut.c,
# and the -lstdc++  & root libraries for linking.

ifeq ($(USE_ROOT),0)
   CCoutput = cc
   CPPLIB   = 
   SNTOOLS_ROOT = 
else
   CCoutput    = c++ $(EXTRA_FLAGS_CPLUS)
   CPPLIB      = -lstdc++
   IROOT       = -I/$(ROOT_DIR)/include -I/opt/local/include/freetype2
   ROOT_LIBDIR = $(LIB_SDIR)
   ROOT_LIBDIR = lib

# include lCint if it exists (root v5), and exclude if not (root v6)
# 9/22/2020; need -ltbb for root 6
   LCINT = -lCint
   LIBCINT_ROOT5 = $(ROOT_DIR)/$(ROOT_LIBDIR)/libCint.so
   ifeq (,$(wildcard $(LIBCINT_ROOT5)))
	LCINT =  -ltbb
   endif	 

   LROOT = -L$(ROOT_DIR)/$(ROOT_LIBDIR) -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -ldl $(LCINT)
   SNTOOLS_ROOT = $(SRC)/sntools_output_root.c 
endif

#echo " LIBINT_ROOT5 = $(LIBCINT_ROOT5) " ;

# Apr 2019: check for python codes
# Jun 2022: use VERSION_LIBPYTHON to decide on C-python interface
ifndef VERSION_LIBPYTHON
  IPY = 
  LPY = 
else
  IPY = -I/$(CONDA_PREFIX)/include/python$(VERSION_LIBPYTHON)
  LPY = -L$(CONDA_PREFIX)/lib -ldl  -lutil -Wl,-rpath,$(CONDA_PREFIX)/lib -lpython$(VERSION_LIBPYTHON)
endif

# --------------------------------------------
# set compile flags.
ifeq ($(BITNESS),i686)
   SNFFLAGS  = -c -static -fno-automatic -fsecond-underscore $(EXTRA_FLAGS_FORTRAN)
   SNCFLAGS  = -c -static -m64 $(EXTRA_FLAGS)
   SNLDFLAGS = -m64 -lg2c -lpthread $(PROFILEFLAG)
else
   SNFFLAGS  = -c -static -fno-automatic -m64 -fsecond-underscore $(EXTRA_FLAGS_FORTRAN)
   SNCFLAGS  = -c -static -m64 $(EXTRA_FLAGS_C)
   SNLDFLAGS = -m64 -lpthread $(PROFILEFLAG)
endif


# From S. Rodney (May 2014):
# OSX returns i386 for uname -p, even when 64-bit, so fall back to uname -m
ifeq ($(OSNAME),Darwin)
   OSXBITNESS = $(shell uname -m)
ifeq ($(OSXBITNESS),x86_64)
   SNFFLAGS  = -c -static -fno-automatic -fsecond-underscore -m64
   SNCFLAGS  = -c -mdynamic-no-pic -m64
   SNLDFLAGS =
else
   SNFFLAGS  = -c -static -fno-automatic -fsecond-underscore -m64
   SNCFLAGS  = -c -mdynamic-no-pic -m64
   SNLDFLAGS =
endif
endif

#GIT_VERSION := "$(shell git describe --dirty --always --tags)"
GIT_VERSION := "$(shell git describe --always --tags)"
SNCFLAGS += -DGIT_SNANA_VERSION=\"$(GIT_VERSION)\"

# ===========================

UTIL = ../util
BIN  = ../bin
OBJ  = ../obj
SRC  = ../src
LIB  = ../lib
TMPANA = $(SRC)/tempana
TMPFIT = $(SRC)/tempfit

# --------------------------------------------

# define fortran pre-processor tools (from obsolete cernlib)
FPP_EXECUTABLES = $(BIN)/fcasplit

# Define snana executables 
SNANA_EXECUTABLES	=		\
	$(BIN)/snlc_sim.exe		\
	$(BIN)/snana.exe		\
	$(BIN)/snlc_fit.exe		\
	$(BIN)/psnid.exe		\
	$(BIN)/kcor.exe			\
	$(BIN)/nearnbr_maxFoM.exe	\
	$(BIN)/nearnbr_apply.exe	\
	$(BIN)/wfit.exe			\
	$(BIN)/SALT2mu.exe		\
	$(BIN)/combine_fitres.exe 	\
	$(BIN)/sntable_dump.exe 	\
	$(BIN)/sntable_combine.exe 	\
	$(BIN)/SIMSED_fudge.exe		\
	$(BIN)/SIMSED_extractSpec.exe	\
	$(BIN)/SIMSED_check.exe		\
	$(BIN)/SIMSED_rebin.exe		\
	$(BIN)/simlib_coadd.exe		\
	$(BIN)/unfold_snpar.exe

# ----------------------
# check for executables that depend explicitly on output format

ifeq ($(USE_HBOOK),1)
	HBOOK_EXECUTABLES = $(BIN)/merge_hbook.exe
else
	HBOOK_EXECUTABLES =
endif

ifeq ($(USE_ROOT),1)
	ROOT_EXECUTABLES = $(BIN)/merge_root.exe
else
	ROOT_EXECUTABLES = 
endif


EXECUTABLES = \
	$(SNANA_EXECUTABLES) 	\
	$(ROOT_EXECUTABLES) 	\
	$(HBOOK_EXECUTABLES)

# -------------------------------------
# Define list of object files for each program.

# start with lists that are used in multiple programs
OBJ_GENMAG_FIT	=		\
	$(OBJ)/genmag_stretch2.o \
	$(OBJ)/genmag_mlcs2k2.o \
	$(OBJ)/genmag_SALT2.o	\
	$(OBJ)/genmag_BAYESN.o	\
	$(OBJ)/genmag_snoopy.o	\
	$(OBJ)/genmag_S11DM15.o	\
	$(OBJ)/genmag_SEDtools.o \
	$(OBJ)/genmag_extinction.o \
	$(OBJ)/genmag_extrap.o


OBJ_GENMAG_SIM =	\
	$(OBJ_GENMAG_FIT) \
	$(OBJ)/genmag_NON1ASED.o \
	$(OBJ)/genmag_NON1AGRID.o \
	$(OBJ)/genmag_LCLIB.o	\
	$(OBJ)/genmag_PySEDMODEL.o \
	$(OBJ)/genmag_SIMSED.o

OBJ_SNTOOLS_SNANA	 = 	\
	$(OBJ)/sntools.o 	\
	$(OBJ)/sntools_gridmap.o \
	$(OBJ)/sntools_nonlinearity.o  	\
	$(OBJ)/sntools_fluxErrModels.o	  \
	$(OBJ)/sntools_fluxErrModels_legacy.o	\
	$(OBJ)/sntools_simeff_obsolete.o	\
	$(OBJ)/sntools_dataformat_fits.o	\
	$(OBJ)/sntools_dataformat_text.o	\
	$(OBJ)/sntools_data.o 			\
	$(OBJ)/sntools_modelgrid_read.o \
	$(OBJ)/sntools_calib.o \
	$(OBJ)/sntools_nearnbr.o \
	$(OBJ)/sntools_output.o  \
	$(OBJ)/yaml.o	\
	$(OBJ)/minuit.o		\
	$(OBJ)/MWgaldust.o	\
	$(OBJ)/sntools_cosmology.o \
	$(OBJ)/multiseason.o

OBJ_SNTOOLS_DATA =  \
	$(OBJ)/sntools_dataformat_fits.o        \
        $(OBJ)/sntools_dataformat_text.o        \
        $(OBJ)/sntools_data.o                   \
	$(OBJ)/sntools_output.o			\
	$(OBJ)/sntools_cosmology.o
 
OBJ_OUTPUT = $(OBJ)/sntools_output.o

# now the program-specific lists

OBJ_SNANA	 = 	\
	$(OBJ)/snana.o	\
	$(OBJ_SNTOOLS_SNANA)

OBJ_SNFIT	 = 		\
	$(OBJ)/snlc_fit.o	\
	$(OBJ_SNTOOLS_SNANA)	\
	$(OBJ)/sntools_genSmear.o \
	$(OBJ)/sntools_devel.o    \
	$(OBJ)/sntools_zPDF_spline.o \
	$(OBJ_GENMAG_FIT)

OBJ_PSNID	=  		\
	$(OBJ)/psnid.o		\
	$(OBJ_SNTOOLS_SNANA)	\
	$(OBJ)/psnid_tools.o 	\
	$(OBJ)/psnid_BEST.o 	\
	$(OBJ)/psnid_BEST2.o

OBJ_SIM	 = 	\
	$(OBJ)/snlc_sim.o	\
	$(OBJ_SNTOOLS_SNANA)	\
	$(OBJ)/sntools_host.o   \
        $(OBJ)/sntools_wgtmap.o   \
	$(OBJ)/sntools_weaklens.o  \
	$(OBJ)/sntools_stronglens.o  \
	$(OBJ)/sntools_wronghost.o   \
	$(OBJ)/sntools_genSmear.o  \
	$(OBJ)/sntools_devel.o  \
	$(OBJ)/sntools_sim_readme.o \
	$(OBJ)/sntools_sim_atmosphere.o \
	$(OBJ)/sntools_trigger.o  \
	$(OBJ)/sntools_modelgrid_gen.o  \
	$(OBJ)/sntools_spectrograph.o \
	$(OBJ)/sntools_genPDF.o \
	$(OBJ)/sntools_genGauss_asym.o \
        $(OBJ)/sntools_genExpHalfGauss.o \
	$(OBJ)/inoue_igm.o \
	$(OBJ_GENMAG_SIM) 

# -------
# list of dependencies for sntools_output.o

SNTOOLS_OUTPUT = \
	$(SRC)/sntools_output.c \
	$(SRC)/sntools_output.h \
	$(SRC)/sntools_output_hbook.c \
	$(SRC)/sntools_output_marz.c \
	$(SRC)/sntools_output_text.c \
	$(SNTOOLS_ROOT)

# ------------------------------------------------------

LIBSNANA  =   $(LIB)/libsnana.a
LIBSNFIT  =   $(LIB)/libsnfit.a

CLEAN_FILES =  $(FPP_EXECUTABLES) $(LIBSNANA)  $(LIBSNFIT)  $(EXECUTABLES)

# --------------------------

all: set_Cpreproc_flags makeDirs $(FPP_EXECUTABLES) $(LIBSNANA)  $(LIBSNFIT)  $(EXECUTABLES)

debug:	EXTRA_FLAGS = -g -fbounds-check  
debug:  all

clean:
	rm -f $(CLEAN_FILES) $(OBJ)/*.o
	rm -rf $(TMPANA)
	rm -rf $(TMPFIT)

# -----------------------------------
# set define statements in sntools_output.c to pick HBOOK and/or ROOT
# based on existance of $CERN_DIR and $ROOT_DIR.
set_Cpreproc_flags :
	$(UTIL)/set_Cpreproc_flags.py

# May 3 2019: create blank dirs that are not in Git
makeDirs:
	mkdir -p $(OBJ)
	mkdir -p $(BIN)
	mkdir -p $(LIB)	

# ------------------------
# fcasplit (April 2014)

$(OBJ)/fcasplit.o : $(SRC)/fcasplit.F
	(cd $(OBJ);  $(FFC)  $(SNCFLAGS) $(SRC)/fcasplit.F ) 

$(BIN)/fcasplit : $(OBJ)/fcasplit.o 
	$(FFC) -o  $@ $(SNLDFLAGS) \
	$(OBJ)/fcasplit.o 
	(cd $(OBJ);  rm fcasplit.o ) 

# -------------------------
$(OBJ)/minuit.o :  minuit.F
	(cd $(OBJ); \
	$(FFC) $(SNFFLAGS) $(SRC)/minuit.F)

# -------------------------
$(OBJ)/yaml.o : $(SRC)/yaml.c $(SRC)/yaml.h
	(cd $(OBJ); $(CC) $(SNCFLAGS)  $(SRC)/yaml.c )

$(OBJ)/cnpy.o : $(SRC)/cnpy.cpp $(SRC)/cnpy.h
	(cd $(OBJ); $(CC) $(SNCFLAGS) -I/$(CONDA_PREFIX)/include  $(SRC)/cnpy.cpp )

# -----------------------------------------------------
#  cosmology fitter (from Alex Connelly, modified by G.Miknaitis & RK)

$(OBJ)/wfit.o : $(SRC)/wfit.c $(SRC)/sntools.c
	(cd $(OBJ);  $(CC)  $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/wfit.c ) 

$(OBJ)/sntools_npz.o : $(SRC)/sntools_npz.c
	(cd $(OBJ); $(CCoutput) $(SNCFLAGS) -I/$(CONDA_PREFIX)/include $(SRC)/sntools_npz.c )

$(BIN)/wfit.exe : $(OBJ)/wfit.o $(OBJ)/cnpy.o $(OBJ)/sntools.o  $(OBJ)/sntools_npz.o $(OBJ)/sntools_output.o $(OBJ)/sntools_cosmology.o
	$(FFC) -o  $@ $(SNLDFLAGS) \
	$(OBJ)/wfit.o  \
	$(OBJ)/sntools.o  \
	$(OBJ)/sntools_output.o  \
	$(OBJ)/sntools_npz.o 	\
	$(OBJ)/cnpy.o 		\
	$(OBJ)/sntools_cosmology.o \
	$(LCERN) $(LROOT) $(CPPLIB) \
	$(LCFITSIO)  \
	-lm $(LGSL) \
	-lz 
	(cd $(OBJ);  rm wfit.o )


# -----------------------------------------------------
#  MCMC cosmology fitter (Parkinson, Marriner)
#  Mar 2012: Finally just pure C code (no more f90)

$(OBJ)/sncosmo_mcmc.o : $(SRC)/sncosmo_mcmc.c $(SRC)/sntools.c
	(cd $(OBJ);  $(CC)  $(SNCFLAGS) $(SRC)/sncosmo_mcmc.c ) 

$(BIN)/sncosmo_mcmc.exe : $(OBJ)/sncosmo_mcmc.o $(OBJ)/sntools.o
	$(FFC) -o  $@ $(SNLDFLAGS) \
	$(OBJ)/sncosmo_mcmc.o  \
	$(OBJ)/sntools.o \
	$(OBJ)/sntools_output.o \
	$(LCERN) $(LROOT) -lm $(LGSL) 
	(cd $(OBJ);  rm sncosmo_mcmc.o ) 

# -------------------------------------------
# SALT2mu : convert SALT2 parameters (c,x1,x0) into distances


$(OBJ)/SALT2mu.o : $(SRC)/SALT2mu.c $(SRC)/sntools.c $(SRC)/sntools_output.c $(SRC)/sntools_genPDF.c  $(SRC)/sntools_genGauss_asym.c $(SRC)/minuit.F
	(cd $(OBJ);  $(CC)  $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(SRC)/SALT2mu.c ) 

$(BIN)/SALT2mu.exe : $(OBJ)/SALT2mu.o  $(OBJ)/sntools.o $(OBJ)/sntools_output.o $(OBJ)/minuit.o $(OBJ)/sntools_gridmap.o $(OBJ)/sntools_genGauss_asym.o $(OBJ)/sntools_genExpHalfGauss.o
	$(FFC) -o  $@ $(SNLDFLAGS) \
	$(OBJ)/SALT2mu.o  \
	$(OBJ)/sntools.o \
	$(OBJ)/sntools_output.o \
	$(OBJ)/sntools_gridmap.o \
	$(OBJ)/sntools_genGauss_asym.o \
	$(OBJ)/sntools_genExpHalfGauss.o \
	$(OBJ)/minuit.o 	\
	$(LCERN) $(LROOT) -lm $(LGSL) $(LCFITSIO) $(CPPLIB) -lpthread
	(cd $(OBJ);  rm SALT2mu.o ) 

# -----------------------------------------------------
# May 2006: 
#   make libsnana.a fortran library for access by C-code
#   Note that "main.o" is excluded from library
#

$(LIB)/libsnana.a:   snana.car
	for item in  snana ; do \
	    rm -f $$item.lnk; \
	    ln -s $(SRC)/$$item.car $$item.lnk; \
	done
	@ rm -rf $(TMPANA)
	mkdir $(TMPANA)
	$(UTIL)/ypatchy.pl - $(TMPANA)/snana.f snana.cra .go
	rm snana.lnk 
	(cd $(TMPANA); \
	    ../$(BIN)/fcasplit -f gfortran snana.f "$(SNFFLAGS)" " " ; \
	    echo "Compiling"; \
	    . ./snana.shfca; \
	)
	echo "Creating Library"; \
	rm -f $(TMPANA)/main.o
	ar cr $@ $(TMPANA)/*.o
	rm -f $(TMPANA)/*.o
	rm -f $(TMPANA)/*.f

# -------------------------------------------------
# May 15 2008:   make libsnfit.a 
#

$(LIB)/libsnfit.a:   snlc_fit.car snana.car
	for item in snana snlc_fit ; do \
	    rm -f $$item.lnk; \
	    ln -s $(SRC)/$$item.car $$item.lnk; \
	done
	@ rm -rf $(TMPFIT)
	mkdir $(TMPFIT)
	$(UTIL)/ypatchy.pl - $(TMPFIT)/snlc_fit.f snlc_fit.cra .go
	rm snana.lnk snlc_fit.lnk
	(cd $(TMPFIT); \
	    ../$(BIN)/fcasplit -f gfortran snlc_fit.f "$(SNFFLAGS)" " " ;  \
	    echo "Compiling"; \
	    .  ./snlc_fit.shfca; \
	)
	echo "Creating Library"; \
	rm -f $(TMPFIT)/main.o
	ar cr $@ $(TMPFIT)/*.o
	rm -f $(TMPFIT)/*.o
	rm -f $(TMPFIT)/*.f

# -------------------------------------------------
# SN analysis shell (reads DATA Files, makes cuts ...)

$(OBJ)/snana.o :  snana.f  
	(cd $(OBJ); $(FFC) $(SNFFLAGS) $(SRC)/snana.f)

snana.f : snana.cra  snana.cra snana.car  
	for item in snana ; do \
	    rm -f $$item.lnk; \
	    ln -s $(SRC)/$$item.car $$item.lnk; \
	done
	$(UTIL)/ypatchy.pl - $@ snana.cra .go
	rm snana.lnk

$(BIN)/snana.exe :  $(OBJ_SNANA)
	$(FFC) -o $@ $(SNLDFLAGS)  \
	$(OBJ_SNANA) \
	$(LCFITSIO)  $(LGSL) \
	$(LCERN) $(LROOT) -lm  $(CPPLIB)
	(cd $(OBJ); rm snana.o )	

# -------------------------------------------------
# SN light curve fitter 

$(OBJ)/snlc_fit.o :  snlc_fit.f 
	(cd $(OBJ); $(FFC) $(SNFFLAGS) $(SRC)/snlc_fit.f)

snlc_fit.f : snlc_fit.cra  snlc_fit.cra  snlc_fit.car snana.car 
	for item in  snana snlc_fit; do \
	    rm -f $$item.lnk; \
	    ln -s $(SRC)/$$item.car $$item.lnk; \
	done
	$(UTIL)/ypatchy.pl - $@ snlc_fit.cra .go
	rm snana.lnk  snlc_fit.lnk 

$(BIN)/snlc_fit.exe : $(OBJ_SNFIT)
	$(FFC) -o $@ $(SNLDFLAGS)  \
	$(OBJ_SNFIT) 	\
	$(LCFITSIO)  $(LGSL) 	\
	$(LCERN)  $(LROOT) -lm $(CPPLIB)
	(cd $(OBJ); rm snlc_fit.o )

# -------------------------------------------------
# PSNID photometric classifier (Jun 2012)

$(OBJ)/psnid_tools.o : $(SRC)/psnid_tools.c $(SRC)/sntools.h \
	$(SRC)/sntools_modelgrid_read.c
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/psnid_tools.c )

$(OBJ)/psnid_BEST.o : $(SRC)/psnid_BEST.c $(SRC)/sntools.h \
	$(SRC)/sntools_modelgrid_read.c
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/psnid_BEST.c )


$(OBJ)/psnid_BEST2.o : $(SRC)/psnid_BEST2.c $(SRC)/sntools.h \
	$(SRC)/sntools_modelgrid_read.c
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/psnid_BEST2.c )


$(OBJ)/psnid.o :  psnid.f
	(cd $(OBJ); $(FFC) $(SNFFLAGS) $(SRC)/psnid.f)

psnid.f : psnid.cra  psnid.cra  psnid.car  snana.car  snlc_fit.cra snlc_fit.car
	for item in  snana  snlc_fit psnid ; do \
	    rm -f $$item.lnk; \
	    ln -s $(SRC)/$$item.car $$item.lnk; \
	done
	$(UTIL)/ypatchy.pl - $@ psnid.cra .go
	rm snana.lnk snlc_fit.lnk psnid.lnk

$(BIN)/psnid.exe : $(OBJ_PSNID) 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ_PSNID) 		\
	$(LCFITSIO) $(LGSL) 	\
	$(LCERN)  $(LROOT) -lm  $(CPPLIB)
	(cd $(OBJ); rm psnid.o)

# -------------------------------------------------
# K corrections and LC templates 

$(OBJ)/kcor.o : $(SRC)/kcor.c $(SRC)/kcor.h  $(SRC)/sntools_spectrograph.c
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO)  $(IGSL) $(SRC)/kcor.c )


$(BIN)/kcor.exe : $(OBJ)/kcor.o $(OBJ)/sntools.o $(OBJ)/MWgaldust.o $(OBJ)/sntools_spectrograph.o $(OBJ)/genmag_SEDtools.o $(OBJ_SNTOOLS_DATA)
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/kcor.o \
	$(OBJ)/sntools.o  	\
	$(OBJ)/MWgaldust.o 	\
	$(OBJ_SNTOOLS_DATA) 	\
	$(OBJ)/sntools_spectrograph.o  	\
	$(OBJ)/genmag_SEDtools.o  	\
	$(LCFITSIO)  $(LGSL) $(LCERN) $(LROOT) -lm $(CPPLIB)
	(cd $(OBJ); rm kcor.o )

# -------------------------------------------------
# kcordump program (removed Nov 29 2022)

# -------------------------------------------------
# combine_fitres (C-version with no array limits)

$(OBJ)/combine_fitres.o :  $(SRC)/combine_fitres.c
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(SRC)/combine_fitres.c )

$(BIN)/combine_fitres.exe : \
	$(OBJ)/combine_fitres.o $(OBJ)/sntools.o  $(OBJ)/sntools_output.o 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/combine_fitres.o \
	$(OBJ)/sntools.o  	\
	$(OBJ_OUTPUT)		\
	$(LGSL) $(LCERN) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	(cd $(OBJ); rm combine_fitres.o )

# --------------------
# sntable_dump (Apr 2013)

$(OBJ)/sntable_dump.o : $(SRC)/sntable_dump.c $(SNTOOLS_OUTPUT)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) \
	$(SRC)/sntable_dump.c )

$(BIN)/sntable_dump.exe :  \
	$(OBJ)/sntable_dump.o $(OBJ)/sntools.o  $(OBJ)/sntools_output.o 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/sntable_dump.o 	\
	$(OBJ)/sntools.o  	\
	$(OBJ_OUTPUT)		\
	$(LGSL) $(LCERN) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	(cd $(OBJ); rm sntable_dump.o )


# --------------------
# sntable_combine (Jan 2017)

$(OBJ)/sntable_combine.o : $(SRC)/sntable_combine.c $(SNTOOLS_OUTPUT)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) \
	$(SRC)/sntable_combine.c )

$(BIN)/sntable_combine.exe :  \
	$(OBJ)/sntable_combine.o $(OBJ)/sntools.o  $(OBJ)/sntools_output.o 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/sntable_combine.o 	\
	$(OBJ)/sntools.o  	\
	$(OBJ_OUTPUT)		\
	$(LGSL) $(LCERN) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	(cd $(OBJ); rm sntable_combine.o )

# --------------------
# nearnbr_maxFoM (Jun 2013)
# Link to libsnana.a to get forclose(LUN) wrapper

$(OBJ)/nearnbr_maxFoM.o : $(SRC)/nearnbr_maxFoM.c $(SNTOOLS_OUTPUT)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) \
	$(SRC)/nearnbr_maxFoM.c )

$(BIN)/nearnbr_maxFoM.exe :  \
	$(OBJ)/nearnbr_maxFoM.o $(OBJ)/sntools.o  $(OBJ)/sntools_output.o 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/nearnbr_maxFoM.o \
	$(OBJ)/sntools.o  	\
	$(OBJ_OUTPUT)		\
	$(LGSL) $(LCERN) $(LROOT) $(LCFITSIO)  -lm $(CPPLIB)
	(cd $(OBJ); rm nearnbr_maxFoM.o )

# --------------------
# nearnbr_apply (Jan 2017)

$(OBJ)/nearnbr_apply.o : $(SRC)/nearnbr_apply.c $(SNTOOLS_OUTPUT)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO)  \
	$(SRC)/nearnbr_apply.c )

$(BIN)/nearnbr_apply.exe :  \
	$(OBJ)/nearnbr_apply.o $(OBJ)/sntools.o  $(OBJ)/sntools_output.o  $(OBJ)/sntools_nearnbr.o 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/nearnbr_apply.o \
	$(OBJ)/sntools.o  	\
	$(OBJ)/sntools_gridmap.o \
	$(OBJ)/sntools_nearnbr.o \
	$(OBJ_OUTPUT)		\
	$(LGSL) $(LCERN) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	(cd $(OBJ); rm nearnbr_apply.o )


# -------------------------------------------------
# compile sntools (.o only)

$(OBJ)/sntools.o : $(SRC)/sntools.c $(SRC)/sntools.h $(SRC)/sndata.h 
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/sntools.c )

# -------------------------------------------------
# compile sntools_genSmear (.o only)
$(OBJ)/sntools_genSmear.o : $(SRC)/sntools_genSmear.c $(SRC)/sntools_genSmear.h $(SRC)/sntools.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/sntools_genSmear.c )

# -------------------------------------------------
# Apr 26 2016: non-linearity package

$(OBJ)/sntools_nonlinearity.o : $(SRC)/sntools_nonlinearity.c $(SRC)/sntools_nonlinearity.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_nonlinearity.c )

# -------------------------------------------------
# Oct 2015: fluxErr model

$(OBJ)/sntools_fluxErrModels.o : $(SRC)/sntools_fluxErrModels.c $(SRC)/sntools_fluxErrModels.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_fluxErrModels.c )

$(OBJ)/sntools_fluxErrModels_legacy.o : $(SRC)/sntools_fluxErrModels_legacy.c $(SRC)/sntools_fluxErrModels_legacy.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_fluxErrModels_legacy.c )

$(OBJ)/sntools_simeff_obsolete.o : $(SRC)/sntools_simeff_obsolete.c $(SRC)/sntools_simeff_obsolete.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_simeff_obsolete.c )

# -------------------------------------------------
# Oct 2020: cosmology util (extracted from sntools.c)

$(OBJ)/sntools_cosmology.o : $(SRC)/sntools_cosmology.c $(SRC)/sntools_cosmology.h $(SRC)/sntools_cosmology.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_cosmology.c )

# -------------------------------------------------
# compile sntools_dataformat_fits[text]

$(OBJ)/sntools_dataformat_fits.o : $(SRC)/sntools_dataformat_fits.c $(SRC)/sntools_dataformat_fits.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_dataformat_fits.c )

$(OBJ)/sntools_dataformat_text.o : $(SRC)/sntools_dataformat_text.c $(SRC)/sntools_dataformat_text.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_dataformat_text.c )

$(OBJ)/sntools_data.o : $(SRC)/sntools_data.c $(SRC)/sntools_data.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_data.c )

# ---------
$(OBJ)/sntools_calib.o : $(SRC)/sntools_calib.c $(SRC)/sntools_calib.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_calib.c )

$(OBJ)/sntools_modelgrid_gen.o : $(SRC)/sntools_modelgrid_gen.c $(SRC)/sntools_modelgrid.h $(SRC)/sndata.h $(SRC)/snlc_sim.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_modelgrid_gen.c )

$(OBJ)/sntools_modelgrid_read.o : $(SRC)/sntools_modelgrid_read.c $(SRC)/sntools_modelgrid.h $(SRC)/sndata.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_modelgrid_read.c )

# -------------------------------------------------
# compile sntools_host (.o only)

$(OBJ)/sntools_host.o : $(SRC)/sntools_host.c  $(SRC)/sntools_host.h $(SRC)/snlc_sim.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_host.c )

# -------------------------------------------------
# compile sntools_wgtmap (.o only)

$(OBJ)/sntools_wgtmap.o : $(SRC)/sntools_wgtmap.c  $(SRC)/sntools_wgtmap.h $(SRC)/snlc_sim.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_wgtmap.c )

# -------------------------------------------------
# compile sntools_devel (.o only)

$(OBJ)/sntools_devel.o : $(SRC)/sntools_devel.c  $(SRC)/sntools_devel.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_devel.c )

# -------------------------------------------------
# compile sntools_sim_readme (.o only)

$(OBJ)/sntools_sim_readme.o : $(SRC)/sntools_sim_readme.c  $(SRC)/sntools_sim_readme.h  $(SRC)/snlc_sim.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_sim_readme.c )

# --------------
# compile sntools_sim_atmosphere (.o only)
$(OBJ)/sntools_sim_atmosphere.o : $(SRC)/sntools_sim_atmosphere.c  $(SRC)/sntools_sim_atmosphere.h $(SRC)/snlc_sim.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_sim_atmosphere.c )

# -------------------------------------------------
# compile sntools_weaklens (.o only)

$(OBJ)/sntools_weaklens.o : $(SRC)/sntools_weaklens.c  $(SRC)/sntools_weaklens.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_weaklens.c )

# -------------------------------------------------
# compile sntools_stronglens (.o only)

$(OBJ)/sntools_stronglens.o : $(SRC)/sntools_stronglens.c  $(SRC)/sntools_stronglens.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_stronglens.c )

# -------------------------------------------------
# compile inoue_igm.o 

$(OBJ)/inoue_igm.o : $(SRC)/inoue_igm.c  $(SRC)/inoue_igm.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/inoue_igm.c )

# -------------------------------------------------
# compile sntools_wronghost (.o only)

$(OBJ)/sntools_wronghost.o : $(SRC)/sntools_wronghost.c  $(SRC)/sntools_wronghost.h $(SRC)/snlc_sim.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/sntools_wronghost.c )

# -------------------------------------------------
# compile sntools_wronghost (.o only)

$(OBJ)/sntools_zPDF_spline.o : $(SRC)/sntools_zPDF_spline.c $(SRC)/sntools_zPDF_spline.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_zPDF_spline.c )


# -------------------------------------------------
# compile sntools_trigger (.o only)
$(OBJ)/sntools_trigger.o : $(SRC)/sntools_trigger.c  $(SRC)/sntools_trigger.h  $(SRC)/snlc_sim.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_trigger.c )


# -------------------------------------------------
# compile sntools_spectrograph (.o only)
$(OBJ)/sntools_spectrograph.o : $(SRC)/sntools_spectrograph.c  $(SRC)/sntools_spectrograph.h  $(SRC)/kcor.h $(SRC)/kcor.c
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_spectrograph.c )

# -------------------------------------------------
# compile sntools_nearnbr (.o only)

$(OBJ)/sntools_nearnbr.o : $(SRC)/sntools_nearnbr.c $(SRC)/sntools_nearnbr.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(SRC)/sntools_nearnbr.c )

$(OBJ)/kdtree.o : $(SRC)/kdtree.c $(SRC)/kdtree.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(SRC)/kdtree.c )

# -------------------------------------------------
# compile sntools_output (.o only)

$(OBJ)/sntools_output.o : $(SNTOOLS_OUTPUT)
	(cd $(OBJ); \
	$(CCoutput) $(SNCFLAGS) $(IROOT) $(ICFITSIO) $(SRC)/sntools_output.c )

# -------------------------------------------------
#   MWgaldust from D.Cinabro (.o only)

$(OBJ)/MWgaldust.o : $(SRC)/MWgaldust.c 
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(SRC)/MWgaldust.c )

# -------------------------------------------------
#  multiseason toolkit (Oct 2014)

$(OBJ)/multiseason.o : $(SRC)/multiseason.c 
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(SRC)/multiseason.c )


# -------------------------------------------------
#  sntools_genPDF  tools (Jun 2020)

$(OBJ)/sntools_genPDF.o : $(SRC)/sntools_genPDF.c  $(SRC)/sntools_genPDF.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(IGSL) $(SRC)/sntools_genPDF.c )

# -------------------------------------------------
#  sntools_gridmap (July 2021)

$(OBJ)/sntools_gridmap.o : $(SRC)/sntools_gridmap.c  $(SRC)/sntools_gridmap.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(IGSL) $(SRC)/sntools_gridmap.c )
# -------------------------------------------------
#  sntools_genGauss_asym  tools (Sep 2016)

$(OBJ)/sntools_genGauss_asym.o : $(SRC)/sntools_genGauss_asym.c  $(SRC)/sntools_genGauss_asym.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(IGSL) $(SRC)/sntools_genGauss_asym.c )

# -------------------------------------------------        
#  sntools_genExpHalfGauss  tools (March 2020)                                         

$(OBJ)/sntools_genExpHalfGauss.o : $(SRC)/sntools_genExpHalfGauss.c  $(SRC)/sntools_genExpHalfGauss.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(IGSL) $(SRC)/sntools_genExpHalfGauss.c )

# ------------------------------------------------
# compile genmag_xxxx functions (.o only

$(OBJ)/genmag_stretch2.o : $(SRC)/genmag_stretch2.c $(SRC)/sntools.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(SRC)/genmag_stretch2.c )

$(OBJ)/genmag_mlcs2k2.o : $(SRC)/genmag_mlcs2k2.c  $(SRC)/sntools.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS)  $(SRC)/genmag_mlcs2k2.c )

$(OBJ)/genmag_SEDtools.o : $(SRC)/genmag_SEDtools.c  $(SRC)/genmag_SEDtools.c 
	(cd $(OBJ); \
	$(CC)  $(SNCFLAGS) $(SRC)/genmag_SEDtools.c )

$(OBJ)/genmag_SALT2.o : $(SRC)/genmag_SALT2.c  $(SRC)/genmag_SEDtools.c
	(cd $(OBJ); \
	$(CC)  $(SNCFLAGS) $(IGSL) $(SRC)/genmag_SALT2.c )

$(OBJ)/genmag_BAYESN.o : $(SRC)/genmag_BAYESN.c  $(SRC)/genmag_BAYESN.c
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_BAYESN.c )
#	$(CC) $(LYAML) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_BAYESN.c )

$(OBJ)/genmag_SIMSED.o : $(SRC)/genmag_SIMSED.c  $(SRC)/genmag_SEDtools.c 
	(cd $(OBJ); \
	$(CC)  $(SNCFLAGS) $(SRC)/genmag_SIMSED.c )

$(OBJ)/genmag_PySEDMODEL.o : $(SRC)/genmag_PySEDMODEL.c  $(SRC)/genmag_SEDtools.c 
	(cd $(OBJ); \
	$(CC)  $(SNCFLAGS) $(IPY) $(SRC)/genmag_PySEDMODEL.c )

$(OBJ)/genmag_snoopy.o : $(SRC)/genmag_snoopy.c  $(SRC)/genmag_snoopy.h $(SRC)/sntools.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_snoopy.c )

$(OBJ)/genmag_S11DM15.o : $(SRC)/genmag_S11DM15.c  $(SRC)/genmag_S11DM15.h $(SRC)/sntools.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(SRC)/genmag_S11DM15.c )

$(OBJ)/genmag_NON1ASED.o : $(SRC)/genmag_NON1ASED.c  $(SRC)/sntools.h  $(SRC)/snlc_sim.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_NON1ASED.c )

$(OBJ)/genmag_NON1AGRID.o : $(SRC)/genmag_NON1AGRID.c  $(SRC)/sntools.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL)  $(SRC)/genmag_NON1AGRID.c )

$(OBJ)/genmag_LCLIB.o : $(SRC)/genmag_LCLIB.c  $(SRC)/sntools.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL)  $(SRC)/genmag_LCLIB.c )

$(OBJ)/genmag_extinction.o : $(SRC)/genmag_extinction.c  $(SRC)/sntools.h $(SRC)/sntools_calib.c
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL)  $(SRC)/genmag_extinction.c )

$(OBJ)/genmag_extrap.o : $(SRC)/genmag_extrap.c  $(SRC)/sntools.h
	(cd $(OBJ); \
	$(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL)  $(SRC)/genmag_extrap.c )

# -------------------------------------------------
# Generic Lightcurve simulator 
# $(SRC)/sntools_genPDF.c

$(OBJ)/snlc_sim.o : $(SRC)/snlc_sim.c \
	$(SRC)/sntools.c $(SRC)/sntools_cosmology.c \
	$(SRC)/sntools_calib.c \
	$(SRC)/sntools_modelgrid_gen.c \
	$(SRC)/sntools_data.c 	\
	$(SRC)/sntools_dataformat_fits.c $(SRC)/sntools_dataformat_text.c \
	$(SRC)/sntools_host.c  \
	$(SRC)/sntools_genPDF.c $(SRC)/sntools_wronghost.c \
	$(SRC)/sntools_genSmear.c $(SRC)/sntools_devel.c \
	$(SRC)/SNcadenceFoM.c  \
	$(SRC)/sim_unit_tests.c \
	$(SRC)/snlc_sim.h $(SRC)/sndata.h 
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(IPY) $(SRC)/snlc_sim.c )

$(BIN)/snlc_sim.exe : $(OBJ_SIM)  $(LIB)/libsnana.a
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ_SIM)	  	\
	$(LIB)/libsnana.a  	\
	$(LCFITSIO) $(LGSL)  	\
	$(LPY)  $(LCERN) $(LROOT)  -lm  $(CPPLIB)
	(cd $(OBJ); rm snlc_sim.o)

# -------------------------------------------------
# program to compact the SIMLIB into one measure per night

$(OBJ)/simlib_coadd.o : $(SRC)/simlib_coadd.c $(SRC)/simlib_tools.c
	(cd $(OBJ); $(CC) $(SNCFLAGS)  $(SRC)/simlib_coadd.c )

$(BIN)/simlib_coadd.exe : \
	$(OBJ)/simlib_coadd.o $(OBJ)/MWgaldust.o
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/simlib_coadd.o  \
	$(OBJ)/sntools.o \
	$(OBJ)/sntools_output.o \
	$(OBJ)/sntools_cosmology.o \
	$(OBJ)/MWgaldust.o \
	$(LGSL) $(LCERN) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	(cd $(OBJ); rm simlib_coadd.o )

# ---------------------------------------
# merge_hbook.exe program (to replace hmerge.exe)
$(OBJ)/merge_hbook.o :  merge_hbook.F
	(cd $(OBJ); $(FFC) $(SNFFLAGS) $(SRC)/merge_hbook.F)

$(BIN)/merge_hbook.exe : $(OBJ)/merge_hbook.o 
	$(FFC) -o $@ $(SNLDFLAGS) $(OBJ)/merge_hbook.o  $(LCERN)
	(cd $(OBJ); rm merge_hbook.o )


# merge_root.exe program (Apr 2013)

$(OBJ)/merge_root.o : $(SRC)/merge_root.c $(SNTOOLS_OUTPUT) 
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) \
	$(SRC)/merge_root.c )

$(BIN)/merge_root.exe :  \
	$(OBJ)/merge_root.o $(OBJ)/sntools.o  $(OBJ)/sntools_output.o
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/merge_root.o 	\
	$(OBJ)/sntools.o  	\
	$(OBJ_OUTPUT)		\
	$(LGSL) $(LCERN) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	(cd $(OBJ); rm merge_root.o )

# -------------
# sim filter-calib

$(OBJ)/filtercal_sim.o : $(SRC)/filtercal_sim.c
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/filtercal_sim.c )

$(BIN)/filtercal_sim.exe : \
	$(OBJ)/filtercal_sim.o  $(OBJ)/sntools.o
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/filtercal_sim.o  \
	$(OBJ)/sntools.o \
	$(LCERN)  -lm $(LGSL)
	(cd $(OBJ); rm filtercal_sim.o )

# -------------
# SIMSED_fudge

$(OBJ)/SIMSED_fudge.o : $(SRC)/SIMSED_fudge.c $(SRC)/genmag_SIMSED.c 
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/SIMSED_fudge.c )

$(BIN)/SIMSED_fudge.exe : \
	$(OBJ)/SIMSED_fudge.o  $(OBJ)/sntools.o  $(OBJ)/sntools_genSmear.o \
	$(OBJ)/genmag_SEDtools.o  $(OBJ)/genmag_SIMSED.o 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/SIMSED_fudge.o  \
	$(OBJ)/sntools.o  \
	$(OBJ_SNTOOLS_DATA) 		\
	$(OBJ)/genmag_SEDtools.o \
	$(OBJ)/genmag_SIMSED.o   \
	$(OBJ)/genmag_extrap.o   \
	$(OBJ)/genmag_SALT2.o   \
	$(OBJ)/sntools_genSmear.o \
	$(OBJ)/sntools_devel.o \
	$(OBJ)/sntools_wgtmap.o \
	$(OBJ)/sntools_gridmap.o \
	$(OBJ)/MWgaldust.o \
	$(LIB)/libsnana.a  \
	-lm $(LCFITSIO) $(LGSL) $(LCERN) $(LROOT) $(CPPLIB)
	(cd $(OBJ); rm SIMSED_fudge.o )

# -------------
# SIMSED_extractSpec
# Extract a single interpolated SED for SIMSED version, 
# epoch, and parameter-list

$(OBJ)/SIMSED_extractSpec.o : $(SRC)/SIMSED_extractSpec.c \
	$(SRC)/genmag_SEDtools.c	
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/SIMSED_extractSpec.c )

$(BIN)/SIMSED_extractSpec.exe : \
	$(OBJ)/SIMSED_extractSpec.o  \
	$(OBJ)/sntools.o  \
	$(OBJ)/genmag_SEDtools.o 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/SIMSED_extractSpec.o  \
	$(OBJ)/sntools.o  \
	$(OBJ)/sntools_output.o  \
	$(OBJ)/sntools_cosmology.o  \
	$(OBJ)/MWgaldust.o  \
	$(OBJ)/genmag_SEDtools.o  \
	$(OBJ)/sntools_wgtmap.o \
	$(OBJ)/sntools_gridmap.o \
	-lm $(LCFITSIO) $(LGSL) $(LCERN) $(LROOT) $(CPPLIB)
	(cd $(OBJ); rm SIMSED_extractSpec.o )

# -------------
# SIMSED_rebin (May 31 2018)

$(OBJ)/SIMSED_rebin.o : $(SRC)/SIMSED_rebin.c \
	$(SRC)/genmag_SEDtools.c	
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/SIMSED_rebin.c )

$(BIN)/SIMSED_rebin.exe : \
	$(OBJ)/SIMSED_rebin.o  \
	$(OBJ)/sntools.o  $(OBJ)/genmag_SEDtools.o 
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/SIMSED_rebin.o  \
	$(OBJ)/sntools.o  \
	$(OBJ)/sntools_output.o  \
	$(OBJ)/sntools_cosmology.o  \
	$(OBJ)/MWgaldust.o  \
	$(OBJ)/genmag_SEDtools.o  \
	-lm $(LCFITSIO) $(LGSL) $(LCERN) $(LROOT) $(CPPLIB)
	(cd $(OBJ); rm SIMSED_rebin.o )

# -------------
# SIMSED_check

$(OBJ)/SIMSED_check.o : $(SRC)/SIMSED_check.c \
	$(SRC)/genmag_SEDtools.c
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/SIMSED_check.c )

$(BIN)/SIMSED_check.exe : \
	$(OBJ)/SIMSED_check.o  \
	$(OBJ)/sntools.o  $(OBJ)/genmag_SEDtools.o $(OBJ)/genmag_SIMSED.o
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/SIMSED_check.o  \
	$(OBJ)/sntools.o  \
	$(OBJ)/sntools_output.o  \
	$(OBJ)/sntools_cosmology.o  \
	$(OBJ)/sntools_wgtmap.o \
	$(OBJ)/sntools_gridmap.o \
	$(OBJ)/MWgaldust.o  \
	$(OBJ)/genmag_SEDtools.o  \
	$(OBJ)/genmag_SIMSED.o  \
	-lm $(LCFITSIO) $(LGSL) $(LCERN) $(LROOT) $(CPPLIB)
	(cd $(OBJ); rm SIMSED_check.o )

# -------------
# simualate effects of atmosphere

$(OBJ)/atmosphere_sim.o : $(SRC)/atmosphere_sim.c 
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/atmosphere_sim.c )

$(BIN)/atmosphere_sim.exe : \
	$(OBJ)/atmosphere_sim.o  $(OBJ)/sntools.o
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/atmosphere_sim.o	\
	$(OBJ)/sntools.o 	\
	$(OBJ_OUTPUT)		\
	-lm $(LCFITSIO) $(LGSL) $(LCERN) $(LROOT) $(CPPLIB)
	(cd $(OBJ); rm atmosphere_sim.o )

# -------------
# unfold_snpar

$(OBJ)/unfold_snpar.o : $(SRC)/unfold_snpar.c
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(SRC)/unfold_snpar.c )

$(BIN)/unfold_snpar.exe : $(OBJ)/unfold_snpar.o  $(OBJ)/sntools.o $(OBJ)/sntools_output.o
	$(FFC) -o $@ $(SNLDFLAGS) \
	$(OBJ)/unfold_snpar.o	\
	$(OBJ)/sntools.o 	\
	$(OBJ_OUTPUT)		\
	$(LGSL) $(LCERN) $(LROOT) $(LCFITSIO) -lm $(CPPLIB) 
	(cd $(OBJ); rm unfold_snpar.o )

####################  END ##########################################
