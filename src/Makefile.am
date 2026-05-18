# ===========================================================================
# Makefile.am for SNANA
# Autotools rewrite of the legacy src/Makefile.
#
# Design:
#   - configure.ac handles all library/flag detection
#   - Object files land in $(top_builddir)/obj  (mirrors original ../obj)
#   - Binaries land in $(top_builddir)/bin       (mirrors original ../bin)
#   - Rules are explicit 1:1 translations of the original Makefile
#   - Single-executable builds supported: make snlc_sim, make wfit, etc.
# ===========================================================================

# ---------------------------------------------------------------------------
# Directory layout — mirrors original Makefile's OBJ=../obj, BIN=../bin
# ---------------------------------------------------------------------------

SRC  = $(top_srcdir)/src
OBJ  = $(top_srcdir)/obj
BIN  = $(top_srcdir)/bin
UTIL = $(top_srcdir)/util

# tack on flag to provide git version to C/fortran codes
AM_CFLAGS = -DGIT_SNANA_VERSION=\"@GIT_VERSION@\"

EXTRA_FLAGS = -D_FORTIFY_SOURCE=2


# ---------------------------------------------------------------------------
# Compiler flags — translated from original SNFFLAGS / SNCFLAGS / SNLDFLAGS
# These are set by configure.ac and passed in via AM_ variables; we
# assemble the final flags here to match the original Makefile exactly.
# ---------------------------------------------------------------------------
SNFFLAGS  = -c -fno-automatic -m64 -fsecond-underscore $(AM_FCFLAGS) $(FCFLAGS) $(DEBUG_FCFLAGS)
SNCFLAGS  = -c -m64 $(AM_CFLAGS) $(CFLAGS) -fcommon $(DEBUG_CFLAGS)
SNCXXFLAGS = -c -m64 $(AM_CXXFLAGS) $(CXXFLAGS) -fcommon $(DEBUG_CXXFLAGS)
SNLDFLAGS = -m64 $(PROFILEFLAG) $(DEBUG_LDFLAGS)

# Python includes/libs (set by configure.ac via @PYTHON_INCLUDES@ / @PYTHON_LIBS@)
IPY = @PYTHON_INCLUDES@
LPY = @PYTHON_LIBS@

# Library flags (set by configure.ac)
LCFITSIO  = $(CFITSIOLIB)
LGSL      = $(GSLLIB) $(GSLCBLASLIB)

# ROOT (conditionally set by configure.ac)
if USE_ROOT
LROOT  = -L$(ROOT_DIR)/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d \
          -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics \
          -lMathCore -lThread -ldl
IROOT  = -I$(ROOT_DIR)/include
else
LROOT  =
IROOT  =
endif

CPPLIB = -lstdc++

# ---------------------------------------------------------------------------
# Create obj/ and bin/ before anything is compiled
# clean-local removes object files and binaries
# ---------------------------------------------------------------------------
all-local: set_Cpreproc_flags $(SNANA_EXECUTABLES)  list-all
	$(MKDIR_P) $(OBJ) $(BIN)

clean-local:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod
	rm -f $(BIN)/*.exe

# alter define flags in a few C-include files based on existance of
# ROOT and GSL_DIR version
set_Cpreproc_flags:
	$(UTIL)/set_Cpreproc_flags.py

# maybe somebody, can figure out how to implement this list after make-all completes
list-all:
	ls -l $(BIN)

# ---------------------------------------------------------------------------
# Tell automake we are managing binaries ourselves.
# Leaving bin_PROGRAMS empty prevents automake from generating
# its own compile/link rules that would conflict with ours.
# ---------------------------------------------------------------------------
bin_PROGRAMS =
noinst_PROGRAMS =

# ---------------------------------------------------------------------------
# Object-file groups — mirrors original OBJ_GENMAG_FIT, OBJ_SIM, etc.
# ---------------------------------------------------------------------------

OBJ_GENMAG_FIT = \
	$(OBJ)/genmag_stretch2.o \
	$(OBJ)/genmag_mlcs2k2.o \
	$(OBJ)/genmag_SALT2.o \
	$(OBJ)/genmag_BAYESN.o \
	$(OBJ)/genmag_snoopy.o \
	$(OBJ)/genmag_S11DM15.o \
	$(OBJ)/genmag_SEDtools.o \
	$(OBJ)/genmag_extinction.o \
	$(OBJ)/genmag_extrap.o

OBJ_GENMAG_SIM = \
	$(OBJ_GENMAG_FIT) \
	$(OBJ)/genmag_NON1ASED.o \
	$(OBJ)/genmag_NON1AGRID.o \
	$(OBJ)/genmag_LCLIB.o \
	$(OBJ)/genmag_PySEDMODEL.o \
	$(OBJ)/genmag_SIMSED.o

OBJ_SNTOOLS_SNANA = \
	$(OBJ)/sntools.o \
	$(OBJ)/sntools_gridmap.o \
	$(OBJ)/sntools_nonlinearity.o \
	$(OBJ)/sntools_fluxErrModels.o \
	$(OBJ)/sntools_fluxErrModels_legacy.o \
	$(OBJ)/sntools_simeff_obsolete.o \
	$(OBJ)/sntools_dataformat_fits.o \
	$(OBJ)/sntools_dataformat_text.o \
	$(OBJ)/sntools_data.o \
	$(OBJ)/sntools_zPDF_spline.o \
	$(OBJ)/sntools_spline.o \
	$(OBJ)/sntools_modelgrid_read.o \
	$(OBJ)/sntools_calib.o \
	$(OBJ)/sntools_output.o \
	$(OBJ)/yaml.o \
	$(OBJ)/MWgaldust.o \
	$(OBJ)/sntools_cosmology.o \
	$(OBJ)/multiseason.o

SNTOOLS_DATA_OBJ = \
	$(OBJ)/sntools_dataformat_fits.o \
	$(OBJ)/sntools_dataformat_text.o \
	$(OBJ)/sntools_data.o \
	$(OBJ)/sntools_output.o \
	$(OBJ)/sntools_cosmology.o

OBJ_OUTPUT = $(OBJ)/sntools_output.o

OBJ_SNANA = \
	$(OBJ)/snana.o \
	$(OBJ_SNTOOLS_SNANA) \
	$(OBJ)/minuit.o

OBJ_SNFIT = \
	$(OBJ)/snlc_fit.o \
	$(OBJ_SNTOOLS_SNANA) \
	$(OBJ)/sntools_genSmear.o \
	$(OBJ)/sntools_devel.o \
	$(OBJ)/minuit.o \
	$(OBJ_GENMAG_FIT)

OBJ_PSNID = \
	$(OBJ)/psnid.o \
	$(OBJ_SNTOOLS_SNANA) \
	$(OBJ)/psnid_tools.o \
	$(OBJ)/psnid_BEST.o \
	$(OBJ)/psnid_BEST2.o \
	$(OBJ)/minuit.o

OBJ_SIM = \
	$(OBJ)/snlc_sim.o \
	$(OBJ_SNTOOLS_SNANA) \
	$(OBJ)/sntools_host.o \
	$(OBJ)/sntools_wgtmap.o \
	$(OBJ)/sntools_weaklens.o \
	$(OBJ)/sntools_stronglens.o \
	$(OBJ)/sntools_wronghost.o \
	$(OBJ)/sntools_genSmear.o \
	$(OBJ)/sntools_devel.o \
	$(OBJ)/sntools_sim_readme.o \
	$(OBJ)/sntools_sim_atmosphere.o \
	$(OBJ)/sntools_trigger.o \
	$(OBJ)/sntools_modelgrid_gen.o \
	$(OBJ)/sntools_spectrograph.o \
	$(OBJ)/sntools_genPDF.o \
	$(OBJ)/sntools_genGauss_asym.o \
	$(OBJ)/sntools_genExpHalfGauss.o \
	$(OBJ)/inoue_igm.o \
	$(OBJ_GENMAG_SIM)

# ---------------------------------------------------------------------------
# Full list of executables (used by all-local and .PHONY aliases)
# ---------------------------------------------------------------------------
SNANA_EXECUTABLES = \
	$(BIN)/snlc_sim.exe \
	$(BIN)/snana.exe \
	$(BIN)/snlc_fit.exe \
	$(BIN)/psnid.exe \
	$(BIN)/kcor.exe \
	$(BIN)/wfit.exe \
	$(BIN)/SALT2mu.exe \
	$(BIN)/combine_fitres.exe \
	$(BIN)/sntable_dump.exe \
	$(BIN)/sntable_combine.exe \
	$(BIN)/SIMSED_extractSpec.exe \
	$(BIN)/SIMSED_check.exe \
	$(BIN)/SIMSED_rebin.exe \
	$(BIN)/simlib_coadd.exe

# ---------------------------------------------------------------------------
# Short-name phony aliases for single-executable builds.
# Usage: make snlc_sim  OR  make wfit  etc.
# ---------------------------------------------------------------------------
.PHONY: snlc_sim snana snlc_fit psnid kcor wfit SALT2mu \
        combine_fitres sntable_dump sntable_combine \
        SIMSED_extractSpec SIMSED_check SIMSED_rebin simlib_coadd

snlc_sim:        $(BIN)/snlc_sim.exe
snana:           $(BIN)/snana.exe
snlc_fit:        $(BIN)/snlc_fit.exe
psnid:           $(BIN)/psnid.exe
kcor:            $(BIN)/kcor.exe
wfit:            $(BIN)/wfit.exe
SALT2mu:         $(BIN)/SALT2mu.exe
combine_fitres:  $(BIN)/combine_fitres.exe
sntable_dump:    $(BIN)/sntable_dump.exe
sntable_combine: $(BIN)/sntable_combine.exe
SIMSED_extractSpec: $(BIN)/SIMSED_extractSpec.exe
SIMSED_check:    $(BIN)/SIMSED_check.exe
SIMSED_rebin:    $(BIN)/SIMSED_rebin.exe
simlib_coadd:    $(BIN)/simlib_coadd.exe

# ===========================================================================
# OBJECT FILE RULES
# Each rule mirrors the original Makefile exactly:
#   (cd $(OBJ); compiler flags $(SRC)/file.c)
# This places the .o file in $(OBJ) naturally since the compiler
# writes output to the current directory when no -o is specified.
# ===========================================================================

# ---------------------------------------------------------------------------
# Fortran
# ---------------------------------------------------------------------------

$(OBJ)/minuit.o: $(SRC)/minuit.F
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(FC) $(SNFFLAGS) $(SRC)/minuit.F)

$(OBJ)/snana.o: $(SRC)/snana.F90
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(FC) $(SNFFLAGS) -DSNANA -DMINUIT $(SRC)/snana.F90)

$(OBJ)/snlc_fit.o: $(SRC)/snlc_fit.F90
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(FC) $(SNFFLAGS) -DSNFIT -DMINUIT $(SRC)/snlc_fit.F90)

$(OBJ)/psnid.o: $(SRC)/psnid.F90
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(FC) $(SNFFLAGS) -DPSNID $(SRC)/psnid.F90)

# ---------------------------------------------------------------------------
# C++ 
# ---------------------------------------------------------------------------

$(OBJ)/sntools_output.o: $(SRC)/sntools_output.cpp $(SRC)/sntools_output.h \
	$(SRC)/sntools_output_marz.c $(SRC)/sntools_output_text.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CXX) $(SNCXXFLAGS) $(IROOT) $(ICFITSIO) $(SRC)/sntools_output.cpp)

$(OBJ)/sntools_npz.o: $(SRC)/sntools_npz.cpp
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CXX) $(SNCXXFLAGS) -I$(CONDA_PREFIX)/include $(SRC)/sntools_npz.cpp)

$(OBJ)/cnpy.o: $(SRC)/cnpy.cpp $(SRC)/cnpy.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CXX) $(SNCXXFLAGS) -I$(CONDA_PREFIX)/include $(SRC)/cnpy.cpp)

# ---------------------------------------------------------------------------
# C — core sntools
# ---------------------------------------------------------------------------

$(OBJ)/sntools.o: $(SRC)/sntools.c $(SRC)/sntools.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/sntools.c)

$(OBJ)/yaml.o: $(SRC)/yaml.c $(SRC)/yaml.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/yaml.c)

$(OBJ)/sntools_cosmology.o: $(SRC)/sntools_cosmology.c $(SRC)/sntools_cosmology.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_cosmology.c)

$(OBJ)/sntools_gridmap.o: $(SRC)/sntools_gridmap.c $(SRC)/sntools_gridmap.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_gridmap.c)

$(OBJ)/sntools_nonlinearity.o: $(SRC)/sntools_nonlinearity.c \
	$(SRC)/sntools_nonlinearity.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_nonlinearity.c)

$(OBJ)/sntools_fluxErrModels.o: $(SRC)/sntools_fluxErrModels.c \
	$(SRC)/sntools_fluxErrModels.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_fluxErrModels.c)

$(OBJ)/sntools_fluxErrModels_legacy.o: $(SRC)/sntools_fluxErrModels_legacy.c \
	$(SRC)/sntools_fluxErrModels_legacy.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_fluxErrModels_legacy.c)

$(OBJ)/sntools_simeff_obsolete.o: $(SRC)/sntools_simeff_obsolete.c \
	$(SRC)/sntools_simeff_obsolete.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_simeff_obsolete.c)

$(OBJ)/sntools_dataformat_fits.o: $(SRC)/sntools_dataformat_fits.c \
	$(SRC)/sntools_dataformat_fits.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_dataformat_fits.c)

$(OBJ)/sntools_dataformat_text.o: $(SRC)/sntools_dataformat_text.c \
	$(SRC)/sntools_dataformat_text.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_dataformat_text.c)

$(OBJ)/sntools_data.o: $(SRC)/sntools_data.c $(SRC)/sntools_data.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_data.c)

$(OBJ)/sntools_calib.o: $(SRC)/sntools_calib.c $(SRC)/sntools_calib.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_calib.c)

$(OBJ)/sntools_modelgrid_read.o: $(SRC)/sntools_modelgrid_read.c \
	$(SRC)/sntools_modelgrid.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_modelgrid_read.c)

$(OBJ)/sntools_modelgrid_gen.o: $(SRC)/sntools_modelgrid_gen.c \
	$(SRC)/sntools_modelgrid.h $(SRC)/sndata.h $(SRC)/snlc_sim.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_modelgrid_gen.c)

$(OBJ)/sntools_zPDF_spline.o: $(SRC)/sntools_zPDF_spline.c $(SRC)/sntools_zPDF_spline.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_zPDF_spline.c)

$(OBJ)/sntools_spline.o: $(SRC)/sntools_spline.c $(SRC)/sntools_spline.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_spline.c)

$(OBJ)/MWgaldust.o: $(SRC)/MWgaldust.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/MWgaldust.c)

$(OBJ)/multiseason.o: $(SRC)/multiseason.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/multiseason.c)

# ---------------------------------------------------------------------------
# C — simulation tools
# ---------------------------------------------------------------------------

$(OBJ)/snlc_sim.o: $(SRC)/snlc_sim.c $(SRC)/sntools.c \
	$(SRC)/sntools_cosmology.c $(SRC)/sntools_calib.c \
	$(SRC)/sntools_modelgrid_gen.c $(SRC)/sntools_data.c \
	$(SRC)/sntools_dataformat_fits.c $(SRC)/sntools_dataformat_text.c \
	$(SRC)/sntools_host.c $(SRC)/sntools_genPDF.c \
	$(SRC)/sntools_wronghost.c $(SRC)/sntools_genSmear.c \
	$(SRC)/sntools_devel.c $(SRC)/snlc_sim.h $(SRC)/sndata.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(IPY) $(SRC)/snlc_sim.c)

$(OBJ)/sntools_host.o: $(SRC)/sntools_host.c $(SRC)/sntools_host.h $(SRC)/snlc_sim.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_host.c)

$(OBJ)/sntools_wgtmap.o: $(SRC)/sntools_wgtmap.c $(SRC)/sntools_wgtmap.h $(SRC)/snlc_sim.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_wgtmap.c)

$(OBJ)/sntools_weaklens.o: $(SRC)/sntools_weaklens.c $(SRC)/sntools_weaklens.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_weaklens.c)

$(OBJ)/sntools_stronglens.o: $(SRC)/sntools_stronglens.c $(SRC)/sntools_stronglens.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_stronglens.c)

$(OBJ)/sntools_wronghost.o: $(SRC)/sntools_wronghost.c \
	$(SRC)/sntools_wronghost.h $(SRC)/snlc_sim.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_wronghost.c)

$(OBJ)/sntools_genSmear.o: $(SRC)/sntools_genSmear.c \
	$(SRC)/sntools_genSmear.h $(SRC)/sntools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/sntools_genSmear.c)

$(OBJ)/sntools_devel.o: $(SRC)/sntools_devel.c $(SRC)/sntools_devel.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_devel.c)

$(OBJ)/sntools_sim_readme.o: $(SRC)/sntools_sim_readme.c \
	$(SRC)/sntools_sim_readme.h $(SRC)/snlc_sim.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_sim_readme.c)

$(OBJ)/sntools_sim_atmosphere.o: $(SRC)/sntools_sim_atmosphere.c \
	$(SRC)/sntools_sim_atmosphere.h $(SRC)/snlc_sim.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/sntools_sim_atmosphere.c)

$(OBJ)/sntools_trigger.o: $(SRC)/sntools_trigger.c \
	$(SRC)/sntools_trigger.h $(SRC)/snlc_sim.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_trigger.c)

$(OBJ)/sntools_spectrograph.o: $(SRC)/sntools_spectrograph.c \
	$(SRC)/sntools_spectrograph.h $(SRC)/kcor.h $(SRC)/kcor.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(SRC)/sntools_spectrograph.c)

$(OBJ)/sntools_genPDF.o: $(SRC)/sntools_genPDF.c $(SRC)/sntools_genPDF.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_genPDF.c)

$(OBJ)/sntools_genGauss_asym.o: $(SRC)/sntools_genGauss_asym.c \
	$(SRC)/sntools_genGauss_asym.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_genGauss_asym.c)

$(OBJ)/sntools_genExpHalfGauss.o: $(SRC)/sntools_genExpHalfGauss.c \
	$(SRC)/sntools_genExpHalfGauss.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/sntools_genExpHalfGauss.c)

$(OBJ)/inoue_igm.o: $(SRC)/inoue_igm.c $(SRC)/inoue_igm.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/inoue_igm.c)

# ---------------------------------------------------------------------------
# C — genmag models
# ---------------------------------------------------------------------------

$(OBJ)/genmag_SEDtools.o: $(SRC)/genmag_SEDtools.c $(SRC)/genmag_SEDtools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/genmag_SEDtools.c)

$(OBJ)/genmag_stretch2.o: $(SRC)/genmag_stretch2.c $(SRC)/sntools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/genmag_stretch2.c)

$(OBJ)/genmag_mlcs2k2.o: $(SRC)/genmag_mlcs2k2.c $(SRC)/sntools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/genmag_mlcs2k2.c)

$(OBJ)/genmag_SALT2.o: $(SRC)/genmag_SALT2.c $(SRC)/genmag_SEDtools.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/genmag_SALT2.c)

$(OBJ)/genmag_BAYESN.o: $(SRC)/genmag_BAYESN.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_BAYESN.c)

$(OBJ)/genmag_snoopy.o: $(SRC)/genmag_snoopy.c $(SRC)/genmag_snoopy.h $(SRC)/sntools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_snoopy.c)

$(OBJ)/genmag_S11DM15.o: $(SRC)/genmag_S11DM15.c $(SRC)/genmag_S11DM15.h $(SRC)/sntools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/genmag_S11DM15.c)

$(OBJ)/genmag_NON1ASED.o: $(SRC)/genmag_NON1ASED.c $(SRC)/sntools.h $(SRC)/snlc_sim.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_NON1ASED.c)

$(OBJ)/genmag_NON1AGRID.o: $(SRC)/genmag_NON1AGRID.c $(SRC)/sntools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_NON1AGRID.c)

$(OBJ)/genmag_LCLIB.o: $(SRC)/genmag_LCLIB.c $(SRC)/sntools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_LCLIB.c)

$(OBJ)/genmag_SIMSED.o: $(SRC)/genmag_SIMSED.c $(SRC)/genmag_SEDtools.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/genmag_SIMSED.c)

$(OBJ)/genmag_PySEDMODEL.o: $(SRC)/genmag_PySEDMODEL.c $(SRC)/genmag_SEDtools.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IPY) $(SRC)/genmag_PySEDMODEL.c)

$(OBJ)/genmag_extinction.o: $(SRC)/genmag_extinction.c \
	$(SRC)/sntools.h $(SRC)/sntools_calib.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_extinction.c)

$(OBJ)/genmag_extrap.o: $(SRC)/genmag_extrap.c $(SRC)/sntools.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/genmag_extrap.c)

# ---------------------------------------------------------------------------
# C — psnid tools
# ---------------------------------------------------------------------------

$(OBJ)/psnid_tools.o: $(SRC)/psnid_tools.c $(SRC)/sntools.h \
	$(SRC)/sntools_modelgrid_read.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/psnid_tools.c)

$(OBJ)/psnid_BEST.o: $(SRC)/psnid_BEST.c $(SRC)/sntools.h \
	$(SRC)/sntools_modelgrid_read.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/psnid_BEST.c)

$(OBJ)/psnid_BEST2.o: $(SRC)/psnid_BEST2.c $(SRC)/sntools.h \
	$(SRC)/sntools_modelgrid_read.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/psnid_BEST2.c)

# ---------------------------------------------------------------------------
# C — standalone program objects
# ---------------------------------------------------------------------------

$(OBJ)/wfit.o: $(SRC)/wfit.c $(SRC)/sntools.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/wfit.c)

$(OBJ)/kcor.o: $(SRC)/kcor.c $(SRC)/kcor.h $(SRC)/sntools_spectrograph.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(ICFITSIO) $(IGSL) $(SRC)/kcor.c)

$(OBJ)/SALT2mu.o: $(SRC)/SALT2mu.c $(SRC)/sntools.c \
	$(SRC)/sntools_output.cpp $(SRC)/sntools_genPDF.c \
	$(SRC)/sntools_genGauss_asym.c $(SRC)/minuit.F
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(SRC)/SALT2mu.c)

$(OBJ)/combine_fitres.o: $(SRC)/combine_fitres.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(SRC)/combine_fitres.c)

$(OBJ)/sntable_dump.o: $(SRC)/sntable_dump.c $(SRC)/sntools_output.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(SRC)/sntable_dump.c)

$(OBJ)/sntable_combine.o: $(SRC)/sntable_combine.c $(SRC)/sntools_output.h
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(ICFITSIO) $(SRC)/sntable_combine.c)

$(OBJ)/SIMSED_extractSpec.o: $(SRC)/SIMSED_extractSpec.c $(SRC)/genmag_SEDtools.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(IGSL) $(SRC)/SIMSED_extractSpec.c)

$(OBJ)/SIMSED_check.o: $(SRC)/SIMSED_check.c $(SRC)/genmag_SEDtools.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/SIMSED_check.c)

$(OBJ)/SIMSED_rebin.o: $(SRC)/SIMSED_rebin.c $(SRC)/genmag_SEDtools.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/SIMSED_rebin.c)

$(OBJ)/simlib_coadd.o: $(SRC)/simlib_coadd.c $(SRC)/simlib_tools.c
	$(MKDIR_P) $(OBJ)
	(cd $(OBJ); $(CC) $(SNCFLAGS) $(SRC)/simlib_coadd.c)

# ===========================================================================
# BINARY LINK RULES
# Each rule mirrors the original Makefile exactly.
# The rm of the primary .o after linking matches the original behavior
# (shared .o files like sntools.o are NOT removed; only the
# program-specific .o is removed to signal "this exe is up to date").
# ===========================================================================

$(BIN)/snana.exe: $(OBJ_SNANA)
	$(MKDIR_P) $(BIN)
	$(FC) -o $@ $(SNLDFLAGS) \
		$(OBJ_SNANA) \
		$(LCFITSIO) $(LGSL) \
		$(LROOT) -lm $(CPPLIB)
	cd $(OBJ); rm -f snana.o; rm -f *.mod

$(BIN)/snlc_fit.exe: $(OBJ_SNFIT)
	$(MKDIR_P) $(BIN)
	$(FC) -o $@ $(SNLDFLAGS) \
		$(OBJ_SNFIT) \
		$(LCFITSIO) $(LGSL) \
		$(LROOT) -lm $(CPPLIB)
	cd $(OBJ); rm -f snlc_fit.o; rm -f *.mod

$(BIN)/psnid.exe: $(OBJ_PSNID)
	$(MKDIR_P) $(BIN)
	$(FC) -o $@ $(SNLDFLAGS) \
		$(OBJ_PSNID) \
		$(LCFITSIO) $(LGSL) \
		$(LROOT) -lm $(CPPLIB)
	cd $(OBJ); rm -f psnid.o; rm -f *.mod

$(BIN)/snlc_sim.exe: $(OBJ_SIM)
	$(MKDIR_P) $(BIN)
	$(CC) -o $@ $(SNLDFLAGS) \
		$(OBJ_SIM) \
		$(LCFITSIO) $(LGSL) \
		$(LPY) $(LROOT) -lm $(CPPLIB)
	cd $(OBJ); rm -f snlc_sim.o

$(BIN)/wfit.exe: $(OBJ)/wfit.o $(OBJ)/cnpy.o $(OBJ)/sntools.o \
	$(OBJ)/sntools_npz.o $(OBJ)/sntools_output.o $(OBJ)/sntools_cosmology.o
	$(MKDIR_P) $(BIN)
	$(CC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/wfit.o \
		$(OBJ)/sntools.o \
		$(OBJ)/sntools_output.o \
		$(OBJ)/sntools_npz.o \
		$(OBJ)/cnpy.o \
		$(OBJ)/sntools_cosmology.o \
		$(LROOT) $(CPPLIB) \
		$(LCFITSIO) -lm $(LGSL) -lz
	cd $(OBJ); rm -f wfit.o

$(BIN)/kcor.exe: $(OBJ)/kcor.o $(OBJ)/sntools.o $(OBJ)/MWgaldust.o \
	$(OBJ)/sntools_spectrograph.o $(OBJ)/genmag_SEDtools.o \
	$(SNTOOLS_DATA_OBJ)
	$(MKDIR_P) $(BIN)
	$(CC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/kcor.o \
		$(OBJ)/sntools.o \
		$(OBJ)/MWgaldust.o \
		$(SNTOOLS_DATA_OBJ) \
		$(OBJ)/sntools_spectrograph.o \
		$(OBJ)/genmag_SEDtools.o \
		$(LCFITSIO) $(LGSL) $(LROOT) -lm $(CPPLIB)
	cd $(OBJ); rm -f kcor.o

$(BIN)/SALT2mu.exe: $(OBJ)/SALT2mu.o $(OBJ)/sntools.o \
	$(OBJ)/sntools_output.o $(OBJ)/minuit.o \
	$(OBJ)/sntools_gridmap.o $(OBJ)/sntools_genGauss_asym.o \
	$(OBJ)/sntools_genExpHalfGauss.o
	$(MKDIR_P) $(BIN)
	$(FC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/SALT2mu.o \
		$(OBJ)/sntools.o \
		$(OBJ)/sntools_output.o \
		$(OBJ)/sntools_gridmap.o \
		$(OBJ)/sntools_genGauss_asym.o \
		$(OBJ)/sntools_genExpHalfGauss.o \
		$(OBJ)/minuit.o \
		$(LROOT) -lm $(LGSL) $(LCFITSIO) $(CPPLIB) -lpthread
	cd $(OBJ); rm -f SALT2mu.o

$(BIN)/combine_fitres.exe: $(OBJ)/combine_fitres.o \
	$(OBJ)/sntools.o $(OBJ)/sntools_output.o
	$(MKDIR_P) $(BIN)
	$(CC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/combine_fitres.o \
		$(OBJ)/sntools.o \
		$(OBJ_OUTPUT) \
		$(LGSL) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	cd $(OBJ); rm -f combine_fitres.o

$(BIN)/sntable_dump.exe: $(OBJ)/sntable_dump.o \
	$(OBJ)/sntools.o $(OBJ)/sntools_output.o
	$(MKDIR_P) $(BIN)
	$(CC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/sntable_dump.o \
		$(OBJ)/sntools.o \
		$(OBJ_OUTPUT) \
		$(LGSL) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	cd $(OBJ); rm -f sntable_dump.o

$(BIN)/sntable_combine.exe: $(OBJ)/sntable_combine.o \
	$(OBJ)/sntools.o $(OBJ)/sntools_output.o
	$(MKDIR_P) $(BIN)
	$(CC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/sntable_combine.o \
		$(OBJ)/sntools.o \
		$(OBJ_OUTPUT) \
		$(LGSL) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	cd $(OBJ); rm -f sntable_combine.o

$(BIN)/SIMSED_extractSpec.exe: $(OBJ)/SIMSED_extractSpec.o \
	$(OBJ)/sntools.o $(OBJ)/genmag_SEDtools.o \
	$(OBJ)/sntools_output.o $(OBJ)/sntools_cosmology.o \
	$(OBJ)/MWgaldust.o $(OBJ)/sntools_wgtmap.o $(OBJ)/sntools_gridmap.o
	$(MKDIR_P) $(BIN)
	$(FC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/SIMSED_extractSpec.o \
		$(OBJ)/sntools.o \
		$(OBJ)/sntools_output.o \
		$(OBJ)/sntools_cosmology.o \
		$(OBJ)/MWgaldust.o \
		$(OBJ)/genmag_SEDtools.o \
		$(OBJ)/sntools_wgtmap.o \
		$(OBJ)/sntools_gridmap.o \
		-lm $(LCFITSIO) $(LGSL) $(LROOT) $(CPPLIB)
	cd $(OBJ); rm -f SIMSED_extractSpec.o

$(BIN)/SIMSED_check.exe: $(OBJ)/SIMSED_check.o \
	$(OBJ)/sntools.o $(OBJ)/genmag_SEDtools.o $(OBJ)/genmag_SIMSED.o \
	$(OBJ)/sntools_output.o $(OBJ)/sntools_cosmology.o \
	$(OBJ)/sntools_wgtmap.o $(OBJ)/sntools_gridmap.o $(OBJ)/MWgaldust.o
	$(MKDIR_P) $(BIN)
	$(FC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/SIMSED_check.o \
		$(OBJ)/sntools.o \
		$(OBJ)/sntools_output.o \
		$(OBJ)/sntools_cosmology.o \
		$(OBJ)/sntools_wgtmap.o \
		$(OBJ)/sntools_gridmap.o \
		$(OBJ)/MWgaldust.o \
		$(OBJ)/genmag_SEDtools.o \
		$(OBJ)/genmag_SIMSED.o \
		-lm $(LCFITSIO) $(LGSL) $(LROOT) $(CPPLIB)
	cd $(OBJ); rm -f SIMSED_check.o

$(BIN)/SIMSED_rebin.exe: $(OBJ)/SIMSED_rebin.o \
	$(OBJ)/sntools.o $(OBJ)/genmag_SEDtools.o \
	$(OBJ)/sntools_output.o $(OBJ)/sntools_cosmology.o $(OBJ)/MWgaldust.o
	$(MKDIR_P) $(BIN)
	$(FC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/SIMSED_rebin.o \
		$(OBJ)/sntools.o \
		$(OBJ)/sntools_output.o \
		$(OBJ)/sntools_cosmology.o \
		$(OBJ)/MWgaldust.o \
		$(OBJ)/genmag_SEDtools.o \
		-lm $(LCFITSIO) $(LGSL) $(LROOT) $(CPPLIB)
	cd $(OBJ); rm -f SIMSED_rebin.o

$(BIN)/simlib_coadd.exe: $(OBJ)/simlib_coadd.o $(OBJ)/MWgaldust.o \
	$(OBJ)/sntools.o $(OBJ)/sntools_output.o $(OBJ)/sntools_cosmology.o
	$(MKDIR_P) $(BIN)
	$(CC) -o $@ $(SNLDFLAGS) \
		$(OBJ)/simlib_coadd.o \
		$(OBJ)/sntools.o \
		$(OBJ)/sntools_output.o \
		$(OBJ)/sntools_cosmology.o \
		$(OBJ)/MWgaldust.o \
		$(LGSL) $(LROOT) $(LCFITSIO) -lm $(CPPLIB)
	cd $(OBJ); rm -f simlib_coadd.o

# ===========================================================================
BUILT_SOURCES =
EXTRA_DIST =
CLEANFILES =
