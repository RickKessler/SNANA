
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

bin_PROGRAMS = sntable_dump sntable_combine combine_fitres.exe SALT2mu wfit.exe kcor.exe snlc_sim.exe SIMSED_fudge.exe SIMSED_extractSpec.exe SIMSED_check.exe SIMSED_rebin.exe simlib_coadd.exe unfold_snpar.exe

sntable_dump_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm $(GSLCBLASLIB) -lpthread -lz -lstdc++ $(ROOTLIBS)

sntable_combine_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm $(GSLCBLASLIB) -lpthread -lz -lstdc++ $(ROOTLIBS)

combine_fitres_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm $(GSLCBLASLIB) -lpthread -lz -lstdc++ $(ROOTLIBS)

wfit_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread -lz -lstdc++ $(ROOTLIBS)

kcor_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread -lz -lstdc++ $(ROOTLIBS)

snlc_sim_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread $(ROOTLIBS) @PYTHON_LIBS@

snlc_sim_exe_CFLAGS = $(AM_CFLAGS) @PYTHON_INCLUDES@

#snlc_fit_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread $(ROOTLIBS) $(FLIBS)

#snlc_fit_exe_LINK = $(FC) -o snlc_fit_exe


SIMSED_fudge_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread -lstdc++ $(ROOTLIBS)

SIMSED_extractSpec_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread -lstdc++ $(ROOTLIBS)

SIMSED_check_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread -lstdc++ $(ROOTLIBS)

SIMSED_rebin_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread -lstdc++ $(ROOTLIBS)

simlib_coadd_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread -lstdc++ $(ROOTLIBS)


unfold_snpar_exe_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm -lpthread -lstdc++ $(ROOTLIBS)


SALT2mu_LDADD =  $(CFITSIOLIB) $(GSLLIB) -lm $(GSLCBLASLIB) -lpthread -lz -lstdc++ $(ROOTLIBS) $(FLIBS)

# Need to force automake to use gfortran as the linker
SALT2mu_LINK = $(FC) -o SALT2mu


SNTOOLS_OUTPUT = \
        sntools_output.cpp

sntable_dump_SOURCES = sntable_dump.c sntools.c sntools_output.cpp

sntable_combine_SOURCES = sntable_combine.c sntools.c sntools_output.cpp

combine_fitres_exe_SOURCES = combine_fitres.c sntools.c sntools_output.cpp

wfit_exe_SOURCES = sntools.c sntools_output.cpp sntools_npz.cpp cnpy.cpp sntools_cosmology.c wfit.c

kcor_exe_SOURCES = kcor.c sntools.c MWgaldust.c sntools_spectrograph.c genmag_SEDtools.c sntools_dataformat_fits.c sntools_dataformat_text.c sntools_data.c sntools_output.cpp sntools_cosmology.c


sntools_snana_SRC = sntools.c sntools_gridmap.c sntools_nonlinearity.c sntools_fluxErrModels.c sntools_fluxErrModels_legacy.c sntools_simeff_obsolete.c sntools_dataformat_fits.c sntools_dataformat_text.c sntools_data.c sntools_modelgrid_read.c sntools_calib.c sntools_output.cpp yaml.c MWgaldust.c sntools_cosmology.c multiseason.c

genmag_fit_SRC = genmag_stretch2.c genmag_mlcs2k2.c genmag_SALT2.c genmag_BAYESN.c genmag_snoopy.c genmag_S11DM15.c genmag_extinction.c genmag_extrap.c

genmag_sim_SRC = genmag_NON1ASED.c genmag_NON1AGRID.c genmag_LCLIB.c genmag_PySEDMODEL.c genmag_SEDtools.c genmag_SIMSED.c $(genmag_fit_SRC)

snlc_sim_exe_SOURCES = snlc_sim.c $(sntools_snana_SRC) sntools_host.c sntools_wgtmap.c sntools_weaklens.c sntools_stronglens.c sntools_wronghost.c sntools_genSmear.c sntools_devel.c sntools_sim_readme.c sntools_sim_atmosphere.c sntools_trigger.c sntools_modelgrid_gen.c sntools_spectrograph.c sntools_genPDF.c sntools_genGauss_asym.c sntools_genExpHalfGauss.c inoue_igm.c $(genmag_sim_SRC)

#snlc_fit_exe_SOURCES = snlc_fit.f $(sntools_snana_SRC) sntools_genSmear.c sntools_devel.c sntools_zPDF_spline.c minuit.F $(genmag_fit_SRC)

SALT2mu_SOURCES = SALT2mu.c sntools.c sntools_output.cpp sntools_gridmap.c sntools_genGauss_asym.c minuit.F sntools_genExpHalfGauss.c

sntools_data_SRC = sntools_dataformat_fits.c sntools_dataformat_text.c sntools_data.c sntools_output.cpp sntools_cosmology.c

SIMSED_fudge_exe_SOURCES = SIMSED_fudge.c genmag_SIMSED.c sntools.c sntools_genSmear.c genmag_SEDtools.c genmag_SALT2.c genmag_extrap.c sntools_devel.c sntools_wgtmap.c sntools_gridmap.c MWgaldust.c $(sntools_data_SRC)

SIMSED_extractSpec_exe_SOURCES = SIMSED_extractSpec.c genmag_SEDtools.c sntools.c sntools_output.cpp sntools_cosmology.c MWgaldust.c sntools_wgtmap.c sntools_gridmap.c

SIMSED_check_exe_SOURCES = SIMSED_check.c sntools.c sntools_output.cpp sntools_cosmology.c sntools_wgtmap.c sntools_gridmap.c MWgaldust.c genmag_SEDtools.c genmag_SIMSED.c

SIMSED_rebin_exe_SOURCES = SIMSED_rebin.c sntools.c genmag_SEDtools.c sntools_output.cpp sntools_cosmology.c MWgaldust.c 

simlib_coadd_exe_SOURCES = simlib_coadd.c MWgaldust.c sntools.c sntools_output.cpp sntools_cosmology.c

unfold_snpar_exe_SOURCES = unfold_snpar.c sntools.c sntools_output.cpp 
