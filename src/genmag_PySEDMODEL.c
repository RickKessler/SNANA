/*****************************************
  Created Sep 2018
  Sep 10 2020: rename BYOSED -> more generic name PySEDMODEL

  C-wrapper to call python function that returns rest-frame SED,
  snd then C functions compute & return observer-frame magnitudes
  to the simulation. Also returns obs-frame spectra if requested.

  Each PySEDMODEL is associated with a separate gensed_[model].py :
     gensed_BYOSED.py : Build Your Own SED  (J.Pierel)
     gensed_SNEMO.py  : SNFactory model (Ben Rose)
     gensed_PYBAYESN.py : BayeSN model (Gautham Narayan, Stephen Thorp, Kaisey Mandel)
     gensed_AGN.py : AGN model (Qifeng Cheng, Konstantin Malanchev)

  Initial motivation is to build underlying "true" SED model to
  test SNIa model training. However, this function could in
  principle be used for SNCC or other transients.

  Invoke intrinsic scatter models as follows:
    GENMAG_SMEAR_MODELNAME: C11
    GENMAG_SMEAR_MODELNAME: G10:[MODELPATH_SAL2]

  Note that C11 model does not require extra args, but the
  G10 model needs to be linked to a SALT2 model.

  Mar 30 2019 RK - fix hc factors in spectra
  Apr 11 2019 RK - check for intrinsic scatter models (e.g., C11, G10 ...)
  Jul 12 2019 RK - store inputs in INPUTS_BYOSED struct, and add
                   internal DUMPFLAG_HOSTPAR in genmag_BYOSED().
  Sep 11 2020 RK - major refactor for BYOSED -> PySEDMODEL
  Feb 20 2021 RK - fix bug in genmag_PySEDMODEL that was preventing
                   BYOSED params from being returned.

  Mar 19 2021 RK - pass "RANSEED <ISEED>" via ARGLIST to
            init_genmag_PySEDMODEL.

  Apr 04 2021: if SED wavelength range does NOT cover a band, return
               mag = 666 (instead of 99) to instruct sim NOT to write this
               invalid flux to data file.

  Jun 13 2022 RK - if no python code, abort by default unless
                   GENMODEL_MSKOPT = 4096

  Aug 24 2022 D.Jones: fix genSpec to produce spectrum at requested phase,
                       and not only at peak.

  Sep 16 2022 RK
   + start new function prepEvent_PySEDMODEL() to receive list of
     all Tobs before calling fetch functions.
   + prep MODEL_NAME_AGN; python code to be added later.

  July 1 2023
    + new function set_lamRanges_PySEDMODEL sets wavelength ranges for each
      model, and aborts if MODEL_NAME is not found. This forces a human
      check for new PySEDMODELs.

 *****************************************/

#include  <stdio.h>
#include  <math.h>
#include  <stdlib.h>
#include  <sys/stat.h>

#include  "sntools.h"           // SNANA community tools
#include  "sntools_genSmear.h"
#include  "sntools_spectrograph.h"
#include  "sntools_devel.h"
#include  "MWgaldust.h"
#include  "genmag_SEDtools.h"
#include  "genmag_SIMSED.h"
#include  "genmag_PySEDMODEL.h"

#ifdef USE_PYTHON
#include  <Python.h>
//#include <numpy/arrayobject.h>
//#include <numpy/ndarrayobject.h>
PyObject *numpy_empty, *numpy_double; // numpy.empty, numpy.double

PyObject *geninit_PySEDMODEL ;

//int init_numpy(){
//  import_array(); // PyError if not successful
//  return 0;
//}

#endif

// ===============================================

#ifdef USE_PYTHON
void handle_python_exception(char *fnam, const char *desc) {
  // Exit with an error message for Python exceptions
  if (PyErr_Occurred() == NULL) {
    return ;
  }
  print_preAbort_banner(fnam);
  PyErr_Print();
  sprintf(c1err,"Python raised an uncaught exception");
  sprintf(c2err,"while %s", desc);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
}
#endif

// ===============================================

// =========================================================
void init_genmag_PySEDMODEL(char *MODEL_NAME, char *PATH_VERSION, int OPTMASK,
			    char *ARGLIST, char *NAMES_HOSTPAR  ) {

  // Read input directory file(s) for parameters characterizing
  // how to build your SED.
  //
  // Inputs:
  //  PATH_VERSION : points to model-param directory;
  //               :  passed from GENMODEL arg of sim-input
  //  OPTMASK      : bit mask of options; interpreted by python code.
  //               : OPTMASK=-1 is a flag to print options.
  //               :  passed from sim-input GENMODEL_MSKOPT: <MSKOPT>
  //
  //  ARGLIST      : string of options passed ONLY from the command-line,
  //                   snlc_sim.exe <myInput> GENMODEL_ARGLIST 'bla bla'
  //                 First ARGLIST element is hard-coded to be
  //                 "RANSEED <ISEED>" from sim-input file.
  //
  // NAMES_HOSTPAR : comma-separate list of names of host parameters.
  //                 Includes RV, AV, variables used in WGTMAP, and
  //                 HOSTLIB_STOREPAR list from sim-input file.
  //                 Duplicates are automatically removed.
  //
  // Mar 19 2021: always pass "RANSEED <ISEED>" via ARGLIST
  //

#ifdef USE_PYTHON
  PyObject *genmod, *genclass, *genmod_base, *gensed_base, *pargs;
#endif

  char *PyMODEL_NAME   = INPUTS_PySEDMODEL.MODEL_NAME ;
  char *PyCLASS_NAME   = INPUTS_PySEDMODEL.PyCLASS_NAME ;
  int  L, ipar, NPAR ;
  int  MEMD   = sizeof(double);
  int  MEMC   = sizeof(char);
  char comma[] = ",";
  char fnam[] = "init_genmag_PySEDMODEL" ;

  // -------------- BEGIN ------------------

  sprintf(BANNER, "%s", fnam);
  print_banner(BANNER);

  sprintf(PyMODEL_NAME, "%s",        MODEL_NAME);
  sprintf(PyCLASS_NAME, "gensed_%s", PyMODEL_NAME) ;

  printf("   %s PATH    = '%s' \n",  PyMODEL_NAME, PATH_VERSION);
  printf("   %s OPTMASK = %d \n",    PyMODEL_NAME, OPTMASK );
  printf("   %s ARGLIST = '%s' \n",  PyMODEL_NAME, ARGLIST );
  printf("   %s HOSTPAR = '%s' \n",  PyMODEL_NAME, NAMES_HOSTPAR );
  fflush(stdout);

  // - - - - - - - - - - -
  // store inputs in global (RK - Jul 12 2019)
  L=strlen(PATH_VERSION)+4;
  INPUTS_PySEDMODEL.PATH    = (char*)malloc(L*MEMC);

  L=strlen(ARGLIST)+4;
  INPUTS_PySEDMODEL.ARGLIST = (char*)malloc(L*MEMC);

  L=strlen(NAMES_HOSTPAR)+4;
  INPUTS_PySEDMODEL.NAMES_HOSTPAR=(char*)malloc(L*MEMC);

  INPUTS_PySEDMODEL.OPTMASK = OPTMASK;

  sprintf(INPUTS_PySEDMODEL.PATH,          "%s", PATH_VERSION);
  sprintf(INPUTS_PySEDMODEL.ARGLIST,       "%s", ARGLIST );
  sprintf(INPUTS_PySEDMODEL.NAMES_HOSTPAR, "%s", NAMES_HOSTPAR );

  // split comma-separated HOSTPAR_NAMES and store each
  // name separately (for debug dumps)
  for(ipar=0; ipar < MXHOSTPAR_PySEDMODEL; ipar++ )
    { INPUTS_PySEDMODEL.NAME_ARRAY_HOSTPAR[ipar] = (char*)malloc(60*MEMC);  }

  splitString(NAMES_HOSTPAR, comma, fnam, MXHOSTPAR_PySEDMODEL,
	      &NPAR, INPUTS_PySEDMODEL.NAME_ARRAY_HOSTPAR );

  // - - - - - - - - - - -
  // print summary of filter info
  filtdump_SEDMODEL();

  // set wavelength range
  set_lamRanges_PySEDMODEL(MODEL_NAME);

  /* xxx mark delete July 1 2023 xxxxxxx
  // Dec 2022: set SED and filtercen limits; here they are hard wired, 
  //  but we should really call a python method to return the values  
  // start with generi default ...
  SEDMODEL.RESTLAMMIN_FILTERCEN =  2000.0 ;
  SEDMODEL.RESTLAMMAX_FILTERCEN = 20000.0 ;
  SEDMODEL.LAMMIN_ALL           =  1000.0 ;
  SEDMODEL.LAMMAX_ALL           = 25000.0 ;
  xxxxxxx end mark xxxxxxxxx */


  // init a few C struct items
  Event_PySEDMODEL.LAST_EXTERNAL_ID = -9;
  Event_PySEDMODEL.LAM  = (double*) malloc( MXLAM_PySEDMODEL*MEMD ) ;
  Event_PySEDMODEL.SED  = (double*) malloc( MXLAM_PySEDMODEL*MEMD ) ;

  SEDMODEL_MWEBV_LAST     = -999.   ;
  SEDMODEL_HOSTXT_LAST.AV = -999.   ;
  SEDMODEL_HOSTXT_LAST.z  = -999.   ;

#ifndef USE_PYTHON

  bool ALLOW_C_ONLY = ( OPTMASK & OPTMASK_ALLOW_C_ONLY ) > 0;

  if ( ALLOW_C_ONLY ) {
    printf("\n no python ==> read SALT2 template with C code. \n");
    // read M0 surface of SALT2 model, corresponding to x1=c=0.
    // This is for debug-comparisons with SALT2 model where x1=c=0.
    
    read_SALT2_template0();
  }
  else {
    sprintf(c1err,"SNANA is not linked to Python-model code.");
    sprintf(c2err,"Pass GENMODEL_MSKOPT += %d to debug with C code and SALT2", 
	    OPTMASK_ALLOW_C_ONLY);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

#endif


#ifdef USE_PYTHON
  printf("\n\t Begin %s python-init from C code ... \n", PyMODEL_NAME );   fflush(stdout);

  Py_Initialize();
  int nResult1 = PyRun_SimpleStringFlags("import numpy", NULL);
  int nResult2 = PyRun_SimpleStringFlags("import os", NULL);
  int nResult3 = PyRun_SimpleStringFlags("import optparse",NULL);
  int nResult4 = PyRun_SimpleStringFlags("import configparser",NULL);
  if ((nResult1 != 0) || (nResult2 != 0) || (nResult3 != 0) || (nResult4 != 0)) {
    sprintf(c1err,"Could not import numpy, os, optparse or configparser");
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // Get numpy empty, we will reuse it
  PyObject *numpy = PyImport_ImportModule("numpy");
  numpy_empty = PyObject_GetAttrString(numpy, "empty");
  numpy_double = PyObject_GetAttrString(numpy, "double");
  handle_python_exception(fnam, "importing numpy and getting numpy.empty & numpy.double");
  Py_DECREF(numpy);

  genmod_base = PyImport_ImportModule("gensed_base");
  if (genmod_base == NULL) {
    handle_python_exception(fnam, "tried to import gensed_base");
  }
  gensed_base = PyObject_GetAttrString(genmod_base, "gensed_base");
  if (gensed_base == NULL) {
    sprintf(c1err,"Could not find gensed_base class in gensed_base module");
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  genmod = PyImport_ImportModule(PyCLASS_NAME);
  handle_python_exception(fnam, "tried to import pySED model");

  if (genmod == NULL) {
    sprintf(c1err,"Could not import module %s", PyCLASS_NAME);
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // xxxx  genclass = PyObject_GetAttrString(genmod, "genmag_BYOSED");
  genclass = PyObject_GetAttrString(genmod, PyCLASS_NAME);
  if (genclass == NULL) {
    sprintf(c1err,"Could not find %s class",PyCLASS_NAME);
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  if (PyObject_IsSubclass(genclass, gensed_base) != 1) {
    handle_python_exception(fnam, "checking model class is a subclass of gensed_base.gensed_base");
    sprintf(c1err,"%s is not a subclass of gensed_base.gensed_base",PyCLASS_NAME);
    sprintf(c2err,"Please check SNANA/src/gensed_base.py and inherit gensed_base from your model");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  pargs  = Py_BuildValue("(siss)",PATH_VERSION,OPTMASK,ARGLIST,NAMES_HOSTPAR);
  geninit_PySEDMODEL = PyObject_CallObject(genclass, pargs);
  if (geninit_PySEDMODEL == NULL) {
    handle_python_exception(fnam, "tried to construct pySED model class");
  }

  Py_DECREF(genmod);
  Py_DECREF(genclass);
  Py_DECREF(pargs);

  printf("\t Finished %s python-init from C code \n", PyMODEL_NAME );
  fflush(stdout);
#endif

  // -----------------------------------------------------
  // set SED par names and allocate arrays for parameters
  Event_PySEDMODEL.PARVAL  = (double*) malloc(MXPAR_PySEDMODEL*MEMD );
  Event_PySEDMODEL.PARNAME = (char**)  malloc(MXPAR_PySEDMODEL*sizeof(char*) );
  Event_PySEDMODEL.NPAR    = 0 ;
  for(ipar=0; ipar < MXPAR_PySEDMODEL; ipar++ )
    { Event_PySEDMODEL.PARNAME[ipar] = (char*) malloc(40*MEMC);  }
#ifdef USE_PYTHON
  NPAR = fetchParNames_PySEDMODEL(Event_PySEDMODEL.PARNAME);
#else
  NPAR = 2;  // test with C code only
  sprintf(Event_PySEDMODEL.PARNAME[0],"TEST0");
  sprintf(Event_PySEDMODEL.PARNAME[1],"TEST1");
#endif

  Event_PySEDMODEL.NPAR = NPAR;
  printf("\t %s parameters to store in data files:\n", PyMODEL_NAME);
  for(ipar=0; ipar < NPAR; ipar++ )
    { printf("\t\t %s \n", Event_PySEDMODEL.PARNAME[ipar] ); }

  // - - - -
  printf("\n\t Done with %s \n", fnam);
  fflush(stdout);

  return ;

} // end init_genmag_PySEDMODEL



/* xxxxxxxxxxxxx mark delete July 1 2023 xxxxxxxx
// =========================================================
void get_MODEL_NAME_PySEDMODEL(char *PATH,char *MODEL_NAME) {

  // For input PATH, return MODEL_NAME

  int i;
  char *ptrModel;
  char fnam[] = "get_MODEL_NAME_PySEDMODEL" ;
  // ----------- BEGIN -----------

  // load all possible model choices
  load_PySEDMODEL_CHOICE_LIST();

  MODEL_NAME[0] = 0;
  for ( i=0; i < NCHOICE_PySEDMODEL; i++  ) {
    ptrModel = PySEDMODEL_CHOICE_LIST[i] ;
    if ( strstr(PATH,ptrModel) ) { sprintf(MODEL_NAME, "%s", ptrModel); }
  }

  if ( strlen(MODEL_NAME) == 0 ) {
    sprintf(c1err,"Could not determine MODEL_NAME from") ;
    sprintf(c2err,"PATH='%s' ", PATH );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  return;
} // end get_MODEL_NAME
xxxxxxxxx end mark xxxxxxxx */



// ===================================================
void set_lamRanges_PySEDMODEL(char *MODEL_NAME) {

  // Created July 1 2023
  // set wavelengh range for center of filter, and for SED
  
  char fnam[] = "set_lamRanges_PySEDMODEL";

  // ----------- BEGIN ---------

  if ( strcmp(MODEL_NAME,MODEL_NAME_AGN) == 0 ) {
    SEDMODEL.LAMMIN_ALL           =  100.0 ;
    SEDMODEL.LAMMAX_ALL           = 25000.0 ;
    SEDMODEL.RESTLAMMIN_FILTERCEN =  500.0 ;
    SEDMODEL.RESTLAMMAX_FILTERCEN = 20000.0 ;
  }

  else if ( strcmp(MODEL_NAME,MODEL_NAME_PYBAYESN) == 0 ) { 
    SEDMODEL.LAMMIN_ALL           =  2000.0 ;
    SEDMODEL.LAMMAX_ALL           = 20000.0 ;
    SEDMODEL.RESTLAMMIN_FILTERCEN =  4000.0 ;
    SEDMODEL.RESTLAMMAX_FILTERCEN = 16000.0 ;
  }

  else if ( strcmp(MODEL_NAME,MODEL_NAME_BYOSED) == 0 ) { 
    SEDMODEL.LAMMIN_ALL           =  1000.0 ;
    SEDMODEL.LAMMAX_ALL           = 25000.0 ;
    SEDMODEL.RESTLAMMIN_FILTERCEN =  2000.0 ;
    SEDMODEL.RESTLAMMAX_FILTERCEN = 20000.0 ;
  }

  else if ( strcmp(MODEL_NAME,MODEL_NAME_SNEMO) == 0 ) { 
    SEDMODEL.LAMMIN_ALL           =  1000.0 ;
    SEDMODEL.LAMMAX_ALL           = 25000.0 ;
    SEDMODEL.RESTLAMMIN_FILTERCEN =  2000.0 ;
    SEDMODEL.RESTLAMMAX_FILTERCEN = 20000.0 ;
  }

  else {
    sprintf(c1err,"Cannot set wavelength ranges for MODEL_NAME=%s", MODEL_NAME);    
    sprintf(c2err,"Check MODEL_NAME");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  printf("\n\t Hard-wired wavelength ranges for %s: \n", MODEL_NAME);
  printf("\t   FILTERCEN: %.1f to %.1f A\n", 
	 SEDMODEL.RESTLAMMIN_FILTERCEN, SEDMODEL.RESTLAMMAX_FILTERCEN );
  printf("\t   SED: %.1f to %.1f A\n", 
	 SEDMODEL.LAMMIN_ALL, SEDMODEL.LAMMAX_ALL);
  fflush(stdout);

  return;

} // end set_lamRanges_PySEDMODEL

// =========================================================
void prepEvent_PySEDMODEL(int EXTERNAL_ID, double zHEL,
			  int NHOSTPAR, double *HOSTPAR_LIST,
			  int NOBS,  double *TOBS_LIST,
			  int NSPEC, double *TSPEC_LIST ) {

  // Created Sep 2022 by R.Kessler
  // Interface to prepare seds for all epochs (over all bands),
  // instead of genmag_PySEDMODEL that is called for each band.
  // Note that genmag_PySEDMODEL is still called the same way,
  // so this function is an option for preparing all SEDs for
  // the event rather than by band.
  //
  // Initial use is for AGN where model must always go forward in time,
  // and thus cannot processed by band.
  //
  // Inputs:
  //   EXTERNAL_ID  :  SNID passed from main program
  //   zHEL  : heliocentric redshift (to allow for z-dependent model)
  //  NHOSTPAR:      number of host params
  //  HOSTPAR_LIST:  host property values
  //  NOBS:          total number of observations, all bands
  //  TOBS_LIST:     observer-frame times w.r.t. PEAKMJD
  //  NSPEC:         number of spectra
  //  TSPEC_LIST     TOBS for each spectrum
  //
  // Mar 27 2024: pass extra list of NSPEC and TSPEC_LIST so that
  //   this works for both broadband and spectrograph observations.

  #ifdef USE_PYTHON
 
  int     NOBS_ALL = NOBS + NSPEC;
  int     MEMI     = (NOBS_ALL+10) * sizeof(int);
  int     MEMD     = (NOBS_ALL+10) * sizeof(double);
  int    *INDEX_SORT = (int   *)   malloc(MEMI);
  double *TOBS_ALL   = (double*)   malloc(MEMD);

  char *MODEL_NAME   = INPUTS_PySEDMODEL.MODEL_NAME ;

  double *arrTrest;
  double z1      = 1.0 + zHEL;
  double z1inv   = 1.0/z1;
  double Tobs_template;
  int o, o_sort, i, ihost, NOBS_STORE;
  char fnam[] = "prepEvent_PySEDMODEL";
 

  PyObject *pHOSTPARS, *pTrest, *prepmeth;
  Py_buffer bufTrest = {NULL, NULL};

  // ------------- BEGIN ------------

  for(o = 0; o < NOBS;  o++ ) { TOBS_ALL[o]      = TOBS_LIST[o]; }
  for(o = 0; o < NSPEC; o++ ) { TOBS_ALL[NOBS+o] = TSPEC_LIST[o]; }

  // sort TOBS_LIST in increasing order
  int ORDER_SORT = +1; // ==> increasing
  sortDouble(NOBS_ALL, TOBS_ALL, ORDER_SORT, INDEX_SORT);

  // for recurring models (e.g., AGN), we need one extra DIA tempalte
  bool RECUR = ( strcmp(MODEL_NAME,MODEL_NAME_AGN) == 0 );
  
  if ( RECUR ) {
    o_sort = INDEX_SORT[0];
    Tobs_template = TOBS_ALL[o_sort] - 365.0 ; // arbitrary value of one year
    NOBS_STORE = NOBS_ALL + 1;
  } else {
    Tobs_template = -1.0E8; // not used
    NOBS_STORE = NOBS_ALL;
  }
  Event_PySEDMODEL.Tobs_template = Tobs_template ;

  // Creating an empty numpy array to handle Trest
  // It would be more clear to use numpy C-API, but we wouldn't introduce numpy
  // as a build dependency.
  pTrest = PyObject_CallFunction(numpy_empty, "(iO)", NOBS_STORE, numpy_double);
  handle_python_exception(fnam, "creating ndarray for trest");
  if (PyObject_GetBuffer(pTrest, &bufTrest, PyBUF_CONTIG) != 0) {
    sprintf(c1err,"pTrest must be a contiguous numpy array");
    sprintf(c2err,"??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  arrTrest = (double*) bufTrest.buf;

  if ( RECUR ) {
    arrTrest[0] = Tobs_template * z1inv;
  }
  for(o=0; o < NOBS_ALL; o++ )  { 
    o_sort = INDEX_SORT[o];
    i = o + (int) RECUR;
    arrTrest[i] = TOBS_ALL[o_sort] * z1inv ;
  }

  // Give Python model Trest array so it is prepared
  // for fetching Trest in an arbitrary order
  prepmeth = PyObject_GetAttrString(geninit_PySEDMODEL, "prepEvent");
  handle_python_exception(fnam, "getting prepEvent method");
  
  pHOSTPARS = PyTuple_New(NHOSTPAR);
  for(ihost=0; ihost < NHOSTPAR; ihost++ ){
    PyTuple_SetItem(pHOSTPARS,ihost,PyFloat_FromDouble(HOSTPAR_LIST[ihost]));
  }

  PyObject_CallFunction(prepmeth, "(OiO)", pTrest, EXTERNAL_ID, pHOSTPARS);
  handle_python_exception(fnam, "calling prepEvent");

  PyBuffer_Release(&bufTrest);
  Py_DECREF(pHOSTPARS);
  Py_DECREF(pTrest);
  Py_DECREF(prepmeth);

  free(INDEX_SORT);

  #else
  // We don't need to do a template for SALT2
  Event_PySEDMODEL.Tobs_template = -1.0E8;
  #endif

  return;

} // end prepEvent_PySEDMODEL

// =========================================================
void genmag_PySEDMODEL(int EXTERNAL_ID, double zHEL, double zCMB, double MU,
		       double MJDOFF, double MWEBV, int NHOSTPAR, double *HOSTPAR_LIST,
		       int IFILT_OBS, int NOBS, double *TOBS_list,
		       double *MAGOBS_list, double *MAG_TEMPLATE,
		       double *MAGERR_list ) {

  // Created Sep 2018
  //
  // Inputs:
  //   EXTERNAL_ID  : set NEWEVT_FLAG logical when this changes
  //   zHEL         : helio redshift
  //   MU           : distance modulus
  //   MJDOFF       : MJD = MJDOFF + TOBS
  //   MWEBV        : E(B-V) for Milky Way
  //   NHOSTPAR     : Number of host params in HOSTPAR_LIST
  //   HOSTPAR_LIST : RV, AV, LOGMASS, SFR ...
  //   IFILT_OBS    : absolute filter index
  //   NOBS         : number of observations
  //   TOBS_list    : list of MJD-PEAKMJD
  //
  // Outputs"
  //   MAGOBS_list   : list of true mags
  //   MAG_TEMPLATE  : mag of template (for DIA)
  //   MAGERR_list   : list of mag errors (place-holder, in case)
  //
  // Feb 20 2021 RK
  //   + fix bug that was preventing call to fetchParVal_PySEDMODEL()
  // 
  // Sep 16 2022 RK - compute new output *MAG_TEMPLATE arg to handle AGN
  //
  // May 5 2023 RK - add MJDOFF arg (for defining AGN t_transition as MJD)
  //

  int   MXLAM      = MXLAM_PySEDMODEL;
  char *MODEL_NAME = INPUTS_PySEDMODEL.MODEL_NAME ;
  char *PyCLASS_NAME = INPUTS_PySEDMODEL.PyCLASS_NAME ;

  double RV_host = HOSTPAR_LIST[0];
  double AV_host = HOSTPAR_LIST[1];
  double z1    = 1.0 + zHEL ;
  double z1inv   = 1.0/z1;
  double *LAM  = Event_PySEDMODEL.LAM;
  double *SED  = Event_PySEDMODEL.SED ;
  bool DO_TEMPLATE = Event_PySEDMODEL.Tobs_template > -1.0E7 ;

  int  ifilt  = IFILTMAP_SEDMODEL[IFILT_OBS] ; // sparse filter index
  double  ZP  = FILTER_SEDMODEL[ifilt].ZP ;    // ZP for flux->mag
  double x0   = pow(10.0,-0.4*MU);             // dimming from dist. mod.
  int    NEWEVT_FLAG = 0 ;
  int    NEWEVT_FLAG_TMP ;
  int    DUMPFLAG_HOSTPAR = 0 ;
  int    FLAG_Finteg;

  int    NLAM, o, ipar ;
  double Tobs, Trest, FLUXSUM_OBS, FspecDUM[2], magobs ;
  char fnam[] = "genmag_PySEDMODEL" ;

   #ifdef USE_PYTHON
  // python declarations here

   #endif

  // ------------ BEGIN -----------

  // init output args
  *MAG_TEMPLATE = 99.0; 
  for(o=0; o < NOBS; o++ )  { MAGOBS_list[o] = 99.0 ; }

  // check of this is a new event, or same event
  // with different epoch
  Event_PySEDMODEL.EXTERNAL_ID = EXTERNAL_ID ;
  if ( EXTERNAL_ID != Event_PySEDMODEL.LAST_EXTERNAL_ID )
    { NEWEVT_FLAG=1; }


  // RK - check internal flag to dump host params
  if ( DUMPFLAG_HOSTPAR && NEWEVT_FLAG ) {
    printf(" xxx ----------------------------------------------\n");
    printf(" xxx %s: dump HOSTPAR_LIST for EXTERNAL_ID=%d \n",
	   fnam, EXTERNAL_ID);
    for(ipar=0; ipar < NHOSTPAR; ipar++ ) {
      printf(" xxx (%2d) %-20s = %f \n", ipar,
	     INPUTS_PySEDMODEL.NAME_ARRAY_HOSTPAR[ipar],
	     HOSTPAR_LIST[ipar] );
    }
    fflush(stdout);
  }

  /*
  printf(" xxx ------------------------------------ \n" ) ;
  printf(" xxx %s: process z=%.3f MU=%.3f RV=%3.1f IFILT_OBS=%d(%s) \n",
	 fnam, zHEL, MU, MWXT_SEDMODEL.RV, IFILT_OBS, cfilt );
  printf(" xxx RVMW=%3.1f  MWEBV=%.3f      RV_host=%3.1f AV=%4.2f \n",
	 MWXT_SEDMODEL.RV, MWEBV,   RV_host, AV_host );
  */


  // store table info for Galactic & host extinction
  fill_TABLE_MWXT_SEDMODEL(MWXT_SEDMODEL.RV, MWEBV);
  fill_TABLE_HOSTXT_SEDMODEL(RV_host, AV_host, zHEL);   // July 2016

  int NOBS_LOCAL = NOBS;
  if ( DO_TEMPLATE ) { NOBS_LOCAL++ ; }

  for(o=0; o < NOBS_LOCAL; o++ ) {

    if ( o < NOBS ) 
      { Tobs = TOBS_list[o]; }
    else
      { Tobs =  Event_PySEDMODEL.Tobs_template; }

    Trest = Tobs * z1inv;
    if (o == 0 )
      { NEWEVT_FLAG_TMP = NEWEVT_FLAG; }
    else
      { NEWEVT_FLAG_TMP = 0; }

    fetchSED_PySEDMODEL(EXTERNAL_ID, NEWEVT_FLAG_TMP, Trest,
			MXLAM, HOSTPAR_LIST, &NLAM, LAM, SED);
    Event_PySEDMODEL.NLAM = NLAM ;

    // integrate redshifted SED to get observer-frame flux in IFILT_OBS band.
    // FLUXSUM_OBS is returned (ignore FspecDUM)
    INTEG_zSED_PySEDMODEL(0, IFILT_OBS, Tobs, zHEL, x0,RV_host,AV_host,
			  NLAM, LAM, SED,
			  &FLUXSUM_OBS, FspecDUM, &FLAG_Finteg ); //<=returned

    // convert calibrated flux into true magnitude
    if ( FLAG_Finteg == 0 )
      { magobs = ZP - 2.5*log10(FLUXSUM_OBS); }
    else if ( FLAG_Finteg == (int)MAG_ZEROFLUX )
      { magobs = MAG_ZEROFLUX ; }
    else if ( FLAG_Finteg == (int)MAG_UNDEFINED )
      { magobs = MAG_UNDEFINED ; }

    if ( o < NOBS ) {
      MAGOBS_list[o] = magobs;  // load output array
      MAGERR_list[o] = 0.01;    // not used
    }
    else {
      *MAG_TEMPLATE = magobs ;
    }

  }

  // for NEW EVENT, store SED parameters so that sim can
  // write them to data files
  // hack

  if ( NEWEVT_FLAG ) {
    fetchParVal_PySEDMODEL(Event_PySEDMODEL.PARVAL);
  }

  // keep track of last ID
  Event_PySEDMODEL.LAST_EXTERNAL_ID = EXTERNAL_ID ;

  return ;


} // end genmag_PySEDMODEL


// ============================================================
//  FETCH UTILITIES TO RETURN EXTRA INFO TO MAIN PROGRAM
// =============================================================


// ================================================
int fetchParNames_PySEDMODEL(char **parNameList) {

  // Pass name of each parameter to calling function, so that
  // these parameters can be included in the data files.
  // Function Returns number of parameters used to create
  // each SED.  **parNameList is a list of parameter names.
  //
  // Called once during init stage.

  char *MODEL_NAME = INPUTS_PySEDMODEL.MODEL_NAME ;
  int NPAR = 0 ;
  char fnam[] = "fetchParNames_PySEDMODEL" ;

#ifdef USE_PYTHON
  // python declarations here

  int ipar;
  PyObject *parnamesmeth,*pNames,*pnamesitem;
  printf("fetching parameter names from Python\n");
  // David: need your python magic to return these string names.

  parnamesmeth  = PyObject_GetAttrString(geninit_PySEDMODEL, "fetchParNames");
  handle_python_exception(fnam, "getting fetchParNames method");
  //xx  parnamesmeth  = PyObject_GetAttrString(geninit_PySEDMODEL,
  //xx					 "fetchParNames_BYOSED");

  pNames  = PyObject_CallObject(parnamesmeth, NULL);
  handle_python_exception(fnam, "fetching paramweter names");
  if (PySequence_Check(pNames) != 1) {
    sprintf(c1err,"fetchParNames is expected to return a sequence of str, for example a list");
    sprintf(c2err,"but it is %s instead", PyUnicode_AsUTF8(PyObject_Type(pNames)));
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  NPAR = PySequence_Length(pNames);
  if (NPAR == -1) {
    handle_python_exception(fnam, "getting parameter sequence length");
  }
  if (NPAR > MXPAR_PySEDMODEL) {
    sprintf(c1err,"fetchParNames returned too long sequence of parameters");
    sprintf(c2err,"%d is fetched, %d is the maximum allowed value", NPAR, MXPAR_PySEDMODEL);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  for(ipar=0; ipar < NPAR; ipar++ ) {
    pnamesitem = PySequence_GetItem(pNames,ipar);
    if (pnamesitem == NULL) {
      handle_python_exception(fnam, "getting parameter from fetched object");
    }
    sprintf(parNameList[ipar],PyUnicode_AsUTF8(pnamesitem));
    Py_DECREF(pnamesitem);
  }

  Py_DECREF(pNames);
  Py_DECREF(parnamesmeth);

#endif


  return(NPAR) ;

} // fetchParNames_PySEDMODEL


void fetchParVal_PySEDMODEL(double *parVal) {

  // return list of parameters to calling function (sim)
  // so that these parameters can be included in the
  // data files.
  //
  // Called once per event.

#ifdef USE_PYTHON
  PyObject *parvalmeth,*pParVal,*pargs;
#endif
  char *MODEL_NAME = INPUTS_PySEDMODEL.MODEL_NAME ;
  double val;
  char **parNameList ;
  int NPAR, ipar;
  char fnam[] = "fetchParVal_PySEDMODEL" ;

  // ------------- BEGIN ------------------

  NPAR = Event_PySEDMODEL.NPAR;
  parNameList = Event_PySEDMODEL.PARNAME;
  // David: need python function to return these values.
#ifdef USE_PYTHON

  parvalmeth  = PyObject_GetAttrString(geninit_PySEDMODEL, "fetchParVals");
  handle_python_exception(fnam, "getting fetchParVals attribute of the class");

  for(ipar=0; ipar < NPAR; ipar++ ) {
    pParVal  = PyObject_CallFunction(parvalmeth, "(s)", parNameList[ipar]);
    handle_python_exception(fnam, "fetching a parameter value");
    val          = PyFloat_AsDouble(pParVal);
    parVal[ipar] = val;
    // printf("   PARVAL    = '%d' \n",  val);
  }

  Py_DECREF(pParVal);
  Py_DECREF(parvalmeth);
#endif

#ifndef USE_PYTHON
  for(ipar=0; ipar < NPAR; ipar++ ) {
    val = (double)Event_PySEDMODEL.EXTERNAL_ID + 0.1*(double)ipar ;
    parVal[ipar] = val;
  }
#endif

  return ;

} // end fetchParVal_PySEDMODEL

// =================================================
void fetchSED_PySEDMODEL(int EXTERNAL_ID, int NEWEVT_FLAG, double Trest, int MXLAM,
			 double *HOSTPAR_LIST, int *NLAM_SED,
			 double *LAM_SED, double *FLUX_SED) {

  // return rest-frame SED to calling function;
  // Inputs:
  //   EXTERNAL_ID  :  SNID passed from main program
  //   NEWEVT_FLAG  :  logical flag: True for new event
  //   Trest        : rest frame epochs (Trest=0 at peak)
  //   MXLAM        : abort iof *NLAM > MXLAM
  //   HOSTPAR_LIST : RV, AV, LOGMAS ...
  //
  // Output
  //  *NLAM_SED  : number of wavelenth bins for SED
  //  *LAM_SED   : array of wavelengths for SED
  //  *FLUX_SED  : SED flux in each wave bin (erg/s/cm^2/A)
  //               Note that this is flux, not dF/dLam
  //

  char *MODEL_NAME = INPUTS_PySEDMODEL.MODEL_NAME ;
  char fnam[] = "fetchSED_PySEDMODEL" ;

  // ------------ BEGIN -----------

  *NLAM_SED = 0 ; // init output

#ifdef USE_PYTHON
  PyObject *pmeth, *pargs, *pargs2, *pLAM, *pFLUX, *plammeth;
  Py_buffer bufLAM = {NULL, NULL};
  Py_buffer bufFLUX = {NULL, NULL};
  int NLAM, ilam, ihost;
  //int numpy_initialized =  init_numpy();

  // python declarations here
  pmeth  = PyObject_GetAttrString(geninit_PySEDMODEL, "_fetchSED");
  plammeth  = PyObject_GetAttrString(geninit_PySEDMODEL, "_fetchSED_LAM");

  pargs = PyTuple_New(5);
  pargs2 = PyTuple_New(sizeof(HOSTPAR_LIST));
  PyTuple_SetItem(pargs,0,PyFloat_FromDouble(Trest));
  PyTuple_SetItem(pargs,1,PyLong_FromLong(MXLAM));
  PyTuple_SetItem(pargs,2,PyLong_FromLong(EXTERNAL_ID));
  PyTuple_SetItem(pargs,3,PyLong_FromLong(NEWEVT_FLAG));

  for(ihost=0; ihost < sizeof(HOSTPAR_LIST); ihost++ ){
    PyTuple_SetItem(pargs2,ihost,PyFloat_FromDouble(HOSTPAR_LIST[ihost]));
  }
  PyTuple_SetItem(pargs,4,pargs2);
  
  pLAM   = PyEval_CallObject(plammeth, NULL);
  handle_python_exception(fnam, "calling _fetchSED_LAM method");
  pFLUX  = PyEval_CallObject(pmeth, pargs);
  handle_python_exception(fnam, "calling _fetchSED method");

  Py_DECREF(pmeth);
  Py_DECREF(plammeth);


  if (PyObject_CheckBuffer(pLAM) != 1) {
    sprintf(c1err,"_fetchSED_LAM must return numpy array");
    sprintf(c2err,"type of return value is %s",PyUnicode_AsUTF8(PyObject_Str(PyObject_Type(pLAM))));
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  if (PyObject_CheckBuffer(pFLUX) != 1) {
    sprintf(c1err,"_fetchSED must return numpy array");
    sprintf(c2err,"type of return value is %s",PyUnicode_AsUTF8(PyObject_Str(PyObject_Type(pFLUX))));
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  };

  if (PyObject_GetBuffer(pLAM, &bufLAM, PyBUF_FULL_RO) != 0) {
    handle_python_exception(fnam, "setting buffer from pLAM");
  }
  if (PyObject_GetBuffer(pFLUX, &bufFLUX, PyBUF_FULL_RO) != 0) {
    handle_python_exception(fnam, "setting buffer from pFLUX");
  }

  if (bufLAM.itemsize != sizeof(double)) {
    sprintf(c1err,"_fetchSED_LAM must return numpy array with np.float64 dtype");
    sprintf(c2err,"itemsize of returned dtype is %d",bufLAM.itemsize);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  if (bufFLUX.itemsize != sizeof(double)) {
    sprintf(c1err,"_fetchSED must return numpy array with np.float64 dtype");
    sprintf(c2err,"itemsize of returned dtype is %d",bufFLUX.itemsize);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  NLAM = bufLAM.len / bufLAM.itemsize;
  if (NLAM != bufFLUX.len / bufFLUX.itemsize) {
    sprintf(c1err,"size of array returned by _fetchSED_LAM doesn't equal to one returned by _fetchSED");
    sprintf(c2err,"NLAM = %d, NFLUX = %d",bufLAM.len / bufLAM.itemsize, bufFLUX.len / bufFLUX.itemsize);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  PyBuffer_ToContiguous(LAM_SED, &bufLAM, bufLAM.len, 'C');
  PyBuffer_ToContiguous(FLUX_SED, &bufFLUX, bufFLUX.len, 'C');

  *NLAM_SED = NLAM;

  PyBuffer_Release(&bufLAM);
  PyBuffer_Release(&bufFLUX);
  Py_DECREF(pLAM);
  Py_DECREF(pFLUX);
  Py_DECREF(pargs);
  Py_DECREF(pargs2);

#endif


#ifndef USE_PYTHON
  // C-code return SALT2 SED interpolated to Trest
  int NLAM, NDAY, ilam, iday, jf0, jf1, jf2 ;
  double TMPDAY[3], TMPSED[3],  F_interp, FSCALE ;

  NLAM = TEMP_SEDMODEL.NLAM ;    *NLAM_SED = NLAM ;
  NDAY = TEMP_SEDMODEL.NDAY ;
  iday = quickBinSearch(Trest, NDAY,TEMP_SEDMODEL.DAY, fnam, fnam);
  if ( iday >= NDAY-2 ) { iday = NDAY-3; }
  TMPDAY[0] = TEMP_SEDMODEL.DAY[iday+0] ;
  TMPDAY[1] = TEMP_SEDMODEL.DAY[iday+1] ;
  TMPDAY[2] = TEMP_SEDMODEL.DAY[iday+2] ;

  // hard-code 0.27 mag offset here to avoid reading SALT2.INFO file
  FSCALE = pow(10.0,-0.4*0.27);

  for(ilam=0; ilam < NLAM; ilam++ ) {
    // interpolate flux to Trest
    jf0  = NLAM*(iday+0) + ilam ;  TMPSED[0] = TEMP_SEDMODEL.FLUX[jf0];
    jf1  = NLAM*(iday+1) + ilam ;  TMPSED[1] = TEMP_SEDMODEL.FLUX[jf1];
    jf2  = NLAM*(iday+2) + ilam ;  TMPSED[2] = TEMP_SEDMODEL.FLUX[jf2];
    F_interp = quadInterp( Trest, TMPDAY, TMPSED, fnam);
    LAM_SED[ilam]  = TEMP_SEDMODEL.LAM[ilam] ;
    FLUX_SED[ilam] = F_interp * FSCALE;
  }
#endif


  // abort if too many wavelength bins
  if (NLAM >= MXLAM ) {
    sprintf(c1err,"NLAM=%d exceeds bound of %d", NLAM, MXLAM);
    sprintf(c2err,"Trest=%.2f ", Trest );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  return;

} // end fetchSED_PySEDMODEL


// =====================================================
void INTEG_zSED_PySEDMODEL(int OPT_SPEC, int ifilt_obs, double Tobs,
			   double zHEL, double x0,
			   double RV_host, double AV_host,
			   int NLAM, double *LAM, double *SED,
			   double *Finteg, double *Fspec, int *FLAG_Finteg) {

  // Created Dec 2018 by R.K.
  // Return integrated obs-frame flux in filter passband
  // Be careful with rest-frame and obs-frame
  //
  // Inputs:
  //   OPT_SPEC     1 -> return obs-frame spectrum withing filter-band.
  //   ifilt_obs    filter index
  //   zHEL         helio redshift
  //   Tobs         MJD - MJDpeak (for intrinsic scatter model)
  //   x0           flux scale from distance modulus
  //   RV_host      RV for host (optional)
  //   AV_host      AV for host
  //   NLAM         Number of SED bins
  //  *LAM          array of rest-frame wavelengths to define SED
  //  *SED          array of rest-fame SED model fluxes
  //
  // Outputs
  //   *Finteg     Integrated flux
  //   *Fspec      obs-frame spectrum (if OPT_SPEC==1)
  //   *FLAG_Finteg
  //        0 -> normal result; no issues
  //       99 -> zero flux
  //      666 -> undefined because model does not cover band wavelength
  //
  // !!! Dec 12 2018: Finteg is tested against SALT2 filter-fluxes,
  //     but Fspec is not tested.
  //
  // Jul 12 2019 RK - few fixes to work with spectrograph.
  //
  // Dec 8 2020: replace local magSmear[] with global GENSMEAR.MAGSMEAR_LIST
  //               (to work properly with G10 and C11 models)
  //
  // Apr 02 2021:
  //   + bug fix : multiply spectral flux by x0 (broadband fluxes were OK)
  //   + add FLAG_Finteg output arg.
  //
  // Sep 21 2022: fix awful bug fom Dec 2020:
  //    use GENSMEAR.MAGSMEAR_LIST[ilamobs] instead of 
  //    undefined/obsolsete magSmear[ilamobs]
  //

  int    ifilt          = IFILTMAP_SEDMODEL[ifilt_obs] ;
  int    NLAMFILT       = FILTER_SEDMODEL[ifilt].NLAM ;
  int   DO_SPECTROGRAPH = ( ifilt_obs == JFILT_SPECTROGRAPH ) ;

  double z1             = 1.0 + zHEL ;
  double Trest          = Tobs/z1 ;
  double minlam_filt    = FILTER_SEDMODEL[ifilt].lammin ;
  double maxlam_filt    = FILTER_SEDMODEL[ifilt].lammax ;
  double lamstep_filt   = FILTER_SEDMODEL[ifilt].lamstep ;
  double minlam_SED     = LAM[0];
  double maxlam_SED     = LAM[NLAM-1];
  double hc8            = (double)hc ;
  double MODELNORM_Finteg  = lamstep_filt / hc8 ;
  double MODELNORM_Fspec   = lamstep_filt ;
  double FLUXSUM_MIN       = 1.0E-30 ;

  int    ilamobs, ilamsed, ISTAT_SMEAR ;
  double TRANS, MWXT_FRAC, HOSTXT_FRAC;
  double LAMOBS, LAMSED, LAMSED_MIN, LAMSED_MAX;
  double LAMSED_STEP, LAMSPEC_STEP, LAMRATIO ;
  double TMPLAM[3], TMPSED[3];
  double lam[MXBIN_LAMFILT_SEDMODEL]; 
  double Fbin_forFlux, Fbin_forSpec, FTMP, arg, FSMEAR ;
  double Finteg_filter=0.0, Finteg_spec=0.0 ;

  //  int  LDMP = 0 ;
  char fnam[]  = "INTEG_zSED_PySEDMODEL" ;

  // ------------- BEGIN -----------


  *Finteg  = 0.0 ; // init output flux for filter
  Fspec[0] = 0.0 ;
  *FLAG_Finteg = 0 ;

  // first make sure that model SED wavelength range covers filter
  if ( !DO_SPECTROGRAPH ) {
    bool MODEL_COVERS_FILTER =
      ( minlam_filt >= minlam_SED*z1 && maxlam_filt <= maxlam_SED*z1 );
    if ( !MODEL_COVERS_FILTER )
      { *FLAG_Finteg = (int)MAG_UNDEFINED;  return ; }
  }


  // - - - - - - -
  // check for intrinsic scatter models in sntools_genSmear.c
  //  (e..g, G10, C11).  Get magSmear at all wavelengths.
  ISTAT_SMEAR = istat_genSmear(); // check for smear model

  if ( ISTAT_SMEAR ) {
    double cdum=0.0, x1dum=0.0 ;
    double parList[4] = { Trest, x1dum, cdum, -9.0 } ;
    for ( ilamobs=0; ilamobs < NLAMFILT; ilamobs++ ) {
      LAMOBS       = FILTER_SEDMODEL[ifilt].lam[ilamobs] ;
      LAMSED       = LAMOBS/z1;   // rest-frame wavelength
      lam[ilamobs] = LAMSED ;
    }

    /*
    printf(" xxx %s:  Trest=%.3f  ifilt_obs=%d \n",
	   fnam,  Trest, ifilt_obs); fflush(stdout);
    */
    get_genSmear(parList, NLAMFILT, lam, GENSMEAR.MAGSMEAR_LIST);
  }

  // - - - - - -

  LAMSED_STEP = lamstep_filt ;

  for ( ilamobs=0; ilamobs < NLAMFILT; ilamobs++ ) {

    TRANS  = FILTER_SEDMODEL[ifilt].transSN[ilamobs] ;

    if ( TRANS < 1.0E-12 && OPT_SPEC==0)   { continue ; }

    MWXT_FRAC  = SEDMODEL_TABLE_MWXT_FRAC[ifilt][ilamobs] ;

    // July 2016: check for host extinction.
    if( RV_host > 1.0E-9 && AV_host > 1.0E-9 )
      { HOSTXT_FRAC = SEDMODEL_TABLE_HOSTXT_FRAC[ifilt][ilamobs] ; }
    else
      { HOSTXT_FRAC = 1.0 ; } // standard SALT2 model has no host extinction

    LAMOBS     = FILTER_SEDMODEL[ifilt].lam[ilamobs] ;
    LAMSED     = LAMOBS / z1 ;           // rest-frame lambda
    LAMSED_MIN = LAMSED_MAX = LAMSED ;   // default is no sub-bins

    // check spectrum options
    if ( OPT_SPEC  && DO_SPECTROGRAPH )  {
      Finteg_spec = 0.0 ;
      // prepare sub-bins since SPECTROGRAPH bins can be large
      LAMSED_MIN = SPECTROGRAPH_SEDMODEL.LAMMIN_LIST[ilamobs]/z1 ;
      LAMSED_MAX = SPECTROGRAPH_SEDMODEL.LAMMAX_LIST[ilamobs]/z1 ;
    } // end OPT_SPEC

    // RK - Jul 12 2019 - make sure spectrograph bin is covered by SED model
    if ( LAMSED_MIN < minlam_SED ) { continue; }
    if ( LAMSED_MAX > maxlam_SED ) { continue; }

    // loop over rest-frame lambda (for SPECTROGRAPH)
    for(LAMSED = LAMSED_MIN; LAMSED <= LAMSED_MAX; LAMSED+=LAMSED_STEP ) {

      // find rest-frame bin for PySEDMODEL ... note that non-uniform
      // bins are allowed, but non-uniform might lead to trouble elsewhere.
      ilamsed = quickBinSearch(LAMSED, NLAM,LAM, fnam, fnam);
      if ( ilamsed >= NLAM-2 ) { ilamsed=NLAM-3; }

      TMPLAM[0]=LAM[ilamsed+0];  TMPSED[0]=SED[ilamsed+0];
      TMPLAM[1]=LAM[ilamsed+1];  TMPSED[1]=SED[ilamsed+1];
      TMPLAM[2]=LAM[ilamsed+2];  TMPSED[2]=SED[ilamsed+2];
      FTMP = quadInterp( LAMSED, TMPLAM, TMPSED, fnam);

      // check option to smear flux with intrinsic scatter (Apr 11 2019)
      if ( ISTAT_SMEAR ) {
	arg     = -0.4*GENSMEAR.MAGSMEAR_LIST[ilamobs];
	FSMEAR  =  pow(TEN,arg)  ;        // fraction change in flux
	FTMP   *=  FSMEAR;                // adjust flux for smearing
      }

      Fbin_forFlux = (FTMP * HOSTXT_FRAC*MWXT_FRAC * LAMSED*TRANS);
      Fbin_forSpec = (FTMP * HOSTXT_FRAC*MWXT_FRAC );
      Finteg_filter  +=  Fbin_forFlux ;

      if ( OPT_SPEC && DO_SPECTROGRAPH ) {
	if ( LAMSED+LAMSED_STEP < LAMSED_MAX )
	  { LAMSPEC_STEP = LAMSED_STEP * z1 ; } // obs-frame lamStep
	else
	  { LAMSPEC_STEP = (LAMSED_MAX-LAMSED)*z1 ; }

	LAMRATIO        = LAMSPEC_STEP/LAMSED_STEP ;
	Finteg_spec    += (Fbin_forSpec * LAMRATIO );

      } // end OPT_SPEC


    } // end loop over LAMSED

    // store spectrum
    if ( OPT_SPEC )
      { Fspec[ilamobs] = (Finteg_spec * x0 * MODELNORM_Fspec) ; }

  } // end ilamobs


  // - - - - - - -
  // store integrated flux in passband
  *Finteg = (Finteg_filter * x0 * MODELNORM_Finteg);

  if ( *Finteg < FLUXSUM_MIN ) { *FLAG_Finteg = (int)MAG_ZEROFLUX; }

  return;

} // end INTEG_zSED_PySEDMODEL


// ====================================================
void genSpec_PySEDMODEL(double Tobs, double zHEL, double MU,
			double MWEBV,                   // (I) galactic
			double RV_host, double AV_host, // (I) host
			int NHOSTPAR, double *HOSTPAR_LIST, // (I) host	
			double *GENFLUX_LIST,           // (O) fluxGen per bin
			double *GENMAG_LIST ) {         // (O) magGen per bin

  // March 2019
  // Return PySEDMODEL spectrum in SPECTROGRAPH bins
  //
  // Aug 26 2021: remove buggy 1/z1 factor

  double hc8   = (double)hc ;
  double z1    = 1.0 + zHEL ;
  double x0    = pow(TEN,-0.4*MU);
  int NBLAM    = SPECTROGRAPH_SEDMODEL.NBLAM_TOT ;
  int ilam, NLAM, FLAG_ignore ;
  int    NEWEVT_FLAG = 0 ;
  double Finteg_ignore, FTMP, MAG, ZP, LAM, Trest ;
  double *SED   = Event_PySEDMODEL.SED ;
  double *FLAM  = Event_PySEDMODEL.LAM;

  char fnam[] = "genSpec_PySEDMODEL" ;


  // --------- BEGIN ------------

  // get the spectrum
  Trest = Tobs/z1;
  fetchSED_PySEDMODEL(Event_PySEDMODEL.EXTERNAL_ID, NEWEVT_FLAG, Trest,
		      MXLAM_PySEDMODEL, HOSTPAR_LIST, &NLAM, FLAM, SED);
  Event_PySEDMODEL.NLAM = NLAM ;


  // init entire spectum to zero.
  for(ilam=0; ilam < NBLAM; ilam++ ) 
    { GENFLUX_LIST[ilam] = 0.0 ; }

  INTEG_zSED_PySEDMODEL(1, JFILT_SPECTROGRAPH, Tobs, zHEL, x0,
			RV_host, AV_host,
			Event_PySEDMODEL.NLAM,
			Event_PySEDMODEL.LAM,
			Event_PySEDMODEL.SED,
			&Finteg_ignore, GENFLUX_LIST, // <= returned
			&FLAG_ignore );

  // convert generated fluxes into mags
  for(ilam=0; ilam < NBLAM; ilam++ ) {
    LAM  = SPECTROGRAPH_SEDMODEL.LAMAVG_LIST[ilam] ;
    ZP   = SPECTROGRAPH_SEDMODEL.ZP_LIST[ilam] ;
    FTMP = ( LAM/ hc8 ) * GENFLUX_LIST[ilam] ;
    if ( ZP > ZPMIN_SPECTROGRAPH && FTMP > 0.0 )   {
      MAG = -2.5*log10(FTMP) + ZP ;
    }
    else  {
      MAG = MAG_UNDEFINED ;  // model undefined
    }
    GENMAG_LIST[ilam] = MAG ;
  } // end ilam loop over SPECTROGRAPH bins


  return ;

} // end genSpec_PySEDMODEL


// =====================================
void  read_SALT2_template0(void) {

  // read M0 surface of SALT2 model, corresponding to x1=c=0.

  char sedcomment[80] ;
  char SALT2_tempate0_file[MXPATHLEN];

  int nflux_nan;
  double Trange[2] = { -40.0,  100.0   } ;
  double Lrange[2] = { 2000.0, 9200.0  } ;
  char   MODELNAME_SALT2[] = "SALT2.JLA-B14" ;
  //  char fnam[] = "read_SALT2_template0" ;

  // ------------- BEGIN -------------

  sprintf(SALT2_tempate0_file,
	  "%s/models/SALT2/%s/salt2_template_0.dat",
	  PATH_SNDATA_ROOT, MODELNAME_SALT2 );

  malloc_SEDFLUX_SEDMODEL(&TEMP_SEDMODEL,0,0,0);
  sprintf(sedcomment,"template-0(%s)", MODELNAME_SALT2);

  rd_sedFlux(SALT2_tempate0_file, sedcomment, Trange, Lrange
	     ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL, 0
	     ,&TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY, &TEMP_SEDMODEL.DAYSTEP
	     ,&TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM, &TEMP_SEDMODEL.LAMSTEP
	     ,TEMP_SEDMODEL.FLUX,  TEMP_SEDMODEL.FLUXERR
	     ,&nflux_nan );

  return ;

} // end read_SALT2_template0
