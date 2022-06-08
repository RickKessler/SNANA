/*****************************************
  Created Sep 2018
  Sep 10 2020: rename BYOSED -> more generic name PySEDMODEL

  C-wrapper to call python function that returns rest-frame SED,
  snd then C functions compute & return observer-frame magnitudes
  to the simulation. Also returns obs-frame spectra if requested.

  Each PySEDMODEL is associated with a separate gensed_[model].py :
     gensed_BYOSED.py : Build Your Own SED  (J.Pierel)
     gensed_SNEMO.py  : SNFactory model (Ben Rose)
     gensed_BAYESN.py : BayeSN model (Gautham Narayan, Stephen Thorp, Kaisey Mandel)

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
PyObject *geninit_PySEDMODEL ;

//int init_numpy(){
//  import_array(); // PyError if not successful
//  return 0;
//}

#endif

// ===============================================

void load_PySEDMODEL_CHOICE_LIST(void) {

  int N=0;
  char fnam[] = "load_PySEDMODEL_CHOICE_LIST" ;

  // generic utility to store all possible PySEDMODEL names.
  // Used by sim, parsing, etc ...
  sprintf(PySEDMODEL_CHOICE_LIST[N], "%s", MODEL_NAME_BYOSED); N++ ;
  sprintf(PySEDMODEL_CHOICE_LIST[N], "%s", MODEL_NAME_SNEMO ); N++ ;
  sprintf(PySEDMODEL_CHOICE_LIST[N], "%s", MODEL_NAME_BAYESN ); N++ ;

  if ( N != NCHOICE_PySEDMODEL ) {
    sprintf(c1err,"Expected %d PySEDMODEL choices");
    sprintf(c2err,"but loaded %d", N);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
} // end load_PySEDMODEL_CHOICE_LIST


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
  PyObject *genmod, *genclass, *pargs;
#endif

  char *PyMODEL_NAME = INPUTS_PySEDMODEL.MODEL_NAME ;
  char *PyFUN_NAME   = INPUTS_PySEDMODEL.PyFUN_NAME ;
  int  L, ipar, NPAR ;
  int  MEMD   = sizeof(double);
  int  MEMC   = sizeof(char);
  char comma[] = ",";
  char fnam[] = "init_genmag_PySEDMODEL" ;

  // -------------- BEGIN ------------------

  sprintf(BANNER, "%s", fnam);
  print_banner(BANNER);

  sprintf(PyMODEL_NAME, "%s",      MODEL_NAME);
  sprintf(PyFUN_NAME, "gensed_%s", PyMODEL_NAME) ;

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

  splitString(NAMES_HOSTPAR, comma, MXHOSTPAR_PySEDMODEL,
	      &NPAR, INPUTS_PySEDMODEL.NAME_ARRAY_HOSTPAR );

  // - - - - - - - - - - -
  // print summary of filter info
  filtdump_SEDMODEL();

  // init a few C struct items
  Event_PySEDMODEL.LAST_EXTERNAL_ID = -9;
  Event_PySEDMODEL.LAM  = (double*) malloc( MXLAM_PySEDMODEL*MEMD ) ;
  Event_PySEDMODEL.SED  = (double*) malloc( MXLAM_PySEDMODEL*MEMD ) ;

  SEDMODEL_MWEBV_LAST     = -999.   ;
  SEDMODEL_HOSTXT_LAST.AV = -999.   ;
  SEDMODEL_HOSTXT_LAST.z  = -999.   ;

#ifndef USE_PYTHON
  printf("\n no python ==> read SALT2 template with C code. \n");
  // read M0 surface of SALT2 model, corresponding to x1=c=0.
  // This is for debug-comparisons with SALT2 model where x1=c=0.
  read_SALT2_template0();

#endif


#ifdef USE_PYTHON
  printf("\t Begin %s python-init from C code ... \n", PyMODEL_NAME );
  Py_Initialize();
  int nResult1 = PyRun_SimpleStringFlags("import numpy", NULL);
  int nResult2 = PyRun_SimpleStringFlags("import os", NULL);
  int nResult3 = PyRun_SimpleStringFlags("import optparse",NULL);
  int nResult4 = PyRun_SimpleStringFlags("import configparser",NULL);

  // xxxx  genmod = PyImport_ImportModule("genmag_BYOSED");
  printf("DEBUG", PyFUN_NAME, "\n");
  genmod = PyImport_ImportModule(PyFUN_NAME);

  if (genmod == NULL) {
    sprintf(c1err,"Could not import module %s", PyFUN_NAME);
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // xxxx  genclass = PyObject_GetAttrString(genmod, "genmag_BYOSED");
  genclass = PyObject_GetAttrString(genmod, PyFUN_NAME);
  if (genclass == NULL) {
    sprintf(c1err,"Could not import PyObject_GetAttrString class");
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  pargs  = Py_BuildValue("(siss)",PATH_VERSION,OPTMASK,ARGLIST,NAMES_HOSTPAR);
  geninit_PySEDMODEL = PyEval_CallObject(genclass, pargs);
  if (geninit_PySEDMODEL == NULL) {
    sprintf(c1err,"Could not run PyEval_CallObject module");
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  Py_DECREF(genmod);
  Py_DECREF(genclass);
  Py_DECREF(pargs);

  printf("\t Finished %s python-init from C code \n", PyMODEL_NAME );
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

  return ;

} // end init_genmag_PySEDMODEL

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

// =========================================================
void genmag_PySEDMODEL(int EXTERNAL_ID, double zHEL, double zCMB, double MU,
		       double MWEBV, int NHOSTPAR, double *HOSTPAR_LIST,
		       int IFILT_OBS, int NOBS, double *TOBS_list,
		       double *MAGOBS_list, double *MAGERR_list ) {

  // Created Sep 2018
  //
  // Inputs:
  //   EXTERNAL_ID  : set NEWEVT_FLAG logical when this changes
  //   zHEL         : helio redshift
  //   MU           : distance modulus
  //   MWEBV        : E(B-V) for Milky Wat
  //   NHOSTPAR     : Number of host params in HOSTPAR_LIST
  //   HOSTPAR_LIST : RV, AV, LOGMASS, SFR ...
  //   IFILT_OBS    : absolute filter index
  //   NOBS         : number of observations
  //   TOBS_list    : list of MJD-PEAKMJD
  //
  // Outputs"
  //   MAGOBS_list   : list of true mags
  //   MAGERR_list   : list of mag errors (place-holder, in case)
  //
  // Feb 20 2021 RK
  //   + fix bug that was preventing call to fetchParVal_PySEDMODEL()

  char *MODEL_NAME = INPUTS_PySEDMODEL.MODEL_NAME ;
  char *PyFUN_NAME = INPUTS_PySEDMODEL.PyFUN_NAME ;

  double RV_host = HOSTPAR_LIST[0];
  double AV_host = HOSTPAR_LIST[1];
  double z1    = 1.0 + zHEL ;
  double *LAM  = Event_PySEDMODEL.LAM;
  double *SED  = Event_PySEDMODEL.SED ;

  int  ifilt  = IFILTMAP_SEDMODEL[IFILT_OBS] ; // sparse filter index
  double  ZP  = FILTER_SEDMODEL[ifilt].ZP ;    // ZP for flux->mag
  double x0   = pow(10.0,-0.4*MU);             // dimming from dist. mod.
  int    NEWEVT_FLAG = 0 ;
  int    NEWEVT_FLAG_TMP ;
  int    DUMPFLAG_HOSTPAR = 0 ;
  int    FLAG_Finteg;

  int    NLAM, o, ipar ;
  double Tobs, Trest, FLUXSUM_OBS, FspecDUM[2], magobs ;
  char   pyFORMAT_STRING_HOSTPAR[100] ;;
  char fnam[] = "genmag_PySEDMODEL" ;

   #ifdef USE_PYTHON
  // python declarations here

   #endif

  // ------------ BEGIN -----------

  for(o=0; o < NOBS; o++ ) { MAGOBS_list[o] = 99.0 ; }

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

  // construct hostpar string to pass to python
  sprintf(pyFORMAT_STRING_HOSTPAR,"diii[" );
  for(ipar=0; ipar < NHOSTPAR; ipar++ ) {
    strcat(pyFORMAT_STRING_HOSTPAR,"d");
    if ( ipar < NHOSTPAR-1 )
      { strcat(pyFORMAT_STRING_HOSTPAR,","); }
    else
      { strcat(pyFORMAT_STRING_HOSTPAR,"]"); }
  }
  //  printf(" xxx pySTRING_HOSTPAR = '%s' \n", pyFORMAT_STRING_HOSTPAR );


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

  for(o=0; o < NOBS; o++ ) {
    Tobs  = TOBS_list[o];
    Trest = Tobs/z1;
    if (o == 0 )
      { NEWEVT_FLAG_TMP = NEWEVT_FLAG; }
    else
      { NEWEVT_FLAG_TMP = 0; }

    fetchSED_PySEDMODEL(EXTERNAL_ID, NEWEVT_FLAG_TMP, Trest,
			MXLAM_PySEDMODEL, HOSTPAR_LIST, &NLAM, LAM, SED,
			pyFORMAT_STRING_HOSTPAR);
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

    MAGOBS_list[o] = magobs;  // load output array
    MAGERR_list[o] = 0.01;    // not used
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
  char pyfun_tmp[80];
  int NPAR = 0 ;
  //  char fnam[] = "fetchParNames_PySEDMODEL" ;

#ifdef USE_PYTHON
  // python declarations here

  int ipar;
  PyObject *parnamesmeth,*pNames,*pNPARmeth,*pNPAR,*pnamesitem;
  PyListObject *arrNames;
  printf("fetching parameter names from Python\n");
  // David: need your python magic to return these string names.

  sprintf(pyfun_tmp, "fetchParNames_%s", MODEL_NAME);

  parnamesmeth  = PyObject_GetAttrString(geninit_PySEDMODEL, pyfun_tmp);
  //xx  parnamesmeth  = PyObject_GetAttrString(geninit_PySEDMODEL,
  //xx					 "fetchParNames_BYOSED");

  sprintf(pyfun_tmp, "fetchNParNames_%s", MODEL_NAME );

  pNPARmeth  = PyObject_GetAttrString(geninit_PySEDMODEL, pyfun_tmp);
  // xxx  pNPARmeth  = PyObject_GetAttrString(geninit_PySEDMODEL,
  // xxx			      "fetchNParNames_BYOSED");

  pNames  = PyEval_CallObject(parnamesmeth, NULL);
  pNPAR  = PyEval_CallObject(pNPARmeth, NULL);

  NPAR = PyLong_AsLong(pNPAR);
  arrNames = (PyListObject *)(pNames);
  for(ipar=0; ipar < NPAR; ipar++ ) {
    pnamesitem = PyList_GetItem(arrNames,ipar);
    sprintf(parNameList[ipar],PyUnicode_AsUTF8(pnamesitem));
  }

  Py_DECREF(arrNames);
  Py_DECREF(pNames);
  Py_DECREF(parnamesmeth);
  Py_DECREF(pNPARmeth);
  Py_DECREF(pNPAR);

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
  char **parNameList, pyfun_tmp[60] ;
  int NPAR, ipar;
  char fnam[] = "fetchParVal_PySEDMODEL" ;

  // ------------- BEGIN ------------------

  NPAR = Event_PySEDMODEL.NPAR;
  parNameList = Event_PySEDMODEL.PARNAME;
  // David: need python function to return these values.
#ifdef USE_PYTHON

  sprintf(pyfun_tmp, "fetchParVals_%s_4SNANA", MODEL_NAME );
  parvalmeth  = PyObject_GetAttrString(geninit_PySEDMODEL, pyfun_tmp);
  // xxx  parvalmeth  = PyObject_GetAttrString(geninit_PySEDMODEL,
  // xxx			       "fetchParVals_BYOSED_4SNANA");

  for(ipar=0; ipar < NPAR; ipar++ ) {
    pargs  = Py_BuildValue("(s)",parNameList[ipar]);
    pParVal  = PyEval_CallObject(parvalmeth, pargs);
    val = PyFloat_AsDouble(pParVal);
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
			 double *LAM_SED, double *FLUX_SED,
			 char *pyFORMAT_STRING_HOSTPAR) {

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
  char pyfun_tmp[60];
  char fnam[] = "fetchSED_PySEDMODEL" ;

  // ------------ BEGIN -----------

  *NLAM_SED = 0 ; // init output

#ifdef USE_PYTHON
  PyObject *pmeth, *pargs, *pargs2, *pNLAM, *pLAM, *pFLUX, *plammeth, *pnlammeth;
  int NLAM, ilam, ihost;
  PyListObject *arrLAM, *arrFLUX;
  PyObject *pylamitem, *pyfluxitem;
  //int numpy_initialized =  init_numpy();

  // python declarations here
  sprintf(pyfun_tmp, "fetchSED_%s", MODEL_NAME );
  pmeth  = PyObject_GetAttrString(geninit_PySEDMODEL, pyfun_tmp);

  // xxx  pmeth  = PyObject_GetAttrString(geninit_PySEDMODEL,
  // xxx			  "fetchSED_BYOSED"); //

  plammeth  = PyObject_GetAttrString(geninit_PySEDMODEL, "fetchSED_LAM");
  pnlammeth = PyObject_GetAttrString(geninit_PySEDMODEL, "fetchSED_NLAM");

  pargs = PyTuple_New(5);
  pargs2 = PyTuple_New(sizeof(HOSTPAR_LIST));
  PyTuple_SetItem(pargs,0,PyFloat_FromDouble(Trest));
  PyTuple_SetItem(pargs,1,PyFloat_FromDouble(MXLAM));
  PyTuple_SetItem(pargs,2,PyFloat_FromDouble(EXTERNAL_ID));
  PyTuple_SetItem(pargs,3,PyFloat_FromDouble(NEWEVT_FLAG));

  for(ihost=0; ihost < sizeof(HOSTPAR_LIST); ihost++ ){
    PyTuple_SetItem(pargs2,ihost,PyFloat_FromDouble(HOSTPAR_LIST[ihost]));
  }

  PyTuple_SetItem(pargs,4,pargs2);
  pNLAM  = PyEval_CallObject(pnlammeth, NULL);
  pLAM  = PyEval_CallObject(plammeth, NULL);
  pFLUX   = PyEval_CallObject(pmeth, pargs);

  Py_DECREF(pmeth);
  Py_DECREF(plammeth);
  Py_DECREF(pnlammeth);

  NLAM = PyFloat_AsDouble(pNLAM);
  Py_DECREF(pNLAM);

  arrLAM  = (PyListObject *)(pLAM);
  arrFLUX = (PyListObject *)(pFLUX);
  for(ilam=0; ilam < NLAM; ilam++ ) {
    // interpolate flux to Trest
    pylamitem  = PyList_GetItem(arrLAM,ilam);
    pyfluxitem = PyList_GetItem(arrFLUX,ilam);

    LAM_SED[ilam]  = PyFloat_AsDouble(pylamitem);
    FLUX_SED[ilam] = PyFloat_AsDouble(pyfluxitem);
  }

  *NLAM_SED = NLAM;

  Py_DECREF(pLAM);
  Py_DECREF(pFLUX);
  Py_DECREF(arrLAM);
  //Py_DECREF(arrFLUX);
  Py_DECREF(pargs);
  Py_DECREF(pargs2);
  //Py_DECREF(pylamitem);
  //Py_DECREF(pyfluxitem);

#endif


#ifndef USE_PYTHON
  // C-code return SALT2 SED interpolated to Trest
  int NLAM, NDAY, ilam, iday, jf0, jf1, jf2 ;
  double TMPDAY[3], TMPSED[3],  F_interp, FSCALE ;

  NLAM = TEMP_SEDMODEL.NLAM ;    *NLAM_SED = NLAM ;
  NDAY = TEMP_SEDMODEL.NDAY ;
  iday = quickBinSearch(Trest, NDAY,TEMP_SEDMODEL.DAY, fnam);
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
  double lam[MXBIN_LAMFILT_SEDMODEL], magSmear[MXBIN_LAMFILT_SEDMODEL];
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
      ilamsed = quickBinSearch(LAMSED, NLAM,LAM, fnam);
      if ( ilamsed >= NLAM-2 ) { ilamsed=NLAM-3; }

      TMPLAM[0]=LAM[ilamsed+0];  TMPSED[0]=SED[ilamsed+0];
      TMPLAM[1]=LAM[ilamsed+1];  TMPSED[1]=SED[ilamsed+1];
      TMPLAM[2]=LAM[ilamsed+2];  TMPSED[2]=SED[ilamsed+2];
      FTMP = quadInterp( LAMSED, TMPLAM, TMPSED, fnam);

      // check option to smear flux with intrinsic scatter (Apr 11 2019)
      if ( ISTAT_SMEAR ) {
	arg     =  -0.4*magSmear[ilamobs] ;
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
  int ilam, FLAG_ignore ;
  double Finteg_ignore, FTMP, MAG, ZP, LAM ;
  //  char fnam[] = "genSpec_PySEDMODEL" ;

  // --------- BEGIN ------------

  // init entire spectum to zero.
  for(ilam=0; ilam < NBLAM; ilam++ ) { GENFLUX_LIST[ilam] = 0.0 ; }

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
    if ( ZP > 0.0 && FTMP > 0.0 )   {
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
	     ,TEMP_SEDMODEL.FLUX,  TEMP_SEDMODEL.FLUXERR );

  return ;

} // end read_SALT2_template0
