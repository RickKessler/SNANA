/*****************************************
  Created Sep 2018
  BYOSED = Build Your Own SED

  Lots of options to mangle and tweak an inital SED sequence
  such as Hsiao. Example options are to add random stretch,
  apply color law, add spectral features, correlate with
  host properties, etc ...

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

 *****************************************/


#include  <stdio.h> 
#include  <math.h>     
#include  <stdlib.h>   
#include  <sys/stat.h>


#include  "sntools.h"           // SNANA community tools
#include  "sntools_genSmear.h"
#include  "sntools_spectrograph.h"
#include  "MWgaldust.h"
#include  "genmag_SEDtools.h"
#include  "genmag_SIMSED.h"
#include  "genmag_BYOSED.h"

#ifdef USE_PYTHON
#include  <Python.h>
//#include <numpy/arrayobject.h>
//#include <numpy/ndarrayobject.h>
PyObject *geninit_BYOSED;

//int init_numpy(){
//  import_array(); // PyError if not successful
//  return 0;
//}

#endif

// =========================================================
void init_genmag_BYOSED(char *PATH_VERSION, int OPTMASK, char *ARGLIST ) {

  // Read input directory file(s) for parameters characterizing 
  // how to build your SED.
  //  
  // Inputs:
  //  PATH_VERSION : points to model-param directory;
  //               :  passed from GENMODEL arg of sim-input
  //  OPTMASK      : bit mask of options; interpreted by python code.
  //               : OPTMASK=-1 is a flag to print options.
  //               :  passed from GENMODEL_MSKOPT arg in sim-input
  //  
  //  ARGLIST      : string of options
  //

#ifdef USE_PYTHON
  PyObject *genmod, *genclass, *pargs;
#endif

  int  MEMD   = sizeof(double);
  int  MEMC   = sizeof(char);
  char fnam[] = "init_genmag_BYOSED" ;

  // -------------- BEGIN ------------------

  sprintf(BANNER, "%s", fnam);
  print_banner(BANNER);
  printf("   BYOSED PATH    = '%s' \n",  PATH_VERSION);
  printf("   BYOSED OPTMASK = %d \n",    OPTMASK );	
  printf("   BYOSED ARGLIST = '%s' \n",  ARGLIST );	
  fflush(stdout);

  // print summary of filter info
  filtdump_SEDMODEL();

  // init a few C struct items
  Event_BYOSED.LAST_EXTERNAL_ID = -9;
  Event_BYOSED.LAM  = (double*) malloc( MXLAM_BYOSED*MEMD ) ;
  Event_BYOSED.SED  = (double*) malloc( MXLAM_BYOSED*MEMD ) ;

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
  printf("\t Begin python-init from C code ... \n");
  Py_Initialize();
  int nResult1 = PyRun_SimpleStringFlags("import numpy", NULL);
  int nResult2 = PyRun_SimpleStringFlags("import os", NULL);
  int nResult3 = PyRun_SimpleStringFlags("import optparse",NULL);
  int nResult4 = PyRun_SimpleStringFlags("import configparser",NULL);
  genmod = PyImport_ImportModule("genmag_BYOSED");
  if (genmod == NULL) {
    sprintf(c1err,"Could not import class genmag_BYOSED");
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  genclass = PyObject_GetAttrString(genmod, "genmag_BYOSED");
  if (genclass == NULL) {
    sprintf(c1err,"Could not import PyObject_GetAttrString module");
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  pargs  = Py_BuildValue("(sis)",PATH_VERSION,OPTMASK,ARGLIST);
  geninit_BYOSED = PyEval_CallObject(genclass, pargs);
  if (geninit_BYOSED == NULL) {
    sprintf(c1err,"Could not run PyEval_CallObject module");
    sprintf(c2err,"2nd message ??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  Py_DECREF(genmod);
  Py_DECREF(genclass);
  Py_DECREF(pargs);

  printf("\t Finished python-init from C code \n");
#endif
  
  // -----------------------------------------------------
  // set SED par names and allocate arrays for parameters
  int NPAR, ipar;
  Event_BYOSED.PARVAL  = (double*) malloc ( MXPAR_BYOSED*MEMD );
  Event_BYOSED.PARNAME = (char**)  malloc ( MXPAR_BYOSED*sizeof(char*) );
  Event_BYOSED.NPAR    = 0 ;
  for(ipar=0; ipar < MXPAR_BYOSED; ipar++ ) 
    { Event_BYOSED.PARNAME[ipar] = (char*) malloc(40*MEMC);  }

#ifdef USE_PYTHON
  NPAR = fetchParNames_BYOSED(Event_BYOSED.PARNAME);
  Event_BYOSED.NPAR = NPAR;
  printf("\t BYOSED parameters to store in data files:\n");
  for(ipar=0; ipar < NPAR; ipar++ ) 
    { printf("\t\t %s \n", Event_BYOSED.PARNAME[ipar] ); }
#endif

  // - - - - 
  printf("\n\t Done with %s \n", fnam);

  return ;

} // end init_BYOSED

// =========================================================
void genmag_BYOSED(int EXTERNAL_ID, double zHEL, double MU, 
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

  double RV_host = HOSTPAR_LIST[0];
  double AV_host = HOSTPAR_LIST[1];
  double FLUXSUM_MIN = 1.0E-30 ;
  double z1    = 1.0 + zHEL ;
  double *LAM  = Event_BYOSED.LAM;
  double *SED  = Event_BYOSED.SED ;

  int  ifilt  = IFILTMAP_SEDMODEL[IFILT_OBS] ; // sparse filter index
  char *cfilt = FILTER_SEDMODEL[ifilt].name ;
  double  ZP  = FILTER_SEDMODEL[ifilt].ZP ;    // ZP for flux->mag
  double x0   = pow(10.0,-0.4*MU);             // dimming from dist. mod.
  int    NEWEVT_FLAG = 0 ;

  int    NLAM, o, ipar ;
  double Tobs, Trest, FLUXSUM_OBS, FspecDUM[2], magobs ; 
  char pySTRING_HOSTPAR[100], dSTRING_HOSTPAR[100], ctmp[20];
  char fnam[] = "genmag_BYOSED" ;

   #ifdef USE_PYTHON
  // python declarations here
  
   #endif

  // ------------ BEGIN -----------

  for(o=0; o < NOBS; o++ ) { MAGOBS_list[o] = 99.0 ; }

  // check of this is a new event, or same event
  // with different epoch
  if ( EXTERNAL_ID != Event_BYOSED.LAST_EXTERNAL_ID )
    { NEWEVT_FLAG=1; }


  // construct hostpar string to pass to python
  dSTRING_HOSTPAR[0] = 0 ;
  for(ipar=0; ipar < NHOSTPAR; ipar++ ) {
    sprintf(ctmp,"%f", HOSTPAR_LIST[ipar] );
    strcat(dSTRING_HOSTPAR,ctmp);
    if ( ipar < NHOSTPAR-1 ) { strcat(dSTRING_HOSTPAR,","); }
  }
  sprintf(pySTRING_HOSTPAR,"diii[%s]", dSTRING_HOSTPAR);
  //  printf(" xxx pySTRING_HOSTPAR = '%s' \n", pySTRING_HOSTPAR );
  
  /* xxx
  printf(" xxx ------------------------------------ \n" ) ;
  printf(" xxx %s: process z=%.3f MU=%.3f RV=%3.1f IFILT_OBS=%d(%s) \n",
	 fnam, zHEL, MU, MWXT_SEDMODEL.RV, IFILT_OBS, cfilt );
  printf(" xxx RVMW=%3.1f  MWEBV=%.3f      RV_host=%3.1f AV=%4.2f \n", 
	 MWXT_SEDMODEL.RV, MWEBV,   RV_host, AV_host );
  xxxx */


  
  // store table info for Galactic & host extinction
  fill_TABLE_MWXT_SEDMODEL(MWXT_SEDMODEL.RV, MWEBV);
  fill_TABLE_HOSTXT_SEDMODEL(RV_host, AV_host, zHEL);   // July 2016

  for(o=0; o < NOBS; o++ ) {
    Tobs  = TOBS_list[o];
    Trest = Tobs/z1;

    fetchSED_BYOSED(EXTERNAL_ID, NEWEVT_FLAG, Trest, 
		    MXLAM_BYOSED, HOSTPAR_LIST, &NLAM, LAM, SED );  
    Event_BYOSED.NLAM = NLAM ;

    // integrate redshifted SED to get observer-frame flux in IFILT_OBS band.
    // FLUXSUM_OBS is returned (ignore FspecDUM)
    INTEG_zSED_BYOSED(0, IFILT_OBS, Tobs, zHEL, x0,RV_host,AV_host, 
		      NLAM, LAM, SED, 
		      &FLUXSUM_OBS, FspecDUM ); // <= returned 

    
    // convert calibrated flux into true magnitude
    if ( FLUXSUM_OBS > FLUXSUM_MIN ) 
      { magobs = ZP - 2.5*log10(FLUXSUM_OBS); }
    else
      { magobs = MAG_ZEROFLUX ; }

    MAGOBS_list[o] = magobs;  // load output array
    MAGERR_list[o] = 0.01;    // not used
  }

  // for NEW EVENT, store SED parameters so that sim can 
  // write them to data files
  // hack
  if ( NEWEVT_FLAG ) { 
    fetchParVal_BYOSED(Event_BYOSED.PARVAL); 
  }

  // keep track of last ID
  Event_BYOSED.LAST_EXTERNAL_ID = EXTERNAL_ID ;

  return ;


} // end genmag_BYOSED


// ============================================================
//  FETCH UTILITIES TO RETURN EXTRA INFO TO MAIN PROGRAM
// =============================================================


// ================================================
int fetchParNames_BYOSED(char **parNameList) {

  // Pass name of each parameter to calling function, so that
  // these parameters can be included in the data files.
  // Function Returns number of parameters used to create
  // each SED.  **parNameList is a list of parameter names.
  //
  // Called once during init stage.

  int NPAR,ipar;
  char fnam[] = "fetchParNames_BYOSED" ;
#ifdef USE_PYTHON
  // python declarations here
  
  PyObject *parnamesmeth,*pNames,*pNPARmeth,*pNPAR,*pnamesitem;
  PyListObject *arrNames;
  printf("fetching parameter names from Python\n");
  // David: need your python magic to return these string names.
  parnamesmeth  = PyObject_GetAttrString(geninit_BYOSED, "fetchParNames_BYOSED");
  pNPARmeth  = PyObject_GetAttrString(geninit_BYOSED, "fetchNParNames_BYOSED");
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

} // fetchParNames_BYOSED


void fetchParVal_BYOSED(double *parVal) {

  // return list of parameters to calling function (sim)
  // so that these parameters can be included in the
  // data files.
  //
  // Called once per event.

#ifdef USE_PYTHON  
  PyObject *parvalmeth,*pParVal,*pargs;
#endif
  double val;
  char **parNameList;
  int NPAR,ipar;
  char fnam[] = "fetchParVal_BYOSED" ;

  // ------------- BEGIN ------------------

  NPAR = Event_BYOSED.NPAR; //fetchParNames_BYOSED(parNameList);
  parNameList = Event_BYOSED.PARNAME;
  // David: need python function to return these values.
#ifdef USE_PYTHON
  parvalmeth  = PyObject_GetAttrString(geninit_BYOSED, "fetchParVals_BYOSED_4SNANA");

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
  
  return ;

} // end fetchParVal_BYOSED

// =================================================
void fetchSED_BYOSED(int EXTERNAL_ID, int NEWEVT_FLAG, double Trest, int MXLAM,
		     double *HOSTPAR_LIST, int *NLAM_SED, double *LAM_SED, double *FLUX_SED) {

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

  char fnam[] = "fetchSED_BYOSED" ;

  // ------------ BEGIN -----------

  *NLAM_SED = 0 ; // init output

#ifdef USE_PYTHON
  PyObject *pmeth, *pargs, *pNLAM, *pLAM, *pFLUX, *plammeth, *pnlammeth;
  int NLAM, ilam;
  PyListObject *arrLAM, *arrFLUX;
  PyObject *pylamitem, *pyfluxitem;
  //int numpy_initialized =  init_numpy();
  
  // python declarations here
  pmeth  = PyObject_GetAttrString(geninit_BYOSED, "fetchSED_BYOSED");
  plammeth  = PyObject_GetAttrString(geninit_BYOSED, "fetchSED_LAM");
  pnlammeth  = PyObject_GetAttrString(geninit_BYOSED, "fetchSED_NLAM");
  pargs  = Py_BuildValue("diii",Trest,MXLAM,EXTERNAL_ID,NEWEVT_FLAG);

  pNLAM  = PyEval_CallObject(pnlammeth, NULL);
  pLAM  = PyEval_CallObject(plammeth, NULL);
  pFLUX   = PyEval_CallObject(pmeth, pargs);

  Py_DECREF(pmeth);
  Py_DECREF(plammeth);
  Py_DECREF(pnlammeth);
  
  NLAM = PyFloat_AsDouble(pNLAM);
  Py_DECREF(pNLAM);
  
  arrLAM = (PyListObject *)(pLAM);
  arrFLUX = (PyListObject *)(pFLUX);
  for(ilam=0; ilam < NLAM; ilam++ ) {
    // interpolate flux to Trest
    pylamitem = PyList_GetItem(arrLAM,ilam);
    pyfluxitem = PyList_GetItem(arrFLUX,ilam);

    LAM_SED[ilam]  = PyFloat_AsDouble(pylamitem);
    FLUX_SED[ilam] = PyFloat_AsDouble(pyfluxitem);
  }

  *NLAM_SED = NLAM;
  
  Py_DECREF(pLAM);
  Py_DECREF(pFLUX);
  Py_DECREF(arrLAM);
  Py_DECREF(arrFLUX);
  Py_DECREF(pargs);
  //Py_DECREF(pylamitem);
  //Py_DECREF(pyfluxitem);

#endif


#ifndef USE_PYTHON
  // C-code return SALT2 SED interpolated to Trest
  int NLAM, NDAY, ilam, iday, jf0, jf1, jf2 ;  
  double TMPDAY[3], TMPSED[3],  F_interp, FSCALE ;
  
  NLAM = TEMP_SEDMODEL.NLAM ;    *NLAM_SED = NLAM ;
  NDAY = TEMP_SEDMODEL.NDAY ;
  iday = quickBinSearch(NDAY, Trest, TEMP_SEDMODEL.DAY, fnam);
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

} // end fetchSED_BYOSED


// =====================================================
void INTEG_zSED_BYOSED(int OPT_SPEC, int ifilt_obs, double Tobs, 
		       double zHEL, double x0, 
		       double RV_host, double AV_host,
		       int NLAM, double *LAM, double *SED,
		       double *Finteg, double *Fspec ) {

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
  //  *SED          array of rest-fame SED fluxes
  //
  // Outputs
  //   *Finteg     Integrated flux 
  //   *Fspec      obs-frame spectrum (if OPT_SPEC==1)
  //
  // !!! Dec 12 2018: Finteg is tested against SALT2 filter-fluxes, 
  //     but Fspec is not tested.
  //

  int    ifilt          = IFILTMAP_SEDMODEL[ifilt_obs] ;
  int    NLAMFILT       = FILTER_SEDMODEL[ifilt].NLAM ;
  int   DO_SPECTROGRAPH = ( ifilt_obs == JFILT_SPECTROGRAPH ) ;
  double meanlam_obs    = FILTER_SEDMODEL[ifilt].mean ;  // mean lambda
  double z1             = 1.0 + zHEL ;
  double Trest          = Tobs/z1 ;
  double meanlam_rest   = meanlam_obs/z1 ;
  double minlam_filt    = FILTER_SEDMODEL[ifilt].minlam ;
  double maxlam_filt    = FILTER_SEDMODEL[ifilt].maxlam ;
  double lamstep_filt   = FILTER_SEDMODEL[ifilt].lamstep ;
  double minlam_SED     = LAM[0];
  double maxlam_SED     = LAM[NLAM-1];
  double hc8            = (double)hc ;
  double MODELNORM_Finteg  = lamstep_filt / hc8 ;
  double MODELNORM_Fspec   = lamstep_filt ;

  int    LABORT, ilamobs, ilamsed, ISTAT_SMEAR ;
  double TRANS, MWXT_FRAC, HOSTXT_FRAC, FLUXSUM = 0.0 ;
  double LAMOBS, LAMSED, LAMSED_MIN, LAMSED_MAX;
  double LAMSED_STEP, LAMSPEC_STEP, LAMRATIO ;
  double TMPLAM[3], TMPSED[3];
  double lam[MXBIN_LAMFILT_SEDMODEL], magSmear[MXBIN_LAMFILT_SEDMODEL];
  double Fbin_forFlux, Fbin_forSpec, FTMP, arg, FSMEAR ; 
  double Finteg_filter=0.0, Finteg_spec=0.0 ;
  char *cfilt  = FILTER_SEDMODEL[ifilt].name ;

  int  LDMP = 0 ;
  char fnam[]  = "INTEG_zSED_BYOSED" ;

  // ------------- BEGIN -----------

  *Finteg  = 0.0 ; // init output flux for filter
  Fspec[0] = 0.0 ;

  // first make sure that SED wavelength range covers filter
  if ( minlam_filt < minlam_SED*z1 || maxlam_filt>maxlam_SED*z1 ) {
    sprintf(c1err,"Invalid obs-frame SED wave range (%.1f - %.1f), zHEL=%.3f",
	    minlam_SED*z1, maxlam_SED*z1, zHEL );
    sprintf(c2err,"ifilt=%d(%s) wave range is %.1f - %.1f",
	    ifilt, cfilt, minlam_filt, maxlam_filt);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  // - - - - - - -
  // check for intrinsic scatter models in sntools_genSmear.c
  //  (e..g, G10, C11).  Get magSmear at all wavelengths.
  ISTAT_SMEAR = istat_genSmear(); // check for smear model
  if ( ISTAT_SMEAR ) {
    for ( ilamobs=0; ilamobs < NLAMFILT; ilamobs++ ) {
      LAMOBS       = FILTER_SEDMODEL[ifilt].lam[ilamobs] ;
      LAMSED       = LAMOBS/z1;   // rest-frame wavelength 
      lam[ilamobs] = LAMSED ;
      magSmear[ilamobs] = 0.0 ;
    }
    get_genSmear( Trest, NLAMFILT, lam, magSmear) ;
  }

  // - - - - - -

  LAMSED_STEP = lamstep_filt ; 

  for ( ilamobs=0; ilamobs < NLAMFILT; ilamobs++ ) {

    TRANS  = FILTER_SEDMODEL[ifilt].transSN[ilamobs] ;

    if ( TRANS < 1.0E-12 && OPT_SPEC==0) 
      { continue ; } // Jul 2013 - skip zeros for leakage

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

    // loop over rest-frame lambda (for SPECTROGRAPH)
    for(LAMSED = LAMSED_MIN; LAMSED <= LAMSED_MAX; LAMSED+=LAMSED_STEP ) {

      // find rest-frame bin for BYOSED ... note that non-uniform
      // bins are allowed, but non-uniform might lead to trouble elsewhere.
      ilamsed = quickBinSearch(NLAM, LAMSED, LAM, fnam);
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
    if ( OPT_SPEC ) { Fspec[ilamobs] = (Finteg_spec * MODELNORM_Fspec) ; }

  } // end ilamobs
  

  // - - - - - - - 
  // store integrated flux in passband 
  *Finteg = (Finteg_filter * x0 * MODELNORM_Finteg);

  return;

} // end INTEG_zSED_BYOSED


// ====================================================
void genSpec_BYOSED(double Tobs, double zHEL, double MU, double MWEBV,
		    double RV_host, double AV_host, // (I)		     
		    double *GENFLUX_LIST,           // (O) fluxGen per bin 
		    double *GENMAG_LIST ) {         // (O) magGen per bin

  // March 2019
  // Return BYOSED spectrum in SPECTROGRAPH bins
  //

  double hc8   = (double)hc ;
  double z1    = 1.0 + zHEL ;
  double x0    = pow(TEN,-0.4*MU);
  double Trest = Tobs/z1 ;
  int NBLAM    = SPECTROGRAPH_SEDMODEL.NBLAM_TOT ;
  int NEWEVT_FLAG = 0;

  int ilam ;
  double Finteg_ignore, FTMP, MAG, ZP, LAM ;
  char fnam[] = "genSpec_BYOSED" ;

  // --------- BEGIN ------------

  // init entire spectum to zero.
  for(ilam=0; ilam < NBLAM; ilam++ ) { GENFLUX_LIST[ilam] = 0.0 ; }

  INTEG_zSED_BYOSED(1, JFILT_SPECTROGRAPH, Tobs, zHEL, x0, 
		    RV_host, AV_host, 
		    Event_BYOSED.NLAM,
		    Event_BYOSED.LAM, 
		    Event_BYOSED.SED, 
		    &Finteg_ignore, GENFLUX_LIST ); // <= returned 


  // convert generated fluxes into mags
  for(ilam=0; ilam < NBLAM; ilam++ ) { 
    LAM  = SPECTROGRAPH_SEDMODEL.LAMAVG_LIST[ilam] ;
    ZP   = SPECTROGRAPH_SEDMODEL.ZP_LIST[ilam] ;
    FTMP = (LAM/(hc8*z1)) * GENFLUX_LIST[ilam] ;
    if ( ZP > 0.0 && FTMP > 0.0 )   { 
      MAG = -2.5*log10(FTMP) + ZP ;
    }
    else  { 
      MAG = MAG_UNDEFINED ;  // model undefined
    }
    GENMAG_LIST[ilam] = MAG ;
  } // end ilam loop over SPECTROGRAPH bins


  return ;

} // end genSpec_BYOSED


// =====================================
void  read_SALT2_template0(void) {

  // read M0 surface of SALT2 model, corresponding to x1=c=0.

  char sedcomment[80] ;
  char SALT2_tempate0_file[MXPATHLEN];

  double Trange[2] = { -40.0,  100.0   } ;
  double Lrange[2] = { 2000.0, 9200.0  } ;
  char   MODELNAME_SALT2[] = "SALT2.JLA-B14" ;
  char fnam[] = "read_SALT2_template0" ;

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
