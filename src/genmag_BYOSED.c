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

 *****************************************/


#include  <stdio.h> 
#include  <math.h>     
#include  <stdlib.h>   
#include  <sys/stat.h>


#include  "sntools.h"           // SNANA community tools
#include  "genmag_SEDtools.h"
#include  "genmag_SIMSED.h"
#include  "sntools_spectrograph.h"
#include  "MWgaldust.h"
#include  "genmag_BYOSED.h"

#ifdef  USE_PYTHON
#include  <Python.h>
//#include <numpy/arrayobject.h>
PyObject *geninit_BYOSED;
#endif

// =========================================================
void init_genmag_BYOSED(char *PATH_VERSION) {

  // Read input directory file(s) for parameters characterizing 
  // how to build your SED.
  //  
  // Inputs:
  //  PATH_VERSION : points to model-param directory


#ifdef USE_PYTHON
  PyObject *genmod, *genclass, *pargs;
#endif

  int  MEMD   = sizeof(double);
  int  MEMC   = sizeof(char);
  char fnam[] = "init_genmag_BYOSED" ;

  // -------------- BEGIN ------------------

  sprintf(BANNER, "%s", fnam);
  print_banner(BANNER);
  printf("   path = '%s' \n",  PATH_VERSION);

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

  pargs  = Py_BuildValue("(s)",PATH_VERSION);
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
		   double MWEBV, double RV_host, double AV_host,
		   int IFILT_OBS, int NOBS, double *TOBS_list, 
		   double *MAGOBS_list, double *MAGERR_list ) {

  // Created Sep 2018
  //
  // Inputs:
  //   EXTERNAL_ID : read next SED when this changes
  //   zHEL        : helio redshift
  //   MU          : distance modulus
  //   MWEBV       : E(B-V) for Milky Wat
  //   IFILT_OBS   : absolute filter index
  //   NOBS        : number of observations
  //   TOBS_list   : list of MJD-PEAKMJD
  //
  // Outputs"
  //   MAGOBS_list   : list of true mags
  //   MAGERR_list   : list of mag errors (place-holder, in case)
  //

  double FLUXSUM_MIN = 1.0E-30 ;
  double z1    = 1.0 + zHEL ;
  double *LAM  = Event_BYOSED.LAM;
  double *SED  = Event_BYOSED.SED ;

  int  ifilt  = IFILTMAP_SEDMODEL[IFILT_OBS] ; // sparse filter index
  char *cfilt = FILTER_SEDMODEL[ifilt].name ;
  double  ZP  = FILTER_SEDMODEL[ifilt].ZP ;    // ZP for flux->mag
  double x0   = pow(10.0,-0.4*MU);             // dimming from dist. mod.

  int    NLAM, o ;
  double Tobs, Trest, FLUXSUM_OBS, FspecDUM[2], magobs ; 

  char fnam[] = "genmag_BYOSED" ;

   #ifdef USE_PYTHON
  // python declarations here
  
   #endif

  // ------------ BEGIN -----------

  for(o=0; o < NOBS; o++ ) { MAGOBS_list[o] = 99.0 ; }



  /* xxx
  printf(" xxx ------------------------------------ \n" ) ;
  printf(" xxx %s: process z=%.3f MU=%.3f RV=%3.1f IFILT_OBS=%d(%s) \n",
	 fnam, zHEL, MU, MWXT_SEDMODEL.RV, IFILT_OBS, cfilt );
  printf(" xxx RVMW=%3.1f  MWEBV=%.3f      RV_host=%3.1f AV=%4.2f \n", 
	 MWXT_SEDMODEL.RV, MWEBV,   RV_host, AV_host );
  xxxx */

  // make sure filter-lambda range is valid
  // xxxx  checkLamRange_SEDMODEL(ifilt,zHEL,fnam);

  // store table info for Galactic & host extinction
  fill_TABLE_MWXT_SEDMODEL(MWXT_SEDMODEL.RV, MWEBV);
  fill_TABLE_HOSTXT_SEDMODEL(RV_host, AV_host, zHEL);   // July 2016

    
  for(o=0; o < NOBS; o++ ) {
    Tobs  = TOBS_list[o];
    Trest = Tobs/z1;
    fetchSED_BYOSED(Trest, MXLAM_BYOSED,   &NLAM, LAM, SED );  
    Event_BYOSED.NLAM = NLAM ;

    // integrate redshifted SED to get observer-frame flux in IFILT_OBS band.
    // FLUXSUM_OBS is returned (ignore FspecDUM)
    INTEG_zSED_BYOSED(0, IFILT_OBS, zHEL, x0, RV_host, AV_host, NLAM, LAM, SED, 
		      &FLUXSUM_OBS, FspecDUM ); // <= returned 
		      
    
    // convert calibrated flux into true magnitude
    if ( FLUXSUM_OBS > FLUXSUM_MIN ) 
      { magobs = ZP - 2.5*log10(FLUXSUM_OBS); }
    else
      { magobs = MAG_ZEROFLUX ; }

    MAGOBS_list[o] = magobs;  // load output array
    MAGERR_list[o] = 0.01;    // not used
  }

  // store SED parameters so that sim can pass to data files
  if ( EXTERNAL_ID != Event_BYOSED.LAST_EXTERNAL_ID )
    { fetchParVal_BYOSED(Event_BYOSED.PARVAL); }

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

  int NPAR=0;
  char fnam[] = "fetchParNames_BYOSED" ;

  // David: need your python magic to return these string names.

  sprintf(parNameList[NPAR],"DUMMY_PAR0" ); NPAR++;
  sprintf(parNameList[NPAR],"DUMMY_PAR1" ); NPAR++;
  sprintf(parNameList[NPAR],"DUMMY_PAR2" ); NPAR++;

  return(NPAR) ;

} // fetchParNames_BYOSED


void fetchParVal_BYOSED(double *parVal) {

  // return list of parameters to calling function (sim)
  // so that these parameters can be included in the
  // data files.
  //
  // Called once per event.
  
  int NPAR=0;
  char fnam[] = "fetchParVal_BYOSED" ;

  // ------------- BEGIN ------------------

  // David: need python function to return these values.

  parVal[0] = 1.111 ; 
  parVal[1] = 2.222 ; 
  parVal[2] = 3.333 ; 

  return ;

} // end fetchParVal_BYOSED

// =================================================
void fetchSED_BYOSED(double Trest, int MXLAM,
		     int *NLAM_SED, double *LAM_SED, double *FLUX_SED) {

  // return rest-frame SED to calling function; 
  // Inputs:
  //   Trest : rest frame epochs (Trest=0 at peak)
  //   MXLAM : abort iof *NLAM > MXLAM
  //
  // Output
  //  *NLAM_SED  : number of wavelenth bins for SED
  //  *LAM_SED   : array of wavelengths for SED
  //  *FLUX_SED  : SED flux in each wave bin

  char fnam[] = "fetchSED_BYOSED" ;

  // ------------ BEGIN -----------

  *NLAM_SED = 0 ; // init output
  
#ifdef USE_PYTHON
  PyObject *pmeth, *pargs, *pValue, *NLAM;
  PyArrayObject *arrReturn;
  // python declarations here
  pmeth  = PyObject_GetAttrString(geninit_BYOSED, "fetchSED_BYOSED");
  pargs  = Py_BuildValue("ddi",z,Tobs,MXLAM);
  pValue   = PyEval_CallObject(pmeth, pargs);

  arrReturn = (PyArrayObject *)(pValue);
  NLAM = PyArray_GETPTR1(arrReturn,0,0)
    
  for(ilam=0; ilam < NLAM; ilam++ ) {
    // interpolate flux to Trest
    LAM_SED[ilam]  = PyArray_GETPTR2(arrReturn,1,ilam) ;
    FLUX_SED[ilam] = PyArray_GETPTR2(arrReturn,2,ilam);
  }

  Py_DECREF(pmeth);
  Py_DECREF(pargs);
  Py_DECREF(pValue);
  //Py_DECREF(pNLAM);
  //Py_DECREF(pFLUX);
  //Py_DECREF(pLAM);
  
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
void INTEG_zSED_BYOSED(int OPT_SPEC, int ifilt_obs, double zHEL, double x0, 
		       double RV_host, double AV_host,
		       int NLAM, double *LAM, double *SED,
		       double *Finteg, double *Fspec ) {

  // Created Dec 2018 by R.K.
  // Return integrated obs-frame flux in filter passband 
  // Be careful with rest-frame and obs-frame
  //s 
  // Inputs:
  //   OPT_SPEC     1 -> return obs-frame spectrum withing filter-band.
  //   ifilt_obs    filter index
  //   zHEL         helio redshift
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
  // !!! Dec 12 2018: Finteg is ready, but Fspec is NOT ready !!!
  //

  int    ifilt          = IFILTMAP_SEDMODEL[ifilt_obs] ;
  int    NLAMFILT       = FILTER_SEDMODEL[ifilt].NLAM ;
  int   DO_SPECTROGRAPH = ( ifilt_obs == JFILT_SPECTROGRAPH ) ;
  double meanlam_obs    = FILTER_SEDMODEL[ifilt].mean ;  // mean lambda
  double z1             = 1.0 + zHEL ;
  double meanlam_rest   = meanlam_obs/z1 ;
  double minlam_filt    = FILTER_SEDMODEL[ifilt].minlam ;
  double maxlam_filt    = FILTER_SEDMODEL[ifilt].maxlam ;
  double lamstep_filt   = FILTER_SEDMODEL[ifilt].lamstep ;
  double minlam_SED     = LAM[0];
  double maxlam_SED     = LAM[NLAM-1];
  double hc8            = (double)hc ;
  double MODELNORM      = lamstep_filt * 1.0E0 / hc8 ;

  int    LABORT, ilamobs, ilamsed;
  double TRANS, MWXT_FRAC, HOSTXT_FRAC, FLUXSUM = 0.0 ;
  double LAMOBS, LAMSED, LAMSED_MIN, LAMSED_MAX, LAMSED_STEP ;
  double TMPLAM[3], TMPSED[3];
  double Fbin_forFlux, Fbin_forSpec, FTMP; 
  double Finteg_filter=0.0, Finteg_spec=0.0 ;
  char *cfilt  = FILTER_SEDMODEL[ifilt].name ;

  int  LDMP = 0 ;
  char fnam[]  = "INTEG_zSED_BYOSED" ;

  // ------------- BEGIN -----------

  *Finteg = 0.0 ; // init output flux for filter

  // first make sure that SED wavelength range covers filter
  if ( minlam_filt < minlam_SED*z1 || maxlam_filt>maxlam_SED*z1 ) {
    sprintf(c1err,"Invalid obs-frame SED wave range (%.1f - %.1f), zHEL=%.3f",
	    minlam_SED*z1, maxlam_SED*z1, zHEL );
    sprintf(c2err,"ifilt=%d(%s) wave range is %.1f - %.1f",
	    ifilt, cfilt, minlam_filt, maxlam_filt);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

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

      Fbin_forFlux = (FTMP * HOSTXT_FRAC*MWXT_FRAC * LAMSED*TRANS);
      Fbin_forSpec = (FTMP * HOSTXT_FRAC*MWXT_FRAC );
      Finteg_filter  +=  Fbin_forFlux ;

    } // end loop over LAMSED

  } // end ilamobs
  

  *Finteg = Finteg_filter * x0 * MODELNORM ;

  return;

} // end INTEG_zSED_BYOSED


// ====================================================
void genSpec_BYOSED(double TOBS, double z, double MWEBV,
		    double RV_host, double AV_host, // (I)		     
		    double *FLUXOBS_LIST,           // (O) fluxGen per bin 
		    double *GENMAG_LIST ) {         // (O) magGen per bin

  // Dec 2018
  // Return BYOSED spectrum in SPECTROGRAPH bins
  //
  // !!! NOT READY !!!
  // --------- BEGIN ------------


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
