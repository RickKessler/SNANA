/**************************************************************
 Generate SN obs-frame mags using SALT2 model.
                                                             
 R. Kessler  Apr 2009 : re-written for sim & fitter; includes errors


 There are three init functions that must be called
 exeternally in the following order:

 1.  init_primary_SALT2(name, NLAM, spec);
     passes the primary spectrum (Vega, BD17, 1/lam^2 ...)

 2.  init_filter_SALT2(ifilt_obs..) is called for each
     observer-frame filter; user passes transmission vs.
     wavelength along with a global lambda-shift.
     Remember to pass Bessell-B band to compute mB at
     the end of each fit.

 3.  init_genmag_SALT2() reads SED templates and computes
     needed info from filters.


             HISTORY
            ~~~~~~~~~~~


 Oct 2020: minor refactor for INTGEG_zSED_SALT2 and SALT2magerr;  
           needed to handle SALT3 or SALT2.

 Dec 28 2020: use function setFlags_ISMODEL_SALT2 to set ISMODEL_SALT2[3]

 Mar 23 2021: use get_LAMTRANS_SEDMODEL util to getch LAM and TRANS

 Apr 27 2021: minor refactor so that default SALT3 models are checked
              in $SNDATA_ROOT/models/SALT3 (not under SALT2)

 May 31 2021: refactor to receive parList_SN and parList_HOST so that
              logMass can be included and passed to get_genSmear.

 Aug 26 2021: remove buggy 1/z1 factor in genSpec_SALT2

 Oct 01 2021: no longer set magerr=5.0 -> avoid LC fit discontinuity.

*************************************/

#include "sntools.h"           // community tools
#include "sntools_genSmear.h"
#include "sntools_spectrograph.h"
#include "sntools_devel.h"
#include "genmag_SEDtools.h"
#include "genmag_extrap.h"
#include "genmag_SALT2.h" 
#include "MWgaldust.h"

// =======================================================
// define mangled functions with underscore (for fortran)

int init_genmag_salt2__(char *model_version, char *model_extrap, 
			int *OPTMASK) {
  int istat;
  istat = init_genmag_SALT2 ( model_version, model_extrap,  *OPTMASK);
  return istat ;
} 


void genmag_salt2__(int *OPTMASK, int *ifilt, 
		    double *parList_SN, double *parList_HOST, double *mwebv,
		    double *z, double *z_forErr, int *nobs, double *Tobs_list, 
		    double *magobs_list, double *magerr_list ) {

  genmag_SALT2(*OPTMASK, *ifilt, parList_SN, parList_HOST, *mwebv,
	       *z, *z_forErr, *nobs, Tobs_list, 
	       magobs_list, magerr_list );
}


double salt2x0calc_(double *alpha, double *beta, double *x1,   
		    double *c, double *dlmag ) {
  double x0;
  x0 = SALT2x0calc(*alpha, *beta, *x1,  *c, *dlmag );
  return x0;
} //

double salt2mbcalc_(double *x0) {
  double mB;
  mB = SALT2mBcalc(*x0);
  return mB ;
} //



int gencovar_salt2__ (
		  int *MATSIZE         // (I) row-len (or col-len)
                  ,int *ifilt_obs      // (I) list of 'matsize' filter indices
                  ,double *epobs       // (I) list of 'matsize' rest days
		  ,double *z            // (I) redshift
		  ,double *parList_SN 
		  ,double *parList_HOST
		  ,double *mwebv
                  ,double *covar       // (O) covariance matrix
                  ) {
  int istat ;
  istat = gencovar_SALT2 ( *MATSIZE, ifilt_obs, epobs, *z, parList_SN,
			   parList_HOST, *mwebv, covar ) ;
  return istat;
}

// external spline function

extern void in2dex_(int *ispline, int *N2D,
		    double *XX, double *YY, double *ZZ,
		    double *XLIM, double *YLIM, double *SS, int *IERR );

extern double ge2dex_ ( int *IND, double *Trest, double *Lrest, int *IERR ) ;

/****************************************************************
  init_genmag_SALT2:
    o reads in the filters from FilterFiles
    o calculates filter mean and AB zeropoint
    o reads the templates from TemplateFiles
    o returns true if successful

   Note: must call init_filter_SEDMODEL() and init_primary_SEDMODEL()
          before calling this function.

  Feb 24, 2009 RSK: read colorCorrection from file

  Oct 13, 2009: if SED.DAT exists, NSURFACE=1; else NSURFACE=2 (nominal)

  Apr 24, 2010:  call  check_sedflux_bins(...) to ensure idential 
                 SED binning for each surface.

  Jun 27, 2010: add OPTMASK argument. Bit 8 (OPTMASK=128) sets 
                NLAMPOW_SEDMODEL=0 to implement legacy style 
                where color*XTMW is outside the integrals.


  OPTMASK bit 8 (128) : legacy option sets NLAMPOW_SEDMODEL=0

  Jul 27, 2010:  new call to load_mBoff_SALT2();

  Aug 04, 2010: call colordump_SALT2 for more wavelengths ...
                enough to span 2000 to 10,000 A.

  Jan 18, 2011: read model from getenv(PRIVATE_MODELPATH_NAME) if this
                env variable exists.

  Mar 3, 2011: move errmap-reading  into read_SALT2errmaps().

  Aug 9 2017: Lrange[1] -> 30000 (was 20000)  
 
  Mar 18 2017: add genmag_SALT2 argument z_forErr; allows passing
               fixed redshift for photo-z fits, analogous to x1_forErr.

  Jun 24 2018: add argument MODEL_EXTRAP to override
                GENMODEL_EXTRAP_LATETIME in SALT2.INFO file

  Aug 02 2019:
    set SALT2_PREFIX_FILENAME to either salt2 (default) or salt3.
    File name prefixes thus correspond to model name.

 Aug 26 2019: implement RELAX_IDIOT_CHECK_SALT2 for P18 to avoid abort.

 Nov 7 2019: for SALT3, remove x1*M1/M0 term in error; see ISMODEL_SALT3.
 Jan 19 2020: in INTEG_zSED_SALT2, fix memory leak related to local magSmear.
 Mar 24 2020: if ISMODEL_SALT3, then read SALT3.INFO 
 Sep 03 2020: check REQUIRE_DOCANA
 Oct 16 2020: refactor SALT2magerr to handle SALT3 vs. SALT2 vartot

 Nov 23 2020: pass and store *SURVEY; used to match
              MAGSHIFT and WAVESHIFT keys in SALT2.INFO file.

 Apr 27 2021: minor refactor to set default SALT2 or SALT3 model location.

 Feb 22 2022: set DEBUG_SALT2 with OPTMASK += 1024
 Dec    2023: integrate 2nd component.

****************************************************************/

int init_genmag_SALT2(char *MODEL_VERSION, char *MODEL_EXTRAP_LATETIME,
		      int OPTMASK) {

  double Trange[2], Lrange[2];
  int  ised, nflux_nan;
  int  retval = 0   ;
  int  ABORT_on_LAMRANGE_ERROR = 0;
  int  ABORT_on_BADVALUE_ERROR = 1;
  char BANNER[120], tmpFile[200], sedcomment[40], version[60]  ;
  char fnam[] = "init_genmag_SALT2" ;

  // -------------- BEGIN --------------

  // extrac OPTMASK options
  sprintf(BANNER, "%s : Initialize %s", fnam, MODEL_VERSION );
  print_banner(BANNER);
 

  if ( NFILT_SEDMODEL == 0 ) {
    sprintf(c1err,"No filters defined ?!?!?!? " );
    sprintf(c2err,"Need to call init_filter_SEDMODEL");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  ABORT_on_LAMRANGE_ERROR = ( OPTMASK & GENMODEL_MSKOPT_SALT2_ABORT_LAMRANGE ) ; 

  init_NEGFLAM_SEDMODEL(OPTMASK);


  DEBUG_SALT2 = ( OPTMASK & GENMODEL_MSKOPT_SALT2_DEBUG );
  //  DEBUG_SALT2 = 1; // xxx REMOVE

  // summarize filter info
  filtdump_SEDMODEL();


  // ==========================================
  // construct path to SALT2 surfaces

  extract_MODELNAME(MODEL_VERSION,             // input
		    SALT2_MODELPATH, version); // returned

  // parse version string to check if SALT2 or SALT3
  setFlags_ISMODEL_SALT2(version); // set ISMODEL_SALT2 and ISMODEL_SALT3
  sprintf(SALT2_PREFIX_FILENAME,"salt%d", IMODEL_SALT);
  
  if ( getenv(PRIVATE_MODELPATH_NAME) != NULL ) {
    sprintf( SALT2_MODELPATH, "%s/%s", 
	     getenv(PRIVATE_MODELPATH_NAME), version );    
  }
  else if ( strlen(SALT2_MODELPATH) > 0 ) {
    // do nothing;
  }
  else {
    // default location under $SNDATA_ROOT
    sprintf( SALT2_MODELPATH, "%s/models/SALT%d/%s",  
	     getenv("SNDATA_ROOT"), IMODEL_SALT, version );
  }

  // Mar 24 2020 INFO file depends on SALT2 or SALT3
  if ( ISMODEL_SALT2 ) 
    { sprintf(SALT2_INFO_FILE,  "SALT2.INFO" ); }
  else if ( ISMODEL_SALT3 ) 
    { sprintf(SALT2_INFO_FILE,  "SALT3.INFO" ); }
  else {
    sprintf(c1err,"Unknown model; expecting SALT2 or SALT3");
    sprintf(c2err,"Check SALT2.INFO/SALT3.INFO") ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }


  // set defaults for two surfaces (nominal SALT2)
  SEDMODEL.NSURFACE   = NSURFACE_SALT2();
  SEDMODEL.FLUXSCALE  = X0SCALE_SALT2; 
  SEDMODEL.MAGERR_FIX = -9.0 ;        // => use calculated errors


  // May 2013:
  // if re-reading the same SALT2 version, then skip re-reading the files.
  // WARNING: SALT2_VERSION is not set on first pass, so it may give
  //          valgrind errors.
  int SKIPREAD = 0 ;
  if ( strcmp(SALT2_VERSION,version) == 0 ) { 
    printf("\t Re-init %s -> skip reading files. \n", version);
    fflush(stdout);
    init_SALT2interp_SEDFLUX();
    init_SALT2interp_ERRMAP();
    SKIPREAD = 1;  // set logical in case we need it later
    return retval ;
  }

  sprintf(SALT2_VERSION,"%s", version); // May 15 2013

  // ============================

  read_SALT2_INFO_FILE(OPTMASK);  

  RELAX_IDIOT_CHECK_SALT2 = (INPUT_SALT2_INFO.RESTLAMMAX_FILTERCEN > 12000.0);

  // ============================
  // set extreme ranges to read anything
  Trange[0] = -20. ;
  Trange[1] = 200. ;
  Lrange[0] = (double)LAMMIN_SEDMODEL ;
  Lrange[1] = (double)LAMMAX_SEDMODEL ;

  SEDMODEL_MWEBV_LAST     = -999.   ;
  SEDMODEL_HOSTXT_LAST.AV = -999.   ;
  SEDMODEL_HOSTXT_LAST.z  = -999.   ;

  SPECTROGRAPH_SEDMODEL.NBLAM_TOT = 0 ; // spectrograph option

  malloc_METADATA_SEDMODEL(SEDMODEL.NSURFACE, 0); 
  malloc_SEDFLUX_SEDMODEL(&TEMP_SEDMODEL,0,0,0);

  // ------- Now read the spectral templates -----------

  for ( ised = 0 ; ised < SEDMODEL.NSURFACE ; ised++ ) {

    sprintf(tmpFile, "%s/%s_template_%d.dat", 
	    SALT2_MODELPATH, SALT2_PREFIX_FILENAME, ised );

    //  printf("  Read Template file: \n\t %s \n", tmpFile);

    sprintf(sedcomment,"SALT%d-%d", IMODEL_SALT, ised);

    rd_sedFlux(tmpFile, sedcomment, Trange, Lrange
	       ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL, 0
	       ,&TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY, &TEMP_SEDMODEL.DAYSTEP
	       ,&TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM, &TEMP_SEDMODEL.LAMSTEP
	       ,TEMP_SEDMODEL.FLUX,  TEMP_SEDMODEL.FLUXERR
	       ,&nflux_nan );
    
    // July 18 2018: check for UV extrap to avoid filter dropouts
    double UVLAM = INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX;
    if ( UVLAM > 0.0 ) { UVLAM_EXTRAPFLUX_SEDMODEL(UVLAM,&TEMP_SEDMODEL); }

    // make sure that DAY and LAM binning is identical for each surface
    check_sedflux_bins(ised, "DAY", 
       TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY[0], TEMP_SEDMODEL.DAYSTEP);
    check_sedflux_bins(ised, "LAM", 
       TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM[0], TEMP_SEDMODEL.LAMSTEP);

    // transfer TEMP_SEDMODEL to permanent storage
    fill_SALT2_TABLE_SED(ised);

  } //  end loop over SED templates



  load_mBoff_SALT2();

  // ========== read error maps with same format as SED flux
  read_SALT2errmaps(Trange,Lrange);  

  // ------------------ 
  // Read color-dispersion vs. wavelength
  read_SALT2colorDisp();
  
  // abort if any ERRMAP has invalid wavelength range (Sep 2019)
  if ( ABORT_on_LAMRANGE_ERROR ) { check_lamRange_SALT2errmap(-1); }
  if ( ABORT_on_BADVALUE_ERROR ) { check_BADVAL_SALT2errmap(-1); }
  

  // fill/calculate color-law table vs. color and rest-lambda
  fill_SALT2_TABLE_COLORLAW();

  // init interp (for splines only)
  init_SALT2interp_SEDFLUX();
  init_SALT2interp_ERRMAP();

  NCALL_DBUG_SALT2 = 0;

  // Summarize CL and errors vs. lambda for Trest = x1 = 0.
  errorSummary_SALT2();
  
  init_extrap_latetime_Ia(MODEL_EXTRAP_LATETIME);

  init_calib_shift_SALT2train(); // Nov 2020

  // determine method to extrapolate in phase (RK, Sep 2023)
  int EXTRAP_METHOD_PREFER = EXTRAP_PHASE_FLAM; 
  set_METHOD_EXTRAP_PHASE(EXTRAP_METHOD_PREFER);


  fflush(stdout) ;


  //  test_SALT2colorlaw1();

  // ===========================================


  printf("\n  %s : Done. \n", fnam );

  //  debugexit("SALT2 init");
  return retval;

} // end of function init_genmag_SALT2

// ***********************************************
void setFlags_ISMODEL_SALT2(char *version) {

  // Created Dec 28 2020
  // Based on input *version, set global flags
  // ISMODEL_SALT2 and ISMODEL_SALT3
  //
  // Models are of the form
  //   [path]/SALT2.XYZ  or
  //   [path]/SALT3.XYZ
  //
  // So check 5 characters before the dot. 
  //
  // Mar 16 2021: abort on null version
  // Apr 27 2021: set IMODEL_SALT = 2 or 3
  // Sep 20 2023: extend version_near_dot array from 60 to MXPATHLEN [fix crash]

  int  index_dot, set=0 ;
  int  LENSALT2 = strlen("SALT2");
  char *dot, version_near_dot[MXPATHLEN];
  char fnam[] = "setFlags_ISMODEL_SALT2" ;

  // ------------- BEGIN ------------

  ISMODEL_SALT2 = false;
  ISMODEL_SALT3 = false;
  IMODEL_SALT   = -9 ;

  if ( strlen(version) < LENSALT2+1 ) {
    sprintf(c1err,"Invalid version = '%s' (name too short)", version);
    sprintf(c2err,"version must contain SALT2.[bla] or SALT2.[bla]");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  dot       = strchr(version, '.');
  index_dot = (int)(dot - version);
  sprintf(version_near_dot, "%s", &version[index_dot-LENSALT2] );

  if ( strstr(version_near_dot,"SALT2") != NULL ) 
    { set=1; ISMODEL_SALT2 = true;  IMODEL_SALT = 2 ; } 
  if ( strstr(version_near_dot,"SALT3") != NULL ) 
    { set=1; ISMODEL_SALT3 = true;  IMODEL_SALT = 3 ; } 

  if ( !set ) {
    printf("\n\t index_dot = %d \n", index_dot);
    sprintf(c1err,"Unable to set ISMODEL_SALT2 or ISMODEL_SALT3");
    sprintf(c2err,"Check GENNODEL: %s (%s)", version, version_near_dot);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  printf("\t ISMODEL = SALT%d\n", IMODEL_SALT);

  fflush(stdout);

  return;

} // end setFlags_ISMODEL_SALT2

// ==============================
int NSURFACE_SALT2(void) {
  
  // Created Dec 2023
  // Return number of SALT2 surfaces. Original model has 2,
  // but there can be 3 for x2 * M2.

  char template_file[200];
  int  gzipFlag, i, N = 0 ;
  FILE *fp;
  char fnam[] = "NSURFACE_SALT2";

  // ------------ BEGIN -------------

  if ( DEBUG_SALT2 == 0 ) { return 2; } // legacy

  for(i=0; i < 3; i++ ) {
    sprintf(template_file, "%s/%s_template_%d.dat", 
	    SALT2_MODELPATH, SALT2_PREFIX_FILENAME, i);
    fp = open_TEXTgz(template_file, "rt", &gzipFlag );
    if ( fp != NULL ) { N++; fclose(fp); }
  }

  // xxx printf(" xxx %s: return NSURFACE = %d \n", fnam, N);
  
  return N;

} // NSURFACE_SALT2

// ***********************************************
void fill_SALT2_TABLE_SED(int ISED) {

  // transfer original TEMP_SEDMODEL contents to permanent
  // (dynamically allocated) SALT2_TABLE.SEDFLUX.
  // Note that SEDMODEL.DAY is not allocated or filled
  // since the DAY-bins here must be uniform.
  //
  // If spline-interp option is set, then apply spline-interp
  // here and store finer-binned SEDs so that faster linear-
  // interpolation can be used inside the integration loops.
  //
  // Jun 9, 2011: load SEDMODEL.LAMMIN[ISED] and SEDMODEL.LAMMAX[ISED]
  // Dec 30, 2013: add PRE-ABORT dump when FRATIO is too large.
  // Jul 30, 2016: call check_uniform_bins( ... )
  // May 05 2022: add abort protection when  F_orig=0.0 
  //          (for SALT2 surfaces computed with Jacobian method)
  //
#define N1DBIN_SPLINE 3

  int 
    jflux_orig, IDAY, ILAM, IDAY_ORIG, ILAM_ORIG,iday, ilam
    ,NDAY_ORIG,  NLAM_ORIG, NDAY_TABLE, NLAM_TABLE
    ,NREBIN_DAY, NREBIN_LAM, INTERP_OPT, EDGE, I8, I8p
    ;

  double 
    xi, DAY, LAM, DIF, FRAC, DAYSTEP_ORIG, LAMSTEP_ORIG
    ,F2D_orig[N1DBIN_SPLINE][N1DBIN_SPLINE]
    ,F_interp, F_orig, FDIF, FSUM, FRATIO
    ,FDAY[N1DBIN_SPLINE],FRATIO_CHECK
    ,*ptrLAM, *ptrDAY
    ;

  char 
     cmsg1[60], cmsg2[60]
    ,fnam[]   = "fill_SALT2_TABLE_SED" 
    ,tagLAM[] = "interp-LAM"
    ,tagDAY[] = "interp-DAY"
    ;

  // ---------------- BEGIN --------------

  NDAY_ORIG  = TEMP_SEDMODEL.NDAY ;
  NLAM_ORIG  = TEMP_SEDMODEL.NLAM ;
  INTERP_OPT = INPUT_SALT2_INFO.SEDFLUX_INTERP_OPT ;
  
  // ----------
  // check for uniform binning (July 2016)
  check_uniform_bins(NDAY_ORIG, TEMP_SEDMODEL.DAY, "DayGrid(SALT2)");

  // ----------
  if ( INTERP_OPT == SALT2_INTERP_SPLINE ) {
    NREBIN_DAY = INPUT_SALT2_INFO.INTERP_SEDREBIN_DAY ; 
    NREBIN_LAM = INPUT_SALT2_INFO.INTERP_SEDREBIN_LAM ; 
  }
  else {
    NREBIN_DAY = NREBIN_LAM = 1;
  }

  NDAY_TABLE = NDAY_ORIG * NREBIN_DAY ; 
  NLAM_TABLE = NLAM_ORIG * NREBIN_LAM ; 

  SEDMODEL.LAMSTEP[ISED] = TEMP_SEDMODEL.LAMSTEP ; // Jul 30 2016

  I8  = sizeof(double);
  I8p = sizeof(double*);

  if ( ISED == 0 ) {

    // load SEDMODEL struct for IFILTSTAT function
    SEDMODEL.LAMMIN_ALL = TEMP_SEDMODEL.LAM[0] ;
    SEDMODEL.LAMMAX_ALL = TEMP_SEDMODEL.LAM[NLAM_ORIG-1] ;

    SALT2_TABLE.NDAY    = NDAY_TABLE ;
    SALT2_TABLE.DAYSTEP = TEMP_SEDMODEL.DAYSTEP/(double)(NREBIN_DAY) ;
    SALT2_TABLE.DAYMIN  = TEMP_SEDMODEL.DAY[0] ;
    SALT2_TABLE.DAYMAX  = TEMP_SEDMODEL.DAY[NDAY_ORIG-1] ;
    SALT2_TABLE.DAY     = (double*)malloc(I8*NDAY_TABLE);

    SALT2_TABLE.NLAMSED = NLAM_TABLE ;
    SALT2_TABLE.LAMSTEP = TEMP_SEDMODEL.LAMSTEP/(double)(NREBIN_LAM) ;
    SALT2_TABLE.LAMMIN  = TEMP_SEDMODEL.LAM[0] ;
    SALT2_TABLE.LAMMAX  = TEMP_SEDMODEL.LAM[NLAM_ORIG-1] ;
    SALT2_TABLE.LAMSED  = (double*)malloc(I8*NLAM_TABLE);

    for ( IDAY=0; IDAY < NDAY_TABLE; IDAY++ ) {  
      xi = (double)IDAY ;    
      SALT2_TABLE.DAY[IDAY] = 
	SALT2_TABLE.DAYMIN + (xi * SALT2_TABLE.DAYSTEP);
    }

    for ( ILAM=0; ILAM < NLAM_TABLE; ILAM++ ) {  
      xi = (double)ILAM ;    
      SALT2_TABLE.LAMSED[ILAM] = 
	SALT2_TABLE.LAMMIN + (xi * SALT2_TABLE.LAMSTEP);
    }
  }


  sprintf(cmsg1,"LAM(MIN,MAX,STEP)=%4.0f,%4.0f,%1.0f",
	  SALT2_TABLE.LAMMIN, SALT2_TABLE.LAMMAX, SALT2_TABLE.LAMSTEP );
  sprintf(cmsg2,"DAY(MIN,MAX,STEP)=%2.0f,%2.0f,%2.1f",
	  SALT2_TABLE.DAYMIN, SALT2_TABLE.DAYMAX, SALT2_TABLE.DAYSTEP );
  printf("  Store SED-%d  %s  %s \n\n", ISED, cmsg1, cmsg2 );


  DAYSTEP_ORIG = TEMP_SEDMODEL.DAYSTEP ;
  LAMSTEP_ORIG = TEMP_SEDMODEL.LAMSTEP ;


  // --------------------------------------
  // allocate memory for SED surface
  SALT2_TABLE.SEDFLUX[ISED] = (double**)malloc(I8p*NDAY_TABLE);
  for ( IDAY=0; IDAY < NDAY_TABLE; IDAY++ ) {
    SALT2_TABLE.SEDFLUX[ISED][IDAY] = (double*)malloc(I8*NLAM_TABLE); 
  }


  // --------------------------------------
  // store SED table.

  for ( IDAY=0; IDAY < NDAY_TABLE; IDAY++ ) {

    // get day-index on original grid
    DAY       = SALT2_TABLE.DAY[IDAY]; // fine-binned DAY
    DIF       = DAY - SALT2_TABLE.DAYMIN + 0.0001 ;
    IDAY_ORIG = (int)(DIF/DAYSTEP_ORIG);

    if ( INTERP_OPT == SALT2_INTERP_SPLINE ) {
      FRAC = (DAY - TEMP_SEDMODEL.DAY[IDAY_ORIG])/DAYSTEP_ORIG;
      if ( FRAC < 0.5 && IDAY_ORIG > 0 ) { IDAY_ORIG-- ; }
      if ( IDAY_ORIG > NDAY_ORIG - N1DBIN_SPLINE ) 
	{ IDAY_ORIG = NDAY_ORIG - N1DBIN_SPLINE ; }
    }

    for ( ILAM=0; ILAM < NLAM_TABLE; ILAM++ ) {

      // for LINEAR option, just take value at node (no interp here)
      if ( INTERP_OPT == SALT2_INTERP_LINEAR ) {
	  jflux_orig  = NLAM_ORIG*IDAY + ILAM ;
	  SALT2_TABLE.SEDFLUX[ISED][IDAY][ILAM] = 
	    TEMP_SEDMODEL.FLUX[jflux_orig] ;
	  continue ; // skip interp stuff below
      }

      // the code below is for the spline option

      // get lam-index on original grid
      LAM       = SALT2_TABLE.LAMSED[ILAM]; // fine-binned lambda
      DIF       = LAM - SALT2_TABLE.LAMMIN + 0.0001 ;
      ILAM_ORIG = (int)(DIF/LAMSTEP_ORIG);

      FRAC  = (LAM - TEMP_SEDMODEL.LAM[ILAM_ORIG])/LAMSTEP_ORIG ;
      if ( FRAC < 0.5  &&  ILAM_ORIG > 0 ) { ILAM_ORIG-- ; }
      if ( ILAM_ORIG > NLAM_ORIG - N1DBIN_SPLINE ) 
	{ ILAM_ORIG = NLAM_ORIG - N1DBIN_SPLINE ; }

      // build temp flux-grid around point to interpolate

      ptrLAM = &TEMP_SEDMODEL.LAM[ILAM_ORIG] ;
      ptrDAY = &TEMP_SEDMODEL.DAY[IDAY_ORIG] ;

      for ( iday=0; iday < N1DBIN_SPLINE; iday++ ) {
	for ( ilam=0; ilam < N1DBIN_SPLINE; ilam++ ) {
	  jflux_orig  = NLAM_ORIG*(IDAY_ORIG+iday) + (ILAM_ORIG+ilam) ;
	  F2D_orig[iday][ilam] = TEMP_SEDMODEL.FLUX[jflux_orig] ;
	}
	// interpolate across lambda to get FDAY
	FDAY[iday] = quadInterp( LAM, ptrLAM, F2D_orig[iday], tagLAM);
      }

      // Now interpolate across DAY
      F_interp = quadInterp( DAY, ptrDAY, FDAY, tagDAY);
      SALT2_TABLE.SEDFLUX[ISED][IDAY][ILAM] = F_interp ;

      // DDDDDDDDDDDDDDDDDDDDDDDDDDDDD
      if ( IDAY == -10 && ILAM < -6 ) {
	printf(" DDDDD ------------------------------------- \n");

	printf(" DDDDD ptrDAY = %5.1f %5.1f %5.1f \n",
	       *(ptrDAY+0),  *(ptrDAY+1), *(ptrDAY+2) );

	printf(" DDDDD FDAY = %le %le %le \n",
	       *(FDAY+0),  *(FDAY+1), *(FDAY+2) );

	printf(" DDDDD DAY[%d]=%6.1f LAM[%d]=%7.1f  F_interp = %le \n",
	       IDAY, DAY, ILAM, LAM, F_interp );
      }
      // DDDDDDDDDDDDDDDDDDDDD


    } // ILAM
  } // IDAY



  // Now an idiot check.
  // Loop over original grid (nodes) and make sure that
  // the finer grid agrees at the nodes.
  // Start at 2nd index for  both DAY and LAM to avoid sharp rise at start

  for ( IDAY_ORIG=1; IDAY_ORIG < NDAY_ORIG; IDAY_ORIG++ ) {
    for ( ILAM_ORIG=1; ILAM_ORIG < NLAM_ORIG; ILAM_ORIG++ ) {

      EDGE = 0;

      if ( IDAY_ORIG <=1 || IDAY_ORIG >= NDAY_ORIG-2 ) { EDGE = 1 ; }
      if ( ILAM_ORIG <=1 || ILAM_ORIG >= NLAM_ORIG-2 ) { EDGE = 1 ; }

      // make looser check at edge-boundary where the interpolation
      // may be a little off.
      if ( EDGE || RELAX_IDIOT_CHECK_SALT2 ) 
	{ FRATIO_CHECK = 1.0E-3 ; }
      else
	{ FRATIO_CHECK = 1.0E-5 ; } // Aug 28 2019: E-6 -> E-5

      IDAY = IDAY_ORIG * NREBIN_DAY ;
      ILAM = ILAM_ORIG * NREBIN_LAM ;

      jflux_orig  = NLAM_ORIG*IDAY_ORIG + ILAM_ORIG ;
      F_orig      = TEMP_SEDMODEL.FLUX[jflux_orig] ;
      F_interp    = SALT2_TABLE.SEDFLUX[ISED][IDAY][ILAM];
      FDIF        = F_interp - F_orig ;
      FSUM        = F_interp + F_orig ;

      if ( RELAX_IDIOT_CHECK_SALT2 && F_orig < 1.0E-25 ) { continue; } 

      if ( FSUM > 0.0 && F_orig != 0.0 ) 
	{ FRATIO = FDIF / FSUM ; }
      else
	{ FRATIO = 0.0 ; }
	    
      if ( fabs(FRATIO) > FRATIO_CHECK ) {
	print_preAbort_banner(fnam);
	printf("  FRATIO = FDIF/FSUM = %f  (FRATIO_CHECK=%le)\n", 
	       FRATIO, FRATIO_CHECK);
	printf("  IDAY=%4d  IDAY_ORIG=%4d  \n", IDAY, IDAY_ORIG);
	printf("  ILAM=%4d  ILAM_ORIG=%4d  \n", ILAM, ILAM_ORIG);
	printf("\n");

	// print 3x3 Flux matrix vs. LAM and DAY 
	int ilam, iday, jflux;	
	
	printf("\t LAM\\DAY");
	for(iday=IDAY_ORIG-1; iday<=IDAY_ORIG+1; iday++ ) 
	  { printf("   %8.1f     ", SALT2_TABLE.DAY[iday*NREBIN_DAY] ); }
	printf("\n");

	for(ilam=ILAM_ORIG-1; ilam <= ILAM_ORIG+1; ilam++ ) {
	  printf("\t %6.1f : ", SALT2_TABLE.LAMSED[ilam*NREBIN_LAM] );
	  for(iday=IDAY_ORIG-1; iday<=IDAY_ORIG+1; iday++ ) {
	    jflux     = NLAM_ORIG*iday + ilam ;
	    printf("%14.6le  ", TEMP_SEDMODEL.FLUX[jflux]);
	  }
	  printf("\n");
	}	

	sprintf(c1err,"Bad SED-%d interp at DAY[%d]=%3.1f  LAM[%d]=%6.1f"
		,ISED
		,IDAY_ORIG, SALT2_TABLE.DAY[IDAY]
		,ILAM_ORIG, SALT2_TABLE.LAMSED[ILAM] );
	sprintf(c2err,"F[interp/orig] = %le / %le = %f",
		F_interp, F_orig, F_interp/F_orig );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
    }
  }

  
  /*
  IDAY = 12;  ILAM = 40; 
  printf("\t  xxxx DAY = SALT2_TABLE.DAY[%d] = %f \n", 
	 IDAY, SALT2_TABLE.DAY[IDAY] );
  printf("\t  xxxx LAM = SALT2_TABLE.LAMSED[%d] = %f \n", 
	 ILAM, SALT2_TABLE.LAMSED[ILAM] );
  printf("\t  xxxx SEDFLUX[%d] = %le \n", 
	 ISED, SALT2_TABLE.SEDFLUX[ISED][IDAY][ILAM] );
  if ( ISED == 1 )  { debugexit("store SED"); } 
  */

  return ;

} // end of fill_SALT2_TABLE_SED





// ***********************************************
void fill_SALT2_TABLE_COLORLAW(void) {

  // Create and fill color law table as a function of
  // color and rest-frame wavelength. The lambda bins
  // are the same bins as for the SED.

  int NLAMSED, NCBIN, ilam, ic, I8, I8p;
  double LAM, CMIN, CMAX, CSTEP, CVAL, xc, CCOR ;

  char fnam[] = "fill_SALT2_TABLE_COLORLAW";
  
  // ------------- BEGIN ------------------

  // first hard-wire the binning for color ...
  // perhaps later read this from somewhere.

  
  SALT2_TABLE.NCBIN  =  401  ;
  SALT2_TABLE.CMIN   = -2.0  ;
  SALT2_TABLE.CMAX   = +2.0  ;
  SALT2_TABLE.CSTEP  =  0.01 ;


  // get local 'NBIN' variables
  NLAMSED  = SALT2_TABLE.NLAMSED ;
  NCBIN    = SALT2_TABLE.NCBIN   ;
  CMIN     = SALT2_TABLE.CMIN    ;
  CMAX     = SALT2_TABLE.CMAX    ;
  CSTEP    = SALT2_TABLE.CSTEP   ;

  printf("  Create ColorLaw Table: COLOR(MIN,MAX,STEP) = %3.1f,%2.1f,%3.2f\n",
	 CMIN,CMAX, CSTEP );

  // allocate memory for table
  I8  = sizeof(double);
  I8p = sizeof(double*);

  SALT2_TABLE.COLOR     = (double *)malloc(I8*NCBIN);
  SALT2_TABLE.COLORLAW  = (double**)malloc(I8p*NCBIN);
  for ( ic=0; ic < NCBIN; ic++ ) {

    xc = (double)ic;
    CVAL = CMIN + xc*CSTEP ;
    SALT2_TABLE.COLOR[ic] = CVAL ;
    SALT2_TABLE.COLORLAW[ic] = (double*)malloc(I8*NLAMSED);

    for ( ilam=0; ilam < NLAMSED; ilam++ ) {
      LAM  = SALT2_TABLE.LAMSED[ilam];
      CCOR = SALT2colorCor(LAM,CVAL);
      SALT2_TABLE.COLORLAW[ic][ilam] = CCOR ;
    } // ilam
  } //  ic

  // sanity checks on color-table

  if ( SALT2_TABLE.COLOR[0] != CMIN ) {
    sprintf(c1err, "SALT2_TABLE.COLOR[0] = %f", SALT2_TABLE.COLOR[0]);
    sprintf(c2err, "but should be CMIN = %f", CMIN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( SALT2_TABLE.COLOR[NCBIN-1] != CMAX ) {
    sprintf(c1err, "SALT2_TABLE.COLOR[%d] = %f", 
	    NCBIN-1, SALT2_TABLE.COLOR[NCBIN-1]);
    sprintf(c2err, "but should be CMAX = %f", CMAX);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


} // end of fill_SALT2_TABLE_colorLaw



// ***********************************************
void read_SALT2errmaps(double Trange[2], double Lrange[2] ) {
  // Mar 2011
  // Read error maps that depend on Trest vs. lambda
  // (move code from init_genmag_SALT2).
  //
  // May 2011: check array bound
  // Jul 2013: add array-bound check on NBTOT = NBLAM*NDAY
  // Jul 23 2020: for SALT3, remove _relative from file names
  // Apr 27 2021: skip reading ERRSCAL map for SALT3
  //
  // Dec 2023: refactor to read covariances for up to NSURFACE surfaces
  //    (no longer hard-wired for 2 surfaces)

  int NSED = SEDMODEL.NSURFACE ;
  int imap, NDAY, NLAM, NBTOT, nflux_nan ;
  double DUMMY[20];

  bool LEGACY = (DEBUG_SALT2 == 0 ) ;
  char tmpFile[200], sedcomment[80], lc_string[20] ;
  char *prefix = SALT2_PREFIX_FILENAME ;
  char fnam[] = "read_SALT2errmaps" ;

  // ----------- BEGIN -----------    

  if ( LEGACY ) { read_SALT2errmaps_legacy(Trange,Lrange); return; }

  printf("\n Read SALT%d ERROR MAPS: \n", IMODEL_SALT );
  fflush(stdout);

  NERRMAP_BADRANGE_SALT2 = 0 ;
  NERRMAP_BADVALUE_SALT2 = 0 ; // July 2020

  // hard-wire filenames for error maps
  if ( ISMODEL_SALT2 ) 
    { sprintf(lc_string,"lc_relative"); }
  else if ( ISMODEL_SALT3 ) 
    { sprintf(lc_string,"lc"); }

  int NERRMAP_COV =  (NSED*NSED - NSED)/2 ;
  NERRMAP_SALT2 = 
    NSED +                    // Number of variance maps
    NERRMAP_COV +             // Number of cov maps
    2;                        // scaling and disp map

  if ( NERRMAP_SALT2 >= MXERRMAP_SALT2 ) {
    sprintf(c1err,"NERRMAP_SALT2=%d exceeds bound MXERRMAP_SALT2=%d",
	    NERRMAP_SALT2, MXERRMAP_SALT2 );
    sprintf(c2err,"NSED = %d ", NSED);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // define file name, comment, and index for each error map
  int ised, ised0=0, ised1=0, icov=0 ;
  for(imap=0; imap < NERRMAP_SALT2; imap++ ) {
    if ( imap < NSED ) {
      ised = imap;
      INDEX_SALT2_ERRMAP.VAR[ised] = imap;
      sprintf(SALT2_ERRMAP_FILES[imap], "%s_%s_variance_%d.dat", 
	      prefix, lc_string, ised );
      sprintf(SALT2_ERRMAP_COMMENT[imap],  "VAR%d", ised );
    }
    else if ( imap < NSED + NERRMAP_COV ) {
      ised0++ ;
      if ( ised0 == NSED ) { ised0--; ised1++; }
      icov++ ;
      INDEX_SALT2_ERRMAP.COVAR[ised0][ised1] = imap;
      INDEX_SALT2_ERRMAP.COVAR[ised1][ised0] = imap;
      sprintf(SALT2_ERRMAP_FILES[imap], "%s_%s_covariance_%d%d.dat", 
	      prefix,lc_string, ised1, ised0 );   
      sprintf(SALT2_ERRMAP_COMMENT[imap],  "COVAR%d%d", ised1, ised0 );
    }
    else if ( imap == NSED + NERRMAP_COV ) {
      sprintf(SALT2_ERRMAP_FILES[imap], "%s_lc_dispersion_scaling.dat", prefix );
      sprintf(SALT2_ERRMAP_COMMENT[imap],  "ERRSCALE" );
      INDEX_SALT2_ERRMAP.ERRSCALE = imap;
    }

    else if ( imap == NSED + NERRMAP_COV + 1 ) {
      sprintf(SALT2_ERRMAP_FILES[imap], "%s_color_dispersion.dat",      prefix );
      sprintf(SALT2_ERRMAP_COMMENT[imap],  "COLOR-DISP" ); 
      INDEX_SALT2_ERRMAP.COLORDISP = imap;
    }

    printf("  Define error map %2d(%s): %s \n",
	   imap, SALT2_ERRMAP_COMMENT[imap], SALT2_ERRMAP_FILES[imap] );
    fflush(stdout);
    
  } // end imap

  
  for ( imap=0; imap < NERRMAP_SALT2; imap++ ) {

    init_BADVAL_SALT2errmap(imap);
     
    if ( imap >= INDEX_SALT2_ERRMAP.COLORDISP ) { continue ; } // read elsewhere

    if ( ISMODEL_SALT3 && imap == INDEX_SALT2_ERRMAP.ERRSCALE ) { continue; }

    sprintf(tmpFile, "%s/%s", SALT2_MODELPATH, SALT2_ERRMAP_FILES[imap] );
    sprintf(sedcomment, "SALT%d-%s", IMODEL_SALT, SALT2_ERRMAP_COMMENT[imap] );

    rd_sedFlux(tmpFile, sedcomment, Trange, Lrange
	       ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL, 0   // inputs
	       ,&SALT2_ERRMAP[imap].NDAY    // outputs
	       ,SALT2_ERRMAP[imap].DAY      // idem ...
	       ,&SALT2_ERRMAP[imap].DAYSTEP
	       ,&SALT2_ERRMAP[imap].NLAM
	       ,SALT2_ERRMAP[imap].LAM
	       ,&SALT2_ERRMAP[imap].LAMSTEP
	       ,SALT2_ERRMAP[imap].VALUE 
	       ,DUMMY, &nflux_nan
	       );

    NLAM = SALT2_ERRMAP[imap].NLAM ;
    SALT2_ERRMAP[imap].LAMMIN  = SALT2_ERRMAP[imap].LAM[0] ;
    SALT2_ERRMAP[imap].LAMMAX  = SALT2_ERRMAP[imap].LAM[NLAM-1] ;

    NDAY = SALT2_ERRMAP[imap].NDAY ;
    SALT2_ERRMAP[imap].DAYMIN  = SALT2_ERRMAP[imap].DAY[0] ;
    SALT2_ERRMAP[imap].DAYMAX  = SALT2_ERRMAP[imap].DAY[NDAY-1] ;

    NBTOT = NLAM*NDAY ;
    if ( NBTOT >= MXBIN_VAR_SALT2 ) {
      sprintf(c1err,"NLAM*NDAY=%d*%d = %d exceeds bound of "
	      "MXBIN_VAR_SALT2=%d",
	      NLAM, NDAY, NBTOT, MXBIN_VAR_SALT2);
      sprintf(c2err,"See '%s'", tmpFile);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }


    // Sep 2019: make sure wave range covers SED wave range
    check_lamRange_SALT2errmap(imap);
    check_dayRange_SALT2errmap(imap);
    check_BADVAL_SALT2errmap(imap);

    fflush(stdout);

  }   //  imap

  //  if ( DEBUG_SALT2 ) { debugexit(fnam); }

} // end of read_SALT2errmaps


// ***********************************************
void read_SALT2errmaps_legacy(double Trange[2], double Lrange[2] ) {
  // Mar 2011
  // Read error maps that depend on Trest vs. lambda
  // (move code from init_genmag_SALT2).
  //
  // May 2011: check array bound
  // Jul 2013: add array-bound check on NBTOT = NBLAM*NDAY
  // Jul 23 2020: for SALT3, remove _relative from file names
  // Apr 27 2021: skip reading ERRSCAL map for SALT3
  //
  int imap, NDAY, NLAM, NBTOT, nflux_nan ;
  double DUMMY[20];

  char tmpFile[200], sedcomment[80], lc_string[20] ;
  char *prefix = SALT2_PREFIX_FILENAME ;
  char fnam[] = "read_SALT2errmaps_legacy" ;

  // ----------- BEGIN -----------    

  printf("\n Read SALT2 ERROR MAPS: \n");
  fflush(stdout);

  NERRMAP_BADRANGE_SALT2 = 0 ;
  NERRMAP_BADVALUE_SALT2 = 0 ; // July 2020

  // ****** LEGACY *********

  // hard-wire filenames for error maps
  if ( ISMODEL_SALT2 ) 
    { sprintf(lc_string,"lc_relative"); }
  else if ( ISMODEL_SALT3 ) 
    { sprintf(lc_string,"lc"); }

  // ****** LEGACY *********

  INDEX_SALT2_ERRMAP.VAR[0]       = 0;
  INDEX_SALT2_ERRMAP.VAR[1]       = 1;
  INDEX_SALT2_ERRMAP.COVAR[0][1]  = 2 ;
  INDEX_SALT2_ERRMAP.COVAR[1][0]  = 2 ;
  INDEX_SALT2_ERRMAP.ERRSCALE     = 3 ;
  INDEX_SALT2_ERRMAP.COLORDISP    = 4 ;
  NERRMAP_SALT2 = 5; // hard wired

  sprintf(SALT2_ERRMAP_FILES[0], "%s_%s_variance_0.dat", prefix, lc_string );
  sprintf(SALT2_ERRMAP_FILES[1], "%s_%s_variance_1.dat", prefix, lc_string );
  sprintf(SALT2_ERRMAP_FILES[2], "%s_%s_covariance_01.dat", prefix,lc_string);
  sprintf(SALT2_ERRMAP_FILES[3], "%s_lc_dispersion_scaling.dat", prefix );
  sprintf(SALT2_ERRMAP_FILES[4], "%s_color_dispersion.dat",      prefix );

  sprintf(SALT2_ERRMAP_COMMENT[0],  "VAR0" );
  sprintf(SALT2_ERRMAP_COMMENT[1],  "VAR1" );
  sprintf(SALT2_ERRMAP_COMMENT[2],  "COVAR" );
  sprintf(SALT2_ERRMAP_COMMENT[3],  "ERRSCALE" );
  sprintf(SALT2_ERRMAP_COMMENT[4],  "COLOR-DISP" ); // 10 chars long

  // ****** LEGACY *********
  
  for ( imap=0; imap < NERRMAP_SALT2; imap++ ) {

    init_BADVAL_SALT2errmap_legacy(imap);
     
    if ( imap >= INDEX_SALT2_ERRMAP.COLORDISP ) { continue ; } // read elsewhere

    if ( ISMODEL_SALT3 && imap==INDEX_SALT2_ERRMAP.ERRSCALE ) { continue; }

    sprintf(tmpFile, "%s/%s", SALT2_MODELPATH, SALT2_ERRMAP_FILES[imap] );
    sprintf(sedcomment, "SALT2-%s", SALT2_ERRMAP_COMMENT[imap] );
    // ****** LEGACY *********

    rd_sedFlux(tmpFile, sedcomment, Trange, Lrange
	       ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL, 0   // inputs
	       ,&SALT2_ERRMAP[imap].NDAY    // outputs
	       ,SALT2_ERRMAP[imap].DAY      // idem ...
	       ,&SALT2_ERRMAP[imap].DAYSTEP
	       ,&SALT2_ERRMAP[imap].NLAM
	       ,SALT2_ERRMAP[imap].LAM
	       ,&SALT2_ERRMAP[imap].LAMSTEP
	       ,SALT2_ERRMAP[imap].VALUE 
	       ,DUMMY, &nflux_nan
	       );

    NLAM = SALT2_ERRMAP[imap].NLAM ;
    SALT2_ERRMAP[imap].LAMMIN  = SALT2_ERRMAP[imap].LAM[0] ;
    SALT2_ERRMAP[imap].LAMMAX  = SALT2_ERRMAP[imap].LAM[NLAM-1] ;

    // ****** LEGACY *********

    NDAY = SALT2_ERRMAP[imap].NDAY ;
    SALT2_ERRMAP[imap].DAYMIN  = SALT2_ERRMAP[imap].DAY[0] ;
    SALT2_ERRMAP[imap].DAYMAX  = SALT2_ERRMAP[imap].DAY[NDAY-1] ;

    NBTOT = NLAM*NDAY ;
    if ( NBTOT >= MXBIN_VAR_SALT2 ) {
      sprintf(c1err,"NLAM*NDAY=%d*%d = %d exceeds bound of "
	      "MXBIN_VAR_SALT2=%d",
	      NLAM, NDAY, NBTOT, MXBIN_VAR_SALT2);
      sprintf(c2err,"See '%s'", tmpFile);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    // ****** LEGACY *********

    // Sep 2019: make sure wave range covers SED wave range
    check_lamRange_SALT2errmap(imap);
    check_dayRange_SALT2errmap(imap);
    check_BADVAL_SALT2errmap(imap);

    fflush(stdout);

  }   //  imap

  // ****** LEGACY *********  
  return; 

} // end of read_SALT2errmaps_legacy


// ***************************************
void init_SALT2interp_SEDFLUX(void) {

  int OPT ;
  //  char fnam[] = "init_SALT2interp_SEDFLUX" ;

  // ---------- BEGIN -----------

  OPT = INPUT_SALT2_INFO.SEDFLUX_INTERP_OPT ;
  if ( OPT != SALT2_INTERP_SPLINE ) { return ; }

  // nothing to do here.

} // end of init_SALT2interp_SEDFLUX


// ***************************************
void init_SALT2interp_ERRMAP(void) {

  // if spline-option is set for error maps,
  // then init splines

  int OPT, imap, ispline, iday, ilam, N2DBIN, jtmp, IERR ;
  double ERRTMP, XM, SS ;


  char fnam[] = "init_SALT2interp_ERRMAP" ;

  // ---------- BEGIN -----------

  OPT = INPUT_SALT2_INFO.ERRMAP_INTERP_OPT ;
  if ( OPT != SALT2_INTERP_SPLINE ) { return ; }

  for ( imap=0; imap < NERRMAP_SALT2; imap++ ) {

    if ( imap >= INDEX_SALT2_ERRMAP.COLORDISP ) { continue ; }

      ispline = SALT2_TABLE.INDEX_SPLINE[1] + imap + 1 ; 
      SALT2_ERRMAP[imap].INDEX_SPLINE = ispline ; 

      SALT2_SPLINE_ARGS.DAYLIM[0] = SALT2_ERRMAP[imap].DAYMIN ;
      SALT2_SPLINE_ARGS.DAYLIM[1] = SALT2_ERRMAP[imap].DAYMAX ;
      SALT2_SPLINE_ARGS.LAMLIM[0] = SALT2_ERRMAP[imap].LAMMIN ;
      SALT2_SPLINE_ARGS.LAMLIM[1] = SALT2_ERRMAP[imap].LAMMAX ;

      // for spline, use every other day and every other lambda bin
      N2DBIN = 0;
      for ( iday=0; iday <  SALT2_ERRMAP[imap].NDAY ; iday+=2 ) {
	for ( ilam=0; ilam <  SALT2_ERRMAP[imap].NLAM ; ilam+=2 ) {
	  N2DBIN++ ;

	  jtmp = SALT2_ERRMAP[imap].NLAM *iday + ilam ;
	  ERRTMP = SALT2_ERRMAP[imap].VALUE[jtmp] ;
          if ( ERRTMP == 0.0 ) ERRTMP = 1.0E-9 ;

	  SALT2_SPLINE_ARGS.DAY[N2DBIN-1] = SALT2_ERRMAP[imap].DAY[iday] ;
	  SALT2_SPLINE_ARGS.LAM[N2DBIN-1] = SALT2_ERRMAP[imap].LAM[ilam] ;
	  SALT2_SPLINE_ARGS.VALUE[N2DBIN-1] = log10(ERRTMP*ERRTMP) ;

	}// ilam
      } // iday
     
      XM  = (double)N2DBIN ;  SS  = XM ;

      in2dex_(&ispline, &N2DBIN
	      , SALT2_SPLINE_ARGS.DAY
	      , SALT2_SPLINE_ARGS.LAM
	      , SALT2_SPLINE_ARGS.VALUE
	      , SALT2_SPLINE_ARGS.DAYLIM
	      , SALT2_SPLINE_ARGS.LAMLIM
	      , &SS, &IERR ) ; 

      printf("\t Init SPLINE %2d  for error map: %d nodes (IERR=%d) \n", 
	     ispline, N2DBIN, IERR);

      if ( IERR > 0 ) {
	sprintf(c1err,"IN2DEX SPLINE-INIT is bad: IERR=%d", IERR );
	sprintf(c2err,"ispline=%d  SS=%le \n", ispline, SS);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

    
  } // imap

} // end of init_SALT2interp_ERRMAP


// **************************************
void getFileName_SALT2colorDisp(char *fileName) {

  int imap;

  imap = INDEX_SALT2_ERRMAP.COLORDISP;
  sprintf(fileName, "%s/%s", SALT2_MODELPATH, SALT2_ERRMAP_FILES[imap] );
} 


// ***********************************************
void read_SALT2colorDisp(void) {

  // Mar 2, 2011
  // Read color dispersion vs. wavelength 
  //
  // If the color_dispersion  file returns Nrow=0, we have
  // the older Guy07 color-dispersion model, so just hard-wire
  // this map using a 3rd-order polynomial fit.
  // This hard-wire allows reading older SALT2 models
  // without needing to update the color-dispersion map.
  //
  // Sep 6, 2019: 
  //   + start ERRMAP at ilam index=0 (not 1)
  //   + call check_lamRange_SALT2errmap(imap);
  //
  // Oct 18 2022 RK - abort if all colorDisp values are zero

  int imap, NLAM, ilam, MXBIN ;
  char tmpFile[MXPATHLEN];

  // define parameters for Guy07 color dispersion model
#define NPOLY_G07 4  // 3rd-order poly fit to G07 color law
  double 
    lam, cDisp, xi,  PTMP
  ,*ptrPoly
  ,G07POLY_NULL[NPOLY_G07] = { 0.0 , 0.0, 0.0, 0.0 }
  ,G07POLY_UB[NPOLY_G07]= {6.2736, -0.43743E-02, 0.10167E-05, -0.78765E-10 }
  ,G07POLY_RI[NPOLY_G07]= {0.53882, -0.19852E-03, 0.18285E-07, -0.81849E-16 }
  ;

  int i;

  char fnam[] = "read_SALT2colorDisp";

  // ---------- BEGIN ------------

  imap = INDEX_SALT2_ERRMAP.COLORDISP ;   // refac Dec 2023

  SALT2_ERRMAP[imap].NLAM   = 0 ;  // init to nothing
  SALT2_ERRMAP[imap].LAMMIN = 0.0 ;
  SALT2_ERRMAP[imap].LAMMAX = 0.0 ;

  if ( INPUT_SALT2_INFO.ERRMAP_KCOR_OPT == 0 ) {
    printf("\n  Ignore color-dispersion (KCOR) errors. \n" );
    return ;
  }

  sprintf(tmpFile, "%s/%s", SALT2_MODELPATH, SALT2_ERRMAP_FILES[imap] );

  if ( MXBIN_VAR_SALT2 < MXBIN_LAMSED_SEDMODEL ) 
    { MXBIN = MXBIN_VAR_SALT2-1 ; }        // size of VALUE array
  else 
    { MXBIN = MXBIN_LAMSED_SEDMODEL-1 ; }  // size of LAM array


  rd2columnFile( tmpFile, MXBIN
		 ,&SALT2_ERRMAP[imap].NLAM
		 ,SALT2_ERRMAP[imap].LAM
		 ,SALT2_ERRMAP[imap].VALUE, 0   );

  NLAM = SALT2_ERRMAP[imap].NLAM ;

  printf("\n  Read color-dispersion vs. lambda from %s \n",
	 SALT2_ERRMAP_FILES[imap] );


  // check that all cDisp values are non-zero, and set ceiling.
  bool found_nonzero = false;
  for (ilam=0; ilam < NLAM; ilam++ ) {
    cDisp = SALT2_ERRMAP[imap].VALUE[ilam];
    if ( cDisp > INPUT_SALT2_INFO.COLOR_DISP_MAX ) 
      { SALT2_ERRMAP[imap].VALUE[ilam] = INPUT_SALT2_INFO.COLOR_DISP_MAX; }
    if ( cDisp > 1.0E-12 ) { found_nonzero=true; } 
  }
  
  if ( !found_nonzero ) {
    sprintf(c1err,"All colorDisp elements are zero ");
    sprintf(c2err,"Check training options." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // if nothing was read, then assume we have the older
  // Guy07 model and use polynominal parametrization
  // to hard-wire the color disp.

  if ( NLAM == 0 ) {

    // use binning of nominal surfaces
    NLAM = SALT2_TABLE.NLAMSED;
    SALT2_ERRMAP[imap].NLAM =  NLAM ;

    for ( ilam=0; ilam < NLAM; ilam++ ) {

      lam = SALT2_TABLE.LAMSED[ilam];

      if      ( lam < 4400.0 ) {  ptrPoly = G07POLY_UB   ; }
      else if ( lam < 5500.0 ) {  ptrPoly = G07POLY_NULL ; }
      else                     {  ptrPoly = G07POLY_RI   ; }

      cDisp = 0.0 ;      
      for ( i=0; i < NPOLY_G07 ; i++ ) {
	xi     = (double)i;
	PTMP   = *(ptrPoly+i) ;
	cDisp += ( PTMP * pow(lam,xi) );
      }
      
      SALT2_ERRMAP[imap].LAM[ilam]   = lam ;
      SALT2_ERRMAP[imap].VALUE[ilam] = cDisp ; 
    }

    printf("  Model is pre-G10 => hard-wire G07 color disp. \n");
  }

  SALT2_ERRMAP[imap].LAMMIN = SALT2_ERRMAP[imap].LAM[0] ;
  SALT2_ERRMAP[imap].LAMMAX = SALT2_ERRMAP[imap].LAM[NLAM-1] ;

  check_lamRange_SALT2errmap(imap);
  check_BADVAL_SALT2errmap(imap); // July 2020

  fflush(stdout);

  return ;

} // end of read_SALT2colorDisp


// =================================
void read_SALT2_INFO_FILE(int OPTMASK) {

  // March 18, 2010 R.Kessler
  // read SALT2.INFO file (or SALT3.INFO), and fill SALT2_INFO structure
  // 
  // Aug  2, 2010: read COLORLAW_VERSION: <version>
  // May  2, 2011: read SEDFLUX_INTERP_OPT 
  // Nov 24, 2011: read MAG_OFFSET
  // Oct 25, 2015: read optional RESTLAM_FORCEZEROFLUX
  // Sep 03, 2020: pass REQUIRE_DOCANA arg
  // Sep 17, 2020: read and use NPAR_POLY from COLORLAW line
  // Nov 10, 2020: read MAGSHIFT and WAVESHIFT keys
  // Nov 12, 2020: pass OPTMASK instead of REQUIRE_DOCANA
  // Aug 05, 2021: read last char for BAND
  // Oct 10, 2022: abort of NSHIFT_CALIB exceeds bount
  char
     infoFile[MXPATHLEN]
    ,c_get[60]
    ,CHAR_ERROPT[3][20] = { "OFF", "Linear", "Spline" }
    ,CHAR_SEDOPT[3][20] = { "OFF", "Linear", "Spline" }
    ,CHAR_OFFON[2][8]   = { "OFF", "ON" }
    ,ctmp[100]
       ;

  bool   REQUIRE_DOCANA, DISABLE_MAGSHIFT, DISABLE_WAVESHIFT;
  double *errtmp, *ptrpar, SHIFT_CALIB;
  double UVLAM = INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX ;
  int     OPT, ipar, NPAR_READ, NPAR_POLY, IVER, i, WHICH, NSHIFT, LEN ;
  FILE  *fp ;

  char fnam[]  = "read_SALT2_INFO_FILE"  ;

  // ------- BEGIN ---------

  REQUIRE_DOCANA     = ( OPTMASK & OPENMASK_REQUIRE_DOCANA );
  DISABLE_MAGSHIFT   = ( OPTMASK & GENMODEL_MSKOPT_SALT2_DISABLE_MAGSHIFT );
  DISABLE_WAVESHIFT  = ( OPTMASK & GENMODEL_MSKOPT_SALT2_DISABLE_WAVESHIFT );

  printf("  Read SALT2 model parameters from  \n\t  %s\n",
	 SALT2_MODELPATH );
  printf("\t OPTMASK = %d  \n", OPTMASK );
  if ( DISABLE_MAGSHIFT ) 
    { printf("\t WARNING: MAGSHIFT keys DISABLED !\n"); }
  if ( DISABLE_WAVESHIFT ) 
    { printf("\t WARNING: WAVESHIFT keys DISABLED !\n"); }
  fflush(stdout);

  sprintf(infoFile, "%s/%s", SALT2_MODELPATH, SALT2_INFO_FILE );

  check_file_docana((int)REQUIRE_DOCANA, infoFile);

  if (( fp = fopen(infoFile, "rt")) == NULL ) {
    print_preAbort_banner(fnam);
    printf("\t SALT2_MODELPATH = '%s' \n", SALT2_MODELPATH);
    sprintf(c1err,"Could not open SALT2 info file:");
    sprintf(c2err," %s", infoFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // init info variables to reasonable  default values

  INPUT_SALT2_INFO.RESTLAMMIN_FILTERCEN   = 2900. ; // Angstroms
  INPUT_SALT2_INFO.RESTLAMMAX_FILTERCEN   = 7000. ;
  INPUT_SALT2_INFO.MAGERR_FLOOR      = 0.005 ; 
  for ( i=0; i< 3; i++ ) {
    INPUT_SALT2_INFO.MAGERR_LAMOBS[i]  = 0.00 ;
    INPUT_SALT2_INFO.MAGERR_LAMOBS[i]  = 0.00 ;
  }

  // SED rebin factors if SEDFLUX_INTERP_OPT=2
  INPUT_SALT2_INFO.INTERP_SEDREBIN_LAM    = 2 ;
  INPUT_SALT2_INFO.INTERP_SEDREBIN_DAY    = 5 ;
  INPUT_SALT2_INFO.SEDFLUX_INTERP_OPT = 2 ; // 1=linear,   2=> Spline
  INPUT_SALT2_INFO.ERRMAP_INTERP_OPT  = 2 ; // 1=>linear,  2=> spline
  INPUT_SALT2_INFO.ERRMAP_KCOR_OPT    = 1 ; // 1=>ON 
  INPUT_SALT2_INFO.ERRMAP_BADVAL_ABORT= 1 ; // 1=>ON (July 2020)
  INPUT_SALT2_INFO.COLORLAW_VERSION   = IVER = 0;
  INPUT_SALT2_INFO.NCOLORLAW_PARAMS   = 4;
  INPUT_SALT2_INFO.COLOR_DISP_MAX     = COLOR_DISP_MAX_DEFAULT ;
  INPUT_SALT2_INFO.MAG_OFFSET         = 0.0 ;

  ptrpar = INPUT_SALT2_INFO.COLORLAW_PARAMS ;
  for ( ipar = 0; ipar < MXCOLORPAR ; ipar++ ) 
    { INPUT_SALT2_INFO.COLORLAW_PARAMS[ipar] = 0.0 ; }

  // Jul 2, 2010: 1st two parameters are reference wavelengths
  INPUT_SALT2_INFO.COLORLAW_PARAMS[ICLPAR_REFLAM_CL0] = B_WAVELENGTH ;
  INPUT_SALT2_INFO.COLORLAW_PARAMS[ICLPAR_REFLAM_CL1] = V_WAVELENGTH ;

  INPUT_SALT2_INFO.RESTLAM_FORCEZEROFLUX[0] = 0.0 ;
  INPUT_SALT2_INFO.RESTLAM_FORCEZEROFLUX[1] = 0.0 ;

  INPUT_SALT2_INFO.NSHIFT_CALIB = 0;

  // - - - - - - - - - -
  // read info variables

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get, "RESTLAMBDA_RANGE:") == 0 ) {
      readdouble(fp, 1, &INPUT_SALT2_INFO.RESTLAMMIN_FILTERCEN );
      readdouble(fp, 1, &INPUT_SALT2_INFO.RESTLAMMAX_FILTERCEN );

      if ( UVLAM > 0.0 )   { 
	double RESTLAMMIN_FILTERCEN_ORIG = INPUT_SALT2_INFO.RESTLAMMIN_FILTERCEN;
	double RESTLAMMIN_FILTERCEN_NEW  = UVLAM_EXTRAP_RESTLAMMIN_FILTERCEN(UVLAM);
	INPUT_SALT2_INFO.RESTLAMMIN_FILTERCEN = RESTLAMMIN_FILTERCEN_NEW ;

	printf("\t RESTLAMMIN_FILTERCEN = %.0f --> %.0f for UVLAM_EXTRAP=%.0f A\n",
	       RESTLAMMIN_FILTERCEN_ORIG, RESTLAMMIN_FILTERCEN_NEW, UVLAM);
	fflush(stdout);
      }
    }


    if ( strcmp(c_get, "COLORLAW_VERSION:") == 0 ) {
      readint(fp, 1, &IVER );
      INPUT_SALT2_INFO.COLORLAW_VERSION = IVER ;

      if ( IVER == 0 ) 
	INPUT_SALT2_INFO.NCOLORLAW_PARAMS  = 4 ;
      else if ( IVER == 1 ) 
	INPUT_SALT2_INFO.NCOLORLAW_PARAMS  = 9 ;
      else {
	sprintf(c1err,"Invalid COLORLAW_VERSION = %d", IVER );
	sprintf(c2err,"Valid versions are 0,1 only");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
    }


    if ( strcmp(c_get, "COLORLAW_PARAMS:") == 0 ||  // new key
	 strcmp(c_get, "COLORCOR_PARAMS:") == 0     // allow old key
	 ) {

      // read LAM_MIN, LAM_MAX, NPAR_POLY
      readdouble(fp, 3, &INPUT_SALT2_INFO.COLORLAW_PARAMS[ICLPAR_LAM_MIN] );
      NPAR_POLY = (int)INPUT_SALT2_INFO.COLORLAW_PARAMS[ICLPAR_NPAR_POLY];

      // read CL poly param
      readdouble(fp,NPAR_POLY,&INPUT_SALT2_INFO.COLORLAW_PARAMS[ICLPAR_POLY]);

    }
    if ( strcmp(c_get, "COLOR_OFFSET:") == 0 ) {
      readdouble(fp, 1, &INPUT_SALT2_INFO.COLOR_OFFSET );
    }
    if ( strcmp(c_get, "COLOR_DISP_MAX:") == 0 ) {
      readdouble(fp, 1, &INPUT_SALT2_INFO.COLOR_DISP_MAX );
    }

    if ( strcmp(c_get, "MAG_OFFSET:") == 0 ) {
      readdouble(fp, 1, &INPUT_SALT2_INFO.MAG_OFFSET );
    }
    if ( strcmp(c_get, "MAGERR_FLOOR:") == 0 ) {
      readdouble(fp, 1, &INPUT_SALT2_INFO.MAGERR_FLOOR );
    }
    if ( strcmp(c_get, "MAGERR_LAMOBS:") == 0 ) {
      readdouble(fp, 3, INPUT_SALT2_INFO.MAGERR_LAMOBS );
    }
    if ( strcmp(c_get, "MAGERR_LAMREST:") == 0 ) {
      readdouble(fp, 3, INPUT_SALT2_INFO.MAGERR_LAMREST );
    }

    if ( strcmp(c_get, "ERRMAP_INTERP_OPT:") == 0 ) {
      readint(fp, 1, &INPUT_SALT2_INFO.ERRMAP_INTERP_OPT );
    }

    if ( strcmp(c_get, "SEDFLUX_INTERP_OPT:") == 0 ) {
      readint(fp, 1, &INPUT_SALT2_INFO.SEDFLUX_INTERP_OPT );
    }

    if ( strcmp(c_get, "ERRMAP_KCOR_OPT:") == 0 ) {
      readint(fp, 1, &INPUT_SALT2_INFO.ERRMAP_KCOR_OPT );
    }

    if ( strcmp(c_get, "ERRMAP_BADVAL_ABORT:") == 0 ) {
      readint(fp, 1, &INPUT_SALT2_INFO.ERRMAP_BADVAL_ABORT );
    }

    if ( strcmp(c_get, "RESTLAM_FORCEZEROFLUX:") == 0 ) {
      readdouble(fp, 2, INPUT_SALT2_INFO.RESTLAM_FORCEZEROFLUX );
    }

    // Nov 10 2020: read optional calib shifts done in training
    WHICH = 0;
    if ( !DISABLE_MAGSHIFT && strcmp(c_get,"MAGSHIFT:" )==0) 
      { WHICH = CALIB_SALT2_MAGSHIFT  ; }
    if ( !DISABLE_WAVESHIFT && strcmp(c_get,"WAVESHIFT:")==0) 
      { WHICH = CALIB_SALT2_WAVESHIFT ; }
    if ( !DISABLE_WAVESHIFT && strcmp(c_get,"LAMSHIFT:")==0)  // May 21 2021
      { WHICH = CALIB_SALT2_WAVESHIFT ; }

    if ( WHICH > 0 ) {
      NSHIFT = INPUT_SALT2_INFO.NSHIFT_CALIB;

      if ( NSHIFT < MXSHIFT_CALIB_SALT2 ) {
	INPUT_SALT2_INFO.SHIFT_CALIB[NSHIFT].WHICH = WHICH ;
	readchar(fp, INPUT_SALT2_INFO.SHIFT_CALIB[NSHIFT].SURVEY_STRING );

	// e.g., if band is SDSS-r, just strip off last 
	// character 'r' and preserve full filter string
	readchar(fp, ctmp );  LEN = strlen(ctmp) ;
	readdouble(fp, 1, &SHIFT_CALIB );

	sprintf(INPUT_SALT2_INFO.SHIFT_CALIB[NSHIFT].BAND, "%c", ctmp[LEN-1] );
	sprintf(INPUT_SALT2_INFO.SHIFT_CALIB[NSHIFT].FILTER_STRING,"%s",ctmp);
	INPUT_SALT2_INFO.SHIFT_CALIB[NSHIFT].SHIFT = SHIFT_CALIB;
      }
      INPUT_SALT2_INFO.NSHIFT_CALIB++ ;
    }

  } // end while

  // check array bount (Oct 2022)
  NSHIFT = INPUT_SALT2_INFO.NSHIFT_CALIB;  
  if ( NSHIFT >= MXSHIFT_CALIB_SALT2 ) {
    sprintf(c1err,"NSHIFT_CALIB=%d exceeds bound of %d",
	    NSHIFT, MXSHIFT_CALIB_SALT2);
    sprintf(c2err,"Check %s file", SALT2_INFO_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }


  // transfer filter-lambda range to SEDMODEL struct
  SEDMODEL.RESTLAMMIN_FILTERCEN = INPUT_SALT2_INFO.RESTLAMMIN_FILTERCEN ;
  SEDMODEL.RESTLAMMAX_FILTERCEN = INPUT_SALT2_INFO.RESTLAMMAX_FILTERCEN ;

  // print INFO to screen

  printf("\n  %s \n", SALT2_INFO_FILE );
  printf("\t RESTLAMBDA_RANGE:  %6.0f - %6.0f A\n"
	 ,INPUT_SALT2_INFO.RESTLAMMIN_FILTERCEN
	 ,INPUT_SALT2_INFO.RESTLAMMAX_FILTERCEN );

  printf("\t Global MAG OFFSET:  %6.3f mag  \n", 
	 INPUT_SALT2_INFO.MAG_OFFSET ); 

  printf("\t COLOR OFFSET:  %6.3f mag  \n", 
	 INPUT_SALT2_INFO.COLOR_OFFSET ); 

  printf("\t COLOR DISP_MAX:  %5.1f mag  \n", 
	 INPUT_SALT2_INFO.COLOR_DISP_MAX ); 

  // dump colorlaw parameters based on version
  printf("\t COLORLAW PARAMS:  \n" );
  printf("\t    B,V_WAVELENGTH = %6.1f , %6.1f \n", 
	 B_WAVELENGTH,  V_WAVELENGTH );

  if ( IVER == 0 ) {
    printf("\t    Polynomial params: %f %f \n", ptrpar[2], ptrpar[3] );
  }
  else if ( IVER == 1 ) {
    printf("\t    INTERP LAMBDA RANGE: %7.1f - %7.1f \n", ptrpar[2], ptrpar[3] );
    int ipar, NPAR = (int)ptrpar[4] ;
    printf("\t    CL-polynomial params: ");
    for(ipar=0; ipar < NPAR; ipar++ ) 
      { printf("%6.3f ", ptrpar[5+ipar] ); }
    printf("\n"); fflush(stdout);
  }

  
  printf("\t MAGERR_FLOOR:  %6.3f mag  \n", INPUT_SALT2_INFO.MAGERR_FLOOR ); 

  errtmp = INPUT_SALT2_INFO.MAGERR_LAMOBS ;
  if ( *errtmp > 0.0 ) {
    printf("\t MAGERR(OBS)  += %6.3f mag for %6.0f  < LAMOBS < %6.0f \n",
	   *(errtmp+0), *(errtmp+1), *(errtmp+2) );
  }

  errtmp = INPUT_SALT2_INFO.MAGERR_LAMREST ;
  if ( *errtmp > 0.0 ) {
    printf("\t MAGERR(REST) = %6.3f mag for %6.0f < LAMREST < %6.0f \n",
	   *(errtmp+0), *(errtmp+1), *(errtmp+2) );
  }

  double  *ptrLam = INPUT_SALT2_INFO.RESTLAM_FORCEZEROFLUX;
  if ( ptrLam[1] > 0.0 ) {
    printf("\t Force flux=0 for %.0f < RESTLAM < %.0f \n",
	   ptrLam[0], ptrLam[1] );
  }

  // ---

  OPT = INPUT_SALT2_INFO.SEDFLUX_INTERP_OPT;
  
  if ( OPT == SALT2_INTERP_LINEAR ) {
    sprintf(ctmp,"%s", CHAR_SEDOPT[OPT] );
  }
  else if ( OPT == SALT2_INTERP_SPLINE ) {
    sprintf(ctmp,"%s  then Linear with LAMSTEP/%d and DAYSTEP/%d"
	    ,CHAR_SEDOPT[OPT]
	    ,INPUT_SALT2_INFO.INTERP_SEDREBIN_LAM
	    ,INPUT_SALT2_INFO.INTERP_SEDREBIN_DAY
	    );
  }
  else {
    sprintf(c1err,"Invalid SEDFLUX_INTERP_OPT = %d", OPT );
    sprintf(c2err,"Check %s file above", SALT2_INFO_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);        
  }
  
  printf("\t SEDFLUX_INTERP_OPT:  %d  (%s) \n", OPT, ctmp);


  // ---
  OPT = INPUT_SALT2_INFO.ERRMAP_INTERP_OPT;
  printf("\t ERRMAP_INTERP_OPT:   %d  (%s) \n", OPT, CHAR_ERROPT[OPT] );
  if ( OPT < 0 || OPT > 2 ) {
    sprintf(c1err,"Invalid ERRMAP_INTERP_OPT = %d", OPT );
    sprintf(c2err,"Check %s file above", SALT2_INFO_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  OPT = INPUT_SALT2_INFO.ERRMAP_KCOR_OPT;
  printf("\t ERRMAP_KCOR_OPT:     %d  (%s) \n", OPT, CHAR_OFFON[OPT] );

  OPT = INPUT_SALT2_INFO.ERRMAP_BADVAL_ABORT ;
  printf("\t ERRMAP_BADVAL_ABORT: %d  (%s) \n", OPT, CHAR_OFFON[OPT] );

  
  NSHIFT = INPUT_SALT2_INFO.NSHIFT_CALIB ;
  char KEY_SHIFT[4][12] = { "", "MAGSHIFT", "WAVESHIFT", "" } ;
  for(i=0; i < NSHIFT; i++ ) {
    WHICH = INPUT_SALT2_INFO.SHIFT_CALIB[i].WHICH ;
    printf("\t Read %s=%9.4f for SURVEY=%s, FILTER=%s\n",
	   KEY_SHIFT[WHICH],
	   INPUT_SALT2_INFO.SHIFT_CALIB[i].SHIFT,
	   INPUT_SALT2_INFO.SHIFT_CALIB[i].SURVEY_STRING,
	   INPUT_SALT2_INFO.SHIFT_CALIB[i].FILTER_STRING );
  }
  if ( NSHIFT > 0 ) { check_surveyDefined_SEDMODEL(); }

  printf("\n");    fflush(stdout);

  return ;

} // end of read_SALT2_INFO_FILE


// =========================================
void check_lamRange_SALT2errmap(int imap) {

  // Sep 6 2019
  // If imap>=0, print ERROR message if ERRMAP wave range
  // does not cover SED wave range. Also increment NERRMAP_BADRANGE_SALT2.
  // If imap < 0 && NERRMAP_BADRANGE_SALT2>0, abort.

  double SED_LAMMIN = SALT2_TABLE.LAMMIN ;
  double SED_LAMMAX = SALT2_TABLE.LAMMAX ;

  double ERRMAP_LAMMIN, ERRMAP_LAMMAX ;

  double tol     = 10.0;
  int    DISABLE = 0 ;
  char fnam[] = "check_lamRange_SALT2errmap" ;

  // ----------- BEGIN -------------

  if ( DISABLE ) { return ; }

  if ( imap < 0 ) {
    if ( NERRMAP_BADRANGE_SALT2 > 0 ) {
      sprintf(c1err,"%d ERRMAPs have invalid wavelength range.",
	      NERRMAP_BADRANGE_SALT2 );
      sprintf(c2err,"grep stdout for 'ERRMAP:'  to see all errors.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    else
      { return ; }
  }

  ERRMAP_LAMMIN = SALT2_ERRMAP[imap].LAMMIN ;
  ERRMAP_LAMMAX = SALT2_ERRMAP[imap].LAMMAX ;

  if ( ERRMAP_LAMMIN-tol > SED_LAMMIN || ERRMAP_LAMMAX+tol < SED_LAMMAX ) {
    NERRMAP_BADRANGE_SALT2++ ;
    printf("\nERRMAP: WARNING for ERRMAP file %d: %s\n", 
	   imap, SALT2_ERRMAP_FILES[imap] );
    printf("ERRMAP:     SED_LAMRANGE:    %.1f to %.1f A\n", 
	   SED_LAMMIN, SED_LAMMAX);
    printf("ERRMAP:     ERRMAP_LAMRANGE: %.1f to %.1f A "
	   "does not cover SED_LAMRANGE\n", 
	   ERRMAP_LAMMIN, ERRMAP_LAMMAX);       
  }

  return;

} // end check_lamRange_SALT2errmap

// =========================================
void check_dayRange_SALT2errmap(int imap) {

  // Sep 6 2019
  // Give error warning if ERRMAP[imap] day-range does not
  // cover SED day range, and increment NERRMAP_BADRANGE_SALT2.

  double SED_DAYMIN    = SALT2_TABLE.DAYMIN ;
  double SED_DAYMAX    = SALT2_TABLE.DAYMAX ;
  double ERRMAP_DAYMIN = SALT2_ERRMAP[imap].DAYMIN ;
  double ERRMAP_DAYMAX = SALT2_ERRMAP[imap].DAYMAX ;
  double tol = 1.1 ;
  int    DISABLE = 0 ;
  //  char fnam[] = "check_dayRange_SALT2errmap" ;

  // ----------- BEGIN -------------

  if ( DISABLE ) { return ; }

  if ( ERRMAP_DAYMIN-tol > SED_DAYMIN || ERRMAP_DAYMAX+tol < SED_DAYMAX ) {
    NERRMAP_BADRANGE_SALT2++ ;
    printf("\nERRMAP: WARNING for ERRMAP file: %s\n", 
	   SALT2_ERRMAP_FILES[imap] );
    printf("ERRMAP:     SED_DAYRANGE:    %.1f to %.1f days\n", 
	   SED_DAYMIN, SED_DAYMAX);
    printf("ERRMAP:     ERRMAP_DAYRANGE: %.1f to %.1f days "
	   "does not cover SED_DAYRANGE\n", 
	   ERRMAP_DAYMIN, ERRMAP_DAYMAX);       
  }
  
  return;

} // end check_dayRange_SALT2errmap


// ==========================================================
void  init_BADVAL_SALT2errmap(int imap) {

  // Created July 26 2020
  // Init stuff to count bad values in error maps.
  // Goal is to quickly catch retraining pathologies.
  //
  int i,i2,  NSED = SEDMODEL.NSURFACE ;
  // --------------- BEGIN ------------

  SALT2_ERRMAP[imap].NBADVAL_NAN   = 0 ;
  SALT2_ERRMAP[imap].NBADVAL_CRAZY = 0 ;
  SALT2_ERRMAP[imap].RANGE_FOUND[0] = +1.0E8 ; 
  SALT2_ERRMAP[imap].RANGE_FOUND[1] = -1.0E8 ;
  SALT2_ERRMAP[imap].RANGE_VALID[0] = -1.0E5 ; 
  SALT2_ERRMAP[imap].RANGE_VALID[1] =  1.0E5 ;

  // - - - - - - - - 
  for(i=0; i < NSED; i++ ) {
  // set valid ranges for SALT2
    if ( imap == INDEX_SALT2_ERRMAP.VAR[i] ) {
      SALT2_ERRMAP[imap].RANGE_VALID[0] =  -10.0 ; 
      SALT2_ERRMAP[imap].RANGE_VALID[1] =   500.0 ;
    }

    for(i2=0; i2 < NSED; i2++ ) {
      if ( imap == INDEX_SALT2_ERRMAP.COVAR[i][i2] ) {
	SALT2_ERRMAP[imap].RANGE_VALID[0] = -30.0 ; 
	SALT2_ERRMAP[imap].RANGE_VALID[1] =  100.0 ;
      } 
    }
  }

  // - - - -
  if ( imap == INDEX_SALT2_ERRMAP.ERRSCALE ) {
    SALT2_ERRMAP[imap].RANGE_VALID[0] =  0.0 ;
    SALT2_ERRMAP[imap].RANGE_VALID[1] =  200. ;
  }
  if ( imap == INDEX_SALT2_ERRMAP.COLORDISP ) {
    SALT2_ERRMAP[imap].RANGE_VALID[0] =  0.0 ;
    SALT2_ERRMAP[imap].RANGE_VALID[1] =  5.0 ;
  }

  // - - - - - - - - 
  // make valid range adjustments for SALT3 
  if ( ISMODEL_SALT3 ) { 
    // for D'Arcy, David, Mi ??
  }

  return ;

} // end init_check_BADVAL_SALT2errmap 

// ==========================================================
void  init_BADVAL_SALT2errmap_legacy(int imap) {

  // ******* LEGACY *********
  // Created July 26 2020
  // Init stuff to count bad values in error maps.
  // Goal is to quickly catch retraining pathologies.
  //
  // Mar 25 2021: increase color-disp crazy range to 5 (was 3)
  // Oct 18 2021: extend valid COVAR range to -30 (was -10)
  //
  // ****** Dec 28 2023: will soon be OBSOLETE ********
  // --------------- BEGIN ------------

  SALT2_ERRMAP[imap].NBADVAL_NAN   = 0 ;
  SALT2_ERRMAP[imap].NBADVAL_CRAZY = 0 ;
  SALT2_ERRMAP[imap].RANGE_FOUND[0] = +1.0E8 ; 
  SALT2_ERRMAP[imap].RANGE_FOUND[1] = -1.0E8 ;
  SALT2_ERRMAP[imap].RANGE_VALID[0] = -1.0E5 ; 
  SALT2_ERRMAP[imap].RANGE_VALID[1] =  1.0E5 ;

  // ******* LEGACY *********

  // - - - - - - - - 
  // set valid ranges for SALT2
  if ( imap == INDEX_SALT2_ERRMAP.VAR[0] ) {
    SALT2_ERRMAP[imap].RANGE_VALID[0] =  -10.0 ; 
    SALT2_ERRMAP[imap].RANGE_VALID[1] =   500.0 ;
  }
  else if ( imap == INDEX_SALT2_ERRMAP.VAR[1] ) {
    SALT2_ERRMAP[imap].RANGE_VALID[0] = -10.0 ; 
    SALT2_ERRMAP[imap].RANGE_VALID[1] =  500.0 ;
  }
  else if ( imap == INDEX_SALT2_ERRMAP.COVAR[0][1] ) {
    SALT2_ERRMAP[imap].RANGE_VALID[0] = -30.0 ; 
    SALT2_ERRMAP[imap].RANGE_VALID[1] =  100.0 ;
  }
  else if ( imap == INDEX_SALT2_ERRMAP.ERRSCALE ) {
    SALT2_ERRMAP[imap].RANGE_VALID[0] =  0.0 ;
    SALT2_ERRMAP[imap].RANGE_VALID[1] =  200. ;
  }
  else if ( imap == INDEX_SALT2_ERRMAP.COLORDISP ) {
    SALT2_ERRMAP[imap].RANGE_VALID[0] =  0.0 ;
    SALT2_ERRMAP[imap].RANGE_VALID[1] =  5.0 ;
  }

  // - - - - - - - - 
  // make valid range adjustments for SALT3 
  if ( ISMODEL_SALT3 ) { 
    // for D'Arcy, David, Mi ??
  }

  // ******* LEGACY *********
  return ;

} // end init_check_BADVAL_SALT2errmap_legacy

// ==========================================================
void  check_BADVAL_SALT2errmap(int imap) {

  // July 2020
  // check errmap values for NaN and crazy values
  // Will likely need this protection as retraining ramps up.
  //   imap >= 0 --> check entire map for bad values
  //   imap <  0 --> abort on any bad values; print stats for each map

  int    iday, ilam, jtmp, NDAY, NLAM ;
  double ERRTMP, NERR_LOCAL = 0 ;
  double *RANGE_FOUND = SALT2_ERRMAP[imap].RANGE_FOUND ;
  double *RANGE_VALID = SALT2_ERRMAP[imap].RANGE_VALID ;
  int    BADVAL_ABORT = INPUT_SALT2_INFO.ERRMAP_BADVAL_ABORT;

  double DAYMIN       = SALT2_ERRMAP[0].DAYMIN ; // fixed map index
  double DAYMAX       = SALT2_ERRMAP[0].DAYMAX ;
  double DAYEDGE_TOLERANCE = 9.0; // crazy check this far from boundary
  double DAYMIN_CRAZY = DAYMIN + DAYEDGE_TOLERANCE - 0.001 ;
  double DAYMAX_CRAZY = DAYMAX - DAYEDGE_TOLERANCE + 0.001;

  double day, lam;
  double DO_CRAZY_CHECK ;
  char fnam[] = "check_BADVAL_SALT2errmap" ;

  // ----------- BEGIN --------------

  if ( imap < 0 ) {
    printf("\n");
    printf("              NBADVAL NBADVAL       "
	   "ERRMAP-value    ERRMAP-value\n");
    printf("    ERRMAP     (NaN)  (crazy^)      "
	   "actual-range^   [valid-range]\n");
    printf("# --------------------------------------"
	   "---------------------------------\n");

    int NBAD_NAN, NBAD_CRAZY;  char *COMMENT;
    for(jtmp = 0; jtmp < NERRMAP_SALT2; jtmp++ ) {
      NBAD_NAN    = SALT2_ERRMAP[jtmp].NBADVAL_NAN; 
      NBAD_CRAZY  = SALT2_ERRMAP[jtmp].NBADVAL_CRAZY; 
      COMMENT     = SALT2_ERRMAP_COMMENT[jtmp];
      RANGE_FOUND = SALT2_ERRMAP[jtmp].RANGE_FOUND ;
      RANGE_VALID = SALT2_ERRMAP[jtmp].RANGE_VALID ;
      printf(" %10s  %5d    %5d    %8.1f - %8.1f  [%8.1f - %8.1f]\n",
	     COMMENT, NBAD_NAN, NBAD_CRAZY, 
	     RANGE_FOUND[0], RANGE_FOUND[1],
	     RANGE_VALID[0], RANGE_VALID[1]);
      fflush(stdout);
    }
    printf("# --------------------------------------"
	   "---------------------------------\n");
    printf("#         ^restricted to %.1f < DAY < %.1f days\n",
	   DAYMIN_CRAZY, DAYMAX_CRAZY);
    printf("\n"); fflush(stdout);

    if ( NERRMAP_BADVALUE_SALT2 > 0 && BADVAL_ABORT ) {
      print_preAbort_banner(fnam);       
      sprintf(c1err,"%d bad ERRMAP values (NaN and/or crazy)", 
	      NERRMAP_BADVALUE_SALT2);
      sprintf(c2err,"Check bad-value stats above");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    
    return ;
  }

  // - - - - -
  NDAY = SALT2_ERRMAP[imap].NDAY ;
  NLAM = SALT2_ERRMAP[imap].NLAM ; 

  if  (imap == INDEX_SALT2_ERRMAP.COLORDISP) { NDAY=1; } // only wave bins

  for ( iday=0; iday <  NDAY ; iday++ ) {
    for ( ilam=0; ilam < NLAM ; ilam++ ) {
      day    = SALT2_ERRMAP[imap].DAY[iday] ;
      lam    = SALT2_ERRMAP[imap].LAM[ilam] ;

      jtmp   = NLAM *iday + ilam ;
      ERRTMP = SALT2_ERRMAP[imap].VALUE[jtmp] ;

      // always check for NaN
      if ( isnan(ERRTMP) ) { 
	SALT2_ERRMAP[imap].NBADVAL_NAN++ ; NERR_LOCAL++; 
	continue ;
      }

      // check crazy values at least 5 days away from 
      // min/max day ... because the edge value really are crazy.
      DO_CRAZY_CHECK = ( day > DAYMIN_CRAZY && day < DAYMAX_CRAZY ) ;
      if ( NDAY == 1 ) { DO_CRAZY_CHECK = true; }
      if ( DO_CRAZY_CHECK ) {
	if ( ERRTMP < RANGE_VALID[0] || ERRTMP > RANGE_VALID[1] ) {
	  SALT2_ERRMAP[imap].NBADVAL_CRAZY++ ; NERR_LOCAL++; 
	}
	if ( ERRTMP < RANGE_FOUND[0] )  { RANGE_FOUND[0] = ERRTMP ; }
	if ( ERRTMP > RANGE_FOUND[1] )  { RANGE_FOUND[1] = ERRTMP ; }
      }

    } // end ilam
  }   // end iday

  NERRMAP_BADVALUE_SALT2 += NERR_LOCAL ;

  int LDMP = 0 ;
  if ( LDMP ) {
    printf(" xxx %s ---------------------------------- \n", fnam);
    printf(" xxx %s: NDAY = %d   NLAM=%d \n", fnam, NDAY, NLAM);
    printf(" xxx %s: imap=%d  RANGE = %.2f to %.2f  NERR_LOCAL=%d \n",
	   fnam, imap, RANGE_VALID[0], RANGE_VALID[1], NERR_LOCAL );
    printf(" xxx %s: RANGE_FOUND = %f to %f \n", 
	   fnam, 
	   SALT2_ERRMAP[imap].RANGE_FOUND[0],
	   SALT2_ERRMAP[imap].RANGE_FOUND[1] );
  }

  return;
} // end check_BADVAL_SALT2errmap



// =========================================
void init_calib_shift_SALT2train(void) {

  // Nov 10 2020
  // apply training-calibration shifts from SALT2.INFO file
  // Works for MAGSHIFT and WAVESHIFT keys in SALT2.INFO file.
  //
  // Notation: 
  //   _calib  -> corresponds to calib shift in the training

  int  NSHIFT_TOT   = INPUT_SALT2_INFO.NSHIFT_CALIB ;
  int  NSHIFT_APPLY = 0 ;
  int  i, which, n_survey, isurvey, ifilt, ifilt_obs, NLAM, ilam, MEMD ;
  char **survey_calib_list, *survey_calib_string, *survey_calib, *band_calib;
  char *filter_name, *filter_calib;
  double shift, magprimary, mag_shift, lam_shift;
  double *lam, *trans, *transREF ;
  bool MATCH ;
  char string_shift[3][12] = { "", "MAGSHIFT", "LAMSHIFT" } ;
  char fnam[] = "init_calib_shift_SALT2train" ;

  // ----------- BEGIN -------------

  if ( NSHIFT_TOT == 0 ) { return; }

  sprintf(BANNER,"Propagate calibration shifts from SALT2 training");
  print_banner(BANNER);

  set_FILTERSTRING(FILTERSTRING);

  for(i=0; i < NSHIFT_TOT; i++ ) {
    which               = INPUT_SALT2_INFO.SHIFT_CALIB[i].WHICH ;
    survey_calib_string = INPUT_SALT2_INFO.SHIFT_CALIB[i].SURVEY_STRING ;
    band_calib          = INPUT_SALT2_INFO.SHIFT_CALIB[i].BAND ;
    filter_calib        = INPUT_SALT2_INFO.SHIFT_CALIB[i].FILTER_STRING;
    shift               = INPUT_SALT2_INFO.SHIFT_CALIB[i].SHIFT ;
   
    mag_shift = lam_shift = 0.0 ;
    if ( which == CALIB_SALT2_MAGSHIFT ) 
      {  mag_shift = shift ; }
    else 
      {  lam_shift = shift;  }	

    // extract array of surveys for comma-sep input:
    // e.g., survey_string = 'CFA3,CFA3S,CFA3K' ->
    // survey_list = 'CFA3', 'CFA3S', 'CFA3K'
    parse_commaSepList(fnam,survey_calib_string, 10, 40,  // passed arguments 
		       &n_survey, &survey_calib_list); // returned arguments

    for(isurvey=0; isurvey < n_survey; isurvey++ ) {
      survey_calib = survey_calib_list[isurvey] ;
      
      for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++ ) {
      
	MATCH  = match_SALT2train(survey_calib, filter_calib, ifilt);
	if ( !MATCH ) { continue; }

	// store current filter trans in separate array so that
	// original array can be modified.
	NLAM = copy_filter_trans_SALT2(ifilt, &lam, &trans, &transREF);

	magprimary    = FILTER_SEDMODEL[ifilt].magprimary ;
	filter_name   = FILTER_SEDMODEL[ifilt].name ;
	ifilt_obs     = INTFILTER(filter_name);
	
	printf("\t Update %s(%2.2d) with %s = %9.4f (survey = %s)\n",
	       filter_name, ifilt_obs, string_shift[which], shift, 
	       survey_calib); 
	fflush(stdout);
	
	NSHIFT_APPLY++ ;
	init_filter_SEDMODEL(ifilt_obs, filter_name, survey_calib, 
			     magprimary+mag_shift, 
			     NLAM, lam, trans, transREF, lam_shift );
	
	free(lam); free(trans); free(transREF);
      } // end ifilt

    } // end isurvey loop
    
  } // end i loop over NSHIFT


  printf("\t --> Apply %d of %d calibration shifts.\n\n",
	 NSHIFT_APPLY, NSHIFT_TOT );  fflush(stdout);

  if ( NSHIFT_APPLY > 0 ) {  filtdump_SEDMODEL(); }

  return ;

} // end init_calib_shift_SALT2train


// =========================================
int copy_filter_trans_SALT2(int ifilt, double **lam, double **trans, 
			    double **transREF) {

  // Created Nov 11 2020
  // malloc arrays and copy filter info for 'ifilt'.

  int  NLAM          = FILTER_SEDMODEL[ifilt].NLAM ;
  int  MEMD          = NLAM * sizeof(double);
  int  ilam ;
  char fnam[] = "copy_filter_trans_SALT2" ;
  // ----------- BEGIN ------------

  *lam           = (double*) malloc(MEMD) ;
  *trans         = (double*) malloc(MEMD) ;
  *transREF      = (double*) malloc(MEMD) ;
  for(ilam=0; ilam < NLAM; ilam++ ) {
    (*lam)[ilam]      = FILTER_SEDMODEL[ifilt].lam[ilam];
    (*trans)[ilam]    = FILTER_SEDMODEL[ifilt].transSN[ilam];
    (*transREF)[ilam] = FILTER_SEDMODEL[ifilt].transREF[ilam];
  }
  
  return NLAM;

} // end copy_filter_trans_SALT2


// ========================================
bool match_SALT2train(char *survey_calib, char *filter_calib, int ifilt) {

  // Created Febr 17 2022
  // return true if input calibration shift matches ifilt.
  //
  // INPUTS:
  //     survey_calib = survey used in the training
  //     filter_calib = filter used in the training  
  //     ifilt = sparse filter of current filter for LC fit
  //
  // Oct 10 2022: fix subtle bug and treat survey_genmag as a 
  //              comma-sep list.
  //

  char *survey_genmag  = FILTER_SEDMODEL[ifilt].survey ;
  char *filter_genmag  = FILTER_SEDMODEL[ifilt].name ; // full filter name
  char  filter_calib_base[60], filter_genmag_base[60]; 
  int   j_band, j_slash ;
  bool  MATCH_SURVEY, MATCH_FILTER ;
  int   LDMP = 0 ;
  char fnam[] = "match_SALT2train";

  // ---------- BEGIN ----------

  /*
  printf(" xxx %s: survey[calib,filt] = [ %s, %s ] \n",
	 fnam, survey_calib, survey_genmag); fflush(stdout);
  */

  int n_survey_genmag=0, i_srvy;
  char **survey_genmag_list, *s;   
  parse_commaSepList(fnam, survey_genmag, 10, 100,
                     &n_survey_genmag, &survey_genmag_list );
  MATCH_SURVEY = false;
  for(i_srvy=0; i_srvy < n_survey_genmag; i_srvy++ ) {
    s = survey_genmag_list[i_srvy] ;
    if ( strcmp(survey_calib,s) == 0 ) { MATCH_SURVEY = true; }
    free(survey_genmag_list[i_srvy]);
  }
  free(survey_genmag_list);


  if ( !MATCH_SURVEY ) { return MATCH_SURVEY; }

  // get original filter string "base" from filter_genmag.
  // E.g., filter_genmag = CFA3K-B/q -> filter_genmag_base = CFA3K-B (not q)
  sprintf(filter_genmag_base, "%s", filter_genmag) ;
  j_slash      = index_charString("/", filter_genmag) ;
  if ( j_slash > 0 ) { filter_genmag_base[j_slash]=0; }

  // repeat for calib filter
  sprintf(filter_calib_base, "%s", filter_calib) ;
  j_slash      = index_charString("/", filter_calib) ;
  if ( j_slash > 0 ) { filter_calib_base[j_slash]=0; }

  // - - - - - - - - - - - 
  // special back-compatability hack for Pantheon+ that specified
  // single character band in SALT2.INFO's MAGSHIFT args
  bool HACK_BAND_CALIB_ONLY = true;
  if ( HACK_BAND_CALIB_ONLY ) {
    int lenf_genmag = strlen(filter_genmag);
    int lenf_calib  = strlen(filter_calib);
    if ( lenf_calib == 1 ) {
      filter_genmag_base[0] = filter_calib_base[0] = 0;
      sprintf(filter_genmag_base,"%c", filter_genmag[lenf_genmag-1]);
      sprintf(filter_calib_base, "%c", filter_calib[lenf_calib-1]);
    }
  }


  // check filter-string match for entire name, not just last char.
  MATCH_FILTER =  ( strcmp(filter_genmag_base,filter_calib_base) == 0 ) ;

  if ( LDMP ) {
    printf(" xxx ----------------------------------------- \n");
    printf(" xxx %s: survey_calib=%s, filter_calib=%s  base=%s\n",
	   fnam, survey_calib, filter_calib , filter_calib_base);

    printf(" xxx %s: survey_genmag=%s  filter_genmag=%s base=%s (ifilt=%d)\n",
	   fnam, survey_genmag, filter_genmag, filter_genmag_base, ifilt);

    printf(" xxx %s:  MATCH_FILTER=%d \n",
	   fnam,  MATCH_FILTER);
    fflush(stdout);
  }

  bool MATCH = MATCH_SURVEY && MATCH_FILTER ;
  return MATCH ;

} // end match_SALT2train


/**********************************************
  SALT-II color correction formula
**********************************************/
double SALT2colorCor(double lam_rest, double c ) {

  // Compute/return flux-correction from color-law.
  // The color-law version 'IVER' is from the SALT2.INFO file.
  //
  // Jul 2, 2010: 
  // replace formula with a call to generic SALT2colorlaw0 function

  double cc ;
  int IVER;
  char fnam[] = "SALT2colorCor" ;

  // -------- BEGIN ---------


  IVER = INPUT_SALT2_INFO.COLORLAW_VERSION ;
  cc   = c - INPUT_SALT2_INFO.COLOR_OFFSET ;

  if ( IVER == 0 ) {
    return SALT2colorlaw0(lam_rest, cc, INPUT_SALT2_INFO.COLORLAW_PARAMS );
  }
  else if ( IVER == 1 ) {
    return SALT2colorlaw1(lam_rest, cc, INPUT_SALT2_INFO.COLORLAW_PARAMS );
  }
  else {
    sprintf(c1err,"Invalid COLORLAW_VERSION = %d", IVER );
    sprintf(c2err,"Valid versions are 0,1 only");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    return(-9.0);
  }

  return(-9.0);

} // end of SALT2colorCor



// ****************************************************************
void genmag_SALT2(
		  int OPTMASK     // (I) bit-mask of options (LSB=0)
		  ,int ifilt_obs  // (I) absolute filter index
		  ,double *parList_SN   // x0, x1, x1_forErr, c, x2
		  ,double *parList_HOST // RV, AV, logMass
		  ,double mwebv   // (I) Galactic extinction: E(B-V)
		  ,double z       // (I) Supernova redshift
		  ,double z_forErr// (I) z used for error calc (Mar 2018)
		  ,int    Nobs         // (I) number of epochs
		  ,double *Tobs_list   // (I) list of Tobs (w.r.t peakMJD) 
		  ,double *magobs_list  // (O) observed mag values
		  ,double *magerr_list  // (O) model mag errors
		  ) {

  /****
  Return observer frame mag in absolute filter index "ifilt_obs" 
  for input SALT2 parameters.

   OPTMASK=1 => return flux instead of mag ; magerr still in mag.
                     (to avoid discontinuity for negative flux)

   OPTMASK=2 => print warning message when model flux < 0

   OPTMASK=4 => set errors to zero  (July 2013)

  May 1, 2011: in addition to major change to use brute-force integration,
               fix bug in calculation of relx1.


  Jan 27, 2013: fix Trest logic for EXTRAPFLAG,
    Trest <  SALT2_TABLE.DAYMIN  -> Trest <=  SALT2_TABLE.DAYMIN+epsT 
    Trest >  SALT2_TABLE.DAYMAX  -> Trest >=  SALT2_TABLE.DAYMAX-epsT

 Dec 27 2014: fix aweful extrapolation bug (EXTRAPFLAG>0);
              was causing constant flux/mag when Trest > model range.

 Oct 25 2015: check option to force flux=0 for specific restLam range.
              See  INPUT_SALT2_INFO.RESTLAM_FORCEZEROFLUX.

 July 2016: add arguments RV_Host and AV_host
 
 Jun 25 2018: check mag-extrap option for late times

 May 31 2021: refactor to receive parList_SN[HOST]

 Dec 28 2023: pass x2 = parlist_SN[4] if SALT3 model includes M2 component.

  ***/

  double x0        = parList_SN[0];
  double x1        = parList_SN[1];
  double c         = parList_SN[2];
  double x1_forErr = parList_SN[3];
  double x2        = parList_SN[4]; // Dec 2023

  double RV_host      = parList_HOST[0];
  double AV_host      = parList_HOST[1];
  double logMass_host = parList_HOST[2];

  double 
    meanlam_obs,  meanlam_rest, ZP, z1
    ,Tobs, Tobs_interp, Trest, Trest_interp, flux, flux_interp
    ,arg, magerr, Finteg, Finteg_errPar, FspecDum[10]
    ,lamrest_forErr, Trest_forErr, z1_forErr, magobs_tmp, magobs
    ;

  double
     fluxmin = 1.0E-30
    ,epsT    = 1.0E-5
    ;

  char *cfilt;
  int  ifilt, epobs, EXTRAPFLAG_DMP = 0;
  int  LDMP_DEBUG,OPT_PRINT_BADFLUX, OPT_RETURN_MAG ;
  int  OPT_RETURN_FLUX, OPT_DOERR ;    

  bool DO_EXTRAP_LOCAL       = false;
  bool EXTRAP_METHOD_MAG     = (EXTRAP_PHASE_METHOD == EXTRAP_PHASE_MAG );
  bool EXTRAP_METHOD_SEDFLUX = (EXTRAP_PHASE_METHOD == EXTRAP_PHASE_SEDFLUX );
  bool EXTRAP_METHOD_FLAM    = (EXTRAP_PHASE_METHOD == EXTRAP_PHASE_FLAM);

  char fnam[] = "genmag_SALT2" ;

  // ----------------- BEGIN -----------------
  
  // parse bit-mask options

  OPT_PRINT_BADFLUX = OPT_RETURN_FLUX =  0; // default flags OFF
  OPT_RETURN_MAG = OPT_DOERR = 1 ;       // default flags on
  LDMP_DEBUG=0;

  if ( (OPTMASK & 1)  ) { OPT_RETURN_MAG    = 0 ; OPT_RETURN_FLUX = 1; }
  if ( (OPTMASK & 2)  ) { OPT_PRINT_BADFLUX = 1 ; }
  if ( (OPTMASK & 4)  ) { OPT_DOERR         = 0 ; } // Jul 2013
  if ( (OPTMASK & 8)  ) { LDMP_DEBUG        = 1 ; }

  // translate absolute filter index into sparse index
  ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ;
  z1    = 1. + z ;

  // filter info for this "ifilt"
  meanlam_obs  = FILTER_SEDMODEL[ifilt].mean ;  // mean lambda
  ZP           = FILTER_SEDMODEL[ifilt].ZP ;
  cfilt        = FILTER_SEDMODEL[ifilt].name ;
  meanlam_rest = meanlam_obs/z1 ;


  // make sure filter-lambda range is valid
  checkLamRange_SEDMODEL(ifilt,z,fnam);

  // store info for Galactic & host extinction
  fill_TABLE_MWXT_SEDMODEL(MWXT_SEDMODEL.RV, mwebv);
  fill_TABLE_HOSTXT_SEDMODEL(RV_host, AV_host, z);   // July 2016

  //determine integer times which sandwich the times in Tobs

  for ( epobs=0; epobs < Nobs; epobs++ ) {

    Tobs = Tobs_list[epobs];
    Trest = Tobs / z1 ;

    DO_EXTRAP_LOCAL = false;
    Trest_interp = Trest; 

    flux = FLUX_UNDEFINED ;

    // - - - - - -
    // check for phase extrapolation
    if ( Trest <= SALT2_TABLE.DAYMIN+epsT )
      { DO_EXTRAP_LOCAL = true;  Trest_interp = SALT2_TABLE.DAYMIN+epsT ; }

    // FLAM extrapolation is done in INTEG_zSED (not here);
    // if !EXTRAP_METHOD_FLAM then enable local extrapolation
    if ( Trest >= SALT2_TABLE.DAYMAX-epsT && !EXTRAP_METHOD_FLAM ) 
      { DO_EXTRAP_LOCAL = true;  Trest_interp = SALT2_TABLE.DAYMAX-epsT ; }

    // check mag-extrap option for late times (June 25 2018)
    double DAYMIN_EXTRAP = INPUT_EXTRAP_LATETIME_Ia.DAYMIN ;
    if ( EXTRAP_METHOD_MAG && Trest > DAYMIN_EXTRAP  )    // legacy method
      { DO_EXTRAP_LOCAL = true; Trest_interp = DAYMIN_EXTRAP ; }

    // - - - - - - -

    // brute force integration
    Tobs_interp = Trest_interp * z1 ;
    INTEG_zSED_SALT2(0,ifilt_obs, z, Tobs_interp, parList_SN, parList_HOST,
		     &Finteg, &Finteg_errPar, FspecDum); // returned
    flux_interp = Finteg ;

    flux = flux_interp;

    // ------------------------
   
    // for Trest < 0, always do SEDFLUX-extrap method. For Trest > 0, check option.
    bool DO_EXTRAP_METHOD_SEDFLUX = ( Trest < 0.0 || EXTRAP_METHOD_SEDFLUX );
    if ( DO_EXTRAP_LOCAL && DO_EXTRAP_METHOD_SEDFLUX ) {
      // default is SED flux extrapolation
      double Trest_edge, Trest_tmp, flux_edge, flux_tmp, Tobs_tmp ;
      double slope_flux, nday_slope = 0.0 ;

      if ( Trest != 0.0 ) { nday_slope  = 3.0 * Trest / fabs(Trest); }

      // measure slope dTrest/dFlux using last nday_slope days of model
      Trest_edge = Trest_interp ;
      Trest_tmp  = Trest_edge - nday_slope ;
      flux_edge  = flux_interp ;
      Tobs_tmp   = Trest_tmp * z1 ;
      INTEG_zSED_SALT2(0,ifilt_obs, z, Tobs_tmp, parList_SN, parList_HOST,
		       &Finteg, &Finteg_errPar, FspecDum); // return
      flux_tmp = Finteg;
      
      slope_flux = -(flux_tmp - flux_edge)/nday_slope ;

      // extrapolate model
      flux = modelflux_extrap( Trest, Trest_edge, 
			       flux_edge, slope_flux, EXTRAPFLAG_DMP ) ;
    }

    // -------------------
    // Oct 25 2015: check option to force flux to zero
    if ( meanlam_rest > INPUT_SALT2_INFO.RESTLAM_FORCEZEROFLUX[0] &&
	 meanlam_rest < INPUT_SALT2_INFO.RESTLAM_FORCEZEROFLUX[1] ) {
      flux = 0.0 ;
    }

    // ------------------------

    if ( flux == FLUX_UNDEFINED ) {
      sprintf(c1err,"Undefined flux at Tobs(%s)=%.2f z=%f", cfilt, Tobs, z);
      sprintf(c2err,"x0=%le x1=%.3f x2=%.2f c=%.3f \n", x0, x1, x2, c);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( flux <= fluxmin || isnan(flux) ) {

      magobs = MAG_ZEROFLUX ;

      if ( OPT_PRINT_BADFLUX ) {
	printf("  genmag_SALT2 Warning:");
	printf(" Flux(%s)<%le at Trest = %6.2f => return mag=%.0f \n",
	       cfilt, fluxmin, Trest, MAG_ZEROFLUX );
      }
    }
    else{
      magobs = ZP - 2.5*log10(flux) + INPUT_SALT2_INFO.MAG_OFFSET ;

      if ( DO_EXTRAP_LOCAL && EXTRAP_METHOD_MAG && Trest > 0.0 ) { 
	// legacy method from original PLASTICC 2018
	magobs_tmp = magobs ;
	magobs = genmag_extrap_latetime_Ia(magobs_tmp, Trest, meanlam_rest); 
      }
    }
    
    // check option to return flux intead of mag;
    // preserves continutity for negative model fluxes.
    if ( OPT_RETURN_FLUX ) 
      { arg = -0.4 * magobs;  magobs = pow(TEN,arg);  }


    // load output array
    magobs_list[epobs] = magobs ;

    // ------------- DEBUG DUMP ONLY ------------------
    //    LDMP_DEBUG = ( ifilt_obs == 1 && magobs > 90.0 ) ;

    if ( LDMP_DEBUG ) {
      printf("\n xxxx ================================================= \n");
      printf(" xxxx genmag_SALT2 dump \n" ) ;
      printf(" xxxx Trest(%s) = %6.2f   LAMrest = %6.0f  z=%6.4f\n", 
	     cfilt, Trest, meanlam_rest, z );
      printf(" xxxx flux=%f   mag=%f   OPT_RETURN_FLUX=%d \n", 
	     flux, magobs_list[epobs], OPT_RETURN_FLUX );
      printf(" xxxx x1=%6.3f  c=%6.3f  Finteg=%9.3le   \n", 
	     x1, c, Finteg ) ;
      printf(" xxxx ZP=%f  mwebv=%f \n", ZP, mwebv);
      printf(" xxxx colorCor = %f\n",  SALT2colorCor(meanlam_rest,c) ) ;
      fflush(stdout);
    }
    // -------- END OF DEBUG DUMP  ------------

    // get the mag error and pass the LDMP_DEBUG flag from above.

    if ( OPT_DOERR ) {
      z1_forErr      = (1.0 + z_forErr);
      Trest_forErr   = Tobs / z1_forErr ;
      lamrest_forErr = meanlam_obs / z1_forErr ;
      magerr = SALT2magerr(Trest_forErr, lamrest_forErr, z_forErr,
			   x1_forErr, x2, Finteg_errPar, LDMP_DEBUG);
    }
    else
      { magerr = 0.0 ; }

    // load magerr onto output list
    magerr_list[epobs] = magerr ;  // load error to output array


  } // end epobs loop over epochs


  return ;

} // end of genmag_SALT2


// *****************************************
double SALT2magerr(double Trest, double lamRest, double z,
		   double x1, double x2, double Finteg_errPar, int LDMP ) {


  // Created Jun 2011 by R.Kessler
  // return mag-error for this epoch and rest-frame <lamRest>.
  //
  // Inputs:
  //   - Trest   : rest-frame epoch relative to peak brightness (days)
  //   - lamRest : <lamObs>/(1+z) = mean wavelength in rest-frame
  //   - z       : redshift
  //   - x1      : stretch parameter, multiplies M1
  //   - x2      : multiplies optional M2
  //
  //   - Finteg_errPar 
  //         : for SALT2, Finteg[1] / Finteg[0]
  //         : for SALT3, (M0 + x1*M1) 
  //
  //
  //   - LDMP : dump-and-exit flag
  //
  // Nov 7 2019: for SALT3 (retraining), set relx1=0. We don't understand
  //             the origin of this term, so scrap it for SALT3.
  // 
  // Oct 16 2020: vartot 
  //   + pass new arg Finteg_noMW
  //   + vartot_flux for SALT3 (relative for SALT2)
  //
  // Oct 01 2021: no longer set magerr=5.0 to avoid discontinuity in LC fit.

  int NSED = SEDMODEL.NSURFACE;
  int i, i2 ;
  double 
     ERRMAP[MXERRMAP_SALT2], Trest_tmp
    ,vartot_rel, vartot_flux, var[10], relsig0, relsig1
    ,covar[6][6], covtmp, rho, errscale, fracerr_snake, fracerr_kcor, fracerr_TOT
    ,magerr_model, magerr, lamObs
    ,ONE = 1.0, relx1=0.0 ;
    ;

  char fnam[] = "SALT2magerr" ;

  // ---------------- BEGIN ---------------
  
  if ( !DEBUG_SALT2 ) {
    // legacy magerr for only 2 components
    return SALT2magerr_legacy(Trest, lamRest, z, x1, Finteg_errPar, LDMP );
  }

  lamObs = lamRest * ( 1. + z );

  // Make sure that Trest is within range of map.

  if ( Trest > SALT2_ERRMAP[0].DAYMAX ) 
    { Trest_tmp = SALT2_ERRMAP[0].DAYMAX ; }
  else if ( Trest < SALT2_ERRMAP[0].DAYMIN ) 
    { Trest_tmp = SALT2_ERRMAP[0].DAYMIN ; }
  else
    { Trest_tmp = Trest ; }

  get_SALT2_ERRMAP(Trest_tmp, lamRest, ERRMAP ) ;

  // strip off the goodies
  for(i=0; i < NSED; i++ ) {
    var[i]    = ERRMAP[INDEX_SALT2_ERRMAP.VAR[i]] ;  // sigma(Si)/S0
    for(i2=i; i2 < NSED; i2++ ) {
      if ( i == i2 ) 
	{ covar[i][i] = var[i]; }
      else {
	covtmp = ERRMAP[INDEX_SALT2_ERRMAP.COVAR[i][i2]] ;
        covar[i][i2]  =	 covar[i2][i]  = covtmp ;
	}
    } 
  }

  errscale = ERRMAP[INDEX_SALT2_ERRMAP.ERRSCALE] ;  // error fudge  

  vartot_rel = vartot_flux = 0.0 ;

  if ( ISMODEL_SALT2 ) {
    // SALT2: fractional error as in  Guy's ModelRelativeError function
    double Fratio = Finteg_errPar;    // Finteg[1]/Finteg[0], no MW 
    relx1         = x1 * Fratio ;
    vartot_rel    = var[0] + var[1]*x1*x1 + (2.0 * x1* covar[0][1] ) ;
    if ( vartot_rel < 0 ) { vartot_rel = 0.01*0.01 ; } // 7/2013: follow JG 
    fracerr_snake = errscale * sqrt(vartot_rel) / fabs(ONE + relx1) ;   
  }
  else if ( ISMODEL_SALT3 ) {
    // Dave and D'Arcy's vartot has flux units (M0 + x1*M1), not relative units
    double flux_train   = Finteg_errPar ;  // M0 + x1*M1; no Gal extinc and c=0    

    if ( DEBUG_SALT2 ) {
      double x_loop[3]    = { 1.0, x1, x2 };
      vartot_flux = 0.0 ;
      for(i=0; i < NSED; i++ ) {
	for(i2=i; i2 < NSED; i2++ ) {
	  covtmp = x_loop[i] * x_loop[i2] * covar[i][i2] ;
	  vartot_flux += covtmp;
	  if ( i != i2 ) { vartot_flux += covtmp; } // add other off-diag element
	}
      }
    } // end DEBUG_SALT2

    if ( !DEBUG_SALT2 ) 
      { vartot_flux = var[0] + var[1]*x1*x1 + (2.0 * x1 * covar[0][1]) ; } // legacy

    if ( vartot_flux < 0   ) { vartot_flux = -vartot_flux ; }  // W.A.G
    if ( flux_train  < 0.0 ) { flux_train  = -flux_train  ; }  // W.A.G
    fracerr_snake = sqrt(vartot_flux) / flux_train ;
  }

  // kcor/color error is the same for SALT2,SALT3
  fracerr_kcor = SALT2colorDisp(lamRest,fnam); 

  // get total fractional  error.
  fracerr_TOT  = sqrt( pow(fracerr_snake,2.0) + pow(fracerr_kcor,2.0) ) ;

  // convert frac-error to mag-error, and load return array
  magerr_model  = (2.5/LNTEN) * fracerr_TOT ;   // exact

  // check for error fudges
  magerr = magerrFudge_SALT2(magerr_model, lamObs, lamRest );


  // ------------- DEBUG DUMP ONLY ------------------
  if ( LDMP ) {
    relsig0 = sqrt(var[0]);
    relsig1 = sqrt(var[1]);
    rho = covar[0][1] / (relsig0*relsig1) ;

    // printf("\n xxxx ================================================= \n");
    printf(" xxxx \t SALT2magerr dump \n" );

    printf(" xxxx Trest=%6.2f  lamRest = %6.0f   z=%6.4f\n", 
	   Trest, lamRest, z );

    printf(" xxx Finteg_errPar = %le \n",
	   Finteg_errPar );

    printf(" xxxx var0=%le  var1=%le  vartot_[rel,flux]=%le,%le  \n", 
	   var[0], var[1], vartot_rel, vartot_flux );

    printf(" xxxx relsig0=%f  relsig1=%f  rho=%f  scale=%f\n", 
	   relsig0, relsig1, rho, errscale );

    printf(" xxxx fracerr[snake,kcor] = %f , %f \n", 
	   fracerr_snake, fracerr_kcor );

    printf(" xxxx fracerr_TOT=%f   x1*S1/S0=%f \n", 
	   fracerr_TOT, relx1 );
    printf(" xxxx magerr(model,final) = %7.3f , %7.3f \n", 
	   magerr_model, magerr );

    //    debugexit("SALT2 MODEL DUMP");

  }

  // -------- END OF DEBUG DUMP  ------------

  return magerr ;

} // end of SALT2magerr


// *****************************************
double SALT2magerr_legacy(double Trest, double lamRest, double z,
			  double x1, double Finteg_errPar, int LDMP ) {


  // Created Jun 2011 by R.Kessler
  // return mag-error for this epoch and rest-frame <lamRest>.
  //
  // Inputs:
  //   - Trest   : rest-frame epoch relative to peak brightness (days)
  //   - lamRest : <lamObs>/(1+z) = mean wavelength in rest-frame
  //   - z       : redshift
  //   - x1      : stretch parameter.
  //   - Finteg_errPar 
  //         : for SALT2, Finteg[1] / Finteg[0]
  //         : for SALT3, (M0 + x1*M1) 
  //
  //
  //   - LDMP : dump-and-exit flag
  //
  // Nov 7 2019: for SALT3 (retraining), set relx1=0. We don't understand
  //             the origin of this term, so scrap it for SALT3.
  // 
  // Oct 16 2020: vartot 
  //   + pass new arg Finteg_noMW
  //   + vartot_flux for SALT3 (relative for SALT2)
  //
  // Oct 01 2021: no longer set magerr=5.0 to avoid discontinuity in LC fit.

  double 
     ERRMAP[MXERRMAP_SALT2], Trest_tmp
    ,vartot_rel, vartot_flux, var0, var1, relsig0, relsig1
    ,covar01, rho, errscale, fracerr_snake, fracerr_kcor, fracerr_TOT
    ,magerr_model, magerr, lamObs
    ,ONE = 1.0, relx1=0.0 ;
    ;

  char fnam[] = "SALT2magerr_legacy" ;

  // ---------------- BEGIN ---------------
  
  lamObs = lamRest * ( 1. + z );

  // ********* LEGACY **********

  // Make sure that Trest is within range of map.

  if ( Trest > SALT2_ERRMAP[0].DAYMAX ) 
    { Trest_tmp = SALT2_ERRMAP[0].DAYMAX ; }
  else if ( Trest < SALT2_ERRMAP[0].DAYMIN ) 
    { Trest_tmp = SALT2_ERRMAP[0].DAYMIN ; }
  else
    { Trest_tmp = Trest ; }

  get_SALT2_ERRMAP(Trest_tmp, lamRest, ERRMAP ) ;

  // strip off the goodies
  var0     = ERRMAP[INDEX_SALT2_ERRMAP.VAR[0]] ;  // sigma(S0)/S0
  var1     = ERRMAP[INDEX_SALT2_ERRMAP.VAR[1]] ;  // sigma(S1)/S0
  covar01  = ERRMAP[INDEX_SALT2_ERRMAP.COVAR[0][1]] ;  // 
  errscale = ERRMAP[INDEX_SALT2_ERRMAP.ERRSCALE] ;  // error fudge  

  vartot_rel = vartot_flux = 0.0 ;

  // ********* LEGACY **********

  if ( ISMODEL_SALT2 ) {
    // SALT2: fractional error as in  Guy's ModelRelativeError function
    double Fratio = Finteg_errPar;    // Finteg[1]/Finteg[0], no MW 
    relx1         = x1 * Fratio ;
    vartot_rel    = var0 + var1*x1*x1 + (2.0 * x1* covar01) ;
    if ( vartot_rel < 0 ) { vartot_rel = 0.01*0.01 ; } // 7/2013: follow JG 
    fracerr_snake = errscale * sqrt(vartot_rel) / fabs(ONE + relx1) ;   
  }
  else if ( ISMODEL_SALT3 ) {
    // Dave and D'Arcy's vartot has flux units (M0+x1*M1), not relative units
    double flux_train   = Finteg_errPar ;  // M0+x1*M1; no Gal extinc and c=0
    vartot_flux = var0 + var1*x1*x1 + (2.0 * x1* covar01) ;

    if ( vartot_flux < 0   ) { vartot_flux = -vartot_flux ; }  // W.A.G
    if ( flux_train  < 0.0 ) { flux_train  = -flux_train  ; }  // W.A.G
    fracerr_snake = sqrt(vartot_flux) / flux_train ;
  }

  // ********* LEGACY **********

  // kcor/color error is the same for SALT2,SALT3
  fracerr_kcor = SALT2colorDisp(lamRest,fnam); 

  // get total fractional  error.
  fracerr_TOT  = sqrt( pow(fracerr_snake,2.0) + pow(fracerr_kcor,2.0) ) ;

  // convert frac-error to mag-error, and load return array
  magerr_model  = (2.5/LNTEN) * fracerr_TOT ;   // exact

  // check for error fudges
  magerr = magerrFudge_SALT2(magerr_model, lamObs, lamRest );

  // ********* LEGACY **********

  // ------------- DEBUG DUMP ONLY ------------------
  if ( LDMP ) {
    relsig0 = sqrt(var0);
    relsig1 = sqrt(var1);
    rho = covar01 / (relsig0*relsig1);

    // printf("\n xxxx ================================================= \n");
    printf(" xxxx \t SALT2magerr dump \n" );

    printf(" xxxx Trest=%6.2f  lamRest = %6.0f   z=%6.4f\n", 
	   Trest, lamRest, z );

    printf(" xxx Finteg_errPar = %le \n",
	   Finteg_errPar );

    printf(" xxxx var0=%le  var1=%le  vartot_[rel,flux]=%le,%le  \n", 
	   var0, var1, vartot_rel, vartot_flux );

    printf(" xxxx relsig0=%f  relsig1=%f  rho=%f  scale=%f\n", 
	   relsig0, relsig1, rho, errscale );

    printf(" xxxx fracerr[snake,kcor] = %f , %f \n", 
	   fracerr_snake, fracerr_kcor );

    printf(" xxxx fracerr_TOT=%f   x1*S1/S0=%f \n", 
	   fracerr_TOT, relx1 );
    printf(" xxxx magerr(model,final) = %7.3f , %7.3f \n", 
	   magerr_model, magerr );

  // ********* LEGACY **********

    //    debugexit("SALT2 MODEL DUMP");

  }
    // -------- END OF DEBUG DUMP  ------------

  return magerr ;

} // end of SALT2magerr_legacy


// **********************************************
double magerrFudge_SALT2(double magerr_model,
			 double meanlam_obs, double meanlam_rest ) {

  // Mar 01, 2011 R.Kessler
  // return fudged magerr for requested fudges.
  // If no fudges are defined or used, return magerr_model
  //
  // May 11 2013: arg, add fudged errors in quadrature to magerr_model
  //              instead of replacing error.
  //          

  double magerr, floor, magerr_add, sqerr ;

  // ------- BEGIN --------

  magerr = magerr_model ; // default error is input model error


  // check global floor on error
  floor = INPUT_SALT2_INFO.MAGERR_FLOOR;
  if ( magerr < floor )  { magerr = floor ; }

  // check obs-frame fudge
  if ( meanlam_obs >= INPUT_SALT2_INFO.MAGERR_LAMOBS[1] &&
       meanlam_obs <= INPUT_SALT2_INFO.MAGERR_LAMOBS[2] ) {
    magerr_add = INPUT_SALT2_INFO.MAGERR_LAMOBS[0] ;
    sqerr  = magerr_model*magerr_model + magerr_add*magerr_add ;
    magerr = sqrt(sqerr);
  }

  // check rest-frame fudge
  if ( meanlam_rest >= INPUT_SALT2_INFO.MAGERR_LAMREST[1] &&
       meanlam_rest <= INPUT_SALT2_INFO.MAGERR_LAMREST[2] ) {
    magerr_add = INPUT_SALT2_INFO.MAGERR_LAMREST[0] ;
    sqerr = magerr_model*magerr_model + magerr_add*magerr_add ;
    magerr = sqrt(sqerr);
  }
       
  return magerr ;

} // end of magerrFudge_SALT2


// **********************************************
void INTEG_zSED_SALT2(int OPT_SPEC, int ifilt_obs, double z, double Tobs, 
		      double *parList_SN, double *parList_HOST,
		      double *Finteg, double *Finteg_errPar, 
		      double *Fspec ) {

  // May 2011
  // obs-frame integration of SALT2 flux.
  // Returns Finteg that includes  filter-trans, SALT2 SEDs, 
  // color-law, and Galactic extinction.
  // This routine samples each filter-transmission
  // grid-point, and should give a better result
  // than integrating over SED-lambda (integSALT2_SEDFLUX),
  // particularly at high redshifts.
  //
  // Finteg_errPar is needed later for error estimate;
  //   for SALT2: errPar = Finteg[1]/Finteg[0]
  //   for SALT3: errPar = Finteg with no MW and normalization to per Ang.
  //
  // Do linear interpolation here.
  // Optional splines are done in fill_SALT2_TABLE_SED
  // where the SEDs are stored with finer binning so that
  // linear interp is adequate here.
  //
  // OPT_SPEC > 0 --> return Fspec spectrum within filter band.
  //                  Array size = size of filter-trans array.
  //
  // ------------ HISTORY --------------
  //
  //
  // Jan 19 2020:
  //   replace local magSmear[ilam] with global GENSMEAR.MAGSMEAR_LIST
  //   so that it works properly with repeat function.
  //
  // Oct 2020: replace Fratio with general Finteg_errPar
  // Mar 23 2021: use get_LAMTRANS_SEDMODEL() to avoid overwriting
  //              FILTER_SEDMODEL with spectrograph info
  //
  // Mar 24 2021: 
  //   + check ALLOW_NEGFLUX_SALT2 to allow or avoid negative spectral flux.
  //
  // May 31 2021: refactor to pass parList_SN and parList_HOST
  // Aug 31 2023: use zero_NEGFLAM_SEDMODEL() util
  // Dec 28 2023: implement x2 component

  int NSED = SEDMODEL.NSURFACE;

  // strip of SN and HOST params
  double x0   = parList_SN[0];
  double x1   = parList_SN[1];
  double c    = parList_SN[2];
  double x2   = parList_SN[4];  // note that skipped [3] element is xx1
  double x_loop[3] = { 1.0, x1, x2 } ; // enable looping over SED-dependent arrays

  double RV_host   = parList_HOST[0];
  double AV_host   = parList_HOST[1];
  double m_host    = parList_HOST[2];

  int  
    ifilt, NLAMFILT, ilamobs, ilamsed, jlam
    ,IDAY, NDAY, nday, iday, ised, ic
    ,ISTAT_GENSMEAR, LABORT_ILAM, LABORT_FRAC, LDMP
    ;

  double
    LAMOBS, LAMSED, z1, LAMDIF, LAMSED_MIN, LAMSED_MAX
    ,LAMFILT_STEP, LAMSED_STEP, LAMSPEC_STEP, LAMRATIO
    ,DAYSTEP, DAYMIN, DAYDIF, Trest, Trest_model
    ,MWXT_FRAC, HOSTXT_FRAC, CCOR, CCOR_LAM0, CCOR_LAM1, CDIF, CNEAR
    ,FRAC_INTERP_DAY, FRAC_INTERP_COLOR, FRAC_INTERP_LAMSED
    ,TRANS, MODELNORM_Fspec, MODELNORM_Finteg, *ptr_FLUXSED[4][4] 
    ,FSED[4], FTMP, FDIF, VAL0, VAL1, mean, arg, FSMEAR, *lam
    ,Finteg_filter[4], Finteg_forErr[4], Finteg_spec[4]
    ,Fbin_forFlux, Fbin_forSpec, Fnorm_SALT3, Fcheck, Ftmp
    ,Flam_filter[4], Flam_err[4], Flam_spec[4], parList_genSmear[10]
    ,hc8 = (double)hc ;

  bool zero_FLAM;

  int  DO_SPECTROGRAPH = ( ifilt_obs == JFILT_SPECTROGRAPH ) ;

  bool   DO_EXTRAP_LOCAL       = false;
  bool   EXTRAP_METHOD_FLAM    = (EXTRAP_PHASE_METHOD == EXTRAP_PHASE_FLAM);
  double DAYMIN_EXTRAP         = INPUT_EXTRAP_LATETIME_Ia.DAYMIN ;

  char *cfilt ;
  char fnam[] = "INTEG_zSED_SALT2" ;

  // ----------- BEGIN ---------------

  *Finteg = *Finteg_errPar = 0.0 ;
  Fspec[0] = 0.0 ; // init only first element

  for(ised=0; ised < 3; ised++ )  { 
    Finteg_filter[ised]  = Finteg_forErr[ised] = 0.0 ;  
    Flam_filter[ised] = Flam_spec[ised] =  Flam_err[ised] = 0.0 ;
  }

  Fnorm_SALT3 = 0.0 ; // for SALT3

  ifilt     = IFILTMAP_SEDMODEL[ifilt_obs] ;
  NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;
  cfilt     = FILTER_SEDMODEL[ifilt].name ;
  z1        = 1. + z ;
  Trest     = Tobs/z1 ;
  Trest_model = Trest ; // use this for SALT2 model flux

  if ( EXTRAP_METHOD_FLAM ) { // Sep 20, 2023
    if ( Trest > DAYMIN_EXTRAP ) 
      { DO_EXTRAP_LOCAL = true ; Trest_model = DAYMIN_EXTRAP ; }
  }


  LAMFILT_STEP = FILTER_SEDMODEL[ifilt].lamstep; 
  LAMSED_STEP  = SALT2_TABLE.LAMSTEP ;    // step size of SALT2 model

  // Compute flux normalization factor.
  // Note that the 1+z factor is missing because the integration 
  // is over observer-lambda instead of lambda-rest.
  MODELNORM_Fspec  = LAMFILT_STEP * SEDMODEL.FLUXSCALE ;
  MODELNORM_Finteg = LAMFILT_STEP * SEDMODEL.FLUXSCALE / hc8 ;

  // for SED find rest-frame 'iday' and DAYFRAC used to 
  // interpolate SED in TREST-space.
  DAYSTEP = SALT2_TABLE.DAYSTEP ;
  DAYMIN  = SALT2_TABLE.DAY[0]  ;
  DAYDIF  = Trest_model - DAYMIN ;
  IDAY    = (int)(DAYDIF/DAYSTEP);  
  NDAY    = SALT2_TABLE.NDAY ;
  DAYDIF  = Trest_model - SALT2_TABLE.DAY[IDAY] ;

  nday    = 2 ; 
  FRAC_INTERP_DAY = DAYDIF/DAYSTEP ;
  
  // get color-index needed for interpolation of table.
  CDIF  = c - SALT2_TABLE.CMIN ;
  ic    = (int)(CDIF / SALT2_TABLE.CSTEP) ;
  // make sure that 'ic' is within bounds.
  if ( ic < 0 ) 
    { ic = 0 ; }
  if ( ic > SALT2_TABLE.NCBIN - 2 ) 
    { ic = SALT2_TABLE.NCBIN - 2 ; }

  CNEAR = SALT2_TABLE.COLOR[ic] ;
  FRAC_INTERP_COLOR = (c - CNEAR)/SALT2_TABLE.CSTEP ;

  // get rest-frame SED pointers
  for(ised=0; ised < NSED ; ised++ ) {
    for ( iday=0; iday<nday; iday++ ) {
      ptr_FLUXSED[ised][iday] = SALT2_TABLE.SEDFLUX[ised][IDAY+iday] ;
    }
  }

  // evaluate optional smearing from function
  // Should be used only for simulation (not for fitting mode)
  ISTAT_GENSMEAR = istat_genSmear();  
  if ( ISTAT_GENSMEAR  ) {
    lam = (double*) malloc(NLAMFILT*sizeof(double) );
    //  printf(" xxx %s: z=%.3f ifilt_obs=%d \n", fnam, z, ifilt_obs); 
    int NLAMTMP = 0 ;

    parList_genSmear[0] = Trest_model;
    parList_genSmear[1] = x1;
    parList_genSmear[2] = c;
    parList_genSmear[3] = m_host ; 

    for ( ilamobs=0; ilamobs < NLAMFILT; ilamobs++ ) {

      get_LAMTRANS_SEDMODEL(ifilt,ilamobs, &LAMOBS, &TRANS);
      LAMSED       = LAMOBS/z1;   // rest-frame wavelength

      // protect undefined red end for low-z (July 2016)
      if ( LAMSED >= SALT2_TABLE.LAMMAX ) { continue ; }  

      lam[ilamobs] = LAMSED ; 
      NLAMTMP++ ;
    }

    get_genSmear(parList_genSmear, NLAMTMP, lam, GENSMEAR.MAGSMEAR_LIST) ;
    free(lam);
  }


  // Loop over obs-filter lambda-bins. XTMW has the same binning,
  // but the color and SED flux must be interpolated.

  for ( ilamobs=0; ilamobs < NLAMFILT; ilamobs++ ) {

    // fetch LAM and TRANS with utility to account for spectrograph
    get_LAMTRANS_SEDMODEL(ifilt,ilamobs, &LAMOBS, &TRANS);

    if ( TRANS < 1.0E-12 && OPT_SPEC==0) 
      { continue ; } // Jul 2013 - skip zeros for leakage

    MWXT_FRAC  = SEDMODEL_TABLE_MWXT_FRAC[ifilt][ilamobs] ;

    // July 2016: check for host extinction.    
    if( RV_host > 1.0E-9 && AV_host > 1.0E-9 ) 
      { HOSTXT_FRAC = SEDMODEL_TABLE_HOSTXT_FRAC[ifilt][ilamobs] ; }
    else 
      { HOSTXT_FRAC = 1.0 ; } // standard SALT2 model has no host extinction

    LAMSED     = LAMOBS / z1 ;  // rest-frame lambda
    LAMSED_MIN = LAMSED_MAX = LAMSED ;  // default is no sub-bins 

    // Jan 2021: bail if outside model range 
    if ( LAMSED <= SALT2_TABLE.LAMMIN ) { continue ; }
    if ( LAMSED >= SALT2_TABLE.LAMMAX ) { continue ; } 

    LDMP = 0; // (OPT_SPEC>0 && ifilt_obs==2 );

    // check spectrum options
    if ( OPT_SPEC > 0 ) {
      for(ised=0; ised<=1; ised++ ) { Finteg_spec[ised] = 0.0 ; }
      if ( DO_SPECTROGRAPH )  {
	// prepare sub-bins since SPECTROGRAPH bins can be large
	LAMSED_MIN = SPECTROGRAPH_SEDMODEL.LAMMIN_LIST[ilamobs]/z1 ; 
	LAMSED_MAX = SPECTROGRAPH_SEDMODEL.LAMMAX_LIST[ilamobs]/z1 ;
      }
    } // end OPT_SPEC

    // loop over rest-frame lambda (for SPECTROGRAPH)
    for(LAMSED = LAMSED_MIN; LAMSED <= LAMSED_MAX; LAMSED+=LAMSED_STEP ) {

      // bail if outside model range 
      if ( LAMSED <= SALT2_TABLE.LAMMIN ) { continue ; }
      if ( LAMSED >= SALT2_TABLE.LAMMAX ) { continue ; } 

      // get rest-frame lambda index and interp-fraction for SED space
      LAMDIF  = LAMSED - SALT2_TABLE.LAMMIN ;
      ilamsed = (int)(LAMDIF/LAMSED_STEP); 
      LAMDIF  = LAMSED - SALT2_TABLE.LAMSED[ilamsed] ;
      FRAC_INTERP_LAMSED = LAMDIF / LAMSED_STEP ; // 0-1


      if ( LDMP ) { 
	printf(" xxx -------------- %s DUMP ------------- \n", fnam ); 
	printf(" xxx LAMOBS=%.1f  LAMSED=%.2f \n", LAMOBS, LAMSED ); 
	printf(" xxx FRAC_INTERP_[CCOR,LAMSED] = %.3f , %.3f \n",	       
	       FRAC_INTERP_COLOR , FRAC_INTERP_LAMSED ); 
	printf(" xxx Tobs=%.3f  Trest_model=%.3f \n", Tobs, Trest_model);
	fflush(stdout);
      }

      LABORT_ILAM = ( ilamsed < 0 || ilamsed >= SALT2_TABLE.NLAMSED );     
      LABORT_FRAC = ( FRAC_INTERP_LAMSED < -1.0E-8 || 
		      FRAC_INTERP_LAMSED > 1.0000000001 ) ;

      if ( LABORT_ILAM || LABORT_FRAC ) {
	mean = FILTER_SEDMODEL[ifilt].mean ;
	print_preAbort_banner(fnam);
	printf("\t IS_ABORT(FRAC_INTERP)=%d  IS_ABORT(ilamsed)=%d \n",
	       LABORT_FRAC, LABORT_ILAM);
	printf("\t LAMOBS = %7.2f  LAMDIF=%7.2f\n",  LAMOBS, LAMDIF);
	printf("\t LAMSED = LAMOBS/(1+z) = %7.2f \n", LAMSED );
	printf("\t LAMSTEP=%4.1f  LAMMIN=%6.1f \n", 
	       LAMSED_STEP, SALT2_TABLE.LAMMIN );
	printf("\t ilamobs=%d   ilamsed= %d \n", 	     
	       ilamobs, ilamsed );
	printf("\t Tobs=%f  Trest=%f  Trest_model=%f \n", Tobs, Trest, Trest_model);
	printf("\t <LAMFILT(%s)> = %7.2f(OBS)  %7.2f(REST) \n", 
	       cfilt, mean, mean/z1);
	for( jlam=ilamsed-2; jlam <= ilamsed+2; jlam++ ) {
	  printf("\t SALT2_TABLE.LAMSED[ilamsed=%d] = %f\n", 
		 jlam, SALT2_TABLE.LAMSED[jlam] ); 
	}
	fflush(stdout);
	if ( LABORT_FRAC )  { 
	  sprintf(c1err,"Invalid FRAC_INTERP_LAMSED=%le ", 
		  FRAC_INTERP_LAMSED ); 
	}
	else if ( LABORT_ILAM ) {
	  sprintf(c1err,"Invalid ilamsed=%d (NLAMSED=%d)", 
		    ilamsed, SALT2_TABLE.NLAMSED);
	}
	sprintf(c2err,"Tobs(%s) = %6.2f at z=%5.3f x1=%.2f  c=%6.3f",
		cfilt, Tobs, z, x1, c);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      // interpolated color correction in 2D space of color and LAMSED
      VAL0  = SALT2_TABLE.COLORLAW[ic+0][ilamsed];
      VAL1  = SALT2_TABLE.COLORLAW[ic+1][ilamsed];
      CCOR_LAM0  = VAL0 + (VAL1-VAL0) * FRAC_INTERP_COLOR ;
      
      VAL0  = SALT2_TABLE.COLORLAW[ic+0][ilamsed+1];
      VAL1  = SALT2_TABLE.COLORLAW[ic+1][ilamsed+1];
      CCOR_LAM1  = VAL0 + (VAL1-VAL0) * FRAC_INTERP_COLOR ;
      
      CCOR = CCOR_LAM0 + (CCOR_LAM1-CCOR_LAM0)*FRAC_INTERP_LAMSED ;

      // interpolate SED Fluxes to LAMSED
      for(ised=0; ised < NSED; ised++ ) {
	for ( iday=0; iday<nday; iday++ ) {
	  
	  VAL0 = ptr_FLUXSED[ised][iday][ilamsed+0] ;
	  VAL1 = ptr_FLUXSED[ised][iday][ilamsed+1] ;
	  FSED[iday] = VAL0 + (VAL1-VAL0)*FRAC_INTERP_LAMSED ;
	  if ( LDMP ) {
	    printf(" xxx ised=%d iday=%d : VAL0,1=%f,%f  FSED=%f \n",
		   ised, iday, VAL0, VAL1, FSED[iday] );
	  }
	} // iday	
	
	FDIF = FSED[1] - FSED[0] ;
	FTMP = FSED[0] + FDIF * FRAC_INTERP_DAY ; 
	
	// check option to smear SALT2 flux with intrinsic scatter
	if ( ISTAT_GENSMEAR ) {
	  arg     =  -0.4*GENSMEAR.MAGSMEAR_LIST[ilamobs] ; 
	  FSMEAR  =  pow(TEN,arg)  ;        // fraction change in flux
	  FTMP   *=  FSMEAR;                // adjust flux for smearing
	}
	
	// update integral for each SED surface
	Fbin_forFlux = (FTMP * CCOR * HOSTXT_FRAC*MWXT_FRAC * LAMSED*TRANS);
	Fbin_forSpec = (FTMP * CCOR * HOSTXT_FRAC*MWXT_FRAC );
	
	// extrapolate past 45 days here ??
	if ( DO_EXTRAP_LOCAL ) {
	  Ftmp = Fbin_forFlux;
	  Fbin_forFlux = genflux_extrap_latetime_Ia(Ftmp, Trest, LAMSED);

	  Ftmp = Fbin_forSpec;
	  Fbin_forSpec = genflux_extrap_latetime_Ia(Ftmp, Trest, LAMSED);
	}
	
	if ( OPT_SPEC ) { 
	  LAMSPEC_STEP = LAMFILT_STEP ; // default for filters

	  // switch to SED bin size for SPECTROGRAPH
	  if ( DO_SPECTROGRAPH ) {
	    if ( LAMSED+LAMSED_STEP < LAMSED_MAX ) 
	      { LAMSPEC_STEP = LAMSED_STEP  ; } // obs-frame lamStep
	    else
	      { LAMSPEC_STEP = (LAMSED_MAX-LAMSED) ; }
	  }

	  LAMRATIO            = LAMSPEC_STEP/LAMFILT_STEP ; // binSize ratio
	  Flam_spec[ised]     = (Fbin_forSpec * LAMRATIO );

	} // end OPT_SPEC

	Flam_filter[ised]     = Fbin_forFlux ;
	Flam_err[ised]        = (Fbin_forFlux/MWXT_FRAC) ;  

      } // ised


      // check option to force negative flux to zero
      Fcheck = 0.0 ;
      for(ised=0; ised < NSED; ised++ )
	{ Fcheck += x_loop[ised] * Flam_filter[ised] ; }

      zero_FLAM = zero_NEGFLAM_SEDMODEL(0, LAMSED, Trest, Fcheck);
      if ( zero_FLAM ) {
	Flam_filter[0] = 1.0E-15; Flam_filter[1] = Flam_filter[2] = 0.0 ; // very small
	Flam_spec[0]   = 1.0E-15; Flam_spec[1]   = Flam_spec[2]   = 0.0 ;
      }

      for(ised=0; ised < NSED; ised++ ) {
	Finteg_filter[ised]  +=  Flam_filter[ised];
	Finteg_forErr[ised]  +=  Flam_err[ised];
	if(OPT_SPEC) { Finteg_spec[ised] +=  Flam_spec[ised]; }
      }

    } // end LAMSED loop 

    if ( OPT_SPEC ) {
      Fspec[ilamobs]  = 0.0 ;
      for(ised=0; ised < NSED; ised++ ) 
	{ Fspec[ilamobs]  += (x_loop[ised] * Finteg_spec[ised]) ; }
      Fspec[ilamobs] *= x0;
      Fspec[ilamobs] *=  (MODELNORM_Fspec) ;
    }
   
    Fnorm_SALT3  += (TRANS * LAMOBS ); 

  } // end ilamobs loop over obs filter

 
  // - - - - - - - - - - 
  // compute total flux in filter
  for(ised=0; ised < NSED; ised++ ) 
    { *Finteg  +=  ( x_loop[ised] * Finteg_filter[ised]) ; }
  *Finteg *= x0 ;
  *Finteg *= (MODELNORM_Finteg) ;

  // - - - - - - -
  // determine Finteg_errPar based on model

  if ( ISMODEL_SALT2 ) {
    if ( Finteg_filter[0] != 0.0 ) 
      { *Finteg_errPar = Finteg_forErr[1] / Finteg_forErr[0] ; }
  }
  else if ( ISMODEL_SALT3 ) {
    // exclude x0 and MODELNORM; instead, normalize to per Angstrom
    // following K21
    for(ised=0; ised < NSED; ised++ )
      { *Finteg_errPar  += ( x_loop[ised] * Finteg_forErr[ised] ); }
    *Finteg_errPar /= Fnorm_SALT3 ;
  }

  return ;

} // end of INTEG_zSED_SALT2


// **********************************************
double SALT2x0calc(
		   double alpha   // (I)
		   ,double beta   // (I)
		   ,double x1     // (I)
		   ,double c      // (I)
		   ,double dlmag  // (I) distance modulus
		   ) {

  // April 2, 2009 R.Kessler
  // Translate luminosity and color parameters into x0.

  double x0, x0inv, arg ;

  // -------- BEGIN ---------
  arg     = 0.4 * ( dlmag - alpha*x1 + beta*c ) ;
  x0inv   = X0SCALE_SALT2 * pow(TEN,arg); 
  x0      = 1./x0inv ;
  return  x0;  

}  // end of SALT2x0calc


// ***********************************
void load_mBoff_SALT2(void) {

  // Created July 27, 2010 by R.Kessler
  // Fill global mBoff_SALT2 used to compute mB.
  // Note that in the simulation mB is not used to 
  // generate fluxes but is used only as a reference.
  //
  // If B filter exists then
  //    mBoff_SALT2 = ZP  - 2.5*log10(S0)
  // where S0 is the zero-surface integral with z = Zat10pc,
  // and then  mB = mBoff_SALT2 - 2.5*log10(x0)
  //
  // Aug 11, 2010: just hard-wire mBoff to have same offset
  //               regardless of filter system.

  //  char fnam[] = "load_mBoff_SALT2" ;

  // -------- BEGIN --------

  // Aug 11, 2010: hard-wire to value based on SNLS VEGA system.
  mBoff_SALT2 = 10.635 ;
  printf("\t mB = %7.4f - 2.5*log10(x0)  \n", mBoff_SALT2 );


} // end of load_mBoff_SALT2

// ***********************************
double SALT2mBcalc(double x0) {

  // April 12, 2009 R.Kessler
  // Translate x0 into mB.
  //
  // July 27, 2010: use mBoff_SALT2 instead of hard-wired 10.63

  double mB;

  // -------- BEGIN ---------

  mB = mBoff_SALT2 - 2.5*log10(x0);
  return mB;

}  // end of SALT2mBcalc


// ***********************************************
void get_SALT2_ERRMAP(double Trest, double Lrest, double *ERRMAP ) {

  /***********************************************
   Apr 14, 2009: 
   return error values from each of the NERRMAP maps.
   Trest         :  (I) rest-frame epoch (days,  T=0 at peak)
   Lrest         :  (I) rest-frame wavelength (A)
   *ERRMAP       :  (O) error-map values: var0, var1, covar01, scale

   Aug 27, 2009: 
      interpolate in both Trest & lambda (instead of just lambda).
      Still not always continuous, but chi2-kinks are much smaller
      than before.

  Jun 2, 2011: renamed from get_SALT2modelerr to get_SALT2_ERRMAP().

  Sep 9 2019: 
    + protect iday_min and ilam_min from being negative. Negative indices
      can occur because ERRMAPs don't always cover SED range.
           
  Apr 27 2021:
     for SALT3, ignore error-scale map and hard wired scale=1

  *****************************************************/

  int imap, jval, iday_min, iday_max, ilam_min, ilam_max ;
  int NLAM, NDAY, IND, IERR ;

  double val, val0, val1, valdif, val_linear, val_spline, tmp;
  double LMIN, LSTEP, LDIF, TMIN, TSTEP, TDIF, val_atlammin, val_atlammax ;

  char fnam[] = "get_SALT2_ERRMAP";

  // ------------ BEGIN --------

  for ( imap=0; imap < NERRMAP_SALT2; imap++ ) {

    if ( imap >= INDEX_SALT2_ERRMAP.COLORDISP ) { continue ; }

    // 4.2021: there is no error-scale map for SALT3, so hard wired scale=1
    if ( ISMODEL_SALT3 && imap == INDEX_SALT2_ERRMAP.ERRSCALE )
      { ERRMAP[imap] = 1.0 ; continue ; }

    LMIN  = SALT2_ERRMAP[imap].LAMMIN ;
    LSTEP = SALT2_ERRMAP[imap].LAMSTEP ;

    TMIN  = SALT2_ERRMAP[imap].DAYMIN ;
    TSTEP = SALT2_ERRMAP[imap].DAYSTEP ;

    NLAM  = SALT2_ERRMAP[imap].NLAM ;
    NDAY  = SALT2_ERRMAP[imap].NDAY ;

    // get indices that sandwhich Trest and Lrest

    iday_min = (int)((Trest - TMIN)/TSTEP) ;
    if ( iday_min >= NDAY-1 ) { iday_min = NDAY - 2 ; }
    if ( iday_min <  0      ) { iday_min = 0; } // Sep 9 2019
    iday_max = iday_min + 1;

    ilam_min = (int)((Lrest - LMIN)/LSTEP) ;
    if ( ilam_min >= NLAM-1 ) { ilam_min = NLAM - 2 ; }
    if ( ilam_min <  0      ) { ilam_min = 0;         }
    ilam_max = ilam_min + 1;
    
    // Aug 27, 2009: 
    // interpolate Trest at LAM-MIN
    jval  = NLAM*iday_min + ilam_min ;
    val0  = SALT2_ERRMAP[imap].VALUE[jval];
    jval  = NLAM*iday_max + ilam_min ;
    val1  = SALT2_ERRMAP[imap].VALUE[jval];
    TDIF  = Trest - SALT2_ERRMAP[imap].DAY[iday_min];
    val_atlammin  = val0 + (val1-val0) * TDIF/TSTEP ;

    // interpolate Trest at LAM-MAX
    jval  = NLAM*iday_min + ilam_max ;
    val0  = SALT2_ERRMAP[imap].VALUE[jval];
    jval  = NLAM*iday_max + ilam_max ;
    val1  = SALT2_ERRMAP[imap].VALUE[jval];
    TDIF  = Trest - SALT2_ERRMAP[imap].DAY[iday_min];
    val_atlammax  = val0 + (val1-val0) * TDIF/TSTEP ;

    // interpolate in lambda space
    LDIF       = Lrest - SALT2_ERRMAP[imap].LAM[ilam_min];
    valdif     = val_atlammax - val_atlammin ;
    val_linear = val_atlammin + (valdif * LDIF/LSTEP) ;
    val        = val_linear ;

    if ( INPUT_SALT2_INFO.ERRMAP_INTERP_OPT == 0 ) 
      { val = 0.0 ; }

    if ( INPUT_SALT2_INFO.ERRMAP_INTERP_OPT == 2 ) {

      IND    = SALT2_ERRMAP[imap].INDEX_SPLINE ; 
      tmp    = ge2dex_ ( &IND, &Trest, &Lrest, &IERR ) 	;
      val_spline  = sqrt(pow(TEN,tmp)) ;

      // Use the sign of the linear interp because
      // the sign in storing the error-squared is lost.

      if ( val_linear < 0.0 ) { val = -val_spline ; }
      else                    { val = +val_spline ; }

      /*
      NCALL_DBUG_SALT2++ ;
      printf(" xxx imap=%d  Trest=%6.2f Lrest=%6.0f : ",
	     imap, Trest, Lrest );
      printf("val(lin,spline)= %le, %le \n", val_linear, val );
      if ( NCALL_DBUG_SALT2 > 50 )  debugexit("SALT2 sline");
      */

    }  // SPLINE option

    ERRMAP[imap] = val ;
    
  } // end of imap loop

} // end of get_SALT2_ERRMAP


// *******************************************************
int gencovar_SALT2(int MATSIZE, int *ifiltobsList, double *epobsList, 
		   double z, double *parList_SN, double *parList_HOST,
		   double mwebv, double *covar ) {

  // Jun 2, 2011 R.Kessler
  // return *covar matrix that depends on ifilt_obs and redshift. 
  // *covar is the covariance in mag^2 units,
  //        FAC * k(lambda)^2 
  // where FAC converts flux-fraction-squared error into mag^2 error,
  // and k(lambda) are from the salt2_color_dispersion.dat file.
  // Input 'matsize' is the size of one row or column;
  // the output *covar size is matsize^2.
  //
  // Jul 3 2013: 
  //  fix aweful bug that was double-counting the kcor error
  //  for the diagonal elements.
  //    COV_TMP += COV_DIAG  -> COV_TMP = COV_DIAG
  //  ARRRRRRRRRGH !!!
  //
  //  July 2016: add new inputs args RV_host & AV_host

  int  icovar, irow, icol, ifilt_obs, ifilt_row, ifilt_col, ifilt ;
  int ISDIAG, LDMP ;

  double x1    = parList_SN[1] ;
  double x2    = parList_SN[4] ;
  double z1    = 1.0 + z;
  double invZ1 = 1.0/z1;

  double 
    COV_TMP,  COV_DIAG, meanlam_obs, meanlam_rest
    ,cDisp[MXFILT_SEDMODEL]
    ,Finteg, Finteg_errPar, FspecDum[10], magerr
    ,Tobs, Trest, Trest_tmp, Trest_row, Trest_col
    ,FAC = 1.17882   //  [ 2.5/ln(10) ]^2
    ;

    char *cfilt, cdum0[40], cdum1[40];

    char fnam[] = "gencovar_SALT2" ;

  // -------------- BEGIN -----------------
  
  icovar = 0 ;

  // init  cDisp to -9 in each filter
  for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {
    ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs ;
    cDisp[ifilt_obs] = -9.0 ;
  }


  for ( irow=0; irow < MATSIZE; irow++ ) {
    for ( icol=0; icol < MATSIZE; icol++ ) {

      Tobs  = *(epobsList+irow) ;  Trest_row  = Tobs * invZ1 ;
      Tobs  = *(epobsList+icol) ;  Trest_col  = Tobs * invZ1 ;
      ISDIAG = ( irow == icol ) ;

      meanlam_rest  = -9.0 ;
      ifilt_row     = *(ifiltobsList+irow);
      ifilt_col     = *(ifiltobsList+icol);
      COV_TMP       = 0.0 ; 
      COV_DIAG      = 0.0 ;

      // get cDisp for this filter; avoid repeating SALT2colorDisp calls
      if ( cDisp[ifilt_row] < -1.0 ) {
	ifilt         = IFILTMAP_SEDMODEL[ifilt_row] ;
	meanlam_obs   = FILTER_SEDMODEL[ifilt].mean ;  // mean lambda
	meanlam_rest  = meanlam_obs * invZ1 ; 
	cDisp[ifilt_row] = SALT2colorDisp(meanlam_rest,fnam);    
      }

      // set covariances only for same passband.
      if ( ifilt_col == ifilt_row ) 
	{ COV_TMP = FAC * pow(cDisp[ifilt_row],2.0);  }

      // check for local dump option
      LDMP  = (COV_TMP != 0.0 || ISDIAG) && 
	(fabs(Trest_row) < -1.0 && fabs(Trest_col) < -1.0); 

      if ( LDMP ) 
	{ printf(" xxx ############ COV_MODEL DUMP ############### \n"); }

      // diagonal-only term
      if ( ISDIAG ) {
	Tobs          = *(epobsList+irow) ;  
	Trest         = Tobs * invZ1 ;
	ifilt         = IFILTMAP_SEDMODEL[ifilt_row] ;
	cfilt         = FILTER_SEDMODEL[ifilt].name ;
	meanlam_obs   = FILTER_SEDMODEL[ifilt].mean ;  // mean lambda
	meanlam_rest  = meanlam_obs * invZ1 ; 


	// make sure that Trest is within the map range
	if ( Trest > SALT2_ERRMAP[0].DAYMAX ) 
	  { Trest_tmp = SALT2_ERRMAP[0].DAYMAX ; }
	else if ( Trest < SALT2_ERRMAP[0].DAYMIN ) 
	  { Trest_tmp = SALT2_ERRMAP[0].DAYMIN ; }
	else
	  { Trest_tmp = Trest ; }

	Trest = Trest_tmp ;
	Tobs  = Trest * z1 ;

	INTEG_zSED_SALT2(0,ifilt_row,z,Tobs, parList_SN, parList_HOST, // (I)
			 &Finteg, &Finteg_errPar, FspecDum); // returned

	magerr = SALT2magerr(Trest, meanlam_rest, z, x1, x2, 
			     Finteg_errPar, LDMP );
	COV_DIAG = magerr*magerr ;
	COV_TMP = COV_DIAG ;
      }
      
      covar[icovar] = COV_TMP ;  // load output array  
      icovar++ ;                   // increment local pointer

      if ( LDMP && COV_TMP != 0.0  ) {

	ifilt         = IFILTMAP_SEDMODEL[ifilt_row] ;
	cfilt         = FILTER_SEDMODEL[ifilt].name ;
	sprintf(cdum0,"%s:Tobs=%7.3f", cfilt, Tobs);

	ifilt         = IFILTMAP_SEDMODEL[ifilt_col] ;
	cfilt         = FILTER_SEDMODEL[ifilt].name ;
	sprintf(cdum1,"%s:Tobs=%7.3f", cfilt, Tobs);
 
	printf(" xxx COV_MAGERR[ %s , %s ] = %le \n", cdum0,cdum1, COV_TMP );

	if ( ISDIAG ) {
	  printf(" xxx ----------------- \n");
	  printf(" xxx COV_DIAGON[ %s , %s ] = %le \n", cdum0,cdum1, COV_DIAG);
	  printf(" xxx meanlam_rest = %f  z=%f  x1=%f  errPar=%f \n",
		 meanlam_rest, z, x1, Finteg_errPar );
	  printf(" xxx ----------------- \n");
	}

	fflush(stdout);
      }

    } // icol
  } //  irow

  return SUCCESS ; 

} // end of gencovar_SALT2

// ***********************************************
double SALT2colorDisp(double lam, char *callFun) {

  // Mar 2011
  // Return color dispersion for input rest-wavelength "lam".
  // Since this function INTERPOLATES, and does NOT extrapolate,
  // this function aborts if 'lam' is outside the valid
  // rest-lambda range. Make sure to check that 'lam' is valid
  // before calling this function.
  //
  // Jan 28 2020:
  //  if UV extrap is used, extrapolate instead of aborting
  //
  int imap, NLAM ;
  double cDisp, LAMMIN, LAMMAX ;
  double *mapLam, *mapDisp ;
  char fnam[] = "SALT2colorDisp" ;

  // ------------ BEGIN --------------

  // strip off goodies into local variables
  imap    = INDEX_SALT2_ERRMAP.COLORDISP ;
  NLAM    = SALT2_ERRMAP[imap].NLAM ;
  LAMMIN  = SALT2_ERRMAP[imap].LAMMIN ;
  LAMMAX  = SALT2_ERRMAP[imap].LAMMAX ;
  mapLam  = SALT2_ERRMAP[imap].LAM ;
  mapDisp = SALT2_ERRMAP[imap].VALUE ;

  if ( NLAM <= 0 ) { cDisp = 0.0 ; return cDisp ; }

  
  if ( INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX > 0.0 && lam < LAMMIN ) 
    { cDisp = mapDisp[0]; return(cDisp);  }

  // first some sanity checks
  if ( lam < LAMMIN || lam > LAMMAX ) {  
    sprintf(c1err,"lam=%f outside lookup range (called from %s)", 
	    lam, callFun );
    sprintf(c2err,"Valid range is %7.1f to %7.1f A ", LAMMIN, LAMMAX);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( NLAM <= 1 ) {
    sprintf(c1err,"Cannot do map-lookup with %d lambda bins (callFun=%s).", 
	    NLAM, callFun);
    sprintf(c2err,"Check %s",  SALT2_ERRMAP_FILES[imap] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // use generic linear interpolator.
  cDisp = interp_1DFUN(OPT_INTERP_LINEAR, lam, 
		       NLAM, mapLam, mapDisp, "cDisp" );


  // return final result
  return cDisp ;
    

} // end of SALT2colorDisp


// ***********************************************
void colordump_SALT2(double lam, double c, char *cfilt) {

  double colorCor;
  char cCor[60];
  // -------- BEGIN --------

  colorCor = SALT2colorCor(lam,c);
  if ( fabs(colorCor) < 100 ) 
    { sprintf(cCor,"%7.3f", colorCor ); }
  else
    { sprintf(cCor,"%9.3le", colorCor ); }

  printf("\t ColorTerm[ lam=%5.0f (%s)  c=%4.1f ] = %s \n", 
	 lam, cfilt, c, cCor);

} // end of colordump_SALT2

// ===============================
void errorSummary_SALT2(void) {

  // summarize errors and CL(lambda) in a list vs. lambda.

  int NLAM, ilam, imap, LLAM ;
  double LAMLIST[100], lam, Trest, ERRMAP[MXERRMAP_SALT2] ;
  double var0, var1, covar01, errscale, S0fracErr, colorCor, c, colorDisp ;   
  char cCor[20];
  //  char fnam[] = "errorSummary_SALT2" ;

  // ---------- BEGIN --------

  NLAM = 0;

  // hard-wire list of lambda values to check
  NLAM++ ; LAMLIST[NLAM] = 2000.0 ; // A
  NLAM++ ; LAMLIST[NLAM] = 2500.0 ;
  NLAM++ ; LAMLIST[NLAM] = 3000.0 ;
  NLAM++ ; LAMLIST[NLAM] = U_WAVELENGTH ;
  NLAM++ ; LAMLIST[NLAM] = 3560.0 ;   // u band
  NLAM++ ; LAMLIST[NLAM] = 3900.0 ;
  NLAM++ ; LAMLIST[NLAM] = B_WAVELENGTH ;
  NLAM++ ; LAMLIST[NLAM] = 4720.0 ;  // g band
  NLAM++ ; LAMLIST[NLAM] = V_WAVELENGTH ;
  NLAM++ ; LAMLIST[NLAM] = 6185.0 ;  // r band
  NLAM++ ; LAMLIST[NLAM] = R_WAVELENGTH ;
  NLAM++ ; LAMLIST[NLAM] = 7500.0 ;  // i band
  NLAM++ ; LAMLIST[NLAM] = 8030.0 ;  // I band
  NLAM++ ; LAMLIST[NLAM] = 8500.0 ;
  NLAM++ ; LAMLIST[NLAM] = 9210.0 ;  // z
  NLAM++ ; LAMLIST[NLAM] = 9940.0 ;  // Y


  Trest = 0.0;
  c = 1.0; // color value

  printf("\n");
  printf("                       e^CL     peak     color  \n" );
  printf("            LAMBDA(A)  (c=1)    dS0/S0   disp   \n" );
  printf("  --------------------------------------------- \n" );

  for ( ilam = 1; ilam <= NLAM; ilam++ ) {

    lam = LAMLIST[ilam];

    // color correcction (note that this is not an error)
    colorCor = SALT2colorCor(lam,c);
    if ( fabs(colorCor) < 100. ) 
      { sprintf(cCor,"%7.3f", colorCor ); }
    else
      { sprintf(cCor,"%9.3le", colorCor ); }

    // fractional flux error with x1=0
    get_SALT2_ERRMAP ( Trest, lam, ERRMAP );
    var0       = ERRMAP[INDEX_SALT2_ERRMAP.VAR[0]] ;  // sigma(S0)/S0
    var1       = ERRMAP[INDEX_SALT2_ERRMAP.VAR[1]] ;  // sigma(S1)/S0
    covar01    = ERRMAP[INDEX_SALT2_ERRMAP.COVAR[0][1]] ;  // 
    errscale   = ERRMAP[INDEX_SALT2_ERRMAP.ERRSCALE] ;  // error fudge
    S0fracErr  = errscale * sqrt(var0);         // dF/F with x1=0
	   
    // color dispersion
    imap = INDEX_SALT2_ERRMAP.COLORDISP ;
    if ( lam >= SALT2_ERRMAP[imap].LAMMIN &&
	 lam <= SALT2_ERRMAP[imap].LAMMAX ) {
      
      colorDisp = interp_1DFUN(OPT_INTERP_LINEAR, lam
			    ,SALT2_ERRMAP[imap].NLAM
			    ,SALT2_ERRMAP[imap].LAM
			    ,SALT2_ERRMAP[imap].VALUE
			    ,"colorDispSummary" );      
    }
    else
      { colorDisp = 0.0 ; }


    LLAM = (int)lam ;
    printf("  LAMINFO:  %6d  %8s   %6.4f   %5.3f  \n", 
	   LLAM, cCor, S0fracErr, colorDisp );

  }  // end of ilam loop



} // end of errorSummary_SALT2

// ================================
void test_SALT2colorlaw1(void) {

#define NCTEST 5
  double c[NCTEST] ;
  double claw[NCTEST] ;
  double colorPar[20];
  double lambda;

  int i, irow ;

  // --------- BEGIN ------------

  // define parameters from Julien's test code (read_color_law.c)
  colorPar[0] = B_WAVELENGTH ;
  colorPar[1] = V_WAVELENGTH ;
  colorPar[2] = 3700.0 ;
  colorPar[3] = 8000.0 ;
  colorPar[4] = 4.0 ;      // nparams
  colorPar[5] = -1.77139; 
  colorPar[6] =  2.38305 ; 
  colorPar[7] = -1.16417 ; 
  colorPar[8] =  0.178494 ;

  c[0] = 0.2 ;
  c[1] = 0.4 ;
  c[2] = 0.6 ;
  c[3] = 0.8 ;
  c[4] = 1.0 ;

  irow = 0;

  for( lambda=2500; lambda<9000; lambda+=100) {

    for ( i=0; i < NCTEST ; i++ ) {
      claw[i] = SALT2colorlaw1(lambda,c[i], colorPar) ;
    }

    irow++ ;
    printf("SN: %4.4d  %6.1f  %f %f %f %f %f\n", 
	   irow, lambda, claw[0], claw[1], claw[2],claw[3], claw[4] );

  }
 
  debugexit("Done testing SALT2colorlaw1");
} // end of test_SALT2colorlaw1


// ========================================================================
// ============== SPECTROGRAPH FUNCTIONS (July 2016) ======================
// ========================================================================


// ==========================
void genSpec_SALT2(double *parList_SN, double *parList_HOST, double mwebv,
                   double z, double Tobs, 
		   double *GENFLUX_LIST, double *GENMAG_LIST ) {

  // July 2016
  // For input SALT2 params, return *GENFLUX_LIST and *GENMAG_LIST.
  // FLUXGEN units are arbitrary since snlc_sim program will add
  // fluctuations based on user-input SNR.
  // 
  // Inputs:
  //   x0,x1,c  : SALT2 params
  //   mwebv    : MW E(B-V)
  //   RV/AV_host : host extinction
  //   Tobs     : T - Tpeak, obs frame (scalar, not array)
  //
  // Output:
  //   *GENFLUX_LIST : flux vs. wavelength bin.
  //                   flux=0 outside SALT2 model range (no aborts).
  //
  //   *GENMAG_LIST :  mag
  //
  // Note that output arrays have length NBLAM (see below)
  //
  // Aug 31 2016: return of Trest is outside defined epoch range.
  //              Can extrapolate mags, but NOT spectra !
  //
  // Mar 29 2019: apply MAG_OFFSET to GENFLUX_LIST
  //
  // Mar 23 2021: call fill_TABLE_MWXT_SEDMODEL
  //
  // Sep 18 2023: check option to allow extrapolation past DAYMAX;
  //              see local EXTRAP_LATETIME bool.
  //
  // Nov 6 2023: fix bug in which MWXT_FRAC had been ignored.
  //             This bug was probably in the SALT3 training papers.
  //
  // Jan 2024: call fill_TABLE_MWXT_SEDMODEL() in case genmag_SALT2 isn't called
  //        because all bands are outside model wavelength range
  //
  // ------------------------------------------

  int    NBLAM      = SPECTROGRAPH_SEDMODEL.NBLAM_TOT ;
  double MAG_OFFSET = INPUT_SALT2_INFO.MAG_OFFSET ;  

  int ilam ;  
  bool   DEFINED;
  double DAYMIN_SALT2  = SALT2_TABLE.DAYMIN+0.01 ;
  double DAYMAX_SALT2  = SALT2_TABLE.DAYMAX-0.01 ;
  double DAYMIN_EXTRAP = INPUT_EXTRAP_LATETIME_Ia.DAYMIN;

  bool DO_EXTRAP_LOCAL = false; 
  bool EXTRAP_METHOD_MAG     = (EXTRAP_PHASE_METHOD == EXTRAP_PHASE_MAG );
  bool EXTRAP_METHOD_SEDFLUX = (EXTRAP_PHASE_METHOD == EXTRAP_PHASE_SEDFLUX );
  bool EXTRAP_METHOD_FLAM    = (EXTRAP_PHASE_METHOD == EXTRAP_PHASE_FLAM);

  double Tobs_SED = Tobs; // Tobs to fetch SED

  double Trest, Finteg, Finteg_errPar, MWXT_FRAC ;
  double FTMP, GENFLUX, ZP, MAG, LAM, LAMREST, z1, FSCALE_ZP;
  double FTMP_DAYMAX, MAG_DAYMAX ;
  double hc8 = (double)hc ;
  int    LDMP;
  char fnam[] = "genSpec_SALT2" ;

  // -------------- BEGIN --------------

  z1 = 1.0 + z;

  fill_TABLE_MWXT_SEDMODEL(MWXT_SEDMODEL.RV, mwebv); 

  // init entire spectum to zero.
  for(ilam=0; ilam < NBLAM; ilam++ ) { GENFLUX_LIST[ilam] = 0.0 ; }

  Trest = Tobs / (1.0+z) ;

  if ( Trest < DAYMIN_SALT2 ) { return ; }

  // bail for generic SEDFLUX-extrap that does not work on spectra/SED
  if ( Trest > DAYMAX_SALT2 && EXTRAP_METHOD_SEDFLUX ) { return; }

  // note that DAYMIN_EXTRAP ~ 45 days, which is inside the nominal
  // SALT2 DAYMAX ~ 50 days.
  if ( EXTRAP_METHOD_MAG && Trest > DAYMIN_EXTRAP ) { // legacy; should be deleted
      Tobs_SED        = DAYMIN_EXTRAP * (1.0+z); 
      DO_EXTRAP_LOCAL = true;
  } 
  
  // - - - - -  
      
  INTEG_zSED_SALT2(1, JFILT_SPECTROGRAPH, z, Tobs_SED, 
		   parList_SN, parList_HOST,
		   &Finteg, &Finteg_errPar,  GENFLUX_LIST ) ;

  FSCALE_ZP = pow(TEN,-0.4*MAG_OFFSET);

  // convert generated fluxes into mags
  for(ilam=0; ilam < NBLAM; ilam++ ) { 
    GENFLUX_LIST[ilam] *= FSCALE_ZP ;  // Mar 29 2019

    MWXT_FRAC  = SEDMODEL_TABLE_MWXT_FRAC[JFILT_SPECTROGRAPH][ilam] ;
    GENFLUX = GENFLUX_LIST[ilam] ;
    LAM     = SPECTROGRAPH_SEDMODEL.LAMAVG_LIST[ilam] ;
    ZP      = SPECTROGRAPH_SEDMODEL.ZP_LIST[ilam] ;
    FTMP    = (LAM/hc8) * GENFLUX * MWXT_FRAC;

    /*xxxx
    printf(" xxx %s: ilam=%d LAM=%.1f ZP=%.3f  GENFLUX=%le\n",
	   fnam, ilam, LAM, ZP, GENFLUX); fflush(stdout);
    */

    DEFINED = ( ZP > 0.0 && FTMP > 0.0 ) ;

    // xxxxxxx mark: this extrap should be deleted when FLAM method is well tested
    if ( DEFINED && DO_EXTRAP_LOCAL ) {
      MAG_DAYMAX  = -2.5*log10(FTMP) + ZP; 
      FTMP_DAYMAX = FTMP ;
      LAMREST     = LAM/(1.0+z);
      MAG  = genmag_extrap_latetime_Ia(MAG_DAYMAX, Trest, LAMREST);
      FTMP = pow(10.0,0.4*(ZP-MAG)) ;
      GENFLUX_LIST[ilam] *= ( FTMP/FTMP_DAYMAX );

      LDMP = fabs(Tobs-9999994.543) < 0.001 && fabs(LAM-8050)<10.0 ;
      if ( LDMP) {
        printf(" xxx ----------------- \n");
        printf(" xxx %s: Tobs=%.3f z=%.3f  LAMREST = %.1f  \n",
               fnam, Tobs, z, LAMREST);
        printf(" xxx %s:    FTMP(edge,extrap) = %le  %le \n",
               fnam, FTMP_DAYMAX, FTMP);
        printf(" xxx %s:    MAG(edge,extrap)= %.3f %.3f \n",
               fnam, MAG_DAYMAX, MAG );
        fflush(stdout);
      } 
      // xxxxxxxx end mark xxxxxxxxx

    }
    else if ( DEFINED ) {
      MAG = -2.5*log10(FTMP) + ZP; 
    }
    else { 
      MAG = MAG_UNDEFINED ;  // model undefined
    }

    GENMAG_LIST[ilam] = MAG ;

  } // end ilam loop over SPECTROGRAPH bins

  return ;

} // end genSpec_SALT2


// ======================================================
int getSpec_band_SALT2(int ifilt_obs, float Tobs_f, float z_f,
                       float x0_f, float x1_f, float c_f, float mwebv_f,
		       float *LAMLIST_f, float *FLUXLIST_f) {

  // Created Nov 2016
  // Special utility used by LC fit program:
  // Return spectrum in band 'ifilt_obs' with passed SALT2 params.
  // Spectrum is returned as LAMLIST_f and FLUXLIST_f.
  // Note that all function args are float, but local 
  // variables are double.

  int ifilt      = IFILTMAP_SEDMODEL[ifilt_obs] ;
  int NBLAM      = FILTER_SEDMODEL[ifilt].NLAM ;
  int MEMD   = NBLAM * sizeof(double);
  int ilam ;
  double LAMOBS, LAMREST, z1, Finteg, Finteg_errPar, Finteg_check, TRANS ;
  double RV_host=-9.0, AV_host=0.0, m_host = -9.0  ;

  double Tobs  = (double)Tobs_f ;
  double z     = (double)z_f ;
  double x0    = (double)x0_f ;
  double x1    = (double)x1_f ;
  double c     = (double)c_f ;
  double *FLUXLIST = (double*) malloc ( MEMD );
  double Trest = Tobs/(1.0 + z) ;

  double parList_SN[3]   = { x0, x1, c };
  double parList_HOST[3] = { RV_host, AV_host, m_host } ;

  char fnam[] = "getSpec_band_SALT2" ;

  // ------------- BEGIN ---------------

  if ( Trest <= SALT2_TABLE.DAYMIN ) { return(0); }
  if ( Trest >= SALT2_TABLE.DAYMAX ) { return(0); }

  INTEG_zSED_SALT2(1, ifilt_obs, z, Tobs,         // (I)
		   parList_SN, parList_HOST,      // (I)
		   &Finteg, &Finteg_errPar, FLUXLIST ) ; // (O)
  
  Finteg_check = 0.0 ;  z1=1.0+z ;
  for(ilam=0; ilam < NBLAM; ilam++ ) {
    
    get_LAMTRANS_SEDMODEL(ifilt, ilam, &LAMOBS, &TRANS);

    LAMLIST_f[ilam]  = (float)LAMOBS ;
    FLUXLIST_f[ilam] = (float)FLUXLIST[ilam];

    // check Finteg; FLUXLIST already includes LAMSTEP

    LAMREST = LAMOBS/z1 ;
    Finteg_check += ( TRANS * LAMREST * FLUXLIST[ilam] ); 
  }
  
  /*
  printf(" xxx Tobs=%5.1f z=%.3f ifiltobs=%2d: Ratio_Finteg=%.3f (%le)\n",
	 Tobs, z, ifilt_obs, Finteg_check/Finteg, Finteg );  
  */

  free(FLUXLIST);

  return(NBLAM);

} // end getSpec_band_SALT2

int getspec_band_salt2__(int *ifilt_obs, float *Tobs, float *z,
			 float *x0, float *x1, float *c, float *mwebv,
			 float *LAMLIST, float *FLUXLIST) {
  int NBLAM;
  NBLAM = getSpec_band_SALT2(*ifilt_obs, *Tobs, *z, 
			     *x0, *x1, *c, *mwebv, LAMLIST, FLUXLIST ) ;
  return(NBLAM);
} 


// ======================================================
//  COLOR LAW FUNCTIONS
//  (moved back here from sntools.c, Sep 2020)
// ======================================================


double SALT2colorlaw0(double lam_rest, double c, double *colorPar ) {

  // Jul 2, 2010
  // Returns 10^[.4*c*C(lambda)] as defined in Guy 2007 SALT2 paper.
  //
  // *colorPar is an array of five parameters:
  // LAMBDA_B, LAMBDA_V, color-offset and two polynomial parameters
  // This code is moved from genamg_SALT2.c to allow more access.
  //
  // Aug 2, 2010: remove C_OFF parameter (previously 3rd colorPar arg)


  // define local args for *colorPar inputs
  double LAM_B, LAM_V, COR0, COR1 ;

  // local args
  double arg, lr, lr2, lr3, numerator, denominator, CLAM  ;

  char fnam[] = "SALT2colorlaw0" ;

  // -------- BEGIN ---------

  // strip off color law parameters
  LAM_B = *(colorPar+0); // mean lambda of B filter
  LAM_V = *(colorPar+1); // mean labmda of V filter
  COR0  = *(colorPar+2); // 1st fitted poly param from training
  COR1  = *(colorPar+3); // 2nd "    "

  // --------------------------------------
  // make a few sanity checks on passed parameters

  sprintf(c2err,"Check colorlaw parameters");

  if ( LAM_B < 4000 || LAM_B > 4500 ) {
    sprintf(c1err, "insane LAM_B = %6.0f", LAM_B );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( LAM_V < 5000 || LAM_V >6000 ) {
    sprintf(c1err, "insane LAM_V = %6.0f", LAM_V );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( COR0 < -0.4 || COR0 > -0.1 ) {
    sprintf(c1err, "insane COR0 = %6.3f", COR0 );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  if ( COR1 < -0.1 || COR1 > +0.1 ) {
    sprintf(c1err, "insane COR1 = %6.3f", COR1 );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // --------------------------------------
  lr = (lam_rest - LAM_B)/( LAM_V - LAM_B );
  lr2 = lr * lr;
  lr3 = lr * lr2 ;

  numerator   = lr + lr2*COR0 + lr3*COR1 ;
  denominator = 1.0 + COR0 + COR1 ;
  CLAM        = numerator/denominator ;
  arg         = 0.4 * c * CLAM ;

  return  pow(10.0,arg);


} // end of  SALT2colorlaw0


// *************************************************************
double SALT2colorlaw1(double lambda, double c, double *colorPar ) {

  // Created Jul 31, 2010 by R.Kessler and J.Guy
  // Returns 10^[.4*c*C(lambda)] as defined in 
  // Guy 2010 SALT2/SNLS3 paper.
  //
  // *colorPar is an array of parameters as defined below.
  // Code is from Julien, with slight modifications for
  // SNANA compatibility.
  // 
  // Sep 2020: use checkval_D to check values.

  double LAM_B, LAM_V, LAM_MIN, LAM_MAX, XN, params[10] ;
  int nparams, i  ;
  double constant = log(10.0)/2.5 ;
  double alpha    = 1.0 ;
  double val      = 0.0 ;
  double rl, rlmin, rlmax, tmp ;
  char fnam[] = "SALT2colorlaw1" ;

  // --------------- BEGIN ------------

  sprintf(c2err,"Check colorlaw parameters");

  // strip of colorPar values into local variables.
  LAM_B   = colorPar[ICLPAR_REFLAM_CL0] ;   // mean lambda of B filter
  LAM_V   = colorPar[ICLPAR_REFLAM_CL1] ;   // mean labmda of V filter
  LAM_MIN = colorPar[ICLPAR_LAM_MIN] ;   // special extrap function below this
  LAM_MAX = colorPar[ICLPAR_LAM_MAX] ;   // idem for upper wavelength
  XN      = colorPar[ICLPAR_NPAR_POLY] ;   // Number of parameters for function

  nparams = (int)XN ;        // number of parameters to follow
  for ( i=0; i < nparams; i++ ) {
    tmp       = colorPar[5+i]; 
    params[i] = tmp;
    if ( fabs(tmp) > 10.0 ) {
      sprintf(c1err, "insane params[%d] = %f", i, tmp );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  // --------------------------------------
  // make a few sanity checks on passed parameters

  checkval_D("CL1-LAM_B",   1, &LAM_B,   4000.0,  4500.0 );
  checkval_D("CL1-LAM_V",   1, &LAM_V,   5000.0,  6000.0 );
  checkval_D("CL1-LAM_MIN", 1, &LAM_MIN, 1000.0,  6000.0 );
  checkval_D("CL1-LAM_MAX", 1, &LAM_MAX, 6000.0, 18000.0 );

  // ------------------------------

  for(i=0; i < nparams; i++)
    { alpha = alpha - params[i]; }

  // compute reduced wavelengths
  rl    = (lambda  - LAM_B) / ( LAM_V - LAM_B );
  rlmin = (LAM_MIN - LAM_B) / ( LAM_V - LAM_B );
  rlmax = (LAM_MAX - LAM_B) / ( LAM_V - LAM_B );

  if(lambda >= LAM_MIN  && lambda <= LAM_MAX ) {
    val = SALT2colorfun_pol(rl,nparams,params,alpha) ;
  }
  else if(lambda < LAM_MIN ) {
    // extrapolate to UV
      double Pmin  = SALT2colorfun_pol(rlmin,nparams,params,alpha);
      double dPmin = SALT2colorfun_dpol(rlmin,nparams,params,alpha);
      val =  Pmin + dPmin* (rl-rlmin);
  } 
  else {
    // extrapolate to IR
    double Pmax  = SALT2colorfun_pol(rlmax,nparams,params,alpha);
    double dPmax = SALT2colorfun_dpol(rlmax,nparams,params,alpha);
    val =  Pmax + dPmax* (rl-rlmax);
  }

  double CL = exp(c*constant*val);
  return CL ;

} // end of SALT2colorlaw1



double SALT2colorfun_dpol(const double rl, int nparams, 
			  const double *params, const double alpha) {
  double v = alpha;
  double rlp = rl;
  int i;
  for(i=0; i<nparams; i++) {
    v += (i+2)*params[i]*rlp;
    rlp *= rl; 
  }
  return v;
}

double SALT2colorfun_pol(const double rl, int nparams, 
			 const double *params, const double alpha) {
  double v   = alpha*rl;
  double rlp = rl*rl;  
  int i;
  for(i =0; i<nparams; i++) {
    v += params[i]*rlp; // v = alpha*rl + rl^2*( sum_i p_i*rl^i)
    rlp *= rl; // rl^(i+2)
  }
  return v;
} 

