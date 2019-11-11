/* =================================================

  March, 2011  R.Kessler

  Inlcude this file in the simulation to 

  - select a host-galaxy with ZTRUE ~ ZSN and randomly selected
    based on a user-defined weight map.

  - determine a SN mag-offset based on galaxy properties.

  - randomly pick SN position weighted by surface brightness.

  - for list of SN apertures (r = 2x sig_PSF), determine galaxy 
    light-fraction inside each aperture and corresponding
    galaxy mag. At each simulated epoch, these GALMAGs are 
    interpolated according to the PSF to determine the
    extra Poisson noise from the host galaxy.

  Particular attention is given to minimize CPU time by
  pre-computing
  - Sersic integrals for each Sersic indexs
  - 2d Gaussian integrals contained in aperture
  - grid of sin(TH) and cos(TH) 

  CPU time on sdssdp62 with 100 radial bins and 36 theta bins :
     3 msec per host (DES  with 1 Sersic profile)
     5 msec per host (LSST with 2 Sersic profiles)
 

  The user-input control is via the sim-input structure
  INPUTS.HOSTLIB_MSKOPT, and the bits are defined in the
  parameters HOSTLIB_MSKOPT_XXX below. To see a list of 
  mask-options,
    grep MSKOPT $SNANA_DIR/src/sntools_host.h 
  

  This HOSTLIB package is designed specifically for the 
  SNANA simulation, and is not intended as a general 
  package for other software environments. 


  HOSTLIB size: there is no limit to the number of HOSTLIB
  entries since memory is allocated incrementally (in chunks
  of MALLOCSIZE_HOSTLIB) as the HOSTLIB is read. The limitation 
  is therefore set by the memory available to your processor.


 TODO_LIST: 
  - photoZ clipping with HOSTLIB_GENRANGE_NSIGZ

  - set DEBUG_WGTFLUX=0 when LSST hostlib has Sersic weights

  - find image of rotated galaxy to check treatment of a_rot.

  - increase Gal noise due to template subtraction 
     (JPAS/Xavier request)


          HISTORY
  ------------------------

  Jul 5, 2011
    replace  get_GALWGT_HOSTLIB() with utility function
    interp_GRIDMAP(). Replace most WGTMAP init with
    call to init_interp_GRIDMAP().


 Sep 22, 2011:
   Select Host galaxies based on RA and DECL;
   see INPUTS.HOSTLIB_GENRANGE_RA[2] and _DECL[2].

 May 22, 2012: for gen-filters that are not defined in the HOSTLIB,
               set MAGOBS=99 to avoid seg fault (see GEN_SNHOST_GALMAG)

 Sep 09, 2012: add lots of fflush(stdout) calls after printf().

 Sep 14, 2012: implement new FIXRAN_RADIUS and FIXRAN_PHI options

 Dec 17, 2012: in GEN_SNHOST_GALMAG(int IGAL) store total host mag
               in  SNHOSTGAL.GALMAG[ifilt_obs][i=0]

 Feb 12, 2013: alway read/store host mags so that they are
               written out, even if gal-noise option is not set.

 Feb 22, 2013: add include statements to compile this module instead
               of including it in the simulation.

 Feb 11 2014:  call init_OUTVAR_HOSTLIB()

 Jan 29 2015: allow RA or RA_GAL or RA_HOST, and same for DEC.
              --> distinguish RA_GAL and RA_SN in SIMGEN-DUMP file
                  and in SNANA tables.

  Feb 5 2015: int GALID -> long long GALID in many places.

  Mar 4 2015: fix aweful index bug set FGAL_TOT in
              double get_GALFLUX_HOSTLIB(..)

  Jun 24 2015: fix check on R/Re at F/FTOT=0.5 to interpolate R/Re 
               instead of using closest integral bin to 0.5.


  July 2015: allow re-using host if MJD is separated from previous use.
             See  USEHOST_GALID(int IGAL)
             See INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL ;

  Aug 18 2015: new routine init_HOSTLIB_ZPHOTEFF to read eff vs. ZTRUE
               to find host photoz. Implemented in  GEN_SNHOST_ZPHOT().

  Jan 14 2017: include "sntools_output.h"
    
  Dec 29 2017: speed up init, redGal_HOSTLIB and zptr_HOSTLIB

  Dec 2018: ugly refactor to change all 1-N loops to 0 to N-1 loops.

  Mar 20 2019: 
    +read VARNAMES list from zHOST effic map and append STOREPAR_LIST.
     See copy_VARNAMES_zHOST_to_HOSTLIB_STOREPAR()

  Apr 06 2019: 
    + fix aweful bug when SN coords are transferred to host;
      CMB redshift was incorrectly computed from original SN
      coords instead of from HOST coords. The crucial fix is
      getting the right distance modulus since the smeared
      redshift is computed correctly later in gen_zsmear().

    + account for VPEC when SN coords are transferred to host.

=========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_trigger.h"
#include "snlc_sim.h" 
#include "sntools_host.h"
#include "sntools_output.h"
#include "sntools_genGauss_asym.h"
#include "MWgaldust.h"
#include "sntools_spectrograph.h"
#include "genmag_SEDtools.h"

// ==================================
void INIT_HOSTLIB(void) {

  //  March 2011, R.Kessler
  //  read and  initialize host-galaxy library (HOSTLIB)
  //  for SNANA simulation.

  int USE ;
  FILE *fp_hostlib ;
  char fnam[] = "INIT_HOSTLIB" ;

  // ---------------- BEGIN -------------

  USE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USE) ;
  if ( USE == 0 ) 
    {  checkAbort_noHOSTLIB() ;    return ;  }
  else
    { checkAbort_HOSTLIB(); }

  TIME_INIT_HOSTLIB[0]  = time(NULL);
  print_banner("INIT_HOSTLIB(): Read host-galaxy library.");

  // check for spectral templates to determin host spectrum
  read_specbasis_HOSTLIB();

  // set inital values for HOSTLIB structure
  initvar_HOSTLIB();

  // check to read external WEIGHT-MAP instead of the HOSTLIB WEIGHT-MAP
  read_wgtmap_HOSTLIB();

  // open hostlib file and  return file pointer
  open_HOSTLIB(&fp_hostlib);

  // read header info : NVAR, VARNAMES ...
  read_head_HOSTLIB(fp_hostlib);

  // check for match among spec templates and hostlib varnames (Jun 2019)
  match_specbasis_HOSTVAR();

  //  debugexit(fnam); // xxx REMOVE

  // read GAL: keys
  read_gal_HOSTLIB(fp_hostlib);

  // summarize SNPARams that were/weren't found
  summary_snpar_HOSTLIB();

  // close HOSTLIB file

  if ( HOSTLIB.GZIPFLAG ) 
    { pclose(fp_hostlib); } // close gzip file
  else
    { fclose(fp_hostlib); } // close normal file stream

  // sort HOSTLIB entries by redshift
  sortz_HOSTLIB();

  // abort if any GALID+ZTRUE pair appears more than once
  check_duplicate_GALID();
  
  // set redshift pointers for faster lookup
  zptr_HOSTLIB();

  // setup optional wgt-map grid
  init_HOSTLIB_WGTMAP();

  // read optional EFF(zPHOT) vs. ZTRUE (Aug 2015)
  init_HOSTLIB_ZPHOTEFF();

  // prepare integral tables for Sersic profile(s).
  init_Sersic_VARNAMES();
  init_Sersic_HOSTLIB();

  // init options for re-using same host
  init_SAMEHOST();

  // init parameters and Gauss2d integrals for galaxy aperture mag
  init_GALMAG_HOSTLIB();

  // prepare comments for README file and/or screen dump
  readme_HOSTLIB();


  TIME_INIT_HOSTLIB[1]  = time(NULL);
  double dT = (TIME_INIT_HOSTLIB[1]-TIME_INIT_HOSTLIB[0]);
  printf("\t HOSTLIB Init time: %.2f seconds \n", dT );
  fflush(stdout);

  //  debugexit(fnam); // xxxx REMOVE
  return ;

} // end of INIT_HOSTLIB


// ====================================
void initvar_HOSTLIB(void) {

  // one-time init of variables used for HOSTLIB

  int ivar, j, igal, ifilt  ;
  //  char fnam[] = "initvar_HOSTLIB" ;

  // ----------- BEGIN -------------

  HOSTLIB.SORTFLAG     = 0 ;
  HOSTLIB.NGAL_READ    = 0 ;
  HOSTLIB.NGAL_STORE   = 0 ;
  HOSTLIB.NVAR_ALL     = 0 ; 
  HOSTLIB.NVAR_STORE   = 0 ; 
  HOSTLIB.MALLOCSIZE_D = 0 ;
  HOSTLIB.MALLOCSIZE_I = 0 ;

  HOSTLIB.NVAR_SNPAR   = 0 ;
  HOSTLIB.VARSTRING_SNPAR[0] = 0 ;

  for ( ivar=0; ivar < MXVAR_HOSTLIB; ivar++ )  { 
    HOSTLIB.IVAR_WGTMAP[ivar] = -9 ; 
    HOSTLIB.IVAR_STORE[ivar]  = -9 ; 
    HOSTLIB.VALMAX[ivar]      = -1.0E18 ;
    HOSTLIB.VALMIN[ivar]      = +1.0E18 ;
    sprintf(HOSTLIB.VARNAME_STORE[ivar],"%s", NULLSTRING );

    HOSTLIB.IS_SNPAR_OPTIONAL[ivar] = 0 ;
    HOSTLIB.IS_SNPAR_STORE[ivar]    = 0 ;
  }

  HOSTLIB.ZGAPMAX       = -9. ;
  HOSTLIB.Z_ATGAPMAX[0] = -9. ;
  HOSTLIB.Z_ATGAPMAX[1] = -9. ;
  HOSTLIB.ZGAPAVG       = 0.0 ;
  HOSTLIB.ZMAX          = -99.0 ; 
  HOSTLIB.ZMIN          = +99.0 ; 

  HOSTLIB.NLINE_COMMENT = 0;

  SERSIC_PROFILE.NFIX   = 0;
  SERSIC_PROFILE.NDEF   = 0 ; 

  SERSIC_TABLE.NBIN_reduced   = 0 ;
  SERSIC_TABLE.TABLEMEMORY    = 0 ;

  HOSTLIB.NFILT_MAGOBS = 0;
  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ )     {
    HOSTLIB.IVAR_MAGOBS[ifilt] = -9 ;
  }
  sprintf(HOSTLIB.filterList, "%s", "" );


  // -----------------------
  // init Sersic shape parameters a0, a1 ...a9 and b0, b1 ... b9
  for ( j=0; j < MXSERSIC_HOSTLIB ; j++ ) {       
    SERSIC_PROFILE.IVAR_a[j]  = -9 ;
    SERSIC_PROFILE.IVAR_b[j]  = -9 ;
    SERSIC_PROFILE.IVAR_w[j]  = -9 ;
    SERSIC_PROFILE.IVAR_n[j]  = -9 ;
    SERSIC_PROFILE.FIXn[j]    = -9.0 ; 
  }

  // init Sersic Table
  for ( j=0; j < NSERSIC_TABLE ; j++ ) {   
    SERSIC_TABLE.inv_n[j] = 999. ;
    SERSIC_TABLE.n[j]     = 999. ; 
    SERSIC_TABLE.bn[j]    = 999. ;
  }

  // ------------------------
  // set array of required and optional keys
  init_REQUIRED_HOSTVAR();

  // now the optional keys
  init_OPTIONAL_HOSTVAR();


  HOSTLIB_WGTMAP.GRIDMAP.NDIM = 0;
  HOSTLIB_WGTMAP.GRIDMAP.NFUN = 0;
  HOSTLIB_WGTMAP.GRIDMAP.NROW = 0;
  HOSTLIB_WGTMAP.WGTMAX    = 0.0 ;
  HOSTLIB_WGTMAP.ISTAT     = 0 ; // init to weight-map NOT read
  HOSTLIB_WGTMAP.NCHECKLIST = 0;
  for ( ivar=0; ivar < MXVAR_WGTMAP_HOSTLIB; ivar++ ) {  
    sprintf(HOSTLIB_WGTMAP.VARNAME[ivar], "%s", NULLSTRING );
  }

  for ( igal = 0; igal < MXCHECK_WGTMAP ; igal++ ) {
    HOSTLIB_WGTMAP.CHECKLIST_IGAL[igal] = -9 ;
  }

  return ;

} // end of initvar_HOSTLIB

// ==========================================
void init_OPTIONAL_HOSTVAR(void) {

  // Feb 2014 [extracted from initvar_HOSTLOB]
  // Define array of optional variables to read from hostlib.
  // --> will read them if they exist; otherwise ignore them
  //     and continue.
  //
  // Jan 30 2015: allow RA or RA_HOST, DEC or DEC_HOST

  int NVAR, j, ifilt, ifilt_obs ;

  char anam[12], bnam[12], wnam[12], nnam[12];
  char varName[40], *cptr ;
  char fnam[] = "init_OPTIONAL_HOSTVAR" ;

  // ----------- BEGIN ---------------

  NVAR = 0;

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_ZPHOT );

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_ZPHOT_ERR );

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_RA ) ;
  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_RA_HOST ) ;
  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_RA_GAL ) ; 

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_DEC );
  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_DEC_HOST ); 
  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_DEC_GAL ); 

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_LOGMASS );

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_LOGMASS_ERR );

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_FIELD ); 

  // allow Sersic shape parameters a0, a1 ...a9 and b0, b1 ... b9
  for ( j=0; j < MXSERSIC_HOSTLIB ; j++ ) {   
    
    Sersic_names(j, anam, bnam, wnam, nnam );

    cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;    NVAR++;  
    sprintf(cptr, "%s", anam );    // major half-light axis

    cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;    NVAR++; 
    sprintf(cptr, "%s", bnam );    // minor half-light axis

    cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;    NVAR++;  
    sprintf(cptr, "%s", wnam );    // wgt = Flux/flutTotal

    cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;   NVAR++;  
    sprintf(cptr, "%s", nnam );    // index 

  }


  cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ; NVAR++; 
  sprintf(cptr,"%s", HOSTLIB_VARNAME_ANGLE );

  // check for observer-frame mags '[filt]_obs' 
  // used to determine galaxy noise  to SN signal-flux
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    magkey_HOSTLIB(ifilt_obs,varName); // returns varname
    cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ; NVAR++;  
    sprintf(cptr,"%s", varName);
  }


  // check for optional use of SN params from HOLSTIB; e.g., c, x1, delta
  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) > 0 ) {

    int N_SNPAR = 0 ;
    char tmpName[10][40] ;

    sprintf(tmpName[N_SNPAR],"%s", GENLC.SHAPEPAR_NAME);    N_SNPAR++ ; 
    sprintf(tmpName[N_SNPAR],"%s", GENLC.COLORPAR_NAME);    N_SNPAR++ ; 
    sprintf(tmpName[N_SNPAR],"%s", GENLC.SHAPEPAR2_NAME);   N_SNPAR++ ; 
    sprintf(tmpName[N_SNPAR],"%s", GENLC.COLORPAR2_NAME);   N_SNPAR++ ; 
    sprintf(tmpName[N_SNPAR],"%s", HOSTLIB_VARNAME_SNMAGSHIFT );  N_SNPAR++ ; 
    HOSTLIB.NVAR_SNPAR   = N_SNPAR ;

    for ( j=0; j < N_SNPAR; j++ ) {
      strcat(HOSTLIB.VARSTRING_SNPAR,tmpName[j]);
      strcat(HOSTLIB.VARSTRING_SNPAR," ");
      cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
      HOSTLIB.IS_SNPAR_OPTIONAL[NVAR] = 1 ;
      HOSTLIB.FOUND_SNPAR_OPTIONAL[NVAR] = 0 ; // init FOUND logicals to 0
      sprintf(cptr, "%s", tmpName[j] );
      NVAR++;  
    }   
  }


  if ( NVAR >= MXVAR_HOSTLIB ) {
    sprintf(c1err, "NVAR_OPTIONAL=%d exceeds bound", NVAR ) ;
    sprintf(c2err, "of MXVAR_HOSTLIB=%d", MXVAR_HOSTLIB ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  HOSTLIB.NVAR_OPTIONAL = NVAR ; // load global NVAR_OPTIONAL

} // end of init_OPTIONAL_HOSTVAR


// ==========================================
void init_REQUIRED_HOSTVAR(void) {

  int NVAR,  LOAD, i ;
  char *cptr, *varName ;
  //  char fnam[] = "init_REQUIRED_HOSTVAR" ;

  // ----------- BEGIN -----------

  NVAR = 0;

  cptr = HOSTLIB.VARNAME_REQUIRED[NVAR] ;  NVAR++; 
  sprintf(cptr, "%s", HOSTLIB_VARNAME_GALID );
  LOAD = load_VARNAME_STORE(cptr) ;

  cptr = HOSTLIB.VARNAME_REQUIRED[NVAR] ;  NVAR++; 
  sprintf(cptr, "%s", HOSTLIB_VARNAME_ZTRUE );
  LOAD = load_VARNAME_STORE(cptr) ;


  // check for required specTemplate coefficients
  for(i = 0; i < HOSTSPEC.NSPECBASIS; i++ ) {
    varName = HOSTSPEC.VARNAME_SPECBASIS[i];
    cptr = HOSTLIB.VARNAME_REQUIRED[NVAR] ;  NVAR++;
    sprintf(cptr, "%s%s", PREFIX_SPECBASIS_HOSTLIB, varName );
    LOAD = load_VARNAME_STORE(cptr) ;
  }

  HOSTLIB.NVAR_REQUIRED = NVAR ;

  // zHOST gets parsed later with trigger maps, but here we check
  // for HOSTLIB dependence to ensure needed HOSTLIB columns are read.
  // Beware of spaghetti code.
  copy_VARNAMES_zHOST_to_HOSTLIB_STOREPAR(); 

  // check user-specified HOSTLIB variables to store in data files.
  // --> add to the VARNAME_REQUIRED list if not already there.
  // See sim-input key HOSTLIB_OUTVARLIST.
  init_OUTVAR_HOSTLIB();

  return ;

} // end of init_REQUIRED_HOSTVAR


// ==============================================
void  checkAbort_HOSTLIB(void) {

  // Apr 8 2019

  int USE_GALCOORDS  = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC);
  int USE_MWEBV_FILE = (INPUTS.OPT_MWEBV == OPT_MWEBV_FILE);  
  char fnam[] = "checkAbort_HOSTLIB" ;

  // ---------------- BEGIN ----------------


  if ( USE_GALCOORDS && USE_MWEBV_FILE )  {
    sprintf ( c1err, "Cannot transfer SN coords to Gal coords");
    sprintf ( c2err, "without specifying MWEBV model (OPT_MWEBV>1)");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  return ;
} // end checkAbort_HOSTLIB

// ==============================================
void  checkAbort_noHOSTLIB(void) {

  // abort if user picked options requiring a HOSTLIB

  char fnam[] = "checkAbort_noHOSTLIB" ;

  if ( INPUTS.USE_HOSTLIB_GENZPHOT ) {
    sprintf ( c1err, "Cannot use HOSTLIB_GENZPHOT_XXX options " );
    sprintf ( c2err," without defining HOSTLIB_FILE" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  double GAMMA_GRID_MIN = INPUTS.BIASCOR_SALT2GAMMA_GRID[0]; 
  double GAMMA_GRID_MAX = INPUTS.BIASCOR_SALT2GAMMA_GRID[1]; 
  if ( GAMMA_GRID_MAX > GAMMA_GRID_MIN ) {
    sprintf ( c1err, "Cannot use sim-input BIASCOR_SALT2GAMMA_GRID");
    sprintf ( c2err," without defining HOSTLIB_FILE" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return ;

} // end checkAbort_noHOSTLIB


// ==============================================
int load_VARNAME_STORE(char *varName) {

  // Feb 2014
  // load HOSTLIB.VARNAME_STORE array with input *varName 
  // if NOT already on the list.
  // Return +1 if varName is loaded here.
  // return -1 if this varName was already loaded.

  int N ;
  if ( IVAR_HOSTLIB(varName,0) < 0 ) {
    N = HOSTLIB.NVAR_STORE ;
    sprintf( HOSTLIB.VARNAME_STORE[N],"%s", varName);
    HOSTLIB.NVAR_STORE++ ;
    return 1 ;
  }
  else  { 
    return -1;
  }


} // end of load_VARNAME_STORE

// ==============================================
void copy_VARNAMES_zHOST_to_HOSTLIB_STOREPAR(void) {

  // Mar 20 2019
  // copy variables from zHOST efficiency map to 
  // INPUTS.HOSTLIB_STOREPAR_LIST --> ensure that
  // all of the HOSTLIB-zHOST parameters are read
  // from the HOSTLIB.

  FILE *fp ;
  //  char fnam[] = "copy_VARNAMES_zHOST_to_HOSTLIB_STOREPAR" ;

  // -------------- BEGIN ------------

  // open zHOST file
  fp = open_zHOST_FILE(-1);
  if ( fp != NULL ) 
    { read_VARNAMES_zHOST(fp); fclose(fp); }
  else 
    { return ; }


  // if we get here, append zHOST varNames to HOSTLIB_STOREPAR_LIST.
  int ivar; 
  int  NVAR = SEARCHEFF_zHOST[0].NVAR ;
  char *plist = INPUTS.HOSTLIB_STOREPAR_LIST ;
  char *varName;

  for(ivar=0; ivar < NVAR; ivar++ ) {
    varName = SEARCHEFF_zHOST[0].VARNAMES_HOSTLIB[ivar];
    if ( strlen(plist) > 0 ) { strcat(plist,","); }
    strcat(plist,varName);

  } // end ivar

  return ;

} // end copy_VARNAMES_zHOST_to_HOSTLIB_STOREPAR

// ====================================================
void  init_OUTVAR_HOSTLIB(void) {

  // Feb 2014
  // strip variables from comma-separated INPUTS.HOSTLIB_STOREPAR_LIST,
  // and store them in VARNAME_STORE and REQUIRED lists. 
  // Be careful NOT to store the same variable twice in case
  // user requests a variable that is already required or
  // in the WGTMAP.

  int   NVAR_STOREPAR, NVAR_OUT, NVAR_REQ, LOAD, ivar, ivar2, ISDUPL;
  char  VARLIST_ALL[MXPATHLEN], VARLIST_LOAD[MXPATHLEN];
  char  varName[60], *varName2;
  //  char  fnam[] = "init_OUTVAR_HOSTLIB"     ;

  // ---------------- BEGIN ----------------

  HOSTLIB_OUTVAR_EXTRA.NOUT = NVAR_OUT = 0 ;

  // note that VARLIST_ALL may contain duplicates
  // VARLIST_LOAD contains only unique variables, and is for print only
  sprintf(VARLIST_ALL, "%s", INPUTS.HOSTLIB_STOREPAR_LIST) ;
  VARLIST_LOAD[0] = 0 ; 


  // bail if nothing is specified.
  if ( strlen(VARLIST_ALL) == 0 ) { return ; }

  // ----------------------------------

  NVAR_REQ = HOSTLIB.NVAR_REQUIRED ;

  // split VARLIST_ALL into individual var names
  NVAR_STOREPAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,VARLIST_ALL);

  for(ivar=0; ivar < NVAR_STOREPAR; ivar++ ) {

    get_PARSE_WORD(0,ivar,varName); // return varName

    // skip duplicates, which could come from zHOST effic map,
    // or user mistake setting up STOREPAR_LIST
    ISDUPL = 0 ;
    for(ivar2=0; ivar2 < ivar; ivar2++ ) {
      varName2 = HOSTLIB_OUTVAR_EXTRA.NAME[ivar2];
      if ( strcmp(varName,varName2) == 0 ) { ISDUPL=1; }
    }
    if  ( ISDUPL ) { continue ; }

    // load global VARNAME_STORE array;
    // LOAD=+1 if not already loaded; LOAD=-1 if already loaded.
    LOAD = load_VARNAME_STORE(varName);

    // load REQUIRED list if VARNAME_STORE (above) was also loaded.
    if ( LOAD > 0 ) { 
      NVAR_REQ++ ;
      sprintf(HOSTLIB.VARNAME_REQUIRED[NVAR_REQ], "%s", varName);
      //printf("\t xxx add '%s' to  VARNAME_REQUIRED list \n", varName);
    }


    // always store variable in OUTVAR list, along with IVAR
    sprintf(HOSTLIB_OUTVAR_EXTRA.NAME[NVAR_OUT], "%s", varName );
    HOSTLIB_OUTVAR_EXTRA.IVAR_STORE[NVAR_OUT] = IVAR_HOSTLIB(varName,1);

    // set USED_IN_WGTMAP flat to zero; fill this later in read_head_HOSTLIB.
    HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[NVAR_OUT] = 0 ;

    NVAR_OUT++ ;  

    if(NVAR_OUT > 1 ) { strcat(VARLIST_LOAD,"," ); }
    strcat(VARLIST_LOAD,varName); // for stdout message.

  } // end ivar loop


  // update global counters
  HOSTLIB_OUTVAR_EXTRA.NOUT = NVAR_OUT ;
  HOSTLIB.NVAR_REQUIRED     = NVAR_REQ ;

  // - - - - - -
  printf("\t Load '%s' for SNTABLE storage. \n", VARLIST_LOAD);
  fflush(stdout);

  return ;

} // end of  init_OUTVAR_HOSTLIB



// =========================================
void Sersic_names(int j, char *a, char *b, char *w, char *n) {

  // Return hostlib VARNAME  for type = 'a', 'b', 'w'  and 'n'
  // and  integer label j.

  char suffix[12] = "_Sersic" ;
  // -------------- BEGIN -----------
  sprintf(a, "a%d%s", j, suffix );
  sprintf(b, "b%d%s", j, suffix );
  sprintf(w, "w%d%s", j, suffix );
  sprintf(n, "n%d%s", j, suffix );
} // end of Sersic  name


// ==================================
void open_HOSTLIB(FILE **fp) {

  // Dec 29 2017: use snana_openTextFile utility to allow gzipped library.

  char libname_full[MXPATHLEN] ;
  char PATH_DEFAULT[MXPATHLEN] ;
  char fnam[] = "open_HOSTLIB" ;

  // ----------- BEGIN ----------

  sprintf(PATH_DEFAULT, "%s/simlib", PATH_SNDATA_ROOT );
  *fp = snana_openTextFile(0,PATH_DEFAULT, INPUTS.HOSTLIB_FILE,
			   libname_full, &HOSTLIB.GZIPFLAG );  // <== returned


  if ( *fp == NULL ) {
    sprintf ( c1err, "Cannot open file :" );
    sprintf ( c2err," '%s' ", libname_full );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  sprintf(HOSTLIB.FILENAME , "%s", libname_full );
  printf("\t Reading %s \n", libname_full );
  fflush(stdout);

}  // end of open_HOSTLIB

// ====================================
void  read_wgtmap_HOSTLIB(void) {

  // Function to read OPTIONAL weight-map to over-ride
  // weight map in the HOSTLIB. If the weight map is read
  // here, then the corresponding weight-map in the HOSTLIB
  // will be ignored. Note that this function must be called
  // before read_head_HOSTLIB().

  char 
    *ptrFile
    ,c_get[40]
    ,fnam[] = "read_wgtmap_HOSTLIB" 
    ;

  FILE *fp ;

  // ------------- BEGIN --------------

  HOSTLIB_WGTMAP.ISTAT = 0 ;

  ptrFile = INPUTS.HOSTLIB_WGTMAP_FILE ;
  if ( IGNOREFILE(ptrFile) )  { return ; }

  if ( (fp = fopen(ptrFile, "rt")) == NULL ) {
    sprintf(c1err,"%s", "Could not find supplemental WGTMAP file:" );
    sprintf(c2err,"%s", ptrFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // if we get here, open and read WGTMAP file.

  printf("\t Read WEIGHT-MAP from supplemental file:\n\t   %s\n", ptrFile );
  fflush(stdout);

  while( (fscanf(fp, "%s", c_get)) != EOF) 
    { parse_WGTMAP_HOSTLIB(fp,c_get);  }

  fclose(fp);

  HOSTLIB_WGTMAP.ISTAT = 1 ;

  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SNMAGSHIFT )
    { printf("\t Implement SNMAGSHIFT in WGTMAP \n"); fflush(stdout); }
  else
    { printf("\t Ignore SNMAGSHIFT in WGTMAP \n"); fflush(stdout); }


} // end of read_wgtmap_HOSTLIB


// ====================================
void parse_WGTMAP_HOSTLIB(FILE *fp, char *string) {

  // Parse WGTMAP variables from file *fp.
  // *string is the current string value to check
  // if this is one of the WGTMAP keys.
  //
  // Mar 14 2019: refactor to use read_GRIDMAP().
  // Apr 12 2019: return of string != VARNAMES_WGTMAP

  int  IDMAP = IDGRIDMAP_HOSTLIB_WGTMAP ;
  long long GALID ;
  int FOUND_VARNAMES;
  int NVAR_WGTMAP, IVAR_STORE, NDIM, NFUN, ivar, N ;

  char LINE[100], *VARNAME ;
  char fnam[] = "parse_WGTMAP_HOSTLIB"  ;

  // ----------- BEGIN -------------

  if ( strcmp(string,"NVAR_WGTMAP:")==0 ) {
    printf("\n WARNING: Should remove obsolete "
	   "NVAR_WGTMAP key from %s\n", INPUTS.HOSTLIB_FILE );
  }

  IVAR_STORE = HOSTLIB.NVAR_STORE ;

  FOUND_VARNAMES = ( strcmp(string,"VARNAMES_WGTMAP:") ==0 );
  if ( !FOUND_VARNAMES) { return ; }

    
  fgets(LINE,100,fp);
  NVAR_WGTMAP = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
  NDIM = NVAR_WGTMAP-2; NFUN=2;
  if ( NDIM < 1 ) {
    sprintf(c1err, "Invalid NDIM=%d for %s", NDIM, string);
    sprintf(c2err, "VARNAMES_WGTMAP must inclulde WGT & SNMAGSHIFT");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  
  // read in names used for weight map
  HOSTLIB_WGTMAP.GRIDMAP.VARLIST[0] = 0 ;
  for ( ivar=0; ivar < NVAR_WGTMAP ; ivar++ ) {
    VARNAME = HOSTLIB_WGTMAP.VARNAME[ivar] ;
    get_PARSE_WORD(0,ivar,VARNAME);
    
    strcat(HOSTLIB_WGTMAP.GRIDMAP.VARLIST,VARNAME);
    strcat(HOSTLIB_WGTMAP.GRIDMAP.VARLIST," ");
    
    // load variable if it's not already loaded
    if ( IVAR_HOSTLIB(VARNAME,0) < 0 && ivar < NDIM ) {
      sprintf(HOSTLIB.VARNAME_STORE[IVAR_STORE], "%s", VARNAME );
      IVAR_STORE++ ;
    }  
  } // end of ivar loop

    // read WGT keys and load GRIDMAP struct.
  read_GRIDMAP(fp,"WGT:", "", IDMAP, NDIM, NFUN, 0, 
	       MXWGT_HOSTLIB, fnam,
	       &HOSTLIB_WGTMAP.GRIDMAP ); // <== return GRIDMAP
  
  HOSTLIB_WGTMAP.WGTMAX = HOSTLIB_WGTMAP.GRIDMAP.FUNMAX[0];
  
  // update global counter
  HOSTLIB.NVAR_STORE = IVAR_STORE ;

  // check for optional WGTMAP_CHECK key to verify
  // WGTMAP interpolation.

  double TMPVAL[10];
  if ( strcmp(string,"WGTMAP_CHECK:") == 0 ) {
    readlong  (fp, 1, &GALID ); // Feb 2015
    readdouble(fp, 2, TMPVAL  );
    N = HOSTLIB_WGTMAP.NCHECKLIST ;
    HOSTLIB_WGTMAP.CHECKLIST_GALID[N] = GALID ;
    HOSTLIB_WGTMAP.CHECKLIST_ZTRUE[N] = TMPVAL[0] ;
    HOSTLIB_WGTMAP.CHECKLIST_WGT[N]   = TMPVAL[1] ;
    HOSTLIB_WGTMAP.CHECKLIST_SNMAG[N] = TMPVAL[2] ;
    HOSTLIB_WGTMAP.NCHECKLIST++ ;  
  }

  return ;

} // end of parse_WGTMAP_HOSTLIB


// ====================================
void  read_specbasis_HOSTLIB(void) {

  // Created Jun 28 2019 by R.Kessler
  // Read supplemental file of spectral templates, used later
  // to construct Host spectrum for each event.
  // This read must be done before reading the HOSTLIB,
  // so that the template names can be matched between
  // this specTemplate file and the VARNAMES in the HOSTLIB.

  FILE *fp;
  int  NBIN_WAVE, NBIN_READ, IFILETYPE;
  int  NVAR, ivar, ICOL_WAVE, NT, MEMD, NUM ;
  int  NVAR_WAVE  = 0 ;
  int  OPT_VARDEF = 0 ;
  int  LEN_PREFIX = strlen(PREFIX_SPECBASIS);
  char *ptrFile, *varName, c_get[60];  
  char TBLNAME[] = "SPECBASIS";
  char fnam[] = "read_specbasis_HOSTLIB";
  
  // --------------- BEGIN -----------------

  HOSTSPEC.NSPECBASIS  = 0 ;
  HOSTSPEC.NBIN_WAVE   = 0 ;
  HOSTSPEC.FLAM_SCALE        = 1.0 ;
  HOSTSPEC.FLAM_SCALE_POWZ1  = 0.0 ;

  ptrFile = INPUTS.HOSTLIB_SPECBASIS_FILE ;
  if ( IGNOREFILE(ptrFile) )  { return ; }

  printf("\n\t Read SPEC-TEMPLATEs from supplemental file:\n" );
  fflush(stdout);

  // - - - - - - - - - - - - - -
  // read until VARNAMES key in case there are supplemental keys
  if ( (fp = fopen(ptrFile, "rt")) == NULL ) {
    sprintf(c1err,"%s", "Could not open SPEC-TEMPLATE file:" );
    sprintf(c2err,"%s", ptrFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  

  int STOP=0;
  while( !STOP ) {
    fscanf(fp, "%s", c_get);
    if ( strcmp(c_get,"VARNAMES:") == 0 ) { STOP=1; }

    if ( strcmp(c_get,"FLAM_SCALE:") == 0 ) 
      { readdouble(fp, 1, &HOSTSPEC.FLAM_SCALE); }    

    if ( strcmp(c_get,"FLAM_SCALE_POWZ1:") == 0 ) 
      { readdouble(fp, 1, &HOSTSPEC.FLAM_SCALE_POWZ1); }    
  }

  fclose(fp);

  // - - - - - - - - - - - - - -
  // now read with standard routines
  TABLEFILE_INIT();
  NBIN_WAVE   = SNTABLE_NEVT(ptrFile,TBLNAME);
  IFILETYPE   = TABLEFILE_OPEN(ptrFile,"read");
  NVAR        = SNTABLE_READPREP(IFILETYPE,TBLNAME);
  MEMD        = (NBIN_WAVE+100) * sizeof(double);

  if ( NBIN_WAVE > MXBIN_SPECBASIS ) {
    sprintf(c1err,"NBIN_WAVE=%d exceeds bound of %d",
	    NBIN_WAVE, MXBIN_SPECBASIS );
    sprintf(c2err,"Reduce NBIN_WAVE or increase MXBIN_SPECBASIS");
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err); 
  }

  // example VARNAMES list to make sure that there is a wavelength column,
  // and count how many template[nn] colummns
  ICOL_WAVE = -9;  NT=0;
  for(ivar=0; ivar < NVAR; ivar++ ) {
    varName = READTABLE_POINTERS.VARNAME[ivar];
    if ( strstr(varName,"wave") != NULL ) { ICOL_WAVE = ivar; }
    if ( strstr(varName,"WAVE") != NULL ) { ICOL_WAVE = ivar; }
    if ( strstr(varName,"lam" ) != NULL ) { ICOL_WAVE = ivar; }
    if ( strstr(varName,"LAM" ) != NULL ) { ICOL_WAVE = ivar; }

    if ( ICOL_WAVE == ivar ) {
      NVAR_WAVE++ ;
      if ( NVAR_WAVE == 1 ) { 
	HOSTSPEC.WAVE_CEN     = (double*) malloc(MEMD);
	HOSTSPEC.WAVE_MIN     = (double*) malloc(MEMD);
	HOSTSPEC.WAVE_MAX     = (double*) malloc(MEMD);
	HOSTSPEC.WAVE_BINSIZE = (double*) malloc(MEMD);
	SNTABLE_READPREP_VARDEF(varName, HOSTSPEC.WAVE_CEN, NBIN_WAVE, 
				OPT_VARDEF); 
      }
    }  

   
    if ( strstr(varName,PREFIX_SPECBASIS) != NULL ) {
      if ( NT < MXSPECBASIS_HOSTLIB ) {
	sscanf(&varName[LEN_PREFIX], "%d",  &NUM);

	HOSTSPEC.ICOL_SPECBASIS[NT] = ivar;
	HOSTSPEC.NUM_SPECBASIS[NT]  = NUM; // store template NUM
	sprintf(HOSTSPEC.VARNAME_SPECBASIS[NT],"%s", varName);

	HOSTSPEC.FLAM_BASIS[NT] = (double*) malloc(MEMD);
	SNTABLE_READPREP_VARDEF(varName,HOSTSPEC.FLAM_BASIS[NT],
				NBIN_WAVE, OPT_VARDEF);
      }
      NT++ ; // always increment number of templates, NT
    }
   

  } // end ivar loop


  HOSTSPEC.FLAM_EVT = (double*) malloc(MEMD);

  // - - - - - - - - - - - - - -
  // abort tests
  if ( ICOL_WAVE < 0 ) {
    sprintf(c1err,"Could not find wavelength column");
    sprintf(c2err,"Check VARNAMES");
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err); 
  }

  if ( NVAR_WAVE != 1 ) {
    sprintf(c1err,"Found %d wavelength columns; expect 1", NVAR_WAVE);
    sprintf(c2err,"Check VARNAMES");
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err);     
  }

  if ( NT >= MXSPECBASIS_HOSTLIB ) {
    sprintf(c1err,"NSPECBASIS=%d exceeds bound of %d", 
	    NT, MXSPECBASIS_HOSTLIB ) ;
    sprintf(c2err,"Remove templates or increase MXSPECBASIS_HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err);     
  }
	
  HOSTSPEC.ICOL_WAVE = ICOL_WAVE ;
  HOSTSPEC.NSPECBASIS = NT;

  // read the entire table, and close it.
  NBIN_READ = SNTABLE_READ_EXEC();


  // Loop over wave bins and
  // + determine WAVE_BINSIZE for each wave bin
  // + truncate NBIN_WAVE so that lam < MAXLAM_SEDMODEL
  double WAVE_BINSIZE, WAVE_MIN, WAVE_MAX, LAM, UNIT;
  double LAM_LAST=0.0, LAM_NEXT=0.0;
  int FIRST, LAST, i, ilam, NBIN_KEEP=0 ;
  for(ilam=0; ilam < NBIN_WAVE; ilam++ ) {

    FIRST = ( ilam == 0 ) ;
    LAST  = ( ilam == NBIN_WAVE-1 ) ;

    LAM  = HOSTSPEC.WAVE_CEN[ilam];
    if ( LAM > LAMMAX_SEDMODEL ) { continue; }

    if ( !FIRST )  { LAM_LAST = HOSTSPEC.WAVE_CEN[ilam-1]; }
    if ( !LAST  )  { LAM_NEXT = HOSTSPEC.WAVE_CEN[ilam+1]; }

    if ( FIRST ) { LAM_LAST = LAM - (LAM_NEXT-LAM) ; }
    if ( LAST  ) { LAM_NEXT = LAM + (LAM-LAM_LAST) ; }

    WAVE_MIN     = LAM - (LAM - LAM_LAST)/2.0;
    WAVE_MAX     = LAM + (LAM_NEXT - LAM)/2.0;
    WAVE_BINSIZE = WAVE_MAX - WAVE_MIN;

    HOSTSPEC.WAVE_MIN[ilam]     = WAVE_MIN ;
    HOSTSPEC.WAVE_MAX[ilam]     = WAVE_MAX ;
    HOSTSPEC.WAVE_BINSIZE[ilam] = WAVE_BINSIZE;

    NBIN_KEEP++ ;

    if ( (ilam < -10) ) {
      printf(" xxx : ilam=%d: LAM[-1,0,+1] = %6.1f , %6.1f, %6.1f  "
	     "BINSIZE=%.1f \n",
	     ilam, LAM_LAST, LAM, LAM_NEXT, WAVE_BINSIZE);
      fflush(stdout);
    }

  } // end ilam loop

  NBIN_WAVE = NBIN_KEEP;
  HOSTSPEC.NBIN_WAVE  = NBIN_WAVE;

  printf("\t Found %d spectral basis vectors and %d wavelength bins\n",
	 NT, NBIN_WAVE);
  fflush(stdout);

  return;

} // end read_specbasis_HOSTLIB

// ============================================
void match_specbasis_HOSTVAR(void) {

  // Created June 2019
  // Match specTemplate names to names in HOSTLIB (HOSTVAR).
  // The specTemplate file has column names template00, template01, etc ...
  // The HOSTLIB VARNAMES must have corresponding list of
  // coeff_template00, coeff_template01, etc ...
  //

  int  NSPECBASIS    = HOSTSPEC.NSPECBASIS ;
  int  ivar_HOSTLIB, i, NERR=0;
  int  LDMP = 0 ;
  char *VARNAME_SPECBASIS, VARNAME_HOSTLIB[40];
  char fnam[] = "match_specbasis_HOSTVAR";

  // ----------------- BEGIN -----------------

  if ( NSPECBASIS == 0 ) { return ; }

  for(i=0; i < NSPECBASIS; i++ ) {
    VARNAME_SPECBASIS = HOSTSPEC.VARNAME_SPECBASIS[i];
    sprintf(VARNAME_HOSTLIB, "%s%s", 
	    PREFIX_SPECBASIS_HOSTLIB, VARNAME_SPECBASIS);

    ivar_HOSTLIB = IVAR_HOSTLIB(VARNAME_HOSTLIB,0);
    if ( ivar_HOSTLIB < 0 ) {
      printf(" ERROR: required '%s' is not in HOSTLIB' \n", VARNAME_HOSTLIB);
      NERR++; 
    }
    else {
      HOSTSPEC.IVAR_HOSTLIB[i] = ivar_HOSTLIB;
    }

    if ( LDMP ) 
      { printf(" xxx %s: ivar_HOSTLIB=%2d for '%s' \n",
	       fnam, ivar_HOSTLIB, VARNAME_SPECBASIS);
      }
  }


  if ( NERR > 0 ) {
    sprintf(c1err,"%d specTemplates have no coeff_template in HOSTLIB:", NERR );
    sprintf(c2err,"Check specTemplate file and HOSTLIB VARNAMES");
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err); 
  }

  return ;

} // end match_specbasis_HOSTVAR

// =========================================
void checkVarName_specTemplate(char *varName) {

  // Created Jun 2019
  // ABORT if varName (from HOSTLIB VARNAMES) is an 
  // un-used template coefficient.
  // 
  int icol;
  char fnam[] = "checkVarName_specTemplate";

  // --------------- BEGIN ------------

  if ( HOSTSPEC.NSPECBASIS <= 0 ) { return; }

  if (strstr(varName,PREFIX_SPECBASIS)         == NULL ) { return ; }
  if (strstr(varName,PREFIX_SPECBASIS_HOSTLIB) == NULL ) { return ; }

  // now we have varName of the form coeff_template[nn]
  icol = ICOL_SPECBASIS(varName,0); // fetch column in specTemplate file
  
  if ( icol < 0 ) {
    sprintf(c1err,"Found unused '%s' in HOSTLIB", varName );
    sprintf(c2err,"Check VARNAMES in specTemplate and HOSTLIB files.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return;

} // end checkVarName_specTemplate

// =========================================================
void genSpec_HOSTLIB(double zhel, double MWEBV, int DUMPFLAG,
		     double *GENFLUX_LIST, double *GENMAG_LIST) {

  // Created Jun 28 2019 by R.Kessler
  // Return host spectrum, including Galactic extinction.
  // If option is set to compute broadband mags, load them
  // into SNHOSTGAL.GALMAG[ifilt_obs][0] 
  //
  // Issues:
  //  - fraction of galaxy in fiber or slit ? Or does ETC include this ?
  //

  int      NBLAM_SPECTRO   = INPUTS_SPECTRO.NBIN_LAM;
  double  *LAMAVG_SPECTRO  = INPUTS_SPECTRO.LAMAVG_LIST;
  double  *ZP_SPECTRO      = SPECTROGRAPH_SEDMODEL.ZP_LIST;
  double   LAMMIN_SPECTRO  = LAMAVG_SPECTRO[0];
  double   LAMMAX_SPECTRO  = LAMAVG_SPECTRO[NBLAM_SPECTRO-1];

  int  NBLAM_BASIS = HOSTSPEC.NBIN_WAVE; 
  int  IGAL        = SNHOSTGAL.IGAL ;
  double z1        = 1.0 + zhel;
  double znorm     = pow(z1,HOSTSPEC.FLAM_SCALE_POWZ1) ;
  double hc8       = (double)hc;

  long long GALID;
  int  ilam, ilam_basis, ilam_last=-9, i, ivar_HOSTLIB, ivar ;
  int  ilam_near, NLAMSUM, LDMP=0;
  double FLAM_TMP, FLAM_SUM, COEFF, FLUX_TMP;
  double MWXT_FRAC, LAMOBS, LAMOBS_LAST=0.0;
  double LAMOBS_BIN, LAMOBS_MIN, LAMOBS_MAX, LAM_BASIS;
  double LAMREST_MIN, LAMREST_MAX ;
  double LAMMIN_TMP, LAMMAX_TMP, LAMBIN_TMP, LAMBIN_CHECK ;
  double ZP, FTMP, MAG, FLUX, ZCHECK;
  char fnam[] = "genSpec_HOSTLIB" ;

  // ------------------ BEGIN --------------

  if ( HOSTSPEC.NSPECBASIS <= 0 ) {
    sprintf(c1err,"Cannot generate host spectrum without spec basis.");
    sprintf(c2err,"Check HOSTLIB_SPECBASIS_FILE key in sim-input.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( DUMPFLAG ) {
    ivar = HOSTLIB.IVAR_GALID;
    GALID = (long long)HOSTLIB.VALUE_ZSORTED[ivar][IGAL] ;    

    ivar = HOSTLIB.IVAR_ZTRUE;
    ZCHECK = HOSTLIB.VALUE_ZSORTED[ivar][IGAL] ; 
    
    printf(" xxx ----------------------------------------------- \n");
    printf(" xxx %s DUMP for  GALID = %lld  (IGAL=%d, ZTRUE=%.4f)\n",
	   fnam, GALID, IGAL, ZCHECK );
  }

  // first construct total rest-frame spectrum using specbasis binning

  for(ilam_basis=0; ilam_basis < NBLAM_BASIS; ilam_basis++ ) {
    FLAM_SUM = 0.0 ;
    LAM_BASIS    = HOSTSPEC.WAVE_CEN[ilam_basis];
    for(i=0; i < HOSTSPEC.NSPECBASIS; i++ ) {
      ivar_HOSTLIB = HOSTSPEC.IVAR_HOSTLIB[i];
      COEFF        = HOSTLIB.VALUE_ZSORTED[ivar_HOSTLIB][IGAL] ; 
      FLAM_TMP     = HOSTSPEC.FLAM_BASIS[i][ilam_basis] ;
      FLAM_SUM    += ( COEFF * FLAM_TMP );

      if ( DUMPFLAG && ilam_basis == 0 && COEFF > 0.0 ) {
	printf(" xxx COEFF(%2d) = %le  (ivar_HOSTLIB=%d)\n", 
	       i, COEFF, ivar_HOSTLIB );
      }
    }

    // global scale for physical units
    HOSTSPEC.FLAM_EVT[ilam_basis] = (FLAM_SUM * HOSTSPEC.FLAM_SCALE * znorm);

    /* xxxxxxxxxxxxxxxx 
    LAMOBS  = z1*HOSTSPEC.WAVE_CEN[ilam_basis]; // obs frame
       // ZP is no good here because ZP is in SPECTROGRAPH bins,
       // not in HOSTSPEC bins.
    if ( DUMPFLAG ) {
      if ( LAMOBS > LAMMIN_SPECTRO && 
	   LAMOBS < LAMMAX_SPECTRO && LAMOBS-LAMOBS_LAST > 20.0 ) {
	
	LAMOBS_MIN   = z1*HOSTSPEC.WAVE_MIN[ilam_basis];
	LAMOBS_MAX   = z1*HOSTSPEC.WAVE_MAX[ilam_basis];
	LAMOBS_BIN   = z1*HOSTSPEC.WAVE_BINSIZE[ilam_basis];
	FLAM_TMP     = FLAM_BASIS[ilam_basis] ;
	ZP = interp_1DFUN(OPT_INTERP_LINEAR, LAMOBS, NBLAM_SPECTRO,
			  LAMAVG_SPECTRO, ZP_SPECTRO, fnam); // ZP NOT RIGHT

	FLUX    = FLAM_TMP * LAMOBS_BIN * LAMOBS/(hc8*z1);
	MAG     = ZP - 2.5*log10(FLUX);
	
	printf(" xxx %7.1f - %7.1f: FLAM(RAW,NORM) = %10.3le, %10.3le   "
	       "MAG=%5.2f (ZP=%5.2f) \n",
	       LAMOBS_MIN, LAMOBS_MAX, FLAM_SUM, FLAM_TMP, MAG, ZP );
	fflush(stdout);        LAMOBS_LAST = LAMOBS ;
      }

    } // end DUMPFLAG
    xxxxxxxxxxxxxxxx */    

  } // end ilam_basis


  // ---------------
  // HOST FLAM is in wavelength bins defined in the specTemplate file.
  // Here we to convert to SPECTROGRAPH bins in GENFLUX_LIST.
  // "ilam" is the index for SPECTROGRAPH, while ilam_basis is for specbasis.

  for(ilam=0; ilam < NBLAM_SPECTRO; ilam++ ) { 
    if ( MWEBV > 1.0E-9 ) 
      { MWXT_FRAC  = SEDMODEL_TABLE_MWXT_FRAC[0][ilam] ; }
    else
      { MWXT_FRAC  = 1.0 ; }

    GENFLUX_LIST[ilam] = 0.0;
    GENMAG_LIST[ilam]  = MAG_UNDEFINED ;

    LAMOBS      = SPECTROGRAPH_SEDMODEL.LAMAVG_LIST[ilam] ;
    LAMOBS_MIN  = SPECTROGRAPH_SEDMODEL.LAMMIN_LIST[ilam] ;
    LAMOBS_MAX  = SPECTROGRAPH_SEDMODEL.LAMMAX_LIST[ilam] ;
    LAMOBS_BIN  = LAMOBS_MAX - LAMOBS_MIN ;

    LDMP = (ilam == -217);

    LAMREST_MIN = LAMOBS_MIN/z1 ;
    LAMREST_MAX = LAMOBS_MAX/z1 ;

    if ( MWXT_FRAC < 1.0E-9 || MWXT_FRAC > 1.000001 ) {
      sprintf(c1err,"Invalid MWXT_FRAC = %f for LAMOBS=%.1f", 
	      MWXT_FRAC, LAMOBS );
      sprintf(c2err,"SEDMODEL_TABLE_MWXT_FRAC was probably not initialized.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // loop over specBasis bins that overlap this SPECTROGRAPH bin,
    // and sum specBasis flux
    FLUX_TMP  = 0.0 ;    LAMMIN_TMP = -9.0 ;   LAMBIN_CHECK=0.0 ;
    ilam_near = -9; NLAMSUM=0;
    ilam_basis = ilam_last-2 ;  if(ilam_basis<0) {ilam_basis=0;}

    while ( LAMMIN_TMP < LAMREST_MAX ) {
      FLAM_TMP     = HOSTSPEC.FLAM_EVT[ilam_basis];
      LAM_BASIS    = HOSTSPEC.WAVE_CEN[ilam_basis];
      LAMBIN_TMP   = HOSTSPEC.WAVE_BINSIZE[ilam_basis];
      LAMMIN_TMP   = HOSTSPEC.WAVE_MIN[ilam_basis];
      LAMMAX_TMP   = HOSTSPEC.WAVE_MAX[ilam_basis];

      ilam_basis++; ilam_last = ilam_basis;

      if ( LAM_BASIS < LAMOBS/z1 ) { ilam_near = ilam_basis; }

      if( LAMMAX_TMP < LAMREST_MIN ) { continue ; }
      if( LAMMIN_TMP > LAMREST_MAX ) { continue ; }

      // compute flux contained in this SPECTROGRAPH bin;
      // i.e. exclude flux that leaks out of thie SPECTROGRAPh bin.
      if ( LAMMIN_TMP < LAMREST_MIN ) { LAMMIN_TMP = LAMREST_MIN; }
      if ( LAMMAX_TMP > LAMREST_MAX ) { LAMMAX_TMP = LAMREST_MAX; }

      LAMBIN_CHECK += (LAMMAX_TMP - LAMMIN_TMP);
      FLUX_TMP     += ( FLAM_TMP * (LAMMAX_TMP - LAMMIN_TMP)*z1 );
      NLAMSUM++ ;

      if ( LDMP ) {
	printf(" xxx ilam=%d ilam_basis=%2d: LAM_BASIS=%.2f to %.2f  "
	       "FLUX=%.2le\n", 
	       ilam, ilam_basis, LAMMIN_TMP, LAMMAX_TMP, FLUX_TMP );
      }

    } // end while loop over rest-frame wavelength
    
    // if NLAMSUM==0, then no basis lambins overlap SPECTROGRAPH;
    // in thise case, average over nearest basis bins
    if ( NLAMSUM == 0 ) {
      double FLAM0 = HOSTSPEC.FLAM_EVT[ilam_near];
      double FLAM1 = HOSTSPEC.FLAM_EVT[ilam_near+1];
      FLAM_TMP     = (FLAM0+FLAM1)/2.0;
      FLUX_TMP     = FLAM_TMP * LAMOBS_BIN;
      LAMBIN_CHECK = LAMOBS_BIN/z1 ;
    }

    // check that sum over basis wave bins = SPECTROGRAPH bin size.
    // Basis wave is rest frame and SPECTROGRAPH bin is obs frame;
    // hence z1 factor needed to compare.

    if ( fabs(LAMBIN_CHECK - LAMOBS_BIN/z1) > 0.001 ) {
      printf("\n PRE-ABORT DUMP: \n");
      printf("   SPECTROGRAPH LAM(OBS) : %.3f to %.3f  (ilam=%d)\n",
	     LAMOBS_MIN, LAMOBS_MAX, ilam );
      printf("   SPECTROGRAPH LAM(Rest): %.3f to %.3f \n",
	     LAMREST_MIN, LAMREST_MAX);
      printf("   NLAMSUM = %d \n", NLAMSUM);
      sprintf(c1err,"Failed LAMBIN_CHECK=%.3f, but expected %.3f",
	      LAMBIN_CHECK, LAMOBS_BIN/z1);
      sprintf(c2err,"zhel=%.4f, MWXT_FRAC=%.3f", zhel, MWXT_FRAC );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // store flux (not FLAM) in SPECTROGRAPH bin
    GENFLUX_LIST[ilam] = FLUX_TMP * MWXT_FRAC ;  

    // convert to mag
    ZP    = SPECTROGRAPH_SEDMODEL.ZP_LIST[ilam] ;
    FTMP  = (LAMOBS/(hc8*z1)) * FLUX_TMP;
    if ( ZP > 0.0 && FTMP > 0.0 ) 
      { MAG = -2.5*log10(FTMP) + ZP; }
    else
      { MAG = MAG_UNDEFINED; }

    GENMAG_LIST[ilam]  = MAG;

  } // end ilam loop

  return;

} // end genSpec_HOSTLIB


// ============================================
void read_head_HOSTLIB(FILE *fp) {

  // Mar 6, 2011
  // read HOSTLIB header keys 
  // -  NVAR: %d
  // -  VARNMAES: <NVAR names>
  // -  NVAR_WGTMAP: %d
  // -  VARNAMES_WGTMAP:  <NVAR_WGTMAP names>
  // -  WGT: %f %f ... repeated NVAR_WGTMAP+2 times
  //
  // Jan 2014: add logic for SNPAR
  //
  // Feb 15 2018: fix bug to allow SNMAGSHIFT in wgtmap for MSKOPT&64
  // Dec 20 2018: use snana_rewind to allow for gzipped library
  // Mar 15 2019: ignore NVAR key, and parse entire line after VARNAMES
  // Mar 28 2019: use MXCHAR_LINE_HOSTLIB
  // Nov 11 2019: set HOSTLIB.IVAR_NBR_LIST

  int MXCHAR = MXCHAR_LINE_HOSTLIB;
  int ivar, ivar_map, IVAR_STORE, i, N, NVAR, NVAR_WGTMAP, FOUND_SNPAR;
  int MATCH, NVAR_STORE_SNPAR, USE, IS_SNPAR, ISTAT_VARNAMES, VBOSE ;
  int NCHAR;
  char  key[40], c_get[40], c_var[40], ctmp[80], wd[20], *cptr ;
  char  LINE[MXCHAR_LINE_HOSTLIB];
  char fnam[] = "read_head_HOSTLIB" ;

  // ------------- BEGIN ---------

  NVAR = NVAR_WGTMAP = 0 ;
  ISTAT_VARNAMES     = 0 ;  // change to 1 after reading
  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );
  NVAR_STORE_SNPAR = 0 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    // stop reading when first GAL: key is reached.
    if ( strcmp(c_get,"GAL:") == 0 )  { 
      // for later: snana_rewind(fp, INPUTS.HOSTLIB_FILE, HOSTLIB.GZIPFLAG );  
      fseek(fp,-4,SEEK_CUR); //rewind to before 1st GAL key
      goto VARCHECK ; 
    }

    if ( strcmp(c_get,"NVAR:")==0 ) { warn_NVAR_KEY(INPUTS.HOSTLIB_FILE); }

    sprintf(key, "%s", "VARNAMES:") ;
    if ( strcmp(c_get,key) == 0 ) {

      fgets(LINE,MXCHAR,fp);  NCHAR = strlen(LINE);

      if (NCHAR >= MXCHAR-5 ) {
	printf("\n PRE-ABORT DUMP\n");
	printf(" LINE = '%s'  (LEN=%d) \n", LINE, NCHAR );
	sprintf(c1err,"LINE likely exceeds bound of %d", MXCHAR);
	sprintf(c2err,"Shorten VARNAMES LINE, "
		"or increase MXCHAR_LINE_HOSTLIB");
        errmsg(SEV_FATAL, 0, fnam, c1err,c2err); 
      }

      NVAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE); 
      HOSTLIB.NVAR_ALL = NVAR;
      if ( NVAR < HOSTLIB.NVAR_REQUIRED || NVAR > MXVAR_HOSTLIB ) {
        sprintf(c1err, "NVAR= %d  is invalid (MXVAR_HOSTLIB=%d)",
                NVAR, MXVAR_HOSTLIB );
	sprintf(c2err, "NVAR must be between %d and %d",
		HOSTLIB.NVAR_REQUIRED, MXVAR_HOSTLIB );
        errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
      }

      ISTAT_VARNAMES = 1 ; 

      // parse the VARNAMES from LINE string
      for ( ivar=0; ivar < NVAR; ivar++ ) {
	get_PARSE_WORD(0,ivar,c_var);
	checkAlternateVarNames(c_var);

	// if coeff_tempalte[nn], make sure it's actually needed
	checkVarName_specTemplate(c_var);

	// load ALL array
	sprintf( HOSTLIB.VARNAME_ALL[ivar], "%s", c_var);

	// check for optional variables to add to the store list
	for ( i=0; i < HOSTLIB.NVAR_OPTIONAL; i++ ) {
	  cptr = HOSTLIB.VARNAME_OPTIONAL[i];
	  if( strcmp_ignoreCase(c_var,cptr) == 0 ) {
	    N = HOSTLIB.NVAR_STORE ;  // global var index (required+optional)
	    sprintf(HOSTLIB.VARNAME_STORE[N],"%s", c_var );
	    if ( HOSTLIB.IS_SNPAR_OPTIONAL[i] )  
	      { 
		HOSTLIB.FOUND_SNPAR_OPTIONAL[i] = 1 ;
		NVAR_STORE_SNPAR++ ; 
		HOSTLIB.IS_SNPAR_STORE[N] = 1 ; // index is for all variables
	      }
	    HOSTLIB.NVAR_STORE++ ;   
	  }
	}   // i       

      }  // end ivar loop

    } // end of VARNAMES


    // look for fixed Sersic index 'n#_Sersic' outside of VARNAMES list
    if ( ISTAT_VARNAMES ) 
      { parse_Sersic_n_fixed(fp,c_get); }
   

    // -----------
    // look for variables to use in weight-map 
    // (unless already read from elsewhere)

    if ( HOSTLIB_WGTMAP.ISTAT == 0 )
      { parse_WGTMAP_HOSTLIB(fp,c_get); }


  } // end of while-fscanf 


  // -----------------------------

 VARCHECK:  // note there is NO more file-reading below

  //-----------------------
  // sanity check on optioanl SNPARams
  USE = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) ;
  FOUND_SNPAR = ( NVAR_STORE_SNPAR>0  ) ;

  if ( USE && FOUND_SNPAR == 0 ) {
    sprintf(c1err, "Found zero SN params for HOSTLIB_MSKOPT&%d option.",
	    HOSTLIB_MSKOPT_USESNPAR );
    sprintf(c2err, "Expected 1 or more of '%s'", 
	    HOSTLIB.VARSTRING_SNPAR) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // make sure the required header was read
  if ( HOSTLIB.NVAR_ALL <= 0 ) {
    sprintf(c1err,"Could not find 'NVAR:' and/or 'VARNAMES:' keys.");
    sprintf(c2err,"Check HOSTLIB file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // make sure that each VARNAME_STORE maps onto
  // a real VARNAME_ALL, and store which 'ALL' index
  // Also load HOSTLIB.IVAR_XXX variables.

  NVAR_WGTMAP =  HOSTLIB_WGTMAP.GRIDMAP.NDIM ;
  for ( IVAR_STORE=0; IVAR_STORE < HOSTLIB.NVAR_STORE ; IVAR_STORE++ ) {

    sprintf(c_var,"%s", HOSTLIB.VARNAME_STORE[IVAR_STORE] );
    for ( ivar_map=0;  ivar_map < NVAR_WGTMAP; ivar_map++ ) {
      if ( strcmp(c_var,HOSTLIB_WGTMAP.VARNAME[ivar_map])==0 ) {
	HOSTLIB.IVAR_STORE[ivar_map]    = IVAR_STORE ; 
	HOSTLIB.IVAR_WGTMAP[IVAR_STORE] = ivar_map ;
      }      
    }

    // now look for match among ALL of the hostlib variables
    // Feb 2014 -> case insensitive check
    MATCH = 0;
    for( ivar=0; ivar < HOSTLIB.NVAR_ALL; ivar++ ) {
      if ( strcmp_ignoreCase(c_var,HOSTLIB.VARNAME_ALL[ivar]) == 0 ) { 
	MATCH = 1 ; 
	HOSTLIB.IVAR_ALL[IVAR_STORE] = ivar ;
	goto DONEMATCH ; 
      }
    } // end ivar

  DONEMATCH:
    if ( MATCH && VBOSE ) {
      IS_SNPAR = HOSTLIB.IS_SNPAR_STORE[IVAR_STORE];
      sprintf(wd,"HOSTLIB");
      if ( IS_SNPAR ) { sprintf(wd,"SNPAR  "); }

      sprintf(ctmp, "Found match to %s Variable", wd);
      printf("\t %s '%s' (IVAR_STORE=%d) \n", ctmp, c_var, IVAR_STORE ); 
      fflush(stdout);
    }
    if ( MATCH == 0 ) {
      sprintf(c1err,"Could not find required HOSTLIB var '%s' (IVAR_STORE=%d)", 
	      c_var, IVAR_STORE );
      sprintf(c2err,"Check HOSTLIB VARNAMES,WGTMAP & HOSTLIB_STOREVAR key.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
   
  } // end IVAR_STORE
  

  // make sure we found what we are looking for by checking 
  // that each var-string returns a valid 'ivar.
  // The 2nd argument to IVAR_HOSTLIB is a flag to
  // abort on any mis-match.

  HOSTLIB.IVAR_GALID   = IVAR_HOSTLIB(HOSTLIB_VARNAME_GALID,1) ; // required
  HOSTLIB.IVAR_ZTRUE   = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZTRUE,1) ; // required 

  // optional
  HOSTLIB.IVAR_ZPHOT       = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZPHOT,   0) ; 
  HOSTLIB.IVAR_ZPHOT_ERR   = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZPHOT_ERR,0);
  HOSTLIB.IVAR_LOGMASS     = IVAR_HOSTLIB(HOSTLIB_VARNAME_LOGMASS, 0) ; 
  HOSTLIB.IVAR_LOGMASS_ERR = IVAR_HOSTLIB(HOSTLIB_VARNAME_LOGMASS_ERR,0);
  HOSTLIB.IVAR_ANGLE       = IVAR_HOSTLIB(HOSTLIB_VARNAME_ANGLE,0) ;   
  HOSTLIB.IVAR_FIELD       = IVAR_HOSTLIB(HOSTLIB_VARNAME_FIELD,0) ;   
  HOSTLIB.IVAR_NBR_LIST    = IVAR_HOSTLIB(HOSTLIB_VARNAME_NBR_LIST,0) ; 

  // Jan 2015: Optional RA & DEC have multiple allowed keys
  int IVAR_RA[3], IVAR_DEC[3] ;
  IVAR_RA[0]   = IVAR_HOSTLIB(HOSTLIB_VARNAME_RA,0);
  IVAR_RA[1]   = IVAR_HOSTLIB(HOSTLIB_VARNAME_RA_HOST,0);
  IVAR_RA[2]   = IVAR_HOSTLIB(HOSTLIB_VARNAME_RA_GAL,0);
  IVAR_DEC[0]  = IVAR_HOSTLIB(HOSTLIB_VARNAME_DEC,0);
  IVAR_DEC[1]  = IVAR_HOSTLIB(HOSTLIB_VARNAME_DEC_HOST,0);
  IVAR_DEC[2]  = IVAR_HOSTLIB(HOSTLIB_VARNAME_DEC_GAL,0) ;
  HOSTLIB.IVAR_RA = HOSTLIB.IVAR_DEC = -9 ;
  for(i = 0; i < 3; i++ ) {
    if ( IVAR_RA[i]  > 0 ) { HOSTLIB.IVAR_RA  = IVAR_RA[i] ; }
    if ( IVAR_DEC[i] > 0 ) { HOSTLIB.IVAR_DEC = IVAR_DEC[i] ; }
  }


  // just make sure that these WGTMAP variables are really defined.
  // Also flag user-STOREPAR [EXTRA] variables that are also in WGTMAP (7.2019)
  int  NVAR_EXTRA  = HOSTLIB_OUTVAR_EXTRA.NOUT ;
  char *varName_WGTMAP, *varName_EXTRA;
  for ( ivar_map=0;  ivar_map < NVAR_WGTMAP; ivar_map++ )  { 
    varName_WGTMAP = HOSTLIB_WGTMAP.VARNAME[ivar_map] ;
    ivar = IVAR_HOSTLIB(varName_WGTMAP,1);  

    for(ivar=0; ivar < NVAR_EXTRA; ivar++ ) {
      varName_EXTRA = HOSTLIB_OUTVAR_EXTRA.NAME[ivar];
      if ( strcmp(varName_EXTRA,varName_WGTMAP) == 0 ) 
	{  HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[ivar] = 1 ; }

    }
  }

  return ;

} // end of read_head_HOSTLIB


// =====================================
void  checkAlternateVarNames(char *varName) {

  // Feb 12 2014
  // If input varName matches an allowed [hard-wired] alternative,
  // reset varName to the official name.
  // Note that the input argument is modified !

  // --------- BEGIN ---------

  if ( strcmp(varName,"ZERR") == 0 ) 
    { sprintf(varName,"%s", HOSTLIB_VARNAME_ZPHOT_ERR); }

  if ( strcmp(varName,"ZPHOTERR") == 0 ) 
    { sprintf(varName,"%s", HOSTLIB_VARNAME_ZPHOT_ERR); }

  if ( strcmp(varName,"LOGMASS_OBS") == 0 ) 
    { sprintf(varName,"%s", HOSTLIB_VARNAME_LOGMASS); }

  if ( strcmp(varName,"LOGMASS_OBS_ERR") == 0 ) 
    { sprintf(varName,"%s", HOSTLIB_VARNAME_LOGMASS_ERR); }

} // end of 	checkAlternateVarNames

// ====================================
void  parse_Sersic_n_fixed(FILE *fp, char  *string) {

  // look for fixed Sersic index defined after VARANMES.
  // Check all possible keys in case n#_Sersic is defined
  // for '#' that has no corresponding  a#_Sersic and  b#_Sersic.

  int j, N  ;
  double xn;
  char key[40], anam[12], bnam[12], wnam[12], nnam[12];

  //  char fnam[] = "parse_Sersic_n_fixed";

  // ------------ BEGIN ------------

  for ( j=0; j< MXSERSIC_HOSTLIB; j++ ) {

    Sersic_names(j, anam, bnam, wnam, nnam );

    sprintf(key,"%s:" , nnam );
    if ( strcmp(string,key) != 0 ) { continue ; }

    N = SERSIC_PROFILE.NFIX ;
    readdouble(fp, 1, &xn ) ;
    SERSIC_PROFILE.FIX_VALUE[N] = xn ;      
    sprintf(SERSIC_PROFILE.FIX_NAME[N], "%s", nnam);
    SERSIC_PROFILE.NFIX++ ; 
  }

} // end of parse_Sersic_n_fixed


// ====================================
void read_gal_HOSTLIB(FILE *fp) {

  // Mar 6, 2011
  // Read HOSTLIB by reading info following each "GAL:" key.
  // Store only the variables requested to be stored.
  //
  // Sep 22, 2011: apply RA+DEC cuts
  // Feb 13 2014: fix safety margin for LOGZCUT
  //
  // Sep 16 2015: call read_galRow_HOSTLIB() to allow reading FIELD string.
  // Feb 16 2016: check how many in GALID_PRIORITY
  //
  // Dec 29 2017: time to read 416,000 galaxies 5 sec ->
  // Nov 11 2019: check NBR_LIST

  char c_get[40], FIELD[MXCHAR_FIELDNAME], NBR_LIST[MXCHAR_NBR_LIST] ;
  char fnam[] = "read_gal_HOSTLIB"  ;
  
  long long GALID, GALID_MIN, GALID_MAX ;
  int  ivar_ALL, ivar_STORE, NVAR_STORE, NGAL, NGAL_READ, MEMC ;
  int  NPRIORITY;

  double xval[MXVAR_HOSTLIB], val, ZTMP, LOGZCUT[2], DLOGZ_SAFETY ;

  // ----------------- BEGIN -----------

  NVAR_STORE = HOSTLIB.NVAR_STORE;

  GALID_MIN = INPUTS.HOSTLIB_GALID_PRIORITY[0] ;
  GALID_MAX = INPUTS.HOSTLIB_GALID_PRIORITY[1] ;
  NPRIORITY = 0 ;

  // define cut windows
  HOSTLIB_CUTS.RAWIN[0] = INPUTS.HOSTLIB_GENRANGE_RA[0] - .0001 ;
  HOSTLIB_CUTS.RAWIN[1] = INPUTS.HOSTLIB_GENRANGE_RA[1] + .0001 ;

  HOSTLIB_CUTS.DECWIN[0] = INPUTS.HOSTLIB_GENRANGE_DEC[0] - .0001 ;
  HOSTLIB_CUTS.DECWIN[1] = INPUTS.HOSTLIB_GENRANGE_DEC[1] + .0001 ;

  // define redshift cut as user GENRANGE_REDSHIFT with
  // an extra safety margin (in logz space) of 3 bins.
  DLOGZ_SAFETY = 3.0 * DZPTR_HOSTLIB ;
  ZTMP       = INPUTS.GENRANGE_REDSHIFT[0] ;
  LOGZCUT[0] = log10(ZTMP-0.01) - DLOGZ_SAFETY ;
  ZTMP       = INPUTS.GENRANGE_REDSHIFT[1] ;
  LOGZCUT[1] = log10(ZTMP+0.01) + DLOGZ_SAFETY ;
  HOSTLIB_CUTS.ZWIN[0] = pow(10.0, LOGZCUT[0] );
  HOSTLIB_CUTS.ZWIN[1] = pow(10.0, LOGZCUT[1] );

  NGAL = -9;


  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"GAL:") == 0 ) {

      NGAL_READ = HOSTLIB.NGAL_READ; // C-like index
      HOSTLIB.NGAL_READ++ ;          // fortran-like index

      read_galRow_HOSTLIB(fp, HOSTLIB.NVAR_ALL, xval, FIELD, NBR_LIST ); 

      if ( HOSTLIB.NGAL_READ > INPUTS.HOSTLIB_MAXREAD ) 
	{ goto DONE_RDGAL ; } 

      if ( passCuts_HOSTLIB(xval) == 0 ) { continue; }

      // count how many priority entries (for print summary below)
      if ( GALID_MIN < GALID_MAX ) {
	ivar_ALL    = HOSTLIB.IVAR_ALL[HOSTLIB.IVAR_GALID] ;
	GALID       = (long long)xval[ivar_ALL];
	if ( GALID >= GALID_MIN && GALID <= GALID_MAX ) { NPRIORITY++ ; }
      }


      // store this galaxy
      NGAL = HOSTLIB.NGAL_STORE ;   HOSTLIB.NGAL_STORE++ ;   

      // check to allocate more storage/memory
      malloc_HOSTLIB(NGAL,NGAL_READ);    

      // strip off variables to store
      for ( ivar_STORE=0; ivar_STORE < NVAR_STORE; ivar_STORE++ ) {
	ivar_ALL = HOSTLIB.IVAR_ALL[ivar_STORE] ;  // column in file
	val      = xval[ivar_ALL] ;

	HOSTLIB.VALUE_UNSORTED[ivar_STORE][NGAL] = val ;

	// keep track of min and max for each variable
	if ( val > HOSTLIB.VALMAX[ivar_STORE] ) 
	  { HOSTLIB.VALMAX[ivar_STORE] = val; }
	if ( val < HOSTLIB.VALMIN[ivar_STORE] ) 
	  { HOSTLIB.VALMIN[ivar_STORE] = val; }

      }

      // store optional FIELD string 
      ivar_STORE = HOSTLIB.IVAR_FIELD ;
      if ( ivar_STORE > 0  ) {
	ivar_ALL   = HOSTLIB.IVAR_ALL[ivar_STORE];
	sprintf(HOSTLIB.FIELD_UNSORTED[NGAL],"%s", FIELD);
      }


      // Nov 11 2019: store optional NBR_LIST string 
      ivar_STORE = HOSTLIB.IVAR_NBR_LIST ;
      if ( ivar_STORE > 0  ) {
	ivar_ALL   = HOSTLIB.IVAR_ALL[ivar_STORE];
	MEMC = strlen(NBR_LIST) * sizeof(char);
	HOSTLIB.NBR_UNSORTED[NGAL] = (char*) malloc(MEMC);
	sprintf(HOSTLIB.NBR_UNSORTED[NGAL],"%s", NBR_LIST);

	// xxxxxxxx .xyz
	if ( NGAL < 20 ) {
	  printf(" xxx %s: NGAL=%2d: GALID=%lld  NBR_LIST = '%s' \n", 
		 GALID, NBR_LIST); fflush(stdout);
	}
	//xxxxxxxxxxxx
      }

      // store NGAL index vs. absolute READ index (for HOSTNBR)
      HOSTLIB.LIBINDEX_READ[NGAL_READ] = NGAL ;

    }
  } // end of while-fscanf


 DONE_RDGAL:

  // load min/max ZTRUE into separate variables
  ivar_STORE   = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZTRUE,1) ;
  HOSTLIB.ZMIN = HOSTLIB.VALMIN[ivar_STORE];
  HOSTLIB.ZMAX = HOSTLIB.VALMAX[ivar_STORE];

  printf("\t Stored %d galaxies from HOSTLIB. \n",
	 HOSTLIB.NGAL_STORE );

  if ( NPRIORITY > 0 ) {
    printf("\t   --> %d galaxies have priority "
	   "(GALID between %lld and %lld)\n",
	   NPRIORITY, GALID_MIN, GALID_MAX );
  }

  // if user did not specify photo-z outlier range,
  // then use HOSTLIB z-range as default.
  if ( INPUTS.HOSTLIB_GENZPHOT_OUTLIER[0] < 0.0 ) {
    INPUTS.HOSTLIB_GENZPHOT_OUTLIER[0] = HOSTLIB.ZMIN ;
    INPUTS.HOSTLIB_GENZPHOT_OUTLIER[1] = HOSTLIB.ZMAX ;
  }
  
  fflush(stdout);

  // abort if requested z-range exceed range of HOSTLIB
  if ( HOSTLIB.ZMIN > INPUTS.GENRANGE_REDSHIFT[0] ||
       HOSTLIB.ZMAX < INPUTS.GENRANGE_REDSHIFT[1] ) {

    sprintf(c1err,"HOSTLIB z-range (%f - %f) does not contain",
	    HOSTLIB.ZMIN, HOSTLIB.ZMAX );
    sprintf(c2err,"GENRANGE_REDSHIFT (%f - %f)"
	    ,INPUTS.GENRANGE_REDSHIFT[0]
	    ,INPUTS.GENRANGE_REDSHIFT[1] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 

  }
  
  return ;

} // end of read_gal_HOSTLIB

// ====================================
int passCuts_HOSTLIB(double *xval ) {

  // Return 1 if cuts are satisfied; zero otherwise.
  int ivar_ALL, LRA ,LRA2;
  double ZTRUE, RA, RA2, DEC;

  // ---------- BEGIN ---------

  // if computing host mags, do not apply cuts
  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_PLUSMAGS)>0 ) 
    { return(1); }

  // ditto for HOST neighbors
  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_PLUSNBR)>0 ) 
    { return(1); }

  // REDSHIFT
  ivar_ALL    = HOSTLIB.IVAR_ALL[HOSTLIB.IVAR_ZTRUE] ; 
  ZTRUE       = xval[ivar_ALL];
  if ( ZTRUE < HOSTLIB_CUTS.ZWIN[0] ) { return(0); }
  if ( ZTRUE > HOSTLIB_CUTS.ZWIN[1] ) { return(0); }

  // RA
  if ( HOSTLIB.IVAR_RA > 0 && HOSTLIB.IVAR_DEC > 0 ) { 
    ivar_ALL    = HOSTLIB.IVAR_ALL[HOSTLIB.IVAR_RA] ;
    RA          = xval[ivar_ALL];
    RA2         = RA + 360.0 ; 
    LRA   = (RA  > HOSTLIB_CUTS.RAWIN[0] && RA  < HOSTLIB_CUTS.RAWIN[1] );
    LRA2  = (RA2 > HOSTLIB_CUTS.RAWIN[0] && RA2 < HOSTLIB_CUTS.RAWIN[1] );
    if ( LRA == 0 && LRA2 == 0 ) { return(0) ; }
    
    // DEC
    ivar_ALL    = HOSTLIB.IVAR_ALL[HOSTLIB.IVAR_DEC] ;
    DEC         = xval[ivar_ALL];   
    if ( DEC  < HOSTLIB_CUTS.DECWIN[0] ) { return(0) ; }
    if ( DEC  > HOSTLIB_CUTS.DECWIN[1] ) { return(0) ; }
  }

  return(1);

} // end passCuts_HOSTLIB

// ==========================================
void read_galRow_HOSTLIB(FILE *fp, int NVAL, double *VALUES, 
			 char *FIELD, char *NBR_LIST ) {

  // Created 9/16.2015
  // If there is no FIELD key, then read NVAL double from fp.
  // If there is a FIELD key, then read doubles and string separately.
  //
  // Dec 29 2017: use fgets for faster read.
  // Mar 28 2019: use MXCHAR_LINE_HOSTLIB
  // Nov 11 2019: check for NBR_LIST (string), analogous to FIELD check

  int MXCHAR = MXCHAR_LINE_HOSTLIB;
  int MXWD = NVAL;
  int ival_FIELD, ival_NBR_LIST, ival, NWD=0, len, NCHAR ;
  char WDLIST[MXVAR_HOSTLIB][40], *ptrWDLIST[MXVAR_HOSTLIB] ;
  char sepKey[] = " " ;
  char tmpWORD[200], tmpLine[MXCHAR_LINE_HOSTLIB], *pos ;
  char fnam[] = "read_galRow_HOSTLIB" ;
  // ---------------- BEGIN -----------------

  // scoop up rest of line with fgets
  fgets(tmpLine, MXCHAR, fp);

  NCHAR = strlen(tmpLine);
  if ( NCHAR >= MXCHAR-5 ) {
    printf("\n PRE-ABORT DUMP\n");
    printf(" LINE = '%s' (LEN=%d) \n", tmpLine, NCHAR );
    sprintf(c1err,"LINE likely exceeds bound of %d", MXCHAR);
    sprintf(c2err,"Shorten HOSTLIB lines, or increase MXCHAR_LINE_HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // remove line feed
  if ( (pos=strchr(tmpLine,'\n') ) != NULL )  { *pos = '\0' ; }

  // split string into words
  for(ival=0; ival < NVAL; ival++ ) { ptrWDLIST[ival] = WDLIST[ival]; }
  splitString2(tmpLine, sepKey, MXWD, &NWD, ptrWDLIST);

  // abort if too few columns, but allow extra columns
  // (e..g, comment or catenated files with extra columns)
  if ( NWD < NVAL ) {
    printf("\n PRE-ABORT DUMP: \n");
    printf("\t NGAL_READ = %d \n",  HOSTLIB.NGAL_READ );
    printf("\t LINE = '%s %s %s  ... ' \n", 
	   WDLIST[0], WDLIST[1], WDLIST[2] );
    sprintf(c1err,"Found %d words after GAL key", NWD);
    sprintf(c2err,"but expected %d words.", NVAL );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  sprintf(FIELD,"NULL"); ival_FIELD = -9;
  if ( HOSTLIB.IVAR_FIELD > 0 ) 
    { ival_FIELD = HOSTLIB.IVAR_ALL[HOSTLIB.IVAR_FIELD] ;  }

  sprintf(NBR_LIST,"NULL"); ival_NBR_LIST = -9;
  if ( HOSTLIB.IVAR_NBR_LIST > 0 ) 
    { ival_NBR_LIST = HOSTLIB.IVAR_ALL[HOSTLIB.IVAR_NBR_LIST] ;  }


  for(ival=0; ival < NVAL; ival++ ) {
    if ( ival == ival_FIELD )  { 
      sprintf(tmpWORD, "%s", WDLIST[ival] );      len = strlen(tmpWORD);
      if ( len > MXCHAR_FIELDNAME ) {
	sprintf(c1err,"strlen(FIELD=%s) = %d exceeds storage array of %d",
		tmpWORD, len, MXCHAR_FIELDNAME);
	sprintf(c2err,"Check HOSTLIB or increase MXCHAR_FIELDNAME");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      sscanf(tmpWORD, "%s", FIELD ); 
    }

    else if ( ival == ival_NBR_LIST )  { 
      sprintf(tmpWORD, "%s", WDLIST[ival] );      len = strlen(tmpWORD);
      if ( len > MXCHAR_NBR_LIST ) {
	sprintf(c1err,"strlen(NBR_LIST=%s) = %d exceeds storage array of %d",
		tmpWORD, len, MXCHAR_NBR_LIST);
	sprintf(c2err,"Check HOSTLIB or increase MXCHAR_NBR_LIST");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      sscanf(tmpWORD, "%s", NBR_LIST ); 
    }

    else {
      sscanf(WDLIST[ival], "%le", &VALUES[ival] ); 
    }

  } // end ival loop

  
  /*
  printf(" xxx NVAL=%d  N1=%d  N2=%d  FIELD=%s  ZTRUE=%.5f\n",
	 NVAL, N1, N2, FIELD, VALUES[2]);
  debugexit("read row"); // xxx
  */

  return ;

}  // read_galRow_HOSTLIB

// ==========================================
void  summary_snpar_HOSTLIB(void) {

  // Feb 2014
  // if use-SNPAR bit of HOSTLIB_MSKOPT is set,
  // then print list of SNPARAMS that were found,
  // and also a list of allowed SNparams that
  // not found. 

  int i, NVAR, LFOUND, NVAR_FOUND=0 ;
  char
     tmpVar[60]
    ,VARLIST_FOUND[MXPATHLEN]
    ,VARLIST_NOTFOUND[MXPATHLEN]
    ;

  //   char fnam[] = "summary_snpar_HOSTLIB" ;

  // ------------ BEGIN ------------

  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) == 0 ) 
    { return ; }

  NVAR = HOSTLIB.NVAR_OPTIONAL; 
  VARLIST_FOUND[0] = 0 ;
  VARLIST_NOTFOUND[0] = 0 ;

  for ( i=0; i < NVAR; i++ ) {
    if ( HOSTLIB.IS_SNPAR_OPTIONAL[i] == 0 ) { continue ; }

    sprintf(tmpVar, "%s ", HOSTLIB.VARNAME_OPTIONAL[i]);
    LFOUND = HOSTLIB.FOUND_SNPAR_OPTIONAL[i] ;
    if ( LFOUND ) 
      { strcat(VARLIST_FOUND, tmpVar);  NVAR_FOUND++; }
    else
      { strcat(VARLIST_NOTFOUND, tmpVar); }    
  }

  printf("     Valid SNPAR Variables Found:     '%s' \n", VARLIST_FOUND );
  printf("     Valid SNPAR Variables NOT Found: '%s' \n", VARLIST_NOTFOUND );
  fflush(stdout);

} // end of   summary_snpar_HOSTLIB

// ==================================
void malloc_HOSTLIB(int NGAL_STORE, int NGAL_READ) {

  // Mar 2011
  // allocate memory with malloc or extend memory with realloc
  // The value of NGAL determines if/when to allocate more memory.
  //
  // Nov 2019: 
  //  + add NGAL_READ argument to malloc HOSTLIB.LIBINDEX_READ
  //  + check IVAR_NBR_LIST

  double XNTOT, XNUPD;
  int ivar, I8, I8p, I4, DO_FIELD, DO_NBR, igal ;
  int LDMP = 0 ;
  char fnam[] = "malloc_HOSTLIB";

  // ------------- BEGIN ----------

  I8  = sizeof(double);
  I8p = sizeof(double*);
  I4  = sizeof(int);

  DO_FIELD    = ( HOSTLIB.IVAR_FIELD    > 0 );
  DO_NBR      = ( HOSTLIB.IVAR_NBR_LIST > 0 );

  if ( NGAL_READ == 0 ) {
    if  ( LDMP ) 
      { printf("\t xxx Initial malloc for HOSTLIB ... \n"); fflush(stdout); }

    HOSTLIB.MALLOCSIZE_D += (I8 * MALLOCSIZE_HOSTLIB) ;
    // allocate memory for each variable to store
    for ( ivar = 0; ivar < HOSTLIB.NVAR_STORE; ivar++ ) {
      HOSTLIB.VALUE_UNSORTED[ivar] = (double *)malloc(HOSTLIB.MALLOCSIZE_D);
    }

    HOSTLIB.MALLOCSIZE_I += (I4 * MALLOCSIZE_HOSTLIB) ;
    HOSTLIB.LIBINDEX_READ = (int *)malloc(HOSTLIB.MALLOCSIZE_I);

    if ( DO_FIELD ) {
      HOSTLIB.FIELD_UNSORTED = 
	(char**)malloc( MALLOCSIZE_HOSTLIB*sizeof(char*) );
      for(igal=0; igal < MALLOCSIZE_HOSTLIB; igal++ ) {
	HOSTLIB.FIELD_UNSORTED[igal] = 
	  (char*)malloc( MXCHAR_FIELDNAME *sizeof(char) );
      }  
    }

    if ( DO_NBR ) {
      HOSTLIB.NBR_UNSORTED = 
	(char**)malloc( MALLOCSIZE_HOSTLIB*sizeof(char*) );
      // do NOT malloc string length here; malloc later
      // when length of string is known.
    }

    return ;
  }


  /* xxx mark deletew
  // check when to extend memory
  XNTOT = (double)NGAL ;
  XNUPD = (double)MALLOCSIZE_HOSTLIB ;
  if ( fmod(XNTOT,XNUPD) == 0.0 ) {
  */

  if ( (NGAL_STORE % MALLOCSIZE_HOSTLIB) == 0 ) {

    if  ( LDMP )  { 
      printf("\t xxx Extend HOSTLIB malloc at NGAL_STORE=%d \n", 
	     NGAL_STORE);        fflush(stdout); 
    }
    
    HOSTLIB.MALLOCSIZE_D += (I8*MALLOCSIZE_HOSTLIB) ;
    for ( ivar = 0; ivar < HOSTLIB.NVAR_STORE; ivar++ ) {
      HOSTLIB.VALUE_UNSORTED[ivar] =
	(double *)realloc(HOSTLIB.VALUE_UNSORTED[ivar], HOSTLIB.MALLOCSIZE_D);
    }

    if ( DO_FIELD ) {
      HOSTLIB.FIELD_UNSORTED = 
	(char**)realloc( HOSTLIB.FIELD_UNSORTED, MALLOCSIZE_HOSTLIB );
      for(igal=NGAL_STORE; igal < NGAL_STORE+MALLOCSIZE_HOSTLIB; igal++ ) {
	HOSTLIB.FIELD_UNSORTED[igal] = 
	  (char*)malloc( MXCHAR_FIELDNAME *sizeof(char) );
      }  
    }

    if ( DO_NBR ) {
      HOSTLIB.NBR_UNSORTED = 
	(char**)realloc( HOSTLIB.NBR_UNSORTED, MALLOCSIZE_HOSTLIB );
      // do NOT malloc size of each string
    }

  } // end if block

  // separate check for READ index
  if ( (NGAL_READ % MALLOCSIZE_HOSTLIB) == 0 ) {
    HOSTLIB.MALLOCSIZE_I += (I4 * MALLOCSIZE_HOSTLIB) ;
    HOSTLIB.LIBINDEX_READ = 
      (int *)realloc(HOSTLIB.LIBINDEX_READ,HOSTLIB.MALLOCSIZE_I);
    /*
    printf(" 2. xxx %s: HOSTLIB.MALLOCSIZE_I = %d \n", 
	   fnam, HOSTLIB.MALLOCSIZE_I); fflush(stdout);
    */
  }
  
} // end of malloc_HOSTLIB



// =============================================
void check_duplicate_GALID(void) {

  // Sort hostlib by GALID, then abort if any 
  // [GALID,ZTRUE] pair appears more than once.
  //
  // Jan 2013: replace CERNLIB sortzv with sortInt wrapper.
  //

  long long GALID, GALID_LAST ;
  double   *ptrGALID ;

  int  *LIBINDEX_UNSORT ;
  int  NGAL, igal, IVAR_GALID, NERR, VBOSE, isort, ORDER_SORT  ;

  double ZTRUE, ZTRUE_LAST, ZDIF, XGAL;
  double ZTOL = 0.00001 ;    // z-tolerance   

  char fnam[] = "check_duplicate_GALID" ;

  // ---------- BEGIN ----------

  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );

  if ( VBOSE ) { 
    printf("\t Check HOSTLIB for duplicate entries. \n"); 
    fflush(stdout);
  }


  NGAL = HOSTLIB.NGAL_STORE ;
  IVAR_GALID = HOSTLIB.IVAR_GALID ;


  // allocate memory for sort-pointers
  LIBINDEX_UNSORT = (int   *)malloc( (NGAL+1) * sizeof(int)    );
  ptrGALID        = (double*)malloc( (NGAL+1) * sizeof(double) );

  // load array of double-precision GALID to be sorted
  for ( igal = 0; igal < NGAL; igal++ ) {
    XGAL = HOSTLIB.VALUE_ZSORTED[IVAR_GALID][igal] ;
    ptrGALID[igal] = XGAL ;
  }

  ORDER_SORT = +1 ; // increasing order
  sortDouble( NGAL, ptrGALID, ORDER_SORT, 
	      LIBINDEX_UNSORT );  // return array of indices

  GALID_LAST  =  -9  ;
  ZTRUE_LAST  = 0.0 ;
  NERR = 0;

  // check current and previous GALID for duplicates
  for ( igal = 0; igal < NGAL ; igal++ ) {

    isort = LIBINDEX_UNSORT[igal] ;
    GALID = get_GALID_HOSTLIB(isort) ;
    ZTRUE = get_ZTRUE_HOSTLIB(isort) ;

    ZDIF = fabs(ZTRUE - ZTRUE_LAST);

    /*
    printf(" xxx igal=%4d  isort=%4d  GALID=%lld  GALID_LAST=%lld   ZTRUE=%.5f\n",
	   igal, isort, GALID, GALID_LAST, ZTRUE); fflush(stdout);
    */

    if ( GALID == GALID_LAST &&   ZDIF < ZTOL ) {
      NERR++ ;
      printf(" !*!*!* Duplicate GALID = %lld  at ZTRUE=%f ZDIF=%f !*!*! \n", 
	     GALID, ZTRUE, ZDIF  );
      fflush(stdout);
    }

    GALID_LAST = GALID ;
    ZTRUE_LAST = ZTRUE  ;
  }

  if ( NERR ) {
    sprintf(c1err,"Found duplicate GALID (galaxy IDs)");
    sprintf(c1err,"Check HOSTLIB.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  free(LIBINDEX_UNSORT);
  free(ptrGALID);

}   // end of check_duplicate_GALID

// ============================
void sortz_HOSTLIB(void) {

  // Mar 2011
  // Use CERNLIB sortz function to sort library
  // by redshift. Note that sortzv accepts a 
  // real*4 array, so be careful.
  //
  // Also compute ZGAPMAX, ZGAPAVG and Z_ATGAPMAZ
  //

  int  NGAL, igal, ival, unsort, VBOSE, DOFIELD;
  int  IVAR_ZTRUE, NVAR_STORE, ORDER_SORT     ;

  double ZTRUE, ZLAST, ZGAP, ZSUM, *ZSORT, VAL ;
  char *FIELD;
  char fnam[] = "sortz_HOSTLIB" ;

  // ------------- BEGIN -------------

  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );

  if ( VBOSE )  { 
    printf("\t Sort HOSTLIB by redshift (%.4f to %.4f) \n",
	   HOSTLIB.ZMIN , HOSTLIB.ZMAX );   
    printf("\t |zSN-zGAL| tolerance: %.3f + %.3f*z + %.4f*z^2 \n"
	   ,INPUTS.HOSTLIB_DZTOL[0]
	   ,INPUTS.HOSTLIB_DZTOL[1]
	   ,INPUTS.HOSTLIB_DZTOL[2] );
    fflush(stdout); 
  }

  DOFIELD = ( HOSTLIB.IVAR_FIELD > 0 ) ;

  NGAL = HOSTLIB.NGAL_STORE ;
  NVAR_STORE = HOSTLIB.NVAR_STORE ;
  IVAR_ZTRUE = HOSTLIB.IVAR_ZTRUE ;

  // allocate memory for sort-pointers
  HOSTLIB.LIBINDEX_UNSORT  = (int*)malloc( (NGAL+1) * sizeof(int) );
  HOSTLIB.LIBINDEX_ZSORT   = (int*)malloc( (NGAL+1) * sizeof(int) );

  // allocate memory for sorted values
  for ( ival=0; ival < NVAR_STORE; ival++ ) {
    HOSTLIB.VALUE_ZSORTED[ival] = 
      (double*)malloc( (NGAL+1) * sizeof(double) ) ;
  }

  if ( DOFIELD  ) {
    HOSTLIB.FIELD_ZSORTED = 
      (char**)malloc( (NGAL+1) * sizeof(char*) ) ;
    for ( igal=0; igal <= NGAL; igal++ ) {
      HOSTLIB.FIELD_ZSORTED[igal] = 
	(char*)malloc( MXCHAR_FIELDNAME * sizeof(char) ) ; 
    }
  }


  // allocate memory and load ZSORT need for sorting routine
  ZSORT = (double*)malloc( (NGAL+1) * sizeof(double) );

  // load ZSORT array
  for ( igal=0; igal < NGAL; igal++ ) {
    ZTRUE = HOSTLIB.VALUE_UNSORTED[IVAR_ZTRUE][igal]; 
    ZSORT[igal] = ZTRUE ;
  }

  ORDER_SORT = +1 ;    // increasing order
  sortDouble( NGAL, ZSORT, ORDER_SORT, HOSTLIB.LIBINDEX_UNSORT ) ;


  HOSTLIB.SORTFLAG = 1 ;
  ZLAST = HOSTLIB.ZMIN ;
  ZSUM = 0.0 ;

  // fill sorted array. 'igal' is the z-sorted index; 
  // 'unsort' is the  un-sorted index matching the original HOSTLIB order.
  for ( igal = 0; igal < NGAL ; igal++ ) {
    unsort = HOSTLIB.LIBINDEX_UNSORT[igal]  ;
    HOSTLIB.LIBINDEX_ZSORT[unsort] = igal;

    for ( ival=0; ival < NVAR_STORE; ival++ ) {
      VAL = HOSTLIB.VALUE_UNSORTED[ival][unsort] ; 
      HOSTLIB.VALUE_ZSORTED[ival][igal] = VAL;
    }

    if ( DOFIELD ) {
      FIELD = HOSTLIB.FIELD_UNSORTED[unsort] ;
      sprintf(HOSTLIB.FIELD_ZSORTED[igal],"%s", FIELD);
    }

    // update redshift variables
    ZTRUE  = get_ZTRUE_HOSTLIB(igal);
    ZGAP   = ZTRUE - ZLAST;
    ZSUM  += (double)ZGAP;

    if ( ZGAP > HOSTLIB.ZGAPMAX ) {
      HOSTLIB.ZGAPMAX = ZGAP ;
      HOSTLIB.Z_ATGAPMAX[0] = ZLAST ;
      HOSTLIB.Z_ATGAPMAX[1] = ZTRUE ;
    }
    ZLAST = ZTRUE;

  } // end of igal loop

  HOSTLIB.ZGAPAVG = ZSUM/(double)(NGAL);

  // free memory for the pointers and the unsorted array.
  free(ZSORT);
  for ( ival=0; ival < NVAR_STORE; ival++ ) 
    { free(HOSTLIB.VALUE_UNSORTED[ival]);  }

  int  OPT_PLUSMAGS  = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_PLUSMAGS);
  int  OPT_PLUSNBR   = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_PLUSNBR);
  if ( !(OPT_PLUSMAGS || OPT_PLUSNBR) ) {
    free(HOSTLIB.LIBINDEX_UNSORT);
    free(HOSTLIB.LIBINDEX_ZSORT);
  }

} // end of sortz_HOSTLIB


// =======================================
void zptr_HOSTLIB(void) {

  // construct z-pointers to z-sorted library 
  // so that later we can quickly find the
  // nearest redshift in library.
  //
  // Sep 11, 2012: switch from linear to logz grid
  // Nov 20, 2015: fix syntax bug: fabsf -> fabs  for zdif
  // Dec 29, 2017: speed up; see igal_start and NPAST
  // Sep 20, 2019: compute NZPTR from define params, and alloate IZPTR
  //

  double DZPTR     = (double)DZPTR_HOSTLIB;
  double LOGZRANGE = (double)LOGZRANGE_HOSTLIB ;

  int iz, igal, igal_start, igal_last, NTMP, NSET, NPAST, NZPTR ;
  double ZTRUE, ZSAVE, LOGZ_GRID, Z_GRID, zdif, zdifmin ;
  char fnam[] = "zptr_HOSTLIB" ;

  // ----------- BEGIN -------------

  HOSTLIB.MINiz = 9999999;
  HOSTLIB.MAXiz = -9;
  igal_last=1;

  NZPTR         = (int)(LOGZRANGE/DZPTR); ;
  HOSTLIB.NZPTR = NZPTR;
  HOSTLIB.IZPTR = (int*) malloc ( NZPTR*sizeof(int) );

  /*
  printf(" xxx %s: NZPTR=%d bins %.3f < logz < %.3f (range=%.3f)\n", 
	 fnam, NZPTR, MAXLOGZ_HOSTLIB, MINLOGZ_HOSTLIB, LOGZRANGE);
  */

  for ( iz = 0; iz < NZPTR ; iz++ ) {

    HOSTLIB.IZPTR[iz] = 0 ;

    LOGZ_GRID    = MINLOGZ_HOSTLIB + DZPTR_HOSTLIB * (double)iz ;
    Z_GRID       = pow(10.0,LOGZ_GRID);

    zdifmin = 9999. ;
    ZSAVE   = -9.0  ;
    NSET    = NPAST = NTMP = 0 ;

    if ( iz == 0 ) 
      { igal_start = 0; }
    else
      { igal_start = igal_last ; }

    // find closest lib entry above grid value
    for ( igal = igal_start ; igal < HOSTLIB.NGAL_STORE; igal++ ) {
      ZTRUE = get_ZTRUE_HOSTLIB(igal);

      if ( ZTRUE < Z_GRID ) { continue ; }
      zdif = fabs(Z_GRID - ZTRUE) ;

      NTMP++ ;

      if ( zdif < zdifmin ) {
        zdifmin  = zdif ;
        HOSTLIB.IZPTR[iz] = igal_last = igal ;
	ZSAVE    = ZTRUE ;
	NSET++ ;
      }
      else if ( NSET > 0 ) {
	NPAST++ ;
      }

      if ( NPAST > 10 ) { igal = HOSTLIB.NGAL_STORE; }

    } // end of igal loop

    if ( HOSTLIB.IZPTR[iz] > 0 ) {
      if ( iz < HOSTLIB.MINiz ) { HOSTLIB.MINiz = iz; }
      if ( iz > HOSTLIB.MAXiz ) { HOSTLIB.MAXiz = iz; }
    }

    /*
    if ( iz >= 0 && iz < 5550 ) {
      igal = HOSTLIB.IZPTR[iz] ;
      printf(" zzzzz iz=%3d HOSTLIB.IZPTR[z=%.3f] = %5d (ZLIB=%.3f)\n", 
             iz, Z_GRID, igal,  ZSAVE );
	     fflush(stdout);
    }
    */

  } // end of iz loop


  return ;

} // end of zptr_HOSTLIB

// =======================================
void init_HOSTLIB_WGTMAP(void) {

  // Wgt map is already read;
  // do some inits, allocate memory and make sanity checks.
  //
  // Jan 27 2017: fix bug mallocing GRIDMAP_HOSTLIB_WGT.FUNVAL;
  //              I8p -> I8p*2
  // 
  // Jun 18 2019: if interp_GRIDMAP fails, print more PRE-ABORT info.
  // Jun 25 2019: check GAMMA_GRID option

  int  i, NDIM, ivar, ivar_STORE,NFUN, NROW, istat ;
  int  NGAL, NCHECK, NN, igal, igal_difmax, LDMPWGT, VBOSE ;

  double GAMMA_GRID_MIN = INPUTS.BIASCOR_SALT2GAMMA_GRID[0]; 
  double GAMMA_GRID_MAX = INPUTS.BIASCOR_SALT2GAMMA_GRID[1]; 
  int    USE_GAMMA_GRID = (GAMMA_GRID_MAX > GAMMA_GRID_MIN );  

  long long GALID, GALID_CHECK ;
  int I8  = sizeof(double);

  double
    VAL, VALMIN, VALMAX, WGT, WGTSUM, WGTSUM_LAST, SNMAGSHIFT
    ,VAL_WGTMAP[MXVAR_HOSTLIB]
    ,ZTRUE, ZTRUE_CHECK, ZDIF, WGT_EXACT, WGT_INTERP, WDIF
    ,WDIF_SUM, SQWDIF_SUM, WDIF_AVG, WDIF_RMS, XN, WDIF_MAX, SQTMP
    ,TMPVAL[2]
    ;

  char cvar[40], *varName ;
  char fnam[] = "init_HOSTLIB_WGTMAP" ;

  // --------- BEGIN -----------

  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );
  printf("\t Interpolate WGTMAP for each galaxy ... \n" ) ;
  fflush(stdout);

  NDIM = HOSTLIB_WGTMAP.GRIDMAP.NDIM ;
  NFUN = HOSTLIB_WGTMAP.GRIDMAP.NFUN ;
  NROW = HOSTLIB_WGTMAP.GRIDMAP.NROW ;
  NGAL = HOSTLIB.NGAL_STORE;
  HOSTLIB_WGTMAP.USE_SALT2GAMMA_GRID = USE_GAMMA_GRID;

  if ( USE_GAMMA_GRID ) {
    printf("\t Implement BIASCOR_SALT2GAMMA_GRID: %.2f to %.2f mag\n",
	   GAMMA_GRID_MIN, GAMMA_GRID_MAX);
    fflush(stdout);
  }

  // ------------------------------------------------------
  // allocate memory for wgt, snmagshift and cumulative weight-sum.
  // Note we need double-precision here.
  // Here the  memory is allocated for each GALID.

  HOSTLIB_WGTMAP.WGT        = (double *)malloc(2*HOSTLIB.MALLOCSIZE_D+I8);
  HOSTLIB_WGTMAP.WGTSUM     = (double *)malloc(2*HOSTLIB.MALLOCSIZE_D+I8);
  HOSTLIB_WGTMAP.SNMAGSHIFT = (double *)malloc(2*HOSTLIB.MALLOCSIZE_D+I8);
  WGTSUM_LAST = 0.0 ;

  for ( igal=0; igal < NGAL; igal++ ) {

    HOSTLIB_WGTMAP.WGT[igal]        = 0.0 ;
    HOSTLIB_WGTMAP.WGTSUM[igal]     = 0.0 ;
    HOSTLIB_WGTMAP.SNMAGSHIFT[igal] = 0.0 ;

    GALID  = get_GALID_HOSTLIB(igal);
    ZTRUE  = get_ZTRUE_HOSTLIB(igal);

    if ( NROW == 0 ) {
      WGT        = 1.0 ;
      SNMAGSHIFT = 0.0 ;
      goto WGTSUM ;
    }

    // strip off variables used for weighting

    for ( ivar=0; ivar < NDIM; ivar++ ) {  // WGTMAP variables
      ivar_STORE   = HOSTLIB.IVAR_STORE[ivar];
      VAL          = HOSTLIB.VALUE_ZSORTED[ivar_STORE][igal] ;
      VAL_WGTMAP[ivar] = VAL ;
    }

    istat = interp_GRIDMAP(&HOSTLIB_WGTMAP.GRIDMAP, VAL_WGTMAP, TMPVAL ) ;
    if ( istat != SUCCESS ) {
      printf("\n PRE-ABORT DUMP: \n");
      printf("\t GALID = %lld \n", GALID);
      for ( ivar=0; ivar < NDIM; ivar++ ) {  // WGTMAP variables
	ivar_STORE   = HOSTLIB.IVAR_STORE[ivar];
	varName      = HOSTLIB.VARNAME_STORE[ivar_STORE] ;
	VAL          = HOSTLIB.VALUE_ZSORTED[ivar_STORE][igal] ;
	VALMIN       = HOSTLIB_WGTMAP.GRIDMAP.VALMIN[ivar];
	VALMAX       = HOSTLIB_WGTMAP.GRIDMAP.VALMAX[ivar];
	printf("\t %s = %f  (WGTMAP range: %f to %f)\n", 
	       varName, VAL, VALMIN, VALMAX ); fflush(stdout);
      }
      sprintf(c1err,"Could not interpolate WGTMAP for GALID=%lld .", GALID);
      sprintf(c2err,"interp_GRIDMAP() returned istat = %d", istat);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    WGT        = TMPVAL[0] / HOSTLIB_WGTMAP.WGTMAX ;
    SNMAGSHIFT = TMPVAL[1] ; 

    
  WGTSUM:

    if(USE_GAMMA_GRID) { SNMAGSHIFT = snmagshift_salt2gamma_HOSTLIB(GALID); }

    // local sum
    WGTSUM = WGTSUM_LAST + WGT;

    // print first 3 weights and last wgt
    if ( VBOSE && (igal <= 1 || igal == NGAL-1) ) {
      printf("\t\t WGT(GALID=%lld) = %f -> WGTSUM = %f \n", 
	     GALID, WGT, WGTSUM ); 
      fflush(stdout);
    }
    
    // load global array for each sum
    HOSTLIB_WGTMAP.WGT[igal]        = WGT    ;
    HOSTLIB_WGTMAP.WGTSUM[igal]     = WGTSUM ;
    HOSTLIB_WGTMAP.SNMAGSHIFT[igal] = SNMAGSHIFT ;
    WGTSUM_LAST =  WGTSUM ;

    // store IGAL if this GALID is on the check-list
    for ( i=0; i < HOSTLIB_WGTMAP.NCHECKLIST; i++ ) {
      GALID_CHECK = HOSTLIB_WGTMAP.CHECKLIST_GALID[i];
      ZTRUE_CHECK = HOSTLIB_WGTMAP.CHECKLIST_ZTRUE[i] ;
      ZDIF        = fabs(ZTRUE - ZTRUE_CHECK) ;
      if ( GALID == GALID_CHECK && ZDIF < 2.0E-4 ) 
	{ HOSTLIB_WGTMAP.CHECKLIST_IGAL[i] = igal ; }
    }


    LDMPWGT = ( igal == -9 ) ; //  INPUTS.HOSTLIB_MAXREAD - 10  );
    if ( LDMPWGT ) {
      sprintf(cvar,"%s", HOSTLIB_WGTMAP.VARNAME[0] );
      printf(" xxx GALID=%lld  %s=%6.2f   WGT=%5.3f  WGTSUM=%f \n", 
	     GALID, cvar, VAL_WGTMAP[0], WGT, WGTSUM );
      fflush(stdout);
    }

  } // end if igal loop


  // --------------------------
  // verify interpolated WGTMAP values against optional list of 
  // exact WGT values specified by the WGTMAP_CHECK keys.

  NCHECK = HOSTLIB_WGTMAP.NCHECKLIST; 
  WDIF_SUM = SQWDIF_SUM = WDIF_MAX = 0.0 ;
  igal_difmax = 0;
  NN = 0 ;

  for ( i=0; i <  NCHECK ; i++ ) {
    igal        = HOSTLIB_WGTMAP.CHECKLIST_IGAL[i] ;
    GALID       = HOSTLIB_WGTMAP.CHECKLIST_GALID[i] ;
    WGT_EXACT   = HOSTLIB_WGTMAP.CHECKLIST_WGT[i] ;
    WGT_INTERP  = HOSTLIB_WGTMAP.WGT[igal] ;
    WGT_INTERP *= HOSTLIB_WGTMAP.WGTMAX ; // back to user's WGT definition
    // if ( WGT_EXACT < 0.1 ) { continue ; } // test only
    NN++ ;
    WDIF        = WGT_INTERP - WGT_EXACT ;
    WDIF_SUM   += WDIF ;
    SQWDIF_SUM += (WDIF*WDIF) ;    

    if ( fabs(WDIF) > WDIF_MAX ) 
      { WDIF_MAX = fabs(WDIF);  igal_difmax = igal ; }
  } // end loop over NCHECK


  // print summary of WGTMAP interp-Exact comparison :
  // mean, rms and largest outlier.
  if ( NN > 0 ) {
    XN = (double)NN ;
    WDIF_AVG = WDIF_SUM / XN ;
    SQTMP    = fabs( (SQWDIF_SUM / XN) - (WDIF_AVG*WDIF_AVG) ) ;
    WDIF_RMS = sqrt(SQTMP);

    printf("\t\t WGTMAP CHECK: <Interp-Exact> = %6.4f +- %6.4f(RMS) \n",
	   WDIF_AVG,WDIF_RMS);

    GALID  = get_GALID_HOSTLIB(igal_difmax);
    printf("\t\t WGTMAP CHECK OUTLIER: GALID=%lld  Interp-Exact=%f\n",
	   GALID, WDIF_MAX );
    fflush(stdout);
  }


  //  debugexit("checklist"); // xxxxxxxxx

} // end of init_HOSTLIB_WGTMAP


// ==============================================
double snmagshift_salt2gamma_HOSTLIB(int GALID) {

  // Created Jun 25 2019
  // Randomaly assign gamma (magshift) as follows:
  //    odd  GALID --> SALT2GAMA_GRID_MIN
  //    even GALID --> SALT2GAMA_GRID_MAX
  // Note that this function is used only to create a
  // biasCor sample for BBC/SALT2mu.

  double snmagshift = 0.0 ;
  double GAMMA_GRID_MIN = INPUTS.BIASCOR_SALT2GAMMA_GRID[0]; 
  double GAMMA_GRID_MAX = INPUTS.BIASCOR_SALT2GAMMA_GRID[1]; 

  // -------------- BEGIN ------------

  if ( (GALID%2) == 1 ) 
    { snmagshift = GAMMA_GRID_MIN ; }
  else
    { snmagshift = GAMMA_GRID_MAX ; }

  return(snmagshift);

} // end snmagshift_salt2gamma_HOSTLIB

// =======================================
void init_HOSTLIB_ZPHOTEFF(void) {

  char fnam[] = "init_HOSTLIB_ZPHOTEFF" ;
  char *ptrFile = INPUTS.HOSTLIB_ZPHOTEFF_FILE ;

  HOSTLIB_ZPHOTEFF.NZBIN = 0 ;
  if ( IGNOREFILE(ptrFile) ) { return ; }

  rd2columnFile(ptrFile, MXBIN_ZPHOTEFF, 
		&HOSTLIB_ZPHOTEFF.NZBIN,
		HOSTLIB_ZPHOTEFF.ZTRUE,
		HOSTLIB_ZPHOTEFF.EFF );

  printf("\t Read EFF(zPHOT) vs. ZTRUE (%d bins) from \n\t\t %s\n",
	 HOSTLIB_ZPHOTEFF.NZBIN, ptrFile);
  fflush(stdout);

  if ( HOSTLIB.IVAR_ZPHOT < 0 ) {
    sprintf(c1err,"Cannot implement HOSTLIB_EFF(zPHOT) option");
    sprintf(c2err,"because ZPHOT is not defined in the hostlib.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


} // end of init_HOSTLIB_ZPHOTEFF


// =======================================
void init_GALMAG_HOSTLIB(void) {

  int NR, j, jth, ifilt, ifilt_obs, IVAR, NMAGOBS ;
  char cvar[12], cfilt[2];
  char fnam[] = "init_GALMAG_HOSTLIB" ;
  double Rmax, TH, THbin ; 

  // --------------- BEGIN --------------

  // always check for gal mags. Store IVAR for each observer-mag
  NMAGOBS = 0 ;
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];

    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] ); 
    magkey_HOSTLIB(ifilt_obs,cvar); // returns cvar

    IVAR = IVAR_HOSTLIB(cvar,0) ;
    if ( IVAR > 0 ) {
      NMAGOBS++ ;
      HOSTLIB.IVAR_MAGOBS[ifilt_obs] = IVAR ;
      strcat(HOSTLIB.filterList,cfilt) ;
    }
  }

  HOSTLIB.NFILT_MAGOBS = NMAGOBS ;  // store in global

  // first check option to compute host mags.

  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG) == 0 ) { 
    printf("\t Skip Host Galaxy Mags and Noise. \n");  fflush(stdout);
    return ; 
  }
  else {
    printf("\t Compute Host Galaxy Mags and Noise. \n"); fflush(stdout);
  }
  
  
  // make sure we have a galaxy profile 
  if ( SERSIC_PROFILE.NDEF <= 0 ) {
    sprintf(c1err,"%s", "Galaxy (Sersic) profile NOT defined" );
    sprintf(c2err,"%s", "==> cannot compute galaxy mag at SN position.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  NR=0;  HOSTLIB.Aperture_PSFSIG[NR] = 0.0 ; 
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 0.03/2.35; // Added 9.28.2015 for WFIRST
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 0.07/2.35; // Added 5.22.2015 for WFIRST
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 0.10/2.35; // idem
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 0.20/2.35; // idem
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 0.40/2.35; // PSF-sigma, arcsec
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 0.80/2.35; 
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 1.30/2.35;
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 2.10/2.35;
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = 5.00/2.35;


  if ( NR != NMAGPSF_HOSTLIB ) {
    sprintf(c1err,"%d defined PSF values", NR );
    sprintf(c2err,"but expected NMAGPSF_HOSTLIB = %d", 
	    NMAGPSF_HOSTLIB );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  
  if ( NMAGOBS <= 0 ) {
    sprintf(c1err, "%s" , "Galaxy mags under SN requested, but no mags are");
    sprintf(c2err, "given in HOSTLIB. Need [filt]%s keys.", 
	    HOSTLIB_MAGOBS_SUFFIX );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // -------------------------------------
  // Set aperture Radii = 2*PSF since effective area is 
  // 4*PI*R^2 = PI*( 2*PSF)^2
  for ( j=0; j <= NMAGPSF_HOSTLIB ; j++ ) { 
    HOSTLIB.Aperture_Radius[j] = 2.0 * HOSTLIB.Aperture_PSFSIG[j];    
  }

  // integrate out to 2 x maximum aperture radius
  Rmax = 2. * HOSTLIB.Aperture_Radius[NMAGPSF_HOSTLIB];
  HOSTLIB.Aperture_Rmax  = Rmax ;
  HOSTLIB.Aperture_Rbin  = Rmax  / (double)NRBIN_GALMAG  ; 
  HOSTLIB.Aperture_THbin = TWOPI / (double)NTHBIN_GALMAG ; 


  // pre-compute each cos(TH) and sin(TH)
  THbin = HOSTLIB.Aperture_THbin ;
  jth   = 0;
  for ( TH = 0.0; TH < .99999*TWOPI; TH += THbin ) {
    HOSTLIB.Aperture_cosTH[jth] = cos(TH);
    HOSTLIB.Aperture_sinTH[jth] = sin(TH);
    jth++ ;
  }


  // init 2d Gaussian integrals
  init_Gauss2d_Overlap();

} // init_GALMAG_HOSTLIB


// =====================================
void init_Gauss2d_Overlap(void) {
  
  // Read and store table of 2d Gaussian integrals that
  // overlap with a circle (aperture) of radius = 1.
  //
  // Oct 27 2014: switch to refactored fitres-read (SNTABLE_XXX)
  // May 02 2019: remove TABLEFILE_CLOSE() call since READ_EXEC does it.

  char 
     TABLENAME[] = "Gass2D"
    ,cvar[20]
    ,Gauss2dFile[MXPATHLEN] 
    ,GaussVAR[NVAR_Gauss2d][20] = { "R", "SIG", "INTEG2D" }
  ;

  char fnam[] = "init_Gauss2d_Overlap" ;

  int IVAR, istat, i, Nr, Ns, VBOSE, IFILETYPE ;
  double 
    *ptr
    ,r, s, Integ2d 
    ,rmin, rmax, rbin, rlast
    ,smin, smax, sbin, slast
    ;

  // ------------- BEGIN ----------

  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );

  // read and store Gauss2d table

  sprintf(Gauss2dFile,"%s", FILENAME_Gauss2d);
  ENVreplace(Gauss2dFile,fnam,1);

  TABLEFILE_INIT();
  IFILETYPE = TABLEFILE_OPEN(Gauss2dFile, "read");
  SNTABLE_READPREP(IFILETYPE,TABLENAME);

  for ( IVAR = 0; IVAR < NVAR_Gauss2d; IVAR++ ) {

    sprintf(cvar,"%s", GaussVAR[IVAR] );
    ptr   = &HOSTLIB.Gauss2dTable[IVAR][0];
    istat = SNTABLE_READPREP_VARDEF(cvar, ptr, MXGauss2dTable, VBOSE );
    
    if ( istat  < 0 ) {
      sprintf(c1err,"Could not find '%s' in table-file", cvar);
      sprintf(c2err,"%s", Gauss2dFile);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }
  }


  // Read the Gauss2d table
  HOSTLIB.NGauss2d = SNTABLE_READ_EXEC();
  // xxx obsolete TABLEFILE_CLOSE(Gauss2dFile);

  // decode binning for easy lookup
  rmin = 1.0E9; rmax = -1.0E9 ; rlast = -999. ;  rbin = 0.0 ;
  smin = 1.0E9; smax = -1.0E9 ; slast = -999. ;  sbin = 0.0 ;

  for ( i=0; i < HOSTLIB.NGauss2d; i++ ) {
    r       = HOSTLIB.Gauss2dTable[0][i];  // radius
    s       = HOSTLIB.Gauss2dTable[1][i];  // sigma
    Integ2d = HOSTLIB.Gauss2dTable[2][i];

    if ( r < rmin ) { rmin = r; }
    if ( r > rmax ) { rmax = r; }
    if ( s < smin ) { smin = s; }
    if ( s > smax ) { smax = s; }

    if ( r != rlast && i > 0 ) { rbin = r - rlast ; }
    if ( s != slast && i > 0 ) { sbin = s - slast ; }

    rlast = r ; slast = s;
    //  printf("\t xxxx r=%f  sig=%f  Integ2d=%f \n", r, sig, Integ2d);
  }

  // determine number of R and SIG bins
  Nr = Ns = 1; 
  if ( rbin>0.0 ) { Nr = (int)((rmax - rmin + 0.00001)/rbin) + 1 ; }
  if ( sbin>0.0 ) { Ns = (int)((smax - smin + 0.00001)/sbin) + 1 ; }

  
  if ( VBOSE ) {
    printf("\t Gauss2d Table radius : %4.2f - %4.2f (%d bins of %4.2f) \n", 
	   rmin, rmax, Nr, rbin );
    printf("\t Gauss2d Table sigma  : %4.2f - %4.2f (%d bins of %4.2f) \n", 
	   smin, smax, Ns, sbin );
    fflush(stdout);
  }

  // store bin info in global struct

  HOSTLIB.Gauss2dRadius[0] = rbin ;
  HOSTLIB.Gauss2dRadius[1] = rmin ;
  HOSTLIB.Gauss2dRadius[2] = rmax ;
  HOSTLIB.NBIN_Gauss2dRadius  = Nr;

  HOSTLIB.Gauss2dSigma[0] = sbin ;
  HOSTLIB.Gauss2dSigma[1] = smin ;
  HOSTLIB.Gauss2dSigma[2] = smax ;
  HOSTLIB.NBIN_Gauss2dSigma  = Ns;

  return ;

} // end of init_Gauss2d_Overlap

// ======================================
double interp_GALMAG_HOSTLIB(int ifilt_obs, double PSF ) {
 
  //
  // Return interpolated GALMAG for this 
  // obs-filter index and PSF-sigma (arcsec).
  // Interpolate pre-computed GALMAG vs. PSF.
  // GALMAG corresponds to flux contained in 4*PI*PSF^2 aperture.
  //
  // Inputs:
  //  ifilt_obs = absolute obs-filter index
  //  PSF       = sigma(PSF) in arcsec

  int NPSF ;
  double GALMAG, PSFmin, PSFmax ;
  double *PTRGRID_GALMAG, *PTRGRID_PSF ;
  char fnam[] = "interp_GALMAG_HOSTLIB" ;

  // ------------- BEGIN -------------
  
  NPSF = NMAGPSF_HOSTLIB ;
  PSFmin = HOSTLIB.Aperture_PSFSIG[1] ;
  PSFmax = HOSTLIB.Aperture_PSFSIG[NPSF] ;

  if ( PSF < PSFmin ||  PSF > PSFmax ) {
    char cfilt[2];
    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );
    sprintf(c1err,"%s-PSFSIG=%f asec outside valid range (%5.2f - %5.2f)",
	    cfilt, PSF, PSFmin, PSFmax );
    sprintf(c2err, "Cannot extrapolate GALMAG. SIMLIB_ID=%d", 
	    GENLC.SIMLIB_ID) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // note that the zero'th element is total aperture,
  // so ignore it for interpolatio.
  PTRGRID_PSF     = &HOSTLIB.Aperture_PSFSIG[1] ;
  PTRGRID_GALMAG  = &SNHOSTGAL.GALMAG[ifilt_obs][1] ;

  GALMAG = interp_1DFUN (1,PSF,NPSF, PTRGRID_PSF, PTRGRID_GALMAG, "GALMAG");

  return(GALMAG) ;

}  // interp_GALMAG_HOSTLIB



// ======================================
void magkey_HOSTLIB(int ifilt_obs, char *key) {
  sprintf(key,"%c%s", FILTERSTRING[ifilt_obs], HOSTLIB_MAGOBS_SUFFIX );
} // end of magkey_HOSTLIB



// =======================================
void init_Sersic_VARNAMES(void) {

  // check that every defined a#,b# also has a defined
  // Sersic index, and vice-versa; abort on any problem.
  //
  // Also load index arrays
  //  SERSIC.IVAR_a[j1] = HOSTLIB.IVAR_a[j2] ;
  //  SERSIC.IVAR_b[j1] = HOSTLIB.IVAR_b[j2] ;
  //
  // where j1 is the sparse SERSIC-array index (1-SERSIC.NDEF)
  // and j2 is the absolute 0-MXSERSIC integer key. Note that 
  // neither j1 or j2 is the physical SERSIC.n[j1] index.
  //
  // The index-matching is done here after read_head_HOSTLIB()
  // so that the VARNAMES and SERSIC keys can go in any order
  // in the HOSTLIB
  //

  int 
    NDEF, NFIX, j, ifix, NFIX_MATCH, NERR
    , IVAR_a, IVAR_b, IVAR_w, IVAR_n    
    ;

  double FIXn;
  char anam[12], bnam[12], wnam[12], nnam[12], tmpnam[12] ;
  char fnam[] = "init_Sersic_VARNAMES" ;

  // ----------- BEGIN ------------

  NDEF = 0;

  // check which VARNAMES correspond to Sersic profiles,
  // and store list. Search all possible profile names
  // and check which ones are on the VARNAMES list.

  for ( j=0; j < MXSERSIC_HOSTLIB; j++ ) {

    Sersic_names(j, anam, bnam, wnam, nnam );

    IVAR_a = IVAR_HOSTLIB(anam,0);
    IVAR_b = IVAR_HOSTLIB(bnam,0);
    IVAR_w = IVAR_HOSTLIB(wnam,0);
    IVAR_n = IVAR_HOSTLIB(nnam,0);

    // make sure that if 'ai_Sersic' is defined, then 'bi_sersic' 
    // is defined and vice-versa. Check Sersic index and weight later.
    if ( IVAR_a > 0 && IVAR_b < 0 ) {
      sprintf(c1err,"%s is defined but %s is NOT ?", anam, bnam);
      sprintf(c2err,"%s", "Check VARNAMES in HOSTLIB.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    if ( IVAR_a < 0 && IVAR_b > 0 ) {
      sprintf(c1err,"%s is defined but %s is NOT ?", bnam, anam);
      sprintf(c2err,"%s", "Check VARNAMES in HOSTLIB.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // load sersic list if both a & b are specified.
    if ( IVAR_a > 0 && IVAR_b > 0 ) {
      SERSIC_PROFILE.IVAR_a[NDEF] = IVAR_a ;
      SERSIC_PROFILE.IVAR_b[NDEF] = IVAR_b ;
      SERSIC_PROFILE.IVAR_w[NDEF] = IVAR_w ;
      SERSIC_PROFILE.IVAR_n[NDEF] = IVAR_n ;

      sprintf(SERSIC_PROFILE.VARNAME_a[NDEF], "%s", anam);
      sprintf(SERSIC_PROFILE.VARNAME_b[NDEF], "%s", bnam);
      sprintf(SERSIC_PROFILE.VARNAME_w[NDEF], "%s", wnam);
      sprintf(SERSIC_PROFILE.VARNAME_n[NDEF], "%s", nnam);   
      NDEF++ ;   
    }
  } // end of j-loop

  SERSIC_PROFILE.NDEF = NDEF ;
  if ( NDEF == 0 ) { return ; }


  // check list of FIXED Sersic indices and make sure
  // that they all match an a#,b# pair.
  NFIX_MATCH = 0;
  NFIX   = SERSIC_PROFILE.NFIX;

 
  for ( ifix=0; ifix < NFIX; ifix++ ) {
    sprintf(tmpnam, "%s", SERSIC_PROFILE.FIX_NAME[ifix]  ) ;

    for ( j=0; j < NDEF; j++ ) {
      sprintf(nnam, "%s", SERSIC_PROFILE.VARNAME_n[j] ) ;

      if ( strcmp(nnam,tmpnam)==0 ) {
	SERSIC_PROFILE.FIXn[j] = SERSIC_PROFILE.FIX_VALUE[ifix]; 
	NFIX_MATCH++ ;
      }
    }
  }


  // abort on any undefined Sersic indices
  if ( NFIX_MATCH != NFIX ) {
    sprintf(c1err,"%d FIXED Sersic indices defined (n#_Sersic)", NFIX);
    sprintf(c2err,"but %d matched a,b profiles.", NFIX_MATCH );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // now make sure that every Sersic component has a defined index;
  // either specified for each host (IVAR_n>0) or globally fixed (FIXn).
  // Abort if defined  for both.

  NERR = 0;
  for ( j=0; j < NDEF; j++ ) {
    IVAR_n = SERSIC_PROFILE.IVAR_n[j] ;
    FIXn   = SERSIC_PROFILE.FIXn[j] ;

    if ( IVAR_n < 0 && FIXn < 0.0 ) {
      NERR++ ;
      printf(" ERROR: Sersic index '%s' is NOT defined for %s/%s \n"
	,SERSIC_PROFILE.VARNAME_n[j]
	,SERSIC_PROFILE.VARNAME_a[j]
	,SERSIC_PROFILE.VARNAME_b[j] );
      fflush(stdout);
    }

    if ( IVAR_n >= 0 && FIXn >= 0.0 ) {
      sprintf(c1err,"Sersic index '%s' is defined twice.", 
	      SERSIC_PROFILE.VARNAME_n[j] );
      sprintf(c2err,"Only 1 definition allowed.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  }

  if ( NERR > 0 ) {
    sprintf(c1err,"%d missing Sersic indices !!!", NERR);
    sprintf(c2err,"Check VARNAMES keys in HOSTLIB.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

} // end of init_Sersic_VARNAMES


// =======================================
void init_Sersic_HOSTLIB(void) {

  // Driver routine to init Sersic integrals
  // need for relative flux calculations/

  int  j, NBIN, ir, VBOSE    ;
  double  Rmax, Rmin, logRbin, logRmin, logRmax, logR, dif, xmem  ;
  //  char fnam[] = "init_Sersic_HOSTLIB" ;

  // --------------- BEGIN ------------

  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );

  // setup reduced R/Re grid in space of log10(R/Re)
  NBIN = NBIN_RADIUS_SERSIC ; 
  Rmax = MAXRADIUS_SERSIC ; 
  Rmin = MINRADIUS_SERSIC ; 
  logRmin = log10(Rmin);
  logRmax = log10(Rmax);
  logRbin = (logRmax - logRmin) / (double)NBIN ;

  SERSIC_TABLE.NBIN_reduced      = NBIN ;
  SERSIC_TABLE.reduced_logRmin   = logRmin ; 
  SERSIC_TABLE.reduced_logRmax   = logRmax ; 
  SERSIC_TABLE.reduced_logRbin   = logRbin ; 
  SERSIC_TABLE.reduced_logR      = (double *)malloc(8*NBIN+8);


  // compute R/Re grid 
  for ( ir=0; ir <= NBIN; ir++ ) {
    logR = logRmin + logRbin * (double)ir ;
    SERSIC_TABLE.reduced_logR[ir]  = logR ;
  }


  // setup index-table with equal binsize for 1/n
  NBIN = NSERSIC_TABLE ;
  SERSIC_TABLE.INVINDEX_MIN = 1.0/SERSIC_INDEX_MAX ; // min 1/n ==> 1/n_max 
  SERSIC_TABLE.INVINDEX_MAX = 1.0/SERSIC_INDEX_MIN ; // max 1/n ==> 1/n_min 

  dif  = SERSIC_TABLE.INVINDEX_MAX - SERSIC_TABLE.INVINDEX_MIN ;
  SERSIC_TABLE.INVINDEX_BIN = dif / ((double)NBIN - 1.0 ) ; // ??


  // read b_n(n) from ascii file (Jun 2015)
  read_Sersic_bn();

  // init Sersic integral table for each 1/n value
  for ( j=0; j < NSERSIC_TABLE; j++ )
    { init_Sersic_integrals(j); }

  
  // print stuff if VBOSE mode is set.
  if ( VBOSE ) {
    printf("\t Init %d log10(R/Re) bins from %6.3f to %6.3f \n",
	   NBIN+1, logRmin, logRmax );
    printf("\t Init %d Sersic integral tables spanning n = %4.2f - %4.2f",
	   NBIN, SERSIC_INDEX_MIN , SERSIC_INDEX_MAX );
    fflush(stdout);
    xmem = 1.0E-3 * (double)SERSIC_TABLE.TABLEMEMORY ;
    printf(" (%4.1f kB) \n", xmem );
    fflush(stdout);

    if ( fabs(INPUTS.HOSTLIB_SERSIC_SCALE-1.0)>1.0E-5  ) {
      printf("\t Sersic(a0,b0) scale:  %.3f \n", 
	     INPUTS.HOSTLIB_SERSIC_SCALE);
    }
  }

  // optional test of Seric-integral interpolation via 1/n
  //  test_Sersic_interp(); // manually comment this in or out.

  return;

} // end of init_Sersic_HOSTLIB


// ==========================
void read_Sersic_bn(void) {

  int  MXROW, Nrow ;
  char fileName[MXPATHLEN];
  char fnam[] = "read_Sersic_bn" ;

  // ------------------ BEGIN -----------------

  sprintf(fileName, "%s", FILENAME_Sersic_bn );
  ENVreplace(fileName,fnam,1);

  MXROW = MXBIN_SERSIC_bn ;
  rd2columnFile(fileName, MXROW, 
		&Nrow,     // return number of rows in file
		SERSIC_TABLE.grid_n, SERSIC_TABLE.grid_bn ); // output

  SERSIC_TABLE.Ngrid_bn = Nrow ;

  printf("\t Read %d rows of b_n(n) table from \n\t   %s\n", 
	 Nrow, fileName);

} // end of read_Sersic_bn

// =====================================
void test_Sersic_interp(void) {

  // Loop over every other 'j' and test interpolation
  // at skipped index.

  int j0, j1, j2, integ, half ;
  double SINT0, SINT1, SINT2, SINTERP1, n1, logR, R; 
  //   char fnam[] =  "test_Sersic_interp" ;

  // ---------- BEGIN -------------

  for  ( j0=0; j0 < NSERSIC_TABLE-2 ; j0 += 2 ) {
    j1 = j0 + 1 ;
    j2 = j0 + 2 ;
    //    printf(" interp bins %d-%d to test bin %d \n", j0, j2, j1 ); 

    n1 = SERSIC_TABLE.n[j1] ; // index to interpolate

    printf(" Test Sersic integral for n = %f (j0,j1,j2=%d,%d,%d)\n", 
	   n1, j0, j1, j2 );
    fflush(stdout);
    half = SERSIC_TABLE.BIN_HALFINTEGRAL[j1];

    for ( integ=5; integ <= 2*half; integ += 20 ) {
      logR     =   *(SERSIC_TABLE.reduced_logR  + integ)  ;
      SINT0    =   *(SERSIC_TABLE.INTEG_CUM[j0] + integ)  ;
      SINT1    =   *(SERSIC_TABLE.INTEG_CUM[j1] + integ)  ;
      SINT2    =   *(SERSIC_TABLE.INTEG_CUM[j2] + integ)  ;
      SINTERP1 = 0.5*(SINT0 + SINT2);

      R = pow(10.0,logR);
      printf("\t R/Re = %6.3f  Sint(interp/exact) = %f/%f = %f \n",
	     R, SINTERP1, SINT1, SINTERP1/SINT1 );
      fflush(stdout);
    }

  } // end of j0 loop

  debugexit("Done  testing Seric 1/n-interpolation");

} // end of test_Sersic_interp



// ==================================
void get_Sersic_info(int IGAL) {

  // Fill SERSIC info a, b, n(and bn), w, wsum.
  // Require at least NDEF-1 weights to be defined;
  // the last weight is simply 1 = sum(other wgts).
  // If NDEF=1 then no weights are required.
  //
  // Sep 9, 2012: abort if sersic index is outside valid range.
  //

  int j, NDEF, NWGT, j_nowgt, IVAR_a, IVAR_b, IVAR_w, IVAR_n ;
  double WGT, WGTSUM, WTOT, n, wsum_last ;
  char fnam[] = "get_Sersic_info" ;

  // -------------- BEGIN -------------

  NDEF = SERSIC_PROFILE.NDEF; 
  SNHOSTGAL.SERSIC_w[0]    = 0.0 ;
  SNHOSTGAL.SERSIC_wsum[0] = 0.0 ;

  for ( j=0; j < NDEF; j++ ) {
    IVAR_a = SERSIC_PROFILE.IVAR_a[j] ;
    IVAR_b = SERSIC_PROFILE.IVAR_b[j] ;
    IVAR_n = SERSIC_PROFILE.IVAR_n[j] ;

    if ( IVAR_n >= 0 ) 
      { n = HOSTLIB.VALUE_ZSORTED[IVAR_n][IGAL] ; }
    else
      { n = SERSIC_PROFILE.FIXn[j] ; }

    SNHOSTGAL.SERSIC_a[j]  = HOSTLIB.VALUE_ZSORTED[IVAR_a][IGAL] ; 
    SNHOSTGAL.SERSIC_b[j]  = HOSTLIB.VALUE_ZSORTED[IVAR_b][IGAL] ; 
    SNHOSTGAL.SERSIC_n[j]  = n ;
    SNHOSTGAL.SERSIC_bn[j] = get_Sersic_bn(n);

    // apply user-scale on size (Mar 28 2018)
    SNHOSTGAL.SERSIC_a[j] *= INPUTS.HOSTLIB_SERSIC_SCALE ;
    SNHOSTGAL.SERSIC_b[j] *= INPUTS.HOSTLIB_SERSIC_SCALE ;

    if ( n < SERSIC_INDEX_MIN || n > SERSIC_INDEX_MAX ) {
      sprintf(c1err,"Sersic index=%f outside valid range (%5.2f-%5.2f)",
	      n, SERSIC_INDEX_MIN, SERSIC_INDEX_MAX ) ;
      sprintf(c2err,"Check GALID = %lld  ZTRUE=%f ", 
	      SNHOSTGAL.GALID, SNHOSTGAL.ZTRUE );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
  }

  // Now the weights:
  // check how many weights are defined, and track the sum.
  
  NWGT = 0;  IVAR_w=0;  j_nowgt = -9 ;     WGTSUM = WGT = 0.0 ;
  for ( j=0; j < NDEF; j++ ) {
    IVAR_w = SERSIC_PROFILE.IVAR_w[j] ;
    if ( IVAR_w > 0 ) { 
      NWGT++ ;
      WGT     = HOSTLIB.VALUE_ZSORTED[IVAR_w][IGAL] ;
      WGTSUM += WGT;
      SNHOSTGAL.SERSIC_w[j]  = WGT ;
    }
    else { j_nowgt = j ; }
  }


  if ( NWGT == NDEF-1 ) {
    SNHOSTGAL.SERSIC_w[j_nowgt]  = 1.0 - WGTSUM ;
  }

  // check debug option for fixed weight with 2 profiles
  if ( NDEF == 2 && DEBUG_WGTFLUX2 > 0.00001 ) {
    SNHOSTGAL.SERSIC_w[0]  = 1.0 - DEBUG_WGTFLUX2 ;
    SNHOSTGAL.SERSIC_w[1]  = DEBUG_WGTFLUX2 ;
    goto WGTSUM ;
  }


  if ( NWGT < NDEF - 1 ) {
    sprintf(c1err,"%s", "Inadequate Sersic weights.");
    sprintf(c2err,"%d Sersic terms, but only %d weights defined.",
	    NDEF, NWGT );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

 WGTSUM:

  // fill cumulative WGTFLUX_SUM array
  wsum_last = 0.0 ;
  for ( j=0; j < NDEF; j++ ) {
    WGT = SNHOSTGAL.SERSIC_w[j] ;
    // xxx bug SNHOSTGAL.SERSIC_wsum[j] = WGT + SNHOSTGAL.SERSIC_w[j-1] ;
    SNHOSTGAL.SERSIC_wsum[j] = WGT + wsum_last; // bug fix
    wsum_last = SNHOSTGAL.SERSIC_wsum[j] ; 
  }

  // finally check that sum of weights are one
  WTOT = SNHOSTGAL.SERSIC_wsum[NDEF-1] ;
  if ( fabs(WTOT-1.0) > 0.0001 ) {
    sprintf(c1err,"Sum of Sersic weights = %f", WTOT);
    sprintf(c2err,"%s", "Check values of w1, w2 ...");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

} // end of get_Sersic_info


// =======================================
void init_Sersic_integrals(int j) {

  // Mar 14, 2011 R.Kessler
  //
  // for Sersic profile with sparse index 'j', compute integral 
  // table as a function of reduced radius R/Re
  //
  // The integrated flux contined in r < R is
  //      Flux(R) = Re^2 * 2pi * integral[ x*exp[-bn*(x^(1/n)-1)] ]dx
  //
  // where x = R/Re is dimensionless. The integral above
  // is tabulated and stored up to SERSIC.reduced_Rmax.
  //
  // Input 'j' is a lookup-table index from 1 to NSERSIC_TABLE
  // (j is NOT the index over hostlib-sersic profiles)

  int    NBIN, ir, MEM ;
  double 
    inv_n, xj
    ,logRbin, logRmin, logRmax, logR0, logR1, R0, R1, R
    ,F1, F0, FTMP, FSUM, FTOT, FHALF, FSLP
    ,bn, arg, renorm, difmin
    ;

  char fnam[] = "init_Sersic_integrals";

  // ----------- BEGIN --------------


  // allocate memory for this cumulative integral table.  
  MEM  = 8*SERSIC_TABLE.NBIN_reduced + 8;
  SERSIC_TABLE.TABLEMEMORY    += MEM ;
  SERSIC_TABLE.INTEG_CUM[j]    = (double *)malloc(MEM);
  SERSIC_TABLE.INTEG_CUM[j][0] = 0.0 ; 


  xj    =  (double)j  ;
  inv_n = SERSIC_TABLE.INVINDEX_MIN + SERSIC_TABLE.INVINDEX_BIN * xj;
  SERSIC_TABLE.inv_n[j]  = inv_n ;
  SERSIC_TABLE.n[j]      = 1./inv_n ;

  //  bn   = 1.9992 * SERSIC.n[j] - 0.3271 ; // not so good for small 'n'
  bn   = get_Sersic_bn(SERSIC_TABLE.n[j]);
  SERSIC_TABLE.bn[j] = bn ; // store in global array

  logRbin = SERSIC_TABLE.reduced_logRbin ;
  logRmin = SERSIC_TABLE.reduced_logRmin ;
  logRmax = SERSIC_TABLE.reduced_logRmax ;
  NBIN    = SERSIC_TABLE.NBIN_reduced ;
  
  /*
  printf("(%2d) Init Sersic Integ for: 1/n=1/%5.3f = %5.3f   bn=%6.3f\n",
	 j, SERSIC_TABLE.n[j], SERSIC_TABLE.inv_n[j], SERSIC_TABLE.bn[j]  );
  */


  FSUM = 0.0 ;

  for ( ir=1; ir <= NBIN; ir++ ) {
    logR1 = logRmin + logRbin * (double)ir ;  // log10 of reduced radius, R/Re
    logR0 = logR1 - logRbin ;

    R0 = pow(TEN,logR0) ;
    R1 = pow(TEN,logR1) ;

    //compute r*I(R) at start and end of bin
    arg = pow(R0,inv_n) - 1.0 ;    F0 = R0 * exp(-bn*arg) ;
    arg = pow(R1,inv_n) - 1.0 ;    F1 = R1 * exp(-bn*arg) ;

    // do trapezoidal integration 
    FTMP  = 0.5 * (R1-R0) * (F0 + F1) ;
    FSUM += FTMP ;

    // store cumulative integral
    SERSIC_TABLE.INTEG_CUM[j][ir] = FSUM ;
  }


  // tack on 2pi Jacobian factor for the total integral  
  FTOT = SERSIC_TABLE.INTEG_CUM[j][NBIN];
  SERSIC_TABLE.INTEG_SUM[j] = (FTOT * TWOPI) ;  // 2pi is for azim. integral

  if ( FTOT <= 0.0 ) {
    sprintf(c1err,"Sersic integral = %f", FTOT );
    sprintf(c2err,"for SERSIC.n[%d] = %5.2f", j, SERSIC_TABLE.n[j] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }
 
  // re-normalize integrals so that max integral = 1.
  // Also find and store R/Re bin with half-integral.


  renorm = 1. / FTOT ;
  SERSIC_TABLE.BIN_HALFINTEGRAL[j] = -9 ;

  difmin = 1.0E30 ;

  for ( ir=1; ir <= NBIN; ir++ ) {

    FTMP  = SERSIC_TABLE.INTEG_CUM[j][ir] ;
    FTMP *= renorm ;
    SERSIC_TABLE.INTEG_CUM[j][ir]  = FTMP ;
    
    // check for bin where integral is closest to 1/2 the total,
    // without going over.
    if ( FTMP < 0.500 ) {
      SERSIC_TABLE.BIN_HALFINTEGRAL[j] = ir ; 
    }
  }   // end ir loop
    

  // now check that half-integral corresponds to r=1
  ir    = SERSIC_TABLE.BIN_HALFINTEGRAL[j] ;

  if ( ir <= 0 ) {
    sprintf(c1err,"%s", "Could not find R/Re bin with half-integral.");
    sprintf(c2err,"Something is messed up for SERSIC.n[%d]= %5.2f",
	    j, SERSIC_TABLE.n[j] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // interpolate integral
  logR0  = logRmin + logRbin * (double)ir ;  // reduced radius, R/Re  
  logR1  = logRmin + logRbin * (double)(ir+1);
  R0     = pow(TEN,logR0);
  R1     = pow(TEN,logR1);
  F0 = SERSIC_TABLE.INTEG_CUM[j][ir] ;
  F1 = SERSIC_TABLE.INTEG_CUM[j][ir+1] ;

  FSLP  = (0.5-F0)/(F1-F0) ;
  R     = R0 + FSLP*(R1-R0);
  FHALF = F0 + FSLP*(F1-F0); // idiot check: should be 0.5

  /*
  printf(" xxx %3.3d : n=%6.3f    %0.3f<FHALF<%.3f   R/Re = %.3f \n", 
	 j, SERSIC_TABLE.n[j], F0, F1, R);
  */

  if ( fabs(R-1.0) > 0.1 ) {
    sprintf(c1err,"integral=%f  at  R/Re=%f", FHALF, R);
    sprintf(c2err,"but expect F=0.5 at R/Re=1  (n=%5.3f)", 
	    SERSIC_TABLE.n[j] ); 
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


} // end of init_Sersic_integrals


// =======================================
double get_Sersic_bn(double n) {

  // Returns Sersic b_n term .
  // For n > 0.36 use approx from Ciotti & Bertin 1999.
  // For n < 0.36 use exact calc from R.Gupta.
  //
  // Jun 2015: replace pow(n,arg) with n*n*...
  //           Fix bug: 20690717750 -> 30690717750
  //
  double bn, n2, n3, n4 ;
  char fnam[] = "get_Sersic_bn" ;

  // -------------- BEGIN ------------

  if ( n > 0.36 ) {
    n2 = n*n;
    n3 = n2*n;
    n4 = n2*n2;
    
    //   approx from Ciotti & Bertin 1999.
    bn  = 0.0 ;
    bn += 2.0 * n - 1./3. ;
    bn += 4. / (405.*n) ;
    bn += 46./( 25515.0 * n2 );
    bn += 131. / ( 1148175. * n3 );
    bn += 2194697 / ( 30690717750. * n4 ); 
  }
  else {
    // interp exact calc
    int opt = 1; // 1=linear interp; 2=quadratic
    bn = interp_1DFUN(opt, n, 
		      SERSIC_TABLE.Ngrid_bn,
		      SERSIC_TABLE.grid_n,
		      SERSIC_TABLE.grid_bn,
		      fnam );
  }

  return bn;
} // end  of bn_CiottiBertin


// =============================================
void  init_SAMEHOST(void) {

  // July 2015
  // init pointers for storing list of PEAKMJD on each host.

  int MINDAYSEP, USEONCE, i2size, igal, NGAL ;
  char fnam[] = "init_SAMEHOST" ;

  // -----------------------------------------

  SAMEHOST.REUSE_FLAG = 0 ;

  // always allocate NUSE for each host
  i2size = sizeof(unsigned short);
  NGAL = HOSTLIB.NGAL_STORE ;
  SAMEHOST.NUSE =  (unsigned short *)malloc ( (NGAL+10) * i2size) ; 
  for(igal=0; igal < NGAL+10; igal++ ) 
    { SAMEHOST.NUSE[igal] = 0; }

  USEONCE   = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEONCE );
  MINDAYSEP = INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL ;

  // trap user-input mistake
  if ( USEONCE && MINDAYSEP < 9999 ) {
    sprintf(c1err ,"Cannot re-use host after %d days", MINDAYSEP );
    sprintf(c2err ,"because USEONCE flag is set (%d-bit of HOSTLIB_MSKOPT)", 
	    HOSTLIB_MSKOPT_USEONCE) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // if useOnce-bit it set, then skip init
  if ( USEONCE ) { 
    printf("\t %s: Use a host only once.\n", fnam );
    return ; 
  }

  if ( USEONCE == 0 && MINDAYSEP > 9999 ) {
    SAMEHOST.REUSE_FLAG = 1 ;  // re-use arbitrary number of times
    printf("\t %s: Re-use host arbitrary number of times.\n", fnam);
    return ;
  }

  // if MINDAYSEP is too big, skip this option
  if ( MINDAYSEP > 9999 ) { return ; }

  // -------- if we get here, use this feature ------------

  int MJD0, MJD1, TTOT;

  MJD0 = (int)INPUTS.GENRANGE_PEAKMJD[0];
  MJD1 = (int)INPUTS.GENRANGE_PEAKMJD[1];
  TTOT = MJD1 - MJD0 ;
  SAMEHOST.PEAKMJD_STORE_OFFSET = MJD0;

  // make sure that max TTOT fits in 2-byte integer
  if ( TTOT > 62000 ) {
    sprintf(c1err ,"PEAKMJD range spans %d days", TTOT);
    sprintf(c2err ,"too big for PEAKDAY_STORE 2-byte storage.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // allocate
  SAMEHOST.PEAKDAY_STORE = 
    (unsigned short **)malloc ( (NGAL+10) * sizeof(unsigned short*) ) ; 

  // set flag for PEAKMJD-sep option
  SAMEHOST.REUSE_FLAG = 2 ;

  // ---------
  printf("\t %s: Allocate PEAKDAY_STORE array to re-use host after %d days\n",
	 fnam, MINDAYSEP);
  fflush(stdout);

} // end of  init_SAMEGAL_HOSTLIB

// =============================================
void readme_HOSTLIB(void) {

  // prepare comments for README file.

  int    NTMP, ivar, NVAR, NROW, LSN2GAL, LGALMAG, j, igal ;
  long long GALID ;
  double RAD, WGT, fixran ;
  char *cptr,  copt[40],  ctmp[20], ZNAME[40] ;
  char fnam[] = "readme_HOSTLIB" ;

  // -------- BEGIN ---------

  RAD = RADIAN ;
  NTMP = 0;

  cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
  sprintf(cptr, "Read %d total HOSTLIB entries from", HOSTLIB.NGAL_READ );
  cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
  sprintf(cptr, "  %s ", HOSTLIB.FILENAME );

  cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
  sprintf(ZNAME,"%s", "GENRANGE_REDSHIFT");
  sprintf(cptr, "Stored & z-sorted %d entries in %s= %5.3f - %5.3f"
          , HOSTLIB.NGAL_STORE, ZNAME
          , INPUTS.GENRANGE_REDSHIFT[0]
          , INPUTS.GENRANGE_REDSHIFT[1] );

  cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
  sprintf(cptr, "HOSTLIB Z-RANGE = %6.4f - %6.4f" , 
	  HOSTLIB.ZMIN, HOSTLIB.ZMAX );

  if ( INPUTS.HOSTLIB_GENRANGE_RA[1] < 361.0 ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "HOSTLIB RA-RANGE = %8.4f - %8.4f"
	    ,INPUTS.HOSTLIB_GENRANGE_RA[0]
	    ,INPUTS.HOSTLIB_GENRANGE_RA[1] );

    
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "HOSTLIB DEC-RANGE = %8.4f - %8.4f"
	    ,INPUTS.HOSTLIB_GENRANGE_DEC[0]
	    ,INPUTS.HOSTLIB_GENRANGE_DEC[1] );
  }


  cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
  sprintf(cptr, "HOSTLIB ZGAPMAX = %6.4f (at Z = %6.4f - %6.4f)",
          HOSTLIB.ZGAPMAX, HOSTLIB.Z_ATGAPMAX[0], HOSTLIB.Z_ATGAPMAX[1] );
  

  cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
  sprintf(cptr, "HOSTLIB <ZGAP>  = %7.5f", HOSTLIB.ZGAPAVG );
  

  // --------------------------------------
  // check HOSTLIB options and abort on inconsistencies.

  sprintf(copt,"HOSTLIB Opt: ");

  LGALMAG = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG ) ;
  if ( LGALMAG ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "%s compute host-noise contribution to SN noise", copt );    
  }

  LSN2GAL = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_Z ) ;
  if ( LSN2GAL ) {
    cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
    sprintf(cptr, "%s change SN redshift to host redshift", copt );    
  }

  LSN2GAL = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC ) ;
  if ( LSN2GAL ) {
    cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
    sprintf(cptr, "%s change SN position to host-SN pos.", copt );    

    if ( HOSTLIB.IVAR_RA < 0 || HOSTLIB.IVAR_DEC < 0 ) {
      sprintf(c1err ,"%s", "User requested SN pos -> host-SN pos");
      sprintf(c2err ,"%s", "but RA and/or DEC are missing from HOSTLIB.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
   
    fixran = INPUTS.HOSTLIB_FIXRAN_RADIUS ;
    if ( fixran > -0.0000001 && fixran < 1.000001 ) {  
      cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
      sprintf(cptr, "\t DEBUG OPT -> "
	      "Fix random number for reduced radius to %5.3f", fixran);       
    }

    fixran = INPUTS.HOSTLIB_FIXRAN_PHI ;
    if ( fixran > -0.0000001 && fixran < 1.000001 ) {
      cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
      sprintf(cptr, "\t DEBUG OPT -> "
	      "Fix random number for rotation angle to %5.3f", fixran );       
    }
  }

  // check for analytica ZPHOT model
  if ( INPUTS.USE_HOSTLIB_GENZPHOT ) {
      cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
      sprintf(cptr, "\t Use analytic ZPHOT[ERR] model " );      
  }

  // ---------------------
  // now store WGTMAP info

  cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
  sprintf(cptr,"Weight MAP variables: " ); 
  NVAR = HOSTLIB_WGTMAP.GRIDMAP.NDIM ; 
  NROW = HOSTLIB_WGTMAP.GRIDMAP.NROW ; 
  if ( NVAR == 0 ) 
    { strcat(cptr,"none"); }
  else {
    for(ivar=0; ivar < NVAR; ivar++ ) {
      strcat(cptr,HOSTLIB_WGTMAP.VARNAME[ivar] );
    }
  }

  if ( NVAR > 0 ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr,"Weight MAP size : %d ", NROW );
    cptr = HOSTLIB.COMMENT[NTMP];    NTMP++ ; 
    igal  = 0 ;
    GALID = get_GALID_HOSTLIB(igal) ;
    WGT   = HOSTLIB_WGTMAP.WGTSUM[igal];
    sprintf(cptr,"Weight of GALID=%lld : %f (1st HOSTLIB entry) ", GALID,WGT);
  }


  // summarize Sersic profiles
  for ( j=0; j < SERSIC_PROFILE.NDEF; j++ ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr,"GALSHAPE profile : %s %s"
	    ,SERSIC_PROFILE.VARNAME_a[j]
	    ,SERSIC_PROFILE.VARNAME_b[j] );
  }

  cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
  sprintf(cptr,"GALPOS generated with %5.4f of total flux.",
	  INPUTS.HOSTLIB_MXINTFLUX_SNPOS);


  // print out list of aperture radii
  LGALMAG = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG ) ;
  if ( LGALMAG ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "GALMAG %s interp-grid for PSFSIG(asec) = ",
	    HOSTLIB.filterList );
    for ( j=1; j <= NMAGPSF_HOSTLIB; j++ ) { // ??
      sprintf(ctmp, "%4.3f ", HOSTLIB.Aperture_PSFSIG[j] );
      strcat(cptr,ctmp);
      //  sprintf(cptr, "%s%4.3f ", cptr, HOSTLIB.Aperture_PSFSIG[j] );
    }


    cptr = HOSTLIB.COMMENT[NTMP] ; NTMP++ ; 
    sprintf(cptr,"%s", "GALMAG aperture integration bins: " );

    cptr = HOSTLIB.COMMENT[NTMP] ; NTMP++ ; 
    sprintf(cptr,"\t %3d x %6.3f arcsec radial bins.", 
	    NRBIN_GALMAG, HOSTLIB.Aperture_Rbin );

    cptr = HOSTLIB.COMMENT[NTMP] ; NTMP++ ; 
    sprintf(cptr,"\t %3d x %6.3f degree azimuthal bins.", 
	    NTHBIN_GALMAG, HOSTLIB.Aperture_THbin/RAD );

  } // end of LGALMAG if-block


  // store total number of comment/readme lines
  if ( NTMP >= MXCOMMENT_HOSTLIB ) {
    sprintf(c1err,"Number of HOSTLIB-comment lines=%d", NTMP);
    sprintf(c2err,"exceeds bound of MXCOMMENT_HOSTLIB = %d", 
	    MXCOMMENT_HOSTLIB);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  HOSTLIB.NLINE_COMMENT = NTMP;

  return ;

} // end of readme_HOSTLIB


// ========================================
double get_ZTRUE_HOSTLIB(int igal) {
  // Returns ZTRUE from redshift-sorted galaxy list
  double ZTRUE;
  int   ivar ;
  char fnam[] = "get_ZTRUE_HOSTLIB";
  // -------------- BEGIN ---------------
  ZTRUE = -9.0 ;
  if ( HOSTLIB.SORTFLAG == 0 ) {
    sprintf(c1err,"Cannot return sorted ZTRUE(igal=%d)", igal);
    sprintf(c2err,"until HOSTLIB is redshift-sorted.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  ivar  = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZTRUE,1);
  ZTRUE = HOSTLIB.VALUE_ZSORTED[ivar][igal] ;
  return ZTRUE ;
}

// ========================================
long long get_GALID_HOSTLIB(int igal) {
  // Returns GALID from redshift-sorted galaxy list
  double dval ;
  long long GALID ;
  int  ivar ;
  char fnam[] = "get_GALID_HOSTLIB";
  // -------------- BEGIN --------------
  if ( HOSTLIB.SORTFLAG == 0 ) {
    sprintf(c1err,"Cannot return sorted GALID(igal=%d)", igal);
    sprintf(c2err,"until HOSTLIB is redshift-sorted.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  ivar  = IVAR_HOSTLIB(HOSTLIB_VARNAME_GALID,1) ;
  dval  = HOSTLIB.VALUE_ZSORTED[ivar][igal] ;
  GALID = (long long)dval ;

  return GALID ;
}


// ========================================
int IVAR_HOSTLIB(char *varname, int ABORTFLAG) {


  // Mar 14, 2011
  // For input variable name return 'IVAR' index
  // to be used  the array HOSTLIB.VALUE_ZSORTED[ivar].
  // 
  // ABORTFLAG = 1  : abort if key not found
  // ABORTFLAG = 0  : return -9 if key nor found
  //
  // Feb 2014: use case-insenstive string check.

  int ivar, NVAR, ICMP ;
  char fnam[] = "IVAR_HOSTLIB" ;

  // ---------- BEGIN ----------

  NVAR = HOSTLIB.NVAR_STORE ; 
  for ( ivar = 0; ivar < NVAR; ivar++ ) {

    //    ICMP = strcmp(varname,HOSTLIB.VARNAME_STORE[ivar]) ;
    ICMP = strcmp_ignoreCase(varname,HOSTLIB.VARNAME_STORE[ivar]) ;
    
    if ( ICMP == 0 ) { return(ivar); }
  }

  // if we get here, abort
  if ( ABORTFLAG ) {
    sprintf(c1err,"Could not find IVAR index for '%s'", varname);
    sprintf(c2err,"Check VARNAMES keys in HOSTLIB file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    return(-9) ; 
  }
  else  { 
    return(-9) ; 
  }

} // end of IVAR_HOSTLIB


// ========================================
int ICOL_SPECBASIS(char *varname, int ABORTFLAG) {


  // June 2019
  // For input variable name return 'ICOL' column in specTemplate file.
  // 
  // ABORTFLAG = 1  : abort if key not found
  // ABORTFLAG = 0  : return -9 if key nor found
  //

  int icol, NCOL, ICMP ;
  char VARNAME_TMP[2][60];
  char fnam[] = "ICOL_SPECBASIS";

  // ---------- BEGIN ----------

  NCOL = HOSTSPEC.NSPECBASIS ;
  for ( icol = 0; icol < NCOL; icol++ ) {

    sprintf(VARNAME_TMP[0],"%s", 
	    HOSTSPEC.VARNAME_SPECBASIS[icol]);
    sprintf(VARNAME_TMP[1],"%s%s", 
	    PREFIX_SPECBASIS_HOSTLIB, HOSTSPEC.VARNAME_SPECBASIS[icol]);

    // check varname
    ICMP = strcmp_ignoreCase(varname,VARNAME_TMP[0] ) ;
    if ( ICMP == 0 ) { return(icol); }

    // check coeff_[varname]
    ICMP = strcmp_ignoreCase(varname,VARNAME_TMP[1] ) ;
    if ( ICMP == 0 ) { return(icol); }

  }

  // if we get here, abort
  if ( ABORTFLAG ) {
    sprintf(c1err,"Could not find specTemplate column '%s'", varname);
    sprintf(c2err,"Check VARNAMES keys in specTemplate file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    return(-9) ; 
  }
  else  { 
    return(-9) ; 
  }

} // end of ICOL_SPECBASIS

// =========================================
void GEN_SNHOST_DRIVER(double ZGEN_HELIO, double PEAKMJD) {

  // Mar 2011
  // Driver function to select host-galaxy from library,
  // and to determine other properties: 
  // photoz, SN pos, flux-noise ...
  // Note that ZGEN_HELIO is the desired host HELIOCENTRIC redshift,
  // which could be different than SN redshift if GENLC.CORRECT_HOSTMATCH=0.
  //
  // Check new FIXRAN_RADIUS and FIXRAN_PHI options
  // July 2015: add PEAKMJD arg to allow re-using host after MINDAYSEP

  int    USE, IGAL, ilist ;
  double fixran ;
  char   fnam[] = "GEN_SNHOST_DRIVER" ;

  // ------------ BEGIN -----------

  // always burn random numbers to stay synced.
  ilist = 1 ; 
  SNHOSTGAL.FlatRan1_GALID     = FlatRan1(ilist) ; // random GAL in smal z-bin
  SNHOSTGAL.FlatRan1_radius[0] = FlatRan1(ilist) ; // ran Sersic profile
  SNHOSTGAL.FlatRan1_radius[1] = FlatRan1(ilist) ; // random integral
  SNHOSTGAL.FlatRan1_phi       = FlatRan1(ilist) ;  // relative to major axis

  // ------------------------------------------------
  // check option to fix randoms (Sep 14, 2012)
  fixran = INPUTS.HOSTLIB_FIXRAN_RADIUS ;
  if ( fixran > -1.0E-9 ) { SNHOSTGAL.FlatRan1_radius[1] = fixran ; }

  fixran = INPUTS.HOSTLIB_FIXRAN_PHI   ;
  if ( fixran > -1.0E-9 ) { SNHOSTGAL.FlatRan1_phi = fixran ; }

  // ------------------------------------------------
  // init SNHOSTGAL values

  init_SNHOSTGAL();

  SNHOSTGAL.ZGEN              = ZGEN_HELIO ; 
  SNHOSTGAL.ZSPEC             = ZGEN_HELIO ;
  SNHOSTGAL.PEAKMJD           = PEAKMJD ;

  // check option to use HOST library
  USE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USE ) ;
  if ( USE == 0 ) { return ; }

  if ( ZGEN_HELIO < HOSTLIB.ZMIN || ZGEN_HELIO > HOSTLIB.ZMAX ) {
    double zCMB = zhelio_zcmb_translator(ZGEN_HELIO,GENLC.RA,GENLC.DEC,"eq",+1);
    sprintf(c1err,"Invalid ZGEN(helio,CMB)=%f ", ZGEN_HELIO, zCMB );
    sprintf(c2err,"HOSTLIB z-range is %6.4f to %6.4f",
	    HOSTLIB.ZMIN, HOSTLIB.ZMAX );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // select host-galaxy from library
  GEN_SNHOST_GALID(ZGEN_HELIO);
  IGAL = SNHOSTGAL.IGAL ;
  if ( IGAL < 0 ) { return ; } // Aug 2015

  // check on host photoz
  GEN_SNHOST_ZPHOT(IGAL);

  // generate SN position at galaxy
  GEN_SNHOST_POS(IGAL);

  // check of redshift needs to be updated (Apr 8 2019)
  TRANSFER_SNHOST_REDSHIFT(IGAL);

  // host-mag withing SN aperture
  GEN_SNHOST_GALMAG(IGAL);

  // load user-specified variables for output file
  LOAD_OUTVAR_HOSTLIB(IGAL); 

  // ---------------------------
  STORE_SNHOST_MISC(IGAL);


  // check debug-dump options
  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_DEBUG ) 
    {  DEBUG_1LINEDUMP_SNHOST(); }

  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_DUMP ) 
    { DUMP_SNHOST(); }


} // end of GEN_SNHOST_DRIVER



// =========================================
void GEN_SNHOST_GALID(double ZGEN) {
  
  // Mar 2011
  // Select weighted hostlib entry for redshift \simeq ZGEN(SN)
  // Fills
  // * SNHOSTGAL.IGAL
  // * SNHOSTGAL.GALID
  // * SNHOSTGAL.ZTRUE (GAL)
  //
  // Jan 2012: check new USEONCE flag
  //
  // Sep 2012: switch to logz binning; see iz_cen
  //
  // Oct 9 2013: ranval_GALID -> ranval_GALID * 0.95  to leave
  //             extra room to find an unused host.
  //
  // Nov 20 2015: compute ZCMB from ZTRUE = ZHELIO
  // Nov 23 2015: new algorithm to select igal_start & igal_end
  //              based on input key HOSTLIB_DZTOL
  // 
  // Dec 18 2015: if we have intentional wrong host, do NOT move SN redshift
  //              to match that of HOST. See GENLC.CORRECT_HOSTMATCH .
  //

  int 
    IZ_CEN, iz_cen, IGAL_SELECT
    ,igal_start, igal_end, igal, LDMP
    ,NSKIP_WGT, NSKIP_USED, NGAL_CHECK, MATCH; 
    ;

  long long GALID_FORCE, GALID ;
  double  ZTRUE, LOGZGEN ,WGT_start, WGT_end, WGT_dif, WGT_select, WGT ;
  double  ztol, dztol, z, z_start, z_end ;

  char fnam[] = "GEN_SNHOST_GALID" ;

  // ---------- BEGIN ------------


  IGAL_SELECT = -9 ; 

  // compute zSN-zGAL tolerance for this ZGEN = zSN
  dztol = INPUTS.HOSTLIB_DZTOL[0]
    +     INPUTS.HOSTLIB_DZTOL[1]*(ZGEN)
    +     INPUTS.HOSTLIB_DZTOL[2]*(ZGEN*ZGEN) ;


  // find start zbin 
  LOGZGEN = log10(ZGEN);

  if ( LOGZGEN < MINLOGZ_HOSTLIB ) {
    sprintf(c1err,"LOGZGEN=%.3f < MINLOGZ_HOSTLIB=%.3f",
	    LOGZGEN, MINLOGZ_HOSTLIB);
    sprintf(c2err,"Increase GENRANGE_REDSHIFT[0] or "
	    "descrease MINLOGZ_HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  IZ_CEN  = (int)((LOGZGEN-MINLOGZ_HOSTLIB)/DZPTR_HOSTLIB) ; 

  // Nov 23 2015: New algorithm to select min/max GALID using
  //              dztol from user

  iz_cen=IZ_CEN;    z=ZGEN ;  ztol=z-dztol ; 
  if ( iz_cen >= HOSTLIB.MAXiz ) 
    { igal_start = HOSTLIB.NGAL_STORE-1; }
  else 
    { igal_start = HOSTLIB.IZPTR[iz_cen+1]; }

  
  LDMP = ( GENLC.CID == -117 ) ; 

  while ( z > ztol && igal_start > 1 )   { 
    igal_start-- ; z = get_ZTRUE_HOSTLIB(igal_start);  
  }


  // - - - - - - - - - 

  iz_cen   = IZ_CEN;  z = ZGEN ; ztol = z+dztol ;
  if ( iz_cen < HOSTLIB.MINiz ) 
    { igal_end = 0; }
  else   
    { igal_end = HOSTLIB.IZPTR[iz_cen-1]; }
  
  while ( z < ztol && igal_end < HOSTLIB.NGAL_STORE-1 ) 
    { igal_end++ ; z = get_ZTRUE_HOSTLIB(igal_end);  }

  // back up one step to stay inside dztol
  igal_start++ ;    igal_end--   ; 
  z_start = get_ZTRUE_HOSTLIB(igal_start); 
  z_end   = get_ZTRUE_HOSTLIB(igal_end); 

  /*
  printf(" xxx ---------- %s DUMP ---------- \n", fnam);
  printf(" xxx ZGEN = %f  dztol=%f \n", ZGEN, dztol );
  printf(" xxx NEW igal(start,cen,end) = %d, %d, %d \n",
	 igal_start, igal_cen, igal_end);
  */


  if ( igal_end < igal_start ) { igal_end = igal_start ; }

  // error checks
  sprintf(c2err,"ZGEN=%f  dztol=%f  NGAL_STORE=%d", 
	  ZGEN, dztol, HOSTLIB.NGAL_STORE);
  if ( igal_start <= 0 || igal_start > HOSTLIB.NGAL_STORE ) {
    sprintf(c1err,"Invalid igal_start=%d ", igal_start );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  if ( igal_start < 0 || igal_start >= HOSTLIB.NGAL_STORE ) {
    sprintf(c1err,"Invalid igal_start=%d ", igal_start );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  if ( igal_end < 0 || igal_end >= HOSTLIB.NGAL_STORE ) {
    sprintf(c1err,"Invalid igal_end=%d ", igal_end );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  SNHOSTGAL.IGAL_SELECT_RANGE[0] = igal_start ;
  SNHOSTGAL.IGAL_SELECT_RANGE[1] = igal_end ;


  // now pick random igal between igal_start and igal_end
  WGT_start  = HOSTLIB_WGTMAP.WGTSUM[igal_start];
  WGT_end    = HOSTLIB_WGTMAP.WGTSUM[igal_end];
  WGT_dif    = WGT_end - WGT_start ;
  WGT_select = WGT_start + WGT_dif*SNHOSTGAL.FlatRan1_GALID * 0.95 ;

  NSKIP_WGT   = NSKIP_USED = NGAL_CHECK = 0 ;
  GALID_FORCE = INPUTS.HOSTLIB_GALID_FORCE ;

  // ---------------------------------------------------------
  // Feb 16 2016
  // Check option to pick out specific GALID_PRIORITY range.
  // Here just pick the first available GALID in range.
  // If there are no more GALID in range, then the nominal loop
  // below will pick from the full GALID range.
  long long GALID_MIN = INPUTS.HOSTLIB_GALID_PRIORITY[0] ;
  long long GALID_MAX = INPUTS.HOSTLIB_GALID_PRIORITY[1] ;
  if ( GALID_MIN < GALID_MAX  ) {
    for ( igal = igal_start; igal <= igal_end; igal++ ) {
      GALID = get_GALID_HOSTLIB(igal) ;
      if ( GALID < GALID_MIN ) { continue ; }
      if ( GALID > GALID_MAX ) { continue ; }
      if ( USEHOST_GALID(igal) == 0 ) { continue ; }
      IGAL_SELECT = igal; 
    }
    if ( IGAL_SELECT >= 0 ) { goto DONE_SELECT_GALID ; }
  } 


  // ---------------------------------------------------
  // start nominal loop over galaxies
  for ( igal = igal_start; igal <= igal_end; igal++ ) {

    NGAL_CHECK++ ;

    // check for forced GALID (Mar 12 2012)
    if ( GALID_FORCE > 0 ) {
      MATCH = (GALID_FORCE == get_GALID_HOSTLIB(igal) ) ; 
      if ( MATCH ) 
	{ IGAL_SELECT = igal; goto DONE_SELECT_GALID ; }
      else
	{ continue ; }
    }

    WGT  = HOSTLIB_WGTMAP.WGTSUM[igal]; 

    if ( WGT <  WGT_select  ) 
      { NSKIP_WGT++; continue ; }

    if ( USEHOST_GALID(igal) == 0 ) 
      { NSKIP_USED++ ; continue ; }

    IGAL_SELECT = igal ;
    goto DONE_SELECT_GALID ;

  }


  // - - - - - - - - - - - - - - - - - - - - - - - - 

 DONE_SELECT_GALID:

  if ( IGAL_SELECT < 0 ) {
    
    // if using same Galaxy with MJD-sep, just return so that
    // this event is rejected.
    if ( INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL < 9999 ) 
      { SNHOSTGAL.IGAL = IGAL_SELECT ; return ; }

    printf("\n PRE-ABORT DUMP: \n") ;
    printf("\t NSKIP_WGT  = %d/%d \n", NSKIP_WGT,  NGAL_CHECK );
    printf("\t NSKIP_USED = %d/%d \n", NSKIP_USED, NGAL_CHECK );
    printf("\t igal(start-end) = %d - %d\n", igal_start, igal_end);
    printf("\t WGT(start-end)  = %f - %f\n",  WGT_start, WGT_end );
    printf("\t WGT_select      = %f \n", WGT_select);

    DUMP_SNHOST();
    fflush(stdout);

    sprintf(c1err,"Could not find HOSTLIB entry for ZGEN=%5.4f .", ZGEN);
    sprintf(c2err,"           ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  // -----------------------------------------

  GALID  = get_GALID_HOSTLIB(IGAL_SELECT);
  ZTRUE  = get_ZTRUE_HOSTLIB(IGAL_SELECT);  // helio z

  GENLC.REDSHIFT_HOST  = ZTRUE ; // Jan 2016
  SNHOSTGAL.IGAL       = IGAL_SELECT ;
  SNHOSTGAL.GALID      = GALID ;
  SNHOSTGAL.ZTRUE      = ZTRUE  ;  
  SNHOSTGAL.ZDIF       = ZGEN - ZTRUE ; // zSN - zGAL (helio)

  if ( fabs(SNHOSTGAL.ZDIF) > dztol ) {
    printf("\n\n ---------- %s PRE-ABORT DUMP ---------- \n", fnam);
    printf("\t CID=%d  ZGEN(SN,GAL) = %.4f, %.4f \n",
	   GENLC.CID, ZGEN, ZTRUE );
    printf("\t igal(start,end) = %6d, %6d   --> IGAL=%d \n",
	   igal_start, igal_end, IGAL_SELECT);
    printf("\t zgal(start,end) = %.4f, %.4f \n" ,
	   z_start, z_end );
    printf("\t MINiz=%d  MAXiz=%d \n", HOSTLIB.MINiz, HOSTLIB.MAXiz);
    sprintf(c1err,"ZDIF=%f exceeds dztol=%f", SNHOSTGAL.ZDIF, dztol);
    sprintf(c2err,"CID=%d", GENLC.CID );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return ;

} // end of GEN_SNHOST_GALID



// ===============================
void init_SNHOSTGAL(void) {

  // init elements of SNHOSTGAL struct; called for each event.

  SNHOSTGAL.ZGEN              = -9.0 ;
  SNHOSTGAL.ZSPEC             = -9.0 ;
  SNHOSTGAL.PEAKMJD           = -9.0 ;
  SNHOSTGAL.WGTMAP_WGT        = 0.0 ;
  SNHOSTGAL.WGTMAP_SNMAGSHIFT = 0.0 ;  // default is no SN mag shift

  SNHOSTGAL.IGAL        = -9 ;
  SNHOSTGAL.GALID       =  INPUTS.HOSTLIB_GALID_NULL ; // Jul 15 2013
  SNHOSTGAL.ZTRUE       = -9.0 ;
  SNHOSTGAL.ZPHOT       = -9.0 ;
  SNHOSTGAL.ZPHOT_ERR   = -9.0 ;
  SNHOSTGAL.LOGMASS     = -9.0 ;
  SNHOSTGAL.LOGMASS_ERR = -9.0 ;  

  SNHOSTGAL.a_SNGALSEP_ASEC   = -999.0 ;
  SNHOSTGAL.b_SNGALSEP_ASEC   = -999.0 ;
  SNHOSTGAL.RA_SNGALSEP_ASEC  = -999.0 ;
  SNHOSTGAL.DEC_SNGALSEP_ASEC = -999.0 ;
  SNHOSTGAL.RA_GAL_DEG        = -999.0 ;
  SNHOSTGAL.DEC_GAL_DEG       = -999.0 ;
  SNHOSTGAL.RA_SN_DEG         = -999.0 ;
  SNHOSTGAL.DEC_SN_DEG        = -999.0 ;

  // always init GALMAG quantities to garbage
  int i, ifilt ;
  for ( i=0; i <= NMAGPSF_HOSTLIB ; i++ ) {
    SNHOSTGAL.GALFRAC[i]          = -9.0 ;     // global
    for ( ifilt = 0; ifilt < MXFILTINDX; ifilt++ ) { 
      SNHOSTGAL.GALMAG[ifilt][i]  = MAG_UNDEFINED ;
      SNHOSTGAL.SB_FLUX[ifilt]    = 0.0 ; 
      SNHOSTGAL.SB_MAG[ifilt]     = MAG_UNDEFINED ;
      SNHOSTGAL.GALMAG_TOT[ifilt] = MAG_UNDEFINED ;
    }
  }

} // end init_SNHOSTGAL


// ===============================
int USEHOST_GALID(int IGAL) {

  // July 2015
  // Returns 1 if this IGAL should be used.
  // Returns 0 if this IGAL should be rejected.
  //
  // Based on user-inputs determining if a host
  // can be re-used.

  int retCode,  NUSE_PRIOR, NUSE, use;
  int MJD_STORE[MXUSE_SAMEGAL], BLOCK_STORE[MXUSE_SAMEGAL], BLOCKSUM;
  int DAY, DAYDIF, ABSDIF;
  int DAYDIF_STORE[MXUSE_SAMEGAL], MINDAYSEP ;
  int LDMP, MJD;
  char fnam[] = "USEHOST_GALID" ;

  // --------------- BEGIN -----------------
  
  retCode = 0 ; // default is reject
  NUSE_PRIOR = SAMEHOST.NUSE[IGAL] ; 
  BLOCKSUM = 0 ;

  if ( SAMEHOST.REUSE_FLAG == 0 ) {
    // option to use each host just once
    if ( NUSE_PRIOR == 0 ) {
      SAMEHOST.NUSE[IGAL]++ ;
      retCode = 1 ;
    }
  }
  else if ( SAMEHOST.REUSE_FLAG == 1 ) {
    // option for arbitrary re-use of host regardless of PEAKMJD

    // don't increment past 1 to avoid array-bound overflow.
    if ( SAMEHOST.NUSE[IGAL] == 0 ) 
      { SAMEHOST.NUSE[IGAL]++ ; } 
    retCode = 1 ;
  }
  else if ( SAMEHOST.REUSE_FLAG == 2 ) {
    // option for re-use of host with PEAKMD-sep constraint
    BLOCKSUM = 0 ;
    MINDAYSEP = INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL ;
    DAY = (int)SNHOSTGAL.PEAKMJD - SAMEHOST.PEAKMJD_STORE_OFFSET ;
    // allocate memory for this IGAL on first usage
    if ( NUSE_PRIOR == 0 ) {
      SAMEHOST.PEAKDAY_STORE[IGAL] = 
	(unsigned short*)malloc( MXUSE_SAMEGAL*sizeof(unsigned short) ) ;
    }

    // check if any previously generated PEAKMJD on this host 
    // is too close in time.
    for ( use=0; use < NUSE_PRIOR; use++ )  {
      DAYDIF = DAY - SAMEHOST.PEAKDAY_STORE[IGAL][use] ;
      if ( DAYDIF >= 0 ) { ABSDIF=DAYDIF; } else { ABSDIF = -DAYDIF; }
      MJD  = 
	SAMEHOST.PEAKDAY_STORE[IGAL][use] + 
	SAMEHOST.PEAKMJD_STORE_OFFSET ;
      BLOCK_STORE[use] = 0 ;
      if ( ABSDIF < MINDAYSEP ) { BLOCKSUM = 1; BLOCK_STORE[use]=1; }

      MJD_STORE[use]    = MJD ;
      DAYDIF_STORE[use] = DAYDIF ;
    }
    
    // if current PEAKMJD was not blocked, then use this host [again]
    if ( BLOCKSUM == 0 ) {
      SAMEHOST.NUSE[IGAL]++ ;
      NUSE = SAMEHOST.NUSE[IGAL] ;
      SAMEHOST.PEAKDAY_STORE[IGAL][NUSE-1] = DAY ;
      retCode = 1 ;
    }
  }
  else {
    sprintf(c1err,"Unknown REUSE_FLAG = %d", SAMEHOST.REUSE_FLAG);
    sprintf(c2err,"check function init_SAMEHOST()");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // --------------------------------------------
  // ----------------- DEBUG DUMP ---------------
  LDMP = ( SAMEHOST.NUSE[IGAL] >= 2 && SAMEHOST.REUSE_FLAG == 2 ) ;
  LDMP = 0 ;
  if ( LDMP ) {
    printf(" xxx ================ CID = %d ================== \n", 
	   GENLC.CID);
    printf(" xxx IGAL=%d  PEAKMJD = %f  BLOCKSUM=%d \n", 
	   IGAL, SNHOSTGAL.PEAKMJD, BLOCKSUM );
    for ( use=0; use < NUSE_PRIOR; use++ )  {
      printf(" xxx\t MJD(prior) = %d  DAYDIF=%5d  BLOCK=%d \n",
	     MJD_STORE[use], DAYDIF_STORE[use], BLOCK_STORE[use] );
    }
    printf(" xxx retCode = %d \n", retCode );
    fflush(stdout);
  } // end of LDMP
  // --------------------------------------------

  NUSE = SAMEHOST.NUSE[IGAL] ;
  if ( NUSE >= MXUSE_SAMEGAL ) {
    printf("\n PRE-ABORT DUMP: \n");
    printf("\t IGAL=%d   PEAKMJD = %f  BLOCKSUM=%d \n", 
	   IGAL, SNHOSTGAL.PEAKMJD, BLOCKSUM );
    for ( use=0; use < NUSE_PRIOR; use++ )  {
      printf("\t MJD(previous) = %d  DAYDIF=%5d  BLOCK=%d \n",
	     MJD_STORE[use], DAYDIF_STORE[use], BLOCK_STORE[use] );
    }

    sprintf(c1err,"NUSE=%d exceeds bound (MXUSE_SAMEGAL)", NUSE);
    sprintf(c2err,"Increase MXUSE_SAMEGAL or add more galaxies to HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return retCode;

}  // end of USEHOST_GALID

// =========================================
void FREEHOST_GALID(int IGAL) {

  // Sep 11 2015
  // If generated event is rejected, and host galaxies are re-used
  // with MJD-separation, this function frees  up the host so that 
  // it can be used again.

  if ( IGAL < 0 ) { return ; }

  if ( SAMEHOST.REUSE_FLAG == 2 ) {

    SAMEHOST.NUSE[IGAL]-- ;
  }

  return ; 

} // end UNUSE_HOST_GALID

// =========================================
void GEN_SNHOST_ZPHOT(int IGAL) {

  // If photoz is in librray, then fill
  // * SNHOSTGAL.ZPHOT
  // * SNHOSTGAL.ZPHOT_ERR
  //
  // and adjust ZPHOT to account for SN-galaxy difference
  // in redshift.
  //
  // Aug 18 2015: check  HOSTLIB_ZPHOTEFF option
  // Feb 23 2017: add option to compute ZPHOT[_ERR] from formulas
  //              See new functions GEN_SNHOST_ZPHOT_from_XXXX
  //
  // Oct 04 2018: fix awful bug with checking HOSTLIB_GENZPHOT usage
  //              
  int j;
  double ZPHOT, ZPHOT_ERR, ZBIAS, tmp ;
  double ZTRUE = SNHOSTGAL.ZTRUE ;
  char fnam[] = "GEN_SNHOST_ZPHOT" ;

  // ----------- BEGIN ----------

  if ( INPUTS.USE_HOSTLIB_GENZPHOT  ) 
    { GEN_SNHOST_ZPHOT_from_CALC(IGAL,&ZPHOT,&ZPHOT_ERR); }
  else
    { GEN_SNHOST_ZPHOT_from_HOSTLIB(IGAL,&ZPHOT,&ZPHOT_ERR); }

  // Apply user-specified bias
  ZBIAS = 0 ;
  for(j=0; j < 4; j++ ) {
    tmp = (double)INPUTS.HOSTLIB_GENZPHOT_BIAS[j];
    if ( tmp != 0.0 ) {
      ZBIAS += tmp * pow(ZTRUE, (double)j );
    }
  }
  ZPHOT += ZBIAS;
  SNHOSTGAL.ZPHOT       = ZPHOT ;
  SNHOSTGAL.ZPHOT_ERR   = ZPHOT_ERR ;
  

  // -------------------------------
  // Aug 18 2015
  // check ZPHOTEFF option for host-zphot efficiency --> 
  // set some ZPHOT values to -9 to reflect 'no-host' cases.
  // Use randome number for azimuthal angle to avoid burning
  // another random and changing the sync.
  double flatRan, EFF;
  int NZBIN = HOSTLIB_ZPHOTEFF.NZBIN ;
  if ( NZBIN > 0 ) {
    flatRan =  SNHOSTGAL.FlatRan1_phi ;
    EFF  = interp_1DFUN(1, SNHOSTGAL.ZTRUE, NZBIN, 
			HOSTLIB_ZPHOTEFF.ZTRUE, HOSTLIB_ZPHOTEFF.EFF, fnam );
    if ( flatRan > EFF )  { SNHOSTGAL.ZPHOT = SNHOSTGAL.ZPHOT_ERR = -9.0 ; }
  }

  // ---------------------------------------------------------
  // negative redshift error is a flag to use the host photo-Z
  if ( INPUTS.GENSIGMA_REDSHIFT < 0.0 ) {
    GENLC.REDSHIFT_CMB_SMEAR     = SNHOSTGAL.ZPHOT ;
    GENLC.REDSHIFT_HELIO_SMEAR   = SNHOSTGAL.ZPHOT ;
    GENLC.REDSHIFT_SMEAR_ERR     = SNHOSTGAL.ZPHOT_ERR ;
  }


  return ;

} // end of  GEN_SNHOST_ZPHOT


// =======================================
void GEN_SNHOST_ZPHOT_from_CALC(int IGAL, double *ZPHOT, double *ZPHOT_ERR) {

  // Created Feb 23 2017
  // Compute ZPHOT and ZPHOT_ERR from Gaussian profiles specified
  // by sim-input  HOSTLIB_GENZPHOT_FUDGEPAR
  //
  // Mar 29 2018: 
  //  + float -> double
  //  + ZPHOT_ERR = sigma_core, not true sigma
  //
  // May 10 2018: 
  //  + if sigz_outlier >= 10, pick from flat distribution over 
  //    entire HOSTLIB z-range
  //
  // June 7 2018: protect against ZPHOT < 0

  double  sigz1_core[3] ;  // sigma/(1+z): a0 + a1*(1+z) + a2*(1+z)^2
  double  sigz1_outlier, prob_outlier, zpeak, sigz_lo, sigz_hi ;
  double  sigma_core, sigma_outlier, HOSTLIB_ZRANGE[2] ;
  double  z, z1, ranGauss, ranProb, ranzFlat, zphotErr, zshift_ran ;

  GENGAUSS_ASYM_DEF ZPHOTERR_ASYMGAUSS ;

  int    OUTLIER_FLAT=0;
  double SIGMA_OUTLIER_FLAT = 9.999 ; // pick random z if sig_outlier> this
  //  char   fnam[] = "GEN_SNHOST_ZPHOT_from_CALC" ;

  // ----------- BEGIN -------------

  *ZPHOT = *ZPHOT_ERR = -9.0 ;
  
  sigz1_core[0] = (double)INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[0] ;
  sigz1_core[1] = (double)INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[1] ;
  sigz1_core[2] = (double)INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[2] ;

  prob_outlier  = (double)INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[3] ;
  sigz1_outlier = (double)INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[4] ;

  // - - - - - -

  z  = GENLC.REDSHIFT_HELIO ;
  z1 = 1.0 + z ;
  sigma_core = (sigz1_core[0] + 
		sigz1_core[1]*z1 + 
		sigz1_core[2]*(z1*z1) ) ;

  sigma_outlier = sigz1_outlier*z1 ;
  
 PICKRAN:


  ranProb  = FlatRan1(1) ;

  if ( ranProb > prob_outlier ) 
    { zphotErr = sigma_core ;  }
  else  { 
    if ( sigz1_outlier > SIGMA_OUTLIER_FLAT ) 
      {  OUTLIER_FLAT=1; }
    else 
      { zphotErr = sigma_outlier ; }
  }

  if ( OUTLIER_FLAT ) {
    HOSTLIB_ZRANGE[0] = INPUTS.HOSTLIB_GENZPHOT_OUTLIER[0] ;
    HOSTLIB_ZRANGE[1] = INPUTS.HOSTLIB_GENZPHOT_OUTLIER[1] ;
    ranzFlat = FlatRan(2, HOSTLIB_ZRANGE); 
    *ZPHOT = ranzFlat;    
  }
  else {
    // to avoid negative ZPHOT (or truncated ZPHOT),
    // check to modify symmetric error into asymmetric error

    zphoterr_asym(z, zphotErr, &ZPHOTERR_ASYMGAUSS);
    zpeak   = ZPHOTERR_ASYMGAUSS.PEAK;
    sigz_lo = ZPHOTERR_ASYMGAUSS.SIGMA[0] ;
    sigz_hi = ZPHOTERR_ASYMGAUSS.SIGMA[1] ;
    *ZPHOT  = zpeak + biGaussRan(sigz_lo, sigz_hi );

    /*
    printf(" xxx zpeak=%.4f sig(-/+)=%.4f/%.4f  ZPHOT=%.3f \n",
	   zpeak, sigz_lo, sigz_hi, *ZPHOT); fflush(stdout);
    */

  }

  if ( *ZPHOT < 0.001 ) { goto PICKRAN; }

  // reported ZPHOT_ERR is sigma_core, or the RMS of asym Gaussian.
  zphoterr_asym(z, sigma_core, &ZPHOTERR_ASYMGAUSS) ;
  *ZPHOT_ERR = ZPHOTERR_ASYMGAUSS.RMS ;
 

  return ;

} // end GEN_SNHOST_ZPHOT_from_CALC


// ===================================================
void zphoterr_asym(double ZTRUE, double ZPHOTERR,
		   GENGAUSS_ASYM_DEF *asymGaussPar ) {

  // Created Jun 8 2018 
  //
  // For inpts ZTRUE and ZPHOTERR,
  // returns  assymetric Gaussian parameters:
  //  asymGaussPar[0] = Zpeak (z with max prob)
  //  asymGaussPar[1] = sigma_plus
  //  asymGaussPar[2] = sigma_minus
  //
  //  MINSIGNA0=3 (below) is how many stdev to require between 
  //  ZMIN and Zpeak:   Zpeak >= ZMIN + 3*sigma_minus
  //
  //  If there is no solution for Zpeak, then reduce
  //  ZPHOTERR_LOCAL to make sure a solution exists.
  //
  // Useful relations for asymmetric Gaussian:
  //  z_mean - z_peak = (sigma_plus - sigma_minus) * sqrt(2/PI)
  //  VARIANCE    = SIGMA_SUM^2 + SIGMA_DIF^2(1-2/PI)
  //        -->  SIGMA_SUM = 0.5*(sigma_plus + sigma_minus)
  //        -->  SIGMA_DIF = 0.5*(sigma_plus - sigma_minus)
  //
  // In principle there is an analytic solution with a 
  // nasty quadratic equation, but here it is solved
  // numerically by looping over sigma_minus and solving
  // for the VARIANCE. Constraint is that 
  //   VARIANCE(asymGauss) = ZPHOTERR^2 = VARIANCE(input).
  //
  // Beware of a few hard-wired numbers in the code

#define PI             3.1415926535
#define SQRT2_over_PI  0.797885
#define SQRT8_over_PI  1.59577

  int    NITER=2, iter, NTRY_ZRATIO=0;
  double ZRATIO_SCALE[2] = { 1.10, 1.003 };
  double ZMIN             = 0.01 ;
  double ZRATIO_MIN_SOLVE = 0.75 ;
  double MINSIGMA0        = 3.0 ; // min Nsigma for (Zpeak-ZMIN)/sigma_minus
  double ZERO             = 0.0 ;

  double ZPHOTERR_LOCAL, SQZPHOTERR, RMS ;
  double sigma_minus, sigma_plus, zpeak, zpeak_atVMIN=0.0;
  double sigma_minus_atVMIN=0.0, sigma_plus_atVMIN=0.0;
  double SIGMA_SUM, SIGMA_DIF, VARIANCE, V1, V2, VRATIO, VDIF, VDIFMIN;
  double ZRATIO, ZRATIO_START, ZRATIO_END, ZRATIO_MIN, ZRATIO_MAX;
  double zratio_tmp, zratio_atVMIN=0.0, VARIANCE_atVMIN=0.0 ;

  char   text[20];
  int    LDMP = 0 ;
  char fnam[] = "zphoterr_asym" ;

  // --------------- BEGIN ---------------

  if ( LDMP ) { printf("\n xxxxx DUMP %s xxxxx \n", fnam ); }

  init_GENGAUSS_ASYM(asymGaussPar, ZERO );

  ZPHOTERR_LOCAL = ZPHOTERR ;
  ZRATIO = (ZTRUE-ZMIN)/ZPHOTERR_LOCAL;

  // if already more than 3 sigma away from ZMIN, 
  // keep symmetric Gauss params and bail.
  if ( ZRATIO > MINSIGMA0 ) {
    zpeak = ZTRUE ;
    sigma_plus = sigma_minus = RMS = ZPHOTERR ;
    goto LOAD_GAUSSPAR;
  }

  // if ZRATIO is too small there is no solution ;
  // in this case force ZPHOTERR to be smaller to force a solution.
  if ( ZRATIO < ZRATIO_MIN_SOLVE ) {
    double zratio_orig = ZRATIO ;
    ZPHOTERR_LOCAL = (ZTRUE-ZMIN)/ZRATIO_MIN_SOLVE ;
    ZRATIO         = ZRATIO_MIN_SOLVE ;

    if ( LDMP ) { 
      printf("# WARNING: ZRATIO=%.3f < %.3f ==> ZPHOTERR=%.3f -> %.3f \n",
	     zratio_orig, ZRATIO, ZPHOTERR, ZPHOTERR_LOCAL);
    }
  }


  SQZPHOTERR = ZPHOTERR_LOCAL*ZPHOTERR_LOCAL;

  if ( LDMP ) {
    printf("# ZTRUE = %.3f  ZPHOTERR=%.3f  (zRatio=%.3f)\n", 
	   ZTRUE, ZPHOTERR_LOCAL, ZRATIO );
    printf("#\n");
    printf("# NITER  zratio   zpeak   sigma(-/+)   VAR/ZPHOTERR^2  \n");
    printf("# ------------------------------------------------------\n");
  }


  ZRATIO_MAX = ZRATIO;
  ZRATIO_MIN = 0.02 ;

  for(iter=0; iter<NITER; iter++ ) {

    NTRY_ZRATIO = 0 ;
    if ( iter == 0 ) {
      ZRATIO_START = ZRATIO_MIN ;
      ZRATIO_END   = ZRATIO_MAX ;
    }
    else {
      ZRATIO_START = zratio_atVMIN * 0.96 ;
      ZRATIO_END   = zratio_atVMIN * 1.04 ;
      if ( ZRATIO_START < ZRATIO_MIN ) { ZRATIO_START = ZRATIO_MIN; }
      if ( ZRATIO_END   > ZRATIO_MAX ) { ZRATIO_END   = ZRATIO_MAX; }
    }
    VDIFMIN = 99999.0 ;
    zratio_tmp = ZRATIO_START / ZRATIO_SCALE[iter] ;

    if ( LDMP ) {
      printf("\t (START iter=%d: %.4f < zratio < %.4f) \n",
	     iter, ZRATIO_START, ZRATIO_END);
    }

    while ( zratio_tmp < ZRATIO_END ) {
      NTRY_ZRATIO++ ;
      zratio_tmp  *= ZRATIO_SCALE[iter] ;
      sigma_minus  = (ZTRUE-ZMIN) * zratio_tmp;
      SIGMA_DIF    = (ZTRUE - ZMIN - MINSIGMA0*sigma_minus)/SQRT8_over_PI ;
      sigma_plus   = 2.0*SIGMA_DIF + sigma_minus ;
      SIGMA_SUM    = 0.5*(sigma_plus + sigma_minus) ;
      if ( sigma_plus < 0.001 ) { continue ; }

      V1           = SIGMA_SUM * SIGMA_SUM ;
      V2           = SIGMA_DIF * SIGMA_DIF * (1.0 - 2.0/PI);
      VARIANCE     = V1 + V2 ;
      zpeak        = ZMIN + (MINSIGMA0 * sigma_minus);
      VRATIO       = VARIANCE/SQZPHOTERR;
      VDIF         = fabs(VRATIO-1.0); // minmize this quantity
      if ( VDIF < VDIFMIN ) { 
	VDIFMIN=VDIF;  
	zratio_atVMIN      = zratio_tmp ;
	sigma_minus_atVMIN = sigma_minus; 
	sigma_plus_atVMIN  = sigma_plus; 
	zpeak_atVMIN       = zpeak ;
	VARIANCE_atVMIN    = VARIANCE ;
      }

      if ( LDMP ) {
	text[0] = 0 ;
	if ( VDIF < 0.02 ) { sprintf(text, "*") ; }
	printf("  %d-%3d  %.4f   %.4f  %.4f/%.4f   %.4f  %s\n",
	       iter, NTRY_ZRATIO, zratio_tmp,
	       zpeak, sigma_minus, sigma_plus, VRATIO, text );
      }
    }
  } // end iter
  
  // - - - - - load values at minium VARIANCE/ZPHOTERR^2

  sigma_minus = sigma_minus_atVMIN ;
  sigma_plus  = sigma_plus_atVMIN  ;
  zpeak       =	zpeak_atVMIN ;
  RMS         = sqrt(VARIANCE_atVMIN);

 LOAD_GAUSSPAR:
  asymGaussPar->PEAK     = zpeak;
  asymGaussPar->SIGMA[0] = sigma_minus ;
  asymGaussPar->SIGMA[1] = sigma_plus  ;
  asymGaussPar->RMS      = RMS ;
  return ;

} // end zphoterr_asym


// =================================================
void GEN_SNHOST_ZPHOT_from_HOSTLIB(int IGAL, double *ZPHOT, double *ZPHOT_ERR) {

  // Created Feb 23 2017:
  // Generate host-ZPHOT from host library.
  //
  // Note that reported ZPHOT += [ zGEN(HELIO) - ZTRUE(HOSTLIB) ]
  // However, ZPHOTERR is not adjusted, which assumes that the
  // library ZTRUE is always close to generated redshift.
  //
  // Mar 28 2018: float -> double
  // May 17 2018: if either zphot<0 or zphoterr<0, set both
  //              to -9. These indicate failed host photo-z fit.
  //
  // Apr 06 2019: initialize outputs (*ZPHOT,*ZPHOT_ERR) to -9.
  //

  int IVAR_ZPHOT, IVAR_ZPHOT_ERR ;
  double ZDIF, zphot_local, zerr_local ;
  //  char fnam[] = "GEN_SNHOST_ZPHOT_from_HOSTLIB" ;

  // ----------- BEGIN -----------

  *ZPHOT  = *ZPHOT_ERR = -9.0; 
  IVAR_ZPHOT     = HOSTLIB.IVAR_ZPHOT ;
  IVAR_ZPHOT_ERR = HOSTLIB.IVAR_ZPHOT_ERR ;

  if ( IVAR_ZPHOT < 0 ) { return ; }

  if ( GENLC.CORRECT_HOSTMATCH ) 
    { ZDIF = GENLC.REDSHIFT_HELIO - SNHOSTGAL.ZTRUE ; }
  else
    { ZDIF = 0.0 ; }    // wrong host

  zphot_local   = HOSTLIB.VALUE_ZSORTED[IVAR_ZPHOT][IGAL] ;
  zphot_local  += ZDIF ;

  zerr_local  = HOSTLIB.VALUE_ZSORTED[IVAR_ZPHOT_ERR][IGAL] ;

  if ( zphot_local < 0.0 || zerr_local < 0.0 ) 
    { zphot_local = zerr_local = -9.0; }
  

  // load output args
  *ZPHOT     = zphot_local;
  *ZPHOT_ERR = zerr_local;

} // end GEN_SNHOST_ZPHOT_from_HOSTLIB

// =======================================
void GEN_SNHOST_POS(int IGAL) {

  // Mar 10, 2011
  // Fill SNHOSTGAL.[position]
  // Input IGAL is the sequential library index
  //
  // Galaxy profile is a sum of Sersic profiles.
  // First randomly pick a profile based on its total weight (flux).
  // Then pick a random location from this profile using 
  // pre-tabulated integrals of flux vs. reduced radius (r=R/Rhalf).
  //
  // if ( LSN2GAL  ) then call gen_MWEBV with new SN coords.
  // 
  // Nov 12 2015: make a few cos(DEC) fixes ; bugs found by Ravi @ANL
  // Nov 16 2015: use exact angSep function to compute SN-host sep.
  //
  // Feb 04 2019: compute DLR, DDLR, and make it work even if
  //              host galaxy RA,DEC are not given.
  //
  
  int  LSN2GAL, LDEBUG, IVAR_RA, IVAR_DEC, IVAR_ANGLE;
  int j, JPROF, k_table, NBIN    ;

  double 
    RA_GAL, DEC_GAL
    ,phi, cphi, sphi, RAD, crot, srot
    ,reduced_logR0, reduced_logR1, reduced_logR, reduced_R
    ,Ran0, Ran1, WGT, RanInteg, dif, bin, fbin
    ,a, b, a_half, b_half, ang, n, inv_n, DTMP, DCOS, SNSEP, DLR    
    ,*ptr, *ptr_r, *ptr_integ0, *ptr_integ1
    ;

  char fnam[] = "GEN_SNHOST_POS" ;

  // -------------- BEGIN -------------

  SNHOSTGAL.phi               =  0.0 ;
  SNHOSTGAL.reduced_R         = -999. ;
  SNHOSTGAL.SERSIC_INDEX      = -999. ;
  SNHOSTGAL.RA_SNGALSEP_ASEC  = -999. ;
  SNHOSTGAL.DEC_SNGALSEP_ASEC = -999. ;
  SNHOSTGAL.SNSEP             = -999. ;
  SNHOSTGAL.DLR               = -999. ;
  SNHOSTGAL.DDLR              = -999. ;

  // bail out if there are no galaxy shape parameters
  if ( SERSIC_PROFILE.NDEF == 0 ) { return ; }

  // extract info for each Sersic term
  get_Sersic_info(IGAL) ;    

  // strip off user options passed via sim-input file
  LSN2GAL = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC ) ;
  LDEBUG  = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_DEBUG ) ;

  // strip off indices
  IVAR_RA     = HOSTLIB.IVAR_RA ;
  IVAR_DEC    = HOSTLIB.IVAR_DEC ;
  IVAR_ANGLE  = HOSTLIB.IVAR_ANGLE ;
  RAD         = RADIAN ;

  // strip off random numbers to randomly generate a host-location
  Ran0  = SNHOSTGAL.FlatRan1_radius[0] ;
  Ran1  = SNHOSTGAL.FlatRan1_radius[1] ;
  phi   = SNHOSTGAL.FlatRan1_phi * TWOPI ; //azimuth angle rel. to major axis

  // Pick reduced radius (r=R/Rhalf),
  // Start by randonly picking (Ran0) which Sersic profile 
  // based on the WGT of each profile.

  JPROF = -9;
  for ( j=0; j < SERSIC_PROFILE.NDEF; j++ ) {
    WGT = SNHOSTGAL.SERSIC_wsum[j];
    if ( WGT >= Ran0 && JPROF < 0 ) { JPROF = j ; }
  }

  // bail if we cannot pick a Sersic profile.
  if ( JPROF < 0 ) {
    ptr = SNHOSTGAL.SERSIC_wsum ; 
    sprintf(c1err,"Could not find random Sersic profile for Ran0=%f", Ran0);
    sprintf(c2err,"SERSIC_wsum = %f %f %f %f",
	    ptr[0], ptr[1], ptr[2], ptr[3] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // get integral-table index for this Sersic_n
  n        = SNHOSTGAL.SERSIC_n[JPROF]; 
  inv_n    = 1.0/n ; 
  dif      = (inv_n - SERSIC_TABLE.INVINDEX_MIN) ;
  bin      = SERSIC_TABLE.INVINDEX_BIN ;
  k_table  = (int)(dif/bin) ;

  if ( k_table < 0 || k_table >= NSERSIC_TABLE ) {
    sprintf(c1err,"Invalid SERSIC_TABLE ptr = %d", k_table );
    sprintf(c2err,"occurs for n=%5.3f  JPROF=%d", n, JPROF);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // pick reduced_R = R/Rhalf from this "JPROF" Sersic profile.
  // To limit very distant SN from the long Sersic tails,
  // 'MXINTFLUX_SNPOS' of the integrated flux is used to
  // determine the reduced radius.
  RanInteg    = Ran1 * INPUTS.HOSTLIB_MXINTFLUX_SNPOS ;
  ptr_r       = SERSIC_TABLE.reduced_logR ;  
  ptr_integ0  = SERSIC_TABLE.INTEG_CUM[k_table] ;
  ptr_integ1  = SERSIC_TABLE.INTEG_CUM[k_table+1] ;
  NBIN        = NBIN_RADIUS_SERSIC + 1 ;

  reduced_logR0  = interp_1DFUN(1, RanInteg, NBIN, ptr_integ0, ptr_r, 
				"reduced_logR0" );
  reduced_logR1  = interp_1DFUN(1, RanInteg, NBIN, ptr_integ1, ptr_r, 
				"reduced_logR1" );

  // interpolate integral tables in 1/n space
  dif  = inv_n -  1./SERSIC_TABLE.n[k_table] ;
  fbin = dif/bin ;		
  reduced_logR = reduced_logR0 + fbin*(reduced_logR1 - reduced_logR0);
  reduced_R    = pow(10.0,reduced_logR);

  /*
  printf(" ---------------------------------------------------- \n");
  printf(" kkkkkk = %d  1/n=%5.3f  dif=%6.3f  bin=%6.3f  fbin=%f \n", 
	 k_table, inv_n, dif, bin, fbin);
  printf(" xxx R0=%5.3f R1=%5.3f R(%d)=%5.3f   PHI=%5.1f deg. RanInteg=%f \n",
	 reduced_R0, reduced_R1, JPROF, reduced_R, phi/RAD, RanInteg );
  */


  // check test-option from sim-input file to fix the gal size
  if ( LDEBUG ) { 
    //    SNHOSTGAL.SERSIC_a[JPROF] = 2.0 ;
    //    SNHOSTGAL.SERSIC_b[JPROF] = 1.0 ;
    HOSTLIB.VALUE_ZSORTED[IVAR_ANGLE][IGAL] = 0.0 ; // a_rot angle
    // phi = 0.0 ; // along major axis only
    //  reduced_R = 0.0 ; // 
  } 


  // get major and minor half-light axes (arcsec) for this 
  // Sersic profile and this galaxy

  if ( HOSTLIB.IVAR_ANGLE < 0 ) {
    sprintf(c1err,"Missing required %s in hostlib", 
	    HOSTLIB_VARNAME_ANGLE);
    sprintf(c2err,"Needed to choose position near host.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  a_half   = SNHOSTGAL.SERSIC_a[JPROF]; // half-light radius, major axis
  b_half   = SNHOSTGAL.SERSIC_b[JPROF]; // half-light radius, minor axis
  ang      = HOSTLIB.VALUE_ZSORTED[IVAR_ANGLE][IGAL] ; // a_rot angle
  cphi = cos(phi);  sphi = sin(phi);

  // Feb 4 2019: 
  // get DLR, distance to half-light ellipse, in direction of SN
  // See Eq 3.7 in Gupta hostmatching thesis
  double top, bottom;
  top    = a_half * b_half;
  bottom = sqrt(a_half*a_half*sphi*sphi + b_half*b_half*cphi*cphi ) ; 
  DLR    = top/bottom ;

  // get coords (arcsec) in ellipse-frame
  a = reduced_R * a_half * cphi ;
  b = reduced_R * b_half * sphi ;

  SNHOSTGAL.a_SNGALSEP_ASEC  = a ;
  SNHOSTGAL.b_SNGALSEP_ASEC  = b ;
  SNHOSTGAL.phi              = phi ;
  SNHOSTGAL.reduced_R        = reduced_R ;
  SNHOSTGAL.SERSIC_INDEX     = n ;
  SNHOSTGAL.DLR              = DLR ;

  // now rotate coords based on major axis rotation angle w.r.t. RA

  crot = cos(ang*RAD);
  srot = sin(ang*RAD);

  SNHOSTGAL.RA_SNGALSEP_ASEC  = a * crot + b * srot ;
  SNHOSTGAL.DEC_SNGALSEP_ASEC = b * crot - a * srot ;
  
  // ----------------------------------------------------
  // Determine Galaxy and SN coords

  if ( IVAR_RA >= 0 ) {
    // get absolute SN coords if we have host sky coords.
    // fetch galaxy coords from library
    RA_GAL   = HOSTLIB.VALUE_ZSORTED[IVAR_RA][IGAL] ;  // degrees
    DEC_GAL  = HOSTLIB.VALUE_ZSORTED[IVAR_DEC][IGAL] ; // degrees
    
    // compute absolute SN position relative to center of host
    // (Nov 2015 - correct for 1/cos(DEC))
    DCOS = cos(DEC_GAL*RAD) ;
    DTMP                 = DEG_ARCSEC * SNHOSTGAL.RA_SNGALSEP_ASEC / DCOS ;
    SNHOSTGAL.RA_SN_DEG  = RA_GAL  + DTMP ;
    
    DTMP                 = DEG_ARCSEC * SNHOSTGAL.DEC_SNGALSEP_ASEC ;
    SNHOSTGAL.DEC_SN_DEG = DEC_GAL + DTMP ;
  }
  else {
    // Feb 2019: no host coords, so fudge host position relative 
    // to already-determined SN pos. Goal is to get SNSEP and DLR.
    //
    SNHOSTGAL.RA_SN_DEG   = GENLC.RA ; // SN coord already selected
    SNHOSTGAL.DEC_SN_DEG  = GENLC.DEC ;

    DCOS    = cos(SNHOSTGAL.DEC_SN_DEG*RAD) ; 
    DTMP    = DEG_ARCSEC * SNHOSTGAL.RA_SNGALSEP_ASEC / DCOS ;
    RA_GAL  = SNHOSTGAL.RA_SN_DEG  - DTMP ;

    DTMP    = DEG_ARCSEC * SNHOSTGAL.DEC_SNGALSEP_ASEC ;
    DEC_GAL = SNHOSTGAL.DEC_SN_DEG - DTMP ;
  }

  // load gal coords into global
  SNHOSTGAL.RA_GAL_DEG  = RA_GAL ;
  SNHOSTGAL.DEC_GAL_DEG = DEC_GAL ;

  // compute SN-host separation in arcsec.
  SNSEP = angSep(SNHOSTGAL.RA_GAL_DEG, SNHOSTGAL.DEC_GAL_DEG,
		 SNHOSTGAL.RA_SN_DEG, SNHOSTGAL.DEC_SN_DEG,
		 (double)3600.);

  SNHOSTGAL.SNSEP = SNSEP ;
  SNHOSTGAL.DDLR  = SNSEP / DLR ; // Feb 2019

  // check user sim-input option to change previously 
  // selected SN coord to SNHOST coord. 
  // Likely use is for input to Image sim where the
  // SN coords must correspond to the galaxy location.

  if ( LSN2GAL ) {
    GENLC.RA   = SNHOSTGAL.RA_SN_DEG ;
    GENLC.DEC  = SNHOSTGAL.DEC_SN_DEG ;
    gen_MWEBV(); // compute MWEBV with SN coords (was skipped in snlc_sim)
  }


} // end of GEN_SNHOST_POS


// =================================
void TRANSFER_SNHOST_REDSHIFT(IGAL) {

  // Apr 2019
  // Check for redshift change; either explicit transfer to use
  // zHOST, or switch to host coords and change zHEL.
  //
  // Apr 20 2019: bail on incorrect host match.
  //

  double ZTRUE        = SNHOSTGAL.ZTRUE ;             // helio redshift
  double zPEC         = GENLC.VPEC/LIGHT_km ;
  double ZTRUE_noVPEC = ZTRUE - zPEC ;
  double RA           = GENLC.RA ;
  double DEC          = GENLC.DEC ;
  int    MSKOPT       = INPUTS.HOSTLIB_MSKOPT ;
  int OVP_Z           = (MSKOPT & HOSTLIB_MSKOPT_SN2GAL_Z) ;
  int OVP_RADEC       = (MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC) ;
  char eq[]           = "eq";
  //  char fnam[]         = "TRANSFER_SNHOST_REDSHIFT" ;

  double zCMB, zHEL;
  // ------------ BEGIN ------------

  if ( !GENLC.CORRECT_HOSTMATCH) { return ; }

  // - - - - - - - - - - - - - - - - - - - - - 
  // check for transferring redshift to host redshift.
  // Here zHEL & zCMB both change
  if ( OVP_Z ) {
    
    if ( INPUTS.VEL_CMBAPEX > 0.0 ) {
      zCMB = zhelio_zcmb_translator(ZTRUE_noVPEC,RA,DEC,eq,+1);
    }
    else {
      zCMB = ZTRUE_noVPEC ; 
    }

    GENLC.REDSHIFT_HELIO = ZTRUE ;  // preserve this
    SNHOSTGAL.ZSPEC      = ZTRUE ; 
    GENLC.REDSHIFT_CMB   = zCMB ;   // store adjusted zCMB
    gen_distanceMag(zCMB, ZTRUE_noVPEC,
		    &GENLC.DLMU, &GENLC.LENSDMU ); // <== returned
  }

  // - - - - - - - - - - - - - - - - - - - - - 
  // check for switching to host coordinates; 
  // zCMB does not change, but zHEL changes.

  if ( OVP_RADEC && !OVP_Z ) {
    zCMB = GENLC.REDSHIFT_CMB  ; // preserve this
    if ( INPUTS.VEL_CMBAPEX > 0.0 ) 
      { zHEL = zhelio_zcmb_translator(zCMB,RA,DEC,eq,-1); }   
    else 
      { zHEL = zCMB ;  }

    gen_distanceMag(zCMB, zHEL, 
		    &GENLC.DLMU, &GENLC.LENSDMU ); // <== returned
    zHEL                 += zPEC;  // add zPEC after computing MU
    GENLC.REDSHIFT_HELIO  = zHEL ;     
    SNHOSTGAL.ZSPEC       = zHEL ;
  }

  return ;

} // end of TRANSFER_SNHOST_REDSHIFT

// ==============================
void GEN_SNHOST_GALMAG(int IGAL) {

  // Compute host mag contained within PSF for each observer-filter. 
  // Assume that each [filt]_obs value corresponds to the total
  // flux of the host; then use the generated SN position to 
  // determine several apertures and compute the flux-fraction 
  // in each aperture. Finally, compute host mags based on the
  // the total-host mag and the flux-fractions.
  //
  // For a given PSF, the effective aperture is
  // taken to be 4*PI*PSF^2.
  //
  // Fill struct
  //   SNHOSTGAL.GALMAG_TOT[ifilt_obs]
  //   SNHOSTGAL.GALMAG[ifilt_obs][iPSF]
  //
  // May 27, 2011: apply Galactic extinction, see MWXT[ifilt_obs]
  //
  // Dec 17, 2012: in NMSGPSF loop, include i=0 to store total mag.
  //
  // Mar 3 2015: compute local surface brightness SBFLUX and SBMAG
  //
  // May 5 2017: use user-input INPUTS.HOSTLIB_SBRADIUS
  //

  double 
     x_SN, y_SN
    ,xgal, ygal, MAGOBS, MAGOBS_LIB, PSF, FGAL, dF
    ,Rmin, Rmax, Rbin, R, Rcen, dm
    ,THmin, THmax, THbin, TH
    ,dRdTH, Jac
    ,GALFRAC_SUM[NMAGPSF_HOSTLIB+1]       // summed over Sersic profile
    ,sigFrac, RcenFrac
    ,GaussOvp[NMAGPSF_HOSTLIB+1] 
    ,AV, LAMOBS_AVG, MWXT[MXFILTINDX]
    ,RVMW = 3.1
    ;

  float lamavg4, lamrms4, lammin4, lammax4  ;
  int ifilt, ifilt_obs, i, IVAR, jbinTH, opt_frame    ;
  //  char fnam[] = "GEN_SNHOST_GALMAG" ;

  // ------------ BEGIN -------------


  for ( i=0; i <= NMAGPSF_HOSTLIB ; i++ ) {
    GALFRAC_SUM[i]   =  0.0 ;     // local 
  }

  // compute MilkyWay Galactic extinction in each filter using
  // extinction at mean wavelength
  opt_frame = GENFRAME_OBS ; 
  AV        = RVMW * GENLC.MWEBV_SMEAR ;
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs   = GENLC.IFILTMAP_OBS[ifilt];
    get_filtlam__(&opt_frame, &ifilt_obs, 
		  &lamavg4, &lamrms4, &lammin4, &lammax4 );
    LAMOBS_AVG = (double)lamavg4 ;
    MWXT[ifilt_obs] = GALextinct ( RVMW, AV, LAMOBS_AVG, 94 );

    // compute & store galaxy mag if they are defined
    IVAR      = HOSTLIB.IVAR_MAGOBS[ifilt_obs] ;

    if ( IVAR >= 0 ) {
      MAGOBS_LIB = HOSTLIB.VALUE_ZSORTED[IVAR][IGAL] ; 
      MAGOBS     = MAGOBS_LIB + MWXT[ifilt_obs] ; // dim  with Gal extinction
      SNHOSTGAL.GALMAG_TOT[ifilt_obs] = MAGOBS ;

      /* xxx
      printf(" xxx ------------------------ \n");
      printf(" xxx %s: load ifilt_obs=%d IGAL=%d, IVAR=%d \n",
	     fnam, ifilt_obs, IGAL, IVAR );
      printf(" xxx %s: MWXT=%.3f   MAGOBS=%.3f -> %.3f \n" ,
	     fnam, MWXT[ifilt_obs], MAGOBS_LIB, MAGOBS);
      xxxx */

    }
  }

  // --------------------------------------------------------
  // check option to compute host noise under SN
  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG) == 0 ) 
    { return ; }

  // get coords in galaxy frame; 
  // a,b  are along major/minor axes so we don't have to
  // worry about the relative rotation w.r.t. RA.
  // This simplification works with assumed circular PSF.

  x_SN = SNHOSTGAL.a_SNGALSEP_ASEC  ; 
  y_SN = SNHOSTGAL.b_SNGALSEP_ASEC  ;

  // max radius (w.r.t. SN) to integrate is twice the 
  // aperture radius,  or 2 x (2*PSFMAX).
  Rmin  = 0.0 ;
  Rmax  = HOSTLIB.Aperture_Rmax ;
  Rbin  = HOSTLIB.Aperture_Rbin ;

  THmin = 0.0 ;
  THmax = TWOPI * .99999 ;
  THbin = HOSTLIB.Aperture_THbin ;

  dRdTH = Rbin * THbin ;

  // start integration loop in polar coords around the SN.
  for ( R = Rmin; R < Rmax; R += Rbin ) {
    Rcen = R  + 0.5*Rbin ;  // center of R-bin w.r.t SN (arcsec)
    Jac  = Rcen * dRdTH ;   // Jacobian factor = r dr dtheta
    
    // store aperture-dependent stuff before THETA loop
    for ( i=1; i <= NMAGPSF_HOSTLIB ; i++ ) { 
      PSF         = HOSTLIB.Aperture_PSFSIG[i] ; // = 1/2 x Aperture_Radius
      RcenFrac    = Rcen / HOSTLIB.Aperture_Radius[i] ; 
      sigFrac     = PSF  / HOSTLIB.Aperture_Radius[i] ; 
      GaussOvp[i] = Gauss2d_Overlap(RcenFrac,sigFrac);   
    }

    jbinTH = 0;
    for ( TH = THmin; TH < THmax; TH += THbin ) {

      // Translate from SN polar coords to galaxy a,b coords
      // Use pre-computed cos(TH) and sin(TH) for speed
      xgal = x_SN + ( Rcen * HOSTLIB.Aperture_cosTH[jbinTH] ) ;
      ygal = y_SN + ( Rcen * HOSTLIB.Aperture_sinTH[jbinTH] ) ;
      FGAL = get_GALFLUX_HOSTLIB(xgal,ygal) ;
      jbinTH++ ;

      // increment GALFRAC for each PSF grid value
      for ( i=1; i <= NMAGPSF_HOSTLIB ; i++ ) { 
	dF       = Jac * FGAL * GaussOvp[i] ;
	GALFRAC_SUM[i] += dF ;
      }

    }  // end of TH loop

  }  // end of R loop


  // convert GALFRAC into mag shift and add to galaxy mag
  // i=0 is for the total mag; i=1-N is for aperture PSFs.


  for ( i=0; i <= NMAGPSF_HOSTLIB ; i++ ) {

    if ( GALFRAC_SUM[i] > 1.000 )  { GALFRAC_SUM[i] = 1.000 ; }

    // load global arrays
    SNHOSTGAL.GALFRAC[i]  = GALFRAC_SUM[i] ;

    if ( i == 0 ) 
      {  dm = 0.0 ; }
    else { 
      dm = -2.5*log10( GALFRAC_SUM[i] );  // aperture mag  - total galmag
    }

    /*
    printf(" GGG xxx GALID=%d GALFRAC[%d]=%f   R=%6.3f arcsec \n", 
	   SNHOSTGAL.GALID, i, GALFRAC_SUM[i], HOSTLIB.Aperture_Radius[i] );
    */

    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      IVAR      = HOSTLIB.IVAR_MAGOBS[ifilt_obs] ;

      if ( IVAR >= 0 ) {
	MAGOBS_LIB = HOSTLIB.VALUE_ZSORTED[IVAR][IGAL] ; 
	MAGOBS     = MAGOBS_LIB + MWXT[ifilt_obs] ; // dim  with Gal extinction
      }
      else {
	MAGOBS = 99.0 ; // allow missing filter in HOSTLIB 
      }

      SNHOSTGAL.GALMAG[ifilt_obs][i] = MAGOBS + dm ;
    }
  } // end i lopop over PSF bins


  // ----------------------
  // Mar 3 2015:
  // compute local surface brightness mag with effective PSF = 0.6''
  // Note that the effective area is 4*PI*PSF^2.

  double psfsig, arg, SB_MAG, SB_FLUX, AREA ;

  psfsig  = INPUTS.HOSTLIB_SBRADIUS/2.0; 
  AREA    = 3.14159*(4.0*psfsig*psfsig) ; // effective noise-equiv area

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs  = GENLC.IFILTMAP_OBS[ifilt];
    SB_MAG     = interp_GALMAG_HOSTLIB(ifilt_obs, psfsig ) ;
    SB_MAG    += 2.5*log10(AREA); // normalize to 1 sq-arcsec
    
    // convert mag to fluxcal
    arg     = -0.4*(SB_MAG - ZEROPOINT_FLUXCAL_DEFAULT) ;
    SB_FLUX = pow(10.0,arg) ;
    
    SNHOSTGAL.SB_FLUX[ifilt_obs] = SB_FLUX ;  
    SNHOSTGAL.SB_MAG[ifilt_obs]  = SB_MAG ;
  }

  return ;

} // end of GEN_SNHOST_GALMAG



// =====================================
double Gauss2d_Overlap(double offset, double sig) {

  // March 2011  R.Kessler
  //
  // Uses lookup table to quickly returns fraction (0-1) of 
  // 2d Gaussian integral that overlaps circle of radius 1.
  // 
  // Inputs:
  // offset : location of the Gaussian center w.r.t. 
  //          the center of the circle. offset=1 =>
  //          Gaussian center is on the perimeter.
  //
  // sig : sigma of Gaussian, assumed to be the same in x & y.
  //
  // For practical usage, pass dimensionless args
  //  offset = r/R
  //  sig    = sigma/R
  //
  // where R is the aperture radius.
  // Since aperture Radius = 2*PSFSIG, sig=0.5 is always
  // passed to this function, so this is really a 1-d interp
  // in 'offset' space. If we ever need to pass different
  // sig values, the Gauss2dIntegral.dat table will need
  // to be updated and the interp below needs to be upgraded.
  //

  double 
    ovp, ovp1, ovp2, ftmp, r1, xnsig, xnbin
    ,rmin, rmax, rbin, smin, smax, sbin 
    ;

  int Nr, Ns, ibin_r, ibin_s, IBIN1_TABLE, IBIN2_TABLE ;
  int MODETEST = 0 ;

  char fnam[] =  "Gauss2d_Overlap" ;
  
  // -------------- BEGIN ---------------

  ovp = 0.0 ;

  // return zero if offfset is more than 6 sigma away 
  // from center of  aperture
  xnsig = offset/sig ;
  if ( xnsig > 6.0 ) { return 0.0 ; }

  // -----------------------
  // check test-option to use delta-function instead of Gaussian
  if  ( MODETEST ) {
    if ( offset <= 1.0 ) 
      { ovp = 1.0 ; }
    else 
      { ovp = 0.0 ; }

    return ovp ;
  }

  // ------------------------

  rbin = HOSTLIB.Gauss2dRadius[0];
  rmin = HOSTLIB.Gauss2dRadius[1];
  rmax = HOSTLIB.Gauss2dRadius[2];
  Nr   = HOSTLIB.NBIN_Gauss2dRadius;

  sbin = HOSTLIB.Gauss2dSigma[0];
  smin = HOSTLIB.Gauss2dSigma[1];
  smax = HOSTLIB.Gauss2dSigma[2];
  Ns   = HOSTLIB.NBIN_Gauss2dSigma;

  // check for crazy values
  if ( offset < 0.0 || offset > rmax ) {
    sprintf(c1err,"invalid offset = %f ", offset);
    sprintf(c2err,"Valid offset (r/R) range is %5.2f - %5.2f", rmin,rmax);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( sig < 0.0 || sig > smax ) {
    sprintf(c1err,"invalid sigma = %f ", sig);
    sprintf(c2err,"Valid sigma/R  range is %5.2f - %5.2f", smin, smax);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( Nr > 1 )     
    { xnbin = (offset-rmin)/rbin ;  ibin_r = (int)xnbin  ; }
  else
    { ibin_r = 0 ; }

  if ( Ns > 1 ) 
    { xnbin = (sig-smin)/sbin ; ibin_s = (int)xnbin ; }
  else
    { ibin_s = 0; }

  // interpolate in r-space only since for now there
  // is just one sigma-bin.

  IBIN1_TABLE = ibin_s + (ibin_r+0)*Ns ; // nearest table bin
  IBIN2_TABLE = ibin_s + (ibin_r+1)*Ns ; 

  if ( IBIN1_TABLE < 0 || 
       IBIN2_TABLE < 0 ||
       IBIN1_TABLE > HOSTLIB.NGauss2d  ||
       IBIN2_TABLE > HOSTLIB.NGauss2d        
       ) {

    printf("\n PRE-ABORT DUMP: (offset=%f, sig=%f) \n", offset, sig);
    printf("\t rmin=%.3f rbin=%.3f  ibin_r=%d \n",
	   rmin, rbin, ibin_r);
    printf("\t smin=%.3f sbin=%.3f  ibin_s=%d \n",
	   smin, sbin, ibin_s );

    sprintf(c1err,"Invalid IBIN[1,2]_TABLE = %d,%d  (valid range: 1 to %d)",
	    IBIN1_TABLE, IBIN2_TABLE, HOSTLIB.NGauss2d );
    sprintf(c2err,"offset=%6.3f sig=%6.3f  ibin_[r,s]=%d,%d",
	    offset, sig, ibin_r, ibin_s );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // get r-interp fraction
  r1   = HOSTLIB.Gauss2dTable[0][IBIN1_TABLE]; // radius
  ftmp = (offset - r1)/rbin ;
  if ( ftmp < -1.0E-9 || ftmp > 1.000000001 ) {
    sprintf(c1err,"r-interp problem: ftmp=%le", ftmp);
    sprintf(c2err,"offset=%f sig=%f r1=%f IBIN1_TABLE=%d",
	    offset, sig, r1, IBIN1_TABLE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  ovp1 = HOSTLIB.Gauss2dTable[2][IBIN1_TABLE]; // integral
  ovp2 = HOSTLIB.Gauss2dTable[2][IBIN2_TABLE];
  ovp  = ovp1 + ftmp*(ovp2 - ovp1);

  /*
  printf(" xxxx offset=%6.3f  sig=%6.3f  IBIN_TABLE=%2.2d-%2.2d  ovp=%f \n",
	 offset, sig, IBIN1_TABLE, IBIN2_TABLE, ovp );
  */

  return ovp ;

} // end of Gauss2d_Overlap



// ===================================
double get_GALFLUX_HOSTLIB(double xgal, double ygal) {

  //
  // Return relative flux at this galaxy location (arcsec)
  // xgal and ygal are w.r.t. galaxy center;
  // xgal is along major axis, ygal along minor axis.
  //
  // Mar 4 2015: fix aweful index bug setting FGAL_TOT.
  //             Bug affects only the host-noise contribution.
  //

  int  j, NBIN ;

  double 
    a, b, w, n, bn, rexp, sqsum,  arg
    ,reduced_R, xx, yy, FSUM_PROFILE, F, FGAL_TOT
    ;

  //  char fnam[] = "get_GALFLUX_HOSTLIB" ;

  // ---------------- BEGIN ---------------

  FSUM_PROFILE = 0.0 ;
  NBIN = SERSIC_TABLE.NBIN_reduced ;

  for ( j=0; j < SERSIC_PROFILE.NDEF; j++ ) {
    
    // strip off info for this galaxy component.
    a   = SNHOSTGAL.SERSIC_a[j] ;
    b   = SNHOSTGAL.SERSIC_b[j] ;
    n   = SNHOSTGAL.SERSIC_n[j] ;
    w   = SNHOSTGAL.SERSIC_w[j] ;
    bn  = SNHOSTGAL.SERSIC_bn[j] ;
    rexp = 1./n ;

    // Flux normalization = total flux over galaxy.
    // The stored INTEG_SUM is integrated over the reduced radius;
    // the pre-factor a*b converts to the total flux over the 
    // physical galaxy size.
    FGAL_TOT = (a*b) * SERSIC_TABLE.INTEG_SUM[NSERSIC_TABLE-1];

    // get scale needed to move input xgal,ygal onto the
    // half-light ellipse; this scale is the reduced radius (R/Re)
    xx        = xgal/a ;
    yy        = ygal/b ;
    sqsum     = xx*xx + yy*yy ; 
    reduced_R = sqrt(sqsum) ;

    // calculate Flux for this Sersic component
    arg  = pow(reduced_R,rexp) - 1.0 ; 
    F    = (w/FGAL_TOT) * exp(-bn*arg) ;

    FSUM_PROFILE += F ;
  }

  return(FSUM_PROFILE) ;

} // end of get_GALFLUX_HOSTLIB

// =============================================
void  STORE_SNHOST_MISC(int IGAL) {

  // Feb 11 2014.
  // store misc variables in SNHOSTGAL structure.
  // Needed to write to data file.

  int NDIM = HOSTLIB_WGTMAP.GRIDMAP.NDIM;
  int ivar, IVAR_STORE, IVAR ;
  char fnam[] = "STORE_SNHOST_MISC" ;

  // ------------ BEGIN ---------
  
  // store wgt and SN mag shift for this host
  SNHOSTGAL.WGTMAP_WGT        = HOSTLIB_WGTMAP.WGT[IGAL] ;
  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SNMAGSHIFT )
    {  SNHOSTGAL.WGTMAP_SNMAGSHIFT = HOSTLIB_WGTMAP.SNMAGSHIFT[IGAL]; }
  else
    {  SNHOSTGAL.WGTMAP_SNMAGSHIFT = 0.0 ; }


  for ( ivar=0; ivar < NDIM; ivar++ ) {
    IVAR_STORE = HOSTLIB.IVAR_STORE[ivar];
    SNHOSTGAL.WGTMAP_VALUES[ivar] = HOSTLIB.VALUE_ZSORTED[IVAR_STORE][IGAL] ;
  }

  // store host mass if it's there. 
  IVAR  = HOSTLIB.IVAR_LOGMASS ;
  if ( IVAR > 0 ) 
    { SNHOSTGAL.LOGMASS = HOSTLIB.VALUE_ZSORTED[IVAR][IGAL] ; }

  IVAR = HOSTLIB.IVAR_LOGMASS_ERR ;
  if ( IVAR > 0 ) 
    { SNHOSTGAL.LOGMASS_ERR = HOSTLIB.VALUE_ZSORTED[IVAR][IGAL] ; }


  // -----------------------------------------------------------
  // Sep 16 2015
  // set GENLC.FIELDNAME if FIELD column is present in hostlib
  int   ep ;
  char *FIELD ;
  if ( HOSTLIB.IVAR_FIELD > 0 ) {
    FIELD = HOSTLIB.FIELD_ZSORTED[IGAL] ;
    if ( strlen(FIELD) == 0 ) {
      sprintf(c1err,"FIELD='' for IGAL=%d  GALID=%lld", 
	      IGAL, get_GALID_HOSTLIB(IGAL) ) ;
      sprintf(c2err,"Check HOSTLIB");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    for(ep=0; ep < MXEPSIM; ep++ ) {    
      sprintf(GENLC.FIELDNAME[ep], "%s", FIELD);
    }
  }

  return ;

} // end of STORE_SNHOST_MISC


// =============================================
void LOAD_OUTVAR_HOSTLIB(int IGAL) {

  // Feb 2014.
  // load HOSTLIB_OUTVAR_EXTRA.VALUE array for output data file.
  // These variables are user-defined by sim-input key
  // HOSTLIB_STOREVAR:  <var1>,<var2>,<var3>, ...

  int NVAR_OUT, ivar, IVAR_STORE ;
  double DVAL ;
  //  char fnam[] = "LOAD_OUTVAR_HOSTLIB" ;

  // ------------- BEGIN ------------

  NVAR_OUT = HOSTLIB_OUTVAR_EXTRA.NOUT ;
  if ( NVAR_OUT == 0 ) { return ; }

  for(ivar=0; ivar < NVAR_OUT; ivar++ ) {
    IVAR_STORE = HOSTLIB_OUTVAR_EXTRA.IVAR_STORE[ivar] ;
    DVAL       = HOSTLIB.VALUE_ZSORTED[IVAR_STORE][IGAL] ;
    HOSTLIB_OUTVAR_EXTRA.VALUE[ivar] = DVAL ;   

    //    printf(" xxx IGAL=%d  ivar=%d IVAR_STORE=%d DVAL=%f \n", 
    //	   IGAL, ivar, IVAR_STORE, DVAL);
  }

} // end of LOAD_OUTVAR_EXTRA_HOSTLIB


// ===============================================
double modelPar_from_SNHOST(double PARVAL_ORIG, char *PARNAME) {

  // Created Jan 28 2014 by R.Kessler
  // Return SN parameter (shape, color, RV, beta ...)
  // corresponding to input *parName. If *parName
  // is not specified in the hostlib, return PARVAL_ORIG
  // i.e., using the SN parameter from the HOSTLIB is optional.

  // Input: PARVAL_ORIG = parameter selection in sim before HOSLOT
  // Input: PARNAME     = name of SN parameter
  //
  // Mar 23 2018: allow SNMAGSHIFT as well.

  int GALID          = SNHOSTGAL.GALID ;
  int IGAL           = SNHOSTGAL.IGAL ;
  int USE_GAMMA_GRID = HOSTLIB_WGTMAP.USE_SALT2GAMMA_GRID;
  int IVAR, USE1, USE2 ;
  int noABORT = 0 ;
  double PARVAL_OUT ;
  //  char fnam[] = "modelPar_from_SNHOST" ;

  // ----------------- BEGIN ----------------
  
  PARVAL_OUT = PARVAL_ORIG ; // default output valie 

  // check for GAMMA_GRID (Jun 25 2019)
  if ( USE_GAMMA_GRID ) {
    if ( strcmp(PARNAME,HOSTLIB_VARNAME_SNMAGSHIFT)== 0 ) {
      PARVAL_OUT = snmagshift_salt2gamma_HOSTLIB(GALID); 
      return(PARVAL_OUT);
    }
  }

  // if USESNPAR option is not set, then return PARVAL_ORIG immediately.
  USE1   = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) ;
  USE2   = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SNMAGSHIFT) ;
  if ( USE1==0 && USE2==0 ) { return PARVAL_OUT ; }

  // get the vaue
  IVAR   = IVAR_HOSTLIB(PARNAME, noABORT );

  if ( IVAR > 0 ) {
    IGAL        = SNHOSTGAL.IGAL ;
    PARVAL_OUT  = HOSTLIB.VALUE_ZSORTED[IVAR][IGAL] ;
  }


  return PARVAL_OUT ;

} // end of  modelPar_from_SNHOST


// ======================================
void DEBUG_1LINEDUMP_SNHOST(void) {

  // dump 1 line-dump for easy readback.
  // <job.exe args> | grep 222222 > host.dump

  printf(" %d  %6.4f %6.4f %4.2f  %7.4f %7.4f %5.3f %6.2f  %f %f \n"    
	 , 222222    // for easy grep
	 , SNHOSTGAL.SERSIC_a[1]
	 , SNHOSTGAL.SERSIC_b[1]
	 , SNHOSTGAL.SERSIC_INDEX
	 , SNHOSTGAL.RA_SNGALSEP_ASEC 
	 , SNHOSTGAL.DEC_SNGALSEP_ASEC
	 , SNHOSTGAL.reduced_R
	 , SNHOSTGAL.phi/.0174533
	 , SNHOSTGAL.GALFRAC[2] 
	 , SNHOSTGAL.GALFRAC[4] 
	 );
  fflush(stdout);

} // end of DEBUG_1LINEDUMP_SNHOST

// ======================
void DUMP_SNHOST(void) {

  // Feb 4 2019: add more stuff (DLR, DDLR, SNSEP, RA,DEC ...)

  double mag, a, b, w, n, bn, fgal, psf, ZDIF;
  int    IVAR_ZPHOT, ir, IVAR, ifilt, ifilt_obs, j, IGAL ;
  char cfilt[2];
  char fnam[] = "DUMP_SNHOST" ;

  // -------------- BEGIN ---------

  ZDIF = SNHOSTGAL.ZGEN - SNHOSTGAL.ZTRUE ;
  IGAL = SNHOSTGAL.IGAL;

  printf("# -------------------------------------------------------- \n");
  printf(" %s for CID = %d \n", fnam, GENLC.CID );
  printf("    SN-GENERATED z=%.4f, RA=%.6f, DEC=%.6f deg (FIELD=%s).\n",
	 SNHOSTGAL.ZGEN, GENLC.RA, GENLC.DEC, SIMLIB_HEADER.FIELD );

  printf("\t => IGAL-range=%d-%d  IGAL=%d  ran(GALID)=%4.3f \n"
	 , SNHOSTGAL.IGAL_SELECT_RANGE[0]
	 , SNHOSTGAL.IGAL_SELECT_RANGE[1]
	 , SNHOSTGAL.IGAL
	 , SNHOSTGAL.FlatRan1_GALID );
  
  printf("\t => GALID=%10lld  ZTRUE=%5.4f  ZGEN-ZTRUE= %f \n"
	 ,SNHOSTGAL.GALID
	 ,SNHOSTGAL.ZTRUE, ZDIF );

  fflush(stdout);

  IVAR_ZPHOT = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZPHOT,0);
  if ( IVAR_ZPHOT > 0 ) {
    printf("\t => ZPHOT = %6.4f +- %6.4f \n"
	   ,SNHOSTGAL.ZPHOT
	   ,SNHOSTGAL.ZPHOT_ERR );
    fflush(stdout);
  }

  if ( SERSIC_PROFILE.NDEF > 0 ) {

    printf("\t => GAL(RA,DEC) = %.6f, %.6f deg. \n",
	   SNHOSTGAL.RA_GAL_DEG,  SNHOSTGAL.DEC_GAL_DEG  );
    printf("\t => SN(RA,DEC)  = %.6f, %.6f deg. \n",
	   SNHOSTGAL.RA_SN_DEG,  SNHOSTGAL.DEC_SN_DEG  );
    printf("\t => SN-Gal sep(RA,DEC) : %6.3f , %6.3f arcsec \n"
	   , SNHOSTGAL.RA_SNGALSEP_ASEC
	   , SNHOSTGAL.DEC_SNGALSEP_ASEC );
    printf("\t => SNSEP=%.3f arcsec, DLR=%.3f  arcsec, d_DLR=%.3f \n", 
	   SNHOSTGAL.SNSEP, SNHOSTGAL.DLR, SNHOSTGAL.DDLR );

    fflush(stdout);

    for ( j=0; j < SERSIC_PROFILE.NDEF; j++ ) {
      a  =  SNHOSTGAL.SERSIC_a[j] ;
      b  =  SNHOSTGAL.SERSIC_b[j] ;
      w  =  SNHOSTGAL.SERSIC_w[j] ;
      n  =  SNHOSTGAL.SERSIC_n[j] ;
      bn =  SNHOSTGAL.SERSIC_bn[j] ;
      printf("\t => Sersic a_half=%7.4f  b_half=%7.4f  "
	     " n=%5.2f (bn=%5.3f)  wgt=%5.3f \n",
	     a,b, n, bn, w );
      fflush(stdout);
    }
  }

  printf("\t => a_sep=%.4f  b_sep=%.4f  phi=%.4f  R_red=%.4f \n",
	 SNHOSTGAL.a_SNGALSEP_ASEC, SNHOSTGAL.b_SNGALSEP_ASEC,
	 SNHOSTGAL.phi, SNHOSTGAL.reduced_R  );


  // ----------------------------------------------------
  // weight map parameters, wgt and SN mag shift
  printf("\t => WGT = %f   SNMAGSHIFT=%f \n",
	 SNHOSTGAL.WGTMAP_WGT, SNHOSTGAL.WGTMAP_SNMAGSHIFT );
  fflush(stdout);

  // ----------------------------------------------------
  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG ) {

    //  first print the PSFSIG grid

    // first print the gal-fractions
    printf("\t => PSFSIG interp-grid: ");      
    fflush(stdout);
    for ( ir  = 1; ir <= NMAGPSF_HOSTLIB; ir++ ) {
      psf = HOSTLIB.Aperture_PSFSIG[ir] ;
      printf("%7.4f ", psf);
    }
    printf("\n");
    fflush(stdout);


    // first print the gal-fractions
    printf("\t => Gal-flux fractions: ");      
    for ( ir  = 1; ir <= NMAGPSF_HOSTLIB; ir++ ) {
      fgal = SNHOSTGAL.GALFRAC[ir] ;
      printf("%7.4f ", fgal);
    }
    printf("\n");   fflush(stdout);

      // next the aperture mags under SN
    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );

      IVAR      = HOSTLIB.IVAR_MAGOBS[ifilt_obs] ;
      if ( IVAR <= 0 ) { continue ; }

      printf("\t =>  %s aperture mags  : ", cfilt );

      for ( ir  = 1; ir <= NMAGPSF_HOSTLIB; ir++ ) {
	mag = SNHOSTGAL.GALMAG[ifilt_obs][ir] ;
	printf("%7.4f ", mag);
      }
      printf("\n");   fflush(stdout);
    } // end of ifilt loop
    
  } // end of GALMAG if-block


} // end of DUMP_SNHOST


// ========================================
void setbit_HOSTLIB_MSKOPT(int MSKOPT) {

  // set MSKOPT bit; if already set do nothing.
  // Note that MSKOPT is a mask with a bit set (not a bit number)

  int USE   = (INPUTS.HOSTLIB_MSKOPT & MSKOPT) ;
  if ( USE <= 0 ) { INPUTS.HOSTLIB_MSKOPT += MSKOPT ; }
  return ;

}  // end of setbit_HOSTLIB_MSKOPT


// =================================
int fetch_HOSTPAR_GENMODEL(int OPT, char *NAMES_HOSTPAR, double*VAL_HOSTPAR) {

  // Created July 2019
  // Inputs:
  //   OPT = 1 --> store NAMES_HOSTPAR
  //   OPT = 2 --> store VAL_HOSTPAR
  //
  // Output:
  //    NAMES_HOSTPAR : comma-separated list of par names (filled if OPT=1)
  //    VAL_HOSTPAR   : array of hostpar values (filled if OPT=2)
  //
  // Function returns number of HOSTPAR params

  int  NPAR=0, ivar ;
  int  NVAR_WGTMAP = HOSTLIB_WGTMAP.GRIDMAP.NDIM ;
  int  NVAR_EXTRA  = HOSTLIB_OUTVAR_EXTRA.NOUT ;
  char comma[] = ",";
  char fnam[] = "fetch_HOSTPAR_GENMODEL";

  // ----------- BEGIN -----------

  if ( OPT == 1 ) {
    // always start with RV and AV    
    sprintf(NAMES_HOSTPAR,"RV,AV"); NPAR=2;

    for ( ivar=0; ivar < NVAR_WGTMAP; ivar++ ) {  
      strcat(NAMES_HOSTPAR,comma);
      strcat(NAMES_HOSTPAR,HOSTLIB_WGTMAP.VARNAME[ivar]); 
      NPAR++ ;

    }

    // tack on user-defined variables from sim-input HOSTLIB_STOREPAR key
    for(ivar=0; ivar < NVAR_EXTRA; ivar++ ) {
      if ( HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[ivar] ) { continue; }
      strcat(NAMES_HOSTPAR,comma);
      strcat(NAMES_HOSTPAR,HOSTLIB_OUTVAR_EXTRA.NAME[ivar]);
      NPAR++ ;
    }
    // add anything used in WGTMAP or HOSTLIB_STOREPAR
  }
  else if ( OPT ==  2 ) {
    VAL_HOSTPAR[0] = GENLC.RV;
    VAL_HOSTPAR[1] = GENLC.AV;
    NPAR=2;

    for ( ivar=0; ivar < NVAR_WGTMAP; ivar++ ) {  
      VAL_HOSTPAR[NPAR] = SNHOSTGAL.WGTMAP_VALUES[ivar] ;
      NPAR++ ;
    }

    for(ivar=0; ivar < NVAR_EXTRA; ivar++ ) {
      if ( HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[ivar] ) { continue; }
      VAL_HOSTPAR[NPAR] = HOSTLIB_OUTVAR_EXTRA.VALUE[ivar];
      NPAR++ ;
    }

  }
  else {
    sprintf(c1err,"Invalid OPT=%d", OPT);
    sprintf(c2err,"Allowed OPT values, 1 or 2");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  

  return(NPAR);

} // end fetch_HOSTPAR_GENMODEL

// ===================================
void rewrite_HOSTLIB(HOSTLIB_APPEND_DEF *HOSTLIB_APPEND) {

  // generic utility to rewrite HOSTLIB using information
  // Passed here vias HOSTLIB_APPEND.
  // 

  char *SUFFIX       = HOSTLIB_APPEND->FILENAME_SUFFIX; // or new HOSTLIB
  int  NLINE_COMMENT = HOSTLIB_APPEND->NLINE_COMMENT;
  char *VARNAMES     = HOSTLIB_APPEND->VARNAMES_APPEND ;

  FILE *FP_ORIG, *FP_NEW;
  char *HLIB_ORIG = INPUTS.HOSTLIB_FILE;
  char  HLIB_TMP[MXPATHLEN], HLIB_NEW[100], DUMPATH[MXPATHLEN];
  char fnam[] = "rewrite_HOSTLIB" ;

  // -------------- BEGIN --------------

  // create name of new hostlib file
  sprintf(HLIB_TMP,"%s%s", HLIB_ORIG, SUFFIX ); 

  // remove path from HLIB_NEW to ensure that new file is created
  // locally and not in somebody else's directory.
  extract_MODELNAME(HLIB_TMP, DUMPATH, HLIB_NEW);

  FP_ORIG = fopen(HLIB_ORIG,"rt");
  FP_NEW  = fopen(HLIB_NEW, "wt");
  if ( !FP_ORIG ) {
    sprintf(c1err,"Could not open original HOSTLIB_FILE");
    sprintf(c2err,"'%s' ", HLIB_ORIG);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  if ( !FP_NEW ) {
    sprintf(c1err,"Could not open new HOSTLIB_FILE");
    sprintf(c2err,"'%s' ", HLIB_NEW);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n");
  printf("  Created '%s' \n", HLIB_NEW);

  int iline;
  for(iline=0; iline < NLINE_COMMENT; iline++ ) 
    { fprintf(FP_NEW,"# %s\n", HOSTLIB_APPEND->COMMENT[iline] ); }

  fprintf(FP_NEW,"# \n");
  fprintf(FP_NEW,"# Below are original HOSTLIB comments and table\n");
  fprintf(FP_NEW,"# with NBR_LIST column appended.\n");
  fprintf(FP_NEW,"# - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(FP_NEW,"# \n");

  // - - - - - - - - - - - - - - - - - 
  // read each original line

  long long GALID, GALID_orig ;
  int   igal_unsort, igal_zsort, ivar, NWD_LINE;
  char *LINE, *LINE_APPEND, *FIRSTWORD, *NEXTWORD, *ptrCR ;

  LINE         = (char*) malloc ( sizeof(char) * MXCHAR_LINE_HOSTLIB );
  LINE_APPEND  = (char*) malloc ( sizeof(char) * 100 );
  FIRSTWORD    = (char*) malloc ( sizeof(char) * 100 );
  NEXTWORD     = (char*) malloc ( sizeof(char) * 100 );

  igal_unsort = 0;
  while ( fgets(LINE, MXCHAR_LINE_HOSTLIB, FP_ORIG) != NULL ) {

    NWD_LINE = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
    LINE_APPEND[0] = 0 ;

    if ( NWD_LINE > 2 ) {
      get_PARSE_WORD(0, 0, FIRSTWORD);
      if ( strcmp(FIRSTWORD,"VARNAMES:") == 0 ) 
	{ sprintf(LINE_APPEND,"%s", VARNAMES); }

      else if ( strcmp(FIRSTWORD,"GAL:") == 0 ) {

	// make sure GALID matches
	ivar       = HOSTLIB.IVAR_GALID;
	igal_zsort = HOSTLIB.LIBINDEX_ZSORT[igal_unsort];
	GALID      = (long long)HOSTLIB.VALUE_ZSORTED[ivar][igal_zsort] ;

	get_PARSE_WORD(0, 1, NEXTWORD); // read GALID
	sscanf(NEXTWORD, "%lld", &GALID_orig);
	if ( GALID != GALID_orig ) {
	  sprintf(c1err,"GALID mis-match for igal_unsort=%d", igal_unsort);
	  sprintf(c2err,"GALID(orig)=%lld, but stored GALID=%lld",
		  GALID_orig, GALID ) ;
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	}

	sprintf(LINE_APPEND, HOSTLIB_APPEND->LINE_APPEND[igal_unsort]);

	/* xxx
	for(ifilt=1; ifilt <= NFILT; ifilt++ ) {
	  sprintf(cval, " %6.3f", MAG_STORE[ifilt][igal_unsort] );
	  strcat(LINE_APPEND,cval);
	}
	xxxxxxx */

	igal_unsort++ ;
      }
    }

    ptrCR = strchr(LINE,'\n'); if(ptrCR){*ptrCR=' ';} // remove <CR>
    fprintf(FP_NEW,"%s %s\n", LINE, LINE_APPEND);
  }


  fclose(FP_ORIG); fclose(FP_NEW);

  return ;

} // end rewrite_HOSTLIB

// =======================================
void malloc_HOSTLIB_APPEND(int NGAL, HOSTLIB_APPEND_DEF *HOSTLIB_APPEND) {

  // malloc and init
  int i;
  char fnam[] = "malloc_HOSTLIB_APPEND" ;

  // ---------------- BEGIN ---------------

  HOSTLIB_APPEND->NLINE_APPEND = NGAL;
  HOSTLIB_APPEND->LINE_APPEND = (char**) malloc( NGAL*sizeof(char*) );

  for(i=0; i < NGAL; i++ ) {
    HOSTLIB_APPEND->LINE_APPEND[i] = (char*) malloc( 100*sizeof(char*) );
    sprintf(HOSTLIB_APPEND->LINE_APPEND[i],"NULL_APPEND");
  }

  HOSTLIB_APPEND->NLINE_COMMENT = 0;

  return;

} // end malloc_HOSTLIB_APPEND

// =========================================
void addComment_HOSTLIB_APPEND(char *COMMENT, 
			       HOSTLIB_APPEND_DEF *HOSTLIB_APPEND) {

  int NL   = HOSTLIB_APPEND->NLINE_COMMENT;
  int MEMC = 120 * sizeof(char);
  
  HOSTLIB_APPEND->COMMENT[NL] = (char*) malloc(MEMC);
  sprintf(HOSTLIB_APPEND->COMMENT[NL],"%s", COMMENT);
  HOSTLIB_APPEND->NLINE_COMMENT++ ;

  return;

}  // end add_HOSTLIB_APPEND_COMMENT

// ===================================
void rewrite_HOSTLIB_plusMags(void) {
  
  // If +HOSTMAGS option is given on command-line, this function
  // is called to compute synthetic mags and re-write the HOSTLIB
  // with [band]_obs columns. 
  // Beware that 'igal' is a redshift-sorted index for the other 
  // HOSTLIB functions, but here we use the original HOSTLIB order.
  // 

  int NGAL     = HOSTLIB.NGAL_STORE ;
  int NFILT    = NFILT_SEDMODEL ; 
  int NBIN_LAM = INPUTS_SPECTRO.NBIN_LAM ;
  int MEMD     = sizeof(double) * NBIN_LAM ;

  int igal_unsort, igal_zsort, ifilt, ifilt_obs, ivar, DUMPFLAG=0 ;
  long long GALID, GALID_orig ;
  double ZTRUE, MWEBV=0.0, mag, *GENFLUX_LIST, *GENMAG_LIST;
  float  *MAG_STORE ;

  HOSTLIB_APPEND_DEF HOSTLIB_APPEND ;
  char cval[20], LINE_APPEND[100] ;
  char fnam[] = "rewrite_HOSTLIB_plusMags" ;

  // internal debug
  long long GALID_DUMP = 205925 ;
  int  NGAL_DEBUG      = -500; 

  // ------------ BEGIN -----------

  print_banner(fnam);

  printf("\t NGAL[READ,STORE] = %d, %d \n", 
	 HOSTLIB.NGAL_READ, HOSTLIB.NGAL_STORE);

  // ----------------------------

  // malloc local arrays
  GENFLUX_LIST = (double*) malloc(MEMD);
  GENMAG_LIST  = (double*) malloc(MEMD);

  MAG_STORE = (float*) malloc( (NFILT+1)*sizeof(float) );

  for(ifilt=0; ifilt <= NFILT; ifilt++ ) 
    { HOSTSPEC.NWARN_INTEG_HOSTMAG[ifilt] = 0 ; }


  malloc_HOSTLIB_APPEND(NGAL, &HOSTLIB_APPEND);

  // ----------------------------  
  if ( NGAL_DEBUG > 0 ) { NGAL = NGAL_DEBUG; }

  for(igal_unsort=0; igal_unsort < NGAL; igal_unsort++ ) {

    igal_zsort = HOSTLIB.LIBINDEX_ZSORT[igal_unsort];
    SNHOSTGAL.IGAL = igal_zsort; // needed by genSpec_HOSTLIB

    ivar  = HOSTLIB.IVAR_GALID;
    GALID = (long long)HOSTLIB.VALUE_ZSORTED[ivar][igal_zsort] ;

    ivar  = HOSTLIB.IVAR_ZTRUE ;
    ZTRUE = HOSTLIB.VALUE_ZSORTED[ivar][igal_zsort] ;
    
    DUMPFLAG = ( GALID == GALID_DUMP );    

    genSpec_HOSTLIB(ZTRUE,          // (I) helio redshift
		    MWEBV,          // (I) Galactic extinction
		    DUMPFLAG,       // (I) dump flag
		    GENFLUX_LIST,   // (O) fluxGen per bin 
		    GENMAG_LIST );  // (O) magGen per bin

    // ignore GENFLUX_LIST & GENMAG_LIST returned by genSpec_HOSTLIB.
    // Instead, integmag_hostSpec below uses global HOSTSPEC.FLAM_EVT 
    // that is loaded in genSpec_HOSTLIB.

    LINE_APPEND[0] = 0;
    for ( ifilt=1; ifilt <= NFILT; ifilt++ ) {
      ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs;
      mag       = integmag_hostSpec(ifilt_obs,ZTRUE,DUMPFLAG);
      MAG_STORE[ifilt] = mag;
      sprintf(cval, " %6.3f", MAG_STORE[ifilt] );
      strcat(LINE_APPEND,cval);
    }
    sprintf(HOSTLIB_APPEND.LINE_APPEND[igal_unsort],"%s", LINE_APPEND);

    if ( (igal_unsort % 10000) == 0 ) {
      printf("\t Processing igal %8d of %8d \n", igal_unsort, NGAL);
      fflush(stdout);
    }

    // xxxxxxxxxxxxxxxxxxxxxxxxxxx
    if ( ZTRUE < -1.0 ) {
      printf(" xxx ------------------------------------- \n");
      printf(" xxx igal_unsort=%2d  GALID=%8d  ZTRUE=%.5f \n", 
	     igal_unsort, GALID, ZTRUE);

      printf("     mag(%s) = ", INPUTS.GENFILTERS);
      for ( ifilt=1; ifilt <= NFILT; ifilt++ ) 
	{ printf("%7.2f", MAG_STORE[ifilt] );   }
      printf("\n"); fflush(stdout);
    } // end igal
    // xxxxxxxxxxxxxxxxxxxx

  } // end igal


  // ------------------------
  // set HOSTLIB_APPEND struct, then rewrite hostlib.

  int L  ;
  char VARNAMES_HOSTMAGS[200], varname_mag[40], *cfilt, msg[80];
  // create additional varnames of mags to append to VARNAMES list
  VARNAMES_HOSTMAGS[0] = 0;
  for(ifilt=1; ifilt <= NFILT; ifilt++ ) {
    cfilt = FILTER_SEDMODEL[ifilt].name;     L = strlen(cfilt);
    sprintf(varname_mag," %c_obs", cfilt[L-1] );
    strcat(VARNAMES_HOSTMAGS,varname_mag);
  }

  sprintf(HOSTLIB_APPEND.FILENAME_SUFFIX,"+HOSTMAGS");
  sprintf(HOSTLIB_APPEND.VARNAMES_APPEND,"%s", VARNAMES_HOSTMAGS);

  // construct comment lines for top of new HOSTLIB
  sprintf(msg,"Append synthetic host mags computed from host spectra.");
  addComment_HOSTLIB_APPEND(msg,&HOSTLIB_APPEND);
				   
  // execute re-write
  rewrite_HOSTLIB(&HOSTLIB_APPEND);

  /* xxxxxxxxxxxxx mark delete Nov 7 2019 xxxxxxxxxxxx
  // --------------------------------------------------------
  // create new HOSTLIB with host mags appended to varnames.
  // Read original hostlib and copy each line so that format
  // is not changed; then append host mags
  char *HLIB_ORIG = INPUTS.HOSTLIB_FILE,  HLIB_NEW[MXPATHLEN];
  char LINE[MXCHAR_LINE_HOSTLIB] ;
  char *ptrCR ;
  char FIRSTWORD[100], NEXTWORD[100] ;
  char tmpFile[MXPATHLEN], NULLPATH[] = "" ;
  FILE *FP_ORIG, *FP_NEW;
  int  NWD_LINE, gzipFlag ;
  sprintf(HLIB_NEW,"%s+HOSTMAGS_OLD", HLIB_ORIG);

  //FP_ORIG = snana_openTextFile(1, NULLPATH, HLIB_ORIG, tmpFile, &gzipFlag);
  FP_ORIG = fopen(HLIB_ORIG,"rt");
  FP_NEW  = fopen(HLIB_NEW, "wt");
  if ( !FP_ORIG ) {
    sprintf(c1err,"Could not open original HOSTLIB_FILE");
    sprintf(c2err,"'%s' ", HLIB_ORIG);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  if ( !FP_NEW ) {
    sprintf(c1err,"Could not open new HOSTLIB_FILE");
    sprintf(c2err,"'%s' ", HLIB_NEW);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n");
  printf("  Created '%s' with synthetic host mags\n", HLIB_NEW);
  fprintf(FP_NEW,
	  "# Append synthetic host mags computed from host spectra.\n\n");

  igal_unsort = 0;
  while ( fgets(LINE, MXCHAR_LINE_HOSTLIB, FP_ORIG) != NULL ) {

    NWD_LINE = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
    LINE_APPEND[0] = 0;

    if ( NWD_LINE > 2 ) {
      get_PARSE_WORD(0, 0, FIRSTWORD);
      if ( strcmp(FIRSTWORD,"VARNAMES:") == 0 ) 
	{ sprintf(LINE_APPEND,"%s", VARNAMES_HOSTMAGS); }
      else if ( strcmp(FIRSTWORD,"GAL:") == 0 ) {
	// make sure GALID matches

	ivar       = HOSTLIB.IVAR_GALID;
	igal_zsort = HOSTLIB.LIBINDEX_ZSORT[igal_unsort];
	GALID      = (long long)HOSTLIB.VALUE_ZSORTED[ivar][igal_zsort] ;

	get_PARSE_WORD(0, 1, NEXTWORD); // read GALID
	sscanf(NEXTWORD, "%lld", &GALID_orig);
	if ( GALID != GALID_orig ) {
	  sprintf(c1err,"GALID mis-match for igal_unsort=%d", igal_unsort);
	  sprintf(c2err,"GALID(orig)=%lld, but stored GALID=%lld",
		  GALID_orig, GALID ) ;
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	}

	for(ifilt=1; ifilt <= NFILT; ifilt++ ) {
	  sprintf(cval, " %6.3f", MAG_STORE[ifilt][igal_unsort] );
	  strcat(LINE_APPEND,cval);
	}
	igal_unsort++ ;
      }
    }

    ptrCR = strchr(LINE,'\n'); if(ptrCR){*ptrCR=' ';} // remove <CR>
    fprintf(FP_NEW,"%s %s\n", LINE, LINE_APPEND);
  }

  fclose(FP_ORIG);   fclose(FP_NEW);
  xxxxxxxxxxxx end mark xxxxxxxxxxxxx  */

  // ------------------------------------
  free(GENFLUX_LIST); free(GENMAG_LIST); free(MAG_STORE);

  exit(0);

  return ;

} // end rewrite_HOSTLIB_plusMags


// ======================================
double integmag_hostSpec(int IFILT_OBS, double z, int DUMPFLAG) {

  // Sep 2019
  // integrate global HOSTSPEC.FLAM_EVT over IFILT_OBS bandpass
  // and return synthetic mag for filter IFILT_OBS.

  int     NBLAM_BASIS   = HOSTSPEC.NBIN_WAVE ;
  double *LAMCEN_BASIS  = HOSTSPEC.WAVE_CEN ;
  double  LAMMIN_BASIS  = HOSTSPEC.WAVE_MIN[0] ;
  double  LAMMAX_BASIS  = HOSTSPEC.WAVE_MAX[NBLAM_BASIS-1] ;

  int    IFILT         = IFILTMAP_SEDMODEL[IFILT_OBS];
  int    NBLAM_FILT    = FILTER_SEDMODEL[IFILT].NLAM ;
  double LAMSTEP_FILT  = FILTER_SEDMODEL[IFILT].lamstep ; 
  double ZP_FILT       = FILTER_SEDMODEL[IFILT].ZP ;
  char  *cfilt         = FILTER_SEDMODEL[IFILT].name ;
  double hc8           = (double)hc;
  double z1            = 1.0 + z;

  double TRANS, LAMOBS, LAMREST, FLAM, FLUX, FSUM=0.0, mag=0.0 ;
  double SUMTRANS_TOT=0.0, SUMTRANS_UNDEFINED=0.0 ;
  int    ilamobs, NBLAM_UNDEFINED=0;
  char   comment[100];
  char   fnam[] = "integmag_hostSpec" ;

  // -------------- BEGIN -------------

  sprintf(comment,"%s(%s)", fnam, cfilt);

  for ( ilamobs=0; ilamobs < NBLAM_FILT; ilamobs++ ) {

    LAMOBS  = FILTER_SEDMODEL[IFILT].lam[ilamobs] ;
    LAMREST = LAMOBS/z1;
    TRANS   = FILTER_SEDMODEL[IFILT].transSN[ilamobs] ;
    SUMTRANS_TOT += TRANS;

    if ( LAMOBS >= LAMMIN_BASIS  &&  LAMOBS <= LAMMAX_BASIS ) {
      // interpolate to get FLAM
      FLAM = interp_1DFUN(OPT_INTERP_LINEAR, LAMREST, NBLAM_BASIS, LAMCEN_BASIS,
			  HOSTSPEC.FLAM_EVT, comment );   
    }
    else {
      FLAM  = 0.0 ;   
      NBLAM_UNDEFINED++ ;     
      SUMTRANS_UNDEFINED += TRANS;
    }

    FSUM += (FLAM * LAMSTEP_FILT * TRANS * LAMOBS) ;

  } // end ilamobs 


  // - - - - - - - - - - - - - - - -
  // apply normalization factors
  FSUM /= hc8 ;

  // - - - - - - - - - - - - - - 
  double FRAC_UNDEFINED = SUMTRANS_UNDEFINED/SUMTRANS_TOT ;
  int IWARN_FILTER_FRAC = (FRAC_UNDEFINED > 0.01) ;

  if ( FSUM > 1.0E-20 && !IWARN_FILTER_FRAC ) {
    mag = ZP_FILT - 2.5*log10(FSUM) ;
  }
  else {
    mag = MAG_UNDEFINED ;
  }

  if ( DUMPFLAG ) {
    printf(" xxx %s: FSUM(%d:%s)=%10.3le  ZP=%.3f  mag=%.3f\n", 
	   fnam, IFILT, cfilt, FSUM, ZP_FILT, mag );
    fflush(stdout);
  }

  
  if ( IWARN_FILTER_FRAC ) {
    HOSTSPEC.NWARN_INTEG_HOSTMAG[IFILT]++ ;
    int NWARN = HOSTSPEC.NWARN_INTEG_HOSTMAG[IFILT] ;
    if ( NWARN == 1 ) {
      printf(" %s WARNING filter %s: %.3f of trans undefined.\n",
	     fnam, cfilt, FRAC_UNDEFINED);
      fflush(stdout);
    }
  }  

  return(mag) ;

} // end integmag_hostSpec


// ===================================
void rewrite_HOSTLIB_plusNbr(void) {

  // Created Nov 2019 by R.Kessler
  // Re-write HOSTLIB with list of nearby IGALs.
  //
  // Beware that 'igal' is a redshift-sorted index for the other 
  // HOSTLIB functions, but here we use the original HOSTLIB order.

  int  NGAL        = HOSTLIB.NGAL_STORE;
  int  IVAR_RA     = HOSTLIB.IVAR_RA ;
  int  IVAR_DEC    = HOSTLIB.IVAR_DEC ;
  int  MEMD        = NGAL * sizeof(double);
  int  MEMI        = NGAL * sizeof(int);

  int   igal_unsort, igal_zsort, igal_DECsort, isort ;

  HOSTLIB_APPEND_DEF HOSTLIB_APPEND;
  char  *LINE_APPEND, MSG[100] ;

  // internal debug
  int  NGAL_DEBUG  = -500;

  char *INPUT_FILE = INPUTS.INPUT_FILE_LIST[0];
  char fnam[] = "rewrite_HOSTLIB_plusNbr" ;

  // --------------- BEGIN ---------------

  print_banner(fnam);

  if ( IVAR_RA < 0 || IVAR_DEC < 0 ) {
    sprintf(c1err,"Must include galaxy coords to find neighbors.");
    sprintf(c2err,"Check VARNAMES in HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  printf("\t NGAL[READ,STORE] = %d, %d \n", 
	 HOSTLIB.NGAL_READ, HOSTLIB.NGAL_STORE);

  malloc_HOSTLIB_APPEND(NGAL, &HOSTLIB_APPEND);

  LINE_APPEND = (char*) malloc (MXCHAR_LINE_HOSTLIB * sizeof(char) ) ;
  HOSTLIB_NBR.SKY_SORTED_DEC          = (double*) malloc(MEMD) ;
  HOSTLIB_NBR.SKY_SORTED_RA           = (double*) malloc(MEMD) ;
  HOSTLIB_NBR.SKY_SORTED_IGAL_zsort   = (int*) malloc(MEMI) ;
  HOSTLIB_NBR.SKY_SORTED_IGAL_DECsort = (int*) malloc(MEMI) ;
  HOSTLIB_NBR.GALID_atNNBR_MAX = -9 ;
  HOSTLIB_NBR.NNBR_MAX         =  0 ;

  // sort by DEC to improve NBR-matching speed
  int  ORDER_SORT = +1 ;
  double *ptrDEC = HOSTLIB.VALUE_ZSORTED[IVAR_DEC] ; 
  double *ptrRA  = HOSTLIB.VALUE_ZSORTED[IVAR_RA] ; 
  sortDouble( NGAL, ptrDEC, ORDER_SORT, HOSTLIB_NBR.SKY_SORTED_IGAL_DECsort);
  

  // load new lists of RA & DEC sorted by DEC
  for(igal_DECsort=0; igal_DECsort < NGAL; igal_DECsort++ ) {
    igal_zsort = HOSTLIB_NBR.SKY_SORTED_IGAL_DECsort[igal_DECsort];
    HOSTLIB_NBR.SKY_SORTED_DEC[igal_DECsort]      = ptrDEC[igal_zsort] ;
    HOSTLIB_NBR.SKY_SORTED_RA[igal_DECsort]       = ptrRA[igal_zsort] ;
    HOSTLIB_NBR.SKY_SORTED_IGAL_zsort[igal_zsort] = igal_DECsort ;
    if ( igal_DECsort < -5 || igal_DECsort > NGAL+5 ) {
      printf(" xxx igal_DECsort=%d  igal_zsort=%6d,  DEC = %9.5f \n",
	     igal_DECsort, igal_zsort, ptrDEC[igal_zsort] ); fflush(stdout);
    }
  }
  
  // ----------------------------
  if ( NGAL_DEBUG > 0 ) { NGAL = NGAL_DEBUG; }

  // init diagnistic counters (filled in get_LINE_APPEND_HOSTLIB_plusNbr)
  monitor_HOSTLIB_plusNbr(0,&HOSTLIB_APPEND); 

  // loop over all galaxies and prepare string to append.
  for(igal_unsort=0; igal_unsort < NGAL; igal_unsort++ ) {

    // search for neighbors and fill line to append
    get_LINE_APPEND_HOSTLIB_plusNbr(igal_unsort, LINE_APPEND);

    sprintf(HOSTLIB_APPEND.LINE_APPEND[igal_unsort],"%s", LINE_APPEND);

  } // end igal_unsort loop over all galaxies


  // - - - - - - - - - - - - 

  sprintf(HOSTLIB_APPEND.FILENAME_SUFFIX, "+HOSTNBR");
  sprintf(HOSTLIB_APPEND.VARNAMES_APPEND, "NBR_LIST" );


  // construct message strings for top of new HOSTLIB
  sprintf(MSG, "Append up to %d host neighbors within %.1f'' radius.",
	  HOSTLIB_NBR.NNBR_WRITE_MAX, HOSTLIB_NBR.SEPNBR_MAX );
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);
 
  sprintf(MSG, "snlc_sim.exe %s +HOSTNBR  SEPNBR_MAX %.1f  NNBR_WRITE_MAX %d",
	  INPUT_FILE, HOSTLIB_NBR.NNBR_WRITE_MAX, HOSTLIB_NBR.SEPNBR_MAX );
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);

  sprintf(MSG,"Added column NBR_LIST = "
	  "comma-sep list of row numbers (not GALID)" );
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);

  // write out monitor info
  monitor_HOSTLIB_plusNbr(1,&HOSTLIB_APPEND);

  // exectute re-write
  rewrite_HOSTLIB(&HOSTLIB_APPEND);

  exit(0);

  return ;

} // end rewrite_HOSTLIB_plusNbr


// ==============================
void get_LINE_APPEND_HOSTLIB_plusNbr(int igal_unsort, char *LINE_APPEND) {


#define MXNNBR_STORE 100         // max number of neighbors to track
  double SEPNBR_MAX      = HOSTLIB_NBR.SEPNBR_MAX ;
  int    NNBR_WRITE_MAX  = HOSTLIB_NBR.NNBR_WRITE_MAX ;
  //  int    MXCHAR_NBR_LIST = HOSTLIB_NBR.MXCHAR_NBR_LIST ;

  double ASEC_PER_DEG  = 3600.0 ;

  int  NGAL        = HOSTLIB.NGAL_STORE;
  int  IVAR_GALID  = HOSTLIB.IVAR_GALID;
  int  IVAR_RA     = HOSTLIB.IVAR_RA ;
  int  IVAR_DEC    = HOSTLIB.IVAR_DEC ;
  int  LDMP        = (igal_unsort < -3);

  double SEP_NBR_LIST[MXNNBR_STORE];
  int    IGAL_LIST[MXNNBR_STORE];
  long long GALID, GALID_NBR, GALID_LIST[MXNNBR_STORE] ;
  double RA_GAL, DEC_GAL, RA_NBR, DEC_NBR, SEP_NBR, SEP_DEC;
  int  igal_zsort, igal2_unsort, igal2_zsort, igal_DECsort ;
  int  NNBR, NTRY, j, ISORT_CHANGE, isort, NPASS_DEC, LSTDOUT ;
  char cval[20], cval2[20], LINE_STDOUT[200];;
  char fnam[] = "get_LINE_APPEND_HOSTLIB_plusNbr";

  // ------------ BEGIN -----------

  igal_zsort   = HOSTLIB.LIBINDEX_ZSORT[igal_unsort];
  igal_DECsort = HOSTLIB_NBR.SKY_SORTED_IGAL_zsort[igal_zsort];
  RA_GAL       = HOSTLIB.VALUE_ZSORTED[IVAR_RA][igal_zsort] ;  
  DEC_GAL      = HOSTLIB.VALUE_ZSORTED[IVAR_DEC][igal_zsort] ; 
  GALID        = (long long)HOSTLIB.VALUE_ZSORTED[IVAR_GALID][igal_zsort] ;
  

  if ( LDMP ) {
    printf("\n xxx ---------- %s DUMP ---------------- \n", fnam);
    printf(" xxx Input igal_unsort=%d  \n", igal_unsort);
    printf(" xxx recover igal_zsort=%d,  igal_DECsort=%d \n", 
	   igal_zsort, igal_DECsort );
    fflush(stdout);
  }


  sprintf(LINE_APPEND,"-1");
  LINE_STDOUT[0] = 0 ;
  NNBR = NTRY = 0 ; NPASS_DEC = 1;  ISORT_CHANGE=1;

  while ( NPASS_DEC > 0 ) {
    NPASS_DEC = 0;

    for(j = -1; j <=1; j+=2 ) { // try both directions

      NTRY++ ;
      isort = igal_DECsort + j*ISORT_CHANGE;

      if ( isort < 0  || isort >= NGAL ) { continue; }
      igal2_zsort   = HOSTLIB_NBR.SKY_SORTED_IGAL_DECsort[isort];
      igal2_unsort  = HOSTLIB.LIBINDEX_UNSORT[igal2_zsort];

      RA_NBR      = HOSTLIB_NBR.SKY_SORTED_RA[isort] ;  
      DEC_NBR     = HOSTLIB_NBR.SKY_SORTED_DEC[isort] ;  
      SEP_DEC     = fabs(DEC_NBR - DEC_GAL)*ASEC_PER_DEG;

      /*      
      if ( SEP_DEC < 1.0E8 ) {
	printf("\t zzz igal2_zsort=%d  DEC_GAL/NBR=%f/%f  SEP=%.2f\n", 
	       igal_zsort, DEC_GAL, DEC_NBR, SEP_DEC );
      }
      */

      if ( SEP_DEC > SEPNBR_MAX ) { continue ; }
      NPASS_DEC++ ;

      /*
      if ( LDMP ) {
	printf("\t j=%2d, isort=%6d  igal2=%6d NPASS=%d  dDEC=%f\n",
	       j, isort, igal2_unsort, NPASS_DEC, SEP_DEC); fflush(stdout);
      }
      */

      SEP_NBR = angSep(RA_GAL, DEC_GAL, RA_NBR, DEC_NBR, ASEC_PER_DEG);
      if ( SEP_NBR > SEPNBR_MAX ) { continue ; }
      GALID_NBR = (long long)HOSTLIB.VALUE_ZSORTED[IVAR_GALID][igal2_zsort] ;
      
      SEP_NBR_LIST[NNBR] = SEP_NBR;
      IGAL_LIST[NNBR]    = igal2_unsort;
      GALID_LIST[NNBR]   = GALID_NBR;
      NNBR++ ;
      if ( NNBR > HOSTLIB_NBR.NNBR_MAX ) { 
	HOSTLIB_NBR.NNBR_MAX = NNBR; 
	HOSTLIB_NBR.GALID_atNNBR_MAX = GALID;
      }

    } // end j loop over smaller/larger DEC

    ISORT_CHANGE++ ;
  } // end while



  /* xxxxx mark delete xxxx
  // - - - - - -
  for(igal2_unsort=0; igal2_unsort < NGAL; igal2_unsort++ ) {
    if ( igal2_unsort == igal_unsort ) { continue; }
    if ( strlen(LINE_APPEND) > MXCHAR_LINE_HOSTLIB - 20 ) { continue; }

    igal2_zsort = HOSTLIB.LIBINDEX_ZSORT[igal2_unsort];
    RA_NBR      = HOSTLIB.VALUE_ZSORTED[IVAR_RA][igal2_zsort] ;  
    DEC_NBR     = HOSTLIB.VALUE_ZSORTED[IVAR_DEC][igal2_zsort] ; 
    if ( fabs(DEC_NBR-DEC_GAL) > SEPMAX_NBR ) { continue ; }
    
    SEP_NBR = angSep(RA_GAL, DEC_GAL, RA_NBR, DEC_NBR,  (double)3600.);
    if ( SEP_NBR > SEPMAX_NBR ) { continue ; }
    GALID_NBR = (long long)HOSTLIB.VALUE_ZSORTED[IVAR_GALID][igal2_zsort] ;
    
    SEP_NBR_LIST[NNBR] = SEP_NBR;
    IGAL_LIST[NNBR]    = igal2_unsort;
    NNBR++ ;

    if ( NNBR > MXNNBR ) {
      sprintf(c1err,"NNBR=%d exceeds bound of MXNNBR", NNBR);
      sprintf(c2err,"Check SEPMAX_NBR = %.2f arcSec and HOSTLIB density.",
	      SEPMAX_NBR );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }
  xxxx */

  // - - - - - - - - - - - - - - - - - - 
  // sort SEP_NBR_LIST inascending order
  int ORDER_SORT = +1 ;    // increasing order
  int UNSORT[MXNNBR_STORE];
  int i, IGAL, TRUNCATE=0;


  if ( LDMP ) { 
    printf("\n xxx %s DUMP igal_unsort=%d  NNBR=%d   NTRY=%d\n", 
	   fnam, igal_unsort, NNBR, NTRY );
  }

  sortDouble( NNBR, SEP_NBR_LIST, ORDER_SORT, UNSORT ) ;
  for(i=0; i < NNBR; i++ ) {
    isort       = UNSORT[i];
    SEP_NBR     = SEP_NBR_LIST[isort];
    IGAL        = IGAL_LIST[isort];
    GALID_NBR   = GALID_LIST[isort];

    if ( i > NNBR_WRITE_MAX ) 
      { TRUNCATE = 1 ; continue; }

    if ( strlen(LINE_APPEND) > MXCHAR_NBR_LIST) 
      { TRUNCATE = 1 ; continue; }

    if ( LDMP ) {
      printf("\t xxx SEP(%2d) = %8.2f  IGAL=%d \n", isort, SEP_NBR, IGAL);
      fflush(stdout);
    }

    if ( i == 0 )  { 
      sprintf(cval, "%d",   IGAL     ); LINE_APPEND[0]=0; 
      sprintf(cval2,"%lld", GALID_NBR); LINE_STDOUT[0]=0; 
    }
    else   { 
      sprintf(cval, ",%d",   IGAL); 
      sprintf(cval2,",%lld", GALID_NBR); 
    }

    // maybe need flag for missing NBR ???
    strcat(LINE_APPEND,cval); 
    strcat(LINE_STDOUT,cval2); 

  } // end i loop over NBR


  // set stdout dump for a few events so that igal_unsiort <-> GALID
  // can be visually checked when appended HOSTLIB is read back.
  LSTDOUT = (igal_unsort < 10 || igal_unsort > NGAL-10);
  if ( LSTDOUT  && strlen(LINE_STDOUT) > 0 ) {
    printf("# ---------------------------------------------------- \n");
    printf(" Crosscheck dump for igal_read=%d, GALID=%lld \n",
	   igal_unsort, GALID); fflush(stdout);
    
    printf("\t READ  NBR list: %s\n", LINE_APPEND );
    printf("\t GALID NBR list: %s\n", LINE_STDOUT );
    fflush(stdout);
  }

  if ( (igal_unsort % 10000) == 0 ) {
    NNBR = HOSTLIB_NBR.NNBR_MAX ; GALID=HOSTLIB_NBR.GALID_atNNBR_MAX;
    printf("\t Processing igal %8d of %8d  (NNBR_MAX=%2d for GALID=%lld)\n", 
	   igal_unsort, NGAL, NNBR, GALID );
    fflush(stdout);
  }

  if ( NNBR < 100 ) { HOSTLIB_NBR.NGAL_PER_NNBR[NNBR]++ ; }
  if ( TRUNCATE   ) { HOSTLIB_NBR.NGAL_TRUNCATE++ ; }

  return ;

} // end get_LINE_APPEND_HOSTLIB_plusNbr

// =========================================
void  monitor_HOSTLIB_plusNbr(int OPT, HOSTLIB_APPEND_DEF *HOSTLIB_APPEND) {

  // OPT=0 --> init
  // OPT=1 --> write info to stdout and to addComment

  int NGAL        = HOSTLIB.NGAL_STORE;
  int nnbr, NGAL_TMP;
  float frac;
  char MSG[100];
  // ------------ BEGIN -------------

  if ( OPT == 0 ) {
    for(nnbr=0; nnbr < 100; nnbr++ ) { HOSTLIB_NBR.NGAL_PER_NNBR[nnbr]=0; }
    HOSTLIB_NBR.NGAL_TRUNCATE = 0 ;
  }
  else {
    printf("\n");
    for(nnbr=0; nnbr <= HOSTLIB_NBR.NNBR_WRITE_MAX; nnbr++ ) {
      NGAL_TMP = HOSTLIB_NBR.NGAL_PER_NNBR[nnbr];
      frac     = (float)NGAL_TMP / (float)NGAL;
      sprintf(MSG, "\t HOSTLIB fraction with %2d NBR: %8.3f %% ",
	     nnbr, 100.0*frac ); 
      printf("%s\n", MSG); fflush(stdout);
      addComment_HOSTLIB_APPEND(MSG,HOSTLIB_APPEND);
    }

    frac = (float)HOSTLIB_NBR.NGAL_TRUNCATE / (float)NGAL ;
    sprintf(MSG,"\t Truncated fraction with > %d NBR: %8.3f %% ",
	   HOSTLIB_NBR.NNBR_WRITE_MAX, 100.*frac); 
    printf("%s\n", MSG); fflush(stdout);
    addComment_HOSTLIB_APPEND(MSG,HOSTLIB_APPEND);
  }

  return ;

} // end monitor_HOSTLIB_plusNbr
