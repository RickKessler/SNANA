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

  LOGMASS is internally converted to LOGMASS_TRUE.
  User can provide LOGMASS_TRUE & LOGMASS_OBS, or provide
  LOGMASS_TRUE and LOGMASS_ERR from which LOGMASS_OBS is
  computed from smearing the true LOGMASS.

 TODO_LIST: 
  - sum GALMAG[iPSF] over galaxy neihbors (Jan 31 2020)
  - fix ZPHOT so that it is computed for each neighbor, rather
     than only the true host (Jan 31 2020)


          HISTORY
  ------------------------


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

 Nov 18 2019:
   + fix problems generating SN position near galaxy.
     1) fix bug computing SN coords from ANGLE and reduced_R
     2) pick random ANGLE weighted by DLR^2 (no longer random over 0:2PI)
        (see new function GEN_SNHOST_ANGLE)
     WARNING: still need to visually verify 2D profile

 Jan 14 2020: 
   + in GEN_SNHOST_DDLR(), a & b are wgted avg among Sersic terms.

 Jan 31 2020: little more work on DDLR_SORT:
    +  make sure sorted (total) host mags have extinction
         (but beware that GALMAG under SN is only from TRUE host)
    +  implement LOGMASS_TRUE, LOGMASS_OBS

 Feb 26 2020: check SWAPZPHOT option to set ZTRUE = ZPHOT
 Feb 27 2020: fix index bug translating NBR_LIST.
 May 24 2020: add option to VPEC, VPEC_ERR from HOSTLIB
 Jul 14 2020: automatically store WGTMAP variables;
              see new function read_VARNAMES_WGTMAP

 Sep 04 2020 : implement REQUIRE_DOCANA 

 Jan 22 2021 : 
   + convert GENRANGE_REDSHIFT to zHEL min/max, then compare to HOSTLIB range.
   + count NSTAR for z < ZMAX_STAR; give warning if NSTAR>0

 Apr 30 2021: abort on NaN; abort if synthetic band is in HOSTLIB.

 Jul 6 2021: require DOCANA in HOSTLIB and WGTMAP

 Aug 11 2021: for SIMLIB model, if HOSTLIB doesn't match that used to
              generate fakes, fix to work.

 Sep 16 2021: abort on duplicate columns in HOSTLIB; e.g., 
              ZERR and ZPHOTERR are both mapped to ZPHOT_ERR -> abort.

 Nov 18 2021: allow WGTMAP to be missing SNMAGSHIFT column

 Dec 30 2021: refactor GEN_SNHOST_GALID() to use binary search for speed.

=========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_cosmology.h"
#include "sntools_trigger.h"
#include "snlc_sim.h" 
#include "sntools_host.h"
#include "sntools_output.h"
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

  // M. Vincenzi Feb 2022
  REFAC_HOSTLIB = false; // ( INPUTS.DEBUG_FLAG == 203 );

  // ---------------- BEGIN -------------

  USE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USE) ;
  if ( USE == 0 ) 
    {  checkAbort_noHOSTLIB() ;  return ;  }
  else
    { checkAbort_HOSTLIB(); }

  if ( INPUTS.README_DUMPFLAG ) { return; } // Feb 21 2022

  TIME_INIT_HOSTLIB[0]  = time(NULL);

  sprintf(BANNER,"%s: Read host-galaxy library (MSKOPT=%d)",
	  fnam, INPUTS.HOSTLIB_MSKOPT);
  print_banner(BANNER);

  OPTMASK_OPENFILE_HOSTLIB = 0;
  if ( INPUTS.REQUIRE_DOCANA ) 
    { OPTMASK_OPENFILE_HOSTLIB = OPENMASK_REQUIRE_DOCANA; }

  // check for spectral templates to determin host spectrum
  read_specTable_HOSTLIB();

  // set inital values for HOSTLIB structure
  initvar_HOSTLIB();

  // check to read external WEIGHT-MAP instead of the HOSTLIB WEIGHT-MAP
  read_HOSTLIB_WGTMAP();

  // open hostlib and start reading
  open_HOSTLIB(&fp_hostlib);     // open and return file pointer
  read_head_HOSTLIB(fp_hostlib); // read header info: VARNAMES, ...
  close_HOSTLIB(fp_hostlib);     // close HOSTLIB to rewind

  // check for match among spec templates and hostlib varnames (Jun 2019)
  match_specTable_HOSTVAR();

  // re-open hostlib for GAL keys so that it works for gzipped files.
  // (cannot rewind gzip file, so close and re-open is only way)
  open_HOSTLIB(&fp_hostlib);     // re-open
  read_gal_HOSTLIB(fp_hostlib);  // read "GAL:" keys
  close_HOSTLIB(fp_hostlib);     // close HOSTLIB

  // summarize SNPARams that were/weren't found
  summary_snpar_HOSTLIB();

  // sort HOSTLIB entries by redshift
  sortz_HOSTLIB();

  // abort if any GALID+ZTRUE pair appears more than once
  check_duplicate_GALID();
  
  // set redshift pointers for faster lookup
  zptr_HOSTLIB();

  // setup optional wgt-map grid
  int IGAL_START = 0,  IGAL_END = HOSTLIB.NGAL_STORE-1;
  init_HOSTLIB_WGTMAP(1, IGAL_START, IGAL_END);
  
  
  // read optional EFF(zPHOT) vs. ZTRUE (Aug 2015)
  init_HOSTLIB_ZPHOTEFF();

  // check for zphot quantiles (Apr 2022)
  init_HOSTLIB_ZPHOT_QUANTILE();

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

  return ;

} // end of INIT_HOSTLIB


// ====================================
void initvar_HOSTLIB(void) {

  // one-time init of variables used for HOSTLIB

  int ivar, j, igal, ifilt  ;
  char fnam[] = "initvar_HOSTLIB" ;

  // ----------- BEGIN -------------

  sprintf(PATH_DEFAULT_HOSTLIB, "%s %s/simlib", 
	  PATH_USER_INPUT, PATH_SNDATA_ROOT ); // Jul 14 2020
 

  NCALL_GEN_SNHOST_DRIVER = 0 ;

  HOSTLIB.SORTFLAG     = 0 ;
  HOSTLIB.NGAL_READ    = 0 ;
  HOSTLIB.NGAL_STORE   = 0 ;
  HOSTLIB.NVAR_ALL     = 0 ; 
  HOSTLIB.NVAR_STORE    = 0 ; 
  HOSTLIB.MALLOCSIZE_D  = 0 ;
  HOSTLIB.MALLOCSIZE_I  = 0 ;
  HOSTLIB.MALLOCSIZE_Cp = 0 ;

  HOSTLIB.NVAR_SNPAR   = 0 ;
  HOSTLIB.VARSTRING_SNPAR[0] = 0 ;

  HOSTLIB.NERR_NAN = 0 ;
  HOSTLIB.NSTAR = 0;

  for ( ivar=0; ivar < MXVAR_HOSTLIB; ivar++ )  { 
    HOSTLIB.IVAR_WGTMAP[ivar] = -9 ; 
    HOSTLIB.IVAR_STORE[ivar]  = -9 ; 
    HOSTLIB.VALMAX[ivar]      = -1.0E18 ;
    HOSTLIB.VALMIN[ivar]      = +1.0E18 ;
    sprintf(HOSTLIB.VARNAME_STORE[ivar],"%s", NULLSTRING );

    HOSTLIB.IS_SNPAR_OPTIONAL[ivar] = 0 ;
    HOSTLIB.IS_SNPAR_STORE[ivar]    = 0 ;
  }
  
  
  malloc_HOSTGAL_PROPERTY();

  HOSTLIB.ZGAPMAX       = -9. ;
  HOSTLIB.Z_ATGAPMAX[0] = -9. ;
  HOSTLIB.Z_ATGAPMAX[1] = -9. ;
  HOSTLIB.ZGAPAVG       =   0.0 ;
  HOSTLIB.ZMAX          = -99.0 ; 
  HOSTLIB.ZMIN          = +99.0 ; 

  HOSTLIB.NLINE_COMMENT = 0;

  SERSIC_PROFILE.NFIX    = 0;
  SERSIC_PROFILE.NPROF   = 0 ; 

  SERSIC_TABLE.NBIN_reduced   = 0 ;
  SERSIC_TABLE.TABLEMEMORY    = 0 ;

  HOSTLIB.NFILT_MAGOBS = 0;
  HOSTLIB.NZPHOT_Q = 0;
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

  HOSTLIB.IVAR_a_DLR  = -9 ;
  HOSTLIB.IVAR_b_DLR  = -9 ;
  // ------------------------
  // set array of required and optional keys
  init_REQUIRED_HOSTVAR();

  // now the optional keys
  init_OPTIONAL_HOSTVAR();

  HOSTLIB_WGTMAP.GRIDMAP.NDIM = 0;
  HOSTLIB_WGTMAP.GRIDMAP.NFUN = 0;
  HOSTLIB_WGTMAP.GRIDMAP.NROW = 0;
  HOSTLIB_WGTMAP.WGTMAX    = 0.0 ;
  HOSTLIB_WGTMAP.READSTAT  = false ; // init to weight-map NOT read
  HOSTLIB_WGTMAP.NCHECKLIST = 0;
  HOSTLIB_WGTMAP.N_SNVAR      =  0 ;
  HOSTLIB_WGTMAP.NBTOT_SNVAR  =  1 ; // at least 1 dummy bin of no WGTMAP
  HOSTLIB_WGTMAP.ibin_SNVAR   = -9 ; 
  HOSTLIB_WGTMAP.OPT_EXTRAP   =  0 ;
  for ( ivar=0; ivar < MXVAR_WGTMAP_HOSTLIB; ivar++ ) {  
    sprintf(HOSTLIB_WGTMAP.VARNAME[ivar], "%s", NULLSTRING );
    HOSTLIB_WGTMAP.NB1D_SNVAR[ivar] = 0 ;
    HOSTLIB_WGTMAP.IS_SNVAR[ivar]   = false ;
    HOSTLIB_WGTMAP.ISPARSE_SNVAR[ivar]   = -9 ;
    HOSTLIB_WGTMAP.INVSPARSE_SNVAR[ivar] = -9 ;
  }

  for ( igal = 0; igal < MXCHECK_WGTMAP ; igal++ ) 
    { HOSTLIB_WGTMAP.CHECKLIST_IGAL[igal] = -9 ;  }
   
  // malloc temp string pointers for splitString function
  for(ivar=0; ivar < MXTMPWORD_HOSTLIB; ivar++ ) 
    { TMPWORD_HOSTLIB[ivar] = (char*)malloc( 100*sizeof(char) ); }

  reset_SNHOSTGAL_DDLR_SORT(MXNBR_LIST);

  
  HOSTLIB.IGAL_FORCE = -9 ;

  return ;

} // end of initvar_HOSTLIB

// ===================================================
void malloc_HOSTGAL_PROPERTY(void) {

  // Created Feb 2022 by M.Vincenzi and R.Kessler
  
  int N_PROP;
  int MEM, nbr, i, index;
  char *varName, *BASENAME;
  char fnam[] = "malloc_HOSTGAL_PROPERTY";
  // ------------ BEGIN ------------

  N_PROP = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING, HOSTGAL_PROPERTY_NAME_LIST);

  // set a global variable 
  N_HOSTGAL_PROPERTY = N_PROP;

  MEM = N_PROP * sizeof(HOSTGAL_PROPERTY_IVAR_DEF);
  HOSTLIB.HOSTGAL_PROPERTY_IVAR = (HOSTGAL_PROPERTY_IVAR_DEF*) malloc(MEM);
  
  MEM = N_PROP * sizeof(HOSTGAL_PROPERTY_VALUE_DEF);
  for (nbr=0; nbr<MXNBR_LIST; nbr++) {
    SNHOSTGAL_DDLR_SORT[nbr].HOSTGAL_PROPERTY_VALUE = 
      (HOSTGAL_PROPERTY_VALUE_DEF*) malloc(MEM);

    for (i=0; i < N_PROP; i++) {
      SNHOSTGAL_DDLR_SORT[nbr].HOSTGAL_PROPERTY_VALUE[i].VAL_TRUE =
	  HOSTLIB_PROPERTY_UNDEFINED ;
      SNHOSTGAL_DDLR_SORT[nbr].HOSTGAL_PROPERTY_VALUE[i].VAL_OBS =
	  HOSTLIB_PROPERTY_UNDEFINED ;
      SNHOSTGAL_DDLR_SORT[nbr].HOSTGAL_PROPERTY_VALUE[i].VAL_ERR =
	  HOSTLIB_PROPERTY_UNDEFINED ;
    }
  }

  for (i=0; i<N_PROP; i++){
    BASENAME = HOSTLIB.HOSTGAL_PROPERTY_IVAR[i].BASENAME;
    get_PARSE_WORD(0,i,BASENAME);
    HOSTLIB.HOSTGAL_PROPERTY_IVAR[i].SCALE_ERR = 1.0; //default
  } 

  // check legacy input that works only for LOGMASS
  if (INPUTS.HOSTLIB_SCALE_LOGMASS_ERR!= 1.0){
    i = getindex_HOSTGAL_PROPERTY(HOSTGAL_PROPERTY_BASENAME_LOGMASS);
    HOSTLIB.HOSTGAL_PROPERTY_IVAR[i].SCALE_ERR = 
      INPUTS.HOSTLIB_SCALE_LOGMASS_ERR;
  }

  // check for user input error scale for arbitrary property
  char *str_error = INPUTS.HOSTLIB_SCALE_PROPERTY_ERR;

  if ( !IGNOREFILE(str_error) ) {
    char **tmp_item_list, basename[40], item[40];
    int N_item;
    double scale;
    if ( INPUTS.HOSTLIB_SCALE_LOGMASS_ERR != 1.0) {
      sprintf(c1err, "Cannot specify both "
	      "HOSTLIB_SCALE_LOGMASS_ERR and HOSTLIB_SCALE_PROPERTY_ERR" ) ;
      sprintf(c2err, "Pick one or the other") ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
    parse_commaSepList("Host property error", str_error, 100, 40,
		       &N_item, &tmp_item_list ); // this is returned 

    for (i=0; i<N_item; i++){
      sprintf(item, "%s",tmp_item_list[i]); // presevre original tmp_item_list 
      extractStringOpt(item, basename);  // beaware that item is altered
      sscanf(item, "%le", &scale ) ; //convert string item to double scale
      index = getindex_HOSTGAL_PROPERTY(basename);
      if (index<0){
	sprintf(c1err, "INVALID basename = %s", basename ) ;
	sprintf(c2err, "Check key=HOSTLIB_SCALE_PROPERTY_ERR, "
		"item=%s", tmp_item_list[i] ) ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }
      HOSTLIB.HOSTGAL_PROPERTY_IVAR[index].SCALE_ERR = scale;
      printf("\t Scale %s error by %.2f\n", basename, scale);
    } 
  } // end of if SCALE_PROPERTY_ERR given

  return;
 
} // end of malloc_HOSTGAL_PROPERTY

// ==================================================
int getindex_HOSTGAL_PROPERTY(char *PROPERTY){
  // created Febr 2022 by M. Vincenzi 
  // for input hostgal property return index to HOSTGAL_PROPERTY_IVAR
  // return -9 is there's no match to input property
  //  
  char fnam[] = "getindex_HOSTGAL_PROPERTY";
  int index=-9;
  int N_PROP=N_HOSTGAL_PROPERTY;
  char *BASENAME;
  int i;

  for (i=0; i<N_PROP; i++){
    BASENAME = HOSTLIB.HOSTGAL_PROPERTY_IVAR[i].BASENAME;
    if (strcmp(BASENAME, PROPERTY)==0){ index=i; };
  }
  
  return index;
}  // end of getindex_HOSTGAL_PROPERTY

// ==========================================
void init_OPTIONAL_HOSTVAR(void) {

  // Feb 2014 [extracted from initvar_HOSTLOB]
  // Define array of optional variables to read from hostlib.
  // --> will read them if they exist; otherwise ignore them
  //     and continue.
  //
  // Jan 30 2015: allow RA or RA_HOST, DEC or DEC_HOST
  // May 23 2020: add VPEC and VPEC_ERR

  int NVAR, j, ifilt, ifilt_obs ;

  char anam[12], bnam[12], wnam[12], nnam[12];
  char varName[40], *cptr ;
  char fnam[] = "init_OPTIONAL_HOSTVAR" ;

  // ----------- BEGIN ---------------

  NVAR = 0;

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_TRUE_MATCH );

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_ZPHOT );
  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_ZPHOT_ERR );

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_VPEC );
  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_VPEC_ERR );

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
  sprintf(cptr,"%s", HOSTLIB_VARNAME_NBR_LIST ); 

  for (j=0; j<N_HOSTGAL_PROPERTY; j++){
    init_OPTIONAL_HOSTVAR_PROPERTY(HOSTLIB.HOSTGAL_PROPERTY_IVAR[j].BASENAME, &NVAR);
    }

  /* xxx Mark delete Febr 2022 after implementing more generic Host properties
     else {  //legacy
    NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ; 
    sprintf(cptr,"%s", HOSTLIB_VARNAME_LOGMASS_TRUE );
    NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
    sprintf(cptr,"%s", HOSTLIB_VARNAME_LOGMASS_ERR );
    NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
    sprintf(cptr,"%s", HOSTLIB_VARNAME_LOGMASS_OBS );
  }
  xxx */

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_GALID2 );

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_ELLIPTICITY );

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", HOSTLIB_VARNAME_SQRADIUS );

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
  for (j=0; j < MXBIN_ZPHOT_Q; j++){
    cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;   NVAR++;
    sprintf(cptr, "%s%d",    HOSTLIB_PREFIX_ZPHOT_Q, j); 

    // allow zero-padding ; e.g. ZPHOT_Q020
    if ( j < 100 ) {
      cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;   NVAR++;
      sprintf(cptr, "%s%3.3d", HOSTLIB_PREFIX_ZPHOT_Q, j); 
    }
  }

  cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;   NVAR++;
  sprintf(cptr, "%s", HOSTLIB_VARNAME_A_DLR);    // index 
  cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;   NVAR++;
  sprintf(cptr, "%s", HOSTLIB_VARNAME_B_DLR);    // index 

  cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ; NVAR++; 
  sprintf(cptr,"%s", HOSTLIB_VARNAME_ANGLE );

  char varName_err[50]; 
  // check for observer-frame mags '[filt]_obs' 
  // used to determine galaxy noise  to SN signal-flux
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    magkey_HOSTLIB(ifilt_obs,varName,varName_err); // returns varname
    cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ; NVAR++; 
    sprintf(cptr,"%s", varName);
    cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ; NVAR++;
    sprintf(cptr,"%s", varName_err);
  }


  // check for optional use of SN params from HOSTLIB; e.g., c, x1, delta
  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) > 0 ) {

    int N_SNPAR = 0 ;
    char tmpName[10][40] ;

    sprintf(tmpName[N_SNPAR],"%s", GENLC.SHAPEPAR_NAME);    N_SNPAR++ ; 
    sprintf(tmpName[N_SNPAR],"%s", GENLC.COLORPAR_NAME);    N_SNPAR++ ; 
    sprintf(tmpName[N_SNPAR],"%s", GENLC.SHAPEPAR2_NAME);   N_SNPAR++ ; 
    sprintf(tmpName[N_SNPAR],"%s", GENLC.COLORPAR2_NAME);   N_SNPAR++ ; 
    sprintf(tmpName[N_SNPAR],"%s", HOSTLIB_VARNAME_SNMAGSHIFT );  N_SNPAR++ ;

    // May 2020: for SALT2 model, allow additional host dust params  
    if ( INDEX_GENMODEL  == MODEL_SALT2 ) {
      sprintf(tmpName[N_SNPAR],"RV");   N_SNPAR++ ;
      sprintf(tmpName[N_SNPAR],"AV");   N_SNPAR++ ;
      sprintf(tmpName[N_SNPAR],"EBV");  N_SNPAR++ ;
    }   

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
  return ;
} // end of init_OPTIONAL_HOSTVAR

void init_OPTIONAL_HOSTVAR_PROPERTY(char *basename, int *NVAR_PROPERTY) {
#define N_SUFFIX_PROP 3 
  char fnam[]="init_OPTIONAL_HOSTVAR_PROPERTY";
  int NVAR = *NVAR_PROPERTY;
  char suffix_list[N_SUFFIX_PROP][20] = {"TRUE", "OBS", "ERR"};
  char *cptr ;

  NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
  sprintf(cptr,"%s", basename);

  int i;
  for (i=0; i<N_SUFFIX_PROP; i++){
    NVAR++; cptr = HOSTLIB.VARNAME_OPTIONAL[NVAR] ;
    char *suffix=suffix_list[i];
    sprintf(cptr,"%s_%s", basename, suffix);
  }

  *NVAR_PROPERTY = NVAR;
  return ;
} //end of init_OPTIONAL_HOSTVAR_PROPERTY


// ==========================================
void init_REQUIRED_HOSTVAR(void) {

  // Feb 23 2021: check for SPECDATA

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

  // check for IDSPECDATA (Feb 2021)
  if ( HOSTSPEC.ITABLE == ITABLE_SPECDATA ) {
    cptr = HOSTLIB.VARNAME_REQUIRED[NVAR] ;  NVAR++;
    sprintf(cptr, "%s", VARNAME_SPECDATA_HOSTLIB ); // IDSPECDATA
    LOAD = load_VARNAME_STORE(cptr) ;
  }

  HOSTLIB.NVAR_REQUIRED = NVAR ;

  // append STOREPAR automatically so that used HOSTLIB variables
  // can be added to SIMGEN-DUMP
  append_HOSTLIB_STOREPAR();

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

  if (INPUTS.HOSTLIB_MNINTFLUX_SNPOS < 0) {
    sprintf ( c1err, "Cannot have negative value for");
    sprintf ( c2err, "HOSTLIB_MNINTFLUX_SNPOS");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if (INPUTS.HOSTLIB_MNINTFLUX_SNPOS > INPUTS.HOSTLIB_MXINTFLUX_SNPOS) {
    sprintf ( c1err, "Cannot have HOSTLIB min sep > max sep");
    sprintf ( c2err, "MNINTFLUX_SNPOS > MXINTFLUX_SNPOS");
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
void append_HOSTLIB_STOREPAR(void) {

  // May 24 2020
  // Append INPUTS.HOSTLIB_STOREPAR_LIST with HOSTLIB variables
  // read to implement HOSTLIB_MSKOPT options. This allows
  // adding clearly-used variables to SIMGEN_DUMP without also 
  // defining HOSTLIB_STOREPAR.
  //
  // Jun 12 2020: set NVAR_zHOST after reading zHOST file

  char *STOREPAR  = INPUTS.HOSTLIB_STOREPAR_LIST ;
  int  ivar, NVAR_zHOST ;
  char *ptrVarName;
  FILE *fp ;
  char fnam[] = "append_HOSTLIB_STOREPAR" ;

  // -------------- BEGIN ------------

  // If zHOST_FILE exists,  copy variables from zHOST efficiency map to 
  // INPUTS.HOSTLIB_STOREPAR_LIST --> ensure that all of the 
  // HOSTLIB-zHOST parameters are read from the HOSTLIB.

  fp = open_zHOST_FILE(-1);
  if ( fp != NULL ) { 
    read_VARNAMES_zHOST(fp); fclose(fp);
    NVAR_zHOST = SEARCHEFF_zHOST[0].NVAR ; // Jun 12 2020
    for(ivar=0; ivar < NVAR_zHOST; ivar++ ) {
      ptrVarName = SEARCHEFF_zHOST[0].VARNAMES_HOSTLIB[ivar];
      catVarList_with_comma(STOREPAR,ptrVarName);
    } // end ivar
  }

  // - - - - - - - 
  // check PDF maps for populations

  fp = fopen(INPUTS.GENPDF_FILE,"rt");
  if ( fp ) {
    int MXVAR = 50, NVAR_SKIP=-1, NVAR, NKEY, *UNIQUE;
    char **VARNAMES;

    UNIQUE   = (int*)malloc(MXVAR*sizeof(int));
    VARNAMES = (char**) malloc( MXVAR*sizeof(char*) );
    for(ivar=0; ivar < MXVAR; ivar++ ) 
      { VARNAMES[ivar] = (char*) malloc( 40*sizeof(char) ); }

    read_VARNAMES_KEYS(fp, MXVAR, NVAR_SKIP, fnam, &NVAR, &NKEY,
		       UNIQUE, VARNAMES ) ;
    for(ivar=0; ivar < NVAR; ivar++ ) {
      ptrVarName = VARNAMES[ivar] ;
      if ( strcmp(ptrVarName,"PROB") == 0 ) { continue; }
      if ( UNIQUE[ivar] ) 
	{ catVarList_with_comma(STOREPAR,ptrVarName); }

      
     /* printf(" xxx %s: found varName[%2d] = '%s' (UNIQUE=%d)\n",
	     fnam, ivar, VARNAMES[ivar], UNIQUE[ivar]); fflush(stdout);
     */

    }
    fclose(fp);
  }

  // Jul 14 2020: check for WGTMAP variables here.
  //  WGTMAP is read later, but possibly after reading HOSTLIB
  //  VARNAMES where it is too late to store WGTMAP var.
  //  Here we are well before reading HOSTLIB.
  char VARLIST_WGTMAP[200];
  int NVAR_WGTMAP = 0 ;
  NVAR_WGTMAP = read_VARNAMES_WGTMAP(VARLIST_WGTMAP);
  if ( NVAR_WGTMAP > 0 ) {
    if ( strlen(STOREPAR) > 0 ) { strcat(STOREPAR,COMMA); }
    strcat(STOREPAR,VARLIST_WGTMAP);
  }
    
  return ;

} // end append_HOSTLIB_STOREPAR

// ====================================================
void  init_OUTVAR_HOSTLIB(void) {

  // Feb 2014
  // strip variables from comma-separated INPUTS.HOSTLIB_STOREPAR_LIST,
  // and store them in VARNAME_STORE and REQUIRED lists. 
  // Be careful NOT to store the same variable twice in case
  // user requests a variable that is already required or
  // in the WGTMAP.
  //
  // Feb 01 2020: 
  //  + call checkAlternateVarnames so that LOGMASS -> LOGMASS_TRUE.
  // Jun 02 2021: abort if strlen(HOSTLIB_STOREPAR_LIST) is too long

  int   NVAR_STOREPAR, NVAR_OUT, NVAR_REQ, LOAD, ivar, ivar2, ISDUPL ;
  char  VARLIST_ALL[MXPATHLEN], VARLIST_LOAD[MXPATHLEN];
  char  varName[60], *varName2;
  int   LENLIST, LDMP = 0 ;
  char  fnam[] = "init_OUTVAR_HOSTLIB"     ;

  // ---------------- BEGIN ----------------

  HOSTLIB_OUTVAR_EXTRA.NOUT = NVAR_OUT = 0 ;

  // check that STOREPAR_LIST doesn't overwrite array bound.
  // HOSTLIB_STOREPAR_LIST is allocated with 2*MXPATHLEN,  but here we
  // abort if len > MXPATHLEN ... the extra x2 memory padding avoids
  // unrelated memory overwrite errors so that we see the correct error.
  LENLIST = strlen(INPUTS.HOSTLIB_STOREPAR_LIST);
  if ( LENLIST > MXPATHLEN  ) {
    sprintf(c1err, "len(HOSTLIB_STOREPAR) = %d exceeds bound of %d",
	   LENLIST, MXPATHLEN) ;
    sprintf(c2err, "Remove some variables from HOSTLIB_STOREPAR");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

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

    if ( LDMP ) {
      printf(" xxx %s: --------------------------- \n", fnam);
      printf(" xxx %s: ivar=%2d varName='%s' \n", fnam, ivar, varName);
    }

    checkAlternateVarNames_HOSTLIB(varName);
    
    // skip duplicates, which could come from zHOST effic map,
    // or user mistake setting up STOREPAR_LIST
    ISDUPL = 0 ;
    for(ivar2=0; ivar2 < ivar; ivar2++ ) {
      varName2 = HOSTLIB_OUTVAR_EXTRA.NAME[ivar2];
      if ( QstringMatch(varName,varName2) ) { ISDUPL = 1 ; }
    }
    if(LDMP)
      { printf(" xxx %s: ISDUPL=%d for varName='%s' \n",fnam,ISDUPL,varName);}

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

    if(LDMP)
      { printf(" xxx %s: LOAD=%d for varName='%s' \n",fnam,LOAD,varName);}

    // always store variable in OUTVAR list, along with IVAR
    sprintf(HOSTLIB_OUTVAR_EXTRA.NAME[NVAR_OUT], "%s", varName );
    HOSTLIB_OUTVAR_EXTRA.IVAR_STORE[NVAR_OUT] = IVAR_HOSTLIB(varName,1);

    // set USED_IN_WGTMAP to zero; fill this later in read_head_HOSTLIB.
    HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[NVAR_OUT] = 0 ;

    NVAR_OUT++ ;  

    catVarList_with_comma(VARLIST_LOAD,varName); // for print only

  } // end ivar loop


  // update global counters
  HOSTLIB_OUTVAR_EXTRA.NOUT = NVAR_OUT ;
  HOSTLIB.NVAR_REQUIRED     = NVAR_REQ ;

  // - - - - - -
  printf("\t Load '%s' for SNTABLE storage. \n", VARLIST_LOAD);
  fflush(stdout);

  return ;

} // end of  init_OUTVAR_HOSTLIB

// ===================================================
bool QstringMatch(char *varName0, char *varName1 ) {
  // May 24 2020
  // check for string match after removing possible question mark
  // at end of string.
  // Examples:
  //
  //  varName0  varName1  return
  //  apple      apple     true
  //  apple      apple?    true
  //  apple?     apple     true
  //  apple?     apple?    true
  //  apple      app       false
  // 

  char tmp0[100], tmp1[100];
  int  len0, len1;
  // --------- BEGIN ---------
  sprintf(tmp0, "%s", varName0);   len0 = strlen(tmp0);
  sprintf(tmp1, "%s", varName1);   len1 = strlen(tmp1);
  if( tmp0[len0-1] == '?' ) { tmp0[len0-1] = 0; }
  if( tmp1[len1-1] == '?' ) { tmp1[len1-1] = 0; }

  if ( strcmp(tmp0,tmp1) == 0 ) 
    { return true; }
  else
    { return false ; }
  
} // end QstringMatch

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
  char fnam[] = "open_HOSTLIB" ;

  // ----------- BEGIN ----------

  *fp = snana_openTextFile(OPTMASK_OPENFILE_HOSTLIB, 
			   PATH_DEFAULT_HOSTLIB, INPUTS.HOSTLIB_FILE,
			   libname_full, &HOSTLIB.GZIPFLAG );  // <== returned

  if ( *fp == NULL ) {
    abort_openTextFile("HOSTLIB_FILE", 
		       PATH_DEFAULT_HOSTLIB, INPUTS.HOSTLIB_FILE, fnam);
  }

  sprintf(HOSTLIB.FILENAME , "%s", libname_full );
  printf("\t Reading %s \n", libname_full );
  fflush(stdout);

}  // end of open_HOSTLIB

void close_HOSTLIB(FILE *fp) {
  // Created July 16 2021

  char fnam[] = "close_HOSTLIB";

  check_EOF(fp, INPUTS.HOSTLIB_FILE, fnam, HOSTLIB.NGAL_READ);

  if ( HOSTLIB.GZIPFLAG ) 
    { pclose(fp); } // close gzip file
  else
    { fclose(fp); } // close normal file stream
} // end close_HOSTLIB

// ====================================
void  read_HOSTLIB_WGTMAP(void) {

  // Function to read OPTIONAL weight-map to over-ride
  // weight map in the HOSTLIB. If the weight map is read
  // here, then the corresponding weight-map in the HOSTLIB
  // will be ignored. Note that this function must be called
  // before read_head_HOSTLIB().
  //
  // July 14 2020: replace PATH_USER_INPUT with PATH_DEFAULT_HOSTLIB

  FILE *fp ;
  int  gzipFlag ;
  char *ptrFile, fileName_full[MXPATHLEN], c_get[200] ;
  char fnam[] = "read_HOSTLIB_WGTMAP"  ;

  // ------------- BEGIN --------------

  HOSTLIB_WGTMAP.READSTAT = false ;

  ptrFile = INPUTS.HOSTLIB_WGTMAP_FILE ;
  if ( IGNOREFILE(ptrFile) )  { return ; }

  fp = snana_openTextFile(OPTMASK_OPENFILE_HOSTLIB, 
			  PATH_DEFAULT_HOSTLIB, ptrFile,
			  fileName_full, &gzipFlag );  // <== returned

  if ( !fp ) {
      abort_openTextFile("HOSTLIB_WGTMAP_FILE", 
			 PATH_DEFAULT_HOSTLIB, ptrFile, fnam);
  }

  // if we get here, open and read WGTMAP file.

  printf("\t Read WEIGHT-MAP from supplemental file:\n\t   %s\n", ptrFile );
  fflush(stdout);

  while( (fscanf(fp, "%s", c_get)) != EOF) 
    { parse_HOSTLIB_WGTMAP(fp,c_get);  }

  if ( gzipFlag ) { pclose(fp); }   else { fclose(fp); }

  HOSTLIB_WGTMAP.READSTAT = true ;

  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SNMAGSHIFT )
    { printf("\t Implement SNMAGSHIFT in WGTMAP \n"); fflush(stdout); }
  else
    { printf("\t Ignore SNMAGSHIFT in WGTMAP \n"); fflush(stdout); }


} // end of read_HOSTLIB_WGTMAP


// ====================================
void parse_HOSTLIB_WGTMAP(FILE *fp, char *string) {

  // Parse WGTMAP variables from file *fp.
  // *string is the current string value to check
  // if this is one of the WGTMAP keys.
  //
  // Mar 14 2019: refactor to use read_GRIDMAP().
  // Apr 12 2019: return if string != VARNAMES_WGTMAP
  // Jun 20 2020: fix what seems like a cut-and-paste bug for N_SNVAR
  // Nov 18 2021: check HOSTLIB_WGTMAP.FOUNDVAR_SNMAGSHIFT

  int  IDMAP = IDGRIDMAP_HOSTLIB_WGTMAP ;
  long long GALID ;
  bool FOUND_VARNAMES, IS_SNVAR, IS_STORED ;
  int NVAR_WGTMAP, IVAR_STORE, NDIM, NFUN, ivar, N, N_SNVAR=0 ;

  char LINE[100], *VARNAME ;
  char fnam[] = "parse_HOSTLIB_WGTMAP"  ;

  // ----------- BEGIN -------------

  if ( strcmp(string,"NVAR_WGTMAP:")==0 ) {
    printf("\n WARNING: Should remove obsolete "
	   "NVAR_WGTMAP key from %s\n", INPUTS.HOSTLIB_FILE );
  }

  IVAR_STORE = HOSTLIB.NVAR_STORE ;

  if ( strcmp(string,"OPT_EXTRAP_WGTMAP:") == 0 ) 
    { HOSTLIB_WGTMAP.OPT_EXTRAP = 1;  } // Jun 11 2021

  FOUND_VARNAMES = ( strcmp(string,"VARNAMES_WGTMAP:") ==0 );
  if ( !FOUND_VARNAMES) { return ; }
 
  fgets(LINE,100,fp);
  NVAR_WGTMAP = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);

  NFUN = 1;    // WGT is required fun-val
  if ( HOSTLIB_WGTMAP.FOUNDVAR_SNMAGSHIFT )  { NFUN++; }
  NDIM = NVAR_WGTMAP-NFUN ; 

  if ( NDIM < 1 ) {
    sprintf(c1err, "Invalid NDIM=%d for %s", NDIM, string);
    sprintf(c2err, "VARNAMES_WGTMAP must inclulde WGT & SNMAGSHIFT");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  
  // - - - - - - -  -
  // read in names used for weight map
  HOSTLIB_WGTMAP.GRIDMAP.VARLIST[0] = 0 ;

  for ( ivar=0; ivar < NVAR_WGTMAP ; ivar++ ) {
    VARNAME = HOSTLIB_WGTMAP.VARNAME[ivar] ;
    get_PARSE_WORD(0,ivar,VARNAME) ;

    checkAlternateVarNames_HOSTLIB(VARNAME); // Jan 31 2020

    // check SN properties (e..g, x1, c) that are not in HOSTLIB
    IS_SNVAR = checkSNvar_HOSTLIB_WGTMAP(VARNAME); // Mar 2020
    HOSTLIB_WGTMAP.IS_SNVAR[ivar] = IS_SNVAR ; 
    if ( IS_SNVAR ) { 
      N = HOSTLIB_WGTMAP.N_SNVAR;
      HOSTLIB_WGTMAP.ISPARSE_SNVAR[N]       = ivar ;
      HOSTLIB_WGTMAP.INVSPARSE_SNVAR[ivar]  = N ;
      HOSTLIB_WGTMAP.N_SNVAR++ ; 
      N_SNVAR = HOSTLIB_WGTMAP.N_SNVAR ;
    }

    strcat(HOSTLIB_WGTMAP.GRIDMAP.VARLIST,VARNAME) ;
    strcat(HOSTLIB_WGTMAP.GRIDMAP.VARLIST," ") ;
    
    // load variable if it's not already loaded, and NOT SN var.
    IS_STORED = (IVAR_HOSTLIB(VARNAME,0) >= 0 ); 
    if ( !IS_STORED && !IS_SNVAR && ivar < NDIM ) {
      sprintf(HOSTLIB.VARNAME_STORE[IVAR_STORE], "%s", VARNAME );
      IVAR_STORE++ ;
    }  
  } // end of ivar loop
  
  // read WGT keys and load GRIDMAP struct.
  read_GRIDMAP(fp, "WGTMAP", "WGT:", "", IDMAP, NDIM, NFUN, 
	       HOSTLIB_WGTMAP.OPT_EXTRAP,
	       MXROW_WGTMAP, fnam,
	       &HOSTLIB_WGTMAP.GRIDMAP ); // <== return GRIDMAP
  
  HOSTLIB_WGTMAP.WGTMAX = HOSTLIB_WGTMAP.GRIDMAP.FUNMAX[0];

  // update global counter
  HOSTLIB.NVAR_STORE = IVAR_STORE ;


  // check for optional WGTMAP_CHECK key to verify
  // WGTMAP interpolation.
  double TMPVAL[10];
  if ( strcmp(string,"WGTMAP_CHECK:") == 0  && N_SNVAR==0  ) {
    readlong  (fp, 1, &GALID ); // Feb 2015
    readdouble(fp, 2, TMPVAL  );
    N = HOSTLIB_WGTMAP.NCHECKLIST ;
    HOSTLIB_WGTMAP.CHECKLIST_GALID[N] = GALID ;
    HOSTLIB_WGTMAP.CHECKLIST_ZTRUE[N] = TMPVAL[0] ;
    HOSTLIB_WGTMAP.CHECKLIST_WGT[N]   = TMPVAL[1] ;
    HOSTLIB_WGTMAP.CHECKLIST_SNMAG[N] = TMPVAL[2] ;
    HOSTLIB_WGTMAP.NCHECKLIST++ ;  
  }

  // - - - - - - -

  prep_SNVAR_HOSTLIB_WGTMAP();

  return ;

} // end of parse_HOSTLIB_WGTMAP

// ============================================
int read_VARNAMES_WGTMAP(char *VARLIST_WGTMAP) {

  // July 14 2020
  // pre-HOSTLIB-read utility to fetch & return list of
  // VARNAMES appearing in WGTMAP so that WGTMAP variables can
  // be automatically stored in data files.
  // Check external WGTMAP file first; if no external WGTMAP file,
  // then check if WGTMAP is embedded in HOSTLIB.
  // Function returns number of WGTMAP variables found.
  //
  // Jan 15 2021: abort if EOF is reached.
  // Nov 18 2021: 
  //   + abort if WGT is not found
  //   + set FOUNDVAR_SNMAGSHIFT 
  //

  int NVAR   = 0 ;
  int MXCHAR = MXPATHLEN;
  int NWD, ivar, gzipFlag, NTMP=0 ;
  FILE *fp ;
  bool IS_SNVAR, FOUNDVAR_WGT=false, FOUNDVAR_SNMAGSHIFT=false;
  char FILENAME_FULL[MXPATHLEN], *WGTMAP_FILE, LINE[MXPATHLEN];
  char c_get[60], VARNAME[60];
  char KEY_VARNAMES[] = "VARNAMES_WGTMAP:" ;
  char KEY_STOP[]     = "GAL:" ; // stop reading when this key is found
  char fnam[]         = "read_VARNAMES_WGTMAP" ;
  int  LDMP = 0 ;
  // ------------- BEGIN ------------

  if ( IGNOREFILE(INPUTS.HOSTLIB_FILE) ) { return(0) ; }

  VARLIST_WGTMAP[0] = 0 ;

  if ( !IGNOREFILE(INPUTS.HOSTLIB_WGTMAP_FILE)  )
    { WGTMAP_FILE = INPUTS.HOSTLIB_WGTMAP_FILE ; } // external WGTMAP 
  else
    { WGTMAP_FILE = INPUTS.HOSTLIB_FILE ; } // ...or check inside HOSTLIB
  
  if ( LDMP ) {
    printf(" xxx \n" ) ;
    printf(" xxx %s: PATH = '%s' \n", fnam, PATH_DEFAULT_HOSTLIB );  
    printf(" xxx %s: WGTMAP_FILE = '%s' \n", fnam, WGTMAP_FILE );
    printf(" xxx \n" ) ;
  }

  // open file containing WGTMAP
  fp = snana_openTextFile(OPTMASK_OPENFILE_HOSTLIB, 
			  PATH_DEFAULT_HOSTLIB, WGTMAP_FILE,
			  FILENAME_FULL, &gzipFlag );  // <== returned
  
  if ( !fp ) {
    sprintf(c1err, "Unable to open WGTMAP file (to read VARNAMES)");
    sprintf(c2err, "'%s'", WGTMAP_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  bool STOP_READ = false;
  while( !STOP_READ ) { 

    if ( fscanf(fp, "%s", c_get) == EOF ) {
      sprintf(c1err,"Reached EOF before finding WGTMAP or GAL key.");
      sprintf(c2err,"Check format for HOSTLIB_FILE.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    NTMP++;
    //    if ( NTMP < 200 ) 
    //{  printf(" xxx %s: c_get(%3d)='%s' \n", fnam, NTMP, c_get);  }

    // avoid reading the entire HOSTLIB file if there is no WGTMAP
    if ( strcmp(c_get,KEY_STOP) == 0 ) { STOP_READ = true; }

    if ( strcmp(c_get,KEY_VARNAMES) == 0 ) {
      STOP_READ = true ;
      fgets(LINE, MXCHAR, fp);
      NWD  = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
     
      for(ivar=0; ivar < NWD; ivar++ ) {
	get_PARSE_WORD(0,ivar,VARNAME);

	if ( strcmp(VARNAME,"WGT") == 0 )
	  { FOUNDVAR_WGT = true; continue; }

	if ( strcmp(VARNAME,HOSTLIB_VARNAME_SNMAGSHIFT) == 0 )
	  { FOUNDVAR_SNMAGSHIFT = true; continue; }

	NVAR++ ;
	IS_SNVAR = checkSNvar_HOSTLIB_WGTMAP(VARNAME);
	if ( !IS_SNVAR ) 
	  { catVarList_with_comma(VARLIST_WGTMAP,VARNAME); }
	if ( LDMP ) 
	  { printf(" xxx %s wgtmap var '%s'  (storeFlag=%d)\n", 
		   fnam, VARNAME, !IS_SNVAR ); }
      }

    } // end reading VARNAMES_WGTMAP line
  } // end STOP_READ

  // - - - - -  
  if ( NVAR > 0 && !FOUNDVAR_WGT ) {
    sprintf(c1err,"Missing required WGT column in WGTMAP");
    sprintf(c2err,"Check WGTMAP");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  HOSTLIB_WGTMAP.FOUNDVAR_SNMAGSHIFT = FOUNDVAR_SNMAGSHIFT;

  if ( gzipFlag ){ pclose(fp); }     else { fclose(fp); }

  return(NVAR) ;

} // end read_VARNAMES_WGTMAP

// ====================================
void prep_SNVAR_HOSTLIB_WGTMAP(void) {

  // Mar 13 2020
  // prepare optional SN variables in WGTMAP that are NOT in HOSTLIB.
  // Beware: indexing is tricky.

  int  N_SNVAR = HOSTLIB_WGTMAP.N_SNVAR ;
  int  ivar_SN, ivar, NBIN, NBTOT, IBIN, IBIN_TMP, NB1D, IB1D, ID  ;
  double VALMIN, VALBIN, VALGRID ;
  char *VARNAME, VARLIST[100] ;
  char fnam[] = "prep_SNVAR_HOSTLIB_WGTMAP" ;
  bool LDMP = false;
  // ------------- BEGIN -----------

  if ( N_SNVAR == 0 ) { return ; }
 

  VARLIST[0] = 0 ;
  for(ivar_SN = 0; ivar_SN < N_SNVAR; ivar_SN++ ) {
    ivar    = HOSTLIB_WGTMAP.ISPARSE_SNVAR[ivar_SN];
    VARNAME = HOSTLIB_WGTMAP.VARNAME[ivar] ;
    strcat(VARLIST,VARNAME);
    strcat(VARLIST," " );
    NBIN    = HOSTLIB_WGTMAP.GRIDMAP.NBIN[ivar];

    HOSTLIB_WGTMAP.NB1D_SNVAR[ivar_SN] = NBIN;
    HOSTLIB_WGTMAP.NBTOT_SNVAR *= NBIN;
    sprintf(HOSTLIB_WGTMAP.VARNAME_SNVAR[ivar_SN],"%s", VARNAME);
    // printf(" xxx %s: NBIN[%d,%s] = %d \n", fnam, ivar_SN, VARNAME, NBIN);
  }
  NBTOT = HOSTLIB_WGTMAP.NBTOT_SNVAR;

  printf("\t %s: SNVAR = '%s' -> %d %dD bins \n",
	 fnam, VARLIST, NBTOT, N_SNVAR ); 
  fflush(stdout);

  // setup bin map from global 1D bin to 1D bin for each SNvar
  int MEMI = NBTOT * sizeof(int) ;
  int MEMD = NBTOT * sizeof(double) ;
  for(ivar_SN = 0; ivar_SN < N_SNVAR; ivar_SN++ ) { 
    HOSTLIB_WGTMAP.IBIN1D_SNVAR[ivar_SN]  = (int   *) malloc(MEMI); 
    HOSTLIB_WGTMAP.VALGRID_SNVAR[ivar_SN] = (double*) malloc(MEMD); 
  }

  for(IBIN=0; IBIN < HOSTLIB_WGTMAP.NBTOT_SNVAR; IBIN++ ) {
    IBIN_TMP = IBIN ;
    NBTOT = HOSTLIB_WGTMAP.NBTOT_SNVAR;
    for(ivar_SN = 0; ivar_SN < N_SNVAR; ivar_SN++ ) {
      NB1D = HOSTLIB_WGTMAP.NB1D_SNVAR[ivar_SN] ;
      NBTOT   /= NB1D ;
      IB1D     = IBIN_TMP/NBTOT ;
      IBIN_TMP = ( IBIN_TMP % NBTOT);
      HOSTLIB_WGTMAP.IBIN1D_SNVAR[ivar_SN][IBIN] = IB1D ;

      // store SN value(s) at each grid point
      ivar    = HOSTLIB_WGTMAP.ISPARSE_SNVAR[ivar_SN];
      VALMIN  = HOSTLIB_WGTMAP.GRIDMAP.VALMIN[ivar];
      VALBIN  = HOSTLIB_WGTMAP.GRIDMAP.VALBIN[ivar];
      VALGRID = VALMIN + VALBIN * (double)IB1D;
      HOSTLIB_WGTMAP.VALGRID_SNVAR[ivar_SN][IBIN] = VALGRID ;
	

      if ( LDMP ) {
	VARNAME = HOSTLIB_WGTMAP.VARNAME_SNVAR[ivar_SN] ;
	printf(" xxx %s: IBIN=%3d -> %s bin = %d\n",
	       fnam, IBIN, VARNAME, IB1D);	     
      }
    }
  }

  // - - - - - - - - - - - 
  // init map to convert multiple 1D indices into a single index
  // for 1D array.
  HOSTLIB_WGTMAP.IDMAP_INDEX_SNVAR = IDGRIDMAP_HOSTLIB_WGTMAP + 1 ;
  ID = HOSTLIB_WGTMAP.IDMAP_INDEX_SNVAR ;
  init_1DINDEX(ID, N_SNVAR, HOSTLIB_WGTMAP.NB1D_SNVAR );

  return ;

} // end prep_SNVAR_HOSTLIB_WGTMAP

// =======================================
void getVal_SNVAR_HOSTLIB_WGTMAP(int ibin, double *VAL_WGTMAP) {

  // Created Mar 14 2020
  // For input SNVar "ibin", return SN value for each SNVar
  // that is NOT in the HOSTLIB.
  //
  // Example with x1(5 bins) and c(3 bins) so that ibin can take
  // and value from 0 to 14. 
  // Ouptut VAL_WGTMAP[0] is the x1 value for the weigt map, 
  // and VAL_WGTMAP[1] is the c value.

  int N_SNVAR     = HOSTLIB_WGTMAP.N_SNVAR;
  int ivar_SN ;
  double VAL;
  bool LDMP = false ;
  //  char fnam[] = "getVal_SNVAR_HOSTLIB_WGTMAP" ;

  // ------------ BEGIN ------------

  for(ivar_SN=0; ivar_SN < N_SNVAR; ivar_SN++ ) {
    VAL = HOSTLIB_WGTMAP.VALGRID_SNVAR[ivar_SN][ibin] ;
    VAL_WGTMAP[ivar_SN] = VAL;
    if ( LDMP ) {
      char *VARNAME = HOSTLIB_WGTMAP.VARNAME_SNVAR[ivar_SN] ;
      printf(" xxx %s: ibin=%3d -> %s = %8.3f \n", 
	     "getVal_SNVAR", ibin, VARNAME, VAL); fflush(stdout);
    }
  }

  return ;

} // end getVal_SNVAR_HOSTLIB_WGTMAP


// ========================================
int getBin_SNVAR_HOSTLIB_WGTMAP(void) {

  // Mar 14 2020
  // Called for each event, return SNVAR ibin for current
  // ptrVal_SNVAR values. Note that here we find ibin
  // such that SNVAR values have closest match to 
  // current generated values in ptrVal_SNVAR.
  //
  // Beware that uniform bins are required in each SNVAR dimension.
  // Jun 20 2020: abort of SN value is outside WGTMAP range.

  int  NBTOT_SNVAR = HOSTLIB_WGTMAP.NBTOT_SNVAR;
  int  N_SNVAR     = HOSTLIB_WGTMAP.N_SNVAR;
  int  ID          = HOSTLIB_WGTMAP.IDMAP_INDEX_SNVAR ;
  int  ivar_SN, ivar, IBIN = 0 ;
  int  IB1D[MXVAR_HOSTLIB];
  double VAL, VALMIN, VALMAX, VALBIN;
  char *VARNAME ;
  char fnam[] = "getBin_SNVAR_HOSTLIB_WGTMAP" ;
  bool LDMP   = false ;

  // ------------- BEGIN ------------
  
  if ( N_SNVAR == 0 ) { return(IBIN); }

  if ( LDMP ) { 
    printf(" xxx -------------------------------------------- \n");
    printf(" xxx    %s DUMP \n", fnam ); 
  }

  for(ivar_SN=0; ivar_SN < N_SNVAR; ivar_SN++ ) {
    ivar    =  HOSTLIB_WGTMAP.ISPARSE_SNVAR[ivar_SN];
    VALBIN  =  HOSTLIB_WGTMAP.GRIDMAP.VALBIN[ivar] ; 
    VALMIN  =  HOSTLIB_WGTMAP.GRIDMAP.VALMIN[ivar] ; 
    VALMAX  =  HOSTLIB_WGTMAP.GRIDMAP.VALMAX[ivar] ; 
    VAL     = *HOSTLIB_WGTMAP.ptrVal_SNVAR[ivar_SN] ;
    VARNAME = HOSTLIB_WGTMAP.VARNAME_SNVAR[ivar_SN]; 
    IB1D[ivar_SN] = (int)(0.5+(VAL-VALMIN)/VALBIN) ; 

    if ( VAL < VALMIN || VAL > VALMAX ) {
      sprintf(c1err,"%s = %.4f is outside WGTMAP range",  VARNAME, VAL);
      sprintf(c2err,"%.4f <= %s <= %.4f in WGTMAP", VALMIN,VARNAME,VALMAX);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( LDMP ) {
      printf(" xxx %3s = %8.4f (ivar=%d, MIN/MAX=%.3f/%.3f "
	     "BIN=%.3f, bin1D=%d) \n", 
	     VARNAME, VAL, ivar, VALMIN, VALMAX, VALBIN, IB1D[ivar_SN] );
    }
  }


  IBIN = get_1DINDEX(ID, N_SNVAR, IB1D);

  if ( LDMP ) { printf(" xxx     IBIN = %d \n", IBIN); fflush(stdout); }

  if ( IBIN < 0 || IBIN >= NBTOT_SNVAR ) {
    sprintf(c1err, "Invalid IBIN=%d (valid range: 0 to %d", IBIN,NBTOT_SNVAR);
    sprintf(c2err, "SNVAR is messed up.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return(IBIN) ;

} // end getBin_SNVAR_HOSTLIB_WGTMAP

// ====================================
void  read_specTable_HOSTLIB(void) {

  // Created Jun 28 2019 by R.Kessler
  // Read supplemental file of spectral templates, used later
  // to construct Host spectrum for each event.
  // This read must be done before reading the HOSTLIB,
  // so that the template names can be matched between
  // this specTemplate file and the VARNAMES in the HOSTLIB.
  //
  // May 22 2020: abort if VARNAMES key not found
  // Feb 12 2021: adapt to also work for SPECDATA 

  FILE *fp;
  int  NBIN_WAVE, NBIN_READ, IFILETYPE;
  int  NVAR, ivar, ICOL_WAVE, MEMD, NUM, gzipFlag, ifile ;
  int  NVAR_WAVE  = 0 ;
  int  OPT_VARDEF = 0 ;
  int  LEN_PREFIX,  *NSPEC ;
  char *ptrFile_list[2], *ptrFile ;
  char *varName, varName_tmp[100], VARNAME_PREFIX[20];
  char c_get[60], fileName_full[MXPATHLEN] ;  
  char TABLENAME_LIST[2][12] = { "BASIS", "DATA" } ;
  char PREFIX_LIST[2][12]    = { PREFIX_SPECBASIS, PREFIX_SPECDATA } ;
  char TBLNAME[20];
  char fnam[] = "read_specTable_HOSTLIB";
  
  // --------------- BEGIN -----------------

  HOSTSPEC.ITABLE       = -9 ; // BASIS or DATA
  HOSTSPEC.TABLENAME[0] =  0 ;
  HOSTSPEC.NSPECBASIS   =  0 ;
  HOSTSPEC.NSPECDATA    =  0 ;
  HOSTSPEC.IDSPECDATA   = -9 ;
  HOSTSPEC.NBIN_WAVE    =  0 ;
  HOSTSPEC.FLAM_SCALE        = 1.0 ;
  HOSTSPEC.FLAM_SCALE_POWZ1  = 0.0 ;

  ptrFile_list[ITABLE_SPECBASIS] = INPUTS.HOSTLIB_SPECBASIS_FILE ;
  ptrFile_list[ITABLE_SPECDATA]  = INPUTS.HOSTLIB_SPECDATA_FILE ;
  for(ifile = 0; ifile < 2; ifile++ ) {
    if ( !IGNOREFILE(ptrFile_list[ifile]) ) {
      HOSTSPEC.ITABLE = ifile;
      sprintf(HOSTSPEC.TABLENAME, "%s", TABLENAME_LIST[ifile] );
      sprintf(VARNAME_PREFIX,"%s", PREFIX_LIST[ifile]);
      LEN_PREFIX = strlen(VARNAME_PREFIX);
      ptrFile = ptrFile_list[ifile];
      sprintf(TBLNAME,"%s", HOSTSPEC.TABLENAME);
    }
  }

  if ( HOSTSPEC.ITABLE < 0 ) { return; } // nothing to do here

  printf("\n\t Read SPEC-%s from supplemental file:\n", HOSTSPEC.TABLENAME );
  fflush(stdout);

  // - - - - - - - - - - - - - -
  // read until VARNAMES key in case there are supplemental keys

  fp = snana_openTextFile(OPTMASK_OPENFILE_HOSTLIB, 
			  PATH_USER_INPUT, ptrFile,
			  fileName_full, &gzipFlag );  // <== returned
  if ( !fp ) {
    sprintf(varName_tmp,"HOSTLIB_SPEC%_FILE", HOSTSPEC.TABLENAME );
    abort_openTextFile(varName_tmp, PATH_USER_INPUT, ptrFile, fnam);
  }

  if ( HOSTSPEC.ITABLE == ITABLE_SPECBASIS )
    { NSPEC = &HOSTSPEC.NSPECBASIS ; }
  else
    { NSPEC = &HOSTSPEC.NSPECDATA ; }


  bool FOUND_VARNAMES = false ;
  while( !FOUND_VARNAMES  && (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"VARNAMES:") == 0 ) { FOUND_VARNAMES=true ; }

    if ( strcmp(c_get,"FLAM_SCALE:") == 0 ) 
      { readdouble(fp, 1, &HOSTSPEC.FLAM_SCALE); }    

    if ( strcmp(c_get,"FLAM_SCALE_POWZ1:") == 0 ) 
      { readdouble(fp, 1, &HOSTSPEC.FLAM_SCALE_POWZ1); }    
  }

  if ( gzipFlag ) { pclose(fp); }  else {  fclose(fp); }

  if ( !FOUND_VARNAMES ) {
    sprintf(c1err,"Could not find VARNAMES in header of");
    sprintf(c2err,"HOSTLIB_SPEC%s_FILE", TBLNAME );
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err); 
  }

  // - - - - - - - - - - - - - -
  // now read spec-table with standard routines
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

  // examine VARNAMES list to make sure that there is a wavelength column,
  // and count how many template[nn] colummns
  ICOL_WAVE = -9;  *NSPEC=0 ;
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

    sprintf(varName_tmp,"%s", varName);
    varName_tmp[LEN_PREFIX] = 0;

    if ( strcmp_ignoreCase(varName_tmp,VARNAME_PREFIX) == 0 ) {
      int N = *NSPEC;
      if ( N < MXSPECBASIS_HOSTLIB ) {
	sscanf(&varName[LEN_PREFIX], "%d",  &NUM);

	HOSTSPEC.ICOL_SPECTABLE[N] = ivar;
	HOSTSPEC.NUM_SPECBASIS[N]  = NUM; // store template NUM
	sprintf(HOSTSPEC.VARNAME_SPECBASIS[N],"%s", varName);

	HOSTSPEC.FLAM_BASIS[N] = (double*) malloc(MEMD);
	SNTABLE_READPREP_VARDEF(varName,HOSTSPEC.FLAM_BASIS[N],
				NBIN_WAVE, OPT_VARDEF);
      }
      (*NSPEC)++ ; // always increment number of templates or data 
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

  if ( *NSPEC >= MXSPECBASIS_HOSTLIB ) {
    sprintf(c1err,"NSPEC%s=%d exceeds bound of %d", 
	    TBLNAME, *NSPEC, MXSPECBASIS_HOSTLIB ) ;
    sprintf(c2err,"Remove templates or increase MXSPECBASIS_HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err);     
  }
	
  HOSTSPEC.ICOL_WAVE = ICOL_WAVE ;

  // read the entire table, and close it.
  NBIN_READ = SNTABLE_READ_EXEC();

  // Loop over wave bins and
  // + determine WAVE_BINSIZE for each wave bin
  // + truncate NBIN_WAVE so that lam < MAXLAM_SEDMODEL
  double WAVE_BINSIZE, WAVE_MIN, WAVE_MAX, LAM ;
  double LAM_LAST=0.0, LAM_NEXT=0.0;
  int FIRST, LAST, ilam, NBIN_KEEP=0 ;
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

  printf("\t Found %d spectral %s vectors and %d wavelength bins\n",
	 *NSPEC, HOSTSPEC.TABLENAME, NBIN_WAVE);
  fflush(stdout);

  return;

} // end read_specTable_HOSTLIB

// ============================================
void match_specTable_HOSTVAR(void) {

  // Created June 2019
  // Match specTemplate names to names in HOSTLIB (HOSTVAR).
  // The specTemplate file has column names specbasis00, specbasis01, etc ...
  // The HOSTLIB VARNAMES must have corresponding list of
  // coeff_specbasis00, coeff_specbasis01, etc ...
  //
  // Output is to load global HOSTSPEC.IVAR_HOSTLIB[i]
  //
  // Feb 23 2021: adapt for SPECDATA

  int  ITABLE        = HOSTSPEC.ITABLE ;
  int  NSPECBASIS    = HOSTSPEC.NSPECBASIS ;
  int  ivar_HOSTLIB, i, imin, imax, NERR=0;
  int  LDMP = 0 ;
  char *VARNAME_SPECBASIS, VARNAME_HOSTLIB[40];
  char fnam[] = "match_specTable_HOSTVAR";

  // ----------------- BEGIN -----------------

  if ( ITABLE == ITABLE_SPECDATA ) { 
    HOSTSPEC.IVAR_HOSTLIB[0] = IVAR_HOSTLIB(VARNAME_SPECDATA_HOSTLIB,0);
    return ;
  }
  
  if ( NSPECBASIS == 0 ) { return ; }

  // - - - -

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
    sprintf(c1err,"%d specTemplates have no coeff_template in HOSTLIB:", 
	    NERR );
    sprintf(c2err,"Check specTemplate file and HOSTLIB VARNAMES");
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err); 
  }

  return ;

} // end match_specTable_HOSTVAR

// =========================================
void checkVarName_specTemplate(char *varName) {

  // Created Jun 2019
  // ABORT if varName (from HOSTLIB VARNAMES) is an 
  // un-used template coefficient.
  // 
  int icol;
  char fnam[] = "checkVarName_specTemplate";

  // --------------- BEGIN ------------

  if ( HOSTSPEC.ITABLE == ITABLE_SPECDATA ) { return; }
  if ( HOSTSPEC.NSPECBASIS <= 0 ) { return; }

  if (strstr(varName,PREFIX_SPECBASIS)         == NULL ) { return ; }
  if (strstr(varName,PREFIX_SPECBASIS_HOSTLIB) == NULL ) { return ; }

  // now we have varName of the form coeff_template[nn]
  icol = ICOL_SPECTABLE(varName,0); // fetch column in specTemplate file
  
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
  // Return true host spectrum (no noise), including Galactic extinction.
  // If option is set to compute broadband mags, load them
  // into SNHOSTGAL.GALMAG[ifilt_obs][0] 
  //
  // Issues:
  //  - fraction of galaxy in fiber or slit ? Or does ETC include this ?
  //
  // Call fill_TABLE_MWXT_SEDMODEL
  //
  // May 6 2021: check ABMAG_FORCE
  // Dec 17 2021: exclude last LAM bin from LAMBIN_CHECK test;
  //             -> avoids mysterious abort.
  //

  int  NBLAM_SPECTRO    = INPUTS_SPECTRO.NBIN_LAM;
  double ABMAG_FORCE    = INPUTS.HOSTLIB_ABMAG_FORCE;
  double ABMAG_SCALE    = 1.0/pow(10.0, 0.4*ABMAG_FORCE);
  bool   DO_ABMAG_FORCE = ( ABMAG_FORCE > -8.0 );

  int  NBLAM_BASIS     = HOSTSPEC.NBIN_WAVE; 
  bool IS_SPECBASIS    = HOSTSPEC.ITABLE == ITABLE_SPECBASIS ;
  bool IS_SPECDATA     = HOSTSPEC.ITABLE == ITABLE_SPECDATA ;
  double FLAM_SCALE_POWZ1 = HOSTSPEC.FLAM_SCALE_POWZ1 ;
  int  IGAL        = SNHOSTGAL.IGAL ;
  double z1        = 1.0 + zhel;
  if ( DO_ABMAG_FORCE ) { z1 = 1.0 ; }
  double znorm     = pow(z1,FLAM_SCALE_POWZ1) ;
  double hc8       = (double)hc;

  long long GALID;
  int  ilam, ilam_basis, ilam_last=-9, i, ivar_HOSTLIB, ivar ;
  int  ilam_near, NLAMSUM, LDMP=0;
  int  ISPEC_MIN, ISPEC_MAX, IDSPEC;
  double FLAM_TMP, FLAM_SUM, COEFF, FLUX_TMP, MWXT_FRAC, LAMOBS;
  double LAMOBS_BIN, LAMOBS_MIN, LAMOBS_MAX, LAM_BASIS;
  double LAMREST_MIN, LAMREST_MAX ;
  double LAMMIN_TMP, LAMMAX_TMP, LAMBIN_TMP, LAMBIN_CHECK ;
  double ZP, FTMP, MAG, ZCHECK;
  char fnam[] = "genSpec_HOSTLIB" ;

  // ------------------ BEGIN --------------

  fill_TABLE_MWXT_SEDMODEL(MWXT_SEDMODEL.RV, MWEBV);

  if ( HOSTSPEC.ITABLE < 0 ) {
    sprintf(c1err,"Cannot generate host spectrum without spec basis or data.");
    sprintf(c2err,"Check HOSTLIB_SPECBASIS[DATA]_FILE key in sim-input.");
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

  if ( IS_SPECBASIS ) {
    // check all basis vectors
    ISPEC_MIN = 0;  ISPEC_MAX = HOSTSPEC.NSPECBASIS;
  }
  else {
    // pick out the one spectrum in IDSPECDATA column of HOSTLIB
    ISPEC_MIN = 0;  ISPEC_MAX = 1;
    ivar_HOSTLIB = HOSTSPEC.IVAR_HOSTLIB[0];
    COEFF        = HOSTLIB.VALUE_ZSORTED[ivar_HOSTLIB][IGAL] ; 
    IDSPEC       = (int)COEFF ;
    HOSTSPEC.IDSPECDATA = IDSPEC ;
    if  ( IDSPEC < 0 || IDSPEC > HOSTSPEC.NSPECDATA ) {
      sprintf(c1err,"Invalid IDSPEC = %d (ivar_HOSTLIB=%d)", 
	      IDSPEC, ivar_HOSTLIB );
      sprintf(c2err,"Valid IDSPEC range is 0 to %d", HOSTSPEC.NSPECDATA-1);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  // construct total rest-frame spectrum using specbasis/specdata binning
  for(ilam_basis=0; ilam_basis < NBLAM_BASIS; ilam_basis++ ) {
    FLAM_SUM     = 0.0 ;
    LAM_BASIS    = HOSTSPEC.WAVE_CEN[ilam_basis]; // basis or data
    for(i=ISPEC_MIN; i < ISPEC_MAX; i++ ) {

      if ( IS_SPECBASIS ) { 
	ivar_HOSTLIB = HOSTSPEC.IVAR_HOSTLIB[i];
	COEFF        = HOSTLIB.VALUE_ZSORTED[ivar_HOSTLIB][IGAL] ; 
	FLAM_TMP     = HOSTSPEC.FLAM_BASIS[i][ilam_basis] ;
	FLAM_SUM    += ( COEFF * FLAM_TMP );
      }
      else {
	FLAM_SUM   = HOSTSPEC.FLAM_BASIS[IDSPEC][ilam_basis] ; 
      }
      if ( DUMPFLAG && ilam_basis == 0 && COEFF > 0.0 ) {
	printf(" xxx COEFF(%2d) = %le  (ivar_HOSTLIB=%d)\n", 
	       i, COEFF, ivar_HOSTLIB );
      }
    }

    // global scale for physical units
    HOSTSPEC.FLAM_EVT[ilam_basis] = (FLAM_SUM * HOSTSPEC.FLAM_SCALE * znorm);

    // May 2021: check force ABMAG (constant, NOT z-dependent)
    if ( DO_ABMAG_FORCE ) {
      FLAM_SUM = 3.631E-20 * LIGHT_A / (LAM_BASIS*LAM_BASIS) ;
      HOSTSPEC.FLAM_EVT[ilam_basis] = FLAM_SUM * ABMAG_SCALE ;
    }

  } // end ilam_basis


  // ---------------
  // HOST FLAM is in wavelength bins defined in the specTemplate file.
  // Here we to convert to SPECTROGRAPH bins in GENFLUX_LIST.
  // "ilam" is the index for SPECTROGRAPH, while ilam_basis is for specbasis.

  double LAMMIN_SPEC = HOSTSPEC.WAVE_MIN[0];
  double LAMMAX_SPEC = HOSTSPEC.WAVE_MIN[NBLAM_BASIS-1] ;

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

    // allow rest-frame spectra to be truncated
    if ( LAMREST_MIN < LAMMIN_SPEC ) { continue; } // Feb 25 2021
    if ( LAMREST_MAX > LAMMAX_SPEC ) { continue; }

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
      // i.e. exclude flux that leaks out of this SPECTROGRAPh bin.
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
    
    // if NLAMSUM==0, then no basis/data lambins overlap SPECTROGRAPH;
    // in thise case, average over nearest basis bins
    if ( NLAMSUM == 0 ) {
      double FLAM0 = HOSTSPEC.FLAM_EVT[ilam_near];
      double FLAM1 = HOSTSPEC.FLAM_EVT[ilam_near+1];
      FLAM_TMP     = (FLAM0+FLAM1)/2.0;
      FLUX_TMP     = FLAM_TMP * LAMOBS_BIN;
      LAMBIN_CHECK = LAMOBS_BIN/z1 ;
    }

    // check that sum over basis/data wave bins = SPECTROGRAPH bin size.
    // Basis wave is rest frame and SPECTROGRAPH bin is obs frame;
    // hence z1 factor needed to compare.

    bool LAST_LAM          = (ilam == NBLAM_SPECTRO-1);
    bool FAIL_LAMBIN_CHECK = ( fabs(LAMBIN_CHECK - LAMOBS_BIN/z1) > 0.001 );
    if ( FAIL_LAMBIN_CHECK && !LAST_LAM ) {
      print_preAbort_banner(fnam);
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

  // Mar 15 2019: ignore NVAR key, and parse entire line after VARNAMES
  // Mar 28 2019: use MXCHAR_LINE_HOSTLIB
  // Nov 11 2019: set HOSTLIB.IVAR_NBR_LIST
  // May 23 2020: 
  //   + read option VPEC_ERR
  //   + add error checking based on MSKOPT
  //
  // Jun 09 2021: use match_varname_HOSTLIB to apply logic for
  //              case-sensitive matching.
  //

  int MXCHAR    = MXCHAR_LINE_HOSTLIB ;
  bool DO_VPEC  = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEVPEC ) ;
  bool DO_RADEC = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC ) ;
  bool DO_SWAPZ = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SWAPZPHOT ) ;
  bool CHECK_DUPLICATE_COLUMNS = false; // set True after GHOSTLIBS are fixed
  int  LEN_PREFIX_ZPHOT_Q = strlen(HOSTLIB_PREFIX_ZPHOT_Q);

  int ivar, ivar_map, IVAR_STORE, i, N, NVAR, NVAR_WGTMAP, FOUND_SNPAR;
  int MATCH, NVAR_STORE_SNPAR, USE, IS_SNPAR, VBOSE ;
  int ivar2, NDUPL;
  bool FOUND_VARNAMES, FOUND_VPECERR;
  int NCHAR;
  char  key[40], c_get[200], c_var[100], ctmp[80], wd[20], *cptr, *cptr2 ;
  char *basename;
  char  LINE[MXCHAR_LINE_HOSTLIB];
  char  fnam[] = "read_head_HOSTLIB" ;

  // ------------- BEGIN ---------

  NVAR = NVAR_WGTMAP = 0 ;
  FOUND_VARNAMES     = false ;  // change to 1 after reading
  FOUND_VPECERR      = false ;
  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );
  NVAR_STORE_SNPAR = 0 ;
  HOSTLIB.FIX_VPEC_ERR = -9.0 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    // stop reading when first GAL: key is reached.
    if ( strcmp(c_get,"GAL:") == 0 )  {
      // rewind doesn't work for gzip file
      // xxxx snana_rewind(fp, INPUTS.HOSTLIB_FILE, HOSTLIB.GZIPFLAG );  

      fseek(fp,-4,SEEK_CUR); //rewind to before 1st GAL key
      goto VARCHECK ; 
    }

    if ( strcmp(c_get,"NVAR:")==0 ) { warn_NVAR_KEY(INPUTS.HOSTLIB_FILE); }

    sprintf(key, "%s", "VARNAMES:") ;
    if ( strcmp(c_get,key) == 0 ) {

      fgets(LINE,MXCHAR,fp);  NCHAR = strlen(LINE);

      if (NCHAR >= MXCHAR-5 ) {
	print_preAbort_banner(fnam);
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

      FOUND_VARNAMES = true ; 

      // parse the VARNAMES from LINE string
      for ( ivar=0; ivar < NVAR; ivar++ ) {
	get_PARSE_WORD(0,ivar,c_var);
	sprintf( HOSTLIB.VARNAME_ORIG[ivar], "%s", c_var); // 9.16.2021

	checkAlternateVarNames_HOSTLIB(c_var);

	// if coeff_template[nn], make sure it's actually needed
	checkVarName_specTemplate(c_var);

	// load ALL array
	sprintf( HOSTLIB.VARNAME_ALL[ivar], "%s", c_var);

	// check for duplicate columns (9.16.2021)  
	// E.g., ZERR and ZPHOTERR are both converted to ZPHOT_ERR.
	if ( CHECK_DUPLICATE_COLUMNS ) {
	  for(ivar2=0; ivar2 < ivar; ivar2++ ) {
	    cptr2 = HOSTLIB.VARNAME_ALL[ivar2];
	    if ( strcmp(c_var,cptr2) == 0 ) {
	      sprintf(c1err,"Found two '%s' columns at ivar=%d and %d", 
		      c_var, ivar2, ivar );
	      sprintf(c2err,"Original column names are %s and %s",
		      HOSTLIB.VARNAME_ORIG[ivar2],
		      HOSTLIB.VARNAME_ORIG[ivar] );
	      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
	    }
	  }
	} // end CHECK_DUPLICATE_COLUMNS

	// check for optional variables to add to the store list
	for ( i=0; i < HOSTLIB.NVAR_OPTIONAL; i++ ) {
	  cptr = HOSTLIB.VARNAME_OPTIONAL[i];
	  if( match_varname_HOSTLIB(c_var,cptr)  ) {
	    N = HOSTLIB.NVAR_STORE ;  // global var index (required+optional)
	    sprintf(HOSTLIB.VARNAME_STORE[N],"%s", c_var );
	    if ( HOSTLIB.IS_SNPAR_OPTIONAL[i] ) { 
	      HOSTLIB.FOUND_SNPAR_OPTIONAL[i] = 1 ;
	      NVAR_STORE_SNPAR++ ; 
	      HOSTLIB.IS_SNPAR_STORE[N] = 1 ; // index is for all variables
	    }

            if ( strstr(c_var,HOSTLIB_PREFIX_ZPHOT_Q) != NULL ) {
	      int   percentile, N_Q = HOSTLIB.NZPHOT_Q;
	      char *VARNAME = HOSTLIB.VARNAME_ZPHOT_Q[N_Q];
	      sscanf(&c_var[LEN_PREFIX_ZPHOT_Q], "%d", &percentile);   
	      // xxx mark delete sprintf(VARNAME, "%s", c_var);

	      // HOSTLIB ZPHOT_Qnn may not include pad zeros, but outpu
	      // VARNAME must include pad zeros.
	      sprintf(VARNAME, "%s%3.3d", HOSTLIB_PREFIX_ZPHOT_Q, percentile);
	      HOSTLIB.PERCENTILE_ZPHOT_Q[N_Q] = percentile ;
              HOSTLIB.NZPHOT_Q++ ;
            }

	    HOSTLIB.NVAR_STORE++ ;   
	  }
	}   // i       

      }  // end ivar loop

    } // end of VARNAMES


    // - - - - - 
    // look for fixed Sersic index 'n#_Sersic' outside of VARNAMES list
    if ( FOUND_VARNAMES ) 
      { parse_Sersic_n_fixed(fp,c_get); }

    // -----------
    // look for variables to use in weight-map 
    // (unless already read from elsewhere)

    if ( !HOSTLIB_WGTMAP.READSTAT  )
      { parse_HOSTLIB_WGTMAP(fp,c_get); }


    // check for fixed VPEC_ERR or VPECERR
    if ( DO_VPEC && !FOUND_VPECERR ) {
      if ( strcmp(c_get,"VPEC_ERR:") == 0 || strcmp(c_get,"VPECERR:")==0 ) {
	FOUND_VPECERR = true ;
	readdouble(fp, 1, &HOSTLIB.FIX_VPEC_ERR);
      }
    }

  } // end of while-fscanf 


  // -----------------------------

 VARCHECK:  // note there is NO more file-reading below


  //-----------------------
  // sanity check on optioanl SNPARams
  USE = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) ;
  FOUND_SNPAR = ( NVAR_STORE_SNPAR>0  ) ;

  if ( USE && FOUND_SNPAR == 0 ) {
    sprintf(c1err, "Found zero SN params for HOSTLIB_MSKOPT(%d) & %d option.",
	    INPUTS.HOSTLIB_MSKOPT, HOSTLIB_MSKOPT_USESNPAR );
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
      // xxxx if ( strcmp_ignoreCase(c_var,HOSTLIB.VARNAME_ALL[ivar])==0){ 
      if ( match_varname_HOSTLIB(c_var,HOSTLIB.VARNAME_ALL[ivar]) ) { 
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
    if ( !MATCH  ) {
      sprintf(c1err,"Could not find required HOSTLIB var '%s' (IVAR_STORE=%d)", 
	      c_var, IVAR_STORE );
      sprintf(c2err,"Check HOSTLIB VARNAMES, WGTMAP & HOSTLIB_STOREPAR key.");
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
  HOSTLIB.IVAR_TRUE_MATCH   = IVAR_HOSTLIB(HOSTLIB_VARNAME_TRUE_MATCH, 0) ; 
  HOSTLIB.IVAR_ZPHOT        = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZPHOT,   0) ; 
  HOSTLIB.IVAR_ZPHOT_ERR    = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZPHOT_ERR,0);
  HOSTLIB.IVAR_ZPHOT_Q0     = IVAR_HOSTLIB_PREFIX(HOSTLIB_PREFIX_ZPHOT_Q, 0);
  HOSTLIB.IVAR_VPEC         = IVAR_HOSTLIB(HOSTLIB_VARNAME_VPEC,    0) ; 
  HOSTLIB.IVAR_VPEC_ERR     = IVAR_HOSTLIB(HOSTLIB_VARNAME_VPEC_ERR,0);

  for (ivar=0;ivar<N_HOSTGAL_PROPERTY;ivar++){
    basename = HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar].BASENAME;
    
    sprintf(c_var, "%s_TRUE", basename);
    HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar].IVAR_TRUE = IVAR_HOSTLIB(c_var, 0);
    
    sprintf(c_var, "%s_OBS", basename);
    HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar].IVAR_OBS = IVAR_HOSTLIB(c_var, 0);
      
    sprintf(c_var, "%s_ERR", basename);
    HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar].IVAR_ERR = IVAR_HOSTLIB(c_var, 0);  
  }

  HOSTLIB.IVAR_ANGLE        = IVAR_HOSTLIB(HOSTLIB_VARNAME_ANGLE,0) ;   
  HOSTLIB.IVAR_FIELD        = IVAR_HOSTLIB(HOSTLIB_VARNAME_FIELD,0) ;   
  HOSTLIB.IVAR_ELLIPTICITY  = IVAR_HOSTLIB(HOSTLIB_VARNAME_ELLIPTICITY,0) ;
  HOSTLIB.IVAR_GALID2       = IVAR_HOSTLIB(HOSTLIB_VARNAME_GALID2,0) ;
  HOSTLIB.IVAR_SQRADIUS     = IVAR_HOSTLIB(HOSTLIB_VARNAME_SQRADIUS,0) ;
  HOSTLIB.IVAR_NBR_LIST     = IVAR_HOSTLIB(HOSTLIB_VARNAME_NBR_LIST,0) ; 
  HOSTLIB.IVAR_a_DLR        = IVAR_HOSTLIB(HOSTLIB_VARNAME_A_DLR, 0) ; 
  HOSTLIB.IVAR_b_DLR        = IVAR_HOSTLIB(HOSTLIB_VARNAME_B_DLR, 0) ; 

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

  // Mmake sure that these WGTMAP variables are really defined.
  // Also flag user-STOREPAR [EXTRA] variables that are also in WGTMAP (7.2019)
  int  NVAR_EXTRA  = HOSTLIB_OUTVAR_EXTRA.NOUT ;
  char *varName_WGTMAP, *varName_EXTRA;
  bool IS_SNVAR;
  for ( ivar_map=0;  ivar_map < NVAR_WGTMAP; ivar_map++ )  { 
    varName_WGTMAP = HOSTLIB_WGTMAP.VARNAME[ivar_map] ;
    IS_SNVAR       = HOSTLIB_WGTMAP.IS_SNVAR[ivar_map];

    // abort if varName_WGTMAP is not defined.
    if ( !IS_SNVAR )  { ivar = IVAR_HOSTLIB(varName_WGTMAP,1); }
  
    for(ivar=0; ivar < NVAR_EXTRA; ivar++ ) {
      varName_EXTRA = HOSTLIB_OUTVAR_EXTRA.NAME[ivar];
      if ( strcmp(varName_EXTRA,varName_WGTMAP) == 0 ) 
	{  HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[ivar] = 1 ; }

    }
  }

  // --------------------------------------------------
  // Misc error checking (added May 2020)

  if ( DO_VPEC ) {
    if ( HOSTLIB.IVAR_VPEC < 0 ) {  
      sprintf(c1err,"Could not find required VPEC column.");
      sprintf(c2err,"Add VPEC column, or remove %d-bit of HOSTLIB_MSKOPT",
	      HOSTLIB_MSKOPT_USEVPEC ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // make sure VPEC_ERR is defined in table column, or in header
    if ( HOSTLIB.IVAR_VPEC_ERR < 0 && HOSTLIB.FIX_VPEC_ERR < 0.0 ) {
      sprintf(c1err,"Found VPEC column, but VPEC_ERR not defined.");
      sprintf(c2err,"Add VPEC_ERR column, or define VPEC_ERR in header.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  if ( DO_RADEC && (HOSTLIB.IVAR_RA<0 || HOSTLIB.IVAR_DEC < 0) ) {
    print_preAbort_banner(fnam);
    printf("\t Valid RA  column names:  %s  %s  %s \n",
	   HOSTLIB_VARNAME_RA, HOSTLIB_VARNAME_RA_HOST, 
	   HOSTLIB_VARNAME_RA_GAL);
    printf("\t Valid DEC column names:  %s  %s  %s \n",
	   HOSTLIB_VARNAME_RA, HOSTLIB_VARNAME_RA_HOST, 
	   HOSTLIB_VARNAME_RA_GAL);
    sprintf(c1err,"Couldn't find required RA or DEC column (see above).");
    sprintf(c2err,"Add required column(s), or remove "
	    "%d-bit from HOSTLIB_MSKOPT", HOSTLIB_MSKOPT_SN2GAL_RADEC ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }


  if ( DO_SWAPZ && (HOSTLIB.IVAR_ZPHOT<0 || HOSTLIB.IVAR_ZPHOT_ERR<0) ) {
    sprintf(c1err,"Couldn't find required ZPHOT[_ERR] column .");
    sprintf(c2err,"Add required column(s), or remove "
	    "%d-bit from HOSTLIB_MSKOPT", HOSTLIB_MSKOPT_SWAPZPHOT ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  return ;

} // end of read_head_HOSTLIB

// ==============================
bool match_varname_HOSTLIB(char *varName0, char *varName1) {

  // Created Jun 10 2021
  // Perform case-insensitive match ... unless _obs is part of
  // varname. For host mag, we do NOT want to match a_obs to A_obs.

  bool match ;
  char fnam[] = "match_varname_HOSTLIB" ;

  // ------------ BEGIN -------------

  if ( strstr(varName0,"_obs") != NULL ) 
    // case-sensitive match for mag obs
    { match = ( strcmp(varName0,varName1) == 0 ) ; }
  else
    // case-insensitive match for all other variables
    { match = ( strcmp_ignoreCase(varName0,varName1) == 0 ) ; }

  return match;

} // end match_varname_HOSTLIB

// ==========================
bool checkSNvar_HOSTLIB_WGTMAP(char *varName) {

  // Created Mar 12 2020
  // Inputs:
  //   varName = variable name

  int N = HOSTLIB_WGTMAP.N_SNVAR;
  bool IS_SNVAR = false;

  // ---------------- BEGIN --------------

  if ( strcmp(varName,"x1")      == 0 ||
       strcmp(varName,"SALT2x1") == 0 ||
       strcmp(varName,"DM15")    == 0 ||
       strcmp(varName,"DELTA")   == 0  ) {
    IS_SNVAR = true ; 
    HOSTLIB_WGTMAP.ptrVal_SNVAR[N] = GENLC.ptr_SHAPEPAR;
  }

  if ( strcmp(varName,"c") == 0  || strcmp(varName,"SALT2c")==0 ) { 
    IS_SNVAR = true ; 
    HOSTLIB_WGTMAP.ptrVal_SNVAR[N] = &GENLC.SALT2c ;
  }

  if ( strcmp(varName,"AV") == 0 ) { 
    IS_SNVAR = true ; 
    HOSTLIB_WGTMAP.ptrVal_SNVAR[N] = &GENLC.AV ;
  }

  if ( strcmp(varName,"RV") == 0 ) { 
    IS_SNVAR = true ; 
    HOSTLIB_WGTMAP.ptrVal_SNVAR[N] = &GENLC.RV ;
  }

  return(IS_SNVAR);

} // end checkSNvarNames_HOSTLIB(int OPT, char *varName) {
 
// =====================================
void  checkAlternateVarNames_HOSTLIB(char *varName) {

  // Feb 12 2014
  // If input varName matches an allowed [hard-wired] alternative,
  // reset varName to the official name.
  // Note that the input argument is modified !

  char *BASENAME;
  int j;
  char fnam[] = "checkAlternateVarNames_HOSTLIB";

  // --------- BEGIN ---------

  if ( strcmp(varName,"ZERR") == 0 ) 
    { sprintf(varName,"%s", HOSTLIB_VARNAME_ZPHOT_ERR); }

  if ( strcmp(varName,"ZPHOTERR") == 0 ) 
    { sprintf(varName,"%s", HOSTLIB_VARNAME_ZPHOT_ERR); }

  if ( strcmp(varName,"VPECERR") == 0 ) 
    { sprintf(varName,"%s", HOSTLIB_VARNAME_VPEC_ERR); }

  for (j=0; j<N_HOSTGAL_PROPERTY; j++){
    BASENAME= HOSTLIB.HOSTGAL_PROPERTY_IVAR[j].BASENAME;
    // printf ("xxxx 1  %s varName=%s \n",fnam, varName);
    if ( strcmp(varName, BASENAME) == 0 ) { sprintf(varName,"%s_TRUE", BASENAME); }
    // printf ("xxxx 2  %s varName=%s \n",fnam, varName);
  }

  /* xxx Mark delete on Febr 2022
     else{
    if ( strcmp(varName,"LOGMASS") == 0 )  // legacy name (Jan 31 2020)
      { sprintf(varName,"%s", HOSTLIB_VARNAME_LOGMASS_TRUE); }
      } 
  xxx */
  
  if ( strcmp(varName,"REDSHIFT") == 0 )  // allowed in GENPDF_FILE (6/2020)
    { sprintf(varName,"%s", HOSTLIB_VARNAME_ZTRUE); }

} // end of   checkAlternateVarNames_HOSTLIB

// ====================================
void  parse_Sersic_n_fixed(FILE *fp, char  *string) {

  //
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
  // Feb 25 2020: set VALMIN & VALMAX for float; skip for ISCHAR.
  // Jan 22 2021: print WARNING if HOSTLIB.NSTAR > 0
  // Apr 30 2021: abort on NaN.

  bool DO_SWAPZPHOT = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SWAPZPHOT) ;
  int  IVAR_ZPHOT    = HOSTLIB.IVAR_ZPHOT ; // ivar_STORE
  int  IVAR_ZTRUE    = HOSTLIB.IVAR_ZTRUE ;
  int  IVAR_FIELD    = HOSTLIB.IVAR_FIELD ;
  int  IVAR_NBR_LIST = HOSTLIB.IVAR_NBR_LIST ;
  char c_get[200], FIELD[MXCHAR_FIELDNAME], NBR_LIST[MXCHAR_NBR_LIST] ;
  char fnam[] = "read_gal_HOSTLIB"  ;
  
  long long GALID, GALID_MIN, GALID_MAX ;
  int  ivar_ALL, ivar_STORE, NVAR_STORE, NGAL, NGAL_READ, MEMC ;
  int  NPRIORITY ;
  bool ISCHAR ;
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

      malloc_HOSTLIB(HOSTLIB.NGAL_STORE,HOSTLIB.NGAL_READ);

      NGAL_READ = HOSTLIB.NGAL_READ; // C-like index
      HOSTLIB.NGAL_READ++ ;          // fortran-like index

      read_galRow_HOSTLIB(fp, HOSTLIB.NVAR_ALL, xval, FIELD, NBR_LIST ); 

      if ( (HOSTLIB.NGAL_READ % 400000)==0  ) {
	printf("\t\t read %6d GAL rows\n", HOSTLIB.NGAL_READ);
	fflush(stdout);
      }

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

      // strip off variables to store
      for ( ivar_STORE=0; ivar_STORE < NVAR_STORE; ivar_STORE++ ) {
	ivar_ALL = HOSTLIB.IVAR_ALL[ivar_STORE] ;  // column in file
	val      = xval[ivar_ALL] ;

	HOSTLIB.VALUE_UNSORTED[ivar_STORE][NGAL] = val ;
	ISCHAR = ISCHAR_HOSTLIB(ivar_STORE);

	// keep track of min and max for each variable
	if ( ISCHAR == false ) {
	  if ( val > HOSTLIB.VALMAX[ivar_STORE] ) 
	    { HOSTLIB.VALMAX[ivar_STORE] = val; }
	  if ( val < HOSTLIB.VALMIN[ivar_STORE] ) 
	    { HOSTLIB.VALMIN[ivar_STORE] = val; }
	}
      }

      // store optional FIELD string 
      if ( IVAR_FIELD > 0  ) {
	sprintf(HOSTLIB.FIELD_UNSORTED[NGAL],"%s", FIELD);
      }

      if ( NGAL < -3 ) {
	ivar_ALL    = HOSTLIB.IVAR_ALL[HOSTLIB.IVAR_GALID] ;
	GALID       = (long long)xval[ivar_ALL];
	printf(" xxx %s: NGAL=%2d: GALID=%lld  NBR_LIST = '%s' \n", 
	       fnam, NGAL, GALID, NBR_LIST); fflush(stdout);
      }

      // Nov 11 2019: store optional NBR_LIST string 
      if ( IVAR_NBR_LIST > 0 ) {
	MEMC       = (1+strlen(NBR_LIST)) * sizeof(char) ;
	HOSTLIB.NBR_UNSORTED[NGAL] = (char*) malloc(MEMC);
	sprintf(HOSTLIB.NBR_UNSORTED[NGAL],"%s", NBR_LIST);
      }

      // Feb 26 2020: check option to replace ZTRUE with ZPHOT
      if ( DO_SWAPZPHOT ) {
	HOSTLIB.VALUE_UNSORTED[IVAR_ZTRUE][NGAL] = 
	  HOSTLIB.VALUE_UNSORTED[IVAR_ZPHOT][NGAL] ;
      }

      // store NGAL index vs. absolute READ index (for HOSTNBR)
      HOSTLIB.LIBINDEX_READ[NGAL_READ] = NGAL ; // NGAL_READ starts at 0

    }
  } // end of while-fscanf


 DONE_RDGAL:

  // load min/max ZTRUE into separate variables
  ivar_STORE   = IVAR_HOSTLIB(HOSTLIB_VARNAME_ZTRUE,1) ;
  HOSTLIB.ZMIN = HOSTLIB.VALMIN[ivar_STORE]; // min zHELIO
  HOSTLIB.ZMAX = HOSTLIB.VALMAX[ivar_STORE]; // max zHELIO

  printf("\t Stored %d galaxies from HOSTLIB (from %d total). \n",
	 HOSTLIB.NGAL_STORE, HOSTLIB.NGAL_READ );

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

  // check warning for stars (i.e, z ~ 0)
  if ( HOSTLIB.NSTAR > 0 ) {
    printf("\n\t *** WARNING: %d entries might be stars (z<%.4f) **** \n\n",
	   HOSTLIB.NSTAR, ZMAX_STAR);
    fflush(stdout);
  }

  // abort if requested z-range exceeds range of HOSTLIB.
  // Note that GENRANGE_REDSHIFT is for zCMB, so subtract/add "dz_safe"
  // to convert GENRANGE_REDSHIFT(CMB) to helio frame.
  double dz_safe  = 0.002 ;
  double ZMIN_GEN = INPUTS.GENRANGE_REDSHIFT[0] - dz_safe; // convert to helio
  double ZMAX_GEN = INPUTS.GENRANGE_REDSHIFT[1] + dz_safe;

  if ( HOSTLIB.ZMIN > ZMIN_GEN || HOSTLIB.ZMAX < ZMAX_GEN ) {

    print_preAbort_banner(fnam);
    printf("  HOSTLIB z-HEL range: %.5f - %.5f \n", 
	    HOSTLIB.ZMIN, HOSTLIB.ZMAX );
    printf("  Gen     z-HEL range: %.5f - %.5f (GENRANGE_REDSHIFT +- %.4f)\n",
	   ZMIN_GEN, ZMAX_GEN, dz_safe );
    printf("  Gen     z-CMB range: %.5f - %.5f (GENRANGE_REDSHIFT)\n",
	   INPUTS.GENRANGE_REDSHIFT[0], INPUTS.GENRANGE_REDSHIFT[1] );

    sprintf(c1err,"HOSTLIB z-HEL range does not contain user "
	    "GENRANGE_REDSHIFT;" );
    sprintf(c2err,"See redshift ranges above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
       
  // abort on NAN
  int NERR_NAN = HOSTLIB.NERR_NAN ;
  if ( NERR_NAN > 0 ) {
    sprintf(c1err,"%d HOSTLIB values are NaN", NERR_NAN );
    sprintf(c2err,"Must fix/remove these NaN values.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  return ;

} // end of read_gal_HOSTLIB

// ====================================
int passCuts_HOSTLIB(double *xval ) {

  // Return 1 if cuts are satisfied; zero otherwise.
  int ivar_ALL, LRA ,LRA2;
  double ZTRUE, RA, RA2, DEC;
  char fnam[] = "passCuts_HOSTLIB" ;

  // ---------- BEGIN ---------

  // bail for any option to rewrite HOSTLIB with appended columns
  if ( INPUTS.HOSTLIB_USE == 2 ) { return(1); }

  // REDSHIFT
  ivar_ALL    = HOSTLIB.IVAR_ALL[HOSTLIB.IVAR_ZTRUE] ; 
  ZTRUE       = xval[ivar_ALL];
  if ( ZTRUE < ZMAX_STAR ) { HOSTLIB.NSTAR++; }      // diagnostic
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

  // Created 9/16/2015
  // If there is no FIELD key, then read NVAL double from fp.
  // If there is a FIELD key, then read doubles and string separately.
  //
  // Dec 29 2017: use fgets for faster read.
  // Mar 28 2019: use MXCHAR_LINE_HOSTLIB
  // Nov 11 2019: check for NBR_LIST (string), analogous to FIELD check
  // May 06 2021: 
  //   + count NaN values (NERR_NAN)
  //   + check user override of GALMAG (ABMAG_FORCE)
  //
  // May 16 2022; fix bug reading too-short char.
  //              Replace WDLIST with global TMPWORD_HOSTLIB.
  //

  int MXCHAR         = MXCHAR_LINE_HOSTLIB;
  double ABMAG_FORCE = INPUTS.HOSTLIB_ABMAG_FORCE;
  int MXWD = NVAL;
  int ival_FIELD, ival_NBR_LIST, ival, NWD=0, len, NCHAR ;
  // xxx char WDLIST[MXVAR_HOSTLIB][100], *ptrWDLIST[MXVAR_HOSTLIB] ;
  char sepKey[] = " " ;
  char tmpWORD[200], tmpLine[MXCHAR_LINE_HOSTLIB], *pos, *varName ;
  char fnam[] = "read_galRow_HOSTLIB" ;
  // ---------------- BEGIN -----------------

  // scoop up rest of line with fgets
  fgets(tmpLine, MXCHAR, fp);

  NCHAR = strlen(tmpLine);
  if ( NCHAR >= MXCHAR-5 ) {
    print_preAbort_banner(fnam);
    printf(" LINE = '%s' (LEN=%d) \n", tmpLine, NCHAR );
    sprintf(c1err,"LINE likely exceeds bound of %d", MXCHAR);
    sprintf(c2err,"Shorten HOSTLIB lines, or increase MXCHAR_LINE_HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // remove line feed
  if ( (pos=strchr(tmpLine,'\n') ) != NULL )  { *pos = '\0' ; }

  // split string into words
  //xx  for(ival=0; ival < NVAL; ival++ ) { ptrWDLIST[ival] = WDLIST[ival]; }
  //xx  splitString2(tmpLine, sepKey, MXWD, &NWD, ptrWDLIST);
  splitString2(tmpLine, sepKey, MXWD, &NWD, TMPWORD_HOSTLIB);

  // abort if too few columns, but allow extra columns
  // (e..g, comment or catenated files with extra columns)
  if ( NWD < NVAL ) {
    print_preAbort_banner(fnam);
    printf("\t NGAL_READ = %d \n",  HOSTLIB.NGAL_READ );
    printf("\t LINE = '%s %s %s  ... ' \n", 
	   TMPWORD_HOSTLIB[0], TMPWORD_HOSTLIB[1], TMPWORD_HOSTLIB[2] );
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
    VALUES[ival] = -9.0 ; 
    varName      = HOSTLIB.VARNAME_ALL[ival] ;

    if ( ival == ival_FIELD )  { 
      sprintf(tmpWORD, "%s", TMPWORD_HOSTLIB[ival] );
      len = strlen(tmpWORD);
      if ( len > MXCHAR_FIELDNAME ) {
	sprintf(c1err,"strlen(FIELD=%s) = %d exceeds storage array of %d",
		tmpWORD, len, MXCHAR_FIELDNAME);
	sprintf(c2err,"Check HOSTLIB or increase MXCHAR_FIELDNAME");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      sscanf(tmpWORD, "%s", FIELD ); 
    }

    else if ( ival == ival_NBR_LIST )  { 
      sprintf(tmpWORD, "%s", TMPWORD_HOSTLIB[ival] );
      len = strlen(tmpWORD);
      if ( len > MXCHAR_NBR_LIST ) {
	sprintf(c1err,"strlen(NBR_LIST=%s) = %d exceeds storage array of %d",
		tmpWORD, len, MXCHAR_NBR_LIST);
	sprintf(c2err,"Check HOSTLIB or increase MXCHAR_NBR_LIST");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      sscanf(tmpWORD, "%s", NBR_LIST ); 

      if ( strstr(NBR_LIST,".") != NULL ) {
	long long  GALID = (long long)VALUES[0];
	//	print_preAbort_banner(fnam);
	sprintf(c1err,"Invalid NBR_LIST='%s'", NBR_LIST);
	sprintf(c2err,"No decimals allowed; GALID=%lld", GALID);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }

    else {
      // float or int (non-string) value
      sscanf(TMPWORD_HOSTLIB[ival], "%le", &VALUES[ival] ); 

      // check option to force override for gal mags (May 2021)
      if ( ABMAG_FORCE > -8.0 ) {	
	int lenvar = strlen(varName);
	int lensuf = strlen(HOSTLIB_SUFFIX_MAGOBS);
	bool MATCHMAG = ( strstr(varName,HOSTLIB_SUFFIX_MAGOBS) != NULL );
	if ( MATCHMAG && lenvar == lensuf+1 ) 
	  { VALUES[ival] = ABMAG_FORCE ; }
      }

      // check for NaN (NaN abort is after reading entire HOSTLIB)
      if ( isnan(VALUES[ival]) )  {
	HOSTLIB.NERR_NAN++ ;
	if ( HOSTLIB.NERR_NAN < 20 ) 
	  { printf("\t ERROR: HOSTLIB %s = NaN \n", varName );	}
      }

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
  //
  // Jan 15 2020:  bug fix: check NGAL_READ==0, not NGAL_STORE==0.

  int ivar, I8, I8p, I4, ICp, MEMC, DO_FIELD, DO_NBR, igal ;
  int LDMP = 0 ;
  //  char fnam[] = "malloc_HOSTLIB";

  // ------------- BEGIN ----------

  I8   = sizeof(double);
  I8p  = sizeof(double*);
  I4   = sizeof(int);
  ICp  = sizeof(char*);

  DO_FIELD    = ( HOSTLIB.IVAR_FIELD    > 0 );
  DO_NBR      = ( HOSTLIB.IVAR_NBR_LIST > 0 );

  if ( NGAL_READ == 0 ) {

    if  ( LDMP ) 
      { printf(" xxx Initial malloc for HOSTLIB ... \n"); fflush(stdout); }

    HOSTLIB.MALLOCSIZE_D  += (I8  * MALLOCSIZE_HOSTLIB) ;
    HOSTLIB.MALLOCSIZE_I  += (I4  * MALLOCSIZE_HOSTLIB) ;
    HOSTLIB.MALLOCSIZE_Cp += (ICp * MALLOCSIZE_HOSTLIB) ;

    HOSTLIB.NGAL_STORE_MALLOC = 0;

    // allocate memory for each variable to store
    for ( ivar = 0; ivar < HOSTLIB.NVAR_STORE; ivar++ ) {
      HOSTLIB.VALUE_UNSORTED[ivar] = (double*)malloc(HOSTLIB.MALLOCSIZE_D);
    }

    HOSTLIB.LIBINDEX_READ = (int *)malloc(HOSTLIB.MALLOCSIZE_I);
    for(igal=0; igal < MALLOCSIZE_HOSTLIB; igal++ ) 
      { HOSTLIB.LIBINDEX_READ[igal] = -9; }

    if ( DO_FIELD ) {
      HOSTLIB.FIELD_UNSORTED = (char**)malloc( HOSTLIB.MALLOCSIZE_Cp );
      MEMC = MXCHAR_FIELDNAME *sizeof(char) ;
      for(igal=0; igal < MALLOCSIZE_HOSTLIB; igal++ )
	{ HOSTLIB.FIELD_UNSORTED[igal] = (char*)malloc( MEMC );  }  
    }

    if ( DO_NBR ) {
      HOSTLIB.NBR_UNSORTED = (char**)malloc(HOSTLIB.MALLOCSIZE_Cp);
      // do NOT malloc string length here; malloc later
      // when length of string is known.

    }

    return ;
  }


  // --------------------------------------------------
  int UPD = ( (NGAL_STORE % MALLOCSIZE_HOSTLIB) == 0 );
  // xxx mark delete  if ( UPD && NGAL_STORE > 0 ) {    
  if ( UPD && NGAL_STORE > HOSTLIB.NGAL_STORE_MALLOC ) {    
    HOSTLIB.MALLOCSIZE_D  += (I8  * MALLOCSIZE_HOSTLIB) ;
    HOSTLIB.MALLOCSIZE_Cp += (ICp * MALLOCSIZE_HOSTLIB) ;
    HOSTLIB.NGAL_STORE_MALLOC = NGAL_STORE; // May 2022

    if  ( LDMP )  { 
      printf("\t xxx Extend HOSTLIB malloc at NGAL_STORE=%d \n", 
	     NGAL_STORE);    fflush(stdout); 
    }

    for ( ivar = 0; ivar < HOSTLIB.NVAR_STORE; ivar++ ) {
      HOSTLIB.VALUE_UNSORTED[ivar] =
	(double*)realloc(HOSTLIB.VALUE_UNSORTED[ivar], HOSTLIB.MALLOCSIZE_D);
    }

    if ( DO_FIELD ) {
      HOSTLIB.FIELD_UNSORTED = 
	(char**)realloc( HOSTLIB.FIELD_UNSORTED, HOSTLIB.MALLOCSIZE_Cp );
      MEMC = MXCHAR_FIELDNAME *sizeof(char) ;
      for(igal=NGAL_STORE; igal < NGAL_STORE+MALLOCSIZE_HOSTLIB; igal++ ) 
	{ HOSTLIB.FIELD_UNSORTED[igal] = (char*)malloc( MEMC ); }
    }

    if ( DO_NBR ) {
      HOSTLIB.NBR_UNSORTED = 
	(char**)realloc( HOSTLIB.NBR_UNSORTED, HOSTLIB.MALLOCSIZE_Cp );
      // do NOT malloc size of each string
    }

  } // end if block

  // - - - - - - - - - - -
  // separate check for READ index
  if ( (NGAL_READ % MALLOCSIZE_HOSTLIB) == 0 ) {
    HOSTLIB.MALLOCSIZE_I  += (I4  * MALLOCSIZE_HOSTLIB) ;
    HOSTLIB.LIBINDEX_READ = 
      (int *)realloc(HOSTLIB.LIBINDEX_READ, HOSTLIB.MALLOCSIZE_I);
    for(igal=NGAL_READ; igal < NGAL_READ+MALLOCSIZE_HOSTLIB; igal++ ) 
      { HOSTLIB.LIBINDEX_READ[igal] = -9; }
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
    if (INPUTS.HOSTLIB_GALID_UNIQUE){
      printf("\t Assign GALID_UNIQUE for re-used hosts and "
	     "randomize HOSTMAG.\n");
    }
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
  // Use CERNLIB sortz function to sort library by redshift. 
  // Note that sortzv accepts a real*4 array, so be careful.
  //
  // Also compute ZGAPMAX, ZGAPAVG and Z_ATGAPMAZ
  // Also compute VPEC stuff to avoid another igal loop.
  //
  // Nov 11 2019: check NBR_LIST
  // May 23 2020: compute a few VPEC quantities for README

  bool DO_VPEC  = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEVPEC ) ;
  int  NGAL, igal, ival, unsort, VBOSE, DO_FIELD, DO_NBR;
  int  IVAR_ZTRUE, NVAR_STORE, ORDER_SORT, MEMC, IVAR_VPEC ;
  double ZTRUE, ZLAST, ZGAP, ZSUM, *ZSORT, VAL ;
  double VPEC, VSUM, VSUMSQ ;
  long long GALID;
  char *ptr_UNSORT ;
  char fnam[] = "sortz_HOSTLIB" ;

  // ------------- BEGIN -------------

  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );

  if ( VBOSE )  { 
    printf("\t Sort HOSTLIB by redshift (%.4f to %.4f) \n",
	   HOSTLIB.ZMIN , HOSTLIB.ZMAX );   
    
    printf("\t |zSN-zGAL| tolerance zpoly: %s \n",
	   INPUTS.HOSTLIB_GENPOLY_DZTOL.STRING );
    
    fflush(stdout); 
  }


  DO_FIELD = ( HOSTLIB.IVAR_FIELD    > 0 ) ;
  DO_NBR   = ( HOSTLIB.IVAR_NBR_LIST > 0 ) ;

  NGAL = HOSTLIB.NGAL_STORE ;
  NVAR_STORE = HOSTLIB.NVAR_STORE ;
  IVAR_ZTRUE = HOSTLIB.IVAR_ZTRUE ;
  IVAR_VPEC  = HOSTLIB.IVAR_VPEC ;

  // allocate memory for sort-pointers
  HOSTLIB.LIBINDEX_UNSORT  = (int*)malloc( (NGAL+1) * sizeof(int) );
  HOSTLIB.LIBINDEX_ZSORT   = (int*)malloc( (NGAL+1) * sizeof(int) );

  // allocate memory for sorted values
  for ( ival=0; ival < NVAR_STORE; ival++ ) {
    HOSTLIB.VALUE_ZSORTED[ival] = 
      (double*)malloc( (NGAL+1) * sizeof(double) ) ;
  }

  if ( DO_FIELD  ) {
    HOSTLIB.FIELD_ZSORTED = (char**)malloc( (NGAL+1) * sizeof(char*) ) ;
    MEMC = MXCHAR_FIELDNAME * sizeof(char) ;
    for ( igal=0; igal <= NGAL; igal++ ) 
      { HOSTLIB.FIELD_ZSORTED[igal] = (char*)malloc(MEMC) ; }
  }

  if ( DO_NBR ) {
    HOSTLIB.NBR_ZSORTED = (char**)malloc( (NGAL+1) * sizeof(char*) ) ;
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

  HOSTLIB.VPEC_RMS = HOSTLIB.VPEC_AVG = 0.0 ;
  HOSTLIB.VPEC_MIN = HOSTLIB.VPEC_MAX = 0.0 ;
  VSUM = VSUMSQ = 0.0;

  // fill sorted array. 'igal' is the z-sorted index; 
  // 'unsort' is the  un-sorted index matching the original HOSTLIB order.
  for ( igal = 0; igal < NGAL ; igal++ ) {
    unsort = HOSTLIB.LIBINDEX_UNSORT[igal]  ;
    HOSTLIB.LIBINDEX_ZSORT[unsort] = igal;

    for ( ival=0; ival < NVAR_STORE; ival++ ) {
      VAL = HOSTLIB.VALUE_UNSORTED[ival][unsort] ; 
      HOSTLIB.VALUE_ZSORTED[ival][igal] = VAL;
    }

    if ( DO_FIELD ) {
      ptr_UNSORT = HOSTLIB.FIELD_UNSORTED[unsort];
      sprintf(HOSTLIB.FIELD_ZSORTED[igal],"%s", ptr_UNSORT);
      free(HOSTLIB.FIELD_UNSORTED[unsort]); // Nov 11 2019
    }
    
    if ( DO_NBR ) {  // Nov 11 2019 
      ptr_UNSORT = HOSTLIB.NBR_UNSORTED[unsort];
      MEMC = (1+strlen(ptr_UNSORT)) * sizeof(char);
      if ( MEMC == 0 ) { MEMC = 4 ; }
      HOSTLIB.NBR_ZSORTED[igal] = (char*)malloc(MEMC) ; 
      sprintf(HOSTLIB.NBR_ZSORTED[igal], "%s", ptr_UNSORT);
      free(HOSTLIB.NBR_UNSORTED[unsort]); 
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

    // store IGAL if GALID_FORCE is set (Dec 29 2021)
    if ( INPUTS.HOSTLIB_GALID_FORCE > 0 ) {
      GALID  = get_GALID_HOSTLIB(igal);
      if ( INPUTS.HOSTLIB_GALID_FORCE == GALID ) 
	{ HOSTLIB.IGAL_FORCE = igal; }
    }

    if ( DO_VPEC ) {
      VPEC = HOSTLIB.VALUE_ZSORTED[IVAR_VPEC][igal];
      if ( VPEC > HOSTLIB.VPEC_MAX ) { HOSTLIB.VPEC_MAX = VPEC; }
      if ( VPEC < HOSTLIB.VPEC_MIN ) { HOSTLIB.VPEC_MIN = VPEC; }
      VSUM += VPEC;  VSUMSQ += (VPEC*VPEC);

      // Feb 24 2021: abort on ZTRUE+VPEC too small
      double zcheck     = ZTRUE + (VPEC/LIGHT_km);
      double ZTRUE_MIN  = INPUTS.GENRANGE_REDSHIFT[0] ;
      if ( zcheck < ZMIN_HOSTLIB && ZTRUE >= ZTRUE_MIN ) {
	print_preAbort_banner(fnam);
	printf("\t ZTRUE(HOSTLIB) = %f \n", ZTRUE);
	printf("\t VPEC(HOSTLIB)  = %f \n", VPEC);
	printf("\t ZMIN_HOSTLIB   = %f \n", ZMIN_HOSTLIB);
	sprintf(c1err,"ZTRUE + VPEC/c = %f  < %f", zcheck, ZMIN_HOSTLIB);
	sprintf(c2err,"Fix HOSTLIB VPEC or increase GENRANGE_REDSHIFT[0] > %f",
		INPUTS.GENRANGE_REDSHIFT[0]);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }

  } // end of igal loop

  // - - - - - - - 

  double d_NGAL = (double)NGAL;

  HOSTLIB.ZGAPAVG = ZSUM/d_NGAL;

  // VPEC stuff
  if ( DO_VPEC  ) {
    HOSTLIB.VPEC_AVG = VSUM/d_NGAL;   
    HOSTLIB.VPEC_RMS = STD_from_SUMS(NGAL, VSUM, VSUMSQ);
  }

  // free memory for the pointers and the unsorted array.
  free(ZSORT);
  for ( ival=0; ival < NVAR_STORE; ival++ ) 
    { free(HOSTLIB.VALUE_UNSORTED[ival]);  }

  int  OPT_PLUSMAGS  = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_PLUSMAGS);
  int  OPT_PLUSNBR   = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_PLUSNBR);
  if ( !(OPT_PLUSMAGS || OPT_PLUSNBR) ) {
    free(HOSTLIB.LIBINDEX_UNSORT);
  }

} // end of sortz_HOSTLIB


// =======================================
void zptr_HOSTLIB(void) {

  // construct z-pointers to z-sorted library so that later we can 
  // quickly find the nearest redshift in library.
  //
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
void init_HOSTLIB_WGTMAP(int OPT_INIT, int IGAL_START, int IGAL_END) {

  // Wgt map is already read and redshift-sorted.
  // Do some inits, allocate memory and make sanity checks.
  //
  // Inputs:
  //   OPT_INIT   : not used
  //   IGAL_START : first IGAL to compute wgt
  //   IGAL_END   : last IGAL to compute wgt
  //
  // - - - - - - - - - - - - - - - 
  // Jan 27 2017: fix bug mallocing GRIDMAP_HOSTLIB_WGT.FUNVAL;
  //              I8p -> I8p*2
  // 
  // Jun 18 2019: if interp_GRIDMAP fails, print more PRE-ABORT info.
  // Jun 25 2019: check GAMMA_GRID option
  // May 03 2022: check new HOSTLIB feature for TRUE column (IVAR_TRUE_MATCH)

  bool IS_SNVAR ;
  int  i, NDIM, ivar, ivar_STORE, NFUN, NROW, ibin, istat ;
  int  NGAL, igal, isparse, IVAL ;
  short int I2MAG;
  bool VBOSE, LDMPWGT ;

  int N_SNVAR         = HOSTLIB_WGTMAP.N_SNVAR ;
  int NBTOT_SNVAR     = HOSTLIB_WGTMAP.NBTOT_SNVAR ;
  int IVAR_TRUE_MATCH = HOSTLIB.IVAR_TRUE_MATCH;

  double GAMMA_GRID_MIN = INPUTS.BIASCOR_SALT2GAMMA_GRID[0]; 
  double GAMMA_GRID_MAX = INPUTS.BIASCOR_SALT2GAMMA_GRID[1]; 
  bool   USE_GAMMA_GRID = (GAMMA_GRID_MAX > GAMMA_GRID_MIN );  

  long long GALID, GALID_CHECK ;

  double
    VAL, VALMIN, VALMAX, WGT, WGTSUM, WGTSUM_LAST, SNMAGSHIFT
    ,VAL_WGTMAP[MXVAR_HOSTLIB], VAL_SNVAR[MXVAR_HOSTLIB]
    ,ZTRUE, ZTRUE_CHECK, ZDIF, TMPVAL[2]
    ;

  char cvar[40], *varName ;
  char fnam[] = "init_HOSTLIB_WGTMAP" ;

  // --------- BEGIN -----------

  VBOSE = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_VERBOSE );

  malloc_HOSTLIB_WGTMAP();

  printf("\t Interpolate WGTMAP for each galaxy (%.2f MB) ... \n",
	 HOSTLIB_WGTMAP.MEMTOT_MB ) ;  fflush(stdout);

  NDIM = HOSTLIB_WGTMAP.GRIDMAP.NDIM ;
  NFUN = HOSTLIB_WGTMAP.GRIDMAP.NFUN ;
  NROW = HOSTLIB_WGTMAP.GRIDMAP.NROW ;
  NGAL = HOSTLIB.NGAL_STORE ;
  HOSTLIB_WGTMAP.USE_SALT2GAMMA_GRID = USE_GAMMA_GRID;

  if ( USE_GAMMA_GRID  ) {
    printf("\t Implement BIASCOR_SALT2GAMMA_GRID: %.2f to %.2f mag\n",
	   GAMMA_GRID_MIN, GAMMA_GRID_MAX);
    fflush(stdout);
  }

  for(ibin=0; ibin < NBTOT_SNVAR ; ibin++ ) {

    if ( N_SNVAR > 0 )  // fetch SN grid values
      { getVal_SNVAR_HOSTLIB_WGTMAP(ibin,VAL_SNVAR); }

    WGTSUM_LAST = 0.0 ;
    for ( igal=IGAL_START ; igal <= IGAL_END ; igal++ ) {
	
      GALID  = get_GALID_HOSTLIB(igal);
      ZTRUE  = get_ZTRUE_HOSTLIB(igal);
      
      if ( NROW == 0 ) 
	{  WGT = 1.0 ;  SNMAGSHIFT = 0.0 ;  goto WGTSUM ;    }
      
      // May 2022: 
      // if TRUE column exists in HOSTLIB, set WGT=0.0 for TRUE_MATCH=0
      // so that these excluded hosts are used only for DLR-matching
      // and not for true-host matches.
      if ( IVAR_TRUE_MATCH > 0 ) { 
	IVAL  = get_VALUE_HOSTLIB(IVAR_TRUE_MATCH,igal); 
	if ( IVAL == 0 ) { WGT = SNMAGSHIFT = 0.0; goto WGTSUM; }
      }

      // strip off variables used for weighting
      for ( ivar=0; ivar < NDIM; ivar++ ) {  // WGTMAP variables
	IS_SNVAR     = HOSTLIB_WGTMAP.IS_SNVAR[ivar]; 
	isparse      = -9 ;
	if ( !IS_SNVAR ) {
	  // get VAL from HOSTLIB
	  ivar_STORE   = HOSTLIB.IVAR_STORE[ivar];
	  VAL          = HOSTLIB.VALUE_ZSORTED[ivar_STORE][igal] ;
	}
	else {
	  // get VAL from SN property
	  isparse =  HOSTLIB_WGTMAP.INVSPARSE_SNVAR[ivar] ;
	  VAL     =  VAL_SNVAR[isparse];
	}
	VAL_WGTMAP[ivar] = VAL ;

      } // end ivar loop
      

      // interpolate to get TMPVAL = WGT and SNMAGSIFT
      istat = interp_GRIDMAP(&HOSTLIB_WGTMAP.GRIDMAP, VAL_WGTMAP, TMPVAL ) ;
      
      if ( istat != SUCCESS ) {
	print_preAbort_banner(fnam);
	printf("\t GALID = %lld  (ibin_SNVAR=%d, igal=%d)\n", 
	       GALID, ibin, igal );
	for ( ivar=0; ivar < NDIM; ivar++ ) {  // WGTMAP variables
	  IS_SNVAR     = HOSTLIB_WGTMAP.IS_SNVAR[ivar]; 
	  if ( IS_SNVAR ) { continue; }
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
      
      // check for random assignment of SNMAGSHIFT (for BiasCor)
      if(USE_GAMMA_GRID) 
	{ SNMAGSHIFT = snmagshift_salt2gamma_HOSTLIB(GALID); }

      // convert mag shift to 2-byte int to reduce memory
      I2MAG = (short int)(SNMAGSHIFT*I2MAGSCALE_HOSTLIB);

      // local sum
      WGTSUM = WGTSUM_LAST + WGT;
      
      // print first 3 weights and last wgt
      if ( VBOSE && ibin==0 && (igal <= 1 || igal == NGAL-1) ) {
	printf("\t   WGT(igal=%d,GALID=%lld) = %f -> WGTSUM = %f \n", 
	       igal, GALID, WGT, WGTSUM ); 
	fflush(stdout);
      }
	
      // load global array for each sum            
      if ( N_SNVAR > 0 ) {
	HOSTLIB_WGTMAP.WGTSUM_SNVAR[ibin][igal]     = WGTSUM ;
	HOSTLIB_WGTMAP.I2SNMAGSHIFT_SNVAR[ibin][igal] = I2MAG ;
      }
      else {
	HOSTLIB_WGTMAP.WGTSUM[igal]     = WGTSUM ;
	HOSTLIB_WGTMAP.I2SNMAGSHIFT[igal] = I2MAG ;
      }

      WGTSUM_LAST  =  WGTSUM ;
	
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

    }   // end if igal loop
  }   // end NBTOT


  // --------------------------
  // verify interpolated WGTMAP values against optional list of 
  // exact WGT values specified by the WGTMAP_CHECK keys.

  runCheck_HOSTLIB_WGTMAP();
  
  return ;

} // end of init_HOSTLIB_WGTMAP


// =========================================
void malloc_HOSTLIB_WGTMAP(void) {

  // Mar 15 2020
  // allocate memory (for each host gal) for wgt, snmagshift and 
  // cumulative weight-sum.

  int  NGAL        = HOSTLIB.NGAL_STORE ;
  int  N_SNVAR     = HOSTLIB_WGTMAP.N_SNVAR ;
  int  NBTOT_SNVAR = HOSTLIB_WGTMAP.NBTOT_SNVAR ;
  int  MEMDD       = NBTOT_SNVAR * sizeof(double*);
  int  MEMD2       = NGAL * sizeof(double);
  int  MEMSS       = NBTOT_SNVAR * sizeof(short int*);
  int  MEMS2       = NGAL * sizeof(short int);

  int  ibin ;
  double MEMTOT = 0.0 ;
  //  char fnam[] = "malloc_HOSTLIB_WGTMAP" ;

  // --------------- BEGIN -----------

  
  if ( N_SNVAR > 0 ) { 

    HOSTLIB_WGTMAP.WGTSUM_SNVAR       = (double**) malloc(MEMDD);
    HOSTLIB_WGTMAP.I2SNMAGSHIFT_SNVAR = (short int **) malloc(MEMSS);

    for(ibin=0; ibin < NBTOT_SNVAR; ibin++ ) {
      HOSTLIB_WGTMAP.WGTSUM_SNVAR[ibin]       = (double*) malloc(MEMD2);
      HOSTLIB_WGTMAP.I2SNMAGSHIFT_SNVAR[ibin] = (short int*) malloc(MEMS2);
      MEMTOT += (double)(MEMD2 + MEMS2);
    }
  }
  else {
    // HOSTLIB vars only
    HOSTLIB_WGTMAP.WGTSUM       = (double *)malloc(MEMD2);
    HOSTLIB_WGTMAP.I2SNMAGSHIFT = (short int *)malloc(MEMS2);
    MEMTOT += (double)(MEMD2 + MEMS2) ;
  }

  HOSTLIB_WGTMAP.MEMTOT_MB = ( MEMTOT*1.0E-6 );

  return;

} // end malloc_HOSTLIB_WGTMAP


// =========================================
void runCheck_HOSTLIB_WGTMAP(void) {

  // Created Mar 11 2020   [moved code out of init_HOSTLIB_WGTMAP]
  // verify interpolated WGTMAP values against optional list of 
  // exact WGT values specified by the WGTMAP_CHECK keys.

  int NCHECK, igal_difmax, i, igal, NN ;
  long long GALID; 
  double WDIF, WDIF_SUM, WDIF_MAX, WDIF_RMS, WDIF_AVG, SQWDIF_SUM ;
  double XN, SQTMP, WGT_EXACT, WGT_INTERP ;
  char fnam[] = "runCheck_HOSTLIB_WGTMAP" ;

  // ----------- BEGIN ------------

  NCHECK = HOSTLIB_WGTMAP.NCHECKLIST; 
  WDIF_SUM = SQWDIF_SUM = WDIF_MAX = 0.0 ;
  igal_difmax = 0;
  NN = 0 ;

  for ( i=0; i <  NCHECK ; i++ ) {
    igal        = HOSTLIB_WGTMAP.CHECKLIST_IGAL[i] ;
    GALID       = HOSTLIB_WGTMAP.CHECKLIST_GALID[i] ;
    WGT_EXACT   = HOSTLIB_WGTMAP.CHECKLIST_WGT[i] ;
    WGT_INTERP  = HOSTLIB_WGTMAP.WGTSUM[igal] - HOSTLIB_WGTMAP.WGTSUM[igal-1];
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

  return ;

} // end  runCheck_HOSTLIB_WGTMAP 


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
void init_HOSTLIB_ZPHOT_QUANTILE(void) {

  // Created Apr 19 2022 by R.Kess;er
  // For nominal usage of zphot quantiles in HOSTLIB, this function
  // only prints information to stdout.
  // If force-Gauss option is set, this function initializes
  // Gaussian integrals at percentile values.

  int  N_Q        = HOSTLIB.NZPHOT_Q ;
  int  USE_QGAUSS = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_ZPHOT_QGAUSS );
  int  q, qbin, percentile;
  char fnam[] = "init_HOSTLIB_ZPHOT_QUANTILE" ;

  // ------------ BEGIN ------------

  // if Gaussian quantile is requested without any quantiles in the
  // HOSTLIB, force 10 quantiles.
  if ( USE_QGAUSS && N_Q == 0 ) {
    N_Q = HOSTLIB.NZPHOT_Q = 11; // 0, 10, 20 ..... 100
    qbin = 100/(N_Q - 1);
    for(q=0; q < N_Q; q++ ) {
      percentile = q * qbin;
      HOSTLIB.PERCENTILE_ZPHOT_Q[q] = percentile ;
      sprintf(HOSTLIB.VARNAME_ZPHOT_Q[q],"%s%3.3d", 
	      PREFIX_ZPHOT_Q, percentile);
    }
  }

  if ( N_Q == 0 ) { return; }

  // check option to force Gaussian quantiles using Gaussian defined by
  // mean=ZPHOT and sigma=ZPHOT_ERR. Here define quantile values for 
  // unit Gaussian

  if ( USE_QGAUSS ) {
    printf("\t Force %d Gaussian quantiles\n", N_Q ); fflush(stdout);

    int IVAR_ZPHOT     = HOSTLIB.IVAR_ZPHOT;
    int IVAR_ZPHOT_ERR = HOSTLIB.IVAR_ZPHOT_ERR;
    if ( IVAR_ZPHOT < 0 || IVAR_ZPHOT_ERR < 0 ) {
      sprintf(c1err,"Cannot compute forced Gaussian quantiles with");
      sprintf(c2err,"IVAR_ZPHOT=%d and IVAR_ZPHOT_ERR=%d",
	      IVAR_ZPHOT, IVAR_ZPHOT_ERR);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
    }

    double *xsig_store = (double*)malloc( sizeof(double) * N_Q ) ;
    double *xdif_store = (double*)malloc( sizeof(double) * N_Q ) ;
    double xsigmin = -6.0, xsigmax=6.0, xsigstep=0.001, xsig, xdif ;
    double gint, gint_target;

    for (q=0; q < N_Q; q++ )  { xdif_store[q] = xsig_store[q] = 999.0; }

    for(xsig=xsigmin; xsig < xsigmax; xsig += xsigstep) {
      gint = GaussIntegral(xsigmin,xsig);

      //  printf(" xxx %s: xsig = %8.5f, gint = %8.5f \n", fnam, xsig, gint);
      for ( q=0; q < N_Q; q++ ) {
	gint_target =  (double)HOSTLIB.PERCENTILE_ZPHOT_Q[q]/100.0 ;
	if ( gint_target == 0.0 ) { gint_target = 0.001; }
	if ( gint_target == 1.0 ) { gint_target = 0.999; }

	xdif = fabs(gint - gint_target);
	if ( xdif < xdif_store[q] ) {
	  xdif_store[q] = xdif ;
	  xsig_store[q] = xsig ;
	}
      }
    } // end xsig loop

    // store Nsigma values in global struct.
    for ( q=0; q < N_Q; q++ ) { HOSTLIB.SIGMA_QGAUSS[q] = xsig_store[q]; }

    int LDMP_QGAUSS = 0 ;
    if ( LDMP_QGAUSS ) {
      for ( q=0; q < N_Q; q++ ) {
	percentile = HOSTLIB.PERCENTILE_ZPHOT_Q[q];
	xsig       = xsig_store[q];
	gint = GaussIntegral(xsigmin,xsig);
	printf(" xxxx Gauss percentile =%3d -> NSIGMA = %8.5f  "
	       "(exact GaussInt=%6.4f)\n", 
	       percentile, xsig, gint); fflush(stdout);
      }
      debugexit(fnam);
    }

  }  // end USE_QGAUSS


  // list zphot quantiles
  char STRING_Q[200], str_q[40];
  sprintf(STRING_Q,"ZPHOT_QUANTILES: ");
  for(q=0; q<N_Q; q++ ) {
    percentile = HOSTLIB.PERCENTILE_ZPHOT_Q[q] ;
    sprintf(str_q,"%d ", percentile);
    strcat(STRING_Q, str_q);
  }

  printf("\t %s\n", STRING_Q);
  fflush(stdout);

  return;

} // init_HOSTLIB_ZPHOT_QUANTILE

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

  int NR, j, jth, ifilt, ifilt_obs, IVAR, NMAGOBS, MATCH_FLAG, IVAR_ERR ;
  char cvar[12], cfilt[2], cvar_err[40];
  char fnam[] = "init_GALMAG_HOSTLIB" ;
  double Rmax, TH, THbin ; 

  // --------------- BEGIN --------------

  // always check for gal mags. Store IVAR for each observer-mag
  NMAGOBS = 0 ;
  MATCH_FLAG = 2; // case-sensitive varname match (Apr 29 2021)

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];

    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] ); 
    magkey_HOSTLIB(ifilt_obs,cvar,cvar_err); // returns cvar

    IVAR = IVAR_HOSTLIB(cvar,MATCH_FLAG) ;
    IVAR_ERR = IVAR_HOSTLIB(cvar_err,MATCH_FLAG) ;
    //printf("xxx cvar_err = %s, ivar_err = %d\n", cvar_err, IVAR_ERR);


    if ( IVAR > 0 ) {
      NMAGOBS++ ;
      HOSTLIB.IVAR_MAGOBS[ifilt_obs] = IVAR ;
      strcat(HOSTLIB.filterList,cfilt) ;
    }
    if ( IVAR_ERR > 0 ) {
      NMAGOBS++ ;
      HOSTLIB.IVAR_MAGOBS_ERR[ifilt_obs] = IVAR_ERR ;
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
  if ( SERSIC_PROFILE.NPROF <= 0 ) {
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
  NR++ ; HOSTLIB.Aperture_PSFSIG[NR] = PSFMAX_SNANA/2.35 ;

  if ( NR != NMAGPSF_HOSTLIB ) {
    sprintf(c1err,"%d defined PSF values", NR );
    sprintf(c2err,"but expected NMAGPSF_HOSTLIB = %d", 
	    NMAGPSF_HOSTLIB );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  
  if ( NMAGOBS <= 0 ) {
    sprintf(c1err, "%s" , "Galaxy mags under SN requested, but no mags are");
    sprintf(c2err, "given in HOSTLIB. Need [filt]%s keys.", 
	    HOSTLIB_SUFFIX_MAGOBS );
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
double interp_GALMAG_HOSTLIB(int ifilt_obs, double PSFSIG ) {
 
  //
  // Return interpolated GALMAG for this 
  // obs-filter index and PSF-sigma (arcsec).
  // Interpolate pre-computed GALMAG vs. PSF.
  // GALMAG corresponds to flux contained in 4*PI*PSF^2 aperture.
  //
  // Inputs:
  //  ifilt_obs = absolute obs-filter index
  //  PSFSIG    = sigma(PSF) in arcsec
  //
  // Mar 23 2021: if PSF is outside range of map, bring PSF withing
  //    range instead of abort.

  int NPSF ;
  double GALMAG, PSFSIGmin, PSFSIGmax, PSFSIG_local ;
  double *PTRGRID_GALMAG, *PTRGRID_PSF ;
  char fnam[] = "interp_GALMAG_HOSTLIB" ;

  // ------------- BEGIN -------------
  
  NPSF = NMAGPSF_HOSTLIB ;
  PSFSIGmin = HOSTLIB.Aperture_PSFSIG[1] ;
  PSFSIGmax = HOSTLIB.Aperture_PSFSIG[NPSF] ;

  PSFSIG_local = PSFSIG;
  if ( PSFSIG < PSFSIGmin ) { PSFSIG_local = PSFSIGmin + 0.0001; }
  if ( PSFSIG > PSFSIGmax ) { PSFSIG_local = PSFSIGmax - 0.0001; }

  // note that the zero'th element is total aperture,
  // so ignore it for interpolatio.
  PTRGRID_PSF     = &HOSTLIB.Aperture_PSFSIG[1] ;
  PTRGRID_GALMAG  = &SNHOSTGAL.GALMAG[ifilt_obs][1] ;

  GALMAG = interp_1DFUN (1, PSFSIG_local,
			 NPSF, PTRGRID_PSF, PTRGRID_GALMAG, "GALMAG");

  return(GALMAG) ;

}  // interp_GALMAG_HOSTLIB



// ======================================
void magkey_HOSTLIB(int ifilt_obs, char *key, char *key_err) {
  sprintf(key,"%c%s", FILTERSTRING[ifilt_obs], HOSTLIB_SUFFIX_MAGOBS );
  sprintf(key_err,"%c%s", FILTERSTRING[ifilt_obs], HOSTLIB_SUFFIX_MAGOBS_ERR);

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
  // where j1 is the sparse SERSIC-array index (1-SERSIC.NPROF)
  // and j2 is the absolute 0-MXSERSIC integer key. Note that 
  // neither j1 or j2 is the physical SERSIC.n[j1] index.
  //
  // The index-matching is done here after read_head_HOSTLIB()
  // so that the VARNAMES and SERSIC keys can go in any order
  // in the HOSTLIB
  //

  int 
    NPROF, NFIX, j, ifix, NFIX_MATCH, NERR
    , IVAR_a, IVAR_b, IVAR_w, IVAR_n    
    ;

  double FIXn;
  char anam[12], bnam[12], wnam[12], nnam[12], tmpnam[12] ;
  char fnam[] = "init_Sersic_VARNAMES" ;

  // ----------- BEGIN ------------

  NPROF = 0;

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
      SERSIC_PROFILE.IVAR_a[NPROF] = IVAR_a ;
      SERSIC_PROFILE.IVAR_b[NPROF] = IVAR_b ;
      SERSIC_PROFILE.IVAR_w[NPROF] = IVAR_w ;
      SERSIC_PROFILE.IVAR_n[NPROF] = IVAR_n ;

      sprintf(SERSIC_PROFILE.VARNAME_a[NPROF], "%s", anam);
      sprintf(SERSIC_PROFILE.VARNAME_b[NPROF], "%s", bnam);
      sprintf(SERSIC_PROFILE.VARNAME_w[NPROF], "%s", wnam);
      sprintf(SERSIC_PROFILE.VARNAME_n[NPROF], "%s", nnam);   
      NPROF++ ;   
    }
  } // end of j-loop

  SERSIC_PROFILE.NPROF = NPROF ;
  if ( NPROF == 0 ) { return ; }


  // check list of FIXED Sersic indices and make sure
  // that they all match an a#,b# pair.
  NFIX_MATCH = 0;
  NFIX   = SERSIC_PROFILE.NFIX;

 
  for ( ifix=0; ifix < NFIX; ifix++ ) {
    sprintf(tmpnam, "%s", SERSIC_PROFILE.FIX_NAME[ifix]  ) ;

    for ( j=0; j < NPROF; j++ ) {
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
  for ( j=0; j < NPROF; j++ ) {
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
  double  Rmax, Rmin, logRbin, logRmin, logRmax, logR, dif, xmem, SCALE  ;
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
  SERSIC_TABLE.INVINDEX_BIN = dif / ((double)NBIN - 1.0 ) ; 


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

    SCALE = INPUTS.HOSTLIB_SCALE_SERSIC_SIZE ;
    if ( fabs(SCALE-1.0)>1.0E-5  ) 
      { printf("\t Sersic(a0,b0) scale:  %.3f \n", SCALE);    }
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
void get_Sersic_info(int IGAL, SERSIC_DEF *SERSIC) {

  // Fill SERSIC info a, b, n(and bn), w, wsum.
  // Also fill a_rot = rotation angle w.r.t. RA.
  // Require at least NPROF-1 weights to be defined;
  // the last weight is simply 1 = sum(other wgts).
  // If NPROF=1 then no weights are required.
  //
  // Feb 2020:
  //  + For FIXSERSIC option, set HOSTLIB.VALUE_ZSORTED
  //    so that fix Sersic params show up in SIMGEN_DUMP file.
  //
  // May 12 2020: 
  //   + allow wgt-truncation slop below E-3 by re-normalizing weights.
  //     (bug found by M.Vincenzi using MICECAT with double Sersic profile)
  //
  int IVAR_ANGLE   = HOSTLIB.IVAR_ANGLE ;
  double FIXa      = INPUTS.HOSTLIB_FIXSERSIC[0] ;
  double FIXb      = INPUTS.HOSTLIB_FIXSERSIC[1] ;
  double FIXn      = INPUTS.HOSTLIB_FIXSERSIC[2] ;
  double FIXANG    = INPUTS.HOSTLIB_FIXSERSIC[3];

  int j, NPROF, NWGT, j_nowgt, IVAR_a, IVAR_b, IVAR_w, IVAR_n ;
  int NPROF_ab0 = 0;
  double WGT, WGTSUM, WTOT, n, wsum_last, a, b ;
  char fnam[] = "get_Sersic_info" ;

  // -------------- BEGIN -------------

  NPROF = SERSIC_PROFILE.NPROF; 
  SERSIC->w[0]    = 0.0 ;
  SERSIC->wsum[0] = 0.0 ;
  SERSIC->NPROF   = NPROF; 

  for ( j=0; j < NPROF; j++ ) {
    IVAR_a = SERSIC_PROFILE.IVAR_a[j] ;
    IVAR_b = SERSIC_PROFILE.IVAR_b[j] ;
    IVAR_n = SERSIC_PROFILE.IVAR_n[j] ;

    // check for option(s) to fix Sersic params
    if ( FIXa > 0.0 ) 
      { HOSTLIB.VALUE_ZSORTED[IVAR_a][IGAL] = FIXa; }
    if ( FIXb > 0.0 ) 
      { HOSTLIB.VALUE_ZSORTED[IVAR_b][IGAL] = FIXb; }
    if ( FIXn > -998.0 && IVAR_n >= 0 ) 
      { HOSTLIB.VALUE_ZSORTED[IVAR_n][IGAL] = FIXn; }
    if ( FIXn > -998.0 && IVAR_n < 0 ) 
      { SERSIC_PROFILE.FIXn[j] = FIXn; }
    if ( FIXANG > -998.0 ) 
      { HOSTLIB.VALUE_ZSORTED[IVAR_ANGLE][IGAL] = FIXANG; }

    // - - - - - - 
    if ( IVAR_n >= 0 ) 
      { n = HOSTLIB.VALUE_ZSORTED[IVAR_n][IGAL] ; }
    else
      { n = SERSIC_PROFILE.FIXn[j] ; }

    SERSIC->a[j]  = HOSTLIB.VALUE_ZSORTED[IVAR_a][IGAL] ; 
    SERSIC->b[j]  = HOSTLIB.VALUE_ZSORTED[IVAR_b][IGAL] ; 
    SERSIC->n[j]  = n ;
    SERSIC->bn[j] = get_Sersic_bn(n);
    SERSIC->a_rot = HOSTLIB.VALUE_ZSORTED[IVAR_ANGLE][IGAL] ; 

    // apply user-scale on size (Mar 28 2018)
    SERSIC->a[j] *= INPUTS.HOSTLIB_SCALE_SERSIC_SIZE ;
    SERSIC->b[j] *= INPUTS.HOSTLIB_SCALE_SERSIC_SIZE ;
    
    if ( n < SERSIC_INDEX_MIN || n > SERSIC_INDEX_MAX ) {
      sprintf(c1err,"Sersic index=%f outside valid range (%5.2f-%5.2f)",
	      n, SERSIC_INDEX_MIN, SERSIC_INDEX_MAX ) ;
      sprintf(c2err,"Check GALID = %lld  ZTRUE=%f ", 
	      SNHOSTGAL.GALID, SNHOSTGAL.ZTRUE );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    if ( SERSIC->a[j] < 1.0E-9 || SERSIC->b[j] < 1.0E-9 ) 
      { NPROF_ab0++ ; }

  } // end NPROF


  // - - - - - - - - - - - - - - -
  // abort if all profile sizes are zero
  if ( NPROF_ab0 == NPROF ) {
    sprintf(c1err,"Profile size is zero for GALID=%lld", SNHOSTGAL.GALID);
    sprintf(c2err,"Check Sersic params.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // -----------------------------------
  // Now the weights:
  // check how many weights are defined, and track the sum.
  
  NWGT = 0;  IVAR_w=0;  j_nowgt = -9 ;     WGTSUM = WGT = 0.0 ;
  for ( j=0; j < NPROF; j++ ) {
    IVAR_w = SERSIC_PROFILE.IVAR_w[j] ;
    if ( IVAR_w > 0 ) { 
      NWGT++ ;
      WGT     = HOSTLIB.VALUE_ZSORTED[IVAR_w][IGAL] ;
      WGTSUM += WGT;
      SERSIC->w[j] = WGT ;
    }
    else { j_nowgt = j ; }
  }


  if ( NWGT == NPROF-1 ) {
    SERSIC->w[j_nowgt]  = 1.0 - WGTSUM ;
  }

  // check debug option for fixed weight with 2 profiles
  if ( NPROF == 2 && DEBUG_WGTFLUX2 > 0.00001 ) {
    SERSIC->w[0]  = 1.0 - DEBUG_WGTFLUX2 ;
    SERSIC->w[1]  = DEBUG_WGTFLUX2 ;
    goto WGTSUM ;
  }


  if ( NWGT < NPROF - 1 ) {
    sprintf(c1err,"Inadequate Sersic weights.");
    sprintf(c2err,"%d Sersic terms, but only %d weights defined.",
	    NPROF, NWGT );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

 WGTSUM:

  // fill cumulative WGTFLUX_SUM array
  wsum_last = 0.0 ;
  for ( j=0; j < NPROF; j++ ) {
    WGT = SERSIC->w[j] ;
    SERSIC->wsum[j] = WGT + wsum_last; 
    wsum_last = SERSIC->wsum[j] ; 
  }

  // finally check that sum of weights are one within E-3,
  // which allows a little  truncation slop in writing the HOSTLIB
  WTOT = SERSIC->wsum[NPROF-1] ;
  if ( fabs(WTOT-1.0) > 1.0E-3 ) {
    print_preAbort_banner(fnam);
    printf("\t NPROF = %d   (NWGT=%d, j_nowgt=%d) \n", 
	   NPROF, NWGT, j_nowgt);
    for ( j=0; j < NPROF; j++ ) { 
      printf("\t WGT[%d] = %12.5le  WGTSUM = %12.5le \n", 
	     j, SERSIC->w[j], SERSIC->wsum[j] ); 
    }
    sprintf(c1err,"Sum of Sersic weights = %f (ne 1) for GALID=%lld", 
	    WTOT, SNHOSTGAL.GALID );
    sprintf(c2err,"%s", "Check values of w1, w2 ...");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // May 12 2020
  // To allow truncation slop in writing HOSTLIB, 
  // if WTOT is slightly off of 1, re-normalize the weights
  if ( fabs(WTOT-1.0) > 1.0E-12 ) {
    double w_scale = 1.0/WTOT;
    for ( j=0; j < NPROF; j++ ) 
      { SERSIC->w[j] *= w_scale ;  SERSIC->wsum[j] *= w_scale ;  }
  }
  return ;

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

  return ;

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

  bool DO_GALMAG       = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG);
  bool DO_SN2GAL_Z     = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_Z);
  bool DO_SN2GAL_COORD = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC);
  bool DO_SWAPZPHOT    = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SWAPZPHOT ) ;
  bool DO_VPEC         = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEVPEC ) ;

  int    NTMP, ivar, NVAR, NROW, j ;
  double RAD, fixran, *fixab, SCALE ;
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
  sprintf(cptr, "HOSTLIB ZGAPMAX = %6.4f (%6.4f-%6.4f)  <ZGAP>=%10.3le",
          HOSTLIB.ZGAPMAX, HOSTLIB.Z_ATGAPMAX[0], HOSTLIB.Z_ATGAPMAX[1],
	  HOSTLIB.ZGAPAVG );
  
  // --------------------------------------
  // check HOSTLIB options and abort on inconsistencies.

  sprintf(copt,"HOSTLIB Opt: ");

  if ( DO_GALMAG ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "%s compute host-noise contribution to SN noise", copt );    
  }

  if ( DO_SN2GAL_Z ) {
    cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
    sprintf(cptr, "%s change SN redshift to host redshift", copt );    
  }

  if ( DO_SN2GAL_COORD ) {
    cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
    sprintf(cptr, "%s change SN position to host-SN pos.", copt );    

    if ( HOSTLIB.IVAR_RA < 0 || HOSTLIB.IVAR_DEC < 0 ) {
      sprintf(c1err ,"%s", "User requested SN pos -> host-SN pos");
      sprintf(c2err ,"%s", "but RA and/or DEC are missing from HOSTLIB.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  if ( DO_SWAPZPHOT ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "%s set ZTRUE = ZPHOT ", copt);
  }

  if ( DO_VPEC ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "%s Use VPEC (RMS=%.0f  min/max=%.0f/%.0f km/sec)", 
	    copt, HOSTLIB.VPEC_RMS, HOSTLIB.VPEC_MIN, HOSTLIB.VPEC_MAX);
  }

  // - - - - -
  SCALE = INPUTS.HOSTLIB_SCALE_SERSIC_SIZE;
  if ( fabs(SCALE-1.0) > 0.00001 ) {
    cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
    sprintf(cptr, "\t Scale Sersic size by %.3f ", SCALE);
  }
  SCALE = INPUTS.HOSTLIB_SCALE_LOGMASS_ERR;
  if ( fabs(SCALE-1.0) > 0.00001 ) {
    cptr = HOSTLIB.COMMENT[NTMP];  NTMP++ ; 
    sprintf(cptr, "\t Scale LOGMASS_ERR by %.3f ", SCALE);
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

  fixab = INPUTS.HOSTLIB_FIXSERSIC ;
  if ( fixab[0] > 0.000001 || fixab[1]>0.0001 || fixab[2]>-998.0 ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "\t DEBUG OPT -> "
	    "Fix Sersic a,b,n = %.2f,%.2f,%.2f  a_rot=%.1f deg \n", 
	    fixab[0], fixab[1], fixab[2], fixab[3] );
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
  }


  // summarize Sersic profiles
  for ( j=0; j < SERSIC_PROFILE.NPROF; j++ ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr,"GALSHAPE profile : %s %s"
	    ,SERSIC_PROFILE.VARNAME_a[j]
	    ,SERSIC_PROFILE.VARNAME_b[j] );
  }

  cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
  sprintf(cptr,"GALPOS generated with %5.4f -- %5.4f of total flux.",
	  INPUTS.HOSTLIB_MNINTFLUX_SNPOS, INPUTS.HOSTLIB_MXINTFLUX_SNPOS);


  // print out list of aperture radii
  if ( DO_GALMAG ) {
    cptr = HOSTLIB.COMMENT[NTMP]; NTMP++ ; 
    sprintf(cptr, "GALMAG %s interp-grid for PSFSIG(asec) = ",
	    HOSTLIB.filterList );
    for ( j=1; j <= NMAGPSF_HOSTLIB; j++ ) { 
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


// ================================
double get_VALUE_HOSTLIB(int ivar, int igal) {
  // Created Nov 2019
  // Returns VALUE for index ivar, from redshift-sorted galaxy list
  double VALUE;
  char fnam[] = "get_VALUE_HOSTLIB";
  // -------------- BEGIN ---------------
  VALUE = -9.0 ;
  if ( HOSTLIB.SORTFLAG == 0 ) {
    sprintf(c1err,"Cannot return sorted VALUE(ivar=%d,igal=%d)", ivar, igal);
    sprintf(c2err,"until HOSTLIB is redshift-sorted.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( ivar >= 0 )  { VALUE = HOSTLIB.VALUE_ZSORTED[ivar][igal] ; }

  return(VALUE) ;

} // end get_VALUE_HOSTLIB


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
int IVAR_HOSTLIB(char *varname, int FLAG) {

  // Mar 14, 2011
  // For input variable name return 'IVAR' index
  // to be used  the array HOSTLIB.VALUE_ZSORTED[ivar].
  // 
  // FLAG += 1  : abort if key not found; else return -9
  // FLAG += 2  : case sensitive compar (default is case insensitive)
  //
  // Feb 2014: use case-insenstive string check.
  // Apr 2021: FLAG +=2 for case-sensitive compare; e.g., b_obs != B_obs
  //
  int ivar, NVAR, ICMP ;
  bool ABORTFLAG = ( FLAG & 1 ) ;
  bool CASE      = ( FLAG & 2 ) ; 
  char fnam[] = "IVAR_HOSTLIB" ;

  // ---------- BEGIN ----------
  
  
  NVAR = HOSTLIB.NVAR_STORE ; 
  for ( ivar = 0; ivar < NVAR; ivar++ ) {

    if ( CASE ) 
      {  ICMP = strcmp(varname,HOSTLIB.VARNAME_STORE[ivar]) ; }
    else
      { ICMP = strcmp_ignoreCase(varname,HOSTLIB.VARNAME_STORE[ivar]) ; }
    
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
int IVAR_HOSTLIB_PREFIX(char *prefix, int FLAG) {

  // Feb 10, 2022
  // For input prefix return 'IVAR' index
  // of first instance containing that prefix
  // to be used  the array HOSTLIB.VALUE_ZSORTED[ivar].
  // Comparisons are case-sensitive. 
  //
  // FLAG += 1  : abort if key not found; else return -9
  //
  int ivar, NVAR, ICMP ;
  bool ABORTFLAG = ( FLAG & 1 ) ;
  bool CASE      = ( FLAG & 2 ) ;
  char fnam[] = "IVAR_HOSTLIB_PREFIX" ;

  // ---------- BEGIN ----------

  NVAR = HOSTLIB.NVAR_STORE ;
  for ( ivar = 0; ivar < NVAR; ivar++ ) {
        if (strstr(HOSTLIB.VARNAME_STORE[ivar],prefix)  != NULL){
	      return(ivar);
      }
  }

  // if we get here, abort
  if ( ABORTFLAG ) {
    sprintf(c1err,"Could not find IVAR index for prefix='%s'", prefix);
    sprintf(c2err,"Check VARNAMES keys in HOSTLIB file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    return(-9) ;
  }
  else  {
    return(-9) ;
  }
} // end of IVAR_HOSTLIB_PREFIX

bool ISCHAR_HOSTLIB(int IVAR) {
  // Feb 25 2020
  // return true if input IVAR corresponds to a string variable.
  // IVAR is stored index, not ALL index.
  if ( IVAR == HOSTLIB.IVAR_FIELD )     { return(true); }
  if ( IVAR == HOSTLIB.IVAR_NBR_LIST )  { return(true); }
  return(false);

} // end ISCHAR_HOSTLIB

// ========================================
int ICOL_SPECTABLE(char *varname, int ABORTFLAG) {


  // June 2019
  // For input variable name return 'ICOL' column in specTemplate file.
  // 
  // ABORTFLAG = 1  : abort if key not found
  // ABORTFLAG = 0  : return -9 if key nor found
  //

  int icol, NCOL, ICMP ;
  char VARNAME_TMP[2][60];
  char fnam[] = "ICOL_SPECTABLE";

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

} // end of ICOL_SPECTABLE

// =========================================
void GEN_SNHOST_DRIVER(double ZGEN_HELIO, double PEAKMJD) {

  // Mar 2011
  // Driver function to select host-galaxy from library,
  // and to determine other properties: 
  // photoz, SN pos, flux-noise ...
  // Note that ZGEN_HELIO is the desired host HELIOCENTRIC redshift,
  // which could be different than SN redshift if GENLC.CORRECT_HOSTMATCH=0.
  //
  //
  // Jan 31 2020: 
  //  little more refactor to use DDLR sorting. Note that 
  //  GEN_SNHOST_ZPHOT(IGAL) is moved after DDLR sorting,
  //  which changes random sync.
  //

  int    USE, IGAL, ilist ;
  double fixran ;
  char   fnam[] = "GEN_SNHOST_DRIVER" ;

  // ------------ BEGIN -----------

  NCALL_GEN_SNHOST_DRIVER++ ;

  // always burn random numbers to stay synced.
  ilist = 1 ; 
  SNHOSTGAL.FlatRan1_GALID     = getRan_Flat1(ilist) ; // random GAL in z-bin
  SNHOSTGAL.FlatRan1_radius[0] = getRan_Flat1(ilist) ; // ran Sersic profile
  SNHOSTGAL.FlatRan1_radius[1] = getRan_Flat1(ilist) ; // random integral
  SNHOSTGAL.FlatRan1_phi       = getRan_Flat1(ilist) ;  // relative to major axis

  // ------------------------------------------------
  // check option to fix randoms (Sep 14, 2012)
  fixran = INPUTS.HOSTLIB_FIXRAN_RADIUS ;
  if ( fixran > -1.0E-9 ) { SNHOSTGAL.FlatRan1_radius[1] = fixran ; }

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
    double zCMB;
    zCMB = zhelio_zcmb_translator(ZGEN_HELIO,GENLC.RA,GENLC.DEC,"eq",+1);
    sprintf(c1err,"Invalid ZGEN(Helio,CMB)=%f,%f ", ZGEN_HELIO, zCMB );
    sprintf(c2err,"HOSTLIB z-range is %6.4f to %6.4f",
	    HOSTLIB.ZMIN, HOSTLIB.ZMAX );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // select host-galaxy from library
  GEN_SNHOST_GALID(ZGEN_HELIO);
  IGAL = SNHOSTGAL.IGAL ;
  if ( IGAL < 0 ) { return ; } // Aug 2015

  // generate SN position at galaxy
  GEN_SNHOST_POS(IGAL);

  // - - - - - - - - - - - - - -
  // check for neighbors
  GEN_SNHOST_NBR(IGAL);
    
  // determine DLR and ordered list
  for(ilist=0; ilist < SNHOSTGAL.NNBR; ilist++ ) 
    { GEN_SNHOST_DDLR(ilist); }

  // sort by DDLR
  SORT_SNHOST_byDDLR();

  // - - - - - - 
  // check on host photoz
  GEN_SNHOST_ZPHOT(IGAL);

  // check for vpec
  GEN_SNHOST_VPEC(IGAL);

  // check if redshift needs to be updated (Apr 8 2019)
  TRANSFER_SNHOST_REDSHIFT(IGAL);

  // check on host properties for each possible host match (Feb 2020)
  int ivar_property;
  for (ivar_property=0; ivar_property<N_HOSTGAL_PROPERTY; ivar_property++){
    GEN_SNHOST_PROPERTY(ivar_property);
  }

  // host-mag within SN aperture
  GEN_SNHOST_GALMAG(IGAL);

  // load user-specified variables for output file
  LOAD_OUTVAR_HOSTLIB(IGAL); 

  // ---------------------------
  STORE_SNHOST_MISC(IGAL,HOSTLIB_WGTMAP.ibin_SNVAR);


  // check debug-dump options
  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_DUMPROW ) 
    {  DUMPROW_SNHOST(); }

  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_DUMP ) 
    { DUMP_SNHOST(); }


} // end of GEN_SNHOST_DRIVER



// =========================================
void GEN_SNHOST_GALID(double ZGEN) {
  
  // Mar 2011
  // Select weighted hostlib entry for heliocentric redshift ZGEN(SN)
  // Fills
  // * SNHOSTGAL.IGAL
  // * SNHOSTGAL.GALID
  // * SNHOSTGAL.ZTRUE (GAL)
  //
  // Nov 23 2019: for MODEL_SIMLIB, force GALID to value in SIMLIB header.
  // Dec 30 2021: minor refactor to make igal loops faster with binary search.

  int  IGAL_JUMP           = 20 ; // spped search for approx start range
  int  IGAL_RANGE_CONVERGE = 5;   // convergence for binary search
  int  NGAL_CHECK_ABORT    = 50;  // avoid infinite loop in binary search
  bool ISMODEL_SIMLIB   = ( INDEX_GENMODEL == MODEL_SIMLIB ) ;
  int    USEONCE        = INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEONCE;
  double DZPTR          = DZPTR_HOSTLIB ;
  double FlatRan1_GALID = SNHOSTGAL.FlatRan1_GALID ;
  int N_SNVAR           = HOSTLIB_WGTMAP.N_SNVAR ;

  int  IZ_CEN, IZ_TOLMIN, IZ_TOLMAX ;
  int  IGAL_SELECT, igal_start, igal_end, igal;
  int  igal0, igal1, igal_middle;
  int  igal_start_init, igal_end_init ;
  int  NSKIP_WGT, NSKIP_USED, NGAL_CHECK, MATCH, ibin_SNVAR=-9; 
  long long GALID ;
  double ZTRUE, LOGZGEN, LOGZTOLMIN, LOGZTOLMAX, LOGZDIF ;
  double WGT_start, WGT_end, WGT_dif, WGT_select, WGT, *ptrWGT ;
  double ztol, dztol, z, z_start, z_end ;

  int  LDMP   = 0; // (GENLC.CID>50000 && GENLC.CID < 50005) ;
  char fnam[] = "GEN_SNHOST_GALID" ;
  bool REFAC;
  // ---------- BEGIN ------------

  IGAL_SELECT = -9 ; 
  
  // compute zSN-zGAL tolerance for this ZGEN = zSN
  dztol = eval_GENPOLY(ZGEN, &INPUTS.HOSTLIB_GENPOLY_DZTOL, fnam) ;
    
  // find start zbin 
  LOGZGEN    = log10(ZGEN);
  LOGZTOLMAX = log10(ZGEN+dztol);
  if ( (ZGEN-dztol) > ZMIN_HOSTLIB )
    { LOGZTOLMIN = log10(ZGEN-dztol); }
  else
    { LOGZTOLMIN = MINLOGZ_HOSTLIB ; }

  if ( LOGZGEN < MINLOGZ_HOSTLIB ) {
    sprintf(c1err,"LOGZGEN=%.3f < MINLOGZ_HOSTLIB=%.3f",
	    LOGZGEN, MINLOGZ_HOSTLIB);
    sprintf(c2err,"Increase GENRANGE_REDSHIFT[0] or "
	    "descrease MINLOGZ_HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // compute approx IZ index at z, z-dztol, z+dztol
  LOGZDIF    = LOGZGEN - MINLOGZ_HOSTLIB;
  IZ_CEN     = (int)( LOGZDIF/DZPTR ) ; 

  LOGZDIF    = LOGZTOLMIN - MINLOGZ_HOSTLIB;
  IZ_TOLMIN  = (int)( LOGZDIF/DZPTR + 0.5 ) ; 

  LOGZDIF    = LOGZTOLMAX - MINLOGZ_HOSTLIB;
  IZ_TOLMAX  = (int)( LOGZDIF/DZPTR + 0.5 ) ; 

  // - - - - - - -
  // select min GALID using dztol from user

  // begin with approx calculation using IZPTR grid
  if ( IZ_CEN >= HOSTLIB.MAXiz ) 
    { igal_start = HOSTLIB.NGAL_STORE-1; }
  else {
    // xxx mark delete igal_start = HOSTLIB.IZPTR[IZ_CEN+1]; 
    igal_start = HOSTLIB.IZPTR[IZ_TOLMIN+1];
  }

  // loop every IGAL_JUMP galaxy to get refined range
  igal_start_init = igal_start;
  z = get_ZTRUE_HOSTLIB(igal_start) ;  ztol=ZGEN-dztol ; 
  while ( z > ztol && igal_start > IGAL_JUMP ) {
    igal_start -= IGAL_JUMP ; 
    z  = get_ZTRUE_HOSTLIB(igal_start);  
  }

  if ( igal_start != igal_start_init ) 
    { igal_start += IGAL_JUMP ;  z = get_ZTRUE_HOSTLIB(igal_start); }


  // final loop one galaxy at a time to find exact igal at z < ztol
  while ( z > ztol && igal_start > 1 )  
    {  igal_start-- ; z = get_ZTRUE_HOSTLIB(igal_start);  }

  // - - - - - - - - - 
  // select max GALID using dztol from user
  if ( IZ_CEN < HOSTLIB.MINiz ) 
    { igal_end = 0; }
  else { 
    // xxx mark    igal_end = HOSTLIB.IZPTR[IZ_CEN-1]; 
    igal_end = HOSTLIB.IZPTR[IZ_TOLMAX-1]; 
  }

  igal_end_init = igal_end;
  z = get_ZTRUE_HOSTLIB(igal_end) ; ztol = ZGEN+dztol ;

  while ( z < ztol && igal_end < HOSTLIB.NGAL_STORE-IGAL_JUMP ) { 
    igal_end += IGAL_JUMP ; 
    z = get_ZTRUE_HOSTLIB(igal_end);  
  }
  
  if ( igal_end != igal_end_init )  { 
    igal_end -= IGAL_JUMP; 
    z = get_ZTRUE_HOSTLIB(igal_end);
  }

  while ( z < ztol && igal_end < HOSTLIB.NGAL_STORE-1 ) 
    { igal_end++ ; z = get_ZTRUE_HOSTLIB(igal_end);  }

  // - - - - - - - - - - - - - - - - - - - - - - - -  -
  // back up one step to stay inside dztol
  igal_start++ ;    igal_end--   ; 
  z_start = get_ZTRUE_HOSTLIB(igal_start); 
  z_end   = get_ZTRUE_HOSTLIB(igal_end); 

  if ( LDMP ) {
    printf(" xxx ---------- %s DUMP ---------- \n", fnam);
    printf(" xxx CID=%d  ZGEN = %f  dztol=%f \n", 
	   GENLC.CID, ZGEN, dztol );
    printf(" xxx igal_start = %d -> %d (loop %d)\n", 
	   igal_start_init, igal_start, igal_start_init-igal_start);
    printf(" xxx igal_end   = %d -> %d (loop %d)\n", 
	   igal_end_init, igal_end, igal_end-igal_end_init );
    fflush(stdout);
  }

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
  if ( N_SNVAR > 0 ) {
    ibin_SNVAR = getBin_SNVAR_HOSTLIB_WGTMAP();
    HOSTLIB_WGTMAP.ibin_SNVAR = ibin_SNVAR ;
    ptrWGT     = HOSTLIB_WGTMAP.WGTSUM_SNVAR[ibin_SNVAR];
  }
  else {
    ptrWGT     = HOSTLIB_WGTMAP.WGTSUM;
  }

  WGT_start  = ptrWGT[igal_start];
  WGT_end    = ptrWGT[igal_end];    
  WGT_dif    = WGT_end - WGT_start ;
  WGT_select = WGT_start + ( WGT_dif * FlatRan1_GALID * 0.95 ) ;

  NSKIP_WGT   = NSKIP_USED = NGAL_CHECK = 0 ;
  // xxx mark  GALID_FORCE = INPUTS.HOSTLIB_GALID_FORCE ;

  if ( ISMODEL_SIMLIB  &&  SIMLIB_HEADER.GALID > 0 ) { 
    INPUTS.HOSTLIB_GALID_PRIORITY[0] = SIMLIB_HEADER.GALID;
    INPUTS.HOSTLIB_GALID_PRIORITY[1] = SIMLIB_HEADER.GALID;
  }

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

  
  // check GALID_FORCE here (Dec 29 2021)
  if ( HOSTLIB.IGAL_FORCE > 0 ) {
    IGAL_SELECT = HOSTLIB.IGAL_FORCE;
    goto DONE_SELECT_GALID ;
  }

  // ---------------------------------------------------


  // perform binary search to restrict igal range to within a few galaxies.
  bool CONVERGE = false;
  igal0 = igal_start;
  igal1 = igal_end;    
  igal_middle = (int)( (igal0 + igal1)/2 );
  
  while ( !CONVERGE ) {
    WGT = ptrWGT[igal_middle] ;
    NGAL_CHECK++ ;
    if ( NGAL_CHECK > NGAL_CHECK_ABORT ) {
      print_preAbort_banner(fnam);
      printf("\t CID=%d  ZGEN=%f\n", GENLC.CID, ZGEN);
      printf("\t WGT_select=%le   current WGT=%le\n", WGT_select, WGT);
      printf("\t igal[start,end]=%d,%d  igal[0,middle,1]=%d,%d,%d\n",
	     igal_start, igal_end,   igal0, igal_middle, igal1);
      sprintf(c1err,"Cannot converge finding igal range.");
      sprintf(c2err,"NGAL_CHECK=%d", NGAL_CHECK );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
    }
    if ( WGT < WGT_select )
      {  igal0 = igal_middle ;  } // select upper half for next iter
    else 
      {  igal1 = igal_middle ;   } // lower half for next iter

    igal_middle = (int)( (igal0 + igal1)/2 );
    CONVERGE = (igal1 - igal0) < IGAL_RANGE_CONVERGE ;
    
    if ( LDMP ) {
      double WGT_ratio = (WGT-WGT_start) / ( WGT_select-WGT_start);
      printf(" xxx igal[0,m,1] = %d, %d, %d   CONVERGE=%d "
	     "WGT_ratio=%.4f \n",
	     igal0, igal_middle, igal1, CONVERGE, WGT_ratio);
    }
    
  } // end while not CONVERGE

    // reset igal_start[end] for brute-force search below.
  igal_start = igal0;
  if ( !USEONCE ) { igal_end = igal1; }

  // - - - - - - - - - - - 
  // Brute force search, one igal at a time.

  for ( igal = igal_start; igal <= igal_end; igal++ ) {
    NGAL_CHECK++ ;
    WGT = ptrWGT[igal];

    if ( WGT <  WGT_select  )  
      { NSKIP_WGT++; continue ; }
    
    if ( USEHOST_GALID(igal) == 0 ) 
      { NSKIP_USED++ ; continue ; }
    
    // select first igal with WGT > WGT_select
    IGAL_SELECT = igal ;
    goto DONE_SELECT_GALID ;
    
  } // end igal loop

  // - - - - - - - - - - - - - - - - - - - - - - - - 

 DONE_SELECT_GALID:

  if ( LDMP ) {
    bool GOOD = ( IGAL_SELECT >= igal_start && IGAL_SELECT <= igal_end);
    printf(" xxx IGAL_SEL=%d  igal_range=%d,%d [GOOD=%d]  NGAL_CHECK=%d\n",
	   IGAL_SELECT, igal_start, igal_end, GOOD, NGAL_CHECK ); 
    fflush(stdout);
  }

  if ( IGAL_SELECT < 0 ) {
    
    // if using same Galaxy with MJD-sep, just return so that
    // this event is rejected.
    if ( INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL < 9999 ) 
      { SNHOSTGAL.IGAL = IGAL_SELECT ; return ; }

    print_preAbort_banner(fnam);
    printf("\t NSKIP_WGT  = %d/%d \n", NSKIP_WGT,  NGAL_CHECK );
    printf("\t NSKIP_USED = %d/%d \n", NSKIP_USED, NGAL_CHECK );
    printf("\t igal(start-end) = %d - %d\n", igal_start, igal_end);
    printf("\t WGT(start-end)  = %f - %f\n",  WGT_start, WGT_end );
    printf("\t WGT_select      = %f \n", WGT_select);

    DUMP_SNHOST();
    fflush(stdout);

    sprintf(c1err,"Could not find HOSTLIB entry for ZGEN=%5.4f .", ZGEN);
    if ( USEONCE ) 
      { sprintf(c2err,"Each HOST used only once -> need larger library."); }
    else
      { sprintf(c2err,"This error is baffling ?!?!?!"); }
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
    print_preAbort_banner(fnam);
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

  SNHOSTGAL.a_SNGALSEP_ASEC   = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.b_SNGALSEP_ASEC   = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.RA_SNGALSEP_ASEC  = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.DEC_SNGALSEP_ASEC = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.RA_GAL_DEG        = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.DEC_GAL_DEG       = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.RA_SN_DEG         = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.DEC_SN_DEG        = HOSTLIB_SNPAR_UNDEFINED ;

  // always init GALMAG quantities to garbage
  int i, ifilt ;
  for ( i=0; i <= NMAGPSF_HOSTLIB ; i++ ) {
    SNHOSTGAL.GALFRAC[i]          = -9.0 ;     // global
    for ( ifilt = 0; ifilt < MXFILTINDX; ifilt++ ) { 
      SNHOSTGAL.GALMAG[ifilt][i]  = MAG_UNDEFINED ;
      SNHOSTGAL.SB_FLUXCAL[ifilt] = 0.0 ; 
      SNHOSTGAL.SB_MAG[ifilt]     = MAG_UNDEFINED ;
    }
  }

  for ( i=0; i < MXNBR_LIST ; i++ ) {
    SNHOSTGAL.IGAL_NBR_LIST[i]  = HOSTLIB_IGAL_UNDEFINED ; //
    SNHOSTGAL.DDLR_NBR_LIST[i]  = HOSTLIB_SNPAR_UNDEFINED ;
    SNHOSTGAL.SNSEP_NBR_LIST[i] = HOSTLIB_SNPAR_UNDEFINED ;
  }

  return ;

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
  int DAY, DAYDIF, ABSDIF, PEAKDAY ;
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
    if ( SAMEHOST.NUSE[IGAL] == 0 )  { SAMEHOST.NUSE[IGAL]++ ; } 
    retCode = 1 ;
  }
  else if ( SAMEHOST.REUSE_FLAG == 2 ) {
    // Jan 2022 RK - there is a rare bug causing NUSE=0 ??
    // option for re-use of host with PEAKMD-sep constraint
    BLOCKSUM = 0 ;
    MINDAYSEP = INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL ;
    DAY = (int)SNHOSTGAL.PEAKMJD - SAMEHOST.PEAKMJD_STORE_OFFSET ;

    // allocate memory for this IGAL on first usage
    // Problem: if this IGAL is rejected, then malloc might happen again later??
    if ( NUSE_PRIOR == 0 ) {
      SAMEHOST.PEAKDAY_STORE[IGAL] = 
	(unsigned short*)malloc( MXUSE_SAMEGAL*sizeof(unsigned short) ) ;
 
      for(use=0; use < MXUSE_SAMEGAL; use++ )  // Jan 12 2022
        { SAMEHOST.PEAKDAY_STORE[IGAL][use] = 0 ; }
    }

    // check if a previous PEAKMJD on this host is too close in time.
    for ( use=0; use < NUSE_PRIOR; use++ )  {
      PEAKDAY =  SAMEHOST.PEAKDAY_STORE[IGAL][use] ;
      MJD     =  PEAKDAY + SAMEHOST.PEAKMJD_STORE_OFFSET ;
      DAYDIF  =  DAY - PEAKDAY ;
      if ( DAYDIF >= 0 ) { ABSDIF=DAYDIF; }    else { ABSDIF = -DAYDIF; }
      BLOCK_STORE[use] = 0 ;
      if ( ABSDIF < MINDAYSEP ) { BLOCKSUM = 1; BLOCK_STORE[use]=1; }

      MJD_STORE[use]    = MJD ;
      DAYDIF_STORE[use] = DAYDIF ;
    }
    
    // if current PEAKMJD was not blocked, then use this host [again]
    if ( BLOCKSUM == 0 ) {
      NUSE  = SAMEHOST.NUSE[IGAL] ;
      SAMEHOST.PEAKDAY_STORE[IGAL][NUSE] = DAY ;
      SAMEHOST.NUSE[IGAL]++ ;
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
    print_preAbort_banner(fnam);
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
  //
  // Jan 12 2022: check that NUSE > 0 before NUSE--.
  //     Beware that there may still be an underlying bug
  //     that rarely caused NUSE=0 -> NUSE = -1
  //
  if ( IGAL < 0 ) { return ; }

  if ( SAMEHOST.REUSE_FLAG == 2 && SAMEHOST.NUSE[IGAL] > 0 ) {
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
  // Jan 31 2020: compute for all neighbors

  int    NNBR  = SNHOSTGAL.NNBR;
  double ZTRUE = SNHOSTGAL.ZTRUE ;
  double ZGEN  = SNHOSTGAL.ZGEN ;
  int j, inbr;
  double ZPHOT, ZPHOT_ERR, ZBIAS, tmp ;
  char fnam[] = "GEN_SNHOST_ZPHOT" ;

  // ----------- BEGIN ----------

  // Compte user-specified bias
  ZBIAS = 0 ;
  for(j=0; j < 4; j++ ) {
    tmp = (double)INPUTS.HOSTLIB_GENZPHOT_BIAS[j];
    if ( tmp != 0.0 )  { ZBIAS += tmp * pow(ZTRUE, (double)j );  }
  }

  for(inbr=0; inbr < NNBR; inbr++ ) {
    if ( INPUTS.USE_HOSTLIB_GENZPHOT  ) 
      { GEN_SNHOST_ZPHOT_from_CALC(ZGEN, &ZPHOT,&ZPHOT_ERR); }
    else
      { GEN_SNHOST_ZPHOT_from_HOSTLIB(inbr, ZGEN, &ZPHOT,&ZPHOT_ERR); }

    SNHOSTGAL_DDLR_SORT[inbr].ZPHOT     = ZPHOT + ZBIAS;
    SNHOSTGAL_DDLR_SORT[inbr].ZPHOT_ERR = ZPHOT_ERR ;
  }


  // final ZPHOT is from host with closest match (Jan 31 2020)
  SNHOSTGAL.ZPHOT     = SNHOSTGAL_DDLR_SORT[0].ZPHOT;
  SNHOSTGAL.ZPHOT_ERR = SNHOSTGAL_DDLR_SORT[0].ZPHOT_ERR ;

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
void GEN_SNHOST_ZPHOT_from_CALC(double ZGEN, double *ZPHOT, double *ZPHOT_ERR) {

  // Created Feb 23 2017
  // Compute ZPHOT and ZPHOT_ERR from Gaussian profiles specified
  // by sim-input  HOSTLIB_GENZPHOT_FUDGEPAR or HOSTLIB_GENZPHOT_FUDGEMAP
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
  // Jan 31 2020: pass ZGEN instead of unused IGAL
  // Nov 15 2021: check for HOSTLIB_GENZPHOT_FUDGEMAP

  int      NzBIN             = INPUTS.HOSTLIB_GENZPHOT_FUDGEMAP.NzBIN ;
  float   *GENZPHOT_FUDGEPAR = INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR ;

  double  sigz1_core[3] ;  // sigma/(1+z): a0 + a1*(1+z) + a2*(1+z)^2
  double  sigz1_outlier, prob_outlier, zpeak, sigz_lo, sigz_hi ;
  double  sigma_core, sigma_outlier, HOSTLIB_ZRANGE[2] ;
  double  z1, ranProb, ranzFlat, zphotErr  ;

  GENGAUSS_ASYM_DEF ZPHOTERR_ASYMGAUSS ;

  int    OUTLIER_FLAT=0;
  double SIGMA_OUTLIER_FLAT = 9.999 ; // pick random z if sig_outlier> this
  char   fnam[] = "GEN_SNHOST_ZPHOT_from_CALC" ;

  // ----------- BEGIN -------------

  *ZPHOT = *ZPHOT_ERR = -9.0 ;
  
  if ( NzBIN > 0 ) {
    // use map of RMS vs. z
    double *z_LIST   = INPUTS.HOSTLIB_GENZPHOT_FUDGEMAP.z_LIST   ;
    double *RMS_LIST = INPUTS.HOSTLIB_GENZPHOT_FUDGEMAP.RMS_LIST ;
    sigz1_core[0] = interp_1DFUN(1, ZGEN, NzBIN, z_LIST, RMS_LIST, fnam);
    sigz1_core[1] = sigz1_core[2] = prob_outlier = sigz1_outlier = 0.0;

    //    printf(" xxx %s: zgen=%.3f -> sigz1 = %.3f \n",
    //	   fnam, ZGEN, sigz1_core[0] ); fflush(stdout);
  }
  else {
    // use 3rd-order polynomial and outlier Gaussian
    sigz1_core[0] = (double)GENZPHOT_FUDGEPAR[0] ;
    sigz1_core[1] = (double)GENZPHOT_FUDGEPAR[1] ;
    sigz1_core[2] = (double)GENZPHOT_FUDGEPAR[2] ;
    prob_outlier  = (double)GENZPHOT_FUDGEPAR[3] ;
    sigz1_outlier = (double)GENZPHOT_FUDGEPAR[4] ;
  }

  // - - - - - - - - - - - - 

  // xxx delete Jan 31 2020  z  = GENLC.REDSHIFT_HELIO ;

  z1 = 1.0 + ZGEN ;

  sigma_core = (sigz1_core[0] + 
		sigz1_core[1]*z1 + 
		sigz1_core[2]*(z1*z1) ) ;

  sigma_outlier = sigz1_outlier*z1 ;
  
 PICKRAN:

  ranProb  = getRan_Flat1(1) ;

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
    ranzFlat = getRan_Flat(2, HOSTLIB_ZRANGE); 
    *ZPHOT = ranzFlat;    
  }
  else {
    // to avoid negative ZPHOT (or truncated ZPHOT),
    // check to modify symmetric error into asymmetric error

    zphoterr_asym(ZGEN, zphotErr, &ZPHOTERR_ASYMGAUSS);
    zpeak   = ZPHOTERR_ASYMGAUSS.PEAK;
    sigz_lo = ZPHOTERR_ASYMGAUSS.SIGMA[0] ;
    sigz_hi = ZPHOTERR_ASYMGAUSS.SIGMA[1] ;
    *ZPHOT  = zpeak + getRan_GaussAsym(sigz_lo, sigz_hi, 0. );

    /*
    printf(" xxx zpeak=%.4f sig(-/+)=%.4f/%.4f  ZPHOT=%.3f \n",
	   zpeak, sigz_lo, sigz_hi, *ZPHOT); fflush(stdout);
    */

  }

  if ( *ZPHOT < 0.001 ) { goto PICKRAN; }

  // reported ZPHOT_ERR is sigma_core, or the RMS of asym Gaussian.
  zphoterr_asym(ZGEN, sigma_core, &ZPHOTERR_ASYMGAUSS) ;
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
void GEN_SNHOST_ZPHOT_from_HOSTLIB(int INBR, double ZGEN, 
				   double *ZPHOT, double *ZPHOT_ERR) {

  // Created Feb 23 2017:
  // Generate host-ZPHOT from host library.
  //
  // Note that reported ZPHOT += [ zGEN(HELIO) - ZTRUE(HOSTLIB) ]
  // However, ZPHOTERR is not adjusted, which assumes that the
  // library ZTRUE is always close to generated redshift.
  //
  // Inputs:
  //   INBR = neighbor index (0 -> true host)
  //   ZGEN = true helio redshift of SN
  //   
  // Outputs;
  //   *ZPHOT = photo z of host
  //   *ZPHOT_ERR = error on above
  //
  //
  // Mar 28 2018: float -> double
  // May 17 2018: if either zphot<0 or zphoterr<0, set both
  //              to -9. These indicate failed host photo-z fit.
  //
  // Apr 06 2019: initialize outputs (*ZPHOT,*ZPHOT_ERR) to -9.
  // Jan 31 2020: 
  //   + pass INBR (neighbor index) instead of IGAL
  //   + pass ZGEN instead of using GENLC.REDSHIFT_HELIO


  int IVAR_ZPHOT, IVAR_ZPHOT_ERR ;
  double ZDIF, ZTRUE, zphot_local, zerr_local ;
  char fnam[] = "GEN_SNHOST_ZPHOT_from_HOSTLIB" ;

  // ----------- BEGIN -----------

  *ZPHOT  = *ZPHOT_ERR = -9.0; 
  IVAR_ZPHOT     = HOSTLIB.IVAR_ZPHOT ;
  IVAR_ZPHOT_ERR = HOSTLIB.IVAR_ZPHOT_ERR ;
  if ( IVAR_ZPHOT < 0 ) { return ; }

  ZDIF = 0.0 ;
  if ( GENLC.CORRECT_HOSTMATCH ) {
    // from legacy wrong-host map
    ZTRUE = SNHOSTGAL_DDLR_SORT[INBR].ZSPEC ; // ZTRUE from hostlib
    ZDIF = ZGEN - ZTRUE ; 
  }

  zphot_local = SNHOSTGAL_DDLR_SORT[INBR].ZPHOT + ZDIF ;
  zerr_local  = SNHOSTGAL_DDLR_SORT[INBR].ZPHOT_ERR ;

  if ( zphot_local < 0.0 || zerr_local < 0.0 ) 
    { zphot_local = zerr_local = -9.0; }
  

  // load output args
  *ZPHOT     = zphot_local;
  *ZPHOT_ERR = zerr_local;

} // end GEN_SNHOST_ZPHOT_from_HOSTLIB

// =========================================
void  GEN_SNHOST_VPEC(int IGAL) {

  // Created May 24 2020
  // Check for VPEC column in HOSTLIB table.

  bool DO_VPEC       = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEVPEC) ;
  int  IVAR_VPEC     = HOSTLIB.IVAR_VPEC ;
  int  IVAR_VPEC_ERR = HOSTLIB.IVAR_VPEC_ERR ;
  double VPEC        = 0.0, ERR = -999.9 ;
  //  char fnam[]        = "GEN_SNHOST_VPEC" ;

  // ------------ BEGIN ----------

  if ( DO_VPEC ) {
    VPEC  = get_VALUE_HOSTLIB(IVAR_VPEC,IGAL); 
    if ( IVAR_VPEC_ERR > 0 ) 
      { ERR = get_VALUE_HOSTLIB(IVAR_VPEC_ERR,IGAL);  }
    else
      { ERR = HOSTLIB.FIX_VPEC_ERR ; }	
  }
  
  SNHOSTGAL.VPEC     = VPEC ;
  SNHOSTGAL.VPEC_ERR = ERR ;

  return ;

} // end GEN_SNHOST_VPEC

// =========================================
void GEN_SNHOST_LOGMASS(void) {

  // *** legacy function as of Feb 2022 ***
  // Created Feb 2020
  // If LOGMASS_OBS is defined in HOSTLIB, do nothing.
  // Otherwise, use LOGMASS_TRUE and LOGMASS_ERR to determine 
  // LOGMASS_OBS.
  // Alex Gagliano 10/12/21 Added if block to set LOGMASS_OBS = LOGMASS_TRUE
  // if LOGMASS_OBS and LOGMASS_ERR not in HOSTLIB

  int  NNBR       = SNHOSTGAL.NNBR;
  int  IVAR_TRUE  = HOSTLIB.IVAR_LOGMASS_TRUE ;
  int  IVAR_OBS   = HOSTLIB.IVAR_LOGMASS_OBS ;
  int  IVAR_ERR   = HOSTLIB.IVAR_LOGMASS_ERR ;
  double SCALE    = INPUTS.HOSTLIB_SCALE_LOGMASS_ERR;
  int i;

  double LOGMASS_TRUE, LOGMASS_OBS, LOGMASS_ERR, GauRan ;
  double rmin=-3.0, rmax=3.0 ;
  char fnam[] = "GEN_SNHOST_LOGMASS" ;

  // ---------- BEGIN -----------
  
  if ( IVAR_TRUE < 0 ) { return; }

  for(i=0; i < NNBR; i++ ) {

    LOGMASS_OBS = -9.0 ;

    if ( IVAR_OBS > 0 ) { 
      LOGMASS_OBS = SNHOSTGAL_DDLR_SORT[i].LOGMASS_OBS ;
    }
    else if ( IVAR_TRUE > 0 && IVAR_ERR > 0 ) {
      LOGMASS_TRUE = SNHOSTGAL_DDLR_SORT[i].LOGMASS_TRUE ;
      LOGMASS_ERR  = SNHOSTGAL_DDLR_SORT[i].LOGMASS_ERR ;
      LOGMASS_ERR *= SCALE ;
      GauRan = getRan_GaussClip(1,rmin,rmax);
      LOGMASS_OBS = LOGMASS_TRUE + GauRan*LOGMASS_ERR ;
    } 
    else {
      LOGMASS_TRUE = SNHOSTGAL_DDLR_SORT[i].LOGMASS_TRUE ;
      LOGMASS_OBS = LOGMASS_TRUE;
    }

    SNHOSTGAL_DDLR_SORT[i].LOGMASS_OBS = LOGMASS_OBS;
      
  }

  return ;

} // end GEN_SNHOST_LOGMASS


// =========================================    
void GEN_SNHOST_PROPERTY(int ivar_property) {

  // Created on Feb 2022 by M. Vincenzi and R. Kessler
  // For input ivar property, [property]_TRUE and [property]_ERR are used to generate 
  // [property]_OBS prop if [property]_OBS is not already in the HOSTLIB
  // If [property]_OBS is in the Hostlib, do nothing
  // [property] can be LOGMASS, LOGSFR, LOGsSFR, COLOR

  int  NNBR       = SNHOSTGAL.NNBR;
  int  IVAR_TRUE  = HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar_property].IVAR_TRUE ;
  int  IVAR_OBS   = HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar_property].IVAR_OBS;
  int  IVAR_ERR   = HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar_property].IVAR_ERR;
  double SCALE = HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar_property].SCALE_ERR;
  int i;
  double VAL_TRUE, VAL_OBS, VAL_ERR, GauRan ;
  double rmin=-3.0, rmax=3.0 ;
  char fnam[] = "GEN_SNHOST_PROPERTY" ;
  // ---------- BEGIN -----------  

  if ( IVAR_TRUE < 0 ) { return; }

  for(i=0; i < NNBR; i++ ) {

    VAL_OBS = -9.0 ;

    if ( IVAR_OBS > 0 ) {
      VAL_OBS = SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar_property].VAL_OBS ;
    }
    else if ( IVAR_TRUE > 0 && IVAR_ERR > 0 ) {
      VAL_TRUE = SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar_property].VAL_TRUE ;
      VAL_ERR  = SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar_property].VAL_ERR ;
      VAL_ERR *= SCALE ;
      SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar_property].VAL_ERR = VAL_ERR;
      GauRan = getRan_GaussClip(1,rmin,rmax);
      VAL_OBS = VAL_TRUE + GauRan*VAL_ERR ;
    }
    else {
      VAL_TRUE = SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar_property].VAL_TRUE ;
      VAL_OBS = VAL_TRUE;
    }

    SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar_property].VAL_OBS = VAL_OBS;
    //    printf("xxxxx %s: i=%d, ivar_property=%d, VAL_TRUE=%f VAL_OBS=%f VAL_ERR=%f\n", fnam,i,ivar_property, VAL_TRUE, VAL_OBS, VAL_ERR);
  }
  return ;
} // end GEN_SNHOST_PROPERTY

 

// =======================================
void GEN_SNHOST_POS(int IGAL) {

  // Mar 10, 2011
  // Fill SNHOSTGAL.[position]
  // Input IGAL is the sequential library index
  //
  // Galaxy profile is a sum of Sersic profiles.
  // First, randomly pick a profile based on its total weight (flux).
  // Next, pick a random location from this profile using 
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
  // Nov 15 2019: 
  //   + fix aweful bug for local a,b coords in ellipse frame.
  //   + call  GEN_SNHOST_ANGLE(IGAL,&phi);
  //
  // Apr 9 2020:
  //  + move IVAR_ANGLE<0 abort trap up to before get_Sersic_info()
  //    so that we get clean abort instead of seg fault.
  //
  // Jan 17 2022
  //   if RanInteg=0, set reduced_R=0 explicitly (e.g., AGN)
  //

  // strip off user options passed via sim-input file
  int LSN2GAL = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC ) ;
  int NPROF   =  SERSIC_PROFILE.NPROF ;
  int IVAR_RA     = HOSTLIB.IVAR_RA ;
  int IVAR_DEC    = HOSTLIB.IVAR_DEC ;
  int IVAR_ANGLE  = HOSTLIB.IVAR_ANGLE ;
  double RAD       = RADIAN ;

  int  j, JPROF, k_table, NBIN ;

  double 
    RA_GAL, DEC_GAL
    ,phi, cphi, sphi, crot, srot, top, bottom
    ,reduced_logR0, reduced_logR1, reduced_logR, reduced_R
    ,Ran0, Ran1, WGT, RanInteg, dif, bin, fbin
    ,a, b, a_half, b_half, a_rot, n, inv_n, DTMP, COSDEC, SNSEP, DLR    
    ,*ptr, *ptr_r, *ptr_integ0, *ptr_integ1
    ;

  int  DEBUG_MODE_SIMLIB = 0 ;
  int  LDMP = 0 ;
  char fnam[] = "GEN_SNHOST_POS" ;

  // -------------- BEGIN -------------

  SNHOSTGAL.phi               =  0.0 ;
  SNHOSTGAL.reduced_R         = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.SERSIC.INDEX      = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.RA_SNGALSEP_ASEC  = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.DEC_SNGALSEP_ASEC = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.SNSEP             = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.DLR               = HOSTLIB_SNPAR_UNDEFINED ;
  SNHOSTGAL.DDLR              = HOSTLIB_SNPAR_UNDEFINED ;

  // bail out if there are no galaxy shape parameters
  if ( NPROF == 0 ) { 

    if ( LSN2GAL ) {
      sprintf(c1err,"Cannot exec HOSTLIB_MSKOPT += %d", 
	      HOSTLIB_MSKOPT_SN2GAL_RADEC);
      sprintf(c2err,"Must define host Sersic profile to move SN near host.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    return ; 
  }

  if ( IVAR_ANGLE < 0 ) {
    sprintf(c1err,"Missing required %s in hostlib", HOSTLIB_VARNAME_ANGLE);
    sprintf(c2err,"Needed to choose position near host.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // extract info for each Sersic term
  get_Sersic_info(IGAL, &SNHOSTGAL.SERSIC) ;    

  a_rot  = SNHOSTGAL.SERSIC.a_rot ;
  crot   = cos(a_rot*RAD);
  srot   = sin(a_rot*RAD); 

  // for SIMLIB model, use already defined RA,DEC in SIMLIB header
  if ( INDEX_GENMODEL == MODEL_SIMLIB && !DEBUG_MODE_SIMLIB ) {
    SIMLIB_SNHOST_POS(IGAL, &SNHOSTGAL.SERSIC, 0 );
    
    // Aug 2021: hostlib should be the same as used to generated fakes;
    // but if not the SN-HOST separation is huge. If SN-host sep
    // is really big, then run through nominal process of shifting
    // the host to the SN
    bool CORRECT_HOST =  
      fabs(SNHOSTGAL.a_SNGALSEP_ASEC) < 5.0 &&
      fabs(SNHOSTGAL.b_SNGALSEP_ASEC) < 5.0 ;
    if ( CORRECT_HOST ) { DLR = 99999.0 ;    goto SNSEP_CALC ; }
  }

  // strip off random numbers to randomly generate a host-location
  Ran0  = SNHOSTGAL.FlatRan1_radius[0] ;
  Ran1  = SNHOSTGAL.FlatRan1_radius[1] ;

  // Start by randonly picking (Ran0) which Sersic profile 
  // based on the WGT of each profile.

  JPROF = -9;
  for ( j=0; j < NPROF; j++ ) {
    WGT = SNHOSTGAL.SERSIC.wsum[j];
    if ( WGT >= Ran0 && JPROF < 0 ) { JPROF = j ; }
  }


  // bail if we cannot pick a Sersic profile.
  if ( JPROF < 0 ) {
    ptr = SNHOSTGAL.SERSIC.wsum ; 
    sprintf(c1err,"Could not find random Sersic profile for Ran0=%f", Ran0);
    sprintf(c2err,"SERSIC_wsum = %f %f %f %f (GALID=%lld)",
	    ptr[0], ptr[1], ptr[2], ptr[3], SNHOSTGAL.GALID );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  a_half   = SNHOSTGAL.SERSIC.a[JPROF]; // half-light radius, major axis
  b_half   = SNHOSTGAL.SERSIC.b[JPROF]; // half-light radius, minor axis
  n        = SNHOSTGAL.SERSIC.n[JPROF]; 
  a_rot    = SNHOSTGAL.SERSIC.a_rot ;   // w.r.t RA, degrees
  crot = cos(a_rot*RAD);  srot = sin(a_rot*RAD); 


  GEN_SNHOST_ANGLE(a_half, b_half, &phi); // return phi angle
  cphi = cos(phi);  sphi = sin(phi);

  // Pick reduced radius (r=R/Rhalf),

  // get integral-table index for this Sersic_n
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
  //
  //G.Narayan - added min int flux to allow min offset from host center

  float MNINTFLUX = INPUTS.HOSTLIB_MNINTFLUX_SNPOS; 
  float MXINTFLUX = INPUTS.HOSTLIB_MXINTFLUX_SNPOS;
  RanInteg    = MNINTFLUX +  Ran1 * (MXINTFLUX - MNINTFLUX) ;

  if ( RanInteg > 0.0 ) {
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
  }
  else {
    reduced_R = 0.0 ; // Jan 17 2022 RK
  }


  /*
  printf(" ---------------------------------------------------- \n");
  printf(" kkkkkk = %d  1/n=%5.3f  dif=%6.3f  bin=%6.3f  fbin=%f \n", 
	 k_table, inv_n, dif, bin, fbin);
  printf(" xxx R0=%5.3f R1=%5.3f R(%d)=%5.3f   PHI=%5.1f deg. RanInteg=%f \n",
	 reduced_R0, reduced_R1, JPROF, reduced_R, phi/RAD, RanInteg );
  */

  // get major and minor half-light axes (arcsec) for this 
  // Sersic profile and this galaxy


  // Feb 4 2019: 
  // get DLR, distance to half-light ellipse, in direction of SN
  // See Eq 3.7 in Gupta hostmatching thesis
  top    = a_half * b_half ;
  bottom = sqrt(a_half*a_half*sphi*sphi + b_half*b_half*cphi*cphi ) ; 
  DLR    = top / bottom ;

  // get SN coords (arcsec) in ellipse-frame
  a = reduced_R * (DLR * cphi) ;   // bug fix, Nov 15 2019
  b = reduced_R * (DLR * sphi) ;

  if ( INPUTS.RESTORE_HOSTLIB_BUGS == true ) { // Nov 15 2019
    a = reduced_R * (a_half * cphi) ;  // restore bug
    b = reduced_R * (b_half * sphi) ;  // idem
  }

  SNHOSTGAL.a_SNGALSEP_ASEC  =  a ;
  SNHOSTGAL.b_SNGALSEP_ASEC  =  b ;
  SNHOSTGAL.phi              =  phi ;
  SNHOSTGAL.reduced_R        =  reduced_R ;
  SNHOSTGAL.SERSIC.INDEX     =  n ;
  SNHOSTGAL.DLR              =  DLR ;

  // now rotate coords based on major axis rotation angle w.r.t. RA

  // note that RA_SNGALSEP is angular sep, which is not the 
  // same as RA-difference.
  SNHOSTGAL.RA_SNGALSEP_ASEC  = (a * crot  +  b * srot) ;
  SNHOSTGAL.DEC_SNGALSEP_ASEC = (b * crot  -  a * srot) ;
  
  // ----------------------------------------------------
  // Determine Galaxy and SN coords

  if ( IVAR_RA >= 0 ) {
    // get absolute SN coords if we have host sky coords.
    // fetch galaxy coords from library
    RA_GAL   = HOSTLIB.VALUE_ZSORTED[IVAR_RA][IGAL] ;  // degrees
    DEC_GAL  = HOSTLIB.VALUE_ZSORTED[IVAR_DEC][IGAL] ; // degrees
    
    // compute absolute SN position relative to center of host
    COSDEC               = cos(DEC_GAL*RAD) ;
    DTMP                 = DEG_ARCSEC * SNHOSTGAL.RA_SNGALSEP_ASEC/COSDEC;
    SNHOSTGAL.RA_SN_DEG  = RA_GAL + DTMP ;
    
    DTMP                 = DEG_ARCSEC * SNHOSTGAL.DEC_SNGALSEP_ASEC ;
    SNHOSTGAL.DEC_SN_DEG = DEC_GAL + DTMP ;
  }
  else {
    // Feb 2019: no host coords, so fudge host position relative 
    // to already-determined SN pos. Goal is to get SNSEP and DLR.
    //
    SNHOSTGAL.RA_SN_DEG   = GENLC.RA ; // SN coord already selected
    SNHOSTGAL.DEC_SN_DEG  = GENLC.DEC ;

    COSDEC  = cos(SNHOSTGAL.DEC_SN_DEG*RAD) ; 
    DTMP    = DEG_ARCSEC * SNHOSTGAL.RA_SNGALSEP_ASEC / COSDEC ;
    RA_GAL  = SNHOSTGAL.RA_SN_DEG  - DTMP ;

    DTMP    = DEG_ARCSEC * SNHOSTGAL.DEC_SNGALSEP_ASEC ;
    DEC_GAL = SNHOSTGAL.DEC_SN_DEG - DTMP ;
  }

  // load gal coords into global
  SNHOSTGAL.RA_GAL_DEG  = RA_GAL ;
  SNHOSTGAL.DEC_GAL_DEG = DEC_GAL ;

 SNSEP_CALC:

  // compute SN-host separation in arcsec.
  SNSEP = angSep(SNHOSTGAL.RA_GAL_DEG, SNHOSTGAL.DEC_GAL_DEG,
		 SNHOSTGAL.RA_SN_DEG,  SNHOSTGAL.DEC_SN_DEG,
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

    // compute MWEBV with SN coords (was skipped in snlc_sim)
    // Note that return MWEBV is not used here.
    double MWEBV = gen_MWEBV(GENLC.RA, GENLC.DEC); 
  }

  // debug mode for SIMLIB model. Use forward-modeled RA,DEC as if
  // they were read from SIMLIB header ... then compare reverse-computed
  // a,b separations with those computed above.
  if ( INDEX_GENMODEL == MODEL_SIMLIB && DEBUG_MODE_SIMLIB ) {
    SIMLIB_HEADER.RA  = SNHOSTGAL.RA_SN_DEG ;
    SIMLIB_HEADER.DEC = SNHOSTGAL.DEC_SN_DEG ;
    SIMLIB_SNHOST_POS(IGAL, &SNHOSTGAL.SERSIC, 1 );
  }


  if ( LDMP ) {
    printf(" xxx %s: RA,DEC(SN) = %f, %f \n",
	   fnam, GENLC.RA, GENLC.DEC);
    printf(" xxx %s: RA,DEC(GAL)= %f, %f \n", 
	   fnam, SNHOSTGAL.RA_GAL_DEG, SNHOSTGAL.DEC_GAL_DEG );
    printf(" xxx %s: Ran0,Ran1 = %.4f, %.4f \n", fnam, Ran0, Ran1);
    printf(" xxx %s: a,b(half)=%.3f,%.3f n=%.3f a_rot=%.3f \n",
	   fnam, a_half, b_half, n, a_rot);
    printf(" xxx %s: a=%.4f b=%.4f  r=%.4f  phi=%.4f  DLR=%.4f \n",
	   fnam, a, b, reduced_R, phi, DLR);
    printf(" xxx %s: JPROF=%d IGAL=%d \n", fnam, JPROF, IGAL);
    fflush(stdout);
  }

  return ;

} // end of GEN_SNHOST_POS


// ================================================
void   GEN_SNHOST_ANGLE(double a, double b, double *ANGLE) {

  // Created Nov 18 2019
  // pick random ANGLE around ellipse.
  // [Fixes bug of random ANGLE between 0 and TWPI]
  //
  // Inputs: a, b = major and minor axis sizes
  // Ouptut: ANGLE (radians)
  //
  // Method 0: integrate r(ANGLE) and invert ... but I can't do the
  //   integral numerically.
  //
  // Method 1: pick random angle [0:2PI] and weight by DLR ...
  //   might be slow from cosine calculation each time, 
  //   or from many repeats with high eccentricities.
  //
  // Method 2: pick random point inside rectangle containing ellipse,
  //   and use rejection method to pick point inside ellipse. Then 
  //   tan(ANGLE)=  y/x. Should be ~3/4 efficient for any eccentricity.
  //  
  // Method 0 is optimal if there is an analytic form for the
  // integral. Here we go with method 2.

  double asq    = a*a;
  double bsq    = b*b;
  int    ilist  = 1 ;
  int    LEGACY = 0 ;
  int    LDMP   = (GENLC.CID == -9) ;
  //  double RAD    = RADIAN ;
  double fixran, phi, FlatRan_x, FlatRan_y, x, y, SUM ;
  char fnam[] = "GEN_SNHOST_ANGLE";

  // ------------------- BEGIN ---------------

  // check option to fix angle
  fixran = INPUTS.HOSTLIB_FIXRAN_PHI   ;
  if ( fixran > -1.0E-9 ) 
    { SNHOSTGAL.FlatRan1_phi = fixran ; LEGACY=1; }

  if ( LEGACY || INPUTS.RESTORE_HOSTLIB_BUGS ) 
    { *ANGLE  = (SNHOSTGAL.FlatRan1_phi * TWOPI) ; return; }


  if ( LDMP ) {  printf(" xxx ------------------------------- \n"); }
  
 PICK:

  // pick random point inside a rectangle containing ellipse
  FlatRan_x = getRan_Flat1(ilist) ;
  FlatRan_y = getRan_Flat1(ilist) ;
  x         = a * (2.0*FlatRan_x - 1.0); // -a to +a
  y         = b * (2.0*FlatRan_y - 1.0); // -b to +b

  // if not inside ellipse, try again.
  SUM = (x*x/asq) + (y*y/bsq);

  if ( LDMP ) {
    printf(" xxx %s: CID=%d x=%f, y=%f, SUM=%f \n", 
	   fnam, GENLC.CID, x, y, SUM); fflush(stdout);
  }

  if ( SUM > 1.0 ) { goto PICK; }

  if ( x == 0.000 ) { x = 1.0E-20; }
  phi  = atan2(y,x); 
  phi += PI;  // convert to 0 to 2PI range

  // load return ANGLE arg, and ignore radial component.
  *ANGLE = phi;

  return ;

} // end GEN_SNHOST_ANGLE


// ========================================
void GEN_SNHOST_NBR(int IGAL) {

  // Created Nov 2019 by R. Kessler
  // If NBR_LIST column exists, parse it and convert row numbers
  // to SNHOSTGAL.IGAL_NBR_LIST
  //
  // Jun 29 2021: change rownum-1 to rownum to fix index bug.

  int  LDMP = 0 ; // ( NCALL_GEN_SNHOST_DRIVER < 20 );
  int  i, ii, NNBR_READ, NNBR_STORE, rowNum, IGAL_STORE, IGAL_ZSORT ;
  int  ROWNUM_LIST[MXNBR_LIST];
  long long GALID ;
  char NBR_LIST[MXCHAR_NBR_LIST] ;
  char NO_NBR[] = "-1" ;
  char fnam[] = "GEN_SNHOST_NBR";

  // ---------------- BEGIN ----------------

  reset_SNHOSTGAL_DDLR_SORT(SNHOSTGAL.NNBR);

  SNHOSTGAL.NNBR = 1; // true host sets default at 1
  SNHOSTGAL.IGAL_NBR_LIST[0] = IGAL;

  // set default DDLR and SNSEP to that of true host
  SNHOSTGAL_DDLR_SORT[0].SNSEP = SNHOSTGAL.SNSEP ;
  SNHOSTGAL_DDLR_SORT[0].DDLR  = SNHOSTGAL.DDLR ;


  // bail if there is no NBR list
  if ( HOSTLIB.IVAR_NBR_LIST < 0 ) { return ; }

  sprintf(NBR_LIST, "%s", HOSTLIB.NBR_ZSORTED[IGAL] );
  GALID      = get_GALID_HOSTLIB(IGAL);

  // bail if this GAL has no NBR
  if ( strcmp(NBR_LIST,NO_NBR) == 0 ) { return ; }

  if ( LDMP ) {
    printf(" xxx ----------------------------------------- \n");
    printf(" xxx %s: parse NBR_LIST for GALID=%lld: \n", fnam, GALID ); 
  }

  if ( INPUTS.HOSTLIB_MAXREAD < MXROW_HOSTLIB ) {
    sprintf(c1err,"Cannot use HOSTLIB_MAXREAD with HOSTLIB that has NBR.");
    sprintf(c2err,"Remove HOSTLIB_MAXREAD or use HOSTLIB without NBR.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  // parse comma-sep list of HOSTLIB row numbers
  splitString2(NBR_LIST, COMMA, MXNBR_LIST , &NNBR_READ, &TMPWORD_HOSTLIB[1]);

  NNBR_READ++;    // include true host
  NNBR_STORE = 1; // start counter on stored neighbors

  ROWNUM_LIST[0] = -9;
  for(i=1; i < NNBR_READ; i++ ) {   // start at 1 to skip true host
    sscanf(TMPWORD_HOSTLIB[i], "%d", &rowNum);

    if ( rowNum > HOSTLIB.NGAL_READ ) {
      sprintf(c1err,"rowNum=%d exceeds NGAL_READ=%d",
	      rowNum, HOSTLIB.NGAL_READ);
      sprintf(c2err,"CID=%d  GALID=%lld NBR_LIST=%s", 
	      GENLC.CID, GALID, HOSTLIB.NBR_ZSORTED[IGAL]);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }
    if ( rowNum < 0 ) {
      sprintf(c1err,"Invalid rowNum=%d i=%d of %d", rowNum, i, NNBR_READ);
      sprintf(c2err,"CID=%d  GALID=%lld NBR_LIST=%s", 
	      GENLC.CID, GALID, HOSTLIB.NBR_ZSORTED[IGAL]);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }

    // Jun 2021 rowNum is no longer a fortran index due to other refactoring.
    // Use two layers of  indexing to get the desired z-sorted IGAL
    IGAL_STORE = HOSTLIB.LIBINDEX_READ[rowNum];
    if ( IGAL_STORE < 0 ) { continue; } // neighbor was cut from sample

    IGAL_ZSORT = HOSTLIB.LIBINDEX_ZSORT[IGAL_STORE]; 
    GALID      = get_GALID_HOSTLIB(IGAL_ZSORT);

    ii = NNBR_STORE; NNBR_STORE++ ;
    SNHOSTGAL.IGAL_NBR_LIST[ii] = IGAL_ZSORT;
    //HOSTLIB.IGAL_NBR_LIST = IGAL_ZSORT; // AG 09/2021 Need to fix later

    ROWNUM_LIST[ii] = rowNum; // for dump
    if( LDMP == 6 ) {
      printf("\t xxx rowNum(%d) = %d  --> GALID = %lld\n", i, rowNum, GALID); 
      fflush(stdout);
    }

  }
  
  SNHOSTGAL.NNBR = NNBR_STORE ;

  if ( LDMP ) {
    int  IVAR_RA     = HOSTLIB.IVAR_RA ;
    int  IVAR_DEC    = HOSTLIB.IVAR_DEC ;
    double RA_NBR, DEC_NBR, RA_REF, DEC_REF;
    RA_REF  = HOSTLIB.VALUE_ZSORTED[IVAR_RA][IGAL] ;  
    DEC_REF = HOSTLIB.VALUE_ZSORTED[IVAR_DEC][IGAL] ;  
    for(i=0; i < NNBR_STORE; i++ ) {
      IGAL_ZSORT = SNHOSTGAL.IGAL_NBR_LIST[i];
      GALID      = get_GALID_HOSTLIB(IGAL_ZSORT);
      RA_NBR     = HOSTLIB.VALUE_ZSORTED[IVAR_RA][IGAL_ZSORT] ;  
      DEC_NBR    = HOSTLIB.VALUE_ZSORTED[IVAR_DEC][IGAL_ZSORT] ; 
      printf(" xxx %2d : IGAL=%6d  ROW=%6d  GALID=%8lld  "
	     "Del(RA,DEC)=%7.4f,%7.4f \n",
	     i, IGAL_ZSORT, ROWNUM_LIST[i], GALID, 
	     RA_NBR-RA_REF, DEC_NBR-DEC_REF);
      fflush(stdout);
    }
  } // end LDMP

  return;

} // end GEN_SNHOST_NBR


// ==================================
void GEN_SNHOST_DDLR(int i_nbr) {

  // Created Nov 2019 by R.Kessler
  // Determine DLR for galaxy with sparse neighbor index i_nbr.
  // Note that i_nbr=0 corresponds to the true host for which
  // the SN was previously overlaid.
  //
  // If there are multiple Sersic terms, a & b are wgted average
  // among Sersic terms (Jan 2020)
  //
  // Nov 7 2020: Protect cosTH = 1 + tiny

  int IVAR_RA     = HOSTLIB.IVAR_RA ;
  int IVAR_DEC    = HOSTLIB.IVAR_DEC ;
  //  int IVAR_ANGLE  = HOSTLIB.IVAR_ANGLE ;
  double RAD      = RADIAN ;
  double ASEC_PER_DEG = 3600.0 ;
  double RA_SN    = SNHOSTGAL.RA_SN_DEG;  // deg
  double DEC_SN   = SNHOSTGAL.DEC_SN_DEG;
  
  double RA_GAL, DEC_GAL, DLR, DDLR, WTOT ;
  double SNSEP, top, bottom, a_rot, a_half, b_half, w,a,b ;
  int    IGAL, NPROF, j ;
  SERSIC_DEF SERSIC;
  char fnam[] = "GEN_SNHOST_DDLR" ;

  // -------------- BEGIN ------------

  SNHOSTGAL.DDLR_NBR_LIST[i_nbr]  = 0.0 ;
  SNHOSTGAL.SNSEP_NBR_LIST[i_nbr] = 0.0 ;

  // bail out if there are no galaxy shape parameters
  if ( SERSIC_PROFILE.NPROF == 0 ) { return ; }
  
  // get IGAL index to access full info.
  IGAL = SNHOSTGAL.IGAL_NBR_LIST[i_nbr] ;

  if ( IVAR_RA >= 0 ) {
    RA_GAL   = HOSTLIB.VALUE_ZSORTED[IVAR_RA][IGAL] ;   // deg
    DEC_GAL  = HOSTLIB.VALUE_ZSORTED[IVAR_DEC][IGAL] ; 
  }
  else {
    // just one host, 
    RA_GAL  =   SNHOSTGAL.RA_GAL_DEG ;
    DEC_GAL =   SNHOSTGAL.DEC_GAL_DEG ;
    if ( i_nbr > 0 ) { return ; }
  }
  // fetch Sersic profile info for this IGAL neighbor
  get_Sersic_info(IGAL, &SERSIC) ; 

  // compute SN-galaxy separation in arcsec.

  SNSEP = angSep(RA_GAL,DEC_GAL,  RA_SN,DEC_SN,  ASEC_PER_DEG);

  // For a & b, take weighted average among Sersic terms (Jan 2020).
  // This method is gut-feeling and not rigorous.
  NPROF = SERSIC.NPROF;
  WTOT  = SERSIC.wsum[NPROF-1] ;
  a_half = b_half = 0.0 ;
  for (j=0; j < NPROF; j++ ) {
    w = SERSIC.w[j];   a = SERSIC.a[j] ;   b = SERSIC.b[j];
    a_half += w*a;    b_half += w*b;
  }
  a_half /= WTOT;  b_half /= WTOT;

  // check option to use specialized a_DLR, b_DLR for DLR calculation (Helen Qu, 1/25/2022)
  if (HOSTLIB.IVAR_a_DLR > 0) {
    a_half = HOSTLIB.VALUE_ZSORTED[HOSTLIB.IVAR_a_DLR][IGAL]; 
    b_half = HOSTLIB.VALUE_ZSORTED[HOSTLIB.IVAR_b_DLR][IGAL]; 
  }

  a_rot    = SERSIC.a_rot ;   // rot angle (deg) w.r.t. RA

  // for DLR calc, move to frame where RA=DEC=0 for galaxy center
  // so that RA,DEC can be treated as cartesian coordinates.

  double VEC_aHALF[2], VEC_SN[2], DOTPROD, LEN_SN ;  
  double cosTH, sqcos, sqsin;

  VEC_aHALF[0] = +a_half * cos(RAD*a_rot) ;
  VEC_aHALF[1] = -a_half * sin(RAD*a_rot) ;
  VEC_SN[0]    = RA_SN  - RA_GAL  ;
  VEC_SN[1]    = DEC_SN - DEC_GAL ; // i.e., DEC=0 for galaxy center
  LEN_SN       = sqrt( VEC_SN[0]*VEC_SN[0] + VEC_SN[1]*VEC_SN[1] ) ;

  DOTPROD = VEC_aHALF[0]*VEC_SN[0] + VEC_aHALF[1]*VEC_SN[1];
  cosTH   = DOTPROD/(LEN_SN*a_half);

  if ( fabs(cosTH) > 1.000001 ) {
    sprintf(c1err,"Invalid cosTH = %le for GALID=%lld", 
	    cosTH, SNHOSTGAL.GALID);
    sprintf(c2err,"LEN_SN = %f, a_half=%f, DOT=%f ",
	    LEN_SN, a_half, DOTPROD);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  
  if (cosTH >  1.00) { cosTH = +1.0 - 1.0E-9 ; }
  if (cosTH < -1.00) { cosTH = -1.0 + 1.0E-9 ; }

  sqcos  = cosTH*cosTH;   sqsin = 1.0 - sqcos;
  top    = ( a_half * b_half ) ;
  bottom = sqrt(a_half*a_half*sqsin + b_half*b_half*sqcos ) ;   
  DLR    = top/bottom; 

  if ( DLR > 0.0 ) 
    {  DDLR = SNSEP/DLR ; }
  else
    { DDLR = 0.0 ; }

  // store in global to be analyzed later
  SNHOSTGAL.DDLR_NBR_LIST[i_nbr]  = DDLR;
  SNHOSTGAL.SNSEP_NBR_LIST[i_nbr] = SNSEP ;

  return ;

} // end GEN_SNHOST_DDLR


// ==============================
void reset_SNHOSTGAL_DDLR_SORT(int MAXNBR) {

  SNHOSTGAL.NNBR = 0;
  int i, ifilt, q;
  for(i=0; i < MAXNBR; i++ ) {    
    SNHOSTGAL_DDLR_SORT[i].GALID = -9 ;
    SNHOSTGAL_DDLR_SORT[i].SNSEP = -9.0 ;
    SNHOSTGAL_DDLR_SORT[i].DDLR  = -9.0 ;  
    SNHOSTGAL_DDLR_SORT[i].RA    = 999.0 ;
    SNHOSTGAL_DDLR_SORT[i].DEC   = 999.0 ;
    SNHOSTGAL_DDLR_SORT[i].TRUE_MATCH = false ;
    SNHOSTGAL_DDLR_SORT[i].GALID2       = -9;
    SNHOSTGAL_DDLR_SORT[i].GALID_UNIQUE = -9;
    SNHOSTGAL_DDLR_SORT[i].ELLIPTICITY  = 999.0;
    SNHOSTGAL_DDLR_SORT[i].SQRADIUS = 999.0;
    for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
      SNHOSTGAL_DDLR_SORT[i].MAG[ifilt]      = -9.0 ;
      SNHOSTGAL_DDLR_SORT[i].MAG_ERR[ifilt]  = -9.0 ;
    }
    for(q=0; q < MXBIN_ZPHOT_Q; q++){
      SNHOSTGAL_DDLR_SORT[i].ZPHOT_Q[q] = -9.0;
    }
  }

  SNHOSTGAL.IMATCH_TRUE = -9 ;   // May 31 2021

} // end reset_SNHOSTGAL_DDLR_SORT

// =======================================================
void SIMLIB_SNHOST_POS(int IGAL, SERSIC_DEF *SERSIC, int DEBUG_MODE) {

  // For SIMLIB model (using MAG column of SIMLIB file),
  // load RA and DEC coordinates in SNHOSTGAL struct,
  // and compute coordinates in a,b frame that is rotated 
  // w.r.t. RA,DEC frame.
  //
  // DEBUG_MODE > 0 --> compare a_sep,b_sep with already computed values.

  int IVAR_RA     = HOSTLIB.IVAR_RA ;
  int IVAR_DEC    = HOSTLIB.IVAR_DEC ;
  double RAD      = RADIAN ;
  double a_rot    = SERSIC->a_rot;
  double crot     = cos(RAD*a_rot);
  double srot     = sin(RAD*a_rot);
  double RA_GAL, DEC_GAL, RA_SN, DEC_SN, RA_SEP, DEC_SEP, COSDEC;
  double a_sep, b_sep;
  char fnam[] = "SIMLIB_SNHOST_POS";

  // ------------- BEGIN ------------

  RA_GAL   = HOSTLIB.VALUE_ZSORTED[IVAR_RA][IGAL] ; 
  DEC_GAL  = HOSTLIB.VALUE_ZSORTED[IVAR_DEC][IGAL] ;
  COSDEC   = cos(DEC_GAL*RAD);

  RA_SN    = SIMLIB_HEADER.RA;
  DEC_SN   = SIMLIB_HEADER.DEC;
  RA_SEP   = 3600.0 * (RA_SN - RA_GAL) * COSDEC ; // angSep, not RA-diff
  DEC_SEP  = 3600.0 * (DEC_SN - DEC_GAL) ;

  SNHOSTGAL.RA_GAL_DEG    = RA_GAL; 
  SNHOSTGAL.DEC_GAL_DEG   = DEC_GAL; 
  SNHOSTGAL.RA_SN_DEG     = RA_SN; 
  SNHOSTGAL.DEC_SN_DEG    = DEC_SN; 
  
  SNHOSTGAL.RA_SNGALSEP_ASEC  = RA_SEP ;
  SNHOSTGAL.DEC_SNGALSEP_ASEC = DEC_SEP ;
  
  // compute SN-GAL separations along a & b axes: needed for GALMAG calc.
  a_sep = RA_SEP  * crot - DEC_SEP * srot ; 
  b_sep = DEC_SEP * crot + RA_SEP  * srot ; 

  if ( DEBUG_MODE == 0 ) {
      SNHOSTGAL.a_SNGALSEP_ASEC = a_sep ;
      SNHOSTGAL.b_SNGALSEP_ASEC = b_sep ;
  }
  else {
    long long GALID = get_GALID_HOSTLIB(IGAL);
    printf(" xxx ------------------------------------------------ \n" ) ;
    printf(" xxx %s  DEBUB  DUMP for GALID = %lld\n", fnam, GALID );
    printf(" xxx RA,DEC = %.3f, %.3f \n", RA_GAL, DEC_GAL);
    printf(" xxx Sersic a, b = %.3f, %.3f   a_rot=%.1f deg \n", 
	   SERSIC->a[0], SERSIC->b[0], a_rot );
    printf(" xxx SEP(RA,DEC): %8.4f , %8.4f \n", RA_SEP, DEC_SEP);
    printf(" xxx original-fwd a,b = %8.4f , %8.4f \n", 
	   SNHOSTGAL.a_SNGALSEP_ASEC, SNHOSTGAL.b_SNGALSEP_ASEC);
    printf(" xxx SIMLIB-test  a,b = %8.4f , %8.4f \n", a_sep, b_sep);

    printf(" xxx difference   a,b = %8.4f , %8.4f   (DEC=%7.2f)\n",
	   a_sep - SNHOSTGAL.a_SNGALSEP_ASEC, 
	   b_sep - SNHOSTGAL.b_SNGALSEP_ASEC, DEC_GAL );
    fflush(stdout);
  }
 
  return ;

} // end SIMLIB_SNHOST_POS

// =================================
void SORT_SNHOST_byDDLR(void) {

  // created Nov 2019
  // Sort galaxy NBRs by DDLR, and load global structure
  // SNHOSTGAL_DDLR_SORT[i], where i=0 has smallest DDLR.
  //
  // At end of function, set SNHOSTGAL.NNBR = number passing DDLR cut.
  //
  // May 20 2020: bug fix for LSN2GAL
  // Oct 25 2021: compute optional GALID_UNIQUE (for LSST broker test)
  // Nov 17 2021: correct host mags by DMUCOR = MU(zSN) - MU(zGAL)
  // Jan 22 2022: set GENLC.CORRECT_HOSTMATCH=False for wrong host match

  int  MSKOPT           = INPUTS.HOSTLIB_MSKOPT ;
  bool LSN2GAL_Z        = (MSKOPT & HOSTLIB_MSKOPT_SN2GAL_Z) ;
  bool LSN2GAL_RADEC    = (MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC);
  int  USE_QGAUSS      =  (MSKOPT & HOSTLIB_MSKOPT_ZPHOT_QGAUSS);
  int  NNBR             = SNHOSTGAL.NNBR ;
  int  IVAR_RA          = HOSTLIB.IVAR_RA;
  int  IVAR_DEC         = HOSTLIB.IVAR_DEC ;
  int  IVAR_ZPHOT       = HOSTLIB.IVAR_ZPHOT; 
  int  IVAR_ZPHOT_ERR   = HOSTLIB.IVAR_ZPHOT_ERR; 
  int  IVAR_Q0          = HOSTLIB.IVAR_ZPHOT_Q0 ;
  int  ORDER_SORT       = +1 ;     // increasing order
  int  LDMP = 0 ; // (GENLC.CID == 9 ) ;

  int  INDEX_UNSORT[MXNBR_LIST];
  int  i, unsort, IGAL, IVAR, IVAR_ERR, ifilt, ifilt_obs, q ,ivar, IVAR_Q;
  int  NNBR_DDLRCUT = 0 ;
  double DDLR, SNSEP, MAG, MAG_ERR, RA_GAL, DEC_GAL, zq ;
  double DMUCOR = 0.0 ;
  char fnam[] = "SORT_SNHOST_byDDLR" ;

  // ------------- BEGIN ---------------

  // sort by DDLR
  sortDouble( NNBR, SNHOSTGAL.DDLR_NBR_LIST, ORDER_SORT, INDEX_UNSORT ) ;

  // check for SN-host distance-diff correction on host gal mags (Nov 2021)
  ifilt_obs = GENLC.IFILTMAP_OBS[0];
  IVAR      = HOSTLIB.IVAR_MAGOBS[ifilt_obs] ;
  if ( IVAR >= 0 && !LSN2GAL_Z ) {
    double HOST_DLMU, LENSDMU, zCMB, zHEL;
    zHEL = SNHOSTGAL.ZTRUE; 
    zCMB = zhelio_zcmb_translator(zHEL, GENLC.RA, GENLC.DEC, "eq",+1);
    gen_distanceMag(zCMB, zHEL,
		    &HOST_DLMU, &LENSDMU ); // <== returned
    DMUCOR = GENLC.DLMU - HOST_DLMU ; // ignore LENSDMU that cancels

    /* 
    printf(" xxx %s: DMUCOR = %.4f(zSN=%.4f) - %.4f(zHOST=%.4f) = %.4f\n",
	   fnam, GENLC.DLMU, GENLC.REDSHIFT_CMB, HOST_DLMU, zCMB, DMUCOR);  
    */
  }

  //  LDMP = ( INDEX_SORT[0] > 0 ) ;

  if ( LDMP ) 
    { printf(" xxx ----------------------------- \n"); }

  // load info sorted by DDLR
  for(i=0; i < NNBR; i++ ) {
    unsort = INDEX_UNSORT[i];
    IGAL   = SNHOSTGAL.IGAL_NBR_LIST[unsort] ;
    DDLR   = SNHOSTGAL.DDLR_NBR_LIST[unsort] ;
    SNSEP  = SNHOSTGAL.SNSEP_NBR_LIST[unsort] ;

    if ( DDLR < INPUTS.HOSTLIB_MAXDDLR ) { NNBR_DDLRCUT++ ; }

    RA_GAL = GENLC.RA;    DEC_GAL = GENLC.DEC;
    if ( i == 0 ) { 
      RA_GAL  += (SNHOSTGAL.RA_GAL_DEG  - SNHOSTGAL.RA_SN_DEG  ) ;
      DEC_GAL += (SNHOSTGAL.DEC_GAL_DEG - SNHOSTGAL.DEC_SN_DEG ) ;
    }

    // load logical for true host
    SNHOSTGAL_DDLR_SORT[i].TRUE_MATCH = false ;
    if ( unsort == 0 ) { // first element of unsorted array is true host
       SNHOSTGAL_DDLR_SORT[i].TRUE_MATCH = true ; 
       SNHOSTGAL.IMATCH_TRUE = i;
       if ( i != 0 ) { GENLC.CORRECT_HOSTMATCH = false; } // Jan 20 2022
       if ( LDMP ) { printf("\t xxx %s: IMATCH_TRUE = %d \n", fnam, i); }
    }

    // load global struct
    SNHOSTGAL_DDLR_SORT[i].DDLR  = DDLR ;
    SNHOSTGAL_DDLR_SORT[i].SNSEP = SNSEP ;
    SNHOSTGAL_DDLR_SORT[i].GALID = get_GALID_HOSTLIB(IGAL);

    // if HOSTLIB coords don't match the SN, then use GAL-SN difference
    // to determine final host coords near SN. This feature allows using
    // HOSTLIB with any set of coordinates, even coords well outside
    // SN fields. If SN->GAL option, then don't change galaxy position.
    SNHOSTGAL_DDLR_SORT[i].RA  = RA_GAL;
    SNHOSTGAL_DDLR_SORT[i].DEC = DEC_GAL ;

    if ( IVAR_RA>0 && IVAR_DEC>0 &&  !LSN2GAL_RADEC ) {
	SNHOSTGAL_DDLR_SORT[i].RA += 
	  ( get_VALUE_HOSTLIB(IVAR_RA,IGAL) - SNHOSTGAL.RA_SN_DEG );

	SNHOSTGAL_DDLR_SORT[i].DEC += 
	  ( get_VALUE_HOSTLIB(IVAR_DEC,IGAL) - SNHOSTGAL.DEC_SN_DEG );
    }

    IVAR = HOSTLIB.IVAR_ZTRUE ;
    SNHOSTGAL_DDLR_SORT[i].ZSPEC = get_VALUE_HOSTLIB(IVAR,IGAL);
    SNHOSTGAL_DDLR_SORT[i].ZSPEC_ERR = 0.0005; // any small value

    if ( IVAR_ZPHOT > 0 ) {
      SNHOSTGAL_DDLR_SORT[i].ZPHOT     = get_VALUE_HOSTLIB(IVAR_ZPHOT,IGAL); 
      SNHOSTGAL_DDLR_SORT[i].ZPHOT_ERR = get_VALUE_HOSTLIB(IVAR_ZPHOT_ERR,IGAL);
    }
    else {
      SNHOSTGAL_DDLR_SORT[i].ZPHOT     = -9.0 ;
      SNHOSTGAL_DDLR_SORT[i].ZPHOT_ERR = -9.0 ;
    }
    
    // - - - - - - -
    if ( IVAR_Q0 > 0 || USE_QGAUSS ) {
      for (q = 0; q < HOSTLIB.NZPHOT_Q; q++) {
	zq = GEN_SNHOST_ZPHOT_QUANTILE(IGAL,q);
	SNHOSTGAL_DDLR_SORT[i].ZPHOT_Q[q] = zq; 
      }
    }

    for (ivar=0; ivar<N_HOSTGAL_PROPERTY; ivar++){
      IVAR = HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar].IVAR_TRUE;
      if ( IVAR > 0 ) {SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar].VAL_TRUE = get_VALUE_HOSTLIB(IVAR,IGAL);}

      IVAR = HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar].IVAR_OBS;
      if ( IVAR > 0 ) {SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar].VAL_OBS = get_VALUE_HOSTLIB(IVAR,IGAL);}
      
      IVAR = HOSTLIB.HOSTGAL_PROPERTY_IVAR[ivar].IVAR_ERR;
      if ( IVAR > 0 ) {SNHOSTGAL_DDLR_SORT[i].HOSTGAL_PROPERTY_VALUE[ivar].VAL_ERR = get_VALUE_HOSTLIB(IVAR,IGAL);}
    }

    // - - - - -
    IVAR = HOSTLIB.IVAR_GALID2;
    if ( IVAR > 0 ) 
      { SNHOSTGAL_DDLR_SORT[i].GALID2 = get_VALUE_HOSTLIB(IVAR,IGAL); } 

    if ( INPUTS.HOSTLIB_GALID_UNIQUE ) {
      // compute GALID_UNIQUE from GALID (it's NOT read from hostlib)
      set_GALID_UNIQUE(i);
    }

    IVAR = HOSTLIB.IVAR_ELLIPTICITY;
    if ( IVAR > 0 )
      { SNHOSTGAL_DDLR_SORT[i].ELLIPTICITY = get_VALUE_HOSTLIB(IVAR,IGAL); }

    IVAR = HOSTLIB.IVAR_SQRADIUS;
    if ( IVAR > 0 )
      { SNHOSTGAL_DDLR_SORT[i].SQRADIUS = get_VALUE_HOSTLIB(IVAR,IGAL); }

    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      IVAR      = HOSTLIB.IVAR_MAGOBS[ifilt_obs] ;
      if ( IVAR > 0 ) {
	MAG       = get_VALUE_HOSTLIB(IVAR,IGAL) ;
	SNHOSTGAL_DDLR_SORT[i].MAG[ifilt_obs] = MAG + DMUCOR ; 
      }

      IVAR_ERR      = HOSTLIB.IVAR_MAGOBS_ERR[ifilt_obs] ;
      if ( IVAR_ERR > 0 ) {
      	MAG_ERR       = get_VALUE_HOSTLIB(IVAR_ERR,IGAL) ;
     	SNHOSTGAL_DDLR_SORT[i].MAG_ERR[ifilt_obs] = MAG_ERR ;
      }

      if ( INPUTS.HOSTLIB_GALID_UNIQUE && IVAR_ERR > 0 ) {
	// for unique GALIDs, fluctuate mags to avoid duplicate mags
	double GauRan, rmin=-3., rmax=3.;
	MAG_ERR       = get_VALUE_HOSTLIB(IVAR_ERR,IGAL) ;
	GauRan = getRan_GaussClip(1,rmin,rmax);
	SNHOSTGAL_DDLR_SORT[i].MAG[ifilt_obs] += MAG_ERR*GauRan;
      }
    }

    if ( LDMP ) { 
      printf("\t xxx %s: i=%d unsort=%d  NNBR_DDLRCUT=%d \n",
	     fnam, i, unsort, NNBR_DDLRCUT);
      printf("\t xxx %s: DDLR=%6.2f  SEP=%8.1f\n",
	     fnam,SNHOSTGAL_DDLR_SORT[i].DDLR, SNHOSTGAL_DDLR_SORT[i].SNSEP);
	     // xxx DDLR, SNSEP ); 
      printf("\t xxx %s: RA_GAL=%f DEC_GAL=%f \n",
	     fnam, SNHOSTGAL_DDLR_SORT[i].RA, SNHOSTGAL_DDLR_SORT[i].DEC);
      printf("\t xxx %s: LOGMASS = %f \n",
	     fnam, SNHOSTGAL_DDLR_SORT[i].LOGMASS_TRUE );
      fflush(stdout);      
    }
  }
  
  // truncate list to those satisfying DDLR cut (Feb 2020)
  SNHOSTGAL.NNBR = NNBR_DDLRCUT ;

  return ;

} // SORT_SNHOST_byDDLR

// ==========================================
double GEN_SNHOST_ZPHOT_QUANTILE(int IGAL, int q) {

  // Created Apr 19 2022
  // For sparse zphot-quantile index q, return zPHOT value.
  // Default method is to return value from HOSTLIB.
  // If HOSTLIB_MSKOPT_ZPHOT_QGAUSS bit of HOSTLIB_MSKOPT is set,
  // then compute quantile from Gaussian using sigma=ZPHOT_ERR .

  int   USE_QGAUSS = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_ZPHOT_QGAUSS );
  int   IVAR_Q0          = HOSTLIB.IVAR_ZPHOT_Q0 ;
  int   IVAR_ZPHOT       = HOSTLIB.IVAR_ZPHOT; 
  int   IVAR_ZPHOT_ERR   = HOSTLIB.IVAR_ZPHOT_ERR; 
  int   IVAR_Q;
  double ZPHOT, ZPHOT_ERR, SIGMA_QGAUSS;
  double zq = 0.0 ;
  char fnam[] = "GEN_SNHOST_ZPHOT_QUANTILE";

  // ------ BEGIN -------

  if ( USE_QGAUSS ) {
    // test feature to use Gauss approx 
    ZPHOT     = get_VALUE_HOSTLIB(IVAR_ZPHOT,    IGAL);
    ZPHOT_ERR = get_VALUE_HOSTLIB(IVAR_ZPHOT_ERR,IGAL);

    SIGMA_QGAUSS = HOSTLIB.SIGMA_QGAUSS[q]; // Nsigma from mean at this quantile
    zq           = ZPHOT + ZPHOT_ERR * SIGMA_QGAUSS;

    /* 
    printf(" xxx %s: q=%d SIGMA_QGAUSS=%.4f  ZPHOT=%.3f +_ %.3f  zq=%.3f \n",
	   fnam, q, SIGMA_QGAUSS, ZPHOT, ZPHOT_ERR, zq); fflush(stdout);
    //    debugexit(fnam);
    */
  } 
  else {
    // default quantiles from HOSTLIB
    IVAR_Q = IVAR_Q0 + q;
    zq  = get_VALUE_HOSTLIB(IVAR_Q,IGAL) ;
  }

  zq += SNHOSTGAL.ZDIF; // shift is  zSN - zGAL

  return zq;

} // end GEN_SNHOST_ZPHOT_QUANTILE

// ================================= 
void set_GALID_UNIQUE(int i){
  // Created Oct 2021 by A.Gagliaon and R.Kessler
  //
  // Assign unique GALID in separate GALID_UNIQUE variable.
  // Initial use is LSST DESC broker test where HOSTLIB galaxies are
  // re-used many times, but a unique GALID is needed for each event.
  // Algorithm:
  //   GALID_UNIQUE = CID * CID_MULTIPLIER + ran[0,CID_MULTIPLIER]
  // and works only if the CID are unique, which is an already
  // existing feature (FORMAT_MASK += 16).
  //
  // Input i is the neighbor index 
  //

  int CID = GENLC.CID;
  int CID_MULTIPLIER = 57;
  long long GALID_UNIQUE;
  int i2, ilist = 1;
  double rand1 ;
  bool match_GALID = true ;
  bool test_DUPLICATE_AVOID = false ; // for unit test only
  char fnam[] = "set_GALID_UNIQUE";

  // ----------- BEGIN ------------

  if ( test_DUPLICATE_AVOID ) { CID_MULTIPLIER = 5; }

  while ( match_GALID ) {
    rand1        = getRan_Flat1(ilist);
    GALID_UNIQUE = CID*CID_MULTIPLIER + (int)(rand1*CID_MULTIPLIER);
    SNHOSTGAL_DDLR_SORT[i].GALID_UNIQUE = GALID_UNIQUE;

    // done for first host on list.
    if ( i == 0 ) { return ; }

    // check if GALID_UNIQUE matches previous neighbor ;
    // if there's a match, pick another GALID
    match_GALID = false;
    for (i2 = 0; i2 < i; i2++){
      if (GALID_UNIQUE == SNHOSTGAL_DDLR_SORT[i2].GALID_UNIQUE){
	match_GALID = true;
	if ( test_DUPLICATE_AVOID ) {
	  printf("xxx %s: DUPLICATE GALID FOR CID = %d, DDLR = %.2f\n", 
		 fnam, CID, SNHOSTGAL_DDLR_SORT[i].DDLR);
	}
      }
    } //end i2 over previous neighbors

  } // end while

} // end set_GALID_UNIQUE

// =================================
void TRANSFER_SNHOST_REDSHIFT(int IGAL) {

  // Apr 2019
  // Check for redshift change; either explicit transfer to use
  // zHOST, or switch to host coords and change zHEL.
  // Add vpec here, based on sim-input VPEC or HOSTLIB VPEC.
  //
  // Apr 20 2019: bail on incorrect host match.
  // May 25 2020: add vpec

  double ZTRUE         = SNHOSTGAL.ZTRUE ;     // helio redshift
  double zPEC_GAUSIG   = GENLC.VPEC/LIGHT_km ; // from GENSIGMA_VPEC key
  double zPEC_HOSTLIB  = SNHOSTGAL.VPEC/LIGHT_km; // from hostlib

  double RA            = GENLC.RA ;
  double DEC           = GENLC.DEC ;
 
  int  MSKOPT          = INPUTS.HOSTLIB_MSKOPT ;
  bool DO_VPEC         = (MSKOPT & HOSTLIB_MSKOPT_USEVPEC ) ;
  bool DO_SN2GAL_Z     = (MSKOPT & HOSTLIB_MSKOPT_SN2GAL_Z) ;
  bool DO_SN2GAL_RADEC = (MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC) ;

  double zCMB, zHEL, zPEC ;
  char eq[]           = "eq";
  char fnam[]         = "TRANSFER_SNHOST_REDSHIFT" ;
  int LDMP = 0; // ( GENLC.CID < -1010 ) ;

  // ------------ BEGIN ------------

  // if wrong host (based on mag), bail
  if ( !GENLC.CORRECT_HOSTMATCH) { return ; }

  // un-do zPEC to get zHEL without vPEC

  zHEL = (1.0+GENLC.REDSHIFT_HELIO)/(1.0+zPEC_GAUSIG) - 1.0 ;

  // - - - - - - - - - - - - - - - - - - - - - 
  // check for transferring SN redshift to host redshift.
  // Here zHEL & zCMB both change
  if ( DO_SN2GAL_Z ) {
    zHEL = ZTRUE ;  // true host z
    if ( INPUTS.VEL_CMBAPEX > 1.0 ) {
      zCMB = zhelio_zcmb_translator(zHEL,RA,DEC,eq,+1);
    }
    else {
      zCMB = zHEL ; 
    }

    GENLC.REDSHIFT_CMB   = zCMB ;   // store adjusted zCMB
    gen_distanceMag(zCMB, zHEL,
		    &GENLC.DLMU, &GENLC.LENSDMU ); // <== returned
  }

  // - - - - - - - - - - - - - - - - - - - - - 
  // check for switching to host coordinates; 
  // zCMB does not change, but zHEL changes.

  if ( DO_SN2GAL_RADEC && !DO_SN2GAL_Z ) {
    zCMB = GENLC.REDSHIFT_CMB  ; // preserve this
    if ( INPUTS.VEL_CMBAPEX > 1.0 ) 
      { zHEL = zhelio_zcmb_translator(zCMB,RA,DEC,eq,-1); }   
    else 
      { zHEL = zCMB; }

    gen_distanceMag(zCMB, zHEL, 
		    &GENLC.DLMU, &GENLC.LENSDMU ); // <== returned

  }

  // - - - - - - - - - 
  // May 25 2020 add zPEC to zHEL
  if ( DO_VPEC ) 
    { zPEC = zPEC_HOSTLIB; }  // from VPEC key in HOSTLIB
  else
    { zPEC = zPEC_GAUSIG ; }  // from sim-input GENSIGMA_VPEC key

  double zHEL_ORIG = GENLC.REDSHIFT_HELIO ;
  GENLC.VPEC = zPEC * LIGHT_km;
  zHEL = (1.0+zHEL)*(1.0+zPEC) - 1.0 ;    // Jan 27 2021 

  GENLC.REDSHIFT_HELIO  = zHEL ;     
  SNHOSTGAL.ZSPEC       = zHEL ;  

  if ( LDMP ) {
    printf(" xxx -------------------------------- \n");
    printf(" xxx %s  DUMP for CID = %d \n", fnam, GENLC.CID);
    printf(" xxx DO[VPEC,SN2GAL_Z,SN2GAL_RADEC = %d %d %d \n",
	   DO_VPEC, DO_SN2GAL_Z, DO_SN2GAL_RADEC );
    printf(" xxx zPEC(Gauss,HOSTLIB->FINAL) = %f, %f -> %f\n",
	   zPEC_GAUSIG, zPEC_HOSTLIB, zPEC );
    printf(" xxx zHEL = %f -> %f \n", zHEL_ORIG, zHEL );
    fflush(stdout);
  }

  return ;

} // end of TRANSFER_SNHOST_REDSHIFT

// ==============================
void GEN_SNHOST_GALMAG(int IGAL) {

  // Compute TRUE host mag contained within PSF for each observer-filter,
  // which is used to compute local surface brightness.
  // Assume that each [filt]_obs value corresponds to the total
  // flux of the host; then use the generated SN position to 
  // determine several apertures and compute the flux-fraction 
  // in each aperture. Finally, compute total host mag for each neighbor.
  //
  // BEWARE: for host flux under SN, only true host is used.
  //         TO be exactly right, should sum all neighbors.
  //
  // For a given PSF, the effective aperture is
  // taken to be 4*PI*PSF^2.
  //
  // Fill struct
  //   SNHOSTGAL_DDLR_SORT[inbr].MAG[ifilt_obs] = MAG ; 
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
  // Nov 25 2019: protect dm for GALFRAC=0
  //
  // Jan 31 2020: refactor to load DDLR_SORT array for MAG.

  int  NNBR       = SNHOSTGAL.NNBR ;

  double 
     x_SN, y_SN
    ,xgal, ygal, MAGOBS, MAGOBS_LIB, PSF, FGAL, dF
    ,Rmin, Rmax, Rbin, R, Rcen, dm
    ,THmin, THmax, THbin, TH
    ,dRdTH, Jac
    ,GALFRAC_SUM[NMAGPSF_HOSTLIB+1]       // summed over Sersic profile
    ,sigFrac, RcenFrac, GALFRAC
    ,GaussOvp[NMAGPSF_HOSTLIB+1] 
    ,AV, LAMOBS_AVG, MWXT[MXFILTINDX]
    ,RVMW = 3.1
    ;

  float lamavg4, lamrms4, lammin4, lammax4  ;
  int ifilt, ifilt_obs, i, inbr, IVAR, jbinTH, opt_frame  ;
  char cfilt[2];
  char fnam[] = "GEN_SNHOST_GALMAG" ;

  // ------------ BEGIN -------------

  for ( i=0; i <= NMAGPSF_HOSTLIB ; i++ ) { GALFRAC_SUM[i] =  0.0 ; }
  

  // compute MilkyWay Galactic extinction in each filter using
  // extinction at mean wavelength.
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
    if ( IVAR < 0 ) { continue; }

    // dim each possible host/neighbor (Jan 31 2020)
    for(inbr=0; inbr < NNBR; inbr++ ) {     
      SNHOSTGAL_DDLR_SORT[inbr].MAG[ifilt_obs] += MWXT[ifilt_obs] ;       
    }

    /*
      printf(" xxx ------------------------ \n");
      printf(" xxx %s: load ifilt_obs=%d IGAL=%d, IVAR=%d \n",
      fnam, ifilt_obs, IGAL, IVAR );
      printf(" xxx %s: MWXT=%.3f   MAGOBS=%.3f -> %.3f \n" ,
      fnam, MWXT[ifilt_obs], MAGOBS_LIB, MAGOBS);
    */
       
  } // end ifilt

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

  if ( x_SN == HOSTLIB_SNPAR_UNDEFINED || y_SN == HOSTLIB_SNPAR_UNDEFINED ) {
    sprintf(c1err,"Undefined SNGALSEP");
    sprintf(c2err,"SNGALSEP(a,b) = %f,%f arcSec", x_SN, y_SN);	
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }


  // max radius (w.r.t. SN) to integrate is twice the 
  // aperture radius,  or 2 x (2*PSFMAX).
  Rmin  = 0.0 ;
  Rmax  = HOSTLIB.Aperture_Rmax ;
  Rbin  = HOSTLIB.Aperture_Rbin ;

  THmin = 0.0 ;
  THmax = TWOPI * .99999 ;
  THbin = HOSTLIB.Aperture_THbin ;

  dRdTH = Rbin * THbin ;

  /* xxxxxx
  long long GALID  = get_GALID_HOSTLIB(IGAL);
  int j=0;
  printf(" xxx %s:  GALID=%lld   a,b(%d) = %.3f , %.3f \n",
	 fnam, GALID, j, SNHOSTGAL.SERSIC.a[j], SNHOSTGAL.SERSIC.b[j] );
  fflush(stdout);
  xxxxxxxxxx */


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
  // Note that here we use TRUE host labelled by IGAL.

  for ( i=0; i <= NMAGPSF_HOSTLIB ; i++ ) {

    if ( GALFRAC_SUM[i] > 1.000 )  { GALFRAC_SUM[i] = 1.000 ; }

    GALFRAC = GALFRAC_SUM[i] ;
    // load global arrays
    SNHOSTGAL.GALFRAC[i]  = GALFRAC ;

    if ( i == 0 ) 
      {  dm = 0.0 ; }
    else if ( GALFRAC == 0.0 ) {
      dm = 30.0 ;
    }
    else { 
      dm = -2.5*log10( GALFRAC );  // aperture mag  - total galmag
    }

    /*
    printf(" GGG xxx GALID=%d GALFRAC[%d]=%f   R=%6.3f asec  dm=%.3f\n", 
	   SNHOSTGAL.GALID, i, GALFRAC, HOSTLIB.Aperture_Radius[i], dm );
    */

    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      IVAR      = HOSTLIB.IVAR_MAGOBS[ifilt_obs] ;
      sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );

      if ( IVAR >= 0 ) {
	MAGOBS_LIB = HOSTLIB.VALUE_ZSORTED[IVAR][IGAL] ; 
	MAGOBS     = MAGOBS_LIB + MWXT[ifilt_obs] ; // dim  with Gal extinction
      }
      else {
	MAGOBS = 99.0 ; // allow missing filter in HOSTLIB 
      }

      if ( isnan(MAGOBS) ) {
	sprintf(c1err,"MAGOBS(%s)=NaN for IGAL=%d", cfilt, IGAL);
	sprintf(c2err,"Check GALID = %lld", SNHOSTGAL.GALID);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      SNHOSTGAL.GALMAG[ifilt_obs][i] = MAGOBS + dm ;
    }
  } // end i lopop over PSF bins


  // ----------------------
  // Mar 3 2015:
  // compute local surface brightness mag with effective PSF = 0.6''
  // Note that the effective area is 4*PI*PSF^2.

  double psfsig, arg, SB_MAG, SB_FLUXCAL, AREA ;

  psfsig  = INPUTS.HOSTLIB_SBRADIUS/2.0; 
  AREA    = 3.14159*(4.0*psfsig*psfsig) ; // effective noise-equiv area

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs  = GENLC.IFILTMAP_OBS[ifilt];
    SB_MAG     = interp_GALMAG_HOSTLIB(ifilt_obs, psfsig ) ;
    SB_MAG    += 2.5*log10(AREA); // normalize to 1 sq-arcsec
    
    // convert mag to fluxcal
    arg     = -0.4*(SB_MAG - ZEROPOINT_FLUXCAL_DEFAULT) ;
    SB_FLUXCAL = pow(10.0,arg) ;
    
    SNHOSTGAL.SB_FLUXCAL[ifilt_obs] = SB_FLUXCAL ;  
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

    print_preAbort_banner(fnam);
    printf("\t offset=%f, sig=%f) \n", offset, sig);
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

  int    j, NBIN ;
  double a, b, w, n, bn, rexp, sqsum,  arg ;
  double reduced_R, xx, yy, FSUM_PROFILE, F, FGAL_TOT    ;
  //  char fnam[] = "get_GALFLUX_HOSTLIB" ;

  // ---------------- BEGIN ---------------

  FSUM_PROFILE = 0.0 ;
  NBIN = SERSIC_TABLE.NBIN_reduced ;

  for ( j=0; j < SERSIC_PROFILE.NPROF; j++ ) {
    
    // strip off info for this galaxy component.
    a   = SNHOSTGAL.SERSIC.a[j] ;
    b   = SNHOSTGAL.SERSIC.b[j] ;
    n   = SNHOSTGAL.SERSIC.n[j] ;
    w   = SNHOSTGAL.SERSIC.w[j] ;
    bn  = SNHOSTGAL.SERSIC.bn[j] ;
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
void  STORE_SNHOST_MISC(int IGAL, int ibin_SNVAR) {

  // Feb 11 2014.
  // store misc variables in SNHOSTGAL structure.
  // Needed to write to data file.

  // Mar 27 2020: init I2SNMAGSHIFT=0 to fix crazy unitialized values.

  int  NDIM           = HOSTLIB_WGTMAP.GRIDMAP.NDIM;
  bool USE_SNVAR      = ( ibin_SNVAR >= 0 ) ;
  bool USE_SNMAGSHIFT = ((INPUTS.HOSTLIB_MSKOPT&HOSTLIB_MSKOPT_SNMAGSHIFT)>0);
  double  SNMAGSHIFT = 0.0, VAL, WGT, WGTSUM, WGTSUM_LAST ;
  int ivar, IVAR_STORE ;
  short int I2SNMAGSHIFT;
  char fnam[] = "STORE_SNHOST_MISC" ;

  // ------------ BEGIN ---------

  WGTSUM = WGTSUM_LAST = 0.0 ;
  I2SNMAGSHIFT = 0 ;
  
  // store wgt and SN mag shift for this host
  if ( USE_SNVAR ) {
    WGTSUM  = HOSTLIB_WGTMAP.WGTSUM_SNVAR[ibin_SNVAR][IGAL];
    if(IGAL>0) {WGTSUM_LAST=HOSTLIB_WGTMAP.WGTSUM_SNVAR[ibin_SNVAR][IGAL-1];}
    if ( USE_SNMAGSHIFT ) 
      { I2SNMAGSHIFT = HOSTLIB_WGTMAP.I2SNMAGSHIFT_SNVAR[ibin_SNVAR][IGAL];}
  }
  else {
    WGTSUM  = HOSTLIB_WGTMAP.WGTSUM[IGAL] ; 
    if(IGAL>0) { WGTSUM_LAST = HOSTLIB_WGTMAP.WGTSUM[IGAL-1] ; }
    if ( USE_SNMAGSHIFT ) 
      { I2SNMAGSHIFT = HOSTLIB_WGTMAP.I2SNMAGSHIFT[IGAL]; }

  }

  WGT        = WGTSUM - WGTSUM_LAST ;
  SNMAGSHIFT = (double)I2SNMAGSHIFT / I2MAGSCALE_HOSTLIB ;

  /*
  printf(" xxx %s: CID=%2d   IGAL=%d, ibin=%d, WGT=%.4f \n",
	 fnam, GENLC.CID, IGAL, ibin_SNVAR, WGT); fflush(stdout);
  */

  SNHOSTGAL.WGTMAP_WGT        = WGT;
  SNHOSTGAL.WGTMAP_SNMAGSHIFT = SNMAGSHIFT ;

  for ( ivar=0; ivar < NDIM; ivar++ ) {
    IVAR_STORE = HOSTLIB.IVAR_STORE[ivar];
    if ( HOSTLIB_WGTMAP.IS_SNVAR[ivar] ) 
      { VAL = -999.0 ; }    
    else 
      { VAL = HOSTLIB.VALUE_ZSORTED[IVAR_STORE][IGAL]; }
    SNHOSTGAL.WGTMAP_VALUES[ivar] = VAL ;
  }


  // -----------------------------------------------------------
  // Sep 16 2015
  // set GENLC.FIELDNAME if FIELD column is present in hostlib

  /* xxxxx mark delete Feb 21 2022 xxxxxxxxx
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
  xxxxxxxxxxx end mark xxxxxxxxxxxx */

  return ;

} // end of STORE_SNHOST_MISC


// =============================================
void LOAD_OUTVAR_HOSTLIB(int IGAL) {

  // Feb 2014.
  // load HOSTLIB_OUTVAR_EXTRA.VALUE array for output data file.
  // These variables are user-defined by sim-input key
  // HOSTLIB_STOREVAR:  <var1>,<var2>,<var3>, ...

  int NNBR = SNHOSTGAL.NNBR ; // includes true host
  int IGAL_NBR, i_NBR;
  int NVAR_OUT, ivar, IVAR_STORE ;
  double DVAL ;
  char fnam[] = "LOAD_OUTVAR_HOSTLIB" ;

  // ------------- BEGIN ------------
 
  NVAR_OUT = HOSTLIB_OUTVAR_EXTRA.NOUT ;
  if ( NVAR_OUT == 0 ) { return ; }

  for(ivar=0; ivar < NVAR_OUT; ivar++ ) {
    IVAR_STORE = HOSTLIB_OUTVAR_EXTRA.IVAR_STORE[ivar] ;
    DVAL       = HOSTLIB.VALUE_ZSORTED[IVAR_STORE][IGAL] ;
    HOSTLIB_OUTVAR_EXTRA.VALUE[ivar][0] = DVAL ;   

    i_NBR    = 1; 
    HOSTLIB_OUTVAR_EXTRA.VALUE[ivar][i_NBR] = NULLDOUBLE ;
    if ( NNBR > 1 ) {
      IGAL_NBR = SNHOSTGAL.IGAL_NBR_LIST[i_NBR];
      DVAL     = HOSTLIB.VALUE_ZSORTED[IVAR_STORE][IGAL_NBR] ;
      HOSTLIB_OUTVAR_EXTRA.VALUE[ivar][i_NBR] = DVAL ;
      
      //	   printf("xxx %s: IGAL_NBR = %d, %s = %f\n", 
      //	  fnam, IGAL_NBR, HOSTLIB_OUTVAR_EXTRA.NAME[ivar], DVAL);
    }

    //    printf(" xxx IGAL=%d  ivar=%d IVAR_STORE=%d DVAL=%f \n", 
    //	   IGAL, ivar, IVAR_STORE, DVAL);
  } // end IVAR loop

  return ;

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

  int  GALID          = SNHOSTGAL.GALID ;
  int  IGAL           = SNHOSTGAL.IGAL ;
  bool USE_GAMMA_GRID = HOSTLIB_WGTMAP.USE_SALT2GAMMA_GRID;
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
void DUMPROW_SNHOST(void) {

  // dump 1 row for easy readback.
  // <job.exe args> | grep 222222 > host.dump

  printf(" %d  %6.4f %6.4f %4.2f  %7.4f %7.4f %5.3f %6.2f  %f %f \n"    
	 , 222222    // for easy grep
	 , SNHOSTGAL.SERSIC.a[0]
	 , SNHOSTGAL.SERSIC.b[0]
	 , SNHOSTGAL.SERSIC.INDEX
	 , SNHOSTGAL.RA_SNGALSEP_ASEC 
	 , SNHOSTGAL.DEC_SNGALSEP_ASEC
	 , SNHOSTGAL.reduced_R
	 , SNHOSTGAL.phi/.0174533
	 , SNHOSTGAL.GALFRAC[2] 
	 , SNHOSTGAL.GALFRAC[4] 
	 );
  fflush(stdout);

} // end of DUMPROW_SNHOST

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

  if ( SERSIC_PROFILE.NPROF > 0 ) {

    printf("\t => GAL(RA,DEC) = %.6f, %.6f deg. \n",
	   SNHOSTGAL.RA_GAL_DEG,  SNHOSTGAL.DEC_GAL_DEG  );
    printf("\t => SN(RA,DEC)  = %.6f, %.6f deg. \n",
	   SNHOSTGAL.RA_SN_DEG,  SNHOSTGAL.DEC_SN_DEG  );
    printf("\t => SN-Gal sep(RA,DEC) : %6.3f , %6.3f arcsec \n"
	   , SNHOSTGAL.RA_SNGALSEP_ASEC
	   , SNHOSTGAL.DEC_SNGALSEP_ASEC );
    printf("\t => SNSEP=%.3f arcsec, DLR=%.3f  arcsec, d_DLR=%.3f \n", 
	   SNHOSTGAL.SNSEP, SNHOSTGAL.DLR, SNHOSTGAL.DDLR );

    if ( INDEX_GENMODEL == MODEL_SIMLIB ) {
      printf("\t => SIMLIB HEADER PARAMS: \n");
      printf("\t   LIBID=%d  GALID=%lld, RA=%.4f  DEC=%.4f  z=%.4f\n",
	     SIMLIB_HEADER.LIBID, SIMLIB_HEADER.GALID,
	     SIMLIB_HEADER.RA, SIMLIB_HEADER.DEC,
	     SIMLIB_HEADER.GENRANGE_REDSHIFT[0] );
    }

    fflush(stdout);

    for ( j=0; j < SERSIC_PROFILE.NPROF; j++ ) {
      a  =  SNHOSTGAL.SERSIC.a[j] ;
      b  =  SNHOSTGAL.SERSIC.b[j] ;
      w  =  SNHOSTGAL.SERSIC.w[j] ;
      n  =  SNHOSTGAL.SERSIC.n[j] ;
      bn =  SNHOSTGAL.SERSIC.bn[j] ;
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

  return ;

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
  // Initial use was for the sim to pass host par to BYOSED.
  //
  // Nov 5 2020: add REDSHIFT to list.

  int  NPAR=0, ivar ;
  int  NVAR_WGTMAP = HOSTLIB_WGTMAP.GRIDMAP.NDIM ;
  int  NVAR_EXTRA  = HOSTLIB_OUTVAR_EXTRA.NOUT ;
  char fnam[] = "fetch_HOSTPAR_GENMODEL";

  // ----------- BEGIN -----------

  if ( OPT == 1 ) {
    // always start with RV and AV    
    // xxx more delete    sprintf(NAMES_HOSTPAR,"RV,AV"); NPAR=2;
    sprintf(NAMES_HOSTPAR,"RV,AV,REDSHIFT");  NPAR=3;

    for ( ivar=0; ivar < NVAR_WGTMAP; ivar++ ) {  
      catVarList_with_comma(NAMES_HOSTPAR,HOSTLIB_WGTMAP.VARNAME[ivar]);
      NPAR++ ;
    }

    // tack on user-defined variables from sim-input HOSTLIB_STOREPAR key
    for(ivar=0; ivar < NVAR_EXTRA; ivar++ ) {
      if ( HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[ivar] ) { continue; }
      catVarList_with_comma(NAMES_HOSTPAR,HOSTLIB_OUTVAR_EXTRA.NAME[ivar]);
      NPAR++ ;
    }
    // add anything used in WGTMAP or HOSTLIB_STOREPAR
  }
  else if ( OPT ==  2 ) {
    VAL_HOSTPAR[0] = GENLC.RV ;
    VAL_HOSTPAR[1] = GENLC.AV ;
    VAL_HOSTPAR[2] = GENLC.REDSHIFT_CMB ;
    NPAR = 3 ; 

    for ( ivar=0; ivar < NVAR_WGTMAP; ivar++ ) {  
      VAL_HOSTPAR[NPAR] = SNHOSTGAL.WGTMAP_VALUES[ivar] ;
      NPAR++ ;
    }

    for(ivar=0; ivar < NVAR_EXTRA; ivar++ ) {
      if ( HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[ivar] ) { continue; }
      VAL_HOSTPAR[NPAR] = HOSTLIB_OUTVAR_EXTRA.VALUE[ivar][0];
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
  // Passed here via *HOSTLIB_APPEND.
  // 
  // July 16 2021: write DOCANA keys to FP_NEW

  char *SUFFIX       = HOSTLIB_APPEND->FILENAME_SUFFIX; // or new HOSTLIB
  int  NLINE_COMMENT = HOSTLIB_APPEND->NLINE_COMMENT;
  char *VARNAMES     = HOSTLIB_APPEND->VARNAMES_APPEND ;
  char *LINE         = (char*) malloc ( sizeof(char) * MXCHAR_LINE_HOSTLIB );

  int  gzipFlag;
  FILE *FP_ORIG, *FP_NEW;
  char *HLIB_ORIG = INPUTS.HOSTLIB_FILE;
  char  HLIB_TMP[MXPATHLEN], HLIB_NEW[100], DUMPATH[MXPATHLEN];
  char fnam[] = "rewrite_HOSTLIB" ;

  // -------------- BEGIN --------------

  // create local string with name of HOSTLIB
  sprintf(HLIB_TMP,"%s%s", HLIB_ORIG, SUFFIX ); 

  // remove path from HLIB_NEW to ensure that new file is created
  // locally and not in somebody else's directory.
  extract_MODELNAME(HLIB_TMP, DUMPATH, HLIB_NEW);
  
  // open orig HOSTLIB without checking for DOCANA so that
  // we don't skip DOCUMENTATION key
  FP_ORIG = open_TEXTgz(HLIB_ORIG, "rt", &gzipFlag );
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

  printf("\n  Created '%s' \n", HLIB_NEW);
  fflush(stdout);

  // - - - - - - - -
  // transfer DOCANA block (July 2021)
  bool DOCANA_END = false ;
  int  iline, NLINE_DOCANA = 0, ipad, NPAD=0 ;
  char *str_tmp, PAD_YAML[12]="";
  char KEY_DOCANA_PADCHECK[] = "PURPOSE:"; // key to measure pad space

  if ( INPUTS.REQUIRE_DOCANA ) {

    // keep re-writing HOSTLIB lines until reaching DOCUMENTATION_END key
    while ( !DOCANA_END ) {
      fgets(LINE, MXCHAR_LINE_HOSTLIB, FP_ORIG);

      // get yaml pad spacing
      str_tmp = strstr(LINE,KEY_DOCANA_PADCHECK);
      if ( str_tmp != NULL ) {
	NPAD = (str_tmp - LINE);
	for(ipad=0; ipad<NPAD; ipad++ ) { strcat(PAD_YAML," "); }
	// printf("\n xxx %s: NPAD=%d  PAD='%s' \n", fnam, NPAD, PAD_YAML);
      }
      // add new comments just before DOCANA end key
      if ( strstr(LINE,KEYNAME2_DOCANA_REQUIRED) != NULL ) { 
	DOCANA_END = true; 

	if ( NPAD == 0 ) {
	  char line_warn[] = "@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-" ;
	  printf("\n%s\n", line_warn);
	  printf(" WARNING: unable to measure yaml pad spacing\n");
	  printf("\t because %s key is missing.\n", KEY_DOCANA_PADCHECK);
	  printf("\t Add HOSTLIB_APPEND_NOTES to DOCANA yaml with "
		 "3 pad spaces as a guess.\n");
	  printf("%s\n\n", line_warn);
	  fflush(stdout);
	  sprintf(PAD_YAML,"   ");
	}
	fprintf(FP_NEW,"\n%sHOSTLIB_APPEND_NOTES:\n", PAD_YAML);
	for(iline=0; iline < NLINE_COMMENT; iline++ ) { 
	  fprintf(FP_NEW,"%s- %s\n",PAD_YAML,HOSTLIB_APPEND->COMMENT[iline]);
	}
      }

      fprintf(FP_NEW,"%s", LINE); 
      if ( DOCANA_END ) { fprintf(FP_NEW,"\n");  }

      NLINE_DOCANA++ ;

      if ( NLINE_DOCANA > MXLINE_DOCANA ) 
	{ abort_docana_tooLong(HLIB_ORIG, fnam);  } 

    }

    printf("  Wrote %d DOCANA lines to new HOSTLIB\n", NLINE_DOCANA);
  } // end REQUIRE_DOCANA


  // - - - - - 

  // - - - - - - - - - - - - - - - - - 
  // read each original line

  long long GALID, GALID_orig ;
  int   igal_unsort, igal_zsort, ivar, NWD_LINE, NLINE_GAL=0;
  char *LINE_APPEND, *FIRSTWORD, *NEXTWORD, *ptrCR ;

  LINE_APPEND  = (char*) malloc ( sizeof(char) * MXCHAR_LINE_APPEND );
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

	NLINE_GAL++ ;

	// make sure GALID matches
	ivar       = HOSTLIB.IVAR_GALID ;
	igal_zsort = HOSTLIB.LIBINDEX_ZSORT[igal_unsort] ;
	GALID      = (long long)HOSTLIB.VALUE_ZSORTED[ivar][igal_zsort] ;

	get_PARSE_WORD(0, 1, NEXTWORD);    // read GALID
	sscanf(NEXTWORD, "%lld", &GALID_orig);
	if ( GALID != GALID_orig ) {
	  sprintf(c1err,"GALID mis-match for igal_unsort=%d", igal_unsort);
	  sprintf(c2err,"GALID(orig)=%lld, but stored GALID=%lld",
		  GALID_orig, GALID ) ;
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	}

	sprintf(LINE_APPEND,"%s",HOSTLIB_APPEND->LINE_APPEND[igal_unsort]);

	/* xxx
	for(ifilt=1; ifilt <= NFILT; ifilt++ ) {
	  sprintf(cval, " %6.3f", MAG_STORE[ifilt][igal_unsort] );
	  strcat(LINE_APPEND,cval);
	}
	xxxxxxx */

	igal_unsort++ ;
      }
    }

    if ( NLINE_GAL >= INPUTS.HOSTLIB_MAXREAD ) { break; }

    ptrCR = strchr(LINE,'\n'); if(ptrCR){*ptrCR=' ';} // remove <CR>
    fprintf(FP_NEW,"%s %s\n", LINE, LINE_APPEND);
  }


  if(gzipFlag ) { pclose(FP_ORIG); }  else { fclose(FP_ORIG); }
  fclose(FP_NEW);
  printf("  Wrote %d GAL rows\n", NLINE_GAL);
  fflush(stdout);

  return ;

} // end rewrite_HOSTLIB

// =======================================
void malloc_HOSTLIB_APPEND(int NGAL, HOSTLIB_APPEND_DEF *HOSTLIB_APPEND) {

  // malloc and init
  int i;
  int MEMC = MXCHAR_LINE_APPEND * sizeof(char*) ;
  char fnam[] = "malloc_HOSTLIB_APPEND" ;

  // ---------------- BEGIN ---------------

  HOSTLIB_APPEND->NLINE_APPEND = NGAL;
  HOSTLIB_APPEND->LINE_APPEND = (char**) malloc( NGAL*sizeof(char*) );
  for(i=0; i < NGAL; i++ ) {
    HOSTLIB_APPEND->LINE_APPEND[i] = (char*) malloc(MEMC);
    sprintf(HOSTLIB_APPEND->LINE_APPEND[i],"NULL_APPEND");
  }

  HOSTLIB_APPEND->NLINE_COMMENT = 0;

  return;

} // end malloc_HOSTLIB_APPEND

// =========================================
void addComment_HOSTLIB_APPEND(char *COMMENT, 
			       HOSTLIB_APPEND_DEF *HOSTLIB_APPEND) {

  int NL   = HOSTLIB_APPEND->NLINE_COMMENT;
  int LENC = strlen(COMMENT) + 10;
  int MEMC = LENC * sizeof(char);
  
  HOSTLIB_APPEND->COMMENT[NL] = (char*) malloc(MEMC);
  sprintf(HOSTLIB_APPEND->COMMENT[NL],"%s", COMMENT);
  HOSTLIB_APPEND->NLINE_COMMENT++ ;

  return;

}  // end add_HOSTLIB_APPEND_COMMENT

// ===================================
void rewrite_HOSTLIB_plusMags(void) {
  
  // If +HOSTMAGS option is given on command-line, this function is
  // called to compute synthetic mag for each band in GENFILTERS key,
  // and re-write the HOSTLIB with extra [band]_obs columns. 
  // Beware that 'igal' is a redshift-sorted index for the other 
  // HOSTLIB functions, but here we use the original HOSTLIB order.
  // 

  int NGAL     = HOSTLIB.NGAL_STORE ;
  int NFILT    = NFILT_SEDMODEL ; 
  int NBIN_LAM = INPUTS_SPECTRO.NBIN_LAM ;
  int MEMD     = sizeof(double) * NBIN_LAM ;

  int igal_unsort, igal_zsort, ifilt, ifilt_obs, ivar, DUMPFLAG=0, LENLINE ;
  long long GALID ;
  double ZTRUE, MWEBV=0.0, mag, *GENFLUX_LIST, *GENMAG_LIST;
  float  *MAG_STORE ;

  HOSTLIB_APPEND_DEF HOSTLIB_APPEND ;
  char cval[20], LINE_APPEND[MXCHAR_LINE_APPEND] ;
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

      LENLINE = strlen(LINE_APPEND) ;
      if ( LENLINE > MXCHAR_LINE_APPEND ) {
	sprintf(c1err,"strlen(LINE_APPEND)=%d exceeds bound (NFILT=%d)", 
		LENLINE, NFILT);
	sprintf(c2err,"MXCHAR_LINE_APPEND = %d", MXCHAR_LINE_APPEND);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    } // end ifilt

    sprintf(HOSTLIB_APPEND.LINE_APPEND[igal_unsort],"%s", LINE_APPEND);

    if ( (igal_unsort % 10000) == 0 ) {
      printf("\t Processing igal %8d of %8d \n", igal_unsort, NGAL);
      fflush(stdout);
    }

    // xxxxxxxxxxxxxxxxxxxxxxxxxxx
    if ( ZTRUE < -1.0 ) {
      printf(" xxx ------------------------------------- \n");
      printf(" xxx igal_unsort=%2d  GALID=%8lld  ZTRUE=%.5f \n", 
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
  char VARNAMES_HOSTMAGS[5*MXFILTINDX], varname_mag[40], *cfilt, msg[80];
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

  double TRANS, LAMOBS, LAMREST, FLAM, FSUM=0.0, mag=0.0 ;
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

  int   igal_unsort, igal_zsort, igal_DECsort  ;

  HOSTLIB_APPEND_DEF HOSTLIB_APPEND;
  char  *LINE_APPEND, MSG[200] ;

  // internal debug
  int  NGAL_DEBUG  = INPUTS.HOSTLIB_MAXREAD ;

  char *INPUT_FILE   = INPUTS.INPUT_FILE_LIST[0];
  char *HOSTLIB_FILE = INPUTS.HOSTLIB_FILE;
  char SUFFIX_HOSTNBR[] = "+HOSTNBR" ;
  char fnam[] = "rewrite_HOSTLIB_plusNbr" ;

  // --------------- BEGIN ---------------

  print_banner(fnam);

  printf("Append up to %d host neighbors within %.1f'' radius.",
	  HOSTLIB_NBR_WRITE.NNBR_WRITE_MAX, HOSTLIB_NBR_WRITE.SEPNBR_MAX );
  fflush(stdout);

  if ( strstr(HOSTLIB_FILE,SUFFIX_HOSTNBR) != NULL ) {
    sprintf(c1err,"HOSTLIB already has NBR_LIST");
    sprintf(c2err,"Check HOSTLIB_FILE='%s'", HOSTLIB_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( IVAR_RA < 0 || IVAR_DEC < 0 ) {
    sprintf(c1err,"Must include galaxy coords to find neighbors.");
    sprintf(c2err,"Check VARNAMES in HOSTLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  printf("\t NGAL[READ,STORE] = %d, %d \n", 
	 HOSTLIB.NGAL_READ, HOSTLIB.NGAL_STORE);

  malloc_HOSTLIB_APPEND(NGAL, &HOSTLIB_APPEND);

  LINE_APPEND = (char*) malloc (MXCHAR_LINE_HOSTLIB * sizeof(char) ) ;
  HOSTLIB_NBR_WRITE.SKY_SORTED_DEC          = (double*) malloc(MEMD) ;
  HOSTLIB_NBR_WRITE.SKY_SORTED_RA           = (double*) malloc(MEMD) ;
  HOSTLIB_NBR_WRITE.SKY_SORTED_IGAL_zsort   = (int*) malloc(MEMI) ;
  HOSTLIB_NBR_WRITE.SKY_SORTED_IGAL_DECsort = (int*) malloc(MEMI) ;
  HOSTLIB_NBR_WRITE.GALID_atNNBR_MAX = -9 ;
  HOSTLIB_NBR_WRITE.NNBR_MAX         =  0 ;


  // sort by DEC to improve NBR-matching speed
  int  ORDER_SORT = +1 ;
  double *ptrDEC = HOSTLIB.VALUE_ZSORTED[IVAR_DEC] ; 
  double *ptrRA  = HOSTLIB.VALUE_ZSORTED[IVAR_RA] ; 

  sortDouble( NGAL, ptrDEC, ORDER_SORT, 
	      HOSTLIB_NBR_WRITE.SKY_SORTED_IGAL_DECsort);
  
  // load new lists of RA & DEC sorted by DEC
  for(igal_DECsort=0; igal_DECsort < NGAL; igal_DECsort++ ) {
    igal_zsort = HOSTLIB_NBR_WRITE.SKY_SORTED_IGAL_DECsort[igal_DECsort];
    HOSTLIB_NBR_WRITE.SKY_SORTED_DEC[igal_DECsort]      = ptrDEC[igal_zsort] ;
    HOSTLIB_NBR_WRITE.SKY_SORTED_RA[igal_DECsort]       = ptrRA[igal_zsort] ;
    HOSTLIB_NBR_WRITE.SKY_SORTED_IGAL_zsort[igal_zsort] = igal_DECsort ;
    if ( igal_DECsort < -5 || igal_DECsort > NGAL+5 ) {
      printf(" xxx igal_DECsort=%d  igal_zsort=%6d,  DEC = %9.5f \n",
	     igal_DECsort, igal_zsort, ptrDEC[igal_zsort] ); fflush(stdout);
    }
  }
  
  // ----------------------------
  if ( NGAL_DEBUG < MXROW_HOSTLIB ) { NGAL = NGAL_DEBUG; }

  // init diagnistic counters (filled in get_LINE_APPEND_HOSTLIB_plusNbr)
  monitor_HOSTLIB_plusNbr(0,&HOSTLIB_APPEND); 

  // loop over all galaxies and prepare string to append.
  for(igal_unsort=0; igal_unsort < NGAL; igal_unsort++ ) {

    // search for neighbors and fill line to append
    get_LINE_APPEND_HOSTLIB_plusNbr(igal_unsort, LINE_APPEND);
    
    fflush(stdout);

    sprintf(HOSTLIB_APPEND.LINE_APPEND[igal_unsort],"%s", LINE_APPEND);

  } // end igal_unsort loop over all galaxies


  // - - - - - - - - - - - - 

  sprintf(HOSTLIB_APPEND.FILENAME_SUFFIX, "%s", SUFFIX_HOSTNBR );
  sprintf(HOSTLIB_APPEND.VARNAMES_APPEND, "NBR_LIST" );


  // construct message strings for top of new HOSTLIB
  sprintf(MSG, "Append up to %d host neighbors within %.1f'' radius.",
	  HOSTLIB_NBR_WRITE.NNBR_WRITE_MAX, HOSTLIB_NBR_WRITE.SEPNBR_MAX );
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);
 
  sprintf(MSG, "snlc_sim.exe %s +HOSTNBR  SEPNBR_MAX %.1f  NNBR_WRITE_MAX %d",
	  INPUT_FILE, HOSTLIB_NBR_WRITE.SEPNBR_MAX,
	  HOSTLIB_NBR_WRITE.NNBR_WRITE_MAX );
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);

  sprintf(MSG,"Added column NBR_LIST = "
	  "comma-sep list of row numbers (not GALID)" );
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);

  // write out monitor info
  monitor_HOSTLIB_plusNbr(1,&HOSTLIB_APPEND);

  // execute re-write
  rewrite_HOSTLIB(&HOSTLIB_APPEND);

  exit(0);

  return ;

} // end rewrite_HOSTLIB_plusNbr


// ==============================
void get_LINE_APPEND_HOSTLIB_plusNbr(int igal_unsort, char *LINE_APPEND) {


#define MXNNBR_STORE 100         // max number of neighbors to track
  double SEPNBR_MAX      = HOSTLIB_NBR_WRITE.SEPNBR_MAX ;
  int    NNBR_WRITE_MAX  = HOSTLIB_NBR_WRITE.NNBR_WRITE_MAX ;
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
  igal_DECsort = HOSTLIB_NBR_WRITE.SKY_SORTED_IGAL_zsort[igal_zsort];
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
      igal2_zsort   = HOSTLIB_NBR_WRITE.SKY_SORTED_IGAL_DECsort[isort];
      igal2_unsort  = HOSTLIB.LIBINDEX_UNSORT[igal2_zsort];

      RA_NBR      = HOSTLIB_NBR_WRITE.SKY_SORTED_RA[isort] ;  
      DEC_NBR     = HOSTLIB_NBR_WRITE.SKY_SORTED_DEC[isort] ;  
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
      if ( NNBR > HOSTLIB_NBR_WRITE.NNBR_MAX ) { 
	HOSTLIB_NBR_WRITE.NNBR_MAX = NNBR; 
	HOSTLIB_NBR_WRITE.GALID_atNNBR_MAX = GALID;
      }

    } // end j loop over smaller/larger DEC

    ISORT_CHANGE++ ;
  } // end while


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
  LSTDOUT = ( igal_unsort < 10 || 
	      igal_unsort > NGAL-10 ||
	      GALID == 293050
	      );
  if ( LSTDOUT  && strlen(LINE_STDOUT) > 0 ) {
    printf("# ---------------------------------------------------- \n");
    printf(" Crosscheck dump for igal_read=%d, GALID=%lld \n",
	   igal_unsort, GALID); fflush(stdout);
    
    printf("\t READ  NBR list: %s\n", LINE_APPEND );
    printf("\t GALID NBR list: %s\n", LINE_STDOUT );
    fflush(stdout);
  }

  if ( (igal_unsort % 10000) == 0 ) {
    NNBR  = HOSTLIB_NBR_WRITE.NNBR_MAX ; 
    GALID = HOSTLIB_NBR_WRITE.GALID_atNNBR_MAX;
    printf("\t Processing igal %8d of %8d  (NNBR_MAX=%2d for GALID=%lld)\n", 
	   igal_unsort, NGAL, NNBR, GALID );
    fflush(stdout);
  }

  if ( NNBR < 100 ) { HOSTLIB_NBR_WRITE.NGAL_PER_NNBR[NNBR]++ ; }
  if ( TRUNCATE   ) { HOSTLIB_NBR_WRITE.NGAL_TRUNCATE++ ; }

  return ;

} // end get_LINE_APPEND_HOSTLIB_plusNbr

// =========================================
void  monitor_HOSTLIB_plusNbr(int OPT, HOSTLIB_APPEND_DEF *HOSTLIB_APPEND) {

  // OPT=0 --> init
  // OPT=1 --> write info to stdout and to addComment

  int NGAL        = HOSTLIB.NGAL_STORE;
  int nnbr, NGAL_TMP;
  float frac;
  char MSG[200];
  char fnam[] = "monitor_HOSTLIB_plusNbr";
  // ------------ BEGIN -------------

  if ( OPT == 0 ) {
    for(nnbr=0; nnbr < 100; nnbr++ ) 
      { HOSTLIB_NBR_WRITE.NGAL_PER_NNBR[nnbr]=0; }
    HOSTLIB_NBR_WRITE.NGAL_TRUNCATE = 0 ;
  }
  else {
    printf("\n");
    for(nnbr=0; nnbr <= HOSTLIB_NBR_WRITE.NNBR_WRITE_MAX; nnbr++ ) {
      NGAL_TMP = HOSTLIB_NBR_WRITE.NGAL_PER_NNBR[nnbr];
      frac     = (float)NGAL_TMP / (float)NGAL;
      sprintf(MSG, "\t HOSTLIB fraction with %2d NBR: %8.3f %% ",
	     nnbr, 100.0*frac ); 
      printf("%s\n", MSG); fflush(stdout);
      addComment_HOSTLIB_APPEND(MSG,HOSTLIB_APPEND);
    }

    frac = (float)HOSTLIB_NBR_WRITE.NGAL_TRUNCATE / (float)NGAL ;
    sprintf(MSG,"\t Truncated fraction with > %d NBR: %8.3f %% ",
	   HOSTLIB_NBR_WRITE.NNBR_WRITE_MAX, 100.*frac); 
    printf("%s\n", MSG); fflush(stdout);
    addComment_HOSTLIB_APPEND(MSG,HOSTLIB_APPEND);
  }
  return ;

} // end monitor_HOSTLIB_plusNbr

// **************************************************
void rewrite_HOSTLIB_plusAppend(char *append_file) {

  // May 4 2021
  // Read columns in *append_file, and append to original HOSTLIB.
  // 

  int NGAL        = HOSTLIB.NGAL_STORE;
  HOSTLIB_APPEND_DEF HOSTLIB_APPEND ;
  int  NROW, NVAR_APPEND, ivar ;
  int  optMask       = 1 ;
  int  NMISSING      = 0;
  char tableName[]   = "HOSTLIB" ;
  char KEY_ALLVAR[]  = "ALL" ;
  char *varName, varList[100], *LINE_APPEND ;
  char fnam[]      = "rewrite_HOSTLIB_plusAppend" ;

  // ------------- BEGIN -------------

  print_banner(fnam);

  // read append file
  NROW = SNTABLE_AUTOSTORE_INIT(append_file, tableName, KEY_ALLVAR, optMask);
  NVAR_APPEND = READTABLE_POINTERS.NVAR_TOT; 

  malloc_HOSTLIB_APPEND(NGAL, &HOSTLIB_APPEND);
  LINE_APPEND = (char*) malloc (MXCHAR_LINE_APPEND * sizeof(char) ) ;

  varList[0] = 0;
  for( ivar = 1 ; ivar < NVAR_APPEND; ivar++ ) {
    varName = READTABLE_POINTERS.VARNAME[ivar];
    strcat(varList,varName);
    strcat(varList," ");
  }


  sprintf(HOSTLIB_APPEND.VARNAMES_APPEND, "%s", varList);
  sprintf(HOSTLIB_APPEND.FILENAME_SUFFIX, "%s", "+APPEND");
  
  int  igal_unsort, igal_zsort, istat_read;
  long long GALID;
  char cGALID[40], cVAL[40];
  double dVAL ;


  for(igal_unsort=0; igal_unsort < NGAL; igal_unsort++ ) {
    igal_zsort = HOSTLIB.LIBINDEX_ZSORT[igal_unsort];
    ivar  = HOSTLIB.IVAR_GALID;
    GALID = (long long)HOSTLIB.VALUE_ZSORTED[ivar][igal_zsort] ;
    sprintf(cGALID,"%lld", GALID);
    
    LINE_APPEND[0] = 0;
    for( ivar = 1 ; ivar < NVAR_APPEND; ivar++ ) {
      varName = READTABLE_POINTERS.VARNAME[ivar];
      SNTABLE_AUTOSTORE_READ(cGALID, varName, &istat_read, &dVAL, cVAL); 
      if ( istat_read == 0 ) {
	sprintf(cVAL, "%.5f ", dVAL);
      }
      else {
	dVAL = -9.0 ;  if ( ivar==1 ) { NMISSING++ ; }
	sprintf(cVAL, "-9.0 ");
      }
      strcat(LINE_APPEND,cVAL);
    }

    sprintf(HOSTLIB_APPEND.LINE_APPEND[igal_unsort],"%s", LINE_APPEND);

  } // end igal


  char MSG[200]; 
  sprintf(MSG,"%s appended %d variables (%s)", 
	  getenv("USER"), NVAR_APPEND, varList);
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);

  sprintf(MSG,"Append file is %s", append_file);
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);

  sprintf(MSG,"Missing %d GALIDs in append file.", NMISSING);
  addComment_HOSTLIB_APPEND(MSG, &HOSTLIB_APPEND);

  // execute re-write
  rewrite_HOSTLIB(&HOSTLIB_APPEND);
 
  free(LINE_APPEND);

  if ( NMISSING > 0 ) {
    printf("\t WARNING: Missing %5d GALIDs (wrote -9 for %s)\n",
	   NMISSING, varList ); 
    fflush(stdout);
  }

  exit(0);

  return ;

} // end rewrite_HOSTLIB_plusAppend

