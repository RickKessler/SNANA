/*******************************************
  snlc_sim:  Mar 2006: Created by R.Kessler
             Apr 2007: Remove SQL dependencies for public usage. 
             Feb 2009: add nonIa option
             Oct 2010: add GRID option (for psnid)
             Mar 2011: install snhost package
             May 2011: install snfitsio (default sim-> fits format)
             Apr 2012: use wavelength-dependent intrinsic smear in 
                       genSmear_models.c
             Jan 2014: separate trigger code into sntools_trigger.c[h]
             Jan 2017: add SPECTROGRAPH 
             Aug 2017: refactor SIMLIB_read

 ---------------------------------------------------------

  Fast Lightcurve simulator (LCSIM). 
  There are no images and no photometry !
  Survey conditions (PSF,SKY,ZPT) are used to calculate
  photometry and errors analytically.

  Usage:
    snlc_sim.exe  <input_file>
                       ^
                       |
    sample input files: $SNDATA_ROOT/analysis/sample_input_files/

********************************************/

#include "fitsio.h"
#include "MWgaldust.h"
#include "sntools.h"
#include "sntools_cosmology.h"
#include "snlc_sim.h"
#include "sntools_devel.h"
#include "sntools_host.h"
#include "sntools_weaklens.h"
#include "sntools_stronglens.h"
#include "sntools_wronghost.h"
#include "sntools_nonlinearity.h"
#include "sntools_fluxErrModels.h"
#include "sntools_fluxErrModels_legacy.h"
#include "sntools_genSmear.h"
#include "sntools_dataformat_fits.h"
#include "sntools_dataformat_text.h"
#include "sntools_kcor.h"
#include "sntools_trigger.h" 
#include "sntools_modelgrid.h"
#include "sntools_spectrograph.h"
#include "sntools_genPDF.h"
#include "sntools_output.h"   // added Jan 11 2017
#include "inoue_igm.h"        // added Jun 27 2019

#include "genmag_ALL.h"
#include "MWgaldust.h" 

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <sys/stat.h>

// include C code
#include "SNcadenceFoM.c"

#define MODELGRID_GEN
#define TWO_RANDONE_STREAMS

// ******************************************
int main(int argc, char **argv) {

  int ilc, istat, i  ;

  // define local structures
  SIMFILE_AUX_DEF  SIMFILE_AUX ;
  FILE *FP;
  char fnam[] = "main"; 

  // ------------- BEGIN --------------

  sprintf(BANNER,"Begin execution of snlc_sim.exe  " );
  print_banner(BANNER);

  //  errmsg(SEV_FATAL, 0, "main", "testing CodeTest", "Remove this" ); 

  TIMERS.t_start = time(NULL);

  init_SNPATH();    // init PATH_SNDATA_ROOT

  init_commandLine_simargs(argc, argv);

  // init fortran variables
  istat = 0;
  init_snvar__(&istat); 

  // one-time init of sim-variables
  init_simvar();

  //  test_igm(); // xxxx

  // read user input file for directions
  get_user_input();

  // init random number generator, and store first random.
  if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID  ) 
    { init_random_seed(INPUTS.ISEED, INPUTS.NSTREAM_RAN); }

  // prepare user input after init_random_seed to allow 
  // random systematic shifts.
  prep_user_input();

  // check for random CID option (after randoms are inited above)
  init_CIDRAN();

  // init based on GENSOURCE
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENRANDOM ) {
    init_RANDOMsource();

    // prepare randome systematic shifts after reading SURVEY from SIMLIB,
    // but before init_RateModel
    prep_RANSYSTPAR() ;

    init_RateModel();

    if ( INPUTS.INIT_ONLY == 1 ) { debugexit("main: QUIT AFTER RATE-INIT"); }
  }
  else if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
#ifdef MODELGRID_GEN
    init_GRIDsource(0); 
#endif
  }  
  else {
    sprintf(c1err,"%s is invalid GENSOURCE", INPUTS.GENSOURCE );
    errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
  }


  // abort on NGEN=0 after printing N per season (init_Rate above)
  if ( INPUTS.NGEN_LC <= 0 && INPUTS.NGENTOT_LC <= 0 ) {
    sprintf(c1err,"NGEN_LC=0 & NGENTOT_LC=0" );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }

  // init version for simulated output
  init_VERSION(INPUTS.GENVERSION);

  // modify SIM  path with genversion subdir; 
  prep_simpath();

  // skip to optional SIMLIB_DUMP after random number init.
  if ( INPUTS.SIMLIB_DUMP >= 0 ) { goto SIMLIB_DUMP ; }

  // init host-galaxy model (in sntools_host.c)
  INIT_HOSTLIB();  

  // init weak and strong lensing 
  if ( INDEX_GENMODEL != MODEL_LCLIB )  { 
    init_lensDMU(INPUTS.WEAKLENS_PROBMAP_FILE, INPUTS.WEAKLENS_DSIGMADZ ); 
    init_stronglens(INPUTS.STRONGLENS_FILE);
  }


  // init model fudges for flux-errors (Feb 2018)
  INIT_FLUXERRMODEL(INPUTS.FLUXERRMODEL_OPTMASK,  
		    INPUTS.FLUXERRMODEL_FILE,
		    INPUTS.FLUXERRMODEL_REDCOV,
		    INPUTS.FLUXERRMAP_IGNORE_DATAERR);


  // init anomalous host-subtraction noise
  INIT_NOISEMODEL_HOST_LEGACY(INPUTS.HOSTNOISE_FILE);

  // init WRONGHOST model
  double zMIN = INPUTS.GENRANGE_REDSHIFT[0] - 0.01 ;
  double zMAX = INPUTS.GENRANGE_REDSHIFT[1] + 0.01 ;
  INIT_WRONGHOST(INPUTS.WRONGHOST_FILE, zMIN, zMAX);


  // init user-specified z-dependence of sim parameters
  init_zvariation();

  // test triger init here ...
  if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID )  { 
    // init search efficiency
    init_SEARCHEFF(GENLC.SURVEY_NAME,INPUTS.APPLY_SEARCHEFF_OPT); 
  } 

  DASHBOARD_DRIVER();

  // initialize model that generates magnitudes.
  if ( INPUTS.USE_KCOR_LEGACY   ) { init_kcor_legacy(INPUTS.KCOR_FILE); }
  if ( INPUTS.USE_KCOR_REFACTOR ) { init_kcor_refactor(); }


  init_genPDF(INPUTS.GENPDF_OPTMASK, NULL,
	      INPUTS.GENPDF_FILE, INPUTS.GENPDF_IGNORE ) ;

  prioritize_genPDF_ASYMGAUSS();

  init_genmodel();
  init_modelSmear(); 
  init_genSpec();     // July 2016: prepare optional spectra

  // check options to rewrite hostlib and quit
  rewrite_HOSTLIB_DRIVER();

  // create/init output sim-files
  init_simFiles(&SIMFILE_AUX);

  // check option to dump rest-frame mags
 SIMLIB_DUMP:
  if ( INPUTS.USEFLAG_DMPTREST ) { DUMP_GENMAG_DRIVER() ; }

  // check for simlib-dump
  SIMLIB_DUMP_DRIVER();


  if ( INPUTS.INIT_ONLY ==2 ) { debugexit("main: QUIT AFTER FULL INIT"); }

  // =================================================
  // start main loop over "ilc"

  print_banner( " Begin Generating Lightcurves. " );
  fflush(stdout);

  TIMERS.t_end_init    = time(NULL); // Mar 15 2020
  TIMERS.t_update_last = TIMERS.t_end_init;
  TIMERS.NGENTOT_LAST  = 0 ;

  // - - - - - 
  for ( ilc = 1; ilc <= INPUTS.NGEN ; ilc++ ) {

    NGENLC_TOT++;

    if ( INPUTS.TRACE_MAIN  ) { dmp_trace_main("01", ilc) ; }

    if ( fudge_SNR() == 2 ) { goto GETMAGS ; }

    init_GENLC();

    if ( INPUTS.TRACE_MAIN  ) { dmp_trace_main("02", ilc) ; }

    if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID ) 
      { fill_RANLISTs(); }      // init list of random numbers for each SN    
    gen_event_driver(ilc); 

    if ( GENLC.STOPGEN_FLAG ) { NGENLC_TOT--;  goto ENDLOOP ; }

    if ( GENLC.NEPOCH < INPUTS.CUTWIN_NEPOCH[0] ) {   // avoid NEPOCH=0
      gen_event_reject(&ilc, &SIMFILE_AUX, "NEPOCH");
      goto GENEFF;
    }

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("03", ilc) ; }

    // build filter-maps & logical after shapepar in case
    // new filters are read (i.e, for non1a)
    gen_filtmap(ilc);   

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("04", ilc) ; }

    // Quit if we get negative LIBID;
    // mainly to stop reading fakes at end of file.

    if ( GENLC.SIMLIB_ID < 0 ) { continue ; }

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("05", ilc) ;  }

    // apply generation cuts
    if ( GENRANGE_CUT() == 0  ) {
      gen_event_reject(&ilc, &SIMFILE_AUX, "GENRANGE");
      goto GENEFF;
    }

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("06", ilc) ; }

  GETMAGS:

    // first check if peakMag-dependent trigger fails (to speed generation)
    if ( gen_TRIGGER_PEAKMAG_SPEC() == 0 ) { 
      gen_event_reject(&ilc, &SIMFILE_AUX, "SEARCHEFF");
      goto GENEFF; 
    }

    // now check zHOST-dependent efficiency (Dec 1 2017)
    if ( gen_TRIGGER_zHOST() == 0 ) { 
      gen_event_reject(&ilc, &SIMFILE_AUX, "SEARCHEFF");
      goto GENEFF; 
    }


    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("07", ilc) ; }
    GENMAG_DRIVER();   // July 2016

    if ( GENMAG_CUT() == 0  ) {
      gen_event_reject(&ilc, &SIMFILE_AUX, "GENMAG");
      goto GENEFF;
    }

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("08", ilc) ; }

    // generate spectra before broadband fluxes in case TEXPOSE
    // is computed from requested SNR; TEXPOSE is then used for
    // synthetic bands.
    GENSPEC_DRIVER(); 

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("09", ilc) ; }

    // convert generated mags into observed fluxes
    GENFLUX_DRIVER(); 

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("10", ilc) ; }

    // Jan 2016: Force ACCEPT if too many repeats with same LIBID 
    // to keep track of LIBIDs with no accepts.
    // In this case, set NOBS=NEPOCH=0 for data file since the
    // event has actually been discarded.
    if ( GENLC.NGEN_SIMLIB_ID >= SIMLIB_MXGEN_LIBID )  { 
      GENLC.ACCEPTFLAG_FORCE = 1; 
      GENLC.NOBS = GENLC.NEPOCH = 0; 
    }

    // check if search finds this SN:
    GENLC.SEARCHEFF_MASK = 3 ;
    if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID  ) {
      MJD_DETECT_DEF MJD_DETECT;
      LOAD_SEARCHEFF_DATA();
      GENLC.SEARCHEFF_MASK = 
	gen_SEARCHEFF(GENLC.CID                 // (I) ID for dump/abort
		      ,&GENLC.SEARCHEFF_SPEC     // (O)
		      ,&GENLC.SEARCHEFF_zHOST    // (O) Mar 2018
		      ,&MJD_DETECT   );          // (O) Oct 2021

      GENLC.MJD_TRIGGER        = (float)MJD_DETECT.TRIGGER ;
      GENLC.MJD_DETECT_FIRST   = (float)MJD_DETECT.FIRST ;
      GENLC.MJD_DETECT_LAST    = (float)MJD_DETECT.LAST ;
    }

    for ( i=1; i<= GENRAN_INFO.NLIST_RAN ; i++ )  
      { GENRAN_INFO.RANLAST[i] = getRan_Flat1(i); }

    // if APPLY opt is set, then require search MASK to keep SN;
    // otherwise keep all SNe    
    int OVP ;
    OVP = ( INPUTS.APPLY_SEARCHEFF_OPT & GENLC.SEARCHEFF_MASK ) ;
    if ( OVP != INPUTS.APPLY_SEARCHEFF_OPT  && GENLC.ACCEPTFLAG_FORCE==0 ) {
      gen_event_reject(&ilc, &SIMFILE_AUX, "SEARCHEFF");
      goto GENEFF ;
    }

    // check for Spectroscopic typing tag (not related to GENSPEC_DRIVER)
    gen_spectype();

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("11", ilc) ; }

    // check option to apply CUT windows
    if ( gen_cutwin() != SUCCESS  && GENLC.ACCEPTFLAG_FORCE==0 ) {
      gen_event_reject(&ilc, &SIMFILE_AUX, "CUTWIN");
      goto GENEFF;
    }

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("12", ilc) ; }

    // check for dump option
    dmp_event(ilc);

    // check option to lock first accepted LIBID
    if ( INPUTS.SIMLIB_IDLOCK == 1 ) 
      { GENLC.SIMLIB_IDLOCK = GENLC.SIMLIB_ID ; }

    // update various counters for accepted event.
    update_accept_counters();

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("13", ilc) ; }

    // update SNDATA files & auxiliary files
    update_simFiles(&SIMFILE_AUX);

    GENLC.ACCEPTFLAG = 1 ;  // Added Dec 2015

    if ( INPUTS.TRACE_MAIN ) { dmp_trace_main("14", ilc) ; }

  GENEFF:

    if ( INPUTS.NGENTOT_LC > 0 ) { screen_update(); }

    GENLC.STOPGEN_FLAG = geneff_calc();  // calc generation effic & error  
    if ( GENLC.STOPGEN_FLAG )  { goto ENDLOOP; }
    
    fflush(stdout);
    
  } // end of ilc loop


 ENDLOOP:

  TIMERS.t_end = time(NULL);

  // print final statistics on generated lightcurves.

  simEnd(&SIMFILE_AUX);

  return(0);

} // end of main



// **************************
void init_commandLine_simargs(int argc, char **argv) {

  // Created April 2012.
  //  + print full command to stdout
  //  + Extract input_file name 
  //  + store user command-line args in ARGV_LIST array
  //  + check KEY_DUMP
  // 
  // Jul 21 2020
  //   + malloc ARGV_LIST
  //   + abort of NARGV_LIST >= MXARGV
  // Jun 02 2021: use strlen(argv[i]) to malloc ARGV

  int i, LENARG;
  char *inFile = INPUTS.INPUT_FILE_LIST[0] ;
  char fnam[] = "init_commandLine_simargs";

  // ---------------- BEGIN --------------

  printf("   Full command: ");

  if ( argc >= 2 ) {
    sprintf(inFile, "%s", argv[1] );
    NARGV_LIST = argc ;

    if ( NARGV_LIST >= MXARGV ) {
      sprintf(c1err,"%d command line args exceeds MXARGV=%d", 
	      NARGV_LIST, MXARGV);
      sprintf(c2err,"Either reduce number of args, or increase MXARGV");
      errmsg(SEV_WARN, 0, fnam, c1err, c2err); 
    }

    for ( i = 0; i < NARGV_LIST ; i++ ) {
      LENARG = strlen(argv[i]) + 10;
	// xxx mark ARGV_LIST[i]=(char*) malloc( MXPATHLEN*sizeof(char) );
      ARGV_LIST[i] = (char*) malloc( LENARG*sizeof(char) );
      sprintf( ARGV_LIST[i], "%s", argv[i] );
      USE_ARGV_LIST[i] = 0 ;
      printf("%s ", argv[i]);
    }
    USE_ARGV_LIST[0] = 1;  // program name
    USE_ARGV_LIST[1] = 1;  // input file
  }
  else
    { sprintf(inFile, "snlc_sim.input" ); } // default name

  printf("\n\n"); fflush(stdout);


  // 7.17.2020: prepare KEY_DUMP option to read a dummy file
  if ( strcmp(inFile,"KEY_DUMP") == 0 ) {
    INPUTS.KEYNAME_DUMPFLAG = true;  
    sprintf(inFile,"%s/SURVEY.DEF", PATH_SNDATA_ROOT); // dummy file to read
  }
  else {
    INPUTS.KEYNAME_DUMPFLAG   = false; 
  }

  return ;

} // end of init_commandLine_simargs


// ***************************
void simEnd(SIMFILE_AUX_DEF *SIMFILE_AUX) {

  char fnam[] = "simEnd" ;

  // after all SN are generated, wrap up.
  
  sprintf(BANNER, " Done generating %d SN lightcurves from %s source.", 
	  NGENLC_TOT, INPUTS.GENSOURCE );
  print_banner(BANNER);
  printf("\t (%d lightcurves requested => %d were written) \n",
	 INPUTS.NGEN, NGENLC_WRITE );

  if ( NGEN_ALLSKIP >= INPUTS.NGEN ) {
    sprintf(c1err,"No generated observations because of invalid passbands.");
    sprintf(c2err,"Check wavelength ranges for filters and %s",
	    INPUTS.GENMODEL);
    //    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    errmsg(SEV_WARN, 0, fnam, c1err, c2err);  // 9.24.2017
  }

  END_FLUXERRMODEL();

  end_simFiles(SIMFILE_AUX);

  if ( NAVWARP_OVERFLOW[0] > 0 ) 
    { printf("%s", WARNING_AVWARP_OVERFLOW ); }

  if ( NGEN_REJECT.NEPOCH == NGENLC_TOT ) {
    sprintf ( c1err, "Every generated event fails NEPOCH>=%d .",
	      (int)INPUTS.CUTWIN_NEPOCH[0] );
    sprintf ( c2err, "GENRANGE_PEAKMJD and SIMLIB probably don't overlap."); 
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n DONE with snlc_sim.\n");
  fflush(stdout);
  
} // end of simEnd


// *************************
int LUPDGEN(int N) {
  // May 27, 2009
  // return 1 to make screen-update
  // Apr 26 2017: float -> double to handle very large N
  // Jun 25 2019: use % instead of fmod.
  // ------------- BEGIN ---------------
  if ( N == 1   ) { return(1); }
  if ( (N % INPUTS.NGEN_SCREEN_UPDATE)== 0 ) { return(1); }
  return 0;
}

// ******************************************
void get_user_input(void) {

  /**********

  Get user input from four different priority levels:
  (higher priority number overrides value from lower number)
  
  1 primary input file
  2 INCLUDE files in primary
  3 INCLUDE files on command line
  4 command line args

  Jun 18 2017: move prep_user_input into main, after ranSeed init
  Jul 30 2018: check input_file_include2
  Jul 20 2020: 
     + check a few command-line args before reading any files.
  Jun 25 2021: update logic to implment priority list above.

  ***********/
  int i, ifile ;
  bool FOUNDKEY[2];

  char fnam[] = "get_user_input"    ;

  // ------------ BEGIN ---------------

  // set hard-wired default values

  set_user_defaults();

  // force reading MXINPUT_FILE_SIM input files because
  // the primary input file may have an INCLUDE file that
  // is added to INPUTS.INPUT_FILE_LIST.
  for(ifile=0; ifile < MXINPUT_FILE_SIM ; ifile++ ) {
    if ( !IGNOREFILE(INPUTS.INPUT_FILE_LIST[ifile])  ) {
       read_input_file(INPUTS.INPUT_FILE_LIST[ifile] ); 
    }
  }

  // -------
  // check for INCLUDE files on command-line override;
  // read these AFTER reading include files from primary input.  
  for ( i = 2; i < NARGV_LIST ; i++ ) {   
    FOUNDKEY[0] = ( keyMatch(ARGV_LIST[i],"INPUT_FILE_INCLUDE", COLON)  );
    FOUNDKEY[1] = ( keyMatch(ARGV_LIST[i],"INPUT_INCLUDE_FILE", COLON)  );
    if ( FOUNDKEY[0] || FOUNDKEY[1] ) {
      char *include_file = ARGV_LIST[i+1] ;
      read_input_file(include_file);
    }
  } 

  // ------------------------------------------------------
  // check for command line overrides after reading all input files
  // -----------------------------------------------------ß

  sim_input_override(); 

  // check that all command-line args were used
  check_argv(); 

  // check for "PERFECT" options
  genperfect_override();


  // -------------------------------------------
  // make a few checks, compute a few flags and print user input to screen
  // -------------------------------------------

  return;

}  // end of get_user_input



// *********************************
void set_user_defaults(void) {

  int ifilt;
  int i, iz;
  float x;
  double zero = 0.0 ;
  char fnam[] = "set_user_defaults" ;
  // --------------- BEGIN ---------------

  INPUTS.USE_KCOR_REFACTOR = 0 ;
  INPUTS.USE_KCOR_LEGACY   = 1 ;

  INPUTS.DASHBOARD_DUMPFLAG = false ;

  INPUTS.TRACE_MAIN = 0;
  INPUTS.DEBUG_FLAG = 0; 

  INPUTS.RESTORE_DES3YR       = false; // Mar 2020
  INPUTS.RESTORE_HOSTLIB_BUGS = false; // Nov 2019
  INPUTS.RESTORE_FLUXERR_BUGS = false; // Jan 2020
  //INPUTS.RESTORE_WRONG_VPEC   = true ; // Oct 26, 2020 (keep wrong VPEC sign)
  INPUTS.RESTORE_WRONG_VPEC   = false ; // Nov 2, 2020 (fix VPEC sign)
  NLINE_RATE_INFO   = 0;

  // don't init zero'th input file since that is the main input file
  for(i=1; i < MXINPUT_FILE_SIM; i++ ) { INPUTS.INPUT_FILE_LIST[i][0] = 0 ; }
  INPUTS.NREAD_INPUT_FILE = 0 ;

  // - - - - - - - 

  NPAR_ZVAR_USR = USE_ZVAR_FILE = 0;
  sprintf(INPUT_ZVARIATION_FILE,"NONE");
  for ( i=0; i < MXPAR_ZVAR; i++ ) {
    INPUT_ZVARIATION[i].FLAG = -9 ;
    sprintf(INPUT_ZVARIATION[i].PARNAME,"BLANK");
    INPUT_ZVARIATION[i].NZBIN = 0 ;
    init_GENPOLY(&INPUT_ZVARIATION[i].POLY);
    for ( iz=0; iz < MXZBIN_ZVAR; iz++ ) {
      INPUT_ZVARIATION[i].ZBIN_VALUE[iz]    = -9.0 ;
      INPUT_ZVARIATION[i].ZBIN_PARSHIFT[iz] = -9.0 ;
    }
  }

  // - - - - - -
  GENLC.NFILTDEF_OBS = 0;

  INPUTS.ISEED       = 1 ;

  INPUTS.RANLIST_START_GENSMEAR = 1 ;

#ifdef ONE_RANDOM_STREAM
  INPUTS.NSTREAM_RAN = 1 ; // for Mac (7.30.2020
#else
  INPUTS.NSTREAM_RAN = 2 ; // June 6 2020 (2nd stream for spectro noise)
#endif

  INPUTS.NGEN_SCALE         =  1.0 ;
  INPUTS.NGEN_SCALE_NON1A   =  1.0 ;

  INPUTS.NGEN_LC    = 0;
  INPUTS.NGENTOT_LC = 0;
  INPUTS.NGEN_SEASON = 0.0 ;
  INPUTS.INIT_ONLY   = 0 ;  // flag

  INPUTS.NGEN       = 0;
  INPUTS.CIDOFF     = 0;
  INPUTS.CIDRAN_MAX = 3000000 ;
  INPUTS.CIDRAN_MIN = 0 ;
  INPUTS.JOBID      = 0;         // for batch only
  INPUTS.NJOBTOT    = 0;         // for batch only
  INPUTS.NSUBSAMPLE_MARK = 0 ;

  // Mar 2020: use updated cosmoparameters defined in sntools.h
  INPUTS.OMEGA_MATTER  =  (double)OMEGA_MATTER_DEFAULT ; 
  INPUTS.OMEGA_LAMBDA  =  (double)OMEGA_LAMBDA_DEFAULT ;
  INPUTS.w0_LAMBDA     =  (double)w0_DEFAULT ;
  INPUTS.wa_LAMBDA     =  (double)wa_DEFAULT ;
  INPUTS.H0            =  (double)H0_SALT2 ;
  INPUTS.MUSHIFT       =   0.0 ;
  INPUTS.HzFUN_FILE[0] = 0 ;

  INPUTS.GENRANGE_RA[0]   = -360.0  ;
  INPUTS.GENRANGE_RA[1]   = +360.0  ;
  INPUTS.GENRANGE_DEC[0]  = -360.0  ;
  INPUTS.GENRANGE_DEC[1]  = +360.0  ;
  INPUTS.MXRADIUS_RANDOM_SHIFT = 0.0 ;

  INPUTS.SOLID_ANGLE =  0.0  ;
  INPUTS.SNTYPE_Ia_SPEC   =  1 ; // default Ia type 
  INPUTS.SNTYPE_Ia_PHOT   =  INPUTS.SNTYPE_Ia_SPEC + OFFSET_TYPE_PHOT ;

  INPUTS.GENTYPE_SPEC   =  -9 ; 
  INPUTS.GENTYPE_PHOT   =  -9 ;

  INPUTS.GENRANGE_REDSHIFT[0] = 0.0 ;
  INPUTS.GENRANGE_REDSHIFT[1] = 0.0 ;
  INPUTS.GENSIGMA_REDSHIFT    = 0.0 ;
  INPUTS.GENBIAS_REDSHIFT     = 0.0 ;
  INPUTS.GENSIGMA_VPEC        = 0.0 ;
  INPUTS.VPEC_ERR             = 0.0 ;
  INPUTS.VEL_CMBAPEX          = CMBapex_v ; // CMB dipole velocity


  init_RATEPAR ( &INPUTS.RATEPAR );
  init_RATEPAR ( &INPUTS.RATEPAR_PEC1A );
  INPUTS.RATEPAR.SEASON_FRAC           = 1.0 ;
  INPUTS.RATEPAR_PEC1A.SEASON_FRAC     = 0.0 ;

  INPUTS.GENRANGE_PEAKMAG[0] = -99990.0 ;
  INPUTS.GENRANGE_PEAKMAG[1] =  99999.0 ;
  INPUTS.GENRANGE_PEAKMJD[0] = 0.0 ;
  INPUTS.GENRANGE_PEAKMJD[1] = 0.0 ;
  INPUTS.MJD_EXPLODE         = 0.0 ;
  INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE  = -9   ;
  INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX         = -9.0 ;
  INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE  = 0.0 ;
  
  INPUTS.OPT_SETPKMJD      = OPTMASK_SETPKMJD_FLUXMAX2; // May 2019
  INPUTS.MJDWIN_SETPKMJD   = 60.0;  // for default Fmax-clump method
  INPUTS.SNRCUT_SETPKMJD   =  5.0;   // idem

  INPUTS.GENSIGMA_PEAKMJD  = 0.0 ;
  INPUTS.GENRANGE_MJD[0]   = 20000.0 ; // wide open 
  INPUTS.GENRANGE_MJD[1]   = 80000.0 ;
  INPUTS.NEWMJD_DIF  = 0.007 ; // default: same "epoch" for obs within 10'

  INPUTS.GENRANGE_TREST[0]   = 0.0 ;
  INPUTS.GENRANGE_TREST[1]   = 0.0 ;

  INPUTS.GENRANGE_TOBS[0]   = 0.0 ;
  INPUTS.GENRANGE_TOBS[1]   = 0.0 ;

  INPUTS.TGRIDSTEP_MODEL_INTERP = 0.0 ;

  sprintf(INPUTS.STRETCH_TEMPLATE_FILE,"BLANK");

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SHAPEPAR, zero ); 
  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_DELTA,    zero ); 
  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_DM15,     zero ); 
  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_STRETCH,  zero ); 

  INPUTS.DOGEN_AV       = 0 ;
  INPUTS.GENRANGE_AV[0] = 0.0 ;
  INPUTS.GENRANGE_AV[1] = 0.0 ;
  INPUTS.GENEXPTAU_AV   = 0.0 ;
  INPUTS.GENGAUSIG_AV   = 0.0 ;
  INPUTS.GENGAUPEAK_AV  = 0.0 ;
  INPUTS.GENRATIO_AV0   = 0.0 ;
  INPUTS.WV07_GENAV_FLAG    =  0;
  INPUTS.WV07_REWGT_EXPAV   = -9.0;

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_RV, zero );
  INPUTS.GENGAUSS_RV.RANGE[0] = 1.0;  // 2.1 ;
  INPUTS.GENGAUSS_RV.RANGE[1] = 5.0;  // 4.1 ;
  INPUTS.GENGAUSS_RV.PEAK     = RV_MWDUST ; // for SN host

  init_GEN_EXP_HALFGAUSS( &INPUTS.GENPROFILE_AV, (double)-9.0 );
  init_GEN_EXP_HALFGAUSS( &INPUTS.GENPROFILE_EBV_HOST, (double)-9.0 );

  // init SALT2 gen ranges
  INPUTS.GENPOP_ASYMGAUSS_FILE[0] = INPUTS.GENPOP_ASYMGAUSS_MODEL[0] = 0 ;

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2c, zero );
  INPUTS.GENGAUSS_SALT2c.PEAK     = -9.0 ;
  INPUTS.GENGAUSS_SALT2c.SIGMA[0] =  1.0 ;
  INPUTS.GENGAUSS_SALT2c.SIGMA[1] =  1.0 ;
  INPUTS.GENGAUSS_SALT2c.RANGE[0] = -0.5 ;
  INPUTS.GENGAUSS_SALT2c.RANGE[1] =  1.0 ;

  INPUTS.DOGEN_SHAPE = INPUTS.DOGEN_COLOR = true ;

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2x1, zero );
  INPUTS.GENGAUSS_SALT2x1.PEAK     = -9.0 ;
  INPUTS.GENGAUSS_SALT2x1.SIGMA[0] =  5.0 ;
  INPUTS.GENGAUSS_SALT2x1.SIGMA[1] =  5.0 ;
  INPUTS.GENGAUSS_SALT2x1.RANGE[0] = -5.0 ;
  INPUTS.GENGAUSS_SALT2x1.RANGE[1] =  5.0 ;

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2ALPHA, zero );
  INPUTS.GENGAUSS_SALT2ALPHA.RANGE[0] =  0.001; 
  INPUTS.GENGAUSS_SALT2ALPHA.RANGE[1] =  0.40 ; 

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2BETA, zero );
  INPUTS.GENGAUSS_SALT2BETA.RANGE[0] =  0.5 ;
  INPUTS.GENGAUSS_SALT2BETA.RANGE[1] =  9.9 ;

  INPUTS.BIASCOR_SALT2GAMMA_GRID[0] = +9.0 ; // min
  INPUTS.BIASCOR_SALT2GAMMA_GRID[1] = -9.0 ; // max

  init_GENPOLY(&INPUTS.SALT2BETA_cPOLY);

  INPUTS.GENALPHA_SALT2     =  0.0 ; // legacy variable
  INPUTS.GENBETA_SALT2      =  0.0 ; // legacy variable

  // init SIMSED stuff
  INPUTS.NPAR_SIMSED = 0;
  INPUTS.NPAR_SIMSED_PARAM    = 0;
  INPUTS.NPAR_SIMSED_GRIDONLY = 0;
  INPUTS.NPAR_SIMSED_param    = 0;
  INPUTS.NPAR_SIMSED_MODEL    = 0;

  INPUTS.USE_BINARY_SIMSED   = 1; // default is to use binary files
  sprintf(INPUTS.PATH_BINARY_SIMSED,"%s", ".");

  INPUTS.IPAR_SIMSED_SHAPE = -9 ;
  INPUTS.IPAR_SIMSED_COLOR = -9 ;
  INPUTS.IPAR_SIMSED_CL    = -9 ;

  INPUTS.NPAIR_SIMSED_COV = 0 ; 
  INPUTS.OPTMASK_SIMSED   = 0 ;
  for ( i=0; i < MXPAR_SIMSED; i++ ) {
    init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SIMSED[i], (double)(0.0) );
    INPUTS.GENGAUSS_SIMSED[i].RANGE[0] = -1.0E8; // default is wide open
    INPUTS.GENGAUSS_SIMSED[i].RANGE[1] = +1.0E8;

    INPUTS.GENGAUSS_SIMSED[i].SIGMA[0] = 1.0E8; // default is flat distribution
    INPUTS.GENGAUSS_SIMSED[i].SIGMA[1] = 1.0E8;
    INPUTS.GENFLAG_SIMSED[i]  = 0 ; // default is continuous interp

    INPUTS.IPARPAIR_SIMSED_COV[i][0] = -9;
    INPUTS.IPARPAIR_SIMSED_COV[i][1] = -9;

    INPUTS.IROWLIST_SIMSED_COV[i] = -9; 
    INPUTS.IPARLIST_SIMSED_COV[i] = -9; 
  }

  // ----

  INPUTS.SMEARFLAG_FLUX    = 1 ;
  INPUTS.SMEARFLAG_ZEROPT  = 0 ; // Mar 10 2018
  INPUTS.SMEARFLAG_HOSTGAL = 0 ; // internal flag as of Mar 2011 (v9_30G)
  INPUTS.MAGMONITOR_SNR    = 0 ;

  // --- Galactic extinction params
  INPUTS.RV_MWCOLORLAW       = RV_MWDUST ;
  INPUTS.OPT_MWCOLORLAW      = OPT_MWCOLORLAW_ODON94 ; // default
  INPUTS.OPT_MWEBV           = OPT_MWEBV_FILE   ;      // default
  INPUTS.APPLYFLAG_MWEBV     = 0 ;    // default: do NOT correct fluxes
  INPUTS.MWEBV_FLAG          = 1 ;    // default is to do MW extinction
  INPUTS.MWEBV_SIGRATIO      =  0.16;  // default ERR(MWEBV)/MWEBV
  INPUTS.MWEBV_SIG           =  0.0 ;  // add  error in quad on MWEBV (mag)
  INPUTS.MWEBV_SHIFT         =  0.0 ;
  INPUTS.MWEBV_SCALE         =  1.0 ;
  INPUTS.GENRANGE_MWEBV[0]   = -9.0 ;
  INPUTS.GENRANGE_MWEBV[1]   = -9.0 ;

  // ----------

  INPUTS.KCORFLAG_COLOR = GENFRAME_OBS ; // default => use obs-frame closest bands

  INPUTS.GENMAG_OFF_GLOBAL  = 0.0 ;
  INPUTS.GENMAG_OFF_NON1A   = 0.0 ;

  for ( ifilt=0 ; ifilt < MXFILTINDX; ifilt++ ) {
    INPUTS.GENMAG_OFF_MODEL[ifilt]  = 0.0 ;
    INPUTS.GENMAG_OFF_ZP[ifilt]     = 0.0 ;

    INPUTS.TMPOFF_ZP[ifilt]        =  0.0 ;
    INPUTS.TMPOFF_MODEL[ifilt]     =  0.0 ;
    INPUTS.TMPOFF_LAMFILT[ifilt]   =  0.0 ;

    INPUTS.IFILT_SMEAR[ifilt]= 0;
    INPUTS.GENMAG_SMEAR_FILTER[ifilt] = 0.0 ;

   
    INPUTS.MJD_TEMPLATE_FILTER[ifilt] = 0.0 ;
  }
  
  INPUTS.USE_MJD_TEMPLATE     = 0 ;
  INPUTS.MJD_TEMPLATE        = 0.0 ;

  INPUTS.SIGMACLIP_MAGSMEAR[0] = -3.0; // default is +-3 sigma clipping
  INPUTS.SIGMACLIP_MAGSMEAR[1] = +3.0;

  INPUTS.NFILT_SMEAR = 0;
  INPUTS.GENMODEL[0] = 0 ;
  INPUTS.MODELPATH[0] = 0 ;

  INPUTS.GENPDF_FILE[0] = 0 ;
  INPUTS.GENPDF_IGNORE[0] = 0 ;
  INPUTS.GENPDF_FLAT[0]   = 0 ;
  INPUTS.GENPDF_OPTMASK   = 0;
  KEYSOURCE_GENPDF        = -9 ;

  INPUTS.GENMODEL_ERRSCALE     = 0.00 ; // .001 -> 0 (Jun 20 2016) 
  INPUTS.GENMODEL_ERRSCALE_OPT = 1;   // use peak error at all epochs
  INPUTS.GENMODEL_ERRSCALE_CORRELATION = 0.0;   // corr with GENMAG_SMEAR
  INPUTS.GENMODEL_MSKOPT             = 0 ; 
  INPUTS.GENMODEL_ARGLIST[0]         = 0 ;
  INPUTS.GENMAG_SMEAR[0]             = 0.0 ;
  INPUTS.GENMAG_SMEAR[1]             = 0.0 ;
  INPUTS.GENMAG_SMEAR_ADDPHASECOR[0] = 0.0 ;
  INPUTS.GENMAG_SMEAR_ADDPHASECOR[1] = 0.0 ;
  INPUTS.GENSMEAR_RANGauss_FIX       = -999.0 ;
  INPUTS.GENSMEAR_RANFlat_FIX        = -999.0 ;

  INPUTS.GENMAG_SMEAR_SCALE[0] = 0 ;
  sprintf(INPUTS.GENMAG_SMEAR_MODELNAME, "NONE") ;
  INPUTS.GENMAG_SMEAR_MODELARG[0] = 0;
  INPUTS.GENMAG_SMEAR_MSKOPT      = 0 ;

  sprintf(INPUTS.STRONGLENS_FILE,       "NONE");
  sprintf(INPUTS.WEAKLENS_PROBMAP_FILE, "NONE");
  INPUTS.WEAKLENS_DMUSCALE = 1.0 ;
  INPUTS.WEAKLENS_DSIGMADZ = 0.0 ;

  INPUTS.NPAR_GENSMEAR_USRFUN     = 0 ;
  for (i=0; i < 100; i++ ) 
    { INPUTS.GENMAG_SMEAR_USRFUN[i] = -9.0 ; }

  INPUTS.DO_MODELSMEAR   =  0 ;
  NSMEARPAR_OVERRIDE     =  0 ; // see sntools_genSmear.h

  INPUTS.GENFILTERS[0] = 0 ;
  INPUTS.NFILTDEF_OBS = 0;

  INPUTS.GENRANGE_DMPEVENT[0]  = 0 ;
  INPUTS.GENRANGE_DMPEVENT[1]  = 0 ;

  INPUTS.GENRANGE_DMPTREST[0]  = 9999. ;
  INPUTS.GENRANGE_DMPTREST[1]  = 9999. ;
  INPUTS.USEFLAG_DMPTREST = 0 ;

  INPUTS.FORMAT_MASK      = 32 ;                   // 2=TEXT  32=FITS
  INPUTS.WRITE_MASK       = WRITE_MASK_SIM_SNANA ; // default
  INPUTS.WRFLAG_MODELPAR  = 1;  // default is yes
  INPUTS.WRFLAG_YAML_FILE = 0;  // batch-sumbit scripts should set this

  INPUTS.NPE_PIXEL_SATURATE = 1000000000; // billion
  INPUTS.PHOTFLAG_SATURATE = 0 ;
  INPUTS.PHOTFLAG_SNRMAX   = 0 ;
  INPUTS.PHOTFLAG_NEARPEAK = 0 ;

  INPUTS.EXPOSURE_TIME_MSKOPT = 7 ; // bits 1,2,3 => ZPT, SKYSIG, READNOISE
  INPUTS.EXPOSURE_TIME = 1 ;
  for ( ifilt=0 ; ifilt < MXFILTINDX; ifilt++ ) {
    INPUTS.EXPOSURE_TIME_FILTER[ifilt] = 1.0;

    GENLC.SCALE_NOISE[ifilt]      = 1.0 ;
    GENLC.SHIFT_ZPTSIMLIB[ifilt]  = 0.0 ;

    INPUTS.IFILTMAP_OBS[ifilt]   = -9;
    GENLC.IFILTMAP_OBS[ifilt]    = -9;
    GENLC.IFILTINVMAP_OBS[ifilt] = -9;

  }

  INPUTS.GENPERFECT = 0;
  GENPERFECT.NVAR = 0;

  INPUTS.PATH_SNDATA_SIM[0] = 0 ;
  sprintf(INPUTS.GENVERSION,      "BLANK" );  
  sprintf(INPUTS.GENPREFIX,       "BLANK" );
  sprintf(INPUTS.GENSOURCE,       "RANDOM" );
  sprintf(INPUTS.GENMODEL,        "BLANK" );
  INPUTS.GENMODEL_EXTRAP_LATETIME[0] = 0 ;
  sprintf(INPUTS.KCOR_FILE, "NONE" );

  sprintf(INPUTS.GENSNXT, "CCM89" );
  INPUTS.OPT_SNXT = OPT_SNXT_CCM89;

  GENLC.IFLAG_GENSOURCE = 0;

  sprintf(INPUTS.SIMLIB_FILE,      "NONE" );
  sprintf(INPUTS.SIMLIB_SURVEY,    "NONE" );
  sprintf(INPUTS.SIMLIB_FIELDLIST, "ALL" );
  INPUTS.SIMLIB_FIELDSKIP_FLAG = 1 ; //-->count skipped fields in NGENTOT

  INPUTS.SIMLIB_IDSTART  =  0 ;
  INPUTS.SIMLIB_MAXRANSTART =  0 ;
  INPUTS.SIMLIB_IDLOCK   = -9 ;
  INPUTS.SIMLIB_MINOBS   =  1 ; 
  INPUTS.SIMLIB_DUMP     = -9 ;
  INPUTS.SIMLIB_NSKIPMJD =  0 ;
  INPUTS.SIMLIB_NREPEAT  =  1 ;
  INPUTS.NSKIP_SIMLIB    =  0 ;
  INPUTS.SIMLIB_MINSEASON = 0.0 ;

  INPUTS.SIMLIB_CADENCEFOM_ANGSEP     = 0.0 ; // (deg); default is all calc.
  INPUTS.SIMLIB_CADENCEFOM_PARLIST[0] = 0.0 ; // 1st param and parList flag

  INPUTS.USE_SIMLIB_GENOPT   = 0 ; // use gen-options in header
  INPUTS.USE_SIMLIB_REDSHIFT = 0;  // use 'REDSHIFT: <z>' in header
  INPUTS.USE_SIMLIB_DISTANCE = 0;  // use 'DISTANCE: <d(Mpc)>' in header
  INPUTS.USE_SIMLIB_PEAKMJD  = 0;  // use 'PEAKMJD: <t0>'  in header
  INPUTS.USE_SIMLIB_MAGOBS   = 0;
  INPUTS.USE_SIMLIB_SPECTRA  = 0;
  INPUTS.USE_SIMLIB_SALT2    = 0;

  INPUTS.SIMLIB_MSKOPT = 0 ;
  GENLC.SIMLIB_IDLOCK  = -9;

  sprintf(INPUTS.HOSTLIB_FILE,          "NONE" );  // input library
  sprintf(INPUTS.HOSTLIB_WGTMAP_FILE,   "NONE" );  // optional wgtmap
  sprintf(INPUTS.HOSTLIB_ZPHOTEFF_FILE, "NONE" );  // optional zphot-eff
  sprintf(INPUTS.HOSTLIB_SPECBASIS_FILE,"NONE" );  //optional host-spec templ
  sprintf(INPUTS.HOSTLIB_SPECDATA_FILE, "NONE" ); 
  INPUTS.HOSTLIB_STOREPAR_LIST[0] = 0 ; // optional vars -> outfile
  INPUTS.HOSTLIB_PLUS_COMMAND[0]  = 0 ;

  INPUTS.HOSTLIB_USE    = 0;
  INPUTS.HOSTLIB_MSKOPT = 0;
  INPUTS.HOSTLIB_MAXREAD     = MXROW_HOSTLIB;
  INPUTS.HOSTLIB_MNINTFLUX_SNPOS = 0.00 ;  // use 0% as the minimum limit for SNPOS
  INPUTS.HOSTLIB_MXINTFLUX_SNPOS = 0.99 ;  // use 99% of total flux for SNPOS
  INPUTS.HOSTLIB_GALID_NULL      = -9;     // value for no host
  INPUTS.HOSTLIB_GALID_PRIORITY[0]  = 0 ;
  INPUTS.HOSTLIB_GALID_PRIORITY[1]  = 0 ;
  INPUTS.HOSTLIB_GALID_UNIQUE = 0;
  INPUTS.HOSTLIB_SBRADIUS           = 0.6*2.0 ; // arcsec
  INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL  = 9999999;  // default is never re-use host
  INPUTS.HOSTLIB_SCALE_SERSIC_SIZE  = 1.0 ;
  INPUTS.HOSTLIB_SCALE_LOGMASS_ERR  = 1.0 ;

  INPUTS.HOSTLIB_GENZPHOT_OUTLIER[0] = -9.0 ;
  INPUTS.HOSTLIB_GENZPHOT_OUTLIER[1] = -9.0 ;
  for(i=0; i<5; i++) { 
    INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[i] =  0.0 ; 
    INPUTS.HOSTLIB_GENZPHOT_BIAS[i]     =  0.0 ; 
  }
  INPUTS.USE_HOSTLIB_GENZPHOT = 0 ; // logical flag

  INPUTS.HOSTLIB_MAXDDLR = 4.0 ;   // cut on SN-host separation

  HOSTLIB_NBR_WRITE.SEPNBR_MAX = 10.0; // +HOSTNBR keeps neighbors within 10''
  HOSTLIB_NBR_WRITE.NNBR_WRITE_MAX  = 10;   // write up to 10 NBRs
  //  HOSTLIB_NBR.MXCHAR_NBR_LIST = 80;   // max string-length of list

  // define polynom function of ztrue for zSN-zGAL tolerance.
  char CPOLY_DZTOL[] = "0.002,0.04" ;
  init_GENPOLY(&INPUTS.HOSTLIB_GENPOLY_DZTOL);
  parse_GENPOLY(CPOLY_DZTOL, "DZTOL", &INPUTS.HOSTLIB_GENPOLY_DZTOL, fnam);
  
  // debug options
  INPUTS.HOSTLIB_GALID_FORCE   = -9;
  INPUTS.HOSTLIB_ABMAG_FORCE   = -9.0 ;
  INPUTS.HOSTLIB_FIXRAN_RADIUS = -9;
  INPUTS.HOSTLIB_FIXRAN_PHI    = -9;
  INPUTS.HOSTLIB_FIXSERSIC[0]  =  0.0   ; // a
  INPUTS.HOSTLIB_FIXSERSIC[1]  =  0.0   ; // b
  INPUTS.HOSTLIB_FIXSERSIC[2]  = -999.0 ; // n
  INPUTS.HOSTLIB_FIXSERSIC[3]  = -999.0 ; // a_rot, deg

  INPUTS.FLUXERRMODEL_OPTMASK = 0 ;
  INPUTS.FLUXERRMODEL_REDCOV[0] = 0;
  sprintf(INPUTS.FLUXERRMODEL_FILE,          "NONE" ); 
  sprintf(INPUTS.FLUXERRMAP_IGNORE_DATAERR,  "NONE" );
  sprintf(INPUTS.HOSTNOISE_FILE,             "NONE" ); 
  sprintf(INPUTS.WRONGHOST_FILE,             "NONE" ); 

  for ( i=0; i < 2; i++ ) {
    x = ((float)i - 0.5) * 999.* 2.0 ;  // default is no cut 
    INPUTS.HOSTLIB_GENRANGE_NSIGZ[i] = x ;
    INPUTS.HOSTLIB_GENRANGE_RA[i]    = (double)x ;
    INPUTS.HOSTLIB_GENRANGE_DEC[i]   = (double)x ;
  }
  

  INPUTS.OPT_FUDGE_SNRMAX = 0 ;
  INPUTS.FUDGE_SNRMAX     = -9.0 ;
  INPUTS.IFILTOBS_FUDGE_SNRMAX = -1 ;
  INPUTS.STRING_FUDGE_SNRMAX[0] = 0 ;

  INPUTS.FORCEVAL_PSF              = -9.0 ;
  INPUTS.FUDGESCALE_PSF            = 1.0 ;
  INPUTS.FUDGESCALE_NOISE_SKY      = 1.0 ;
  INPUTS.FUDGESCALE_NOISE_READ     = 1.0 ;
  INPUTS.FUDGESCALE_NOISE_TEMPLATE = 1.0 ;
  INPUTS.FUDGESHIFT_ZPT       = 0.0 ;
  INPUTS.FUDGESCALE_FLUXERR   = 1.0 ;
  INPUTS.FUDGESCALE_FLUXERR2  = 1.0 ;
  INPUTS.FUDGEOPT_FLUXERR     = 0 ;
  INPUTS.FUDGE_MAGERR         = 0.0 ; // always applied
  INPUTS.FUDGE_ZPTERR         = 0.0 ; // same, but respects SMEARFLAG_ZPT

  for(ifilt=0; ifilt<MXFILTINDX; ifilt++)  { 
    INPUTS.FUDGE_MAGERR_FILTER[ifilt]        = 0.0; 
    INPUTS.FUDGE_ZPTERR_FILTER[ifilt]        = 0.0; 
    INPUTS.FUDGESHIFT_ZPT_FILTER[ifilt]      = 0.0 ;
    INPUTS.FUDGESHIFT_LAM_FILTER[ifilt]      = 0.0 ;
    INPUTS.FUDGESCALE_FLUXERR_FILTER[ifilt]  = 1.0 ;
    INPUTS.FUDGESCALE_FLUXERR2_FILTER[ifilt] = 1.0 ;
  }

  INPUTS.APPLY_SEARCHEFF_OPT    = 0 ;  // evaluate, but NOT apply

  INPUTS.EFFERR_STOPGEN         = 0.0002 ; // stop when effic error <= this

  // ------
  INPUTS.APPLY_CUTWIN_OPT = 0;

  INPUTS_SEARCHEFF.FIX_EFF_PIPELINE = -9.0 ;
  INPUTS_SEARCHEFF.FUNEFF_DEBUG     = 0 ; // 1->100% eff, 2-> hackFun
  INPUTS_SEARCHEFF.NMAP_DETECT      = 0 ;
  INPUTS_SEARCHEFF.NMAP_PHOTPROB    = 0 ;
  INPUTS_SEARCHEFF.NMAP_SPEC        = 0 ;
  INPUTS_SEARCHEFF.NMAP_zHOST       = 0 ;
  INPUTS_SEARCHEFF.MAGSHIFT_SPECEFF = 0.0 ;
  INPUTS_SEARCHEFF.MINOBS       = 2 ;  // at least 2 obs for search trigger
  INPUTS_SEARCHEFF.PHOTFLAG_DETECT  = 0 ;
  INPUTS_SEARCHEFF.PHOTFLAG_TRIGGER = 0 ;
  INPUTS_SEARCHEFF.OPTMASK_OPENFILE = 0 ;
  sprintf(INPUTS_SEARCHEFF.USER_SPEC_FILE, "NONE");
  sprintf(INPUTS_SEARCHEFF.USER_zHOST_FILE,"NONE");

  //  sprintf(INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE,  "DEFAULT" ); 
  //  sprintf(INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE,    "DEFAULT" ); 

  sprintf(INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE,  "NONE" ); 
  sprintf(INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE,    "NONE" ); 
 
  INPUTS_SEARCHEFF.USER_SPECEFF_SCALE  = 1.0 ; // May 2018
  INPUTS_SEARCHEFF.IFLAG_SPEC_EFFZERO  = 0 ;
  INPUTS_SEARCHEFF.IFLAG_zHOST_EFFZERO = 0 ;
  INPUTS_SEARCHEFF.IVERSION_zHOST      = 0 ;

  INPUTS.EPCUTWIN_LAMREST[0] =  2000.0 ;
  INPUTS.EPCUTWIN_LAMREST[1] = 22000.0 ;

  INPUTS.EPCUTWIN_SNRMIN[0] =  -9999. ;
  INPUTS.EPCUTWIN_SNRMIN[1] =  1.0E9 ;

  INPUTS.CUTWIN_TRESTMIN[0] = -9999. ;
  INPUTS.CUTWIN_TRESTMIN[1] = +9999. ;

  INPUTS.CUTWIN_TRESTMAX[0] = -9999. ;
  INPUTS.CUTWIN_TRESTMAX[1] = +9999. ;

  INPUTS.CUTWIN_TGAPMAX[0] = -999999. ;
  INPUTS.CUTWIN_TGAPMAX[1] = +999999. ;

  INPUTS.CUTWIN_T0GAPMAX[0] = -999999. ;
  INPUTS.CUTWIN_T0GAPMAX[1] = +999999. ;

  INPUTS.CUTWIN_REDSHIFT_TRUE[0] = -9999. ;
  INPUTS.CUTWIN_REDSHIFT_TRUE[1] = +9999. ;

  INPUTS.CUTWIN_REDSHIFT_FINAL[0] = -9999. ;
  INPUTS.CUTWIN_REDSHIFT_FINAL[1] = +9999. ;

  INPUTS.CUTWIN_HOST_ZPHOT[0] = -9999. ;
  INPUTS.CUTWIN_HOST_ZPHOT[1] = +9999. ;

  INPUTS.NCUTWIN_SATURATE   = 0 ;
  INPUTS.NCUTWIN_NOSATURATE = 0 ;

  INPUTS.NCUTWIN_SNRMAX = 0;
  INPUTS.OVERRIDE_CUTWIN_SNRMAX = 0;
  for ( i=0; i < MXCUTWIN_SNRMAX; i++ ) {
    INPUTS.CUTWIN_SNRMAX[i][0] = -9999. ;
    INPUTS.CUTWIN_SNRMAX[i][1] = +9999. ;
    INPUTS.CUTWIN_SNRMAX_NFILT[i]    = 0 ;
    INPUTS.CUTWIN_SNRMAX_TREST[i][0] = -9999. ;
    INPUTS.CUTWIN_SNRMAX_TREST[i][1] = +9999. ;
  }

  // set internal cut 0 for entire light curve so that we
  // get SNRMAX for each filter.
  INPUTS.CUTWIN_SNRMAX_TREST[0][0] = -200.0 ;
  INPUTS.CUTWIN_SNRMAX_TREST[0][1] = +500.0 ;

  INPUTS.CUTWIN_NEPOCH[0] =  1. ;
  INPUTS.CUTWIN_NEPOCH[1] = -9999. ; // SNR requirement

  INPUTS.CUTWIN_NOBSDIF[0] =  0 ;
  INPUTS.CUTWIN_NOBSDIF[1] = 9999 ; 

  INPUTS.CUTWIN_MJDDIF[0] =  0.0 ;
  INPUTS.CUTWIN_MJDDIF[1] = +999999. ; 

  INPUTS.CUTWIN_MWEBV[0] =  0.0 ;
  INPUTS.CUTWIN_MWEBV[1] =  5.0 ;  // Jun e25 2018:  1 -> 5 

  INPUTS.CUTWIN_PEAKMAG[0] =  -30 ;
  INPUTS.CUTWIN_PEAKMAG[1] =  999.0 ; 
  INPUTS.CUTWIN_PEAKMAG_ALL[0] = -30 ;
  INPUTS.CUTWIN_PEAKMAG_ALL[1] =  999.0 ; 

  
  INPUTS.CUTWIN_EPOCHS_SNRMIN      = 999.0 ; 
  INPUTS.CUTWIN_EPOCHS_TRANGE[0]   = 0.0 ;
  INPUTS.CUTWIN_EPOCHS_TRANGE[1]   = 99999.0 ;
  INPUTS.CUTWIN_EPOCHS_FILTERS[0]  = 0 ;
  INPUTS.CUTWIN_EPOCHS_NFILT       = 0 ; 

  INPUTS.NCUTWIN_PEAKMAG_BYFIELD = 0;
  for(i=0; i < MXCUTWIN_PEAKMJD_BYFIELD; i++ ) {
    INPUTS.CUTWIN_PEAKMAG_BYFIELD[i][0] = 0;
    INPUTS.CUTWIN_PEAKMAG_BYFIELD[i][1] = 0;
    INPUTS.CUTWIN_BYFIELDLIST[i][0] = 0 ;
  }


  INPUTS.NON1ASED.PATH[0] = 0 ;
  INPUTS.NON1A_MODELFLAG = -9 ;
  INPUTS.NON1ASED.STOP   = 0 ;
  INPUTS.NON1ASED.NINDEX = 0 ; 
  INPUTS.NON1ASED.NKEY   = 2;
  INPUTS.NON1ASED.NNON1A = 0 ;
  INPUTS.NON1ASED.NPEC1A = 0 ;
  sprintf ( INPUTS.NON1ASED.KEYLIST[1], "INDEX" );
  sprintf ( INPUTS.NON1ASED.KEYLIST[2], "WGT"   );
  for ( i=0; i < MXNON1A_TYPE; i++ ) {
    INPUTS.NON1ASED.INDEX[i]       =  0;
    INPUTS.NON1ASED.WGT[i]         = NULLFLOAT ;
    INPUTS.NON1ASED.MAGOFF[i]      = NULLFLOAT ;
    INPUTS.NON1ASED.MAGSMEAR[i][0] = NULLFLOAT ;
    INPUTS.NON1ASED.MAGSMEAR[i][1] = NULLFLOAT ;
    INPUTS.NON1ASED.SNTAG[i]       = NULLINT ;
    INPUTS.NON1ASED.ISPEC1A[i]     = 0;
    INPUTS.NON1ASED.FLUXSCALE[i]   = 1.0 ;
  }


  // default rise time shifts (relative to model) are zero
  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_RISETIME_SHIFT, zero );

  // repeat for fall-times
  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_FALLTIME_SHIFT, zero );

  // default is to NOT prompt user before clearing (removing) old version 
  INPUTS.CLEARPROMPT      = 0 ;
  INPUTS.REQUIRE_DOCANA   = 1 ;      // set true, July 17 2021
  INPUTS.NVAR_SIMGEN_DUMP = -9 ;    // note that 0 => list variables & quit
  INPUTS.IFLAG_SIMGEN_DUMPALL = 0 ; // dump only SN written to data file.
  INPUTS.PRESCALE_SIMGEN_DUMP = 1 ; // prescale

#ifdef MODELGRID_GEN
  sprintf(GRIDGEN_INPUTS.FORMAT, "%s", "BLANK" );
  for ( i = 0; i <= NPAR_GRIDGEN; i++ ) 
    {  GRIDGEN_INPUTS.NBIN[i] = 0;  }
#endif

  sprintf(INPUTS.NONLINEARITY_FILE,"NONE");

  INPUTS.LCLIB_FILE[0] = 0 ;
  INPUTS.LCLIB_TEMPLATE_EPOCHS[0] = 0 ;
  LCLIB_CUTS.NCUTWIN      = 0 ;
  LCLIB_DEBUG.DAYSCALE             = 1.0 ;
  LCLIB_DEBUG.TOBS_OFFSET_RANGE[0] = 0.0 ;
  LCLIB_DEBUG.TOBS_OFFSET_RANGE[1] = 0.0 ;
  LCLIB_DEBUG.ZERO_TEMPLATE_FLUX   =   0 ;
  LCLIB_DEBUG.FORCE_NREPEAT        =  -9 ;

  set_user_defaults_SPECTROGRAPH();
  set_user_defaults_RANSYSTPAR();

  return ;

}  // end of set_user_defaults

// *******************************************
void set_user_defaults_SPECTROGRAPH(void) {

  // Jun 1 2020: NLAMSIGMA  -> 3.0 (was 5.0)

  // set default spectrograph options
  INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK         =  0 ;
  INPUTS.SPECTROGRAPH_OPTIONS.DOFLAG_SPEC     =  0 ;
  INPUTS.SPECTROGRAPH_OPTIONS.NLAMSIGMA       =  3.0;
  INPUTS.SPECTROGRAPH_OPTIONS.SCALE_LAMSIGMA  =  1. ;
  INPUTS.SPECTROGRAPH_OPTIONS.SCALE_SNR       =  1. ;  // scale on SNR
  INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE   =  1. ;  // scale Texpose
  INPUTS.SPECTROGRAPH_OPTIONS.ILAM_SPIKE      = -9 ;   // lambda bin

  NPEREVT_TAKE_SPECTRUM =  0 ; // Mar 14 2017
  INPUTS.TAKE_SPECTRUM_DUMPCID    = -9 ;
  INPUTS.TAKE_SPECTRUM_HOSTFRAC   =  0.0 ;
  INPUTS.TAKE_SPECTRUM_HOSTSNFRAC =  0.0 ;
  INPUTS.TAKE_SPECTRUM_TEMPLATE_TEXPOSE_SCALE =  1.0 ;

  INPUTS.WARP_SPECTRUM_STRING[0] = 0 ;
  INPUTS.NWARP_TAKE_SPECTRUM = 0 ;
  INPUTS.NHOST_TAKE_SPECTRUM = 0 ;

} // end set_user_defaults_SPECTROGRAPH

// *******************************************
void set_user_defaults_RANSYSTPAR(void) {

  int ifilt ; 

  INPUTS.RANSYSTPAR.USE = 0 ;

  INPUTS.RANSYSTPAR.RANSEED_GEN = 0 ;
  
  // coherent (all bands) scale of flux-errors
  INPUTS.RANSYSTPAR.SIGSCALE_FLUXERR  = 0.0 ; // true & measured 
  INPUTS.RANSYSTPAR.SIGSCALE_FLUXERR2 = 0.0 ; // measured only
  INPUTS.RANSYSTPAR.SIGSCALE_MWEBV    = 0.0 ;
  INPUTS.RANSYSTPAR.SIGSHIFT_MWRV     = 0.0 ;

  INPUTS.RANSYSTPAR.SIGSHIFT_REDSHIFT = 0.0 ; // PA 2020
  INPUTS.RANSYSTPAR.GENMODEL_WILDCARD[0] = 0; // PA 2020
  INPUTS.RANSYSTPAR.GENPDF_FILE_WILDCARD[0] = 0; // Nov 2021

  INPUTS.RANSYSTPAR.SIGSHIFT_OMEGA_MATTER  = 0.0 ;
  INPUTS.RANSYSTPAR.SIGSHIFT_W0            = 0.0 ;

  INPUTS.RANSYSTPAR.RANGESHIFT_OMEGA_MATTER[0]  = 0.0 ;
  INPUTS.RANSYSTPAR.RANGESHIFT_OMEGA_MATTER[1]  = 0.0 ;
  INPUTS.RANSYSTPAR.RANGESHIFT_W0[0]            = 0.0 ;
  INPUTS.RANSYSTPAR.RANGESHIFT_W0[1]            = 0.0 ;

  for(ifilt=0; ifilt < MXFILTINDX; ifilt++ )  { 
    INPUTS.RANSYSTPAR.SIGSHIFT_ZP[ifilt]      = 0.0 ; 
    INPUTS.RANSYSTPAR.SIGSHIFT_LAMFILT[ifilt] = 0.0 ; 
  }

  return ;

} // end set_user_defaults_RANSYSTPAR


// ***********************************
int read_input_file(char *input_file) {

  // Created Jul 2020 by R.Kessler
  //
  // Re-write reading of input so that:
  //  + use same key-parsing for file and command-line overrides
  //  + allow command-line overrids keys to have optional colon
  //     (to match YAML syntax for refactored submit script)

  int MSKOPT = MSKOPT_PARSE_WORDS_FILE + MSKOPT_PARSE_WORDS_IGNORECOMMENT ;
  int  iwd, NWD_FILE, NWD_READ, LENWD, INIT_FLAG_STRING ;
  FILE *fp;
  char tmpWord[MXPATHLEN];  
  char  stringSource[] = "sim-input file" ;
  char fnam[] = "read_input_file" ;

  // ---------- BEGIN ----------

  // unpack ENV, and make sure that input file exists
  ENVreplace(input_file,fnam,1);
  if ( (fp = fopen(input_file, "rt"))==NULL ) {       
    sprintf ( c1err, "Cannot open input file :" );
    sprintf ( c2err," '%s' ", input_file);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  fclose(fp);

  INPUTS.NREAD_INPUT_FILE++ ;
  printf(" --------------------------------------------------------\n");
  INIT_FLAG_STRING = 0 ;
  if ( INPUTS.KEYNAME_DUMPFLAG ) { INIT_FLAG_STRING = -1; }   
  NstringMatch(INIT_FLAG_STRING, STRINGMATCH_INIT, stringSource); // 7.2020

  printf("  Read user input file %d: %s \n", 
	 INPUTS.NREAD_INPUT_FILE, input_file );

  NWD_FILE = store_PARSE_WORDS(MSKOPT,input_file);

  // copy input file contents to another [WORDLIST] buffer because
  // store_PARSE_WORDS may be used later to parse some of the input
  INPUTS.NWORDLIST = NWD_FILE ;
  INPUTS.WORDLIST  = (char**) malloc ( sizeof(char*) * NWD_FILE );
  for(iwd=0; iwd < NWD_FILE; iwd++ ) {    
    get_PARSE_WORD(0,iwd,tmpWord);
    LENWD = strlen(tmpWord) ;
    if ( LENWD < MXCHARWORD_PARSE_WORDS ) { LENWD = MXCHARWORD_PARSE_WORDS; }
    INPUTS.WORDLIST[iwd] = (char*) malloc ( sizeof(char) * (LENWD+10) );
    sprintf(INPUTS.WORDLIST[iwd], "%s", tmpWord );
  }


  for(iwd=0; iwd < NWD_FILE; iwd++ ) {
    NWD_READ = parse_input_key_driver(&INPUTS.WORDLIST[iwd], KEYSOURCE_FILE);
    iwd += NWD_READ ;
  }

  // ------------------- 
  char msg_die[100] ;
  sprintf(msg_die, "DEBUG DIE for refac %s: NWD = %d", fnam, NWD_FILE );
  //  debugexit(msg_die) ;

  if ( INPUTS.KEYNAME_DUMPFLAG ) { happyend(); }

  return(SUCCESS) ;

} // end read_input_file

// ============================
int parse_input_key_driver(char **WORDS, int keySource ) {

  // Jul 20 2020
  // Examine input string(s) WORDS by comparing against all 
  // possible keys. Input keySource indicates if key is
  // file or command-line.
  //
  // Function returns number of words read after the WORDS[0] key.
  // i.e., WORDS[0] is not included in the return count.
  //
  // Apr 6 2021: parse FLUXERRMODEL_REDCOV that was left out of refactor
  // Jun 02 2021: add calls to check_arg_len
  // Jun 17 2021: for hostlib, check "NBR" keys too.
  // Jun 25 2021: check INCLUDE files only if reading file;
  //              command line INCLUDE files are parsed elsewhere.
  //
  // Oct 06 2021: read missing "WV07_REWGT_EXPAV" 
  //

  bool IS_ARG  = (keySource == KEYSOURCE_ARG );
  int j, ITMP, NFILTDEF, NPAR, NFILT, N = 0 ;
  double TMPVAL[2];
  bool ISKEY_INCLUDE, ISKEY_HOSTLIB, ISKEY_SIMLIB, ISKEY_RATE ;
  bool ISKEY_EBV, ISKEY_AV;
  FILE *fpNull = NULL ;
  char strPoly[60], ctmp[60], *parName ;
  char fnam[] = "parse_input_key_driver" ;
  
  // ------------- BEGIN -----------

  // printf(" xxx %s: WORDS = '%s' \n", fnam, WORDS[0] );

  ISKEY_HOSTLIB = (strstr(WORDS[0],"HOSTLIB_") != NULL || 
		   strstr(WORDS[0],"NBR"     ) != NULL );

  ISKEY_SIMLIB  = (strstr(WORDS[0],"SIMLIB_") != NULL );

  ISKEY_EBV     = (strstr(WORDS[0],"_EBV") != NULL );
  ISKEY_AV      = (strstr(WORDS[0],"_AV" ) != NULL );

  ISKEY_RATE    = (strstr(WORDS[0],"DNDZ") != NULL || 
		   strstr(WORDS[0],"DNDB") != NULL  ) ;

  ISKEY_INCLUDE = 
    ( keyMatchSim(2, "INPUT_FILE_INCLUDE", WORDS[0], keySource) ||
      keyMatchSim(2, "INPUT_INCLUDE_FILE", WORDS[0], keySource) );

  // - - - - - - -

  if ( ISKEY_INCLUDE ) {
    if ( keySource == KEYSOURCE_ARG ) { N++; return(N); }
    if ( strlen(INPUTS.INPUT_FILE_LIST[1]) == 0 )  // 1st include file
      { N++; sscanf(WORDS[N], "%s", INPUTS.INPUT_FILE_LIST[1]); }
    else if ( strlen(INPUTS.INPUT_FILE_LIST[2])==0 )  // 2nd include file
      { N++; sscanf(WORDS[N], "%s", INPUTS.INPUT_FILE_LIST[2] ); }
    else if ( strlen(INPUTS.INPUT_FILE_LIST[3])==0 )  // 3rd include file
      { N++; sscanf(WORDS[N], "%s", INPUTS.INPUT_FILE_LIST[3] ); }
    else {
      sprintf(c1err,"Cannot specify %d INPUT_INCLUDE_FILE keys",
	      MXINPUT_FILE_SIM );
      sprintf(c2err,"%d or fewer INCLUDE files allowed.",
	      MXINPUT_FILE_SIM-1 ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }    
  }
  else if ( keyMatchSim(1, "USE_KCOR_REFACTOR", WORDS[0], keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.USE_KCOR_REFACTOR ) ; 
  }

  else if ( keyMatchSim(1, "DASHBOARD", WORDS[0], keySource) ) {
    INPUTS.DASHBOARD_DUMPFLAG = true ; // restore, Mar 9 2021
    N++ ; // no argument, but increment word count to avoid command-line abort
  }

  else if ( keyMatchSim(1, "TRACE_MAIN", WORDS[0], keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.TRACE_MAIN ) ; 
  }
  else if ( keyMatchSim(1, "DEBUG_FLAG", WORDS[0], keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.DEBUG_FLAG) ; 
  }
  else if ( keyMatchSim(1, "RESTORE_DES3YR", WORDS[0], keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &ITMP);  
    if (ITMP) { INPUTS.RESTORE_DES3YR = true; }    
  }
  else if ( keyMatchSim(1, "RESTORE_WRONG_VPEC", WORDS[0], keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &ITMP);  
    INPUTS.RESTORE_WRONG_VPEC = ( ITMP > 0 );
  }
  else if ( ISKEY_HOSTLIB ) {
    N += parse_input_HOSTLIB(WORDS, keySource);
  }
  // - - - - -
  else if ( keyMatchSim(1, "FLUXERRMODEL_FILE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.FLUXERRMODEL_FILE );
  }
  else if ( keyMatchSim(1, "FLUXERRMODEL_OPTMASK", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.FLUXERRMODEL_OPTMASK );
  }

  else if ( keyMatchSim(1, "FLUXERRMODEL_REDCOV", WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], 200);
    char *STR_REDCOV = INPUTS.FLUXERRMODEL_REDCOV;
    N++; sscanf(WORDS[N], "%s", ctmp);
    strcat(STR_REDCOV,WORDS[0] );   // store key name 
    strcat(STR_REDCOV," ");         // blank space 
    strcat(STR_REDCOV,ctmp );   // store argument
    strcat(STR_REDCOV," ");         // blank space  
    if ( IGNOREFILE(ctmp) ) { sprintf(STR_REDCOV,"NONE"); }

  }

  else if ( keyMatchSim(1, "FLUXERRMAP_IGNORE_DATAERR", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.FLUXERRMAP_IGNORE_DATAERR );
  }
  else if ( keyMatchSim(1, "HOSTNOISE_FILE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.HOSTNOISE_FILE );
  }
  else if ( strstr(WORDS[0],"ZVARIATION_") != NULL ) {
    N += parse_input_ZVARIATION(WORDS,keySource);
  }
  // - - - - -
  else if ( keyMatchSim(1, "SNTYPE_Ia",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &ITMP);
    INPUTS.SNTYPE_Ia_SPEC = ITMP ;
    INPUTS.SNTYPE_Ia_PHOT = ITMP + OFFSET_TYPE_PHOT ;
  }
  else if ( keyMatchSim(1, "SNTYPES_Ia",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SNTYPE_Ia_SPEC );
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SNTYPE_Ia_PHOT );
  }
  else if ( keyMatchSim(1, "GENTYPE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &ITMP );
    INPUTS.GENTYPE_SPEC = ITMP ;
    INPUTS.GENTYPE_PHOT = ITMP + OFFSET_TYPE_PHOT ; 
  }
  else if ( keyMatchSim(1, "GENTYPES",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.GENTYPE_SPEC );
    N++;  sscanf(WORDS[N], "%d", &INPUTS.GENTYPE_PHOT );
  }
  else if ( keyMatchSim(1, "NONLINEARITY_FILE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.NONLINEARITY_FILE );
  }
  // - - - - -
  else if ( ISKEY_SIMLIB ) {
    N += parse_input_SIMLIB(WORDS, keySource);
  }
  // - - - - 
  else if ( keyMatchSim(1, "NGEN_LC",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.NGEN_LC );
  }
  else if ( keyMatchSim(1, "NGENTOT_LC",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.NGENTOT_LC );
  }
  else if ( keyMatchSim(1, "NGEN_SEASON",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.NGEN_SEASON );
  }
  else if ( keyMatchSim(1, "NGEN_SCALE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.NGEN_SCALE );
  }
  else if ( keyMatchSim(1, "NGEN_SCALE_NON1A NGEN_SCALE_NONIA",
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.NGEN_SCALE_NON1A );
  }
  else if ( keyMatchSim(1, "NSUBSAMPLE_MARK",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.NSUBSAMPLE_MARK );
  }
  else if ( keyMatchSim(1, "CIDOFF",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.CIDOFF );
  }
  else if ( keyMatchSim(1, "CIDRAN_MAX",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.CIDRAN_MAX );
  }
  else if ( keyMatchSim(1, "CIDRAN_MIN",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.CIDRAN_MIN );
  }
  else if ( keyMatchSim(1, "FORMAT_MASK",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.FORMAT_MASK );
  }
  else if ( keyMatchSim(1, "WRFLAG_MODELPAR",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.WRFLAG_MODELPAR );
  }
  else if ( keyMatchSim(1, "WRFLAG_YAML_FILE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.WRFLAG_YAML_FILE );
  }
  // - - - -
  else if ( keyMatchSim(1, "NPE_PIXEL_SATURATE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.NPE_PIXEL_SATURATE );
  }
  else if ( keyMatchSim(1, "PHOTFLAG_SATURATE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.PHOTFLAG_SATURATE );
  }
  else if ( keyMatchSim(1, "PHOTFLAG_SNRMAX", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.PHOTFLAG_SNRMAX );
  }
  else if ( keyMatchSim(1, "PHOTFLAG_NEARPEAK", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.PHOTFLAG_NEARPEAK );
  }
  else if ( keyMatchSim(1, "PHOTFLAG_DETECT", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS_SEARCHEFF.PHOTFLAG_DETECT );
  }
  else if ( keyMatchSim(1, "PHOTFLAG_TRIGGER", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS_SEARCHEFF.PHOTFLAG_TRIGGER );
  }
  // - - - - - -
  else if ( strstr(WORDS[0],"SIMGEN_DUMP") != NULL ) {
    N += parse_input_SIMGEN_DUMP(WORDS,keySource);
  }
  // - - - - 
  else if ( keyMatchSim(1, "GENVERSION",  WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], 200);
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENVERSION );
    sprintf(INPUTS.GENPREFIX,"%s", INPUTS.GENVERSION);
  }
  else if ( keyMatchSim(1, "GENPREFIX",  WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], 200);
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENPREFIX );
  }
  else if ( keyMatchSim(0, "CLEARPROMPT",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.CLEARPROMPT );
  }
  else if ( keyMatchSim(1, "REQUIRE_DOCANA",  WORDS[0], keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &ITMP );
    // only allow command-line override; not allowed to read from input file
    // to avoid people forgetting
    if ( keySource == KEYSOURCE_ARG ) 
      { INPUTS.REQUIRE_DOCANA = ITMP; }
    else {
      sprintf(c1err,"REQUIRE_DOCANA key not allowed in sim-input file.");
      sprintf(c2err,"Try 'REQUIRE_DOCANA %d' as command-line arg", ITMP);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }
  else if ( keyMatchSim(1, "GENSOURCE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENSOURCE );
  }

  else if ( keyMatchSim(1, "GENMODEL_MSKOPT GENMODEL_OPTMASK",  
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.GENMODEL_MSKOPT );
  }
  else if ( keyMatchSim(1, "GENMODEL_ARGLIST",  WORDS[0],keySource) ) {
    N += parse_input_GENMODEL_ARGLIST(WORDS,keySource);
  }
  else if ( keyMatchSim(1, "GENMODEL",  WORDS[0],keySource) ) {
    N += parse_input_GENMODEL(WORDS,keySource);
  }
  // - - - -
  else if ( keyMatchSim(1, "GENPDF_FILE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENPDF_FILE );
    KEYSOURCE_GENPDF = keySource ; // for prioritization w.r.t. asymGauss
  }
  else if ( keyMatchSim(1, "GENPDF_IGNORE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENPDF_IGNORE );
  }
  else if ( keyMatchSim(1, "GENPDF_FLAT",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENPDF_FLAT );
  }
  else if ( keyMatchSim(1, "GENPDF_OPTMASK",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.GENPDF_OPTMASK );
  }
  else if ( keyMatchSim(1,"GENMODEL_EXTRAP_LATETIME",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENMODEL_EXTRAP_LATETIME );
  }
  // - - - - - PATHs  - - - -
  else if ( keyMatchSim(1, "PATH_USER_INPUT",  WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN);
    N++;  sscanf(WORDS[N], "%s", PATH_USER_INPUT );
  }
  else if ( keyMatchSim(1, "PATH_SNDATA_SIM",  WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN);
    N++;  sscanf(WORDS[N], "%s", INPUTS.PATH_SNDATA_SIM );
  }
  if ( keyMatchSim(1, "PATH_NON1ASED PATH_NONIASED", WORDS[0], keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN);
    N++;  sscanf(WORDS[N], "%s", INPUTS.NON1ASED.PATH );
    return(N);
  }
  else if ( 
	   // check first 5 chars to avoid other strings with NON1A after 1st char
	   strncmp(WORDS[0],"NONIA",5) == 0  ||
	   strncmp(WORDS[0],"NON1A",5) == 0  ||
	   strncmp(WORDS[0],"PECIA",5) == 0  ||
	   strncmp(WORDS[0],"PEC1A",5) == 0 ) {
    N += parse_input_NON1ASED(WORDS,keySource);
  }
  else if ( keyMatchSim(1,"MJD_EXPLODE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.MJD_EXPLODE );
    if(INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE < 0 ) 
      { INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE=0; }
    INPUTS.GENRANGE_PEAKMJD[0] = INPUTS.MJD_EXPLODE ;
    INPUTS.GENRANGE_PEAKMJD[1] = INPUTS.MJD_EXPLODE ;
  }
  else if ( keyMatchSim(1,"OPTMASK_T0SHIFT_EXPLODE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE );
  }
  else if ( keyMatchSim(1,"UVLAM_EXTRAPFLUX", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX );
  }
  else if ( keyMatchSim(1,"MINSLOPE_EXTRAPMAG_LATE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE );
  }
  // - - - -
  else if ( keyMatchSim(1,"RANSEED", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &ITMP);
    INPUTS.ISEED = ITMP; // set unsigned int
  }
  else if ( keyMatchSim(1,"NSTREAM_RAN", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.NSTREAM_RAN );
  } 
  else if ( keyMatchSim(1,"RANLIST_START_GENSMEAR", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.RANLIST_START_GENSMEAR );
  }
  // - - - - DNDZ stuff - - - - 
  else if ( ISKEY_RATE ) {
    N += parse_input_RATEPAR(WORDS, keySource, "NOMINAL",
			     &INPUTS.RATEPAR);
    N += parse_input_RATEPAR(WORDS, keySource, "PEC1A",  
			     &INPUTS.RATEPAR_PEC1A);    
  }
  // - - - -  GENRANGEs - - - - 
  else if ( keyMatchSim(1,"GENRANGE_RA", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_RA[0] );
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_RA[1] );
  }
  else if ( keyMatchSim(1,"GENRANGE_DEC GENRANGE_DECL", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_DEC[0] );
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_DEC[1] );
  }

  else if ( keyMatchSim(1,"MXRADIUS_RANDOM_SHIFT", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.MXRADIUS_RANDOM_SHIFT );
  }

  else if ( strstr(WORDS[0],"SOLID_ANGLE") != NULL ) {
    N += parse_input_SOLID_ANGLE(WORDS,keySource);
  }
  else if ( keyMatchSim(1,"GENRANGE_REDSHIFT", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_REDSHIFT[0] );
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_REDSHIFT[1] );
  }
  else if ( keyMatchSim(1,"GENSIGMA_REDSHIFT", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENSIGMA_REDSHIFT );
  }
  else if ( keyMatchSim(1,"GENBIAS_REDSHIFT", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENBIAS_REDSHIFT );
  }
  else if ( keyMatchSim(1,"GENSIGMA_VPEC", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENSIGMA_VPEC );
  }
  else if ( keyMatchSim(1,"VPEC_ERR", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.VPEC_ERR );
  }
  else if ( keyMatchSim(1,"VEL_CMBAPEX", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.VEL_CMBAPEX );
  }
  else if ( keyMatchSim(1,"WRONGHOST_FILE", WORDS[0],keySource) ) {
    // wrong host -> wrong redshift
    N++;  sscanf(WORDS[N], "%s", &INPUTS.WRONGHOST_FILE);
  }
  else if ( keyMatchSim(1,"GENRANGE_PEAKMAG", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_PEAKMAG[0] );
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_PEAKMAG[1] );
  }
  else if ( keyMatchSim(1,"GENRANGE_MJD", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_MJD[0] );
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_MJD[1] );
  }
  else if ( keyMatchSim(1,"GENRANGE_PEAKMJD", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_PEAKMJD[0] );
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_PEAKMJD[1] );
  }
  else if ( keyMatchSim(1,"GENSIGMA_PEAKMJD GENSIGMA_SEARCH_PEAKMJD", 
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENSIGMA_PEAKMJD );
  }
  else if ( keyMatchSim(1,"OPT_SETPKMJD OPT_SETPEAKMJD", 
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.OPT_SETPKMJD );
  }
  else if ( keyMatchSim(1,"NEWMJD_DIF", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.NEWMJD_DIF );
  }
  else if ( keyMatchSim(1,"GENRANGE_TREST", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_TREST[0] );
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_TREST[1] );
  }
  else if ( keyMatchSim(1,"GENRANGE_TOBS", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_TOBS[0] );
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_TOBS[1] );
  }
  else if ( keyMatchSim(1,"TGRIDSTEP_MODEL_INTERP", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.TGRIDSTEP_MODEL_INTERP );
  }
  // - - - - - SIMSED - - - - 
  else if ( strstr(WORDS[0],"SIMSED") != NULL ) {
    N += parse_input_SIMSED(WORDS,keySource); 
  }
  else if ( keyContains_SIMSED_PARAM(WORDS[0]) ) {
    NPAR = INPUTS.NPAR_SIMSED ;
    ITMP = NPAR-1;
    parName = INPUTS.PARNAME_SIMSED[ITMP];
    N += parse_input_GENGAUSS(parName, WORDS, keySource,
			      &INPUTS.GENGAUSS_SIMSED[ITMP] );
    INPUTS.GENGAUSS_SIMSED[ITMP].FUNINDEX = ITMP ;
  }
  // - - - - - - - - - - 
  // read risetime-shift info
  else if ( strstr(WORDS[0],"TIME_SHIFT") != NULL ) {
    N += parse_input_GENGAUSS("RISETIME_SHIFT", WORDS, keySource,
			      &INPUTS.GENGAUSS_RISETIME_SHIFT );
  
    N += parse_input_GENGAUSS("FALLTIME_SHIFT", WORDS, keySource,
			      &INPUTS.GENGAUSS_FALLTIME_SHIFT );  
  }  
  // - - - -  non-SALT2 SHAPE PDFs - - - - 
  else if ( strstr(WORDS[0],"DELTA") != NULL ) {
    N += parse_input_GENGAUSS("DELTA", WORDS, keySource,
			      &INPUTS.GENGAUSS_DELTA );  
  }
  else if ( strstr(WORDS[0],"DM15") != NULL ) {
    N += parse_input_GENGAUSS("DM15", WORDS, keySource,
			      &INPUTS.GENGAUSS_DM15 );  
  }
  else if ( strstr(WORDS[0],"STRETCH") != NULL ) {
    N += parse_input_GENGAUSS("STRETCH", WORDS, keySource,
			      &INPUTS.GENGAUSS_STRETCH );  
  } 
  else if ( keyMatchSim(1, "STRETCH_TEMPLATE_FILE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.STRETCH_TEMPLATE_FILE );
  }
  // - - - - SALT2 population params - - - - - -

  else if ( keyMatchSim(1, "GENPOP_ASYMGAUSS",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENPOP_ASYMGAUSS_FILE);
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENPOP_ASYMGAUSS_MODEL);
    parse_input_GENPOP_ASYMGAUSS();
  }

  else if ( strstr(WORDS[0],"SALT2c") != NULL ) {
    N += parse_input_GENGAUSS("SALT2c", WORDS, keySource,
			      &INPUTS.GENGAUSS_SALT2c );  
  } 
  else if ( strstr(WORDS[0],"SALT2x1") != NULL ) {
    N += parse_input_GENGAUSS("SALT2x1", WORDS, keySource,
			      &INPUTS.GENGAUSS_SALT2x1 );  
  } 

  else if ( strstr(WORDS[0],"BIASCOR_SALT2ALPHA_GRID") != NULL ) {
    // load asymGauss params for BBC alpha grid
    N++;  sscanf(WORDS[N], "%le", &TMPVAL[0] ) ;
    N++;  sscanf(WORDS[N], "%le", &TMPVAL[1] ) ;
    INPUTS.GENGAUSS_SALT2ALPHA.USE      = true ;
    INPUTS.GENGAUSS_SALT2ALPHA.PEAK     = 0.5*(TMPVAL[0]+TMPVAL[1]);
    INPUTS.GENGAUSS_SALT2ALPHA.RANGE[0] = TMPVAL[0];
    INPUTS.GENGAUSS_SALT2ALPHA.RANGE[1] = TMPVAL[1];
    INPUTS.GENGAUSS_SALT2ALPHA.SIGMA[0] = 1.0E8 ;
    INPUTS.GENGAUSS_SALT2ALPHA.SIGMA[1] = 1.0E8 ;
    INPUTS.GENGAUSS_SALT2ALPHA.NGRID    = 2 ;
    prepIndex_GENGAUSS(WORDS[0], &INPUTS.GENGAUSS_SALT2BETA);    
  }
  else if ( strstr(WORDS[0],"BIASCOR_SALT2BETA_GRID") != NULL ) {
    // load asymGauss params for BBC beta grid 
    N++;  sscanf(WORDS[N], "%le", &TMPVAL[0] ) ;
    N++;  sscanf(WORDS[N], "%le", &TMPVAL[1] ) ;
    INPUTS.GENGAUSS_SALT2BETA.USE      = true ;
    INPUTS.GENGAUSS_SALT2BETA.PEAK     = 0.5*(TMPVAL[0]+TMPVAL[1]);
    INPUTS.GENGAUSS_SALT2BETA.RANGE[0] = TMPVAL[0];
    INPUTS.GENGAUSS_SALT2BETA.RANGE[1] = TMPVAL[1];
    INPUTS.GENGAUSS_SALT2BETA.SIGMA[0] = 1.0E8 ;
    INPUTS.GENGAUSS_SALT2BETA.SIGMA[1] = 1.0E8 ;
    INPUTS.GENGAUSS_SALT2BETA.NGRID     = 2 ;
    prepIndex_GENGAUSS(WORDS[0], &INPUTS.GENGAUSS_SALT2BETA);
  }

  else if ( strstr(WORDS[0],"SALT2ALPHA") != NULL ) {
    N += parse_input_GENGAUSS("SALT2ALPHA", WORDS, keySource,
			      &INPUTS.GENGAUSS_SALT2ALPHA );  
  } 
  else if ( strstr(WORDS[0],"SALT2BETA") != NULL ) {
    N += parse_input_GENGAUSS("SALT2BETA", WORDS, keySource,
			      &INPUTS.GENGAUSS_SALT2BETA );  
  }  
  // - - - - legacy SALT2 alpha,beta keys - - - - - - 
  else if ( keyMatchSim(1, "GENALPHA_SALT2",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENALPHA_SALT2 );
    sprintf(INPUTS.GENGAUSS_SALT2ALPHA.NAME, "SALT2ALPHA" );
    INPUTS.GENGAUSS_SALT2ALPHA.PEAK     = INPUTS.GENALPHA_SALT2 ; 
    INPUTS.GENGAUSS_SALT2ALPHA.USE      = true ;
  }
  else if ( keyMatchSim(1, "GENBETA_SALT2",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENBETA_SALT2 );
    sprintf(INPUTS.GENGAUSS_SALT2BETA.NAME, "SALT2BETA" );
    INPUTS.GENGAUSS_SALT2BETA.PEAK     = INPUTS.GENBETA_SALT2 ; 
    INPUTS.GENGAUSS_SALT2BETA.USE      = true ;
  }
  // - - - - misc SALT2 features - - - - 
  else if ( keyMatchSim(1, "BIASCOR_SALT2GAMMA_GRID",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.BIASCOR_SALT2GAMMA_GRID[0] );
    N++;  sscanf(WORDS[N], "%le", &INPUTS.BIASCOR_SALT2GAMMA_GRID[1] );
  } 
  else if ( keyMatchSim(1, "SALT2BETA_cPOLY",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", strPoly);
    parse_GENPOLY(strPoly, "SALT2BETA_cPOLY", &INPUTS.SALT2BETA_cPOLY, fnam);
  }
  // - - - - - host RV and AV  - - - -
  else if ( strstr(WORDS[0],PARNAME_RV) != NULL ) {
    N += parse_input_GENGAUSS(PARNAME_RV, WORDS, keySource,
			      &INPUTS.GENGAUSS_RV );  
  }

  // - - - - EBV for HOST - - - - 

  else if ( ISKEY_EBV ) {
    N += parse_input_EXP_HALFGAUSS("EBV", WORDS, keySource,
				   &INPUTS.GENPROFILE_EBV_HOST );
    N += parse_input_EXP_HALFGAUSS("EBV_HOST", WORDS, keySource,
				   &INPUTS.GENPROFILE_EBV_HOST );

  }
  else if ( ISKEY_AV ) {
    N += parse_input_EXP_HALFGAUSS("AV", WORDS, keySource,
				   &INPUTS.GENPROFILE_AV );
    N += parse_input_EXP_HALFGAUSS("AV_HOST", WORDS, keySource,
				   &INPUTS.GENPROFILE_AV );
  }

  // - - - - - WV07 options - - - - 
  else if ( keyMatchSim(1, "GENAV_WV07  WV07_GENAV_FLAG", 
			WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%d", &INPUTS.WV07_GENAV_FLAG ) ;
  }
  else if ( keyMatchSim(1, "WV07_REWGT_EXPAV GENAV_REWGT_EXPAV", 
			WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%le", &INPUTS.WV07_REWGT_EXPAV ) ;
  }

  // - - - - MW RV - - - - -

  else if ( keyMatchSim(1, "RV_MWCOLORLAW RVMW", WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%le", &INPUTS.RV_MWCOLORLAW ) ;
  }
  else if ( keyMatchSim(1, "OPT_MWCOLORLAW", WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%d", &INPUTS.OPT_MWCOLORLAW ) ;
  }
  else if ( keyMatchSim(1, "OPT_MWEBV", WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%d", &ITMP );
    INPUTS.OPT_MWEBV = abs(ITMP);
    if ( ITMP < 0 ) { INPUTS.APPLYFLAG_MWEBV=1; } // correct FLUXCAL
    if ( ITMP== 0 ) { INPUTS.APPLYFLAG_MWEBV=0; } // turn off with override

    // Oct 26 2021: no longer allow correcting FLUXCAL for MWEBV
    if ( INPUTS.OPT_MWEBV < 0 ) {
      sprintf(c1err,"Correcting FLUXCAL for MWEBV (OPT_MWEBV=%d)"
	      "no longer allowed.", INPUTS.OPT_MWEBV);
      sprintf(c2err,"OPT_MWEBV must be > 0");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }


  else if ( keyMatchSim(1, "GENSIGMA_MWEBV_RATIO", WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%le", &INPUTS.MWEBV_SIGRATIO ) ;
  }
  else if ( keyMatchSim(1, "GENSIGMA_MWEBV", WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%le", &INPUTS.MWEBV_SIG ) ;
  }
  else if ( keyMatchSim(1, "GENSHIFT_MWEBV FUDGESHIFT_MWEBV", 
			WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%le", &INPUTS.MWEBV_SHIFT ) ;
  }
  else if ( keyMatchSim(1, "GENSCALE_MWEBV FUDGESCALE_MWEBV", 
			WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%le", &INPUTS.MWEBV_SCALE ) ;
  }
  else if ( keyMatchSim(1, "GENRANGE_MWEBV", WORDS[0],keySource) ) {
    N++; sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_MWEBV[0] ) ;
    N++; sscanf(WORDS[N], "%le", &INPUTS.GENRANGE_MWEBV[1] ) ;
  }
  // - - - - host extinction color law model - - - -
  else if ( keyMatchSim(1, "EXTINC_HOSTGAL  COLORLAW_HOSTGAL",  
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.GENSNXT );
  }
  // - - - - -

  else if ( keyMatchSim(1, "GENMAG_OFF_ZP  GENMAG_OFF_ZP", 
			WORDS[0],keySource) ) {
    NFILTDEF = INPUTS.NFILTDEF_OBS ;
    if ( NFILTDEF == 0 ) {
      sprintf(c1err,"Filters NOT specified: cannot read ZP/AB offsets.");
      sprintf(c2err,"Must define filters BEFORE ZP/AB offsets.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    for(j=0; j < NFILTDEF; j++ ) 
      { N++;  sscanf(WORDS[N], "%f", INPUTS.TMPOFF_ZP[j] ); }
  }


  else if ( keyMatchSim(1, "GENMAG_OFF_GLOBAL", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENMAG_OFF_GLOBAL );
  }
  else if ( keyMatchSim(1, "GENMAG_OFF_NON1A", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENMAG_OFF_NON1A );
  }

  else if ( keyMatchSim(1, "GENMAG_OFF_MODEL", WORDS[0],keySource) ) {
    NFILTDEF = INPUTS.NFILTDEF_OBS ;
    if ( NFILTDEF == 0 ) {
      sprintf(c1err,"Filters NOT specified: cannot read MODEL offsets.");
      sprintf(c2err,"Must define filters BEFORE MODEL offsets.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    for(j=0; j < NFILTDEF; j++ ) 
      { N++;  sscanf(WORDS[N], "%f", INPUTS.TMPOFF_MODEL[j] ); }
  }
  else if ( keyMatchSim(1, "GENMODEL_ERRSCALE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENMODEL_ERRSCALE );
  }
  else if ( keyMatchSim(1, "GENMAG_SMEAR", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", ctmp);
    split2floats(ctmp, COMMA, INPUTS.GENMAG_SMEAR );
  }
  else if ( keyMatchSim(1, "GENMAG_SMEAR_ADDPHASECOR", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENMAG_SMEAR_ADDPHASECOR[0] );
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENMAG_SMEAR_ADDPHASECOR[1] );
  }
  else if ( keyMatchSim(1, "GENMAG_SMEAR_USRFUN", WORDS[0],keySource) ) {
    sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"USRFUN");
    INPUTS.NPAR_GENSMEAR_USRFUN     = 8 ; // fix hard-wired param             
    for(j=0; j < INPUTS.NPAR_GENSMEAR_USRFUN; j++ ) 
      { N++ ; sscanf(WORDS[N], "%le", INPUTS.GENMAG_SMEAR_USRFUN ); }    	

    // turn off this option if first USRFUN parameter is negative                    
    if ( INPUTS.GENMAG_SMEAR_USRFUN[0] <= -1.0E-7 ) {
      sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"NONE");
      INPUTS.NPAR_GENSMEAR_USRFUN  = 0 ;
    }
  }
  else if ( strstr(WORDS[0],"GENMAG_SMEAR_SCALE") != NULL ) {
    N += parse_input_GENMAG_SMEAR_SCALE(WORDS,keySource); 
  }
  else if ( keyMatchSim(1, "GENMAG_SMEAR_MSKOPT",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.GENMAG_SMEAR_MSKOPT );
  }

  if ( keyMatchSim(1, "GENMAG_SMEAR_MODELNAME", WORDS[0], keySource) ) {
    char *modelName = INPUTS.GENMAG_SMEAR_MODELNAME  ;
    N++ ; sscanf(WORDS[N], "%s", modelName );
    if ( strcmp(modelName,"G10FUDGE") == 0 )
      { N++ ; sscanf(WORDS[N], "%le", &INPUTS.GENMAG_SMEAR_USRFUN[0] ); }

    if ( strstr(modelName,":")!= NULL)  { parse_GENMAG_SMEAR_MODELNAME(); }
  }
  // - - - - strong and weak lens - - - - 
  else if ( keyMatchSim(1, "STRONGLENS_FILE",  WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN);
    N++;  sscanf(WORDS[N], "%s", &INPUTS.STRONGLENS_FILE );
  }
  else if ( keyMatchSim(1, "WEAKLENS_PROBMAP_FILE  LENSING_PROBMAP_FILE",  
			WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN);
    N++;  sscanf(WORDS[N], "%s", &INPUTS.WEAKLENS_PROBMAP_FILE );
  }
  else if ( keyMatchSim(1, "WEAKLENS_DMUSCALE  LENSING_DMUSCALE",
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.WEAKLENS_DMUSCALE );
  }

  else if ( keyMatchSim(1, "WEAKLENS_DSIGMADZ LENSING_DSIGMADZ",
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.WEAKLENS_DSIGMADZ );
  }
  else if ( keyMatchSim(5, "GENMAG_SMEARPAR_OVERRIDE", WORDS[0],keySource) ) {
    int NVAL;      char key[60], parName[60];
    double tmpList[MXSMEARPAR_OVERRIDE];
    N++ ; sscanf(WORDS[N], "%s", key);  
    NVAL = nval_genSmear_override(key, parName); // return NVAL,parName
    for(j=0; j<NVAL; j++) 
      { N++; sscanf(WORDS[N],"%le", &tmpList[j] ); } // read tmpList 
    store_genSmear_override(parName,NVAL,tmpList);
  }
  else if ( keyMatchSim(1, "GENSMEAR_RANGauss_FIX GENSMEAR_RANGAUSS_FIX",  
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENSMEAR_RANGauss_FIX );
  }
  else if ( keyMatchSim(1, "GENSMEAR_RANFlat_FIX GENSMEAR_RANFLAT_FIX",  
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.GENSMEAR_RANFlat_FIX );
  }
  else if ( keyMatchSim(1, "SIGMACLIP_MAGSMEAR", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.SIGMACLIP_MAGSMEAR[0] );
    N++;  sscanf(WORDS[N], "%f", &INPUTS.SIGMACLIP_MAGSMEAR[1] );
  }
  else if ( keyMatchSim(4, "GENMAG_SMEAR_FILTER", WORDS[0],keySource) ) {
    char filterList[50];  float smear;  int ifilt;
    N++ ; sscanf(WORDS[N] , "%s", filterList );
    N++ ; sscanf(WORDS[N] , "%f", smear );
    NFILT = INPUTS.NFILT_SMEAR ; 
      INPUTS.NFILT_SMEAR += 
	PARSE_FILTLIST(filterList,&INPUTS.IFILT_SMEAR[NFILT+1]);
      for ( j=NFILT+1; j <= INPUTS.NFILT_SMEAR; j++ ) {
        ifilt = INPUTS.IFILT_SMEAR[j];
        INPUTS.GENMAG_SMEAR_FILTER[ifilt] = smear ;
      }
  }
  else if ( keyMatchSim(1, "GENMODEL_ERRSCALE_OPT", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.GENMODEL_ERRSCALE_OPT );
  }
  else if ( keyMatchSim(1, "GENMODEL_ERRSCALE_CORRELATION", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.GENMODEL_ERRSCALE_CORRELATION );
  }
  else if ( strncmp(WORDS[0],"COVMAT_SCATTER",14) == 0 ) {
    char ckey[40];    sprintf(ckey,"%s", WORDS[0] );
    N++ ; sscanf(WORDS[N], "%s", ctmp);
    PARSE_COVMAT_SCATTER(ckey,ctmp) ;
  }
  // - - - filters - - - -
  else if ( keyMatchSim(1, "GENFILTERS", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%s", INPUTS.GENFILTERS );
    INPUTS.NFILTDEF_OBS = PARSE_FILTLIST( INPUTS.GENFILTERS,INPUTS.IFILTMAP_OBS);
  }
  else if ( keyMatchSim(1, "GENRANGE_DMPEVENT", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%d", &INPUTS.GENRANGE_DMPEVENT[0] );
    N++ ; sscanf(WORDS[N], "%d", &INPUTS.GENRANGE_DMPEVENT[1] );
  }
  else if ( keyMatchSim(1, "GENRANGE_DMPTREST", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_DMPTREST[0] );
    N++ ; sscanf(WORDS[N], "%f", &INPUTS.GENRANGE_DMPTREST[1] );
  }
  else if ( keyMatchSim(1, "SMEARFLAG_FLUX", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%d", &INPUTS.SMEARFLAG_FLUX );
  }
  else if ( keyMatchSim(1, "SMEARFLAG_ZEROPT", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%d", &INPUTS.SMEARFLAG_ZEROPT );
  }
  else if ( keyMatchSim(1, "MAGMONITOR_SNR", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%d", &INPUTS.MAGMONITOR_SNR );
  }
  else if ( keyMatchSim(1, "EXPOSURE_TIME_MSKOPT", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%d", &INPUTS.EXPOSURE_TIME_MSKOPT );
  }
  else if ( keyMatchSim(1, "KCOR_FILE", WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN);
    N++ ; sscanf(WORDS[N], "%s", INPUTS.KCOR_FILE );
  }
  else if ( keyMatchSim(1, "OMEGA_MATTER", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%le", &INPUTS.OMEGA_MATTER );
  }
  else if ( keyMatchSim(1, "OMEGA_LAMBDA", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%le", &INPUTS.OMEGA_LAMBDA );
  }
  else if ( keyMatchSim(1, "W0_LAMBDA w0_LAMBDA", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%le", &INPUTS.w0_LAMBDA );
  }
  else if ( keyMatchSim(1, "Wa_LAMBDA wa_LAMBDA", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%le", &INPUTS.wa_LAMBDA );
  }
  else if ( keyMatchSim(1, "H0", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%le", &INPUTS.H0 );
  }
  else if ( keyMatchSim(1, "MUSHIFT", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%le", &INPUTS.MUSHIFT );
  }
  else if ( keyMatchSim(1, "HzFUN_FILE", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%s", INPUTS.HzFUN_FILE );
  }
  // - - - 
  else if ( keyMatchSim(1, "FUDGE_SNRMAX", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%s", INPUTS.STRING_FUDGE_SNRMAX );
    INPUTS.OPT_FUDGE_SNRMAX = 1;
  }
  else if ( keyMatchSim(1, "FUDGE2_SNRMAX", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N], "%s", INPUTS.STRING_FUDGE_SNRMAX );
    INPUTS.OPT_FUDGE_SNRMAX = 2 ;
  }
  else if ( keyMatchSim(1, "FORCEVAL_PSF",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.FORCEVAL_PSF );
  }
  else if ( keyMatchSim(1, "FUDGESCALE_PSF",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.FUDGESCALE_PSF );
  }
  else if ( keyMatchSim(1, "FUDGESCALE_SKYNOISE",  WORDS[0],keySource)||
	    keyMatchSim(1, "FUDGESCALE_NOISE_SKY", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.FUDGESCALE_NOISE_SKY );
  }
  else if ( keyMatchSim(1, "FUDGESCALE_READNOISE",  WORDS[0],keySource)||
	    keyMatchSim(1, "FUDGESCALE_NOISE_READ", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.FUDGESCALE_NOISE_READ );
  }
  else if ( keyMatchSim(1, "FUDGESCALE_NOISE_TEMPLATE",WORDS[0],keySource) ){
    N++;  sscanf(WORDS[N], "%f", &INPUTS.FUDGESCALE_NOISE_TEMPLATE );
  }
  else if ( keyMatchSim(1, "FUDGEOPT_FLUXERR",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.FUDGEOPT_FLUXERR );
  }
  //  - - - - - - - 
  else if ( strstr(WORDS[0], "EXPOSURE_TIME") != NULL ) {
    N += parse_input_KEY_PLUS_FILTER(WORDS, keySource, "EXPOSURE_TIME",
				     &INPUTS.EXPOSURE_TIME, 
				     INPUTS.EXPOSURE_TIME_FILTER) ;
  }
  else if ( strstr(WORDS[0], "FUDGE_MAGERR") != NULL ) {
    N += parse_input_KEY_PLUS_FILTER(WORDS, keySource, "FUDGE_MAGERR",
				     &INPUTS.FUDGE_MAGERR, 
				     INPUTS.FUDGE_MAGERR_FILTER);
  }
  else if ( strstr(WORDS[0], "FUDGE_ZPTERR") != NULL ) {
    N += parse_input_KEY_PLUS_FILTER(WORDS, keySource, "FUDGE_ZPTERR",
				     &INPUTS.FUDGE_ZPTERR, 
				     INPUTS.FUDGE_ZPTERR_FILTER);
  }
  else if ( strstr(WORDS[0], "FUDGESCALE_FLUXERR") != NULL ) {   
    N += parse_input_KEY_PLUS_FILTER(WORDS, keySource, "FUDGESCALE_FLUXERR",
				     &INPUTS.FUDGESCALE_FLUXERR, 
				     INPUTS.FUDGESCALE_FLUXERR_FILTER );
  }
  else if ( strstr(WORDS[0], "FUDGESCALE_FLUXERR2") != NULL ) {   
    N += parse_input_KEY_PLUS_FILTER(WORDS, keySource, "FUDGESCALE_FLUXERR2",
				     &INPUTS.FUDGESCALE_FLUXERR2, 
				     INPUTS.FUDGESCALE_FLUXERR2_FILTER );
  }
  else if ( strstr(WORDS[0], "FUDGESHIFT_ZPT") != NULL ) {   
    N += parse_input_KEY_PLUS_FILTER(WORDS, keySource, "FUDGESHIFT_ZPT",
				     &INPUTS.FUDGESHIFT_ZPT, 
				     INPUTS.FUDGESHIFT_ZPT_FILTER);
  }
  else if ( strstr(WORDS[0], "FUDGESHIFT_LAM") != NULL ) {   
    N += parse_input_KEY_PLUS_FILTER(WORDS, keySource, "FUDGESHIFT_LAM",
				     &INPUTS.FUDGESHIFT_LAM, 
				     INPUTS.FUDGESHIFT_LAM_FILTER);
  }
  else if ( strstr(WORDS[0], "MJD_TEMPLATE") != NULL ) {      
    N += parse_input_KEY_PLUS_FILTER(WORDS, keySource, "MJD_TEMPLATE",
				     &INPUTS.MJD_TEMPLATE, 
				     INPUTS.MJD_TEMPLATE_FILTER);
  }
  // - - - - - - - - 
  else if ( strstr(WORDS[0], "RANSYSTPAR") != NULL ) {  
    N += parse_input_RANSYSTPAR(WORDS, keySource );
  }
  // - - - 
  else if ( keyMatchSim(1, "GENPERFECT",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.GENPERFECT );
  }
  else if ( keyMatchSim(1, "MAGSHIFT_SPECEFF",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS_SEARCHEFF.MAGSHIFT_SPECEFF );
  }
  else if ( keyMatchSim(1, "SEARCHEFF_PIPELINE_LOGIC_FILE", WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN);
    N++;  sscanf(WORDS[N], "%s", INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE );
  }
  else if ( keyMatchSim(1, "SEARCHEFF_PIPELINE_FILE SEARCHEFF_PIPELINE_EFF_FILE",  
			WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE );
  }
  else if ( keyMatchSim(1, "SEARCHEFF_SPEC_FILE  SEARCHEFF_SPECEFF_FILE", 
			WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS_SEARCHEFF.USER_SPEC_FILE );
  }
  else if ( keyMatchSim(1, "SEARCHEFF_SPEC_SCALE SEARCHEFF_SPECEFF_SCALE",
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS_SEARCHEFF.USER_SPECEFF_SCALE );
  }
  else if ( keyMatchSim(1, "SEARCHEFF_zHOST_FILE", WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS_SEARCHEFF.USER_zHOST_FILE );
  }
  else if ( keyMatchSim(1, "APPLY_SEARCHEFF_OPT",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.APPLY_SEARCHEFF_OPT );
  }
  else if ( keyMatchSim(1, "MINOBS_SEARCH",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS_SEARCHEFF.MINOBS );
  }
  else if ( strstr(WORDS[0],"LCLIB") != NULL ) {
    N += parse_input_LCLIB(WORDS,keySource);
  }
  else if ( strstr(WORDS[0],"CUTWIN") != NULL ) {
    N += parse_input_CUTWIN(WORDS,keySource);
  }
  else if ( keyMatchSim(1, "EFFERR_STOPGEN",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.EFFERR_STOPGEN );
  }
  // - - -  specrtrograph - - - -
  else if ( keyMatchSim(1, "SPECTROGRAPH_OPTMASK",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK );
  }
  else if ( keyMatchSim(1, "SPECTROGRAPH_SCALE_TEXPOSE",WORDS[0],keySource)) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE );
  }
  // - - - - TAKE_SPECTRUM - - - - -
  else if ( keyMatchSim(1, "TAKE_SPECTRUM_HOSTFRAC",  WORDS[0], keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.TAKE_SPECTRUM_HOSTFRAC );
  }
  else if ( keyMatchSim(1, "TAKE_SPECTRUM_HOSTSNFRAC",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.TAKE_SPECTRUM_HOSTSNFRAC );
  }
  else if ( keyMatchSim(1, "TAKE_SPECTRUM_PRESCALE",  WORDS[0],keySource) ) {
    char FIELDLIST[60];
    STRING_DICT_DEF *DICT = &INPUTS.DICT_SPECTRUM_FIELDLIST_PRESCALE ;
    N++;  sscanf(WORDS[N], "%s", FIELDLIST);
    parse_string_prescales(FIELDLIST, DICT);
  }
  else if ( keyMatchSim(1, "TAKE_SPECTRUM_DUMPCID",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.TAKE_SPECTRUM_DUMPCID );
  }
  else if ( strstr(WORDS[0],"TAKE_SPECTRUM") ) {  
    N += parse_input_TAKE_SPECTRUM(WORDS, keySource, fpNull );
  }
  else if ( keyMatchSim(10, "WARP_SPECTRUM",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.WARP_SPECTRUM_STRING );
  }

#ifdef MODELGRID_GEN
  else if ( strstr(WORDS[0],"GRID") != NULL ) {
    N += parse_input_GRIDGEN(WORDS,keySource);
  }
#endif
  // - - - - COMMAND-LINE ONLY - - - - - 
  else if ( keyMatchSim(1, "INIT_ONLY",  WORDS[0],keySource) ) {
    // flag to quit after computing number of SN
    N++;  sscanf(WORDS[N], "%d", &INPUTS.INIT_ONLY );
  }
  else if ( keyMatchSim(1, "JOBID",  WORDS[0],keySource) ) {
    // JOBID for batch job
    N++;  sscanf(WORDS[N], "%d", &INPUTS.JOBID );
  }
  else if ( keyMatchSim(1, "NJOBTOT",  WORDS[0],keySource) ) {
    // total number of batch jobs
    N++;  sscanf(WORDS[N], "%d", &INPUTS.NJOBTOT );
  }
  // - - - - ABORT ON OBSOLETE KEYS - - - - 
  else {
    parse_input_OBSOLETE(WORDS,keySource);
  }

  return(N);

} // end parse_input_key_driver


// *******************************************
int parse_input_RATEPAR(char **WORDS, int keySource, char *WHAT, 
			RATEPAR_DEF *RATEPAR ) {

  // Created July 21 2020 [refactor]
  // read input file for DNDZ or DNDB keys.
  // *WHAT = 'NOMINAL' or 'PEC1A'
  //
  // Oct 16 2020: abort on multiple rate models. Note that same model
  //              key name can appear more than once; e.g., POWERLAW2.

  bool  IS_NOMINAL = (strcmp(WHAT,"NOMINAL") == 0 ) ;
  bool  IS_PEC1A   = (strcmp(WHAT,"PEC1A"  ) == 0 ) ;

  bool FOUND_PRIMARY_KEY, CONTINUE ;
  int  N=0, n, j, NLOCAL, nread, NMODEL_LIST ;
  double l=0.0, b=0.0, bmax, R=0.0, TMPVAL ;
  char KEYNAME[40], TMPNAME[60];
  char fnam[] = "parse_input_RATEPAR" ;

  // ------------ BEGIN ------------

  sprintf(KEYNAME, "%s", WORDS[0] );
  CONTINUE = false ;
  if ( strstr(KEYNAME,"DNDZ") != NULL ) { CONTINUE = true ; }
  if ( strstr(KEYNAME,"DNDB") != NULL ) { CONTINUE = true ; }
  if ( !CONTINUE  && !INPUTS.KEYNAME_DUMPFLAG ) { return(N) ; }

  // check a few misc keys
  if ( IS_NOMINAL ) {

    if ( keyMatchSim(1, "DNDZ_ZEXP_REWGT", KEYNAME, keySource) ) {
      N++; sscanf(WORDS[N], "%le", &RATEPAR->DNDZ_ZEXP_REWGT ); 
    }
    else if ( keyMatchSim(1, "DNDZ_ZPOLY_REWGT", KEYNAME, keySource) ) {
      N += read_genpoly(KEYNAME, &WORDS[N+1], 3, &RATEPAR->DNDZ_ZPOLY_REWGT);
    }
    else if ( keyMatchSim(1, "DNDZ_SCALE", KEYNAME, keySource) ) {
      N++ ; sscanf(WORDS[N], "%le", &RATEPAR->DNDZ_SCALE[0] ); 
      N++ ; sscanf(WORDS[N], "%le", &RATEPAR->DNDZ_SCALE[1] ); 
    }
    else if ( keyMatchSim(1, "DNDZ_ALLSCALE", KEYNAME, keySource) ) {
      N++; sscanf(WORDS[N], "%le", &RATEPAR->DNDZ_ALLSCALE ); 
    }
    else if ( keyMatchSim(1, "DNDZ_SCALE_NON1A", KEYNAME, keySource) ||
	      keyMatchSim(1, "DNDZ_SCALE_NONIA", KEYNAME, keySource) ) {
      N++; sscanf(WORDS[N], "%le", &RATEPAR->DNDZ_SCALE[1] ); 
    }

    // return if any key above was read
    if ( N > 0 ) { return(N); }
  }

  // -----------------------------------

  FOUND_PRIMARY_KEY = valid_DNDZ_KEY(WHAT, keySource, KEYNAME ) ;

  /*
  printf(" xxx %s: FOUND=%d  IS_NOM/PEC1A = %d/%d \n",
	 fnam, FOUND_PRIMARY_KEY,  IS_NOMINAL, IS_PEC1A );
  */

  // --------------------

  if ( FOUND_PRIMARY_KEY ) {

    N++; sscanf(WORDS[N], "%s", TMPNAME ); 

    // abort on multiple rate models
    if ( !IGNOREFILE(RATEPAR->NAME) ) { 
      if ( strcmp(RATEPAR->NAME,TMPNAME) != 0 ) {
	sprintf(c1err,"cannot mix multiple rate models with %s", KEYNAME);
	sprintf(c2err,"%s and %s (pick one)", RATEPAR->NAME, TMPNAME);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
    }
    
    sprintf(RATEPAR->NAME, "%s", TMPNAME);

    // make sure that PEC1A rateModel is defined only for NON1A mode
    if ( IS_PEC1A &&  INPUTS.NON1A_MODELFLAG < 0 ) {
      sprintf(c1err,"%s key allowed only for NON1A model.", KEYNAME);
      sprintf(c2err,"Check sim-input file." );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( strcmp(RATEPAR->NAME,"ABMODEL") == 0 ) {
	RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_AB ;
	RATEPAR->NMODEL_ZRANGE = 1 ;  
	N++; sscanf(WORDS[N], "%le", &RATEPAR->MODEL_PARLIST[1][0] ); 
	N++; sscanf(WORDS[N], "%le", &RATEPAR->MODEL_PARLIST[1][1] ); 
    }
    else if ( strcmp(RATEPAR->NAME,"POWERLAW") == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_POWERLAW ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      NLOCAL = RATEPAR->NMODEL_ZRANGE ;
      N++; sscanf(WORDS[N], "%le", &RATEPAR->MODEL_PARLIST[NLOCAL][0] ); 
      N++; sscanf(WORDS[N], "%le", &RATEPAR->MODEL_PARLIST[NLOCAL][1] ); 
    }
    else if ( strcmp(RATEPAR->NAME,"POWERLAW2") == 0 ) {
      RATEPAR->NMODEL_ZRANGE++ ;
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_POWERLAW2 ;
      NLOCAL = RATEPAR->NMODEL_ZRANGE ;

      N++; nread=sscanf(WORDS[N], "%le", &RATEPAR->MODEL_PARLIST[NLOCAL][0] ); 
      N++; nread=sscanf(WORDS[N], "%le", &RATEPAR->MODEL_PARLIST[NLOCAL][1] ); 
      N++; nread=sscanf(WORDS[N], "%le", &RATEPAR->MODEL_ZRANGE[NLOCAL][0] ); 
      N++; nread=sscanf(WORDS[N], "%le", &RATEPAR->MODEL_ZRANGE[NLOCAL][1] ); 
      if(nread!=1) { abort_bad_input(KEYNAME, WORDS[N], 3, fnam); }
    }
    else if ( strstr(RATEPAR->NAME,RATEMODELNAME_CCS15) != NULL ) {
      parse_multiplier(RATEPAR->NAME,RATEMODELNAME_CCS15, &TMPVAL);
      sprintf(RATEPAR->NAME,"%s", RATEMODELNAME_CCS15); // strip off scale
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_CCS15 ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      RATEPAR->MODEL_PARLIST[1][0] = TMPVAL; // rate-scale
    }
    else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_PISN) == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_PISN ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
    }
    else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_TDE) == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_TDE ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      // read rate at z=0
      N++; sscanf(WORDS[N], "%le", &RATEPAR->MODEL_PARLIST[1][0] ); 
    }
    else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_MD14) == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_MD14 ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      // read rate at z=0
      N++; sscanf(WORDS[N], "%le", &RATEPAR->MODEL_PARLIST[1][0] ); 
    }
    else if ( strcmp(RATEPAR->NAME,"ZPOLY") == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_ZPOLY ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      N += read_genpoly(RATEPAR->NAME, &WORDS[N+1], 3, 
			&RATEPAR->MODEL_ZPOLY);
    }
    else if ( strcmp(RATEPAR->NAME,"HUBBLE") == 0 ) {
      // Jun 20 2016: set powerlaw model with alpha=0 to avoid abort later
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_POWERLAW ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      RATEPAR->MODEL_PARLIST[1][0] =  1.0 ; // crazy rate per Mpc^3 per year  
      RATEPAR->MODEL_PARLIST[1][1] =  0.0 ; // alpha for (1+z)^alpha 
      RATEPAR->MODEL_ZRANGE[1][0]  =  0.0 ;  // Zmin
      RATEPAR->MODEL_ZRANGE[1][1]  =  ZMAX_SNANA ;  // Zmax
    }
    else if ( strcmp(RATEPAR->NAME,"FLAT") == 0 ) {
      RATEPAR->INDEX_MODEL   = INDEX_RATEMODEL_FLAT ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
    }
    else if ( strcmp(RATEPAR->NAME,"COSBPOLY") == 0 ) {
      RATEPAR->INDEX_MODEL   = INDEX_RATEMODEL_COSBPOLY ;
      RATEPAR->NMODEL_ZRANGE = 0 ;
      N += read_genpoly(RATEPAR->NAME, &WORDS[N+1], 5,
			&RATEPAR->MODEL_BPOLY );     
      // get max rate (vs. b) for weighting
      for(b=0.0; b < 90.0; b+=1.0 ) {
	R = GALrate_model(l,b, RATEPAR);
	if ( R > RATEPAR->RATEMAX ) { RATEPAR->RATEMAX=R; bmax=b; }
      }
    }
    else if ( strcmp(RATEPAR->NAME,"BPOLY") == 0 ) {
      RATEPAR->INDEX_MODEL   = INDEX_RATEMODEL_BPOLY ;
      RATEPAR->NMODEL_ZRANGE = 0 ;
      N += read_genpoly(RATEPAR->NAME, &WORDS[N+1], 5,
			&RATEPAR->MODEL_BPOLY );

      // get max rate (vs. b) for weighting
      for(b=0.0; b < 90.0; b+=1.0 ) {
	R = GALrate_model(l,b, RATEPAR);
	if ( R > RATEPAR->RATEMAX ) { RATEPAR->RATEMAX=R; bmax=b; }
      }

    }
    else {
      sprintf(c1err,"'%s %s' is invalid", KEYNAME, RATEPAR->NAME );
      sprintf(c2err,"Check sim-input file '%s'", 
	      INPUTS.INPUT_FILE_LIST[0] ); 
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    
  } // DNDZ

  return(N);

} // end parse_input_RATEPAR


// ====================================
bool valid_DNDZ_KEY(char *WHAT, int keySource, char *KEYNAME ) {

  // Created Aug 2016
  // Return True if this is a valid DNDZ key.
  //
  // Inputs:
  //   *WHAT     = 'NOMINAL' or 'PEC1A'
  //   keySource  = KEYSOURCE_FILE ==> read from file with colon after key
  //   keySource  = KEYSOURCE_ARG  ==> command-line override
  //   KEYNAME   = name of key to test
  //
  // Jul 21 2020: use keyMatchSim(...) and keySource

  bool FROM_FILE = ( keySource == KEYSOURCE_FILE ) ;
  int  ISNOMINAL = strcmp(WHAT,"NOMINAL") == 0 ;
  int  ISPEC1A   = strcmp(WHAT,"PEC1A"  ) == 0 ;
  char PRIMARY_KEYLIST[10][20], KEYTEST[20] ;
  int  ikey, NKEY ;
  bool FOUND = false ;
  char fnam[] = "valid_DNDZ_KEY" ;

  // ----------- BEGIN -----------


  NKEY = 0 ;
  if ( ISNOMINAL ) { 
    sprintf(PRIMARY_KEYLIST[NKEY],"DNDZ");        NKEY++ ;
    sprintf(PRIMARY_KEYLIST[NKEY],"DNDZ_NON1A");  NKEY++ ;
    sprintf(PRIMARY_KEYLIST[NKEY],"DNDZ_NONIA");  NKEY++ ;
    sprintf(PRIMARY_KEYLIST[NKEY],"DNDB");        NKEY++ ; // Nov 26 2017
  }
  else if ( ISPEC1A ) {
    sprintf(PRIMARY_KEYLIST[NKEY],"DNDZ_PEC1A");  NKEY++ ;
    sprintf(PRIMARY_KEYLIST[NKEY],"DNDZ_PECIA");  NKEY++ ;
  }
  else {
    sprintf(c1err,"Invalid WHAT = '%s'", WHAT);
    sprintf(c2err,"Must be 'NOMINAL'  or 'PEC1A'" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // check for primary key
  for(ikey=0; ikey < NKEY; ikey++ ) {
    sprintf(KEYTEST,"%s", PRIMARY_KEYLIST[ikey] );
    if ( keyMatchSim(2, KEYTEST, KEYNAME, keySource) ) { FOUND=true; }
  } // end ikey

  return(FOUND);

} // end valid_DNDZ_KEY 




// *******************************************
void parse_input_FIXMAG(char *string) {

  // Created Sep 12 2015
  // if string is a number, then set 
  // INPUTS.FIXMAG[0] = INPUTS.FIXMAG[1] = number
  //
  // if string is mag0:mag1 then set
  //   INPUTS.FIXMAG[0] = mag0
  //   INPUTS.FIXMAG[1] = mag1
  //

  int i, LEN, jsplit ;
  char c1[2], strMag0[20], strMag1[20], strTmp[20];
  //  char fnam[] = "parse_input_FIXMAG" ;

  // ------------- BEGIN ---------------

  LEN = strlen(string);
  jsplit = 0 ;

  strMag0[0] = 0 ;
  strMag1[0] = 0 ;

  for(i=0; i < LEN; i++ ) {
    sprintf(c1, "%c", string[i] );
    if ( *c1 == ':' ) { jsplit = i;  continue ; }
    
    if ( jsplit == 0 )  {
      sprintf(strTmp,"%s", strMag0);
      strcat(strMag0,c1);
    }
    //     { sprintf(strMag0,"%s%s", strMag0, c1); }
    else {
      sprintf(strTmp,"%s", strMag1);
      strcat(strMag1,c1);
    }
    //      { sprintf(strMag1,"%s%s", strMag1, c1); }
  }
 
  if ( strlen(strMag1) == 0 )  { sprintf(strMag1,"%s", strMag0); }

  sscanf(strMag0, "%le", &INPUTS.FIXMAG[0] ); 
  sscanf(strMag1, "%le", &INPUTS.FIXMAG[1] ); 

  return ;

} // end parse_input_FIXMAG


// ==============================================================
void  parse_input_GENZPHOT_OUTLIER(char *string) {

  // Created May 10 2018
  //
  // Check 5th element of HOSTLIB_GENZPHOT_FUDGEPAR,
  // which is the sigma of 2nd Gaussian. This element can
  // be a number, or a FLAT string with z-range in ().
  //
  // If *string is a float, then set
  // INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[5] = string.
  //
  // If *string = FLAT(z1:z2) then set
  // INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[4] = 99
  // INPUTS.HOSTLIB_GENZPHOT_OUTLIER[0:1] = z1:z2
  //

  int   NSPLIT;
  int   LDMP = 0 ;
  float *LOADPAR = &INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[4];
  float *zrange  = INPUTS.HOSTLIB_GENZPHOT_OUTLIER ;
  char string_z[40];
  char string_zsplit[2][20], *ptr_z[2];
  char colon[] = ":" ;
  char fnam[] = "parse_input_GENZPHOT_OUTLIER" ;

  // --------------- BEGIN -----------

  if ( strstr(string,"FLAT") == NULL ) {
    sscanf(string, "%f", LOADPAR ); 
  }
  else {
    extractStringOpt(string,string_z);
    *LOADPAR = 99.0 ;  // large sigma --> flat
    ptr_z[0] = string_zsplit[0]; 
    ptr_z[1] = string_zsplit[1]; 
    splitString(string_z, colon, 3,  // inputs    
		&NSPLIT, ptr_z );    // outputs             

    // load zphot outlier range 
    sscanf(ptr_z[0], "%f", &zrange[0] );
    sscanf(ptr_z[1], "%f", &zrange[1] );

    if ( LDMP ) {
      printf(" xxx string_z = '%s'  zrange = %.2f to %.2f \n",  
	     string_z,
	   INPUTS.HOSTLIB_GENZPHOT_OUTLIER[0],
	   INPUTS.HOSTLIB_GENZPHOT_OUTLIER[1] );
    }
  }

  if ( LDMP ) {
    printf(" xxx LOADPAR=%f \n", *LOADPAR);
    debugexit(fnam); 
  }

  return ;

} // end parse_input_GENZPHOT_OUTLIER

// =========================================
int parse_input_NON1ASED(char **WORDS, int keySource) { 

  int  NKEY, key, N=0, NN;
  float sigTmp[2] ;
  bool IS_NON1AKEY, IS_PEC1AKEY, FOUND_KEY[4], DO_NEXT_READ ; 
  char ckey[40], ctmp[40];
  char fnam[] = "parse_input_NON1ASED" ;

  // ----------- BEGIN -----------

  
  if ( INPUTS.NON1A_MODELFLAG != MODEL_NON1ASED ) { return(N); }

  // below works only for keys inside input file


  // check for list of keys to define columms
  IS_NON1AKEY = 
    NstringMatch(1,WORDS[0],"NON1A_KEYS:")    ||
    NstringMatch(1,WORDS[0],"NONIA_KEYS:")    ||
    NstringMatch(1,WORDS[0],"NON1ASED_KEYS:") ||
    NstringMatch(1,WORDS[0],"NONIASED_KEYS:") ;
  
  if ( IS_NON1AKEY ) {
    if ( keySource == KEYSOURCE_ARG ) {
      sprintf(c1err,"Cannot read NON1A keys on command line");
      sprintf(c2err,"NON1A_KEYS and NON1ASED belong in the input file");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    N++; sscanf(WORDS[N], "%d", &NKEY) ;
    if ( NKEY < 2 ) {
      sprintf(c1err,"NON1A_KEYS = %d, but must be >=2.", NKEY);
      sprintf(c1err,"Check sim-input file." );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    INPUTS.NON1ASED.NKEY = NKEY;
    for ( key = 1; key <= NKEY; key++ ) {
      N++; sscanf(WORDS[N], "%s", ckey) ;
      sprintf(INPUTS.NON1ASED.KEYLIST[key], "%s", ckey );
    }

    // check required keys
    FOUND_KEY[1] = (strcmp(INPUTS.NON1ASED.KEYLIST[1],"INDEX")==0 );
    FOUND_KEY[2] = (strcmp(INPUTS.NON1ASED.KEYLIST[2],"WGT"  )==0 ) ;
    if ( !FOUND_KEY[1] ) {
      sprintf(c1err,"First NON1ASED_KEY is '%s' , but must be 'INDEX' ",
	      INPUTS.NON1ASED.KEYLIST[1] ) ;
      sprintf(c2err,"Check sim-input file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    if ( !FOUND_KEY[2] ) {
      sprintf(c1err,"Second NON1ASED_KEY is '%s' , but must be 'WGT' ", 
	      INPUTS.NON1ASED.KEYLIST[2]) ;
      sprintf(c2err,"Check sim-input file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    
    return(N);

  } // end of NON1ASED_KEYS



  if ( NstringMatch(1,WORDS[0],"NON1A_STOP:") ||
       NstringMatch(1,WORDS[0],"NON1ASED_STOP:") ) 
    { INPUTS.NON1ASED.STOP = 1 ;  return(N); }

  // check for individual NONIA[SED] keys for each template

  IS_NON1AKEY  = 
    NstringMatch(MXNON1A_TYPE,WORDS[0],"NON1A:")    ||
    NstringMatch(MXNON1A_TYPE,WORDS[0],"NONIA:")    ||
    NstringMatch(MXNON1A_TYPE,WORDS[0],"NONIASED:") ||
    NstringMatch(MXNON1A_TYPE,WORDS[0],"NON1ASED:") ;
   
  IS_PEC1AKEY  = 
    NstringMatch(MXNON1A_TYPE,WORDS[0],"PEC1A:")    || 
    NstringMatch(MXNON1A_TYPE,WORDS[0],"PECIA:")    ||
    NstringMatch(MXNON1A_TYPE,WORDS[0],"PECIASED:") ||
    NstringMatch(MXNON1A_TYPE,WORDS[0],"PEC1ASED:") ;


  DO_NEXT_READ = ( (IS_NON1AKEY || IS_PEC1AKEY) && !INPUTS.NON1ASED.STOP ) ;

    
  if ( DO_NEXT_READ ) { 

      INPUTS.NON1ASED.NINDEX++ ; NN = INPUTS.NON1ASED.NINDEX;
      if ( IS_PEC1AKEY ) { 
	INPUTS.NON1ASED.NPEC1A++ ;
	INPUTS.NON1ASED.ISPEC1A[NN] = 1; 
      }
      else { 
	INPUTS.NON1ASED.NNON1A++ ; 
	if ( INPUTS.NON1ASED.NPEC1A > 0 ) {
	  sprintf(c1err, "All NON1A keys must appear before PEC1A keys");
	  sprintf(c2err, "Check sim-input file (and INCLUDE file)");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
	}
      }

      for ( key=1; key <= INPUTS.NON1ASED.NKEY ; key++ ) {

	sprintf(ckey, "%s", INPUTS.NON1ASED.KEYLIST[key] );
	N++; sscanf(WORDS[N], "%s", ctmp);

	// store float value for README file
	sscanf(ctmp,"%f", &INPUTS.NON1ASED.KEYVAL[NN][key] );

	if ( strcmp(ckey, "INDEX") == 0 ) 
	  { sscanf(ctmp, "%d", &INPUTS.NON1ASED.INDEX[NN]); }
	else if ( strcmp(ckey, "WGT") == 0 ) 
	  { sscanf(ctmp, "%f", &INPUTS.NON1ASED.WGT[NN]); }
	else if ( strcmp(ckey, "MAGOFF") == 0 ) 
	  { sscanf(ctmp, "%f", &INPUTS.NON1ASED.MAGOFF[NN]); }
	else if ( strcmp(ckey, "MAGSMEAR") == 0 ) {
	  split2floats(ctmp, COMMA, sigTmp ); 
	  INPUTS.NON1ASED.MAGSMEAR[NN][0] = sigTmp[0] ;
	  INPUTS.NON1ASED.MAGSMEAR[NN][1] = sigTmp[1] ;
	  // if non1a magsmear is nonzero, then set global GENMAG_SMEAR
	  // so that random numbers are generated for smearing.
	  if ( sigTmp[0] > 0.0  &&  INPUTS.GENMAG_SMEAR[0] == 0.0 ) 
	    { INPUTS.GENMAG_SMEAR[0] = INPUTS.GENMAG_SMEAR[1] = 0.001; }
	}
	else if ( strcmp(ckey, "SNTYPE") == 0 ) 
	  { sscanf(ctmp, "%d", &INPUTS.NON1ASED.SNTAG[NN]); }
	else if ( strcmp(ckey, "SNTAG") == 0 )   // same as SNTYPE
	  { sscanf(ctmp, "%d", &INPUTS.NON1ASED.SNTAG[NN]); }
	else {
	  sprintf(c1err, "Invalid NON1A key: '%s'", ckey );
	  sprintf(c2err, "Valid keys are INDEX WGT MAGOFF MAGSMEAR" );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
	}
      }

      return(N);
  }

  return(N);

} // end parse_input_NON1ASED

// =================================================
int parse_input_KEY_PLUS_FILTER(char **WORDS, int keySource, char *KEYCHECK,
				float *VALUE_GLOBAL, float *VALUE_FILTERLIST){
				

  // Created Jul 2020 [refactor]
  // WORDS is list of strings read from input file.
  // KEYCHECK is the key to check.
  //
  // If INPUT_STRING == KEYCHECK,        then fill VALUE_GLOBAL.
  // If INPUT_STRING == KEYCHECK_FILTER, then fill VALUE_FILTERLIST
  //
  // Sim-input is of the form
  //    INPUT_STRING:  VALUE_GLOBAL   # same value for all bands
  //        or
  //    INPUTS_STRING_FILTER:  <filterStringList>  <values>
  //
  // Examples:
  //  +  FUDGESHIFT_LAM_FILTER: gr 10
  //      --> g & r bands have 10 A shift
  //  +  FUDGESHIFT_LAM_FILTER: g,r,i,z  -6.4,3.2,6.9,0.4
  //      --> g,r,i,z bands each have respective shifts
  // 
  // Note that outputs are float, not double.
  // Functions returns number of filled values.
  //
  // Jun 22 2021: allow comma-sep list to specify multiple bands
  //              with one key
  //
  float ftmp, val_list[MXFILTINDX];
  char cfilt[MXFILTINDX], cval[200], **cval_list, KEY[80] ;
  int  MAXKEY = 10;
  int  NFILT, NVAL, ifilt_obs, ifilt, ifilt_list[MXFILTINDX];
  char fnam[] = "parse_input_KEY_PLUS_FILTER" ;

  // ----------- BEGIN -----------

  // read from command-line args   
  sprintf(KEY,"%s", KEYCHECK);

  if ( keyMatchSim(1, KEY, WORDS[0], keySource) ) {
    sscanf(WORDS[1] , "%f", VALUE_GLOBAL ); 
    for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ ) 
      { VALUE_FILTERLIST[ifilt] = *VALUE_GLOBAL ; }
    return(1);
  }

  sprintf(KEY,"%s_FILTER", KEYCHECK);
  if ( keyMatchSim(MAXKEY, KEY, WORDS[0], keySource) ) {

    sscanf(WORDS[1] , "%s", cfilt ); // e.g., griz or g,r,i,z
    sscanf(WORDS[2] , "%s", cval  ); // value or comma-sep list of values

    NFILT = PARSE_FILTLIST(cfilt, ifilt_list );  // return ifilt_list

    if ( strstr(cfilt,COMMA) != NULL ) {
      // parse comman sep list of float values
      parse_commaSepList(fnam, cval, MXFILTINDX, 20, 
			 &NVAL, &cval_list); // <== returned
      if ( NVAL != NFILT ) {
	sprintf(c1err,"NFILT=%d does not match NVAL=%d",
		NFILT, NVAL);
	sprintf(c2err,"Check input %s %s %s", KEY, cfilt, cval);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
      } 
      for(ifilt=0; ifilt < NFILT; ifilt++) {
	sscanf(cval_list[ifilt], "%f", &val_list[ifilt] );
	free(cval_list[ifilt]);
      }
      free(cval_list);
    }
    else {
      // parse single float val for all filters
      sscanf(cval , "%f", &ftmp );
      for(ifilt=0; ifilt<NFILT; ifilt++) { val_list[ifilt] = ftmp; }
    }

    // load output values 
    for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
      ifilt_obs = ifilt_list[ifilt] ;
      VALUE_FILTERLIST[ifilt_obs] = val_list[ifilt];
    }      
    return(2);
  }

  return(0);

} // end parse_input_KEY_PLUS_FILTER


// ======================================
int parse_input_SIMLIB(char **WORDS, int keySource ) {

  // Created July 2020
  // parse keys containing SIMLIB_
  int ITMP, N=0;
  char fnam[] = "parse_input_SIMLIB" ;

  // ---------------- BEGIN ------------

  if ( keyMatchSim(1, "SIMLIB_FILE",  WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS.SIMLIB_FILE );
  }

  else if ( keyMatchSim(1, "SIMLIB_SURVEY",  WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS.SIMLIB_SURVEY );
  }

  else if ( keyMatchSim(1, "SIMLIB_FIELDLIST",  WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], 200 );
    N++;  sscanf(WORDS[N], "%s", INPUTS.SIMLIB_FIELDLIST );

    char *FIELDLIST       = INPUTS.SIMLIB_FIELDLIST;
    STRING_DICT_DEF *DICT = &INPUTS.DICT_SIMLIB_FIELDLIST_PRESCALE ;
    parse_string_prescales(FIELDLIST, DICT); // store pre-scales in *DICT
    INPUTS.SIMLIB_FIELDSKIP_FLAG = 1 ; // count skipped fields in NGENTOT
  }
  else if ( keyMatchSim(1, "SIMLIB_IDSTART",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_IDSTART );
  }
  else if ( keyMatchSim(1, "SIMLIB_IDLOCK",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_IDLOCK );
  }
  else if ( keyMatchSim(1, "SIMLIB_MAXRANSTART",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_MAXRANSTART );
  }
  else if ( keyMatchSim(1, "SIMLIB_MINOBS",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_MINOBS );
  }
  else if ( keyMatchSim(1, "SIMLIB_MINSEASON",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_MINSEASON );
  }
  else if ( keyMatchSim(1, "SIMLIB_DUMP",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_DUMP );
  }
  else if ( keyMatchSim(1, "SIMLIB_NREPEAT",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_NREPEAT );
  }
  else if ( keyMatchSim(1, "SIMLIB_NSKIPMJD",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_NSKIPMJD );
  }
  else if ( keyMatchSim(1, "SIMLIB_IDSKIP",  WORDS[0],keySource) ) {
    int NSKIP = INPUTS.NSKIP_SIMLIB;
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_IDSKIP[NSKIP] );
    INPUTS.NSKIP_SIMLIB++ ;
  }
  else if ( keyMatchSim(1, "USE_SIMLIB_REDSHIFT",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.USE_SIMLIB_REDSHIFT );
    INPUTS.USE_SIMLIB_GENOPT=1;
  }
  else if ( keyMatchSim(1, "USE_SIMLIB_DISTANCE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.USE_SIMLIB_DISTANCE );
    INPUTS.USE_SIMLIB_GENOPT=1;
  }
  else if ( keyMatchSim(1, "USE_SIMLIB_PEAKMJD",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.USE_SIMLIB_PEAKMJD );
    INPUTS.USE_SIMLIB_GENOPT=1;
  }
  else if ( keyMatchSim(1, "USE_SIMLIB_MAGOBS",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.USE_SIMLIB_MAGOBS );
    INPUTS.USE_SIMLIB_GENOPT=1;
  }
  else if ( keyMatchSim(1, "USE_SIMLIB_SPECTRA",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.USE_SIMLIB_SPECTRA );
    INPUTS.USE_SIMLIB_GENOPT=1;
  }
  else if ( keyMatchSim(1, "USE_SIMLIB_SALT2",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.USE_SIMLIB_SALT2 );
    INPUTS.USE_SIMLIB_GENOPT=1;
  }
  else if ( keyMatchSim(1, "SIMLIB_MSKOPT",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.SIMLIB_MSKOPT );
  }

  return(N);
} // end parse_input_SIMLIB

// ======================================
int parse_input_HOSTLIB(char **WORDS, int keySource ) {

  // Created July 2020
  // parse keys starting with HOSTLIB
  // Oct 16 2020: check IGNOREFILE(HOSTLIB_FILE)
  // Dec 02 2020: fix bug setting MSKOPT to allow command-line override.
  // May 04 2021: restore +HOSTMAGS and +HOSTNBR 
  // Jun 02 2021: add calls to check_arg_len
  // Jun 17 2021: restore SEPNBR_MAX & NNBR_WRITE_MAX
  // Jul 01 2021: read forgotten INPUTS.HOSTLIB_MAXDDLR
  // Jul 16 2021: fix override logic by setting INPUTS.HOSTLIB_USE=1
  //               only if it is not already set.

  int  j, ITMP, N=0, nread ;
  char fnam[] = "parse_input_HOSTLIB" ;

  // ------------ BEGIN ------------

  if ( keyMatchSim(1, "HOSTLIB_FILE", WORDS[0], keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS.HOSTLIB_FILE ) ; 
    if ( IGNOREFILE(INPUTS.HOSTLIB_FILE) ) 
      { INPUTS.HOSTLIB_MSKOPT = INPUTS.HOSTLIB_USE = 0;  }
    else { 
      setbit_HOSTLIB_MSKOPT(HOSTLIB_MSKOPT_USE); 
      if ( !INPUTS.HOSTLIB_USE )  { INPUTS.HOSTLIB_USE = 1; }
    }
  }
  else if ( keyMatchSim(1, "HOSTLIB_WGTMAP_FILE", WORDS[0], keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS.HOSTLIB_WGTMAP_FILE ) ; 
  }
  else if ( keyMatchSim(1, "HOSTLIB_ZPHOTEFF_FILE", WORDS[0], keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS.HOSTLIB_ZPHOTEFF_FILE ) ; 
  }

  else if ( keyMatchSim(1, "HOSTLIB_SPECBASIS_FILE", WORDS[0], keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS.HOSTLIB_SPECBASIS_FILE ) ; 
  }
  else if ( keyMatchSim(1, "HOSTLIB_SPECDATA_FILE", WORDS[0], keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN );
    N++;  sscanf(WORDS[N], "%s", INPUTS.HOSTLIB_SPECDATA_FILE ) ; 
  }
  else if ( keyMatchSim(1, "HOSTLIB_MSKOPT", WORDS[0], keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.HOSTLIB_MSKOPT );
    setbit_HOSTLIB_MSKOPT(HOSTLIB_MSKOPT_USE) ;
  }

  else if ( keyMatchSim(1, "+HOSTMAGS", WORDS[0], keySource ) ) {
    INPUTS.HOSTLIB_MSKOPT += HOSTLIB_MSKOPT_PLUSMAGS ;
    N += FLAG_NWD_ZERO; // flag that key has no argument 
    setbit_HOSTLIB_MSKOPT(HOSTLIB_MSKOPT_USE) ;
    INPUTS.HOSTLIB_USE = 2; // set rewrite flag
    sprintf(INPUTS.HOSTLIB_PLUS_COMMAND,"%s", WORDS[0]);
  }

  else if ( keyMatchSim( 1, "+HOSTNBR", WORDS[0], keySource ) ) {
    INPUTS.HOSTLIB_MSKOPT += HOSTLIB_MSKOPT_PLUSNBR ;
    N += FLAG_NWD_ZERO; // flag that key has no argument
    setbit_HOSTLIB_MSKOPT(HOSTLIB_MSKOPT_USE) ;
    INPUTS.HOSTLIB_USE = 2; // set rewrite flag
    sprintf(INPUTS.HOSTLIB_PLUS_COMMAND,"%s", WORDS[0]);
  }
  else if ( keyMatchSim( 1, "SEPNBR_MAX", WORDS[0], keySource ) ) {
    N++; sscanf(WORDS[N], "%le", &HOSTLIB_NBR_WRITE.SEPNBR_MAX );
  }
  else if ( keyMatchSim( 1, "NNBR_WRITE_MAX", WORDS[0], keySource ) ) {
    N++; sscanf(WORDS[N], "%d", &HOSTLIB_NBR_WRITE.NNBR_WRITE_MAX );
  }

  else if ( keyMatchSim( 1, "+HOSTAPPEND", WORDS[0], keySource ) ) {
    INPUTS.HOSTLIB_MSKOPT += HOSTLIB_MSKOPT_APPEND ;
    N++ ; sscanf(WORDS[N], "%s", INPUTS.HOSTLIB_APPEND_FILE );
    setbit_HOSTLIB_MSKOPT(HOSTLIB_MSKOPT_USE) ;
    INPUTS.HOSTLIB_USE = 2; // set rewrite flag
    sprintf(INPUTS.HOSTLIB_PLUS_COMMAND,"%s", WORDS[0]);
  }

  else if ( keyMatchSim(1,"HOSTLIB_GENZPHOT_FUDGEPAR",WORDS[0],keySource)){
    // read first 4 elements as float
    for(j=0; j < 4; j++ )  {
      N++; nread = sscanf(WORDS[N],"%f",&INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[j]);
      if ( nread != 1 ) { abort_bad_input(WORDS[0], WORDS[N], j, fnam); }
    }
    // read 5th element as string
    N++; parse_input_GENZPHOT_OUTLIER(WORDS[N]);
  }
  else if ( keyMatchSim(1, "HOSTLIB_GENZPHOT_BIAS",WORDS[0],keySource) ) {
    for(j=0; j < 4; j++ )  { 
      N++; nread = sscanf(WORDS[N],"%f",&INPUTS.HOSTLIB_GENZPHOT_BIAS[j]) ; 
      if ( nread != 1 ) { abort_bad_input(WORDS[0], WORDS[N], j, fnam); }
    }
  }
  else if ( keyMatchSim(1, "HOSTLIB_DZTOL",WORDS[0],keySource) ) {
    N += read_genpoly(WORDS[0], &WORDS[1], 2, &INPUTS.HOSTLIB_GENPOLY_DZTOL);
  }
  else if ( keyMatchSim(1, "HOSTLIB_SCALE_SERSIC_SIZE  HOSTLIB_SCALE_SERSIC",
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_SCALE_SERSIC_SIZE ) ;   
  }
  else if ( keyMatchSim(1, "HOSTLIB_SCALE_LOGMASS_ERR",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_SCALE_LOGMASS_ERR ) ;   
  }
  else if ( keyMatchSim(1, "HOSTLIB_STOREVAR  HOSTLIB_STOREPAR",
			WORDS[0],keySource) ) {
    check_arg_len(WORDS[0], WORDS[1], MXPATHLEN);
    N++;  sscanf(WORDS[N], "%s", INPUTS.HOSTLIB_STOREPAR_LIST ) ; 
  }  
  else if ( keyMatchSim(1, "HOSTLIB_MAXREAD",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.HOSTLIB_MAXREAD ) ;
  }
  else if ( keyMatchSim(1, "HOSTLIB_GALID_NULL",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.HOSTLIB_GALID_NULL ) ;
  }
  else if ( keyMatchSim(1, "HOSTLIB_GALID_PRIORITY",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.HOSTLIB_GALID_PRIORITY[0] ) ;
    N++;  sscanf(WORDS[N], "%d", &INPUTS.HOSTLIB_GALID_PRIORITY[1] ) ;
  }
  else if ( keyMatchSim(1, "HOSTLIB_GALID_UNIQUE", 
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.HOSTLIB_GALID_UNIQUE ) ;
  }
  else if ( keyMatchSim(1, "HOSTLIB_MINDAYSEP_SAMEGAL",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL);
  }
  else if ( keyMatchSim(1, "HOSTLIB_MNINTFLUX_SNPOS",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.HOSTLIB_MNINTFLUX_SNPOS);
  }
  else if ( keyMatchSim(1, "HOSTLIB_MXINTFLUX_SNPOS",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.HOSTLIB_MXINTFLUX_SNPOS);
  }
  else if ( keyMatchSim(1, "HOSTLIB_GENRANGE_NSIGZ",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.HOSTLIB_GENRANGE_NSIGZ[0] );
    N++;  sscanf(WORDS[N], "%f", &INPUTS.HOSTLIB_GENRANGE_NSIGZ[1] );
  }
  else if ( keyMatchSim(1, "HOSTLIB_GENRANGE_RA",WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_GENRANGE_RA[0] );
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_GENRANGE_RA[1] );
  }
  else if ( keyMatchSim(1, "HOSTLIB_GENRANGE_DEC HOSTLIB_GENRANGE_DECL", 
			WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_GENRANGE_DEC[0] );
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_GENRANGE_DEC[1] );
  }
  else if ( keyMatchSim(1, "HOSTLIB_MAXDDLR", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.HOSTLIB_MAXDDLR );
  }
  else if ( keyMatchSim(1, "HOSTLIB_SBRADIUS", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_SBRADIUS );
  }
  else if ( keyMatchSim(1, "HOSTLIB_GALID_FORCE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.HOSTLIB_GALID_FORCE );
  }
  else if ( keyMatchSim(1, "HOSTLIB_ABMAG_FORCE", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_ABMAG_FORCE );
  }
  else if ( keyMatchSim(1, "HOSTLIB_FIXRAN_RADIUS", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_FIXRAN_RADIUS );
  }
  else if ( keyMatchSim(1, "HOSTLIB_FIXRAN_PHI", WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_FIXRAN_PHI );
  }
  else if ( keyMatchSim(1, "HOSTLIB_FIXSERSIC", WORDS[0],keySource) ) {
    for(j=0; j < 4; j++ ) {
      N++; nread = sscanf(WORDS[N], "%le", &INPUTS.HOSTLIB_FIXSERSIC[0] ); 
      if ( nread != 1 ) { abort_bad_input(WORDS[0], WORDS[N], j, fnam); }
    }
  }


  return(N);

} // end parse_input_HOSTLIB

// ================================================
int  parse_input_LCLIB(char **WORDS, int keySource ) {

  // Created July 2020
  // Feb 11 2021: Fix address bug reading LCLIB_CUTWIN

  int NCUT, N=0;
  char fnam[] = "parse_input_LCLIB";

  // ---------- BEGIN ----------

  if ( keyMatchSim(MXCUT_LCLIB, "LCLIB_CUTWIN",  WORDS[0],keySource) ) {
    // --- LCLIB cut windows to select library events  
    NCUT = LCLIB_CUTS.NCUTWIN ;
    N++ ; sscanf(WORDS[N], "%s",  LCLIB_CUTS.PARNAME[NCUT] );
    N++ ; sscanf(WORDS[N], "%le", &LCLIB_CUTS.CUTWIN[NCUT][0] );
    N++ ; sscanf(WORDS[N], "%le", &LCLIB_CUTS.CUTWIN[NCUT][1] );
    LCLIB_CUTS.NCUTWIN++ ;
  }
  else if ( keyMatchSim(1, "LCLIB_DEBUG_DAYSCALE",  WORDS[0],keySource) ) {
    N++;  sscanf(WORDS[N], "%le", &LCLIB_DEBUG.DAYSCALE );
  }
  else if ( keyMatchSim(1, "LCLIB_DEBUG_TOBS_OFFSET", WORDS[0],keySource)) {
    N++;  sscanf(WORDS[N], "%le", &LCLIB_DEBUG.TOBS_OFFSET_RANGE[0] );
    N++;  sscanf(WORDS[N], "%le", &LCLIB_DEBUG.TOBS_OFFSET_RANGE[1] );
  }
  else if ( keyMatchSim(1, "LCLIB_DEBUG_ZERO_TEMPLATE_FLUX", WORDS[0],keySource)) {
    N++;  sscanf(WORDS[N], "%d", &LCLIB_DEBUG.ZERO_TEMPLATE_FLUX );
  }
  else if ( keyMatchSim(1, "LCLIB_DEBUG_FORCE_NREPEAT", WORDS[0],keySource)) {
    N++;  sscanf(WORDS[N], "%d", &LCLIB_DEBUG.FORCE_NREPEAT);
  }

  return(N);

} // end parse_input_LCLIB


// ================================================
int  parse_input_CUTWIN(char **WORDS, int keySource ) {

  int  NCUT, N=0;
  char fnam[] = "parse_input_CUTWIN" ;

  // ---------- BEGIN ----------

 if ( keyMatchSim(1, "APPLY_CUTWIN_OPT", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%d", &INPUTS.APPLY_CUTWIN_OPT );
  }
 else if ( keyMatchSim(1, "EPCUTWIN_LAMREST", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.EPCUTWIN_LAMREST[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.EPCUTWIN_LAMREST[1] );
 } 
 else if ( keyMatchSim(1, "EPCUTWIN_SNRMIN", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.EPCUTWIN_SNRMIN[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.EPCUTWIN_SNRMIN[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_REDSHIFT CUTWIN_REDSHIFT_TRUE", 
		       WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_REDSHIFT_TRUE[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_REDSHIFT_TRUE[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_REDSHIFT_FINAL", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_REDSHIFT_FINAL[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_REDSHIFT_FINAL[1] );
 } 
 else if ( keyMatchSim(1, "CUTWIN_HOST_PHOTOZ CUTWIN_HOST_ZPHOT",
		       WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_HOST_ZPHOT[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_HOST_ZPHOT[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_TRESTMIN", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_TRESTMIN[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_TRESTMIN[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_TRESTMAX", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_TRESTMAX[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_TRESTMAX[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_TGAPMAX", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_TGAPMAX[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_TGAPMAX[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_T0GAPMAX", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_T0GAPMAX[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_T0GAPMAX[1] );
 }
 else if ( keyMatchSim(MXCUTWIN_SNRMAX, "CUTWIN_SNRMAX", WORDS[0],keySource)) {
   INPUTS.NCUTWIN_SNRMAX++ ;  // should refactor later to start at zero, not 1
   NCUT = INPUTS.NCUTWIN_SNRMAX ;
   N++ ; sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_SNRMAX[NCUT][0]       );
   N++ ; sscanf(WORDS[N], "%s", INPUTS.CUTWIN_SNRMAX_FILTERS[NCUT]   );
   N++ ; sscanf(WORDS[N], "%s", INPUTS.CUTWIN_SNRMAX_LOGIC[NCUT]     );
   N++ ; sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_SNRMAX_TREST[NCUT][0] );
   N++ ; sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_SNRMAX_TREST[NCUT][1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_NEPOCH", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_NEPOCH[0] );   // NOBS cut
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_NEPOCH[1] );   // SNRMIN cut
 }
 else if ( keyMatchSim(1, "CUTWIN_NOBSDIF", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%d", &INPUTS.CUTWIN_NOBSDIF[0] ); 
   N++;  sscanf(WORDS[N], "%d", &INPUTS.CUTWIN_NOBSDIF[1] ); 
 }
 else if ( keyMatchSim(1, "CUTWIN_MJDDIF", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_MJDDIF[0] ); 
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_MJDDIF[1] ); 
 }
 else if ( keyMatchSim(1, "CUTWIN_NOBS_SATURATE", WORDS[0],keySource)) {
   NCUT = INPUTS.NCUTWIN_SATURATE;
   N++ ; sscanf(WORDS[N], "%d", &INPUTS.CUTWIN_SATURATE_NOBS[NCUT][0] );
   N++ ; sscanf(WORDS[N], "%d", &INPUTS.CUTWIN_SATURATE_NOBS[NCUT][1] );
   N++ ; sscanf(WORDS[N], "%s",  INPUTS.CUTWIN_SATURATE_FILTERS[NCUT] );
   INPUTS.NCUTWIN_SATURATE++ ;
 }
 else if ( keyMatchSim(1, "CUTWIN_NOBS_NOSATURATE", WORDS[0],keySource)) {
   NCUT = INPUTS.NCUTWIN_NOSATURATE;
   N++ ; sscanf(WORDS[N], "%d", &INPUTS.CUTWIN_NOSATURATE_NOBS[NCUT][0] );
   N++ ; sscanf(WORDS[N], "%d", &INPUTS.CUTWIN_NOSATURATE_NOBS[NCUT][1] );
   N++ ; sscanf(WORDS[N], "%s",  INPUTS.CUTWIN_NOSATURATE_FILTERS[NCUT] );
   INPUTS.NCUTWIN_NOSATURATE++ ;
 }
 else if ( keyMatchSim(1, "CUTWIN_MWEBV", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_MWEBV[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_MWEBV[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_PEAKMAG", WORDS[0],keySource)) {
   // PEAKMAG cut on brightest epoch only
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_PEAKMAG[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_PEAKMAG[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_PEAKMAG_ALL", WORDS[0],keySource)) {
   // require ALL filters to satisfy PEAKMAG cut
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_PEAKMAG_ALL[0] );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_PEAKMAG_ALL[1] );
 }
 else if ( keyMatchSim(1, "CUTWIN_PEAKMAG_BYFIELD", WORDS[0],keySource)) {
   INPUTS.NCUTWIN_PEAKMAG_BYFIELD++ ;  NCUT = INPUTS.NCUTWIN_PEAKMAG_BYFIELD;
   N++;  sscanf(WORDS[N], "%f", INPUTS.CUTWIN_PEAKMAG_BYFIELD[NCUT][0] );
   N++;  sscanf(WORDS[N], "%f", INPUTS.CUTWIN_PEAKMAG_BYFIELD[NCUT][1] );
   N++;  sscanf(WORDS[N], "%s", INPUTS.CUTWIN_BYFIELDLIST[NCUT] );
 }
 else if ( keyMatchSim(1, "CUTWIN_EPOCHS_SNRMIN", WORDS[0],keySource)) {
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_EPOCHS_SNRMIN );
   N++;  sscanf(WORDS[N], "%f", &INPUTS.CUTWIN_EPOCHS_TRANGE[1] );
   N++;  sscanf(WORDS[N], "%s", INPUTS.CUTWIN_EPOCHS_FILTERS );
 }

 return(N);

} // end parse_input_CUTWIN

// ==================================================
#ifdef MODELGRID_GEN
int parse_input_GRIDGEN(char **WORDS, int keySource) {

  // Created July 2020
  // read/parse input keys to generate grid of population params.

  int N=0, *IPTR ;
  char fnan[] = "parse_input_GRIDGEN";

  // ------------- BEGIN -------------

  // check if source is GRID
  if ( strcmp(INPUTS.GENSOURCE,"GRID") != 0 )  { return(N); }

  if ( keyMatchSim(1,"GRID_FORMAT", WORDS[0], keySource) ) {
    N++ ; sscanf(WORDS[N], "%s", GRIDGEN_INPUTS.FORMAT );
  }
  else if ( keyMatchSim(1,"NGRID_LOGZ", WORDS[0], keySource) ) {
    IPTR = &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_LOGZ] ;
    N++ ; sscanf(WORDS[N], "%d", IPTR);
  }
  else if ( keyMatchSim(1,"NGRID_SHAPEPAR NGRID_LUMIPAR", WORDS[0], keySource) ) {
    IPTR = &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_SHAPEPAR] ;
    N++ ; sscanf(WORDS[N], "%d", IPTR);
  }
  else if ( keyMatchSim(1,"NGRID_COLORPAR", WORDS[0], keySource) ) {
    IPTR = &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_COLORPAR] ;
    N++ ; sscanf(WORDS[N], "%d", IPTR);
  }
  else if ( keyMatchSim(1,"NGRID_COLORLAW", WORDS[0], keySource) ) {
    IPTR = &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_COLORLAW] ;
    N++ ; sscanf(WORDS[N], "%d", IPTR);
  }
  else if ( keyMatchSim(1,"NGRID_TREST", WORDS[0], keySource) ) {
    IPTR = &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_TREST] ;
    N++ ; sscanf(WORDS[N], "%d", IPTR);
  }

  return(N);

}  // end parse_input_GRIDGEN
#endif


// ==============================================================
int parse_input_SOLID_ANGLE(char **WORDS, int keySource) {

  // Created July 2020 [refactor]
  //
  // Parse sim-input for
  // SOLID_ANGLE:  <value>
  //    or
  // SOLID_ANGLE(FIELDLIST):  <value>
  //
  // The latter sets FIELDSKIP_FLAG=0 so that skipped fields
  // are NOT counted in NGENTOT --> reads SIMLIB as if the
  // skipped fields were never there and assumes SOLID_ANGLE
  // corresponds to FIELDLIST subset.

  int N=0;
  char *FIELDLIST = INPUTS.SIMLIB_FIELDLIST ;
  char KEYLOCAL[100], FIELDLIST_ADD[100];
  int  ADDFIELD = 0 ;
  int  LDMP   = 0 ; 
  char fnam[] = "parse_input_SOLID_ANGLE" ;


  // ------------ BEGIN ---------------

  sprintf(KEYLOCAL, "%s", WORDS[0]);

  // extract contents of optional ()
  extractStringOpt(KEYLOCAL,FIELDLIST_ADD);

  if ( LDMP ) {
    printf(" xxx %s: input key = '%s' \n", fnam, WORDS[0] );
    printf(" xxx %s: KEYLOCAL  = '%s' \n", fnam, KEYLOCAL );
    printf(" xxx %s: FIELDLIST_ADD  = '%s' \n", fnam, FIELDLIST_ADD );
    fflush(stdout);
  }

  if ( keyMatchSim(1, "SOLID_ANGLE", KEYLOCAL, keySource) ) {
    N++ ; sscanf(WORDS[N] , "%f", &INPUTS.SOLID_ANGLE ); 
  }
  else {
    return(N);
  }
    
  if ( LDMP ) {
    printf(" xxx %s: SOLID_ANGLE  = %f \n", fnam, INPUTS.SOLID_ANGLE ); 
    fflush(stdout);
  }

  ADDFIELD = strlen(FIELDLIST_ADD);
  if ( ADDFIELD && strcmp(FIELDLIST,"ALL") != 0 ) {
    sprintf(c1err, "Cannot declare FIEDLIST in two places.");
    sprintf(c2err, "Check SOLID_ANGLE and SIMLIB_FIELDLIST keys.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }
  if ( ADDFIELD ) {
    sprintf(FIELDLIST,"%s", FIELDLIST_ADD);
    INPUTS.SIMLIB_FIELDSKIP_FLAG = 0 ;
  }

  /*
  printf(" xxx %s: KEY= '%s'  SKIP_FLAG=%d \n",
	 fnam, KEYLOCAL, INPUTS.SIMLIB_FIELDSKIP_FLAG );
  */

  return(N) ;

} // end parse_input_SOLID_ANGLE


// ==============================
int parse_input_ZVARIATION(char **WORDS, int keySource) {

  // Created Sep 30 2020
  // Parse ZVARIATION file, or polynomial ... latter using read_genpoly.

  int  N=0, NPAR ;
  char *parName, polyVarName[60] ;
  char fnam[] = "parse_input_ZVARIATION" ;

  // ------------ BEGIN -----------

  if ( keyMatchSim(1,"ZVARIATION_FILE", WORDS[0], keySource)  ) {  
    N++; sscanf(WORDS[N], "%s", INPUT_ZVARIATION_FILE ); 
    if ( !IGNOREFILE(INPUT_ZVARIATION_FILE) ) { USE_ZVAR_FILE = 1 ; }
  }

  if ( keyMatchSim(10, "ZVARIATION_POLY", WORDS[0], keySource) ) { 
      NPAR = NPAR_ZVAR_USR ;
      INPUT_ZVARIATION[NPAR].FLAG = FLAG_ZPOLY_ZVAR ;
      INPUT_ZVARIATION[NPAR].NZBIN = 0 ;
      NPAR_ZVAR_USR++ ;  
      
      parName = INPUT_ZVARIATION[NPAR].PARNAME ;
      N++; sscanf(WORDS[N], "%s", parName );

      // read either comma-sep poly of arbitrary order,
      // or read legacy space-sep poly with order = POLYORDER_ZVAR.
      N += read_genpoly(WORDS[0], &WORDS[2], POLYORDER_ZVAR, 
			&INPUT_ZVARIATION[NPAR].POLY);
    }

  return(N);

} // end parse_input_ZVARIATION

// ==============================
int parse_input_ZVARIATION_LEGACY(char **WORDS, int keySource) {

  // Sep 2020: LEGACY -> does not use read_genpoly(...)

  int N=0, j, NPAR ;
  char tmpLine[60], cpoly[60];
  char fnam[] = "parse_input_ZVARIATION_LEGACY" ;

  // ------------ BEGIN -----------

  if ( keyMatchSim(1,"ZVARIATION_FILE", WORDS[0], keySource)  ) {  
    N++; sscanf(WORDS[N], "%s", INPUT_ZVARIATION_FILE ); 
    if ( !IGNOREFILE(INPUT_ZVARIATION_FILE) ) { USE_ZVAR_FILE = 1 ; }
    return(N) ;
  }

  if ( keyMatchSim(10, "ZVARIATION_POLY", WORDS[0], keySource) ) { 
      NPAR = NPAR_ZVAR_USR ;
      INPUT_ZVARIATION[NPAR].FLAG = FLAG_ZPOLY_ZVAR ;
      INPUT_ZVARIATION[NPAR].NZBIN = 0 ;
      NPAR_ZVAR_USR++ ;  

      N++; sscanf(WORDS[N], "%s", INPUT_ZVARIATION[NPAR].PARNAME ) ;
      N++; sscanf(WORDS[N], "%s", tmpLine); 
      // if there is a comma, use new GENPOLY struct; else read legacy format
      if ( strstr(tmpLine,COMMA) ) {
	// store comma-separate poly coefficients in cpoly
	sprintf(cpoly, "%s", tmpLine);	
      }
      else {
	// read LEGACY format of space-separated 3rd order poly
	double zpoly[POLYORDER_ZVAR+2];
	sscanf(tmpLine, "%le", &zpoly[0]); 
	for(j=1; j <= POLYORDER_ZVAR; j++ ) 
	  { N++; sscanf(WORDS[N], "%le", &zpoly[j] ); }

	sprintf(cpoly,"%f,%f,%f,%f", zpoly[0],zpoly[1],zpoly[2],zpoly[3]);
      }
      parse_GENPOLY(cpoly, "z", &INPUT_ZVARIATION[NPAR].POLY, fnam);
      return(N);
    }

  return(N);

} // end parse_input_ZVARIATION_LEGACY

// ======================================
int parse_input_RANSYSTPAR(char **WORDS, int keySource ) {

  // Jul 20 2020 [refactor]
  // Parse RANSYSTPAR_XXX keys.
  // Function returns number of WORDS that are read.
  // 
  // Oct 30 2020: 
  //   + check for filter list in parentheses for SIGZP and SIGLAMFILT.
  //

  int  NFILTDEF = INPUTS.NFILTDEF_OBS ;  
  int  N = 0 ;
  int  i, igrp, ifilt, ifilt_obs, NFILTGROUP, NFILT_PER_GROUP[MXFILTINDX] ;
  double tmpVal;
  float *ptrShift ;
  char KEYNAME[60], FILTGROUP_STRING[MXFILTINDX+10], **FILTGROUP_LIST ;
  char cfilt[2] ;
  char KEYNAMES_FILTER[] = 
    "RANSYSTPAR_SIGSHIFT_ZP  RANSYSTPAR_SIGSHIFT_LAMFILT" ;
  char fnam[] = "parse_input_RANSYSTPAR" ;

  // ---------- BEGIN ----------

  // if WORDS[0] of of the for RANSYSTPAR_XYZ(BLA), then return
  // KEYNAME = "RANSYSTPAR_XYZ" and SYST_OPT = "BLA"
  sprintf(KEYNAME, "%s", WORDS[0]);   
  extractStringOpt(KEYNAME,FILTGROUP_STRING) ;

  // - - - - - - 

  if ( keyMatchSim(1, "RANSYSTPAR_RANSEED_GEN", 
		   KEYNAME, keySource) ) {
    N++;  sscanf(WORDS[N], "%d", &INPUTS.RANSYSTPAR.RANSEED_GEN );
  }
  else if ( keyMatchSim(1, "RANSYSTPAR_SIGSCALE_FLUXERR", 
		   KEYNAME, keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.SIGSCALE_FLUXERR );
  }
  else if ( keyMatchSim(1,"RANSYSTPAR_SIGSCALE_FLUXERR2", 
			KEYNAME, keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.SIGSCALE_FLUXERR2 );
  }
  else if ( keyMatchSim(1,"RANSYSTPAR_SIGSCALE_MWEBV", 
			KEYNAME, keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.SIGSCALE_MWEBV );
  }
  else if ( keyMatchSim(1,"RANSYSTPAR_SIGSHIFT_MWRV", 
			KEYNAME, keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.SIGSHIFT_MWRV );
  }

  else if ( keyMatchSim(1,"RANSYSTPAR_SIGSHIFT_REDSHIFT", // P.Armstrong, 2020
			KEYNAME, keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.SIGSHIFT_REDSHIFT );
  }

  else if ( keyMatchSim(1,"RANSYSTPAR_GENMODEL_WILDCARD",  // P.Armstron, 2020
			KEYNAME, keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.RANSYSTPAR.GENMODEL_WILDCARD ); 
  }
  else if ( keyMatchSim(1,"RANSYSTPAR_GENPDF_FILE_WILDCARD", 
			KEYNAME, keySource) ) {
    N++;  sscanf(WORDS[N], "%s", INPUTS.RANSYSTPAR.GENPDF_FILE_WILDCARD ); 
  }

  else if ( keyMatchSim(1,"RANSYSTPAR_SIGSHIFT_OMEGA_MATTER",
			KEYNAME,keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.SIGSHIFT_OMEGA_MATTER );
  }

  else if ( keyMatchSim(1,"RANSYSTPAR_RANGESHIFT_OMEGA_MATTER",
			KEYNAME,keySource) ) {
    N++;  sscanf(WORDS[N],"%f",&INPUTS.RANSYSTPAR.RANGESHIFT_OMEGA_MATTER[0]);
    N++;  sscanf(WORDS[N],"%f",&INPUTS.RANSYSTPAR.RANGESHIFT_OMEGA_MATTER[1]);
  }
  else if ( keyMatchSim(1,"RANSYSTPAR_SIGSHIFT_W0",
			KEYNAME,keySource) ) {
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.SIGSHIFT_W0 );
  }
  else if ( keyMatchSim(1,"RANSYSTPAR_RANGESHIFT_W0",
			KEYNAME,keySource)){
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.RANGESHIFT_W0[0] );
    N++;  sscanf(WORDS[N], "%f", &INPUTS.RANSYSTPAR.RANGESHIFT_W0[1] );
  }
  else if ( keyMatchSim(1, KEYNAMES_FILTER, KEYNAME,keySource) ) {

    if ( NFILTDEF == 0 ) {
      sprintf(c1err,"Must define GENFILTERS before %s", KEYNAME );
      sprintf(c2err,"Check key order in sim-input file");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    
    if ( strstr(KEYNAME,"ZP") != NULL ) 
      { ptrShift = INPUTS.RANSYSTPAR.SIGSHIFT_ZP ; }
    else if ( strstr(KEYNAME,"LAMFILT") != NULL ) 
      { ptrShift = INPUTS.RANSYSTPAR.SIGSHIFT_LAMFILT ; }
    else 
      { return 0 ; }


    // if no filter argument, then assume all filters are specified.
    if ( strlen(FILTGROUP_STRING) == 0 ) { 
      for(ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
	sprintf(cfilt, "%c", INPUTS.GENFILTERS[ifilt] ); 
	catVarList_with_comma(FILTGROUP_STRING,cfilt); 
      }
    }

    parse_commaSepList("RANSYSTPAR-OPT", FILTGROUP_STRING, 
		       MXFILTINDX,MXFILTINDX, &NFILTGROUP, &FILTGROUP_LIST);

    for(ifilt=0; ifilt < NFILTGROUP; ifilt++ ) 
      { NFILT_PER_GROUP[ifilt] = strlen(FILTGROUP_LIST[ifilt]); }
 
    for(igrp=0; igrp < NFILTGROUP; igrp++ ) {
      tmpVal = NOFLOAT ;       
      N++; sscanf(WORDS[N], "%le", &tmpVal ) ;
      // printf(" xxx %s: igrp=%d/%d   N=%d tmpVal=%f\n",
      //    fnam, igrp, NFILTGROUP,   N, tmpVal);

      if ( tmpVal == NOFLOAT ) {
	sprintf(c1err,"Expecting to read float, but read '%s' ", WORDS[N] ) ;
	sprintf(c2err,"Check igrp=%d after %s", igrp, KEYNAME );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
      }
      
      for(i=0; i < NFILT_PER_GROUP[igrp]; i++ ) {
	sprintf(cfilt, "%c", FILTGROUP_LIST[igrp][i] );
	ifilt = index_charString(cfilt, INPUTS.GENFILTERS);
	ptrShift[ifilt] = tmpVal ;
      }
  
    }    // end igrp loop 
    //  debugexit(fnam) ;
  }

  if ( N > 0 ) {  INPUTS.RANSYSTPAR.USE = 1; }

  return(N) ;

} // end  parse_input_RANSYSTPAR

// ==============================================================
int parse_input_GENMODEL(char **WORDS, int keySource) {

  // Created July 2020 [refactor]
  // Parse GENMODEL key.

  int  N=0;
  int  jnam;
  char *GENMODEL = INPUTS.GENMODEL ;
  char ctmp[60], *NAME0 ;
  bool LDMP   = false ;
  char fnam[] = "parse_input_GENMODEL" ;
  
  // ---------- BEGIN ------------


  N++ ; sscanf(WORDS[N], "%s", GENMODEL ); 

  // check for path + model 
  extract_MODELNAME(GENMODEL,                   // input path/model
		    INPUTS.MODELPATH, INPUTS.MODELNAME); // returned

  // check for NONIA 
  INPUTS.NON1A_MODELFLAG = get_NON1A_MODELFLAG(INPUTS.MODELNAME);
  if ( INPUTS.NON1A_MODELFLAG == MODEL_NON1AGRID ) 
    { N++ ; sscanf(WORDS[N] , "%s", INPUTS.NON1AGRID_FILE );  }


  // check LCLIB (Galactic transients)
  NAME0 = GENMODEL_NAME[MODEL_LCLIB][0];
  if ( strcmp(GENMODEL,NAME0)==0 ) {
      N++ ; sscanf(WORDS[N] , "%s", INPUTS.LCLIB_FILE ); 
      N++ ; sscanf(WORDS[N] , "%s", INPUTS.LCLIB_TEMPLATE_EPOCHS ); 
  }

  // - - - - - - - - - - - - - - - - - - -
  // check python models: BYOSED & SNEMO
  NAME0 = GENMODEL_NAME[MODEL_BYOSED][0];
  if ( strcmp(GENMODEL,NAME0)==0 ) {
    N++ ; sscanf(WORDS[N], "%s", INPUTS.MODELPATH );
    IS_PySEDMODEL = true;
  }

  NAME0 = GENMODEL_NAME[MODEL_SNEMO][0];  // Sep 2020
  if ( strcmp(GENMODEL,NAME0)==0 ) {
    N++ ; sscanf(WORDS[N], "%s", INPUTS.MODELPATH ); 
    IS_PySEDMODEL = true;
  }

  // - - - - - - - - - - - - - - - - - - -
  // check fixmag/FIXMAG
  for(jnam=0; jnam < MXNAME_PER_MODEL; jnam++ ) {
    NAME0  = GENMODEL_NAME[MODEL_FIXMAG][jnam];
    if ( strcmp(GENMODEL,NAME0) == 0 ) {
      N++ ; sscanf(WORDS[N], "%s", ctmp ); 
      parse_input_FIXMAG(ctmp);	

      // check if rest-frame or obs-frame
      if ( strstr(GENMODEL,"mag") != NULL ) 
	{ INPUTS.GENFRAME_FIXMAG = GENFRAME_OBS; }
      else if ( strstr(GENMODEL,"MAG") != NULL ) 
	{ INPUTS.GENFRAME_FIXMAG = GENFRAME_REST; }
    }
  }

  // allow mlcs in place of mlcs2k2
  if ( strcmp(GENMODEL,"mlcs")== 0) { sprintf(GENMODEL, "mlcs2k2"); }    

  // - - - - - - - - - - - - - - - - 
  if ( LDMP ) {
    printf("\n xxx ----------------------------------------- \n");
    printf(" xxx %s DUMP: \n", fnam );
    printf(" xxx Input GENMODEL:  '%s' \n", GENMODEL);
    printf(" xxx INPUTS.MODELPATH = '%s' \n", INPUTS.MODELPATH );
    printf(" xxx INPUTS.MODELNAME = '%s' \n", INPUTS.MODELNAME );
    debugexit(fnam);
  }


  return(N) ;

} // end parse_input_GENMODEL


// ==============================================================
int parse_input_GENMODEL_ARGLIST(char **WORDS, int keySource) {

  // Created July 2020 [refactor]
  // parse GENMODEL_ARGLIST 'XX YY ZZ'
  // which can have arbitrary number of elements between quotes.
  // Initial use is to pass command-line overrides to BYOSED model.

  int  N=0;
  char fnam[] = "parse_input_GENMODEL_ARGLIST" ;
  
  // ---------- BEGIN ------------

  if ( keySource == KEYSOURCE_FILE ) {
    sprintf(c1err,"GENMODEL_ARGLIST valid only as command-line arg.");
    sprintf(c2err,"snlc_sim.exe <inFile> "
	    "GENMODEL_ARGLIST 'KEY1 VAL1 KEY2 VAL2 ...' ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // read command-line arg
  // (C code scoops everything within quotes)
  N++ ;  sprintf(INPUTS.GENMODEL_ARGLIST, "%s", WORDS[N] );

  
  return(N) ;

} // end parse_input_GENMODEL_ARGLIST


// ==================================================================
int parse_input_SIMGEN_DUMP(char **WORDS,int keySource) {

  // Created Jul 21 2020 [refactor]
  // Read/parse SIMGEN_DUMP[ALL] info
  //
  // Apr 16 2021: check SIMGEN_DUMPALL SWITCH 

  int  ivar, NVAR=0, N=0 ;
  bool LRD = false, LRD_COMMA_SEP=false, LRD_SPACE_SEP=false ;
  char *varName ;
  char fnam[] = "parse_input_SIMGEN_DUMP";

  // ------------- BEGIN ------------

  if ( keyMatchSim(1, "PRESCALE_SIMGEN_DUMP SIMGEN_DUMP_PRESCALE", 
		   WORDS[0], keySource) ) {
    N++ ; sscanf(WORDS[N] , "%d", &INPUTS.PRESCALE_SIMGEN_DUMP ); 
  }
  else if ( keyMatchSim(1, "SIMGEN_DUMP", WORDS[0], keySource) ) {
    LRD = true ;
  }
  else if ( keyMatchSim(1, "SIMGEN_DUMPALL", WORDS[0], keySource) ) {
    LRD = true ;
    INPUTS.IFLAG_SIMGEN_DUMPALL = 1 ;

    // Apr 16 2021: check option to switch from DUMP to DUMPALL
    if ( strcmp(WORDS[1],"SWITCH") == 0 ) {
      INPUTS.IFLAG_SIMGEN_DUMPALL = 1 ;
      return(N+1) ; 
    }
  }

  // - - - - - - - 

  if ( !LRD ) { return(N); }

  // check of comma sep or space-sep format.
  if ( strstr(WORDS[N+1],COMMA) != NULL ) 
    { LRD_COMMA_SEP = true;  }
  else 
    { LRD_SPACE_SEP = true; }


  // - - - - -
  if ( LRD_SPACE_SEP ) {         
    // original space-separate list of variables after NVAR integer
    // e.g., SIMGEN_DUMP: 4 CID ZCMB RA DEC
    N++ ; sscanf(WORDS[N] , "%d", &NVAR);

    if(NVAR <=0 ) { INPUTS.IFLAG_SIMGEN_DUMPALL=0; }

    for(ivar=0; ivar < NVAR; ivar++ ) {
      varName = INPUTS.VARNAME_SIMGEN_DUMP[ivar] ;
      N++ ; sscanf(WORDS[N], "%s", varName );
    }
  } // end LRD_SPACE_SEP

  if ( LRD_COMMA_SEP ) {
    // Apr 2021: optional comma-sep list that doesn't need NVAR
    // e.g., SIMGEN_DUMP: CID,ZCMB,RA,DEC

    char *ptrSplit[MXSIMGEN_DUMP];
    int MXCHARWD = MXCHARWORD_PARSE_WORDS;
    int MXVAR    = MXSIMGEN_DUMP ;
      
    N++ ;

    for(ivar=0; ivar < MXVAR; ivar++ ) 
      { ptrSplit[ivar] = INPUTS.VARNAME_SIMGEN_DUMP[ivar]; }

    splitString(WORDS[N], COMMA, MXVAR, &NVAR, ptrSplit);

  } // end LRD_COMMA_SEP


  // - - - - -

  INPUTS.NVAR_SIMGEN_DUMP = NVAR; // load global

  // check for alternate varnames
  for(ivar=0; ivar < NVAR; ivar++ ) {
    varName = INPUTS.VARNAME_SIMGEN_DUMP[ivar] ;      
    checkAlternateVarNames_HOSTLIB(varName);
  }

  return(N);

} // end parse_input_SIMGEN_DUMP

// ======================================================
int parse_input_GENMAG_SMEAR_SCALE(char **WORDS, int keySource ) {

  // July 2020: refactor
  // Examine argument of GENMAG_SMEAR_SCALE.
  //
  // Argument is either a global scale or a polynominal function:   
  //                                 
  // For these example user inputs                              
  //  GENMAG_SMEAR_SCALE: 1.3            ! global scale         
  //  GENMAG_SMEAR_SCALE(SALT2c) 0.9,0.3 ! scale = 0.9 + 0.3*c 
  // 
  // This function loads INPUTS.GENMAG_SMEAR_SCALE with              
  //     'NOVAR  1.3'
  //     'SALT2c 0.9,0.3'
  //
  
  int N=0;
  int LDMP = 0 ;
  char *KEYNAME = WORDS[0];
  char *STRVAL  = WORDS[1];
  char VARNAME[20], KEYTMP[100];
  char fnam[] = "parse_input_GENMAG_SMEAR_SCALE";

  // ----------- BEGIN ------------

  sprintf(KEYTMP, "%s", KEYNAME);

  // extract optional varname from () of key name
  extractStringOpt(KEYTMP,VARNAME); // return VARNAME 

  if ( !keyMatchSim(1, "GENMAG_SMEAR_SCALE",  KEYTMP, keySource) ) 
    { return(N); }

  N++ ; // do not explicitly read into local var; just use STRVAL

  if ( strlen(VARNAME) == 0 ) { sprintf(VARNAME,"NOVAR"); }

  sprintf(INPUTS.GENMAG_SMEAR_SCALE,"%s  %s", VARNAME, STRVAL);

  if ( LDMP ) {
    printf("\n xxx %s DEBUG DUMP -------------- \n", fnam );
    printf(" xxx input KEYNAME = '%s' \n", KEYNAME);
    printf(" xxx read STRVAL   = '%s' \n", STRVAL);
    printf(" xxx KEYTMP = '%s'   VARNAME = '%s' \n",  KEYTMP, VARNAME);
    fflush(stdout);
    //    debugexit(fnam);
  }

  return(N);

} // end parse_input_GENMAG_SMEAR_SCALE


// ==============================================================
void parse_GENMAG_SMEAR_MODELNAME(void) {

  // Split GENMAG_SMEAR_MODELSTRING by colon;
  // Left side of colon is MODELNAME, and right side of colon 
  // is a model argument stored as GENMAG_SMEAR_MODELARG.

  int  MEMC = MXCHAR_FILENAME * sizeof(char) ;
  char colon[] = ":" ;
  int  NSPLIT ;
  char *inString, *ptrSplit[2];
  char fnam[] = "parse_GENMAG_SMEAR_MODELNAME" ;

  // -------------- BEGIN -----------------

  inString    = (char*) malloc(MEMC);
  ptrSplit[0] = (char*) malloc(MEMC);
  ptrSplit[1] = (char*) malloc(MEMC);
  

  sprintf(inString,"%s", INPUTS.GENMAG_SMEAR_MODELNAME);

  splitString(inString, colon, 2,      // inputs               
	      &NSPLIT, ptrSplit );      // outputs             
  
  sprintf(INPUTS.GENMAG_SMEAR_MODELNAME, "%s", ptrSplit[0] );
  sprintf(INPUTS.GENMAG_SMEAR_MODELARG,  "%s", ptrSplit[1] );
  ENVreplace(INPUTS.GENMAG_SMEAR_MODELARG,fnam,1);

  free(inString); 
  free(ptrSplit[0]) ;  free(ptrSplit[1]) ;

  return ;

} // end parse_GENMAG_SMEAR_MODELNAME


// *****************************************************
int parse_input_TAKE_SPECTRUM(char **WORDS, int keySource, FILE *fp) {

  // Created July 2020 [refactored]
  // 
  // Works with SPECTROGRAPH to specify spectra w.r.t peakMJD,
  // instead of at aboslute MJD (simlib) values.
  //
  // Read from file pointer if *fp != NULL; e.g., from SIMLIB header.
  // Else sscanf(WORDS ...
  //
  // WORDS[0] = TAKE_SPECTRUM   or
  //            TAKE_SPECTRUM(FIELDNAME)
  //
  // WORDS[1] = TREST([Tmin]:[Tmax])  or
  //            TOBS([Tmin]:[Tmax])   or
  //            HOST
  //
  // WORDS[2] = SNR_ZPOLY([a0],[a1],[a2],,,,)  or
  //            TEXPOSE_ZPOLY([a0],[a1],[a2],,,,) 
  //       = a0 + a1*z + a2*z^2 + ...
  //
  // If WARP_SPECTRUM_STRING is not blank, parse this option
  // to apply mis-calibration vs. wavelength using polyFun of wavelength.
  //
  // Function returns number or words read.
  //
  // Note that colon denotes a range, while comma separates a list.
  //
  // Mar 23 2019: 
  //  + use parse_GENPOLY to allow arbitrary poly-order.
  //  + pass WARP_SPECTRUM_STRING
  //
  // June 28 2019: allow HOST argument
  //
  // May 29 2020: 
  //   + if reading SIMLIB header, skip epochs outside GENRANGE_TREST
  //
  // Jul 30 2020: fix ABORT logic for TAKE_SPECTRUM keys.
  // Feb 24 2021; fix TREST cut to work with GENLC.REDSHIFT or header z.
  // Apr 01 2021: parse MJD and FIELD
  // May 11 2021: abort if SNR_ZPOLY is used with HOST
  //      [because SNR calc works only with SN]
  //
  // May 27 2021: check for FIELD arg in WORDS(0)
  // Jun 21 2021: FIELD check only for keySource == KEYSOURCE_ARG
  
  bool READ_fp = (fp != NULL);
  int  NTAKE = NPEREVT_TAKE_SPECTRUM ;
  GENPOLY_DEF *GENLAMPOLY_WARP  = &INPUTS.TAKE_SPECTRUM[NTAKE].GENLAMPOLY_WARP ;
  GENPOLY_DEF *GENZPOLY_TEXPOSE = &INPUTS.TAKE_SPECTRUM[NTAKE].GENZPOLY_TEXPOSE;
  GENPOLY_DEF *GENZPOLY_SNR     = &INPUTS.TAKE_SPECTRUM[NTAKE].GENZPOLY_SNR ;
  char *WARP_SPECTRUM_STRING    =  INPUTS.WARP_SPECTRUM_STRING;
  int  N=0;  // track number of WORDS read
  
  float *ptrRange, *ptrLam, *ptrStep ;
  char string0[80], string1[80], string2[80], string3[80]; 
  char *ptrSplit[4], strValues[4][20], *ptrFrame ;
  int  NSPLIT, i ;
  bool IS_REST=false, IS_OBS=false, IS_MJD=false;
  bool IS_HOST = false;
  char stringTmp[80], stringOpt[200];
  char fnam[] = "parse_input_TAKE_SPECTRUM" ;

  // ----------- BEGIN -----------

  /* xxxxxxx mark delete 9/08/2021 xxxxxx
     can't make test since SIMLIB and input file are parsed same way
  if ( keySource == KEYSOURCE_FILE && INPUTS.USE_SIMLIB_SPECTRA ) {
    sprintf(c1err,"Cannot mix TAKE_SPECTRUM keys in sim-input & SIMLIB.") ;
    sprintf(c2err,"Remove one of these TAKE_SPECTRUM sources.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  xxxxxxxxxx */

  // init TAKE_SPECTRUM structure 
  for(i=0; i < 2; i++)  { 
    INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_RANGE[i]  = -99.0;
    INPUTS.TAKE_SPECTRUM[NTAKE].SNR_LAMRANGE[i] = 0.0 ;  
  }

  init_GENPOLY(GENLAMPOLY_WARP  ) ; 
  init_GENPOLY(GENZPOLY_TEXPOSE ) ; 
  init_GENPOLY(GENZPOLY_SNR     ) ; 

  INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_EPOCH  = 0 ;
  INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_LAMBDA = 0 ;
  INPUTS.TAKE_SPECTRUM[NTAKE].OPT_TEXPOSE      = 0 ;
  INPUTS.TAKE_SPECTRUM[NTAKE].FIELD[0]  = 0 ;

  // before reading TAKE_SPECTRUM options from file, 
  // parse optional WARP_SPECTRUM_STRING that was previously read.
  if ( strlen(WARP_SPECTRUM_STRING) > 0 ) {
    char warpOpt[100];  sprintf(warpOpt,"%s", WARP_SPECTRUM_STRING);
    extractStringOpt(warpOpt,stringOpt); // return stringOpt
    if ( strcmp(warpOpt,"LAMPOLY") != 0 ) {
      sprintf(c1err, "%s is invalid WARP_SPECTRUM arg", WARP_SPECTRUM_STRING);
      sprintf(c2err, "Expected LAMPOLY(,,,)");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    parse_GENPOLY(stringOpt, "wave", GENLAMPOLY_WARP, fnam);
    INPUTS.NWARP_TAKE_SPECTRUM++ ;
  }
  

  // for input file & command line arg, 
  // check for FIELD arg in (); e.g., TAKE_SPECTRUM(X1)
  if ( keySource == KEYSOURCE_ARG ) {
    sprintf(string0, "%s", WORDS[0]);
    char *FIELD = INPUTS.TAKE_SPECTRUM[NTAKE].FIELD ;
    extractStringOpt(string0,FIELD); // return FIELD
  }

  // ----------------------------------------------
  // read 1st arg and parse as either TREST or TOBS

  if ( READ_fp ) 
    { N++; readchar(fp, string1); }
  else
    { N++ ; sscanf(WORDS[N], "%s", string1); }

  sprintf(stringTmp, "%s", string1);
  extractStringOpt(stringTmp,stringOpt); // return stringOpt; 
  if ( strcmp(stringTmp,"TREST") == 0 ) {
    INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_EPOCH = GENFRAME_REST ;
    sprintf(INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_FRAME,"REST");
    IS_REST = true ;
  }

  else if ( strcmp(stringTmp,"TOBS") == 0 ) {
    INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_EPOCH = GENFRAME_OBS ;
    sprintf(INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_FRAME,"OBS");
    IS_OBS = true ;
  }

  else if ( strcmp(stringTmp,"MJD") == 0 ) {
    INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_EPOCH = GENFRAME_MJD ;
    sprintf(INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_FRAME,"MJD"); 
    IS_MJD = true ;
  }

  else if ( strcmp(stringTmp,"HOST") == 0 ) {
    INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_EPOCH = GENFRAME_HOST ;
    sprintf(INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_FRAME,"HOST");
    IS_HOST = true ;
    INPUTS.NHOST_TAKE_SPECTRUM++ ;
  }
  else if ( strcmp(stringTmp,"TEMPLATE_TEXPOSE_SCALE") == 0 ) {
    sscanf(stringOpt, "%f", 
	   &INPUTS.TAKE_SPECTRUM_TEMPLATE_TEXPOSE_SCALE);
    return(N) ;
  }
  else if ( strcmp(stringTmp,"NONE") == 0 ) {
    // turn off all spectra with command line arg: "TAKE_SPECTRUM NONE"
    INPUTS.NHOST_TAKE_SPECTRUM = 0;
    NPEREVT_TAKE_SPECTRUM = 0 ;
    SPECTROGRAPH_USEFLAG  = 0 ;
    return(N) ;
  }
  else {
    sprintf(c1err, "Cannot parse '%s' after TAKE_SPECTRUM key.",string1);
    sprintf(c2err, "Expecting TREST or TOBS string" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
 
  // get epoch range from colon-seperated values in stringOpt
  ptrRange = INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_RANGE ;
  ptrFrame = INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_FRAME ;


  if ( IS_HOST ) {
    ptrRange[0] = ptrRange[1] = 9999.0 ;  
  }
  else {
    // SN spectrum

    for(i=0; i < 4; i++ ) 
      { ptrSplit[i] = strValues[i];    ptrRange[i] = -9.0 ; }

    splitString(stringOpt, COLON, 4,      // inputs               
		&NSPLIT, ptrSplit );      // outputs             

    if ( NSPLIT < 1 || NSPLIT > 3 ) {
      sprintf(c1err, "\n   Found %d colon-separated values in '%s'", 
	      NSPLIT, string1); 
      sprintf(c2err, "but expected 1 or 2 or 3 values.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }   

    // load TREST_RANGE or TOBS_RANGE or MJD_RANGE
    for(i=0; i < NSPLIT; i++ ) 
      { sscanf( strValues[i] ,  "%f", &ptrRange[i] );  }

    // if only one time given, set range to delta function.
    if ( NSPLIT == 1 ) { ptrRange[1] = ptrRange[0];  }


  } // end SN spectrum


  // - - - - - - - - - - - - -  -
  // if reading SIMLIB header, apply Trest cut to 
  // skip epochs outside GENRANGE_TREST
  if ( !IS_HOST && INPUTS.USE_SIMLIB_SPECTRA ) {
    double z, Tmin, Tmax, Trest=9999.;  int OPT_FRAME;

    if ( INPUTS.USE_SIMLIB_REDSHIFT ) 
      { z = SIMLIB_HEADER.GENRANGE_REDSHIFT[0] ; }
    else 
      { z  = GENLC.REDSHIFT_CMB; }
    Tmin = INPUTS.GENRANGE_TREST[0] ;
    Tmax = INPUTS.GENRANGE_TREST[1] ;
    OPT_FRAME = INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_EPOCH ;
    if ( OPT_FRAME == GENFRAME_OBS ) 
      { Trest = ptrRange[0] / (1.0+z); }
    else if ( OPT_FRAME == GENFRAME_REST) 
      { Trest = ptrRange[0] ; }
    else {
      sprintf(c1err, "Invalid OPT_FRAME = %d (%s)",   OPT_FRAME, ptrFrame);
      sprintf(c2err, "epoch range: %f to %f \n", 
	      ptrRange[0], ptrRange[1] ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
    }

    if ( Trest < Tmin ) { return(N); }
    if ( Trest > Tmax ) { return(N); }

  }  // end if block

  // -------------------------------
  //  parse string2 for SNR or TEXPOSE
  // -------------------------------

  if ( READ_fp ) 
    { N++; readchar(fp, string2); }
  else
    { N++ ; sscanf(WORDS[N], "%s", string2); }

  sprintf(stringTmp, "%s", string2);
  extractStringOpt(stringTmp,stringOpt); // return stringOpt
  if ( strcmp(stringTmp,"SNR_ZPOLY") == 0 || strcmp(stringTmp,"SNR") == 0 ) {
    INPUTS.TAKE_SPECTRUM[NTAKE].OPT_TEXPOSE = 2;
    parse_GENPOLY(stringOpt, "SNR", GENZPOLY_SNR, fnam);

    if ( IS_HOST ) {
      sprintf(c1err,"Cannot use SNR_ZPOLY with HOST spectrum.");
      sprintf(c2err,"Only TEXPOSE_ZPOLY allowed for HOST");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }
  else if ( strcmp(stringTmp,"TEXPOSE_ZPOLY") == 0 ) {
    INPUTS.TAKE_SPECTRUM[NTAKE].OPT_TEXPOSE = 1;
    parse_GENPOLY(stringOpt, "TEXPOSE", GENZPOLY_TEXPOSE, fnam);
  }
  else {
    sprintf(c1err, "Cannot parse '%s' after TAKE_SPECTRUM key.",string1);
    sprintf(c2err, "Expecting SNR_ZPOLY or TEXPOSE_ZPOLY string" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }



  // - - - - - - - - - - - - 
  // for SNR option, also read SNR_LAMREST or SNR_LAMOBS

  if ( INPUTS.TAKE_SPECTRUM[NTAKE].OPT_TEXPOSE == 2 ) {

    if ( READ_fp ) 
      { N++; readchar(fp, string3); }
    else
      { N++ ; sscanf(WORDS[N], "%s", string3); }
    
    sprintf(stringTmp, "%s", string3);
    extractStringOpt(stringTmp,stringOpt); // return stringOpt

    splitString(stringOpt, COLON, 5,      // inputs               
		&NSPLIT, ptrSplit );      // outputs             

    if ( strcmp(stringTmp,"SNR_LAMREST") == 0 ||
	 strcmp(stringTmp,"LAMREST_SNR") == 0 ) 
      { INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_LAMBDA = GENFRAME_REST ; }
    else if ( strcmp(stringTmp,"SNR_LAMOBS") == 0 ||
	      strcmp(stringTmp,"LAMOBS_SNR") == 0   ) 
      { INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_LAMBDA = GENFRAME_OBS ; }
    else {
      sprintf(c1err,"Invalid key '%s' for LAMBDA-RANGE.", string3);
      sprintf(c2err,"Must be SNR_LAMREST or SNR_LAMOBS");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( NSPLIT != 2 ) {
      sprintf(c1err, "\n   Found %d colon-separated values in '%s'", 
	      NSPLIT, string3);
      sprintf(c2err, "but expected %d values.", 2);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    ptrLam = INPUTS.TAKE_SPECTRUM[NTAKE].SNR_LAMRANGE ;
    sscanf( strValues[0] , "%f", &ptrLam[0] );  // load LAMRANGE
    sscanf( strValues[1] , "%f", &ptrLam[1] ); 

  }  

  // - - - - - -

  int LDMP = 0 ;
  if ( LDMP ) {
    printf(" xxx ----------------------- \n");
    printf(" xxx NPEREVT_TAKE_SPECTRUM = %d \n", NPEREVT_TAKE_SPECTRUM);
    printf(" xxx EPOCH_RANGE = %6.2f to %6.2f \n", 
	   INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_RANGE[0],
	   INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_RANGE[1] );

    fflush(stdout);
  }

  if ( IS_MJD && ptrRange[2] > 0.01 )  
    { expand_TAKE_SPECTRUM_MJD(ptrRange) ; }
  else
    { NPEREVT_TAKE_SPECTRUM++ ; }

  // - - - - 
  if ( NPEREVT_TAKE_SPECTRUM >= MXPEREVT_TAKE_SPECTRUM ) {
    sprintf(c1err, "%d TAKE_SPECTRUM keys exceeds bound", 
	    NPEREVT_TAKE_SPECTRUM );
    sprintf(c2err,"Check MXPEREVT_TAKE_SPECTRUM in snlc_sim.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return(N);

} // end parse_input_TAKE_SPECTRUM


// ***********************************************
void expand_TAKE_SPECTRUM_MJD(float *MJD_RANGE) {

  // Created Apr 5 2021
  // Called for e.g.,
  //   TAKE_SPECTRUM: MJD(55000:55100:10)
  //
  // to internally expand INPUTS_TAKE_SPECTRUM struct to include
  // all 11 MJD: 55000, 55010, 55020, ... 55100.
  // Note that OPT_TEXPOSE must be 1 (TEXPOSE_ZPOLY).
 
  int   NTAKE_ORIG  = NPEREVT_TAKE_SPECTRUM;
  int   NTAKE       = NPEREVT_TAKE_SPECTRUM;
  int   OPT_TEXPOSE = INPUTS.TAKE_SPECTRUM[NTAKE].OPT_TEXPOSE;
  float MJD, MJD_STEP, MJD_MIN, MJD_MAX ;

  GENPOLY_DEF GENPOLY_WARP  ;
  GENPOLY_DEF GENPOLY_TEXPOSE ;
  GENPOLY_DEF GENPOLY_SNR     ;

  char fnam[] = "expand_TAKE_SPECTRUM_MJD";

  // ------- BEGIN -------

  MJD_MIN  = MJD_RANGE[0];  MJD_MAX=MJD_RANGE[1];
  MJD_STEP = MJD_RANGE[2];

  if ( OPT_TEXPOSE != 1 ) {
    sprintf(c1err,"Invalid OPT_TEXPOSE = %d", OPT_TEXPOSE);	    
    sprintf(c2err,"Must use TEXPOSE_ZPOLY after "
	    "TAKE_SPECTRUM: MJD(min:max:step)");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // store copies of polynom functions
  copy_GENPOLY(&INPUTS.TAKE_SPECTRUM[NTAKE].GENLAMPOLY_WARP,  &GENPOLY_WARP);
  copy_GENPOLY(&INPUTS.TAKE_SPECTRUM[NTAKE].GENZPOLY_TEXPOSE, &GENPOLY_TEXPOSE);
  copy_GENPOLY(&INPUTS.TAKE_SPECTRUM[NTAKE].GENZPOLY_SNR,     &GENPOLY_SNR);

  //  NTAKE--; // subtract 1 since it will be added again here
  for (MJD=MJD_MIN; MJD <= MJD_MAX; MJD += MJD_STEP ) {

    if ( NTAKE < MXPEREVT_TAKE_SPECTRUM ) {
      INPUTS.TAKE_SPECTRUM[NTAKE].OPT_TEXPOSE    = OPT_TEXPOSE ;
      INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_RANGE[0] = MJD;
      INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_RANGE[1] = MJD;
      INPUTS.TAKE_SPECTRUM[NTAKE].OPT_FRAME_EPOCH = GENFRAME_MJD ;
      sprintf(INPUTS.TAKE_SPECTRUM[NTAKE].EPOCH_FRAME,"MJD"); 

      copy_GENPOLY(&GENPOLY_WARP, 
		   &INPUTS.TAKE_SPECTRUM[NTAKE].GENLAMPOLY_WARP );
      copy_GENPOLY(&GENPOLY_TEXPOSE, 
		   &INPUTS.TAKE_SPECTRUM[NTAKE].GENZPOLY_TEXPOSE);
      copy_GENPOLY(&GENPOLY_SNR, 
		   &INPUTS.TAKE_SPECTRUM[NTAKE].GENZPOLY_SNR );
    }

    NTAKE++ ;
  }

  if ( NTAKE >= MXPEREVT_TAKE_SPECTRUM ) {
    sprintf(c1err, "%d TAKE_SPECTRUM keys exceeds bound of %d", 
	    NTAKE, MXPEREVT_TAKE_SPECTRUM );
    sprintf(c2err,"Check MXPEREVT_TAKE_SPECTRUM in snlc_sim.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  NPEREVT_TAKE_SPECTRUM = NTAKE ;
  return;

} // end expand_TAKE_SPECTRUM_MJD



// *****************************************
int parse_input_SIMSED(char **WORDS, int keySource) {

  int MXPAR = MXPAR_SIMSED ;
  int NPAR, N=0;
  char PARNAME[60], SIMSED_STRINGOPT[80];
  char fnam[] = "parse_input_SIMSED" ;

  // ----------- BEGIN -----------


  if ( keyMatchSim(1, "SIMSED_USE_BINARY", WORDS[0], keySource) ) {
    N++ ; sscanf(WORDS[N], "%d", &INPUTS.USE_BINARY_SIMSED );
  }
  else if ( keyMatchSim(1, "SIMSED_PATH_BINARY", WORDS[0], keySource) ) {
    N++ ; sscanf(WORDS[N], "%s", INPUTS.PATH_BINARY_SIMSED );
  }
  else if ( keyMatchSim(MXPAR, "SIMSED_PARAM", WORDS[0], keySource) ) {
    N += parse_input_SIMSED_PARAM(WORDS); 
  }
  else if ( keyMatchSim(1, "SIMSED_SHAPEPAR", WORDS[0], keySource) ) {
    N += parse_input_SIMSED_PARAM(WORDS); 
    INPUTS.IPAR_SIMSED_SHAPE = INPUTS.NPAR_SIMSED-1 ;
  }
  else if ( keyMatchSim(1, "SIMSED_COLORPAR", WORDS[0], keySource) ) {
    N += parse_input_SIMSED_PARAM(WORDS); 
    INPUTS.IPAR_SIMSED_COLOR = INPUTS.NPAR_SIMSED-1 ;
  }
  else if ( keyMatchSim(1, "SIMSED_COLORLAW", WORDS[0], keySource) ) {
    N += parse_input_SIMSED_PARAM(WORDS); 
    INPUTS.IPAR_SIMSED_CL = INPUTS.NPAR_SIMSED-1 ;
  }

  else if ( keyMatchSim(MXPAR,"SIMSED_GRIDONLY", WORDS[0],keySource) ) {
    N++ ; sscanf(WORDS[N],"%s", PARNAME);

    // in case we have VARNAME(STRINGOPT), extract optional STRINGOPT
    extractStringOpt(PARNAME,SIMSED_STRINGOPT);

    if ( strcmp(PARNAME,"SEQUENTIAL") == 0  ) {
      INPUTS.NPAR_SIMSED = 0 ;  // no interp pars => treat as baggage
      INPUTS.OPTMASK_SIMSED    = OPTMASK_GEN_SIMSED_GRIDONLY ;
    }
    else {
      // PARNAME is one of the SIMSED params, but GRIDONLY
      NPAR = INPUTS.NPAR_SIMSED ;
      sprintf(INPUTS.PARNAME_SIMSED[NPAR], "%s", PARNAME );
      INPUTS.GENFLAG_SIMSED[NPAR] = 
	OPTMASK_GEN_SIMSED_PARAM + OPTMASK_GEN_SIMSED_GRIDONLY ;
      sprintf(INPUTS.KEYWORD_SIMSED[NPAR],"SIMSED_GRIDONLY");

      INPUTS.NPAR_SIMSED++ ; 
      INPUTS.NPAR_SIMSED_GRIDONLY++ ; 

      parse_input_SIMSED_SUBSET(PARNAME,SIMSED_STRINGOPT) ;
    }

  } // end SIMSED_GRIDONLY 

  else if ( keyMatchSim(1,"SIMSED_REDCOR SIMSED_COV", WORDS[0], keySource) ) {
    N += parse_input_SIMSED_COV(WORDS, keySource);
  }
  
  // - - - -
  NPAR       = INPUTS.NPAR_SIMSED ;
  if ( NPAR >= MXPAR_SIMSED ) {
    sprintf(c1err, "NPAR_SIMSED=%d exceeds array bound", NPAR );
    sprintf(c2err, "Check SIMSED_PARAM keywords in input file." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  
  return(N);

} // end parse_input_SIMSED

// *****************************************
int parse_input_SIMSED_PARAM(char **WORDS) {
  // July 2010 [refactor]
  // call this function when any of the
  // SIMSED_[PARAM/SHAPEPAR/COLORPAR] keys is found.
  int NPAR, N=0;
  char fnam[] = "parse_input_SIMSED_PARAM" ;
  // ------------- BEGIN -----------
  NPAR = INPUTS.NPAR_SIMSED ;
  N++ ; sscanf(WORDS[N], "%s", INPUTS.PARNAME_SIMSED[NPAR] );
  INPUTS.GENFLAG_SIMSED[NPAR] = OPTMASK_GEN_SIMSED_PARAM; //continuous interp
  sprintf(INPUTS.KEYWORD_SIMSED[NPAR],"SIMSED_PARAM"); 
  INPUTS.NPAR_SIMSED++ ; 
  INPUTS.NPAR_SIMSED_PARAM++ ; 
  return(N) ;
} // end of parse_input_SIMSED_PARAM


// *******************************************      
void parse_input_SIMSED_SUBSET(char *parName, char *stringOpt) {

  // if stringOpt = 2,4,5,6 then split string and load
  // each grid index into global array.

  // USE parse_commaSepList(...)

  int NINDEX_SUBSET = 0;  // number of indices for subset
  int MXINDEX_SUBSET = MXPAR_SIMSED;
  int i, INDEX_SUBSET ;
  char **ptrSplit;
  char fnam[] = "parse_input_SIMSED_SUBSET";
  // ------ BEGIN ------

  INPUTS.NINDEX_SUBSET_SIMSED_GRIDONLY = 0 ;
  if ( strlen(stringOpt) == 0 ) { return ; }

  if ( !IS_INDEX_SIMSED(parName) )  {
    sprintf(c1err,"Cannot specify index subset for parName = %s", parName);
    sprintf(c2err,"parName string must include INDEX/INDX/index/indx");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  parse_commaSepList("GRINDONLY_SUBSET", stringOpt, MXINDEX_SUBSET, 8,
		     &NINDEX_SUBSET, &ptrSplit) ;

  for(i=0; i < NINDEX_SUBSET; i++ ) {
    sscanf(ptrSplit[i], "%d", &INDEX_SUBSET);
    INPUTS.INDEX_SUBSET_SIMSED_GRIDONLY[i] = INDEX_SUBSET ;
  }

  INPUTS.NINDEX_SUBSET_SIMSED_GRIDONLY = NINDEX_SUBSET ;

} // end parse_input_SIMSED_SUBSET

// *******************************************      
bool keyContains_SIMSED_PARAM(char *KEYNAME) {

  // check if input *KEYBANE contains most recently 
  // read SIMSED parameter;
  // E.g. KEYNAME = "GENPEAK_MNI:" and most recent SIMSED param is MNI
  //   --> returns true.

  int NPAR    = INPUTS.NPAR_SIMSED ;
  char *parName ;
  char fnam[] = "keyContains_SIMSED_PARAM" ;

  // ----------- BEGIN ------------

  if ( NPAR == 0 ) return(false);

  // require continuous distribution
  int GENFLAG = INPUTS.GENFLAG_SIMSED[NPAR-1];
  if ( (GENFLAG & OPTMASK_GEN_SIMSED_PARAM ) == 0 ) { return(false); }

  parName = INPUTS.PARNAME_SIMSED[NPAR-1];
  if ( strstr(KEYNAME,parName) != NULL ) 
    { return(true); }
  else
    { return(false); }

} // end keyContains_SIMSED_PARAM

// *******************************************
int parse_input_SIMSED_COV(char **WORDS, int keySource) {

  // Created July 2020 [refactor]
  // Store information for off-diagonal covariances used
  // to generate correlated random Gaussians.
  // Called if "SIMSED_REDCOR(var1,var2):" is found.
  // varList = var1,var2 needs to be parsed.
  //
  // Inputs:
  //   WORDS     = array of sim-input words, starting with valid key
  //   keySource = KEUSOURCE_FILE or KEYSOURCE_ARG
  //
  // WORDS[0] is a KEY with pair of varnames in parentheses;
  // e.g.  SIMSED_COV(MNI,COSANGLE):
  //

#define OPTREAD_REDCOR_SIMSED 1
#define OPTREAD_COV_SIMSED    2

  int N=0,  NPAIR = INPUTS.NPAIR_SIMSED_COV; // == Noffdiag/2
  double sigLo, sigHi, SIGLIST[2], RHO;
  int  OPT, MXNAME=4, NNAME, ipar, IPAR[2];
  int  IPAR0, IPAR1, IPARTMP, j, NFOUND=0 ;
  char KEYNAME[100], varList[100] ;
  char *ptrName[2], varNames[2][40];
  char fnam[] = "parse_input_SIMSED_COV" ;

  // -------------- BEGIN --------------

  // extract varList pair from parentheses of key
  sprintf(KEYNAME, "%s", WORDS[0]); 
  extractStringOpt(KEYNAME,varList);  
  
  // check keyname to see if next value is REDCOV or COV
  if ( strstr(KEYNAME,"SIMSED_REDCOV") != NULL ) 
    { OPT = OPTREAD_REDCOR_SIMSED; }
  else
    { OPT = OPTREAD_COV_SIMSED; }

  if ( NPAIR >= MXPAR_SIMSED ) { goto COUNT; }

  // read either REDCOR or COV value
  N++; sscanf(WORDS[N], "f", &INPUTS.COVPAIR_SIMSED_COV[NPAIR]) ;  

  // break varList into two SIMSED varnames
  ptrName[0] = varNames[0] ;
  ptrName[1] = varNames[1] ;
  splitString(varList, COMMA, MXNAME, &NNAME, ptrName);

  // find which IPAR
  for(ipar=0; ipar < INPUTS.NPAR_SIMSED ; ipar++ ) {
    for(j=0; j < 2; j++ ) {
      if ( strcmp(INPUTS.PARNAME_SIMSED[ipar],varNames[j]) == 0 )
	{ IPAR[j] = ipar; NFOUND++; }
    }
  } // end loop over all SIMSED params
    
  if ( NFOUND != 2 ) {
    sprintf(c1err, "Found %d matches for REDCOR_SIMSED", NFOUND);
    sprintf(c2err, "Check '%s' in the sim-input file", varList);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  IPAR0 = IPAR[0];  IPAR1=IPAR[1];

  // abort on diagonal terms
  if ( IPAR0 == IPAR1 ) {
    sprintf(c1err, "Cannot define diagonal COV here" );
    sprintf(c2err, "Use GENSIGMA_%s instead.", varNames[0] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  for(j=0; j<2; j++ ) {
    IPARTMP = IPAR[j];

    INPUTS.IPARPAIR_SIMSED_COV[NPAIR][j] = IPARTMP ;

    // require symmetric Gaussian
    sigLo = INPUTS.GENGAUSS_SIMSED[IPARTMP].SIGMA[0];
    sigHi = INPUTS.GENGAUSS_SIMSED[IPARTMP].SIGMA[1];
    SIGLIST[j] = sigLo;

    if ( sigLo != sigHi ) {
      sprintf(c1err, "Cannot apply correlation to asymmetric Gaussian");
      sprintf(c2err, "Check input '%s' sigLo,sigHi = %.5f, %.5f\n", 
	      INPUTS.GENGAUSS_SIMSED[IPARTMP].NAME,
	      sigLo, sigHi );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  // if input is COV, divide by sigma to get reduced corr.
  if ( OPT == OPTREAD_REDCOR_SIMSED ) {
    INPUTS.COVPAIR_SIMSED_COV[NPAIR] *= ( SIGLIST[0]*SIGLIST[1] );
  }

  // check for valid reduced cov
  RHO = INPUTS.COVPAIR_SIMSED_COV[NPAIR] / ( SIGLIST[0]*SIGLIST[1] );
  if ( fabs(RHO) > 1.0000001 ) {
    sprintf(c1err, "Invalid REDCOR = %f", RHO );
    sprintf(c2err, "between '%s' and '%s' ", varNames[0], varNames[1]);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

 COUNT:
  INPUTS.NPAIR_SIMSED_COV++ ;

  return(N) ;

} // end parse_input_SIMSED_COV

// ************************************************** 
void parse_input_OBSOLETE(char **WORDS, int keySource ) {

  // Created July 2020
  // ABORT on obsolete keys, giving message about correct key to use.
  //
  // Note that abort call is via
  //    legacyKey_abort(fnam, key_obsolete, key_new)
  // to print name of valid key.
  //
  // Dec 23 2020: abort if CONFIG file key is found.

  char fnam[] = "parse_input_OBSOLETE";

  // ---------- BEGIN -------------

  if ( keyMatchSim(1, "NGRID_LUMIPAR", WORDS[0], keySource) )
    { legacyKey_abort(fnam, "NGRID_LUMIPAR:", "NGRID_SHAPEPAR:"); }

  else if ( keyMatchSim(1, "GEN_SNDATA_SIM", WORDS[0], keySource) )
    { legacyKey_abort(fnam, "GEN_SNDATA_SIM:", "FORMAT_MASK:"); }

  else if ( keyMatchSim(1, "EXTINC_MILKYWAY", WORDS[0], keySource) )
    { legacyKey_abort(fnam, "EXTINC_MILKYWAY:", "OPT_MWEBV:"); }

  else if ( keyMatchSim(1, "RVMW", WORDS[0], keySource) )
    { legacyKey_abort(fnam, "RVMW:", "RV_MWCOLORLAW:"); }

  else if ( keyMatchSim(1, "HUMAN_SEARCHEFF_OPT", WORDS[0], keySource) )
    { legacyKey_abort(fnam, "HUMAN_SEARCHEFF_OPT:", "SEARCHEFF_SPEC_FILE:");}


  else if ( keyMatchSim(1, "SPECTYPE", WORDS[0], keySource) )
    { legacyKey_abort(fnam, "SPECTYPE:", ""); }

  else if ( keyMatchSim(1, "SMEARFLAG_HOSTGAL", WORDS[0], keySource) )
    { legacyKey_abort(fnam, "SMEARFLAG_HOSTGAL:", ""); }

 
  bool IS_CONFIG_FILE = 
    ( keyMatchSim(1, "BATCH_INFO", WORDS[0], keySource)  ||
      keyMatchSim(1, "CONFIG",     WORDS[0], keySource) ) ;
  if ( IS_CONFIG_FILE ) {
    sprintf(c1err,"'%s' key is for config/master file -->", WORDS[0]);
    sprintf(c2err,"This file is for submit_batch_jobs, not snlc_sim.exe");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  

  return;

} // end parse_input_OBSOLETE

// **************************************************
void  parse_input_GENPOP_ASYMGAUSS(void) {

  // Created Jun 2021
  // parse separate file with asym-Gauss keys for SALT2 color and stretch.
  // If model is G10_HIZ, then search for
  //
  // MODEL_NAME: G10_HIZ
  // GENPEAK_SALT2x1: xxx
  // GENSIGMA_SALT2x1: xxx xxx
  // GENRANGE_SALT2x1: xxx xxx
  // [and repeat for SALT2c]
  // MODEL_END:

  FILE *fp;
  char *ptrFile  = INPUTS.GENPOP_ASYMGAUSS_FILE ;
  char *ptrModel = INPUTS.GENPOP_ASYMGAUSS_MODEL ;
  int   MSKOPT_PARSE = MSKOPT_PARSE_WORDS_FILE ;
  int   keySource    = KEYSOURCE_FILE ;
  char  KEY_MODEL_NAME[] = "MODEL_NAME:" ;
  char  KEY_MODEL_END[]  = "MODEL_END:" ;
  bool  FOUND_MODEL = false, FOUND_END=false;
  bool  FOUND_c=false, FOUND_x1=false;
  bool  ISKEY_MODEL, ISKEY_END;

  int   iwd, NWD_TOT, NRD=0 ;
  char **WORDS, KEY_TMP[MXPATHLEN], MODEL_TMP[40] ;
  char fnam[] = "parse_input_GENPOP_ASYMGAUSS";

  // ----------- BEGIN ------------
  
  ENVreplace(ptrFile, fnam, 1);

  printf("  Read GENPOP_ASYMGAUSS model '%s' from \n     %s\n",
	 ptrModel, ptrFile); fflush(stdout);

  NWD_TOT = store_PARSE_WORDS(MSKOPT_PARSE, ptrFile);

  // transfer store_PARSSE_WORDS to local list to avoid
  // conflict with other parsing functions that use store_PARSE_WORDS
  WORDS = (char**) malloc( NWD_TOT * sizeof(char*) );
  for ( iwd=0; iwd < NWD_TOT; iwd++ ) {
    WORDS[iwd]    = (char*) malloc( 100*sizeof(char) ); 
    get_PARSE_WORD(0, iwd, WORDS[iwd] );
  }

  // - - - - - -
  for(iwd=0; iwd < NWD_TOT; iwd++ ) {
    
    sprintf(KEY_TMP, "%s", WORDS[iwd]);

    //    printf(" xxx read KEY_TMP[%d] = '%s' \n", iwd, KEY_TMP);
    ISKEY_MODEL = ( strcmp(KEY_TMP,KEY_MODEL_NAME) == 0 ) ;
    ISKEY_END   = ( strcmp(KEY_TMP,KEY_MODEL_END)  == 0 ) ;

    if  ( ISKEY_MODEL ) {
      if (  strcmp(ptrModel,WORDS[iwd+1])==0 )  { FOUND_MODEL = true; }
    }

    if ( !FOUND_MODEL ) { continue; }

    if ( ISKEY_END ) { FOUND_END = true; break; }

    // parse sim-input keys within model block
    //    printf(" xxx read iwd=%d of %d(%d): '%s' \n", 
    //	   iwd, NWD_TOT, PARSE_WORDS.NWD, KEY_TMP);  

    if ( strstr(WORDS[iwd],"SALT") == NULL ) { continue; }

    NRD = parse_input_GENGAUSS("SALT2c", &WORDS[iwd], keySource,
				&INPUTS.GENGAUSS_SALT2c );  
    if ( NRD > 0 ) { FOUND_c = true; }

    NRD = parse_input_GENGAUSS("SALT2x1", &WORDS[iwd], keySource,
				&INPUTS.GENGAUSS_SALT2x1 );  
    if ( NRD > 0 ) { FOUND_x1 = true; }

  }

  // - - -- 
  bool FOUND_ALL = FOUND_MODEL && FOUND_END && FOUND_c && FOUND_x1;
  if ( !FOUND_ALL ) {
    print_preAbort_banner(fnam);
    printf("  FOUND_MODEL = %d for '%s' \n", FOUND_MODEL, ptrModel);
    printf("  FOUND_END   = %d \n", FOUND_END);
    printf("  FOUND_[x1,c] = %d, %d \n", FOUND_x1, FOUND_c);

    sprintf(c1err,"Missed reading something in GENPOP file");
    sprintf(c2err,"Check FOUND_XXX logicals above");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // free local memory for file contents
  for ( iwd=0; iwd < NWD_TOT; iwd++ ) { free(WORDS[iwd]); }
  free(WORDS);

  //  debugexit(fnam);

  return ;

} // end parse_input_GENPOP_ASYMGAUSS


// *************************************************
int get_NON1A_MODELFLAG(char *GENMODEL) {

  // Created Mar 2016
  // Return MODELFLAG = 
  //  -9 -> not a nonIa model
  // MODEL_NON1ASED  for NON1A or NONIA or NON1ASED
  // MODEL_NON1AGRID for NON1AGRID or NONIAGRID
  //
  // Jan 11 2017: fix but and return -9 if NOT NON1a model.
  // Jul 03 2017: check for LCLIB
  // Sep 05 2018: allow for names like NON1ASED.ABC
  // Apr 18 2019: no more check for LCLIB

  char *ptr;
  char GENMODEL_LOCAL[80];
  int  ikey, NKEY=0;
  char MODELKEY_LIST[10][20];
  int  MODELFLAG_LIST[10];
  int  MODELFLAG = -9 ;

  // ------------- BEGIN -------------

  if ( strstr(GENMODEL,"NON") == NULL )  { return -9 ; }

  // copy GENMODEL into local var to make sure GENMODEL isn't modified.
  sprintf(GENMODEL_LOCAL, "%s", GENMODEL);

  // if there is a dot in the model name, remove everything
  // after the dot.
  ptr = strchr(GENMODEL_LOCAL,'.');
  if ( ptr != NULL ) { *ptr = '\0' ; }     

  // build list of model name keys to check
  sprintf(MODELKEY_LIST[NKEY], "NON1A");   
  MODELFLAG_LIST[NKEY] = MODEL_NON1ASED ;
  NKEY++ ;

  sprintf(MODELKEY_LIST[NKEY], "NONIA");   
  MODELFLAG_LIST[NKEY] = MODEL_NON1ASED ;
  NKEY++ ;

  sprintf(MODELKEY_LIST[NKEY], "NON1ASED");   
  MODELFLAG_LIST[NKEY] = MODEL_NON1ASED ;
  NKEY++ ;

  sprintf(MODELKEY_LIST[NKEY], "NONIASED");   
  MODELFLAG_LIST[NKEY] = MODEL_NON1ASED ;
  NKEY++ ;

  sprintf(MODELKEY_LIST[NKEY], "NON1AGRID");   
  MODELFLAG_LIST[NKEY] = MODEL_NON1AGRID ;
  NKEY++ ;

  sprintf(MODELKEY_LIST[NKEY], "NONIAGRID");   
  MODELFLAG_LIST[NKEY] = MODEL_NON1AGRID ;
  NKEY++ ;

  // - - - - 

  for(ikey=0; ikey < NKEY; ikey++ ) {
    if ( strcmp(GENMODEL_LOCAL,MODELKEY_LIST[ikey]) == 0 ) 
      { MODELFLAG = MODELFLAG_LIST[ikey] ; }
  }


  //  printf(" xxxx MODELFLAG(%s) = %d \n", GENMODEL, MODELFLAG);

  return(MODELFLAG);

} // end get_NON1A_MODELFLAG

// *********************************
void sim_input_override(void) {

  // July 2020 [refactor]
  // Parse command-line overrides using parse_input_key_driver util.
  //
  // argv[0] is the name of the program.
  // argv[1] is the name of the input file.
  // start parsing at argv[2].

  int IWD_START = 2,  iwd, NWD_READ, iwd_use ;
  char fnam[] = "sim_input_override" ;
  // --------- BEGIN ----------

  for(iwd = IWD_START; iwd < NARGV_LIST; iwd++ ) {
    NWD_READ = parse_input_key_driver(&ARGV_LIST[iwd],KEYSOURCE_ARG);

    /*
    printf(" xxx %s: iwd_use range is %d to %d (NWD_READ=%d)\n", 
	   fnam, iwd, iwd+NWD_READ, NWD_READ );
    */

    // set USE flag to mark valid command-line inputs
    if ( NWD_READ > 0 && NWD_READ < FLAG_NWD_ZERO ) {
      for(iwd_use = iwd; iwd_use < (iwd+NWD_READ+1); iwd_use++ )  { 
	USE_ARGV_LIST[iwd_use] = 1; 
      }
    }
    else if ( NWD_READ == FLAG_NWD_ZERO ) {
      USE_ARGV_LIST[iwd] = 1;  // Jun 17 2021
      continue ;
    }

    iwd += NWD_READ ;
  }

  // debugexit(fnam); // xxx REMOVE
  return ;

} // end sim_input_override


// **********************************************
void prep_user_input(void) {

  /**************************
   Dec 2007 R.Kessler
   after reading all user-input, construct additional 
   input variables.
  

  May 21, 2013: NGEN *= INPUTS.NGEN_SCALE

  May 26, 2013: check SALT2mu_FILE for alpha,beta values

  Oct 16 2013: apply INPUTS.NGEN_SCALE to INPUTS.NGEN_LC and NGENTOT_LC
               so that screen dump/readme includes this factor.

  Apr 21 2014: 
    fix bug initializing  GENMODEL_NAME; indx=0 to MXMODEL_INDEX-1
    instead of 0 to MXMODEL_INDEX that overwrote memory.
      [bug found by B. Patel].

  Oct 30 2014: remove reference to _non1a model.

  May 11 2015: call INIT_FUDGE_SNRMAX

  Mar 19 2016: 
    + prepare new input INPUTS.GENMAG_OFF_NON1A
    + replace most floats with double


  Jun 19 2017:
   + remove !MODEL_FIXMAG requirement for prep_genmag_offset().

  Aug 16 2017: abort of GENRANGE_REDSHIFT[1] > ZMAX_SNANA

  Sep 25 2017: move CUTWIN stuff into prep_user_CUTWIN()
  Jan 23 2018: check WRFLAG_COMPACT
  Sep 05 2018: set INPUTS.NON1ASED.PATH
  Oct 04 2018: set INPUTS.USE_HOSTLIB_GENZPHOT
  Feb 15 2019: turn off INPUTS.GENMAG_SMEAR_MODELNAME for SIMSED model.
  Dec 01 2019: allow COH scatter gmodel for NON1ASED and SIMSED models.
  Feb 06 2020: set DOGEN_AV for GRIDGEN
  Oct 16 2020: call prep_user_cosmology()
  Feb 21 2021: abort on FORMAT_MASK +=1, or legacy VERBOSE 
  oct 14 2021: set spectra bit of WRITE_MASK if spectrograph is used.

  *******************/

  int ifilt, ifilt_obs, IBLANK,i, j, indx, lentmp, NGEN, OPT, NTMP  ;
  int DOCHECK_FORMAT_MASK=1, ISRATE_LCLIB, INDEX_RATEMODEL;
  int USE_SMEAR_MODELNAME, ISCOH_SMEAR;

  float tmp_F, XGEN_F ; 
  char  *PTR_GENMODEL, *PTR_SMEAR_MODELNAME, vtmp[40]   ;
  char *input_file = INPUTS.INPUT_FILE_LIST[0] ;
  char fnam[] = "prep_user_input" ;

  // ---------- BEGIN ----------

  if ( INPUTS.USE_KCOR_REFACTOR == 2 )  { INPUTS.USE_KCOR_LEGACY = 0 ; }


  if ( INPUTS.DEBUG_FLAG != -1024 ) 
    { INPUTS_SEARCHEFF.OPTMASK_OPENFILE += OPENMASK_REQUIRE_DOCANA ; }

  // Feb 2015: replace ENV names in inputs
  ENVreplace(INPUTS.KCOR_FILE,fnam,1);  
  ENVreplace(INPUTS.SIMLIB_FILE,fnam,1);
  ENVreplace(INPUTS.HOSTLIB_FILE,fnam,1);
  ENVreplace(INPUTS.HOSTLIB_WGTMAP_FILE,fnam,1);
  ENVreplace(INPUTS.HOSTLIB_ZPHOTEFF_FILE,fnam,1);
  ENVreplace(INPUTS.HOSTLIB_SPECBASIS_FILE,fnam,1);
  ENVreplace(INPUTS.HOSTLIB_SPECDATA_FILE,fnam,1);
  ENVreplace(INPUTS.FLUXERRMODEL_FILE,fnam,1 );
  ENVreplace(INPUTS.HOSTNOISE_FILE,fnam,1 );
  ENVreplace(INPUTS.WRONGHOST_FILE,fnam,1 );
  ENVreplace(INPUTS.PATH_BINARY_SIMSED,fnam,1);
  ENVreplace(INPUTS.NON1ASED.PATH,fnam,1);
  ENVreplace(INPUTS.NON1AGRID_FILE,fnam,1);
  ENVreplace(INPUTS.NONLINEARITY_FILE,fnam,1);
  ENVreplace(INPUTS.WEAKLENS_PROBMAP_FILE,fnam,1);
  ENVreplace(INPUTS.STRONGLENS_FILE,fnam,1);
  ENVreplace(INPUTS.LCLIB_FILE,fnam,1);
  ENVreplace(INPUTS.MODELPATH,fnam,1);
  ENVreplace(PATH_USER_INPUT,fnam,1);
  ENVreplace(INPUTS.GENPDF_FILE,fnam,1);

  if ( strlen(INPUTS.PATH_SNDATA_SIM) > 0 ) {
    add_PATH_SNDATA_SIM(INPUTS.PATH_SNDATA_SIM);
    ENVreplace(INPUTS.PATH_SNDATA_SIM,fnam,1);
    sprintf(PATH_SNDATA_SIM, "%s", INPUTS.PATH_SNDATA_SIM);
  }

  ENVreplace(INPUTS_SEARCHEFF.USER_SPEC_FILE,fnam,1);
  ENVreplace(INPUTS_SEARCHEFF.USER_zHOST_FILE,fnam,1);
  ENVreplace(INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE,fnam,1);
  ENVreplace(INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE,fnam,1);
  ENVreplace(INPUT_ZVARIATION_FILE,fnam,1);

  INDEX_GENMODEL = 0 ;
  SUBINDEX_GENMODEL = 0 ; // May 2 2021

  sprintf(GENLC.SHAPEPAR_NAME,    "NULL");
  sprintf(GENLC.SHAPEPAR_GENNAME, "NULL");
  sprintf(GENLC.COLORPAR_NAME,    "NULL");
  sprintf(GENLC.SHAPEPAR2_NAME,   "NULL");
  sprintf(GENLC.COLORPAR2_NAME,   "NULL");
  sprintf(GENLC.DISTANCE_NAME,    "NULL");


  // Now identify MODEL_INDEX
  for ( indx = 1; indx < MXMODEL_INDEX; indx++ ) {
    for ( j=0; j < MXNAME_PER_MODEL ;  j++ ) {
      PTR_GENMODEL = GENMODEL_NAME[indx][j] ;
      if ( strcmp(PTR_GENMODEL,"NULL") == 0 ) { continue ; }

      lentmp = strlen(PTR_GENMODEL);
      if ( strncmp(INPUTS.MODELNAME,PTR_GENMODEL,lentmp) == 0 ) { 
	INDEX_GENMODEL = indx ; 

	// for SALT; check subindex for SALT2,SALT3,SALT4 ...
	if ( INDEX_GENMODEL == MODEL_SALT2 ) { SUBINDEX_GENMODEL=j; }
      }
    }
  }

  if ( INDEX_GENMODEL <= 0 ) {
    sprintf(c1err,"'%s' is not a valid genmag-model", INPUTS.MODELNAME);
    sprintf(c2err,"Check GENMODEL keyword in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  check_model_default( INDEX_GENMODEL ); // adjust INPUTS.GENMODEL if needed


  ISRATE_LCLIB = 0 ; 
  INDEX_RATEMODEL = INPUTS.RATEPAR.INDEX_MODEL ;
  if ( INDEX_RATEMODEL == INDEX_RATEMODEL_COSBPOLY ) { ISRATE_LCLIB=1; }
  if ( INDEX_RATEMODEL == INDEX_RATEMODEL_BPOLY    ) { ISRATE_LCLIB=1; }


  // finally, set MODELPATH, except for NON1A
  if ( strlen(INPUTS.MODELPATH) == 0  && INPUTS.NON1A_MODELFLAG < 0 ) {
    PTR_GENMODEL = GENMODEL_NAME[INDEX_GENMODEL][SUBINDEX_GENMODEL] ; 
    sprintf(INPUTS.MODELPATH,"%s/models/%s/%s", 
	    PATH_SNDATA_ROOT, PTR_GENMODEL, INPUTS.MODELNAME );
  }

  // check for optional over-ride of model path;
  // getenv(PRIVATE_MODELPATH_NAME is equivalent of 
  // $SNDATA_ROOT/models/[SALT2,mlcs2k2...]
  // Jan 12 2021: abort if SNANA_MODELPATH is set.
  if ( getenv(PRIVATE_MODELPATH_NAME) != NULL ) {

    sprintf(c1err,"ENV %s is no longer valid", PRIVATE_MODELPATH_NAME);
    sprintf(c2err,"Include full model path in GENMODEL arg.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 

    //sprintf(INPUTS.MODELPATH,"%s/%s", 
    //	    getenv(PRIVATE_MODELPATH_NAME), INPUTS.MODELNAME );
  }


  // disable HOSTLIB_MSKOPT if HOSTLIB_FILE isn't set.
  if ( IGNOREFILE(INPUTS.HOSTLIB_FILE) ) { INPUTS.HOSTLIB_MSKOPT=0; }

  // -------------------------------
  // set INPUTS.GENSOURCE

  sprintf(vtmp, "%s", INPUTS.GENSOURCE);
  if ( strcmp(vtmp,"RANDOM") == 0 ) 
    { GENLC.IFLAG_GENSOURCE = IFLAG_GENRANDOM ; }
  else if ( strcmp(vtmp,"GRID") == 0 ) 
    { GENLC.IFLAG_GENSOURCE = IFLAG_GENGRID  ; }

  // set OPT_SNXT
  if ( strcmp(INPUTS.GENSNXT,"CCM89") == 0 ) 
    { INPUTS.OPT_SNXT = OPT_SNXT_CCM89; }
  else if ( strcmp(INPUTS.GENSNXT,"SJPAR") == 0 ) 
    { INPUTS.OPT_SNXT = OPT_SNXT_SJPAR; }
  else {
    sprintf(c1err,"'%s' is invalid option for EXTINC_HOSTGAL",INPUTS.GENSNXT );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }


  // check HOSTLIB_GENZPHOT options
  int USE = 0 ;
  if ( INPUTS.HOSTLIB_GENZPHOT_OUTLIER[0] > 0.0 ) { USE=1; }
  if ( INPUTS.HOSTLIB_GENZPHOT_OUTLIER[1] > 0.0 ) { USE=1; }
  for(i=0; i<5; i++) { 
    if ( INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[i] >  0.0 ) { USE=1; }
    if ( INPUTS.HOSTLIB_GENZPHOT_BIAS[i]     != 0.0 ) { USE=1; }
  }
  INPUTS.USE_HOSTLIB_GENZPHOT = USE ;

  // -----------------------------------------------
  // load generic "LUMIPAR" variables based on model
  // Also load INDEX_GENMODEL index

  GENFRAME_OPT = NULLINT ;
  GENLC.ptr_SHAPEPAR = &GENLC.NOSHAPE ;  // default

  PTR_SMEAR_MODELNAME = INPUTS.GENMAG_SMEAR_MODELNAME ;
  USE_SMEAR_MODELNAME = !IGNOREFILE(PTR_SMEAR_MODELNAME);
  ISCOH_SMEAR         = ( strstr(PTR_SMEAR_MODELNAME,"COH") != NULL );

  if ( INPUTS.SIMLIB_DUMP > 0 ) {
    // set params to read entire SIMLIB once; then quit. (Mar 2021)
    INPUTS.SIMLIB_MSKOPT |= SIMLIB_MSKOPT_QUIT_NOREWIND;
    INPUTS.NGENTOT_LC = 10000; INPUTS.NGEN_LC=0;
  }

  // - - - - - -  - - - 

  if ( INDEX_GENMODEL  == MODEL_STRETCH ) {
    
    GENFRAME_OPT    = GENFRAME_REST ;
    sprintf(GENLC.SHAPEPAR_NAME,    "STRETCH");
    sprintf(GENLC.SHAPEPAR_GENNAME, "STRETCH");
    GENLC.ptr_SHAPEPAR = &GENLC.STRETCH ;
    copy_GENGAUSS_ASYM( &INPUTS.GENGAUSS_STRETCH, &INPUTS.GENGAUSS_SHAPEPAR );
  }

  else if ( INDEX_GENMODEL  == MODEL_MLCS2k2  ) {

    GENFRAME_OPT    = GENFRAME_REST;
    sprintf(GENLC.DISTANCE_NAME,    "DLMAG" );
    sprintf(GENLC.SHAPEPAR_NAME,    "DELTA" );
    sprintf(GENLC.SHAPEPAR_GENNAME, "DELTA" );
    sprintf(GENLC.COLORPAR_NAME,  "%s",  PARNAME_AV  );
    sprintf(GENLC.COLORPAR2_NAME, "%s",  PARNAME_RV  );
    GENLC.ptr_SHAPEPAR = &GENLC.DELTA ;
    copy_GENGAUSS_ASYM( &INPUTS.GENGAUSS_DELTA, &INPUTS.GENGAUSS_SHAPEPAR );

    if ( INPUTS.H0 == H0_SALT2 ) { INPUTS.H0 = (double)H0_MLCS; }
  }

  else if ( INDEX_GENMODEL == MODEL_SNOOPY  ) {
    
    GENFRAME_OPT    = GENFRAME_REST; 
    sprintf(GENLC.DISTANCE_NAME,    "DLMAG"    );
    sprintf(GENLC.SHAPEPAR_NAME,    "STRETCH"  );
    sprintf(GENLC.SHAPEPAR_GENNAME, "STRETCH"  );
    sprintf(GENLC.COLORPAR_NAME,  "%s",   PARNAME_AV );
    sprintf(GENLC.COLORPAR2_NAME, "%s",   PARNAME_RV );

    GENLC.ptr_SHAPEPAR = &GENLC.STRETCH ;
    copy_GENGAUSS_ASYM( &INPUTS.GENGAUSS_STRETCH, &INPUTS.GENGAUSS_SHAPEPAR );
  }

  else if ( INDEX_GENMODEL == MODEL_S11DM15  ) {

    GENFRAME_OPT    = GENFRAME_OBS; 
    sprintf(GENLC.DISTANCE_NAME,     "DLMAG" );
    sprintf(GENLC.SHAPEPAR_NAME,     "DM15"  );
    sprintf(GENLC.SHAPEPAR_GENNAME,  "DM15"  );
    sprintf(GENLC.COLORPAR_NAME,  "%s", PARNAME_AV );
    sprintf(GENLC.COLORPAR2_NAME, "%s", PARNAME_RV );
    GENLC.ptr_SHAPEPAR = &GENLC.DM15 ;
    copy_GENGAUSS_ASYM( &INPUTS.GENGAUSS_DM15, &INPUTS.GENGAUSS_SHAPEPAR );
  }

  else if ( INDEX_GENMODEL == MODEL_SALT2 ) {
   
    GENFRAME_OPT   = GENFRAME_OBS ;
    sprintf(GENLC.DISTANCE_NAME,    "mB"     );
    sprintf(GENLC.SHAPEPAR_NAME,    "x1"     );
    sprintf(GENLC.SHAPEPAR_GENNAME, "SALT2x1");
    sprintf(GENLC.COLORPAR_NAME,    "c"      );
    sprintf(GENLC.SHAPEPAR2_NAME,   "alpha"  );
    sprintf(GENLC.COLORPAR2_NAME,   "beta"   );
    GENLC.ptr_SHAPEPAR = &GENLC.SALT2x1 ;
    // if there is just one BETA bin for SNGRID, set RANGE = MEAN
    double BETA = INPUTS.GENGAUSS_SALT2BETA.PEAK ;
    int    NBIN = GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_COLORLAW] ;
    if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  && NBIN == 1 ) { 
      INPUTS.GENGAUSS_SALT2BETA.RANGE[0] = BETA ;
      INPUTS.GENGAUSS_SALT2BETA.RANGE[1] = BETA ;
    }

    copy_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2x1, &INPUTS.GENGAUSS_SHAPEPAR );

  }
  else if ( INDEX_GENMODEL == MODEL_SIMSED ) {
   
    // allow only COH model for intrinsic scatter
    if ( USE_SMEAR_MODELNAME && !ISCOH_SMEAR ) {
      sprintf(c1err,"Invalid GENMAG_SMEAR_MODELNAME: %s",PTR_SMEAR_MODELNAME);
      sprintf(c2err,"Only COH and COH([sig]) allowed with NON1ASED");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    GENFRAME_OPT   = GENFRAME_OBS;
    INPUTS.GENGAUSS_SHAPEPAR.RANGE[0] = -1.0 ; // anything non-zero
    INPUTS.GENGAUSS_SHAPEPAR.RANGE[1] = +1.0 ;
    prep_user_SIMSED(); 
  }

  else if ( IS_PySEDMODEL ) {

    GENFRAME_OPT   = GENFRAME_OBS;

  }

  else  if ( INDEX_GENMODEL  == MODEL_NON1ASED ) {    
    GENFRAME_OPT    = GENFRAME_OBS ;

    // allow only COH model for intrinsic scatter
    if ( USE_SMEAR_MODELNAME && !ISCOH_SMEAR ) {
      sprintf(c1err,"Invalid GENMAG_SMEAR_MODELNAME: %s",PTR_SMEAR_MODELNAME);
      sprintf(c2err,"Only COH and COH([sig]) allowed with NON1ASED");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if( strlen(INPUTS.MODELPATH) > 0 ) 
      { sprintf(INPUTS.NON1ASED.PATH, "%s", INPUTS.MODELPATH);   }
  }
  else if ( INDEX_GENMODEL  == MODEL_NON1AGRID ) {
    GENFRAME_OPT = GENFRAME_OBS ;
    sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"NONE"); // no SNIa scatter model
  }
  else if ( INDEX_GENMODEL == MODEL_LCLIB ) {
    GENFRAME_OPT = GENFRAME_OBS ;
    sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"NONE");    // no SNIa scatter model

    LCLIB_INFO.HOSTLIB_MSKOPT = INPUTS.HOSTLIB_MSKOPT ; 
    INPUTS.HOSTLIB_MSKOPT = INPUTS.HOSTLIB_USE = 0;
    INPUTS.GENSIGMA_VPEC = INPUTS.VPEC_ERR = INPUTS.VEL_CMBAPEX = 0.0 ;
    INPUTS.MJD_EXPLODE   = 0.0 ;
    INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE = -9 ;

    for(j=0; j<2; j++ ) {
      INPUTS.GENRANGE_REDSHIFT[j] = 0.0 ;
      INPUTS.GENRANGE_PEAKMJD[j] = 
	INPUTS.GENRANGE_MJD[0] ;      // delta func
      INPUTS.GENRANGE_TREST[j]   =
	INPUTS.GENRANGE_MJD[j]-INPUTS.GENRANGE_MJD[0];
    }

    
    if ( INPUTS.SIMLIB_DUMP > 0 ) {
      INPUTS.GENRANGE_PEAKMJD[0] = INPUTS.GENRANGE_MJD[0] ; 
      INPUTS.GENRANGE_PEAKMJD[1] = INPUTS.GENRANGE_MJD[1] ; 
    }

    INPUTS.SIMLIB_MXREPEAT = INPUTS.SIMLIB_NREPEAT ;
    
    if ( ISRATE_LCLIB == 0 ) {
      sprintf(c1err,"Must use DNDB rate model with LCLIB");
      sprintf(c2err,"%s rate-model is invalid.", INPUTS.RATEPAR.NAME);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( INPUTS.GENRANGE_MJD[0] < 20001.0 ) {
      sprintf(c1err,"Must specify GENRANGE_MJD for LCLIB model");
      sprintf(c2err,"Check sim-input file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }
  }

  else if ( INDEX_GENMODEL == MODEL_FIXMAG ) {
    
    GENFRAME_OPT  = GENFRAME_OBS;
    INPUTS.GENGAUSS_SHAPEPAR.USE      =  true ;
    INPUTS.GENGAUSS_SHAPEPAR.PEAK     =  0.0 ;
    INPUTS.GENGAUSS_SHAPEPAR.RANGE[0] = -0.01 ;  // made up numbers
    INPUTS.GENGAUSS_SHAPEPAR.RANGE[1] = +0.01 ;  // 
    INPUTS.GENGAUSS_SHAPEPAR.SIGMA[0] =  0.0 ;
    INPUTS.GENGAUSS_SHAPEPAR.SIGMA[1] =  0.0 ;

    // turn off Galactic extinction and intrinsic smearing
    INPUTS.MWEBV_FLAG     = 0  ;
    INPUTS.MWEBV_SIGRATIO = 0.0 ;
    INPUTS.MWEBV_SIG      = 0.0 ;
    INPUTS.MWEBV_SHIFT    = 0.0 ;
    sprintf(INPUTS.GENMAG_SMEAR_MODELNAME, "NONE") ;


    if ( INPUTS.GENFRAME_FIXMAG == GENFRAME_OBS ) {
      // set redshift to 0 for obs-frame to avoid adding DLMAG
      INPUTS.GENRANGE_REDSHIFT[0] = 0.0 ;
      INPUTS.GENRANGE_REDSHIFT[1] = 0.0 ;
      INPUTS.GENSIGMA_REDSHIFT    = 0.0 ; // Apr 20 2018
    }
    else {
      // abort if redshift range not specified
      if ( INPUTS.GENRANGE_REDSHIFT[1] == 0.0 ) {
	sprintf(c1err,"Must specify GENRANGE_REDSHIFT");
	sprintf(c2err,"for GENMODEL: %s", INPUTS.GENMODEL);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }

    // turn off all mag offsets
    INPUTS.GENMAG_OFF_GLOBAL = 0.0 ;
    INPUTS.GENMAG_OFF_NON1A  = 0.0 ;
    for ( ifilt  = 0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs  = INPUTS.IFILTMAP_OBS[ifilt];
      INPUTS.GENMAG_OFF_MODEL[ifilt_obs] = 0.0 ;
      INPUTS.GENMAG_OFF_ZP[ifilt_obs]    = 0.0 ;
    }
  }

  else if ( INDEX_GENMODEL == MODEL_SIMLIB ) {
    // Nov 2019
    GENFRAME_OPT    = GENFRAME_OBS ; 

    char CPOLY_DZTOL[] = "0.05" ;
    parse_GENPOLY(CPOLY_DZTOL, "DZTOL", &INPUTS.HOSTLIB_GENPOLY_DZTOL,fnam);

    GENLC.ptr_SHAPEPAR = &GENLC.SALT2x1 ; 
    sprintf(INPUTS_SEARCHEFF.USER_zHOST_FILE, "NONE" );
    sprintf(INPUTS_SEARCHEFF.USER_SPEC_FILE,  "NONE" );

  }
  else {
    sprintf(c1err,"%s is not a valid genmag-model", INPUTS.MODELNAME);
    sprintf(c2err,"Check GENMODEL keyword in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // construct mapping between sparse filter index and obs-filter index.

  GENLC.NFILTDEF_OBS = INPUTS.NFILTDEF_OBS ;

  for ( ifilt  = 0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs  = INPUTS.IFILTMAP_OBS[ifilt];
    GENLC.IFILTMAP_OBS[ifilt]         = ifilt_obs ;  
    GENLC.IFILTINVMAP_OBS[ifilt_obs]  = ifilt ;
  }


  // make sure to turn off SIMSED params if not SIMSED model.
  // This allows leaving SIMSED parameter defs in the sim-input
  // file when switching to other models.
  if ( INDEX_GENMODEL != MODEL_SIMSED ) { INPUTS.NPAR_SIMSED = 0 ; }

  // ===========================================
  
  // check options to work only for NON1A (with sim_SNmix)

  if ( INPUTS.NON1A_MODELFLAG > 0 || INDEX_GENMODEL == MODEL_SIMSED ) { 
      INPUTS.NGEN_SCALE  *= INPUTS.NGEN_SCALE_NON1A ;  

      if ( INPUTS.GENMAG_OFF_NON1A != 0.0 )  
	{  GENLC.GENMAG_OFF_GLOBAL = INPUTS.GENMAG_OFF_NON1A ; }
  }


  // check MJD_EXPLODE option only for NON1A or SIMSED
  if ( INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE >=0 ) {
    int OK=0 ;
    if ( INDEX_GENMODEL == MODEL_SIMSED   ) { OK=1; }
    if ( INDEX_GENMODEL == MODEL_NON1ASED ) { OK=1; }
    if ( OK==0 ) {
      sprintf(c1err,"MJD_EXPLODE = %.3f invalid for GENMODEL=%s",
	      INPUTS.MJD_EXPLODE, INPUTS.GENMODEL );
      sprintf(c2err,
	      "Define Explosion time for NON1ASED or SIMSED models.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }


  INPUTS.NGEN_LC    *= INPUTS.NGEN_SCALE ;
  INPUTS.NGENTOT_LC *= INPUTS.NGEN_SCALE ;

  if ( INPUTS.NGEN_LC > 0 ) 
    { NGEN = INPUTS.NGEN_LC ; }
  else if ( INPUTS.NGENTOT_LC > 0 ) 
    { NGEN = INPUTS.NGENTOT_LC ; }
  else if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID  ) {
    NGEN = -9;
    sprintf(c1err,"NGEN_LC=0 & NGENTOT_LC=0" );
    errmsg(SEV_WARN, 0, fnam, c1err, ""); 
  }
  else if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
    NGEN = 1000 ;  // just to set screen-update; will be set later
  }

  if ( INPUTS.NGEN_LC > 0 && INPUTS.NGENTOT_LC > 0 ) {
    sprintf(c1err,"NGEN_LC = %d and NGENTOT_LC = %d :", 
	    INPUTS.NGEN_LC, INPUTS.NGENTOT_LC );
    sprintf(c2err,"Only one of the above can be non-zero.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  XGEN_F = (float)NGEN ;
  NGEN  = (int)(XGEN_F+0.5) ;

  INPUTS.NGEN = NGEN;
  set_screen_update(NGEN);

  // -------- Galactic extinction ---------------
  if ( INPUTS.GENRANGE_MWEBV[0]  >= 0.0 ) { INPUTS.MWEBV_FLAG = 2; }
  
  // get string describing Galactic reddening color law
  OPT  = INPUTS.OPT_MWCOLORLAW ;
  text_MWoption("COLORLAW", OPT, INPUTS.STR_MWCOLORLAW ); // return STR
  if ( OPT == 0 ) { INPUTS.MWEBV_FLAG = 0; } // turn off

  // get string describing option to modify MWEBV_SFD
  OPT  = INPUTS.OPT_MWEBV ;
  text_MWoption(PARNAME_EBV, OPT, INPUTS.STR_MWEBV ); // return STR
  if ( OPT == 0 ) { INPUTS.MWEBV_FLAG = 0; } // turn off
  
  // --------------------------------------------------
  // check exposure time for each filter
  for ( ifilt=0; ifilt < MXFILTINDX ; ifilt++ ) {
    tmp_F = INPUTS.EXPOSURE_TIME_FILTER[ifilt] ;
    if ( tmp_F != 1.0 ) { 
      GENLC.SCALE_NOISE[ifilt]       = sqrtf ( tmp_F ); 
      GENLC.SHIFT_ZPTSIMLIB[ifilt]   = 2.5 * log10f ( tmp_F );
    }

    // set global logical if any band uses MJD_TEMPLATE (Sep 2017)
    if ( INPUTS.MJD_TEMPLATE_FILTER[ifilt] > 1.0 ) 
      { INPUTS.USE_MJD_TEMPLATE = 1 ; }

    // Jun 2021: set lamshift per filter, sparse indices
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;  
    INPUTS.TMPOFF_LAMFILT[ifilt] = INPUTS.FUDGESHIFT_LAM_FILTER[ifilt_obs];
  }
	 

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENRANDOM ) {
    if ( IGNOREFILE(INPUTS.SIMLIB_FILE) ) {
      sprintf(c1err,"SIMLIB_FILE not specified." );
      sprintf(c2err,"Check %s/simlib/  area", PATH_SNDATA_ROOT );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  prep_user_CUTWIN(); // 9.25.2017


  // turn off all sim-ouput if we are just going to dump
  // stuff to the screen and then quit

  if ( INPUTS.NVAR_SIMGEN_DUMP == 0  ||
       INPUTS.SIMLIB_DUMP  >= 0        ) 
    { INPUTS.FORMAT_MASK = 0;  DOCHECK_FORMAT_MASK=0; }



  // 9.28.2020: find kcor file and update INPUTS.KCOR_FILE if needed
  char PATH_KCOR_LIST[2*MXPATHLEN], kcorFile[MXPATHLEN];
  sprintf(kcorFile, "%s", INPUTS.KCOR_FILE);
  sprintf(PATH_KCOR_LIST, "%s %s/kcor",  PATH_USER_INPUT, PATH_SNDATA_ROOT );
  find_pathfile(kcorFile, PATH_KCOR_LIST, INPUTS.KCOR_FILE, fnam ); 

  prep_user_cosmology(); // Oct 16 2020 

  // --------------------------------------
  //----------- PRINT SUMMARY -------------
  // --------------------------------------

  if ( INPUTS.DASHBOARD_DUMPFLAG ) { return ; }

  print_banner(" prep_user_inputs summary ");

  printf("\t SIMLIB file        : %s   (start LIBID=%d)\n", 
	 INPUTS.SIMLIB_FILE, INPUTS.SIMLIB_IDSTART );

  printf("\t Generation Version : %s \n", INPUTS.GENVERSION );
  printf("\t Generation source  : %s \n", INPUTS.GENSOURCE  );
  printf("\t Generation model   : %s \n", INPUTS.MODELNAME  );

  printf("\t Number of LC to Generate: %d \n", INPUTS.NGEN_LC );

  if ( INPUTS.NGEN_LC >= MXCID_SIM ) {
    sprintf(c1err,"NGEN_LC = %d exceeds max allowed.", INPUTS.NGEN_LC);
    sprintf(c2err,"Do you really need that many SNe ?");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  
  if ( IGNOREFILE(INPUTS.HzFUN_FILE) ) {
    printf("\t OMEGA_(MATTER,LAMBDA)= %5.3f, %5.3f,    w0= %5.2f   H0=%5.1f \n"
	   ,INPUTS.OMEGA_MATTER, INPUTS.OMEGA_LAMBDA
	   ,INPUTS.w0_LAMBDA, INPUTS.H0  );
  }
  else {
    printf("\t HzFUN from file: %s \n", INPUTS.HzFUN_FILE);
  }

  printf("\t KCOR  file : %s \n", INPUTS.KCOR_FILE );


  printf("\t Observer Gen-FILTERS  :  %s ", INPUTS.GENFILTERS );

  printf(" \n" );

  printf("\t Random number seed: %d  (NSTREAM=%d)\n", 
	 INPUTS.ISEED, INPUTS.NSTREAM_RAN );

  printf("\t Gen-Range for RA(deg)  : %8.3f to %8.3f \n", 
	 INPUTS.GENRANGE_RA[0], INPUTS.GENRANGE_RA[1] );
  printf("\t Gen-Range for DEC(deg): %8.3f to %8.3f \n", 
	 INPUTS.GENRANGE_DEC[0], INPUTS.GENRANGE_DEC[1] );

  printf("\t Gen-Range for ZCMB : %6.3f to %6.3f  "
	 "(sigma=%7.4f, bias=%.5f) \n"
	 , INPUTS.GENRANGE_REDSHIFT[0], INPUTS.GENRANGE_REDSHIFT[1] 
	 , INPUTS.GENSIGMA_REDSHIFT,    INPUTS.GENBIAS_REDSHIFT );

  if ( INPUTS.GENRANGE_REDSHIFT[1] > ZMAX_SNANA ) {
    sprintf(c1err,"GENRANGE_REDSHIFT: %.3f %.3f",
	    INPUTS.GENRANGE_REDSHIFT[0], INPUTS.GENRANGE_REDSHIFT[1]);
    sprintf(c2err,"extends beyond ZMAX_SNANA= %.3f", ZMAX_SNANA);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\t Gen-Range for PEAKMJD  : %8.1f to %8.1f  \n", 
	 INPUTS.GENRANGE_PEAKMJD[0], INPUTS.GENRANGE_PEAKMJD[1] );

  printf("\t Gen-Range for Trest    : %8.1f to %8.1f  days \n", 
	 INPUTS.GENRANGE_TREST[0], INPUTS.GENRANGE_TREST[1] );

  if ( INPUTS.GENRANGE_TREST[0] >= INPUTS.GENRANGE_TREST[1] ) {
    sprintf(c1err,"Invalid GENRANGE_TREST: %f %f", 
	    INPUTS.GENRANGE_TREST[0], INPUTS.GENRANGE_TREST[1] );
    sprintf(c2err,"Check sim-input key");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  printf("\t Gen-Range for SHAPEPAR  : %8.1f to %8.1f  \n", 
	 INPUTS.GENGAUSS_SHAPEPAR.RANGE[0], 
	 INPUTS.GENGAUSS_SHAPEPAR.RANGE[1] );

  printf("\t Gen-Range for AV  : %4.2f to %4.2f  (dN/dAv = exp(-AV/%4.2f) \n", 
	 INPUTS.GENRANGE_AV[0], INPUTS.GENRANGE_AV[1], INPUTS.GENEXPTAU_AV );

  printf("\t Gen-Mean  for RV  : %4.2f  \n", INPUTS.GENGAUSS_RV.PEAK );
  printf("\t Gen-sigma for RV  : %4.2f , %4.2f (lower , upper ) \n",
	 INPUTS.GENGAUSS_RV.SIGMA[0], INPUTS.GENGAUSS_RV.SIGMA[1] );
  printf("\t Gen-Range for RV  : %4.2f to %4.2f \n",
	 INPUTS.GENGAUSS_RV.RANGE[0], INPUTS.GENGAUSS_RV.RANGE[1] );


  prep_genmag_offsets(); 

  prep_GENPDF_FLAT();    // Jul 2020: GENPDF

  // check for required input;
  // Allow exceptions for GRID option and INIT_ONLY

  if ( !INPUTS.INIT_ONLY && GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID ) {
    sprintf(c2err,"Check input file: %s ", input_file );

    IBLANK = strcmp(INPUTS.GENVERSION,"BLANK");
    if ( IBLANK == 0) {
      sprintf(c1err,"Must fill GENVERSION: argument in input file");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    IBLANK = strcmp(INPUTS.GENSOURCE,"BLANK");
    if ( IBLANK == 0) {
      sprintf(c1err,"GENSOURCE: must either be RANDOM or a DATA-VERSION");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }  // end of NOT-GRID if-block



  if ( INPUTS.GENSIGMA_REDSHIFT < 0.0 && INPUTS.HOSTLIB_USE == 0 ) {
    sprintf(c1err,"GENSIGMA_REDSHIFT = %f => use HOST-PHOTOZ", 
	    INPUTS.GENSIGMA_REDSHIFT);
    sprintf(c2err,"but HOSTLIB not selected.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }



  // for DMPTREST option: turn off
  // mag-smearing, redshift-smearing, HEL->CMB
  if ( INPUTS.GENRANGE_DMPTREST[0] < 9000.0 ) {
    INPUTS.USEFLAG_DMPTREST = 1;

    INPUTS.NFILT_SMEAR        = 0 ;      // no intrinsic scatter
    INPUTS.GENMAG_SMEAR[0]    = 0.0;
    INPUTS.GENMAG_SMEAR[1]    = 0.0;
    sprintf(INPUTS.GENMAG_SMEAR_MODELNAME, "NONE");

    INPUTS.GENSIGMA_REDSHIFT = 0.0 ; // no redshift smearing
    INPUTS.GENBIAS_REDSHIFT  = 0.0 ;
    INPUTS.VEL_CMBAPEX       = 0.0 ;

    INPUTS.HOSTLIB_USE = INPUTS.HOSTLIB_MSKOPT = 0 ;
    sprintf(INPUTS.HOSTLIB_FILE,  "NONE" );  // no HOSTLIB
  }

  // init all of the WRFLAGs to zero
  WRFLAG_VBOSE     = 0 ;
  WRFLAG_TEXT      = 0 ;
  WRFLAG_MODEL     = 0 ;
  WRFLAG_BLINDTEST = 0 ;
  WRFLAG_CIDRAN    = 0 ;
  WRFLAG_FITS      = 0 ;
  WRFLAG_FILTERS   = 0 ;
  WRFLAG_COMPACT   = 0 ;

  // check for whether to write FULL, TERSE, FITS, etc ,
  // EXCEPT for the GRID-GEN option (for psnid ...), 
  if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID  ) {
    WRFLAG_VBOSE     = ( INPUTS.FORMAT_MASK  & WRMASK_VBOSE     ) ;
    WRFLAG_TEXT      = ( INPUTS.FORMAT_MASK  & WRMASK_TEXT      ) ;
    WRFLAG_MODEL     = ( INPUTS.FORMAT_MASK  & WRMASK_MODEL     ) ;
    WRFLAG_BLINDTEST = ( INPUTS.FORMAT_MASK  & WRMASK_BLINDTEST ) ;
    WRFLAG_CIDRAN    = ( INPUTS.FORMAT_MASK  & WRMASK_CIDRAN    ) ;
    WRFLAG_FITS      = ( INPUTS.FORMAT_MASK  & WRMASK_FITS      ) ;
    WRFLAG_FILTERS   = ( INPUTS.FORMAT_MASK  & WRMASK_FILTERS   ) ;
    WRFLAG_COMPACT   = ( INPUTS.FORMAT_MASK  & WRMASK_COMPACT   ) ;
  }
  if ( WRFLAG_BLINDTEST ) { INPUTS.WRITE_MASK  = WRITE_MASK_LCMERGE ; }
  if ( WRFLAG_COMPACT   ) { INPUTS.WRITE_MASK += WRITE_MASK_COMPACT ; }
  if ( INPUTS.MAGMONITOR_SNR) { 
    SNDATA.MAGMONITOR_SNR = INPUTS.MAGMONITOR_SNR ;
    sprintf(SNDATA.VARNAME_SNRMON, "SIM_SNRMAG%2.2d", SNDATA.MAGMONITOR_SNR);
    INPUTS.WRITE_MASK    += WRITE_MASK_SIM_SNRMON;
  }

  if ( INPUTS.WRFLAG_MODELPAR ) 
    { INPUTS.WRITE_MASK += WRITE_MASK_SIM_MODELPAR; }

  // abort if no valid format is given
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENRANDOM  && 
       DOCHECK_FORMAT_MASK  ) {

    if ( WRFLAG_VBOSE ) {
      sprintf(c1err,"No longer support FORMAT_MASK += 1 (legacy VERBOSE)");
      sprintf(c2err,"Check manual for valid FORMAT_MASK");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( !WRFLAG_TEXT  && !WRFLAG_FITS  ) {
      sprintf(c1err,"Invalid FORMAT_MASK=%d. Check manual. ", 
	      INPUTS.FORMAT_MASK) ;
      sprintf(c2err,"Must specify TEXT (%d)  or FITS (%d)", 
	      WRMASK_TEXT, WRMASK_FITS );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG ) {
    INPUTS.SMEARFLAG_HOSTGAL += SMEARMASK_HOSTGAL_PHOT; // internal logical flag
  }

  if ( !IGNOREFILE(INPUTS.HOSTNOISE_FILE)  ) {
    INPUTS.SMEARFLAG_HOSTGAL += SMEARMASK_HOSTGAL_IMAGE; 
  }

  // if we are fudging SNRMAX, then turn off zero-point smearing
  // and host-gal noise
  if ( INPUTS.OPT_FUDGE_SNRMAX > 0  ) {
    INPUTS.SMEARFLAG_ZEROPT    = 0 ;
    INPUTS.SMEARFLAG_HOSTGAL   = 0 ;
    INPUTS.APPLY_SEARCHEFF_OPT = 0 ;    // Feb 2012
    INPUTS.MWEBV_SIGRATIO      = 0.0 ;  // Feb 2012

    // turn off GALMAG bit-option
    if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG ) 
      { INPUTS.HOSTLIB_MSKOPT -= HOSTLIB_MSKOPT_GALMAG ; }
  }

  // Feb 2012:
  // make sure smear-INTERP filters (SMEAR = -1) 
  // are at the end of the list.
  float SMEAR_F ;
  int INTERP ;
  INTERP = 0;
  for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;  
    SMEAR_F   = INPUTS.GENMAG_SMEAR_FILTER[ifilt_obs] ;
    
    if ( SMEAR_F < -0.001 ) { INTERP = 1; }

    if ( INTERP  && (SMEAR_F > 1.0E-6) ) {
      sprintf(c1err,"SMEAR-filter '%c' must be listed before ",
	      FILTERSTRING[ifilt_obs] );
      sprintf(c2err,"SMEAR-INTERP filters.") ;
      errmsg(SEV_WARN, 0, fnam, c1err, c2err ); 
    }

    //    printf(" xxx SMEAR(%c) = %6.3f \n", FILTERSTRING[ifilt_obs], SMEAR );
  }


  // init optional scatter matrix
  INIT_COVMAT_SCATTER();

  // init optional fudge to get specific SNRMAX
  INIT_FUDGE_SNRMAX();

  // init non-linearity function (Apr 2016)
  INIT_NONLIN(INPUTS.NONLINEARITY_FILE);

  // pass trigger epoch duration to sntools_trigger.
  INPUTS_SEARCHEFF.TIME_SINGLE_DETECT = INPUTS.NEWMJD_DIF ;

  // if there are no PEC1A SEDs, make sure PEC1A rate is off
  NTMP = INPUTS.RATEPAR_PEC1A.NMODEL_ZRANGE ;
  if ( INPUTS.NON1ASED.NPEC1A == 0  && NTMP>0 ) { 
    printf("\t Found no PEC1A SEDs --> turn off DNDZ_PEC1A rate model\n");
    INPUTS.RATEPAR_PEC1A.NMODEL_ZRANGE = 0 ; 
  }
  

  if ( ISRATE_LCLIB &&  INDEX_GENMODEL != MODEL_LCLIB ) { 
    sprintf(c1err,"Must use DNDB rate-model with GENMODEL=LCLIB.");
    sprintf(c2err,"GENMODEL = %s is not a valid Galactic model.", 
	    INPUTS.GENMODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  } 


  // Feb 2018:
  // if using FLUXERRMODEL, make sure none of the legacy options are used
  if ( !IGNOREFILE(INPUTS.FLUXERRMODEL_FILE)  ) {
    if ( !IGNOREFILE(INPUTS.HOSTNOISE_FILE)  ) {
      sprintf(c1err,"Cannot mix FLUXERRMODEL_FILE with HOSTNOISE_FILE");
      sprintf(c2err,"Pick one, not both.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( INPUTS.FUDGEOPT_FLUXERR  ) {
      sprintf(c1err,"Cannot mix FLUXERRMODEL_FILE with FUDGEOPT_FLUXERR");
      sprintf(c2err,"Pick one, not both.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    
    /* NOT_YET ...
    // Aug 26 2020: check for option to require DOCANA
    if ( INPUTS.REQUIRE_DOCANA ) 
      { INPUTS.FLUXERRMODEL_OPTMASK += MASK_REQUIRE_DOCANA_FLUXERRMAP; }
    */
  }

  if ( INPUTS.RESTORE_DES3YR ) {
    INPUTS.RESTORE_HOSTLIB_BUGS = true ;
    INPUTS.RESTORE_FLUXERR_BUGS = true ;
    printf("\t Restore bugs for DES3YR analysis.\n");
  }
    

  printf("\n");

  return ;

} // end of prep_user_input().


// ===================================
void prep_user_cosmology(void) {

  // Created Oct 2020
  // Call init_HzFUN_INFO to either store user-input cosmology params,
  // or to read z,H(z) from 2-column input file.
  //
 
  double cosPar[NCOSPAR_HzFUN];
  char  *HzFUN_FILE = INPUTS.HzFUN_FILE ;
  int ipar;
  char fnam[] = "prep_user_cosmology";

  // ------- BEGIN ----------

  cosPar[ICOSPAR_HzFUN_H0] = INPUTS.H0;
  cosPar[ICOSPAR_HzFUN_OM] = INPUTS.OMEGA_MATTER ;
  cosPar[ICOSPAR_HzFUN_OL] = INPUTS.OMEGA_LAMBDA ;
  cosPar[ICOSPAR_HzFUN_w0] = INPUTS.w0_LAMBDA ;
  cosPar[ICOSPAR_HzFUN_wa] = INPUTS.wa_LAMBDA ;
  
  // init structure that gets passed later to cosmology functions
  int VBOSE = 1;
  init_HzFUN_INFO(VBOSE, cosPar, HzFUN_FILE, 
		  &INPUTS.HzFUN_INFO ); // <== returned 

  return;

} // end prep_user_cosmology

// ===================================
void prep_user_CUTWIN(void) {

  int   MXFILT_CUT, NFILT_CUT, i ;
  char  CFILT[MXFILTINDX], CLOGIC[20] ;
  char  fnam[] = "prep_user_CUTWIN" ;

  // --------------- BEGIN ----------------

  // Parse CUTWIN_SNRMAX_LOGIC (only if cuts are requested)
  for ( i=1 ; i <= INPUTS.NCUTWIN_SNRMAX ;  i++ ) {
    sprintf(CLOGIC, "%s", INPUTS.CUTWIN_SNRMAX_LOGIC[i] );
    sprintf(CFILT,  "%s", INPUTS.CUTWIN_SNRMAX_FILTERS[i] );
    MXFILT_CUT =  strlen(CFILT) ;   

    if ( strcmp(CLOGIC,"AND") == 0 )
      { NFILT_CUT = MXFILT_CUT ; }
    else if ( strcmp(CLOGIC,"OR") == 0 ) 
      { NFILT_CUT = 1 ; }
    else {
      sscanf(CLOGIC , "%d", &NFILT_CUT ); 
      if ( NFILT_CUT < 1 || NFILT_CUT > MXFILT_CUT ) {
	sprintf(c1err,"Invalid NFILT_CUT=%d for FILTERS='%s'",
		NFILT_CUT, CFILT);
	sprintf(c2err,"Check CUTWIN_SNRMAX keywords in input file");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }
    // load structure
    INPUTS.CUTWIN_SNRMAX_NFILT[i] = NFILT_CUT ;
  }


  // Parse EPOCHS_FILTERS to select fast transients
  int NFILT = strlen(INPUTS.CUTWIN_EPOCHS_FILTERS);
  int IFILT_OBS ;
  char cfilt[2];
  INPUTS.CUTWIN_EPOCHS_NFILT = NFILT;
  if ( NFILT > 0 ) {
    sprintf(CFILT,"%s", INPUTS.CUTWIN_EPOCHS_FILTERS);
    for(i=0; i < NFILT; i++ ) {
      sprintf(cfilt, "%c", CFILT[i] ) ;
      IFILT_OBS = INTFILTER(cfilt);
      INPUTS.CUTWIN_EPOCHS_IFILTLIST[i] = IFILT_OBS;
    }
  } // end NFILT>0 if block


  return ;

} // end prep_user_CUTWIN


// ******************************
void prep_user_SIMSED(void) {

  // Created Feb 26 2018
  // Prepare Cholesky decomp for correlations;
  // otherwise do nothing.

  int NPAIR = INPUTS.NPAIR_SIMSED_COV; 
  int NPAR  = INPUTS.NPAR_SIMSED ;
  int ipair, j, ipartmp, ipar, IPAR0, IPAR1;
  int NTMP, NIPAR[MXPAR_SIMSED] ;
  double COV;
  char *NAME ;
  char fnam[] = "prep_user_SIMSED";

  // ------------- BEGIN ---------------

  
  if ( NPAIR == 0 ) { return; }

  sprintf(BANNER, "%s: prepare SIMSED Covariances.", fnam);
  print_banner(BANNER);

  if ( NPAIR >= MXPAR_SIMSED ) {
    sprintf ( c1err, "NPAIR(COV)=%d exceeds bound of %d",
	      NPAIR, MXPAR_SIMSED);
    sprintf( c2err, "Check SIMSED_REDCOR and SIMSED_COV keys");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  
  // idiot check on cov terms; make sure each off-diagonal
  // is filled correctly
  for(ipair=0; ipair<MXPAR_SIMSED; ipair++ ) { NIPAR[ipair]=0; }
  for(ipair=0; ipair < NPAIR; ipair++ ) {
    for(j=0; j < 2; j++ ) {
      ipartmp = INPUTS.IPARPAIR_SIMSED_COV[ipair][j] ;
      NIPAR[ipartmp]++ ;
    }
  }
  // make sure that all NIPAR are the same or zero.
  NTMP = 0 ;
  for(ipar=0; ipar < NPAR; ipar++ ) {
    if ( NIPAR[ipar] == 0 ) { continue; }
    if ( NTMP == 0  ){ NTMP = NIPAR[ipar]; }
    if ( NIPAR[ipar] != NTMP ) {
      sprintf( c1err, "Invalid redCorr matrix.");
      sprintf( c2err, "Incorrect number of off-diagonal %s-elements",
	       INPUTS.PARNAME_SIMSED[ipar] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - 
  // determine list of IPAR for COV matrix.

  INPUTS.NROW_SIMSED_COV=0;
  int irow, irow1, MATCH, NROW=0;

  for(ipair=0; ipair < NPAIR; ipair++ ) {
    for(j=0; j < 2; j++ ) {

      ipartmp = INPUTS.IPARPAIR_SIMSED_COV[ipair][j] ;
      MATCH   = 0 ;
      for(irow=0; irow < NROW; irow++ ) {
	if ( ipartmp == INPUTS.IPARLIST_SIMSED_COV[irow] ) { MATCH=1; }
      }
      if ( MATCH == 0 ) { 
	INPUTS.IPARLIST_SIMSED_COV[NROW]    = ipartmp;
	INPUTS.IROWLIST_SIMSED_COV[ipartmp] = NROW ;
	NROW++ ;
      }
    } 
  } // end loop over COV pairs 
  INPUTS.NROW_SIMSED_COV = NROW;

  //  printf("\t Size of SIMSED COV matrix: %d x %d \n", NROW, NROW);

  // - - - - - - - - - - - - - - - - - - - - - - - - - 
  int MEMD, NMAT, ISPAIR0, ISPAIR1 ;
  double SIG0, SIG1, RHO, COV_OFFDIAG ;

  MEMD = NROW*NROW * sizeof(double) ;
  INPUTS.SIMSED_DECOMP.MATSIZE  = NROW;
  INPUTS.SIMSED_DECOMP.COVMAT1D = (double *) malloc ( MEMD);

  // load COVMAT1D
  NMAT = 0 ;
  for(irow=0; irow < NROW; irow++ ) {
    IPAR0 = INPUTS.IPARLIST_SIMSED_COV[irow] ;
    SIG0  = INPUTS.GENGAUSS_SIMSED[IPAR0].SIGMA[0];
    NAME  = INPUTS.PARNAME_SIMSED[IPAR0] ;
    printf("    SIGMA(%10.10s)=%8.4f : ", NAME, SIG0 );

    for(irow1=0; irow1 < NROW; irow1++ ) {
      IPAR1 = INPUTS.IPARLIST_SIMSED_COV[irow1] ;
      SIG1  = INPUTS.GENGAUSS_SIMSED[IPAR1].SIGMA[0];
      COV_OFFDIAG = 0.0 ;

      for(ipair=0; ipair < NPAIR; ipair++ ) {	
	ISPAIR0 = (IPAR0 == INPUTS.IPARPAIR_SIMSED_COV[ipair][0] );
	ISPAIR1 = (IPAR1 == INPUTS.IPARPAIR_SIMSED_COV[ipair][1] );
	if ( ISPAIR0 && ISPAIR1 ) 
	  { COV_OFFDIAG = INPUTS.COVPAIR_SIMSED_COV[ipair]; }

	ISPAIR0 = (IPAR1 == INPUTS.IPARPAIR_SIMSED_COV[ipair][0] );
	ISPAIR1 = (IPAR0 == INPUTS.IPARPAIR_SIMSED_COV[ipair][1] );
	if ( ISPAIR0 && ISPAIR1 ) 
	  { COV_OFFDIAG = INPUTS.COVPAIR_SIMSED_COV[ipair]; }
      }

      if ( IPAR0 == IPAR1 ) 
	{ COV = SIG0*SIG0 ; }    // diagonal
      else 
	{ COV = COV_OFFDIAG ; }  // off-diag

      INPUTS.SIMSED_DECOMP.COVMAT1D[NMAT] = COV;
      RHO = COV / ( SIG0*SIG1 );
      printf("%8.4f ", RHO );

      NMAT++ ;
    }
    printf("\n"); fflush(stdout);
  }


  // -------------------------------------------
  // store Cholesky decomp

  init_Cholesky(+1, &INPUTS.SIMSED_DECOMP );
  
  printf("\n");
  //  debugexit(fnam);

  return ;
	
} // end prep_user_SIMSED

// ******************************
void prep_dustFlags(void) {

  // Created Jun 12 2020
  // Set INPUTS.DOGEN_AV for host-galaxy dust.
  // Beware to call this function after calling init_genPDF().
  //
  // Set INPUTS.DOGEN_AV = 1 for analytic EXP+HALFGAUSS function
  // Set INPUTS.DOGEN_AV = 2 for AV  map in GENPDF_FILE.
  // Set INPUTS.DOGEN_AV = 4 for EBV map in GENPDF_FILE.
  //
  // Abort if profile is defined twice (analytic and map),
  // or if only one of AV or RV profile is set.

  int  DO_WV07=0, DO_GRID=0, DO_RV=0, DO_AV=0 ;
  bool LOGPARAM;
  char fnam[] = "prep_dustFlags" ;

  // ------------ BEGIN --------------

  if ( INPUTS.GENGAUSS_RV.USE        ) { DO_RV  = 1; }
  if ( IDMAP_GENPDF(PARNAME_RV,&LOGPARAM) >= 0 ) { DO_RV += 2; }

  // check for WV07 option
  if ( INPUTS.WV07_REWGT_EXPAV > -1.0E-9 ) 
    { INPUTS.WV07_GENAV_FLAG = DO_WV07 = 1; } 


  // setUseFlag_GEN_EXP_HALFGAUSS(&INPUTS.GENPROFILE_AV,       PARNAME_AV);
  // setUseFlag_GEN_EXP_HALFGAUSS(&INPUTS.GENPROFILE_EBV_HOST, PARNAME_EBV);

  // if AV and EBV_HOST useflags are both set, then abort
  if (INPUTS.GENPROFILE_AV.USE && INPUTS.GENPROFILE_EBV_HOST.USE) {
    sprintf(c1err,"Not allowed to generate both AV and EBV_HOST");
    sprintf(c2err,"Check input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if (GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID ) {
    char *gridName = SNGRID_WRITE.NAME[IPAR_GRIDGEN_COLORPAR] ;
    DO_GRID = (strcmp(gridName,PARNAME_AV)==0) ;
  }


  // check for AV function. 
  // for AV in genPDF map.
  if ( INPUTS.GENPROFILE_AV.USE       ) { DO_AV = 1; }
  if ( INPUTS.GENPROFILE_EBV_HOST.USE ) { DO_AV = 1; }
  if ( DO_WV07  || DO_GRID            ) { DO_AV = 1; }
  if ( IDMAP_GENPDF(PARNAME_AV, &LOGPARAM ) >= 0 ) { DO_AV +=2; }
  if ( IDMAP_GENPDF(PARNAME_EBV,&LOGPARAM ) >= 0 ) { DO_AV +=4; }
  INPUTS.DOGEN_AV = DO_AV ; // store global for gen_modelPar_dust()

  // make sure that AV and RV are each defined once and only once.
  if ( DO_AV == 3 || DO_AV == 5 ) {
    sprintf(c1err,"AV/EBV profile defined twice (expFun and GENPDF_FILE)");
    sprintf(c2err,"Only one AV profile allowed.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( DO_RV == 3 ) {
    sprintf(c1err,"RV profile defined twice (expFun and GENPDF_FILE)");
    sprintf(c2err,"Only one RV profile allowed.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // next, make sure that AV and RV are both defined, or that
  // neither are defined.
  if ( ( DO_RV && !DO_AV ) || ( !DO_RV && DO_AV ) ) {
    sprintf(c1err,"DO_AV=%d and DO_RV=%d", DO_AV, DO_RV);
    sprintf(c2err,"Must specify both AV & RV, or specify neither.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return;

} // end prep_dustFlags


// *******************************************
void  prep_GENPDF_FLAT(void) {

  // Jul 8 2020
  // if GENPDF_FLAT option is set, split string to examine
  // each variable, and set approproate GENGAUSS params for
  // flat distribution.
  // E.g., GENPDF_FLAT = SALT2x1(-4:4),SALT2c(-0.4:0.6),RV(1:5)
  //   -> x1 PDF is flat from -4 to +4
  //   -> c PDF is flat from -0.4 to +0.6
  //   -> RV PDF is flat from 1 to 5

  int  NVAR, NDUM, ivar;
  char *ptrStringVar[MXVAR_GENPDF], stringVar[40], stringOpt[40] ;
  char *ptrRange[2], *varName;
  double PEAK, RANGE[2];
  double SIGMA_FLAT[2] = { 1.0E6, 1.0E6 }, TAU_FLAT=1.0E6 ;
  bool   USE_RANGE ;
  int  MEMC = 40*sizeof(char);
  char fnam[] = "prep_GENPDF_FLAT" ;

  // ------------ BEGIN -------------

  //  for FLAT distributions, put on the GENPDF_IGNORE list
  if ( strlen(INPUTS.GENPDF_FLAT) == 0 ) { return; }
  
  // split comma-sep GENPDF_FLAT string 
  for(ivar=0; ivar < MXVAR_GENPDF; ivar++ ) 
    { ptrStringVar[ivar] = (char*) malloc( MEMC ); }

  ptrRange[0] = (char*) malloc(MEMC);

  splitString(INPUTS.GENPDF_FLAT, COMMA, MXVAR_GENPDF, 
	      &NVAR, ptrStringVar);  // <== output

  for(ivar=0; ivar < NVAR; ivar++ ) {

    USE_RANGE = false;
    // note that stringVar = XYZ or XYZ([min]:[max])
    sprintf(stringVar, "%s", ptrStringVar[ivar]);

    // strip out optional range
    extractStringOpt(stringVar,stringOpt);
    varName = stringVar;

    // add varName to ignore list
    catVarList_with_comma(INPUTS.GENPDF_IGNORE,varName);

    // if stringOpt is set, split it by colon to get range
    if ( strlen(stringOpt) > 0 ) {
      splitString(stringOpt, COLON, MXVAR_GENPDF, 
		  &NDUM, ptrRange);  // <== output      
      if ( NDUM != 2 ) {
	sprintf(c1err,"Invalid NDUM=%d (should be 2)", NDUM);
	sprintf(c2err,"Check '%s' in GENPDF_FLAT = '%s' ",
		ptrStringVar[ivar], INPUTS.GENPDF_FLAT);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      sscanf(ptrRange[0], "%le", &RANGE[0] );
      sscanf(ptrRange[1], "%le", &RANGE[1] );
      PEAK = 0.5 * ( RANGE[0] + RANGE[1] );
      USE_RANGE = true;
    }

    if ( !USE_RANGE ) {
      sprintf(c1err,"Must provide gen-range for '%s'", varName);
      sprintf(c2err,"Check '%s' in GENPDF_FLAT = '%s' ",
	      ptrStringVar[ivar], INPUTS.GENPDF_FLAT);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    printf("\t Force flat PDF for %s \n", varName );

    if ( strcmp(varName,"SALT2x1") == 0 ) {   
      set_GENGAUSS_ASYM(PEAK, SIGMA_FLAT, RANGE, &INPUTS.GENGAUSS_SALT2x1);
      copy_GENGAUSS_ASYM(&INPUTS.GENGAUSS_SALT2x1,&INPUTS.GENGAUSS_SHAPEPAR);
    }
    else if ( strcmp(varName,"SALT2c") == 0 ) {
      set_GENGAUSS_ASYM(PEAK, SIGMA_FLAT, RANGE, &INPUTS.GENGAUSS_SALT2c);
    }
    else if ( strcmp(varName,PARNAME_RV) == 0 ) {
      set_GENGAUSS_ASYM(PEAK, SIGMA_FLAT, RANGE, &INPUTS.GENGAUSS_RV);
    }
    else if ( strcmp(varName,PARNAME_AV) == 0 ) {
      set_GEN_EXPON(TAU_FLAT, RANGE, &INPUTS.GENPROFILE_AV);
    }
    else if ( strcmp(varName,PARNAME_EBV) == 0 ) {
      set_GEN_EXPON(TAU_FLAT, RANGE, &INPUTS.GENPROFILE_EBV_HOST);
    }
    else {
      sprintf(c1err,"Unable to generate flat PDF for '%s'", varName);
      sprintf(c2err,"Check sim-input key GENPDF_FLAT");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }


  } // end ivar

  fflush(stdout);

  return ;

} // end  prep_GENPDF_FLAT


// **************************************
void  prep_RANSYSTPAR(void) {

  // Created Jun 2017
  // Prepare optional systematic offsets using random numbers.
  // This allows user to specify one set of Gaussian sigmas
  // (using RANSYSTPAR_XXX params) and then changing RANSEED
  // results in random set of systematic offsets.
  //
  // Do NOT try to sync randoms here; burn randoms only if required.
  //
  // Nov 9 2020: refactor filter-dependent RANSYSTPAR (see manual)

  int   ifilt, ifilt_obs, NSET=0; 
  int   NFILTDEF = INPUTS.NFILTDEF_OBS ;
  int   ILIST_RAN=1;
  float tmp, tmpSigma, *tmpRange, Range ;
  float SIGSCALE_MIN = -1.0E-6, SIGSCALE_MAX = 0.2 ;
  double gmin = -3.0, gmax=+3.0; // Gaussian clip params
  char cfilt[2], *wildcard ;
  char fnam[] = "prep_RANSYSTPAR" ;

  // ---------- BEGIN -----------

  if ( INPUTS.RANSYSTPAR.USE == 0 ) { return ; }

  sprintf(BANNER,"%s: Prepare Random set of Systematic Errors", fnam );
  print_banner(BANNER);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //   Start with variations that are synched among sub-samples
  //   (e.g., among separate simulation for each survey)
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  printf("\t* First Sync-Syst Random : %f \n", getRan_Flat1(ILIST_RAN) );

  // Galactic extinction
  tmpSigma = INPUTS.RANSYSTPAR.SIGSCALE_MWEBV ;
  checkval_F("SIGSCALE_MWEBV", 1, &tmpSigma, SIGSCALE_MIN, SIGSCALE_MAX);
  if ( tmpSigma != 0.0 ) {   
    NSET++; tmp = 1.0 + tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
    INPUTS.MWEBV_SCALE = tmp;
    printf("\t FUDGESCALE_MWEBV  = %.2f \n", tmp );
  }

  tmpSigma = INPUTS.RANSYSTPAR.SIGSHIFT_MWRV ;
  if ( tmpSigma != 0.0 ) { 
    NSET++; tmp = tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
    INPUTS.RV_MWCOLORLAW += tmp ;
    printf("\t RV_MWCOLORLAW  = %.3f \n", INPUTS.RV_MWCOLORLAW );
  }
  
  // Redshift P.Armstrong 2020
  tmpSigma = INPUTS.RANSYSTPAR.SIGSHIFT_REDSHIFT;
  if ( tmpSigma != 0.0 ) {
    NSET++; tmp = tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
    INPUTS.GENBIAS_REDSHIFT = tmp ;
    printf("\t GENBIAS_REDSHIFT  = %f \n", INPUTS.GENBIAS_REDSHIFT );
   }

  // - - - - - 
  // check wild card files
  wildcard = INPUTS.RANSYSTPAR.GENMODEL_WILDCARD; 
  if ( strlen(wildcard) > 0 ) 
    { pick_RANSYSTFILE_WILDCARD(wildcard,INPUTS.GENMODEL); }

  wildcard = INPUTS.RANSYSTPAR.GENPDF_FILE_WILDCARD;
  if ( strlen(wildcard) > 0 ) 
    { pick_RANSYSTFILE_WILDCARD(wildcard,INPUTS.GENPDF_FILE); }

  /* xxxxxxxx mark delete Nov 2 2021 xxxxxxxxxx
  if ( strlen(wildcard) > 0 ) {
      ENVreplace(wildcard,fnam,1);
      n_files = glob_file_list(wildcard, &genmodel_list);
      rand_num = getRan_Flat1(ILIST_RAN);
      // generate random ifile-index between 0 and n_files-1
      ifile_ran = (int)(rand_num * (double)n_files); 
      printf("\t Select GENMODEL %d of %d\n", ifile_ran, n_files);
      sprintf(INPUTS.GENMODEL, "%s", genmodel_list[ifile_ran]);
      if ( ifile_ran < 0 || ifile_ran >= n_files ) {
        sprintf(c1err,"Invalid ifile_ran = %d for GENMODEL_WILDCARD", 
		ifile_ran);
        sprintf(c2err,"Expected ifile_ran between 0 and %d", 
		n_files - 1);
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
  }

  // - - - - - 
  // GENPDF_FILE Wildcard R.Kessler, Nov 2 2021
  wildcard = INPUTS.RANSYSTPAR.GENPDF_FILE_WILDCARD;
  if ( strlen(wildcard) > 0 ) {
      ENVreplace(wildcard,fnam,1);
      n_files = glob_file_list(wildcard, &genmodel_list);
      rand_num = getRan_Flat1(ILIST_RAN);
      // generate random ifile-index between 0 and n_files-1
      ifile_ran = (int)(rand_num * (double)n_files); 
      printf("\t Select GENPDF_FILE %d of %d\n", ifile_ran, n_files);
      sprintf(INPUTS.GENPDF_FILE, "%s", genmodel_list[ifile_ran]);
      if ( ifile_ran < 0 || ifile_ran >= n_files ) {
        sprintf(c1err,"Invalid ifile_ran = %d for GENPDF_FILE_WILDCARD", 
		ifile_ran);
        sprintf(c2err,"Expected ifile_ran between 0 and %d", 
		n_files - 1);
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
  }
  xxxxxxxxxx end mark xxxxxxxxxxx */

  // cosmology params (Aug 2019)
  tmpSigma = INPUTS.RANSYSTPAR.SIGSHIFT_OMEGA_MATTER ;
  if ( tmpSigma != 0.0 ) { 
    NSET++; tmp = tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
    INPUTS.OMEGA_MATTER += tmp ;
    printf("\t OMEGA_MATTER  = %.3f \n", INPUTS.OMEGA_MATTER );
  }

  tmpRange = INPUTS.RANSYSTPAR.RANGESHIFT_OMEGA_MATTER ;
  if ( tmpRange[1] > tmpRange[0] ) { 
    NSET++; Range = tmpRange[1] - tmpRange[0];
    tmp = tmpRange[0] + Range * getRan_Flat1(ILIST_RAN);
    INPUTS.OMEGA_MATTER += tmp ;
    printf("\t OMEGA_MATTER  = %.3f \n", INPUTS.OMEGA_MATTER );
  }


  tmpSigma = INPUTS.RANSYSTPAR.SIGSHIFT_W0 ;
  if ( tmpSigma != 0.0 ) { 
    NSET++; tmp = tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
    INPUTS.w0_LAMBDA += tmp ;
    printf("\t w0_LAMBDA  = %.3f \n", INPUTS.w0_LAMBDA );
  }

  tmpRange = INPUTS.RANSYSTPAR.RANGESHIFT_W0 ;
  if ( tmpRange[1] > tmpRange[0] ) { 
    NSET++; Range = tmpRange[1] - tmpRange[0];
    tmp = tmpRange[0] + Range * getRan_Flat1(ILIST_RAN);
    INPUTS.w0_LAMBDA += tmp ;
    printf("\t w0_LAMBDA  = %.3f \n", INPUTS.w0_LAMBDA );
  }

  printf("\t* Last  Sync-Syst Random : %f "
	 "(should be same each survey)\n", getRan_Flat1(ILIST_RAN) );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //     Now the unsynched variations
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // need to unsync the randoms among surveys. Add 137*IDSURVEY
  // to ISEED and then re-init the randoms with new SEED.

  int IDUM      = GENLC.IDSURVEY ;
  int ISEED_OLD = INPUTS.ISEED ;
  int ISEED_NEW = ISEED_OLD + 137*IDUM ;
  INPUTS.ISEED  = ISEED_NEW ;
  init_random_seed(INPUTS.ISEED, INPUTS.NSTREAM_RAN);
  printf("\t* ISEED = %d --> %d \n", ISEED_OLD, ISEED_NEW );

  printf("\t* First Unsync-Syst Random : %f "
	 "(should differ each survey)\n", getRan_Flat1(ILIST_RAN) );

  // start with fluxerr fudging; SIGSCALE is sigma on fractional change
  tmpSigma = INPUTS.RANSYSTPAR.SIGSCALE_FLUXERR  ;
  checkval_F("SIGSCALE_FLUXERR", 1, &tmpSigma, SIGSCALE_MIN, SIGSCALE_MAX);
  if ( tmpSigma != 0.0 ) {   
    NSET++; tmp = 1.0 + tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
    INPUTS.FUDGESCALE_FLUXERR = tmp;
    printf("\t FUDGESCALE_FLUXERR(true&measured) = %.3f \n", tmp );
    for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) 
      { INPUTS.FUDGESCALE_FLUXERR_FILTER[ifilt] = tmp; }
  }

  tmpSigma = INPUTS.RANSYSTPAR.SIGSCALE_FLUXERR2 ;
  checkval_F("SIGSCALE_FLUXERR2", 1, &tmpSigma, SIGSCALE_MIN, SIGSCALE_MAX);
  if ( tmpSigma != 0.0 ) {   
    NSET=1; tmp = 1.0 + tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
    INPUTS.FUDGESCALE_FLUXERR2 = tmp;
    printf("\t FUDGESCALE_FLUXERR2(measured) = %.3f \n", tmp );
    for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) 
      { INPUTS.FUDGESCALE_FLUXERR2_FILTER[ifilt] = tmp; }
  }

  // ZP error
  for(ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt];
    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );
    tmpSigma = INPUTS.RANSYSTPAR.SIGSHIFT_ZP[ifilt];
    if ( tmpSigma != 0.0 ) {
      NSET++ ;  tmp = tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
      INPUTS.TMPOFF_ZP[ifilt]         = tmp;
      INPUTS.GENMAG_OFF_ZP[ifilt_obs] = tmp ;
      printf("\t ZPerr(%s) = %7.4f  (SIG=%.3f) \n", 
	     cfilt, tmp, tmpSigma );
    }
  }

  // LAMshift error
  for(ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt];
    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );
    tmpSigma = INPUTS.RANSYSTPAR.SIGSHIFT_LAMFILT[ifilt];
    if ( tmpSigma != 0 ) {
      NSET++ ;  tmp = tmpSigma * getRan_GaussClip(ILIST_RAN,gmin,gmax);
      INPUTS.FUDGESHIFT_LAM_FILTER[ifilt_obs] = tmp ;
      INPUTS.TMPOFF_LAMFILT[ifilt]            = tmp ;
      printf("\t LAMSHIFT(%s) = %6.2f A  (SIG=%.1f A)\n", 
	     cfilt, tmp, tmpSigma );
    }
  }


  printf("   %d Systematic Errors have been set. \n", NSET);

  // - - - - - - - - - - - -
  // check option to reset randoms with fixed random seed
  // so that there are no stat fluctuations between
  // GENVERSIONs with different systematics.
  int RANSEED_GEN = INPUTS.RANSYSTPAR.RANSEED_GEN;
  if ( RANSEED_GEN > 0 ) {
    printf("   Re-init randoms with RANSEED = %d\n", RANSEED_GEN ) ;
    init_random_seed(RANSEED_GEN, INPUTS.NSTREAM_RAN); 
  }

  printf("\n");

  return ;

} // end prep_RANSYSTPAR

// *********************************************
void pick_RANSYSTFILE_WILDCARD(char *wildcard, char *randomFile) {

  // Created Nov 2 2021 by R.kessler
  // For input *wildcard, get list of files and select a random file.
  // Load random file name into ouput arg *randomFile.
  // 

  int   ILIST_RAN=1;
  int   NJOBTOT = INPUTS.NJOBTOT ; // >0 for batch job
  int   JOBID   = INPUTS.JOBID   ; // for batch job
  int i, n_files, ifile;
  double rand_num ;
  char **genmodel_list ;
  char fnam[] = "pick_RANSYSTFILE_WILDCARD";
  // ------------- BEGIN ----------

  randomFile[0] = 0;

  ENVreplace(wildcard,fnam,1);
  n_files = glob_file_list(wildcard, &genmodel_list);

  // pick sequential ifile for batch job; else random file
  if ( NJOBTOT > 0 ) {
    ifile = JOBID-1;  // batch job; JOBID starts at 1; ifile starts at 0
  }
  else {
    // interactive
    rand_num = getRan_Flat1(ILIST_RAN);
    ifile    = (int)(rand_num * (double)n_files); 
  }

  printf("\t* Select GENPDF_FILE %d of %d from\n\t %s \n", 
	 ifile, n_files, wildcard);

  if ( ifile < 0 || ifile >= n_files ) {
    sprintf(c1err,"Invalid ifile = %d for GENPDF_FILE_WILDCARD", ifile );
    sprintf(c2err,"Expected ifile between 0 and %d", n_files-1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  sprintf(randomFile, "%s", genmodel_list[ifile]);  

  for(i=0; i < n_files; i++ ) { free(genmodel_list[i]); }
  free(genmodel_list);

  return;

} // end pick_RANSYSTFILE_WILDCARD 

// *******************************************
void prep_genmag_offsets(void) {

  int NFILTDEF = INPUTS.NFILTDEF_OBS ;
  int ifilt, ifilt_obs ;
  float tmp;

  char fnam[] = "prep_genmag_offsets" ;

  // -------------- BEGIN ---------------

  // write genmag ZP offsets and check for crazy values
  printf("\t ZP offsets (%s) : ", INPUTS.GENFILTERS);

  for ( ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt];
    tmp       = INPUTS.TMPOFF_ZP[ifilt] ;
    INPUTS.GENMAG_OFF_ZP[ifilt_obs] = tmp ;
    printf(" %5.3f", tmp );
  }
  printf("\n"); fflush(stdout);


  // write genmag MODEL offsets and check for crazy values
  printf("\t MODEL mag offsets (%s) : ", INPUTS.GENFILTERS);

  for ( ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt];
    tmp       = INPUTS.TMPOFF_MODEL[ifilt] ;
    INPUTS.GENMAG_OFF_MODEL[ifilt_obs] = tmp ;
    printf(" %5.2f", tmp );
    if ( tmp < -5.0 || tmp > 20.0 ) {
      sprintf(c1err,"Crazy GENMAG_OFF_MODEL[%c] = %f ", 
	      FILTERSTRING[ifilt_obs], tmp);
      errmsg(SEV_FATAL, 0, fnam, c1err, "  "); 
    }
  }
  printf("\n"); fflush(stdout); 

  return ;

} // end of prep_genmag_offsets

// ********************************
void genmag_offsets(void) {

  // Created Mar 16 2016
  // Apply filter-dependent mag-offset to each epoch,
  // and load peakmag per filter.
  //
  // Sep 18 2017: add lensing correction here instead of absorbing
  //              it into distance modulus.
  //
  // Sep 24 2017: check for template genmag
  // May 20 2020: check IMGNUM >= 0 before add lens-mag offset

  int epoch, ifilt_obs, IMGNUM ;
  double MAGOFF, genmag8 ;
  //  char fnam[] = "genmag_offsets";

  // ------------- begin -------------

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {

    ifilt_obs = GENLC.IFILT_OBS[epoch] ;

    // Aug 2012: always init USE_EPOCH[epoch]=0 (fixes IDLOCK bug)     
    GENLC.OBSFLAG_WRITE[epoch] = false ; 
    
    if ( !GENLC.DOFILT[ifilt_obs]  ) { continue ; }
    
    // mag offset is :
    MAGOFF = 
      + GENLC.GENMAG_OFF_GLOBAL             // user-defined global offset
      + INPUTS.GENMAG_OFF_MODEL[ifilt_obs]  // user-defined model offs
      - INPUTS.GENMAG_OFF_ZP[ifilt_obs]     // user-defined ZP offsets
      + GENLC.LENSDMU                       // lensing correction
      + INPUTS.MUSHIFT                      // user distance offset (Oct 2020)
      + GENLC.SALT2gammaDM                  // gamma from SN-host corr
    ;


    if ( INPUTS_STRONGLENS.USE_FLAG )  { 
      IMGNUM = GENSL.IMGNUM; 
      if ( IMGNUM>=0 ) { MAGOFF += GENSL.MAGSHIFT_LIST[IMGNUM];  }
    }

    // ------
    // apply mag-offset to each epoch-mag, unless mag is
    // set to flag value (UNDEFINED or ZEROFLUX)
    // Note that mag=99 --> zero flux is well defined.
    genmag8 = GENLC.genmag_obs[epoch];
    
    if ( genmag8 < 50.0 && genmag8 > -2.0 ) 
      {  genmag8 = GENLC.genmag_obs[epoch] + MAGOFF ; }

    if ( genmag8 > 600.0 ) 
      { genmag8 = MAG_UNDEFINED; } // avoid crazy values
   
    GENLC.genmag_obs[epoch] = genmag8 ;        // load global struct

    // -----------------------------
    // keep peak mags separately (before smearing)
    if ( GENLC.OBSFLAG_PEAK[epoch]  )  
      {  GENLC.peakmag_obs[ifilt_obs] = genmag8 ; }

    if ( GENLC.OBSFLAG_TEMPLATE[epoch] ) 
      { GENLC.genmag_obs_template[ifilt_obs] = genmag8; }

  } // end of epoch loop

  return ;

} // end genmag_offsets


// ********************************************
void prep_simpath(void) {

  // Created Jan 7, 2008 by R.Kessler
  //
  // * Create name of path for sim (PATH_SNDATA_SIM)
  // * Check that name of GENVERSION & path is not too long. 
  // * Use system call to create new subdir

  // May 2008: just create PATH_SNDATA_SIM; do NOT create subdir
  // Aug 14 2019: change suffix len from 20 to 7 to better trap
  //              too-long filename (see lensuffix)
  //

  char fnam[] = "prep_simpath" ;

  char tmp_path[2*MXPATHLEN] ;

  int lenpath, lenprefix, lensuffix, lenfile ;

    // ---------- BEGIN ----------

  sprintf(tmp_path, "%s/%s", 
	  PATH_SNDATA_SIM, INPUTS.GENVERSION );

  // check string lengths to avoid memory overwrites
  lenpath    = strlen(tmp_path);
  lenprefix  = strlen(INPUTS.GENPREFIX);
  lensuffix  = 7 ;      // e.g., '.README'
  lenfile    = lenpath + lenprefix + lensuffix ; 

  if ( lenprefix >= MXLEN_VERSION_PREFIX ) {
    sprintf(c1err,"GENPREFIX string len = %d exceeds array bound of %d",
	    lenprefix, MXLEN_VERSION_PREFIX);
    sprintf(c2err,"See input GENPREFIX: %s", INPUTS.GENPREFIX);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( lenfile >= MXPATHLEN  ) {
    sprintf(c1err, "Estimated filename length= %d is too long", lenfile);
    sprintf(c2err, "LEN(path,prefix,suffix) = %d, %d, %d : MXPATHLEN=%d",
	    lenpath, lenprefix, lensuffix, MXPATHLEN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  sprintf(PATH_SNDATA_SIM, "%s", tmp_path );

  return ;

}  // end of prep_simpath


// ************************************
void genperfect_override(void) {

  /****
   9/22/2007: 
   modify INPUTS variables to generate 'perfect'
   lightcurves with x1000 photostats, no zeropt-smearing,
   no galaxy noise and no host extinction.
  
   12/05/2007: turn off host-galaxy noise
   4/19/2009: update exposure time for each filter
   4/27/2009: include GENMAG_SMEAR -> 0

   4/29/2009: allow bitmask to select options.
            INPUT.GENPERFECT =  1 => set all bits 
            INPUT.GENPERFECT = -1 => list bit-options, then quit
 
   9/03/2009: turn off filter-dependent mag-smearing by 
              setting INPUTS.NFILT_SMEAR=0

  Nov 9, 2009: set GENSIGMA_RV=0 under mandatory category

  Jan 07, 2011: dump GENPERFECT value (bit #) instead of just bit #
                to aboid confusion among users

  Apr 06, 2012: set GENMAG_SMEAR_MODELNAME = "" for no mag-smear

  Oct 22, 2021: add GENPERFECT_HOSTLIB bit A Gagliano

  ****/

  int NVAR, ifilt, MASK, OVP, bit, MSKTMP ;
  double *dptr ;
  float  *fptr ;
  int    *iptr ;

#define BITPERFECT_ALL        0  // set all bits below
#define BITPERFECT_EXPTIME    1  // MASK=2: exposure time
#define BITPERFECT_MAGSMEAR   2  // MASK=4: intrinsic mag-smearing        
#define BITPERFECT_XTMW       3  // MASK=8: Milky Way extinction
#define BITPERFECT_AV         4  // MASK=16: host-galaxy extinction
#define BITPERFECT_HOSTLIB    5  // MASK=32: Hostlib A Gagliano 
#define MSKPERFECT_ALL       63  // mask with all bits

  char fnam[] = "genperfect_override" ;

  // ---------- BEGIN ---------------

  MASK = INPUTS.GENPERFECT ;
  if ( MASK == 0 )  return ;


  // check option to set all bits.  
  // Passing argument 1 sets all bits so that
  // it behaves as in previous versions.

  OVP = MASK & (1 <<  BITPERFECT_ALL ) ;
  if ( OVP > 0 && MASK > 0 ) MASK = MSKPERFECT_ALL ;

  // for negative mask, list bit-options and exit program.

  if ( MASK < 0 ) {

    bit = BITPERFECT_ALL ; MSKTMP = (1 <<  bit ) ;
    printf("  GENPERFECT %2d (bit %d) => set all bits below. \n", 
	   MSKTMP,bit);

    bit = BITPERFECT_EXPTIME ; MSKTMP = (1 <<  bit ) ;
    printf("  GENPERFECT %2d (bit %d) => x10^5 EXPOSURE TIME \n", 
	   MSKTMP, bit);

    bit = BITPERFECT_MAGSMEAR ; MSKTMP = (1 <<  bit ) ;
    printf("  GENPERFECT %2d (bit %d) => no intrinsic mag-smearing \n", 
	   MSKTMP, bit );

    bit = BITPERFECT_XTMW ; MSKTMP = (1 <<  bit ) ;
    printf("  GENPERFECT %2d (bit %d) => no Milky Way extinction \n", 
	   MSKTMP, bit );

    bit = BITPERFECT_AV ; MSKTMP = (1 <<  bit ) ;
    printf("  GENPERFECT %2d (bit %d) => no host-galaxy extinction (AV=0) \n", 
	   MSKTMP, bit );

    bit = BITPERFECT_HOSTLIB ; MSKTMP = (1 <<  bit ) ;
    printf("  GENPERFECT %2d (bit %d) => no HOSTLIB \n",
           MSKTMP, bit );


    happyend();

  }

  NVAR = 0;
  GENPERFECT.NVAR = NVAR ;


  // start with optional variables

  OVP = MASK & (1 <<  BITPERFECT_XTMW) ;
  if ( OVP > 0 ) {

    NVAR++ ;  iptr = &INPUTS.MWEBV_FLAG ;
    sprintf(GENPERFECT.parnam[NVAR], "MWEBV_FLAG" );
    GENPERFECT.parval[NVAR][0] = (float)*iptr ;
    *iptr = 0 ;
    GENPERFECT.parval[NVAR][1] = (float)*iptr ;
    GENPERFECT.partype[NVAR]   = 1 ;

    NVAR++ ;  iptr = &INPUTS.OPT_MWEBV ;
    sprintf(GENPERFECT.parnam[NVAR], "OPT_MWEBV" );
    GENPERFECT.parval[NVAR][0] = (float)*iptr ;
    *iptr = 0 ;
    GENPERFECT.parval[NVAR][1] = (float)*iptr ;
    GENPERFECT.partype[NVAR]   = 1 ;

    NVAR++ ;  iptr = &INPUTS.OPT_MWCOLORLAW ;
    sprintf(GENPERFECT.parnam[NVAR], "OPT_MWCOLORLAW" );
    GENPERFECT.parval[NVAR][0] = (float)*iptr ;
    *iptr = 0 ;
    GENPERFECT.parval[NVAR][1] = (float)*iptr ;
    GENPERFECT.partype[NVAR]   = 1 ;

  }


  OVP = MASK & (1 <<  BITPERFECT_EXPTIME ) ;
  if ( OVP > 0 ) {
    NVAR++ ;  iptr = &INPUTS.SMEARFLAG_ZEROPT ;
    sprintf(GENPERFECT.parnam[NVAR], "SMEARFLAG_ZEROPT" );
    GENPERFECT.parval[NVAR][0] = (float)*iptr ;
    *iptr = 0 ;
    GENPERFECT.parval[NVAR][1] = (float)*iptr ;
    GENPERFECT.partype[NVAR]   = 1 ;

    NVAR++ ;  fptr = &INPUTS.EXPOSURE_TIME ;
    sprintf(GENPERFECT.parnam[NVAR], "EXPOSURE_TIME" ) ;
    GENPERFECT.parval[NVAR][0] = *fptr ;
    *fptr = 100000. ;
    GENPERFECT.parval[NVAR][1] = *fptr ;
    GENPERFECT.partype[NVAR]   = 2 ;

    // update exposure time for each filter separately
    for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ )
      { INPUTS.EXPOSURE_TIME_FILTER[ifilt] *= *fptr ; }

  }


  OVP = MASK & (1 <<  BITPERFECT_AV ) ;
  if ( OVP > 0 ) {

    NVAR++ ;  dptr = &INPUTS.GENEXPTAU_AV ;
    sprintf(GENPERFECT.parnam[NVAR], "GENEXPTAU_AV" ) ;
    GENPERFECT.parval[NVAR][0] = *dptr ;
    *dptr = 0.0 ;
    GENPERFECT.parval[NVAR][1] = *dptr ;
    GENPERFECT.partype[NVAR]   = 2 ;

    NVAR++ ;  dptr = &INPUTS.GENGAUSIG_AV ;
    sprintf(GENPERFECT.parnam[NVAR], "GENGAUSIG_AV" ) ;
    GENPERFECT.parval[NVAR][0] = *dptr ;
    *dptr = 0.0 ;
    GENPERFECT.parval[NVAR][1] = *dptr ;
    GENPERFECT.partype[NVAR]   = 2 ;

    INPUTS.GENRANGE_AV[0] = 0.0 ;  // do this just to be safe
    INPUTS.GENRANGE_AV[1] = 0.0 ;

  }


  OVP = MASK & (1 <<  BITPERFECT_MAGSMEAR ) ;
  if ( OVP > 0 ) {
    NVAR++ ;  fptr = &INPUTS.GENMODEL_ERRSCALE ;
    sprintf(GENPERFECT.parnam[NVAR], "GENMODEL_ERRSCALE" ) ;
    GENPERFECT.parval[NVAR][0] = *fptr ;
    *fptr = 0.0 ;
    GENPERFECT.parval[NVAR][1] = *fptr ;
    GENPERFECT.partype[NVAR]   = 2 ;

    NVAR++ ;  fptr = INPUTS.GENMAG_SMEAR ;
    sprintf(GENPERFECT.parnam[NVAR], "GENMAG_SMEAR" ) ;
    GENPERFECT.parval[NVAR][0] = *fptr ;
    *fptr = 0.0 ;
    GENPERFECT.parval[NVAR][1] = *fptr ;
    GENPERFECT.partype[NVAR]   = 2 ;

    NVAR++ ;  iptr = &INPUTS.NFILT_SMEAR ;
    sprintf(GENPERFECT.parnam[NVAR], "NFILT_SMEAR" ) ;
    GENPERFECT.parval[NVAR][0] = (float)*iptr ;
    *iptr = 0 ;
    GENPERFECT.parval[NVAR][1] = (float)*iptr ;
    GENPERFECT.partype[NVAR]   = 1 ;

    sprintf(INPUTS.GENMAG_SMEAR_MODELNAME, "NONE");
  }

  OVP = MASK & (1 <<  BITPERFECT_HOSTLIB ) ;
  if ( OVP > 0 ){
   NVAR++ ;  iptr = &INPUTS.HOSTLIB_USE ;
   sprintf(GENPERFECT.parnam[NVAR], "HOSTLIB_USE" );
   GENPERFECT.parval[NVAR][0] = (float)*iptr ;
   *iptr = 0 ;
   GENPERFECT.parval[NVAR][1] = (float)*iptr ;
   GENPERFECT.partype[NVAR]   = 1 ;
 }

  // now do the mandatory ones.

  NVAR++ ;  dptr = &INPUTS.GENSIGMA_REDSHIFT ;
  sprintf(GENPERFECT.parnam[NVAR], "GENSIGMA_REDSHIFT");
  GENPERFECT.parval[NVAR][0] = *dptr ;
  *dptr = 1.0E-5 ;
  GENPERFECT.parval[NVAR][1] = *dptr ;
  GENPERFECT.partype[NVAR]   = 2 ;

  NVAR++ ;  fptr = &INPUTS.GENSIGMA_VPEC ;
  sprintf(GENPERFECT.parnam[NVAR], "GENSIGMA_VPEC");
  GENPERFECT.parval[NVAR][0] = *fptr ;
  *fptr = 0.0 ;
  GENPERFECT.parval[NVAR][1] = *fptr ;
  GENPERFECT.partype[NVAR]   = 2 ;


/* xxx mark delete
  NVAR++ ;  iptr = &INPUTS.HOSTLIB_USE ;
  sprintf(GENPERFECT.parnam[NVAR], "HOSTLIB_USE" );
  GENPERFECT.parval[NVAR][0] = (float)*iptr ;
  *iptr = 0 ;
  GENPERFECT.parval[NVAR][1] = (float)*iptr ;
  GENPERFECT.partype[NVAR]   = 1 ;
*/

  NVAR++ ;  fptr = &INPUTS.GENSIGMA_PEAKMJD ;
  sprintf(GENPERFECT.parnam[NVAR], "GENSIGMA_PEAKMJD" ) ;
  GENPERFECT.parval[NVAR][0] = *fptr ;
  *fptr = 0.1 ;
  INPUTS.OPT_SETPKMJD=0;
  GENPERFECT.parval[NVAR][1] = *fptr ;
  GENPERFECT.partype[NVAR]   = 2 ;


  NVAR++ ;  // fptr = &INPUTS.GENSIGMA_RV[0] ; float -> double 5/03/2012
  sprintf(GENPERFECT.parnam[NVAR], "GENSIGMA_RV" ) ;
  GENPERFECT.parval[NVAR][0] = 0.0 ; // *fptr ;
  *fptr = 0.0 ;
  GENPERFECT.parval[NVAR][1]  = 0.0 ;  // *fptr ;
  GENPERFECT.partype[NVAR]    = 2 ;
  INPUTS.GENGAUSS_RV.SIGMA[0] = 0.0 ;  
  INPUTS.GENGAUSS_RV.SIGMA[1] = 0.0 ;  // do this just to be safe


  GENPERFECT.NVAR = NVAR ;


} // end of genperfect_override




// *******************************************
void init_RateModel(void) {

  // Created Nov 26 2017
  int i;
  char fnam[] = "init_RateModel" ;

  // -------------- BEGIN ---------------

  if ( INDEX_GENMODEL == MODEL_SIMLIB ) { return; } // Oct 2020

  if ( INPUTS.GENRANGE_REDSHIFT[1] > 1.0E-8 )     
    { init_DNDZ_Rate(); } // Rate vs. redshift
  else
    { init_DNDB_Rate(); } // Rate vs. Gal coords l & b
  
  // --------------------------------------
  // now print info to stdout

  for ( i=0; i< NLINE_RATE_INFO; i++ ) 
    { printf("%s\n", LINE_RATE_INFO[i] ); }

  // July 30 2020
  // for INIT_ONLY mode, write explcit NGENTOT_CALC: <NGENTOT>
  // so that submit-batch script can parse a clean key and
  // avoid parsing awkwared human-readable text.
  if ( INPUTS.INIT_ONLY == 1 ) { 
    int NSN     = (int)(INPUTS.RATEPAR.SEASON_COUNT+0.5) ;
    int NPEC1A  = (int)(INPUTS.RATEPAR_PEC1A.SEASON_COUNT+0.5) ;
    int NGENTOT_CALC = NSN + NPEC1A ;
    printf("NGENTOT_RATECALC: %d\n", NGENTOT_CALC);
  }

  // debugexit(fnam); // xxx REMOVE

  return ;

} // end init_RateModel

// *******************************************
void init_DNDZ_Rate(void) {

  /****
   compute VOLUME and exposure time for GENRANGE.REDSHIFT
   and GENRANGE.PEAKMJD
   Just print result to screen to use for rate calc.

  Dec 2016: add optional CC rate from Strolger 2015/CANDELS

  Jan 19 2017: check TOTsum>0 to protect FRAC_PEC1A from NaN.
  Nov 27 2017: 
    + replace 'SN per season' with 'EVENTS per season'
    + rename function init_DNDZ_Rate (was VT_redshift)

  July 25 2019: abort if there is no rate model from DNDZ key.
  Mar  12 2021: abort of NGENTOT_LC > MXCID_SIM

  *****/

  double Z0, Z1, Z_AVG, ZVtmp[2], ZVint[2], ZVOL, VT;
  double dOmega, dOmega_user, dphi, dth, sin0, sin1;
  double delMJD, Tyear, Tcomoving, SNsum, PEC1Asum, TOTsum;
  double ztmp, rtmp, rtmp1, rtmp2, FRAC_PEC1A;

  int i, iz, OPT_DVDZ ;
  char cH0[8], cnum[20];
  char fnam[] = "init_DNDZ_Rate" ;

  // ----------- BEGIN ------------

  Z0 = INPUTS.GENRANGE_REDSHIFT[0] ;
  Z1 = INPUTS.GENRANGE_REDSHIFT[1] ;

  OPT_DVDZ = 0;
  ZVint[0] = dVdz_integral( OPT_DVDZ, Z0, &INPUTS.HzFUN_INFO );
  ZVint[1] = dVdz_integral( OPT_DVDZ, Z1, &INPUTS.HzFUN_INFO );

  // compute solid angle for stripe 82
  dphi   = INPUTS.GENRANGE_RA[1]   - INPUTS.GENRANGE_RA[0] ;
  dth    = INPUTS.GENRANGE_DEC[1]  - INPUTS.GENRANGE_DEC[0] ;

  sin1 = sin(INPUTS.GENRANGE_DEC[1]*RADIAN);
  sin0 = sin(INPUTS.GENRANGE_DEC[0]*RADIAN);

  dOmega_user = INPUTS.SOLID_ANGLE ;
  if ( dOmega_user == 0.0 ) { 
    dOmega  = dphi*RADIAN * ( sin1 - sin0 ) ;  // calculated dOmega
    INPUTS.SOLID_ANGLE = dOmega ;
  }
  dOmega  = INPUTS.SOLID_ANGLE ; 

  // get volume in (MPc/h)^3
  ZVOL = dOmega * ( ZVint[1] - ZVint[0] ) ;

  // compute average <z> wgted by volume
  OPT_DVDZ = 1;  // z-wgted option
  ZVtmp[0] = dVdz_integral( OPT_DVDZ, Z0, &INPUTS.HzFUN_INFO);  
  ZVtmp[1] = dVdz_integral( OPT_DVDZ, Z1, &INPUTS.HzFUN_INFO);

  /*
  printf(" xxx %s: ZVint = %f, %f \n", fnam, ZVint[0], ZVint[1] );
  printf(" xxx %s: ZVtmp = %f, %f \n", fnam, ZVtmp[0], ZVtmp[1] );
  */

  if ( Z0 == Z1 ) 
    { Z_AVG = Z1 ; }
  else
    { Z_AVG    = (ZVtmp[1] - ZVtmp[0]) / (ZVint[1] - ZVint[0]); }


  // get time in years
  delMJD    = INPUTS.GENRANGE_PEAKMJD[1] - INPUTS.GENRANGE_PEAKMJD[0] ;
  Tyear     = delMJD / 365.0 ;
  Tcomoving = Tyear / ( 1.0 + Z_AVG ) ;
  VT        = ZVOL * Tcomoving ;

  sprintf(cH0, "h%d" , (int)INPUTS.H0 );


  // ---------------------------------------------------
  // compute expected number of SN per season (SNsum):
  //
  //             /
  //            |            dV   Tyear
  //    SNsum = |  Rate(z) * -- * -----  dz
  //            |            dz    1+z
  //           /

  SNsum    = SNcount_model ( Z0, Z1, &INPUTS.RATEPAR ) ; // Ia or CC
  PEC1Asum = SNcount_model ( Z0, Z1, &INPUTS.RATEPAR_PEC1A ) ;

  if ( INDEX_GENMODEL == MODEL_NON1ASED ) {
    if ( SNsum > 0.0 && INPUTS.NON1ASED.NNON1A == 0 ) {
      sprintf(c1err,"DNDZ is specified, but found no NON1A SEDs.");
      sprintf(c2err,"Remove DNDZ key or add NON1A SEDs.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  if ( PEC1Asum > 0.0 && INPUTS.NON1ASED.NPEC1A == 0 ) {
    sprintf(c1err,"DNDZ_PEC1A is specified, but found no PEC1A SEDs.");
    sprintf(c2err,"Remove DNDZ_PEC1A key or add PEC1A SEDs.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // -----------------------------------------------------
  // store info in LINE_RATE_INFO array so that it can
  // be printed to screen and also to README file.

  i = -1 ;

  i++; sprintf(LINE_RATE_INFO[i], 
	       "\n *********************************** " );

  i++; sprintf(LINE_RATE_INFO[i],
	       "   SIMULATED VOLUME, TIME, RATE(%s)", 
	       INPUTS.RATEPAR.NAME );

  i++; sprintf(LINE_RATE_INFO[i], 
	 "\t Survey dOmega  = %8.4e steradians  (%8.5f PI) ", 
	 dOmega, 2.0*dOmega/TWOPI );

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t Redshift range = %7.4f - %7.4f", Z0, Z1 );

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t <redsfhit>     = %7.4f  (volume-weighted) ", Z_AVG ) ;

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t Survey Volume  = %8.4e  sr*(MPc/%s)^3 ", ZVOL, cH0 );

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t Survey Time    = %7.4f  years/season ", Tyear );

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t Co-moving Time = %7.4f  years/season  [ T/(1+<z>) ] ", 
	       Tcomoving );

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t Co-moving V*T  = %8.4e  sr*(MPc/%s)^3 * yr / season ", VT, cH0 );


  // June 4, 2012 - dNdz and rate info moved from README_doc()

  int  DNDZFLAG = 0;
  char ctmp_z[80], ctmp[100], ctmp_pec1a[80], *NAME;
  int  IMODEL_AB, IMODEL_PLAW, IMODEL_PLAW2, IMODEL_ZPOLY, IMODEL_FLAT ;
  int  IMODEL_CCS15, IMODEL_MD14, IMODEL_PISN, IMODEL_TDE, IMODEL_HUBBLE ;
  int  IFLAG_REWGT_ZEXP, IFLAG_REWGT_ZPOLY ;
  
  // model flags 
  NAME = INPUTS.RATEPAR.NAME ;
  IMODEL_AB      = ( strcmp(NAME,"ABMODEL"  )  == 0 ) ;
  IMODEL_PLAW    = ( strcmp(NAME,"POWERLAW" )  == 0 ) ;
  IMODEL_PLAW2   = ( strcmp(NAME,"POWERLAW2")  == 0 ) ;
  IMODEL_ZPOLY   = ( strcmp(NAME,"ZPOLY"    )  == 0 ) ;
  IMODEL_FLAT    = ( strcmp(NAME,"FLAT"     )  == 0 ) ;
  IMODEL_HUBBLE  = ( strcmp(NAME,"HUBBLE"   )  == 0 ) ;
  IMODEL_CCS15   = ( strcmp(NAME,RATEMODELNAME_CCS15)  == 0 ) ;
  IMODEL_MD14    = ( strcmp(NAME,RATEMODELNAME_MD14)   ==0 );

  IMODEL_PISN  = ( (strcmp(NAME,"IPP") ==0)  ||
		   (strcmp(NAME,RATEMODELNAME_PISN)==0) );

  IMODEL_TDE  = ( (strcmp(NAME,"EPM")==0) ||
		  (strcmp(NAME,RATEMODELNAME_TDE) ==0) ) ;

  // re-wgt flags
  IFLAG_REWGT_ZEXP  = ( INPUTS.RATEPAR.DNDZ_ZEXP_REWGT != 0.0    ) ;
  IFLAG_REWGT_ZPOLY = (INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT.ORDER > 0 ) ;

  // ----------------------------

  if ( IMODEL_AB ) {
    DNDZFLAG = 1;
    i++; sprintf(LINE_RATE_INFO[i],
		 "\t A+B MODEL: A(delayed) = %7.2e    B(prompt) = %7.2e", 
		 INPUTS.RATEPAR.MODEL_PARLIST[1][0], 
		 INPUTS.RATEPAR.MODEL_PARLIST[1][1] );
  }
  else if ( IMODEL_PLAW || IMODEL_PLAW2 ) {
    DNDZFLAG = 1;
    for ( iz=1; iz <= INPUTS.RATEPAR.NMODEL_ZRANGE; iz++ ) {
      sprintf(ctmp,"%7.2e*(1+z)^%4.2f",
	      INPUTS.RATEPAR.MODEL_PARLIST[iz][0], 
	      INPUTS.RATEPAR.MODEL_PARLIST[iz][1] );
      sprintf(ctmp_z,"%3.1f < z < %3.1f", 
	      INPUTS.RATEPAR.MODEL_ZRANGE[iz][0], 
	      INPUTS.RATEPAR.MODEL_ZRANGE[iz][1] );

      i++; sprintf(LINE_RATE_INFO[i],
		   "\t POWERLAW MODEL:  %s  (%s) ", ctmp, ctmp_z );
    }
  }
  else if ( IMODEL_ZPOLY ) {
    i++; 
    sprintf(LINE_RATE_INFO[i],"\t dN/dz = ZPOLY(%s)", 
	    INPUTS.RATEPAR.MODEL_ZPOLY.STRING) ;    
  }
  else if ( IMODEL_FLAT ) {
    i++; 
    sprintf(LINE_RATE_INFO[i], "\t dN/dz = FLAT  ~  (MJDrage * dOmega) ");  

    // Oct 22 2015
    // there is no rate for FLAT option, so make up SNsum so that
    // sim_SNmix will work with FLAT option.
    SNsum = delMJD * (Z1-Z0) * (dOmega * (41000./(2.*TWOPI)) ) ;
  }
  else if ( IMODEL_HUBBLE ) {
    i++; 
    sprintf(LINE_RATE_INFO[i], "\t dN/dz = constant/volume (HUBBLE)  ");  
  }
  else if ( IMODEL_CCS15 ) {
    DNDZFLAG = 1;
    i++; sprintf(LINE_RATE_INFO[i],
		 "\t dN/dz from %4.2f x Strolger15(CANDELS): ", 
		 INPUTS.RATEPAR.MODEL_PARLIST[1][0] );
  }
  else if ( IMODEL_PISN ) {
    DNDZFLAG = 1;
    i++; sprintf(LINE_RATE_INFO[i],
		 "\t dN/dz from PISN: " );
  }
  else if ( IMODEL_TDE ) {
    DNDZFLAG = 1;
    i++; sprintf(LINE_RATE_INFO[i],
		 "\t dN/dz from TDE: " );
  }
  else if ( IMODEL_MD14 ) {
    DNDZFLAG = 1;
    i++; sprintf(LINE_RATE_INFO[i],
		 "\t dN/dz = SFR(MD14,rV=%9.2le):  ", 
		 INPUTS.RATEPAR.MODEL_PARLIST[1][0] ) ;
    
  }
  else {
    sprintf(c1err,"Unknown rate model.");
    sprintf(c2err,"Check DNDZ key in sim-input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  // -------------------------
  // check re-wgt options.

  if ( IFLAG_REWGT_ZEXP ) {
    i++; 
    sprintf(LINE_RATE_INFO[i],
	    "\t Reweight dN/dz by z^(%4.2f) ", 
	    INPUTS.RATEPAR.DNDZ_ZEXP_REWGT);
  }
  else if ( IFLAG_REWGT_ZPOLY ) {
    i++; 
    sprintf(LINE_RATE_INFO[i],"\t Reweight dN/dz by ZPOLY(%s)",
	    INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT.STRING);
  }


  // if DNDZFLAG is set then print MODEL-RATE in 0.2 z-bins
  // from 0 to Z1

  if ( DNDZFLAG ) {
    ztmp = Z0 ;   ctmp_pec1a[0] = 0 ;
    while ( ztmp <= Z1 ) {
      rtmp1 = genz_wgt(ztmp,&INPUTS.RATEPAR) ;   
      rtmp2 = genz_wgt(ztmp,&INPUTS.RATEPAR_PEC1A) ; 
      rtmp  = rtmp1 + rtmp2 ;
      FRAC_PEC1A = rtmp2/rtmp ;
      if ( rtmp2 > 0.0 ) { sprintf(ctmp_pec1a,"(%.3f PEC1a)", FRAC_PEC1A); }
      sprintf(ctmp_z,"\t    MODEL-RATE(z=%4.2f) = %8.3e/Mpc^3/yr   %s ", 
	      ztmp, rtmp, ctmp_pec1a );
      i++; sprintf(LINE_RATE_INFO[i],"%s", ctmp_z);
      ztmp += 0.2 ;
    }
  }

  TOTsum  = SNsum + PEC1Asum ;
  if ( TOTsum > 0.0 ) 
    {  FRAC_PEC1A = PEC1Asum / TOTsum ; }
  else
    {  FRAC_PEC1A = 0.0 ; }

  INPUTS.RATEPAR.SEASON_COUNT          = SNsum ;
  INPUTS.RATEPAR_PEC1A.SEASON_COUNT    = PEC1Asum ;
  INPUTS.RATEPAR.SEASON_FRAC           = 1.0 - FRAC_PEC1A ;
  INPUTS.RATEPAR_PEC1A.SEASON_FRAC     = FRAC_PEC1A ;


  if ( TOTsum > 1.0E-12 ) {
    if ( TOTsum < 1.0 ) 
      { sprintf(cnum,"%le",   TOTsum) ; }
    else if ( TOTsum < 1.0E4 ) 
      { sprintf(cnum,"%5.0f",   TOTsum) ; }
    else                 
      { sprintf(cnum,"%10.3le", TOTsum) ; }

    i++; sprintf(LINE_RATE_INFO[i],
		 "\t Number of EVENTS per season = %s ", cnum );
    if ( FRAC_PEC1A > 0.0 ) {
      sprintf(LINE_RATE_INFO[i],"%s  (%.3f PEC1A)",
	      LINE_RATE_INFO[i], FRAC_PEC1A ) ;
    }
  }

  // Oct 26 2015: check season option 
  if ( INPUTS.NGEN_SEASON > 0.0 ) {
    double SCALE = INPUTS.NGEN_SCALE ;
    double XTOT  = (double)INPUTS.NGEN_SEASON * TOTsum * SCALE;
    if ( XTOT > (double)MXCID_SIM ) {
      sprintf(c1err,"NGENTOT_LC = %lld is too large", (long long)XTOT);
      sprintf(c2err,"MXCID_SIM  = %lld", MXCID_SIM);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( INPUTS.NGENTOT_LC > 0 ) {
      sprintf(c1err,"Cannot have both NGENTOT_LC>0 and NGEN_SEASON>0");
      sprintf(c2err,"At least one must be set to zero.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    INPUTS.NGENTOT_LC = (int)(XTOT);
    INPUTS.NGEN_LC    = 0 ;
    INPUTS.NGEN       = INPUTS.NGENTOT_LC ;
    i++; sprintf(LINE_RATE_INFO[i],
		 "\t NGEN_SEASON=%.2f --> NGENTOT_LC=%d ", 
		 INPUTS.NGEN_SEASON, INPUTS.NGENTOT_LC );
    set_screen_update(INPUTS.NGENTOT_LC) ;
  }

  i++; LINE_RATE_INFO[i][0] = 0 ;

  NLINE_RATE_INFO = i+1 ;


  return ;

}  // end of init_DNDZ_Rate

// ******************************
void init_DNDB_Rate(void) {

  // Created Nov 27 2017
  // Analog to init_DNDZ_Rate, summarize DNDB rate, 
  // but nothing is computed here.

  int i, j ;
  double *PARLIST     = INPUTS.RATEPAR.MODEL_PARLIST[1] ;
  int INDEX_RATEMODEL = INPUTS.RATEPAR.INDEX_MODEL ; 
  char varName[8] = "";
  char fnam[] = "init_DNDB_Rate" ;

  // --------------- BEGIN ---------------

  i = -1 ;

  i++; sprintf(LINE_RATE_INFO[i], 
	       "\n *********************************** " );

  i++; sprintf(LINE_RATE_INFO[i],
	       "   SIMULATED RATE(%s)",  INPUTS.RATEPAR.NAME );

  if ( INDEX_GENMODEL == MODEL_LCLIB ) {
    
    if ( INDEX_RATEMODEL == INDEX_RATEMODEL_COSBPOLY )  
      { sprintf(varName,"cosb"); }
    else if ( INDEX_RATEMODEL == INDEX_RATEMODEL_BPOLY ) 
      { sprintf(varName,"b"); }
    else {
      sprintf(c1err,"Unknown INDEX_RATEMODEL=%d for '%s'",
	      INDEX_RATEMODEL, INPUTS.RATEPAR.NAME );
      sprintf(c2err,"Check Galactic rate model.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    

    i++ ;
    sprintf(LINE_RATE_INFO[i],"\t dN/d%s = POLY(%s)",
	    varName, INPUTS.RATEPAR.MODEL_BPOLY.STRING);
    
  } // end MODEL_LCLIB

  // - - - - - - - 

  i++; sprintf(LINE_RATE_INFO[i],
	       "\t Number of EVENTS per season = %d  (NGENTOT_LC)", 
	       INPUTS.NGENTOT_LC );

  i++; LINE_RATE_INFO[i][0] = 0 ;

  NLINE_RATE_INFO = i+1 ;

  // load season count for INIT_ONLY=1 (Feb 2021)
  INPUTS.RATEPAR.SEASON_COUNT = (double)INPUTS.NGENTOT_LC ;

  return;


} // end init_DNDB_Rate

// ***************************
void init_simvar(void) {

  // One-time init of sim-variables; mostly counters.
  // Nov 24, 2017: init GENLC.MWEBV[_ERR] here instead of in init_GENLC().

  int ifilt, type, N ;
  float xmb ;
  char fnam[] = "init_simvar";

  // ----------- BEGIN -----------
 
  N = store_PARSE_WORDS(-1,""); // May 26 2021
  set_GENMODEL_NAME();

  init_GaussIntegral();
  ENVreplace("init", fnam, 1);

  GENLC.STOPGEN_FLAG = 0 ;
  GENLC.ACCEPTFLAG   = GENLC.ACCEPTFLAG_LAST = 0 ;
  
  set_FILTERSTRING(FILTERSTRING);
  set_EXIT_ERRCODE(EXIT_ERRCODE_SIM);

  //  sprintf(FILTERSTRING,"%s", FILTERSTRING_DEFAULT );

  NGENLC_TOT = NGENLC_WRITE = NGENSPEC_TOT = NGENSPEC_WRITE = 0;
  NGENFLUX_DRIVER = 0 ;

  NGEN_REJECT.GENRANGE  = 0;
  NGEN_REJECT.GENMAG    = 0;
  NGEN_REJECT.HOSTLIB   = 0;
  NGEN_REJECT.SEARCHEFF = 0;
  NGEN_REJECT.CUTWIN    = 0;
  NGEN_REJECT.NEPOCH    = 0;

  GENLC.MWEBV           = 0.0 ;
  GENLC.MWEBV_ERR       = 0.0 ;

  GENLC.NTYPE_SPEC_CUTS = 0;
  GENLC.NTYPE_SPEC      = 0 ;
  GENLC.NTYPE_PHOT_CUTS = 0;
  GENLC.NTYPE_PHOT      = 0 ;

  GENLC.NON1ASED.ISPARSE    = -9 ;  //  (non1a sparse index: 1-NINDEX)

  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
    GENLC.SNRMAX_FILT[ifilt]  = -9.0 ;   
    for(type=0; type < NTYPE_FLUXNOISE; type++ ) { 
      GENLC.MONITOR_REDCOV_FLUXNOISE[ifilt][type].RHO_SUM = 0.0 ; 
      GENLC.MONITOR_REDCOV_FLUXNOISE[ifilt][type].NSUM    = 0 ; 
    }
  }

  GENLC.FUDGE_SNRMAX_FLAG = 0 ;

  xmb = 1.0E-6 * (float)sizeof(INPUTS);
  printf("   sizeof(INPUTS) = %7.3f MB \n", xmb);

  xmb = 1.0E-6 * (float)sizeof(GENLC);
  printf("   sizeof(GENLC)  = %7.3f MB \n", xmb );


  // set all intrinsic scatter elements to zero
  ZERO_COVMAT_SCATTER();
  GENLC.TELESCOPE[0][0] = 0 ;

  NGEN_ALLSKIP = 0 ;

  NFUN_GENGAUSS_ASYM = 0 ;

  // init SIMLIB_HEADER counters (to help flag ABORT)
  SIMLIB_HEADER.NREPEAT      = 0 ;

  SIMLIB_HEADER.NFOUND_TOT   = 0 ;
  SIMLIB_HEADER.NFOUND_RA    = 0 ;
  SIMLIB_HEADER.NFOUND_DEC   = 0 ;
  SIMLIB_HEADER.NFOUND_MJD   = 0 ;
  SIMLIB_HEADER.NFOUND_FIELD = 0 ;
  SIMLIB_HEADER.NFOUND_GENCUTS = 0 ;

  SIMLIB_OBS_GEN.PIXSIZE[0] = -9.0 ; // ?? why is this here

  // init strong lens struct.
  GENSL.INIT_FLAG = GENSL.NIMG = GENSL.IDLENS = 0;
  GENSL.IMGNUM = -1;
  GENSL.PEAKMJD_noSL = -9.0 ;

  SPECTROGRAPH_USEFLAG = 0; // Jan 2021

  init_SNDATA_GLOBAL() ;  // Feb 2021

  // init dictionary of pre-scales vs. field for SIMLIB and TAKE_SPECTRUM
  init_string_dict(&INPUTS.DICT_SIMLIB_FIELDLIST_PRESCALE,
		   "SIMLIB_FIELDLIST_PRESCALES", 2*MXFIELD_OVP);
  init_string_dict(&INPUTS.DICT_SPECTRUM_FIELDLIST_PRESCALE, 
		   "SPECTRUM_FIELDLIST_PRESCALES", 2*MXFIELD_OVP);

  return ;

} // end of init_simvar


// ************************************
void  set_GENMODEL_NAME(void) {

  // April 18 2019
  // Called from init_simvar so that MODEL_NAMEs can be used
  // parsing GENMODEL without hard-wiring names.
  //
  // May 31 2019: allow SALT3 or SALT2

  int indx, j;

  // ----------- BEGIN ---------
  for ( indx=0; indx < MXMODEL_INDEX; indx++ ) {
    for (j=0; j < MXNAME_PER_MODEL; j++ ) {
      sprintf(GENMODEL_NAME[indx][j],"%s", "NULL");
    }
  }

  IS_PySEDMODEL = false ;

  // hard-wire list of valid GENMODELs
  // Note that some models allow for two different user-strings

  sprintf(GENMODEL_NAME[MODEL_STRETCH][0], "%s", "stretch" );
  sprintf(GENMODEL_NAME[MODEL_STRETCH][1], "%s", "stretch2" );

  sprintf(GENMODEL_NAME[MODEL_SALT2][0],   "%s", "SALT2"   );
  sprintf(GENMODEL_NAME[MODEL_SALT2][1],   "%s", "SALT3"   ); // May 30 2019

  sprintf(GENMODEL_NAME[MODEL_MLCS2k2][0], "%s", "mlcs2k2" );
  sprintf(GENMODEL_NAME[MODEL_MLCS2k2][1], "%s", "mlcs"    );

  sprintf(GENMODEL_NAME[MODEL_SNOOPY][0],  "%s", "snoopy"  );

  sprintf(GENMODEL_NAME[MODEL_S11DM15][0], "%s", "S11DM15" );

  sprintf(GENMODEL_NAME[MODEL_BYOSED][0],  "%s", "BYOSED"  ); // pyModel

  sprintf(GENMODEL_NAME[MODEL_SNEMO][0],   "%s", "SNEMO"   ); // pyModel

  sprintf(GENMODEL_NAME[MODEL_SIMSED][0],  "%s", "SIMSED"  );

  sprintf(GENMODEL_NAME[MODEL_FIXMAG][0],  "%s", "FIXMAG"  );
  sprintf(GENMODEL_NAME[MODEL_FIXMAG][1],  "%s", "RANMAG"  );
  sprintf(GENMODEL_NAME[MODEL_FIXMAG][2],  "%s", "fixmag"  );
  sprintf(GENMODEL_NAME[MODEL_FIXMAG][3],  "%s", "ranmag"  );

  sprintf(GENMODEL_NAME[MODEL_SIMLIB][0],  "%s", "SIMLIB"  );

  sprintf(GENMODEL_NAME[MODEL_NON1ASED][0], "%s", "NONIA"   );
  sprintf(GENMODEL_NAME[MODEL_NON1ASED][1], "%s", "NON1A"   );
  sprintf(GENMODEL_NAME[MODEL_NON1ASED][2], "%s", "NONIASED"   );
  sprintf(GENMODEL_NAME[MODEL_NON1ASED][3], "%s", "NON1ASED"   );

  sprintf(GENMODEL_NAME[MODEL_NON1AGRID][0],"%s", "NONIAGRID"  );
  sprintf(GENMODEL_NAME[MODEL_NON1AGRID][1],"%s", "NON1AGRID"  );

  sprintf(GENMODEL_NAME[MODEL_LCLIB][0], "%s", "LCLIB" );

  return ;

} // end set_GENMODEL_NAME




// ************************************
void  init_GENLC(void) {

  // called for each SN.
  // Oct    2011: return if SIMLIB_IDLOCK is set
  // Jan    2017: init GENLC.SNTYPE=0 instead of -9
  // Apr 6, 2017: move more stuff before USE_SAME_SIMLIB call.
  // Nov 06 2017: init SEARCHEFF_RANDOMS structure
  // Nov 22 2017: use NEP_RESET (instead of MXEPSIM) to limit
  //               wasted CPU on initializing (based on gprof)
  // Jul 20 2019: add skip for repeated strong lens images

  int epoch, ifilt, ifilt_obs, i, obs, imjd, NEP_RESET ;
  char fnam[] = "init_GENLC" ;

  // -------------- BEGIN ---------------

  GENLC.PEAKMAG_TRIGGER_FLAG = 0 ;
  GENLC.ACCEPTFLAG_LAST  = GENLC.ACCEPTFLAG ;
  GENLC.ACCEPTFLAG       = 0 ;
  GENLC.ACCEPTFLAG_FORCE = 0 ;

  GENLC.SEARCHEFF_MASK = 0 ;
  GENLC.SEARCHEFF_SPEC = 0.0 ;
  GENLC.SEARCHEFF_zHOST= 0.0 ;
  GENLC.CUTBIT_MASK    = 0 ;

  GENLC.NOBS           = 0 ;
  GENLC.NOBS_UNDEFINED = 0 ;
  GENLC.NOBS_SATURATE[0]  = 0 ; // NOBS NOT saturated
  GENLC.NOBS_SATURATE[1]  = 0 ; // NOBS saturated

  if ( NGENLC_TOT == 1 ) 
    { NEP_RESET = MXEPSIM ; }
  else
    { NEP_RESET = GENLC.NEPOCH + 1; }

  // reset data struct for trigger evaluation (Jan 2014)
  SEARCHEFF_DATA.NOBS        =  0 ;
  SEARCHEFF_DATA.REDSHIFT    = -9.0 ;
  SEARCHEFF_DATA.PEAKMJD     = -9.0 ;
  SEARCHEFF_DATA.DTPEAK_MIN  = -9.0 ;
  SEARCHEFF_DATA.SALT2mB     = -9.0 ;

  //  for(obs=0; obs < MXOBS_TRIGGER; obs++ ) {
  for(obs=0; obs < NEP_RESET; obs++ ) {
    SEARCHEFF_DATA.IFILTOBS[obs]    = -9   ;
    SEARCHEFF_DATA.MJD[obs]         = -9.0 ;
    SEARCHEFF_DATA.MAG[obs]         =  MAG_UNDEFINED ;
    SEARCHEFF_DATA.SNR[obs]         = -9.0 ;
    SEARCHEFF_DATA.detectFlag[obs]  =  0   ;
    SEARCHEFF_DATA.PHOTPROB[obs]    =  -9.0 ; 
    SEARCHEFF_RANDOMS.FLAT_PIPELINE[obs] = -9.0 ;
    SEARCHEFF_RANDOMS.FLAT_PHOTPROB[obs]  = -999.0 ;
    SEARCHEFF_RANDOMS.GAUSS_PHOTPROB[obs] = -999.0 ;
  }
  

  GENLC.MJD_TRIGGER      = -9.0 ;
  GENLC.MJD_DETECT_FIRST = -9.0 ;
  GENLC.MJD_DETECT_LAST  = -9.0 ;

  // Aug 9 2014: moved from gen_cutwin()
  GENLC.NOBS_SNR    = 0 ;
  GENLC.NOBS_MJDDIF = 0 ;
  GENLC.TRESTMIN   = +99999. ;
  GENLC.TRESTMAX   = -99999. ;
  GENLC.TGAPMAX    = -99999. ;  // max rest-frame gap among all filters
  GENLC.T0GAPMAX   = -99999. ;  // idem, near peak
 
  GENLC.SL_MAGSHIFT = 0.0 ;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // for repeated strong lens images, skip init (Jul 2019)
  if ( INPUTS_STRONGLENS.USE_FLAG  && GENSL.IMGNUM < GENSL.NIMG-1 )  
    { return; }

  /* 
      printf("\t xxx %s SKIP for SL IMGNUM=%d and NIMG=%d \n", 
	     fnam, GENSL.IMGNUM, GENSL.NIMG); fflush(stdout);
  */
  
 

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


  GENLC.METHOD_TYPE    = NULLINT ;
  GENLC.SNTYPE         = 0; 
  
  GENLC.CORRECT_HOSTMATCH = 1 ;  // default is correct match
  GENLC.SALT2gammaDM  = 0.0 ;

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    GENLC.DOFILT[ifilt_obs] = 1;  // do all user-filters by default.
  }

  // - - - - -
  // init shape-par pameters to NULLFLOAT;
  // only the selected model will get over-written

  GENLC.STRETCH    = NULLFLOAT ;
  GENLC.DELTA      = NULLFLOAT ;
  GENLC.DM15       = NULLFLOAT ;

  GENLC.SALT2x0    = NULLFLOAT ;
  GENLC.SALT2x1    = NULLFLOAT ;
  GENLC.SALT2c     = NULLFLOAT ;
  GENLC.SALT2beta  = NULLFLOAT ;
  GENLC.SALT2alpha = NULLFLOAT ;
  GENLC.SALT2mB    = NULLFLOAT ;

  GENLC.AV      = NULLFLOAT ;
  GENLC.RV      = NULLFLOAT ;

  GENLC.SNRMAX_GLOBAL           = -9.0 ;
  GENLC.IEPOCH_SNRMAX_GLOBAL    = -9 ;
  GENLC.IEPOCH_NEARPEAK         = -9 ;
  GENLC.DTPEAK_MIN              = 99999. ;

  GENLC.GENMAG_OFF_GLOBAL = 0.0 ;
  GENLC.LENSDMU = 0.0 ;
  GENLC.DLMU    = 0.0 ;

  // ------- Aug 2016 ----------
  // for NON1AGRID model, set random numbers here so that
  // gen_redshift knows to pick z from NON1A or PEC1A z distribution.
  // hard-wire 3 sigma clip for magSmear param.
  if ( INDEX_GENMODEL == MODEL_NON1AGRID  ) {
    double rmin=-3.0, rmax=3.0 ;
    int    ILIST_RAN = 2 ;
    GENLC.NON1AGRID_RANGauss = getRan_GaussClip(ILIST_RAN, rmin, rmax);
    GENLC.NON1AGRID_RANFlat  = getRan_Flat1(ILIST_RAN); 
  }


  // init filter-dependent stuff.
  for ( ifilt_obs=0; ifilt_obs < MXFILTINDX; ifilt_obs++ ) {
    GENLC.peakmag_obs[ifilt_obs]   = NULLFLOAT ;    

    GENLC.IEPOCH_PEAK[ifilt_obs]     = -9 ; 
    GENLC.IEPOCH_TEMPLATE[ifilt_obs] = -9 ; 
    GENLC.IEPOCH_SNRMAX[ifilt_obs]   = -9 ;
  
    GENLC.genmag_obs_template[ifilt_obs] = 99.0 ; // zero flux in template

    if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) 
      { GENLC.SIMLIB_USEFILT_ENTRY[ifilt_obs] = 1 ; }
    else
      { GENLC.SIMLIB_USEFILT_ENTRY[ifilt_obs] = 0 ; }

    
    GENLC.SNRMAX_FILT[ifilt_obs]    = -9.0 ;  
    GENLC.SNRMAX_SORTED[ifilt_obs]  = -9.0 ;  

    SEARCHEFF_RANDOMS.FLAT_SPEC[ifilt_obs] = -9.0 ;

    GENLC.NOBS_FILTER[ifilt_obs] = 0 ;
    GENLC.NOBS_SATURATE_FILTER[0][ifilt_obs] = 0 ;
    GENLC.NOBS_SATURATE_FILTER[1][ifilt_obs] = 0 ;

    GENLC.FLUXCOR_MWEBV_MAP[ifilt_obs]   = 1.0 ;
    GENLC.MAGCOR_MWEBV_MAP[ifilt_obs]    = 0.0 ;
    GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs]   = 0.0 ;
  }
  


  /*
  printf(" xxx NGENLC_TOT=%d  NEP_RESET=%d   (NEPOCH=%d)\n",
	 NGENLC_TOT, NEP_RESET, GENLC.NEPOCH );
  */


  GENLC.MAGSMEAR_COH[0] = 0.0 ;
  GENLC.MAGSMEAR_COH[1] = 0.0 ;
     
  for ( epoch = 0; epoch < NEP_RESET; epoch++ ) {
    GENLC.IFILT_OBS[epoch] = NULLINT ;

    GENLC.OBSFLAG_WRITE[epoch] = false ;

    sprintf(GENLC.FIELDNAME[epoch],"%s", FIELD_NONAME );

    GENLC.MJD[epoch]          = NULLFLOAT ;
    GENLC.epoch_obs[epoch]   = NULLFLOAT ;
    GENLC.epoch_rest[epoch]  = NULLFLOAT ;
    GENLC.epoch_obs_range[0] = +999999.0 ;
    GENLC.epoch_obs_range[1] = -999999.0 ;

    GENLC.flux[epoch]         = NULLFLOAT ;
    GENLC.fluxerr_data[epoch] = NULLFLOAT ;
    GENLC.mag[epoch]          = NULLFLOAT ;
    GENLC.mag_err[epoch]      = NULLFLOAT ;

    GENLC.genmag_obs[epoch]   = NULLFLOAT ;
    GENLC.generr_obs[epoch]   = NULLFLOAT ; // Apr 2013

    GENLC.genmag_rest[epoch]  = NULLFLOAT ;
    GENLC.generr_rest[epoch]  = 0.000 ;
    GENLC.genmag_rest2[epoch]  = NULLFLOAT ;
    GENLC.generr_rest2[epoch]  = NULLFLOAT ;

    GENLC.kcorval8[epoch]     = NULLFLOAT ;
    GENLC.warpcolval8[epoch]  = NULLFLOAT ;
    GENLC.AVwarp8[epoch]      = NULLFLOAT ;
    sprintf(GENLC.kcornam[epoch],    "NULL" );
    sprintf(GENLC.warpcolnam[epoch], "NULL" );

    GENLC.OBSFLAG_GEN[epoch] = true ; // default is to generate all obs
    GENLC.OBSFLAG_PEAK[epoch]      = false ;    
    GENLC.OBSFLAG_TEMPLATE[epoch]  = false ;    
    GENLC.SNR_CALC[epoch]    = 0.0 ;
    GENLC.SNR_MON[epoch]     = 0.0 ;

  } // end of epoch loop


  for ( i=0; i < MXZRAN; i++ )
    {  GENLC.REDSHIFT_RAN[i] = -9. ; }


  // init GENSPEC stuff
  if ( INPUTS.SPECTROGRAPH_OPTIONS.DOFLAG_SPEC ) {
    GENSPEC.NMJD_TOT = 0 ;
    GENSPEC.IMJD_HOST = -9 ;
    for(imjd=0; imjd < MXSPEC; imjd++ )  { GENSPEC_INIT(1,imjd); }
  }

  // keep track of last coord to skip parts of gen_MWEBV
  if ( NGENLC_TOT == 1 ) {
    // no previous RA,DEC on 1st event
    GENLC.RA_LAST         = -999.0 ;
    GENLC.DEC_LAST        = -999.0 ;
  }
  else {
    GENLC.RA_LAST   = GENLC.RA ;
    GENLC.DEC_LAST  = GENLC.DEC ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // skip rest of init if we are going to use the same LIBID
  if ( USE_SAME_SIMLIB_ID(1) != 0 ) { return ; } // Dec 2015
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  GENLC.NGEN_SIMLIB_ID = 0 ;
  GENLC.NEPOCH   = 0 ;
  GENLC.NEWMJD   = 0 ;

  GENLC.CID       = NULLINT ;
  GENLC.PEAKMJD   = NULLFLOAT ;
  GENLC.DTSEASON_PEAK   = 9999999.9 ;
  GENLC.MJD_EXPLODE = NULLFLOAT ;
  GENLC.ISOURCE_PEAKMJD = -9 ;
  GENLC.SDSS_SIM  = 0 ;
  GENLC.RA        = -999. ;
  GENLC.DEC       = -999. ;

  GENLC.random_shift_RA = GENLC.random_shift_DEC = 0.0 ;
  GENLC.random_shift_RADIUS = GENLC.random_shift_PHI = 0.0 ;

  GENLC.REDSHIFT_CMB          = NULLFLOAT ;
  GENLC.REDSHIFT_HELIO        = NULLFLOAT ;  
  GENLC.REDSHIFT_CMB_SMEAR    = NULLFLOAT ;
  GENLC.REDSHIFT_HELIO_SMEAR  = NULLFLOAT ;
  GENLC.REDSHIFT_SMEAR_ERR    = NULLFLOAT ;
  GENLC.REDSHIFT_FLAG         = REDSHIFT_FLAG_NONE ;
  GENLC.VPEC = GENLC.VPEC_SMEAR = 0.0 ;

  sprintf(GENLC.SNTYPE_NAME, "UNKNOWN" ); 

  for(i=0; i < MXOBS_SIMLIB; i++ ) {
    SIMLIB_OBS_GEN.IFILT_OBS[i]  = -9 ;
    SIMLIB_OBS_GEN.MJD[i]        = -99. ;
    SIMLIB_OBS_GEN.CCDGAIN[i]    = -99. ;
    SIMLIB_OBS_GEN.SKYSIG[i]     = -99. ;
    SIMLIB_OBS_GEN.READNOISE[i]  = -99. ;
    SIMLIB_OBS_GEN.PSFSIG1[i]    = -99. ;
    SIMLIB_OBS_GEN.PSFSIG2[i]    = -99. ;
    SIMLIB_OBS_GEN.PSFRATIO[i]   = -99. ;
    SIMLIB_OBS_GEN.NEA[i]        = -99. ; // Feb 2021
    SIMLIB_OBS_GEN.ZPTADU[i]     = -99. ;
    SIMLIB_OBS_GEN.ZPTERR[i]     = -99. ;
    SIMLIB_OBS_GEN.APPEND_PHOTFLAG[i] = 0 ;

    SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[i]     =  0.0 ;
    SIMLIB_OBS_GEN.TEMPLATE_READNOISE[i]  =  0.0 ;
    SIMLIB_OBS_GEN.TEMPLATE_ZPT[i]        = -99. ;

    // init a few SIMLIB_OBS_RAW quantities
    SIMLIB_OBS_RAW.APPEND_PHOTFLAG[i] = 0 ;
    SIMLIB_OBS_RAW.PTR_FIELDNAME[i] = SIMLIB_OBS_RAW.FIELDNAME[i] ;
    SIMLIB_OBS_RAW.PTR_BAND[i]      = SIMLIB_OBS_RAW.BAND[i] ;
    SIMLIB_OBS_RAW.FIELDNAME[i][0] = 0 ;
    SIMLIB_OBS_RAW.BAND[i][0] = 0 ;
  }

  return ;
 
}  // end of init_GENLC


// ***********************************
void init_modelSmear(void) {

  // Created Jul 1, 2010
  // one-time initialization for intrinsic SNIa model smearing
  // of the SN mags/fluxes. Note that the 'modelSmear' 
  // is defined independent of the SN model (MLCS,SALT2,DM15...).
  //
  // Jul 19, 2010: set DO_MODELSMEAR flag if any smear-option is set.
  //
  // Apr 06, 2012:
  //  Call new init_modelSmear_XXX functions, including
  //  new option GENMAG_SMEAR_MODELNAME
  //
  // Dec 2012 TODO: for NON1A, allow only global GENMAG_SMEAR.
  //   
  // Nov 14 2013: new model G10FUDGE to allow changing sigma_coh.
  //
  // May 01 2014: major re-org of SALT2-smearing model.
  //
  // Oct 09 2018: 
  //  + implement INPUTS.GENMAG_SMEAR_SCALE. See new SMEAR_SCALE 
  //    argument passed to   init_genSmear_FLAGS(SMEAR_SCALE);
  // 
  // Apr 11 2019:
  //  + adapt so that G10 model works for BYOSED
  //
  // May 27 2019: 
  //  + check for call to init_obs_atFLUXMAX() to prep for 
  //    PEAKMJD estimate
  //
  // Nov 30 2019: pass SMEAR_MODELNAME to init_genmag_COH
  //
  // Feb 11 2020: call init_genSmear_phaseCor(magSmear,expTau);
  //

  double GENMODEL_ERRSCALE   = (double)INPUTS.GENMODEL_ERRSCALE ;
  char  *SMEAR_SCALE_STRING  = INPUTS.GENMAG_SMEAR_SCALE;
  int    SMEAR_MSKOPT        = INPUTS.GENMAG_SMEAR_MSKOPT ;
  int    OPT, j, USE_SALT2smear ;
  double LAMRANGE[2], SIGCOH,  PARLIST_SETPKMJD[10];
  double magSmear, expTau;
  char *ptrName, key[40], NAM3[8]; 
  char MODELPATH_SALT2[MXPATHLEN];
  char fnam[] = "init_modelSmear"  ;

  // --------- BEGIN ----------

  sprintf(BANNER,"%s: init intrinsic SN smearing with model=%s",
	  fnam, INPUTS.GENMAG_SMEAR_MODELNAME);
  print_banner(BANNER);

  // print smearing mode
  if ( IFLAG_GENSMEAR == IFLAG_GENSMEAR_FILT ) {
    printf("\t Smear-mode: interpolate from central filter wavelengths.\n") ;
  } 
  else if ( IFLAG_GENSMEAR == IFLAG_GENSMEAR_LAM ) {
    printf("\t Smear-mode: complete wavelength dependence. \n");
  }
  else {
    sprintf(c1err,"Invalid FLAG_GENSMEAR = %d", IFLAG_GENSMEAR);
    sprintf(c2err,"Allowed values are %d or %d only",
	    IFLAG_GENSMEAR_FILT, IFLAG_GENSMEAR_LAM );
  }

  // ------------

  if (  INPUTS.GENMAG_SMEAR[0]  > 0. || GENMODEL_ERRSCALE  > 0.  ) 
    { INPUTS.DO_MODELSMEAR = 1 ; }

  // check passband magsmear 
  init_genSmear_filters();

  // internal inits for genSmear 
  init_genSmear_FLAGS(SMEAR_MSKOPT,SMEAR_SCALE_STRING);

  ptrName = INPUTS.GENMAG_SMEAR_MODELNAME ;

  if ( IGNOREFILE(ptrName) ) {  goto SKIP_GENSMEAR ;  }

  INPUTS.DO_MODELSMEAR  = 1 ;

    sprintf(key,"GENMAG_SMEAR_MODELNAME") ;
  if ( GENMODEL_ERRSCALE > 1.0E-9 ) {
    print_preAbort_banner(fnam);
    printf("  %s = %s \n" , key, ptrName);
    printf("  GENMODEL_ERRSCALE = %le \n", GENMODEL_ERRSCALE );
    sprintf(c1err, "Cannot set %s and GENMODEL_ERRSCALE.", key);
    sprintf(c2err, "Only one of the above is allowed.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // -------------------------------------
  // check for Guy10 wavelength-dependent color smear model.

  USE_SALT2smear = 0 ;
  SIGCOH = INPUTS.GENMAG_SMEAR_USRFUN[0]; // May 30 2018

  if ( INDEX_GENMODEL==MODEL_SALT2  || INDEX_GENMODEL==MODEL_BYOSED ) {
    OPT = 2;
    get_LAMRANGE_SEDMODEL( OPT, &LAMRANGE[0], &LAMRANGE[1] );

    if ( strstr(ptrName,"G10")      != NULL ) { USE_SALT2smear = 1; }
    if ( strstr(ptrName,"JLA")      != NULL ) { USE_SALT2smear = 1; }
    if ( strstr(ptrName,"SALT2")    != NULL ) { USE_SALT2smear = 1; }
    if ( strstr(ptrName,"G10Fig8")  != NULL ) { SIGCOH = 0.0 ; }

    if ( strstr(ptrName,"G10FUDGE") != NULL ) 
      {   SIGCOH  = INPUTS.GENMAG_SMEAR_USRFUN[0] ; }
    

    if ( INDEX_GENMODEL == MODEL_SALT2 ) { 
      sprintf(MODELPATH_SALT2,"%s", INPUTS.MODELPATH ); 
    }
    else if ( INDEX_GENMODEL == MODEL_BYOSED ) {
      sprintf(MODELPATH_SALT2,"%s", INPUTS.GENMAG_SMEAR_MODELARG); //BYOSED
      SIGCOH = 0.10; // hard-wired, Oct 31 2019 

      if ( USE_SALT2smear && strlen(MODELPATH_SALT2) == 0 ) { 
	sprintf(c1err,"For BYOSED model, missing SALT2PATH argument in ");
	sprintf(c2err,"GENMAG_SMEAR_MODELNAME: G10:[SALT2PATH]");  
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      } 
    }

    /*
    printf(" xxx MODELNAME='%s' USE=%d  SIGCOH=%f \n",
	   ptrName, USE_SALT2smear, SIGCOH);  debugexit(fnam);   */
  }
  else {
    // just guess with safety margin.
    LAMRANGE[0] = GENLC.RESTLAM_MODEL[0] - 2000.0 ;
    LAMRANGE[1] = GENLC.RESTLAM_MODEL[1] + 3000.0 ;
    if ( LAMRANGE[0] < 500.0 ) { LAMRANGE[0] = 500.0 ; }
  }


  // ---------------------------------------------
  // strncpy(NAM3,ptrName,3)  returns garbage, so use brute force
  for (j=0; j < 3; j++)  { sprintf( &NAM3[j], "%c", *(ptrName+j) ) ; }


  if ( strcmp(ptrName,"USRFUN") == 0 )  {  
      init_genSmear_USRFUN(INPUTS.NPAR_GENSMEAR_USRFUN, 
			   INPUTS.GENMAG_SMEAR_USRFUN, LAMRANGE ) ; 
  }
  // ----
  else if ( USE_SALT2smear ) {
    init_genSmear_SALT2(MODELPATH_SALT2, ptrName, SIGCOH, 
			INPUTS.GENRANGE_REDSHIFT);
  }

  // --------
  else if ( strcmp(ptrName,"Chotard11") == 0 ) 
    {  init_genSmear_Chotard11(0) ; }

  else if ( strcmp(ptrName,"C11") == 0 )     // same as Chotard11
    {  init_genSmear_Chotard11(0) ; }

  else if ( strcmp(ptrName,"C11_0") == 0 )     // same as Chotard11
    {  init_genSmear_Chotard11(0) ; }

  else if ( strcmp(ptrName,"C11_1") == 0 )     // 100% corr in far UV
    {  init_genSmear_Chotard11(1) ; }

  else if ( strcmp(ptrName,"C11_2") == 0 )     // 100% anti corr in far UV
    {  init_genSmear_Chotard11(2) ; }

  // -------

  else if ( strcmp(NAM3,"VCR") == 0 )     // Velocity-color relation
    {  init_genSmear_VCR(ptrName,INDEX_GENMODEL) ; }

  // --------

  else if ( strcmp(ptrName,"CCM89") == 0 )   // modify color law with CCM89
    {  init_genSmear_CCM89(LAMRANGE) ; }

  else if ( strstr(ptrName,"COH") != NULL ) 
    {  init_genSmear_COH(ptrName) ; }

  else if ( strcmp(ptrName,"BIMODAL_UV") == 0 ) 
    {  init_genSmear_biModalUV() ; }

  else if ( strstr(ptrName,"OIR.") != NULL ) 
    {  init_genSmear_OIR(ptrName) ; }

  else if ( strstr(ptrName,"COVSED.") != NULL ) 
    {  init_genSmear_COVSED(ptrName,0); }

  else if ( strcmp(ptrName,"PRIVATE") == 0 ) 
    {  init_genSmear_private() ; }

  else {
    sprintf(c1err,"Invalid smear model: '%s'", ptrName);
    sprintf(c2err,"Check GENMAG_SMEAR_MODELNAME key.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // -------------------------------
  if ( INPUTS.DO_MODELSMEAR  == 0 ) 
    { printf("\t ==> No model smearing options selected. \n" ) ; }
  else{ 
    // print magSmear sigma at various wavelengths
    //    dump_modelSmearSigma();
  }


 SKIP_GENSMEAR:

  // phase-dependent smearing is independent of GENMAG_SMEAR_MODELNAME
  magSmear = INPUTS.GENMAG_SMEAR_ADDPHASECOR[0];
  expTau   = max(INPUTS.GENMAG_SMEAR_ADDPHASECOR[1],1.0E-12);
  init_genSmear_phaseCor(magSmear,expTau) ;

  
  // May 2019: init method to estimate peakmjd for data files
  if ( INPUTS.GENSIGMA_PEAKMJD > 1.0E-9 ) // legacy Gauss smear
    { INPUTS.OPT_SETPKMJD = 0; }  
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID ) // do nothing for GRID
    { INPUTS.OPT_SETPKMJD=0;  INPUTS.GENSIGMA_PEAKMJD=0.0;  }

  OPT                 = INPUTS.OPT_SETPKMJD;
  PARLIST_SETPKMJD[0] = INPUTS.MJDWIN_SETPKMJD ;
  PARLIST_SETPKMJD[1] = INPUTS.SNRCUT_SETPKMJD ;
  PARLIST_SETPKMJD[2] = 3.0 ;          // hard-wired SNRCUT_BACKUP
  init_obs_atFLUXMAX(OPT,PARLIST_SETPKMJD,1);

  return ;

} // end of  init_modelSmear



// *****************************************
void dump_modelSmearSigma(void) {

  // Nov 23 2013
  // determine magSmear sigma in wavelength bins, and print to screen.

  int m = SNHOSTGAL.IMATCH_TRUE;
  int NLAM, NRANGEN, igen, ilam ;
  double LAMMIN, LAMMAX, LAMBIN, LAMARRAY[100], LAM, TREST, LOGMASS;
  double **MAGSMEAR, MAGARRAY[100], AVG, RMS, MEDIAN ;
  double parList[10];
  char fnam[] = "dump_modelSmearSigma" ;

  // --------------- BEGIN --------------

  sprintf(BANNER,"%s: Determine sigma(magSmear) vs. wavelength", fnam);
  print_banner(BANNER);

  LAMMIN = 2000.0 ;
  LAMMAX = GENLC.RESTLAM_MODEL[1] + 3000.0 ;
  LAMBIN = 500.0 ;

  TREST = 0.0 ;  LOGMASS = -9.0 ;
  GENLC.SALT2c = GENLC.AV = GENLC.SHAPEPAR = 0.0 ;

  NLAM = 0 ;
  for(LAM=LAMMIN; LAM < LAMMAX ; LAM += LAMBIN ) {
    LAMARRAY[NLAM] = LAM ;
    NLAM++ ;
  }

  NRANGEN = 5000; // number of random generations to compute sigma

  MAGSMEAR = (double**) malloc( NLAM * sizeof(double*) );
  for(ilam=0; ilam < NLAM; ilam++ ) 
    { MAGSMEAR[ilam] = (double*) malloc( NRANGEN * sizeof(double) ); }


  parList[0] = TREST ;
  parList[1] = GENLC.SHAPEPAR ;
  parList[2] = GENLC.SALT2c ;
  if ( m >= 0 ) { LOGMASS = SNHOSTGAL_DDLR_SORT[m].LOGMASS_TRUE ; }
  parList[3] = LOGMASS ;
    

  for(igen=0; igen < NRANGEN; igen++ ) {
    fill_RANLISTs();      // init list of random numbers 
    genran_modelSmear(); // load randoms for genSmear
    get_genSmear(parList, NLAM, LAMARRAY, MAGARRAY); // return MAGARRAY

    for(ilam=0; ilam < NLAM; ilam++ ) 
      { MAGSMEAR[ilam][igen] = MAGARRAY[ilam]  ; }

  } // end of igen loop
  

  for(ilam=0; ilam < NLAM; ilam++ ) {
    LAM = LAMARRAY[ilam] ;
    arrayStat(NRANGEN, MAGSMEAR[ilam], &AVG, &RMS, &MEDIAN);
    printf("\t LAM = %6.0f : magSmear(AVG,RMS) = %6.3f , %5.3f \n", 
	   LAM, AVG, RMS); fflush(stdout);
  }


  for(ilam=0; ilam < NLAM; ilam++ ) { free(MAGSMEAR[ilam]); }
  free(MAGSMEAR);

  debugexit("Graceful Exit") ;

} // end of dump_modelSmearSigma(void) {


// ***************************************
void init_genSmear_filters(void) {

  // Created April 2012 by R.Kessler
  // Move code form init_modelSmear().
  //
  // Init filter-dependent model smearing.
  // This makes sense for rest-frame models like mlcs and snoopy,
  // but for observer-models (SALT2) this model is redshift dependent.

  int   USESMEAR_LOCAL, j, ifilt ;
  float tmpSmear ;
  //  char  fnam[] = "init_genSmear_filters" ;

  // ----------- BEGIN ---------

  USESMEAR_LOCAL = 0 ;
  for ( j=1; j <= INPUTS.NFILT_SMEAR; j++ ) {
    ifilt    = INPUTS.IFILT_SMEAR[j]; 
    tmpSmear = INPUTS.GENMAG_SMEAR_FILTER[ifilt] ;
    if ( tmpSmear > 0.0 ) {
      USESMEAR_LOCAL = 1;
      INPUTS.DO_MODELSMEAR   = 1;
      printf("\t sigma_smear(%c) = %4.2f mag \n", 
	     FILTERSTRING[ifilt], INPUTS.GENMAG_SMEAR_FILTER[ifilt] ) ;
    }
  }

  if ( USESMEAR_LOCAL == 0 )  { INPUTS.NFILT_SMEAR = 0; }

  return ;

} // end of init_genSmear_filters


// ************************************
void  init_genSpec(void) {

  // Created Fall 2016.
  // one-time init for generating spectra.
  // See GENSPEC_DRIVER for execution.
  //
  // Jul 12 2019: allow BYOSED

  char *modelName = GENMODEL_NAME[INDEX_GENMODEL][0] ; // generic model name
  int OPTMASK     = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;
  char fnam[]     =  "init_genSpec" ;

  // -------------- BEGIN ----------------

  if ( SPECTROGRAPH_USEFLAG == 0 ) { return ; }
  INPUTS.SPECTROGRAPH_OPTIONS.DOFLAG_SPEC = 1 ;

  // check user option to turn off spectra
  if ( (OPTMASK & SPECTROGRAPH_OPTMASK_NOSPEC) > 0 ) 
    { INPUTS.SPECTROGRAPH_OPTIONS.DOFLAG_SPEC = 0; }

  if ( INDEX_GENMODEL == MODEL_SALT2     ||
       INDEX_GENMODEL == MODEL_NON1ASED  ||
       INDEX_GENMODEL == MODEL_SIMSED    || 
       IS_PySEDMODEL         ) {
    int     NB = INPUTS_SPECTRO.NBIN_LAM ;
    double *L0 = INPUTS_SPECTRO.LAMMIN_LIST ;
    double *L1 = INPUTS_SPECTRO.LAMMAX_LIST ;

    INIT_SPECTROGRAPH_SEDMODEL(modelName, NB, L0, L1); 
  }
  else if ( INDEX_GENMODEL == MODEL_FIXMAG ) { 
    // do nothing for spectra
  }
  else {
    sprintf(c1err,"Spectrograph option not available for");    
    sprintf(c2err,"INDEX_GENMODEL=%d (%s)", INDEX_GENMODEL, modelName);
    errmsg(SEV_WARN, 0, fnam, c1err, c2err ); 
  }

  // - - - - - - - - - - - - - - - - - - - -
  // check user options via sim-input key SPECTROGRAPH_OPTMASK: <mask>
  char tmpText[]=  "SPECTROGRAPH OPTION: " ;
  if ( (OPTMASK & SPECTROGRAPH_OPTMASK_LAMSIGMAx0)>0 )  { 
    INPUTS.SPECTROGRAPH_OPTIONS.SCALE_LAMSIGMA = 0.0; 
    printf("\t %s LAMSMEAR = 0 \n", tmpText );
  }
  if ( (OPTMASK & SPECTROGRAPH_OPTMASK_LAMSIGMAx2)>0 ) { 
    INPUTS.SPECTROGRAPH_OPTIONS.SCALE_LAMSIGMA = 2.0; 
    printf("\t %s LAMSMEAR *= 2 \n", tmpText );
  }
  if ( (OPTMASK & SPECTROGRAPH_OPTMASK_LAMSPIKE)>0 ) { 
    int ILAM_SPIKE = INPUTS_SPECTRO.NBIN_LAM/2 ;
    INPUTS.SPECTROGRAPH_OPTIONS.ILAM_SPIKE = ILAM_SPIKE ;
    printf("\t %s SED Flux only in ILAM_SPIKE=%d (%.1f A) \n", 
	   tmpText, ILAM_SPIKE, INPUTS_SPECTRO.LAMAVG_LIST[ILAM_SPIKE] );
  }
  if ( (OPTMASK & SPECTROGRAPH_OPTMASK_SNRx100)>0 ) { 
    INPUTS.SPECTROGRAPH_OPTIONS.SCALE_SNR = 100.0; 
    printf("\t %s SNR *= 100 \n", tmpText );
  }

  // this knob is tunable, not a fixed value from mask
  double *s = &INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE ;
  if ( *s != 1.0 ) {
    printf("\t %s Texpose *= %.3f \n", tmpText, *s );
  }

  if ( (OPTMASK & SPECTROGRAPH_OPTMASK_noTEMPLATE)>0 ) { 
    printf("\t %s no template noise \n", tmpText );
  }

  if ( (OPTMASK & SPECTROGRAPH_OPTMASK_NOSPEC)>0 ) { 
    printf("\t %s no spectra (but keep synthetic bands) \n", tmpText );
  }

  // - - - - - - - - - - - - - 
  // check option to scale lambda-resolution
  double SCALE_LAMSIGMA =  INPUTS.SPECTROGRAPH_OPTIONS.SCALE_LAMSIGMA ;
  int    ilam, NBLAM    =  INPUTS_SPECTRO.NBIN_LAM ;
  if ( SCALE_LAMSIGMA != 1.0 ) {
    for(ilam=0; ilam < NBLAM; ilam++ ) 
      { INPUTS_SPECTRO.LAMSIGMA_LIST[ilam] *= SCALE_LAMSIGMA ; }
  }

  GENSPEC.TEXPOSE_TEMPLATE = 0.0 ; // May 2021

  // - - - - - - - - - - - - - - - - - - - -
  // malloc arrays vs. lambda 
  int ispec, MXLAM = NBLAM ;
  int MEMD = MXLAM * sizeof(double);
  for(ispec=0; ispec < MXSPEC; ispec++ ) {
    GENSPEC.LAMMIN_LIST[ispec]           = (double*) malloc(MEMD) ;
    GENSPEC.LAMMAX_LIST[ispec]           = (double*) malloc(MEMD) ;
    GENSPEC.GENMAG_LIST[ispec]           = (double*) malloc(MEMD) ;
    GENSPEC.GENSNR_LIST[ispec]           = (double*) malloc(MEMD) ;
    GENSPEC.GENFLUX_LIST[ispec]          = (double*) malloc(MEMD) ;
    GENSPEC.GENFLUX_LAMSMEAR_LIST[ispec] = (double*) malloc(MEMD) ;
    GENSPEC.OBSFLUX_LIST[ispec]          = (double*) malloc(MEMD) ;
    GENSPEC.OBSFLUXERR_LIST[ispec]       = (double*) malloc(MEMD) ;
    GENSPEC.OBSFLUXERRSQ_LIST[ispec]     = (double*) malloc(MEMD) ;
    GENSPEC.GENFLAM_LIST[ispec]          = (double*) malloc(MEMD) ;
    GENSPEC.FLAM_LIST[ispec]             = (double*) malloc(MEMD) ;
    GENSPEC.FLAMERR_LIST[ispec]          = (double*) malloc(MEMD) ;
    GENSPEC.FLAMWARP_LIST[ispec]         = (double*) malloc(MEMD) ;
  } 

  if ( INPUTS.TAKE_SPECTRUM_HOSTSNFRAC > 0.000001 ) {   // Mar 2 2021
    GENSPEC.GENFLUX_PEAK          = (double*) malloc(MEMD) ;
    GENSPEC.GENMAG_PEAK           = (double*) malloc(MEMD) ;
  }

  // check if any spectrum has a warp
  if ( INPUTS.NWARP_TAKE_SPECTRUM > 0 ) 
    {  GENSPEC.USE_WARP = 1 ; }
  else
    {  GENSPEC.USE_WARP = 0 ; }


  GENSPEC.RANGauss_NOISE_TEMPLATE = (double*) malloc(MEMD) ;
  for(ilam=0; ilam < MXLAMSMEAR_SPECTROGRAPH; ilam++ ) 
    { GENSPEC.RANGauss_NOISE_TEMPLATE[ilam] = 0.0 ; }
 

  return ;

} // end init_genSpec


// *****************************************************
void rewrite_HOSTLIB_DRIVER(void) {

  char fnam[] = "rewrite_HOSTLIB_DRIVER" ;

  // ------------- BEGIN -------------

  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_PLUSMAGS)>0 ) 
    { rewrite_HOSTLIB_plusMags(); }

  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_PLUSNBR)>0 ) 
    { rewrite_HOSTLIB_plusNbr(); }

  if ( (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_APPEND )>0 ) 
    { rewrite_HOSTLIB_plusAppend(INPUTS.HOSTLIB_APPEND_FILE); }

} // end rewrite_HOSTLIB_DRIVER

// *****************************************************
void GENSPEC_DRIVER(void) {

  // Created July 2016
  // Driver for generating spectra (units: erg/s/cm^2/A)
  //
  // Sep  1 2016: add flat spectra for FIXRAN model.
  // Oct 16 2016: apply Gaussian LAMRES smearing
  // Jan 14 2021: abort if NMJD>0 but there is no SPECTROGRAPH instrument.
  // Feb 24 2021: increment NMJD_PROC only if NBLAM_VALID > 0
  // May 24 2021: check prescale for SN spectra

  int    NMJD = GENSPEC.NMJD_TOT  ;
  double MJD; 
  double SNR_LAMMIN, SNR_LAMMAX;
  
  int    i, imjd ;
  char fnam[] = "GENSPEC_DRIVER" ;

  // ------------ BEGIN ------------

  GENSPEC.NMJD_PROC = 0;

  // bail if no  spectra are requested
  if ( NMJD == 0 ) { return; }

  // if there is no SPECTROGRPAH instrument, abort
  if ( !SPECTROGRAPH_USEFLAG ) {
    sprintf(c1err,"Cannot generate %d spectra for CID=%d", 
	    NMJD, GENLC.CID );
    sprintf(c2err,"because SPECTROGRAPH is not defined in kcor file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // Jan 2018: allow some LIBIDs to NOT have a spectrograph key
  if ( NPEREVT_TAKE_SPECTRUM == 0 && SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH == 0 ) 
    { return ; }

  if ( NPEREVT_TAKE_SPECTRUM > 0 && NPEREVT_TAKE_SPECTRUM != NMJD ) {
    print_preAbort_banner(fnam);
    printf("  NPEREVT_TAKE_SPECTRUM = %d \n", NPEREVT_TAKE_SPECTRUM );
    printf("  NMJD_TOT = %d \n", NMJD);
    sprintf(c1err,"Cannot mix TAKE_SPECTRUM keys in sim-input file");
    sprintf(c2err,"with SPECTROGRAPH keys in SIMLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // - - - - - - - - - - - - - - - - 

  if ( INPUTS.TAKE_SPECTRUM_HOSTSNFRAC > 1.0E-8 ) 
    { GENSPEC_TRUE(ISPEC_PEAK); }

  int    NFIELD_OVP = SIMLIB_HEADER.NFIELD_OVP ;
  int    ifield, imjd_order[MXSPEC];
  double SNR, LAMMIN=100.0, LAMMAX=25000. ; 

  GENSPEC_MJD_ORDER(imjd_order); // check if nearPeak is done first

  for(ifield=0; ifield < NFIELD_OVP; ifield++ ) 
    { GENSPEC.FLATRAN_LIST[ifield]  = getRan_Flat1(1); } // for optional prescale

  // - - - - -
  for(i=0; i < NMJD; i++ ) {

    imjd = imjd_order[i];
    MJD  = GENSPEC.MJD_LIST[imjd] ;

    //    printf(" xxx %s: i=%d imjd=%d  MJD=%.2f \n",
    //	   fnam, i, imjd, GENSPEC.MJD_LIST[imjd] ); fflush(stdout);

    if ( !DO_GENSPEC(imjd) ) { continue; }

    SNR_LAMMIN = INPUTS.TAKE_SPECTRUM[imjd].SNR_LAMRANGE[0] ;
    SNR_LAMMAX = INPUTS.TAKE_SPECTRUM[imjd].SNR_LAMRANGE[1] ;
    if ( INPUTS.USE_SIMLIB_SPECTRA && SNR_LAMMIN > 0.1 ) {
      // clip spectra wavelength range to match that for SNR range.
      LAMMIN = SNR_LAMMIN ;
      LAMMAX = SNR_LAMMAX ;
    }
    
    GENSPEC_INIT(2,imjd);   // 2-> event-dependent init

    // compute true GENMAG and FLUXGEN in each lambda bin
    GENSPEC_TRUE(imjd); 

    // July 2019: option to add host contamination to SN spectrum
    GENSPEC_HOST_CONTAMINATION(imjd);

    // apply optional fudges for test or debug
    GENSPEC_FUDGES(imjd); 

    // if TAKE_SPECTRUM is defined by SNR, compute TEXPOSE
    GENSPEC_TEXPOSE_TAKE_SPECTRUM(imjd);

    // smear fluxes from Poisson noise and wavelength
    // Returned SNR is over entire wavelength range, but not used here.
    SNR = GENSPEC_SMEAR(imjd,LAMMIN,LAMMAX);

    // Feb 2 2017: convert flux to FLAM (dF/dlam)
    GENSPEC_FLAM(imjd) ;

    if ( GENSPEC.NBLAM_VALID[imjd] > 0 ) {
      GENSPEC.NMJD_PROC++ ; // total Nspec for this event.
      NGENSPEC_TOT++ ;      // total NSpec over all Light curve
    }
    else {
      GENSPEC.SKIP[imjd] = true; // 4.2021
    }

  } // end imjd
  

  return ;

} // end GENSPEC_DRIVER


// *********************************************
bool DO_GENSPEC(int imjd) {

  // Created May 27 2021
  // Return True to process spectrum for this imjd.

  STRING_DICT_DEF *DICT = &INPUTS.DICT_SPECTRUM_FIELDLIST_PRESCALE ;
  bool   DO             = false ;
  double MJD            = GENSPEC.MJD_LIST[imjd] ;

  char  *FIELD_REQUIRE  = INPUTS.TAKE_SPECTRUM[imjd].FIELD ;
  bool   CHECK_FIELD    = ( strlen(FIELD_REQUIRE) > 0 );
  bool   CHECK_PS       = ( DICT->N_ITEM > 0 ) ;
  int    NFIELD_OVP     = SIMLIB_HEADER.NFIELD_OVP ; 

  bool MATCH_FIELD, ACCEPT_FIELD, ACCEPT_PS ;
  int  OPT_DICT = 1 ;   // 0=exact match, 1=partial string match
  int  ifield ;
  double preScale, r1 ;
  char *field_tmp ; 
  char fnam[] = "DO_GENSPEC" ;

  // ------------ BEGIN ---------------

  // skip if outside Trest range
  if ( GENSPEC.SKIP[imjd] ) { return(false); } 

  // if there are no fields, skip FIELD-logic
  if ( NFIELD_OVP == 0 )   { return(true); }

  // - - - - - - - - - -
  // check each field for SELECT and PRESCALE(PS)
  // Beware the invalid FIELD in user input is not trapped
  // since we don't know a-priori which fields are in the SIMLIB.

  ACCEPT_PS = ACCEPT_FIELD = false;

  for(ifield=0; ifield < NFIELD_OVP; ifield++ ) {
    field_tmp    = SIMLIB_HEADER.FIELDLIST_OVP[ifield] ;

    // check field-dependent spectrum (host and SN)
    if ( CHECK_FIELD ) {
      MATCH_FIELD  = ( strcmp(field_tmp,FIELD_REQUIRE) == 0 );
      if ( MATCH_FIELD)  { ACCEPT_FIELD = true ; }   
    }
    else {
      ACCEPT_FIELD = true; 
    }

    // check field-dependent pre-scales for SN spectra (MJD>0);
    // host spectra always generated.
    if ( CHECK_PS && MJD > 0.0 ) {
      r1       = GENSPEC.FLATRAN_LIST[ifield] ;
      preScale = get_string_dict(OPT_DICT, field_tmp, DICT);
      if ( preScale < 0.0    ) { preScale  = 1.0  ; }
      if ( r1 < 1.0/preScale ) { ACCEPT_PS = true ; } 
    }
    else {
      ACCEPT_PS = true ;
    }

  } // end ifield

  // - - - - - - - - 

  DO = (ACCEPT_FIELD && ACCEPT_PS) ;

  // Jul 1 2021: set SKIP logical to account for pre-scale logic above.
  //     -> used later for writing data.
  if ( !DO ) { GENSPEC.SKIP[imjd] = true ; }

  return DO;

} // end DO_GENSPEC

// *************************************************
void  GENSPEC_LAMOBS_RANGE(int INDX, double *LAMOBS_RANGE) {

  // For input INDX (pointing to INPUTS.TAKE_SPECTRUM),
  // return obs-frame wavelength range in LAMOBS_RANGE.

  int    OPT_TEXPOSE      = INPUTS.TAKE_SPECTRUM[INDX].OPT_TEXPOSE ;
  int    OPT_FRAME_LAMBDA = INPUTS.TAKE_SPECTRUM[INDX].OPT_FRAME_LAMBDA ;
  double LAMMIN           = INPUTS.TAKE_SPECTRUM[INDX].SNR_LAMRANGE[0];
  double LAMMAX           = INPUTS.TAKE_SPECTRUM[INDX].SNR_LAMRANGE[1];
  double z    = GENLC.REDSHIFT_CMB ;
  double z1   = 1.0 + z ;

  //  char fnam[] = "GENSPEC_LAMOBS_RANGE" ;

  // ----------- BEGIN -----------

  LAMOBS_RANGE[0] = LAMOBS_RANGE[1] = -9.0 ;

  // make sure this TAKE_SPECTRUM uses the option to compute
  // Texpose based on specified SNR in lamrange.
  if ( OPT_TEXPOSE != 2 ) { return; }

  if ( OPT_FRAME_LAMBDA == GENFRAME_REST )  { 
    LAMOBS_RANGE[0] = LAMMIN * z1 ;  
    LAMOBS_RANGE[1] = LAMMAX * z1 ; 
  }
  else  { 
    LAMOBS_RANGE[0] = LAMMIN ;
    LAMOBS_RANGE[1] = LAMMAX ; 
  }

  return ;

} // end GENSPEC_LAMOBS_RANGE


// *************************************************
double  GENSPEC_PICKMJD(int OPT_MJD, int INDX, double z, 
			double *TOBS, double *TREST) {

  // For TAKE_SPECTRUM option , pick MJD.
  //
  // Inputs:
  //  OPT_MJD = 0 --> compute MJD at center of epoch range
  //          = 1 --> pick random MJD within epoch range
  //
  //  INDX = pointer to INPUTS.TAKE_SPECTRUM.
  //  z    = redshift
  //
  // Ouput:
  //   TOBS = MJD - GENLC.PEAKMJD
  //   TREST = TOBS/(1+z)
  //   GENSPEC_PICKMJD = returned MJD value
  //
  // Aug 23 2017: fix aweful bug setting EPOCH when OPT==0
  //              Somehow it was OK for LEGACY SIMLIB routines.
  //
  // Juj 28 2019: return 9999 for HOST spectrum

  int  OPT_FRAME = INPUTS.TAKE_SPECTRUM[INDX].OPT_FRAME_EPOCH ;
  int  ILIST_RAN = ILIST_RANDOM_SPECTROGRAPH ;
  double z1     = 1.0 + z ;
  double EPOCH_RANGE[2], EPOCH, MJD, Tobs, Trest ;
  char fnam[] = "GENSPEC_PICKMJD" ;

  // ------------ BEGIN ------------

  if ( OPT_FRAME == GENFRAME_HOST ) 
    {  *TOBS = *TREST = 9999.0;  MJD = -9.0 ; return(MJD); }
  
  EPOCH_RANGE[0]  = INPUTS.TAKE_SPECTRUM[INDX].EPOCH_RANGE[0] ;
  EPOCH_RANGE[1]  = INPUTS.TAKE_SPECTRUM[INDX].EPOCH_RANGE[1] ;


  if ( OPT_MJD == 0 ) 
    { EPOCH = 0.5 * (EPOCH_RANGE[1] + EPOCH_RANGE[0] ) ; }
  else
    { EPOCH = getRan_Flat(ILIST_RAN, EPOCH_RANGE); }  // pick random epoch

  z1 = 1.0 + z ; 
  if ( OPT_FRAME == GENFRAME_OBS ) 
    { Tobs = EPOCH ;    Trest = EPOCH/z1; }
  else if ( OPT_FRAME == GENFRAME_REST ) 
    { Tobs = EPOCH*z1 ; Trest = EPOCH ; }

  /*
  printf(" xxx  z1=%.3f  OPT=%d  INDX=%d  EPRANGE=(%5.1f,%5.1f) "
	 "EP=%.1f  Tobs=%.1f\n", 
	 z1, OPT_MJD, INDX,
	 EPOCH_RANGE[0], EPOCH_RANGE[1], EPOCH, Tobs ); fflush(stdout);
  */


  // load output arguments

  if ( OPT_FRAME == GENFRAME_MJD ) {
    // pick pre-defined MJD regardless of SN phase
    MJD    = EPOCH; 
    *TOBS  = MJD - GENLC.PEAKMJD ;
    *TREST = *TOBS / z1 ;

  }
  else {
    // pick MJD based on SN phase
    MJD    = GENLC.PEAKMJD + Tobs;
    *TOBS  = Tobs ;
    *TREST = Trest ;
  }

  /*
  printf(" xxx %s: CID=%d  INDX=%d  MJD=%.3f \n", 
	 fnam, GENLC.CID, INDX, MJD); 
  */

  return(MJD) ;

} // end GENSPEC_PICKMJD

// *************************************************
void  GENSPEC_MJD_ORDER(int *imjd_order) {

  // Created Mar 2017
  // determine order to generate spectra.
  // Default is by increasing MJD (as in the SIMLIB).
  // But for TAKE_SPECTRUM option with scaled template 
  // exposure time, the spectrum closest to max must
  // be processed first in order to determine 
  // TEXPOSE(TEMPLATE) for all other epochs.
  //
  // Always process HOST spectra first to allow for host 
  // contamination in the SN spectra.
  //
  //
 
  int  imjd;
  float SCALE = INPUTS.TAKE_SPECTRUM_TEMPLATE_TEXPOSE_SCALE ;
  char fnam[] = "GENSPEC_MJD_ORDER" ;

  // -------------- BEGIN -------------

  GENSPEC.IMJD_NEARPEAK = -9;

  if ( SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH > 0 ) {

    // SPECTROGRPAH keys in SIMLIB file
    for(imjd=0; imjd < GENSPEC.NMJD_TOT; imjd++ ) 
      { imjd_order[imjd] = imjd; }
  }
  else {
    // TAKE_SPECTRUM
    // epoch closet to peak is first in list;
    // after that it doesn't matter.
    if ( NPEREVT_TAKE_SPECTRUM == 0 ) {
      sprintf(c1err,"TEMPLATE_TEXPOSE_SCALE = %.2f", SCALE);
      sprintf(c2err,"but no TAKE_SPECTRUM keys found.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    // find epoch closest to peak
    double TREST, TRESTmin=99999.0 ;
    int ntmp, imjd_nearPeak=-9,  nhost=0, IS_HOST[MXSPEC] ;
    for(imjd=0; imjd < GENSPEC.NMJD_TOT; imjd++ ) {

      IS_HOST[imjd] = 0;

      // store hosts first (note there can be more than one host spec)
      if ( GENSPEC.MJD_LIST[imjd] < 0.0 )  
	{ imjd_order[nhost] = imjd;  nhost++ ; IS_HOST[imjd] = 1; }

      TREST = fabs(GENSPEC.TREST_LIST[imjd]) ;
      if ( TREST < TRESTmin ) { TRESTmin = TREST;  imjd_nearPeak=imjd; }
    }

    // set nearest-peak as first SN spectrum.
    imjd_order[nhost] = imjd_nearPeak ;  ntmp=nhost ;
    GENSPEC.IMJD_NEARPEAK = imjd_nearPeak ;

    // store remaining spectra
    for(imjd=0; imjd < GENSPEC.NMJD_TOT; imjd++ ) {
      if ( imjd == imjd_nearPeak ) { continue; }
      if ( IS_HOST[imjd]         ) { continue; }
      ntmp++ ; imjd_order[ntmp] = imjd;
    }

  } // end if block


  return ;

} // end GENSPEC_MJD_ORDER

// *************************************************
void GENSPEC_INIT(int OPT, int imjd) {

  // Init imjd & lambda-arrays.
  // OPT=1 --> init everything for one-time init
  // OPT=2 --> init some stuff each event.

  int  NBLAM = INPUTS_SPECTRO.NBIN_LAM ;
  int  ilam ;
  //  char fnam[] = "GENSPEC_INIT" ;

  // ------------ BEGIN ----------

  GENSPEC.NBLAM_VALID[imjd]   =  0 ;
  GENSPEC.NBLAM_TOT[imjd]     = INPUTS_SPECTRO.NBIN_LAM ;

  if ( OPT == 1 ) {
    GENSPEC.MJD_LIST[imjd]      = -9.0 ;
    GENSPEC.TOBS_LIST[imjd]     = -9.0 ;
    GENSPEC.TREST_LIST[imjd]    = -9.0 ;
    GENSPEC.TEXPOSE_LIST[imjd]        = -9.0 ;
    GENSPEC.OPT_TEXPOSE_LIST[imjd]    = -9 ;
    GENSPEC.INDEX_TAKE_SPECTRUM[imjd] = -9 ;

    GENSPEC.SNR_REQUEST_LIST[imjd] = -9.0 ;
    GENSPEC.SNR_COMPUTE_LIST[imjd] = -99.0 ;
    GENSPEC.IS_HOST[imjd]          = false;
  }


  // init arrays
  for(ilam=0; ilam < NBLAM; ilam++ ) {
    GENSPEC.GENMAG_LIST[imjd][ilam]           =  0.0 ;
    GENSPEC.GENSNR_LIST[imjd][ilam]           =  0.0 ;
    GENSPEC.GENFLUX_LIST[imjd][ilam]          =  0.0 ;
    GENSPEC.GENFLUX_LAMSMEAR_LIST[imjd][ilam] =  0.0 ;

    GENSPEC.OBSFLUX_LIST[imjd][ilam]          =  0.0 ;
    GENSPEC.OBSFLUXERR_LIST[imjd][ilam]       =  0.0 ;
    GENSPEC.OBSFLUXERRSQ_LIST[imjd][ilam]     =  0.0 ;

    GENSPEC.GENFLAM_LIST[imjd][ilam]          =  0.0 ; // true dF/dlam
    GENSPEC.FLAM_LIST[imjd][ilam]             =  0.0 ; // dF/dlam
    GENSPEC.FLAMERR_LIST[imjd][ilam]          =  0.0 ; // error on above
    GENSPEC.FLAMWARP_LIST[imjd][ilam]         =  1.0 ;
  }

  return ;

} // end GENSPEC_INIT


void GENSPEC_OBSFLUX_INIT(int imjd, int ILAM_MIN, int ILAM_MAX) {

  int ilam;
  for(ilam = ILAM_MIN; ilam <= ILAM_MAX; ilam++ ) {
    GENSPEC.OBSFLUX_LIST[imjd][ilam]          =  0.0 ;
    GENSPEC.OBSFLUXERR_LIST[imjd][ilam]       =  0.0 ;
    GENSPEC.OBSFLUXERRSQ_LIST[imjd][ilam]     =  0.0 ;
  }
  GENSPEC.NBLAM_VALID[imjd] = 0;

} // end  GENSPEC_OBSFLUX_INIT

// *************************************************
void GENSPEC_TRUE(int imjd) {

  // Load true GENMAG and FLUXGEN for each lambda bin at this
  // "imjd".
  // Fill arrays
  //    GENSPEC.GENMAG_LIST[imjd][ilam] 
  //    GENSPEC.GENFLUX_LIST[imjd][ilam]
  //
  // Dec 04 2018: call genSpec_BYOSED
  // Jun 28 2019: call genSpec_HOST if IS_HOST is true
  // Mar 02 2021: generate PEAK spectrum if imjd == ISPEC_PEAK

  int  NBLAM      = INPUTS_SPECTRO.NBIN_LAM ;
  double TOBS, GENMAG, ZP, ARG, FLUXGEN, MAGOFF ;
  double *ptrGENFLUX, *ptrGENMAG ;

  double x0 = GENLC.SALT2x0 ;
  double x1 = GENLC.SALT2x1 ;
  double c  = GENLC.SALT2c ;
  double RV = GENLC.RV;
  double AV = GENLC.AV ;
  double  logMass  = -9.0 ;

  double parList_SN[4]   = { x0, x1, c, x1 } ;
  double parList_HOST[3] = { RV, AV, logMass } ;
  int IS_HOST, ilam ;
  int DUMPFLAG=0;
  char fnam[] = "GENSPEC_TRUE" ;

  // --------------- BEGIN ----------------

  if ( imjd == ISPEC_PEAK ) {
    TOBS       = 0.0 ;  
    IS_HOST    = 0 ;
    ptrGENFLUX = GENSPEC.GENFLUX_PEAK ;
    ptrGENMAG  = GENSPEC.GENMAG_PEAK ;
  }
  else {
    TOBS       = GENSPEC.TOBS_LIST[imjd] ;
    IS_HOST    = GENSPEC.IS_HOST[imjd] ;
    ptrGENFLUX = GENSPEC.GENFLUX_LIST[imjd] ;
    ptrGENMAG  = GENSPEC.GENMAG_LIST[imjd] ;
  }

  // - - - - - - - - - - - 

  if ( IS_HOST ) {    
    genSpec_HOSTLIB(GENLC.REDSHIFT_HELIO,         // (I) helio redshift
		    GENLC.MWEBV,                  // (I) Galactic extinction
		    DUMPFLAG,                     // (I)
		    ptrGENFLUX,       // (O) true fluxGen per bin 
		    ptrGENMAG );      // (O) magGen per bin
    return;
  }

  // below is a SN spectrum, so check which model.

  if ( INDEX_GENMODEL == MODEL_SALT2 )  {
    genSpec_SALT2(parList_SN, parList_HOST,
		  GENLC.MWEBV,             // Galactic		 
		  GENLC.REDSHIFT_HELIO, TOBS,
		  ptrGENFLUX,        // (O) fluxGen per bin (erg/s/cm^2)
		  ptrGENMAG          // (O) magGen per bin
		  );
  }
  else if ( INDEX_GENMODEL == MODEL_FIXMAG )  {
    for(ilam=0; ilam < NBLAM; ilam++ ) {
      GENMAG  = GENLC.FIXMAG ;   
      ZP      = ZEROPOINT_FLUXCAL_DEFAULT ;
      ARG     = -0.4*(GENMAG-ZP);
      FLUXGEN = pow(TEN,ARG);
      
      ptrGENMAG[ilam]  = GENMAG ;
      ptrGENFLUX[ilam] = FLUXGEN ;  // FLUXCAL units
    }
  }
  else  if ( INDEX_GENMODEL  == MODEL_NON1ASED ) {    
    // works with NON1ASED, but NOT with NON1AGRID
    MAGOFF = 
      GENLC.GENMAG_OFF_GLOBAL + GENLC.MAGSMEAR_COH[0] + GENLC.MAGSMEAR_COH[1];

    getSpec_SEDMODEL(ISED_NON1A,
		     GENLC.MWEBV, GENLC.RV, GENLC.AV, // (I)   
		     GENLC.REDSHIFT_HELIO,            // (I) redshift
		     GENLC.DLMU,                      // (I) dist mod
		     TOBS, MAGOFF,                    // (I) Tobs, magoff 
		     ptrGENFLUX,           // (O) fluxGen per bin 
		     ptrGENMAG             // (O) magGen per bin
		     ); 

  }
  else  if ( INDEX_GENMODEL  == MODEL_NON1AGRID ) {    
    sprintf(c1err,"Spectrograph option not available for NON1AGRID;");    
    sprintf(c2err,"try NON1ASED instead.");
    errmsg(SEV_WARN, 0, fnam, c1err, c2err ); 
  }
  else if ( INDEX_GENMODEL == MODEL_SIMSED ) {

    if ( NPEREVT_TAKE_SPECTRUM > 0 ) {
      sprintf(c1err,
	      "\n\tSIMSED model is not yet compatiable with SPECTROGRAPH");
      sprintf(c2err,"Try NON1ASED or SALT2 model instead");
      errmsg(SEV_WARN, 0, fnam, c1err, c2err ); 
    }
    

    MAGOFF = GENLC.GENMAG_OFF_GLOBAL ;
    getSpec_SEDMODEL(ISED_SEDMODEL,
		     GENLC.MWEBV, GENLC.RV, GENLC.AV,  // (I)  
		     GENLC.REDSHIFT_HELIO,        // (I) redshift
		     GENLC.DLMU,
		     TOBS, MAGOFF,                // (I) Tobs, magoff 
		     ptrGENFLUX,         // (O) fluxGen per bin 
		     ptrGENMAG           // (O) magGen per bin
		     ); 
  }
  else if ( IS_PySEDMODEL ) {
    
    genSpec_PySEDMODEL(TOBS, 
		       GENLC.REDSHIFT_HELIO,            // (I) helio redshift
		       GENLC.DLMU,                      // (I) dist. mod.
		       GENLC.MWEBV, GENLC.RV, GENLC.AV, // (I)		     
		       ptrGENFLUX,      // (O) fluxGen per bin 
		       ptrGENMAG        // (O) magGen per bin
		       );		
   
  }
  else { 
    /*  don't abort since init_genSpec gives warning.  
	sprintf(c1err,"Cannot make spectrum for model=%s", modelName);
	sprintf(c2err,"Check manual for valid models");
	errmsg(SEV_WARN, 0, fnam, c1err, c2err ); 
    */
  }


  // Aug 2021: check option to integrate flux within each band
  //            and check peak mags.

  int VERIFY_PEAKMAG = 0 ;
  int ifilt, ifilt_obs;
  double PEAKMAG ;
  if ( TOBS == 0.0 && VERIFY_PEAKMAG ) {
    printf("\n %s: Verify PEAKMAGs for CID=%d  z=%.4f: \n", 
	   fnam, GENLC.CID, GENLC.REDSHIFT_CMB);
    for(ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      GENSPEC_VERIFY_PEAKMAG(ifilt_obs, ptrGENFLUX);
    }
  }

  return ;

} // end GENSPEC_TRUE

// ==================================
void GENSPEC_VERIFY_PEAKMAG(int ifilt_obs, double *GENFLUX_LIST) {

  // Created Aug 26 2021
  // Diagnostic:
  // Convert ptr_GENMAG in each wave bin to flux, sum flux
  // over filter transmission, then convert back to broadband mag.
  // Compare with already-computed PEAKMAG.

  int  NLAMSPEC       = SPECTROGRAPH_SEDMODEL.NBLAM_TOT ;
  double *LAMAVG_LIST = SPECTROGRAPH_SEDMODEL.LAMAVG_LIST ;
  double *LAMMAX_LIST = SPECTROGRAPH_SEDMODEL.LAMMAX_LIST ;
  double *LAMMIN_LIST = SPECTROGRAPH_SEDMODEL.LAMMIN_LIST ;
  double ZP_SNANA = ZEROPOINT_FLUXCAL_DEFAULT ;
  double hc8    = (double)hc ;
  double z1     = 1.0 + GENLC.REDSHIFT_CMB;
  double PEAKMAG, PEAKMAG_VERIFY, LAMOBS, TRANS, flam, flux, flux_sum=0.0 ;
  double *FLAM_LIST, lamstep, lambin, ZP ; 
  int  ifilt, NLAMFILT, ilam;
  char *cfilt;
  int  OPT_INTERP = 1; // linear
  char fnam[] = "GENSPEC_VERIFY_PEAKMAG";

  // --------- BEGIN ------------
  PEAKMAG = GENLC.peakmag_obs[ifilt_obs] ;
  if ( PEAKMAG < 1.0 || PEAKMAG > 30.0 )  { return ; }

  // compute true Flam in each spectrograph bin
  FLAM_LIST = (double*)malloc( NLAMSPEC * sizeof(double) );
  for(ilam=0; ilam < NLAMSPEC; ilam++ ) {
    flux   = GENFLUX_LIST[ilam] ;
    lambin = LAMMAX_LIST[ilam] - LAMMIN_LIST[ilam];
    flam   = flux / lambin ;
    FLAM_LIST[ilam] = flam ;
  }

  ifilt     = IFILTMAP_SEDMODEL[ifilt_obs] ;
  NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;
  cfilt     = FILTER_SEDMODEL[ifilt].name ;
  lamstep   = FILTER_SEDMODEL[ifilt].lamstep ;
  ZP        = FILTER_SEDMODEL[ifilt].ZP ;

  // loop over filter wave bins
  for ( ilam=0; ilam < NLAMFILT; ilam++ ) {
    get_LAMTRANS_SEDMODEL(ifilt, ilam, &LAMOBS, &TRANS );
    if ( TRANS < 1.0E-12 ) { continue; }

    flam = interp_1DFUN(OPT_INTERP, LAMOBS, NLAMSPEC,
			LAMAVG_LIST, FLAM_LIST, fnam );

    flux_sum += ( flam * LAMOBS * TRANS );
  } // end ilam

  flux_sum *= (lamstep/hc8) ;
  PEAKMAG_VERIFY = ZP - 2.5*log10(flux_sum);

  printf(" %s: Peakmag(%s) = %.3f / %.3f (orig/specVerify) \n",
	 fnam, cfilt, PEAKMAG, PEAKMAG_VERIFY );
  fflush(stdout);

  free(FLAM_LIST);

  return ;

} // end GENSPEC_VERIFY_PEAKMAG

// *****************************************
void GENSPEC_HOST_CONTAMINATION(int imjd) {

  // Created July 11 2019
  // Check option to add host contamination
  // Feb 26 2021: check for HOST/SN fraction: HOSTSNFRAC
  // Jun 07 2021: fix IMJD_HOST

  int    IS_HOST    = GENSPEC.IS_HOST[imjd];
  double HOSTFRAC   = (double)INPUTS.TAKE_SPECTRUM_HOSTFRAC;
  double HOSTSNFRAC = (double)INPUTS.TAKE_SPECTRUM_HOSTSNFRAC;
  int    IMJD_HOST  = GENSPEC.IMJD_HOST ;
  int    NBLAM      = INPUTS_SPECTRO.NBIN_LAM ;
  bool   ALLOW_HOST_ZEROFLUX = true; // for host

  int ilam, ilam2, NOPT=0 ;
  bool   IS_HOST_ZEROFLUX = false;
  double FLAM_PEAK, FLAM_HOST, FLAM_TOT, FLAM_SN ;
  double arg, MAGSHIFT, SCALE_FLAM_HOST=1.0 ;
  double FSUM_PEAK, FSUM_HOST, LAMAVG, LAMMIN, LAMMAX, LAMBIN ;
  char fnam[] = "GENSPEC_HOST_CONTAMINATION" ;

  // ------------- BEGIN --------------

  if ( IS_HOST )           { return; }
  if ( HOSTFRAC < 1.0E-8 && HOSTSNFRAC < 1.0E-8 ) { return; }

  // check that imjd=0 is indeed a host spectrum to add
  IS_HOST = GENSPEC.IS_HOST[IMJD_HOST];
  if ( !IS_HOST ) {
    sprintf(c1err,"IMJD=%d is not a host spectrum", IMJD_HOST);
    sprintf(c2err,"Cannot add host contamination to SN spectrum.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( HOSTFRAC > .001 ) {
    SCALE_FLAM_HOST = HOSTFRAC ; 
    NOPT++ ;
  }

  if ( HOSTSNFRAC > 0.001 ) {
    FSUM_PEAK = FSUM_HOST = 0.0 ;
    for(ilam=0; ilam < NBLAM; ilam++ ) {
      FLAM_PEAK   = GENSPEC.GENFLUX_PEAK[ilam];
      FLAM_HOST   = GENSPEC.GENFLUX_LIST[IMJD_HOST][ilam];

      if ( isnan(FLAM_HOST) ) {
        print_preAbort_banner(fnam);
        for(ilam2=ilam-2; ilam2 <= ilam+2; ilam2++ ) {
          if ( ilam2 >= 0 && ilam2 < NBLAM ) {
	    LAMAVG      = INPUTS_SPECTRO.LAMAVG_LIST[ilam2] ;
            printf("\t FLAM_HOST[%4d] = %le  at LAM=%.1f\n",
                   ilam2, GENSPEC.GENFLUX_LIST[IMJD_HOST][ilam2], LAMAVG );
          }
        }
        LAMAVG      = INPUTS_SPECTRO.LAMAVG_LIST[ilam] ;
        sprintf(c1err,"FLAM_HOST=NaN at LAM=%.1f A (ilam=%d), IDSPEC=%d",
                LAMAVG, ilam, HOSTSPEC.IDSPECDATA );
        sprintf(c2err,"Check SPECDATA.");
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

      LAMMIN      = INPUTS_SPECTRO.LAMMIN_LIST[ilam] ; 
      LAMMAX      = INPUTS_SPECTRO.LAMMAX_LIST[ilam] ;    
      LAMBIN      = LAMMAX - LAMMIN ;
      FSUM_PEAK  += (FLAM_PEAK*LAMBIN);
      FSUM_HOST  += (FLAM_HOST*LAMBIN);
    }

    IS_HOST_ZEROFLUX = FSUM_HOST < 1.0E-40 ;
    if ( IS_HOST_ZEROFLUX && !ALLOW_HOST_ZEROFLUX ) {
      print_preAbort_banner(fnam);
      double MJD = GENSPEC.MJD_LIST[imjd];
      printf("   imjd=%d  IMJD_HOST=%d \n", imjd, IMJD_HOST);
      printf("   zCMB = %.3f \n", GENLC.REDSHIFT_CMB);
      printf("   MJD  = %.3f, TOBS=%.2f\n", MJD, MJD-GENLC.PEAKMJD);
      printf("   GALID = %lld \n", SNHOSTGAL.GALID );
      sprintf(c1err,"Cannot implement HOSTSNFRAC=%f because",
              HOSTSNFRAC );
      sprintf(c2err,"FSUM_HOST = %le \n", FSUM_HOST );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    // HOSTSNFRAC \equiv (SCALE*FSUM_HOST/FSUM_PEAK)
    if ( !IS_HOST_ZEROFLUX ) 
      { SCALE_FLAM_HOST = HOSTSNFRAC * FSUM_PEAK / FSUM_HOST ; }

    NOPT++ ;
  }

  if ( NOPT > 1 ) {
    sprintf(c1err,"Cannot use both HOSTFRAC and HOSTSNFRAC.");
    sprintf(c1err,"Pick one (and see manual).");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // - - - - - -
  for(ilam=0; ilam < NBLAM; ilam++ ) {
    FLAM_SN   = GENSPEC.GENFLUX_LIST[imjd][ilam];
    FLAM_HOST = GENSPEC.GENFLUX_LIST[IMJD_HOST][ilam];
    FLAM_TOT  = FLAM_SN + (FLAM_HOST*SCALE_FLAM_HOST);
    
    GENSPEC.GENFLUX_LIST[imjd][ilam] = FLAM_TOT ;

    arg      =  FLAM_TOT/FLAM_SN ;
    MAGSHIFT = -2.5*log10(arg);
    GENSPEC.GENMAG_LIST[imjd][ilam] += MAGSHIFT ;
  } // end ilam

  return;

} // end GENSPEC_HOST_CONTAMINATION

// *************************************************
void GENSPEC_TEXPOSE_TAKE_SPECTRUM(int imjd) {

  // Compute TEXPOSE from requested SNR.
  // For synthetic filters, update ZPT and SKYSIG.
  // May 27 2020: check option to extrapolate TEXPOSE beyond defined range
  // Dec 09 2020: tol_converge -> 0.03 (was 0.02) to fix barely-missed
  //              convergence for SNLS spectrum.
  //

  bool LDMP       = (GENLC.CID == INPUTS.TAKE_SPECTRUM_DUMPCID );
  int  OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;  
  bool DO_TEXTRAP = ( (OPTMASK & SPECTROGRAPH_OPTMASK_TEXTRAP)>0 );

  int  NBLAM            = INPUTS_SPECTRO.NBIN_LAM ;  
  int  INDX             = GENSPEC.INDEX_TAKE_SPECTRUM[imjd] ;
  int  OPT_FRAME_EPOCH  = INPUTS.TAKE_SPECTRUM[INDX].OPT_FRAME_EPOCH ;
  int  OPT_FRAME_LAMBDA = INPUTS.TAKE_SPECTRUM[INDX].OPT_FRAME_LAMBDA ;
  int  OPT_TEXPOSE      = INPUTS.TAKE_SPECTRUM[INDX].OPT_TEXPOSE ;
  char *FRAME_EPOCH     = INPUTS.TAKE_SPECTRUM[INDX].EPOCH_FRAME ;
  double LAMMIN         = INPUTS.TAKE_SPECTRUM[INDX].SNR_LAMRANGE[0];
  double LAMMAX         = INPUTS.TAKE_SPECTRUM[INDX].SNR_LAMRANGE[1];

  double MJD         = GENSPEC.MJD_LIST[imjd]; 
  double TEXPOSE_MIN = INPUTS_SPECTRO.TEXPOSE_MIN ;
  double TEXPOSE_MAX = INPUTS_SPECTRO.TEXPOSE_MAX ;
  double SCALE       = 
    (double)INPUTS.TAKE_SPECTRUM_TEMPLATE_TEXPOSE_SCALE ; 
  GENPOLY_DEF *GENZPOLY_SNR = 
    &INPUTS.TAKE_SPECTRUM[INDX].GENZPOLY_SNR;

  bool IS_NEARPEAK   = ( imjd == GENSPEC.IMJD_NEARPEAK ); // SN near peak
  bool IS_HOST       = GENSPEC.IS_HOST[imjd];             // host

  double z           = GENLC.REDSHIFT_CMB ;  
  double z1          = 1.0 + z ;
  double TOBS        = MJD - GENLC.PEAKMJD ;
  double TREST       = TOBS/z1;

  double LAMMIN_OBS, LAMMAX_OBS ;
  double SNR=0.0, SNR_REQUEST, SNR_RATIO, ZPT, PSFSIG;
  double TEXPOSE_REQUEST, TEXPOSE, TEXPOSE_T, SKYSIG, SKYSIG_T;
  char fnam[] = "GENSPEC_TEXPOSE_TAKE_SPECTRUM" ;
  
  // ------------ BEGIN --------------

  if ( INDX < 0 )         { return ; } // not from TAKE_SPECTRUM key
  if ( OPT_TEXPOSE != 2 ) { return ; } // not SNR option

  SNR_REQUEST = eval_GENPOLY(z,GENZPOLY_SNR,fnam);
  GENSPEC.SNR_REQUEST_LIST[imjd] = SNR_REQUEST ;

  // extract min/max wavelength to determine SNR
  LAMMIN_OBS = GENSPEC.LAMOBS_SNR_LIST[imjd][0] ;
  LAMMAX_OBS = GENSPEC.LAMOBS_SNR_LIST[imjd][1] ;

  // get ILAM_MIN & ILAM_MAX to speed-up repeated calls
  // to GENSPEC_OBSFLUX_INIT  
  int ilam, ILAM_MIN=99999,  ILAM_MAX=-9;  
  for(ilam=0; ilam < NBLAM; ilam++ ) {
    if ( INPUTS_SPECTRO.LAMMAX_LIST[ilam] < LAMMIN_OBS ) { continue; }
    if ( INPUTS_SPECTRO.LAMMIN_LIST[ilam] > LAMMAX_OBS ) { continue; }
    if ( ILAM_MIN > 99990 ) { ILAM_MIN=ilam; }
    ILAM_MAX = ilam;
  }

  // set flag to NOT apply noise
  INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK |= SPECTROGRAPH_OPTMASK_noNOISE ;

  // start with TEXPOSE_MIN and MAX
  double SNR_MIN, SNR_MAX ;
  GENSPEC.TEXPOSE_LIST[imjd] = TEXPOSE_MIN ;
  GENSPEC_OBSFLUX_INIT(imjd,ILAM_MIN,ILAM_MAX);
  SNR_MIN = GENSPEC_SMEAR(imjd,LAMMIN_OBS,LAMMAX_OBS);

  GENSPEC.TEXPOSE_LIST[imjd] = TEXPOSE_MAX ;
  GENSPEC_OBSFLUX_INIT(imjd,ILAM_MIN,ILAM_MAX);
  SNR_MAX = GENSPEC_SMEAR(imjd,LAMMIN_OBS,LAMMAX_OBS);

  if ( LDMP ) {
    printf("\n xxx --------- %s DUMP for CID=%d ----------- \n", 
	   fnam, GENLC.CID );
    printf(" xxx imjd=%d  MJD=%.3f \n", imjd, MJD);
    printf(" xxx SNR_REQUEST = %.1f , z=%.3f  TOBS=%.2f  TREST=%.2f\n", 
	   SNR_REQUEST, z, TOBS, TREST);

    printf(" xxx LAMOBS = %.1f to %.1f \n", 
	   LAMMIN_OBS, LAMMAX_OBS);
    printf(" xxx LAMRANGE(stored) = %.1f to %.1f \n",
	   LAMMIN, LAMMAX);
    printf(" xxx SNR_MIN(%d sec)=%.3f, SNR_MAX(%d sec)=%.3f \n",
	    (int)TEXPOSE_MIN, SNR_MIN, (int)TEXPOSE_MAX, SNR_MAX);
    printf(" xxx OPT_FRAME_EPOCH=%d(%s)  OPT_FRAME_LAMBDA=%d \n",
	   OPT_FRAME_EPOCH, FRAME_EPOCH, OPT_FRAME_LAMBDA );
    fflush(stdout);
  }

  if ( SNR_REQUEST <= SNR_MIN ) {
    if ( DO_TEXTRAP ) { 
      SNR_RATIO                      = SNR_REQUEST/SNR_MIN;
      TEXPOSE_REQUEST                = TEXPOSE_MIN * (SNR_RATIO*SNR_RATIO);
      GENSPEC.TEXPOSE_LIST[imjd]     = TEXPOSE_REQUEST ;
      GENSPEC.SNR_COMPUTE_LIST[imjd] = SNR_REQUEST ;
    }
    else {
      GENSPEC.TEXPOSE_LIST[imjd]     = TEXPOSE_MIN ;
      GENSPEC.SNR_COMPUTE_LIST[imjd] = SNR_MIN ;
    }
    goto CLEANUP ;
  }


  if ( SNR_REQUEST >= SNR_MAX ) {
    if ( DO_TEXTRAP ) {  // option to extrapolate TEXPOSE
      SNR_RATIO                      = SNR_REQUEST/SNR_MAX;
      TEXPOSE_REQUEST                = TEXPOSE_MAX * (SNR_RATIO*SNR_RATIO);
      GENSPEC.TEXPOSE_LIST[imjd]     = TEXPOSE_REQUEST ;
      GENSPEC.SNR_COMPUTE_LIST[imjd] = SNR_REQUEST ;
    }
    else {
      GENSPEC.TEXPOSE_LIST[imjd]     = TEXPOSE_MAX ;
      GENSPEC.SNR_COMPUTE_LIST[imjd] = SNR_MAX ;
    }

    /*
    printf(" xxx TEXPOSE_TAKE_SPEC: DO_TEXTRAP=%d "
	   "SNR(MAX,REQ)=%5.1f,%5.1f  TEXPOSE(MAX,REQ)=%6.0f,%6.0f \n",
	   DO_TEXTRAP,	SNR_MAX, SNR_REQUEST,
	   TEXPOSE_MAX, TEXPOSE_REQUEST );
    */

    goto CLEANUP ;
  }

  // assume SNR is linear function of sqrt(TEXPOSE);
  // iterate until convergence.
  //    SNR ~ SNR0 + (SNR1-SNR0) * sqrt[(T-T0)/(T1-T0)]
  //    T -> T0 + (T1-T0) * [ (SNR-SNR0)/(SNR1-SNR0) ]^2 
  // 

  double T0=TEXPOSE_MIN, T1=TEXPOSE_MAX, TDIF_LAST=T1-T0;
  double SNR0=SNR_MIN, SNR1=SNR_MAX;
  double argSNR, last_tol, last_TEXPOSE, sgn=0.0 ;
  double tol=99.0 , tol_converge = 0.03 ;
  int    NITER=0, MAXITER=40, FLAG_TEXPOSE=0 ;  

  TEXPOSE = 0.0;

  while ( fabs(tol) > tol_converge ) {

    if ( NITER==0 ) 
      { FLAG_TEXPOSE=0 ; }
    else if ( sgn == 0.0 ) 
      { FLAG_TEXPOSE=1 ; }
    else
      { FLAG_TEXPOSE=2 ; }

    last_tol = tol;
    last_TEXPOSE = TEXPOSE ;

    if ( FLAG_TEXPOSE == 0 ) {
      if ( IS_NEARPEAK ) {
	// first epoch is always closest to peak
	argSNR  = (SNR_REQUEST - SNR0)/(SNR1 - SNR0);
	TEXPOSE = T0 + TDIF_LAST * (argSNR*argSNR) ;
	TEXPOSE *= 2.0 ; // guess from added template noise
      }
      else  { 
	// best guess is that TEXPOSE is similar to TEXPOSE_TEMPLATE
	TEXPOSE = GENSPEC.TEXPOSE_TEMPLATE*1.2 ; 
      }
    }
    else if ( FLAG_TEXPOSE == 1 ) {
      // keep scaling by SNR-ratio until we go too far
      argSNR  = (SNR_REQUEST/SNR) ;
      TEXPOSE = TEXPOSE * argSNR ;
    }
    else if ( FLAG_TEXPOSE == 2 ) {
      TEXPOSE = last_TEXPOSE * ( 1.0 - 0.02*sgn) ;
    }

    // avoid slipping outside defined range
    if ( TEXPOSE <= TEXPOSE_MIN ) { TEXPOSE = TEXPOSE_MIN ; }
    if ( TEXPOSE >= TEXPOSE_MAX ) { TEXPOSE = TEXPOSE_MAX ; }
    
    // for GENSPEC_SMEAR, check option to set global template 
    // exposure time by scaling T_expose for epoch nearest peak.
    if ( IS_NEARPEAK && SCALE > 0.01 )  { 
      TEXPOSE_T = TEXPOSE * SCALE ; 
      if ( TEXPOSE_T <= TEXPOSE_MIN ) { TEXPOSE_T = TEXPOSE_MIN; }
      if ( TEXPOSE_T >= TEXPOSE_MAX ) { TEXPOSE_T = TEXPOSE_MAX; }
      GENSPEC.TEXPOSE_TEMPLATE = TEXPOSE_T ;
    }
    
    // set global used in GENSPEC_SMEAR
    GENSPEC.TEXPOSE_LIST[imjd] = TEXPOSE ;     
    GENSPEC_OBSFLUX_INIT(imjd,ILAM_MIN,ILAM_MAX);
    SNR     = GENSPEC_SMEAR(imjd,LAMMIN_OBS,LAMMAX_OBS);
    tol     = (SNR/SNR_REQUEST) - 1.0 ;

    // save SNR_COMPUTE for output
    GENSPEC.SNR_COMPUTE_LIST[imjd] = SNR ;   

    if ( LDMP ) {
      printf(" xxx iter=%2d: Texpose(S/T)->%6.1f/%.1f  SNR=%6.2f  tol=%7.3f "
	     "(sgn=%4.1f)\n",
	     NITER,  TEXPOSE, GENSPEC.TEXPOSE_TEMPLATE, SNR, tol, sgn);
    }

    // if tolerance changes sign, then just add small amount
    // to TEXPOSE to avoid bouncing below/above SNR_REQUEST
    if ( tol/last_tol < 0.0 && NITER>0 && FLAG_TEXPOSE==1 ) {
      sgn = tol/fabs(tol);
      FLAG_TEXPOSE = 2 ;
    }

    NITER++;
    if ( NITER >= MAXITER ) {
      print_preAbort_banner(fnam);
      printf("\t SNR_MIN(%d sec)=%.1f \n", (int)TEXPOSE_MIN, SNR_MIN ) ;
      printf("\t SNR_MAX(%d sec)=%.1f \n", (int)TEXPOSE_MAX, SNR_MAX ) ;

      sprintf(c1err,"Could not converge after NITER=%d (CID=%d)", 
	      NITER, GENLC.CID );
      sprintf(c2err,"SNR_REQUEST=%.1f  z=%.3f  TOBS=%.2f",
	      SNR_REQUEST, z, TOBS);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }
    

 CLEANUP:

  if ( LDMP ) {
    printf(" xxx FINAL TEXPOSE(S,T)=%.2f,%.2f   SNR_COMPUTE=%.2f  \n",
	   GENSPEC.TEXPOSE_LIST[imjd] ,
	   GENSPEC.TEXPOSE_TEMPLATE,
	   GENSPEC.SNR_COMPUTE_LIST[imjd] );
    fflush(stdout);
  }

  INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK -= SPECTROGRAPH_OPTMASK_noNOISE ;
  GENSPEC_OBSFLUX_INIT(imjd,ILAM_MIN,ILAM_MAX);

  // ------------------------------------------------------
  // Jun 19 2017: 
  // use computed TEXPOSE to update ZPT and SKYSIG for SYNTHETIC filters.

  if ( GENLC.NFILTDEF_SPECTROGRAPH == 0 ) { return ; }
  
  // first find SIMLIB epoch matching this MJD
  int ep, ifilt, ifilt_obs;
  int ISTORE_SIMLIB_GEN=-9, ISTORE_SIMLIB_RAW=-9, ISTORE ;
  double VAL_STORE[10];

  for(ep=0; ep < GENLC.NEPOCH; ep++ ) {

    if ( SIMLIB_OBS_GEN.INDX_TAKE_SPECTRUM[ep] != INDX ) { continue ; }

    ifilt_obs  = SIMLIB_OBS_GEN.IFILT_OBS[ep] ;
    ifilt      = GENLC.IFILTINV_SPECTROGRAPH[ifilt_obs] ; 
    if ( ifilt_obs < 0 ) { continue ; }

    ISTORE_SIMLIB_GEN = ep ;
    ISTORE_SIMLIB_RAW = SIMLIB_OBS_GEN.ISTORE_RAW[ep] ;
    TEXPOSE           = GENSPEC.TEXPOSE_LIST[imjd] ;
    TEXPOSE_T         = GENSPEC.TEXPOSE_TEMPLATE ;
    
    VAL_STORE[0] = MJD ;
    VAL_STORE[1] = TEXPOSE ;
    VAL_STORE[2] = (int)INDX ;
    
    get_SPECTROGRAPH_ZPTPSFSKY(ISTORE_SIMLIB_RAW, ifilt, 
			       TEXPOSE, TEXPOSE_T,   
			       &ZPT, &PSFSIG, &SKYSIG, &SKYSIG_T ); 
    
    ISTORE = ISTORE_SIMLIB_GEN ;
    SIMLIB_OBS_GEN.ZPTADU[ISTORE]   = ZPT ;
    SIMLIB_OBS_GEN.SKYSIG[ISTORE]   = SKYSIG ; 
    SIMLIB_OBS_GEN.CCDGAIN[ISTORE]  = 1.0 ;      // unity GAIN 
    SIMLIB_OBS_GEN.PSFSIG1[ISTORE]  = PSFSIG ;   // sigma1
    SIMLIB_OBS_GEN.NEA[ISTORE]      = NoiseEquivAperture(PSFSIG, 0.0, 0.0);
   
    // set optional template noise for passbands made from spectrograph
    if ( TEXPOSE_T > 0.001 ) {
      SIMLIB_TEMPLATE.USEFLAG |= 2 ;
      SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[ISTORE]    = SKYSIG_T ; 
      SIMLIB_OBS_GEN.TEMPLATE_ZPT[ISTORE]       = ZPT ;	
    }

    if ( LDMP ) { 
      printf(" xxx FINAL SKYSIG-%c (S,T) = %8.1f, %8.1f  [TOBS=%6.2f] \n",
	     FILTERSTRING[ifilt_obs], SKYSIG, SKYSIG_T, 
	     GENSPEC.TOBS_LIST[imjd] );  
      fflush(stdout);
    }
      
  }  // end ep loop over epochs


  return ;

} // end GENSPEC_TEXPOSE_TAKE_SPECTRUM

// *************************************************
double GENSPEC_SMEAR(int imjd, double LAMMIN, double LAMMAX ) {
  
  // apply noise to flux and smear over wavelength bins.
  // Return SNR over input wavelength range.
  //
  // Mar  6 2020: skip bin if SNR_true = 0
  // May 10 2021: 
  //  Major refactor to use SNR from measured wave bins, not true wave bin.
  //  Fix bug from v10_78 where GRAN_T no longer followed correlated option.

  int    NBLAM = INPUTS_SPECTRO.NBIN_LAM ;
  int    MEMD  = NBLAM * sizeof(double);

  int    ilam, ILAM_MIN=99999, ILAM_MAX=-9, NBLAM_USE=0 ;
  double GENFLUX, GENFLUXERR, GENFLUXERR_T, GENMAG, LAMAVG ;
  double *SNR_TRUE_LIST,   SNR_TRUE, *ERRFRAC_T_LIST, ERRFRAC_T ; 

  double  TEXPOSE_S  = GENSPEC.TEXPOSE_LIST[imjd] ;
  double  TEXPOSE_T  = GENSPEC.TEXPOSE_TEMPLATE ;
  double  SCALE_SNR  = INPUTS.SPECTROGRAPH_OPTIONS.SCALE_SNR ;
  double  SNR_SPEC ;

  int  OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;  
  bool ALLOW_TEXTRAP = ( (OPTMASK & SPECTROGRAPH_OPTMASK_TEXTRAP)>0 );
  bool IS_HOST       = GENSPEC.IS_HOST[imjd];
  char   fnam[] = "GENSPEC_SMEAR" ;

  // ------------- BEGIN -------------

  // - - - - - -

  SNR_TRUE_LIST   = (double*) malloc( MEMD ) ;
  ERRFRAC_T_LIST  = (double*) malloc( MEMD ) ;
 
  for(ilam=0; ilam < NBLAM; ilam++ ) {

    SNR_TRUE_LIST[ilam] = -9.0 ;

    LAMAVG = INPUTS_SPECTRO.LAMAVG_LIST[ilam] ;
    if ( LAMAVG < LAMMIN ) { continue; }
    if ( LAMAVG > LAMMAX ) { continue; }

    if ( ILAM_MIN > 99990 ) { ILAM_MIN=ilam; }
    ILAM_MAX = ilam;

    GENFLUX = GENSPEC.GENFLUX_LIST[imjd][ilam] ;
    GENMAG  = GENSPEC.GENMAG_LIST[imjd][ilam] ;

    // skip unphysical fluxes
    if ( GENFLUX <= 0.0  ) { continue ; }
    if ( GENMAG  > 600.0 ) { continue ; } // Mar 2019

    // get true SNR in this lambda bin
    SNR_TRUE =
      getSNR_spectrograph(ilam, TEXPOSE_S, TEXPOSE_T, ALLOW_TEXTRAP, GENMAG,
			  &ERRFRAC_T);  // template frac of error

    SNR_TRUE_LIST[ilam]  = SNR_TRUE ;
    ERRFRAC_T_LIST[ilam] = ERRFRAC_T ;

    // apply lambda smear to distribute GENFLUX over lambda bins 
    GENSPEC_LAMSMEAR(imjd, ilam, GENFLUX );

    NBLAM_USE++ ;

  } // end ilam  

  if ( NBLAM_USE == 0 ) { goto DONE ; }

  // - - - - - - - - - - - - - - 
  // after smearing flux in neighbor bins, loop again over wavelegth
  // and apply Poisson noise.

  double OBSFLUX, OBSFLUX_SMEAR, OBSFLUXERR, OBSFLUXERR_T, *GAURAN_T;
  double ERRSQ, SUM_FLUX, SUM_ERRSQ, SUM_ERR ;
  SUM_FLUX = SUM_ERRSQ = SNR_SPEC = 0.0 ;

  for(ilam = ILAM_MIN; ilam <= ILAM_MAX ; ilam++ ) {

    SNR_TRUE  = SNR_TRUE_LIST[ilam];
    ERRFRAC_T = ERRFRAC_T_LIST[ilam];

    if ( SNR_TRUE < 1.0E-18 ) { continue; } 
    if ( SCALE_SNR != 1.00 ) { SNR_TRUE *= SCALE_SNR ;  }

    OBSFLUX       = GENSPEC.OBSFLUX_LIST[imjd][ilam] ;
    OBSFLUXERR    = OBSFLUX / SNR_TRUE ;

    if ( OBSFLUXERR < 1.0E-50 ) { continue; }

    // compute random flucution of spectrograph flux
    // be careful with correlated template noise in each spectrum
    GAURAN_T = &GENSPEC.RANGauss_NOISE_TEMPLATE[ilam];
    OBSFLUX_SMEAR = 
      GENSPEC_OBSFLUX_RANSMEAR(imjd, OBSFLUXERR,ERRFRAC_T, GAURAN_T);
					     

    // obdate observed flux and store it
    OBSFLUX += OBSFLUX_SMEAR ;

    // store results
    ERRSQ                                 = OBSFLUXERR * OBSFLUXERR ;
    GENSPEC.OBSFLUX_LIST[imjd][ilam]      = OBSFLUX ;
    GENSPEC.OBSFLUXERR_LIST[imjd][ilam]   = OBSFLUXERR ;
    GENSPEC.OBSFLUXERRSQ_LIST[imjd][ilam] = ERRSQ ;

    GENSPEC.NBLAM_VALID[imjd]++ ; 
    SUM_FLUX    += GENSPEC.OBSFLUX_LIST[imjd][ilam] ;
    SUM_ERRSQ   += ERRSQ ;
  } // end ilam loop

  if ( SUM_ERRSQ > 0.0 ) {
    SUM_ERR  = sqrt(SUM_ERRSQ);
    SNR_SPEC = SUM_FLUX / SUM_ERR ;
  }

  /*
  printf("\t xxx %s: SNR = %le / %le = %le \n", 
	 fnam, SUM_FLUX, sqrt(SUM_ERRSQ), SNR_SPEC); fflush(stdout);
  */

 DONE:
  free(SNR_TRUE_LIST);
  free(ERRFRAC_T_LIST);

  return(SNR_SPEC) ;

} // end GENSPEC_SMEAR

// *********************************************
void  GENSPEC_LAMSMEAR(int imjd, int ilam, double GenFlux ) {

  // Use Gaussian lambda-resolution to smear flux(imjd,ilam)
  // over nearby lambda bins. Do NOT apply Poisson noise here.
  //
  // Inputs
  //  + imjd         = sparse MJD index for spectra
  //  + ilam         = wavelength index
  //  + GenFlux      = flux
  //
  //
  // Store smeared/observed flux in  arrays
  //   + GENSPEC.OBSFLUX_LIST 
  //
  // May 10 2021:
  //  + refactor/simplify to NOT apply Poisson noise. Noise is added
  //    later so that SNR properties are from smeared wave bin,
  //    and not from true wave bin.
  //

  int OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;
  int noNOISE    = ( OPTMASK & SPECTROGRAPH_OPTMASK_noNOISE    ) ;

  double NSIGLAM, LAMAVG, LAMSIGMA, LAMBIN, LAMSIG0, LAMSIG1;
  double GINT, tmp_GenFlux ;
  int    NBIN2, ilam2, ilam_tmp, NBLAM, NRAN, LDMP=0 ;
  char fnam[] = "GENSPEC_LAMSMEAR" ;

  // ----------- BEGIN ---------------

  NBLAM    = INPUTS_SPECTRO.NBIN_LAM ;

  LAMAVG   = INPUTS_SPECTRO.LAMAVG_LIST[ilam] ;
  LAMSIGMA = INPUTS_SPECTRO.LAMSIGMA_LIST[ilam] ;
  NSIGLAM  = INPUTS.SPECTROGRAPH_OPTIONS.NLAMSIGMA ;

  if ( noNOISE > 0  ) { LAMSIGMA = 0.0 ; }

  // for LAMBIN, avoid edge bin which can be artificially small
  // leading to excessively large NBIN2.
  ilam_tmp = ilam;  if ( ilam == NBLAM-1 ) { ilam_tmp = NBLAM-2; }
  LAMBIN   = INPUTS_SPECTRO.LAMBIN_LIST[ilam_tmp] ;
  NBIN2    = (int)(NSIGLAM*LAMSIGMA/LAMBIN + 0.5) ;  
  NRAN     = 0 ;

  /*
  if( ilam > NBLAM-6 && imjd==0 ) {
    printf(" xxx %s: ilam=%d LAMAVG=%7.1f  NBIN2=%d \n",
	   fnam, ilam, LAMAVG, NBIN2);
  }
  //xxxxxxx */

  // loop over neighbor bins to smear flux over lambda bins
  for(ilam2=ilam-NBIN2; ilam2 <= ilam+NBIN2; ilam2++ ) {
    if ( ilam2 <  0     ) { continue ; }
    if ( ilam2 >= NBLAM ) { continue ; }
    
    // don't bother loading extended lambda bins for lam-res
    if ( INPUTS_SPECTRO.ISLAM_EXTEND_LIST[ilam2] ) { continue; }    

    GINT = 1.0 ;
    if ( LAMSIGMA > 0.0 ) {
      LAMSIG0 = (INPUTS_SPECTRO.LAMMIN_LIST[ilam2]-LAMAVG)/LAMSIGMA ;
      LAMSIG1 = (INPUTS_SPECTRO.LAMMAX_LIST[ilam2]-LAMAVG)/LAMSIGMA ;
      GINT = GaussIntegral(LAMSIG0,LAMSIG1); 
    }
    else if ( ilam != ilam2 ) // LAMSIGMA=0; dump all flux in same bin
      { GINT = 0.0 ; }

    // true flux in this lambda bin
    tmp_GenFlux      = ( GINT*GenFlux ) ;

    // increment sum of obsFlux 
    // OBSFLUX_LIST is modified later to include Poisson noise;
    // GENFLUX_LAMSMEAR_LIST is not modified.
    GENSPEC.OBSFLUX_LIST[imjd][ilam2]          += tmp_GenFlux ;  
    GENSPEC.GENFLUX_LAMSMEAR_LIST[imjd][ilam2] += tmp_GenFlux ;

    NRAN++ ;

  } // end ilam2
    
  // - - - - - -
  if ( NRAN >= MXLAMSMEAR_SPECTROGRAPH ) {
    print_preAbort_banner(fnam);    
    printf("\t NSIGLAM  = %f \n", NSIGLAM);
    printf("\t LAMSIGMA = %f \n", LAMSIGMA );
    printf("\t LAMBIN   = %f \n", LAMBIN);
    printf("\t NBIN2    = %d \n", NBIN2 );
    sprintf(c1err,"NLAMSMEAR = %d exceeds bound of %d",
	    NRAN, MXLAMSMEAR_SPECTROGRAPH );
    sprintf(c2err,"ilam=%d LAMAVG=%.2f", ilam, LAMAVG );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  return ;

} // end GENSPEC_LAMSMEAR

// *************************************************
double GENSPEC_OBSFLUX_RANSMEAR(int imjd, double OBSFLUXERR, double ERRFRAC_T, 
				double *GAURAN_T) {

  // Created May 10 2021
  // Compute and return random fluctuation for spectrograph flux.
  // Inputs:
  //   + OBSFLUXERR = uncertainty 
  //   + ERRFRAC_T  = fraction of error from template -> correlated
  //
  // For correlated template noise, *GAURAN_T is set on first
  // MJD nearest peak; then re-used for other SN spectra.
  // If ERRFRAC_T = 0, *GAURAN_T is set to zero.

  int OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;
  int onlyTNOISE = ( OPTMASK & SPECTROGRAPH_OPTMASK_onlyTNOISE ) ;
  int noTNOISE   = ( OPTMASK & SPECTROGRAPH_OPTMASK_noTEMPLATE ) ;
  int noNOISE    = ( OPTMASK & SPECTROGRAPH_OPTMASK_noNOISE    ) ;

  int NSTREAM      = GENRAN_INFO.NSTREAM ;
  int ISTREAM_RAN  = ISTREAM_RANDOM_SPECTROGRAPH ;
  int ILIST_RAN    = ILIST_RANDOM_SPECTROGRAPH ; // mark obsolete, Jun 4 2020

  int IMJD_NEARPEAK = GENSPEC.IMJD_NEARPEAK ;

  double RanFlux_S, RanFlux_T, GRAN_S, GRAN_T ;
  double FluxErr_S, FluxErr_T;

  double OBSFLUX_RANSMEAR = 0.0 ;
  char fnam[] = "GENSPEC_OBSFLUX_RANSMEAR";

  // ---------- BEGIN ----------

  if ( noNOISE ) {
    RanFlux_S = RanFlux_T = 0.0 ; 
  }
  else {

    FluxErr_T = ERRFRAC_T * OBSFLUXERR;
    FluxErr_S = sqrt(OBSFLUXERR*OBSFLUXERR - FluxErr_T*FluxErr_T);

    if ( ERRFRAC_T < 1.0E-7 ) {
      *GAURAN_T = 0.0; 
    }
    else if ( imjd == IMJD_NEARPEAK &&  FluxErr_T > 0.0 ) { 
      if ( NSTREAM == 2 ) 
	{ *GAURAN_T = unix_getRan_Gauss(ISTREAM_RAN); }
      else
	{ *GAURAN_T = getRan_Gauss(ILIST_RAN); }
    }
    GRAN_T    = *GAURAN_T ;

    if ( NSTREAM ==  2 ) 
      { GRAN_S = unix_getRan_Gauss(ISTREAM_RAN); }
    else
      { GRAN_S = getRan_Gauss(ILIST_RAN); }

    // random noise from search 
    RanFlux_S = FluxErr_S * GRAN_S ;
    if ( onlyTNOISE ) { RanFlux_S = 0.0 ; }
    
    // correlated random noise from template
    RanFlux_T = FluxErr_T * GRAN_T ;
    if ( noTNOISE ) { RanFlux_T = 0.0 ; }
    
    /* xxxxxxxx
    printf(" xxx ERRFRAC_T=%f   GRAN_T=%f RanFlux_T = %f \n", 
	   ERRFRAC_T, GRAN_T, RanFlux_T);
    debugexit(fnam); 
    xxxxxx */
  }

  // - - - - - -
  // add fluctuation from search (S) and correlated template (T)
  OBSFLUX_RANSMEAR = RanFlux_S + RanFlux_T ; 


  return OBSFLUX_RANSMEAR  ;

} // end GENSPEC_OBSFLUX_RANSMEAR

// *************************************************
double GENSPEC_SMEAR_LEGACY(int imjd, double LAMMIN, double LAMMAX ) {
  
  // apply noise to flux and smear over wavelength bins.
  // Return SNR over input wavelength range.
  //
  // Mat 6 2020: skip bin if SNR_true = 0

  int    NBLAM = INPUTS_SPECTRO.NBIN_LAM ;
  int    ilam, ILAM_MIN=99999, ILAM_MAX=-9, NBLAM_USE=0 ;
  double GENFLUX, GENFLUXERR, GENFLUXERR_T, GENMAG, SNR_true ;
  double ERRFRAC_T, LAMAVG ; 

  double  TEXPOSE_S  = GENSPEC.TEXPOSE_LIST[imjd] ;
  double  TEXPOSE_T  = GENSPEC.TEXPOSE_TEMPLATE ;
  double  SCALE_SNR  = INPUTS.SPECTROGRAPH_OPTIONS.SCALE_SNR ;
  double  SNR_SPEC ;

  int  OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;  
  bool ALLOW_TEXTRAP = ( (OPTMASK & SPECTROGRAPH_OPTMASK_TEXTRAP)>0 );
  //  char   fnam[] = "GENSPEC_SMEAR_LEGACY" ;

  // ------------- BEGIN -------------

  // ****** MARK OBSOLETE *******

  for(ilam=0; ilam < NBLAM; ilam++ ) {

    LAMAVG = INPUTS_SPECTRO.LAMAVG_LIST[ilam] ;
    if ( LAMAVG < LAMMIN ) { continue; }
    if ( LAMAVG > LAMMAX ) { continue; }

    if ( ILAM_MIN > 99990 ) { ILAM_MIN=ilam; }
    ILAM_MAX = ilam;

    GENFLUX = GENSPEC.GENFLUX_LIST[imjd][ilam] ;
    GENMAG  = GENSPEC.GENMAG_LIST[imjd][ilam] ;

    // skip unphysical fluxes
    if ( GENFLUX <= 0.0  ) { continue ; }
    if ( GENMAG  > 600.0 ) { continue ; } // Mar 2019

    // get true SNR in this lambda bin
    SNR_true = getSNR_spectrograph(ilam, TEXPOSE_S, TEXPOSE_T, ALLOW_TEXTRAP,
				   GENMAG, 
				   &ERRFRAC_T);  // template frac of error

    // ****** MARK OBSOLETE *******

    if ( SNR_true < 1.0E-18 ) { continue; } // May 2020

    if ( SCALE_SNR != 1.00 ) { SNR_true *= SCALE_SNR ;  }
    
    GENFLUXERR    = GENFLUX/SNR_true ;       // sigma on FLUXGEN
    GENFLUXERR_T  = GENFLUXERR * ERRFRAC_T;  // template contribution

    // apply lambda smear to distribute GENFLUX over lambda bins 
    GENSPEC_LAMSMEAR_LEGACY(imjd, ilam, GENFLUX, GENFLUXERR, GENFLUXERR_T );

    NBLAM_USE++ ;

  } // end ilam  

  if ( NBLAM_USE == 0 ) { return(0.0); }

  // ****** MARK OBSOLETE *******

  // - - - - - - - - - - - - - - 
  // convert sum(ERRSQ) -> ERR in each lambda bin and
  // count number of wavelength bins with flux;
  double ERRSQ, SUM_FLUX, SUM_ERRSQ, SUM_ERR ;
  SUM_FLUX = SUM_ERRSQ = SNR_SPEC = 0.0 ;

  for(ilam=ILAM_MIN; ilam <=ILAM_MAX ; ilam++ ) {
    ERRSQ = GENSPEC.OBSFLUXERRSQ_LIST[imjd][ilam] ;
    if ( ERRSQ <= 1.0E-100 ) { continue ; }

    GENSPEC.OBSFLUXERR_LIST[imjd][ilam] = sqrt(ERRSQ);
    GENSPEC.NBLAM_VALID[imjd]++ ; 

    SUM_FLUX    += GENSPEC.OBSFLUX_LIST[imjd][ilam] ;
    SUM_ERRSQ   += ERRSQ ;
  }

  // ****** MARK OBSOLETE *******

  if ( SUM_ERRSQ > 0.0 ) {
    SUM_ERR  = sqrt(SUM_ERRSQ);
    SNR_SPEC = SUM_FLUX / SUM_ERR ;
  }

  /*
  printf("\t xxx %s: SNR = %le / %le = %le \n", 
	 fnam, SUM_FLUX, sqrt(SUM_ERRSQ), SNR_SPEC); fflush(stdout);
  */

  return(SNR_SPEC) ;

} // end GENSPEC_SMEAR_LEGACY

// *********************************************
void  GENSPEC_LAMSMEAR_LEGACY(int imjd, int ilam, double GenFlux, 
			      double GenFluxErr, double GenFluxErr_T ) {

  // Use Gaussian lambea-resolution to smear flux(imjd,ilam)
  // over nearby lambda bins.  In each lambda bin, use
  // GenFluxErr to add Poisson noise.
  // Inputs
  //  + imjd         = sparse MJD index for spectra
  //  + ilam         = wavelength index
  //  + GenFlux      = flux
  //  + GenFluxErr   = total flux error, including template noise
  //  + GenFluxErr_T = template noise contribution
  //
  //
  // Store smeared/observed flux in  arrays
  //   + GENSPEC.OBSFLUX_LIST 
  //   + GENSPEC.OBSFLUXERRSQ_LIST 
  //
  // Jan 17 2018: use ILIST_RANDOM_SPECTROGRAPH
  // Oct 25 2019: fix bug setting GRAN_T if there is no template.
  // Jun 01 2020: move NRAN abort outside loop with more info

  int OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;
  int onlyTNOISE = ( OPTMASK & SPECTROGRAPH_OPTMASK_onlyTNOISE ) ;
  int noTNOISE   = ( OPTMASK & SPECTROGRAPH_OPTMASK_noTEMPLATE ) ;
  int noNOISE    = ( OPTMASK & SPECTROGRAPH_OPTMASK_noNOISE    ) ;

  int NSTREAM      = GENRAN_INFO.NSTREAM ;
  int ISTREAM_RAN  = ISTREAM_RANDOM_SPECTROGRAPH ;
  int ILIST_RAN    = ILIST_RANDOM_SPECTROGRAPH ; // mark obsolete, Jun 4 2020

  double NSIGLAM, LAMAVG, LAMSIGMA, LAMBIN, LAMSIG0, LAMSIG1;
  double GINT, SUM_GINT, GINT_SQRT, GRAN_S, GRAN_T ;
  double tmp_GenFlux, tmp_GenFluxErr, tmp_GenFluxErr_S, tmp_GenFluxErr_T ;
  double tmp_RanFlux_S, tmp_RanFlux_T, RANGauss_NOISE_TEMPLATE;
  double GenFluxErr_S, OBSFLUX, OBSFLUXERR ;
  int    NBIN2, ilam2, ilam_tmp, NBLAM, NRAN, LDMP=0 ;
  char fnam[] = "GENSPEC_LAMSMEAR_LEGACY" ;

  // ----------- BEGIN ---------------

  NBLAM    = INPUTS_SPECTRO.NBIN_LAM ;

  LAMAVG   = INPUTS_SPECTRO.LAMAVG_LIST[ilam] ;
  LAMSIGMA = INPUTS_SPECTRO.LAMSIGMA_LIST[ilam] ;
  NSIGLAM  = INPUTS.SPECTROGRAPH_OPTIONS.NLAMSIGMA ;

  if ( noNOISE > 0  ) { LAMSIGMA = 0.0 ; }

  // get search error by subtracting template contribution
  GenFluxErr_S = sqrt( GenFluxErr*GenFluxErr - GenFluxErr_T*GenFluxErr_T);

  // for LAMBIN, avoid edge bin which can be artificially small
  // leading to excessively large NBIN2.
  ilam_tmp = ilam;  if ( ilam == NBLAM-1 ) { ilam_tmp = NBLAM-2; }
  LAMBIN   = INPUTS_SPECTRO.LAMBIN_LIST[ilam_tmp] ;
  NBIN2    = (int)(NSIGLAM*LAMSIGMA/LAMBIN + 0.5) ;  
  SUM_GINT = 0.0 ; 
  NRAN     = 0 ;

  /*
  if( ilam > NBLAM-6 && imjd==0 ) {
    printf(" xxx %s: ilam=%d LAMAVG=%7.1f  NBIN2=%d \n",
	   fnam, ilam, LAMAVG, NBIN2);
  }
  //xxxxxxx */

  // loop over neighbor bins to smear flux over lambda bins
  for(ilam2=ilam-NBIN2; ilam2 <= ilam+NBIN2; ilam2++ ) {
    if ( ilam2 <  0     ) { continue ; }
    if ( ilam2 >= NBLAM ) { continue ; }
    
    // don't bother loading extended lambda bins for lam-res
    if ( INPUTS_SPECTRO.ISLAM_EXTEND_LIST[ilam2] ) { continue; }    

    GINT = 1.0 ;
    if ( LAMSIGMA > 0.0 ) {
      LAMSIG0 = (INPUTS_SPECTRO.LAMMIN_LIST[ilam2]-LAMAVG)/LAMSIGMA ;
      LAMSIG1 = (INPUTS_SPECTRO.LAMMAX_LIST[ilam2]-LAMAVG)/LAMSIGMA ;
      GINT = GaussIntegral(LAMSIG0,LAMSIG1); 
    }
    else if ( ilam != ilam2 ) // LAMSIGMA=0; dump all flux in same bin
      { GINT = 0.0 ; }

    SUM_GINT += GINT ; // for debug only; should be 1 except near edges    
    GINT_SQRT = sqrt(GINT);

    // true flux in this lambda bin
    tmp_GenFlux      = ( GINT*GenFlux ) ;
    tmp_GenFluxErr   = GenFluxErr   * GINT_SQRT; // error scaled to flux-frac
    tmp_GenFluxErr_S = GenFluxErr_S * GINT_SQRT; 
    tmp_GenFluxErr_T = GenFluxErr_T * GINT_SQRT; 

    if ( noNOISE > 0 ) 
      { tmp_RanFlux_S = tmp_RanFlux_T = 0.0 ; }
    else {

      GRAN_S = GRAN_T = RANGauss_NOISE_TEMPLATE = 0.0 ;
      if ( GENSPEC.NMJD_PROC==0 && tmp_GenFluxErr_T > 0.0 ) { 
	if ( NSTREAM == 2 ) 
	  { RANGauss_NOISE_TEMPLATE = unix_getRan_Gauss(ISTREAM_RAN); }
	else
	  { RANGauss_NOISE_TEMPLATE = getRan_Gauss(ILIST_RAN); }
      }

      if ( NSTREAM ==  2 ) 
	{ GRAN_S = unix_getRan_Gauss(ISTREAM_RAN); }
      else
	{ GRAN_S = getRan_Gauss(ILIST_RAN); }

      GRAN_T = RANGauss_NOISE_TEMPLATE ;

      // random noise from search 
      tmp_RanFlux_S = tmp_GenFluxErr_S * GRAN_S ;
      if ( onlyTNOISE ) { tmp_RanFlux_S = 0.0 ; }
      
      // correlated random noise from template
      tmp_RanFlux_T = tmp_GenFluxErr_T * GRAN_T ;
      if ( noTNOISE ) { tmp_RanFlux_T = 0.0 ; }
    }

    
    LDMP = (ilam > NBLAM-6 && imjd == -10 );
    if ( LDMP ) {
      printf(" xxx ilam=%3d(%d)  imjd=%d  noise(S,T) = %10.3le , %10.3le\n",
	     ilam, ilam2, imjd, tmp_RanFlux_S, tmp_RanFlux_T );
      //      printf(" xxx \t (%f, %f) \n", GenFluxErr, GenFluxErr_T  );
    }
    
    // add noise to true flux
    OBSFLUX    = tmp_GenFlux + tmp_RanFlux_S + tmp_RanFlux_T ;
    OBSFLUXERR = tmp_GenFluxErr ; // naive obs-error is true error

    // increment sum of obsFlux and sum of error-squared.

    GENSPEC.OBSFLUX_LIST[imjd][ilam2]      += OBSFLUX ;  
    GENSPEC.OBSFLUXERRSQ_LIST[imjd][ilam2] += (OBSFLUXERR*OBSFLUXERR) ;
    GENSPEC.GENFLUX_LAMSMEAR_LIST[imjd][ilam2] += tmp_GenFlux ;

    NRAN++ ;

  } // end ilam2
    
  // - - - - - -
  if ( NRAN >= MXLAMSMEAR_SPECTROGRAPH ) {
    print_preAbort_banner(fnam);    
    printf("\t NSIGLAM  = %f \n", NSIGLAM);
    printf("\t LAMSIGMA = %f \n", LAMSIGMA );
    printf("\t LAMBIN   = %f \n", LAMBIN);
    printf("\t NBIN2    = %d \n", NBIN2 );
    sprintf(c1err,"NLAMSMEAR = %d exceeds bound of %d",
	    NRAN, MXLAMSMEAR_SPECTROGRAPH );
    sprintf(c2err,"ilam=%d LAMAVG=%.2f", ilam, LAMAVG );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  return ;

} // end GENSPEC_LAMSMEAR_LEGACY

// *************************************************
void GENSPEC_FUDGES(int imjd) {

  // Apply spectrograph fudges for tests or debugging.
  // Note that SCALE_SNR is applied in GENSPEC_DRIVER().

  int    ILAM_SPIKE     =  INPUTS.SPECTROGRAPH_OPTIONS.ILAM_SPIKE ;
  int    NBLAM          =  INPUTS_SPECTRO.NBIN_LAM ;
  int ilam, DOFUDGE=0 ;
  //  char fnam[] = "GENSPEC_FUDGES";

  // ------------- BEGIN ---------

  if ( ILAM_SPIKE     >  0    ) { DOFUDGE = 1; }
  if ( DOFUDGE == 0 ) { return ; }

  for(ilam=0; ilam < NBLAM; ilam++ ) {
    
    // check option for spectral flux to be in only one lambda bin
    if ( ILAM_SPIKE > 0 && ilam != ILAM_SPIKE ) {  
      GENSPEC.GENFLUX_LIST[imjd][ilam] =  0.0 ; 
      GENSPEC.GENMAG_LIST[imjd][ilam]  = 40.0 ;
    }

  } // end ilam

  return;

} // end GENSPEC_FUDGES


// ************************************
void  GENSPEC_FLAM(int imjd) {

  // Convert true flux in each wave bin to dF/dlam,
  // and apply optional calibration warp from WARP_SPECTRUM key(s).
  // WARP applies only to measured dF/dlam, NOT to true dF/dlam.

  int  ilam, DO_WARP ;
  double GENFLUX, FLUX, FLUXERR, FLAM, FLAMERR;
  double LAMMIN, LAMMAX, LAMBIN, LAMAVG, WARP ;
  GENPOLY_DEF *GENLAMPOLY_WARP = &INPUTS.TAKE_SPECTRUM[imjd].GENLAMPOLY_WARP ;
  int  ORDER  = GENLAMPOLY_WARP->ORDER; 
  int  NBLAM  = INPUTS_SPECTRO.NBIN_LAM ;
  char fnam[] = "GENSPEC_FLAM" ;
  int  LDMP   = 0 ; 

  // ------------ BEGIN ------------

  DO_WARP = ( ORDER >= 0 ) ;

  if ( LDMP ) {
    printf(" xxx ------------------------------------------ \n");
    printf(" xxx imjd=%d  WARP_ORDER=%d  DO_WARP=%d \n",
	   imjd, ORDER, DO_WARP );  fflush(stdout);
  }

  for(ilam=0; ilam < NBLAM; ilam++ ) {

    LAMMIN = INPUTS_SPECTRO.LAMMIN_LIST[ilam] ; 
    LAMMAX = INPUTS_SPECTRO.LAMMAX_LIST[ilam] ;    
    LAMBIN = LAMMAX - LAMMIN ;
    LAMAVG = (LAMMAX + LAMMIN)/2.0 ;

    GENFLUX  = GENSPEC.GENFLUX_LIST[imjd][ilam] ; // true flux
    FLUX     = GENSPEC.OBSFLUX_LIST[imjd][ilam] ;
    FLUXERR  = GENSPEC.OBSFLUXERR_LIST[imjd][ilam] ;

    FLAM    = FLUX / LAMBIN ;      // dF/dlam
    FLAMERR = FLUXERR / LAMBIN ;   

    WARP = 1.0;
    if ( DO_WARP ) { WARP = eval_GENPOLY(LAMAVG, GENLAMPOLY_WARP, fnam); }

    if ( LDMP && ilam < -5 ) {
      printf(" xxx \t LAM=%.1f (BIN=%.1f)  --> WARP = %8.5f  F=%le\n", 
	     LAMAVG, LAMBIN, WARP, FLUX );
      fflush(stdout);
    }

    GENSPEC.LAMMIN_LIST[imjd][ilam]   = LAMMIN ;  // Apr 2 2021
    GENSPEC.LAMMAX_LIST[imjd][ilam]   = LAMMAX ;
    GENSPEC.FLAM_LIST[imjd][ilam]     = FLAM * WARP ;
    GENSPEC.FLAMERR_LIST[imjd][ilam]  = FLAMERR ;
    GENSPEC.GENFLAM_LIST[imjd][ilam]  = GENFLUX/LAMBIN ;
    GENSPEC.FLAMWARP_LIST[imjd][ilam] = WARP ;
  }
  

} // end GENSPEC_FLAM

// ************************************
int fudge_SNR ( void ) {

  // Created Aug 22, 2009 by R.Kessler
  //
  // Fudge EXPOSURE_TIME_FILTER to force
  // SNRMAX = INPUTS.FUDGE_SNRMAX or INPUTS.FUDGE2_SNRMAX
  // in each filter. Strategy is to alternate between first
  // generation with EXPOSURE_TIME_FILT[ifilt]= SCALE_ITER1, 
  // and then set EXPOSURE_TIME based on SNRMAX and repeat event.
  //
  // Mar 9, 2011: update to handle both 
  //  1. FUDGE_SNRMAX  => adjust exposure time to change both 
  //                      ZPT and skynoise
  //  2. FUDGE2_SNRMAX => adjust sky noise only (don't change SN flux)
  //                         
  //
  // HISTORY
  // May 27, 2011: set GENLC.NOBS = 0 to avoid NOBS doubling.
  //
  // Jan 17 2014: float -> double
  //
  // May 2015: check INPUTS.IFILTOBS_FUDGE_SNRMAX and use SNRMAX_TRUE
  //           instead of measure SNRMAX

  int ifilt, ifilt_obs, ep, ITER1, ITER2, NFILT ;
  int LDMP = 0 ;

  double SNRMAX, SNRMAX_TRUE, SCALE_NOISE, SHIFT_ZPT ;
  double Tobs, DF, DFINV, DSNR1, DSNR2, DSQSNR1, DSQSNR2, DTMP ;
  double RATIO, SQRATIO ;

  char cfilt[2];
  char fnam[] = "fudge_SNR";

  // ----------- BEGIN -----------

  if ( INPUTS.OPT_FUDGE_SNRMAX == 0 )     { return 0; }

  GENLC.NOBS = 0;  // always reset NOBS since init_GENLC() may be skipped.

  ITER1 = ITER2 = 0 ;
  NFILT = GENLC.NFILTDEF_OBS;
  SCALE_NOISE = SHIFT_ZPT = SQRATIO = -9.0 ;

  if ( GENLC.FUDGE_SNRMAX_FLAG == 2 || GENLC.FUDGE_SNRMAX_FLAG == 0 ) {
    ITER1 = 1; 
    GENLC.FUDGE_SNRMAX_FLAG = 1 ;   // global flag for uniform exposure times
  }
  else if ( GENLC.FUDGE_SNRMAX_FLAG == 1 ) {
    ITER2 = 1; 
    GENLC.FUDGE_SNRMAX_FLAG = 2 ;   // flag for fudged exposure times
  }
  else {
    sprintf(c1err,"Invalid GENLC.FUDGE_SNRMAX_FLAG=%d",
	    GENLC.FUDGE_SNRMAX_FLAG ) ;
    sprintf(c2err,"Something is really messed up.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( LDMP && ITER2 ) { printf(" xxx ---------------------------- \n"); }

  // ----------------------------------------

  for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    SNRMAX      = GENLC.SNRMAX_FILT[ifilt_obs];  // measured (with fluctuations)
    SNRMAX_TRUE = GENLC.SNRMAX_TRUE[ifilt_obs];  // true flux / trueERR

    if ( INPUTS.OPT_FUDGE_SNRMAX == 1 ) {  // adjust EXPOSURE_TIME
      if ( ITER1 ) 
	{ RATIO = 1.0 ; }
      else 
	{ RATIO = INPUTS.FUDGE_SNRMAX / SNRMAX_TRUE ; }     

      SQRATIO     = RATIO * RATIO ;
      SCALE_NOISE = RATIO;
      SHIFT_ZPT   = 5.0 * log10(RATIO);  
    }
    else if ( INPUTS.OPT_FUDGE_SNRMAX == 2 ) { // adjust sigSKY only
      SQRATIO     = 1.0 ; // don't change signal with this option
      SHIFT_ZPT   = 0.0 ;
      if ( ITER1 ) 
	{ SCALE_NOISE = 1.0 ; }
      else { 
	// do SKY-NOISE-SCALE calc. in double precision
	DF      = GENLC.FLUXpe_at_SNRMAX[ifilt_obs] ;
	DFINV   = 1.0 / (DF + 1.0E-40) ;  // avoid 1/0
	DSNR1   = SNRMAX_TRUE;                  // 1st-iter SNRMAX
	DSNR2   = (double)INPUTS.FUDGE_SNRMAX;  // target SNRMAX
	DSQSNR1 = DSNR1 * DSNR1 ;
	DSQSNR2 = DSNR2 * DSNR2 ;
	DTMP    = (1./DSQSNR2 - DFINV) / (1./DSQSNR1 - DFINV) ;
	SCALE_NOISE = sqrt(DTMP) ;  
      }
    }

    INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] = SQRATIO ;
    GENLC.SCALE_NOISE[ifilt_obs]           = SCALE_NOISE ;
    GENLC.SHIFT_ZPTSIMLIB[ifilt_obs]       = SHIFT_ZPT ;
    
    if ( LDMP && ITER2 ) {
      sprintf(cfilt, "%c",  FILTERSTRING[ifilt_obs] );
      printf(" xxx ITER-%d : SNRMAX(%s)=%6.1f/%6.1f   RATIO^2=%9.3f  EXPO->%9.2f \n",
	     GENLC.FUDGE_SNRMAX_FLAG, cfilt, SNRMAX, SNRMAX_TRUE, SQRATIO,
	     INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] );
      /*
      printf("\t SCALE_NOISE=%7.3f  SHIFT_ZPT=%6.3f \n"
	     , GENLC.SCALE_NOISE[ifilt_obs]
	     , GENLC.SHIFT_ZPTSIMLIB[ifilt_obs] 
	     ) ;  */
    }

  } // end of ifilt loop


  // ----------------------------------
  // May 2015
  // If IFILTOBS_FUDGE_SNRMAX > 0, then set each filter's exposure time
  // to the exposure time of this filter.

  
  int IFILTOBS_FIX = INPUTS.IFILTOBS_FUDGE_SNRMAX ;
  if ( IFILTOBS_FIX > 0 ) {
    SQRATIO     = INPUTS.EXPOSURE_TIME_FILTER[IFILTOBS_FIX] ;
    SCALE_NOISE = GENLC.SCALE_NOISE[IFILTOBS_FIX] ;
    SHIFT_ZPT   = GENLC.SHIFT_ZPTSIMLIB[IFILTOBS_FIX] ;
    
    for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] = SQRATIO ;
      GENLC.SCALE_NOISE[ifilt_obs]           = SCALE_NOISE ;
      GENLC.SHIFT_ZPTSIMLIB[ifilt_obs]       = SHIFT_ZPT ;
    }
  }


  // --------------------------------------------
  // set and apply exposure times
  if ( ITER2 ) {
    for ( ep=1; ep <= GENLC.NEPOCH; ep++ ) {

      if ( !GENLC.OBSFLAG_GEN[ep] ) { continue; }

      ifilt_obs  = GENLC.IFILT_OBS[ep] ;
      Tobs       = GENLC.MJD[ep] - GENLC.PEAKMJD ;  // for dump only

      SCALE_NOISE = GENLC.SCALE_NOISE[ifilt_obs] ;
      SHIFT_ZPT   = GENLC.SHIFT_ZPTSIMLIB[ifilt_obs];

      // Mar 1 2014: template scaling below is NOT tested ... beware !
      SIMLIB_OBS_GEN.SKYSIG[ep]     *= SCALE_NOISE ;
      SIMLIB_OBS_GEN.READNOISE[ep]  *= SCALE_NOISE ;
      SIMLIB_OBS_GEN.ZPTADU[ep]     += SHIFT_ZPT  ;

      SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[ep]     *= SCALE_NOISE ;
      SIMLIB_OBS_GEN.TEMPLATE_READNOISE[ep]  *= SCALE_NOISE ;
      SIMLIB_OBS_GEN.TEMPLATE_ZPT[ep]        += SHIFT_ZPT  ;

      if ( LDMP && fabs(Tobs) < 1.0 ) {
	sprintf(cfilt, "%c",  FILTERSTRING[ifilt_obs] );
	printf(" Peak-SKYSIG[%s] = %9.3f  ZPT=%6.3f  SCALE_NOISE=%9.3f\n"
	       ,cfilt
	       ,SIMLIB_OBS_GEN.SKYSIG[ep]
	       ,SIMLIB_OBS_GEN.ZPTADU[ep] 
	       ,SCALE_NOISE	       );
      }      

    } // end of ep
  }  // end of ITER2

  return ( GENLC.FUDGE_SNRMAX_FLAG );


} // end of fudge_SNR


// *******************************
void gen_event_driver(int ilc) {

  /**************************************************
    generate the following:
      - CID
      - REDSHIFT_CMB
      - PEAKMJD
      - AV
      - RISETIME_SHIFT
      - MWEBV
      - GENMAG_OFF_GLOBAL

  These quantities are either read from data files
  to mimic real data, or they are generated randomly.



  Apr 14 2016: move gen_shapepar() in here [from main] so that
               all generated quantities are before SIMLIB_read().

  Jun 11 2016: 
    call override_SNPARVAL_from_SNHOST() to ensure that all 
    HOSTLIB parameters are read before overriding.
    Fixes bug where HOSTLIB x1 & c were applied from previous event.

  Mar 1 2017: allow call to gen_AV for SIMSED model.

  Jul 19 2017: set GENLC.GENMAG_OFF_GLOBAL

  Oct 18 2017: compute and store GENLC.epoch_obs_range

  Jan 6 2018: refactor so that MU is computed after SIMLIB_READ,
              and has the zHEL dependence.

  Apr 9 2019: move gen_modelPar() call down after SNHOST_DRIVER()
              so that SN->Gal transfers work for coords or redshift.

  Apr 22 2019: 
   + GENLC.GENMAG_OFF_GLOBAL += instead of =. Adds to offset
     selected in pick_NON1ASED().
  
  Jul 2019: add strong lens option

  *********************************************************/

  int    NEPMIN = (int)INPUTS.CUTWIN_NEPOCH[0];
  double z, z1, Tobs, Trest,  zHOST, Tobs_min, Tobs_max, MWEBV  ;
  int    CID, ifilt, ifilt_obs, NEP, iep, USE_HOSTCOORD ;
  char   fnam[] = "gen_event_driver" ;

  // -------- BEGIN ----------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENRANDOM  ) {

    // pick PEAKMJD and CMB-redshift BEFORE SIMLIB_READ_DRIVER
    // so that we keep only those epochs given by
    // GENRANGE_TREST ... avoids wasting memory.
    // Note that z_helio is computed after SIMLIB_READ because
    // we need RA/DEC for z_helio.

    GENLC.PEAKMJD = gen_peakmjd(); 

    // set GENLC.RATEPAR struct (Aug 2016)
    set_RATEPAR(ilc, &INPUTS.NON1ASED);

    GENLC.REDSHIFT_CMB = gen_redshift_cmb();

    // check for strong lens multiple images before reading SIMLIB
    // so that SIMLIB-MJDRANGE to read is based on all images
    if ( INPUTS_STRONGLENS.USE_FLAG ) {
      gen_event_stronglens(ilc,1); 
      if ( GENSL.REPEAT_FLAG ) { goto  LOAD_TOBS ; }
    }

    // read entry from libray after generated PEAKMJD and redshift ;
    // see comment above.

    SIMLIB_READ_DRIVER();

    GENLC.CID   = GENLC.CIDOFF + ilc ; 

    if ( GENLC.NEPOCH < NEPMIN ) { return ; }

    GENLC.REDSHIFT_HELIO = gen_redshift_helio(); // needs RA,DEC from SIMLIB

    gen_distanceMag( GENLC.REDSHIFT_CMB,            // input
		     GENLC.REDSHIFT_HELIO,          // input 
		     &GENLC.DLMU, &GENLC.LENSDMU ); // returned 

    // - - - - -   Tricky MWEBV LOGIC (Spaghetti alert) - - - - - - -
    // get MWEBV here if SNHOST_DRIVER will NOT change the SN coords; 
    // otherwise skip it here and let SNHOST_DRIVER generate MWEBV 
    // after it picks new coords.
    // Don't call gen_MWEBV after SNHOST_DRIVER because it needs
    // MWEBV for the galaxy mags.
    USE_HOSTCOORD = ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SN2GAL_RADEC ) ;
    if ( !USE_HOSTCOORD  )  { MWEBV = gen_MWEBV(GENLC.RA,GENLC.DEC); }


    // Get redshift of host: equals SN helio redshift, or that of wrongHost.
    zHOST = gen_zHOST(GENLC.CID, GENLC.REDSHIFT_HELIO, 
		      &GENLC.CORRECT_HOSTMATCH );
    GENLC.REDSHIFT_HOST = zHOST ;  // helio frame

    // Mar 14 2020: move modelPar before GEN_SNHOST so that WGTMAP
    // can depend on SN params (c,x1) without adding c,x1 in HOSTLIB
    if ( HOSTLIB_WGTMAP.N_SNVAR > 0  )
      { gen_modelPar(ilc, OPT_FRAME_REST);  }

    // Fetch host-galaxy using HOSTLIB (except for LCLIB model)
    // Note that SNHOST_DRIVER can change GENLC.REDSHIFT_CMB 
    // and DLMAG to match that of the HOST
    // Similarly, GENLC.REDSHIFT_HOST is changed to be the true zhost
    GEN_SNHOST_DRIVER(zHOST, GENLC.PEAKMJD); 

    // Jun 12 2020 
    //  if no SN par in WGTMAP, generate SN params after picking host
    //  (spagetti alert here)
    if ( HOSTLIB_WGTMAP.N_SNVAR==0  )
      { gen_modelPar(ilc, OPT_FRAME_REST);  }

    // pick model params AFTER redshift/host selection (4.09.2019),
    // and note that GEN_SNHOST can modify GENLC.DLMU that is used
    // for model params. e.g.: for SALT2, params are c,x1,alpha,beta.
    gen_modelPar(ilc, OPT_FRAME_OBS ); 

    // - - - - - - - 
    // get host galaxy extinction for rest-frame models and for NON1A
    gen_modelPar_dust(GENFRAME_OPT);

    // check for SN params in HOSTLIB
    override_modelPar_from_SNHOST(); 

    // Now get smeared [measured] redshift and PEAKMJD
    // Smeared quantities are written to SNDATA files,
    // but are not used for anything else.

    if ( INPUTS.GENSIGMA_REDSHIFT >= 0.0 )
      { gen_zsmear( INPUTS.GENSIGMA_REDSHIFT ); }  

    // global mag offset + z-dependence 
    GENLC.GENMAG_OFF_GLOBAL += (double)INPUTS.GENMAG_OFF_GLOBAL
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENMAG_OFF_GLOBAL");

    // option to randomly shift coords to avoid spatial overlap 
    gen_random_coord_shift();

    // option to magnify SN and generate multiple LCs
    gen_event_stronglens(ilc,2);

  } 

  else if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID ) {
#ifdef MODELGRID_GEN
    gen_GRIDevent(ilc);
    gen_modelPar(ilc, OPT_FRAME_REST) ;
    gen_modelPar(ilc, OPT_FRAME_OBS ) ;
#endif
    return ;
  }
  
  // --------------------------
  if ( WRFLAG_CIDRAN > 0 ) {
    CID = INPUTS.CIDRAN_LIST[GENLC.CID-INPUTS.CIDOFF];
    GENLC.CIDRAN = CID ;

    if ( CID < 0 ) {
      sprintf(c1err,"Invalid CIDRAN=%d", CID);
      sprintf(c2err,"CID=%d CIDOFF=%d", GENLC.CID, INPUTS.CIDOFF);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  // ----------------------------------------------------------
  // misc. filter-dependent stuff:
  // + tack on peakMJDs at very end of list.

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];

    // skip this filter if it is not in this SIMLIB entry.
    if ( GENLC.SIMLIB_USEFILT_ENTRY[ifilt_obs] == 0 &&
	 INPUTS.CUTWIN_NEPOCH[0] > 0 )
      { continue ; }

    // check option to add artificial template epoch
    float MJD_TEMPLATE = INPUTS.MJD_TEMPLATE_FILTER[ifilt_obs] ;
    if ( MJD_TEMPLATE > 1.0 ) {
      GENLC.NEPOCH++ ; NEP   = GENLC.NEPOCH ;
      GENLC.MJD[NEP]         = MJD_TEMPLATE ;
      GENLC.IFILT_OBS[NEP]   = ifilt_obs ;
      GENLC.OBSFLAG_TEMPLATE[NEP]  = true ;
      GENLC.OBSFLAG_GEN[NEP]       = false ;
      GENLC.IEPOCH_TEMPLATE[ifilt_obs] = NEP ; 
    }

    
    // always add artificial PEAKMJD epoch to be last
    GENLC.NEPOCH++ ; NEP = GENLC.NEPOCH ;
    GENLC.MJD[NEP]       = GENLC.PEAKMJD ;
    GENLC.IFILT_OBS[NEP] = ifilt_obs ;
    GENLC.OBSFLAG_PEAK[NEP]    = true ;
    GENLC.OBSFLAG_GEN[NEP]     = false ;
    GENLC.IEPOCH_PEAK[ifilt_obs] = NEP ; 

  }  // ifilt_obs loop

    
  // ----------------------------------------------------------

  // pick random shift in rise & fall-times
  genshift_risefalltimes();

 LOAD_TOBS:

  // Compute epochs relative to peak
  z = GENLC.REDSHIFT_HELIO ;  z1 = 1.0 + z;
  Tobs_min = 1.0E9;  Tobs_max= -1.0E-9;
  for ( iep=1; iep <= GENLC.NEPOCH ; iep++ ) {
    Tobs = GENLC.MJD[iep] - GENLC.PEAKMJD ;
    Trest = Tobs/z1 ;
    GENLC.epoch_obs[iep]  = Tobs  ;
    GENLC.epoch_rest[iep] = Trest ;
    if ( Tobs < Tobs_min )  { Tobs_min = Tobs; }
    if ( Tobs > Tobs_max )  { Tobs_max = Tobs; }
  }
  GENLC.epoch_obs_range[0] = Tobs_min ;  // global range including all bands
  GENLC.epoch_obs_range[1] = Tobs_max ;  // --> needed for LCLIB

  return ;

}  // end of gen_event_driver


// *******************************************
void override_modelPar_from_SNHOST(void) {

  // Jun 2016
  // check if any SN parameters should be overwritten by
  // value from HOSTLIB, to enable SNpar-host correlations
  //
  // Mar 23 2018: allow SNMAGSHIFT or USESNPAR
  // May 23 2019: adjust amplitude for SALT2gammaDM
  // May 28 2020: fix call for RV, and add AV & EBV options too.

  double GAMMA_GRID_MIN = INPUTS.BIASCOR_SALT2GAMMA_GRID[0];
  double GAMMA_GRID_MAX = INPUTS.BIASCOR_SALT2GAMMA_GRID[1];
  int USE1, USE2, USE3 ;
  double DM_HOSTCOR, shape, PKMJD, RV, AV, EBV ;
  //  char fnam[] = "override_modelPar_from_SNHOST" ;

  // ---------------- BEGIN ------------------

  // if USESNPAR option is not set, then return PARVAL_ORIG immediately.        
  USE1   = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) ;
  USE2   = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_SNMAGSHIFT) ;
  USE3   = (GAMMA_GRID_MAX > GAMMA_GRID_MIN );
  if ( USE1==0 && USE2==0 && USE3==0 ) { return ; }

  // - - - - - -  load SNMAGSHIFT - - - -   
  // now check if SNMAGSHIFT is one of the hostlib variables
  DM_HOSTCOR = modelPar_from_SNHOST(SNHOSTGAL.WGTMAP_SNMAGSHIFT,
				    HOSTLIB_VARNAME_SNMAGSHIFT);

  // check option to pick shape param from HOSTLIB 
  shape  = modelPar_from_SNHOST(GENLC.SHAPEPAR,GENLC.SHAPEPAR_NAME);
  *GENLC.ptr_SHAPEPAR = shape ; 
  GENLC.SHAPEPAR      = shape ;

  // check option to get peakMJD from hostlib
  PKMJD         = GENLC.PEAKMJD ;
  GENLC.PEAKMJD = modelPar_from_SNHOST(PKMJD,"PEAKMJD") ;

  RV       = GENLC.RV ;
  GENLC.RV = modelPar_from_SNHOST(RV,PARNAME_RV);

  AV       = GENLC.AV ;
  GENLC.AV = modelPar_from_SNHOST(AV,PARNAME_AV);

  // check for EBV in HOSTLOB ... update AV
  if ( GENLC.RV > 0.001 ) {
    EBV = GENLC.AV/GENLC.RV ;
    GENLC.AV  = GENLC.RV * modelPar_from_SNHOST(EBV,PARNAME_EBV);
  }

  if ( INDEX_GENMODEL  == MODEL_SALT2 ) {
    double a = GENLC.SALT2alpha ;
    double b = GENLC.SALT2beta  ;
    double c = GENLC.SALT2c ;
    // note that generic shape param (x1) is modified above.

    GENLC.SALT2c =   // optional overwrite from HOSTLIB
      modelPar_from_SNHOST(c, GENLC.COLORPAR_NAME);   
    GENLC.SALT2alpha = // optional overwrite from HOSTLIB
      modelPar_from_SNHOST(a, GENLC.SHAPEPAR2_NAME);
    GENLC.SALT2beta = // optional overwrite from HOSTLIB
      modelPar_from_SNHOST(b, GENLC.COLORPAR2_NAME);    

    // re-compute x0 and mB since x1,c have changed (Jun 15 2016)
    GENLC.SALT2x0 = SALT2x0calc(GENLC.SALT2alpha, GENLC.SALT2beta, 
				GENLC.SALT2x1, GENLC.SALT2c, 
				GENLC.DLMU );   
    GENLC.SALT2mB = SALT2mBcalc(GENLC.SALT2x0) ;

    // May 23 2019: adjust amplitude for SNMAGSHIFT_HOSTCOR
    if ( DM_HOSTCOR != 0.0 ) {  GENLC.SALT2gammaDM = DM_HOSTCOR ; }
  }

  return ;

} // end override_modelPar_from_HOSTLIB


// *******************************************
void gen_random_coord_shift(void) {

  // Created Apr 26 2021
  // Check option to randomly shift RA,DEC within a circle of MXRADIUS
  // of original coordinate. Make sure to shift SN and all host galaxy 
  // candidates. Units in degrees.
  // Initial motivation is to avoid spatial duplicates for
  // testing alert brokers in LSST.
  //
  // BEWARE: MXRADIUS is assumed to be very small (<<1 deg)
  //   and thus MWEBV and zCMB(zHEL) are NOT updated !!!
  //

  double MXRADIUS = (double)INPUTS.MXRADIUS_RANDOM_SHIFT;
  if ( MXRADIUS < 1.0E-12 ) { return ; }

  double RAD       = RADIAN ;
  double ran_r, ran_phi, r, phi, ANGSEP ;
  double shift_RA=0.0, shift_DEC=0.0;
  double RA = GENLC.RA, DEC=GENLC.DEC;
  double COSDEC = cos(DEC*RAD);
  int    m;
  bool CHECK_ANGSEP = true;
  int  LDMP = 0 ;
  char fnam[] = "gen_random_coord_shift" ;

  // -------------- BEGIN ------------

  // pick random radius and randam azimuth angle
  ran_r   = getRan_Flat1(1);   // random 0-1
  ran_phi = getRan_Flat1(1); 

  phi = TWOPI    * ran_phi ;
  r   = MXRADIUS * sqrt(ran_r);

  shift_RA    = r * cos(phi) / COSDEC ;
  shift_DEC   = r * sin(phi) ;

  // load globals for SIMGEN_DUMP file
  GENLC.random_shift_RA     = shift_RA;
  GENLC.random_shift_DEC    = shift_DEC;
  GENLC.random_shift_RADIUS = r   ;
  GENLC.random_shift_PHI    = phi ;

  // - - - - -
  // apply shifts

  // start with SN ...
  GENLC.RA  += shift_RA;  
  GENLC.DEC += shift_DEC;
 
  // now the host(s)
  for(m=0; m < SNHOSTGAL.NNBR; m++ ) {
    SNHOSTGAL_DDLR_SORT[m].RA  += shift_RA ;
    SNHOSTGAL_DDLR_SORT[m].DEC += shift_DEC ;
  }

  if ( CHECK_ANGSEP ) {
    ANGSEP = angSep(GENLC.RA,GENLC.DEC,  RA,DEC, (double)1.0 ) ;
    double ratio = ANGSEP/r;
    double dif_arcsec = 3600.0*(ANGSEP-r);
    if ( fabs(ratio-1.0) > 1.0E-3 ) {
      print_preAbort_banner(fnam);
      printf(" CID = %d \n", GENLC.CID);
      printf(" Original RA,DEC = %.3f, %.3f deg \n", RA, DEC);
      printf(" random phi    = %f radians\n", phi);
      printf(" random radius = %f degrees \n", r);
      printf(" shift(RA,DEC) = %f, %f \n", shift_RA, shift_DEC);
      printf(" actual ANGSEP = %f \n", ANGSEP);
      printf(" ANGSEP(actual)/ANGSEP(expect) = %le \n", ratio);
      printf(" ANGSEP(actual)-ANGSEP(expect) = %le arcSec\n", dif_arcsec);

      sprintf(c1err,"Problem with random coord shift.");
      sprintf(c2err,"See preAbort info above.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  return;

} // end gen_random_coord_shift

// *******************************************
void gen_event_stronglens(int ilc, int istage) {

  // Created July 2019 by R.Kessler
  // Generate multiple SL images, each with
  // time delay, magnification, angle shift.
  // 
  // Input :
  //  ilc  :  internal LC index used to set CID
  //
  //  istage=1 --> PEAKMJD and redshift are generated; 
  //               SIMLIB not read, and thus RA,DEC are not known.
  //
  //  istage=2 --> SIMLIB has been read (RA,DEC known)
  //

  int    INIT_FLAG = GENSL.INIT_FLAG;
  int    NIMG      = GENSL.NIMG;
  int    IMGNUM    = GENSL.IMGNUM;
  double TRESTMIN  = INPUTS.GENRANGE_TREST[0];
  double TRESTMAX  = INPUTS.GENRANGE_TREST[1];
  int    MEMD      = MXIMG_STRONGLENS * sizeof(double);
  double RAD       = RADIAN;
  int    LDMP      = 0; // (ilc < 4) ; 

  double zLENS, zSN=-9.0, z1, hostpar[10];
  double PEAKMJD, tdelay_min=1.0E9, tdelay_max=-1.0E9;
  double tdelay=0.0,  magnif=0.0, magshift=0.0;
  double XIMG=0.0, YIMG=0.0;
  double cosDEC, ANGSEP_TRUE ;
  int    NEXTLENS=0, IDLENS=0, blend_flag, img, NGEN_MIN, ep ;
  char fnam[] = "gen_event_stronglens";

  // ------------- BEGIN ------------------

  
  GENSL.REPEAT_FLAG  =  0 ;
  if ( !INPUTS_STRONGLENS.USE_FLAG ) { return; }

  if ( WRFLAG_CIDRAN ) {
    sprintf(c1err,"Cannot use CIDRAN option with strong lens model.");
    sprintf(c2err,"Remove %d from FORMAT_MASK", WRMASK_CIDRAN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  GENLC.CID = GENLC.CIDOFF + ilc ; 

  // ------------------
  if ( INIT_FLAG == 0 ) {
    GENSL.TDELAY_LIST   = (double*) malloc(MEMD);
    GENSL.MAGNIF_LIST   = (double*) malloc(MEMD);
    GENSL.MAGSHIFT_LIST = (double*) malloc(MEMD);
    GENSL.XIMG_LIST     = (double*) malloc(MEMD);
    GENSL.YIMG_LIST     = (double*) malloc(MEMD);
    GENSL.INIT_FLAG     = 1;
  }

  if ( INPUTS.USE_SIMLIB_REDSHIFT ) {
    sprintf(c1err,"Cannot use USE_SIMLIB_REDSHIFT option with strong lens");
    sprintf(c2err,"Check sim-input file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( INPUTS.USE_SIMLIB_PEAKMJD ) {
    sprintf(c1err,"Cannot use USE_SIMLIB_PEAKMJD option with strong lens");
    sprintf(c2err,"Check sim-input file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // -----------------------
  if ( istage == 2 ) {
    if ( NIMG == 0 ) { return; } // May 2020
    if ( !GENSL.REPEAT_FLAG ) {
      // store original coords on first image
      GENSL.RA_noSL    = GENLC.RA;
      GENSL.DEC_noSL   = GENLC.DEC ;
    }

    XIMG          = GENSL.XIMG_LIST[IMGNUM] ;   // arcSec
    YIMG          = GENSL.YIMG_LIST[IMGNUM] ;   // arcSec
    cosDEC        = cos(RAD*GENSL.DEC_noSL) ;
    GENLC.RA      = GENSL.RA_noSL  + (XIMG/3600.0)/cosDEC ;
    GENLC.DEC     = GENSL.DEC_noSL + (YIMG/3600.0) ;

    if ( fabs(GENLC.RA) > 400.0 || fabs(GENLC.DEC) > 400.0 ) {
      sprintf(c1err,"Insane RA,DEC = %f, %f", GENLC.RA, GENLC.DEC);
      sprintf(c2err,"IDLENS=%d, X,Yimg=%.2f,%.2f arcSec", 
	      GENSL.IDLENS, XIMG, YIMG );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  
    goto DONE ;
  }


  // istage=1
  // check if next lens is to be generated
  NEXTLENS = ( IMGNUM == NIMG-1 );

  if ( NEXTLENS ) {

    zSN       = GENLC.REDSHIFT_CMB;
    z1        = 1.0 + zSN;

    get_stronglens(zSN, hostpar, LDMP,  // <== inputs
		   &IDLENS, &zLENS, &blend_flag, // <== returned
		   &GENSL.NIMG,         // <== returned
		   GENSL.TDELAY_LIST, GENSL.MAGNIF_LIST, 
		   GENSL.XIMG_LIST, GENSL.YIMG_LIST );
    
    GENSL.IDLENS       = IDLENS;
    GENSL.zSN          = zSN ;
    GENSL.zLENS        = zLENS;
    GENSL.BLEND_FLAG   = blend_flag ;
    GENSL.IMGNUM       = -1;
    if ( GENSL.NIMG == 0 ) 
      { GENSL.IDLENS = -9; GENSL.zLENS = -9.0;   goto DONE ;  }

    // get min and max delay for reading enough of the cadence
    tdelay_min = +1.0E9 ;
    tdelay_max = -1.0E9 ;
    for(img=0; img < GENSL.NIMG; img++ ) {
      tdelay = GENSL.TDELAY_LIST[img];
      if ( tdelay < tdelay_min ) { tdelay_min = tdelay; }
      if ( tdelay > tdelay_max ) { tdelay_max = tdelay; }
    }

    // set MJD range to cover all time delays; used to read SIMLIB cadence
    PEAKMJD            = GENLC.PEAKMJD ; // without lens
    GENSL.PEAKMJD_noSL = PEAKMJD ;
    GENSL.MJDMIN       = PEAKMJD + z1*TRESTMIN + tdelay_min - 0.1;
    GENSL.MJDMAX       = PEAKMJD + z1*TRESTMAX + tdelay_max + 0.1;

    // check if NGEN needs to be increased to generate all of the lenses.
    // This affects only the last event.
    // Beware that actual number of generated events may exceed
    // requested NGENTOT_LC by a few.
    NGEN_MIN = ilc + GENSL.NIMG - 1 ;
    if ( NGEN_MIN > INPUTS.NGEN && GENSL.INIT_FLAG != 777 ) 
      { INPUTS.NGEN = NGEN_MIN; GENSL.INIT_FLAG=777; }

  }

  //  - - - - - - - - - - - - - - - -
  // get peakmjd for this image
  GENSL.IMGNUM++ ;  IMGNUM = GENSL.IMGNUM; 
  tdelay        = GENSL.TDELAY_LIST[IMGNUM];
  PEAKMJD       = GENSL.PEAKMJD_noSL + tdelay;
  GENLC.PEAKMJD = PEAKMJD ;

  // update MJD for epoch flagged as peakMJD
  for(ep=0; ep <= GENLC.NEPOCH; ep++ ) {
    if ( GENLC.OBSFLAG_PEAK[ep] ) { GENLC.MJD[ep] = GENLC.PEAKMJD ; }
  }


  // convert magnifation to magshift
  magnif   = GENSL.MAGNIF_LIST[IMGNUM];
  magshift = -2.5*log10(magnif);
  GENSL.MAGSHIFT_LIST[IMGNUM] = magshift ;
  GENLC.SL_MAGSHIFT = magshift; // for simgen-dump

  // restore same SN redshift 
  zSN = GENSL.zSN;
  GENLC.REDSHIFT_CMB = GENSL.zSN;

  // set REPEAT flag for other sim functions
  GENSL.REPEAT_FLAG  = ( GENSL.IMGNUM > 0 ) ;

  // sanity checks
  if ( PEAKMJD < 2.0E4 || PEAKMJD > 1.0E5 ) {
    sprintf(c1err,"Insane PEAKMJD = %f (IMGNUM=%d of %d) ", 
	    PEAKMJD, IMGNUM, GENSL.NIMG );
    sprintf(c2err,"IDLENS=%d, tdelay = %f ", GENSL.IDLENS, tdelay);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

 DONE:



  // - - - - - - - - 
  if ( LDMP ) {

    if ( istage == 1 && IMGNUM == 0 ) {
      printf(" xxx ========================================"
	     "=============================== \n");
    }

    printf(" xxx ------------ %s DUMP CID=%d istage=%d ------------- \n", 
	   fnam, GENLC.CID, istage );
    if ( istage == 1 ) {
       printf(" xxx NEXTLENS=%d  REPEAT=%d  IMGNUM=%d of %d  zSN=%.4f \n",
	      NEXTLENS, GENSL.REPEAT_FLAG, IMGNUM, GENSL.NIMG, zSN );
       printf(" xxx PEAKMJD_noSL=%.2f, PEAKMJD=%.2f  TDELAY=%.2f \n",
	      GENSL.PEAKMJD_noSL, GENLC.PEAKMJD, tdelay);
       if ( IMGNUM == 0 ) {
	 printf(" xxx PEAKMJD[min,max] = %.3f to %.3f \n",
		GENSL.MJDMIN, GENSL.MJDMAX );    
	 printf(" xxx tdelay[min,max] = %.2f to %.2f \n", 
		tdelay_min, tdelay_max);
       }
    }
    else {
      printf(" xxx RA: %.5f -> %.5f deg ", GENSL.RA_noSL, GENLC.RA);
      printf("   DEC: %.5f -> %.5f deg \n", GENSL.DEC_noSL, GENLC.DEC);

      ANGSEP_TRUE = angSep(GENLC.RA, GENLC.DEC, 
			   GENSL.RA_noSL, GENSL.DEC_noSL, (double)3600.0 );
      printf(" xxx XIMG,YIMG = %.3f, %.3f  ANGSEP=%.3f arcsec \n",
	     XIMG, YIMG, ANGSEP_TRUE );

    }

    fflush(stdout);
  }

  return ;

} // end gen_event_stronglens



// *******************************************
void gen_event_reject(int *ILC, SIMFILE_AUX_DEF *SIMFILE_AUX,
		      char *REJECT_STAGE ) {


  // Sep 11 2015
  // Take care of stuff for rejecting generated event.
  // * Decrement ilc counter if non-zero
  // * do stuff based on REJECT_STAGE
  //
  // Mar 18 2018: add separate category for NEPOCH 

  int ilc_orig, ilc;
  char fnam[] = "gen_event_reject" ;

  // ----------- BEGIN --------------

  ilc = ilc_orig = *ILC ;
  FREEHOST_GALID(SNHOSTGAL.IGAL);

  if ( strcmp(REJECT_STAGE,"GENRANGE") == 0 ) {
    ilc-- ;
    NGEN_REJECT.GENRANGE++ ;
  }
  else if ( strcmp(REJECT_STAGE,"GENMAG") == 0 ) {
    if ( INPUTS.NGEN_LC > 0 ) { ilc-- ; }
    NGEN_REJECT.GENMAG++ ;
  }
  else if ( strcmp(REJECT_STAGE,"SEARCHEFF") == 0 ) {
    if ( INPUTS.NGEN_LC > 0 ) { ilc-- ; }
    NGEN_REJECT.SEARCHEFF++ ;
    if ( doReject_SIMGEN_DUMP("SEARCHEFF") ) 
      { wr_SIMGEN_DUMP(2,SIMFILE_AUX); }
    
  }
  else if ( strcmp(REJECT_STAGE,"CUTWIN") == 0 ) {
    if ( INPUTS.NGEN_LC > 0 ) { ilc-- ; }
    NGEN_REJECT.CUTWIN++ ;
     if ( doReject_SIMGEN_DUMP("CUTWIN") ) 
      { wr_SIMGEN_DUMP(2,SIMFILE_AUX); }
  }
  else if ( strcmp(REJECT_STAGE,"NEPOCH") == 0 ) {  // Mar 17 2018
    if ( INPUTS.NGEN_LC > 0 ) { ilc-- ; }
    NGEN_REJECT.NEPOCH++ ;
     if ( doReject_SIMGEN_DUMP("NEPOCH") ) 
      { wr_SIMGEN_DUMP(2,SIMFILE_AUX); }
  }
  else {
    sprintf(c1err,"Undefined REJECT_STAGE = '%s'", REJECT_STAGE);
    sprintf(c1err,"at ilc = %d" , *ILC);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  
  int LDMP=0;
  if ( LDMP && *ILC != ilc ) {
    printf(" xxx %s: ILC=%d  ilc=%d  CID=%d  LIBID=%d\n",
	   fnam, *ILC, ilc, GENLC.CID, GENLC.SIMLIB_ID );
    printf(" xxx %s: NREJEVT[GENRANGE,GENMAG|SEARCH,CUT,NEP] = "
	   " %d %d | %d %d %d\n",
	   fnam, 
	   NGEN_REJECT.GENRANGE, NGEN_REJECT.GENMAG,
	   NGEN_REJECT.SEARCHEFF, NGEN_REJECT.CUTWIN, 
	   NGEN_REJECT.NEPOCH );
    fflush(stdout);
  }

  *ILC = ilc ;

  return;

} // end gen_event_reject


// *********************************
double gen_MWEBV(double RA, double DEC) {

  // compute the following globals:
  // GENLC.MWEBV        --> 'measured' value (SFD map) passed to analysis
  // GENLC.MWEBV_ERR    --> uncertainty passed to analysis.
  // GENLC.MWEBV_SMEAR  --> true value to simulate
  //
  // and function returns GENLC.MWEBV_SMEAR
  //
  // Note that the logic here is different than usual;
  // the SMEAR value is the true value rather than the
  // measured value.
  //
  // May 5 2014: use getRan_Flat for random EBV
  // Nov 26, 2016: for GENRANGE_MWEBV, apply error.
  // Jun 21, 2017: apply INPUTS.MWEBV_SCALE
  // Nov 24, 2017: modify MWEBV if coordinates have changed (see NEWCOORD)
  //
  // Jun 25 2018: 
  //  + do NOT require MWEBV<=5 because MWEBV_ERR was not re-scaled,
  //     resulting in MWEBV_SMEAR < 0
  //  + apply +-3 sigma clip to Gauss smear.
  //
  // Jul 21 2018:
  //  + if APPLYFLAG_MWEBV is true, compute FLUXCOR_MWEBV[ifilt]
  //
  // Feb 03 2021: pass RA,DEC as args so it can be called fron LCLIB


  int    RANFLAG, OPT, NEWCOORD ;
  double EBV, MWXT_GaussRan ;
  double SigmaClip[2] = { -3.0, 3.0 } ; // 3 sigma clip
  int    LDMP = 0 ;
  char fnam[] = "gen_MWEBV";

  // -------------- BEGIN ----------------

  // always burn random number to remain synced
  // MWXT_GaussRan = GaussRan(1);  
  MWXT_GaussRan = getRan_GaussClip(1,SigmaClip[0],SigmaClip[1]);

  if ( INPUTS.MWEBV_FLAG  == 0 ) { return(0.0) ; }  


  RANFLAG = ( INPUTS.GENRANGE_MWEBV[0] >= 0.0 ) ; // random MWEBV

  if ( RANFLAG ) {
    // special map-generation option 
    EBV             = getRan_Flat(1, INPUTS.GENRANGE_MWEBV );
    GENLC.MWEBV     = EBV ;
  }
  else {
    // check option to modify MWEBV using map based on sky location    

    OPT = INPUTS.OPT_MWEBV ; // user option for EBV map

    // if EBV is from SIMLIB, and it is <= 0 (i.e., missing), 
    // then switch to SFD98
    if ( OPT == OPT_MWEBV_FILE && GENLC.MWEBV <= 0.0 ) 
      { OPT = OPT_MWEBV_SFD98 ; }

    // modify MWEBV only if coordinates have changed
    // (since fetching MWEBV is slow)

    NEWCOORD = (RA != GENLC.RA_LAST || DEC != GENLC.DEC_LAST) ;

    if ( NEWCOORD ) 
      { modify_MWEBV_SFD(OPT, RA,DEC, &GENLC.MWEBV, &GENLC.MWEBV_ERR); }
   
  }


  // re-compute smearing-error only for FILE option
  // Note defaults are ERR1=0 and ERR2 = MWEBV/6.
  if ( INPUTS.OPT_MWEBV == OPT_MWEBV_FILE ) {
    double SQERR1, SQERR2, ERR1, ERR2 ;
    ERR1 = (double)INPUTS.MWEBV_SIG ;  // fixed error (default is 0)
    ERR2 = (double)(INPUTS.MWEBV_SIGRATIO * GENLC.MWEBV); 
    SQERR1 = ERR1 * ERR1 ;
    SQERR2 = ERR2 * ERR2 ;
    GENLC.MWEBV_ERR = sqrt(SQERR1 + SQERR2);
  }
  else if ( INPUTS.MWEBV_SIG > 0.0 ) {
    GENLC.MWEBV_ERR = INPUTS.MWEBV_SIG ; // for systematic test
  }
  else if ( INPUTS.MWEBV_SIGRATIO == 0.0 && INPUTS.MWEBV_SIG==0.0 ) {
    GENLC.MWEBV_ERR = 0.0 ;  // turn off MWEBV error (Dec 2018)
  }

  // compute smeared (true) value to generate extinction
  GENLC.MWEBV_SMEAR =
    GENLC.MWEBV +                        // value from extinc. map
    GENLC.MWEBV_ERR * MWXT_GaussRan  +   // smearing from uncert.
    INPUTS.MWEBV_SHIFT                   // systematic shift
    ;
  
  // apply user scale (June 2017)
  GENLC.MWEBV_SMEAR *= INPUTS.MWEBV_SCALE ;

  if ( LDMP ) {
    printf(" xxx --------------------------------------- \n");
    printf(" xxx CID=%2d NEWC=%d  RA=%.6f  DEC=%.6f  "
	   "--> True MWEBV=%.4f +- %.4f \n",
	   GENLC.CID, NEWCOORD, RA, DEC, GENLC.MWEBV, GENLC.MWEBV_ERR );
    printf(" xxx MWEBV_SMEAR = %.4f  (OPT_MWEBV=%d) \n", 
	   GENLC.MWEBV_SMEAR, INPUTS.OPT_MWEBV );
    printf(" xxx RA_LAST=%f  DEC_LAST=%f \n", GENLC.RA_LAST, GENLC.DEC_LAST);
  }

  // - - - - - - - - - - - - - - - - - 
  // Jul 21 2018: checking option to compute MWEBV-FLUXCOR for each band.
  //   Used to correct output fluxes for MWEBV (e..g, for public challenge)

  if ( INPUTS.APPLYFLAG_MWEBV ) {

    int ifilt, ifilt_obs ;
    double MCOR_MAP, MCOR_TRUE, FCOR_MAP, LAMOBS;
    double RV       = INPUTS.RV_MWCOLORLAW ;
    double AV_MAP   = RV * GENLC.MWEBV;        // reported from map
    double AV_TRUE  = RV * GENLC.MWEBV_SMEAR ; // yes, smear=true
    int    OPT      = INPUTS.OPT_MWCOLORLAW ;
    
    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs  = GENLC.IFILTMAP_OBS[ifilt];
      LAMOBS     = (double)INPUTS.LAMAVG_OBS[ifilt_obs];
      MCOR_MAP   = GALextinct( RV, AV_MAP,  LAMOBS, OPT );
      MCOR_TRUE  = GALextinct( RV, AV_TRUE, LAMOBS, OPT );
      FCOR_MAP   = pow(TEN,0.4*MCOR_MAP); // note that FCOR > 1

      GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs]   = MCOR_TRUE;
      GENLC.MAGCOR_MWEBV_MAP[ifilt_obs]    = MCOR_MAP;
      GENLC.FLUXCOR_MWEBV_MAP[ifilt_obs]   = FCOR_MAP ;

      if ( LDMP ) {
	printf(" xxx ifilt_obs=%d --> MAGCOR[MAP,TRUE]=%5.3f,%5.3f \n",
	       ifilt_obs, MCOR_MAP, MCOR_TRUE ); fflush(stdout);
      }
    }
  }

  return(GENLC.MWEBV_SMEAR) ;

} // end of gen_MWEBV

// ****************************************
void gen_filtmap(int ilc) {

  // Created May 23, 2009 by R.Kessler
  // Moved from gen_event() so that gen_non1a() can get
  // called in between.
  // 
  // Sep 8, 2009: set ifilt_rest3
  //
  // Feb 13,2012: test lamz inside GENLC.RESTLAM_MODEL and
  //              remove call to  IFILTSTAT_SEDMODEL
  //
  // Feb 01, 2014: check ALLSKIP flag and increment NGEN_ALLSKIP
  // Sep 26 2017: give WARNING only if SIMLIB.NOBS > 0

  int  colopt, irank, istat, ALLSKIP, NDOFILT_ZERO=0;
  int  ifilt, ifilt_obs, ifilt_rest1, ifilt_rest2, ifilt_rest3, ep  ;

  float  z4, lamdif4[10];
  double z, ztmp, lamz ;
  char fnam[] = "gen_filtmap" ;

  // ------------- BEGIN ----------

  z   = GENLC.REDSHIFT_HELIO ;
  z4  = (float)z;
  colopt  = INPUTS.KCORFLAG_COLOR;
  colopt += 8; // flag to NOT ABORT if nearest-filt is not found.
  
  ALLSKIP = 1 ;  

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];

    // skip this filter if it is not in this SIMLIB entry.
    if ( GENLC.SIMLIB_USEFILT_ENTRY[ifilt_obs] == 0  ) {
      GENLC.DOFILT[ifilt_obs] = 0; 
      NDOFILT_ZERO++ ;
      continue ;
    }

    // skip this filter if it's outside the valid rest-frame
    // wavelength range of the model (Feb 2012)
    lamz = INPUTS.LAMAVG_OBS[ifilt_obs]/(1.0+z);
    if ( lamz < GENLC.RESTLAM_MODEL[0] || lamz > GENLC.RESTLAM_MODEL[1] ) {
      GENLC.DOFILT[ifilt_obs] = 0; 
      NDOFILT_ZERO++ ; 
      NSKIP_FILTER[ifilt_obs]++ ; 
      continue ;
    }


    // check if this filter is valid for SALT2/SIMSED model in observer-frame.
    if ( INDEX_GENMODEL == MODEL_SALT2  || 
	 INDEX_GENMODEL == MODEL_SIMSED	||
	 INDEX_GENMODEL == MODEL_NON1ASED   // 9.24.2017
	 ) {
      istat = IFILTSTAT_SEDMODEL(ifilt_obs, z) ; 
      if ( istat == 0 ) {
	GENLC.DOFILT[ifilt_obs] = 0 ;
	NDOFILT_ZERO++ ;
	NSKIP_FILTER[ifilt_obs]++ ; 
	continue ;
      }

      //  printf(" xxxx ifilt_obs=%d  z=%f  istat=%d \n",ifilt_obs,z8,istat);
    }

    if ( GENFRAME_OPT  == GENFRAME_REST  ) {

      irank = 1;
      ifilt_rest1 =
	nearest_ifilt_rest__(&colopt, &ifilt_obs, &irank, &z4,&lamdif4[irank]);

      irank = 2;
      ifilt_rest2 = 
	nearest_ifilt_rest__(&colopt, &ifilt_obs, &irank, &z4,&lamdif4[irank]);

      irank = 3;
      ifilt_rest3 = 
	nearest_ifilt_rest__(&colopt, &ifilt_obs, &irank, &z4,&lamdif4[irank]);
      
      GENLC.IFILTMAP_REST1[ifilt_obs] = ifilt_rest1 ;
      GENLC.IFILTMAP_REST2[ifilt_obs] = ifilt_rest2 ;
      GENLC.IFILTMAP_REST3[ifilt_obs] = ifilt_rest3 ;

      GENLC.LAMDIF_REST1[ifilt_obs] = lamdif4[1];
      GENLC.LAMDIF_REST2[ifilt_obs] = lamdif4[2];
      GENLC.LAMDIF_REST3[ifilt_obs] = lamdif4[3];
      
      // skip obs-filter if rest-frame filters are not defined.
      if ( ifilt_rest1 < 0 || ifilt_rest2 < 0 ) {
	GENLC.DOFILT[ifilt_obs] = 0; 	
	NDOFILT_ZERO++ ;
	NSKIP_FILTER[ifilt_obs]++ ;	
	continue ; 
      } 
    }


    // if DOFILT logical is still set, then update the
    // valid redshift range for this obs-filter.

    if ( GENLC.DOFILT[ifilt_obs]  ) {
      ztmp = ZVALID_FILTER[0][ifilt_obs] ;
      if ( z  < ztmp ) { ZVALID_FILTER[0][ifilt_obs] = z; }

      ztmp = ZVALID_FILTER[1][ifilt_obs] ;
      if ( z  > ztmp ) { ZVALID_FILTER[1][ifilt_obs] = z; }

      ALLSKIP = 0;  // at least one simulated epoch
    }


  }  // ifilt/ifilt_obs loop


  if ( ALLSKIP && SIMLIB_OBS_GEN.NOBS>0  ) {

    if ( NGEN_ALLSKIP < 100  ) {
      printf("  %s: WARNING %3d: no bands for "
	     "CID=%d, LIBID=%d, zCMB=%.4f  SIMLIB.NOBS=%d\n", 
	     fnam, NGEN_ALLSKIP, 
	     GENLC.CID, GENLC.SIMLIB_ID, GENLC.REDSHIFT_CMB,
	     SIMLIB_OBS_GEN.NOBS );
      fflush(stdout);
    }

    NGEN_ALLSKIP++ ; // checked in simEnd().
  }


  // Dec 22 2019: if any filter was excluded, set ISOBS[ep]=0
  if ( NDOFILT_ZERO > 0 ) {
    for ( ep = 1; ep <= GENLC.NEPOCH; ep++ )  {  
      ifilt_obs = GENLC.IFILT_OBS[ep] ;    
      if ( !GENLC.DOFILT[ifilt_obs] )  { GENLC.OBSFLAG_GEN[ep] = false ; }
    }
  }


  return ;

} // end of gen_filtmap

// ****************************************
void genshift_risefalltimes(void) {

  // Jun 12 2020: 
  //  check USE flags to set shifts, and make sure to burn randoms
  

  double shift ;
  char fnam[] = "genshift_risefalltimes";

  // ---------- BEGIN ------------

  shift = getRan_GENGAUSS_ASYM(&INPUTS.GENGAUSS_RISETIME_SHIFT);
  if ( INPUTS.GENGAUSS_RISETIME_SHIFT.USE ) 
    { GENLC.RISETIME_SHIFT = (float)shift ; }
  
  shift = getRan_GENGAUSS_ASYM(&INPUTS.GENGAUSS_FALLTIME_SHIFT);
  if ( INPUTS.GENGAUSS_FALLTIME_SHIFT.USE ) 
    { GENLC.FALLTIME_SHIFT = (float)shift ; }

} // end of genshift_risefalltimes



// *****************************************
void gen_modelPar(int ilc, int OPT_FRAME ) {

  /*********
   generate shape/luminosity & color parameters for model:
   stretch, delta, dm15 ... Also computes SALT2x0 and mB.

  Feb 27 2018: refactor with functions gen_shapepar_SALT2[SIMSED]

  Apr 09 2019: change function name, gen_shapepar() -> gen_modelPar().

  Nov 25 2019: for SIMLIB model, set SALT2x1[c]

  Mar 11 2020: pass OPT_FRAME = rest or obs.

  Jul 23 2020: DOSHAPE = F for LCLIB model.
  Aug 31 2020: DOSHAPE = F for BYOSED model
  Oct 20 2020: skip get_random_genPDF for SIMLIB model
  Dec 23 2020: skip get_random_genPDF if x1 is from HOSTLIB (ISx1_HOSTLIB)
  **********/

  bool ISFRAME_REST      = ( OPT_FRAME == OPT_FRAME_REST );
  bool ISFRAME_OBS       = ( OPT_FRAME == OPT_FRAME_OBS );
  bool ISMODEL_SALT2     = ( INDEX_GENMODEL == MODEL_SALT2  );
  bool ISMODEL_SIMSED    = ( INDEX_GENMODEL == MODEL_SIMSED );
  bool ISMODEL_FIXMAG    = ( INDEX_GENMODEL == MODEL_FIXMAG );
  bool ISMODEL_SIMLIB    = ( INDEX_GENMODEL == MODEL_SIMLIB );
  bool ISMODEL_NON1ASED  = ( INDEX_GENMODEL == MODEL_NON1ASED );
  bool ISMODEL_NON1A     = ( INPUTS.NON1A_MODELFLAG > 0 );
  bool ISMODEL_LCLIB     = ( INDEX_GENMODEL == MODEL_LCLIB ) ;

  // Check if x1/shape is extracted elsewhere (e.g., SIMLIB header or HOSTLIB).
  // If x1 is from HOSTLIB, SKIP only if SALT2x1 asymGauss is NOT defined
  // in order to preserve random sync.

  /* xxx mark delete 2021 xxxxxxx
  bool SHAPE_ASYMGAUSS   = (INPUTS.GENGAUSS_SHAPEPAR.USE);
  bool SHAPE_HOSTLIB = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) ;
  bool SKIP_SHAPE     = SHAPE_SIMLIB || (SHAPE_HOSTLIB && !SHAPE_ASYMGAUSS );
  bool DOSHAPE = !( SKIP_SHAPE || ISMODEL_SIMSED || ISMODEL_NON1A || 
		    ISMODEL_LCLIB || IS_PySEDMODEL || ISMODEL_SIMLIB );
  xxxxxxx end mark xxxxx*/

  bool SHAPE_SIMLIB = (SIMLIB_HEADER.GENGAUSS_SALT2x1.USE) ;
  bool DOSHAPE      = INPUTS.DOGEN_SHAPE && !SHAPE_SIMLIB ; // Aug 11 2021

  double ZCMB = GENLC.REDSHIFT_CMB ; // for z-dependent populations
  double shape;
  GENGAUSS_ASYM_DEF GENGAUSS_ZVAR ;
  char fnam[] = "gen_modelPar";

  //------------ BEGIN function ------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { 
    bool DOFUN = ( ISMODEL_NON1ASED || ISMODEL_SIMSED);
    if ( !DOFUN) { return ; }
  }

  // ---------------------------------------
  // evaluate shape with z-dependence on population, 

  if ( DOSHAPE && ISFRAME_REST ) {

    char *snam = GENLC.SHAPEPAR_GENNAME ;
    GENGAUSS_ZVAR = 
      get_zvariation_GENGAUSS(ZCMB, snam, &INPUTS.GENGAUSS_SHAPEPAR);

    // pick random shape value from populatoin at this redshift
    shape = get_random_genPDF(snam, &GENGAUSS_ZVAR);

    // load shape value into global GENLC struct
    GENLC.SHAPEPAR      = shape ; 
    *GENLC.ptr_SHAPEPAR = shape ; // load model-spefic shape variables 
  }

  // - - - - - - - - - - - - - - - - 

  if ( ISMODEL_SALT2 ) {
    gen_modelPar_SALT2(OPT_FRAME);
  }
  else if ( ISMODEL_NON1ASED  )  {
    if ( ISFRAME_REST ) 
      { pick_NON1ASED(ilc, &INPUTS.NON1ASED, &GENLC.NON1ASED); }
  }
  else if ( ISMODEL_SIMSED ) {
    // generate all of the SIMSED params, whatever they are.

    gen_modelPar_SIMSED(OPT_FRAME);

  } // end of SIMSED if-block

  else  if ( ISMODEL_FIXMAG ) {	  
    if ( ISFRAME_REST ) 
      { GENLC.NOSHAPE   = getRan_Flat ( 2, INPUTS.FIXMAG );  }

    if ( ISFRAME_OBS ) {
      if ( INPUTS.GENFRAME_FIXMAG == GENFRAME_REST ) 
	{ GENLC.NOSHAPE += GENLC.DLMU; }
    }
  }
  else if ( ISMODEL_SIMLIB ) {

    // this is a fragile hack to link SALT2 params with SIMLIB model.
    // Purpose is to allow cutting on SALT2x1[c] in SIMLIB header.
    if ( ISFRAME_REST ) {
      shape               = SIMLIB_HEADER.GENGAUSS_SALT2x1.RANGE[0] ;
      GENLC.SHAPEPAR      = shape ; 
      *GENLC.ptr_SHAPEPAR = shape ; 
    }
  }
  
  return ;

}  // end of gen_modelPar


//***************************************
void  gen_modelPar_SALT2(int OPT_FRAME) {

  // Created Feb 26 2018
  // Generated c, x1, alpha, beta
  // Mar 11 2020: pass OPT_FRAME argument.

  bool ISFRAME_REST  = ( OPT_FRAME == OPT_FRAME_REST );
  bool ISFRAME_OBS   = ( OPT_FRAME == OPT_FRAME_OBS  );

  /* xxx mark delete Aug 11 2021 
  bool GETc_ASYMGAUSS   = (INPUTS.GENGAUSS_SALT2c.USE);
  bool GETc_HOSTLIB     = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) ;
  bool SKIPc            = GETc_SIMLIB || (GETc_HOSTLIB && !GETc_ASYMGAUSS );
  xxxxxxxxxx mark delete xxxxxxx */
  bool GETc_SIMLIB  = (SIMLIB_HEADER.GENGAUSS_SALT2c.USE) ; // each event
  bool SKIPc        = (GETc_SIMLIB || !INPUTS.DOGEN_COLOR); // Aug 11 2021

  double   ZCMB = GENLC.REDSHIFT_CMB ; // for z-dependent populations
  GENGAUSS_ASYM_DEF  GENGAUSS_ZVAR ;
  char fnam[] = "gen_modelPar_SALT2";

  // ---------- BEGIN -----------

  // for SALT2, the color term is analogous to shapepar
  // so generate the 'c' and beta term here.

  if ( ISFRAME_REST ) {

    if ( !SKIPc ) {

      GENGAUSS_ZVAR = 
	get_zvariation_GENGAUSS(ZCMB,"SALT2c",&INPUTS.GENGAUSS_SALT2c);

      GENLC.SALT2c = 
	get_random_genPDF("SALT2c", &GENGAUSS_ZVAR );
    }

    GENGAUSS_ZVAR = 
      get_zvariation_GENGAUSS(ZCMB,"SALT2ALPHA",&INPUTS.GENGAUSS_SALT2ALPHA);
    GENLC.SALT2alpha = 
      getRan_GENGAUSS_ASYM(&GENGAUSS_ZVAR) ;   

    GENGAUSS_ZVAR = 
      get_zvariation_GENGAUSS(ZCMB,"SALT2BETA",&INPUTS.GENGAUSS_SALT2BETA);
    GENLC.SALT2beta = 
      getRan_GENGAUSS_ASYM(&GENGAUSS_ZVAR) ;   

    // 2/29/2016: optional  beta(c) polynomial 
    // 3/23/2020: refactor using GENPOLY tools
    if( INPUTS.SALT2BETA_cPOLY.ORDER >= 0 ) {
      double c = GENLC.SALT2c ;
      GENLC.SALT2beta += eval_GENPOLY(c, &INPUTS.SALT2BETA_cPOLY, fnam);
    }

    if ( GENLC.SALT2alpha < 0.0 || GENLC.SALT2beta < 0.0 ) {
      sprintf(c1err,"Invalid alpha, beta = %f, %f", 
	      GENLC.SALT2alpha, GENLC.SALT2beta  );
      sprintf(c2err,"Check sim-inputs for alpha and beta");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  } // end ISFRAME_REST

  if ( ISFRAME_OBS ) {
    // now compute x0 parameter from MU and the alpha,beta params.
    GENLC.SALT2x0 = SALT2x0calc(GENLC.SALT2alpha, GENLC.SALT2beta, 
				GENLC.SALT2x1, GENLC.SALT2c, GENLC.DLMU );
    
    GENLC.SALT2mB = SALT2mBcalc(GENLC.SALT2x0);

    /*
    printf(" xxx %s: c=%.3f  x1=%.2f  mB = %f  (a,b=%.2f,%.2f) \n", 
	   fnam, GENLC.SALT2c, GENLC.SALT2x1, GENLC.SALT2mB,
	   GENLC.SALT2alpha, GENLC.SALT2beta );
    */
  }

  return;

} // end gen_modelPar_SALT2

//***************************************
void  gen_modelPar_SIMSED(int OPT_FRAME) {

  // Created Feb 26 2018
  // Move code from gen_modelPar() to here, and add code to
  // pick correlated Guass randoms among arbitrary number of
  // variables.
  //
  // Dec 20 2018: set ISIMSED_SEQUENTIAL instead of PARVAL[0]
  //     (part of refactor for SIMSED loops)
  //
  // Apr 28 2019: return for GRIDGEN, after computing DLMU
  // Mar 11 2020: pass OPT_FRAME

  bool ISFRAME_REST    = ( OPT_FRAME == OPT_FRAME_REST );
  bool ISFRAME_OBS     = ( OPT_FRAME == OPT_FRAME_OBS  );
  int     NPAR      = INPUTS.NPAR_SIMSED;
  int     NROW_COV  = INPUTS.NROW_SIMSED_COV;
  double  ZCMB      = GENLC.REDSHIFT_CMB ; // for z-dependent populations

  int  ipar, ipar_model, genflag, opt_interp, opt_gridonly ;
  int  irow, irow_COV, NRANGEN_ITER=0 ;
  double ARG, parVal, parVal_old ;
  double PEAK, PMIN, PMAX ;
  double GAURAN[MXPAR_SIMSED], CORRVAL[MXPAR_SIMSED] ;
  GENGAUSS_ASYM_DEF GENGAUSS_ZVAR ;
  char *parName;
  char fnam[] = "gen_modelPar_SIMSED" ;
  int LDMP_COV = 0 ;

  // ----------- BEGIN ------------

  if ( ISFRAME_OBS ) {
    // use SALT2x0 parameter for SIMSED ... it just converts
    // MU into a flux-scale.
    ARG = -0.4 * GENLC.DLMU ;
    GENLC.SALT2x0 = pow(TEN , ARG );
  }

  if ( !ISFRAME_REST ) { return; }
  // everything below is rest frame.

#ifdef MODELGRID_GEN
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return; }
#endif

  parVal = parVal_old = 0.0 ;

  if ( NPAR <= 0 &&  INPUTS.OPTMASK_SIMSED == 0 ) {
    sprintf(c1err,"NPAR_SIMSED = %d", NPAR );
    errmsg(SEV_FATAL, 0, fnam, c1err, "  " ); 
  }

  // ----------------------------------------------
  // check for correlated parameters
  if ( NROW_COV > 0 ) {
  PICK_RANCOV: 
    NRANGEN_ITER++ ;
    if ( NRANGEN_ITER > 100 ) {
      sprintf(c1err,"Could not generate bounded correlated randoms");
      sprintf(c2err,"after %d tries", NRANGEN_ITER);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    
    for(irow=0; irow < NROW_COV; irow++ ) { 
      GAURAN[irow]  = getRan_Gauss(1);   // Gauss Random
      CORRVAL[irow] = 0.0 ;
    }

    
    getRan_GaussCorr(&INPUTS.SIMSED_DECOMP,GAURAN,CORRVAL); // return CORRVAL

    for(irow=0; irow < NROW_COV; irow++ ) {
      ipar = INPUTS.IPARLIST_SIMSED_COV[irow];
      PEAK = INPUTS.GENGAUSS_SIMSED[ipar].PEAK;
      PMIN = INPUTS.GENGAUSS_SIMSED[ipar].RANGE[0];
      PMAX = INPUTS.GENGAUSS_SIMSED[ipar].RANGE[1];
      CORRVAL[irow] += PEAK ;
      if ( CORRVAL[irow] > PMAX ) { goto PICK_RANCOV; }
      if ( CORRVAL[irow] < PMIN ) { goto PICK_RANCOV; }
    } // end irow loop

  } // end check on correlated randoms
 


  if ( LDMP_COV ) {
    printf(" xxx -------------------------------- \n");
    printf(" xxx NRANGEN_ITER = %d \n", NRANGEN_ITER );
    for(irow=0; irow < NROW_COV; irow++ ) {
      ipar = INPUTS.IPARLIST_SIMSED_COV[irow];
      parName      = INPUTS.PARNAME_SIMSED[ipar] ;
      printf(" xxx GAURAN= %6.3f  CORRVAL=%6.3f  (%s)\n",
	     GAURAN[irow], CORRVAL[irow], parName ); fflush(stdout);
    }
  }  


  // ------------------------------------------------
  for ( ipar = 0; ipar < NPAR ; ipar++ ) {
    ipar_model   = GENLC.SIMSED_IPARMAP[ipar]  ;
    genflag      = INPUTS.GENFLAG_SIMSED[ipar] ;
    opt_interp   = ( genflag & OPTMASK_GEN_SIMSED_PARAM    ) ;
    opt_gridonly = ( genflag & OPTMASK_GEN_SIMSED_GRIDONLY ) ;
    parName      = INPUTS.PARNAME_SIMSED[ipar] ;
    irow_COV     = INPUTS.IROWLIST_SIMSED_COV[ipar] ;
    GENLC.SIMSED_PARVAL[ipar]  = -9.0 ;

    // skip baggage params
    if ( (genflag & OPTMASK_GEN_SIMSED_param)  != 0 ) { continue ; }

    if ( opt_interp ) {

      if (  opt_gridonly ) {
	parVal = pick_gridval_SIMSED(ipar_model); 
      }
      else if ( irow_COV >= 0 ) {
	// correlated random, Gaussian only
	parVal = CORRVAL[irow_COV] ; 
      }
      else {
	// uncorrelated random, more general profile
	GENGAUSS_ZVAR = 
	  get_zvariation_GENGAUSS(ZCMB,parName,
				  &INPUTS.GENGAUSS_SIMSED[ipar]);	
	parVal = getRan_GENGAUSS_ASYM( &GENGAUSS_ZVAR );
      }

    } // opt_interp
    
    //  printf(" xxx %s: PARVAL[%d] = %f \n", fnam, ipar, parVal);
    GENLC.SIMSED_PARVAL[ipar]  = parVal ;

  }  // end ipar

  // ------------------
  // set PARVAL[0] to the SED index ... the IGEN index.
  if ( INPUTS.OPTMASK_SIMSED == OPTMASK_GEN_SIMSED_GRIDONLY ) 
  { ISIMSED_SEQUENTIAL  = (double)NGENLC_TOT ; } // IGEN = ISED

  
  return;

} // end gen_modelPar_SIMSED

// ***********************************
double pick_gridval_SIMSED(int ipar_model) {

  // Created Oct 23 2020
  // pick SIMSED value from grid using 1 of 2 options:
  //  1) default option is pick random SED within GENRANGE_XXX
  //  2) select from user-specified indices

  bool PICK_SUBSET  = (INPUTS.NINDEX_SUBSET_SIMSED_GRIDONLY > 0) ;
  int  nsed=0, ised_ran=-9, itmp=-9 ;
  double flatRan  = getRan_Flat1(2); // random between 0-1
  double PARVAL, PARVAL_TMP ;
  int    LDMP = 0 ;
  char fnam[]  = "pick_gridval_SIMSED" ;

  // ---------- BEGIN ----------

  if ( PICK_SUBSET )  {
    // pick random index from user-defined subset;
    // SIMSED_GRIDONLY:  INDEX_NAME(a,b,c,...)
    nsed = INPUTS.NINDEX_SUBSET_SIMSED_GRIDONLY ; 
    if ( ipar_model != SEDMODEL.IPAR_NON1A_INDEX ) {
      char *parName = SEDMODEL.PARNAMES[ipar_model];
      sprintf(c1err,"Cannot select GRIDONLY subset for %s", parName) ;
      sprintf(c2err,"Check SIMSED_XXX keys in sim-input file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    itmp    = (int) ( flatRan * (double)nsed ) ;
    PARVAL  = (double)INPUTS.INDEX_SUBSET_SIMSED_GRIDONLY[itmp] ;
  }
  else { 
    // pick random sed. Use GENGAUSS struct to enable 
    // selecting parameter range.
    PARVAL_TMP = getRan_GENGAUSS_ASYM(&INPUTS.GENGAUSS_SIMSED[ipar_model]);
    PARVAL     = nearest_gridval_SIMSED(ipar_model,PARVAL_TMP);
  }

  // - - - - -
  if ( LDMP ) {
    printf(" xxx %s DUMP ---------------------- \n", fnam );
    printf(" xxx %s: ipar_model = %d\n", fnam, ipar_model);
    printf(" xxx %s: PICK_SUBSET=%d  nsed=%d itmp=%d ised_ran=%d\n",	   
	   fnam, PICK_SUBSET, nsed, itmp, ised_ran);
    printf(" xxx %s: PARVAL = %f \n", fnam, PARVAL );
    fflush(stdout);
  }

  return(PARVAL);

} // end pick_gridval_SIMSED

// ***********************************
void gen_modelPar_dust(int OPT_FRAME) {

  // Created Jun 12 2020 [code moved from gen_event_driver]
  // Generate AV and RV for host galaxy.
  // 

  bool ISREST    = ( OPT_FRAME == GENFRAME_REST );
  char fnam[] = "gen_modelPar_dust" ;

  // ----------- BEGIN ------------

  if ( !INPUTS.DOGEN_AV ) { return; }

  GENLC.RV = gen_RV() ; 
  GENLC.AV = gen_AV() ;  //DJB March 20 2020:  EBV.      

  //  printf(" xxx %s: AV=%f  RV=%f \n", fnam, GENLC.AV, GENLC.RV);

  return;

} // end gen_modelPar_dust

// ***********************************
void pick_NON1ASED(int ilc,
		   INPUTS_NON1ASED_DEF *INP_NON1ASED,
		   GENLC_NON1ASED_DEF *GEN_NON1ASED) {

  /*********************
   March 2016,  R.Kessler
   Analog of gen_modelPar(), here we pick a NON1ASED template index
   based on input lightcurve GEN index  "ilc". 
 
  Aug 26 2017:  use already determined SED_FILE[ispgen] rather
                than re-constructing sedFile string.

  Jan 2019: set GENLC.SIMTYPE= user GENTYPE if specified in input file.

  ******/

  int LDMP = 0 ;
  int isp, MXGEN, ispgen, NINDEX, ifilt, ifilt_obs, NEW_ISPGEN  ;
  float MAGOFF, MAGOFF_USER, MAGSMEAR[2] ;
  char *name;
  char fnam[] = "pick_NON1ASED";

  // --------  BEGIN ----------

  NINDEX = INP_NON1ASED->NINDEX; 
  if ( NINDEX <= 0 ) { return ; }
  
  GENLC.TEMPLATE_INDEX      = -9 ;  //  (non1a class/index)
  NEW_ISPGEN = 0;

  // find sparse non1a index for this 'ilc'
  ispgen = -9 ;
  for ( isp = NINDEX; isp >= 1; isp-- ) {
    MXGEN = INP_NON1ASED->MXGEN[isp];
    if ( ilc <= MXGEN ) { ispgen = isp; }
  }

#ifdef MODELGRID_GEN
  // special case for GRID generation; 
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
    ispgen = (int)GENLC.SHAPEPAR ;  // NONIA sparse index 
  }
#endif

  if ( ispgen < 1 || ispgen > NINDEX ) {
    sprintf(c1err,"Invalid ispgen=%d (NINDEX=%d)", ispgen, NINDEX );
    sprintf(c2err,"Check NON1A keys." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // check if this is a new ISPGEN
  if ( ispgen != GEN_NON1ASED->ISPARSE )  { NEW_ISPGEN = 1 ; }

  GEN_NON1ASED->ISPARSE = ispgen ;
  GEN_NON1ASED->NGENTOT[ispgen]++ ;
  GENLC.SIMTYPE        = INP_NON1ASED->SNTAG[ispgen]; 
  if ( INPUTS.GENTYPE_SPEC > 0 ) { GENLC.SIMTYPE = INPUTS.GENTYPE_SPEC; }
  GENLC.TEMPLATE_INDEX = INP_NON1ASED->INDEX[ispgen]; 
  name   = INP_NON1ASED->LIST_NAME[GENLC.TEMPLATE_INDEX] ;

  // ----------------------------
  if ( NEW_ISPGEN == 0 ) {
    // same NON1A template extend CID range.       
    GEN_NON1ASED->CIDRANGE[ispgen][1] = GENLC.CID;
  }
  else {

    sprintf(BANNER,"Prepare '%s' Simulation for %s (INDEX=%3.3d)", 
	    MODELNAME_NON1ASED, name, GENLC.TEMPLATE_INDEX );
    print_banner(BANNER);
    
    // new template --> start next CID range
    GEN_NON1ASED->CIDRANGE[ispgen][0] = GENLC.CID;
    GEN_NON1ASED->CIDRANGE[ispgen][1] = GENLC.CID;

    init_genmag_NON1ASED( ispgen, INP_NON1ASED);

    printf("   Changes for %s : \n", name);
    
    MAGOFF      = INP_NON1ASED->MAGOFF[ispgen] ;
    MAGSMEAR[0] = INP_NON1ASED->MAGSMEAR[ispgen][0] ; 
    MAGSMEAR[1] = INP_NON1ASED->MAGSMEAR[ispgen][1] ; 
      
    if ( MAGOFF != NULLFLOAT ) {
      printf("\t GENMAG_OFF_MODEL -> %5.2f \n", MAGOFF );
      for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
	ifilt_obs  = GENLC.IFILTMAP_OBS[ifilt] ;
	MAGOFF_USER = INPUTS.TMPOFF_MODEL[ifilt] ; // Feb 2016
	INPUTS.GENMAG_OFF_MODEL[ifilt_obs] = MAGOFF_USER ;
      }
    }
    
    if ( MAGSMEAR[0] != NULLFLOAT ) {
      printf("\t GENMAG_SMEAR -> %5.2f,%5.2f \n", MAGSMEAR[0],MAGSMEAR[1] );
      INPUTS.GENMAG_SMEAR[0] = MAGSMEAR[0] ;
      INPUTS.GENMAG_SMEAR[1] = MAGSMEAR[1] ;
    }


#ifdef MODELGRID_GEN
    // Aug 20 2013: turn off mag smearing for GRID, but leave 
    // INPUTS.NON1A_MAGSMEAR as is so that its value is stored 
    // in the GRID.
    if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID )  { 
      INPUTS.GENMAG_SMEAR[0] = MAGSMEAR[0] ;  
      INPUTS.GENMAG_SMEAR[1] = MAGSMEAR[1] ;
    }
#endif

  }

  // --------------------------------------
  // Nov 17 2016: define global mag off
  // Jul 19 2017: always set global mag off
  GENLC.GENMAG_OFF_GLOBAL = 
    INPUTS.GENMAG_OFF_NON1A + 
    INP_NON1ASED->MAGOFF[ispgen] ;
  
  if ( LDMP ) {
    printf(" xxx %s: CID=%d  ilc=%d  isp=%d  TMPL=%d  MAG_OFF_GLOBAL=%.2f \n", 
	   fnam, GENLC.CID, ilc, ispgen, GENLC.TEMPLATE_INDEX,
	   GENLC.GENMAG_OFF_GLOBAL );
  }

  return ;

} // end of pick_NON1ASED


// ***********************************
void set_RATEPAR(int ilc, INPUTS_NON1ASED_DEF *INP_NON1ASED ) {

  /*********************
   Created Aug 13 2016
   Set global GENLC.RATEPAR pointer to either the default,
   or to RATEPAR_PEC1A for the PEC1A subset of NON1A(CC).

   Needed to use the correct rate-vs-redshift model for
   NON1A(CC) and PecIa in the same sim-job. 

  ******/

  int  isp, MXGEN, ispgen, NINDEX, ISPEC1A ;
  //  char fnam[] = "set_RATEPAR" ;

  // --------  BEGIN ----------

  // set default RATEPAR
  GENLC.RATEPAR = &INPUTS.RATEPAR ;

  if ( INPUTS.NON1A_MODELFLAG <= 0 ) { return ; }

  // ----------------------------------------
  if ( INDEX_GENMODEL == MODEL_NON1ASED )  {
    
    NINDEX = INP_NON1ASED->NINDEX; 
    if ( NINDEX <= 0 ) { return ; }
  
    // Find sparse non1a index for this 'ilc'
    // Beware of fragile code:  
    //   this isp loop is duplicated later in pick_NON1ASED.
    ispgen = -9 ;
    for ( isp = NINDEX; isp >= 1; isp-- ) {
      MXGEN = INP_NON1ASED->MXGEN[isp];
      if ( ilc <= MXGEN ) { ispgen = isp; }
    }
    ISPEC1A = INP_NON1ASED->ISPEC1A[ispgen];
  }
  else {
    // NON1AGRID model:
    // Use NON1AGRID_RANFlat to determine if random template id will
    // be pecIa. Assumes that NON1A templates have been sorted
    // so that all the PEC1A are at the end of the list.
    double FRAC_PEC1A = INPUTS.RATEPAR_PEC1A.SEASON_FRAC ;
    ISPEC1A = ( GENLC.NON1AGRID_RANFlat  > 1.0 - FRAC_PEC1A  ) ;
  }

  if ( ISPEC1A ) { GENLC.RATEPAR = &INPUTS.RATEPAR_PEC1A ;  }

} // set_RATEPAR




// *****************************************
void genran_modelSmear(void) {

 
  // Oct 2007: 
  // generate random number per filter according for
  // smearing-model. Fill 
  //    GENLC.GENSMEAR_RAN_FILTER[0]
  //    GENLC.GENSMEAR_RAN_FILTER[ifilt_obs]
  // Note that actual mag-smear is computed in genmag_model
  // using the randoms generated here.
  // This function ensure the proper coherence or randomness
  // among filters and epochs.
  //
  //
  // Mar 24 2016: for GENGRID, return after all randoms are initalized to zero.
  // Jul 20 2019: return for repeat strong lens
 
  int    ifilt ;
  int    ILIST_RAN = 2 ; // list to use for genSmear randoms
  double rr8, rho, RHO, rmax, rmin, rtot ;
  //  char fnam[] = "genran_modelSmear" ;

  // -------------- BEGIN --------

  if ( GENSL.REPEAT_FLAG ) { return; } // same randoms for strong lens repeat

  // init all smearing to zero

  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ )  {  
    GENLC.GENSMEAR_RANGauss_FILTER[ifilt] = 0.0 ;  
  }

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return ; }

  // always generate randoms to stay synced, even if mag smear is zero.

  GENRAN_INFO.NSTORE_RAN[ILIST_RAN] = INPUTS.RANLIST_START_GENSMEAR ; // reset

  rmin = INPUTS.SIGMACLIP_MAGSMEAR[0] ;
  rmax = INPUTS.SIGMACLIP_MAGSMEAR[1] ;


  // for global coherent smearing
  
  GENLC.GENSMEAR_RANGauss_FILTER[0] = getRan_GaussClip(1,rmin,rmax);

  // Jun 26 2019: check option for asymmetric smear
  double siglo = (double)INPUTS.GENMAG_SMEAR[0] ;
  double sighi = (double)INPUTS.GENMAG_SMEAR[1] ;
  if ( fabs(siglo-sighi) > 1.0E-6 ) {
    sighi /= siglo ;      siglo /= siglo ;
    GENLC.GENSMEAR_RANGauss_FILTER[0] = getRan_GaussAsym(siglo,sighi, 0.);
  }

  // Now check for model with 100% correlation within a passband,
  // but an independent Gaussian fluctuation in each passband.
  // Note that random # is ADDED here in case both
  // models are specified.
  // Also note that random numbers are always burned
  // to keep randoms synced whether or not this option is used.
  // xxx check to remove ERRSCALE_CORRELATION and change random sync.

  rho  = INPUTS.GENMODEL_ERRSCALE_CORRELATION ;
  RHO  = sqrt( 1.0 - rho*rho );
  
  for ( ifilt=1; ifilt < MXFILTINDX; ifilt++ ) {
    rr8 = getRan_GaussClip(ILIST_RAN,rmin,rmax);
    rtot =  rho * GENLC.GENSMEAR_RANGauss_FILTER[0] +  RHO * rr8 ; 
    GENLC.GENSMEAR_RANGauss_FILTER[ifilt]  = rtot ;      
  }

  if ( !IGNOREFILE(INPUTS.GENMAG_SMEAR_MODELNAME) )  {
    load_genSmear_randoms(GENLC.CID, rmin, rmax, 
			  INPUTS.GENSMEAR_RANGauss_FIX);
  }

  // set randoms for instrinsic scatter matrix (July 2011)
  GENLC.COVMAT_SCATTER_GRAN[0] =  getRan_Gauss(1);
  GENLC.COVMAT_SCATTER_GRAN[1] =  getRan_Gauss(1);
  GENLC.COVMAT_SCATTER_GRAN[2] =  getRan_Gauss(1);
  GEN_COVMAT_SCATTER ( GENLC.COVMAT_SCATTER_GRAN,    // <== input 
		       GENLC.COVMAT_SCATTER );       // <== output
  
  // -----------------------------------
  // pass SN parameters to genSmear function
  double shape, color, z;
  shape = GENLC.SHAPEPAR ;
  z     = GENLC.REDSHIFT_CMB ;
  if ( INDEX_GENMODEL == MODEL_SALT2 ) { color = GENLC.SALT2c ; }
  else                                 { color = GENLC.AV ; }

  SETSNPAR_genSmear(shape, color, z); // Jan 2014

}  // end of genran_modelSmear



// ***********************************
void wr_SIMGEN_FILTERS( char *PATH_FILTERS ) {


  // Created Dec 31, 2010 by R.Kessler
  // Make subdir [VERSION]/[VERSION].FILTERS
  // and write out filter response text files inside.
  // These filter-text files are NOT used by any 
  // SNANA program, but are there for non-SNANA programs.
  //
  // Dec 20 2013: write trans as %.6f instead of %.4f
  // Aug 19 2014: pass output arg PATH_FILTERS

  char cmd[MXPATHLEN], cfilt[4] ;
  char filtFile[MXPATHLEN], filtName[40], surveyName[80]  ;
  double  magprim8, lam8[MXLAMSIM], TransSN8[MXLAMSIM], TransREF8[MXLAMSIM];  
  int  ifilt, ifilt_obs, NLAM, ilam, MASKFRAME, isys ;
  FILE *fp_filt ;

  //  char  fnam[] = "wr_SIMGEN_FILTERS" ;

  // ---------- BEGIN ------------

  sprintf(BANNER,"Write Filter Transmision files " );
  print_banner(BANNER);

  
  // construct path to filter --> set output arg.
  sprintf(PATH_FILTERS, "%s/%s.FILTERS", 
	  PATH_SNDATA_SIM, INPUTS.GENVERSION );
  

  sprintf(cmd,"mkdir -m g+wr %s", PATH_FILTERS );
  isys = system(cmd) ;
  printf("    mkdir %s\n", PATH_FILTERS );
  fflush(stdout);

  MASKFRAME = GENFRAME_OBS ;

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;

    sprintf(filtName,"NULL");
    sprintf(surveyName,"NULL");
    magprim8 = lam8[0] = TransSN8[0] = TransREF8[0] = 0.0 ; 
    NLAM = 0 ;

    // note that both the SN and REF trans are returned,
    // but for simulations they are always the same.

    get_filttrans__(&MASKFRAME, &ifilt_obs, surveyName, filtName,
		    &magprim8, &NLAM, lam8, TransSN8, TransREF8, 80, 40 );

    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
    sprintf(filtFile,"%s/%s.dat", PATH_FILTERS, cfilt );
   
    fp_filt = fopen(filtFile, "wt") ;

    for ( ilam=0; ilam < NLAM; ilam++ ) {
      fprintf(fp_filt,"%7.1f  %.6f \n", lam8[ilam], TransSN8[ilam] );
    }

    fclose(fp_filt);

  }  // end of ifilt loop

  return ;

} // end of wr_SIMGEN_FILTERS

// ***********************************************
void wr_SIMGEN_YAML(SIMFILE_AUX_DEF *SIMFILE_AUX) {
  
  FILE *fp ;
  char *ptrFile  = SIMFILE_AUX->YAML ;
  char fnam[] = "wr_SIMGEN_YAML" ;
  double t_gen   = (TIMERS.t_end - TIMERS.t_end_init); // total time after init

  // ------------ BEGIN ---------------
    if ( (fp = fopen(ptrFile, "wt")) == NULL ) {       
      sprintf ( c1err, "Cannot open SIMGEN YAML file :" );
      sprintf ( c2err," '%s' ", ptrFile );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    fprintf(fp, "SURVEY:          %s\n",    GENLC.SURVEY_NAME   );
    fprintf(fp, "IDSURVEY:        %d\n",    GENLC.IDSURVEY );
    fprintf(fp, "NGENLC_TOT:      %d\n",    NGENLC_TOT     );
    fprintf(fp, "NGENLC_WRITE:    %d\n",    NGENLC_WRITE   );
    fprintf(fp, "NGENSPEC_WRITE:  %d\n",    NGENSPEC_WRITE );
    fprintf(fp, "CPU_MINUTES:     %.2f\n",  t_gen/60.0     );
    fprintf(fp, "ABORT_IF_ZERO:   %d\n",    NGENLC_WRITE   );
    fclose(fp);

    return;

} // end wr_SIMGEN_YAML

// ***********************************************
void wr_SIMGEN_DUMP(int OPT_DUMP, SIMFILE_AUX_DEF *SIMFILE_AUX) {

  /***
   Created Feb 20, 2009 by R.Kessler

   Write SIMGEN variables to column-formated text file
   (same format as fitter "FITRES" file).
   Useful for quick analysis of simgen quantities.

   Example keyword in sim-input file:
 
    SIMGEN_DUMP: 3  Z RA PEAKMAG_g

  will write the redshift and peakmags for each SN to
  [VERSION]_SIMGEN.DUMP

  OPT_DUMP =  1  => init file, write header
  OPT_DUMP >  2  => update 
  OPT_DUMP =  3  => close file (end of job)

       HISTORY

  Mar 16, 2017: write each line to strig SIMFILE_AUX->OUTLINE,
                then make single write per line instead of per value.

  Aug 14 2017: check prescale, and write a little extra info at top of file.

  Dec 01 2017: add hash (#) in front of all comments.

  May 14 2019: free(SIMFILE_AUX->OUTLINE)

  Apr 25 2021: fix bug by exiting after PREP_SIMGEN_DUMP(0).

  ****/

  int   NVAR, ivar, IDSPEC, imjd, index, FIRST ; 

  long long i8, ir8 ;
  int    i4 ;
  float  r4 ; 
  double r8 ;
  char  *ptrFile, *pvar, *str  ;

  FILE *fp ;
  char cval[40] ;
  char fnam[] = "wr_SIMGEN_DUMP" ;

  // --------------- BEGIN ----------

  if ( INPUTS.NVAR_SIMGEN_DUMP < 0 ) return ;

  if ( OPT_DUMP < 1 || OPT_DUMP > 3 ) {
    sprintf ( c1err, "Invalid OPT_DUMP = %d ", OPT_DUMP );
    sprintf ( c2err, "OPT_DUMP must be 1,2, or 3.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  ptrFile = SIMFILE_AUX->DUMP ;

  if ( OPT_DUMP == 1 ) {

    sprintf(BANNER,"Init SIMGEN_DUMP file " );
    print_banner(BANNER );
    
    if ( INPUTS.NVAR_SIMGEN_DUMP == 0 ) 
      { PREP_SIMGEN_DUMP(0); debugexit(fnam); }  // list variables, then quit.
    else
      { PREP_SIMGEN_DUMP(1); }    // prepare valid list

    NVAR = INPUTS.NVAR_SIMGEN_DUMP ; // update NVAR

    // allocate memory to hold one line of output
    // (for faster writing)
    SIMFILE_AUX->OUTLINE = (char *) malloc( 50 + sizeof(char)*NVAR*20 );
    // open file and write header
    if ( (SIMFILE_AUX->FP_DUMP = fopen(ptrFile, "wt")) == NULL ) {       
      sprintf ( c1err, "Cannot open SIMGEN dump file :" );
      sprintf ( c2err," '%s' ", ptrFile );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    printf("\t open %s\n", ptrFile );
    fflush(stdout);

    fp = SIMFILE_AUX->FP_DUMP ;

    // check all user-variables before writing NVAR to dump-file;
    // abort if invalid variable is specified.
    for ( ivar = 0; ivar < NVAR ; ivar++ ) {
      pvar = INPUTS.VARNAME_SIMGEN_DUMP[ivar] ;
      INDEX_SIMGEN_DUMP[ivar] = MATCH_INDEX_SIMGEN_DUMP(pvar);
      if ( strstr(pvar,"WIDTH") ) { GENLC.NWIDTH_SIMGEN_DUMP++; }
    } // end of ivar loop over user variables
    if ( GENLC.NWIDTH_SIMGEN_DUMP>0 ) { init_lightCurveWidth(); }

    // - - - - - - - - - 
    // now write header info to dump file.
    fprintf(fp, "#\n"  );
    fprintf(fp, "# Simulation SUMMARY: one row per event. \n");
    fprintf(fp, "#  MODEL:     %s \n", INPUTS.GENMODEL   );
    fprintf(fp, "#  PRESCALE:  %d \n", INPUTS.PRESCALE_SIMGEN_DUMP);

    fprintf(fp, "#  SELECTION: " );
    if  ( INPUTS.IFLAG_SIMGEN_DUMPALL )
      { fprintf(fp,"NONE (write every generated event)\n");  }
    else
      { fprintf(fp,"Pass Trigger + Cuts\n");  }

    fprintf(fp, "\n") ;

#ifdef TEXTFILE_NVAR
    fprintf(fp, "NVAR: %d \n", NVAR);    
#endif

    sprintf(SIMFILE_AUX->OUTLINE,"VARNAMES: ");
    for ( ivar = 0; ivar < NVAR ; ivar++ ) {
      pvar = INPUTS.VARNAME_SIMGEN_DUMP[ivar] ;
      strcat(SIMFILE_AUX->OUTLINE, " ") ;
      strcat(SIMFILE_AUX->OUTLINE, pvar) ;
    }
    fprintf(fp, "%s\n\n", SIMFILE_AUX->OUTLINE );

    return;  
  } // end of OPT_DUMP==1


  // update dump file for each SN.

  double XN, XNPS=(double)INPUTS.PRESCALE_SIMGEN_DUMP ;

  if ( OPT_DUMP == 2 ) {

    FIRST = (NEVT_SIMGEN_DUMP==0 ) ; // used for SPECTROGRAPH info
    NEVT_SIMGEN_DUMP++ ;  XN=(double)NEVT_SIMGEN_DUMP ;

    // check pre-scale (Aug 2017)
    if ( fmod(XN,XNPS) != 0 ) { return; }

    fp = SIMFILE_AUX->FP_DUMP ;

    // on first event, write comment info for TAKE_SPECTRUM keys.
    // Can't do this during init because GENSPEC arrays not fille until now.
    if ( FIRST ) {
      for(imjd=0; imjd < NPEREVT_TAKE_SPECTRUM; imjd++ ) {
	//	if ( imjd == 0 ) { fprintf(fp,"\n"); }
	IDSPEC = imjd + 1 ; // fortran-like index for ID
	index  = GENSPEC.INDEX_TAKE_SPECTRUM[imjd] ;
	fprintf(fp,"# SPEC%2.2d_T%-4.4s : %6.1f to %6.1f \n",
		IDSPEC, 
		INPUTS.TAKE_SPECTRUM[index].EPOCH_FRAME,
		INPUTS.TAKE_SPECTRUM[index].EPOCH_RANGE[0],
		INPUTS.TAKE_SPECTRUM[index].EPOCH_RANGE[1] ) ;
      }
      fprintf(fp, "\n" );    fflush(fp);
    }


    sprintf(SIMFILE_AUX->OUTLINE, "SN: " );

    NVAR = INPUTS.NVAR_SIMGEN_DUMP ; // update NVAR
    for ( ivar=0; ivar < NVAR; ivar++ ) {

      pvar  = INPUTS.VARNAME_SIMGEN_DUMP[ivar] ;
      index = INDEX_SIMGEN_DUMP[ivar] ;

      if ( index < 0 || index > NVAR_SIMGEN_DUMP ) {
	sprintf(c1err,"invalid index=%d for var='%s' ivar=%d", 
		index, pvar, ivar );
	errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
      }

      r4   = *SIMGEN_DUMP[index].PTRVAL4 ;
      r8   = *SIMGEN_DUMP[index].PTRVAL8 ;
      i4   = *SIMGEN_DUMP[index].PTRINT4 ;
      i8   = *SIMGEN_DUMP[index].PTRINT8 ;
      str  =  SIMGEN_DUMP[index].PTRCHAR ;  // 7.30.2014
      
      if ( r4 != SIMGEN_DUMMY.VAL4 )  
	{ sprintf(cval," %.5le",  r4 ); }

      else if ( r8 != SIMGEN_DUMMY.VAL8 )  { 
	ir8 = (long long)r8 ;

	if ( strstr(pvar,"MJD") != NULL ) 
	  {  sprintf(cval," %.3f", r8 ); }
	else if ( strstr(pvar,"RA") != NULL ) 
	  {  sprintf(cval," %.6f", r8 ); }
	else if ( strstr(pvar,"DEC") != NULL ) 
	  {  sprintf(cval," %.6f", r8 ); }
	else if ( (r8 - ir8) == 0.0 ) // it's really an integer
	  {  sprintf(cval," %lld", ir8 ); }
	else
	  { sprintf(cval," %.5le", r8 );  }
      }
      else if ( i4 != SIMGEN_DUMMY.IVAL4 )  
	{  sprintf(cval," %d",  i4 ); }

      else if ( i8 != SIMGEN_DUMMY.IVAL8 )  
	{  sprintf(cval," %lld",  i8 ); }

      else if ( strcmp(str,SIMGEN_DUMMY.CVAL) != 0 )
	{  sprintf(cval," %s",  str ); }   // 7.30.2014

      else {
	sprintf(c1err,"no value for variable %d (%s)", ivar, pvar);
	errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
      }

      strcat(SIMFILE_AUX->OUTLINE,cval);

    } // end of ivar loop

    fprintf(fp, "%s\n", SIMFILE_AUX->OUTLINE );
    fflush(fp);

  } // end of OPT_DUMP=2 if-block


  if ( OPT_DUMP == 3 ) {
    free(SIMFILE_AUX->OUTLINE);
    fclose(SIMFILE_AUX->FP_DUMP);
    printf("  %s\n", ptrFile );
    return ;
  }


} // end of wr_SIMGEN_DUMP


// ******************************************
int MATCH_INDEX_SIMGEN_DUMP(char *varName ) {

  // Mar 16 2014
  // Return SIMGEN_DUMP index corresponding to input *varName.
  // First check for exact match; if none found then check for
  // case-insensitive match. This logic allows for similar
  // strings such as SNRMAX_a and SNRMAX_A, yet also allows
  // for user-sloppiness; e.g., snrmax instead of SNRMAX.
  //
  // If there are 0 matches, or more than 1 match, ABORT.
  //

  int index, MATCH_INDEX, jcmp ;
  int NMATCH, NMATCH_EXACT, NMATCH_ignoreCase ;
  char fnam[] = "MATCH_INDEX_SIMGEN_DUMP" ;

  // ----------- BEGIN -----------

  NMATCH = NMATCH_EXACT = NMATCH_ignoreCase = 0 ;
  MATCH_INDEX = -9 ;

  // first try for exact match:
  for ( index=0; index < NVAR_SIMGEN_DUMP; index++ ) {
    jcmp = strcmp(varName,SIMGEN_DUMP[index].VARNAME) ;
    if ( jcmp == 0 ) {  NMATCH_EXACT++ ;  MATCH_INDEX = index ; }
  } 
  NMATCH = NMATCH_EXACT ;


  // if we don't find exact match, look for case-insensitive match
  if ( NMATCH_EXACT == 0 ) {
    for ( index=0; index < NVAR_SIMGEN_DUMP; index++ ) {
      jcmp = strcmp_ignoreCase(varName,SIMGEN_DUMP[index].VARNAME) ; 
      if ( jcmp == 0 ) {  NMATCH_ignoreCase++ ; MATCH_INDEX = index ;  }
    } 
    NMATCH = NMATCH_ignoreCase ;    
  }

  // ------------------------
  // Error checking
  
  if ( NMATCH == 0 ) {
    print_preAbort_banner(fnam);
    PREP_SIMGEN_DUMP(0); // print list of valid varnames, then quit
    sprintf(c1err,"Undefined SIMGEN_DUMP variable: '%s'", varName);
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }
  
  if ( NMATCH > 1 ) {
    print_preAbort_banner(fnam);
    PREP_SIMGEN_DUMP(0); // print list of valid varnames, then quit
    sprintf(c1err,"SIMGEN_DUMP variable '%s'", varName);
    sprintf(c2err,"defined %d times ??", NMATCH);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  if ( MATCH_INDEX < 0 ) {
    sprintf(c1err,"MATCH_INDEX = %d  for  varName = '%s'", 
	    MATCH_INDEX, varName);
    sprintf(c2err,"is really messed up.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  
  // if we get here, return valid matched index.
  return MATCH_INDEX ;

} // end of MATCH_INDEX_SIMGEN_DUMP

// ******************************************
void PREP_SIMGEN_DUMP(int OPT_DUMP) {

  //
  // Created May 9, 2009 by R.Kessler
  // Prepare list of all potential variable, and assign pointer. 
  // Note that this list is hard-wired, and that some variables
  // have more than one name (i.e, T0 and PEAKMJD are the same)
  // 
  // OPT_DUMP = 0 => dump variable list, return
  // OPT_DUMP = 1 => prepare list (no output), then return

  // Nov 2, 2009: remove mlcs/salt2 restrictions so that we can
  //              combine SIMGEN_DUMP files
  //
  // Jul 6, 2010; add EFFMASK -> &GENLC.CUTBIT_MASK ;
  // Nov 10 2010: add DM15
  // Jun 28,2011: add ZSMEAR = observed redshift
  //
  // Feb 12, 2014: allow any variable on HOSTLIB_STOREVAR list
  // Feb 18, 2014: add MAGSMEAR_COH
  // Feb 25, 2014: allow SALT2 or S2 prefix for SALT2 params
  // Mar 16, 2014: fix SIMSED logic to avoud multiple entries with
  //               same variables; see SKIP1 logic.
  //
  // Jul 30 2014: add "FIELD" string and init SIMGEN_DUMMY.CVAL
  //
  // Jan 29 2015: allow RA_SN and DEC_SN to distinguish host coords.
  //
  // Dec 27 2015: add MJD_TRIGGER
  // Jun 23 2016: add SBMAG_[band]
  // Mar 01 2017: allow AV & RV for SIMSED model.
  // Mar 14 2017: tack on TAKE_SPECTRUM info
  // mar 28 2017: add IDSURVEY
  // Oct 16 2019: move MAGSMEAR_COH after SKIP1 
  // Apr 28 2020: allow list of var names for SALT2c,x1,x0 (see strList_)
  // Jul 24 2020: for OPT_DUMP=0, do NOT quit

  int i, ifilt, ifilt_obs, ifilt_rest, ipar, imap, ivar, NTMP ;
  char *cptr ;
  char fnam[] = "PREP_SIMGEN_DUMP";

  // --------------- BEGIN -------------

  NVAR_SIMGEN_DUMP = 0 ;
  NEVT_SIMGEN_DUMP = 0 ;

  SIMGEN_DUMMY.VAL4  = -987.0 ;
  SIMGEN_DUMMY.VAL8  = -987.0 ;
  SIMGEN_DUMMY.IVAL4 = -987   ;
  SIMGEN_DUMMY.IVAL8 = -987   ;
  sprintf(SIMGEN_DUMMY.CVAL,"xxx");

  // init pointers to DUMMY values
  for ( i=0; i< MXSIMGEN_DUMP; i++ ) {
    SIMGEN_DUMP[i].PTRVAL4 = &SIMGEN_DUMMY.VAL4 ;
    SIMGEN_DUMP[i].PTRVAL8 = &SIMGEN_DUMMY.VAL8 ;
    SIMGEN_DUMP[i].PTRINT4 = &SIMGEN_DUMMY.IVAL4 ;
    SIMGEN_DUMP[i].PTRINT8 = &SIMGEN_DUMMY.IVAL8 ;
    SIMGEN_DUMP[i].PTRCHAR =  SIMGEN_DUMMY.CVAL ;
  }


  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"CID");
  if ( WRFLAG_CIDRAN == 0 )  
    {  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.CID ;  }
  else  
    {  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.CIDRAN ; }
  NVAR_SIMGEN_DUMP++ ;


  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"LIBID");
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.SIMLIB_ID ;
  NVAR_SIMGEN_DUMP++ ;


  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NGEN_LIBID");
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.NGEN_SIMLIB_ID ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"IDSURVEY");
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.IDSURVEY ;
  NVAR_SIMGEN_DUMP++ ;

  // --------- coords ----------
  // Allow _SN suffix to avoid collision with host coords.
  // Also, allow either DEC or DECL.

  char RASTRINGS[2][12] = { "RA", "RA_SN" } ;
  for(i=0; i < 2 ; i++ )  {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr, "%s", RASTRINGS[i]);
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.RA ;
    NVAR_SIMGEN_DUMP++ ;
  }

  char DECSTRINGS[4][12] = { "DECL", "DEC", "DECL_SN", "DEC_SN" } ;
  for(i=0; i<4; i++ )  {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr, "%s", DECSTRINGS[i]);
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.DEC ;
    NVAR_SIMGEN_DUMP++ ;
  }

  // random coord shift (Apr 2021)
  if ( INPUTS.MXRADIUS_RANDOM_SHIFT > 1.0E-12 ) {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"shift_RA"); 
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.random_shift_RA ;
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"shift_DEC"); 
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.random_shift_DEC ;
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"shift_RADIUS"); 
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.random_shift_RADIUS ;
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"shift_PHI"); 
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.random_shift_PHI ;
    NVAR_SIMGEN_DUMP++ ;
  }


  // ---------
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"FIELD");  
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRCHAR = SIMLIB_HEADER.FIELD ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NFIELD_OVP");  
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &SIMLIB_HEADER.NFIELD_OVP ;
  NVAR_SIMGEN_DUMP++ ;

  // Galactic  extinction
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MWEBV");  // simulated MWEBV applied to SN
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.MWEBV_SMEAR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MWEBVMAP");  // true MWEBV from map
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.MWEBV ;
  NVAR_SIMGEN_DUMP++ ;


  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MWEBV_TRUE");  // true MWEBV from map
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.MWEBV ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MWEBVERR");  // reported MWEBV error for fitting
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.MWEBV_ERR ;
  NVAR_SIMGEN_DUMP++ ;

  // redshift variables
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"ZCMB_SMEAR");  // reported redshift in data file
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.REDSHIFT_CMB_SMEAR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"ZCMB");  // true zcmb
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.REDSHIFT_CMB ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"ZHELIO");  // true ZHELIO
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.REDSHIFT_HELIO ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"ZFLAG");  // integer ZFLAG 
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.REDSHIFT_FLAG ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"VPEC"); // true VPEC
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.VPEC ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"VPEC_SMEAR");  // measured VPEC
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.VPEC_SMEAR ;
  NVAR_SIMGEN_DUMP++ ;

  // --------------------------------------
  // legacy redshift variables to avoid aborting on old input files
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"Z");
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.REDSHIFT_CMB ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GENZ");
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.REDSHIFT_CMB ;
  NVAR_SIMGEN_DUMP++ ;
  // --------------------------------------

  // distance modulus
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MU");
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.DLMU ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"DLMAG") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.DLMU ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"LENSDMU") ; // from weak lensing
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.LENSDMU ;
  NVAR_SIMGEN_DUMP++ ;

  // strong lens magnification 
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SL_MAGSHIFT") ; // from strong lens
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SL_MAGSHIFT ;
  NVAR_SIMGEN_DUMP++ ;
  // ... or ...
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"STRONGLENS_MAGSHIFT") ; // from strong lens
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SL_MAGSHIFT ;
  NVAR_SIMGEN_DUMP++ ;
  
  // host Z stuff
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALID") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT8 = &SNHOSTGAL.GALID ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALZTRUE") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.ZTRUE ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALZPHOT") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.ZPHOT ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALZERR") ;  // legacy name for ZPHOT ERROR
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.ZPHOT_ERR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALZPHOTERR") ;  // ZPHOT ERROR
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.ZPHOT_ERR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALSNSEP") ;   // host-SN separation, arcsec
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL_DDLR_SORT[0].SNSEP ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALSNDLR") ;   // 2/2019: DLR from Sako 2014, Gupta 2016
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL_DDLR_SORT[0].DLR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALSNDDLR") ;   //2/2019:  d_DLR = SNSEP/DLR
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL_DDLR_SORT[0].DDLR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALZDIF") ;   // zSN - zGAL, Nov 2015
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.ZDIF ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALSNDM") ; // SN mag shift
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.WGTMAP_SNMAGSHIFT ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALWGT") ; 
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.WGTMAP_WGT ;
  NVAR_SIMGEN_DUMP++ ;

  // allow any variable used in HOSTLIB_WGTMAP
  for ( imap=0; imap < HOSTLIB_WGTMAP.GRIDMAP.NDIM; imap++ ) {   
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr, "%s", HOSTLIB_WGTMAP.VARNAME[imap]) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.WGTMAP_VALUES[imap] ;
    NVAR_SIMGEN_DUMP++ ;
  }

  // allow any variable on extra HOSTLIB_STOREVAR list
  NTMP = NVAR_SIMGEN_DUMP ; 
  for(ivar=0; ivar < HOSTLIB_OUTVAR_EXTRA.NOUT; ivar++ ) {
    if ( HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[ivar] ) { continue; } 
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;

    sprintf(cptr,"%s", HOSTLIB_OUTVAR_EXTRA.NAME[ivar] );
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &HOSTLIB_OUTVAR_EXTRA.VALUE[ivar][0];
    NVAR_SIMGEN_DUMP++ ;
  } // end of ivar loop


  // Aug 2014: allow host mags and SBFLUX

  if ( INPUTS.HOSTLIB_USE ) {
    for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt]; 
      if ( ifilt > 20 ) { continue; } // avoid array-bound overflow
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr,"HOSTMAG_%c", FILTERSTRING[ifilt_obs] ) ;
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.GALMAG[ifilt_obs][0];
      NVAR_SIMGEN_DUMP++ ;
    }


    for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt]; 
      if ( ifilt > 20 ) { continue; } // avoid array-bound overflow
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr,"SBFLUX_%c", FILTERSTRING[ifilt_obs] ) ;
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.SB_FLUXCAL[ifilt_obs];
      NVAR_SIMGEN_DUMP++ ;
    }

    for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt]; 
      if ( ifilt > 20 ) { continue; } // avoid array-bound overflow
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr,"SBMAG_%c", FILTERSTRING[ifilt_obs] ) ;
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.SB_MAG[ifilt_obs] ;
      NVAR_SIMGEN_DUMP++ ;
    }

  }

  // -----------
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"PEAKMJD") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.PEAKMJD ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"DTSEASON_PEAK") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.DTSEASON_PEAK ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"PEAKMJD_SMEAR") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.PEAKMJD_SMEAR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MJD0") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.PEAKMJD ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MJD_TRIGGER") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.MJD_TRIGGER ;
  NVAR_SIMGEN_DUMP++ ;

  // write true/calculated peakMag (not observed peak mag)
  for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt]; 
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"MAGT0_%c", FILTERSTRING[ifilt_obs] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.peakmag_obs[ifilt_obs] ;
    NVAR_SIMGEN_DUMP++ ;
  }

  // Aug 2017: allow traditional PEAKMAG_x for true peakmag (not observed)
  for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt]; 
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"PEAKMAG_%c", FILTERSTRING[ifilt_obs] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.peakmag_obs[ifilt_obs] ;
    NVAR_SIMGEN_DUMP++ ;
  }


  // Aug 2017: true LC width (days)
  GENLC.NWIDTH_SIMGEN_DUMP = 0;
  for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt]; 
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"WIDTH_%c", FILTERSTRING[ifilt_obs] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.WIDTH[ifilt_obs] ;
    NVAR_SIMGEN_DUMP++ ;
  }


  if ( GENFRAME_OPT == GENFRAME_REST ) {
    for ( ifilt=0; ifilt < GENLC.NFILTDEF_REST; ifilt++ ) {
      ifilt_rest = GENLC.IFILTMAP_REST[ifilt]; 
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr,"magT0_%c", FILTERSTRING[ifilt_rest] ) ;
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8=&GENLC.peakmag_rest[ifilt_rest];
      NVAR_SIMGEN_DUMP++ ;
    }
  }


  if ( INDEX_GENMODEL == MODEL_SIMSED && INPUTS.WRFLAG_MODELPAR ) {
    for( ipar=0;  ipar < INPUTS.NPAR_SIMSED; ipar++ ) {
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr, "%s", INPUTS.PARNAME_SIMSED[ipar] ) ;
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SIMSED_PARVAL[ipar];
      NVAR_SIMGEN_DUMP++ ; 
    }
  }


  if ( IS_PySEDMODEL ) {
    for( ipar=0 ;  ipar < Event_PySEDMODEL.NPAR; ipar++ ) {
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr, "%s", Event_PySEDMODEL.PARNAME[ipar] );
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &Event_PySEDMODEL.PARVAL[ipar]; 
      NVAR_SIMGEN_DUMP++ ; 
    }
  }

  if ( INDEX_GENMODEL == MODEL_LCLIB  && INPUTS.WRFLAG_MODELPAR) {
    for( ipar=0 ;  ipar <= LCLIB_INFO.NPAR_MODEL; ipar++ ) {
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr, "%s", LCLIB_INFO.PARNAME_MODEL[ipar] );
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = 
	&LCLIB_EVENT.PARVAL_MODEL[ipar]; 
      NVAR_SIMGEN_DUMP++ ; 
    }


    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"GLAT");  
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.GLAT ;
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"b");  
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.GLAT ;
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"GLON");  
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.GLON ;
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"l");  
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.GLON ;
    NVAR_SIMGEN_DUMP++ ;
  }


  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,PARNAME_AV) ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.AV ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,PARNAME_RV) ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.RV ;
  NVAR_SIMGEN_DUMP++ ;

  if ( INDEX_GENMODEL == MODEL_SIMSED ) { goto SKIP1; }

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"DELTA") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.DELTA ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"DM15") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.DM15 ;
  NVAR_SIMGEN_DUMP++ ;

  // add STRETCH for B18-style snoopy model that uses stretch
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"STRETCH") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.STRETCH ;
  NVAR_SIMGEN_DUMP++ ;

  char strList_alpha[3][20] = { "SALT2alpha", "S2alpha", "SIM_alpha" };
  for(i=0; i < 3; i++ ) {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s", strList_alpha[i] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2alpha ;
    NVAR_SIMGEN_DUMP++ ;
  }

  char strList_beta[3][20] = { "SALT2beta", "S2beta", "SIM_beta" };
  for(i=0; i < 3; i++ ) {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s", strList_beta[i] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2beta ;
    NVAR_SIMGEN_DUMP++ ;
  }


  char strList_x0[3][20] = { "S2x0", "SALT2x0", "SIM_x0" };
  for(i=0; i < 3; i++ ) {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s", strList_x0[i] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2x0 ;
    NVAR_SIMGEN_DUMP++ ;
  }

  char strList_x1[3][20] = { "S2x1", "SALT2x1", "SIM_x1" };
  for(i=0; i < 3; i++ ) {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s", strList_x1[i] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2x1 ;
    NVAR_SIMGEN_DUMP++ ;
  }

  char strList_c[3][20] = { "S2c", "SALT2c", "SIM_c" };
  for(i=0; i < 3; i++ ) {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s", strList_c[i] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2c ;
    NVAR_SIMGEN_DUMP++ ;
  }

  char strList_mb[5][20] = { "S2mb", "SALT2mb", "SALT2mB", "SIM_mb", "SIM_mB" };
  for(i=0; i < 5; i++ ) {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s", strList_mb[i] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2mB ;
    NVAR_SIMGEN_DUMP++ ;
  }


  // check COVMAT_SCATTER
  if ( INPUTS.NCOVMAT_SCATTER > 0 ) {
    int iscat;
    for ( iscat=0; iscat < 3; iscat++ ) {
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr,"SCAT%s", GENLC.COVMAT_SCATTER_NAME[iscat] );
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.COVMAT_SCATTER[iscat];
      NVAR_SIMGEN_DUMP++ ;   
    }
  }


  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SALT2gammaDM") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2gammaDM ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SNMAGSHIFT") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2gammaDM ;
  NVAR_SIMGEN_DUMP++ ;

  // - - - - - - - - - - - - - - - -
  // params from VCR intrinsic scatter model
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"VSI") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENSMEAR_VCR.GENRAN_VSI ;
  NVAR_SIMGEN_DUMP++ ;

  for(i=0; i<5; i++ ) {
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"COLORSHIFT%d", i) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENSMEAR_VCR.GENRAN_COLORSHIFT[i];
    NVAR_SIMGEN_DUMP++ ;
  }

  // - - - - - - - - - - - - - - - -
  // rise and fall shifts
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TRISE_SHIFT") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.RISETIME_SHIFT ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TFALL_SHIFT") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.FALLTIME_SHIFT ;
  NVAR_SIMGEN_DUMP++ ;

 SKIP1:
  
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MAGSMEAR_COH") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.MAGSMEAR_COH[0] ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"MAGSMEAR_COH2") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.MAGSMEAR_COH[1] ;
  NVAR_SIMGEN_DUMP++ ;

  // typings
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GENTYPE") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.SIMTYPE ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SNTYPE") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.SNTYPE ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TYPE") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.SNTYPE ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NON1A_INDEX") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.TEMPLATE_INDEX ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NONIA_INDEX") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.TEMPLATE_INDEX ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NOBS") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.NOBS ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NOBS_UNDEFINED") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.NOBS_UNDEFINED ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NOBS_SATURATE") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.NOBS_SATURATE[1] ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NOBS_NOSATURATE") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.NOBS_SATURATE[0] ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NOBSDIF") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.NOBS_MJDDIF ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"NEPOCH") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.NEPOCH ;
  NVAR_SIMGEN_DUMP++ ;

  // --------------------------------------------------------
  //  TAKE_SPECTRUM info: TOBS, SNR, TEXPOSE, LAMMIN, LAMMAX
  // --------------------------------------------------------
  PREP_SIMGEN_DUMP_TAKE_SPECTRUM();

  // ----------------------------------------------
  //
  // Start analysis variables (i.e,, not used in generation)
  //
  // -----------------------------------------------
  NVAR_SIMGEN_DUMP_GENONLY =   NVAR_SIMGEN_DUMP ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"REDSHIFT_FINAL") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.REDSHIFT_CMB_SMEAR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TRESTMIN") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENLC.TRESTMIN ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TMIN") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENLC.TRESTMIN ;  
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TRESTMAX") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENLC.TRESTMAX ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TMAX") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENLC.TRESTMAX ;
  NVAR_SIMGEN_DUMP++ ;
  // ---

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SNRMAX") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SNRMAX_SORTED[1] ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SNRMAX1") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SNRMAX_SORTED[1] ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SNRMAX2") ;  // 2nd SNRMAX among filters
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SNRMAX_SORTED[2] ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SNRMAX3") ;  // 3rd SNRMAX among filters
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SNRMAX_SORTED[3] ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SNRMAX4") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SNRMAX_SORTED[4] ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SNRMAX5") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SNRMAX_SORTED[5] ;
  NVAR_SIMGEN_DUMP++ ;
  // ----- SNRMAX by filter
  for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt]; 
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"SNRMAX_%c", FILTERSTRING[ifilt_obs] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SNRMAX_FILT[ifilt_obs] ;
    NVAR_SIMGEN_DUMP++ ;
  }

  // -----------

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TGAPMAX") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENLC.TGAPMAX ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"T0GAPMAX") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENLC.T0GAPMAX ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"TIME_ABOVE_SNRMIN") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENLC.TIME_ABOVE_SNRMIN ;
  NVAR_SIMGEN_DUMP++ ;


  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"CUTMASK") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.CUTBIT_MASK ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SIM_EFFMASK") ;  // legacy name
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.SEARCHEFF_MASK ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SIM_SEARCHEFF_MASK") ;  // repeat to match manual
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRINT4 = &GENLC.SEARCHEFF_MASK ;
  NVAR_SIMGEN_DUMP++ ;

  if ( NVAR_SIMGEN_DUMP >= MXSIMGEN_DUMP ) {
    sprintf(c1err,"NVAR_SIMGEN_DUMP=%d exceeds array bound of %d.",
	    NVAR_SIMGEN_DUMP, MXSIMGEN_DUMP );
    errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
  }

  //check option to LIST variables and quit
  if ( OPT_DUMP == 0 ) {

    printf("\n\n");
    print_banner("Allowed SIMGEN_DUMP variables \n");

    printf(" From SIM-GENERATION: \n");
    for( i = 0; i < NVAR_SIMGEN_DUMP_GENONLY; i++ ) 
      { printf(" %s", SIMGEN_DUMP[i].VARNAME ) ; }

    printf("\n\n From CUTWIN-ANALYSIS: \n");
    for( i = NVAR_SIMGEN_DUMP_GENONLY; i <NVAR_SIMGEN_DUMP; i++ ) 
      {  printf(" %s", SIMGEN_DUMP[i].VARNAME ) ; }

    printf("\n\n");
    printf("\t Example of sim-input file syntax is \n");
    printf("\t SIMGEN_DUMP:  5  CID Z RA DEC SNRMAX \n");
  }

  return;

} // end of PREP_SIMGEN_DUMP


// **********************************************
void PREP_SIMGEN_DUMP_TAKE_SPECTRUM(void) {

  // Mar 2017
  // if there are TAKE_SPECTRUM keys defined by SNR_LAMREST[OBS],
  // then automatically append info to SIMGEN_DUMP file so that
  // user only needs to include the usual SIMGEN_DUMP key in the
  // sim-input file. User does NOT add TAKE_SPECTRUM variables
  // since these are added here by default.

  int i, imjd, IDSPEC, OPT_TEXPOSE, NOPT_SNR=0;
  int idump_inp, idump_list;
  int NVAR_SIMGEN_DUMP_START = NVAR_SIMGEN_DUMP ;
  int NVAR_SIMGEN_DUMP_END   = NVAR_SIMGEN_DUMP ;
  char PREFIX[12], *cptr ;
  //  char fnam[] = "PREP_SIMGEN_DUMP_TAKE_SPECTRUM" ;

  // ------------- BEGIN ----------

  // logic is tricky because GENSPEC.OPT_TEXPOSE is not
  // set until first event, and we don't know the mapping
  // between INPUTS.TAKE_SPECTRUM and GENSPEC. 
  // If any TAKE_SPECTRUM key uses SNR option, then add
  // all TAKE_SPECTRUM epochs to list.

  // begin by counting how many TAKE_SPECTRUM keys use
  // the SNR option
  for(i=0; i < NPEREVT_TAKE_SPECTRUM ; i++ ) {
    OPT_TEXPOSE = INPUTS.TAKE_SPECTRUM[i].OPT_TEXPOSE;
    if (OPT_TEXPOSE==2 ) { NOPT_SNR++; }
  }

  if ( NOPT_SNR == 0 ) { return ; }

  for(i=0; i < NPEREVT_TAKE_SPECTRUM ; i++ ) {

    imjd   = i ;  
    IDSPEC = imjd + 1;  // start count at 1 to match data files

    sprintf(PREFIX,"SPEC%2.2d", IDSPEC);
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s_TREST", PREFIX) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENSPEC.TREST_LIST[imjd] ;    
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s_TEXPOSE", PREFIX) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENSPEC.TEXPOSE_LIST[imjd] ; 
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s_SNR_REQUEST", PREFIX) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENSPEC.SNR_REQUEST_LIST[imjd] ; 
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s_SNR_COMPUTE", PREFIX) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL4 = &GENSPEC.SNR_COMPUTE_LIST[imjd] ; 
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s_LAMOBS_MIN", PREFIX) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = 
      &GENSPEC.LAMOBS_SNR_LIST[imjd][0] ;
    NVAR_SIMGEN_DUMP++ ;

    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s_LAMOBS_MAX", PREFIX) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = 
      &GENSPEC.LAMOBS_SNR_LIST[imjd][1] ; 
    NVAR_SIMGEN_DUMP++ ;
  } // end imjd loop

  // finally, tack on the template exposure time.

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SPEC_TEXPOSE_TEMPLATE" ) ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENSPEC.TEXPOSE_TEMPLATE ;
  NVAR_SIMGEN_DUMP++ ;

  NVAR_SIMGEN_DUMP_END = NVAR_SIMGEN_DUMP;
  

  // add these to SIMGEN_DUMP list as if the user had
  // added them after the SIMGEN_DUMP key.

  if ( NVAR_SIMGEN_DUMP_END > NVAR_SIMGEN_DUMP_START ) {
    for ( idump_list = NVAR_SIMGEN_DUMP_START-1 ; 
	  idump_list<  NVAR_SIMGEN_DUMP_END  ;
	  idump_list++ ) {
      idump_inp = INPUTS.NVAR_SIMGEN_DUMP ;
      sprintf(INPUTS.VARNAME_SIMGEN_DUMP[idump_inp], "%s", 
	      SIMGEN_DUMP[idump_list].VARNAME);
      INPUTS.NVAR_SIMGEN_DUMP++ ; 
    }
  }

  return ;

} // end PREP_SIMGEN_DUMP_TAKE_SPECTRUM

// **********************************************
int doReject_SIMGEN_DUMP(char *rejectStage) {

  // Feb 2014
  // Called for SN that fail
  //    rejectStage = 'SEARCHEFF'   or
  //    rejectStage = 'CUTWIN'
  // to decide if this SN should be updated in the SIMGEN_DUMP.
  // Default is to update SIMGEN_DUMP only for SN that are 
  // written to the data files, but check for options to 
  // include some or all of the rejected SN.
  //

  int DOFLAG, ISCUTWIN ;

  // ------------------ BEGIN ---------

  DOFLAG = 0 ; // default is do NOT udpdate SIMGEN_DUMP file.

  // check option to write ALL rejected SN.
  if ( INPUTS.IFLAG_SIMGEN_DUMPALL ) { DOFLAG = 1; }


  // check option to write SN that fail CUTWIN
  ISCUTWIN =  ( strcmp(rejectStage,"CUTWIN") == 0 ) ;
  if ( ISCUTWIN && INPUTS.APPLY_CUTWIN_OPT == 3 ) { DOFLAG = 1; }

  return DOFLAG ;

} // end of  doReject_SIMGEN_DUMP

// ******************************
void dmp_trace_main(char *string, int ilc) {
  printf(" >>>>> TRACE_MAIN-%s  at ilc=%d <<<<<< \n", string, ilc );
  fflush(stdout);
} // end of dmp_trace_main




// ************************************************
void dmp_event(int ilc) {

  // debug option to dump event.
  // Jan 9 2016: include ilc cut here instead of in main.

  int i, ifilt_obs, ifilt_rest1, ifilt_rest2;

  // -------- BEGIN ---------

  if ( ilc < INPUTS.GENRANGE_DMPEVENT[0] ) { return ; }
  if ( ilc > INPUTS.GENRANGE_DMPEVENT[1] ) { return ; }


  sprintf(BANNER, " Dump Gen-Event %d ", ilc);
  print_banner ( BANNER );

  printf(" RA=%f, DEC=%f,  PEAKMJD=%8.1f, REDSHIFT(CMB)=%6.4f \n",
	 GENLC.RA, GENLC.DEC, GENLC.PEAKMJD, GENLC.REDSHIFT_CMB );


  printf("  MJD      Trest     restmag1 restmag2     obsmag zpt \n");
  for ( i=1; i<= GENLC.NEPOCH; i++ ) {

    ifilt_obs   = GENLC.IFILT_OBS[i] ;

    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

    ifilt_rest1 = GENLC.IFILTMAP_REST1[ifilt_obs];
    ifilt_rest2 = GENLC.IFILTMAP_REST2[ifilt_obs];

    printf(" %9.3f %6.2f   %6.2f(%c) %6.2f(%c)  %6.2f(%c) %6.2f \n",
	   GENLC.MJD[i]
	   ,GENLC.epoch_rest[i] 
	   ,GENLC.genmag_rest[i],  FILTERSTRING[ifilt_rest1]
	   ,GENLC.genmag_rest2[i], FILTERSTRING[ifilt_rest2]
	   ,GENLC.genmag_obs[i],   FILTERSTRING[ifilt_obs]
	   ,SIMLIB_OBS_GEN.ZPTADU[i]	   );
  }


}  // end of dmp_event


// ***********************************
double gen_peakmjd(void) {

  // Jan 2012.
  // return PEAKMJD randomly selected within user range.
  // Avoid optional GENSKIP_PEAKMJD gaps from SIMLIB header.
  //
  //
  // May 2014: 
  //   float rangen -> double getRan_Flat  
  //   float GENRANGE_PEAKMJD -> double GENRANGE_PEAKMJD
  //   float PKMJD, MJD[2] -> double
  //
  // Sep 17 2015: check option to snap to nearest PEAKMJD read from file.
  // Jan 25 2017: if last PEAKMJD was read from SIMLIB header, then
  //              just return current PEAKMJD in case SIMLIB_NREPEAT
  //              key is set.
  //

  double PKMJD=-9.0 , MJD[2];
  int    NSKIP_RANGE, NSKIP_MJD, i ;
  char   fnam[] = "gen_peakmjd" ;

  // ------------- BEGIN --------------

  if ( INDEX_GENMODEL == MODEL_SIMLIB ) { return(PKMJD); }

  if ( GENLC.ISOURCE_PEAKMJD == ISOURCE_PEAKMJD_SIMLIB ) 
    { return(GENLC.PEAKMJD); }

  NSKIP_RANGE   = SIMLIB_GLOBAL_HEADER.NGENSKIP_PEAKMJD ;
  NSKIP_MJD     = 0 ;
  
 PICKRAN_PEAKMJD:
  
  if ( NSKIP_MJD > 1000 ) {
    sprintf(c1err,"Could not find valid PEAKMJD after 1000 tries");
    sprintf(c2err,"Check GENRANGE_PEAKMJD and GENSKIP_PEAKMJD");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  PKMJD = getRan_Flat (1,INPUTS.GENRANGE_PEAKMJD );

  // check for MJD-ranges to skip
  for ( i=0; i < NSKIP_RANGE ; i++ ) {
    MJD[0] = SIMLIB_GLOBAL_HEADER.GENSKIP_PEAKMJD[i][0] ;
    MJD[1] = SIMLIB_GLOBAL_HEADER.GENSKIP_PEAKMJD[i][1] ;
    if ( PKMJD > MJD[0] && PKMJD < MJD[1] )  { 
      NSKIP_MJD++ ;
      goto PICKRAN_PEAKMJD ; 
    }
  }

  GENLC.ISOURCE_PEAKMJD = ISOURCE_PEAKMJD_RANDOM ;
  return PKMJD ;

}   // end of gen_peakmjd


double gen_peakmjd_smear(void) {

  // May 2019
  // Determine PEAKMJD estimate for data file.
  // Either Gaussian smear, or search for max-flux.
  // May 31 2019: generate random Gaussian here.

  double PEAKMJD_SMEAR = GENLC.PEAKMJD; // default
  double smear;
  int    o=-9, obs_atFLUXMAX[MXFILTINDX];
  //  char fnam[] = "gen_peakmjd_smear" ;

  // ----------------- BEGIN -----------------

  // always burn Gaussian random, regardless of option.
  GENLC.PEAKMJD_RANGauss = getRan_Gauss(1); 

  // check option of Gaussian smear.
  // Note that PEAKMJD_RANGauss is selected earlier in gen_event_driver()
  // in order to preserve random sync, but at some point PEAKMJD_RANGauss
  // should be generated here.
  if ( INPUTS.GENSIGMA_PEAKMJD > 1.0E-9 ) {
    smear = GENLC.PEAKMJD_RANGauss * INPUTS.GENSIGMA_PEAKMJD;
    PEAKMJD_SMEAR = GENLC.PEAKMJD + smear ;
  }

  // check option to associate MJD  with epoch at max flux.
  // Note that get_obs_atFLUXMAX is also used by fitting
  // programs.
  if ( INPUTS.OPT_SETPKMJD > 0 ) {
    int NOBS = GENLC.NOBS ;
    get_obs_atFLUXMAX(SNDATA.CCID, NOBS, 
		      &SNDATA.FLUXCAL[1], &SNDATA.FLUXCAL_ERRTOT[1],
		      &SNDATA.MJD[1], &GENLC.IFILT_OBS[1],
		      obs_atFLUXMAX ) ;
    o = obs_atFLUXMAX[0] ; 
    if ( o >= 0 ) 
      { PEAKMJD_SMEAR = SNDATA.MJD[o+1] ; }
    else
      { PEAKMJD_SMEAR = -9.0 ; }

  }

  int LDMP = (GENLC.CID == -4530) ;
  if ( LDMP ) {
    printf(" xxx ---------------------------------------- \n");
    printf(" xxx GENSIGMA_PEAKMJD=%.2f, RANGauss=%.3f, OPT_SETPKMJD=%d\n",
	   INPUTS.GENSIGMA_PEAKMJD, GENLC.PEAKMJD_RANGauss, 
	   INPUTS.OPT_SETPKMJD);
    printf(" xxx CID=%d: PEAKMJD(true,smear)=%.1f, %.1f  (o=%d)\n",
	   GENLC.CID, GENLC.PEAKMJD, PEAKMJD_SMEAR, o); fflush(stdout);
  }

  return(PEAKMJD_SMEAR);

} // end gen_peakmjd_smear

// *********************************
double gen_redshift_cmb ( void) {


  /************
  Generate cmb redshift from Hubble or from data distribution
  Return cmb redshift.

  Sep 23, 2011: float -> double 
  Dec 21, 2015: pass zmin & zmax to genz_hubble()
  Apr 12, 2017: remove dlmag arg.

  *****************/

  double zmin = INPUTS.GENRANGE_REDSHIFT[0];
  double zmax = INPUTS.GENRANGE_REDSHIFT[1];
  double zcmb = -9.0 ; 
  //  char fnam[] = "gen_redshift_cmb" ;

  // --------------- BEGIN ------------

  if ( INDEX_GENMODEL == MODEL_SIMLIB ) { return(zcmb); }
  if ( INDEX_GENMODEL == MODEL_LCLIB  ) { return(zmin); }

  // check for delta-function in redshift
  if ( zmin == zmax ) {
    zcmb = INPUTS.GENRANGE_REDSHIFT[0] ; 
  }
  else {

    // pick from pure Hubble distribution with reweight based on rate model
    zcmb = genz_hubble(zmin,zmax, GENLC.RATEPAR ) ;  
  }

  return(zcmb);

}  // end of gen_redshiftt_cmb



// *******************************************
double gen_redshift_helio(void) {
  //
  // Created Jan 2018 by R.Kessler
  // Convert zcmb to zhelio and apply peculiar velocity (vpec in km/s)
  // Function returns z_helio, and also stores GENLC.VPEC
  //
  // Jan 27 2021: fix to use exact zpec formula instead of approximation.
  // Feb 15 2021: if VEL_CMBAPEX==0, continue instead of returning zCMB

  double zCMB = GENLC.REDSHIFT_CMB ;
  double RA   = GENLC.RA;
  double DEC  = GENLC.DEC ;
  bool   USE_HOSTLIB_VPEC = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEVPEC );
  double vpec, zhelio, dzpec ;
  //  char fnam[] = "gen_redshift_helio" ;

  // ----------- BEGIN ------------

  // check (v10_31) legacy option to keep zhelio = zcmb

  if ( INPUTS.VEL_CMBAPEX < 1.0 ) 
    { zhelio = zCMB ; }  
  else 
    { zhelio = zhelio_zcmb_translator(zCMB, RA,DEC, "eq", -1); }

  // apply v_pec
  if ( USE_HOSTLIB_VPEC ) {
    // May 2020: do nothing; will be done later in HOSTLIB call.
    vpec = 0.0;  
  }
  else {
    // pick random vpec
    vpec = ((double)INPUTS.GENSIGMA_VPEC) * getRan_Gauss(2) ;    
  }


  GENLC.VPEC = vpec; 
  if ( vpec != 0.0 ) {
    dzpec = vpec/LIGHT_km ;
    zhelio = (1.0+zhelio)*(1.0+dzpec) - 1.0 ;    // Jan 27 2021 
  }

  return zhelio ;

} // end of gen_redshift_helio


// ==================================
void gen_redshift_LCLIB(void) {

  // Aug 22 2018: 
  // for LCLIB model with REDSHIFT, store redshift and MU.
  // Also check for HOSTLIB photo-z
  //

  int LDMP =  0 ;
  double granz, ZERR;
  double RA        = GENLC.RA;
  double DEC       = GENLC.DEC ;
  double ZCMB_TRUE = LCLIB_EVENT.REDSHIFT;  // see genmag_LCLIB.c
  double ZHEL_TRUE = zhelio_zcmb_translator(ZCMB_TRUE, RA,DEC, "eq", -1);
  //  char fnam[] = "gen_redshift_LCLIB" ;

  // --------------- BEGIN --------------

  if ( INDEX_GENMODEL != MODEL_LCLIB ) { return ; }

  if ( LDMP ) {
    printf("    xxx ------------------------------------- \n");
    printf(" 0. xxx ZCMB_TRUE = %.4f  (CID=%d) \n", 
	   ZCMB_TRUE, GENLC.CID ) ;
    printf(" 1. xxx  ZPHOT(Gauss .05) = %.4f += %.4f \n",
	   SNHOSTGAL.ZPHOT, SNHOSTGAL.ZPHOT_ERR ); 
  }

  GENLC.REDSHIFT_HELIO = ZHEL_TRUE ;
  GENLC.REDSHIFT_CMB   = ZCMB_TRUE ;
  if ( LCLIB_INFO.IPAR_REDSHIFT > 0  ) 
    { GEN_SNHOST_DRIVER(ZHEL_TRUE, GENLC.PEAKMJD); }
  else
    { SNHOSTGAL.ZPHOT = SNHOSTGAL.ZPHOT_ERR  = 0.0 ; }

  if ( LDMP ) {
    printf(" 2. xxx  ZPHOT(HOSTLIB) = %.4f += %.4f  \n",
	   SNHOSTGAL.ZPHOT, SNHOSTGAL.ZPHOT_ERR ); 
  }

  granz = getRan_GaussClip(1, (double)-3.0, (double)+3.0) ;
  ZERR  = INPUTS.GENSIGMA_REDSHIFT ;
  SNHOSTGAL.ZSPEC             = ZHEL_TRUE ;
  SNHOSTGAL.ZSPEC_ERR         = ZERR ;
  GENLC.TEMPLATE_INDEX        = LCLIB_EVENT.ID ; 
  GENLC.REDSHIFT_HOST         = ZHEL_TRUE ;
  GENLC.REDSHIFT_CMB_SMEAR    = ZCMB_TRUE + ZERR*granz ; 
  GENLC.REDSHIFT_HELIO_SMEAR  = ZHEL_TRUE + ZERR*granz ; 
  GENLC.REDSHIFT_SMEAR_ERR    = ZERR ;
  SNHOSTGAL.SNSEP             = 0.0 ;

  // get update DLMU and LENSDMU; these quantities are NOT
  // used for anything, so this is just to avoid crazy
  // MU-vs-ZCMB in the output.
  gen_distanceMag( GENLC.REDSHIFT_CMB,            // input
		   GENLC.REDSHIFT_HELIO,          // input 
		   &GENLC.DLMU, &GENLC.LENSDMU ); // returned 

  /* xxxxxxxxx
  printf(" xxx zSpec = %.4f +- %.4f  (granz=%.4f) \n",
	 SNHOSTGAL.ZSPEC, SNHOSTGAL.ZSPEC_ERR, granz);
	 xxxxxxx  */

  return ;

} // end   gen_redshift_LCLIB

// *******************************************
void gen_zsmear(double zerr) {

  // Created Nov 7, 2009:
  // generate smeared redshift and fill GENLC struct
  // Note that randoms are generated to stay synced.
  // Also generate smeared vpec -> GENLC.VPEC_SMEAR.
  //
  // Sep 23, 2011: float -> double
  // Apr 26 2017: 
  // + skip zerr calc if zerr>0.999 to avoid abort on ZGEN_MIN
  //    (recall zerr=1 --> REDSHIFT_FINAL[_ERR] = -9
  //
  // Oct 18 2017:  add INPUTS.GENBIAS_REDSHIFT to zsmear
  // Jan 05 2018:  compute GENLC.VPEC_SMEAR
  // Apr 20 2018:  bail for FIXMAG model where z=zerr=0
  // Oct 26 2020:  implement RESTORE_WRONG_VPEC logic

  int    i, NZRAN ;
  double zsmear, zerr_loc;
  double ZGEN_MIN = 0.0001;
  double zhelio, RA, DEC ;
  char   fnam[] = "gen_zsmear";

  // ---------- BEGIN ----------

  if ( INDEX_GENMODEL == MODEL_LCLIB ) { 
    // set all redshift info to zero
    GENLC.REDSHIFT_CMB_SMEAR    = GENLC.REDSHIFT_CMB    = 0.0 ;
    GENLC.REDSHIFT_HELIO_SMEAR  = GENLC.REDSHIFT_HELIO  = 0.0 ;
    GENLC.REDSHIFT_SMEAR_ERR    = 0 ;
    return ; 
  }


  // generate list of randoms if not already done
  if ( GENLC.REDSHIFT_RAN[0] < -1.0 ) {
    for ( i=0; i < MXZRAN; i++ ) 
      { GENLC.REDSHIFT_RAN[i] = getRan_Gauss(1); }
  }


  if ( zerr == 0.0 || zerr > 0.999 ) {
    GENLC.REDSHIFT_HELIO_SMEAR   = GENLC.REDSHIFT_HELIO ;
    GENLC.REDSHIFT_SMEAR_ERR     = zerr ;
    goto ZCMB_SMEAR ;
  }

  // don't work with ZERR that is larger than the redshift.
  if ( zerr < GENLC.REDSHIFT_CMB ) 
    { zerr_loc = zerr; }
  else
    { zerr_loc = GENLC.REDSHIFT_CMB ; }

  NZRAN= 0;

 ZSMEAR:
  if ( NZRAN >= MXZRAN ) {
    print_preAbort_banner(fnam);
    printf("\t z[cmb,hel] = %f , %f \n", 
	   GENLC.REDSHIFT_CMB, GENLC.REDSHIFT_HELIO );
    printf("\t VEL_CMBAPEX = %.2f \n", INPUTS.VEL_CMBAPEX); 
    printf("\t VPEC = %.1f  (GALID=%lld) \n", GENLC.VPEC, SNHOSTGAL.GALID );
    sprintf(c1err,"Could not generate random Z after %d tries", NZRAN);
    sprintf(c2err,"Each Z < ZGEN_MIN = %f ", ZGEN_MIN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  /*
  printf(" xxx zerr_loc=%f  RAN=%f  RA,DEC=%f,%f\n",
	 zerr_loc, GENLC.REDSHIFT_RAN[NZRAN],
	 GENLC.RA, GENLC.DEC ); fflush(stdout); // xxxxx
  */

  zsmear  = GENLC.REDSHIFT_RAN[NZRAN] * zerr_loc ; 
  zsmear += INPUTS.GENBIAS_REDSHIFT ;   // add user-defined bias
  GENLC.REDSHIFT_HELIO_SMEAR  = GENLC.REDSHIFT_HELIO + zsmear ;
  GENLC.REDSHIFT_SMEAR_ERR    = zerr ;

  // Dec 18 2015: set zspec for host
  SNHOSTGAL.ZSPEC     += zsmear ;  // smear redshift 
  SNHOSTGAL.ZSPEC_ERR  = zerr ;
  
  // avoid z-smearing to pathological value
  if ( GENLC.REDSHIFT_HELIO_SMEAR < ZGEN_MIN ) 
    { NZRAN++ ; goto ZSMEAR ; }


 ZCMB_SMEAR:
  // translate observed zhelio_smear back to zcmb_smear 
  zhelio = GENLC.REDSHIFT_HELIO_SMEAR ;
  RA     = GENLC.RA ;  
  DEC    = GENLC.DEC ;

  if ( INPUTS.VEL_CMBAPEX < 1.0 ) 
    { GENLC.REDSHIFT_CMB_SMEAR = zhelio ; } // legacy option, ZCMB = ZHELIO
  else 
    { GENLC.REDSHIFT_CMB_SMEAR = 
	zhelio_zcmb_translator(zhelio,RA,DEC, "eq", +1); 
    }

 
  // --------------------------------------------
  // Determine VPEC correction using Gaussian-random number.
  // Note that correction has oppoisite sign of true value.
  // If sim-input VPEC_ERR is >= true scatter, this is a flag
  // to NOT apply a correction and set measured VPEC_SMEAR=0.
  // If VPEC is from HOSTLIB, a correction is made.

  bool   USE_HOSTLIB_VPEC = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEVPEC );
  bool   APPLY_VPEC_SMEAR = 
    ( USE_HOSTLIB_VPEC || INPUTS.VPEC_ERR < INPUTS.GENSIGMA_VPEC );
  double GAURAN_VPEC = GENLC.REDSHIFT_RAN[MXZRAN-1];
  double VPEC_ERR, SGN_VPEC ;
  GENLC.VPEC_SMEAR   = 0.0; // default is no vpec estimate

  if ( APPLY_VPEC_SMEAR ) {
    if ( USE_HOSTLIB_VPEC )
      { VPEC_ERR = SNHOSTGAL.VPEC_ERR; } // from HOSTLIB
    else
      { VPEC_ERR = INPUTS.VPEC_ERR; }     // from sim-input file
	  
    if ( INPUTS.RESTORE_WRONG_VPEC ) 
      { SGN_VPEC = -1.0; }  // Jan 2018: incorrect sign convention
    else
      { SGN_VPEC = +1.0; }  // Oct 2020: correct sign convention

    GENLC.VPEC_SMEAR =  (SGN_VPEC*GENLC.VPEC) + (VPEC_ERR * GAURAN_VPEC) ; 
  }

  return ;

} // end of gen_zsmear


// =====================================================
void gen_distanceMag(double zCMB, double zHEL, double *MU, double *LENSDMU) {

  // Created Apr 2017
  // for input CMB redshift zCMB, returns
  //  MU = true distance modulus
  //  LENSDMU = distance correction
  //
  //  9/17/2017: 
  //    + allow new model with symmetric Gaussian
  //    + return true MU instead of corrected MU
  //
  // Jan 5 2018: pass new arg zHEL
  //

  double lensDMU, ran1 ;
  //  char fnam[] = "gen_distanceMag" ;

  // -------------- BEGIN ------------


  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID ) 
    { ran1 = 0.5; } // no randoms for GRID mode
  else
    { ran1 = getRan_Flat1(2); }  // normal gen: always burn random: May 7 2017

  if ( IGNOREFILE(INPUTS.WEAKLENS_PROBMAP_FILE) ) 
    { lensDMU = 0.0 ; }
  else 
    { lensDMU = gen_lensDMU(zCMB,ran1);  }


  if ( INPUTS.WEAKLENS_DSIGMADZ > 1.0E-8 ) {
    lensDMU = zCMB * INPUTS.WEAKLENS_DSIGMADZ * getRan_Gauss(1) ;
  }

  lensDMU *= INPUTS.WEAKLENS_DMUSCALE ; // user-scale

  // load return args
  *MU      = gen_dLmag(zCMB,zHEL);
  *LENSDMU = lensDMU ;

  return;

} // end gen_distanceMag

// ********************************
double gen_dLmag(double zCMB, double zHEL ) {

  // Returns lumi-distance mag with input cosmology
  // Note that INPUTS.H0 ~ 70 km/s/Mpc, and H0 -> ~2E-18

  double mu ;
  if ( zCMB <= 1.0E-10 ) { return 0.0 ; }  // avoid inf, May 2013
  mu = dLmag (zCMB, zHEL, &INPUTS.HzFUN_INFO);

  return(mu) ;

} // end of gen_dLmag


// *********************************
double genz_hubble ( double zmin, double zmax, RATEPAR_DEF *RATEPAR ) {


  /************
   return random SN redshift between Z0,Z1 = INPUTS.REDSHIFT[0,1]
   using dV/dz/(1+z) Hubble distribution. 
   Note that time-dilation factor 1/(1+z) for SN.

   For special tests, if INPUTS.DNDZ_ZEXP_REWGT  is non-zero, 
   then re-wgt
      dN/dz -> dN/dz * z^ZEXP_REWGT

    Example: since dN/dz ~ z^2 for z << 1,
    setting ZEXP_REWGT = -2.0 will result in a flat
    z-distribution.  Setting ZEXP_REWGT = -1 results
    in linear z-distribution.

   Global variable ZGENWGT_MAX is set on the first pass;
   this is the max weight based on looping over all z.
  
   Dec 17, 2006: fix dumb bug and include time-dilation factor

   Jan 22, 2011: fix dumb bug and force NZ to be at least 3 to
                 avoid divide-by-zero when zmax is very small

   Sep 23, 2011: float -> double 

   May 5, 2014: float rangen -> double getRan_Flat
                and remaining legacy floats -> double.

   Feb 5 2015: move 'ilist=1' after PICKRAN instead of before so that
               it works with FLATRAN option.

   Dec 21 2015: 
      pass zmin & zmax as args instead of using INPUTS.GENRANGE_REDSHIFT.
      Allows picking from narrower redshift windows.

   Jan 7 2016: if NEWZRANGE is true, then recompute ZGENWGT_MAX.
               --> more efficient generation when GENRANGE_REDSHIFT
                   is different for each simlib entry.

  Aug 12 2016; refactor to pass RATEPAR struct.

  Jan 26 2018: if USE_SIMLIB_DISTANCE, return wgt=1

  Feb 19 2018:
    set ZGENWGT_MAX=1 if 
      (USE_FLAT || USE_SIMLIB_DISTANCE || USE_SIMLIB_REDSHIFT)

 Nov 24 2019: if zmin == zmax, return immediately

  *****************/

  double z, zran, z_atmax, dz, w, wgt, wran1 ; // xxx H0, OM, OL, W0, 
  double zrange[2] = { zmin, zmax } ;
  int iz, NZ, ISFLAT, ISPOLY, ilist, NEWZRANGE, FIRST ;
  char fnam[] = "genz_hubble" ;

  // --------------- BEGIN ------------

  if ( zmin <= 1.0E-9 || zmax <= 1.0E-9 ) {
    sprintf(c1err,"Invalid zmin,zmax = %le, %le at CID=%d LIBID=%d", 
	    zmin, zmax, GENLC.CID, GENLC.SIMLIB_ID );
    sprintf(c2err,"Both must be > %le (10pc)", ZAT10PC );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( zmin == zmax ) { return(zmin); }

  ISFLAT = ( strcmp(RATEPAR->NAME,"FLAT" ) == 0 ) ;
  ISPOLY = ( strcmp(RATEPAR->NAME,"ZPOLY") == 0 ) ;

  int Iget_zFromSIMLIB = 
    (INPUTS.USE_SIMLIB_DISTANCE || INPUTS.USE_SIMLIB_REDSHIFT );
  FIRST = ( RATEPAR->ZGENWGT_MAX == 0.0 && Iget_zFromSIMLIB==0 ) ;

  // =======================================

  if ( ISFLAT ) { RATEPAR->ZGENWGT_MAX = 1.0 ; goto PICKRAN ;  }

  // get max wgt if zmin or zmax have changed.
  NEWZRANGE = ( ( zmin != RATEPAR->ZGENMIN_STORE) || 
		( zmax != RATEPAR->ZGENMAX_STORE) );

  // a + a1*z + a2*z^2 + a3*z^3
  // a1 + 2*a2*z + 3*a3 * z^2 = 0

  if ( NEWZRANGE ) {

    RATEPAR->ZGENWGT_MAX = z_atmax = -9.0 ;

    NZ = (int)(zmax * 100) ;
    if ( NZ < 3 ) { NZ = 3 ; }

    dz = (zmax - zmin) / (double)(NZ-1) ;

    for ( iz=1; iz<=NZ; iz++ ) {
      z = zmin + dz * (double)(iz-1) ;

      if ( ISPOLY ) {
	// xxx mark  wgt = eval_GENPOLY(z, &RATEPAR->MODEL_ZPOLY, fnam) ; 
	w = eval_GENPOLY(z, &RATEPAR->MODEL_ZPOLY, fnam) ; 
      }
      else {
	w = dVdz (z, &INPUTS.HzFUN_INFO);
	w /= (1.0+z);          // Dec 2006: time dilation factor
	w *= genz_wgt(z,RATEPAR) ;
      }

      if ( w > RATEPAR->ZGENWGT_MAX ) {
	RATEPAR->ZGENWGT_MAX = w ;  
	z_atmax     = z;
      }

    } // end of iz loop

    if ( RATEPAR->ZGENWGT_MAX == 0.0 ) {
      print_preAbort_banner(fnam);
      printf("\t RATE MODEL = '%s' \n", RATEPAR->NAME );
      sprintf(c1err,"Max dN/dz*wgt = 0 ?!?!?");
      sprintf(c2err,"zmin=%f zmax=%f", zmin, zmax);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( FIRST ) {
      printf("  Found Max dN/dz * wgt = %e at z = %8.3f \n", 
	     RATEPAR->ZGENWGT_MAX, z_atmax); fflush(stdout);
    }

  }  // end of ZGENWGT_MAX if-block


  // =======================================

 PICKRAN:

  ilist = 1 ;  // select which random list

  // pick random redshift
  zran   = getRan_Flat ( ilist, zrange ); // random redshift
  wran1  = getRan_Flat1( ilist );  // for weight

  if ( zran < zmin || zran > zmax ) { goto PICKRAN ; }

  // compute wgt = dN/dz for this z
  if ( ISFLAT )  { 
    // flat redshift distribution
    wgt = 1.0 ; 
  }
  else if ( ISPOLY ) {
    // user-specified polynomial function of redshift
    wgt = eval_GENPOLY(zran, &RATEPAR->MODEL_ZPOLY, fnam) ; 
    wgt /= RATEPAR->ZGENWGT_MAX ;
  }
  else if ( INPUTS.USE_SIMLIB_DISTANCE ) {
    wgt = 1.0 ;
  }
  else {
    // physical distribution
    w    = dVdz (zran, &INPUTS.HzFUN_INFO);
    w   /= (1.0+zran);  
    w   *= genz_wgt(zran,RATEPAR);
    if ( w > RATEPAR->ZGENWGT_MAX )  { RATEPAR->ZGENWGT_MAX = w ; }
    wgt  =  w / RATEPAR->ZGENWGT_MAX ;
  }

  if ( wran1 > wgt ) { goto PICKRAN ; }

  // store zmin & zmax for next time
  RATEPAR->ZGENMIN_STORE = zmin ;
  RATEPAR->ZGENMAX_STORE = zmax ;
  
  return(zran) ;

}  // end of genz_hubble


// *******************************************
double genz_wgt(double z, RATEPAR_DEF *RATEPAR ) {

  /*********
   Created Dec 17, 2007 by R.Kessler

   Compute weight factor on Hubble distribution.
   Options are 

    - rate model
    - z^zwgt

   Feb 03, 2014: ZPOLY now has z^3 term.

   Aug 30 2017: use global scale (DNDZ_ALLSCALE) for all models.
                Needed for SIMSED since these models can be either
                Ia or NON1A.
   June 14 2019:
     Apply DNDZ_SCALE[1] to SIMSED models in addition to NON1A models.

  *******/

  double w, zexp, wpoly ;
  char fnam[] = "genz_wgt" ;

  // ------ BEGIN ---------

  if ( RATEPAR->NMODEL_ZRANGE <= 0 ) { return(0.0); }

  w = 1.0 ;

  // Jul 2007: check for rate-evolution models
  w *= SNrate_model(z, RATEPAR);  

  // check for user re-wgt functions

  // first a simple exponential of z
  zexp = RATEPAR->DNDZ_ZEXP_REWGT ;
  if ( zexp != 0.0 ) w *= pow(z,zexp);

  //check polynominal re-weight
  wpoly = eval_GENPOLY(z, &RATEPAR->DNDZ_ZPOLY_REWGT, fnam) ; 
  w *= wpoly ;

  // global rate scale for all models.
  w *= RATEPAR->DNDZ_ALLSCALE ;

  // check DNDZ scale (Apr 19 2017)
  if ( INPUTS.NON1A_MODELFLAG > 0 || INDEX_GENMODEL == MODEL_SIMSED )
    { w *= RATEPAR->DNDZ_SCALE[1]; } // NON1A scale
  else
    { w *= RATEPAR->DNDZ_SCALE[0]; } // Ia scale

  return(w) ;

} // end of genz_wgt


// *******************************************
void  init_RATEPAR ( RATEPAR_DEF *RATEPAR ) {

  int i;
  char fnam[] = "init_RATEPAR" ;
  // ------------- BEGIN -------------

  RATEPAR->DNDZ_ZEXP_REWGT     = 0.0 ;
  RATEPAR->DNDZ_SCALE[0]       = 1.0 ; // Apr 19 2017
  RATEPAR->DNDZ_SCALE[1]       = 1.0 ; // Apr 19 2017
  RATEPAR->DNDZ_ALLSCALE       = 1.0 ; // Aug 30 2017
  RATEPAR->RATEMAX = 0.0 ;

  init_GENPOLY(&RATEPAR->MODEL_ZPOLY);
  init_GENPOLY(&RATEPAR->MODEL_BPOLY);

  // init REWGT to 1.000
  parse_GENPOLY("1", "DNDZ_ZPOLY_REWGT", &RATEPAR->DNDZ_ZPOLY_REWGT, fnam);

  sprintf(RATEPAR->NAME, "NONE"); 
  RATEPAR->NMODEL_ZRANGE  = 0 ;
  RATEPAR->INDEX_MODEL    = 0 ;

  for ( i=0; i <= MXRATEPAR_ZRANGE; i++ ) {
    RATEPAR->MODEL_PARLIST[i][0]  = 0.0 ;  // rate param 1
    RATEPAR->MODEL_PARLIST[i][1]  = 0.0 ;  // rate param 2
    RATEPAR->MODEL_ZRANGE[i][0] = 0.0 ;  // Zmin
    RATEPAR->MODEL_ZRANGE[i][1] = ZMAX_SNANA ;  // Zmax
  }


  RATEPAR->ZGENWGT_MAX = 0.0;
  RATEPAR->ZGENMIN_STORE =  RATEPAR->ZGENMAX_STORE = 0.0 ; 

} // end init_RATEPAR


// *******************************************
double SNcount_model(double zMIN, double zMAX, RATEPAR_DEF *RATEPAR ) {

  // Created Aug 2016
  // return expected number of SN between zMIN & zMAX,
  // with model params passed via RATEPAR.
  // Solid angle & time window are passed via global INPUTS.
  //
  // Sep 4 2016: add missing return(SNsum)

  double SNsum, dz, ztmp, vtmp, rtmp, tmp ;
  int NBZ, iz;
  char fnam[] = "SNcount_model" ;

  // ------------- BEGIN ------------

  SNsum = 0.0 ;

  NBZ = (int)( (zMAX - zMIN ) * 1000. ) ;
  if ( NBZ < 10 ) { NBZ = 10; }
  dz  = ( zMAX - zMIN ) / (double)NBZ ;

  for ( iz=1; iz <= NBZ; iz++ ) {
    ztmp   = zMIN + dz * ((double)iz - 0.5 ) ;
    vtmp   = dVdz (ztmp, &INPUTS.HzFUN_INFO);
    rtmp   = genz_wgt(ztmp,RATEPAR) ;   // rate * user-rewegt fudge
    tmp    = rtmp * vtmp / ( 1.0 + ztmp );
    SNsum += tmp ;
  }

  // tack on factors for solid angle and time-window.
  double dOmega = INPUTS.SOLID_ANGLE ;
  double delMJD = INPUTS.GENRANGE_PEAKMJD[1] - INPUTS.GENRANGE_PEAKMJD[0] ;
  double Tyear  = delMJD / 365.0 ;

  SNsum *= ( dOmega * Tyear * dz) ;

  return(SNsum);

} // end SNcount_model

// *******************************************
double SNrate_model(double z, RATEPAR_DEF *RATEPAR ) {

  /***
      Created July 2007 by R.Kessler
      
      returns calculated SN rate at redshift 'z' using model 
      where A= delayed component and B = prompt component.

      Check A+B or POWERLAW 

    Feb 20, 2009: allow negative powerlaw (B < 0)
    Dec 02, 2011: allow POWERLAW2 with z-dependent rate params.
    Feb 06, 2014: new ZPOLY model
    Aug 12, 2016: refactor to pass RATEPAR struct
    Dec     2016: add Strolger 2015 CC rate
    Mar 12  2021: fix bug for INDEX_RATEMODEL_ZPOLY 

  ***/

  double sfr, sfrint, h, H0, OM, OL, w0, wa, z1, rate ;
  double A, B, k, zmin, zmax, z2,z3,z4,z5, arg ;
  double MD14parList[8], R0;
  double zero=0.0 ;
  int iz, FOUND_iz, ISPOW ;
  char *cptr ;
  char fnam[] = "SNrate_model" ;

  // ------------- BEGIN ----------------

  // extract rate parameters from user input.
  // Make sure to extract from appropriate redshift range.

  A = B = rate = 0.0 ;
  FOUND_iz = 0;
  for ( iz = 1; iz <= RATEPAR->NMODEL_ZRANGE; iz++ ) {
    zmin = RATEPAR->MODEL_ZRANGE[iz][0];
    zmax = RATEPAR->MODEL_ZRANGE[iz][1];

    if ( z >= zmin && z < zmax ) {
      A  = RATEPAR->MODEL_PARLIST[iz][0];
      B  = RATEPAR->MODEL_PARLIST[iz][1];
      FOUND_iz += 1 ;
    }
  }

  // abort if we cannot find the z-dependent rate params.
  if ( FOUND_iz != 1 ) {
    sprintf(c1err,"Found %d sets of rate params for z=%f", FOUND_iz, z);
    sprintf(c2err,"Check DNDZ key(s) in sim-input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  ISPOW = 0 ;
  if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_POWERLAW  ) { ISPOW=1; }
  if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_POWERLAW2 ) { ISPOW=1; }

  // ----------------------------------
  cptr = RATEPAR->NAME ;
  z1 = 1. + z;  
  h  = INPUTS.H0/100.0 ;

  if ( strcmp(cptr,"HUBBLE") == 0 ) {
    rate = 3.0E-5 ;  // approx SNIa rate at z ~ 0
  }
  else if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_FLAT ) {
    rate = 1.0E-8 ; // just to avoid abort. Rate here is meaningless.
  }
  else if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_CCS15 ) {
    // define MD14 params corresponding to green line in
    // Fig 6 of Strolger 15 --> best fit to measured R_CC(z)
    // from A,B,C,D parameters for Eq 9.
    MD14parList[0] = 0.015;  // A
    MD14parList[1] = 1.50;   // B
    MD14parList[2] = 5.00;   // C
    MD14parList[3] = 6.10;   // D
    k=0.0091 ;
    rate  = k * (h*h) * SFRfun_MD14(z,MD14parList);
    rate *= RATEPAR->MODEL_PARLIST[1][0] ; // user-defined scale
  }
  else if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_MD14 ) {
    MD14parList[0] = 0.015;  // A
    MD14parList[1] = 2.90;   // B
    MD14parList[2] = 2.70;   // C
    MD14parList[3] = 5.60;   // D
    R0   = RATEPAR->MODEL_PARLIST[1][0] ; // user-define rate at z=0
    rate = R0 * SFRfun_MD14(z,MD14parList)/SFRfun_MD14(zero,MD14parList);

  }
  else if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_PISN ) {

    // R= -0.0508z^5 + 0.8312z^4 -4.42z^3 + 6.55z^2 + 6.38z + 1.98 (/yr/Gpc^3)
    z2=z*z; z3=z*z2; z4=z*z3; z5=z*z4; 
    rate = 1.98 + 6.38*z + 6.558*z2 - 4.42*z3 + 0.8312*z4 - 0.0508*z5;
    rate /= 1.0E9;  // convert Gpc^3 to Mpc^3
  }
  else if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_TDE ) {
    R0   = RATEPAR->MODEL_PARLIST[1][0] ; // user-define rate at z=0
    arg  = (-0.5*z/0.6) ;
    rate = R0 * pow(10.0,arg);
  }
  else if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_AB ) {
    if ( A > 1.0E-13 || B > 0.01 ) {
      sprintf(c1err,"Invalid parameters for A+B rate model");
      sprintf(c2err,"A=%e  B=%e", A, B );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    // get star formation rate 
    sfr    = SFRfun_BG03(INPUTS.H0,z);

    // determine integrated stellar mass 
    sfrint = SFR_integral( z, &INPUTS.HzFUN_INFO);

    // compute rate from 2-component model
    rate = A * sfrint  +  B * sfr ;

  }
  else if ( ISPOW ) {
    if ( A > 10.0 || B > 10. ) {
      sprintf(c1err,"Invalid parameters for powerlaw rate model");
      sprintf(c2err,"alpha=%e  beta=%e", A, B );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    rate = A * pow(z1,B);
  }
  else if ( RATEPAR->INDEX_MODEL == INDEX_RATEMODEL_ZPOLY ) {
    rate = eval_GENPOLY(z, &RATEPAR->MODEL_ZPOLY, fnam) ; 
  }
  else {
    sprintf(c1err,"Invalid model: '%s'", cptr);
    sprintf(c2err,"ISPOW = %d, INDEX_MODEL=%d",
	    ISPOW, RATEPAR->INDEX_MODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // -------------

  return(rate) ;

}  // end of SNrate_model()

// *************************************
double GALrate_model(double l, double b, RATEPAR_DEF *RATEPAR ) {

  // Created Nov 26 2017
  // Return Galactic rate for input b and l
  //
  // May 25 2018: Fix order of l,b arguments.
  // Sep 30 2020: switch to using polyEval or eval_GENPOLY

  double Rate=0.0, b_val ;
  double BPOW[MXPOLY_GALRATE+1], COSBPOW[MXPOLY_GALRATE+1], Rtest=0.0 ;
  int i;
  char fnam[] = "GALrate_model" ;

  // -------------- BEGIN ---------------

  if ( strcmp(RATEPAR->NAME,"COSBPOLY") == 0 ) {
    b_val = cos(b*RADIAN);
  }
  else if ( strcmp(RATEPAR->NAME,"BPOLY") == 0 ) {
    b_val     = fabs(b);
  }

  else {
    sprintf(c1err,"Unknown GALrate model: '%s' ", RATEPAR->NAME );
    sprintf(c2err,"Check GENMODEL key");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // - - - - - - - - - - 
  // evaluate polynomial functon
  Rate = eval_GENPOLY(b_val, &RATEPAR->MODEL_BPOLY, fnam) ;

  return(Rate);

} // end GALrate_model


// ***********************************
double gen_AV(void) {

  // Mar 20, 2020: refactor by D.Brout and R.Kessler
  // 
  // Select AV from exponential + halfGauss distribution,
  //    dN/dAv = exp(-av/tau) + exp(-0.5*av^2/sig^2)
  //      or
  // select EBV_HOST from same distribiution and then
  // AV = EVB_HOST*RV
  //
  // 

  double RV       =  GENLC.RV ;
  double AV       = -9.0;
  double EBV_HOST =  0.0;
  double epsilon  =  1.0E-12 ;
  GENGAUSS_ASYM_DEF  GENGAUSS_NULL ;
  char  fnam[] = "gen_AV" ;

  // ------------ BEGIN -------------

  // preserve old option to generate WV07 extinction model (RK)
  if ( INPUTS.WV07_GENAV_FLAG )  
    { AV = GENAV_WV07(); goto DONE ; }

  if ( INPUTS.GENPROFILE_AV.USE ) {
    copy_GEN_EXP_HALFGAUSS(&INPUTS.GENPROFILE_AV,&GENLC.GENPROFILE_AV);

    GENLC.GENPROFILE_AV.EXP_TAU = INPUTS.GENEXPTAU_AV 
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENEXPTAU_AV");// legacy
    GENLC.GENPROFILE_AV.EXP_TAU = INPUTS.GENPROFILE_AV.EXP_TAU
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENTAU_AV");

    GENLC.GENPROFILE_AV.SIGMA = INPUTS.GENGAUSIG_AV
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENEXPSIG_AV");// legacy
    GENLC.GENPROFILE_AV.SIGMA = INPUTS.GENPROFILE_AV.SIGMA
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENSIG_AV");

    GENLC.GENPROFILE_AV.RATIO = INPUTS.GENPROFILE_AV.RATIO
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENRATIO_AV0" ); 

    AV = getRan_GEN_EXP_HALFGAUSS(&GENLC.GENPROFILE_AV);
    
  }
 
  if ( INPUTS.GENPROFILE_EBV_HOST.USE ) {
    copy_GEN_EXP_HALFGAUSS(&INPUTS.GENPROFILE_EBV_HOST,
			   &GENLC.GENPROFILE_EBV_HOST);

    GENLC.GENPROFILE_EBV_HOST.EXP_TAU = INPUTS.GENPROFILE_EBV_HOST.EXP_TAU
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENTAU_EBV_HOST");

    GENLC.GENPROFILE_EBV_HOST.SIGMA = INPUTS.GENPROFILE_EBV_HOST.SIGMA
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENSIG_EBV_HOST");

    GENLC.GENPROFILE_EBV_HOST.RATIO = INPUTS.GENPROFILE_EBV_HOST.RATIO
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENRATIO_EBV0_HOST" );

    EBV_HOST = getRan_GEN_EXP_HALFGAUSS(&GENLC.GENPROFILE_EBV_HOST);
    AV       = EBV_HOST * RV ;
    
  }

  // Jun 2020: check for map in GENPDF_FILE
  if ( INPUTS.DOGEN_AV == 2 ) {   // flag to get AV from map
    GENGAUSS_NULL.USE = false ;
    AV = get_random_genPDF(PARNAME_AV, &GENGAUSS_NULL); 
  }
  if ( INPUTS.DOGEN_AV == 4 ) {  // flag to get EBV from map
    GENGAUSS_NULL.USE = false ;
    EBV_HOST = get_random_genPDF(PARNAME_EBV, &GENGAUSS_NULL); 
    AV       = EBV_HOST * RV ;
  }

 DONE: 
  if ( AV < -epsilon  || RV < epsilon ) {
    sprintf(c1err,"Crazy dust params: AV = %f and RV = %f ; ", AV, RV);
    sprintf(c2err,"Check dust params in sim-input and/or GENPDF_FILE.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }


  return(AV) ;

}  // end of gen_AV


// ***********************************
double GENAV_WV07(void) {

  // Created Sep 12,2007 by R.Kessler
  // return AV from distribution used by ESSENCE-WV07
  // Mar 2013: float -> double
  // Apr 2018: 
  //   + expf -> exp
  //   + check REWGT_EXPAV option

  double AV ;
  double tau = 0.4, sqsigma=0.01;
  double REWGT_AEXP = INPUTS.WV07_REWGT_EXPAV ;
  double AEXP, BEXP, arg_A, arg_B, W0, W;

  // ----------- BEGIN -----------

  AEXP = 1./tau;
  BEXP = 1./sqrt(sqsigma * 2. * 3.14159) ;

  if ( REWGT_AEXP > -1.0E-9 ) { AEXP *= REWGT_AEXP; } // April 2018

  W0 = AEXP + BEXP ; // weight at AV=0

  // pick random AV on defined interval

 PICKAV:
  AV = getRan_Flat ( 1 ,INPUTS.GENRANGE_AV );

  // compute relative wgt
  arg_A = -AV/tau ;                  // broad exponential
  arg_B = -0.5 * AV*AV / sqsigma ;   // sharp half-Gaussian core

  W = AEXP * exp(arg_A) + BEXP * exp(arg_B)  ;
  W /= W0;

  if ( W < getRan_Flat1(1) ) { goto PICKAV ; }

  return(AV) ;

}  // end of GENAV_WV07

// ***********************************
double gen_RV(void) {

  // Aug 2006: select RV from bifurcated gaussian
  // Mar 2008: after 1000 tries, abort
  // Apr 20, 2017: GENMEAN_RV -> GENPEAK_RV
  // Jun 12, 2020: refactor to prepare for using genPDF map file.

  double RV ;
  double  ZCMB  = GENLC.REDSHIFT_CMB ; // for z-dependent populations
  GENGAUSS_ASYM_DEF GENGAUSS_ZVAR ;
  char fnam[] = "gen_RV" ;

  // ------------ BEGIN -------------
  
  GENGAUSS_ZVAR = 
    get_zvariation_GENGAUSS(ZCMB, PARNAME_RV, &INPUTS.GENGAUSS_RV); 

  RV = get_random_genPDF(PARNAME_RV, &GENGAUSS_ZVAR); 

  return RV ;

}  // end of gen_RV


// ************************************
void init_RANDOMsource(void) {

  char fnam[] = "init_RANDOMsource" ;

  // ------------- BEGIN ----------

  GENLC.CIDOFF = INPUTS.CIDOFF ; 
  sprintf(BANNER," %s : CIDOFF=%d ",  fnam, GENLC.CIDOFF);
  print_banner ( BANNER );

  // initialize simlib/cadence library
  SIMLIB_INIT_DRIVER();

}  // end of init_RANDOMsource

// *************************************************
void SIMLIB_INIT_DRIVER(void) {

  // Jan 31 2021: 
  //   for INIT_ONLY flag, return after initGlobalHeader in case
  //   the global header has rate info such as SOLID_ANGLE.
  
  //  char fnam[] = "SIMLIB_INIT_DRIVER" ;

  // --------------- BEGIN --------------

  // xxx  if ( INPUTS.INIT_ONLY == 1 ) { return; } // July 30 2020

  SIMLIB_initGlobalHeader();        // generic init of SIMLIB_GLOBAL_HEADER

  SIMLIB_readGlobalHeader_TEXT();   // open and read global header

  SIMLIB_prepGlobalHeader();   

  if ( INPUTS.INIT_ONLY == 1 ) { return; } 

  SIMLIB_findStart();    // find first LIBID to start reading

} // end SIMLIB_INIT


// *************************************************
void SIMLIB_initGlobalHeader(void) {

  // Created Aug 2017
  // initialize global SIMLIB header arrays.

  SIMLIB_GLOBAL_HEADER.SURVEY_NAME[0]     = 0 ;
  SIMLIB_GLOBAL_HEADER.SUBSURVEY_NAME[0]  = 0 ;
  SIMLIB_GLOBAL_HEADER.FILTERS[0]         = 0 ;
  SIMLIB_GLOBAL_HEADER.USERNAME[0]        = 0 ;
  SIMLIB_GLOBAL_HEADER.PIXSIZE            = 0.0 ;
  SIMLIB_GLOBAL_HEADER.SOLID_ANGLE        = 0.0 ;
  SIMLIB_GLOBAL_HEADER.NGENSKIP_PEAKMJD   = 0 ;
  sprintf(SIMLIB_GLOBAL_HEADER.TELESCOPE, "UNKNOWN"); // cannot be blank string

  sprintf(SIMLIB_GLOBAL_HEADER.SKYSIG_UNIT, "%s", 
	  SIMLIB_SKYSIG_SQPIX );
  sprintf(SIMLIB_GLOBAL_HEADER.PSF_UNIT,    "%s", 
	  SIMLIB_PSF_PIXEL_SIGMA );  // default
  SIMLIB_GLOBAL_HEADER.NEA_PSF_UNIT = false;

  SIMLIB_FLUXERR_COR.USE     = 0 ;
  SIMLIB_TEMPLATE.USEFLAG    = 0 ;
  SIMLIB_TEMPLATE.NFIELD_OVP = 0;
  SIMLIB_TEMPLATE.TEXPOSE_SPECTROGRAPH = 0.0 ;

  int i;
  for ( i=0; i < MXFILTINDX; i++ ) { 

    // set global template params to zero (first field index=0)
    SIMLIB_TEMPLATE.SKYSIG[0][i]    = 0.0 ;
    SIMLIB_TEMPLATE.READNOISE[0][i] = 0.0 ;
    SIMLIB_TEMPLATE.ZPT[0][i]       = 0.0 ;  
  }

  SIMLIB_GLOBAL_HEADER.FLUXERR_ADD_FILTERS[0] = 0 ;
  SIMLIB_GLOBAL_HEADER.NFLUXERR_COR           = 0 ;

  SIMLIB_GLOBAL_HEADER.NPE_PIXEL_SATURATE = INPUTS.NPE_PIXEL_SATURATE ;
  SIMLIB_GLOBAL_HEADER.PHOTFLAG_SATURATE  = INPUTS.PHOTFLAG_SATURATE ;
  SIMLIB_GLOBAL_HEADER.PHOTFLAG_SNRMAX    = INPUTS.PHOTFLAG_SNRMAX ;
  SIMLIB_GLOBAL_HEADER.PHOTFLAG_NEARPEAK  = INPUTS.PHOTFLAG_NEARPEAK ;

  return ;

} // end SIMLIB_initGlobalHeader


// *************************************************
void SIMLIB_readGlobalHeader_TEXT(void) {

  // Re-factored Aug 2017
  // Open SIMLIB file and read global header into
  // SIMLIB_GLOBAL_HEADER structure.
  //
  // Sep 3 2020: check REQUIRE_DOCANA

  char PATH_DEFAULT[2*MXPATHLEN];
  char *OPENFILE      = INPUTS.SIMLIB_OPENFILE;
  int  REQUIRE_DOCANA = INPUTS.REQUIRE_DOCANA ; 
  int  OPENMASK       = OPENMASK_VERBOSE ;
  if (REQUIRE_DOCANA) { OPENMASK += OPENMASK_REQUIRE_DOCANA; }
  char c_get[80];
  int  NTMP, NFILT;
  char fnam[] = "SIMLIB_readGlobalHeader_TEXT" ;

  // ---------- BEGIN ----------

  print_banner(fnam);

  sprintf(PATH_DEFAULT, "%s %s/simlib",  PATH_USER_INPUT, PATH_SNDATA_ROOT );
  fp_SIMLIB = snana_openTextFile(OPENMASK,PATH_DEFAULT, INPUTS.SIMLIB_FILE, 
				 OPENFILE, &INPUTS.SIMLIB_GZIPFLAG );
  
  if ( fp_SIMLIB == NULL ) {
    abort_openTextFile("SIMLIB_FILE", PATH_DEFAULT, INPUTS.SIMLIB_FILE, fnam);
  }

  // - - - - - - - - - - - - - - - -
  // read header keywords. Stop when we reach "BEGIN"

  while( (fscanf(fp_SIMLIB, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"SURVEY:") == 0 ) {
      readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.SURVEY_NAME );

      // 9.15.2021: check optoin to override survey name,
      // E.g., lowz and hiz from LSST can be split into LOWZ and LSST.
      char *S = INPUTS.SIMLIB_SURVEY ;
      if ( !IGNOREFILE(S) ) 
	{ sprintf(SIMLIB_GLOBAL_HEADER.SURVEY_NAME,"%s", S); }
      
    }
    else if ( strcmp(c_get,"SUBSURVEY_LIST:")==0  ) {
      readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.SUBSURVEY_NAME );
    }
    else if ( strcmp(c_get,"TELESCOPE:") == 0 ) {
      readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.TELESCOPE );
    }
    else if ( strcmp(c_get,"PIXSIZE:") == 0 ) {      
      readdouble(fp_SIMLIB, 1, &SIMLIB_GLOBAL_HEADER.PIXSIZE );
    }
    else if ( strcmp(c_get,"SOLID_ANGLE:") == 0 ) { 
      readdouble(fp_SIMLIB, 1, &SIMLIB_GLOBAL_HEADER.SOLID_ANGLE );
    }
    else if ( strcmp(c_get,"NLIBID:") == 0 ) { 
      readint(fp_SIMLIB, 1, &SIMLIB_GLOBAL_HEADER.NLIBID );
      SIMLIB_GLOBAL_HEADER.NLIBID_VALID = SIMLIB_GLOBAL_HEADER.NLIBID;
    }
    else if ( strcmp(c_get,"PSF_UNIT:") == 0 ) {
      readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.PSF_UNIT );
    }
    else if ( strcmp(c_get,"SKYSIG_UNIT:") == 0 ) {
      readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.SKYSIG_UNIT );
    }
    else if ( strcmp(c_get,"FILTERS:") == 0 ) {
      readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.FILTERS );
    }
    else if ( strcmp(c_get,"USER:") == 0 ) {
      readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.USERNAME );
    }
    else if ( strcmp(c_get,"NPE_PIXEL_SATURATE:") == 0 ) {
      if ( SIMLIB_GLOBAL_HEADER.NPE_PIXEL_SATURATE == 1000000000 ) 
	{ readint(fp_SIMLIB, 1, &SIMLIB_GLOBAL_HEADER.NPE_PIXEL_SATURATE); }
    }    
    else if ( strcmp(c_get,"PHOTFLAG_SATURATE:") == 0 || 
	      strcmp(c_get,"PHOTMASK_SATURATE:") == 0 ) {
      if ( SIMLIB_GLOBAL_HEADER.PHOTFLAG_SATURATE == 0 ) 
	{ readint(fp_SIMLIB, 1, &SIMLIB_GLOBAL_HEADER.PHOTFLAG_SATURATE ); }
    }
    else if ( strcmp(c_get,"PHOTFLAG_SNRMAX:") == 0 || 
	      strcmp(c_get,"PHOTMASK_SNRMAX:") == 0 ) {
      if ( SIMLIB_GLOBAL_HEADER.PHOTFLAG_SNRMAX == 0 ) 
	{ readint(fp_SIMLIB, 1, &SIMLIB_GLOBAL_HEADER.PHOTFLAG_SNRMAX ); }
    }
    else if ( strcmp(c_get,"PHOTFLAG_NEARPEAK:") == 0 || 
	      strcmp(c_get,"PHOTMASK_NEARPEAK:") == 0 ) {
      if ( SIMLIB_GLOBAL_HEADER.PHOTFLAG_NEARPEAK==0 ) 
	{ readint(fp_SIMLIB, 1, &SIMLIB_GLOBAL_HEADER.PHOTFLAG_NEARPEAK ); }
    }

    else if ( strcmp(c_get, "FLUXERR_ADD:" ) == 0 ) { 
      if ( INPUTS.OPT_FUDGE_SNRMAX == 0 ) {
	readchar(fp_SIMLIB,SIMLIB_GLOBAL_HEADER.FLUXERR_ADD_FILTERS);
	NFILT = strlen(SIMLIB_GLOBAL_HEADER.FLUXERR_ADD_FILTERS);
	readdouble(fp_SIMLIB, NFILT, SIMLIB_GLOBAL_HEADER.FLUXERR_ADD_VALUES);
      }
    }

    else if ( strcmp(c_get,"FLUXERR_COR:") == 0 ) { 
      NTMP = SIMLIB_GLOBAL_HEADER.NFLUXERR_COR ;
      if ( NTMP < MXFLUXERR_COR_SIMLIB ) {
	readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.FLUXERR_COR_FILTERS[NTMP]);
	NFILT = strlen(SIMLIB_GLOBAL_HEADER.FLUXERR_COR_FILTERS[NTMP]);
	readfloat(fp_SIMLIB, 1,
		  &SIMLIB_GLOBAL_HEADER.FLUXERR_COR_LOG10SNR[NTMP]);
	readfloat(fp_SIMLIB, NFILT,
		  SIMLIB_GLOBAL_HEADER.FLUXERR_COR_SCALE[NTMP]);
      }
      SIMLIB_GLOBAL_HEADER.NFLUXERR_COR++ ;
    } 


    else if ( strcmp(c_get,"GENSKIP_PEAKMJD:") == 0 ) {
      NTMP  = SIMLIB_GLOBAL_HEADER.NGENSKIP_PEAKMJD ; 
      readdouble(fp_SIMLIB, 2, SIMLIB_GLOBAL_HEADER.GENSKIP_PEAKMJD[NTMP] );
      SIMLIB_GLOBAL_HEADER.NGENSKIP_PEAKMJD++ ;
    }

    else if ( strcmp(c_get,"BEGIN") == 0 ) 
      { return ; }  // DONE READING GLOBAL HEADER ==> BAIL

    else if ( strcmp(c_get,"LIBID:")   == 0 ) {
      sprintf(c1err,"Found 1st LIBID before BEGIN keyword.");
      sprintf(c2err,"Check simlib file: %s \n", OPENFILE );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  } // end while

  return ;

} // end SIMLIB_readGlobalHeader_TEXT



// *************************************************
void SIMLIB_prepGlobalHeader(void) {

  // Created Aug 2017
  // Transfer SIMLIB_GLOBAL_HEADER to GENLC struct, 
  // and print stuff to stdout.
  // This is called after reading global header from 
  // whatever format.
  //
  // Mar 7 2018: replace filtindx_ with INTFILTER
  // Dec 1 2020: abort if GENLC.IDSURVEY < 0

  int i, NTMP, ifilt, ifilt_obs ;
  char cfilt[4], *FILTERS, *TEL ;
  char *SURVEY, *SUBSURVEY ;
  char fnam[] = "SIMLIB_prepGlobalHeader" ;

  // -------------- BEGIN ------------

  print_banner(fnam);


  SURVEY    = SIMLIB_GLOBAL_HEADER.SURVEY_NAME ;
  SUBSURVEY = SIMLIB_GLOBAL_HEADER.SUBSURVEY_NAME ;
  
  sprintf(GENLC.SURVEY_NAME,    "%s", SURVEY);
  sprintf(GENLC.SUBSURVEY_NAME, "%s", SUBSURVEY);
  if ( !IGNOREFILE(SUBSURVEY) ) 
    { sprintf(GENLC.SUBSURVEY_NAME,"%s",SURVEY); }
  printf("\t SIMLIB Survey    : %s \n", SURVEY );

  // get integer IDSURVEY from SURVEY string
  read_SURVEYDEF();   
  GENLC.IDSURVEY = get_IDSURVEY(GENLC.SURVEY_NAME);
  if ( GENLC.IDSURVEY < 0 ) {
    sprintf(c1err,"Invalid 'SURVEY: %s' in SIMLIB header", GENLC.SURVEY_NAME);
    sprintf(c2err,"Check valid SURVEY names in $SNDATA_ROOT/SURVEY.DEF" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  TEL = SIMLIB_GLOBAL_HEADER.TELESCOPE ;
  sprintf(GENLC.TELESCOPE[0],       "%s", TEL );
  printf("\t SIMLIB telescope : %s \n", TEL );


  double PIXSIZE = SIMLIB_GLOBAL_HEADER.PIXSIZE ;
  SIMLIB_OBS_GEN.PIXSIZE[0] = PIXSIZE ;
  if ( PIXSIZE > 0.000001 ) 
    {  printf("\t SIMLIB pixel size: %5.3f asec \n", PIXSIZE );  }

  int NPE_SAT = SIMLIB_GLOBAL_HEADER.NPE_PIXEL_SATURATE;
  if ( NPE_SAT < 999999999  ) 
    { printf("\t SIMLIB Pixel Saturation: %d photoelectrons. \n", NPE_SAT); }

  fflush(stdout) ;

  // - - -  prepare filter maps - - - - 
  for ( i=0; i < MXFILTINDX; i++ ) {
    GENLC.IFILTMAP_SIMLIB[i]    = -9 ;
    GENLC.IFILTINVMAP_SIMLIB[i] = -9 ;
  }
  FILTERS = SIMLIB_GLOBAL_HEADER.FILTERS ;
  NTMP = strlen(FILTERS);    GENLC.NFILTDEF_SIMLIB = NTMP;

  for ( i=0; i < NTMP; i++ ) {
    sprintf(cfilt, "%c", FILTERS[i] );
    ifilt_obs = INTFILTER(cfilt);
    GENLC.IFILTMAP_SIMLIB[i]            = ifilt_obs ;
    GENLC.IFILTINVMAP_SIMLIB[ifilt_obs] = i;
  }
  printf("\t SIMLIB Filters   : %s \n", FILTERS );

  // make sure GENFILTERS is a subset of SURVEY filters
  GENFILTERS_CHECK();
  fflush(stdout);

  char *USERNAME = SIMLIB_GLOBAL_HEADER.USERNAME ;
  if ( strlen(USERNAME) > 0 ) 
    {  printf("\t SIMLIB created by: %s \n", USERNAME );  }


  // Use SOLID_ANGLE from simlib header only if
  // INPUTS.SOLID_ANGLE is still zero ... 
  // allow user-input to over-ride SIMLIB global header
  double OMEGA = SIMLIB_GLOBAL_HEADER.SOLID_ANGLE ;
  if ( OMEGA > 1.0E-12 &&  INPUTS.SOLID_ANGLE == 0.0 ) {
    INPUTS.SOLID_ANGLE = (float)OMEGA ;
    printf("\t SIMLIB Solid Angle  : %f steridians \n", OMEGA);
  }


  NTMP  = SIMLIB_GLOBAL_HEADER.NGENSKIP_PEAKMJD ;
  for(i=0; i < NTMP; i++ ) {
    printf("\t SIMLIB MJD-SKIP RANGE: %6.0f - %6.0f \n"
	   ,SIMLIB_GLOBAL_HEADER.GENSKIP_PEAKMJD[NTMP][0]
	   ,SIMLIB_GLOBAL_HEADER.GENSKIP_PEAKMJD[NTMP][1] );    
  }

  fflush(stdout);


  char *unit ;
  unit = SIMLIB_GLOBAL_HEADER.PSF_UNIT ;
  if ( strcmp(unit,SIMLIB_PSF_PIXEL_SIGMA) == 0  ) {  }
  else if ( strcmp(unit,SIMLIB_PSF_ARCSEC_FWHM ) == 0 ) {  }
  else if ( strcmp(unit,SIMLIB_PSF_NEA_PIXEL  )  == 0 ) 
    { SIMLIB_GLOBAL_HEADER.NEA_PSF_UNIT = true; }
  else if ( strcmp(unit,SIMLIB_PSF_NEA_ARCSECSQ ) == 0 ) 
    { SIMLIB_GLOBAL_HEADER.NEA_PSF_UNIT = true; }
  else {
    sprintf(c1err,"Invalid 'PSF_UNIT: %s' in header of ", unit);
    sprintf(c2err,"%s", INPUTS.SIMLIB_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  unit = SIMLIB_GLOBAL_HEADER.SKYSIG_UNIT ;
  if ( strcmp(unit,SIMLIB_SKYSIG_SQPIX) == 0 ) {  }
  else if ( strcmp(unit,SIMLIB_SKYSIG_SQASEC ) == 0 ) {  }
  else {
    sprintf(c1err,"Invalid 'SKYSIG_UNIT: %s' in header of ", unit);
    sprintf(c2err,"%s", INPUTS.SIMLIB_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  printf("\t SIMLIB PSF unit:    %s \n", SIMLIB_GLOBAL_HEADER.PSF_UNIT );
  printf("\t SIMLIB SKYSIG unit: %s \n", SIMLIB_GLOBAL_HEADER.SKYSIG_UNIT );
  fflush(stdout);


  // abort if mandatory arg is missing

  sprintf(c2err,"See top of %s", INPUTS.SIMLIB_FILE );

  if ( strlen(SIMLIB_GLOBAL_HEADER.SURVEY_NAME) == 0 ) {
    sprintf(c1err," SURVEY not specified in SIMLIB header.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( strlen(SIMLIB_GLOBAL_HEADER.FILTERS) == 0 ) {
    sprintf(c1err,"'FILTERS  not specified in SIMLIB header.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  GENLC.SIMLIB_IDMAX = 0 ;
  SIMLIB_HEADER.NWRAP = 0 ;

  // ---------------------- ERROR FUDGES -------------------
  // check option to fudge the flux-errors (Oct 2015)
  init_fluxErrModel_legacy(INPUTS.SIMLIB_FILE, INPUTS.FUDGEOPT_FLUXERR);

  // error-fudge maps

  // add fixed FLUXERR in quadrature with each filter.
  FILTERS = SIMLIB_GLOBAL_HEADER.FLUXERR_ADD_FILTERS ;
  NTMP    = strlen(FILTERS);
  for ( ifilt=0; ifilt < NTMP; ifilt++ ) {
    sprintf(cfilt,"%c", FILTERS[ifilt] );   
    ifilt_obs = INTFILTER(cfilt);
    GENLC.SIMLIB_FLUXERR_ADDPAR[ifilt_obs] = 
      (float)SIMLIB_GLOBAL_HEADER.FLUXERR_ADD_VALUES[ifilt] ;
  }


  // scale FLUXERR as a function of log10(SNR) and filter
  // [code moved/refactored from SIMLIB_read_fluxerrCor]
  SIMLIB_prep_fluxerrScale_LEGACY();


  return ;

} // end SIMLIB_prepGlobalHeader


// *************************************************
void  SIMLIB_prep_fluxerrScale_LEGACY(void) {

  // Created Aug 2017
  // prepare fluxerr scale fudges specified in global header
  // of SIMLIB file. Transfer sparse contents of
  // SIMLIB_GLOBAL_HEADER into SIMLIB_FLUXERR_COR struct.
  //
  // Feb 2018: 
  //  + if using FLUXERRMODEL_FILE, ignore errFudge map inside SIMLIB
  // 

  int  NFLUXERR_COR = SIMLIB_GLOBAL_HEADER.NFLUXERR_COR ;
  int  NFTMP, size, ifilt, ifilt_obs, IFILTLIST_MAP[MXFILTINDX] ;

  char cfilt[2], *FILTERS;
  char fnam[] = "SIMLIB_prep_fluxerrScale_LEGACY" ;

  // --------------- BEGIN -------------


  if ( IGNOREFILE(INPUTS.FLUXERRMODEL_FILE) == 0 ) 
    { INPUTS.SIMLIB_MSKOPT |= SIMLIB_MSKOPT_IGNORE_FLUXERR_COR ; }


  if ( (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_IGNORE_FLUXERR_COR)>0 ) {
    NFLUXERR_COR = 0;
    printf("\t SIMLIB IGNORE FLUXERR_COR map. \n");
    fflush(stdout);
  } 
  
  if ( NFLUXERR_COR == 0 ) { return ; }

  if ( NFLUXERR_COR >= MXFLUXERR_COR_SIMLIB ) {
    sprintf(c1err,"NFLUXERR_COR=%d exceeds bound MXFLUXERR_COR=%d",
	    NFLUXERR_COR, MXFLUXERR_COR_SIMLIB );
    sprintf(c2err,"Check FLUXERR_COR keys in file %s", 
	    INPUTS.SIMLIB_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }



  SIMLIB_FLUXERR_COR.USE = 1 ;
  printf("\t SIMLIB Found FLUXERR_COR map. \n");
  fflush(stdout);

  // malloc memory
  NFTMP = GENLC.NFILTDEF_SIMLIB ;

  size = sizeof(int);
  SIMLIB_FLUXERR_COR.MAPSIZE = (int*)malloc( size*(NFTMP+1) );

  size = sizeof(double*);
  SIMLIB_FLUXERR_COR.LOG10SNR = (double**)malloc( size*(NFTMP+1) );
  SIMLIB_FLUXERR_COR.SCALE    = (double**)malloc( size*(NFTMP+1) );

  size = sizeof(double);
  for ( ifilt=0; ifilt <= NFTMP; ifilt++ ) {

      SIMLIB_FLUXERR_COR.LOG10SNR[ifilt] = 
	(double*)malloc(size*MXFLUXERR_COR_SIMLIB) ;

      SIMLIB_FLUXERR_COR.SCALE[ifilt]   = 
	(double*)malloc(size*MXFLUXERR_COR_SIMLIB) ;

      SIMLIB_FLUXERR_COR.MAPSIZE[ifilt]  = 0 ; 
  }


  int icor, i, M ;
  double SCALE, LOG10SNR ;

  for(icor=0; icor < NFLUXERR_COR; icor++ ) {

    FILTERS = SIMLIB_GLOBAL_HEADER.FLUXERR_COR_FILTERS[icor];
    NFTMP = PARSE_FILTLIST( FILTERS, IFILTLIST_MAP );

    // note that 'i' is sparse filter index for 
    // SIMLIB_GLOBAL_HEADER.FLUXERR_COT_FILTERS, while ifilt
    // is sparse index for SURVEY filters.

    for ( i = 0; i < NFTMP; i++ ) {
      ifilt_obs = IFILTLIST_MAP[i] ; 
      ifilt     = GENLC.IFILTINVMAP_SIMLIB[ifilt_obs] ;
      sprintf(cfilt,  "%c", FILTERSTRING[ifilt_obs] ) ;

      LOG10SNR = SIMLIB_GLOBAL_HEADER.FLUXERR_COR_LOG10SNR[icor];
      SCALE    = (double)SIMLIB_GLOBAL_HEADER.FLUXERR_COR_SCALE[icor][i];

      // abort on crazy SCALE value
      if ( SCALE < 0.1 || SCALE > 10.0 ) {
	sprintf(c1err,"Crazy ERR-SCALE(%s) = %f for  LOG10(SNR)=%f",
		cfilt, SCALE, LOG10SNR);
	sprintf(c2err,"Check FLUXERR_COR keys in file %s", 
		INPUTS.SIMLIB_FILE);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
      }

      M = SIMLIB_FLUXERR_COR.MAPSIZE[ifilt] ;
      SIMLIB_FLUXERR_COR.LOG10SNR[ifilt][M] = LOG10SNR ;
      SIMLIB_FLUXERR_COR.SCALE[ifilt][M]    = SCALE ;

      SIMLIB_FLUXERR_COR.MAPSIZE[ifilt]++ ;

    } // end i loop over sparse filter list
  } // end icor loop

  
} // end SIMLIB_prep_fluxerrScale_LEGACY



// *************************************
void SIMLIB_findStart(void) {

  // Created Aug 2017
  //   [some code moved out of SIMLIB_open]
  //
  // Move to start of first LIBID. Possible options are
  //  + first LIBID in file (do nothing here)
  //  + move to IDLOCK    
  //  + move to IDSTART
  //  + use JOBID to compute IDSTART, then move there.
  //
  // Jun 9 2018: 
  //  + for batch job and NLIBID in header, auto-compute NSKIP_LIBID.
  //  + use fgets to skip LIBIDs much more quickly
  //
  // Aug 19 2018: use JOBID and NJBOTOT to compute NSKIP so that
  //    entire SIMLIB range is covered with small NGENTOT.
  //
  // Dec 09 2020: return immediately if QUIT_REWIND option is set to
  //     read SIMLIB once and stop. See SIMLIB_MSKOPT += 4.

  int IDSTART  = INPUTS.SIMLIB_IDSTART ;
  int IDLOCK   = INPUTS.SIMLIB_IDLOCK ;
  int NLIBID   = SIMLIB_GLOBAL_HEADER.NLIBID ;
  int JOBID    = INPUTS.JOBID ; // batch JOBID
  int NJOBTOT  = INPUTS.NJOBTOT ; // totoal number of batch jobs

  int  MSKOPT = INPUTS.SIMLIB_MSKOPT ;
  bool QUIT_NOREWIND = ( (MSKOPT & SIMLIB_MSKOPT_QUIT_NOREWIND)>0 );
  if ( QUIT_NOREWIND ) { return; }

  int NOPT, IDSEEK, MAXRANSTART, NSKIP_LIBID=-9, NSKIP_EXTRA=0 ;
  int DOSKIP, NREAD, NREPEAT, MXREPEAT, NTMP, NLIBID_EXTRA ;
  time_t t0=time(NULL), t1=time(NULL); 
  double XSKIP, XTMP, flatRan ;
  char LINE[100];
  char fnam[] = "SIMLIB_findStart" ;

  // -------------------- BEGIN -------------

  DOSKIP = 0 ;    IDSEEK = -9 ;

  // check user-input SIMLIB_MAXRANSTART
  MAXRANSTART = INPUTS.SIMLIB_MAXRANSTART ; // default
  if ( MAXRANSTART > 0 ) {
    flatRan     = unix_getRan_Flat1(0) ;
    XSKIP       = (double)MAXRANSTART * flatRan ;
    NSKIP_LIBID = (int)XSKIP + 1 ;
    DOSKIP = 1;
    printf("\t SIMLIB MAXRANSTART after  %d  LIBIDs ", 
	   NSKIP_LIBID ); 
  }

  // for batch job, autom-compute NSKIP 
  if ( NJOBTOT > 0  &&  NLIBID > 100 ) { 
    flatRan     = unix_getRan_Flat1(0) ;
    XTMP        = (double)NLIBID / (double)NJOBTOT;    
    NTMP        = (int)XTMP ;
    NLIBID_EXTRA = NTMP - INPUTS.NGENTOT_LC ; // Number of extra LIBIDs 

    // if NGENTOT_LC doesn't use the whole SIMLIB range,
    // pick a random start point 
    if ( INPUTS.NGENTOT_LC>0 && NLIBID_EXTRA > 5 ) {
      NSKIP_EXTRA = (int)( flatRan * (double)NLIBID_EXTRA );
    }
 
    NSKIP_LIBID = (JOBID-1) * (int)XTMP  + NSKIP_EXTRA;

    DOSKIP = 1;
    printf("\t SIMLIB BATCH-MODE START at %d of %d LIBIDs ", 
	   NSKIP_LIBID, NLIBID ); 
  }


  if ( IDSTART > 0  )  { 
    IDSEEK = IDSTART - 1 ; 
    printf("\t SIMLIB IDSTART:  %d  (start at this LIBID) ", IDSTART );
    DOSKIP = 1;
  }

  if ( IDLOCK  >= 0 ) {
    IDSEEK = IDLOCK ;
    printf("\t SIMLIB IDLOCK:   %d  (will use only this LIBID) ", IDLOCK );
    DOSKIP = 1;
  }

  
  // count number of valid options. Only 0 or 1 allowed.
  NOPT=0;
  if ( IDSTART     > 0 ) { NOPT++ ; }
  if ( IDLOCK      > 0 ) { NOPT++ ; }
  if ( NSKIP_LIBID > 0 ) { NOPT++ ; }

  if ( NOPT   == 0 ) { return ; } // do nothing ==> start at first LIBID

  if ( NOPT > 1 ) {
    sprintf ( c1err, "Cannot combine SIMLIB-start options:");
    sprintf ( c2err," (IDSTART=%d  IDLOCK=%d JOBID=%d)", 
	      IDSTART, IDLOCK, JOBID );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

 
  if ( DOSKIP ) { t0 = time(NULL) ; }


  GENLC.SIMLIB_ID = NREAD = 0;

  NREPEAT  = INPUTS.SIMLIB_NREPEAT  ; INPUTS.SIMLIB_NREPEAT  = 1 ;
  MXREPEAT = INPUTS.SIMLIB_MXREPEAT ; INPUTS.SIMLIB_MXREPEAT = 1 ;
  fflush(stdout);

  
  // skip fixed number of LIBIDs
  while ( NREAD < NSKIP_LIBID ) {
    fgets(LINE, 40, fp_SIMLIB) ;
    if ( strstr(LINE,"END_LIBID:") != NULL ) { NREAD++; }
  }
  

  // search for specific LIBID
  while ( IDSEEK > SIMLIB_HEADER.LIBID ) {   
    SIMLIB_READ_DRIVER();
    if ( SIMLIB_HEADER.NWRAP > 0 ) {
      sprintf(c1err,"Wrapped around SIMLIB without finding start.");
      sprintf(c2err,"Try again with smaller "
	      "SIMLIB_IDSTART or MAXRANSTART value.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }
  }


 
  if ( DOSKIP ) { 
    t1 = time(NULL); 
    printf(" (%d seconds)\n", (int)(t1-t0) ); fflush(stdout); 
  }

  // -----------------------------------------------
  // restore REPEAT
  INPUTS.SIMLIB_NREPEAT  = NREPEAT ;
  INPUTS.SIMLIB_MXREPEAT = MXREPEAT ;

  GENLC.SIMLIB_IDLOCK = INPUTS.SIMLIB_IDLOCK ; // set lock on LIBID


  return ;

}  // end of SIMLIB_findStart


// **********************************************
void  GENFILTERS_CHECK(void) {

  // Created Mar 2013.
  // check that each GENFILTER (INPUTS.GENFILTERS) is in the 
  // list of filters from the simlib 
  // If not, then abort.

  int  NFILT_GEN, NFILT_SIMLIB, ifilt, ifilt_obs,  NFILT_UNDEFINED ;
  char cfilt[4], CFILT_UNDEFINED[MXFILTINDX];
  char fnam[] = "GENFILTERS_CHECK" ;

  // ---------- BEGIN ----------------

  NFILT_GEN    = strlen(INPUTS.GENFILTERS);
  NFILT_SIMLIB = GENLC.NFILTDEF_SIMLIB ;
  
  NFILT_UNDEFINED = 0 ;
  CFILT_UNDEFINED[0] = 0 ;

  for(ifilt=0; ifilt < NFILT_GEN ; ifilt++ ) {
    sprintf(cfilt, "%c ",  INPUTS.GENFILTERS[ifilt] );
    ifilt_obs = filtindx_( cfilt, strlen(cfilt) ); 

    if ( GENLC.IFILTINVMAP_SIMLIB[ifilt_obs] < 0 ) {
      NFILT_UNDEFINED++ ;
      strcat(CFILT_UNDEFINED, cfilt);
    }

    /*
    printf(" xxx check GENFILTER = '%s' -> ifilt_obs = %d (inv=%d)\n", 
	   cfilt, ifilt_obs, ifiltinvmap ); fflush(stdout);
    */
  }

  if ( NFILT_UNDEFINED > 0 ) {
    sprintf(c1err,"%d undefined GENFILTERS (%s)", 
	    NFILT_UNDEFINED, CFILT_UNDEFINED);
    sprintf(c2err,"not part of SIMLIB filters (%s)",
	    SIMLIB_GLOBAL_HEADER.FILTERS );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

} // end of GENFILTERS_CHECK


// *************************************************
int SIMLIB_read_templateNoise(char *FIELD, char *whatNoise, char **wdlist) {

  // Feb 2014
  // For *whatNoise = "SKY" or "CCD" or "ZPT", read correlated template 
  // noise for each defined filter and store in 
  // SIMLIB_TEMPLATE_SKY[CCD]SIG array.
  // The template noise values must be listed in the same order
  // as for the FILTERS key.
  // Note that this is different from the T: keys where
  // the template noise is UN-correlated among epochs.
  // **wdlist is list of words passed from file;
  // Function returns number of words read (0 or NFILT)
  //
  //     HISTORY
  //
  // Mar 02 2018: check option to ignore template noise
  // Sep 02 2021: 
  //   + refactor to read from wdlist instead of reading file from
  //     fp_SIMLIB; return number or words read
  //

  double noise, *ptrNoise, validNoise_min, validNoise_max;
  int  NRD=0, NFILT, ifilt, ifilt_obs, IFIELD, IGNORE  ;
  char fnam[] = "SIMLIB_read_templateNoise" ;
  
  // ------------- BEGIN -----------

  IGNORE = (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_IGNORE_TEMPLATENOISE);
  if ( IGNORE ) { return NRD ; }

  // determine local FIELD index.
  IFIELD = IFIELD_OVP_SIMLIB(0,FIELD);

  // if we did not find IFIELD, then increment for new field
  if ( IFIELD < 0 ) {
    SIMLIB_TEMPLATE.NFIELD_OVP++ ;
    IFIELD = SIMLIB_TEMPLATE.NFIELD_OVP ;
  }
  
  
  if ( IFIELD < 0 || IFIELD >= MXFIELD_OVP ) {
    sprintf(c1err,"Invalid IFIELD=%d for FIELD=%s, LIBID=%d", 
	    IFIELD, FIELD, SIMLIB_HEADER.LIBID );
    sprintf(c2err,"IFIELD must be 0 to %d", MXFIELD_OVP-1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // load FIELDNAME to keep track of overlaps.
  sprintf(SIMLIB_TEMPLATE.FIELDNAME[IFIELD],"%s", FIELD); 

  // -------------------------------------
  validNoise_min = validNoise_max = 0.0 ;
  ptrNoise  = &validNoise_min  ; // avoid compile warning

  if ( strcmp(whatNoise,"SKYSIG") == 0 )  { 
    ptrNoise = SIMLIB_TEMPLATE.SKYSIG[IFIELD] ; 
    validNoise_min = 0.0;  validNoise_max = 1.0E5 ;
  }

  else if ( strcmp(whatNoise,"CCDSIG") == 0 )  { 
    ptrNoise = SIMLIB_TEMPLATE.READNOISE[IFIELD] ; 
    validNoise_min = 0.0;  validNoise_max = 100. ;
  }

  else if ( strcmp(whatNoise,"ZPT") == 0 )  { 
    ptrNoise = SIMLIB_TEMPLATE.ZPT[IFIELD] ; 
    validNoise_min = 0.0;  validNoise_max = 50. ;
  }

  else {
    sprintf(c1err,"Invalid whatNoise = '%s' ", whatNoise);
    sprintf(c2err,"Valid whatNoise = SKY or CCD or ZPT");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  NFILT = GENLC.NFILTDEF_SIMLIB ;

  if ( NFILT <= 0 ) {
    sprintf(c1err,"Cannot read template noise before reading filters.");
    sprintf(c2err,"Check FILTERS key in simlib header");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  char varName[60] ;

  for(ifilt= 0 ; ifilt < NFILT; ifilt++ ) {
    // xxx mark    readdouble(fp_SIMLIB, 1, &noise);
    sscanf(wdlist[ifilt], "%le", &noise ); NRD++ ;
    ifilt_obs = GENLC.IFILTMAP_SIMLIB[ifilt] ;
    ptrNoise[ifilt_obs] = noise ;

    sprintf(varName,"TEMPLATE_%s(%c-%d)", 
	    whatNoise, FILTERSTRING[ifilt_obs], ifilt );
    checkval_D(varName, 1, &noise, validNoise_min, validNoise_max );

    /*    
    printf("\t xxx %s(ifilt,ifiltobs=%2d,%2d) = %f \n",
	   whatNoise, ifilt, ifilt_obs, noise );  fflush(stdout);  
    */
  }

  SIMLIB_TEMPLATE.USEFLAG |= 1;

  return NRD ;
  
} // end of SIMLIB_read_templateNoise

// *************************************************
int IFIELD_OVP_SIMLIB(int OPT, char *FIELD) {

  // for input *FIELD, return IFIELD index that matches
  // list in SIMLIB_TEMPLATE.FIELDNAME.
  // Return -9 if there is no match.
  //
  // OPT=0 --> read stage, just check FIELD name
  // OPT=1 --> use stage, so check global ZPT

  int  IFIELD, NFIELD, i  ;
  //  char fnam[] = "IFIELD_OVP_SIMLIB" ;
  // ---------------- BEGIN ------------

  IFIELD = -9;

  if ( strcmp(FIELD,ALLFIELDS) == 0 )  { IFIELD=0; }

  if ( OPT == 1 ) {
    // set global IFIELD if global ZPT is set. But do NOT return
    // to allow override below.
    int ifield    = 0;
    int ifilt_obs = GENLC.IFILTMAP_SIMLIB[0] ;
    if ( SIMLIB_TEMPLATE.ZPT[ifield][ifilt_obs] > 10.0 ) { IFIELD = 0 ; }
  }

  // -------------------------
  // if we get here then check list of fields.
  // start at ifield=1 since 0 is reserved for global ALLFIELDS

  NFIELD = SIMLIB_TEMPLATE.NFIELD_OVP ;
  for ( i=1; i <= NFIELD; i++ ) {
    if ( strcmp(SIMLIB_TEMPLATE.FIELDNAME[i],FIELD) == 0 ) 
      { IFIELD = i; }
  } 

  return IFIELD ;

} // end of IFIELD_OVP_SIMLIB



// ************************************
double get_SIMLIB_fluxerrScale_LEGACY(int ifiltobs, double SNR) {

  // Dec 5, 2011 R.Kessler
  // Return error correction on flux based on FLUXERR_COR map
  // in the simlib file. If there is no map for this filter
  // then return 1.
  //
  // Feb 24, 2012: if SNR is outside map, return value at min or max edge.
  // Aug 15, 2017: arrays are 0 to M-1 (instead of 1 to M)
  //               ERRCOR -> SCALE

  double SCALE, LOGSNR, LOGSNR_MIN, LOGSNR_MAX;
  int    ifilt, M, LDMP ;
  int    OPT_INTERP = 1;  // 1=linear, 2=quadratic
  char   fnam[] = "get_SIMLIB_fluxerrScale_LEGACY" ;

  // -------------- BEGIN ---------------

  LDMP = 0 ;

  SCALE = 1.0 ;
  if ( SIMLIB_FLUXERR_COR.USE == 0 ) { goto END; }

  ifilt  = GENLC.IFILTINVMAP_SIMLIB[ifiltobs] ;
  M      = SIMLIB_FLUXERR_COR.MAPSIZE[ifilt];
  
  if ( M <= 0 ) { goto END; }

  // return value at map edge if SNR is outside the range of the map.
  LOGSNR     =  log10(SNR);
  LOGSNR_MIN =  SIMLIB_FLUXERR_COR.LOG10SNR[ifilt][0] ;
  LOGSNR_MAX =  SIMLIB_FLUXERR_COR.LOG10SNR[ifilt][M-1] ;
  if ( LOGSNR < LOGSNR_MIN ) 
    {  SCALE = SIMLIB_FLUXERR_COR.SCALE[ifilt][0] ;    goto END;  }

  if ( LOGSNR > LOGSNR_MAX ) 
    {  SCALE = SIMLIB_FLUXERR_COR.SCALE[ifilt][M-1] ;    goto END;  }

  // LOG(SNR) is defined by the map ... find nearest bin and
  // get error correction-scale (SCALE) from linear interpolation.

  SCALE = interp_1DFUN(OPT_INTERP, LOGSNR, M
		       ,SIMLIB_FLUXERR_COR.LOG10SNR[ifilt] 
		       ,SIMLIB_FLUXERR_COR.SCALE[ifilt] 
		       ,fnam );
 END:  
  if ( LDMP ) {
    printf(" xxxx ERR SCALE(%c) = %6.4f for SNR=%f  LOG10(SNR)=%f\n",
	   FILTERSTRING[ifiltobs], SCALE, SNR, LOGSNR);
  }

  return(SCALE);

} // end of get_SIMLIB_fluxerrScale_LEGACY


// ************************************
void SIMLIB_READ_DRIVER(void) {

  // Created Aug 2017
  // Re-factored SIMLIB-reader to separate reading and processing
  // (i.e., untangle spagetti code). This refactor allows reading
  // from traditional TEXT file, and also from a binary file such
  // as ROOT or FITS.
  //
  // Nov 22 2019: init REPEAT=0

  int  REPEAT = 0 ;
  char fnam[] = "SIMLIB_READ_DRIVER" ;

  // ------------------ BEGIN ------------------

  // Begin refactored code.
  // read next cadence and load SIMLIB_OBS_RAW array

  // check for option to repeat Cadence to save time 
  GENLC.NGEN_SIMLIB_ID++ ;
  REPEAT = USE_SAME_SIMLIB_ID(2);

  if ( !REPEAT ) {  // process next cadence

    // read next cadence from SIMLIB/Cadence file (any format)
    SIMLIB_readNextCadence_TEXT(); 

    // expand SPECTROGRAPH keys (MJD,TEXPOSE) to include 
    // each synthetic band as if they were in the SIMLIB file
    SIMLIB_addCadence_SPECTROGRAPH(); 
  }

  // TAKE_SPECTRUM uses PEAKMJD and thus must always be
  // computed, even if cadence is repeated.
  SIMLIB_TAKE_SPECTRUM(); 

  // transfer OBS_RAW list to GEN_OBS list; cuts and MJD-sorting
  SIMLIB_prepCadence(REPEAT);

  return ;

} // end SIMLIB_READ_DRIVER



// ==============================================
void  SIMLIB_readNextCadence_TEXT(void) {

  // Created Aug 7 2017
  // Read next cadence from TEXT SIMLIB file,
  // and store in the following structures:
  //   + SIMLIB_HEADER 
  //   + SIMLIB_OBS_RAW 
  //   + SIMLIB_TEMPLATE
  //
  // Notes:
  //  + SIMLIB text file has already been opened, with global 
  //    file handle fp_SIMLIB
  //  + Defining filter as sum of other filters is NO longer
  //    supported. i.e., ZPT='1+2+3+4' is no longer valid.
  // 
  // to do: apply GENRANGE_MJD cut ? And keep track of NOBS
  //
  // Dec 1 2016: 
  //  fix to work properly when LCLIB NREPEAT=0; see NOBS_FOUND_ALL
  //    
  // Jan 3 2018: use parse_SIMLIB_IDplusNEXPOSE() to read IDEXPT & NEXPOSE
  // Sep 17 2019: rewind on EOF so that END_OF_SIMLIB: key is optional.
  // May 15 2020: don't read SPECTROGRPAH key unless
  //     SPECTROGRAPH_USEFLAG is set.
  //
  // Feb 05 2021: fix to handle mutliple FIELDs
  // Feb 25 2021: c_get[80] -> c_get[200] to allow for long comment strings
  //
  // Sep 01 2021: refactor to read entire line with fgets;
  //              hopefully speeds up Cori batch jobs.
  //
#define MXWDLIST_SIMLIB 20  // max number of words per line to read

  int ISMODEL_SIMLIB =  (INDEX_GENMODEL == MODEL_SIMLIB);
  int ID, NOBS_EXPECT, NOBS_FOUND, NOBS_FOUND_ALL, ISTORE=0 ;
  int APPEND_PHOTFLAG, ifilt_obs, DONE_READING, NWD, iwd, IWD ;
  int NTRY, USEFLAG_LIBID, USEFLAG_MJD, OPTLINE, NTMP, NFIELD ;
  int NOBS_SKIP, SKIP_FIELD, SKIP_APPEND, OPTLINE_REJECT, NMAG_notZeroFlux;
  int OPTMASK, noTEMPLATE ;
  double TEXPOSE, TSCALE ;
  bool  FOUND_SPECTROGRAPH, FOUND_EOF, FOUND_ENDKEY, SKIP_MJD ;
  bool  ISKEY, ISKEY_S, ISKEY_TEMPLATE;
  double PIXSIZE, TEXPOSE_S, MJD, MAG ;
  char wd0[200], wd1[200], ctmp[80], *BAND, cline[400], *pos ;
  char WDLIST[MXWDLIST_SIMLIB][200], *ptrWDLIST[MXWDLIST_SIMLIB];
  char *FIELD = SIMLIB_HEADER.FIELD;
  char *TEL   = SIMLIB_HEADER.TELESCOPE ;
  char sepKey[] = " ";
  char fnam[] = "SIMLIB_readNextCadence_TEXT" ;

  // ------------ BEGIN --------------

  SIMLIB_HEADER.NWRAP = NTRY = 0 ; // reset wrap before START 
  SIMLIB_OBS_RAW.NOBS = SIMLIB_OBS_RAW.NOBS_READ = 0 ;
  SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH = 0 ;

  for(iwd=0; iwd < MXWDLIST_SIMLIB; iwd++ ) { ptrWDLIST[iwd] = WDLIST[iwd]; }


 START:

  init_SIMLIB_HEADER();
  NOBS_EXPECT = NOBS_FOUND = NOBS_FOUND_ALL = ISTORE = 0 ;
  USEFLAG_LIBID =USEFLAG_MJD = 0 ;
  DONE_READING = NOBS_SKIP = SKIP_FIELD = SKIP_APPEND = APPEND_PHOTFLAG = 0 ;
  NMAG_notZeroFlux = 0 ;
  SIMLIB_LIST_forSORT.MJD_LAST = -9.0 ;
  SIMLIB_TEMPLATE.NFIELD_OVP = 0;
  NTRY++ ;

  if ( NTRY >= MXREAD_SIMLIB ) {
    sprintf(c1err,"NTRY=%d exceeds bound of MXREAD_SIMLIB=%d",
	    NTRY, MXREAD_SIMLIB );
    sprintf(c2err,"Check SIMLIB file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  // - - - - - - - start reading SIMLIB - - - - - - - - - 
  int NLINE=0;

  while ( !DONE_READING ) {

    NLINE++ ;
    cline[0] = 0 ;   FOUND_EOF = false ;
    if ( fgets(cline, 380, fp_SIMLIB) == NULL ) { FOUND_EOF = true; }

    // skip comment character
    if ( commentchar(cline) ) { continue; }

    // remove line feed    
    if ( (pos=strchr(cline,'\n') ) != NULL )  { *pos = '\0' ; }

    // note that splitString2 is fast, but destroys cline
    splitString2(cline, sepKey, MXWDLIST_SIMLIB, &NWD, ptrWDLIST);

    // check end of file, or end of simlib keywork -> rewind
    FOUND_ENDKEY = ( strcmp(WDLIST[0],"END_OF_SIMLIB:") == 0 );
    if ( FOUND_EOF || FOUND_ENDKEY ) {
      // check SIMLIB after 5 passes to avoid infinite loop
      ENDSIMLIB_check();
      if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENRANDOM ) {
	snana_rewind(fp_SIMLIB, INPUTS.SIMLIB_OPENFILE,
		     INPUTS.SIMLIB_GZIPFLAG);
	SIMLIB_HEADER.NWRAP++ ; 
	SIMLIB_HEADER.LIBID = SIMLIB_ID_REWIND ; 
	NOBS_FOUND = NOBS_FOUND_ALL = USEFLAG_LIBID = USEFLAG_MJD = 0 ;
	ISTORE = 0;
      }
      continue ;
    }  // end REWIND

    ISKEY_S        = ( strcmp(WDLIST[0], "S:") == 0 ) ;
    ISKEY_TEMPLATE = ( strstr(WDLIST[0], "TEMPLATE_") != NULL ) ;

    for(iwd=0; iwd < NWD; iwd++ ) {

      wd0[0] = wd1[0] = 0;
      sprintf(wd0,"%s", WDLIST[iwd] );
      if ( NWD > 1 ) { sprintf(wd1,"%s", WDLIST[iwd+1] ); }

      if ( strcmp(wd0,"LIBID:") == 0 ) {
	sscanf(wd1, "%d", &ID ); 
	SIMLIB_HEADER.LIBID = ID ; 
	sprintf(SIMLIB_HEADER.LIBNAME, "LIB%5.5d", ID );
	USEFLAG_LIBID = ACCEPT_FLAG ;
	NFIELD = 0 ;
	iwd++ ; continue ;
      }

      if ( strcmp(wd0,"RA:") == 0 ) 
	{ sscanf(wd1, "%le", &SIMLIB_HEADER.RA );  iwd++ ; continue; }
      else if ( strcmp(wd0,"DEC:") == 0 ) 
	{ sscanf(wd1, "%le", &SIMLIB_HEADER.DEC ); iwd++ ; continue; }
      else if ( strcmp(wd0,"DECL:") == 0 ) 
	{ sscanf(wd1, "%le", &SIMLIB_HEADER.DEC ); iwd++ ; continue; }
      
      else if ( strcmp(wd0,"SUBSURVEY:")==0 ) 
	{ sscanf(wd1,"%s",SIMLIB_HEADER.SUBSURVEY_NAME); iwd++; continue;}
      
      else if ( strcmp(wd0,"FIELD:") == 0 )  { 

	char tmp_field[40];
	sscanf(wd1, "%s", tmp_field ); 
	sprintf(SIMLIB_HEADER.FIELDLIST_OVP[NFIELD], "%s", tmp_field);

	NFIELD++ ;   SIMLIB_HEADER.NFIELD_OVP = NFIELD;
	if ( NFIELD == 1 ) 
	  { sprintf(FIELD, "%s", tmp_field) ; }
	else
	  { strcat(FIELD,"+"); strcat(FIELD,tmp_field); }
	
	SKIP_FIELD = ( SKIP_SIMLIB_FIELD(FIELD) &&
		       (INPUTS.SIMLIB_FIELDSKIP_FLAG ==0 ) ) ;

	iwd++ ; continue;
      }
      
      else if ( strcmp(wd0,"PIXSIZE:") == 0 )  
	{ sscanf(wd1,"%le", &SIMLIB_HEADER.PIXSIZE ); iwd++; continue;}
      
      else if ( strcmp(wd0,"TELESCOPE:") == 0 ) 
	{ sscanf(wd1,"%s",SIMLIB_HEADER.TELESCOPE ); iwd++; continue;}
  
      else if ( strcmp(wd0,"MWEBV:") == 0  )
	{ sscanf(wd1,"%le", &SIMLIB_HEADER.MWEBV ); iwd++ ; continue; }

      else if ( strcmp(wd0,"CCD:") == 0 || strcmp(wd0,"CCDNUM:") == 0 ) 
	{ sscanf(wd1,"%d",&SIMLIB_HEADER.CCDNUM);  iwd++ ; continue; }

      // read optional header keys for FAKEID option
      else if ( strcmp(wd0,"GALID:") == 0 )  
	{ sscanf(wd1, "%lld", &SIMLIB_HEADER.GALID ); iwd++ ; continue; }
      else if ( strcmp(wd0,"FAKEID:") == 0 )  
	{ sscanf(wd1, "%d", &SIMLIB_HEADER.FAKEID ); iwd++ ; continue; }

      // check for APPEND option --> last MJDs are not sorted.
      else if ( strcmp(wd0,"APPEND:") == 0 ) { 
	sscanf(wd1, "%d", &APPEND_PHOTFLAG ); 
	SKIP_APPEND = 0 ; // to do: replace with user input flag
	iwd++ ; continue;
      }
    
      // - - - -
      // note that NOBS can exceed MXEPSIM ... 
      // as long as it's less than MXOBS_SIMLIB
      if ( strcmp(wd0,"NOBS:") == 0 )   {  
	sscanf(wd1, "%d", &NOBS_EXPECT ); 
	SIMLIB_HEADER.NOBS = NOBS_EXPECT ;
	
	if ( NOBS_EXPECT >= MXOBS_SIMLIB ) {
	  sprintf(c1err,"NOBS=%d exceeds bound for LIBID=%d.", 
		  NOBS_EXPECT, ID);
	  sprintf(c2err,"Check bound: MXOBS_SIMLIB = %d", 
		  MXOBS_SIMLIB );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
	}
	iwd++; continue ;
      }
    
      // check for correlated template noise ; fill SIMLIB_TEMPLATE
      if ( ISKEY_TEMPLATE && iwd==0 ) {
	int iwd_start = iwd;
	if ( strcmp(wd0,"TEMPLATE_SKYSIG:") == 0 )
	  { iwd += SIMLIB_read_templateNoise(FIELD,"SKYSIG", &ptrWDLIST[1] ); } 
	else if ( strcmp(wd0,"TEMPLATE_CCDSIG:") == 0 )
	  { iwd += SIMLIB_read_templateNoise(FIELD,"CCDSIG", &ptrWDLIST[1] ); } 
	else if ( strcmp(wd0,"TEMPLATE_ZPT:") == 0 )
	  { iwd += SIMLIB_read_templateNoise(FIELD,"ZPT", &ptrWDLIST[1] ); }    

	// read spectrograph exposure time for template
	OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;
	noTEMPLATE = ( OPTMASK & SPECTROGRAPH_OPTMASK_noTEMPLATE ) ;
	ISKEY      = ( strcmp(wd0,"TEMPLATE_TEXPOSE_SPECTROGRAPH:")==0);
	if ( ISKEY && (noTEMPLATE==0) ) {
	  TSCALE=INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE ;
	  sscanf(wd1,"%le", &TEXPOSE ); iwd++ ;
	  SIMLIB_TEMPLATE.TEXPOSE_SPECTROGRAPH = TEXPOSE * TSCALE ;
	}
	// ?? if ( iwd > iwd_start ) { continue; }
      } // end ISKEY_TEMPLATE

      // check for optional GENRANGE or cut keys in header.
      parse_SIMLIB_GENRANGES( &ptrWDLIST[iwd] ); 
    
      // ------------------------------------------------
      // check for epochs
      OPTLINE = 0;
      if ( ISKEY_S ) {
	NOBS_FOUND_ALL++ ;
	if ( USEFLAG_LIBID == ACCEPT_FLAG ) { OPTLINE = OPTLINE_SIMLIB_S; }
      }
    
      FOUND_SPECTROGRAPH = 
	( SPECTROGRAPH_USEFLAG && strcmp(wd0,"SPECTROGRAPH:")==0 );
      if ( FOUND_SPECTROGRAPH && USEFLAG_LIBID==ACCEPT_FLAG )
	{ OPTLINE = OPTLINE_SIMLIB_SPECTROGRAPH ; }

      // always check reasons to reject (header cuts, FIELD, APPEND ...)
      OPTLINE_REJECT = ( USEFLAG_LIBID == REJECT_FLAG || 
			 SKIP_FIELD || SKIP_APPEND ) ;

      if ( OPTLINE && OPTLINE_REJECT )  {    
	// xxx MJD line in already rejected LIBID --> read rest of line 
	// xxx mark delete	fgets(cline, 180, fp_SIMLIB) ;
	if ( SKIP_FIELD ) { NOBS_SKIP++ ; }
      }
      else if ( OPTLINE == OPTLINE_SIMLIB_S )  { 
	NOBS_FOUND++ ; 	IWD = iwd;  SKIP_MJD = false;

	/* xxxxxx mark delete 
	if ( INPUTS.DEBUG_FLAG == 903 ) {
	  if ( MJD < GENLC.MJD_RANGE[0] || MJD > GENLC.MJD_RANGE[1] ) {
	    SKIP_MJD = true ;
	  }
	}
	if ( SKIP_MJD ) { iwd=NWD; break; }
	xxxxxx */

	SIMLIB_OBS_RAW.OPTLINE[ISTORE] = OPTLINE ;

	IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.MJD[ISTORE]);

	IWD++; sscanf(WDLIST[IWD], "%s", ctmp );
	parse_SIMLIB_IDplusNEXPOSE(ctmp, 
				   &SIMLIB_OBS_RAW.IDEXPT[ISTORE],
				   &SIMLIB_OBS_RAW.NEXPOSE[ISTORE] );

	IWD++; sscanf(WDLIST[IWD], "%s" , SIMLIB_OBS_RAW.BAND[ISTORE]     );
	IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.CCDGAIN[ISTORE]  );
	IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.READNOISE[ISTORE]);
	IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.SKYSIG[ISTORE]   );

	if ( SIMLIB_GLOBAL_HEADER.NEA_PSF_UNIT ) 
	  { IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.NEA[ISTORE] ); }
	else {
	  IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.PSFSIG1[ISTORE] );
	  IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.PSFSIG2[ISTORE] );
	  IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.PSFRATIO[ISTORE]);

	  // if NEA is here, but user forgets "PSF_UNIT: NEA_PIXEL" in header,
	  // this trap will hopefully abort.
	  checkval_D("PSF1(readNextCadence)", 1, 
		     &SIMLIB_OBS_RAW.PSFSIG1[ISTORE], 0.0, 30.0 ) ;
	}
	IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.ZPTADU[ISTORE]   ); 
	checkval_D("ZPT(readNextCadence)", 1, 
		   &SIMLIB_OBS_RAW.ZPTADU[ISTORE], 5.0, 50.0 ) ;

	IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.ZPTERR[ISTORE]   );  
	IWD++; sscanf(WDLIST[IWD], "%le", &SIMLIB_OBS_RAW.MAG[ISTORE]      );
	iwd = NWD; 

	if ( INPUTS.FORCEVAL_PSF > 0.001 )  // Sep 2020
	  { SIMLIB_OBS_RAW.PSFSIG1[ISTORE] = INPUTS.FORCEVAL_PSF;  }

	// check MAG column for SIMLIB model (Nov 2019)
	MAG = SIMLIB_OBS_RAW.MAG[ISTORE];
	if ( MAG < MAG_ZEROFLUX-0.001 ) { NMAG_notZeroFlux++ ; }

	// convert filter-char to integer		   
	BAND      = SIMLIB_OBS_RAW.BAND[ISTORE] ;
	ifilt_obs = INTFILTER(BAND);
	SIMLIB_OBS_RAW.IFILT_OBS[ISTORE] = ifilt_obs ;

	// update few header items for each epoch since
	// these item can be changed at any epoch.
	sprintf(SIMLIB_OBS_RAW.FIELDNAME[ISTORE], "%s", FIELD );
	sprintf(SIMLIB_OBS_RAW.TELESCOPE[ISTORE], "%s", TEL);
	PIXSIZE = SIMLIB_HEADER.PIXSIZE ;
	SIMLIB_OBS_RAW.PIXSIZE[ISTORE] = PIXSIZE ;

	// set 'not from spectrograph' values
	SIMLIB_OBS_RAW.IFILT_SPECTROGRAPH[ISTORE]   = -9 ;
	SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[ISTORE] =  0.0 ; 
	SIMLIB_OBS_RAW.INDX_TAKE_SPECTRUM[ISTORE]   = -9 ; 

	ISTORE++ ;
      }
      else if ( OPTLINE == OPTLINE_SIMLIB_SPECTROGRAPH  )  { 
	
	NOBS_FOUND++ ;

	sscanf(WDLIST[iwd+1], "%le", &MJD );
	sscanf(WDLIST[iwd+2], "%le", &TEXPOSE_S );
	iwd = NWD;

	// increment sparse list so that SPECTROGRAPH entries
	// can be found later (in SIMLIB_addCadence_SPECTROGRAPH)
	// without looping thru entire list.
	NTMP = SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH ;
	SIMLIB_OBS_RAW.OBSLIST_SPECTROGRAPH[NTMP] = ISTORE;
	SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH++ ;
	
	// store few things at this ISTORE location
	SIMLIB_OBS_RAW.OPTLINE[ISTORE]    = OPTLINE ;
	SIMLIB_OBS_RAW.MJD[ISTORE]        = MJD ;
	SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[ISTORE] = TEXPOSE_S ;      
	SIMLIB_OBS_RAW.BAND[ISTORE][0] = 0 ;
	sprintf(SIMLIB_OBS_RAW.FIELDNAME[ISTORE], "%s", FIELD );

	ISTORE++ ;

      } // end OPTLINE == OPTLINE_SIMLIB_SPECTROGRAPH

      // check APPEND_PHOTFLAG to NOT sort this MJD (Jan 2018)
      if ( OPTLINE && (OPTLINE_REJECT==0) )  {    
	if ( APPEND_PHOTFLAG > 0 ) { SIMLIB_HEADER.NOBS_APPEND++ ; }
	SIMLIB_OBS_RAW.APPEND_PHOTFLAG[ISTORE] = APPEND_PHOTFLAG ;
      }

      // after first OBS is found we are done with header.
      // -> check for random RA,DEC shift
      // -> apply header cuts on ID, redshift, etc...
      if ( NOBS_FOUND_ALL == 1 ) { 
	SIMLIB_randomize_skyCoords();
	USEFLAG_LIBID = keep_SIMLIB_HEADER(); 
	if ( USEFLAG_LIBID != ACCEPT_FLAG && SIMLIB_HEADER.NWRAP==0 )
	  { SIMLIB_GLOBAL_HEADER.NLIBID_VALID-- ; }	
      }

      // stop reading when we reach the end of this LIBID
      if ( strcmp(wd0, "END_LIBID:") == 0 ) { 
	if ( USEFLAG_LIBID == ACCEPT_FLAG ) 
	  { DONE_READING = 1 ; }
	else
	  { goto START ; }      // read another 
      }

    } // end while wd loop with 
  }   // end !DONE_READING
  
  // ---------------------

  /*
  printf(" xxx %s  NOBS(EXPECT,FOUND,SKIP,not0) = %d,%d,%d,%d   FIELD=%s \n",
	 fnam, NOBS_EXPECT, NOBS_FOUND, NOBS_SKIP, NMAG_notZeroFlux,
	 FIELD );
  */

  if ( INPUTS.SIMLIB_FIELDSKIP_FLAG==0 && ( NOBS_EXPECT==NOBS_SKIP) ) 
    { goto START ; }


  SIMLIB_OBS_RAW.NOBS      = ISTORE ;      // can change with SPECTROGRAPH
  SIMLIB_OBS_RAW.NOBS_READ = ISTORE ;  // won't change for this cadence


  NOBS_EXPECT -= NOBS_SKIP ;
  if ( NOBS_EXPECT != NOBS_FOUND ) {
    sprintf(c1err,"Found %d observations in LIBID %d", 
	    NOBS_FOUND,  SIMLIB_HEADER.LIBID );
    sprintf(c2err,"But expected NOBS = %d from simlib header (NOBS_SKIP=%d)", 
	    NOBS_EXPECT, NOBS_SKIP );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }


  if ( ISMODEL_SIMLIB  && NMAG_notZeroFlux==0 ) {
    sprintf(c1err,"All MAG=%.1f (i.e., FLUX=0) for LIBID=%d .",
	    MAG_ZEROFLUX, SIMLIB_HEADER.LIBID);
    sprintf(c2err,"For SIMLIB model, some valid MAG values are required.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  return ;

} // end SIMLIB_readNextCadence_TEXT


// ==============================================
void SIMLIB_randomize_skyCoords(void) {

  // Created Aug 24 2018
  //
  // Check option to randomly shift 
  //   SIMLIB_HEADER.RA 
  //   SIMLIB_HEADER.DEC
  // To avoid fixed set of RA,DEC from SIMLIB.
  //
  // Initial function is just a do-nothing shell;
  // waiting to get C function that picks randonly
  // within healPix area.

  //  char fnam[] = "SIMLIB_randomize_skyCoords" ;

  // --------------- BEGIN -------------

  return ;

} // end SIMLIB_randomize_skyCoords


// ==============================================
int keep_SIMLIB_HEADER(void) {

  // Created Aug 2017
  // Return ACCEPT_FLAG if SIMLIB_HEADER values pass cuts.
  // Return REJECT_FLAG if SIMLIB_HEADER fails cuts, or REWIND flag is set
  //
  // Nov 28 2019: few checks for SIMLIB model.
  // May 30 2020: increment SIMLIB_HEADER.NFOUND_GENCUTS after all cuts.

  int  ID      = SIMLIB_HEADER.LIBID ;
  int  NOBS    = SIMLIB_HEADER.NOBS ;
  int  IDLOCK  = GENLC.SIMLIB_IDLOCK ;
  int  NSKIP   = INPUTS.NSKIP_SIMLIB ;
  int  NWRAP   = SIMLIB_HEADER.NWRAP ;
  double  zCMB = GENLC.REDSHIFT_CMB ;
  bool  ISMODEL_SIMLIB = ( INDEX_GENMODEL == MODEL_SIMLIB );
  bool  ISMODEL_SALT2  = ( INDEX_GENMODEL == MODEL_SALT2  );
  double z_min     = SIMLIB_HEADER.GENRANGE_REDSHIFT[0];
  double pkmjd_min = SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[0];

  int  iskip, icheck ;
  double *ptrGen ;
  char fnam[] = "keep_SIMLIB_HEADER" ;
  int LTRACE  = 0; // (SIMLIB_HEADER.LIBID == 4 );

  // ----------- BEGIN ---------------

  if(LTRACE) {
    printf(" xxx %s: ----------- LIBID=%d  NOBS_RAW=%d ------------ \n", 
	   fnam, SIMLIB_HEADER.LIBID, NOBS );
  }

  // for SIMLIB model, make sure required REDSHIFT and PEAKMJD
  // are in the header.
  if ( ISMODEL_SIMLIB ) {
    if ( pkmjd_min < 1000.0 ) {
      sprintf(c1err,"Missing 'PEAKMJD:' key in SIMLIB header (ID=%d)",
	      SIMLIB_HEADER.LIBID); 
      sprintf(c2err,"PEAKMJD key is required for SIMLIB model.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
    }
    if ( z_min < 0 ) {
      sprintf(c1err,"Missing  'REDSHIFT:' key in SIMLIB header.");
      sprintf(c2err,"REDSHIFT key needed is required for SIMLIB model.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
    }
  }

  if(LTRACE) {printf(" xxx %s: 0  IDLOCK=%d \n", fnam, IDLOCK);}
  if ( IDLOCK >= 0 && IDLOCK != ID ) 
    { return(REJECT_FLAG); }

  if(LTRACE) {printf(" xxx %s: 1  NSKIP=%d \n", fnam, NSKIP );}
  for ( iskip=0; iskip < NSKIP; iskip++ ) {
    if ( ID == INPUTS.SIMLIB_IDSKIP[iskip] ) 
      { return(REJECT_FLAG); }
  }

  if(LTRACE) {printf(" xxx %s: 2  NWRAP=%d \n", fnam, NWRAP );}
  if ( NWRAP > 0 && (SIMLIB_HEADER.LIBID == SIMLIB_ID_REWIND) )
    { return(REJECT_FLAG); }

  if(LTRACE) {printf(" xxx %s: 3  NOBS=%d \n", fnam, NOBS );}
  if ( NOBS < INPUTS.SIMLIB_MINOBS ) 
    { return(REJECT_FLAG); }

  // check SIMLIB GENRANGES and user-input ranges
  if(LTRACE) {printf(" xxx %s: 4  RA=%f\n", fnam, SIMLIB_HEADER.RA );}
  if( SIMLIB_HEADER.RA < INPUTS.GENRANGE_RA[0] ) { return(REJECT_FLAG) ; }
  if( SIMLIB_HEADER.RA > INPUTS.GENRANGE_RA[1] ) { return(REJECT_FLAG) ; }
  SIMLIB_HEADER.NFOUND_RA++ ;    

  if(LTRACE) {printf(" xxx %s: 5  DEC=%f\n", fnam, SIMLIB_HEADER.DEC );}
  if( SIMLIB_HEADER.DEC < INPUTS.GENRANGE_DEC[0] ) { return(REJECT_FLAG) ; }
  if( SIMLIB_HEADER.DEC > INPUTS.GENRANGE_DEC[1] ) { return(REJECT_FLAG) ; }
  SIMLIB_HEADER.NFOUND_DEC++ ;    


  // Check optional CUTWIN_REDSHIFT so that each LIBID
  // can have its own z-range (originally for WFIRST study)
  // Note that this is a cut, not a range to generate.
  if(LTRACE) {printf(" xxx %s: 6  zCMB=%f\n", fnam, zCMB );}
  if ( zCMB < SIMLIB_HEADER.CUTWIN_REDSHIFT[0] )  { return(REJECT_FLAG); }
  if ( zCMB > SIMLIB_HEADER.CUTWIN_REDSHIFT[1] )  { return(REJECT_FLAG); }


  // check redshift window
  ptrGen = SIMLIB_HEADER.GENRANGE_REDSHIFT ;
  if(LTRACE) {
    printf(" xxx %s: 8 check zrange: %.3f to %.3f\n", 
	   fnam, ptrGen[0], ptrGen[1] );
  }
  if ( ptrGen[0] > 0.0 ) {
    icheck = check_SIMLIB_GENRANGE(INPUTS.GENRANGE_REDSHIFT, ptrGen);
    if ( icheck < 0 ) { return(REJECT_FLAG); }
    if ( ptrGen[1] > ptrGen[0] ) { SIMLIB_HEADER.REGEN_FLAG = 1; }
  }

  // check PEAKMJD
  ptrGen = SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE ;
  if(LTRACE) 
    {printf(" xxx %s: 9 PEAKMJD.RANGE=%.3f,%.3f\n", fnam,ptrGen[0],ptrGen[1]);}
  if ( ptrGen[0] > 1000. ) {
    icheck = check_SIMLIB_GENRANGE(INPUTS.GENRANGE_PEAKMJD, ptrGen);
    if ( icheck < 0 ) { return(REJECT_FLAG); }    
    if ( ptrGen[1] > ptrGen[0] ) 
      { SIMLIB_HEADER.REGEN_FLAG = 1; }
    else { 
      GENLC.PEAKMJD         = ptrGen[0]; 
      GENLC.ISOURCE_PEAKMJD = ISOURCE_PEAKMJD_SIMLIB; 
    }
  }


  if ( ISMODEL_SALT2 ) {
    // check SALT2c
    if(LTRACE) {printf(" xxx %s: 10 check SALT2c \n", fnam );}
    ptrGen = SIMLIB_HEADER.GENGAUSS_SALT2c.RANGE ;
    if ( fabs(ptrGen[0]) < 90.0 ) {
      icheck = check_SIMLIB_GENRANGE(INPUTS.GENGAUSS_SALT2c.RANGE, ptrGen);
      if ( icheck < 0 ) { return(REJECT_FLAG); }    
      SIMLIB_HEADER.REGEN_FLAG = 1; 
    }
    
    // check SALT2x1
    if(LTRACE) {printf(" xxx %s: 11 check SALT2x1\n", fnam );}
    ptrGen = SIMLIB_HEADER.GENGAUSS_SALT2x1.RANGE ;
    if ( fabs(ptrGen[0]) < 90.0 ) {
      icheck = check_SIMLIB_GENRANGE(INPUTS.GENGAUSS_SALT2x1.RANGE, ptrGen);
      if ( icheck < 0 ) { return(REJECT_FLAG); }    
      SIMLIB_HEADER.REGEN_FLAG = 1; 
    }
  } // end ISMODEL_SALT2


  // Nov 26 2017: check NREPEAT for LCLIB
  set_SIMLIB_NREPEAT();
  if(LTRACE) 
    {printf(" xxx %s: 12 SIMLIB_NREPEAT=%d\n", fnam, INPUTS.SIMLIB_NREPEAT );}
  if ( INPUTS.SIMLIB_NREPEAT == 0 ) { return(REJECT_FLAG); }
  

  if(LTRACE) {printf(" xxx %s: 99 END\n", fnam ); debugexit(fnam); }
  
  // if we get here, keep this ID
  SIMLIB_HEADER.NFOUND_GENCUTS++ ;

  return(ACCEPT_FLAG) ;

} // end keep_SIMLIB_HEADER


// =======================================================
void SIMLIB_addCadence_SPECTROGRAPH(void) {

  // Created Aug 23 2017
  // For each SPECTROGRAPH key in the SIMLIB, which specifies
  // MJD and TEXPOSE, add additional SIMLIB_OBS_RAW observation
  // for each synthetic band as if that band had been included
  // in the SIMLIB file.
  // The original SPECTROGRAPH entry is masked out by setting
  // its OPTLINE = -9, but then appended to the end of the 
  // SIMLIB_OBS_RAW list along with each synthetic band.
  //
  // 9.18/2017: set PIXSIZE

  int NOBS       = SIMLIB_OBS_RAW.NOBS ;
  int NOBS_SPEC  = SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH ;
  int NFILT_SPEC = GENLC.NFILTDEF_SPECTROGRAPH  ; // Numb of synthetic filters
  int ifilt, ispec, OBSRAW, ISTORE, NOBS_ADD, APP ;
  double MJD, TEXPOSE, PIXSIZE ;
  char *FIELD, *TEL ;
  //  char fnam[] = "SIMLIB_addCadence_SPECTROGRAPH" ;

  // --------------- BEGIN --------------

  // always mask off original SPECTROGRAPH def in case we bail below.
  for ( ispec=0; ispec < NOBS_SPEC; ispec++ ) {
    OBSRAW = SIMLIB_OBS_RAW.OBSLIST_SPECTROGRAPH[ispec];
    SIMLIB_OBS_RAW.OPTLINE[OBSRAW] = -9 ;
  }

  // if SPCTROGRAPH is to be ignored, then bail.
  if ( SPECTROGRAPH_USEFLAG == 0 ) { return ; }
  if ( NOBS_SPEC == 0 )            { return ; }


  if ( NFILT_SPEC == 0 ) { NFILT_SPEC=1; }   // at least the spectrum
  ISTORE = NOBS;  NOBS_ADD=0;

  for ( ispec=0; ispec < NOBS_SPEC; ispec++ ) {

    // get MJD and epxosure time 
    OBSRAW  = SIMLIB_OBS_RAW.OBSLIST_SPECTROGRAPH[ispec];
    MJD     = SIMLIB_OBS_RAW.MJD[OBSRAW] ;
    TEXPOSE = SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[OBSRAW] ;
    FIELD   = SIMLIB_OBS_RAW.FIELDNAME[OBSRAW] ;
    APP     = SIMLIB_OBS_RAW.APPEND_PHOTFLAG[OBSRAW] ;
    TEL     = INPUTS_SPECTRO.INSTRUMENT_NAME ;
    PIXSIZE = SIMLIB_HEADER.PIXSIZE ;

    // add each synthetic filter to SIMLIB OBS_RAW as if it
    // were defined in the SIMLIB file.

    for(ifilt=0 ; ifilt < NFILT_SPEC; ifilt++ ) {
      
      SIMLIB_OBS_RAW.INDX_TAKE_SPECTRUM[ISTORE]   = ispec ;
      SIMLIB_OBS_RAW.IFILT_SPECTROGRAPH[ISTORE]   = ifilt ; 
      SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[ISTORE] = TEXPOSE ;      
      SIMLIB_OBS_RAW.OPTLINE[ISTORE]    = OPTLINE_SIMLIB_SPECTROGRAPH ;
      SIMLIB_OBS_RAW.MJD[ISTORE]        = MJD ;
      SIMLIB_OBS_RAW.PIXSIZE[ISTORE]    = PIXSIZE ;
      
      SIMLIB_OBS_RAW.BAND[ISTORE][0] = 0 ;
      sprintf(SIMLIB_OBS_RAW.FIELDNAME[ISTORE], "%s", FIELD );
      sprintf(SIMLIB_OBS_RAW.TELESCOPE[ISTORE], "%s", TEL);

      if ( APP > 0 ) { SIMLIB_HEADER.NOBS_APPEND++ ; }
      SIMLIB_OBS_RAW.APPEND_PHOTFLAG[ISTORE] = APP ;

      // the quantities below will be evaluated later in 
      // store_SIMLIB_SPECTROGRAPH (as part of 'prepCadence')
      SIMLIB_OBS_RAW.IDEXPT[ISTORE]     = -9 ;
      SIMLIB_OBS_RAW.NEXPOSE[ISTORE]    =  0 ;
      SIMLIB_OBS_RAW.CCDGAIN[ISTORE]    = -9.0 ;
      SIMLIB_OBS_RAW.READNOISE[ISTORE]  = -9.0 ;
      SIMLIB_OBS_RAW.SKYSIG[ISTORE]     = -9.0 ;
      SIMLIB_OBS_RAW.PSFSIG1[ISTORE]    = -9.0 ;
      SIMLIB_OBS_RAW.PSFSIG2[ISTORE]    = -9.0 ;
      SIMLIB_OBS_RAW.PSFRATIO[ISTORE]   = -9.0 ;
      SIMLIB_OBS_RAW.NEA[ISTORE]        = -9.0 ;
      SIMLIB_OBS_RAW.ZPTADU[ISTORE]     = -9.0 ;
      SIMLIB_OBS_RAW.ZPTERR[ISTORE]     = -9.0 ;
      SIMLIB_OBS_RAW.MAG[ISTORE]        = 99.0 ;
      ISTORE++ ;  NOBS_ADD++ ;

    } // end ifilt loop
  } // end i Loop over SPECTROGRAPH keys

  
  // -------------------------
  // update NOBS
  SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH = NOBS_ADD;

  // make sure that NOBS is deterministic in case we neeed
  // to recompute this later.
  SIMLIB_OBS_RAW.NOBS = 
    SIMLIB_OBS_RAW.NOBS_READ + SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH ;

  SIMLIB_HEADER.NOBS = SIMLIB_OBS_RAW.NOBS ;

  return ;

} // end SIMLIB_addCadence_SPECTROGRAPH 


// ==============================================
void  SIMLIB_TAKE_SPECTRUM(void) {
  
  // add spectrograph entries w.r.t peakMJD
  // These OBS are not in the SIMLIB, so here we add
  // them as if they were in the SIMLIB.
  //
  // Aug 23 2017: remove NLINE_MJD argument for refactor
  // Mar 02 2019: remove LEGACY flag, and fix bug setting VAL_STORE 
  //              when NFILT==0.
  // Jun 28 2019: update to work with HOST spectrum
  // Apr 01 2021: set GENSPEC.SKIP
  // Apr 26 2021: abort if final NOBS > MXOBS_SIMLIB
  // Apr 28 2021: TEXPOSE = -1 if outside TREST range

  int i, OPT, ifilt, OBSRAW, NOBS_ADD, OPT_FRAME, IS_HOST, IS_TREST ;
  int NSPEC = NPEREVT_TAKE_SPECTRUM ;
  int NFILT = GENLC.NFILTDEF_SPECTROGRAPH ;  // all spectroscopic filters

  float Trest_min, Trest_max, Trest_pad ;
  double EPOCH[2], MJD_REF;
  double z, TOBS, TREST ;
  double TEXPOSE=0.0, VAL_STORE[8] ;
  char fnam[] = "SIMLIB_TAKE_SPECTRUM" ;
  // -------------- BEGIN -------------

  if ( NSPEC == 0 ) { return ; }

  // in case of repeated cadence, reset NOBS to leave out 
  // previous TAKE_SPECTRUM obs.
  SIMLIB_OBS_RAW.NOBS = SIMLIB_HEADER.NOBS =
    (SIMLIB_OBS_RAW.NOBS_READ + SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH) ;
  

  OBSRAW = SIMLIB_OBS_RAW.NOBS ;  NOBS_ADD=0;
 
  z = GENLC.REDSHIFT_CMB ;

  Trest_pad = GENRANGE_TOBS_PAD/(1.0+z) ;
  Trest_min = INPUTS.GENRANGE_TREST[0] - Trest_pad;
  Trest_max = INPUTS.GENRANGE_TREST[1] + Trest_pad ;

  for(i=0; i < NSPEC; i++ ) {  

    OPT_FRAME = INPUTS.TAKE_SPECTRUM[i].OPT_FRAME_EPOCH ;
    IS_HOST   = (OPT_FRAME == GENFRAME_HOST);
    EPOCH[0] = INPUTS.TAKE_SPECTRUM[i].EPOCH_RANGE[0] ;
    EPOCH[1] = INPUTS.TAKE_SPECTRUM[i].EPOCH_RANGE[1] ;

    MJD_REF =  GENSPEC_PICKMJD( 0, i, z, &TOBS, &TREST) ; // for MJD sorting
    GENSPEC.SKIP[i] = false;

    IS_TREST = ( TREST > Trest_min && TREST < Trest_max );

    /* xxx
    if ( SIMLIB_HEADER.LIBID == 22 ) {
      printf(" xxx %s: i=%d of %d  TOBS=%.1f \n", fnam, i, NSPEC, TOBS);
    }
    xxxx */

    // Now get exposure time
    OPT  = INPUTS.TAKE_SPECTRUM[i].OPT_TEXPOSE ;
    if ( OPT == 1 ) {
      GENPOLY_DEF *GENZPOLY = &INPUTS.TAKE_SPECTRUM[i].GENZPOLY_TEXPOSE;
      TEXPOSE = eval_GENPOLY(z,GENZPOLY,fnam);
    }
    else if ( OPT == 2 ) {
      // use max TEXPOSE for now; later correct for user SNR
      TEXPOSE = INPUTS_SPECTRO.TEXPOSE_MAX ;
    }
    else {
      sprintf(c1err,"Invalid OPT_TEXPOSE[%d] = %d", i, OPT);
      sprintf(c2err,"Should be 1(TEXPOSE) or 2(SNR)" );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
    }

    /*
    printf(" xxx i=%d : load MJD=%.3f  TEXPOSE=%.0f  EP=%.2f  z=%.3f \n",
	   i, MJD, TEXPOSE, EPOCH_RAN, z ); fflush(stdout);
    */

    // skip SN spectra outside Trest range of model
    if ( !IS_HOST && !IS_TREST )
      { GENSPEC.SKIP[i] = true ; TEXPOSE = -1.0 ; }

    //    printf(" xxx %s: SKIP[%d] = %d  MJD_REF=%.3f \n", 
    //	   fnam, i, GENSPEC.SKIP[i], MJD_REF );

    if ( NFILT==0 ) {  // no synthetic filters; just spectra
      VAL_STORE[0] = MJD_REF   ; // same MJD_REF each event for MJD-sorting
      VAL_STORE[1] = TEXPOSE ;
      VAL_STORE[2] = (double)i ;   
      store_GENSPEC(VAL_STORE) ; 
    }

    // check synthetic filters
    for(ifilt=0; ifilt < NFILT; ifilt++ ) {  

      if ( OBSRAW < MXOBS_SIMLIB ) {
	SIMLIB_OBS_RAW.OPTLINE[OBSRAW]  = OPTLINE_SIMLIB_SPECTROGRAPH ; 
	SIMLIB_OBS_RAW.MJD[OBSRAW]      = MJD_REF ;
	SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[OBSRAW]  = TEXPOSE ;
	SIMLIB_OBS_RAW.IFILT_SPECTROGRAPH[OBSRAW]    = ifilt;
	SIMLIB_OBS_RAW.INDX_TAKE_SPECTRUM[OBSRAW]    = i ;
	SIMLIB_OBS_RAW.PIXSIZE[OBSRAW] = SIMLIB_HEADER.PIXSIZE ;

	// For field, give survey name (since we don;t know)
	// TO DO: determine FIELD based on MJD range
	sprintf(SIMLIB_OBS_RAW.FIELDNAME[OBSRAW],"%s", 
		SIMLIB_GLOBAL_HEADER.SURVEY_NAME );
	
	sprintf(SIMLIB_OBS_RAW.TELESCOPE[OBSRAW], "%s", 
		INPUTS_SPECTRO.INSTRUMENT_NAME );
      }

      OBSRAW++ ;   NOBS_ADD++ ;	

    } // end ifilt

  } // end i-loop over TAKE_SPECTRUM


  SIMLIB_OBS_RAW.NOBS += NOBS_ADD ;
  SIMLIB_HEADER.NOBS  += NOBS_ADD ;

  if ( SIMLIB_HEADER.NOBS >= MXOBS_SIMLIB ) {
    sprintf(c1err,"NOBS=%d exceeds MXOBS_SIMLIB=%d",
	    SIMLIB_HEADER.NOBS, MXOBS_SIMLIB);
    sprintf(c2err,"Check number of synthetic bands from SPECTROGRAPH");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  return;

} // end  SIMLIB_TAKE_SPECTRUM


// ===================================
void  SIMLIB_prepCadence(int REPEAT_CADENCE) {

  // Created Aug 2017
  //
  // After calling SIMLIB_readNext,
  // apply selection and MJD-sorting on SIMLIB_OBS_RAW entries 
  // and fill SIMLIB_OBS_GEN to use for generation.
  // Also set some GENLC arrays for current LIBID/Cadence.
  //
  // REPEAT_CADENCE=0 --> this is a new cadence
  // REPEAT_CADENCE=1 --> same cadence as before, but new event
  //
  // Mar 11 2018: fix PIXSIZE for FWHM_ARCSEC units
  // Mar 15 2018: MAG -> MAG_UNDEFINED
  // Jul 11 2018: fix bug by setting PIXSIZE outside UNIT if-block
  // Dec 16 2019: if MAG = MAG_ZEROFLUX, then skip undefined mag-check.
  // Jan 23 2020: check for FUDGE_ZPTERR
  // Feb 22 2020: fetch SCALE_SKYSIG_T for template noise scale.
  // Dec 08 2020: increase max PSF limit fro 9.9 to 20 (for LSST)
  // Feb 28 2021: check PSF_UNIT for NEA 

  int NOBS_RAW    = SIMLIB_OBS_RAW.NOBS; // xxx SIMLIB_HEADER.NOBS ;
  int NEW_CADENCE = (REPEAT_CADENCE == 0 ) ;
  int ISTORE,  OPTLINE, OBSRAW ;
  double PIXSIZE, FUDGE_ZPTERR, NEA, PSF[3], TREST ;
  double z1       = 1.0 + GENLC.REDSHIFT_CMB ;
  bool IS_SPECTRO, BAD_MJD;
  char *UNIT, *BAND ;  
  char fnam[] = "SIMLIB_prepCadence" ;

  // --------------- BEGIN ----------------

  GENLC.NEPOCH = 0 ; // Mar 17 2015: needed for LIBID repeats

  if ( NOBS_RAW == 0 ) { return; } // 9.03.2021

  // Apr 13, 2016: 
  // check if anything needs to be re-generated based on header info
  if ( regen_SIMLIB_GENRANGES() < 0 ) { return ; }

  // transfer some SIMLIB_HEADER info to GENLC struct
  GENLC.RA         = SIMLIB_HEADER.RA ;
  GENLC.DEC        = SIMLIB_HEADER.DEC ;

  if ( INPUTS.OPT_MWEBV == OPT_MWEBV_FILE ) 
    { GENLC.MWEBV = SIMLIB_HEADER.MWEBV ; }
  GENLC.SIMLIB_ID  = SIMLIB_HEADER.LIBID ;
  sprintf(GENLC.FIELDNAME[0], "%s", SIMLIB_HEADER.FIELD);
  sprintf(GENLC.TELESCOPE[0], "%s", SIMLIB_HEADER.TELESCOPE);

  // --------------------------------------------------------------
  // do a few things for a NEW cadence, 
  // but not for a repeated/re-used cadence
  //   1) prepare duplicate MJDs for sorting
  //   2) sanity checks on values
  //   3) check change of units for PSF and SKYSIG

  if ( NEW_CADENCE ) { 

    for(ISTORE=0; ISTORE < NOBS_RAW ; ISTORE++ ) {

      // 1a. prep duplicate MJDs
      SIMLIB_prepMJD_forSORT(ISTORE);

      OPTLINE  = SIMLIB_OBS_RAW.OPTLINE[ISTORE] ;
      if ( OPTLINE != OPTLINE_SIMLIB_S ) { continue ; }

      // 2. sanity checks to catch crazy [nan] values. Second arg is NVAL=1
      checkval_D("ZPTAVG", 1, &SIMLIB_OBS_RAW.ZPTADU[ISTORE],  6.0, 50.0) ;
      checkval_D("ZPTERR", 1, &SIMLIB_OBS_RAW.ZPTERR[ISTORE],  0.0,  5.0 ) ;
      if ( !SIMLIB_GLOBAL_HEADER.NEA_PSF_UNIT ) {
	checkval_D("PSF1",   1, &SIMLIB_OBS_RAW.PSFSIG1[ISTORE], 0.0, 30.0 ) ;
	checkval_D("PSF2",   1, &SIMLIB_OBS_RAW.PSFSIG2[ISTORE], 0.0, 30.0 ) ;
	checkval_D("PSFrat", 1, &SIMLIB_OBS_RAW.PSFRATIO[ISTORE],0.0,  1.0 ) ;
      }
      checkval_D("SKYSIG", 1, &SIMLIB_OBS_RAW.SKYSIG[ISTORE],  0.0,  1.0E5);

      PIXSIZE = SIMLIB_OBS_RAW.PIXSIZE[ISTORE] ; 

      // 3a. unit check for optional units of PSF and SKYSIG
      UNIT = SIMLIB_GLOBAL_HEADER.PSF_UNIT ;
      if ( strcmp(UNIT,SIMLIB_PSF_ARCSEC_FWHM ) == 0 ) {
	// convert FWHM(arcsec) back to Sigma(pixels)
	SIMLIB_OBS_RAW.PSFSIG1[ISTORE] /= (PIXSIZE * FWHM_SIGMA_RATIO);
	SIMLIB_OBS_RAW.PSFSIG2[ISTORE] /= (PIXSIZE * FWHM_SIGMA_RATIO);
      }
      else if ( SIMLIB_GLOBAL_HEADER.NEA_PSF_UNIT ) {
	// convert NEA back to Sigma(pixels)
	NEA = SIMLIB_OBS_RAW.NEA[ISTORE];
	if ( strcmp(UNIT,SIMLIB_PSF_NEA_ARCSECSQ ) == 0 )
	  { NEA /= (PIXSIZE * PIXSIZE);  SIMLIB_OBS_RAW.NEA[ISTORE]=NEA; } 
	SIMLIB_OBS_RAW.PSFSIG1[ISTORE]  = sqrt(NEA/(2.0*TWOPI)) ;
	SIMLIB_OBS_RAW.PSFSIG2[ISTORE]  = 0.0;
	SIMLIB_OBS_RAW.PSFRATIO[ISTORE] = 0.0;
      }

      // if PSF is NOT already NEA, compute NEA (Feb 2021)
      if ( !SIMLIB_GLOBAL_HEADER.NEA_PSF_UNIT ) {
	PSF[0] = SIMLIB_OBS_RAW.PSFSIG1[ISTORE];
	PSF[1] = SIMLIB_OBS_RAW.PSFSIG2[ISTORE];
	PSF[2] = SIMLIB_OBS_RAW.PSFRATIO[ISTORE];
	NEA    = NoiseEquivAperture(PSF[0], PSF[1], PSF[2] );
	SIMLIB_OBS_RAW.NEA[ISTORE] = NEA ; // pixels
      }
      
      // 3b. check for optional units of SKYSIG
      UNIT = SIMLIB_GLOBAL_HEADER.SKYSIG_UNIT ;
      if ( strcmp(UNIT,SIMLIB_SKYSIG_SQASEC ) == 0 ) {
	// convert SKYSIG per sqrt(asec^2) back into pixel
	SIMLIB_OBS_RAW.SKYSIG[ISTORE] *= PIXSIZE ;
      }
      
    } // end ISTORE loop over everything read from cadence

    // -------------------------------------------
    // 1b. sort MJDs for new LIBID only 

    SIMLIB_sortbyMJD();

  } // end NEW_CADENCE 

  // - - - - - - - - - - - - - - - - - - - 
  // chop MJD range into seasons to allow user options
  store_SIMLIB_SEASONS();

  // - - - - - - -
  int isort, ifilt, IFILT_OBS, NEXPOSE, KEEP, NEP, NEP_NEWMJD ;
  int  IFLAG_SYNFILT, IFLAG_TEMPLATE, IFIELD, APP ;
  double MJD, CCDGAIN, RDNOISE, SKYSIG, ZPT[2], MAG ;
  double SKYSIG_T, RDNOISE_T, ZPT_T ;
  double SHIFT_ZPT, SCALE_SKYSIG, SCALE_SKYSIG_T, SCALE_RDNOISE, SCALE_PSF ;
  double MJD_DIF, MJD_LAST_KEEP, DT, DUMMY_STORE[3] ;
  char   *TEL, *FIELD, cfilt[2];

  set_SIMLIB_MJDrange(2,GENLC.MJD_RANGE); // 9.03.2021 refac

  // init stuff before loop over MJDs
  NEP=NEP_NEWMJD=0;  MJD_LAST_KEEP=-9.0;  

  // transfer OBS_RAW to OBS_GEN; latter has cuts and is sorted by MJD

  for ( isort = 0; isort < NOBS_RAW; isort++ ) {

    OBSRAW   = SIMLIB_LIST_forSORT.INDEX_SORT[isort]; 
    OPTLINE  = SIMLIB_OBS_RAW.OPTLINE[OBSRAW] ;

    if ( OPTLINE < 0 ) { continue; }

    IS_SPECTRO = (OPTLINE == OPTLINE_SIMLIB_SPECTROGRAPH );
    if ( IS_SPECTRO ) {
      ifilt = SIMLIB_OBS_RAW.IFILT_SPECTROGRAPH[OBSRAW]; //sparse synth index
      DUMMY_STORE[0] = DUMMY_STORE[1] = -9.0 ;
      store_SIMLIB_SPECTROGRAPH(ifilt, DUMMY_STORE, OBSRAW) ;
    }

    // strip off SIMLIB variables for this sorted MJD
    MJD        = SIMLIB_OBS_RAW.MJD[OBSRAW] ;
    CCDGAIN    = SIMLIB_OBS_RAW.CCDGAIN[OBSRAW] ;
    RDNOISE    = SIMLIB_OBS_RAW.READNOISE[OBSRAW] ;
    SKYSIG     = SIMLIB_OBS_RAW.SKYSIG[OBSRAW] ;
    PSF[0]     = SIMLIB_OBS_RAW.PSFSIG1[OBSRAW] ;
    PSF[1]     = SIMLIB_OBS_RAW.PSFSIG2[OBSRAW] ;
    PSF[2]     = SIMLIB_OBS_RAW.PSFRATIO[OBSRAW] ;
    NEA        = SIMLIB_OBS_RAW.NEA[OBSRAW] ;
    ZPT[0]     = SIMLIB_OBS_RAW.ZPTADU[OBSRAW] ;
    ZPT[1]     = SIMLIB_OBS_RAW.ZPTERR[OBSRAW] ;
    PIXSIZE    = SIMLIB_OBS_RAW.PIXSIZE[OBSRAW];
    MAG        = SIMLIB_OBS_RAW.MAG[OBSRAW] ;
    IFILT_OBS  = SIMLIB_OBS_RAW.IFILT_OBS[OBSRAW] ;
    NEXPOSE    = SIMLIB_OBS_RAW.NEXPOSE[OBSRAW] ;
    BAND       = SIMLIB_OBS_RAW.BAND[OBSRAW];
    FIELD      = SIMLIB_OBS_RAW.FIELDNAME[OBSRAW];
    TEL        = SIMLIB_OBS_RAW.TELESCOPE[OBSRAW];   
    APP        = SIMLIB_OBS_RAW.APPEND_PHOTFLAG[OBSRAW];   

    RDNOISE_T  = SIMLIB_OBS_RAW.TEMPLATE_READNOISE[OBSRAW] ;
    SKYSIG_T   = SIMLIB_OBS_RAW.TEMPLATE_SKYSIG[OBSRAW] ;
    ZPT_T      = SIMLIB_OBS_RAW.TEMPLATE_ZPT[OBSRAW] ;

    // idiot check on MJD might catch invalid SIMLIB entries
    BAD_MJD = ( MJD < 10000.0 || MJD > 2.0E5 );
    if ( BAD_MJD && !IS_SPECTRO ) {
      sprintf(c1err,"Invalid MJD[%d]=%f for LIBID=%d", 
	      OBSRAW, MJD, SIMLIB_HEADER.LIBID);

      int last = OBSRAW-1 ;
      sprintf(c2err,"Previous MJD / band = %f / %s", 
	      SIMLIB_OBS_RAW.MJD[last], SIMLIB_OBS_RAW.BAND[last] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    // Jan 2020: check for ZPTERR fudge from FUDGE_ZPTERR key
    FUDGE_ZPTERR = INPUTS.FUDGE_ZPTERR_FILTER[IFILT_OBS];
    if ( FUDGE_ZPTERR > 0.0000001 ) { ZPT[1] = FUDGE_ZPTERR; }

    // compute a few things from OBS_RAW
    if ( INPUTS.SMEARFLAG_ZEROPT == 0 ) { ZPT[1] = 0.0 ; }

    if ( MAG != MAG_ZEROFLUX ) 
      { if ( MAG < 5.0 || MAG > 50 ) { MAG = MAG_UNDEFINED ; } }
    if ( MAG < MAG_ZEROFLUX    ) { SIMLIB_HEADER.NOBS_SIM_MAGOBS++ ; }
    sprintf(cfilt,  "%c", FILTERSTRING[IFILT_OBS] ) ;

    ifilt     = GENLC.IFILTINVMAP_SIMLIB[IFILT_OBS] ;
    if ( ifilt < 0 || ifilt >= GENLC.NFILTDEF_SIMLIB ) 
      {  ABORT_SIMLIB_FILTER(isort);  }
    //      {  ABORT_SIMLIB_FILTER(OPTLINE,MJD,cfilt);  }

    // check if this MJD is kept.
    KEEP = keep_SIMLIB_OBS(isort);

    if ( KEEP == 0 ) { continue ; }

    // mark this filter as 'used'
    GENLC.SIMLIB_USEFILT_ENTRY[IFILT_OBS] = 1;

    // get optional fudge scales
    get_SIMLIB_SCALES(IFILT_OBS, &SHIFT_ZPT, &SCALE_SKYSIG, &SCALE_SKYSIG_T, 
		      &SCALE_RDNOISE);
    SCALE_PSF = INPUTS.FUDGESCALE_PSF ;

    NEP++; 
    GENLC.NEPOCH = NEP ;

    if ( NEP >= MXEPSIM ) {
      sprintf(c1err, "NEPOCH = %d exceeds array bound at LIBID=%d", 
	      NEP, GENLC.SIMLIB_ID );
      sprintf(c2err, "Either increase MXEPSIM or reduce library size.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    // determine epoch for trigger; either new epoch, or same as before.
    MJD_DIF = MJD - MJD_LAST_KEEP ;
    if ( fabs(MJD_DIF) > INPUTS.NEWMJD_DIF ) {
      NEP_NEWMJD++; 
      GENLC.NEWMJD = NEP_NEWMJD ;
      GENLC.EPOCH_RANGE_NEWMJD[NEP_NEWMJD][0] = NEP ;
    }
    GENLC.EPOCH_RANGE_NEWMJD[NEP_NEWMJD][1] = NEP ;

    
    // store a few things in GENLC struct
    GENLC.IFILT_OBS[NEP]     = IFILT_OBS;
    GENLC.genmag_obs[NEP]    = MAG ;
    GENLC.MJD[NEP]           = MJD ;
    sprintf( GENLC.FIELDNAME[NEP], "%s", FIELD );
    sprintf( GENLC.TELESCOPE[NEP], "%s", TEL   );

    // store min/max TREST (Apr 2021) for efficiency studies.
    // If event is accepted, TRESTMIN/MAX are overwritten in gen_cutwin(),
    // so these values here are for rejected events.
    TREST = (GENLC.MJD[NEP] - GENLC.PEAKMJD )/z1 ;
    if ( TREST < GENLC.TRESTMIN ) { GENLC.TRESTMIN = TREST; }
    if ( TREST > GENLC.TRESTMAX ) { GENLC.TRESTMAX = TREST; }

    // store SIMLIB_GEN quanties vs. NEP that include fudges                  
    SIMLIB_OBS_GEN.MJD[NEP]         = MJD ;
    SIMLIB_OBS_GEN.IFILT_OBS[NEP]   = IFILT_OBS ;
    SIMLIB_OBS_GEN.PIXSIZE[NEP]     = PIXSIZE ;
    SIMLIB_OBS_GEN.CCDGAIN[NEP]     = CCDGAIN ;
    SIMLIB_OBS_GEN.READNOISE[NEP]   = RDNOISE * SCALE_RDNOISE ;
    SIMLIB_OBS_GEN.SKYSIG[NEP]      = SKYSIG * SCALE_SKYSIG ;
    SIMLIB_OBS_GEN.PSFSIG1[NEP]     = PSF[0] * SCALE_PSF ;
    SIMLIB_OBS_GEN.PSFSIG2[NEP]     = PSF[1] * SCALE_PSF ;
    SIMLIB_OBS_GEN.PSFRATIO[NEP]    = PSF[2] ;    // ratio         
    SIMLIB_OBS_GEN.NEA[NEP]         = NEA * (SCALE_PSF*SCALE_PSF); // Feb 2021
    SIMLIB_OBS_GEN.ZPTADU[NEP]      = ZPT[0] + SHIFT_ZPT ;
    SIMLIB_OBS_GEN.ZPTERR[NEP]      = ZPT[1] ;
    SIMLIB_OBS_GEN.MAG[NEP]         = MAG ;
    SIMLIB_OBS_GEN.NEXPOSE[NEP]     = NEXPOSE ;
    SIMLIB_OBS_GEN.ISTORE_RAW[NEP]  = OBSRAW ;
    SIMLIB_OBS_GEN.APPEND_PHOTFLAG[NEP] = APP ;
    sprintf(SIMLIB_OBS_GEN.TELESCOPE[NEP], "%s", TEL   );
    sprintf(SIMLIB_OBS_GEN.FIELDNAME[NEP], "%s", FIELD );

    SIMLIB_OBS_GEN.TEXPOSE_SPECTROGRAPH[NEP] =
      SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[OBSRAW] ;
    SIMLIB_OBS_GEN.INDX_TAKE_SPECTRUM[NEP] =
      SIMLIB_OBS_RAW.INDX_TAKE_SPECTRUM[OBSRAW] ;

    // - --  store correlated template info ... - - - -                       
    IFLAG_SYNFILT  = GENLC.IFLAG_SYNFILT_SPECTROGRAPH[IFILT_OBS] ;
    IFLAG_TEMPLATE = SIMLIB_TEMPLATE.USEFLAG ;

    // check template noise for synthetic filters from spectrograph           
    if ( (IFLAG_TEMPLATE & 2)>0 && IFLAG_SYNFILT==1 ) {
      SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[NEP]    = SKYSIG_T ;
      SIMLIB_OBS_GEN.TEMPLATE_READNOISE[NEP] = RDNOISE_T ;
      SIMLIB_OBS_GEN.TEMPLATE_ZPT[NEP]       = ZPT_T ;
    }


    // check template noise for normal filters   
    if ( (IFLAG_TEMPLATE & 1)>0 && IFLAG_SYNFILT == 0 ) {

      IFIELD = IFIELD_OVP_SIMLIB(1,FIELD);
      if ( IFIELD < 0 || IFIELD >= MXFIELD_OVP ) {
	sprintf(c1err,"Invalid IFIELD=%d for template FIELD='%s'",
		IFIELD,FIELD);
	sprintf(c2err, "Check LIBID=%d  %s-band  OBSRAW=%d of %d", 
		GENLC.SIMLIB_ID,cfilt, OBSRAW, NOBS_RAW );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;
      }  

      SKYSIG_T = SIMLIB_TEMPLATE.SKYSIG[IFIELD][IFILT_OBS] ;
      SKYSIG_T *= (SCALE_SKYSIG * SCALE_SKYSIG_T) ; // Feb 2020
      SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[NEP]    = SKYSIG_T ;	
      SIMLIB_OBS_RAW.TEMPLATE_SKYSIG[OBSRAW] = SKYSIG_T ;
      
      RDNOISE_T  = SIMLIB_TEMPLATE.READNOISE[IFIELD][IFILT_OBS] ;
      RDNOISE_T *= (SCALE_RDNOISE * SCALE_SKYSIG_T) ;
      SIMLIB_OBS_GEN.TEMPLATE_READNOISE[NEP]    = RDNOISE_T ;	
      SIMLIB_OBS_RAW.TEMPLATE_READNOISE[OBSRAW] = RDNOISE_T ;

      ZPT_T = SIMLIB_TEMPLATE.ZPT[IFIELD][IFILT_OBS] + SHIFT_ZPT ;
      SIMLIB_OBS_GEN.TEMPLATE_ZPT[NEP]    = ZPT_T ;	
      SIMLIB_OBS_RAW.TEMPLATE_ZPT[OBSRAW] = ZPT_T ;

    } // end template 

    MJD_LAST_KEEP = MJD ;
  
    // keep track of MJD nearest peak
    DT = MJD - GENLC.PEAKMJD ; 
    if ( fabs(DT) < fabs(GENLC.DTPEAK_MIN) ) 
      { GENLC.IEPOCH_NEARPEAK = NEP ; GENLC.DTPEAK_MIN   = DT ;  }


  } // end isort loop (REFACTOR)


  return ;

} // end SIMLIB_prepCadence


// ====================================================	
void store_SIMLIB_SPECTROGRAPH(int ifilt, double *VAL_STORE, int ISTORE) {

  // Read info following SPECTROGRAPH key in SIMLIB,
  // and store in SIMLIB_OBS_RAW struct.
  //
  // Inputs:
  //  ifilt     --> sparse spectrograph-filter index
  //  VAL_STORE --> MJD,TEXPOSE,INDX_TAKE_SPECTRUM  if MJD > 0
  //  ISTORE    --> storage index for SIMLIB_OBS_RAW
  //
  // When ifilt=0, read rest of SIMLIB line
  //
  //
  // Aug 22 2017: remove *fp argument and use fp_SIMLIB for LEGACY only.
  //              Change is made to work with refactored SIMLIB routines.
  //
  //
  // Aug 23 2017: as part of refactor move lots of local code 
  //              into new functino get_SPECTROGRAPH_ZPTPSFSKY
  //

  int    LDMP = ( ifilt < 0 ) ;
  int    ifilt_obs, IFLAG_TEMPLATE, LCUT;
  int    INDX_TAKE_SPECTRUM ;
  double MJD, TEXPOSE_S, TEXPOSE_T ;
  double ZPT, SKYSIG_S,  SKYSIG_T, PSFSIG  ;

  double TSCALE = INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE ;

  char   cfilt[2] ;
  char   fnam[] = "store_SIMLIB_SPECTROGRAPH" ;

  // --------------- BEGIN ---------------

  // strip off MJD and TEXPOSE from 1st time these
  // quantites were read (i.e., for ISTORE at ifilt=0)
  if ( ifilt == 0 ) {

    if ( VAL_STORE[0] < 0.0 ) {

      MJD       = SIMLIB_OBS_RAW.MJD[ISTORE] ;
      TEXPOSE_S = SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[ISTORE] ;
      INDX_TAKE_SPECTRUM  = SIMLIB_OBS_RAW.INDX_TAKE_SPECTRUM[ISTORE];

      VAL_STORE[0] = MJD ;
      VAL_STORE[1] = TEXPOSE_S ;
      VAL_STORE[2] = (double)INDX_TAKE_SPECTRUM ;
    }
    else {
      // use passed arguments from TAKE_SPECTRUM
      MJD                     = VAL_STORE[0] ;
      TEXPOSE_S               = VAL_STORE[1] ;
      INDX_TAKE_SPECTRUM      = (int)VAL_STORE[2] ;
    }

    SIMLIB_OBS_RAW.MJD[ISTORE]                  = MJD ;
    SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[ISTORE] = TEXPOSE_S * TSCALE ;
    SIMLIB_OBS_RAW.INDX_TAKE_SPECTRUM[ISTORE]   = INDX_TAKE_SPECTRUM ; 

    store_GENSPEC(VAL_STORE) ; // Mar 20 2017

    SIMLIB_OBS_RAW.MJD[ISTORE]  = VAL_STORE[0]; // in case MJD changes


  } // end ifilt=0

  
  // note that ISTORE-ifilt is the SPECTROGRAPH index with ifilt=0
  MJD       = SIMLIB_OBS_RAW.MJD[ISTORE-ifilt] ;
  TEXPOSE_S = SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[ISTORE-ifilt] ; // search
  INDX_TAKE_SPECTRUM = SIMLIB_OBS_RAW.INDX_TAKE_SPECTRUM[ISTORE-ifilt] ;

  TEXPOSE_T = SIMLIB_TEMPLATE.TEXPOSE_SPECTROGRAPH ;          // template
  IFLAG_TEMPLATE = ( TEXPOSE_T > 0.01 ) ;

  LCUT = 
    ( TEXPOSE_T > INPUTS_SPECTRO.TEXPOSE_MIN ) && 
    ( TEXPOSE_T < INPUTS_SPECTRO.TEXPOSE_MAX ) ;

  if ( IFLAG_TEMPLATE && LCUT==0 ) {
    sprintf(c1err,
	    "\n\t Invalid TEMPLATE_TEXPOSE_SPECTROGRAPH = %.1f in SIMLIB",
            TEXPOSE_T);
    sprintf(c2err,"Valid range is %.2f to %.2f sec",
            INPUTS_SPECTRO.TEXPOSE_MIN, INPUTS_SPECTRO.TEXPOSE_MAX);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
  }


  // load successive filters 
  SIMLIB_OBS_RAW.MJD[ISTORE]                  = MJD ;
  SIMLIB_OBS_RAW.TEXPOSE_SPECTROGRAPH[ISTORE] = TEXPOSE_S ;
  SIMLIB_OBS_RAW.INDX_TAKE_SPECTRUM[ISTORE]   = INDX_TAKE_SPECTRUM ;

  if ( LDMP) {
    printf(" xxx ------------- %s DUMP ------------------- \n", fnam); 
    printf(" xxx ifilt=%d  ISTORE=%d  MJD=%.3f\n", 
	   ifilt, ISTORE, MJD ); 
    fflush(stdout);
  }

  // bail if no synthetic filters.
  if ( GENLC.NFILTDEF_SPECTROGRAPH == 0 )
    { SIMLIB_OBS_RAW.IFILT_OBS[ISTORE]  = 0;  return ; }

  ifilt_obs = GENLC.IFILTDEF_SPECTROGRAPH[ifilt] ;
  sprintf(cfilt,  "%c", FILTERSTRING[ifilt_obs] ) ;

  // negative TEXPOSE is for outside TREST range of SN model
  if ( TEXPOSE_S < 0.0 ) 
    { SIMLIB_OBS_RAW.IFILT_OBS[ISTORE] = ifilt_obs;  return ; }

  if ( GENLC.IFLAG_SYNFILT_SPECTROGRAPH[ifilt_obs] == 0 ) {
    sprintf(c1err,"Invalid %s-filter is not a SYN-FILTER ?", cfilt );  
    sprintf(c2err,"Check SYN-FILTERs in kcor input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  get_SPECTROGRAPH_ZPTPSFSKY(ISTORE, ifilt, TEXPOSE_S, TEXPOSE_T,   // inputs
			     &ZPT, &PSFSIG, &SKYSIG_S, &SKYSIG_T ); // returned


 STORE_SIMLIB_RAW:

  SIMLIB_OBS_RAW.IFILT_OBS[ISTORE]  = ifilt_obs ;
  sprintf(SIMLIB_OBS_RAW.BAND[ISTORE],"%s", cfilt );

  SIMLIB_OBS_RAW.CCDGAIN[ISTORE]    = 1.0 ; // unity GAIN 
  SIMLIB_OBS_RAW.READNOISE[ISTORE]  = 0.0 ;
  SIMLIB_OBS_RAW.SKYSIG[ISTORE]     = SKYSIG_S ; // noise per pixel or arcsec^2
  SIMLIB_OBS_RAW.PSFSIG1[ISTORE]    = PSFSIG ;   // sigma1
  SIMLIB_OBS_RAW.PSFSIG2[ISTORE]    = 0.0 ;  // sigma2
  SIMLIB_OBS_RAW.PSFRATIO[ISTORE]   = 0.0 ;  // ratio of two Gaussians
  SIMLIB_OBS_RAW.NEA[ISTORE]        = NoiseEquivAperture(PSFSIG, 0.0, 0.0);
  SIMLIB_OBS_RAW.ZPTADU[ISTORE]     = ZPT ;  // TEXPOSE zero point (pe)
  SIMLIB_OBS_RAW.ZPTERR[ISTORE]     = 0.0 ;  // no error; included in SKYSIG
  SIMLIB_OBS_RAW.MAG[ISTORE]        = 99.0; 

  // set optional template noise for passbands made from spectrograph
  if ( IFLAG_TEMPLATE ) {
    SIMLIB_TEMPLATE.USEFLAG |= 2;
    SIMLIB_OBS_RAW.TEMPLATE_SKYSIG[ISTORE]    = SKYSIG_T ; 
    SIMLIB_OBS_RAW.TEMPLATE_READNOISE[ISTORE] = 0.0 ;
    SIMLIB_OBS_RAW.TEMPLATE_ZPT[ISTORE]       = ZPT ;
    
    /*
    printf(" xxx ifilt_obs=%2d : SKYSIG(S/T) = %.2f/%.2f = %.3f "
	   "TEXPOSE(S,T)=%.1f,%.1f \n", 
	   ifilt_obs, SKYSIG_S, SKYSIG_T, SKYSIG_S/SKYSIG_T,
	   TEXPOSE_S, TEXPOSE_T );    */
  }

  return ;

} // end store_SIMLIB_SPECTROGRAPH

// =================================================
void get_SPECTROGRAPH_ZPTPSFSKY(int OBSRAW, int ifilt,
				double TEXPOSE_S, double TEXPOSE_T,
				double *ZPT, double *PSFSIG, 
				double *SKYSIG_S, double *SKYSIG_T) {

  // Created Aug 12 2017
  // For input exposure times (TEXPOSE_S,TEXPOSE_T), 
  // compute and return
  //  + ZPT       (zero point for search)
  //  + PSFSIG    (sigma of PSF)
  //  + SKYSIG_S  (search SKY sigma)
  //  + SKYSIG_T  (template SKY sigma)
  //
  // Mar 11 2018: fix PIXSIZE to depend on OBSRAW

  // get min/max wavelentgh of syn filter
  double LAMMIN  = INPUTS_SPECTRO.SYN_FILTERLIST_LAMMIN[ifilt] ;
  double LAMMAX  = INPUTS_SPECTRO.SYN_FILTERLIST_LAMMAX[ifilt] ;
  double MJD     = SIMLIB_OBS_RAW.MJD[OBSRAW]; // for debug
  double PIXSIZE = SIMLIB_OBS_RAW.PIXSIZE[OBSRAW]; 
  int ifilt_obs  = GENLC.IFILTDEF_SPECTROGRAPH[ifilt] ;
  int NBT        = INPUTS_SPECTRO.NBIN_TEXPOSE ;
  int IFLAG_TEMPLATE = ( TEXPOSE_T > 0.01 ) ;
  int UNIT_IS_ARCSEC = 
    ( strcmp(SIMLIB_GLOBAL_HEADER.SKYSIG_UNIT,SIMLIB_SKYSIG_SQASEC )==0) ;

  int l,  NBLAM ;
  int OPT_INTERP = 1; // linear
  int LDMP = 0 ;
  double LCEN, LBINMIN, LBINMAX, ZSUM, ZPLAM ;
  double SQSKY_S, SQSKY_T, SQSUM_S, SQSUM_T ;
  double PSFSIG_PIXEL, PSFSIG_ARCSEC, PSF_ARCSECFWHM, NEA ;

  char fnam[] = "get_SPECTROGRAPH_ZPTPSFSKY" ;
  char STR_ZPLAM[]   = "get_ZPTPSFSKY(ZPLAM)";
  char STR_SQSIG_S[] = "get_ZPTPSFSKY(SQSIG_S)";
  char STR_SQSIG_T[] = "get_ZPTPSFSKY(SQSIG_T)";

  // ---------------- BEGIN -------------

  // init output args
  *ZPT = *SKYSIG_S = *SKYSIG_T = *PSFSIG = 0.0 ;

  // bail if filter is not requested by user.

  if ( GENLC.IFILTINVMAP_OBS[ifilt_obs] < 0 ) { return ; }


  ZSUM = SQSUM_S = SQSUM_T = 0.0 ;   NBLAM=0 ;
  for(l=0; l < INPUTS_SPECTRO.NBIN_LAM; l++ ) {

    LBINMIN = INPUTS_SPECTRO.LAMMIN_LIST[l] ;  // of spectral bin
    LBINMAX = INPUTS_SPECTRO.LAMMAX_LIST[l] ;    
    LCEN    = 0.5 * (LBINMIN + LBINMAX) ; // bin center

    // keep bin if its center is withing the user-defined
    // wavelength range in the kcor-input file.
    if ( LCEN < LAMMIN ) { continue ; }
    if ( LCEN > LAMMAX ) { continue ; }

    // interpolate ZPLAM(pe) for this exposure time, at this lambda
    ZPLAM = interp_1DFUN(OPT_INTERP, TEXPOSE_S, NBT,
			 INPUTS_SPECTRO.TEXPOSE_LIST, 
			 INPUTS_SPECTRO.ZP[l], STR_ZPLAM);
    ZSUM    += pow(TEN, 0.4*ZPLAM) ;

    // interpolate SQSIGSKY
    SQSKY_S = interp_1DFUN(OPT_INTERP, TEXPOSE_S, NBT,
			   INPUTS_SPECTRO.TEXPOSE_LIST, 
			   INPUTS_SPECTRO.SQSIGSKY[l], STR_SQSIG_S);
    SQSUM_S += SQSKY_S ; 

    if ( IFLAG_TEMPLATE ) {
      SQSKY_T = interp_1DFUN(OPT_INTERP, TEXPOSE_T, NBT,
			     INPUTS_SPECTRO.TEXPOSE_LIST, 
			     INPUTS_SPECTRO.SQSIGSKY[l], STR_SQSIG_T);
      SQSUM_T += SQSKY_T ; 
    }

    NBLAM++ ; 

    if ( LDMP ==5 ) { 
      double L0 = INPUTS.LAMMIN_OBS[ifilt_obs];
      double L1 = INPUTS.LAMMAX_OBS[ifilt_obs];
      printf(" xxx LCEN=%.1f (%.1f:%.1f) ZPLAM=%.3f \n", 
	     LCEN, L0, L1, ZPLAM );
    }
  }  // end l loop over lambda bins


  if ( ZSUM == 0.0 ) {
    print_preAbort_banner(fnam);
    printf("   for SYN_FILTER='%c' (ifilt_obs=%d)\n",
	   FILTERSTRING[ifilt_obs], ifilt_obs ) ; 
    printf("   LAMRANGE = %.1f to %.1f  at  MJD=%.3f \n",
	    INPUTS.LAMMIN_OBS[ifilt_obs], 
	    INPUTS.LAMMAX_OBS[ifilt_obs], MJD );

    sprintf(c1err,"No spectro-bin center overlaps SYN_FILTER='%c'",
	    FILTERSTRING[ifilt_obs] );
    sprintf(c2err,"--> FLUX=0.  Try wider SYN_FILTER.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  // get ZPT of all spectrograph slices summed into one band
  *ZPT = 2.5*log10(ZSUM) ;
 
  // define artificial PSF to get sky noise per pixel 
  PSF_ARCSECFWHM = 1.0 ; // 1'' FWHM   
  PSFSIG_ARCSEC  = PSF_ARCSECFWHM / 2.355 ;
  PSFSIG_PIXEL   = PSFSIG_ARCSEC / PIXSIZE ;
  if ( UNIT_IS_ARCSEC )  { 
    *PSFSIG = PSF_ARCSECFWHM ;   // later gets converted back to SIGMA
    NEA    = 2.0 * TWOPI * (PSFSIG_ARCSEC * PSFSIG_ARCSEC ) ;
  } 
  else  { 
    *PSFSIG = PSFSIG_PIXEL ; 
     NEA    = 2.0 * TWOPI * (PSFSIG_PIXEL  * PSFSIG_PIXEL ) ;
  }


  // get skyNoise per pixel (or per sqArcSec)
  *SKYSIG_S = sqrt(SQSUM_S/NEA) ; // noise per pixel or arcsec^2
  if ( IFLAG_TEMPLATE )  { 
    *SKYSIG_T  = sqrt(SQSUM_T/NEA);  // template sky noise
    *SKYSIG_T *= ( TEXPOSE_S/TEXPOSE_T ); // scale template to search image
  }

  if ( LDMP ) {
    printf(" xxx -------------------------------------- \n");
    printf(" xxx synBand %c : TEXPOSE[S,T] = %.1f , %.1f \n",
	   FILTERSTRING[ifilt_obs], TEXPOSE_S, TEXPOSE_T );
    printf(" xxx synBand %c : %.0f-%.0f  "
	   "MJD=%.3f ZPT=%.3f  PSF=%.2f  SKYSIG=%.2f,%.2f\n",
	   FILTERSTRING[ifilt_obs], 
	   INPUTS.LAMMIN_OBS[ifilt_obs],
	   INPUTS.LAMMAX_OBS[ifilt_obs],
	   MJD, *ZPT, *PSFSIG, *SKYSIG_S, *SKYSIG_T ); fflush(stdout);
  }
  
  return ;

} // end get_SPECTROGRAPH_ZPTPSFSKY


// =================================================
int keep_SIMLIB_OBS(int isort) {

  // Created Aug 2017
  // Return 1 of this raw-OBS is kept.
  // Return 0 to reject this OBS.
  //
  // Inputs:
  //   OBS      = MJD-sort index; needed to get absolute OBS index
  //
  // Sep 24 2017: check field.
  // Sep 02 2021: remove obsolete REPEAT arg

  int  KEEP=1, NOKEEP=0; 
  int  ifilt, ifilt_obs, OBS ;
  int  LTRACE= 0 ; 
  double  MJD, MJDrange[2];
  char *FIELD ;
  char fnam[] = "keep_SIMLIB_OBS" ;

  // ------------ BEGIN --------------

  
  if (LTRACE) {
    printf(" xxx ---------------------------- \n");
    printf(" xxx 0 isort=%d  \n",isort); 
  }
  OBS  = SIMLIB_LIST_forSORT.INDEX_SORT[isort] ; // absolute OBS_RAW index

  FIELD = SIMLIB_OBS_RAW.FIELDNAME[OBS] ;
  MJD   = SIMLIB_OBS_RAW.MJD[OBS] ;

  // compute & store SEASON info; return MJDrange to keep SIMLIB entries
  if ( INPUTS.DEBUG_FLAG == 903 ) { // legacy
    set_SIMLIB_MJDrange(2,GENLC.MJD_RANGE);
  }

  if ( MJD < GENLC.MJD_RANGE[0] ) { return(NOKEEP); }
  if ( MJD > GENLC.MJD_RANGE[1] ) { return(NOKEEP); }

  if (LTRACE) {printf(" xxx 1 \n"); fflush(stdout); }
  // check option to skip MJDs
  if ( INPUTS.SIMLIB_NSKIPMJD > 0 ) {
    float xskip = (float)(INPUTS.SIMLIB_NSKIPMJD+1) ;
    float xline = (float)isort ;
    float xmod  = fmodf(xline,xskip);
    if ( xmod > 0.0 ) { return(NOKEEP); }
  }


  // check option to skip field 
  if ( SKIP_SIMLIB_FIELD(FIELD) ) { return(NOKEEP); }

  // apply PEAKMJD cut-window for SIMLIB_DUMP option (Feb 2013)
  // Skip this check for nominal use because OBS before/after
  // PEAKMJD can be used.
  if (LTRACE) {
    printf(" xxx 2 MJD=%f (CUTWIN_PEAKMJD: %.1f to %.1f)\n",
	   MJD, INPUTS.GENRANGE_PEAKMJD[0],INPUTS.GENRANGE_PEAKMJD[1] );
    fflush(stdout);
  }
  if ( INPUTS.SIMLIB_DUMP >= 0 ) {      
    if ( MJD < INPUTS.GENRANGE_PEAKMJD[0] ) { return(NOKEEP); }
    if ( MJD > INPUTS.GENRANGE_PEAKMJD[1] ) { return(NOKEEP); }
  }

  // always apply MJD cut on MJD window (Aug 2017)
  if (LTRACE) {
    printf(" xxx 3 MJD=%f (CUTWIN_MJD: %.1f to %.1f) \n",
	   MJD, INPUTS.GENRANGE_MJD[0], INPUTS.GENRANGE_MJD[1] );
    fflush(stdout);
  }
  if ( MJD < INPUTS.GENRANGE_MJD[0] )  { return(NOKEEP); }
  if ( MJD > INPUTS.GENRANGE_MJD[1] )  { return(NOKEEP); }
  

  // bail if this filter is not requested by user
  ifilt_obs = SIMLIB_OBS_RAW.IFILT_OBS[OBS] ;
  ifilt     = GENLC.IFILTINVMAP_OBS[ifilt_obs]; 
  if (LTRACE) {
    printf(" xxx 4 ifilt=%d of %d\n",ifilt, GENLC.NFILTDEF_OBS ); 
    fflush(stdout); 
  }
  if ( ifilt < 0                   ) { return(NOKEEP); }
  if ( ifilt >= GENLC.NFILTDEF_OBS ) { return(NOKEEP); }

  // if we get here, KEEP epoch
  return(KEEP);

} // end keep_SIMLIB_OBS

// ==============================================
void set_SIMLIB_MJDrange(int OPT, double *MJDrange) {

  // 
  // Return MJDrange
  //
  // Input: 
  //   OPT = 1 --> called before processing SIMLIB header
  //   OPT = 2 --> called after processing SIMLIB header
  //     (OPT=1 isn't used since processing before SIMLIB header
  //          has too many problems)
  //
  // Output:  MJDrange[0:1] = range of MJD to keep SIMLIB entries
  //
  // Jan  6 2018: zHEL -> zCMB
  // Mar 18 2018:
  //  + NLINE starts at 0 (not 1) to fix bug that was removing
  //    first epoch in season for ENTIRE_SEASON option.
  //
  // Mar 27 2018: check ENTIRE_SEASON & ENTIRE_SURVEY options
  // Aug 06 2018: refactor by moving SEASON computation elsewhere
  // Jul 20 2019: update MJD range for strong lens.
  // Sep 03 2021: remove obsolete sameFlag arg; add OPT arg

  int  MSKOPT = INPUTS.SIMLIB_MSKOPT;
  int  KEEP_ENTIRE_SEASON = 
    (MSKOPT & SIMLIB_MSKOPT_ENTIRE_SEASON );
  int  KEEP_ENTIRE_SURVEY = 
    (MSKOPT & SIMLIB_MSKOPT_ENTIRE_SURVEY );
  int REPEAT_UNTIL_ACCEPT = 
    (MSKOPT & SIMLIB_MSKOPT_REPEAT_UNTIL_ACCEPT );

  double PEAKMJD = GENLC.PEAKMJD ;
  double z       = GENLC.REDSHIFT_CMB ;
  double z1      = 1. + z ;
  double Tpad    = GENRANGE_TOBS_PAD ;
  double Tmin    = PEAKMJD + (z1 * INPUTS.GENRANGE_TREST[0]) - Tpad ;
  double Tmax    = PEAKMJD + (z1 * INPUTS.GENRANGE_TREST[1]) + Tpad ;
  double TMPmin, TMPmax;
  int    ISEASON, LDMP = 0 ; 
  char fnam[] = "set_SIMLIB_MJDrange" ;

  // -------------- BEGIN ------------

  // printf(" xxx %s: Tmin / Tmax = %.3f / %.3f \n", fnam, Tmin, Tmax);

  MJDrange[0] = -9.0 ;
  MJDrange[1] = -9.0 ;

  if ( GENSL.NIMG > 0 ) { Tmin = GENSL.MJDMIN; Tmax=GENSL.MJDMAX;  }

  if ( KEEP_ENTIRE_SEASON && KEEP_ENTIRE_SURVEY ) {
    sprintf(c1err,"Cannot specify ENTIRE_SEASON (%d) and ENTIRE_SURVEY(%d)",
	    SIMLIB_MSKOPT_ENTIRE_SEASON, SIMLIB_MSKOPT_ENTIRE_SURVEY) ;
    sprintf(c2err,"Pick one or neither in SIMLIB_MSKOPT key.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  // - - - - - - - - 
  // if LIBID header hasn't been processed yet, just set MJD range
  // based on GENRANGE_TREST. If KEEP_ENTIRE_SEASON/SURVEY is set,
  // set MJD_RANGE to be wide open so that it gets set properly
  // in a 2nd call later.
  if ( OPT == 1 ) {
    if ( KEEP_ENTIRE_SEASON   || 
	 KEEP_ENTIRE_SURVEY   ||
	 REPEAT_UNTIL_ACCEPT  || 
	 INPUTS.SIMLIB_NREPEAT > 1 ) {
      MJDrange[0] = 1000.0 ;  // not ready to compute, so wide open
      MJDrange[1] = 100000.0 ;
    } else {
      MJDrange[0] = Tmin ;
      MJDrange[1] = Tmax ;
    }
    return ;
  }
  
  // - - - - - - - - 
  if ( INPUTS.SIMLIB_DUMP >= 0 ) {
    MJDrange[0] = INPUTS.GENRANGE_PEAKMJD[0] ;
    MJDrange[1] = INPUTS.GENRANGE_PEAKMJD[1] ;
    return ;
  }

  if ( KEEP_ENTIRE_SEASON ) {
    // option: keep all MJDs in season if Trest range overlaps
    // any part of season. e.g. useful for blind analysis challenge.
    for(ISEASON=0; ISEASON < SIMLIB_HEADER.NSEASON; ISEASON++ ) {
      TMPmin = SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][0] ;
      TMPmax = SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][1] ;
      if ( Tmax > TMPmin ) 
	{ MJDrange[1] = TMPmax; }

      if ( Tmin < TMPmax && MJDrange[0] < 0.0 ) 
	{ MJDrange[0] = TMPmin; MJDrange[1] = TMPmax; }
    } // end ISEASON loop
    
  }
  else if ( KEEP_ENTIRE_SURVEY ) {
    // keep every MJD in the entire survey (Mar 27 2018)
    MJDrange[0] = SIMLIB_HEADER.MJDRANGE_SURVEY[0] ;
    MJDrange[1] = SIMLIB_HEADER.MJDRANGE_SURVEY[1] ;
  }
  else {
    // default: keep MJDs only within user-defined TREST range
    MJDrange[0] = Tmin;
    MJDrange[1] = Tmax;
  }

  // if user Tmin/Tmax range is completely outside SIMLIB MJD-range,
  // return with MJDrange=-9,-9.
  // 
  if ( Tmax < SIMLIB_HEADER.MJDRANGE_SURVEY[0]+Tpad ||
       Tmin > SIMLIB_HEADER.MJDRANGE_SURVEY[1]-Tpad ) 
    { 
      MJDrange[0] = MJDrange[1] = -9.0 ;
      return ; 
    }

  
  // sanity tests
  int ERRFLAG = ( MJDrange[0] < 0.0 || 
		  MJDrange[1] < 0.0 || 
		  MJDrange[0] > MJDrange[1] 
		  );

  // ------------ check dump option ---------------
  LDMP = 0 ;
  if ( LDMP || ERRFLAG ) {
    printf("\n xxx --------- SEASON DUMP for LIBID=%d --------------- \n",
	   GENLC.SIMLIB_ID );
    
    for(ISEASON=0; ISEASON < SIMLIB_HEADER.NSEASON; ISEASON++ ) {
      if ( ISEASON < MXSEASON_SIMLIB ) {
	printf(" xxx ISEASON=%2d : MJDrange = %.3f to %.3f \n",
	       ISEASON, 
	       SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][0],
	       SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][1]  );
      }
    }
    printf(" xxx user Tmin, Tmax = %.1f to %.1f \n", Tmin, Tmax);

    printf(" xxx MJDRANGE_SURVEY = %.3f to %.3f \n",
	   SIMLIB_HEADER.MJDRANGE_SURVEY[0],
	   SIMLIB_HEADER.MJDRANGE_SURVEY[1] );
       
    printf(" xxx KEEP_ENTIRE_[SEASON,SURVEY] = %d, %d \n",
	   KEEP_ENTIRE_SEASON, KEEP_ENTIRE_SURVEY);

    printf(" xxx PEAKMJD=%.3f  z=%.3f \n", PEAKMJD, z);
    printf(" xxx MJDrange(keep) = %.3f to %.3f  (%d days)\n", 
	   MJDrange[0], MJDrange[1], (int)(MJDrange[1]-MJDrange[0]) );
    fflush(stdout);
  } // end LDMP

  if ( ERRFLAG ) {
    sprintf(c1err,"Invalid MJDrange = %.3f to %.3f",
	    MJDrange[0], MJDrange[1] );
    sprintf(c2err,"Check SIMLIB and SEASON logic");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  return ;

} // end set_SIMLIB_MJDrange


// ==============================================
void store_SIMLIB_SEASONS(void) {

  // Compute and store SIMLIB_HEADER.[seasonInfo]
  //
 
  double MJD, MJD_LAST;
  int    isort, NLINE_MJD, NLINE, ISEASON, ISGAP, FIRST ;  
  char fnam[] = "store_SIMLIB_SEASONS" ;

  // -------------- BEGIN ------------

  SIMLIB_HEADER.NSEASON = 1;
  MJD_LAST = -9.0 ;
  NLINE_MJD = SIMLIB_LIST_forSORT.NMJD;

  for ( isort = 0; isort < NLINE_MJD; isort++ ) {

    NLINE = SIMLIB_LIST_forSORT.INDEX_SORT[isort] ;

    // strip off  MJD
    MJD   = SIMLIB_OBS_RAW.MJD[NLINE] ;

    // check first MJD in SIMLIB
    FIRST = ( MJD_LAST < 0.0 ) ;

    if ( FIRST )  { 
      ISEASON=0; 
      SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][0]  = MJD; 
      SIMLIB_HEADER.MJDRANGE_SURVEY[0]           = MJD; 
    }
    SIMLIB_HEADER.MJDRANGE_SURVEY[1]             = MJD; 

    // check for large time-gap which defines new season
    ISGAP = ( MJD - MJD_LAST > TGAP_SEASON_SIMLIB ) ;
    
    if ( FIRST==0 && ISGAP ) { 
      SIMLIB_HEADER.NSEASON++ ; 
      ISEASON = SIMLIB_HEADER.NSEASON-1 ;
      if ( ISEASON < MXSEASON_SIMLIB ) 
	{ SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][0] = MJD; } // MJDmin
    }
    ISEASON = SIMLIB_HEADER.NSEASON-1 ;
    if ( ISEASON < MXSEASON_SIMLIB ) 
      { SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][1] = MJD; } // MJDmax
    SIMLIB_OBS_RAW.ISEASON[NLINE] = ISEASON ;
    MJD_LAST = MJD ;
    
  } // end isort loop


  // compute season lengths
  for(ISEASON=0; ISEASON < SIMLIB_HEADER.NSEASON; ISEASON++ ) {
    SIMLIB_HEADER.TLEN_SEASON[ISEASON] =
      SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][1] -
      SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][0] ;
  }

  if ( SIMLIB_HEADER.NSEASON >= MXSEASON_SIMLIB ) {
    sprintf(c1err,"Invalid NSEASON=%d  exceeds bound of %d",
	    SIMLIB_HEADER.NSEASON, MXSEASON_SIMLIB ) ;
    sprintf(c2err,"LIBID=%d ", GENLC.SIMLIB_ID );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  // Aug 2018: check for minimum  Season length requirement 
  remove_short_SIMLIB_SEASON();

  // check time in season |MJD_season_edge - PEAKMJD|
  // DTSEASON_PEAK > 0 if in season; negative if out of season.
  // This variable is useful for selecting NEVT_TOTAL for efficiencies.

  int j;
  double DT, MJD_MIN, MJD_MAX, DT_MIN=999999.9;
  double DT_TEST[2], z1 = 1.0 + GENLC.REDSHIFT_CMB ;

  //  printf(" xxx %s: ------- LIBID=%d ---------- \n", 
  //	 fnam, SIMLIB_HEADER.LIBID);

  for(ISEASON=0; ISEASON < SIMLIB_HEADER.NSEASON; ISEASON++ ) {
    MJD_MIN = SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][0];
    MJD_MAX = SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][1];
    DT_TEST[0] = GENLC.PEAKMJD - MJD_MIN ;
    DT_TEST[1] = MJD_MAX - GENLC.PEAKMJD ;
    for ( j=0; j < 2; j++ ) {
      DT = DT_TEST[j];
      if ( fabs(DT) < fabs(DT_MIN) )  {  DT_MIN = DT; }
    }
  } // end ISEASON
  GENLC.DTSEASON_PEAK =  DT_MIN ;


  // ------------ check dump option ---------------
  int LDMP = 0 ;
  if ( LDMP ) {
    printf("\n xxx --------- SEASON DUMP for LIBID=%d --------------- \n",
	   GENLC.SIMLIB_ID );
    
    for(ISEASON=0; ISEASON < SIMLIB_HEADER.NSEASON; ISEASON++ ) {
      if ( ISEASON < MXSEASON_SIMLIB ) {
	printf(" xxx ISEASON=%2d : MJDrange = %.3f to %.3f \n",
	       ISEASON, 
	       SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][0],
	       SIMLIB_HEADER.MJDRANGE_SEASON[ISEASON][1]  );
      }
    }
  } // end LDMP


  return ;

} // end store_SIMLIB_SEASONS

// =========================================
void remove_short_SIMLIB_SEASON(void) {

  // Aug 2018
  // If 1st or last season length is < INPUTS.SIMLIB_MINSEASON,
  // remove short SEASON as if it didn't exist.
  //
  // Oct 15 2018
  //  + Fix bug with   SIMLIB_HEADER.NSEASON = NSEASON ;
  //

  int    NSEASON   = SIMLIB_HEADER.NSEASON ;
  double MINSEASON = INPUTS.SIMLIB_MINSEASON ;
  char fnam[] = "remove_short_SIMLIB_SEASON" ;
  int  LDMP=0;
  int  iseason, FLAG_REMOVE[MXSEASON_SIMLIB] ;

  double TLEN ;

  int NSEASON_ORIG = NSEASON ;
  double MJDRANGE_SEASON_ORIG[MXSEASON_SIMLIB][2];
  double TLEN_SEASON_ORIG[MXSEASON_SIMLIB];

  // ---------------- BEGIN -------------

  if ( MINSEASON < 0.001 ) { return ; }

  for(iseason=0; iseason < NSEASON; iseason++ )  { 
    FLAG_REMOVE[iseason] = 0;
    TLEN = SIMLIB_HEADER.TLEN_SEASON[iseason] ;
    if ( TLEN < MINSEASON )  { FLAG_REMOVE[iseason] = 1; }
  }


  // store original seasons in separate array
  for(iseason=0; iseason < NSEASON_ORIG; iseason++ ) {
     MJDRANGE_SEASON_ORIG[iseason][0] = 
       SIMLIB_HEADER.MJDRANGE_SEASON[iseason][0] ;
     MJDRANGE_SEASON_ORIG[iseason][1] = 
       SIMLIB_HEADER.MJDRANGE_SEASON[iseason][1] ;

     TLEN_SEASON_ORIG[iseason] =
       SIMLIB_HEADER.TLEN_SEASON[iseason] ;
  }

  // over-write SIMLIB_HEADER.[seasonInfo]
  NSEASON=0 ;
  for(iseason=0; iseason < NSEASON_ORIG; iseason++ ) {
    if ( FLAG_REMOVE[iseason] ) { continue ; }

    SIMLIB_HEADER.MJDRANGE_SEASON[NSEASON][0] = 
      MJDRANGE_SEASON_ORIG[iseason][0];
    SIMLIB_HEADER.MJDRANGE_SEASON[NSEASON][1] = 
      MJDRANGE_SEASON_ORIG[iseason][1];

    SIMLIB_HEADER.TLEN_SEASON[NSEASON] = TLEN_SEASON_ORIG[iseason] ;
    NSEASON++ ;
  }

  if ( NSEASON > 0 ) {
    SIMLIB_HEADER.MJDRANGE_SURVEY[0] = 
      SIMLIB_HEADER.MJDRANGE_SEASON[0][0] - 0.1;
    SIMLIB_HEADER.MJDRANGE_SURVEY[1] = 
      SIMLIB_HEADER.MJDRANGE_SEASON[NSEASON-1][1] + 0.1;
  }
  else {
    SIMLIB_HEADER.MJDRANGE_SURVEY[0] = -9.0 ;
    SIMLIB_HEADER.MJDRANGE_SURVEY[1] = -9.0 ;
  }

  // -------------------------------------
  if ( LDMP ) {
    char cstatus[12];
    printf("\n xxx ================================== \n");
    printf(" xxx DUMP  %s  for  LIBID=%d  CID=%d\n", 
	   fnam, SIMLIB_HEADER.LIBID, GENLC.CID );
    for(iseason=0; iseason < NSEASON_ORIG; iseason++ ) {
      if( FLAG_REMOVE[iseason] )
	{ sprintf(cstatus,"Removed"); }
      else
	{ sprintf(cstatus,"Kept"); }
      printf(" xxx Original SEASON %d: %.1f to %.1f (%.1f days)  %s \n"
	     ,iseason
	     ,MJDRANGE_SEASON_ORIG[iseason][0]
	     ,MJDRANGE_SEASON_ORIG[iseason][1]
	     ,TLEN_SEASON_ORIG[iseason]
	     ,cstatus );
    } // end iseason

    printf(" xxx NSEASON(final)=%d  (%.3f to %.3f) \n",
	   NSEASON, 
	   SIMLIB_HEADER.MJDRANGE_SURVEY[0],
	   SIMLIB_HEADER.MJDRANGE_SURVEY[1] );
    
    fflush(stdout);
  }

  SIMLIB_HEADER.NSEASON = NSEASON ; 

  return ;

} // end  remove_short_SIMLIB_SEASON


// ==============================================
void SIMLIB_prepMJD_forSORT(int ISTORE) {

  // Created Aug 21 2017
  // Increment SIMLIB_LIST_forSORT to use later for
  // sorting MJD. This routing does NOT do the actual
  // MJD sorting, but prepares for duplicate MJDs to
  // preserve SIMLIB-order and avoid random re-order.
  // For each duplicate MJD, add 0.000001 so that 
  // sorting will preserve order.

  double MJD, MJD_DIF, MJD_LAST;
  //  char fnam[] = "SIMLIB_prepMJD_forSORT" ;

  // --------------BEGIN --------------

  MJD      = SIMLIB_OBS_RAW.MJD[ISTORE] ;
  MJD_LAST = SIMLIB_LIST_forSORT.MJD_LAST ;
  MJD_DIF  = MJD - MJD_LAST ;
  SIMLIB_LIST_forSORT.MJD[ISTORE] = MJD; 
  if ( fabs(MJD_DIF) < 0.0001 ) {
    SIMLIB_LIST_forSORT.MJD[ISTORE] = 
      SIMLIB_LIST_forSORT.MJD[ISTORE-1] + 0.00001 ;  // about 1 second
  }
  SIMLIB_LIST_forSORT.MJD_LAST    = MJD ;

  return ;

} // end  SIMLIB_prepMJD_forSORT


// ============================
void SIMLIB_sortbyMJD(void) {

  // Created Jan 17 2018
  // Sort SIMLIB MJDs. Tricky part is that entries with
  // APPEND_PHOTFLAG>0 are NOT to be sorted, but left at
  // the end of the list. The unsorted MJDs allow re-generating
  // the exact same light curves for sorted MJDs, and then 
  // appending new Light curve data. This can happen when
  // observations are used to define new sets of observations.
  // Initial use: WFIRST-imaging triggers IFC(spectra),
  // which overlaps more imaging that appends observations.

  int NOBS_RAW    = SIMLIB_OBS_RAW.NOBS; // xxx SIMLIB_HEADER.NOBS ;
  int NOBS_APPEND = SIMLIB_HEADER.NOBS_APPEND ;
  int NOBS_SORT   = NOBS_RAW - NOBS_APPEND ;
  int  ORDER_SORT = +1 ; // increasing order
  int  isort;
  char fnam[] = "SIMLIB_sortbyMJD" ;

  // ------------- BEGIN --------------

  /* xxx mark delete Sep 3 2021; 
     NOBS==0 is allowed with MDJRANGE cut while reading SIMLIB
  if ( NOBS_SORT < 1 ) {
    sprintf(c1err,"Invalid NOBS_SORT=%d for LIBID=%d CID=%d", 
	    NOBS_SORT, SIMLIB_HEADER.LIBID, GENLC.CID );
    sprintf(c2err,"NOBS_RAW=%d  NOBS_APPEND=%d", 
	    NOBS_RAW, NOBS_APPEND);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }
  xxxx*/

  if ( NOBS_SORT > 0 ) {
    sortDouble( NOBS_SORT, SIMLIB_LIST_forSORT.MJD, ORDER_SORT, 
		SIMLIB_LIST_forSORT.INDEX_SORT ) ;  // return array
  }

  // if there are unsorted MJDs, add them to the back fo the sort-list
  if ( NOBS_SORT < NOBS_RAW ) {
    for ( isort = NOBS_SORT-1; isort < NOBS_RAW; isort++ ) {
      SIMLIB_LIST_forSORT.INDEX_SORT[isort] = isort ; 
    }
  }

  SIMLIB_LIST_forSORT.NMJD = NOBS_RAW ;
  return ;

} // end of  SIMLIB_sortbyMJD(

// ====================================================
void store_GENSPEC(double *VAL_STORE) {

  // load GENSPEC array
  // Input VAL_STORE :
  //  MJD
  //  TEXPOSE
  //  INDX_TAKE_SPECTRUM  if MJD > 0
  //
  // Aug 23 2017: set VAL_STORE[0] = MJD in case MJD changes
  // Mar 01 2019: remove LEGACY flag option
  // Jun 28 2019: store GENSPEC.IS_HOST

  double MJD     = VAL_STORE[0];
  double TEXPOSE = VAL_STORE[1];
  int    INDX    = (int)VAL_STORE[2];  // for TAKE_SPECTRUM 
  double z       = GENLC.REDSHIFT_CMB;
  double z1      = 1.0 + z;

  double TOBS, TREST ;
  int OPT_TEXPOSE, imjd ;
  char fnam[] = "store_GENSPEC" ;

  // -------------- BEGIN ---------------

  imjd        = INDX ;

  if ( NPEREVT_TAKE_SPECTRUM > 0 ) 
    {  OPT_TEXPOSE = INPUTS.TAKE_SPECTRUM[INDX].OPT_TEXPOSE ; }
  else
    { OPT_TEXPOSE = 0 ; }


  if ( OPT_TEXPOSE == 0 )  { 
    TOBS = MJD - GENLC.PEAKMJD ;  TREST = TOBS/z1;
    GENSPEC.INDEX_TAKE_SPECTRUM[imjd]   = -9 ; 
    GENSPEC.INV_TAKE_SPECTRUM[INDX]     = -9 ;
  }
  else if ( OPT_TEXPOSE > 0 ) {
    MJD = GENSPEC_PICKMJD(1,INDX,z, &TOBS, &TREST) ; // return random MJD
    GENSPEC.INDEX_TAKE_SPECTRUM[imjd] = INDX ;
    GENSPEC.INV_TAKE_SPECTRUM[INDX]   = imjd;
    GENSPEC_LAMOBS_RANGE(INDX,GENSPEC.LAMOBS_SNR_LIST[imjd] ); 
  }

  // adjust input MJD in case MJD changes. (Aug 23, 2017)
  VAL_STORE[0] = MJD; 

  GENSPEC.MJD_LIST[imjd]          = MJD ;
  GENSPEC.TOBS_LIST[imjd]         = TOBS ;
  GENSPEC.TREST_LIST[imjd]        = TREST; 
  GENSPEC.TEXPOSE_LIST[imjd]      = TEXPOSE ;      
  GENSPEC.OPT_TEXPOSE_LIST[imjd]  = OPT_TEXPOSE ;

  if ( MJD > 0.0 ) 
    { GENSPEC.IS_HOST[imjd] = 0; }  // SN spectrum
  else
    { GENSPEC.IS_HOST[imjd] = 1; GENSPEC.IMJD_HOST=imjd; }  // HOST spectrum

  GENSPEC.NMJD_TOT++ ;

  /*
  printf(" xxx %s: NMJD_TOT=%d  MJD=%.3f  TOBS=%.1f  (INDX=%d,OPT=%d)\n",
	 fnam, GENSPEC.NMJD_TOT, MJD, TOBS, INDX,OPT_TEXPOSE ); fflush(stdout);
  */

  return ;

} // end store_GENSPEC

// ====================================================
void init_SIMLIB_HEADER(void) {

  // Created Apr 2016 to initialize SIMLIB_HEADER struct.
  // Called for each event.

  int i;
  //  char fnam[] = "init_SIMLIB_HEADER" ;
  // ------------ BEGIN ------------

  SIMLIB_HEADER.REGEN_FLAG  = 0 ;
  SIMLIB_HEADER.NOBS        = NULLINT ;
  SIMLIB_HEADER.NOBS_APPEND     = 0 ;
  SIMLIB_HEADER.NOBS_SIM_MAGOBS = 0 ; // Nov 2019
  SIMLIB_HEADER.RA          = 9999. ;
  SIMLIB_HEADER.DEC         = 9999. ;
  SIMLIB_HEADER.FAKEID      = NULLINT ;
  SIMLIB_HEADER.GALID       = NULLINT ;

  SIMLIB_HEADER.GENRANGE_REDSHIFT[0] = -9.0 ;
  SIMLIB_HEADER.GENRANGE_REDSHIFT[1] = -9.0 ;

  SIMLIB_HEADER.CUTWIN_REDSHIFT[0] = -9.0 ; // wide open
  SIMLIB_HEADER.CUTWIN_REDSHIFT[1] = 99.0 ;

  
  init_GENGAUSS_ASYM( &SIMLIB_HEADER.GENGAUSS_PEAKMJD, (double)999. ) ;
  init_GENGAUSS_ASYM( &SIMLIB_HEADER.GENGAUSS_SALT2x1, (double)999. ) ;
  init_GENGAUSS_ASYM( &SIMLIB_HEADER.GENGAUSS_SALT2c,  (double)999. ) ;  

  sprintf(SIMLIB_HEADER.FIELD,"%s", FIELD_NONAME );
  SIMLIB_HEADER.NFIELD_OVP = 0 ;
  SIMLIB_HEADER.SUBSURVEY_NAME[0] = 0 ;
  // doesn't work  sprintf(SIMLIB_HEADER.SUBSURVEY_NAME, "NULL" );

  SIMLIB_HEADER.NSEASON = 1 ;
  for(i=0; i < MXSEASON_SIMLIB; i++ ) { SIMLIB_HEADER.USE_SEASON[i] = 0;}


  sprintf(SIMLIB_HEADER.TELESCOPE,  "%s", SIMLIB_GLOBAL_HEADER.TELESCOPE) ;
  SIMLIB_HEADER.PIXSIZE  = SIMLIB_GLOBAL_HEADER.PIXSIZE ;

  if ( INPUTS.USE_SIMLIB_SPECTRA ) { NPEREVT_TAKE_SPECTRUM = 0 ; }

  return ;

} // end init_SIMLIB_HEADER


// ================================================
void parse_SIMLIB_IDplusNEXPOSE(char *inString, int *IDEXPT, int *NEXPOSE) {

  // Created Jan 3 2018
  // If inString has no *   -->  IDEXPT = inString and NEXPOSE=1
  // If inString = ID*NEXP  -->  IDEXPT = ID and NEXPOSE=NEXP

  int  IDTMP, NTMP, NRD;
  char star[] = "*" ;
  char WDLIST[2][20], *ptrWDLIST[2];
  //  char fnam[] = "parse_SIMLIB_IDplusNEXPOSE" ;

  // ----------- BEGIN ------------

  NTMP = 1;  // default

  if ( strchr(inString,'*') == NULL ) {
    // no star, just read IDEXPT
    sscanf(inString , "%d", &IDTMP ); 
  }
  else {
    // found star, read both ID and NEXPOSE
    ptrWDLIST[0] = WDLIST[0] ;
    ptrWDLIST[1] = WDLIST[1] ;
    splitString(inString, star, 3,  &NRD, ptrWDLIST ); 
    sscanf( WDLIST[0] , "%d", &IDTMP ); 
    sscanf( WDLIST[1] , "%d", &NTMP ); 
  }

  // load output arguments
  *IDEXPT  = IDTMP ;
  *NEXPOSE = NTMP ;

  /*
  printf(" xxx %s: inString='%s'  NRD=%d   ID=%d  NEXPOSE=%d \n",
	 fnam, inString, NRD, IDTMP, NTMP);
  debugexit(fnam);
  */

  return ;

} // end parse_SIMLIB_IDplusNEXPOSE

// ================================================
void parse_SIMLIB_GENRANGES(char **WDLIST ) {

  // Created Aug 2017
  // Read optional GENRANGE_XXX keys.
  //
  // Inputs:
  //   **WDLIST : key + args
  //
  // Jan 4 2018: read optional DISTANCE key, and convert to zCMB
  // May 29 2020: check for TAKE_SPECTRUM key(s) in header.
  // May 31 2020: parse SALT2 params only of USE_SIMLIB_SALT2 flag is set
  // Nov 12 2020: check RDFLAG_SPECTRA before reading TAKE_SPECTRUM keys
  // Sep 01 2021: refactor to use WDLIST input instead of fp_SIMLIB input
  //
  //  TO-DO: if there are no GENRANGE keys on first LIBID, then don't
  //     bother checking remaining LIBIDs
  //

  char *KEY  = WDLIST[0];

  bool USE_MODEL_SIMLIB = (INDEX_GENMODEL == MODEL_SIMLIB);
  bool RDFLAG_REDSHIFT  = (INPUTS.USE_SIMLIB_REDSHIFT || USE_MODEL_SIMLIB);
  bool RDFLAG_PEAKMJD   = (INPUTS.USE_SIMLIB_PEAKMJD  || USE_MODEL_SIMLIB);
  bool RDFLAG_DISTANCE  = (INPUTS.USE_SIMLIB_DISTANCE || USE_MODEL_SIMLIB);
  bool RDFLAG_SPECTRA   = (INPUTS.USE_SIMLIB_SPECTRA );
  bool RDFLAG_SALT2     = (INPUTS.USE_SIMLIB_SALT2    || USE_MODEL_SIMLIB);
  bool UPDATE_MJDrange = false;
  int  LTMP ;
  double TMPVAL, TMPRANGE[2], dist, MU ;

  char fnam[] = "parse_SIMLIB_GENRANGES" ;

  // ------------ BEGIN ----------------

  // Nov 2015: check optional CUTWIN_REDSHIFT so that each LIBID
  //           can have its own z-range (originally for WFIRST study)
  //    Note that this is a cut, not a range to generate.
  if (  strcmp(KEY,"CUTWIN_REDSHIFT:") == 0  ||
	strcmp(KEY,"REDSHIFT_RANGE:" ) == 0  ) 
    {
      sscanf(WDLIST[1], "%le", &SIMLIB_HEADER.CUTWIN_REDSHIFT[0]);
      sscanf(WDLIST[2], "%le", &SIMLIB_HEADER.CUTWIN_REDSHIFT[1]);
    }
  
  // -------------------------------------------------
  // ------------ READ GENRANGEs ---------------------
  // -------------------------------------------------


  // check for redshift value or range 
  LTMP = 0 ;
  if ( strcmp(KEY,"REDSHIFT:")==0 && RDFLAG_REDSHIFT ) {
    sscanf(WDLIST[1], "%le", &TMPVAL);
    TMPRANGE[0] = TMPRANGE[1] = TMPVAL;  LTMP=1;
  }

  if ( strcmp(KEY,"DISTANCE:")==0 && RDFLAG_DISTANCE ) {
    sscanf(WDLIST[1], "%le", &dist);
    MU = 5.0*log10(dist/1.0E-5);  // 10pc = 1.0E-5 Mpc
    TMPVAL = zcmb_dLmag_invert(MU, &INPUTS.HzFUN_INFO); // returns zCMB
    TMPRANGE[0] = TMPRANGE[1] = TMPVAL;  LTMP=1;
  }

  if ( strcmp(KEY,"GENRANGE_REDSHIFT:")==0 ) {
    sscanf(WDLIST[1], "%le", &TMPRANGE[0] ) ;
    sscanf(WDLIST[2], "%le", &TMPRANGE[1] ) ;
    LTMP=1;
  }

  if ( LTMP ) {
    SIMLIB_HEADER.GENRANGE_REDSHIFT[0] = TMPRANGE[0] ;
    SIMLIB_HEADER.GENRANGE_REDSHIFT[1] = TMPRANGE[1] ;
    SIMLIB_HEADER.REGEN_FLAG = 1;
  }

  // check for optional PEAKMJD value or range in SIMLIB; 
  LTMP=0;
  if ( strcmp(KEY,"PEAKMJD:")==0  && RDFLAG_PEAKMJD )  { 
    sscanf(WDLIST[1], "%le", &TMPVAL ) ;
    TMPRANGE[0] = TMPRANGE[1] = TMPVAL;  LTMP=1;
  }
  if ( strcmp(KEY,"GENRANGE_PEAKMJD:") == 0   || 
       strcmp(KEY,"GENRANGE_PKMJD:")   == 0  ) {
    sscanf(WDLIST[1], "%le", &TMPRANGE[0] ) ;
    sscanf(WDLIST[2], "%le", &TMPRANGE[1] ) ;
    LTMP=1 ;
  }
  if ( LTMP ) {
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.USE      = true;      
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.PEAK     = 0.5*(TMPRANGE[0]+TMPRANGE[1]);
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[0] = TMPRANGE[0] ;
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[1] = TMPRANGE[1] ;
    SIMLIB_HEADER.REGEN_FLAG = 1;
    UPDATE_MJDrange = true ;
  }

  // check for PEAKMJD sigma
  if ( strcmp(KEY,"GENSIGMA_PEAKMJD:")==0  ) {
    sscanf(WDLIST[1], "%le", &TMPVAL ) ;
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.SIGMA[0] = TMPVAL ;
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.SIGMA[1] = TMPVAL ;
  }

  // - - - - - - - - - - - - - - - - - 
  if ( RDFLAG_SALT2 ) {
    LTMP=0;
    // check for SALT2c range & sigma
    if ( strcmp(KEY,"SALT2c:") == 0 ) {
      sscanf(WDLIST[1], "%le", &TMPVAL ) ;
      TMPRANGE[0] = TMPRANGE[1] = TMPVAL ; LTMP=1;
    }
    else if ( strcmp(KEY,"GENRANGE_SALT2c:") == 0 ) {
      sscanf(WDLIST[1], "%le", &TMPRANGE[0] ) ;
      sscanf(WDLIST[2], "%le", &TMPRANGE[1] ) ;
      LTMP=1 ;
    }
    if ( LTMP ) {
      SIMLIB_HEADER.GENGAUSS_SALT2c.USE      = true ;
      SIMLIB_HEADER.GENGAUSS_SALT2c.PEAK     = 0.5*(TMPRANGE[0]+TMPRANGE[1]);
      SIMLIB_HEADER.GENGAUSS_SALT2c.RANGE[0] = TMPRANGE[0] ;
      SIMLIB_HEADER.GENGAUSS_SALT2c.RANGE[1] = TMPRANGE[1] ;
      SIMLIB_HEADER.REGEN_FLAG = 1;
    }
    
    if ( strcmp(KEY,"GENSIGMA_SALT2c:") == 0 ) {
      sscanf(WDLIST[1], "%le", &TMPVAL ) ;
      SIMLIB_HEADER.GENGAUSS_SALT2c.SIGMA[0] = TMPVAL ;
      SIMLIB_HEADER.GENGAUSS_SALT2c.SIGMA[1] = TMPVAL ;
    } 
    
    // check for SALT2x1 range & sigma
    LTMP=0;
    if ( strcmp(KEY,"SALT2x1:")==0 ) {
      sscanf(WDLIST[1], "%le", &TMPVAL ) ; 
      TMPRANGE[0] = TMPRANGE[1] = TMPVAL ; LTMP=1;  
    }
    else if ( strcmp(KEY,"GENRANGE_SALT2x1:")==0 ) {
      sscanf(WDLIST[1], "%le", &TMPRANGE[0] ) ;   
      sscanf(WDLIST[2], "%le", &TMPRANGE[1] ) ;   
      LTMP=1;
    }
    if ( LTMP ) {
      SIMLIB_HEADER.GENGAUSS_SALT2x1.USE      = true ;
      SIMLIB_HEADER.GENGAUSS_SALT2x1.PEAK     = 0.5*(TMPRANGE[0]+TMPRANGE[1]);
      SIMLIB_HEADER.GENGAUSS_SALT2x1.RANGE[0] = TMPRANGE[0] ;
      SIMLIB_HEADER.GENGAUSS_SALT2x1.RANGE[1] = TMPRANGE[1] ;
      SIMLIB_HEADER.REGEN_FLAG = 1;
    }
    
    if ( strcmp(KEY,"GENSIGMA_SALT2x1:") == 0 ) {
      sscanf(WDLIST[1], "%le", &TMPVAL ) ;    
      SIMLIB_HEADER.GENGAUSS_SALT2x1.SIGMA[0] = TMPVAL ;
      SIMLIB_HEADER.GENGAUSS_SALT2x1.SIGMA[1] = TMPVAL ;
    } 
  }  // end RDFLAG_SALT2


  // - - - - - - - - - 
  // May 29 2020 : check for TAKE_SPECTRUM keys
  if ( RDFLAG_SPECTRA ) {
    if ( strcmp(KEY,"TAKE_SPECTRUM:") == 0 ) {
      parse_input_TAKE_SPECTRUM( WDLIST, KEYSOURCE_FILE, NULL ); 
    }
    else if ( strcmp(KEY,"WARP_SPECTRUM:") == 0 ) {
      sscanf(WDLIST[1], "%s", INPUTS.WARP_SPECTRUM_STRING );
    }

  } // end RDFLAG_SPECTRA

  return ;

} // end parse_SIMLIB_GENRANGES


int check_SIMLIB_GENRANGE(double *GENRANGE_ORIG, double *GENRANGE_NEW) {

  // If GENRANGE_NEW has no overlap with GENRANGE_ORIG, return REJECT;
  // otherwise return ACCEPT.
  // If any part of GENRANGE_NEW is outisde GENRANGE_ORIG,
  // modify GENRANGE_NEW so that it is contained inside GENRANGE_ORIG.
  // Note that GENRANGE_NEW can change, but not GENRANGE_ORIG.

  int REJECT = REJECT_FLAG ;

  // --------------- BEGIN ------------

  if ( GENRANGE_NEW[0] >= GENRANGE_ORIG[1] ) { return(REJECT); }
  if ( GENRANGE_NEW[1] <= GENRANGE_ORIG[0] ) { return(REJECT); }
  
  // now make sure that GENRANGE_NEW is contained inside GENRANGE_ORIG
  if ( GENRANGE_NEW[0] < GENRANGE_ORIG[0] ) 
    { GENRANGE_NEW[0] = GENRANGE_ORIG[0];  }

  if ( GENRANGE_NEW[1] > GENRANGE_ORIG[1] ) 
    { GENRANGE_NEW[1] = GENRANGE_ORIG[1];  }

  return(ACCEPT_FLAG);

} // end check_SIMLIB_GENRANGE


// ============================================
int regen_SIMLIB_GENRANGES(void) {

  // Created Apr 13 2016 by R.Kessler
  // called after reading each LIBID header,
  // to check if anything needs to be re-generated.
  //
  // Functions returns +1 to keep LIBID, -1 to reject.
  //
  // Initial motivation is to re-select PEAKMJD, redshift, 
  // color and stretch for ABC method. 
  //
  // Jan 26 2018: fix bug and require HEADER-PEAKMJD>1000 to use it.

  double tmpVal ;
  int LTRACE = 0 ; 
  int REJECT = REJECT_FLAG ;
  double t0_old, z_old, x1_old, c_old ;
  double t0_new, z_new, x1_new, c_new ;
  char fnam[] = "regen_SIMLIB_GENRANGES" ;

  // ------------ BEGIN -------------

  if(LTRACE) {
    printf(" xxx -----------CID=%d  LIBID=%d -------------- \n",
	   GENLC.CID, SIMLIB_HEADER.LIBID ) ;
    printf(" xxx %s: 0 \n", fnam); 
  }

  // check if there is anything to regenerate
  if ( SIMLIB_HEADER.REGEN_FLAG == 0 ) { return(ACCEPT_FLAG); }

  if ( SIMLIB_HEADER.FAKEID > 0 ) 
    { GENLC.CID = SIMLIB_HEADER.FAKEID ; }


  if ( LTRACE ) {
    t0_old = GENLC.PEAKMJD ;
    z_old  = GENLC.REDSHIFT_CMB ; 
    x1_old = GENLC.SALT2x1 ;  
    c_old  = GENLC.SALT2c ;
  }

  //check redshift update 
  double *zWIN = SIMLIB_HEADER.GENRANGE_REDSHIFT ;
  if(LTRACE) {printf(" xxx %s: 1 zWIN=%f,%f\n", fnam,zWIN[0],zWIN[1] ); }
  if ( zWIN[0] > 1.0E-5 ) {
    tmpVal = genz_hubble( zWIN[0], zWIN[1], &INPUTS.RATEPAR );
    if(LTRACE) {printf(" xxx %s: 1b ZGEN=%f\n", fnam,tmpVal ); }
    if ( tmpVal < INPUTS.GENRANGE_REDSHIFT[0] ) { return(REJECT); }
    if ( tmpVal > INPUTS.GENRANGE_REDSHIFT[1] ) { return(REJECT); }
    GENLC.REDSHIFT_CMB     = tmpVal ;
    gen_distanceMag(GENLC.REDSHIFT_CMB, GENLC.REDSHIFT_CMB, // to fix later
		    &GENLC.DLMU, &GENLC.LENSDMU);
  }

  // check PEAKMJD update
  if(LTRACE) {printf(" xxx %s: 2 PEAKMJD=%f\n", fnam,GENLC.PEAKMJD ); }
  if ( SIMLIB_HEADER.GENGAUSS_PEAKMJD.USE  ) {
    tmpVal         = getRan_GENGAUSS_ASYM(&SIMLIB_HEADER.GENGAUSS_PEAKMJD) ;
    GENLC.PEAKMJD  = tmpVal ;  
    if(LTRACE) 
      { printf(" xxx %s: 2b PEAKMJD(HEADER) -> %f\n", fnam,GENLC.PEAKMJD ); }
    if ( tmpVal < INPUTS.GENRANGE_PEAKMJD[0] ) { return(REJECT); }
    if ( tmpVal > INPUTS.GENRANGE_PEAKMJD[1] ) { return(REJECT); }
  }


  if(LTRACE) {printf(" xxx %s: 3 SALT2x1=%f\n", fnam,GENLC.SALT2x1 ); }
  if ( SIMLIB_HEADER.GENGAUSS_SALT2x1.USE ) {
    tmpVal         = getRan_GENGAUSS_ASYM(&SIMLIB_HEADER.GENGAUSS_SALT2x1);
    GENLC.SALT2x1  = tmpVal ;
    if(LTRACE) 
      { printf(" xxx %s: 3b SALT2x1(HEADER) -> %f\n", fnam,GENLC.SALT2x1); }
    if ( tmpVal < INPUTS.GENGAUSS_SALT2x1.RANGE[0] ) { return(REJECT); }
    if ( tmpVal > INPUTS.GENGAUSS_SALT2x1.RANGE[1] ) { return(REJECT); }
  }

  if(LTRACE) {printf(" xxx %s: 4 SALT2c=%f\n", fnam,GENLC.SALT2c ); }
  if ( SIMLIB_HEADER.GENGAUSS_SALT2c.USE ) {
    tmpVal        = getRan_GENGAUSS_ASYM(&SIMLIB_HEADER.GENGAUSS_SALT2c) ;
    GENLC.SALT2c  = tmpVal ;
    if(LTRACE) 
      { printf(" xxx %s: 4b SALT2c(HEADER) -> %f\n", fnam,GENLC.SALT2c); }
    if ( tmpVal < INPUTS.GENGAUSS_SALT2c.RANGE[0] ) { return(REJECT); }
    if ( tmpVal > INPUTS.GENGAUSS_SALT2c.RANGE[1] ) { return(REJECT); }
  }
  

  if ( LTRACE ) {
    t0_new = GENLC.PEAKMJD ;
    z_new  = GENLC.REDSHIFT_CMB ; 
    x1_new = GENLC.SALT2x1 ;  
    c_new  = GENLC.SALT2c ;

    printf(" xxx t0 = %.2f --> %.2f (%.1f to %.1f)\n", t0_old, t0_new, 
	   SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[0], 
	   SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[1] );
    printf(" xxx z  = %7.4f --> %7.4f \n", z_old,  z_new );
    printf(" xxx x1 = %7.4f --> %7.4f \n", x1_old, x1_new );
    printf(" xxx c  = %7.4f --> %7.4f \n", c_old,  c_new );
    fflush(stdout);
  }

  return(ACCEPT_FLAG);

} //  end regen_SIMLIB_GENRANGES


// ========================
void  ABORT_SIMLIB_FILTER(int isort) {

  // Created Feb 2017
  // Missing filter in SIMLIB, so here give proper abort message
  // depending on regular filter or SPECTROGRAPH filter.

  int NOBS_RAW = SIMLIB_OBS_RAW.NOBS; // xxxx SIMLIB_HEADER.NOBS ;
  int OBSRAW   = SIMLIB_LIST_forSORT.INDEX_SORT[isort]; 
  int OPTLINE  = SIMLIB_OBS_RAW.OPTLINE[OBSRAW] ;
  int IFILT_OBS= SIMLIB_OBS_RAW.IFILT_OBS[OBSRAW] ;
  double MJD   = SIMLIB_OBS_RAW.MJD[OBSRAW] ;
  char cfilt[2];
  char fnam[]  = "ABORT_SIMLIB_FILTER" ;

  // ----------- BEGIN ------------


  sprintf(cfilt,  "%c", FILTERSTRING[IFILT_OBS] ) ;

  print_preAbort_banner(fnam);
  printf("   isort=%d  OBSRAW= %d of %d \n", isort, OBSRAW, NOBS_RAW );
  printf("   LIBID=%d  MJD=%9.3f  OPTLINE=%d\n", 
	 GENLC.SIMLIB_ID, MJD, OPTLINE ) ;


  if ( OPTLINE == OPTLINE_SIMLIB_SPECTROGRAPH ) {
    printf("   SPECTROGRAPH FILTERS: %s \n",
	   INPUTS_SPECTRO.SYN_FILTERLIST_BAND );
    sprintf(c1err, "missing Spectrograph filter='%s' ",  cfilt );
    sprintf(c2err, "after 'FILTERS:'  key in SIMLIB header.");
  }
  else {
    sprintf(c1err, "Unknown SIMLIB filter='%s' ",  cfilt );
    sprintf(c2err, "Expecting filter from : '%s'", 
	    SIMLIB_GLOBAL_HEADER.FILTERS );
  }
  
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 


} // end ABORT_SIMLIB_FILTER

// ================================================
int USE_SAME_SIMLIB_ID(int IFLAG) {

  // Dec 16 2015: 
  // Return 1 to re-use the same SIMLIB_ID as last event.
  // Return 0 otherwise to read another LIBID.
  //
  // IFLAG used for debug:
  // IFLAG=1 if called from init_GENLC
  // IFLAG=2 if called from SIMLIB_READ_DRIVER
  //
  // Apr 17 2017: fmodf -> fmod and XNTOT, XNUPD -> double
  //             (to handle very large NGENTOT_LC)
  //

  int IDLOCK = GENLC.SIMLIB_IDLOCK ;
  int OVP, OLD_SAMEFLAG=0, NEW_SAMEFLAG=0 ;
  int LDMP = 0 ; // (GENLC.SIMLIB_ID >= 0 && IFLAG>0 );
  //  char fnam[] = "USE_SAME_SIMLIB_ID" ;

  // -------------- BEGIN --------------

  if ( IFLAG==1 ) { SIMLIB_HEADER.NREPEAT++ ; }

  // cannot repeast on 1st reading
  if ( SIMLIB_HEADER.NREPEAT<=1 ) { return(0); }
  
  // check option to keep using same LIBID
  if ( IDLOCK > 1 && GENLC.SIMLIB_ID > 0 ) { return 1; }

  if ( SIMLIB_HEADER.NREPEAT <= INPUTS.SIMLIB_NREPEAT ) { 
    NEW_SAMEFLAG=1; 
  }
  else if ( IFLAG == 2 ) { 
    SIMLIB_HEADER.NREPEAT=1 ;    // reset
  }

  if ( LDMP ) {
    printf(" xxx ------------------------------------ \n");
    printf(" xxx IFLAG=%d CID=%6d  LIBID=%3d  NREPEAT=%d  NGENLC_TOT=%d \n",
	   IFLAG, GENLC.CID, GENLC.SIMLIB_ID,
	   SIMLIB_HEADER.NREPEAT, NGENLC_TOT );
    printf(" xxx --> SAMEFLAG[OLD,NEW] = %d, %d   IDLOCK=%d\n",
	   OLD_SAMEFLAG, NEW_SAMEFLAG, IDLOCK );    
    fflush(stdout);
  }

  if ( NEW_SAMEFLAG ) { return 1; }

  // - - - - - - - - - - - - - - - - 
  // check option to re-use same LIBID until an event is accepted
  // (for ABC sims)
  OVP = (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_REPEAT_UNTIL_ACCEPT );
  if ( OVP > 0  &&  GENLC.ACCEPTFLAG_LAST == 0 ) 
    { return(1); }
  // - - - - - - - - - - - - - - - - 

  // if we get here, read next LIBID
  return 0 ;

} // end USE_SAME_SIMLIB_ID


// ==============================================
void set_SIMLIB_NREPEAT(void) {

  // Nov 2017
  // set SIMLIB_NREPEAT for LCLIB model, using b-dependent rate.
  // Do nothing for other models.
  //
  // Jun 28 2018: for SIMLIB_DUMP, set NREPEAT=1

  int    NREPEAT_LO, NREPEAT ;
  double relRate, XNREPEAT, XDIF, RA, DEC ;
  char fnam[] = "set_SIMLIB_NREPEAT" ;

  // ----------------- BEGIN -----------------

  if ( INDEX_GENMODEL != MODEL_LCLIB   ) { return ; }

  RA  = SIMLIB_HEADER.RA ;
  DEC = SIMLIB_HEADER.DEC ;
  slaEqgal (RA, DEC,                    // input
	    &GENLC.GLON, &GENLC.GLAT ); // output

  if ( INPUTS.SIMLIB_DUMP >= 0 ) 
    { SIMLIB_HEADER.NREPEAT = 1; return ; }

  if ( SIMLIB_HEADER.NREPEAT != 1 )  
    { return ; }


  relRate  = GALrate_model(GENLC.GLON, GENLC.GLAT, &INPUTS.RATEPAR);
  relRate /= INPUTS.RATEPAR.RATEMAX ;
  XNREPEAT = (double)INPUTS.SIMLIB_MXREPEAT * relRate;
  
  // use random number to pick the low or high integer of XNREPEAT.
  NREPEAT_LO = (int)XNREPEAT ;
  XDIF       = XNREPEAT - (double)NREPEAT_LO ;

  if ( XDIF < getRan_Flat1(1) ) 
    { NREPEAT = NREPEAT_LO ; }
  else
    { NREPEAT = NREPEAT_LO+1 ; }

  int LDMP = 0 ;

  if ( LDMP ) {
    printf(" xxx ------------------------------------------ \n");
    printf(" xxx CID=%d    RA=%f   DEC=%f   RATEMAX=%le\n",	
	   GENLC.CID, RA, DEC, INPUTS.RATEPAR.RATEMAX ) ;
    printf(" xxx b=%f  relRate=%.3f   XNREPEAT=%.3f  (XDIF=%.3f) NREPEAT=%d\n",
	   GENLC.GLAT, relRate, XNREPEAT, XDIF, NREPEAT );
  }

  INPUTS.SIMLIB_NREPEAT = NREPEAT ;

  return ;

} // end set_SIMLIB_NREPEAT

// ================================================
int SKIP_SIMLIB_FIELD(char *field) {

  // Created Mar 2015.
  // Refactored Feb 2021 to enable optional preScale per FIELD
  //
  // Returns 1 if this field should be skipped.
  // Returns 0 if this field should be kept.
  //

  double preScale;
  int iPS, iTEST ;
  char *FIELDLIST = INPUTS.SIMLIB_FIELDLIST ;
  int   OPT_DICT  = 1 ; // --> do partial match
  char fnam[] = "SKIP_SIMLIB_FIELD" ;
  
  if ( strcmp(FIELDLIST,"ALL") == 0 )  { return 0 ; }

  // fetch prescale for this field
  preScale = get_string_dict(OPT_DICT, field, 
			     &INPUTS.DICT_SIMLIB_FIELDLIST_PRESCALE);
  iPS = (int)preScale ;

  //  printf(" xxx %s: field=%s -> PS = %d \n", fnam, field, iPS);
  // fflush(stdout);

  if ( iPS < 0 ) {
    return 1 ;    // field not specified -> skip
  }

  // if we get here, check prescale.
  // Make sure to pick iTEST that is the same for all epochs.
  // Here we pick LIBID and add NWRAP so that each wrap-around
  // of the SIMLIB picks a different set of LIBIDs.
  // To avoid artifacts from harmonics, do NOT use NGENLC_TOT
  // or NGENLC_WRITE because these are corrleated with iTEST.

  iTEST = (SIMLIB_HEADER.LIBID + SIMLIB_HEADER.NWRAP) ;
  if ( ( iTEST % iPS ) == 0 )  // beware of harmonic effects XXX
    { return 0; }  // keep field
  else
    { return 1 ; }  // reject due to prescale

} // end of SKIP_SIMLIB_FIELD

// ================================================
int  parse_SIMLIB_ZPT(char *cZPT, double *ZPT, 
		      char *cfiltList, int *ifiltList) {

  // Created Dec 2011
  // Parse string *cZPT to return either float ZPT
  // or list of filters to sum and *ZPT = -9.0;
  // Examples:
  //  cZPT = 28.654   => *ZPT = 28.654 and cfiltList = unchanged 
  //  cZPT = +28.654  => *ZPT = 28.654 and cfiltList = unchanged 
  //  cZPT = 2+3+4+5  => *ZPT = -9     and cfiltList = 2345  
  //
  // Returned function value is strlen(cfiltList)
  //
  // Mar 1 2014: float *ZPT -> double

  int NFILT, jplus, LEN, I, i, ifilt, ifilt_obs;
  char *ptrPlus, cfilt[4], ctmp[4] ;

  char fnam[] = "parse_SIMLIB_ZPT";

  // -------------- BEGIN -------------

  // init output args
  NFILT = 0;
  *ZPT = -999.0;

  // check if there is a '+' sign in the cZPT string.
  ptrPlus  = strstr(cZPT,"+");
  jplus    = ptrPlus - cZPT ;  // = 1 few if '+' in there
      
  if ( jplus < 1 || jplus > 100 ) {
    *ZPT = atof(cZPT);
    sprintf(ctmp,"%s ", cfiltList );
    ifiltList[0] = filtindx_( ctmp, strlen(ctmp) ); 
  }
  else  { 

    *ZPT = 19.999 ; // 19.999
    cfiltList[0] = 0 ;  // reset input value
    // strip off the filter chars between the plus signs
    LEN = strlen(cZPT);
    i   = -1 ;

    for ( I = 0 ; I < LEN; I++ ) {

      sprintf(cfilt, "%c", cZPT[I] );
      if ( strcmp(cfilt,"+") == 0 ) { continue; }

      strcat(cfiltList, cfilt); // append filter list
      i++ ;

      // check that this is a valid filter.
      // Note that we need a stupid blank space for filtindx_ function
      sprintf(ctmp,"%s ", cfilt);
      ifilt_obs    = filtindx_( ctmp, strlen(ctmp) ); 
      ifiltList[i] = ifilt_obs ;
      ifilt        = GENLC.IFILTINVMAP_OBS[ifilt_obs];
      if ( ifilt < 0 || ifilt >= GENLC.NFILTDEF_OBS ) {
	sprintf(c1err,"Invalid filter='%s'(%d) in cZPT='%s' ", 
		cfilt, ifilt_obs, cZPT );
	sprintf(c2err,"Check simlib file %s", INPUTS.SIMLIB_FILE );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }

      
    } // end loop over LEN
  }
      
  // get number of filters in list for return function arg.
  NFILT = strlen(cfiltList);

  /*
  printf(" ZZZZZ cZPT(%s) = '%s'  ZPT=%f   NFILT=%d  jplus=%d\n", 
	 cfiltList, cZPT, *ZPT, NFILT, jplus );
  */

  return NFILT ;

} // end of parse_SIMLIB_ZPT


// ====================================
void get_SIMLIB_SCALES( int ifilt_obs
		       ,double *SHIFT_ZPT
		       ,double *SCALE_SKYSIG
		       ,double *SCALE_SKYSIG_T
		       ,double *SCALE_READNOISE
		       ) {

  // Apr 26, 2011
  // return SIMLIB scales based on MSKOPT
  //
  // Feb 2012: init to INPUTS.FUDGESCALE_XXX
  // May 2013: READNOISE is scaled separately from SKYNOISE
  //
  // init to fudged values
  //
  // Jun 19 2017: use filter-dependent ZPT shift.
  // Feb 22 2020: add *SCALE_SKYSIG_T argument
  //
  // --------- BEGIN ----------
  
  *SCALE_SKYSIG    = (double)INPUTS.FUDGESCALE_NOISE_SKY ;
  *SCALE_SKYSIG_T  = (double)INPUTS.FUDGESCALE_NOISE_TEMPLATE ;
  *SCALE_READNOISE = (double)INPUTS.FUDGESCALE_NOISE_READ ;
  *SHIFT_ZPT       = (double)INPUTS.FUDGESHIFT_ZPT_FILTER[ifilt_obs] ;
  

  if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 0) )
    { *SHIFT_ZPT  += GENLC.SHIFT_ZPTSIMLIB[ifilt_obs]; }

  if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 1) )  
    {  *SCALE_SKYSIG   *= GENLC.SCALE_NOISE[ifilt_obs] ; }

  if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 2) )
    { *SCALE_READNOISE  *= GENLC.SCALE_NOISE[ifilt_obs] ; }


} // end of get_SIMLIB_SCALES

// ***********************************
void ENDSIMLIB_check(void) {

  // Created Apr 29, 2009 by R.Kessler
  // called after reading entire simlib to check
  // if anything passes MJD, RA and DECL cuts.
  // If not, abort with message instead of wrapping
  // around in an infinite loop.
  //
  // May 30 2020: abort if SIMLIB_HEADER.NFOUND_GENCUTS == 0
  // Dec 09 2020: leave message for QUIT_REWIND

  int  MSKOPT = INPUTS.SIMLIB_MSKOPT ;
  bool QUIT_NOREWIND = ( (MSKOPT & SIMLIB_MSKOPT_QUIT_NOREWIND)>0 );
  char fnam[] = "ENDSIMLIB_check";

  // ------- BEGIN ---------

  if ( INPUTS.SIMLIB_DUMP >= 0 ) { return ; }

  // check option to quit generating after reading SIMLIB once
  if ( QUIT_NOREWIND && SIMLIB_HEADER.NWRAP == 0 ) { 
    GENLC.STOPGEN_FLAG = 1; 
    printf("\n\t STOP generation after reading SIMLIB one time.\n");
    printf("\t    (see option with SIMLIB_MSKOPT += %d) \n",
	   SIMLIB_MSKOPT_QUIT_NOREWIND);
    printf("\n");
    fflush(stdout);
  }

  // don't do error checking until a few wrap-arounds.
  if ( SIMLIB_HEADER.NWRAP < 5 ) { return ; }

  if ( SIMLIB_HEADER.NFOUND_RA == 0 ) {
    sprintf(c1err,"Could not find SIMLIB RA within");
    sprintf(c2err,"GENRANGE_RA : %f - %f  deg" ,
	    INPUTS.GENRANGE_RA[0], INPUTS.GENRANGE_RA[1]  );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( SIMLIB_HEADER.NFOUND_DEC == 0 ) {
    sprintf(c1err,"Could not find SIMLIB DEC within");
    sprintf(c2err,"GENRANGE_DEC : %f - %f  deg" ,
	    INPUTS.GENRANGE_DEC[0], INPUTS.GENRANGE_DEC[1] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( strcmp(INPUTS.SIMLIB_FIELDLIST,"ALL") != 0 ) {
    if ( SIMLIB_HEADER.NFOUND_FIELD == 0 ) {
      sprintf(c1err,"Could not find SIMLIB FIELD '%s' ", 
	      INPUTS.SIMLIB_FIELDLIST );
      sprintf(c2err,"Check sim-input file and SIMLIB file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  if ( SIMLIB_HEADER.NFOUND_GENCUTS == 0 ) {
    sprintf(c1err,"Could not find SIMLIB HEADER passing GENRANGE_XXX cuts.");
    sprintf(c2err,"Check USE_SIMLIB_XXX flags and header values.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


} // end of ENDSIMLIB_check


// *********************************
void init_zvariation(void) {

  // Created May 11, 2009 by R.Kessler
  //
  // Read/parse ZVARIATION_FILE and load ZVARIATION struct.
  // To allow another z-dependent parameter,
  //   * increment MXPAR_ZVAR in snlc_sim.h
  //   * udpate PARDEF_ZVAR array in snlc_sim.h
  //   * apply get_zvariation() call at appropriate place in snlc_sim.c
  //
  // Feb 22 2013: move hard-wired PARDEF_ZVAR list from snlc_sim.h to here.
  //
  // Mar 12 2014: if ZPOLY keys have already been read from sim-input file,
  //              skip reading separate ZVAR file.
  //
  // Mar 16, 2014: major overhaul to define z-dependence on any
  //               population param; see PREFIX_GENGAUSS.
  //
  // Mar 24, 2014: Fixed auxilliary ZPAR file 
  //     INPUT_ZVARIATION[i].FLAG bug (JLM)
  //
  // Apr 23 2015: allow 'END:' key in ZVAR-file.
  // Oct 27 2017: add SKEWNORMAL[0,1,2]
  //
  // Apr 20 2017:
  //  +  always allow host AV & RV (not just for mlcs)
  //  + add "EXPSIG","EXPTAU" to list of GEN[prefixes]
  //
  // Jun 2 2017: add 2nd-peak params to PREFIX_GENGAUSS
  // Jul 19 2017: add GENMAG_OFF_GLOBAL to update list
  // Jul 24 2017: update error reporting of NZBIN exceeding boung.
  // Mar 21 2020: DJB added more parameters to check for zvar. 

  char *ptrZfile, *ptrparname, *ptrPar, *ptrPoly, fileName_full[MXPATHLEN] ;
  char GENPREFIX[60], c_get[60], method[20], parName[60], cpoly[60] ;

#define NPREFIX_GENGAUSS 12
  char PREFIX_GENGAUSS[NPREFIX_GENGAUSS][20] 
    = 
    { "PEAK", "MEAN", "SKEW[0]", "SKEW[1]", "SIGMA[0]", "SIGMA[1]",
      "PROB2", "PEAK2", "SIGMA2[0]", "SIGMA2[1]" ,
      "EXPSIG", "EXPTAU"
    } ;

  double *ptrzval, *ptrzshift, ZTMP, ZMIN, ZMAX, ZGEN[2], shift ;
  int i, i2, NZ, MATCH, FLAG, ipar, IPAR_USR=-9 ;
  int IZVAR_DEJA, IZVAR_FILE, gzipFlag ;
  FILE *fpz;
  char    fnam[] = "init_zvariation" ;

  // -------------- BEGIN ----------------

  IZVAR_DEJA = ( NPAR_ZVAR_USR > 0 ) ; // already read from sim-input file
  IZVAR_FILE = USE_ZVAR_FILE ;
  ptrZfile   = INPUT_ZVARIATION_FILE ;
  ptrparname = fnam; // avoid compile warning
  FLAG       = -9 ;

  if ( IZVAR_DEJA == 0 && IZVAR_FILE == 0 ) { return ; }

  // allow only one input option
  if ( IZVAR_DEJA && IZVAR_FILE ) {
    sprintf ( c1err,"Cannot specify ZVARIATION_FILE and ZPOLY keys");
    sprintf ( c2err,"in sim-input file. Pick one or the other." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;
  }

  // -------------------------------------------
  // hard-wire list of valid parameter names

  NPAR_ZVAR_TOT = 0 ;

  // allow any GENPREFIX to modify population params.

  for (i=0; i < NPREFIX_GENGAUSS; i++ ) {

    sprintf(GENPREFIX,"GEN%s", PREFIX_GENGAUSS[i] );

    // always allow host AV & RV
    sprintf(parName,"%s_RV", GENPREFIX );
    update_PARDEF_ZVAR( parName );
    
    sprintf(parName,"%s_AV", GENPREFIX );
    update_PARDEF_ZVAR( parName );
    
    if ( INPUTS.NPAR_SIMSED > 0 ) { 
      for ( ipar=0; ipar < INPUTS.NPAR_SIMSED; ipar++ ) {
	ptrPar = INPUTS.PARNAME_SIMSED[ipar] ;
	sprintf(parName,"%s_%s", GENPREFIX, ptrPar);      
	update_PARDEF_ZVAR( parName );
      }
    } 
    else if ( INDEX_GENMODEL == MODEL_SALT2 ) { 

      sprintf(parName,"%s_SALT2x1", GENPREFIX );
      update_PARDEF_ZVAR( parName );

      sprintf(parName,"%s_SALT2c", GENPREFIX );
      update_PARDEF_ZVAR( parName );

      sprintf(parName,"%s_SALT2ALPHA", GENPREFIX );
      update_PARDEF_ZVAR( parName );

      sprintf(parName,"%s_SALT2BETA", GENPREFIX );
      update_PARDEF_ZVAR( parName );
    }
    else {

      sprintf(parName,"%s_DELTA", GENPREFIX );
      update_PARDEF_ZVAR( parName );

      sprintf(parName,"%s_DM15", GENPREFIX );
      update_PARDEF_ZVAR( parName );

      sprintf(parName,"%s_STRETCH", GENPREFIX );
      update_PARDEF_ZVAR( parName );
    }
    
  } // end i loop over PREFIX_GENGAUSS


  // add a few miscellaneous variables.
  if ( INPUTS.NPAR_SIMSED == 0 ) { 
    update_PARDEF_ZVAR( "GENEXPTAU_AV"       );  // Mar 2013
    update_PARDEF_ZVAR( "GENTAU_AV"          );  // new naming March 2020
    update_PARDEF_ZVAR( "GENGAUSIG_AV"       );
    update_PARDEF_ZVAR( "GENSIG_AV"          );
    update_PARDEF_ZVAR( "GENGAUPEAK_AV"      );
    update_PARDEF_ZVAR( "GENPEAK_AV"         );
    update_PARDEF_ZVAR( "GENRATIO_AV0"       );
    update_PARDEF_ZVAR( "GENTAU_EBV_HOST"    );
    update_PARDEF_ZVAR( "GENSIG_EBV_HOST"    );
    update_PARDEF_ZVAR( "GENPEAK_EBV_HOST"   );
    update_PARDEF_ZVAR( "GENRATIO_EBV0_HOST" );

    update_PARDEF_ZVAR( "VSI"                ); // Si velocity for VCR model
    update_PARDEF_ZVAR( "GENMAG_OFF_GLOBAL"  ); // added July 2017
  }

  // ---------------------------------
  // check for PARDEF-dump option
  if ( strcmp(ptrZfile,"0") == 0 ) {
    printf("\n\n");
    print_banner("Valid parameters to simulate Z-dependence:");
    for ( i=0; i < MXPAR_ZVAR; i++ )  { printf(" %s", PARDEF_ZVAR[i] ); }
    printf("\n");    happyend();
  }


  // if ZPOLY keys have already been read from sim-input file,
  // skip reading separate ZVAR file.
  if ( IZVAR_DEJA ) { goto ZVAR_SUMMARY ; }


  // -------------------------------------------
  // open file

  fpz = snana_openTextFile(1, PATH_USER_INPUT, ptrZfile,
			   fileName_full, &gzipFlag );
  
  if ( !fpz ) {
    abort_openTextFile("ZVARIATION_FILE", 
		       PATH_USER_INPUT, ptrZfile, fnam);
  }

  print_banner(" Read Z-DEPENDENCE for SIM PARAMETERS ");

  while( (fscanf(fpz, "%s", c_get)) != EOF) {

    // allow END key to ignore anything after
    if ( strcmp(c_get,"END:") == 0 ) { fclose(fpz) ;  goto ZVAR_SUMMARY ; }
    
    if ( strcmp(c_get,"PARAMETER:") == 0 ) {
      IPAR_USR = NPAR_ZVAR_USR;
      ptrparname = INPUT_ZVARIATION[IPAR_USR].PARNAME ;
      readchar(fpz, ptrparname);
      FLAG = INPUT_ZVARIATION[IPAR_USR].FLAG ;
      NPAR_ZVAR_USR++ ;
    }

    if ( strcmp(c_get,"ZPOLY:") == 0 ) {
      if ( FLAG > 0 ) {
	sprintf(c1err, "Found ZPOLY but ZVARIATION[%d].FLAG=%d already set", 
		IPAR_USR, FLAG);
	sprintf(c2err,"for %s . Check file: %s", ptrparname, ptrZfile );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      readchar(fpz,cpoly);
      if ( strstr(cpoly,",") == NULL ) {
	sprintf(c1err,"Must use comma-separate poly coefficients");
	sprintf(c2err,"for ZPOLY option with param = '%s'", ptrparname);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      parse_GENPOLY(cpoly, "z", &INPUT_ZVARIATION[IPAR_USR].POLY, fnam);
    }

    if ( strcmp(c_get,"ZBIN:") == 0 ) {
      INPUT_ZVARIATION[IPAR_USR].FLAG = FLAG_ZBINS_ZVAR ; //JLM
      NZ = INPUT_ZVARIATION[IPAR_USR].NZBIN;
      INPUT_ZVARIATION[IPAR_USR].NZBIN++ ;  

      if ( NZ < MXZBIN_ZVAR ) {
	ptrzval   = &INPUT_ZVARIATION[IPAR_USR].ZBIN_VALUE[NZ];
	ptrzshift = &INPUT_ZVARIATION[IPAR_USR].ZBIN_PARSHIFT[NZ];
	readdouble(fpz, 1, ptrzval );
	readdouble(fpz, 1, ptrzshift );
	
	ZMIN   = INPUT_ZVARIATION[IPAR_USR].ZBIN_VALUE[0];
	ZMAX   = INPUT_ZVARIATION[IPAR_USR].ZBIN_VALUE[NZ];
      }
    }

  } // end of read loop

  // - - - - -

  NZ = INPUT_ZVARIATION[NPAR_ZVAR_USR-1].NZBIN;
  if ( NZ >= MXZBIN_ZVAR ) {
    sprintf(c1err,"NZBIN=%d for %s exceeds bound MXZBIN_ZVAR=%d", 
	    NZ,ptrparname, MXZBIN_ZVAR );
    sprintf ( c2err,"Check file: '%s' ", ptrZfile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  fclose(fpz) ;

  // -------------------

 ZVAR_SUMMARY:

  // print summary, and make sure that GEN-Z range is covered
  // by ZBIN option

  for ( i=0; i < NPAR_ZVAR_USR; i++ ) {
    ptrparname = INPUT_ZVARIATION[i].PARNAME;
    ptrPoly    = INPUT_ZVARIATION[i].POLY.STRING ;
    FLAG       = INPUT_ZVARIATION[i].FLAG ;
    NZ         = INPUT_ZVARIATION[i].NZBIN;
    ZMIN       = INPUT_ZVARIATION[i].ZBIN_VALUE[0];
    ZMAX       = INPUT_ZVARIATION[i].ZBIN_VALUE[NZ-1];

    // set method comment based on flag

    if ( FLAG == FLAG_ZPOLY_ZVAR ) 
      { sprintf(method,"ZPOLY"); }
    else if ( FLAG == FLAG_ZBINS_ZVAR ) {
      { sprintf(method,"ZBINS"); }
    }
    else {
      sprintf(c1err,"Invalid FLAG=%d for %s", FLAG, ptrparname );
      errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
    }

    printf("   %s for %s(%s)   method=%s\n", 
	   fnam, ptrparname, ptrPoly, method );

    // ========== IDIOT CHECKS ===============

    // make sure that this parameter is valid
    MATCH = 0;
    for ( i2=0; i2 < MXPAR_ZVAR; i2++ ) 
      { if ( strcmp(ptrparname,PARDEF_ZVAR[i2]) == 0 ) MATCH=1; }

    if ( MATCH == 0 ) {
      sprintf(c1err,"Unrecognized SIM-parameter '%s'", ptrparname );
      sprintf(c2err,"Cannot apply Z-dependence.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // do NOT allow extrapolation past ZBIN range

    if ( FLAG == FLAG_ZBINS_ZVAR ) {
      ZGEN[0] = INPUTS.GENRANGE_REDSHIFT[0];
      ZGEN[1] = INPUTS.GENRANGE_REDSHIFT[1];
      if ( ZMIN > ZGEN[0] ||  ZMAX < ZGEN[1] ) {
	sprintf(c1err," ZBIN range %5.3f - %5.3f for %s",
		ZMIN, ZMAX, ptrparname );
	sprintf(c2err,"does not cover ZGEN-range %5.3f - %5.3f", 
		ZGEN[0], ZGEN[1] );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }

    // make sure ZVARIATION is zero at z=0
    ZTMP  = 0.0 ;
    shift = get_zvariation(ZTMP,ptrparname);
    if ( shift != 0.0 ) {
      sprintf(c1err,"%s-shift = %f at z=0", ptrparname, shift);
      sprintf(c2err,"but shift(z=0) must be zero");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  }  // end if i-loop over NPAR_ZVAR_USR

  return ;

} // end of init_zvariation


// *********************************
void update_PARDEF_ZVAR(char *parName) {


  char fnam[] = "update_PARDEF_ZVAR" ;

  // --------------- BEGIN -------------

  NPAR_ZVAR_TOT++ ;

  if ( NPAR_ZVAR_TOT >= MXPAR_ZVAR ) {

    int ipar ;
    print_preAbort_banner(fnam);
    for(ipar=1; ipar < NPAR_ZVAR_TOT; ipar++ ) {
      printf("\t ZVARIATION param %2d: '%s' \n", 
	     ipar, PARDEF_ZVAR[ipar] );
      fflush(stdout);
    }

    sprintf(c1err,"NPAR_ZVAR_TOT = %d exceeds bound of ", NPAR_ZVAR_TOT);
    sprintf(c2err,"MXPAR_ZVAR = %d", MXPAR_ZVAR );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  sprintf(PARDEF_ZVAR[NPAR_ZVAR_TOT], "%s", parName );

  return ;

} // end of  update_PARDEF_ZVAR


// *********************************
void cp_zvariation(char *outFile_zvar) {

  // read & copy Zvariation file to SIM/[VERSION] area
  // Aug 19 2014: pass outFile_zvar

  FILE *fpz, *fpz2 ;
  char *ptrZfile, *ptrZfile2,  cline[84] ;
  char fnam[] = "cp_zvariation" ;

  // ------------ BEGIN ---------------

  if ( USE_ZVAR_FILE == 0 ) { return ; }

  ptrZfile  = INPUT_ZVARIATION_FILE;
  ptrZfile2 = outFile_zvar ;  // target file for copy

  // open files

  if ( (fpz = fopen(ptrZfile, "rt"))==NULL ) {   
    sprintf ( c1err, "Cannot open ZVARIATION_FILE (NPAR_ZVAR_USR=%d)",
	      NPAR_ZVAR_USR );
    sprintf ( c2err," '%s' ", ptrZfile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  if ( (fpz2 = fopen(ptrZfile2, "wt"))==NULL ) {   
    sprintf ( c1err, "Cannot open target ZVARIATION_FILE " );
    sprintf ( c2err," '%s' ", ptrZfile2 );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  while( fgets(cline, 80, fpz)  != NULL ) 
    { fprintf(fpz2, "%s", cline); }
  
  fclose(fpz);
  fclose(fpz2);

  return ;

} // end of cp_zvariation

// *********************************
double get_zvariation(double z, char *parname) {

  // Created May 11, 2009 by R.Kessler
  // Determine additive z-dependent shift for parname.
  // Note that init_zvariation() must be called first.
  //
  // Sep 23, 2011: float -> double

  double shift,  ZMIN, ZMAX, ZTMP, FMIN, FMAX  ;
  int  ipar, ipar_match, MATCH, FLAG, NZ, iz, izmin, izmax ;
  char fnam[] = "get_zvariation" ;

  // ------------ BEGIN -------------

  shift = 0.0 ;

  if ( NPAR_ZVAR_USR == 0 ) { return(shift);  }

  // for first few events, check  that parnmae is on the
  // PARDEF list, even if it's not on the user-list;
  // this is to catch cases where a parname is passed that
  // would never get used.


  MATCH = 0;
  if ( NGENLC_WRITE < 3 ) {
    for ( ipar=0; ipar < MXPAR_ZVAR; ipar++ ) 
      { if ( strcmp(parname,PARDEF_ZVAR[ipar]) == 0 ) MATCH=1; }

    if ( MATCH == 0 ) {
      sprintf(c1err,"Undefined parname = '%s' at z=%f", parname, z );
      sprintf(c2err,"Check file: '%s'", INPUT_ZVARIATION_FILE );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }


  // look for index that matches this IPAR
  ipar_match = -9 ;
  for ( ipar=0; ipar < NPAR_ZVAR_USR; ipar++ ) {
    if ( strcmp(parname,INPUT_ZVARIATION[ipar].PARNAME) == 0  ) 
      { ipar_match = ipar; }
  }

  if ( ipar_match < 0 ) { return(shift) ; }

  // --- Now we have a match; compute z-dependent shift

  ipar = ipar_match ;
  FLAG = INPUT_ZVARIATION[ipar].FLAG;

  if ( FLAG == FLAG_ZPOLY_ZVAR ) {
    shift = eval_GENPOLY(z,&INPUT_ZVARIATION[ipar].POLY,fnam);
  }
  else if ( FLAG == FLAG_ZBINS_ZVAR ) {  // linear interpolate ZBINS

    // find zinbs that bound input redshift
    NZ = INPUT_ZVARIATION[ipar].NZBIN;
    izmin = 0 ;
    for ( iz=0; iz < NZ; iz++ ) {
      ZTMP = INPUT_ZVARIATION[ipar].ZBIN_VALUE[iz];
      if ( ZTMP <= z ) izmin = iz ;
    }
    if ( izmin <= 0 ) {
      sprintf(c1err,"izmin=0 for %s and z=%5.3f", parname, z );
      errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
    }

    izmax = izmin + 1;
    ZMIN = INPUT_ZVARIATION[ipar].ZBIN_VALUE[izmin];
    ZMAX = INPUT_ZVARIATION[ipar].ZBIN_VALUE[izmax];
    FMIN = INPUT_ZVARIATION[ipar].ZBIN_PARSHIFT[izmin];
    FMAX = INPUT_ZVARIATION[ipar].ZBIN_PARSHIFT[izmax];
    shift = FMIN + (FMAX-FMIN)*(z-ZMIN)/(ZMAX-ZMIN) ;
  }
  else {
    sprintf(c1err,"Invalid FLAG=%d  for  %s", FLAG, parname);
    errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
  }

  return(shift) ;

} // end of get_zvariation


// *****************************************
GENGAUSS_ASYM_DEF get_zvariation_GENGAUSS(double z, char *parName, 
					  GENGAUSS_ASYM_DEF *GENGAUSS_INP) {

  // Mar 2014
  // Return GENGAUSS structure with z-variation.
  // This allows applying z-variation to the population parameters
  // instead of simply the generated parameters.
  //
  // Inputs:
  //  * z = redshift
  //  * parName: name of parameter described by asymmetric Gaussian
  //  * GENGAUSS_INP = struct with PEAK, SIGMA, RANGE, SKEW
  //
  // Aug 30 2016: use copy_ASYMMGAUSS 
  // Sep 01 2016: add SKEWNORMAL, and reFactor with for loops
  // Jun 02 2017: check for 2nd Gaussian peak

  int i;
  char PARNAME[60];
  GENGAUSS_ASYM_DEF GENGAUSS_OUT ;

  char fnam[] = "get_zvariation_GENGAUSS" ;

  // ---------- BEGIN ---------

  copy_GENGAUSS_ASYM(GENGAUSS_INP, &GENGAUSS_OUT) ;

  // check for z-variation in the peak
  sprintf(PARNAME, "GENMEAN_%s", parName );
  GENGAUSS_OUT.PEAK += get_zvariation(z,PARNAME) ;

  sprintf(PARNAME, "GENPEAK_%s", parName );
  GENGAUSS_OUT.PEAK += get_zvariation(z,PARNAME) ;

  // check for variation in the skew
  for(i=0; i < 2; i++ ) {
    sprintf(PARNAME, "GENSKEW[%d]_%s", i, parName );
    GENGAUSS_OUT.SKEW[i] += get_zvariation(z,PARNAME) ;
  }

  // check for variation in the SIGMAs
  for(i=0; i < 2; i++ ) {
    sprintf(PARNAME, "GENSIGMA[%d]_%s", i, parName );
    GENGAUSS_OUT.SIGMA[i] += get_zvariation(z,PARNAME) ;
  }


  // idiot checks 
  sprintf(PARNAME, "GENSIGMA_%s", parName );
  checkVal_GENGAUSS(PARNAME, &GENGAUSS_OUT.SIGMA[0], fnam );

  sprintf(PARNAME, "GENRANGE_%s", parName );
  checkVal_GENGAUSS(PARNAME, &GENGAUSS_OUT.RANGE[0], fnam );

  // -------------------------------------------
  // check for 2nd peak (Jun 2017)

  sprintf(PARNAME, "GENPROB2_%s", parName );
  GENGAUSS_OUT.PROB2 += get_zvariation(z,PARNAME) ;

  sprintf(PARNAME, "GENPEAK2_%s", parName );
  GENGAUSS_OUT.PEAK2 += get_zvariation(z,PARNAME) ;

  // check for variation in the SIGMAs
  for(i=0; i < 2; i++ ) {
    sprintf(PARNAME, "GENSIGMA2[%d]_%s", i, parName );
    GENGAUSS_OUT.SIGMA2[i] += get_zvariation(z,PARNAME) ;
  }


  return GENGAUSS_OUT ;

} // end of  get_zvariation_GENGAUSS" ;


// ***********************************
void init_CIDRAN(void) {

  // Oct 21, 2009 : init random CID option.
  //
  // CIDOFF is the offset to select from list
  // so that different samples (with different CIDOFF)
  // are guaranteed to have different CIDs.
  //
  // Jan 9, 2011: Declare logical USED array to quickly
  //              check which CIDs have already been used.
  //              Saves lots of init time for large jobs.
  //
  // Sep 09, 2012: replace local USED[] array with global *USEDCID
  //               and free memory when finished.
  //
  // Dec 4, 2012: allocate memory for CIDRAN_LIST
  //
  // Apr 6 2017: replace USEDCID array with exec_cidmask function
  //             that burns 1 bit per CID instead of 2 bytes per CID.
  //
  // Apr 10 2017: new algorithm to pick CID closest to CIDRAN to 
  //              avoid abort when there are too few CIDRANs 
  //              remaining.
  //
  // Aug 11 2017: after generating CID list, re-init randoms with NEWSEED
  //
  // Jul 22 2018: check minimum CID from user input (see MNCID)
  //
  // Jan 04 2021; 
  //   + remove useless calc of NEWSEED
  //
  // Jul 21 2021: 
  //   + init all CIDRAN_LIST[i] = -9
  //   + fix bug implementing CIDRAN_MIN (see CIDRAN_OFF)
  //   + abort if CIDOFF < CIDRAN_MIN
  //

  int NPICKRAN_ABORT ; // abort after this many tries
  int i, i2, j, NPICKRAN, NSTORE_ALL, NSTORE, CIDRAN, CIDTMP, CIDADD;
  int CIDMAX, CIDMIN, CIDRAN_OFF, USED, NTRY, NCALL_random=0, LDMP=0 ;
  int L_STDOUT, NEWSEED ;

  double r8, XN8, XNUPD8;

  char fnam[] = "init_CIDRAN";

  // ------------ BEGIN ------------

  fflush(stdout);
  if ( WRFLAG_CIDRAN  <= 0 ) { return ; }


  CIDMAX     = INPUTS.CIDRAN_MAX ; // -> local var      
  CIDMIN     = INPUTS.CIDRAN_MIN ; 

  CIDRAN_OFF = INPUTS.CIDOFF - INPUTS.CIDRAN_MIN ;
  NSTORE_ALL = CIDRAN_OFF + INPUTS.NGEN_LC + INPUTS.NGENTOT_LC ;

  NSTORE     = INPUTS.NGEN_LC + INPUTS.NGENTOT_LC ; // for this sim-job

  if ( CIDRAN_OFF < 0 ) {
    sprintf(c1err,"Invalid  CIDOFF=%d < CIDRAN_MIN=%d", 
	    INPUTS.CIDOFF, INPUTS.CIDRAN_MIN );
    sprintf(c2err,"CIDOFF must be >= CIDRAN_MIN");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  if ( NSTORE == 0 ) {
    sprintf(c1err,"NGEN_LC = NGENTOT_LC = 0 ?? ");
    sprintf(c2err,"Cannot initialize CIDRAN");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  sprintf(BANNER,"init_CIDRAN: initialize %d RANDOM CIDs (CIDOFF=%d)",
	  NSTORE_ALL, INPUTS.CIDOFF );

  print_banner( BANNER );
  printf("  ");

  // --------------------------------------
  // error checking on bounds.

  if ( CIDMAX > MXCID_SIM ) {
    sprintf(c1err, "CIDRAN_MAX=%d exceeds limit of MXCID_SIM=%d",
	    CIDMAX, MXCID_SIM );
    sprintf(c2err, "Reduce CIDRAN_MAX to be <= %d", MXCID_SIM );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( NSTORE_ALL >= CIDMAX ) {
    sprintf(c1err,"NSTORE_ALL=%d exceeds bound of CIDRAN_MAX = %d", 
	    NSTORE_ALL, CIDMAX );
    sprintf(c2err,"NSTORE_ALL= %d(CIDOFF) + %d(NGEN_LC) + %d(NGENTOT_LC)",
	    INPUTS.CIDOFF, INPUTS.NGEN_LC, INPUTS.NGENTOT_LC );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  exec_cidmask(0,CIDMAX); // allocate bit-mask memory to store all CIDs

  // allocate memory for big CIDRAN_LIST and print mem usage
  // Beware that allocation could be as much as 400 MB,
  // more than the static size of snlc_sim.exe.
  float XMEM = (float)(NSTORE+2) * (4.0 / 1.0E6) ;
  printf("\t CIDRAN_MIN/MAX = %d/%d   CIDOFF=%d   NSTORE_JOB=%d\n", 
	 CIDMIN, CIDMAX, INPUTS.CIDOFF, NSTORE );
  printf("\t CIDRAN_LIST memory allocation: %6.3f MB \n", XMEM);
  fflush(stdout);

  INPUTS.CIDRAN_LIST = (int*)malloc( (NSTORE+2)*sizeof(int) );
  for(i=0; i < NSTORE+2; i++ ) { INPUTS.CIDRAN_LIST[i] = -9; }

  // start main loop

  if ( NSTORE_ALL < 200000 ) 
    { XNUPD8 = 1.0E4; }     // screen update
  else if ( NSTORE_ALL < 2000000 )   // 2 million
    { XNUPD8 = 1.0E5 ; }
  else if ( NSTORE_ALL < 20000000 )  // 20 million
    { XNUPD8 = 1.0E6 ; }
  else
    { XNUPD8 = 5.0E6 ; }

  // set abort to 5% of NSTORE (Apr 10 1017)
  NPICKRAN_ABORT = (NSTORE_ALL/20);
  if ( NPICKRAN_ABORT < 100 ) { NPICKRAN_ABORT = 100; }

  // - - - - - - - - - - - - - - - - - - 

  for( i = 0; i <= NSTORE_ALL ; i++ ) {

    NPICKRAN = 0 ;
    
  PICKRAN:
    NPICKRAN++ ;

    if ( NPICKRAN >= NPICKRAN_ABORT ) {
      sprintf(c1err,"Could not find unused random CID after %d tries.",
	      NPICKRAN);
      sprintf(c2err,"Something is wrong at istore=%d of %d ", 
	      i, NSTORE_ALL) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // use unix 'random' instead of snana 'rangen()'.
    // because rangen() uses a finite list of randoms.
    r8 = unix_getRan_Flat1(0);   
    NCALL_random++ ;
    CIDRAN = CIDMIN + (int)( r8 * (double)(CIDMAX-CIDMIN) ) ;

    if ( CIDRAN >= CIDMAX  )       { goto PICKRAN ; }
    if ( CIDRAN <= CIDMIN  )       { goto PICKRAN ; }

    // Find closest unused CID to CIDRAN
    USED=1; i2 = 0 ; CIDADD = -9 ;
    while ( USED ) {

      // loop over +- 1 to try both directions from CIDRAN
      NTRY = 0 ;
      for(j=-1; j <= +1; j+=2 ) {
	CIDTMP = CIDRAN + (j*i2) ;
	if ( CIDTMP < CIDMAX  &&  CIDTMP > CIDMIN ) {
	  NTRY++ ;
	  if ( exec_cidmask(2,CIDTMP)==0 && CIDADD<0 ) { CIDADD = CIDTMP; } 
	}
      }

      if ( NTRY == 0 ) {  // check for error
	
	for(j=1; j <= CIDMAX; j++ ) {

	  if ( j < 100 || (CIDMAX-j)<100 ) {
	    printf("  xxx CID= %3d --> cidmask=%d \n",
		   j, exec_cidmask(2,j) );
	  }
	}
	printf(" xxx global MXCIDMASK = %d   NCIDMASK_LIST=%d\n", 
	       MXCIDMASK, NCIDMASK_LIST );
	printf(" xxx global CIDMASK_LIST = %d, %d, %d \n",
	       CIDMASK_LIST[0], CIDMASK_LIST[1], CIDMASK_LIST[2] ) ;

	sprintf(c1err,"Unable to try CIDTMP = %d or %d  (CIDRAN=%d,i2=%d)",
		CIDRAN-i2, CIDRAN+i2, CIDRAN, i2 );
	sprintf(c2err,"CIDMAX=%d  NPICKRAN=%d  i=%d of %d  CID=%d", 
		CIDMAX, NPICKRAN, i, NSTORE_ALL, GENLC.CID );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      i2++ ;
      if ( CIDADD > 0 ) {
	exec_cidmask(1,CIDADD); // set "USED" bit

	if ( i >= CIDRAN_OFF ) 
	  { INPUTS.CIDRAN_LIST[i-CIDRAN_OFF] = CIDADD ;  }

	/* xxxxxxxxxxxxx mark delete Jul 21 2021 xxxx
	if ( i >= INPUTS.CIDOFF ) 
	  { INPUTS.CIDRAN_LIST[i-INPUTS.CIDOFF] = CIDADD ; }
	xxxxxxxxx */

	USED = 0 ;
      }

    } // end USED block

    // - - - - - - - -  -
    // check for screen update 
    XN8    = (double)i ;
    L_STDOUT = ( fmod(XN8,XNUPD8) == 0 && i > 0 );
    if ( L_STDOUT ) {
      if ( XNUPD8 < 0.999E6 )  {  printf("%6d ", i);         }
      else          	       {  printf("%dM ", i/1000000); }
      fflush(stdout);
    }

    if ( LDMP ) {
      printf("\n xxx i=%d  i2=%d  CIDRAN=%d CIDADD=%d \n",
	     i, i2, CIDRAN, CIDADD);
    }

    //    printf(" xxx CIDRAN[%6.6d] = %6d \n", i, CIDRAN ) ;

  } // end i-loop

  printf("\n\t Done allocating %d random CIDs: %d for current sim-job \n", 
	 NSTORE_ALL, NSTORE); fflush(stdout);
  printf("\t Called random() function %d times.\n", NCALL_random );

  return ;

} // end of init_CIDRAN



// **************************
void sort_CIDRAN(void) {

  // Created Nov 30, 2011 by R.Kessler
  // sort random CIDs so that they will be generated in 
  // monotonically increasing order. This is needed to simplify
  // the merging of [FITS] samples with random CIDs.
  // Make sure to order only the randoms after CIDOFF,
  // INPUTS.CIDRAN_LIST[CIDOFF+1 ...]
  //
  // Jun 9 2013: do NOT sort for NON1a, otherwise the NON1A CIDs
  //             are ordered.
  //
  // ---------------- BEGIN --------------

  int
    NCID, CIDOFF
    ,ORDER_SORT   
    ,*INDEX_SORT
    ,*CIDRAN_TMPLIST
    ,isort, i
    ;

  if ( INDEX_GENMODEL == MODEL_NON1ASED  ) { return ; }
  if ( INDEX_GENMODEL == MODEL_NON1AGRID ) { return ; } // Apr 3 2016 

  CIDOFF = INPUTS.CIDOFF ;
  NCID   = INPUTS.NGEN_LC + INPUTS.NGENTOT_LC ;

  // allocate temp array to store sorted indices
  CIDRAN_TMPLIST = (int*)malloc( (NCID+1) * sizeof(int));
  INDEX_SORT     = (int*)malloc( (NCID+1) * sizeof(int));

  // transfer random CIDs to temporary list
  for ( i=1; i <= NCID; i++ ) 
    { CIDRAN_TMPLIST[i] = INPUTS.CIDRAN_LIST[i] ;  }

  ORDER_SORT = +1 ; // increasing order
  sortInt(NCID, &CIDRAN_TMPLIST[1], ORDER_SORT, &INDEX_SORT[1] );


  //  re-create CIDRAN_LIST with sorted CIDs 
  for ( i=1; i <= NCID; i++ ) {
    isort = INDEX_SORT[i] + 1 ; // fortran-like index
    INPUTS.CIDRAN_LIST[i] = CIDRAN_TMPLIST[isort] ;
  }

  free(CIDRAN_TMPLIST);
  free(INDEX_SORT);

} // end of sort_CIDRAN



// *****************************************
void checkpar_SIMSED(void) {

  // Dec 20, 2009:
  // make sanity checks, mainly that user sim-input
  // SIMSED-parameters match those of the model.
  // Also make sure that generation PEAK, RANGE and SIGMA 
  // are specified for each SIMSED parameter.
  //
  // If no baggage params are specified, then add them all to list.
  // If SIMGEN_DUMP is used, add all SIMSED params to SIMGEN_DUMP list.
  //
  //
  // Jun 15, 2013: for SIMSED model, do NOT check for MEAN and SIGMA
  //
  // July 28 2017:
  //  For GRIDONLY and FLAT distribution, fix the range to be 
  //  MIN-0.5*BIN to MAX+0.5*BIN so that each bin is equally
  //  populated.
  //
  // Feb 19 2018:
  //  Fix bug from last July. Use ipar_model instead of ipar_user
  //  for SIMSED lookups.
  //
  // Aug 19 2020: for GRIDONLY, set  INPUTS.GENGAUSS_SIMSED[].USE=true
  //                (bug fix)
  //

  int  GENFLAG, ipar, ipar_model, ipar_user, ipar_tmp
    ,IPAR_MODEL_USED[MXPAR_SIMSED]
    ,NBPAR, NTMP, N, NPAR_MODEL, idump, ISGRIDONLY, ISFLAT
    ;

  double  tmpdif, range_model[2], range_user[2], bin_user, SIGMA[2] ;

  char
    fnam[] = "checkpar_SIMSED"
    ,parname_model[40]  // defined in modes/SIMSED
    ,parname_user[40]   // defined in sim-input file
    ,*ptrDumpVar
    ,*ptrSIMSEDVar
    ;

  // -------------- BEGIN ---------------

  // get total number of parameters in this SIMSED model
  INPUTS.NPAR_SIMSED_MODEL = NPAR_SEDMODEL();
  NPAR_MODEL               = INPUTS.NPAR_SIMSED_MODEL;

  for ( ipar_tmp = 0; ipar_tmp < INPUTS.NPAR_SIMSED_MODEL; ipar_tmp++ )
    { IPAR_MODEL_USED[ipar_tmp] = 0; }

  for ( ipar_user = 0; ipar_user < INPUTS.NPAR_SIMSED; ipar_user++ ) {

    sprintf(parname_user,"%s", INPUTS.PARNAME_SIMSED[ipar_user] );
    GENFLAG  = INPUTS.GENFLAG_SIMSED[ipar_user];   
    SIGMA[0] = INPUTS.GENGAUSS_SIMSED[ipar_user].SIGMA[0] ;
    SIGMA[1] = INPUTS.GENGAUSS_SIMSED[ipar_user].SIGMA[1] ;

    // find which model parameter matches user parameter
    ipar_model = -9 ;
    for ( ipar_tmp = 0; ipar_tmp < NPAR_MODEL; ipar_tmp++ ) {
      fetch_parInfo_SEDMODEL(ipar_tmp, parname_model, &NBPAR, range_model );

      if ( strcmp(parname_user,parname_model) == 0 ) {
	ipar_model = ipar_tmp ;
	IPAR_MODEL_USED[ipar_model] = 1;
      }
    }

    if ( ipar_model < 0 ) {
      sprintf(c1err,"invalid user SIMSED param: '%s'", parname_user);
      sprintf(c2err,"Check sim-input file and SED.INFO file." );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }


    GENLC.SIMSED_IPARMAP[ipar_user] = ipar_model ; // Nov 18, 2011

    fetch_parInfo_SEDMODEL(ipar_model, parname_model, &NBPAR, range_model );

    // July 28 2017
    // if GRIDONLY option and flat distribution, set range to
    // min-0.5*bin to max+0.5*bin so that each bin is sample
    // with equal prob
    ISGRIDONLY = (GENFLAG & OPTMASK_GEN_SIMSED_GRIDONLY);
    ISFLAT     = (SIGMA[0] > 1.0E7 );
    if ( ISGRIDONLY && ISFLAT ) {

      bin_user      = SEDMODEL.PARVAL_BIN[ipar_model] ;
      range_user[0] = SEDMODEL.PARVAL_MIN[ipar_model] - 0.5*bin_user ;
      range_user[1] = SEDMODEL.PARVAL_MAX[ipar_model] + 0.5*bin_user;

      sprintf(INPUTS.GENGAUSS_SIMSED[ipar_user].NAME, "%s",
              SEDMODEL.PARNAMES[ipar_model] );
      INPUTS.GENGAUSS_SIMSED[ipar_user].USE      = true  ; // Aug 19 2020
 
      INPUTS.GENGAUSS_SIMSED[ipar_user].RANGE[0] = range_user[0] ;
      INPUTS.GENGAUSS_SIMSED[ipar_user].RANGE[1] = range_user[1] ;
    }

    // skip the tests below for baggage or GRIDONLY parameters
    if ( GENFLAG != OPTMASK_GEN_SIMSED_PARAM    ) { continue ; }

    // continue here only for continuous distribution
    range_user[0] = INPUTS.GENGAUSS_SIMSED[ipar_user].RANGE[0] ;
    range_user[1] = INPUTS.GENGAUSS_SIMSED[ipar_user].RANGE[1] ;

    // extend range by epsilon to allow boundary values.
    // Note this doesn't help for delta-function range.
    tmpdif = range_model[1] - range_model[0];
    range_model[0] -= (tmpdif * 1.E-8) ;
    range_model[1] += (tmpdif * 1.E-8) ;

    sprintf(c2err,"Check sim-input keyword GENRANGE_%s  (ipar_user=%d)", 
	    parname_model, ipar_user );

    // make sure that user-range does not exceed defined range of model.
    if ( range_user[0] < range_model[0] ) {
      sprintf(c1err,"User %s-min = %6.3f < model min-range=%6.3f", 
	      parname_model, range_user[0], range_model[0] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( range_user[1] > range_model[1] ) {
      sprintf(c1err,"User %s-max = %6.3f > model max-range=%6.3f", 
	      parname_model, range_user[1], range_model[1] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }


    // make sure that all information is specified.
    sprintf(c2err,"Check sim-input file.");

    if ( INDEX_GENMODEL != MODEL_SIMSED ) {
      if ( INPUTS.GENGAUSS_SIMSED[ipar_user].PEAK == -999. ) {
	sprintf(c1err,"GENMEAN_%s not specified", parname_model );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
    }

    if ( INPUTS.GENGAUSS_SIMSED[ipar_user].RANGE[0] == -999. ) {
      sprintf(c1err,"GENRANGE_%s not specified", parname_model );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    if ( INPUTS.GENGAUSS_SIMSED[ipar_user].RANGE[1] == -999. ) {
      sprintf(c1err,"GENRANGE_%s not specified", parname_model );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( INDEX_GENMODEL != MODEL_SIMSED ) {
      if ( INPUTS.GENGAUSS_SIMSED[ipar_user].SIGMA[0] == -999. ) {
	sprintf(c1err,"GENSIGMA_%s not specified", parname_model );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
      if ( INPUTS.GENGAUSS_SIMSED[ipar_user].SIGMA[1] == -999. ) {
	sprintf(c1err,"GENSIGMA_%s not specified", parname_model );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
    }

  } // end ipar loop


  // Set un-used MODEL parameters as baggage params.

  NTMP = INPUTS.NPAR_SIMSED_param;
  if ( NTMP == 0 ) {

    for ( ipar_model = 0; ipar_model < NPAR_MODEL; ipar_model++ ) {
      if ( IPAR_MODEL_USED[ipar_model] == 0 ) {     

	N = INPUTS.NPAR_SIMSED ;
	INPUTS.GENFLAG_SIMSED[N] = OPTMASK_GEN_SIMSED_param ;
	sprintf(INPUTS.KEYWORD_SIMSED[N],"SIMSED_param");
	fetch_parInfo_SEDMODEL(ipar_model, INPUTS.PARNAME_SIMSED[N],
			  &NBPAR, range_model );

	GENLC.SIMSED_IPARMAP[N] = ipar_model ; // Nov 18, 2011

	printf("\t Load baggage parameter '%s' \n",
	       INPUTS.PARNAME_SIMSED[N] );

	INPUTS.NPAR_SIMSED++ ;
	INPUTS.NPAR_SIMSED_param++ ;
	if ( INPUTS.NPAR_SIMSED >= MXPAR_SIMSED ) {
	  sprintf(c1err,"NPAR_SIMSED=%d exceeds MXPAR_SIMSED=%d",
		  INPUTS.NPAR_SIMSED, MXPAR_SIMSED ) ;
	  sprintf(c2err,"Need to increase MXPAR_SIMSED");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
	}

      }
    }

  } // NTMP=0


  // Add the SIMSED parameters to the SIMGEN_DUMP list
  // if the list is specified.

  NTMP = INPUTS.NVAR_SIMGEN_DUMP ;
  if ( NTMP > 0  && INPUTS.WRFLAG_MODELPAR ) {

    // first check that no SIMSED PARAM/param is specified ... 
    // otherwise abort.
    for ( idump=0; idump < INPUTS.NVAR_SIMGEN_DUMP; idump++ ) {
      for ( ipar = 0; ipar < INPUTS.NPAR_SIMSED; ipar++ ) {
	ptrDumpVar    = INPUTS.VARNAME_SIMGEN_DUMP[idump] ;
	ptrSIMSEDVar  = INPUTS.PARNAME_SIMSED[ipar] ;
	if ( strcmp(ptrDumpVar,ptrSIMSEDVar) == 0 ) {
	  sprintf(c1err,"SIMSED Param '%s' cannot be in SIMGEN_DUMP list",
		  ptrDumpVar);
	  sprintf(c2err,"because SIMSED Params are added automatically.");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
	}

      } // ipar
    } // idump

    // Now load all SIMSED params onto SIMGEN_DUMP list
    for ( ipar = 0; ipar < INPUTS.NPAR_SIMSED; ipar++ ) {
      idump = INPUTS.NVAR_SIMGEN_DUMP ;
      ptrDumpVar    = INPUTS.VARNAME_SIMGEN_DUMP[idump] ;
      ptrSIMSEDVar  = INPUTS.PARNAME_SIMSED[ipar] ;
      sprintf(ptrDumpVar,"%s", ptrSIMSEDVar );     
      INPUTS.NVAR_SIMGEN_DUMP++ ;
    }    
  }


  // -----------------------------------
  // For the "SIMSED_GRIDONLY: SEQUENTIAL" option, 
  // set NGEN_LC equal to the number of SEDs so that 
  // each SED is processed once.
  if ( INPUTS.OPTMASK_SIMSED == OPTMASK_GEN_SIMSED_GRIDONLY ) {
    N = NSED_SEDMODEL();
    INPUTS.NGEN = N ;
    INPUTS.NGENTOT_LC = N;
    INPUTS.NGEN_LC    = 0;

    set_screen_update(N);
    /*
    if ( INPUTS.NGEN_LC    > 0 ) { INPUTS.NGEN_LC    = N; }
    if ( INPUTS.NGENTOT_LC > 0 ) { INPUTS.NGENTOT_LC = N; }
    */
  }


  //  debugexit("add baggage params");
  
} // end of checkpar_SIMSED


// **********************************
int GENRANGE_CUT(void) {

  // Cut on user-specified GEN-ranges.
  // returns 1 if generated quantities pass cuts; 0 otherwise
  // Note that these cuts do not effect the rate calculation
  //
  // Aug 24 2015: if using HOSTLIB, require value ZTRUE.
  //              Abort if too many hosts aren't found.
  //
  // Aug 22 2018: skip REDSHIFT cut for LCLIB model.

  int LTRACE = 0; 
  int istat ;
  char fnam[] = "GENRANGE_CUT" ;

  // ----------- BEGIN ------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
    istat = 1 ; return istat ;
  }

  istat = 0 ;

  if(LTRACE) {
    printf(" xxx %s: trace CID=%d LIBID=%d\n", 
	   fnam, GENLC.CID, GENLC.SIMLIB_ID ); fflush(stdout);
  }

  if(LTRACE) { printf(" xxx %s: 0 check RA=%f \n", fnam, GENLC.RA); }
  if ( GENLC.RA < INPUTS.GENRANGE_RA[0] ) { return istat; }
  if ( GENLC.RA > INPUTS.GENRANGE_RA[1] ) { return istat; }

  if(LTRACE) { printf(" xxx %s: 1 check DEC=%f \n", fnam, GENLC.DEC); }
  if ( GENLC.DEC < INPUTS.GENRANGE_DEC[0] )  { return istat; }
  if ( GENLC.DEC > INPUTS.GENRANGE_DEC[1] )  { return istat; }

  if ( INDEX_GENMODEL != MODEL_LCLIB ) {
    if(LTRACE) {printf(" xxx %s: 2 check zCMB=%f \n",
		       fnam, GENLC.REDSHIFT_CMB);}
    if ( GENLC.REDSHIFT_CMB < INPUTS.GENRANGE_REDSHIFT[0] ) { return istat; }
    if ( GENLC.REDSHIFT_CMB > INPUTS.GENRANGE_REDSHIFT[1] ) { return istat; }
  }

  if(LTRACE) { printf(" xxx %s: 3 check PEAKMJD=%f \n", fnam, GENLC.PEAKMJD); }
  if ( GENLC.PEAKMJD < INPUTS.GENRANGE_PEAKMJD[0] )  { return istat; }
  if ( GENLC.PEAKMJD > INPUTS.GENRANGE_PEAKMJD[1] )  { return istat; }

  if ( INPUTS.HOSTLIB_USE && SNHOSTGAL.ZTRUE < 0.0 )  {  
    // if number of missing host-gals exceeds NGEN, then abort
    // to avoid infinite loop
    NGEN_REJECT.HOSTLIB++ ;
    int BIGRATIO = ( NGEN_REJECT.HOSTLIB > 2*NGENLC_WRITE );
    if ( BIGRATIO && NGENLC_WRITE > 10 ) {
      float x = (float)NGEN_REJECT.HOSTLIB / (float)NGENLC_WRITE ;
      sprintf(c1err,"%d events rejected because HOST cannot be found:",
	      NGEN_REJECT.HOSTLIB );
      sprintf(c2err,"NREJECT/NGENLC_WRITE=%d/%d = %.1f --> "
	      "need bigger HOSTLIB",
	      NGEN_REJECT.HOSTLIB, NGENLC_WRITE, x );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }
    return istat; 
  }


  if(LTRACE) { printf(" xxx %s: 4 SUCCESS \n", fnam); fflush(stdout); }

  istat = 1 ;
  return istat ;
 
} // end of GENRANGE_CUT


// **********************************
int GENMAG_CUT(void) {

  // Created Mar 16 2016
  // Cut on PEAKMAG here after generation but before
  // the light curve is digitized.
  // Returns 1 if passing cut
  // Returns 0 if failing cut

  int ifilt_obs, ifilt, PASS ;
  double PEAKMAG;
  //  char fnam[] = "GENMAG_CUT";

  // --------------- BEGIN ----------------

  PASS = 0 ;

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];    
    PEAKMAG = GENLC.peakmag_obs[ifilt_obs] ;
    if ( PEAKMAG < INPUTS.GENRANGE_PEAKMAG[0] ) { continue; }
    if ( PEAKMAG > INPUTS.GENRANGE_PEAKMAG[1] ) { continue; }
    PASS = 1;
  }

  /*
  if ( PASS == 0 ) {
    printf(" xxx fail PEAKMAG=%f for z=%f \n",
	   PEAKMAG, GENLC.REDSHIFT_CMB );
  }
  */

  return(PASS);

} // end GENMAG_CUT


// ***************************************
int geneff_calc(void) {

  // Feb 9, 2008: compute generation efficiency & error.
  //              Mainly useful for generating efficiency maps.
  //
  // Jul 17, 2009: return 0 for normal status; return 1 to STOP because
  //               effic-error is below user-allowed limit
  // 
  // Dec 22 2015: modify STOPGEN logic to be more robust. 
  //              Min allowed EFFERR -> EFFERR_STOPGEN * sqrt(1000/NGEN_LC)
  //
  float XN1, XN0, EFF, EFF_ERR, SQERR ;
  int ISTOP;

  char fnam[] = "geneff_calc" ;

  // -------- BEGIN ------

  ISTOP = 0;

  XN1 = (float)NGENLC_WRITE ;

  // normalization (XN0) include generated light curves that were
  // rejected by the search & by CUTWIN-cuts, but EXCLUDEs those
  // rejected by generation-ranges.

  XN0 = (float)(NGENLC_TOT - NGEN_REJECT.GENRANGE );

  if ( XN0 > 0.0 ) {
    EFF     = XN1 / XN0 ;
    SQERR   = XN1 * (XN0-XN1) / ( XN0*XN0*XN0 ) ;
    EFF_ERR = sqrtf( SQERR ) ;

    // effic error cannot be zero ...
    if ( EFF_ERR == 0.0 ) 
      EFF_ERR = 1.0 / XN0 ;

  }
  else {
    EFF     = 0.0 ; 
    EFF_ERR = 1.0 ;  // set large error to avoid STOPGEN

    // if NGEN_REJECT.GENRANGE is 1000, and we still have
    // not accepted anything abort with error message.

    if ( NGEN_REJECT.GENRANGE >= 1000 ) {
      sprintf(c1err,"Generated %d light curves, but NONE pass GENRANGEs",
	      NGEN_REJECT.GENRANGE );
      sprintf(c2err,"Check GENRANGEs for RA, DECL, REDSHIFT, PEAKMJD") ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }
  }

  GENLC.GENEFF    = EFF ;
  GENLC.GENEFFERR = EFF_ERR ;

  return ISTOP ;

} // end of geneff_calc

// ***************************************
int gen_cutwin(void) {

  // apply cut-windows; returns SUCCESS or ERROR.
  // ERROR just means that this SN fails cuts.
  //
  // Jan 12, 2013: replace CERNLIB sortzv with sortFloat wrapper
  //
  // Jan 31, 2014: abort on any undefined filter in CUTWIN_SNRMAX.
  //
  // Mar 1, 2014: for SNRMAX, use flux_errstat instead flux_errtot
  //              since flux_errstat is written to data file.
  //
  //
  // May 11 2015: fill GENLC.SNRMAX_TRUE[ifilt_obs]
  //
  // Aug 10 2017: malloc and free TLIST and INDEX_SORT
  //
  // Sep 25 2017: check CUTWIN_EPOCHS_xxx to select fast transients.
  //
  // May 2018: check SATURATION  cuts
  // Jun 22 2018: fux bug initializing SNRMAX_FILT[icut][ifilt]
  // Aug 19 2018: apply HOST_ZPHOT cut

  int 
    NEP, ep, i, icut, ifilt, ifilt_obs 
    ,LTMP, LSNR, NTMP, NFCUT, IFILTOBS_LIST[MXFILTINDX]
    ,NTLIST, USE_TGAP, USE_T0GAP, LDMP, opt_frame, MEM
    ;

  float 
    SNR, trueSNR
    ,SNRMAX_FILT[MXCUTWIN_SNRMAX][MXFILTINDX] 
    ,SNRMAX_TMP[MXFILTINDX] 
    ,lamrest[MXFILTINDX] 
    ,Trest, Tlast, TGAP, T0GAPOVP[2], flux, fluxerr
    ,lamobs, lamrms, lammin, lammax, z1, zhel
    ,*TLIST    
    ;

  double MJD, MJD_last, MJDDIF, GAIN, MJD_SNRMIN_FIRST, MJD_SNRMIN_LAST;

  // TLIST sorting variables
  int   DOCUTS, ORDER_SORT, isort ;
  int  *INDEX_SORT ;

  char cfilt[2];
  char fnam[] = "gen_cutwin" ;

  // ------------ BEGIN -------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return SUCCESS; }

  // bail if no cuts and no dump are specified:

  if ( INPUTS.APPLY_CUTWIN_OPT > 0  || 
       INPUTS.NVAR_SIMGEN_DUMP > 0  ||
       INPUTS.OPT_FUDGE_SNRMAX > 0    )  
    { DOCUTS = 1; }
  else 
    { DOCUTS = 0 ;    return SUCCESS ;  }

  // if we get here, calcualte CUT-WINDOW variables.
  
  LDMP = ( GENLC.CID == -991 ) ;

  if ( LDMP ) {
    printf("\n xxx ----------------------------------------------- \n");
    printf(" xxx START gen_cutwin() dump for CID=%d  LIBID=%d\n", 
	   GENLC.CID, GENLC.SIMLIB_ID );
    printf(" xxx NEP=%d  zHel=%.4f  PEAKMJD=%.3f  MWEBV=%.3f\n",   
	   GENLC.NEPOCH, GENLC.REDSHIFT_HELIO, 
	   GENLC.PEAKMJD, GENLC.MWEBV );
    fflush(stdout);
  }

  zhel  = GENLC.REDSHIFT_HELIO ;
  z1 = ( 1. + zhel );

  NEP = GENLC.NEPOCH;

  // INDEX_SORT is used for epochs and also for filters,
  // so make sure the MEM size can handle all filters.
  MEM=NEP; if (NEP < MXFILTINDX) { MEM=MXFILTINDX+1; }

  TLIST      = malloc( MEM * sizeof(float) ) ;
  INDEX_SORT = malloc( MEM * sizeof(int)   ) ;

  for ( ifilt_obs=0; ifilt_obs < MXFILTINDX; ifilt_obs++ ) {
    SNRMAX_TMP[ifilt_obs]           = -9.0 ;
    GENLC.SNRMAX_FILT[ifilt_obs]    = -9.0 ;  
    GENLC.SNRMAX_SORTED[ifilt_obs]  = -9.0 ;  
    GENLC.SNRMAX_TRUE[ifilt_obs]    = -9.0 ;  
    for ( icut=0; icut <= MXCUTWIN_SNRMAX; icut++ )
      { SNRMAX_FILT[icut][ifilt_obs] = -9.0 ; }
  }
  GENLC.TIME_ABOVE_SNRMIN = 0.0 ;


  // get lamrest for each observer filter
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    opt_frame = GENFRAME_OBS;
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];    
    get_filtlam__(&opt_frame, &ifilt_obs,  &lamobs,&lamrms,&lammin,&lammax);
    lamrest[ifilt_obs] = lamobs / z1 ;
  }

  // make redshift cut on true and/or final(measured) 
  int LZTRUE=0, LZFINAL=0, LZPHOT=0 ;
  double ZTRUE, ZFINAL, ZPHOT ;
  ZTRUE  = GENLC.REDSHIFT_CMB;
  ZFINAL = GENLC.REDSHIFT_CMB_SMEAR ;
  ZPHOT  = SNHOSTGAL.ZPHOT ; 

  LZTRUE = ( ZTRUE >= INPUTS.CUTWIN_REDSHIFT_TRUE[0] && 
	     ZTRUE <= INPUTS.CUTWIN_REDSHIFT_TRUE[1] ) ;
  LZFINAL = ( ZFINAL >= INPUTS.CUTWIN_REDSHIFT_FINAL[0] && 
	      ZFINAL <= INPUTS.CUTWIN_REDSHIFT_FINAL[1] ) ;
  LZPHOT  = ( ZPHOT >= INPUTS.CUTWIN_HOST_ZPHOT[0] && 
	      ZPHOT <= INPUTS.CUTWIN_HOST_ZPHOT[1] ) ;

  if ( LZTRUE && LZFINAL && LZPHOT ) 
    { GENLC.CUTBIT_MASK |= (1 << CUTBIT_REDSHIFT); }

  // start loop over epochs
  NTLIST    = 0 ;
  MJD_last  = -99999. ;
  MJD_SNRMIN_FIRST = MJD_SNRMIN_LAST = 0.0 ;

  for ( ep = 1; ep <= NEP ; ep++ ) {

    Trest  = (float)GENLC.epoch_rest[ep] ;
    MJD    = GENLC.MJD[ep];
    MJDDIF = MJD - MJD_last ;

    // require valid epoch to continue

    if ( !GENLC.OBSFLAG_GEN[ep] ) { continue ; }
    
    if ( Trest < INPUTS.CUTWIN_TRESTMIN[0] ) { continue ; }
    if ( Trest > INPUTS.CUTWIN_TRESTMAX[1] ) { continue ; }

    ifilt_obs = GENLC.IFILT_OBS[ep] ;
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );

    // keep track of the true/generated SNR
    trueSNR     = GENLC.trueSNR[ep] ;
    if ( trueSNR > GENLC.SNRMAX_TRUE[ifilt_obs] ) {
      GENLC.SNRMAX_TRUE[ifilt_obs] = trueSNR ;
    }

    if ( LDMP ) {
      printf(" xxx ep=%d   MJD(%s)= %f   Trest=%f  \n", 
	     ep, cfilt, MJD, Trest );   
      fflush(stdout); 
    }

    // apply rest-frame lambda cut to all epochs 
    // used in TREST & SNRMAX cuts
    if ( lamrest[ifilt_obs] < INPUTS.EPCUTWIN_LAMREST[0] ) continue ;
    if ( lamrest[ifilt_obs] > INPUTS.EPCUTWIN_LAMREST[1] ) continue ;

    // apply SNRMIN_EP cut to all epochs
    flux    = GENLC.flux[ep] ;
    fluxerr = GENLC.fluxerr_data[ep] ;
    SNR     = flux / fluxerr ;

    if ( SNR < INPUTS.EPCUTWIN_SNRMIN[0] ) continue ;

    NTLIST++;
    TLIST[NTLIST] = Trest; // used later for sorting

    if ( Trest < GENLC.TRESTMIN ) { GENLC.TRESTMIN = Trest ; }
    if ( Trest > GENLC.TRESTMAX ) { GENLC.TRESTMAX = Trest ; }

    if ( Trest >= INPUTS.CUTWIN_TRESTMIN[0] &&
	 Trest <= INPUTS.CUTWIN_TRESTMIN[1] )
      { GENLC.CUTBIT_MASK |= (1 << CUTBIT_TRESTMIN) ; }

    if ( Trest >= INPUTS.CUTWIN_TRESTMAX[0] &&
	 Trest <= INPUTS.CUTWIN_TRESTMAX[1] )
      { GENLC.CUTBIT_MASK  |= (1 << CUTBIT_TRESTMAX) ; }

    // check epochs above special SNR cut (not the global EPCUTWIN_SNRMIN)
    if ( SNR >= INPUTS.CUTWIN_NEPOCH[1] ) { GENLC.NOBS_SNR++; }

    // evaluate SNRMAX for each Trest window & filter

    for ( icut=0; icut <= INPUTS.NCUTWIN_SNRMAX; icut++ ) {

      LTMP = ( Trest >= INPUTS.CUTWIN_SNRMAX_TREST[icut][0])
	&&   ( Trest <= INPUTS.CUTWIN_SNRMAX_TREST[icut][1]) ;

      LSNR = SNR  >  SNRMAX_FILT[icut][ifilt_obs] ;

      if (  LSNR   && LTMP ) {
	SNRMAX_FILT[icut][ifilt_obs] = SNR ; 
	if ( icut == 0 ) {
	  GAIN = SIMLIB_OBS_GEN.CCDGAIN[ep] ;
	  GENLC.FLUXADU_at_SNRMAX[ifilt_obs] = flux ; 
	  GENLC.FLUXpe_at_SNRMAX[ifilt_obs]  = flux * GAIN ;
	}
      }
    } // icut

    // count NOBS with at least MJDDIF time separation
    if ( MJDDIF >= INPUTS.CUTWIN_MJDDIF[0] ) 
      { GENLC.NOBS_MJDDIF++ ; }

    // count time above SNRMIN
    if ( SNR > INPUTS.CUTWIN_EPOCHS_SNRMIN &&
	 strstr(INPUTS.CUTWIN_EPOCHS_FILTERS,cfilt)!=NULL  ) {
      if ( MJD_SNRMIN_FIRST < 0.01 ) { MJD_SNRMIN_FIRST = MJD; }
      MJD_SNRMIN_LAST = MJD ;
    }

    MJD_last = MJD;

  } // end of ep loop

  // - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ( GENLC.NOBS_SNR >= INPUTS.CUTWIN_NEPOCH[0] ) 
    { GENLC.CUTBIT_MASK |= (1<<CUTBIT_NOBS_SNR) ; }

  if ( GENLC.NOBS_MJDDIF >= INPUTS.CUTWIN_NOBSDIF[0] ) 
    { GENLC.CUTBIT_MASK |= (1<<CUTBIT_NOBS_MJDDIF) ; }

  GENLC.TIME_ABOVE_SNRMIN = (float)(MJD_SNRMIN_LAST - MJD_SNRMIN_FIRST);
  if ( GENLC.TIME_ABOVE_SNRMIN < INPUTS.CUTWIN_EPOCHS_TRANGE[1] ) 
    { GENLC.CUTBIT_MASK |= (1<<CUTBIT_TIME_ABOVE) ; }


  // now sort SNRMAX for each filter and store in GENLC.SNRMAX_SORTED
  // First get SNRMAX-list sorted by filter.
  // These are computed just for the SIMGEN_DUMP option,
  // and are NOT used for any of the CUTWIN tests.

  icut = 0 ;
  NTMP = GENLC.NFILTDEF_OBS;
  for ( ifilt=0; ifilt < NTMP ; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );

    SNRMAX_TMP[ifilt+1]          = SNRMAX_FILT[icut][ifilt_obs] ;
    GENLC.SNRMAX_FILT[ifilt_obs] = SNRMAX_FILT[icut][ifilt_obs] ;

    if ( LDMP ) {
      printf(" xxx CID=%d  SNRMAX[%s] = %6.2f \n", 
	     GENLC.CID, cfilt, SNRMAX_FILT[icut][ifilt_obs] );
      fflush(stdout);
    }
  }


  ORDER_SORT = -1; // decreasing order
  sortFloat( NTMP, &SNRMAX_TMP[1], ORDER_SORT, &INDEX_SORT[1] );

  for ( ifilt = 1; ifilt <= NTMP ; ifilt++ ) {
    isort = INDEX_SORT[ifilt] + 1 ; // fortran-line index
    SNR   = SNRMAX_TMP[isort] ;
    GENLC.SNRMAX_SORTED[ifilt] = SNR ;
  }


  // -------------
  // Check if SNRMAX passes cut in each passband for each test
  // Set default LSNR=TRUE; then set false if any SNRMAX cut fails. 

    LSNR = 1;

    for ( icut=1; icut <= INPUTS.NCUTWIN_SNRMAX; icut++ ) {

      // get filter-range to test for this cut
      NFCUT = PARSE_FILTLIST( INPUTS.CUTWIN_SNRMAX_FILTERS[icut], 
			      IFILTOBS_LIST ); // <== returned

      NTMP = 0;
      for ( i=0; i < NFCUT; i++ ) {

	ifilt_obs = IFILTOBS_LIST[i] ;

	// abort if ifilt_obs is undefined
	ifilt = GENLC.IFILTINVMAP_OBS[ifilt_obs]; 
	if ( ifilt < 0 ) {
	  sprintf(c1err,"Invalid filter '%c' in", FILTERSTRING[ifilt_obs]);
	  sprintf(c2err,"CUTWIN_SNRMAX: %.1f  %s   etc ...",
		  INPUTS.CUTWIN_SNRMAX[icut][0],
		  INPUTS.CUTWIN_SNRMAX_FILTERS[icut] );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
	}

	SNR = SNRMAX_FILT[icut][ifilt_obs] ;
	if ( SNR > INPUTS.CUTWIN_SNRMAX[icut][0] )  { NTMP++ ; }

      } // end of ifilt_obs loop


      if ( NTMP < INPUTS.CUTWIN_SNRMAX_NFILT[icut] ) { LSNR = 0 ; }


    } // end of icut loop


    if ( LSNR == 1 ) 
      { GENLC.CUTBIT_MASK  |= (1<<CUTBIT_SNRMAX) ; }

    // ========================
    // To get max gap, need sort TLIST
    // INDEX_SORT is returned from CERNLIB function sortzv_

    ORDER_SORT = +1; // increasing order
    sortFloat(NTLIST, &TLIST[1], ORDER_SORT, &INDEX_SORT[1] );

    Tlast = -999999. ;


    for ( i=1; i<= NTLIST; i++ ) {

      isort = INDEX_SORT[i] + 1; // fortran-like index
      Trest = TLIST[isort];

      // valid TGAP requires both Trest & Tlast to be within usable 
      // light curve defined by TRESTMIN & TRESTMAX cuts

      if ( Trest < INPUTS.CUTWIN_TRESTMIN[0] ) goto SKIPGAP;
      if ( Trest > INPUTS.CUTWIN_TRESTMAX[1] ) goto SKIPGAP;
      if ( Tlast < INPUTS.CUTWIN_TRESTMIN[0] ) goto SKIPGAP;
      if ( Tlast > INPUTS.CUTWIN_TRESTMAX[1] ) goto SKIPGAP;

      TGAP = Trest - Tlast; 
      USE_TGAP = USE_T0GAP = 0;

      // use T0-gap only of it covers the lower range of TRESTMAX-window,
      // or the upper range of the TRESTMIN-window.

      T0GAPOVP[0] = INPUTS.CUTWIN_TRESTMIN[1] ; // max TRESTMIN window 
      T0GAPOVP[1] = INPUTS.CUTWIN_TRESTMAX[0] ; // min TRESTMAX window
      if ( Trest >= T0GAPOVP[0] && Trest <= T0GAPOVP[1] ) USE_T0GAP = 1 ;
      if ( Tlast >= T0GAPOVP[0] && Tlast <= T0GAPOVP[1] ) USE_T0GAP = 1 ;
      if ( Trest >= T0GAPOVP[1] && Tlast <= T0GAPOVP[0] ) USE_T0GAP = 1 ;

      // use T-gap anywhere along light curve
      USE_TGAP = 1 ;

      if ( USE_T0GAP == 1 && TGAP > GENLC.T0GAPMAX )  
	GENLC.T0GAPMAX = TGAP ;

      if ( USE_TGAP == 1 && TGAP > GENLC.TGAPMAX )  
	GENLC.TGAPMAX = TGAP ;

      /*
      if ( LDMP > 0 ) {
	printf(" xxx sorted Trest(%2d) = %7.3f  TGAP=%6.3f ", 
	       i, Trest, TGAP );
	printf(" USE_T0GAP=%d USE_TGAP=%d\n", USE_T0GAP, USE_TGAP );
      }
      */

    SKIPGAP:
      Tlast = Trest;
    } // end i-loop over TLIST




    if ( LDMP > 0 ) {
      printf(" xxx CID=%d : TGAPMAX= %6.2f   T0GAPMAX=%6.2f\n", 
	     GENLC.CID, GENLC.TGAPMAX, GENLC.T0GAPMAX );
      printf(" xxx CID=%d : TRESTMIN=%f TRESTMAX=%f \n",
	     GENLC.CID, GENLC.TRESTMIN, GENLC.TRESTMAX );
      fflush(stdout);
    }

    // set TGAP cutbits

    if ( GENLC.TGAPMAX <= INPUTS.CUTWIN_TGAPMAX[1] && 
	 GENLC.TGAPMAX >= INPUTS.CUTWIN_TGAPMAX[0]  ) 
      { GENLC.CUTBIT_MASK  |= (1<<CUTBIT_TGAPMAX) ; }

    if ( GENLC.T0GAPMAX <= INPUTS.CUTWIN_T0GAPMAX[1] && 
	 GENLC.T0GAPMAX >= INPUTS.CUTWIN_T0GAPMAX[0]  ) 
      { GENLC.CUTBIT_MASK  |= (1<<CUTBIT_T0GAPMAX) ; }

    // Set MWEBV cutbit

    if ( GENLC.MWEBV <= INPUTS.CUTWIN_MWEBV[1] && 
	 GENLC.MWEBV >= INPUTS.CUTWIN_MWEBV[0]  ) 
      { GENLC.CUTBIT_MASK  |= (1<<CUTBIT_MWEBV) ; }

    // Set PEAKMAG cutbit if any filter satisfies (Sep 2012)
    int NPASS_PEAKMAG=0, NFAIL_PEAKMAG_ALL=0;
    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {      
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];    
      if ( gen_cutwin_PEAKMAG(1,ifilt_obs) == SUCCESS ) // OR logic
	{ NPASS_PEAKMAG++; } 
      if ( gen_cutwin_PEAKMAG(2,ifilt_obs) != SUCCESS )  // AND logic
	{ NFAIL_PEAKMAG_ALL++; } 

    }
    if ( NPASS_PEAKMAG > 0 && NFAIL_PEAKMAG_ALL == 0 ) 
      {  GENLC.CUTBIT_MASK  |= (1<<CUTBIT_PEAKMAG) ; }

    // May 2018: check SATURATION  cuts 
    int MINOBS, MAXOBS, NOBS_TOT, NOBS, IPASS=1 ;
    for(icut=0; icut < INPUTS.NCUTWIN_SATURATE; icut++ ) {
      MINOBS = INPUTS.CUTWIN_SATURATE_NOBS[icut][0] ;
      MAXOBS = INPUTS.CUTWIN_SATURATE_NOBS[icut][1] ;
      NFCUT = PARSE_FILTLIST( INPUTS.CUTWIN_SATURATE_FILTERS[icut], 
			      IFILTOBS_LIST ); // <== returned   
      for ( i=0; i < NFCUT; i++ ) {
	ifilt_obs = IFILTOBS_LIST[i] ;
	NOBS   = GENLC.NOBS_SATURATE_FILTER[1][ifilt_obs];
	if ( NOBS < MINOBS || NOBS > MAXOBS ) { IPASS=0; }
      }
    } // end icut loop


    // repeat for NOBS_NOSATURATE
    for(icut=0; icut < INPUTS.NCUTWIN_NOSATURATE; icut++ ) {
      MINOBS = INPUTS.CUTWIN_NOSATURATE_NOBS[icut][0] ;
      MAXOBS = INPUTS.CUTWIN_NOSATURATE_NOBS[icut][1] ;
      NFCUT = PARSE_FILTLIST( INPUTS.CUTWIN_NOSATURATE_FILTERS[icut], 
			      IFILTOBS_LIST ); // <== returned      
      for ( i=0; i < NFCUT; i++ ) {
	ifilt_obs = IFILTOBS_LIST[i] ;
	NOBS_TOT  = GENLC.NOBS_FILTER[ifilt_obs] ;

	// require enough NOBS_TOT in this band in case
	// the model is undefined and we have NOBS_TOT=0.
	if ( NOBS_TOT > MINOBS ) {
	  NOBS      = GENLC.NOBS_SATURATE_FILTER[0][ifilt_obs];
	  if ( NOBS < MINOBS || NOBS > MAXOBS ) { IPASS=0; }
	}
	  /* 
	printf(" xxx ifilt_obs=%d  NOBS=%d (%d-%d)  IPASS=%d\n", 
	ifilt_obs, NOBS, MINOBS, MAXOBS, IPASS ); */
      }

    } // end icut loop

    //    printf(" xxx IPASS(NOSATURATE) = %d \n", IPASS);
    if ( IPASS ) { GENLC.CUTBIT_MASK  |= (1<<CUTBIT_SATURATE) ; }

    // -------------------------------------------

    if ( LDMP  ) {
      printf(" xxxxx CID=%d  CUTBIT_MASK = %d  (requires %d) \n",
	     GENLC.CID, GENLC.CUTBIT_MASK, ALLBIT_CUTMASK );
      for ( ifilt_obs  = IFILTOBS_LIST[0]; 
	    ifilt_obs <= IFILTOBS_LIST[1]; ifilt_obs++ ) {
	SNR = SNRMAX_FILT[1][ifilt_obs] ;
	printf(" xxxx CID = %d, ifilt_obs = %d, SNRMAX(cut1) = %f \n", 
	       GENLC.CID, ifilt_obs, SNR);
	SNR = SNRMAX_FILT[2][ifilt_obs] ;
	printf(" xxxx CID = %d, ifilt_obs = %d, SNRMAX(cut2) = %f \n", 
	       GENLC.CID, ifilt_obs, SNR);
	fflush(stdout);
      } // end of ifilt_obs loop
    }

    // free memory
    free(TLIST);  free(INDEX_SORT);

    if ( LDMP ) {
      printf(" xxx END OF gen_cutwin() dump for CID=%d \n\n", GENLC.CID );
    }

    // fail event if SNRMAX-fudge is requested, but exposure time is normal
    // Next time around with adjusted EXPOSURE_TIME, the event of processed
    if ( GENLC.FUDGE_SNRMAX_FLAG == 1 ) 
      { return ERROR; }

    if ( INPUTS.APPLY_CUTWIN_OPT == 0 )   
      { return SUCCESS ; }


    // if we get here, apply cuts.

    if ( GENLC.CUTBIT_MASK == ALLBIT_CUTMASK ) 
      {  return SUCCESS ;  }
    else
      { return ERROR ; }


} // end of gen_cutwin


// ******************************************
int gen_cutwin_PEAKMAG(int OPT, int ifilt_obs) {

  // Created Sep 2015
  // Apply user-specified PEAKMAG cut(s);
  // return SUCCESS or ERROR.
  //
  // OPT=1 --> test with CUTWIN_PEAKMAG
  // OPT=2 --> test with CUTWIN_PEAKMAG_ALL
  //
  // First apply global CUTWIN_PEAKMAG cut.
  // Then check field-dependent CUTWIN_PEAKMAG_BYFIELD cuts
  // to apply stricter cut.
  //
  // Aug 31 2018: for LCLIB, apply same cut to template mag.

  int PASS = 1;
  float MLO, MHI ;
  double MCOR_TRUE_MW  = GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs]; 
  float  peakmag       = (float)GENLC.peakmag_obs[ifilt_obs] ;
  float  mag_T         = (float)GENLC.genmag_obs_template[ifilt_obs];

// ------------- BEGIN ---------------

  // correct MWCOR if option is set. Default MCOR is 0
  peakmag -=  (float)MCOR_TRUE_MW ;
  mag_T   -=  (float)MCOR_TRUE_MW ;

  if ( OPT == 1 ) {
    // OR logic; only 1 filter must pass cut
     MLO     = INPUTS.CUTWIN_PEAKMAG[0] ; 
     MHI     = INPUTS.CUTWIN_PEAKMAG[1] ; 
  }
  else {
    // AND logic; all filters must pass cut
     MLO     = INPUTS.CUTWIN_PEAKMAG_ALL[0] ; 
     MHI     = INPUTS.CUTWIN_PEAKMAG_ALL[1] ; 
  }

  if ( peakmag < MLO ) { PASS = 0 ; }
  if ( peakmag > MHI ) { PASS = 0 ; }

  // for LCLIB model, check template mags too
  if ( INDEX_GENMODEL == MODEL_LCLIB ) {
    if ( mag_T < MLO ) { PASS = 0 ; }
    if ( mag_T > MHI ) { PASS = 0 ; }
  }

  // skip FIELD-dependent cuts for AND logic
  if ( OPT !=1 ) { goto DONE; }

  // now check field-dependent PEAKMAG cut(s)
  int  i, N;
  char *field, *fieldList;
  N      = INPUTS.NCUTWIN_PEAKMAG_BYFIELD;
  field  = GENLC.FIELDNAME[1]; // get field from header/1st epoch

  for(i=1; i<=N; i++ ) {
    fieldList = INPUTS.CUTWIN_BYFIELDLIST[i];
    if ( strstr(fieldList,field) != NULL ) {
      if ( peakmag < INPUTS.CUTWIN_PEAKMAG_BYFIELD[i][0] ) { PASS = 0 ; }
      if ( peakmag > INPUTS.CUTWIN_PEAKMAG_BYFIELD[i][1] ) { PASS = 0 ; }

      /*
      printf(" xxx compare '%s' to '%s' (PASS=%d)\n", 
      field, fieldList, PASS );  */
    }
  }

  // return appropriate function code
 DONE:

  if( PASS ) 
    { return SUCCESS; }
  else
    { return ERROR; }


} // end  gen_cutwin_PEAKMAG

// ******************************************
void  LOAD_SEARCHEFF_DATA(void) {

  // Created Jan , 2014
  // Load structures
  //   - SEARCHEFF_DATA 
  //   - SEARCHEFF_RANDOMS
  // to be used by  SEARCHEFF_xxx functions to evaluate
  // pipeline and SPEC efficiencies.
  //
  // Aug 24 2014: 
  //   major fix, use SNR_CALC instead of measured SNR.
  //   But for SDSS, keep using measured SNR so that spec-effic is OK.
  //
  // Jun 23 2016: load HOSTMAG and HOSTSB
  // Jan 03 2018: load NPE_SAT
  // Jan 15 2018: load GENLC.FIELDNAME[0], not FIELDNAME[1]
  // Jun 18 2018: load SNRMAX
  //
  // Dec 22 1019: 
  //  + loop over GENLC.NEPOCH intead of GENLC.NEWMJD -->
  //    more efficiency, and beware change of random sync.
  //
  // Feb 05 2021: load each overlap field and NFIELD_OVP

  bool ISCORR_PHOTRPBOB = (INPUTS_SEARCHEFF.NREDUCED_CORR_PHOTPROB > 0);
  int  NMAP_PHOTPROB    = INPUTS_SEARCHEFF.NMAP_PHOTPROB;

  int ep, NOBS,  NRANTMP=0;
  double flux, flux_err, SNR_CALC, SNR_MEAS, SNR, oldRan ;
  //  char fnam[] = "LOAD_SEARCHEFF_DATA";

  // --------------- BEGIN ----------------

  SEARCHEFF_DATA.CID        = GENLC.CID ;
  SEARCHEFF_DATA.REDSHIFT   = GENLC.REDSHIFT_HELIO ;
  SEARCHEFF_DATA.PEAKMJD    = GENLC.PEAKMJD ;
  SEARCHEFF_DATA.DTPEAK_MIN = GENLC.DTPEAK_MIN ; // closest T-Tpeak
  SEARCHEFF_DATA.SALT2mB    = GENLC.SALT2mB ;
  SEARCHEFF_DATA.SNRMAX     = GENLC.SNRMAX_GLOBAL ;


  // load field(s) and be careful about overlaps (e.g., X1+X3)

  // xxx  sprintf(SEARCHEFF_DATA.FIELDNAME, "%s", GENLC.FIELDNAME[0] );
  int ifield, NFIELD_OVP = SIMLIB_HEADER.NFIELD_OVP ;
  sprintf(SEARCHEFF_DATA.FIELDNAME, "%s", SIMLIB_HEADER.FIELD );
  SEARCHEFF_DATA.NFIELD_OVP = NFIELD_OVP ;
  for(ifield=0; ifield < NFIELD_OVP; ifield++ ) {
    sprintf(SEARCHEFF_DATA.FIELDLIST_OVP[ifield], "%s",
	    SIMLIB_HEADER.FIELDLIST_OVP[ifield] );
  }
  // - - - - - - - -
  NOBS = 0 ;

  for(ep=1; ep <= GENLC.NEPOCH; ep++ ) {

    if ( !GENLC.OBSFLAG_GEN[ep] ) { continue; } 

    SNR_CALC = GENLC.SNR_CALC[ep] ; // Aug 24, 2014

    flux      = GENLC.flux[ep] ;
    flux_err  = GENLC.fluxerr_data[ep] ;
    SNR_MEAS  = -9.0 ;
    if ( flux_err > 0.0 ) { SNR_MEAS = flux / flux_err ; }
   
    SNR = SNR_CALC ;
      
    // for SDSS, continue using wrong SNR based on measured flux
    // so that the spec-efficiency function is still correct.
    if ( strcmp(GENLC.SURVEY_NAME,"SDSS") == 0 ) { SNR = SNR_MEAS; }
    
    NOBS = ep-1; // Dec 22 2019
    SEARCHEFF_DATA.IFILTOBS[NOBS]  = GENLC.IFILT_OBS[ep] ;     
    SEARCHEFF_DATA.MJD[NOBS]       = GENLC.MJD[ep] ;
    SEARCHEFF_DATA.MAG[NOBS]       = GENLC.genmag_obs[ep] ; 
    SEARCHEFF_DATA.SNR[NOBS]       = SNR ;
    SEARCHEFF_DATA.NPE_SAT[NOBS]   = GENLC.npe_above_sat[ep];
    
    oldRan = SEARCHEFF_RANDOMS.FLAT_PIPELINE[NOBS] ;
    if ( oldRan < -0.001 ) 
      { SEARCHEFF_RANDOMS.FLAT_PIPELINE[NOBS] = getRan_Flat1(1);  NRANTMP++ ; }
    

    if ( NMAP_PHOTPROB > 0 ) {
      // load Gaussian randoms for correlated PHOTPROB
      if ( ISCORR_PHOTRPBOB ) {
	oldRan = SEARCHEFF_RANDOMS.GAUSS_PHOTPROB[NOBS] ;
	if ( oldRan < -998.0 ) 
	  { SEARCHEFF_RANDOMS.GAUSS_PHOTPROB[NOBS] = getRan_Gauss(1); }
      }
      else {
	// load flat randoms for uncorrelated PHOTPROB
	oldRan = SEARCHEFF_RANDOMS.FLAT_PHOTPROB[NOBS] ;
	if ( oldRan < -998.0 ) 
	  { SEARCHEFF_RANDOMS.FLAT_PHOTPROB[NOBS] = getRan_Flat1(1); }	
      }
    } // end NMAP_PHOTPROB

  } // end ep loop over epochs

  SEARCHEFF_DATA.NOBS =  GENLC.NEPOCH ;


  // load SPEC-EFF randoms and filter-dependent quantities

  int ifilt, ifilt_obs ;

  for ( ifilt=0; ifilt <= MXFILTINDX; ifilt++ ) {

    if ( SEARCHEFF_RANDOMS.FLAT_SPEC[ifilt] < -0.01 ) 
      { SEARCHEFF_RANDOMS.FLAT_SPEC[ifilt]  = getRan_Flat1(1); }

    if ( ifilt == MXFILTINDX ) { continue ; } // avoid array overwrite
    SEARCHEFF_DATA.PEAKMAG[ifilt] = MAG_UNDEFINED ;
    SEARCHEFF_DATA.HOSTMAG[ifilt] = MAG_UNDEFINED ;
    SEARCHEFF_DATA.SBMAG[ifilt]   = MAG_UNDEFINED ;
    
  }
  
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      /*
      printf(" xxx ifilt_obs=%d: MAG[PEAK,HOST,SB] = %.3f, %.3f, %.3f \n",
	     ifilt_obs, GENLC.peakmag_obs[ifilt_obs],
	     SNHOSTGAL.GALMAG[ifilt_obs][0], SNHOSTGAL.SB_MAG[ifilt_obs] );
      */
      SEARCHEFF_DATA.PEAKMAG[ifilt_obs] =  GENLC.peakmag_obs[ifilt_obs] ;
      SEARCHEFF_DATA.HOSTMAG[ifilt_obs] =  SNHOSTGAL.GALMAG[ifilt_obs][0] ;
      SEARCHEFF_DATA.SBMAG[ifilt_obs] =    SNHOSTGAL.SB_MAG[ifilt_obs];
  }  // ifilt

  return ;

} // end of LOAD_SEARCHEFF_DATA



// ******************************************
void gen_spectype(void) {

  // Created July 2011
  // Major modif. Mar 2013
  // Set GENLC.SNTYPE based on spec vs. phot typing.
  //
  // Jun 7, 2013: 
  //  - for photoid, add OFFSET_TYPE_PHOT to SNTYPE.
  //  - set GENLC.SNTYPE for CC as well as for Ia
  //
  //
  // Jun 2 2018: if IFLAG_SPEC_EFFZERO, then set PHOTID flag

  int L_PHOTID, ispgen ;
  //  char fnam[] = "gen_spectype" ;

  // ---------- BEGIN --------------

  // if SPECTYPE has not been specified, then
  // leave GENLC.SNTYPE as-is => 100% spec-typing.
  // However, for SPECTYPE_SEARCH_FLAG, then use search 
  // eff for spec-typing.


  L_PHOTID =  ( INPUTS_SEARCHEFF.NMAP_SPEC > 0  &&
		GENLC.SEARCHEFF_MASK  != 3  ) ;

  // un 2018: if user forces EFF_SPEC=0 without a map, then we have PHOTID
  if ( INPUTS_SEARCHEFF.IFLAG_SPEC_EFFZERO ) { L_PHOTID = 1; }

  ispgen = GENLC.NON1ASED.ISPARSE ;
  

  // set SNTYPE 

  if ( INDEX_GENMODEL == MODEL_FIXMAG ) {
    GENLC.SNTYPE   = MODEL_FIXMAG ;
  }
  else if ( LGEN_SNIA )  { 
    GENLC.SNTYPE   = INPUTS.SNTYPE_Ia_SPEC ; 
  }
  else if ( ispgen >=0 ) {
    GENLC.SNTYPE = INPUTS.NON1ASED.SNTAG[ispgen];
  }

  // user-defined GENTYPE over-rides any other type (Apr 11 2017)
  if ( INPUTS.GENTYPE_SPEC > 0 ) {
    GENLC.SNTYPE = INPUTS.GENTYPE_SPEC ; 
  }

 
  if ( L_PHOTID ) {
    // photometric candidate
    GENLC.METHOD_TYPE = METHOD_TYPE_PHOT ;
    GENLC.NTYPE_PHOT++ ; 
    if ( LGEN_SNIA )  
      { GENLC.SNTYPE = INPUTS.SNTYPE_Ia_PHOT ; }
    else
      { GENLC.SNTYPE += OFFSET_TYPE_PHOT ; }

    // check over-ride from GENTYPE (Apr 11 2017)
    if ( INPUTS.GENTYPE_PHOT > 0 ) 
      { GENLC.SNTYPE = INPUTS.GENTYPE_PHOT ; }

    setz_unconfirmed();
  }
  else { 
    // spec-confirmed
    GENLC.REDSHIFT_FLAG = REDSHIFT_FLAG_SNSPEC ; 
    GENLC.METHOD_TYPE   = METHOD_TYPE_SPEC ;
    GENLC.NTYPE_SPEC++ ; 
  }


} // end gen_spectype


// ***************************************
void  setz_unconfirmed(void) {

  // Use GENLC.SEARCHEFF_MASK (see sntoolc_trigger.c) 
  // to determine if this UNCONFIRMED SN has a host-galaxy 
  // redshift ... and if so, apply optional WRONGHOST model.
  // If there is no zHOST then set measured redshift and
  // its error to HOST PHOTOZ, or to -9.
  //
  //
  // July 26 2017: 
  //  Fix bug when there is no zSPEC(SN or host), but there is
  //  zPhotHost (SNHOSTGAL.ZPHOT).
  //
  // Mar 17 2018: z -> HOST_PHOTOZ if no host-zSPEC
  //
  // May 30 2018: 
  //  + GENLC.REDSHIFT_CMB[HELIO] -> -9 of ZHOST < 0
  //    This doesn't fix any bug, but avoids wierd small-negative 
  //    ZCMB & ZHELIO that look like a bug.
  //

  double RA        = GENLC.RA ;
  double DEC       = GENLC.DEC ;
  double ZPHOT     = SNHOSTGAL.ZPHOT ;
  double ZPHOT_ERR = SNHOSTGAL.ZPHOT_ERR ;
  double ZCMB_ORIG = GENLC.REDSHIFT_CMB_SMEAR ;
  int    FOUND_zHOST ;
  int LDMP = 0 ;
  char eq[]   = "eq" ;
  char fnam[] = "setz_unconfirmed" ;

  // ---------- BEGIN -------------
  
  FOUND_zHOST = ( (GENLC.SEARCHEFF_MASK & APPLYMASK_SEARCHEFF_zHOST) > 0 ) ;
  
  if ( LDMP ) {
    printf(" xxx -------------------------- \n");
    printf(" xxx %s: CID=%d  Found_zHOST=%d   CORRECT_HOSTMATCH=%d \n",
	   fnam, GENLC.CID, FOUND_zHOST, GENLC.CORRECT_HOSTMATCH );
  }

  if ( FOUND_zHOST ) {

    GENLC.REDSHIFT_FLAG = REDSHIFT_FLAG_HOSTSPEC ; 

    // update only if we have the wrong host
    if ( !GENLC.CORRECT_HOSTMATCH  ) {
      GENLC.REDSHIFT_FLAG        = REDSHIFT_FLAG_WRONGHOST ; 
      GENLC.REDSHIFT_HELIO_SMEAR = SNHOSTGAL.ZSPEC ;  

      if ( INPUTS.VEL_CMBAPEX > 1.0 ) {
	GENLC.REDSHIFT_CMB_SMEAR =  
	  zhelio_zcmb_translator(SNHOSTGAL.ZSPEC, RA,DEC,eq, +1);
      }
      else {
	GENLC.REDSHIFT_CMB_SMEAR = SNHOSTGAL.ZSPEC ;
      }

      if ( LDMP ) {
	printf(" xxx zcmb = %.3f(SN) -> %.3f(WRONGHOST) \n", 
	       ZCMB_ORIG, GENLC.REDSHIFT_CMB_SMEAR); fflush(stdout);
      }
    }
    
  }
  else {
    // no host galaxy redshift (and no SN redshift)
    SNHOSTGAL.ZSPEC            = -9.0 ;
    SNHOSTGAL.ZSPEC_ERR        = -9.0 ;

    if ( ZPHOT < 0.0 ) { // if no host photoz ...
      GENLC.REDSHIFT_FLAG        = REDSHIFT_FLAG_NONE ; 
      GENLC.REDSHIFT_CMB_SMEAR   = -9.0 ; 
      GENLC.REDSHIFT_HELIO_SMEAR = -9.0 ; 
      GENLC.REDSHIFT_SMEAR_ERR   = -9.0 ;
    } 
    else {
      // use host zPHOT  here if it is defined
      GENLC.REDSHIFT_FLAG        = REDSHIFT_FLAG_HOSTPHOT ; 
      GENLC.REDSHIFT_HELIO_SMEAR = ZPHOT ;
      GENLC.REDSHIFT_SMEAR_ERR   = ZPHOT_ERR ;
      GENLC.REDSHIFT_CMB_SMEAR = zhelio_zcmb_translator(ZPHOT,RA,DEC,eq, +1);
    }
  }

  return ;

} // end of  setz_unconfirmed



// **********************************
int gen_smearMag ( int epoch, int VBOSE) {

  // Dec 2011
  // Convert flux(ADU) and ZP into observed mag and error.

  double flux, flux_errstat, flux_tmp, zpt ;
  double mag, mag_err , mag_tmp, genmag ;

  //  char fnam[] = "gen_smearMag" ;

  // -------------- BEGIN --------------

  flux         = GENLC.flux[epoch];
  flux_errstat = GENLC.fluxerr_data[epoch];
  genmag       = GENLC.genmag_obs[epoch];
  zpt          = SIMLIB_OBS_GEN.ZPTADU[epoch] ;

  // init output
  GENLC.mag[epoch]      = NULLFLOAT ;     
  GENLC.mag_err[epoch]  = NULLFLOAT ;     

  if ( flux > 0.0  )
    { mag = -2.5 * log10(flux) + zpt ; }
  else
    { mag = MAG_NEGFLUX ; }


  // now get the mag error by subtracting 1 sigma from observed flux
  flux_tmp  =  flux - flux_errstat ; 

  if ( flux_tmp > 0.0  )
    { mag_tmp = -2.5 * log10(flux_tmp) + zpt ; }
  else
    { mag_tmp = MAG_NEGFLUX ; }

  mag_err = mag_tmp - mag ;


  // Jan 22, 2010: preserve mag=666 to indicate filter-problem
  if ( genmag > 600.0 )  
    { mag = MAG_UNDEFINED ; mag_err  = MAGERR_UNDEFINED; }
 
  GENLC.mag[epoch]      = mag ;
  GENLC.mag_err[epoch]  = mag_err ;

  if ( GENLC.npe_above_sat[epoch] > 0 ) {
    GENLC.mag[epoch]                   = MAG_SATURATE ;
    GENLC.mag_err[epoch]               = 5.0 ;
  }

  return SUCCESS ;

} // end of gen_smearMag


// ****************************************************
int npe_above_saturation ( int epoch, double flux_pe) {

  // Created jan 3 2018
  // For input epoch index and source flux (p.e.),
  // return number of photo-electrons above saturation
  // in central pixel. 
  //   npe < 0 ==> not saturated
  //   npe > 0 ==> saturated
  //
  // Used only if NPE_PIXEL_SATURATION key is given in the SIMLIB header.
  // 
  // Oct 28 2021: reutrn -999 if INPUTS.EXPOSURE_TIME>10
  //               (otherwise all epochs saturate and fail trigger)
  //
  int npe = -999 ; // default is not saturated
  int NEXPOSE      = SIMLIB_OBS_GEN.NEXPOSE[epoch] ;
  int LDMP ;

  double ccdgain    = SIMLIB_OBS_GEN.CCDGAIN[epoch] ;
  double skysig_adu = SIMLIB_OBS_GEN.SKYSIG[epoch] ;
  double psfsig1    = SIMLIB_OBS_GEN.PSFSIG1[epoch] ; // pixels
  double PIX        = SIMLIB_OBS_GEN.PIXSIZE[epoch];  // pix size, arcsec
  double PI         = 3.1415926535 ;

  double NPE_SAT = SIMLIB_GLOBAL_HEADER.NPE_PIXEL_SATURATE ;
  double skysig_pe, sky_pe, areaFrac, areaCor, AREAPIX, SQPSFSIG ;
  double fluxtot_pe, fluxtot_pe_avg;
  char fnam[] = "npe_above_saturation" ;

  // -------------- BEGIN ------------

  // if no saturation value is given, then never saturate
  if ( NPE_SAT > 999999999 ) { return(npe); }

  // skip saturation check if using large exposure_time (Oct 28 2021)     
  if ( INPUTS.EXPOSURE_TIME > 10.0 ) { return(npe); } //.xyz                         

  skysig_pe = (skysig_adu * ccdgain) ;  // convert ADU -> pe
  sky_pe    = skysig_pe * skysig_pe ;   // sky counts per pix, Npe

  // compute fraction of flux contained in 1 pixel at the
  // center of the PSF.  For pixel size << PSF, flux-fraction
  // is about PIXSIZE^2/2*PI*sigma^2. For safety, also include
  // Taylor-expansion correction in areaCor factor below.
  AREAPIX    = PIX*PIX ;
  SQPSFSIG   = psfsig1*psfsig1 * AREAPIX ; // arcSec
  areaCor    =   AREAPIX/(4.0*PI*SQPSFSIG);  // from Taylor expansion
  areaFrac   = ( AREAPIX/(2.0*PI*SQPSFSIG) ) * ( 1.0 - areaCor );

  // flux in central pixel
  fluxtot_pe = (flux_pe*areaFrac) + sky_pe ;
  
  // divide by number of exposues
  fluxtot_pe_avg = fluxtot_pe / (double)NEXPOSE ;

  // comput number of photoelectrons above saturation
  npe = (int)( fluxtot_pe_avg - NPE_SAT);


  LDMP = ( epoch <= -4) ;
  if ( LDMP ) {
    int ifilt_obs = GENLC.IFILT_OBS[epoch] ; 
    double genmag = GENLC.genmag_obs[epoch];
    printf("\n");
    printf(" xxx ------------------------------------ \n");
    printf(" xxx %s DUMP \n", fnam);
    printf(" xxx INPUTS: epoch=%d  flux_pe=%.1f  NEXPOSE=%d   NPE_SAT=%.0f\n",
	   epoch, flux_pe, NEXPOSE, NPE_SAT);
    printf(" xxx genmag = %.2f   ifiltobs=%d \n", genmag, ifilt_obs);
    printf(" xxx areaCor = %f  areaFrac = %f \n", areaCor, areaFrac );
    printf(" xxx sky,fluxtot, fluxtot_avg = %.1f, %.1f, %.1f pe \n", 
	   sky_pe, fluxtot_pe, fluxtot_pe_avg );
    printf(" xxx npe above saturation :  %d \n", npe ) ;
    if ( epoch ==4 ) { debugexit(fnam); }
  }


  return(npe);

} // end npe_above_saturation


// **************************************
void snlc_to_SNDATA(int FLAG) {


  // load SNDATA struct using info from GENLC structure.
  // FLAG = 0 => load everything.
  // FLAG = 1 => load header info only (for fits format)
  //
  // Set CUTBIT_MASK bits
  //
  // Mar 02 2018: scale SKYSIG_T to same ZP as SKYSIG_S
  //              Only affects output to data files.
  // Jun 01 2018: for LCLIB, set all redshifts and errors to zero
  //
  // Jul 21 2018: apply MW correction on fluxes (default is no cor)
  //
  // Dec 10 2018: load BYOSED info
  // Jul 20 2019: load strong lens info into SNDATA.SIM_SL_XXX
  // Feb 28 2021: load NEA
  // ---------------------------------------

  int PHOTFLAG_DETECT  = INPUTS_SEARCHEFF.PHOTFLAG_DETECT;
  int PHOTFLAG_TRIGGER = INPUTS_SEARCHEFF.PHOTFLAG_TRIGGER ;

  int epoch, ifilt, ifilt_obs, NFILT, NPAR, ipar, MSKTMP, istat, zFLAG ;
  double ZP_S, ZP_T, SKYSIG_T_scale, arg;
  double MCOR_MAP_MW, MCOR_TRUE_MW ;
  
  char ctmp[200], ccid[12], cfilt[2];
  char *cptr, *tmpName;
  char fnam[] = "snlc_to_SNDATA" ;

  // --------------- BEGIN -------------

  // always start with header info

  
  sprintf(SNDATA.SNANA_VERSION, "%s", SNANA_VERSION_CURRENT);
  sprintf(SNDATA.SURVEY_NAME,   "%s", SIMLIB_GLOBAL_HEADER.SURVEY_NAME );
  sprintf(SNDATA.SUBSURVEY_NAME,"%s", SIMLIB_HEADER.SUBSURVEY_NAME );  
  sprintf(SNDATA.DATATYPE,      "%s", DATATYPE_SIM_SNANA);

  SNDATA.SIMLIB_MSKOPT  = INPUTS.SIMLIB_MSKOPT ;
  sprintf(SNDATA.SIMLIB_FILE,    "%s", INPUTS.SIMLIB_FILE );
  sprintf(SNDATA.KCOR_FILE,      "%s", INPUTS.KCOR_FILE   );

  sprintf(SNDATA.HOSTLIB_FILE,   "%s", INPUTS.HOSTLIB_FILE );
  sprintf(SNDATA.SIM_MODEL_NAME, "%s", INPUTS.MODELNAME );
  sprintf(SNDATA_FILTER.LIST,    "%s", INPUTS.GENFILTERS );
  SNDATA.SIM_MODEL_INDEX    = INDEX_GENMODEL ;

  SNDATA.SIMOPT_FLUXERR     = INPUTS.FUDGEOPT_FLUXERR; // Oct 30 2015
  SNDATA.SIMOPT_MWCOLORLAW  = INPUTS.OPT_MWCOLORLAW ;
  SNDATA.SIMOPT_MWEBV       = INPUTS.OPT_MWEBV ;
  SNDATA.SIM_MWRV           = INPUTS.RV_MWCOLORLAW ;
  SNDATA.APPLYFLAG_MWEBV    = INPUTS.APPLYFLAG_MWEBV ;
    
  SNDATA.SIM_NOBS_UNDEFINED   = GENLC.NOBS_UNDEFINED;

  SNDATA.MJD_TRIGGER          = GENLC.MJD_TRIGGER ;
  SNDATA.MJD_DETECT_FIRST     = GENLC.MJD_DETECT_FIRST;
  SNDATA.MJD_DETECT_LAST      = GENLC.MJD_DETECT_LAST;

  SNDATA.NEA_PSF_UNIT         = SIMLIB_GLOBAL_HEADER.NEA_PSF_UNIT;

  // Jul 18,2011: write INPUTS GENFILTERS only
  // advantage: does not exceed Kcor bound of there are more than
  //    10 obs-filters defined in the SIMLIB
  // disadvantage:  cannot combine simulations that generate
  //   different filter-subsets; i.e., ugr and ugri
  //

  NFILT = INPUTS.NFILTDEF_OBS ;
  SNDATA_FILTER.NDEF = NFILT;
  for (ifilt=0; ifilt < NFILT; ifilt++ ) 
    { SNDATA_FILTER.MAP[ifilt] = INPUTS.IFILTMAP_OBS[ifilt]; }


  // load SIMSED info (for SIMSED model)
  if ( INDEX_GENMODEL == MODEL_SIMSED ) {
    NPAR = INPUTS.NPAR_SIMSED ;
    SNDATA.NPAR_SIMSED = NPAR ;
    for ( ipar=0; ipar < NPAR; ipar++ ) {

      sprintf(SNDATA.SIMSED_PARNAME[ipar],"%s",
	      INPUTS.PARNAME_SIMSED[ipar] );

      sprintf(SNDATA.SIMSED_KEYWORD[ipar],"%s(%s)" 
	      ,INPUTS.KEYWORD_SIMSED[ipar]
	      ,INPUTS.PARNAME_SIMSED[ipar] );

      SNDATA.SIMSED_PARVAL[ipar]  = GENLC.SIMSED_PARVAL[ipar];
    }
  } // model = MODEL_SIMSED


  // load PySEDMODEL info  (Dec 2018)
  if ( IS_PySEDMODEL ) {
    NPAR = Event_PySEDMODEL.NPAR ;
    SNDATA.NPAR_PySEDMODEL = NPAR ; 
    for ( ipar=0; ipar < NPAR; ipar++ ) {
      tmpName = Event_PySEDMODEL.PARNAME[ipar];
      sprintf(SNDATA.PySEDMODEL_PARNAME[ipar], "%s", tmpName ) ;
      sprintf(SNDATA.PySEDMODEL_KEYWORD[ipar], 
	      "%s_PARAM(%s)", INPUTS_PySEDMODEL.MODEL_NAME, tmpName); 
      SNDATA.PySEDMODEL_PARVAL[ipar]  = Event_PySEDMODEL.PARVAL[ipar];
    }
  } 


  // load LCLIB info 
  if ( INDEX_GENMODEL == MODEL_LCLIB ) {
    NPAR = LCLIB_INFO.NPAR_MODEL_STORE ; 
    SNDATA.NPAR_LCLIB = NPAR ; 
    for ( ipar=0; ipar < NPAR; ipar++ ) {
      tmpName = LCLIB_INFO.PARNAME_MODEL[ipar];
      sprintf(SNDATA.LCLIB_PARNAME[ipar], "%s", tmpName ) ;
      sprintf(SNDATA.LCLIB_KEYWORD[ipar],"LCLIB_PARAM(%s)", tmpName);
      SNDATA.LCLIB_PARVAL[ipar]  = LCLIB_EVENT.PARVAL_MODEL[ipar];
    }
  } // model = MODEL_LCLIB



  // load strong lens (SL) info
  SNDATA.SIM_SL_FLAG = INPUTS_STRONGLENS.USE_FLAG;
  if ( SNDATA.SIM_SL_FLAG ) {
    int IMGNUM = GENSL.IMGNUM ;
    SNDATA.SIM_SL_IDLENS    = GENSL.IDLENS;
    SNDATA.SIM_SL_zLENS     = GENSL.zLENS;
    SNDATA.SIM_SL_NIMG      = GENSL.NIMG ;
    SNDATA.SIM_SL_IMGNUM    = GENSL.IMGNUM ;
    if ( IMGNUM >= 0 ) {
      SNDATA.SIM_SL_TDELAY    = GENSL.TDELAY_LIST[IMGNUM] ;
      SNDATA.SIM_SL_MAGSHIFT  = GENSL.MAGSHIFT_LIST[IMGNUM] ;
      SNDATA.SIM_SL_XIMG      = GENSL.XIMG_LIST[IMGNUM] ;
      SNDATA.SIM_SL_YIMG      = GENSL.YIMG_LIST[IMGNUM] ;
    }
    else {
      SNDATA.SIM_SL_TDELAY    = 0.0 ;
      SNDATA.SIM_SL_MAGSHIFT  = 0.0 ;
      SNDATA.SIM_SL_XIMG      = 0.0 ;
      SNDATA.SIM_SL_YIMG      = 0.0 ;
    }
  }

  hostgal_to_SNDATA(FLAG,0);  // for header init only.

  if ( INPUTS.GENGAUSS_SALT2ALPHA.NGRID > 1 ) 
    { SNDATA.SIM_BIASCOR_MASK += 1; }
  if ( INPUTS.GENGAUSS_SALT2BETA.NGRID > 1 ) 
    { SNDATA.SIM_BIASCOR_MASK += 2; }


  // #####################################

  if ( FLAG == 1 ) { return ; }

  // #####################################

  VERSION_INFO.N_SNLC       = NGENLC_WRITE ;
  VERSION_INFO.GENFRAME_SIM = GENFRAME_OPT ;


  // do smearing up front to see what's going on.
 
  SNDATA.WRFLAG_BLINDTEST = (WRFLAG_BLINDTEST>0) ;
  SNDATA.WRFLAG_PHOTPROB  = INPUTS_SEARCHEFF.NMAP_PHOTPROB > 0 ;
  SNDATA.WRFLAG_SKYSIG_T  = SIMLIB_TEMPLATE.USEFLAG;

  if ( GENLC.NEPOCH >= MXEPOCH ) {
    print_preAbort_banner(fnam);
    printf("\t LIBID=%d  z=%.3f  PEAKMJD=%.3f\n", 
	   GENLC.SIMLIB_ID, GENLC.REDSHIFT_CMB, GENLC.PEAKMJD );

    sprintf(c1err,"NEPOCH=%d exceeds bound of MXEPOCH=%d",
	    GENLC.NEPOCH, MXEPOCH );
    sprintf(c2err,"Check TREST range, PEAKMJD range and SIMLIB range");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  SNDATA.NEPOCH        = GENLC.NEPOCH ;
  SNDATA.RA            = GENLC.RA ;
  SNDATA.DEC           = GENLC.DEC ;
  SNDATA.SNTYPE        = GENLC.SNTYPE;
  SNDATA.SIM_TYPE_INDEX  = GENLC.SIMTYPE ;
  sprintf(SNDATA.SIM_TYPE_NAME, "%s", GENLC.SNTYPE_NAME );
  
  if ( WRFLAG_BLINDTEST ) 
    { SNDATA.FAKE  = FAKEFLAG_LCSIM_BLINDTEST ; }
  else
    { SNDATA.FAKE  = FAKEFLAG_LCSIM ; }

  if ( WRFLAG_CIDRAN == 0 ) 
    { SNDATA.CID  = GENLC.CID ; }
  else 
    { SNDATA.CID  = GENLC.CIDRAN ; }
  

  sprintf(SNDATA.CCID,      "%d", SNDATA.CID ) ;
  sprintf(SNDATA.IAUC_NAME, "%s", "NULL" );

  SNDATA.SUBSAMPLE_INDEX = GENLC.SUBSAMPLE_INDEX ;

  // Hard-wire legacy search type used by the SDSS
  // Users really should set "SNTYPE_Ia:  120" in the SDSS sim-input file.
  if ( strcmp(SNDATA.SURVEY_NAME,"SDSS") == 0  ) {
    SNDATA.SEARCH_TYPE   = 120 ;  // type Ia
    if ( INDEX_GENMODEL == MODEL_NON1ASED ) { SNDATA.SEARCH_TYPE  = 110 ; }
    if ( INDEX_GENMODEL == MODEL_NON1AGRID) { SNDATA.SEARCH_TYPE  = 110 ; }
    if ( SNDATA.SNTYPE < 0 ) 
      { SNDATA.SNTYPE = SNDATA.SEARCH_TYPE ; }
  }

  SNDATA.NOBS                = GENLC.NOBS ; // 5/26/2011

  SNDATA.REDSHIFT_FINAL      = GENLC.REDSHIFT_CMB_SMEAR ;
  SNDATA.REDSHIFT_HELIO      = GENLC.REDSHIFT_HELIO_SMEAR ;
  SNDATA.REDSHIFT_FINAL_ERR  = GENLC.REDSHIFT_SMEAR_ERR ;
  SNDATA.REDSHIFT_HELIO_ERR  = GENLC.REDSHIFT_SMEAR_ERR ;
  SNDATA.VPEC                = GENLC.VPEC_SMEAR;
  SNDATA.VPEC_ERR            = INPUTS.VPEC_ERR;

  SNDATA.SIM_SEARCHEFF_MASK  = GENLC.SEARCHEFF_MASK ;
  SNDATA.SIM_SEARCHEFF_SPEC  = GENLC.SEARCHEFF_SPEC ;
  SNDATA.SIM_SEARCHEFF_zHOST = GENLC.SEARCHEFF_zHOST ;

  // assign photoz to REDSHIFT_FINAL if Zspec error is >= 1
  if ( GENLC.REDSHIFT_SMEAR_ERR > 0.999  ) {
    SNDATA.REDSHIFT_FINAL      = SNHOSTGAL.ZPHOT ;
    SNDATA.REDSHIFT_FINAL_ERR  = SNHOSTGAL.ZPHOT_ERR ;

    SNDATA.REDSHIFT_HELIO      = SNHOSTGAL.ZPHOT ;
    SNDATA.REDSHIFT_HELIO_ERR  = SNHOSTGAL.ZPHOT_ERR ;
  }

  if ( INPUTS.NGENTOT_LC > 1000000 || 
       INPUTS.NGEN_LC    > 1000000 || 
       WRFLAG_CIDRAN     > 0 ) {
    sprintf(ccid, "%8.8d", SNDATA.CID );
  }
  else {
    sprintf(ccid, "%6.6d", SNDATA.CID );
  }
  sprintf(SNDATA.snfile_output, "%s_SN%s.DAT",
	  INPUTS.GENPREFIX, ccid );
  sprintf(SNDATA.SNFILE_OUTPUT, "%s/%s",
	  PATH_SNDATA_SIM, SNDATA.snfile_output );

  sprintf(SNDATA.MAGREF, "%s", GENLC.primary );

  // load SIM_XXX info

  if ( INPUTS.NON1A_MODELFLAG > 0 ) 
    { sprintf(ctmp,"MODEL = %s", GENLC.SNTEMPLATE ); }
  else
    {  sprintf(ctmp,"MODEL = %s", INPUTS.MODELNAME ); }
  sprintf(SNDATA.SIM_COMMENT, "SN Type = %s , %s", GENLC.SNTYPE_NAME, ctmp );

  SNDATA.SIM_LIBID           = GENLC.SIMLIB_ID;
  SNDATA.SIM_NGEN_LIBID      = GENLC.NGEN_SIMLIB_ID;
  SNDATA.SIM_REDSHIFT_CMB    = GENLC.REDSHIFT_CMB ;
  SNDATA.SIM_REDSHIFT_HELIO  = GENLC.REDSHIFT_HELIO ;
  SNDATA.SIM_REDSHIFT_HOST   = GENLC.REDSHIFT_HOST ; // Jan 2016

  // 4.19.2019: store redshift flag to indicate where redshift is from
  zFLAG = GENLC.REDSHIFT_FLAG;
  if ( zFLAG >= 0 ) {
    SNDATA.SIM_REDSHIFT_FLAG   = zFLAG ;
    sprintf(SNDATA.SIM_REDSHIFT_COMMENT,"%s", STRING_REDSHIFT_FLAG[zFLAG] );
  }

  SNDATA.SIM_VPEC         = GENLC.VPEC ;
  SNDATA.SIM_DLMU         = GENLC.DLMU ; 
  SNDATA.SIM_LENSDMU      = GENLC.LENSDMU ;
  SNDATA.SIM_RA           = GENLC.RA ;
  SNDATA.SIM_DEC          = GENLC.DEC ;
  SNDATA.SIM_AVTAU        = GENLC.AVTAU ;
  SNDATA.SIM_MWEBV        = GENLC.MWEBV_SMEAR ; // use smeared MWEBV
  SNDATA.SIM_PEAKMJD      = GENLC.PEAKMJD ;

  SNDATA.PIXSIZE          = SIMLIB_OBS_GEN.PIXSIZE[1];
  SNDATA.CCDNUM[0]        = SIMLIB_HEADER.CCDNUM ;
  SNDATA.SIM_AV           = GENLC.AV ;
  SNDATA.SIM_RV           = GENLC.RV ;

  if ( INDEX_GENMODEL == MODEL_STRETCH ) 
    {  SNDATA.SIM_STRETCH  = GENLC.STRETCH ; }

  if ( INDEX_GENMODEL == MODEL_MLCS2k2 ) 
    { SNDATA.SIM_DELTA    = GENLC.DELTA ; }

  if ( INDEX_GENMODEL == MODEL_SNOOPY ) 
    { SNDATA.SIM_STRETCH = GENLC.STRETCH ; }

  if ( INDEX_GENMODEL == MODEL_S11DM15 ) 
    { SNDATA.SIM_DM15 = GENLC.DM15 ; }


  SNDATA.SIM_SALT2x0          = GENLC.SALT2x0 ; 
  SNDATA.SIM_SALT2x1          = GENLC.SALT2x1 ; 
  SNDATA.SIM_SALT2c           = GENLC.SALT2c  ; 
  SNDATA.SIM_SALT2mB          = GENLC.SALT2mB ; 
  SNDATA.SIM_SALT2alpha       = GENLC.SALT2alpha ; 
  SNDATA.SIM_SALT2beta        = GENLC.SALT2beta ; 
  SNDATA.SIM_SALT2gammaDM     = GENLC.SALT2gammaDM ;
  SNDATA.SIM_TEMPLATE_INDEX   = GENLC.TEMPLATE_INDEX ; 
  SNDATA.SIM_MAGSMEAR_COH     = GENLC.MAGSMEAR_COH[0] ;
  SNDATA.SIM_RISETIME_SHIFT   = GENLC.RISETIME_SHIFT ;
  SNDATA.SIM_FALLTIME_SHIFT   = GENLC.FALLTIME_SHIFT ;
  SNDATA.SIM_TRESTMIN         = GENLC.TRESTMIN ;
  SNDATA.SIM_TRESTMAX         = GENLC.TRESTMAX ;

  // set GALID here in case HOSTLIB_USE=0 in hostgal_to_SNDATA
  if ( SNHOSTGAL.GALID > 0 ) { 
    SNDATA.HOSTGAL_NMATCH[0] = 1 ; 
    SNDATA.HOSTGAL_NMATCH[1] = 1 ; 
  }
  SNDATA.HOSTGAL_OBJID[0]   = SNHOSTGAL.GALID ;
  
  
  // set HOSTLIB variables
  ifilt_obs=0 ;  hostgal_to_SNDATA(FLAG,ifilt_obs);

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];  

    MCOR_TRUE_MW  = GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs]; 
    MCOR_MAP_MW   = GENLC.MAGCOR_MWEBV_MAP[ifilt_obs];

    SNDATA.SIM_PEAKMAG[ifilt_obs]     
      = (float)( GENLC.peakmag_obs[ifilt_obs] - MCOR_TRUE_MW) ;

    SNDATA.SIM_TEMPLATEMAG[ifilt_obs]       // LCLIB source mag in template
      = (float)GENLC.genmag_obs_template[ifilt_obs] - MCOR_MAP_MW ;

    SNDATA.SIM_EXPOSURE_TIME[ifilt_obs]  = 
      INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] ; 


    hostgal_to_SNDATA(FLAG,ifilt_obs);

  } // end of ifilt if-block



  SNDATA.MWEBV      = GENLC.MWEBV;      // reported MWEBV from Shlagel map
  SNDATA.MWEBV_ERR  = GENLC.MWEBV_ERR ; // reported error
  SNDATA.APPLYFLAG_MWEBV = INPUTS.APPLYFLAG_MWEBV ;
    
  // load NEWMJD info
  SNDATA.NEWMJD = GENLC.NEWMJD ;
  for ( epoch = 1; epoch <= GENLC.NEWMJD ; epoch++ ) {
    SNDATA.EPOCH_RANGE_NEWMJD[epoch][0] = 
      GENLC.EPOCH_RANGE_NEWMJD[epoch][0] ;
    SNDATA.EPOCH_RANGE_NEWMJD[epoch][1] = 
      GENLC.EPOCH_RANGE_NEWMJD[epoch][1] ;
  }


  // covmat-scatter info
  int iscat ;
  if ( INPUTS.NCOVMAT_SCATTER > 0 ) {
    SNDATA.SIMFLAG_COVMAT_SCATTER = 1 ;
    for ( iscat=0; iscat < 3 ; iscat++ ) {
      SNDATA.SIM_COVMAT_SCATTER[iscat] = 
	GENLC.COVMAT_SCATTER[iscat] ;
      cptr = SNDATA.SIM_COVMAT_SCATTER_NAME[iscat];
      sprintf(cptr,"%s", GENLC.COVMAT_SCATTER_NAME[iscat] );
    }
  }

  // load epoch info


  for ( epoch = 1; epoch <= GENLC.NEPOCH ; epoch++ ) {

    SNDATA.OBSFLAG_WRITE[epoch] = GENLC.OBSFLAG_WRITE[epoch] ;

    if ( !GENLC.OBSFLAG_GEN[epoch] ) { continue; }

    ifilt_obs    = GENLC.IFILT_OBS[epoch];
      
    SNDATA.FILTINDX[epoch]  = ifilt_obs ;
    sprintf(SNDATA.FILTCHAR[epoch], "%c",  FILTERSTRING[ifilt_obs] );

    MCOR_TRUE_MW  = GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs]; 
    MCOR_MAP_MW   = GENLC.MAGCOR_MWEBV_MAP[ifilt_obs];

    SNDATA.SIMEPOCH_TREST[epoch]  = GENLC.epoch_rest[epoch] ;
    SNDATA.SIMEPOCH_TOBS[epoch]   = GENLC.epoch_obs[epoch] ;
    SNDATA.SIMEPOCH_MAG[epoch]    = GENLC.genmag_obs[epoch] - MCOR_TRUE_MW ;
    SNDATA.SIMEPOCH_MODELMAGERR[epoch] = GENLC.generr_rest[epoch] ;
    SNDATA.SIMEPOCH_MAGSMEAR[epoch] = GENLC.magsmear8[epoch] ;
    SNDATA.SIMEPOCH_FLUXCAL_HOSTERR[epoch] = GENLC.NOISE_HOSTGAL_PHOT[epoch];

    SNDATA.MJD[epoch]          = GENLC.MJD[epoch];

    sprintf(SNDATA.TELESCOPE[epoch], "%s", GENLC.TELESCOPE[epoch] );
    sprintf(SNDATA.FIELDNAME[epoch], "%s", GENLC.FIELDNAME[epoch] );

    SNDATA.GAIN[epoch]      =  SIMLIB_OBS_GEN.CCDGAIN[epoch] ;
    SNDATA.READNOISE[epoch] =  SIMLIB_OBS_GEN.READNOISE[epoch] ;

    // -------------------------------
    // PHOTFLAG_DETECT is the optional mask to set for each detection;
    // detectFlag = 0,1 for undetected/detected. Note the epoch offset
    // bewteen the SNDATA and SEARCHEFF_DATA structures
    MSKTMP = SEARCHEFF_DATA.detectFlag[epoch-1];
    if ( (MSKTMP & 1)>0  &&  PHOTFLAG_DETECT>0 ) 
      { SNDATA.PHOTFLAG[epoch] += PHOTFLAG_DETECT ; } 

    if ( (MSKTMP & 2)>0  &&  PHOTFLAG_TRIGGER>0 ) 
      { SNDATA.PHOTFLAG[epoch] += PHOTFLAG_TRIGGER ; } 


    // Jan 17 2018: check appended observations
    int APP = SIMLIB_OBS_GEN.APPEND_PHOTFLAG[epoch] ;
    if ( APP > 0 ) 
      { SNDATA.PHOTFLAG[epoch] += APP ; }

    if ( GENLC.npe_above_sat[epoch] > 0 ) 
      { SNDATA.PHOTFLAG[epoch] += SIMLIB_GLOBAL_HEADER.PHOTFLAG_SATURATE; }

    if ( GENLC.IEPOCH_SNRMAX_GLOBAL == epoch ) 
      { SNDATA.PHOTFLAG[epoch] += SIMLIB_GLOBAL_HEADER.PHOTFLAG_SNRMAX ; }

    if ( GENLC.IEPOCH_NEARPEAK == epoch ) 
      { SNDATA.PHOTFLAG[epoch] += SIMLIB_GLOBAL_HEADER.PHOTFLAG_NEARPEAK ; }

    // load photProb
    SNDATA.PHOTPROB[epoch] = SEARCHEFF_DATA.PHOTPROB[epoch-1] ;

    // - - - - - - - - -  -
    double diff = SNDATA.MJD[epoch] - SEARCHEFF_DATA.MJD[epoch-1] ;
    if ( diff != 0.0 ) {
      sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
      print_preAbort_banner(fnam);
      printf("   SNDATA.MJD[%d]=%.4f   SEARCHEFF_DATA.MJD[%d]=%.4f\n",
	     epoch, SNDATA.MJD[epoch], 
	     epoch-1, SEARCHEFF_DATA.MJD[epoch-1] );
      sprintf(c1err,"Index problem with SEARCHEFF_DATA struct.");
      sprintf(c2err,"CID=%d  epoch=%d of %d,  %s-band, MJD-diff=%f", 
	      GENLC.CID, epoch, GENLC.NEPOCH, cfilt, diff);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
      
    // ------------------------
    // store conditions
    SNDATA.PSF_SIG1[epoch]  = SIMLIB_OBS_GEN.PSFSIG1[epoch] ;
    SNDATA.PSF_SIG2[epoch]  = SIMLIB_OBS_GEN.PSFSIG2[epoch] ;
    SNDATA.PSF_RATIO[epoch] = SIMLIB_OBS_GEN.PSFRATIO[epoch] ;

    if ( SNDATA.NEA_PSF_UNIT ) {
      SNDATA.PSF_NEA[epoch] = SIMLIB_OBS_GEN.NEA[epoch] ; // Feb 2021
    }

    // Feb 23, 2012: store search-run zpt instead of template run  zpt
    SNDATA.ZEROPT[epoch]      = SIMLIB_OBS_GEN.ZPTADU[epoch] ;
    SNDATA.ZEROPT_ERR[epoch]  = SIMLIB_OBS_GEN.ZPTERR[epoch] ;
    SNDATA.ZEROPT_SIG[epoch]  = SIMLIB_OBS_GEN.ZPTERR[epoch] ;

    // mar 18 2018: store SNR at fixed mag to monitor data quality 
    SNDATA.SIMEPOCH_SNRMON[epoch] = GENLC.SNR_MON[epoch];

    // ----------------

    SNDATA.NPE_ABOVE_SAT[epoch]     = GENLC.npe_above_sat[epoch];
    SNDATA.FLUX[epoch]              = GENLC.flux[epoch];
    SNDATA.FLUX_ERRTOT[epoch]       = GENLC.fluxerr_data[epoch]; 
    SNDATA.FLUX_ERRTEMPLATE[epoch]  = GENLC.template_err[epoch] ;

    SNDATA.MAG[epoch]          = GENLC.mag[epoch];
    SNDATA.MAG_ERRPLUS[epoch]  = GENLC.mag_err[epoch];
    SNDATA.MAG_ERRMINUS[epoch] = GENLC.mag_err[epoch];

    SNDATA.SIMEPOCH_WARPCOLVAL[epoch]  = GENLC.warpcolval8[epoch] ;

    sprintf(SNDATA.SIMEPOCH_WARPCOLNAM[epoch], "%s",
	    GENLC.warpcolnam[epoch] ) ;

    SNDATA.SIMEPOCH_AVWARP[epoch]   = GENLC.AVwarp8[epoch] ;
    SNDATA.SIMEPOCH_KCORVAL[epoch]  = GENLC.kcorval8[epoch] ;
    sprintf(SNDATA.SIMEPOCH_KCORNAM[epoch], "%s",
	    GENLC.kcornam[epoch] ) ;
    
    // --> fill SNDATA.FLUXCAL
    int OPT_ZPERR = 2; // --> do NOT add ZP_sig since it's already added
    istat = fluxcal_SNDATA ( epoch, "log10", OPT_ZPERR ) ; 
    
    // check option to make MWEBV-analysis correction on the
    // reported flux, flux-error and mag
    if ( INPUTS.APPLYFLAG_MWEBV ) { MWEBVfluxCor_to_SNDATA(epoch); }

    // store SKY sigmas in ADU/pixel
    // Mar 2 2018: scale output SIG_T to same ZP as search
    ZP_S = SIMLIB_OBS_GEN.ZPTADU[epoch] ;
    ZP_T = SIMLIB_OBS_GEN.TEMPLATE_ZPT[epoch] ;
    arg  = 0.4*(ZP_S-ZP_T);
    SKYSIG_T_scale = pow(TEN,arg);  // Mar 2 2018

    SNDATA.SKY_SIG[epoch]    = SIMLIB_OBS_GEN.SKYSIG[epoch] ;
    SNDATA.SKY_SIG_T[epoch]  = 
      SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[epoch] * SKYSIG_T_scale ;

  } // end epoch loop


  // May 2019: 
  // estimate PEAKMJD after all of the FLUXCAL[ERR] are  evaluated.
  GENLC.PEAKMJD_SMEAR   = gen_peakmjd_smear(); 
  SNDATA.SEARCH_PEAKMJD = GENLC.PEAKMJD_SMEAR;

  return ;

}  // end of snlc_to_SNDATA


// **************************************************
void MWEBVfluxCor_to_SNDATA(int epoch) {

  // Correct FLUXCAL for MWEBV: 
  // i.e., FLUXCAL -> FLUXCAL* ( 1 + Delta )
  // and also increase FLUXCAL_ERRTOT.

  int    ifilt_obs     = GENLC.IFILT_OBS[epoch];
  double FCOR_MAP_MW   = GENLC.FLUXCOR_MWEBV_MAP[ifilt_obs] ;
  double MCOR_MAP_MW   = GENLC.MAGCOR_MWEBV_MAP[ifilt_obs]  ;
  //  double MCOR_TRUE_MW  = GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs] ;
  
  double FLUX    = (double)SNDATA.FLUXCAL[epoch] ;
  double FLUXERR = (double)SNDATA.FLUXCAL_ERRTOT[epoch] ;
  double FLUXERR_ORIG = FLUXERR;
  double MAG     = (double)SNDATA.MAG[epoch] ;
  double ratio_MWEBV, FLUXERR_ADD, arg;

  double   MAG_T  = SNDATA.SIM_TEMPLATEMAG[ifilt_obs] ;
  double   FLUX_T = pow(10.0, 0.4*(27.5-MAG_T) ); 
  int LDMP = (MAG_T < 12 && ifilt_obs == -3);
  char fnam[] = "MWEBVfluxCor_to_SNDATA" ;

  // --------------- BEGIN --------------

  // bail on saturated epochs
  if( GENLC.npe_above_sat[epoch] > 0 ) { return ; }

  FLUX      *= FCOR_MAP_MW ;
  FLUXERR   *= FCOR_MAP_MW ;
  MAG       -= MCOR_MAP_MW ; 
    
  // add extra error based on MWEBV uncertainty
  // here we use approximation with ratios.
  ratio_MWEBV = GENLC.MWEBV_ERR / GENLC.MWEBV_SMEAR ;
  arg         = 0.4*(MCOR_MAP_MW*ratio_MWEBV);
  FLUXERR_ADD = (pow(10.0,arg)-1.0) * fabs(FLUX) ;
  FLUXERR     = sqrt(FLUXERR*FLUXERR + FLUXERR_ADD*FLUXERR_ADD);

  // update global FLUX and MAG
  SNDATA.FLUXCAL[epoch]         = (float)FLUX ;
  SNDATA.FLUXCAL_ERRTOT[epoch]  = (float)FLUXERR ;
  SNDATA.MAG[epoch]             = (float)MAG ;

  if ( LDMP ) {
    printf(" xxx ----------------------------------------- \n");
    printf(" xxx %s DUMP for CID=%d ifiltobs=%d xxx \n", 
	   fnam, GENLC.CID, ifilt_obs);
    printf(" xxx MAG=%.3f  FLUX=%le   MWEBV(meas)=%.3f +- %.3f \n",	   
	   MAG, FLUX, GENLC.MWEBV_SMEAR, GENLC.MWEBV_ERR);
    printf(" xxx MAG_T=%.3f   FLUX_T=%le \n" , MAG_T, FLUX_T );
    printf(" xxx FCOR_MAP = %.4f   MCOR_MAP=%.4f\n", 
	   FCOR_MAP_MW, MCOR_MAP_MW);
    printf(" xxx FLUXERR(orig)  = %.2f (%le x FLUX) \n",
	   FLUXERR_ORIG, FLUXERR_ORIG/FLUX );
    printf(" xxx FLUXERR(add)   = %.2f (%le x FLUX) \n",
	   FLUXERR_ADD, FLUXERR_ADD/FLUX );
    printf(" xxx FLUXERR(final) = %.2f (%le x FLUX) \n",
	   FLUXERR, FLUXERR_ADD/FLUX );
  }

  return;

} // end MWEBVfluxCor_to_SNDATA

// **************************************************
void hostgal_to_SNDATA(int IFLAG, int ifilt_obs) {

  // Mar 2012
  
  // IFLAG = 0 --> nominal loading
  // IFLAG = 1 --> load header info only (for FITS-file init)

  // Compute and global hostgal properties (i.e., epoch-independent)
  // and load SNDATA structure,
  //   SNDATA.SIM_GALFRAC[ifilt_obs]      => reference quantity
  //   SNDATA.HOSTGAL_SB_FLUXCAL[ifilt_obs]  => data-like quantity
  //   SNDATA.HOSTGAL_USEMASK
  //   SNDATA.HOSTGAL_MAG[ifilt_obs]   ==> data-like quanity
  // where SB = surface brightness.
  //
  // Dec 17, 2012: fill HOSTGAL_NFILT_MAGOBS and others with ifilt_obs=0
  // Feb 12, 2014: fill SNDATA.SIM_HOSTGAL_xxx, and add IFLAG arg
  // Jun 02, 2018: load zphot info for LCLIB
  // Jan 29 2020: USE_REFACTOR -> true to get multiple hosts
  // Sep 09 2021: 
  //   + Fill nbr dimension of SNDATA.SIM_HOSTLIB_PARVAL Alex Gagliano


  int    NPAR, ipar, nbr, OVP, ifilt, NMATCH, m ;
  double psfsig, mag_GAL, mag_SN, mag_dif, fgal ;
  char  *name ;
  char fnam[] = "hostgal_to_SNDATA" ;

  // --------------- BEGIN ------------

  // for LCLIB, load photo-z info and return
  if ( ifilt_obs == 0 &&  INDEX_GENMODEL==MODEL_LCLIB ) {
    SNDATA.HOSTGAL_SPECZ[0]          = SNHOSTGAL.ZSPEC ;
    SNDATA.HOSTGAL_SPECZ_ERR[0]      = SNHOSTGAL.ZSPEC_ERR ;
    SNDATA.HOSTGAL_PHOTOZ[0]         = SNHOSTGAL.ZPHOT ;
    SNDATA.HOSTGAL_PHOTOZ_ERR[0]     = SNHOSTGAL.ZPHOT_ERR ;
    return ;
  }


  if ( INPUTS.HOSTLIB_USE == 0 ) { return ; }

  ifilt = GENLC.IFILTINVMAP_OBS[ifilt_obs] ; // sparse index

  if ( IFLAG == 1 ) {

    // one-time init before output files are initialized.
    // set counters and par-names.

    SNDATA.SIM_HOSTLIB_MSKOPT = 
      INPUTS.HOSTLIB_MSKOPT ; // needed in sntools_fitsio

    NPAR = HOSTLIB_OUTVAR_EXTRA.NOUT ;
    SNDATA.NPAR_SIM_HOSTLIB = NPAR ;
    for(ipar=0; ipar < NPAR ; ipar++ ) {
      name = HOSTLIB_OUTVAR_EXTRA.NAME[ipar] ;
      sprintf(SNDATA.SIM_HOSTLIB_PARNAME[ipar],"%s", name);

      // set key name for ascii output
      sprintf(SNDATA.SIM_HOSTLIB_KEYWORD[ipar],"SIM_HOSTLIB(%s)", name);
    }
    SNDATA.HOSTLIB_NFILT_MAGOBS   = HOSTLIB.NFILT_MAGOBS ;

    return ;
  }  // end IFLAG==1


  NMATCH = SNHOSTGAL.NNBR ;
  if ( NMATCH > MXHOSTGAL ) { NMATCH = MXHOSTGAL; }

  if ( ifilt_obs == 0 ) {

    // NEW(Nov 2019): test multiple host matches with NBR_LIST in HOSTLIB
    SNDATA.HOSTGAL_NMATCH[0] = SNDATA.HOSTGAL_NMATCH[1] = NMATCH ;
    for(m=0; m < NMATCH; m++ ) {
      SNDATA.HOSTGAL_OBJID[m]      = SNHOSTGAL_DDLR_SORT[m].GALID;
      SNDATA.HOSTGAL_PHOTOZ[m]     = SNHOSTGAL_DDLR_SORT[m].ZPHOT;
      SNDATA.HOSTGAL_PHOTOZ_ERR[m] = SNHOSTGAL_DDLR_SORT[m].ZPHOT_ERR;
      
      if ( SNHOSTGAL_DDLR_SORT[m].TRUE_MATCH == true ) {
	SNDATA.HOSTGAL_SPECZ[m]      = SNHOSTGAL.ZSPEC ;
	SNDATA.HOSTGAL_SPECZ_ERR[m]  = SNHOSTGAL.ZSPEC_ERR ;
      }
      else {
	SNDATA.HOSTGAL_SPECZ[m]      = SNHOSTGAL_DDLR_SORT[m].ZSPEC;
	SNDATA.HOSTGAL_SPECZ_ERR[m]  = SNHOSTGAL_DDLR_SORT[m].ZSPEC_ERR;
      }
      
      SNDATA.HOSTGAL_RA[m]           = SNHOSTGAL_DDLR_SORT[m].RA ;
      SNDATA.HOSTGAL_DEC[m]          = SNHOSTGAL_DDLR_SORT[m].DEC ;
      SNDATA.HOSTGAL_DDLR[m]         = SNHOSTGAL_DDLR_SORT[m].DDLR ;
      SNDATA.HOSTGAL_SNSEP[m]        = SNHOSTGAL_DDLR_SORT[m].SNSEP ;
      SNDATA.HOSTGAL_LOGMASS_TRUE[m] = SNHOSTGAL_DDLR_SORT[m].LOGMASS_TRUE;
      SNDATA.HOSTGAL_LOGMASS_OBS[m]  = SNHOSTGAL_DDLR_SORT[m].LOGMASS_OBS ;
      SNDATA.HOSTGAL_LOGMASS_ERR[m]  = SNHOSTGAL_DDLR_SORT[m].LOGMASS_ERR ;
      // Added for LSST but may be of more general use
      // Alex Gagliano 09/2021
      SNDATA.HOSTGAL_OBJID2[m]       = SNHOSTGAL_DDLR_SORT[m].GALID2;
      SNDATA.HOSTGAL_ELLIPTICITY[m]  = SNHOSTGAL_DDLR_SORT[m].ELLIPTICITY;
      SNDATA.HOSTGAL_SQRADIUS[m]     = SNHOSTGAL_DDLR_SORT[m].SQRADIUS;

      SNDATA.HOSTGAL_OBJID_UNIQUE[m] = SNHOSTGAL_DDLR_SORT[m].GALID_UNIQUE;
    }
    
  
    NPAR = SNDATA.NPAR_SIM_HOSTLIB ;
    for(nbr=0; nbr < MXHOSTGAL; nbr++){  // loop over neighbors 
      for(ipar=0; ipar < NPAR ; ipar++ ) {
	SNDATA.SIM_HOSTLIB_PARVAL[ipar][nbr] = 
	  HOSTLIB_OUTVAR_EXTRA.VALUE[ipar][nbr] ;
      }
    }
    SNDATA.SIM_HOSTLIB_GALID = SNHOSTGAL.GALID; // store true GALID, Feb 2020

    return ;

  } // end of ifilt_obs==0
  

  // transfer total host mag (Feb 2013)

  SNDATA.HOSTGAL_MAG[0][ifilt] = (float)SNHOSTGAL.GALMAG[ifilt_obs][0]; 
  SNDATA.HOSTGAL_USEMASK |= 1 ; // flag to write host mag


  for(m=0; m < NMATCH; m++ ) {
    SNDATA.HOSTGAL_MAG[m][ifilt] = 
      (float)SNHOSTGAL_DDLR_SORT[m].MAG[ifilt_obs] ;
    SNDATA.HOSTGAL_MAGERR[m][ifilt] =
      (float)SNHOSTGAL_DDLR_SORT[m].MAG_ERR[ifilt_obs] ;
  }

  SNDATA.HOSTGAL_SB_FLUXCAL[ifilt] = (float)SNHOSTGAL.SB_FLUXCAL[ifilt_obs];
  SNDATA.HOSTGAL_SB_MAG[ifilt]     = (float)SNHOSTGAL.SB_MAG[ifilt_obs];

  OVP = (INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_PHOT) ;
  if ( OVP > 0 ) {
    SNDATA.HOSTGAL_USEMASK |= 4 ; // flag to write surface brightness
    psfsig   = 1./2.355 ;     // typical PSF in arcsec
    mag_GAL  = interp_GALMAG_HOSTLIB(ifilt_obs,psfsig );
    mag_SN   = (double)SNDATA.SIM_PEAKMAG[ifilt_obs] ;
    mag_dif  = mag_GAL - mag_SN ;
    if ( mag_SN < 98.0 ) {
      // gal/SNpeak flux-fraction in 1" aperture  
      fgal  = pow(10.0,-0.4*mag_dif); 
    }
    else {
      fgal = -9.0 ; // undefined
    }
    SNDATA.SIM_GALFRAC[ifilt_obs] = (float)fgal; 

  }

  return ;

} // end of hostgal_to_SNDATA


// **********************************
void genmag_boost(void) {

  /****
       Start with GENLC.*rest* quantities, and apply generic
       transformations to get GENLC.*obs* quantities.
       Apply 
       - local extinction using AV
       - K-correction from rest to obs frame

       Note that rest-frame mags are not necessarity in SDSS system.
       K-correction transforms from rest-system to SDSS mags.

     Dec 27 2019: skip AV part if AV < 1E-9 rather than if AV==0.

  *****/

  int     OPT_SNXT = INPUTS.OPT_SNXT ;
  double  AVwarp[4], AV, RV, z, Trest,  x[10];
  double  mag[4], lamdif[4], mag_obs, kcor    ;

  int ifilt_obs, ifilt_rest1, ifilt_rest2, ifilt_rest3, epoch, NZ ;
  int ifilt_rest_tmp[4] ;   

  char fnam[] = "genmag_boost" ;

  // ------------------ BEGIN ------------

  // do NOT apply extinction if AV=0

  if ( GENLC.AV  < 1.0E-9 ) { goto KCOR ; }

#ifdef MODELGRID_GEN
  // skip boost for SNOOPY model with just 1 logz-bin
  NZ = GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_LOGZ];
  if ( NZ == 1 &&  INDEX_GENMODEL == MODEL_SNOOPY ) { return ; }
#endif


  // apply extinction in rest frame mags/filters
  AV = (double)GENLC.AV ;
  RV = (double)GENLC.RV ;
  z  = 0.0;

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {  

    Trest = GENLC.epoch_rest[epoch]; 

    ifilt_obs   = GENLC.IFILT_OBS[epoch];

    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

    ifilt_rest1 = GENLC.IFILTMAP_REST1[ifilt_obs];
    ifilt_rest2 = GENLC.IFILTMAP_REST2[ifilt_obs];
    ifilt_rest3 = GENLC.IFILTMAP_REST3[ifilt_obs];     
       
    // start with nearest filter.
    x[1] = get_snxt8__( &OPT_SNXT, &ifilt_rest1, &Trest, &AV, &RV );
    GENLC.genmag_rest[epoch] += x[1] ;
    
    // now do 2nd nearest filter
    x[2] = get_snxt8__( &OPT_SNXT, &ifilt_rest2, &Trest, &AV, &RV );
    GENLC.genmag_rest2[epoch] += x[2] ;
    
    // 3rd nearest filter
    x[3] = get_snxt8__( &OPT_SNXT, &ifilt_rest3, &Trest, &AV, &RV );
    GENLC.genmag_rest3[epoch] += x[3] ;
    
    
  } // end of epoch loop


 KCOR:

  z     = GENLC.REDSHIFT_HELIO ;

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {  

    Trest     = GENLC.epoch_rest[epoch]; 

    ifilt_obs   = GENLC.IFILT_OBS[epoch];
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

    ifilt_rest1 = GENLC.IFILTMAP_REST1[ifilt_obs] ;
    ifilt_rest2 = GENLC.IFILTMAP_REST2[ifilt_obs] ;
    ifilt_rest3 = GENLC.IFILTMAP_REST3[ifilt_obs] ;

    lamdif[1]  = GENLC.LAMDIF_REST1[ifilt_obs];
    lamdif[2]  = GENLC.LAMDIF_REST2[ifilt_obs];
    lamdif[3]  = GENLC.LAMDIF_REST3[ifilt_obs];

    ifilt_rest_tmp[1] = ifilt_rest1 ;
    ifilt_rest_tmp[2] = ifilt_rest2 ;
    ifilt_rest_tmp[3] = ifilt_rest3 ;

    mag[1]  = GENLC.genmag_rest[epoch]  ;
    mag[2]  = GENLC.genmag_rest2[epoch] ;    
    mag[3]  = GENLC.genmag_rest3[epoch] ;

    kcor = kcorfun8_ ( &ifilt_obs, &ifilt_rest_tmp[1], 
		       &mag[1], &lamdif[1],  &Trest, &z, &AVwarp[1] );

    if ( isnan( AVwarp[2]) ) {
      sprintf(c1err,"AVwarp=nan for T8=%5.1f  mag8[1,2]=%6.2f,%6.2f",
	      Trest, mag[1], mag[2] );
      sprintf(c2err,"ifilt_rest[1,2]=%d,%d (%c,%c) "
	      ,ifilt_rest1, ifilt_rest2
	      ,FILTERSTRING[ifilt_rest1]
	      ,FILTERSTRING[ifilt_rest2]
	      );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( mag[1] < -10. ) {
      NAVWARP_OVERFLOW[0]++ ;             // grand total
      NAVWARP_OVERFLOW[ifilt_rest1]++ ;   // filter-dependent sum
    }

    // add up all the contributions

    mag_obs = mag[1] + kcor + GENLC.DLMU ; 

    // Jan 22, 2010: 
    // if the spectral warping is crazy, set mag_obs to crazy value
    if ( Trest > 0 ) {
      if ( fabs(mag[1]) > 90.0 ) { mag_obs = MAG_UNDEFINED ; }
      if ( fabs(mag[2]) > 90.0 ) { mag_obs = MAG_UNDEFINED ; }
    }
    else {
      if ( fabs(mag[1]) > 90.0 ) { mag_obs = 99; }
      if ( fabs(mag[2]) > 90.0 ) { mag_obs = 99; }
    }

    if ( fabs(Trest) < -90.01 ) {
      printf(" BOOST: %c(%c) -> %c : Trest=%6.1f (MJD=%7.1f) iep=%d \n"
	     ,FILTERSTRING[ifilt_rest1]
	     ,FILTERSTRING[ifilt_rest2]
	     ,FILTERSTRING[ifilt_obs]
	     ,Trest, GENLC.MJD[epoch], epoch ) ;

      printf("\t M%c = %6.2f(M%c) + %6.3f(kcor) + %7.3f(mu) = %7.3f \n"
	     ,FILTERSTRING[ifilt_obs]
	     ,mag[1]
	     ,FILTERSTRING[ifilt_rest1]
	     ,kcor, GENLC.DLMU,  mag_obs );

    }

    // store Kcor info to include in SNDATA files

    sprintf(GENLC.kcornam[epoch], "K_%c%c"
	    ,FILTERSTRING[ifilt_rest1]
	    ,FILTERSTRING[ifilt_obs] 
	    ) ;
    sprintf(GENLC.warpcolnam[epoch], "%c-%c"
	    ,FILTERSTRING[ifilt_rest1]
	    ,FILTERSTRING[ifilt_rest2] 
	    ) ;

    GENLC.kcorval8[epoch]        = kcor ;
    GENLC.warpcolval8[epoch]     = mag[1] - mag[2] ;
    GENLC.AVwarp8[epoch]         = AVwarp[2] ;
    GENLC.ifilt_AVwarp[epoch][0] = ifilt_rest1 ; // closest filter
    GENLC.ifilt_AVwarp[epoch][1] = ifilt_rest2 ; // 2nd closest

      // store true obs mag in GENLC structure

    GENLC.genmag_obs[epoch] = mag_obs;

    // load model mag-err for observer frame using the model-error
    // from the nearest rest-frame filter (ifilt_rest1).

    GENLC.generr_obs[epoch] = GENLC.generr_rest[epoch] ;

  } // end of epoch loop

  return ;

}  // end of genmag_boost


// ***********************************
void genmag_MWXT_fromKcor(void) {


  /***********
  Created June 11, 2008 by R.Kessler
  [major update: Sep 2013]

  For rest-frame models only, compute and add MilkyWay extinction 
  for all observer-frame mags. Also make correction if user 
  RVMW != RVMW(kcorTable). The correction is an approximation based
  on the MWXT-difference at the central filter wavelengths, 
  and should be OKay for small corrections. However, for large 
  RV differences (order 1 or 2) the approximation may not work so well,
  particularly in the bluer bands.


  May 15, 2009:  check GENLC.DOFILT[ifilt_obs] logical and skip 
                 calc if observer-mag > 30.

  May 10, 2012: allow RV != RV_MWDUST = 3.1.
                Compute difference from INPUTS.RVMW and 3.1 at
                central wavelength of observer-frame filter.
                Works the same for rest-frame and obs-frame models.

  Feb 13, 2013: return immediately for SALT2 model.

  Sep 22, 2013: 
    - return immediatly on obs-frame model.
    - move RV-correction into get_mwxt8 function; 
      see extra RV and OPT_COLORLAW args to get_mwxt8.

  **********/

  double  z, Trest, mwebv, AVwarp, MWXT, RV    ;
  int  epoch, ifilt_obs, OPT ;
  //  char fnam[] = "genmag_MWXT_fromKcor" ;

  // ---------- BEGIN -----------

  if ( GENFRAME_OPT == GENFRAME_OBS ) { return ; }

  z      = GENLC.REDSHIFT_HELIO ;
  mwebv  = GENLC.MWEBV_SMEAR ;       //  smeared MWEBV
  RV     = INPUTS.RV_MWCOLORLAW ;
  OPT    = INPUTS.OPT_MWCOLORLAW ;

  // ------------------------------------------------

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {  

    // don't bother with really large mags that have ~ zero flux
    if ( GENLC.genmag_obs[epoch]  > 30.0 ) { continue ; }

    ifilt_obs  = GENLC.IFILT_OBS[epoch];
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

    Trest   = GENLC.epoch_rest[epoch]; 

    AVwarp  = GENLC.AVwarp8[epoch] ;
    MWXT  = get_mwxt8__(&ifilt_obs, &Trest, &z, &AVwarp, &mwebv, &RV, &OPT);
      
    // increment observer-frame magnitude.
    GENLC.genmag_obs[epoch] += MWXT ;      


  } // end of epoch loop



} // end of genmag_MWXT


// ************************************
int NEPFILT_GENLC(
		int opt   // (I) +1 => retrieve from GENLC,  -1 => load GENLC
		,int ifilt_obs  // (I) 
		) {

  // Returns number of 'ifilt_obs' epochs stored in GENLC & GENFILT
  // Fill GENFILT structure for 'ifilt_obs;
  // loads epochs in order for each filter so that epoch-array
  // can be passed to genmag_xxx functions.

  int  epoch, ifilt_tmp, NEP ;
  double Trest1, Trest2, Tdif;

  char fnam[] = "NEPFILT_GENLC" ;

  // --------------- BEGIN -------------

  NEP = 0;

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {

    ifilt_tmp = GENLC.IFILT_OBS[epoch];

    if ( ifilt_tmp == ifilt_obs ) {
      NEP++ ;

      if ( NEP >= MXEPSIM_PERFILT ) {
	sprintf(c1err,"NEP=%d exceeds bound for ifilt_obs=%d(%c).",  
		NEP, ifilt_obs, FILTERSTRING[ifilt_obs] );
	sprintf(c2err,"%s", "Check MXEPSIM_PERFILT");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      if ( opt == +1 ) {

	GENFILT.Trest[ifilt_obs][NEP]   = GENLC.epoch_rest[epoch];
	GENFILT.Tobs[ifilt_obs][NEP]    = GENLC.epoch_obs[epoch];

	GENFILT.genmag_obs[ifilt_obs][NEP]    = GENLC.genmag_obs[epoch];
	GENFILT.genmag_rest[ifilt_obs][NEP]   = GENLC.genmag_rest[epoch];
	GENFILT.genmag_rest2[ifilt_obs][NEP]  = GENLC.genmag_rest2[epoch];
	GENFILT.genmag_rest3[ifilt_obs][NEP]  = GENLC.genmag_rest3[epoch];

	GENFILT.genmag_smear[ifilt_obs][NEP]  = GENLC.magsmear8[epoch];

	GENFILT.generr_obs[ifilt_obs][NEP]    = GENLC.generr_obs[epoch];
	GENFILT.generr_rest[ifilt_obs][NEP]   = GENLC.generr_rest[epoch];
	GENFILT.generr_rest2[ifilt_obs][NEP]  = GENLC.generr_rest2[epoch];

	//printf(" xxx load Trest(%c=%d) = %f at  epoch=%d,  FILTEPOCH=%d \n", 
	// FILTERSTRING[ifilt_obs], ifilt_obs, GENLC.epoch_rest[epoch], epoch,NEP);
      }
      else if ( opt == -1 ) {

	// first make sure that Trest matches up

	Trest1	= GENFILT.Trest[ifilt_obs][NEP] ;
	Trest2  = GENLC.epoch_rest[epoch];
	Tdif    = fabs(Trest1-Trest2) ;
	if ( Tdif > 0.001 ) {
	  sprintf(c1err,"GENFILT[ifilt_obs=%d][NEP=%d] = %8.3f", 
		  ifilt_obs, NEP, Trest1);
	  sprintf(c2err,"GENLC.epoch_rest[ep=%d]=%8.3f  (Tdif=%8.3f) \n", 
		  epoch, Trest2, Tdif );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	}

	GENLC.genmag_obs[epoch]   = GENFILT.genmag_obs[ifilt_obs][NEP] ;
	GENLC.genmag_rest[epoch]  = GENFILT.genmag_rest[ifilt_obs][NEP] ; 
	GENLC.genmag_rest2[epoch] = GENFILT.genmag_rest2[ifilt_obs][NEP] ; 
	GENLC.genmag_rest3[epoch] = GENFILT.genmag_rest3[ifilt_obs][NEP] ; 

	GENLC.magsmear8[epoch]   = GENFILT.genmag_smear[ifilt_obs][NEP] ;

	GENLC.generr_obs[epoch]   = GENFILT.generr_obs[ifilt_obs][NEP] ;
	GENLC.generr_rest[epoch]  = GENFILT.generr_rest[ifilt_obs][NEP] ;
	GENLC.generr_rest2[epoch] = GENFILT.generr_rest2[ifilt_obs][NEP] ;
      }
    }

  }

  return NEP ;


}  // end of NEPFILT_GENLC


// ************************************
void init_genmodel(void) {

  /**********
  History:


 Feb 13, 2012: 
    make calls to get_LAMRANGE_[model] to fill  GENLC.RESTLAM_MODEL[2].
    Init GENLC.RESTLAM_MODEL to be wide open before calling get_LAMRANGE_XXX.

 Apr 7, 2012: set IFLAG_GENSMEAR

 Jul 20 2018: for LCLIB model, set GENRANGE_PEAKMJD= start of season,
              and GENRANGE_TREST = 0 to MJDLAST

 Nov 23 2020: pass SURVEY arg to init_genmag_SALT2.

  ************/

  char *GENMODEL        = INPUTS.GENMODEL;
  char *GENMODEL_EXTRAP = INPUTS.GENMODEL_EXTRAP_LATETIME ;
  char  covFile[] = ""  ;
  char *ARGLIST_PySEDMODEL ;
  int istat, OPTMASK,  ifilt, ifilt_obs, ifilt_rest ;
  float scale_covar_flt ;
  char fnam[] = "init_genmodel" ;

  //--------- BEGIN --------

  // init a few things
  for ( ifilt_obs=0; ifilt_obs < MXFILTINDX; ifilt_obs++ ) {
    NSKIP_FILTER[ifilt_obs] = 0;
    ZVALID_FILTER[0][ifilt_obs] =  +999. ;
    ZVALID_FILTER[1][ifilt_obs] =  -999. ;
  }

  GENLC.RESTLAM_MODEL[0] =  1000. ;
  GENLC.RESTLAM_MODEL[1] = 30000.0 ;

  // default is to generate intrinsic smearing at center of filter only
  IFLAG_GENSMEAR = IFLAG_GENSMEAR_FILT ;

  LGEN_SNIA = 0 ;

  // =========================

  if ( INDEX_GENMODEL == MODEL_STRETCH ) {

    init_genmag_stretch2(GENMODEL, GENLC.FILTLIST_REST );
    LGEN_SNIA = 1 ;
  }

  else if ( INDEX_GENMODEL == MODEL_FIXMAG ) {
    GENLC.SIMTYPE  = MODEL_FIXMAG ;
  }
  else if ( INDEX_GENMODEL == MODEL_SIMLIB ) {
    printf("\n Read SIM_MAGOBS from SIMLIB file. \n");
  }
  else if ( INDEX_GENMODEL == MODEL_MLCS2k2 ) {

    scale_covar_flt = 1.0;  // 5/02/2009: scale error at mag-smear stage
    istat = init_genmag_mlcs2k2(GENMODEL, covFile, scale_covar_flt, 
				GENLC.FILTLIST_REST );
    get_LAMRANGE_mlcs2k2(&GENLC.RESTLAM_MODEL[0], &GENLC.RESTLAM_MODEL[1] );
    LGEN_SNIA = 1 ;
  }

  else if ( INDEX_GENMODEL == MODEL_SNOOPY ) {

    OPTMASK = 0;
    init_genmag_snoopy(GENMODEL, OPTMASK, GENLC.FILTLIST_REST );

    get_LAMRANGE_snoopy(&GENLC.RESTLAM_MODEL[0], &GENLC.RESTLAM_MODEL[1] );
    LGEN_SNIA = 1 ;
  }

  else if ( INDEX_GENMODEL == MODEL_S11DM15 ) {

    // init generic part of any SEDMODEL (filter & primary ref)
    init_genSEDMODEL();

    OPTMASK = 0 ;
    istat = init_genmag_S11DM15(GENMODEL,OPTMASK);
    LGEN_SNIA = 1 ;
  }

  else if ( INDEX_GENMODEL == MODEL_SALT2 ) {

    // init generic part of any SEDMODEL (filter & primary ref)
    init_genSEDMODEL();

    // model-specific init
    // xxxx    OPTMASK = 0;
    OPTMASK = INPUTS.GENMODEL_MSKOPT; 

    /* NOT_YET
    if ( INPUTS.REQUIRE_DOCANA  ) { OPTMASK |= OPENMASK_REQUIRE_DOCANA; }
    */

    istat = init_genmag_SALT2(GENMODEL, GENMODEL_EXTRAP, OPTMASK) ;

    get_LAMRANGE_SEDMODEL(1,&GENLC.RESTLAM_MODEL[0],&GENLC.RESTLAM_MODEL[1] );

    if ( istat != 0 ) {
      sprintf(c1err,"init_genmag_SALT2 failed.");
      errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
    }

    // set flag to generate intrinsic smear vs. wavelength
    // (inside genmag_SALT2).
    IFLAG_GENSMEAR = IFLAG_GENSMEAR_LAM ;
    LGEN_SNIA = 1 ;
  }

  else if ( INDEX_GENMODEL == MODEL_SIMSED ) {

    // init generic part of any SEDMODEL (filter & primary ref)
    init_genSEDMODEL();

    // model-specific init

    OPTMASK = INPUTS.GENMODEL_MSKOPT ;
    if( INPUTS.USE_BINARY_SIMSED > 0 ) 
      { OPTMASK += OPTMASK_INIT_SIMSED_BINARY; }
    if( INPUTS.JOBID             > 0 ) 
      { OPTMASK += OPTMASK_INIT_SIMSED_BATCH; }

    istat = init_genmag_SIMSED (INPUTS.GENMODEL
			       ,INPUTS.PATH_BINARY_SIMSED
			       ,GENLC.SURVEY_NAME
			       ,INPUTS.KCOR_FILE
			       ,OPTMASK );

    get_LAMRANGE_SEDMODEL(1,&GENLC.RESTLAM_MODEL[0], &GENLC.RESTLAM_MODEL[1] );
    
    if ( istat == 0 ) {
      sprintf(c1err,"init_genmag_SIMSED failed.");
      errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
    }

    // check that SIMSED parameters in sim-input file match
    // those in the SIMSED model.
    checkpar_SIMSED();
    LGEN_SNIA = 0 ;  // July 2017
  }

  else if ( IS_PySEDMODEL ) {

    OPTMASK  = INPUTS.GENMODEL_MSKOPT;

    ARGLIST_PySEDMODEL = (char*)malloc(400*sizeof(char) );
    sprintf(ARGLIST_PySEDMODEL,"RANSEED %d  %s",
	    INPUTS.ISEED, INPUTS.GENMODEL_ARGLIST );

    int NPAR; char NAMES_HOSTPAR[200];  double VAL_HOSTPAR=0.0 ;
    NPAR = fetch_HOSTPAR_GENMODEL(1, NAMES_HOSTPAR, &VAL_HOSTPAR);

    // init generic part of any SEDMODEL (filter & primary ref)
    init_genSEDMODEL();

    init_genmag_PySEDMODEL(INPUTS.GENMODEL, INPUTS.MODELPATH, 
			   OPTMASK, ARGLIST_PySEDMODEL, NAMES_HOSTPAR);
  }

  else if ( INDEX_GENMODEL == MODEL_NON1ASED ) {

    INPUTS.NON1ASED.NGENTOT   = INPUTS.NGEN ;
    INPUTS.NON1ASED.CIDOFF    = INPUTS.CIDOFF ;
    INPUTS.NON1ASED.IFLAG_GEN = GENLC.IFLAG_GENSOURCE ;
    GENLC.NON1ASED.IFLAG_GEN  = GENLC.IFLAG_GENSOURCE ;
    
    GENLC.NON1ASED.FRAC_PEC1A = INPUTS.RATEPAR_PEC1A.SEASON_FRAC ;
    prep_NON1ASED( &INPUTS.NON1ASED, &GENLC.NON1ASED );

    init_genSEDMODEL();       // pass filters and primary ref
    init_genmag_NON1ASED(-9,&INPUTS.NON1ASED); // do one-time inits for SEDs
  }
  else if ( INDEX_GENMODEL == MODEL_NON1AGRID ) {
    double FRAC_PEC1A = INPUTS.RATEPAR_PEC1A.SEASON_FRAC ;
    init_genmag_NON1AGRID(INPUTS.NON1AGRID_FILE, FRAC_PEC1A );
  }
  else if ( INDEX_GENMODEL == MODEL_LCLIB ) {
    
    double TOBS_RANGE_MAX = INPUTS.GENRANGE_MJD[1] - INPUTS.GENRANGE_MJD[0];
    int    IFLAG_RANSTART = (INPUTS.JOBID>0) ;
    OPTMASK = INPUTS.GENMODEL_MSKOPT ;
    init_genmag_LCLIB(INPUTS.LCLIB_FILE, INPUTS.LCLIB_TEMPLATE_EPOCHS,
		      GENLC.SURVEY_NAME, INPUTS.GENFILTERS,
		      TOBS_RANGE_MAX, INPUTS.NGEN, IFLAG_RANSTART, OPTMASK );
    
    //  Aug 22 2017:
    //  Hack ugly spaghetti code to run INIT_HOSTLIB if LCLIB has a 
    //  REDSHIFT parameter (e.g., AGN). Note that redshift range is read 
    //  from  the LCLIB header, not from sim-input file.
    //  SHOULD AVOID USING LCLIB MODEL for EXTRA-GALACTIC models.
    if ( LCLIB_INFO.IPAR_REDSHIFT >= 0 && LCLIB_INFO.HOSTLIB_MSKOPT>0 )  { 
      INPUTS.GENRANGE_REDSHIFT[0] = LCLIB_INFO.REDSHIFT_RANGE[0] ;
      INPUTS.GENRANGE_REDSHIFT[1] = LCLIB_INFO.REDSHIFT_RANGE[1] ;
      INPUTS.HOSTLIB_MSKOPT       = LCLIB_INFO.HOSTLIB_MSKOPT ;
      INPUTS.HOSTLIB_USE          = 1 ;
      INIT_HOSTLIB(); 
    }
    else {
      sprintf(INPUTS.HOSTLIB_FILE,"NONE");              // no host  
      sprintf(INPUTS_SEARCHEFF.USER_zHOST_FILE,"ZERO"); // force EFF_zHOST=0
    }

  }     // end LCLIB if block

  else {
    sprintf(c1err,"%s is not a valid genmag-model", INPUTS.MODELNAME);
    sprintf(c2err,"Check GENMODEL keyword in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // convert char-filter list into integer map and inverse map

  if ( GENFRAME_OPT == GENFRAME_REST ) {

    printf("\n Rest-frame MODEL FILTERS : '%s' \n", GENLC.FILTLIST_REST );

    GENLC.NFILTDEF_REST = 
      PARSE_FILTLIST(GENLC.FILTLIST_REST, GENLC.IFILTMAP_REST) ;

    // create inverse map need for genmag_mlcs2k2().
    for ( ifilt=0 ; ifilt < GENLC.NFILTDEF_REST ; ifilt++ ) {
      ifilt_rest = GENLC.IFILTMAP_REST[ifilt] ;
      GENLC.IFILTINVMAP_REST[ifilt_rest] = ifilt;
    }

    if ( INDEX_GENMODEL == MODEL_MLCS2k2 ) { init_covar_mlcs2k2(); }

  }

  if ( LGEN_SNIA ) 
    { GENLC.SIMTYPE  = INPUTS.SNTYPE_Ia_SPEC ; }

  if ( INPUTS.GENTYPE_SPEC > 0 ) 
    { GENLC.SIMTYPE  = INPUTS.GENTYPE_SPEC ; }
 
  prep_dustFlags();


  return ;

}  // end of init_genmodel


// *********************************
void init_genSEDMODEL(void) {

  // Created Oct 15, 2009 by R.Kessler
  // Perform generic part of initialization for observer-frame
  // SED model.
  //
  // Read K-cor files to fetch primary SED and filter-responses.
  // Note  that the K-corrections tables are NOT used.
  //
  // Dec 19, 2011: abort if zmin <= 1.E-9
  // Dec 14, 2013: include safety margin on zmin,zmax to allow for
  //               zcmb -> zhelio because we generate zcmb, but
  //               generate with zhelio.
  //
  // Jan 23, 2014: call extraFilters_4genSmear( ...) to determine
  //               extra filters needed by genSmear functions.
  //               See NFILT_extra and IFILTOBS_extra.
  //
  // --------------------
  char  fnam[] = "init_genSEDMODEL" ;
  char   filtName[40], surveyName[80], cfilt[2]    ;
  //  char  *survey = GENLC.SURVEY_NAME ;

  int 
    ifilt, ifilt_obs
    ,NLAM, NFILT, NZBIN, NFILT_extra, LEN=40
    ,DONEFILT[MXFILTINDX]
    ,IFILTLOCAL[MXFILTINDX]
    ,IFILTOBS_extra[MXFILTINDX]
    ,MASKFRAME[MXFILTINDX]      
    ;
  
  double  lamshift, magprim, zmin, zmax, dz, z1min, z1max,logzdif  ;

  // ------------ BEGIN -------------

  reset_SEDMODEL();

  // allocate temporary arrays to hold primary SED and filter trans.
  int MEMD = MXLAMSIM * sizeof(double);
  genSEDMODEL.lam         = (double*) malloc( MEMD );
  genSEDMODEL.primaryFlux = (double*) malloc( MEMD );
  genSEDMODEL.TransSN     = (double*) malloc( MEMD );
  genSEDMODEL.TransREF    = (double*) malloc( MEMD );

  // -----------
  sprintf(GENLC.primary,"%s", "NULL");
  get_primary__(GENLC.primary, &NLAM, 
		genSEDMODEL.lam, genSEDMODEL.primaryFlux, LEN) ;
  //		strlen(GENLC.primary) );

  init_primary_SEDMODEL( GENLC.primary, NLAM, 
			 genSEDMODEL.lam, genSEDMODEL.primaryFlux );


  init_MWXT_SEDMODEL(INPUTS.OPT_MWCOLORLAW, INPUTS.RV_MWCOLORLAW);

  // init z-range; include safety margin for zcmb->zhelio 
  // and for VPEC
  z1min    = 1. + INPUTS.GENRANGE_REDSHIFT[0] ;
  z1max    = 1. + INPUTS.GENRANGE_REDSHIFT[1] ;
  dz       = INPUTS.VEL_CMBAPEX/LIGHT_km  ;
  dz      += INPUTS.GENSIGMA_VPEC*5.0/LIGHT_km ;  // 5sigma safety margin
  zmin     = INPUTS.GENRANGE_REDSHIFT[0] - dz*z1min ;
  zmax     = INPUTS.GENRANGE_REDSHIFT[1] + dz*z1max ;

  if ( zmin <= 1.0E-9 ) {
    sprintf(c1err, "Invalid zmin=%le (must be > 0)", zmin);
    sprintf(c2err, "Check GENRANGE_REDSHIFT key in sim-input file." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  logzdif =  log10(zmax) - log10(zmin) ;
  NZBIN = (int)( 50. * logzdif ) + 1 ;
  init_redshift_SEDMODEL(NZBIN, zmin, zmax);

  NFILT = INPUTS.NFILTDEF_OBS ;
  for ( ifilt = 0; ifilt < MXFILTINDX; ifilt++ ) {
    DONEFILT[ifilt]   = 0 ;
    IFILTLOCAL[ifilt] = 0 ;
    MASKFRAME[ifilt]  = GENFRAME_OBS ;
    if ( ifilt < NFILT ) 
      { IFILTLOCAL[ifilt] = INPUTS.IFILTMAP_OBS[ifilt]; }
  }


  // check for extra filters needed by the intrinsic scatter model.
  extraFilters_4genSmear(INPUTS.GENMAG_SMEAR_MODELNAME,
			 &NFILT_extra, IFILTOBS_extra) ;

  printf("\n");

  for(ifilt=0; ifilt < NFILT_extra ; ifilt++ ) {
    ifilt_obs           = IFILTOBS_extra[ifilt] ;
    IFILTLOCAL[NFILT]   = ifilt_obs ;
    MASKFRAME[NFILT]    = GENFRAME_OBS + GENFRAME_REST ; // check both
    NFILT++ ;
    sprintf(cfilt,  "%c", FILTERSTRING[ifilt_obs] ); 
    printf("  %s : add %s filter needed by %s \n", 
	   fnam, cfilt, INPUTS.GENMAG_SMEAR_MODELNAME ); fflush(stdout);
  }


  // pass filter info to SEDMODEL generator

  for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
      
    ifilt_obs = IFILTLOCAL[ifilt] ;	
        
    // avoid initializing same filter twice.
    if ( DONEFILT[ifilt_obs] > 0 ) { continue ; }

    // get_filttrans returns separate trans for SN and REF,
    // but these are always the same for sims.

    sprintf(filtName,"NULL") ;
    sprintf(surveyName,"NULL") ;

    get_filttrans__(&MASKFRAME[ifilt],     // (I) obs/rest-frame mask
		    &ifilt_obs,            // (I) absolute filter index
		    surveyName,            // (I) name(s) of survey (11.2020)
		    filtName,              // (O) full name of filter
		    &magprim,              // (O) mag of primary ref 
		    &NLAM,                 // (O) Number of lambda bins
		    genSEDMODEL.lam,       // (O) lambda array 
		    genSEDMODEL.TransSN,   // (O) filter trans 
		    genSEDMODEL.TransREF,  // (O) idem
		    80,40 ); 

    if ( NLAM > MXLAMSIM ) {
      sprintf(cfilt,  "%c", FILTERSTRING[ifilt_obs] );
      sprintf(c1err,"NLAM(%s) = %d exceeds bound of MXLAMSIM=%d",
	      cfilt, NLAM, MXLAMSIM );
      sprintf(c2err,"Check filter trans files.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    lamshift = 0.0 ;
    init_filter_SEDMODEL(ifilt_obs, filtName, surveyName, magprim, NLAM, 
			 genSEDMODEL.lam, 
			 genSEDMODEL.TransSN, 
			 genSEDMODEL.TransREF, 
			 lamshift );

    DONEFILT[ifilt_obs] = 1 ;
  }


  // free temp memory
  free(genSEDMODEL.lam);
  free(genSEDMODEL.primaryFlux);
  free(genSEDMODEL.TransSN);
  free(genSEDMODEL.TransREF);

} // end of init_genSEDMODEL


// *******************************
void check_model_default(int index_model ) {

  // if default model is used, then change INPUTS.MODELNAME
  // to that listed in [model].default file.

  int   j, lentmp, ICMP;
  FILE *fp;
  char *cptr, *PTR_GENMODEL ;
  char  model[40],suffix[20], defaultFile[MXPATHLEN], defaultPath[MXPATHLEN] ;
  char  model_lc[40] = { "  " } ;    
  char  fnam[] = "check_model_default"  ;


  // ----------------- BEGIN ---------------

  if ( index_model == MODEL_NON1ASED  ) return ;
  if ( index_model == MODEL_NON1AGRID ) return ;
  if ( index_model == MODEL_FIXMAG    ) return ;
  if ( index_model == MODEL_SIMLIB    ) return ;
  if ( index_model == MODEL_LCLIB     ) return ;
  if ( IS_PySEDMODEL                  ) return ;

  sprintf(model, "%s", INPUTS.MODELNAME );

  PTR_GENMODEL = GENMODEL_NAME[index_model][0] ; // generic model name

  // if there is no dot (.) in the model name,
  // then assume that we have the generic name and
  // read the [model].default file to get the model version.

  cptr = strstr ( model, "." ) ; 
  if ( cptr == NULL ) {

    // check if model is lower/upper case to determine if
    // .default  suffix is lower/upper case.

    lentmp = strlen(PTR_GENMODEL);
    for (j=0; j<lentmp; j++) { model_lc[j] = tolower(PTR_GENMODEL[j]); }
    model_lc[lentmp] = '\0' ;
    ICMP = strcmp(model_lc,PTR_GENMODEL) ;

    if ( ICMP == 0 ) 
      { sprintf(suffix, "%s", "default"); }
    else
      { sprintf(suffix, "%s", "DEFAULT"); }


    // construct full name of file with default model version
    sprintf(defaultPath,"%s/models/%s", 
	    PATH_SNDATA_ROOT, PTR_GENMODEL );
    sprintf(defaultFile,"%s/%s.%s", 
	    defaultPath, PTR_GENMODEL, suffix);

    if ( (fp = fopen(defaultFile, "rt")) == NULL ) {
      sprintf(c1err,"Default version file does not exist:");
      sprintf(c2err,"%s",  defaultFile );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // read default model-version from default-file
    readchar(fp, model ) ;

    fclose(fp);

    // re-load global GENMODEL with default model version
    sprintf(INPUTS.GENMODEL, "%s/%s", defaultPath, model);
    sprintf(INPUTS.MODELNAME, "%s", model ); // copy full name 
  }


} // end of check_model_default

// ******************************
void init_kcor_legacy(char *kcorFile) {

  /********
   Created May 18, 2008 by R.Kessler

   Read K-corrections and filter transmissions (rdkcor).
   We do not yet know GENFRAME_OPT ('rest' or 'obs'),
   but always read the K-cors to be safe, since even
   obs-frame models will likely need the filter transmissions.

   Also do miscelaneous chores like fetching
   AB (zeropoint) offets.

   Sep 1, 2010: fix evil bug: tmpoff[20] -> tmpoff[MXFILTINDX]

   May 13, 2013: few fixes for MODEL_FIXMAG.

   Sep 22 2013: for rest-frame models call get_kcor_mwpar(RV,OPT_COLORLAW)

   Jun 20 2017: pass FILTER_LAMSHIFT argument to set_survey_

   Nov 2 2017: add AB offsets to user offsets; see tmpoff_kcor
   Apr 24 2019: for FIXMAG model, return before doing rest-frame stuff.

   Jul 26 2019: 
     + get_kcor_mwpar() -> get_kcor_info(), and NKCOR is returned. 
     + for rest-frame model, abort if NKCOR==0.

   Sep 28 2020: 
     + use new function find_pathfile(...) to find kcor_file by searching
       current, SNDATA_ROOT and PATH_USER_INPUT.

   Aug 18 2021: enable for FIXMAG model

  *********/

  int ISMODEL_FIXMAG = ( INDEX_GENMODEL == MODEL_FIXMAG );
  int ISMODEL_SIMLIB = ( INDEX_GENMODEL == MODEL_SIMLIB );
  int ierrstat, ifilt, ifilt_obs, OPT, NKCOR=0 ;
  float tmpoff_kcor[MXFILTINDX] ;
  char   copt[40], xtDir[MXPATHLEN], cfilt[4];
  char fnam[] = "init_kcor_legacy" ;

  // -------------- BEGIN --------------

  // set fortran arrays needed for K-correction lookup
  // The GENLC variables are set in SIMLIB_open
 
  set_survey__( GENLC.SURVEY_NAME  
		,&GENLC.NFILTDEF_OBS
		,GENLC.IFILTMAP_OBS
		,INPUTS.TMPOFF_LAMFILT
		,strlen(GENLC.SURVEY_NAME)   
		);

  rdkcor_(kcorFile, &ierrstat, strlen(kcorFile) );

  if ( ierrstat != 0 ) {
    sprintf(c1err, "Could not open kcor file: '%s'", kcorFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
  } 


  // get zeropoint (AB) offsets from kcor/calib file.
  sprintf(copt,"ABOFF") ;
  get_filtmap__(copt, tmpoff_kcor, strlen(copt) );


  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs  = GENLC.IFILTMAP_OBS[ifilt] ;

    if ( !ISMODEL_FIXMAG ) 
      { INPUTS.GENMAG_OFF_ZP[ifilt_obs]  += tmpoff_kcor[ifilt];  }

    // get mean and RMS of filter-wavelength 
    OPT = GENFRAME_OBS ;
    get_filtlam__(&OPT, &ifilt_obs
		  ,&INPUTS.LAMAVG_OBS[ifilt_obs]
		  ,&INPUTS.LAMRMS_OBS[ifilt_obs]   
		  ,&INPUTS.LAMMIN_OBS[ifilt_obs]  
		  ,&INPUTS.LAMMAX_OBS[ifilt_obs]   
		  );
  }  // end ifilt loop


  // misc. inits
  for ( ifilt=0; ifilt<MXFILTINDX; ifilt++ ) { 
    NAVWARP_OVERFLOW[ifilt] = 0; 
    GENLC.IFLAG_SYNFILT_SPECTROGRAPH[ifilt] = 0 ;
  }

  // xxx mark delete  if ( ISMODEL_FIXMAG ) { return ; } // 4.24.2019
  if ( ISMODEL_SIMLIB ) { return ; } // 11.22.2019

  // init MW extinction (moved from end of init_genmodel on Aug 30 2010)
  if ( GENFRAME_OPT == GENFRAME_REST )  {

    // fetch kcor-info:
    // Number of K-cor rables, and params used to compute MW extinct in 
    // kcor files: needed to re-compute MWXT if RV is changed.
    get_kcor_info__(&NKCOR, &GENLC.kcor_RVMW, &GENLC.kcor_OPT_MWCOLORLAW );

    double RVDIF ;
    RVDIF = fabs ( INPUTS.RV_MWCOLORLAW - GENLC.kcor_RVMW) ;
    if ( RVDIF > 0.001 ) {
      printf("\n");
      printf("   User RVMW=%.3f != RVMW(kcorTable)=%.3f  \n",
	     INPUTS.RV_MWCOLORLAW, GENLC.kcor_RVMW) ;
      printf("   -> Compute MWXT correction with OPT_MWCOLORLAW=%d\n", 
	     GENLC.kcor_OPT_MWCOLORLAW );
      fflush(stdout);
    }

    if ( NKCOR == 0 ) {
      sprintf(c1err,"zero k-correction tables found. ");
      sprintf(c2err,"Rest-frame model requires valid k-cor tables.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( INPUTS.OPT_SNXT  == OPT_SNXT_CCM89 ) {
      init_xthost__(&INPUTS.OPT_SNXT);
    }
    else if ( INPUTS.OPT_SNXT  == OPT_SNXT_SJPAR ) {
      print_banner("Init HOST-GALAXY Extinction Parameters");
      sprintf(xtDir, "%s ",INPUTS.MODELPATH); // leave blank space for fortran
      rdxtpar_(xtDir, strlen(xtDir) );      
    }



  } // end rest-frame


  // check for optional spectrograph information (July 2016)
  sprintf(copt,"SPECTROGRAPH") ;
  int NF = get_filtmap__(copt, tmpoff_kcor, strlen(copt) );
  GENLC.NFILTDEF_SPECTROGRAPH = NF ;
  GENLC.FILTERLIST_SPECTROGRAPH[0] = 0 ;
  for(ifilt=0; ifilt < NF; ifilt++ ) {
    ifilt_obs = (int)tmpoff_kcor[ifilt] ;
    GENLC.IFILTDEF_SPECTROGRAPH[ifilt]     = ifilt_obs ;
    GENLC.IFILTINV_SPECTROGRAPH[ifilt_obs] = ifilt ;
    GENLC.IFLAG_SYNFILT_SPECTROGRAPH[ifilt_obs] = 1 ;
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
    strcat(GENLC.FILTERLIST_SPECTROGRAPH,cfilt);
  }

  // check for optional SPECTROGRPH info
  read_spectrograph_fits(kcorFile) ;

  if ( SPECTROGRAPH_USEFLAG ) {

    dump_INPUTS_SPECTRO(4,"");

    extend_spectrograph_lambins();

    printf("   Found %d synthetic spectrograph filters (%s) \n",
	   GENLC.NFILTDEF_SPECTROGRAPH, GENLC.FILTERLIST_SPECTROGRAPH );
    fflush(stdout);
  }

  return ;

} // end of init_kcor_legacy

// ********************************** 
void init_kcor_refactor(void) {

  // Oct 2019
  // Begin translating fortran kcor-read codes into C.
  // Call functions in sntools_kcor.c[h]

  int ISMODEL_FIXMAG    = ( INDEX_GENMODEL == MODEL_FIXMAG );
  int    ifilt, ifilt_obs, NBIN, NFILTDEF ;
  double MAGOBS_SHIFT[MXFILTINDX], MAGREST_SHIFT[MXFILTINDX];
  double ZPOFF;
  char fnam[] = "init_kcor_refactor" ;

  // ------------ BEGIN -------------

  KCOR_INFO.NCALL_READ = 0 ;
  for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) 
    { MAGOBS_SHIFT[ifilt] = MAGREST_SHIFT[ifilt] = 0.0 ; }

  READ_KCOR_DRIVER(INPUTS.KCOR_FILE, SIMLIB_GLOBAL_HEADER.FILTERS,
		   MAGREST_SHIFT, MAGOBS_SHIFT );


  NFILTDEF = KCOR_INFO.FILTERMAP_OBS.NFILTDEF;

  for(ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
    ifilt_obs = KCOR_INFO.FILTERMAP_OBS.IFILTDEF[ifilt];

    if ( !ISMODEL_FIXMAG ) {
      // ZPOFF.DAT means filter system isn't right; here we simulate
      // incorrect mag system so that data and sim are corrected the
      // same way with offsets in ZPOFF.DAT file (which are stored
      // in kcor file(
      ZPOFF     = KCOR_INFO.FILTERMAP_OBS.PRIMARY_ZPOFF_FILE[ifilt];
      INPUTS.GENMAG_OFF_ZP[ifilt_obs]  += ZPOFF ;
    }
    
    NBIN = KCOR_INFO.FILTERMAP_OBS.NBIN_LAM[ifilt];
    INPUTS.LAMAVG_OBS[ifilt_obs] = KCOR_INFO.FILTERMAP_OBS.LAMMEAN[ifilt];
    INPUTS.LAMRMS_OBS[ifilt_obs] = KCOR_INFO.FILTERMAP_OBS.LAMRMS[ifilt];
    INPUTS.LAMMIN_OBS[ifilt_obs] = 
      (double)KCOR_INFO.FILTERMAP_OBS.LAM[ifilt][0];
    INPUTS.LAMMAX_OBS[ifilt_obs] = 
      (double)KCOR_INFO.FILTERMAP_OBS.LAM[ifilt][NBIN-1];

  }

  debugexit(fnam);

  return ;

} // end init_kcor_refactor

// *********************************************
int gen_TRIGGER_PEAKMAG_SPEC(void) {

  // Nov 2017
  // Utility to speed generation when the spec-efficieny map
  // depends ONLY on peakMag or peak-color. Idea here is to
  // call GENMAG_DRIVER only for peakmag (not entire LC)
  // and check if spec-trigger is satisfied. This will skip
  // generating mags for entire light curve, and also skip
  // checking for detection at each epoch.
  //
  // Return 1 if trigger is satisfied (nominal epochs restored)
  // Return 0 if trigger fails.       (nominal epochs NOT restored)
  //
  // Jun 2 2018: bail out for LCLIB model because PEAKMAG is ill-defined.

  int  LFIND_SPEC=0, NEP_PEAKONLY=0, iep  ;
  int  MEMI  = (GENLC.NEPOCH+1) * sizeof(int);
  int  MEMD  = (GENLC.NEPOCH+1) * sizeof(double);
  int  DOSPEC = (INPUTS.APPLY_SEARCHEFF_OPT & 2) ;
  double EFF ;

  struct {
    int NEPOCH, *IFILT_OBS, *ISPEAK;
    double *MJD, *TOBS, *TREST ;
  } GENLC_ORIG ;


  //  char fnam[] = "gen_TRIGGER_PEAKMAG_SPEC" ;

  // -------------- BEGIN ----------------

  if ( DOSPEC == 0 ) 
    { return(1); }

  if ( SEARCHEFF_SPEC_INFO.FLAG_PEAKMAG_ONLY == 0 ) 
    { return(1); }

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) 
    { return(1); }

  if ( INDEX_GENMODEL == MODEL_LCLIB ) // Jun 2 2018
    { return(1); }

  GENLC_ORIG.NEPOCH      = GENLC.NEPOCH;
  GENLC_ORIG.IFILT_OBS   = (int*)malloc ( MEMI ) ; 
  GENLC_ORIG.ISPEAK      = (int*)malloc ( MEMI ) ;
  GENLC_ORIG.MJD         = (double*)malloc( MEMD ) ;
  GENLC_ORIG.TOBS        = (double*)malloc( MEMD ) ;
  GENLC_ORIG.TREST       = (double*)malloc( MEMD ) ;

  for(iep=1; iep <= GENLC_ORIG.NEPOCH ; iep++ ) {
    GENLC_ORIG.ISPEAK[iep]    = GENLC.OBSFLAG_PEAK[iep] ;
    GENLC_ORIG.IFILT_OBS[iep] = GENLC.IFILT_OBS[iep] ;
    GENLC_ORIG.MJD[iep]       = GENLC.MJD[iep];
    GENLC_ORIG.TOBS[iep]      = GENLC.epoch_obs[iep];  
    GENLC_ORIG.TREST[iep]     = GENLC.epoch_rest[iep] ;

    if ( GENLC_ORIG.ISPEAK[iep] == 0 ) { continue ; }
    NEP_PEAKONLY++ ;
    GENLC.OBSFLAG_PEAK[NEP_PEAKONLY]       = GENLC_ORIG.ISPEAK[iep] ;
    GENLC.IFILT_OBS[NEP_PEAKONLY]    = GENLC_ORIG.IFILT_OBS[iep] ;
    GENLC.MJD[NEP_PEAKONLY]          = GENLC_ORIG.MJD[iep] ;
    GENLC.epoch_obs[NEP_PEAKONLY]   = GENLC_ORIG.TOBS[iep] ;
    GENLC.epoch_rest[NEP_PEAKONLY]  = GENLC_ORIG.TREST[iep] ;
  }


  GENLC.NEPOCH = NEP_PEAKONLY;
  GENLC.PEAKMAG_TRIGGER_FLAG = 1 ; // global flag
  GENMAG_DRIVER(); 
  LOAD_SEARCHEFF_DATA();
  LFIND_SPEC = gen_SEARCHEFF_SPEC(GENLC.CID, &EFF) ;  // return EFF 
    
  GENLC.NEPOCH = GENLC_ORIG.NEPOCH ; // always needed for init_GENLC

  // check to restore ALL epochs
  if ( LFIND_SPEC ) {
    for(iep=1; iep <= GENLC.NEPOCH ; iep++ ) {
      GENLC.OBSFLAG_PEAK[iep]  = GENLC_ORIG.ISPEAK[iep];
      GENLC.IFILT_OBS[iep]     = GENLC_ORIG.IFILT_OBS[iep] ;
      GENLC.MJD[iep]           = GENLC_ORIG.MJD[iep];
      GENLC.epoch_obs[iep]     = GENLC_ORIG.TOBS[iep]  ;
      GENLC.epoch_rest[iep]    = GENLC_ORIG.TREST[iep] ;
    }
  }


  free(GENLC_ORIG.IFILT_OBS);
  free(GENLC_ORIG.ISPEAK);
  free(GENLC_ORIG.MJD);
  free(GENLC_ORIG.TOBS);
  free(GENLC_ORIG.TREST);

  GENLC.PEAKMAG_TRIGGER_FLAG = 0 ;
  return(LFIND_SPEC) ;

} // end gen_TRIGGER_PEAKMAG_SPEC


// ==========================================
int gen_TRIGGER_zHOST(void) {

  // Created Dec 1 2017
  // Evaluate zHOST trigger and return 1=PASS, 0=FAIL.
  // Intended to reject events before digitization to
  // speed simulation.
  //
  // Mar 15 2018: bail if zHOST is NOT required by user
  // Jun 19 2018: bail of CUTWIN_SNRMAX is specified.

  double EFF, IPASS;
  int  NEPOCH_ORIG   = GENLC.NEPOCH;
  int  APPLY_OPT     = INPUTS.APPLY_SEARCHEFF_OPT ;
  int  REQUIRE_SPEC  = (APPLY_OPT & APPLYMASK_SEARCHEFF_SPEC) ;
  int  REQUIRE_zHOST = (APPLY_OPT & APPLYMASK_SEARCHEFF_zHOST);
  //  char fnam[] = "gen_TRIGGER_zHOST" ;

  // -------------- BEGIN ------------------

  // bail if spec-confirmation is required
  if ( REQUIRE_SPEC ) { return(1); }

  // bail if zHOST is not required  (Mar 15, 2018)
  if ( !REQUIRE_zHOST ) { return(1); }

  // bail if zHOST map is NOT defined.
  if ( INPUTS_SEARCHEFF.NMAP_zHOST == 0 )     { return(1); }

  // bail if CUTWIN_SNRMAX is specified (June 18 2018)
  if ( INPUTS_SEARCHEFF.CUTWIN_SNRMAX_zHOST[0] > 0.0 ) { return(1); }

  // bail if generating GRID
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  )     { return(1); }

  // load data to load randoms used for efficiencies.
  GENLC.NEPOCH=1;  LOAD_SEARCHEFF_DATA();

  // evaluate zHOST trigger
  IPASS = gen_SEARCHEFF_zHOST(GENLC.CID,&EFF);

  // restore GENLC.NEPOCH
  GENLC.NEPOCH = NEPOCH_ORIG; 

  return(IPASS) ;

} // end gen_TRIGGER_zHOST

// *********************************************
void GENMAG_DRIVER(void) {

  // Created July 2016
  // [for spectra refactor, move code from main to here]
  // Driver routine to generate true mags at each epoch & passband.
  //
  // Aut 17 2017: call get_lightCurveWidth

  int ifilt, ifilt_obs, DOFILT ;
  //  char fnam[] = "GENMAG_DRIVER" ;

  // -------------- BEGIN ---------------

  genran_modelSmear(); // randoms for intrinsic scatter

  // this loop is to generate ideal mag in each filter.
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;

    DOFILT = GENLC.DOFILT[ifilt_obs] ;

    if ( DOFILT == 0 ) { continue ; }
    
    genmodel(ifilt_obs,1); 

    if ( GENFRAME_OPT == GENFRAME_REST ) {
      genmodel(ifilt_obs,2);      // 2nd nearest filter
      genmodel(ifilt_obs,3);      // 3rd nearest filter
    } 
  } // ifilt
 

  // spaghetti hack to pass LCLIB redshift and compute HOSTLIB photo-z
  gen_redshift_LCLIB();

  // if generation is in rest frame, apply generic "boost" 
  // to observer frame
  if ( GENFRAME_OPT == GENFRAME_REST ) {
    genmag_boost();           // warp rest-frame SED and apply KCOR      
    genmag_MWXT_fromKcor() ;  //  apply MilkyWay extinction
  }

  // Mar 2016: apply filter-dependent MAGOFF and load peakMag per band
  genmag_offsets();
 
  
  // estimate light-curve Width in each band (Aug 17 2017)
  compute_lightCurveWidths();

  return ;

} // end GENMAG_DRIVER



// *********************************************
void GENFLUX_DRIVER(void) {

  // Created Dec 2019
  // Driver routine to generate observed fluxes and uncertainties
  // from true/generated mags.
  //
  // Dec 2019: begin refactor to allow off-diagonal covariances;
  //           e.g., correlations for anomalous host noise.
  //

  int NEPOCH = GENLC.NEPOCH ;
  int MEM    = (NEPOCH+1)*sizeof(FLUXNOISE_DEF);
  int epoch, icov;
  int VBOSE_CALC  = 0 ; 
  int VBOSE_FUDGE = 0 ;
  int VBOSE_APPLY = 0 ;
  char fnam[] = "GENFLUX_DRIVER" ;

  // -------------- BEGIN ---------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return; }

  GENLC.FLUXNOISE = (FLUXNOISE_DEF*) malloc(MEM);
  NGENFLUX_DRIVER++ ;

  // generate randoms for each epopch and filter
  // Avoid calling randoms twice when running both legacy gen_smearFlux
  // and refactored code here.
  gen_fluxNoise_randoms(); 

  for(icov=0; icov < NREDCOV_FLUXERRMODEL; icov++ )
    { COVINFO_FLUXERRMODEL[icov].NOBS = 0 ; }

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {

    GENLC.flux[epoch]         = NULLFLOAT ; 
    GENLC.fluxerr_data[epoch] = NULLFLOAT ;     
    GENLC.FLUXNOISE[epoch].IFILT_OBS = -888 ;

    if ( !GENLC.OBSFLAG_GEN[epoch]  )  { continue ; }
    gen_fluxNoise_calc(epoch,VBOSE_CALC, &GENLC.FLUXNOISE[epoch]);

    // check noise fudge-options; diagonal COV only
    gen_fluxNoise_fudge_diag(epoch, VBOSE_FUDGE, &GENLC.FLUXNOISE[epoch]);

    if ( VBOSE_CALC ) 
      { dumpLine_fluxNoise("NEW", epoch, &GENLC.FLUXNOISE[epoch] );  }

  }

  // check for optional flux covariance in FLUXERRMODEL_FILE maps
  for(icov=0; icov < NREDCOV_FLUXERRMODEL; icov++ )
    { gen_fluxNoise_fudge_cov(icov); }


  // apply random noise to each flux
  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {
    if ( !GENLC.OBSFLAG_GEN[epoch]  )  { continue ; }
    gen_fluxNoise_apply(epoch, VBOSE_APPLY, &GENLC.FLUXNOISE[epoch] );

    set_GENFLUX_FLAGS(epoch);
  }

  // monitor covariance separately for S,T,F components, and each band.
  if ( (INPUTS.FLUXERRMODEL_OPTMASK & MASK_MONITORCOV_FLUXERRMODEL)>0 )
    {  monitorCov_fluxNoise(); }

  free(GENLC.FLUXNOISE);

  return ;

} // end GENFLUX_DRIVER


// *****************************************
void gen_fluxNoise_randoms(void) {

  // Created Mar 1, 2014
  // Generate Gaussian randoms for each epoch (for flux noise)
  // and for each filter (coherent template noise).
  // These randoms are used in gen_smearFlux().
  //
  // Aug 27 2014: template randoms also depend on field-overlap
  //              Note that this breaks random sequence.
  //
  // Feb 14 2018: set GENLC.RANGauss_NOISE_ZP[ep] 
  //

  double RAN1, RAN2;
  int ep, ifilt, ifilt_obs, ifield;
  //  char fnam[] = "gen_fluxNoise_randoms" ;

  // -------------- BEGIN --------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return ; }

  // one random per epoch
  for ( ep = 1; ep <= GENLC.NEPOCH; ep++ )  {  

    ifilt_obs = GENLC.IFILT_OBS[ep] ;    
    GENLC.RANGauss_NOISE_SEARCH[ep] = -99999. ;  
    GENLC.RANGauss_NOISE_FUDGE[ep]  = -99999. ;  
    GENLC.RANGauss_NOISE_ZP[ep]     = -99999. ;  

    // skip un-used epochs so that randoms stay synced with
    // previous (10_33g) snana version.
    if ( !GENLC.OBSFLAG_GEN[ep]  ) { continue ; }

    // load randoms into global
    RAN1 = getRan_Gauss(1) ;  
    RAN2 = getRan_Gauss(1) ;  
    GENLC.RANGauss_NOISE_SEARCH[ep] = RAN1;
    GENLC.RANGauss_NOISE_ZP[ep]     = RAN2; // Jan 2020; soon to be obsolete
    GENLC.RANGauss_NOISE_FUDGE[ep]  = RAN2; // for refactored GENFLUX_DRIVER

  } // end ep loop


  // one random per filter (for template noise) and field overlap
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_SIMLIB; ifilt++ ) {
    ifilt_obs =  GENLC.IFILTMAP_SIMLIB[ifilt] ;

    // init to crazy values
    for(ifield=0; ifield < MXFIELD_OVP; ifield++ ) {      
      GENLC.RANGauss_NOISE_TEMPLATE[ifield][ifilt_obs] = -99999. ; 
    }

    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }
   
    for(ifield=0; ifield < MXFIELD_OVP; ifield++ ) {      
      GENLC.RANGauss_NOISE_TEMPLATE[ifield][ifilt_obs] = getRan_Gauss(1) ; 
    } 
    
  }

  return ;

}   // end of gen_fluxNoise_randoms


// *************************************
void gen_fluxNoise_calc(int epoch, int vbose, FLUXNOISE_DEF *FLUXNOISE) {

  // Created Dec 27 2019
  // Calculate Poisson noise (ie., sigma_flux) for epoch index 'ep'
  // and store calculated errors in FLUXNOISE struct.
  // Units are p.e.
  // Do not include error fudges here, and do not apply errors here.

  int ifilt_obs      = GENLC.IFILT_OBS[epoch] ;
  //  int ifilt          = GENLC.IFILTINVMAP_OBS[ifilt_obs] ; // sparse index

  double  mjd        = SIMLIB_OBS_GEN.MJD[epoch] ;
  double  pixsize    = SIMLIB_OBS_GEN.PIXSIZE[epoch] ;
  double  ccdgain    = SIMLIB_OBS_GEN.CCDGAIN[epoch] ;
  double  skysig     = SIMLIB_OBS_GEN.SKYSIG[epoch] ;
  double  readnoise  = SIMLIB_OBS_GEN.READNOISE[epoch] ;
  double  psfsig1    = SIMLIB_OBS_GEN.PSFSIG1[epoch] ; // pixels
  double  psfsig2    = SIMLIB_OBS_GEN.PSFSIG2[epoch] ;
  double  psfratio   = SIMLIB_OBS_GEN.PSFRATIO[epoch] ;
  double  nea        = SIMLIB_OBS_GEN.NEA[epoch] ;
  double  zpt        = SIMLIB_OBS_GEN.ZPTADU[epoch] ;
  double  zpterr     = SIMLIB_OBS_GEN.ZPTERR[epoch] ;
  double  template_skysig     = SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[epoch] ; 
  double  template_readnoise  = SIMLIB_OBS_GEN.TEMPLATE_READNOISE[epoch] ; 
  double  template_zpt        = SIMLIB_OBS_GEN.TEMPLATE_ZPT[epoch] ;

  double  genmag     = GENLC.genmag_obs[epoch];
  double  genmag_T   = GENLC.genmag_obs_template[ifilt_obs]; // for LCLIB

  int    OVP, NERR=0;

  double ZPTDIF_ADU, NADU_over_FLUXCAL, Npe_over_FLUXCAL, NADU_over_Npe ;
  double sqsig_noZ, sqsig_true, sqsig_data, sqsig_mon, flux_T, arg;
  double fluxsn_adu, fluxsn_pe, fluxmon_pe=0.0 ;
  double area_bg, psfsig_arcsec ;
  double sqerr_ccd_pe, skysig_tmp_pe, sqerr_sky_pe, sqerr_zp_pe;
  double fluxgal_pe, galmag=0.0 ;
  double template_sqerr_sky_pe, template_sqerr_ccd_pe, template_sqerr_pe;
  double SNR_CALC_S, SNR_CALC_ST, SNR_CALC_SZT, SQSIG_CALC, SNR_MON ;
  char band[2];
  char fnam[] = "gen_fluxNoise_calc" ;

  // ------------- begin ---------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return ; }

  sprintf(band, "%c", FILTERSTRING[ifilt_obs] );

  NERR=0;
  if ( zpt     < 10.0   ) { NERR++ ; }
  if ( psfsig1 < 0.0001 ) { NERR++ ; }
  if ( skysig  < 0.0001 ) { NERR++ ; } 
  if ( NERR > 0 ) {
    sprintf(c1err,"%d invalid observing conditions for ep=%d, band=%s",
	    NERR, epoch, band);
    sprintf(c2err,"mjd=%.3f zpt=%.2f, psf=%.3f, skysig=%.2f", 
	    mjd, zpt, psfsig1, skysig);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // compute zp in photo-electrons 
  ZPTDIF_ADU        = zpt - ZEROPOINT_FLUXCAL_DEFAULT ;
  NADU_over_FLUXCAL = pow( TEN , 0.4*ZPTDIF_ADU) ;
  Npe_over_FLUXCAL  = NADU_over_FLUXCAL * ccdgain;
  NADU_over_Npe     = NADU_over_FLUXCAL/Npe_over_FLUXCAL ;  

  flux_T   = 0.0 ;
  if ( genmag_T < 90.0 ) {
    arg     = 0.4 * ( zpt - genmag_T );
    flux_T  = pow(10.0,arg) / NADU_over_Npe ; 
  }


  // use search-run zero-point to convert mag -> flux.
  arg           = 0.4 * ( zpt - genmag );
  fluxsn_adu    = pow(10.0,arg);         // flux in ADUs
  fluxsn_pe     = fluxsn_adu * ccdgain ; // flux in pe
  
  // compute optional signal for monitor mag
  if ( INPUTS.MAGMONITOR_SNR > 10 ) {
    double magmon = (double)INPUTS.MAGMONITOR_SNR ;
    arg           = 0.4 * ( zpt - magmon );
    fluxmon_pe    = ccdgain * pow(10.0,arg); 
  }

  // get effective aperture  area (pixels)
  area_bg       = nea; // Feb 28 2021
  psfsig_arcsec = pixsize * sqrt(area_bg/(2.0*TWOPI))  ; 

  // get total sky noise for search run; includes sky & CCD noise
  skysig_tmp_pe = skysig * ccdgain ;  // convert ADU -> pe per pixel
  sqerr_sky_pe = area_bg * (skysig_tmp_pe*skysig_tmp_pe);
  sqerr_ccd_pe = area_bg * (readnoise*readnoise) ;  // in Npe

  // non-linearity ??

  // add sky-noise from template, integrated over effective aperture 
  template_sqerr_sky_pe = template_sqerr_ccd_pe = 0.0 ;
  template_sqerr_pe = 0.0 ;

  if ( SIMLIB_TEMPLATE.USEFLAG && template_skysig > 0.0 ) {

    if ( template_zpt < 10.0 ) {
      sprintf(c1err,"Invalid template_zpt(%c)=%f for  LIBID=%d at MJD=%.3f", 
	      FILTERSTRING[ifilt_obs], template_zpt, GENLC.SIMLIB_ID, mjd );
      sprintf(c2err,"Need TEMPLATE_ZPT to scale template noise.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    skysig_tmp_pe = template_skysig * ccdgain ;  // convert sigma in ADU -> pe.
    template_sqerr_sky_pe = area_bg * (skysig_tmp_pe*skysig_tmp_pe) ;
    template_sqerr_ccd_pe = area_bg * (template_readnoise*template_readnoise);

    // scale template noise to the search image
    double zparg, zfac;
    zparg = 0.8*(zpt - template_zpt); // Feb 2 2017 bug fix
    zfac  = pow(TEN, zparg);
    template_sqerr_sky_pe *= zfac ;
    template_sqerr_ccd_pe *= zfac ;
    template_sqerr_pe = template_sqerr_sky_pe + template_sqerr_ccd_pe ;

  }

  // -------------------------------------------
  // galaxy noise from photo-stats
  fluxgal_pe = 0.0;
  OVP = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_PHOT ;
  if ( OVP > 0 ) {
    // get galmag over NEA
    galmag        = interp_GALMAG_HOSTLIB(ifilt_obs, psfsig_arcsec );
    arg           = 0.4 * ( zpt - galmag );
    fluxgal_pe    = ccdgain * pow(10.0,arg);   // effec-aper flux in pe.
  }


  // check option ZP smearing 
  sqerr_zp_pe = 0.0 ;
  if ( INPUTS.SMEARFLAG_ZEROPT > 0 ) {
    double relerr, err;
    relerr  = pow(TEN, 0.4*zpterr) - 1.0 ;
    err     = (fluxsn_pe-flux_T) * relerr ;    
    sqerr_zp_pe = err*err;
  }

  // --------------------
  // add up SN flux error (photo-electrons^2) in quadrature
  // Do not include ZPerr nor correlated template noise here.
  sqsig_noZ
    = fluxsn_pe        // signal stat-error
    + fluxgal_pe       // square of host galaxy stat-error
    + sqerr_sky_pe     // sky-err from search run
    + sqerr_ccd_pe     // CCD read noise (added Dec 13, 2010)
    ;

  sqsig_true = sqsig_noZ ;
  sqsig_data = sqsig_noZ ;

  // Check options to include ZPerr in true and reported flux-error.
  // Note that correlated template noise is not included here.

  
  if ( (INPUTS.SMEARFLAG_ZEROPT & 1) > 0 ) 
    { sqsig_true += sqerr_zp_pe; }
  
  
  if ( (INPUTS.SMEARFLAG_ZEROPT & 2) > 0 ) 
    { sqsig_data += sqerr_zp_pe; }  // reported error includes zperr
  

  SNR_MON = 0.0 ;
  if ( INPUTS.MAGMONITOR_SNR > 10 ) {
    sqsig_mon = (sqsig_data - fluxsn_pe + fluxmon_pe + template_sqerr_pe);
    SNR_MON = fluxmon_pe / sqrt(sqsig_mon);
  }

  // calculated ERROR and SNR are for error fudges 
  SQSIG_CALC    = sqsig_data + template_sqerr_pe ;
  SNR_CALC_SZT  = fluxsn_pe / sqrt(sqsig_data + template_sqerr_pe);
  SNR_CALC_ST   = fluxsn_pe / sqrt(sqsig_noZ  + template_sqerr_pe);
  SNR_CALC_S    = fluxsn_pe / sqrt(sqsig_noZ); 

  // - - - - - - - - - - - - - - -   
  // load info in output structure

  FLUXNOISE->SQSIG_SRC       = fluxsn_pe ; // image source noise
  FLUXNOISE->SQSIG_TSRC      = flux_T ;    // template source noise (LCLIB)
  FLUXNOISE->SQSIG_SKY       = sqerr_sky_pe + sqerr_ccd_pe ;
  FLUXNOISE->SQSIG_TSKY      = template_sqerr_pe ; // template sky noise
  FLUXNOISE->SQSIG_ZP        = sqerr_zp_pe;  
  FLUXNOISE->SQSIG_HOST_PHOT = fluxgal_pe ;

  FLUXNOISE->SQSIG_CALC_TRUE[TYPE_FLUXNOISE_S]    = sqsig_noZ ;
  FLUXNOISE->SQSIG_CALC_TRUE[TYPE_FLUXNOISE_SZ]   = sqsig_true ;
  FLUXNOISE->SQSIG_CALC_TRUE[TYPE_FLUXNOISE_T]    = template_sqerr_pe;
  FLUXNOISE->SQSIG_CALC_TRUE[TYPE_FLUXNOISE_Z]    = sqerr_zp_pe;
  FLUXNOISE->SQSIG_CALC_TRUE[TYPE_FLUXNOISE_F]    = 0.0 ;


  FLUXNOISE->SQSIG_CALC_TRUE[TYPE_FLUXNOISE_SUM]  = 
    FLUXNOISE->SQSIG_CALC_TRUE[TYPE_FLUXNOISE_SZ] +
    FLUXNOISE->SQSIG_CALC_TRUE[TYPE_FLUXNOISE_T]  ;

  int itype;
  for(itype=0; itype < NTYPE_FLUXNOISE ; itype++ )  { 
    FLUXNOISE->SQSIG_FUDGE_TRUE[itype] = 0.0 ;
    FLUXNOISE->SQSIG_FINAL_TRUE[itype] = FLUXNOISE->SQSIG_CALC_TRUE[itype]; 
  }

  FLUXNOISE->SQSIG_CALC_DATA   = sqsig_data + template_sqerr_pe;
  FLUXNOISE->SQSIG_FUDGE_DATA  = 0.0 ;
  FLUXNOISE->SQSIG_FINAL_DATA  = sqsig_data + template_sqerr_pe;

  FLUXNOISE->SNR_CALC_S          = SNR_CALC_S ;
  FLUXNOISE->SNR_CALC_ST         = SNR_CALC_ST ;
  FLUXNOISE->SNR_CALC_SZT        = SNR_CALC_SZT ;

  FLUXNOISE->SNR_CALC_MON        = SNR_MON  ;
  FLUXNOISE->SNR_FINAL_MON       = SNR_MON  ;

  FLUXNOISE->NEA                 = area_bg ;
  FLUXNOISE->GALMAG_NEA          = galmag ;
  FLUXNOISE->Npe_over_FLUXCAL    = Npe_over_FLUXCAL;
  FLUXNOISE->NADU_over_Npe       = NADU_over_Npe ;

  FLUXNOISE->IFILT_OBS = ifilt_obs;
  sprintf(FLUXNOISE->BAND,"%s",band);

  return;

} // end gen_fluxNoise_calc

// ********************************************************
void  gen_fluxNoise_fudge_diag(int epoch, int VBOSE, FLUXNOISE_DEF *FLUXNOISE){

  // Created Dec 27, 2019
  // Compute diagonal error fudges, if specified 
  // (ignore off-diag correlations among epochs)
  //
  // Mar 12 2020:
  //  Fix bug applying INPUTS.FUDGESCALE_FLUXERR_FILTER(2) because
  //  it was using undefined SQSCALE.
  //
  // Apr 14 2021: abort of SQSIG_F<0 (happens if err scale < 1)

  int    ifilt_obs  = GENLC.IFILT_OBS[epoch] ;
  char   *FIELD     = GENLC.FIELDNAME[epoch];

  double  MJD       = SIMLIB_OBS_GEN.MJD[epoch] ;
  double  SKYSIG    = SIMLIB_OBS_GEN.SKYSIG[epoch] ;
  double  PSFSIG1   = SIMLIB_OBS_GEN.PSFSIG1[epoch] ; // pixels
  double  ZPADU     = SIMLIB_OBS_GEN.ZPTADU[epoch] ;
  double  PIXSIZE   = SIMLIB_OBS_GEN.PIXSIZE[epoch];

  long long GALID           = SNHOSTGAL.GALID ;
  double SBmag              = SNHOSTGAL.SB_MAG[ifilt_obs];
  double SNSEP              = SNHOSTGAL.SNSEP ;

  double Npe_over_FLUXCAL   = FLUXNOISE->Npe_over_FLUXCAL;
  double NEA                = FLUXNOISE->NEA;
  
  double SQSIG_CALC         = FLUXNOISE->SQSIG_CALC_DATA;
  double SIG_CALC           = sqrt(SQSIG_CALC);
  double SNR_CALC_SZT       = FLUXNOISE->SNR_CALC_SZT ;
  double SNR_CALC_ST        = FLUXNOISE->SNR_CALC_ST ;
  //  double SNR_CALC_S         = FLUXNOISE->SNR_CALC_S ;

  double SQSIG_SRC          = FLUXNOISE->SQSIG_SRC; // = SN flux, p.e.
  double SQSIG_SKY          = FLUXNOISE->SQSIG_SKY;
  double SQSIG_HOST_PHOT    = FLUXNOISE->SQSIG_HOST_PHOT;
  double GALMAG_NEA         = FLUXNOISE->GALMAG_NEA ;

  int    NTYPE = NTYPE_FLUXNOISE;
  int    OVP, itype ;
  double SCALE, SQSCALE, magerr_tmp, fluxErr_tmp ;
  double SQSIG_TRUE[NTYPE_FLUXNOISE];
  double SQSIG_DATA, SQSIG_TMP, SNR_MON, SQSIG_SCALED, SQSIG_F ;
  char   BAND[2];
  char   fnam[] = "gen_fluxNoise_fudge_diag" ;

  // ------------ BEGIN ----------

  sprintf(BAND, "%c", FILTERSTRING[ifilt_obs] );
  if ( SBmag > 32.0 ) { SBmag = 32.0; }    // to limit fluxerrmap size
  
  SQSIG_DATA = FLUXNOISE->SQSIG_FINAL_DATA; 
  for(itype=0; itype < NTYPE; itype++ ) { 
    SQSIG_TRUE[itype] = FLUXNOISE->SQSIG_FINAL_TRUE[itype]; 
  }

  // - - - - - - - - - - - - - - 

  // optional magErr fudge in quadrature (default=0)
  magerr_tmp = INPUTS.FUDGE_MAGERR_FILTER[ifilt_obs]; 
  if ( magerr_tmp > 1.0E-9 ) {
    fluxErr_tmp   = SQSIG_SRC* ( 1.0 - pow(10.0,-0.4*magerr_tmp)) ;
    SQSIG_TMP     = fluxErr_tmp * fluxErr_tmp ;
    SQSIG_TRUE[TYPE_FLUXNOISE_S]  += SQSIG_TMP ;
    SQSIG_DATA                    += SQSIG_TMP ;
  }

  // Optional errscale fudge (default=1) applied to true and reported errors.
  SCALE  = INPUTS.FUDGESCALE_FLUXERR_FILTER[ifilt_obs] ;
  if ( fabs(SCALE-1.0) > 1.0E-9 ) {
    // SCALE      = SCALE * SCALE ; xxx bug removed, Mar 12 2020
    SQSCALE      = SCALE * SCALE ;
    SQSIG_DATA  *= SQSCALE ;  
    for(itype=0; itype < NTYPE; itype++ ) { SQSIG_TRUE[itype] *= SQSCALE ; }
  }
  

  // Optional errscale fudge (default=1) applied only to reported errors
  SCALE  = INPUTS.FUDGESCALE_FLUXERR2_FILTER[ifilt_obs] ;
  if ( fabs(SCALE-1.0) > 1.0E-9 ) {
    SQSCALE        = SCALE * SCALE ;
    SQSIG_DATA    *= SQSCALE ;  
  }
  
 
  // Feb 2018: fudge error from FLUXERRMODEL. Should replace _legacy codes.
  if ( NMAP_FLUXERRMODEL > 0 ) {
    double ERRPARLIST[MXPAR_FLUXERRMAP];
    double LOGSNR = log10(SNR_CALC_SZT); 
    double PSF_FWHM = (PSFSIG1/PIXSIZE)*2.3548; // sigma(pix) -> FWHM(arcsec)
    int OPT = 0;
    if ( LOGSNR < -0.9 ) { LOGSNR = -0.9 ; }
    ERRPARLIST[IPAR_FLUXERRMAP_MJD]    = MJD;
    ERRPARLIST[IPAR_FLUXERRMAP_PSF]    = PSF_FWHM;  // FWHM, arcsec
    ERRPARLIST[IPAR_FLUXERRMAP_SKYSIG] = SKYSIG;    // ADU/pixel
    ERRPARLIST[IPAR_FLUXERRMAP_ZP]     = ZPADU;     // observed ZP, ADU
    ERRPARLIST[IPAR_FLUXERRMAP_LOGSNR] = LOGSNR ;
    ERRPARLIST[IPAR_FLUXERRMAP_SBMAG]  = SBmag ;
    ERRPARLIST[IPAR_FLUXERRMAP_GALMAG] = GALMAG_NEA ;
    ERRPARLIST[IPAR_FLUXERRMAP_SNSEP]  = SNSEP ;
    
    // pass FLUXCAL units to fluxErrModel in case of additive term.
    double FLUXCALERR_in  = SIG_CALC/Npe_over_FLUXCAL ;
    double FLUXCALERR_TRUE;  // generated error
    double FLUXCALERR_DATA ; // reported error in data file

    get_FLUXERRMODEL(OPT, FLUXCALERR_in, BAND, FIELD,      // (I)
		     NPAR_FLUXERRMAP_REQUIRE, ERRPARLIST,  // (I)
		     &FLUXCALERR_TRUE, &FLUXCALERR_DATA) ; // (O)
    
    SCALE   = FLUXCALERR_TRUE/FLUXCALERR_in ;  SQSCALE = SCALE*SCALE ;

    // don't modify true SQSIG(S); instead, add extra error to FUDGE term.
    SQSIG_TMP    = SQSIG_TRUE[TYPE_FLUXNOISE_S] ; 
    SQSIG_SCALED = SQSCALE * SQSIG_TMP ;
    SQSIG_F      = SQSIG_SCALED - SQSIG_TMP;
    SQSIG_TRUE[TYPE_FLUXNOISE_F] = SQSIG_F;

    if ( SQSIG_F < 0.0 ) {
      print_preAbort_banner(fnam);
      printf("  MJD  = %le \n", MJD);
      printf("  BAND = %s \n", BAND);
      printf("  PSFSIG1=%.3f  SBmag=%.3f\n", PSFSIG1, SBmag );

      sprintf(c1err,"Invalid SQSIG_FUDGE = %f < 0", SQSIG_F);
      sprintf(c2err,"fluxErrScale probably < 1");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    SQSIG_TRUE[TYPE_FLUXNOISE_T] *= SQSCALE ; // scale template noise

    SCALE   = FLUXCALERR_DATA/FLUXCALERR_in ;  SQSCALE = SCALE*SCALE ;
    SQSIG_DATA *= SQSCALE ;

    // keep track of NOBS per covariance matrix
    FLUXNOISE->INDEX_REDCOV = -9 ;
    if ( NREDCOV_FLUXERRMODEL > 0 ) {
      int ICOV = INDEX_REDCOV_FLUXERRMODEL(BAND,FIELD,2,fnam);      
      COVINFO_FLUXERRMODEL[ICOV].NOBS++ ;
      FLUXNOISE->INDEX_REDCOV = ICOV;
    }

  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //            BELOW ARE LEGACY OPTIONS
  // - - - - - - - - - - - - - - - - - - - - - - - - - - 


  // @@@@@@@@@@@@ LEGACY ERR-FUDGE FROM SIMLIB HEADER @@@@@@@@@@@@@@
  // Used only for SDSS.
  if ( SIMLIB_FLUXERR_COR.USE  ) {
    double  ERR_CAL, ERR_pe, XT, sqadderr_pe=0.0, SQSIG_ST_NEW, SQSIG_ST_ORIG;
    double  SQSIG_S = SQSIG_TRUE[TYPE_FLUXNOISE_S];
    double  SQSIG_T = SQSIG_TRUE[TYPE_FLUXNOISE_T];
    SQSIG_ST_ORIG = SQSIG_ST_NEW = (SQSIG_S + SQSIG_T);      
    ERR_CAL  = GENLC.SIMLIB_FLUXERR_ADDPAR[ifilt_obs];
    XT       = INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] ; 
    // translate ERR_CAL from FLUXCAL back to p.e.
    if ( ERR_CAL > 1.0E-6 && XT==1.0 ) {
      ERR_pe      = ERR_CAL * Npe_over_FLUXCAL ;
      sqadderr_pe = (ERR_pe * ERR_pe) ;
      SQSIG_TRUE[TYPE_FLUXNOISE_S]  += sqadderr_pe;
      SQSIG_TRUE[TYPE_FLUXNOISE_SZ] += sqadderr_pe;
      SQSIG_DATA                    += sqadderr_pe;
      SQSIG_ST_NEW                  += sqadderr_pe;
    }

    // to be consistent with legacy code, SNR_CALC here must
    // include sqadderr_pe
    double SNR_CALC = SNR_CALC_ST * sqrt(SQSIG_ST_ORIG/SQSIG_ST_NEW);
    SCALE     = get_SIMLIB_fluxerrScale_LEGACY(ifilt_obs, SNR_CALC);    
    SQSCALE   = SCALE * SCALE ;
    SQSIG_DATA *= SQSCALE ;
    for(itype=0; itype < NTYPE; itype++ ) { SQSIG_TRUE[itype] *= SQSCALE ; }

  }
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  // @@@@@@@@@@@@@@@@@ LEGACY ERRFUDGE @@@@@@@@@@@@@@@@@@@@@@@@@
  // Aug 2014: anomolous host-subtraction noise (HOSTNOISE_FILE)
  //           Note the dependence on both band and field.
  //  Should use newer FLUXERRMODEL_FILE
  OVP = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_IMAGE ;
  if ( OVP ) {    
    double noisePar[10];
    double sqImageNoise_pe = 0.0 ;
    double HOSTNOISE_pe = 0.0, HOSTNOISE_FLUXCAL = 0.0 ;
    double HOSTNOISE_ERRSCALE, SQ0, SQ1 ;

    GEN_NOISEMODEL_HOST_LEGACY(BAND,FIELD,(int)GALID,GALMAG_NEA,SBmag,SNSEP, 
			       noisePar);  // <== return this array

    HOSTNOISE_FLUXCAL  = noisePar[0] ; // add this noise, per pixel
    HOSTNOISE_ERRSCALE = noisePar[1] ; // scale on added sky noise

    // convert extra noise in FLUXCAL unit to Npe (per pixel)
    HOSTNOISE_pe    = HOSTNOISE_FLUXCAL * Npe_over_FLUXCAL ;
    SQ0 = HOSTNOISE_pe * HOSTNOISE_pe ;
    SQ1 = HOSTNOISE_ERRSCALE * HOSTNOISE_ERRSCALE ;

    sqImageNoise_pe  = 
      (NEA * SQ0)  +                       // quadrature model
      (SQSIG_SKY+SQSIG_HOST_PHOT)*(SQ1-1.0)  ;    // err-scale;

    // update TRUE error only; leave reported error as is.
    SQSIG_TRUE[TYPE_FLUXNOISE_S] += sqImageNoise_pe;    
  }
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




  // @@@@@@@@@@@@@@@@@ LEGACY ERROR SCALE @@@@@@@@@@@@@@@@@@@@@@@@@
  if ( INPUTS.FUDGEOPT_FLUXERR > 0 ) {
    SCALE  = scale_fluxErrModel_legacy(BAND,FIELD,MJD,ZPADU,SKYSIG,PSFSIG1);
    SQSCALE     = SCALE * SCALE ;
    SQSIG_DATA *= SQSCALE ;      // scale reported error only
  }  // @@@@@@@@@@@@@@@@@ LEGACY ERROR SCALE @@@@@@@@@@@@@@@@@@@@@@@@@


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  FLUXNOISE->SQSIG_FINAL_DATA = SQSIG_DATA;
  for(itype=0; itype < NTYPE; itype++ ) { 
    FLUXNOISE->SQSIG_FINAL_TRUE[itype] = SQSIG_TRUE[itype];
  }

  FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_SUM]  = 
    FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_S] +
    FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_T] +
    FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_F] ; 


  // store monitor SNR, corrected for error fudges
  SNR_MON = FLUXNOISE->SNR_CALC_MON ;
  if ( SNR_MON > 1.0E-9 ) {
    SQSCALE = FLUXNOISE->SQSIG_FINAL_DATA/FLUXNOISE->SQSIG_CALC_DATA ;
    FLUXNOISE->SNR_FINAL_MON = SNR_MON / sqrt(SQSCALE) ;
  }

  return ;

} // end gen_fluxNoise_fudge_diag

// ******************************
void gen_fluxNoise_fudge_cov(int icov) {

  // Created Jan 10 2020
  // Compute correlated randoms among epochs in "icov" group.
  // Modify GAURAN_F = GENLC.RANGauss_NOISE_FUDGE[epoch] so that
  //     GAURAN_F * sqrt(SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_F])
  //
  // is the correlated flux-shift.
  //
  // Beware that epochs in "icov" group are a subset of the
  // GENLC.NEPOCH total epochs, so watch indices.

  int  NOBS = COVINFO_FLUXERRMODEL[icov].NOBS ;
  int  MEMD0 = NOBS*sizeof(double);
  int  MEMD1 = NOBS*sizeof(double*);
  int  NEPOCH  = GENLC.NEPOCH ;
  //  int  TYPE_S  = TYPE_FLUXNOISE_S ;  // search noise
  int  TYPE_F  = TYPE_FLUXNOISE_F ;  // fudge noise
  //  double COVFAC = 1.0E-3;

  int  ep, iep0, iep1, indx_1D, IFILT_OBS, INDEX_REDCOV ;
  int  N0=0, N1=0, o, *epMAP ;
  double *covFlux_1D, **covCholesky_2D ;
  double SQSIG_FUDGE[2] ;
  double SIGxSIG, REDCOV, COV ;
  int LDMP = 0 ; 
  char fnam[] = "gen_fluxNoise_fudge_cov" ;

  // --------- BEGIN --------

  if ( NOBS <= 1 ) { return ; }

  if ( LDMP ) 
    { printf("\n xxx ------------- DUMP on for %s --------------- \n", fnam); }

  // create matrix for each cov matrix

  covFlux_1D     = (double*)malloc (NOBS*NOBS * sizeof(double) ) ;
  covCholesky_2D = (double**)malloc(MEMD1) ;
  for(o=0; o < NOBS; o++ ) { covCholesky_2D[o] = (double*) malloc(MEMD0); }

  epMAP = (int*) malloc( NOBS * sizeof(int) );

  // - - - - - - 
  for(iep0=1; iep0 <= NEPOCH; iep0++ ) {

    if ( !GENLC.OBSFLAG_GEN[iep0]  ) { continue ; }
    IFILT_OBS    = GENLC.IFILT_OBS[iep0] ;
    INDEX_REDCOV = GENLC.FLUXNOISE[iep0].INDEX_REDCOV;
    if ( INDEX_REDCOV != icov ) { continue; }
    N0++ ;  N1=0;

    epMAP[N0-1] = iep0 ; // preserve mapping between GENLC and COV arrays

    for(iep1=1; iep1 <= iep0; iep1++ ) {

      if ( !GENLC.OBSFLAG_GEN[iep1]  ) { continue ; }

      IFILT_OBS = GENLC.IFILT_OBS[iep1] ;
      INDEX_REDCOV = GENLC.FLUXNOISE[iep1].INDEX_REDCOV;
      if ( INDEX_REDCOV != icov ) { continue; }
      N1++ ;

      REDCOV  = COVINFO_FLUXERRMODEL[INDEX_REDCOV].REDCOV;
      if ( N0 == N1 ) { REDCOV = 1.0 ; }

      SQSIG_FUDGE[0] = GENLC.FLUXNOISE[iep0].SQSIG_FINAL_TRUE[TYPE_F];
      SQSIG_FUDGE[1] = GENLC.FLUXNOISE[iep1].SQSIG_FINAL_TRUE[TYPE_F];
      SIGxSIG  = sqrt(SQSIG_FUDGE[0] * SQSIG_FUDGE[1]) ;
      COV      = SIGxSIG * REDCOV ;
            
      indx_1D = (N1-1)*NOBS + (N0-1) ;
      covFlux_1D[indx_1D] = COV;

      indx_1D = (N0-1)*NOBS + (N1-1) ;
      covFlux_1D[indx_1D] = COV;

      if ( LDMP ) {
	printf(" xxx iep[0,1]=%2d,%2d  N[0,1]=%2d,%2d  indx_1D=%3d  "
	       "COV=%.1f\n",
	       iep0, iep1, N0, N1, indx_1D, COV); fflush(stdout);
	fflush(stdout);
      }


    }  // end iep2
  }  // end iep


  // - - - - - - - - 
  // sanity checks
  if ( N0 != NOBS || N1 != NOBS ) {
    sprintf(c1err,"Invalid %d x %d matrix; expected %d^2 ",
	    N0, N1, NOBS);
    sprintf(c2err,"CID=%d  NEPOCH=%d  icov=%d(%s)",
	    GENLC.CID, GENLC.NEPOCH, 
	    icov, COVINFO_FLUXERRMODEL[icov].BANDSTRING);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  if(LDMP) { dumpCovMat_fluxNoise(icov, NOBS, covFlux_1D); }

  // - - - - - - - - - - - - - - - - - - - - -
  // do Cholesky decomp ...
  int obs0, obs1;
  double *flux_scatter, GAURAN, GAURAN_ORIG, GAURAN_NEW, SIG_F ;
  gsl_matrix_view chk;

  flux_scatter = (double*) malloc( MEMD0 );
  chk  = gsl_matrix_view_array ( covFlux_1D, NOBS, NOBS); 
  gsl_linalg_cholesky_decomp ( &chk.matrix )  ;
  
  for(obs0=0; obs0 < NOBS; obs0++ ) {  
    flux_scatter[obs0] = 0.0 ;

    for(obs1=0; obs1 < NOBS; obs1++ ) {
      covCholesky_2D[obs0][obs1] = 0.0 ;
      if ( obs0 <= obs1 ) {
	covCholesky_2D[obs0][obs1] = gsl_matrix_get(&chk.matrix,obs0,obs1); 
      }

      covCholesky_2D[obs1][obs0] = 0.0 ;
      if ( obs1 <= obs0 ) {
	covCholesky_2D[obs1][obs0] = gsl_matrix_get(&chk.matrix,obs1,obs0); 
      }
     
      ep = epMAP[obs1];  GAURAN  = GENLC.RANGauss_NOISE_FUDGE[ep] ; 
      flux_scatter[obs0] += (GAURAN * covCholesky_2D[obs1][obs0]) ;
    }
  }

  
  // modify independent RANGauss_FUDGE
  for(o=0; o < NOBS; o++ ) {
    ep             = epMAP[o];
    SIG_F          = sqrt(GENLC.FLUXNOISE[ep].SQSIG_FINAL_TRUE[TYPE_F]);
    GAURAN_ORIG    = GENLC.RANGauss_NOISE_FUDGE[ep] ; 
    GAURAN_NEW     = flux_scatter[o]/SIG_F ;
    GENLC.RANGauss_NOISE_FUDGE[ep] = GAURAN_NEW ;
  } // end o loop


  /* xxxxxx
  for (i = 0 ; i < 3 ; i++) {
    for (j = 0 ; j < 3 ; j++) {
      double tmp = INPUTS.COVMAT_CHOLESKY_SCATTER[j][i]; 
      SCATTER_VALUES[i] += tmp * normalvector[j];      
    }
  }
  xxxxxxxxx*/

  // free matrix memory.
  free(covFlux_1D);

  for(o=0; o < NOBS; o++ ) { free(covCholesky_2D[o]); }
  free(covCholesky_2D);

  free(epMAP);  free(flux_scatter);

  return ;

} // end of  gen_fluxNoise_fudge_cov

// *********************a****************
void gen_fluxNoise_apply(int epoch, int vbose, FLUXNOISE_DEF *FLUXNOISE) {
  
  // Created Dec 2019
  // apply random flux shifts 

  char  *FIELD     = GENLC.FIELDNAME[epoch] ;

  int   ifilt_obs    = FLUXNOISE->IFILT_OBS;
  char  *BAND        = FLUXNOISE->BAND;
  double fluxTrue    = FLUXNOISE->SQSIG_SRC ;          // p.e.
  double fluxgal     = FLUXNOISE->SQSIG_HOST_PHOT ;
  
  double SQSIG_S     = FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_S];
  double SQSIG_SZ    = FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_SZ];
  double SQSIG_T     = FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_T];
  double SQSIG_F     = FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_F];
  //  double SQSIG_Z     = FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_Z];
  double SIG_S       = sqrt(SQSIG_S);  // search image
  double SIG_SZ      = sqrt(SQSIG_SZ);  // search image + zp err
  double SIG_T       = sqrt(SQSIG_T);  // template
  double SIG_F       = sqrt(SQSIG_F);  // fudge noise
  //  double SIG_Z       = sqrt(SQSIG_Z);  // ZP noise, included in _S

  int    ifield, OVP ;
  double GAURAN_SZ=0.0, GAURAN_T=0.0, GAURAN_F=0.0, GAURAN_S=0.0 ;
  double SHIFT_SZ,  SHIFT_T,  SHIFT_F, SHIFT_S ;
  double fluxObs, flux_T ;
  double SQSIG_TMP, SIG_TMP, SCALE_TMP, SNR_CALC ;
  char fnam[] = "gen_fluxNoise_apply" ;

  // ----------- BEGIN -----------

  if ( ifilt_obs < 0 || ifilt_obs > MXFILTINDX ) {
    sprintf(c1err,"Undefined IFILT_OBS=%d for ep=%d CID=%d .",
	    ifilt_obs, epoch, GENLC.CID);
    sprintf(c2err,"Probably a code bug.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  ifield        = IFIELD_OVP_SIMLIB(1,FIELD);
  if ( ifield < 0 ) { ifield=0; }   // Nov 2016: needed for GAURAN_TEMPLATE


  // strip off previuously generated Gaussian randoms
  GAURAN_SZ = GAURAN_F = GAURAN_T = 0.0 ;
  if ( INPUTS.SMEARFLAG_FLUX > 0 ) {     
    GAURAN_S   = GENLC.RANGauss_NOISE_SEARCH[epoch] ; 
    GAURAN_SZ  = GENLC.RANGauss_NOISE_SEARCH[epoch] ; 
    GAURAN_F   = GENLC.RANGauss_NOISE_FUDGE[epoch] ; 
    //    GAURAN_Z   = GENLC.RANGauss_NOISE_ZP[epoch] ; 
    GAURAN_T   = GENLC.RANGauss_NOISE_TEMPLATE[ifield][ifilt_obs] ;
  }

  SHIFT_SZ = SIG_SZ * GAURAN_SZ ; // independent part of search image
  SHIFT_S  = SIG_S  * GAURAN_S  ; // for monitor only
  SHIFT_T  = SIG_T  * GAURAN_T ;  // 100% correlated noise from template
  SHIFT_F  = SIG_F  * GAURAN_F ;  // fudged noise, maybe with correlations

  // store each shift to monitor later.
  FLUXNOISE->FLUX_SHIFT_TRUE[TYPE_FLUXNOISE_S]  = SHIFT_S ;
  FLUXNOISE->FLUX_SHIFT_TRUE[TYPE_FLUXNOISE_SZ] = SHIFT_SZ ;
  FLUXNOISE->FLUX_SHIFT_TRUE[TYPE_FLUXNOISE_T]  = SHIFT_T ;
  FLUXNOISE->FLUX_SHIFT_TRUE[TYPE_FLUXNOISE_F]  = SHIFT_F ;

  // - - - - - - - - - - -  -
  // Sum the shifts to get observed flux. Note that 
  //   * SIG_F != 0 only for FLUXERRMODEL_FILE option.
  fluxObs  = fluxTrue + (SHIFT_SZ + SHIFT_T + SHIFT_F) ;

  // - - - - - - - - - - -  -


  // check option for random template noise instead of default correlated noise
  OVP = (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_RANDOM_TEMPLATENOISE);
  if ( OVP ) {
    SQSIG_TMP = SQSIG_SZ + SQSIG_T ;
    SIG_TMP   = sqrt(SQSIG_TMP)   ;
    fluxObs   = fluxTrue + (SIG_TMP*GAURAN_SZ);
  }

  // Adjust reported error to be based on observed flux instead of true flux.  
  double sqerr_ran ;
  if ( fluxObs > 0 ) 
    { sqerr_ran = (fluxObs - fluxTrue); }
  else
    { sqerr_ran = -fluxTrue; }
  
  // update reported error in data file
  double SIG_FINAL_TRUE  = sqrt(FLUXNOISE->SQSIG_FINAL_DATA); // local use

  SQSIG_TMP = FLUXNOISE->SQSIG_FINAL_DATA + sqerr_ran;
  FLUXNOISE->SQSIG_FINAL_DATA = SQSIG_TMP ;
  FLUXNOISE->SIG_FINAL_DATA   = sqrt(SQSIG_TMP) ;
  FLUXNOISE->SQSIG_RAN        = sqerr_ran ;



  // check option to ignore source & host error in reported error
  // (SMP-like)
  if ( (INPUTS.SMEARFLAG_FLUX & 2) > 0 ) {
    SQSIG_TMP  = FLUXNOISE->SQSIG_FINAL_DATA - (fluxTrue + fluxgal) ;
    FLUXNOISE->SQSIG_FINAL_DATA  = SQSIG_TMP ;
    FLUXNOISE->SIG_FINAL_DATA    = sqrt(SQSIG_TMP) ;    
  }


  // --------------------------------------------
  // Check optional template flux to subtract (for LCLIB model).
  // Beware that coherent template fluctuations are not included,
  // so deep templates are assumed.
  // This template-flux subtraction is done at the very end so that
  // search-soure noise is included.
  flux_T   = FLUXNOISE->SQSIG_TSRC ;
  if ( flux_T > 1.0E-9 ) {
    fluxObs       -= flux_T ;  // obs flux; can be pos or neg
    fluxTrue      -= flux_T ;  // true flux without fluctuations

    // update SNR_CALC
    SCALE_TMP      = ( fluxTrue / ( fluxTrue + flux_T) ) ;
    FLUXNOISE->SNR_CALC_SZT *= SCALE_TMP ;
    FLUXNOISE->SNR_CALC_ST  *= SCALE_TMP ;
    FLUXNOISE->SNR_CALC_S   *= SCALE_TMP ;
  }


  // - - - - - - - - - - - - - - - - - - - 
  // Jan 2018: check for saturation. NPE > 0 --> saturation
  int npe_above_sat = npe_above_saturation(epoch,fluxTrue+fluxgal);
  GENLC.npe_above_sat[epoch] = npe_above_sat ;
  if ( npe_above_sat > 0 ) {
    fluxObs = 0.0 ;     SIG_TMP = FLUXCALERR_SATURATE ;
    FLUXNOISE->SQSIG_FINAL_DATA = SIG_TMP * SIG_TMP ;
    FLUXNOISE->SIG_FINAL_DATA   = SIG_TMP;
  }

  // ---------------------------------------------
  // load global GENLC array.
  // ---------------------------------------------


  double NADU_over_Npe       = FLUXNOISE->NADU_over_Npe;
  //  double Npe_over_FLUXCAL    = FLUXNOISE->Npe_over_FLUXCAL;
  double legacy_flux         = GENLC.flux[epoch];
  double legacy_fluxerr_data = GENLC.fluxerr_data[epoch];

  GENLC.flux[epoch]         = fluxObs * NADU_over_Npe ;  
  GENLC.fluxerr_data[epoch] = FLUXNOISE->SIG_FINAL_DATA * NADU_over_Npe;


  // store true SNR without fluctuations (used later to force SNRMAX)
  GENLC.trueSNR[epoch] =  fluxTrue/SIG_FINAL_TRUE;

  // store coherent template error.
  GENLC.template_err[epoch] = NADU_over_Npe * SIG_T ;

  // store SNR of fixed monitor mag
  GENLC.SNR_MON[epoch]  = FLUXNOISE->SNR_FINAL_MON ;
  
  // keep track of epoch with max SNR (Jun 2018)
  SNR_CALC              = FLUXNOISE->SNR_CALC_SZT; 
  GENLC.SNR_CALC[epoch] = SNR_CALC ;
  if (SNR_CALC > GENLC.SNRMAX_GLOBAL) 
    { GENLC.SNRMAX_GLOBAL = SNR_CALC;  GENLC.IEPOCH_SNRMAX_GLOBAL = epoch;  }

  //  if (SNR_CALC > GENLC.SNRMAX) 
  //{ GENLC.SNRMAX_GLOBAL = SNR_CALC;  GENLC.IEPOCH_SNRMAX_GLOBAL = epoch;  }

  
  if ( vbose ) {
    double flux         = GENLC.flux[epoch];
    double fluxerr_data = GENLC.fluxerr_data[epoch];
    double ratio_flux   = flux/legacy_flux;
    double ratio_err    = fluxerr_data/legacy_fluxerr_data;
    double ratio_tol    = 0.005;
    char starFlux[2]=" ", starErr[2]=" " ;
    if ( fabs(ratio_flux-1.0)>ratio_tol ) { sprintf(starFlux,"*"); }
    if ( fabs(ratio_err -1.0)>ratio_tol ) { sprintf(starErr, "*"); }

    printf(" xxx %s(%3d-%s) NEW/OLD flux=%7.2f/%7.2f=%7.4f%s  "
	   "err=%6.2f/%6.2f=%.4f%s\n"
	   ,"apply", epoch, BAND
	   ,flux, legacy_flux, ratio_flux, starFlux
	   ,fluxerr_data, legacy_fluxerr_data, ratio_err, starErr );
    fflush(stdout);
	   
  }

  // check for really crazy flux values
  check_crazyFlux(epoch, FLUXNOISE);
  
  if ( epoch == -7 ) 
    { dumpEpoch_fluxNoise_apply(fnam,epoch,FLUXNOISE); }


  return ;


} // end gen_fluxNoise_apply


// ********************************************
void  check_crazyFlux(int ep, FLUXNOISE_DEF *FLUXNOISE) {

  // Jan 2020
  // Abort if flux (in ADU) is way too large (i.e., crazy)
  //   [part of GENFLUX_DRIVER refactor] 

  int     ifilt_obs    = FLUXNOISE->IFILT_OBS;  
  double  ZPADU        = SIMLIB_OBS_GEN.ZPTADU[ep] ;
  double  zsn          = GENLC.REDSHIFT_HELIO ;
  double  flux         =  GENLC.flux[ep];         // ADU
  double  fluxerr      =  GENLC.fluxerr_data[ep]; // ADU
  double  mag_smear    =  GENLC.magsmear8[ep];

  double arg, pow_arg, xt, crazyFlux, crazyFlux_neg ;
  char fnam[] = "check_crazyFlux" ;

  // ----------- BEGIN -----------

  // use 1/z^2 dependence on mag to set bounds for crazy flux abort;
  // account for exposure time (xt) and SIMLIB zeropoint (zptfac)
  // Also add 10 sigma of noise to allow for fluctuations.

  xt = INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] ;
  if ( xt < 1.0 ) { xt = 1.0 ; }

  arg        = 0.4 * ( ZPADU - 31.0 );  
  pow_arg     = pow(10.0,arg);  
  if ( zsn > 1.0E-9 ) 
    { crazyFlux  = (4.E4 * pow_arg * xt) / (zsn*zsn) ; }
  else
    { crazyFlux = 1.0E14 ; } // for LCLIB (July 2018)

  crazyFlux += (10.*fluxerr);


  if ( mag_smear < 0.0 )  { // adjust for intrinsic smearing
    arg        = -0.4*mag_smear;  pow_arg = pow(TEN,arg); 
    crazyFlux *= pow_arg; 
  }


  if ( GENLC.SL_MAGSHIFT < 0.0 ) { // adjust for strong lens magnification
    arg = -0.4*GENLC.SL_MAGSHIFT ;  pow_arg = pow(TEN,arg); 
    crazyFlux *= pow_arg; 
  }



  if ( GENLC.FUDGE_SNRMAX_FLAG == 2 && INPUTS.FUDGE_SNRMAX > 1.0 ) 
    { crazyFlux *= INPUTS.FUDGE_SNRMAX; }

  if ( INDEX_GENMODEL == MODEL_SIMSED ) 
    { crazyFlux *= 10.0; }    // allow for really bright objects (Aug 2017)

  if ( INDEX_GENMODEL == MODEL_LCLIB ) 
    { crazyFlux *= 100.0; }  


  // determine NEGATIVE crazy flux 
  crazyFlux_neg = -crazyFlux ; // Mar 20 2021

  // - - - - - - - - - - - - - 
  if ( flux > crazyFlux || flux < crazyFlux_neg ) {
    print_preAbort_banner(fnam);
    dumpEpoch_fluxNoise_apply(fnam, ep, FLUXNOISE);
    if ( flux > 0.0 ) {
      sprintf(c1err, "flux=%le exceeds CRAZYFLUX=%le", flux, crazyFlux); 
    }
    else {
      sprintf(c1err, "flux=%le is below CRAZYFLUX=%le", flux, crazyFlux_neg); 
    }

    sprintf(c2err, "See dumpEpoch details");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  return ;

} // end check_crazyFlux


// ********************************************************
void dumpLine_fluxNoise(char *fnam, int ep, FLUXNOISE_DEF *FLUXNOISE) {
  // Dump util for refactoring GENFLUX_DRIVER
  double sqsig_calc       = FLUXNOISE->SQSIG_CALC_DATA ;
  double sqsig_final_true = FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_SUM] ;
  double sqsig_final_data = FLUXNOISE->SQSIG_FINAL_DATA ;
  double sqsig_src        = FLUXNOISE->SQSIG_SRC ;
  char  *band             = FLUXNOISE->BAND;
  // --------- begin --------
  printf(" xxx %s: sig(%3d-%s) = %9.3f -> %9.3f(true)/%9.3f(data) "
	 " F_SN=%.1f\n",
	 fnam, ep, band,
	 sqrt(sqsig_calc),
	 sqrt(sqsig_final_true), 
	 sqrt(sqsig_final_data), 
	 sqsig_src ); 

  return;

} // end dumpLine_fluxNoise

// ***********************
void dumpCovMat_fluxNoise(int icov, int NOBS, double *COV) {

  int i0, i1, J1D, NOBS_dump = NOBS;
  double COVTMP ;
  //  char fnam[] = "dumpCovMat_fluxNoise" ;
  if ( NOBS_dump > 8 ) { NOBS_dump = 8; }
  printf("\n Dump COV matrix for icov=%d  (%s) \n", 
	 icov, COVINFO_FLUXERRMODEL[icov].BANDSTRING);
  
  for(i0=0; i0 < NOBS_dump; i0++ ) {
    for(i1=0; i1 < NOBS_dump; i1++ ) {
      J1D = i0*NOBS + i1;
      COVTMP = COV[J1D];
      printf("%10.1f ", COVTMP); fflush(stdout);
    }
    printf("\n");
  }

  //  debugexit(fnam);

} // end dumpCovMat_fluxNoise



// *********************************************
void monitorCov_fluxNoise(void) {

  // Jan 2020
  // evaluate reduced correlation for each band,
  // and separately for Search, Template, Fudge components.
  // REDCORs are computed for each event, and then averaged
  // over all events. 

  int NEPOCH = GENLC.NEPOCH ;
  //  int TYPE_T = TYPE_FLUXNOISE_T ;

#define NTYPE_CHECK 3
  int TYPE_LIST[NTYPE_CHECK];
  char TYPE_STRING[] = "STF";
  int ifilt_obs, ifilt, ep0, ep1, NSUM, i, TYPE ;
  double SHIFT[2], SIG[2], RHO, RHO_SUM, RHO_AVG;
  //  char fnam[] = "monitorCov_fluxNoise" ;

  // ------------ BEGIN -----------

  TYPE_LIST[0] = TYPE_FLUXNOISE_SZ ;
  TYPE_LIST[1] = TYPE_FLUXNOISE_T ;
  TYPE_LIST[2] = TYPE_FLUXNOISE_F ;

  // - - - - - - - - - - - - - - - - - - - - -
  for(ep0=1; ep0 <= NEPOCH; ep0++ ) {
    if ( !GENLC.OBSFLAG_GEN[ep0] ) { continue; }

    for(ep1=1; ep1 < ep0; ep1++ ) {
      if ( !GENLC.OBSFLAG_GEN[ep1] ) { continue; }

      if ( GENLC.IFILT_OBS[ep0] != GENLC.IFILT_OBS[ep1] ) { continue; }
      ifilt_obs = GENLC.IFILT_OBS[ep0];

      for(i=0; i < NTYPE_CHECK; i++ ) {
	TYPE = TYPE_LIST[i];
	SHIFT[0]  = GENLC.FLUXNOISE[ep0].FLUX_SHIFT_TRUE[TYPE];
	SHIFT[1]  = GENLC.FLUXNOISE[ep1].FLUX_SHIFT_TRUE[TYPE];
	SIG[0]    = sqrt(GENLC.FLUXNOISE[ep0].SQSIG_FINAL_TRUE[TYPE]) ;
	SIG[1]    = sqrt(GENLC.FLUXNOISE[ep1].SQSIG_FINAL_TRUE[TYPE]) ;
	RHO       = (SHIFT[0]*SHIFT[1]) / (SIG[0]*SIG[1]) ;
      
	if ( SIG[0] < 1.0E-9 ) { continue; }
	if ( SIG[1] < 1.0E-9 ) { continue; }

	/*
	  printf(" xxx %s: ifilt_obs=%d SHIFT=%7.3f,%7.3f  SIG=%.3f,%.3f\n", 
	  fnam, ifilt_obs, SHIFT[0], SHIFT[1], SIG[0], SIG[1] );
	*/

	GENLC.MONITOR_REDCOV_FLUXNOISE[ifilt_obs][TYPE].RHO_SUM += RHO;
	GENLC.MONITOR_REDCOV_FLUXNOISE[ifilt_obs][TYPE].NSUM++ ;
      }
    }
  } // end ep
   

  // -------------------
  char cfilt[2];

  if ( (NGENFLUX_DRIVER % 50) != 0 ) { return; }
  
  printf(" xxx \n");
  printf(" xxx ------------ NCALL = %d------------------------- \n",
	 NGENFLUX_DRIVER );

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;

    if ( !GENLC.DOFILT[ifilt_obs] ) { continue; }
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
    printf(" xxx RHOFLUX(%s) = ", cfilt );

    for(i=0; i < NTYPE_CHECK; i++ ) {
      TYPE    = TYPE_LIST[i];
      RHO_SUM = GENLC.MONITOR_REDCOV_FLUXNOISE[ifilt_obs][TYPE].RHO_SUM;
      NSUM    = GENLC.MONITOR_REDCOV_FLUXNOISE[ifilt_obs][TYPE].NSUM ;
      RHO_AVG = RHO_SUM / (double)NSUM;
      printf(" %6.3f(%c)", RHO_AVG, TYPE_STRING[i] );
    } // end NTYPE_CHECK loop
    printf("\n"); fflush(stdout);

  } // end ifilt

  return ;

} // end monitorCov_fluxNoise(int ifilt_obs) {

// *****************************************
void dumpEpoch_fluxNoise_apply(char *fnam, int ep, FLUXNOISE_DEF *FLUXNOISE) {

  // Created Dec 2019
  // Complete fluxNoise dump for epoch 'ep'.
  // This is part of the GENFLUX_DRIVER refactor.

  
  int  ifilt_obs = FLUXNOISE->IFILT_OBS;
  char *band     = FLUXNOISE->BAND;
  double NADU_over_Npe       = FLUXNOISE->NADU_over_Npe;
  double Npe_over_FLUXCAL    = FLUXNOISE->Npe_over_FLUXCAL;

  double flux_data    = GENLC.flux[ep]/NADU_over_Npe;   // Npe
  double flux_true    = FLUXNOISE->SQSIG_SRC;           // Npe

  double SQSIG_DATA   = FLUXNOISE->SQSIG_FINAL_DATA;  // Npe
  double SQSIG_TRUE   = FLUXNOISE->SQSIG_FINAL_TRUE[TYPE_FLUXNOISE_SUM];
  double fluxerr_data = sqrt(SQSIG_DATA);
  double fluxerr_true = sqrt(SQSIG_TRUE);

  double genmag       = GENLC.genmag_obs[ep];
  double Trest        = GENLC.epoch_rest[ep]; 
  double Tobs         = GENLC.epoch_obs[ep]; 

  char fnam_local[] = "dumpEpoch_fluxNoise_apply";

  // ----------- BEGIN --------------

  printf("\n");
  printf(" xxx ------------------------------------------------------- \n");
  printf(" xxx %s called from %s \n", fnam_local, fnam);

  printf(" xxx CID=%d  MJD=%.3f  ifilt_obs=%d(%s)  LIBID=%d\n",
	 GENLC.CID, GENLC.MJD[ep], ifilt_obs, band, GENLC.SIMLIB_ID );

  printf(" xxx z=%.3f  mu=%.3f  Trest=%6.2f  Tobs=%.2f  genmag(%c)=%6.1f  \n", 
	 GENLC.REDSHIFT_CMB, GENLC.DLMU, Trest, Tobs, 
	 FILTERSTRING[ifilt_obs], genmag );
  
  printf(" xxx GEN(AV,RV) = %7.3f , %7.3f  SHAPEPAR=%7.3f  "
	 "(c=%7.3f,mB=%.3f)\n", 
	 GENLC.AV, GENLC.RV, GENLC.SHAPEPAR, GENLC.SALT2c, GENLC.SALT2mB );

  if ( GENFRAME_OPT  == GENFRAME_REST ) {
    printf(" xxx Kcor  %s = %le   AVwarp=%7.3f\n",
	   GENLC.kcornam[ep], GENLC.kcorval8[ep], GENLC.AVwarp8[ep] );
  }
  else if ( INDEX_GENMODEL  == MODEL_SIMSED ) {
      printf(" xxx SIMSED PARAMS %s,%s = %f, %f  (x0=%le)\n"
	     ,INPUTS.PARNAME_SIMSED[0]
	     ,INPUTS.PARNAME_SIMSED[1]
	     ,GENLC.SIMSED_PARVAL[0]
	     ,GENLC.SIMSED_PARVAL[1], GENLC.SALT2x0 );
  }

  else if ( INDEX_GENMODEL == MODEL_LCLIB ) {
    printf(" xxx LCLIB EVENT ID = %lld \n", LCLIB_EVENT.ID);
  }

  printf(" xxx \n");    fflush(stdout);
  // - - - - - - - - - - - -  
  printf(" xxx FLUXNPE(obs ) = %le +_  %le\n",  flux_data, fluxerr_data );
  printf(" xxx FLUXNPE(true) = %le +_  %le\n",  flux_true, fluxerr_true );

  printf(" xxx FLUXCAL(obs ) = %le, +_ %le \n", 
	 flux_data/Npe_over_FLUXCAL,  fluxerr_data/Npe_over_FLUXCAL );
  printf(" xxx FLUXCAL(true) = %le, +_ %le \n", 
	 flux_true/Npe_over_FLUXCAL,  fluxerr_true/Npe_over_FLUXCAL );

  printf(" xxx SKYSIG(S,T) = %.3f, %.3f ADU/pix   PSFSIG = %.3f pixels \n",
	 SIMLIB_OBS_GEN.SKYSIG[ep], SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[ep],
	 SIMLIB_OBS_GEN.PSFSIG1[ep] );
  printf(" xxx ZPT(S,T) = %.3f, %.3f  (ADU)   GAIN=%f\n",	
	 SIMLIB_OBS_GEN.ZPTADU[ep], SIMLIB_OBS_GEN.TEMPLATE_ZPT[ep],
	 SIMLIB_OBS_GEN.CCDGAIN[ep] );

  printf(" xxx GAURAN(S,T,ZP) = %f, %f, %f \n",	 
	 GENLC.RANGauss_NOISE_SEARCH[ep], 
	 GENLC.RANGauss_NOISE_TEMPLATE[0][ifilt_obs],
	 GENLC.RANGauss_NOISE_ZP[ep] ) ;

  printf(" xxx SIG_pe(SRC,TSRC, SKY,TSKY) = %.2f,%.2f   %.2f,%.2f \n",
	 sqrt(FLUXNOISE->SQSIG_SRC), sqrt(FLUXNOISE->SQSIG_TSRC),
	 sqrt(FLUXNOISE->SQSIG_SKY), sqrt(FLUXNOISE->SQSIG_TSKY) );

  printf(" xxx SIG_pe(ZP,HOST) = %.2f, %.2f    SQSIG_RAN=%.2f\n",
	 sqrt(FLUXNOISE->SQSIG_ZP), sqrt(FLUXNOISE->SQSIG_HOST_PHOT),
	 FLUXNOISE->SQSIG_RAN );

  printf(" xxx SIG_pe(CALC, FINAL_TRUE, FINAL_DATA) = %.2f, %.2f, %.2f \n",
	 sqrt(FLUXNOISE->SQSIG_CALC_DATA),
	 fluxerr_data, fluxerr_true );

  printf(" xxx\n");  fflush(stdout);


  return ;

} // end dumpEpoch_fluxNoise_apply

// *************************************
void set_GENFLUX_FLAGS(int epoch) {

  // Created Dec 27 2019
  // Called from GENFLUX_DRIVER to set flags for
  // saturation, undefined, etc ...

  int  ifilt_obs, indx;
  bool IS_ERRPOS, IS_UNDEFINED, IS_SATURATE ;
  double obsmag, genmag, fluxerr;
  char fnam[] = "set_GENFLUX_FLAGS" ;

  // ---------- BEGIN -------------

  ifilt_obs = GENLC.IFILT_OBS[epoch] ;
  obsmag       = GENLC.mag[epoch];
  genmag       = GENLC.genmag_obs[epoch] ;
  fluxerr      = GENLC.fluxerr_data[epoch];

  IS_ERRPOS    = (fluxerr > 0 );
  IS_UNDEFINED = (genmag == MAG_UNDEFINED) ; // model undefined
  IS_SATURATE  = (obsmag == MAG_SATURATE ) ;    


  if ( IS_UNDEFINED ) 
    { GENLC.NOBS_UNDEFINED++ ; } // model is undefined 
  
  if ( IS_ERRPOS ) {
    if ( !IS_UNDEFINED ) { 
      GENLC.OBSFLAG_WRITE[epoch] = true ; 
      GENLC.NOBS++ ;
      GENLC.NOBS_FILTER[ifilt_obs]++ ;
    }
    
    if ( IS_SATURATE ) 
      { indx = INDEX_SATURATE ;  }
    else 
      { indx = INDEX_NOTSATURATE ;  }

    GENLC.NOBS_SATURATE[indx]++ ; 
    GENLC.NOBS_SATURATE_FILTER[indx][ifilt_obs]++ ; 

  } // end IS_ERRPOS
  
  return ;

} // end set_GENFLUX_FLAGS

// **************************************
void compute_lightCurveWidths(void) {

  // call gen_lightCurveWidth for each band,
  // and store GENLC.WIDTH[ifilt_obs].
  // These quantities are for SIMGEN_DUMP file
  // (varNames WIDTH_[band] ).

  int ifilt, ifilt_obs, N, ep, NEP[MXFILTINDX], ERRFLAG ;
  int MEM = (GENLC.NEPOCH+1) * sizeof(double) ;
  int OPTMASK = 1 ;
  double *TLIST[MXFILTINDX], *MAGLIST[MXFILTINDX] ;
  char FUNCALL[100], cfilt[2] ;
  char fnam[] = "compute_lightCurveWidths" ;
  
  // --------------- BEGIN -------------

  if ( GENLC.NWIDTH_SIMGEN_DUMP == 0 ) { return ; }

  // first allocate memory to store LC separately in each band.
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    GENLC.WIDTH[ifilt_obs] = -9.0 ; // init to undefined
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }
    TLIST[ifilt_obs]   = (double*)malloc(MEM);
    MAGLIST[ifilt_obs] = (double*)malloc(MEM);
    NEP[ifilt_obs]     = 0 ;
  }

  // load LC for each band

  for(ep=1; ep <= GENLC.NEPOCH; ep++ ) {
    ifilt_obs = GENLC.IFILT_OBS[ep] ;    
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }
    N = NEP[ifilt_obs];
    TLIST[ifilt_obs][N]   = GENLC.epoch_obs[ep] ;
    MAGLIST[ifilt_obs][N] = GENLC.genmag_obs[ep] ;
    NEP[ifilt_obs]++ ;
  }
 

  // compute Width, then free memory
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

    // construct error-message string
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
    sprintf(FUNCALL, "%s(%s)", fnam, cfilt);
    
    //      printf(" xxx ------------------ \n");
    //      printf(" xxx PEAKMJD = %.3f  LIBID=%d \n",
    //	     GENLC.PEAKMJD, GENLC.SIMLIB_ID );
  GENLC.WIDTH[ifilt_obs] = 
    get_lightCurveWidth(OPTMASK, NEP[ifilt_obs], 
			TLIST[ifilt_obs], MAGLIST[ifilt_obs], &ERRFLAG,FUNCALL);
    
    // free memory for this band
    free(TLIST[ifilt_obs]);  free(MAGLIST[ifilt_obs]);
  }

  return ;

} // end compute_lightCurveWidths


// *********************************************
void genmodel(
	      int ifilt_obs  // observer filter index 
	      ,int inear      // 1=>nearest rest-filter; 2=2nd nearest
	      ) {

  /*********
   Generate magnitude at each epoch according to model.
   Output is GENLC.genmag_obs[ifilt][epoch] for obs-frame model,
   or GENLC.genmag_rest[ifilt][epoch] for rest frame model.
   Note that for rest-frame models, genmag_boost transforms
   rest-frame mags to observer frame.

  ***********/

  int istat, ifilt_tmp, ifilt_rest, iep ;
  int NEPFILT, NGRID, NEPFILT_SAVE ;
  int index, isp, OPTMASK, NPAR; 

  double 
    z, zz, mu, stretch[2], delta, dm15, av, RV, AV, tmpdif
    ,errtmp, mwebv, Tep, Tshift, dummy
    ,Tmodel[MXEPSIM], TGRID[MXEPSIM]    
    ,*ptr_genmag, *ptr_generr, *ptr_epoch
    ;

  char  model[40], cfilt_obs[2], cfilt_rest[2]    ;
  char fnam[] = "genmodel" ;

  // -------- BEGIN ---------
 
  // create temp structure with epoch-list for this filter
  NEPFILT = NEPFILT_GENLC(1,ifilt_obs);

  // make local copy of model
  sprintf(model, "%s", INPUTS.MODELNAME );
  sprintf(cfilt_obs,  "%c", FILTERSTRING[ifilt_obs] );

  if ( GENFRAME_OPT == GENFRAME_REST ) {
    z      = ZAT10PC ; 
    av     = 0.0; 
    mu     = 0.0;
    mwebv  = -9.0; // not used for mlcs
  }
  else {
    z      = GENLC.REDSHIFT_HELIO ;
    av     = GENLC.AV ;
    mwebv  = GENLC.MWEBV_SMEAR ; // smeared value, not sigma
    mu     = GENLC.DLMU ;
  }


  zz = 1. + z ;
  stretch[0]   = GENLC.STRETCH ;
  stretch[1]   = GENLC.STRETCH ;
  delta        = GENLC.DELTA;  // for MLCS model
  dm15         = GENLC.DM15;
  RV           = GENLC.RV ;
  AV           = GENLC.AV ;
  ifilt_rest = -9;
  ptr_genmag = ptr_generr = &dummy; // avoid compile warning

  if ( GENFRAME_OPT == GENFRAME_REST ) {

    ptr_epoch    = &GENFILT.Trest[ifilt_obs][1] ;

    if ( inear == 1 ) {
      ptr_genmag   = &GENFILT.genmag_rest[ifilt_obs][1];
      ptr_generr   = &GENFILT.generr_rest[ifilt_obs][1];
      ifilt_rest   = GENLC.IFILTMAP_REST1[ifilt_obs];
    }
    else if ( inear == 2 ) {
      ptr_genmag   = &GENFILT.genmag_rest2[ifilt_obs][1];
      ptr_generr   = &GENFILT.generr_rest2[ifilt_obs][1];
      ifilt_rest   = GENLC.IFILTMAP_REST2[ifilt_obs];
    }
    else if ( inear == 3 ) {

      // to save CPU, only do the 3rd filter if it is nearly as 
      // close to filt1 as the 2nd closest filter.

      tmpdif  = GENLC.LAMDIF_REST3[ifilt_obs] - GENLC.LAMDIF_REST2[ifilt_obs];
      if ( tmpdif > 300.0 ) return ;  // beware hard-wired cut in Angstroms

      if ( tmpdif < 0.0 ) {
	sprintf(c1err,"GENLC.LAMDIF_REST(2,3)[ifilt_obs=%d] = %6.0f,%6.0f"
		,ifilt_obs
		,GENLC.LAMDIF_REST2[ifilt_obs]
		,GENLC.LAMDIF_REST3[ifilt_obs] );

	sprintf(c2err,"But LAMDIF_REST3 should be larger than LAMDIF_REST2");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      ptr_genmag   = &GENFILT.genmag_rest3[ifilt_obs][1];
      ptr_generr   = &GENFILT.generr_rest3[ifilt_obs][1];
      ifilt_rest   = GENLC.IFILTMAP_REST3[ifilt_obs];
    }

    sprintf(cfilt_rest, "%c", FILTERSTRING[ifilt_rest] );

  }
  else {
    // observer-frame model
    ptr_genmag   = &GENFILT.genmag_obs[ifilt_obs][1] ;
    ptr_epoch    = &GENFILT.Tobs[ifilt_obs][1] ;
    ptr_generr   = &GENFILT.generr_obs[ifilt_obs][1];    
  } 


  // - - - - - - 
  // Aug 20, 2008
  // convert ptr_epoch to account for perturbing the
  // light curve shape such as changing the rise time.
  // Note that shift is positive for positive risetime shifts,
  // and negative for positive falltime shifts.
  for ( iep=0; iep <= NEPFILT ; iep++ ) {
    Tep    = ptr_epoch[iep] ;
    Tshift = genmodel_Tshift(Tep,z); // generate epoch shift    
    Tmodel[iep] = Tep + Tshift ;
  }
  // change epoch-pointer to local array
  ptr_epoch = Tmodel ;


  // check option to interpolate model on an Epoch grid.
  // e.g., to mimic interpolated fakes on images.
  if ( INPUTS.TGRIDSTEP_MODEL_INTERP > 0.001 ) {
    NEPFILT_SAVE  = NEPFILT ;
    NGRID     = setEpochGrid( Tmodel[0], Tmodel[NEPFILT-2], TGRID );
    NEPFILT   = NGRID ;
    ptr_epoch = TGRID ;  // TGRID returned from setEpochGrid
  }
  else
    { NGRID = NEPFILT_SAVE = 0 ; }


  // -------------------------------------------
  // -------------------------------------------
  // -------------------------------------------

  if (  INDEX_GENMODEL  == MODEL_STRETCH ) {

    // this model can generate rest or observer frame

    istat = 
      genmag_stretch2 
	(
	 stretch[0]              // (I) rise-time stretch
	 ,stretch[1]             // (I) fall-time stretch
	 ,ifilt_rest          // (I) absolute rest-filter index
	 ,NEPFILT             // (I) number of epochs passed
	 ,ptr_epoch           // (I) Tobs-Tpeak (days)
	 ,ptr_genmag          // (O) perfect magnitudes
	 ,ptr_generr          // (O) ideal rest mag-errs
	 );


    // cap errors at .25 mag to avoid ridiculous light curves
    // from blown-up U,R,I errors (blown up for fitter)

    for ( iep=0; iep<NEPFILT; iep++ ) {
      errtmp = ptr_generr[iep];
      if ( errtmp > .25 ) { errtmp = .25 ; }
      ptr_generr[iep] = errtmp ;
    }

  }

  else if ( INDEX_GENMODEL  == MODEL_FIXMAG ) {
    GENLC.FIXMAG = GENLC.NOSHAPE ;
    for ( iep = 0 ; iep < NEPFILT; iep++ ) {
      ptr_genmag[iep] = GENLC.FIXMAG ;
      ptr_generr[iep] = 0.0 ;
    }
    sprintf(GENLC.SNTYPE_NAME, "%s", INPUTS.MODELNAME ); 
    GENLC.TEMPLATE_INDEX = MODEL_FIXMAG ;  //anything but zero
  }
  else if ( INDEX_GENMODEL == MODEL_SIMLIB ) {

    // gen mags are already loaded since ptr_genmag points to 
    // GENFILT.genmag_obs

    /*
    printf(" xxx %s: %d \n", SIMLIB_HEADER.NOBS_SIM_MAGOBS);

    if ( SIMLIB_HEADER.NOBS_SIM_MAGOBS == 0 ) {
      sprintf(c1err,"All SIMLIB SIM_MAGOBS=%.1f --> all true flux = 0", 
	      MAG_ZEROFLUX);
      sprintf(c2err,"Probably have wrong SIMLIB");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    */

  }
  else if ( INDEX_GENMODEL  == MODEL_MLCS2k2 ) {

    // define function to translage ifilt_rest into ifilt_model(ifilt_rest]
    ifilt_tmp  = GENLC.IFILTINVMAP_REST[ifilt_rest] ;

    genmag_mlcs2k2( 
		   ifilt_tmp          // (I) 0-8 => UBVRIYJHK
		   ,delta             // (I) shape parameter 
		   ,NEPFILT           // (I) number of epochs
		   ,ptr_epoch         // (I) Trest-Tpeak (days)
		   ,ptr_genmag        // (O) ideal rest mags
		   ,ptr_generr        // (O) ideal rest mag-errs
		   ) ;

  }

  else if ( INDEX_GENMODEL  == MODEL_SALT2 ) {

    // 1: bit0 => return flux instead of mag;
    // 2: bit1 => print warning when flux < 0
    // 8: bit3 => SALT2 dump

    OPTMASK      = 0 ; // return mag
                      
    /*
    printf(" xxx passing ifilt_obs=%d x0=%f x1=%f c=%f \n",
	   ifilt_obs,GENLC.SALT2x0, GENLC.SALT2x1, GENLC.SALT2c );
    */
    // apply scatter matrix
    double S2x0, S2x1, S2c, S2mB, tmp, logMass=-9.0 ;
    S2mB = GENLC.SALT2mB + GENLC.COVMAT_SCATTER[0] ;
    S2x1 = GENLC.SALT2x1 + GENLC.COVMAT_SCATTER[1] ;
    S2c  = GENLC.SALT2c  + GENLC.COVMAT_SCATTER[2] ;
    tmp  = -0.4 * GENLC.COVMAT_SCATTER[0] ;
    S2x0 = GENLC.SALT2x0 * pow(10.0,tmp);
    int m = SNHOSTGAL.IMATCH_TRUE;
    if ( m >= 0 ) { logMass = SNHOSTGAL_DDLR_SORT[m].LOGMASS_TRUE; }

    double parList_SN[4]   = { S2x0, S2x1, S2c, S2x1 } ;
    double parList_HOST[3] = { RV, AV, logMass } ;
 
    ptr_generr[0] = 4.44; // xxx REMOVE

    genmag_SALT2 (
		  OPTMASK         // (I) bit-mask options
		  ,ifilt_obs      // (I) obs filter index 
		  ,parList_SN     // (I) SN params: x0, x1, x1forERR, c
		  ,parList_HOST   // (I) host params: RV, AV, logMass
		  ,mwebv          // (I) Galactic E(B-V)
		  ,z,z            // (I) redshift, and z used for error
		  ,NEPFILT        // (I) number of epochs
		  ,ptr_epoch         // (I) obs-frame time (days)
		  ,ptr_genmag        // (O) mag vs. Tobs
		  ,ptr_generr        // (O) mag-errs
		  ) ;    

    /* xxxxxxxxxx
    if ( ifilt_obs == INPUTS.DEBUG_FLAG ) {
      double Trest_tmp = ptr_epoch[0]/(1.0+z);
      printf(" xxx %s: ifilt_obs=%d  z=%.5f  Trest=%.3f magerr=%.3f \n",
	     fnam, ifilt_obs, z, Trest_tmp, ptr_generr[0] ) ;
    }
    xxxxxxx */

  }

  else if ( INDEX_GENMODEL  == MODEL_SIMSED ) {

    // 1: bit0=> return flux instead of mag;
    // 2: bit1 => print warning when flux < 0

    NPAR         = INPUTS.NPAR_SIMSED;

    OPTMASK      =  0              // 0->return mags   
      + INPUTS.OPTMASK_SIMSED      // 4 => sequential GRIDONLY option
      ;

    //    if ( GENLC.CID==19201 && ifilt_obs==1 ) { OPTMASK += 8 ; }
    
    genmag_SIMSED (
		  OPTMASK            // (I) bit-mask options
		  ,ifilt_obs         // (I) obs filter index 
		  ,GENLC.SALT2x0     // (I) x0 (flux-scale)
		  ,NPAR              // (I) Number of SIMSED pars + baggage
		  ,INPUTS.GENFLAG_SIMSED // (I) SIMSED flag per param
		  ,GENLC.SIMSED_IPARMAP  // (I) ipar map
		  ,GENLC.SIMSED_PARVAL   // (I) SIMSED gen params (lumipar)
		  ,GENLC.RV, GENLC.AV        // (I) added Mar 1 2017
		  ,mwebv             // (I) Galactic E(B-V)
		  ,z                 // (I) redshift
		  ,NEPFILT           // (I) number of epochs
		  ,ptr_epoch         // (I) rest-frame time (days)
		  ,ptr_genmag        // (O) mag vs. Tobs
		  ,ptr_generr        // (O) ideal rest mag-errs
		  ,&GENLC.TEMPLATE_INDEX // (O) SED index
		  ) ;

    if ( INPUTS.OPTMASK_SIMSED == OPTMASK_GEN_SIMSED_GRIDONLY ) 
      { GENLC.CID = GENLC.TEMPLATE_INDEX ; }

  }

  else if ( IS_PySEDMODEL ) {

    // python-based SED model: BYOSED or SNEMO
    int NHOSTPAR; char *NAMES_HOSTPAR = NULL; 
    double VAL_HOSTPAR[MXHOSTPAR_PySEDMODEL];
    NHOSTPAR = fetch_HOSTPAR_GENMODEL(2, NAMES_HOSTPAR, VAL_HOSTPAR);

    genmag_PySEDMODEL(
		      GENLC.CID
		      ,GENLC.REDSHIFT_HELIO  // (I) heliocentric redshift 
		      ,GENLC.REDSHIFT_CMB    // (I) CMB redshift
		      ,GENLC.DLMU            // (I) distance modulus
		      ,mwebv               // (I) E(B-V) for Milky Way
		      ,NHOSTPAR            // (I) number of host params to pass
		      ,VAL_HOSTPAR         // (I) host property values
		      ,ifilt_obs           // (I) filter index
		      ,NEPFILT             // (I) number of epochs
		      ,ptr_epoch           // (I) rest-frame time (days)
		      ,ptr_genmag        // (O) mag vs. Tobs
		      ,ptr_generr        // (O) ideal rest mag-errs
		      );		  
    
  }

  else if ( INDEX_GENMODEL  == MODEL_SNOOPY ) {

    ifilt_tmp = GENLC.IFILTINVMAP_REST[ifilt_rest] ;
    genmag_snoopy (
		   ifilt_tmp             // (I) SNoopY internal filter index
		   ,stretch[0]           // (I) shape parameter 
		   ,NEPFILT              // (I) number of epochs
		   ,ptr_epoch            // (I) Trest-Tpeak (days)
		   ,ptr_genmag           // (O) ideal rest mags
		   ,ptr_generr           // (O) rest mag-errors
		   );

  }

  else if ( INDEX_GENMODEL  == MODEL_S11DM15 ) {

    genmag_S11DM15(ifilt_obs, dm15, AV, RV, mwebv, z, mu,
		   NEPFILT, ptr_epoch, ptr_genmag, ptr_generr);

  }

  else if (  INDEX_GENMODEL  == MODEL_NON1ASED ) {

    index    = GENLC.TEMPLATE_INDEX; 
    isp      = GENLC.NON1ASED.ISPARSE ;
    sprintf(GENLC.SNTYPE_NAME,  "%s", INPUTS.NON1ASED.LIST_TYPE[index] );
    sprintf(GENLC.SNTEMPLATE,   "%s", INPUTS.NON1ASED.LIST_NAME[index] );
    sprintf(GENLC.NON1ASED.TYPE[isp], "%s", GENLC.SNTYPE_NAME );

    genmag_NON1ASED (
		  index                // (I) NON1A class
		  ,ifilt_obs           // (I) obs-filter index
		  ,mwebv               // (I) Galactic E(B-V)
		  ,z                   // (I) redshift
		  ,GENLC.DLMU          // (I) distance modulus
		  ,RV, AV             // (I) RV and AV of host
		  ,NEPFILT             // (I) number of epochs
		  ,ptr_epoch           // (I) Tobs (days)
		  ,ptr_genmag          // (O) obs mags
		  ,ptr_generr          // (O) obs mag-errs
		  );
  }
  else if (  INDEX_GENMODEL  == MODEL_NON1AGRID ) {  // Mar 2016

    /*
    printf(" xxx pass NON1AGRID args: \n");
    printf(" xxx ifilt_obs=%d  mwebv=%f z=%f \n", ifilt_obs, mwebv8, z8);
    printf(" xxx RV=%f AV=%f  NEP=%d \n", RV8, AV8, NEPFILT);
    fflush(stdout);
    */
    
    // use same randoms for each band
    double ranWgt   = GENLC.NON1AGRID_RANFlat ;
    double ranSmear = GENLC.NON1AGRID_RANGauss ;
    
    genmag_NON1AGRID (
		   ifilt_obs           // (I) obs-filter index
		   ,mwebv               // (I) Galactic E(B-V)
		   ,GENLC.REDSHIFT_CMB  // (I) CMB redshift, NOT zhelio
		   ,RV, AV             // (I) RV and AV of host
		   ,ranWgt, ranSmear    // (I) to pick random index & magSmear
		   ,NEPFILT             // (I) number of epochs
		   ,ptr_epoch           // (I) Trest-Tpeak (days)
		   ,ptr_genmag          // (O) obs mags
		   ,ptr_generr          // (O) obs mag-errs
		   ,&GENLC.MAGSMEAR_COH[0] // (O) store random magSmear
		      );

    GENLC.TEMPLATE_INDEX = fetchInfo_NON1AGRID("NON1A_INDEX");
    GENLC.SNTYPE         = fetchInfo_NON1AGRID("NON1A_ITYPE_USER");
  }
  else if ( INDEX_GENMODEL  == MODEL_LCLIB ) {  // July 2017

    set_TobsRange_LCLIB(GENLC.epoch_obs_range);

    sprintf(GENLC.SNTYPE_NAME, "%s", LCLIB_INFO.NAME_MODEL ); 
    double *ptr_template = &GENLC.genmag_obs_template[ifilt_obs];
    double LAMAVG        = INPUTS.LAMAVG_OBS[ifilt_obs];
    double TobsPeak ;

    genmag_LCLIB(GENLC.CID, ifilt_obs, &GENLC.RA, &GENLC.DEC, 
		 &mwebv, LAMAVG,
		 NEPFILT, ptr_epoch, 
		 ptr_genmag, ptr_template,    // return args
		 &TobsPeak                    // Tobs with peak brightness
		 );

    // for non-recurring events, set PKMJD like any other transient
    if ( LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_NONRECUR ) {
      GENLC.PEAKMJD  = INPUTS.GENRANGE_PEAKMJD[0] + TobsPeak ;
    }

    /*
    printf(" xxx genmag_LCLIB returns zphot = %.3f +- %.3f  (ztrue=%.3f)\n",
	   SNHOSTGAL.ZPHOT, SNHOSTGAL.ZPHOT_ERR, GENLC.REDSHIFT_CMB );
    */

  }

  // !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=
  //         end of model if-block
  // !=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=

  // ---------------------------------------
  // check option to interpolate model mags (see setEpochGrid above)
  if ( INPUTS.TGRIDSTEP_MODEL_INTERP > 0.001 ) {
    NEPFILT = NEPFILT_SAVE ;
    interpEpochGrid(NEPFILT, Tmodel, NGRID,              // I
		    ptr_epoch, ptr_genmag, ptr_generr) ; // I/O
  } 

  // ======================================
  // store peak mag in rest-frame

  if ( GENFRAME_OPT == GENFRAME_REST ) {
    iep =  NEPFILT - 1 ; // last epoch is always peak
    GENLC.peakmag_rest[ifilt_rest] = ptr_genmag[iep] ;
  }

  // apply intrinsic model mag-smearing AFTER model-mag generation
  // Note that ptr_genmag is both an input and output arg.
  genmodelSmear(NEPFILT, ifilt_obs, ifilt_rest, z, 
		 ptr_epoch, ptr_genmag, ptr_generr ); 

  NEPFILT = NEPFILT_GENLC(-1, ifilt_obs);

  return ;

}  // end of genmodel


// ********************************************
double  genmodel_Tshift(double T, double z) {

  // Created Jun 2017
  // Return model phase shift to smear rise/fall times.
  // Input T = phase (days) w.r.t. peak
  // Works only for rest-frame models.

  double Tshift; 
  double zz = 1.0 + z ;
  // --------------- BEGIN --------------

  if ( GENFRAME_OPT == GENFRAME_REST ) {
    if ( T < 0.0 ) {
      Tshift = -(T/18.0) * GENLC.RISETIME_SHIFT * zz ;
    }
    else {
      Tshift = -(T/15.0) * GENLC.FALLTIME_SHIFT * zz ;
    }    
  }
  else {
    // no shift for OBS-frame models (e.g., SALT2, NON1ASED)
    Tshift = 0.0 ;
  }

  return(Tshift) ;

} // end genmodel_Tshift

// ********************************************
int setEpochGrid(double TMIN, double TMAX, double *TGRID ) {

  // Created Oct 23 2015
  // Inputs: TMIN and TMAX
  // Ouput: TGRID[i] = fixed grid
  // Function outut = size of grid

  double t, TSTEP, TSTART, TEND;
  int  NGRID, iTround ;
  //  char fnam[] = "setEpochGrid" ;

  // ---------------- BEGIN ----------------

  TSTEP = INPUTS.TGRIDSTEP_MODEL_INTERP ;
  NGRID = 0 ;
  if ( TSTEP <= 0.0 ) { return(NGRID) ; }

  iTround = (int)TMIN ;
  TSTART  = (double)iTround -2.0 ;
  if ( TSTART > 0.0 ) { TSTART = -TSTEP; }

  TEND    = TMAX + TSTEP ;
  if ( TEND < TSTEP ) { TEND = TSTEP + 1.0 ; }

  /*
  if( TEND < 20.0 ) {
    printf(" xxx %s: TMIN,TMAX -> TSTART,TEND = %.2f, %.2f -> %.2f, %.2f \n", 
	   fnam, TMIN, TMAX, TSTART, TEND );
  }
  */

  for(t=TSTART; t <= TEND; t += TSTEP ) {
    TGRID[NGRID] = t;
    NGRID++ ;
  }

  return(NGRID);

} // end of setEpochGrid


void interpEpochGrid(
		     int NEP_LC,           // (I) Nep for LC
		     double *TList_LC,     // (I) Tlist for LC
		     int NGRID,            // (I) Number of grid points
		     double *TList,        // (I/O) T at each grid point
		     double *magList,      // (I/O) mag at each grid point
		     double *magerrList )  // (I/O) magerr at each grid
{ 

  // Oct 23 2015
  // Use input grid of mags (maglist vs. Tlist) to interpolate
  // mags on the light curve epochs (TList_LC).
  // Note that TList, magList and magerrList are input
  // with GRID values and overwritten with the LightCurve values.
  //
  // This function is associated with user-input key 
  // TGRIDSTEP_MODEL_INTERP  to generate LC mags on a grid
  // and then interpolate.

  double *GRID_T, *GRID_MAG, *GRID_MAGERR;
  double T, mag, magerr;
  int MEMD, ep;
  int OPT_INTERP=1;  // linear interp
  char fnam[] = "interpEpochGrid" ;

  // -------------- BEGIN -------------------
  
  // create arrays to store input GRID
  MEMD        = NGRID * sizeof(double);
  GRID_T      = (double*) malloc( MEMD );
  GRID_MAG    = (double*) malloc( MEMD );
  GRID_MAGERR = (double*) malloc( MEMD );

  for(ep=0; ep<NGRID; ep++ ) {
    GRID_T[ep]      = TList[ep];
    GRID_MAG[ep]    = magList[ep];
    GRID_MAGERR[ep] = magerrList[ep];    

    TList[ep]      = 999. ; // init output
    magList[ep]    = 999. ; // idem 
    magerrList[ep] = 999. ; // idem
  }

  // interpolate for each LC epoch and overwrite input arrays
  for(ep=0; ep < NEP_LC; ep++ ) {    

    T      = TList_LC[ep];

    mag    = interp_1DFUN(OPT_INTERP, T, 
			  NGRID, GRID_T, GRID_MAG, fnam);

    magerr = interp_1DFUN(OPT_INTERP, T,
			  NGRID, GRID_T, GRID_MAGERR, fnam);

    TList[ep]      = T ;
    magList[ep]    = mag;
    magerrList[ep] = magerr ;

    /*
    printf(" xxx ep=%2d  T=%.2f  mag= %.2f +- %.2f \n", 
	   ep, T, mag, magerr ); // xxxxxxxxxxxxxxxxxxxxxx
    */

  }

  free(GRID_T);
  free(GRID_MAG);
  free(GRID_MAGERR);

} // end of interpEpochGrid

// ********************************************
void genmodelSmear(int NEPFILT, int ifilt_obs, int ifilt_rest,  double z, 
		   double *ptr_epoch, double *ptr_genmag, double *ptr_generr ) {

  //
  // Created Jul 1, 2010 by R.Kessler
  // Move model-mag smearing  code from end of genmag_smear to here.
  // Also add new option to specify smearing vs. lambda
  //
  // Inputs:
  //
  // NEPFILT     : Number of epochs for this filter
  // ifilt_obs   : observer filter 
  // ifilt_rest  : rest-frame filter index (for rest-frame model only)
  // z           : redshift
  //
  // *ptr_epoch  : pointer to epochs (rest or observer)
  //
  // *ptr_genmag : pointer to list of generated mags
  //               NOTE THAT *ptr_genmag values are CHANGED
  //               based on smearing model !!!
  //
  // Jul 19, 2010: return if DO_MODELSMEAR == 0.
  // Oct 16, 2010: return for GRID option
  //
  // Apr 2012: cleanup with call to genSmear_ERRSCALE 
  //           Also call new functino get_genSmear(...)
  //
  // Jun 14 2016: load GENLC.MAGSMEAR_COH at end of function
  //

  int
    iep, opt_frame, ifilt_local
    ,USE_GENMODEL_ERRSCALE, USE_GENSMEAR_MODEL       
    ,ONE = 1
    ;

  double 
    ran_COH, ran_FILT
    ,magSmear=0.0, magSmear_model, magSmear_tmp
    ,smearsig=0.0, smearsig_fix, smearsig_model
    ,Tep, Tpeak, Trest, lamrest, Z1, parList[10]
    ;

  float lamavg4, lamrms4, lammin4, lammax4 ;
  bool  USE_SMEARSIG = false;
  char cfilt[2];
  char fnam[] = "genmodelSmear" ;

  // -------------- BEGIN ------------

  if ( !INPUTS.DO_MODELSMEAR  ) 
    { return ; }

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) 
    { return ; }

  if ( GENFRAME_OPT == GENFRAME_REST  ) 
    { ifilt_local = ifilt_rest ;  Z1 = 1.0 ; }
  else 
    { ifilt_local = ifilt_obs ;   Z1 = 1.0 + z ; }


  sprintf(cfilt,  "%c", FILTERSTRING[ifilt_local] ); // for err msg only

  // sanity check that epoch at peak really as Tpeak = 0
  Tpeak = ptr_epoch[NEPFILT - 1];

  if ( fabs(Tpeak) > 0.0001 && NEPFILT > 0 ) {
    sprintf(c1err,"Tpeak(%s) = %f but expected T=0 at peak", 
	    cfilt, Tpeak );
    sprintf(c2err,"PEAKMJD(GEN) = %f  LIBID=%d  NEPFILT=%d \n", 
	    GENLC.PEAKMJD, GENLC.SIMLIB_ID, NEPFILT );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // init all smearing to zero
  for ( iep = 0; iep <= NEPFILT; iep++ )  {  
    GENFILT.genmag_smear[ifilt_obs][iep] = 0.0 ; 
  }


  // ----------------------------------------------
  // start with coherent smearing (GENMAG_SMEAR)
  // that is common to all passbands.

  if ( INPUTS.GENMAG_SMEAR[0] > 0.0001 ) {

    ran_COH = GENLC.GENSMEAR_RANGauss_FILTER[0] ; // retrieve random #

    magSmear = INPUTS.GENMAG_SMEAR[0] * ran_COH ;
    GENLC.MAGSMEAR_COH[0] = magSmear ;  // store global

    for ( iep = 1; iep <= NEPFILT; iep++ )   { 
      ptr_genmag[iep-1]                    += magSmear ;
      GENFILT.genmag_smear[ifilt_obs][iep] += magSmear ;
    }
  }


  // check old model of using the fit-model errors
  if ( INPUTS.GENMODEL_ERRSCALE > 1.0E-6  ) 
    { USE_GENMODEL_ERRSCALE = 1 ; }
  else
    { USE_GENMODEL_ERRSCALE = 0 ; }


  USE_GENSMEAR_MODEL = 
    ( (istat_genSmear() > 0) && (IFLAG_GENSMEAR == IFLAG_GENSMEAR_FILT) ) ;

  smearsig_fix   = 0.0 ;
  smearsig_model = 0.0 ;

  // first get random number for this passband based on 
  // rest- or observer-frame model

  ran_FILT = GENLC.GENSMEAR_RANGauss_FILTER[ifilt_local] ; 

  // now get filter-dependent smearing
  if ( INPUTS.NFILT_SMEAR > 0 )  { 
    smearsig_fix = INPUTS.GENMAG_SMEAR_FILTER[ifilt_local]; 
    USE_SMEARSIG = true ;
  }

  // loop over epochs and apply same mag-smear 

  for ( iep=1; iep <= NEPFILT; iep++ ) {

    magSmear = smearsig_model = 0.0 ;
    Tep   = ptr_epoch[iep - 1] ; // rest or obs.
    Trest = Tep / Z1 ;          // Z1=1 (rest) or 1+z (obs)
     
    if ( USE_GENMODEL_ERRSCALE  ) { // practially obsolete
      smearsig_model = genSmear_ERRSCALE(ptr_generr,iep, NEPFILT ); 
      USE_SMEARSIG = true ;
    }

    // add two sources of smear in quadrature: FILTER and ERRSCALE
    if ( USE_SMEARSIG ) {
      smearsig = sqrt( pow(smearsig_model,2.) + pow(smearsig_fix,2.) );      
      magSmear = smearsig * ran_FILT ;  // apply filter-dependent random smear
    }

    // model-dependent smearing; note that get_genSmear returns
    // a randomly generated magSmear rather than a sigma-smear.
    // Also, this magSmear over-writes previous magSmear since 
    // the two cannot both be set.
    if ( USE_GENSMEAR_MODEL ) {
      // convert mean lambda of filter in obs frame to rest frame
      opt_frame = GENFRAME_OBS ;
      get_filtlam__(&opt_frame, &ifilt_obs, 
		    &lamavg4, &lamrms4, &lammin4, &lammax4 );
      lamrest = (double)lamavg4 / ( 1.0 + z );   
      
      parList[0] = Trest;
      parList[1] = GENLC.SALT2x1;
      parList[2] = GENLC.SALT2c ;
      parList[3] = SNHOSTGAL_DDLR_SORT[0].LOGMASS_TRUE ;
      get_genSmear(parList, ONE,  &lamrest, &magSmear_model);
      magSmear = magSmear_model ;
    }
   

    //* xxxxxxxxxxxxxxxx
    if  ( fabs(Tep) < .01 && ifilt_obs == -2 ) {
      printf("\t xxx ------------------------------------- \n");
      printf("\t xxx ifllt[obs,rest] = %d %d   T=%5.2f  mag=%6.3f \n", 
	     ifilt_obs, ifilt_rest, Tep, *(ptr_genmag+iep-1) );
      printf("\t xxx magSmear=%7.4f  sig[model,fix]=%6.3f,%6.3f  rr=%6.3f \n",
		magSmear, smearsig_model, smearsig_fix, ran_FILT );
    }
    // xxxxxxxxxxxxxx */


    // insanity check:
    if ( fabs(magSmear) > 30.0 ) {
      printf("\n xxx ---------------------------------------- \n");
      printf(" xxx smearsig(model,fix,tot) = %.3f,%.3f,%.3f  magSmear=%f\n", 
	     smearsig_model, smearsig_fix, smearsig, magSmear );
      printf(" xxx lamrest = %.1f   Trest=%.2f\n", lamrest, Trest);
      printf(" xxx USE_GENSMEAR_MODEL = %d \n", USE_GENSMEAR_MODEL );
      sprintf(c1err,"magSmear = %f is for filt=%s T=%6.1f", 
	      magSmear, cfilt, Tep );
      sprintf(c2err,"magsig(model)=%4.1f  GaussRan=%5.3f  errscale=%4.2f", 
	      smearsig, ran_FILT, INPUTS.GENMODEL_ERRSCALE );
      errmsg(SEV_WARN, 0, fnam, c1err, c2err );
    }

    // Aug 15, 2009: clip mag-smear at +- 3 mag to avoid catastrophe
    if ( magSmear < -3.0 ) { magSmear = -3.0 ; }
    if ( magSmear > +3.0 ) { magSmear = +3.0 ; }

    // check flag (negative smearsig) to interpolate smearing
    if ( smearsig_fix < -1.0E-6 ) {
      magSmear = genmodelSmear_interp(ifilt_obs,iep) ;
    }

    // add small phase-dependent scatter (Feb 11 2020)
    get_genSmear_phaseCor(GENLC.CID, Trest, &magSmear_tmp);
    magSmear += magSmear_tmp ; 

    if ( iep == -3 ) {  
      printf(" 66666  %d/%d  %2d  %7.2f  %8.4f \n", 
	     GENLC.CID, GENLC.SIMLIB_ID, ifilt_obs, Trest, magSmear); 
    }

    // store magSmear in global
    ptr_genmag[iep-1]                    += magSmear ;
    GENFILT.genmag_smear[ifilt_obs][iep] += magSmear ;

  } // ep loop

  if ( istat_genSmear() > 0 ) {
    GENLC.MAGSMEAR_COH[0] = GENSMEAR.MAGSMEAR_COH[0];
    GENLC.MAGSMEAR_COH[1] = GENSMEAR.MAGSMEAR_COH[1]; // optional 2nd term
  } 

  return ;

} // end of genmodelSmear


// ********************************
double genSmear_ERRSCALE(double *ptr_generr, int iep, int NEPFILT ) {

  // April 2012 [code moved from genmodelSmear]
  // Implement ERRSCALE model where intrinsic smear
  // is based on the genmag model error.
  // This model is really obsolete and is not recommended.
  //
  // ptr_generr is an array of fitter-model errors.
  // 'iep' is the current epoch index.

  double errtmp, smearsig_model, GENMODEL_ERRSCALE  ;
  int iep_modelerr;
  char fnam[] = "genSmear_ERRSCALE" ;

  // ----------- BEGIN ----------

  iep_modelerr = -9 ;

  if ( INPUTS.GENMODEL_ERRSCALE_OPT == 1 )       
    { iep_modelerr = NEPFILT ; }  // use peak model error at all epochs
  else if ( INPUTS.GENMODEL_ERRSCALE_OPT == 2 )
    { iep_modelerr = iep ; } // use epoch-dependent model error    
  else {
    sprintf(c1err,"Invalid  GENMODEL_ERRSCALE_OPT = %d",
	    INPUTS.GENMODEL_ERRSCALE_OPT  );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }

  // make a few sanity checks

  if ( iep_modelerr <= 0 || iep_modelerr > NEPFILT ) {
    sprintf(c1err,"Insane iep_modelerr = %d", iep_modelerr);
    sprintf(c2err,"at iep=%d ", iep );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  GENMODEL_ERRSCALE  = (double)INPUTS.GENMODEL_ERRSCALE; 
  errtmp             = *(ptr_generr + iep_modelerr - 1 ); // model error
  smearsig_model     = errtmp * GENMODEL_ERRSCALE ;

  return smearsig_model ; 

} // end of genSmear_ERRSCALE
 
// ****************************************************
double  genmodelSmear_interp(int ifilt_interp, int iep) {

  // Created Feb 3, 2012 by R.Kessler
  // interpolate smearing for this observer-frame ifilt_interp.
  // 'iep' is the epoch index.
  // Interpolate with filters that have closest wavelength.
  //
  // Initial use is for SNACC teeth since the nominal
  // filter mags must already be set.

  int ifilt, ifilt_obs,  ifilt1_obs, ifilt2_obs ;

  double 
    smear 
    ,LAM_INTERP
    ,LAM, LAM1, LAM2, LAMDIF, LAMDIF_MIN
    ,SMEAR1, SMEAR2
    ;

  char fnam[] = "genmodelSmear_interp" ;

  // -------------- BEGIN --------------

  LAM_INTERP = (double)INPUTS.LAMAVG_OBS[ifilt_interp] ;


  ifilt1_obs = ifilt2_obs = -9;

  // get closest filter
  LAMDIF_MIN = 99999999. ;
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];

    // only include filters with positive smear
    if ( INPUTS.GENMAG_SMEAR_FILTER[ifilt_obs] < 0.000001 ) { continue ; }

    LAM       = INPUTS.LAMAVG_OBS[ifilt_obs] ;
    LAMDIF    = fabs(LAM - LAM_INTERP) ;
    if ( LAMDIF < LAMDIF_MIN ) { 
      LAMDIF_MIN = LAMDIF ; 
      ifilt1_obs = ifilt_obs ; 
    }
  }     // ifilt


  // get 2nd closest filter
  LAMDIF_MIN = 99999999. ;
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    if ( ifilt_obs == ifilt1_obs ) { continue ; }

    // only include filters with positive smear
    if ( INPUTS.GENMAG_SMEAR_FILTER[ifilt_obs] < .0000001 ) { continue ; }

    LAM       = INPUTS.LAMAVG_OBS[ifilt_obs] ;
    LAMDIF    = fabs(LAM - LAM_INTERP) ;
    if ( LAMDIF < LAMDIF_MIN ) { 
      LAMDIF_MIN = LAMDIF ; 
      ifilt2_obs = ifilt_obs ; 
    }
  }     // ifilt


  // idiot check.
  if ( ifilt2_obs < 1 ) {
    sprintf(c1err,"Could not find nearest filters to %c (LAM=%6.1f)",
	    FILTERSTRING[ifilt_interp], INPUTS.LAMAVG_OBS[ifilt_interp] );
    sprintf(c2err,"ifilt1_obs=%d  ifilt2_obs=%d", 
	    ifilt1_obs, ifilt2_obs);
    errmsg( SEV_FATAL,0 , fnam, c1err, c2err);
  }

  // interpolate
  LAM1   = (double)INPUTS.LAMAVG_OBS[ifilt1_obs] ;
  LAM2   = (double)INPUTS.LAMAVG_OBS[ifilt2_obs] ;
  LAMDIF = (LAM_INTERP -  LAM1);
  SMEAR1 = GENFILT.genmag_smear[ifilt1_obs][iep] ;  
  SMEAR2 = GENFILT.genmag_smear[ifilt2_obs][iep] ;  
  
  smear = SMEAR1 + (SMEAR2-SMEAR1) * LAMDIF/(LAM2-LAM1);

  if ( iep == -3 ) {
    
    printf(" xxx -------------------------------------------------- \n");
    printf(" xxx LAMDA[%c,%c -> %c] =  %7.1f , %7.1f  -> %7.1f   \n"
	   ,FILTERSTRING[ifilt1_obs]
	   ,FILTERSTRING[ifilt2_obs]
	   ,FILTERSTRING[ifilt_interp]
	   ,LAM1, LAM2, LAM_INTERP );

    printf(" xxx SMEAR[%c,%c -> %c] = %7.4f  , %7.4f  -> %7.4f \n" 
	   ,FILTERSTRING[ifilt1_obs]
	   ,FILTERSTRING[ifilt2_obs]
	   ,FILTERSTRING[ifilt_interp]
	   ,SMEAR1, SMEAR2, smear );
  }



  return smear ;

}  // end of smearsig_interp




// ==================================================
void ZERO_COVMAT_SCATTER(void)
{
  //  Created July 2011 by R.Biswas 
  //  April 12, 2012 (RK) init the Cholesky matrix too.

  // ------------ BEGIN --------------

  /// Set INPUTS. covariance matrix to 0.  
  INPUTS.NCOVMAT_SCATTER = 0;

  int i,j ;
  for (i = 0 ; i < 3 ; i++ ){
    for (j = 0 ; j < 3 ; j++ ){
      INPUTS.COVMAT_SCATTER[i][j]          = 0. ;
      INPUTS.COVMAT_SCATTER_SQRT[i][j]     = 0. ;
      INPUTS.COVMAT_SCATTER_REDUCED[i][j]  = 0. ;
      INPUTS.COVMAT_CHOLESKY_SCATTER[i][j] = 0.0 ; // RK
      
    }
  }
}//End of ZERO_COVMAT_SCATTER

void PARSE_COVMAT_SCATTER(char *CKEY , char *CVAL) 
{
  /// Created July 2011 by R.Biswas 
  /// Parses a character array CKEY specifying the element of 
  /// the covariance matrix and a character array CVAL 
  /// having the associated value to assign the covariance matrix.

  /// CKEY     (I) char string specifying COVMAT_SCATTER element 
  ///			of the form COVMAT_SCATTER[0][0]:
  ///			or COVMAT_SCATTER_SQRT[0][0]:
  ///			or COVMAT_SCATTER_REDUCED[0][0]:
  /// CVAL     (I) char string giving value of the the covariance
  ///			matrix entry

  char fname [100];
  
  // --------------- BEGIN -----------------

  strcpy(fname,"PARSE_COVMAT_SCATTER");


  //convert the char value into a double
  double cvalcont;
  cvalcont = atof(CVAL);	


  char covmatentries[1000];

  int m , n ; 		
  //PARSE routine looks for something like 
  //COVMAT_SCATTER[m][n] only
			
  strcpy(covmatentries,"COVMAT_SCATTER[0][0]:");
  int lencovmat= strlen(covmatentries);
  lencovmat = lencovmat+1;

  for ( m = 0 ; m < 3 ; m++ ){
    for ( n = 0 ; n < 3 ; n++ ){ 
      sprintf(covmatentries,"COVMAT_SCATTER[%d][%d]:", m ,n );
      if (0==strncmp(covmatentries,CKEY,lencovmat)){
	INPUTS.COVMAT_SCATTER[m][n]= cvalcont;
	INPUTS.NCOVMAT_SCATTER++ ;
	return;
      }
    }
  }
  
  //PARSE routine looks for something like 
  //COVMAT_SCATTER_SQRT[m][n] only

  strcpy(covmatentries , "COVMAT_SCATTER_SQRT[0][0]:");
  lencovmat= strlen(covmatentries);
  lencovmat = lencovmat+1;
  
  for ( m = 0 ; m < 3 ; m++ ){
    for ( n = 0 ; n < 3 ; n++ ){ 
      sprintf(covmatentries,"COVMAT_SCATTER_SQRT[%d][%d]:", m ,n );
      if (0==strncmp(covmatentries,CKEY,lencovmat)){
	//code to abort if one enters COVMAT_SCATTER_SQRT
	//elements that are off diagonal
	if (n != m) {
	  sprintf(c1err,
		  " off-diagonal COVMAT_SCATTER_SQRT entered: %d %d", m , n );
	  errmsg( SEV_FATAL,0 , fname, c1err, c2err);
	}
	INPUTS.COVMAT_SCATTER_SQRT[m][n]= cvalcont;
	INPUTS.NCOVMAT_SCATTER++ ;
	return;
      }
    }
  }

  //Code to parse off diagonal reduced scatter 
  //terms
  for ( m = 0 ; m < 3 ; m++ ){
    for ( n = 0 ; n < 3 ; n++ ){ 
      sprintf(covmatentries,"COVMAT_SCATTER_REDUCED[%d][%d]:", m ,n );
      if (0==strncmp(covmatentries,CKEY,lencovmat)){
	//code to abort if one enters COVMAT_SCATTER_SQRT
	//elements that are off diagonal
	if (n == m) {
	  sprintf(c1err,
		  " diagonal COVMAT_SCATTER_REDUCED entered: %d %d", m , n );
	  errmsg( SEV_FATAL,0 , fname, c1err, c2err);
	}
	INPUTS.COVMAT_SCATTER_REDUCED[m][n]= cvalcont;
	INPUTS.NCOVMAT_SCATTER++ ;
	return;
      }
    }
  }
  // If nothing found return
  return;
 }//END OF PARSE_COVMAT_SCATTER

// =======================================
void INIT_COVMAT_SCATTER( void )
{
  // Created July 2011 by R.Biswas
  // todo: 
  // - use c1err and c2err instead of errormsg and c2err (done)
  // - move comments from PARSE_COVMAT to here; 
  // - fill GENLC.COVMAT_SCATTER_README with cov matrix text string
  // - print the above
  //

  char fnam[28] = "INIT_COVMAT_SCATTER" ;
  int m, n, LDMP, i,j ;
  double epsilon = 0.00000000001; //tolerance level for differences
  double xsig, xred, xx;

  // ----------------- BEGIN --------------------

  if ( INPUTS.NCOVMAT_SCATTER <= 0 ) { return ; }

  print_banner("INIT_COVMAT_SCATTER:");

  LDMP = 1;

  // asign name to each scatter element based on model
  if ( INDEX_GENMODEL == MODEL_SALT2 ) {
    sprintf(GENLC.COVMAT_SCATTER_NAME[0],"%s", "mb" );
    sprintf(GENLC.COVMAT_SCATTER_NAME[1],"%s", "x1" );
    sprintf(GENLC.COVMAT_SCATTER_NAME[2],"%s", "c"  );
  }
  else {
    sprintf(c1err,"COVMAT_SCATTER not defined for GENMODEL=%s",
	    INPUTS.MODELNAME );
    sprintf(c2err,"%s", "Remove COVMAT_SCATTER* keys");
    errmsg( SEV_FATAL,0 , fnam, c1err, c2err);
  }



  //Abort if the same element of covmat
  //has been entered twice
  strcpy(c1err," Same element of covariance matrix entered twice");
  //char comb[10];
  for (m =0; m<3; m++){
    for (n=0; n<3;n++){
      sprintf(c2err, 
	      "INPUTS.COVMAT_SCATTER[%d][%d]", m , n);
      float t1 = fabs(INPUTS.COVMAT_SCATTER[m][n]);
      float t2 = fabs(INPUTS.COVMAT_SCATTER_SQRT[m][n]);
      float t3 = fabs(INPUTS.COVMAT_SCATTER_REDUCED[m][n]);
      if( (t1>0) && (t2>0)  )
	{ errmsg( SEV_FATAL,0 , fnam, c1err, c2err); }

      if ( (t1>0) && (t3>0)  ) 
	{ errmsg( SEV_FATAL,0 , fnam, c1err, c2err ); }
    }
  }
  

  //If INPUTS.COVMAT_SCATTER are negative abort
  strcpy(c1err," diagonal entries of covariance matrix negative");
  for (m = 0 ; m < 3 ; m++){
    sprintf(c2err, "INPUTS.COVMAT_SCATTER[%d][%d]", m , m);
    if (INPUTS.COVMAT_SCATTER[m][m] < -epsilon ) 
      { errmsg( SEV_FATAL,0 , fnam, c1err , c2err); }
  }

  //If values are given in terms of sqrt or reduced
  //scatter, convert to covmat_scatter
  for (m = 0 ; m < 3 ; m++){

    if (INPUTS.COVMAT_SCATTER_SQRT[m][m] > epsilon ){
      xsig = INPUTS.COVMAT_SCATTER_SQRT[m][m]; // sqrt(DIAG) = sigma
      INPUTS.COVMAT_SCATTER[m][m] = xsig * xsig;

      if ( LDMP ) 
	{ printf("\t Read SQRT[COV(%d,%d)] = %g \n",  m, m, xsig ); }
      
    }
  }


  //If values are given in terms of a reduced 
  //covariance matrix, convert to covmat_scatter
  for (m = 0; m < 3 ; m++){
    for (n = 0; n <3; n++){
      if (m == n) continue;

      xred = INPUTS.COVMAT_SCATTER_REDUCED[m][n] ;

      if (fabs(xred) > epsilon) {
	xx = INPUTS.COVMAT_SCATTER[m][m] * INPUTS.COVMAT_SCATTER[n][n] ;
	INPUTS.COVMAT_SCATTER[m][n] = xred * sqrt(xx);

	if ( LDMP ) 
	  { printf("\t Read RHO(%d,%d) = %g\n", m, n, xred ); }
      }

    }  // n
  }  // m
	

  // Test that intrinsic matrix is symmetric in case both 
  // off-diagonal elements have been provided.

  strcpy(c1err,"Symmetric off-diagonal entries entered and different for ");
  double diff ;
  int    Ldiff, Lmn, Lnm ;
  //	int mm, nn ;
  for (m = 0; m < 3 ; m++){
    for (n = 0;  n < 3 ; n++){
      sprintf(c2err, "INPUTS.COVMAT_SCATTER[%d][%d]", m , n);
      if (n <= m) continue;
      diff =  fabs(INPUTS.COVMAT_SCATTER[m][n] - 
	           INPUTS.COVMAT_SCATTER[n][m]);

      Ldiff = diff > epsilon ;
      Lmn   = fabs(INPUTS.COVMAT_SCATTER[m][n]) > epsilon ;
      Lnm   = fabs(INPUTS.COVMAT_SCATTER[n][m]) > epsilon ;

      if ( Ldiff && Lmn && Lnm ) 
	{  errmsg( SEV_FATAL,0 , fnam, c1err, c2err);  }

    }
  }    	



  //If only one of the asymmetric off-diagonal element
  //entered, put in the remaining terms.

  if(fabs(INPUTS.COVMAT_SCATTER[0][1])>epsilon){
    INPUTS.COVMAT_SCATTER[1][0] = INPUTS.COVMAT_SCATTER[0][1];
  } else {
    INPUTS.COVMAT_SCATTER[1][0] = INPUTS.COVMAT_SCATTER[1][0];
  }
  
  if (fabs(INPUTS.COVMAT_SCATTER[2][0])>epsilon){
    INPUTS.COVMAT_SCATTER[0][2] = INPUTS.COVMAT_SCATTER[2][0];
  } else {
    INPUTS.COVMAT_SCATTER[2][0] = INPUTS.COVMAT_SCATTER[0][2];
  }	
  if (fabs(INPUTS.COVMAT_SCATTER[2][1])>epsilon){
    INPUTS.COVMAT_SCATTER[1][2] = INPUTS.COVMAT_SCATTER[2][1];
  }else {
    INPUTS.COVMAT_SCATTER[2][1] = INPUTS.COVMAT_SCATTER[1][2];
  }

  //Have INPUTS.COVMAT_SCATTER completely. 

  //Check that submatrices have non-singular dets
  strcpy(c1err,"Submatrices are not pos-definite for combination ");

  for (m =0; m < 3 ; m++){
    n = m + 1; 
    
    if (3 == n ) {n = 0;}
    
    diff = (INPUTS.COVMAT_SCATTER[m][m]*
	    INPUTS.COVMAT_SCATTER[n][n]) - 
      (INPUTS.COVMAT_SCATTER[m][n]*
       INPUTS.COVMAT_SCATTER[m][n]);

    sprintf(c2err, "%d\t%d", m , n);
    
    if (diff < 0) {
      print_preAbort_banner(fnam);
      printf("\t m=%d    n=%d    diff=%le\n", m, n , diff);
      errmsg( SEV_FATAL,0 , fnam, c1err, c2err);
    }
  }

  int N = -1;
  char *cptr ;
  N++ ; cptr = GENLC.COVMAT_SCATTER_README[N] ;
  sprintf(cptr,"%s", "\t Cov matrix:");


  N++ ; cptr = GENLC.COVMAT_SCATTER_README[N] ;
  sprintf(cptr,"\t\t %8s(%d) %8s(%d) %8s(%d) "
	  ,GENLC.COVMAT_SCATTER_NAME[0], 0
	  ,GENLC.COVMAT_SCATTER_NAME[1], 1
	  ,GENLC.COVMAT_SCATTER_NAME[2], 2 );

  N++ ; cptr = GENLC.COVMAT_SCATTER_README[N] ;
  sprintf(cptr,"%s", "\t -----------------------------------------------");
  
  
  for (m = 0 ; m < 3 ; m++){
    N++ ; cptr = GENLC.COVMAT_SCATTER_README[N] ;
    sprintf(cptr,"\t %-4s(%d) :  %10.6f  %10.6f  %10.6f"
	    ,GENLC.COVMAT_SCATTER_NAME[m], m
	    ,INPUTS.COVMAT_SCATTER[m][0]
	    ,INPUTS.COVMAT_SCATTER[m][1]
	    ,INPUTS.COVMAT_SCATTER[m][2]  );

  }
  
  N++ ; cptr = GENLC.COVMAT_SCATTER_README[N] ;
  cptr[0] = 0 ;
	
  for ( i=0; i <= N; i++ ) {
    printf("%s\n", GENLC.COVMAT_SCATTER_README[i] );
  }

  //printvals();
  //Do the Cholesky Decomposition once and for all
  double a_data[] = { 
    INPUTS.COVMAT_SCATTER[0][0],
    INPUTS.COVMAT_SCATTER[0][1],
    INPUTS.COVMAT_SCATTER[0][2],
    INPUTS.COVMAT_SCATTER[1][0],
    INPUTS.COVMAT_SCATTER[1][1],
    INPUTS.COVMAT_SCATTER[1][2],
    INPUTS.COVMAT_SCATTER[2][0],
    INPUTS.COVMAT_SCATTER[2][1],
    INPUTS.COVMAT_SCATTER[2][2]
  };

  gsl_matrix_view chk;
  chk  = gsl_matrix_view_array (a_data, 3, 3); 
  gsl_linalg_cholesky_decomp ( &chk.matrix)  ;
  
  for (i =0; i<3; i++){
    for (j = 0; j < 3 ; j++){      
      if (i <= j ) { 
	INPUTS.COVMAT_CHOLESKY_SCATTER[i][j] =  
	  gsl_matrix_get(&chk.matrix,i,j); 
      }
      else {
	INPUTS.COVMAT_CHOLESKY_SCATTER[i][j] = 0.0 ;
      }
    }
    
  }
  //printvals();
  

  return;	
	
} // END OF INIT_COVMAT_SCATTER



// ========================================================
int GEN_COVMAT_SCATTER ( double *randoms, double *SCATTER_VALUES  ) {

  // Created July 2011 by R.Biswas
  // Ensure Initialization of scatter vector;
  //
  // *randoms        (I) : 3 gaussian random numbers
  // *SCATTER_VALUES (O) : ..
  //

  SCATTER_VALUES[0] = 0. ;
  SCATTER_VALUES[1] = 0. ;
  SCATTER_VALUES[2] = 0. ;
  
  if ( INPUTS.NCOVMAT_SCATTER <= 0 ) { return 0 ; }

  //Draw 3 vector from a standard (uncorrelated multinormal distribution)
  int i , j;
  //double sigma =1.0;
  double normalvector[3];
  
  //Form 3 vector of standard normal RV
  for (i =0; i <3; i++ ){
    normalvector[i] = randoms[i];    
  }

  //Matrix Multiply
  //scatter_values = ch^T normalvector
  for (i = 0 ; i < 3 ; i++){
    for (j = 0 ; j < 3 ; j++){
      //transpose cholesky matrix
      double tmp = INPUTS.COVMAT_CHOLESKY_SCATTER[j][i]; 
      SCATTER_VALUES[i] += tmp * normalvector[j];      
    }
  }
  return 0;

} //END OF GEN_COVMAT_SCATTER


// *******************************
void INIT_FUDGE_SNRMAX(void) {


  // Created May 2015 by R.Kessler
  // Parse INPUTS.STRING_FUDGE_SNRMAX
  
  //  char fnam[] = "INIT_FUDGE_SNRMAX" ;
  char *ptrStr, f1[2], strSNRMAX[20] ;

  // -------------- BEGIN -------------

  if ( INPUTS.OPT_FUDGE_SNRMAX == 0 ) { return ; }

  ptrStr = INPUTS.STRING_FUDGE_SNRMAX ;
  if ( strstr(ptrStr,"=") != NULL ) {
    sprintf(f1, "%c", ptrStr[0] ) ; // single-char band
    sprintf(strSNRMAX, "%s", (ptrStr+2) );
    INPUTS.IFILTOBS_FUDGE_SNRMAX = INTFILTER(f1);
  }
  else {
    sprintf(strSNRMAX, "%s", ptrStr );
  }

  // read SNRMAX part of string into float
  sscanf(strSNRMAX, "%f", &INPUTS.FUDGE_SNRMAX); 

} // end of INIT_FUDGE_SNRMAX

// ********************************************
void init_covar_mlcs2k2(void) {

  /*******
   Initialize covariance matrix used to generate MLCS model.
   First version (9/2006) assumes 100% correlation among epochs
   within a bandpass, and zero correlation between bandpasses.
   The correlation model may get more sophisticate later.

   Note that this uses "gencovar_mlcs2k2" to fetch pieces
   of the real MLCS smatrix.


  History
  ~~~~~~~~
  Dec 30, 2007:  fill GENLC.PEAKMAGERR_MODEL[ifilt] 

  May 28,  2009:  don't write covFile because this is called
                  before the [VERSION]/ subdir is created,
                  and because nobody uses this COV file.

  *******/

  int 
    NFILT, MINEP, MAXEP, NEP, NMAT
    ,iep, ifilt,ifilt_rest, iep1, iep2, ifilt1, ifilt2
    ;

  double epoch8, sqerr8    ;
  float MODELERR[MXFILTINDX][MXEPCOV] ;
  float err1, err2, errtmp,T1 ;
  char fnam[] = "init_covar_mlcs2k2" ;

  // --------------- BEGIN ---------


  //  if ( INPUTS.GENMODEL_ERRSCALE <= 0.000001 ) return

  NMAT = 1 ;      // size of covar matrix to fetch
  MINEP = mlcs2k2_Tmin() ;
  MAXEP = mlcs2k2_Tmax() ;
  NEP   = MAXEP - MINEP + 1 ;

  sprintf(BANNER,"Init covariance matrix for sim: %d to %d days.", 
	  MINEP, MAXEP);
  print_banner(BANNER);

  if ( NEP >= MXEPCOV ) {
    sprintf(c1err,"NEP = %d exceeds array bound (MXEPCOV=%d)", NEP, MXEPCOV);
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }



  NFILT   = GENLC.NFILTDEF_REST; 

  if ( NFILT > MXFILT_COVAR ) {
    sprintf(c1err,"NFILT=%d exceeds COVAR array bound MXFILT=%d", 
	    NFILT, MXFILT_COVAR );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }

  for ( ifilt = 0 ; ifilt < NFILT ; ifilt++ ) {

    ifilt_rest = GENLC.IFILTMAP_REST[ifilt];

    for ( iep = MINEP; iep <= MAXEP ; iep++ ) {
      epoch8 = (double)iep ;
      gencovar_mlcs2k2(NMAT, &ifilt, &epoch8, &sqerr8 );  // get sqerr8
      errtmp  = (float)sqrt(sqerr8);
      MODELERR[ifilt][iep-MINEP] = errtmp ;

      T1 = (float)(iep);
      
      // store errors at peak for each filter.
      if ( fabsf(T1) < 0.5 )
	{ GENLC.PEAKMAGERR_MODEL[ifilt_rest] = errtmp ; }


    } // end of iep loop
  }  // end of ifilt loop


  // now construct full covar matrix using error at each epoch

  for ( ifilt1=0; ifilt1 < NFILT; ifilt1++ ) {
    for ( iep1=0; iep1 < NEP; iep1++ ) {
      for ( ifilt2=0; ifilt2 < NFILT; ifilt2++ ) {
	for ( iep2=0; iep2 < NEP; iep2++ ) {
   
	  GENLC.COVAR[ifilt1][iep1][ifilt2][iep2] = 0.0 ;  // init to zero

	  // set 100% correlation only within same passband
	  if ( ifilt1 == ifilt2 ) {
	    err1 = MODELERR[ifilt1][iep1] ;
	    err2 = MODELERR[ifilt2][iep2] ;
	    GENLC.COVAR[ifilt1][iep1][ifilt2][iep2] = err1 * err2 ;
	  }

	}
      }
    }
  }


}  // end of init_covar_mlcs2k2



// *********************************************
void readme_doc(int iflag_readme) {

  // Created 2007 by R.Kessler
  // fill VERSION_INFO.README_DOC structure with 
  // "VERSION_INFO.NLINE_README" lines of README content.
  // Note that the contents are printed from elsewhere
  //
  // iflag_readme = 1 => fill init-part
  // iflag_readme = 2 => fill post-sim part
  //
  // Jul 24, 2009: re-write mag-smearing stuff, 
  //               and include new GENMAG_SMEAR_FILTER
  //
  // Dec 20, 2009: print out SIMSED parameters and ranges
  // Jul 01, 2010: print lambda-dependent mag smearing stuff
  //
  // Oct 29, 2010: print nonIa info for both non1a and NON1A;
  //
  // Mar 2011 : write FUDGE_SNRMAX, FUDGE2_SNRMAX if they are non-zero.
  //
  // Jul 19, 2011: clarify APPLY_SEARCHEFF_OPT
  //
  // Mar 14, 2012: write SMEAR_FUNPAR
  //
  // Jan 19 2016: remove redundant NON1A dump at beginning.
  // Feb 01 2017: call readme_doc_NON1ASED
  // Jan 16 2019: always print KCOR file (before, only printed for SNIA)
  // May 27 2019: print PEAKMJD-estimate method
  // Feb 16 2020: write SIMULATION key at top for Pippin.
  // Feb 20 2020: write PHOTPROB info
  // Jul 31 2020: for batch jobs (JOBID>0), write keys for monitor task
  // Aug 26 2020: use DOCANA structure at top of file

  char ctmp[MXPATHLEN], cfilt[2], cwd[MXPATHLEN] ;
  char *cptr;
  char conoff[2][4] = { "OFF" , "ON" } ;

  int i, j, j2, ifilt_obs, ifilt_rest, ifilt, itmp, imap, iopt, ipar ;
  int NLINE, NOV, NON1A_non1a, OVP1, OVP2 ;

  double XN, XNERR ;
  float xt, xtprod, val, ZMIN, ZMAX, shift[2];

  char  fnam[] = "readme_doc" ;

  // ------------ BEGIN readme_doc() ---------------

  i=0 ;

  print_banner ( " Fill comments for README doc-file" );

  if ( iflag_readme == 2 ) goto AFTERSIM ;

  //--- brief description

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s \n", KEYNAME_DOCANA_REQUIRED ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SURVEY:       %s\n",  GENLC.SURVEY_NAME);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  GENMODEL:     %s \n", INPUTS.MODELNAME);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  HOST_MACHINE: %s \n", getenv("HOST") );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf( cptr, "  USERNAME:  %s \n", getenv("USER") );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SNDATA_ROOT:  %s \n", PATH_SNDATA_ROOT );
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SNANA_DIR:     %s \n", PATH_SNANA_DIR );
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SNANA_VERSION: %s \n", SNANA_VERSION_CURRENT );

  // write current directory (Sep 5 2013)
  if ( getcwd(cwd,MXPATHLEN) != NULL ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"  CWD:   %s \n", cwd );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s \n", KEYNAME2_DOCANA_REQUIRED ); 

  // -----------------------------

  // indicate changes if "GENPERFECT" is requested.
  readme_doc_GENPERFECT(&i);

  // ---- full description

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n FULL_DESCRIPTION: \n" ); 

  readme_doc_SIMLIB(&i);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Generation VERSION: %s \n", INPUTS.GENVERSION ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Generation source : %s \n", INPUTS.GENSOURCE ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Generation FILTERS: %s \n", INPUTS.GENFILTERS ); 

  if ( GENFRAME_OPT == GENFRAME_REST ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Rest-frame FILTERS: %s  for  model=%s\n", 
	    GENLC.FILTLIST_REST, INPUTS.MODELNAME );

  }
  // ------------
  // print info for skipped filters

  readme_doc_filterWarn(&i);

  NON1A_non1a = ( INDEX_GENMODEL == MODEL_NON1ASED ) ;
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t KCOR lookup tables: %s \n", INPUTS.KCOR_FILE ); 
    

  if ( GENFRAME_OPT == GENFRAME_REST ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t KCOR lookup  time : Trest" );
    strcat ( cptr, "\n" );

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t AVwarp option for SED : color from " );
    if ( INPUTS.KCORFLAG_COLOR == 1 ) 
      { strcat ( cptr, "closest passbands \n" ); }
    else if ( INPUTS.KCORFLAG_COLOR == 2 ) 
      { strcat ( cptr, "Jha table \n" ); }
    else
      { strcat ( cptr, "UNKNOWN\n" ); }
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( INPUTS.SMEARFLAG_FLUX>0 ) { j=1; } else { j=0; }
  sprintf(cptr, "\t Flux-smearing is %s \n", conoff[j] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Reported flux-uncertainty includes " ); 
  j2 = (INPUTS.SMEARFLAG_FLUX & 2);
  if ( j2 == 0 ) 
    { strcat(cptr,"SKY+GALAXY+SOURCE\n"); }
  else
    { strcat(cptr,"SKY only\n"); } // SMP-like errors


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if( INPUTS.SMEARFLAG_ZEROPT > 0 ) { j=1; } else { j=0; }
  sprintf(cptr, "\t Zeropt-smearing is %s \n", conoff[j] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;  j=0;
  OVP1 = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_PHOT ; 
  if ( OVP1 ) { j = 1; }
  sprintf(cptr, "\t Host-galaxy shot-noise  is %s \n", conoff[j] );


  i++; cptr = VERSION_INFO.README_DOC[i] ;  j=0;
  OVP2 = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_IMAGE ; 
  if ( OVP2 ) { j = 1; }
  sprintf(cptr, "\t Host-galaxy image-noise  is %s \n", conoff[j] );


  readme_doc_MWXT(&i);


  // ------- RATE INFO ---------

  for ( j=0; j < NLINE_RATE_INFO; j++ ) {
     i++; cptr = VERSION_INFO.README_DOC[i] ;
     sprintf(cptr,"%s\n", LINE_RATE_INFO[j] );
  }

  // ------- generation --------

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  GENERATION RANGES: \n");

  if ( INPUTS.USE_SIMLIB_GENOPT ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t ==> Override gen-values with SIMLIB header values.\n");
  }

  if ( INPUTS.USE_SIMLIB_REDSHIFT == 0 ) 
    { sprintf(ctmp,"%s distribution.", INPUTS.RATEPAR.NAME ); }
  else
    { sprintf(ctmp,"SIMLIB values." ); }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,
	  "\t Generate Redshift : %6.3f to %6.3f  using  %s \n"
	  ,INPUTS.GENRANGE_REDSHIFT[0]
	  ,INPUTS.GENRANGE_REDSHIFT[1] 
	  ,ctmp	  );


  if ( INPUTS.GENBIAS_REDSHIFT != 0.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Generate redshift bias : %f \n", 
	    INPUTS.GENBIAS_REDSHIFT );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( INPUTS.GENSIGMA_REDSHIFT >= 0.0 ) {
    sprintf(cptr,"\t REDSHIFT_FINAL is ZCMB_GEN smeared by : %8.5f \n", 
	    INPUTS.GENSIGMA_REDSHIFT);
  }
  else
    { sprintf(cptr,"\t REDSHIFT_FINAL is host-photoZ \n"); }

  
  if ( WRONGHOST.NLIST > 0.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t WRONGHOST model frac: %.3f + %.3f*z + %.5f*z^2 \n"
	    ,WRONGHOST.PROB_POLY[0]
	    ,WRONGHOST.PROB_POLY[1]
	    ,WRONGHOST.PROB_POLY[2]     );
  }

  if ( INPUTS.VEL_CMBAPEX == 0.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,
	    "\t v(cmb)=0 -> REDSHIFT_HELIO=REDSHIFT_CMB (legacy option)\n");
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  bool  USE_HOSTLIB_VPEC = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEVPEC );
  if ( USE_HOSTLIB_VPEC  ) {
    sprintf(cptr,"\t Peculiar Velocity (VPEC) HOSTLIB RMS : %.1f km/sec\n", 
	    HOSTLIB.VPEC_RMS );
  }
  else {
    sprintf(cptr,"\t Peculiar Velocity (VPEC) Gauss sigma: %.1f km/sec\n", 
	    INPUTS.GENSIGMA_VPEC );
  }
  if ( INPUTS.RESTORE_WRONG_VPEC ) {
    i++ ; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t\t (RESTORE wrong VPEC sign convention)\n") ;
  }


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %s   ZP   offsets : ", INPUTS.GENFILTERS);
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    sprintf(ctmp,"%6.3f ", INPUTS.GENMAG_OFF_ZP[ifilt_obs] ); 
    strcat(cptr,ctmp); 
  }
  strcat(cptr,"\n"); 


  float magoff;
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %s  MODEL offsets : ", INPUTS.GENFILTERS);
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    magoff = INPUTS.GENMAG_OFF_MODEL[ifilt_obs] + INPUTS.GENMAG_OFF_GLOBAL ;
    sprintf(ctmp,"%6.3f ", magoff ); 
    strcat(cptr,ctmp); 
  }
  strcat(cptr,"\n"); 


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %s  exposure times: ", INPUTS.GENFILTERS);
  xtprod = 1.0 ;
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    xt = INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] ;
    xtprod *= xt ;

    if ( xt > 10.0 )   sprintf(ctmp,"%6.1f ", xt ); 
    else               sprintf(ctmp,"%6.4f ", xt ); 

    strcat(cptr,ctmp); 
  } // end if ifilt
  strcat(cptr,"\n"); 

  // indicate which quantities were scaled by EXPOSURE_TIME  
  if ( xtprod != 1.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t %s  exposure MSKOPT=%d => ", 
	    INPUTS.GENFILTERS, INPUTS.EXPOSURE_TIME_MSKOPT );
    
    if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 0) )
      { strcat(cptr, "ZPT  "); }
    if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 1) )
      { strcat(cptr, "SKYSIG  "); }
    if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 2) )
      { strcat(cptr, "READNOISE  "); }

    strcat(cptr,"\n");
  }

  // - - - - - - -  -
  // optional MJD_TEMPLAT (Sep 2017)
  if ( INPUTS.USE_MJD_TEMPLATE ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t MJD_TEMPLATE: " );
    
    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      val       = INPUTS.MJD_TEMPLATE_FILTER[ifilt_obs] ;
      if ( val > 1.0 ) {
	sprintf(ctmp,"%6.1f(%c) ", val, FILTERSTRING[ifilt_obs] ); 
	strcat(cptr,ctmp); 
      }
    } // end if ifilt
    strcat(cptr,"\n"); 
  }

  // - - - - - - -
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t RA       : %6.2f to %6.2f  deg\n", 
	  INPUTS.GENRANGE_RA[0], INPUTS.GENRANGE_RA[1] );

  // - - - -  PEAKMJD stuff - - - - - 
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t PEAKMJD  : %8.1f to %8.1f   \n", 
	  INPUTS.GENRANGE_PEAKMJD[0], INPUTS.GENRANGE_PEAKMJD[1] );

  if ( INPUTS.GENSIGMA_PEAKMJD > 1.0E-9 ) {
    sprintf(ctmp,"Gauss smear, sigma=%5.2f days", INPUTS.GENSIGMA_PEAKMJD);
  }
  else if ( (INPUTS.OPT_SETPKMJD & OPTMASK_SETPKMJD_FLUXMAX2)>0 ) {
    sprintf(ctmp,"Fmax-clump, MJDWIN=%.1f, SNRCUT>%.1f(3.0)",
	    INPUTS.MJDWIN_SETPKMJD, INPUTS.SNRCUT_SETPKMJD );
  }
  else if ( (INPUTS.OPT_SETPKMJD & OPTMASK_SETPKMJD_FLUXMAX)>0 ) {
    sprintf(ctmp,"naive max flux.");
  }
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t PEAKMJD-estimate  : %s\n", ctmp);

  // - - - - - - - - - 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Trest    : %8.2f to %8.2f  days \n", 
	  INPUTS.GENRANGE_TREST[0], INPUTS.GENRANGE_TREST[1] );

  
  if ( INPUTS.TGRIDSTEP_MODEL_INTERP > 0.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t TGRIDSTEP  : %.1f days (for model-mag interp) \n",
	    INPUTS.TGRIDSTEP_MODEL_INTERP );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t RISETIME-SHIFT(days) SIGMA(lo,hi) : %3.1f , %3.1f  (Mean= %3.1f) \n"
	  ,INPUTS.GENGAUSS_RISETIME_SHIFT.SIGMA[0]
	  ,INPUTS.GENGAUSS_RISETIME_SHIFT.SIGMA[1] 
	  ,INPUTS.GENGAUSS_RISETIME_SHIFT.PEAK
	  );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t FALLTIME-SHIFT(days) SIGMA(lo,hi) : %3.1f , %3.1f  (Mean= %3.1f) \n"
	  ,INPUTS.GENGAUSS_FALLTIME_SHIFT.SIGMA[0]
	  ,INPUTS.GENGAUSS_FALLTIME_SHIFT.SIGMA[1] 
	  ,INPUTS.GENGAUSS_FALLTIME_SHIFT.PEAK
	  );

  // print shape-par info for all models except SIMSED
  if ( INPUTS.NPAR_SIMSED == 0  ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(ctmp,"\t Shape-par(%s)", GENLC.SHAPEPAR_NAME );
    sprintf_GENGAUSS(cptr, ctmp, &INPUTS.GENGAUSS_SHAPEPAR);
  }

  // SIMSED parameters
  readme_doc_SIMSED(&i);

  // ---- SALT2 params
  readme_doc_SALT2params(&i);

  // ---- GENPDF populations (Aug 2021)
  readme_doc_GENPDF(&i);

  // ---- FIXMAG params
  readme_doc_FIXMAG(&i);


  // --------------------------------
  // dump host extinction params

  if ( INPUTS.GENPROFILE_AV.USE )
    { readme_doc_hostxt(&i, &INPUTS.GENPROFILE_AV); }
  else if ( INPUTS.GENPROFILE_EBV_HOST.USE ) 
    { readme_doc_hostxt(&i, &INPUTS.GENPROFILE_EBV_HOST); }
  else {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n  Host Extinction Parameters: NONE  (AV=0) \n");    
  }



  // ====================================
  // Z-dependent parameters (if requested)

  if ( NPAR_ZVAR_USR > 0 ) {
    
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\n  Z-dependent SN-parameter shifts: \n");

    for ( ipar=0; ipar < NPAR_ZVAR_USR; ipar++ ) {
      sprintf(ctmp,"%s", INPUT_ZVARIATION[ipar].PARNAME );
      ZMIN = INPUTS.GENRANGE_REDSHIFT[0] ; 
      ZMAX = INPUTS.GENRANGE_REDSHIFT[1] ; 
      shift[0] = get_zvariation(ZMIN,ctmp) ;
      shift[1] = get_zvariation(ZMAX,ctmp) ;

      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr," %20s-shift = %6.3f(ZMIN=%5.3f), %6.3f(ZMAX=%5.3f) \n",
	      ctmp, shift[0],ZMIN, shift[1], ZMAX );
    }  // end of NPAR_ZVAR_USR loop
  } else {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\n  Z-dependent SN-parameter shifts:  None. \n");
  }


  // ===========================================

  readme_doc_magSmear(&i);

  readme_doc_nonLin(&i); // May 27 2016

  // ============================================
  // AVWARP-overflows vs. rest-frame filter

  i++; cptr = VERSION_INFO.README_DOC[i] ;

  sprintf(ctmp, "\n  AVWARP_OVERFLOWS: ");
  if ( NAVWARP_OVERFLOW[0] == 0 ) {
    strcat(ctmp," NONE. ");
  }
  else {
    
    for ( ifilt=0; ifilt < GENLC.NFILTDEF_REST; ifilt++ ) {
      ifilt_rest = GENLC.IFILTMAP_REST[ifilt] ;
      sprintf(cfilt,"%c", FILTERSTRING[ifilt_rest] ); 
      NOV = NAVWARP_OVERFLOW[ifilt_rest] ;
      if ( NOV > 0 ) sprintf(ctmp, "%s %d(%s)", ctmp, NOV, cfilt);
    } 
  }
  strcat(ctmp,"\n");
  strcat(cptr,ctmp);
  sprintf(WARNING_AVWARP_OVERFLOW,"\n  WARNING: %s", ctmp); 

  readme_doc_TAKE_SPECTRUM(&i);

  // ----- cosmology parameters

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Cosmology Parameters: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t H0 = %6.2f km/s per MPc \n", INPUTS.H0 );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Omega_{M,L} = %6.3f, %6.3f   w0,wa = %5.2f,%5.3f  \n",
	  INPUTS.OMEGA_MATTER, INPUTS.OMEGA_LAMBDA, 
	  INPUTS.w0_LAMBDA, INPUTS.wa_LAMBDA );



  // ------ software Search efficiency

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n --------------------------------------------------- \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  Software-Pipeline Search Efficiency (MINOBS=%d) from \n", 
	  INPUTS_SEARCHEFF.MINOBS );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s", SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[0] ) ;
  for ( imap=0; imap < INPUTS_SEARCHEFF.NMAP_DETECT; imap++ ) {
    NLINE = SEARCHEFF_DETECT[imap].NLINE_README ; 
    for ( itmp = 0; itmp < NLINE; itmp++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr, "%s", SEARCHEFF_DETECT[imap].README[itmp] ) ;
    }    
  }

  for ( imap=0; imap < INPUTS_SEARCHEFF.NMAP_PHOTPROB; imap++ ) {
    NLINE = SEARCHEFF_PHOTPROB[imap].NLINE_README ; 
    for ( itmp = 0; itmp < NLINE; itmp++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr, "%s", SEARCHEFF_PHOTPROB[imap].README[itmp] ) ;
    }    
  }

  // print detection logic
  NLINE = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ; 
  for ( itmp = 0; itmp < NLINE; itmp++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[itmp] ) ;
  }    
  

  // ------ Spec Search efficiency

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s", "\n  Spectroscopic Efficiency : \n" );
  for ( iopt=0; iopt < SEARCHEFF_SPEC_INFO.NLINE_README; iopt++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr, "%s \n", SEARCHEFF_SPEC_INFO.README[iopt] ) ;    
  }
  

  if ( INPUTS_SEARCHEFF.NMAP_zHOST ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", "\n  Unconfirmed zHOST Efficiency map from \n" );
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t %s \n", INPUTS_SEARCHEFF.zHOST_FILE );
  }
  else {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", "\n  Unconfirmed zHOST Efficiency : 100% \n" );
  }
  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  %s \n", COMMENT_README_TRIGGER);

  // print SNTYPE values for SPEC and PHOT Ia-subsets 
  // For NONIA there is only one type specified with the NONIA keys.
  if ( LGEN_SNIA ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"  SNTYPE(Ia) = %d(SPEC)  and %d(PHOT) \n", 
	    INPUTS.SNTYPE_Ia_SPEC, INPUTS.SNTYPE_Ia_PHOT);    
  }


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr," --------------------------------------------------- \n");


  // ------ software cuts -----------
  if ( INPUTS.APPLY_CUTWIN_OPT  ) {  readme_doc_CUTWIN(&i) ; }


  // ---- HOST-GALAXY 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  HOSTLIB Summary: " );

  if ( INPUTS.HOSTLIB_USE == 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"None. \n" );
  }
  else {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n" );
    NLINE = HOSTLIB.NLINE_COMMENT ;
    for ( itmp=0; itmp < NLINE; itmp++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t %s \n", HOSTLIB.COMMENT[itmp] );
    }
  }


  // -----  FUDGES on observing conditions ------
  readme_doc_FUDGES(&i);

  // write list of MAPs so that pipelines can check time-stamps, etc ...
  readme_doc_mapFileList(&i);


  // ======================================
  VERSION_INFO.NLINE_README_INIT = i;   // can dump to here after init
  if ( iflag_readme == 1 ) return ;
  // ======================================

 AFTERSIM:
  i = VERSION_INFO.NLINE_README_INIT ;

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n ============ END OF SIMULATION SUMMARY ============== \n");

  // ---- random seed and first/last random

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Random Number Sync: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t RANDOM SEED: %d   (RANLIST_START_GENSMEAR: %d)\n", 
	  INPUTS.ISEED, INPUTS.RANLIST_START_GENSMEAR );


  int ilist;
  sumstat_RANLISTs(2);
  for ( ilist=1; ilist <= GENRAN_INFO.NLIST_RAN; ilist++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"   FIRST/LAST Random Number (List=%d): %f %f  "
	    "AVG(wrap) = %.1f +_ %.1f \n", ilist, 
	    GENRAN_INFO.RANFIRST[ilist], GENRAN_INFO.RANLAST[ilist],
	    GENRAN_INFO.NWRAP_AVG[ilist], GENRAN_INFO.NWRAP_RMS[ilist]	);
  }

  // ---- statistics

  double t_gen   = (TIMERS.t_end - TIMERS.t_end_init); // total time after init
  double R_gen   = (double)NGENLC_TOT / t_gen ;  // NGEN/sec
  double R_write = (double)NGENLC_WRITE/t_gen ;  // NWRITE/sec

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Generation Statistics (gen CPU=%.1f minutes): \n", 
	  t_gen/60.);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Generated %5d simulated light curves "
	  "(%.f/sec) \n",  NGENLC_TOT, R_gen );
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Wrote     %5d simulated light curves to SNDATA files "
	  "(%.f/sec). \n",  NGENLC_WRITE, R_write );

  if ( NGENSPEC_WRITE > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Wrote     %5d simulated spectra to SNDATA files \n",
	    NGENSPEC_WRITE );
  }

  // spectroscopic tags

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  Spectroscopic-type: %d -> %d (before -> after cuts)\n",
	  GENLC.NTYPE_SPEC, GENLC.NTYPE_SPEC_CUTS);
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  Photometric-type:   %d -> %d (before -> after cuts)\n",
	  GENLC.NTYPE_PHOT, GENLC.NTYPE_PHOT_CUTS);


  if ( !IGNOREFILE(INPUTS.WRONGHOST_FILE) ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    double frac = 0.0 ;
    int N0      = NGENLC_WRITE ;
    int N1      = GENLC.NTYPE_PHOT_WRONGHOST ;
    if ( GENLC.NTYPE_PHOT_CUTS > 0 ) { frac = (double)N1 / (double)N0 ; }
    sprintf(cptr,"  Wrong-host fraction: %d/%d = %.4f\n", N1,N0,frac);
  }

  // ----------------

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  Rejection Statistics: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %5d rejected by NEPOCH<%d \n",
	  NGEN_REJECT.NEPOCH, (int)INPUTS.CUTWIN_NEPOCH[0] ) ;


  double MAGMIN = INPUTS.GENRANGE_PEAKMAG[0];
  double MAGMAX = INPUTS.GENRANGE_PEAKMAG[1];
  if ( MAGMIN > 0.0 && MAGMAX < 9999. ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t %5d rejected by GENRANGE_PEAKMAG(%.1f to %.1f) \n",  
	    NGEN_REJECT.GENMAG, MAGMIN, MAGMAX );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %5d rejected by GENRANGEs \n",  
	  NGEN_REJECT.GENRANGE );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %5d rejected by SEARCH-TRIGGER \n",  
	  NGEN_REJECT.SEARCHEFF );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %5d rejected by CUTWIN-SELECTION \n",  
	  NGEN_REJECT.CUTWIN );

  // ---

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SEARCH+CUTS Efficiency: %7.4f +- %7.4f \n", 
	  GENLC.GENEFF, GENLC.GENEFFERR);

  // give warning if generation stops early
  if ( GENLC.STOPGEN_FLAG  ) {

    bool QUIT_NOREWIND=((INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_QUIT_NOREWIND)>0);
  
    if ( QUIT_NOREWIND ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,
	      "\n  WARNING: GENERATION STOPPED AFTER ONE PASS THRU SIMLIB\n");
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t  AS REQUESTED BY SIM-INPUT SIMLIB_MSKOPT += %d\n",
	      SIMLIB_MSKOPT_QUIT_NOREWIND );
    }
    else {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\n  WARNING: GENERATION STOPPED WHEN ERROR\n");
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t   on ERR(EFF) <= EFFERR_STOPGEN(=%f)\n",  
	      INPUTS.EFFERR_STOPGEN );
      
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t  => YOU HAVE %d FEWER LIGHT CURVES THAN REQUESTED\n",
	      INPUTS.NGEN_LC - NGENLC_WRITE );
    }
  }


  // give SN-stats with cuts
  if ( NLINE_RATE_INFO > 0 ) {
    XN   = (INPUTS.RATEPAR.SEASON_COUNT  + INPUTS.RATEPAR_PEC1A.SEASON_COUNT) ;
    XN  *= GENLC.GENEFF ; // multiply by cut-efficiency
    if ( NGENLC_WRITE > 0 ) 
      { XNERR = XN/sqrt((double)NGENLC_WRITE); }
    else
      { XNERR = 0.0 ; }
    
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n  Number of SNe per season AFTER CUTS : "
	    "%6.0f +- %5.0f \n", XN, XNERR );
  }


  // NON1ASED
  if ( INPUTS.NON1ASED.NINDEX > 0 ) { readme_doc_NON1ASED(&i); }

  // end marker

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n\t ===== END OF README FILE ====== \n");

  // ------
  if ( i >= MXDOCLINE ) {
    sprintf ( c1err, "%d README.DOC lines exceeds array bound of %d",
	      i, MXDOCLINE );
    sprintf ( c2err,"Increase parameter MXDOCLINE and re-make ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  VERSION_INFO.NLINE_README  = i;


}  // end of readme_doc


// ********************************************
void readme_doc_CUTWIN(int *iline) {
  
  int i, icut ;
  char *cptr, ctmp[80] ;
  // ---------- BEGIN ---------

  i = *iline ;

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  SOFTWARE CUTS: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t EPOCH CUT: %5.0f Lambda(rest) < %5.0f A \n",
	  INPUTS.EPCUTWIN_LAMREST[0], INPUTS.EPCUTWIN_LAMREST[1] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t EPOCH CUT: SNR >= %2.0f  \n",
	  INPUTS.EPCUTWIN_SNRMIN[0] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t TrestMIN < %5.1f  &&  TrestMAX > %4.1f days \n",
	  INPUTS.CUTWIN_TRESTMIN[1], INPUTS.CUTWIN_TRESTMAX[0] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Max TGAP(rest) <= %5.1f  days \n",
	  INPUTS.CUTWIN_TGAPMAX[1] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Max T0GAP(rest) <= %5.1f  days \n",
	  INPUTS.CUTWIN_T0GAPMAX[1] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t NOBS(MJDDIF > %3.1f) >= %2d \n",
	  INPUTS.CUTWIN_MJDDIF[0], INPUTS.CUTWIN_NOBSDIF[0] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t NEPOCH(SNR > %4.1f) >= %2.0f \n",
	  INPUTS.CUTWIN_NEPOCH[1], INPUTS.CUTWIN_NEPOCH[0] ) ;

  for ( icut=1; icut <= INPUTS.NCUTWIN_SNRMAX; icut++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SNRMAX > %4.1f for %d of the '%s' filters "
	    "(%5.1f < Trest < %5.1f) \n"
	    , INPUTS.CUTWIN_SNRMAX[icut][0]
	    , INPUTS.CUTWIN_SNRMAX_NFILT[icut]
	    , INPUTS.CUTWIN_SNRMAX_FILTERS[icut]
	    , INPUTS.CUTWIN_SNRMAX_TREST[icut][0]
	    , INPUTS.CUTWIN_SNRMAX_TREST[icut][1]
	    );
  }


  for ( icut=0; icut < INPUTS.NCUTWIN_SATURATE; icut++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t NOBS_SATURATE = %d to %d for each %s filter. \n"
	    , INPUTS.CUTWIN_SATURATE_NOBS[icut][0]
	    , INPUTS.CUTWIN_SATURATE_NOBS[icut][1]
	    , INPUTS.CUTWIN_SATURATE_FILTERS[icut]	    );
  }
  for ( icut=0; icut < INPUTS.NCUTWIN_NOSATURATE; icut++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t NOBS_NOSATURATE = %d to %d for each %s filter. \n"
	    , INPUTS.CUTWIN_NOSATURATE_NOBS[icut][0]
	    , INPUTS.CUTWIN_NOSATURATE_NOBS[icut][1]
	    , INPUTS.CUTWIN_NOSATURATE_FILTERS[icut]	    );
  }


  if ( INPUTS.NVAR_SIMGEN_DUMP > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(ctmp,"\t SIMGEN_DUMP file includes" );
    if ( INPUTS.APPLY_CUTWIN_OPT  == 1 ) 
      { sprintf(cptr,"%s SNe passing software cuts. \n", ctmp ); }
    else 
      { sprintf(cptr,"%s ALL SNe; use CUTMASK for software cuts.\n", ctmp); }
    
  }

  if ( INPUTS.CUTWIN_PEAKMAG[1] < 998.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t PEAKMAG(All bands) :  %.1f - %.1f\n", 
	    INPUTS.CUTWIN_PEAKMAG_ALL[0], INPUTS.CUTWIN_PEAKMAG_ALL[1] ) ;
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t PEAKMAG(any filter) < %3.1f \n", 
	  INPUTS.CUTWIN_PEAKMAG[1] ) ;
  
  for ( icut=1; icut <= INPUTS.NCUTWIN_PEAKMAG_BYFIELD; icut++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t PEAKMAG(%s): %.1f to %1.f \n", 
	    INPUTS.CUTWIN_BYFIELDLIST[icut],
	    INPUTS.CUTWIN_PEAKMAG_BYFIELD[icut][0],
	    INPUTS.CUTWIN_PEAKMAG_BYFIELD[icut][1] );
  }
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t MWEBV <= %3.1f \n", INPUTS.CUTWIN_MWEBV[1] ) ;
  
  // optional cut on time above SNRMIN
  if ( INPUTS.CUTWIN_EPOCHS_NFILT > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Time < %.1f days for SNR(%s) > %.1f \n",
	    INPUTS.CUTWIN_EPOCHS_TRANGE[1],
	    INPUTS.CUTWIN_EPOCHS_FILTERS,
	    INPUTS.CUTWIN_EPOCHS_SNRMIN );
  }

  // ------------
  *iline  = i ;

}  // end of readme_doc_CUTWIN


// ********************************************
void  readme_doc_TAKE_SPECTRUM(int *iline) {

  // Mar 14 2017: 
  // write TAKE_SPECTRUM info to readme file.
  //
  // Mar 23 2019: use GENPOLY typedefs and write WARP info

  int N = NPEREVT_TAKE_SPECTRUM ;
  int i, j , OPT_FRAME_EPOCH, OPT_TEXPOSE, IS_HOST ;
  float *ptrEP, *ptrLAM ;
  char *cptr, Tname[20], name2[20], zpolyString[100], warpString[100] ;
  char *ptrFIELD, fieldString[60] ;
  GENPOLY_DEF *GENZPOLY_SNR, *GENZPOLY_TEXPOSE, *GENLAMPOLY_WARP; 

  // ------- BEGIN -------

  if ( N == 0 ) { return ; }

  i = *iline ;
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  TAKE_SPECTRUM: \n");

  
  for(j=0; j < N; j++ ) {

    ptrEP    = INPUTS.TAKE_SPECTRUM[j].EPOCH_RANGE; 
    ptrLAM   = INPUTS.TAKE_SPECTRUM[j].SNR_LAMRANGE ;
    ptrFIELD = INPUTS.TAKE_SPECTRUM[j].FIELD ;

    GENZPOLY_SNR     = &INPUTS.TAKE_SPECTRUM[j].GENZPOLY_SNR ;
    GENZPOLY_TEXPOSE = &INPUTS.TAKE_SPECTRUM[j].GENZPOLY_TEXPOSE ;
    GENLAMPOLY_WARP  = &INPUTS.TAKE_SPECTRUM[j].GENLAMPOLY_WARP ;

    OPT_FRAME_EPOCH = INPUTS.TAKE_SPECTRUM[j].OPT_FRAME_EPOCH ;
    OPT_TEXPOSE = INPUTS.TAKE_SPECTRUM[j].OPT_TEXPOSE ;
    Tname[0] = name2[0] = zpolyString[0] = warpString[0] = IS_HOST = 0 ;
 
    if ( OPT_FRAME_EPOCH == GENFRAME_REST ) 
      { sprintf(Tname,"TREST"); }
    else if ( OPT_FRAME_EPOCH == GENFRAME_OBS )
      { sprintf(Tname,"TOBS"); }
    else if ( OPT_FRAME_EPOCH == GENFRAME_HOST )
      { sprintf(Tname,"HOST");  IS_HOST=1; }

    if ( OPT_TEXPOSE == 1 ) { 
      sprintf(name2, "TEXPOSE" ); 
      sprintf(zpolyString,"%s", GENZPOLY_TEXPOSE->STRING) ;
    }
    else { 
      // sprintf(name2, "SNR(%.0f:%.0f)", ptrLAM[0], ptrLAM[1] ) ; 
      sprintf(name2, "SNR" );
      sprintf(zpolyString,"%s", GENZPOLY_SNR->STRING) ;
    }

    if ( GENLAMPOLY_WARP->ORDER > 0 ) 
      { sprintf(warpString,"WARP:%s", GENLAMPOLY_WARP->STRING ); }

    if ( strlen(ptrFIELD) > 0 ) 
      { sprintf(fieldString,"FIELD=%s", ptrFIELD); }
    else
      { fieldString[0] = 0 ; }

    // - - - - 
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    if ( IS_HOST ) {
      sprintf(cptr,"    %s                      %s-zPOLY:%s  %s  %s\n",
	      Tname, name2, zpolyString, warpString, fieldString );
    }
    else {
      sprintf(cptr,"    %s = %5.1f to %5.1f  "
	      "  %s-zPOLY:%s  %s  %s\n",
	      Tname, ptrEP[0], ptrEP[1],
	      name2, zpolyString, warpString, fieldString );
    }

  } // end j loop over spectra

  // print optional prescale string
  STRING_DICT_DEF *DICT = &INPUTS.DICT_SPECTRUM_FIELDLIST_PRESCALE ;
  if ( DICT->N_ITEM > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"    SPECTRUM_PRESCALE: %s\n", DICT->STRING);
  }

  *iline = i;
  return ;

} // end readme_doc_TAKE_SPECTRUM


// ********************************************
void readme_doc_FUDGES(int *iline) {

  // Created May 2014 (pulled from readme_doc)
  
  int i, NLINE_FUDGE, j ;
  char *cptr ;
  char fudgeLine[10][100] ;

  // ------- BEGIN -------
  i = *iline ;

  NLINE_FUDGE = 0 ;

  // create local comment lines based on used fudges.

  if ( INPUTS.FUDGESCALE_NOISE_SKY != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for SIMLIB NOISE(SKY) : %5.2f ", 
	  INPUTS.FUDGESCALE_NOISE_SKY );
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESCALE_NOISE_READ != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for SIMLIB NOISE(CCD-read) : %5.2f ", 
	  INPUTS.FUDGESCALE_NOISE_READ );
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESCALE_PSF != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for SIMLIB PSF      : %5.2f ", 
	  INPUTS.FUDGESCALE_PSF);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESHIFT_ZPT != 0.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-shift for SIMLIB ZPT      : %5.2f ", 
	  INPUTS.FUDGESHIFT_ZPT);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGE_MAGERR != 0.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge MAGERR (add in quad)      : %6.4f ", 
	  INPUTS.FUDGE_MAGERR);
    NLINE_FUDGE++ ;
  }
  if ( INPUTS.FUDGE_ZPTERR != 0.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge ZPTERR in SIMLIB          : %6.4f ", 
	  INPUTS.FUDGE_ZPTERR);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESCALE_FLUXERR != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for FLUX-ERROR      : %5.2f ", 
	  INPUTS.FUDGESCALE_FLUXERR);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.OPT_FUDGE_SNRMAX == 1 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t FUDGE_SNRMAX    : %s  (adjust exposure time)",  
	    INPUTS.STRING_FUDGE_SNRMAX);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.OPT_FUDGE_SNRMAX == 2 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t FUDGE_SNRMAX    : %s  (adjust sky noise)",  
	    INPUTS.STRING_FUDGE_SNRMAX);
    NLINE_FUDGE++ ;
  }

  // ----------------------
  // make README comment based on whether any of the fudges are used.

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if( NLINE_FUDGE == 0 ) 
    {  sprintf(cptr,"\n  Fudges on SIMLIB Seeing Conditions: NONE. \n"); }
  else
    {  sprintf(cptr,"\n  Fudges on SIMLIB Seeing Conditions: \n"); }


  // print summary of USED fudges only.
  for(j=0; j < NLINE_FUDGE; j++ ) {
    i++ ; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "%s\n", fudgeLine[j]);
  }


  *iline = i;

  return ;

} // end readme_doc_FUDGES

// ******************************************
void readme_doc_mapFileList(int *iline) {

  // write each map file so that higher level pipelines
  // can grep out these files and compare time stamps
  // against the sim data time stamp.

  int i;
  // -------- BEGIN ---------

  i = *iline;

  i++; 
  sprintf(VERSION_INFO.README_DOC[i], "\n");

  readme_doc_mapFile(&i, "SIMLIB_FILE:", INPUTS.SIMLIB_FILE);
  readme_doc_mapFile(&i, "KCOR_FILE:",   INPUTS.KCOR_FILE);
  readme_doc_mapFile(&i, "HOSTLIB_WGTMAP_FILE:", 
		     INPUTS.HOSTLIB_WGTMAP_FILE);
  readme_doc_mapFile(&i, "HOSTLIB_ZPHOTEFF_FILE:", 
		     INPUTS.HOSTLIB_ZPHOTEFF_FILE);
  readme_doc_mapFile(&i, "HOSTLIB_SPECBASIS_FILE:", 
		     INPUTS.HOSTLIB_SPECBASIS_FILE);
  readme_doc_mapFile(&i, "HOSTLIB_SPECDATA_FILE:", 
		     INPUTS.WRONGHOST_FILE);
  readme_doc_mapFile(&i, "WRONGHOST_FILE:", 
		     INPUTS.WRONGHOST_FILE);
  readme_doc_mapFile(&i, "FLUXERRMODEL_FILE:",
		     INPUTS.FLUXERRMODEL_FILE);
  readme_doc_mapFile(&i, "NONLINEARITY_FILE:",
		     INPUTS.NONLINEARITY_FILE);
  readme_doc_mapFile(&i, "ZVARIATION_FILE:",
		     INPUT_ZVARIATION_FILE );
  readme_doc_mapFile(&i, "WEAKLENS_PROBMAP_FILE" ,
		     INPUTS.WEAKLENS_PROBMAP_FILE );
  readme_doc_mapFile(&i, "SEARCHEFF_PIPELINE_LOGIC_FILE:" ,
		     INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE);
  readme_doc_mapFile(&i, "SEARCHEFF_PIPELINE_EFF_FILE:" ,
		     INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE);
  readme_doc_mapFile(&i, "SEARCHEFF_SPEC_FILE:" ,
		     INPUTS_SEARCHEFF.USER_SPEC_FILE);
  readme_doc_mapFile(&i, "SEARCHEFF_zHOST_FILE:" ,
		     INPUTS_SEARCHEFF.USER_zHOST_FILE);

  *iline = i;

  return ;

}  // end readme_doc_mapFile_list

void readme_doc_mapFile(int *iline, char *KEY, char *FILENAME) {

  int i;
  char *cptr  ;
  char KEY_MAP[] = "MAP:" ;

  // -------------- BEGIN ------------
  i = *iline ;
  
  if ( !IGNOREFILE(FILENAME) ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s %s %s\n", KEY_MAP, KEY, FILENAME);
  }


  *iline = i;
  return ;

} // end readme_doc_mapFile

// ********************************************
void readme_doc_GENPERFECT(int *iline) {

  // April 2014

  int i, itmp ;
  char *cptr, *PARNAME ;
  double DPARVAL_ORIG, DPARVAL_USED ;
  int    IPARVAL_ORIG, IPARVAL_USED ;
  // --------- BEGIN ---------

  i = *iline ;

  if ( GENPERFECT.NVAR <= 0 ) { return ; }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n PERFECT LIGHTCURVES REQUESTED: \n");
  
  for ( itmp=1; itmp <= GENPERFECT.NVAR; itmp++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    
    PARNAME      = GENPERFECT.parnam[itmp] ;
    DPARVAL_ORIG = GENPERFECT.parval[itmp][0] ;
    DPARVAL_USED = GENPERFECT.parval[itmp][1] ;

    if ( GENPERFECT.partype[itmp] == 1 ) {
      IPARVAL_ORIG = (int)DPARVAL_ORIG ;
      IPARVAL_USED = (int)DPARVAL_USED ;
      sprintf(cptr,"\t %-28.28s : %d -> %d \n",
	      PARNAME, IPARVAL_ORIG, IPARVAL_USED );
    }
    else {
      sprintf(cptr,"\t %-28.28s : %.5f -> %.5f \n" ,
	      PARNAME, DPARVAL_ORIG, DPARVAL_USED );
    }
  }

  *iline = i;

}  // end of readme_doc_GENPERFECT

// ********************************************
void readme_doc_NON1ASED(int *iline) {

  // Created Feb 1 2017    [code moved out of readme_doc()]
  //
  // Print SEDs which have NGENTOT>0.

  int i, j, isp, index, NINDEX ;
  int NGENTOT, NGENWR, NGENUSR, CID0, CID1 ;
  float eff;
  char *cptr,  *ptrtype, ctmp[100], cline[200] ;
  char fnam[] = "readme_doc_NON1ASED" ;

  // --------------- BEGIN --------------

  i = *iline ;

  NINDEX = INPUTS.NON1ASED.NINDEX ;

  sprintf(cline,
	  " ---------------------------------------------------"
	  "-------------------------\n");

  // summarize user inputs
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n\n  NON1A USER INPUTS vs. INDEX: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr," INDEX(TYPE) " );
  for ( j = 2; j <= INPUTS.NON1ASED.NKEY; j++ )
    { sprintf(cptr,"%s %8s", cptr, INPUTS.NON1ASED.KEYLIST[j] ); }
  strcat(cptr,"\n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s", cline);

  for ( isp=1; isp <= NINDEX; isp++ ) {
    index   = INPUTS.NON1ASED.INDEX[isp];
    ptrtype = GENLC.NON1ASED.TYPE[isp];
    NGENTOT = GENLC.NON1ASED.NGENTOT[isp];
    if ( NGENTOT == 0 ) { continue ; }

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr," %3.3d (%5s) ", index, ptrtype );
    
    for ( j = 2; j <= INPUTS.NON1ASED.NKEY; j++ )
      { sprintf(cptr,"%s %8.3f", cptr, INPUTS.NON1ASED.KEYVAL[isp][j] ); }
    
    sprintf(cptr,"%s  %s\n", cptr,  INPUTS.NON1ASED.LIST_NAME[index] );
  } // end isp loop
  
  // ----------------------------------------
  // now summarize stats, effic, CID-range, peakmags ...
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  NON1A-SUMMARY vs. INDEX: \n");
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"                  NGEN    NGEN   SEARCH   \n");
  
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr," INDEX(TYPE)      written total  Effic    CID-range \n" );

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", cline);

    for ( isp=1; isp <= NINDEX; isp++ ) {
      index     = INPUTS.NON1ASED.INDEX[isp];
      NGENTOT   = GENLC.NON1ASED.NGENTOT[isp];
      NGENWR    = GENLC.NON1ASED.NGENWR[isp];
      NGENUSR   = INPUTS.NON1ASED.NGEN[isp];
      ptrtype   = GENLC.NON1ASED.TYPE[isp];
      CID0      = GENLC.NON1ASED.CIDRANGE[isp][0];
      CID1      = GENLC.NON1ASED.CIDRANGE[isp][1];

      if ( NGENTOT > 0 ) { eff = (float)NGENWR/(float)NGENTOT ; }
      else               { eff = 0.0 ; continue ; }
      
      // glue together rest-frame peakmags

      ctmp[0] = 0 ;
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr," %3.3d (%10s)  %4d  %5d   %5.3f %6d -%6d\n", 
	      index, ptrtype, NGENWR, NGENTOT, eff, CID0, CID1 );

      if ( NGENUSR != NGENWR  && INPUTS.NGEN_LC > 0 ) {
	sprintf(c1err,"NGEN[NONIA index=%d]=%d, but expected NGEN=%d",
		index, NGENWR, NGENUSR );
	errmsg(SEV_WARN, 0, fnam, c1err, ""); 
      }

    }  // end of isp loop

  *iline = i;

  return ;

} // end readme_doc_NON1ASED


// ********************************************
void readme_doc_SIMLIB(int *iline) {

  // add SIMLIB info to readme file.
  int i, j, NSKIP ;
  char *cptr, ctmp[100];

  // ------------- BEGIN ---------------

  i = *iline ;

  // - - - - - - 
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( GENLC.SIMLIB_IDLOCK >= 0 ) 
    { sprintf(ctmp,"Lock LIBID=%d", GENLC.SIMLIB_IDLOCK ); }
  else
    { sprintf(ctmp,"start at LIBID=%d", INPUTS.SIMLIB_IDSTART ); }

  if ( INPUTS.SIMLIB_NREPEAT > 1 ) // Apr 26 2017
    { sprintf(ctmp,"NREPEAT=%d", INPUTS.SIMLIB_NREPEAT ); } 

  sprintf(cptr,"\t SIMLIB filename  : %s (%s) \n", 
	  INPUTS.SIMLIB_FILE, ctmp ); 

  // - - - - - - 
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t SIMLIB SURVEY    : %s  (TELESCOPE=%s, MINOBS=%d) \n", 
	  GENLC.SURVEY_NAME, GENLC.TELESCOPE[0], INPUTS.SIMLIB_MINOBS  ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t SIMLIB UNITS     : %s for PSF,  %s for SKYSIG \n", 
	  SIMLIB_GLOBAL_HEADER.PSF_UNIT, SIMLIB_GLOBAL_HEADER.SKYSIG_UNIT  ); 

  if ( INPUTS.SIMLIB_MSKOPT != 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB MSKOPT   : %d \n", INPUTS.SIMLIB_MSKOPT );
  }

  if ( strcmp(INPUTS.SIMLIB_FIELDLIST,"ALL") != 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\t SIMLIB FIELDLIST : %s\n", INPUTS.SIMLIB_FIELDLIST); 
  }
  


  double MINSEASON = INPUTS.SIMLIB_MINSEASON;
  if ( MINSEASON > 0.01 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB min season len : %.0f days", MINSEASON);
  }
 

  int NPE_SAT  = SIMLIB_GLOBAL_HEADER.NPE_PIXEL_SATURATE ;
  int PHOTFLAG = SIMLIB_GLOBAL_HEADER.PHOTFLAG_SATURATE ;
  if ( NPE_SAT < 999999999 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB Saturation : %d pe in central pixel "
	    "-> PHOTFLAG=%d\n",  NPE_SAT, PHOTFLAG);
  }

  if ( INPUTS.FUDGEOPT_FLUXERR ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t FUDGEOPT_FLUXERR : %d \n",
	    INPUTS.FUDGEOPT_FLUXERR );
  }

  if ( INPUTS.SIMLIB_IDLOCK > 1 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB ID LOCKED to LIBID : %d \n", 
	    GENLC.SIMLIB_IDLOCK ); 
  }
  if ( INPUTS.SIMLIB_IDLOCK == 1 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB ID LOCKED to first accepted LIBID : %d \n", 
	    GENLC.SIMLIB_IDLOCK ); 
  }

  NSKIP = INPUTS.NSKIP_SIMLIB ;
  if ( NSKIP > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB IDs skipped : " );
    for ( j=0; j < NSKIP ; j++ ) {
      sprintf(ctmp," %d", INPUTS.SIMLIB_IDSKIP[j] );
      strcat(cptr, ctmp);
    }
    strcat ( cptr, "\n" );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t NEWMJD_DIF       : "
	  "%5.2f  minutes (defines trigger epoch)\n", 
	  24.*60.*INPUTS.NEWMJD_DIF ); 

  *iline = i;


  //  int    KEEP_ENTIRE_SEASON = 
  //    (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_ENTIRE_SEASON );

} // end of readme_doc_SIMLIB


// *******************************
void  readme_doc_magSmear(int *iline) {

  int i, j, ifilt, onoff;
  char *cptr, ctmp[80] ;
  char conoff[2][4] = { "OFF" , "ON" } ;
  //  char fnam[] = "readme_doc_magSmear" ;

  // ----------------- BEGIN -----------------
  i = *iline ;

  // intrinsic smearing params

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Intrinsic MAG-smearing models "
	  "(sigma clip %4.1f to %4.1f) : \n"
	  ,INPUTS.SIGMACLIP_MAGSMEAR[0]
	  ,INPUTS.SIGMACLIP_MAGSMEAR[1]
	  );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 1: Coherent MAG-smearing (GENMAG_SMEAR) : %6.3f  \n", 
	  INPUTS.GENMAG_SMEAR[0] );


  onoff = 0;
  if ( INPUTS.GENMODEL_ERRSCALE > 0.0 ) { onoff=1; }
  if ( INPUTS.NFILT_SMEAR       > 0   ) { onoff=1; }


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 2: passband MAG-smearing is %s \n", conoff[onoff] );

  if ( onoff > 0 ) {
    if ( INPUTS.GENMODEL_ERRSCALE_OPT == 1 )  
      { sprintf(ctmp,"PEAK"); }
    else           
      { sprintf(ctmp,"EPOCH-DEPENDENT"); }


    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Scale %s model errors (GENMODEL_ERRSCALE) : %6.3f  \n", 
	    ctmp, INPUTS.GENMODEL_ERRSCALE );

    if ( INPUTS.NFILT_SMEAR > 0 ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t Fixed ");
      for ( j=1; j <= INPUTS.NFILT_SMEAR; j++ ) {
	ifilt = INPUTS.IFILT_SMEAR[j]; 
	sprintf(cptr, "%s%c", cptr, FILTERSTRING[ifilt] ) ;
      }
      strcat(cptr," mag-smear:");

      for ( j=1; j <= INPUTS.NFILT_SMEAR; j++ ) {
	ifilt = INPUTS.IFILT_SMEAR[j]; 
	sprintf(cptr, "%s %4.2f", cptr, INPUTS.GENMAG_SMEAR_FILTER[ifilt] ) ;
      }
      strcat(cptr,"\n");

    }

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Correlation between models 1 & 2 : %5.2f \n",
	    INPUTS.GENMODEL_ERRSCALE_CORRELATION );

  } // NFILT_SMEAR

  // --------------------------------
  // restlambda-dependent smearing using  GENMAG_SMEAR_MODELNAME

  onoff = 0;   ctmp[0] = 0 ;
  if ( istat_genSmear() > 0 ) 
    { onoff=1;     sprintf(ctmp,"%s", INPUTS.GENMAG_SMEAR_MODELNAME); }
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 3: %s model-smear is %s  \n", 
	  ctmp, conoff[onoff] );

  // print override parameter names xxxgenmag
  for(j=0; j < NSMEARPAR_OVERRIDE; j++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\t    %s-Override :  %d values for %s\n" 
	    , ctmp
	    , GENMAG_SMEARPAR_OVERRIDE[j].NVAL
	    , GENMAG_SMEARPAR_OVERRIDE[j].NAME );
  }

  // -------------------------
  // model 4: intrinsic scatter matrix

  if ( INPUTS.NCOVMAT_SCATTER > 0 ) { onoff=1;} else { onoff=0;}

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 4: intrinsic scatter matrix is %s \n", 
	  conoff[onoff] );

  char *ptrScat; 
  for ( j=0; j < 8; j++ ) {
    ptrScat = GENLC.COVMAT_SCATTER_README[j] ;
    if ( strlen(ptrScat) > 0 ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s\n", ptrScat ); 
    }
  }

  
  // model 5: SMEAR_USRFUN
  if ( INPUTS.NPAR_GENSMEAR_USRFUN > 0 )  { onoff=1;} else { onoff=0;}
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 5: GENMAG_SMEAR_USRFUN is  %s \n", 
	  conoff[onoff] );
  
  if ( onoff == 1 ) {
    i++ ; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t FUNPAR[1-%d] = ", INPUTS.NPAR_GENSMEAR_USRFUN );
    for ( j=0; j < INPUTS.NPAR_GENSMEAR_USRFUN; j++ )  { 
      sprintf(ctmp,"%6.2f ", INPUTS.GENMAG_SMEAR_USRFUN[j] );  
      cptr = strcat(cptr,ctmp);
    }
    cptr = strcat(cptr,"\n");
  }



  *iline = i;
  return ;

} // end of readme_doc_magSmear

// ***************************************
void  readme_doc_nonLin(int *iline) {

  int i,j;
  char *cptr ;
  // ------------- BEGIN ----------------

  if ( NONLIN_README.NLINE == 0 ) { return ; }

  i = *iline ;

  for(j=0; j < NONLIN_README.NLINE; j++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s\n", NONLIN_README.LINE[j] );
  }

  *iline = i;

  return ;

} // end readme_doc_nonLin

// ***************************************
void  readme_doc_SIMSED(int *iline) {

  int i, ipar, iflag ;
  char *cptr, ctmp[80], ctmp2[80] ;

  i = *iline ;

  if ( INPUTS.OPTMASK_SIMSED == OPTMASK_GEN_SIMSED_GRIDONLY ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;  
    sprintf(cptr,"\n   %s model option is GRIDONLY & SEQUENTIAL \n", 
	    INPUTS.MODELNAME );
  } 
  else if ( INPUTS.NPAR_SIMSED > 0  ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;  
    sprintf(cptr,"\n   %s model parameters: \n", 
	    INPUTS.MODELNAME );
  }


  for ( ipar = 0; ipar < INPUTS.NPAR_SIMSED; ipar++ ) {

    iflag = INPUTS.GENFLAG_SIMSED[ipar];

    if ( (iflag & 1) == 0 ) 
      { continue ; } // baggage

    if ( (iflag & 4) > 0 ) 
      { sprintf(ctmp, "GRIDONLY" ); }
    else
      { sprintf(ctmp, "contin." ); }

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(ctmp2,"    %s(%s)", INPUTS.PARNAME_SIMSED[ipar], ctmp );
    sprintf_GENGAUSS(cptr, ctmp2, &INPUTS.GENGAUSS_SIMSED[ipar] );

  }


  *iline = i;

} // end of  readme_doc_SIMSED

// *******************************
void  readme_doc_MWXT(int *iline) {

  int i ;
  char *cptr ;

  i = *iline ;
   
  if ( INPUTS.MWEBV_FLAG ) {

    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t MilkyWay extinction  is ON  \n" );

    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t    Color law: %s  (OPT_MWCOLORLAW=%d) \n",  
	    INPUTS.STR_MWCOLORLAW, INPUTS.OPT_MWCOLORLAW );


    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    if ( INPUTS.GENRANGE_MWEBV[0] >= 0.0 ) {
      sprintf(cptr, "\t    E(B-V) randomly picked between %.3f and %.3f\n",  
	      INPUTS.GENRANGE_MWEBV[0], INPUTS.GENRANGE_MWEBV[1] );
    }
    else {
      sprintf(cptr, "\t    E(B-V): %s   (OPT_MWEBV=%d)\n",  
	      INPUTS.STR_MWEBV, INPUTS.OPT_MWEBV );
    }

    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t    sigma(MWEBV) = %.2f*MWEBV + %.2f  \n", 
	    INPUTS.MWEBV_SIGRATIO, INPUTS.MWEBV_SIG );

    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t    shift(MWEBV) = %5.3f mag \n", 
	    INPUTS.MWEBV_SHIFT );

    if ( INPUTS.APPLYFLAG_MWEBV ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ; 
      sprintf(cptr, "\t    Correct dataFile FLUXCAL for MWEBV \n");
    }
    
  }
  else {
    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t MilkyWay extinction  is OFF \n" );
  }

  *iline = i ;

} // end of readme_doc_MWXT

// *******************************
void readme_doc_filterWarn(int *iline) {

  double ZMIN, ZMAX;
  int i, ifilt, ifilt_obs, NSKIP;
  char *cptr, cfilt[2];

  i = *iline ;

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];

    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] ); 
    NSKIP = NSKIP_FILTER[ifilt_obs] ;
    ZMIN  = ZVALID_FILTER[0][ifilt_obs] ;
    ZMAX  = ZVALID_FILTER[1][ifilt_obs] ;

    if ( NSKIP > 0 && ZMIN > 100. ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t FILTER-WARNING: %s not generated for any z-range. \n", 
	      cfilt);
    }

    if ( NSKIP > 0 && ZMIN < 10. ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t FILTER-WARNING: %s generated only for %5.3f < z < %5.3f\n",
	      cfilt, ZMIN, ZMAX );       
    }

  }

  *iline = i;

} // end of readme_doc_filterWarn


// *******************************
void readme_doc_hostxt(int *iline, GEN_EXP_HALFGAUSS_DEF *GENPROFILE) {

  // add host extinction info to README 
  // Apr 4 2020: refactor and update to allow EBV or AV spec.

  char *VARNAME = GENPROFILE->NAME;
  int i;
  char *cptr;
  char cEXPON[40], cGAUSS[40] ;

  // ------------ BEGIN --------------

  i = *iline ;

  // write host extinct info

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Host Extinction Parameters: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf_GENGAUSS(cptr, "\t RV ", &INPUTS.GENGAUSS_RV);
  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Gen-Range for %s  : %4.2f to %4.2f  (model=%s) \n", 
	  VARNAME,GENPROFILE->RANGE[0], GENPROFILE->RANGE[1], INPUTS.GENSNXT );
  
 
  
  double TAU   = GENPROFILE->EXP_TAU;
  double RATIO = GENPROFILE->RATIO;
  double SIG   = GENPROFILE->SIGMA;

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( TAU < 50.0 ) {
    
    sprintf(cEXPON,"ExpNONE");
    sprintf(cGAUSS,"GaussNONE");

    if ( TAU > 1.0E-9 ) 
      { sprintf(cEXPON,"exp(-%s/%4.2f)", VARNAME, TAU );   }
    if ( INPUTS.GENGAUSIG_AV > 1.0E-9 ) 
      { sprintf(cGAUSS,"%4.2f x Gauss(%s,sig=%4.2f)", RATIO, VARNAME,SIG );  }
    
    sprintf(cptr,"\t dN/d%s = %s + %s \n", VARNAME, cEXPON, cGAUSS ); 
  }
  else  { 
    sprintf(cptr,"\t dN/d%s = flat \n", VARNAME ); 
  }
  

  *iline = i ;

} // end of readme_doc_hostxt


// **************************************
void readme_doc_FIXMAG(int *iline ) {

  int i ;
  char *cptr ;

  if ( INDEX_GENMODEL != MODEL_FIXMAG ) { return ; }

  i = *iline ;

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr, "\t %s Range: %.2f to %.2f",
	  INPUTS.GENMODEL, INPUTS.FIXMAG[0], INPUTS.FIXMAG[1] );

  *iline = i;

} // end readme_doc_FIXMAG


// **************************************
void readme_doc_SALT2params(int *iline ) {

  // Aug 11 2021: write SALT2 only if asymGauss fun is used

  int i ;
  char *cptr, string[40] ;
  char star[2];

  // ------------ BEGIN ---------

  if ( INDEX_GENMODEL != MODEL_SALT2 ) { return ; }

  i = *iline ;

  if ( INPUTS.GENGAUSS_SALT2c.USE ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf_GENGAUSS(cptr, "\t SALT2c", &INPUTS.GENGAUSS_SALT2c);
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(string, "\t Alpha" );
  sprintf_GENGAUSS(cptr, string, &INPUTS.GENGAUSS_SALT2ALPHA );
  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( INPUTS.SALT2BETA_cPOLY.ORDER >= 0 ) { 
    sprintf(cptr,"\t Beta(c) = %s \n", INPUTS.SALT2BETA_cPOLY.STRING);
  }
  else {
    sprintf(string, "\t Beta ");
    sprintf_GENGAUSS(cptr, string, &INPUTS.GENGAUSS_SALT2BETA);
  }


  *iline = i ;

  return; 

} // end of readme_doc_SALT2params
 

// **************************************
void readme_doc_GENPDF(int *iline ) {

  // Aug 11 2021: write SALT2 only if asymGauss fun is used

  int i, imap ;
  char *cptr, string[40] ;
  char star[2];

  // ------------ BEGIN ---------

  if ( INDEX_GENMODEL != MODEL_SALT2 ) { return ; }

  i = *iline ;

  for (imap=0; imap < NMAP_GENPDF; imap++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\t GENPDF: %s(%s)\n", 
	    GENPDF[imap].MAPNAME,  GENPDF[imap].GRIDMAP.VARLIST );
  }

  *iline = i ;

  return; 

} // end of readme_doc_GENPDF 

// **************************************
void sprintf_GENGAUSS(char *string, char *name, 
		      GENGAUSS_ASYM_DEF *genGauss) {

  // write genGauss info to string
  // Mar 16 2015: include new SKEW parameter.

  double s0, s1;
  char cPEAK[80], cSIGMA[80], cSKEW[80], cRANGE[80];

  // ------------------- BEGIN ---------------

  sprintf(cPEAK,"Peak=%.2f", genGauss->PEAK);

  s0 =  genGauss->SIGMA[0];
  s1 =  genGauss->SIGMA[1];
  if ( s0 < 1.0E5 || s1 < 1.0E5 ) {
    sprintf(cSIGMA,"SIG-+= %.3f,%.3f", 
	    genGauss->SIGMA[0], genGauss->SIGMA[1] );
  }
  else
    { sprintf(cSIGMA,"SIG-+ > E5,E5"); }


  sprintf(cSKEW,"SKEW=%.2f,%.2f",  
	  genGauss->SKEW[0], genGauss->SKEW[1] );
  
  sprintf(cRANGE,"BND=%.2f,%.2f", 
	  genGauss->RANGE[0], genGauss->RANGE[1] );
  
  sprintf(string,"%s: %s  %s  %s  %s\n"
	  ,name, cPEAK, cSIGMA, cSKEW, cRANGE  );

  return ;

} // sprintf_GENGAUSS




// ***************************************************
void update_accept_counters(void) {

  // Created June 2017
  // Called from main for accepted event, so here we
  // increment various counters.
  // Most code moved from main for cleanup.

  int isp ;
  //  char fnam[] = "update_accept_counters" ;

  // -------------- BEGIN ---------------

  // increment number of generated SN that are written
  NGENLC_WRITE++ ;
  NGENSPEC_WRITE += GENSPEC.NMJD_PROC ; // Jan 14 2021
  
  // increment stats based on typing method
  if ( GENLC.METHOD_TYPE == METHOD_TYPE_SPEC )    
    { GENLC.NTYPE_SPEC_CUTS++ ; }   // spec-tags after cuts
  if ( GENLC.METHOD_TYPE == METHOD_TYPE_PHOT ) { 
    GENLC.NTYPE_PHOT_CUTS++ ;    // phot-tags after cuts   
    if ( !GENLC.CORRECT_HOSTMATCH ) { GENLC.NTYPE_PHOT_WRONGHOST++; }
  }
  // check screen update
  if ( INPUTS.NGEN_LC > 0 )  { screen_update(); }
      

  // increment NGENWR per non1a index
  if ( INPUTS.NON1ASED.NINDEX > 0  ) {
    isp = GENLC.NON1ASED.ISPARSE;
    GENLC.NON1ASED.NGENWR[isp]++; 
  }

  // Jun 2017: check SUBSAMPLE_MARK option
  double XN0, XN1;
  int NTMP = INPUTS.NSUBSAMPLE_MARK ;
  if ( NTMP > 1 ) {
    XN0   = (double)NTMP ; 
    XN1   = (double)NGENLC_WRITE ;
    GENLC.SUBSAMPLE_INDEX = (int)fmod( XN1, XN0 ) ;
  }
  else {
    GENLC.SUBSAMPLE_INDEX = -9 ;
  }

  return ;

} // end update_accept_counters

// ***************************************************
void init_simFiles(SIMFILE_AUX_DEF *SIMFILE_AUX) {

  // May 28, 2009 R.Kessler
  // Init auxiliary files create by sim program.
  // Does NOT include special dump files where
  // the sim quits after the dump.
  //
  // Jun 2011: check option to call WR_SNFITSIO_INIT
  // Feb 12, 2014: always call snlc_to_SNDATA(1) instead of only
  //               for FITS format.
  // Feb 06, 2021: Remove .IGNORE file (no longer required)
  // Oct 14 2021: set spectra bit of INPUTS.WRITE_MASK 

  int i, isys ;
  char headFile[MXPATHLEN];
  char cmd[2*MXPATHLEN], prefix[2*MXPATHLEN];
  char fnam[] = "init_simFiles" ;

  // ------------ BEGIN -------------

  // always construct readme lines
  readme_doc(1);    

  // init DUMP file regardless of SNDATA file status

  if ( INPUTS.FORMAT_MASK <= 0 ) {
    sprintf(SIMFILE_AUX->DUMP,  "%s.DUMP",  INPUTS.GENVERSION );
    wr_SIMGEN_DUMP(1,SIMFILE_AUX);  // always make DUMP file if requested
    return ;
  }

  // clear out old GENVERSION files; 2nd arg is PROMPT flag
  clr_VERSION(INPUTS.GENVERSION,INPUTS.CLEARPROMPT);

  // create new subdir for simulated SNDATA files.
  // Note that -p is not used to avoid bad behavior.
  sprintf(cmd,"mkdir -m g+wr %s", PATH_SNDATA_SIM );
  isys = system(cmd);


  // create full names for auxilliary files,
  // whether they are used or not.
  sprintf(prefix,"%s/%s", PATH_SNDATA_SIM, INPUTS.GENVERSION );

  // mandatory
  sprintf(SIMFILE_AUX->LIST,       "%s.LIST",        prefix );
  sprintf(SIMFILE_AUX->README,     "%s.README",      prefix );
  // xxx  sprintf(SIMFILE_AUX->IGNORE,     "%s.IGNORE",      prefix );


  // optional
  sprintf(SIMFILE_AUX->DUMP,       "%s.DUMP",        prefix );
  sprintf(SIMFILE_AUX->ZVAR,       "%s.ZVARIATION",  prefix );
  sprintf(SIMFILE_AUX->GRIDGEN,    "%s.GRID",        prefix );

  // Aug 10 2020: for batch mode, write YAML file locally so that
  //              it is easily found by batch script.
  sprintf(SIMFILE_AUX->YAML,  "%s.YAML",  INPUTS.GENVERSION ); // Aug 10, 2020


  // create mandatory files.
  SIMFILE_AUX->FP_LIST   = fopen(SIMFILE_AUX->LIST,   "wt") ;  
  SIMFILE_AUX->FP_README = fopen(SIMFILE_AUX->README, "wt") ;  

  // dump out the README file
  for ( i = 1; i <= VERSION_INFO.NLINE_README_INIT; i++ )
    { fprintf(SIMFILE_AUX->FP_README, "%s", VERSION_INFO.README_DOC[i] ); }

  fflush(SIMFILE_AUX->FP_README);

  // if FITRES DUMP-file is requested, open and init header
  wr_SIMGEN_DUMP(1,SIMFILE_AUX);
  
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
#ifdef MODELGRID_GEN
    init_GRIDsource(1); 
    wr_GRIDfile(1,SIMFILE_AUX->GRIDGEN);
#endif
  }

  snlc_to_SNDATA(1) ;  // 1 => load header only

  // check option for fits format (Jun 2011)
  if ( WRFLAG_FITS ) { 

    // Oct 14 2021 - check to set write-mask bit for spectra
    if ( SPECTROGRAPH_USEFLAG ) {
      int OPTMASK = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK;
      int REFAC  = (OPTMASK & SPECTROGRAPH_OPTMASK_FITS_REFAC);
      int LEGACY = (OPTMASK & SPECTROGRAPH_OPTMASK_FITS_LEGACY); // default
      if ( REFAC ) 
	{ INPUTS.WRITE_MASK += WRITE_MASK_SPECTRA ; }
      else 
	{ INPUTS.WRITE_MASK += WRITE_MASK_SPECTRA_LEGACY ; }
    } 

    // abort of any text-option is defined along with fits format
    if ( WRFLAG_TEXT  ) {
      sprintf(c1err, "Cannot mix TEXT and FITS format ; FORMAT_MASK=%d",
	      INPUTS.FORMAT_MASK );
      sprintf(c2err,"WRFLAG[TEXT,FITS] = %d, %d",
	      WRFLAG_TEXT, WRFLAG_FITS );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    WR_SNFITSIO_INIT(PATH_SNDATA_SIM
		     , INPUTS.GENVERSION
		     , INPUTS.GENPREFIX
		     , INPUTS.WRITE_MASK
		     , INPUTS.NSUBSAMPLE_MARK
		     , headFile); // <== return arg for list file

    fprintf(SIMFILE_AUX->FP_LIST,"%s", headFile);
  }

  // write filter responses for non-SNANA programs
  if ( WRFLAG_FILTERS ) 
    { wr_SIMGEN_FILTERS(SIMFILE_AUX->PATH_FILTERS); }
 
  return ;

} // end of init_simFiles

// ***********************************
void update_simFiles(SIMFILE_AUX_DEF *SIMFILE_AUX) {

  // May 28, 2009: control output from here
  // May 24, 2011: call snfitsio_update()
  // Jan 23, 2014: remove call to append_SNDATA_MODEL();
  //
  // May 27, 2019: 
  //  + call wr_SIMGEN_DUMP after snlc_to_SNDATA to allow for
  //    things like PEAKMJD_SMEAR

  int  CID    ;
  char fnam[] = "update_simFiles";

  // ------------ BEGIN -------------

#ifdef MODELGRID_GEN
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
    update_GRIDarrays();
    wr_GRIDfile(2,"");
    return ;
  }
#endif

  if ( WRFLAG_CIDRAN == 0 ) 
    { CID = GENLC.CID ; }
  else
    { CID = GENLC.CIDRAN ; }

  if ( CID > MXCID_SIM  ) {
    sprintf(c1err,"CID=%d exceeds MXCID_SIM=%d", CID, MXCID_SIM );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }


  // init SNDATA strucure
  init_SNDATA_EVENT() ; 

  // load SNDATA structure
  snlc_to_SNDATA(0) ;

  // always fill DUMP file even if SNDATA files are not written
  wr_SIMGEN_DUMP(2,SIMFILE_AUX);

  if ( INPUTS.FORMAT_MASK <= 0 ) { return ; }

  if ( WRFLAG_FITS ) { 
    WR_SNFITSIO_UPDATE(); 
    return ;
  }

  WR_SNTEXTIO_DATAFILE(SNDATA.SNFILE_OUTPUT);

  // below are the text-output options

  // update LIST file
  fprintf(SIMFILE_AUX->FP_LIST, "%s\n", SNDATA.snfile_output);

  return ;

} // end of update_simFiles


// ***********************************
void end_simFiles(SIMFILE_AUX_DEF *SIMFILE_AUX) {

  // May 28, 2009
  // close out auxiliary simfiles, and leave 1-line 
  // summary to stdout for each aux file.
  //
  // Oct 16, 2010: print GRIDGEN summary if used
  // Jun 09, 2012:  replace GRIDGEN stuff call to wr_GRIDfile(3);
  // May 14, 2019: add call to wr_SIMGEN_DUMP(3,SIMFILE_AUX);
  int i, N1, N2 ;
  
  // ------------ BEGIN -------------

  iter_summary_genPDF();

  // fill post-sim part of readme
  readme_doc(2); 

  // always dump entire readme contents to screen
  print_banner("DUMP README CONTENTS TO SCREEN\n");
  for ( i = 1; i<= VERSION_INFO.NLINE_README; i++ )
    printf("%s", VERSION_INFO.README_DOC[i] );

  // ==========================================
  // continue only if SNDATA files are written for each SN

  if ( INPUTS.FORMAT_MASK <= 0 ) return ;

  // ==========================================

  printf("\n  AUXILIARY FILES: \n");
  printf("  %s \n", SIMFILE_AUX->LIST );
  printf("  %s \n", SIMFILE_AUX->README );

  if ( INPUTS.WRFLAG_YAML_FILE > 0 ) 
    { printf("  %s \n", SIMFILE_AUX->YAML ); }  // for batch mode, Aug 10 2020

  if ( WRFLAG_FILTERS ) 
    { printf("  %s \n", SIMFILE_AUX->PATH_FILTERS ); }  // it's a subdir

  fflush(stdout);

  // dump post-sim part of readme to README file.
  N1 = VERSION_INFO.NLINE_README_INIT + 1 ; 
  N2 = VERSION_INFO.NLINE_README ; 
  for ( i = N1; i <= N2; i++ )
    {  fprintf(SIMFILE_AUX->FP_README, "%s", VERSION_INFO.README_DOC[i] ); }

  // close files.
  fclose(SIMFILE_AUX->FP_LIST);
  fclose(SIMFILE_AUX->FP_README);

  // check optional auxiliary files.

  // close out SIMGEN_DUMP file if it exists
  if ( INPUTS.NVAR_SIMGEN_DUMP > 0 ) {
    wr_SIMGEN_DUMP(3,SIMFILE_AUX);
  }


  // copy ZVARATION file to SIM/[VERSION]
  if ( USE_ZVAR_FILE ) {
    cp_zvariation(SIMFILE_AUX->ZVAR);  
    printf("  %s\n", SIMFILE_AUX->ZVAR);
  }

  // Aug 10 2020: in batch mode, write few stats to YAML formatted file
  if ( INPUTS.WRFLAG_YAML_FILE > 0 ) {  wr_SIMGEN_YAML(SIMFILE_AUX); } 

#ifdef MODELGRID_GEN
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID ) {
    printf("  %s\n", SIMFILE_AUX->GRIDGEN );
    wr_GRIDfile(3,"");
  }
#endif

  fflush(stdout);
  if ( WRFLAG_FITS) { WR_SNFITSIO_END(); }

} // end of end_simFiles


// ===========================
void set_screen_update(int NGEN) {

  if ( NGEN < 300 ) 
    { INPUTS.NGEN_SCREEN_UPDATE = 1 ; }
  else if ( NGEN < 1000 ) 
    { INPUTS.NGEN_SCREEN_UPDATE = 10 ; }
  else if ( NGEN < 20000 ) 
    { INPUTS.NGEN_SCREEN_UPDATE = 100 ; }
  else
    { INPUTS.NGEN_SCREEN_UPDATE = 500 ; }

}

// ===========================
void screen_update(void) {

  int CID, NGEN;
  bool QUIT_NOREWIND = (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_QUIT_NOREWIND) >0;
  int  NLIBID        = SIMLIB_GLOBAL_HEADER.NLIBID_VALID ;
  char ctmp[40];

  // July 26 2021: print GenRate every 100 sec.

  // ---------- BEGIN --------------

  if ( !WRFLAG_CIDRAN  ) 
    { CID = GENLC.CID ; }
  else
    { CID = GENLC.CIDRAN ; }

  // - - - - -

  if ( INPUTS.NGEN_LC > 0 ) {

    if ( LUPDGEN(NGENLC_WRITE)  ) {
      printf("\t Finished writing %6d of %d (CID=%8d, NEP=%3d) \n", 
	     NGENLC_WRITE, INPUTS.NGEN_LC, CID, GENLC.NEPOCH );
      fflush(stdout);
    }

  } else {

    if ( LUPDGEN(NGENLC_TOT) ) {

      NGEN = INPUTS.NGENTOT_LC ; // expected number to generate at end

      if ( QUIT_NOREWIND && NLIBID > 0 ) {
	// user set option to make 1 and only 1 pass thru SIMLIB
	printf("\t Finished generating %8d of %d valid LIBIDs \n", 
	       NGENLC_TOT, NLIBID );
      }
      else {
	// most likely to end up here

	ctmp[0] = 0 ;
	if ( INPUTS.NGEN_SCREEN_UPDATE > 50 && NGENLC_TOT > 100 ) {
	  // time function has 1 sec resolution, so wait at least
	  // 100 sec to update GenRate so we get 1% precision.
	  time_t t_now         = time(NULL);
	  time_t t_update_last = TIMERS.t_update_last ;
	  if ( t_now - t_update_last > 100 ) {
	    int    NGENTOT_LAST  = TIMERS.NGENTOT_LAST ;
	    double XNDIF         = (double)(NGENLC_TOT - NGENTOT_LAST) ;
	    double TDIF          = (double)(t_now - t_update_last) ;
	    int GenRate          = (int)( XNDIF/TDIF ) ;
	    sprintf(ctmp,"GenRate=%d/sec", GenRate);
	    TIMERS.t_update_last = t_now ;
	    TIMERS.NGENTOT_LAST  = NGENLC_TOT ;
	  }
	}

	printf("\t Finished generating %8d of %d (CID=%d) %s\n", 
	       NGENLC_TOT, INPUTS.NGENTOT_LC, CID, ctmp );
      }

      fflush(stdout);
    } // end LUPDGEN

  }

  return ;

} // end of screen_update


// ============================================
void prioritize_genPDF_ASYMGAUSS(void) {

  // Created Aug 11 2021
  // When GENPDF_FILE arg is read, the map contents are not known
  // until it is parsed. Here (after parsing) we check priority 
  // when both genPDF and GENGAUSS are defined. Command line arg 
  // has priority over input file; lower priority option is
  // disabled as if it were never read.
  // Abort if both options have same priority.
  // Priority is based on KEYSOURCE = 1(FILE) or 2(command line arg)
  //
  // Finally, check for reasons to set DOGEN_SHAPE[COLOR] to false.

  int  FUNTYPE_ASYMGAUSS     = 1;
  int  FUNTYPE_EXPHALFGAUSS  = 2;
  char FUNNAME[3][20] = { "", "AsymGauss", "Exp+HalfGauss" };

  int  IDMAP_LIST[10], FUNTYPE_LIST[10], FUNTYPE, IDMAP, i, NCHECK = 0;
  bool IS_LOGPAR, USE_GENPDF, USE_FUN;
  int  KEYSOURCE_FUN ; // asymGauss or exp_halfGauss function
  char PARNAME_LIST[10][20], *PARNAME;
  GENGAUSS_ASYM_DEF     *ptr_ASYMGAUSS_LIST[10];
  GEN_EXP_HALFGAUSS_DEF *ptr_EXPHALFGAUSS_LIST[10];
  char KEYSOURCE_STR[3][20] = { "", "INPUT-FILE", "COMMAND-LINE" };

  char fnam[] = "prioritize_genPDF_ASYMGAUSS";

  // ------------- BEGIN -------------
  
  // store asymGauss functions
  sprintf(PARNAME_LIST[NCHECK],"SALT2x1");
  IDMAP_LIST[NCHECK]        = IDMAP_GENPDF(PARNAME_LIST[NCHECK], &IS_LOGPAR);
  ptr_ASYMGAUSS_LIST[NCHECK] = &INPUTS.GENGAUSS_SALT2x1 ;
  FUNTYPE_LIST[NCHECK] = FUNTYPE_ASYMGAUSS;
  NCHECK++ ;

  sprintf(PARNAME_LIST[NCHECK],"SALT2c");
  IDMAP_LIST[NCHECK]        = IDMAP_GENPDF(PARNAME_LIST[NCHECK], &IS_LOGPAR);
  ptr_ASYMGAUSS_LIST[NCHECK] = &INPUTS.GENGAUSS_SALT2c ;
  FUNTYPE_LIST[NCHECK] = FUNTYPE_ASYMGAUSS;
  NCHECK++ ;

  sprintf(PARNAME_LIST[NCHECK],"RV");
  IDMAP_LIST[NCHECK]        = IDMAP_GENPDF(PARNAME_LIST[NCHECK], &IS_LOGPAR);
  ptr_ASYMGAUSS_LIST[NCHECK] = &INPUTS.GENGAUSS_RV ;
  FUNTYPE_LIST[NCHECK] = FUNTYPE_ASYMGAUSS;
  NCHECK++ ;

  // store expHalfGauss functions
  sprintf(PARNAME_LIST[NCHECK],"EBV");
  IDMAP_LIST[NCHECK]        = IDMAP_GENPDF(PARNAME_LIST[NCHECK], &IS_LOGPAR);
  ptr_EXPHALFGAUSS_LIST[NCHECK] = &INPUTS.GENPROFILE_EBV_HOST ;
  FUNTYPE_LIST[NCHECK] = FUNTYPE_EXPHALFGAUSS;
  NCHECK++ ;

  sprintf(PARNAME_LIST[NCHECK],"EBV_HOST");
  IDMAP_LIST[NCHECK]        = IDMAP_GENPDF(PARNAME_LIST[NCHECK], &IS_LOGPAR);
  ptr_EXPHALFGAUSS_LIST[NCHECK] = &INPUTS.GENPROFILE_EBV_HOST ;
  FUNTYPE_LIST[NCHECK] = FUNTYPE_EXPHALFGAUSS;
  NCHECK++ ;

  sprintf(PARNAME_LIST[NCHECK],"AV");
  IDMAP_LIST[NCHECK]        = IDMAP_GENPDF(PARNAME_LIST[NCHECK], &IS_LOGPAR);
  ptr_EXPHALFGAUSS_LIST[NCHECK] = &INPUTS.GENPROFILE_EBV_HOST ;
  FUNTYPE_LIST[NCHECK] = FUNTYPE_EXPHALFGAUSS;
  NCHECK++ ;

  // - - - - -
  for(i=0; i < NCHECK; i++ ) {
    PARNAME      = PARNAME_LIST[i];
    USE_GENPDF   = ( IDMAP_LIST[i] >= 0 );
    FUNTYPE      = FUNTYPE_LIST[i];

    if ( FUNTYPE == FUNTYPE_ASYMGAUSS ) {
      USE_FUN       = ptr_ASYMGAUSS_LIST[i]->USE ; 
      KEYSOURCE_FUN = ptr_ASYMGAUSS_LIST[i]->KEYSOURCE ;
    }
    else {
      USE_FUN       = ptr_EXPHALFGAUSS_LIST[i]->USE ; 
      KEYSOURCE_FUN = ptr_EXPHALFGAUSS_LIST[i]->KEYSOURCE ;
    }
   
    if ( USE_GENPDF && USE_FUN ) {

      // both analytic-function and GENPDF-map are specified;
      // either abort, or implement priority.

      if ( KEYSOURCE_FUN == KEYSOURCE_GENPDF ) {
	sprintf(c1err,"Ambiguous method to generate '%s' ; ", PARNAME);
	sprintf(c2err,"GENPDF and %s are both from %s",
		FUNNAME[FUNTYPE], KEYSOURCE_STR[KEYSOURCE_GENPDF] );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	
      }
      else if ( KEYSOURCE_FUN < KEYSOURCE_GENPDF ) {
	// disable function as if it were never read
	printf("\t GENPDF overrides %s for %s\n", FUNNAME[FUNTYPE],PARNAME);
	if ( FUNTYPE == FUNTYPE_ASYMGAUSS )  {
	  init_GENGAUSS_ASYM(ptr_ASYMGAUSS_LIST[i], 0.0 ); 
	}
	else {
	  init_GEN_EXP_HALFGAUSS(ptr_EXPHALFGAUSS_LIST[i], 0.0 ); 
	}
      }
      else if ( KEYSOURCE_FUN > KEYSOURCE_GENPDF ) {
	// disable GENPDF map by changing VARNAME
	// so that IDMAP_GENPDF cannot find a match.
	printf("\t %s overrides GENPDF for %s\n", 
	       FUNNAME[FUNTYPE], PARNAME); 
	sprintf(GENPDF[IDMAP].VARNAMES[0],"DISABLE-%s", PARNAME);
      }
      
    } // end both methods

  } // end i loop over possible variables to generate

  fflush(stdout);

  // - - - - - - - -
  
  bool GETPAR_HOSTLIB     = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USESNPAR) ;

  bool GETc_ASYMGAUSS   = INPUTS.GENGAUSS_SALT2c.USE;
  bool GETc_GENPDF      = (IDMAP_GENPDF("SALT2c", &IS_LOGPAR) >= 0);
  bool GETc             = GETc_ASYMGAUSS || GETc_GENPDF ;
  bool SKIPc            = (GETPAR_HOSTLIB && !GETc );
  if ( SKIPc ) { INPUTS.DOGEN_COLOR = false; }

  bool ISMODEL_SIMSED    = ( INDEX_GENMODEL == MODEL_SIMSED );
  bool ISMODEL_SIMLIB    = ( INDEX_GENMODEL == MODEL_SIMLIB );
  bool ISMODEL_NON1A     = ( INPUTS.NON1A_MODELFLAG > 0 );
  bool ISMODEL_LCLIB     = ( INDEX_GENMODEL == MODEL_LCLIB ) ;
  bool NOSHAPE = ( ISMODEL_SIMSED || ISMODEL_SIMLIB|| ISMODEL_NON1A || ISMODEL_LCLIB || IS_PySEDMODEL );

  bool GETx_ASYMGAUSS   = (INPUTS.GENGAUSS_SHAPEPAR.USE);
  bool GETx_GENPDF      = (IDMAP_GENPDF("SALT2x1", &IS_LOGPAR) >= 0);
  bool GETx             = GETx_ASYMGAUSS || GETx_GENPDF ;
  bool SKIPx            = NOSHAPE || ( GETPAR_HOSTLIB && !GETx ) ;
  if ( SKIPx ) { INPUTS.DOGEN_SHAPE = false; }

  return;

} // end prioritize_genPDF_ASYMGAUSS

// ***********************************
void DASHBOARD_DRIVER(void) {

  // Created July 2019
  // Dash-board dump of files and key options.
  // Not too detailed, but enough to get a global overview.

  char fileName_orig[MXPATHLEN];  // fileName before ENVreplace
  char fnam[] = "DASHBOARD_DRIVER";

  // ---------------- BEGIN ----------------

  if ( !INPUTS.DASHBOARD_DUMPFLAG ) { return ; }

  print_banner(fnam);


  ENVrestore(INPUTS.GENMODEL,fileName_orig);
  printf("GENMODEL:        %s \n", fileName_orig);

  // ------- SIMLIB -------

  INPUTS.SIMLIB_DUMP = 1 ;  SIMLIB_DUMP_DRIVER();
  ENVrestore(INPUTS.SIMLIB_FILE,fileName_orig);
  printf("SIMLIB_FILE:            %s\n", fileName_orig);
  printf("\t %d observation sequences  %d < MJD < %d \n",  
	 NREAD_SIMLIB, 
	 (int)SIMLIB_DUMP_AVGALL.MJDMIN, (int)SIMLIB_DUMP_AVGALL.MJDMAX);

  ENVrestore(INPUTS.KCOR_FILE,fileName_orig);
  printf("KCOR_FILE:              %s\n", fileName_orig);

  ENVrestore(INPUTS.HOSTLIB_FILE,fileName_orig);
  printf("HOSTLIB_FILE:           %s\n", fileName_orig);
  if ( HOSTLIB.NGAL_READ > 0 ) {
    printf("\t %d galaxies, %.2f < z < %.2f \n", 
	   HOSTLIB.NGAL_READ, HOSTLIB.ZMIN, HOSTLIB.ZMAX );
  }

  ENVrestore(INPUTS.HOSTLIB_WGTMAP_FILE,fileName_orig);
  printf("HOSTLIB_WGTMAP_FILE:    %s\n", fileName_orig);

  ENVrestore(INPUTS.HOSTLIB_ZPHOTEFF_FILE,fileName_orig);
  printf("HOSTLIB_ZPHOTEFF_FILE:  %s\n", fileName_orig);

  printf("HOSTLIB_SPECBASIS_FILE: %s\n", INPUTS.HOSTLIB_SPECBASIS_FILE);
  printf("HOSTLIB_SPECDATA_FILE:  %s\n", INPUTS.HOSTLIB_SPECDATA_FILE);
  printf("WRONGHOST_FILE:         %s\n", INPUTS.WRONGHOST_FILE);
  printf("FLUXERRMODEL_FILE:      %s\n", INPUTS.FLUXERRMODEL_FILE);

  ENVrestore(INPUTS.NONLINEARITY_FILE,fileName_orig);
  printf("NONLINEARITY_FILE:      %s\n", fileName_orig );

  ENVrestore(INPUT_ZVARIATION_FILE,fileName_orig);
  printf("ZVARIATION_FILE:        %s\n", fileName_orig );

  ENVrestore(INPUTS.WEAKLENS_PROBMAP_FILE,fileName_orig);
  printf("WEAKLENS_PROBMAP_FILE:  %s\n", fileName_orig);

  ENVrestore(INPUTS.STRONGLENS_FILE,fileName_orig);
  printf("STRONGLENS_FILE:        %s\n", fileName_orig);

  ENVrestore(INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE,fileName_orig);
  printf("SEARCHEFF_PIPELINE_LOGIC_FILE: %s\n", fileName_orig );
  printf("\t Logic: %s  (epoch-sep > %.3f days)\n", 
	 SEARCHEFF_LOGIC.INPUT_STRING, INPUTS_SEARCHEFF.TIME_SINGLE_DETECT);

  ENVrestore(INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE,fileName_orig);
  printf("SEARCHEFF_PIPELINE_EFF_FILE:   %s\n", fileName_orig);
  if ( INPUTS_SEARCHEFF.NMAP_DETECT > 0 ) {
    printf("\t NMAP_EFFDETECT=%d, NMAP_PHOTPROB=%d \n", 
	   INPUTS_SEARCHEFF.NMAP_DETECT, INPUTS_SEARCHEFF.NMAP_PHOTPROB);
  }

  ENVrestore(INPUTS_SEARCHEFF.USER_SPEC_FILE,fileName_orig);
  printf("SEARCHEFF_SPEC_FILE:    %s\n", fileName_orig );  
  if ( INPUTS_SEARCHEFF.NMAP_SPEC > 0 ) 
    { printf("\t NMAP_SPECEFF = %d \n", INPUTS_SEARCHEFF.NMAP_SPEC); }

  ENVrestore(INPUTS_SEARCHEFF.USER_zHOST_FILE,fileName_orig);
  printf("SEARCHEFF_zHOST_FILE:   %s\n", fileName_orig );  
  if ( INPUTS_SEARCHEFF.NMAP_zHOST > 0 ) 
    { printf("\t NMAP_zHOST = %d \n", INPUTS_SEARCHEFF.NMAP_zHOST); }

  happyend();
  return;

} // end DASHBOARD_DRIVER

// ***********************************
void SIMLIB_DUMP_DRIVER(void) {

  // Dump summary information for SIMLIB 
  //
  // History:
  //
  // May 2 2017;
  //   +  add RA & DEC to fpdmp1 (dump for every MJD)
  //   +  bail when NREAD > NGENTOT_LC
  //
  // Aug 4 2017: 
  //  + refactor so that SIMLIB_DUMP[MXREAD_SIMLIB] ->  SIMLIB_DUMP
  //    and more than 100 MB of memory is eliminated.
  //  + use new utils zero_SIMLIB_DUMP and update_SIMLIB_DUMP_AVGALL
  //
  // Sep 22 2017:  
  //  + set INPUTS.SIMLIB_NREPEAT=1 to override user value
  //
  //
  // Oct 16 2017:  MJDGAP_IGNORE-> 50 (was 100)
  // Oct 23 2017:  add MJD_DIF to SEQuence file (time since last obs)
  //
  // Jun 18 2018: compute avergae weight for LCLIB model
  //
  // July 16 2018: include BAND: key in output table for parsing
  //
  // Jul 02 2021: malloc large MJDLIST arrays to avoid stack issues.
  // ----------------------------------

#define MXSIMLIB_DUMP_STDOUT 50   // max simlib entries to screen-dump
#define ZPTERR_MAX  0.05      // exclude entries with this much ZPT-variation

  int QUIET =  INPUTS.DASHBOARD_DUMPFLAG;

  int 
    ID, IDLAST, NREAD, LDMP_LOCAL
    ,LDMP_SEQ_TEXT, LDMP_OBS_TEXT, LDMP_ROOT
    ,ifilt, ifilt_obs, iep, icut, Ncut, NVAR, NROW_MJD, NROW
    ;

  char  cfilt[2];

  FILE *fpdmp0, *fpdmp1 ;
  
  double
    MJD, MJD_LAST, GAPMAX, GAPAVG, MJDWIN, FRAC, *ptrmjd
    ,ZPT_pe, PSF, SKYSIG_ADU, SKYSIG_pe
    ,ZPTERR, M5SIG
    ,ZPT_SIMLIB, PSF_SIMLIB, SNR_maglimit
    ,FSKY_pe, SKYMAG, RA, DEC, MWEBV
    ,XNobs, TMP, TMP0,  TMP1, wgt_LCLIB, wgtsum_LCLIB=0.0
    ,GLOBAL_RANGE_RA[2]
    ,GLOBAL_RANGE_DEC[2]
    ,GAIN_SIMLIB, PIXSIZE_SIMLIB
    ;


  float  MJDMIN4, MJDMAX4, RA4, DEC4 ;

  // SNcadenceFoM args

  double 
     *MJDLIST_ALL
    ,*MJDLIST[MXFILTINDX]    // MJD list per band
    ,*M5SIGLIST[MXFILTINDX]   // idem for M5sigma    
    ;

  double  MJDGAP_IGNORE = 50.0 ; // ignore Gaps (days) longer than this
  int Nobs;
  char ctmp[40], FIELDNAME[60] ;

#define MXSIMLIB_CUTCHECK 10

  int    NCUT_SIMLIB;
  double SIMLIB_GENRANGE[MXSIMLIB_CUTCHECK][2] ;
  int    NSIMLIB_CUTFAIL[MXSIMLIB_CUTCHECK] ;
  char   SIMLIB_CUTVARNAME[MXSIMLIB_CUTCHECK][20] ;
  float *PTR_SIMLIB_CUTVAR[MXSIMLIB_CUTCHECK][2] ;

  char  fnam[] = "SIMLIB_DUMP_DRIVER" ;

  // ------------ BEGIN  ----------


  LDMP_LOCAL = 1 * (1-QUIET) ;

  if ( INPUTS.SIMLIB_DUMP < 0 ) { return ; }

  // check dump options (Aug 8 2013)

  if ( INPUTS.SIMLIB_DUMP == 0 ) {
    sprintf(c1err,"SIMLIB_DUMP = 0 is no longer valid.");
    sprintf(c2err,"Must set bit0,1 to dump per LIBID, per MJD");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  LDMP_SEQ_TEXT =  LDMP_OBS_TEXT = LDMP_ROOT = 0 ; 
  if ( !QUIET ) { 
    if ( INPUTS.SIMLIB_DUMP & SIMLIB_DUMPMASK_SEQ ) { LDMP_SEQ_TEXT = 1; }
    if ( INPUTS.SIMLIB_DUMP & SIMLIB_DUMPMASK_OBS ) { LDMP_OBS_TEXT = 1; }
  }

#ifdef USE_ROOT
  if ( INPUTS.SIMLIB_DUMP & 4 ) { LDMP_ROOT  = 1; }
#endif

  if ( INPUTS.DASHBOARD_DUMPFLAG == 0 ) { print_banner("SIMLIB_DUMP"); }

  /*
  printf("\t size of SIMLIB_DUMP struct: %6.2f MB \n",
	sizeof(SIMLIB_DUMP)/1.0E6 );
  */

  // store SIMLIB cuts to check

  init_GENLC();

  icut=0;
  icut++;
  SIMLIB_GENRANGE[icut][0] = INPUTS.GENRANGE_RA[0] ;
  SIMLIB_GENRANGE[icut][1] = INPUTS.GENRANGE_RA[1] ;
  sprintf(SIMLIB_CUTVARNAME[icut],"RA");
  PTR_SIMLIB_CUTVAR[icut][0] = &RA4 ;
  PTR_SIMLIB_CUTVAR[icut][1] = &RA4 ;

  icut++;
  SIMLIB_GENRANGE[icut][0] = INPUTS.GENRANGE_DEC[0] ;
  SIMLIB_GENRANGE[icut][1] = INPUTS.GENRANGE_DEC[1] ;
  sprintf(SIMLIB_CUTVARNAME[icut],"DEC");
  PTR_SIMLIB_CUTVAR[icut][0] = &DEC4 ;
  PTR_SIMLIB_CUTVAR[icut][1] = &DEC4 ;


  icut++;
  SIMLIB_GENRANGE[icut][0] = INPUTS.GENRANGE_PEAKMJD[0] - 50. ;
  SIMLIB_GENRANGE[icut][1] = INPUTS.GENRANGE_PEAKMJD[1] + 50. ;
  sprintf(SIMLIB_CUTVARNAME[icut],"PEAKMJD");
  PTR_SIMLIB_CUTVAR[icut][0] = &MJDMIN4 ;
  PTR_SIMLIB_CUTVAR[icut][1] = &MJDMAX4 ;

  NCUT_SIMLIB = icut;

  for ( icut=0; icut < MXSIMLIB_CUTCHECK; icut++ )
    { NSIMLIB_CUTFAIL[icut] = 0 ; }


  // open cuts to see entire library
  INPUTS.GENRANGE_TREST[0] = -99999. ;
  INPUTS.GENRANGE_TREST[1] = +99999. ;

  INPUTS.GENRANGE_RA[0] = -99999. ;
  INPUTS.GENRANGE_RA[1] = +99999. ;

  INPUTS.GENRANGE_DEC[0] = -99999. ;
  INPUTS.GENRANGE_DEC[1] = +99999. ;

  INPUTS.SIMLIB_NREPEAT = 1 ;

  // init local counters
  NREAD_SIMLIB = 0;
  ID = NREAD = 0; IDLAST = -1 ;

  // init global min/max values
  GLOBAL_RANGE_RA[0]   = +9999999. ;  
  GLOBAL_RANGE_DEC[0]  = +9999999. ;  
  GLOBAL_RANGE_RA[1]   = -9999999. ;  
  GLOBAL_RANGE_DEC[1]  = -9999999. ;  
  
  zero_SIMLIB_DUMP(&SIMLIB_DUMP_AVGALL);
  zero_SIMLIB_DUMP(&SIMLIB_DUMP_NAVGALL);

  NROW_MJD = NROW = 0 ;

  fpdmp0 = fpdmp1 = NULL ;

  // malloc local arrays
  int MEMD = sizeof(double);
  MJDLIST_ALL = (double*) malloc( MXEPSIM * MEMD);
  for (ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
    MJDLIST[ifilt]   = (double*)malloc( MXEPSIM_PERFILT * MEMD);
    M5SIGLIST[ifilt] = (double*)malloc( MXEPSIM_PERFILT * MEMD);
  }

  // =======================================
  // open Dump SIMLIB to fitres-style file with 1 line per LIB

  if ( LDMP_SEQ_TEXT ) {
    sprintf(SIMLIB_DUMPFILE_SEQ, "SIMLIB_DUMP_SUMMARY_%s-%s.TEXT", 
	    GENLC.SURVEY_NAME, SIMLIB_GLOBAL_HEADER.FILTERS );

    fpdmp0 = fopen(SIMLIB_DUMPFILE_SEQ, "wt") ;    
    if ( !fpdmp0 ) {
      sprintf(c1err,"Could not open DUMP_LIBID file to write:");
      sprintf(c2err,"%s", SIMLIB_DUMPFILE_SEQ );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }


    NVAR = 10 + 4*GENLC.NFILTDEF_OBS ; 

    // xxxx    fprintf(fpdmp0,"NVAR: %d \n", NVAR );
    fprintf(fpdmp0,"VARNAMES: ROW LIBID RA DEC FIELD MWEBV GAPMAX GAPAVG "
	    "NOBS MJDMIN MJDMAX ");

    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;
      sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
      fprintf(fpdmp0,"N_%s ZPT_%s PSF_%s M5SIG_%s ", 
	      cfilt, cfilt, cfilt, cfilt );
            
    }
    fprintf(fpdmp0,"\n\n");

  } // LDMP_LIBID


  // open full table for each MJD/epoch
  if ( LDMP_ROOT ) {
    SIMLIB_DUMP_openTable(LDMP_OBS_TEXT, LDMP_ROOT);
    SIMLIB_DUMP_makeTable(LDMP_OBS_TEXT, LDMP_ROOT);
  }

  if ( LDMP_OBS_TEXT ) {
    sprintf(SIMLIB_DUMPFILE_OBS,"SIMLIB_DUMP_OBS_%s-%s.TEXT", 
	    GENLC.SURVEY_NAME, SIMLIB_GLOBAL_HEADER.FILTERS );

    fpdmp1 = fopen(SIMLIB_DUMPFILE_OBS, "wt") ; 

    NVAR = 11 ;
    // xxxx mark    fprintf(fpdmp1,"NVAR: %d \n", NVAR );
    fprintf(fpdmp1,"VARNAMES: ROW "
	    "LIBID RA DEC MJD BAND ZP_pe SKYMAG PSF M5SIG MJD_DIF\n");
    fprintf(fpdmp1,"\n");
  }



  // stdout table header

  if ( !QUIET ) {
    printf("\n  LIBID  MJD-range    NEPOCH(all,%s)    GAPMAX(frac) <GAP> \n", 
	   INPUTS.GENFILTERS );
    printf(" ------------------------------------------------------------- \n");
    fflush(stdout);
  }
    if ( INPUTS.NGENTOT_LC <= 0 ) { INPUTS.NGENTOT_LC = 9999999; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  while ( NREAD < MXREAD_SIMLIB ) {

    SIMLIB_READ_DRIVER();

    if ( SIMLIB_HEADER.NWRAP   > 0 ) { goto DONE_READ ; }
    if ( NREAD > INPUTS.NGENTOT_LC ) { goto DONE_READ ; } // May 2 2017

    ID    = GENLC.SIMLIB_ID; 
    RA    = GENLC.RA ;
    DEC   = GENLC.DEC ;
    MWEBV = GENLC.MWEBV;

    RA4  = (float)GENLC.RA ;   // need float version for pointers
    DEC4 = (float)GENLC.DEC ;

    // keep track of global min/max
    if ( RA   < GLOBAL_RANGE_RA[0]   ) { GLOBAL_RANGE_RA[0]   = RA ; }
    if ( RA   > GLOBAL_RANGE_RA[1]   ) { GLOBAL_RANGE_RA[1]   = RA ; }
    if ( DEC  < GLOBAL_RANGE_DEC[0]  ) { GLOBAL_RANGE_DEC[0]  = DEC ; }
    if ( DEC  > GLOBAL_RANGE_DEC[1]  ) { GLOBAL_RANGE_DEC[1]  = DEC ; }

    NREAD++ ;
    NREAD_SIMLIB = NREAD;  // set global 

    // store summary info for each SIMLIB entry

    zero_SIMLIB_DUMP(&SIMLIB_DUMP_AVG1);

    // for LCLIB, sum weights
    if ( INDEX_GENMODEL == MODEL_LCLIB ) {
      wgt_LCLIB  = GALrate_model(GENLC.GLON, GENLC.GLAT, &INPUTS.RATEPAR);
      wgt_LCLIB /= INPUTS.RATEPAR.RATEMAX ;
      wgtsum_LCLIB += wgt_LCLIB ;
    }

    MJDMIN4 =  99999. ;
    MJDMAX4 = -99999. ;
    sprintf(FIELDNAME, "%s", SIMLIB_OBS_GEN.FIELDNAME[1]) ;

    for ( iep = 1; iep <= GENLC.NEPOCH; iep++ ) {
      ifilt_obs = GENLC.IFILT_OBS[iep] ;

      MJD_LAST = -9.0; if(iep>1) { MJD_LAST=GENLC.MJD[iep-1]; }
      MJD = GENLC.MJD[iep];

      ZPTERR    = SIMLIB_OBS_GEN.ZPTERR[iep]; 
      if ( ZPTERR > ZPTERR_MAX ) { continue;} // exclude extreme variations

      if ( MJD < MJDMIN4 ) { MJDMIN4 = MJD ; }
      if ( MJD > MJDMAX4 ) { MJDMAX4 = MJD ; }

      SIMLIB_DUMP_AVG1.NEPFILT[0]         += 1.0 ;
      SIMLIB_DUMP_AVG1.NEPFILT[ifilt_obs] += 1.0 ;

      Nobs = (int)SIMLIB_DUMP_AVG1.NEPFILT[0] ;
      MJDLIST_ALL[Nobs] = MJD ; 
      Nobs = (int)SIMLIB_DUMP_AVG1.NEPFILT[ifilt_obs] ;
      MJDLIST[ifilt_obs][Nobs] = MJD ; 

      GAIN_SIMLIB    = SIMLIB_OBS_GEN.CCDGAIN[iep] ;
      ZPT_SIMLIB     = SIMLIB_OBS_GEN.ZPTADU[iep] ;
      SKYSIG_ADU     = SIMLIB_OBS_GEN.SKYSIG[iep] ;
      PSF_SIMLIB     = SIMLIB_OBS_GEN.PSFSIG1[iep] ;
      PIXSIZE_SIMLIB = SIMLIB_OBS_GEN.PIXSIZE[1] ;

      // convert PSF(pixels,sigma) into PSF(arcsec,FWHM)
      PSF = PSF_SIMLIB * PIXSIZE_SIMLIB * FWHM_SIGMA_RATIO ; 

      // convert ZPT(ADU) into ZPT(p.e.)
      ZPT_pe = ZPT_SIMLIB + 2.5*log10f(GAIN_SIMLIB);

      SKYSIG_pe  = SKYSIG_ADU * GAIN_SIMLIB ;

      // convert SKYSIG (ADU/pixel) into Perry mag/arcsec^2
      TMP     = pow(SKYSIG_pe,2.0) ;            // sky level in p.e., per pixel
      FSKY_pe = TMP / pow(PIXSIZE_SIMLIB,2.0) ;  // sky flux in p.e. per arcsec
      if ( FSKY_pe > 1.0E-9 ) {
	SKYMAG = ZPT_pe - 2.5*log10(FSKY_pe); // flux -> mag conversion
      }
      else { SKYMAG = 0.0; }

      
      // calculate 5 sigma limiting mag
      SNR_maglimit = 5.0 ;
      M5SIG = MAGLIMIT_calculator(ZPT_pe,PSF,SKYMAG, SNR_maglimit);

      // store M5SIG list for this entry 
      Nobs = (int)SIMLIB_DUMP_AVG1.NEPFILT[ifilt_obs] ;
      M5SIGLIST[ifilt_obs][Nobs] = M5SIG ; 

      SIMLIB_DUMP_AVG1.GAIN[ifilt_obs]       += GAIN_SIMLIB ;
      SIMLIB_DUMP_AVG1.ZPT[ifilt_obs]        += ZPT_pe ;
      SIMLIB_DUMP_AVG1.PSF[ifilt_obs]        += PSF ;
      SIMLIB_DUMP_AVG1.SKYSIG_ADU[ifilt_obs] += SKYSIG_ADU ;
      SIMLIB_DUMP_AVG1.SKYSIG_pe[ifilt_obs]  += SKYSIG_pe ;
      SIMLIB_DUMP_AVG1.SKYMAG[ifilt_obs]     += SKYMAG;
      SIMLIB_DUMP_AVG1.M5SIG[ifilt_obs]      += M5SIG ;


      if ( LDMP_OBS_TEXT ) {
	NROW_MJD++ ;
	sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
	fprintf(fpdmp1,"ROW: %3d %4d %.4f %.4f %.3f  "
		"%s  %.3f  %.3f  %.3f %.3f %.3f \n",
		NROW_MJD, ID, RA, DEC, MJD, 
		cfilt, ZPT_pe, SKYMAG, PSF,  M5SIG, MJD-MJD_LAST );
      }



    } // end of 'ep' epoch loop for this simlib entry


    // get max temporal gap for unsorted list of MJDs
    Nobs   = (int)SIMLIB_DUMP_AVG1.NEPFILT[0] ;
    ptrmjd = &MJDLIST_ALL[1] ;
    MJDGAP(Nobs, ptrmjd, MJDGAP_IGNORE, &GAPMAX, &GAPAVG );
    SIMLIB_DUMP_AVG1.GAPMAX[0]  = GAPMAX ;
    SIMLIB_DUMP_AVG1.GAPAVG[0]  = GAPAVG ;
    
    if ( LDMP_SEQ_TEXT ) {
      NROW++ ;
      fprintf(fpdmp0,"ROW: %4d %4d %7.3f %7.3f %s %6.3f %5.0f %3.1f "
	      "%3d %.2f %.2f ", 
	      NROW, ID, RA, DEC, FIELDNAME, MWEBV, GAPMAX, GAPAVG, 
	      Nobs, MJDMIN4, MJDMAX4 );
      fflush(fpdmp0);
    }

    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;
      sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );

      Nobs   = (int)SIMLIB_DUMP_AVG1.NEPFILT[ifilt_obs] ;
      ptrmjd = &MJDLIST[ifilt_obs][1] ;
      MJDGAP(Nobs, ptrmjd, MJDGAP_IGNORE, &GAPMAX, &GAPAVG );
      SIMLIB_DUMP_AVG1.GAPMAX[ifilt_obs] = GAPMAX ;
      SIMLIB_DUMP_AVG1.GAPAVG[ifilt_obs] = GAPAVG ;

      XNobs = (double)Nobs;
      if ( XNobs == 0.0 ) XNobs = 1.0E12;
      SIMLIB_DUMP_AVG1.GAIN[ifilt_obs]       /= XNobs ;
      SIMLIB_DUMP_AVG1.ZPT[ifilt_obs]        /= XNobs ;
      SIMLIB_DUMP_AVG1.PSF[ifilt_obs]        /= XNobs ;
      SIMLIB_DUMP_AVG1.SKYSIG_ADU[ifilt_obs] /= XNobs ;
      SIMLIB_DUMP_AVG1.SKYSIG_pe[ifilt_obs]  /= XNobs ;
      SIMLIB_DUMP_AVG1.SKYMAG[ifilt_obs]     /= XNobs ;
      SIMLIB_DUMP_AVG1.M5SIG[ifilt_obs]      /= XNobs ;


      if ( LDMP_SEQ_TEXT ) {
	fprintf(fpdmp0,"%2d %6.2f %5.2f %5.2f ", Nobs
		, SIMLIB_DUMP_AVG1.ZPT[ifilt_obs]
		, SIMLIB_DUMP_AVG1.PSF[ifilt_obs]
		, SIMLIB_DUMP_AVG1.M5SIG[ifilt_obs]  );
	fflush(fpdmp0);
      }

      // load global array
      SIMLIB_DUMP_AVG1.FOM[ifilt_obs] = -9.0 ;

    } // end of ifilt loop over filters

    SIMLIB_DUMP_AVG1.MJDMIN = MJDMIN4 ;
    SIMLIB_DUMP_AVG1.MJDMAX = MJDMAX4 ;
    
    if ( MJDMIN4 < SIMLIB_DUMP_AVGALL.MJDMIN  ) 
      { SIMLIB_DUMP_AVGALL.MJDMIN = MJDMIN4; }
    if ( MJDMAX4 > SIMLIB_DUMP_AVGALL.MJDMAX  ) 
      { SIMLIB_DUMP_AVGALL.MJDMAX = MJDMAX4; }

    // check user-cuts that would have failed
    for ( icut=1; icut <= NCUT_SIMLIB; icut++ ) {
      TMP0 = *PTR_SIMLIB_CUTVAR[icut][0];
      TMP1 = *PTR_SIMLIB_CUTVAR[icut][1];
      if ( TMP0 < SIMLIB_GENRANGE[icut][0] ) { NSIMLIB_CUTFAIL[icut]++ ; }
      if ( TMP1 > SIMLIB_GENRANGE[icut][1] ) { NSIMLIB_CUTFAIL[icut]++ ; }
    }

    if ( LDMP_SEQ_TEXT ) { fprintf(fpdmp0,"\n" ); }
    
    // screen dump
    if ( LDMP_LOCAL  && NREAD <= MXSIMLIB_DUMP_STDOUT ) {

      GAPMAX = SIMLIB_DUMP_AVG1.GAPMAX[0];
      GAPAVG = SIMLIB_DUMP_AVG1.GAPAVG[0];
      MJDWIN = SIMLIB_DUMP_AVG1.MJDMAX - SIMLIB_DUMP_AVG1.MJDMIN;
      FRAC   = GAPMAX/MJDWIN;

      
      printf("  %4.4d   %5.0f-%5.0f  %3d,"
	     ,ID, MJDMIN4, MJDMAX4, (int)SIMLIB_DUMP_AVG1.NEPFILT[0] );

      for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
	ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;
	printf("%2d ", (int)SIMLIB_DUMP_AVG1.NEPFILT[ifilt_obs] );
      }

      printf("  %5.1f(%3.2f) %4.1f \n", GAPMAX, FRAC, GAPAVG );

      fflush(stdout);

    } // end of LDMP_LOCAL - screen dump

    update_SIMLIB_DUMP_AVGALL(1);
    if ( LDMP_ROOT) {
      SNTABLE_FILL(TABLEID_SIMLIB_DUMP); 
    }

  } // end of READ loop



  // ======================================
  // summmarize average quantities for each filter

 DONE_READ:

  if ( QUIET ) { return; }

  if ( LDMP_SEQ_TEXT ) { 
    fclose(fpdmp0); 
    printf("\n One-line dump per LIBID-sequence written to \n\t '%s' \n", 
	   SIMLIB_DUMPFILE_SEQ );
  }

  if ( LDMP_OBS_TEXT   ) { 
    fclose(fpdmp1); 
    printf("\n One-line dump per OBS written to \n\t '%s' \n", 
	   SIMLIB_DUMPFILE_OBS );
  }

  if ( LDMP_ROOT ) 
    {  TABLEFILE_CLOSE(SIMLIB_DUMPFILE_ROOT); }

  // compute grand-average quantities over all LIBIDs (for each filter)
  update_SIMLIB_DUMP_AVGALL(2);


  printf("\n Done reading %d SIMLIB entries. \n", NREAD );
  
  printf("\n LIBRARY AVERAGES PER FILTER:  \n");
  printf("                   <PSF>  \n" );
  printf("          <ZPT-pe> FWHM  <SKYSIG>  <SKYMAG>                    \n");
  printf("      FLT  (mag)  (asec) (pe/pix)  (asec^-2) <m5sig> <Nep> <GAP>\n");
  printf(" --------------------------------------------------------------- \n");

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {

    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;

    TMP = SIMLIB_DUMP_AVGALL.SKYSIG_pe[ifilt_obs] ; 
    if ( TMP < 10.0 )
      { sprintf(ctmp,"%8.3f", TMP);    }
    else 
      { sprintf(ctmp,"%8.1f", TMP); }

    //      FLT   ZPT     PSF    SKYSIG  SKYMAG  M5SIG  NEP   <GAP> 
    printf("BAND:  "
	   "%c  %6.2f   %5.3f  %8.8s   %5.2f   %6.2f  %5.2f  %4.1f\n",
	   FILTERSTRING[ifilt_obs]
	   ,SIMLIB_DUMP_AVGALL.ZPT[ifilt_obs]
	   ,SIMLIB_DUMP_AVGALL.PSF[ifilt_obs]
	   ,ctmp  // AVG_SKYSIG_pe[ifilt_obs]
	   ,SIMLIB_DUMP_AVGALL.SKYMAG[ifilt_obs]
	   ,SIMLIB_DUMP_AVGALL.M5SIG[ifilt_obs] 
	   ,SIMLIB_DUMP_AVGALL.NEPFILT[ifilt_obs] 
	   ,SIMLIB_DUMP_AVGALL.GAPAVG[ifilt_obs]
	   // ccc	   ,SIMLIB_DUMP_AVGALL.FOM[ifilt_obs]
	   );
    fflush(stdout);    
  }


  // ==========================
  // library ranges (min/max)
  printf("\n LIBRARY MIN-MAX RANGES:  \n");

  printf("\t RA:   %8.3f  to  %8.3f deg \n", 
	 GLOBAL_RANGE_RA[0], GLOBAL_RANGE_RA[1] );
  printf("\t DEC: %8.3f  to  %8.3f deg \n", 
	 GLOBAL_RANGE_DEC[0], GLOBAL_RANGE_DEC[1] );
  printf("\t MJD:  %8.1f  to  %8.1f  (GENRANGE_PEAKMJD = %5.0f-%5.0f)\n", 
	 SIMLIB_DUMP_AVGALL.MJDMIN, SIMLIB_DUMP_AVGALL.MJDMAX ,
	 INPUTS.GENRANGE_PEAKMJD[0], INPUTS.GENRANGE_PEAKMJD[1] );

  
  fflush(stdout);

  // ==============
  // check for user-cuts that fail some of the simlibs

  printf("\n");
  for ( icut=1; icut <= NCUT_SIMLIB; icut++ ) {
    Ncut = NSIMLIB_CUTFAIL[icut] ;
    if ( Ncut > 0 ) {
      printf("  CUT-WARNING: %4d SIMLIBS will fail user-cut on '%s' \n",
	     Ncut, SIMLIB_CUTVARNAME[icut] );
    }
  }

  // print average wgt for LCLIB (Jun 2018)
  if ( INDEX_GENMODEL == MODEL_LCLIB ) {
    double avgwgt = wgtsum_LCLIB / (double)NREAD;
    printf("\t Average LCLIB-wgt = %f/%d = %f \n",
	   wgtsum_LCLIB, NREAD, avgwgt);
  }

  // ============
  happyend();

  return ;

} // end of SIMLIB_DUMP_DRIVER


// =========================================================
void update_SIMLIB_DUMP_AVGALL(int OPT) {

  // Created Aug 2017
  // Update SIMLIB_DUMP_AVGALL to average over entire library.
  // OPT=1 --> increment sum over LIBID
  // OPT=2 --> take avg

  int ifilt_obs, ifilt ;
  double VAL, XN  ;
  char  cfilt[2];
  //  char fnam[] = "update_SIMLIB_DUMP_AVGALL" ;

  // ------------ BEGIN -----------

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {

    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt] ;
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );

    if ( OPT == 1 ) {  
      // increment sums for AVGALL

      VAL  = SIMLIB_DUMP_AVG1.ZPT[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.ZPT[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.ZPT[ifilt_obs]  += VAL; 
      }
      
      VAL  = SIMLIB_DUMP_AVG1.PSF[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.PSF[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.PSF[ifilt_obs]  += VAL; 
      }
      
      VAL  = SIMLIB_DUMP_AVG1.SKYSIG_ADU[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.SKYSIG_ADU[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.SKYSIG_ADU[ifilt_obs]  += VAL; 
      }
      
      VAL  = SIMLIB_DUMP_AVG1.SKYSIG_pe[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.SKYSIG_pe[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.SKYSIG_pe[ifilt_obs]  += VAL; 
      }
      
      VAL  = SIMLIB_DUMP_AVG1.SKYMAG[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.SKYMAG[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.SKYMAG[ifilt_obs]  += VAL; 
      }
      
      VAL  = SIMLIB_DUMP_AVG1.M5SIG[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.M5SIG[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.M5SIG[ifilt_obs]  += VAL; 
      }
      
      
      VAL  = SIMLIB_DUMP_AVG1.GAPAVG[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.GAPAVG[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.GAPAVG[ifilt_obs]  += VAL; 
      }
      VAL  = SIMLIB_DUMP_AVG1.GAPMAX[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.GAPMAX[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.GAPMAX[ifilt_obs]  += VAL; 
      }
      
      VAL  = SIMLIB_DUMP_AVG1.FOM[ifilt_obs] ; 
      if ( VAL > .001 )  { 
	SIMLIB_DUMP_NAVGALL.FOM[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.FOM[ifilt_obs]  += VAL; 
      }
      
      VAL  = SIMLIB_DUMP_AVG1.NEPFILT[ifilt_obs] ; 
      if ( VAL > -8 )  { 
	SIMLIB_DUMP_NAVGALL.NEPFILT[ifilt_obs] += 1.0 ;
	SIMLIB_DUMP_AVGALL.NEPFILT[ifilt_obs]  += VAL; 
      }
    }
    else if ( OPT == 2 ) {
      // AVG = SUM/N

      XN  = SIMLIB_DUMP_NAVGALL.ZPT[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.ZPT[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.ZPT[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.PSF[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.PSF[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.PSF[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.SKYSIG_ADU[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.SKYSIG_ADU[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.SKYSIG_ADU[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.SKYSIG_pe[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.SKYSIG_pe[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.SKYSIG_pe[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.SKYMAG[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.SKYMAG[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.SKYMAG[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.M5SIG[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.M5SIG[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.M5SIG[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.GAPMAX[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.GAPMAX[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.GAPMAX[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.GAPAVG[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.GAPAVG[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.GAPAVG[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.FOM[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.FOM[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.FOM[ifilt_obs] = VAL/XN ;

      XN  = SIMLIB_DUMP_NAVGALL.NEPFILT[ifilt_obs] + 1.0E-9 ;
      VAL = SIMLIB_DUMP_AVGALL.NEPFILT[ifilt_obs]; 
      SIMLIB_DUMP_AVGALL.NEPFILT[ifilt_obs] = VAL/XN ;

    }


  } // end ifilt loop

  return ;

} // end update_SIMLIB_DUMP_AVGALL

// =========================================================
void zero_SIMLIB_DUMP(SIMLIB_DUMP_DEF *SIMLIB_DUMP) {

  int ifilt ; 

  SIMLIB_DUMP->MJDMIN =  9999999. ;
  SIMLIB_DUMP->MJDMAX = -9999999. ;

  for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
    SIMLIB_DUMP->NEPFILT[ifilt]      =  0   ;
    SIMLIB_DUMP->NEPFILT[ifilt]      =  0   ;
    SIMLIB_DUMP->ZPT[ifilt]          =  0.0 ;
    SIMLIB_DUMP->PSF[ifilt]          =  0.0 ;
    SIMLIB_DUMP->SKYSIG_ADU[ifilt]   =  0.0 ;
    SIMLIB_DUMP->SKYSIG_pe[ifilt]    =  0.0 ;
    SIMLIB_DUMP->SKYMAG[ifilt]       =  0.0 ;
    SIMLIB_DUMP->M5SIG[ifilt]        =  0.0 ;
    SIMLIB_DUMP->GAPMAX[ifilt]       =  0.0 ;
    SIMLIB_DUMP->GAPAVG[ifilt]       =  0.0 ;
  }

} // end zero_SIMLIB_DUMP

// =========================================================
void SIMLIB_DUMP_openTable(int LDMP_MJD_TEXT,int LDMP_ROOT) {

  // Created Aug 4 2017 
  // Check option to open table(s):
  // LDMP_MJD_TEXT --> open TEXT-format table
  // LDMP_ROOT     --> open ROOT table for entire SIMLIB
  //

  int NOUT=0, IFILETYPE; 
  char PREFIX[MXPATHLEN], openOpt[40] ;
  //  char fnam[] = "SIMLIB_DUMP_openTable" ;

  // -------------- BEGIN -------------

  if ( LDMP_MJD_TEXT==0 && LDMP_ROOT==0 ) { return ; }

  TABLEFILE_INIT(); 

  sprintf(PREFIX,"SIMLIB_DUMP_%s-%s", 
	  GENLC.SURVEY_NAME, SIMLIB_GLOBAL_HEADER.FILTERS );

#ifdef USE_ROOT
  if ( LDMP_ROOT )  {
    sprintf(SIMLIB_DUMPFILE_ROOT,"%s.ROOT", PREFIX);
    sprintf(openOpt,"root new");
    IFILETYPE = TABLEFILE_OPEN(SIMLIB_DUMPFILE_ROOT,openOpt);
    NOUT++ ;
  }
#endif

  return ;

} // end SIMLIB_DUMP_openTable

// =========================================================
void SIMLIB_DUMP_makeTable(int LDMP_MJD_TEXT,int LDMP_ROOT) {

  // Created Aug 4 2017
  // Create table (text and/or root) for SIMLIB DUMP
  //
  // !!! NOT READY !!!
  // 
  char varName[40];
  char BLOCKVAR[] = "VAR";
  //  char fnam[] = "SIMLIB_DUMP_makeTable" ;

  // --------------- BEGIN --------------

  if ( LDMP_MJD_TEXT==0 && LDMP_ROOT==0 ) { return ; }

  SNTABLE_CREATE(TABLEID_SIMLIB_DUMP, TABLENAME_SIMLIB_DUMP, "KEY" ) ;


  // define columns; first identifier must be a string
  sprintf(varName, "ROW:C*12"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 SIMLIB_HEADER.LIBNAME, varName, 1 ) ;

  sprintf(varName, "LIBID:I"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.LIBID, varName, 1 ) ;

  sprintf(varName, "CCDNUM:I"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.CCDNUM, varName, 1 ) ;

  sprintf(varName, "FIELD:C*20"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 SIMLIB_HEADER.FIELD, varName, 1 ) ;

  sprintf(varName, "RA:D"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.RA, varName, 1 ) ;

  sprintf(varName, "DEC:D"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.DEC, varName, 1 ) ;
 
  sprintf(varName, "MWEBV:D"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.MWEBV, varName, 1 ) ;

  sprintf(varName, "PIXSIZE:D"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.PIXSIZE, varName, 1 ) ;

  sprintf(varName, "REDSHIFT_MIN:D"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.GENRANGE_REDSHIFT[0], varName, 1 ) ;
  sprintf(varName, "REDSHIFT_MAX:D"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.GENRANGE_REDSHIFT[1], varName, 1 ) ;

  sprintf(varName, "PEAKMJD_MIN:D"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[0], varName, 1 );
  sprintf(varName, "PEAKMJD_MAX:D"); 
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[1], varName, 1 );

  // ---------------------------------
  // --- start epoch-vectors ---------
  // ---------------------------------


  sprintf(varName,"NOBS[%d]:I", MXOBS_SIMLIB );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_HEADER.NOBS, varName, 1 );

  sprintf(varName,"MJD(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.MJD[1], varName, 1 );

  sprintf(varName,"IDEXPT(NOBS):I" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.IDEXPT[1], varName, 1 );

  sprintf(varName,"FIELD(NOBS):C*12" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 SIMLIB_OBS_RAW.PTR_FIELDNAME[1], varName, 1 );

  sprintf(varName,"BAND(NOBS):C*4" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 SIMLIB_OBS_RAW.PTR_BAND[1], varName, 1 );

  sprintf(varName,"CCDGAIN(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.CCDGAIN[1], varName, 1 );

  sprintf(varName,"READNOISE(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.READNOISE[1], varName, 1 );

  sprintf(varName,"SKYSIG(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.SKYSIG[1], varName, 1 );

  sprintf(varName,"PSFSIG1(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.PSFSIG1[1], varName, 1 );
  sprintf(varName,"PSFSIG2(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.PSFSIG2[1], varName, 1 );
  sprintf(varName,"PSFRATIO(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.PSFRATIO[1], varName, 1 );

  sprintf(varName,"ZPTADU(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.ZPTADU[1], varName, 1 );
  sprintf(varName,"ZPTERR(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.ZPTERR[1], varName, 1 );

  sprintf(varName,"MAG(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.MAG[1], varName, 1 );

  // - - -  template noise - - - - 
  sprintf(varName,"TEMPLATE_SKYSIG(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.TEMPLATE_SKYSIG[1], varName, 1 );

  sprintf(varName,"TEMPLATE_ZPT(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.TEMPLATE_ZPT[1], varName, 1 );

  return;

} // end SIMLIB_DUMP_makeTable

// ===================================
void MJDGAP(int NLIST, double *MJDLIST, double MJDGAP_IGNORE,
	    double *GAPMAX, double *GAPAVG ) {

  // Created May 14, 2009 by R.Kessler
  // Return GAPMAX and GAPAVG for input *MJDLIST
  // and NLIST is the size of *MJDLIST
  // Note that MJDLIST can be unsorted.
  // GAPS > MJDGAP_IGNORE (days) are ignored.
  //
  // Jun 10, 2009: set *GAPAVG=0 when NGAP=0 (fixes nan bug)
  //
  // Oct 14, 2010: fix aweful bug by passing float MJDLIST4 to 
  //               SORTZV instead of passing double MJDLIST
  //
  // Jan 21, 2011: count GAPSUM (for avg) only if GAP < 100 days;
  //               Assume that >100 day gaps are seasonal.
  //
  // Jan 31, 2011: init *GAPMAX = *GAPAVG = 0 in case NLIST=1
  //
  // Jun 07, 2011: pass MJDGAP_IGNORE to replace hard-wired '100'
  //
  // Jan 12, 2013: replace CERNLIB sortzv with floatDouble wrapper.
  //               No more MJDLIST4 since we can sort in double now.
  //
  // May 20 2017: INDEX_SORT[1000] -> INDEX_SORT[MXREAD_SIMLIB],
  //
  // Jan 25 2018: INDEX_SORT[MXREAD_SIMLIB] -> INDEX_SORT[MXOBS_SIMLIB]

  double  MJD, MJDLAST, GAP, GAPSUM ;
  int     INDEX_SORT[MXOBS_SIMLIB], ORDER_SORT, i, isort, NGAP ;

  // --------------- BEGIN ---------

  *GAPAVG = 0.0;
  *GAPMAX = 0.0;
  if ( NLIST <= 1 ) return ;

  // first sort the MJDs to be chronological

  ORDER_SORT = +1 ; // increasing
  sortDouble(NLIST, MJDLIST, ORDER_SORT, INDEX_SORT );


  // note that we start at 2nd MJD on list

   GAPSUM = 0.0  ;
  *GAPMAX = -9.0 ;
  *GAPAVG = -9.0 ;
  NGAP = 0;
 
  for ( i=1; i < NLIST; i++ ) {

    isort    = INDEX_SORT[i-1]  ;
    MJDLAST  = MJDLIST[isort];

    isort = INDEX_SORT[i] ;
    MJD   = MJDLIST[isort];

    GAP = MJD - MJDLAST ;

    // average gap only if MJD changes and if gap is < MJDGAP_IGNORE
    if ( GAP > 0.001 && GAP < MJDGAP_IGNORE ) {
      GAPSUM += GAP;
      NGAP++;
    }

    if ( GAP > *GAPMAX ) *GAPMAX = GAP ;
  }

  if ( NGAP > 0 ) 
    { *GAPAVG = GAPSUM/(double)(NGAP); }
  else
    {  *GAPAVG = 0.0 ; }

  //  printf(" MAXGAP = %f for NLIST=%d MJDs\n", *GAPMAX, NLIST );


} // end of MJDGAP



// ************************************
double SIMLIB_angsep_min(int NSTORE, double RA, double DEC, 
			  double *RA_STORE, double *DEC_STORE) {

  // OBSOLETE (Aug 2017) ??
  // return min angular separation (degrees) between RA,DEC and the
  // passed array of NSTORE RA,DEC values.

  int i;

  double     
     RA_TMP, DEC_TMP
    ,RA_RAD, DEC_RAD
    ,ANGSEP, ANGSEP_MIN 
    ,DOTPROD
    ,X, X_TMP
    ,Y, Y_TMP
    ,Z, Z_TMP
    ;

  char fnam[] = "SIMLIB_angsep_min" ;

  // ------------- BEGIN -----------

  if ( NSTORE <= 0 ) { return  (TWOPI / RADIAN) ;  }

  // get polar coord of ref
  X = cos(RA*RADIAN) * cos(DEC*RADIAN);
  Y = sin(RA*RADIAN) * cos(DEC*RADIAN);
  Z = sin(DEC*RADIAN);


  ANGSEP_MIN = 99999. ;

  for ( i = 0; i < NSTORE; i++ ) {
    RA_TMP   = *(RA_STORE+i);     // phi
    DEC_TMP  = *(DEC_STORE+i);   // theta

    RA_RAD   = RADIAN * RA_TMP ;
    DEC_RAD  = RADIAN * DEC_TMP ; 
    X_TMP = cos(RA_RAD) * cos(DEC_RAD);
    Y_TMP = sin(RA_RAD) * cos(DEC_RAD);
    Z_TMP = sin(DEC_RAD);

    DOTPROD = .99999999*(X*X_TMP + Y*Y_TMP + Z*Z_TMP) ;
    
    if ( fabs(DOTPROD) > 1.000000 ) {
      sprintf(c1err,"DOTPROD = %f", DOTPROD);
      sprintf(c2err,"RA,DEC=%f,%f  STORE[RA,DEC]=%f,%f", 
	      RA, DEC, RA_TMP, DEC_TMP );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    ANGSEP  = acos(DOTPROD) ;
    if ( ANGSEP < ANGSEP_MIN )
      {  ANGSEP_MIN = ANGSEP ; }
  }

  // return min angular sep in degrees.

  return(ANGSEP_MIN / RADIAN) ;


} // end of SIMLIB_ANGSEP_MIN

// ***********************************
void DUMP_GENMAG_DRIVER(void) {

  // compute mag in 1-day bins and write to text file.
  //
  // Jan 31 2017: major overhaul and refactor.
  //   + write 1 TEXT tableFile instead of PAW-readable file per band.
  //   + allow delta-function in Trest and/or shapePar
  //
  //

#define MXSHAPEPAR 8

  float TEST_DELTA[MXSHAPEPAR] =  // mlcs2k2
    { -0.4, -0.2, 0.0, +0.2, +0.4, +0.6,  +1.0, +1.5 } ;
  float TEST_SALT2x1[MXSHAPEPAR]   =  // SALt2
    { -5.0, -3.0, -1.0, 0.0, +1.0, +2.0, +3.0, +5.0 } ;
  float TEST_DM15[MXSHAPEPAR]   =  // SNOOPY
    { +1.9, +1.7, +1.5, 1.3, 1.2, 1.1, 1.0,   0.85 } ;
  float TEST_STRETCH[MXSHAPEPAR] =  // stretch
    { 0.5, 0.7, 0.8, 0.9, 1.0,  1.1, 1.2,  1.3 } ;
  float TEST_COSANGLE[MXSHAPEPAR] =  // SIMSED with COSANGLE param
    { -0.9669, -0.70, -0.433, -0.233, 0.233,  0.433, 0.70,  0.9669 } ;
  float TEST_NON1ASED = 1.0 ;

  float 
    *ptr_shapepar
    ,Tmin, Tmax, lamdif, magtmp, magerrtmp, z4, shapepar
    ,GENMAG[MXEPSIM][MXSHAPEPAR]
    ,GENMAGERR[MXEPSIM][MXSHAPEPAR]
    ,DM15_CALC[MXFILTINDX][MXSHAPEPAR]
    ;

  int 
    epoch, ifilt_rest, ifilt_obs, ifilt, NROW
    ,N, NSHAPEPAR, ishape, irank, colopt
    ;

  double 
    mag8, *ptr_genmag8, MU8, ARG8, z8
    ,magerr8, *ptr_generr8, Trest8
      ;

  FILE *FPDMP ;
  char dmpFile[200], cfilt[2] ;
  char fnam[]  = "DUMP_GENMAG_DRIVER"    ;

  // --------- BEGIN ---------

  Tmin = INPUTS.GENRANGE_DMPTREST[0];
  Tmax = INPUTS.GENRANGE_DMPTREST[1];
  if ( Tmin > Tmax ) { 
    sprintf(c1err,"Invalid Tmin=%.3f exceeds Tmax=%.3f", Tmin, Tmax);
    sprintf(c2err,"Check GENRANGE_DMPTREST key in sim-input file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  NSHAPEPAR = MXSHAPEPAR;

  z8 = INPUTS.GENRANGE_REDSHIFT[0] ;

  if ( z8 > 1.0E-8 ) 
    { GENLC.REDSHIFT_CMB       = z8; }
  else
    { GENLC.REDSHIFT_CMB       = ZAT10PC ; }  // 10 pc

  GENLC.REDSHIFT_HELIO = GENLC.REDSHIFT_CMB ; // ??

  gen_distanceMag(GENLC.REDSHIFT_CMB, GENLC.REDSHIFT_HELIO,
		  &GENLC.DLMU, &GENLC.LENSDMU);

  GENLC.AV             =  0.0 ;
  GENLC.SALT2c         =  0.0 ;
  GENLC.NEPOCH         =  0   ;
  MU8 = GENLC.DLMU ;

  //  init_RANLIST();

  // --------------------------
  // set epoch-array to 1 day bins, unless delta-function is requested
  if ( Tmin == Tmax ) {
    GENLC.NEPOCH = 1 ;
    GENLC.epoch_rest[1] = (double)INPUTS.GENRANGE_DMPTREST[0] ;
    GENLC.epoch_obs[1]  = (double)INPUTS.GENRANGE_DMPTREST[0] ;
  }
  else {
    for ( epoch = (int)Tmin; epoch <= (int)Tmax; epoch++ ) {
      GENLC.NEPOCH++ ;   N = GENLC.NEPOCH;      
      if ( N >= MXEPSIM ) {
	sprintf(c1err,"Nepoch =%d exceeds array bound ", N);
	sprintf(c2err,"(MXEPSIM = %d)", MXEPSIM);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }      
      GENLC.epoch_rest[N] = (double)epoch;
      GENLC.epoch_obs[N]  = (double)epoch;
    }
  }


  // --------------------------
  // set shapepar array list; either hard-wired list above
  // or user-defined delta-function
  float Smin = INPUTS.GENGAUSS_SHAPEPAR.RANGE[0] ;
  float Smax = INPUTS.GENGAUSS_SHAPEPAR.RANGE[1] ;
  if ( Smin == Smax ) {
    NSHAPEPAR = 1 ;
    TEST_DELTA[0] = TEST_DM15[0] = TEST_STRETCH[0] = TEST_SALT2x1[0] = Smin ;
  }

  // --------------------------------------------
  // open one output dump file (FITRES format)
  sprintf(dmpFile, "DUMP_GENMAG_%s.TEXT", INPUTS.MODELNAME );
  if ( (FPDMP = fopen(dmpFile, "wt"))==NULL ) {       
    sprintf ( c1err, "Cannot open output file :" );
    sprintf ( c2err," '%s' ", dmpFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  fprintf(FPDMP,"# GENMODEL:  %s \n", INPUTS.MODELNAME );
  fprintf(FPDMP,"# REDSHIFT_HELIO: %.4f \n", GENLC.REDSHIFT_HELIO );
  fprintf(FPDMP,"#\n");
  fprintf(FPDMP,"NVAR:  6\n");
  fprintf(FPDMP,"VARNAMES: ROW BAND TREST SHAPEPAR GENMAG GENMAGERR\n");

  printf("   Open DUMP-file: %s \n", dmpFile );   
  NROW=0;


  if ( INDEX_GENMODEL == MODEL_STRETCH ) 
    {  ptr_shapepar = TEST_STRETCH; }
  else if ( INDEX_GENMODEL == MODEL_SALT2 ) 
    {  ptr_shapepar = TEST_SALT2x1 ; }
  else if ( INDEX_GENMODEL == MODEL_SIMSED ) 
    {  ptr_shapepar = TEST_COSANGLE ; }
  else if ( INDEX_GENMODEL == MODEL_MLCS2k2 ) 
    {  ptr_shapepar = TEST_DELTA ; }
  else if ( INDEX_GENMODEL == MODEL_SNOOPY ) 
    {  ptr_shapepar = TEST_DM15 ; }
  else if ( INDEX_GENMODEL == MODEL_NON1ASED ) 
    {  ptr_shapepar = &TEST_NON1ASED ;  NSHAPEPAR=1 ; }
  else {
    ptr_shapepar = TEST_STRETCH ; // avoid compile warning
    sprintf(c1err,"Do not recognize model = '%s'", INPUTS.MODELNAME );
    errmsg(SEV_FATAL, 0, fnam, c1err , "" ); 
  }

  // ------------------


  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {

    ifilt_obs  = GENLC.IFILTMAP_OBS[ifilt] ;
    GENLC.DOFILT[ifilt_obs] = 1; // Jan 31 2017; needed for restFrame model

    for ( epoch=0; epoch <= GENLC.NEPOCH; epoch++ )
      { GENLC.IFILT_OBS[epoch] = ifilt_obs; }

    ptr_genmag8   = &GENFILT.genmag_obs[ifilt_obs][0] ;
    ptr_generr8   = &GENFILT.generr_obs[ifilt_obs][0] ;
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
   
    if ( GENFRAME_OPT == GENFRAME_REST ) {
      colopt = INPUTS.KCORFLAG_COLOR ; 
      z4 = (float)GENLC.REDSHIFT_HELIO ;

      irank = 1;  
      GENLC.IFILTMAP_REST1[ifilt_obs] = 
	nearest_ifilt_rest__(&colopt, &ifilt_obs, &irank, &z4, &lamdif );
      ifilt_rest = GENLC.IFILTMAP_REST1[ifilt_obs] ;

      irank = 2;  
      GENLC.IFILTMAP_REST2[ifilt_obs] = 
	nearest_ifilt_rest__(&colopt, &ifilt_obs, &irank, &z4, &lamdif );

      irank = 3;  
      GENLC.IFILTMAP_REST3[ifilt_obs] = 
	nearest_ifilt_rest__(&colopt, &ifilt_obs, &irank, &z4, &lamdif );
      
    }

    // try NSHAPEPAR shape/luminosity variations

    for ( ishape=0; ishape < NSHAPEPAR; ishape++ ) {

      GENLC.DELTA          = TEST_DELTA[ishape] ;
      GENLC.DM15           = TEST_DM15[ishape] ;
      GENLC.STRETCH        = TEST_STRETCH[ishape] ;
      GENLC.SALT2x1        = TEST_SALT2x1[ishape] ;

      if ( INDEX_GENMODEL == MODEL_SALT2 ) {
	GENLC.SALT2alpha = (double)INPUTS.GENGAUSS_SALT2ALPHA.PEAK ;
	GENLC.SALT2beta  = (double)INPUTS.GENGAUSS_SALT2BETA.PEAK ;
	GENLC.SALT2x0 = SALT2x0calc(GENLC.SALT2alpha, GENLC.SALT2beta, 
				    GENLC.SALT2x1, GENLC.SALT2c, MU8 );
	/*
	printf(" xxxx x0=%f  x1=%f  a=%f  b=%f \n", 
	GENLC.SALT2x0, GENLC.SALT2x1, GENLC.SALT2alpha, GENLC.SALT2beta ); */
      }

      else  if ( INDEX_GENMODEL == MODEL_SIMSED ) {
	GENLC.SIMSED_PARVAL[1]  = TEST_COSANGLE[ishape] ;
	ARG8 = -0.4 * MU8 ;
	GENLC.SALT2x0 = pow( TEN , ARG8 );
      }
      else  if ( INDEX_GENMODEL == MODEL_SNOOPY ) {
	// Jan 31 2017
      }
      else  if ( INDEX_GENMODEL == MODEL_NON1ASED ) {
	pick_NON1ASED(1, &INPUTS.NON1ASED, &GENLC.NON1ASED); 
      }
      else {
	// unknown
	GENLC.SALT2x0 = -9.0;
      }

      genmodel(ifilt_obs, 1);
      if ( GENFRAME_OPT == GENFRAME_REST ) {
	genmodel(ifilt_obs,2);      // 2nd nearest filter
	genmodel(ifilt_obs,3);      // 3rd nearest filter
	genmag_boost();
      } 

      for ( epoch=1; epoch <= GENLC.NEPOCH; epoch++ ) {
	mag8    = ptr_genmag8[epoch] ;
	mag8   += GENLC.GENMAG_OFF_GLOBAL ;
	magerr8 = ptr_generr8[epoch] ;
	GENMAG[epoch][ishape]    = (float)mag8 ;
	GENMAGERR[epoch][ishape] = (float)magerr8 ;
      }

    } // end of ishape loop


    for ( epoch=1; epoch <= GENLC.NEPOCH; epoch++ ) {
      Trest8 = GENLC.epoch_rest[epoch] ;

      for ( ishape=0; ishape<NSHAPEPAR; ishape++ ) {
	magtmp    =  GENMAG[epoch][ishape] ;
	magerrtmp =  GENMAGERR[epoch][ishape] ;
	shapepar  =  ptr_shapepar[ishape];

	NROW++ ; 
	fprintf(FPDMP,"ROW: %8d  %s  %6.2f  %6.2f  %.3f  %.3f \n",
		NROW, cfilt, Trest8, shapepar,  magtmp, magerrtmp );

	if ( fabsf(Trest8) < .001 ) {
	  DM15_CALC[ifilt_obs][ishape]  = -magtmp ;
	}
	if ( fabsf(Trest8-15.0) < .001 ) {
	  DM15_CALC[ifilt_obs][ishape] +=  magtmp ;
	  fprintf(FPDMP,"DM15_CALC[%s] = %.3f\n",
		  cfilt, DM15_CALC[ifilt_obs][ishape] );
	}

      } // end ishape
    }  // end epoch 
  }  // end of ifilt loop


  // -------------------------------------------------

  fclose(FPDMP);
  printf("   Finished dump --> exit program. \n\n");
  exit(1);

  return ;

}  // end of DUMP_GENMAG_DRIVER


// ************************
void test_fortran(void) {

  char kcorFile[80] = "  " ;
  int IERR, ifilt_rest, ifilt_obs, iz;
  double t8, z8, av8, kcor8 ;

  // ----------- BEGIN -----------

  print_banner ( " TEST FORTRAN INTERFACE " );

  init_snvar__(&IERR);
  sprintf(kcorFile, "%s", INPUTS.KCOR_FILE );
  rdkcor_(kcorFile, &IERR, strlen(kcorFile) );

  av8 = 0.0;
  z8  = 0.20;
  t8  = 0.0;

  ifilt_rest = 7 ;
  ifilt_obs  = 3 ;


  for ( iz=0; iz<=10; iz++ ) {
    z8 = 0.04 * (double)iz ;
    printf(" --------------------------------------- \n");
    kcor8  = get_kcor8__( &ifilt_rest, &ifilt_obs, &t8, &z8, &av8 );
    printf(" xxx z=%4.2f  Tr=%5.2f AV=%4.2f   K_(%d->%d) = %f \n",
	   z8, t8, av8, ifilt_rest, ifilt_obs, kcor8 );
  }


  exit(1);
  
}  // end of test_fortran

// ******************************
void test_zcmb_dLmag_invert(void) {

  char fnam[] = "test_zcmb_dLmag_invert" ;
  double MU, zCMB;
  // ----------- RETURN ------------
  for(MU=32.0; MU < 49.0; MU+=1.0 ) {
    zCMB = zcmb_dLmag_invert(MU, &INPUTS.HzFUN_INFO);
  }
  debugexit(fnam);
  return ;
} // end test_zcmb_dLmag_invert

// ******************************
void test_PARSE_WORDS(void) {

  char fnam[] = "test_PARSE_WORDS" ;
  char dumString[] = "I walk to the beach";
  char dumFile[]   = "SIMGEN_SALT2_SDSS_FITS.input" ;
  //char dumFile[]   = "DES_SVA2_Y2.HOSTLIB" ;
  char tmpWord[80];
  int NWD, iwd;

  store_PARSE_WORDS(-1,"");

  NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE,dumFile);

  for(iwd=0; iwd < NWD; iwd++ ) {
    get_PARSE_WORD(0,iwd,tmpWord);
    printf(" xxx word[%2d] = '%s' \n", iwd, tmpWord);
  }

  NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,dumString);
  for(iwd=0; iwd < NWD; iwd++ ) {
    get_PARSE_WORD(0,iwd,tmpWord);
    printf(" xxx word[%2d] = '%s' \n", iwd, tmpWord);
  }

  debugexit(fnam);

  return ;
} // end test_PARSE_WORDS


// *********************
void test_ran(void) {
  int i, NTMP=0 ;
  double x, y; 
  double x0 = 0.0 ;
  // ------------- BEGIN --------
  for ( i = 1; i<=100000; i++ ) {
    // now smear with sigma=1
    NTMP++ ;
    if ( NTMP==100 ) {  fill_RANLISTs();  NTMP = 0 ; }
    x = getRan_Gauss(1);
    if ( x >=  0.0 ) 
      { y = 0.5 + GaussIntegral(x0,x);     }
    else 
      { y = 0.5 - GaussIntegral(x0,fabs(x) );     }
    printf(" 7777777 %f %f \n", x, y );
  }
  exit(1);
}  // end of test_ran

// ***********************
void test_igm(void) {

  int ilam;
  double lrest, lobs, z, z_a, z_b, tau_a, tau_b;
  char PATH_IGM_PARAM[MXPATHLEN];
  //  char fnam[] = "test_igm";

  // --------------- BEGIN --------------

  print_banner ( " Initialize igm correction for Spectra " );

  sprintf(PATH_IGM_PARAM,"%s/models/inoue_igm", PATH_SNDATA_ROOT);
  sprintf(DLA_FILE,"%s/DLAcoeff.txt", PATH_IGM_PARAM);
  sprintf(LAF_FILE,"%s/LAFcoeff.txt", PATH_IGM_PARAM);

  read_Inoue_coeffs();

  z_a = 2.0;
  z_b = 5.0;

  for (ilam=500; ilam <=6000; ilam+=500 ) {
    lrest = (double)ilam ;
    lobs  = lrest*(1.+z);

    z     = z_a ;
    tau_a = tLSLAF(z,lobs) + tLCLAF(z,lobs) + tLSDLA(z,lobs) + tLCDLA(z,lobs); 

    z     = z_b ;
    tau_b = tLSLAF(z,lobs) + tLCLAF(z,lobs) + tLSDLA(z,lobs) + tLCDLA(z,lobs); 

    printf(" xxx lrest=%6.0f tau(z=%.1f)=%le   tau(z=%.1f)=%le\n", 
	   lrest, z_a, tau_a,  z_b, tau_b) ;
    fflush(stdout);
  }
  
  exit(1);

} // end test_igm 
