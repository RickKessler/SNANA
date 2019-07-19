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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>

#include <sys/stat.h>

#include "fitsio.h"
#include "sntools.h"
#include "MWgaldust.h"
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
#include "sntools_fitsio.h"
#include "sntools_trigger.h" 
#include "sntools_grid.h"
#include "sntools_spectrograph.h"
#include "sntools_genGauss_asym.h"
#include "sntools_output.h"   // added Jan 11 2017
#include "inoue_igm.h"        // added Jun 27 2019

#include "genmag_ALL.h"
#include "MWgaldust.h" 

// include C code
#include "SNcadenceFoM.c"

#define SNGRIDGEN


// ******************************************
int main(int argc, char **argv) {

  int ilc, istat, i, isp, ifilt ;

  // define local structures
  SIMFILE_AUX_DEF  SIMFILE_AUX ;
  char fnam[] = "main"; 

  // ------------- BEGIN --------------

  sprintf(BANNER,"Begin execution of snlc_sim.exe  " );
  print_banner(BANNER);

  //  errmsg(SEV_FATAL, 0, "main", "testing CodeTest", "Remove this" ); 

  t_start = time(NULL);

  parse_commandLine_simargs(argc, argv);

  // init PATH_SNDATA_ROOT
  init_SNPATH();

  // init fortran variables
  istat = 0;
  init_snvar__(&istat); 

  // init sim-variables
  init_simvar();

  //  test_igm(); // xxxx

  // read user input file for directions
  if ( get_user_input() != SUCCESS ) { madend(1) ; }

  // init random number generator, and store first random.
  if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID  ) 
    { init_simRandoms();  }

  // prepare user input after init_simRandoms to allow 
  // random systematic shifts.
  prep_user_input();

  // check for random CID option (after randoms are inited above)
  init_CIDRAN();

  // init based on GENSOURCE
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENRANDOM ) {
    init_RANDOMsource();
    init_RateModel();

    if ( INPUTS.INIT_ONLY == 1 ) { debugexit("main: QUIT AFTER RATE-INIT"); }
  }
  else if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
#ifdef SNGRIDGEN
    init_GRIDsource(0); 
#endif
  }  
  else {
    sprintf(c1err,"%s is invalid GENSOURCE", INPUTS.GENSOURCE );
    errmsg(SEV_FATAL, 0, "main", c1err, "" ); 
  }


  // prepare randome systematic shifts after reading SURVEY from SIMLIB
  prep_RANSYSTPAR() ;

  // abort on NGEN=0 after printing N per season (init_DNDZ_Rate above)
  if ( INPUTS.NGEN_LC <= 0 && INPUTS.NGENTOT_LC <= 0 ) {
    sprintf(c1err,"NGEN_LC=0 & NGENTOT_LC=0" );
    errmsg(SEV_FATAL, 0, "main", c1err, ""); 
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
  INIT_FLUXERRMODEL(INPUTS.FLUXERRMODEL_OPTMASK,  INPUTS.FLUXERRMODEL_FILE,
		    INPUTS.FLUXERRMAP_IGNORE_DATAERR);


  // init anomalous host-subtraction noise
  INIT_NOISEMODEL_HOST_LEGACY(INPUTS.HOSTNOISE_FILE);

  // init WRONGHOST model
  double zMIN = INPUTS.GENRANGE_REDSHIFT[0] - 0.01 ;
  double zMAX = INPUTS.GENRANGE_REDSHIFT[1] + 0.01 ;
  INIT_WRONGHOST(INPUTS.WRONGHOST_FILE, zMIN, zMAX);

  // init user-specified z-dependence of sim parameters
  init_zvariation();

  // xxx mark delete May 27 2019:  init_peakmjdList();

  // initialize model that generates magnitudes.
  init_kcor(INPUTS.KCOR_FILE);

  init_genmodel();
  init_modelSmear(); 

  init_genSpec();     // July 2016: prepare optional spectra

  if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID )  { 
    // init search efficiency
    init_SEARCHEFF(GENLC.SURVEY_NAME,INPUTS.APPLY_SEARCHEFF_OPT); 
  } 
 
  // create/init output sim-files
  init_simFiles(&SIMFILE_AUX);

  print_banner( " Begin Generating Lightcurves. " );
  fflush(stdout);


  // check option to dump rest-frame mags
 SIMLIB_DUMP:
  if ( INPUTS.USEFLAG_DMPTREST ) { DUMP_GENMAG_DRIVER() ; }


  // check for simlib-dump
  SIMLIB_DUMP_DRIVER();


  if ( INPUTS.INIT_ONLY ==2 ) { debugexit("main: QUIT AFTER FULL INIT"); }

  // =================================================
  // start main loop over "ilc"


  for ( ilc = 1; ilc <= INPUTS.NGEN ; ilc++ ) {

    NGENLC_TOT++;

    if ( INPUTS.TRACE_MAIN  ) { dmp_trace_main("01", ilc) ; }

    if ( fudge_SNR() == 2 ) { goto GETMAGS ; }

    init_GENLC();

    if ( INPUTS.TRACE_MAIN  ) { dmp_trace_main("02", ilc) ; }

    if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID ) 
      { init_RANLIST(); }      // init list of random numbers for each SN    

    gen_event_driver(ilc); 

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
    GENFLUX_DRIVER();  // July 2016

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
      LOAD_SEARCHEFF_DATA();
      GENLC.SEARCHEFF_MASK = 
	gen_SEARCHEFF(GENLC.CID                 // (I) ID for dump/abort
		      ,&GENLC.SEARCHEFF_SPEC     // (O)
		      ,&GENLC.SEARCHEFF_zHOST    // (O) Mar 2018
		      ,&GENLC.MJD_TRIGGER ) ;    // (O)
    }

    for ( i=1; i<=NLIST_RAN; i++ )  { RANLAST[i] = FlatRan1(i); }

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
    if ( GENLC.STOPGEN_FLAG == 1 )  { goto ENDLOOP; }
    
    fflush(stdout);
    
  } // end of ilc loop


 ENDLOOP:

  t_end = time(NULL);

  // print final statistics on generated lightcurves.

  simEnd(&SIMFILE_AUX);

  return(0);

} // end of main



// **************************
void parse_commandLine_simargs(int argc, char **argv) {

  // Created April 2012.
  // Store user command-line args.

  int i;
  char *inFile = INPUTS.INPUT_FILE_LIST[0] ;
  char fnam[] = "parse_commandLine_simargs";

  // ---------------- BEGIN --------------

  printf("   Full command: ");

  if ( argc >= 2 ) {
    sprintf(inFile, "%s", argv[1] );
    NARGV_LIST = argc ;

    for ( i = 0; i < NARGV_LIST ; i++ ) {
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

  return ;

} // end of parse_commandLine_simargs


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

  //end_skewNormal();  // Sep 2016

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

  /* xxxxxxx mark delete Jun 25 2019 xxxxxxxx
  double XN, XNUPD;
  XN    = (double)N ; 
  XNUPD = (double)INPUTS.NGEN_SCREEN_UPDATE ;
  if ( fmod( XN, XNUPD ) == 0 ) return 1 ;
  xxxxxxxxxxxx */


  return 0;
}

// ******************************************
int get_user_input(void) {

  /**********

  Get user input from four different priority levels:
  (higher priority number overrides value from lower number)
  
  1. hard-wired default values in set_user_defaults()
  2. read from defaults file in $SNDATA_ROOT/models/snlc_sim.defaults
  3. read from user-input file "input_file"
  3.5 read 2nd (optional) user-input file 
  4. read from command-line override

   (changed from rd_input on 12/6/2007)

  Feb 4, 2008: fix buggy error message when default_file does not
               exist.

  May 24, 2010: change size of default_file from 120 to 200

  Jun 16, 2010: add option to read 2nd input file : input_file_include

  Apr 02, 2012: before reading 2nd input-include file, check command-line
                args to allow override. 

  Jun 20 2016: remove reading of default snlc_sim.defaults file.

  Jun 18 2017: move prep_user_input into main, after ranSeed init
  Jul 30 2018: check input_file_include2

  ***********/
  FILE *fp_input ;  // name of user input file 
  int i ;
  char *input_file;
  char fnam[] = "get_user_input"    ;

  // ------------ BEGIN ---------------

  // set hard-wired default values

  set_user_defaults();

  // ------------------------------------------------


  read_input(INPUTS.INPUT_FILE_LIST[0]);

  // ------------------------------
  // check option to read 2nd input-include file
  // First check command-line args here for override of 2nd input-include...
  // This check is needed because the override list is read below after
  // the 2nd input file is alread read.
  for ( i = 2; i < NARGV_LIST ; i++ ) {   
    if ( strcmp(ARGV_LIST[i],"INPUT_FILE_INCLUDE") == 0 ) { 
      sprintf(INPUTS.INPUT_FILE_LIST[1], "%s", ARGV_LIST[i+1] ) ; 
      INPUTS.INPUT_FILE_LIST[2][0] = 0 ; // erase 2nd include file
    }
  } 

  // check 1st include file
  if ( strlen(INPUTS.INPUT_FILE_LIST[1]) > 1 ) 
    { read_input(INPUTS.INPUT_FILE_LIST[1] );  }

  // check 2nd include file
  if ( strlen(INPUTS.INPUT_FILE_LIST[2]) > 1 ) 
    { read_input(INPUTS.INPUT_FILE_LIST[2] );  }

  // --------------------------------------------
  // check for command line overrides
  // --------------------------------------------

  sim_input_override();
  check_argv();

  // check for "PERFECT" options
  genperfect_override();


  // -------------------------------------------
  // make a few checks, compute a few flags and print user input to screen
  // -------------------------------------------

  
  return SUCCESS;

}  // end of get_user_input



// *********************************
void set_user_defaults(void) {

  int ifilt;
  int i, iz;
  float x;
  double zero = 0.0 ;
  char fnam[] = "set_user_defaults" ;
  // --------------- BEGIN ---------------

  INPUTS.TRACE_MAIN = 0;
  INPUTS.DEBUG_FLAG = 0 ;
  INPUTS.OPT_DEVEL_BBC7D = 0 ;
  NLINE_RATE_INFO   = 0;

  // don't init zero'th input file since that is the main input file
  for(i=1; i < MXINPUT_FILE_SIM; i++ ) { INPUTS.INPUT_FILE_LIST[i][0] = 0 ; }
  INPUTS.NREAD_INPUT_FILE = 0 ;

  sprintf(INPUTS.COMMENT,"BLANK");

  // - - - - - - - 

  NPAR_ZVAR_USR = USE_ZVAR_FILE = 0;
  sprintf(INPUT_ZVARIATION_FILE,"BLANK");
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

  INPUTS.ISEED      = 1 ;
  INPUTS.RANLIST_START_GENSMEAR = 1 ;

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

  INPUTS.OMEGA_MATTER  =  0.3 ;
  INPUTS.OMEGA_LAMBDA  =  0.7 ;
  INPUTS.W0_LAMBDA     = -1.0 ;
  INPUTS.H0            = 70.0 ;

  INPUTS.GENRANGE_RA[0]   = -360.0  ;
  INPUTS.GENRANGE_RA[1]   = +360.0  ;
  INPUTS.GENRANGE_DEC[0]  = -360.0  ;
  INPUTS.GENRANGE_DEC[1]  = +360.0  ;

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

  sprintf(INPUTS.GENPAR_SELECT_FILE,"BLANK");

  sprintf(INPUTS.STRETCH_TEMPLATE_FILE,"BLANK");

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SHAPEPAR, zero ); 

  INPUTS.DO_AV          = 0 ;
  INPUTS.GENRANGE_AV[0] = 0.0 ;
  INPUTS.GENRANGE_AV[1] = 0.0 ;
  INPUTS.GENEXPTAU_AV   = 0.0 ;
  INPUTS.GENGAUSIG_AV   = 0.0 ;
  INPUTS.GENGAUPEAK_AV  = 0.0 ;
  INPUTS.GENRATIO_AV0   = 0.0 ;
  INPUTS.WV07_GENAV_FLAG    =  0;
  INPUTS.WV07_REWGT_EXPAV   = -9.0;

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_RV, zero );
  INPUTS.GENGAUSS_RV.RANGE[0] = 2.1 ;
  INPUTS.GENGAUSS_RV.RANGE[1] = 4.1 ;
  INPUTS.GENGAUSS_RV.PEAK     = RV_MWDUST ; // for SN host

  // init SALT2 gen ranges
  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2c, zero );
  INPUTS.GENGAUSS_SALT2c.PEAK     =  0.0 ;
  INPUTS.GENGAUSS_SALT2c.SIGMA[0] =  1.0 ;
  INPUTS.GENGAUSS_SALT2c.SIGMA[1] =  1.0 ;
  INPUTS.GENGAUSS_SALT2c.RANGE[0] = -0.5 ;
  INPUTS.GENGAUSS_SALT2c.RANGE[1] =  1.0 ;

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2x1, zero );
  INPUTS.GENGAUSS_SALT2x1.PEAK     =  0.0 ;
  INPUTS.GENGAUSS_SALT2x1.SIGMA[0] =  5.0 ;
  INPUTS.GENGAUSS_SALT2x1.SIGMA[1] =  5.0 ;
  INPUTS.GENGAUSS_SALT2x1.RANGE[0] = -5.0 ;
  INPUTS.GENGAUSS_SALT2x1.RANGE[1] =  5.0 ;

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2ALPHA, zero );
  INPUTS.GENGAUSS_SALT2ALPHA.PEAK     =  0.11 ;
  INPUTS.GENGAUSS_SALT2ALPHA.RANGE[0] =  0.001; 
  INPUTS.GENGAUSS_SALT2ALPHA.RANGE[1] =  0.40 ; 

  init_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2BETA, zero );
  INPUTS.GENGAUSS_SALT2BETA.PEAK     =  3.20 ;
  INPUTS.GENGAUSS_SALT2BETA.RANGE[0] =  0.5 ;
  INPUTS.GENGAUSS_SALT2BETA.RANGE[1] =  9.9 ;

  INPUTS.BIASCOR_SALT2GAMMA_GRID[0] = +9.0 ; // min
  INPUTS.BIASCOR_SALT2GAMMA_GRID[1] = -9.0 ; // max

  INPUTS.SALT2BETA_cPOLY[0] = 0.0 ;
  INPUTS.SALT2BETA_cPOLY[1] = 0.0 ;
  INPUTS.SALT2BETA_cPOLY[2] = 0.0 ;

  INPUTS.SALT2mu_FILE[0] = 0 ;  // May 2013

  INPUTS.GENALPHA_SALT2     =  0.0 ; // legacy variable
  INPUTS.GENBETA_SALT2      =  0.0 ; // legacy variable
  INPUTS.LEGACY_colorXTMW_SALT2 = 0 ;

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

  INPUTS.KCORFLAG_STRETCH   = 0 ;
  INPUTS.KCORFLAG_COLOR     = GENFRAME_OBS ; // default => use obs-frame closest bands

  INPUTS.GENMAG_OFF_GLOBAL  = 0.0 ;
  INPUTS.GENMAG_OFF_NON1A   = 0.0 ;

  for ( ifilt=0 ; ifilt < MXFILTINDX; ifilt++ ) {
    INPUTS.GENMAG_OFF_MODEL[ifilt]  = 0.0 ;
    INPUTS.GENMAG_OFF_ZP[ifilt]     = 0.0 ;

    INPUTS.TMPOFF_ZP[ifilt]        =  0.0 ;
    INPUTS.TMPOFF_MODEL[ifilt]     =  0.0 ;
    INPUTS.TMPOFF_LAMSHIFT[ifilt]  = 0.0 ;

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
  INPUTS.GENMODEL_ERRSCALE     = 0.00 ; // .001 -> 0 (Jun 20 2016) 
  INPUTS.GENMODEL_ERRSCALE_OPT = 1;   // use peak error at all epochs
  INPUTS.GENMODEL_ERRSCALE_CORRELATION = 0.0;   // corr with GENMAG_SMEAR
  INPUTS.GENMODEL_MSKOPT             = 0 ; 
  INPUTS.GENMODEL_ARGLIST[0]         = 0 ;
  INPUTS.GENMAG_SMEAR[0]             = 0.0;
  INPUTS.GENMAG_SMEAR[1]             = 0.0;
  INPUTS.GENSMEAR_RANGauss_FIX       = -999.0 ;
  INPUTS.GENSMEAR_RANFlat_FIX        = -999.0 ;

  INPUTS.GENMAG_SMEAR_SCALE = 1.0;
  sprintf(INPUTS.GENMAG_SMEAR_MODELNAME, "NONE") ;
  sprintf(INPUTS.GENMAG_SMEAR_MODELARG,  "") ;

  sprintf(INPUTS.STRONGLENS_FILE,       "NONE");
  sprintf(INPUTS.WEAKLENS_PROBMAP_FILE, "NONE");
  INPUTS.WEAKLENS_DMUSCALE = 1.0 ;
  INPUTS.WEAKLENS_DSIGMADZ = 0.0 ;

  INPUTS.NPAR_GENSMEAR_USRFUN     = 0 ;
  for (i=0; i < 100; i++ ) 
    { INPUTS.GENMAG_SMEAR_USRFUN[i] = -9.0 ; }

  INPUTS.DO_MODELSMEAR   =  0 ;
  NSMEARPAR_OVERRIDE     =  0 ; // see sntools_genSmear.h

  INPUTS.GENRANGE_TYPE[0] = 0;
  INPUTS.GENRANGE_TYPE[1] = 0;
  INPUTS.GENRANGE_CID[0]  = 0;
  INPUTS.GENRANGE_CID[1]  = 0;

  INPUTS.GENFILTERS[0] = 0 ;
  INPUTS.NFILTDEF_OBS = 0;

  INPUTS.GENRANGE_DMPEVENT[0]  = 0 ;
  INPUTS.GENRANGE_DMPEVENT[1]  = 0 ;

  INPUTS.GENRANGE_DMPTREST[0]  = 9999. ;
  INPUTS.GENRANGE_DMPTREST[1]  = 9999. ;
  INPUTS.USEFLAG_DMPTREST = 0 ;

  INPUTS.FORMAT_MASK      = 32 ;                   // 2=TEXT  32=FITS
  INPUTS.WRITE_MASK       = WRITE_MASK_SIM_SNANA ; // default
  INPUTS.WRFLAG_MODELPAR  = 1; // default is yes

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
  sprintf(INPUTS.GENVERSION,      "%s", "BLANK" );  
  sprintf(INPUTS.GENPREFIX,       "%s", "BLANK" );
  sprintf(INPUTS.GENSOURCE,       "%s", "RANDOM" );
  sprintf(INPUTS.GENMODEL,        "%s", "BLANK" );
  INPUTS.GENMODEL_EXTRAP_LATETIME[0] = 0 ;
  sprintf(INPUTS.KCOR_FILE,   "%s",  "BLANK" );

  sprintf(INPUTS.GENSNXT, "%s", "CCM89" );
  OPT_SNXT = OPT_SNXT_CCM89;

  GENLC.IFLAG_GENSOURCE = 0;

  sprintf(INPUTS.SIMLIB_FILE,      "%s", "BLANK" );
  sprintf(INPUTS.SIMLIB_FIELDLIST, "%s", "ALL" );
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

  INPUTS.SIMLIB_MSKOPT = 0 ;
  GENLC.SIMLIB_IDLOCK  = -9;

  sprintf(INPUTS.HOSTLIB_FILE,          "NONE" );  // input library
  sprintf(INPUTS.HOSTLIB_WGTMAP_FILE,   "NONE" );  // optional wgtmap
  sprintf(INPUTS.HOSTLIB_ZPHOTEFF_FILE, "NONE" );  // optional zphot-eff
  sprintf(INPUTS.HOSTLIB_SPECBASIS_FILE, "NONE" );  // optional host-spec templates
  INPUTS.HOSTLIB_STOREPAR_LIST[0] = 0 ; // optional vars -> outfile

  INPUTS.HOSTLIB_USE    = 0;
  INPUTS.HOSTLIB_MSKOPT = 0;
  INPUTS.HOSTLIB_MAXREAD     = 1000000000 ; // default is 1 billion
  INPUTS.HOSTLIB_MXINTFLUX_SNPOS = 0.99 ;  // use 99% of total flux for SNPOS
  INPUTS.HOSTLIB_GALID_NULL      = -9;     // value for no host
  INPUTS.HOSTLIB_GALID_PRIORITY[0] = 0 ;
  INPUTS.HOSTLIB_GALID_PRIORITY[1] = 0 ;
  INPUTS.HOSTLIB_SBRADIUS          = 0.6*2.0 ; // arcsec
  INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL = 9999999;  // default is never re-use host
  INPUTS.HOSTLIB_SERSIC_SCALE      = 1.0 ;

  INPUTS.HOSTLIB_GENZPHOT_OUTLIER[0] = -9.0 ;
  INPUTS.HOSTLIB_GENZPHOT_OUTLIER[1] = -9.0 ;
  for(i=0; i<5; i++) { 
    INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[i] =  0.0 ; 
    INPUTS.HOSTLIB_GENZPHOT_BIAS[i]     =  0.0 ; 
  }
  INPUTS.USE_HOSTLIB_GENZPHOT = 0 ; // logical flag

  // Nov 23 2015
  // define polynom function of ztrue for zSN-zGAL tolerance.
  // Default is close to what we had before, which is quite big
  INPUTS.HOSTLIB_DZTOL[0] = 0.002 ;
  INPUTS.HOSTLIB_DZTOL[1] = 0.040 ;
  INPUTS.HOSTLIB_DZTOL[2] = 0.0   ;

  // debug options
  INPUTS.HOSTLIB_GALID_FORCE   = -9;
  INPUTS.HOSTLIB_FIXRAN_RADIUS = -9;
  INPUTS.HOSTLIB_FIXRAN_PHI    = -9;

  INPUTS.FLUXERRMODEL_OPTMASK=0 ;
  sprintf(INPUTS.FLUXERRMODEL_FILE,          "NONE" ); 
  sprintf(INPUTS.FLUXERRMAP_IGNORE_DATAERR,  "NONE");
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

  INPUTS.FUDGESCALE_PSF       = 1.0 ;
  INPUTS.FUDGESCALE_SKYNOISE  = 1.0 ;
  INPUTS.FUDGESCALE_READNOISE = 1.0 ;
  INPUTS.FUDGESHIFT_ZPT       = 0.0 ;
  INPUTS.FUDGESCALE_FLUXERR   = 1.0 ;
  INPUTS.FUDGESCALE_FLUXERR2  = 1.0 ;
  INPUTS.FUDGEOPT_FLUXERR     = 0 ;
  INPUTS.FUDGE_MAGERR         = 0.0 ;

  for(ifilt=0; ifilt<MXFILTINDX; ifilt++)  { 
    INPUTS.FUDGE_MAGERR_FILTER[ifilt]        = 0.0; 
    INPUTS.FUDGESHIFT_ZPT_FILTER[ifilt]      = 0.0 ;
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
  sprintf(INPUTS_SEARCHEFF.USER_SPEC_FILE, "NONE");
  sprintf(INPUTS_SEARCHEFF.USER_zHOST_FILE,"NONE");

  sprintf(INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE,  "DEFAULT" ); 
  sprintf(INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE,    "DEFAULT" ); 
 
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
  INPUTS.CLEARPROMPT = 0;

  INPUTS.NVAR_SIMGEN_DUMP = -9 ; // note that 0 => list variables & quit
  INPUTS.IFLAG_SIMGEN_DUMPALL = 0 ; // dump only SN written to data file.
  INPUTS.PRESCALE_SIMGEN_DUMP = 1 ; // prescale

#ifdef SNGRIDGEN
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

  // set default spectrograph options
  INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK         =  0 ;
  INPUTS.SPECTROGRAPH_OPTIONS.DOFLAG_SPEC     =  0 ;
  INPUTS.SPECTROGRAPH_OPTIONS.NLAMSIGMA       =  5.0;
  INPUTS.SPECTROGRAPH_OPTIONS.SCALE_LAMSIGMA  =  1. ;
  INPUTS.SPECTROGRAPH_OPTIONS.SCALE_SNR       =  1. ;  // scale on SNR
  INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE   =  1. ;  // scale Texpose
  INPUTS.SPECTROGRAPH_OPTIONS.ILAM_SPIKE      = -9 ;   // lambda bin

  NPEREVT_TAKE_SPECTRUM =  0 ; // Mar 14 2017
  INPUTS.TAKE_SPECTRUM_DUMPCID  = -9 ;
  INPUTS.TAKE_SPECTRUM_HOSTFRAC =  0.0 ;
  INPUTS.TAKE_SPECTRUM_TEMPLATE_TEXPOSE_SCALE =  1.0 ;
  INPUTS.NWARP_TAKE_SPECTRUM = 0 ;
  INPUTS.NHOST_TAKE_SPECTRUM = 0 ;

} // end set_user_defaults_SPECTROGRAPH

// *******************************************
void set_user_defaults_RANSYSTPAR(void) {

  int ifilt; 

  INPUTS.RANSYSTPAR.USE = 0 ;

  // coherent (all bands) scale of flux-errors
  INPUTS.RANSYSTPAR.FUDGESCALE_FLUXERR  = 1.0 ; // true & measured 
  INPUTS.RANSYSTPAR.FUDGESCALE_FLUXERR2 = 1.0 ; // measured only
  INPUTS.RANSYSTPAR.FUDGESCALE_MWEBV    = 1.0 ;
  INPUTS.RANSYSTPAR.FUDGESHIFT_MWRV     = 0.0 ;

  for(ifilt=0; ifilt < MXFILTINDX; ifilt++ )  { 
    INPUTS.RANSYSTPAR.GENMAG_OFF_ZP[ifilt]   = 0.0 ; 
    INPUTS.RANSYSTPAR.FILTER_LAMSHIFT[ifilt] = 0.0 ; 
  }

  return ;

} // end set_user_defaults_RANSYSTPAR


// ***********************************
int read_input(char *input_file) {

  // Jan 2016: skip lines with '#' by reading entire line after each '#'
  // Mar 1 2017: read AV,RV stuff for SIMSED model. See READ_AVRV jump.
  // 
  // Jul 30 2018: pass input file as argument instead of fp
  // Apr 09 2019:
  //  + use new uniqueMatch() util to check for UNIQUE string-match
  //    and abort on duplicate key.
  //  + use keyMatch() util to check for match where the same
  //    key is allowed multiple times.
  //
  // Jun 26 2019: 
  //  + for NON1ASED, allow comma-separated sigmas for MAGSMEAR column

  FILE *fp;

  char 
    c_get[80], ctmp[80], ctmp2[80], tmpLine[200]
    ,ckey[60], cval[40], cpoly[60], warp_spectrum_string[100]
    ,*parname, *modelName   
    ,*ptrTmp = tmpLine
    ,*ptrComment = INPUTS.COMMENT
    ;

  float ftmp, sigTmp[2];
  int L_NON1ASED, L_NON1AKEY, L_PEC1AKEY, L_TMP, ITMP, NWD ;
  int itmp, N, NTMP, j, i, ifilt, ifilt_obs, opt_tmp;
  int key, NKEY, ovp1, ovp2, ITYPE ;
  int ifilt_list[MXFILTINDX];
  int iArg = -9;
  char  stringSource[] = "sim-input file" ;
  char  comma[] = ","; 
  char  fnam[] = "read_input"  ;

  // -------------- BEGIN -------------

  ENVreplace(input_file,fnam,1);
  if ( (fp = fopen(input_file, "rt"))==NULL ) {       
    sprintf ( c1err, "Cannot open input file :" );
    sprintf ( c2err," '%s' ", input_file);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  INPUTS.NREAD_INPUT_FILE++ ;
  printf(" --------------------------------------------------------\n");
  uniqueMatch("INIT",stringSource); // 4.2019
  printf("  Read user input file %d: %s \n", 
	 INPUTS.NREAD_INPUT_FILE, input_file );

  warp_spectrum_string[0] = ckey[0] = cval[0] = tmpLine[0] = 0 ;


  NWD = 0 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    NWD++ ;

    // if comment key is found, read remainder of line into dummy string
    // so that anything after comment key is ignored (even a valid key)
    // xxx mark delete if ( c_get[0] == '#' || c_get[0] == '!' || c_get[0] == '%' ) 
    if ( commentchar(c_get) ) 
      { ptrTmp = fgets(tmpLine, 80, fp) ; continue ; }

    if ( strcmp(c_get,"INPUT_FILE_INCLUDE:") == 0 ) {
      if ( strlen(INPUTS.INPUT_FILE_LIST[1]) == 0 )  // 1st include file
	{ readchar ( fp, INPUTS.INPUT_FILE_LIST[1]); }
      else if ( strlen(INPUTS.INPUT_FILE_LIST[2])==0 )  // 2nd include file
	{ readchar ( fp, INPUTS.INPUT_FILE_LIST[2] ); } // Jul 30 2018
      else {
	sprintf(c1err,"Cannot specify 3 INPUT_FILE_INCLUDE keys");
	sprintf(c2err,"Only 1 or 2 allowed.") ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }

    if ( uniqueMatch(c_get,"TRACE_MAIN:") ) 
      { readint ( fp, 1, &INPUTS.TRACE_MAIN ); continue; }

    if ( uniqueMatch(c_get,"DEBUG_FLAG:")  ) 
      { readint ( fp, 1, &INPUTS.DEBUG_FLAG );  continue ; }

    if ( uniqueMatch(c_get,"OPT_DEVEL_BBC7D:")  ) 
      { readint ( fp, 1, &INPUTS.OPT_DEVEL_BBC7D );  continue ; }

    // --- HOSTLIB stuff
    if ( uniqueMatch(c_get,"HOSTLIB_FILE:")   ) {
      readchar ( fp, INPUTS.HOSTLIB_FILE );
      setbit_HOSTLIB_MSKOPT(HOSTLIB_MSKOPT_USE) ;
      INPUTS.HOSTLIB_USE = 1;    // logical 
      continue ; 
    }

    if ( uniqueMatch(c_get,"HOSTLIB_WGTMAP_FILE:")  )
      {  readchar ( fp, INPUTS.HOSTLIB_WGTMAP_FILE ); continue ; }
    
    if ( uniqueMatch(c_get,"HOSTLIB_ZPHOTEFF_FILE:")  )
      { readchar ( fp, INPUTS.HOSTLIB_ZPHOTEFF_FILE ); continue ; }

    if ( uniqueMatch(c_get,"HOSTLIB_SPECBASIS_FILE:")  )
      { readchar ( fp, INPUTS.HOSTLIB_SPECBASIS_FILE ); continue ; }
    
    if ( uniqueMatch(c_get,"HOSTLIB_MSKOPT:")  ) {
      readint ( fp, 1, &itmp );
      // remove USE-bit if set
      if ( itmp & HOSTLIB_MSKOPT_USE ) { itmp -= HOSTLIB_MSKOPT_USE; }
      INPUTS.HOSTLIB_MSKOPT += itmp ;
      continue ; 
    }

    if ( uniqueMatch(c_get,"HOSTLIB_GENZPHOT_FUDGEPAR:")  ) {
      readfloat ( fp, 4, INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR );
      readchar(fp, cval);
      parse_input_GENZPHOT_OUTLIER(cval);
      continue ; 
    }
    if ( uniqueMatch(c_get,"HOSTLIB_GENZPHOT_BIAS:")  ) 
      { readfloat ( fp, 4, INPUTS.HOSTLIB_GENZPHOT_BIAS ); continue ; }

    if ( uniqueMatch(c_get,"HOSTLIB_DZTOL:")  ) 
      { readdouble ( fp, 3, INPUTS.HOSTLIB_DZTOL );  continue ; }

    if ( uniqueMatch(c_get,"HOSTLIB_SERSIC_SCALE:")  )
      { readdouble ( fp, 1, &INPUTS.HOSTLIB_SERSIC_SCALE ); continue ; }
   

    // for hostlib variables, allow STOREPAR or STOREVA
    if ( uniqueMatch(c_get,"HOSTLIB_STOREVAR:")  )
      { readchar ( fp, INPUTS.HOSTLIB_STOREPAR_LIST ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_STOREPAR:")  )
      { readchar ( fp, INPUTS.HOSTLIB_STOREPAR_LIST ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_MAXREAD:")  )
      { readint ( fp, 1, &INPUTS.HOSTLIB_MAXREAD ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_GALID_NULL:")  )
      { readint ( fp, 1, &INPUTS.HOSTLIB_GALID_NULL ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_GALID_PRIORITY:")  ) 
      { readint ( fp, 2, INPUTS.HOSTLIB_GALID_PRIORITY ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_MINDAYSEP_SAMEGAL:")  ) 
      { readint ( fp, 1, &INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_MXINTFLUX_SNPOS:")  ) 
      { readfloat ( fp, 1, &INPUTS.HOSTLIB_MXINTFLUX_SNPOS ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_GENRANGE_NSIGZ:")  )
      { readfloat ( fp, 2, INPUTS.HOSTLIB_GENRANGE_NSIGZ ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_GENRANGE_RA:")  ) 
      { readdouble ( fp, 2, INPUTS.HOSTLIB_GENRANGE_RA ); continue ; }
    
    if ( uniqueMatch(c_get,"HOSTLIB_GENRANGE_DEC:")  )
      { readdouble ( fp, 2, INPUTS.HOSTLIB_GENRANGE_DEC ); continue ; }
   
    if ( uniqueMatch(c_get,"HOSTLIB_GENRANGE_DECL:")  )
      { readdouble ( fp, 2, INPUTS.HOSTLIB_GENRANGE_DEC ); continue ; }
    
    if ( uniqueMatch(c_get,"HOSTLIB_SBRADIUS:")  )
      { readdouble ( fp, 1, &INPUTS.HOSTLIB_SBRADIUS ); continue ; }
    
    if ( uniqueMatch(c_get,"HOSTLIB_GALID_FORCE:")  )
      { readint( fp, 1, &INPUTS.HOSTLIB_GALID_FORCE ); continue ; }
    
    if ( uniqueMatch(c_get,"HOSTLIB_FIXRAN_RADIUS:")  )
      { readdouble( fp, 1, &INPUTS.HOSTLIB_FIXRAN_RADIUS); continue ; }
    
    if ( uniqueMatch(c_get,"HOSTLIB_FIXRAN_PHI:")  )
      { readdouble( fp, 1, &INPUTS.HOSTLIB_FIXRAN_PHI); continue ; }
   
    if ( uniqueMatch(c_get,"FLUXERRMODEL_FILE:")   )
      { readchar ( fp, INPUTS.FLUXERRMODEL_FILE ); continue ; }
    
    if ( uniqueMatch(c_get,"FLUXERRMODEL_OPTMASK:")   )
      { readint ( fp, 1, &INPUTS.FLUXERRMODEL_OPTMASK ); continue ; }
    
    if ( uniqueMatch(c_get,"FLUXERRMAP_IGNORE_DATAERR:")  )
      { readchar ( fp, INPUTS.FLUXERRMAP_IGNORE_DATAERR ); continue ; }   

    // --- anomalous host-subtraction noise ------
    if ( uniqueMatch(c_get,"HOSTNOISE_FILE:")   ) {
      readchar ( fp, INPUTS.HOSTNOISE_FILE );
      if ( INPUTS.HOSTLIB_USE == 0  ){
	sprintf(c1err,"Cannot define HOSTNOISE_FILE ");
	sprintf(c2err,"without first defining HOSTLIB_FILE");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      continue ; 
    }

    // -------
    // check options for z-variation; either separate file
    // or ZPOLY key(s) in sim-input file.
    
    if ( uniqueMatch(c_get,"ZVARIATION_FILE:")  ) {  
      readchar ( fp, INPUT_ZVARIATION_FILE ); 
      if ( !IGNOREFILE(INPUT_ZVARIATION_FILE) ) { USE_ZVAR_FILE = 1 ; }
      continue ; 
    }

    if ( keyMatch(c_get,"ZVARIATION_POLY:")  ) { 

      N = NPAR_ZVAR_USR ;
      INPUT_ZVARIATION[N].FLAG = FLAG_ZPOLY_ZVAR ;
      INPUT_ZVARIATION[N].NZBIN = 0 ;
      NPAR_ZVAR_USR++ ;  

      readchar(fp, INPUT_ZVARIATION[N].PARNAME ) ;
      readchar(fp, tmpLine);  
      // if there is a comma, use new GENPOLY struct; else ready legacy format
      if ( strstr(tmpLine,comma) ) {
	// store comma-separate poly coefficients in cpoly
	sprintf(cpoly, "%s", tmpLine);	
      }
      else {
	// read LEGACY format of space-separated 3rd order poly
	double zpoly[POLYORDER_ZVAR+2];
	sscanf(tmpLine, "%le", &zpoly[0]); 
	readdouble(fp, POLYORDER_ZVAR, &zpoly[1] );
	sprintf(cpoly,"%f,%f,%f,%f", zpoly[0],zpoly[1],zpoly[2],zpoly[3]);
      }
      parse_GENPOLY(cpoly, &INPUT_ZVARIATION[N].POLY, fnam);
      // xxx delete readdouble(fp,POLYORDER_ZVAR+1,INPUT_ZVARIATION[N].ZPOLY);
      continue ; 
    }

    // ------

    if ( uniqueMatch(c_get,"SNTYPE_Ia:")  ) {
      readint ( fp, 1, &ITYPE ) ;
      INPUTS.SNTYPE_Ia_SPEC = ITYPE ;
      INPUTS.SNTYPE_Ia_PHOT = ITYPE + OFFSET_TYPE_PHOT ;
      continue ; 
    }
    if ( uniqueMatch(c_get,"SNTYPES_Ia:")  ) {
      readint ( fp, 1, &INPUTS.SNTYPE_Ia_SPEC);
      readint ( fp, 1, &INPUTS.SNTYPE_Ia_PHOT);
      continue ; 
    }

    if ( uniqueMatch(c_get,"GENTYPE:")  ) {
      readint ( fp, 1, &ITYPE ) ;
      INPUTS.GENTYPE_SPEC = ITYPE ;
      INPUTS.GENTYPE_PHOT = ITYPE + OFFSET_TYPE_PHOT ;
      continue ;
    }
    if ( uniqueMatch(c_get,"GENTYPES:")  ) {
      readint ( fp, 1, &INPUTS.GENTYPE_SPEC); 
      readint ( fp, 1, &INPUTS.GENTYPE_PHOT);
      continue ;
    }

    if ( uniqueMatch(c_get,"COMMENT:")  )
      { ptrComment = fgets ( INPUTS.COMMENT, 80, fp) ; continue ; }   

    if ( uniqueMatch(c_get,"NONLINEARITY_FILE:")  ) 
      { readchar ( fp, INPUTS.NONLINEARITY_FILE ); continue ; }    

    if ( uniqueMatch(c_get,"SIMLIB_FILE:")  ) 
      { readchar ( fp, INPUTS.SIMLIB_FILE ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_FIELDLIST:")  ) {
      readchar ( fp, INPUTS.SIMLIB_FIELDLIST ); 
      INPUTS.SIMLIB_FIELDSKIP_FLAG = 1 ; // count skipped fields in NGENTOT
      continue ;
    }
    if ( uniqueMatch(c_get,"SIMLIB_IDSTART:")  )
      { readint ( fp, 1, &INPUTS.SIMLIB_IDSTART ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_MAXRANSTART:")  )
      { readint ( fp, 1, &INPUTS.SIMLIB_MAXRANSTART ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_IDLOCK:")  )
      { readint ( fp, 1, &INPUTS.SIMLIB_IDLOCK ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_MINOBS:")  ) 
      { readint ( fp, 1, &INPUTS.SIMLIB_MINOBS ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_MINSEASON:")  )
      { readdouble ( fp, 1, &INPUTS.SIMLIB_MINSEASON ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_DUMP:")  ) 
      { readint ( fp, 1, &INPUTS.SIMLIB_DUMP ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_NREPEAT:")  ) 
      { readint ( fp, 1, &INPUTS.SIMLIB_NREPEAT );  continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_NSKIPMJD:")  )
      { readint ( fp, 1, &INPUTS.SIMLIB_NSKIPMJD ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_IDSKIP:")  ) {
      N = INPUTS.NSKIP_SIMLIB ;
      readint ( fp, 1, &INPUTS.SIMLIB_IDSKIP[N] );
      INPUTS.NSKIP_SIMLIB++ ;
      continue ; 
    }
    if ( uniqueMatch(c_get,"SIMLIB_CADENCEFOM_ANGSEP:")  ) 
      { readfloat ( fp, 1, &INPUTS.SIMLIB_CADENCEFOM_ANGSEP ); continue ; }
    
    if ( uniqueMatch(c_get,"SIMLIB_CADENCEFOM_PARLIST:")  ) {
      readdouble(fp, NPAR_SNCADENCEFOM,INPUTS.SIMLIB_CADENCEFOM_PARLIST );
      continue ; 
    }
    if ( uniqueMatch(c_get,"USE_SIMLIB_REDSHIFT:")   ) { 
      readint ( fp, 1, &INPUTS.USE_SIMLIB_REDSHIFT ); 
      INPUTS.USE_SIMLIB_GENOPT=1;
      continue ; 
    }
    if ( uniqueMatch(c_get,"USE_SIMLIB_DISTANCE:")   ) { 
      readint ( fp, 1, &INPUTS.USE_SIMLIB_DISTANCE ); 
      INPUTS.USE_SIMLIB_GENOPT=1;
      continue ; 
    }
    if ( uniqueMatch(c_get,"USE_SIMLIB_PEAKMJD:")   )  { 
      readint ( fp, 1, &INPUTS.USE_SIMLIB_PEAKMJD ); 
      INPUTS.USE_SIMLIB_GENOPT=1;
      continue ;
    }

    if ( uniqueMatch(c_get,"SIMLIB_MSKOPT:")  ) {
      readint ( fp, 1, &ITMP );   INPUTS.SIMLIB_MSKOPT+=ITMP ; 
      continue ; 
    }
    
    //  - - - 

    if ( uniqueMatch(c_get,"NGEN_LC:")  ) 
      { readint ( fp, 1, &INPUTS.NGEN_LC ); continue ; }

    if ( uniqueMatch(c_get,"NGENTOT_LC:")  ) 
      { readint ( fp, 1, &INPUTS.NGENTOT_LC ); continue ; }    

    if ( uniqueMatch(c_get,"NGEN_SEASON:")  ) 
      { readfloat ( fp, 1, &INPUTS.NGEN_SEASON ); continue ; }    

    if ( uniqueMatch(c_get,"NGEN_SCALE:")  )
      { readfloat ( fp, 1, &INPUTS.NGEN_SCALE ); continue ; }

    if ( uniqueMatch(c_get,"NGEN_SCALE_NON1A:")  )
      { readfloat ( fp, 1, &INPUTS.NGEN_SCALE_NON1A ); continue ; }

    if ( uniqueMatch(c_get,"NSUBSAMPLE_MARK:")  ) 
      { readint ( fp, 1, &INPUTS.NSUBSAMPLE_MARK ); continue ; }

    if ( uniqueMatch(c_get,"CIDOFF:")  ) 
      { readint ( fp, 1, &INPUTS.CIDOFF ); continue ; }
    
    if ( uniqueMatch(c_get,"CIDRAN_MAX:")  ) 
      { readint ( fp, 1, &INPUTS.CIDRAN_MAX ); continue ; }
    
    if ( uniqueMatch(c_get,"CIDRAN_MIN:")  ) 
      { readint ( fp, 1, &INPUTS.CIDRAN_MIN ); continue ; }    

    if ( uniqueMatch(c_get,"FORMAT_MASK:")  ) 
      { readint ( fp, 1, &INPUTS.FORMAT_MASK ); continue ; }
    
    if ( uniqueMatch(c_get,"WRFLAG_MODELPAR:")  )
      { readint ( fp, 1, &INPUTS.WRFLAG_MODELPAR ); continue ; }
    
    // - - - -  -

    if ( uniqueMatch(c_get,"NPE_PIXEL_SATURATE:")  ) 
      { readint ( fp, 1, &INPUTS.NPE_PIXEL_SATURATE ); continue ; }
    
    if ( uniqueMatch(c_get,"PHOTFLAG_SATURATE:")  )
      { readint ( fp, 1, &INPUTS.PHOTFLAG_SATURATE ); continue ; }
    
    if ( uniqueMatch(c_get,"PHOTFLAG_SNRMAX:")  ) 
      { readint ( fp, 1, &INPUTS.PHOTFLAG_SNRMAX ); continue ; }
    
    if ( uniqueMatch(c_get,"PHOTFLAG_NEARPEAK:")  ) 
      { readint ( fp, 1, &INPUTS.PHOTFLAG_NEARPEAK ); continue ; }

    if ( uniqueMatch(c_get,"PHOTFLAG_DETECT:")  ) 
      { readint ( fp, 1, &INPUTS_SEARCHEFF.PHOTFLAG_DETECT ); continue ; }
    
    if ( uniqueMatch(c_get,"PHOTFLAG_TRIGGER:")  )
      { readint ( fp, 1, &INPUTS_SEARCHEFF.PHOTFLAG_TRIGGER ); continue ; }   

    // -----
    // check for SIMGEN_DUMP

    if ( strstr(c_get,"SIMGEN_DUMP") != NULL ) {
      parse_input_SIMGEN_DUMP(fp, &iArg, c_get ); 
    }

    // -----------

    if ( uniqueMatch(c_get,"GENVERSION:")   ) {
      readchar ( fp, INPUTS.GENVERSION );
      sprintf(INPUTS.GENPREFIX,"%s", INPUTS.GENVERSION); // default
      continue ; 
    }
    if ( uniqueMatch(c_get,"GENPREFIX:")   ) 
      { readchar ( fp, INPUTS.GENPREFIX ); continue ; }

    if ( uniqueMatch(c_get,"CLEARPROMPT:")  )
      { readint ( fp, 1, &INPUTS.CLEARPROMPT ); continue ; }

    if ( uniqueMatch(c_get,"GENSOURCE:")  )
      { readchar ( fp, INPUTS.GENSOURCE ); continue ; }

    if ( uniqueMatch(c_get,"GENMODEL_MSKOPT:")  ) 
      { readint ( fp, 1, &INPUTS.GENMODEL_MSKOPT ); continue ; }

    if ( uniqueMatch(c_get,"GENMODEL_ARGLIST:")  )
      { parse_input_GENMODEL_ARGLIST(fp,&iArg); continue ; }

    if ( uniqueMatch(c_get,"GENMODEL:") > 0  ) {

      parse_input_GENMODEL(fp,&iArg);  //4.18.2019
      continue ;

    } // end GENMODEL

    if ( uniqueMatch(c_get,"GENMODEL_EXTRAP_LATETIME:") ) 
      { readchar ( fp, INPUTS.GENMODEL_EXTRAP_LATETIME ); continue ; }

    if ( uniqueMatch(c_get,"PATH_SNDATA_SIM:")  )
      { readchar(fp, INPUTS.PATH_SNDATA_SIM ); continue ; }

    if ( uniqueMatch(c_get,"PATH_NON1ASED:")  )
      { readchar(fp, INPUTS.NON1ASED.PATH ); continue ; }
    
    if ( uniqueMatch(c_get,"PATH_NONIASED:")  )
      { readchar(fp, INPUTS.NON1ASED.PATH ); continue ; }
   
    L_NON1ASED = ( INPUTS.NON1A_MODELFLAG == MODEL_NON1ASED ) ;


    L_NON1AKEY = ( strcmp(c_get,"NON1A_KEYS:")    ==0 || 
		   strcmp(c_get,"NONIA_KEYS:")    ==0 ||
		   strcmp(c_get,"NON1ASED_KEYS:") ==0 || 
		   strcmp(c_get,"NONIASED_KEYS:") ==0      ) ;

    if ( L_NON1ASED && L_NON1AKEY ) {
      readint ( fp, 1, &NKEY );
      if ( NKEY < 2 ) {
	sprintf(c1err,"NON1A_KEYS = %d, but must be >=2.", NKEY);
	sprintf(c1err,"Check sim-input file." );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      INPUTS.NON1ASED.NKEY = NKEY;
      for ( key = 1; key <= NKEY; key++ ) {
	readchar(fp, ckey );
	sprintf(INPUTS.NON1ASED.KEYLIST[key], "%s", ckey );
      }

      // check required keys
      ovp1 = strcmp(INPUTS.NON1ASED.KEYLIST[1],"INDEX") ;
      ovp2 = strcmp(INPUTS.NON1ASED.KEYLIST[2],"WGT"  ) ;
      if ( ovp1 != 0 ) {
	sprintf(c1err,"First NON1ASED_KEY is '%s' , but must be 'INDEX' ",
		INPUTS.NON1ASED.KEYLIST[1] ) ;
	sprintf(c2err,"Check sim-input file.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
      if ( ovp2 != 0 ) {
	sprintf(c1err,"Second NON1ASED_KEY is ' %s' , but must be 'WGT' ", 
		INPUTS.NON1ASED.KEYLIST[2]) ;
	sprintf(c2err,"Check sim-input file.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

    } // end of NON1ASED_KEYS

    if ( strcmp(c_get,"NON1A_STOP:")    ==0 ) {  INPUTS.NON1ASED.STOP = 1 ; }
    if ( strcmp(c_get,"NON1ASED_STOP:") ==0 ) {  INPUTS.NON1ASED.STOP = 1 ; }
    L_NON1ASED  = (INPUTS.NON1A_MODELFLAG == MODEL_NON1ASED) ;
    L_NON1AKEY  = ( strcmp(c_get,"NON1A:")    ==0 || 
		    strcmp(c_get,"NONIA:")    ==0 ||
		    strcmp(c_get,"NONIASED:") ==0 ||
		    strcmp(c_get,"NON1ASED:") ==0 ) ;

    L_PEC1AKEY  = ( strcmp(c_get,"PEC1A:")    ==0 || 
		    strcmp(c_get,"PECIA:")    ==0 ||
		    strcmp(c_get,"PECIASED:") ==0 ||
		    strcmp(c_get,"PEC1ASED:") ==0 ) ;

    L_TMP = 
      (L_NON1ASED) && (L_NON1AKEY || L_PEC1AKEY) && 
      (INPUTS.NON1ASED.STOP==0) ;

    if ( L_TMP ) { 

      INPUTS.NON1ASED.NINDEX++ ; N = INPUTS.NON1ASED.NINDEX;
      if ( L_PEC1AKEY ) { 
	INPUTS.NON1ASED.NPEC1A++ ;
	INPUTS.NON1ASED.ISPEC1A[N] = 1; 
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
	readchar(fp, ctmp );

	INPUTS.NON1ASED.KEYVAL[N][key] = ftmp ;

	if ( strcmp(ckey, "INDEX") == 0 ) 
	  { sscanf(ctmp, "%d", &INPUTS.NON1ASED.INDEX[N]); }
	else if ( strcmp(ckey, "WGT") == 0 ) 
	  { sscanf(ctmp, "%f", &INPUTS.NON1ASED.WGT[N]); }
	else if ( strcmp(ckey, "MAGOFF") == 0 ) 
	  { sscanf(ctmp, "%f", &INPUTS.NON1ASED.MAGOFF[N]); }
	else if ( strcmp(ckey, "MAGSMEAR") == 0 ) {
	  split2floats(ctmp, comma, sigTmp ); 
	  INPUTS.NON1ASED.MAGSMEAR[N][0] = sigTmp[0] ;
	  INPUTS.NON1ASED.MAGSMEAR[N][1] = sigTmp[1] ;
	  // if non1a magsmear is nonzero, then set global GENMAG_SMEAR
	  // so that random numbers are generated for smearing.
	  if ( sigTmp[0] > 0.0  &&  INPUTS.GENMAG_SMEAR[0] == 0.0 ) 
	    { INPUTS.GENMAG_SMEAR[0] = INPUTS.GENMAG_SMEAR[1] = 0.001; }
	}
	else if ( strcmp(ckey, "SNTYPE") == 0 ) 
	  { sscanf(ctmp, "%d", &INPUTS.NON1ASED.SNTAG[N]); }
	else if ( strcmp(ckey, "SNTAG") == 0 )   // same as SNTYPE
	  { sscanf(ctmp, "%d", &INPUTS.NON1ASED.SNTAG[N]); }
	else {
	  sprintf(c1err, "Invalid NON1A key: '%s'", ckey );
	  sprintf(c2err, "Valid keys are INDEX WGT MAGOFF MAGSMEAR" );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
	}
      }

    }

 
    if ( uniqueMatch(c_get,"MJD_EXPLODE:")  ) {
      readdouble(fp, 1, &INPUTS.MJD_EXPLODE );
      if ( INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE < 0 ) 
	{ INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE = 0 ; }
      INPUTS.GENRANGE_PEAKMJD[0] = INPUTS.MJD_EXPLODE ;
      INPUTS.GENRANGE_PEAKMJD[1] = INPUTS.MJD_EXPLODE ;
      continue ;
    }
    if ( uniqueMatch(c_get,"OPTMASK_T0SHIFT_EXPLODE:")  ) {
      readint(fp, 1, &INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE); 
      continue ;
    }

    if ( uniqueMatch(c_get,"UVLAM_EXTRAPFLUX:")  ) {
      readdouble(fp, 1, &INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX); 
      continue ; 
    }

    if ( uniqueMatch(c_get,"MINSLOPE_EXTRAPMAG_LATE:")  )  {
      readdouble(fp, 1, &INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE); 
      continue ; 
    }
    
    if ( uniqueMatch(c_get,"RANSEED:")  ) { 
      readint ( fp, 1, &ITMP ); // read regular int
      INPUTS.ISEED = ITMP ;     // set unsigned int
      continue ; 
    }

    /*
    if ( uniqueMatch(c_get,"RANSEED_ADD:")  ) { 
      readint ( fp, 1, &ITMP ); // read regular int
      INPUTS.ISEED_ADD = ITMP ;     // set unsigned int
      } */

    if ( uniqueMatch(c_get,"RANLIST_START_GENSMEAR:") )
      { readint(fp, 1, &INPUTS.RANLIST_START_GENSMEAR ); continue ; }    

    if ( uniqueMatch(c_get,"GENRANGE_RA:")  ) 
      { readfloat ( fp, 2, INPUTS.GENRANGE_RA ); continue ; }
    
    if ( uniqueMatch(c_get,"GENRANGE_DECL:")  )
      { readfloat ( fp, 2, INPUTS.GENRANGE_DEC ); continue ; }
    
    if ( uniqueMatch(c_get,"GENRANGE_DEC:")  )
      { readfloat ( fp, 2, INPUTS.GENRANGE_DEC ); continue ; }

    if ( strstr(c_get,"SOLID_ANGLE") != NULL ) {
      parse_input_SOLID_ANGLE(fp, &iArg, c_get ); 
      continue ; 
    }

    if ( uniqueMatch(c_get,"GENRANGE_REDSHIFT:")  ) 
      { readdouble ( fp, 2, INPUTS.GENRANGE_REDSHIFT ); continue ; }
    
    if ( uniqueMatch(c_get,"GENSIGMA_REDSHIFT:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENSIGMA_REDSHIFT ); continue ; }
    
    if ( uniqueMatch(c_get,"GENBIAS_REDSHIFT:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENBIAS_REDSHIFT ); continue ; }    

    if ( uniqueMatch(c_get,"GENSIGMA_VPEC:")  )
      { readfloat ( fp, 1, &INPUTS.GENSIGMA_VPEC ); continue ; }
    
    if ( uniqueMatch(c_get,"VPEC_ERR:")  )
      { readfloat ( fp, 1, &INPUTS.VPEC_ERR ); continue ; }    

    if ( uniqueMatch(c_get,"VEL_CMBAPEX:")  ) 
      { readfloat ( fp, 1, &INPUTS.VEL_CMBAPEX ); continue ; }
 
    // - - - - -- 
    // wrongHost model params

    if ( uniqueMatch(c_get,"WRONGHOST_FILE:")  ) 
      { readchar( fp, INPUTS.WRONGHOST_FILE); continue ; }

    // ---- read/parse RATEPAR struct -----------
    read_input_RATEPAR(fp, "NOMINAL", c_get, &INPUTS.RATEPAR );
    read_input_RATEPAR(fp, "PEC1A",   c_get, &INPUTS.RATEPAR_PEC1A );
    // - - - - 

    if ( uniqueMatch(c_get,"GENRANGE_PEAKMAG:")  )
      { readdouble ( fp, 2, INPUTS.GENRANGE_PEAKMAG ); continue ; }

    if ( uniqueMatch(c_get,"GENRANGE_MJD:")  )
      { readdouble ( fp, 2, INPUTS.GENRANGE_MJD ); continue ; }

    if ( uniqueMatch(c_get,"GENRANGE_PEAKMJD:")  ) 
      { readdouble ( fp, 2, INPUTS.GENRANGE_PEAKMJD ); continue ; }

    if ( uniqueMatch(c_get,"GENSIGMA_SEARCH_PEAKMJD:")  )  // legacy key
      { readfloat ( fp, 1, &INPUTS.GENSIGMA_PEAKMJD ); continue ; }
    if ( uniqueMatch(c_get,"GENSIGMA_PEAKMJD:")  ) 
      { readfloat ( fp, 1, &INPUTS.GENSIGMA_PEAKMJD ); continue ; }
    if ( uniqueMatch(c_get,"OPT_SETPKMJD:")  ) 
      { readint ( fp, 1, &INPUTS.OPT_SETPKMJD ); continue ; }

    if ( uniqueMatch(c_get,"NEWMJD_DIF:")  ) 
      { readfloat ( fp, 1, &INPUTS.NEWMJD_DIF ); continue ; }

    if ( uniqueMatch(c_get,"GENRANGE_TREST:")  ) 
      { readfloat ( fp, 2, INPUTS.GENRANGE_TREST ); continue ; }

    if ( uniqueMatch(c_get,"GENRANGE_TOBS:")  ) 
      { readfloat ( fp, 2, INPUTS.GENRANGE_TOBS ); continue ; }

    if ( uniqueMatch(c_get,"TGRIDSTEP_MODEL_INTERP:")  ) 
      { readfloat ( fp, 1, &INPUTS.TGRIDSTEP_MODEL_INTERP ); continue ; }

    // contraint input file
    if ( uniqueMatch(c_get,"GENPAR_SELECT_FILE:")  ) 
      { readchar ( fp, INPUTS.GENPAR_SELECT_FILE ); continue ; }

    // ----- SIMSED parameters ---------

    if ( strstr(c_get,"SIMSED") != NULL ) 
      { read_input_SIMSED(fp,c_get); }

    N       = INPUTS.NPAR_SIMSED ;
    opt_tmp = INPUTS.GENFLAG_SIMSED[N-1];
    if ( N >= MXPAR_SIMSED ) {
      sprintf(c1err, "NPAR_SIMSED=%d exceeds array bound", N);
      sprintf(c2err, "Check SIMSED_PARAM keywords in input file." );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( N > 0 && (opt_tmp & OPTMASK_SIMSED_PARAM ) > 0 ) {
      parname = INPUTS.PARNAME_SIMSED[N-1];
      read_input_GENGAUSS(fp, c_get, parname,  &INPUTS.GENGAUSS_SIMSED[N-1] );
      INPUTS.GENGAUSS_SIMSED[N-1].FUNINDEX = N-1; 
    }

    // Jun 14, 2013
    // to avoid key confusion (e.g., GENRANGE_AV for mlcs or for SIMSED?),
    // skip lots of key checks if SIMSED model is set.
    if ( INPUTS.NPAR_SIMSED > 0 ) { goto READ_AVRV ; }


    // - - - - - - - - - - 
    // read risetime-shift info
    read_input_GENGAUSS(fp, c_get, "RISETIME_SHIFT",  
			&INPUTS.GENGAUSS_RISETIME_SHIFT );

    read_input_GENGAUSS(fp, c_get, "FALLTIME_SHIFT",  
			&INPUTS.GENGAUSS_FALLTIME_SHIFT );    

    // read shape par. Imput key depends on model,

    read_input_GENGAUSS(fp, c_get, "DELTA",    &INPUTS.GENGAUSS_DELTA    );
    read_input_GENGAUSS(fp, c_get, "DM15",     &INPUTS.GENGAUSS_DM15     );
    read_input_GENGAUSS(fp, c_get, "STRETCH",  &INPUTS.GENGAUSS_STRETCH  );

    if ( uniqueMatch(c_get,"STRETCH_TEMPLATE_FILE:" )   )
      { readchar(fp, INPUTS.STRETCH_TEMPLATE_FILE ); continue ; }

    // read SALT2 gen parameters
    read_input_GENGAUSS(fp, c_get, "SALT2c",     &INPUTS.GENGAUSS_SALT2c     );
    read_input_GENGAUSS(fp, c_get, "SALT2x1",    &INPUTS.GENGAUSS_SALT2x1    );
    read_input_GENGAUSS(fp, c_get, "SALT2ALPHA", &INPUTS.GENGAUSS_SALT2ALPHA );
    read_input_GENGAUSS(fp, c_get, "SALT2BETA",  &INPUTS.GENGAUSS_SALT2BETA  );

    if ( uniqueMatch(c_get,"BIASCOR_SALT2GAMMA_GRID:")   ) 
      { readdouble(fp, 2, INPUTS.BIASCOR_SALT2GAMMA_GRID );  continue ; }

    if ( uniqueMatch(c_get,"SALT2BETA_cPOLY:")   ) 
      { readdouble(fp, 3, INPUTS.SALT2BETA_cPOLY );  continue ; }

    if ( uniqueMatch(c_get,"SALT2mu_FILE:" )   )
      { readchar(fp, INPUTS.SALT2mu_FILE ); continue ; }

    // read legacy SALT2 alpha & beta key, but load modern variable
    if ( uniqueMatch(c_get,"GENALPHA_SALT2:")  ) {
      readfloat ( fp, 1, &INPUTS.GENALPHA_SALT2 );
      INPUTS.GENGAUSS_SALT2ALPHA.PEAK     = INPUTS.GENALPHA_SALT2 ; 
      continue ; 
    }
    if ( uniqueMatch(c_get,"GENBETA_SALT2:")  ) {
      readfloat ( fp, 1, &INPUTS.GENBETA_SALT2 );
      INPUTS.GENGAUSS_SALT2BETA.PEAK     = INPUTS.GENBETA_SALT2 ; 
      continue ; 
    }

    if ( uniqueMatch(c_get,"LEGACY_colorXTMW_SALT2:")    )
      { readint(fp, 1, &INPUTS.LEGACY_colorXTMW_SALT2); continue ; }

    // ---

  READ_AVRV:

    read_input_GENGAUSS(fp, c_get, "RV",  &INPUTS.GENGAUSS_RV );

    if ( uniqueMatch(c_get,"GENRANGE_AV:")  ) 
      { readdouble ( fp, 2, INPUTS.GENRANGE_AV ); continue ; }

    // allow old or new key for AV tau
    if ( uniqueMatch(c_get,"GENTAU_AV:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENEXPTAU_AV );  continue ; }
    if ( uniqueMatch(c_get,"GENEXPTAU_AV:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENEXPTAU_AV ); continue ; }

    if ( uniqueMatch(c_get,"GENSIG_AV:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENGAUSIG_AV ); continue ; }
    if ( uniqueMatch(c_get,"GENGAUSIG_AV:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENGAUSIG_AV ); continue ; }

    if ( uniqueMatch(c_get,"GENGAUPEAK_AV:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENGAUPEAK_AV ); continue ; }

    if ( uniqueMatch(c_get,"GENRATIO_AV0:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENRATIO_AV0 ); continue ; }

    if ( uniqueMatch(c_get,"GENAV_WV07:")  ) // legacy key name
      { readint ( fp, 1, &INPUTS.WV07_GENAV_FLAG ); continue ; }
    if ( uniqueMatch(c_get,"WV07_GENAV_FLAG:")  ) 
      { readint ( fp, 1, &INPUTS.WV07_GENAV_FLAG ); continue ; }
    if ( uniqueMatch(c_get,"WV07_REWGT_EXPAV:")  ) 
      { readdouble ( fp, 1, &INPUTS.WV07_REWGT_EXPAV ); continue ; }


  READ_MWEBV:

    if ( uniqueMatch(c_get,"RVMW:")   ) // allow legacy key (Sep 2013)
      { readdouble ( fp, 1, &INPUTS.RV_MWCOLORLAW ); continue ; }
    if ( uniqueMatch(c_get,"RV_MWCOLORLAW:")   )  // new key (Sep 2013)
      { readdouble ( fp, 1, &INPUTS.RV_MWCOLORLAW ); continue ; }

    if ( uniqueMatch(c_get,"OPT_MWCOLORLAW:")   ) 
      { readint ( fp, 1, &INPUTS.OPT_MWCOLORLAW ); continue ; }

    if ( uniqueMatch(c_get,"OPT_MWEBV:")  ) {
      readint ( fp, 1, &opt_tmp) ;
      INPUTS.OPT_MWEBV = abs(opt_tmp);
      if ( opt_tmp < 0 ) { INPUTS.APPLYFLAG_MWEBV=1; } // correct FLUXCAL
      continue ; 
    }

    // next are for systematic tests
    if ( uniqueMatch(c_get,"GENSIGMA_MWEBV_RATIO:")  )
      { readdouble ( fp, 1, &INPUTS.MWEBV_SIGRATIO ); continue ; }

    if ( uniqueMatch(c_get,"GENSIGMA_MWEBV:")  )
      { readdouble ( fp, 1, &INPUTS.MWEBV_SIG ); continue ; }

    if ( uniqueMatch(c_get,"GENSHIFT_MWEBV:")  ||
	 uniqueMatch(c_get,"FUDGESHIFT_MWEBV:")  )
      { readdouble ( fp, 1, &INPUTS.MWEBV_SHIFT ); continue ; }

    if ( uniqueMatch(c_get,"GENSCALE_MWEBV:")  ||
	 uniqueMatch(c_get,"FUDGESCALE_MWEBV:")  )
      { readdouble ( fp, 1, &INPUTS.MWEBV_SCALE ); continue ; }

    if ( uniqueMatch(c_get,"GENRANGE_MWEBV:")  )
      { readdouble ( fp, 2, INPUTS.GENRANGE_MWEBV ); continue ; }

    // ---- host -----

    if ( uniqueMatch(c_get,"EXTINC_HOSTGAL:")  ) 
      { readchar ( fp, INPUTS.GENSNXT );  continue ; }

    // ---

    if ( uniqueMatch(c_get,"GENMAG_OFF_AB:")  ||   // legacy key
	 uniqueMatch(c_get,"GENMAG_OFF_ZP:")  ) {
      if ( INPUTS.NFILTDEF_OBS == 0 ) {
	sprintf(c1err,"Filters NOT specified: cannot read AB offsets.");
	sprintf(c2err,"Please define filters BEFORE AB offsets.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      readfloat ( fp, INPUTS.NFILTDEF_OBS, INPUTS.TMPOFF_ZP);
      continue ;
    }


    if ( uniqueMatch(c_get,"GENMAG_OFF_GLOBAL:")  ) 
      { readfloat ( fp, 1, &INPUTS.GENMAG_OFF_GLOBAL ); continue ; }
    if ( uniqueMatch(c_get,"GENMAG_OFF_NON1A:")  ) 
      { readfloat ( fp, 1, &INPUTS.GENMAG_OFF_NON1A ); continue ; }

    if ( uniqueMatch(c_get,"GENMAG_OFF_MODEL:")  ) {
      if ( INPUTS.NFILTDEF_OBS == 0 ) {
	sprintf(c1err,"Filters NOT specified: cannot read MODEL offsets.");
	sprintf(c2err,"Please define filters BEFORE MODEL offsets.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      readfloat ( fp, INPUTS.NFILTDEF_OBS, INPUTS.TMPOFF_MODEL );
      continue ; 
    }


    if ( uniqueMatch(c_get,"GENMODEL_ERRSCALE:")  ) 
      {  readfloat ( fp, 1, &INPUTS.GENMODEL_ERRSCALE ); continue ; }
    if ( uniqueMatch(c_get,"GENMAG_SMEAR:")  ) {
      readchar(fp, ctmp); 
      split2floats(ctmp, comma, INPUTS.GENMAG_SMEAR );
      continue ;
    }
    //xxx mark delete   readfloat ( fp, 1, &INPUTS.GENMAG_SMEAR ); continue ; }


    if ( uniqueMatch(c_get,"GENMAG_SMEAR_USRFUN:")  ) { 
      INPUTS.NPAR_GENSMEAR_USRFUN     = 8 ; // fix hard-wired param
      readdouble ( fp, INPUTS.NPAR_GENSMEAR_USRFUN, 
		   INPUTS.GENMAG_SMEAR_USRFUN );
      sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"USRFUN");

      // turn off this option of first parameter is negative
      if ( INPUTS.GENMAG_SMEAR_USRFUN[0] <= -1.0E-7 ) {  
	sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"NONE"); 
	INPUTS.NPAR_GENSMEAR_USRFUN  = 0 ;
      }
      continue ; 
    }

    if ( uniqueMatch(c_get,"GENMAG_SMEAR_SCALE:")   ) 
      { readfloat(fp, 1, &INPUTS.GENMAG_SMEAR_SCALE ); continue ; }

    if ( uniqueMatch(c_get,"GENMAG_SMEAR_MODELNAME:")   ) {
      modelName = INPUTS.GENMAG_SMEAR_MODELNAME  ;
      readchar(fp, modelName );
      if ( strcmp(modelName,"G10FUDGE") == 0 ) 
	{ readdouble ( fp, 1, &INPUTS.GENMAG_SMEAR_USRFUN[0] );  }

      if ( strstr(modelName,":")!= NULL)  { parse_GENMAG_SMEAR_MODELNAME(); }

      continue ; 
    }

    // updated keys for strong & weak lensing (July 2019)
    if ( uniqueMatch(c_get,"STRONGLENS_FILE:")   ) 
      { readchar(fp, INPUTS.STRONGLENS_FILE); continue ; }

    if ( uniqueMatch(c_get,"WEAKLENS_PROBMAP_FILE:")   ) 
      { readchar(fp, INPUTS.WEAKLENS_PROBMAP_FILE); continue ; }
    if ( uniqueMatch(c_get,"WEAKLENS_DMUSCALE:")   ) 
      { readfloat(fp, 1, &INPUTS.WEAKLENS_DMUSCALE); continue ; }
    if ( uniqueMatch(c_get,"WEAKLENS_DSIGMADZ:")   ) 
      { readfloat(fp, 1, &INPUTS.WEAKLENS_DSIGMADZ); continue ; }

    // legacy keys for weak lensing (May 2017)
    if ( uniqueMatch(c_get,"LENSING_PROBMAP_FILE:")   ) 
      { readchar(fp, INPUTS.WEAKLENS_PROBMAP_FILE); continue ; }
    if ( uniqueMatch(c_get,"LENSING_DMUSCALE:")   ) 
      { readfloat(fp, 1, &INPUTS.WEAKLENS_DMUSCALE); continue ; }
    if ( uniqueMatch(c_get,"LENSING_DSIGMADZ:")   ) 
      { readfloat(fp, 1, &INPUTS.WEAKLENS_DSIGMADZ); continue ; }

    int NVAL;
    char key[60], parName[60];
    double tmpList[MXSMEARPAR_OVERRIDE];
    // allow multiple of these keys
    if ( keyMatch(c_get,"GENMAG_SMEARPAR_OVERRIDE:") ) { 
      readchar(fp, key);
      NVAL = nval_genSmear_override(key, parName); // return NVAL,parName
      readdouble(fp, NVAL, tmpList);               // read tmpList
      store_genSmear_override(parName,NVAL,tmpList); 
      continue ;
    }


    if ( uniqueMatch(c_get,"GENSMEAR_RANGauss_FIX:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENSMEAR_RANGauss_FIX ); continue ; }
    if ( uniqueMatch(c_get,"GENSMEAR_RANGAUSS_FIX:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENSMEAR_RANGauss_FIX ); continue ; }

    if ( uniqueMatch(c_get,"GENSMEAR_RANFlat_FIX:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENSMEAR_RANFlat_FIX ); continue ; }
    if ( uniqueMatch(c_get,"GENSMEAR_RANFLAT_FIX:")  ) 
      { readdouble ( fp, 1, &INPUTS.GENSMEAR_RANFlat_FIX ); continue ; }

    if ( uniqueMatch(c_get,"SIGMACLIP_MAGSMEAR:")  ) 
      { readfloat ( fp, 2, INPUTS.SIGMACLIP_MAGSMEAR ); continue ; }


    if ( uniqueMatch(c_get,"GENMAG_SMEAR_FILTER:")  ) {
      readchar(fp, ctmp );  // list of filters
      readchar(fp, ctmp2);  // float value or "INTERP" flag

      if ( strcmp(ctmp2,"INTERP") == 0 ) 
	{ ftmp = -1.0 ; } // any negative number -> interpolate smear
      else
	{ sscanf(ctmp2, "%f", &ftmp); }

      if ( ftmp != 0.0 ) {
	N = INPUTS.NFILT_SMEAR ;
	INPUTS.NFILT_SMEAR += 
	  PARSE_FILTLIST(ctmp, &INPUTS.IFILT_SMEAR[N+1] );
	for ( j=N+1; j <= INPUTS.NFILT_SMEAR; j++ ) {
	  ifilt = INPUTS.IFILT_SMEAR[j];
	  INPUTS.GENMAG_SMEAR_FILTER[ifilt] = ftmp ;
	}
      }
      continue ; 
    }

    if ( uniqueMatch(c_get,"GENMODEL_ERRSCALE_OPT:")  ) 
      {  readint(fp,1,&INPUTS.GENMODEL_ERRSCALE_OPT ); continue ; }

    if ( uniqueMatch(c_get,"GENMODEL_ERRSCALE_CORRELATION:")  ) 
      { readfloat(fp,1,&INPUTS.GENMODEL_ERRSCALE_CORRELATION ); continue;}

    // check intrinsic scatter

    if ( strncmp(c_get,"COVMAT_SCATTER",14) == 0 ) {
      sprintf(ckey,"%s", c_get);
      readchar(fp, cval);
      PARSE_COVMAT_SCATTER(ckey,cval) ;
      continue ; 
    }
    // -----

    if ( uniqueMatch(c_get,"GENRANGE_TYPE:")  ) 
      { readint ( fp, 2, INPUTS.GENRANGE_TYPE ); continue ; }
    if ( uniqueMatch(c_get,"GENRANGE_CID:")  ) 
      { readint ( fp, 2, INPUTS.GENRANGE_CID ); continue ; }


    // check old option to specify filter-index range;
    // Users should really switch to using 'GENFILTERS: abcd'

    if ( uniqueMatch(c_get,"GENRANGE_FILT:")  ) {
      sprintf(c1err,"GENRANGE_FILT: keyword no longer supported.");
      sprintf(c2err,"Must use GENFILTERS: ");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // Aug 13, 2007: new option to read filter-string so that
    // users need not worry about indices.

    if ( uniqueMatch(c_get,"GENFILTERS:")  ) {
      readchar ( fp, INPUTS.GENFILTERS );
      // must parse this NOW to get NFILT_OBS
      INPUTS.NFILTDEF_OBS =
	PARSE_FILTLIST( INPUTS.GENFILTERS, INPUTS.IFILTMAP_OBS  );
      continue ; 
    }


    if ( uniqueMatch(c_get,"GENRANGE_DMPEVENT:")  ) 
      { readint ( fp, 2, INPUTS.GENRANGE_DMPEVENT ); continue ; }

    if ( uniqueMatch(c_get,"GENRANGE_DMPTREST:")   ) 
      { readfloat ( fp, 2, INPUTS.GENRANGE_DMPTREST ); continue ; }

    if ( uniqueMatch(c_get,"SMEARFLAG_FLUX:")  ) 
      { readint ( fp, 1, &INPUTS.SMEARFLAG_FLUX ); continue ; }
    if ( uniqueMatch(c_get,"SMEARFLAG_ZEROPT:")  ) 
      { readint ( fp, 1, &INPUTS.SMEARFLAG_ZEROPT ); continue ; }

    if ( uniqueMatch(c_get,"MAGMONITOR_SNR:")  ) 
      { readint ( fp, 1, &INPUTS.MAGMONITOR_SNR ); continue ; }

    if ( uniqueMatch(c_get,"GENKCOR_STRETCH:")  ) 
      { readint ( fp, 1, &INPUTS.KCORFLAG_STRETCH ); }  

    if ( uniqueMatch(c_get,"KCORFLAG_COLOR:")  ) 
      { readint ( fp, 1, &INPUTS.KCORFLAG_COLOR ); }  

    if ( uniqueMatch(c_get,"EXPOSURE_TIME_MSKOPT:")  ) 
      { readint ( fp, 1, &INPUTS.EXPOSURE_TIME_MSKOPT ); } 

    if ( uniqueMatch(c_get,"KCOR_FILE:") ) 
      { readchar ( fp, INPUTS.KCOR_FILE ); continue ; }

    if ( uniqueMatch(c_get,"OMEGA_MATTER:")  ) 
      { readdouble ( fp, 1, &INPUTS.OMEGA_MATTER ) ; continue ; }

    if ( uniqueMatch(c_get,"OMEGA_LAMBDA:")  ) 
      { readdouble ( fp, 1, &INPUTS.OMEGA_LAMBDA ); continue ; }

    if ( uniqueMatch(c_get,"W0_LAMBDA:")  )
      { readdouble ( fp, 1, &INPUTS.W0_LAMBDA ); continue ; }

    if ( uniqueMatch(c_get,"H0:")  ) 
      { readdouble ( fp, 1, &INPUTS.H0 ); continue ; }

    // ------

    if ( uniqueMatch(c_get,"FUDGE_SNRMAX:")  ) { 
      readchar ( fp,  INPUTS.STRING_FUDGE_SNRMAX ); 
      INPUTS.OPT_FUDGE_SNRMAX = 1; 
      continue ; 
    }
    if ( uniqueMatch(c_get,"FUDGE2_SNRMAX:")  )  { 
      readchar( fp,  INPUTS.STRING_FUDGE_SNRMAX ); 
      INPUTS.OPT_FUDGE_SNRMAX = 2; 
      continue ; 
    }

    if ( uniqueMatch(c_get,"FUDGESCALE_PSF:")  ) 
      { readfloat ( fp, 1, &INPUTS.FUDGESCALE_PSF ); continue ; }

    if ( uniqueMatch(c_get,"FUDGESCALE_SKYNOISE:")  ) 
      { readfloat ( fp, 1, &INPUTS.FUDGESCALE_SKYNOISE ); continue ; }

    if ( uniqueMatch(c_get,"FUDGESCALE_READNOISE:")  ) 
      { readfloat ( fp, 1, &INPUTS.FUDGESCALE_READNOISE ); continue ; }

    if ( uniqueMatch(c_get,"FUDGEOPT_FLUXERR:")  )   
      { readint ( fp, 1, &INPUTS.FUDGEOPT_FLUXERR ); continue ; }


    parse_input_KEY_PLUS_FILTER(fp, &iArg, c_get, "EXPOSURE_TIME",
				&INPUTS.EXPOSURE_TIME, 
				INPUTS.EXPOSURE_TIME_FILTER);

    parse_input_KEY_PLUS_FILTER(fp, &iArg, c_get, "FUDGE_MAGERR",
				&INPUTS.FUDGE_MAGERR, 
				INPUTS.FUDGE_MAGERR_FILTER);

    parse_input_KEY_PLUS_FILTER(fp, &iArg, c_get, "FUDGESCALE_FLUXERR",
				&INPUTS.FUDGESCALE_FLUXERR, 
				INPUTS.FUDGESCALE_FLUXERR_FILTER );

    parse_input_KEY_PLUS_FILTER(fp, &iArg, c_get, "FUDGESCALE_FLUXERR2",
				&INPUTS.FUDGESCALE_FLUXERR2, 
				INPUTS.FUDGESCALE_FLUXERR2_FILTER );

    parse_input_KEY_PLUS_FILTER(fp, &iArg, c_get, "FUDGESHIFT_ZPT",
				&INPUTS.FUDGESHIFT_ZPT, 
				INPUTS.FUDGESHIFT_ZPT_FILTER);

    parse_input_KEY_PLUS_FILTER(fp, &iArg, c_get, "MJD_TEMPLATE",
				&INPUTS.MJD_TEMPLATE, 
				INPUTS.MJD_TEMPLATE_FILTER);

    // check for RANSYSTPAR params for systematic shifts (Jun 2017)
    if ( strstr(c_get,"RANSYSTPAR") != NULL ) 
      { parse_input_RANSYSTPAR(fp, &iArg, c_get ); }


    // ------

    // -- search eff
    if ( uniqueMatch(c_get,"GENPERFECT:")  ) 
      { readint ( fp, 1, &INPUTS.GENPERFECT ); continue ; }

    if ( uniqueMatch(c_get,"MAGSHIFT_SPECEFF:")  ) 
      { readdouble(fp,1,&INPUTS_SEARCHEFF.MAGSHIFT_SPECEFF ); continue;}

    if ( uniqueMatch(c_get,"SEARCHEFF_PIPELINE_LOGIC_FILE:")  ) 
      { readchar(fp,INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE ); continue;}

    // allow two different keys here
    if ( uniqueMatch(c_get,"SEARCHEFF_PIPELINE_FILE:")  ) 
      { readchar ( fp, INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE );continue ; }
    if ( uniqueMatch(c_get,"SEARCHEFF_PIPELINE_EFF_FILE:")  ) 
      { readchar ( fp, INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE ); continue ; }

    if ( uniqueMatch(c_get,"SEARCHEFF_SPEC_FILE:")  ) 
      { readchar ( fp, INPUTS_SEARCHEFF.USER_SPEC_FILE ); continue ; }
    if ( uniqueMatch(c_get,"SEARCHEFF_SPEC_SCALE:")  ) 
      { readdouble(fp,1,&INPUTS_SEARCHEFF.USER_SPECEFF_SCALE ); continue ; }

    if ( uniqueMatch(c_get,"SEARCHEFF_zHOST_FILE:")  ) 
      { readchar ( fp, INPUTS_SEARCHEFF.USER_zHOST_FILE ); continue ; }
    
    if ( uniqueMatch(c_get,"APPLY_SEARCHEFF_OPT:")  ) 
      { readint ( fp, 1, &INPUTS.APPLY_SEARCHEFF_OPT ); continue ; }

    if ( uniqueMatch(c_get,"MINOBS_SEARCH:")  ) 
      {  readint ( fp, 1, &INPUTS_SEARCHEFF.MINOBS ); continue ; }

    // --- LCLIB cut windows to select library events
    if ( uniqueMatch(c_get,"LCLIB_CUTWIN:")   ) {
      N = LCLIB_CUTS.NCUTWIN ;
      readchar(fp,      LCLIB_CUTS.PARNAME[N] );
      readdouble(fp, 2, LCLIB_CUTS.CUTWIN[N] );
      LCLIB_CUTS.NCUTWIN++ ;
      continue ; 
    }

    if ( uniqueMatch(c_get,"LCLIB_DEBUG_DAYSCALE:")   ) 
      { readdouble(fp, 1, &LCLIB_DEBUG.DAYSCALE ); continue ; }
    if ( uniqueMatch(c_get,"LCLIB_DEBUG_TOBS_OFFSET:")   ) 
      { readdouble(fp, 2, LCLIB_DEBUG.TOBS_OFFSET_RANGE ); continue ; }
    if ( uniqueMatch(c_get,"LCLIB_DEBUG_ZERO_TEMPLATE_FLUX:")   ) 
      { readint(fp, 1, &LCLIB_DEBUG.ZERO_TEMPLATE_FLUX ); continue ; }
    if ( uniqueMatch(c_get,"LCLIB_DEBUG_FORCE_NREPEAT:")   ) 
      { readint(fp, 1, &LCLIB_DEBUG.FORCE_NREPEAT ); continue ; }

    // -- cut-windows
    if ( uniqueMatch(c_get,"APPLY_CUTWIN_OPT:")  ) 
      { readint ( fp, 1, &INPUTS.APPLY_CUTWIN_OPT ); continue ; }

    if ( uniqueMatch(c_get,"EPCUTWIN_LAMREST:")  ) 
      { readfloat ( fp, 2, INPUTS.EPCUTWIN_LAMREST ); continue ; }

    if ( uniqueMatch(c_get,"EPCUTWIN_SNRMIN:")   ) 
      { readfloat ( fp, 2, INPUTS.EPCUTWIN_SNRMIN ); continue ; }

    if ( uniqueMatch(c_get,"CUTWIN_REDSHIFT:")  )  // legacy 
      { readfloat ( fp, 2, INPUTS.CUTWIN_REDSHIFT_TRUE ); continue ; }
    if ( uniqueMatch(c_get,"CUTWIN_REDSHIFT_TRUE:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_REDSHIFT_TRUE ); continue ; }
    if ( uniqueMatch(c_get,"CUTWIN_REDSHIFT_FINAL:")  )  // measured
      { readfloat ( fp, 2, INPUTS.CUTWIN_REDSHIFT_FINAL ); continue ; }

    if ( uniqueMatch(c_get,"CUTWIN_HOST_PHOTOZ:")  )  // host photo-z
      { readfloat ( fp, 2, INPUTS.CUTWIN_HOST_ZPHOT ); continue ; }
    if ( uniqueMatch(c_get,"CUTWIN_HOST_ZPHOT:")  )  // idem
      { readfloat ( fp, 2, INPUTS.CUTWIN_HOST_ZPHOT ); continue ; }

    if ( uniqueMatch(c_get,"CUTWIN_TRESTMAX:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_TRESTMAX ); continue ; }
    if ( uniqueMatch(c_get,"CUTWIN_TRESTMIN:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_TRESTMIN ); continue ; }

    if ( uniqueMatch(c_get,"CUTWIN_TGAPMAX:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_TGAPMAX ); continue ; }
    if ( uniqueMatch(c_get,"CUTWIN_T0GAPMAX:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_T0GAPMAX ); continue ; }

    if ( keyMatch(c_get,"CUTWIN_SNRMAX:")  ) {
      INPUTS.NCUTWIN_SNRMAX++ ; 
      N = INPUTS.NCUTWIN_SNRMAX ;
      readfloat ( fp, 1, &INPUTS.CUTWIN_SNRMAX[N][0] );
      readchar  ( fp, INPUTS.CUTWIN_SNRMAX_FILTERS[N] );
      readchar  ( fp, INPUTS.CUTWIN_SNRMAX_LOGIC[N] );
      readfloat ( fp, 2, &INPUTS.CUTWIN_SNRMAX_TREST[N][0] );
      continue ;
    }
    if ( uniqueMatch(c_get,"CUTWIN_NEPOCH:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_NEPOCH ); continue ; }

    if ( uniqueMatch(c_get,"CUTWIN_NOBSDIF:")  ) 
      { readint ( fp, 2, INPUTS.CUTWIN_NOBSDIF ); continue ; }
    if ( uniqueMatch(c_get,"CUTWIN_MJDDIF:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_MJDDIF ); continue ; }


    if ( uniqueMatch(c_get,"CUTWIN_NOBS_SATURATE:")  ) { // minobs maxobs filters
      N = INPUTS.NCUTWIN_SATURATE;  
      readint  ( fp, 2, INPUTS.CUTWIN_SATURATE_NOBS[N] ); 
      readchar ( fp,    INPUTS.CUTWIN_SATURATE_FILTERS[N] ); 
      INPUTS.NCUTWIN_SATURATE++ ;
      continue ; 
    }
    if ( uniqueMatch(c_get,"CUTWIN_NOBS_NOSATURATE:")  ) { 
      N = INPUTS.NCUTWIN_NOSATURATE;  
      readint  ( fp, 2, INPUTS.CUTWIN_NOSATURATE_NOBS[N] ); 
      readchar ( fp,    INPUTS.CUTWIN_NOSATURATE_FILTERS[N] ); 
      INPUTS.NCUTWIN_NOSATURATE++ ;
      continue ; 
    }

    if ( uniqueMatch(c_get,"CUTWIN_MWEBV:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_MWEBV ); continue ; }

    if ( uniqueMatch(c_get,"CUTWIN_PEAKMAG:")  ) 
      { readfloat ( fp, 2, INPUTS.CUTWIN_PEAKMAG ); continue ; }
    if ( uniqueMatch(c_get,"CUTWIN_PEAKMAG_ALL:")  )  // require ALL filters
      { readfloat ( fp, 2, INPUTS.CUTWIN_PEAKMAG_ALL ); continue ; }

    if ( uniqueMatch(c_get,"CUTWIN_PEAKMAG_BYFIELD:")  ) {
      INPUTS.NCUTWIN_PEAKMAG_BYFIELD++ ;  N= INPUTS.NCUTWIN_PEAKMAG_BYFIELD;
      readfloat ( fp, 2, INPUTS.CUTWIN_PEAKMAG_BYFIELD[N] );
      readchar  ( fp, INPUTS.CUTWIN_BYFIELDLIST[N] );
      continue ; 
    }

    // CUTWIN_EPOCHS_SNRMIN: 5 20 iz # SNR>5 for < 20 days, i or z
    if ( uniqueMatch(c_get,"CUTWIN_EPOCHS_SNRMIN:")  )  { 
      readfloat ( fp, 1, &INPUTS.CUTWIN_EPOCHS_SNRMIN ); 
      readfloat ( fp, 1, &INPUTS.CUTWIN_EPOCHS_TRANGE[1] ); 
      readchar  ( fp, INPUTS.CUTWIN_EPOCHS_FILTERS ); 
      continue ; 
    }

    if ( uniqueMatch(c_get,"EFFERR_STOPGEN:")  ) 
      { readfloat ( fp, 1, &INPUTS.EFFERR_STOPGEN ); continue ; }

    // - - - - - 
    if ( uniqueMatch(c_get,"SPECTROGRAPH_OPTMASK:")  ) 
      { readint ( fp, 1, &INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ); continue ; }

    if ( uniqueMatch(c_get,"SPECTROGRAPH_SCALE_TEXPOSE:")  ) { 
      readdouble ( fp, 1, &INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE ); 
      continue ; 
    }

    // - - - - - 
    if ( keyMatch(c_get,"WARP_SPECTRUM:")  ) 
      { readchar(fp,warp_spectrum_string);  GENSPEC.USE_WARP=1; continue ; }

    if ( keyMatch(c_get,"TAKE_SPECTRUM:")  ) 
      {  parse_input_TAKE_SPECTRUM(fp,warp_spectrum_string); continue ; }

    if ( keyMatch(c_get,"TAKE_SPECTRUM_HOSTFRAC:")  ) 
      {	readfloat(fp, 1, &INPUTS.TAKE_SPECTRUM_HOSTFRAC ); continue ; }

    if ( keyMatch(c_get,"TAKE_SPECTRUM_DUMPCID:")  ) 
      {	readint(fp, 1, &INPUTS.TAKE_SPECTRUM_DUMPCID ); continue ; }

    // check NGRID values for GRID option only

#ifdef SNGRIDGEN
    if ( strcmp(INPUTS.GENSOURCE,"GRID") == 0 )  {

      if ( uniqueMatch(c_get,"GRID_FORMAT:" )   ) 
	{ readchar(fp, GRIDGEN_INPUTS.FORMAT ); continue ; }

      if ( uniqueMatch(c_get,"NGRID_LOGZ:" )   ) 
	{ readint(fp,1,&GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_LOGZ]); continue ; }

      if ( uniqueMatch(c_get,"NGRID_LUMIPAR:" )   )  { // legacy key	
	readint(fp,1,&GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_SHAPEPAR]); 
	continue;
      }
      if ( uniqueMatch(c_get,"NGRID_SHAPEPAR:" )   ) 	{ 
	readint ( fp, 1, &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_SHAPEPAR] ); 
	continue ; 
      }

      if ( uniqueMatch(c_get,"NGRID_COLORPAR:" )   )  { 
	readint ( fp, 1, &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_COLORPAR] ); 
	continue ; 
      }
      if ( uniqueMatch(c_get,"NGRID_COLORLAW:" )   ) 	{ 
	readint ( fp, 1, &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_COLORLAW] ); 
	continue ; 
      }
      if ( uniqueMatch(c_get,"NGRID_TREST:" )   )  { 
	readint ( fp, 1, &GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_TREST] ); 
	continue ; 
      }
    }
#endif

    // abort on legacy keys ( Jun 2016)
    if ( uniqueMatch(c_get,"NGRID_LUMIPAR:" )   )   // legacy key
      { legacyKey_abort(fnam, "NGRID_LUMIPAR:", "NGRID_SHAPEPAR:"); }

    if ( uniqueMatch(c_get,"GEN_SNDATA_SIM:" )   )   // legacy key
      { legacyKey_abort(fnam, "GEN_SNDATA_SIM:", "FORMAT_MASK:"); }

    if ( uniqueMatch(c_get,"EXTINC_MILKYWAY:" )   )   // legacy key
      { legacyKey_abort(fnam, "EXTINC_MILKYWAY:", "OPT_MWEBV:"); }

    if ( uniqueMatch(c_get,"RVMW:" )   )   // legacy key
      { legacyKey_abort(fnam, "RVMW:", "RV_MWCOLORLAW:"); }

    if ( uniqueMatch(c_get,"HUMAN_SEARCHEFF_OPT:" )   )   // legacy key
      { legacyKey_abort(fnam, "HUMAN_SEARCHEFF_OPT:", "SEARCHEFF_SPEC_FILE:");}

    if ( uniqueMatch(c_get,"SPECTYPE:" )   )   // legacy key
      { legacyKey_abort(fnam, "SPECTYPE:", ""); }

    if ( uniqueMatch(c_get,"SMEARFLAG_HOSTGAL:" )   )   // legacy key
      { legacyKey_abort(fnam, "SMEARFLAG_HOSTGAL:", ""); }

  }  // end fscanf loop


  fclose(fp);
  printf("\n");
  return(SUCCESS) ;

} // end of read_input


// *******************************************
void  read_input_RATEPAR(FILE *fp, char *WHAT, char *KEYNAME, 
			 RATEPAR_DEF *RATEPAR ) {

  // Created Aug 2016
  // read input file for DNDZ keys.
  // *WHAT = 'NOMINAL' or 'PEC1A'
  // 
  // Jan 27 2017: check for turning off rate model (for PEC1A)
  // Jul 29 2017: for FLAT, set RATEPAR->NMODEL_ZRANGE = 1 
  // Mar 30 2019: parse CC_S15*[scale] to allow rate-scale for this model.

  int  ISNOMINAL = strcmp(WHAT,"NOMINAL") == 0 ;
  int  ISPEC1A   = strcmp(WHAT,"PEC1A"  ) == 0 ;

  double l=0.0, b=0.0, bmax, R=0.0, TMPVAL ;
  //  char PRIMARY_KEYLIST[10][20], c_get[80] ;
  int  FOUND_PRIMARY_KEY, NLOCAL, CONTINUE ;
  char fnam[] = "read_input_RATEPAR" ;

  // ----------- BEGIN ----------

  CONTINUE = 0 ;
  if ( strstr(KEYNAME,"DNDZ") != NULL ) { CONTINUE = 1 ; }
  if ( strstr(KEYNAME,"DNDB") != NULL ) { CONTINUE = 1 ; }
  if ( CONTINUE == 0 ) { return ; }

  // check a few misc keys
  if ( ISNOMINAL ) {
    if ( strcmp(KEYNAME,"DNDZ_ZEXP_REWGT:")==0 ) 
      {  readdouble ( fp, 1, &RATEPAR->DNDZ_ZEXP_REWGT ); return; }

    if ( strcmp(KEYNAME,"DNDZ_ZPOLY_REWGT:")==0 ) 
      {  readdouble ( fp, 4, RATEPAR->DNDZ_ZPOLY_REWGT ); return ; }
    
    if ( strcmp(KEYNAME,"DNDZ_SCALE:")==0 ) 
      { readdouble ( fp, 2, RATEPAR->DNDZ_SCALE ); return ; }

    if ( strcmp(KEYNAME,"DNDZ_ALLSCALE:")==0 ) 
      { readdouble ( fp, 1, &RATEPAR->DNDZ_ALLSCALE ); return ; }

    // check legacy key ...
    if ( strcmp(KEYNAME,"DNDZ_SCALE_NON1A:")==0 ||
	 strcmp(KEYNAME,"DNDZ_SCALE_NONIA:")==0 )  { 
      readdouble ( fp, 1, &RATEPAR->DNDZ_SCALE[1] ); return ; 
    }
  }

  FOUND_PRIMARY_KEY = valid_DNDZ_KEY(WHAT, 1, KEYNAME ) ;

  // --------------------

  if ( FOUND_PRIMARY_KEY ) {

    readchar ( fp, RATEPAR->NAME ); 

    if ( IGNOREFILE(RATEPAR->NAME) ) { return ; }

    // make sure that PEC1A rateModel is defined only for NON1A mode
    if ( ISPEC1A &&  INPUTS.NON1A_MODELFLAG < 0 ) {
      sprintf(c1err,"%s key allowed only for NON1A model.", KEYNAME);
      sprintf(c2err,"Check sim-input file." );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( strcmp(RATEPAR->NAME,"ABMODEL") == 0 ) {
	RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_AB ;
	RATEPAR->NMODEL_ZRANGE = 1 ;  
	readdouble ( fp, 2, RATEPAR->MODEL_PARLIST[1] ); 
    }
    else if ( strcmp(RATEPAR->NAME,"POWERLAW") == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_POWERLAW ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      readdouble ( fp, 2, RATEPAR->MODEL_PARLIST[1] ); 
    }     
    else if ( strcmp(RATEPAR->NAME,"POWERLAW2") == 0 ) {
      RATEPAR->NMODEL_ZRANGE++ ;
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_POWERLAW2 ;
      NLOCAL = RATEPAR->NMODEL_ZRANGE ;
      readdouble ( fp, 2, RATEPAR->MODEL_PARLIST[NLOCAL] ); 
      readdouble ( fp, 2, RATEPAR->MODEL_ZRANGE[NLOCAL] ); 
    }
    /* xxxx mark delete Mar 30 2019 xxxxxxxxxx 
    else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_CCS15) == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_CCS15 ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
    }
    xxxxxxxxx */
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
      readdouble ( fp, 1, &RATEPAR->MODEL_PARLIST[1][0] ); // rate at z=0
    }
    else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_MD14) == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_MD14 ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      readdouble ( fp, 1, &RATEPAR->MODEL_PARLIST[1][0] ); // rate at z=0
    }
    else if ( strcmp(RATEPAR->NAME,"ZPOLY") == 0 ) {
      RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_ZPOLY ;
      RATEPAR->NMODEL_ZRANGE = 1 ;
      readdouble ( fp, 4, RATEPAR->MODEL_PARLIST[1] ); 
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
      readdouble ( fp, MXPOLY_GALRATE+1, RATEPAR->MODEL_PARLIST[1] ); 

      // get max rate (vs. b) for weighting

      for(b=0.0; b < 90.0; b+=1.0 ) {
	R = GALrate_model(l,b, RATEPAR);
	if ( R > RATEPAR->RATEMAX ) 
	  { RATEPAR->RATEMAX=R; bmax=b; }
      }
      //  printf(" xxx max GALrate=%f at b=%f \n", RATEPAR->RATEMAX,bmax);
    }
    else if ( strcmp(RATEPAR->NAME,"BPOLY") == 0 ) {
      RATEPAR->INDEX_MODEL   = INDEX_RATEMODEL_BPOLY ;
      RATEPAR->NMODEL_ZRANGE = 0 ;
      readdouble ( fp, MXPOLY_GALRATE+1, RATEPAR->MODEL_PARLIST[1] ); 

      // get max rate (vs. b) for weighting
      for(b=0.0; b < 90.0; b+=1.0 ) {
	R = GALrate_model(l,b, RATEPAR);
	if ( R > RATEPAR->RATEMAX ) 
	  { RATEPAR->RATEMAX=R; bmax=b; }
      }
      //  printf(" xxx max GALrate=%f at b=%f \n", RATEPAR->RATEMAX,bmax);
    }
    else {
      sprintf(c1err,"'%s %s' is invalid", KEYNAME, RATEPAR->NAME );
      sprintf(c2err,"Check sim-input file '%s'", 
	      INPUTS.INPUT_FILE_LIST[0] ); 
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
    
  } // DNDZ


  return ;

} // end read_input_RATEPAR

// ====================================
int  valid_DNDZ_KEY(char *WHAT, int fromFile, char *KEYNAME ) {

  // Created Aug 2016
  // Return True if this is a valid DNDZ key.
  //
  // Inputs:
  //   *WHAT     = 'NOMINAL' or 'PEC1A'
  //   fromFile  = 1 ==> read from file with colon after key
  //   fromFile  = 0 ==> command-line override, so no colon
  //   KEYNAME   = name of key to test
    
  int  ISNOMINAL = strcmp(WHAT,"NOMINAL") == 0 ;
  int  ISPEC1A   = strcmp(WHAT,"PEC1A"  ) == 0 ;
  char PRIMARY_KEYLIST[10][20], KEYTEST[20] ;
  int  ikey, NKEY, FOUND ;
  char fnam[] = "valid_DNDZ_KEY" ;

  // ----------- BEGIN -----------
  FOUND = 0 ;

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
    if ( fromFile ) { strcat(KEYTEST,":"); }
    if ( strcmp(KEYNAME,KEYTEST) == 0 )  { FOUND = 1 ; }
  }

  return(FOUND);

} // end valid_DNDZ_KEY 


// ************************************************************
void  sscanf_RATEPAR(int *i, char *WHAT, RATEPAR_DEF *RATEPAR) {

  // Aug 2016
  // read/fill RATEPAT from command-line override
  // *WHAT = 'NOMINAL' or 'PEC1A'
  //
  // June 18 2018:
  //  + fix bug and set RATEPAR->INDEX_MODEL for POWERLAW, POWERLAW2.
  //

  int iLoc = *i;
  int j, FOUND, NLOCAL=0;
  double TMPVAL;
  char *KEY ;
  char fnam[] = "sscanf_RATEPAR" ;

  // -------------- BEGIN -----------------

  KEY   = ARGV_LIST[iLoc] ;

  if ( strstr(KEY,"DNDZ") == NULL ) { return ; }

  FOUND = valid_DNDZ_KEY(WHAT, 0, KEY ) ;

  if ( FOUND == 0 ) { return ; }

  iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%s", RATEPAR->NAME ); 


  NLOCAL++ ;
  if ( strcmp(RATEPAR->NAME,"ABMODEL") == 0 ) {
    RATEPAR->NMODEL_ZRANGE = 1 ;
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_PARLIST[1][0] ); 
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_PARLIST[1][1] ); 
  }


  else if ( strcmp(RATEPAR->NAME,"POWERLAW") == 0 ) {
    RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_POWERLAW ;
    RATEPAR->NMODEL_ZRANGE = 1 ;
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_PARLIST[1][0] ); 
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_PARLIST[1][1] ); 
  }
  
  else if ( strcmp(RATEPAR->NAME,"POWERLAW2") == 0 ) {
    RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_POWERLAW2 ;
    RATEPAR->NMODEL_ZRANGE++ ;   NLOCAL = RATEPAR->NMODEL_ZRANGE ;
    
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_PARLIST[NLOCAL][0] ); 
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_PARLIST[NLOCAL][1] ); 
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_ZRANGE[NLOCAL][0] ); 
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_ZRANGE[NLOCAL][1] ); 
  }

  else if ( strcmp(RATEPAR->NAME,"HUBBLE") == 0 ) {
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
  /* xxxxxxxxxxxx mark delete Mar 30 2019 xxxxxxxx
  else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_CCS15) == 0 ) {
    RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_CCS15 ;
    RATEPAR->NMODEL_ZRANGE = 1 ;
  }
  xxxxxxxxxxxxxxxxx */
  else if ( strstr(RATEPAR->NAME,RATEMODELNAME_CCS15) != NULL ) {
    parse_multiplier(RATEPAR->NAME,RATEMODELNAME_CCS15, &TMPVAL);
    sprintf(RATEPAR->NAME,"%s", RATEMODELNAME_CCS15); // strip off scale
    RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_CCS15 ;
    RATEPAR->NMODEL_ZRANGE = 1 ;
    RATEPAR->MODEL_PARLIST[1][0] = TMPVAL ;
  }
  else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_PISN) == 0 ) {
    RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_PISN ;
    RATEPAR->NMODEL_ZRANGE = 1 ;
  }
  else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_TDE) == 0 ) {
    RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_TDE ;
    RATEPAR->NMODEL_ZRANGE = 1 ;
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_PARLIST[1][0] ); 
  }
  else if ( strcmp(RATEPAR->NAME,RATEMODELNAME_MD14) == 0 ) {
    RATEPAR->INDEX_MODEL = INDEX_RATEMODEL_MD14 ;
    RATEPAR->NMODEL_ZRANGE = 1 ;
    iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		    &RATEPAR->MODEL_PARLIST[1][0] ); 
  }
  else if ( strcmp(RATEPAR->NAME,"ZPOLY") == 0 ) {
    RATEPAR->NMODEL_ZRANGE = 1 ;
    for(j=0; j<4; j++ ) {
      iLoc++ ; sscanf(ARGV_LIST[iLoc] , "%le", 
		      &RATEPAR->MODEL_PARLIST[1][j] ); 
    }
  }
  else if ( IGNOREFILE(RATEPAR->NAME) ) {
    // if model = NULL, NONE, BLANK, etc ..., turn model off.
    RATEPAR->NMODEL_ZRANGE = 0 ;
  }
  else {
    sprintf(c1err,"'%s %s' is invalid", KEY, RATEPAR->NAME);
    sprintf(c2err,"Check sim-input file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  
  *i = iLoc ;

  return ;

} // end sscanf_RATEPAR


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
  char fnam[] = "parse_input_FIXMAG" ;

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

// ==============================================================
void parse_input_RANSYSTPAR(FILE *fp, int *iArg, char *KEYNAME ) {

  // Created Jun 2017
  // parse RANSYSTPAR_XXX params

  int  i = *iArg ;
  int  ifilt, NFILTDEF = INPUTS.NFILTDEF_OBS ;  
  char fnam[] = "parse_input_RANSYSTPAR" ;
  
  // ---------- BEGIN ------------


  if ( NFILTDEF == 0 ) {
    sprintf(c1err, "Cannot read RANSYSTPAR because NFILTDEF=0");
    sprintf(c2err, "Make sure %s appears after GENFILTERS: key.",KEYNAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  INPUTS.RANSYSTPAR.USE = 1;

  if ( i < 0 ) {
    // read from file
    if ( strcmp(KEYNAME,"RANSYSTPAR_FUDGESCALE_FLUXERR:")==0 ) {
      readfloat ( fp, 1, &INPUTS.RANSYSTPAR.FUDGESCALE_FLUXERR);
    }
    else if ( strcmp(KEYNAME,"RANSYSTPAR_FUDGESCALE_FLUXERR2:")==0 ) {
      readfloat ( fp, 1, &INPUTS.RANSYSTPAR.FUDGESCALE_FLUXERR2);
    }
    else if ( strcmp(KEYNAME,"RANSYSTPAR_FUDGESCALE_MWEBV:")==0 ) {
      readfloat ( fp, 1, &INPUTS.RANSYSTPAR.FUDGESCALE_MWEBV );
    }
    else if ( strcmp(KEYNAME,"RANSYSTPAR_FUDGESHIFT_MWRV:")==0 ) {
      readfloat ( fp, 1, &INPUTS.RANSYSTPAR.FUDGESHIFT_MWRV );
    }

    // filter-dependent params
    else if ( strcmp(KEYNAME,"RANSYSTPAR_GENMAG_OFF_ZP:")==0 ) {
      readfloat ( fp, NFILTDEF, INPUTS.RANSYSTPAR.GENMAG_OFF_ZP );
    }
    else if ( strcmp(KEYNAME,"RANSYSTPAR_FILTER_LAMSHIFT:")==0 ) {
      readfloat ( fp, NFILTDEF, INPUTS.RANSYSTPAR.FILTER_LAMSHIFT );
    }

  }
  else {
    // read command-line arg
    if ( strcmp(ARGV_LIST[i],"RANSYSTPAR_FUDGESCALE_FLUXERR") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", 
		   &INPUTS.RANSYSTPAR.FUDGESCALE_FLUXERR ); 
    }
    else if ( strcmp(ARGV_LIST[i],"RANSYSTPAR_FUDGESCALE_FLUXERR2") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", 
		   &INPUTS.RANSYSTPAR.FUDGESCALE_FLUXERR2 ); 
    }
    else if ( strcmp(ARGV_LIST[i],"RANSYSTPAR_FUDGESCALE_MWEBV") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", 
		   &INPUTS.RANSYSTPAR.FUDGESCALE_MWEBV ); 
    }
    else if ( strcmp(ARGV_LIST[i],"RANSYSTPAR_FUDGESHIFT_MWRV") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", 
		   &INPUTS.RANSYSTPAR.FUDGESHIFT_MWRV ); 
    }

    else if ( strcmp(ARGV_LIST[i],"RANSYSTPAR_GENMAG_OFF_ZP") == 0 ) {
      for ( ifilt = 0; ifilt < NFILTDEF; ifilt++ ) {
	i++ ; sscanf(ARGV_LIST[i] , "%f", 
		     &INPUTS.RANSYSTPAR.GENMAG_OFF_ZP[ifilt] ); 
      }
    }
    else if ( strcmp(ARGV_LIST[i],"RANSYSTPAR_FILTER_LAMSHIFT") == 0 ) {
      for ( ifilt = 0; ifilt < NFILTDEF; ifilt++ ) {
	i++ ; sscanf(ARGV_LIST[i] , "%f", 
		     &INPUTS.RANSYSTPAR.FILTER_LAMSHIFT[ifilt] ); 
      }
    }

    *iArg = i ;
  }

  return ;

} // end parse_input_RANSYSTPAR

// ==============================================================
void parse_input_SOLID_ANGLE(FILE *fp, int *iArg, char *KEYNAME ) {

  // Created Sep 27 2017
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

  char *FIELDLIST = INPUTS.SIMLIB_FIELDLIST ;
  char KEYLOCAL[100], FIELDLIST_ADD[100];
  int  ADDFIELD = 0 ;
  int  i = *iArg;
  char fnam[] = "parse_input_SOLID_ANGLE" ;


  // ------------ BEGIN ---------------

  sprintf(KEYLOCAL, "%s", KEYNAME);

  //  printf(" xxx %s: KEYNAME = '%s' \n", fnam, KEYNAME);

  // extract contents of optional ()
  extractStringOpt(KEYLOCAL,FIELDLIST_ADD);


  if ( i < 0 ) {
    // read from file
    if ( strcmp(KEYLOCAL,"SOLID_ANGLE:") != 0 ) { return ; }
    readfloat(fp, 1, &INPUTS.SOLID_ANGLE);
  }
  else {
    // read from command line
    if ( strcmp(KEYLOCAL,"SOLID_ANGLE") != 0 ) { return ; }
    i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.SOLID_ANGLE ); 
  }

  *iArg = i;

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

  return ;

} // end parse_input_SOLID_ANGLE

// ==============================================================
void parse_input_SIMGEN_DUMP(FILE *fp, int *iArg, char *KEYNAME ) {

  // Created Aug 2017
  // parse SIMGEN_DUMP params, including prescale.
  // 
  // Nov 8 2017: fix bug and read SIMGEN_DUMP[ALL] from command-line arg.
  // Nov 27 2017: 
  //  + if user passes 'SIMGEN_DUMP 0' then set NVAR_SIMGEN_DUMP=-9.
  //    (note that NVAR_SIMGEN_DUMP=0 --> list variables and quit)
  //
  int  i = *iArg, itmp, LRD=0 ;
  char fnam[] = "parse_input_SIMGEN_DUMP" ;
  
  // ---------- BEGIN ------------

  if ( i < 0 ) {
    // read from file
    if ( strcmp(KEYNAME,"SIMGEN_DUMP:")==0 ) {
      // dump SN that are written to data files
      { readint ( fp, 1, &INPUTS.NVAR_SIMGEN_DUMP ); }
      for ( itmp=0; itmp < INPUTS.NVAR_SIMGEN_DUMP; itmp++ ) 
	{ readchar(fp, INPUTS.VARNAME_SIMGEN_DUMP[itmp] );  }
    }

    if ( strcmp(KEYNAME,"SIMGEN_DUMPALL:")==0 ) {
      { readint ( fp, 1, &INPUTS.NVAR_SIMGEN_DUMP ); }
      INPUTS.IFLAG_SIMGEN_DUMPALL = 1; // dump EVERY generated SN
      for ( itmp=0; itmp < INPUTS.NVAR_SIMGEN_DUMP; itmp++ )
	{ readchar(fp, INPUTS.VARNAME_SIMGEN_DUMP[itmp] ); 	}
    }

    if ( strcmp(KEYNAME,"PRESCALE_SIMGEN_DUMP:")==0 ) 
      { readint ( fp, 1, &INPUTS.PRESCALE_SIMGEN_DUMP ); }
    if ( strcmp(KEYNAME,"SIMGEN_DUMP_PRESCALE:")==0 )   // allow alternate key
      { readint ( fp, 1, &INPUTS.PRESCALE_SIMGEN_DUMP ); }
  }
  else {
    // read command-line arg
    if ( strcmp(ARGV_LIST[i],"PRESCALE_SIMGEN_DUMP") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.PRESCALE_SIMGEN_DUMP ); 
    }
    if ( strcmp(ARGV_LIST[i],"SIMGEN_DUMP_PRESCALE") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.PRESCALE_SIMGEN_DUMP ); 
    }


    if ( strcmp(ARGV_LIST[i],"SIMGEN_DUMP") == 0 )
      { LRD=1; }
    if ( strcmp(ARGV_LIST[i],"SIMGEN_DUMPALL") == 0 ) 
      { LRD=1; INPUTS.IFLAG_SIMGEN_DUMPALL=1 ; }

    if ( LRD ) {      
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.NVAR_SIMGEN_DUMP);

      if ( INPUTS.NVAR_SIMGEN_DUMP <= 0 ) 
	{ INPUTS.IFLAG_SIMGEN_DUMPALL=0; }

      for(itmp=0; itmp < INPUTS.NVAR_SIMGEN_DUMP; itmp++ ) {
	i++ ; sscanf(ARGV_LIST[i], "%s", INPUTS.VARNAME_SIMGEN_DUMP[itmp] );
      }
    }

    *iArg = i ;
  }

  return ;

} // end parse_input_SIMGEN_DUMP

// ==================================================
int parse_input_KEY_PLUS_FILTER(FILE *fp, int *iArg, 
				char *INPUT_STRING, char *KEYCHECK,
				float *VALUE_GLOBAL,
				float *VALUE_FILTERLIST) {

  // Created Jun 2017
  // INPUT_STRING is current string read from input file.
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
  // Note that outputs are float, not double.

  float ftmp;
  char cfilt[MXFILTINDX], KEY[80] ;
  int  i, NTMP, ifilt_obs, ifilt, ifilt_list[MXFILTINDX];
  char fnam[] = "parse_input_KEY_PLUS_FILTER" ;

  // ----------- BEGIN -----------

  i = *iArg;

  if ( i < 0 ) { 
    // read from file *fp
    sprintf(KEY,"%s:", KEYCHECK);
    if ( strcmp(INPUT_STRING,KEY)==0 ) { 
      readfloat ( fp, 1, VALUE_GLOBAL );  
      for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ )
	{ VALUE_FILTERLIST[ifilt] = *VALUE_GLOBAL; }
      return(1);
    }
    
    sprintf(KEY,"%s_FILTER:", KEYCHECK);
    if ( strcmp(INPUT_STRING,KEY)==0 ) { 
      readchar  ( fp, cfilt );
      readfloat ( fp, 1, &ftmp );
      NTMP = PARSE_FILTLIST(cfilt, ifilt_list );  // return ifilt_list
      
      /*
      printf(" xxx read %s cfilt=%s ftmp=%.2f NTMP=%d \n",
      KEY, cfilt, ftmp, NTMP ); */
      
      for ( ifilt=0; ifilt < NTMP; ifilt++ ) {
	ifilt_obs = ifilt_list[ifilt] ;
	VALUE_FILTERLIST[ifilt_obs] = ftmp ;
      }
      return(NTMP);
    }
  }
  else {
    // read from command-line args   
    sprintf(KEY,"%s", KEYCHECK);
    if ( strcmp( ARGV_LIST[i], KEY ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", VALUE_GLOBAL ); 
      //      printf(" xxx %s:  %s = %f \n", fnam, KEY, *VALUE_GLOBAL);
      *iArg = i ;
      for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ ) 
	{ VALUE_FILTERLIST[ifilt] = *VALUE_GLOBAL ; }
      return(1);
    }

    sprintf(KEY,"%s_FILTER", KEYCHECK);
    if ( strcmp( ARGV_LIST[i], KEY ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", cfilt );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &ftmp );
      *iArg = i ;
      NTMP = PARSE_FILTLIST(cfilt, ifilt_list );  // return ifilt_obs
      for ( ifilt=0; ifilt < NTMP; ifilt++ ) {
	ifilt_obs = ifilt_list[ifilt] ;
	VALUE_FILTERLIST[ifilt_obs] = ftmp;
      }      
      return(NTMP);
    }


  } // end i if-block

  return(0);

} // end parse_input_KEY_PLUS_FILTER

// ==============================================================
void parse_GENMAG_SMEAR_MODELNAME(void) {

  // Split GENMAG_SMEAR_MODELSTRING by colon;
  // store right side of colin in MODELARG.

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

  free(inString); free(ptrSplit[0]) ;  free(ptrSplit[1]) ;

  return ;

} // end parse_GENMAG_SMEAR_MODELNAME

// ==============================================================
void parse_input_GENMODEL(FILE *fp, int *iArg ) {

  // Created April 18 2019 
  //   [code moved from read_input, part of refactor ]
  // Parse GENMODEL key.

  int  i = *iArg ;
  int  jnam;
  char *GENMODEL = INPUTS.GENMODEL ;
  char ctmp[60], *NAME0, *NAME1;
  int  LDMP   = 0 ;
  char fnam[] = "parse_input_GENMODEL" ;
  
  // ---------- BEGIN ------------

  if ( i < 0 ) {
    // read from file 
    readchar ( fp, GENMODEL );      
    
    // check for path + model (Feb 22 2017)
    extract_MODELNAME(GENMODEL,                   // input path/model
		      INPUTS.MODELPATH, INPUTS.MODELNAME); // returned
    
    // check for NONIA model to help parsing below
    INPUTS.NON1A_MODELFLAG = get_NON1A_MODELFLAG(INPUTS.MODELNAME);
    if ( INPUTS.NON1A_MODELFLAG == MODEL_NON1AGRID ) 
      { readchar(fp, INPUTS.NON1AGRID_FILE ); }

    NAME0 = GENMODEL_NAME[MODEL_LCLIB][0];
    if ( strcmp(GENMODEL,NAME0)==0 ) {
      readchar(fp, INPUTS.LCLIB_FILE );  
      readchar(fp, INPUTS.LCLIB_TEMPLATE_EPOCHS );  
    }

    NAME0 = GENMODEL_NAME[MODEL_BYOSED][0];
    if ( strcmp(GENMODEL,NAME0)==0 ) 
      { readchar(fp, INPUTS.MODELPATH ); }

    // for FIXMAG model, read value of fixed mag
    for(jnam=0; jnam < MXNAME_PER_MODEL; jnam++ ) {
      NAME0  = GENMODEL_NAME[MODEL_FIXMAG][jnam];
      if ( strcmp(GENMODEL,NAME0)==0 ) {
	readchar(fp, ctmp );
	parse_input_FIXMAG(ctmp);
	if ( strstr(GENMODEL,"mag") != NULL ) 
	  { INPUTS.GENFRAME_FIXMAG = GENFRAME_OBS; }
	else if ( strstr(GENMODEL,"MAG") != NULL ) 
	  { INPUTS.GENFRAME_FIXMAG = GENFRAME_REST; }
      }
    }
    
    // allow mlcs in place of mlcs2k2
    if ( strcmp(GENMODEL,"mlcs")==0 ) { sprintf(GENMODEL,"mlcs2k2"); }

    
  }
  else {
    // read command-line arg

    i++ ; sscanf(ARGV_LIST[i] , "%s", GENMODEL ); 

    // check for path + model 
    extract_MODELNAME(GENMODEL,                   // input path/model
		      INPUTS.MODELPATH, INPUTS.MODELNAME); // returned

    // check for NONIA model to help parsing below
    INPUTS.NON1A_MODELFLAG = get_NON1A_MODELFLAG(INPUTS.GENMODEL);
    if ( INPUTS.NON1A_MODELFLAG == MODEL_NON1AGRID ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.NON1AGRID_FILE );  }

    NAME0 = GENMODEL_NAME[MODEL_LCLIB][0];
    if ( strcmp(GENMODEL,NAME0)==0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.LCLIB_FILE ); 
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.LCLIB_TEMPLATE_EPOCHS ); 
    }

    NAME0 = GENMODEL_NAME[MODEL_BYOSED][0];
    if ( strcmp(GENMODEL,NAME0)==0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%s", INPUTS.MODELPATH ); 
    }

    int jnam;
    for(jnam=0; jnam < MXNAME_PER_MODEL; jnam++ ) {
      NAME0  = GENMODEL_NAME[MODEL_FIXMAG][jnam];
      if ( strcmp(GENMODEL,NAME0) == 0 ) {
	i++ ; sscanf(ARGV_LIST[i], "%s", ctmp ); 
	parse_input_FIXMAG(ctmp);	
	if ( strstr(GENMODEL,"mag") != NULL ) 
	  { INPUTS.GENFRAME_FIXMAG = GENFRAME_OBS; }
	else if ( strstr(GENMODEL,"MAG") != NULL ) 
	  { INPUTS.GENFRAME_FIXMAG = GENFRAME_REST; }
      }
    }

    // allow mlcs in place of mlcs2k2
    if ( strcmp(GENMODEL,"mlcs")== 0) { sprintf(GENMODEL, "mlcs2k2"); }    

    *iArg = i ;
  }


  // - - - - - - - - - - - - - - - - 
  if ( LDMP ) {
    printf("\n xxx -------------------------------- \n");
    printf(" xxx %s DUMP: \n", fnam );
    printf(" xxx Input GENMODEL:  '%s' \n", GENMODEL);
    printf(" xxx INPUTS.MODELPATH = '%s' \n", INPUTS.MODELPATH );
    printf(" xxx INPUTS.MODELNAME = '%s' \n", INPUTS.MODELNAME );
    debugexit(fnam);
  }


  return ;

} // end parse_input_GENMODEL


// ==============================================================
void parse_input_GENMODEL_ARGLIST(FILE *fp, int *iArg ) {

  // Created March 2019
  // parse GENMODEL_ARGLIST 'XX YY ZZ'
  // which can have arbitrary number of elements between quotes.
  // Initial use is to pass command-line overrides to BYOSED model.

  int  i = *iArg ;
  char fnam[] = "parse_input_GENMODEL_ARGLIST" ;
  
  // ---------- BEGIN ------------

  if ( i < 0 ) {
    // read from file 
    sprintf(c1err,"GENMODEL_ARGLIST valid only as command-line arg.");
    sprintf(c2err,"snlc_sim.exe <inFile> "
	    "GENMODEL_ARGLIST 'KEY1 VAL1 KEY2 VAL2 ...' ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    // readchar ( fp, &INPUTS.GEMMODEL_ARGLIST);
  }
  else {
    // read command-line arg
    // (C code scoops everything within quotes)
    i++ ;  sprintf(INPUTS.GENMODEL_ARGLIST, "%s", ARGV_LIST[i]);
    *iArg = i ;
  }

  return ;

} // end parse_input_GENMODEL_ARGLIST



// ==================================================
void parse_input_TAKE_SPECTRUM(FILE *fp, char *WARP_SPECTRUM_STRING) {

  // Created Mar 14 2017 
  // 
  // Works with SPECTROGRAPH to specify spectra w.r.t peakMJD,
  // instead of at aboslute MJD (simlib) values.
  //
  // *string1 = TREST([Tmin]:[Tmax])  or
  //            TOBS([Tmin]:[Tmax])   or
  //            HOST
  //
  // *string2 = SNR_ZPOLY([a0],[a1],[a2],,,,)  or
  //            TEXPOSE_ZPOLY([a0],[a1],[a2],,,,) 
  //       = a0 + a1*z + a2*z^2 + ...
  //
  // If WARP_SPECTRUM_STRING is not blank, parse this option
  // to apply mis-calibration vs. wavelength using polyFun of wavelength.
  //
  // Note that colon denotes a range, while comma separates a list.
  //
  // Mar 23 2019: 
  //  + use parse_GENPOLY to allow arbitrary poly-order.
  //  + pass WARP_SPECTRUM_STRING
  //
  // June 28 2019: allow HOST argument

  int  N = NPEREVT_TAKE_SPECTRUM ;
  GENPOLY_DEF *GENLAMPOLY_WARP  = &INPUTS.TAKE_SPECTRUM[N].GENLAMPOLY_WARP ;
  GENPOLY_DEF *GENZPOLY_TEXPOSE = &INPUTS.TAKE_SPECTRUM[N].GENZPOLY_TEXPOSE;
  GENPOLY_DEF *GENZPOLY_SNR     = &INPUTS.TAKE_SPECTRUM[N].GENZPOLY_SNR ;
  float *ptrF ;
  char string1[80], string2[80], string3[80]; 
  char *ptrSplit[4], strValues[4][20] ;
  int  NSPLIT, i, o ;
  int  IS_HOST = 0;
  char colon[] = ":", comma[] = "," ;
  char stringTmp[80], stringOpt[200];
  char fnam[] = "parse_input_TAKE_SPECTRUM" ;

  // ----------- BEGIN -----------

  // init TAKE_SPECTRUM structure 
  for(i=0; i < 2; i++)  { 
    INPUTS.TAKE_SPECTRUM[N].EPOCH_RANGE[i]  = -99.0;
    INPUTS.TAKE_SPECTRUM[N].SNR_LAMRANGE[i] = 0.0 ;  
  }

  init_GENPOLY(GENLAMPOLY_WARP  ) ; 
  init_GENPOLY(GENZPOLY_TEXPOSE ) ; 
  init_GENPOLY(GENZPOLY_SNR     ) ; 

  INPUTS.TAKE_SPECTRUM[N].OPT_FRAME_EPOCH  = 0 ;
  INPUTS.TAKE_SPECTRUM[N].OPT_FRAME_LAMBDA = 0 ;
  INPUTS.TAKE_SPECTRUM[N].OPT_TEXPOSE      = 0 ;

  
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
    parse_GENPOLY(stringOpt, GENLAMPOLY_WARP, fnam);
    INPUTS.NWARP_TAKE_SPECTRUM++ ;
  }
  


  // ----------------------------------------------
  // read 1st arg and parse as either TREST or TOBS
  readchar(fp,string1);
  sprintf(stringTmp, "%s", string1);
  extractStringOpt(stringTmp,stringOpt); // return stringOpt; 
  if ( strcmp(stringTmp,"TREST") == 0 ) {
    INPUTS.TAKE_SPECTRUM[N].OPT_FRAME_EPOCH = GENFRAME_REST ;
    sprintf(INPUTS.TAKE_SPECTRUM[N].EPOCH_FRAME,"REST");
  }
  else if ( strcmp(stringTmp,"TOBS") == 0 ) {
    INPUTS.TAKE_SPECTRUM[N].OPT_FRAME_EPOCH = GENFRAME_OBS ;
    sprintf(INPUTS.TAKE_SPECTRUM[N].EPOCH_FRAME,"OBS");
  }
  else if ( strcmp(stringTmp,"HOST") == 0 ) {
    INPUTS.TAKE_SPECTRUM[N].OPT_FRAME_EPOCH = GENFRAME_HOST ;
    sprintf(INPUTS.TAKE_SPECTRUM[N].EPOCH_FRAME,"HOST");
    IS_HOST = 1;
    INPUTS.NHOST_TAKE_SPECTRUM++ ;
  }
  else if ( strcmp(stringTmp,"TEMPLATE_TEXPOSE_SCALE") == 0 ) {
    sscanf(stringOpt, "%f", 
	   &INPUTS.TAKE_SPECTRUM_TEMPLATE_TEXPOSE_SCALE);
    return ;
  }
  else {
    sprintf(c1err, "Cannot parse '%s' after TAKE_SPECTRUM key.",string1);
    sprintf(c2err, "Expecting TREST or TOBS string" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

 

  // get epoch range from colon-seperated values in stringOpt
  ptrF = INPUTS.TAKE_SPECTRUM[N].EPOCH_RANGE ;

  if ( IS_HOST ) {
    ptrF[0] = ptrF[1] = 9999.0 ;  
  }
  else {
    for(i=0; i < 4; i++ ) { ptrSplit[i] = strValues[i]; }

    splitString(stringOpt, colon, 4,      // inputs               
		&NSPLIT, ptrSplit );      // outputs             

    if ( NSPLIT != 2 ) {
      sprintf(c1err, "\n   Found %d colon-separated values in '%s'", 
	      NSPLIT, string1);
      sprintf(c2err, "but expected %d values.", 2);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }   
    sscanf( strValues[0] , "%f", &ptrF[0] );  // load TREST_RANGE or TOBS_RANGE
    sscanf( strValues[1] , "%f", &ptrF[1] ); 
  }


  // -------------------------------
  //  parse string2 for SNR or TEXPOSE
  // -------------------------------

  readchar(fp,string2);
  sprintf(stringTmp, "%s", string2);
  extractStringOpt(stringTmp,stringOpt); // return stringOpt
  if ( strcmp(stringTmp,"SNR_ZPOLY") == 0 ) {
    INPUTS.TAKE_SPECTRUM[N].OPT_TEXPOSE = 2;
    parse_GENPOLY(stringOpt, GENZPOLY_SNR, fnam);
  }
  else if ( strcmp(stringTmp,"TEXPOSE_ZPOLY") == 0 ) {
    INPUTS.TAKE_SPECTRUM[N].OPT_TEXPOSE = 1;
    parse_GENPOLY(stringOpt, GENZPOLY_TEXPOSE, fnam);
  }
  else {
    sprintf(c1err, "Cannot parse '%s' after TAKE_SPECTRUM key.",string1);
    sprintf(c2err, "Expecting SNR_ZPOLY or TEXPOSE_ZPOLY string" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }



  // - - - - - - - - - - - - 
  // for SNR option, also read SNR_LAMREST or SNR_LAMOBS
  
  if ( INPUTS.TAKE_SPECTRUM[N].OPT_TEXPOSE == 2 ) {
    readchar(fp,string3);
    sprintf(stringTmp, "%s", string3);
    extractStringOpt(stringTmp,stringOpt); // return stringOpt

    splitString(stringOpt, colon, 5,      // inputs               
		&NSPLIT, ptrSplit );      // outputs             

    if ( strcmp(stringTmp,"SNR_LAMREST") == 0 ) 
      { INPUTS.TAKE_SPECTRUM[N].OPT_FRAME_LAMBDA = GENFRAME_REST ; }
    else if ( strcmp(stringTmp,"SNR_LAMOBS") == 0 ) 
      { INPUTS.TAKE_SPECTRUM[N].OPT_FRAME_LAMBDA = GENFRAME_OBS ; }
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

    ptrF = INPUTS.TAKE_SPECTRUM[N].SNR_LAMRANGE ;
    sscanf( strValues[0] , "%f", &ptrF[0] );  // load LAMRANGE
    sscanf( strValues[1] , "%f", &ptrF[1] ); 
  }    

  // - - - - - -

  int LDMP = 0 ;
  if ( LDMP ) {
    printf(" xxx ----------------------- \n");
    printf(" xxx NPEREVT_TAKE_SPECTRUM = %d \n", NPEREVT_TAKE_SPECTRUM);
    printf(" xxx EPOCH_RANGE = %6.2f to %6.2f \n", 
	   INPUTS.TAKE_SPECTRUM[N].EPOCH_RANGE[0],
	   INPUTS.TAKE_SPECTRUM[N].EPOCH_RANGE[1] );

    fflush(stdout);
  }

  NPEREVT_TAKE_SPECTRUM++ ;

  if ( NPEREVT_TAKE_SPECTRUM >= MXPEREVT_TAKE_SPECTRUM ) {
    sprintf(c1err, "%d TAKE_SPECTRUM keys exceeds bound", 
	    NPEREVT_TAKE_SPECTRUM );
    sprintf(c2err,"Check MXPEREVT_TAKE_SPECTRUM in snlc_sim.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return;

} // end parse_input_TAKE_SPECTRUM


// *****************************************
void read_input_SIMSED(FILE *fp, char *KEY) {

  // Created Feb 2018
  // parse SIMSED_XXX keys
  //
  // Apr 28 2019: IPAR_SIMSED_XXX -= 1 to account for earlier
  //              refactor where IPAR starts at zero.

  int N;
  char ctmp[80];
  char fnam[] = "read_input_SIMSED" ;

  // --------------- BEGIN ----------------

  if ( strcmp(KEY,"SIMSED_USE_BINARY:" ) == 0 )
    { readint ( fp, 1, &INPUTS.USE_BINARY_SIMSED ); return; }

  if ( strcmp(KEY,"SIMSED_PATH_BINARY:" ) == 0 )
    { readchar ( fp, INPUTS.PATH_BINARY_SIMSED );  return; }

  // check for SIMSED_PARAM keys
  if ( strcmp(KEY,"SIMSED_PARAM:" ) == 0 ) 
    { read_input_SIMSED_PARAM(fp); return ; }

  if ( strcmp(KEY,"SIMSED_SHAPEPAR:" ) == 0 )  { 
    read_input_SIMSED_PARAM(fp); 
    INPUTS.IPAR_SIMSED_SHAPE = INPUTS.NPAR_SIMSED-1 ;
    return ;
  }

  if ( strcmp(KEY,"SIMSED_COLORPAR:" ) == 0 )  { 
    read_input_SIMSED_PARAM(fp); 
    INPUTS.IPAR_SIMSED_COLOR = INPUTS.NPAR_SIMSED-1 ;
  }
  if ( strcmp(KEY,"SIMSED_COLORLAW:" ) == 0 )  { 
    read_input_SIMSED_PARAM(fp); 
    INPUTS.IPAR_SIMSED_CL = INPUTS.NPAR_SIMSED-1 ;
    return ;
  }
  
  if ( strcmp(KEY,"SIMSED_GRIDONLY:" ) == 0 ) {
    
    readchar(fp, ctmp );
    if ( strcmp(ctmp,"SEQUENTIAL") == 0  ) {
      INPUTS.NPAR_SIMSED = 0 ;  // no interp pars => treat as baggage
      INPUTS.OPTMASK_SIMSED    = OPTMASK_SIMSED_GRIDONLY ;
    }
    else {
      N = INPUTS.NPAR_SIMSED ;
      sprintf(INPUTS.PARNAME_SIMSED[N], "%s", ctmp );
      INPUTS.GENFLAG_SIMSED[N] = 
	OPTMASK_SIMSED_PARAM + OPTMASK_SIMSED_GRIDONLY ;
      sprintf(INPUTS.KEYWORD_SIMSED[N],"SIMSED_GRIDONLY");

      INPUTS.NPAR_SIMSED++ ; 
      INPUTS.NPAR_SIMSED_GRIDONLY++ ; 
    }
    return ;
  } // end SIMSED_SEQUENTIAL
  

  // - - - - - - - - - - - 
  // check for SIMSED_REDCOR(varname1,varname2):  <redCor>  
  char string[60], stringOpt[60];
  sprintf(string,"%s", KEY);
  extractStringOpt(string,stringOpt);

  if ( strcmp(string,"SIMSED_REDCOR:") == 0 ) 
    { read_input_SIMSED_COV(fp, 1, stringOpt);  }
  if ( strcmp(string,"SIMSED_COV:") == 0 ) 
    { read_input_SIMSED_COV(fp, 2, stringOpt);  }

  //  debugexit(fnam);

  return ;

} // end read_input_SIMSED


// *******************************************
void read_input_SIMSED_COV(FILE *fp, int OPT, char *varList ) {

  // Feb 18 2018
  // Store information for off-diagonal covariances used
  // to generate correlated random Gaussians.
  // Called if "SIMSED_REDCOR(var1,var2):" is found.
  // varList = var1,var2 needs to be parsed.
  //
  // Inputs:
  //   fp = file pointer to read
  //   OPT = 1 for reduced correlation, =2 for COV
  //   varList = comma-separated list of two correlated varNames
  //

#define OPTREAD_REDCOR_SIMSED 1
#define OPTREAD_COV_SIMSED    2

  int NPAIR = INPUTS.NPAIR_SIMSED_COV; // == Noffdiag/2
  double sigLo, sigHi, SIGLIST[2], RHO;
  int  MXNAME=4, NNAME, ipar, IPAR[2], IPAR0, IPAR1, IPARTMP, j, NFOUND=0 ;
  char comma[] = "," ;
  char *ptrName[2], varNames[2][40];
  char fnam[] = "read_input_SIMSED_COV" ;

  // -------------- BEGIN --------------

  if ( NPAIR >= MXPAR_SIMSED ) { goto COUNT; }

  readfloat(fp, 1, &INPUTS.COVPAIR_SIMSED_COV[NPAIR]) ;  // either REDCOR or COV

    // break stringOpt into two SIMSED varnames
  ptrName[0] = varNames[0] ;
  ptrName[1] = varNames[1] ;
  splitString(varList, comma, MXNAME, &NNAME, ptrName);

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

  return ;

} // end read_input_SIMSED_COV


// *****************************************
void read_input_SIMSED_PARAM(FILE *fp) {

  // Jun 2013
  // call this function when any of the
  // SIMSED_[PARAM/SHAPEPAR/COLORPAR] keys is found.

  //  char fnam[] = "read_input_SIMSED_PARAM" ;

  // ------------- BEGIN -----------

  int N ;
  N = INPUTS.NPAR_SIMSED ;
  readchar(fp, INPUTS.PARNAME_SIMSED[N] );
  INPUTS.GENFLAG_SIMSED[N] =  OPTMASK_SIMSED_PARAM ; // continuous interp
  sprintf(INPUTS.KEYWORD_SIMSED[N],"SIMSED_PARAM");
  
  INPUTS.NPAR_SIMSED++ ; 
  INPUTS.NPAR_SIMSED_PARAM++ ; 

  return ;

} // end of read_input_SIMSED_PARAM

// **************************************************
void read_input_GENGAUSS(FILE *fp, char *string, char *varName, 
			 GENGAUSS_ASYM_DEF *genGauss ) {

  // May 3, 2012
  // Utility to read GENGAUSS struct from input file.
  // for input arg *string that has been read from file fp,
  // check if it matches any of the GENXXX_[varName] keys;
  // if so then load output genGauss struct.
  // Input *varName is the name of the variable: e.g. 'SALT2c'
  //
  // Mar 16 2014: allow  GENMEAN or GENPEAK
  // Apr 21 2016: check for NGRID; 
  // Aug 30 2016: check for SKEWNORMAL, and fill genGauss->NAME
  // Mar 29 2017: check for 2nd peak

  int FOUND=0 ;
  char KEYNAME[80];
  char fnam[] = "read_input_GENGAUSS" ;

  // ----------------- BEGIN ----------------

  sprintf(KEYNAME,  "GENPEAK_%s:",  varName);
  if ( strcmp(string, KEYNAME )==0 ) 
    { readdouble ( fp, 1, &genGauss->PEAK ); FOUND=1; }

  sprintf(KEYNAME,  "GENMEAN_%s:",  varName);  
  if ( strcmp(string, KEYNAME )==0 ) 
    { readdouble ( fp, 1, &genGauss->PEAK ); FOUND=1; }

  sprintf(KEYNAME, "GENSIGMA_%s:", varName);
  if ( strcmp(string, KEYNAME )==0 ) {
    readdouble ( fp, 2, genGauss->SIGMA ); FOUND=1; 
    checkVal_GENGAUSS(KEYNAME, genGauss->SIGMA, fnam );
  }

  sprintf(KEYNAME,  "GENSKEW_%s:",  varName);
  if ( strcmp(string, KEYNAME )==0 ) 
    { readdouble ( fp, 2, genGauss->SKEW ); FOUND=1; }


  sprintf(KEYNAME, "GENSKEWNORMAL_%s:",  varName);
  if ( strcmp(string, KEYNAME )==0 ) 
    { readdouble ( fp, 3, genGauss->SKEWNORMAL );  FOUND=1; }

  sprintf(KEYNAME, "GENRANGE_%s:", varName);
  if ( strcmp(string, KEYNAME )==0 ) {
    readdouble ( fp, 2, genGauss->RANGE ); FOUND=1; 
    checkVal_GENGAUSS(KEYNAME, genGauss->RANGE, fnam );
  }

  // read NGRID 
  sprintf(KEYNAME, "GENGRID_%s:",  varName); // NGRID = number of grid bins
  if ( strcmp(string, KEYNAME )==0 ) {
    int NGRID ; readint ( fp, 1, &NGRID ); FOUND=1;
    genGauss->NGRID = NGRID ; 
  }

  // check for 2nd peak
  sprintf(KEYNAME,  "GENPROB2_%s:",  varName);
  if ( strcmp(string, KEYNAME )==0 ) 
    { readdouble ( fp, 1, &genGauss->PROB2 ); FOUND=1; }
  sprintf(KEYNAME,  "GENPEAK2_%s:",  varName);
  if ( strcmp(string, KEYNAME )==0 ) 
    { readdouble ( fp, 1, &genGauss->PEAK2 ); FOUND=1; }
  sprintf(KEYNAME, "GENSIGMA2_%s:", varName);
  if ( strcmp(string, KEYNAME )==0 ) {
    readdouble ( fp, 2, genGauss->SIGMA2 ); FOUND=1; 
    checkVal_GENGAUSS(KEYNAME, genGauss->SIGMA2, fnam );
  }

  // - - - - - -
  if ( FOUND ) { prepIndex_GENGAUSS(varName, genGauss); }

  return ;

} // end of read_input_GENGAUSS

void prepIndex_GENGAUSS(char *varName, GENGAUSS_ASYM_DEF *genGauss ) {

  // Created Sep 2 2016
  // Store NAME and increment index.
  // Called by readFile routine and command-line read function.

  char *ptrName = genGauss->NAME;
  //  char fnam[] = "prepIndex_GENGAUSS" ;
  // ---------- BEGIN ---------

  /*
  printf(" xxx FOUND GENGAUSS(%s): peak=%.2f  "
	 "range=%9.3le,%9.3le  sigma=%9.3le,%9.3le \n",
	 varName, genGauss->PEAK, 
	 genGauss->RANGE[0], genGauss->RANGE[1],
	 genGauss->SIGMA[0], genGauss->SIGMA[1] );  fflush(stdout);
  */

  // if genGauss name is not set, then set name and FUNINDEX
  if ( strcmp(ptrName,varName) != 0 ) {
    genGauss->FUNINDEX = NFUN_GENGAUSS_ASYM ;
    NFUN_GENGAUSS_ASYM++;  
    sprintf(genGauss->NAME, "%s", varName);  
  }

  // copy each GENGAUSS_ASYM struct into master list in case
  // some kind of global operation or init is needed.

  int FUNINDEX = genGauss->FUNINDEX ;
  copy_GENGAUSS_ASYM( genGauss, &GENGAUSS_ASYM_LIST[FUNINDEX] );

} // end prepIndex_GENGAUSS

// ************************************************
void sscanf_GENGAUSS(int *i, char *varName, GENGAUSS_ASYM_DEF *genGauss ) {

  // May 3, 2012
  // Utility to read GENGAUSS struct from command line.
  // *varName is the name of the variable, such as SALT2c, DELTA, RV ...
  // Load output genGauss struct.
  //
  // Mar 16 2014: allow  GENMEAN or GENPEAK
  // Apr 21 2016: allow  GENGRID

  int j, FOUND=0 ;
  char KEYNAME[80];
  char fnam[] = "sscanf_GENGAUSS" ;
  
  // ------------ BEGIN ----------

  j = *i; // current index of ARGV_LIST

  sprintf(KEYNAME,  "GENPEAK_%s", varName);
  if ( strcmp( ARGV_LIST[j], KEYNAME ) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->PEAK ); 
  }

  sprintf(KEYNAME, "GENMEAN_%s", varName);
  if ( strcmp( ARGV_LIST[j], KEYNAME ) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->PEAK ); 
  }


  sprintf(KEYNAME, "GENSIGMA_%s", varName);
  if ( strcmp( ARGV_LIST[*i], KEYNAME) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SIGMA[0] ); 
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SIGMA[1] ); 

    // idiot check
    checkVal_GENGAUSS(KEYNAME, genGauss->SIGMA, fnam );
  }


  sprintf(KEYNAME, "GENSKEW_%s", varName);
  if ( strcmp( ARGV_LIST[*i], KEYNAME ) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SKEW[0] ); 
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SKEW[1] ); 
  }

  sprintf(KEYNAME, "GENSKEWNORMAL_%s",  varName);
  if ( strcmp( ARGV_LIST[*i], KEYNAME ) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SKEWNORMAL[0] ); 
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SKEWNORMAL[1] ); 
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SKEWNORMAL[2] ); 
  }

  sprintf(KEYNAME, "GENRANGE_%s", varName);
  if ( strcmp( ARGV_LIST[*i], KEYNAME ) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->RANGE[0] ); 
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->RANGE[1] ); 

    // idiot check
    checkVal_GENGAUSS(KEYNAME, genGauss->RANGE, fnam );
  }

  // Apr 2016: check for GENGRID_XXX
  sprintf(KEYNAME, "GENGRID_%s", varName);
  if ( strcmp( ARGV_LIST[*i], KEYNAME ) == 0 ) {
    int NGRID; 
    j++ ; sscanf(ARGV_LIST[j] , "%d", &NGRID );
    // xxx mark delete  if ( NGRID>=2 ) { genGauss->NGRID = NGRID ; }
    genGauss->NGRID = NGRID ; 
  }

  // Jun 2017: check for 2nd Gaussian
  sprintf(KEYNAME,  "GENPROB2_%s", varName);
  if ( strcmp( ARGV_LIST[j], KEYNAME ) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->PROB2 ); 
  }
  sprintf(KEYNAME,  "GENPEAK2_%s", varName);
  if ( strcmp( ARGV_LIST[j], KEYNAME ) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->PEAK2 ); 
  }
  sprintf(KEYNAME, "GENSIGMA2_%s", varName);
  if ( strcmp( ARGV_LIST[*i], KEYNAME) == 0 ) {
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SIGMA2[0] ); 
    j++ ; sscanf(ARGV_LIST[j] , "%le", &genGauss->SIGMA2[1] ); 

    // idiot check
    checkVal_GENGAUSS(KEYNAME, genGauss->SIGMA, fnam );
  }


  // check if any key was found.
  FOUND = ( j > *i ) ;
  if ( FOUND ) { prepIndex_GENGAUSS(varName, genGauss); }

  *i = j;

  return ;

} // end of sscanf_GENGAUSS


// *************************************
void checkVal_GENGAUSS(char *varName, double *val, char *fromFun ) {

  // Mar 16 2014
  // abort if any GENGAUSS values are invalid, 
  // such as negative GENSIGMA

  double lo, hi ;
  char fnam[] = "checkVal_GENGAUSS" ;

  // ------------- BEGIN ----------

  if ( strstr(varName,"RANGE") != NULL ) {
    lo  = val[0] ;    hi  = val[1] ;
    if ( hi < lo  ) {
      sprintf(c1err, "Invalid %s range pass from %s", varName, fromFun); 
      sprintf(c2err, "Specified  %s  range is %f to %f\n",
	      varName, lo, hi);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    } 
  }


  if ( strstr(varName,"SIGMA") != NULL ) {
    lo  = val[0] ;    hi  = val[1] ;
    if ( lo < 0.0 || hi < 0.0  ) {
      sprintf(c1err, "Invalid %s = %.3f, %.3f passed from %s", 
	      varName, lo, hi, fromFun); 
      sprintf(c2err, "%s must be both positive in sim-input file.", 
	      varName);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    } 
  }

  return ;

}  // end of checkVal_GENGAUSS



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

  /* xxxx mark delete Apr 18 2019 xxxxxxx
  //  if ( strcmp(GENMODEL,"LCLIB" ) == 0 )  { return(MODEL_LCLIB); }
  xxxx end mark xxxxxx*/

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

  // Check for command-line overrides.
  // argv[0] is the name of the program.
  // argv[1] is the name of the input file.
  // start parsing at argv[2].
  //
  // Dec  3 2013: fix to work with DNDZ POWERLAW2
  // Aug 13 2016: refactor to call sscanf_RATEPAR().
  // Mar 10 2017: start add 'goto INCREMENT_COUNTER' to avoid conflicts.
  // Dec 27 2018: check LCLIB_CUTWIN

  int  i, ilast, iuse, ifilt, ifilt_obs, N, j, itmp ;
  int  ipar, opt_tmp, itype, NLOCAL_DNDZ ;
  int  NTMP, ifilt_list[MXFILTINDX] ;
  char  ctmp[80], ctmp2[80], cfilt[2], parname[60] ;
  float ftmp, old, tmpSmear[2];
  char *modelName ;
  char comma[] = "," ;
  char fnam[] = "sim_input_override" ;
  FILE *fpNull ;

  // ------------ BEGIN -----------

  i = 2; ilast = 2 ;
  NLOCAL_DNDZ = 0 ;

  while ( i < NARGV_LIST ) {
    printf("  PROCESS COMMAND LINE ARG: %s \n", ARGV_LIST[i] );

    if ( strcmp( ARGV_LIST[i], "INPUT_FILE_INCLUDE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.INPUT_FILE_LIST[1] ); 
      INPUTS.INPUT_FILE_LIST[2][0] = 0 ; // erase 2nd include file
    }

    if ( strcmp( ARGV_LIST[i], "COMMENT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.COMMENT ); 
    }

    if ( strcmp( ARGV_LIST[i], "TRACE_MAIN" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.TRACE_MAIN ); 
    }
    if ( strcmp( ARGV_LIST[i], "OPT_DEVEL_BBC7D" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.OPT_DEVEL_BBC7D ); 
    }
    if ( strcmp( ARGV_LIST[i], "DEBUG_FLAG" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.DEBUG_FLAG ); 
    }

    if ( strcmp( ARGV_LIST[i], "NONLINEARITY_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.NONLINEARITY_FILE ); 
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.SIMLIB_FILE ); 
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_FIELDLIST" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.SIMLIB_FIELDLIST ); 
    }


    if ( strcmp( ARGV_LIST[i], "ZVARIATION_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUT_ZVARIATION_FILE ); 
      if ( !IGNOREFILE(INPUT_ZVARIATION_FILE) )
	{ USE_ZVAR_FILE = 1 ; }
    }

    if ( strcmp( ARGV_LIST[i], "ZVARIATION_POLY" ) == 0 ) {
      N = NPAR_ZVAR_USR ;
      INPUT_ZVARIATION[N].FLAG = FLAG_ZPOLY_ZVAR ;
      INPUT_ZVARIATION[N].NZBIN = 0 ;  
      NPAR_ZVAR_USR++ ;  
      i++ ; sscanf(ARGV_LIST[i] , "%s",  INPUT_ZVARIATION[N].PARNAME );
      i++ ; sscanf(ARGV_LIST[i] , "%s", ctmp );
      parse_GENPOLY(ctmp, &INPUT_ZVARIATION[N].POLY, fnam );

      /* xxxx mark delete Apr 10 2019 xxxxxxxxx
      for(j=0; j <= POLYORDER_ZVAR; j++ ) {
	i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUT_ZVARIATION[N].ZPOLY[j] );
      }
      xxxxxxxxxxxxxxxxxx */

    }  // ZVARIATION_POLY


    // ------ HOSTLIB -----
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.HOSTLIB_FILE ); 
      INPUTS.HOSTLIB_USE = 1;  
      setbit_HOSTLIB_MSKOPT(HOSTLIB_MSKOPT_USE) ;

      // check if we are just turning off the HOSTLIB
      if ( IGNOREFILE(INPUTS.HOSTLIB_FILE)  )
	{ INPUTS.HOSTLIB_USE = INPUTS.HOSTLIB_MSKOPT = 0 ; }

    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_WGTMAP_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.HOSTLIB_WGTMAP_FILE ); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_ZPHOTEFF_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.HOSTLIB_ZPHOTEFF_FILE ); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_SPECBASIS_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.HOSTLIB_SPECBASIS_FILE ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_MSKOPT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.HOSTLIB_MSKOPT ); 
      setbit_HOSTLIB_MSKOPT(HOSTLIB_MSKOPT_USE) ;
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GENZPHOT_FUDGEPAR" ) == 0 ) {
      for(j=0; j<4; j++ ) {
	i++ ; sscanf(ARGV_LIST[i], "%f", 
		     &INPUTS.HOSTLIB_GENZPHOT_FUDGEPAR[j] );	
      }      
      i++ ; sscanf(ARGV_LIST[i], "%s", ctmp); 
      parse_input_GENZPHOT_OUTLIER(ctmp); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GENZPHOT_BIAS" ) == 0 ) {
      for(j=0; j<4; j++ ) {
	i++ ; sscanf(ARGV_LIST[i] , "%f", 
		     &INPUTS.HOSTLIB_GENZPHOT_BIAS[j] );
      }      
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_MAXREAD" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.HOSTLIB_MAXREAD ); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GALID_NULL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.HOSTLIB_GALID_NULL ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GALID_PRIORITY" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.HOSTLIB_GALID_PRIORITY[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.HOSTLIB_GALID_PRIORITY[1] ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_MINDAYSEP_SAMEGAL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.HOSTLIB_MINDAYSEP_SAMEGAL ); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_MXINTFLUX_SNPOS" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.HOSTLIB_MXINTFLUX_SNPOS ); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GENRANGE_NSIGZ" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.HOSTLIB_GENRANGE_NSIGZ[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.HOSTLIB_GENRANGE_NSIGZ[1] ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GENRANGE_RA" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_GENRANGE_RA[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_GENRANGE_RA[1] ); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GENRANGE_DEC" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_GENRANGE_DEC[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_GENRANGE_DEC[1] ); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GENRANGE_DECL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_GENRANGE_DEC[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_GENRANGE_DEC[1] ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_DZTOL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_DZTOL[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_DZTOL[1] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_DZTOL[2] ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_STOREPAR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.HOSTLIB_STOREPAR_LIST ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_SERSIC_SCALE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_SERSIC_SCALE ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_SBRADIUS" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_SBRADIUS ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_GALID_FORCE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.HOSTLIB_GALID_FORCE ); 
    }

    if ( strcmp( ARGV_LIST[i], "HOSTLIB_FIXRAN_RADIUS" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_FIXRAN_RADIUS ); 
    }
    if ( strcmp( ARGV_LIST[i], "HOSTLIB_FIXRAN_PHI" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.HOSTLIB_FIXRAN_PHI ); 
    }

    if ( strcmp( ARGV_LIST[i], "FLUXERRMODEL_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.FLUXERRMODEL_FILE ); 
    }
    if ( strcmp( ARGV_LIST[i], "FLUXERRMODEL_OPTMASK" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.FLUXERRMODEL_OPTMASK ); 
    }
    if ( strcmp( ARGV_LIST[i], "FLUXERRMAP_IGNORE_DATAERR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.FLUXERRMAP_IGNORE_DATAERR ); 
    }

    // anomalous host-subtraction noise (Aug 2014)
    if ( strcmp( ARGV_LIST[i], "HOSTNOISE_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.HOSTNOISE_FILE ); 
    }

    // --------------

    if ( strcmp( ARGV_LIST[i], "SNTYPE_Ia" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &itype );
      INPUTS.SNTYPE_Ia_SPEC = itype ;
      INPUTS.SNTYPE_Ia_PHOT = itype + OFFSET_TYPE_PHOT ;
    }
    if ( strcmp( ARGV_LIST[i], "SNTYPES_Ia" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SNTYPE_Ia_SPEC );
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SNTYPE_Ia_PHOT );
    }

    if ( strcmp( ARGV_LIST[i], "GENTYPE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &itype );
      INPUTS.GENTYPE_SPEC = itype ;
      INPUTS.GENTYPE_PHOT = itype + OFFSET_TYPE_PHOT ;
    }
    if ( strcmp( ARGV_LIST[i], "GENTYPES" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENTYPE_SPEC );
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENTYPE_PHOT );
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_IDSTART" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SIMLIB_IDSTART ); 
    }
    if ( strcmp( ARGV_LIST[i], "SIMLIB_MAXRANSTART" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SIMLIB_MAXRANSTART ); 
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_IDLOCK" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SIMLIB_IDLOCK ); 
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_MINOBS" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SIMLIB_MINOBS ); 
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_MINSEASON" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.SIMLIB_MINSEASON ); 
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_DUMP" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SIMLIB_DUMP ); 
    }
    if ( strcmp( ARGV_LIST[i], "SIMLIB_CADENCEFOM_ANGSEP" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", 
		   &INPUTS.SIMLIB_CADENCEFOM_ANGSEP ); 
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_CADENCEFOM_PARLIST" ) == 0 ) {
      for ( j=0; j < NPAR_SNCADENCEFOM; j++ ) { 
	  i++ ; sscanf(ARGV_LIST[i] ,"%le",
		       &INPUTS.SIMLIB_CADENCEFOM_PARLIST[j] ); 
      }
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_NREPEAT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SIMLIB_NREPEAT ); 
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_NSKIPMJD" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SIMLIB_NSKIPMJD ); 
    }

    if ( strcmp( ARGV_LIST[i], "USE_SIMLIB_REDSHIFT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.USE_SIMLIB_REDSHIFT ); 
      INPUTS.USE_SIMLIB_GENOPT=1;
    }
    if ( strcmp( ARGV_LIST[i], "USE_SIMLIB_DISTANCE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.USE_SIMLIB_DISTANCE ); 
      INPUTS.USE_SIMLIB_GENOPT=1;
    }
    if ( strcmp( ARGV_LIST[i], "USE_SIMLIB_PEAKMJD" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.USE_SIMLIB_PEAKMJD ); 
      INPUTS.USE_SIMLIB_GENOPT=1;
    }

    if ( strcmp( ARGV_LIST[i], "SIMLIB_MSKOPT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SIMLIB_MSKOPT );       
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "FUDGEOPT_FLUXERR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.FUDGEOPT_FLUXERR ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "FUDGESCALE_PSF" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.FUDGESCALE_PSF ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "FUDGESCALE_SKYNOISE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.FUDGESCALE_SKYNOISE ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "FUDGESCALE_READNOISE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.FUDGESCALE_READNOISE ); 
      goto INCREMENT_COUNTER; 
    }


    // check filter-dependent lists
    parse_input_KEY_PLUS_FILTER(fpNull, &i, ARGV_LIST[i], "EXPOSURE_TIME",
				&INPUTS.EXPOSURE_TIME, 
				INPUTS.EXPOSURE_TIME_FILTER);

    parse_input_KEY_PLUS_FILTER(fpNull, &i, ARGV_LIST[i], "FUDGE_MAGERR",
				&INPUTS.FUDGE_MAGERR, 
				INPUTS.FUDGE_MAGERR_FILTER);

    parse_input_KEY_PLUS_FILTER(fpNull, &i, ARGV_LIST[i], "FUDGESCALE_FLUXERR",
				&INPUTS.FUDGESCALE_FLUXERR, 
				INPUTS.FUDGESCALE_FLUXERR_FILTER );

    parse_input_KEY_PLUS_FILTER(fpNull, &i, ARGV_LIST[i],"FUDGESCALE_FLUXERR2",
				&INPUTS.FUDGESCALE_FLUXERR2, 
				INPUTS.FUDGESCALE_FLUXERR2_FILTER );

    parse_input_KEY_PLUS_FILTER(fpNull, &i, ARGV_LIST[i], "FUDGESHIFT_ZPT",
				&INPUTS.FUDGESHIFT_ZPT, 
				INPUTS.FUDGESHIFT_ZPT_FILTER);

    parse_input_KEY_PLUS_FILTER(fpNull, &i, ARGV_LIST[i], "MJD_TEMPLATE",
				&INPUTS.MJD_TEMPLATE, 
				INPUTS.MJD_TEMPLATE_FILTER);
    
    if ( strstr(ARGV_LIST[i],"RANSYSTPAR") != NULL ) 
      { parse_input_RANSYSTPAR(fpNull, &i, ARGV_LIST[i] ); }


    // ----

    if ( strcmp( ARGV_LIST[i], "INIT_ONLY" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.INIT_ONLY ); 
      goto INCREMENT_COUNTER ; 
    }

    if ( strcmp( ARGV_LIST[i], "NGEN_LC" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.NGEN_LC ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "NGENTOT_LC" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.NGENTOT_LC ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "NGEN_SEASON" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.NGEN_SEASON ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "NGEN_SCALE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.NGEN_SCALE ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "NGEN_SCALE_NON1A" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.NGEN_SCALE_NON1A ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "NSUBSAMPLE_MARK" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.NSUBSAMPLE_MARK ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "CIDOFF" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.CIDOFF ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "CIDRAN_MAX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.CIDRAN_MAX ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "CIDRAN_MIN" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.CIDRAN_MIN ); 
      goto INCREMENT_COUNTER; 
    }

    // JOBID & NJOBTOT can be read only from command-line;
    // NOT from sim-input file.
    if ( strcmp( ARGV_LIST[i], "JOBID" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.JOBID ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "NJOBTOT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.NJOBTOT ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "GENVERSION" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.GENVERSION ); 
      sprintf(INPUTS.GENPREFIX,"%s", INPUTS.GENVERSION);
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "GENPREFIX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.GENPREFIX ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "PATH_SNDATA_SIM" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.PATH_SNDATA_SIM );
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "CLEARPROMPT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.CLEARPROMPT ); 
      goto INCREMENT_COUNTER; 
    }


    if ( strcmp( ARGV_LIST[i], "GENSOURCE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.GENSOURCE ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "GENMODEL_MSKOPT" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENMODEL_MSKOPT ); }
    if ( strcmp( ARGV_LIST[i], "GENMODEL_ARGLIST" ) == 0 ) { 
      parse_input_GENMODEL_ARGLIST(fpNull,&i);
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "GENMODEL" ) == 0 ) {
      parse_input_GENMODEL(fpNull,&i);
      goto INCREMENT_COUNTER; 
    } // end GENMODEL key


    if ( strcmp( ARGV_LIST[i], "FORMAT_MASK" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.FORMAT_MASK ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "WRFLAG_MODELPAR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.WRFLAG_MODELPAR ); 
      goto INCREMENT_COUNTER; 
    }


    if ( strcmp( ARGV_LIST[i], "NPE_PIXEL_SATURATE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.NPE_PIXEL_SATURATE ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "PHOTFLAG_SATURATE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.PHOTFLAG_SATURATE ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "PHOTFLAG_SNRMAX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.PHOTFLAG_SNRMAX ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "PHOTFLAG_NEARPEAK" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.PHOTFLAG_NEARPEAK ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "PHOTFLAG_DETECT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS_SEARCHEFF.PHOTFLAG_DETECT);
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "PHOTFLAG_TRIGGER" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS_SEARCHEFF.PHOTFLAG_TRIGGER);
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "NON1A"    ) == 0 || 
	 strcmp( ARGV_LIST[i], "NON1ASED" ) == 0 ) {
      INPUTS.NON1ASED.NINDEX++ ; N = INPUTS.NON1ASED.NINDEX;
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.NON1ASED.INDEX[N] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.NON1ASED.WGT[N] );
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "PATH_NON1ASED" ) == 0 ||
	 strcmp( ARGV_LIST[i], "PATH_NONIASED" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.NON1ASED.PATH );
      goto INCREMENT_COUNTER; 
    }


    if ( strcmp(ARGV_LIST[i], "MJD_EXPLODE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.MJD_EXPLODE );  
      if(INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE < 0 ) 
	{ INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE=0; }
      INPUTS.GENRANGE_PEAKMJD[0] = INPUTS.MJD_EXPLODE ;
      INPUTS.GENRANGE_PEAKMJD[1] = INPUTS.MJD_EXPLODE ;
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "OPTMASK_T0SHIFT_EXPLODE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", 
		   &INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE );  
      goto INCREMENT_COUNTER; 
    }


    if ( strcmp( ARGV_LIST[i], "UVLAM_EXTRAPFLUX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le",&INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX );
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "MINSLOPE_EXTRAPMAG_LATE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le",
		   &INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE );
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "RANSEED" ) == 0 )  { 
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.ISEED );  
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "RANLIST_START_GENSMEAR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.RANLIST_START_GENSMEAR ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "GENRANGE_RA" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_RA[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_RA[1] ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "GENRANGE_DECL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_DEC[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_DEC[1] ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "GENRANGE_DEC" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_DEC[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_DEC[1] ); 
      goto INCREMENT_COUNTER; 
    }


    if ( strstr(ARGV_LIST[i],"SOLID_ANGLE") != NULL ) { 
      parse_input_SOLID_ANGLE(fpNull, &i, ARGV_LIST[i] ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "GENRANGE_REDSHIFT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_REDSHIFT[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_REDSHIFT[1] ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "GENSIGMA_REDSHIFT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENSIGMA_REDSHIFT ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "GENBIAS_REDSHIFT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENBIAS_REDSHIFT ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "GENSIGMA_VPEC" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENSIGMA_VPEC ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "VPEC_ERR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.VPEC_ERR ); 
      goto INCREMENT_COUNTER; 
    }
    if ( strcmp( ARGV_LIST[i], "VEL_CMBAPEX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.VEL_CMBAPEX ); 
      goto INCREMENT_COUNTER; 
    }

    // - - -  wrongHost model params
    if ( strcmp( ARGV_LIST[i], "WRONGHOST_FILE")==0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%s", INPUTS.WRONGHOST_FILE );
      goto INCREMENT_COUNTER; 
    }

    // - - - - - -  DNDZ stuff - - - - - - - - - - - -

    if ( strcmp( ARGV_LIST[i], "DNDZ_ZEXP_REWGT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", 
		   &INPUTS.RATEPAR.DNDZ_ZEXP_REWGT ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "DNDZ_ZPOLY_REWGT" ) == 0 ) {
      for(j=0; j<4; j++ ) 
	{ i++ ; sscanf(ARGV_LIST[i] , "%le", 
		       &INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT[j]);  }
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "DNDZ_SCALE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.RATEPAR.DNDZ_SCALE[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.RATEPAR.DNDZ_SCALE[1] ); 
      goto INCREMENT_COUNTER; 
    }

    if ( strcmp( ARGV_LIST[i], "DNDZ_ALLSCALE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.RATEPAR.DNDZ_ALLSCALE ); 
    }

    // LEGACY KEY
    if ( strcmp( ARGV_LIST[i], "DNDZ_SCALE_NON1A" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.RATEPAR.DNDZ_SCALE[1] ); 
    }

    // -----------------------------------------------
    // check for rate model following DNDZ key
    sscanf_RATEPAR(&i, "NOMINAL",  &INPUTS.RATEPAR );
    sscanf_RATEPAR(&i, "PEC1A",    &INPUTS.RATEPAR_PEC1A );

    // ------------------------------------------------
    // read RISETIME shift params

    sscanf_GENGAUSS(&i, "RISETIME_SHIFT", &INPUTS.GENGAUSS_RISETIME_SHIFT ) ;
    sscanf_GENGAUSS(&i, "FALLTIME_SHIFT", &INPUTS.GENGAUSS_FALLTIME_SHIFT ) ;

    // ---

    if ( strcmp( ARGV_LIST[i], "GENRANGE_PEAKMAG" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_PEAKMAG[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_PEAKMAG[1] ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENRANGE_MJD" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_MJD[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_MJD[1] ); 
    }
    if ( strcmp( ARGV_LIST[i], "GENRANGE_PEAKMJD" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_PEAKMJD[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_PEAKMJD[1] ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENSIGMA_SEARCH_PEAKMJD" ) == 0 ) { //legacy
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENSIGMA_PEAKMJD); 
    }
    if ( strcmp( ARGV_LIST[i], "GENSIGMA_PEAKMJD" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENSIGMA_PEAKMJD); 
    }
    if ( strcmp( ARGV_LIST[i], "OPT_SETPKMJD" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.OPT_SETPKMJD); 
    }

    if ( strcmp( ARGV_LIST[i], "NEWMJD_DIF" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.NEWMJD_DIF ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENRANGE_TREST" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_TREST[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_TREST[1] ); 
    }

    if ( strcmp( ARGV_LIST[i], "TGRIDSTEP_MODEL_INTERP" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.TGRIDSTEP_MODEL_INTERP );
    }

    if ( strcmp( ARGV_LIST[i], "GENPAR_SELECT_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.GENPAR_SELECT_FILE );
    }

    // ------------------------------------------------
    // read STRETCH parameters
    
    sscanf_GENGAUSS(&i, "STRETCH", &INPUTS.GENGAUSS_STRETCH   ) ;

    if ( strcmp( ARGV_LIST[i], "STRETCH_TEMPLATE_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.STRETCH_TEMPLATE_FILE ); 
    }

    // ------------------------------------------------
    // read DELTA parameters
    sscanf_GENGAUSS(&i, "DELTA",  &INPUTS.GENGAUSS_DELTA  ) ;

    // ------------------------------------------------
    // read DM15 parameters
    sscanf_GENGAUSS(&i, "DM15",  &INPUTS.GENGAUSS_DM15  ) ;

    // ------------------------------------------------
    // read RV parameters
    sscanf_GENGAUSS(&i, "RV",  &INPUTS.GENGAUSS_RV  ) ;

    if ( strcmp( ARGV_LIST[i], "GENRANGE_AV" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_AV[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_AV[1] ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENTAU_AV" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENEXPTAU_AV );   }
    if ( strcmp( ARGV_LIST[i], "GENEXPTAU_AV" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENEXPTAU_AV );   }

    if ( strcmp( ARGV_LIST[i], "GENSIG_AV" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENGAUSIG_AV );   }
    if ( strcmp( ARGV_LIST[i], "GENGAUSIG_AV" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENGAUSIG_AV );   }

    if ( strcmp( ARGV_LIST[i], "GENGAUPEAK_AV" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENGAUPEAK_AV );   }

    if ( strcmp( ARGV_LIST[i], "GENRATIO_AV0" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRATIO_AV0 );   }

    sscanf_GENGAUSS(&i, "SALT2c",     &INPUTS.GENGAUSS_SALT2c     ) ;
    sscanf_GENGAUSS(&i, "SALT2x1",    &INPUTS.GENGAUSS_SALT2x1    ) ;
    sscanf_GENGAUSS(&i, "SALT2ALPHA", &INPUTS.GENGAUSS_SALT2ALPHA );
    sscanf_GENGAUSS(&i, "SALT2BETA",  &INPUTS.GENGAUSS_SALT2BETA  );

    if ( strcmp ( ARGV_LIST[i], "BIASCOR_SALT2GAMMA_GRID" ) == 0 )  { 
      i++; sscanf(ARGV_LIST[i] , "%le", &INPUTS.BIASCOR_SALT2GAMMA_GRID[0] );
      i++; sscanf(ARGV_LIST[i] , "%le", &INPUTS.BIASCOR_SALT2GAMMA_GRID[1] );
    }
      

    if ( strcmp ( ARGV_LIST[i], "SALT2BETA_cPOLY" ) == 0 ) {
       i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.SALT2BETA_cPOLY[0] );  
       i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.SALT2BETA_cPOLY[1] );  
       i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.SALT2BETA_cPOLY[2] );  
    }

    // May 2013: check option to read alpha,beta from SALT2mu-fitres file
    if ( strcmp( ARGV_LIST[i], "SALT2mu_FILE" ) == 0 )  
      { i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.SALT2mu_FILE );   }
    
    // legacy variables ...
    if ( strcmp( ARGV_LIST[i], "GENALPHA_SALT2" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENALPHA_SALT2 ); 
      INPUTS.GENGAUSS_SALT2ALPHA.PEAK = INPUTS.GENALPHA_SALT2; 
    }
    if ( strcmp( ARGV_LIST[i], "GENBETA_SALT2" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENBETA_SALT2 ); 
      INPUTS.GENGAUSS_SALT2BETA.PEAK = INPUTS.GENBETA_SALT2; 
    }
    if ( strcmp( ARGV_LIST[i], "LEGACY_colorXTMW_SALT2" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.LEGACY_colorXTMW_SALT2 ); 
    }

    // --------------

    if ( strcmp( ARGV_LIST[i], "GENAV_WV07" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.WV07_GENAV_FLAG ); }
    if ( strcmp( ARGV_LIST[i], "WV07_GENAV_FLAG" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.WV07_GENAV_FLAG );  }
    if ( strcmp( ARGV_LIST[i], "WV07_REWGT_EXPAV" ) == 0 )  
      { i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.WV07_REWGT_EXPAV ); }

    if ( strcmp( ARGV_LIST[i], "RVMW" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.RV_MWCOLORLAW ); 
    }
    if ( strcmp( ARGV_LIST[i], "RV_MWCOLORLAW" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.RV_MWCOLORLAW ); 
    }

    if ( strcmp( ARGV_LIST[i], "OPT_MWCOLORLAW" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.OPT_MWCOLORLAW ); 
    }
    if ( strcmp( ARGV_LIST[i], "OPT_MWEBV" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &opt_tmp );
      INPUTS.OPT_MWEBV = abs(opt_tmp) ;
      if ( opt_tmp < 0 )
	{ INPUTS.APPLYFLAG_MWEBV=1; } // correct FLUXCAL
      else
	{ INPUTS.APPLYFLAG_MWEBV=0; } // turn off correction
    }


    if ( strcmp( ARGV_LIST[i], "GENSIGMA_MWEBV_RATIO" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.MWEBV_SIGRATIO ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENSIGMA_MWEBV" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.MWEBV_SIG ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENSHIFT_MWEBV" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.MWEBV_SHIFT ); 
    }
    if ( strcmp( ARGV_LIST[i], "GENSCALE_MWEBV" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.MWEBV_SCALE ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENRANGE_MWEBV" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_MWEBV[0] ); 
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENRANGE_MWEBV[1] ); 
    }

    if ( strcmp( ARGV_LIST[i], "EXTINC_HOSTGAL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.GENSNXT ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENMAG_OFF_AB" ) == 0 ||   // legacy key
	 strcmp( ARGV_LIST[i], "GENMAG_OFF_ZP" ) == 0 ) {
      for ( ifilt = 0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
	i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.TMPOFF_ZP[ifilt] ); 
      }
    }

    if ( strcmp( ARGV_LIST[i], "GENMAG_OFF_GLOBAL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENMAG_OFF_GLOBAL ); 
    }
    if ( strcmp( ARGV_LIST[i], "GENMAG_OFF_NON1A" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENMAG_OFF_NON1A ); 
    }

    if ( strcmp( ARGV_LIST[i], "GENMAG_OFF_MODEL" ) == 0 ) {
      for ( ifilt = 0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
	i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.TMPOFF_MODEL[ifilt] ); 
      }
    }

    if ( strcmp( ARGV_LIST[i], "SIGMACLIP_MAGSMEAR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.SIGMACLIP_MAGSMEAR[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.SIGMACLIP_MAGSMEAR[1] );
    }

    if ( strcmp( ARGV_LIST[i], "GENMAG_SMEAR_FILTER" ) == 0 ) {

      i++ ; sscanf(ARGV_LIST[i] , "%s", ctmp );
      i++ ; sscanf(ARGV_LIST[i] , "%s", ctmp2 );

      if ( strcmp(ctmp2,"INTERP") == 0 ) 
	{ ftmp = -1.0 ; } // any negative number -> interpolate smear
      else
	{ sscanf(ctmp2, "%f", &ftmp); }

      // update for any value, even zero to allow override.

      N = INPUTS.NFILT_SMEAR ;
      INPUTS.NFILT_SMEAR += 
	PARSE_FILTLIST(ctmp, &INPUTS.IFILT_SMEAR[N+1] );
      for ( j=N+1; j <= INPUTS.NFILT_SMEAR; j++ ) {
	ifilt = INPUTS.IFILT_SMEAR[j];
	INPUTS.GENMAG_SMEAR_FILTER[ifilt] = ftmp ;
      }

    } // GENMAG_SMEAR_FILTER


    if ( strcmp( ARGV_LIST[i], "GENSMEAR_RANGaussFIX" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i],"%le",&INPUTS.GENSMEAR_RANGauss_FIX );  }
    if ( strcmp( ARGV_LIST[i], "GENSMEAR_RANGAUSSFIX" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i],"%le",&INPUTS.GENSMEAR_RANGauss_FIX );  }

    if ( strcmp( ARGV_LIST[i], "GENSMEAR_RANFlatFIX" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i],"%le",&INPUTS.GENSMEAR_RANFlat_FIX );  }
    if ( strcmp( ARGV_LIST[i], "GENSMEAR_RANFLATFIX" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i],"%le",&INPUTS.GENSMEAR_RANFlat_FIX );  }

    if ( strcmp( ARGV_LIST[i], "GENMAG_SMEAR_SCALE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%f", &INPUTS.GENMAG_SMEAR_SCALE);
    }

    if ( strcmp( ARGV_LIST[i], "GENMAG_SMEAR_MODELNAME" ) == 0 ) {
      modelName = INPUTS.GENMAG_SMEAR_MODELNAME ;
      i++ ; sscanf(ARGV_LIST[i] , "%s", modelName);
      if ( strcmp(modelName,"G10FUDGE") == 0  )  { // read sigma_coh
	i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENMAG_SMEAR_USRFUN[0]);
      }
      if ( strstr(modelName,":")!= NULL)  { parse_GENMAG_SMEAR_MODELNAME();}
    }

    if ( strstr( ARGV_LIST[i], "SIGMA_COH" ) != NULL ) {
      parse_SIGCOH_SALT2(ARGV_LIST[i],ARGV_LIST[i+1], 
			 &GENSMEAR_SALT2.SIGCOH_LAM );
      INPUTS.GENMAG_SMEAR_USRFUN[0] = -8.0;  // flag that SIGCOH_LAM is set
      i++ ;
    }

    // updated keys for strong & weak lensing (July 2019)
    if ( strcmp( ARGV_LIST[i], "STRONGLENS_FILE" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.STRONGLENS_FILE ); }

    if ( strcmp( ARGV_LIST[i], "WEAKLENS_PROBMAP_FILE" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.WEAKLENS_PROBMAP_FILE ); }
    if ( strcmp( ARGV_LIST[i], "WEAKLENS_DMUSCALE" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.WEAKLENS_DMUSCALE ); }
    if ( strcmp( ARGV_LIST[i], "WEAKLENS_DSIGMADZ" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.WEAKLENS_DSIGMADZ );  }

    // legacy keys for weak lensing 
    if ( strcmp( ARGV_LIST[i], "LENSING_PROBMAP_FILE" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.WEAKLENS_PROBMAP_FILE ); }
    if ( strcmp( ARGV_LIST[i], "LENSING_DMUSCALE" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.WEAKLENS_DMUSCALE ); }
    if ( strcmp( ARGV_LIST[i], "LENSING_DSIGMADZ" ) == 0 ) 
      { i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.WEAKLENS_DSIGMADZ );  }

    // - - - - - - 
    // check override params for genSmear models (Mar 22 2014)
    int  NVAL;   char parName[60], ckey[60];   double tmpList[100];
    if ( strcmp( ARGV_LIST[i], "GENMAG_SMEARPAR_OVERRIDE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", ckey);
      NVAL = nval_genSmear_override(ckey, parName); // return NVAL,parName
      for(j=0; j < NVAL; j++ ) 
	{ i++ ; sscanf(ARGV_LIST[i],"%le", &tmpList[j] ); }
      
      store_genSmear_override(parName,NVAL,tmpList);
    }
    // - - - - - - - - 
   
    if ( strcmp( ARGV_LIST[i], "GENMODEL_ERRSCALE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENMODEL_ERRSCALE );
    }
    if ( strcmp( ARGV_LIST[i], "GENMAG_SMEAR" ) == 0 ) {      
      i++ ; sscanf(ARGV_LIST[i] , "%s", ctmp) ;
      split2floats(ctmp, comma, tmpSmear);

    
      // Sep 2014: if NON1a magSmear is already set, 
      //           then do not use GENMAG_SMEAR.
      if ( INPUTS.NON1ASED.MAGSMEAR[1][0] <= 0.0 ) {
	INPUTS.GENMAG_SMEAR[0] = tmpSmear[0]; 
	INPUTS.GENMAG_SMEAR[1] = tmpSmear[1]; 
      }
	
    }

    if ( strcmp( ARGV_LIST[i], "GENMAG_SMEAR_USRFUN" ) == 0 ) {

      INPUTS.NPAR_GENSMEAR_USRFUN = 8 ;
      for ( j=0; j < INPUTS.NPAR_GENSMEAR_USRFUN; j++ ) {
	i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.GENMAG_SMEAR_USRFUN[j] );
	sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"USRFUN");
      }
      
      // turn off if first element is negative
      if ( INPUTS.GENMAG_SMEAR_USRFUN[0] < -0.001 ) {  
	sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"NONE"); 
	INPUTS.NPAR_GENSMEAR_USRFUN = 0 ;
      }

    }

    if ( strcmp( ARGV_LIST[i], "GENMODEL_ERRSCALE_OPT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENMODEL_ERRSCALE_OPT );
    }
    if ( strcmp( ARGV_LIST[i], "GENMODEL_ERRSCALE_CORRELATION" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", 
		   &INPUTS.GENMODEL_ERRSCALE_CORRELATION );
    }


    if ( strcmp( ARGV_LIST[i], "GENRANGE_TYPE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENRANGE_TYPE[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENRANGE_TYPE[1] );
    }

    if ( strcmp( ARGV_LIST[i], "GENRANGE_CID" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENRANGE_CID[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENRANGE_CID[1] );
    }

    if ( strcmp( ARGV_LIST[i], "GENRANGE_FILT" ) == 0 ) {
      sprintf(c1err,"GENRANGE_FILT: keyword no longer supported.");
      sprintf(c2err,"Must use GENFILTERS: ");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( strcmp( ARGV_LIST[i], "GENFILTERS" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.GENFILTERS );
      INPUTS.NFILTDEF_OBS =
	PARSE_FILTLIST( INPUTS.GENFILTERS, INPUTS.IFILTMAP_OBS  );
    }

    if ( strcmp( ARGV_LIST[i], "GENRANGE_DMPEVENT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENRANGE_DMPEVENT[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENRANGE_DMPEVENT[1] );
    }
    if ( strcmp( ARGV_LIST[i], "GENRANGE_DMPTREST" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_DMPTREST[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.GENRANGE_DMPTREST[1] );
    }

    if ( strcmp( ARGV_LIST[i], "SMEARFLAG_FLUX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SMEARFLAG_FLUX );
    }
    if ( strcmp( ARGV_LIST[i], "SMEARFLAG_ZEROPT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SMEARFLAG_ZEROPT );
    }
    if ( strcmp( ARGV_LIST[i], "MAGMONITOR_SNR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.MAGMONITOR_SNR );
    }

    if ( strcmp( ARGV_LIST[i], "KCORFLAG_COLOR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.KCORFLAG_COLOR );
    }

    if ( strcmp( ARGV_LIST[i], "EXPOSURE_TIME_MSKOPT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.EXPOSURE_TIME_MSKOPT );
    }

    if ( strcmp( ARGV_LIST[i], "KCOR_FILE" ) == 0  || 
	 strcmp( ARGV_LIST[i], "KCOR"      ) == 0  ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.KCOR_FILE  );
    }

    if ( strcmp( ARGV_LIST[i], "OMEGA_MATTER" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.OMEGA_MATTER  );
    }
    if ( strcmp( ARGV_LIST[i], "OMEGA_LAMBDA" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.OMEGA_LAMBDA  );
    }
    if ( strcmp( ARGV_LIST[i], "W0_LAMBDA" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.W0_LAMBDA  );
      
    }
    if ( strcmp( ARGV_LIST[i], "H0" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.H0  );
    }

    // --- search eff
    if ( strcmp( ARGV_LIST[i], "MAGSHIFT_SPECEFF" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS_SEARCHEFF.MAGSHIFT_SPECEFF );
    }

    if ( strcmp( ARGV_LIST[i], "SEARCHEFF_PIPELINE_LOGIC_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", 
		   INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE );
    } 

    if ( strcmp( ARGV_LIST[i], "SEARCHEFF_PIPELINE_FILE"     ) == 0 ||
	 strcmp( ARGV_LIST[i], "SEARCHEFF_PIPELINE_EFF_FILE" ) == 0 ) {	
      i++ ; sscanf(ARGV_LIST[i] , "%s", 
		   INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE );
    } 

    if ( strcmp( ARGV_LIST[i], "SEARCHEFF_SPEC_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i],"%s",INPUTS_SEARCHEFF.USER_SPEC_FILE );     
    } 
    if ( strcmp( ARGV_LIST[i], "SEARCHEFF_SPEC_SCALE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i],"%le",&INPUTS_SEARCHEFF.USER_SPECEFF_SCALE);
    } 

    if ( strcmp( ARGV_LIST[i], "SEARCHEFF_zHOST_FILE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i],"%s", INPUTS_SEARCHEFF.USER_zHOST_FILE ); 
    } 

    if ( strcmp( ARGV_LIST[i], "APPLY_SEARCHEFF_OPT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.APPLY_SEARCHEFF_OPT );
    }
    if ( strcmp( ARGV_LIST[i], "MINOBS_SEARCH" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS_SEARCHEFF.MINOBS );
    } 

    // -- CUT WINDOWS
    if ( strcmp( ARGV_LIST[i], "APPLY_CUTWIN_OPT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.APPLY_CUTWIN_OPT );
    }

    if ( strcmp( ARGV_LIST[i], "EPCUTWIN_LAMREST" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.EPCUTWIN_LAMREST[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.EPCUTWIN_LAMREST[1] );
    }

    if ( strcmp( ARGV_LIST[i], "EPCUTWIN_SNRMIN" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.EPCUTWIN_SNRMIN[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.EPCUTWIN_SNRMIN[1] );
    }

    if ( strcmp( ARGV_LIST[i], "CUTWIN_TRESTMAX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_TRESTMAX[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_TRESTMAX[1] );
    }


    if ( strcmp( ARGV_LIST[i], "CUTWIN_REDSHIFT" )      == 0 ||
	 strcmp( ARGV_LIST[i], "CUTWIN_REDSHIFT_TRUE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_REDSHIFT_TRUE[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_REDSHIFT_TRUE[1] );
    }
    if ( strcmp( ARGV_LIST[i], "CUTWIN_REDSHIFT_FINAL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_REDSHIFT_FINAL[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_REDSHIFT_FINAL[1] );
    }
    if ( strcmp( ARGV_LIST[i], "CUTWIN_HOST_PHOTOZ" ) == 0 ||
	 strcmp( ARGV_LIST[i], "CUTWIN_HOST_ZPHOT"  ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_HOST_ZPHOT[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_HOST_ZPHOT[1] );
    }


    if ( strcmp( ARGV_LIST[i], "CUTWIN_TRESTMIN" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_TRESTMIN[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_TRESTMIN[1] );
    }

    if ( strcmp( ARGV_LIST[i], "CUTWIN_TGAPMAX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_TGAPMAX[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_TGAPMAX[1] );
    }
    if ( strcmp( ARGV_LIST[i], "CUTWIN_T0GAPMAX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_T0GAPMAX[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_T0GAPMAX[1] );
    }

    if ( strcmp( ARGV_LIST[i], "CUTWIN_SNRMAX" ) == 0 ) {

      if ( INPUTS.OVERRIDE_CUTWIN_SNRMAX == 0 ) 
	{ INPUTS.NCUTWIN_SNRMAX = 1; } // reset
      else
	{ INPUTS.NCUTWIN_SNRMAX++ ; }  //  keep adding

      INPUTS.OVERRIDE_CUTWIN_SNRMAX = 1 ;
      N = INPUTS.NCUTWIN_SNRMAX ;
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_SNRMAX[N][0]       );
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.CUTWIN_SNRMAX_FILTERS[N]   );
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.CUTWIN_SNRMAX_LOGIC[N]     );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_SNRMAX_TREST[N][0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_SNRMAX_TREST[N][1] );
    }
    if ( strcmp( ARGV_LIST[i], "CUTWIN_NEPOCH" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_NEPOCH[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_NEPOCH[1] );
    }

    if ( strcmp( ARGV_LIST[i], "CUTWIN_NOBSDIF" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.CUTWIN_NOBSDIF[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.CUTWIN_NOBSDIF[1] );
    }
    if ( strcmp( ARGV_LIST[i], "CUTWIN_MJDDIF" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_MJDDIF[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_MJDDIF[1] );
    }

    if ( strcmp( ARGV_LIST[i], "CUTWIN_MWEBV" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_MWEBV[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_MWEBV[1] );
    }

    if ( strcmp( ARGV_LIST[i], "CUTWIN_PEAKMAG" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_PEAKMAG[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_PEAKMAG[1] );
    }
    if ( strcmp( ARGV_LIST[i], "CUTWIN_PEAKMAG_ALL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_PEAKMAG_ALL[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_PEAKMAG_ALL[1] );
    }

    if ( strcmp( ARGV_LIST[i], "CUTWIN_EPOCHS_SNRMIN" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_EPOCHS_SNRMIN );
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.CUTWIN_EPOCHS_TRANGE[1] );
      i++ ; sscanf(ARGV_LIST[i] , "%s",  INPUTS.CUTWIN_EPOCHS_FILTERS );
    }

    if ( strcmp( ARGV_LIST[i], "EFFERR_STOPGEN" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.EFFERR_STOPGEN );
    }

    if ( strcmp( ARGV_LIST[i], "SPECTROGRAPH_OPTMASK" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK );
    }

    if ( strcmp( ARGV_LIST[i], "SPECTROGRAPH_SCALE_TEXPOSE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", 
		   &INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE );
    }

    if ( strcmp( ARGV_LIST[i], "TAKE_SPECTRUM_DUMPCID" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.TAKE_SPECTRUM_DUMPCID );
    }

    if ( strcmp( ARGV_LIST[i], "TAKE_SPECTRUM_HOSTFRAC" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%f", &INPUTS.TAKE_SPECTRUM_HOSTFRAC );
    }


    if ( strstr(ARGV_LIST[i],"SIMGEN_DUMP") != NULL ) 
      { parse_input_SIMGEN_DUMP(fpNull, &i, ARGV_LIST[i] ); }


    // FUDGES

    if ( strcmp( ARGV_LIST[i], "FUDGE_SNRMAX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.STRING_FUDGE_SNRMAX );
      INPUTS.OPT_FUDGE_SNRMAX = 1; 
    }
    if ( strcmp( ARGV_LIST[i], "FUDGE2_SNRMAX" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.STRING_FUDGE_SNRMAX );
      INPUTS.OPT_FUDGE_SNRMAX = 2 ; 
    }

    if ( strcmp( ARGV_LIST[i], "GENPERFECT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.GENPERFECT );
    }

    // - - - - -
    if ( strcmp(ARGV_LIST[i], "LCLIB_DEBUG_DAYSCALE") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &LCLIB_DEBUG.DAYSCALE );
    }
    if ( strcmp(ARGV_LIST[i], "LCLIB_DEBUG_TOBS_OFFSET") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &LCLIB_DEBUG.TOBS_OFFSET_RANGE[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &LCLIB_DEBUG.TOBS_OFFSET_RANGE[1] );
    }
    if ( strcmp(ARGV_LIST[i], "LCLIB_DEBUG_ZERO_TEMPLATE_FLUX") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &LCLIB_DEBUG.ZERO_TEMPLATE_FLUX );
    }
    if ( strcmp(ARGV_LIST[i], "LCLIB_DEBUG_FORCE_NREPEAT") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &LCLIB_DEBUG.FORCE_NREPEAT );
    }
    if ( strcmp(ARGV_LIST[i], "LCLIB_CUTWIN") == 0 ) {
      N = LCLIB_CUTS.NCUTWIN ;
      i++ ; sscanf(ARGV_LIST[i] , "%s",  LCLIB_CUTS.PARNAME[N]   );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &LCLIB_CUTS.CUTWIN[N][0] );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &LCLIB_CUTS.CUTWIN[N][1] );
      LCLIB_CUTS.NCUTWIN++ ;
    }


#ifdef SNGRIDGEN
    // GRID input
    if ( strcmp( ARGV_LIST[i], "GRID_FORMAT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", GRIDGEN_INPUTS.FORMAT );
    }
    if ( strcmp( ARGV_LIST[i], "NGRID_LOGZ" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &N );
      GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_LOGZ] = N;
    }
    if ( strcmp( ARGV_LIST[i], "NGRID_LUMIPAR" ) == 0 ) {  // legacy key
      i++ ; sscanf(ARGV_LIST[i] , "%d", &N );
      GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_SHAPEPAR] = N;
    }
    if ( strcmp( ARGV_LIST[i], "NGRID_SHAPEPAR" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &N );
      GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_SHAPEPAR] = N;
    }
    if ( strcmp( ARGV_LIST[i], "NGRID_TREST" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &N );
      GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_TREST] = N;
    }
#endif

    // SIMSED args

    if ( strcmp( ARGV_LIST[i], "SIMSED_USE_BINARY" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.USE_BINARY_SIMSED );
    }

    if ( strcmp( ARGV_LIST[i], "SIMSED_PATH_BINARY" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.PATH_BINARY_SIMSED );
    }

    if ( strcmp( ARGV_LIST[i], "SIMSED_GRIDONLY" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", ctmp );
      if ( strcmp(ctmp,"SEQUENTIAL") == 0  ) {
	INPUTS.NPAR_SIMSED = 0 ;  // no interp pars => treat as baggage
	INPUTS.OPTMASK_SIMSED = OPTMASK_SIMSED_GRIDONLY ;
	// xxx mark delete  INPUTS.GENFLAG_SIMSED[0] = OPTMASK_SIMSED_GRIDONLY ;
      }
      else {
	sprintf(c1err,"Invalid command-line option");
	sprintf(c2err,"'SIMSED_GRIDONLY: %s'", ctmp);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }

    for ( ipar=0; ipar < INPUTS.NPAR_SIMSED; ipar++ ) {
      opt_tmp = INPUTS.GENFLAG_SIMSED[ipar];
      if ( (opt_tmp & OPTMASK_SIMSED_param ) > 0 ) 
	{ continue; }  // allow only interp params 
      sprintf(parname,"%s", INPUTS.PARNAME_SIMSED[ipar] );
      sscanf_GENGAUSS(&i, parname, &INPUTS.GENGAUSS_SIMSED[ipar] );
    } // ipar loop



    // -------------------
  INCREMENT_COUNTER:

    if ( i > ilast ) {
      for ( iuse = ilast; iuse <= i; iuse++ ) { USE_ARGV_LIST[iuse] = 1; }
    }

    i++ ;  ilast=i;

  }
  return ;

} // end of sim_input_override


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
   + remove !MODEL_FIXMAG requirement for prep_genmag_offsets().

  Aug 16 2017: abort of GENRANGE_REDSHIFT[1] > ZMAX_SNANA

  Sep 25 2017: move CUTWIN stuff into prep_user_CUTWIN()
  Jan 23 2018: check WRFLAG_COMPACT
  Sep 05 2018: set INPUTS.NON1ASED.PATH
  Oct 04 2018: set INPUTS.USE_HOSTLIB_GENZPHOT
  Feb 15 2019: turn off INPUTS.GENMAG_SMEAR_MODELNAME for SIMSED model.

  *******************/

  int ifilt, ifilt_obs, IBLANK,i, j, indx, lentmp, NGEN, OPT, NTMP  ;
  int DOCHECK_FORMAT_MASK=1, ISRATE_LCLIB, INDEX_RATEMODEL;

  float tmp_F, XGEN_F ; 
  double mu, z;
  char  *PTR_GENMODEL, vtmp[40]   ;
  char *input_file = INPUTS.INPUT_FILE_LIST[0] ;
  char fnam[] = "prep_user_input" ;

  // ---------- BEGIN ----------

  // Feb 2015: replace ENV names in inputs
  ENVreplace(INPUTS.KCOR_FILE,fnam,1);  
  ENVreplace(INPUTS.SIMLIB_FILE,fnam,1);
  ENVreplace(INPUTS.HOSTLIB_FILE,fnam,1);
  ENVreplace(INPUTS.HOSTLIB_WGTMAP_FILE,fnam,1);
  ENVreplace(INPUTS.HOSTLIB_ZPHOTEFF_FILE,fnam,1);
  ENVreplace(INPUTS.FLUXERRMODEL_FILE,fnam,1 );
  ENVreplace(INPUTS.HOSTNOISE_FILE,fnam,1 );
  ENVreplace(INPUTS.WRONGHOST_FILE,fnam,1 );
  ENVreplace(INPUTS.SALT2mu_FILE,fnam,1 );
  ENVreplace(INPUTS.PATH_BINARY_SIMSED,fnam,1);
  ENVreplace(INPUTS.NON1ASED.PATH,fnam,1);
  ENVreplace(INPUTS.GENPAR_SELECT_FILE,fnam,1);
  ENVreplace(INPUTS.NON1AGRID_FILE,fnam,1);
  ENVreplace(INPUTS.NONLINEARITY_FILE,fnam,1);
  ENVreplace(INPUTS.WEAKLENS_PROBMAP_FILE,fnam,1);
  ENVreplace(INPUTS.STRONGLENS_FILE,fnam,1);
  ENVreplace(INPUTS.LCLIB_FILE,fnam,1);
  ENVreplace(INPUTS.MODELPATH,fnam,1);

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
      if ( strncmp(INPUTS.MODELNAME,PTR_GENMODEL,lentmp) == 0 ) 
	{ INDEX_GENMODEL = indx ; }
    }
  }

  if ( INDEX_GENMODEL <= 0 ) {
    sprintf(c1err,"%s is not a valid genmag-model", INPUTS.MODELNAME);
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
    PTR_GENMODEL = GENMODEL_NAME[INDEX_GENMODEL][0] ; // generic model name
    sprintf(INPUTS.MODELPATH,"%s/models/%s/%s", 
	    PATH_SNDATA_ROOT, PTR_GENMODEL, INPUTS.MODELNAME );
  }

  // check for optional over-ride of model path;
  // getenv(PRIVATE_MODELPATH_NAME is equivalent of 
  // $SNDATA_ROOT/models/[SALT2,mlcs2k2...]
  if ( getenv(PRIVATE_MODELPATH_NAME) != NULL ) {
    sprintf(INPUTS.MODELPATH,"%s/%s", 
	    getenv(PRIVATE_MODELPATH_NAME), INPUTS.MODELNAME );
  }


  // -------------------------------
  // set INPUTS.GENSOURCE

  sprintf(vtmp, "%s", INPUTS.GENSOURCE);
  if ( strcmp(vtmp,"RANDOM") == 0 ) 
    { GENLC.IFLAG_GENSOURCE = IFLAG_GENRANDOM ; }
  else if ( strcmp(vtmp,"GRID") == 0 ) 
    { GENLC.IFLAG_GENSOURCE = IFLAG_GENGRID  ; }

  // set OPT_SNXT
  if ( strcmp(INPUTS.GENSNXT,"CCM89") == 0 ) 
    { OPT_SNXT = OPT_SNXT_CCM89; }
  else if ( strcmp(INPUTS.GENSNXT,"SJPAR") == 0 ) 
    { OPT_SNXT = OPT_SNXT_SJPAR; }
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
    sprintf(GENLC.COLORPAR_NAME,    "AV"    );
    sprintf(GENLC.COLORPAR2_NAME,   "RV"    );
    GENLC.ptr_SHAPEPAR = &GENLC.DELTA ;
    copy_GENGAUSS_ASYM( &INPUTS.GENGAUSS_DELTA, &INPUTS.GENGAUSS_SHAPEPAR );
  }

  else if ( INDEX_GENMODEL == MODEL_SNOOPY  ) {
    
    GENFRAME_OPT    = GENFRAME_REST; 
    sprintf(GENLC.DISTANCE_NAME,    "DLMAG" );
    sprintf(GENLC.SHAPEPAR_NAME,    "DM15"  );
    sprintf(GENLC.SHAPEPAR_GENNAME, "DM15"  );
    sprintf(GENLC.COLORPAR_NAME,    "AV"    );
    sprintf(GENLC.COLORPAR2_NAME,   "RV"    );
    GENLC.ptr_SHAPEPAR = &GENLC.DM15 ;
    copy_GENGAUSS_ASYM( &INPUTS.GENGAUSS_DM15, &INPUTS.GENGAUSS_SHAPEPAR );
  }

  else if ( INDEX_GENMODEL == MODEL_S11DM15  ) {

    GENFRAME_OPT    = GENFRAME_OBS; 
    sprintf(GENLC.DISTANCE_NAME,     "DLMAG" );
    sprintf(GENLC.SHAPEPAR_NAME,     "DM15"  );
    sprintf(GENLC.SHAPEPAR_GENNAME,  "DM15"  );
    sprintf(GENLC.COLORPAR_NAME,     "AV"    );
    sprintf(GENLC.COLORPAR2_NAME,    "RV"    );
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

    read_SALT2mu_AlphaBeta(INPUTS.SALT2mu_FILE);
    copy_GENGAUSS_ASYM( &INPUTS.GENGAUSS_SALT2x1, &INPUTS.GENGAUSS_SHAPEPAR );
  }
  else if ( INDEX_GENMODEL == MODEL_SIMSED ) {
   
    GENFRAME_OPT   = GENFRAME_OBS;
    INPUTS.GENGAUSS_SHAPEPAR.RANGE[0] = -1.0 ; // anything non-zero
    INPUTS.GENGAUSS_SHAPEPAR.RANGE[1] = +1.0 ;
    prep_user_SIMSED(); // Feb 26 2018
    sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"NONE"); // no SNIa scatter model
  }

  else if ( INDEX_GENMODEL == MODEL_BYOSED ) {

    GENFRAME_OPT   = GENFRAME_OBS;

  }

  else  if ( INDEX_GENMODEL  == MODEL_NON1ASED ) {    
    GENFRAME_OPT    = GENFRAME_OBS ;
    sprintf(INPUTS.GENMAG_SMEAR_MODELNAME,"NONE"); // no SNIa scatter model
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
      // abort of redshift range not specified
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
  text_MWoption("EBV", OPT, INPUTS.STR_MWEBV ); // return STR
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
    if ( INPUTS.MJD_TEMPLATE_FILTER[ifilt] > 1.0 ) { INPUTS.USE_MJD_TEMPLATE = 1 ; }
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


  // check for host AV
  if ( INPUTS.WV07_REWGT_EXPAV > -1.0E-9 ) { INPUTS.WV07_GENAV_FLAG=1; }

  int DO_RV    = INPUTS.GENGAUSS_RV.PEAK > 1.0E-9 ;  
  int DO_AVTAU = INPUTS.GENEXPTAU_AV     > 1.0E-9 ;
  int DO_AVSIG = INPUTS.GENGAUSIG_AV     > 1.0E-9 ;
  int DO_AV    = INPUTS.GENRANGE_AV[1]   > 1.0E-9 ;
  INPUTS.DO_AV = (DO_AV && DO_RV && 
		  ( DO_AVTAU || DO_AVSIG || INPUTS.WV07_GENAV_FLAG) ) ;

  if ( INPUTS.DO_AV==0 && DO_AV ) {
    printf("\n PRE-ABORT DUMP: \n");
    printf("\t GENPEAK_RV   = %f \n", INPUTS.GENGAUSS_RV.PEAK );
    printf("\t GENEXPTAU_AV = %f \n", INPUTS.GENEXPTAU_AV );
    printf("\t GENGAUSIG_AV = %f \n", INPUTS.GENGAUSIG_AV );
    sprintf(c1err,"GENRANGE_AV = %.3f to %.3f", 
	    INPUTS.GENRANGE_AV[0], INPUTS.GENRANGE_AV[1] );
    sprintf(c2err,"But cannot generate AV>0; see above param-dump");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // --------------------------------------
  //----------- PRINT SUMMARY -------------
  // --------------------------------------

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

  printf("\t OMEGA_(MATTER,LAMBDA)= %5.3f, %5.3f,    W= %5.2f   H0=%5.1f \n"
	 ,INPUTS.OMEGA_MATTER
	 ,INPUTS.OMEGA_LAMBDA
	 ,INPUTS.W0_LAMBDA
	 ,INPUTS.H0 	 );


  printf("\t KCOR  file : %s \n", INPUTS.KCOR_FILE );


  printf("\t Observer Gen-FILTERS  :  %s ", INPUTS.GENFILTERS );

  printf(" \n" );

  printf("\t Random number seed: %d \n", INPUTS.ISEED );

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

  // check for required input;
  // Allow exceptions for GRID option

  if ( GENLC.IFLAG_GENSOURCE != IFLAG_GENGRID ) {
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
    INPUTS.WRITE_MASK    += WRITE_MASK_SIM_SNRMON;
  }

  if ( INPUTS.WRFLAG_MODELPAR ) { INPUTS.WRITE_MASK += WRITE_MASK_SIM_MODELPAR; }

  // abort if no valid format is given
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENRANDOM  && 
       DOCHECK_FORMAT_MASK  ) {
    if (WRFLAG_VBOSE==0 && WRFLAG_TEXT==0 && WRFLAG_FITS==0 ) {
      sprintf(c1err,"Invalid FORMAT_MASK=%d. Must specify ", 
	      INPUTS.FORMAT_MASK) ;
      sprintf(c2err,"TEXT (%d)  or FITS (%d)", 
	      WRMASK_TEXT, WRMASK_FITS);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  if ( INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_GALMAG ) {
    INPUTS.SMEARFLAG_HOSTGAL += SMEARMASK_HOSTGAL_PHOT; // internal logical flag
  }

  if ( IGNOREFILE(INPUTS.HOSTNOISE_FILE) == 0 ) {
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
  if ( IGNOREFILE(INPUTS.FLUXERRMODEL_FILE) == 0 ) {
    if ( IGNOREFILE(INPUTS.HOSTNOISE_FILE) == 0 ) {
      sprintf(c1err,"Cannot mix FLUXERRMODEL_FILE with HOSTNOISE_FILE");
      sprintf(c2err,"Pick one, not both.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( INPUTS.FUDGEOPT_FLUXERR  ) {
      sprintf(c1err,"Cannot mix FLUXERRMODEL_FILE with FUDGEOPT_FLUXERR");
      sprintf(c2err,"Pick one, not both.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  
  printf("\n");

  return ;

} // end of prep_user_input().

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
  double *COVMAT_1D, SIG0, SIG1, RHO, COV_OFFDIAG ;

  // allocate memory for COVMAT and CHOLESKY 
  COVMAT_1D                  = (double *) malloc ( NROW*NROW * sizeof(double) );
  INPUTS.CHOLESKY_SIMSED_COV = (double**) malloc ( NROW * sizeof(double*) );
  for(irow=0; irow < NROW; irow++ ) {
    MEMD = NROW * sizeof(double);
    INPUTS.CHOLESKY_SIMSED_COV[irow] = (double*) malloc ( MEMD );
  }


  // load COVMAT_1D
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
	{ COVMAT_1D[NMAT] = SIG0*SIG0 ; }    // diagonal
      else 
	{ COVMAT_1D[NMAT] = COV_OFFDIAG ; }  // off-diag


      RHO = COVMAT_1D[NMAT] / ( SIG0*SIG1 );
      printf("%8.4f ", RHO );

      NMAT++ ;
    }
    printf("\n"); fflush(stdout);
  }

  // -------------------------------------------
  // store Cholesky decomp
  gsl_matrix_view chk; 
  chk = gsl_matrix_view_array ( COVMAT_1D, NROW, NROW ); 
  gsl_linalg_cholesky_decomp ( &chk.matrix)  ;  
  
  for (irow=0; irow<NROW; irow++){
    for (irow1 = 0; irow1 < NROW ; irow1++){      
      if (irow <= irow1 ) { 
	INPUTS.CHOLESKY_SIMSED_COV[irow][irow1] = 
	  gsl_matrix_get(&chk.matrix,irow,irow1) ;
      }
      else {
	INPUTS.CHOLESKY_SIMSED_COV[irow][irow1] = 0.0 ; // Apr 7 2018
      }
    }    
  }
  
  printf("\n");
  //  debugexit(fnam);

  return ;
	
} // end prep_user_SIMSED

// *******************************************
void  read_SALT2mu_AlphaBeta(char *SALT2mu_FILE) {

  // May 2013
  // If SALT2mu_FILE is defined then read alpha0 and beta0
  // from SALT2mu fitres file. Store alpha,beta values in
  //
  //  INPUTS.GENGAUSS_SALT2ALPHA.PEAK
  //  INPUTS.GENGAUSS_SALT2BET.PEAK
  //
  // so that they are used for generation.
  // 
  // -------

  FILE *fp ;
  double a,b;
  char c_get[100], cdum[12];
  char fnam[] = "read_SALT2mu_AlphaBeta" ;
  int NRD;
  // ------------ BEGIN -----------
  
  if ( strlen(SALT2mu_FILE) == 0 ) { return ; }
  
  
  if ( (fp = fopen(SALT2mu_FILE, "rt"))==NULL ) {       
      sprintf ( c1err, "Cannot open SALT2mu file to read alpha,beta: " );
      sprintf ( c2err," '%s' ", SALT2mu_FILE );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  sprintf(BANNER, "Read SALT2 alpha, beta from SALT2mu file: ");
  print_banner(BANNER);
  printf("\t %s\n", SALT2mu_FILE);

  a = b = -999. ;  // init alpha,beta to unphysical values
  NRD = 0;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    NRD++ ; // for debug only
    if ( strcmp(c_get,"alpha0") == 0 ) {
      readchar(fp, cdum); // skip '=' sign
      readdouble(fp, 1, &a);
    }

    if ( strcmp(c_get,"beta0") == 0 ) {
      readchar(fp, cdum); // skip '=' sign
      readdouble(fp, 1, &b);
    }

    // when we hit the first SN key, we are done reading
    if ( strcmp(c_get,"SN:") == 0 ) { goto DONE_READING ; }

  }  // end while

 DONE_READING:

  if ( a < -99 || b < -99 ) {
    sprintf ( c1err, "Could not find alpha0 and/or beta0");
    sprintf ( c2err," alpha = %f  beta=%f ", a, b );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\t Will generate SALT2 with alpha=%6.3f  beta=%6.3f \n", a, b );
  INPUTS.GENGAUSS_SALT2ALPHA.PEAK = a ; 
  INPUTS.GENGAUSS_SALT2BETA.PEAK  = b ;


} // end of read_SALT2mu_AlphaBeta


// **************************************
void  prep_RANSYSTPAR(void) {

  // Created Jun 2017
  // Prepare optional systematic offsets using random numbers.
  // This allows user to specify one set of Gaussian sigmas
  // (using RANSYSTPAR_XXX params) and then changing RANSEED
  // results in random set of systematic offsets.
  //
  // Do NOT try to sync randoms here; burn randoms only if required.

  int   i, ifilt, ifilt_obs, NSET=0 ;
  int   NFILTDEF = INPUTS.NFILTDEF_OBS ;
  int   ILIST_RAN=1;
  float tmp, tmpSigma ;
  double gmin = -3.0, gmax=+3.0; // Gaussian clip params
  char cfilt[2];
  char fnam[] = "prep_RANSYSTPAR" ;

  // ---------- BEGIN -----------

  if ( INPUTS.RANSYSTPAR.USE == 0 ) { return ; }

  sprintf(BANNER,"%s: Prepare Random set of Systematic Errors", fnam );
  print_banner(BANNER);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //   Start with variations that are synched among sub-samples
  //   (e.g., among separate simulation for each survey)
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  printf("\t* First Sync-Syst Random : %f \n", FlatRan1(ILIST_RAN) );

  // Galactic extinction
  tmpSigma = INPUTS.RANSYSTPAR.FUDGESCALE_MWEBV - 1.0 ;
  if ( tmpSigma != 0.0 ) {   
    NSET=1; tmp = 1.0 + tmpSigma * GaussRanClip(ILIST_RAN,gmin,gmax);
    INPUTS.MWEBV_SCALE = tmp;
    printf("\t FUDGESCALE_MWEBV  = %.2f \n", tmp );
  }

  tmpSigma = INPUTS.RANSYSTPAR.FUDGESHIFT_MWRV ;
  if ( tmpSigma != 0.0 ) { 
    NSET=1; tmp = tmpSigma * GaussRanClip(ILIST_RAN,gmin,gmax);
    INPUTS.RV_MWCOLORLAW += tmp ;
    printf("\t RV_MWCOLORLAW  = %.3f \n", INPUTS.RV_MWCOLORLAW );
  }

  printf("\t* Last  Sync-Syst Random : %f "
	 "(should be same each survey)\n", FlatRan1(ILIST_RAN) );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //     Now the unsynched variations
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // need to unsync the randoms among surveys. Add 137*IDSURVEY
  // to ISEED and then re-init the randoms with new SEED.

  int IDUM      = GENLC.IDSURVEY ;
  int ISEED_OLD = INPUTS.ISEED ;
  int ISEED_NEW = ISEED_OLD + 137*IDUM ;
  INPUTS.ISEED  = ISEED_NEW ;
  init_simRandoms();
  printf("\t* ISEED = %d --> %d \n", ISEED_OLD, ISEED_NEW );

  printf("\t* First Unsync-Syst Random : %f "
	 "(should differ each survey)\n", FlatRan1(ILIST_RAN) );

  // start with fluxerr fudging; only the deviation from 1 varies.
  tmpSigma = INPUTS.RANSYSTPAR.FUDGESCALE_FLUXERR - 1.0 ;
  if ( tmpSigma != 0.0 ) {   
    NSET=1; tmp = 1.0 + tmpSigma * GaussRanClip(ILIST_RAN,gmin,gmax);
    INPUTS.FUDGESCALE_FLUXERR = tmp;
    printf("\t FUDGESCALE_FLUXERR(true&measured) = %.3f \n", tmp );
    for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) 
      { INPUTS.FUDGESCALE_FLUXERR_FILTER[ifilt] = tmp; }
  }

  tmpSigma = INPUTS.RANSYSTPAR.FUDGESCALE_FLUXERR2 - 1.0 ;
  if ( tmpSigma != 0.0 ) {   
    NSET=1; tmp = 1.0 + tmpSigma * GaussRanClip(ILIST_RAN,gmin,gmax);
    INPUTS.FUDGESCALE_FLUXERR2 = tmp;
    printf("\t FUDGESCALE_FLUXERR2(measured) = %.3f \n", tmp );
    for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) 
      { INPUTS.FUDGESCALE_FLUXERR2_FILTER[ifilt] = tmp; }
  }

  // ZP error
  for(ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt];
    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );
    tmpSigma = INPUTS.RANSYSTPAR.GENMAG_OFF_ZP[ifilt];
    if ( tmpSigma != 0.0 ) {
      NSET++ ;  tmp = tmpSigma * GaussRanClip(ILIST_RAN,gmin,gmax);
      INPUTS.TMPOFF_ZP[ifilt]         = tmp;
      INPUTS.GENMAG_OFF_ZP[ifilt_obs] = tmp ;
      printf("\t ZPerr(%s) = %7.4f \n", cfilt, tmp );
    }
  }

  // LAMshift error
  for(ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt];
    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );
    tmpSigma = INPUTS.RANSYSTPAR.FILTER_LAMSHIFT[ifilt];
    if ( tmpSigma != 0 ) {
      NSET++ ;  tmp = tmpSigma * GaussRanClip(ILIST_RAN,gmin,gmax);
      INPUTS.TMPOFF_LAMSHIFT[ifilt] = tmp;
      printf("\t LAMSHIFT(%s) = %6.2f A \n", cfilt, tmp );
    }
  }

  
  printf("   %d Systematic Errors have been set. \n\n", NSET);

  return ;

} // end prep_RANSYSTPAR

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

  int epoch, ifilt_obs;
  double MAGOFF, genmag8 ;
  char fnam[] = "genmag_offsets";

  // ------------- begin -------------

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {

    ifilt_obs = GENLC.IFILT_OBS[epoch] ;

    // Aug 2012: always init USE_EPOCH[epoch]=0 (fixes IDLOCK bug)     
    GENLC.USE_EPOCH[epoch] = 0 ; 
    
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }
    
    // mag offset is :
    MAGOFF = 
      + GENLC.GENMAG_OFF_GLOBAL              // user-defined global offset
      + INPUTS.GENMAG_OFF_MODEL[ifilt_obs]   // user-defined model offs
      - INPUTS.GENMAG_OFF_ZP[ifilt_obs]      // user-defined ZP offsets
      + GENLC.LENSDMU                        // lensing correction
    ;

    if ( INPUTS.OPT_DEVEL_BBC7D ) 
      { MAGOFF += GENLC.SALT2gammaDM ;   }


    // ------
    // apply mag-offset to each epoch-mag, unless mag is
    // set to flag value (UNDEFINED or ZEROFLUX)
    // Note that mag=99 --> zero flux is well defined.
    genmag8 = GENLC.genmag8_obs[epoch];
    
    if ( genmag8 < 50.0 && genmag8 > -2.0 ) 
      {  genmag8 = GENLC.genmag8_obs[epoch] + MAGOFF ; }

    if ( genmag8 > 600.0 ) 
      { genmag8 = MAG_UNDEFINED; } // avoid crazy values
   
    GENLC.genmag8_obs[epoch] = genmag8 ;        // load global struct

    // -----------------------------
    // keep peak mags separately (before smearing)
    if ( GENLC.ISPEAK[epoch]  )  
      {  GENLC.peakmag8_obs[ifilt_obs] = genmag8 ; }

    if ( GENLC.ISTEMPLATE[epoch] ) 
      { GENLC.genmag8_obs_template[ifilt_obs] = genmag8; }

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


  char fnam[] = "prep_simpath" ;

  char tmp_path[5*MXPATHLEN] ;

  int lenpath, lenprefix, lenfile ;

    // ---------- BEGIN ----------

  sprintf(tmp_path, "%s/%s", 
	  PATH_SNDATA_SIM, INPUTS.GENVERSION );

  // check string lengths to avoid memory overwrites
  lenpath    = strlen(tmp_path);
  lenprefix  = strlen(INPUTS.GENPREFIX);
  lenfile    = lenpath + lenprefix + 20 ; // allow file-name extensions

  if ( lenprefix >= MXVERLEN ) {
    sprintf(c1err,"GENPREFIX string len = %d exceeds array bound of %d",
	    lenprefix, MXVERLEN);
    sprintf(c2err,"See input GENPREFIX: %s", INPUTS.GENPREFIX);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( lenfile >= MXPATHLEN  ) {
    sprintf(c1err, "Estimated filename length= %d is too long", lenfile);
    sprintf(c2err, "LEN(path,prefix) = %d, %d : MXPATHLEN=%d",
	    lenpath, lenprefix, MXPATHLEN );
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
#define MSKPERFECT_ALL       31  // mask with all bits

  //  char fnam[] = "genperfect_override" ;

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



  NVAR++ ;  iptr = &INPUTS.HOSTLIB_USE ;
  sprintf(GENPERFECT.parnam[NVAR], "HOSTLIB_USE" );
  GENPERFECT.parval[NVAR][0] = (float)*iptr ;
  *iptr = 0 ;
  GENPERFECT.parval[NVAR][1] = (float)*iptr ;
  GENPERFECT.partype[NVAR]   = 1 ;


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

  if ( INPUTS.GENRANGE_REDSHIFT[1] > 1.0E-8 )     
    { init_DNDZ_Rate(); } // Rate vs. redshift
  else
    { init_DNDB_Rate(); } // Rate vs. Gal coords l & b
  
  // --------------------------------------
  // now print info to stdout

  for ( i=0; i< NLINE_RATE_INFO; i++ ) 
    { printf("%s\n", LINE_RATE_INFO[i] ); }

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

  *****/

  double 
    H0, OM, OL, W0, Z0, Z1, Z_AVG, ZVtmp[2]
    ,ZVint[2], ZVOL, VT
    ,dOmega, dOmega_user, dphi, dth, sin0, sin1
    ,delMJD, Tyear, Tcomoving, SNsum, PEC1Asum, TOTsum
    ,ztmp, rtmp, rtmp1, rtmp2, FRAC_PEC1A
    ;

  int i, iz ;
  char ch[8], cnum[20];
  char fnam[] = "init_DNDZ_Rate" ;

  // ----------- BEGIN ------------

  H0 = INPUTS.H0 ;
  OM = INPUTS.OMEGA_MATTER ;
  OL = INPUTS.OMEGA_LAMBDA ;
  W0 = INPUTS.W0_LAMBDA  ;
  Z0 = INPUTS.GENRANGE_REDSHIFT[0] ;
  Z1 = INPUTS.GENRANGE_REDSHIFT[1] ;

  ZVint[0] = dVdz_integral ( H0, OM, OL, W0, Z0, 0 );
  ZVint[1] = dVdz_integral ( H0, OM, OL, W0, Z1, 0 );

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

  // compute average <z> wgted by volume (last arge = 1)
  ZVtmp[0] = dVdz_integral ( H0, OM, OL, W0, Z0, 1 );
  ZVtmp[1] = dVdz_integral ( H0, OM, OL, W0, Z1, 1 );

  if ( Z0 == Z1 ) 
    { Z_AVG = Z1 ; }
  else
    { Z_AVG    = (ZVtmp[1]-ZVtmp[0]) / (ZVint[1] - ZVint[0]); }


  // get time in years
  delMJD    = INPUTS.GENRANGE_PEAKMJD[1] - INPUTS.GENRANGE_PEAKMJD[0] ;
  Tyear     = delMJD / 365.0 ;
  Tcomoving = Tyear / ( 1.0 + Z_AVG ) ;
  VT        = ZVOL * Tcomoving ;

  sprintf(ch, "h%d" , (int)H0 );

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
	 "\t Survey Volume  = %8.4e  sr*(MPc/%s)^3 ", ZVOL, ch );

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t Survey Time    = %7.4f  years/season ", Tyear );

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t Co-moving Time = %7.4f  years/season  [ T/(1+<z>) ] ", 
	       Tcomoving );

  i++; sprintf(LINE_RATE_INFO[i],
	 "\t Co-moving V*T  = %8.4e  sr*(MPc/%s)^3 * yr / season ", VT, ch );


  // June 4, 2012 - dNdz and rate info moved from README_doc()

  int  DNDZFLAG = 0;
  char ctmp_z[80], ctmp[100], ctmp_pec1a[80], *NAME;
  int  IMODEL_AB, IMODEL_PLAW, IMODEL_PLAW2, IMODEL_ZPOLY, IMODEL_FLAT ;
  int  IMODEL_CCS15, IMODEL_MD14, IMODEL_PISN, IMODEL_TDE ;
  int  IFLAG_REWGT_ZEXP, IFLAG_REWGT_ZPOLY ;
  
  // model flags 
  NAME = INPUTS.RATEPAR.NAME ;
  IMODEL_AB      = ( strcmp(NAME,"ABMODEL"  )  == 0 ) ;
  IMODEL_PLAW    = ( strcmp(NAME,"POWERLAW" )  == 0 ) ;
  IMODEL_PLAW2   = ( strcmp(NAME,"POWERLAW2")  == 0 ) ;
  IMODEL_ZPOLY   = ( strcmp(NAME,"ZPOLY"    )  == 0 ) ;
  IMODEL_FLAT    = ( strcmp(NAME,"FLAT"     )  == 0 ) ;
  IMODEL_CCS15   = ( strcmp(NAME,RATEMODELNAME_CCS15)  == 0 ) ;
  IMODEL_MD14    = ( strcmp(NAME,RATEMODELNAME_MD14)   ==0 );

  IMODEL_PISN  = ( (strcmp(NAME,"IPP") ==0)  ||
		   (strcmp(NAME,RATEMODELNAME_PISN)==0) );

  IMODEL_TDE  = ( (strcmp(NAME,"EPM")==0) ||
		  (strcmp(NAME,RATEMODELNAME_TDE) ==0) ) ;

  // re-wgt flags
  IFLAG_REWGT_ZEXP  = ( INPUTS.RATEPAR.DNDZ_ZEXP_REWGT != 0.0 ) ;
  IFLAG_REWGT_ZPOLY = ( INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT[1] != 0.0  ||
			INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT[2] != 0.0     );

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
    i++;    iz=1;
    sprintf(LINE_RATE_INFO[i],
	    "\t dN/dz = %.2f + %.2f*z + %.2f*z^2 + %.2f*z^3"
	    ,INPUTS.RATEPAR.MODEL_PARLIST[iz][0]
	    ,INPUTS.RATEPAR.MODEL_PARLIST[iz][1]
	    ,INPUTS.RATEPAR.MODEL_PARLIST[iz][2]
	    ,INPUTS.RATEPAR.MODEL_PARLIST[iz][3] );
  }
  else if ( IMODEL_FLAT ) {
    i++;    iz=1;
    sprintf(LINE_RATE_INFO[i], "\t dN/dz = FLAT  ~  (MJDrage * dOmega) ");  

    // Oct 22 2015
    // there is no rate for FLAT option, so make up SNsum so that
    // sim_SNmix will work with FLAT option.
    SNsum = delMJD * (Z1-Z0) * (dOmega * (41000./(2.*TWOPI)) ) ;
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
    sprintf(LINE_RATE_INFO[i],
	    "\t Reweight dN/dz by %.2f + %.2f*z + %.2f*z^2 + %.2f*z^3 "
	    ,INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT[0]
	    ,INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT[1]
	    ,INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT[2]
	    ,INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT[3]   );
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
    INPUTS.NGENTOT_LC = (int)(INPUTS.NGEN_SEASON * TOTsum * SCALE);
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
    
    i++ ; sprintf(LINE_RATE_INFO[i], 
	    "\t dN/d%s = %.2f ", varName, PARLIST[0] );
    for(j=1; j <= MXPOLY_GALRATE; j++ ) {
      if( PARLIST[j] != 0.0 ) {
	i++ ; sprintf(LINE_RATE_INFO[i],
		      "\t    + %14.6le * %s^%d ", PARLIST[j], varName, j );
      }
    }

  }

  // - - - - - - - 

  // this line is needed by sim_SNmix.pl
  i++; sprintf(LINE_RATE_INFO[i],
	       "\t Number of EVENTS per season = %d  (NGENTOT_LC)", 
	       INPUTS.NGENTOT_LC );

  i++; LINE_RATE_INFO[i][0] = 0 ;

  NLINE_RATE_INFO = i+1 ;

  return;


} // end init_DNDB_Rate

// ***************************
void init_simvar(void) {

  // May 26, 2009: one-time init of sim-variables; mostly counters.
  // Nov 24, 2017: init GENLC.MWEBV[_ERR] here instead of in init_GENLC().

  int ifilt, itype;
  float xmb ;

  // ----------- BEGIN -----------
 
  set_GENMODEL_NAME();

  // xxx  GENLC.MODELPATH[0] = 0 ;
  // xxx  GENLC.MODELNAME[0] = 0

  init_GaussIntegral();

  GENLC.STOPGEN_FLAG = 0 ;
  GENLC.ACCEPTFLAG   = GENLC.ACCEPTFLAG_LAST = 0 ;
  
  set_FILTERSTRING(FILTERSTRING);
  set_EXIT_ERRCODE(EXIT_ERRCODE_SIM);

  //  sprintf(FILTERSTRING,"%s", FILTERSTRING_DEFAULT );

  NGENLC_TOT = NGENLC_WRITE = NGENSPEC_WRITE = 0;

  NGEN_REJECT.GENRANGE  = 0;
  NGEN_REJECT.GENMAG    = 0;
  NGEN_REJECT.GENPAR_SELECT_FILE = 0 ;
  NGEN_REJECT.HOSTLIB   = 0;
  NGEN_REJECT.SEARCHEFF = 0;
  NGEN_REJECT.CUTWIN    = 0;
  NGEN_REJECT.NEPOCH    = 0;

  /* xxx
  for(itype=0; itype < MXTYPE; itype++ ) {
    NGENTYPE_TOT[itype] = 0 ;
    NGENTYPE_WRITE[itype] = 0 ;
  }
  xxxxx  */

  GENLC.MWEBV     = 0.0 ;
  GENLC.MWEBV_ERR = 0.0 ;

  GENLC.NTYPE_SPEC_CUTS = 0;
  GENLC.NTYPE_SPEC = 0 ;
  GENLC.NTYPE_PHOT_CUTS = 0;
  GENLC.NTYPE_PHOT = 0 ;

  GENLC.NON1ASED.ISPARSE    = -9 ;  //  (non1a sparse index: 1-NINDEX)

  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ )
    { GENLC.SNRMAX_FILT[ifilt]  = -9.0 ;   }

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

  SIMLIB_OBS_GEN.PIXSIZE[0] = -9.0 ; // ?? why is this here

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

  sprintf(GENMODEL_NAME[MODEL_BYOSED][0],  "%s", "BYOSED"  );

  sprintf(GENMODEL_NAME[MODEL_SIMSED][0],  "%s", "SIMSED"  );

  sprintf(GENMODEL_NAME[MODEL_FIXMAG][0],  "%s", "FIXMAG"  );
  sprintf(GENMODEL_NAME[MODEL_FIXMAG][1],  "%s", "RANMAG"  );
  sprintf(GENMODEL_NAME[MODEL_FIXMAG][2],  "%s", "fixmag"  );
  sprintf(GENMODEL_NAME[MODEL_FIXMAG][3],  "%s", "ranmag"  );

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
void init_simRandoms(void) {

  // Create Sep 2016 by R.Kessler & E.Jennings
  // Move init stuff from main, and check for skewNormal.

  int i;
  int ISEED = INPUTS.ISEED ;
  //  char fnam[] = "init_simRandoms" ;

  // ----------- BEGIN ----------------
  srandom(ISEED);
  init_RANLIST(); 
  for ( i=1; i <= NLIST_RAN; i++ )  
    { RANFIRST[i] = FlatRan1(i); }
  
  // ---------------- skewNormal stuff -------------------
  //
  init_skewNormal(ISEED);  // one-time init, to set seed in python

} // end init_simRandoms

// ************************************
void  init_GENLC(void) {

  // called for each SN.
  // Oct    2011: return if SIMLIB_IDLOCK is set
  // Jan    2017: init GENLC.SNTYPE=0 instead of -9
  // Apr 6, 2017: move more stuff before USE_SAME_SIMLIB call.
  // Nov 06 2017: init SEARCHEFF_RANDOMS structure
  // Nov 22 2017: use NEP_RESET (instead of MXEPSIM) to limit
  //               wasted CPU on initializing (based on gprof)
  //
  int epoch, ifilt, ifilt_obs, i, imjd, NEP_RESET ;
  char fnam[] = "init_GENLC" ;

  // -------------- BEGIN ---------------

  GENLC.NOBS           = 0 ;
  GENLC.NOBS_REMOVE    = 0 ;
  GENLC.NOBS_UNDEFINED = 0 ;
  GENLC.NOBS_SATURATE[0]  = 0 ; // NOBS NOT saturated
  GENLC.NOBS_SATURATE[1]  = 0 ; // NOBS saturated

  GENLC.SEARCHEFF_MASK = 0 ;
  GENLC.SEARCHEFF_SPEC = 0.0 ;
  GENLC.SEARCHEFF_zHOST= 0.0 ;
  GENLC.CUTBIT_MASK    = 0 ;
  GENLC.METHOD_TYPE    = NULLINT ;
  GENLC.SNTYPE         = 0; 

  GENLC.PEAKMAG_TRIGGER_FLAG = 0 ;
  GENLC.ACCEPTFLAG_LAST  = GENLC.ACCEPTFLAG ;
  GENLC.ACCEPTFLAG       = 0 ;
  GENLC.ACCEPTFLAG_FORCE = 0 ;
  
  GENLC.CORRECT_HOSTMATCH = 1 ;  // default is correct match
  GENLC.SALT2gammaDM  = 0.0 ;

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    GENLC.DOFILT[ifilt_obs] = 1;  // do all user-filters by default.
  }


  // - - - - -
  // init shape-parameters to NULLFLOAT;
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

  GENLC.TRESTMIN   = NULLFLOAT ;
  GENLC.TRESTMAX   = NULLFLOAT ;

  GENLC.SNRMAX_GLOBAL    = -9.0 ;
  GENLC.IEPOCH_SNRMAX    = -9 ;
  GENLC.IEPOCH_NEARPEAK  = -9 ;
  GENLC.DTPEAK_MIN       = 99999. ;

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
    GENLC.NON1AGRID_RANGauss = GaussRanClip(ILIST_RAN, rmin, rmax);
    GENLC.NON1AGRID_RANFlat  = FlatRan1(ILIST_RAN); 
  }


  // init filter-dependent stuff.
  for ( ifilt_obs=0; ifilt_obs < MXFILTINDX; ifilt_obs++ ) {
    GENLC.peakmag8_obs[ifilt_obs]   = NULLFLOAT ;    

    GENLC.IEPOCH_PEAK[ifilt_obs]     = -9 ; 
    GENLC.IEPOCH_TEMPLATE[ifilt_obs] = -9 ; 
    GENLC.genmag8_obs_template[ifilt_obs] = 99.0 ; // zero flux in template

    if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) 
      { GENLC.SIMLIB_USEFILT_ENTRY[ifilt_obs] = 1 ; }
    else
      { GENLC.SIMLIB_USEFILT_ENTRY[ifilt_obs] = 0 ; }

    
    GENLC.SNRMAX_FILT[ifilt_obs]    = -9.0 ;  
    GENLC.SNRMAX_SORTED[ifilt_obs]  = -9.0 ;  

    SEARCHEFF_RANDOMS.SPEC[ifilt_obs] = -9.0 ;

    GENLC.NOBS_FILTER[ifilt_obs] = 0 ;
    GENLC.NOBS_SATURATE_FILTER[0][ifilt_obs] = 0 ;
    GENLC.NOBS_SATURATE_FILTER[1][ifilt_obs] = 0 ;

    GENLC.FLUXCOR_MWEBV_MAP[ifilt_obs]   = 1.0 ;
    GENLC.MAGCOR_MWEBV_MAP[ifilt_obs]    = 0.0 ;
    GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs]   = 0.0 ;
  }
  

  if ( NGENLC_TOT == 1 ) 
    { NEP_RESET = MXEPSIM ; }
  else
    { NEP_RESET = GENLC.NEPOCH + 1; }

  /*
  printf(" xxx NGENLC_TOT=%d  NEP_RESET=%d   (NEPOCH=%d)\n",
	 NGENLC_TOT, NEP_RESET, GENLC.NEPOCH );
  */

  // reset data struct for trigger evaluation (Jan 2014)
  SEARCHEFF_DATA.NOBS        =  0 ;
  SEARCHEFF_DATA.REDSHIFT    = -9.0 ;
  SEARCHEFF_DATA.PEAKMJD     = -9.0 ;
  SEARCHEFF_DATA.DTPEAK_MIN  = -9.0 ;
  SEARCHEFF_DATA.SALT2mB     = -9.0 ;
  int obs;
  //  for(obs=0; obs < MXOBS_TRIGGER; obs++ ) {
  for(obs=0; obs < NEP_RESET; obs++ ) {
    SEARCHEFF_DATA.IFILTOBS[obs]    = -9   ;
    SEARCHEFF_DATA.MJD[obs]         = -9.0 ;
    SEARCHEFF_DATA.MAG[obs]         = -9.0 ;
    SEARCHEFF_DATA.SNR[obs]         = -9.0 ;
    SEARCHEFF_DATA.detectFlag[obs]  =  0   ;
    SEARCHEFF_DATA.PHOTPROB[obs]    =  0.0 ;
    SEARCHEFF_RANDOMS.PIPELINE[obs] = -9.0 ;
    SEARCHEFF_RANDOMS.PHOTPROB[obs] = -999.0 ;
  }
  

  GENLC.MAGSMEAR_COH = 0.0 ;
     
  for ( epoch = 0; epoch < NEP_RESET; epoch++ ) {
    // xxx mark delete    GENLC.cc[epoch]   = NULLINT ;
    GENLC.IFILT_OBS[epoch] = NULLINT ;

    GENLC.USE_EPOCH[epoch] = 0 ;

    sprintf(GENLC.FIELDNAME[epoch], "NULL" );

    GENLC.MJD[epoch]          = NULLFLOAT ;
    GENLC.epoch8_obs[epoch]   = NULLFLOAT ;
    GENLC.epoch8_rest[epoch]  = NULLFLOAT ;
    GENLC.epoch8_obs_range[0] = +999999.0 ;
    GENLC.epoch8_obs_range[1] = -999999.0 ;

    GENLC.flux[epoch]         = NULLFLOAT ;
    GENLC.flux_errstat[epoch] = NULLFLOAT ;
    GENLC.mag[epoch]          = NULLFLOAT ;
    GENLC.mag_err[epoch]      = NULLFLOAT ;

    GENLC.genmag8_obs[epoch]   = NULLFLOAT ;
    GENLC.generr8_obs[epoch]   = NULLFLOAT ; // Apr 2013

    GENLC.genmag8_rest[epoch]  = NULLFLOAT ;
    GENLC.generr8_rest[epoch]  = 0.000 ;
    GENLC.genmag8_rest2[epoch]  = NULLFLOAT ;
    GENLC.generr8_rest2[epoch]  = NULLFLOAT ;

    GENLC.kcorval8[epoch]     = NULLFLOAT ;
    GENLC.warpcolval8[epoch]  = NULLFLOAT ;
    GENLC.AVwarp8[epoch]      = NULLFLOAT ;
    sprintf(GENLC.kcornam[epoch],    "NULL" );
    sprintf(GENLC.warpcolnam[epoch], "NULL" );

    GENLC.ISPEAK[epoch]      = 0;    
    GENLC.ISTEMPLATE[epoch]  = 0;    
    GENLC.SNR_CALC[epoch]    = 0.0 ;
    GENLC.SNR_MON[epoch]     = 0.0 ;

  } // end of epoch loop

  // Aug 9 2014: moved from gen_cutwin()
  GENLC.NOBS_SNR    = 0 ;
  GENLC.NOBS_MJDDIF = 0 ;
  GENLC.TRESTMIN   = +99999. ;
  GENLC.TRESTMAX   = -99999. ;
  GENLC.TGAPMAX    = -99999. ;  // max rest-frame gap among all filters
  GENLC.T0GAPMAX   = -99999. ;  // idem, near peak

  for ( i=0; i < MXZRAN; i++ )
    {  GENLC.REDSHIFT_RAN[i] = -9. ; }


  // init GENSPEC stuff
  if ( INPUTS.SPECTROGRAPH_OPTIONS.DOFLAG_SPEC ) {
    GENSPEC.NMJD_TOT = 0 ;
    for(imjd=0; imjd < MXSPEC; imjd++ )  { GENSPEC_INIT(1,imjd); }
  }

  // keep track of last coord to skip parts of gen_MWEBV
  GENLC.RA_LAST   = GENLC.RA ;
  GENLC.DEC_LAST  = GENLC.DEC ;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // skip rest of init if we are going to use the same LIBID
  if ( USE_SAME_SIMLIB_ID(1) != 0 ) { return ; } // Dec 2015
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  GENLC.NGEN_SIMLIB_ID = 0 ;
  GENLC.NEPOCH   = 0 ;
  GENLC.NEWMJD   = 0 ;

  GENLC.CID       = NULLINT ;
  GENLC.PEAKMJD   = NULLFLOAT ;
  GENLC.MJD_EXPLODE = NULLFLOAT ;
  GENLC.ISOURCE_PEAKMJD = -9 ;
  GENLC.SDSS_SIM  = 0 ;
  GENLC.RA        = -999. ;
  GENLC.DEC       = -999. ;


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
    SIMLIB_OBS_GEN.ZPTADU[i]     = -99. ;
    SIMLIB_OBS_GEN.ZPTSIG[i]     = -99. ;
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

  double GENMODEL_ERRSCALE = (double)INPUTS.GENMODEL_ERRSCALE ;
  double SMEAR_SCALE       = (double)INPUTS.GENMAG_SMEAR_SCALE;
  int    OPT, j, USE_SALT2smear ;
  double LAMRANGE[2], SIGCOH,  PARLIST_SETPKMJD[10];
  char *ptrName, key[40], NAM3[8], PATHPLUSMODEL[MXPATHLEN] ;
  char MODELPATH_SALT2[MXPATHLEN];
  char fnam[] = "init_modelSmear"  ;

  // --------- BEGIN ----------

  GENLC.NRANGauss_GENSMEAR = 0 ;
  GENLC.NRANFlat_GENSMEAR  = 0 ;


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

  // check model names
  init_genSmear_FLAGS(SMEAR_SCALE); // internal inits

  sprintf(key,"GENMAG_SMEAR_MODELNAME") ;
  ptrName = INPUTS.GENMAG_SMEAR_MODELNAME ;
  if ( IGNOREFILE(ptrName) ) { goto INIT_SETPKMJD ; }

  INPUTS.DO_MODELSMEAR  = 1 ;

  if ( GENMODEL_ERRSCALE > 1.0E-9 ) {
    printf("\n PRE-ABORT DUMP: \n");
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
    
    if ( INDEX_GENMODEL==MODEL_SALT2 ) 
      { sprintf(MODELPATH_SALT2,"%s", INPUTS.MODELPATH ); }
    else
      { sprintf(MODELPATH_SALT2,"%s",INPUTS.GENMAG_SMEAR_MODELARG);} //BYOSED
    
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
    init_genSmear_SALT2(MODELPATH_SALT2, ptrName, SIGCOH);
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

  else if ( strcmp(ptrName,"COH") == 0 ) 
    {  init_genSmear_COH() ; }

  else if ( strcmp(ptrName,"BIMODAL_UV") == 0 ) 
    {  init_genSmear_biModalUV() ; }

  else if ( strcmp(ptrName,"PRIVATE") == 0 ) 
    {  init_genSmear_private() ; }

  else {
    sprintf(c1err,"Invalid smear model: '%s'", ptrName);
    sprintf(c2err,"Check GENMAG_SMEAR_MODELNAME key.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // fetch & store number of randoms to generate for each SN.

  int *NG = &GENLC.NRANGauss_GENSMEAR ;
  int *NF = &GENLC.NRANFlat_GENSMEAR ;

  get_NRAN_genSmear(NG, NF); // returns NG and NF    

  if ( *NG > MXFILTINDX ) {
    sprintf(c1err, "NRANGauss=%d exceeds bound of %d", *NG, MXFILTINDX);
    sprintf(c2err, "Check %s and its parameters.", ptrName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  if ( *NF > MXFILTINDX ) {
    sprintf(c1err, "NRANFlat=%d exceeds bound of %d", *NF, MXFILTINDX);
    sprintf(c2err, "Check %s and its parameters.", ptrName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  
  printf("   Number of %s Gaussian  randoms per SN: %d \n",  ptrName, *NG);
  printf("   Number of %s Flat(0-1) randoms per SN: %d \n",  ptrName, *NF);
  printf("   MagSmear scale: %.3f \n", SMEAR_SCALE);
  fflush(stdout);

  // -------------------------------
  if ( INPUTS.DO_MODELSMEAR  == 0 ) 
    { printf("\t ==> No model smearing options selected. \n" ) ; }
  else{ 
    // print magSmear sigma at various wavelengths
    //    dump_modelSmearSigma();
  }

  
  // May 2019: init method to estimate peakmjd for data files
 INIT_SETPKMJD:
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

  int NLAM, NRANGEN, igen, ilam ;
  double LAMMIN, LAMMAX, LAMBIN, LAMARRAY[100], LAM, TREST;
  double **MAGSMEAR, MAGARRAY[100], AVG, RMS ;
  char fnam[] = "dump_modelSmearSigma" ;

  // --------------- BEGIN --------------

  sprintf(BANNER,"%s: Determine sigma(magSmear) vs. wavelength", fnam);
  print_banner(BANNER);

  LAMMIN = 2000.0 ;
  LAMMAX = GENLC.RESTLAM_MODEL[1] + 3000.0 ;
  LAMBIN = 500.0 ;

  TREST = 0.0 ; 
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


  for(igen=0; igen < NRANGEN; igen++ ) {
    init_RANLIST();      // init list of random numbers 
    genran_modelSmear(); // load randoms for genSmear
    get_genSmear(TREST, NLAM, LAMARRAY, MAGARRAY); // return MAGARRAY


    for(ilam=0; ilam < NLAM; ilam++ ) 
      { MAGSMEAR[ilam][igen] = MAGARRAY[ilam]  ; }

  } // end of igen loop
  

  for(ilam=0; ilam < NLAM; ilam++ ) {
    LAM = LAMARRAY[ilam] ;
    arrayStat(NRANGEN, MAGSMEAR[ilam], &AVG, &RMS);
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
       INDEX_GENMODEL == MODEL_BYOSED         ) {
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

  // this know is tunable, not a fixed value from mask
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

  // - - - - - - - - - - - - - - - - - - - -
  // malloc arrays vs. lambda 
  int ispec, MXLAM = NBLAM ;
  int MEMD = MXLAM * sizeof(double);
  for(ispec=0; ispec < MXSPEC; ispec++ ) {
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

  // check if any spectrum has a warp
  if ( INPUTS.NWARP_TAKE_SPECTRUM > 0 ) 
    {  GENSPEC.USE_WARP = 1 ; }
  else
    {  GENSPEC.USE_WARP = 0 ; }


  // allocate arrays for Gaussian template noise
  for(ilam=0; ilam < MXLAMSMEAR_SPECTROGRAPH; ilam++ ) {
    GENSPEC.RANGauss_NOISE_TEMPLATE[ilam] = (double*) malloc(MEMD) ;
  }

  
  return ;

} // end init_genSpec

// *****************************************************
void GENSPEC_DRIVER(void) {

  // Created July 2016
  // Driver for generating spectra.
  //
  // Sep  1 2016: add flat spectra for FIXRAN model.
  // Oct 16 2016: apply Gaussian LAMRES smearing

  double MJD, MJD_LAST, MJD_DIF, TEXPOSE, TOBS, TREST ;
  double z = GENLC.REDSHIFT_CMB,  z1 = 1.0 + z;
  int  ep, i, imjd, NMJD ;
  int  INDX, OPT_TEXPOSE ;
  char fnam[] = "GENSPEC_DRIVER" ;

  // ------------ BEGIN ------------

  GENSPEC.NMJD_PROC = 0;

  if ( INPUTS.SPECTROGRAPH_OPTIONS.DOFLAG_SPEC == 0 ) { return ; }

  // Jan 2018: allow some LIBIDs to NOT have a spectrograph key
  if ( NPEREVT_TAKE_SPECTRUM == 0 && SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH == 0 ) 
    { return ; }


  GENSPEC.NBLAM_TOT = INPUTS_SPECTRO.NBIN_LAM ;

  MJD_LAST = MJD_DIF = -9.0 ;
  NMJD = 0 ;


  NMJD = GENSPEC.NMJD_TOT  ;
  if ( NPEREVT_TAKE_SPECTRUM > 0 && NPEREVT_TAKE_SPECTRUM != NMJD ) {
    printf("\n PRE-ABORT DUMP: \n");
    printf("  NPEREVT_TAKE_SPECTRUM = %d \n", NPEREVT_TAKE_SPECTRUM );
    printf("  NMJD_TOT = %d \n", NMJD);
    sprintf(c1err,"Cannot mix TAKE_SPECTRUM keys in sim-input file");
    sprintf(c2err,"with SPECTROGRAPH keys in SIMLIB");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // - - - - - - - - - - - - - - - - 
  int imjd_order[MXSPEC];
  double SNR, LAMMIN=100.0, LAMMAX=25000. ; 

  GENSPEC_MJD_ORDER(imjd_order); // check if nearPeak is done first

  for(i=0; i < NMJD; i++ ) {

    imjd = imjd_order[i];

    GENSPEC_INIT(2,imjd);   // 2-> event-dependent init

    // compute true GENMAG and FLUXGEN in each lambda bin
    GENSPEC_TRUE(imjd); 

    // July 2019: option to add host contamination
    GENSPEC_HOST_CONTAMINATION(imjd);

    // apply optional fudges for test or debug
    GENSPEC_FUDGES(imjd); 

    // if TAKE_SPECTRUM is defined by SNR, compute TEXPOSE
    GENSPEC_TEXPOSE_TAKE_SPECTRUM(imjd);

    // smear fluxes from Poisson noise and wavelength
    SNR = GENSPEC_SMEAR(imjd,LAMMIN,LAMMAX);

    // Feb 2 2017: convert flux to FLAM (dF/dlam)
    GENSPEC_FLAM(imjd) ;

    GENSPEC.NMJD_PROC++ ; // total Nspec for this event.
    NGENSPEC_WRITE++ ;    // total NSpec over all Light curve

  } // end imjd
  


  return ;

} // end GENSPEC_DRIVER


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

  char fnam[] = "GENSPEC_LAMOBS_RANGE" ;

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
    { EPOCH = FlatRan(ILIST_RAN, EPOCH_RANGE); }  // pick random epoch

  z1 = 1.0 + z ; 
  if ( OPT_FRAME == GENFRAME_OBS ) 
    { Tobs = EPOCH ;    Trest = EPOCH/z1; }
  else
    { Tobs = EPOCH*z1 ; Trest = EPOCH ; }

  /*
  printf(" xxx  z1=%.3f  OPT=%d  INDX=%d  EPRANGE=(%5.1f,%5.1f) "
	 "EP=%.1f  Tobs=%.1f\n", 
	 z1, OPT_MJD, INDX,
	 EPOCH_RANGE[0], EPOCH_RANGE[1], EPOCH, Tobs ); fflush(stdout);
  */


  // load output arguments
  MJD = GENLC.PEAKMJD + Tobs;
  *TOBS  = Tobs ;
  *TREST = Trest ;

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

  if ( OPT == 1 ) {
    GENSPEC.MJD_LIST[imjd]      = -9.0 ;
    GENSPEC.TOBS_LIST[imjd]     = -9.0 ;
    GENSPEC.TREST_LIST[imjd]    = -9.0 ;
    GENSPEC.TEXPOSE_LIST[imjd]        = -9.0 ;
    GENSPEC.OPT_TEXPOSE_LIST[imjd]    = -9 ;
    GENSPEC.INDEX_TAKE_SPECTRUM[imjd] = -9 ;

    GENSPEC.SNR_REQUEST_LIST[imjd] = -9.0 ;
    GENSPEC.SNR_COMPUTE_LIST[imjd] = -99.0 ;
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

  int  NBLAM      = INPUTS_SPECTRO.NBIN_LAM ;
  int  IS_HOST    = GENSPEC.IS_HOST[imjd];
  double TOBS     = GENSPEC.TOBS_LIST[imjd];
  double TREST    = GENSPEC.TREST_LIST[imjd];

  double GENMAG, ZP, ARG, FLUXGEN, MAGOFF ;
  int ilam ;
  char fnam[] = "GENSPEC_TRUE" ;

  // --------------- BEGIN ----------------


  if ( IS_HOST ) {    
    genSpec_HOSTLIB(GENLC.REDSHIFT_HELIO,         // (I) helio redshift
		    GENLC.MWEBV,                  // (I) Galactic extinction
		    GENSPEC.GENFLUX_LIST[imjd],   // (O) fluxGen per bin 
		    GENSPEC.GENMAG_LIST[imjd] );  // (O) magGen per bin
    return;
  }

  // below is a SN spectrum, so check which model.

  if ( INDEX_GENMODEL == MODEL_SALT2 )  {
    genSpec_SALT2(GENLC.SALT2x0, GENLC.SALT2x1, GENLC.SALT2c, 
		  GENLC.MWEBV,             // Galactic
		  GENLC.RV, GENLC.AV,      // host
		  GENLC.REDSHIFT_HELIO, TOBS,
		  GENSPEC.GENFLUX_LIST[imjd], // (O) fluxGen per bin 
		  GENSPEC.GENMAG_LIST[imjd]   // (O) magGen per bin
		  );
  }
  else if ( INDEX_GENMODEL == MODEL_FIXMAG )  {
    for(ilam=0; ilam < NBLAM; ilam++ ) {
      GENMAG  = GENLC.FIXMAG ;   
      ZP      = ZEROPOINT_FLUXCAL_DEFAULT ;
      ARG     = -0.4*(GENMAG-ZP);
      FLUXGEN = pow(TEN,ARG);
      
      GENSPEC.GENMAG_LIST[imjd][ilam]  = GENMAG ;
      GENSPEC.GENFLUX_LIST[imjd][ilam] = FLUXGEN ;
    }
  }
  else  if ( INDEX_GENMODEL  == MODEL_NON1ASED ) {    
    // works with NON1ASED, but NOT with NON1AGRID
    MAGOFF = GENLC.GENMAG_OFF_GLOBAL + GENLC.MAGSMEAR_COH ;
    getSpec_SEDMODEL(ISED_NON1A,
		     GENLC.MWEBV, GENLC.RV, GENLC.AV, // (I)		     
		     GENLC.REDSHIFT_HELIO,            // (I) redshift
		     GENLC.DLMU,                      // (I) dist mod
		     TOBS, MAGOFF,                    // (I) Tobs, magoff 
		     GENSPEC.GENFLUX_LIST[imjd],   // (O) fluxGen per bin 
		     GENSPEC.GENMAG_LIST[imjd]     // (O) magGen per bin
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
		     GENLC.MWEBV, GENLC.RV, GENLC.AV, // (I)		
		     GENLC.REDSHIFT_HELIO,            // (I) redshift
		     GENLC.DLMU,
		     TOBS, MAGOFF,                   // (I) Tobs, magoff 
		     GENSPEC.GENFLUX_LIST[imjd],     // (O) fluxGen per bin 
		     GENSPEC.GENMAG_LIST[imjd]       // (O) magGen per bin
		     ); 
  }
  else if ( INDEX_GENMODEL == MODEL_BYOSED ) {
    
    
    genSpec_BYOSED(TOBS, 
		   GENLC.REDSHIFT_HELIO,            // (I) helio redshift
		   GENLC.DLMU,                      // (I) dist. mod.
		   GENLC.MWEBV, GENLC.RV, GENLC.AV, // (I)		     
		   GENSPEC.GENFLUX_LIST[imjd],      // (O) fluxGen per bin 
		   GENSPEC.GENMAG_LIST[imjd]        // (O) magGen per bin
		   );		
   
  }
  else { 
    /*  don't abort since init_genSpec gives warning.  
	sprintf(c1err,"Cannot make spectrum for model=%s", modelName);
	sprintf(c2err,"Check manual for valid models");
	errmsg(SEV_WARN, 0, fnam, c1err, c2err ); 
    */
  }

  return ;

} // end GENSPEC_TRUE


// *****************************************
void GENSPEC_HOST_CONTAMINATION(int imjd) {

  // Created July11 2019
  // Check option to add host contamination
  int    IS_HOST  = GENSPEC.IS_HOST[imjd];
  double HOSTFRAC = (double)INPUTS.TAKE_SPECTRUM_HOSTFRAC;
  int    IMJD_HOST = 0;
  int    NBLAM = INPUTS_SPECTRO.NBIN_LAM ;

  int ilam;
  double FLAM_SN, FLAM_HOST, FLAM_TOT, arg, MAGSHIFT ;
  char fnam[] = "GENSPEC_HOST_CONTAMINATION" ;

  // ------------- BEGIN --------------

  if ( IS_HOST )           { return; }
  if ( HOSTFRAC < 1.0E-8 ) { return; }

  // check that imjd=0 is indeed a host spectrum to add
  IS_HOST = GENSPEC.IS_HOST[IMJD_HOST];
  if ( !IS_HOST ) {
    sprintf(c1err,"IMJD=%d is not a host spectrum", IMJD_HOST);
    sprintf(c2err,"Cannot add host contamination to SN spectrum.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  for(ilam=0; ilam < NBLAM; ilam++ ) {
    FLAM_SN   = GENSPEC.GENFLUX_LIST[imjd][ilam];
    FLAM_HOST = GENSPEC.GENFLUX_LIST[IMJD_HOST][ilam];
    FLAM_TOT  = FLAM_SN + (FLAM_HOST*HOSTFRAC);
    
    GENSPEC.GENFLUX_LIST[imjd][ilam] = FLAM_TOT ;

    arg      =  FLAM_TOT/FLAM_SN ;
    MAGSHIFT = -2.5*log10(arg);
    GENSPEC.GENMAG_LIST[imjd][ilam] += MAGSHIFT ;
  }

  return;

} // end GENSPEC_HOST_CONTAMINATION

// *************************************************
void GENSPEC_TEXPOSE_TAKE_SPECTRUM(int imjd) {

  // Compute TEXPOSE from requested SNR.
  // For synthetic filters, update ZPT and SKYSIG.

  int  LDMP   = (GENLC.CID == INPUTS.TAKE_SPECTRUM_DUMPCID );

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
  double z           = GENLC.REDSHIFT_CMB ;  
  double z1          = 1.0 + z ;
  double TOBS        = MJD - GENLC.PEAKMJD ;
  double TREST       = TOBS/z1;

  double LAMMIN_OBS, LAMMAX_OBS ;
  double SNR, SNR_REQUEST, TEXPOSE, TEXPOSE_T, SKYSIG, SKYSIG_T;
  double ZPT, ZPT_T, PSFSIG ;
  char fnam[] = "GENSPEC_TEXPOSE_TAKE_SPECTRUM" ;
  
  // ------------ BEGIN --------------

  if ( INDX < 0 )         { return ; } // not from TAKE_SPECTRUM key
  if ( OPT_TEXPOSE != 2 ) { return ; } // not SNR option

  SNR_REQUEST = eval_GENPOLY(z,GENZPOLY_SNR,fnam);
  GENSPEC.SNR_REQUEST_LIST[imjd] = SNR_REQUEST ;

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
    printf(" xxx SNR_MIN(%d sec)=%.1f, SNR_MAX(%d sec)=%.1f \n",
	    (int)TEXPOSE_MIN, SNR_MIN, (int)TEXPOSE_MAX, SNR_MAX);
    printf(" xxx OPT_FRAME_EPOCH=%d(%s)  OPT_FRAME_LAMBDA=%d \n",
	   OPT_FRAME_EPOCH, FRAME_EPOCH, OPT_FRAME_LAMBDA );
    fflush(stdout);
  }

  if ( SNR_REQUEST <= SNR_MIN ) {
    GENSPEC.TEXPOSE_LIST[imjd]     = TEXPOSE_MIN ;
    GENSPEC.SNR_COMPUTE_LIST[imjd] = SNR_MIN ;
    goto CLEANUP ;
  }
  if ( SNR_REQUEST >= SNR_MAX ) {
    GENSPEC.TEXPOSE_LIST[imjd]     = TEXPOSE_MAX ;
    GENSPEC.SNR_COMPUTE_LIST[imjd] = SNR_MAX ;
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
  double tol=99.0 , tol_converge = 0.02 ;
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
      if ( GENSPEC.NMJD_PROC == 0 ) {
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
    if ( TEXPOSE <= TEXPOSE_MIN ) { TEXPOSE = TEXPOSE_MIN; }
    if ( TEXPOSE >= TEXPOSE_MAX ) { TEXPOSE = TEXPOSE_MAX; }
    
    // for GENSPEC_SMEAR, check option to set global template 
    // exposure time by scaling T_expose for epoch nearest peak.
    if ( GENSPEC.NMJD_PROC == 0 && SCALE > 0.01 )  { 
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
      printf("\n PRE-ABORT DUMP (CID=%d): \n", GENLC.CID);
      printf("\t SNR_MIN(%d sec)=%.1f \n", (int)TEXPOSE_MIN, SNR_MIN ) ;
      printf("\t SNR_MAX(%d sec)=%.1f \n", (int)TEXPOSE_MAX, SNR_MAX ) ;

      sprintf(c1err,"Could not converge after NITER=%d", NITER);
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
   
    // set optional template noise for passbands made from spectrograph
    if ( TEXPOSE_T > 0.001 ) {
      SIMLIB_TEMPLATE.USEFLAG |= 2;
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

  int    NBLAM = INPUTS_SPECTRO.NBIN_LAM ;
  int    ilam, ILAM_MIN=99999, ILAM_MAX=-9, NBLAM_USE=0 ;
  double GENFLUX, GENFLUXERR, GENFLUXERR_T, GENMAG, SNR_true ;
  double ERRFRAC_T, LAMAVG ; 

  double  TEXPOSE_S  = GENSPEC.TEXPOSE_LIST[imjd] ;
  double  TEXPOSE_T  = GENSPEC.TEXPOSE_TEMPLATE ;
  double  SCALE_SNR  = INPUTS.SPECTROGRAPH_OPTIONS.SCALE_SNR ;
  double  SNR_SPEC ;

  char   fnam[] = "GENSPEC_SMEAR" ;

  // ------------- BEGIN -------------


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
    SNR_true = getSNR_spectrograph(ilam, TEXPOSE_S, TEXPOSE_T, GENMAG, 
				   &ERRFRAC_T);  // template frac of error

    if ( SCALE_SNR != 1.00 ) { SNR_true *= SCALE_SNR ;  }
    
    GENFLUXERR    = GENFLUX/SNR_true ;       // sigma on FLUXGEN
    GENFLUXERR_T  = GENFLUXERR * ERRFRAC_T;  // template contribution

    // apply lambda smear to distribute GENFLUX over lambda bins 
    GENSPEC_LAMSMEAR(imjd, ilam, GENFLUX, GENFLUXERR, GENFLUXERR_T );

    NBLAM_USE++ ;

  } // end ilam  

  if ( NBLAM_USE == 0 ) { return(0.0); }

  // - - - - - - - - - - - - - - 
  // convert sum(ERRSQ) -> ERR in each lambda bin and
  // count number of wavelength bins with flux;
  double ERRSQ, SUM_FLUX, SUM_ERRSQ, SUM_ERR ;
  SUM_FLUX = SUM_ERRSQ = 0.0 ;

  for(ilam=ILAM_MIN; ilam <=ILAM_MAX ; ilam++ ) {
    ERRSQ = GENSPEC.OBSFLUXERRSQ_LIST[imjd][ilam] ;
    if ( ERRSQ <= 1.0E-100 ) { continue ; }

    GENSPEC.OBSFLUXERR_LIST[imjd][ilam] = sqrt(ERRSQ);
    GENSPEC.NBLAM_VALID[imjd]++ ; 

    SUM_FLUX    += GENSPEC.OBSFLUX_LIST[imjd][ilam] ;
    SUM_ERRSQ   += ERRSQ ;
  }

  SUM_ERR  = sqrt(SUM_ERRSQ);
  SNR_SPEC = SUM_FLUX / SUM_ERR ;


  /*
  printf("\t xxx %s: SNR = %le / %le = %le \n", 
	 fnam, SUM_FLUX, sqrt(SUM_ERRSQ), SNR_SPEC); fflush(stdout);
  */

  return(SNR_SPEC) ;

} // end GENSPEC_SMEAR

// *********************************************
void  GENSPEC_LAMSMEAR(int imjd, int ilam, double GenFlux, 
		       double GenFluxErr, double GenFluxErr_T ) {

  // Use Gaussian lambea-resolution to smear flux(imjd,ilam)
  // over nearby lambda bins.  In each lambda bin, use
  // GenFluxErr to add Poisson noise.
  // Inputs
  //  + imjd = sparse MJD index for spectra
  //  + ilam = wavelength index
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
  //

  int OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;
  int onlyTNOISE = ( OPTMASK & SPECTROGRAPH_OPTMASK_onlyTNOISE ) ;
  int noTNOISE   = ( OPTMASK & SPECTROGRAPH_OPTMASK_noTEMPLATE ) ;
  int noNOISE    = ( OPTMASK & SPECTROGRAPH_OPTMASK_noNOISE    ) ;
  int ILIST_RAN  = ILIST_RANDOM_SPECTROGRAPH ;
  double NSIGLAM, LAMAVG, LAMSIGMA, LAMBIN, LAMSIG0, LAMSIG1;
  double GINT, SUM_GINT, GINT_SQRT, GRAN_S, GRAN_T ;
  double tmp_GenFlux, tmp_GenFluxErr, tmp_GenFluxErr_S, tmp_GenFluxErr_T ;
  double tmp_RanFlux_S, tmp_RanFlux_T;
  double GenFluxErr_S, OBSFLUX, OBSFLUXERR ;
  int    NBIN2, ilam2, ilam_tmp, NBLAM, NRAN;
  char fnam[] = "GENSPEC_LAMSMEAR" ;

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

  NBIN2    = (int)(NSIGLAM*LAMSIGMA/LAMBIN) ;  
  SUM_GINT = 0.0 ; 
  NRAN     = 0 ;

  // loop over neighbor bins to smear flux over lambda bins
  for(ilam2=ilam-NBIN2; ilam2 <= ilam+NBIN2; ilam2++ ) {
    if ( ilam2 <  0     ) { continue ; }
    if ( ilam2 >= NBLAM ) { continue ; }
    
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

      if ( GENSPEC.NMJD_PROC==0 && tmp_GenFluxErr_T > 0.0 ) 
	{ GENSPEC.RANGauss_NOISE_TEMPLATE[NRAN][ilam] = GaussRan(ILIST_RAN);}

      GRAN_S = GaussRan(ILIST_RAN);
      GRAN_T = GENSPEC.RANGauss_NOISE_TEMPLATE[NRAN][ilam] ;

      // random noise from search 
      tmp_RanFlux_S = tmp_GenFluxErr_S * GRAN_S ;
      if ( onlyTNOISE ) { tmp_RanFlux_S = 0.0 ; }
      
      // correlated random noise from template
      tmp_RanFlux_T = tmp_GenFluxErr_T * GRAN_T ;
      if ( noTNOISE ) { tmp_RanFlux_T = 0.0 ; }
    }

    /* xxxxxxxx
    if ( (ilam > 98 && ilam < 103) && ilam2==ilam ) {
      printf(" xxx ilam=%3d imjd=%d  noise(S,T) = %10.3le , %10.3le \n",
	     ilam, imjd, tmp_RanFlux_S, tmp_RanFlux_T );
    }
    xxxxxxxxxxx */

    // add noise to true flux
    OBSFLUX    = tmp_GenFlux + tmp_RanFlux_S + tmp_RanFlux_T ;
    OBSFLUXERR = tmp_GenFluxErr ; // naive obs-error is true error

    // increment sum of obsFlux and sum of error-squared.
    GENSPEC.OBSFLUX_LIST[imjd][ilam2]      += OBSFLUX ;
    GENSPEC.OBSFLUXERRSQ_LIST[imjd][ilam2] += (OBSFLUXERR*OBSFLUXERR) ;
    GENSPEC.GENFLUX_LAMSMEAR_LIST[imjd][ilam2] += tmp_GenFlux ;

    NRAN++ ;

    if ( NRAN >= MXLAMSMEAR_SPECTROGRAPH ) {
      sprintf(c1err,"NLAMSMEAR = %d exceeds bound of %d",
	      NRAN, MXLAMSMEAR_SPECTROGRAPH );
      sprintf(c2err,"ilam=%d LAMAVG=%.2f", ilam, LAMAVG );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  } // end ilam2

  return ;

} // end GENSPEC_LAMSMEAR

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

    if ( LDMP && fabs(LAMAVG-7000.0) < 20.0 ) {
      printf(" xxx \t LAM=%.1f  --> WARP = %8.5f \n", LAMAVG, WARP);
      fflush(stdout);
    }

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

      if ( GENLC.ISPEAK[ep]     ) { continue ; }
      if ( GENLC.ISTEMPLATE[ep] ) { continue ; }

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

  Oct 18 2017: compute and store GENLC.epoch8_obs_range

  Jan 6 2018: refactor so that MU is computed after SIMLIB_READ,
              and has the zHEL dependence.

  Apr 9 2019: move gen_modelPar() call down after SNHOST_DRIVER()
              so that SN->Gal transfers work for coords or redshift.

  Apr 22 2019: 
   + GENLC.GENMAG_OFF_GLOBAL += instead of =. Adds to offset
     selected in pick_NON1ASED().

  *********************************************************/

  int    NEPMIN = (int)INPUTS.CUTWIN_NEPOCH[0];
  double smear, z, z1, Tobs, Trest,  zHOST, Tobs_min, Tobs_max  ;
  int    ifilt, ifilt_obs, NEP, iep, USE_HOSTCOORD ;
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

    // check for strong lens multiple images (July 2019)
    gen_event_stronglens();

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
    if ( USE_HOSTCOORD == 0 )  { gen_MWEBV(); }

    // Get redshift of host: equals SN helio redshift, or that of wrongHost.
    zHOST = gen_zHOST(GENLC.CID, GENLC.REDSHIFT_HELIO, 
		      &GENLC.CORRECT_HOSTMATCH );
    GENLC.REDSHIFT_HOST = zHOST ;  // helio frame

    // Fetch host-galaxy using HOSTLIB (except for LCLIB model)
    // Note that SNHOST_DRIVER can change GENLC.REDSHIFT_CMB 
    // and DLMAG to match that of the HOST
    // Similarly, GENLC.REDSHIFT_HOST is changed to be the true zhost
    GEN_SNHOST_DRIVER(zHOST, GENLC.PEAKMJD); 

    // pick model params AFTER redshift/host selection (4.09.2019)
    gen_modelPar(ilc); 

    // - - - - - - - 
    // get host galaxy extinction for rest-frame models and for NON1A
    int ISREST  = ( GENFRAME_OPT   == GENFRAME_REST );
    int ISNON1A = ( INDEX_GENMODEL == MODEL_NON1ASED  ||
		    INDEX_GENMODEL == MODEL_NON1AGRID  );
    int ISMISC  = ( INDEX_GENMODEL == MODEL_SALT2   ||
		    INDEX_GENMODEL == MODEL_SIMSED  ||
		    INDEX_GENMODEL == MODEL_S11DM15 ||
		    INDEX_GENMODEL == MODEL_BYOSED 	     );

    if ( (ISREST || ISNON1A || ISMISC) && INPUTS.DO_AV ) {
      GENLC.RV = gen_RV() ;
      GENLC.AV = gen_AV() ;
    }

    override_modelPar_from_SNHOST(); // Jun 2016      

    // Now get smeared [measured] redshift and PEAKMJD
    // Smeared quantities are written to SNDATA files,
    // but are not used for anything else.

    if ( INPUTS.GENSIGMA_REDSHIFT >= 0.0 )
      { gen_zsmear( INPUTS.GENSIGMA_REDSHIFT ); }  

    // global mag offset + z-dependence 
    GENLC.GENMAG_OFF_GLOBAL += (double)INPUTS.GENMAG_OFF_GLOBAL
      + get_zvariation(GENLC.REDSHIFT_CMB,"GENMAG_OFF_GLOBAL");


  } 

  else if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID ) {
#ifdef SNGRIDGEN
    gen_GRIDevent(ilc);
    gen_modelPar(ilc) ;
#endif
    return ;
  }

  
  // --------------------------
  if ( WRFLAG_CIDRAN > 0 ) 
    { GENLC.CIDRAN = INPUTS.CIDRAN_LIST[GENLC.CID-INPUTS.CIDOFF] ; }

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
      GENLC.ISTEMPLATE[NEP]  = 1 ;
      GENLC.IEPOCH_TEMPLATE[ifilt_obs] = NEP ; 
    }

    
    // always add artificial PEAKMJD epoch to be last
    GENLC.NEPOCH++ ; NEP = GENLC.NEPOCH ;
    GENLC.MJD[NEP]       = GENLC.PEAKMJD ;
    GENLC.IFILT_OBS[NEP] = ifilt_obs ;
    GENLC.ISPEAK[NEP]    = 1 ;
    GENLC.IEPOCH_PEAK[ifilt_obs] = NEP ; 

  }  // ifilt_obs loop


  // ----------------------------------------------------------

  // pick random shift in rise & fall-times
  genshift_risefalltimes();

  // Compute epochs relative to peak
  z = GENLC.REDSHIFT_HELIO ;  z1 = 1.0 + z;
  Tobs_min = 1.0E9;  Tobs_max= -1.0E-9;
  for ( iep=1; iep <= GENLC.NEPOCH ; iep++ ) {
    Tobs = GENLC.MJD[iep] - GENLC.PEAKMJD ;
    Trest = Tobs/z1 ;
    GENLC.epoch8_obs[iep]  = Tobs  ;
    GENLC.epoch8_rest[iep] = Trest ;
    if ( Tobs < Tobs_min )  { Tobs_min = Tobs; }
    if ( Tobs > Tobs_max )  { Tobs_max = Tobs; }
  }
  GENLC.epoch8_obs_range[0] = Tobs_min ;  // global range including all bands
  GENLC.epoch8_obs_range[1] = Tobs_max ;  // --> needed for LCLIB

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

  double GAMMA_GRID_MIN = INPUTS.BIASCOR_SALT2GAMMA_GRID[0];
  double GAMMA_GRID_MAX = INPUTS.BIASCOR_SALT2GAMMA_GRID[1];
  int USE1, USE2, USE3 ;
  double DM_HOSTCOR, shape, PKMJD, RV, arg ;
  char fnam[] = "override_modelPar_from_SNHOST" ;

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
  GENLC.RV = modelPar_from_SNHOST(RV,GENLC.COLORPAR_NAME);


  if ( INDEX_GENMODEL  == MODEL_SALT2 ) {
    double a = GENLC.SALT2alpha ;
    double b = GENLC.SALT2beta  ;
    double c = GENLC.SALT2c ;

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
    if ( DM_HOSTCOR != 0.0 ) {
      GENLC.SALT2gammaDM = DM_HOSTCOR ;

      // xxxxxx temp hack until BBC7D is developed xxxxxxxxxxx
      if ( !INPUTS.OPT_DEVEL_BBC7D ) {	
	arg = -0.4*DM_HOSTCOR;
	GENLC.SALT2mB += DM_HOSTCOR;
	GENLC.SALT2x0 *= pow(TEN,arg);
      }
      // xxxxxxxxxxxxxxxxxxxxxx

    }
  }

  return ;

} // end override_modelPar_from_HOSTLIB


// *******************************************
void gen_event_stronglens(void) {

  double zSN = GENLC.REDSHIFT_CMB;
  double zLENS, hostpar[10];
  double mu_list[MXIMG_STRONGLENS];
  double angSep_list[MXIMG_STRONGLENS], phi_list[MXIMG_STRONGLENS];
  int    IDLENS, blend_flag, Nimage;
  char fnam[] = "gen_event_stronglens";

  // ------------- BEGIN ------------------

  IDLENS = get_stronglens(zSN, hostpar, &zLENS, &blend_flag,
			  &Nimage, mu_list, angSep_list, phi_list );

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
  else if ( strcmp(REJECT_STAGE,"GENPAR_SELECT") == 0 ) {
    if ( INPUTS.NGEN_LC > 0 ) { ilc-- ; }
    NGEN_REJECT.GENPAR_SELECT_FILE++ ;
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

  /*
  printf(" xxx REJECT_STAGE = '%s',  ilc=%3d -> %3d \n",
	 REJECT_STAGE, ilc_orig, ilc); fflush(stdout);
  */

  *ILC = ilc ;

} // end gen_event_reject


// *********************************
void gen_MWEBV(void) {

  // compute the following globals:
  // GENLC.MWEBV        --> 'measured' value (SFD map) passed to analysis
  // GENLC.MWEBV_ERR    --> uncertainty passed to analysis.
  // GENLC.MWEBV_SMEAR  --> true value to simulate
  //
  // Note that the logic here is different than usual;
  // the SMEAR value is the true value rather than the
  // measured value.
  //
  // May 5 2014: use FlatRan for random EBV
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
  
  int    RANFLAG, OPT, NEWCOORD ;
  double EBV, RA, DEC, MWXT_GaussRan ;
  double SigmaClip[2] = { -3.0, 3.0 } ; // 3 sigma clip
  int    LDMP = 0 ;
  char fnam[] = "gen_MWEBV";

  // -------------- BEGIN ----------------

  // always burn random number to remain synced
  // MWXT_GaussRan = GaussRan(1);  
  MWXT_GaussRan = GaussRanClip(1,SigmaClip[0],SigmaClip[1]);

  if ( INPUTS.MWEBV_FLAG  == 0 ) { return ; }  

  RA   = GENLC.RA ;
  DEC  = GENLC.DEC  ;
  RANFLAG = ( INPUTS.GENRANGE_MWEBV[0] >= 0.0 ) ; // random MWEBV

  if ( RANFLAG ) {
    // special map-generation option 
    EBV             = FlatRan(1, INPUTS.GENRANGE_MWEBV );
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

  return ;

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


  int 
    colopt, irank, istat, ALLSKIP
    ,ifilt, ifilt_obs 
    ,ifilt_rest1, ifilt_rest2, ifilt_rest3
    ;

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
      continue ;
    }

    // skip this filter if it's outside the valid rest-frame
    // wavelength range of the model (Feb 2012)
    lamz = INPUTS.LAMAVG_OBS[ifilt_obs]/(1.0+z);
    if ( lamz < GENLC.RESTLAM_MODEL[0] || lamz > GENLC.RESTLAM_MODEL[1] ) {
      GENLC.DOFILT[ifilt_obs] = 0; 
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
	NSKIP_FILTER[ifilt_obs]++ ;	
	continue ; 
      } 
    }


    // if DOFILT logical is still set, then update the
    // valid redshift range for this obs-filter.

    if ( GENLC.DOFILT[ifilt_obs] == 1 ) {
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

  return ;

} // end of gen_filtmap

// ****************************************
void genshift_risefalltimes(void) {

  double shift ;

  // ---------- BEGIN ------------


  shift = exec_GENGAUSS_ASYM(&INPUTS.GENGAUSS_RISETIME_SHIFT);
  GENLC.RISETIME_SHIFT = (float)shift ;

  shift = exec_GENGAUSS_ASYM(&INPUTS.GENGAUSS_FALLTIME_SHIFT);
  GENLC.FALLTIME_SHIFT = (float)shift ;


} // end of genshift_risefalltimes




// *****************************************
void gen_modelPar(int ilc) {

  /*********
   generate shape/luminosity parameter for model:
   stretch, delta, dm15 ... Also computes SALT2x0 and mB.

  Jan 10 2014: call get_zvariation() for SIMSED model.

  Mar 16, 2014: major overhaul to apply z-variatoin to any population
                parameter, rather than applying only to the MEAN/PEAK.

  Sep 12 2015: set fixmag

  Feb 27 2018: refactor with functions gen_shapepar_SALT2[SIMSED]

  Apr 09 2019: change function name, gen_shapepar() -> gen_modelPar().

  **********/
  int ISMODEL_SALT2     = ( INDEX_GENMODEL == MODEL_SALT2  );
  int ISMODEL_SIMSED    = ( INDEX_GENMODEL == MODEL_SIMSED );
  int ISMODEL_FIXMAG    = ( INDEX_GENMODEL == MODEL_FIXMAG );
  int ISMODEL_NON1ASED  = ( INDEX_GENMODEL == MODEL_NON1ASED );

  double ZCMB = GENLC.REDSHIFT_CMB ; // for z-dependent populations
  double shape;
  GENGAUSS_ASYM_DEF GENGAUSS_ZVAR ;
  char fnam[] = "gen_modelPar";

  //------------ BEGIN function ------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { 
    int DOFUN = ( ISMODEL_NON1ASED || ISMODEL_SIMSED);
    if ( !DOFUN) { return ; }
  }

  // ---------------------------------------
  // evaluate shape with z-dependence on population, 

  if ( !ISMODEL_SIMSED && INPUTS.NON1A_MODELFLAG<0 ) {

    char *snam = GENLC.SHAPEPAR_GENNAME ;
    GENGAUSS_ZVAR = 
      get_zvariation_GENGAUSS(ZCMB, snam, &INPUTS.GENGAUSS_SHAPEPAR);

    // pick random shape value from populatoin at this redshift
    shape  = exec_GENGAUSS_ASYM(&GENGAUSS_ZVAR) ;

    // load shape value into global GENLC struct
    GENLC.SHAPEPAR      = shape ; 
    *GENLC.ptr_SHAPEPAR = shape ; // load model-spefic shape variables 
  }


  if ( ISMODEL_SALT2 ) {
    gen_modelPar_SALT2();
  }
  else if ( ISMODEL_NON1ASED  )  {
    pick_NON1ASED(ilc, &INPUTS.NON1ASED, &GENLC.NON1ASED);
  }
  else if ( ISMODEL_SIMSED ) {
    // generate all of the SIMSED params, whatever they are.

    gen_modelPar_SIMSED();

  } // end of SIMSED if-block

  else  if ( ISMODEL_FIXMAG ) {	  
    GENLC.NOSHAPE   = FlatRan ( 2, INPUTS.FIXMAG );  
    if ( INPUTS.GENFRAME_FIXMAG == GENFRAME_REST ) 
      { GENLC.NOSHAPE += GENLC.DLMU; }
  }

  return ;

}  // end of gen_modelPar


//***************************************
void  gen_modelPar_SALT2(void) {

  // Created Feb 26 2018

  double   ZCMB = GENLC.REDSHIFT_CMB ; // for z-dependent populations
  GENGAUSS_ASYM_DEF  GENGAUSS_ZVAR ;
  char fnam[] = "gen_modelPar_SALT2";

  // ---------- BEGIN -----------

    // for SALT2, the color term is analogous to shapepar
    // so generate the 'c' and beta term here.

  GENGAUSS_ZVAR = 
    get_zvariation_GENGAUSS(ZCMB,"SALT2c",&INPUTS.GENGAUSS_SALT2c);
  GENLC.SALT2c = 
    exec_GENGAUSS_ASYM(&GENGAUSS_ZVAR) ;

  GENGAUSS_ZVAR = 
    get_zvariation_GENGAUSS(ZCMB,"SALT2ALPHA",&INPUTS.GENGAUSS_SALT2ALPHA);
  GENLC.SALT2alpha = 
    exec_GENGAUSS_ASYM(&GENGAUSS_ZVAR) ;

  GENGAUSS_ZVAR = 
    get_zvariation_GENGAUSS(ZCMB,"SALT2BETA",&INPUTS.GENGAUSS_SALT2BETA);
  GENLC.SALT2beta = 
    exec_GENGAUSS_ASYM(&GENGAUSS_ZVAR) ;


  // 2/29/2016: optional  beta(c) polynomial
  if( INPUTS.SALT2BETA_cPOLY[0] > 0.0 ) {
    double c = GENLC.SALT2c ;
    GENLC.SALT2beta += (
			INPUTS.SALT2BETA_cPOLY[0] +
			INPUTS.SALT2BETA_cPOLY[1]*c +
			INPUTS.SALT2BETA_cPOLY[2]*(c*c)  -
			INPUTS.GENGAUSS_SALT2BETA.PEAK 	 );
  }


  // now compute x0 parameter from MU and the alpha,beta params.
  GENLC.SALT2x0 = SALT2x0calc(GENLC.SALT2alpha, GENLC.SALT2beta, 
			      GENLC.SALT2x1, GENLC.SALT2c, GENLC.DLMU );
  
  GENLC.SALT2mB = SALT2mBcalc(GENLC.SALT2x0);

  return;

} // end gen_modelPar_SALT2

//***************************************
void  gen_modelPar_SIMSED(void) {

  // Created Feb 26 2018
  // Move code from gen_modelPar() to here, and add code to
  // pick correlated Guass randoms among arbitrary number of
  // variables.
  //
  // Dec 20 2018: set ISIMSED_SEQUENTIAL instead of PARVAL[0]
  //     (part of refactor for SIMSED loops)
  //
  // Apr 28 2019: return for GRIDGEN, after computing DLMU

  int     NPAR      = INPUTS.NPAR_SIMSED;
  int     NROW_COV  = INPUTS.NROW_SIMSED_COV;
  double  ZCMB      = GENLC.REDSHIFT_CMB ; // for z-dependent populations

  int  ipar, ipar1, ipar_model, genflag, opt_interp, opt_gridonly ;
  int  irow, irow1, irow_COV, NRANGEN_ITER=0 ;
  double ARG, ranGauss, parVal, parVal_old ;
  double tmpRan, tmpMat, PEAK, PMIN, PMAX, SIGMA ;
  double GAURAN[MXPAR_SIMSED], CORRVAL[MXPAR_SIMSED] ;
  GENGAUSS_ASYM_DEF GENGAUSS_ZVAR ;
  char *parName;
  char fnam[] = "gen_modelPar_SIMSED" ;
  int LDMP_COV = 0 ;

  // ----------- BEGIN ------------

  // use SALT2x0 parameter for SIMSED ... it just converts
    // MU into a flux-scale.
  ARG = -0.4 * GENLC.DLMU ;
  GENLC.SALT2x0 = pow(TEN , ARG );

#ifdef SNGRIDGEN
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
      GAURAN[irow]  = GaussRan(1);   // Gauss Random
      CORRVAL[irow] = 0.0 ;
    }
    for(irow=0; irow < NROW_COV; irow++ ) {
      ipar = INPUTS.IPARLIST_SIMSED_COV[irow];
      PEAK = INPUTS.GENGAUSS_SIMSED[ipar].PEAK;
      PMIN = INPUTS.GENGAUSS_SIMSED[ipar].RANGE[0];
      PMAX = INPUTS.GENGAUSS_SIMSED[ipar].RANGE[1];
      CORRVAL[irow] = PEAK ;
      for(irow1=0; irow1 < NROW_COV; irow1++ ) {
	tmpMat = INPUTS.CHOLESKY_SIMSED_COV[irow1][irow] ;
	tmpRan = GAURAN[irow1] ;
	CORRVAL[irow] += ( tmpMat * tmpRan) ;
      }
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
    opt_interp   = ( genflag & OPTMASK_SIMSED_PARAM    ) ;
    opt_gridonly = ( genflag & OPTMASK_SIMSED_GRIDONLY ) ;
    parName      = INPUTS.PARNAME_SIMSED[ipar] ;
    irow_COV     = INPUTS.IROWLIST_SIMSED_COV[ipar] ;
    GENLC.SIMSED_PARVAL[ipar]  = -9.0 ;

    // skip baggage params
    if ( (genflag & OPTMASK_SIMSED_param)  != 0 ) { continue ; }

    if ( opt_interp ) {

      if ( irow_COV >= 0 ) {
	// correlated random, Gaussian only
	parVal = CORRVAL[irow_COV] ; 
      }
      else {
	// uncorrelated random, more general profile
	GENGAUSS_ZVAR = 
	  get_zvariation_GENGAUSS(ZCMB,parName,
				  &INPUTS.GENGAUSS_SIMSED[ipar]);	
	parVal = exec_GENGAUSS_ASYM( &GENGAUSS_ZVAR );
      }

      if ( opt_gridonly  ) {
	parVal_old = parVal ;
	parVal = nearest_gridval_SIMSED(ipar_model,parVal_old);
      }
    } // opt_interp
    
    //    printf(" xxx %s: PARVAL[%d] = %f \n", fnam, ipar, parVal);
    GENLC.SIMSED_PARVAL[ipar]  = parVal ;

  }  // end ipar

  // ------------------
  // set PARVAL[0] to the SED index ... the IGEN index.
  if ( INPUTS.OPTMASK_SIMSED == OPTMASK_SIMSED_GRIDONLY ) 
  { ISIMSED_SEQUENTIAL  = (double)NGENLC_TOT ; } // IGEN = ISED

  
  return;

} // end gen_modelPar_SIMSED



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

#ifdef SNGRIDGEN
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
    // xxx mark delete 5/2019 init_genmag_NON1ASED(GENLC.TEMPLATE_INDEX,sedFile); 

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


#ifdef SNGRIDGEN
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
void genran_obsNoise(void) {

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

  int ep, ifilt, ifilt_obs, ifield;
  char fnam[] = "genran_obsNoise" ;

  // -------------- BEGIN --------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return ; }

  // one random per epoch
  for ( ep = 1; ep <= GENLC.NEPOCH; ep++ )  {  

    ifilt_obs = GENLC.IFILT_OBS[ep] ;    
    GENLC.RANGauss_NOISE_SEARCH[ep] = -99999. ;  
    GENLC.RANGauss_NOISE_ZP[ep]     = -99999. ;  

    // skip un-used epochs so that randoms stay synced with
    // previous (10_33g) snana version.
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }
    if ( GENLC.ISPEAK[ep]             ) { continue ; }
    if ( GENLC.ISTEMPLATE[ep]         ) { continue ; }
    
    // load random into global
    GENLC.RANGauss_NOISE_SEARCH[ep] = GaussRan(1) ;  // from list 1
    GENLC.RANGauss_NOISE_ZP[ep]     = GaussRan(1) ; 
  }


  // one random per filter (for template noise) and field overlap

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_SIMLIB; ifilt++ ) {

    ifilt_obs =  GENLC.IFILTMAP_SIMLIB[ifilt] ;

    // init to crazy values
    for(ifield=0; ifield < MXFIELD_OVP_SIMLIB; ifield++ ) {      
      GENLC.RANGauss_NOISE_TEMPLATE[ifield][ifilt_obs] = -99999. ; 
    }

    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }
   
    // Aug 27 2014: the extra field loop here changes random sequence.
    for(ifield=0; ifield < MXFIELD_OVP_SIMLIB; ifield++ ) {      
      GENLC.RANGauss_NOISE_TEMPLATE[ifield][ifilt_obs] = GaussRan(1) ; 
    } 
    
  }

}   // end of genran_obsNoise


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
  // May 2, 2009: set smearing for all models, rest or observer
  //              Fix bug: <=MXFILTINDX should be <MXFILTINDX
  //
  // Jul 24, 2009: use GENSMEAR_RAN_FILTER[0] = global
  //                   GENSMEAR_RAN_FILTER[ifilt] per filter
  //
  // March 23, 2010: apply sigma-clipping to avoid crazy
  //                 mag-variations
  // 
  // May 26, 2010:  fix dumb bug by setting GENLC.GENSMEAR_RAN_FILTER[0] = 
  //                clipped rr8 instead of gran().
  //
  // Jul 1, 2010: rename gen_lumismear -> genran_modelSmear
  //              (generate randoms for model smearing)
  //
  // Jul 27, 2011: set GENLC.COVMAT_SCATTER_GRAN and 
  //
  // Mar 12, 2012: call intrinsic_magSmear_SETRAN
  //
  // Feb 26 2013: 
  //   For GENSMEAR_RAN use separate RANLIST so that different randoms
  //   can be used without changing the main random sequence.
  //   See user input RANLIST_START_GENSMEAR so that user can 
  //   vary random intrinsic scatter with the same SN.
  //
  // May 19, 2013: generate separate list for FILTERs and MODEL
  //               so that they are not correlated
  //
  // Jan 21 2014: update to handle flat-randoms as well as Gauss randoms.
  //
  // Mar 24 2016: for GENGRID, return after all randoms are initalized to zero.

  int ifilt ;
  int    ILIST_RAN = 2 ; // list to use for genSmear randoms
  double rr8, rho, RHO, rmax, rmin, rtot, RANFIX ;
  char fnam[] = "genran_modelSmear" ;

  // -------------- BEGIN --------

  // init all smearing to zero

  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ )  {  
    GENLC.GENSMEAR_RANGauss_FILTER[ifilt] = 0.0 ;  
    GENLC.GENSMEAR_RANGauss_MODEL[ifilt]  = 0.0 ;  
    GENLC.GENSMEAR_RANFlat_MODEL[ifilt]   = 0.0 ;  
  }

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return ; }

  // always generate randoms to stay synced, even if mag smear is zero.

  NSTORE_RAN[ILIST_RAN] = INPUTS.RANLIST_START_GENSMEAR ; // reset

  rmin = INPUTS.SIGMACLIP_MAGSMEAR[0] ;
  rmax = INPUTS.SIGMACLIP_MAGSMEAR[1] ;


  // for global coherent smearing
  
  GENLC.GENSMEAR_RANGauss_FILTER[0] = GaussRanClip(1,rmin,rmax);

  // Jun 26 2019: check option for asymmetric smear
  double siglo = (double)INPUTS.GENMAG_SMEAR[0] ;
  double sighi = (double)INPUTS.GENMAG_SMEAR[1] ;
  if ( fabs(siglo-sighi) > 1.0E-6 ) {
    sighi /= siglo ;      siglo /= siglo ;
    GENLC.GENSMEAR_RANGauss_FILTER[0] = biGaussRan(siglo,sighi);
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
    rr8 = GaussRanClip(ILIST_RAN,rmin,rmax);
    rtot =  rho * GENLC.GENSMEAR_RANGauss_FILTER[0] +  RHO * rr8 ; 
    GENLC.GENSMEAR_RANGauss_FILTER[ifilt]  = rtot ;      
  }



  // generate Guassian randoms for intrinsic scatter [genSmear] model
  for ( ifilt=1; ifilt < MXFILTINDX; ifilt++ ) {
    GENLC.GENSMEAR_RANGauss_MODEL[ifilt] = GaussRanClip(ILIST_RAN,rmin,rmax);
  }

  // repeat for 0-1 [flat] randoms (Jan 2014 - changes random sync)
  for ( ifilt=1; ifilt < MXFILTINDX; ifilt++ ) 
    {  GENLC.GENSMEAR_RANFlat_MODEL[ifilt]  = FlatRan1(ILIST_RAN);  }


  // check option to fix genSmear randoms for debugging
  RANFIX = (double)INPUTS.GENSMEAR_RANGauss_FIX ;  
  if ( RANFIX > -99.0 ) {
    printf(" XXXX DEBUG FLAG: "
	   "Fix all GENSMEAR Gauss-randoms to %6.3f  for  CID=%d\n", 
	   RANFIX, GENLC.CID ) ;
    for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
      GENLC.GENSMEAR_RANGauss_FILTER[ifilt] = RANFIX ;
      GENLC.GENSMEAR_RANGauss_MODEL[ifilt]  = RANFIX ;
    }
  }


  // set randoms for instrinsic scatter matrix (July 2011)
  GENLC.COVMAT_SCATTER_GRAN[0] =  GaussRan(1);
  GENLC.COVMAT_SCATTER_GRAN[1] =  GaussRan(1);
  GENLC.COVMAT_SCATTER_GRAN[2] =  GaussRan(1);
  GEN_COVMAT_SCATTER ( GENLC.COVMAT_SCATTER_GRAN,    // <== input 
		       GENLC.COVMAT_SCATTER );       // <== output
  

  // -----------------------------------
  // pass randoms to genSmear function
  int NRAN ;
  NRAN = GENLC.NRANGauss_GENSMEAR ;
  SETRANGauss_genSmear(NRAN, &GENLC.GENSMEAR_RANGauss_MODEL[1] );

  NRAN = GENLC.NRANFlat_GENSMEAR ;
  SETRANFlat_genSmear(NRAN, &GENLC.GENSMEAR_RANFlat_MODEL[1] );

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

  char 
     cmd[MXPATHLEN]
    ,cfilt[4]
    ,filtFile[MXPATHLEN]
    ,filtName[40]   // full name returned from get_filttrans
    ;

  double 
    magprim8
    ,lam8[MXLAMSIM]
    ,TransSN8[MXLAMSIM]
    ,TransREF8[MXLAMSIM] 
    ;
  
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
    magprim8 = lam8[0] = TransSN8[0] = TransREF8[0] = 0.0 ; 
    NLAM = 0 ;

    // note that both the SN and REF trans are returned,
    // but for simulations they are always the same.

    get_filttrans__(&MASKFRAME, &ifilt_obs, filtName, &magprim8, &NLAM, 
		    lam8, TransSN8, TransREF8, strlen(filtName) );

    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );
    sprintf(filtFile,"%s/%s.dat", PATH_FILTERS, cfilt );
   
    fp_filt = fopen(filtFile, "wt") ;

    for ( ilam=0; ilam < NLAM; ilam++ ) {
      fprintf(fp_filt,"%7.1f  %.6f \n", lam8[ilam], TransSN8[ilam] );
    }

    fclose(fp_filt);

  }  // end of ifilt loop


} // end of wr_SIMGEN_FILTERS

// ***********************************
void wr_SIMGEN_DUMP(int OPT_DUMP, SIMFILE_AUX_DEF *SIMFILE_AUX) {

  /***
   Created Feb 20, 2009 by R.Kessler
   Major re-write 5/09/2009

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

  Feb 13 2014: case insensitive string checks.

  Jul 30 2014; allow string ; see *str

  Aug 19 2014: pass SIMFILE_AUX struct.

  Feb 05, 2015: for r*8, write specific formats for MJD, RA, DEC.

  Mar 16, 2017: write each line to strig SIMFILE_AUX->OUTLINE,
                then make single write per line instead of per value.

  Aug 14 2017: check prescale, and write a little extra info at top of file.

  Dec 01 2017: add hash (#) in front of all comments.

  May 14 2019: free(SIMFILE_AUX->OUTLINE)

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
      { PREP_SIMGEN_DUMP(0); }    // list variables, then quit.
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
      // xxx mark deletesprintf(SIMFILE_AUX->OUTLINE,"%s %s", 
      // xxx SIMFILE_AUX->OUTLINE, pvar ) ;
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
    sprintf(c1err,"Undefined SIMGEN_DUMP variable: '%s'", varName);
    errmsg(SEV_WARN, 0, fnam, c1err, ""); 
    madend(0);  //  give ugly face, but do NOT abort
    PREP_SIMGEN_DUMP(0); // print list of valid varnames, then quit
  }
  
  if ( NMATCH > 1 ) {
    sprintf(c1err,"SIMGEN_DUMP variable '%s'", varName);
    sprintf(c2err,"defined %d times ??", NMATCH);
    errmsg(SEV_WARN, 0, fnam, c1err, c2err ); 
    madend(0);  //  give ugly face, but do NOT abort
    PREP_SIMGEN_DUMP(0); // print list of valid varnames, then quit
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
  // OPT_DUMP = 0 => dump variable list then quit.
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

  int i, ifilt, ifilt_obs, ifilt_rest, ipar, imap, ivar, ivar2 ;
  int NTMP, FOUND ;
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


  // ---------
  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"FIELD");  // warning; does not properly treat field-overlap
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRCHAR = GENLC.FIELDNAME[1] ;
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
  sprintf(cptr,"LENSDMU") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.LENSDMU ;
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
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.SNSEP ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALSNDLR") ;   // 2/2019: DLR from Sako 2014, Gupta 2016
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.DLR ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"GALSNDDLR") ;   //2/2019:  d_DLR = SNSEP/DLR
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.DDLR ;
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

    /* xxxx mark delete July 12 2019 xxxxxxxxxx
    // avoid declaring a variable twice.
    cptr = HOSTLIB_OUTVAR_EXTRA.NAME[ivar] ;     FOUND = 0 ;
    for(ivar2=0; ivar2<NTMP; ivar2++ ) 
      { if ( strcmp(SIMGEN_DUMP[ivar2].VARNAME,cptr) == 0 ) { FOUND=1; }    }
    if ( FOUND ) { continue ; }
    xxxxxxxx end mark delete xxxxxxxxxxx*/

    if ( HOSTLIB_OUTVAR_EXTRA.USED_IN_WGTMAP[ivar] ) { continue; } 
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"%s", HOSTLIB_OUTVAR_EXTRA.NAME[ivar] );
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &HOSTLIB_OUTVAR_EXTRA.VALUE[ivar];
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
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &SNHOSTGAL.SB_FLUX[ifilt_obs];
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
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.peakmag8_obs[ifilt_obs] ;
    NVAR_SIMGEN_DUMP++ ;
  }

  // Aug 2017: allow traditional PEAKMAG_x for true peakmag (not observed)
  for ( ifilt=0; ifilt < INPUTS.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt]; 
    cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
    sprintf(cptr,"PEAKMAG_%c", FILTERSTRING[ifilt_obs] ) ;
    SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.peakmag8_obs[ifilt_obs] ;
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
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8=&GENLC.peakmag8_rest[ifilt_rest];
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


  if ( INDEX_GENMODEL == MODEL_BYOSED ) {
    for( ipar=0 ;  ipar < Event_BYOSED.NPAR; ipar++ ) {
      cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
      sprintf(cptr, "%s", Event_BYOSED.PARNAME[ipar] );
      SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &Event_BYOSED.PARVAL[ipar]; 
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
  sprintf(cptr,"AV") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.AV ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"RV") ;
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

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"S2alpha") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2alpha ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SALT2alpha") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2alpha ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"S2beta") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2beta ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SALT2beta") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2beta ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"S2x0") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2x0 ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SALT2x0") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2x0 ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"S2x1") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2x1 ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SALT2x1") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2x1 ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"S2c") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2c ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SALT2c") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2c ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"S2mb") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2mB ;
  NVAR_SIMGEN_DUMP++ ;

  cptr = SIMGEN_DUMP[NVAR_SIMGEN_DUMP].VARNAME ;
  sprintf(cptr,"SALT2mb") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.SALT2mB ;
  NVAR_SIMGEN_DUMP++ ;

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
  sprintf(cptr,"MAGSMEAR_COH") ;
  SIMGEN_DUMP[NVAR_SIMGEN_DUMP].PTRVAL8 = &GENLC.MAGSMEAR_COH ;
  NVAR_SIMGEN_DUMP++ ;

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
    happyend();
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

  int i, imjd, IDSPEC, indx, OPT_TEXPOSE, NOPT_SNR=0;
  int idump_inp, idump_list;
  int NVAR_SIMGEN_DUMP_START = NVAR_SIMGEN_DUMP ;
  int NVAR_SIMGEN_DUMP_END   = NVAR_SIMGEN_DUMP ;
  char PREFIX[12], *cptr ;
  char fnam[] = "PREP_SIMGEN_DUMP_TAKE_SPECTRUM" ;

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
	   ,GENLC.epoch8_rest[i] 
	   ,GENLC.genmag8_rest[i],  FILTERSTRING[ifilt_rest1]
	   ,GENLC.genmag8_rest2[i], FILTERSTRING[ifilt_rest2]
	   ,GENLC.genmag8_obs[i],   FILTERSTRING[ifilt_obs]
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
  //   float rangen -> double FlatRan  
  //   float GENRANGE_PEAKMJD -> double GENRANGE_PEAKMJD
  //   float PKMJD, MJD[2] -> double
  //
  // Sep 17 2015: check option to snap to nearest PEAKMJD read from file.
  // Jan 25 2017: if last PEAKMJD was read from SIMLIB header, then
  //              just return current PEAKMJD in case SIMLIB_NREPEAT
  //              key is set.
  //

  double PKMJD, MJD[2];
  int    NSKIP_RANGE, NSKIP_MJD, i ;
  char   fnam[] = "gen_peakmjd" ;

  // ------------- BEGIN --------------

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

  PKMJD = FlatRan (1,INPUTS.GENRANGE_PEAKMJD );


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
  char fnam[] = "gen_peakmjd_smear" ;

  // ----------------- BEGIN -----------------

  // always burn Gaussian random, regardless of option.
  GENLC.PEAKMJD_RANGauss = GaussRan(1); 

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
  double zcmb ; 
  //  char fnam[] = "gen_redshift_cmb" ;

  // --------------- BEGIN ------------
  // check for delta-function in redshift
  if ( INPUTS.GENRANGE_REDSHIFT[0] == INPUTS.GENRANGE_REDSHIFT[1] ) {
    zcmb = INPUTS.GENRANGE_REDSHIFT[0] ; 
  }
  else {

    // pick from pure Hubble distribution with reweight based on rate model
    zcmb = genz_hubble(zmin,zmax, GENLC.RATEPAR ) ;  
  }

  return(zcmb);

}  // end of gen_redsfhit_cmb



// *******************************************
double gen_redshift_helio(void) {
  //
  // Created Jan 2018 by R.Kessler
  // Convert zcmb to zhelio and apply peculiar velocity (vpec in km/s)
  // Function returns z_helio, and also stores GENLC.VPEC
  
  double zCMB = GENLC.REDSHIFT_CMB ;
  double RA   = GENLC.RA;
  double DEC  = GENLC.DEC ;
  double vpec, zhelio, dzpec ;
  char fnam[] = "gen_redshift_helio" ;

  // ----------- BEGIN ------------

  // check (v10_31) legacy option to keep zhelio = zcmb
  if ( INPUTS.VEL_CMBAPEX == 0.0 ) { return zCMB ; }

  zhelio = zhelio_zcmb_translator(zCMB, RA,DEC, "eq", -1);

  // apply v_pec
  vpec = (double)INPUTS.GENSIGMA_VPEC * GaussRan(2) ;
  GENLC.VPEC = vpec; 
  if ( vpec != 0.0 ) {
    dzpec = vpec/LIGHT_km ;
    zhelio += dzpec ;
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
  char fnam[] = "gen_redshift_LCLIB" ;

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

  granz = GaussRanClip(1, (double)-3.0, (double)+3.0) ;
  ZERR  = INPUTS.GENSIGMA_REDSHIFT ;
  SNHOSTGAL.ZSPEC             = ZHEL_TRUE ;
  SNHOSTGAL.ZSPEC_ERR         = ZERR;
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
      { GENLC.REDSHIFT_RAN[i] = GaussRan(1); }
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
    sprintf(c1err,"Could not generate random Z after %d tries", NZRAN);
    sprintf(c2err,"Each Z < ZGEN_MIN = %f  z[cmb,hel]=%f,%f)", 
	    ZGEN_MIN, GENLC.REDSHIFT_CMB, GENLC.REDSHIFT_HELIO );
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

  if ( INPUTS.VEL_CMBAPEX == 0.0 ) 
    { GENLC.REDSHIFT_CMB_SMEAR = zhelio ; } // legacy option, ZCMB = ZHELIO
  else
    { GENLC.REDSHIFT_CMB_SMEAR = 
	zhelio_zcmb_translator(zhelio,RA,DEC, "eq", +1); 
    }

 
  // --------------------------------------------
  // Determine VPEC correction using Gaussian-random number.
  // Note that correction has oppoisite sign true value.
  // If VPEC_ERR is >= to true scatter, this is a flag
  // to NOT apply a correction and set measured VPEC_SMEAR=0.

  double GAURAN_VPEC = GENLC.REDSHIFT_RAN[MXZRAN-1];
  if ( INPUTS.VPEC_ERR < INPUTS.GENSIGMA_VPEC )
    { GENLC.VPEC_SMEAR = -GENLC.VPEC + (INPUTS.VPEC_ERR * GAURAN_VPEC) ; }
  else
    { GENLC.VPEC_SMEAR = 0.0 ; }  // do NOT apply correction (e.g., high-z)

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

  double mu, lensDMU, ran1 ;
  char fnam[] = "gen_distanceMag" ;

  // -------------- BEGIN ------------


  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID ) 
    { ran1 = 0.5; } // no randoms for GRID mode
  else
    { ran1 = FlatRan1(2); }  // normal gen: always burn random: May 7 2017

  if ( IGNOREFILE(INPUTS.WEAKLENS_PROBMAP_FILE) ) 
    { lensDMU = 0.0 ; }
  else 
    { lensDMU = gen_lensDMU(zCMB,ran1);  }


  if ( INPUTS.WEAKLENS_DSIGMADZ > 1.0E-8 ) {
    lensDMU = zCMB * INPUTS.WEAKLENS_DSIGMADZ * GaussRan(1) ;
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
  double H0 = INPUTS.H0 / ( 1.0E6 * PC_km) ;
  double OM = INPUTS.OMEGA_MATTER ;
  double OL = INPUTS.OMEGA_LAMBDA ;
  double w0 = INPUTS.W0_LAMBDA ;

  if ( zCMB <= 1.0E-10 ) { return 0.0 ; }  // avoid inf, May 2013
  mu = dLmag (H0,OM,OL,w0, zCMB, zHEL);

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

   May 5, 2014: float rangen -> double FlatRan
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

  *****************/

  double z, zran, z_atmax, dz, H0, OM, OL, W0, w, wgt ;
  double wran1, zrange[2] ;
  int iz, NZ, ISFLAT, ISPOLY, ilist, NEWZRANGE, FIRST ;
  char fnam[] = "genz_hubble" ;

  // --------------- BEGIN ------------

  if ( zmin <= 1.0E-9 || zmax <= 1.0E-9 ) {
    sprintf(c1err,"Invalid zmin,zmax = %le, %le at CID=%d LIBID=%d", 
	    zmin, zmax, GENLC.CID, GENLC.SIMLIB_ID );
    sprintf(c2err,"Both must be > %le (10pc)", ZAT10PC );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  zrange[0] = zmin ;
  zrange[1] = zmax ;

  H0   = INPUTS.H0 / ( 1.0E6 * PC_km) ;
  OM   = INPUTS.OMEGA_MATTER ;
  OL   = INPUTS.OMEGA_LAMBDA ;
  W0   = INPUTS.W0_LAMBDA ;   

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
	w = polyEval(4, RATEPAR->MODEL_PARLIST[1], z);
      }
      else {
	w = dVdz ( H0, OM, OL, W0, z );
	w /= (1.0+z);          // Dec 2006: time dilation factor
	w *= genz_wgt(z,RATEPAR) ;
      }

      if ( w > RATEPAR->ZGENWGT_MAX ) {
	RATEPAR->ZGENWGT_MAX = w ;  
	z_atmax     = z;
      }

    } // end of iz loop

    if ( FIRST ) {
      printf("  Found Max dN/dz * wgt = %e at z = %8.3f \n", 
	     RATEPAR->ZGENWGT_MAX, z_atmax); fflush(stdout);
    }

  }  // end of ZGENWGT_MAX if-block


  // =======================================

 PICKRAN:

  ilist = 1 ;  // select which random list

  // pick random redshift
  zran   = FlatRan ( ilist, zrange ); // random redshift
  wran1  = FlatRan1( ilist );  // for weight

  if ( zran < zmin || zran > zmax ) { goto PICKRAN ; }

  // compute wgt = dN/dz for this z
  if ( ISFLAT )  { 
    // flat redshift distribution
    wgt = 1.0 ; 
  }
  else if ( ISPOLY ) {
    // user-specified polynomial function of redshift
    wgt = polyEval(4, RATEPAR->MODEL_PARLIST[1], zran) / RATEPAR->ZGENWGT_MAX ;
  }
  else if ( INPUTS.USE_SIMLIB_DISTANCE ) {
    wgt = 1.0 ;
  }
  else {
    // physical distribution
    w    = dVdz ( H0, OM, OL, W0, zran );
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
  wpoly = polyEval(4, RATEPAR->DNDZ_ZPOLY_REWGT, z) ;
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
  //  char fnam[] = "init_RATPAR" ;

  // ------------- BEGIN -------------

  RATEPAR->DNDZ_ZEXP_REWGT     = 0.0 ;
  RATEPAR->DNDZ_SCALE[0]       = 1.0 ; // Apr 19 2017
  RATEPAR->DNDZ_SCALE[1]       = 1.0 ; // Apr 19 2017
  RATEPAR->DNDZ_ALLSCALE       = 1.0 ; // Aug 30 2017
  RATEPAR->RATEMAX = 0.0 ;

  RATEPAR->DNDZ_ZPOLY_REWGT[0] = 1.0 ;
  for(i=1; i<4; i++ ) { RATEPAR->DNDZ_ZPOLY_REWGT[i] = 0.0 ;  }

  sprintf(RATEPAR->NAME, "HUBBLE"); 
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
  double H0, OM, OL, W0 ;
  int NBZ, iz;
  char fnam[] = "SNcount_model" ;

  // ------------- BEGIN ------------

  SNsum = 0.0 ;

  H0   = INPUTS.H0 ; 
  OM   = INPUTS.OMEGA_MATTER ;
  OL   = INPUTS.OMEGA_LAMBDA ;
  W0   = INPUTS.W0_LAMBDA ;   

  NBZ = (int)( (zMAX - zMIN ) * 1000. ) ;
  if ( NBZ < 10 ) { NBZ = 10; }
  dz  = ( zMAX - zMIN ) / (double)NBZ ;

  for ( iz=1; iz <= NBZ; iz++ ) {
    ztmp   = zMIN + dz * ((double)iz - 0.5 ) ;
    vtmp   = dVdz ( H0, OM, OL, W0, ztmp );
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
  ***/

  double sfr, sfrint, h, H0, OM, OL, W, z1, rate ;
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
    H0=INPUTS.H0 ;
    OM=INPUTS.OMEGA_MATTER; OL=INPUTS.OMEGA_LAMBDA ;  W=INPUTS.W0_LAMBDA;

    // get star formation rate    
    sfr    = SFRfun(H0,z);

    // determine integrated stellar mass 
    sfrint = SFR_integral ( H0, OM, OL, W, z ) ;

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
    // put something here to avoid abort; this isn't really a rate
    rate = 3.0E-5 * polyEval(4, RATEPAR->MODEL_PARLIST[1], z);
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

  double Rate=0.0, cosb, bb;
  double BPOW[MXPOLY_GALRATE+1], COSBPOW[MXPOLY_GALRATE+1] ;
  int i;
  char fnam[] = "GALrate_model" ;

  // -------------- BEGIN ---------------

  if ( strcmp(RATEPAR->NAME,"COSBPOLY") == 0 ) {
    cosb = cos(b*RADIAN);
    COSBPOW[0] = 1.0 ;
    for(i=0;  i <= MXPOLY_GALRATE; i++ ) {
      if( i > 0 ) { COSBPOW[i] = COSBPOW[i-1] * cosb ; }
      Rate += ( COSBPOW[i] * RATEPAR->MODEL_PARLIST[1][i] ) ;
    }
  }
  else if ( strcmp(RATEPAR->NAME,"BPOLY") == 0 ) {
    bb      = fabs(b);
    BPOW[0] = 1.0 ;    
    for(i=0;  i <= MXPOLY_GALRATE; i++ ) {
      if( i > 0 ) { BPOW[i] = BPOW[i-1] * bb ; }
      Rate += ( BPOW[i] * RATEPAR->MODEL_PARLIST[1][i] ) ;
    }
  }

  else {
    sprintf(c1err,"Unknown GALrate model: '%s' ", RATEPAR->NAME );
    sprintf(c2err,"Check GENMODEL key");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  return(Rate);

} // end GALrate_model

// ***********************************
double gen_AV(void) {

  // select AV from exponential distribution,
  // dN/dAv = exp(-av/tau) + exp(.5*av^2/sig^2)
  //
  // Oct 26, 2012: tau > 50 -> flat distribution
  // Mar 15, 2013: include Gaussian core (requested by Rodney for HST)
  //                      float -> double
  // 
  // Feb 02, 2016: for pure Gaussian, repeat until avmin < AV < avmax;
  //               and abort after MAXTRY_ABORT tries
  //
  // Feb 19, 2016: added missing 'goto DONE' in DOFUN_GAUSS (D.Jones)
  //
  // Feb 24, 2016: adjusted gaussian to allow non-0 centroid using
  //               the GENGAUPEAK_AV keyword (D. Jones)
  //
  // Jan 14 2017: 
  //   while fixing warnings from -Wall, found and fixed bug setting
  //   un-initialized "peak":   peakGauss = INPUTS.GENGAUPEAK_AV ;
  //
  double  tau, sig, ratio, peakGauss, expmin, expmax, expdif ;
  double avmin, avmax, AV, ran_EXPON, ran_GAUSS, ran_WGT ;
  int DOFUN_EXPON, DOFUN_GAUSS ;
  char fnam[] = "gen_AV" ;

  // ------------ BEGIN -------------

  AV    = 0.0 ;
  DOFUN_EXPON = DOFUN_GAUSS = 0 ;

  // always burn randoms to stay synced.
  ran_EXPON = FlatRan1(1) ;
  ran_GAUSS = GaussRan(1);    
  ran_WGT   = FlatRan1(1) ;  

  if ( INPUTS.WV07_GENAV_FLAG )  { AV = GENAV_WV07(); return(AV); }

  // pick from exponential
  GENLC.AVTAU = INPUTS.GENEXPTAU_AV 
    + get_zvariation(GENLC.REDSHIFT_CMB,"GENEXPTAU_AV");

  GENLC.AVSIG = INPUTS.GENGAUSIG_AV 
    + get_zvariation(GENLC.REDSHIFT_CMB,"GENEXPSIG_AV");

  peakGauss = INPUTS.GENGAUPEAK_AV; 

  GENLC.AV0RATIO = INPUTS.GENRATIO_AV0
    + get_zvariation(GENLC.REDSHIFT_CMB,"GENRATIO_AV0" );

  tau   = GENLC.AVTAU ;
  sig   = GENLC.AVSIG ;
  ratio = GENLC.AV0RATIO ;  // Gauss/expon ratio

  avmin = INPUTS.GENRANGE_AV[0] ;
  avmax = INPUTS.GENRANGE_AV[1] ;

  // sanity check
  if ( (avmax > avmin) && (tau==0.0 && sig==0.0) ) {
    sprintf(c1err,"GENRANGE_AV = %.3f to %.3f", avmin, avmax);
    sprintf(c2err,"but GENEXPTAU_AV=0 and GENGAUSIG_AV=0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // check for trivial cases

  if ( tau   <= 0.001 && sig <= 0.001 )  // delta function at AVMIN
    { AV = avmin; goto DONE ; }

  if ( avmin == avmax )       // delta-function at AVMIN
    { AV = avmin; goto DONE ; }  

  if ( tau > 50. && ratio == 0.0 )    // flat distribution
    { AV = avmin + (avmax - avmin) * ran_EXPON ; goto DONE ; }


  // check for pure exponential or pure Gaussian
  if ( tau > 0.0 &&  sig <= 0.0 ) { DOFUN_EXPON = 1 ; }
  if ( sig > 0.0 &&  tau <= 0.0 ) { DOFUN_GAUSS = 1 ; }


  // check for mixed.
  if ( tau > 0.0 && sig > 0.0 && ratio > 0.0 ) {
    double WGT_EXPON, WGT_GAUSS, WGT_SUM ;
    WGT_EXPON = 1.0 / tau ;
    WGT_GAUSS = ratio / sqrt(TWOPI * sig);
    WGT_SUM = WGT_EXPON + WGT_GAUSS ;
    if ( ran_WGT < WGT_EXPON/WGT_SUM ) 
      { DOFUN_EXPON = 1; }
    else
      { DOFUN_GAUSS = 1; }
  }

  // pure exponential 
  if ( DOFUN_EXPON ) {
    expmin = expf(-avmin/tau) ;  // note that expmin > expmax !!!
    expmax = expf(-avmax/tau) ;
    expdif = expmin - expmax ;    
    AV = -tau * log( expmin - expdif*ran_EXPON ) ;
    goto DONE ;
  }

  // pure Guassian (AV > 0 only)                   
  if ( DOFUN_GAUSS ) {
    int MAXTRY_ABORT = 10000  ;  // abort after this many tries   
    int itry = 0 ;
    AV = -9999.0 ;
    while ( AV < avmin || AV > avmax ) {
      if ( itry > MAXTRY_ABORT ) {
	sprintf(c1err,"Can't find Gauss-AV between %.2f and %.2f",
                avmin, avmax);
        sprintf(c2err,"after %d tries (sigma=%.2f)\n", itry, sig );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }
      if ( itry  > 1 ) { ran_GAUSS = GaussRan(1); }
      if (peakGauss > 0) 
	{ AV = sig * ran_GAUSS + peakGauss; }
      else
	{ AV = sig * fabs(ran_GAUSS); }
 
      itry++ ;
    }
    goto DONE ;
  }

  // if we get here then abort on confusion.
  sprintf(c1err,"Could not determine AV from Expon. or Gaussian ??");
  sprintf(c2err,"tau=%f  sig=%f  ratio=%f", tau, sig, ratio);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;

 
 DONE:
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
  BEXP = 1./sqrtf(sqsigma * 2. * 3.14159) ;

  if ( REWGT_AEXP > -1.0E-9 ) { AEXP *= REWGT_AEXP; } // April 2018

  W0 = AEXP + BEXP ; // weight at AV=0

  // pick random AV on defined interval

 PICKAV:
  AV = FlatRan ( 1 ,INPUTS.GENRANGE_AV );

  // compute relative wgt
  arg_A = -AV/tau ;                  // broad exponential
  arg_B = -0.5 * AV*AV / sqsigma ;   // sharp half-Gaussian core

  W = AEXP * exp(arg_A) + BEXP * exp(arg_B)  ;
  W /= W0;

  if ( W < FlatRan1(1) ) { goto PICKAV ; }

  return(AV) ;

}  // end of GENAV_WV07

// ***********************************
double gen_RV(void) {

  // Aug 2006: select RV from bifurcated gaussian
  // Mar 2008: after 1000 tries, abort
  // May 03, 2012: simplify further using exec_GENGAUSS_ASYM
  // Jan 28, 2014: call SNPARVAL_from_SNHOST() function.
  // Jun 11, 2016: move SNPARVAL_from_SNOST call to override_SNPARVAL_xxx
  // Apr 20, 2017: GENMEAN_RV -> GENPEAK_RV

  double RV ;
  // char fnam[] = "gen_RV" ;

  // ------------ BEGIN -------------

  RV = exec_GENGAUSS_ASYM(&INPUTS.GENGAUSS_RV)  // pass structure
    + get_zvariation(GENLC.REDSHIFT_CMB,"GENPEAK_RV") ;

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

  char fnam[] = "SIMLIB_INIT_DRIVER" ;

  // --------------- BEGIN --------------

  
  SIMLIB_initGlobalHeader();        // generic init of SIMLIB_GLOBAL_HEADER

  SIMLIB_readGlobalHeader_TEXT();   // open and read global header

  SIMLIB_prepGlobalHeader();        

  //  SIMLIB_openLegacy(); // xxx mark delete 

  SIMLIB_findStart();    // find first LIBID to start reading

  //  debugexit(fnam);

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
	  SIMLIB_PSF_PIXSIG );

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

  char PATH_DEFAULT[MXPATHLEN];
  char *OPENFILE = INPUTS.SIMLIB_OPENFILE;
  char c_get[80];
  int  NTMP, NFILT, gzipFlag ;
  char fnam[] = "SIMLIB_readGlobalHeader_TEXT" ;

  // ---------- BEGIN ----------

  print_banner(fnam);

  sprintf(PATH_DEFAULT, "%s/simlib",  PATH_SNDATA_ROOT );
  fp_SIMLIB = snana_openTextFile(1,PATH_DEFAULT, INPUTS.SIMLIB_FILE, 
				 OPENFILE, &INPUTS.SIMLIB_GZIPFLAG );
  
  if ( fp_SIMLIB == NULL ) {
    sprintf ( c1err, "Cannot open file SIMLIB_FILE" );
    sprintf ( c2err," '%s' ", INPUTS.SIMLIB_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // - - - - - - - - - - - - - - - -
  // read header keywords. Stop when we reach "BEGIN"

  while( (fscanf(fp_SIMLIB, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"SURVEY:") == 0 ) {
      readchar(fp_SIMLIB, SIMLIB_GLOBAL_HEADER.SURVEY_NAME );
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
  //
  int i, NTMP, ifilt, ifilt_obs ;
  char cfilt[4], *FILTERS, *TEL ;
  char fnam[] = "SIMLIB_prepGlobalHeader" ;

  // -------------- BEGIN ------------

  print_banner(fnam);


  char *SURVEY    = SIMLIB_GLOBAL_HEADER.SURVEY_NAME ;
  char *SUBSURVEY = SIMLIB_GLOBAL_HEADER.SUBSURVEY_NAME ;
  
  sprintf(GENLC.SURVEY_NAME,    "%s", SURVEY);
  sprintf(GENLC.SUBSURVEY_NAME, "%s", SUBSURVEY);
  if ( IGNOREFILE(SUBSURVEY) == 0 ) 
    { sprintf(GENLC.SUBSURVEY_NAME,"%s",SURVEY); }
  printf("\t SIMLIB Survey    : %s \n", SURVEY );

  // get integer IDSURVEY from SURVEY string
  read_SURVEYDEF();   GENLC.IDSURVEY = get_IDSURVEY(GENLC.SURVEY_NAME);

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
    // xxx mark delete  ifilt_obs = filtindx_( cfilt, strlen(cfilt) ); 
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
  // allower user-input to over-ride 
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
  if ( strcmp(unit,SIMLIB_PSF_PIXSIG) == 0  ) {  }
  else if ( strcmp(unit,SIMLIB_PSF_ASECFWHM ) == 0 ) {  }
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
  int IDSTART  = INPUTS.SIMLIB_IDSTART ;
  int IDLOCK   = INPUTS.SIMLIB_IDLOCK ;
  int NLIBID   = SIMLIB_GLOBAL_HEADER.NLIBID ;
  int JOBID    = INPUTS.JOBID ; // batch JOBID
  int NJOBTOT  = INPUTS.NJOBTOT ; // totoal number of batch jobs

  int NOPT, IDSEEK, MAXRANSTART, NSKIP_LIBID=-9, NSKIP_EXTRA=0 ;
  int DOSKIP, NREAD, NREPEAT, MXREPEAT, NTMP, NLIBID_EXTRA ;
  time_t t0, t1; 
  double XSKIP, XTMP, flatRan ;
  char LINE[100];
  char fnam[] = "SIMLIB_findStart" ;

  // -------------------- BEGIN -------------

  DOSKIP = 0 ;    IDSEEK = -9 ;

  // check user-input SIMLIB_MAXRANSTART
  MAXRANSTART = INPUTS.SIMLIB_MAXRANSTART ; // default
  if ( MAXRANSTART > 0 ) {
    flatRan     = unix_random() ;
    XSKIP       = (double)MAXRANSTART * flatRan ;
    NSKIP_LIBID = (int)XSKIP + 1 ;
    DOSKIP = 1;
    printf("\t SIMLIB MAXRANSTART after  %d  LIBIDs ", 
	   NSKIP_LIBID ); 
  }

  // for batch job, autom-compute NSKIP 
  if ( NJOBTOT > 0  &&  NLIBID > 100 ) { 

    flatRan     = unix_random() ;
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
void SIMLIB_read_templateNoise(char *FIELD, char *whatNoise) {

  // Feb 2014
  // For *whatNoise = "SKY" or "CCD" or "ZPT", read correlated template 
  // noise for each defined filter and store in 
  // SIMLIB_TEMPLATE_SKY[CCD]SIG array.
  // The template noise values must be listed in the same order
  // as for the FILTERS key.
  // Note that this is different from the T: keys where
  // the template noise is UN-correlated among epochs.
  //
  // Aug 2014: new FIELD arg and major re-org to handle overlapping fields.
  //
  // Sep 29 2015: check that noise is within validNoiseRange to
  //              avoid reading crazy values 
  //
  // Mar 02 2018: check option to ignore template noise

  double noise, *ptrNoise, validNoise_min, validNoise_max;
  int NFILT, ifilt, ifilt_obs, IFIELD, IGNORE  ;
  char fnam[] = "SIMLIB_read_templateNoise" ;
  
  // ------------- BEGIN -----------

  IGNORE = (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_IGNORE_TEMPLATENOISE);
  if ( IGNORE ) { return ; }

  // determine local FIELD index.
  IFIELD = IFIELD_OVP_SIMLIB(0,FIELD);

  // if we did not find IFIELD, then increment for new field
  if ( IFIELD < 0 ) {
    SIMLIB_TEMPLATE.NFIELD_OVP++ ;
    IFIELD = SIMLIB_TEMPLATE.NFIELD_OVP ;
  }
  
  
  if ( IFIELD < 0 || IFIELD >= MXFIELD_OVP_SIMLIB ) {
    sprintf(c1err,"Invalid IFIELD=%d for FIELD=%s, LIBID=%d", 
	    IFIELD, FIELD, SIMLIB_HEADER.LIBID );
    sprintf(c2err,"IFIELD must be 0 to %d", MXFIELD_OVP_SIMLIB-1);
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
    readdouble(fp_SIMLIB, 1, &noise);
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

  return ;
  
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
  char fnam[] = "IFIELD_OVP_SIMLIB" ;
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

  int  REPEAT ;
  char fnam[] = "SIMLIB_READ_DRIVER" ;

  // ------------------ BEGIN ------------------


  // Begin refactored code.
  // read next cadence and load SIMLIB_OBS_RAW array

  // check for option to repeat Cadence to save time 
  GENLC.NGEN_SIMLIB_ID++ ;
  REPEAT = USE_SAME_SIMLIB_ID(2);

  if ( REPEAT == 0 ) {  // process next cadence

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

  int ID, NOBS_EXPECT, NOBS_FOUND, NOBS_FOUND_ALL, ISTORE, ifilt_obs, ifilt;
  int APPEND_PHOTFLAG;
  int NTRY, USEFLAG_LIBID, USEFLAG_MJD, OPTLINE, NWD, NTMP ;
  int   NOBS_SKIP, SKIP_FIELD, SKIP_APPEND, OPTLINE_REJECT  ;
  double PIXSIZE, TEXPOSE_S, MJD ;
  char c_get[80], ctmp[80], *BAND, *UNIT, cline[200] ;
  char *FIELD = SIMLIB_HEADER.FIELD;
  char *TEL   = SIMLIB_HEADER.TELESCOPE ;
  char fnam[] = "SIMLIB_readNextCadence_TEXT" ;

  // ------------ BEGIN --------------

  SIMLIB_HEADER.NWRAP = NTRY = 0 ; // reset wrap before START 
  SIMLIB_OBS_RAW.NOBS = SIMLIB_OBS_RAW.NOBS_READ = 0 ;
  SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH = 0 ;
  NWD=0;


 START:

  init_SIMLIB_HEADER();
  NOBS_EXPECT = NOBS_FOUND = NOBS_FOUND_ALL = USEFLAG_LIBID =USEFLAG_MJD = 0 ;
  NOBS_SKIP = SKIP_FIELD = SKIP_APPEND = APPEND_PHOTFLAG = 0 ;
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

  while( (fscanf(fp_SIMLIB, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"END_OF_SIMLIB:") == 0 ) {

      // check SIMLIB after 5 passes to avoid infinite loop
      if ( SIMLIB_HEADER.NWRAP >= 5 )  { ENDSIMLIB_check(); }
      
      if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENRANDOM ) {
	snana_rewind(fp_SIMLIB, INPUTS.SIMLIB_OPENFILE,
		     INPUTS.SIMLIB_GZIPFLAG);
	SIMLIB_HEADER.NWRAP++ ; 
	SIMLIB_HEADER.LIBID = SIMLIB_ID_REWIND ; 
	NOBS_FOUND = NOBS_FOUND_ALL = USEFLAG_LIBID = USEFLAG_MJD = 0 ;
      }
    }  // end simlib if-block
 
    if ( strcmp(c_get,"LIBID:") == 0 ) {
      readint ( fp_SIMLIB, 1, &ID );
      SIMLIB_HEADER.LIBID = ID ; 
      sprintf(SIMLIB_HEADER.LIBNAME, "LIB%5.5d", ID );
      USEFLAG_LIBID = ACCEPT_FLAG ;
    }

    if ( strcmp(c_get,"RA:") == 0 ) 
      { readdouble(fp_SIMLIB, 1, &SIMLIB_HEADER.RA ); }
    if ( strcmp(c_get,"DEC:") == 0 ) 
      { readdouble(fp_SIMLIB, 1, &SIMLIB_HEADER.DEC ); }
    if ( strcmp(c_get,"DECL:") == 0 ) 
      { readdouble(fp_SIMLIB, 1, &SIMLIB_HEADER.DEC ); }

    if ( strcmp(c_get,"SUBSURVEY:")==0 ) 
      { readchar( fp_SIMLIB, SIMLIB_HEADER.SUBSURVEY_NAME); }

    if ( strcmp(c_get,"FIELD:") == 0 )  { 
      readchar ( fp_SIMLIB, FIELD ); 
      SKIP_FIELD = ( SKIP_SIMLIB_FIELD(FIELD) &&
		     (INPUTS.SIMLIB_FIELDSKIP_FLAG==0) ) ;
    }

    if ( strcmp(c_get,"PIXSIZE:") == 0 )  
      {  readdouble ( fp_SIMLIB, 1, &SIMLIB_HEADER.PIXSIZE ); }

    if ( strcmp(c_get,"TELESCOPE:") == 0 ) 
      { readchar ( fp_SIMLIB, SIMLIB_HEADER.TELESCOPE ); }
  
    if ( strcmp(c_get,"MWEBV:") == 0  )
      { readdouble ( fp_SIMLIB, 1, &SIMLIB_HEADER.MWEBV ); }

    if ( strcmp(c_get,"CCD:") == 0 || strcmp(c_get,"CCDNUM:") == 0 )  
      {  readint ( fp_SIMLIB, 1, &SIMLIB_HEADER.CCDNUM);  }

    // read optional header keys for FAKEID option
    if ( strcmp(c_get,"GALID:") == 0 )  
      { readint ( fp_SIMLIB, 1, &SIMLIB_HEADER.GALID ); }
    if ( strcmp(c_get,"FAKEID:") == 0 )  
      { readint ( fp_SIMLIB, 1, &SIMLIB_HEADER.FAKEID ); }

    // check for APPEND option --> last MJDs are not sorted.
    if ( strcmp(c_get,"APPEND:") == 0 ) { 
      readint ( fp_SIMLIB, 1, &APPEND_PHOTFLAG ); 
      SKIP_APPEND = 0 ; // to do: replace with user input flag
    }
    

    // note that NOBS can exceed MXEPSIM ... 
    // as long as it's less than MXOBS_SIMLIB
    if ( strcmp(c_get,"NOBS:") == 0 )   {  
      readint ( fp_SIMLIB, 1, &NOBS_EXPECT );
      SIMLIB_HEADER.NOBS = NOBS_EXPECT ;

      if ( NOBS_EXPECT >= MXOBS_SIMLIB ) {
	sprintf(c1err,"NOBS=%d exceeds bound for LIBID=%d.", NOBS_EXPECT, ID);
	sprintf(c2err,"Check bound: MXOBS_SIMLIB = %d", MXOBS_SIMLIB );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
      }
    }
    
    // check for correlated template noise ; fill SIMLIB_TEMPLATE
    if ( strcmp(c_get,"TEMPLATE_SKYSIG:") == 0 )
      { SIMLIB_read_templateNoise(FIELD,"SKYSIG"); }
    else if ( strcmp(c_get,"TEMPLATE_CCDSIG:") == 0 )
      { SIMLIB_read_templateNoise(FIELD,"CCDSIG"); }
    else if ( strcmp(c_get,"TEMPLATE_ZPT:") == 0 )
      { SIMLIB_read_templateNoise(FIELD,"ZPT"); }

    // read spectrograph exposure time for template
    int OPTMASK, noTEMPLATE, ISKEY ;
    double TEXPOSE, TSCALE ;
    OPTMASK    = INPUTS.SPECTROGRAPH_OPTIONS.OPTMASK ;
    noTEMPLATE = ( OPTMASK & SPECTROGRAPH_OPTMASK_noTEMPLATE ) ;
    ISKEY      = ( strcmp(c_get,"TEMPLATE_TEXPOSE_SPECTROGRAPH:")==0);
    if ( ISKEY && (noTEMPLATE==0) ) {
      TSCALE=INPUTS.SPECTROGRAPH_OPTIONS.SCALE_TEXPOSE ;
      readdouble(fp_SIMLIB, 1, &TEXPOSE );
      SIMLIB_TEMPLATE.TEXPOSE_SPECTROGRAPH = TEXPOSE * TSCALE ;
    }

    // check for optional GENRANGE or cut keys in header.
    parse_SIMLIB_GENRANGES(fp_SIMLIB,c_get);
    
    // ------------------------------------------------
    // check for epochs
    OPTLINE = 0;
    if ( strcmp(c_get,"S:") == 0 ) {
      NOBS_FOUND_ALL++ ;
      if ( USEFLAG_LIBID == ACCEPT_FLAG ) { OPTLINE = OPTLINE_SIMLIB_S;  }
    }

    if ( strcmp(c_get,"SPECTROGRAPH:") == 0 && USEFLAG_LIBID==ACCEPT_FLAG )
      { OPTLINE = OPTLINE_SIMLIB_SPECTROGRAPH ;  }

    // always check reasons to reject (header cuts, FIELD, APPEND ...)
    OPTLINE_REJECT = ( USEFLAG_LIBID==REJECT_FLAG || 
		       SKIP_FIELD || SKIP_APPEND ) ;

    if ( OPTLINE && OPTLINE_REJECT )  {    
      // MJD line in already rejected LIBID --> read rest of line 
      fgets(cline, 100, fp_SIMLIB) ;
      if ( SKIP_FIELD ) { NOBS_SKIP++ ; }
    }
    else if ( OPTLINE == OPTLINE_SIMLIB_S )  { 
      ISTORE = NOBS_FOUND ;
      SIMLIB_OBS_RAW.OPTLINE[ISTORE] = OPTLINE ;
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.MJD[ISTORE] );

      readchar( fp_SIMLIB, ctmp );
      parse_SIMLIB_IDplusNEXPOSE(ctmp,
				 &SIMLIB_OBS_RAW.IDEXPT[ISTORE],
				 &SIMLIB_OBS_RAW.NEXPOSE[ISTORE] );

      readchar   ( fp_SIMLIB,     SIMLIB_OBS_RAW.BAND[ISTORE]     );
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.CCDGAIN[ISTORE]  );
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.READNOISE[ISTORE]);
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.SKYSIG[ISTORE]   );
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.PSFSIG1[ISTORE]  );
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.PSFSIG2[ISTORE]  );
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.PSFRATIO[ISTORE] );
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.ZPTADU[ISTORE]   );  
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.ZPTSIG[ISTORE]   );  
      readdouble ( fp_SIMLIB, 1, &SIMLIB_OBS_RAW.MAG[ISTORE]      );

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

      NOBS_FOUND++ ;   

    }
    else if ( OPTLINE == OPTLINE_SIMLIB_SPECTROGRAPH  )  { 

      ISTORE = NOBS_FOUND ;
      readdouble( fp_SIMLIB, 1, &MJD );
      readdouble( fp_SIMLIB, 1, &TEXPOSE_S );

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

      NOBS_FOUND++ ;

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
    }

    // stop reading when we reach the end of this LIBID
    if ( strcmp(c_get, "END_LIBID:") == 0 ) { 
      if ( USEFLAG_LIBID == ACCEPT_FLAG ) 
	{ goto DONE_READING ; }
      else
	{ goto START ; }      // read another 
    }

  } // end while loop with fscanf

 DONE_READING:

  /*
  printf(" xxx %s  NOBS(EXPECT,FOUND,SKIP) = %d,%d,%d   FIELD=%s \n",
	 fnam, NOBS_EXPECT, NOBS_FOUND, NOBS_SKIP, FIELD );
  */

  if ( INPUTS.SIMLIB_FIELDSKIP_FLAG==0 && ( NOBS_EXPECT==NOBS_SKIP) ) 
    { goto START ; }


  SIMLIB_OBS_RAW.NOBS      = NOBS_FOUND ;  // can change with SPECTROGRAPH
  SIMLIB_OBS_RAW.NOBS_READ = NOBS_FOUND ;  // won't change for this cadence

  NOBS_EXPECT -= NOBS_SKIP ;
  if ( NOBS_EXPECT != NOBS_FOUND ) {
    sprintf(c1err,"Found %d observations in LIBID %d", 
	    NOBS_FOUND,  SIMLIB_HEADER.LIBID );
    sprintf(c2err,"But expected NOBS = %d from simlib header (NOBS_SKIP=%d)", 
	    NOBS_EXPECT, NOBS_SKIP );
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

  char fnam[] = "SIMLIB_randomize_skyCoords" ;

  // --------------- BEGIN -------------

  return ;

} // end SIMLIB_randomize_skyCoords


// ==============================================
int keep_SIMLIB_HEADER(void) {

  // Created Aug 2017
  // Return ACCEPT_FLAG if SIMLIB_HEADER values pass cuts.
  // Return REJECT_FLAG if SIMLIB_HEADER fails cuts, or REWIND flag is set
  //

  int  ID      = SIMLIB_HEADER.LIBID ;
  int  NOBS    = SIMLIB_HEADER.NOBS ;
  int  IDLOCK  = GENLC.SIMLIB_IDLOCK ;
  int  NSKIP   = INPUTS.NSKIP_SIMLIB ;
  int  NWRAP   = SIMLIB_HEADER.NWRAP ;
  double  zCMB = GENLC.REDSHIFT_CMB ;
  int  iskip, icheck ;
  double *ptrGen ;
  char fnam[] = "keep_SIMLIB_HEADER" ;
  int LTRACE = 0 ;

  // ----------- BEGIN ---------------

  if(LTRACE) {
    printf(" xxx %s: ----------- ID=%d ------------ \n", 
	   fnam, GENLC.SIMLIB_ID );
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
  if(LTRACE) {printf(" xxx %s: 8\n", fnam );}
  ptrGen = SIMLIB_HEADER.GENRANGE_REDSHIFT ;
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


  // Nov 26 2017: check NREPEAT for LCLIB
  if(LTRACE) {printf(" xxx %s: 12\n", fnam );}
  //  printf(" 1. xxx %s: HEADER.NREPEAT=%d   INPUTS.NREPEAT=%d\n", 
  //	 fnam, SIMLIB_HEADER.NREPEAT, INPUTS.SIMLIB_NREPEAT );
  set_SIMLIB_NREPEAT();
  if ( INPUTS.SIMLIB_NREPEAT == 0 ) { return(REJECT_FLAG); }
  

  if(LTRACE) {printf(" xxx %s: 99 END\n", fnam ); debugexit(fnam); }
  
  // if we get here, keep this ID
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
  char fnam[] = "SIMLIB_addCadence_SPECTROGRAPH" ;

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
      SIMLIB_OBS_RAW.ZPTADU[ISTORE]     = -9.0 ;
      SIMLIB_OBS_RAW.ZPTSIG[ISTORE]     = -9.0 ;
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

  int i, OPT, ifilt, ifilt_obs, OBSRAW, NOBS_ADD, OPT_FRAME, IS_HOST ;
  int NSPEC = NPEREVT_TAKE_SPECTRUM ;
  int NFILT = GENLC.NFILTDEF_SPECTROGRAPH ;  // all spectroscopic filters
  double EPOCH[2], MJD, MJD_REF;
  double z, z1, TOBS, TREST, TOBS_REF ;
  double TEXPOSE, VAL_STORE[8] ;
  char *FIELD ;
  char fnam[] = "SIMLIB_TAKE_SPECTRUM" ;


  // -------------- BEGIN -------------

  if ( NSPEC == 0 ) { return ; }

  // in case of repeated cadence, reset NOBS to leave out 
  // previous TAKE_SPECTRUM obs.
  SIMLIB_OBS_RAW.NOBS = SIMLIB_HEADER.NOBS =
    SIMLIB_OBS_RAW.NOBS_READ + SIMLIB_OBS_RAW.NOBS_SPECTROGRAPH ;
  

  OBSRAW = SIMLIB_OBS_RAW.NOBS ;  NOBS_ADD=0;
 
  z = GENLC.REDSHIFT_CMB ;
  for(i=0; i < NSPEC; i++ ) {  

    OPT_FRAME = INPUTS.TAKE_SPECTRUM[i].OPT_FRAME_EPOCH ;
    IS_HOST   = (OPT_FRAME == GENFRAME_HOST);
    EPOCH[0] = INPUTS.TAKE_SPECTRUM[i].EPOCH_RANGE[0] ;
    EPOCH[1] = INPUTS.TAKE_SPECTRUM[i].EPOCH_RANGE[1] ;

    MJD_REF =  GENSPEC_PICKMJD( 0, i, z, &TOBS, &TREST) ; // for MJD sorting

    // make sure that TREST is within valid range
    if ( !IS_HOST ) {
      float *ptrTcut = INPUTS.GENRANGE_TREST ;
      if ( TREST <= ptrTcut[0] ||	 TREST >= ptrTcut[1]  ) {
	sprintf(c1err,"Invalid TREST=%.2f in TAKE_SPECTRUM key.",TREST);
	sprintf(c2err,"User set 'GENRANGE_TREST:  %.2f  %.2f' ", 
		ptrTcut[0], ptrTcut[1] );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
      }
    }

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


    if ( NFILT==0 ) {  // no synthetic filters; just spectra
      VAL_STORE[0] = MJD_REF   ; // same MJD_REF each event for MJD-sorting
      VAL_STORE[1] = TEXPOSE ;
      VAL_STORE[2] = (double)i ;   
      store_GENSPEC(VAL_STORE) ; 
    }

    for(ifilt=0; ifilt < NFILT; ifilt++ ) {  

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
      
      OBSRAW++ ;   NOBS_ADD++ ;	

    } // end ifilt

  } // end i-loop over TAKE_SPECTRUM


  SIMLIB_OBS_RAW.NOBS += NOBS_ADD ;
  SIMLIB_HEADER.NOBS  += NOBS_ADD ;


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

  int NEW_CADENCE = (REPEAT_CADENCE == 0 ) ;
  int ORDER_SORT, ISTORE,  OPTLINE, OBSRAW ;
  double zcmb, RA, DEC, vpec, MJDrange[2], PIXSIZE, DUM, PSF_ORIG ;
  double SNR, SNRMAX=-9.0 ;

  int NOBS_APPEND = SIMLIB_HEADER.NOBS_APPEND ;
  int NOBS_RAW    = SIMLIB_HEADER.NOBS ;

  char *UNIT, *BAND ;
  char fnam[] = "SIMLIB_prepCadence" ;

  // --------------- BEGIN ----------------

  GENLC.NEPOCH = 0 ; // Mar 17 2015: needed for LIBID repeats

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
      checkval_D("ZPTSIG", 1, &SIMLIB_OBS_RAW.ZPTSIG[ISTORE],  0.0, 5.0 ) ;
      checkval_D("PSF1",   1, &SIMLIB_OBS_RAW.PSFSIG1[ISTORE], 0.0, 9.9 );
      checkval_D("PSF2",   1, &SIMLIB_OBS_RAW.PSFSIG2[ISTORE], 0.0, 9.9 );
      checkval_D("PSFrat", 1, &SIMLIB_OBS_RAW.PSFRATIO[ISTORE],0.0, 1.0 );
      checkval_D("SKYSIG", 1, &SIMLIB_OBS_RAW.SKYSIG[ISTORE],  0.0, 1.0E5 );

      PIXSIZE = SIMLIB_OBS_RAW.PIXSIZE[ISTORE] ; 

      // 3a. unit check for optional units of PSF and SKYSIG
      UNIT = SIMLIB_GLOBAL_HEADER.PSF_UNIT ;
      if ( strcmp(UNIT,SIMLIB_PSF_ASECFWHM ) == 0 ) {
	PSF_ORIG = SIMLIB_OBS_RAW.PSFSIG1[ISTORE] ;
	// convert FWHM(arcsec) back to Sigma(pixels)
	SIMLIB_OBS_RAW.PSFSIG1[ISTORE] /= (PIXSIZE * FWHM_SIGMA_RATIO);
	SIMLIB_OBS_RAW.PSFSIG2[ISTORE] /= (PIXSIZE * FWHM_SIGMA_RATIO);
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
  double MJD, CCDGAIN, RDNOISE, SKYSIG, PSF[3], ZPT[2], MAG ;
  double SKYSIG_T, RDNOISE_T, ZPT_T ;
  double SHIFT_ZPT, SCALE_SKYSIG, SCALE_RDNOISE, SCALE_PSF ;
  double MJD_DIF, MJD_LAST_KEEP, DT, DUMMY_STORE[3] ;
  char   *TEL, *FIELD, cfilt[2];

  // init stuff before loop over MJDs
  NEP=NEP_NEWMJD=0;  MJD_LAST_KEEP=-9.0;  

  // transfer OBS_RAW to OBS_GEN; latter has cuts and is sorted by MJD

  for ( isort = 0; isort < NOBS_RAW; isort++ ) {

    OBSRAW   = SIMLIB_LIST_forSORT.INDEX_SORT[isort]; 
    OPTLINE  = SIMLIB_OBS_RAW.OPTLINE[OBSRAW] ;

    if ( OPTLINE < 0 ) { continue; }

    if ( OPTLINE == OPTLINE_SIMLIB_SPECTROGRAPH ) {
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
    ZPT[0]     = SIMLIB_OBS_RAW.ZPTADU[OBSRAW] ;
    ZPT[1]     = SIMLIB_OBS_RAW.ZPTSIG[OBSRAW] ;
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

    // compute a few things from OBS_RAW
    if ( INPUTS.SMEARFLAG_ZEROPT == 0 ) { ZPT[1] = 0.0 ; }
    if ( MAG < 5.0 || MAG > 50 ) { MAG = MAG_UNDEFINED ; }
    sprintf(cfilt,  "%c", FILTERSTRING[IFILT_OBS] ) ;

    ifilt     = GENLC.IFILTINVMAP_SIMLIB[IFILT_OBS] ;
    if ( ifilt < 0 || ifilt >= GENLC.NFILTDEF_SIMLIB ) 
      {  ABORT_SIMLIB_FILTER(isort);  }
    //      {  ABORT_SIMLIB_FILTER(OPTLINE,MJD,cfilt);  }

    // check if this MJD is kept.
    KEEP = keep_SIMLIB_OBS(isort,REPEAT_CADENCE);

    if ( GENLC.SIMLIB_ID < -99 ) {
      printf(" xxx isort=%3d OBSRAW=%3d  IFILT_OBS=%2d(%s) "
	     "MJD=%.3f  KEEP=%d\n",
	     isort, OBSRAW, IFILT_OBS, BAND, MJD, KEEP ); // xxxxxx
      fflush(stdout);
    }

    if ( KEEP == 0 ) { continue ; }

    // mark this filter as 'used'
    GENLC.SIMLIB_USEFILT_ENTRY[IFILT_OBS] = 1;

    // get optional fudge scales
    get_SIMLIB_SCALES(IFILT_OBS, &SHIFT_ZPT, &SCALE_SKYSIG, &SCALE_RDNOISE);
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
    GENLC.genmag8_obs[NEP]   = MAG ;
    GENLC.MJD[NEP]           = MJD ;
    sprintf( GENLC.FIELDNAME[NEP], "%s", FIELD );
    sprintf( GENLC.TELESCOPE[NEP], "%s", TEL   );

    // store SIMLIB_GEN quanties vs. NEP that include fudges                  
    SIMLIB_OBS_GEN.MJD[NEP]         = MJD ;
    SIMLIB_OBS_GEN.IFILT_OBS[NEP]   = IFILT_OBS ;
    SIMLIB_OBS_GEN.PIXSIZE[NEP]     = PIXSIZE ;
    SIMLIB_OBS_GEN.CCDGAIN[NEP]     = CCDGAIN ;
    SIMLIB_OBS_GEN.READNOISE[NEP]   = RDNOISE * SCALE_RDNOISE ;
    SIMLIB_OBS_GEN.SKYSIG[NEP]      = SKYSIG * SCALE_SKYSIG ;
    SIMLIB_OBS_GEN.PSFSIG1[NEP]     = PSF[0] * SCALE_PSF ;
    SIMLIB_OBS_GEN.PSFSIG2[NEP]     = PSF[1] * SCALE_PSF ;
    SIMLIB_OBS_GEN.PSFRATIO[NEP]    = PSF[2] ;  // ratio                      
    SIMLIB_OBS_GEN.ZPTADU[NEP]      = ZPT[0] + SHIFT_ZPT ;
    SIMLIB_OBS_GEN.ZPTSIG[NEP]      = ZPT[1] ;
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
      if ( IFIELD < 0 || IFIELD >= MXFIELD_OVP_SIMLIB ) {
	sprintf(c1err,"Invalid IFIELD=%d for template FIELD=%s",IFIELD,FIELD);
	sprintf(c2err, "Check LIBID=%d  %s-band", GENLC.SIMLIB_ID,cfilt );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;
      }  

      SKYSIG_T = SIMLIB_TEMPLATE.SKYSIG[IFIELD][IFILT_OBS] * SCALE_SKYSIG ;
      SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[NEP]    = SKYSIG_T ;	
      SIMLIB_OBS_RAW.TEMPLATE_SKYSIG[OBSRAW] = SKYSIG_T ;
      
      RDNOISE_T = SIMLIB_TEMPLATE.READNOISE[IFIELD][IFILT_OBS]*SCALE_RDNOISE ;
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

  if ( GENLC.NFILTDEF_SPECTROGRAPH == 0 ) 
    { ifilt_obs=0 ; goto  STORE_SIMLIB_RAW ; }

  ifilt_obs = GENLC.IFILTDEF_SPECTROGRAPH[ifilt] ;
  sprintf(cfilt,  "%c", FILTERSTRING[ifilt_obs] ) ;


  /* xxxxxxxxx mark delete xxxxxxxx
  printf(" xxx %s: ISTORE=%d  MJD=%.3f   ifilt=%d  ifilt_obs=%d \n",
	 fnam, ISTORE, MJD, TEXPOSE_S, ifilt, ifilt_obs ); 
  fflush(stdout);
  xxxxxxxxxxxxxxxxx */


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
  SIMLIB_OBS_RAW.ZPTADU[ISTORE]     = ZPT ;  // TEXPOSE zero point (pe)
  SIMLIB_OBS_RAW.ZPTSIG[ISTORE]     = 0.0 ;  // no error; included in SKYSIG
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
    printf("\n PRE-ABORT DUMP \n");
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
int keep_SIMLIB_OBS(int isort, int REPEAT) {

  // Created Aug 2017
  // Return 1 of this raw-OBS is kept.
  // Return 0 to reject this OBS.
  //
  // Inputs:
  //   OBS      = MJD-sort index; needed to get absolute OBS index
  //   REPEAT   = 1 if this cadence is repeated
  //
  // Sep 24 2017: check field.

  int  KEEP=1, NOKEEP=0; 
  int  ifilt, ifilt_obs, OBS ;
  int  LTRACE=0; 
  double  MJD, MJDrange[2];
  char *FIELD ;
  char fnam[] = "keep_SIMLIB_OBS" ;

  // ------------ BEGIN --------------

  
  if (LTRACE) {
    printf(" xxx ---------------------------- \n");
    printf(" xxx 0 isort=%d REPEAT=%d \n",isort,REPEAT); 
  }
  OBS  = SIMLIB_LIST_forSORT.INDEX_SORT[isort] ; // absolute OBS_RAW index


  FIELD = SIMLIB_OBS_RAW.FIELDNAME[OBS] ;
  MJD   = SIMLIB_OBS_RAW.MJD[OBS] ;

  // compute & store SEASON info; return MJDrange to keep SIMLIB entries
  set_SIMLIB_MJDrange(REPEAT,MJDrange);

  if ( MJD < MJDrange[0] ) { return(NOKEEP); }
  if ( MJD > MJDrange[1] ) { return(NOKEEP); }

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
  // xxx mark delete   ifilt     = GENLC.IFILTINVMAP_SIMLIB[ifilt_obs] ; 
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
void set_SIMLIB_MJDrange(int sameFlag, double *MJDrange) {

  // 
  // Return MJDrange
  //
  // Input: sameFlag=T if LIBID has not changed.
  //    (sameFlag not used; obsolete)
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
  //

  int    KEEP_ENTIRE_SEASON = 
    (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_ENTIRE_SEASON );
  int    KEEP_ENTIRE_SURVEY = 
    (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_ENTIRE_SURVEY );

  double PEAKMJD = GENLC.PEAKMJD ;
  double z       = GENLC.REDSHIFT_CMB ;
  double z1      = 1. + z ;
  double Tpad = 0.1 ;
  double Tmin = PEAKMJD + (z1* INPUTS.GENRANGE_TREST[0] ) - Tpad ;
  double Tmax = PEAKMJD + (z1* INPUTS.GENRANGE_TREST[1] ) + Tpad ;
  double TMPmin, TMPmax;
  int    ISEASON; 
  char fnam[] = "set_SIMLIB_MJDrange" ;

  // -------------- BEGIN ------------

  MJDrange[0] = -9.0 ;
  MJDrange[1] = -9.0 ;

  if ( KEEP_ENTIRE_SEASON && KEEP_ENTIRE_SURVEY ) {
    sprintf(c1err,"Cannot specify ENTIRE_SEASON (%d) and ENTIRE_SURVEY(%d)",
	    SIMLIB_MSKOPT_ENTIRE_SEASON, SIMLIB_MSKOPT_ENTIRE_SURVEY) ;
    sprintf(c2err,"Pick one or neither in SIMLIB_MSKOPT key.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

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
  int LDMP = 0 ;
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
  char fnam[] = "SIMLIB_prepMJD_forSORT" ;

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

  int NOBS_RAW    = SIMLIB_HEADER.NOBS ;
  int NOBS_APPEND = SIMLIB_HEADER.NOBS_APPEND ;
  int NOBS_SORT   = NOBS_RAW - NOBS_APPEND ;
  int  ORDER_SORT = +1 ; // increasing order
  int  isort;
  char fnam[] = "SIMLIB_sortbyMJD" ;

  // ------------- BEGIN --------------

  if ( NOBS_SORT < 1 ) {
    sprintf(c1err,"Invalid NOBS_SORT=%d for LIBID=%d CID=%d", 
	    NOBS_SORT, SIMLIB_HEADER.LIBID, GENLC.CID );
    sprintf(c2err,"NOBS_RAW=%d  NOBS_APPEND=%d", 
	    NOBS_RAW, NOBS_APPEND);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  sortDouble( NOBS_SORT, SIMLIB_LIST_forSORT.MJD, ORDER_SORT, 
	      SIMLIB_LIST_forSORT.INDEX_SORT ) ;  // return array
  
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
    { GENSPEC.IS_HOST[imjd] = 1; }  // HOST spectrum

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
  char fnam[] = "init_SIMLIB_HEADER" ;
  // ------------ BEGIN ------------

  SIMLIB_HEADER.REGEN_FLAG  = 0 ;
  SIMLIB_HEADER.NOBS        = NULLINT ;
  SIMLIB_HEADER.NOBS_APPEND = 0 ;
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

  sprintf(SIMLIB_HEADER.FIELD,"NULL");
  SIMLIB_HEADER.SUBSURVEY_NAME[0] = 0 ;
  // doesn't work  sprintf(SIMLIB_HEADER.SUBSURVEY_NAME, "NULL" );

  SIMLIB_HEADER.NSEASON = 1 ;
  for(i=0; i < MXSEASON_SIMLIB; i++ ) { SIMLIB_HEADER.USE_SEASON[i] = 0;}


  sprintf(SIMLIB_HEADER.TELESCOPE,  "%s", SIMLIB_GLOBAL_HEADER.TELESCOPE) ;
  SIMLIB_HEADER.PIXSIZE  = SIMLIB_GLOBAL_HEADER.PIXSIZE ;


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
  char fnam[] = "parse_SIMLIB_IDplusNEXPOSE" ;

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
void parse_SIMLIB_GENRANGES(FILE *fp_SIMLIB, char *KEY) {

  // Created Aug 2017
  // Read optional GENRANGE_XXX keys.
  //
  // Inputs:
  //   fp_SIMLIB : file pointer to SIMLIB file
  //   KEY       : current simlib string to check
  //
  // Jan 4 2018: read optional DISTANCE key, and convert to zCMB

  double RA  = SIMLIB_HEADER.RA ;
  double DEC = SIMLIB_HEADER.DEC ;
  double TMPVAL, TMPRANGE[2], dist, MU, zHel;
  int    LTMP ;

  char fnam[] = "parse_SIMLIB_GENRANGES" ;

  // ------------ BEGIN ----------------


  // Nov 2015: check optional CUTWIN_REDSHIFT so that each LIBID
  //           can have its own z-range (originally for WFIRST study)
  //    Note that this is a cut, not a range to generate.
  if (  strcmp(KEY,"CUTWIN_REDSHIFT:") == 0  ||
	strcmp(KEY,"REDSHIFT_RANGE:" ) == 0  ) 
    { readdouble( fp_SIMLIB, 2, SIMLIB_HEADER.CUTWIN_REDSHIFT );   }
  
  // -------------------------------------------------
  // ------------ READ GENRANGEs ---------------------
  // -------------------------------------------------

  // check for redshift value or range 
  LTMP = 0 ;
  if ( strcmp(KEY,"REDSHIFT:")==0 && INPUTS.USE_SIMLIB_REDSHIFT ) {
    readdouble ( fp_SIMLIB, 1, &TMPVAL );
    TMPRANGE[0] = TMPRANGE[1] = TMPVAL;  LTMP=1;
  }
  if ( strcmp(KEY,"DISTANCE:")==0 && INPUTS.USE_SIMLIB_DISTANCE ) {
    readdouble ( fp_SIMLIB, 1, &dist );  // Lumi-distance (D_L), Mpc (Jan 2018)
    MU = 5.0*log10(dist/1.0E-5);  // 10pc = 1.0E-5 Mpc
    double H0 = INPUTS.H0 / ( 1.0E6 * PC_km) ;
    double OM = INPUTS.OMEGA_MATTER;
    double OL = INPUTS.OMEGA_LAMBDA;
    double w0 = INPUTS.W0_LAMBDA ;
    TMPVAL = zcmb_dLmag_invert(H0, OM, OL, w0, MU); // returns zCMB
    TMPRANGE[0] = TMPRANGE[1] = TMPVAL;  LTMP=1;
  }
  if ( strcmp(KEY,"GENRANGE_REDSHIFT:")==0 ) {
    readdouble ( fp_SIMLIB, 2, TMPRANGE );  LTMP=1;
  }
  if ( LTMP ) {
    SIMLIB_HEADER.GENRANGE_REDSHIFT[0] = TMPRANGE[0] ;
    SIMLIB_HEADER.GENRANGE_REDSHIFT[1] = TMPRANGE[1] ;
    SIMLIB_HEADER.REGEN_FLAG = 1;
  }

  // check for optional PEAKMJD value or range in SIMLIB; 
  LTMP=0;
  if ( strcmp(KEY,"PEAKMJD:")==0  && INPUTS.USE_SIMLIB_PEAKMJD )  { 
    readdouble ( fp_SIMLIB, 1, &TMPVAL ); 
    TMPRANGE[0] = TMPRANGE[1] = TMPVAL;  LTMP=1;
  }
  if ( strcmp(KEY,"GENRANGE_PEAKMJD:") == 0   || 
       strcmp(KEY,"GENRANGE_PKMJD:")   == 0  ) {
    readdouble ( fp_SIMLIB, 2, TMPRANGE ) ;  LTMP=1 ;
  }
  if ( LTMP ) {
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.PEAK     = 0.5*(TMPRANGE[0]+TMPRANGE[1]);
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[0] = TMPRANGE[0] ;
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.RANGE[1] = TMPRANGE[1] ;
    SIMLIB_HEADER.REGEN_FLAG = 1;
  }

  // check for PEAKMJD sigma
  if ( strcmp(KEY,"GENSIGMA_PEAKMJD:")==0  ) {
    readdouble ( fp_SIMLIB, 1, &TMPVAL ); 
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.SIGMA[0] = TMPVAL ;
    SIMLIB_HEADER.GENGAUSS_PEAKMJD.SIGMA[1] = TMPVAL ;
  }


  // check for SALT2c range & sigma
  if ( strcmp(KEY,"GENRANGE_SALT2c:") == 0 ) {
    readdouble ( fp_SIMLIB, 2, TMPRANGE );
    SIMLIB_HEADER.GENGAUSS_SALT2c.PEAK     = 0.5*(TMPRANGE[0]+TMPRANGE[1]);
    SIMLIB_HEADER.GENGAUSS_SALT2c.RANGE[0] = TMPRANGE[0] ;
    SIMLIB_HEADER.GENGAUSS_SALT2c.RANGE[1] = TMPRANGE[1] ;
    SIMLIB_HEADER.REGEN_FLAG = 1;
  }
  if ( strcmp(KEY,"GENSIGMA_SALT2c:") == 0 ) {
    readdouble ( fp_SIMLIB, 1, &TMPVAL ); 
    SIMLIB_HEADER.GENGAUSS_SALT2c.SIGMA[0] = TMPVAL ;
    SIMLIB_HEADER.GENGAUSS_SALT2c.SIGMA[1] = TMPVAL ;
  } 

  // check for SALT2x1 range & sigma
  if ( strcmp(KEY,"GENRANGE_SALT2x1:")==0 ) {
    readdouble ( fp_SIMLIB, 2, TMPRANGE );   
    SIMLIB_HEADER.GENGAUSS_SALT2x1.PEAK     = 0.5*(TMPRANGE[0]+TMPRANGE[1]);
    SIMLIB_HEADER.GENGAUSS_SALT2x1.RANGE[0] = TMPRANGE[0] ;
    SIMLIB_HEADER.GENGAUSS_SALT2x1.RANGE[1] = TMPRANGE[1] ;
    SIMLIB_HEADER.REGEN_FLAG = 1;
  }
  if ( strcmp(KEY,"GENSIGMA_SALT2x1:") == 0 ) {
    readdouble ( fp_SIMLIB, 1, &TMPVAL ); 
    SIMLIB_HEADER.GENGAUSS_SALT2x1.SIGMA[0] = TMPVAL ;
    SIMLIB_HEADER.GENGAUSS_SALT2x1.SIGMA[1] = TMPVAL ;
  } 

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


int regen_SIMLIB_GENRANGES(void) {

  // Created Apr 13 2016 by R.Kessler
  // called after reading each LIBID header,
  // to check if anything needs to be re-generated.
  //
  // Functions returns +1 to keep event, -1 to reject.
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
    printf(" xxx -----------CID=%d LIBID=%d -------------- \n",
	   GENLC.CID, GENLC.SIMLIB_ID ) ;
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
  double *ZWIN = SIMLIB_HEADER.GENRANGE_REDSHIFT ;
  if(LTRACE) {printf(" xxx %s: 1 ZWIN=%f,%f\n", fnam,ZWIN[0],ZWIN[1] ); }
  if ( ZWIN[0] > 1.0E-5 ) {
    tmpVal = genz_hubble( ZWIN[0], ZWIN[1], &INPUTS.RATEPAR );
    if(LTRACE) {printf(" xxx %s: 1b ZGEN=%f\n", fnam,tmpVal ); }
    if ( tmpVal < INPUTS.GENRANGE_REDSHIFT[0] ) { return(REJECT); }
    if ( tmpVal > INPUTS.GENRANGE_REDSHIFT[1] ) { return(REJECT); }
    GENLC.REDSHIFT_CMB     = tmpVal ;
    gen_distanceMag(GENLC.REDSHIFT_CMB, GENLC.REDSHIFT_CMB, // to fix later
		    &GENLC.DLMU, &GENLC.LENSDMU);
  }

  // check PEAKMJD update
  if(LTRACE) {printf(" xxx %s: 2 PEAKMJD=%f\n", fnam,GENLC.PEAKMJD ); }
  if ( SIMLIB_HEADER.GENGAUSS_PEAKMJD.PEAK > 1000.0 ) {
    tmpVal         = exec_GENGAUSS_ASYM(&SIMLIB_HEADER.GENGAUSS_PEAKMJD) ;
    GENLC.PEAKMJD  = tmpVal ;  
    if ( tmpVal < INPUTS.GENRANGE_PEAKMJD[0] ) { return(REJECT); }
    if ( tmpVal > INPUTS.GENRANGE_PEAKMJD[1] ) { return(REJECT); }
  }


  if(LTRACE) {printf(" xxx %s: 3 SALT2x1=%f\n", fnam,GENLC.SALT2x1 ); }
  if ( fabs(SIMLIB_HEADER.GENGAUSS_SALT2x1.PEAK) < 90. ) {
    tmpVal         = exec_GENGAUSS_ASYM(&SIMLIB_HEADER.GENGAUSS_SALT2x1);
    GENLC.SALT2x1  = tmpVal ;
    if ( tmpVal < INPUTS.GENGAUSS_SALT2x1.RANGE[0] ) { return(REJECT); }
    if ( tmpVal > INPUTS.GENGAUSS_SALT2x1.RANGE[1] ) { return(REJECT); }
  }

  if(LTRACE) {printf(" xxx %s: 4 SALT2c=%f\n", fnam,GENLC.SALT2c ); }
  if ( fabs(SIMLIB_HEADER.GENGAUSS_SALT2c.PEAK) < 90. ) {
    tmpVal        = exec_GENGAUSS_ASYM(&SIMLIB_HEADER.GENGAUSS_SALT2c) ;
    GENLC.SALT2c  = tmpVal ;
    if ( tmpVal < INPUTS.GENGAUSS_SALT2c.RANGE[0] ) { return(REJECT); }
    if ( tmpVal > INPUTS.GENGAUSS_SALT2c.RANGE[1] ) { return(REJECT); }
  }


  // re-compute x0 and mB for SALT2 model
  if ( INDEX_GENMODEL == MODEL_SALT2 ) {
    GENLC.SALT2x0 = SALT2x0calc(GENLC.SALT2alpha, GENLC.SALT2beta, 
				GENLC.SALT2x1, GENLC.SALT2c, GENLC.DLMU ) ;
    GENLC.SALT2mB = SALT2mBcalc(GENLC.SALT2x0) ;
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

  int NOBS_RAW = SIMLIB_HEADER.NOBS ;
  int OBSRAW   = SIMLIB_LIST_forSORT.INDEX_SORT[isort]; 
  int OPTLINE  = SIMLIB_OBS_RAW.OPTLINE[OBSRAW] ;
  int IFILT_OBS= SIMLIB_OBS_RAW.IFILT_OBS[OBSRAW] ;
  double MJD   = SIMLIB_OBS_RAW.MJD[OBSRAW] ;
  char cfilt[2];
  char fnam[]  = "ABORT_SIMLIB_FILTER" ;

  // ----------- BEGIN ------------


  sprintf(cfilt,  "%c", FILTERSTRING[IFILT_OBS] ) ;

  printf("\n PRE-ABORT DUMP: \n" );
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
  char fnam[] = "USE_SAME_SIMLIB_ID" ;

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
    printf(" xxx --> SAMEFLAG[OLD,NEW] = %d, %d \n",
	   OLD_SAMEFLAG, NEW_SAMEFLAG );    
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

  if ( XDIF < FlatRan1(1) ) 
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
  // Returns 1 if this field should be skipped.
  // Returns 0 if this field should be kept.
  //

  if ( strcmp(INPUTS.SIMLIB_FIELDLIST,"ALL") == 0 ) 
    { return 0 ; }

  if ( strstr(INPUTS.SIMLIB_FIELDLIST,field) == NULL ) 
    { return 1; }
  else 
    { return 0; }

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
  //
  // --------- BEGIN ----------
  
  *SCALE_SKYSIG    = (double)INPUTS.FUDGESCALE_SKYNOISE ;
  *SCALE_READNOISE = (double)INPUTS.FUDGESCALE_READNOISE ;
  *SHIFT_ZPT       = (double)INPUTS.FUDGESHIFT_ZPT_FILTER[ifilt_obs] ;


  if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 0) )
    { *SHIFT_ZPT  += GENLC.SHIFT_ZPTSIMLIB[ifilt_obs]; }


  if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 1) )
    { *SCALE_SKYSIG  *= GENLC.SCALE_NOISE[ifilt_obs] ; }


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

  char fnam[] = "ENDSIMLIB_check";

  // ------- BEGIN ---------

  if ( INPUTS.SIMLIB_DUMP >= 0 ) return ;

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

  char *ptrZfile, *ptrparname, *ptrPar;
  char GENPREFIX[60], c_get[60], method[20], parName[60], cpoly[60] ;

#define NPREFIX_GENGAUSS 15
  char PREFIX_GENGAUSS[NPREFIX_GENGAUSS][20] 
    = 
    { "PEAK", "MEAN", "SKEW[0]", "SKEW[1]", "SIGMA[0]", "SIGMA[1]",
      "SKEWNORMAL[0]", "SKEWNORMAL[1]", "SKEWNORMAL[2]",
      "PROB2", "PEAK2", "SIGMA2[0]", "SIGMA2[1]" ,
      "EXPSIG", "EXPTAU"
    } ;

  double *ptrpoly, *ptrzval, *ptrzshift, ZTMP, ZMIN, ZMAX, ZGEN[2], shift ;
  int i, i2, NZ, MATCH, NTMP, FLAG, ipar, IPAR_USR;
  int IZVAR_DEJA, IZVAR_FILE ;
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
    update_PARDEF_ZVAR( "GENGAUSIG_AV"       );
    update_PARDEF_ZVAR( "GENGAUPEAK_AV"      );
    update_PARDEF_ZVAR( "GENRATIO_AV0"       );
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
  if ( (fpz = fopen(ptrZfile, "rt"))==NULL ) {   
    sprintf ( c1err, "Cannot open ZVARIATION_FILE " );
    sprintf ( c2err," '%s' ", ptrZfile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
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
      parse_GENPOLY(cpoly,&INPUT_ZVARIATION[IPAR_USR].POLY, fnam);
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

    printf("\t Init Z-variation (zvar) for '%s' (method = %s)\n", 
	   ptrparname, method );

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
    printf("\n PRE-ABORT DUMP: \n") ;
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

  double shift, tmp, ZMIN, ZMAX, ZTMP, zpow, FMIN, FMAX  ;
  int  ipar, ipar_match, MATCH, FLAG, i, NZ, iz, izmin, izmax ;
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
      sprintf(c1err,"Undefined parname = '%s' at Z=%f", parname, z );
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

    /************* mark delete Apr 10 2019 xxxxxxxxxx
    zpow = 1.0 ;
    for ( i=0; i <= POLYORDER_ZVAR ; i++ ) {  // polynomial function
      tmp = INPUT_ZVARIATION[ipar].ZPOLY[i] ;
      shift += (tmp * zpow) ;
      zpow *= z ; // avoid using slow pow function for powers of z    
    }
    xxxxxxxx end mark xxxxxxxxxxxxx*/

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


  // check for variation in skewNormal (Sep 2016)
  for(i=0; i < 3; i++ ) {
    sprintf(PARNAME, "GENSKEWNORMAL[%d]_%s", i, parName );
    GENGAUSS_OUT.SKEWNORMAL[i] += get_zvariation(z,PARNAME) ;
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

  int NPICKRAN_ABORT ; // abort after this many tries
  int i, i2, j, NPICKRAN, NSTORE_ALL, NSTORE, CIDRAN, CIDTMP, CIDADD;
  int CIDMAX, CIDMIN, USED, NTRY, NCID_MALLOC, NCALL_random=0, LDMP=0 ;
  int L_STDOUT;

  double r8, XN8, XNUPD8;

  char fnam[] = "init_CIDRAN";

  // ------------ BEGIN ------------

  fflush(stdout);
  if ( WRFLAG_CIDRAN  <= 0 ) { return ; }

  CIDMAX     = INPUTS.CIDRAN_MAX ; // -> local var      
  CIDMIN     = INPUTS.CIDRAN_MIN ; 
  NSTORE_ALL = INPUTS.CIDOFF + INPUTS.NGEN_LC + INPUTS.NGENTOT_LC ;
  NSTORE     = INPUTS.NGEN_LC + INPUTS.NGENTOT_LC ; // for this sim-job

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
    r8 = unix_random();   
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

      if ( NTRY == 0 ) {
	
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

	sprintf(c1err,"Unable to try CIDTMP = %d or %d  (CIDRAN=%d)",
		CIDRAN-i2, CIDRAN+i2, CIDRAN );
	sprintf(c2err,"CIDMAX=%d  NPICKRAN=%d  i=%d of %d  CID=%d", 
		CIDMAX, NPICKRAN, i, NSTORE_ALL, GENLC.CID );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      i2++ ;
      if ( CIDADD > 0 ) {
	exec_cidmask(1,CIDADD); // set "USED" bit
	if ( i >= INPUTS.CIDOFF ) 
	  { INPUTS.CIDRAN_LIST[i-INPUTS.CIDOFF] = CIDADD ; }
	USED = 0 ;
      }

    } // end USED block

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

  // Aug 2017: re-init randoms with new SEED that is different
  //           for each sim_SNmix job
  int NEWSEED = INPUTS.ISEED + INPUTS.CIDOFF ;
  srandom(NEWSEED);

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
    ISGRIDONLY = (GENFLAG & OPTMASK_SIMSED_GRIDONLY);
    ISFLAT     = (SIGMA[0] > 1.0E7 );
    if ( ISGRIDONLY && ISFLAT ) {

      bin_user      = SEDMODEL.PARVAL_BIN[ipar_model] ;
      range_user[0] = SEDMODEL.PARVAL_MIN[ipar_model] - 0.5*bin_user ;
      range_user[1] = SEDMODEL.PARVAL_MAX[ipar_model] + 0.5*bin_user;

      INPUTS.GENGAUSS_SIMSED[ipar_user].RANGE[0] = range_user[0] ;
      INPUTS.GENGAUSS_SIMSED[ipar_user].RANGE[1] = range_user[1] ;
    }

    // skip the tests below for baggage or GRIDONLY parameters
    if ( GENFLAG != OPTMASK_SIMSED_PARAM    ) { continue ; }

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
	INPUTS.GENFLAG_SIMSED[N] = OPTMASK_SIMSED_param ;
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
  if ( INPUTS.OPTMASK_SIMSED == OPTMASK_SIMSED_GRIDONLY ) {
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

  int LTRACE=0;
  int istat ;
  char fnam[] = "GENRANGE_CUT" ;

  // ----------- BEGIN ------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
    istat = 1 ; return istat ;
  }

  istat = 0 ;

  if(LTRACE) { printf(" xxx %s: 0 check RA=%f \n", fnam, GENLC.RA); }
  if ( GENLC.RA < INPUTS.GENRANGE_RA[0] ) { return istat; }
  if ( GENLC.RA > INPUTS.GENRANGE_RA[1] ) { return istat; }

  if(LTRACE) { printf(" xxx %s: 1 check RA=%f \n", fnam, GENLC.DEC); }
  if ( GENLC.DEC < INPUTS.GENRANGE_DEC[0] )  { return istat; }
  if ( GENLC.DEC > INPUTS.GENRANGE_DEC[1] )  { return istat; }

  if ( INDEX_GENMODEL != MODEL_LCLIB ) {
    if(LTRACE) {printf(" xxx %s: 2 check zCMB=%f \n",fnam,GENLC.REDSHIFT_CMB);}
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
    PEAKMAG = GENLC.peakmag8_obs[ifilt_obs] ;
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
    ,Trest, Tlast, TGAP, T0GAPOVP[2], flux, fluxerr, peakmag
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
  
  LDMP = ( GENLC.CID < -5555 ) ;

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

    Trest  = (float)GENLC.epoch8_rest[ep] ;
    MJD    = GENLC.MJD[ep];
    MJDDIF = MJD - MJD_last ;

    // require valid epoch to continue

    if ( GENLC.USE_EPOCH[ep] == 0 ) { continue ; }
    
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
    fluxerr = GENLC.flux_errstat[ep] ;
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
  float  peakmag       = (float)GENLC.peakmag8_obs[ifilt_obs] ;
  float  mag_T         = (float)GENLC.genmag8_obs_template[ifilt_obs];

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
  field  = GENLC.FIELDNAME[0]; // get field from header/1st epoch

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

  int imjd, EPMIN, EPMAX, ep, NOBS;
  double flux, flux_err, SNR_CALC, SNR_MEAS, SNR, oldRan, tmpRan ;

  // --------------- BEGIN ----------------

  SEARCHEFF_DATA.CID        = GENLC.CID ;
  SEARCHEFF_DATA.REDSHIFT   = GENLC.REDSHIFT_HELIO ;
  SEARCHEFF_DATA.PEAKMJD    = GENLC.PEAKMJD ;
  SEARCHEFF_DATA.DTPEAK_MIN = GENLC.DTPEAK_MIN ; // closest T-Tpeak
  SEARCHEFF_DATA.SALT2mB    = GENLC.SALT2mB ;
  SEARCHEFF_DATA.SNRMAX     = GENLC.SNRMAX_GLOBAL ;


  sprintf(SEARCHEFF_DATA.FIELDNAME, "%s", GENLC.FIELDNAME[0] );

  NOBS = 0 ;

  for ( imjd = 1; imjd <= GENLC.NEWMJD ; imjd++ ) {

    EPMIN = GENLC.EPOCH_RANGE_NEWMJD[imjd][0] ;
    EPMAX = GENLC.EPOCH_RANGE_NEWMJD[imjd][1] ;

    for ( ep=EPMIN; ep <= EPMAX; ep++ ) {

      if ( GENLC.ISPEAK[ep]     ) { continue ; } // Sep 24 2017
      if ( GENLC.ISTEMPLATE[ep] ) { continue ; }

      SNR_CALC = GENLC.SNR_CALC[ep] ; // Aug 24, 2014

      flux      = GENLC.flux[ep] ;
      flux_err  = GENLC.flux_errstat[ep] ;
      SNR_MEAS  = -9.0 ;
      if ( flux_err > 0.0 ) { SNR_MEAS = flux / flux_err ; }
   
      SNR = SNR_CALC ;
      
      // for SDSS, continue using wrong SNR based on measured flux
      // so that the spec-efficiency function is still correct.
      if ( strcmp(GENLC.SURVEY_NAME,"SDSS") == 0 ) { SNR = SNR_MEAS; }

      SEARCHEFF_DATA.IFILTOBS[NOBS]  = GENLC.IFILT_OBS[ep] ;     
      SEARCHEFF_DATA.MJD[NOBS]       = GENLC.MJD[ep] ;
      SEARCHEFF_DATA.MAG[NOBS]       = GENLC.genmag8_obs[ep] ; 
      SEARCHEFF_DATA.SNR[NOBS]       = SNR ;
      SEARCHEFF_DATA.NPE_SAT[NOBS]   = GENLC.npe_above_sat[ep];

      oldRan = SEARCHEFF_RANDOMS.PIPELINE[NOBS] ;
      if ( oldRan < -0.001 ) 
	{ SEARCHEFF_RANDOMS.PIPELINE[NOBS] = FlatRan1(1);  }

      oldRan = SEARCHEFF_RANDOMS.PHOTPROB[NOBS] ;
      if ( oldRan < -998.0 && INPUTS_SEARCHEFF.NMAP_PHOTPROB >0 )  { 
	if ( INPUTS_SEARCHEFF.NREDUCED_CORR_PHOTPROB == 0 ) 
	  { tmpRan = FlatRan1(1); } // flat randoms [0,1] for uncorrelated
	else
	  { tmpRan = GaussRan(1); } // need Gauss-randoms for correlation
	SEARCHEFF_RANDOMS.PHOTPROB[NOBS] = tmpRan ;
      }
      
      NOBS++ ;
    }
  }
  SEARCHEFF_DATA.NOBS =  NOBS ;


  // load SPEC-EFF randoms and filter-dependent quantities

  int ifilt, ifilt_obs ;

  for ( ifilt=0; ifilt <= MXFILTINDX; ifilt++ ) {

    if ( SEARCHEFF_RANDOMS.SPEC[ifilt] < -0.01 ) 
      { SEARCHEFF_RANDOMS.SPEC[ifilt]  = FlatRan1(1); }

    if ( ifilt == MXFILTINDX ) { continue ; } // avoid array overwrite
    SEARCHEFF_DATA.PEAKMAG[ifilt] = MAG_UNDEFINED ;
    SEARCHEFF_DATA.HOSTMAG[ifilt] = MAG_UNDEFINED ;
    SEARCHEFF_DATA.SBMAG[ifilt]   = MAG_UNDEFINED ;
    
  }


  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      /*
      printf(" xxx ifilt_obs=%d: MAG[PEAK,HOST,SB] = %.3f, %.3f, %.3f \n",
	     ifilt_obs, GENLC.peakmag8_obs[ifilt_obs],
	     SNHOSTGAL.GALMAG[ifilt_obs][0], SNHOSTGAL.SB_MAG[ifilt_obs] );
      */
      SEARCHEFF_DATA.PEAKMAG[ifilt_obs] =  GENLC.peakmag8_obs[ifilt_obs] ;
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
  char fnam[] = "gen_spectype" ;

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
    // xxx mark delete    GENLC.SNTYPE   = (int)INPUTS.FIXMAG[0] ;
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
      GENLC.REDSHIFT_CMB_SMEAR =  
	zhelio_zcmb_translator(SNHOSTGAL.ZSPEC, RA,DEC,eq, +1);

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
      // use zPHOT  here if it is defined
      GENLC.REDSHIFT_FLAG        = REDSHIFT_FLAG_HOSTPHOT ; 
      GENLC.REDSHIFT_HELIO_SMEAR = ZPHOT ;
      GENLC.REDSHIFT_SMEAR_ERR   = ZPHOT_ERR ;
      GENLC.REDSHIFT_CMB_SMEAR = zhelio_zcmb_translator(ZPHOT,RA,DEC,eq, +1);
    }
  }

  return ;

} // end of  setz_unconfirmed


// ***************************************
int gen_smearFlux ( int epoch, int VBOSE ) {


  /**********

   Convert generated  "genmag" into observed flux that
   includes Poisson fluctuations.

   Each genmag is converted into flux in ADU and then into 
   photo-electrons. Error is computed from photo-statistics 
   of signal, sky, readnoise. An effective aperture of 4*pi*PSF^2 
   is used fpr sky- and read-noise.

   The flux is smeared based on fluxcal_err which combines
   photostat error and mag-zeropt error. However, the 
   calculated flux error "flux_err" includes only the
   photostat error and not the mag-zeropt error.


 Mar 1, 2014:
    - Add coherent template noise from SKY and CCD.
    - Use pre-computed randoms from
         GAURAN_Search   = GENLC.RANGauss_NOISE_SEARCH[epoch] ; 
         GAURAN_Template = GENLC.RANGauss_NOISE_TEMPLATE[ifilt_obs] ;
      The simlib 'T: ' option is now ignored -> abort in SIMLIB_READ_DRIVER

    - Use new SIMLIB_GEN array to get SIMLIB conditions 

    - cleanup and variable name changes; remove 'srun' and 'trun' suffix.

 
  Aug 21 2014: compute anomalous host-subtraction noise "sqhostNoise_pe".
               See sim-input key HOSTNOISE_FILE.

  Mar 3 2015: pass SBmag to GEN_NOISEMODEL_HOST().

  Sep 2017: for LCLIB model, subtract template flux. Note that 
            Search-flux can be negative !

  Oct 6 2017: update SNR_CALC after subtracting optional template flux

  Oct 31 2017: if sqsum < 0, set sqsum = small positive number
              (prevents abort)

  Feb 14 2018: minor refactor/cleanup to prepare for FLUXERRMODEL.
               "errtot" is used along with errstat.

  Feb 28 2018: 
   For LCLIB model and SMEARFLAG_ZEROPT, fix bug by subtracting
   flux_T before taking relErr. Also note that flux_T is computed
   earlier than before. Bug found by Gautham.

 Mar 03 2018: GENLC.NOISE_XXX -> FLUXCAL units instead of p.e.

 Mar 05 2018: check option to ignore source & galaxy in reported
              errors (INPUTS.SMEARFLAG_FLUX & 2) 

 Mar 18 2018: compute GENLC.SNR_MON for mag = INPUTS.MAGMONITOR_SNR

  **********************************/

  // define args stripped from GENLC structure

  int ifilt_obs, ifilt, ifield, LDMP=0, GALID, OVP ;

  double 
    genmag, genmag_T
    ,ccdgain, readnoise, skysig, psfsig1, psfsig2, psfratio, zpt, zptsig  
    ,ZPTDIF_pe, Npe_over_FLUXCAL, ZPTDIF_ADU, NADU_over_FLUXCAL, NADU_over_Npe
    ,template_readnoise, template_skysig, template_zpt
    ,zsn, psfsig_arcsec, psfFWHM_arcsec, galmag, GALMAG, SBmag
    ;

  // Note three types of errors:
  //   errS   = search error only  (source+sky+host)
  //   errSZ  = search error + ZPsmear
  //   errSZT = search error + ZPsmear + template error
  
  // local variables   
  double 
    arg, mjd, fluxObs_adu
    ,fluxsn_pe,  fluxsn_pe_err, fluxgal_pe, fluxmon_pe
    ,fluxsn_adu, fluxsn_adu_errS, fluxsn_adu_errSZ, fluxsn_adu_errSZT
    ,fluxsn_adu_errST, fluxsn_adu_errZ, flux_T 
    ,sqskyerr_pe, template_sqskyerr_pe, template_sqerr_pe
    ,sqccderr_pe, template_sqccderr_pe
    ,sqImageNoise_pe, fluxsn_adu_errData
    ,template_pe_err, template_adu_err, fluxsn_adu_errReal
    ,sqadderr_pe, area_bg
    ,mag_smear, magerr_tmp, fluxErr_tmp, sqerr_tmp, skysig_tmp_pe
    ,sqsum, sqerr, sqerr1, sqerr2, sqerr_ran
    ,err1, err2, relerr, errtmp, errtot, errstat, pixsize
    ,xt, zptfac, crazyflux, Trest, Tobs
    ,fluxerrCor, SNR_CALC, SNR_MON, SNR
    ,GAURAN_Search, GAURAN_ZP, GAURAN_Template
    ,snsep, HOSTNOISE_FLUXCAL, HOSTNOISE_ERRSCALE, HOSTNOISE_pe
    ,noisePar[10]
    ;

  double scale_fluxErr = 1.0 ;
  char field[MXCHAR_FIELDNAME], band[4];
  char fnam[] = "gen_smearFlux" ;

  // ----------------- BEGIN --------------

  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) { return SUCCESS; }

  // strip GENLC info into local variables.

  ifilt_obs      = GENLC.IFILT_OBS[epoch] ;
  ifilt = GENLC.IFILTINVMAP_OBS[ifilt_obs] ; // sparse index
  sprintf(band, "%c", FILTERSTRING[ifilt_obs] );

  genmag         = GENLC.genmag8_obs[epoch];
  mag_smear      = GENLC.magsmear8[epoch];
  zsn            = GENLC.REDSHIFT_HELIO ;
  sprintf(field,"%s", GENLC.FIELDNAME[epoch] );
  ifield         = IFIELD_OVP_SIMLIB(1,field);
  if ( ifield < 0 ) { ifield=0; } // Nov 2016: needed for GAURAN_TEMPLATE

  mjd        = SIMLIB_OBS_GEN.MJD[epoch] ;
  pixsize    = SIMLIB_OBS_GEN.PIXSIZE[epoch] ;
  ccdgain    = SIMLIB_OBS_GEN.CCDGAIN[epoch] ;
  skysig     = SIMLIB_OBS_GEN.SKYSIG[epoch] ;
  readnoise  = SIMLIB_OBS_GEN.READNOISE[epoch] ;
  psfsig1    = SIMLIB_OBS_GEN.PSFSIG1[epoch] ; // pixels
  psfsig2    = SIMLIB_OBS_GEN.PSFSIG2[epoch] ;
  psfratio   = SIMLIB_OBS_GEN.PSFRATIO[epoch] ;
  zpt        = SIMLIB_OBS_GEN.ZPTADU[epoch] ;
  zptsig     = SIMLIB_OBS_GEN.ZPTSIG[epoch] ;
  
  GALID   = SNHOSTGAL.GALID ;
  GALMAG  = SNHOSTGAL.GALMAG[ifilt_obs][0] ; // full galaxy mag
  SBmag   = SNHOSTGAL.SB_MAG[ifilt_obs]  ;   // in 1 sq arcsec
  if ( SBmag > 32.0 ) { SBmag = 32.0; }      // to limit fluxerrmap size
  snsep   = SNHOSTGAL.SNSEP ;                // SN-host sep, arcsec

  template_skysig     = SIMLIB_OBS_GEN.TEMPLATE_SKYSIG[epoch] ; 
  template_readnoise  = SIMLIB_OBS_GEN.TEMPLATE_READNOISE[epoch] ; 
  template_zpt        = SIMLIB_OBS_GEN.TEMPLATE_ZPT[epoch] ;

  genmag_T = GENLC.genmag8_obs_template[ifilt_obs]; // Sep 2017/LCLIB
  flux_T   = 0.0 ;
  if ( genmag_T < 90.0 ) {
    arg     = 0.4 * ( zpt - genmag_T );
    flux_T  = pow(10.0,arg);        // flux in ADUs    
  }

  xt = INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] ; // for error check only
  if ( xt < 1.0 ) { xt = 1.0 ; }

  // compute zp in photo-electrons in case we need to 
  // convert a FLUXCAL quantity back into pe.
  ZPTDIF_ADU = zpt - ZEROPOINT_FLUXCAL_DEFAULT ;
  NADU_over_FLUXCAL = pow( TEN , 0.4*ZPTDIF_ADU) ;
  Npe_over_FLUXCAL  = NADU_over_FLUXCAL * ccdgain;
  NADU_over_Npe     = NADU_over_FLUXCAL/Npe_over_FLUXCAL ;  

  // init output
  GENLC.flux[epoch]         = NULLFLOAT ; 
  GENLC.flux_errstat[epoch] = NULLFLOAT ;     
  GENLC.flux_errtot[epoch]  = NULLFLOAT ;     

  // bail on bad input or non-existant peak epoch
  int SKIPIT = 0;
  if ( GENLC.ISPEAK[epoch]      )  { SKIPIT = 1 ; }
  if ( GENLC.ISTEMPLATE[epoch]  )  { SKIPIT = 1 ; }
  if ( zpt     < 10.0   )          { SKIPIT = 1 ; }
  if ( psfsig1 < 0.0001 )          { SKIPIT = 1 ; }
  if ( skysig  < 0.0001 )          { SKIPIT = 1 ; }
  if ( SKIPIT ) {    return ERROR ;  }


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

  // ################ GET FLUX STAT-ERROR ##################

  // get effective aperture  area (pixels)
  area_bg = NoiseEquivAperture(psfsig1, psfsig2, psfratio );
  psfsig_arcsec = pixsize * sqrt(area_bg/(2.0*TWOPI))  ; 
  psfFWHM_arcsec = 2.3548 * psfsig_arcsec ;

  // get total sky noise for search run; includes sky & CCD noise
  skysig_tmp_pe = skysig * ccdgain ;  // convert ADU -> pe
  sqskyerr_pe = area_bg * (skysig_tmp_pe*skysig_tmp_pe);
  sqccderr_pe = area_bg * (readnoise*readnoise) ;  // in Npe


  // Apri 26 2016: optional non-linearity bias to signal (NONLINEARITY_FILE)
  double scale_nonLin = GET_NONLIN(band, fluxsn_pe, sqskyerr_pe, genmag);
  if ( fabs(scale_nonLin-1.0) > 1.0E-7 ) {
    fluxsn_pe   *= scale_nonLin ;
    fluxsn_adu  *= scale_nonLin ;
    sqskyerr_pe *= scale_nonLin ;
    sqccderr_pe *= scale_nonLin ;
  }

  // add sky-noise from template, integrated over effective aperture 
  template_pe_err = template_sqskyerr_pe = template_sqccderr_pe = 0.0 ;
  template_sqerr_pe = 0.0 ;

  if ( SIMLIB_TEMPLATE.USEFLAG && template_skysig > 0.0 ) {
    if ( template_zpt < 10.0 ) {
      sprintf(c1err,"Invalid template_zpt(%c)=%f for  LIBID=%d at MJD=%.3f", 
	      FILTERSTRING[ifilt_obs], template_zpt, GENLC.SIMLIB_ID, mjd );
      sprintf(c2err,"Need TEMPLATE_ZPT to scale template noise.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    skysig_tmp_pe = template_skysig * ccdgain ;  // convert sigma in ADU -> pe.
    template_sqskyerr_pe = area_bg * (skysig_tmp_pe*skysig_tmp_pe) ;
    template_sqccderr_pe = area_bg * (template_readnoise*template_readnoise);

    // scale template noise to the search image
    double zparg, zfac;
    zparg = 0.8*(zpt - template_zpt); // Feb 2 2017 bug fix
    zfac  = pow(TEN, zparg);
    template_sqskyerr_pe *= zfac ;
    template_sqccderr_pe *= zfac ;

    template_sqerr_pe = template_sqccderr_pe + template_sqskyerr_pe ;
    template_pe_err   = sqrt(template_sqerr_pe) ;
  }


  // @@@@@@@@@@@@ LEGACY ERRFDUGE @@@@@@@@@@@@@@
  // Feb 2012: add FLUXERR_ADD in quadrature
  // Jan 2014: require EXPOSURE_TIME=1 to implement; avoid crazy error
  double  ERR_CAL, ERR_pe, XT ;
  ERR_CAL  = GENLC.SIMLIB_FLUXERR_ADDPAR[ifilt_obs];
  XT       = INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] ; // Jan 2014
  sqadderr_pe = 0.0 ;
  if ( SIMLIB_FLUXERR_COR.USE  && ERR_CAL >1.E-6  &&   XT == 1.0 ) {
    // translate ERR_CAL from FLUXCAL back to p.e.
    ERR_pe      = ERR_CAL * Npe_over_FLUXCAL ;
    sqadderr_pe = (ERR_pe * ERR_pe) ;
  }
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



  // -------------------------------------------
  // galaxy noise from photo-stats
  fluxgal_pe = 0.0;

  OVP = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_PHOT ;
  if ( OVP > 0 ) {
    // get galmag and SN-host sep:  4*PI*PSF^2 = area
    galmag        = interp_GALMAG_HOSTLIB(ifilt_obs, psfsig_arcsec );
    arg           = 0.4 * ( zpt - galmag );
    fluxgal_pe    = ccdgain * pow(10.0,arg);   // effec-aper flux in pe.
  }

  // @@@@@@@@@@@@@@@@@ LEGACY ERRFUDGE @@@@@@@@@@@@@@@@@@@@@@@@@
  // Aug 2014: anomolous host-subtraction noise (HOSTNOISE_FILE)
  //           Note the dependence on both band and field.
  OVP = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_IMAGE ;
  sqImageNoise_pe = 0.0 ;
  HOSTNOISE_pe    = HOSTNOISE_FLUXCAL = 0.0 ;
  if ( OVP ) {    
    GEN_NOISEMODEL_HOST_LEGACY(band,field,GALID,GALMAG,SBmag,snsep, 
			       noisePar);  // <== return this array

    HOSTNOISE_FLUXCAL  = noisePar[0] ; // add this noise, per pixel
    HOSTNOISE_ERRSCALE = noisePar[1] ; // scale on added sky noise

    // convert extra noise in FLUXCAL unit to Npe (per pixel)
    HOSTNOISE_pe    = HOSTNOISE_FLUXCAL * Npe_over_FLUXCAL ;

    double SQ0 = HOSTNOISE_pe * HOSTNOISE_pe ;
    double SQ1 = HOSTNOISE_ERRSCALE * HOSTNOISE_ERRSCALE ;

    sqImageNoise_pe  = 
      (area_bg * SQ0)  +                       // quadrature model
      (sqskyerr_pe+fluxgal_pe)*(SQ1-1.0)  ;    // err-scale;
  }
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // --------------------
  // add up SN flux error (photo-electrons^2) in quadrature
  // Do not include correlated template noise here.
  sqsum 
    = fluxsn_pe        // signal stat-error
    + sqskyerr_pe      // sky-err from search run
    + sqccderr_pe      // CCD read noise (added Dec 13, 2010)
    + sqadderr_pe      // from FLUXERR_ADD key in simlib
    + fluxgal_pe       // square of host galaxy stat-error
    ;

  // net search stat error in photoelectrons (no template)
  fluxsn_pe_err = sqrt ( sqsum ) ;  

  // -------------------------------------------
  // compute calculated SNR to use in lookup maps;
  // include coh template noise, but ignore anomalous image noise.
  SNR_CALC      = fluxsn_pe / sqrt(sqsum + template_sqskyerr_pe );
 
  if ( INPUTS.MAGMONITOR_SNR > 10 ) {
    double sqsum_mon = sqsum - fluxsn_pe + fluxmon_pe ;
    SNR_MON = fluxmon_pe / sqrt(sqsum_mon + template_sqskyerr_pe );
  }
  else
    { SNR_MON = 0.0 ; }

  // @@@@@@@@@@@@@@@@@ LEGACY ERRFUDGE @@@@@@@@@@@@@@@@@@@@@@@@@
  // check optional correction from simlib map
  fluxerrCor     = get_SIMLIB_fluxerrScale_LEGACY(ifilt_obs, SNR_CALC );    
  fluxsn_pe_err *= fluxerrCor;
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  // Nov 2016: apply user-optional magErr fudge in quadrature (default=0)
  magerr_tmp = INPUTS.FUDGE_MAGERR_FILTER[ifilt_obs]; 
  if ( magerr_tmp > 1.0E-9 ) {
    fluxErr_tmp   = fluxsn_pe * ( 1.0 - pow(10.0,-0.4*magerr_tmp)) ;
    sqerr_tmp     = fluxErr_tmp * fluxErr_tmp ;
    fluxsn_pe_err = sqrt(fluxsn_pe_err*fluxsn_pe_err + sqerr_tmp);
  }
  // Optional errscale fudge (default=1) applied to true and reported errors.
  fluxsn_pe_err   *= INPUTS.FUDGESCALE_FLUXERR_FILTER[ifilt_obs] ;
  template_pe_err *= INPUTS.FUDGESCALE_FLUXERR_FILTER[ifilt_obs] ;

  // -----------------------
  // get search error in ADU
  fluxsn_adu_errS   = fluxsn_pe_err   / ccdgain ;
  template_adu_err  = template_pe_err / ccdgain ;

  // store individual noise contributions in FLUXCAL units
  GENLC.NOISE_SN[epoch]        = sqrt(fluxsn_pe ) /Npe_over_FLUXCAL ;
  GENLC.NOISE_SKY[epoch]       = sqrt(sqskyerr_pe)/Npe_over_FLUXCAL ;
  GENLC.NOISE_TEMPLATE[epoch]  = sqrt(template_sqskyerr_pe)/Npe_over_FLUXCAL ;
  GENLC.NOISE_CCD[epoch]       = sqrt(sqccderr_pe)/Npe_over_FLUXCAL ;
  GENLC.NOISE_HOSTGAL_PHOT[epoch]  = sqrt(fluxgal_pe)/Npe_over_FLUXCAL ;
  GENLC.NOISE_HOSTGAL_IMAGE[epoch] = sqrt(sqImageNoise_pe)/Npe_over_FLUXCAL ;
  GENLC.NOISE_TOTAL[epoch]         = fluxsn_pe_err/Npe_over_FLUXCAL ;


  // ############### ADD ZEROPOINT ERROR ###################
  // note that zeropoint smearing applies only to fluxcal, and NOT 
  // to fluxsn_adu. This smearing is random for each epoch.
  // NOTE: maybe the zpt error should be correlated among epochs ?

  if ( INPUTS.SMEARFLAG_ZEROPT > 0 ) {
    relerr  = pow(TEN, 0.4*zptsig) - 1.0 ;
    err1    = fluxsn_adu_errS  ;    
    err2    = (fluxsn_adu-flux_T) * relerr ;    
    sqerr1  = err1 * err1 ;    sqerr2  = err2 * err2 ;
    fluxsn_adu_errSZ  = sqrt ( sqerr1 + sqerr2 );
    fluxsn_adu_errZ   = err2; 
  }
  else {
    sqerr1 = sqerr2 = fluxsn_adu_errZ = 0.0 ;
    fluxsn_adu_errSZ = fluxsn_adu_errS ;  // photostat only
  }

  // ############### SMEAR OBSERVED FLUX #####################
  // Now we have nominal flux in "template" ADUs, and its error.
  // Apply gaussian fluctuation to get "observed" fluxObs


  // get Gaussian randoms separately for SEARCH and TEMPLATE.
  // All SEARCH randoms are uncorrelated, but TEMPLATE randoms
  // depend on field/band, but are the same for each epoch.
  if ( INPUTS.SMEARFLAG_FLUX > 0 )  {     
    GAURAN_Search   = GENLC.RANGauss_NOISE_SEARCH[epoch] ; 
    GAURAN_ZP       = GENLC.RANGauss_NOISE_ZP[epoch] ; 
    GAURAN_Template = GENLC.RANGauss_NOISE_TEMPLATE[ifield][ifilt_obs] ;
  }
  else {  
    GAURAN_Search   = 0.0 ; 
    GAURAN_ZP       = 0.0 ; 
    GAURAN_Template = 0.0 ;
  }


  // Feb 2018: fudge error from FLUXERRMODEL. Should replace _legacy codes.
  if ( NMAP_FLUXERRMODEL > 0 ) {
    double ERRPARLIST[MXPAR_FLUXERRMAP], fluxerr_sim, fluxerr_data ;
    double LOGSNR = log10(SNR_CALC);    
    int OPT = 0;
    if ( LOGSNR < -0.9 ) { LOGSNR = -0.9 ; }
    ERRPARLIST[IPAR_FLUXERRMAP_MJD]    = mjd;
    ERRPARLIST[IPAR_FLUXERRMAP_PSF]    = psfFWHM_arcsec;  // FWHM, arcsec
    ERRPARLIST[IPAR_FLUXERRMAP_SKYSIG] = skysig;         // ADU/pixel
    ERRPARLIST[IPAR_FLUXERRMAP_ZP]     = zpt;            // observed ZP, ADU
    ERRPARLIST[IPAR_FLUXERRMAP_LOGSNR] = LOGSNR ;
    ERRPARLIST[IPAR_FLUXERRMAP_SBMAG]  = SBmag ;
    ERRPARLIST[IPAR_FLUXERRMAP_GALMAG] = GALMAG ;
    ERRPARLIST[IPAR_FLUXERRMAP_SNSEP]  = snsep ;
    
    // pass FLUXCAL units to fluxErrModel in case of additive term.
    double FLUXCALERR_in = fluxsn_adu_errS/NADU_over_FLUXCAL ;
    double FLUXCALERR_REAL, FLUXCALERR_DATA ;

    get_FLUXERRMODEL(OPT, FLUXCALERR_in, band, field,      // (I)
		     NPAR_FLUXERRMAP_REQUIRE, ERRPARLIST,  // (I)
		     &FLUXCALERR_REAL, &FLUXCALERR_DATA) ; // (O)

    // convert FLUXCAL back to pe
    fluxsn_adu_errReal = FLUXCALERR_REAL * NADU_over_FLUXCAL ;
    fluxsn_adu_errData = FLUXCALERR_DATA * NADU_over_FLUXCAL ;
    scale_fluxErr      = fluxsn_adu_errData/fluxsn_adu_errS ;
  }
  else {
    // leagcy option
    sqsum = 
      (fluxsn_adu_errS * fluxsn_adu_errS) + 
      (sqImageNoise_pe / (ccdgain*ccdgain)) ;

    fluxsn_adu_errReal = sqrt(sqsum);
  }

  fluxObs_adu  = fluxsn_adu 
    + ( fluxsn_adu_errReal * GAURAN_Search   ) 
    + ( fluxsn_adu_errZ    * GAURAN_ZP       ) // Feb 2018
    + ( template_adu_err   * GAURAN_Template ) // coherent template err
    ;


  OVP = (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_RANDOM_TEMPLATENOISE);
  if ( OVP > 0 ) {
    // treat template noise as random instead of coherent
    // note that fluxObs_adu is overwritten here so that
    // the coherent-template fluctuation is dropped.
    sqerr  = fluxsn_adu_errSZ*fluxsn_adu_errSZ ;
    sqerr += (template_adu_err * template_adu_err);
    errtmp = sqrt(sqerr);
    fluxObs_adu  = fluxsn_adu + ( errtmp * GAURAN_Search )  ;
  }


  // use 1/z^2 dependence on mag to set bounds for crazy flux abort;
  // account for exposure time (xt) and SIMLIB zeropoint (zptfac)
  // Also add 10 sigma of noise to allow for fluctuations.

  arg        = 0.4 * ( zpt - 31.0 );  
  zptfac     = pow(10.0,arg);  
  if ( zsn > 1.0E-9 ) 
    { crazyflux  = 2.*(1.E4 * zptfac * xt) / (zsn*zsn) ; }
  else
    { crazyflux = 1.0E14 ; } // for LCLIB (July 2018)

  crazyflux += (10.*fluxsn_adu_errS) ;

  if ( mag_smear < 0 ) 
    { crazyflux *= pow(TEN,-0.4*mag_smear); }

  if ( GENLC.FUDGE_SNRMAX_FLAG == 2 && INPUTS.FUDGE_SNRMAX > 1.0 ) 
    { crazyflux *= INPUTS.FUDGE_SNRMAX; }

  if ( INDEX_GENMODEL == MODEL_SIMSED ) 
    { crazyflux *= 10.0; }    // allow for really bright objects (Aug 2017)

  if ( INDEX_GENMODEL == MODEL_LCLIB ) 
    { crazyflux *= 100.0; }  

  
  if ( fluxObs_adu > crazyflux || fluxObs_adu < -xt*1.0E9 ) {

    printf("\n\n PRE-ABORT COMMENTS for CID=%d: \n", GENLC.CID );
    printf("\t LIBID = %d   MJD=%f  zHel=%6.4f  MU=%7.3f \n",
	   GENLC.SIMLIB_ID, mjd,
	   GENLC.REDSHIFT_HELIO, GENLC.DLMU );

    Trest = GENLC.epoch8_rest[epoch]; 
    Tobs  = Trest * ( 1.0 + GENLC.REDSHIFT_HELIO );
    printf("\t Trest=%6.2f  Tobs=%.2f  genmag(%c)=%6.1f  \n", 
	   Trest, Tobs, FILTERSTRING[ifilt_obs], genmag );
     
    printf("\t ccdgain=%6.2f \n",  ccdgain ); 
    printf("\t zpt(srun) = %f   zptfac=%f \n", 
	   zpt, zptfac );
    printf("\t GAURAN(photostat)=%f  magsmear=%6.2f\n", 
	   GAURAN_Search, mag_smear  );
    printf("\t sqerr[1,2] = %9.3le , %9.3le \n", sqerr1, sqerr2 );
    printf("\t ERR_CAL = %f \n", ERR_CAL );
    printf("\t SKYerr_pe(srun,trun) = %f, %f \n",
	  sqrt(sqskyerr_pe), sqrt(template_sqskyerr_pe) );
    printf("\t CCDerr_pe(srun,trun) = %f, %f \n",
	  sqrt(sqccderr_pe), sqrt(template_sqccderr_pe) );
    printf("\t err_pe(flux,host,image) = %f, %f, %f \n",
	   sqrt(fluxsn_pe), sqrt(fluxgal_pe), sqrt(sqImageNoise_pe) );

    printf("\t HOSTNOISE(FLUXCAL,pe) = %f, %f \n",
	   HOSTNOISE_FLUXCAL, HOSTNOISE_pe);

    printf("\t zptsig=%f   \n",  zptsig);
    printf("\t GEN(AV,RV) = %7.3f , %7.3f  SHAPEPAR=%7.3f  (c=%7.3f)\n", 
	   GENLC.AV, GENLC.RV, GENLC.SHAPEPAR, GENLC.SALT2c );
    printf("\t Gauss random number: %f \n", 
	   GENLC.GENSMEAR_RANGauss_FILTER[0] );

    if ( GENFRAME_OPT  == GENFRAME_REST ) {
      printf("\t Kcor  %s = %le   AVwarp=%7.3f\n"
	     , GENLC.kcornam[epoch]
	     , GENLC.kcorval8[epoch] 
	     , GENLC.AVwarp8[epoch] 
	     );
    }
    if ( INDEX_GENMODEL  == MODEL_SIMSED ) {
      printf("\t SIMSED PARAMS %s,%s = %f, %f  (x0=%le)\n"
	     ,INPUTS.PARNAME_SIMSED[1]
	     ,INPUTS.PARNAME_SIMSED[2]
	     ,GENLC.SIMSED_PARVAL[1]
	     ,GENLC.SIMSED_PARVAL[2], GENLC.SALT2x0 );
    }

    if ( INDEX_GENMODEL == MODEL_LCLIB ) {
      printf("\t LCLIB EVENT ID = %lld \n", LCLIB_EVENT.ID);
    }

    printf("\t CRAZYFLUX = %9.3le\n", crazyflux);

    sprintf(c1err, "Too large %c(%d)-flux = %9.3le  fluxsn_adu=%9.3le", 
	    FILTERSTRING[ifilt_obs], ifilt_obs, fluxObs_adu, fluxsn_adu );
    sprintf(c2err, "fluxsn_adu_errS=%9.3le  fluxsn_adu_errSZ=%9.3le ", 
	    fluxsn_adu_errS, fluxsn_adu_errSZ );

    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    // errmsg(SEV_WARN, 0, fnam, c1err, c2err) ; 

  }


  // Adjust reported error to be based on observed flux instead of
  // true flux.
  if ( fluxObs_adu > 0 ) 
    { sqerr_ran = (fluxObs_adu - fluxsn_adu)/ccdgain; }
  else
    { sqerr_ran = -fluxsn_adu/ccdgain; }

  sqerr1 = (fluxsn_adu_errSZ * fluxsn_adu_errSZ) ;
  fluxsn_adu_errSZ = sqrt(sqerr1 + sqerr_ran);

  // Add template error to reported error (SZT).
  sqerr1 = (fluxsn_adu_errSZ * fluxsn_adu_errSZ) ;
  sqerr2 = (template_adu_err * template_adu_err) ; // add coh template error
  sqsum  = sqerr1 + sqerr2 ;
  if ( sqsum < 1.0E-9 ) { sqsum = 1.01E-9 ; }

  if ( sqsum < 1.0E-9 ) {
    printf(" \n PRE-ABORT DUMP: \n");
    printf("   sqerr1=%f  sqerr2=%f  sqerr_ran=%f\n",	   
	   sqerr1, sqerr2, sqerr_ran);
    printf("   GAURAN(Search,Template) = %f, %f \n",
	   GAURAN_Search, GAURAN_Template );
    printf("   fluxADU(obs,true) = %f, %f \n",
	   fluxObs_adu, fluxsn_adu);
    printf("   Trest=%6.2f  genmag=%6.2f  zHEL=%.4f \n", 
	   GENLC.epoch8_rest[epoch], genmag, GENLC.REDSHIFT_HELIO );
    sprintf(c1err,"Invalid sqerr=%le for fluxsn_adu_errSZT", sqsum);
    sprintf(c2err,"LIBID=%d  MJD=%.3f  band=%c",
	    GENLC.SIMLIB_ID, mjd, FILTERSTRING[ifilt_obs]);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;
  }
  fluxsn_adu_errSZT = sqrt(sqsum);



  // Compute reported without ZP error
  sqerr1  = fluxsn_adu_errS  * fluxsn_adu_errS ;
  sqerr2  = template_adu_err * template_adu_err ;
  sqsum   = sqerr1 + sqerr2 + sqerr_ran ;
  fluxsn_adu_errST = sqrt(sqsum);   // S & T, but no ZP error

  // store errors for GENLC storage
  errstat = fluxsn_adu_errST ;   // reported total error without ZP err
  errtot  = fluxsn_adu_errSZT  ; // includes ZP error


  // @@@@@@@@@@@@@@@@@ LEGACY ERRFUDGE @@@@@@@@@@@@@@@@@@@@@@@@@
  // Scale reported (not true) errors to better match data
  if ( INPUTS.FUDGEOPT_FLUXERR > 0 ) {
    scale_fluxErr = 
      scale_fluxErrModel_legacy(band,field,mjd,zpt,skysig,psfsig1);
  }
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  errtot  *= scale_fluxErr ;
  errstat *= scale_fluxErr ;      

  // check fudge on measured flux-errors (does not affect true errors)
  errtot    *= INPUTS.FUDGESCALE_FLUXERR2_FILTER[ifilt_obs] ;
  errstat   *= INPUTS.FUDGESCALE_FLUXERR2_FILTER[ifilt_obs] ; 

 
  // check option to ignore source & host error in reported error
  // (SMP-like)
  if ( (INPUTS.SMEARFLAG_FLUX & 2) > 0 ) {
    double sqerr_pe, err_remove, err_orig=errstat;
    sqerr_pe   = fluxsn_pe + fluxgal_pe ;
    err_remove = sqrt(sqerr_pe) * NADU_over_Npe ;  
    errstat    = sqrt(errstat*errstat - err_remove*err_remove);
    errtot     = sqrt(errtot *errtot  - err_remove*err_remove);

    /*
    if ( SBmag < 22. ) {
      printf(" xxx err = %.1f -> %.1f  : err_remove=%.1f \n",
	     err_orig, errstat, err_remove );
    }
    */
  }

  if ( genmag > 600.0 ) 
    { fluxObs_adu = errtot = errstat = template_adu_err = -9.0 ; }


  // --------------------------------------------
  // Check optional template flux to subtract (for LCLIB model).
  // Beware that coherent template fluctuations are not included,
  // so deep templates are assumed.
  // This template-flux subtraction is done at the very end so that
  // search-soure noise is included in the Poisson noise.
  if ( genmag_T < 90.0 ) {
    fluxObs_adu       -= flux_T ;  // can be pos or neg
    fluxsn_adu        -= flux_T ;  // true flux without fluctuations
    SNR_CALC *= ( fluxsn_adu / ( fluxsn_adu + flux_T) ) ;
  }

  // - - - - - - - - - - - - - - - - - - - 
  // Jan 2018: check for saturation. NPE > 0 --> saturation
  int npe_above_sat = npe_above_saturation(epoch,fluxsn_pe+fluxgal_pe );
  GENLC.npe_above_sat[epoch] = npe_above_sat ;
  if ( npe_above_sat > 0 ) {
    fluxObs_adu = 0.0 ;
    errstat = errtot = FLUXCALERR_SATURATE ;
  }
  // - - - - - - - - - - - - - - - - - - - 

  // load global GENLC array.
  GENLC.flux[epoch]         = fluxObs_adu ;        // flux with fluctuation

  // store stat error for the flux in ADU;  the FLUXCAL and FLUXCAL_ERR 
  // (computed later in fluxcal_SNDATA) include the ZPTSIG term.
  GENLC.flux_errstat[epoch] = errstat ;    // without ZP smear

  // separately store total error that include ZPTSIG term.
  GENLC.flux_errtot[epoch]  = errtot ; 

  
  // keep track of epoch with max SNR (Jun 2018)
  SNR = fluxObs_adu/errstat ;
  if (SNR > GENLC.SNRMAX_GLOBAL) 
    { GENLC.SNRMAX_GLOBAL = SNR;  GENLC.IEPOCH_SNRMAX = epoch;  }

  // store generated/true SNR without fluctuations
  GENLC.trueSNR[epoch] =  fluxsn_adu/fluxsn_adu_errSZT  ;

  // Aug 24 2014: store ideal/calculated SNR
  //  Beware that fluxErrModel scale is not included here ; see below.
  GENLC.SNR_CALC[epoch] = SNR_CALC ;

  // Mar 18 2018: store SNR of fixed monitor mag
  GENLC.SNR_MON[epoch]  = SNR_MON ;

  // store coherent template error.
  GENLC.template_err[epoch] = template_adu_err ;

  // debug dump
  //  LDMP = ( GENLC.SIMLIB_ID == 121 && ifilt_obs==5  ) ;
  if ( LDMP ) {

    printf("\n # ------------------------------------------ \n");
    printf(" XXX NOISE DUMP  \n");
    printf(" XXX CID=%d at MJD=%f  ifilt_obs=%d  LIBID=%d\n",
	   GENLC.CID, mjd, ifilt_obs, GENLC.SIMLIB_ID );
 
    printf(" flux(obs,true) = %f, %f  ADU\n", fluxObs_adu, fluxsn_adu);
    printf(" err(S,SZ,SZT) = %f, %f, %f \n",
	   fluxsn_adu_errS, fluxsn_adu_errSZ, fluxsn_adu_errSZT);
    printf(" sqerr_ran = %f   GAIN=%f\n", sqerr_ran, ccdgain ) ;

    printf(" FLUX = %f +- %f  ADU \n", 
	   GENLC.flux[epoch], GENLC.flux_errstat[epoch] );

    printf(" psfsig[1,2] = %f, %f pixels  2/1 = %f   \n",
	   psfsig1, psfsig2, psfratio );
    printf(" Effective area = %f pixels \n", area_bg );
    printf(" skysig(S,T) = %.3f , %.3f (ADU/pix) \n", 
	   skysig, template_skysig) ;
    printf(" zpt(S,T) = %.3f , %3f \n", zpt, template_zpt);
    printf(" scale_fluxErr=%.3f   FUDGEOPT_FLUXERR=%d \n",
	   scale_fluxErr, INPUTS.FUDGEOPT_FLUXERR );

    printf(" Gain=%f \n", ccdgain);
    printf(" GUARAN = %f, %f \n", GAURAN_Search, GAURAN_Template );
    printf(" Noise(src)  = %f  p.e. \n",   
	   GENLC.NOISE_SN[epoch]*Npe_over_FLUXCAL );
    printf(" Noise(sky)  = %f  p.e. \n",  
	   GENLC.NOISE_SKY[epoch]*Npe_over_FLUXCAL );
    printf(" Noise(tmpl) = %f  p.e.  USEFLAG=%d \n",  
	   GENLC.NOISE_TEMPLATE[epoch]*Npe_over_FLUXCAL, 
	   SIMLIB_TEMPLATE.USEFLAG );
    printf(" Noise(host-phot) = %f  p.e.   (GALMAG=%.2f, SBmag=%.2f)\n",   
	   GENLC.NOISE_HOSTGAL_PHOT[epoch]*Npe_over_FLUXCAL, GALMAG, SBmag );
    printf(" Noise(host-image)= %f  p.e. \n",   
	   GENLC.NOISE_HOSTGAL_IMAGE[epoch]*Npe_over_FLUXCAL );
    printf(" Noise(CCD sky+tmpl) = %f p.e. \n", 
	   GENLC.NOISE_CCD[epoch]*Npe_over_FLUXCAL  );
    printf(" Noise(tot)  = %f  p.e. \n", 
	   GENLC.NOISE_TOTAL[epoch]*Npe_over_FLUXCAL );
    printf(" fluxsn_adu_errS = %f ADU \n",  fluxsn_adu_errS );
    printf(" errtot=%.3f  errstat=%.3f \n", errtot, errstat );
    debugexit(fnam); 
  }

  return SUCCESS ;

}  // end of gen_smearFlux


// **********************************
int gen_smearMag ( int epoch, int VBOSE) {

  // Dec 2011
  // Convert flux(ADU) and ZP into observed mag and error.

  double flux, flux_errtot, flux_errstat, flux_tmp, zpt ;
  double mag, mag_err , mag_tmp, genmag ;

  //  char fnam[] = "gen_smearMag" ;

  // -------------- BEGIN --------------

  flux         = GENLC.flux[epoch];
  flux_errtot  = GENLC.flux_errtot[epoch];
  flux_errstat = GENLC.flux_errstat[epoch];
  genmag       = GENLC.genmag8_obs[epoch];
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
    //    GENLC.genmag8_obs[epoch]  = MAG_SATURATE ; // xxx delete Sep 3 2018
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
    double genmag = GENLC.genmag8_obs[epoch];
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
  // Feb 23, 2012: store search-run ZPT instead of ZPT for template run
  //
  // Aug 23, 2012: write file name with %8.8 of NGEN > 1 million
  //               or if the CIDRAN flag is set.
  //               Goes with new MXGEN_SIM = 3 million.
  //
  // Mar 3 2014: SMEARFLAG_ZEROPT=2 --> report total error.
  // Nov 23 2014: set PHOTFLAG
  // Sep 04 2016: SNDATA.MAG computed in fluxcal_SNDATA().
  // Apr 17 2017: SKY_SIG back to ADU (not FLUXCAL)
  // Aug 10 2017: abort of NEPOCH>=MXEPOCH to protect SNDATA arrays.
  // Aug 26 2017: ctmp[60] -> ctmp[200] to handle long model names
  // Sep 08 2017: load LCLIB info
  // Sep 20 2017: fix refactor bug setting SURVEY and SUBSURVEY
  // Mar 02 2018: scale SKYSIG_T to same ZP as SKYSIG_S
  //              Only affects output to data files.
  // Jun 01 2018: for LCLIB, set all redshifts and errors to zero
  //
  // Jul 21 2018: apply MW correction on fluxes (default is no cor)
  //
  // Dec 10 2018: load BYOSED info
  // ---------------------------------------

  int PHOTFLAG_DETECT  = INPUTS_SEARCHEFF.PHOTFLAG_DETECT;
  int PHOTFLAG_TRIGGER = INPUTS_SEARCHEFF.PHOTFLAG_TRIGGER ;

  int epoch, ifilt, ifilt_obs, NFILT, NPAR, ipar, MSKTMP, istat, zFLAG ;
  double ZP_S, ZP_T, SKYSIG_T_scale, arg;
  double MCOR_MAP_MW, MCOR_TRUE_MW ;
  
  char ctmp[200], ccid[12];
  char *cptr, *tmpName;
  char fnam[] = "snlc_to_SNDATA" ;

  // --------------- BEGIN -------------

  // always start with header info
  SNDATA.SIMLIB_MSKOPT  = INPUTS.SIMLIB_MSKOPT ;
  sprintf(SNDATA.SIMLIB_FILE,    "%s", INPUTS.SIMLIB_FILE );


  sprintf(SNDATA.SURVEY_NAME,   "%s", SIMLIB_GLOBAL_HEADER.SURVEY_NAME );
  sprintf(SNDATA.SUBSURVEY_NAME,"%s", SIMLIB_HEADER.SUBSURVEY_NAME );  

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


  // load BYOSED info  (Dec 2018)
  if ( INDEX_GENMODEL == MODEL_BYOSED ) {
    NPAR = Event_BYOSED.NPAR ;
    SNDATA.NPAR_BYOSED = NPAR ; 
    for ( ipar=0; ipar < NPAR; ipar++ ) {
      tmpName = Event_BYOSED.PARNAME[ipar];
      sprintf(SNDATA.BYOSED_PARNAME[ipar], "%s", tmpName ) ;
      sprintf(SNDATA.BYOSED_KEYWORD[ipar], "BYOSED_PARAM(%s)", tmpName);
      SNDATA.BYOSED_PARVAL[ipar]  = Event_BYOSED.PARVAL[ipar];
    }
  } // model = MODEL_BYOSED


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

  hostgal_to_SNDATA(FLAG,0);  // for header init only.


  // #####################################

  if ( FLAG == 1 ) { return ; }

  // #####################################

  VERSION_INFO.N_SNLC       = NGENLC_WRITE ;
  VERSION_INFO.GENFRAME_SIM = GENFRAME_OPT ;


  // do smearing up front to see what's going on.
 
  SNDATA.WRFLAG_BLINDTEST = WRFLAG_BLINDTEST;

  if ( GENLC.NEPOCH >= MXEPOCH ) {
    printf("\n PRE-ABORT DUMP: \n");
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
    { SNDATA.SIM_DM15 = GENLC.DM15 ; }

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
  SNDATA.SIM_MAGSMEAR_COH     = GENLC.MAGSMEAR_COH ;
  SNDATA.SIM_RISETIME_SHIFT   = GENLC.RISETIME_SHIFT ;
  SNDATA.SIM_FALLTIME_SHIFT   = GENLC.FALLTIME_SHIFT ;
  SNDATA.SIM_TRESTMIN         = GENLC.TRESTMIN ;
  SNDATA.SIM_TRESTMAX         = GENLC.TRESTMAX ;

  // set GALID here in case HOSTLIB_USE=0 in hostgal_to_SNDATA
  if ( SNHOSTGAL.GALID>0 ) 
    { SNDATA.HOSTGAL_NMATCH[0]=1 ; SNDATA.HOSTGAL_NMATCH[1]=1 ; }
  SNDATA.HOSTGAL_OBJID[0]   = SNHOSTGAL.GALID ;

  // st HOSTLIB variables
  ifilt_obs=0 ;  hostgal_to_SNDATA(FLAG,ifilt_obs);

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];  

    MCOR_TRUE_MW  = GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs]; 
    MCOR_MAP_MW   = GENLC.MAGCOR_MWEBV_MAP[ifilt_obs];

    SNDATA.SIM_PEAKMAG[ifilt_obs]     
      = (float)( GENLC.peakmag8_obs[ifilt_obs] - MCOR_TRUE_MW) ;

    SNDATA.SIM_TEMPLATEMAG[ifilt_obs]       // LCLIB source mag in template
      = (float)GENLC.genmag8_obs_template[ifilt_obs] - MCOR_MAP_MW ;

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

    SNDATA.USE_EPOCH[epoch] = GENLC.USE_EPOCH[epoch] ;

    if ( GENLC.ISPEAK[epoch]      ) { continue ; }
    if ( GENLC.ISTEMPLATE[epoch]  ) { continue ; }

    ifilt_obs    = GENLC.IFILT_OBS[epoch];
      
    SNDATA.FILTINDX[epoch]  = ifilt_obs ;
    sprintf(SNDATA.FILTCHAR[epoch], "%c",  FILTERSTRING[ifilt_obs] );

    MCOR_TRUE_MW  = GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs]; 
    MCOR_MAP_MW   = GENLC.MAGCOR_MWEBV_MAP[ifilt_obs];

    SNDATA.SIMEPOCH_TREST[epoch]  = GENLC.epoch8_rest[epoch] ;
    SNDATA.SIMEPOCH_TOBS[epoch]   = GENLC.epoch8_obs[epoch] ;
    SNDATA.SIMEPOCH_MAG[epoch]    = GENLC.genmag8_obs[epoch] - MCOR_TRUE_MW ;
    SNDATA.SIMEPOCH_MODELMAGERR[epoch] = GENLC.generr8_rest[epoch] ;
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

    if ( GENLC.IEPOCH_SNRMAX == epoch ) 
      { SNDATA.PHOTFLAG[epoch] += SIMLIB_GLOBAL_HEADER.PHOTFLAG_SNRMAX ; }

    if ( GENLC.IEPOCH_NEARPEAK == epoch ) 
      { SNDATA.PHOTFLAG[epoch] += SIMLIB_GLOBAL_HEADER.PHOTFLAG_NEARPEAK ; }

    // load photProb
    SNDATA.PHOTPROB[epoch] = SEARCHEFF_DATA.PHOTPROB[epoch-1] ;

    // - - - - - - - - -  -
    double diff = SNDATA.MJD[epoch] - SEARCHEFF_DATA.MJD[epoch-1] ;
    if ( diff != 0.0 ) {
      printf("\n PRE-ABORT DUMP: \n");
      printf("   SNDATA.MJD[%d]=%.4f   SEARCHEFF_DATA.MJD[%d]=%.4f\n",
	     epoch, SNDATA.MJD[epoch], 
	     epoch-1, SEARCHEFF_DATA.MJD[epoch-1] );
      sprintf(c1err,"Index problem with SEARCHEFF_DATA struct.");
      sprintf(c2err,"CID=%d  epoch=%d of %d,  MJD-diff=%f", 
	      GENLC.CID, epoch, GENLC.NEPOCH, diff);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
      
    // ------------------------
    // store conditions
    SNDATA.PSF_SIG1[epoch]  = SIMLIB_OBS_GEN.PSFSIG1[epoch] ;
    SNDATA.PSF_SIG2[epoch]  = SIMLIB_OBS_GEN.PSFSIG2[epoch] ;
    SNDATA.PSF_RATIO[epoch] = SIMLIB_OBS_GEN.PSFRATIO[epoch] ;

    // Feb 23, 2012: store search-run zpt instead of template run  zpt
    SNDATA.ZEROPT[epoch]      = SIMLIB_OBS_GEN.ZPTADU[epoch] ;
    SNDATA.ZEROPT_ERR[epoch]  = SIMLIB_OBS_GEN.ZPTSIG[epoch] ;
    SNDATA.ZEROPT_SIG[epoch]  = SIMLIB_OBS_GEN.ZPTSIG[epoch] ;

    // mar 18 2018: store SNR at fixed mag to monitor data quality 
    SNDATA.SIMEPOCH_SNRMON[epoch] = GENLC.SNR_MON[epoch];

    // ----------------

    SNDATA.NPE_ABOVE_SAT[epoch]     = GENLC.npe_above_sat[epoch];
    SNDATA.FLUX[epoch]              = GENLC.flux[epoch];
    SNDATA.FLUX_ERRTOT[epoch]       = GENLC.flux_errstat[epoch]; 
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
    
    istat = fluxcal_SNDATA ( epoch, "log10" ) ; // --> fill SNDATA.FLUXCAL
    
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
  double MCOR_TRUE_MW  = GENLC.MAGCOR_MWEBV_TRUE[ifilt_obs] ;
  
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
  //   SNDATA.HOSTGAL_SB_FLUX[ifilt_obs]  => data-like quantity
  //   SNDATA.HOSTGAL_USEMASK
  //   SNDATA.HOSTGAL_MAG[ifilt_obs]   ==> data-like quanity
  // where SB = surface brightness.
  //
  // Dec 17, 2012: fill HOSTGAL_NFILT_MAGOBS and others with ifilt_obs=0
  // Feb 12, 2014: fill SNDATA.SIM_HOSTGAL_xxx, and add IFLAG arg
  // Jun 02, 2018: load zphot info for LCLIB

  int    NPAR, ipar, OVP, ifilt ;
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



  if ( ifilt_obs == 0 ) {
    if ( SNHOSTGAL.GALID>0 ) 
      { SNDATA.HOSTGAL_NMATCH[0]=1 ; SNDATA.HOSTGAL_NMATCH[1]=1 ; }
    SNDATA.HOSTGAL_OBJID[0]          = SNHOSTGAL.GALID ;
    SNDATA.HOSTGAL_PHOTOZ[0]         = SNHOSTGAL.ZPHOT ;
    SNDATA.HOSTGAL_PHOTOZ_ERR[0]     = SNHOSTGAL.ZPHOT_ERR ;

    SNDATA.HOSTGAL_SPECZ[0]          = SNHOSTGAL.ZSPEC ;
    SNDATA.HOSTGAL_SPECZ_ERR[0]      = SNHOSTGAL.ZSPEC_ERR ;
    SNDATA.HOSTGAL_RA[0]     = SNDATA.RA + 
      (SNHOSTGAL.RA_GAL_DEG-SNHOSTGAL.RA_SN_DEG);
    SNDATA.HOSTGAL_DEC[0]    = SNDATA.DEC + 
      (SNHOSTGAL.DEC_GAL_DEG-SNHOSTGAL.DEC_SN_DEG);
    SNDATA.HOSTGAL_SNSEP[0]          = SNHOSTGAL.SNSEP ;
    SNDATA.HOSTGAL_DDLR[0]           = SNHOSTGAL.DDLR ;
    SNDATA.HOSTGAL_LOGMASS[0]        = SNHOSTGAL.LOGMASS ;
    SNDATA.HOSTGAL_LOGMASS_ERR[0]    = SNHOSTGAL.LOGMASS_ERR ;

    NPAR = SNDATA.NPAR_SIM_HOSTLIB ;
    for(ipar=0; ipar < NPAR ; ipar++ ) {
      SNDATA.SIM_HOSTLIB_PARVAL[ipar] = HOSTLIB_OUTVAR_EXTRA.VALUE[ipar] ;
    }

    return ;

  } // end of ifilt_obs==0
  

  // transfer total host mag (Feb 2013)

  SNDATA.HOSTGAL_MAG[0][ifilt] = (float)SNHOSTGAL.GALMAG[ifilt_obs][0] ;
  SNDATA.HOSTGAL_USEMASK |= 1 ; // flag to write host mag

  SNDATA.HOSTGAL_SB_FLUX[ifilt] = (float)SNHOSTGAL.SB_FLUX[ifilt_obs];
  SNDATA.HOSTGAL_SB_MAG[ifilt]  = (float)SNHOSTGAL.SB_MAG[ifilt_obs];

  OVP = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_PHOT ;
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

    Nov 30, 2006: compute extinction (x8) with get_snxt8 function.

    Dec 4, 2006: use fortran NEAREST_IFILT_REST and pass option

    Jun 7, 2007: take wgted average of two nearest filters
                 Weights wgt1 and wgt2 are inversely proportional
                 to lambda_obs/(1+z) - lambda_rest

   Mar 5, 2008: if ifilt_rest2 is redder than ifilt_obs (for z ~ 0.02)
                then do NOT get kcor8_2 to avoid fatal error. 
                Needed to simulate LOWZ.

   Jun 11, 2008: remove MWXT calculation from here and move 
                 it to genmag_MWXT().

   May 23, 2009: add istat output arg to get_avwarp8(),
                 and increment AVwarp-overflows vs. filter

   Sep 8, 2009: use new kcorfun8 function to do both AVwarp
                and K-corrections. Eventually will do wgted
                avg K-cor with both colors.  Get rid of lots
                of old unecessary code.

   Jan 22, 2010: 
    if mag8[1]  or mag8[2] > is crazy, set mag_obs = 666 as a flag
    to set flux and fluxerr = -9 in gen_observerSmear().

    Remove obsolete 'istat' dependencies.

   Jun 11, 2010: if Trest < 0 then set undefined mags to 99 instead
                 of 666 so that pre-explosion epochs are written out.

  *****/

  double 
    AVwarp8[4], AV8, RV8, Z8, T8,  x8[10]
    ,mag8[4], lamdif8[4], mag8_obs, kcor8
    ;

  int ifilt_obs, ifilt_rest1, ifilt_rest2, ifilt_rest3, epoch, NZ ;
  int ifilt_rest_tmp[4] ;   

  char fnam[] = "genmag_boost" ;

  // ------------------ BEGIN ------------

  // do NOT apply extinction if AV=0

  if ( GENLC.AV  == 0.0 ) { goto KCOR ; }

#ifdef SNGRIDGEN
  // skip boost for SNOOPY model with just 1 logz-bin
  NZ = GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_LOGZ];
  if ( NZ == 1 &&  INDEX_GENMODEL == MODEL_SNOOPY ) { return ; }
#endif


  // apply extinction in rest frame mags/filters
  AV8 = (double)GENLC.AV ;
  RV8 = (double)GENLC.RV ;
  Z8  = 0.0;

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {  


    T8 = GENLC.epoch8_rest[epoch]; 
    if ( INPUTS.KCORFLAG_STRETCH == 1 ) T8 /= GENLC.STRETCH ;

    ifilt_obs   = GENLC.IFILT_OBS[epoch];

      if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

      ifilt_rest1 = GENLC.IFILTMAP_REST1[ifilt_obs];
      ifilt_rest2 = GENLC.IFILTMAP_REST2[ifilt_obs];
      ifilt_rest3 = GENLC.IFILTMAP_REST3[ifilt_obs];     
       
      // start with nearest filter.
      x8[1] = get_snxt8__( &OPT_SNXT, &ifilt_rest1, &T8, &AV8, &RV8 );
      GENLC.genmag8_rest[epoch] += x8[1] ;

      // now do 2nd nearest filter
      x8[2] = get_snxt8__( &OPT_SNXT, &ifilt_rest2, &T8, &AV8, &RV8 );
      GENLC.genmag8_rest2[epoch] += x8[2] ;

      // 3rd nearest filter
      x8[3] = get_snxt8__( &OPT_SNXT, &ifilt_rest3, &T8, &AV8, &RV8 );
      GENLC.genmag8_rest3[epoch] += x8[3] ;


  } // end of epoch loop


 KCOR:

  Z8     = GENLC.REDSHIFT_HELIO ;

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {  

    T8      = GENLC.epoch8_rest[epoch]; 
    if ( INPUTS.KCORFLAG_STRETCH == 1 ) { T8 /= GENLC.STRETCH ; }

    ifilt_obs   = GENLC.IFILT_OBS[epoch];
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

    ifilt_rest1 = GENLC.IFILTMAP_REST1[ifilt_obs] ;
    ifilt_rest2 = GENLC.IFILTMAP_REST2[ifilt_obs] ;
    ifilt_rest3 = GENLC.IFILTMAP_REST3[ifilt_obs] ;

    lamdif8[1]  = GENLC.LAMDIF_REST1[ifilt_obs];
    lamdif8[2]  = GENLC.LAMDIF_REST2[ifilt_obs];
    lamdif8[3]  = GENLC.LAMDIF_REST3[ifilt_obs];

    ifilt_rest_tmp[1] = ifilt_rest1 ;
    ifilt_rest_tmp[2] = ifilt_rest2 ;
    ifilt_rest_tmp[3] = ifilt_rest3 ;

    mag8[1]  = GENLC.genmag8_rest[epoch]  ;
    mag8[2]  = GENLC.genmag8_rest2[epoch] ;    
    mag8[3]  = GENLC.genmag8_rest3[epoch] ;

    kcor8 = kcorfun8_ ( &ifilt_obs, &ifilt_rest_tmp[1], 
			&mag8[1], &lamdif8[1],  &T8, &Z8, &AVwarp8[1] );

    if ( isnan( AVwarp8[2]) ) {
      sprintf(c1err,"AVwarp=nan for T8=%5.1f  mag8[1,2]=%6.2f,%6.2f",
	      T8, mag8[1], mag8[2] );
      sprintf(c2err,"ifilt_rest[1,2]=%d,%d (%c,%c) "
	      ,ifilt_rest1, ifilt_rest2
	      ,FILTERSTRING[ifilt_rest1]
	      ,FILTERSTRING[ifilt_rest2]
	      );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if ( mag8[1] < -10. ) {
      NAVWARP_OVERFLOW[0]++ ;             // grand total
      NAVWARP_OVERFLOW[ifilt_rest1]++ ;   // filter-dependent sum
    }

    // add up all the contributions

    mag8_obs = mag8[1] + kcor8 + GENLC.DLMU ; 

    // Jan 22, 2010: 
    // if the spectral warping is crazy, set mag_obs to crazy value
    if ( T8 > 0 ) {
      if ( fabs(mag8[1]) > 90.0 ) { mag8_obs = MAG_UNDEFINED ; }
      if ( fabs(mag8[2]) > 90.0 ) { mag8_obs = MAG_UNDEFINED ; }
    }
    else {
      if ( fabs(mag8[1]) > 90.0 ) mag8_obs = 99;
      if ( fabs(mag8[2]) > 90.0 ) mag8_obs = 99;
    }

    if ( fabs(T8) < -90.01 ) {
      printf(" BOOST: %c(%c) -> %c : Trest=%6.1f (MJD=%7.1f) iep=%d \n"
	     ,FILTERSTRING[ifilt_rest1]
	     ,FILTERSTRING[ifilt_rest2]
	     ,FILTERSTRING[ifilt_obs]
	     ,T8, GENLC.MJD[epoch], epoch ) ;

      printf("\t M%c = %6.2f(M%c) + %6.3f(kcor) + %7.3f(mu) = %7.3f \n"
	     ,FILTERSTRING[ifilt_obs]
	     ,mag8[1]
	     ,FILTERSTRING[ifilt_rest1]
	     ,kcor8, GENLC.DLMU,  mag8_obs );

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

    GENLC.kcorval8[epoch]        = kcor8 ;
    GENLC.warpcolval8[epoch]     = mag8[1] - mag8[2] ;
    GENLC.AVwarp8[epoch]         = AVwarp8[2] ;
    GENLC.ifilt_AVwarp[epoch][0] = ifilt_rest1 ; // closest filter
    GENLC.ifilt_AVwarp[epoch][1] = ifilt_rest2 ; // 2nd closest

      // store true obs mag in GENLC structure

    GENLC.genmag8_obs[epoch] = mag8_obs;

    // load model mag-err for observer frame using the model-error
    // from the nearest rest-frame filter (ifilt_rest1).

    GENLC.generr8_obs[epoch] = GENLC.generr8_rest[epoch] ;


  } // end of epoch loop

 
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
    if ( GENLC.genmag8_obs[epoch]  > 30.0 ) { continue ; }

    ifilt_obs  = GENLC.IFILT_OBS[epoch];
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

    Trest   = GENLC.epoch8_rest[epoch]; 
    if ( INPUTS.KCORFLAG_STRETCH == 1 ) { Trest /= GENLC.STRETCH ; }

    AVwarp  = GENLC.AVwarp8[epoch] ;
    MWXT  = get_mwxt8__(&ifilt_obs, &Trest, &z, &AVwarp, &mwebv, &RV, &OPT);
      
    // increment observer-frame magnitude.
    GENLC.genmag8_obs[epoch] += MWXT ;      


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

	GENFILT.Trest8[ifilt_obs][NEP]   = GENLC.epoch8_rest[epoch];
	GENFILT.Tobs8[ifilt_obs][NEP]    = GENLC.epoch8_obs[epoch];

	GENFILT.genmag8_obs[ifilt_obs][NEP]    = GENLC.genmag8_obs[epoch];
	GENFILT.genmag8_rest[ifilt_obs][NEP]   = GENLC.genmag8_rest[epoch];
	GENFILT.genmag8_rest2[ifilt_obs][NEP]  = GENLC.genmag8_rest2[epoch];
	GENFILT.genmag8_rest3[ifilt_obs][NEP]  = GENLC.genmag8_rest3[epoch];

	GENFILT.genmag8_smear[ifilt_obs][NEP]  = GENLC.magsmear8[epoch];

	GENFILT.generr8_obs[ifilt_obs][NEP]    = GENLC.generr8_obs[epoch];
	GENFILT.generr8_rest[ifilt_obs][NEP]   = GENLC.generr8_rest[epoch];
	GENFILT.generr8_rest2[ifilt_obs][NEP]  = GENLC.generr8_rest2[epoch];

	//printf(" xxx load Trest(%c=%d) = %f at  epoch=%d,  FILTEPOCH=%d \n", 
	// FILTERSTRING[ifilt_obs], ifilt_obs, GENLC.epoch8_rest[epoch], epoch,NEP);
      }
      else if ( opt == -1 ) {

	// first make sure that Trest matches up

	Trest1	= GENFILT.Trest8[ifilt_obs][NEP] ;
	Trest2  = GENLC.epoch8_rest[epoch];
	Tdif    = fabs(Trest1-Trest2) ;
	if ( Tdif > 0.001 ) {
	  sprintf(c1err,"GENFILT[ifilt_obs=%d][NEP=%d] = %8.3f", 
		  ifilt_obs, NEP, Trest1);
	  sprintf(c2err,"GENLC.epoch8_rest[ep=%d]=%8.3f  (Tdif=%8.3f) \n", 
		  epoch, Trest2, Tdif );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	}

	GENLC.genmag8_obs[epoch]   = GENFILT.genmag8_obs[ifilt_obs][NEP] ;
	GENLC.genmag8_rest[epoch]  = GENFILT.genmag8_rest[ifilt_obs][NEP] ; 
	GENLC.genmag8_rest2[epoch] = GENFILT.genmag8_rest2[ifilt_obs][NEP] ; 
	GENLC.genmag8_rest3[epoch] = GENFILT.genmag8_rest3[ifilt_obs][NEP] ; 

	GENLC.magsmear8[epoch]   = GENFILT.genmag8_smear[ifilt_obs][NEP] ;

	GENLC.generr8_obs[epoch]   = GENFILT.generr8_obs[ifilt_obs][NEP] ;
	GENLC.generr8_rest[epoch]  = GENFILT.generr8_rest[ifilt_obs][NEP] ;
	GENLC.generr8_rest2[epoch] = GENFILT.generr8_rest2[ifilt_obs][NEP] ;
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

  ************/

  char *GENMODEL        = INPUTS.GENMODEL;
  char *GENMODEL_EXTRAP = INPUTS.GENMODEL_EXTRAP_LATETIME ;

  char  covFile[] = ""  ;
  char *ARGLIST ;
  int istat, OPTMASK, NZ, ifilt, ifilt_obs, ifilt_rest ;
  double dummy[10];
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
    // xxx mark delete  GENLC.SIMTYPE  = (int)INPUTS.FIXMAG[0] ;
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
#ifdef SNGRIDGEN
    NZ = GRIDGEN_INPUTS.NBIN[IPAR_GRIDGEN_LOGZ];    
    if ( NZ == 1 ) { OPTMASK = 1; }
#endif
    init_genmag_snoopy(GENMODEL, OPTMASK, dummy, GENLC.FILTLIST_REST );

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
    OPTMASK = 0; 
    if ( INPUTS.LEGACY_colorXTMW_SALT2 ) { OPTMASK += 128 ; }

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
    if( INPUTS.USE_BINARY_SIMSED > 0 ) { OPTMASK += 1; }

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

  else if ( INDEX_GENMODEL == MODEL_BYOSED ) {

    OPTMASK  = INPUTS.GENMODEL_MSKOPT;
    ARGLIST  = INPUTS.GENMODEL_ARGLIST;
    int NPAR; char NAMES_HOSTPAR[200];  double VAL_HOSTPAR=0.0 ;
    NPAR = fetch_HOSTPAR_GENMODEL(1, NAMES_HOSTPAR, &VAL_HOSTPAR);

    // init generic part of any SEDMODEL (filter & primary ref)
    init_genSEDMODEL();
    init_genmag_BYOSED( INPUTS.MODELPATH, OPTMASK, ARGLIST, NAMES_HOSTPAR );
  }

  else if ( INDEX_GENMODEL == MODEL_NON1ASED ) {

    INPUTS.NON1ASED.NGENTOT  = INPUTS.NGEN ;
    INPUTS.NON1ASED.CIDOFF   = INPUTS.CIDOFF ;
    GENLC.NON1ASED.IFLAG_GEN = GENLC.IFLAG_GENSOURCE ;
    
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
    //  REDSHIFT parameter. Note that redshift range is read from 
    //  the LCLIB header, not from sim-input file.
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

    if ( INDEX_GENMODEL == MODEL_MLCS2k2 )   init_covar_mlcs2k2(); 

  }

  if ( LGEN_SNIA ) 
    { GENLC.SIMTYPE  = INPUTS.SNTYPE_Ia_SPEC ; }

  if ( INPUTS.GENTYPE_SPEC > 0 ) 
    { GENLC.SIMTYPE  = INPUTS.GENTYPE_SPEC ; }
 
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
  char  filtName[40], cfilt[2]    ;

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

    get_filttrans__(&MASKFRAME[ifilt],     // (I) obs/rest-frame mask
		    &ifilt_obs,            // (I) absolute filter index
		    filtName,              // (O) full name of filter
		    &magprim,              // (O) mag of primary ref 
		    &NLAM,                 // (O) Number of lambda bins
		    genSEDMODEL.lam,       // (O) lambda array 
		    genSEDMODEL.TransSN,   // (O) filter trans 
		    genSEDMODEL.TransREF,  // (O) idem
		    LEN ); 

    if ( NLAM > MXLAMSIM ) {
      sprintf(cfilt,  "%c", FILTERSTRING[ifilt_obs] );
      sprintf(c1err,"NLAM(%s) = %d exceeds bound of MXLAMSIM=%d",
	      cfilt, NLAM, MXLAMSIM );
      sprintf(c2err,"Check filter trans files.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    lamshift = 0.0 ;
    init_filter_SEDMODEL(ifilt_obs, filtName, magprim, NLAM, 
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
  if ( index_model == MODEL_LCLIB     ) return ;
  if ( index_model == MODEL_BYOSED    ) return ;

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
void init_kcor(char *kcorFile) {

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

  *********/

  int ISMODEL_FIXMAG = ( INDEX_GENMODEL == MODEL_FIXMAG );
  int ierrstat, ifilt, ifilt_obs, OPT ;
  float tmpoff_kcor[MXFILTINDX] ;
  float *ptr ;
  char   copt[40], xtDir[MXPATHLEN], cfilt[4], *NAME;
  char fnam[] = "init_kcor" ;

  // -------------- BEGIN --------------

  // set fortran arrays needed for K-correction lookup
  // The GENLC variables are set in SIMLIB_open
 
  set_survey__( GENLC.SURVEY_NAME  
		,&GENLC.NFILTDEF_OBS
		,GENLC.IFILTMAP_OBS
		,INPUTS.TMPOFF_LAMSHIFT
		,strlen(GENLC.SURVEY_NAME)   
		);

  // read K-cor and mag tables (vs. Z, epoch, AV)
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
    ptr = &INPUTS.GENMAG_OFF_ZP[ifilt_obs] ;

    if ( !ISMODEL_FIXMAG ) {
      *ptr += tmpoff_kcor[ifilt]; // add AB off to user offset
    }

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

  if ( ISMODEL_FIXMAG ) { return ; } // 4.24.2019

  // init MW extinction (moved from end of init_genmodel on Aug 30 2010)
  if ( GENFRAME_OPT == GENFRAME_REST )  {

    // Sep 2013: fetch params used to compute MW extinct in kcor files.
    //           Needed to re-compute MWXT if RV is changed.
    get_kcor_mwpar__(&GENLC.kcor_RVMW, &GENLC.kcor_OPT_MWCOLORLAW );

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

    if ( OPT_SNXT  == OPT_SNXT_CCM89 ) {
      init_xthost__(&OPT_SNXT);
    }
    else if ( OPT_SNXT  == OPT_SNXT_SJPAR ) {
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
    printf("   Found %d synthetic spectrograph filters (%s) \n",
	   GENLC.NFILTDEF_SPECTROGRAPH, GENLC.FILTERLIST_SPECTROGRAPH );
    fflush(stdout);
  }

  return ;

} // end of init_kcor
 

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

  int  LFIND_SPEC=0, NEP_PEAKONLY=0, iep,ifilt ;
  int  MEMI  = (GENLC.NEPOCH+1) * sizeof(int);
  int  MEMD  = (GENLC.NEPOCH+1) * sizeof(double);
  int  DOSPEC = (INPUTS.APPLY_SEARCHEFF_OPT & 2) ;
  double EFF ;

  struct {
    int NEPOCH, *IFILT_OBS, *ISPEAK;
    double *MJD, *TOBS, *TREST ;
  } GENLC_ORIG ;


  char fnam[] = "gen_TRIGGER_PEAKMAG_SPEC" ;

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
    GENLC_ORIG.ISPEAK[iep]    = GENLC.ISPEAK[iep] ;
    GENLC_ORIG.IFILT_OBS[iep] = GENLC.IFILT_OBS[iep] ;
    GENLC_ORIG.MJD[iep]       = GENLC.MJD[iep];
    GENLC_ORIG.TOBS[iep]      = GENLC.epoch8_obs[iep];  
    GENLC_ORIG.TREST[iep]     = GENLC.epoch8_rest[iep] ;

    if ( GENLC_ORIG.ISPEAK[iep] == 0 ) { continue ; }
    NEP_PEAKONLY++ ;
    GENLC.ISPEAK[NEP_PEAKONLY]       = GENLC_ORIG.ISPEAK[iep] ;
    GENLC.IFILT_OBS[NEP_PEAKONLY]    = GENLC_ORIG.IFILT_OBS[iep] ;
    GENLC.MJD[NEP_PEAKONLY]          = GENLC_ORIG.MJD[iep] ;
    GENLC.epoch8_obs[NEP_PEAKONLY]   = GENLC_ORIG.TOBS[iep] ;
    GENLC.epoch8_rest[NEP_PEAKONLY]  = GENLC_ORIG.TREST[iep] ;
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
      GENLC.ISPEAK[iep]      = GENLC_ORIG.ISPEAK[iep];
      GENLC.IFILT_OBS[iep]   = GENLC_ORIG.IFILT_OBS[iep] ;
      GENLC.MJD[iep]         = GENLC_ORIG.MJD[iep];
      GENLC.epoch8_obs[iep]  = GENLC_ORIG.TOBS[iep]  ;
      GENLC.epoch8_rest[iep] = GENLC_ORIG.TREST[iep] ;
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
  char fnam[] = "gen_TRIGGER_zHOST" ;

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
  int NEPCUT = INPUTS.CUTWIN_NEPOCH[0];
  char fnam[] = "GENMAG_DRIVER" ;

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

  // Created July 2016
  // [for spectra refactor, move code from main to here]
  // Driver routine to generate observed fluxes and uncertainties
  // from true/generated mags.
  //
  // May 2018: set GENLC.NOBS_SATURATE
  // Oct 30 2018: fix bug setting L_SATURATE

  int epoch, istat, ifilt_obs, ISPEAK, L_ERRPOS, L_UNDEFINED, L_SATURATE ;
  int VBOSE_SMEAR = 0;
  double genmag, obsmag;
  char fnam[] = "GENFLUX_DRIVER" ;

  // -------------- BEGIN ---------------

  genran_obsNoise();   // randoms for instrumental noise

  for ( epoch = 1; epoch <= GENLC.NEPOCH; epoch++ ) {

    ifilt_obs = GENLC.IFILT_OBS[epoch] ;
    if ( GENLC.DOFILT[ifilt_obs] == 0 ) { continue ; }

    if (   GENLC.ISPEAK[epoch]      ) { continue; }
    if (   GENLC.ISTEMPLATE[epoch]  ) { continue; }

    // convert 'genmag' into Possion-smeared mag and flux
    istat =  gen_smearFlux ( epoch, VBOSE_SMEAR );
    istat =  gen_smearMag  ( epoch, VBOSE_SMEAR );

    obsmag      = GENLC.mag[epoch];
    genmag      = GENLC.genmag8_obs[epoch] ;
    L_ERRPOS    = (GENLC.flux_errstat[epoch] > 0 );
    L_UNDEFINED = (genmag == MAG_UNDEFINED) ; // model undefined
    L_SATURATE  = (obsmag == MAG_SATURATE ) ;    

    if ( L_UNDEFINED ) 
	{ GENLC.NOBS_UNDEFINED++ ; } // model is undefined 

    if ( L_ERRPOS ) {
      if ( !L_UNDEFINED ) { 
	GENLC.USE_EPOCH[epoch] = 1 ; 
	GENLC.NOBS++ ;
	GENLC.NOBS_FILTER[ifilt_obs]++ ;
      }
  
      if ( L_SATURATE ) { 
	GENLC.NOBS_SATURATE[1]++ ; 
	GENLC.NOBS_SATURATE_FILTER[1][ifilt_obs]++ ; 
      }
      else {
	// number of unsaturated epochs
	GENLC.NOBS_SATURATE[0]++; 
	GENLC.NOBS_SATURATE_FILTER[0][ifilt_obs]++ ; 
      }
    } // end L_ERRPOS

  }  // end of epoch loop
  

  return ;

} // end GENFLUX_DRIVER


// **************************************
void compute_lightCurveWidths(void) {

  // call gen_lightCurveWidth for each band,
  // and store GENLC.WIDTH[ifilt_obs].
  // These quantities are for SIMGEN_DUMP file
  // (varNames WIDTH_[band] ).

  int ifilt, ifilt_obs, N, ep, NEP[MXFILTINDX], ERRFLAG ;
  int MEM = (GENLC.NEPOCH+1) * sizeof(double) ;
  int OPTMASK = 1 ;
  double Width, *TLIST[MXFILTINDX], *MAGLIST[MXFILTINDX] ;
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
    TLIST[ifilt_obs][N]   = GENLC.epoch8_obs[ep] ;
    MAGLIST[ifilt_obs][N] = GENLC.genmag8_obs[ep] ;
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
   Output is GENLC.genmag8_obs[ifilt][epoch] for obs-frame model,
   or GENLC.genmag8_rest[ifilt][epoch] for rest frame model.
   Note that for rest-frame models, genmag_boost transforms
   rest-frame mags to observer frame.

  ***********/

  int
    istat, ifilt_tmp, ifilt_rest, iep
    ,NEPFILT, NGRID, NEPFILT_SAVE
    ,index, isp, OPTMASK, NPAR
    ;

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

    ptr_epoch    = &GENFILT.Trest8[ifilt_obs][1] ;

    if ( inear == 1 ) {
      ptr_genmag   = &GENFILT.genmag8_rest[ifilt_obs][1];
      ptr_generr   = &GENFILT.generr8_rest[ifilt_obs][1];
      ifilt_rest   = GENLC.IFILTMAP_REST1[ifilt_obs];
    }
    else if ( inear == 2 ) {
      ptr_genmag   = &GENFILT.genmag8_rest2[ifilt_obs][1];
      ptr_generr   = &GENFILT.generr8_rest2[ifilt_obs][1];
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

      ptr_genmag   = &GENFILT.genmag8_rest3[ifilt_obs][1];
      ptr_generr   = &GENFILT.generr8_rest3[ifilt_obs][1];
      ifilt_rest   = GENLC.IFILTMAP_REST3[ifilt_obs];
    }

    sprintf(cfilt_rest, "%c", FILTERSTRING[ifilt_rest] );

  }
  else {
    // observer-frame model
    ptr_genmag   = &GENFILT.genmag8_obs[ifilt_obs][1] ;
    ptr_epoch    = &GENFILT.Tobs8[ifilt_obs][1] ;
    ptr_generr   = &GENFILT.generr8_obs[ifilt_obs][1];    
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

    /* xxx mark delete April 24 2019 xxxxxxxx
    if ( strcmp(INPUTS.MODELNAME,"FIXMAG") == 0 ) {
      sprintf(GENLC.SNTYPE_NAME, "FIXMAG" ); 
      GENLC.TEMPLATE_INDEX = MODEL_FIXMAG ;  //anything but zero
    }
    else {
      sprintf(GENLC.SNTYPE_NAME, "RANMAG" ); 
      GENLC.TEMPLATE_INDEX = (int)INPUTS.FIXMAG[0] ;
    }
    xxxxxxxxxxx end mark xxxxxxxxxxxxx */

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
    double S2x0, S2x1, S2c, S2mB, tmp ;
    S2mB = GENLC.SALT2mB + GENLC.COVMAT_SCATTER[0] ;
    S2x1 = GENLC.SALT2x1 + GENLC.COVMAT_SCATTER[1] ;
    S2c  = GENLC.SALT2c  + GENLC.COVMAT_SCATTER[2] ;
    tmp  = -0.4 * GENLC.COVMAT_SCATTER[0] ;
    S2x0 = GENLC.SALT2x0 * pow(10.0,tmp);


    genmag_SALT2 (
		  OPTMASK         // (I) bit-mask options
		  ,ifilt_obs      // (I) obs filter index 
		  ,S2x0           // (I) x0 term
		  ,S2x1,S2x1      // (I) x1 and x1 used for error calc.
		  ,S2c            // (I) SN color: E(B-V)
		  ,mwebv          // (I) Galactic E(B-V)
		  ,RV, AV         // (I) host extinc params (July 2016)
		  ,z,z            // (I) redshift, and z used for error
		  ,NEPFILT        // (I) number of epochs
		  ,ptr_epoch         // (I) obs-frame time (days)
		  ,ptr_genmag        // (O) mag vs. Tobs
		  ,ptr_generr        // (O) mag-errs
		  ) ;    
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

    if ( INPUTS.OPTMASK_SIMSED == OPTMASK_SIMSED_GRIDONLY ) 
      { GENLC.CID = GENLC.TEMPLATE_INDEX ; }

  }

  else if ( INDEX_GENMODEL  == MODEL_BYOSED ) {

    int NHOSTPAR; char *NAMES_HOSTPAR; double VAL_HOSTPAR[MXHOSTPAR_BYOSED];
    NHOSTPAR = fetch_HOSTPAR_GENMODEL(2, NAMES_HOSTPAR, VAL_HOSTPAR);

    genmag_BYOSED(
		  GENLC.CID
		  ,z, GENLC.DLMU       // (I) helio-z and distance modulus
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

    ifilt_tmp = GENLC.IFILTINVMAP_REST[ifilt_rest] + 1 ;
    genmag_snoopy (
		   ifilt_tmp             // (I) SNoopY internal filter index
		   ,dm15                 // (I) shape parameter 
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
		   ,&GENLC.MAGSMEAR_COH // (O) store random magSmear
		      );

    GENLC.TEMPLATE_INDEX = fetchInfo_NON1AGRID("NON1A_INDEX");
    GENLC.SNTYPE         = fetchInfo_NON1AGRID("NON1A_ITYPE_USER");
  }
  else if ( INDEX_GENMODEL  == MODEL_LCLIB ) {  // July 2017

    set_TobsRange_LCLIB(GENLC.epoch8_obs_range);

    sprintf(GENLC.SNTYPE_NAME, "%s", LCLIB_INFO.NAME_MODEL ); 
    double *ptr_template = &GENLC.genmag8_obs_template[ifilt_obs];
    double LAMAVG        = INPUTS.LAMAVG_OBS[ifilt_obs];
    double ZTRUE, TobsPeak ;

    genmag_LCLIB(GENLC.CID, ifilt_obs, GENLC.RA, GENLC.DEC, 
		 mwebv, LAMAVG,
		 NEPFILT, ptr_epoch, 
		 ptr_genmag, ptr_template,    // return args
		 &TobsPeak                    // Tobs with peak brightness
		 );

    // for non-recurring events, set PKMJD like any other transient
    if ( LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_NONRECUR ) {
      GENLC.PEAKMJD  = INPUTS.GENRANGE_PEAKMJD[0] + TobsPeak ;
      // xxx mark delete  GENLC.PEAKMJD_SMEAR         = GENLC.PEAKMJD ;
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
    GENLC.peakmag8_rest[ifilt_rest] = ptr_genmag[iep] ;
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
  // ifilt_obs   : observer filter index
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
    ,magSmear, magSmear_model
    ,smearsig, smearsig_fix, smearsig_model
    ,Tep, Tpeak, Trest, lamrest, Z1
    ;

  float lamavg4, lamrms4, lammin4, lammax4 ;

  char cfilt[2];
  char fnam[] = "genmodelSmear" ;

  // -------------- BEGIN ------------

  if ( INPUTS.DO_MODELSMEAR  == 0 ) 
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
    GENFILT.genmag8_smear[ifilt_obs][iep] = 0.0 ; 
  }


  // ----------------------------------------------
  // start with coherent smearing (GENMAG_SMEAR)
  // that is common to all passbands.

  if ( INPUTS.GENMAG_SMEAR[0] > 0.0001 ) {

    ran_COH = GENLC.GENSMEAR_RANGauss_FILTER[0] ; // retrieve random #

    magSmear = INPUTS.GENMAG_SMEAR[0] * ran_COH ;
    GENLC.MAGSMEAR_COH = magSmear ;  // store global

    for ( iep = 1; iep <= NEPFILT; iep++ )   { 
      *(ptr_genmag+iep-1)                   += magSmear ;
      GENFILT.genmag8_smear[ifilt_obs][iep] += magSmear ;
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
  if ( INPUTS.NFILT_SMEAR > 0 ) 
    { smearsig_fix = INPUTS.GENMAG_SMEAR_FILTER[ifilt_local]; }

  // loop over epochs and apply same mag-smear 

  for ( iep=1; iep <= NEPFILT; iep++ ) {

    smearsig_model = 0.0 ;
    Tep   = *(ptr_epoch + iep - 1); // rest or obs.
    Trest = Tep / Z1 ;  // Z1=1 (rest) or 1+z (obs)
     
    if ( USE_GENMODEL_ERRSCALE  ) 
      { smearsig_model = genSmear_ERRSCALE(ptr_generr,iep, NEPFILT ); }

    // add two sources of smear in quadrature: FILTER and ERRSCALE
    smearsig = sqrt( pow(smearsig_model,2.) + pow(smearsig_fix,2.) );      
    magSmear = smearsig * ran_FILT ;  // apply filter-dependent random smear


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
      get_genSmear(Trest, ONE,  &lamrest, &magSmear_model);
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

      printf("\n xxx smearsig(model,fix,tot) = %f , %f , %f \n", 
	     smearsig_model, smearsig_fix, smearsig );
      
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
    
    ptr_genmag[iep-1]                     += magSmear ;
    GENFILT.genmag8_smear[ifilt_obs][iep] += magSmear ;

  } // ep loop

  if ( istat_genSmear() > 0 ) 
    { GENLC.MAGSMEAR_COH = MAGSMEAR_COH ; } // see sntools_genSmear.c

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
  SMEAR1 = GENFILT.genmag8_smear[ifilt1_obs][iep] ;  
  SMEAR2 = GENFILT.genmag8_smear[ifilt2_obs][iep] ;  
  
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

  char fname[28] = "INIT_COVMAT_SCATTER" ;
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
    errmsg( SEV_FATAL,0 , fname, c1err, c2err);
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
	{ errmsg( SEV_FATAL,0 , fname, c1err, c2err); }

      if ( (t1>0) && (t3>0)  ) 
	{ errmsg( SEV_FATAL,0 , fname, c1err, c2err ); }
    }
  }
  

  //If INPUTS.COVMAT_SCATTER are negative abort
  strcpy(c1err," diagonal entries of covariance matrix negative");
  for (m = 0 ; m < 3 ; m++){
    sprintf(c2err, "INPUTS.COVMAT_SCATTER[%d][%d]", m , m);
    if (INPUTS.COVMAT_SCATTER[m][m] < -epsilon ) 
      { errmsg( SEV_FATAL,0 , fname, c1err , c2err); }
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
	{  errmsg( SEV_FATAL,0 , fname, c1err, c2err);  }

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
  strcpy(c1err,"Submatrices are not positive definite for the combination ");

  for (m =0; m < 3 ; m++){
    n = m + 1; 
    
    if (3 == n ) {n = 0;}
    
    diff = (INPUTS.COVMAT_SCATTER[m][m]*
	    INPUTS.COVMAT_SCATTER[n][n]) - 
      (INPUTS.COVMAT_SCATTER[m][n]*
       INPUTS.COVMAT_SCATTER[m][n]);

    sprintf(c2err, "%d\t%d", m , n);
    
    if (diff < 0) {
      printf("PRE-ABORT DUMP \n\t %d \t %d \t %g\n", m, n , diff);
      errmsg( SEV_FATAL,0 , fname, c1err, c2err);
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
int GEN_COVMAT_SCATTER ( double *randoms, double *SCATTER_VALUES  ) 
{
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

  char ctmp[MXPATHLEN], cfilt[2], cwd[MXPATHLEN] ;
  char *cptr;
  char conoff[2][4] = { "OFF" , "ON" } ;

  int i, j, j2, ifilt_obs, ifilt_rest, ifilt, itmp, imap, iopt, ipar ;
  int NLINE, NOV, NON1A_non1a, OVP1, OVP2, OVP3 ;

  double XN, XNERR ;
  float xt, xtprod, val, ZMIN, ZMAX, shift[2];

  char  fnam[] = "readme_doc" ;

  // ------------ BEGIN readme_doc() ---------------

  i=0 ;

  print_banner ( " Fill comments for README doc-file" );

  if ( iflag_readme == 2 ) goto AFTERSIM ;

  //--- brief description

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr," USER COMMENT: %s \n\n", INPUTS.COMMENT);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr," BRIEF_DESCRIPTION: simulate %s SNe with GENMODEL = %s \n", 
	  INPUTS.GENSOURCE, INPUTS.MODELNAME );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  HOST MACHINE: %s \n", getenv("HOST") );

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
    sprintf(cptr,"  Current Dir:  %s \n", cwd );
  }

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
    if ( INPUTS.KCORFLAG_STRETCH == 1 ) strcat ( cptr, "/STRETCH" );
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
  j = INPUTS.SMEARFLAG_ZEROPT ;
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
  sprintf(cptr,"\t Peculiar Velocity Gaussian sigma: %.1f km/sec\n", 
	  INPUTS.GENSIGMA_VPEC );


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
  // xxx mark delete May 2019  if ( fabsf(xtprod-1.0) > .001 ) {
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

  // ---- FIXMAG params
  readme_doc_FIXMAG(&i);


  // --------------------------------
  // dump host extinction params

  readme_doc_hostxt(&i);


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
    // xxx mark delete   sprintf(ctmp, "%s NONE. ", ctmp);
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
  // xxx mark delete  sprintf(ctmp,"%s\n", ctmp);
  strcat(cptr,ctmp);
  // xxx mark delete  sprintf(cptr,"%s", ctmp);
  sprintf(WARNING_AVWARP_OVERFLOW,"\n  WARNING: %s", ctmp); 

  readme_doc_TAKE_SPECTRUM(&i);

  // ----- cosmology parameters

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Cosmology Parameters: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t H0 = %6.2f km/s per MPc \n", INPUTS.H0 );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Omega_{M,L} = %5.2f, %5.2f   w = %5.2f  \n",
	  INPUTS.OMEGA_MATTER, INPUTS.OMEGA_LAMBDA, INPUTS.W0_LAMBDA );



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
    for ( itmp = 1; itmp <= NLINE; itmp++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s", SEARCHEFF_DETECT[imap].README[itmp] ) ;
    }    
  }

  // print detection logic
  NLINE = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ; 
  for ( itmp = 1; itmp <= NLINE; itmp++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[itmp] ) ;
  }    
  

  // ------ Spec Search efficiency

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s", "\n  Spectroscopic Efficiency : \n" );
  for ( iopt=1; iopt <= SEARCHEFF_SPEC_INFO.NLINE_README; iopt++ ) {
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
  


  // write APPLY-mask
  OVP1 = (INPUTS.APPLY_SEARCHEFF_OPT & 1);
  OVP2 = (INPUTS.APPLY_SEARCHEFF_OPT & 2);
  OVP3 = (INPUTS.APPLY_SEARCHEFF_OPT & 4);

  if ( OVP1 && OVP2 ) 
    { sprintf(ctmp,"%s", "Apply both PIPELINE+SPEC efficiencies"); }
  else if ( OVP1 && OVP2 == 0 ) 
    { sprintf(ctmp,"%s", "Apply only PIPELINE efficiency (not SPEC-eff)"); }
  else if ( OVP1 == 0 && OVP2 ) 
    { sprintf(ctmp,"%s", "Apply only SPEC efficiency (not PIPELINE-eff)"); }
  else if ( OVP1 == 0 && OVP2 == 0 ) 
    { sprintf(ctmp,"%s", "Do NOT Apply trigger efficiency "); }
      
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  APPLY_SEARCHEFF_OPT:  %d => %s \n", 
	  INPUTS.APPLY_SEARCHEFF_OPT, ctmp);

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
    for ( itmp=1; itmp <= NLINE; itmp++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t %s \n", HOSTLIB.COMMENT[itmp] );
    }
  }


  // -----  FUDGES on observing conditions ------
  readme_doc_FUDGES(&i);

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
  for ( ilist=1; ilist <= NLIST_RAN; ilist++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t FIRST/LAST Random Number (List=%d): %f %f  \n", 
	    ilist, RANFIRST[ilist], RANLAST[ilist] );
  }

  // ---- statistics

  double t_tot   = (t_end-t_start) ; // total proc time, sec
  double R_gen   = (double)NGENLC_TOT / t_tot ;  // NGEN/sec
  double R_write = (double)NGENLC_WRITE/t_tot ;  // NWRITE/sec

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Generation Statistics (total CPU=%.1f minutes): \n", 
	  t_tot/60.);

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

  if ( GENPAR_SELECT.NSTORE > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t %5d rejected by GENPAR_SELECT_FILE \n",  
	    NGEN_REJECT.GENPAR_SELECT_FILE );
  }

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
  if ( GENLC.STOPGEN_FLAG == 1  ) {

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n  WARNING: GENERATION STOPPED WHEN ERROR\n");
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t   on ERR(EFF) <= EFFERR_STOPGEN(=%f)\n",  
	    INPUTS.EFFERR_STOPGEN );

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t  => YOU HAVE %d FEWER LIGHT CURVES THAN REQUESTED\n",
	    INPUTS.NGEN_LC - NGENLC_WRITE );
  }


  // give SN-stats with cuts
  XN   = (INPUTS.RATEPAR.SEASON_COUNT  + INPUTS.RATEPAR_PEC1A.SEASON_COUNT) ;
  XN  *= GENLC.GENEFF ; // multiply by cut-efficiency
  if ( NGENLC_WRITE > 0 ) 
    { XNERR = XN/sqrt((double)NGENLC_WRITE); }
  else
    { XNERR = 0.0 ; }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Number of SNe per season AFTER CUTS : %6.0f +- %5.0f \n", 
	  XN, XNERR );


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
  GENPOLY_DEF *GENZPOLY_SNR, *GENZPOLY_TEXPOSE, *GENLAMPOLY_WARP; 

  // ------- BEGIN -------

  if ( N == 0 ) { return ; }

  i = *iline ;
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  TAKE_SPECTRUM: \n");

  
  for(j=0; j < N; j++ ) {

    ptrEP    = INPUTS.TAKE_SPECTRUM[j].EPOCH_RANGE; 
    ptrLAM   = INPUTS.TAKE_SPECTRUM[j].SNR_LAMRANGE ;

    GENZPOLY_SNR     = &INPUTS.TAKE_SPECTRUM[j].GENZPOLY_SNR ;
    GENZPOLY_TEXPOSE = &INPUTS.TAKE_SPECTRUM[j].GENZPOLY_TEXPOSE ;
    GENLAMPOLY_WARP  = &INPUTS.TAKE_SPECTRUM[j].GENLAMPOLY_WARP ;

    // .xyz
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

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    if ( IS_HOST ) {
      sprintf(cptr,"    %s                      %s-zPOLY:%s  %s\n",
	      Tname, name2, zpolyString, warpString );
    }
    else {
      sprintf(cptr,"    %s = %5.1f to %5.1f  "
	      "  %s-zPOLY:%s  %s\n",
	      Tname, ptrEP[0], ptrEP[1],
	      name2, zpolyString, warpString );
    }

  } // end j loop over spectra

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

  if ( INPUTS.FUDGESCALE_SKYNOISE != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for SIMLIB NOISE(SKY) : %5.2f ", 
	  INPUTS.FUDGESCALE_SKYNOISE);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESCALE_READNOISE != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for SIMLIB NOISE(CCD-read) : %5.2f ", 
	  INPUTS.FUDGESCALE_READNOISE);
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

  int i, j, NSKIP, isp, index, NINDEX ;
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
  char fnam[] = "readme_doc_magSmear" ;

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
  if ( INPUTS.GENMODEL_ERRSCALE > 0.0 ) onoff=1; 
  if ( INPUTS.NFILT_SMEAR       > 0   ) onoff=1; 


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
  if ( istat_genSmear() > 0 ) { 
    onoff=1; 
    sprintf(ctmp,"%s", INPUTS.GENMAG_SMEAR_MODELNAME);
  }
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

  if ( INPUTS.OPTMASK_SIMSED == OPTMASK_SIMSED_GRIDONLY ) {
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
void readme_doc_hostxt(int *iline) {

  // add host extinction info to README 

  int i;
  int LDMP_HOSTXT, USE_AV ;
  char *cptr, cEXPON[40], cGAUSS[40] ;

  // ------------ BEGIN --------------

  i = *iline ;

  LDMP_HOSTXT = 
    ( GENFRAME_OPT   == GENFRAME_REST ) ||
    ( INPUTS.NON1A_MODELFLAG > 0 )  ;

  USE_AV = 
    (INPUTS.GENRANGE_AV[1] > 1.0E-9 ) &&
    (INPUTS.GENEXPTAU_AV   > 1.0E-9  || INPUTS.GENGAUSIG_AV > 1.0E-9 ) ;

  if ( LDMP_HOSTXT && (USE_AV==0)  ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n  Host Extinction Parameters: NONE  (AV=0) \n");    
  }

  if ( LDMP_HOSTXT && USE_AV  ) {

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n  Host Extinction Parameters: \n");

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf_GENGAUSS(cptr, "\t RV ", &INPUTS.GENGAUSS_RV);

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Gen-Range for AV  : %4.2f to %4.2f  (model=%s) \n", 
	    INPUTS.GENRANGE_AV[0], INPUTS.GENRANGE_AV[1], INPUTS.GENSNXT );
  

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    if ( INPUTS.GENEXPTAU_AV < 50.0 ) {

      cEXPON[0] = 0 ;
      cGAUSS[0] = 0 ;
      if ( INPUTS.GENEXPTAU_AV > 1.0E-9 ) 
	{ sprintf(cEXPON,"exp(-AV/%4.2f)/tau", 
		  INPUTS.GENEXPTAU_AV ); 
	}
      if ( INPUTS.GENGAUSIG_AV > 1.0E-9 ) 
	{ sprintf(cGAUSS,"%4.2f x Gauss(AV,sig=%4.2f)",
		  INPUTS.GENRATIO_AV0, INPUTS.GENGAUSIG_AV ); 
	}

      sprintf(cptr,"\t dN/dAv = %s + %s \n", cEXPON, cGAUSS ); 
    }
    else  { 
      sprintf(cptr,"\t dN/dAv = flat \n");  // Rodney request (Oct 2012)
    }

  }  // end of SALT2/color/dust if-block


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

  int i ;
  char *cptr, string[40] ;
  char star[2];

  // ------------ BEGIN ---------

  if ( INDEX_GENMODEL != MODEL_SALT2 ) { return ; }

  if ( strlen(INPUTS.SALT2mu_FILE) > 0 )   { sprintf(star,"*"); }
  else { sprintf(star," "); }

  i = *iline ;

  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf_GENGAUSS(cptr, "\t SALT2c", &INPUTS.GENGAUSS_SALT2c);
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(string, "\t Alpha%s", star);
  sprintf_GENGAUSS(cptr, string, &INPUTS.GENGAUSS_SALT2ALPHA );
  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( INPUTS.SALT2BETA_cPOLY[0] > 0.0 ) { 
    sprintf(cptr,"\t Beta = %.3f + %.3f*c + %.3f*c^2 "
	    , INPUTS.SALT2BETA_cPOLY[0] 
	    , INPUTS.SALT2BETA_cPOLY[1] 
	    , INPUTS.SALT2BETA_cPOLY[2]  );
  }
  else {
    sprintf(string, "\t Beta%s ", star);
    sprintf_GENGAUSS(cptr, string, &INPUTS.GENGAUSS_SALT2BETA);
  }

  if ( strlen(INPUTS.SALT2mu_FILE) > 0 )   { 
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t * Alpha,Beta read from a SALT2mu_FILE. \n");
  }

  *iline = i ;

} // end of readme_doc_SALT2params
 

// **************************************
void sprintf_GENGAUSS(char *string, char *name, 
		      GENGAUSS_ASYM_DEF *genGauss) {

  // write genGauss info to string
  // Mar 16 2015: include new SKEW parameter.
  // Aug 30 2016: check SKEWNORMAL

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


  if( genGauss->SKEWNORMAL[1] != 0.0 ) {
    double *ptrSKEW = genGauss->SKEWNORMAL ;
    sprintf(string,"%s: SKEWNORMAL(MEAN,SIG,SKEW) = %.3f, %.3f, %.3f\n",
	    name, ptrSKEW[0], ptrSKEW[1], ptrSKEW[2] );
  }
  else {
    sprintf(cSKEW,"SKEW=%.2f,%.2f",  
	    genGauss->SKEW[0], genGauss->SKEW[1] );
    
    sprintf(cRANGE,"BND=%.2f,%.2f", 
	    genGauss->RANGE[0], genGauss->RANGE[1] );
    
    sprintf(string,"%s: %s  %s  %s  %s\n"
	    ,name, cPEAK, cSIGMA, cSKEW, cRANGE  );
  }

  return ;

} // sprintf_GENGAUSS




// ***************************************************
void update_accept_counters(void) {

  // Created June 2017
  // Called from main for accppted event, so here we
  // increment various counters.
  // Most code moved from main for cleanup.

  int isp ;
  char fnam[] = "update_accept_counters" ;

  // -------------- BEGIN ---------------

  // increment number of generated SN that are written
  NGENLC_WRITE++ ;
 
  
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
  //
  // Feb 12, 2014: always call snlc_to_SNDATA(1) instead of only
  //               for FITS format.

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

  // create new subdir for simulated SNDATA files

  sprintf(cmd,"mkdir -m g+wr %s", PATH_SNDATA_SIM );
  isys = system(cmd);

  // create full names for auxilliary files,
  // whether they are used or not.

  sprintf(prefix,"%s/%s", PATH_SNDATA_SIM, INPUTS.GENVERSION );

  // mandatory
  sprintf(SIMFILE_AUX->LIST,       "%s.LIST",        prefix );
  sprintf(SIMFILE_AUX->README,     "%s.README",      prefix );
  sprintf(SIMFILE_AUX->IGNORE,     "%s.IGNORE",      prefix );


  // optional
  sprintf(SIMFILE_AUX->DUMP,       "%s.DUMP",        prefix );
  sprintf(SIMFILE_AUX->ZVAR,       "%s.ZVARIATION",  prefix );
  sprintf(SIMFILE_AUX->GRIDGEN,    "%s.GRID",        prefix );


  // create mandatory files.
  SIMFILE_AUX->FP_LIST   = fopen(SIMFILE_AUX->LIST,   "wt") ;  
  SIMFILE_AUX->FP_README = fopen(SIMFILE_AUX->README, "wt") ;  
  SIMFILE_AUX->FP_IGNORE = fopen(SIMFILE_AUX->IGNORE, "wt") ;   
  fclose(SIMFILE_AUX->FP_IGNORE);


  // dump out the README file
  for ( i = 1; i <= VERSION_INFO.NLINE_README_INIT; i++ )
    { fprintf(SIMFILE_AUX->FP_README, "%s", VERSION_INFO.README_DOC[i] ); }

  fflush(SIMFILE_AUX->FP_README);

  // if FITRES DUMP-file is requested, open and init header
  wr_SIMGEN_DUMP(1,SIMFILE_AUX);
  
  if ( GENLC.IFLAG_GENSOURCE == IFLAG_GENGRID  ) {
#ifdef SNGRIDGEN
    init_GRIDsource(1); 
    wr_GRIDfile(1,SIMFILE_AUX->GRIDGEN);
#endif
  }

  snlc_to_SNDATA(1) ;  // 1 => load header only

  // check option for fits format (Jun 2011)
  if ( WRFLAG_FITS ) { 

    // abort of any text-option is defined along with fits format
    if ( WRFLAG_TEXT || WRFLAG_VBOSE ) {
      sprintf(c1err, "Cannot mix TEXT and FITS format ; FORMAT_MASK=%d",
	      INPUTS.FORMAT_MASK );
      sprintf(c2err,"WRFLAG[TEXT,VBOSE,FITS] = %d, %d, %d",
	      WRFLAG_TEXT, WRFLAG_VBOSE, WRFLAG_FITS );
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

  int  NEWMJD, CID    ;
  char fnam[] = "update_simFiles";

  // ------------ BEGIN -------------


#ifdef SNGRIDGEN
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
    sprintf(c1err,"CID=%d exceeds MXCID_SIM=%d", GENLC.CID, MXCID_SIM );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }


  // init SNDATA strucure
  init_SNDATA() ; 

  // load SNDATA structure
  snlc_to_SNDATA(0) ;

  // always fill DUMP file even if SNDATA files are not written
  wr_SIMGEN_DUMP(2,SIMFILE_AUX);

  if ( INPUTS.FORMAT_MASK <= 0 ) { return ; }

  if ( WRFLAG_FITS ) { 
    WR_SNFITSIO_UPDATE(); 
    return ;
  }

  // below are the text-output options

  NEWMJD = SNDATA.NEWMJD ; // save NEWMJD value

  // check option to suppress verbose output
  if ( WRFLAG_VBOSE == 0 ) { SNDATA.NEWMJD = 0; }

      
  // alawys write header, and write NEWMJD epochs if verbose format 
  wr_SNDATA ( INPUTS.WRITE_MASK, 0 );

  SNDATA.NEWMJD = NEWMJD ;   // restore NEWMJD

  // check option to append TERSE output at end of SNDATA file
  if ( WRFLAG_TEXT > 0 ) { 
    append_SNPHOT_TEXT() ; 
    append_SNSPEC_TEXT() ;   // July 2016
  }

  // update LIST file
  fprintf(SIMFILE_AUX->FP_LIST,"%s\n", SNDATA.snfile_output);

} // end of upd_simFiles


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
    // xxx mark delete    fclose(SIMFILE_AUX->FP_DUMP);
    // xxx mark delete    printf("  %s \n", SIMFILE_AUX->DUMP );
  }


  // copy ZVARATION file to SIM/[VERSION]
  if ( USE_ZVAR_FILE ) {
    cp_zvariation(SIMFILE_AUX->ZVAR);  
    printf("  %s\n", SIMFILE_AUX->ZVAR);
  }

#ifdef SNGRIDGEN
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

  int CID;

  if ( WRFLAG_CIDRAN == 0 ) 
    { CID = GENLC.CID ; }
  else
    { CID = GENLC.CIDRAN ; }

  if ( INPUTS.NGEN_LC > 0 ) {

    if ( LUPDGEN(NGENLC_WRITE)  ) {
      printf("\t Finished writing %6d of %d (CID=%8d, NEP=%3d) \n", 
	     NGENLC_WRITE, INPUTS.NGEN_LC, CID, GENLC.NEPOCH );
      fflush(stdout);
    }

  } else {

    if ( LUPDGEN(NGENLC_TOT) ) {
      printf("\t Finished generating %8d of %d (CID=%6d) \n", 
	     NGENLC_TOT, INPUTS.NGENTOT_LC, CID );
      fflush(stdout);
    }

  }

} // end of screen_update



// ***********************************
void append_SNPHOT_TEXT(void) {

  // Created May 31, 2009
  // Re-open SNDATA file and append terse light curve
  // info with "MJD FILT FLUX FLUXERR" per line
  // Useful to pass to other collaborators not using SNANA
  //
  // Oct  24, 2009: fix dumb bug filling fluxcal. Before was
  //                loading flux(in ADU)
  //
  // Dec 29, 2009: indicate epoch when search-trigger is satisfied.
  //               (see GENLC.MJD_TRIGGER)
  //
  // Jun 11, 2010: if sim_mag=99, then set mag=99 and magerr=5 ...
  //               for pre-explosion epochs
  //
  // Apr 28, 2011: write SIM_MAG with f8.4 instead of f8.3
  //
  // May 30, 2011: write FLUXCAL as %11.4le insteadl %10.3le
  //
  // Jan 23 2014: include columns for rest-filt and rest-mag 
  //              (for rest-models only)
  //
  // Mar 3, 2014: replace remaining GENLC with SNDATA elements.
  //
  // Aug 8 2014: SIM_MAG key -> SIM_MAGOBS key
  //
  // Sep 22 2015: write PSF column
  //
  // Mar 8 2016: if template-error is used, write both SIG_SKY and SIG_SKY_T
  //
  // Nov 4 2016: baseline NVAR_TEXT=12 (was 11) to include SIM_MAGOBS
  //
  // Mar 16 2017: use tmpLine buffer to reduce number of writes
  //
  // Oct 07, 2017: remove 0.4 day buffer on writing DETECTION info
  //
  // Jan 17 2018: add permanent PHOTFLAG; remove SNR,MAG,MAG_ERR
  //
  // Mar 13 2018: write PHOTPROB
  // Mar 18 2018: write optional SNR_MON
  // Feb 08 2019: 
  //   + use strcat in instead of sprintf(tmpLine,"%sxx", tmpLine ...)

  FILE *fp;
  int ISMODEL_FIXMAG    = ( INDEX_GENMODEL == MODEL_FIXMAG );

  int  NVAR_TEXT, ep, ifilt_obs, ifilt_rest, ADD_PHOTFLAG=1 ;
  int LTRIGCHECK, IFLAG, PHOTFLAG, APPEND_MAGREST ;
    
  double mjd, flux, fluxerr, searcherr, templerr, mag, magerr;
  double mag_rest, ZPT, PSF, SNR, SIM_MAGOBS, SNR_MON, PHOTPROB ;

  char  cfilt[2], cfield[20], OBSKEY[12], cval[40];
  char  tmpLine[400], tmpVal[80], tmpVarName[40];
  char  fnam[] = "append_SNPHOT_TEXT" ;

  // ------------ BEGIN --------------

  fp = fopen(SNDATA.SNFILE_OUTPUT, "at") ;

  APPEND_MAGREST = ( WRFLAG_BLINDTEST == 0             && 
		     GENFRAME_OPT     == GENFRAME_REST &&
		     !ISMODEL_FIXMAG		     );

  // determine number of variables to write out
  NVAR_TEXT =  8;         // Mar 13 2018

  if ( WRFLAG_BLINDTEST  ) 
    { NVAR_TEXT -= 4 ; } // exclude sim info

  if ( APPEND_MAGREST    ) 
    { NVAR_TEXT += 2 ; } // rest-filter and mag

  if ( SIMLIB_TEMPLATE.USEFLAG ) 
    { NVAR_TEXT += 1 ; } // template error

  if ( ADD_PHOTFLAG )
    { NVAR_TEXT += 1 ; } 

  if ( INPUTS_SEARCHEFF.NMAP_PHOTPROB ) 
    { NVAR_TEXT += 1 ; } 

  if ( INPUTS.MAGMONITOR_SNR  )  // SNR-monitor for fixed mag
    { NVAR_TEXT += 1 ; }

  if ( GENLC.MJD_TRIGGER > 0.99E6 ) 
    { LTRIGCHECK = 0; }
  else 
    { LTRIGCHECK = 1; }

    
  fprintf(fp,"\n\n# ============================================ \n");
  fprintf(fp,"# TERSE LIGHT CURVE OUTPUT: \n#\n");
  fprintf(fp,"NOBS: %d \n", GENLC.NOBS - GENLC.NOBS_REMOVE );
  fprintf(fp,"NVAR: %d \n", NVAR_TEXT );

  // write VARLIST strings into tmpLine memory
  sprintf(tmpLine,"VARLIST:  MJD  FLT FIELD   FLUXCAL   FLUXCALERR ");

  if ( ADD_PHOTFLAG > 0  )
    { strcat(tmpLine,"PHOTFLAG "); }

  if ( INPUTS_SEARCHEFF.NMAP_PHOTPROB ) 
    { strcat(tmpLine,"PHOTPROB "); }

  if ( SIMLIB_TEMPLATE.USEFLAG )  
    { strcat(tmpLine,"SKY_SIG SKY_SIG_T " ); }

  strcat(tmpLine,"  ZPT") ;
  strcat(tmpLine,"  PSF") ;  // added Sep 22 2015


  if ( WRFLAG_BLINDTEST == 0 ) 
    { strcat(tmpLine,"  SIM_MAGOBS"); }

  if ( INPUTS.MAGMONITOR_SNR  ) { 
    sprintf(tmpVarName, "SIM_SNRMAG%2.2d", INPUTS.MAGMONITOR_SNR ); 
    strcat(tmpLine," " );    strcat(tmpLine, tmpVarName );
  }

  if ( APPEND_MAGREST ) 
    { strcat(tmpLine,"  FLT_REST  SIM_MAGREST "); }

  // write VARLIST line to file
  fprintf(fp,"%s\n", tmpLine);

    
      
  for ( ep=1; ep <= GENLC.NEPOCH; ep++ ) {

    if ( GENLC.USE_EPOCH[ep] == 0  ) { continue ; }

    mjd         = GENLC.MJD[ep];
    ifilt_obs   = GENLC.IFILT_OBS[ep] ;
    sprintf(cfilt,  "%c", FILTERSTRING[ifilt_obs]);
    sprintf(cfield, "%s", GENLC.FIELDNAME[ep] ) ;

    flux        = SNDATA.FLUXCAL[ep] ;
    fluxerr     = SNDATA.FLUXCAL_ERRTOT[ep] ;
    templerr    = SNDATA.SKY_SIG_T[ep] ;
    searcherr   = SNDATA.SKY_SIG[ep] ;

    PHOTFLAG    = SNDATA.PHOTFLAG[ep] ;
    PHOTPROB    = SNDATA.PHOTPROB[ep] ; // Mar 13 2018

    mag         = SNDATA.MAG[ep] ;
    magerr      = SNDATA.MAG_ERRPLUS[ep] ;   // errplus = errminus
    SIM_MAGOBS  = SNDATA.SIMEPOCH_MAG[ep] ;   // true mag
    SNR_MON     = SNDATA.SIMEPOCH_SNRMON[ep];

    ZPT         = SNDATA.ZEROPT[ep];
    PSF         = SNDATA.PSF_SIG1[ep];
    SNR         = flux/fluxerr ;

    // check for trigger marker
    if ( LTRIGCHECK && mjd > GENLC.MJD_TRIGGER+1.0E-5  ) {
      fprintf(fp,"DETECTION:  %d obs satisfy %s  at MJD_TRIGGER=%.4f \n", 
	      SEARCHEFF_LOGIC.NMJD, SEARCHEFF_LOGIC.INPUT_STRING,
	      GENLC.MJD_TRIGGER);
      LTRIGCHECK = 0;
    }

    sprintf(OBSKEY,"OBS");      
    
    sprintf(tmpLine, "%s: %10.3f  %s %4s  %11.4le  %10.3le", 
	    OBSKEY, mjd, cfilt, cfield, flux, fluxerr );

    if ( ADD_PHOTFLAG > 0  ) { 
      sprintf(tmpVal," %4d", PHOTFLAG);
      strcat(tmpLine,tmpVal);
    }

    if ( INPUTS_SEARCHEFF.NMAP_PHOTPROB > 0 )  {
      sprintf(tmpVal," %5.3f", PHOTPROB);
      strcat(tmpLine,tmpVal);
    }  

    if ( SIMLIB_TEMPLATE.USEFLAG )  { 
      sprintf(tmpVal," %.3le %.3le", searcherr, templerr ); 
      strcat(tmpLine,tmpVal);
    }

    sprintf(tmpVal,"  %6.3f", ZPT); strcat(tmpLine,tmpVal);
    sprintf(tmpVal,"  %5.2f", PSF); strcat(tmpLine,tmpVal);

    if ( WRFLAG_BLINDTEST == 0 ) {
      sprintf(tmpVal, " %8.4f ", SIM_MAGOBS ); 
      strcat(tmpLine,tmpVal); 
    }

    if ( INPUTS.MAGMONITOR_SNR ) { 
      sprintf(tmpVal, " %6.1f ", SNR_MON );  
      strcat(tmpLine,tmpVal); 
    }

    if ( APPEND_MAGREST ) { 
      mag_rest    = GENLC.genmag8_rest[ep] ;
      ifilt_rest  = GENLC.IFILTMAP_REST1[ifilt_obs] ;
      sprintf(tmpVal," %c %.4f", FILTERSTRING[ifilt_rest], mag_rest );
      strcat(tmpLine,tmpVal) ;
    }

    fprintf(fp, "%s\n", tmpLine);


  }  //  end of 'ep' epoch loop

  //  fprintf(fp, "END: \n");
  fprintf(fp, "END_PHOTOMETRY: \n");
  fflush(fp);
  fclose(fp);

  return ;

} // end of append_SNPHOT_TEXT


// ***************************************
void append_SNSPEC_TEXT(void) {

  // Created July 2016
  // Write spectra to end of text file.
  //
  // output Format
  //
  // NSPECTRA:  <NMJD>
  //
  // SPECTRUM: 1
  // SPECTRUM_MJD:          <mjd>   
  // SPECTRUM_TEXPOSE:      <texpose>
  // SPECTRUM_SNR_COMPUTE:  <computed snr>
  // SPECTRUM_LAMOBS_SNR:   <lamobs range for SNR_COMPUTE>
  // SPECTRUM_NLAM:    <NLAM>
  // SPECTRUM_HOST_CONTAM:  <frac>
  //
  // Mar 24 2019: add SIM_WARP column for mis-calibration

  int FORMAT_LAMCEN = (INPUTS_SPECTRO.FORMAT_MASK & 1) ;
  int    IDSPEC, imjd, ilam, NMJD, NBLAM_TOT, NBLAM_VALID, NBLAM_WR ;
  int    NVAR, IS_HOST ;
  double L0, L1, LCEN, FLAM, FLAMERR, GENFLAM, GENMAG, GENSNR, WARP ;
  FILE *fp ;
  char tmpLine[400], varList_lam[40], varList_warp[20];
  char fnam[] = "append_SNSPEC_TEXT" ;

  // ------------ BEGIN -----------

  if ( INPUTS.SPECTROGRAPH_OPTIONS.DOFLAG_SPEC == 0 ) { return ; }

  fp = fopen(SNDATA.SNFILE_OUTPUT, "at") ;

  NMJD      = GENSPEC.NMJD_TOT ;
  NBLAM_TOT = GENSPEC.NBLAM_TOT ;

  NVAR=4;  varList_lam[0] = varList_warp[0] = 0 ;
  if ( FORMAT_LAMCEN ) {
    // with fixed lambda bins, just record central lambda per bin
    NVAR++;    sprintf(varList_lam,"LAMCEN");
  }
  else {
    // with non-uniform lambda bins, record MIN and MAX per bin
    NVAR+=2 ;  sprintf(varList_lam,"LAMMIN LAMMAX");
  }

  if ( GENSPEC.USE_WARP )  { NVAR++;  sprintf(varList_warp,"SIM_WARP"); } 

  // write header info
  fprintf(fp,"\n# ============================================= \n");
  fprintf(fp,"NSPECTRA:  %d \n\n", NMJD);
  fprintf(fp,
	  "NVAR_SPEC: %d \n"
	  "VARNAMES_SPEC: %s  "
	  "FLAM  FLAMERR   SIM_GENFLAM SIM_GENMAG  %s\n\n",
	  NVAR, varList_lam, varList_warp );

  for(imjd=0; imjd < NMJD; imjd++ ) {
    IDSPEC = imjd + 1 ;  // start at 1
    NBLAM_VALID = GENSPEC.NBLAM_VALID[imjd] ;
    IS_HOST     = GENSPEC.IS_HOST[imjd];

    if ( NBLAM_VALID == 0 ) { return; } // suppress legacy bug (Aug 23 2017)

    fprintf(fp,"SPECTRUM_ID:       %d  \n", IDSPEC ) ; 

    fprintf(fp,"SPECTRUM_MJD:      %9.3f            ", GENSPEC.MJD_LIST[imjd] );
    if ( IS_HOST ) 
      { fprintf(fp, "# HOST \n"); }
    else
      { fprintf(fp, "# Tobs = %8.3f \n", GENSPEC.TOBS_LIST[imjd] ); }    

    fprintf(fp,"SPECTRUM_TEXPOSE:  %9.1f            "
	    "# seconds\n",  
	    GENSPEC.TEXPOSE_LIST[imjd] );

    fprintf(fp,"SPECTRUM_SNR_COMPUTE:  %9.3f        "
	    "# from user SNR request\n",  
	    GENSPEC.SNR_COMPUTE_LIST[imjd] );

    fprintf(fp,"SPECTRUM_LAMOBS_SNR:   %7.1f %7.1f  "
	    "# range for SNR_COMPUTE\n",  
	    GENSPEC.LAMOBS_SNR_LIST[imjd][0],
	    GENSPEC.LAMOBS_SNR_LIST[imjd][1] );

    fprintf(fp,"SPECTRUM_NLAM:     %d (of %d)         "
	    "# Number of valid wave bins\n", 
	    NBLAM_VALID, NBLAM_TOT );
    
    if ( INPUTS.NHOST_TAKE_SPECTRUM > 0 && !IS_HOST ) {
      fprintf(fp,"SPECTRUM_HOSTFRAC:   %.2f               "
	      "# host-frac flux contamination\n", 
	      INPUTS.TAKE_SPECTRUM_HOSTFRAC );
    }

    NBLAM_WR = 0 ;

    for(ilam=0; ilam < NBLAM_TOT; ilam++ ) {
      GENFLAM    = GENSPEC.GENFLAM_LIST[imjd][ilam];
      GENMAG     = GENSPEC.GENMAG_LIST[imjd][ilam];
      FLAM       = GENSPEC.FLAM_LIST[imjd][ilam];
      FLAMERR    = GENSPEC.FLAMERR_LIST[imjd][ilam];
      WARP       = GENSPEC.FLAMWARP_LIST[imjd][ilam];

      if ( FLAMERR <= 0.0 ) { continue ; } // skip unphysical values

      L0      = INPUTS_SPECTRO.LAMMIN_LIST[ilam];
      L1      = INPUTS_SPECTRO.LAMMAX_LIST[ilam];
      LCEN    = INPUTS_SPECTRO.LAMAVG_LIST[ilam];

      if ( FORMAT_LAMCEN ) 
	{  sprintf(varList_lam,"%8.2f", LCEN); }
      else
	{  sprintf(varList_lam,"%8.2f %8.2f", L0, L1); }

      sprintf(tmpLine,
	      "SPEC: %s  %10.3le  %10.3le   %10.3le  %.2f",
	      varList_lam, FLAM, FLAMERR, GENFLAM, GENMAG );

      if ( GENSPEC.USE_WARP ) 
	{ sprintf(tmpLine,"%s %6.3f", tmpLine, WARP); }

      fprintf(fp,"%s \n", tmpLine);
      NBLAM_WR++ ; 
    }

    fprintf(fp,"SPECTRUM_END:  \n\n" );

    if ( NBLAM_WR != NBLAM_VALID ) {
      sprintf(c1err,"Wrote %d lamBins, but expected %d",
	      NBLAM_WR, NBLAM_VALID);
      sprintf(c2err,"CID=%d  SpecID=%d MJD=%.3f", 
	      GENLC.CID, IDSPEC, GENSPEC.MJD_LIST[imjd] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  }  // end imjd loop
  
  fflush(fp);  fclose(fp);

  return ;

} // end append_SNSPEC_TEXT




// ***********************************
void SIMLIB_DUMP_DRIVER(void) {

  // Dump summary information for SIMLIB 
  //
  // History:
  // Jan 29, 2010: move SIMLIB_DUMP struct to snlc_sim.h
  //                to avoid valgrind errors
  //
  // Sep 27, 2010: 
  //   - add few more columns: SKYMAG and M5SIG
  //   - get CADENCEFOM from SNcadenceFoM()
  //   - include FoM_[filt] in .DUMP file
  //
  // Oct 14, 2010: Replace NREAD-1 with NREAD for first argument to 
  //               SIMLIB_angsep_min(...) ; 
  //               Can't recall why NREAD-1 was used ???
  //
  //
  // Jan 13, 2011: remove stupid ';' from ZPT_pe definition so that
  //               ZPT_pe is really in pe instead of ADU
  //
  // Jan 21, 2011: print <GAP> in filter-summary
  // Jan 26, 2011: print SKYSIG in pe/pix instead of ADU/pix
  //
  // Mar 4, 2011: 
  //              There is a bug calling
  //                SIMLIB_angsep_min(NREAD, RA,DECL, RA_STORE, DECL_STORE);
  //              because RA_STORE and DECL_STORE are not set on 1st call.
  //              ==> need to fix at some point, but for now it's not called.
  //
  // Jun 07, 2011: define MJDGAP_IGNORE = 100.0 ;
  //               -> pass to MJDGAP and  to SNcadenceFOM().
  //
  // Aug 8, 2013:  bits0,1 of SIMLIB_DUMP -> dump per LIBID, per MJD
  //
  // Mar 16 2015: change dump-file name to NOT use name of simlib file;
  //              instead use survey-name and filter-list.
  //              See dmpFile0.
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
  // ----------------------------------

#define MXSIMLIB_DUMP_STDOUT 50   // max simlib entries to screen-dump
#define ZPTSIG_MAX  0.05      // exclude entries with this much ZPT-variation

  int 
    ID, IDLAST, NREAD, LDMP_LOCAL
    ,LDMP_SEQ_TEXT, LDMP_OBS_TEXT, LDMP_ROOT
    ,ifilt, ifilt_obs, iep, i, icut, Ncut, NVAR, NROW_MJD, NROW
    ;

  char  cfilt[2];

  FILE *fpdmp0, *fpdmp1 ;
  
  double
    MJD, MJD_LAST, GAPMAX, GAPAVG, MJDWIN, FRAC, *ptrmjd
    ,ZPT_pe, PSF, SKYSIG_ADU, SKYSIG_pe
    ,XNEP, ZPTSIG, M5SIG, FOM
    ,ZPT_SIMLIB, PSF_SIMLIB, SNR_maglimit
    ,FSKY_pe, SKYMAG, RA, DEC, MWEBV
    ,XNobs, TMP, TMP0,  TMP1, wgt_LCLIB, wgtsum_LCLIB=0.0
    ,GLOBAL_RANGE_RA[2]
    ,GLOBAL_RANGE_DEC[2]
    ,GAIN_SIMLIB, PIXSIZE_SIMLIB
    ;


  float  MJDMIN4, MJDMAX4, RA4, DEC4 ;

  int 
    NLIB_ZPT,NLIB_PSF, NLIB_SKYSIG ,NLIB_SKYMAG ,NLIB_M5SIG
    ,NLIB_NEP, NLIB_GAP, NLIB_FOM
    ;


  // SNcadenceFoM args
  double 
     MJDLIST_ALL[MXEPSIM]   // all filters together 
    ,MJDLIST[MXFILTINDX][MXEPSIM_PERFILT]    // MJD list for one lib entry
    ,M5SIGLIST[MXFILTINDX][MXEPSIM_PERFILT]   // idem for M5sigma    
    ,*ptrList
    ;
  
  double  MJDGAP_IGNORE = 50.0 ; // ignore Gaps (days) longer than this

  int Nobs, OPTMASK  ;
  char SIMLIB_IDENTIFIER[40] ;
  char ctmp[40], FIELDNAME[60] ;

#define MXSIMLIB_CUTCHECK 10

  int    NCUT_SIMLIB;
  double SIMLIB_GENRANGE[MXSIMLIB_CUTCHECK][2] ;
  int    NSIMLIB_CUTFAIL[MXSIMLIB_CUTCHECK] ;
  char   SIMLIB_CUTVARNAME[MXSIMLIB_CUTCHECK][20] ;
  float *PTR_SIMLIB_CUTVAR[MXSIMLIB_CUTCHECK][2] ;

  char  fnam[] = "SIMLIB_DUMP_DRIVER" ;

  // ------------ BEGIN  ----------


  LDMP_LOCAL = 1 ;

  if ( INPUTS.SIMLIB_DUMP < 0 ) { return ; }

  // check dump options (Aug 8 2013)

  if ( INPUTS.SIMLIB_DUMP == 0 ) {
    sprintf(c1err,"SIMLIB_DUMP = 0 is no longer valid.");
    sprintf(c2err,"Must set bit0,1 to dump per LIBID, per MJD");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  LDMP_SEQ_TEXT =  LDMP_OBS_TEXT = LDMP_ROOT = 0 ; 
  if ( INPUTS.SIMLIB_DUMP & SIMLIB_DUMPMASK_SEQ ) { LDMP_SEQ_TEXT = 1; }
  if ( INPUTS.SIMLIB_DUMP & SIMLIB_DUMPMASK_OBS ) { LDMP_OBS_TEXT = 1; }

#ifdef USE_ROOT
  if ( INPUTS.SIMLIB_DUMP & 4 ) { LDMP_ROOT  = 1; }
#endif

  print_banner("SIMLIB_DUMP");

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

    fprintf(fpdmp0,"NVAR: %d \n", NVAR );
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
    fprintf(fpdmp1,"NVAR: %d \n", NVAR );
    fprintf(fpdmp1,"VARNAMES: ROW "
	    "LIBID RA DEC MJD BAND ZP_pe SKYMAG PSF M5SIG MJD_DIF\n");
    fprintf(fpdmp1,"\n");
  }



  // stdout table header

  printf("\n  LIBID  MJD-range    NEPOCH(all,%s)    GAPMAX(frac) <GAP> \n", 
	 INPUTS.GENFILTERS );
  printf(" ------------------------------------------------------------- \n");

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

      //   printf(" xxx MJD[iep=%d] = %f \n", iep, MJD );
      ZPTSIG    = SIMLIB_OBS_GEN.ZPTSIG[iep]; 
      if ( ZPTSIG > ZPTSIG_MAX ) { continue;} // exclude extreme variations


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

  int ifilt_obs, ifilt, IVAL;
  double VAL, XN  ;
  char  cfilt[2];
  char fnam[] = "update_SIMLIB_DUMP_AVGALL" ;

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
  char PREFIX[MXPATHLEN], DUMPFILE[MXPATHLEN], openOpt[40] ;
  char fnam[] = "SIMLIB_DUMP_openTable" ;

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
  char fnam[] = "SIMLIB_DUMP_makeTable" ;

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
  sprintf(varName,"ZPTSIG(NOBS):D" );
  SNTABLE_ADDCOL(TABLEID_SIMLIB_DUMP, BLOCKVAR,
		 &SIMLIB_OBS_RAW.ZPTSIG[1], varName, 1 );

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
  // Oct 6, 2009: open genmagerr_xxx file and dump errors.
  // Feb 10, 2010: add TEST_COSANGLE for SIMSED model
  // Apr 19, 2010: fix SALT2 bug by setting GENLC.SALT2alpha[beta]
  // Apr 23, 2010: calculate DM15 and print along with SHAPEPAR string
  // Feb 26, 2013: comment out init_RANLIST()
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
    GENLC.epoch8_rest[1] = (double)INPUTS.GENRANGE_DMPTREST[0] ;
    GENLC.epoch8_obs[1]  = (double)INPUTS.GENRANGE_DMPTREST[0] ;
  }
  else {
    for ( epoch = (int)Tmin; epoch <= (int)Tmax; epoch++ ) {
      GENLC.NEPOCH++ ;   N = GENLC.NEPOCH;      
      if ( N >= MXEPSIM ) {
	sprintf(c1err,"Nepoch =%d exceeds array bound ", N);
	sprintf(c2err,"(MXEPSIM = %d)", MXEPSIM);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }      
      GENLC.epoch8_rest[N] = (double)epoch;
      GENLC.epoch8_obs[N]  = (double)epoch;
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

    ptr_genmag8   = &GENFILT.genmag8_obs[ifilt_obs][0] ;
    ptr_generr8   = &GENFILT.generr8_obs[ifilt_obs][0] ;
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
      Trest8 = GENLC.epoch8_rest[epoch] ;

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

  double H0 = INPUTS.H0 / ( 1.0E6 * PC_km) ;
  double OM = INPUTS.OMEGA_MATTER;
  double OL = INPUTS.OMEGA_LAMBDA;
  double w0 = INPUTS.W0_LAMBDA ;
  char fnam[] = "test_zcmb_dLmag_invert" ;
  double MU, zCMB;
  // ----------- RETURN ------------
  for(MU=32.0; MU < 49.0; MU+=1.0 ) {
    zCMB = zcmb_dLmag_invert(H0, OM, OL, w0, MU); 
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
  double x, y, xx, yy, sig, xlen ;
  double x0 = 0.0 ;
  // ------------- BEGIN --------
  for ( i = 1; i<=100000; i++ ) {
    // now smear with sigma=1
    NTMP++ ;
    if ( NTMP==100 ) {  init_RANLIST();  NTMP = 0 ; }
    x = GaussRan(1);
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
  char fnam[] = "test_igm";

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
