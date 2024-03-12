/******************************************************
  
   Created Dec 22 2021 [code moved out of snlc_sim]
   Write DOCUMENTATION block to [VERSION].README    

   Jan 17 2022: remove legacy functions

******************************************************/

#include "sntools.h"
#include "sntools_cosmology.h"
#include "sntools_stronglens.h"
#include "snlc_sim.h"
#include "sntools_host.h"
#include "sntools_wronghost.h"
#include "sntools_nonlinearity.h"
#include "sntools_genSmear.h"
#include "sntools_trigger.h" 
#include "sntools_genPDF.h"
#include "sntools_sim_readme.h"

// *************************************************
void README_DOCANA_DRIVER(int iflag_readme) {

  // Created Dec 22 2021 by R.Kessler
  //
  // Complete refactor for output README.
  // Replace casual comments with yaml-compliant format
  // under DOCUMENTATION block in several parts:
  // 1. OVERVIEW        # survey, model, host, user ...
  // 2. INPUTS_KEYS     # summarize all user input keys
  // 3. INPUTS_NOTES    # compute rates, etc ...
  // 4. OUTPUT_SUMMARY  # stats, cpu time ...
  //
  // Input arg:
  //   iflag_readme = 0 => one-time init 
  //   iflag_readme = 1 => write first 3 parts
  //   iflag_readme = 2 => write OUTPUT_SUMMARY
  //
  // The iflag_readme scheme enables viewing most of the README
  // while job is running; e.g., can check INPUT_KEYS before 
  // long job finishes.
  //
  // Jun 21 2023: abort if any README_DOC line is too long

  int i;
  char fnam[] = "README_DOCANA_DRIVER";

  // ----------- BEGIN -----------

  if ( iflag_readme == 0 ) {
    VERSION_INFO.NLINE_README      = 0;
    VERSION_INFO.NLINE_README_INIT = 0;

    README_KEYPLUSARGS_init(&README_KEYS_COSMO);
    README_KEYPLUSARGS_init(&README_KEYS_GENMODEL);
    README_KEYPLUSARGS_init(&README_KEYS_CID);
    README_KEYPLUSARGS_init(&README_KEYS_SIMLIB);
    README_KEYPLUSARGS_init(&README_KEYS_HOSTLIB);
    README_KEYPLUSARGS_init(&README_KEYS_RATEMODEL);
    README_KEYPLUSARGS_init(&README_KEYS_LENS);
    README_KEYPLUSARGS_init(&README_KEYS_SKY);
    README_KEYPLUSARGS_init(&README_KEYS_MWEBV);
    README_KEYPLUSARGS_init(&README_KEYS_NON1ASED);
    README_KEYPLUSARGS_init(&README_KEYS_SIMSED);
    README_KEYPLUSARGS_init(&README_KEYS_LCLIB);
    README_KEYPLUSARGS_init(&README_KEYS_FILTER);
    README_KEYPLUSARGS_init(&README_KEYS_FLUXERRMODEL);
    README_KEYPLUSARGS_init(&README_KEYS_GENMAG_OFF);
    README_KEYPLUSARGS_init(&README_KEYS_GENMAG_SMEAR);
    README_KEYPLUSARGS_init(&README_KEYS_TAKE_SPECTRUM);
    README_KEYPLUSARGS_init(&README_KEYS_RANSYSTPAR);
    README_KEYPLUSARGS_init(&README_KEYS_ZVARIATION);
    README_KEYPLUSARGS_init(&README_KEYS_GRIDGEN);
    README_KEYPLUSARGS_init(&README_KEYS_CUTWIN);
    README_KEYPLUSARGS_init(&README_KEYS_COVMAT_SCATTER);
    README_KEYPLUSARGS_init(&README_KEYS_SIMGEN_DUMP);
    return;
  }


  i = VERSION_INFO.NLINE_README ;

  sprintf(BANNER,"%s: Prepare README content (iflag=%d)", 
	  fnam, iflag_readme);
  print_banner(BANNER);

  if ( iflag_readme == 1 ) { 
    i++; 
    sprintf(VERSION_INFO.README_DOC[i] ,"%s", KEYNAME_DOCANA_REQUIRED ); 

    README_DOCANA_OVERVIEW(&i);
    README_DOCANA_INPUT_KEYS(&i);
    README_DOCANA_INPUT_NOTES(&i);
  }
  else {
    README_DOCANA_SED_TRUE(&i); // write mag-diff warnings

    README_DOCANA_GENTYPE_MAP(&i);  // write map of GENTYPE -> TYPE_NAME

    README_DOCANA_OUTPUT_SUMMARY(&i);

    i++; 
    sprintf(VERSION_INFO.README_DOC[i],"%s",KEYNAME2_DOCANA_REQUIRED); 
  }

  if ( i > MXDOCLINE ) {    
    sprintf(c1err,"%d DOCANA lines exceeds bound of MXDOCLINE=%d", i);
    sprintf(c2err,"Consider increasing MXDOCLINE");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  VERSION_INFO.NLINE_README = i;  
  if ( iflag_readme == 1 ) 
    { VERSION_INFO.NLINE_README_INIT = i; }


  // - - - - - - - - - - -                                                      
  // abort if any string length is too long (Jun 2023) 
  char *line; int len;
  for(i=0; i < VERSION_INFO.NLINE_README; i++ ) {
    line = VERSION_INFO.README_DOC[i];
    len  = strlen(line);
    if ( len > 2*MXPATHLEN ) {
      print_preAbort_banner(fnam);
      printf(" README_DOC line =\n\t '%s' \n", line );
      sprintf(c1err,"DOCANA line len=%d exceeds bound of 2*MXPATHLEN=%d",
              len, 2*MXPATHLEN);
      sprintf(c2err,"Reduce line len or increase MXPATHLEN");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;
    }
  }

  return;

} // end README_DOCANA_DRIVER


void README_DOCANA_OVERVIEW(int *iline) {

  // Prepare OVERVIEW block for README.
  // Feb 21 2022: write SNIA flag as True or False
  // Mar 31 2022: write TIME_START (UT)

  int i = *iline;
  char pad[] = "    ", *cptr, cwd[MXPATHLEN], TF[12];
  char TIME_START[100] ;
  char *SURVEY         = GENLC.SURVEY_NAME;
  char *SUBSURVEY_LIST = SIMLIB_GLOBAL_HEADER.SUBSURVEY_LIST ;
  char fnam[] = "README_DOCANA_OVERVIEW";

  // ----------- BEGIN ------------

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  %s:", DOCANA_OVERVIEW ); 

  get_TIME_START_readme_docana(TIME_START);
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sTIME_START:   %s",  pad, TIME_START);  
    
  // - - - - -
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sSURVEY:       %s",  pad, SURVEY);

  if ( !IGNOREFILE(SUBSURVEY_LIST)  ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sSUBSURVEY_LIST:  %s",  pad, SUBSURVEY_LIST);
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sGENMODEL:     %s", pad, INPUTS.GENMODEL);

  // Feb 21 2022: write python True or False logical for SNIA (or not
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( LGEN_SNIA ) { sprintf(TF,"True"); }  else{sprintf(TF,"False"); }
  sprintf(cptr,"%sSNIA:         %s", pad, TF); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sHOST_MACHINE: %s", pad, getenv("HOSTNAME") );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sUSERNAME:     %s",  pad, getenv("USER") );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sSNDATA_ROOT:  %s", pad, PATH_SNDATA_ROOT );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sSNANA_DIR:     %s", pad, PATH_SNANA_DIR );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sSNANA_VERSION: %s", pad, SNANA_VERSION_CURRENT );

  // write current directory (Sep 5 2013)
  if ( getcwd(cwd,MXPATHLEN) != NULL ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sCWD:   %s", pad, cwd );
  }

  int ifile;
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sINPUT_FILE:", pad );
  for(ifile=0; ifile < INPUTS.NREAD_INPUT_FILE; ifile++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    // printf(" xxx %s: inpfile='%s'\n", fnam,INPUTS.INPUT_FILE_LIST[ifile]);
    ENVrestore(INPUTS.INPUT_FILE_LIST[ifile], ORIG_FILE_README);
    sprintf(cptr,"%s- %s", pad, ORIG_FILE_README );
  }

  *iline = i;
  return ;
} // end README_DOCANA_OVERVIEW


void get_TIME_START_readme_docana(char *TIME_START) {

  // Created Mar 31 2022
  // return TIME_START string with either
  //  + override valude from $TIME_START env
  //  + actual current time

  char *tptr = getenv("SNANA_TIME_START");

  if ( tptr != NULL ) {
    // use ENV to sync times for many jobs
    sprintf(TIME_START, "%s", tptr );
  }
  else if ( strlen(INPUTS.TIME_START) > 0 ) {
    sprintf(TIME_START,"%s", INPUTS.TIME_START);
  }
  else {
    // compute actual current time
    // note that month is 0-11, so add 1 to tm_mon
    time_t tkey ;    struct tm * ptm;
    time(&tkey);     ptm = gmtime ( &tkey );
    sprintf(TIME_START,"%4.4d-%2.2d-%2.2d  %2.2d:%2.2d  # UT",
	    1900+ptm->tm_year, ptm->tm_mon+1, ptm->tm_mday ,
	    ptm->tm_hour, ptm->tm_min );
  }

  return ;

} // end get_TIME_START_readme_docana


void  README_DOCANA_INPUT_KEYS(int *iline) {
  int i = *iline;
  char *cptr, pad[] = "    ", noComment[]="" ;
  int nval1=1, nval2=2, lenkey=20;
  double *dptr, dval, dval_list[10];
  char fnam[] = "README_DOCANA_INPUT_KEYS";

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  %s:", DOCANA_INPUT_KEYS);

  // output genversion, format, nevt ...
  readme_docana_output(&i, pad);

  // source model
  readme_docana_genmodel(&i, pad);

  // instrument: filter, kcor, simlib, noise ...
  readme_docana_instr(&i, pad);

  // hostlib 
  readme_docana_hostlib(&i,pad);

  // search eff maps
  readme_docana_searcheff(&i, pad);

  // redshift and vpec 
  readme_docana_redshift(&i,pad);

  // epoch info (MJD range, Trest range .. )
  readme_docana_epoch(&i,pad);

  // - - - - -  MWEBV - - - - - -
  readme_docana_mwebv(&i,pad);

  // population params for color & stretch
  readme_docana_modelPar(&i,pad);

  // rate model (DNDZ, DNDB ...)
  readme_docana_rate(&i, pad) ;

  // misc
  readme_docana_misc(&i,pad);

  // CUTWIN cut-windows
  readme_docana_cutwin(&i, pad) ;

  *iline = i;
  return;
} // end README_DOCANA_INPUT_KEYS

void README_DOCANA_INPUT_NOTES(int *iline) {
  int  OVP, j, i = *iline;
  char *cptr, onoff[8], pad[] = "  ", dash[]="  -";

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "");
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%s:", pad, DOCANA_INPUT_NOTES );  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  onoff_readme_docana(INPUTS.SMEARFLAG_FLUX,onoff);
  sprintf(cptr, "%s Poisson noise is %s ", dash, onoff);

  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s Reported flux-uncertainty includes ", dash ); 
  OVP = (INPUTS.SMEARFLAG_FLUX & 2);
  if ( OVP == 0 ) 
    { strcat(cptr,"SKY+GALAXY+SOURCE"); }
  else
    { strcat(cptr,"SKY only"); } 
  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  OVP = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_IMAGE ; 
  onoff_readme_docana(OVP,onoff);
  sprintf(cptr, "%s SB-dependent flux scatter is %s",dash, onoff);

  for ( j=0; j < NLINE_RATE_INFO; j++ ) {
     i++; cptr = VERSION_INFO.README_DOC[i] ;
     sprintf(cptr,"%s %s", dash, LINE_RATE_INFO[j] );
  }


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s H0 = %6.2f km/s/Mpc ", dash, INPUTS.H0 );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s Omega_{M,L} = %6.3f, %.3f     w0,wa = %5.2f,%5.3f",
	  dash, INPUTS.OMEGA_MATTER, INPUTS.OMEGA_LAMBDA, 
	  INPUTS.w0_LAMBDA, INPUTS.wa_LAMBDA );

  for(j=0; j < 2; j++ )  {
    char *COMMENT = COMMENT_README_SEARCHEFF[j];
    if ( strlen(COMMENT) > 0 ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s %s ", dash, COMMENT);
    }
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s PIPELINE requires %d detections among %s "
	  "(MJD dif > %.4f days) ", 
	  dash, SEARCHEFF_LOGIC.NMJD, SEARCHEFF_LOGIC.INPUT_STRING,
	  INPUTS.NEWMJD_DIF);


  // write warning for option to read entire SIMLIB once and not rewind.
  int  MSKOPT = SIMLIB_MSKOPT_QUIT_NOREWIND;
  bool QUIT_NOREWIND=( (INPUTS.SIMLIB_MSKOPT & MSKOPT)>0 );
  
  if ( QUIT_NOREWIND ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,
	    "%s WARNING: STOP GENERATION AFTER ONE PASS THRU SIMLIB "
	    "(MSKOPT+=%d)",  dash, MSKOPT);
  }
 
  *iline = i;
  return;
} // end README_DOCANA_INPUT_NOTES


// =========================================
void README_DOCANA_SED_TRUE(int *iline) {

  // Created Oct 2023
  // For option to generate true SED for each observation (sim-input LAMBIN_SED_TRUE),
  // write out alarm values to readme; fraction of observations (per band) in which
  // synthetic mag from SED is discrepant from original genmag using fine-binned SED.

  int  OVP, j, i = *iline;
  int    ifilt, ifilt_obs, N0, N1;
  double frac;

  double LAMBIN_SED_TRUE = INPUTS.SPECTROGRAPH_OPTIONS.LAMBIN_SED_TRUE;
  double magdif_alarm = INPUTS.SPECTROGRAPH_OPTIONS.ALARM_PAR_SED_TRUE[0];
  double mag_alarm    = INPUTS.SPECTROGRAPH_OPTIONS.ALARM_PAR_SED_TRUE[1];

  char cfilt[2];
  char *cptr, pad[] = "  " ;
  char fnam[] = "README_DOCANA_SED_TRUE";

  // ----------- BEGIN ------------

  if ( LAMBIN_SED_TRUE < 0.01 ) { return; }

  readme_docana_comment(&i, "");
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sWARNINGS_SED_TRUE:   "
	  "# fraction of genmag<%.1f with |synmag-genmag|>%.2f", 
	  pad, mag_alarm, magdif_alarm );  

  for(ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = INPUTS.IFILTMAP_OBS[ifilt];
    N0 =  INPUTS.SPECTROGRAPH_OPTIONS.N_ALARM_SED_TRUE[ifilt_obs][0];
    N1 =  INPUTS.SPECTROGRAPH_OPTIONS.N_ALARM_SED_TRUE[ifilt_obs][1];
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );

    if ( N0 > 0 ) 
      { frac = (double)N1 / (double)N0; }
    else
      { frac = 0.0 ; }

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s  %s:  %8.5f   #  %6d / %6d", pad, cfilt, frac, N1, N0);
    
  }
 
  *iline = i;
  return;

} // end README_DOCANA_SED_TRUE


// ===================================================
void README_DOCANA_GENTYPE_MAP(int *iline) {

  // Created Mar 2024
  // Write the following to readme docana section:
  // GENTYPE_TO_NAME:
  //   10:  Ia     SNIa  scale
  //   20:  nonIa  II    scale
  //   32:  nonIa  Ibc   scale
  //   etc ... 
  // Motivation is for classifiers (e.g.. scone, snn) to 
  // automatically read this map from sim-output readme

  int   GENTYPE, i = *iline;
  double SCALE;
  char *TYPE_NAME, *cptr, *STR;
  char STR_Ia[]    = "Ia";
  char STR_notIa[] = "nonIa" ;
  char fnam[]=  "README_DOCANA_GENTYPE_MAP";

  // --------- BEGIN ---------

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"#" );  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  GENTYPE_TO_NAME:  # GENTYPE:  [non]Ia   type  ");  
  
  for (GENTYPE=0; GENTYPE < MXNON1A_TYPE; GENTYPE++ ) {
    TYPE_NAME = INPUTS.GENTYPE_TO_NAME_MAP[GENTYPE];

    SCALE = INPUTS.RATEPAR.DNDZ_ALLSCALE ;
    if ( LGEN_SNIA ) {
      STR = STR_Ia ; 
      SCALE *= INPUTS.RATEPAR.DNDZ_SCALE[0];  // doesn't work with submit_batch_jobs
    } 
    else {
       STR = STR_notIa; 
       SCALE *= INPUTS.RATEPAR.DNDZ_SCALE[0];
    } 

    if ( strlen(TYPE_NAME) > 0 ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"    %d:  %-6s   %s", 
	      GENTYPE, STR, TYPE_NAME );      
    }
  }

  *iline = i;

  return;

} // end README_DOCANA_GENTYPE_MAP

// ========================================
void README_DOCANA_OUTPUT_SUMMARY(int *iline) {
  int  OVP, j, i = *iline;
  char *cptr, *onoff, pad[] = "    ", dash[]="    -";
  char comment[60];
  char *SUBSURVEY_LIST = SIMLIB_GLOBAL_HEADER.SUBSURVEY_LIST; // comma-sep list
  double XN, XNERR;
  double NGEN_PER_SEASON=0.0, NACC_PER_SEASON=0.0, NACCERR_PER_SEASON=0.0;
  char fnam[] = "README_DOCANA_OUTPUT_SUMMARY" ;

  // ------------- BEGIN ------------

  readme_docana_comment(&i, "");
 
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  %s:", DOCANA_OUTPUT_SUMMARY );  

  // first and last random number per random list
  int ilist;
  sumstat_RANLISTs(2);
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sRANDOM_SYNC: ", pad);
  for ( ilist=1; ilist <= GENRAN_INFO.NLIST_RAN; ilist++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s List=%d  FIRST=%f  LAST=%f   "
	    "AVG(wrap) = %.1f +_ %.1f ", 
	    dash, ilist, 
	    GENRAN_INFO.RANFIRST[ilist], GENRAN_INFO.RANLAST[ilist],
	    GENRAN_INFO.NWRAP_AVG[ilist], GENRAN_INFO.NWRAP_RMS[ilist]	);
  }

  readme_docana_comment(&i, "");

  // ---- statistics
  double t_gen, R_gen=0.0, R_write=0.0 ;

  t_gen   = (TIMERS.t_end - TIMERS.t_end_init); // total time after init
  
  if ( t_gen > 0.0 ) {
    R_gen   = (double)NGENLC_TOT / t_gen ;  // NGEN/sec
    R_write = (double)NGENLC_WRITE/t_gen ;  // NWRITE/sec
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sCPU_MINUTES:       %.2f  ",  pad, t_gen/60.);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sNGENLC_TOT:        %d    # (%.f/sec)", 
	  pad, NGENLC_TOT, R_gen );

  if ( GENSL.NGENLC_LENS_TOT > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sNGENLC_LENS_TOT:  %d  ", 
	  pad, GENSL.NGENLC_LENS_TOT );

  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sNGENLC_WRITE:      %d    # (%.f/sec)", 
	  pad, NGENLC_WRITE, R_write );

  // Jan 2022
  // if there are sub-surveys (see global SIMLIB header), then 
  // write stats for each sub-survey
  if ( !IGNOREFILE(SUBSURVEY_LIST) ) {
    int NTOT, NWR, ID, n_subsurvey;  
    char **subsurvey_tmpList, *s, skey[60];
    parse_commaSepList("SUBSURVEY_LIST", SUBSURVEY_LIST, MXIDSURVEY,60,
		       &n_subsurvey, &subsurvey_tmpList);

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sNGENLC_SUBSURVEY:", pad);
    for(j=0; j < n_subsurvey; j++ ) {
      comment[0] = 0;
      if(j==0) { sprintf(comment,"# NTOT NWRITE"); }
      s    = subsurvey_tmpList[j] ;
      sprintf(skey, "%s:", s) ; // e.g., 'CSP:'
      ID   = get_IDSURVEY(s);  // int index in $SNDATA_ROOT/SURVEY.DEF
      NTOT = NGENLC_TOT_SUBSURVEY[ID];
      NWR  = NGENLC_WRITE_SUBSURVEY[ID];
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s  %-12s  %5d %5d   %s", 
	      pad, skey, NTOT, NWR, comment);   
    }
  }


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sNGENSPEC_WRITE:    %d  ", 
	  pad, NGENSPEC_WRITE );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sEFF(SEARCH+CUTS): %7.4f +- %7.4f", 
	  pad, GENLC.GENEFF, GENLC.GENEFFERR);


  // --- number of events per season: calc and accepted
  if ( NLINE_RATE_INFO > 0 ) {
    NGEN_PER_SEASON =  
      INPUTS.RATEPAR.SEASON_COUNT + INPUTS.RATEPAR_PEC1A.SEASON_COUNT ;

    NACC_PER_SEASON   = NGEN_PER_SEASON * GENLC.GENEFF ;
    if ( NGENLC_WRITE > 0 ) 
      { NACCERR_PER_SEASON = NACC_PER_SEASON/sqrt((double)NGENLC_WRITE); }   
  
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sNGEN_PER_SEASON:   %.0f       "
	    "# NSN(GEN) in GENRANGE(z,MJD,dOmega)", pad, NGEN_PER_SEASON );

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sNACC_PER_SEASON:   %.0f +_ %.0f  "
	    "# NSN(ACCEPT) after trigger+cuts", 
	    pad, NACC_PER_SEASON, NACCERR_PER_SEASON);
  }

  // - - - - - 
  // Reject stats


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sNREJECT:  [%d,%d,%d,  %d,%d]   "
	  "# [NEP<%d,GENRANGE,PEAKMAG,  SEARCHEFF,CUTWIN] ", pad,
	  NGEN_REJECT.NEPOCH, 
	  NGEN_REJECT.GENRANGE,
	  NGEN_REJECT.GENMAG,
	  NGEN_REJECT.SEARCHEFF,
	  NGEN_REJECT.CUTWIN,
	  (int)INPUTS.CUTWIN_NEPOCH[0] ) ;

  // check for wrong host info
  if ( !IGNOREFILE(INPUTS.WRONGHOST_FILE) ) {
    int N_WRONGHOST=0 ;  float FRAC=0.0 ;

    N_WRONGHOST = GENLC.NTYPE_PHOT_WRONGHOST;
    if ( NGENLC_WRITE>0 ) { FRAC = (float)N_WRONGHOST/ (float)NGENLC_WRITE; }

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sNWRONGHOST_WRITE:   %d  "
	    "  # frac = %.4f",  pad, N_WRONGHOST, FRAC );
  }


  // - - - - -
  // Jun 29 2022
  // write fractions of host-less and multi-host events in redshift bins
  if ( INPUTS.HOSTLIB_USE && NGENLC_WRITE > 0 )
    { readme_docana_hostmatch(&i, pad); }


  *iline = i;
  return;

} // end README_DOCANA_OUTPUT_SUMMARY


// ========================================
void readme_docana_hostmatch(int *iline, char *pad) {
  
  // Created Jun29 2022
  // print host-nmatch fractions for nmatch=0 and nmatch=2

  int i = *iline;
  int NBINz = WRITE_HOSTMATCH.NRANGE_REDSHIFT;

  int  N, N0,N2, iz; 
  char str_n0_list[100], str_n2_list[100];
  char str_frac0_list[100], str_frac2_list[100], str_z_list[100];
  char str_z[20], str_tmp[20], *cptr ;
  double frac0, frac2, z0, z1;

  char fnam[] = "readme_docana_hostmatch" ;

  // ------------ BEGIN ----------

  str_n0_list[0]    = str_n2_list[0] = 0 ;
  str_frac0_list[0] = str_frac2_list[0] = str_z_list[0] = 0 ;

  for(iz=0; iz < NBINz; iz++ ) {
    N     = WRITE_HOSTMATCH.NGENLC[iz]; 
    N0    = WRITE_HOSTMATCH.NGENLC_NO_HOST[iz]; 
    N2    = WRITE_HOSTMATCH.NGENLC_MULTI_HOST[iz]; 
    z0    = WRITE_HOSTMATCH.REDSHIFT_RANGE[iz][0] ;
    z1    = WRITE_HOSTMATCH.REDSHIFT_RANGE[iz][1] ;

    frac0 = frac2 = 0.0 ;
    if ( N > 0 ) {
      frac0 = (double)N0 / (double)N ;
      frac2 = (double)N2 / (double)N ;
    }

    /*
    printf(" xxx %s: N=%d  N0=%d N2=%d N_WRITE=%d frac[0,2] = %f, %f \n",
	   fnam, N, N0, N2, NGENLC_WRITE, frac0, frac2); fflush(stdout);
    */

    sprintf(str_z,  " %.2f-%.2f", z0, z1);
    catVarList_with_comma(str_z_list, str_z);   

    sprintf(str_tmp," %3d", N0);
    catVarList_with_comma(str_n0_list, str_tmp);
    
    sprintf(str_tmp," %3d", N2);
    catVarList_with_comma(str_n2_list, str_tmp);

    sprintf(str_tmp," %.3f", frac0);
    catVarList_with_comma(str_frac0_list, str_tmp);
    
    sprintf(str_tmp," %.3f", frac2);
    catVarList_with_comma(str_frac2_list, str_tmp);
  }
  
  // - - - - - -
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sHOST_NMATCH:   # accepted", pad );
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s- 0 : [ %s ]  # z-ranges: %s", 
	  pad, str_n0_list, str_z_list);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s- 2 : [ %s ] ", 
	  pad, str_n2_list);

  // - - - - - -
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sHOST_NMATCH_FRACTIONS:  # accepted ", pad );
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s- 0 : [ %s ]  # z-ranges: %s", 
	  pad, str_frac0_list, str_z_list);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s- 2 : [ %s ] ", 
	  pad, str_frac2_list);

  // - - - -

  *iline = i ;
  return ;
} // end readme_docana_hostmatch

// ========================================
void onoff_readme_docana(int FLAG, char *onoff) {
  char STR_ON[]  = "ON";
  char STR_OFF[] = "OFF";
  char *STR ;
  if ( FLAG == 0 )     { STR = STR_OFF; }
  else                 { STR = STR_ON ; }    

  sprintf(onoff,"%s", STR);
  return;

} // end onoff_readme_docana


void readme_docana_comment(int *iline, char *comment) {
  int i = *iline;
  i++; sprintf(VERSION_INFO.README_DOC[i],"# %s", comment);
  *iline = i;
} // end readme_docana_comment


void readme_docana_genmodel(int *iline, char *pad) {
  
  // store GENMODEL, intrinsic scatter model etc ...
  // Sep 19 2022: force writing cosmology parameter keys

  int i = *iline;
  char *cptr, noComment[]="" ;
  int nval1=1, nval2=2, lenkey=24, o ;  
  double dval, *dptr;

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Source model");

  readme_docana_load_list(&i, pad, &README_KEYS_GENMODEL);

  readme_docana_load_list(&i, pad, &README_KEYS_NON1ASED);

  readme_docana_load_list(&i, pad, &README_KEYS_SIMSED);

  readme_docana_load_list(&i, pad, &README_KEYS_LCLIB);

  dval = (double)INPUTS.GENMAG_OFF_GLOBAL ;
  VERSION_INFO_load(&i, pad, "GENMAG_OFF_GLOBAL:", noComment,
		    lenkey, false, nval1, &dval, -20.0,20.0, 0.0); 

  dval = (double)INPUTS.GENMAG_SMEAR[0] ;
  VERSION_INFO_load(&i, pad, "GENMAG_SMEAR:", "coherent scatter (mag)",
		    lenkey, false, nval1, &dval, -20.0,20.0, 0.0); 

  dval = (double)INPUTS.GENMODEL_ERRSCALE ;
  VERSION_INFO_load(&i, pad, "GENMODEL_ERRSCALE:", noComment,
		    lenkey, false, nval1, &dval, 0.0, 20.0, 0.0); 

  readme_docana_load_list(&i, pad, &README_KEYS_GENMAG_OFF);
  readme_docana_load_list(&i, pad, &README_KEYS_GENMAG_SMEAR);

  // - - - - - - -
  readme_docana_comment(&i, "Cosmology inputs");
  dval = (double)INPUTS.OMEGA_MATTER ;
  VERSION_INFO_load(&i, pad, "OMEGA_MATTER:", noComment,
		    lenkey, false, nval1, &dval, 0.0, 2.0, -9.0);
  dval = (double)INPUTS.OMEGA_LAMBDA ;
  VERSION_INFO_load(&i, pad, "OMEGA_LAMBDA:", noComment,
		    lenkey, false, nval1, &dval, 0.0, 2.0, -9.0);
  dval = (double)INPUTS.w0_LAMBDA ;
  VERSION_INFO_load(&i, pad, "w0_LAMBDA:", noComment,
		    lenkey, false, nval1, &dval, -3.0, 3.0, -9.0);
  dval = (double)INPUTS.wa_LAMBDA ;
  VERSION_INFO_load(&i, pad, "wa_LAMBDA:", noComment,
		    lenkey, false, nval1, &dval, -3.0, 3.0, -9.0);

  dval = (double)INPUTS.MUSHIFT ;
  VERSION_INFO_load(&i, pad, "MUSHIFT:", noComment,
		    lenkey, false, nval1, &dval, -3.0, 3.0, -9.0);

  // - - - - 
  // check extra cosmology keys
  if ( README_KEYS_COSMO.NKEY > 0 ) {
    readme_docana_load_list(&i, pad, &README_KEYS_COSMO);
  }

  if ( README_KEYS_LENS.NKEY > 0 ) {
    readme_docana_load_list(&i, pad, &README_KEYS_LENS);
  }

  *iline = i;
  return;

} // end readme_docana_genmodel


// ========================================
void readme_docana_instr(int *iline, char *pad) {

  // instrumental: filters, kcor/calib file, simlib, noise
  int i = *iline;
  int nval1=1, nval2=2, lenkey = 20;
  double dval;
  char noComment[] = "";
  char *cptr;
  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Instrumental inputs");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s %s ",
	  pad, lenkey, "GENFILTERS:", INPUTS.GENFILTERS);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS.KCOR_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s %s ",
	  pad, lenkey, "KCOR_FILE:", ORIG_FILE_README);

  readme_docana_load_list(&i, pad, &README_KEYS_SIMLIB);

  dval = (double)INPUTS.SMEARFLAG_FLUX ;
  VERSION_INFO_load(&i, pad, "SMEARFLAG_FLUX:", "1->add Poisson noise",
		    lenkey, true, nval1, &dval, 0.0,100.0, -9.0); 

  dval = (double)INPUTS.SMEARFLAG_ZEROPT ;
  VERSION_INFO_load(&i, pad, "SMEARFLAG_ZEROPT:", 
		    "+=1->apply scatter, +=2->add to FLUXERRCAL", 
		    lenkey, true, nval1, &dval, 0.0,100.0, -9.0); 

  dval = (double)INPUTS.FUDGE_SNRMAX ;
  VERSION_INFO_load(&i, pad, "FUDGE_SNRMAX:", noComment,
		    lenkey, true, nval1, &dval, 0.0,1.0E5, -9.0); 

  readme_docana_load_list(&i, pad, &README_KEYS_FILTER);

  readme_docana_load_list(&i, pad, &README_KEYS_FLUXERRMODEL);

  if ( README_KEYS_TAKE_SPECTRUM.NKEY > 0 ) {
    readme_docana_comment(&i, "Spectrograph inputs");
    readme_docana_load_list(&i, pad, &README_KEYS_TAKE_SPECTRUM);
  }

  *iline = i;
  return;

} // end readme_docana_instr


void readme_docana_hostlib(int *iline, char *pad) {

  // Created Dec 22 2021
  // Write HOSTLIB_XXX keys and values to DOCUMENTATION block
  //
  // TO DO:
  //  HOSTLIB_FIXSERSIC ??

  int i = *iline;
  int lenkey = 24; // to format printing of key
  double dval, dval_list[20] ;
  char *cptr, *sval, *key;
  int n, nval;
  char  noComment[] = "" ;
  char  fnam[] = "readme_docana_hostlib" ;

  // ----------- BEGIN ------------

  if ( !IGNOREFILE(INPUTS.WRONGHOST_FILE) ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    ENVrestore(INPUTS.WRONGHOST_FILE, ORIG_FILE_README);
    sprintf(cptr,"%s%-*s  %s ", 
	    pad, lenkey, "WRONGHOST_FILE:", ORIG_FILE_README);
    *iline = i;
  }

  if ( IGNOREFILE(INPUTS.HOSTLIB_FILE) ) { return; }

  readme_docana_comment(&i, "HOSTLIB inputs");

  readme_docana_load_list(&i, pad, &README_KEYS_HOSTLIB);


  // - - - - -
  *iline = i;
  return;

} // end readme_docana_hostlib


void readme_docana_modelPar(int *iline, char *pad) {

  // load population params for stretch, color, alpha, beta ...

  int i = *iline;
  char *cptr, *s, noComment[]="" ;
  int nval1=1, nval2=2, lenkey=24;
  double *dptr, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Population and rate-model params");

  s = INPUTS.GENPDF.MAP_FILE;
  if ( !IGNOREFILE(s) ) {
    ENVrestore(s,ORIG_FILE_README);
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sGENPDF_FILE:  %s ", pad, ORIG_FILE_README);
  }

  if ( INDEX_GENMODEL == MODEL_SALT2 ) {
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_SALT2x1);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_SALT2c);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_SALT2ALPHA);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_SALT2BETA);
  }
  else {
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_DM15);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_DELTA);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_STRETCH);
  }


  int GENAV_WV07 = INPUTS.WV07_GENAV_FLAG ;
  if ( GENAV_WV07 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sGENAV_WV07:  %d ", pad, GENAV_WV07);
  }

  readme_docana_load_expHalfGauss(&i, pad, &INPUTS.GENPROFILE_AV);
  readme_docana_load_expHalfGauss(&i, pad, &INPUTS.GENPROFILE_EBV_HOST);
  readme_docana_load_asymGauss   (&i, pad, &INPUTS.GENGAUSS_RV );

  dptr = INPUTS.BIASCOR_SALT2GAMMA_GRID ;
  VERSION_INFO_load(&i, pad, "BIASCOR_SALT2GAMMA_GRID:", noComment,
		    lenkey, false, nval2, dptr, -1.0,1.0, 9.0); 

  *iline = i;
  return;

} // end readme_docana_modelPar


void readme_docana_rate(int *iline, char *pad) {

  int i = *iline;
  int nval1=1, nval2=2, lenkey=24, m, NMODEL ;
  char *cptr, *s;
  double *dptr, dval;
  // ----------- BEGIN ------------

  readme_docana_load_list(&i, pad, &README_KEYS_RATEMODEL);

  if ( INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT.ORDER > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    s = INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT.STRING;
    sprintf(cptr,"%s%-*s %s # %s",
	    pad, lenkey, "DNDZ_ZPOLY_REWGT:", s, 
	    "dNdz *= polyFun(z)" );
  }

  dptr = &INPUTS.RATEPAR.DNDZ_ZEXP_REWGT ;
  VERSION_INFO_load(&i, pad, "DNDZ_ZEXP_REWGT:", "dN/dz *= z^REWGT",
		    lenkey, false, nval1, dptr, -10.0,10.0, 0.0); 

  dptr = &INPUTS.RATEPAR.DNDZ_ALLSCALE ;
  VERSION_INFO_load(&i, pad, "DNDZ_ALLSCALE:", "dN/dz *= ALLSCALE",
		    lenkey, false, nval1, dptr, 0.0,1.0E4, 1.0); 

  // for DNDZ_SCALE, logic in VERSION_INFO_load doesn't work
  // so need to make explicit check here.
  dptr = INPUTS.RATEPAR.DNDZ_SCALE ;
  if ( dptr[0] != 1.0 || dptr[1] != 1.0 ) {
    VERSION_INFO_load(&i, pad, "DNDZ_SCALE:", 
		      "dN/dz(SNIa,NON1A) *= SCALE[0,1]",
		      lenkey, false, nval2, dptr, 0.0,1.0E4, -9.0 );
  }

  *iline = i;
  return;

} // end readme_docana_rate


void readme_docana_cutwin(int *iline, char *pad) {

  int i = *iline;
  int icut, NCUTWIN = README_KEYS_CUTWIN.NKEY ;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval;
  char *KEY, *ARG, KEY_LAST[100]="" ;

  // ----------- BEGIN ------------

  if ( NCUTWIN == 0 ) { return; }

  readme_docana_comment(&i,"CUTWIN inputs");

  dval = (double)INPUTS.APPLY_CUTWIN_OPT ;
  VERSION_INFO_load(&i, pad, "APPLY_CUTWIN_OPT:", noComment,
                      lenkey, true, nval1, &dval, 0.0, 100, 0.0 );
  
  readme_docana_load_list(&i, pad, &README_KEYS_CUTWIN);

  *iline = i;
  return;

} // end readme_docana_cutwin


void readme_docana_redshift(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval;

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Redshift inputs");

  dptr = INPUTS.GENRANGE_REDSHIFT;
  VERSION_INFO_load(&i, pad, "GENRANGE_REDSHIFT:", noComment, 
		    lenkey, false, nval2, dptr, 0.0,10.0, -1.0); 

  dval = (double)INPUTS.GENSIGMA_REDSHIFT;
  VERSION_INFO_load(&i, pad, "GENSIGMA_REDSHIFT:", noComment, 
		    lenkey, false, nval1, &dval, 0.0,10.0, -1.0); 

  dval = (double)INPUTS.GENSIGMA_VPEC ;
  VERSION_INFO_load(&i, pad, "GENSIGMA_VPEC:", "true vpec scatter (km/sec)",
		    lenkey, false, nval1, &dval, 0.0,9000.0, -1.0); 

  dval = (double)INPUTS.VPEC_ERR ;
  VERSION_INFO_load(&i, pad, "VPEC_ERR:", 
		    "vpec error after correction (km/sec)",
		    lenkey, false, nval1, &dval, 0.0,9000.0, -1.0); 

  dval = (double)INPUTS.VEL_CMBAPEX ;
  VERSION_INFO_load(&i, pad, "VEL_CMBAPEX:", "km/sec",
		    lenkey, false, nval1, &dval, 0.0,400.0, CMBapex_v); 
  

  *iline = i;
  return;

} // end readme_docana_redshift

void readme_docana_epoch(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Epoch & sky coverage inputs");

  dptr = INPUTS.GENRANGE_MJD;
  VERSION_INFO_load(&i, pad, "GENRANGE_MJD:", noComment, 
		    lenkey, false, nval2, dptr, 21000.0,79000.0, -1.0); 

  dptr = INPUTS.GENRANGE_PEAKMJD;
  VERSION_INFO_load(&i, pad, "GENRANGE_PEAKMJD:", noComment, 
		    lenkey, false, nval2, dptr, 1.0E3,1.0E5, -1.0); 

  dval = (double)INPUTS.GENSIGMA_PEAKMJD;
  VERSION_INFO_load(&i, pad, "GENSIGMA_PEAKMJD:", noComment, 
		    lenkey, false, nval1, &dval, 0.0,10.0, 0.0); 

  dval_list[0] = (double)INPUTS.GENRANGE_TREST[0];
  dval_list[1] = (double)INPUTS.GENRANGE_TREST[1];
  VERSION_INFO_load(&i, pad, "GENRANGE_TREST:", noComment, 
		    lenkey, false, nval2, dval_list, -1.0E3,1.0E4, 0.111); 

  dptr = INPUTS.GENRANGE_RA ;
  VERSION_INFO_load(&i, pad, "GENRANGE_RA:", noComment, 
		    lenkey, false, nval2, dptr, -359.0, 360.0, -999.0); 

  dptr = INPUTS.GENRANGE_DEC ;
  VERSION_INFO_load(&i, pad, "GENRANGE_DEC:", noComment, 
		    lenkey, false, nval2, dptr, -359.0, 360.0, -999.0); 

  dval = (double)INPUTS.SOLID_ANGLE ;
  VERSION_INFO_load(&i, pad, "SOLID_ANGLE:", noComment, 
		    lenkey, false, nval1, &dval, 0.0,20.0, 0.0); 

  *iline = i;
  return;

} // end readme_docana_epoch

void readme_docana_misc(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="", *ptrFile ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Misc inputs");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s %s", pad, lenkey, "GENSOURCE:", INPUTS.GENSOURCE);


  ptrFile = PATH_USER_INPUT;
  if ( !IGNOREFILE(ptrFile) ) {
    ENVrestore(ptrFile,ORIG_FILE_README);
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s%-*s %s", pad, lenkey, "PATH_USER_INPUT:", 
	    ORIG_FILE_README);
  }


  ptrFile = INPUTS.PATH_SNDATA_SIM;
  if ( !IGNOREFILE(ptrFile) ) {
    ENVrestore(ptrFile,ORIG_FILE_README);
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s%-*s %s", pad, lenkey, "PATH_SNDATA_SIM:", 
	    ORIG_FILE_README);
  }

  dval = (double)INPUTS.ISEED_ORIG ;
  VERSION_INFO_load(&i, pad, "RANSEED:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 

  dval = (double)INPUTS.DEBUG_FLAG ;
  VERSION_INFO_load(&i, pad, "DEBUG_FLAG:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 

  dptr = INPUTS.GENRANGE_PEAKMAG;
  VERSION_INFO_load(&i, pad, "GENRANGE_PEAKMAG:", noComment, 
		    lenkey, false, nval2,dptr, 0.0,40.0, -999.0); 

  readme_docana_load_list(&i, pad, &README_KEYS_CID);

  readme_docana_load_list(&i, pad, &README_KEYS_RANSYSTPAR);

  readme_docana_load_list(&i, pad, &README_KEYS_ZVARIATION);

  readme_docana_load_list(&i, pad, &README_KEYS_GRIDGEN);

  readme_docana_load_list(&i, pad, &README_KEYS_SIMGEN_DUMP);

  *iline = i;
  return;

} // end readme_docana_misc


void readme_docana_mwebv(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Galactic extinction");

  readme_docana_load_list(&i, pad, &README_KEYS_MWEBV);

  *iline = i;
  return;

} // end readme_docana_mwebv

void readme_docana_searcheff(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];
  // ----------- BEGIN ------------

  readme_docana_comment(&i, "SEARCHEFF/detections");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s  %s", 
	  pad, lenkey, "SEARCHEFF_PIPELINE_LOGIC_FILE:", ORIG_FILE_README);

  dval = (double)INPUTS.NEWMJD_DIF ;
  VERSION_INFO_load(&i, pad, "NEWMJD_DIF:", 
		    "day-sep if > 1 detections required",
		    lenkey, false, nval1, &dval, 0.0,2000.0, -1.0); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s  %s", 
	  pad, lenkey, "SEARCHEFF_PIPELINE_EFF_FILE:", ORIG_FILE_README);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS_SEARCHEFF.USER_SPEC_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s %s", 
	  pad, lenkey, "SEARCHEFF_SPEC_FILE:", ORIG_FILE_README);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS_SEARCHEFF.USER_zHOST_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s %s", 
	  pad, lenkey, "SEARCHEFF_zHOST_FILE:", ORIG_FILE_README);

  dval = (double)INPUTS.APPLY_SEARCHEFF_OPT ;
  VERSION_INFO_load(&i, pad, "APPLY_SEARCHEFF_OPT:", 
		    "+= 1,2,4 => pipe,spec,zhost", 
		    lenkey, true, nval1, &dval, 0.0,2000.0, -1.0); 

  *iline = i;
  return;

} // end readme_docana_searcheff


void readme_docana_output(int *iline, char *pad) {

  // Store GENVERSTION, FORMAT_MASK, etc ...
  // Aug 2022: add CIDRAN_SKIPLIST

  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Output data");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s %s", pad, lenkey, "GENVERSION:", INPUTS.GENVERSION);

  dval = (double)INPUTS.NGENTOT_LC ;
  VERSION_INFO_load(&i, pad, "NGENTOT_LC:", noComment,
		    lenkey, true, nval1, &dval, 1.0,1.0E8, 0.0); 

  dval = (double)INPUTS.NGEN_SEASON ;
  VERSION_INFO_load(&i, pad, "NGEN_SEASON:", noComment,
		    lenkey, true, nval1, &dval, 0.0,1.0E8, 0.0); 

  dval = (double)INPUTS.FORMAT_MASK ;
  VERSION_INFO_load(&i, pad, "FORMAT_MASK:", 
		    " += 2,32,16 -> TEXT, FITS, randomCID", 
		    lenkey, true, nval1, &dval, 0.0,2000.0, -1.0); 

  // - - - -- types - - - - -
  dval = (double)GENLC.SIM_GENTYPE ;
  VERSION_INFO_load(&i, pad, "GENTYPE:", "true type", 
		    lenkey, true, nval1, &dval, 0.0,2000.0, -1.0); 

  dval_list[0] = (double)INPUTS.SNTYPE_Ia_SPEC;
  dval_list[1] = (double)INPUTS.SNTYPE_Ia_PHOT;
  VERSION_INFO_load(&i, pad, "SNTYPE:", "spec Type, photID type", 
		    lenkey, true, nval2, dval_list, 0.0,2000.0, -1.0); 

  /* xxx mark delete Nov 12 2023 (check below under Misc) xxxx
  dval = (double)INPUTS.CIDOFF ;
  VERSION_INFO_load(&i, pad, "CIDOFF:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 
  xxx */

  dval = (double)INPUTS.CIDRAN_MIN ;
  VERSION_INFO_load(&i, pad, "CIDRAN_MIN:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 

  dval = (double)INPUTS.CIDRAN_MAX ;
  VERSION_INFO_load(&i, pad, "CIDRAN_MAX:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 
 
  if ( INPUTS.NCIDRAN_SKIPLIST > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sCIDRAN_SKIPLIST:  %s", pad, INPUTS.CIDRAN_SKIPLIST_STRING);
  }

  *iline = i;
  return;

} // end readme_docana_output

void readme_docana_template(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------
  // copy this and change 'template' to something else
  *iline = i;
  return;

} // end readme_docana_template

void  readme_docana_load_list(int *iline, char *pad,
                              README_KEYPLUSARGS_DEF *README_KEYS) {

  // Created Dec 2021
  // Load list of NKEY keys. If KEY is unique, write it as
  //   KEY ARG
  // If there are duplicate keys, write them as
  //   KEY
  //   - ARG0
  //   - ARG1
  //   etc ...
  // Each KEY_LIST item is assumed to include a colon.
  // Beware that duplicates are assumed to be sequential in the list.
  // Non-sequential duplicates will both be written as "KEY ARG".
  // 

  int i = *iline;
  int NKEY = README_KEYS->NKEY;
  int k, lenkey=24 ;
  bool NEW, UNIQUE;
  char *KEY, *ARG, *KEY_NEXT, KEY_LAST[100] = "", *cptr ;

  // ------------- BEGIN -------------

  for(k=0; k < NKEY; k++ ) {
    KEY = README_KEYS->KEY_LIST[k]; 
    ARG = README_KEYS->ARG_LIST[k];
    NEW = (strcmp(KEY,KEY_LAST) != 0) ; // it's a new key

    if ( k < NKEY-1 ) {
      KEY_NEXT = README_KEYS->KEY_LIST[k+1];
      UNIQUE = NEW && (strcmp(KEY,KEY_NEXT) != 0 ); 
    }
    else
      { UNIQUE = NEW; } 

    if ( NEW && UNIQUE ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s%-*s %s", pad, lenkey, KEY, ARG); // write "KEY ARG"
    }
    else if ( NEW ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s%-*s", pad, lenkey, KEY); // write "KEY"
    }

    if ( !UNIQUE ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s- %s", pad, ARG); // write "- ARG"
    }

    sprintf(KEY_LAST,"%s", KEY);
  }

  *iline = i;
  return;

} // end readme_docana_load_list


void readme_docana_load_asymGauss(int *iline, char *pad,
				 GENGAUSS_ASYM_DEF *GENGAUSS) {

  // Utility to load asymGauss params for DOCANA
  // Nov 2023: noprint range for GENSIGMA -> 1E9 (was 1E4)

  int i = *iline ;
  int nval1=1, nval2=2, lenkey=24 ;
  bool USE      = GENGAUSS->USE ;
  char *VARNAME = GENGAUSS->NAME;
  char noComment[] = "" ;

  char keyName[60];
  double *dptr, dval ;
  // ------------ BEGIN -------------
 
  if ( !USE ) { return; }

  // start with required elements

  sprintf(keyName,"GENPEAK_%s:", VARNAME);
  dptr = &GENGAUSS->PEAK ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval1, dptr, -1.E4,1.E4, -9.0); // val,min,max noprint

  sprintf(keyName,"GENSIGMA_%s:", VARNAME);
  dptr = GENGAUSS->SIGMA ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval2, dptr, -1.E9,1.E9, -9.0); // val,min,max noprint

  sprintf(keyName,"GENRANGE_%s:", VARNAME);
  dptr = GENGAUSS->RANGE ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval2, dptr, -1.E4,1.E4, -9.0); // val,min,max noprint

  // optional elements

  sprintf(keyName,"GENGRID_%s:", VARNAME);
  dval = (double)GENGAUSS->NGRID ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, true,
		    nval1, &dval, 0,5, 0.0); // val,min,max noprint

  sprintf(keyName,"GENPEAK2_%s:", VARNAME);
  dptr = &GENGAUSS->PEAK2 ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval1, dptr, -1.E4,1.E4, 0.0); // val,min,max noprint

  sprintf(keyName,"GENPROB2_%s:", VARNAME);
  dptr = &GENGAUSS->PROB2 ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval1, dptr, -1.E4,1.E4, 0.0); // val,min,max noprint

  sprintf(keyName,"GENSIGMA2_%s:", VARNAME);
  dptr = GENGAUSS->SIGMA2 ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval2, dptr, -1.E4,1.E4, 0.0); // val,min,max noprint

  if(  GENGAUSS->PEAKRANGE[1] >  GENGAUSS->PEAKRANGE[0] ) {
    sprintf(keyName,"PEAKRANGE_%s:", VARNAME);
    dptr = GENGAUSS->PEAKRANGE ;
    VERSION_INFO_load(&i, pad, keyName,
		      noComment, lenkey, false,
		      nval2, dptr, -1.E4,1.E4, 0.0); // val,min,max noprint
  }

  *iline = i;
  return;
} // end readme_docana_load_asymGauss


void readme_docana_load_expHalfGauss(int *iline, char *pad,
                                     GEN_EXP_HALFGAUSS_DEF *EXP_HALFGAUASS) {

  // load params for Exp + halfGauss function
  int i = *iline ;
  int nval1=1, nval2=2, lenkey=24 ;
  bool USE      = EXP_HALFGAUASS->USE;
  char *VARNAME = EXP_HALFGAUASS->NAME;
  char noComment[] = "" ;

  char keyName[60];
  double *dptr, dval;

  // ------------- BEGIN -----------
	 
  if ( !USE ) { return; }

  sprintf(keyName,"GENTAU_%s:", VARNAME);
  dptr = &EXP_HALFGAUASS->EXP_TAU ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval1, dptr, 0.0, 1.0E5, -9.0); // val,min,max noprint

  sprintf(keyName,"GENRANGE_%s:", VARNAME);
  dptr = EXP_HALFGAUASS->RANGE ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval2, dptr, 0.0,100.0, -9.0); // val,min,max noprint
  
  // optonal half-Gauss params
  if ( EXP_HALFGAUASS->RATIO == 0.0 ) { *iline=i; return; }

  sprintf(keyName,"GENGAUPEAK_%s:", VARNAME);
  dptr = &EXP_HALFGAUASS->PEAK ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    "peak of half-Gauss", lenkey, false, 
		    nval1, dptr, 0.0,100.0, -9.0); // val,min,max noprint

  sprintf(keyName,"GENSIGMA_%s:", VARNAME);
  dptr = &EXP_HALFGAUASS->SIGMA ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    "sigma of half-Gauss", lenkey, false, 
		    nval1, dptr, 0.0,100.0, -9.0); // val,min,max noprint

  sprintf(keyName,"GENRATIO_%s:", VARNAME);
  dptr = &EXP_HALFGAUASS->RATIO ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    "Gauss(0)/Expon(0)", lenkey, false, 
		    nval1, dptr, 0.0,100.0, -9.0); // val,min,max noprint

  *iline = i;
} // end readme_docana_load_expHalfGauss


// *********************************************
void VERSION_INFO_load(int *iline, char *pad, char *keyName,  char *comment,
		       int lenkey, bool isint, int nval,
		       double *val_list, double valmin, double valmax, 
		       double val_noprint ) {

  // created Dec 22 2021
  // Load info to VERSION_INFO.README_DOC if:
  //   valmin <= val[0] <= valmax && val[0] != val_noprint
  // "nval" is the number of values following keyname.
  // include *comment if not null.

  int  ival, i = *iline ;
  double val0   = val_list[0];
  double val ;
  bool noprint  = fabs(val0 - val_noprint) < 1.0E-8;
  bool passcut  = (val0 >= valmin   &&   val0 <= valmax );
  bool do_print = passcut && !noprint ;
  char *cptr, cval[40] ;
  char fnam[] = "VERSION_INFO_load" ;

  // ------------ BEGIN ------------

  if ( !do_print ) { return; }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s ",  pad, lenkey, keyName );

  for(ival=0; ival < nval; ival++ ) {
    val = val_list[ival];
    if ( isint || val == 0.0 ) 
      { sprintf(cval,"%d  ", (int)val ); }
    else if ( fabs(val) > 0.01 ) 
      { sprintf(cval,"%.4f  ", val ); }
    else
      { sprintf(cval,"%.4le  ", val ); }

    strcat(cptr,cval);
  }
    
  if ( comment[0] != 0 ) { strcat(cptr,"# "); strcat(cptr,comment); }

  *iline = i;
  return ;
} // end VERSION_INFO_load


void README_KEYPLUSARGS_init(README_KEYPLUSARGS_DEF *README_KEYS) {
  README_KEYS->NKEY = 0;
  README_KEYS->MALLOC1 = false;
}

// =============================================================
void README_KEYPLUSARGS_load(int MXKEY, int NWD, char **WORDS, int keySource,
			     README_KEYPLUSARGS_DEF *README_KEYS,
			     char *callFun ) {

  // Store NWD WORDS in README_KEYS strut.
  // Inputs:
  //   MXKEY : max number of keys to store
  //   NWD   : number of args after key
  //   WORDS : list of words;  KEY=WORDS[0], ARGS=WORDS[1:N]
  //   callFun: calling function, for error message
  //
  // Output:
  //   README_KEYS : structure to load
  //
  // Jan 11 2022: pass keySource arg, and check override only for
  //              command-line override.
  //
  // Mar 26 2023: allow special override for GENMODEL key to appear
  //              in two sim-input files.
  //
  bool MALLOC1 = README_KEYS->MALLOC1 ;
  int  NKEY    = README_KEYS->NKEY;  
  int  MEMC1   = MXKEY * sizeof(char*);
  int  MEMC_KEY, MEMC_ARG;
  int iwd, lenkey, k, LENKEY=0, LENARG=0 ;
  char *KEY, *ARG, *ARG_TMP, *KEY_TMP ;
  char BLANK[] = " ";
  char fnam[] = "README_KEYPLUSARGS_load" ; 

  // ------------ BEGIN ----------

  if ( NKEY == 0 && !MALLOC1 ) {
    // allocate pointer for all MXKEY possible keys
    README_KEYS->KEY_LIST = (char**) malloc(MEMC1);
    README_KEYS->ARG_LIST = (char**) malloc(MEMC1);
    README_KEYS->MALLOC1 = true;
  }

  // compute lenth of string needed to store KEY and  WORDS
  LENKEY = strlen(WORDS[0]) + 10;

  LENARG = NWD + 10 ;
  for(iwd=1; iwd<=NWD; iwd++ )  { LENARG += strlen(WORDS[iwd]) ; } 

  MEMC_KEY = ( LENKEY ) * sizeof(char);
  MEMC_ARG = ( LENARG ) * sizeof(char);

  KEY = (char*) malloc(MEMC_KEY); ;  KEY[0]=0;
  ARG = (char*) malloc(MEMC_ARG); ;  ARG[0]=0;

  // Store key in local variable
  // if there is no colon after key, add it to word for command-line
  // overrides without colon
  sprintf(KEY, "%s", WORDS[0]);
  if ( strchr(KEY,':') == NULL ) { strcat(KEY,COLON); }
  
  // load args into local string
  for(iwd=1; iwd<=NWD; iwd++ )  { strcat(ARG,WORDS[iwd]); strcat(ARG,BLANK);  }

  // printf(" xxx %s: check '%s' = '%s' \n", fnam, KEY, ARG); fflush(stdout);

  // For command-line override, check previous keys to override.
  // Mar 2023: special override allows GENMODEL key to appear in
  //    main sim-input file, and to be overwritten in an INCLUDE file.
  bool OVERRIDE_COMMAND = ( keySource == KEYSOURCE_ARG );
  bool OVERRIDE_SPECIAL = ( strcmp(KEY,"GENMODEL:") == 0 ) ;

  if ( OVERRIDE_COMMAND || OVERRIDE_SPECIAL ) {
    for(k=0; k < NKEY; k++ ) {
      KEY_TMP = README_KEYS->KEY_LIST[k];
      ARG_TMP = README_KEYS->ARG_LIST[k];

      if ( strcmp(KEY,KEY_TMP) == 0 ) {
	// check if the arg override needs more memory than original key
	if ( LENARG > strlen(ARG_TMP) ) {
	  README_KEYS->ARG_LIST[k] = (char *)realloc(README_KEYS->ARG_LIST[k],MEMC_ARG);
	  ARG_TMP = README_KEYS->ARG_LIST[k];
	}
	sprintf(ARG_TMP,"%s", ARG);
	free(KEY); free(ARG);
	return ;
      }
    } // end k
  }

  // ---- store new KEY and ARG

  // allocate mem for this key and arg
  README_KEYS->KEY_LIST[NKEY] = (char*) malloc(MEMC_KEY);
  README_KEYS->ARG_LIST[NKEY] = (char*) malloc(MEMC_ARG);

  sprintf(README_KEYS->KEY_LIST[NKEY], "%s", KEY);
  sprintf(README_KEYS->ARG_LIST[NKEY], "%s", ARG);

  // increment number of stored keys in this structure.
  README_KEYS->NKEY++ ;

  free(KEY); free(ARG);

  return;

} // README_KEYPLUSARGS_load
