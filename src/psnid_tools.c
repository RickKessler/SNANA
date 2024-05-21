/* =======================================
 Created July 24, 2012 by R.Kessler

 Miscellaneous functions used by psnid_[METHOD] routines.
 The functions here are not specific to any particular 
 method. Functions include reading templates, dump-functions
 and numerical recipes.


    History
   ~~~~~~~~~~

 Oct 25 2012 MS - fixed typo in SNIA -> NONIA index

 Oct 11 2012 RK - parse new PSNID_INPUTS.OPT_ZPRIOR and .CUTWIN_ZERR

 Nov 29 2012 MS - migrating PSNID variables to namelist

 Dec 02 2012 RK - fix compile bug found by Ishida

 Feb 08, 2013 RK - add psnid_store_data() function to store data
                   in DATA_PSNID_DOFIT structure.

 Feb 27, 2013 RK - improve NONIa summary in
                   psnid_dump_nonIa_templates(void) 

 July 28, 2013 DB - Added PBayes & Fitprob namelist cuts

 Sep 7, 2013 RK - check template_nonIa_list for optional subset of
                  templates to use.

 Mar 14 2016:  move a few utilities into sntools_gridread.c ,
    psnid_init_nonIa_wgts      -> renorm_wgts_SNGRID
    psnid_dump_nonIa_templates -> dump_SNGRID

 Aug 27 2017: 
   + few updates to allow using only NON1A templates without SNIA templates.

 Oct 23 2020: call init_HzFUN_INFO

========================================= */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_sf_gamma.h>

#include "sntools.h"            // general snana stuff
#include "sntools_cosmology.h"  // cosmology functions (10.2020)
#include "MWgaldust.h"    // GALextinct is here

#define MODELGRID_READ // use only the read utilities in sngridtools.c

#include "fitsio.h"
#include "sntools_modelgrid.h"

#include "psnid_tools.h"  // psnid tools (after including sngrindtools)


// ====================================================================
void PSNID_USER_INPUT(int NVAR, double *input_array, char *input_string ) {

  // Receive input values, then strip them off and load
  // PSNID_INPUTS structure in psnid_tools.h.

  int ivar, i, NWD, iwd, OPT, iter, VBOSE ;
  double dval, cosPar[10] ;
  char *ptrtok, cwd[60][200] ;
  char fnam[] = "PSNID_USER_INPUT" ;

  int LDMP = 1 ;

  // ------------- BEGIN -----------

  ivar = -1 ;

  // load &SNLCINP doubles
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.H0 = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.OMAT = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.OLAM = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.W0 = dval ;

  // Oct 2020: create HzFUN_INFO struct
  cosPar[ICOSPAR_HzFUN_H0] = PSNID_INPUTS.H0 ;
  cosPar[ICOSPAR_HzFUN_OM] = PSNID_INPUTS.OMAT ;
  cosPar[ICOSPAR_HzFUN_OL] = PSNID_INPUTS.OLAM ;
  cosPar[ICOSPAR_HzFUN_w0] = PSNID_INPUTS.W0 ;
  cosPar[ICOSPAR_HzFUN_wa] = 0.0 ;
  VBOSE = 1;
  init_HzFUN_INFO(VBOSE, cosPar, "", &PSNID_INPUTS.HzFUN_INFO);
  PSNID_INPUTS.ANISOTROPY_INFO.USE_FLAG = false; // Feb 2023

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.DEBUG_FLAG = (int)dval ;

  // load &PSNIDINP doubles ...
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.AV_TAU = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.AV_SMEAR = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.AV_PRIOR_STR = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.WGT_ZPRIOR = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.OPT_ZPRIOR = (int)dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.CUTWIN_ZERR[0] = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.CUTWIN_ZERR[1] = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.MCMC_NSTEP = (int)dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.COLOR_MIN = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.COLOR_MAX = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.NCOLOR = (int)dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.DMU_MIN = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.DMU_MAX = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.DMU_NON1A_MIN = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.DMU_NON1A_MAX = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.NDMU = (int)dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.TEMPLERR_SCALE_SNIA = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.TEMPLERR_SCALE_NONIA = dval ;

  // - - - - -
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.PBAYES_IA_CUT = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.PBAYES_IBC_CUT = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.PBAYES_II_CUT = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.PBAYES_PEC1A_CUT = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.PBAYES_MODEL1_CUT = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.PBAYES_MODEL2_CUT = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.PBAYES_MODEL3_CUT = dval ;
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.PBAYES_MODEL4_CUT = dval ;
  
  // - - - - -

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.FITPROB_IA_CUT = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.FITPROB_IBC_CUT = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.FITPROB_II_CUT = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.FITPROB_PEC1A_CUT = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.FITPROB_MODEL1_CUT = dval ;
    ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.FITPROB_MODEL2_CUT = dval ;
    ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.FITPROB_MODEL3_CUT = dval ;
    ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.FITPROB_MODEL4_CUT = dval ;

  // - - - - -
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.CHISQMIN_OUTLIER = dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.NREJECT_OUTLIER = (int)dval ;

  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.OPT_SIMCHEAT = (int)dval ;

  for(iter=0; iter < MXITER_PSNID; iter++ ) {
    ivar++ ; dval = input_array[ivar];
    PSNID_INPUTS.TMAX_START[iter] = dval ;
    ivar++ ; dval = input_array[ivar];
    PSNID_INPUTS.TMAX_STOP[iter]  = dval ;
    ivar++ ; dval = input_array[ivar];
    PSNID_INPUTS.TMAX_STEP[iter]  = dval ;
  }
    
  // --------
  ivar++ ; dval = input_array[ivar];
  PSNID_INPUTS.OPT_RATEPRIOR = (int)dval ;
  for(i=0; i<3; i++ ) { 
    ivar++ ; dval = input_array[ivar] ;
    PSNID_INPUTS.ZRATEPRIOR_SNIA[i] = dval ; 
  }
  for(i=0; i<3; i++ ) { 
    ivar++ ; dval = input_array[ivar] ;
    PSNID_INPUTS.ZRATEPRIOR_NONIA[i] = dval ; 
  }

  // -----------------------------------------------
  // break the input string into separate words
  //  printf(" xxx input_string = '%s' \n", input_string);


  ptrtok = strtok(input_string," ");
  NWD = 0;
  while ( ptrtok != NULL  ) {
    sprintf(cwd[NWD], "%s", ptrtok );
    NWD++ ;
    ptrtok = strtok(NULL, " ");
  }
  // load strings ...


  for ( iwd=0 ; iwd < NWD ; iwd++ ) {

    if ( strcmp(cwd[iwd],"NMLFILE:") == 0 ) 
      { sprintf(PSNID_INPUTS.NMLFILE, "%s", cwd[iwd+1]) ; }

    if ( strcmp(cwd[iwd],"MODELNAME_MAGERR:") == 0 ) 
      { sprintf(PSNID_INPUTS.MODELNAME_MAGERR, "%s", cwd[iwd+1]) ; }

    if ( strcmp(cwd[iwd],"FILTLIST_PEAKMAG_STORE:") == 0 ) 
      { sprintf(PSNID_INPUTS.CFILTLIST_PEAKMAG_STORE, "%s", cwd[iwd+1]) ; }
             
  }  // end of iwd
 

  // error checking
  if ( ivar+1 != NVAR ) {
    sprintf(c1err,"Expected NVAR = %d from input_array", NVAR);
    sprintf(c2err,"but unpacked %d values.", ivar+1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // optional 2nd dump since it's already dumped by the nml dump
  // from psnid.cra

  if( LDMP ) {

    printf("\n %s DUMP  PSNID_INPUTS struct (C code): \n", fnam ) ;

    printf("\t Input nmlFile    = '%s' \n", PSNID_INPUTS.NMLFILE ) ;
  
    printf("\t Input AV[TAU,SMEAR,PRIOR_STR] = %f , %f, %f \n",
	   PSNID_INPUTS.AV_TAU, PSNID_INPUTS.AV_SMEAR, 
	   PSNID_INPUTS.AV_PRIOR_STR );

    printf("\t Input OPT_ZPRIOR = %d \n", PSNID_INPUTS.OPT_ZPRIOR);

    printf("\t Input WGT_ZPRIOR = %f \n", PSNID_INPUTS.WGT_ZPRIOR);

    printf("\t Input CUTWIN_ZERR = %f,%f \n", 
	   PSNID_INPUTS.CUTWIN_ZERR[0], PSNID_INPUTS.CUTWIN_ZERR[1] );

    // ----
    OPT = PSNID_INPUTS.OPT_RATEPRIOR ;
    printf("\t Input OPT_RATEPRIOR = %d \n", OPT);
    if ( OPT > 0 ) {
      printf("\t    dN/dz(SNIA) = %.2le * (1+z)^%4.2f,  constant for z>%.2f\n"
	     ,PSNID_INPUTS.ZRATEPRIOR_SNIA[0] 
	     ,PSNID_INPUTS.ZRATEPRIOR_SNIA[1]
	     ,PSNID_INPUTS.ZRATEPRIOR_SNIA[2] );
      printf("\t    dN/dz(SNCC) = %.2le * (1+z)^%4.2f,  constant for z>%.2f\n"
	     ,PSNID_INPUTS.ZRATEPRIOR_NONIA[0] 
	     ,PSNID_INPUTS.ZRATEPRIOR_NONIA[1]
	     ,PSNID_INPUTS.ZRATEPRIOR_NONIA[2] );
    }


    // ---

    printf("\t Input MCMC_NSTEP = %d \n", PSNID_INPUTS.MCMC_NSTEP);

    printf("\t Input COLOR_MIN, COLOR_MAX, NCOLOR = %f %f %d \n"
	   ,PSNID_INPUTS.COLOR_MIN
	   ,PSNID_INPUTS.COLOR_MAX
	   ,PSNID_INPUTS.NCOLOR ) ;

    printf("\t Input DMU_MIN, DMU_MAX, NDMU = %f %f %d \n",
    	   PSNID_INPUTS.DMU_MIN, PSNID_INPUTS.DMU_MAX, PSNID_INPUTS.NDMU);

    printf("\t Input MODELNAME_MAGERR = '%s' \n",
	   PSNID_INPUTS.MODELNAME_MAGERR );

    printf("\t Input TEMPLERR_SCALE_SNIA = %f \n",
	   PSNID_INPUTS.TEMPLERR_SCALE_SNIA);

    printf("\t Input TEMPLERR_SCALE_NONIA = %f \n",
	   PSNID_INPUTS.TEMPLERR_SCALE_NONIA);

    // - - - - 
    printf("\t Input PBAYES_[IA,IBC,II,PEC1A]_CUT = "
	   "%.3f, %.3f, %.3f, %.3f \n"
	   ,PSNID_INPUTS.PBAYES_IA_CUT
           ,PSNID_INPUTS.PBAYES_IBC_CUT
           ,PSNID_INPUTS.PBAYES_II_CUT
           ,PSNID_INPUTS.PBAYES_PEC1A_CUT );

    printf("\t Input PBAYES_[MODEL(1,2,3,4)]_CUT = "
	   "%.3f, %.3f, %.3f, %.3f \n"
	   ,PSNID_INPUTS.PBAYES_MODEL1_CUT
           ,PSNID_INPUTS.PBAYES_MODEL2_CUT
           ,PSNID_INPUTS.PBAYES_MODEL2_CUT
           ,PSNID_INPUTS.PBAYES_MODEL4_CUT );

    // - - - - 
    printf("\t Input FITPROB_[IA,IBC,II,PEC1A]_CUT = "
	   "%.3f, %.3f, %.3f, %.3f \n"
           ,PSNID_INPUTS.FITPROB_IA_CUT
           ,PSNID_INPUTS.FITPROB_IBC_CUT
           ,PSNID_INPUTS.FITPROB_II_CUT
           ,PSNID_INPUTS.FITPROB_PEC1A_CUT );

        printf("\t Input FITPROB_[MODEL(1,2,3,4)]_CUT = "
	   "%.3f, %.3f, %.3f, %.3f \n" 
           ,PSNID_INPUTS.FITPROB_MODEL1_CUT
           ,PSNID_INPUTS.FITPROB_MODEL2_CUT
           ,PSNID_INPUTS.FITPROB_MODEL3_CUT
           ,PSNID_INPUTS.FITPROB_MODEL4_CUT );

    // - - - - 

    printf("\t Input CHISQMIN_OUTLIER, NREJECT_OUTLIER = %f %d \n",
           PSNID_INPUTS.CHISQMIN_OUTLIER,
           PSNID_INPUTS.NREJECT_OUTLIER);

    printf("\n" ) ;
    fflush(stdout);

  }

} // end of PSNID_USER_INPUT

void psnid_user_input__( int *NVAR, double *input_array, char *input_string){
  PSNID_USER_INPUT(*NVAR, input_array, input_string ) ;
}
  

// ==================================
void PSNID_INIT_VAR(void) {

  // Misc. inits 
  int ifilt, itype ;

  FLAG_PSNID_INIT_VAR = 1;

  printf("   Init PSNID variables. \n");
  
  PSNID_INPUTS.CFILTLIST[0] = 0;
  PSNID_INPUTS.CFILTLIST_PEAKMAG_STORE[0] = 0;
  PSNID_INPUTS.NFILT  = 0 ;
  PSNID_INPUTS.WRSTAT_TABLE   = 0 ;
  PSNID_INPUTS.LSIM           = 0 ;

  for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
    PSNID_INPUTS.USEFILT[ifilt] = 0 ;
    PSNID_INPUTS.USEFILT_PEAKMAG_STORE[ifilt] = 0 ;
  }

  PSNID_FITRES.NVAR   = 0 ;
  NLCDUMP_PSNID = 0 ;

  for(itype=0 ; itype <= MXTYPEINDX_PSNID; itype++ ) {   
    USEFLAG_TEMPLATES_PSNID[itype] = 0 ;
    TEMPLATES_FILE_PSNID[itype][0]=0;
  }
  
  sprintf(TEMPLATETYPE_PSNID[TYPEINDX_SNIA_PSNID],  "SNIa" );
  sprintf(TEMPLATETYPE_PSNID[TYPEINDX_NONIA_PSNID], "NONIa");


  // define public area for PSNID templates 
  sprintf(PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") ) ;
  sprintf(PATH_TEMPLATES_PSNID, "%s/models/psnid", PATH_SNDATA_ROOT);
  sprintf(PSNID_INPUTS.SNANA_VERSION, "%s", SNANA_VERSION_CURRENT);

  // init data-storage array
  DATA_PSNID_DOFIT.NOBS = 0 ;

} // end of PSNID_INIT_VAR

void psnid_init_var__(void) {
  PSNID_INIT_VAR();
}


// ==================================
void PSNID_INIT_FILTERS(char *filtlist_fit) {

  // Set USEFILT[ifilt_obs] to check later
  // which filters to use in the fit. Also check that
  // filters are valid; abort on error.

  int ifilt, ifilt_obs ;
  char cfilt[2], ctmp[40], *cfiltlist;
  char fnam[] = "PSNID_INIT_FILTERS" ;

  // ------------- BEGIN ----------------

  sprintf(PSNID_INPUTS.CFILTLIST, "%s", filtlist_fit);

  // define absolute filter indices in SNANA.
  set_FILTERSTRING(FILTERSTRING);
  //  sprintf(FILTERSTRING,"%s", FILTERSTRING_DEFAULT );


  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ )  { 

    // init filter use-mask
    PSNID_INPUTS.USEFILT[ifilt]  = 0; 

    // init GRID offsets
    PSNID_INPUTS.USEFILT_GRIDOFF[TYPEINDX_SNIA_PSNID][ifilt]  = 0; 
    PSNID_INPUTS.USEFILT_GRIDOFF[TYPEINDX_NONIA_PSNID][ifilt] = 0; 
  }

  
  PSNID_INPUTS.NFILT = strlen(PSNID_INPUTS.CFILTLIST);

  // idiot check
  if ( PSNID_INPUTS.NFILT <= 0 ) {
    sprintf(c1err,"Invalid NFILT = %d", PSNID_INPUTS.NFILT);
    sprintf(c2err,"Check &PSNIDINP variable FILTLIST_FIT = '%s'", 
	    PSNID_INPUTS.CFILTLIST);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // set USEFILT for each filter so we can quickly check
  // which ones to include in the psnid fit

  for ( ifilt=0; ifilt < PSNID_INPUTS.NFILT; ifilt++ ) {
    sprintf(cfilt,"%c", PSNID_INPUTS.CFILTLIST[ifilt]);
    ifilt_obs = INTFILTER(cfilt); // get absolute integer index
    PSNID_INPUTS.IFILTLIST[ifilt]         = ifilt_obs ;
    PSNID_INPUTS.IFILTLIST_INV[ifilt_obs] = ifilt ;

    if ( ifilt_obs <= 0 ) {
      sprintf(c1err,"Invalid filter = '%s'", cfilt);
      sprintf(c2err,"Check &PSNIDINP variable FILTLIST_FIT = '%s'",
	      PSNID_INPUTS.CFILTLIST);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    PSNID_INPUTS.USEFILT[ifilt_obs]   = 1;

    // check list of filters to store peak mag
    cfiltlist = PSNID_INPUTS.CFILTLIST_PEAKMAG_STORE ;
    if ( strstr(cfiltlist,cfilt) != NULL ) {
      PSNID_INPUTS.USEFILT_PEAKMAG_STORE[ifilt_obs]   = 1 ;
      sprintf(ctmp,"and store peak-mag in output");
    }
    else  {
      PSNID_INPUTS.USEFILT_PEAKMAG_STORE[ifilt_obs]   = 0 ;
      ctmp[0] = 0 ;
    }

    printf("\t PSNID will fit filter %s (%2d)  %s \n", 
	   cfilt, ifilt_obs, ctmp );
    fflush(stdout);
  }


  // now check which filters to store peak-mag

  fflush(stdout);

} // end of  PSNID_INIT_FILTERS

void psnid_init_filters__(char *filtlist_fit) {
  PSNID_INIT_FILTERS(filtlist_fit);
}


// *************************************************
void  PSNID_INIT_XTMW(int OPT_MWCOLORLAW) {

  // Created Apr 2013 by R.Kessler
  // For AV=1, compute MilkyWay Galactic extinction at 
  // the central wavelength of each filter.  
  // Fill PSNID_INPUTS.XTMW_at_AV1[ifilt] 
  //
  // Sep 18 2013: pass OPT_MWCOLORLAW as argument for GALextinct

  double LAMAVG, RV, AV, XTMW, NXTMW ; ;
  int ifilt ;
  char cfilt[2] ;
  //  char fnam[] = "PSNID_INIT_XTMW" ;

  // ----------------- BEGIN ------------

  RV    = RV_MWDUST ;
  AV    = 1.0 ;
  NXTMW = 0;

  for ( ifilt=0; ifilt < PSNID_INPUTS.NFILT; ifilt++ ) {
    PSNID_INPUTS.XTMW_at_AV1[ifilt] = 0.0 ; // init in case no LAMAVG info
    sprintf(cfilt,"%c", PSNID_INPUTS.CFILTLIST[ifilt]);
    
    LAMAVG = (double)SNGRID_PSNID[TYPEINDX_SNIA_PSNID].FILTER_LAMAVG[ifilt];
    if ( LAMAVG <= 0.0 ) { continue; }

    XTMW   = GALextinct ( RV, AV, LAMAVG, OPT_MWCOLORLAW );
    PSNID_INPUTS.XTMW_at_AV1[ifilt] = XTMW ; // store in global
    NXTMW++ ;

    // print values for first 10 filters
    if ( ifilt < 10 ) 
      { printf("\t XTMW(%s, AV=1,RV=%3.1f) = %6.3f \n", cfilt, RV, XTMW); }
  }


  if ( NXTMW == 0 ) {
    printf("\t WARNING: LAMAVG not defined -> no Galactic extinction.\n");
  }

} // end of PSNID_INIT_XTMW

void psnid_init_xtmw__(int *OPT_MWCOLORLAW) {
  PSNID_INIT_XTMW(*OPT_MWCOLORLAW);
}

// ===========================================
void PSNID_READ_TEMPLATES(int TYPEINDX, char *file ) {

  // read templates from FITS-formatted file.
  // Check local dir first; if not there then check 
  // $SNDATA_ROOT/models/psnid.

  FILE *fp;
  int OPT_READ = 1; // 1-> verbose mode
  int gzipFlag ;
  char FILE[200], *ptrFile ;
  char fnam[] = "PSNID_READ_TEMPLATES";

  // ------------- BEGIN --------------

  if ( FLAG_PSNID_INIT_VAR != 1 ) {
    sprintf(c1err,"PSNID_INIT_VAR was not called.");
    sprintf(c2err,"Check functions called.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }

  sprintf(TEMPLATES_FILE_PSNID[TYPEINDX], "%s", file);
  ptrFile = TEMPLATES_FILE_PSNID[TYPEINDX] ;

  // open as text file just to see if it exists.
  fp = snana_openTextFile(1,PATH_TEMPLATES_PSNID,       // (I) public area
			  ptrFile,                    // (I) filename
			  FILE,                // (O) full filename 
			  &gzipFlag );         // (O) gzip flag

  if ( fp == NULL ) {
    sprintf(c1err,"Could not open template file:");
    sprintf(c2err,"%s", ptrFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  fclose(fp);

  // Now we know it exists;
  // read FITS grid and load SNGRID  in sngridtools.h
  fits_read_SNGRID(OPT_READ, FILE, &SNGRID_PSNID[TYPEINDX] );
  USEFLAG_TEMPLATES_PSNID[TYPEINDX] = 1 ;
    
  // make sure there are templates for each fitted filter (5/01/2014)
  psnid_check_filtlist_fit(TYPEINDX);

  printf("\n");

} // end of PSNID_READ_TEMPLATES

void psnid_read_templates__(int *TYPEINDX, char *file) {
  PSNID_READ_TEMPLATES( *TYPEINDX, file);
}


// ======================================================
void psnid_check_filtlist_fit(int TYPEINDX) {

  // Created 5/01/2014 by R.Kessler
  // Check each filter to fit and make sure there is 
  // a template define for this filter; abort if any
  // template-filters are missing.

  int   NFILT           = PSNID_INPUTS.NFILT ;
  char *FITFILTER_LIST  = PSNID_INPUTS.CFILTLIST ;  // filters to fit
  char *GRIDFILTER_LIST = SNGRID_PSNID[TYPEINDX].FILTERS ;
  int  ifilt, NBAND_MISS ;
  char band[4];
  char fnam[] = "psnid_check_filtlist_fit" ;

  // -------------- BEGIN -----------

  printf("\n");
  NBAND_MISS = 0;

  for(ifilt=0; ifilt < NFILT; ifilt++ ) {
    sprintf(band,"%c", FITFILTER_LIST[ifilt]) ;
    if ( strstr(GRIDFILTER_LIST,band) == NULL ) {
      printf(" ERROR: PSNID TEMPLATES MISSING FIT BAND = '%s' \n", band);
      NBAND_MISS++ ;
    }
  }

  if ( NBAND_MISS > 0 ) {
    sprintf(c1err,"FILTLIST_FIT = %s in &PSNIDINP", FITFILTER_LIST);    
    sprintf(c2err,"but template filters include '%s'", GRIDFILTER_LIST );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
  }

} // end of psnid_check_filtlist_fit

// ======================================================
void PSNID_INIT_TEMPLATES(char *templates_nonIa_list, 
			  char *templates_nonIa_ignore) {

  // inits and crosschecks after all templates have been read.
  // Input argument is 
  //   - optional subset of NONIA-template names to use
  //   - optional subset of NONIA-template names to ignore.
  // Example:
  //  templates_nonIa_ignore = 'SN000018  SN002000'
  //

  int ifilt, ifiltobs ;
  int USE_SNIA  = USEFLAG_TEMPLATES_PSNID[TYPEINDX_SNIA_PSNID];
  int USE_NON1A = USEFLAG_TEMPLATES_PSNID[TYPEINDX_NONIA_PSNID];
  
  char *flt1, *flt2;
  char fnam[] = "PSNID_INIT_TEMPLATES" ;

  // ------------- BEGIN ---------------

  // psnid_dump_SNGRID(); // debug dump only

  // abort if filters don't match for Ia and nonIa
  if ( USE_SNIA && USE_NON1A ) {
    flt1 = SNGRID_PSNID[TYPEINDX_SNIA_PSNID].FILTERS ;
    flt2 = SNGRID_PSNID[TYPEINDX_NONIA_PSNID].FILTERS ;

    if ( strcmp(flt1,flt2) != 0 ) {
      sprintf(c1err,"Filter mismatch: %s(Ia) != %s(NONIA)", 
	      flt1, flt2);
      sprintf(c2err,"Remake GRID-template file or choose another.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }
  }

  // abort if Trest range does not match ...

  // loop over FIT filters and set USEFILT_GRIDOFF to be the offset
  // Note that user FILTLIST_FIT is not necessarily in the
  // same order as the filter-list in the template files.

  int OFF, NBIN, TYPEINDX; 

  for ( TYPEINDX=1; TYPEINDX <=2; TYPEINDX++ ) {
    OFF  = 0;
    NBIN = SNGRID_PSNID[TYPEINDX].NBIN[IPAR_GRIDGEN_TREST] ;
    if ( USEFLAG_TEMPLATES_PSNID[TYPEINDX] == 0 ) { continue ; }
    for ( ifilt=0; ifilt < PSNID_INPUTS.NFILT; ifilt++ ) {
      ifiltobs = PSNID_INPUTS.IFILTLIST[ifilt];
      PSNID_INPUTS.USEFILT_GRIDOFF[TYPEINDX][ifiltobs] = OFF ;
      OFF += NBIN ;
    }
  }
  
  psnid_init_nonIa_ignore(templates_nonIa_list,templates_nonIa_ignore);


  // normalize wgts so that sum = 1
  SNGRID_PSNID[TYPEINDX_NONIA_PSNID].FRAC_PEC1A=0.0; // Aug 15 2016
  renorm_wgts_SNGRID(&SNGRID_PSNID[TYPEINDX_NONIA_PSNID]) ;

  // print info about each nonIa template
  dump_SNGRID(&SNGRID_PSNID[TYPEINDX_NONIA_PSNID]) ;

  // create huge mag -> fluxcal table
  psnid_mktable_pogson2fluxcal(); 

  //  psnid_dumpLC_SNGRID();

}    // end of  PSNID_INIT_TEMPLATES


void psnid_init_templates__(char *templates_nonIa_list, 
			    char *templates_nonIa_ignore) {
  PSNID_INIT_TEMPLATES(templates_nonIa_list,templates_nonIa_ignore);
}


// ==================================================
void  psnid_init_nonIa_ignore(char *templates_nonIa_list,
			      char *templates_nonIa_ignore) {

  // check list of nonIa templates to ignore.
  // Set ITYPE ->  -ITYPE as a flag to ignore it.
  // If non-existant template is given, abort.
  //
  // Sep 7 2013: add arg templates_nonIa_list

  char *ptrtok, name[200];
  int  IMATCH, TYPEINDX ;
  char fnam[] = "psnid_init_nonIa_ignore" ;

  // -------------- BEGIN -------------


  TYPEINDX = TYPEINDX_NONIA_PSNID ;

  if ( strlen(templates_nonIa_ignore) < 2 ) { goto CHECK_LIST ; }

  ptrtok = strtok(templates_nonIa_ignore," ");
  while ( ptrtok != NULL  ) {

    sprintf(name, "%s", ptrtok );
    if ( strlen(name) > 1 ) {
      IMATCH = psnid_match_NONIA_NAME(name);
	
      if ( IMATCH > 0 ) 
	{ SNGRID_PSNID[TYPEINDX].NON1A_ITYPE_AUTO[IMATCH] *= -1 ; }
      else {
	sprintf(c1err,"Could not find NONIA_NAME = '%s'", name);
	sprintf(c2err,"Check namelist variable TEMPLATES_NONIA_IGNORE");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
    }
    ptrtok = strtok(NULL, " ");
  }
  

 CHECK_LIST:

  // check optional subset; ignore those not on the list

  if ( strlen(templates_nonIa_list) < 2 ) { return ; }

  int IPAR, NTMPL, i, NLIST ;
  int USEME[MXGRIDGEN] ;

  // identify those on the NONIA_LIST; fill USEME array.
  NLIST = 0;
  for(i=0; i < MXGRIDGEN; i++ ) { USEME[i] = 0; }

  ptrtok = strtok(templates_nonIa_list," ");
  while ( ptrtok != NULL  ) {

    sprintf(name, "%s", ptrtok );
    if ( strlen(name) > 1 ) {
      IMATCH = psnid_match_NONIA_NAME(name);

      if ( IMATCH == 0 ) {
	sprintf(c1err,"Could not find NONIA_NAME = '%s'", name);
	sprintf(c2err,"Check namelist variable TEMPLATES_NONIA_LIST");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      NLIST++ ;
      USEME[IMATCH] = 1 ;
    }
    ptrtok = strtok(NULL, " ");
  }

  // Now turn off templates that are not in the NONIA_LIST,
  // but only if there is something on the list.
  if ( NLIST > 0 ) {
    IPAR     = IPAR_GRIDGEN_SHAPEPAR ;
    TYPEINDX = TYPEINDX_NONIA_PSNID ;
    NTMPL    = SNGRID_PSNID[TYPEINDX].NBIN[IPAR];
    
    for ( i=1; i <= NTMPL ; i++ ) {
      if ( USEME[i] == 0 ) 
	{ SNGRID_PSNID[TYPEINDX].NON1A_ITYPE_AUTO[i] *= -1 ; }
    }
  }
  
} // end of psnid_init_nonIa_ignore



// ===================================
int psnid_match_NONIA_NAME(char *name) {

  // return nonIa index whose NONIA_NAME matches input name.

  int i, NTMPL, TYPEINDX, IPAR ;
  char *NAME ;
  //  char fnam[] = "psnid_match_NONIA_NAME" ;

  // ------------- BEGIN -------------

  IPAR     = IPAR_GRIDGEN_SHAPEPAR ;
  TYPEINDX = TYPEINDX_NONIA_PSNID ;
  NTMPL    = SNGRID_PSNID[TYPEINDX].NBIN[IPAR];

  for ( i=1; i <= NTMPL ; i++ ) {
    NAME = SNGRID_PSNID[TYPEINDX].NON1A_NAME[i] ;
    if ( strcmp(NAME,name) == 0 ) {
      return i ;
    }
  }

  // if we get here, return error code.
  return -9;

} // end of psnid_match_NONIA_NAME

// ===========================================
void get_template_lc(int TYPEINDX, int iz, int ic1, int ic2, 
		     int ishape, int ifiltobs, int optDump,
		     int *NEPOCH,  double *Trest,   // return
		     double *mag, double *magerr) // return
{

  // Return template light curve for inputs
  //   TYPEINDX = 1,2 for Ia,NONIa
  //   iz = logz index
  //   ic1 = color index (c or AV)
  //   ic2 = color law index (RV or beta)
  //   ishape = shape param (Delta, DM15, x1, NONIA-index)
  //   ifiltobs = absolute filter index
  // Returns
  //   NEPOCH = number of epochs
  //   Trest  = list of Trest (days)
  //   mag    = list of observer-frame mags
  //
  // Note that observer-frame mags are returned even though
  // you pass rest-frame epochs ; this is because the underlying
  // grid is uniform in TREST but not in Tobs.
  //

  int ILC, IPTROFF, IOFF_FILT, NEP, ep, ep0, LDMP ;
  short I2TMP, *I2PTR_MAG, *I2PTR_ERR  ;
  float *R4PTR_EP ;
  char fnam[] = "get_template_lc" ;

  // ------------ BEGIN --------------

  // start with sanity checks

  psnid_binCheck( TYPEINDX, IPAR_GRIDGEN_LOGZ,     iz   );
  psnid_binCheck( TYPEINDX, IPAR_GRIDGEN_COLORPAR, ic1  );
  psnid_binCheck( TYPEINDX, IPAR_GRIDGEN_COLORLAW, ic2  );
  psnid_binCheck( TYPEINDX, IPAR_GRIDGEN_SHAPEPAR, ishape);

  LDMP = ( optDump != 0 && NLCDUMP_PSNID < MXLCDUMP_PSNID ) ;

  if ( PSNID_INPUTS.USEFILT[ifiltobs] == 0 ) {
    sprintf(c1err,"FIlter = '%c' (%d) is not defined for fit.",
	    FILTERSTRING[ifiltobs], ifiltobs);
    sprintf(c2err,"Check &PSNIDINP namelist string FILTLIST_FIT");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }

  ILC = 1  
    + (iz  -1)  * SNGRID_PSNID[TYPEINDX].ILCOFF[IPAR_GRIDGEN_LOGZ]
    + (ic1 -1)  * SNGRID_PSNID[TYPEINDX].ILCOFF[IPAR_GRIDGEN_COLORPAR]
    + (ic2 -1)  * SNGRID_PSNID[TYPEINDX].ILCOFF[IPAR_GRIDGEN_COLORLAW]
    + (ishape-1)* SNGRID_PSNID[TYPEINDX].ILCOFF[IPAR_GRIDGEN_SHAPEPAR]
    ; 

  IPTROFF    =   SNGRID_PSNID[TYPEINDX].PTR_GRIDGEN_LC[ILC] ;

  // make sure that IPTROFF is valid.
  // divide size by 4 because there are 2 2-byte words per row
  long int MXPTR;
  MXPTR     =  SNGRID_PSNID[TYPEINDX].SIZEOF_GRIDGEN/4 ;
  if ( IPTROFF < 0 || IPTROFF > MXPTR ) {
    print_preAbort_banner(fnam);
    printf("\t iz=%d  ic1=%d  ic2=%d  ishape=%d -> ILC=%d \n", 
           iz, ic1, ic2, ishape, ILC );     
    sprintf(c1err,"Invalid IPTROFF = %d for %s-GRID",
            IPTROFF,   TEMPLATETYPE_PSNID[TYPEINDX] );
    sprintf(c2err,"Must be between 0 and %ld", MXPTR);  // fix compile bug, RK Dec 2 2012
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }

  I2PTR_MAG  =  &SNGRID_PSNID[TYPEINDX].I2GRIDGEN_LCMAG[IPTROFF] ;
  I2PTR_ERR  =  &SNGRID_PSNID[TYPEINDX].I2GRIDGEN_LCERR[IPTROFF] ;
  R4PTR_EP   =  &SNGRID_PSNID[TYPEINDX].VALUE[IPAR_GRIDGEN_TREST][0] ;

  // make sure that 1st word is BEGIN-LC marker
  I2TMP = I2PTR_MAG[0] ;
  if ( I2TMP != MARK_GRIDGEN_LCBEGIN ) {
    sprintf(c1err,"First I*2 word of ILC=%d is %d .", ILC, I2TMP );
    sprintf(c2err,"But expected %d", MARK_GRIDGEN_LCBEGIN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // make sure that 2nd word is 8-bits of ILC
  I2TMP = I2PTR_MAG[1];
  if ( I2TMP != ( ILC & 127 ) ){
    sprintf(c1err,"2nd I*2 word of ILC=%d is %d .", ILC, I2TMP );
    sprintf(c2err,"But expected %d", (ILC & 127) );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  IOFF_FILT = NPADWD_LCBEGIN-1 
    + PSNID_INPUTS.USEFILT_GRIDOFF[TYPEINDX][ifiltobs];

  NEP = SNGRID_PSNID[TYPEINDX].NBIN[IPAR_GRIDGEN_TREST] ;

  if ( LDMP ) {
    NLCDUMP_PSNID++ ;
    printf("\n# DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD \n");
    printf("  DUMP TEMPLATE LC for %s  %c-band\n"
	   ,TEMPLATETYPE_PSNID[TYPEINDX] 
	   ,FILTERSTRING[ifiltobs]   );

    psnid_dumpLine(TYPEINDX, IPAR_GRIDGEN_LOGZ,     iz   );
    psnid_dumpLine(TYPEINDX, IPAR_GRIDGEN_COLORPAR, ic1  );
    psnid_dumpLine(TYPEINDX, IPAR_GRIDGEN_COLORLAW, ic2  );
    psnid_dumpLine(TYPEINDX, IPAR_GRIDGEN_SHAPEPAR, ishape );
    printf("\t ILC = %d   PTR_GENGRID_LC=%d   IOFF_FILT=%d\n", 
	   ILC, IPTROFF, IOFF_FILT );
    fflush(stdout);
  }

  // load output light curve for this filter.
  *NEPOCH = NEP ;
  for (ep=1; ep <= NEP ; ep++ ) {
    ep0 = ep - 1;
    Trest[ep0]  = (double)R4PTR_EP[ep] ;
    mag[ep0]    = (double)I2PTR_MAG[IOFF_FILT+ep] / MAGPACK_GRIDGEN ;
    magerr[ep0] = (double)I2PTR_ERR[IOFF_FILT+ep] / MAGPACK_GRIDGEN ;

    if ( LDMP ) {
      printf("\t DUMP: ep=%2d: Trest=%5.1f  mag=%6.3f +- %6.3f \n",
	     ep, Trest[ep0], mag[ep0], magerr[ep0] );
      fflush(stdout) ;
    }
  }


} // end of get_template_lc


// ==============================================
void psnid_dumpLC_SNGRID(void) {

  // few test dumps of templates for debugging
  double Trest[MXEP_PSNID], MAG[MXEP_PSNID], MAGERR[MXEP_PSNID];
  int    optDump = 1;
  int    iz, ic1, ic2, ishape, ifiltobs, NEP ;

  char fnam[] = "psnid_dumpLC_SNGRID" ;

  // ---------- BEGIN -------------

  iz = 1; ic1=1;ic2=1;  ishape=1;  ifiltobs=2;
  get_template_lc(1,iz,ic1,ic2,ishape,ifiltobs, optDump,
		  &NEP, Trest, MAG, MAGERR);


  iz = 78; ic1=10; ic2=1;  ishape=9;  ifiltobs=4;
  get_template_lc(1,iz,ic1,ic2,ishape,ifiltobs, optDump,
		  &NEP, Trest, MAG, MAGERR);

  debugexit(fnam);

} // end of psnid_dumpLC_SNGRID


// ====================================================
void psnid_binCheck(int TYPEINDX, int IPAR, int ibin) {

  int NBIN;
  char fnam[] = "psnid_binCheck";

  // ----------- BEGIN ---------

  NBIN = SNGRID_PSNID[TYPEINDX].NBIN[IPAR];

  if ( ibin < 1 || ibin > NBIN ) {
    sprintf(c1err,"Invalid ibin=%d for TYPEINDEX=%d(%s)  IPAR=%d",
	    ibin, TYPEINDX, SNGRID_PSNID[TYPEINDX].NAME[IPAR], IPAR );
    sprintf(c2err,"Valid bin range is 1 - %d", NBIN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

} // end of psnid_binCheck

// ==============================================
void psnid_dumpLine(int TYPEINDX, int IPAR, int ibin) {

  // dump one line summary for this IPAR and bin.

  //  char fnam[] = "psnid_dumpLine" ;

  printf("\t %-12s = %f  (ibin = %3d) \n"
	 ,SNGRID_PSNID[TYPEINDX].NAME[IPAR]
	 ,SNGRID_PSNID[TYPEINDX].VALUE[IPAR][ibin]
  	 ,ibin
	 );

} // end of psnid_dumpLine


// =================================================
void psnid_dumpBins_SNGRID(void)
{
  int i;
  //  char fnam[] = "psnid_dumpBins_SNGRID" ;

  // --------------- BEGIN ---------------

  printf("\n XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
  printf(" XXXXXXXXXX    SNGRID DUMP   XXXXXXXXXXX\n");

  
  // bins
  printf("\t NBIN_LOGZ(Ia,NONIa)     = %d , %d \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NBIN[IPAR_GRIDGEN_LOGZ]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_LOGZ]
	 );
  printf("\t NBIN_COLORPAR(Ia,NONIa) = %d , %d \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NBIN[IPAR_GRIDGEN_COLORPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_COLORPAR]
	 );
  printf("\t NBIN_COLORLAW(Ia,NONIa) = %d , %d \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NBIN[IPAR_GRIDGEN_COLORLAW]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_COLORLAW]
	 );
  printf("\t NBIN_SHAPEPAR(Ia,NONIa)  = %d , %d \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NBIN[IPAR_GRIDGEN_SHAPEPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_SHAPEPAR]
	 );
  printf("\t NBIN_FILTER(Ia,NONIa)   = %d , %d \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NBIN[IPAR_GRIDGEN_FILTER]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_FILTER]
	 );
  printf("\t NBIN_TREST(Ia,NONIa)    = %d , %d \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NBIN[IPAR_GRIDGEN_TREST]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_TREST]
	 );
  printf("\n");

  // binsizes
  printf("\t BINSIZE_LOGZ(Ia,NONIa)     = %f , %f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].BINSIZE[IPAR_GRIDGEN_LOGZ]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].BINSIZE[IPAR_GRIDGEN_LOGZ]
	 );
  printf("\t BINSIZE_COLORPAR(Ia,NONIa) = %f , %f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].BINSIZE[IPAR_GRIDGEN_COLORPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].BINSIZE[IPAR_GRIDGEN_COLORPAR]
	 );
  printf("\t BINSIZE_COLORLAW(Ia,NONIa) = %f , %f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].BINSIZE[IPAR_GRIDGEN_COLORLAW]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].BINSIZE[IPAR_GRIDGEN_COLORLAW]
	 );
  printf("\t BINSIZE_SHAPEPAR(Ia,NONIa)  = %f , %f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].BINSIZE[IPAR_GRIDGEN_SHAPEPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].BINSIZE[IPAR_GRIDGEN_SHAPEPAR]
	 );
  printf("\t BINSIZE_TREST(Ia,NONIa)    = %f , %f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].BINSIZE[IPAR_GRIDGEN_TREST]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].BINSIZE[IPAR_GRIDGEN_TREST]
	 );
  printf("\n");


  // parameter name/range
  printf("\t LOGZ_NAME(Ia,NONIa)  = %s , %s \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NAME[IPAR_GRIDGEN_LOGZ]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NAME[IPAR_GRIDGEN_LOGZ]
	 );
  printf("\t LOGZ_RANGE(Ia,NONIa) = %6.3f to %6.3f  , %6.3f to %6.3f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMIN[IPAR_GRIDGEN_LOGZ]
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMAX[IPAR_GRIDGEN_LOGZ]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMIN[IPAR_GRIDGEN_LOGZ]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMAX[IPAR_GRIDGEN_LOGZ]
	 );
  printf("\t    Z_RANGE(Ia,NONIa) = %6.3f to %6.3f  , %6.3f to %6.3f \n"
	 ,pow(10,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMIN[IPAR_GRIDGEN_LOGZ])
	 ,pow(10,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMAX[IPAR_GRIDGEN_LOGZ])
	 ,pow(10,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMIN[IPAR_GRIDGEN_LOGZ])
	 ,pow(10,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMAX[IPAR_GRIDGEN_LOGZ])
	 );
  printf("\n");

  printf("\t COLORPAR_NAME(Ia,NONIa)  = %s , %s \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NAME[IPAR_GRIDGEN_COLORPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NAME[IPAR_GRIDGEN_COLORPAR]
	 );
  printf("\t COLORPAR_RANGE(Ia,NONIa) = %6.3f to %6.3f  , %6.3f to %6.3f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMIN[IPAR_GRIDGEN_COLORPAR]
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMAX[IPAR_GRIDGEN_COLORPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMIN[IPAR_GRIDGEN_COLORPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMAX[IPAR_GRIDGEN_COLORPAR]
	 );
  printf("\n");

  printf("\t COLORLAW_NAME(Ia,NONIa)  = %s , %s \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NAME[IPAR_GRIDGEN_COLORLAW]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NAME[IPAR_GRIDGEN_COLORLAW]
	 );
  printf("\t COLORLAW_RANGE(Ia,NONIa) = %6.3f to %6.3f  , %6.3f to %6.3f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMIN[IPAR_GRIDGEN_COLORLAW]
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMAX[IPAR_GRIDGEN_COLORLAW]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMIN[IPAR_GRIDGEN_COLORLAW]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMAX[IPAR_GRIDGEN_COLORLAW]
	 );
  printf("\n");

  printf("\t SHAPEPAR_NAME(Ia,NONIa) = %s , %s \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NAME[IPAR_GRIDGEN_SHAPEPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NAME[IPAR_GRIDGEN_SHAPEPAR]
	 );  
  printf("\t SHAPE_RANGE(Ia,NONIa)   = %6.3f to %6.3f  , %6.3f to %6.3f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMIN[IPAR_GRIDGEN_SHAPEPAR]
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMAX[IPAR_GRIDGEN_SHAPEPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMIN[IPAR_GRIDGEN_SHAPEPAR]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMAX[IPAR_GRIDGEN_SHAPEPAR]
	 );
  printf("\n");

  printf("\t TREST_NAME(Ia,NONIa) = %s , %s \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].NAME[IPAR_GRIDGEN_TREST]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NAME[IPAR_GRIDGEN_TREST]
	 );
  printf("\t TREST_RANGE(Ia,NONIa) = %6.3f to %6.3f  , %6.3f to %6.3f \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMIN[IPAR_GRIDGEN_TREST]
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMAX[IPAR_GRIDGEN_TREST]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMIN[IPAR_GRIDGEN_TREST]
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].VALMAX[IPAR_GRIDGEN_TREST]
	 );
  printf("\n");

  printf("\t FILTER_NAME(Ia,NONIa) = %s , %s \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].FILTERS
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].FILTERS
	 );
  printf("\t SURVEY(Ia,NONIa)      = %s , %s \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].SURVEY
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].SURVEY
	 );
  printf("\t MODEL(Ia,NONIa)       = %s , %s \n"
	 ,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].MODEL
	 ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].MODEL
	 );
  printf("\n");

  // non-Ia types
  for (i=1; i<=SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_SHAPEPAR]; i++) {
        printf("\t NONIA lumin %2d   NONIA_INDEX = %2d   NONIA_NAME = %40s   NONIA_CTYPE = %20s\n"
    	   ,i
    	   ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NON1A_INDEX[i]
    	   ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NON1A_NAME[i]
    	   ,SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NON1A_CTYPE[i]);
  }

  printf(" XXXXXXXXXX    SNGRID DUMP END    XXXXXXXXXXX\n");
  printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");


  fflush(stdout);

} // end of  psnid_dumpBins_SNGRID


/************************************************************************/
void psnid_store_data(char *CCID, int NOBS, int *IFILT,
		      double *MJD, double *FLUXDATA, double *FLUXERR, 
		      double *FLUXSIM, 
		      double *REDSHIFT, double *REDSHIFT_ERR,
		      double MWEBV, double MWEBV_ERR, int SIM_NON1A_INDEX)
/*************************************************************************/
// Feb 18, 2013 RK - store light curve in structure
//                   called from  PSNID_[METHOD]_DOFIT interface.
// Mar 10 2017 RK - pass and store SIM_NON1A_INDEX
{

  int obs, i, MEM_D, MEM_I ;

  // ----------- BEGIN ----------

  // allocate memory; but free current memory if already allocated.
  if ( DATA_PSNID_DOFIT.NOBS > 0 ) {
    free(DATA_PSNID_DOFIT.MJD) ;
    free(DATA_PSNID_DOFIT.FLUXDATA) ;
    free(DATA_PSNID_DOFIT.FLUXERR) ;
    free(DATA_PSNID_DOFIT.FLUXSIM) ;
    free(DATA_PSNID_DOFIT.DUMERR0) ;
    free(DATA_PSNID_DOFIT.IFILT) ;
    free(DATA_PSNID_DOFIT.IFILTOBS) ;
  }

  MEM_D = NOBS * sizeof(double);
  MEM_I = NOBS * sizeof(int); 
  DATA_PSNID_DOFIT.MJD       = (double*)malloc(MEM_D) ;
  DATA_PSNID_DOFIT.FLUXDATA  = (double*)malloc(MEM_D) ;
  DATA_PSNID_DOFIT.FLUXERR   = (double*)malloc(MEM_D) ;
  DATA_PSNID_DOFIT.FLUXSIM   = (double*)malloc(MEM_D) ;
  DATA_PSNID_DOFIT.DUMERR0   = (double*)malloc(MEM_D) ;
  DATA_PSNID_DOFIT.IFILT     = (int   *)malloc(MEM_I) ;
  DATA_PSNID_DOFIT.IFILTOBS  = (int   *)malloc(MEM_I) ;

// ------------------
  sprintf(DATA_PSNID_DOFIT.CCID, "%s", CCID);
  DATA_PSNID_DOFIT.NOBS = NOBS ;

  for ( i=0; i < 3 ; i++ ) {
    DATA_PSNID_DOFIT.REDSHIFT[i]     = REDSHIFT[i] ;
    DATA_PSNID_DOFIT.REDSHIFT_ERR[i] = REDSHIFT_ERR[i] ;
  }

  DATA_PSNID_DOFIT.MWEBV     = MWEBV ;
  DATA_PSNID_DOFIT.MWEBV_ERR = MWEBV_ERR ;

  DATA_PSNID_DOFIT.SIM_NON1A_INDEX = SIM_NON1A_INDEX ;

  for(obs=0; obs < NOBS; obs++ ) {
    DATA_PSNID_DOFIT.MJD[obs]       = MJD[obs] ;
    DATA_PSNID_DOFIT.FLUXDATA[obs]  = FLUXDATA[obs] ;
    DATA_PSNID_DOFIT.FLUXERR[obs]   = FLUXERR[obs] ;
    DATA_PSNID_DOFIT.FLUXSIM[obs]   = FLUXSIM[obs] ;
    DATA_PSNID_DOFIT.DUMERR0[obs]   = 0.0 ;
    DATA_PSNID_DOFIT.IFILT[obs]     = IFILT[obs] ;
    DATA_PSNID_DOFIT.IFILTOBS[obs]  = 
      PSNID_INPUTS.IFILTLIST[IFILT[obs]]; // absolute index
  }

} // end of psnid_store_data



/************************************************************************/
void psnid_dumpInput_data(char *CCID, int NOBS, int *IFILTOBS, 
			  double *MJD, double *FLUX, double *FLUXERR,
			  double *REDSHIFT, double *REDSHIFT_ERR,
			  double MWEBV, double MWEBVERR)
/*************************************************************************/
// Note to Masao (RK, Feb 8 2013): should remove args to this function 
// and use DATA_PSNID_DOFIT struct filled in psnid_store_data().

{
  int i, iobs, ifiltobs ;
  char cfilt[2];

  //  char fnam[] = "psnid_dumpInput_data" ;

  // --------------- BEGIN ---------------

  printf("\n XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
  printf(" XXXXXXXXXX   BEGIN: psnid_dumpInput_data   XXXXXXXXXX\n");
  fflush(stdout);

  for (i=0; i<3; i++) {
    printf("        REDSHIFT[%d] = %9.6f  REDSHIFT_ERR[%d] = %9.6f\n",
	   i, REDSHIFT[i], i, REDSHIFT_ERR[i]);
    fflush(stdout);
  }
  printf("              MWEBV = %9.6f         MWEBVERR = %9.6f\n",
	 MWEBV, MWEBVERR);
  fflush(stdout);

  for ( iobs=0; iobs < NOBS; iobs++ ) {
    ifiltobs = IFILTOBS[iobs] ;
    //    sprintf(cfilt,"%c", FILTERSTRING[ifiltobs] ) ;
    sprintf(cfilt,"%c", PSNID_INPUTS.CFILTLIST[ifiltobs]) ;

    printf("\t xxx iobs=%3d  CFILTLIST[IFILTOBS[iobs]]=%s  "
	      "IFILTOBS[iobs]=%d  IFILTLIST[IFILTOBS[iobs]]=%d  "
	   "MJD=%9.3f   FLUX=%8.3f +- %8.3f    S/N = %8.3f\n",
	   iobs, cfilt, IFILTOBS[iobs], PSNID_INPUTS.IFILTLIST[ifiltobs],
	   MJD[iobs], FLUX[iobs], FLUXERR[iobs],
	   FLUX[iobs]/FLUXERR[iobs]);
    fflush(stdout);

  }

  printf(" XXXXXXXXXX   END: psnid_dumpInput_data   XXXXXXXXXX\n");
  printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
  fflush(stdout);

  return;
}
// end of psnid_dumpInput_data



// ===================================================
void psnid_pogson2mJy(double mag, double mage,
		      double *flux, double *fluxe) {

  // convert mag and mage into flux in mJy.

  int LMAG, LERR ;

  LMAG = (mag  <= PSNID_GOODMAG_HI    && mag  >= PSNID_GOODMAG_LO    ) ;
  LERR = (mage <= PSNID_GOODMAGERR_HI && mage >= PSNID_GOODMAGERR_LO ) ;

  if ( LMAG && LERR ) {
    *flux   = 3.631e9*pow(100.,-0.2*mag);
    *fluxe  = *flux * mage;
  } else {
    *flux   = -9.0;
    *fluxe  = -9.0;
  }

  return;
}
// end of psnid_pogson2mJy



// ======================================
void psnid_mktable_pogson2fluxcal(void) {

  // Created Sep 8 2013 by R.Kessler
  // To save CPU, create massive lookup table for mag -> fluxcal
 
  int    NBIN, MEM, i ;
  double MAGRANGE, MAGBIN, XMB, xi, mag, mage, fluxcal, fluxcale ;

  // -------------- BEGIN ----------

  SIZEOF_PSNID_FLUXCAL_LOOKUP = 0 ;

  MAGRANGE = PSNID_GOODMAG_HI - PSNID_GOODMAG_LO ;
  NBIN     = (int)(MAGRANGE * PSNID_MAGSCALE_for_LOOKUP);
  MAGBIN   = MAGRANGE / (double)NBIN ;

  MEM  = NBIN * sizeof(double) ;
  PSNID_FLUXCAL_LOOKUP = (double*)malloc(MEM);

  XMB = (double)(MEM)/1.0E6;

  printf("\n Create MAG->FLUXCAL Lookup with %d bins (%.2f MB); magBin=%f\n",
	 NBIN, XMB, MAGBIN ); fflush(stdout);

  mage = 0.01 ;  // any positive value
  for(i=0; i<NBIN; i++ ) {
    xi = (double)i;
    mag = PSNID_GOODMAG_LO + ( xi * MAGBIN ) ;
    psnid_pogson2fluxcal(mag, mage, &fluxcal, &fluxcale );
    PSNID_FLUXCAL_LOOKUP[i] = fluxcal ;
  }


  // load global SIZE telling the pogson2fluxcal 
  // function to use the newly-created lookup.

  SIZEOF_PSNID_FLUXCAL_LOOKUP =  MEM ; // store in global


} // psnid_mktable_pogson2fluxcal

// ===================================================
void psnid_pogson2fluxcal(double mag, double mage,
			  double *fluxcal, double *fluxcale)  {

  // convert mag and mage into FLUXCAL and flux-error (snana flux units)
  // Sep 8 2013 RK 
  //   - check option to use fast lookup table 
  //   - fix fluxcale calc with mage -> mage * .921
  //
  // Aug 28, 2014: for undefined model set flux=fluxe = 0 (instead of -9)
  //

  int    LMAG, LERR, ibin ;
  double arg;
  double ZEROPOINT_FLUXCAL = 27.5 ;

  LMAG = (mag  <= PSNID_GOODMAG_HI    && mag  >= PSNID_GOODMAG_LO    ) ;
  LERR = (mage <= PSNID_GOODMAGERR_HI && mage >= PSNID_GOODMAGERR_LO ) ;

  if (LMAG && LERR ) {

    if ( SIZEOF_PSNID_FLUXCAL_LOOKUP == 0 ) {
      // exact calculation 
      arg        = mag - ZEROPOINT_FLUXCAL ;
      *fluxcal   = pow(10.0,-0.4*arg);
    }
    else {
      // use lookup to avoid 'pow' function
      arg  = (mag - PSNID_GOODMAG_LO) * PSNID_MAGSCALE_for_LOOKUP ;
      ibin = (int)arg ;
      *fluxcal = PSNID_FLUXCAL_LOOKUP[ibin];

      /* xxxxxxxx
      double F ;
      arg  = mag - ZEROPOINT_FLUXCAL ;   F  = pow(10.0,-0.4*arg);
      printf(" xxx F(lookup)/F(true) = %le/%le = %f \n",
	     *fluxcal, F, *fluxcal/F );	  xxxxxxxxx */
    }

    // *fluxcale  = (*fluxcal) * (mage) ; // OLD approx.
    *fluxcale  = (*fluxcal) * (mage * .921) ; // .921 = ln(10)/2.5

  } else {
    // model is undefined; assume zero
    *fluxcal   = 0.0;
    *fluxcale  = 0.0;
  }

  return;
}
// end of psnid_pogson2fluxcal


/**************************************************************************/
/**************************************************************************/
/**************                                           *****************/
/**************              Numerical Recipes            *****************/
/**************                                           *****************/
/**************************************************************************/
/**************************************************************************/

/***************************/
/*****    nrutil.c    ******/
/***************************/

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

#define NR_END 1
#define FREE_ARG char*


/* Numerical Recipes standard error handler */
void nrerror(char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

unsigned long long *llvector(long nl, long nh)
/* allocate an unsigned long long vector with subscript range v[nl..nh] */
{
	unsigned long long *v;

	v=(unsigned long long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long long)));
	if (!v) nrerror("allocation failure in llvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate an int 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	int ***t;

	/* allocate pointers to pointers to rows */
	t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int**)));
	if (!t) nrerror("allocation failure 1 in i3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int*)));
	if (!t[nrl]) nrerror("allocation failure 2 in i3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in i3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in d3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

double ****d4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nel, long neh)
/* allocate a double 4tensor with range t[nrl..nrh][ncl..nch][ndl..ndh][nel..neh] */
{
	long i,j,k,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,ndep2=neh-nel+1;
	double ****t;

	/* allocate pointers to pointers to rows */
	t=(double ****) malloc((size_t)((nrow+NR_END)*sizeof(double***)));
	if (!t) nrerror("allocation failure 1 in d4tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double**)));
	if (!t[nrl]) nrerror("allocation failure 2 in d4tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double **) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double*)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in d4tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	
	/* allocate rows and set pointers to them */
	t[nrl][ncl][ndl]=(double *) malloc((size_t)((nrow*ncol*ndep*ndep2+NR_END)*sizeof(double)));
	if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in d4tensor()");
	t[nrl][ncl][ndl] += NR_END;
	t[nrl][ncl][ndl] -= nel;

   for(k=ndl+1;k<=ndh;k++) t[nrl][ncl][k]=t[nrl][ncl][k-1]+ndep2;
   for(j=ncl+1;j<=nch;j++) {
		t[nrl][j]=t[nrl][j-1]+ndep;
		t[nrl][j][ndl]=t[nrl][j-1][ndl]+ndep*ndep2;
		for(k=ndl+1;k<=ndh;k++) t[nrl][j][k]=t[nrl][j][k-1]+ndep2;}
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		t[i][ncl][ndl]=t[i-1][ncl][ndl]+ncol*ndep*ndep2;
		for(k=ndl+1;k<=ndh;k++) t[i][ncl][k]=t[i][ncl][k-1]+ndep2;
		for(j=ncl+1;j<=nch;j++) {
		   t[i][j]=t[i][j-1]+ndep;
		   t[i][j][ndl]=t[i][j-1][ndl]+ndep*ndep2;
		   for(k=ndl+1;k<=ndh;k++) t[i][j][k]=t[i][j][k-1]+ndep2;}
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_llvector(unsigned long long *v, long nl, long nh)
/* free an unsigned long vector allocated with llvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free an int i3tensor allocated by i3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_d4tensor(double ****t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh, long nel, long neh)
/* free a double d4tensor allocated by d4tensor() */
{
   free((FREE_ARG) (t[nrl][ncl][ndl]+nel-NR_END));
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

#else /* ANSI */
/* traditional - K&R */

#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((unsigned int)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
long nch,ncl,nrh,nrl;
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
	if (!m)	nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((unsigned int)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

double ***d3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((unsigned int)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in d3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((unsigned
	int)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

double ****d4tensor(nrl,nrh,ncl,nch,ndl,ndh,nel,neh)
long nrl,nrh,ncl,nch,ndl,ndh,nel,neh;
/* allocate a double 4tensor with range t[nrl..nrh][ncl..nch][ndl..ndh][nel..neh] */
{
	long i,j,k,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,ndep2=neh-nel+1;
	double ****t;

	/* allocate pointers to pointers to rows */
	t=(double ****) malloc((unsigned int)((nrow+NR_END)*sizeof(double***)));
	if (!t) nrerror("allocation failure 1 in d4tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double ***) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double**)));
	if (!t[nrl]) nrerror("allocation failure 2 in d4tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double **) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(double*)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in d4tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	
	/* allocate rows and set pointers to them */
	t[nrl][ncl][ndl]=(double *) malloc((unsigned int)((nrow*ncol*ndep*ndep2+NR_END)*sizeof(double)));
	if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in d4tensor()");
	t[nrl][ncl][ndl] += NR_END;
	t[nrl][ncl][ndl] -= nel;

   for(k=ndl+1;k<=ndh;k++) t[nrl][ncl][k]=t[nrl][ncl][k-1]+ndep2;
   for(j=ncl+1;j<=nch;j++) {
		t[nrl][j]=t[nrl][j-1]+ndep;
		t[nrl][j][ndl]=t[nrl][j-1][ndl]+ndep*ndep2;
		for(k=ndl+1;k<=ndh;k++) t[nrl][j][k]=t[nrl][j][k-1]+ndep2;}
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		t[i][ncl][ndl]=t[i-1][ncl][ndl]+ncol*ndep*ndep2;
		for(k=ndl+1;k<=ndh;k++) t[i][ncl][k]=t[i][ncl][k-1]+ndep2;
		for(j=ncl+1;j<=nch;j++) {
		   t[i][j]=t[i][j-1]+ndep;
		   t[i][j][ndl]=t[i][j-1][ndl]+ndep*ndep2;
		   for(k=ndl+1;k<=ndh;k++) t[i][j][k]=t[i][j][k-1]+ndep2;}
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(v,nl,nh)
float *v;
long nh,nl;
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(v,nl,nh)
double *v;
long nh,nl;
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
long nch,ncl,nrh,nrl;
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
long nch,ncl,nrh,nrl;
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
long nch,ncl,nrh,nrl;
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
float ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_d3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
double ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a double d3tensor allocated by d3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_d4tensor(t,nrl,nrh,ncl,nch,ndl,ndh,nel,neh)
double ****t;
long nrl,nrh,ncl,nch,ndl,ndh,nel,neh;
/* free a double d4tensor allocated by d4tensor() */
{
   	free((FREE_ARG) (t[nrl][ncl][ndl]+nel-NR_END));
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}


#endif /* ANSI */


/*************************/
/*****    hunt.c    ******/
/*************************/


/* xxxxx  does not work with c++ ??
void hunt(xx,n,x,jlo)
double xx[], x ;
int n, *jlo ;
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void hunt(double xx[], int n, double x, int *jlo)
{
   int jm,jhi,inc,ascnd;

   ascnd=(xx[n] > xx[1]);
   if (*jlo <= 0 || *jlo > n) {
      *jlo=0;
      jhi=n+1;
   } else {
      inc=1;
      if ( (x >= xx[*jlo]) == ascnd) {
         if (*jlo == n) return;
         jhi=(*jlo)+1;
         while ( (x >= xx[jhi]) == ascnd) {
            *jlo=jhi;
            inc += inc;
            jhi=(*jlo)+inc;
            if (jhi > n) {
               jhi=n+1;
               break;
            }
         }
      } else {
         if (*jlo == 1) {
            *jlo=0;
            return;
         }
         jhi=(*jlo);
         *jlo -= 1;
         while ( (x < xx[*jlo]) == ascnd) {
            jhi=(*jlo);
            inc += inc;
            *jlo=jhi-inc;
            if (*jlo < 1) {
               *jlo=0;
               break;
            }
         }
      }
   }
   while (jhi-(*jlo) != 1) {
      jm=(jhi+(*jlo)) >> 1;
      if ( (x > xx[jm]) == ascnd)
         *jlo=jm;
      else
         jhi=jm;
   }
}



/***************************/
/*****    indexx.c    ******/
/***************************/

/*
void indexx(n,arrin,indx)
int n,indx[];
double arrin[];
*/
void indexx(int n, double arrin[], int indx[])
{
   int l,j,ir,indxt,i;
   double q;

   for (j=1;j<=n;j++) indx[j]=j;
   l=(n >> 1) + 1;
   ir=n;
   for (;;) {
      if (l > 1)
         q=arrin[(indxt=indx[--l])];
      else {
         q=arrin[(indxt=indx[ir])];
         indx[ir]=indx[1];
         if (--ir == 1) {
            indx[1]=indxt;
            return;
         }
      }
      i=l;
      j=l << 1;
      while (j <= ir) {
         if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
         if (q < arrin[indx[j]]) {
            indx[i]=indx[j];
            j += (i=j);
         }
         else j=ir+1;
      }
      indx[i]=indxt;
   }
}


/*************************/
/*****    ran1.c    ******/
/*************************/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


/*************************/
/*****    ran2.c    ******/
/*************************/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


/*************************/
/****    gasdev.c    *****/
/*************************/

//#include <math.h>

float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


/********************      End of Numerical Recipes     ***********************/
/******************************************************************************/
