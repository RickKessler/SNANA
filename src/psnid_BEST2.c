/***********************************************
   Created Jun 2019 by R.Kessler and M.Sako

   
 ***********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_sf_gamma.h>

#include "fitsio.h"
#include "sntools.h"      // snana stuff
#include "sntools_grid.h"
#include "sntools_output.h"
#include "sntools_nearnbr.h"
#include "psnid_tools.h"


// Delcare functions
// Upper case name are called by fortran/snana program psnid
// Lower case names are internal.

void PSNID_BEST2_INIT(void);
void psnid_best2_init__(void) 
{  PSNID_BEST2_INIT(); }


void PSNID_BEST2_INIT_SNTABLE(int OPT, char *TEXTFORMAT, int LSIM) ;
void psnid_best2_init_sntable__(int *OPT, char *TEXTFORMAT, int *LSIM) 
{ PSNID_BEST2_INIT_SNTABLE(*OPT,TEXTFORMAT,*LSIM); }


void PSNID_BEST2_UPDATE_OUTPUT(char *CCID) ;
void psnid_best2_update_output__(char *CCID) 
{  PSNID_BEST2_UPDATE_OUTPUT(CCID) ; }


int PSNID_BEST2_DOFIT(char *CCID, int NOBS, int *IFILTOBS, 
		      double *MJD, double *FLUX, double *FLUXERR, 
		      double *REDSHIFT, double *REDSHIFT_ERR, 
		      double MWEBV, double MWEBVERR, int SIM_NON1A_INDEX);
int psnid_best2_dofit__(char *CCID, int *NOBS, int *IFILTOBS, 
			double *MJD, double *FLUX, double *FLUXERR, 
			double *REDSHIFT, double *REDSHIFT_ERR, 
			double *MWEBV, double *MWEBVERR, int *SIM_NON1A_INDEX) 
{
  return PSNID_BEST2_DOFIT(CCID, *NOBS, IFILTOBS, MJD, FLUX, FLUXERR, 
			   REDSHIFT, REDSHIFT_ERR, *MWEBV, *MWEBVERR,
			   *SIM_NON1A_INDEX );
}


int PSNID_BEST2_GET_BESTTYPE(void);
int psnid_best2_get_besttype__(void) {
  return PSNID_BEST2_GET_BESTTYPE();
}


void SNLCPLOT_PSNID_BEST2(int iplot) ;    // standard LC+fit plots
void snlcplot_psnid_best2__(int *iplot) 
{  SNLCPLOT_PSNID_BEST2(*iplot); }

void MONPLOT_PSNID_BEST2(int iplot) ;    // arbitrary monitor plots
void monplot_psnid_best2__(int *iplot)
{ MONPLOT_PSNID_BEST2(*iplot); }



// ===================================================
// ========= BEGIN FUNCTIONS =========================
// ===================================================

void PSNID_BEST2_INIT(void) {
  
  char fnam[] = "PSNID_BEST2_INIT" ;

  // --------------- BEGIN --------------

  printf(" xxx Hello from %s\n", fnam); fflush(stdout);

  return ;

} // end PSNID_BEST2_INIT



void PSNID_BEST2_INIT_SNTABLE(int OPT, char *TEXTFORMAT, int LSIM) {

  char fnam[] = "PSNID_BEST2_INIT_SNTABLE" ;

  // ---------------- BEGIN -----------------

    printf(" xxx Hello from %s: OPT=%d  FMT=%s\n", 
	   fnam, OPT, TEXTFORMAT) ; fflush(stdout);

  return;

} // end PSNID_BEST2_INIT_SNTABLE


void PSNID_BEST2_UPDATE_OUTPUT(char *CCID) {

  //  char fnam[] = "PSNID_BEST2_UPDATE_OUTPUT" ;

  // ------------ BEGIN -------------

  return ;

} // end PSNID_BEST2_UPDATE_OUTPUT


int PSNID_BEST2_DOFIT(char *CCID, int NOBS, int *IFILTOBS, 
		      double *MJD, double *FLUX, double *FLUXERR, 
		      double *REDSHIFT, double *REDSHIFT_ERR, 
		      double MWEBV, double MWEBVERR, int SIM_NON1A_INDEX) {

  int ERRFLAG = 0;
  char fnam[] = "PSNID_BEST2_DOFIT" ;

  // ----------- BEGIN ------------

  printf(" xxx Hello: %s for CCID = %s \n", fnam, CCID);

  return(ERRFLAG);

} // end PSNID_BEST2_DOFIT


int PSNID_BEST2_GET_BESTTYPE(void) {

  int TYPE = 1;
  char fnam[] = "PSNID_BEST2_GET_BESTTYPE" ;

  // ----------- BEGIN -----------

  printf(" xxx Hello: %s returns TYPE=%d \n", fnam, TYPE);

  return(TYPE);

} // end PSNID_BEST2_GET_BESTTYPE


void MONPLOT_PSNID_BEST2(int iplot) {

  //  char fnam[] = "MONPLOT_PSNID_BEST2";
  // ----------- BEGIN -------------

  return ;

} // end MONPLOT_PSNID_BEST2


void SNLCPLOT_PSNID_BEST2(int iplot) {

  //  char fnam[] = "SNLCPLOT_PSNID_BEST2";

  // ----------- BEGIN -------------

  return ;

} // end SNLCPLOT_PSNID_BEST2

// END:
