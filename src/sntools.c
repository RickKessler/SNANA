// sntools.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

#include <gsl/gsl_rng.h>      // for Poisson generator
#include <gsl/gsl_randist.h>  // idem

#include "sntools.h"

#include "eispack.h"
#include "eispack.c"

/*********************************************************
**********************************************************

               Utilities for SN Analysis
            (see declarations in sntools.h)

**********************************************************
**********************************************************/


void init_obs_atFLUXMAX(int OPTMASK, double *PARLIST, int VBOSE) {

  // May 24 2019
  // Initialize inputs to estimate PEAKMJD for each event,
  // using only flux and fluxerr.
  //

  char cmethod[100];
  char fnam[] = "init_obs_atFLUXMAX" ;

  // ------------- BEGIN ------------

  if ( OPTMASK == 0 ) { return ; }

  INPUTS_OBS_atFLUXMAX.OPTMASK        = OPTMASK ;
  INPUTS_OBS_atFLUXMAX.MJDWIN         = PARLIST[0];
  INPUTS_OBS_atFLUXMAX.SNRCUT         = PARLIST[1];
  INPUTS_OBS_atFLUXMAX.SNRCUT_BACKUP  = PARLIST[2];

  if ( VBOSE ) {
    if ( (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX) > 0 )  { 
      sprintf(cmethod,"max-flux");  
    }
    else if ( (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX2) > 0 )  { 
      sprintf(cmethod, "Fmax-clump [SNRCUT=%.1f(%.1f), MJDWIN=%.1f]"
	      ,INPUTS_OBS_atFLUXMAX.SNRCUT        
	      ,INPUTS_OBS_atFLUXMAX.SNRCUT_BACKUP
	      ,INPUTS_OBS_atFLUXMAX.MJDWIN	      );    
    }
    else if ( (OPTMASK & OPTMASK_SETPKMJD_TRIGGER) > 0 )  { 

    }

    printf("\n %s: %s \n", fnam, cmethod);    fflush(stdout);

  } // end VBOSE

  return ;

} // end init_obs_atFLUXMAX


void get_obs_atFLUXMAX(char *CCID, int NOBS, 
		       float *FLUX_LIST, float *FLUXERR_LIST,
		       double *MJD_LIST, int *IFILTOBS_LIST, 
		       int *OBS_atFLUXMAX) {


  // Created May 2019
  // Find observation index (o) at max flux.
  // Must call init_obs_atFLUXMAX first.
  // Always do brute-force search on first iteration.
  // If FLUXMAX2 option is set in OPTMASK (Fmax in clump),
  // do 2nd iteration  over a restricted MJDWIN region 
  // that has the most observations with SNR > SNRCUT_USER. 
  // The goal here is to reject crazy-outler observations that can
  // introduce false PEAKMJDs. If there are no obs passing
  // SNR > SNRCUT_USER, then the process is repeated with
  // an SNR-cut of SNRCUT_BACJUP ... this is to ensure that we 
  // always get a PEAKMJD estimate.
  //
  // Be careful with two different MJD windows:
  // 1) Local MJDSTEP_SNRCUT=10:
  //    NSNRCUT (number of SNR detection) is evaluated in each
  //    10 day window
  // 2) User input MJDWIN_USER (several x 10 days)
  //    This larger window is used to find max number of SNR-detections
  //    from which Fmax is evaluated. SNRCUT_SUM is the sum over
  //    NWIN_COMBINE 10-day windows.
  //
  //  While only the 2nd window (MJDWIN_USER) is used in the end, 
  //  the smaller 10-day window is to improve speed of calculation.
  //
  //  OPTMASK = 8  --> return naive Fmax over all observatins
  //  OPTMASK = 16 --> use Fmax-clump method describe above.
  //
  // Inputs:
  //   CCID         : SN id, for error message
  //   NOBS         : number of observations
  //   FLUX_LIST    : list of NOBS fluxes
  //   FLUXERR_LIST : list of NOBS flux-uncertainties
  //   MJD_LIST     : list of NOBS MJDs
  //   IFILTOBS_LIST: list of NOBS absolute filter indices
  //
  // OUTPUT:
  //   OBS_atFLUXMAX:  obs-index for max flux in each filter band
  //
  // Jul 2 2019: fix bug computing MJDMIN/MAX to work with unsorted MJD_LIST.
  //

  int    OPTMASK         = INPUTS_OBS_atFLUXMAX.OPTMASK ;
  if ( OPTMASK == 0 ) { return ; }

  double MJDWIN_USER     = INPUTS_OBS_atFLUXMAX.MJDWIN ;
  double SNRCUT_USER     = INPUTS_OBS_atFLUXMAX.SNRCUT ;
  double SNRCUT_BACKUP   = INPUTS_OBS_atFLUXMAX.SNRCUT_BACKUP;
  double MJDSTEP_SNRCUT  = 10.0 ; // hard wired param

  int USE_MJDatFLUXMAX  = (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX );
  int USE_MJDatFLUXMAX2 = (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX2);
  int USE_BACKUP_SNRCUT, ITER, NITER, IMJD, IMJDMAX=0 ;
  int NOBS_SNRCUT=0, NSNRCUT_MAXSUM=0;
  int IFILTOBS, o, omin, omax, omin2, omax2, NOTHING ;
  int MALLOC=0 ;
  double SNR, SNRCUT=0.0, SNRMAX=0.0, FLUXMAX[MXFILTINDX] ;
  double MJD, MJDMIN, MJDMAX, FLUX, FLUXERR;

  int   *NSNRCUT = NULL;      // Number of obs in each 10-day window
  int   *oMIN_SNRCUT=NULL, *oMAX_SNRCUT=NULL ;
  int    NWIN_COMBINE, MXWIN_SNRCUT, MEMI;
  int    LDMP = 0; // t(strcmp(CCID,"3530")==0 ) ;
  char fnam[] = "get_obs_atFLUXMAX" ;

  // ------------ BEGIN -------------

  NITER  = 1 ;         // always do max flux on 1st iter
  omin2  = omax2 = -9 ;
  NSNRCUT_MAXSUM = 0 ;
  USE_BACKUP_SNRCUT = 0;

  // find MJDMIN,MAX (obs may not be time-ordered)
  MJDMIN = +999999.0 ;
  MJDMAX = -999999.0 ;
  for(o=0; o < NOBS; o++ ) {
    MJD = MJD_LIST[o];
    if ( MJD < MJDMIN ) { MJDMIN = MJD; }
    if ( MJD > MJDMAX ) { MJDMAX = MJD; }
  }


  NWIN_COMBINE = (int)(MJDWIN_USER/MJDSTEP_SNRCUT + 0.01) ;
  MXWIN_SNRCUT = (int)((MJDMAX-MJDMIN)/MJDSTEP_SNRCUT)+1 ;

  if ( NWIN_COMBINE < 0 ) {  }

  if ( MXWIN_SNRCUT < 0 ) {
    sprintf(c1err,"Crazy MXWIN_SNRCUT = %d",  MXWIN_SNRCUT);
    sprintf(c2err,"MJDMIN/MAX=%.2f/%.2f  MJDSTEP_SNRCUT=%.2f  NOBS=%d",
	    MJDMIN, MJDMAX, MJDSTEP_SNRCUT, NOBS);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  MEMI  = sizeof(int) * MXWIN_SNRCUT ;

 START:

  if ( LDMP ) {
    printf(" xxx ---------  START CID=%s -------------- \n",
	   CCID ); fflush(stdout);
  }


  if ( USE_MJDatFLUXMAX ) {
    // for naive fluxmax, start with lower SNRCUT_BACKUP
    SNRCUT  = SNRCUT_BACKUP;
    USE_BACKUP_SNRCUT = 1;
  }
  else if ( USE_MJDatFLUXMAX2 ) {
    NITER   = 2 ;
    // xxx makr delete    MJDMIN  = MJD_LIST[0];
    IMJDMAX = 0;
    SNRCUT  = SNRCUT_USER;
    if ( USE_BACKUP_SNRCUT ) { SNRCUT = SNRCUT_BACKUP; }

    if ( MALLOC == 0 ) {
      NSNRCUT     = (int*)malloc(MEMI);
      oMIN_SNRCUT = (int*)malloc(MEMI);
      oMAX_SNRCUT = (int*)malloc(MEMI);
      MALLOC = 1; 
    }
    // initialize quantities in each 10-day bin
    for(o=0; o < MXWIN_SNRCUT; o++ ) {
      NSNRCUT[o]     =  0 ;
      oMIN_SNRCUT[o] = oMAX_SNRCUT[o] = -9 ;
    }

  } // end if-block over USE_MJDatFLUXMAX2


  //  - - - - - - - - - - - - -- - - - - 
  
  for(ITER=1; ITER <= NITER; ITER++ ) {

    for(IFILTOBS=0; IFILTOBS < MXFILTINDX; IFILTOBS++ ) {
      FLUXMAX[IFILTOBS]       = -9.0 ;
      OBS_atFLUXMAX[IFILTOBS] = -9 ;
    }

    if ( ITER == 1 )
      { omin = 0;  omax = NOBS-1; }
    else 
      { omin = omin2; omax=omax2; }


    if ( omin<0 || omax<0 || omin>= NOBS || omax>= NOBS ) {
      printf("\n PRE-ABORT DUMP: \n");
      printf("\t NSNRCUT_MAXSUM = %d \n", NSNRCUT_MAXSUM);
      printf("\t NSNRCUT[]      = %d, %d, %d, %d, %d ... \n",
	     NSNRCUT[0],NSNRCUT[1],NSNRCUT[2],NSNRCUT[3],NSNRCUT[4]);
      printf("\t NOBS_SNRCUT    = %d \n", NOBS_SNRCUT );
      printf("\t SNRMAX         = %.2f \n", SNRMAX );
      printf("\t IMJDMAX        = %d  \n",  IMJDMAX);      
      sprintf(c1err,"omin, omax = %d,%d   ITER=%d", omin, omax, ITER);
      sprintf(c2err,"CID=%s  NOBS=%d", CCID, NOBS);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    
    if ( LDMP ) {
      printf(" xxx ITER=%d : omin,omax=%3d-%3d   MJDWIN=%.1f-%.1f"
	     " SNRCUT=%.1f \n", 
	     ITER,omin,omax, MJDMIN, MJDMAX, SNRCUT ); 
      fflush(stdout);
    }

    NOBS_SNRCUT=0;   SNRMAX = 0.0 ;
    for(o = omin; o <= omax; o++ ) {
      MJD      = MJD_LIST[o];
      FLUX     = (double)FLUX_LIST[o];
      FLUXERR  = (double)FLUXERR_LIST[o];
      IFILTOBS = IFILTOBS_LIST[o];

      if ( IFILTOBS < 1 || IFILTOBS >= MXFILTINDX ) {
	sprintf(c1err,"Invalid IFILTOBS=%d for FLUX[%d]=%f", 
		IFILTOBS, o, FLUX );
	sprintf(c2err,"NOBS=%d", NOBS);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      IMJD = (int)MJD;    
      if ( IMJD    < 40000   ) { continue; }
      if ( FLUXERR < 1.0E-9  ) { continue ; }
      SNR = FLUX/FLUXERR ;
      
      if ( SNR > SNRMAX ) { SNRMAX = SNR; } // diagnostic 
      if ( SNR < SNRCUT ) { continue ; }

      if ( FLUX > FLUXMAX[IFILTOBS] ) { // max flux in each filter
	FLUXMAX[IFILTOBS] = FLUX ;
	OBS_atFLUXMAX[IFILTOBS] = o;
      }

      if ( FLUX > FLUXMAX[0] ) {  // global FLUXMAX
	FLUXMAX[0] = FLUX ;
	OBS_atFLUXMAX[0] = o;
      }


      // count number of SNR>cut epochs in sliding 10-day windows
      if ( USE_MJDatFLUXMAX2 && ITER==1 ) {
      
	IMJD = (int)((MJD - MJDMIN)/MJDSTEP_SNRCUT) ;
	IMJDMAX = IMJD ;
	   
	if ( IMJD < 0 || IMJD >= MXWIN_SNRCUT ) {
	  sprintf(c1err,"Invalid IMJD=%d (must be 0 to %d)", 
		  IMJD, MXWIN_SNRCUT-1);
	  sprintf(c2err,"Something is really messed up.");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
	}							  
	NSNRCUT[IMJD]++ ;
	if ( oMIN_SNRCUT[IMJD] < 0 ) { oMIN_SNRCUT[IMJD] = o; }
	if ( oMAX_SNRCUT[IMJD] < o ) { oMAX_SNRCUT[IMJD] = o; }
      } // end FmaxClump 

      NOBS_SNRCUT++ ;
      
    } // end o-loop over observations

    // =======================

    if ( LDMP ) {
      printf(" xxx ITER=%d : SNRMAX=%.2f NOBS_SNRCUT=%d \n",
	     ITER, SNRMAX, NOBS_SNRCUT ); fflush(stdout);
    }

    // if no obs pass SNRCUT on ITER=1, try again with SNRCUT_BACKUP         
    NOTHING = ( ITER==1 && NOBS_SNRCUT==0 ) ;
    if ( NOTHING ) {
      if ( USE_BACKUP_SNRCUT ) 
	{ return; }
      else
	{ USE_BACKUP_SNRCUT = 1; goto START; }
    }


    int iwin, iwin_shift, iwin_max, iwin_start;
    int oMIN_TMP, oMAX_TMP, NSNRCUT_SUM ;

    if ( USE_MJDatFLUXMAX2 && ITER==1 ) {
      //   check sliding combined windows for max NSNRCUT
      NSNRCUT_MAXSUM = 0 ;
      iwin_max = IMJDMAX - (NWIN_COMBINE-1) ;
      if ( iwin_max < 0 ) { iwin_max = 0 ; }
    
      for( iwin_start = 0; iwin_start <= iwin_max; iwin_start++ ) {

	// combine multiple 10-day windows to get a MJDWIN-day window
	NSNRCUT_SUM = 0; oMIN_TMP = oMAX_TMP = -9;

	for(iwin_shift = 0; iwin_shift<NWIN_COMBINE; iwin_shift++ ) {
	  iwin = iwin_start + iwin_shift ;
	  if ( iwin >= MXWIN_SNRCUT ) { continue; }
	  if ( NSNRCUT[iwin] > 0 ) {
	    NSNRCUT_SUM += NSNRCUT[iwin] ;
	    if ( oMIN_TMP < 0 ) { oMIN_TMP = oMIN_SNRCUT[iwin]; }
	    oMAX_TMP = oMAX_SNRCUT[iwin];
	  }
	}

	if ( NSNRCUT_SUM > NSNRCUT_MAXSUM ) {
	  NSNRCUT_MAXSUM = NSNRCUT_SUM ;
	  omin2 = oMIN_TMP ;
	  omax2 = oMAX_TMP ;
	}

      } // end loop over iwin_start
    } 

  } // end ITER loop

  if ( MALLOC ) { free(NSNRCUT);  free(oMIN_SNRCUT); free(oMAX_SNRCUT);  }

  return;

} // end get_obs_atFLUXMAX


void init_obs_atfluxmax__(int *OPTMASK, double *PARLIST, int *VBOSE)
{ init_obs_atFLUXMAX(*OPTMASK, PARLIST, *VBOSE); }

void get_obs_atfluxmax__(char *CCID, int *NOBS, float *FLUX, float *FLUXERR,
			 double *MJD, int *IFILTOBS, int *EP_atFLUXMAX) 
{
  get_obs_atFLUXMAX(CCID,*NOBS,FLUX,FLUXERR,MJD,IFILTOBS,EP_atFLUXMAX);
}


// ==============================================
int keyMatch(char *string,char *key ) {
  if ( strcmp(string,key)==0 ) 
    { return(1); }
  else
    { return(0); }
} // end keyMatch

int   uniqueMatch(char *string,char *key ) {
  // April 9 2019
  // utility to check for string match, and abort on
  // duplicate string
  char *msgSource = STRING_UNIQUE.SOURCE_of_STRING ;
  char fnam[] = "uniqueMatch";

  if ( strcmp(string,"INIT") == 0 ) {
    // interpet *key  as "source of string" to store
    printf("  Initialize %s for %s\n", fnam, key); fflush(stdout);
    sprintf(STRING_UNIQUE.SOURCE_of_STRING, "%s", key);
    STRING_UNIQUE.NLIST = 0;
    return(0);
  }

  if ( strcmp(string,key) == 0 )
    { checkStringUnique(string,msgSource,fnam);  return(1); }
  else
    { return(0); }
}  // end uniqueMatch


int uniqueOverlap (char *string,char *key ) {

  // check Overlap string match up to length of key.
  // Example: 
  //   *key    =  'file='
  //   *string = 'file=anything' will return true.
  // 
  //  The overlap match must start at beginning of *string.
  //  Note that this function uses uniqueMatch util.

  int lenkey = strlen(key);
  int match ;
  char tmpString[1000];
  char fnam[] = "uniqueOverlap" ;
  // ---------- BEGIN -----------

  if ( strcmp(string,"INIT") == 0 ) {
    // interpet *key  as "source of string" to store
    printf("  Initialize %s for %s\n", fnam, key); fflush(stdout);
    sprintf(STRING_UNIQUE.SOURCE_of_STRING, "%s", key);
    STRING_UNIQUE.NLIST = 0;
    return(0);
  }

  strncpy(tmpString,string,lenkey); tmpString[lenkey]='\0';
  match = uniqueMatch(tmpString,key);
  return(match);

} // end uniqueOverlap
 
void  checkStringUnique(char *string, char *msgSource, char *callFun) {

  // Apr 2019
  // Utility to store strings and check that each new string
  // is unique. If not unique, abort with error message.
  //
  // Inputs:
  //  *string : string to check
  //  *msgSrouce : message about source; e.g., "sim-input file"
  //  *callFun   : name of calling function
  //
  // The latter two args are used only for error message.
  //
  // string = "INIT" -> initialize arrays.

  int i, NLIST;
  char *tmpString;
  char fnam[] = "checkStringUnique" ;
  // --------------- BEGIN ------------

  if ( strcmp(string,"INIT") == 0 ) { 
    printf("  Initialize %s for %s\n", fnam, msgSource); fflush(stdout);
    STRING_UNIQUE.NLIST = 0;
    return ; 
  }

  NLIST = STRING_UNIQUE.NLIST ;
  for(i=0; i < NLIST; i++ ) {
    tmpString = STRING_UNIQUE.STRING[i];
    if ( strcmp(tmpString,string) == 0 ) {
      sprintf(c1err,"Duplicate '%s' not allowed in %s.", string, msgSource);
      sprintf(c2err,"Calling function is %s .", callFun);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  sprintf(STRING_UNIQUE.STRING[NLIST],"%s", string);
  STRING_UNIQUE.NLIST++ ;

  if ( STRING_UNIQUE.NLIST >= MXLIST_STRING_UNIQUE ) {
    sprintf(c1err,"%d unique %s strings exceeds bound.", 
	    STRING_UNIQUE.NLIST, msgSource );
    sprintf(c2err,"Check MXLIST_STRING_UNIQUE=%d", MXLIST_STRING_UNIQUE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return;

} // end checkStringUnique

// ==========================================
void init_lightCurveWidth(void) {
  // ------------ BEGIN ---------------
  // mallac arrays so that they can be realloc later as NOBS increases
  LCWIDTH.LAST_NOBS = 10 ;
  int NTMP = LCWIDTH.LAST_NOBS ;
  LCWIDTH.TLIST_SORTED    = (double*)malloc( NTMP*sizeof(double) );
  LCWIDTH.MAGLIST_SORTED  = (double*)malloc( NTMP*sizeof(double) );
  LCWIDTH.FLUXLIST_SORTED = (double*)malloc( NTMP*sizeof(double) );
  LCWIDTH.INDEX_SORT      = (int   *)malloc( NTMP*sizeof(int)    );
  return ;

} // end init_lightCurveWidth

void   init_lightcurvewidth__(void) { init_lightCurveWidth(); }

double get_lightcurvewidth__(int *OPTMASK, int *NOBS, double *TLIST,
			     double *MAGLIST, int *ERRFLAG, char *FUNCALL ) {
  double width;
  width = get_lightCurveWidth(*OPTMASK,*NOBS,TLIST,MAGLIST,
			      ERRFLAG, FUNCALL);
  return(width);
}


double get_lightCurveWidth(int OPTMASK_LCWIDTH, int NOBS, 
			   double *TLIST, double *MAGLIST,
			   int *ERRFLAG, char *FUNCALL ) {

  // Created Aug 2017
  // For input light curve, return width (days) based on OPTMASK.
  // Inputs:
  //   OPTMASK_LCWIDTH  : bit-mask of options
  //   NOBS     : number of observations
  //   TLIST    : time (days) for each Obs
  //   MAGLIST  : mag for each Obs
  //   FUNCALL  : name of calling function (for error msg only)
  //
  //   Output
  //      ERRFLAG  : user warnings, error, etc.  ERRFLAG=0 for success.
  //      Function returns width in days. 
  //             If ERRFLAG != 0, Width -> -9.
  //
  // Feb 1 2019: add OPT_RISE (4-bit)
  //

  int  OPT_FWHM     = ( (OPTMASK_LCWIDTH & 1) > 0 ) ; // FWHM 
  int  OPT_MINMAX   = ( (OPTMASK_LCWIDTH & 2) > 0 ) ; // TMAX-TMIN
  int  OPT_RISE     = ( (OPTMASK_LCWIDTH & 4) > 0 ) ; // 10-100% rise time

  int  ERRFLAG_LOCAL = 0 ;
  int  obs, isort, obsFmax, FOUND ;
  int  LDUMP_DEBUG =  0 ; // brute-force prints for debug, then abort
  int  ORDER_SORT  = +1 ; // sort with increasing epochs.

  double Width = -9.0, MAG, FLUX, ARG, TMAX, FLUXMAX, T, F, T0,T1, F0,F1; 
  double THALF_lo, THALF_hi, T_FWHM, FHALF, FMIN  ;
  char fnam[] = "get_lightCurveWidth" ;

  // ------------------- BEGIN -------------------

  // error checks
  if ( OPTMASK_LCWIDTH == 0 ) {
    sprintf(c1err,"Invalid OPTMASK=%d passed by %s", 
	    OPTMASK_LCWIDTH, FUNCALL);
    sprintf(c2err,"Check valid options in function %s", fnam);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // extend local 'sort' arrays with realloc when NOBS increases
  if ( NOBS > LCWIDTH.LAST_NOBS  ) {
    int MEMD = NOBS * sizeof(double);
    int MEMI = NOBS * sizeof(int);
    LCWIDTH.TLIST_SORTED    = (double*)realloc(LCWIDTH.TLIST_SORTED,   MEMD );
    LCWIDTH.MAGLIST_SORTED  = (double*)realloc(LCWIDTH.MAGLIST_SORTED, MEMD );
    LCWIDTH.FLUXLIST_SORTED = (double*)realloc(LCWIDTH.FLUXLIST_SORTED,MEMD );
    LCWIDTH.INDEX_SORT      = (int   *)realloc(LCWIDTH.INDEX_SORT,     MEMI );
  }

  // sort TLIST (since T=0 is usually at the end of the list)
  sortDouble(NOBS, TLIST, ORDER_SORT, LCWIDTH.INDEX_SORT );

  if ( LDUMP_DEBUG ) { printf("\n DEBUG DUMP of SORTED Tobs,mag: \n");  }

  // create local time-ordered lists for easier computation,
  // and store fluxes and epoch of peak flux
  TMAX=-999.0; FLUXMAX = 0.0 ; obsFmax=-9;
  for(obs=0; obs<NOBS; obs++ ) {
      isort = LCWIDTH.INDEX_SORT[obs];      
      T     = TLIST[isort];
      MAG   = MAGLIST[isort] ;
      LCWIDTH.TLIST_SORTED[obs]   = T; 
      LCWIDTH.MAGLIST_SORTED[obs] = MAG ;

      ARG = 0.4*(ZEROPOINT_FLUXCAL_DEFAULT - MAG);
      FLUX = pow(TEN,ARG);
      LCWIDTH.FLUXLIST_SORTED[obs] = FLUX ;
      if ( FLUX > FLUXMAX ) { TMAX=T; FLUXMAX=FLUX; obsFmax=obs; }

      if ( LDUMP_DEBUG && MAG < 40.0 ) {
	printf("\t obs=%3d isort=%3d  T=%7.2f  MAG=%.2f  F=%10.6f\n",
	       obs, isort, T, MAG, FLUX ) ;
	fflush(stdout);
      }
  }


  // ----------- start -----------

  if ( OPT_MINMAX ) {
    // Tmin - Tmax option is for illustration only
    if ( NOBS > 1 ) 
      { Width = LCWIDTH.TLIST_SORTED[NOBS-1] - LCWIDTH.TLIST_SORTED[0]; }

  }
  else if ( OPT_FWHM ) {
     
    THALF_lo=THALF_hi = -9999.0; T_FWHM = -9.0 ;  
    FHALF=FLUXMAX/2.0 ; FMIN=0.1*FLUXMAX ;

    // find half max before peak
    FOUND=0;
    for (obs=obsFmax; obs>=0; obs-- ) {
      T=LCWIDTH.TLIST_SORTED[obs] ; F=LCWIDTH.FLUXLIST_SORTED[obs] ; 
      if ( F > FHALF ) { FOUND=0; THALF_lo=-9999.0 ; }
      if ( F < FHALF && F>FMIN && FOUND==0 ) {
	FOUND=1;
	T0 = T; T1=LCWIDTH.TLIST_SORTED[obs+1] ; 
	F0 = F; F1=LCWIDTH.FLUXLIST_SORTED[obs+1] ; 
	if ( LDUMP_DEBUG ) {
	  printf("\t xxx BEFORE PEAK T0=%.2f T1=%.2f  F0=%.3f F1=%.3f \n",
		 T0, T1, F0, F1 ); 
	}
	if ( F1 != F0 ) { THALF_lo = T0 + (T1-T0)*(FHALF-F0)/(F1-F0); }
      }
    }
    // next find half max after peak ... in case of 2nd bump,
    // identify last half-max epoch.
    FOUND=0;
    for (obs=obsFmax; obs<NOBS; obs++ ) {
      T=LCWIDTH.TLIST_SORTED[obs] ; F=LCWIDTH.FLUXLIST_SORTED[obs] ;
      if ( F > FHALF ) { FOUND=0; THALF_hi=-9.0 ; }
      if ( F < FHALF && F>FMIN && FOUND==0 ) {
	FOUND=1;
	T1 = T; T0=LCWIDTH.TLIST_SORTED[obs-1] ; 
	F1 = F; F0=LCWIDTH.FLUXLIST_SORTED[obs-1] ;
	if ( LDUMP_DEBUG ) {
	  printf("\t xxx AFTER PEAK T0=%.2f T1=%.2f  F0=%.3f F1=%.3f \n",
		 T0, T1, F0, F1 ); 
	}
	if ( F1 != F0 ) { THALF_hi = T0 + (T1-T0)*(FHALF-F0)/(F1-F0); }
      }
    }

    if ( THALF_lo > -9998. && THALF_hi > 0.0 ) 
      { T_FWHM = THALF_hi - THALF_lo; }

    Width = T_FWHM ;
    
  } // end OPT_FWHM
  
  else if ( OPT_RISE ) {
    
    // use existing 'HALF' variables, even though it's 10% of maxFlux.
    THALF_lo = -9999.0 ;  
    FHALF = 0.10*FLUXMAX ; FMIN=0.02*FLUXMAX ;

    FOUND=0;
    for (obs=obsFmax; obs>=0; obs-- ) {
      T=LCWIDTH.TLIST_SORTED[obs] ; F=LCWIDTH.FLUXLIST_SORTED[obs] ; 
      if ( F > FHALF ) { FOUND=0; THALF_lo=-9999.0 ; }
      if ( F < FHALF && F>FMIN && FOUND==0 ) {
	FOUND=1;
	T0 = T; T1=LCWIDTH.TLIST_SORTED[obs+1] ; 
	F0 = F; F1=LCWIDTH.FLUXLIST_SORTED[obs+1] ; 
	if ( LDUMP_DEBUG ) {
	  printf("\t xxx BEFORE PEAK T0=%.2f T1=%.2f  F0=%.3f F1=%.3f \n",
		 T0, T1, F0, F1 ); 
	}
	if ( F1 != F0 ) { THALF_lo = T0 + (T1-T0)*(FHALF-F0)/(F1-F0); }
      }
    }

    // Width = TMAX - T(10% of peakFlux)
    if ( THALF_lo > -9998. ) { Width = (TMAX - THALF_lo); }

  } // end OPT_RISE


  *ERRFLAG = ERRFLAG_LOCAL ; // load output arg
  LCWIDTH.LAST_NOBS = NOBS ;   // update last NOBS

  if ( LDUMP_DEBUG ) {
    printf("  Inputs:  OPTMASK=%d  FUNCALL=%s \n", OPTMASK_LCWIDTH, FUNCALL);
    printf("  Output:  WIDTH=%f  ERRFLAG=%d \n", Width, ERRFLAG_LOCAL);
    if ( OPT_FWHM ) {
      printf(" TMAX=%.3f : THALF(lo,hi,sum) = %.3f, %.3f, %.3f \n",
	     TMAX, THALF_lo, THALF_hi, T_FWHM );
    }

    debugexit(fnam); // friendly abort
  }

  return(Width);

} // end get_lightCurveWidth




// ==================================================
void  update_covMatrix(char *name, int OPTMASK, int MATSIZE,
		       double (*covMat)[MATSIZE], double EIGMIN, 
		       int *istat_cov ) {


  // May 2016
  // Part of re-factor to move code out of read_data so that
  // it can be used when reading simdata_biasCor also.
  // 
  // If covMat is not invertible, fix it.
  // Note that input covMat can change !!! 

  // Inputs:
  //  *name = name of SN, used for error msg only.
  //   OPTMASK: 
  //        += 1 -> abort on bad matrix
  //        += 2 -> allow some diag cov entries to be zero
  //        += 4 -> dump info for bad cov           
  //  
  //   MATSIZE = matrix size
  //   covMat  = matrix, MATSIZE x MATSIZE
  //
  // Outputs
  // *covMat    = fixed covMatrix 
  // *istat_cov =  0 if no change; 
  //            = -1 if covMat has changed.
  //            = -9 if any covMat element is -9 ==> undefined 
  //
  // Nov 9 2018:  
  //  + replace fortran rs_ with C version of rs using eispack.c[h] 
  //  + move function from SALT2mu.c into sntools.c
  //


  int nm = MATSIZE ;
  //  int matz = 1;
  bool matz  ;
  int l, k, m, ierr, ipar, NBAD_EIGVAL ;
  int LDMP = 0; // (strcmp(name,"4249392") == 0);
  int ABORT_ON_BADCOV, ALLOW_ZERODIAG ;


  double eigval[MATSIZE];
  double eigvec[MATSIZE][MATSIZE];
  double eigvalOrig[MATSIZE], covOrig[MATSIZE][MATSIZE];
  double diag[3], EIGEN_MIN ;

  char fnam[] = "update_covMatrix" ;

  // ----------------- BEGIN -----------------

  *istat_cov = 0;

  // check for -9 entries --> undefined.
  // If found, reset cov=0 and return istat_cov=-9
  for (l=0; l<MATSIZE ; ++l)  {         
    for (k=0; k<MATSIZE; ++k)  { 
      if ( fabs(covMat[l][k]+9.0)<1.0E-6 ) {
	*istat_cov = -9 ;
	covMat[l][k] = 0.0 ;
      }
    }
  }
  if ( *istat_cov == -9 ) { return ; }

  // ------------------------
  // check OPTMASK bits
  
  ABORT_ON_BADCOV = ( (OPTMASK & 1) > 0 ) ;
  ALLOW_ZERODIAG  = ( (OPTMASK & 2) > 0 ) ;
  LDMP            = ( (OPTMASK & 4) > 0 ) ;

  // Find eigenvalues and eigenvectors, convention below
  // err[j][i]*eigvec[0][j] = eigval[0]*eigvec[0][i]
  // xxx rs_(&nm,&nm, &covMat[0][0], eigval, &matz, &eigvec[0][0], fv1,fv2, &ierr);

  if(LDMP){ printf("\t 1. xxx %s \n", fnam); fflush(stdout); }

  matz = true;
  ierr = rs(nm, &covMat[0][0], eigval, &matz, &eigvec[0][0] );

  if(LDMP){ printf("\t 2. xxx %s \n", fnam); fflush(stdout); }

  EIGEN_MIN = 0.0 ;
  if ( ALLOW_ZERODIAG ) { EIGEN_MIN = -1.0E-6 ; }

  NBAD_EIGVAL = 0 ;
  for(ipar=0; ipar < MATSIZE; ipar++ )  { 
    if (eigval[ipar] <= EIGEN_MIN ){NBAD_EIGVAL += 1;}  
  }

  // Check for illegal error matrix
  if (ierr == 0 && NBAD_EIGVAL == 0 ) { return ; }

  // check option to abort
  if ( ABORT_ON_BADCOV ) {
    printf("\n PRE-ABORT DUMP: \n");
    for(l=0; l < MATSIZE; l++ ) {
      printf("\t par%d: Eigval=%f  EigVec=%f\n", 
	     l, eigval[l], eigvec[l][l] );
    }
    sprintf(c1err,"%s matrix not positive", name);
    sprintf(c2err,"ierr(rs)=%d  NBAD_EIGVAL=%d", ierr, NBAD_EIGVAL);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  

  // ------------ FIX BAD COV -------------

  *istat_cov = -1 ;

  for (l=0;l<MATSIZE;++l)  { diag[l] = sqrt( covMat[l][l]) ; }
        
  for (k=0;k<MATSIZE;++k)  { 
    eigvalOrig[k] = eigval[k] ;
    if (eigval[k]<EIGMIN) {eigval[k]=EIGMIN;}  
  }
    
  for (l=0;l<MATSIZE;++l)  {
    for (m=0;m<MATSIZE;++m)  {
      covOrig[l][m] = covMat[l][m] ;
      covMat[l][m]=0.0;
      for (k=0;k<MATSIZE;++k) {
	covMat[l][m] += eigvec[k][l]*eigval[k]*eigvec[k][m];
      }
    }
  }


  // --------------------------------------
  // check option to dump info for cov fix
  if ( LDMP ) {
    printf("xxx -------- %s DUMP for SN=%s ------------- \n",  fnam, name);
    printf("xxx rsError=%i  Eigenval(mB,x1,c) = %f %f %f \n",
	   ierr, eigvalOrig[0], eigvalOrig[1], eigvalOrig[2] );
    fflush(stdout);

    printf("xxx  COV_ORIG(mB,x1,c)       -->      COV_FIXED(mB,x1,c) \n");
    for(l=0; l < MATSIZE; l++ ) {
      printf("xxx  ");
      for (m=0;m<MATSIZE;++m) { printf("%8.5f ", covOrig[l][m] ); }
      printf("       ");
      for (m=0;m<MATSIZE;++m) { printf("%8.5f ", covMat[l][m] ); }
      printf("\n");
    }

    printf("\n"); fflush(stdout);

  } // end LDMP


  return ;

} // end update_covMat

void update_covmatrix__(char *name, int *OPTMASK, int *MATSIZE,
			 double (*covMat)[*MATSIZE], double *EIGMIN, 
			 int *istat_cov ) {  
  update_covMatrix(name, *OPTMASK, *MATSIZE, covMat,
		   *EIGMIN, istat_cov); 
} 


int store_PARSE_WORDS(int OPT, char *FILENAME) {

  // Read FILENAME (ascii) and store each word to be fetched later
  // with next_PARSE_WORDS.
  // 
  // OPT  = -1 --> one-time init
  // OPT +=  1 --> parse file FILE
  // OPT +=  2 --> FILENAME is a string to parse, parse by space or comma
  // OPT +=  4 --> ignore comma in parsing string: space-sep only
  //                 
  // Function returns number of stored words separated by either
  // space or comma.
  //
  // Jun 26 2018: free(tmpLine) --> fix awful memory leak
  // May 09 2019: refactor to allow option for ignoring comma in strings.

  int DO_STRING    = ( (OPT & MSKOPT_PARSE_WORDS_STRING) > 0 );
  int DO_FILE      = ( (OPT & MSKOPT_PARSE_WORDS_FILE) > 0 );
  int CHECK_COMMA  = ( (OPT & MSKOPT_PARSE_WORDS_IGNORECOMMA) == 0 );
  int LENF = strlen(FILENAME);

  int NWD, MXWD, iwdStart=0, GZIPFLAG, iwd;
  char LINE[MXCHARLINE_PARSE_WORDS], *pos, sepKey[4] = " ";
  FILE *fp;
  char fnam[] = "store_PARSE_WORDS" ;
  
  // ------------- BEGIN --------------------

  if ( LENF == 0  ) { PARSE_WORDS.NWD = 0 ; return(0); }

  // if input file (or string) is the same as before,
  // just return since everything is still stored.

  if ( LENF > 0 && strcmp(PARSE_WORDS.FILENAME,FILENAME)==0 ) 
    { return(PARSE_WORDS.NWD); }

  /*
  printf(" xxx %s: OPT=%2d  BUFSIZE=%d    FILENAME='%s'\n", 
	 fnam, OPT, PARSE_WORDS.BUFSIZE, FILENAME ); fflush(stdout);
  */

  if ( OPT < 0 ) {
    PARSE_WORDS.BUFSIZE = PARSE_WORDS.NWD = 0 ;
    PARSE_WORDS.FILENAME[0] = 0 ;
    NWD = 0 ;
  }
  else if ( DO_STRING ) {
    // FILENAME is a string to parse
   
    char *tmpLine = (char*) malloc( (LENF+10)*sizeof(char) );
    sprintf(tmpLine, "%s", FILENAME); // Mar 13 2019
    if ( CHECK_COMMA && strchr(tmpLine,',') != NULL ) 
      { sprintf(sepKey,","); }

    malloc_PARSE_WORDS() ;    MXWD = PARSE_WORDS.BUFSIZE ;

    splitString2(tmpLine, sepKey, MXWD,
		 &NWD, &PARSE_WORDS.WDLIST[0] ); // <== returned
    PARSE_WORDS.NWD = NWD ;

    if ( NWD >= MXWD ) {
      sprintf(c1err,"NWD=%d exceeds bound.", NWD);
      sprintf(c2err,"Check PARSE_WORDS.BUFSIZE=%d", MXWD);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
    free(tmpLine);

  }
  else if ( DO_FILE ) {
    // read text file
    fp = open_TEXTgz(FILENAME,"rt", &GZIPFLAG );
    if ( !fp ) {
      sprintf(c1err,"Could not open text file");
      sprintf(c2err,"%s", FILENAME);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    NWD = PARSE_WORDS.NWD = 0 ;
    while( fgets(LINE, MXCHARLINE_PARSE_WORDS, fp)  != NULL ) {
      malloc_PARSE_WORDS();
      if ( (pos=strchr(LINE,'\n') ) != NULL )  { *pos = '\0' ; }
      if ( PARSE_WORDS.NWD < MXWORDFILE_PARSE_WORDS ) 
	{ iwdStart = PARSE_WORDS.NWD; }
      splitString2(LINE, sepKey, MXWORDLINE_PARSE_WORDS, 
		   &NWD, &PARSE_WORDS.WDLIST[iwdStart] ); // <== returned
      PARSE_WORDS.NWD += NWD;
    }
    NWD = PARSE_WORDS.NWD ;

    fclose(fp);
  }
  else {
    sprintf(c1err,"Invalid OPT=%d with FILENAME='%s'", OPT, FILENAME);
    sprintf(c2err,"grep MSKOPT_PARSE $SNANA_DIR/src/sntools.c");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  if ( NWD >= MXWORDFILE_PARSE_WORDS ) {
    sprintf(c1err,"NWD=%d exceeds bound, MXWORDFILE_PARSE_WORDS=%d",
	    NWD, MXWORDFILE_PARSE_WORDS );
    sprintf(c2err,"Check '%s' ", PARSE_WORDS.FILENAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // abort on tabs
  int NTAB=0;
  for(iwd=0; iwd < NWD; iwd++ ) 
    {  if ( PARSE_WORDS.WDLIST[iwd][0] == '\t' ) { NTAB++ ; }   }
  if ( NTAB > 0 ) {
    sprintf(c1err,"Found %d invalid tabs.", NTAB);
    sprintf(c2err,"Check '%s' ", FILENAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }


  if ( LENF < MXPATHLEN ) 
    { sprintf(PARSE_WORDS.FILENAME, "%s", FILENAME); }
  else
    { PARSE_WORDS.FILENAME[0] = 0 ; }


  return(NWD);;

} // end store_PARSE_WORDS

int store_parse_words__(int *OPT, char *FILENAME) 
{ return store_PARSE_WORDS(*OPT, FILENAME); }


void malloc_PARSE_WORDS(void) {
  int ADDBUF    = ADDBUF_PARSE_WORDS ;
  int MXCHARWD  = MXCHARWORD_PARSE_WORDS ;
  int NWD       = PARSE_WORDS.NWD ;
  int iwd, WD0, WD1, BUFSIZE, IFLAG=0 ;
  //  char fnam[] = "malloc_PARSE_WORDS" ;

  // ------------- BEGIN ----------------

  /*
  printf(" xxx %s: BUFSIZE = %d  (NWD=%d) \n", 
  fnam, PARSE_WORDS.BUFSIZE, NWD ); */

  if ( PARSE_WORDS.BUFSIZE == 0 ) {
    PARSE_WORDS.WDLIST    = (char**) malloc( ADDBUF*sizeof(char*) );
    PARSE_WORDS.BUFSIZE  += ADDBUF ;
    WD0 = 0;  WD1 = PARSE_WORDS.BUFSIZE ;
    IFLAG=1;
  }
  else if ( NWD > PARSE_WORDS.BUFSIZE - MXWORDLINE_PARSE_WORDS ) {
    WD0                   = PARSE_WORDS.BUFSIZE ;
    PARSE_WORDS.BUFSIZE  += ADDBUF ;
    WD1                   = PARSE_WORDS.BUFSIZE ;

    BUFSIZE               = PARSE_WORDS.BUFSIZE ;
    PARSE_WORDS.WDLIST    = 
      (char**) realloc(PARSE_WORDS.WDLIST, BUFSIZE*sizeof(char*) );
    IFLAG=2;    
  }

  else {
    IFLAG=3;
    return ;
  }

  
  for(iwd=WD0; iwd < WD1; iwd++ ) {
    PARSE_WORDS.WDLIST[iwd]  = (char*) malloc( MXCHARWD*sizeof(char) );
    PARSE_WORDS.WDLIST[iwd][0] = 0 ;
  }

  return ;

} // end malloc_PARSE_WORDS


void get_PARSE_WORD(int langFlag, int iwd, char *word) {

  // Code-language flag:
  // langFlag=0 ==> called by C code  ==> do NOT leave pad space
  // langFlag=1 ==> called by fortran ==> leave pad space

  int NWD = PARSE_WORDS.NWD ;
  char fnam[] = "get_PARSE_WORD" ;

  if ( iwd >= NWD ) {
    sprintf(c1err,"iwd=%d exceeds NWD_STORE=%d", iwd, NWD);
    sprintf(c2err,"Check '%s' ", PARSE_WORDS.FILENAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  
  // leave extra blank space so that fortran can find length
  sprintf(word, "%s ", PARSE_WORDS.WDLIST[iwd] );

  if ( langFlag==0 ) { trim_blank_spaces(word); }  // remove <CR>


  /*
  printf(" xxx %s return word[%2d] = '%s' \n", 
	 fnam, iwd, word ); 
  */
  
} // end get_PARSE_WORD

void get_parse_word__(int *langFlag, int *iwd, char *word) 
{ get_PARSE_WORD(*langFlag, *iwd, word); }


// ******************************************
void parse_multiplier(char *inString, char *key, double *multiplier) {

  // Mar 2019
  //
  // Examples:
  //  inString(I)   key(I)   multiplier(O)
  // ----------------------------------------------
  //   CC_S15        CC_S15       1.0
  //   .3*CC_S15     CC_S15       0.3
  //   CC_S15*.3     CC_S15       0.3
  //   BLABLA        CC_S15       0.0

  int    MEMC  = 40*sizeof(char);
  int    LDMP = 0;
  int    NSPLIT ;
  char  *splitValues[2], *cnum = NULL ;
  char star[] = "*" ;
  char fnam[] = "parse_multiplier" ;

  // ----------------- BEGIN ----------------

  if ( strstr(inString,key) == NULL ) { 
    *multiplier = 0.0;
  }
  else if ( strcmp(inString,key) == 0 ) {
    *multiplier = 1.0;
  }
  else {
    // do the parsing; start by splitting string around *
    splitValues[0] = (char*)malloc(MEMC);  splitValues[1] = (char*)malloc(MEMC);
    splitString(inString, star, 3,      // inputs
		&NSPLIT, splitValues );      // outputs         
    
    if ( strcmp(splitValues[0],key) == 0 ) 
      { cnum = splitValues[1]; }
    else if ( strcmp(splitValues[1],key) == 0 ) 
      { cnum = splitValues[0]; }
    else {
      sprintf(c1err,"Cannot find key='%s' after splitting string='%s'",
	      key, inString);
      sprintf(c2err,"Something is messed up.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    sscanf(cnum, "%le", multiplier ); // load output arg

    free(splitValues[0]);  free(splitValues[1]);
  }

  if ( LDMP ) {
    printf(" xxx ----------------------------- \n");
    printf(" xxx %s DUMP: \n", fnam);
    printf(" xxx inString = '%s'  key='%s' \n", inString, key);
    printf(" xxx multiplier = %le \n", *multiplier );
    debugexit(fnam);
  }

  return ;

} // end parse_multiplier

// ===========================================
void init_GENPOLY(GENPOLY_DEF *GENPOLY) {
  int o;
  GENPOLY->ORDER = -9;
  for(o=0; o < MXORDER_GENPOLY; o++ ) {
    GENPOLY->COEFF_RANGE[o][0] = 0.0 ;
    GENPOLY->COEFF_RANGE[o][1] = 0.0 ;
  }
} // end init_GENPOLY

void parse_GENPOLY(char *string, GENPOLY_DEF *GENPOLY, char *callFun) {

  // Mar 23 2019
  // parse input string and load GENPOLY structure with 
  // polynomial of arbitrary order (up to 20).
  //
  // Examples:
  //  string = "1.0,.034,1.0E-4" 
  //     --> load 2nd order polynominal
  //  string = "0.9:1.1,.034,1.0E-4,2.2E-8" 
  //     --> load 3rd order polynominal, and a0 coefficient
  //         range is 0.9 to 1.1 for random selection.
  //
  // Input *callFun is for error or debug msg only.
  //
  // Be carefuly to parse both commas and colons
  //

  int MEMC    = sizeof(char);
  int MXSPLIT = MXORDER_GENPOLY;
  int LDMP    = 0 ;

  int o, NSPLIT, NRANGE, ORDER ;
  double DVAL0, DVAL1 ;
  char *splitValue[MXORDER_GENPOLY], *splitRange[2];
  char *tmpVal, colon[] = ":", comma[] = "," ;
  char fnam[] = "parse_GENPOLY" ;

  // ----------- BEGIN ------------

  sprintf(GENPOLY->STRING,"%s", string);

  for(o=0; o < MXSPLIT; o++ ) 
    { splitValue[o] = (char*)malloc( 40*MEMC); }

  splitRange[0] = (char*)malloc( 40*MEMC);
  splitRange[1] = (char*)malloc( 40*MEMC);

  splitString(string, comma, MXSPLIT,      // inputs
              &NSPLIT, splitValue );      // outputs         

  ORDER = NSPLIT-1;
  GENPOLY->ORDER = ORDER ;

  // loop over each poly coeff and check for colon separating range
  for(o=0; o <= ORDER ; o++ ) {
    tmpVal = splitValue[o] ;
    if ( strstr(tmpVal,colon) == NULL ) {
      sscanf(tmpVal, "%le", &DVAL0 ); DVAL1=DVAL0;
    }
    else {
      splitString(tmpVal, colon, 4,
		  &NRANGE, splitRange );  // outputs
      if ( NRANGE != 2 ) {
	sprintf(c1err,"NRANGE=%d for order=%d. Expect 2 args"
		" around one colon.", 	NRANGE, o);
	sprintf(c2err,"Check %s element of %s (callFun=%s)", 
		tmpVal, string, callFun );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      sscanf(splitRange[0], "%le", &DVAL0 ); 
      sscanf(splitRange[1], "%le", &DVAL1 ); 
    }
    
    if ( DVAL1 < DVAL0 ) {
      sprintf(c1err,"Invalid range for poly-order %d (callFun=%s)", 
	      o, callFun);
      sprintf(c2err,"%s has 2nd value smaller than 1st.", tmpVal);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    GENPOLY->COEFF_RANGE[o][0] = DVAL0;
    GENPOLY->COEFF_RANGE[o][1] = DVAL1;

  } // end loop over poly terms

  for(o=0; o < MXSPLIT; o++ )  { free(splitValue[o]); }
  free(splitRange[0]);
  free(splitRange[1]);

  if ( LDMP ) {
    printf("\n xxx ======== %s DUMP =========== \n", fnam );
    printf(" xxx calling function: %s \n", callFun);
    printf(" xxx input string = '%s' \n", string);
    printf(" xxx NORDER = %d \n", ORDER );
    for(o=0; o <= ORDER; o++ ) {
      DVAL0 = GENPOLY->COEFF_RANGE[o][0];
      DVAL1 = GENPOLY->COEFF_RANGE[o][1];
      if ( DVAL1 > DVAL0 ) 
	{ printf(" xxx A%d : %le to %le \n", o, DVAL0, DVAL1); }
      else
	{ printf(" xxx A%d : %le \n", o, DVAL0 ); }
      fflush(stdout);
    }
    
    //    debugexit(fnam);
  } // end LDMP

  return ;

} // end parse_GENPOLY

double eval_GENPOLY(double VAL, GENPOLY_DEF *GENPOLY, char *callFun) {

  // Mar 2019
  // evaluate polynominal for input value 'val'.
  // Note that random numbers are used for ranges.
  // Avoid using slow 'pow' function.

  int o, ORDER = GENPOLY->ORDER;
  double VALPOLY = 0.0 ;
  double VALPOW  = 1.0 ;
  double COEFF_RANGE[2], RANCOEFF;
  char fnam[] = "eval_GENPOLY";
  int LDMP = 0 ;
  // ---------- begin ------------

  for(o=0; o <= ORDER; o++ ) {
    COEFF_RANGE[0] = GENPOLY->COEFF_RANGE[o][0];
    COEFF_RANGE[1] = GENPOLY->COEFF_RANGE[o][1];
    RANCOEFF = FlatRan ( 2, COEFF_RANGE ) ;
    // RANCOEFF = COEFF_RANGE[0];
    VALPOLY += RANCOEFF * VALPOW ;
    VALPOW = VALPOW * VAL;
  }

  if ( ORDER < 0 ) {
    sprintf(c1err,"Undefined POLYNOMIAL passed from fun=%s", callFun);
    sprintf(c2err,"VAL=%le  POLYSTRING='%s' ", VAL, GENPOLY->STRING);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( LDMP ) {
    printf(" xxx --------------------------------------------- \n");
    printf(" xxx %s DUMP for VAL=%le  and POLY='%s\n", 
	   fnam, VAL, GENPOLY->STRING );
    printf(" xxx ORDER=%d   VALPOLY = %le \n", ORDER, VALPOLY );
  }

  return(VALPOLY);

} // end eval_GENPOLY

// ====================================
const gsl_rng_type *T_Poisson ;
gsl_rng *r_Poisson ;

int getRan_Poisson(double mean){

  // Created DEc 27 2017
  // Return random Poisson integer from input mean.


  unsigned int k=0 ;

  // change random seed with gsl_rgn_set ??

  if ( mean < 1.0E-12 ) {
    gsl_rng_env_setup();
    T_Poisson = gsl_rng_default;
    r_Poisson = gsl_rng_alloc(T_Poisson);   
    return(0);
  }

  k = gsl_ran_poisson (r_Poisson,mean);
  
  //  gsl_rng_free(r);

  return( (int)k );

} // end getRan_Poisson

void get_SNANA_VERSION(char *snana_version) 
{ sprintf(snana_version, "%s", SNANA_VERSION_CURRENT); } 
void get_snana_version__(char *snana_version) 
{  get_SNANA_VERSION(snana_version); }



void INIT_SNANA_DUMP(char *STRING) {

  // Created Sep 21 2017
  // Parse input *STRING and load DUMP_STRING_INFO structure.
  // This string can be set in the analysis programs via 
  // &SNLCINP input  DUMP_STRING = 'BLA BLA'. For simulation ... (not yet).
  //
  // Example:
  // DUMP_STRING = 
  //  'FUN FLUXERRCALC  FILTERS iz  CIDLIST 5001,5002
  //      MJDRANGE 53662 53663 ABORT'
  //
  // Will trigger the dump in subroutine FLUXERRCALC for filters
  // i & z, CIDs 5001 & 5002, and within MJDRANGE given above.
  // CIDLIST is interpreated as a list of comma-separated strings.
  // The ABORT flag triggers an abort after all CIDs and FILTERS
  // have been processed.  The job will finish to completion if:
  // a) no ABORT flag, or b) CIDLIST includes an invalid CID, or
  // c) FILTERS includes an invalid filter. For completed jobs,
  // search the stdout for dump info.
  //
  // Default values:
  //   + CIDLIST   --> none
  //   + FILTERS   --> none
  //   + MJDRANGE  --> wide open (0 to 9999999)
  //
  //
  // List of functions integrated with this utility:
  //      file          function         date
  //     ----------------------------------------------
  //    snana.car     FLUXERRCALC     2017-09-21
  //        ... more later ...
  //

  int i;
  char fnam[] = "INIT_SNANA_DUMP" ;

  // --------------- BEGIN --------------

  DUMP_STRING_INFO.NCID      = 0 ;
  DUMP_STRING_INFO.NFILT     = 0 ;
  DUMP_STRING_INFO.ABORTFLAG = 0 ;
  for(i=0; i < MXCID_DUMP; i++ ) 
    { DUMP_STRING_INFO.NFILT_DONE[i] = 0 ; }

  DUMP_STRING_INFO.FUNNAME[0]  = 0 ;
  DUMP_STRING_INFO.FILTLIST[0] = 0 ;

  DUMP_STRING_INFO.MJDRANGE[0] = 0.0 ;
  DUMP_STRING_INFO.MJDRANGE[1] = 9999999. ;

  // bail if required key is not present
  if ( strstr(STRING,"FUN") == NULL ) { return ; }

  // - - - - - - - - - - - - - - - - 
#define MXWORD_DUMP_STRING 20
  char BLANK[] = " " ;
  int MXSPLIT = MXWORD_DUMP_STRING ;
  int Nsplit, iwd, NCID=0 ;
  char *ptrSplit[MXWORD_DUMP_STRING];
  char *wordList[MXWORD_DUMP_STRING];
  int  MEMC = 60 * sizeof(char);
  // allocate wordList and assign ptrSplit
  for(i=0; i < MXWORD_DUMP_STRING; i++) {
    wordList[i] = (char*) malloc ( MEMC ) ;
    ptrSplit[i] = wordList[i];
  }
  splitString(STRING, BLANK, MXSPLIT,
	      &Nsplit, ptrSplit) ;   // returned

  char *cwd0=wordList[0], *cwd1=wordList[0], *cwd2=wordList[0];
  for(iwd=0; iwd < Nsplit; iwd++ ) {
    cwd0 = wordList[iwd] ;
    if ( iwd < (Nsplit-1) ) { cwd1 = wordList[iwd+1] ; }
    if ( iwd < (Nsplit-2) ) { cwd2 = wordList[iwd+2] ; }

    if ( strcmp(cwd0,"FUN") == 0 ) 
      { sscanf(cwd1, "%s", DUMP_STRING_INFO.FUNNAME );  }

    if ( strcmp(cwd0,"FILTLIST") == 0 || strcmp(cwd0,"FILTERS")==0 )  { 
      sscanf(cwd1, "%s", DUMP_STRING_INFO.FILTLIST );  
      DUMP_STRING_INFO.NFILT = strlen(DUMP_STRING_INFO.FILTLIST);
    }

    if ( strcmp(cwd0,"MJDRANGE") ==0 ) {
      sscanf(cwd1, "%le", &DUMP_STRING_INFO.MJDRANGE[0] );  
      sscanf(cwd2, "%le", &DUMP_STRING_INFO.MJDRANGE[1] );  
    }

    if ( strcmp(cwd0,"ABORT") ==0 ) 
      { DUMP_STRING_INFO.ABORTFLAG=1; }  // no argument

    if ( strcmp(cwd0,"CID") ==0 || strcmp(cwd0,"CIDLIST")==0 ) { 
      char *ptrCCID[MXCID_DUMP],  comma[] = "," ;
      for(i=0;i<MXCID_DUMP;i++) {ptrCCID[i] = DUMP_STRING_INFO.CCIDLIST[i];}
      splitString(cwd1, comma, MXCID_DUMP, &NCID, ptrCCID) ;  
      DUMP_STRING_INFO.NCID = NCID;
    }

  } // end iwd loop

  // - - - - - -  print summary - - - - - - - - 

  char BANNER[100];
  sprintf(BANNER,"%s: dump instructions", fnam );
  print_banner(BANNER);
  printf("\t Function to Dump: '%s' \n", DUMP_STRING_INFO.FUNNAME  );
  printf("\t Filters  to Dump: '%s' \n", DUMP_STRING_INFO.FILTLIST );

  printf("\t CID list to Dump: ");
  for(i=0; i < NCID; i++ ) { printf("%s ", DUMP_STRING_INFO.CCIDLIST[i] ); }
  printf("\n");

  printf("\t Abort after Dump: %d \n",   DUMP_STRING_INFO.ABORTFLAG );
  printf("\n");

  fflush(stdout);
  //  debugexit(fnam); // xxx REMOVE

  return ;

} // end INIT_SNANA_DUMP

// ==========================================================
int CHECK_SNANA_DUMP(char *FUNNAME, char *CCID, char *BAND, double MJD) {

  int IFLAG = 0 ;
  int NCID  = DUMP_STRING_INFO.NCID ;
  int NFILT = DUMP_STRING_INFO.NFILT ;
  //  char fnam[] = "CHECK_SNANA_DUMP" ;

  // ------------ BEGIN ----------------

  if ( NCID == 0 ) { return(IFLAG); }

  // check function name
  if ( strcmp( DUMP_STRING_INFO.FUNNAME,FUNNAME) != 0 ) 
    { return(IFLAG); }

  // check valid band
  if ( strstr(DUMP_STRING_INFO.FILTLIST,BAND) == NULL ) 
    { return(IFLAG); } 

  // check valid CCID
  int i, ICID=-9;
  for(i=0; i < NCID; i++ ) {
    if ( strcmp(DUMP_STRING_INFO.CCIDLIST[i],CCID)==0 ) { ICID=i; }
  }
  if ( ICID < 0 ) { return(IFLAG); }

  // check valid MJD
  if ( MJD < DUMP_STRING_INFO.MJDRANGE[0] ) { return(IFLAG); }
  if ( MJD > DUMP_STRING_INFO.MJDRANGE[1] ) { return(IFLAG); }

  // - - - - - - - - - - -
  DUMP_STRING_INFO.NFILT_DONE[ICID]++  ;

  // if we get here, set IFLAG > 0
  IFLAG = DUMPFLAG_noABORT ; // default

  // if ABORTFLAG is set, abort only when all CID and FILTERS have
  // been dumped.
  if ( DUMP_STRING_INFO.ABORTFLAG ) {
    int NFILT_DONE, NOTDONE = 0 ;
    for(i=0; i < NCID; i++ ) {
      NFILT_DONE = DUMP_STRING_INFO.NFILT_DONE[i];
      if ( NFILT_DONE < NFILT ) { NOTDONE=1; }
    }
    if ( NOTDONE==0 ) { IFLAG = DUMPFLAG_withABORT; }
  }

  return(IFLAG);

} // end CHECK_SNANA_DUMP

void ABORT_SNANA_DUMP(void) {

  // user should call this routine if CHECK_SNANA_DUMP returns 2,
  // but call this after doing the last dump.
  // The dump here is a friendly dump, not Evil-Abort face dump.

  printf("\n Done with all user-requested DUMPs for: \n");
  printf("\t %d CIDs, FILTERS=%s, MJDRANGE=%.3f-%.3f\n",
	 DUMP_STRING_INFO.NCID, DUMP_STRING_INFO.FILTLIST,
	 DUMP_STRING_INFO.MJDRANGE[0], DUMP_STRING_INFO.MJDRANGE[1] 
	 );
  printf(" User requests ABORT after all DUMPs. Bye Bye. \n");
  fflush(stdout);
  exit(66);

} // end ABORT_SNANA_DUMP


// - - - - - - - -
void init_snana_dump__(char *STRING) { INIT_SNANA_DUMP(STRING); }

int  check_snana_dump__(char *FUNNAME, char *CCID, char *BAND, double *MJD) 
{  return ( CHECK_SNANA_DUMP(FUNNAME,CCID,BAND,*MJD) ) ; } 

void abort_snana_dump__(void) { ABORT_SNANA_DUMP(); }


// ==========================================
void read_SURVEYDEF(void) {

  // Created Aug 19 2016 [moved frmo SALT2mu.c on July 2017]
  // Read SURVEY.DEF file and store each defined survey name in
  // SURVEY_INFO.SURVEYDEF_LIST[IDSURVEY].
  //

  char SURVEYDEF_FILE[MXPATHLEN], c_get[60], nameTmp[40];
  int  idTmp ;
  FILE *fp ;
  char fnam[] = "read_SURVEYDEF" ;

  // ------------ BEGIN -------------

  sprintf(SURVEYDEF_FILE,"%s/SURVEY.DEF", PATH_SNDATA_ROOT);
  fp = fopen(SURVEYDEF_FILE,"rt");
  if ( !fp ) {
    sprintf(c1err,"Could not open SURVEY.DEF file");
    sprintf(c2err,"%s", SURVEYDEF_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // init each survey name to NULL
  for(idTmp=0; idTmp < MXIDSURVEY; idTmp++ ) { 
    sprintf(SURVEY_INFO.SURVEYDEF_LIST[idTmp],"NULL"); 
  }

  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"SURVEY:") == 0 ) {
      readchar(fp, nameTmp ); readint(fp, 1, &idTmp );
      if ( idTmp >= MXIDSURVEY ) {
	sprintf(c1err,"IDSURVEY=%d(%s) exceeds MXIDSURVEY=%d", 
		idTmp, nameTmp, MXIDSURVEY);
	sprintf(c2err,"check %s", SURVEYDEF_FILE);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      sprintf(SURVEY_INFO.SURVEYDEF_LIST[idTmp],"%s", nameTmp);
    }
  }

  fclose(fp);
  return ;

} // end read_SURVEYDEF

// ==========================================                                   
int get_IDSURVEY(char *SURVEY) {
  // return integer IDSURVEY for input string *SURVEY
  int ID;
  //  char fnam[] = "get_IDSURVEY" ;     
  // -------------- BEGIN -------------  
  for(ID=0; ID < MXIDSURVEY; ID++ ) {
    if ( strcmp(SURVEY_INFO.SURVEYDEF_LIST[ID],SURVEY) == 0 )
      { return ID; }
  }

  return(-9);
} // end get_IDSURVEY 

// ============================================
int  exec_cidmask(int mode, int CID) {

  // Created Apri 2017 (code transfered from snana.car)
  // utility to store integer CIDs in bit mask (to save memory),
  // and to check if any CID-bit is set.  Used by simulation
  // and analysis to flag duplicates.

  // Inputs:
  //   mode =  0 --> allocate memory with MXCID=CID
  //   mode =  1 --> set CID'th bit
  //   mode =  2 --> return 1 if CID'th bit is set; 0 otherwise 
  //   mode = -1 --> free memory
  //
  // Apr 25 2017: bug fix, fmodf -> fmod for double

  int J, JBIT, NINT ;
  double dCID = (double)CID;
  char fnam[] = "exec_cidmask" ;

  // --------------- BEGIN -----------------

  if ( mode == 0 ) {
    NINT = (CID+1)/32 + 2 ;
    float xMB  = (float)(NINT*4)/1.0E6 ;
    printf("   Allocate %.2f MB for CIDMASK array "
	   "(to check duplicates)\n", xMB);  fflush(stdout);
    CIDMASK_LIST = (unsigned int*) malloc ( NINT * sizeof(unsigned int) ) ;
    MXCIDMASK = CID;
    NCIDMASK_LIST = NINT;
    for(J=0; J < NINT; J++ ) { CIDMASK_LIST[J]=0; }
    return(1);
  }
  else if ( mode == 1 ) {
    if ( CID > MXCIDMASK ) { return(0); }
    J = (CID/32) + 1 ;    JBIT = (int)fmod(dCID,32.0) ;
    CIDMASK_LIST[J] |= (1 << JBIT) ;
    /*
    printf(" xxx %s: CID=%2d --> J=%d JBIT=%2d  CIDMASK=%u,%u \n",
	   fnam, CID, J, JBIT, CIDMASK_LIST[1], CIDMASK_LIST[2] ); */
    return(1);
  }
  else if ( mode == 2 ) {
    J = (CID/32) + 1 ;    JBIT = (int)fmod(dCID,32.0) ;
    if ( (CIDMASK_LIST[J] & (1 << JBIT)) > 0 )
      { return(1); }
    else
      { return(0); }
  }
  else if ( mode == -1 ) {
    free(CIDMASK_LIST);
    return(1);
  }
  else {
    sprintf(c1err,"Invalid mode=%d for CID=%d", mode, CID);
    c2err[0] = 0 ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  return(1);

} // end exec_cidmask

int  exec_cidmask__(int *mode, int *CID) {  
  int istat = exec_cidmask(*mode,*CID); 
  return(istat);
}

void test_cidmask(void) {

  // test exec_cidmask function
#define NLIST 5
  int CIDLIST[NLIST] = { 4, 22, 23, 38, 78 } ;
  int i, CID, LSET ;
  int MXCID = 10000;
  char fnam[] = "test_cidmask";

  // ---------- BEGIN ----------
  
  exec_cidmask(0,MXCID);
  for(i=0; i < NLIST; i++ ) {
    CID = CIDLIST[i];
    exec_cidmask(1,CID);
  }

  for(CID=0; CID<100 ; CID++ ) {
    LSET = exec_cidmask(2,CID);
    printf(" xxx CID = %3d --> LSET=%d \n", CID, LSET);
    
  }

  fflush(stdout);
  debugexit(fnam);

  return ;

} // end test_cidmask


void extract_MODELNAME(char *STRING, char *MODELPATH, char *MODELNAME) {

  // Created Feb 22 2016
  // For input STRING = BLA
  //    return MODELPATH='' and MODELNAME=BLA
  // For input STRING = /project/models/SALT2/SALT2.ABC
  //    return MODELPATH=STRING and MODELNAME=SALT2.ABC
  //


  int i,i2, lastSlash, LENSTR, ENVstat ;
  char fnam[] = "extract_MODELNAME" ;

  // ------------ BEGIN -----------

  MODELPATH[0] = MODELNAME[0] = 0 ;

  ENVstat   = ENVreplace(STRING,fnam,1);    // check for ENV
  LENSTR    = strlen(STRING);
  lastSlash = -9;
  for(i=0; i < LENSTR; i++ ) {
    if ( STRING[i] == '/' ) { lastSlash=i; }
  }

  if ( lastSlash < 0 ) {
    sprintf(MODELNAME, "%s", STRING);
  }
  else {
    sprintf(MODELPATH, "%s", STRING);

    for(i=lastSlash+1; i < LENSTR; i++ ) {
      i2 = i - lastSlash - 1; 
      sprintf(&MODELNAME[i2], "%c", STRING[i] ); 
    } // end i loop over string chars

  } // end lastSlah if

  return ;

} // end extract_MODELNAME

void extract_modelname__(char *STRING, char *MODELPATH, char *MODELNAME) {
  extract_MODELNAME(STRING, MODELPATH, MODELNAME);
}

void FILTER_REMAP_INIT(char *remapString, char *VALID_FILTERLIST,
		       int *NFILT_REMAP,
                       int *IFILTLIST_REMAP, char *FILTLIST_REMAP) {

  // Created Feb 2017
  // Example:
  //   Input remapString = 'abc->U def->B gh->V jklm->R'
  //   
  // Outputs
  //   NFILT_REMAP = 4
  //   IFILTLIST_REMAP = array of 4 indices corresponding to UBVR
  //   FILTLIST_REMAP  = 'UBVR'
  //
  // Inputs  VALID_FILTERLIST is a string of valid input filters
  // to check that the remapString does not contain invalid filters.
  //

  int  NFILT_ORIG, NFILT, ifilt, lenStr;
  int  USEFILT_ORIG[MXFILTINDX], USEFILT_REMAP[MXFILTINDX];
  char *localString, *ptrtok, ctmp[80] ;
  char fnam[] = "FILTER_REMAP_INIT" ;

  // ------------ BEGIN ------------

  // init output function args
  *NFILT_REMAP =  FILTLIST_REMAP[0]  = 0 ;

  // init global
  FILTER_REMAP.NMAP  = 0 ;

  lenStr = strlen(remapString);
  if ( lenStr == 0 ) { return ; }

  // - - - - -
  set_FILTERSTRING(FILTERSTRING);

  localString = (char*) malloc( lenStr* sizeof(char) );
  for(ifilt=0; ifilt<MXFILTINDX; ifilt++) { 
    USEFILT_ORIG[ifilt]  = 0; 
    USEFILT_REMAP[ifilt]  = 0; 

    FILTER_REMAP.IFILTOBS_MAP[ifilt] = -9 ;
    FILTER_REMAP.IFILT_MAP[ifilt]    = -9 ;
  }

  printf("\n %s: \n", fnam);
  printf(" remapString = \n  '%s' \n", remapString);


  // split remapString into pieces, where each piece is a map
  sprintf(localString, "%s", remapString); 
  ptrtok = strtok(localString," ") ; // split string
  NFILT  = NFILT_ORIG = 0;
  while ( ptrtok != NULL ) {
      sprintf(ctmp, "%s", ptrtok );
      if ( NFILT < MXFILT_REMAP )
	{  sprintf(FILTER_REMAP.MAPSTRING[NFILT], "%s", ctmp); }
      NFILT++ ;	
      ptrtok = strtok(NULL, " ");
  }
  free(localString) ;
  sprintf(c2err,"See remapString"); // generic abort comment.

  if ( NFILT >= MXFILT_REMAP ) {
    sprintf(c1err,"NFILT_REMAP=%d exeeds bound MXFILT_REMAP=%d.",
		NFILT, MXFILT_REMAP);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  // analyze each mapstring and load IFILTOBS_MAP
  char *ptrMap, BAND[2], band[2] ;  
  int lenMap, i, IFILTOBS, IFILT_REMAP, ifiltobs;
  for(IFILT_REMAP=0; IFILT_REMAP < NFILT; IFILT_REMAP++ ) {  
    ptrMap = FILTER_REMAP.MAPSTRING[IFILT_REMAP] ;
    lenMap = strlen(ptrMap);

    // start with last char, which is the mapped ifiltobs
    sprintf(BAND, "%c", ptrMap[lenMap-1] ) ;
    IFILTOBS     = INTFILTER(BAND);  // mapped filter
    FILTLIST_REMAP         = strcat(FILTLIST_REMAP,BAND);
    IFILTLIST_REMAP[IFILT_REMAP] = IFILTOBS ;
    // abort if remapped filter is used more than once
    if ( USEFILT_REMAP[IFILTOBS] ) {
      sprintf(c1err,"Remapped filter = %s used more than once.", BAND);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }
    USEFILT_REMAP[IFILTOBS] += 1;

    for(i=0; i < lenMap-3 ; i++ ) {
      sprintf(band, "%c", ptrMap[i] );
      ifiltobs = INTFILTER(band); // original filter

      // abort if original filter is used more than once
      if ( USEFILT_ORIG[ifiltobs] ) {
	sprintf(c1err,"original filter = %s used more than once.", band);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
      }
      // abort if original filter is not in the VALID_FILTERLIST      
      if ( strstr(VALID_FILTERLIST,band) == NULL ) {
	printf("\n PRE-ABORT DUMP:\n Defined bands: '%s'\n", 
	       VALID_FILTERLIST );
	sprintf(c1err,"Invalid filter=%s is not defined ", band);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
      }

      USEFILT_ORIG[ifiltobs] += 1;
      FILTER_REMAP.IFILTOBS_MAP[ifiltobs] = IFILTOBS;
      FILTER_REMAP.IFILT_MAP[ifiltobs]    = IFILT_REMAP ;

      NFILT_ORIG++ ;
      printf("\t %s: %s->%s (%2.2d->%2.2d)  IFILT_REMAP=%d\n",
	     ptrMap, band, BAND, ifiltobs, IFILTOBS, IFILT_REMAP);
      fflush(stdout);

    }  // i loop
  } // ifilt
 

  *NFILT_REMAP = NFILT ;
  printf("  Done mapping %d original filters to %d mapped filters.\n",
	 NFILT_ORIG, NFILT);
  fflush(stdout);

  return;

} // end FILTER_REMAP_INIT

void filter_remap_init__(char *remapString, char *VALID_FILTERLIST,
			 int *NFILT_REMAP,
                         int *IFILTLIST_REMAP, char *FILTLIST_REMAP) {
  FILTER_REMAP_INIT(remapString, VALID_FILTERLIST,
		    NFILT_REMAP, IFILTLIST_REMAP, FILTLIST_REMAP);
}

void FILTER_REMAP_FETCH(int IFILTOBS_ORIG,
			int *IFILTOBS_REMAP, int *IFILT_REMAP) {

  // Created Feb 2017
  // For input IFILTOBS_ORIG, returns 
  //   + absolute IFILTOBS_REMAP,
  //   + sparse   IFILT_REMAP indices.

  //  char fnam[] = "FILTER_REMAP_FETCH" ;

  // ------------- BEGIN -----------
  
  *IFILT_REMAP    = FILTER_REMAP.IFILT_MAP[IFILTOBS_ORIG] ;
  *IFILTOBS_REMAP = FILTER_REMAP.IFILTOBS_MAP[IFILTOBS_ORIG]  ;

  return ;

} // end FILTER_REMAP_FETCH

void filter_remap_fetch__(int *IFILTOBS_ORIG,
			  int *IFILTOBS_REMAP, int *IFILT_REMAP) {
  FILTER_REMAP_FETCH(*IFILTOBS_ORIG, IFILTOBS_REMAP, IFILT_REMAP) ;
}



// =============================================
void init_GaussIntegral(void) {

  int ix, NX;
  double xmin = 0.0;
  double xmax, xbin, xi;
  double XMAX = 10.0 ;  // store integrals up to 10 sigma
  GAUSS_INTEGRAL_STORAGE.XMAX = XMAX ;

  char fnam[] = "init_GaussIntegral";

  // --------- BEGIN -----------

  GAUSS_INTEGRAL_STORAGE.INIT_FLAG = 1;

  
  NX = 500 ;
  GAUSS_INTEGRAL_STORAGE.NBIN_XMAX = NX ;

  xbin = XMAX / (double)NX ;

  if ( NX <= 0 || NX >= MXSTORE_RAN ) {
    sprintf(c1err,"Invalid NBIN_XMAX=%d", NX);
    sprintf(c2err,"Valid range is < %d", MXSTORE_RAN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  for(ix=0; ix < NX; ix++ ) {
    xi   = (double)ix;
    xmax = xmin + xbin * (xi - 0.5) ;
    GAUSS_INTEGRAL_STORAGE.XMAX_LIST[ix] = xmax ;
    GAUSS_INTEGRAL_STORAGE.GINT_LIST[ix] = GaussIntegral(xmin,xmax);   
  }

  // done with init
  GAUSS_INTEGRAL_STORAGE.INIT_FLAG = 0;


#define NDUMP 5
  int LDMP=0, idump;
  if ( LDMP ) {
    double xmin_dump[NDUMP] = { 0.0, -1.0, -1.0, -2.0, 1.0 } ;
    double xmax_dump[NDUMP] = { 1.0,  1.0,  0.0, -1.0, 2.0 } ;
    for(idump=0 ; idump < NDUMP ; idump++ ) {
      xmin= xmin_dump[idump];  xmax=xmax_dump[idump];
      printf(" xxx DUMP %s(%4.1f,%4.1f) = %f \n",  
	     fnam, xmin, xmax, GaussIntegral(xmin,xmax) );
    }
    debugexit(fnam);
  }

  return ;

} // end init_GaussIntegral

// =============================================
double GaussIntegral(double nsig1, double nsig2) {

  // Created Oct 2016
  // return Gaussian integral between nsig1 and nsig2
  //
  // Mar 19 2019: if nsig[1,2] > 10, set GINT=0 to avoid abort.

  int    NBIN, ibin ;
  double SUM, DENOM, SIGBIN, xsig, xi, arg ;
  char fnam[] = "GaussIntegral" ;

  // --------- BEGIN ----------

  SUM    = 0.0 ;

  if ( nsig1 == nsig2 ) { return SUM; }

  if ( GAUSS_INTEGRAL_STORAGE.INIT_FLAG  ) {
    DENOM = sqrt(TWOPI);
    // brute-force integral during init stage
    NBIN  = 1000*(int)(nsig2-nsig1);
    if ( NBIN < 5 ) { NBIN=5; }
    SIGBIN = (nsig2-nsig1)/(double)NBIN;
    for(ibin=0; ibin < NBIN; ibin++ ) {
      xi     = (double)ibin ;
      xsig   = nsig1 + SIGBIN*(xi+0.5);
      arg    = (xsig*xsig)/2.0 ;
      SUM   += exp(-arg);
    }
    SUM *= SIGBIN/DENOM ;
  }

  else {
    // interpolate stored integrals from 0 to xmax
    int    OPT_INTERP=1;
    double GINT1, GINT2, x1, x2, xsign1=1.0, xsign2=1.0 ;
    double *ptr_xmax = GAUSS_INTEGRAL_STORAGE.XMAX_LIST ;
    double *ptr_gint = GAUSS_INTEGRAL_STORAGE.GINT_LIST ;
    double  NSIGMAX  = 0.995*GAUSS_INTEGRAL_STORAGE.XMAX ;

    x1 = fabs(nsig1);  if(nsig1!=0.0) { xsign1 = x1/nsig1; }
    x2 = fabs(nsig2);  if(nsig2!=0.0) { xsign2 = x2/nsig2; }
    NBIN = GAUSS_INTEGRAL_STORAGE.NBIN_XMAX ;

    if ( NBIN <= 0 || NBIN >= MXSTORE_RAN ) {
      sprintf(c1err,"Invalid NBIN_XMAX=%d", NBIN);
      sprintf(c2err,"Check if init_GaussIntegral was called");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }


    if (x1 > NSIGMAX ) { x1 = NSIGMAX; }
    if (x2 > NSIGMAX ) { x2 = NSIGMAX; }

    GINT1 = interp_1DFUN(OPT_INTERP,x1,NBIN, ptr_xmax, ptr_gint, fnam);
    GINT2 = interp_1DFUN(OPT_INTERP,x2,NBIN, ptr_xmax, ptr_gint, fnam);

    SUM   = xsign2*GINT2 - xsign1*GINT1;
  }

  /*
    if ( xsign1 > 0.0 && xsign2 > 0.0 ) {
    printf(" xxx SUM2/SUM = %f/%f = %f \n",
    SUM2, SUM, SUM2/SUM ); 
    debugexit(fnam); // xxxx
    }*/

  return(SUM) ;

} // end GaussIntegral

// =============================================
double angSep( double RA1,double DEC1, 
	       double RA2,double DEC2, double  scale) {

  // Copied from DIFFIMG on Nov 16 2015
  //
  // Oct 9, 2012 R. Kessler
  // for input coords of point 1 (RA1,DEC1) and point 2 (RA2,DEC2),
  // return angular separation. Inputs are in degrees and output
  // is in degrees x scale ->
  // * scale = 1    -> output is in degrees
  // * scale = 60   -> output is in arcmin
  // * scale = 3600 -> output is in arcsec

  double X1,Y1,Z1, X2, Y2, Z2, DOTPROD, sep ;
  double RAD = RADIAN ;

  // ------------- BEGIN ------------------

  X1 = cos(RA1*RAD) * cos(DEC1*RAD);
  Y1 = sin(RA1*RAD) * cos(DEC1*RAD);
  Z1 = sin(DEC1*RAD);

  X2 = cos(RA2*RAD) * cos(DEC2*RAD);
  Y2 = sin(RA2*RAD) * cos(DEC2*RAD);
  Z2 = sin(DEC2*RAD);

  DOTPROD = (1.0-1.0E-15)*(X1*X2 + Y1*Y2 + Z1*Z2);
  sep = acos(DOTPROD)/RAD ; // angular sep, degrees

  return (sep * scale) ;

} // end of angSep


int ENVreplace(char *fileName, char *callFun, int ABORTFLAG) {

  // Feb 2015 [major overhaul Mar 30 2019]
  // if input *fileName starts with $XXX/yyy 
  // then replace $XXX with getenv("XXX").
  // Note that input fileName is modified.
  //
  // Inputs:
  //   *fileName   : file to check for ENV
  //   *callFun    : name of calling function to print on error
  //    ABORTFLAG  : non-zero -> abort on error; else return null fileName
  //
  //  Return SUCCESS or ERROR
  //

  int LL, i, FOUNDSLASH, SEV, NFILE; ;
  char firstChar[2], c1[2], ENVname[MXPATHLEN], suffix[MXPATHLEN] ;
  char fnam[] = "ENVreplace" ;  
  
  // ------------- BEGIN -------------

  if ( strcmp(fileName,"init") == 0 || strcmp(fileName,"INIT") ==0 ) 
    { ENVreplace_store.NFILE = 0 ; return(SUCCESS); }

  sprintf(firstChar,"%c", fileName[0] );
  suffix[0]=0;
  ENVname[0]=0;

  NFILE = ENVreplace_store.NFILE;
  if ( NFILE < MXFILE_ENVreplace-1 ) 
    { sprintf(ENVreplace_store.FILENAME_ORIG[NFILE],"%s", fileName); }

  if ( *firstChar == '$' ) { 

    FOUNDSLASH = 0 ;
    LL = strlen(fileName);
    for(i=1; i < LL ; i++ ) {
      sprintf(c1,"%c", fileName[i] );
      if( *c1 == '/' ) {  FOUNDSLASH++ ;  }

      if ( FOUNDSLASH > 0  ) 
	{ strcat(suffix,c1); }
      else
	{ sprintf(ENVname, "%s%s", ENVname, c1) ; }
    }

    
    /*
    printf(" xxx --------------------------------------------- \n");
    printf(" xxx fileName = '%s' \n", fileName);
    printf(" xxx ENVname  = '%s' \n", ENVname );
    printf(" xxx suffix   = '%s' \n", suffix );
    */

    if ( getenv(ENVname) == NULL ) {
      if ( ABORTFLAG ) { SEV = SEV_FATAL; } else { SEV = SEV_WARN; }
      sprintf(c1err,"getenv(%s) failed (callFun=%s)", ENVname, callFun);
      sprintf(c2err,"check fileName='%s'", fileName);
      errmsg(SEV, 0, fnam, c1err, c2err);     
      return(ERROR);
    }

    sprintf(fileName, "%s%s", getenv(ENVname), suffix);
  } 

  if ( NFILE < MXFILE_ENVreplace-1 ) 
    { sprintf(ENVreplace_store.FILENAME_ENVreplace[NFILE],"%s", fileName); }
  ENVreplace_store.NFILE++ ;

  return(SUCCESS) ; 

} // end of ENVreplace

// void envreplace_(char *fileName) { ENVreplace(fileName); }

void ENVrestore(char *fileName_ENVreplace, char *fileName_orig) {

  // Created July 21, 2019
  // return fileName_orig before ENVreplace was called.

  int NFILE = ENVreplace_store.NFILE ;
  int ifile ;
  char *FILENAME_ENVreplace, *FILENAME_ORIG ;
  // ------------- BEGIN -----------

  sprintf(fileName_orig, "%s", fileName_ENVreplace); 

  for(ifile=0; ifile < NFILE; ifile++ ) {
    FILENAME_ENVreplace = ENVreplace_store.FILENAME_ENVreplace[ifile] ;
    FILENAME_ORIG       = ENVreplace_store.FILENAME_ORIG[ifile];
    if ( strcmp(fileName_ENVreplace,FILENAME_ENVreplace) == 0 ) {
      sprintf(fileName_orig,"%s", FILENAME_ORIG );
    }
  }

  return;
} // end ENVrestore


int intrac_() {  return ((int) isatty(0)); }  // needed by fortran minuit (from intrac.c)

void warn_oldInputs(char *varName_old, char *varName_new) {

  int NWARN, i;

  NWARN = OLD_INPUTS.NWARN ;
  if ( strcmp(varName_old,"init") == 0 ) 
    { OLD_INPUTS.NWARN = 0 ; return ; }

  if ( strcmp(varName_old,"list") == 0 ) {
    if ( NWARN == 0 ) { return ; }

    return ; // remove this when new system goes live
    printf("\n");
    for(i=0; i < NWARN; i++ ) {
      printf("  WARNING: REPLACE OLD INPUT  %s  with %s \n"
	     ,OLD_INPUTS.VARNAME_OLD[i]
	     ,OLD_INPUTS.VARNAME_NEW[i] );
      fflush(stdout);
    }
    return ;
  } // end list loop


  sprintf(OLD_INPUTS.VARNAME_OLD[NWARN], "%s", varName_old );
  sprintf(OLD_INPUTS.VARNAME_NEW[NWARN], "%s", varName_new );
  OLD_INPUTS.NWARN++ ;

} // end of warn_oldInputs

void warn_oldinputs__(char *varName_old, char* varName_new) 
{ warn_oldInputs(varName_old, varName_new); }

// ===================================
void set_FILTERSTRING(char *FILTERSTRING) {
  // Feb 2013: return list with every filter defined.
  sprintf(FILTERSTRING,"%s", FILTERSTRING_DEFAULT );
}


void set_EXIT_ERRCODE(int ERRCODE) {   EXIT_ERRCODE = ERRCODE; }
void set_exit_errcode__(int *ERRCODE) { set_EXIT_ERRCODE(*ERRCODE); }

// ====================================================
int IGNOREFILE(char *fileName) {

  // May 2014
  // Return 1 if fileName is "NULL" or "NONE".
  // These are valid names to override an existing fileName 
  // with a command to ignore it.
  // Return 0 otherwise.
  //
  // Dec 18 2017: include ' '

  if ( strcmp_ignoreCase(fileName,"null")   == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"none")   == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"blank")  == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"NULL")   == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"NONE")   == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"BLANK")  == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName," ")      == 0 ) { return 1 ; }
  if ( strlen(fileName) == 0 ) { return 1; } // Sep 2016

  return  0 ;

}  // end of IGNOREFILE

int ignorefile_(char *fileName) { return IGNOREFILE(fileName); }

// ==================================================
int strcmp_ignoreCase(char *str1, char *str2) {

  // Feb 2014: analog of strcmp, but convert to lower case
  //           before testing so that result is case-insensitive.
  //
  // return strcmp on lower-case comparison.

  int len1, len2, j ;
  char str1_lc[100], str2_lc[100];

  // --------- BEGIN ----------
  len1 = strlen(str1) ;
  len2 = strlen(str2) ;
  if ( len1 != len2 ) { return -1 ; }

  for(j=0; j<len1; j++ ) { str1_lc[j] = tolower ( str1[j] ) ;  }
  for(j=0; j<len2; j++ ) { str2_lc[j] = tolower ( str2[j] ) ;  }

  str1_lc[len1] = '\0' ;
  str2_lc[len2] = '\0' ;

  return strcmp(str1_lc,str2_lc) ;
  
} // end of strcmp_ignoreCase

// ==================================================
void invertMatrix(int N, int n, double *Matrix ) {

  // Jan 2013
  // *Matrix is a 2D matrix of dimension NxN
  // but only the n x n subset is to be inverted.
  //
  // input *Matrix is overwritten with its inverse;
  // if you want to save the original matrix, save
  // it before calling this function.

  int s;
  int i1, i2, J ;
  
  // Define all the used matrices
  gsl_matrix * m         = gsl_matrix_alloc (n, n);
  gsl_matrix * inverse   = gsl_matrix_alloc (n, n);
  gsl_permutation * perm = gsl_permutation_alloc (n);

  // Fill the matrix m

  J = 0;
  for ( i1=0; i1 < N; i1++ ) {
    for ( i2=0; i2 < N; i2++ ) {
      
      if ( i1 < n && i2 < n )  { 
	gsl_matrix_set(m,i1,i2, Matrix[J] );
      }

      J++ ;
    }
  }

  // Make LU decomposition of matrix m
  gsl_linalg_LU_decomp (m, perm, &s);

  // Invert the matrix m
  gsl_linalg_LU_invert (m, perm, inverse);

  // load inverse into Matrix
  J = 0;
  for ( i1=0; i1 < N; i1++ ) {
    for ( i2=0; i2 < N; i2++ ) {
      
      if ( i1 < n && i2 < n )  { 
	Matrix[J] = gsl_matrix_get(inverse,i1,i2) ;
      }

      J++ ;
    }
  }

  gsl_matrix_free(m) ;
  gsl_matrix_free(inverse) ;
  gsl_permutation_free(perm) ;

}  // end of invertMatrix

void invertmatrix_(int *N, int *n, double *Matrix ) {
  invertMatrix(*N, *n, Matrix ) ;
}


void randominit_(int *ISEED) {  srandom(*ISEED) ; } 


void sortDouble(int NSORT, double *ARRAY, int ORDER, 
		int *INDEX_SORT) {

  // Created Jan 12 2013 by R.Kessler
  // Wrapper to sort NSORT elements of ARRARY (replaces CERNLIB's SORTZV)
  // ORDER = -1/+1  -> decreasing/increasing order
  // Output is *INDEX_SORT (input ARRAY is NOT changed)

  size_t stride_t = 1;
  size_t n_t, *index_t ;
  int i;

  // ---------------- BEGIN ----------------

  n_t = NSORT ;
  index_t = (size_t*)malloc(n_t * sizeof(size_t) );

  // native function is double
  gsl_sort_index(index_t, ARRAY, stride_t, n_t);

  // load output arg.
  for(i=0; i < NSORT; i++ ) { INDEX_SORT[i] = index_t[i]; }

  // check to flip from increasing to decreasing order.
  if ( ORDER < 0 ) 
    { reverse_INDEX_SORT(NSORT, INDEX_SORT); }
  
  free(index_t);

} // end of sortDouble

// mangled sort function for fortran
void sortDouble_(int *NSORT, double *ARRAY, int *ORDER, 
		 int *INDEX_SORT) {
  sortDouble(*NSORT, ARRAY, *ORDER, INDEX_SORT) ;

  int i;
  // add 1 to INDEX_SORT so that index starts at 1 instead of 0
  for ( i=0; i < *NSORT; i++ ) { INDEX_SORT[i] += 1; }
}


void sortFloat(int NSORT, float *ARRAY, int ORDER, 
	       int *INDEX_SORT) {

  // Created Jan 12 2013 by R.Kessler
  // Wrapper to sort NSORT elements of ARRARY (replaces CERNLIB's SORTZV)
  // ORDER = -1/+1  -> decreasing/increasing order
  // Output is *INDEX_SORT (input ARRAY is NOT changed)

  size_t stride_t = 1;
  size_t n_t, *index_t ;
  int i;

  // ---------------- BEGIN ----------------

  n_t = NSORT ;
  index_t = (size_t*)malloc(n_t * sizeof(size_t) );

  gsl_sort_float_index(index_t, ARRAY, stride_t, n_t );

  for(i=0; i < NSORT; i++ ) { INDEX_SORT[i] = index_t[i]; }
  
  // check to flip from increasing to decreasing order.
  if ( ORDER < 0 )  { 
    reverse_INDEX_SORT(NSORT, INDEX_SORT); 
  }
  
  free(index_t);

} // end of sortFloat

// mangled sort function for fortran
void sortfloat_(int *NSORT, float *ARRAY, int *ORDER, 
		int *INDEX_SORT) {
  sortFloat(*NSORT, ARRAY, *ORDER, INDEX_SORT) ;

  int i;
  // add 1 to INDEX_SORT so that index starts at 1 instead of 0
  for ( i=0; i < *NSORT; i++ ) { INDEX_SORT[i] += 1; }
}



void sortInt(int NSORT, int *ARRAY, int ORDER, int *INDEX_SORT) {

  // Created Jan 12 2013 by R.Kessler
  // Wrapper to sort NSORT elements of ARRARY (replaces CERNLIB's SORTZV)
  // ORDER = -1/+1  -> decreasing/increasing order
  // Output is *INDEX_SORT (input ARRAY is NOT changed)

  size_t stride_t = 1;
  size_t n_t, *index_t ;
  int i;

  // ---------------- BEGIN ----------------

  n_t = NSORT ;
  index_t = (size_t*)malloc(n_t * sizeof(size_t) );

  gsl_sort_int_index(index_t, ARRAY, stride_t, n_t);
  
  for(i=0; i < NSORT; i++ ) { INDEX_SORT[i] = index_t[i]; }

  // check to flip from increasing to decreasing order.
  if ( ORDER < 0 ) 
    { reverse_INDEX_SORT(NSORT, INDEX_SORT); }
  

  free(index_t);

} // end of sortInt

// mangled sort function for fortran
void sortint_(int *NSORT, int *ARRAY, int *ORDER, 
	      int *INDEX_SORT) {
  sortInt(*NSORT, ARRAY, *ORDER, INDEX_SORT) ;

  int i;
  // add 1 to INDEX_SORT so that index starts at 1 instead of 0
  for ( i=0; i < *NSORT; i++ ) { INDEX_SORT[i] += 1; }
}


void reverse_INDEX_SORT(int NSORT, int *INDEX_SORT) {

  // revserse the order of INDEX_SORT.
  // Note that input INDEX_SORT array is changed.

  int *INDEX_TMP, i ;


  INDEX_TMP = (int*)malloc( NSORT * sizeof(int) );
    
  for ( i=0; i < NSORT; i++ ) 
    {  INDEX_TMP[i] = INDEX_SORT[i];  }

  for ( i=0; i < NSORT; i++ )  
    {  INDEX_SORT[i] = INDEX_TMP[NSORT-i-1]; }


  free(INDEX_TMP);
    
}  // end of revserse_INDEX_SORT

// ======================================================
void print_KEYwarning(int ISEV, char *key_old, char *key_new) {

  // print warning message to use the new key instead of the old key
  // ISEV = SEV_WARN  -> give warning but do not abort.
  // ISEV = SEV_FATAL -> abort

  char fnam[] = "print_KEYwarning";
  char line[] = "%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%" ;

  if ( ISEV == SEV_FATAL ) {
    sprintf(c1err,"key = '%s' is no longer valid.", key_old);
    sprintf(c2err,"Must use new key = '%s'", key_new);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }
  else {
    printf("\n");
    printf("%s%s\n", line, line);
    printf("\n");
    printf(" WARNING: use new '%s' key instead of old '%s' key' \n", 
	   key_new, key_old);
    printf("          Old '%s' key will soon become obsolete.\n", 
	   key_old);
    printf("\n");
    printf("%s%s\n", line, line);
    printf("\n");
  }

} // end of print_KEYwarning


// ***************************************************
double SNR_calculator(double ZPT, double PSF, double SKYMAG, double MAG,
		      double *FLUX_and_ERR ) {

  // Created May 2014 by R. Kessler
  // Return signal-to-noise (SNR) for input MAG and 
  // input observing conditions,
  //   ZPT:    Npe = 10**[-0.4*(mag-ZPT)]
  //   PSF:    FWHM, arcsec
  //   SKYMAG: mag/arcsec^2
  //
  // Jun 1 2014: return FLUX and its error in FLUX_and_ERR array.

  double SQSNR, SNR, tA, F, FSKY, arg ;
  double OMEGA = 1.51; // Pararam to define effective sky-noise area
                       //  A = ( OMEGA * PSF )^2
  // ------------- BEGIN -----------

  tA    = ( OMEGA * PSF ) * ( OMEGA * PSF );  // area
  arg   = -0.4*(MAG    - ZPT);   F    = pow(10.0,arg);
  arg   = -0.4*(SKYMAG - ZPT);   FSKY = pow(10.0,arg);
  SQSNR = (F*F)/( F + tA*FSKY );
  SNR   = sqrt(SQSNR);

  // load flux and its error to allow external code to get coadded SNR
  FLUX_and_ERR[0] = F ;
  FLUX_and_ERR[1] = sqrt( F + tA*FSKY );

  return SNR ;

} // end of SNRMAG

// ************************************************************
double MAGLIMIT_calculator(double ZPT, double PSF, double SKYMAG, double SNR){

  // October 2010 by P. Belov & S.Glazov
  //
  // compute and return point-source limiting mag corresponding to SNR = S/N
  // ZPT:    Npe = 10**[-0.4*(mag-ZPT)]
  // PSF:    FWHM, arcsec
  // SKYMAG: mag/arcsec^2
  // SNR   : S/N for limiting mag
  //
  // Aug 7, 2013 RK - fix bugs computing MLIM
  // May 17 2014 RK - moved from simulation code

  double MLIM;
  double tMLIM;  // Temp value for MLIM
  double EA    = 1.00; // The fraction of source-flux
  double OMEGA = 1.51; // Pararam to define effective sky-noise area
                       //  A = ( OMEGA * PSF )^2
  double eps; // epsilon
  double acc; // accuracy
  double arg;

  int cnt, cnt_orig; // counter in order to avoid infinite loop
  double tA, t1, t2; // Temp variables

  char fnam[] = "MAGLIMIT" ;

  // -------------- BEGIN --------------

  // check for crazy args
  if ( ZPT <= 1.01 || PSF <= 0.001 || SKYMAG <= 1.0 ) {
    MLIM = 0.0 ;
    return MLIM;
  }

  eps = 10000.0  ;
  acc =     0.001;  // required accuracy for convergence
  cnt = cnt_orig =  100 ;  
  t2  =  5.0 * log10(SNR/EA);
  tA  = ( OMEGA * PSF ) * ( OMEGA * PSF );  // area
  t1  = 2.5 * log10( tA );
  MLIM = 20.0 ; // very rough guess.

  // here we get the value for mlimit by iterative procedure
  while ( ( eps > acc ) && ( cnt > 0 ) ) {

    if ( cnt == cnt_orig ) 
      // initial guess is with no signal-noise
      {  tMLIM  = 0.5 * ( ZPT + SKYMAG - t1 - t2 ); }
    else
      {  tMLIM = MLIM ; } // previous iteration

    arg    = 0.4*(SKYMAG-tMLIM) ;  // RK added this Feb 2011


    MLIM = 0.5*(tMLIM + ZPT - t2) - 1.25*log10(1.0 +  tA*pow(10.0,-arg) );
    eps    = fabs( MLIM - tMLIM );

    /*
    printf(" xxx cnt=%3d: MLIM=%f  tMLIM=%f  arg=%f \n",
    cnt, MLIM, tMLIM, arg); */
    cnt   -= 1 ;
  }

  if ( cnt == 0 ) {    
    sprintf(c1err,"MLIM=%.3f value has not converged: eps=%f", 
	    MLIM, eps );
    sprintf(c2err,"ZP=%.3f  PSF=%.3f  SKYMAG=%.2f SNR=%.1f ",
	    ZPT, PSF, SKYMAG, SNR );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  /*
  double SNR_CHECK = SNR_calculator(ZPT,PSF,SKYMAG, MLIM);
  printf(" xxxx SNR_CHECK/SNR = %.4f/%.4f = %f \n",
	SNR_CHECK, SNR, SNR_CHECK/SNR ); fflush(stdout);
  */


  return MLIM;

} // end of SIMLIB_maglimit


void set_SNDATA(char *key, int NVAL, char *stringVal, double *parVal ) {

  // Created May 2012
  // Load SNDATA structure so that snfitsio can be called by
  // fortran (snana.exe) to write data in fits format.
  //
  // Key is the variable name  such as REDSHIFT of SNID.
  // *stringVal is set if the variable is a string;
  // *parVal is set for double, float or int, but note that
  // *parVal is always passed as double.
  // NVAL is the number of values passed.
  //
  // Jun 15, 2012: add HOSTGAL info
  // May 26, 2013: add HOSTGLA_OBJID (long long)
  // Feb     2014: add LOGMASS[_ERR]
  // Aug  6, 2014: add SIM_MAGOBS
  // Aug  7, 2014: add NXPIX, NYPIX, XPIX, YPIX
  // Sep 03, 2014: check HOSTGAL_MAG, and HOSTGAL_SB (filter-dependent)
  // Jul 11, 2015: add REDSHIFT_SN[_ERR]
  //
  // Oct 02, 2015: localString[2000] -> *localString + malloc to avoid
  //                array bound problems. See new logical USE_stringVal.
  //
  // Feb 17 2017: add SUBSURVEY
  // Mar 29 2017: allow FLT or BAND

  int i, iep, ifilt_obs, NVAR, USE_stringVal, LEN_stringVal ;
  int LOAD_TEL, LOAD_FLT, LOAD_FIELD;
  char  *ptrtok, ctmp[20], *localString ;
  char  fnam[] =  "set_SNDATA" ;

  // ------------ BEGIN ------------

  // check if using stringVal
  LOAD_TEL   = ( strcmp(key,"TELESCOPE") == 0  );

  LOAD_FLT   = ( (strcmp(key,"FLT")   == 0)  ||
		 (strcmp(key,"BAND")  == 0)  ) ;

  LOAD_FIELD = ( strcmp(key,"FIELD")     == 0  );
  USE_stringVal = ( LOAD_FLT || LOAD_FIELD || LOAD_TEL ) ;
		   
  // ----------------------------------------
  // global header info
  if ( strcmp(key,"SURVEY") == 0 ) 
    {  sprintf(SNDATA.SURVEY_NAME, "%s", stringVal);  }

  else if ( strcmp(key,"SUBSURVEY") == 0 ) 
    {  sprintf(SNDATA.SUBSURVEY_NAME, "%s", stringVal);  }

  else if ( strcmp(key,"FILTERS") == 0 ) {  
    sprintf(SNDATA_FILTER.LIST, "%s", stringVal);  

    set_FILTERSTRING(FILTERSTRING);
    SNDATA_FILTER.NDEF = 
      PARSE_FILTLIST(SNDATA_FILTER.LIST, SNDATA_FILTER.MAP );
    
  }
  else if ( strcmp(key,"SNANA_DIR") == 0 ) 
    {  sprintf(PATH_SNANA_DIR, "%s", stringVal);  }

  else if ( strcmp(key,"PRIVATE_KEYWORD") == 0 )  { 
    SNDATA.NVAR_PRIVATE++ ;
    NVAR = SNDATA.NVAR_PRIVATE ;
    if ( NVAR >= MXVAR_PRIVATE ) {
      sprintf(c1err,"NVAR_PRIVATE exceeds bound of MXVAR_PRIVATE=%d",
	      MXVAR_PRIVATE);
      sprintf(c2err,"See MXVAR_PRIVATE in sndata.h");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    sprintf(SNDATA.PRIVATE_KEYWORD[NVAR],"%s", stringVal );
  }

  else if ( strcmp(key,"PRIVATE_VALUE") == 0 )  { 
    for(i=1; i <= NVAL ; i++ ) 
      {  SNDATA.PRIVATE_VALUE[i] = parVal[i-1] ; }
  }

  // header info for each data file
  else if ( strcmp(key,"SNID") == 0 ) 
    {  sprintf(SNDATA.CCID, "%s", stringVal);  }
  else if ( strcmp(key,"IAUC") == 0 ) 
    {  sprintf(SNDATA.IAUC_NAME, "%s", stringVal);  }
  
  else if ( strcmp(key,"RA") == 0 )
    {  SNDATA.RA = parVal[0] ;  }
  else if ( strcmp(key,"DEC") == 0 )
    {  SNDATA.DEC = parVal[0] ;  }
  else if ( strcmp(key,"DECL") == 0 )
    {  SNDATA.DEC = parVal[0] ;  }

  else if ( strcmp(key,"SNTYPE") == 0 )
    {  SNDATA.SNTYPE = (int)parVal[0] ;  }

  else if ( strcmp(key,"NEPOCH") == 0 ) {  
    SNDATA.NEPOCH = (int)parVal[0] ;  
    SNDATA.NOBS   = (int)parVal[0] ;  
  }

  else if ( strcmp(key,"SEARCH_TYPE") == 0 )  // SDSS only
    {  SNDATA.SEARCH_TYPE = (int)parVal[0] ;  }

  else if ( strcmp(key,"PEAKMJD") == 0 )
    {  SNDATA.SEARCH_PEAKMJD = parVal[0] ;  }

  else if ( strcmp(key,"PIXSIZE") == 0 )
    {  SNDATA.PIXSIZE = parVal[0] ;  }

  else if ( strcmp(key,"CCDNUM") == 0 )
    {  SNDATA.CCDNUM[0] = (int)parVal[0] ;  }

  else if ( strcmp(key,"NXPIX") == 0 )
    {  SNDATA.NXPIX = (int)parVal[0] ; }

  else if ( strcmp(key,"NYPIX") == 0 )
    {  SNDATA.NYPIX = (int)parVal[0] ;  }

  else if ( strcmp(key,"MASK_FLUXCOR_SNANA") == 0 )
    {  SNDATA.MASK_FLUXCOR = (int)parVal[0] ;  }

  else if ( strcmp(key,"REDSHIFT_FINAL") == 0 )
    {  SNDATA.REDSHIFT_FINAL = parVal[0] ;  }
  else if ( strcmp(key,"REDSHIFT_FINAL_ERR") == 0 )
    {  SNDATA.REDSHIFT_FINAL_ERR = parVal[0] ;  }

  else if ( strcmp(key,"REDSHIFT_HELIO") == 0 )
    {  SNDATA.REDSHIFT_HELIO = parVal[0] ;  }
  else if ( strcmp(key,"REDSHIFT_HELIO_ERR") == 0 )
    {  SNDATA.REDSHIFT_HELIO_ERR = parVal[0] ;  }

  else if ( strcmp(key,"VPEC") == 0 )
    {  SNDATA.VPEC = parVal[0] ;  }
  else if ( strcmp(key,"VPEC_ERR") == 0 )
    {  SNDATA.VPEC_ERR = parVal[0] ;  }


  // - - -  host - - - 

  else if ( strcmp(key,"HOSTGAL_NMATCH") == 0 )
    {  SNDATA.HOSTGAL_NMATCH[0] = (int)parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL_NMATCH2") == 0 )
    {  SNDATA.HOSTGAL_NMATCH[1] = (int)parVal[0] ;  }

  else if ( strcmp(key,"HOSTGAL_OBJID") == 0 )
    {  SNDATA.HOSTGAL_OBJID[0] = (long long)parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_OBJID") == 0 )
    {  SNDATA.HOSTGAL_OBJID[1] = (long long)parVal[0] ;  }

  else if ( strcmp(key,"HOSTGAL_PHOTOZ") == 0 )
    {  SNDATA.HOSTGAL_PHOTOZ[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL_PHOTOZ_ERR") == 0 )
    {  SNDATA.HOSTGAL_PHOTOZ_ERR[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_PHOTOZ") == 0 )
    {  SNDATA.HOSTGAL_PHOTOZ[1] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_PHOTOZ_ERR") == 0 )
    {  SNDATA.HOSTGAL_PHOTOZ_ERR[1] = parVal[0] ;  }

  else if ( strcmp(key,"HOSTGAL_SPECZ") == 0 )
    {  SNDATA.HOSTGAL_SPECZ[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL_SPECZ_ERR") == 0 )
    {  SNDATA.HOSTGAL_SPECZ_ERR[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_SPECZ") == 0 )
    {  SNDATA.HOSTGAL_SPECZ[1] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_SPECZ_ERR") == 0 )
    {  SNDATA.HOSTGAL_SPECZ_ERR[1] = parVal[0] ;  }

  else if ( strcmp(key,"HOSTGAL_RA") == 0 )
    {  SNDATA.HOSTGAL_RA[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_RA") == 0 )
    {  SNDATA.HOSTGAL_RA[1] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL_DEC") == 0 )
    {  SNDATA.HOSTGAL_DEC[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_DEC") == 0 )
    {  SNDATA.HOSTGAL_DEC[1] = parVal[0] ;  }

  else if ( strcmp(key,"HOSTGAL_SNSEP") == 0 )
    {  SNDATA.HOSTGAL_SNSEP[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_SNSEP") == 0 )
    {  SNDATA.HOSTGAL_SNSEP[1] = parVal[0] ;  }

  else if ( strcmp(key,"HOSTGAL_DDLR") == 0 )
    {  SNDATA.HOSTGAL_DDLR[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_DDLR") == 0 )
    {  SNDATA.HOSTGAL_DDLR[1] = parVal[0] ;  }

  else if ( strcmp(key,"HOSTGAL_LOGMASS") == 0 )
    {  SNDATA.HOSTGAL_LOGMASS[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL_LOGMASS_ERR") == 0 )
    {  SNDATA.HOSTGAL_LOGMASS_ERR[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_LOGMASS") == 0 )
    {  SNDATA.HOSTGAL_LOGMASS[1] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_LOGMASS_ERR") == 0 )
    {  SNDATA.HOSTGAL_LOGMASS_ERR[1] = parVal[0] ;  }


  else if ( strcmp(key,"HOSTGAL_sSFR") == 0 )
    {  SNDATA.HOSTGAL_sSFR[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL_sSFR_ERR") == 0 )
    {  SNDATA.HOSTGAL_sSFR_ERR[0] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_sSFR") == 0 )
    {  SNDATA.HOSTGAL_sSFR[1] = parVal[0] ;  }
  else if ( strcmp(key,"HOSTGAL2_sSFR_ERR") == 0 )
    {  SNDATA.HOSTGAL_sSFR_ERR[1] = parVal[0] ;  }

  // - - -  filter-dependent HOST properties - - - - 
  else if ( strcmp(key,"HOSTGAL_MAG") == 0 ) {
    SNDATA.HOSTGAL_USEMASK |= 1 ;
    for(i=0; i<NVAL; i++) { SNDATA.HOSTGAL_MAG[0][i]=(float)parVal[i];}
  }
  else if ( strcmp(key,"HOSTGAL2_MAG") == 0 ) {
    for(i=0; i<NVAL; i++) { SNDATA.HOSTGAL_MAG[1][i]=(float)parVal[i];}
  }
  else if ( strcmp(key,"HOSTGAL_MAGERR") == 0 ) {
    for(i=0; i<NVAL; i++) { SNDATA.HOSTGAL_MAGERR[0][i]=(float)parVal[i];}
  }
  else if ( strcmp(key,"HOSTGAL2_MAGERR") == 0 ) {
    for(i=0; i<NVAL; i++) { SNDATA.HOSTGAL_MAGERR[1][i]=(float)parVal[i];}
  }

  // match-independent host properties 
  else if ( strcmp(key,"HOSTGAL_CONFUSION") == 0 )
    { SNDATA.HOSTGAL_CONFUSION = (float)parVal[0] ; }

  else if ( strcmp(key,"HOSTGAL_SB_FLUXCAL") == 0 ) {
    SNDATA.HOSTGAL_USEMASK |= 4 ;
    for(i=0; i < NVAL ; i++ ) 
      {  SNDATA.HOSTGAL_SB_FLUX[i] = (float)parVal[i] ;  }
  }
  // - - - -

  else if ( strcmp(key,"MWEBV") == 0 )
    {  SNDATA.MWEBV = parVal[0] ;  }
  else if ( strcmp(key,"MWEBV_ERR") == 0 )
    {  SNDATA.MWEBV_ERR = parVal[0] ;  }

  // - - - -  EPOCH-DEPENDENT variables - - - - 

  else if ( strcmp(key,"MJD") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { 
      SNDATA.MJD[i] = parVal[i-1] ; 
      SNDATA.USE_EPOCH[i] = 1;
    }
  }
  else if ( strcmp(key,"FLUXCAL") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.FLUXCAL[i] = parVal[i-1] ; }
  }
  else if ( strcmp(key,"FLUXCALERR")     == 0 ||
	    strcmp(key,"FLUXCAL_ERRTOT") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.FLUXCAL_ERRTOT[i] = parVal[i-1] ; }
  }

  else if ( strcmp(key,"MAG") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.MAG[i] = parVal[i-1] ; }
  }
  else if ( strcmp(key,"MAGERR") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { 
      SNDATA.MAG_ERRPLUS[i]  = parVal[i-1] ; 
      SNDATA.MAG_ERRMINUS[i] = parVal[i-1] ; 
    }
  }

  else if ( strcmp(key,"PHOTFLAG") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.PHOTFLAG[i] = (int)parVal[i-1] ; }
  }

  else if ( strcmp(key,"PHOTPROB") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.PHOTPROB[i] = parVal[i-1] ; }
  }

  else if ( strcmp(key,"PSF_SIG1") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.PSF_SIG1[i] = parVal[i-1] ; }
  }
  else if ( strcmp(key,"PSF_SIG2") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.PSF_SIG2[i] = parVal[i-1] ; }
  }
  else if ( strcmp(key,"PSF_RATIO") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.PSF_RATIO[i] = parVal[i-1] ; }
  }

  else if ( strcmp(key,"SKY_SIG") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.SKY_SIG[i] = parVal[i-1] ; }
  }
  else if ( strcmp(key,"SKY_SIG_T") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.SKY_SIG_T[i] = parVal[i-1] ; }
  }

  else if ( strcmp(key,"RDNOISE") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.READNOISE[i] = parVal[i-1] ; }
  }

  else if ( strcmp(key,"GAIN") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.GAIN[i] = parVal[i-1] ; }
  }

  else if ( strcmp(key,"ZEROPT") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.ZEROPT[i] = parVal[i-1] ; }
  }
  else if ( strcmp(key,"ZEROPT_ERR") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.ZEROPT_ERR[i] = parVal[i-1] ; }
  }
  else if ( strcmp(key,"XPIX") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.XPIX[i] = parVal[i-1] ; }
  }
  else if ( strcmp(key,"YPIX") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.YPIX[i] = parVal[i-1] ; }
  }

  else if ( strcmp(key,"SIM_MAGOBS") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) { SNDATA.SIMEPOCH_MAG[i] = parVal[i-1] ; }
  }

  else if ( strcmp(key,"SIM_FLUXCAL_HOSTERR") == 0 ) {
    for(i=1; i <= NVAL ; i++ ) 
      { SNDATA.SIMEPOCH_FLUXCAL_HOSTERR[i] = parVal[i-1] ; }
  }

  else if ( USE_stringVal ) {

    LEN_stringVal = strlen(stringVal);
    localString   = (char*) malloc ( 10 + LEN_stringVal * sizeof(char) );

    /*
    printf(" xxx LEN(%s) = %d \n", stringVal, LEN_stringVal ); 
    fflush(stdout); */

    sprintf(localString, "%s", stringVal);
    ptrtok = strtok(localString," ") ; // split string
    iep=0;
    while ( ptrtok != NULL ) {
      sprintf(ctmp, "%s ", ptrtok );
      iep++ ;	
      if ( LOAD_TEL ) {
	sprintf(SNDATA.TELESCOPE[iep], "%s", ctmp );
      }
      else if ( LOAD_FLT ) {
	ifilt_obs = INTFILTER(ctmp);
	sprintf(SNDATA.FILTCHAR[iep], "%s", ctmp);
      }
      else if ( LOAD_FIELD ) {
	sprintf(SNDATA.FIELDNAME[iep], "%s", ctmp );
      }
      else {
	sprintf(c1err,"USE_stringVal=TRUE for key='%s'", key);
	sprintf(c2err,"but cannot find array to load");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      ptrtok = strtok(NULL, " ");
    }

    free(localString) ;
  }
  else {
    sprintf(c1err,"Unknown key = '%s' ", key);
    sprintf(c2err,"stringVal='%s'  parVal=%f", stringVal, parVal[0] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return ;

} // end of set_SNDATA


// mangled function for fortran
void set_sndata__(char *key, int *NVAL, char *stringVal, double *parVal ) {
  set_SNDATA(key, *NVAL, stringVal, parVal);
}

// ************************************************************

double NoiseEquivAperture(double PSFSIG1, double PSFSIG2, 
			  double PSFratio) {

/* ===========================================================
   Returns effective sky aperture, or noise-equivalent area
   
   = 1/ integral[PSF^2 rdr dtheta]
   
   for double-gaussian PSF where

   SIG1  : sigma of inner gaussian       (psf_sigma1 in psField)
   SIG2  : sigma of outer gaussian       (psf_sigma2 in psField)
   ratio : ratio of gaussians at origin  (pdf_b      in psField)

   PSF(r) is defined locally as
   
              exp(-r^2/2\sig1^2)            exp(-r^2/2\sig2^2)
   PSF = A1 * -------------------   +  A2 * ------------------
               2 * PI * sig1^2               2 * PI * sig1^2 
 
  where A1 + A2 = 1.00 for normalized PSF.
  With this definition,
 
            A2/sig2^2                    1
   ratio = -----------  ,  A1 = ------------------------
            A1/sig1^2            1 + ratio*sig2^2/sig1^2
 
 
  and effective aperture  = 
 
            4 * PI * ( sig1^2 + sig2^2 )
    =  ------------------------------------------
         1 + [ A1 * sig2/sig1 + A2*sig1/sig2]^2
 
  Special cases/checks:
   * sig2 = A2 = 0 : Aperture -> 4 * PI * sig1^2
   * sig1 = sig2   : Aperture -> 4 * PI * sig1^2 
   * sig2 = 2*sig1 : Aperture -> 4 * PI * sig1^2 / ( 1 - A2*6/5 )


  Note:  moved here from snana.car (SKY_APERTURE) on Feb 27, 2012,
         and REAL*4 -> double.
 
   ===================== */

  // local var


  //  double PI = 3.14159265 ;
  double PI = TWOPI/2.0 ;
  double A, A1, A2, SQSIG1, SQSIG2, SQSUM, TMP ;

  // ---------- BEGIN -----------

  A = 0.0 ;

  SQSIG1 = PSFSIG1 * PSFSIG1 ;
  SQSIG2 = PSFSIG2 * PSFSIG2 ;
  SQSUM  = SQSIG1 + SQSIG2 ;

  // check for single-gaussian PSF

  if ( PSFratio < 1.0E-5 || PSFSIG2 < 0.0001 ) {
    A = 4.0 * PI * SQSIG1 ;
    return A ;
  }

  TMP  = PSFratio * SQSIG2 / SQSIG1  ;
  A1   = 1.0 / ( 1.0 + TMP ) ;
  A2   = 1.0 - A1 ;
  TMP  = A1 * PSFSIG2/PSFSIG1  +  A2 * PSFSIG1/PSFSIG2 ;
  A    = 4.0 * PI * SQSUM / ( 1.0 + TMP*TMP ) ;
  return A ;

}  // end of   NoiseEquivAperture 

double noiseequivaperture_(double *PSFSIG1, double *PSFSIG2, 
			   double *PSFratio){
  return NoiseEquivAperture(*PSFSIG1, *PSFSIG2, *PSFratio) ;
}

// =================================================
double PROB_Chi2Ndof(double chi2, int Ndof ) {

  // Nov 17, 2011, R.Biswas
  // Return probability that chi2 > chi2 with Ndof
  // Note that this function should replaces CERNLIB's PROB function.
  //
  // Feb 3, 2012: protect againsta Ndof < 1 -> returns P=0
  //

  double N, chi2red, P;
  N = (double)Ndof/2.0 ;
  chi2red = chi2/2.0 ;

  if ( Ndof < 1 ) 
    { P = 0.0 ; }   
  else if ( chi2 < 0.0 ) 
    { P = 1.0 ; }
  else
    { P = gsl_sf_gamma_inc_Q(N, chi2red); }

  return P ;

} // end of PROB_Chi2Ndof

double prob_chi2ndof__(double *chi2, int *Ndof) {
  return PROB_Chi2Ndof(*chi2, *Ndof);
} 

// =================================================
double get_SIMEFFMAP(int OPTMASK, int NVAR, double *GRIDVALS) {

  // Created July 2011 by R.Kessler
  // return EFFICIENCY from interpolating map.
  // Must call init_SIMEFFMAP() before calling this function.
  // GRIDVALS are the actaul values for each variable
  // on the grid.
  // OPTMASK bit 1 (LSB) => return EFF/EFFMAX (default returns EFF)

  int  ivar, istat, iflag;
  int LDMP = 0 ;

  double 
    EFF
    ,tmpval
    ,TMPVAL
    ,TMPVAL_LIST[MXGENVAR_SIMEFFMAP]
    ,MINVAL_4LOG10 = 1.0E-12
    ,arg
    ;
  char fnam[] =  "get_SIMEFFMAP";

  //  --------------- BEGIN -------------

  if ( NVAR != SIMEFFMAP.NGENVAR || SIMEFFMAP.NGENVAR == 0 ) {
    sprintf(c1err,"Passed NVAR=%d", NVAR);
    sprintf(c2err,"but NGENVAR = %d", SIMEFFMAP.NGENVAR );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // check for LIN, LOG and INV scales

  if ( LDMP ) {
    printf(" 0. xxxxxxxx --------------------------------- \n");
  }

  TMPVAL = 0.0 ;
  for ( ivar=1; ivar <= SIMEFFMAP.NGENVAR; ivar++ ) {
    iflag  = SIMEFFMAP.IFLAGSCALE[ivar];
    tmpval = GRIDVALS[ivar-1];

    if ( iflag == 1 )        // LINEAR
      { TMPVAL = tmpval ; }  
    else if ( iflag == 10 ) {     // LOG10
      arg = tmpval ;
      if ( tmpval < MINVAL_4LOG10 )  { arg = MINVAL_4LOG10 ; }
      TMPVAL = log10(arg) ; 
    }  
    else if ( iflag == -1 )     // INVERSE
      { TMPVAL = 1./tmpval ; }  
    else {
      sprintf(c1err,"Invalid IFLAGSCALE = %d for %s = %f",
	      iflag, SIMEFFMAP.VARNAME[ivar], tmpval );
      sprintf(c2err,"%s", "Something is really screwed up.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // if tmpval is  outside the min/max bound,
    // then pull it to the edge.

    if ( TMPVAL < SIMEFFMAP.VARMIN[ivar] ) {
      TMPVAL = SIMEFFMAP.VARMIN[ivar] ; 
    }

    if ( TMPVAL > SIMEFFMAP.VARMAX[ivar] )  { 
      TMPVAL = SIMEFFMAP.VARMAX[ivar] ;  //  RANGE * 1.0E-6 ; 
    }

    TMPVAL_LIST[ivar-1] = TMPVAL ;

    if ( LDMP ) {
      printf(" 1. xxxxxxx ivar=%d : iflag=%d  VAL(orig,interp) = %le ,%le \n",
	     ivar, iflag, tmpval, TMPVAL );
      fflush(stdout); 
    }

  }

  // do the interpolation
  
  istat = interp_GRIDMAP(&SIMEFF_GRIDMAP, TMPVAL_LIST, &EFF ); // return EFF

  if ( istat != SUCCESS ) {
    printf("\n PRE-ABORT DUMP \n" );
    
    for ( ivar=1; ivar <= SIMEFFMAP.NGENVAR; ivar++ ) {
      printf("\t %-12s = %f\n", 
	     SIMEFFMAP.VARNAME[ivar], *(GRIDVALS+ivar-1) );
    }

    sprintf(c1err,"interp_GRIDMAP returned istat=%d", istat);
    sprintf(c2err,"%s", "--> could not interplate");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }



  if ( OPTMASK & 1 ) { EFF /= SIMEFFMAP.EFFMAX ; }


  return EFF;

} // end of get_SIMEFFMAP


// =================================================
int init_SIMEFFMAP(char *file, char *varnamesList) {

  // Created July 2011 by R.Kessler
  // initialize effic. map stored in input *file.
  // Return list of blank-seprated variable names in *varnamesList.
  // For example, if the variables are redshift (Z), extinction (AV)
  // and mlcs DELTA, then  varnamesList = "Z AV DELTA" .
  //

  int NFUN=1;
  int NGENVAR, NVAR_READ, NEFF_READ, NBINTOT, ivar, gzipFlag ;
  double TMPVAL[MXGENVAR_SIMEFFMAP];
  double EFF;

  FILE *fp;
  char 
    c_get[200]
    ,PATH_SIMEFF[MXPATHLEN]
    ,FILE[MXPATHLEN]
    ,*cptr
    ,fnam[] = "init_SIMEFFMAP" 
    ;

  // -------------- BEGIN --------------

  varnamesList[0] = 0 ;      // init output
  SIMEFFMAP.NGENVAR = 0 ;       // init global 
  SIMEFFMAP.EFFMAX  = 0.0 ;

  // construct official PATH_SIMEFF
  sprintf(PATH_SIMEFF,"%s/models/simeff", getenv("SNDATA_ROOT") );

  // check user's directory, then check PATH_SIMEFF
  fp = snana_openTextFile(1,PATH_SIMEFF, file, FILE, &gzipFlag );
  
  // abort if there is no SIMEFF file.
  if ( fp == NULL ) {
    sprintf(c1err,"%s", "Could not open SIMEFF file");
    sprintf(c2err,"%s", file );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  sprintf(BANNER,"%s", "init_SIMEFF: read efficiency map from ");
  print_banner(BANNER);
  printf("\t %s\n", FILE );

  NVAR_READ =  NEFF_READ =  NGENVAR   = 0 ;
  NBINTOT   = 1 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"NGENVAR:") == 0 )  {
      readint(fp, 1, &NGENVAR);
      SIMEFFMAP.NGENVAR = NGENVAR ;
      if ( NGENVAR >= MXGENVAR_SIMEFFMAP ) {
	sprintf(c1err,"NGENVAR=%d exceeds bound of", NGENVAR);
	sprintf(c2err,"MXGENVAR_SIMEFFMAP = %d", MXGENVAR_SIMEFFMAP );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    } // NGENVAR: key


    if ( strcmp(c_get,"GENVAR:") == 0 )  {
      NVAR_READ++ ;
      readchar(fp, SIMEFFMAP.VARSCALE[NVAR_READ] );
      readchar(fp, SIMEFFMAP.VARNAME[NVAR_READ] );
      readint (fp, 1, &SIMEFFMAP.NBIN[NVAR_READ] );
      readdouble(fp, 1, &SIMEFFMAP.VARMIN[NVAR_READ] );
      readdouble(fp, 1, &SIMEFFMAP.VARMAX[NVAR_READ] );

      NBINTOT *= SIMEFFMAP.NBIN[NVAR_READ] ;
      SIMEFFMAP.NBINTOT = NBINTOT;
      sprintf(varnamesList,"%s %s ", 
	      varnamesList, SIMEFFMAP.VARNAME[NVAR_READ] );

      cptr = SIMEFFMAP.VARSCALE[NVAR_READ];
      if ( strcmp(cptr,"LIN") == 0 ) 
	{ SIMEFFMAP.IFLAGSCALE[NVAR_READ] = 1; }
      else if ( strcmp(cptr,"LOG") == 0 ) 
	{ SIMEFFMAP.IFLAGSCALE[NVAR_READ] = 10; } 
      else if ( strcmp(cptr,"INV") == 0 ) 
	{ SIMEFFMAP.IFLAGSCALE[NVAR_READ] = -1; } 
      else {
	sprintf(c1err,"Invalid scale = '%s'", cptr);
	sprintf(c2err,"%s", "valid scales are LIN, LOG and  INV");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }	       

      printf("\t Found MAP variable %s %-8s : %2d bins from %7.3f to %7.3f\n"
	     ,SIMEFFMAP.VARSCALE[NVAR_READ]
	     ,SIMEFFMAP.VARNAME[NVAR_READ]
	     ,SIMEFFMAP.NBIN[NVAR_READ]
	     ,SIMEFFMAP.VARMIN[NVAR_READ]
	     ,SIMEFFMAP.VARMAX[NVAR_READ]
	     );

      // allocate temp memory when all of the GENVAR keys are read
      if ( NVAR_READ == NGENVAR ) { malloc_SIMEFFMAP(+1); }

    } // GENVAR: key


    if ( strcmp(c_get,"EFF:") == 0 )  {

      NEFF_READ++ ;
      readdouble(fp, NGENVAR+1, TMPVAL ) ;

      for ( ivar=0; ivar < NGENVAR; ivar++ ) {
	SIMEFFMAP.TMPVAL[ivar][NEFF_READ-1] = TMPVAL[ivar];
      }

      EFF = TMPVAL[NGENVAR];

      SIMEFFMAP.TMPEFF[NEFF_READ-1] = EFF ;
      if ( EFF > SIMEFFMAP.EFFMAX )
	{ SIMEFFMAP.EFFMAX  = EFF ; }

    }


  } // fscanf


  //  printf("\t EFF(MAX) = %6.4f \n", SIMEFFMAP.EFFMAX);

  // init multi-dimensional interpolation


  init_interp_GRIDMAP(IDGRIDMAP_SIMEFFMAP, "SIMEFF",
		      NBINTOT, NGENVAR, NFUN, 0,
		      SIMEFFMAP.TMPVAL, &SIMEFFMAP.TMPEFF,
		      &SIMEFF_GRIDMAP ); // <== return this struct

  // free temp memory
  malloc_SIMEFFMAP(-1);


  // make sure we read what was expected.
  if ( NVAR_READ != NGENVAR ) {
    sprintf(c1err, "Expected NGENVAR = %d", NGENVAR );
    sprintf(c2err, "but found %d  'GENVAR:' keys.", NVAR_READ);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( NEFF_READ != NBINTOT ) {
    sprintf(c1err, "Expected %d 'EFF:' keys", NBINTOT );
    sprintf(c2err, "but found %d  keys.", NEFF_READ);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return NGENVAR ;

} // end of init_SIMEFFMAP

// mangled routines for fortrans
int init_simeffmap__(char *file, char *varnamesList) {
  return init_SIMEFFMAP(file, varnamesList);
}
double get_simeffmap__(int *OPTMASK, int *NVAR, double *GRIDVALS) {
  return get_SIMEFFMAP(*OPTMASK, *NVAR, GRIDVALS);
}

// ===========================================
void malloc_SIMEFFMAP(int flag) {

  // Created July 2011 by R.Kessler
  // flag > 0 => allocate
  // flag < 0 => free 

  int  NVAR, NBINTOT, MEM, ivar, I8, I8p, I4 ;
  double XMEM;
  //  char fnam[] = "malloc_SIMEFFMAP" ;

  // --------------- BEGIN ---------------

  NVAR    = SIMEFFMAP.NGENVAR; 
  NBINTOT = SIMEFFMAP.NBINTOT ;

  MEM = (NVAR+1)*(NBINTOT) * 8;
  XMEM = 1.0E-6 * (double)MEM ;

  I8p = sizeof(double*);
  I8  = sizeof(double );
  I4  = sizeof(int);

  if ( flag > 0 ) {
    printf("\t Allocate %6.3f MB of temp SIMEFFMAP memory (NBIN=%d). \n", 
	   XMEM, NBINTOT );

    SIMEFFMAP.TMPEFF = (double  *)malloc(NBINTOT*I8); // EFF memory

    SIMEFFMAP.TMPVAL = (double **)malloc(NVAR*I8p);   // memory for each var
    for ( ivar=0; ivar < NVAR; ivar++ ) {
      SIMEFFMAP.TMPVAL[ivar] = (double *)malloc(NBINTOT*I8);
    }


  }
  else if ( flag < 0 ) {
    printf("\t Free temp SIMEFFMAP memory. \n" );
    free(SIMEFFMAP.TMPEFF);
    for ( ivar=0; ivar < NVAR; ivar++ ) {
      free(SIMEFFMAP.TMPVAL[ivar]);
    }
    free(SIMEFFMAP.TMPVAL) ;
  }

} // end of malloc_SIMEFFMAP


// ====================================
int getInfo_PHOTOMETRY_VERSION(char *VERSION      // (I) photometry version 
			       ,char *DATADIR     // (I/O) dir with data files 
			       ,char *LISTFILE    // (O) name of list file
			       ,char *READMEFILE  // (O) name of readme file
			       ) {

  // Created Mar 4, 2011 by R.Kessler
  // For input photometry VERSION, returns 
  // DATADIR:     SNANA data directory 
  // LISTFILE:    file with list of data files in DATADIR
  // READMEFILE:  file describing data sample
  //
  // If DATADIR is non-null, then it is passed from PRIVATE_DATA_PATH,
  // so use this directory to  determine LISTFILE and READMEFILE.
  // Do NOT change DATADIR in this case. If DATADIR=='', then set
  // to default $SNDATA_ROOT/lcmerge.
  //
  // Function return argument is SUCCESS or ERROR
  // 
  // Dec 2 2012: rename function, 
  //   getInfo_SNANA_VERSION -> getInfo_PHOTOMETRY_VERSION
  //   to avoid confusion with snana version (i.e, v10_19)
  //
  // Nov 11 2014: allow data to be in a folder under lcmerge/[VERSION]
  // Feb 10, 2015: same fix for DATADIR and DATADIR/[VERSION]
  // Feb 26, 2015: for datadir, also check SNDATA_ROOT/SIM
  //
  // Nov 18 2017: 
  //  + check for user-define sim-paths to be compatible with
  //    sim-input option PATH_SNDATA_SIM.. See PATHLIST below.
  //    
  // Sep 19 2018: 
  //  remove local & obsolete MXDIR_CHECK and instead use global
  //  MXPATH_SNDATA_SIM. Abort if NDIR_CHECK >= MXPATH_SNDATA_SIM .
  //
  // Sep 12 2019: 
  //  + abort if DATADIR corresponds to $SNDATA_ROOT/SIM or any SIM path

  char 
    SNDATA_ENV[20] = "SNDATA_ROOT"
    ,SNDATA_ROOT[MXPATHLEN]
    ,tmpDir[MXPATH_SNDATA_SIM][MXPATHLEN]
    ,tmpFile[MXPATH_SNDATA_SIM][MXPATHLEN]
    ,fnam[] = "getInfo_PHOTOMETRY_VERSION"
    ;

  int idir, ifound, NFOUND, NDIR_CHECK;
  int idirFOUND[MXPATH_SNDATA_SIM];
  int LDMP = 0;
  FILE *fp ;

  // ---------- BEGIN -----------

  // init outputs to NULLSTRING value

  sprintf(LISTFILE,   "%s", NULLSTRING );
  sprintf(READMEFILE, "%s", NULLSTRING );

  // first make sure that required env is set.
  if ( getenv(SNDATA_ENV) == NULL ) {
    sprintf(c1err,"env variable '$%s' is not set.",  SNDATA_ENV );
    sprintf(c2err,"%s", "      ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  sprintf(SNDATA_ROOT, "%s", getenv(SNDATA_ENV) ) ;


  // define list of directories to check for data
  idir=0;

  if ( strlen(DATADIR) > 0 ) { 
    // private user dir
    sprintf(tmpDir[idir], "%s" ,          DATADIR ); 
    sprintf(tmpFile[idir],"%s/%s.LIST",   tmpDir[idir], VERSION  );
    idir++ ;

    sprintf(tmpDir[idir], "%s/%s",        DATADIR, VERSION);
    sprintf(tmpFile[idir],"%s/%s.LIST",   tmpDir[idir], VERSION  );
    idir++ ;
    
  }
  else {

    // default locations under SNDATA_ROOT
    sprintf(tmpDir[idir], "%s/lcmerge" ,  SNDATA_ROOT );
    sprintf(tmpFile[idir],"%s/%s.LIST",   tmpDir[idir], VERSION  );    
    idir++ ;

    sprintf(tmpDir[idir], "%s/lcmerge/%s",  SNDATA_ROOT, VERSION);
    sprintf(tmpFile[idir],"%s/%s.LIST",     tmpDir[idir], VERSION  );
    idir++ ;
  }

  /* xxxxx mark delete Sep 12 2019 (see below) xxxxxxxxxxx
  // always tack on SIM dir to check
  sprintf(tmpDir[idir],  "%s/SIM/%s" , SNDATA_ROOT,  VERSION );
  sprintf(tmpFile[idir], "%s/%s.LIST", tmpDir[idir], VERSION  );
  idir++ ;
  xxxxxxxxxxxxx end mark xxxxxxxxxxxxxx  */

  // - - - - - - - - - - - - - - - - - - - - - 
  // Nov 18 2017: check for user-defined SIM-output dirs
  int ipath, NPATH;
  char **PATHLIST;

  PATHLIST = (char**) malloc ( MXPATH_SNDATA_SIM * sizeof(char*) );
  for(ipath=0; ipath < MXPATH_SNDATA_SIM; ipath++ ) 
    { PATHLIST[ipath] = (char*) malloc ( MXPATHLEN * sizeof(char) ); }
  NPATH = getList_PATH_SNDATA_SIM(PATHLIST);

  // tack on default SIM dir (Sep 2019)
  int IPATH_SIM_DEFAULT = NPATH;
  sprintf(PATHLIST[NPATH], "%s/SIM", SNDATA_ROOT, VERSION); NPATH++ ;

  if ( LDMP ) 
    { printf(" xxx DATADIR = '%s' \n", DATADIR); fflush(stdout); }

  for(ipath = 0 ; ipath < NPATH; ipath++ ) {

    if ( idir < MXPATH_SNDATA_SIM ) {
      sprintf(tmpDir[idir],  "%s/%s" , PATHLIST[ipath],  VERSION );
      sprintf(tmpFile[idir], "%s/%s.LIST", tmpDir[idir], VERSION  );

      if ( LDMP) {
	printf(" xxx check PATHLIST[%d] = '%s' \n", ipath, PATHLIST[ipath] );
	fflush(stdout); 
      }

      // Sep 12 2019: abort if DATADIR corresponds to any SIM path
      if ( strcmp(DATADIR,PATHLIST[ipath])== 0 ) {
	printf("\n PRE-ABORT DUMP: \n");
	printf("\t PRIVATE_DATA_PATH = '%s' \n", DATADIR);

	if ( ipath == IPATH_SIM_DEFAULT ) {
	  sprintf(c1err,"PRIVATE_DATA_PATH cannot be the same as");
	  sprintf(c2err,"$SNDATA_ROOT/SIM");
	}
	else {
	  sprintf(c1err,"PRIVATE_DATA_PATH cannot match any path in");
	  sprintf(c2err,"$SNDATA_ROOT/SIM/%s", PATH_SNDATA_SIM_LIST );
	}
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

    }
    idir++ ;
  }

  for(ipath=0; ipath < MXPATH_SNDATA_SIM; ipath++ ) 
    { free(PATHLIST[ipath]) ; }
  free(PATHLIST);

  // -------
  NDIR_CHECK = idir;

  if ( NDIR_CHECK >= MXPATH_SNDATA_SIM ) {
    sprintf(c1err,"NDIR_CHECK=%d exceeds bound", NDIR_CHECK);
    sprintf(c2err,"of MXPATH_SNDATA_SIM=%d", MXPATH_SNDATA_SIM);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // set DATADIR to directory with LIST file
  NFOUND = 0 ;

  for(idir=0; idir < NDIR_CHECK; idir++ ) {
    if ( (fp = fopen(tmpFile[idir], "rt")) != NULL )  { 
      sprintf(DATADIR, "%s", tmpDir[idir]); 
      idirFOUND[NFOUND] = idir;
      fclose(fp);
      NFOUND++ ;
    }
  } // end of idir


  if ( NFOUND == 0 ) {
    printf("\n\n PRE-ABORT DUMP: \n");

    for(idir=0; idir < NDIR_CHECK; idir++ ) {
      printf("   data not in '%s' \n", tmpDir[idir] );
    }	   
    sprintf(c1err,"Could not find SNANA version '%s'",  VERSION );
    sprintf(c2err,"Check directories above");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( NFOUND > 1 ) {
    printf("\n PRE-ABORT DUMP: \n");
    for(ifound=0; ifound < NFOUND; ifound++ ) {
      idir = idirFOUND[ifound];
      printf("   Found %s \n", tmpDir[idir] );
    }
    sprintf(c1err,"Found %d VERSION='%s'", NFOUND, VERSION);
    sprintf(c2err,"Only one unique version allowed.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  // always use DATADIR to construct name of list-file and readme file.
  sprintf(LISTFILE,   "%s/%s.LIST",   DATADIR, VERSION );
  sprintf(READMEFILE, "%s/%s.README", DATADIR, VERSION );

  return SUCCESS ;

}  // end of function


int getinfo_photometry_version__(char *VERSION, char *DATADIR, 
				 char *LISTFILE, char *READMEFILE) {
  return getInfo_PHOTOMETRY_VERSION(VERSION,DATADIR,LISTFILE,READMEFILE);
}


// ==========================================
FILE *openFile_PATH_SNDATA_SIM(char *mode) {

  // Open file for 
  //  + reading (mode='read')
  //  + append  (mode='append')
  //
  // modeArg[2] -> modeArg[4] (fix Mac issue)
  //

  char fileName[MXPATHLEN], SNDATA_ROOT[MXPATHLEN] ;
  char modeArg[4];
  FILE *fp ;

  //  char fnam[] = "openFile_PATH_SNDATA_SIM" ;

  // ------------- BEGIN --------------

  sprintf(SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") ) ;

  // hard-wire name of file with list of alternate PATH_SNDATA_SIM
  sprintf(fileName, "%s/SIM/%s", SNDATA_ROOT, PATH_SNDATA_SIM_LIST );
  sprintf(modeArg, "%ct", mode[0] );
  fp = fopen(fileName,modeArg);
  //  printf("\n Open %s in %s-mode (%s)\n", fileName, mode, modeArg);
  return(fp) ;

} // end openFile_PATH_SNDATA_SIM


void add_PATH_SNDATA_SIM(char *PATH) {

  // Created Nov 18 2017
  // Add PATH to file in $SNDATA_ROOT/SIM so that analysis codes 
  // know where else to check for simulated data files.  If PATH 
  // name is already written, don't write anything.

  FILE *fp;
  int  EXIST=0;
  int  NPATH_DEJA=2; // accounts for default /lcmerge and /SIM 
  char PATH_ENVreplace[MXPATHLEN];
  char ftmp[MXPATHLEN];
  char fnam[] = "add_PATH_SNDATA_SIM" ;

  // -------------- BEGIN ---------------

  // first read file to see if PATH already exists
  fp = openFile_PATH_SNDATA_SIM("read");
  sprintf(PATH_ENVreplace,"%s", PATH);  
  ENVreplace(PATH_ENVreplace,fnam,1);

  if ( fp ) {
    while( (fscanf(fp, "%s", ftmp)) != EOF) {
      ENVreplace(ftmp,fnam,1);
      NPATH_DEJA++ ;
      if ( strcmp(PATH_ENVreplace,ftmp) == 0  ) { EXIST=1; }
    }
    fclose(fp);
  }

  // if PATH does not already exist, append it to list file
  if ( EXIST == 0 ) {
    if ( NPATH_DEJA >= MXPATH_SNDATA_SIM ) {
      printf("\n PRE-ABORT DUMP: \n");
      printf("   NPATH_DEJA = %d (includes /lcmerge and /SIM) \n", NPATH_DEJA);
      printf("   MXPATH_SNDATA_SIM = %d \n", MXPATH_SNDATA_SIM);

      sprintf(c1err,"Too many paths defined in ");
      sprintf(c2err,"$SNDATA_ROOT/SIM/%s", PATH_SNDATA_SIM_LIST);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }
    fp = openFile_PATH_SNDATA_SIM("append");
    fprintf(fp, "%s\n", PATH);
    fclose(fp);
  }

  return;

} // end add_PATH_SNDATA_SIM


int  getList_PATH_SNDATA_SIM(char **pathList) {

  // Created Nov 2017
  // Return list (pathList) of user-defined PATH_SNDATA_SIM .
  // Function returns NPATH, number on list.
  // This list is used by analysis codes to look for
  // simulated data files so that user doesn't have
  // to worry about writing sim files to alternate
  // directories
  //
  // Mar 30 2019: refactor so that undefined paths are skipped
  //              without aborting.
  //
  int ENVstat, NPATH=0;
  FILE *fp ;
  char path[MXPATHLEN] ;
  char fnam[] = "getList_PATH_SNDATA_SIM" ;

  // ------------ BEGIN ---------------- 

  // open list file in read mode
  fp = openFile_PATH_SNDATA_SIM("read");
  if ( !fp ) { return(0); }   // return if it does not exist.
 
  // scoop up each PATH in list file.
  while( (fscanf(fp, "%s", path)) != EOF) {

    // do NOT abort on invalid ENV
    ENVstat = ENVreplace(path,fnam,0); 
    if ( ENVstat != SUCCESS ) { continue ; }

    if (NPATH<MXPATH_SNDATA_SIM) { sprintf(pathList[NPATH], "%s", path); }
    NPATH++ ;
  }
  fclose(fp);

  // abort if too many paths.
  if ( NPATH >= MXPATH_SNDATA_SIM ) {
    sprintf(c1err,"NPATH=%d exceeds bound of MXPATH_SNDATA_SIM=%d",
	    NPATH, MXPATH_SNDATA_SIM );
    sprintf(c2err,"Check opened file above." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  return(NPATH) ;

} // end getList_PATH_SNDATA_SIM


// =============================================================
void arrayStat(int N, double *array, double *AVG, double *RMS) {

  // For input *array return *AVG and *RMS

  int i;
  double XN, avg, sqsum, rms, tmpdif ;

  // ----------- BEGIN ------------
  *AVG = *RMS = 0.0 ;
  if ( N <= 0 ) { return ; }

  avg  = rms = sqsum = 0.0 ;
  XN   = (double)N ;

  for ( i=0; i < N; i++ ) { avg += *(array+i) ; }
  avg /= XN ; 

  for ( i=0; i < N ; i++ ) {
    tmpdif =  ( *(array+i) - avg ) ;
    sqsum += (tmpdif*tmpdif) ;
  }
  rms = sqrt(sqsum/XN) ;
  
  // load output array.
  *AVG = avg ;
  *RMS = rms ;

} // end of arrayStat


double RMSfromSUMS(int N, double SUM, double SQSUM) {

  // Created Aug 2017
  // Compute RMS from sums

  double RMS = 0.0 ;
  double XN  = (double)N;
  if ( N == 0 ) { return(RMS); }

  double ARG = SQSUM/XN - pow((SUM/XN),2.0) ;
  if ( ARG > 0.0 ) { RMS = sqrt(ARG); }

  return(RMS);

} // end RMSfromSUMS

// =============================================
void remove_quote(char *string) {

  // remove quote(s) from string
  // Input *string is returned without quotes.

  char q[]=  "'" ;
  char qq[] = "\"" ;
  if ( strstr(string,q) == NULL && strstr(string,qq) == NULL ) 
    { return ; }


  int i, i2, lens = strlen(string);
  char *string_orig = (char*) malloc( (lens+1) * sizeof(char) );  
  char ctmp[2] ;
  i2=0;  sprintf(string_orig,"%s", string);  string[0]=0; 
  for(i=0; i<lens; i++ ) {
    sprintf(ctmp, "%c", string_orig[i] ) ;
    if ( strcmp(ctmp,q)  == 0 ) { continue ; }
    if ( strcmp(ctmp,qq) == 0 ) { continue ; }
    sprintf(&string[i2], "%c", string_orig[i] );
    i2++ ;
  }

  free(string_orig);

  return ;

} // end remove_quote

void extractStringOpt(char *string, char *stringOpt) {

  // Created Aug 23 2016                          
  //                                                         
  // For input *string = 'blabla(stringOpt)moreBlaBla' 
  // return 
  //   *string    = 'blablamoreBlaBla' 
  //   *stringOpt = 'stringOpt'
  //  
  // i.e., remove () from *string and return contents of () 
  // in *stringOpt.  If there are no (), then *string is
  // returned unchanged and *stringOpt="".   

  int  i,lens = strlen(string);
  int  L=0, R=0 ;                // Left & Right parentheses logicals 
  char *stringLocal, ctmp[2] ;
  //  char fnam[] = "extractStringOpt";

  // --------------- BEGIN --------------    

  stringOpt[0]=0;
  if ( strstr(string,"(") == NULL ) { return ; }

  stringLocal = (char*) malloc ( (lens+2)*sizeof(char) );
  sprintf(stringLocal, "%s", string);  string[0]=0; ;

  for(i=0; i < lens; i++ ) {
    sprintf(ctmp, "%c", stringLocal[i] ) ;
    if ( *ctmp == '(' ) { L=1; R=0; continue ; }
    if ( *ctmp == ')' ) { R=1; L=0; continue ; }

    if ( L==1 && R==0 )
      { strcat(stringOpt,ctmp); }
    else
      { strcat(string,ctmp); }

  }

  free(stringLocal);
  return ;

} // end extractStringOpt


// ===============================================
void  remove_string_termination(char *STRING, int LENTOT) {

  // Created Jan 2014 by R.Kessler
  // remove string termination and replace with padding
  // up to length LENTOT. Note that inputs *STRING is modified.
  int ic, LEN;
  // -------------- BEGIN -------------
  LEN = strlen(STRING);
  for(ic=LEN; ic < LENTOT; ic++ )  { STRING[ic] = ' ' ;  }
} 

// ====================================
void trim_blank_spaces(char *string) {

  // April 2013
  // return string without blank spaces.
  // Assume that first blank space after char is end of string
  // Examples:
  //   "BLA1   BLA2   " -> "BLA1".
  //   "   BLA1   BLA2" -> "BLA1".
  // Note that input string is overwritten !
  // 
  // strlen() is NOT used in case there is no '\0' termination.
  //
  // Jan 10 2017: check for termination char '\0' so that it now
  //              works for properly terminated strings with no
  //              extra blank spaces.
  //
  // Mar 13 2019:
  //  Check \r and \n for <CR> or line-feed. 
  //  

  int MXchar, i, FOUNDCHAR, ISCHAR, ISBLANK, ISTERM ;
  char *tmpString, c1[2] ;
  //  char fnam[] = "trim_blank_spaces" ;

  // -------------- BEGIN ---------

  MXchar  = 1000 ;

  // alloate temporary/huge array since we don't know the length
  // and cannot use strlen.
  tmpString = (char*) malloc( sizeof(char) * MXchar ) ;
  tmpString[0]=0;

  if ( strlen(string) == 0 ) { return; } // Aug 2014

  //  printf(" xxx trim string='%s' \n",string); fflush(stdout);

  // transfer non-null chars in string to tmpString
  FOUNDCHAR = 0 ;

  for ( i=0; i < MXchar-1; i++ ) {
    sprintf(c1, "%c", string[i] );
    ISBLANK = ( string[i] == ' '  )  ;
    ISCHAR  = ( ISBLANK == 0  );
    ISTERM  = ( string[i] == '\0' || string[i] == '\n' || string[i]=='\r' ) ;

    /*
    printf(" xxx ----------------- \n");
    printf(" xxx %s: '%s'(%d): ISCHAR=%d  ISBLANK=%d FOUNDCH=%d  ISTERM=%d \n",
	   fnam, string, i, ISCHAR, ISBLANK, FOUNDCHAR, ISTERM);
    printf(" xxx %s: tmpString = '%s'  will; add c1='%s' \n",  
	   fnam, tmpString, c1);
    */

    if ( ISCHAR ) { FOUNDCHAR++ ; } // once set, stays set

    if ( ISBLANK && FOUNDCHAR ) { goto DONE ; }
    if ( ISTERM               ) { goto DONE ; }

    if ( ISCHAR ) 
      {  sprintf(tmpString,"%s%s", tmpString, c1); }
  }


  // overwrite input arg.
 DONE:

  sprintf(string, "%s", tmpString);
  free(tmpString);
  
} // end of trim_blank_spaces

void splitString(char *string, char *sep, int MXsplit,
		 int *Nsplit, char **ptrSplit) {

  // Created July 2016
  //
  // Inputs:
  //    *string  : string to split (preserved)
  //    *sep     : separator, e.g., ',' or ' ' or '+'
  //    MXsplit  : abort if Nsplit >= MXsplit
  //
  // Output :
  //  *Nsplit     : number of split elements return in ptrSplit
  //  **ptrSplit  : array of pointers to split elements
  // ---------------

  int LEN, N;
  char *localString, *ptrtok ;
  char fnam[] = "splitString" ;

  // ------------ BEGIN ---------------
  LEN         = strlen(string);
  localString = (char*) malloc( (LEN+10) * sizeof(char) );
  sprintf(localString, "%s", string);
  ptrtok      = strtok(localString,sep) ; // split string
  N=0;

  while ( ptrtok != NULL  ) {
    if ( N < MXsplit ) 
      { sscanf(ptrtok,"%s", ptrSplit[N] );  }
    ptrtok = strtok(NULL, sep);
    N++ ;
  }

  if ( N > MXsplit ) {
    printf("\n PRE-ABORT DUMP: \n");
    printf("  string to split: '%s' \n", string);
    printf("  split separator: '%s' \n", sep);
    sprintf(c1err,"Nsplit = %d ", N );
    sprintf(c2err,"Exceeds bound MXsplit=%d", MXsplit);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  *Nsplit = N ; // load output arg
  free(localString);

  return ;

} // end splitString

void splitString2(char *string, char *sep, int MXsplit,
		  int *Nsplit, char **ptrSplit) {

  // ----
  // Use strtok_r ... much faster than strtok, 
  // but input string is destroyed.
  // ----
  //
  // Created July 2016                                
  //                                                                  
  // Inputs:                                                           
  //    *string  : string to split (destroyed)
  //    *sep     : separator, e.g., ',' or ' ' or '+'            
  //    MXsplit  : abort if Nsplit >= MXsplit                    
  //                                   
  // Output : 
  //  *Nsplit     : number of split elements return in ptrSplit 
  //  **ptrSplit  : array of pointers to split elements    
  //
  // Dec 27 2017: avoid <CR> in case of fgets scooping up extra
  //              blank spaces.
  // ---------------                                             

  int   N;
  char *localString, *token ;
  //  char fnam[] = "splitString2" ;

  // ------------ BEGIN ---------------

  localString = string ;
  N=0;

  while (( token = strtok_r(localString, sep, &localString ))) {

    /*
    printf(" xxx %s: token = '%s' L=%d \n", 
    fnam, token, strlen(token));  */

    if ( token[0] != '\0'  && token[0] != '\n' ) {
      if ( N < MXsplit ) { sprintf(ptrSplit[N],"%s", token ); }
      N++ ;
    }
  }
  *Nsplit = N ; // load output arg       

  return ;

}  // end of splitString2

void split2floats(char *string, char *sep, float *fval) {

  // Created Jun 26 2019
  // for *string = 'xxx[sep]yyy,
  // returns fval[0]=xxx and fval[1]=yyy.
  // Example:
  //   Input string   = 1.3,4.6
  //   Output fval[0] = 1.3
  //   Output fval[1] = 4.6
  //
  // Example:
  //   Input string   =  1.3
  //   Output fval[0] =  1.3
  //   Output fval[1] =  1.3
  //
  int Nsplit ;
  char cnum[2][40], *cptr[2];  cptr[0]=cnum[0]; cptr[1]=cnum[1];
  char fnam[] = "split2floats" ;
  // ---------------- BEGIN --------------------

  fval[0] = fval[1] = -9.0 ;

  if ( strstr(string,sep) == NULL ) 
    { sscanf(string, "%f", &fval[0] ); fval[1]=fval[0];  return ;  }
  

  // split the string by the sep input
  splitString(string, sep, 2, &Nsplit, cptr);
  if ( Nsplit != 2 ) {
    sprintf(c1err,"Invalid Nsplit=%d (expected 2)", Nsplit);
    sprintf(c2err,"Input string='%s'  sep='%s' ", string, sep);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  sscanf(cnum[0], "%f", &fval[0] );
  sscanf(cnum[1], "%f", &fval[1] );

  return;
} // end split2floats


// ********************************************************
void read_GRIDMAP(FILE *fp, char *KEY_ROW, char *KEY_STOP, 
		  int IDMAP, int NDIM, int NFUN, int OPT_EXTRAP, int MXROW,
                  char *callFun, GRIDMAP *GRIDMAP_LOAD ) {

  // Mar 2019
  // Utility to read mutil-D map from file and call
  // init_inter_GRIDMAP to create & store GRIDMAP_LOAD.
  // Beware that uniform map-bins are strictly enforce;
  // ABORTs on non-uniform bin.
  //
  // Inputs:
  //   *fp        : already-opened file to read
  //  KEY_ROW     : NVAR columns follows this row-key
  //  KEY_STOP    : stop reading when this key is reached;
  //              " default is to stop reading on blank line.
  //  IDMAP       : integer ID of GRIDMAP_LOAD
  //  NDIM        : number of dimensions of map
  //  NFUN        : number of functions of map
  //  OPT_EXTRAP  : flag for extrapolation outside map range
  //  MXROW       : abort if NROW > MXROW 
  //  callFun     : name of calling function (for error messages)
  //
  // Output:
  //    GRIDMAP_LOAD
  //
  //
  // Apr 12 2019: abort if 10 or more rows read without valid key

  int   READ_NEXTLINE = 1 ;
  int   NROW_READ     = 0 ;
  int   NVARTOT = NDIM + NFUN;
  char  *VARLIST = GRIDMAP_LOAD->VARLIST ; // for comment only
  double DUMVAL = -999.0 ;

  int   MEMD   = sizeof(double);
  int   MEMVAR = NVARTOT  * sizeof(double*);
  int   MEMROW = MXROW    * MEMD ;

  int  NBADBIN = 0 ;
  int  NLINE   = 0 ;
  double **TMPMAP2D ;  // [0:NVARTOT-1][MXROW-1]
  double *TMPVAL, *TMPVAL_LAST, *DIFVAL_LAST, DDIF, DIF;

  int   ivar, NWD, ISKEY_ROW, EXTRA_WORD_OK ;
  int   LDIF1, LDIF2, ivar2, NROW_SKIP=0 ;
  char  LINE[200], word[40], MAPNAME[100] ;
  char fnam[] = "read_GRIDMAP" ;
 
  // ----------- BEGIN -------------

  // create generic MAPNAME using row key and IDMAP
  //  sprintf(MAPNAME,"%s%3.3d", KEY_ROW, IDMAP );
  sprintf(MAPNAME,"%s", KEY_ROW );

  // allocate arrays to monitor uniform binning.
  TMPVAL      = (double*) malloc(NVARTOT * MEMD );
  TMPVAL_LAST = (double*) malloc(NVARTOT * MEMD );
  DIFVAL_LAST = (double*) malloc(NVARTOT * MEMD );
  for(ivar=0; ivar<NVARTOT; ivar++) {
    TMPVAL[ivar] = DUMVAL;
    TMPVAL_LAST[ivar] = DUMVAL;
    DIFVAL_LAST[ivar] = DUMVAL;
  }

  // alloate temp 2D array to read map
  TMPMAP2D = (double**) malloc(MEMVAR);
  for(ivar=0; ivar<NVARTOT; ivar++) {TMPMAP2D[ivar]=(double*)malloc(MEMROW);} 


  while ( READ_NEXTLINE ) {
    LINE[0] = 0 ;
    fgets(LINE,200,fp);  NLINE++ ;
    NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);

    // abort if we read too many lines without finding any valid row keys
    if ( NLINE > 20 && NROW_READ==0 ) {
      sprintf(c1err,"Found no '%s' keys after reading %d lines.",
	      KEY_ROW, NLINE);
      sprintf(c2err,"NDIM=%d, NFUN=%d, callFun=%s", 
	      NDIM, NFUN, callFun );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    // Skip blank line.
    // However, stop reading only after reading at least one valid row;
    // this allows blank line(s) between VARNAMES and first map row.
    if ( NROW_READ > 0  && NWD == 0 )  { READ_NEXTLINE=0; }
    if ( NWD == 0 ) { continue ; }

    get_PARSE_WORD(0,0,word);  

    ISKEY_ROW = 0 ;
    if ( strcmp(word,KEY_ROW) ==0 ) { ISKEY_ROW = 1; }
    if ( strcmp(word,KEY_STOP)==0 ) { READ_NEXTLINE=0; continue; }

    if ( ISKEY_ROW ) {
      
      NROW_SKIP = 0 ;
      // allow comment string on same line as grid data
      EXTRA_WORD_OK = 1 ;
      if ( NWD-1 > NVARTOT ) {
	get_PARSE_WORD(0,NVARTOT+1,word);
	EXTRA_WORD_OK = ( word[0] == '#' ) ;
      }
      //  printf(" xxx extra word = '%s'  OK=%d \n",word, EXTRA_WORD_OK);

      if ( (NWD-1 < NVARTOT) || (!EXTRA_WORD_OK) ) {
	sprintf(c1err,"Expected NVARTOT=%d words after '%s' key,",
		NVARTOT, KEY_ROW);
	sprintf(c2err,"but found %d. (NDIM=%d, NFUN=%d, callFun=%s)", 
		NWD-1, NDIM, NFUN, callFun );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }

      for(ivar=0; ivar < NVARTOT; ivar++ ) {
	get_PARSE_WORD(0,1+ivar,word);
	sscanf ( word, "%le", &TMPVAL[ivar] );
	TMPMAP2D[ivar][NROW_READ] = TMPVAL[ivar];

	// check for uniform binning
	DIF = TMPVAL[ivar] - TMPVAL_LAST[ivar];
	if ( DIF > 0.0  && ivar < NDIM && TMPVAL_LAST[ivar]!=DUMVAL ) { 
	  DDIF  = DIF - DIFVAL_LAST[ivar] ;
	  LDIF1 = ( fabs(DDIF/DIF) > .001 ) ; 
	  LDIF2 = ( DIFVAL_LAST[ivar] > 0.0 ) ;
	  if ( LDIF1 && LDIF2 ) {
	    NBADBIN++ ;
	    printf(" ERROR: non-uniform bin at '%s'=", VARLIST );
	    for(ivar2=0; ivar2 < NVARTOT; ivar2++ ) 
	      { printf("%.3f ", TMPVAL[ivar2] ); }
	    printf(" (row %d)\n", NROW_READ );	    fflush(stdout);
	  }
	  DIFVAL_LAST[ivar] = DIF; 
	} // end DIF>0
	// end uniform bin check

	TMPVAL_LAST[ivar] = TMPVAL[ivar];
      } // end ivar loop

      NROW_READ++ ;
      if ( NROW_READ >= MXROW ) {
	sprintf(c1err,"NROW_READ=%d exceeds MXROW=%d", NROW_READ, MXROW);
	sprintf(c2err,"NDIM=%d  NFUN=%d  callFun=%s", 
		NDIM, NFUN, callFun );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }

    } // end KEY_ROW
    else {
      // 4.2019: abort if too many rows have invalid key
      NROW_SKIP++ ;
      if ( NROW_SKIP >= 10 ) { 
	printf("\n PRE-ABORT DUMP: \n");
	printf("   Last line read: %s\n", LINE);
	sprintf(c1err,"Read %d rows without valid row-key, "
		"stop-key, or blank line.", NROW_SKIP );
	sprintf(c2err,"KEY_ROW='%s'  KEY_STOP='%s'  callFun=%s", 
		KEY_ROW, KEY_STOP, callFun );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }
    }

  } // end while


  // -------------------------------------------------
  // ABORT on non-uniform bins
  if ( NBADBIN > 0 ) {
    sprintf(c1err,"%d non-uniform bin errors detected", NBADBIN);
    sprintf(c2err,"Check %s map. ", KEY_ROW );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  // ----------------
  printf("    Load GRIDMAP-%3.3d '%s(%s)'  NROW=%d \n",
	 IDMAP, MAPNAME, VARLIST, NROW_READ); fflush(stdout);

  init_interp_GRIDMAP(IDMAP, MAPNAME, NROW_READ, NDIM, NFUN, OPT_EXTRAP,
		      TMPMAP2D, 
		      &TMPMAP2D[NDIM],
		      GRIDMAP_LOAD  );       // <== returned

  // free temp memory
  for(ivar=0; ivar < NVARTOT; ivar++ )  { free(TMPMAP2D[ivar]); }
  free(TMPMAP2D);
  free(TMPVAL); free(TMPVAL_LAST); free(DIFVAL_LAST);
  return ;

} // end read_GRIDMAP

// ==============================================================
void init_interp_GRIDMAP(int ID, char *MAPNAME, int MAPSIZE, 
			 int NDIM, int NFUN, int OPT_EXTRAP,
			 double **GRIDMAP_INPUT, double **GRIDFUN_INPUT,
			 GRIDMAP *gridmap ) {

  // Created July 2011 by R.Kessler
  // Return struct *gridmap to assist in mult-dimensional interp .
  // This struct contains all the binning info for each dimension.
  //
  // Arguments
  // (I) ID        reference id
  // (I) MAPNAME   human-readable name for error message
  // (I) MAPSIZE   total number of bins in gridmap
  // (I) NDIM      number of dimensions = number of variables
  // (I) NFUN      Number of functions on same GRID
  // (I) **GRIDMAP_INPUT[idim][i=0 to MAPSIZE-1] 
  // (I) **GRIDFUN_INPUT[ifun][i=0 to MAPSIZE-1] = function values
  // (O) *gridmap  structure to return 
  //         (to be passed later to interp_GRIDMAP)
  //
  // Jun 15 2016: 
  //  + add char MAPNAME argument for error messages.
  //  + print better message for non-uniform binning.
  //  + note that binning check only works for 1D map ... fix later.
  // 
  // Feb 12 2018: 
  //   + malloc gridmap here, instead of externally
  //   + refactor so that all local indices are 0 to N-1 (not 1-N)
  //
  // Mar 13 2018:
  //   + fix bug malloc-ing FUNVAL : I8p -> I8p * NFUN

  int idim, ifun, i, NBIN, igrid_tmp, igrid_1d[100] ;
  double VAL, VALMIN, VALMAX, VALBIN, LASTVAL, RANGE, DIF ;
  double FUNVAL, RANGE_CHECK, RATIO;
  char fnam[] = "init_interp_GRIDMAP" ;

  int I4  = sizeof(int);
  int I8  = sizeof(double);
  int I8p = sizeof(double*);
  // --------- BEGIN ------------
  
  gridmap->NBIN      = (int     *)malloc(I4*NDIM+I4);
  gridmap->VALMIN    = (double  *)malloc(I8*NDIM+I8);
  gridmap->VALMAX    = (double  *)malloc(I8*NDIM+I8);
  gridmap->VALBIN    = (double  *)malloc(I8*NDIM+I8);
  gridmap->RANGE     = (double  *)malloc(I8*NDIM+I8);
  gridmap->FUNVAL    = (double **)malloc(I8p*NFUN);
  gridmap->FUNMIN    = (double  *)malloc(I8*NFUN);
  gridmap->FUNMAX    = (double  *)malloc(I8*NFUN);
  gridmap->INVMAP    = (int     *)malloc(I4*MAPSIZE+I4);
  for(ifun=0; ifun < NFUN; ifun++ ) 
    {  gridmap->FUNVAL[ifun] = (double *)malloc(I8*MAPSIZE);  }


  VALBIN = 0.0 ;
  for ( idim=0; idim < NDIM; idim++ ) {

    NBIN    =  0 ;
    VALMIN  = +1.0E12 ;
    VALMAX  = -1.0E12 ;
    LASTVAL = -999. ;

    for ( i=0; i < MAPSIZE; i++ ) {
      VAL = GRIDMAP_INPUT[idim][i] ;

      if ( VAL > VALMAX  ) { VALMAX = VAL ; }
      if ( VAL < VALMIN  ) { VALMIN = VAL ; }
      if ( VAL > LASTVAL && LASTVAL != -999. ) 
	{ VALBIN = VAL - LASTVAL ; }

      LASTVAL = VAL ;
    } 

    RANGE = VALMAX - VALMIN ;
    NBIN = (int)( (RANGE+0.001*VALBIN) / VALBIN ) + 1;

    // load output struct
    gridmap->ID           = ID ;
    gridmap->NDIM         = NDIM ;
    gridmap->NBIN[idim]   = NBIN ;
    gridmap->VALMIN[idim] = VALMIN ;
    gridmap->VALMAX[idim] = VALMAX ;
    gridmap->VALBIN[idim] = VALBIN ;
    gridmap->RANGE[idim]  = RANGE ;
    gridmap->NFUN         = NFUN ;
    gridmap->NROW         = MAPSIZE ;
    gridmap->OPT_EXTRAP   = OPT_EXTRAP ;

    // make sure that VALBIN x integer = RANGE
    RANGE_CHECK = (double)(NBIN-1) * VALBIN;
    RATIO       = RANGE_CHECK/RANGE ;

    if ( fabs(RATIO-1.0) > 1.0E-4 ) {
      printf("\n PRE-ABORT DUMP:\n");
      printf("\t VALMAX - VALMIN  = %le (%le to %le)\n", 
	     RANGE, VALMIN, VALMAX );
      printf("\t (NBIN-1)*BINSIZE = %le (%d x %le) \n",
	     RANGE_CHECK, NBIN-1, VALBIN);
      printf("\t Ratio-1 = %le \n", RATIO-1. );

      sprintf(c1err,"Non-uniform binning for idim=%d", idim);
      sprintf(c2err,"Check map = '%s' ", MAPNAME );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

  } // idim

  // now store function value at each node along with
  // mapping between mult-D indices and 1D index
  // Load gridmap->FUNVAL  and gridmap->IGRIDMAP

  init_1DINDEX(ID, NDIM, &gridmap->NBIN[0] ) ;

  for ( ifun=0; ifun < NFUN; ifun++ ) { 
    gridmap->FUNMIN[ifun] = +999999.0 ;
    gridmap->FUNMAX[ifun] = -999999.0 ;
  }


  for ( i=0; i < MAPSIZE; i++ )  {  

    for ( ifun=0; ifun < NFUN; ifun++ ) { 
      FUNVAL = GRIDFUN_INPUT[ifun][i] ;
      gridmap->FUNVAL[ifun][i] = FUNVAL ;
      if(FUNVAL < gridmap->FUNMIN[ifun]) { gridmap->FUNMIN[ifun] = FUNVAL; }
      if(FUNVAL > gridmap->FUNMAX[ifun]) { gridmap->FUNMAX[ifun] = FUNVAL; }
    }

      for ( idim=0; idim < NDIM; idim++ ) {      
	VAL    = GRIDMAP_INPUT[idim][i] ;
	VALMIN = gridmap->VALMIN[idim] ;
	VALBIN = gridmap->VALBIN[idim] ;
	DIF    = VAL - VALMIN ;
	if ( VALBIN < 1.0E-9 ) 
	  { igrid_1d[idim] = 0 ; }
	else
	  { igrid_1d[idim] = (int)((DIF+1.0E-9)/VALBIN); }  //  + 1 ; }

      }

      igrid_tmp = get_1DINDEX(ID, NDIM, &igrid_1d[0] ) ; 

      if ( igrid_tmp < 0 || igrid_tmp >= MAPSIZE ) {

	printf("\n PRE-ABORT DUMP for MAPNAME=%s: \n", MAPNAME );
	for ( idim=0; idim < NDIM; idim++ ) {  
	  VAL    = GRIDMAP_INPUT[idim][i] ;
	  VALMIN = gridmap->VALMIN[idim] ;
	  VALBIN = gridmap->VALBIN[idim] ;
	  printf("   idim=%d : VAL=%10.3f  BIN=%10.3f  MIN=%10.3f  "
		 "igrid_1d=%d \n",
		 idim, VAL, VALBIN, VALMIN, igrid_1d[idim] );
	}
	printf("\t Probably have non-uniform binning.\n");
	       
	sprintf(c1err,"Invalid igrid_tmp=%d (ID=%d, MAPSIZE=%d)", 
		igrid_tmp, ID, MAPSIZE );
	sprintf(c2err,"original NDIM=%d  igrid=%d ", NDIM, i ) ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }

      gridmap->INVMAP[igrid_tmp] = i ;

  } // end loop over MAPSIZE
  

} // end of init_interp_GRIDMAP


int interp_GRIDMAP(GRIDMAP *gridmap, double *data, double *interpFun ) {

  // Created Jul 3, 2011
  // Do multi-dimensional interpolation.
  // (I) *gridmap is returned from init_GRIDMAP
  // (I) *data is the multi-dimensional data point to interpolate
  // (O) *interpFun is the interpolated function value, or array
  //                of function values for multiple functions
  //
  // Function returns  0  if *data is withing the grid;
  // Function returns -1  if *data is outside the grid.
  //
  // Note that init_interp_GRIDMAP must be called first
  // to initialize the gridmap structure.
  //
  // Jul 25, 2011; remove   WGT_SUM <= 0 test since it can be
  //               zero if the function is zero nearby.
  //
  // Aug 28, 2011: 
  //     fix bug for when TMPVAL is within "EPSILON" of TMPMAX
  //
  // Mar 14 2019: 
  //  + check OPT_EXTRAP option
  //  + return SUCCESS or ERROR instead of hard-coded values.

  int 
    ivar, ifun, NFUN, NVAR, ID, igrid, MSK, NBIN
    ,NCORNERS, icorner, igrid_tmp, igrid_1D, g
    ,igrid_cell[100], igrid_var[100], IGRID_VAR[100]
    ;

  double  
    WGT_SUM[100], CORNER_WGTSUM, CORNER_WGT
    ,TMPVAL, TMPDIF, TMPMIN, TMPMAX, TMPBIN, TMPRANGE, xgrid, XNBIN
    ,GRIDFRAC[100], FUNVAL[100]
    ;

  double EPSILON = 1.0E-8 ;

  int  LDMP=0;
  char fnam[] = "interp_GRIDMAP" ;

  // ---------- BEGIN ------------

  ID   = gridmap->ID ;
  NVAR = gridmap->NDIM ;
  NFUN = gridmap->NFUN ;

  for  ( ifun=0; ifun < NFUN; ifun++ )   {  
    *(interpFun + ifun) = 0.0 ; 
    WGT_SUM[ifun] = 0.0 ;
  }
  CORNER_WGTSUM = 0.0 ;

  if ( LDMP ) 
    { printf(" xxxxx ------------- START DUMP ----------------- \n"); }

  // get central index and grid-frac in each dimension
  for ( ivar=0; ivar < NVAR; ivar++ ) {
    TMPVAL   = data[ivar] ;
    TMPMIN   = gridmap->VALMIN[ivar] ;
    TMPMAX   = gridmap->VALMAX[ivar] ;
    TMPBIN   = gridmap->VALBIN[ivar] ;
    TMPRANGE = TMPMAX - TMPMIN ;

    // check extrap option
    if ( gridmap->OPT_EXTRAP ) {
      if ( TMPVAL < TMPMIN ) { TMPVAL = TMPMIN + (TMPRANGE*1.0E-12); }
      if ( TMPVAL > TMPMAX ) { TMPVAL = TMPMAX - (TMPRANGE*1.0E-12); }	
    }

    /*
    if ( TMPVAL < TMPMIN  || TMPVAL > TMPMAX ) {
      printf(" %s ERROR: TMPVAL=%le not between %le and %le \n",
	     fnam, TMPVAL, TMPMIN, TMPMAX);
      fflush(stdout);
      return(ERROR);
    } 
    */

    if ( TMPVAL < TMPMIN ) { return(ERROR) ; }
    if ( TMPVAL > TMPMAX ) { return(ERROR) ; }


    TMPDIF  = TMPVAL - TMPMIN ;
    if ( TMPBIN == 0.0 )
      { XNBIN = 0.0 ; igrid = 0; }
    else if ( (TMPMAX - TMPVAL)/TMPRANGE < EPSILON  )  { 
      XNBIN = (TMPDIF - TMPRANGE*EPSILON)/TMPBIN ;
      igrid = (int)(XNBIN) ; 
    }
    else {
      XNBIN = (TMPDIF + TMPRANGE*EPSILON ) / TMPBIN ;
      igrid = (int)(XNBIN); //  + 1; 
    }

    xgrid   = (double)igrid ;

    // store relative cell location: 0-1
    if ( TMPBIN > 0.0 ) 
      {  GRIDFRAC[ivar]  = TMPDIF/TMPBIN - xgrid ; }
    else
      {  GRIDFRAC[ivar]  = 1.0 ; }

    IGRID_VAR[ivar] = igrid ; // store central bin  for each var

    if ( LDMP ) {
      printf(" xxxx VAL=%f  BIN=%f  XNBIN=%f  igrid=%2d GRIDFRAC=%le \n",
	     TMPVAL, TMPBIN, XNBIN, igrid, GRIDFRAC[ivar] );
      fflush(stdout);
    }

  } // ivar


  // determine the grid points at the corners of the
  // NVAR-dimentional cell containing *galpar.
  // Then take weighted average of WGTMAP at each corner.

  double XN ;
  XN = (double)NVAR ;

  NCORNERS = (int)pow(2.0,XN);
  for ( icorner=0; icorner < NCORNERS; icorner++ ) {
   
    if ( LDMP ) {
      printf(" xxx --------- Next icorner = %d / %d ------------ \n", 
	     icorner, NCORNERS ); 
      fflush(stdout);
    }

    CORNER_WGT = 1.0 ;
    for ( ivar=0; ivar < NVAR; ivar++ ) {
      //      MSK = 1 << (ivar-1);
      MSK = 1 << (ivar);
      igrid_cell[ivar] = (icorner & MSK)/MSK ; // 0 or 1 only
      igrid_var[ivar]  = IGRID_VAR[ivar] + igrid_cell[ivar];

      if ( LDMP  ) {
	printf("\t xxxxxx ivar=%d : cell=%d  igrid_var=%d \n",
	       ivar, igrid_cell[ivar], igrid_var[ivar] );
	fflush(stdout);
      }

      // make sure that igrid_var is valid
      NBIN = gridmap->NBIN[ivar] ;   g = igrid_var[ivar];
      if ( g < 0 || g >= NBIN  ) {
	//	TMPVAL  = *(data+ivar-1) ;
	TMPVAL  = data[ivar] ;
	TMPMIN  = gridmap->VALMIN[ivar] ;
	TMPMAX  = gridmap->VALMAX[ivar] ;
	sprintf(c1err,
		"Invalid igrid_var[ivar=%d]=%d  (NBIN=%d  cell=%d icorner=%d)",
		ivar, g, NBIN, igrid_cell[ivar], icorner ) ;
	sprintf(c2err, "VAL=%f  VALMIN/MAX = %f / %f", 
		TMPVAL, TMPMIN, TMPMAX);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      
      if ( igrid_cell[ivar] ) 
	{ CORNER_WGT *= GRIDFRAC[ivar] ; }
      else
	{ CORNER_WGT *= (1.0 - GRIDFRAC[ivar]) ; }
    }

    CORNER_WGTSUM += CORNER_WGT ;

    // translate 1d indices for each variable into absolute lookup index
    igrid_tmp    = get_1DINDEX( ID, NVAR, &igrid_var[0]);
    igrid_1D     = gridmap->INVMAP[igrid_tmp] ;        

    for  ( ifun=0; ifun < NFUN; ifun++ )   {  
      FUNVAL[ifun]       = gridmap->FUNVAL[ifun][igrid_1D];
      WGT_SUM[ifun]     += (CORNER_WGT * FUNVAL[ifun]) ;
    }
    
    if ( LDMP ) {
      printf(" xxx CORNER_[WGT,SUM](%d)=%6.4f,%6.4f  FUNVAL=%f  WGT_SUM=%f\n",
	     icorner, CORNER_WGT, CORNER_WGTSUM, FUNVAL[0], WGT_SUM[0] );
      fflush(stdout);
    }
    
  } // corner

  if ( CORNER_WGTSUM <= 0.0 ) {
    sprintf(c1err,"Could not compute CORNER_WGT for gridmap ID=%d", 
	    gridmap->ID );
    sprintf(c2err,"%s", "data = ") ;
    for ( ivar=0; ivar < NVAR; ivar++ ) 
      { sprintf(c2err,"%s %f", c2err, *(data+ivar) ) ; }
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  for  ( ifun=0; ifun < NFUN; ifun++ )   {  
    *(interpFun+ifun) = WGT_SUM[ifun]   / CORNER_WGTSUM ; 
  }


  if ( LDMP ) {
    printf("  xxxx data=%6.2f , %6.2f  interpFun=%f (WGT_SUM=%f,%f)\n",
	   data[0], data[1], interpFun[0], WGT_SUM[0], CORNER_WGTSUM );
    printf("  xxxx DUMP DONE. \n");
    fflush(stdout);
    
  }

  return(SUCCESS) ;

} // end of interp_GRIDMAP


// ================================================
int  get_1DINDEX(int ID, int NDIM, int *indx ) {

  // Created April 2011 (initial use for hostlib weight-map)
  //
  // Return 1d index for NDIM-dimensional grid.
  // *indx is an array of indices for each dimension.
  // Each *indx value 0 to N-1 
  //
  // Note: must call init_1DINDEX(ID ...) once per ID
  // before calling this function.
  // 
  // If the number of grid-points in each dimension is
  // N1, N2 ... N_NDIM, then the returned index is an
  // integer from 1 to N1*N2* ... N_NDIM.
  // For example, for a 3 dimensional array [4][5][4],
  // the returned index is from 1 to 4*5*4 = 80.
  //
  // Feb 25, 2013: ABORT if *indx exceeds NPT
  // Feb 12, 2018: indx is 0 to N-1 (no longer 1-N)

  char fnam[] = "get_1DINDEX" ;
  int i, offset, INDEX_1D, index_1d, NPT ;

  //------------ BEGIN -------------

  //  printf(" xxxx %s called with ID = %d \n", fnam, ID) ;

  // xxx mark delete  if ( OFFSET_1DINDEX[ID][0] == 0 ) {
  if ( NPT_PERDIM_1DINDEX[ID][0] == 0 ) {
    sprintf(c1err,"ID=%d  is not defined.", ID );
    sprintf(c2err,"%s", "Must first call init_1DINDEX()");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  INDEX_1D = 0;

  for ( i=0; i < NDIM; i++ ) {   
    offset   = OFFSET_1DINDEX[ID][i];
    index_1d =  indx[i] ;       // index in this dimension
    INDEX_1D += (index_1d ) * offset ; // global 1D index

    /*
    printf(" xxx %s: i=%d index_1d=%3d  INDEX_1D=%6d\n",
	   fnam, i, index_1d, INDEX_1D); fflush(stdout);
    */

    // make sure that index does not exceed NPT
    NPT =    NPT_PERDIM_1DINDEX[ID][i] ;
    if ( index_1d >= NPT ) {
      sprintf(c1err,"index_1d=%d exceeds NPT=%d (ID=%d)", 
	      index_1d, NPT, ID );
      sprintf(c2err,"for idim = %d of %d", i, NDIM) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

  }
  return INDEX_1D ;

} // end of get_1DINDEX


void init_1DINDEX(int ID, int NDIM, int *NPT_PERDIM ) {

  // Apr 2011
  // init offsets needed to quickly compute 1d index
  // for multi-dimensional array or grid.
  //
  // ID        = reference ID for this mapping
  // NDIM      = number of dimensions
  // *NPT_PERDIM  = max number of elements in each dimension 
  //
  //
  // Feb 25 2013: store NPT_PERDIM
  // Feb 12 2018: refactor with indices starting at zero
  //

  int LDMP = 0 ;
  int i, NPT, NPT_LAST, OFFSET_LAST, OFFSET ;
  char fnam[] = "init_1DINDEX" ;

  // --------- BEGIN ----------

  if ( ID < 1 || ID >= MXMAP_1DINDEX ) {
    sprintf(c1err,"Invalid ID=%d", ID );
    sprintf(c2err,"Valid ID range is %d to %d", 1, MXMAP_1DINDEX-1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( NDIM < 0 || NDIM >= MXDIM_1DINDEX ) {
    sprintf(c1err,"Invalid NDIM=%d", NDIM );
    sprintf(c2err,"Valid NDIM range is %d to %d", 1, MXDIM_1DINDEX-1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  for ( i=0; i < NDIM; i++ ) {

    NPT_PERDIM_1DINDEX[ID][i]  = NPT_PERDIM[i]; 
    NPT                        = NPT_PERDIM[i]; 
    OFFSET_1DINDEX[ID][i] = OFFSET  = 0 ;

    if ( i > 0 ) {
      NPT_LAST    = NPT_PERDIM[i - 1];  
      OFFSET_LAST = OFFSET_1DINDEX[ID][i-1] ; 
      OFFSET      = OFFSET_LAST * NPT_LAST ;
      OFFSET_1DINDEX[ID][i] = OFFSET;
    }
    else {
      OFFSET_1DINDEX[ID][i] = OFFSET  = 1 ;
    }
      
    
    if ( LDMP ) {
      printf(" xxxx OFFSET_1DINDEX[ID=%d][ivar=%2d] = %7d   "
	     " NPT_PERDIM=%d\n",
	     ID, i, OFFSET, NPT_PERDIM[i] );
    }

  } // end NDIM


} // end of init_1DINDEX

void clear_1DINDEX(int ID) {
  //  printf("  Clear 1DINDEX for ID=%d \n", ID  );
  OFFSET_1DINDEX[ID][0] = 0;
  NPT_PERDIM_1DINDEX[ID][0] = 0 ;
}


// mangled functions for fortran
void clear_1dindex__(int *ID){
  clear_1DINDEX(*ID);
}
void init_1dindex__(int *ID, int *NDIM, int *NPT_PERDIM ) {
  init_1DINDEX(*ID, *NDIM, NPT_PERDIM);
}
int get_1dindex__(int *ID, int *NDIM, int *indx ) {
  return  get_1DINDEX(*ID, *NDIM, indx );
}
  
// =====================================
void warn_NVAR_KEY(char *fileName) {
  printf("   WARNING: Should remove obsolete NVAR key from %s\n", fileName);
} 

// =====================================
int file_timeDif(char *file1, char *file2) {

  // Mar  29, 2011
  // Return time-stamp difference between file1 and file2.
  // Dif > 0 if  file1 is  more  recent
  // Dif < 0 if  file1 is  older.
  //
  
  char fnam[] = "file_timeDif" ;

  struct stat statbuf1;
  struct stat statbuf2;

  int j1, j2, t1,  t2 ;

  // ------------ BEGIN ----------------

  j1 = stat(file1, &statbuf1); // returns 0 if file exists
  j2 = stat(file2, &statbuf2);

  if ( j1 !=0 ) {
    sprintf(c1err,"istat returned %d  for", j1);
    sprintf(c2err,"file1='%s' ", file1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( j2 !=0 ) {
    sprintf(c1err,"istat returned %d  for", j2);
    sprintf(c2err,"file2='%s' ", file2);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  t1 = statbuf1.st_mtime ;
  t2 = statbuf2.st_mtime ;
  return t1 - t2 ;

} // end of  file_timeDif

// =====================================
int nrow_read(char *file, char *callFun) {

  // Created Jan 19 2016 by R.Kessler
  // read and return number of rows in ascii *file.
  // *calFun is the name of the calling function to 
  // print in case of error.
  // Useful to estimate malloc size.
  //
  // Dec 29 2017: use open_TEXTgz to read gzipped files.

  int NROW = 0, GZIPFLAG ;
  FILE *fp;
  char line[1000];
  char fnam[] = "nrow_read" ;

  fp = open_TEXTgz(file, "rt", &GZIPFLAG );

  if ( fp == NULL ) {
    sprintf(c1err,"Could not open file '%s'", file);
    sprintf(c2err,"Called from function %s", callFun );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  while ( fgets (line, 1000, fp) !=NULL  ) {  NROW++ ; }
  fclose(fp);

  return(NROW);

} // end nrow_read

// =====================================
int rd2columnFile(
		   char *file   // (I) file to read
		   , int MXROW  // (I) max size of column arrays
		   , int *Nrow  // (O) number of rows read from file
		   , double *column1, double *column2 // (O) column data
		   ) {

  // open  *file and read/return 2 data columns.
  // Abort if Nrow exceeds MXROW or if file does not exist.
  // Returns SUCCESS flag upon completion
  // Skips rows that have comments starting with #, !, %
  //
  //
  // Dec 29 2017: use open_TEXTgz() to allow for gzipped file.

  FILE *fp ;
  int  n, GZIPFLAG ;
  double tmp1, tmp2;

  char line[MXPATHLEN], tmpline[MXPATHLEN], s1[40], s2[40] ;
  char *ptrtok ;
  char fnam[] = "rd2columnFile" ;

  // -------------- BEGIN ---------

  fp = open_TEXTgz(file, "rt", &GZIPFLAG );
  if ( fp == NULL ) {
    sprintf(c1err,"Could not open file");
    sprintf(c2err,"%s", file);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  *Nrow =  n = 0;

  while ( fgets (line, 180, fp) !=NULL  ) {

    // skip blank lines
    if ( strlen(line) <= 2 ) { continue ; }
    
    // break of tmpline into blank-separated strings
    sprintf(tmpline,"%s", line);
    ptrtok = strtok(tmpline," ");
    sscanf ( ptrtok, "%s", s1 );
    ptrtok = strtok(NULL, " ");

    // skip comment-lines
    if ( commentchar(s1) == 1 ) { continue ; }

    // get 2nd string 
    sscanf ( ptrtok, "%s", s2 );
    ptrtok = strtok(NULL, " ");

    // strip off the numerical values
    sscanf(s1, "%le" , &tmp1 ) ;
    sscanf(s2, "%le" , &tmp2 ) ;

    n++ ;
    if ( n > MXROW ) {
      sprintf(c1err,"Nrow exceeds MXROW=%d for file", MXROW );
      sprintf(c2err,"%s", file);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    *(column1 + n - 1) = tmp1 ;
    *(column2 + n - 1) = tmp2 ;

  } // end of while loop


  *Nrow = n ;
  fclose(fp);

  return SUCCESS ;

}  // rd2columnFile

// define mangle routine for fortran
int rd2columnfile_(char *file, int *MXROW, int *Nrow,
		   double *col1, double *col2) {
  return rd2columnFile(file,*MXROW, Nrow, col1, col2);
}


// ==========================================
int commentchar(char *str) {

  // returns 1 if input string s is a comment character
  // such as #  %  ! 

  char c1[2];
  sprintf(c1,"%c", str[0]);

  if ( strcmp(c1,"#") == 0 ) return  1 ;
  if ( strcmp(c1,"!") == 0 ) return  1 ;
  if ( strcmp(c1,"%") == 0 ) return  1 ;
  if ( strcmp(c1,"@") == 0 ) return  1 ;

  return  0;
}

// ============================================
void fillbins(int OPT, char *name, int NBIN, float *RANGE, 
	      float *BINSIZE, float *GRIDVAL ) {

  // Created Octg 16, 2010 by RSK
  // For input NBIN and *RANGE (min,max),
  // return *BINSIZE  and *GRIDVAL array
  // OPT=1 : grid goes from MIN to MAX
  // OPT=2 : grid goes from MIN+BINSIZE/2 to MAX-BINSIZE/2
  //         (i.e, bin-centers like HBOOK histograms)
  //
  // Initial use is for sngridtools.c
  // if *name is NOT null, then print each gridval


  float VALMIN, VALMAX, DIF, xi, VAL, BIN, BINOFF ;
  int i, LDMP ;

  char fnam[] = "fillbins" ;

  // --------------- BEGIN --------------

  if ( NBIN <= 0 ) {
    sprintf(c1err,"invalid NBIN=%d for '%s' ", NBIN, name );
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
  }

  VALMIN = *(RANGE+0); 
  VALMAX = *(RANGE+1); 
  DIF    =  VALMAX - VALMIN ;

  if ( OPT == 2 || NBIN == 1 ) {
    BIN    =  DIF / (float)(NBIN);
    BINOFF = 0.5;
  }
  else {
    BIN    =  DIF / (float)(NBIN-1);
    BINOFF = 0.0 ;
  }

  *BINSIZE = BIN ;

  LDMP = strlen(name);
  if ( LDMP ) printf(" %2d GRID-BINS for %8s : ", NBIN, name );


  // now fill grid value at each grid center

  for ( i = 0; i < NBIN; i++ ) {
    xi = (float)i;
    VAL = VALMIN + BIN * (xi + BINOFF) ;
    *(GRIDVAL+i) = VAL;

    if ( LDMP ) printf("%6.2f ", VAL);
  }

  if ( LDMP ) printf("\n");


} // end of fillbins



// ***********************************
double FlatRan ( int ilist, double *range ) {

  // return random number between range[0] and range[1]
  // 'ilist' selects which random list.

  double xran, dif, rantmp;
  double x0, x1;

  xran     = FlatRan1(ilist);
  x0       = range[0] ;
  x1       = range[1] ;
  dif      = x1 - x0;
  rantmp   = x0  + xran * dif ;

  return rantmp;

}  // end of rangen


// ****************************************
double biGaussRan(double siglo, double sighi ) {

  // Return random number from bifurcate gaussian
  // with sigmas = "siglo" and "sighi" and peak = 0.0
  //
  // Jan 2012: always pick random number to keep randoms synced.

  double rr, rg, sigsum, plo, biran    ;

    // ---------- BEGIN ------------


  // pick random number to decide which half of the gaussian we are on
  rr = FlatRan1(1) ;  // pick random number between 0 and 1

  biran = 0.;

  if ( siglo == 0.0 && sighi == 0.0 ) { return biran;  } 

  sigsum = siglo + sighi ;
  plo    = siglo / sigsum ;    // prob of picking lo-side gaussian

  if ( sigsum <= 0.0 ) { return biran;  }


  // pick random gaussian number; force it positive
  rg = GaussRan(1);
  if ( rg < 0.0 ) { rg = -1.0 * rg ; }  // force positive random

  if ( rr < plo ) 
    { biran = -rg * siglo;  }
  else
    { biran = +rg * sighi;  }


  return biran ;

} // end of biguassran

// **********************************************
double skewGaussRan(double xmin, double xmax, 
		    double siglo, double sighi, 
		    double skewlo, double skewhi) {

  //  Mar 16 2014
  //
  // Select random number from skewed Gaussian defined by
  //
  //      SG(x) = exp [ x^2 / (2 * SIG^2) ]
  // where
  //   SIGLO = siglo + |x|*skewlo
  //   SIGHI = sighi + |x|*skewhi
  //
  //
  // Notes:
  //   *  'skew' is not the usual definition of a skewed Gaussian.
  //   *  xmin, xmax are bounds such that peak of distribution is 
  //      xpeak=0; therefore pass xmin = XMIN-XPEAK and xmax = XMAX-XPEAK.  
  //   * xmin must be negative, and xmax must be positive.
  //   * siglo and sighi must both be positive
  //   * skew can be positive or negative, but if SIGLO or SIGHI
  //     become negative, code aborts.
  //
  // Method:
  // Since we cannot do exact solution, select bi_Gaussian random 'x' 
  // from Bounding BiGaussian Function [BBGF] defined with the skewed
  // sigma = |xmax| or |xmin|,  depending on the sign of skew.
  // Then assign probability Prob = SG(x)/BBGF(x) and verify that
  // 0 <= P <= 1. Keep selecting from BBGF until random(0-1) < Prob.
  //

  double SIGHI_BBGF, SIGLO_BBGF, x, P_BBGF, P_SG, Prob, ran1 ;
  char fnam[] = "skewGaussRan" ;
  int NTRY;
  int MXTRY = 100 ;

  // ------------ BEGIN -----------

  // define BBGF = Bounding Bi-Gaussian Function
  if ( skewlo > 0.0 ) 
    { SIGLO_BBGF = fabs(xmin) ; }
  else 
    { SIGLO_BBGF = siglo ; }


  if ( skewhi > 0.0 ) 
    { SIGHI_BBGF = fabs(xmax) ; }
  else 
    { SIGHI_BBGF = sighi ; }


  sprintf(c2err,"BND=%6.3f,%6.3f  sig=%6.3f,%6.3f  skew=%.3f,%.3f", 
	  xmin,xmax, siglo,sighi, skewlo,skewhi );
  


  NTRY = 0 ;
  P_SG = P_BBGF = Prob = ran1 = -9.0 ;

 PICKRAN:
  NTRY++ ;

  /* xxx
  printf(" xxx ------------------------------- \n");
  printf(" xxx NTRY = %d (xmin,xmax=%.3f,%.3f) \n", NTRY, xmin, xmax );
  printf(" xxx x=%.3f  Prob= %.4le/%.4le = %.4le   (ran1=%f) \n",
	 x, P_SG, P_BBGF, Prob, ran1); 
  fflush(stdout);
  xxx */


  if ( NTRY > MXTRY ) {
    sprintf(c1err,"Could not find random from skewed Gaussian after %d tries", 
	    NTRY );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  x = biGaussRan(SIGLO_BBGF, SIGHI_BBGF);
  if ( x < xmin ) { goto PICKRAN ; }
  if ( x > xmax ) { goto PICKRAN ; }

  P_BBGF = skewGauss(x, SIGLO_BBGF, SIGHI_BBGF, skewlo, skewhi );
  P_SG   = skewGauss(x, siglo,      sighi,      skewlo, skewhi );

  if ( P_BBGF <= 1.0E-9 ) {
    sprintf(c1err,"Invalid P_BBGF = %le for x=%.3f", P_BBGF, x);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  Prob = P_SG / P_BBGF ;
    
  if ( Prob < 0.0 || Prob > 1.0 ) {
    sprintf(c1err,"Invalid P_SG/P_BBGF = %.4le/%.4le = %.4f for x=%.3f", 
	    P_SG, P_BBGF, Prob, x);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // apply weight for this x value

  ran1 = FlatRan1(2); // 2-> 2nd list of randoms
  if ( ran1 < Prob ) 
    { return x ; }
  else
    { goto PICKRAN ; }


} // end of skewGaussRan 

double skewGauss(double x, double siglo,double sighi, 
		 double skewlo, double skewhi ) {

  // March 2014
  // See function definition at top of skewGaussRan().
  // Note that function peak value (mode) must correspond to x=0.

  double sqx, sig, sqsig, arg, funval ;
  char fnam[] = "skewGauss";

  // ---------------  BEGIN ----------------

  if ( x < 0.0 ) 
    { sig = siglo + skewlo * fabs(x); }
  else
    { sig = sighi + skewhi * x; }

  if ( sig < 0.0 ) {
    sprintf(c1err, "sig = %f < 0 for x=%f", sig, x);
    sprintf(c2err, "sig(lo,hi)=%.3f,%.3f  skew(lo,hi)=%.3f,%.3f",
	    siglo, sighi, skewlo, skewhi);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  sqx    =   x * x ;
  sqsig  = sig * sig ;
  arg    = 0.5 * sqx/sqsig ;
  funval = exp(-arg);

  return funval ;

} // end of skewGauss



// **********************************************
void init_RANLIST(void) {

  // Dec 1, 2006 RSK
  // Load RANLIST8 array with lots of random numbers (uniform from 0-1).
  // If SEED is fixed, then a fixed list is generated
  // for each SN. [note that the list is different for each SN ...
  // but the lists repeate themselves when the simulation reruns.]
  //
  // Feb 26 2013: fill separate RANLIST8_GENSMEAR list so that
  //              intrinsic-scatter randoms can be changed without
  //              changing synced randoms for main generation.
  //
  // Jun 9 2018: use unix_random() call.

  int ilist, istore;
  char fnam[] = "init_RANLIST" ;

  // ---------------- BEGIN ----------------

  NLIST_RAN = 0 ;

  NLIST_RAN++ ;   // main generation
  NLIST_RAN++ ;   // genSmear
  NLIST_RAN++ ;   // GENSPEC_DRIVER (Jan 2018)


  if ( NLIST_RAN > MXLIST_RAN ) {
    sprintf(c1err,"NLIST_RAN=%d exceeds bound.", NLIST_RAN);
    sprintf(c2err,"Check MXLIST_RAN = %d", MXLIST_RAN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  for (ilist = 1; ilist <= NLIST_RAN; ilist++ ) {
    NSTORE_RAN[ilist] = 0 ;
    for ( istore=0; istore < MXSTORE_RAN; istore++ ) {
      RANSTORE8[ilist][istore] = unix_random();
    }
  }

}  // end of init_RANLIST

// **********************************
double unix_random(void) {
  // Created Jun 9 2018
  // Return random between 0 and 1
  long int i8 = random(); 
  double   r8 = (double)i8 / (double)RAND_MAX ;  // 0 < r8 < 1 
  return(r8);
}
double unix_random__(void) { return( unix_random() ); }

// ***********************************
double GaussRan(int ilist) {
  // return Gaussian random number
  double R,  V1, V2, FAC, G ;
  // --------------- BEGIN ----------------
 BEGIN:
  V1 = 2.0 * FlatRan1(ilist) - 1.0;
  V2 = 2.0 * FlatRan1(ilist) - 1.0;
  R  = V1*V1 + V2*V2 ;
  if ( R >= 1.0 ) { goto BEGIN ; }
  FAC = sqrt(-2.*log(R)/R) ;
  G = V2 * FAC ;

  return G ;
}  // end of Gaussran


double GaussRanClip(int ilist, double ranGmin, double ranGmax ) {
  // Created Aug 2016
  double ranG ;
 PICK_RANGAUSS:
  ranG = GaussRan(ilist);
  if ( ranG < ranGmin || ranG > ranGmax ) { goto PICK_RANGAUSS; }
  return(ranG);

} // end GaussRanClip

// *********************************
double FlatRan1(int ilist) {

  int  N ;
  double   x8;
  char fnam[] = "FlatRan1" ;

  // return random number between 0 and 1
  // Feb 2013: pass argument 'ilist' to pick random list.

  if ( ilist < 1 || ilist > NLIST_RAN ) {
    sprintf(c1err,"Invalid ilist = %d", ilist);
    sprintf(c2err,"Valid ilist is 1 to %d", NLIST_RAN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // check to wrap around with random list.
  if ( NSTORE_RAN[ilist] >= MXSTORE_RAN ) { NSTORE_RAN[ilist] = 0;  }

  // use current random in list
  N  = NSTORE_RAN[ilist] ;
  x8 = RANSTORE8[ilist][N] ;

  // increment for next usage.
  NSTORE_RAN[ilist]++;  

  return x8;

}  // end of FlatRan1


double gaussran_(int *ilist) { return GaussRan(*ilist); }
double flatran1_(int *ilist) { return FlatRan1(*ilist); }


// ********************************************************
double interp_SINFUN(double VAL, double *VALREF, double *FUNREF,
		     char *ABORT_COMMENT) {

  // Created April 11, 2012
  // sin-function interpolation.
  // Interpolate function(val) where VALREF[0] <= val <= VALREF[1]
  // and FUNREF[0-1] are the function values at the edges.
  // Use sin function to interpolate so that derivative=0
  // at the edges.
  // Initial usage is to pass wavelength for 'val' for
  // interpolating functions of intrinsic SN variations.
  //

  double VAL_MEAN, VAL_DIF, FUN_MEAN, FUN_DIF, FUN, ARG, S ;
  double PI = TWOPI/2.0 ;
  char fnam[] = "interp_SINFUN" ;

  // ----------------- BEGIN ----------------

  if ( VAL < VALREF[0] || VAL > VALREF[1] ) {
    sprintf(c1err,"Invalid VAL=%f passed from '%s'", VAL, ABORT_COMMENT);
    sprintf(c2err,"VALREF = %f to %f", VALREF[0], VALREF[1]) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  VAL_MEAN = ( VALREF[1] + VALREF[0] ) * 0.5 ;
  VAL_DIF  = ( VALREF[1] - VALREF[0] );
  FUN_MEAN = ( FUNREF[1] + FUNREF[0] ) * 0.5 ;
  FUN_DIF  = ( FUNREF[1] - FUNREF[0] );
  
  ARG = PI * (VAL - VAL_MEAN)/VAL_DIF ;
  S   = sin(ARG);
  FUN = FUN_MEAN + ( 0.5 * FUN_DIF * S ) ;

  return FUN ;

} // end interp_SINFUN


// **************************************************************
double interp_1DFUN(
		    int OPT         // (I) 1=linear, 2=quadratic interp
		    ,double val     // (I) interp at this point (e.g., lambda)
		    ,int NBIN       // (I) number of function bins  passed
		    ,double *VAL_LIST  // (I) list of grid values
		    ,double *FUN_LIST  // (I) list of function values
		    ,char *abort_comment // (I) comment in case of abort
		    ) {

  // Created April 2011 by R.Kessler
  // 1D interp-utility.
  // Interpolate FUN_LIST(VAL_LIST) at the point *val.
  // Note that the input function need NOT be specified
  // on a uniform grid.
  // 
  // NBIN=3   => do the old interp8   functionality
  // NBIN > 3 => do the old lamInterp functionality


  int  IBIN, ibin0, ibin2 ;
  double 
    val0, val1, frac
    ,fun0, fun1, fun
    ,dif0, dif1
    ,*ptrVAL, *ptrFUN
    ;

  char fnam[] = "interp_1DFUN" ;

  // ------------- BEGIN -----------------

  // do binary search to quickly find which bin contains 'val'
  IBIN = quickBinSearch(val, NBIN,VAL_LIST, abort_comment );

  if ( IBIN < 0 || IBIN >= NBIN-1 ) {
    sprintf(c1err,"quickBinSearch returned invalid IBIN=%d (NBIN=%d)", 
	    IBIN, NBIN );
    sprintf(c1err,"Check '%s'", abort_comment);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  if ( OPT == OPT_INTERP_LINEAR ) {
    val0 = *(VAL_LIST + IBIN ) ;
    val1 = *(VAL_LIST + IBIN + 1) ;
    fun0 = *(FUN_LIST + IBIN ) ;
    fun1 = *(FUN_LIST + IBIN + 1) ;
    frac = (val - val0)/(val1-val0) ;
    fun  = fun0 + frac*(fun1-fun0);
    return fun ;
  }
  else if ( NBIN < 3 ) {
    sprintf(c1err,"Cannot do quadratic interp with %d bins.", NBIN);
    sprintf(c1err,"Check '%s'", abort_comment);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  else if ( OPT == OPT_INTERP_QUADRATIC ) {
    // quadratic interp
    dif0   = val - VAL_LIST[IBIN];
    dif1   = VAL_LIST[IBIN+1] - val;
    
    if ( dif0 <= dif1 && IBIN > 0 )  
      { ibin0 = IBIN - 1; }
    else 
      { ibin0 = IBIN ; }

    // if there are only three bins, then there is no choice
    // for which bins to use.
    if ( NBIN == 3 ) { ibin0 = 0; } 

    // make sure that all three bins are defined
    ibin2 = ibin0+2;
    if ( ibin0 < 0 || ibin2 > NBIN-1 ) {
      sprintf(c1err,"Invalid quadInterp bins: %d to %d", ibin0, ibin2);
      sprintf(c2err,"Must be within defined bins 0 to %d (val=%f)", 
	      NBIN-1, val );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    ptrVAL = &VAL_LIST[ibin0];
    ptrFUN = &FUN_LIST[ibin0];
    fun = quadInterp ( val, ptrVAL, ptrFUN, abort_comment);
    return fun;

  }
  else { 
    sprintf(c1err,"Invalid OPT=%d for '%s'", OPT, abort_comment);
    sprintf(c2err,"%s", "=> Cannot interpolate.") ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    return(-999.);
  }

  return(-999.);

} // end of interp_1DFUN


// **************************************************************
// mangle interp function for fortran.
double interp_1dfun__(int *OPT, double *val, int *NBIN
		    ,double *VAL_LIST, double *FUN_LIST
		    ,char *abort_comment ) {
  return interp_1DFUN(*OPT, *val, *NBIN, VAL_LIST, FUN_LIST, abort_comment);
}

// **************************************************************
double quadInterp ( double VAL, double VAL_LIST[3], double FUN_LIST[3],
		    char *abort_comment ) {
 
/***

  Created April 2011 by R.Kessler, 

  Do quadratic interp using the 3 function values.
  It is recommended that
       a_lam[0] <= lambda  <= a_lam[2]

  or you might get pathological solutions
  
  We solve 

     FUN_LIST[0] = A + B ( VAL_LIST[0] - LBAR )^2
     FUN_LIST[1] = A + B ( VAL_LIST[1] - LBAR )^2
     FUN_LIST[2] = A + B ( VAL_LIST[2] - LBAR )^2

 for A, B and LBAR,
 and then convert to

       quadInterp = C0 + C1*lambda + C2*(lambda**2)

 to avoid divide-by-zero when the 2nd order C2=0

 May 5, 2011:
   no need evaluate C[3]. For a little extra speed,
   compute   FUN = A + B(VAL-LBAR)^2

 Dec 30 2013: fix check on pure linear relation by allowing
              slope of E-12

 ****/

   double 
     A, B, LBAR   // defines quadratic function 
     , v0, v1, v2   // three values of input VAL_LIST
     , f0, f1, f2   // three values of input function 
     , sqv0, sqv1, sqv2
     , dprod, fprod  // internal
     , sqdval10, sqdval21, sqdval02
     , dval10,   dval21,   dval02
     , top, dval, dum
     , slope
     , dif
     , fun            // fun(lambda) => return value 
     , C[3]
        ;


   int LDMP;

   /* ------------------------ BEGIN -------------------- */

   f0 = FUN_LIST[0];
   f1 = FUN_LIST[1];
   f2 = FUN_LIST[2];

   v0 = VAL_LIST[0];
   v1 = VAL_LIST[1];
   v2 = VAL_LIST[2];

   // check trivial case of flat function
   if ( f0 == f1 && f1 == f2 ) {  fun = f1; return fun ; }

   // check for pure linear relation 
   if ( f1 != f0 ) {
     dum = (f2-f1)/(f1-f0) ;
     if ( fabs(dum-1.0) < 1.0E-12 ) {  // f2-f1 \simeq f1 - f0
       dval  = v1 - v0 ;
       slope = (f1-f0) / dval ;
       dval  = VAL - v0 ; 
       fun   = f0 + slope * dval ;       
       return fun;     
     }
   }

   // if we get here, do quadratic interpolation

   sqv0 = v0 * v0 ;
   sqv1 = v1 * v1 ;
   sqv2 = v2 * v2 ;

   sqdval10 = sqv1 - sqv0 ;
   sqdval21 = sqv2 - sqv1 ;
   sqdval02 = sqv0 - sqv2 ;

   dval10   = v1 - v0 ; 
   dval21   = v2 - v1 ; 
   dval02   = v0 - v2 ; 
   dprod = dval21 * dval10 * dval02;


   fprod = f0*dval21 + f1*dval02 + f2*dval10;

   if ( dprod == 0.0 ) 
      B = 0;
   else
      B = -fprod / dprod;


   top = f0 * sqdval21 + f1*sqdval02 + f2*sqdval10;
   if ( fprod == 0.0 )
      LBAR = 0.0 ;
   else
      LBAR = 0.5 * top / fprod;
  
   dval = v0 - LBAR;
   A = f0 - B*dval*dval;

   // compute function value

   dif  = VAL - LBAR ;
   fun  = A + B * (dif*dif) ; 

   /*
   C[0] = A + B * LBAR * LBAR;
   C[1] = top / dprod;
   C[2] = B;
   fun  = C[0] + C[1]*VAL + C[2]*(VAL*VAL) ;
   */

   LDMP = 0 ; // ( VAL == 4040.0 ) ; 
   if ( LDMP ) {
     printf(" xxxx -------- NEW interp_1DFUN DUMP ------------------ \n") ;
     printf(" xxxx v0,v1,v2 = %f %f %f \n", v0, v1, v2 );
     printf(" xxxx f0,f1,f2 = %f %f %f \n", f0, f1, f2 );
     printf(" xxxx quadInterp = %f  at val=%f \n", fun, VAL );
     printf(" xxxx C0, C1, C2 = %f %f %f \n", 
	    C[0], C[1], C[2] );
     printf(" xxx A=%le  B=%le  LBAR=%f \n", A, B, LBAR);
     //     debugexit("quad interp");
   }

   return fun ;

} // end of quadInterp


// ===================================================
int quickBinSearch(double VAL, int NBIN, double *VAL_LIST,
		   char *abort_comment) {

  // April 2011.
  // Return integer bin [0 < IBIN < NBIN-1] such that 
  // *(VAL_LIST+IBIN) contains VAL.
  // Use binary search to quickly find IBIN when NBIN is very large.
  //

  char fnam[] = "quickBinSearch" ;
  int  LDMP, NITER, ibin_min, ibin_max, ibin, ibin1, ibin2, ISTEP ;
  double    MINVAL, MAXVAL, VAL1, VAL2 ;

  // ------------- BEGIN --------------

  LDMP = 0 ; // ( fabs(VAL-7000.) < 0.01 );

  MINVAL = VAL_LIST[0] ;
  MAXVAL = VAL_LIST[NBIN-1];

  if ( VAL < MINVAL || VAL > MAXVAL )  {
    sprintf(c1err,"VAL = %le outside '%s' range",  
	    VAL, abort_comment );
    sprintf(c2err,"defined for %le < VAL < %le", 
	    MINVAL, MAXVAL );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  NITER    = 0;  // for efficiency testing only
  ibin_min = 0; 
  ibin_max = NBIN-1 ;

 NEXTSTEP:

  ISTEP    = (ibin_max-ibin_min)/2 ;  // start with just two lambda bins

  // if we are down to fewer than 6 bins to check,
  // skip binary search and just do brute force search
  // one bin at a time.
  if ( ISTEP < 6 ) { ISTEP = 1; }

  /*
  printf("\t xxxx ibin = %3d to %3d  STEPS of %d \n", 
	 ibin_min, ibin_max, ISTEP);
  */

  for ( ibin=ibin_min; ibin <= ibin_max; ibin+=ISTEP ) {

    if ( ibin >= NBIN - 1 ) { continue ; }

    ibin1 = ibin ;
    ibin2 = ibin + ISTEP;
    if ( ibin2 >= NBIN-1 ) { ibin2 = NBIN-1 ; }

    VAL1 = VAL_LIST[ibin1] ;
    VAL2 = VAL_LIST[ibin2] ;
    NITER++ ;

    // abort if NITER gets larger than NBIN
    if ( NITER > NBIN ) {
      printf("\n PRE-ABORT DUMP: \n");
      printf("\t ibin1=%d     ibin2=%d  (ISTEP=%d) \n", 
	     ibin1,  ibin2, ISTEP );
      printf("\t VAL1/VAL2 =%6.0f/%6.0f A\n", VAL1,  VAL2);
      printf("\t MINVAL/MAXVAL(entire grid) = %6.0f/%6.0f \n", MINVAL, MAXVAL);

      sprintf(c1err,"Could not find '%s' bin after NITER=%d (NBIN=%d)", 
	      abort_comment, NITER, NBIN);
      sprintf(c2err,"Something is wrong here.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    // just do linear interpolation when we find the right bin
    if ( VAL >= VAL1 && VAL <= VAL2 ) {
      if ( ISTEP == 1 ) 
	{ return ibin1 ; }
      else {
	ibin_min = ibin1 ;
	ibin_max = ibin2 ;
	goto NEXTSTEP ;
      }
    
    }  // VAL if-block
  }  // ibin loop


  // if we get here, then something is really wrong.
  sprintf(c1err,"Something is REALLY messed up:");
  sprintf(c2err,"Could not find '%s' bin for VAL=%le", abort_comment, VAL );
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  return(-9) ;

} // end of quickBinSearch


// ************************************************
double polyEval(int N, double *coef, double x) {

  // evaluate polynomial function sum_1^N coef[i] * x^i
  // and avoid using the slow pow function.

  double F, xpow ;
  int i;

  F = 0. ;   xpow = 1.0 ;

  for(i=0; i<N; i++ ) {
    F += coef[i] * xpow ;
    xpow *= x ;
  }

  return F ;

} // end of polyEval


// ***************************************************
double modelmag_extrap(
		       double T          // Time to extrapolate to
		       ,double Tref      // nearest time with model mag
		       ,double magref   // mag at Tref
		       ,double magslope  // dmag/dT  at Tref
		       ,int LDMP         // print to screen if true
		       ) {

  /****
    Created Apr 12, 2009 by R.Kessler

    Extrapolates model magnitide from Tref to time T
    (obs or rest-frame), where T is outside range of model.
    Note that T < Tref assumes on the rise, and T > Tref 
    assumes on the fall.
    
    Aug 7, 2014: fix bug and return extrapolated mag 
  ****/

  double mag, arg ;

  // -------- BEGIN ----------

  if ( T > Tref ) {

    // linear extrap at late times
    mag = magref + (T-Tref)*magslope ;
  }
  else {

    // exponential fall for early times
    arg = (T - Tref) * magslope/magref ;
    mag = magref * exp(arg);
  }

  return mag ; // Aug 7 2014

} // end of modelmag_extrap


// ***************************************************
double modelflux_extrap(
		       double T          // (I) Time to extrapolate to
		       ,double Tref      // (I) nearest time with model flux
		       ,double fluxref   // (I) flux at Tref
		       ,double fluxslope // (I) dflux/dT  at Tref
		       ,int LDMP         // (I) print to screen if true
		       ) {

  /****
    Created Apr 12, 2009 by R.Kessler

    Extrapolates model flux from Tref to time T
    (rest-frame), where T is outside range of model.
    T and Tref are in days, defined so that T=0 at peak brightness.

    T < Tref => 
    Rises like  F0*(T-T0)^2 x (T-T1) where the quadratic term 
    is a simple fireball model, and the 2nd (linear) term allows 
    matching both the function and its derivative at T=Tref.

    T > Tref => on the tail; use linear mag-extrap, which
                is a power-law flux-extrap.  If fluxref < 0,
                then just hard-wire flux=1E-30 since model
                is unphysical and extrapolation makes no sense.
  
  ****/

  double flux, slope, arg, F0, FP, T1 ;
  double DTref0, SQDTref0, CUBEDTref0 ,DT0, DT1 ;
  double T0  = -20.0 ; // time of explosion, wher T=0 at peak

  //  char fnam[] = "modelflux_extrap" ;

  // -------- BEGIN ----------

  // always use same sign to ensure that flux-> 0 at T->infinity
  slope = fabs(fluxslope) ;

  if ( T > Tref ) {

    // linear mag-extrap at late times => power-law extrap for flux
    if ( fluxref > 0.0 ) {
      arg  = -(T-Tref)*slope / (LNTEN * fluxref) ;
      flux = fluxref * pow(TEN,arg) ;
    }
    else { 
      flux = 1.0E-30 ; 
    }
  }
  else {

    if ( T < T0 ) { flux = 0.0 ; return flux ; }

    // solve for F0 and T1 to describe pre-max rise
    FP         = fluxslope ;
    DTref0     = Tref - T0;
    SQDTref0   = DTref0*DTref0;
    CUBEDTref0 = SQDTref0*DTref0;

    F0  = (FP * DTref0 - 2.*fluxref)/CUBEDTref0 ;
    T1  = Tref - fluxref / (F0*SQDTref0) ;

    DT0 = T - T0 ;
    DT1 = T - T1 ;
    flux = F0 * DT0*DT0 * DT1 ;
  }

  if ( LDMP > 0 ) {
    printf("FLUX_EXTRAP: F(T=%5.1f)=%10.3le", T, flux );
    printf(" [Fref(Tref=%3.0f)=%10.3le & dF/dT=%10.3le] \n",
	   Tref, fluxref, fluxslope );
  }

  return flux ;

} // end of modelflux_extrap



// ====================================================
int rd_sedFlux(
	     char *sedFile         // (I) name of SED file to read
	     ,char *sedcomment     // (I) comment to print in front of SED
	     ,double DAYrange[2]   // (I) rest-frame range (days) to accept
	     ,double LAMrange[2]   // (I) lambda-range to accept (A)
	     ,int    MXDAY         // (I) bound on *DAY_LIST
	     ,int    MXLAM         // (I) bound on *LAM_LIST
	     ,int    OPTMASK       // (I) see OPTMASK bits below
	     ,int    *NDAY         // (O) # of epochs within Trange
	     ,double *DAY_LIST     // (O) list of epochs (days)
	     ,double *DAY_STEP     // (O) day-step
	     ,int    *NLAM         // (O) # of lambda bins within LAMrange
	     ,double *LAM_LIST     // (O) list of lambdas
	     ,double *LAM_STEP     // (O) lambda-step
	     ,double *FLUX_LIST    // (O) flux list
	     ,double *FLUXERR_LIST // (O) flux error list
	    ) {

  /***************
   Created Apr 8, 2009 by R.Kessler
   General function to read SED flux with format

      DAY-0   WAVELENGTH-0   FLUX-0
      DAY-0   WAVELENGTH-1   FLUX-1
      DAY-0   WAVELENGTH-2   FLUX-2
      etc ...
 
   Binning in DAY and WAVELENGTH must be uniform, 
   although there is no explicit check here.

   OPTMASK += 1 --> read FLUXERR (4th column of SEDFILE)
   OPTMASK += 2 --> allow non-uniform DAY bins 

   The return-arg lengths are
     - DAY_LIST has length NDAY
     - LAM_LIST has length NLAM and 
     - FLUX_LIST has length NDAY * NLAM

    To extract
    for ( jflux=0; jflux < NDAY * NLAM ; jflux++ ) 
       iep  = jflux/NLAM
       ilam = jflux - NLAM * iep

       Flux   = *(FLUX_LIST+jflux);
       epoch  = *(DAY_LIST+iep);
       lambda = *(LAM_LIST+ilam);


     HISTORY

  Jan 19 2017: if just one bin, return DAY_STEP=0 or LAM_STEP=0.
               Fix divide-by-zero bug.

  Aug 7 2017: little refactoring and abort if DAYSTEP changes.
  Aug 11 2017: 
     + use splitString2 to simplify reading & parsing
     + read 80 chars per line instead of 180 ... hopefully faster
     + pass OPTMASK arg instead of IERRFLAG. See mask def above
  
  Dec 29 2017: use open_TEXTgz() to read gzipped file.

  May 29 2018: 
   + prepare to suppress DAYs in which FLUX=0 at all wavelenghts.
     See NONZEROFLUX counter. Not implemented

  Jan 4 2019:
   +  move error checking earlier, right after loop.
   + add error check for NDAY=NLAM=0 (maybe catch tabs)
        
  **********/

#define  MXWORD_RDFLUX 4  // max values per line to read

  FILE *fpsed;

  char txterr[20], line[200] ;
  //  char *ptrtok, s1[60], s2[60], s3[60], s4[60], tmpline[200] ;
  char *ptrStringVal[4], StringVal[4][40];
  char space[] = " ";
  char fnam[]  = "rd_sedFlux" ;

  double day, lam, day_last, lam_last, lam_expect, flux, fluxerr, XN ;
  double daystep_last, daystep, daystep_dif ;
  double lamstep_last, lamstep, lamstep_dif ;
  int iep, ilam, iflux, ival, LAMFILLED, OKBOUND_LAM, OKBOUND_DAY, NBIN ;
  int NRDLINE, NRDWORD, GZIPFLAG, FIRST_NONZEROFLUX, NONZEROFLUX ;
  int OPT_READ_FLUXERR, OPT_FIX_DAYSTEP, OPT_FIX_LAMSTEP ;

  // define tolerances for binning uniformity (Aug 2017)
  double DAYSTEP_TOL = 0.5E-3; // tolerance on DAYSTEP uniformity
  double LAMSTEP_TOL = 0.01;   // tolerance on LAMSTEP uniformity

  // ------------- BEGIN -------------

  // init return args
  *NDAY = *NLAM = 0 ;

  // set flags that *NDAY and *NLAM are within bounds
  OKBOUND_DAY = OKBOUND_LAM = 1; 

  // open SED file.

  fpsed = open_TEXTgz(sedFile, "rt", &GZIPFLAG );
  if ( fpsed == NULL ) {
    sprintf(c1err,"Cannot open SED file: " );
    sprintf(c2err,"  '%s' ", sedFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // check OPTMASK args
  OPT_READ_FLUXERR = 0  ;              // default: ignore fluxerr column
  OPT_FIX_DAYSTEP = OPT_FIX_LAMSTEP = 1; // default => require fixed bin size

  if ( (OPTMASK & 1) > 0 ) { OPT_READ_FLUXERR = 1; }
  if ( (OPTMASK & 2) > 0 ) { OPT_FIX_DAYSTEP  = 0; }

  if ( OPT_READ_FLUXERR  ) 
    { sprintf(txterr, "and errors"); }
  else
    { txterr[0]=0; }

  printf("  Read  %s  SED %s from : \n    %s \n", 
	 sedcomment, txterr, sedFile );

  fflush(stdout);

  iflux = iep = ilam = 0 ;

  LAMFILLED = NRDLINE = NONZEROFLUX = FIRST_NONZEROFLUX = 0;
  daystep_last = daystep=-9.0;  day_last = -999999. ;
  lamstep_last = lamstep=-9.0;  lam_last = -999999. ;

  for(ival=0; ival < MXWORD_RDFLUX; ival++ ) 
    { ptrStringVal[ival] = StringVal[ival];   }

  while ( fgets (line, 80, fpsed ) != NULL  ) {

    NRDLINE++ ;   

    // skip blank lines
    if ( strlen(line) <= 2 ) { continue ; }
    if ( commentchar(line) ) { continue ; }

    splitString2(line, space, MXWORD_RDFLUX,  // input line is destroyed
		 &NRDWORD, ptrStringVal ) ;  // returned

    sscanf(StringVal[0], "%le" , &day  ) ;
    sscanf(StringVal[1], "%le" , &lam  ) ;
    sscanf(StringVal[2], "%le" , &flux ) ;
    if ( OPT_READ_FLUXERR ) { sscanf(StringVal[3], "%le" , &fluxerr ) ;  }

    if ( day < DAYrange[0] ) { continue ; }
    if ( day > DAYrange[1] ) { continue ; }

    if ( lam < LAMrange[0] ) { continue ; }
    if ( lam > LAMrange[1] ) { continue ; }

    if ( day > day_last ) {
      if ( iep >= MXDAY ) 
	{ OKBOUND_DAY = 0 ; }
      else
	{ DAY_LIST[iep] = day; }
      
      // check that each DAYSTEP is the same (Aug 2017)
      if ( iep>0 && OKBOUND_DAY && OPT_FIX_DAYSTEP ) {
	daystep      = DAY_LIST[iep] - DAY_LIST[iep-1] ;
	daystep_dif  = fabs(daystep - daystep_last);
	if ( iep>1 && daystep_dif > DAYSTEP_TOL ) {
	  sprintf(c1err,"daystep[iep=%d] = (%.4f) - (%.4f) = %le",
		  iep, DAY_LIST[iep], DAY_LIST[iep-1], daystep);
	  sprintf(c2err,"daystep_last   = (%.4f) - (%.4f) = %le",
		  DAY_LIST[iep-1], DAY_LIST[iep-2], daystep_last);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
	}
	daystep_last = daystep ;
      }

      iep++ ; *NDAY = iep ;  NONZEROFLUX = 0 ; 

    /* xxx not now
      // do NOT increment this day if all fluxes are zero.
      // Once a valid flux is found, alway store fluxes, 
      // even if non-zero. Thus only UV region is suppressed if zero.
      if ( NONZEROFLUX || FIRST_NONZEROFLUX ) 
	{ iep++ ; *NDAY = iep ;  NONZEROFLUX = 0 ; }
      else
	{ iflux -= *NLAM; }
    */

    }
    else if ( day < day_last ) {
      sprintf(c1err,"day_last = %f   but day=%f ", day_last, day );
      sprintf(c2err,"New day must increment." ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
    }

    if ( lam < lam_last )  { ilam = 0; LAMFILLED = 1; }
 
    if ( LAMFILLED == 0 ) {

      // increment NLAM and LAM_LIST on first epoch
      if ( ilam >= MXLAM ) 
	{ OKBOUND_LAM = 0 ; }
      else
	{ LAM_LIST[ilam] = lam; }

      *NLAM = ilam+1 ;

      // check that each LAMSTEP is the same (Aug 2017)
      if ( ilam>0 && OKBOUND_LAM ) {
	lamstep      = LAM_LIST[ilam] - LAM_LIST[ilam-1] ;
	lamstep_dif  = fabs(lamstep - lamstep_last);
	if ( ilam>1 && lamstep_dif > LAMSTEP_TOL ) {
	  sprintf(c1err,"lamstep[ilam=%d] = (%.4f) - (%.4f) = %le",
		  ilam, LAM_LIST[ilam], LAM_LIST[ilam-1], lamstep);
	  sprintf(c2err,"lamstep_last   = (%.4f) - (%.4f) = %le",
		  LAM_LIST[ilam-1], LAM_LIST[ilam-2], lamstep_last);
	  //	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
	  errmsg(SEV_WARN, 0, fnam, c1err, c2err ) ;
	}
	lamstep_last = lamstep ;
      }

    } else {

      // make sure that lambdas repeat exactly
      lam_expect = LAM_LIST[ilam] ;
      if ( lam != lam_expect && OKBOUND_LAM ) {
	printf("\n PRE-ABORT info: \n");
	printf("\t DAYrange = %7.2f to %7.2f  (NBIN=%d) \n", 
	       DAYrange[0], DAYrange[1], *NDAY  );
	printf("\t LAMrange = %7.1f to %7.1f  (NBIN=%d) \n", 
	       LAMrange[0], LAMrange[1], *NLAM );
	printf("\t iep=%d  ilam=%d \n", iep, ilam);

	sprintf(c1err,"Found LAM=%6.1f at day = %f ", lam, day);
	sprintf(c2err,"But expected LAM=%7.2f (ilam=%d)", lam_expect, ilam );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
    }

    ilam++ ; 

    // load flux after OKBOUND arrays are set
    if( OKBOUND_DAY && OKBOUND_LAM )  { 
      FLUX_LIST[iflux] = flux ;  
      if ( flux > 1.0E-12 ) { NONZEROFLUX++ ; FIRST_NONZEROFLUX = 1 ;}
      if ( OPT_READ_FLUXERR ) { FLUXERR_LIST[iflux] = fluxerr ;  }
      iflux++ ;  
    }

    day_last = day;
    lam_last = lam;

  }  // end of read-while loop


  // - - - - - - - - - 
  // error checking
  if ( OKBOUND_DAY == 0 ) {
    sprintf(c1err,"NDAY=%d exceeds bound of MXDAY=%d", *NDAY, MXDAY );
    sprintf(c2err,"%s", "Increase EPOCH  binsize or increase MXDAY");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( OKBOUND_LAM == 0 ) {
    sprintf(c1err,"NLAM=%d exceeds bound of MXLAM=%d", *NLAM, MXLAM );
    sprintf(c2err,"%s", "Increase LAMBDA  binsize or increase MXLAM");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // more error checking related to tabs in SED file:
  if ( *NDAY == 0 || *NLAM == 0 ) {
    sprintf(c1err, "Invalid NDAY=%d and NLAM=%d", *NDAY, *NLAM);
    sprintf(c2err, "Check sed file (make sure there are no tabs)");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // - - - - - - - - - 
  NBIN = *NDAY;
  if ( NBIN > 1 ) {
    XN   = (double)(NBIN-1);
    *DAY_STEP = ( DAY_LIST[NBIN-1] - DAY_LIST[0] ) / XN ;
  }
  else 
    { *DAY_STEP = 0.0 ; }


  NBIN = *NLAM ;
  if ( NBIN > 1 ) {
    XN   = (double)(NBIN-1);
    *LAM_STEP = ( LAM_LIST[NBIN-1] - LAM_LIST[0] ) / XN ;
  }
  else
    { *LAM_STEP = 0.0 ; }


  if ( GZIPFLAG==0 ) 
    { fclose(fpsed); } // normal file stream
  else
    { pclose(fpsed); } // for gzip file

  return SUCCESS ;

} // end of rd_sedFlux


// ****************************************************
float effective_aperture ( float PSF_sigma, int VBOSE ) {

  // Returns effective aperture to determine background
  // for PSF defined as Gaussian with sigma = PSF_sigma
  // in pixels.
  // Aperture area = 1 / sum PSF^2
  // PSF = exp(-R^2/2*sigma^2)

  int NPIX, ix, iy ;
  double 
    x, y
    ,arg
    ,PSF, SQPSF
    ,SQR, SQSIG
    ,SUM_PSF, SUM_SQPSF
    ,area
    ,factor
    ;

  // -------------- BEGIN -------------



  // get number of pixels to sum
  NPIX = (int)(20.0 * PSF_sigma) ;
 
  SQSIG = (double)(PSF_sigma * PSF_sigma);

  SUM_PSF       = 0.0 ;
  SUM_SQPSF     = 0.0 ;

  for ( ix = -NPIX; ix <= NPIX; ix++ ) {
    x = (double)ix ;
    for ( iy= -NPIX; iy <= NPIX; iy++ ) {
      y          = (double)iy ;
      SQR        = x*x + y*y ;
      arg        = -0.5 * SQR / SQSIG;
      PSF        = exp(arg) ;
      SQPSF      = PSF * PSF ;
      SUM_PSF   += PSF ;
      SUM_SQPSF += SQPSF ;
    }
  }

  area = SUM_PSF * SUM_PSF / SUM_SQPSF ;

  if ( VBOSE == 1 ) {
    factor = area / (SQSIG * 3.14159265);
    printf(" APERTURE(PSF_sigma=%6.3f) = %10.2f = %8.3f * PI * PSFSIG^2 \n",
	   PSF_sigma, area, factor );
  }

  return area ;

} // end of effective_aperture




// ****************************
double SFR_integral( double H0, double OM, double OL, double W, double Z ) {

  /***
   Integrate SFR(t) from 0 to current time.
   For convenience, integrate  over 'a' instead of redshift
   [since z goes to infinity]
 
        /a
       |   SFR(a') 
     c |  ----------- da'
       |  a' * H(a')
       /0


  ***/

  int    ia, NABIN = 100 ;
  double AMIN, AMAX, ABIN, atmp, ztmp, xa, tmp  ;
  double sum, sfr, aH ;

  double SECONDS_PER_YEAR = 3600. * 24. * 365. ;

  // ---------- BEGIN ------------

  AMIN = 0.0 ;
  AMAX = 1. / (1. + Z) ;
  ABIN = (AMAX - AMIN) / (double)NABIN ;

  sum = 0.0 ;

  for ( ia=1; ia <= NABIN; ia++ ) {
    xa   = (double)ia ;
    atmp = AMIN + ABIN * ( xa - 0.5 ) ;
    ztmp = (1. / atmp) - 1.0 ;

    sfr = SFRfun(H0,ztmp);
    aH  = atmp * Hzfun(H0,OM,OL,W,ztmp) ;
    sum += sfr / aH ;
  }

  // convert H (km/s/Mpc) to H(/year)

  tmp = (1.0E6 * PC_km) / SECONDS_PER_YEAR ;
  sum *= (ABIN * tmp) ;
  return sum ;

}  // end of function SFR_integral



// ****************************
double SFRfun(double H0, double z) {

  /***
      Compute sfr(z) = (a+b*z)/(1+(z/c)^d)*h solar_mass/yr/Mpc^3
      where presumably h = H0/(100 km/s/Mpc)
  ***/

  // -- Baldry and Glazebrook IMF --
  double a = 0.0118 ;
  double b = 0.08 ;
  double c = 3.3 ;
  double d = 5.2 ;

  double tmp1, tmp2, zc, h, SFRLOC;

  // ------------ BEGIN --------------

  zc   = z/c ;
  tmp1 = a + b*z ;
  tmp2 = 1.0 + pow(zc,d) ;
  h    = H0 / 100. ;

  SFRLOC = h * tmp1 / tmp2 ;

  return(SFRLOC) ;

}

// *******************************************
double SFRfun_MD14(double z, double *params) {

  // Created Dec 2016 by R.Kessler
  // use function from  Madau & Dickoson 2014,
  // that was also used in Strolger 2015 for CC rate.
  // This function intended for CC rate, so beware
  // using for other purposes (e.g., no H0 factor here).

  double A = params[0];
  double B = params[1];
  double C = params[2];
  double D = params[3];
  double z1     = 1.0 + z;
  double top    = A*pow(z1,C);
  double bottom = 1.0 + pow( (z1/B), D );

  return( top / bottom );  // also see Eq 8+9 of Strolger 2015          

} // end SFRfun_MD14

// *******************************************
double dVdz_integral
( 
  double H0    // (I) km/s per MPc
  ,double OM    // (I) Omega_matter
  ,double OL    // (I) Omega_lamba
  ,double W     //  (I) w = rho/p
  ,double Zmax  //  (I) integrate up to this redshift
  ,int wgtopt   //  (I) weight integral by z^wgtopt
 ) {

  //
  // return integral of dV/dz = r(z)^2/H(z) dz
  // wgtopt = 0:  returns standard volume integral
  // wgtopt = 1:  returns  z-wgted integral

  double sum, tmp, dz, Ztmp, wz, xz ;
  int NZbin, iz;

  // ---- BEGIN ----------

  // compute exact integral

  NZbin = (int)( Zmax * 1000.0 ) ;
  if ( NZbin < 10 ) { NZbin = 10 ; }
  dz   = Zmax / (float)NZbin ;   // integration binsize
  sum = 0.0;

  for ( iz=1; iz <= NZbin; iz++ ) {
    xz   = (double)iz ;
    Ztmp = dz * (xz - 0.5) ;
    tmp  = dVdz ( H0, OM, OL, W, Ztmp );

    wz = pow(Ztmp,(double)wgtopt);
    sum += wz * tmp;

  }

  sum *= dz ;

  //  printf(" xxxx dVdz_integral = %e (approx=%e) \n", sum, sumtmp  );

  return sum ;

}  // end of dVdz_integral


double dvdz_integral__(double *H0, double *OM, double *OL, double *W,
		       double *Zmax, int *wgtopt) {
  return dVdz_integral(*H0,*OM,*OL,*W,*Zmax,*wgtopt);
}




// **********************************
double dVdz 
( 
  double H0    // (I) km/s per MPc
 ,double OM    // (I) Omega_matter
 ,double OL    // (I) Omega_lamba
 ,double W    //  (I) w = rho/p
 ,double Z    //  (I) redshift
 ) {

  // returns dV/dz = r(z)^2 / H(z)

  double r, H, tmp ;
  double zero = 0.0;

  r = Hzinv_integral ( H0, OM, OL, W, zero, Z );
  H = Hzfun ( H0, OM, OL, W, Z );

  tmp = LIGHT_km * r * r / H ;

  return tmp;

}  // end of dVdz


// ******************************************
double Hzinv_integral 
( 
 double H0     // (I) km/s per MPc
 ,double OM     // (I) Omega_matter
 ,double OL     // (I) Omega_lamba
 ,double W      //  (I) w = rho/p
 ,double Zmin   //  (I) min integ bin
 ,double Zmax   //  (I) integrate up to this redshift
 ) {

  // 
  // Jun 2016: bug fix, (float)NZbin -> (double)NZbin

  int iz, NZbin ;
  double dz, Hz, xz, Ztmp, sum, Hzinv, KAPPA, SQRT_KAPPA ; 

  // ------ return integral c*r(z) = int c*dz/H(z) -------------
  // Note that D_L = (1+z)*Hzinv_integral

  sum = 0.0;

  NZbin = (int)( (Zmax-Zmin) * 1000.0 ) ;
  if ( NZbin < 10 ) { NZbin = 10 ; }
  dz  = (Zmax-Zmin) / (double)NZbin ;      // integration binsize

  for ( iz=1; iz <= NZbin; iz++ ) {
    xz   = (double)iz ;
    Ztmp = Zmin + dz * (xz - 0.5) ;
    Hz   = Hzfun ( H0, OM, OL, W, Ztmp );
    sum += (1.0/Hz) ;
  }

  // remove H0 factor from inetgral before checking curvature.

  sum *= (dz * H0) ;

  // check for curvature
  KAPPA      = 1.0 - OM - OL ; 
  SQRT_KAPPA = sqrt(fabs(KAPPA));

  if ( KAPPA < -0.00001 ) 
    { Hzinv = sin( SQRT_KAPPA * sum ) / SQRT_KAPPA ; }
  else if ( KAPPA > 0.00001 ) 
    { Hzinv = sinh( SQRT_KAPPA * sum ) / SQRT_KAPPA ; }
  else
    { Hzinv = sum ; }


  // return Hzinv with c/H0 factor
  return (Hzinv * LIGHT_km / H0 ) ;

} // end of Hzinv_integral


// ******************************************
double Hainv_integral 
( 
 double H0     // (I) km/s per MPc
 ,double OM     // (I) Omega_matter
 ,double OL     // (I) Omega_lamba
 ,double W      //  (I) w = rho/p
 ,double amin   //  (I) min integ bin
 ,double amax   //  (I) integrate up to this redshift
 ) {

  // May 29, 2008: same as Hzinv_integral, but integrate over a
  // instead of over z.
  // dz/E(z) :  z=1/a-1   dz = -da/a^2

  int ia, Nabin ;
  double da, Hz, xa, atmp, Ztmp, sum, Hzinv, KAPPA, SQRT_KAPPA ; 

  // ------ return integral c*r(z) = int c*dz/H(z) -------------
  // Note that D_L = (1+z)*Hzinv_integral

  sum = 0.0;

  Nabin = (int)( (amax-amin) * 1000.0 ) ;
  if ( Nabin < 10 ) { Nabin = 10 ; }
  da   = (amax-amin) / (double)Nabin ;   // integration binsize

  for ( ia=1; ia <= Nabin; ia++ ) {
    xa   = (double)ia ;
    atmp = amin + da * (xa - 0.5) ;
    Ztmp = 1./atmp - 1.0 ;
    Hz   = Hzfun ( H0, OM, OL, W, Ztmp );
    sum += 1.0/( Hz * atmp * atmp) ;
  }

  // remove H0 factor from inetgral before checking curvature.

  sum *= (da * H0) ;

  // check for curvature
  KAPPA      = 1.0 - OM - OL ; 
  SQRT_KAPPA = sqrt(fabs(KAPPA));

  if ( KAPPA < -0.00001 ) 
    { Hzinv = sin( SQRT_KAPPA * sum ) / SQRT_KAPPA ; }
  else if ( KAPPA > 0.00001 ) 
    { Hzinv = sinh( SQRT_KAPPA * sum ) / SQRT_KAPPA ; }
  else
    { Hzinv = sum ; }


  // return Hzinv with c/H0 factor
  return (Hzinv * LIGHT_km / H0 ) ;

} // end of Hainv_integral


// ******************************************
double Hzfun ( double H0, double OM, double OL, double W, double Z ) {

  //  int iz, NZbin ;
  double sqHz, Hz, ZZ, Z2, Z3, ZL, WW, KAPPA ;

  // ------ returns H(z) -------------

  KAPPA = 1.0 - OM - OL ;  // curvature

  ZZ  = 1.0 + Z ;
  Z2  = ZZ * ZZ ;  // avoid pow fun
  Z3  = Z2 * ZZ ; 

  WW  = 3.0 * (1.0 + W) ;
  ZL  = pow(ZZ,WW) ;

  sqHz = OM*Z3  +  OL*ZL  + KAPPA*Z2 ;

  Hz = H0 * sqrt ( sqHz ) ;

  return Hz ;

} // end of Hzfun


// ******************************************
double dLmag ( double H0, double OM, double OL, double W, 
	       double zCMB, double zHEL ) {

  // returns luminosity distance in mags:
  //   dLmag = 5 * log10(DL/10pc)
  //
  // BEWARE that input H0 is 1/seconds (not km/s/Mpc!)
  // Jan 5 2015: pass zCMB and zHEL
  //
  double rz, dl, arg, mu ;
  double zero = 0.0 ;

  // ----------- BEGIN -----------

  rz     = Hzinv_integral ( H0, OM, OL, W, zero, zCMB );
  dl     = ( 1.0 + zHEL ) * rz ;
  arg    = (double)10.0 * PC_km / dl ;
  mu     = -5.0 * log10( arg );

  return mu ;

}  // end of dLmag

double dlmag_ (double *H0, double *OM, double *OL, double *W,
               double *zCMB, double *zHEL ) {
  double mu = dLmag(*H0,*OM,*OL,*W,*zCMB,*zHEL);
  return(mu);
}


double zcmb_dLmag_invert( double H0, double OM, double OL, double W, double MU ) {

  // Created Jan 4 2018
  // for input distance modulus (MU), solve for zCMB.
  // Beware that H0 unit is km/s/pc (not per Mpc)

  double zCMB, zCMB_start, dmu, DMU, mutmp, DL ;
  double DMU_CONVERGE = 1.0E-4 ;
  int    NITER=0;
  char fnam[] = "zcmb_dLmag_invert" ;

  // ---------- BEGIN ----------

  // use naive Hubble law to estimate zCMB_start
  DL    = pow( 10.0,(MU/5.0) ) * 1.0E-5 ; // Mpc
  zCMB_start = (70.0*DL)/LIGHT_km ;
  zCMB_start *= exp(-zCMB_start/6.0); // very ad-hoc estimate

  zCMB = zCMB_start ;
  DMU = 9999.0 ;
  while ( DMU > DMU_CONVERGE ) {
    mutmp  = dLmag(H0,OM,OL,W, zCMB, zCMB); // MU for trial zCMB
    dmu    = mutmp - MU ;             // error on mu
    DMU    = fabs(dmu);
    zCMB  *= exp(-dmu/2.0); 

    NITER++ ;
    if ( NITER > 500 ) {
      sprintf(c1err,"Could not solve for zCMB after NITER=%d", NITER);
      sprintf(c2err,"MU=%f  dmu=%f  ztmp=%f", MU, dmu, zCMB);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
  } // end dz                                                                   

  int LDMP=0 ;
  if ( LDMP ) {
    printf(" xxx --------------------------------------------- \n");
    printf(" xxx MU=%.4f -> DL = %.2f Mpc  zCMB(start) = %.5f \n",
	   MU, DL, zCMB_start );
    printf(" xxx zCMB -> %.5f  DMU=%f after %d iterations \n",
	   zCMB, DMU, NITER);
  }

  return(zCMB);

} // end zcmb_dLmag_invert


// ***************************************************
void clr_VERSION ( char *version, int prompt ) {

  /*****************
   remove files from 

      $SNDATA_ROOT/photometry/version*
      $SNDATA_ROOT/lcmerge/version*
      $SNDATA_ROOT/SIM/version*

     if prompt = 1 then prompt user first before removing.

   Aug 16, 2007: double char-string memory from 100 to 200.

   Oct 7, 2007: include new _SIM path as well.

   Nov 12, 2007: replace VERSION.LIST with VERSION.* to remove
                 .README file as well as .LIST file

  May 29, 2008: remove entire SIM/VERSION directory to avoid
                'arglist too long' for rm * command.

 Sep 30, 2010: replace gets() with scanf() => no more compile warnings !

 Jan 11 2017: remove user-query to remove old version; just do it.

  ***/

  char fnam[] = "clr_VERSION" ;
  char cmd[200];
  char listFile[200];
  char vprefix[200];
  int  RMFILE_PHOTOMETRY, RMFILE_LCMERGE, RMFILE_SIM, DO_RMFILE, isys ;
  FILE *fp;

  // ----------- BEGIN ------------

  print_banner(fnam);

  // check if this version exists by checking listFile

  DO_RMFILE = 0;  // init remove flag to false
  RMFILE_PHOTOMETRY = 0 ;
  RMFILE_LCMERGE = 0 ;
  RMFILE_SIM     = 0;

  sprintf(listFile, "%s/%s.LIST", PATH_SNDATA_LCMERGE, version);
  if ( (fp = fopen(listFile, "rt"))==NULL ) 
    printf("\t LCMERGE Version %s does not exist. \n", version );
  else { 
    printf("\t LCMERGE Version %s exists. \n", version );
    RMFILE_LCMERGE = 1;  DO_RMFILE = 1 ; fclose(fp); 
  }


  sprintf(listFile, "%s/%s.LIST", PATH_SNDATA_SIM, version);
  if ( (fp = fopen(listFile, "rt"))==NULL ) 
    printf("\t SIM Version %s does not exist. \n", version );
  else {
    printf("\t SIM Version %s exists. \n", version );
    RMFILE_SIM = 1;  DO_RMFILE = 1 ; fclose(fp); 
  }


  sprintf(listFile, "%s/%s.LIST", PATH_SNDATA_PHOTOMETRY, version);
  if ( (fp = fopen(listFile, "rt"))==NULL ) 
    printf("\t PHOTOMETRY Version %s does not exist. \n", version );
  else {
    printf("\t PHOTOMETRY Version %s exists. \n", version );
    RMFILE_PHOTOMETRY = 1; DO_RMFILE = 1 ; fclose(fp); 
  }



  if ( DO_RMFILE == 0 ) return ;


  if ( DO_RMFILE == 1 ) {

    printf("\t Removing version %s files ... ", version );

    // make sure to remove ONLY "version" files and not
    // "version*" files in case different versions have
    // similar names.
    // This means that the .LIST and _XXX files are removed separately.

    if ( RMFILE_PHOTOMETRY == 1 ) {
      sprintf(vprefix,"%s/%s", PATH_SNDATA_PHOTOMETRY, version );
      sprintf(cmd,"rm %s.*  %s_*", vprefix, vprefix );
      fflush(stdout);  isys = system ( cmd );
    }

    if ( RMFILE_LCMERGE == 1 ) {
      sprintf(vprefix,"%s/%s", PATH_SNDATA_LCMERGE, version );
      sprintf(cmd,"rm %s.*  %s_*", vprefix, vprefix );
      fflush(stdout);    
      isys = system ( cmd );
    }

    if ( RMFILE_SIM == 1 ) {
      sprintf(cmd,"rm -rf  %s/ ", PATH_SNDATA_SIM );
      printf("\n execute command: %s \n", cmd);
      fflush(stdout);    
      isys = system ( cmd );
    }

    printf(" Done. \n" );
  }

}  // end of clr_VERSION


// ***************************************************
int  init_VERSION ( char *version ) {

  // init VERSION_INFO structure

  // -------------- BEGIN ----------------

  sprintf(VERSION_INFO.NAME, "%s" , version );

  VERSION_INFO.N_SNLC     = 0;
  VERSION_INFO.N_SNFILE   = 0;
  VERSION_INFO.NEPOCH_TOT = 0;

  return SUCCESS ;

}  // end of init_VERSION



// ********************************************
void ld_null (float *ptr, float value) {

  // load *ptr with "value" if *ptr = NULLFLOAT ;
  // If *ptr is not null, then leave it with current value.

  if ( *ptr == NULLFLOAT ) *ptr = value ;

}




// ********************************************
int  init_SNPATH(void) {

  // Feb 19 2015: 
  //   checkg getenv() == NULL rather than if it returns "(null)".
  //   Set SURVEYNAME to "" instead of SDSS.

  char fnam[] = "init_SNPATH" ;

  // ------------ BEGIN ----------


  if ( getenv("SNDATA_ROOT") == NULL ) {
    sprintf(c1err,"Environment variable SNDATA_ROOT not defined.");
    sprintf(c2err,"See SNANA installation instructions.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  sprintf( PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );
  printf (" SNDATA_ROOT = %s \n", PATH_SNDATA_ROOT);


  // ENV is defined, but now check that the directory is really there.
  // Abort if $SNDATA_ROOT/ does not exist
  int istat;
  struct stat statbuf ;
  istat = stat( PATH_SNDATA_ROOT, &statbuf);
  if ( istat != 0 ) {
    sprintf(c1err,"$SNDATA_ROOT = '%s'", PATH_SNDATA_ROOT);
    sprintf(c2err,"does not exist ?!?!? Check your $SNDATA_ROOT .");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( getenv("SNANA_DIR") == NULL ) {
    sprintf(c1err,"Environment variable SNANA_DIR not defined.");
    sprintf(c2err,"See SNANA installation instructions.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  sprintf( PATH_SNANA_DIR, "%s", getenv("SNANA_DIR") );
  printf (" SNANA_DIR   = %s \n", PATH_SNANA_DIR);

  // - - - - -

  sprintf(PATH_SNDATA_PHOTOMETRY, "%s/photometry", PATH_SNDATA_ROOT);
  sprintf(PATH_SNDATA_LCMERGE,    "%s/lcmerge",    PATH_SNDATA_ROOT);
  sprintf(PATH_SNDATA_SIM,        "%s/SIM",        PATH_SNDATA_ROOT);
  SNDATA.SURVEY_NAME[0]=0;
  SNDATA.SUBSURVEY_NAME[0] = 0 ;

  fflush(stdout);
  return(SUCCESS);

}   // end of init_SNPATH


// *******************************************
int init_SNDATA ( void ) {

  // initialize SNDATA.xxxx elements
  // Note that only one SN index is initialized per call.
  //
  int i_epoch, ifilt, i, igal ;
  // --------- BEGIN -----------------

  sprintf(FLUXUNIT, "ADU");

  sprintf(SNDATA.IAUC_NAME,      "UNKNOWN" );
  sprintf(SNDATA.AUXHEADER_FILE, "UNKNOWN" );

  SNDATA.NLINES_AUXHEADER = 0;
  for ( i=0; i < 20; i++ ) 
    { sprintf(SNDATA.AUXHEADER_LINES[i],"GARBAGE AUXHEADER_LINE"); }

  SNDATA.RA     = NULLFLOAT ;
  SNDATA.DEC    = NULLFLOAT ;
  SNDATA.FAKE   = NULLINT ;
  SNDATA.MWEBV  = NULLFLOAT ;
  SNDATA.WRFLAG_BLINDTEST = 0 ; 
  SNDATA.SNTYPE = -999;

  SNDATA.NEPOCH = 0;
  SNDATA.NEWMJD = 0;


  // default mag settings (1/25/2007)
  sprintf(SNDATA.MAGTYPE, "LOG10");
  sprintf(SNDATA.MAGREF,  "AB");

  // init redshift info

  SNDATA.REDSHIFT_HELIO        = NULLFLOAT ;
  SNDATA.REDSHIFT_HELIO_ERR    = NULLFLOAT ;
  SNDATA.REDSHIFT_FINAL        = NULLFLOAT ;
  SNDATA.REDSHIFT_FINAL_ERR    = NULLFLOAT ;
  SNDATA.VPEC = SNDATA.VPEC_ERR = 0.0 ;

  // init HOSTGAL info
  SNDATA.HOSTGAL_NMATCH[0] = 0;
  SNDATA.HOSTGAL_NMATCH[1] = 0;
  for(igal=0; igal<MXHOSTGAL; igal++ ) {  
    SNDATA.HOSTGAL_OBJID[igal]       = 0;
    SNDATA.HOSTGAL_PHOTOZ[igal]      = -9.0 ;
    SNDATA.HOSTGAL_PHOTOZ_ERR[igal]  = -9.0 ;
    SNDATA.HOSTGAL_SNSEP[igal]       = -9.0 ;
    SNDATA.HOSTGAL_RA[igal]          = -999.0 ;
    SNDATA.HOSTGAL_DEC[igal]         = -999.0 ;
    SNDATA.HOSTGAL_DDLR[igal]        = -9.0 ;
    SNDATA.HOSTGAL_LOGMASS[igal]     = -9.0 ;
    SNDATA.HOSTGAL_LOGMASS_ERR[igal] = -9.0 ;
    SNDATA.HOSTGAL_sSFR[igal]        = -9.0 ;
    SNDATA.HOSTGAL_sSFR_ERR[igal]    = -9.0 ;
  }
  SNDATA.HOSTGAL_USEMASK = 0 ;


  // init SEARCH parameters
  SNDATA.SEARCH_TYPE     = NULLINT ;
  SNDATA.SEARCH_PEAKMJD  = NULLFLOAT ;
  
  // init sim parameters (used for simulation only)
  sprintf(SNDATA.SIM_MODEL_NAME, "NULL" );
  sprintf(SNDATA.SIM_COMMENT, "NULL" );
  SNDATA.SIM_MODEL_INDEX    = NULLINT ; // model class; e.g., SIMSED
  SNDATA.SIM_TEMPLATE_INDEX = NULLINT ; // specific template or SED

  SNDATA.SIM_NOBS_UNDEFINED = 0 ;
  SNDATA.SIM_SEARCHEFF_MASK = 0 ;
  SNDATA.SIM_LIBID      = -9 ;
  SNDATA.SIM_NGEN_LIBID =  0 ;
  SNDATA.SIM_SL_FLAG    = 0 ;
  SNDATA.SIMLIB_FILE[0] = 0 ;
  SNDATA.SIMLIB_MSKOPT  = 0 ;

  SNDATA.SIM_REDSHIFT_HELIO = NULLFLOAT ;
  SNDATA.SIM_REDSHIFT_CMB   = NULLFLOAT ;
  SNDATA.SIM_REDSHIFT_HOST  = NULLFLOAT ;
  SNDATA.SIM_REDSHIFT_FLAG        = -9 ;
  SNDATA.SIM_REDSHIFT_COMMENT[0]  =  0 ;

  SNDATA.SIM_DLMU     = NULLFLOAT ;
  SNDATA.SIM_RA       = NULLFLOAT ;
  SNDATA.SIM_DEC      = NULLFLOAT ;
  SNDATA.SIM_PEAKMJD  = NULLFLOAT ;
  SNDATA.SIM_AVTAU    = NULLFLOAT ;
  SNDATA.SIM_AV       = NULLFLOAT ;
  SNDATA.SIM_RV       = NULLFLOAT ;
  SNDATA.SIM_STRETCH  = NULLFLOAT ;
  SNDATA.SIM_DM15     = NULLFLOAT ;
  SNDATA.SIM_DELTA    = NULLFLOAT ;

  SNDATA.SIM_MWRV    = NULLFLOAT ;
  SNDATA.SIM_MWEBV   = NULLFLOAT ;
  SNDATA.SIMOPT_MWCOLORLAW = NULLINT ;
  SNDATA.SIMOPT_MWEBV      = NULLINT ;

  SNDATA.SIM_SALT2alpha = NULLFLOAT ;
  SNDATA.SIM_SALT2beta  = NULLFLOAT ;
  SNDATA.SIM_SALT2x0  = NULLFLOAT ;
  SNDATA.SIM_SALT2x1  = NULLFLOAT ;
  SNDATA.SIM_SALT2c   = NULLFLOAT ;
  SNDATA.SIM_SALT2mB  = NULLFLOAT ;

  SNDATA.SIM_MAGSMEAR_COH = 0.0 ;

  SNDATA.SIM_RISETIME_SHIFT = 0.0 ;
  SNDATA.SIM_FALLTIME_SHIFT = 0.0 ;

  SNDATA.SIM_TRESTMIN = 0.0 ;
  SNDATA.SIM_TRESTMAX = 0.0 ;
  SNDATA.SIMFLAG_COVMAT_SCATTER = 0 ;

  SNDATA.NPAR_SIMSED = 0;
  SNDATA.NPAR_LCLIB  = 0;

  SNDATA_FILTER.NDEF = 0;

  SNDATA.NVAR_PRIVATE = 0 ;

  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
    SNDATA_FILTER.MAP[ifilt] = 0;
    SNDATA.SIM_PEAKMAG[ifilt]      = NULLFLOAT ;
    SNDATA.SIM_TEMPLATEMAG[ifilt]  = NULLFLOAT ;
    SNDATA.SIM_GALFRAC[ifilt]      = NULLFLOAT ;
    SNDATA.NPRESN[ifilt]                = NULLINT ;
    SNDATA.HOSTGAL_SB_FLUX[ifilt]       = NULLFLOAT ;
    SNDATA.HOSTGAL_SB_FLUXERR[ifilt]    = NULLFLOAT ;
    SNDATA.HOSTGAL_SB_MAG[ifilt]        = 99.0 ;

    for(igal=0; igal<MXHOSTGAL; igal++ ) {
      SNDATA.HOSTGAL_MAG[igal][ifilt]        = 999.0     ;
      SNDATA.HOSTGAL_MAGERR[igal][ifilt]     = 999.0     ;
    }
    SNDATA.SIM_EXPOSURE_TIME[ifilt]     = 1.0  ;
  }

  SNDATA.PIXSIZE     = NULLFLOAT ;
  SNDATA.NXPIX       = -9 ;  
  SNDATA.NYPIX       = -9 ;  
  SNDATA.CCDNUM[0]   = -9 ;
  SNDATA.CCDNUM[1]   = -9 ;

  SNDATA.SUBSAMPLE_INDEX = -9 ;
  SNDATA.MASK_FLUXCOR    =  0 ;

  //  -------------------------------------------
  // epoch info

  for ( i_epoch = 0; i_epoch < MXEPOCH; i_epoch++ ) {

    SNDATA.QMASK[i_epoch]        = NULLINT ;
    SNDATA.SEARCH_RUN[i_epoch]   = NULLINT ;
    SNDATA.TEMPLATE_RUN[i_epoch] = NULLINT ;

    SNDATA.MJD[i_epoch]          = (double)NULLFLOAT ;

    // xxx mark delete     SNDATA.IDCCD[i_epoch]        = NULLINT ;
    sprintf ( SNDATA.FIELDNAME[i_epoch], "NULL" );

    SNDATA.IDTEL[i_epoch] = NULLINT ;
    sprintf(SNDATA.TELESCOPE[i_epoch], "BLANK" );
    sprintf(SNDATA.DATE[i_epoch],      "BLANK" );
    SNDATA.IDATE[i_epoch]          = NULLINT ;

    SNDATA.FILTINDX[i_epoch]       = NULLINT ;
    SNDATA.SEARCH_FIELD[i_epoch]   = NULLINT ;
    SNDATA.TEMPLATE_FIELD[i_epoch] = NULLINT ;

    SNDATA.GAIN[i_epoch]           = NULLFLOAT ;
    SNDATA.READNOISE[i_epoch]      = NULLFLOAT ;

    SNDATA.XPIX[i_epoch]         = NULLFLOAT ;
    SNDATA.YPIX[i_epoch]         = NULLFLOAT ;
    SNDATA.EDGEDIST[i_epoch]     = NULLFLOAT ;

    SNDATA.SKY_SIG[i_epoch]      = NULLFLOAT ;
    SNDATA.SKY_SIG_T[i_epoch]    = NULLFLOAT ;
    SNDATA.PSF_SIG1[i_epoch]     = NULLFLOAT ;
    SNDATA.PSF_SIG2[i_epoch]     = NULLFLOAT ;
    SNDATA.PSF_RATIO[i_epoch]    = NULLFLOAT ;

      
    SNDATA.FLUXCAL[i_epoch]         = NULLFLOAT ;
    SNDATA.FLUXCAL_ERRTOT[i_epoch]  = NULLFLOAT ;

    SNDATA.MAG[i_epoch]           = NULLFLOAT ;
    SNDATA.MAG_ERRPLUS[i_epoch]   = NULLFLOAT ;
    SNDATA.MAG_ERRMINUS[i_epoch]  = NULLFLOAT ;

    SNDATA.SKYSUB_ERR[i_epoch]    = NULLFLOAT ;
    SNDATA.GALSUB_ERR[i_epoch]    = NULLFLOAT ;

    SNDATA.ZEROPT[i_epoch]         = NULLFLOAT ;
    SNDATA.ZEROPT_ERR[i_epoch]     = NULLFLOAT ;
    SNDATA.ZEROPT_SIG[i_epoch]     = NULLFLOAT ;

    SNDATA.PHOTFLAG[i_epoch]       = 0   ;
    SNDATA.PHOTPROB[i_epoch]       = 0.0 ;

    sprintf(SNDATA.SIMEPOCH_WARPCOLNAM[i_epoch],"NULL");
    sprintf(SNDATA.SIMEPOCH_KCORNAM[i_epoch],"NULL");
    SNDATA.SIMEPOCH_MAGSMEAR[i_epoch] = 0.0 ;

  }  //  end i_epoch init loop

  return SUCCESS ;

}   // end of init_SNDATA



// ********************************************
int wr_SNDATA ( int IFLAG_WR, int IFLAG_DBUG  ) {

  /*******
    Created Mar 8, 2006  R.Kessler
    write out all info to file for this "isn".
    This output file is used for analysis.

    cid_debug => dump into to screen for this cid

    wrflag = 
           = WRITE_MASK_LCMERGE    => write everything
           = WRITE_MASK_SIM_SNANA  => write to /SIM area instead

    The _PHOTOMETRY flag is used to write input photometry files;
    the _LCMERGE flag is used by sndata_build to build merged
    files for analysis.

    vbose = 1 => print one-line statement to screen 
    vbose = 0 => no screen dump

   ================
  Mar 22, 2011: write  SIM_GALFRAC(ifilt_obs)

  Mar 28, 2011: for list of filter-dependent quanities,
                put <CR> after 10 variables to prevent char-overflow
                for parsing routines. Fixes problem with 56 filters.
                See NTMP usage.

  Apr 27, 2011: fprintf SNDATA.ORIGIN_SEARCH_PEAKMJD  right after
                SEARCH_PEAKMJD

  Jul 27, 2011: write SIM_COVMAT_SCATTER info (if used)

  Mar 08, 2012: if SNDATA.HOSTGAL_SBFLAG is set then write
                surface brightness for data or MC.

  Jun 13, 2012: write optional AUXHEADER_LINES instead of BOSS-specific info

  Dec 17, 2012: use filtlist = SNDATA_FILTER.LIST instead of computing 
                filtlist from integer list.

  Feb 3 2014: for sim, always write redshift info, even if -9.

  Feb 7, 2014: replace legacy SEARCH_PEAKMJD key with PEAKMJD.
  Feb 12, 2014: write new SIM_HOSTLIB keys 

  Sep 08 2017:  for LCLIB model, write SNDATA.SIM_TEMPLATEMAG

  May 23 2019: write SIM_MAGSHIFT_HOSTCOR

  **************/

  int 
    cid, NEWMJD, inext, imjd
    ,EPMIN, EPMAX, iep
    ,wstat, istat, fake
    ,ifilt,ifilt_obs, ifilt_start
    ,NFILT_TOT, NFILT_EP
    ,LWRITE_FINAL, LWRITE_SIMFLAG
    ,JTMP, ipar, NTMP, iscat
    ;


  char 
    outfile[MXPATHLEN]
    ,filtlist[MXFILTINDX]
    ,ctmp[100]
    ;

  FILE *fp;

  int   *iptr;
  float *fptr;
  char  *cptr;
  float Zspec, Zspec_err, FTMP ;

  double mjd ;
  char fnam[] = "wr_SNDATA";
  char blank[2] = " " ;

  // ------------- BEGIN -----------------

  LWRITE_FINAL = 0 ;
  if ( (IFLAG_WR & WRITE_MASK_LCMERGE)>0 || 
       (IFLAG_WR & WRITE_MASK_SIM_SNANA)>0 )  {
    sprintf(outfile , "%s", SNDATA.SNFILE_OUTPUT );
    LWRITE_FINAL = 1 ;
  }


    
  else {
    sprintf(c1err,"IFLAG_WR=%d is invalid",  IFLAG_WR );
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
  }

  if ( SNDATA.FAKE > 0 && SNDATA.WRFLAG_BLINDTEST == 0 ) 
    {  LWRITE_SIMFLAG = 1 ; }
  else
    {  LWRITE_SIMFLAG = 0 ; }


  cid    = SNDATA.CID ;
  NEWMJD = SNDATA.NEWMJD;


  if ( (fp = fopen(outfile, "wt"))==NULL ) {
    sprintf(c1err,"Cannot open output file: " );
    sprintf(c2err,"%s", outfile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    fclose(fp);  return ERROR;
  }

  if ( IFLAG_DBUG > 0 )  
    { printf( " WRITE %s\n", outfile); fflush(stdout); }

  // sort epochs in order of MJD
  sort_epochs_bymjd();

  // construct filterlist string.
  sprintf(filtlist, "%s", SNDATA_FILTER.LIST);



  NFILT_TOT   = SNDATA_FILTER.NDEF; 
  ifilt_start = SNDATA_FILTER.MAP[0] ;  // start filter index for data

  // first spit back the input

  if ( cid == IFLAG_DBUG ) 
    { printf(" xxxx %s : write HEADER info for  CID=%d \n", fnam, cid ); }

  if ( strlen(SNDATA.SUBSURVEY_NAME) == 0 ) {
    fprintf(fp, "SURVEY: %s   \n", SNDATA.SURVEY_NAME  ); 
  }
  else  { 
    fprintf(fp, "SURVEY: %s(%s)   \n", 
	    SNDATA.SURVEY_NAME, SNDATA.SUBSURVEY_NAME  ); 
  }

  fprintf(fp, "SNID:   %d   \n", cid  );


  if ( LWRITE_FINAL == 1  )
    { fprintf(fp, "IAUC:    %s \n", SNDATA.IAUC_NAME ); }

  if( SNDATA.WRFLAG_BLINDTEST == 0 ) {
    fprintf(fp, "PHOTOMETRY_VERSION: %s \n", VERSION_INFO.NAME );    
  }


  fprintf(fp, "SNTYPE:  %d \n", SNDATA.SNTYPE ) ;
  fprintf(fp, "FILTERS: %s \n", filtlist );
  fprintf(fp, "RA:      %f  deg \n", SNDATA.RA );
  fprintf(fp, "DEC:     %f  deg \n", SNDATA.DEC );
  fprintf(fp, "PIXSIZE: %7.4f  arcsec \n", SNDATA.PIXSIZE );

  if ( SNDATA.CCDNUM[0] >=0 ) 
    { fprintf(fp, "CCDNUM: %d  \n", SNDATA.CCDNUM[0] ); }

  if ( NEWMJD > 0 ) { fprintf(fp, "NEPOCH:  %d \n", NEWMJD ); }


  // now write extra info detrermined from this program

  fake = SNDATA.FAKE ;
  fprintf(fp, "FAKE:    %d   ", fake );
  if ( fake == FAKEFLAG_DATA ) 
    { fprintf(fp,"(=> data) \n" ); }
  else if ( fake == FAKEFLAG_LCSIM ) 
    { fprintf(fp,"(=> simulated LC with snlc_sim.exe) \n" ); }
  else if ( fake == FAKEFLAG_LCSIM_BLINDTEST ) 
    { fprintf(fp,"(=> BLIND-TEST simulation) \n" ); }
  else
    { fprintf(fp,"(=> unknown data source) \n" ); }


  fptr = &SNDATA.MWEBV ;
  fprintf(fp, "MWEBV:      %7.4f    MW E(B-V) \n", *fptr );

  fptr = &SNDATA.MWEBV_ERR ;
  fprintf(fp, "MWEBV_ERR:  %7.4f    error on MWEBV \n", *fptr );

  if ( SNDATA.APPLYFLAG_MWEBV ) {
    fprintf(fp, "MWEBV_APPLYFLAG: %d    # FLUXCAL corrected for MWEBV \n",
	   SNDATA.APPLYFLAG_MWEBV ) ;
  }

  // write things for real data only, and only if values are set

  if ( SNDATA.FAKE == 0 ) {
    iptr = &SNDATA.NPRESN[ifilt_start] ;
    if ( *iptr != NULLINT ) 
      istat = wr_filtband_int ( fp, "NEPOCH_PRESN:",  
				NFILT_TOT, iptr, filtlist, 0 );
  }


  // ----------------
  // write redshift info for all surveys

  wstat = WRSTAT ( IFLAG_WR, SNDATA.REDSHIFT_FINAL );

  // fill REDSHIFT_SPEC if SN is typed
  if ( SNDATA.SNTYPE > 0 ) {
    Zspec     = SNDATA.REDSHIFT_FINAL; 
    Zspec_err = SNDATA.REDSHIFT_FINAL_ERR;
  }
  else {
    Zspec     = -9.0 ; 
    Zspec_err =  9.0;
  }

  if( SNDATA.WRFLAG_BLINDTEST ) {
    wstat = 0;
    fprintf(fp,"REDSHIFT_SPEC:    %.6f +- %.6f  \n", Zspec, Zspec_err );
  }


  if ( wstat == 1 || LWRITE_SIMFLAG ) {

    fprintf(fp,"REDSHIFT_HELIO:   %.6f +- %.6f  (Helio, z_best) \n"
	    ,SNDATA.REDSHIFT_HELIO, SNDATA.REDSHIFT_HELIO_ERR );

    fprintf(fp,"REDSHIFT_FINAL:   %.6f +- %.6f  (CMB) \n"
	    ,SNDATA.REDSHIFT_FINAL, SNDATA.REDSHIFT_FINAL_ERR );
  
    fprintf(fp,"\n");
  }

  fprintf(fp,"VPEC:      %7.2f  # v_pec correction\n", SNDATA.VPEC );
  fprintf(fp,"VPEC_ERR:  %7.2f  # error on above  \n", SNDATA.VPEC_ERR );

  // ---------------------------------------
  FTMP = SNDATA.SEARCH_PEAKMJD ; 
  if ( FTMP != NULLFLOAT ) { fprintf(fp, "PEAKMJD:  %10.3f \n", FTMP );  }

  
  fprintf(fp, "\n" ); 


  // ---------------------------------------
  // write host-gal info
  wr_HOSTGAL(fp) ;

  // ---------------------------------------
  
  // merge additional header info ... likely from host-galaxy file
  if ( LWRITE_FINAL == 1 ) {
    fprintf(fp, " \n" );
    header_merge( fp, SNDATA.AUXHEADER_FILE );
  } 

 

  // --------------------------------------------------
  // write optional AUXHEADER_LINES (June 2012)
  // at end of header, but before the SIM_XXX keys
  int j ;
  for(j=0; j < SNDATA.NLINES_AUXHEADER; j++ ) 
    {  fprintf(fp,"%s \n", SNDATA.AUXHEADER_LINES[j]); }


  // --------------------------------------------------
  fflush(fp);

  if ( SNDATA.WRFLAG_BLINDTEST ) { goto START_EPOCHS ; }
  
  // write SIM_XXX variables if FAKE > 0
  if ( LWRITE_SIMFLAG > 0  ) {
    fprintf(fp, "\n" );

    fprintf(fp, "SIM_MODEL_NAME:  %s \n", 
	    SNDATA.SIM_MODEL_NAME  ) ;

    fprintf(fp, "SIM_MODEL_INDEX:  %d \n", 
	    SNDATA.SIM_MODEL_INDEX  ) ;

    fprintf(fp, "SIM_TYPE_INDEX:   %d \n",
	    SNDATA.SIM_TYPE_INDEX );

    fprintf(fp, "SIM_TYPE_NAME:    %s \n",
	    SNDATA.SIM_TYPE_NAME );

    iptr = &SNDATA.SIM_TEMPLATE_INDEX ;
    //      fprintf(fp, "SIM_NON1a:      %d   (NONIA index) \n", *iptr ) ;
    fprintf(fp, "SIM_TEMPLATE_INDEX: %d   \n", *iptr ) ;

    fprintf(fp, "SIM_COMMENT:  %s  \n", SNDATA.SIM_COMMENT  ) ;

    iptr = &SNDATA.SIM_LIBID ; 
    fprintf(fp, "SIM_LIBID:  %d  \n", *iptr ) ;

    iptr = &SNDATA.SIM_NGEN_LIBID ;   // added Dec 2015
    fprintf(fp, "SIM_NGEN_LIBID:  %d  \n", *iptr ) ;

    if ( SNDATA.SIMLIB_MSKOPT != 0 ) {
      iptr = &SNDATA.SIMLIB_MSKOPT ;
      fprintf(fp, "SIMLIB_MSKOPT:  %d  \n", *iptr ) ;
    }

    iptr = &SNDATA.SIM_NOBS_UNDEFINED ;   // March 2017
    fprintf(fp, "SIM_NOBS_UNDEFINED:  %d  \n", *iptr ) ;

    fptr = &SNDATA.SIM_REDSHIFT_HELIO ; 
    fprintf(fp, "SIM_REDSHIFT_HELIO:  %.5f  \n", *fptr ) ;

    fptr = &SNDATA.SIM_REDSHIFT_CMB ; 
    fprintf(fp, "SIM_REDSHIFT_CMB:    %.5f  \n", *fptr ) ;

    fptr = &SNDATA.SIM_REDSHIFT_HOST ; 
    fprintf(fp, "SIM_REDSHIFT_HOST:   %.5f  \n", *fptr ) ;

    iptr = &SNDATA.SIM_REDSHIFT_FLAG ;
    cptr = SNDATA.SIM_REDSHIFT_COMMENT;
    fprintf(fp,"SIM_REDSHIFT_FLAG:   %d  # %s\n", *iptr, cptr);

    fptr = &SNDATA.SIM_VPEC ; 
    fprintf(fp, "SIM_VPEC:  %.1f (km/sec) \n", *fptr ) ;

    // - - - -  SIM_HOSTLIB_XXX  (Feb 2014)

    char key[60];
    int NPAR = SNDATA.NPAR_SIM_HOSTLIB; 
    fprintf(fp, "SIM_HOSTLIB_NPAR: %d \n", NPAR);
    for(ipar=0; ipar < NPAR; ipar++ ) {
      sprintf(key,"%s:", SNDATA.SIM_HOSTLIB_KEYWORD[ipar] );
      fprintf(fp, "%-28.28s  %.3f \n", key, SNDATA.SIM_HOSTLIB_PARVAL[ipar] );
    }
    
    // - - - - 

    fptr = &SNDATA.SIM_DLMU ;
    fprintf(fp, "SIM_DLMU:      %.4f  mag   [ -5*log10(10pc/dL) ]\n", *fptr);

    fptr = &SNDATA.SIM_LENSDMU ;
    fprintf(fp, "SIM_LENSDMU:   %.4f  mag \n", *fptr ) ;

    fptr = &SNDATA.SIM_RA ; 
    fprintf(fp, "SIM_RA:        %f deg  \n", *fptr ) ;
    fptr = &SNDATA.SIM_DEC ; 
    fprintf(fp, "SIM_DEC:       %f deg  \n", *fptr ) ;

    fptr = &SNDATA.SIM_MWRV ; 
    fprintf(fp, "SIM_MWRV:   %6.3f   (MilkyWay RV) \n", *fptr ) ;
    fptr = &SNDATA.SIM_MWEBV ; 
    fprintf(fp, "SIM_MWEBV:  %7.4f   (MilkyWay E(B-V)) \n", *fptr ) ;

    iptr = &SNDATA.SIMOPT_MWCOLORLAW ; 
    fprintf(fp, "SIMOPT_MWCOLORLAW:  %d   (MW color-law option) \n",  *iptr);

    iptr = &SNDATA.SIMOPT_MWEBV ; 
    fprintf(fp, "SIMOPT_MWEBV:       %d   (modify MWEBV_SFD)\n",*iptr);


    if ( SNDATA.SIM_AV != NULLFLOAT ) {

      fptr = &SNDATA.SIM_AVTAU ; 
      fprintf(fp, "SIM_AVTAU:     %8.3f   (dN/dAV = exp(-AV/AVTAU)) \n", 
	      *fptr ) ;

      fptr = &SNDATA.SIM_AV ; 
      fprintf(fp, "SIM_AV:        %6.3f  mag  (host extinction at 5510 A)\n", 
	      *fptr ) ;
      fptr = &SNDATA.SIM_RV ; 
      fprintf(fp, "SIM_RV:        %6.3f   (CCM89 extinction parameter) \n", 
	      *fptr ) ;      
    }
    
    fptr = &SNDATA.SIM_MAGSMEAR_COH ; 
    fprintf(fp, "SIM_MAGSMEAR_COH:     %6.3f  \n", *fptr ) ;      


    // gal/SN flux-fraction
    fprintf(fp, "SIM_GALFRAC: "); NTMP = 0;
    for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
      ifilt_obs = SNDATA_FILTER.MAP[ifilt];
      fprintf(fp," %6.3f",SNDATA.SIM_GALFRAC[ifilt_obs] ) ; 
      NTMP++ ;
      if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
    }
    fprintf(fp,"  (%s F_gal/F_SNpeak for PSF=1'')\n", filtlist);
    


    fprintf(fp, "SIM_PEAKMAG: "); NTMP=0;
    for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
      ifilt_obs = SNDATA_FILTER.MAP[ifilt];
      fprintf(fp," %6.2f",SNDATA.SIM_PEAKMAG[ifilt_obs] ) ;
      NTMP++ ;
      if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; } 
    }
    fprintf(fp,"  (%s obs)\n", filtlist);


    if ( SNDATA.SIM_MODEL_INDEX == MODEL_LCLIB ) {
      fprintf(fp, "SIM_TEMPLATEMAG: "); NTMP=0;
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs = SNDATA_FILTER.MAP[ifilt];
	fprintf(fp," %6.2f",SNDATA.SIM_TEMPLATEMAG[ifilt_obs] ) ;
	NTMP++ ;
	if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; } 
      }
      fprintf(fp,"  (%s)\n", filtlist);      
    }

    fprintf(fp, "SIM_EXPOSURE: ");  NTMP=0;
    for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
      ifilt_obs = SNDATA_FILTER.MAP[ifilt];
      fprintf(fp," %6.1f",SNDATA.SIM_EXPOSURE_TIME[ifilt_obs] ) ; 
      NTMP++ ;
      if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
    }
    fprintf(fp,"  (%s obs)\n", filtlist);


    fptr = &SNDATA.SIM_PEAKMJD ; 
    fprintf(fp, "SIM_PEAKMJD:   %f  days \n", *fptr ) ;

    // write luminosity parameter for valid par only

    sprintf(ctmp,"STRETCH lumi-par"); // init comment about stretch

    if ( SNDATA.SIM_DELTA != NULLFLOAT ) {
      fptr = &SNDATA.SIM_DELTA ;
      fprintf(fp, "SIM_DELTA:     %6.3f  mag  (MLCS lumi-par) \n", *fptr ) ;
      sprintf(ctmp,"approx STRETCH for KCOR loopup");
    }
    if ( SNDATA.SIM_DM15 != NULLFLOAT ) {
      fptr = &SNDATA.SIM_DM15 ;
      fprintf(fp, "SIM_DM15:      %6.3f  mag  (DM15 lumi-par) \n", *fptr ) ;
      sprintf(ctmp,"approx STRETCH for KCOR lookup");
    }

    if ( SNDATA.SIM_SALT2alpha != NULLFLOAT ) {
      fptr = &SNDATA.SIM_SALT2alpha ;
      fprintf(fp, "SIM_SALT2alpha:  %7.3f   \n", *fptr ) ;

      fptr = &SNDATA.SIM_SALT2beta ;
      fprintf(fp, "SIM_SALT2beta:   %7.3f   \n", *fptr ) ;

      fptr = &SNDATA.SIM_SALT2gammaDM ; 
      fprintf(fp, "SIM_SALT2gammaDM: %6.3f  \n", *fptr ) ; 
    }

    if ( SNDATA.SIM_SALT2x0 != NULLFLOAT ) {
      fptr = &SNDATA.SIM_SALT2x0 ;
      fprintf(fp, "SIM_SALT2x0:   %8.4e   \n", *fptr ) ;
    }
    if ( SNDATA.SIM_SALT2x1 != NULLFLOAT ) {
      fptr = &SNDATA.SIM_SALT2x1 ;
      fprintf(fp, "SIM_SALT2x1:   %7.4f   \n", *fptr ) ;
    }
    if ( SNDATA.SIM_SALT2c != NULLFLOAT ) {
      fptr = &SNDATA.SIM_SALT2c ;
      fprintf(fp, "SIM_SALT2c:    %7.4f   \n", *fptr ) ;
    }

    if ( SNDATA.SIM_SALT2mB != NULLFLOAT ) {
      fptr = &SNDATA.SIM_SALT2mB ;
      fprintf(fp, "SIM_SALT2mB:   %7.4f   \n", *fptr ) ;
    }

    if ( SNDATA.SIM_STRETCH > 0 ) {
      fptr = &SNDATA.SIM_STRETCH ;
      fprintf(fp, "SIM_STRETCH:   %6.3f  \n", *fptr ) ;
    }

    // covmat-scatter info
    if ( SNDATA.SIMFLAG_COVMAT_SCATTER ) {
      for ( iscat=0; iscat < 3; iscat++ ) {
	fprintf(fp, "SIM_SCATTER[%s]:  %8.4f \n"
		, SNDATA.SIM_COVMAT_SCATTER_NAME[iscat]
		, SNDATA.SIM_COVMAT_SCATTER[iscat] ) ;
      }
    }


    iptr = &SNDATA.SIM_SEARCHEFF_MASK ;
    sprintf(ctmp,"bits 1,2=> found by software,spec" );
    fprintf(fp, "SIM_SEARCHEFF_MASK:  %d  (%s) \n", *iptr, ctmp );
    fptr = &SNDATA.SIM_SEARCHEFF_SPEC ;
    sprintf(ctmp,"spectro-search efficiency (ignores pipelines)");
    fprintf(fp, "SIM_SEARCHEFF_SPEC:  %6.4f  (%s) \n", *fptr, ctmp );

    if ( SNDATA.SIM_TRESTMIN != NULLFLOAT ) {
      fptr = &SNDATA.SIM_TRESTMIN ; 
      fprintf(fp, "SIM_TRESTMIN:  %7.2f   days \n", *fptr ) ;
      fptr = &SNDATA.SIM_TRESTMAX ; 
      fprintf(fp, "SIM_TRESTMAX:  %7.2f   days \n", *fptr ) ;
    }

    fptr = &SNDATA.SIM_RISETIME_SHIFT ;
    fprintf(fp, "SIM_RISETIME_SHIFT:   %3.1f days \n", *fptr ) ;
    fptr = &SNDATA.SIM_FALLTIME_SHIFT ;
    fprintf(fp, "SIM_FALLTIME_SHIFT:   %3.1f days \n", *fptr ) ;


    // SIMSED info
    if ( SNDATA.NPAR_SIMSED > 0 ) {
      fprintf(fp,"\n");
      fprintf(fp,"SIMSED_NPAR: %d \n", SNDATA.NPAR_SIMSED );
      for ( ipar = 0; ipar < SNDATA.NPAR_SIMSED; ipar++ ) {
	fprintf(fp,"%s:  %f\n"
		,SNDATA.SIMSED_KEYWORD[ipar]
		,SNDATA.SIMSED_PARVAL[ipar] );
      }
    }

    // BYOSED info (Dec 2018) 
    if ( SNDATA.NPAR_BYOSED > 0 ) {
      fprintf(fp,"\n");
      fprintf(fp,"BYOSED_NPAR: %d \n", SNDATA.NPAR_BYOSED );
      for ( ipar = 0; ipar < SNDATA.NPAR_BYOSED; ipar++ ) {
	fprintf(fp,"%s:  %f \n"
		,SNDATA.BYOSED_KEYWORD[ipar]
		,SNDATA.BYOSED_PARVAL[ipar] );
      }
    }


    // LCLIB info (Sep 8 2017)
    if ( SNDATA.NPAR_LCLIB > 0 ) {
      fprintf(fp,"\n");
      fprintf(fp,"LCLIB_NPAR: %d \n", SNDATA.NPAR_LCLIB );
      for ( ipar = 0; ipar < SNDATA.NPAR_LCLIB; ipar++ ) {
	fprintf(fp,"%s:  %f \n"
		,SNDATA.LCLIB_KEYWORD[ipar]
		,SNDATA.LCLIB_PARVAL[ipar] );
      }
    }

    
    // strong lens info (July 20 2019)
    if ( SNDATA.SIM_SL_FLAG ) {
      fprintf(fp,"\n");
      fprintf(fp,"SIM_STRONGLENS_ID:        %d   \n", SNDATA.SIM_SL_IDLENS );
      fprintf(fp,"SIM_STRONGLENS_z:         %.3f \n", SNDATA.SIM_SL_zLENS  );
      fprintf(fp,"SIM_STRONGLENS_DELAY:     %.3f  # days \n", 
	      SNDATA.SIM_SL_TDELAY );
      fprintf(fp,"SIM_STRONGLENS_XIMG:      %.3f  # arcsec\n",
	      SNDATA.SIM_SL_XIMG );
      fprintf(fp,"SIM_STRONGLENS_YIMG:      %.3f  # arcsec\n",
	      SNDATA.SIM_SL_YIMG );
      fprintf(fp,"SIM_STRONGLENS_MAGSHIFT:  %.3f \n", SNDATA.SIM_SL_MAGSHIFT );
      fprintf(fp,"SIM_STRONGLENS_NIMG:      %d   \n", SNDATA.SIM_SL_NIMG    );
      fprintf(fp,"SIM_STRONGLENS_IMGNUM:    %d   \n", SNDATA.SIM_SL_IMGNUM  );
    }
    
    fprintf(fp,"\n");

    // Jun 2017: SUBSAMPLE_INDEX
    if ( SNDATA.SUBSAMPLE_INDEX >= 0 ) {
      fprintf(fp,"SIM_SUBSAMPLE_INDEX: %d \n", SNDATA.SUBSAMPLE_INDEX);
    }

    
  } // end of FAKE > 0 if-block

  
  fflush(fp);

 START_EPOCHS:

  // check option to skip epoch info (i.e, for TERSE output)
  if ( NEWMJD <= 0 ) {
    fclose ( fp );
    return SUCCESS ;
  }

    
  // epoch info. Note that "epoch" is the sorted epoch index

  for ( inext = 1; inext <= SNDATA.NEWMJD; inext++ ) {

    imjd = SNDATA.UNSORTED_EPOCH[inext];  // unsorted index

    EPMIN = SNDATA.EPOCH_RANGE_NEWMJD[imjd][0] ;  
    EPMAX = SNDATA.EPOCH_RANGE_NEWMJD[imjd][1] ;  

    // make sure that the number of epochs here corresponds
    // to the number of filters per epoch

    if ( NFILT_TOT < EPMAX - EPMIN + 1 ) {

      printf("\n");
      // prepare list of filters for screen dump
      for ( iep   = EPMIN; iep <= EPMAX; iep++ ) {
	ifilt_obs = SNDATA.FILTINDX[iep];
	mjd       = SNDATA.MJD[iep];
	printf("  ifilt_obs(iep=%d) = %d (%c)   MJD=%9.3f\n", 
	       iep, ifilt_obs, FILTERSTRING[ifilt_obs], mjd );
      }

      if ( SNDATA.FAKE > 0 ) 
	printf("  (SIMLIB ID=%d) \n", SNDATA.SIM_LIBID );

      sprintf(c1err,"Epoch range %d to %d  (CID=%d NEWMJD=%d)",
	      EPMIN, EPMAX, cid, imjd );
      sprintf(c2err,"But only %d filters defined (%s)", 
	      NFILT_TOT, filtlist);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    
    }

    VERSION_INFO.NEPOCH_TOT++ ;


    fprintf(fp,"\n# ----------------------------------------------- \n");
    fprintf(fp, "  EPOCH: %d  \n", inext );  // write sorted epoch !!!

    fprintf(fp, "  MJD: %f  \n" , SNDATA.MJD[EPMIN] );

    fprintf(fp, "  TELESCOPE: %s    ", SNDATA.TELESCOPE[EPMIN] );

    // check for optional things

    JTMP = SNDATA.SEARCH_RUN[EPMIN] ;
    if ( JTMP != NULLINT ) fprintf(fp, "SEARCH_RUN: %d  ", JTMP);

    JTMP = SNDATA.TEMPLATE_RUN[EPMIN] ;
    if ( JTMP != NULLINT ) fprintf(fp, "  TEMPLATE_RUN: %d  ", JTMP ) ;

    JTMP = SNDATA.QMASK[EPMIN];
    if ( JTMP  != NULLINT ) fprintf(fp, "QMASK: %d   ", JTMP );

    fprintf(fp,"\n  ");

    /* xxxxxxxx mark delete xxxxxxxxxxxxx
    JTMP = SNDATA.IDCCD[EPMIN];
    if ( JTMP  != NULLINT ) fprintf(fp, "IDCCD: %d   ", JTMP );
    xxxxxxxxx */

    cptr = SNDATA.FIELDNAME[EPMIN] ;
    if ( strcmp(cptr,"NULL") != 0 )  fprintf(fp, "FIELD: %s   ", cptr );

    fprintf(fp,"\n");

    if ( SNDATA.IDATE[EPMIN] > 20040000 ) 
      fprintf(fp,"  PROCESS_DATE: %s \n", SNDATA.DATE[EPMIN] );


    /* xxxxxxxx mark delete Jun 23 2019 xxxxxxx
    FTMP = SNDATA.CLOUDCAM_AVG[EPMIN];
    if ( FTMP != NULLFLOAT ) {
      fprintf(fp, "  CLOUDCAM_AVG: %7.2f   ", SNDATA.CLOUDCAM_AVG[EPMIN] );
      fprintf(fp, "  CLOUDCAM_SIG: %5.2f \n", SNDATA.CLOUDCAM_SIG[EPMIN] );
    }
    FTMP = SNDATA.MOONDIST[EPMIN] ;
    if ( FTMP != NULLFLOAT ) {
      fprintf(fp, "  MOONDIST:    %7.2f deg   ", SNDATA.MOONDIST[EPMIN] );
      fprintf(fp, "  MOONPHASE:  %5.2f \n",    SNDATA.MOONPHASE[EPMIN] );
    }
    FTMP = SNDATA.AIRMASS[EPMIN] ;
    if ( FTMP != NULLFLOAT ) fprintf(fp, "  AIRMASS:  %6.3f \n", FTMP );
    xxxxxxxxx */

    // now write info vs. fitler-band

    NFILT_EP = 0;
    fprintf( fp, "\n  %16s ", "PASSBAND:" );
    for ( iep=EPMIN; iep <= EPMAX; iep++ ) {
      ifilt_obs = SNDATA.FILTINDX[iep]; 
      fprintf(fp,"     %c   ", FILTERSTRING[ifilt_obs] );
      NFILT_EP++;
    }
    fprintf( fp, "\n" );


    iptr = &SNDATA.SEARCH_FIELD[EPMIN] ;
    if ( *iptr != NULLINT ) 
      istat = wr_filtband_int(fp, "SEARCH_FIELD:", NFILT_EP, iptr, blank,1 );

    iptr = &SNDATA.TEMPLATE_FIELD[EPMIN] ;
    if ( *iptr != NULLINT ) 
      istat = wr_filtband_int(fp, "TEMPLATE_FIELD:", NFILT_EP, iptr, blank,1);

    iptr = &SNDATA.PHOTFLAG[EPMIN] ;
    if ( *iptr != NULLINT  && fake == 0 ) 
      istat = wr_filtband_int ( fp, "PHOTFLAG:", NFILT_EP, iptr, blank,1) ;

    // xxx mark delete     istat = wr_filtband_int ( fp, "PHOTOMETRYFLAG:", NFILT_EP, iptr, blank,1) ;


    if ( SNDATA.WRFLAG_BLINDTEST ) goto FLUXCAL ;

    wstat = WRSTAT ( IFLAG_WR, SNDATA.GAIN[EPMIN]  );
    if ( wstat == 1 ) {
      fptr = &SNDATA.GAIN[EPMIN] ;
      istat = wr_filtband_float ( fp, "GAIN:", NFILT_EP, fptr, "e/ADU", 3 ) ;

      fptr = &SNDATA.READNOISE[EPMIN] ;
      istat = wr_filtband_float ( fp, "RDNOISE:", NFILT_EP, fptr, "e-", 3 ) ;
    }


    FTMP = SNDATA.XPIX[EPMIN] ;    
    if ( SNDATA.FAKE == 0  &&  FTMP != NULLFLOAT ) {

      fptr = &SNDATA.XPIX[EPMIN] ;
      istat = wr_filtband_float ( fp, "XPIXEL:", NFILT_EP, fptr, "(pixels)", 3 ) ;
      fptr = &SNDATA.YPIX[EPMIN] ;
      istat = wr_filtband_float ( fp, "YPIXEL:", NFILT_EP, fptr, "(pixels)", 3 ) ;
      fptr = &SNDATA.EDGEDIST[EPMIN] ;
      istat = wr_filtband_float ( fp, "EDGEDIST:", NFILT_EP, fptr, "(pixels)", 3 ) ;
    }


    FTMP = SNDATA.SKY_SIG[EPMIN];
    if ( FTMP != NULLFLOAT ) {
      fptr = &SNDATA.SKY_SIG[EPMIN] ;
      istat = wr_filtband_float ( fp, "SKY_SIG:", 
				  NFILT_EP, fptr, "ADU/pix", 2 ) ;
    }


    fptr = &SNDATA.PSF_SIG1[EPMIN] ;
    if ( *fptr != NULLFLOAT )
      istat = wr_filtband_float ( fp, "PSF_SIG1:", NFILT_EP, fptr, "pixels", 3 ) ;

    FTMP = SNDATA.PSF_SIG2[EPMIN] ;
    if ( FTMP > 0.0 ) {

      fptr = &SNDATA.PSF_SIG2[EPMIN] ;
      istat = wr_filtband_float ( fp, "PSF_SIG2:", NFILT_EP, fptr, "pixels", 3 ) ;

      fptr = &SNDATA.PSF_RATIO[EPMIN] ;      
      istat = wr_filtband_float ( fp, "PSF_RATIO:", NFILT_EP, fptr, "(at origin)", 3 ) ;
    }

    /* xxxxxxxxxxxx mark delete Jun 24 2019 xxxxxxxx
    // write FLUX in ADU (or uJy)
    fprintf(fp," \n" );
    fptr = &SNDATA.FLUX[EPMIN] ;
    istat = wr_filtband_float ( fp, "FLUX:", NFILT_EP, fptr, FLUXUNIT, 2 ) ;
    fptr = &SNDATA.FLUX_ERRTOT[EPMIN] ;
    istat = wr_filtband_float (fp,"FLUX_ERRTOT:",NFILT_EP,fptr,FLUXUNIT,3);
    xxxxxxxxxxxxxx */

    // write calibrate fluxes: 

  FLUXCAL:

    fprintf(fp," \n" );
    fptr = &SNDATA.FLUXCAL[EPMIN] ;
    istat = wr_filtband_float ( fp, "FLUXCAL:", NFILT_EP, fptr, 
				"10^(11-.4*m)",2) ;
    fptr = &SNDATA.FLUXCAL_ERRTOT[EPMIN] ;
    istat = wr_filtband_float ( fp, "FLUXCAL_ERRTOT:", NFILT_EP, fptr, 
				blank, 3 ) ;

    // write magnitudes

    if ( SNDATA.WRFLAG_BLINDTEST ) goto END_OF_EPOCH ;

    fprintf(fp," \n" );
    fptr = &SNDATA.MAG[EPMIN] ;
    istat = wr_filtband_float ( fp, "MAG:", NFILT_EP, fptr, " ", 4 ) ;
    fptr = &SNDATA.MAG_ERRPLUS[EPMIN] ;
    istat = wr_filtband_float ( fp, "MAG_ERRPLUS:", NFILT_EP, fptr, blank, 4 ) ;
    fptr = &SNDATA.MAG_ERRMINUS[EPMIN] ;
    istat = wr_filtband_float ( fp, "MAG_ERRMINUS:", NFILT_EP, fptr, blank, 4 ) ;

    ctmp[0]=0;
    // write zero points
    fptr = &SNDATA.ZEROPT[EPMIN] ;
    istat = wr_filtband_float ( fp, "ZEROPT:", NFILT_EP, fptr, ctmp, 4 ) ;
    fptr = &SNDATA.ZEROPT_ERR[EPMIN] ;
    istat = wr_filtband_float ( fp, "ZEROPT_ERR:", NFILT_EP, fptr, ctmp, 4 ) ;

    // write subtraction error for sky & galaxy (SDSS data only)

    FTMP = SNDATA.SKYSUB_ERR[EPMIN] ;
    if ( FTMP != NULLFLOAT ) {
      fptr = &SNDATA.SKYSUB_ERR[EPMIN] ;
      istat = wr_filtband_float ( fp, "SKYSUB_ERR:", NFILT_EP, fptr,  
				  "(fluxcal)", 3 ) ;
    }

    FTMP = SNDATA.SKYSUB_ERR[EPMIN] ;
    if ( FTMP != NULLFLOAT ) {
      fptr = &SNDATA.GALSUB_ERR[EPMIN] ;
      istat = wr_filtband_float ( fp, "GALSUB_ERR:", NFILT_EP, fptr,  
				  "(fluxcal)", 3 ) ;
    }


  // write SIMEPOCH_XXX variables if FAKE > 0
    if ( LWRITE_SIMFLAG > 0  ) {
      wr_SIMKCOR(fp,EPMIN,EPMAX);
    }

  END_OF_EPOCH:

    fprintf(fp, "\n  END_OF_EPOCH: %d \n", inext );

  } // end of  epoch" loop


  fprintf(fp," \n  END_OF_SN: %d \n", cid );
  
  fclose ( fp );
 
  return SUCCESS;

} // end of wr_SNDATA




// ***************************
void wr_HOSTGAL(FILE *fp) {


  // Dec 17 2012 - write HOSTGAL info to current ASCII data file (*fp)
  // May 16,2013 - write no more than 10 per line to avoid lines that
  //               are too long.
  // Dec 18 2015 - write specz
  // Nov 13 2019 - fix to work with NGAL>1

  int ifilt, ifilt_obs, NTMP, igal, NGAL ;
  char PREFIX[20] = "HOSTGAL";
  char filtlist[MXFILTINDX], ctmp[100] ;
  
  // --------------- BEGIN --------------

  
  sprintf(filtlist,"%s", SNDATA_FILTER.LIST );

  NGAL = SNDATA.HOSTGAL_NMATCH[0];
  if ( NGAL > MXHOSTGAL ) { NGAL = MXHOSTGAL ; }

  fprintf(fp, "%s_NMATCH:    %d  \n",  
	  PREFIX, SNDATA.HOSTGAL_NMATCH[0] );
  fprintf(fp, "%s_NMATCH2:   %d  \n",  
	  PREFIX, SNDATA.HOSTGAL_NMATCH[1] );

  for(igal=0; igal < NGAL; igal++ ) {

    if ( igal > 0 ) { sprintf(PREFIX,"HOSTGAL%d", igal+1); }

    fprintf(fp, "%s_OBJID:    %lld  \n",  
	    PREFIX, SNDATA.HOSTGAL_OBJID[igal] );

    fprintf(fp, "%s_PHOTOZ:   %6.4f  +- %6.4f \n", PREFIX,
	    SNDATA.HOSTGAL_PHOTOZ[igal], 
	    SNDATA.HOSTGAL_PHOTOZ_ERR[igal]);

    fprintf(fp, "%s_SPECZ:    %6.4f  +- %6.4f \n", PREFIX,
	  SNDATA.HOSTGAL_SPECZ[igal], SNDATA.HOSTGAL_SPECZ_ERR[igal] ); 
  
    fprintf(fp, "%s_RA:       %.6f    # deg \n", 
	    PREFIX, SNDATA.HOSTGAL_RA[igal] );
    fprintf(fp, "%s_DEC:      %.6f    # deg \n", 
	    PREFIX, SNDATA.HOSTGAL_DEC[igal] );

    fprintf(fp, "%s_SNSEP:    %6.3f    # arcsec \n", 
	    PREFIX, SNDATA.HOSTGAL_SNSEP[igal] );
    fprintf(fp, "%s_DDLR:     %6.3f    # SNSEP/DLR  \n", 
	    PREFIX, SNDATA.HOSTGAL_DDLR[igal] );
    
    if ( igal==0 ) {
      fprintf(fp, "HOSTGAL_CONFUSION:  %6.3f  \n", 
	      SNDATA.HOSTGAL_CONFUSION );
    }

    if ( SNDATA.HOSTGAL_LOGMASS[igal] > 0.0 ) {
      fprintf(fp, "%s_LOGMASS:  %6.3f +- %6.3f   # log10(Mgal/Msolar)\n", 
	      PREFIX, 
	      SNDATA.HOSTGAL_LOGMASS[igal], 
	      SNDATA.HOSTGAL_LOGMASS_ERR[igal] );

      fprintf(fp, "%s_sSFR:  %6.3e +- %6.3e  \n",
	      PREFIX, 
	      SNDATA.HOSTGAL_sSFR[igal], 
	      SNDATA.HOSTGAL_sSFR_ERR[igal] );
    }

    // if MAGOBS has been read for any filter, then write MAGOBS
    // for all filters.

    if ( SNDATA.HOSTLIB_NFILT_MAGOBS > 0 ) {
      fprintf(fp, "%s_MAG:    ", PREFIX ); NTMP=0;    
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs = SNDATA_FILTER.MAP[ifilt] ;
	fprintf(fp,"%6.2f ", SNDATA.HOSTGAL_MAG[igal][ifilt] );
	NTMP++ ;
	if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
      }
      fprintf(fp,"# %s\n", filtlist) ;
    }

    if ( SNDATA.HOSTLIB_NFILT_MAGOBS > 0 ) {
      fprintf(fp, "%s_MAGERR: ", PREFIX ); NTMP=0;    
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs = SNDATA_FILTER.MAP[ifilt] ;
	fprintf(fp,"%6.2f ", SNDATA.HOSTGAL_MAGERR[igal][ifilt] );
	NTMP++ ;
	if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
      }
      fprintf(fp,"# %s\n", filtlist) ;
    }
    
    fprintf(fp,"\n");

  } // end igal loop

  // ---------- surface brightness -----------

  sprintf(ctmp,"%s/asec^2",filtlist );
  if ( (SNDATA.HOSTGAL_USEMASK & 4) > 0 ) {
    fprintf(fp,"HOSTGAL_SB_FLUXCAL:    " ); NTMP=0 ;
    for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
      ifilt_obs = SNDATA_FILTER.MAP[ifilt];
      fprintf(fp," %6.2f",SNDATA.HOSTGAL_SB_FLUX[ifilt] ) ;
      NTMP++ ;
      if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
    }
    fprintf(fp,"  # %s\n", ctmp );
  }


  if ( (SNDATA.HOSTGAL_USEMASK & 8) > 0 ) {
    fprintf(fp,"HOSTGAL_SB_FLUXCAL_ERR:    " ); NTMP=0 ;
    for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
      ifilt_obs = SNDATA_FILTER.MAP[ifilt];
      fprintf(fp," %6.2f",SNDATA.HOSTGAL_SB_FLUXERR[ifilt_obs] ) ;
      NTMP++ ;
      if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
    }
    fprintf(fp," # %s \n", ctmp );
  }


} // end of wr_HOSTGAL


// ********************************
void wr_SIMKCOR(FILE *fp, int EPMIN, int EPMAX) {

  // Mar 2019: remove tabs
  float *fptr ;
  int i;
  char tmpstring[80];
  //  char fnam[] = "wr_SIMKCOR" ;

  // -------------------- BEGIN ------------

    fprintf(fp, "\n" );

    fptr = &SNDATA.SIMEPOCH_TREST[EPMIN] ;
    fprintf(fp,"   SIM_TREST:   %7.3f  rest-frame days \n", *fptr ) ;
    fptr = &SNDATA.SIMEPOCH_TOBS[EPMAX] ;
    fprintf(fp,"   SIM_TOBS:    %7.3f  obs-frame days \n",  *fptr ) ;

    sprintf(tmpstring,"SIM_MAG:            ") ;
    for ( i=EPMIN; i <= EPMAX; i++ ) {
      sprintf(tmpstring,"%s %7.3f", 
	      tmpstring, SNDATA.SIMEPOCH_MAG[i]);
    }
    fprintf(fp,"   %s \n", tmpstring);

    // Jun 21, 2009: write model-mag error
    sprintf(tmpstring,"SIM_MODELMAGERR:    ") ;
    for ( i=EPMIN; i <= EPMAX; i++ ) {
      sprintf(tmpstring,"%s %7.3f", 
	      tmpstring, SNDATA.SIMEPOCH_MODELMAGERR[i]);
    }
    fprintf(fp,"   %s \n", tmpstring);

    // Feb 2, 2009: write intrinsic mag-smearing 
    sprintf(tmpstring,"SIM_MAGSMEAR:      ") ;
    for ( i=EPMIN; i <= EPMAX; i++ ) {
      sprintf(tmpstring,"%s %7.3f", tmpstring, 
	      SNDATA.SIMEPOCH_MAGSMEAR[i]);
    }
    fprintf(fp,"   %s \n", tmpstring);


    // skip K-cor stuff for observer-frame model
    if ( VERSION_INFO.GENFRAME_SIM == 2 ) return ;

    sprintf(tmpstring,"SIM_WARPCOL_SYMBOL: ") ;
    for ( i=EPMIN; i <= EPMAX; i++ ) {
      sprintf(tmpstring,"%s %7s", 
	      tmpstring, SNDATA.SIMEPOCH_WARPCOLNAM[i]);
    }
    fprintf(fp,"   %s \n", tmpstring);


    sprintf(tmpstring,"SIM_WARPCOL_VALUE:  ") ;
    for ( i=EPMIN; i <= EPMAX; i++ ) {
      sprintf(tmpstring,"%s %7.3f", 
	      tmpstring, SNDATA.SIMEPOCH_WARPCOLVAL[i]);
    }
    fprintf(fp,"   %s \n", tmpstring);

    sprintf(tmpstring,"SIM_AVWARP:         ") ;
    for ( i=EPMIN; i <= EPMAX; i++ ) {
      sprintf(tmpstring,"%s %7.3f", 
	      tmpstring, SNDATA.SIMEPOCH_AVWARP[i]);
    }
    fprintf(fp,"   %s \n", tmpstring);


    sprintf(tmpstring,"SIM_KCOR_SYMBOL:    ") ;
    for ( i=EPMIN; i <= EPMAX; i++ ) {
      sprintf(tmpstring,"%s %7s", 
	      tmpstring, SNDATA.SIMEPOCH_KCORNAM[i]);
    }
    fprintf(fp,"   %s \n", tmpstring);


    sprintf(tmpstring,"SIM_KCOR_VALUE:     ") ;
    for ( i=EPMIN; i <= EPMAX; i++ ) {
      sprintf(tmpstring,"%s %7.3f", tmpstring, 
	      SNDATA.SIMEPOCH_KCORVAL[i]);
    }
    fprintf(fp,"   %s \n", tmpstring);


}  // end of wr_SNDATA


// ****************************************************
int WRSTAT ( int wrflag, float value ) {

  // return 1 to write; 0 to suppress write

  int istat, OVP ;

  // ---------- BEGIN ----------

  istat = 0 ;

  OVP = (wrflag & WRITE_MASK_LCMERGE) ;
  if ( OVP > 0 ) { istat = 1; }

  OVP = (wrflag & WRITE_MASK_SIM_SNANA) ;
  if ( OVP > 0 &&  value != NULLFLOAT ) { istat = 1; }

  /* xxxxxxxx mark delete Jan 23 2018 xxxxxxxxxxxxxxx
  xxx if ( wrflag == WRITE_MASK_PHOTOMETRY &&  value != NULLFLOAT ) 
    { istat = 1; }
  xxxxxxxxxxxxxxxxxxxxxxxxxx */

  return istat ;

}


// **********************************
int header_merge(FILE *fp, char *auxheader_file) {

  // Jun 19, 2009
  // Merge contents of auxheader_file into existing file with
  // pointer fp

  FILE *fp_aux;
  char cline[MXPATHLEN];

  // -------- BEGIN -----------

  if ( (fp_aux = fopen(auxheader_file, "rt"))==NULL ) return SUCCESS ;

  while( (fgets(cline, 100, fp_aux)) != NULL) 
    { fprintf(fp,"%s", cline ); }
  
  fprintf(fp,"\n");

  fclose(fp_aux);
  return SUCCESS ;

} // end of header_merge

// ******************************************************
int  fluxcal_SNDATA ( int iepoch, char *magfun ) {


  /*********
    fill SNDATA.FLUXCAL[iepoch][ifilt] 
               and
         SNDATA.FLUXCAL_ERRTOT[iepoch][ifilt] 


   = FLUX * 10**{-0.4*ZP}    if  magfun = "log10"

   = SINH (  )               if  magfun = "asinh"


  Mar 2 2013: float -> double and FLUXCAL_SCALE -> ZEROPOINT_FLUXCAL

  Mar 14, 2014:  remove requirement of asinh mag < 28; 
                 keep extreme negative fluxes

  Sep 5 2016: magfun = asinh is now obsolete

  Jan 3 2018: check for saturation

  *********/


  double
     mag, mag_err, mag_tmp
    ,ZP, ZP_err, ZP_scale, ZP_sig
    ,flux, flux_err
    ,fluxcal, fluxcal_err
    ,ferrp, ferrm
    ,tmperr, arg
    ,sqerrtmp, relerr
    ;

  int VALID_MAGFUN, IFILT, LTMP;

  char fnam[] = "fluxcal_SNDATA" ;

  // ------------- BEGIN ----------------


  VALID_MAGFUN = 0;

  ZP        = SNDATA.ZEROPT[iepoch];
  ZP_err    = SNDATA.ZEROPT_ERR[iepoch];
  ZP_sig    = SNDATA.ZEROPT_SIG[iepoch];
  flux      = SNDATA.FLUX[iepoch];
  flux_err  = SNDATA.FLUX_ERRTOT[iepoch];
  IFILT     = SNDATA.FILTINDX[iepoch]; // absolute obs filter indx

  fluxcal     = NULLDOUBLE ;
  fluxcal_err = NULLDOUBLE ;

  mag = 99.0 ; mag_err=99.0 ;

  // make some idiot checks

  if ( flux_err <   0.  ) { return ERROR; }

  // Jan 2018 : for saturated epoch, set fluxcal = flux
  if ( SNDATA.NPE_ABOVE_SAT[iepoch] > 0 ) {
    SNDATA.FLUXCAL[iepoch]        = flux ;
    SNDATA.FLUXCAL_ERRTOT[iepoch] = flux_err ;
    return(SUCCESS) ;
  }

  // for asinh mags, get calibrated flux from mag

  if ( strcmp(magfun,"asinh") == 0 ) {

    VALID_MAGFUN = 1 ;

    fluxcal = asinhinv ( mag, IFILT ) ;

    mag_tmp = mag - mag_err ;
    ferrp  = asinhinv ( mag_tmp, IFILT) - fluxcal;

    mag_tmp    = mag+mag_err ;
    ferrm      = fluxcal - asinhinv ( mag_tmp, IFILT ) ;

    fluxcal_err = (ferrp+ferrm)/2.0;

    if ( fluxcal_err > 10000. ) { fluxcal_err = 9999. ; }

    // xxxxxxxxxxxxxx dump + and - errors for test
    double mjd ;
    mjd = SNDATA.MJD[iepoch] ;
    LTMP = 0;  if ( fabs(mjd - 53626.37 ) < 0.1 ) { LTMP = 1; }
    if ( LTMP == -9 ) {
      printf(" MJD=%9.3f  magerr=%5.2f  "
	     "fluxerr(%c) = avg(+%5.1f -%5.1f) = %5.1f\n",
	     mjd, mag_err, FILTERSTRING[IFILT],  ferrp, ferrm, fluxcal_err );
    }
    // xxxxxxxxxxxx

  }

  // for log10 mag, get calibrated flux from FLUXCAL

  if ( strcmp(magfun,"log10") == 0 ) {
    VALID_MAGFUN = 1 ;
    arg      = -0.4 * (ZP - ZEROPOINT_FLUXCAL_DEFAULT) ;
    ZP_scale = pow(TEN,arg) ;
 
    if ( flux_err >= 0.0 && ZP > 10.0  && ZP_sig >= 0.0 ) {

      fluxcal     = flux     * ZP_scale ;
      fluxcal_err = flux_err * ZP_scale ;

      // add ZP error here for FLUXCAL; 
      // note that flux in ADU does not have this ZP error.
      relerr      = powf(TEN,0.4*ZP_sig) - 1.0 ;
      tmperr      = fluxcal * relerr ;
      sqerrtmp    = fluxcal_err*fluxcal_err + tmperr*tmperr ;
      fluxcal_err = sqrt(sqerrtmp) ;
    }
  }

  if ( VALID_MAGFUN == 0 ) {
    sprintf(c1err,"Invalid mag function: %s ", magfun );
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
  }


  // load structure

  SNDATA.FLUXCAL[iepoch]        = fluxcal ;
  SNDATA.FLUXCAL_ERRTOT[iepoch] = fluxcal_err ;

  if ( IFILT == -666 ) 
    printf("Ep=%d flt=%d, ZP = %f +- %f, fluxcal=%f +- %f \n",
	   iepoch, IFILT, ZP, ZP_sig, fluxcal, fluxcal_err );


  return(SUCCESS) ;

}  // end of  fluxcal_SNDATA



// ******************************************* 
double asinhinv(double mag, int ifilt) {

  // Invert SDSS mag to get calibrated flux.

  // define SDSS softening parameter vs. filter
  double bb[6]     = {0.0, 1.4e-10, 0.9e-10, 1.2e-10, 1.8e-10, 7.4e-10};
  double arg, b, fluxCal, magoff, fluxScale;

    // ------------ BEGIN ------------

  magoff  = -2.5/log(10.0);
  b       = bb[ifilt];
  arg     = mag/magoff - log(b);

  fluxScale = pow(TEN,0.4*ZEROPOINT_FLUXCAL_DEFAULT);
  fluxCal   = fluxScale * (2*b) * sinh(arg);

  // xxx delete Mar 2013:  fluxcal = FLUXCAL_SCALE * 2*b * sinh(arg);
 
  return fluxCal;

} // end of asinhinv



// **********************************************
int sort_epochs_bymjd ( void ) {

  /*******
   Created  May 18, 2006
   Fill SNDATA[isn].UNSORTED_EPOCH[epoch]
   used by wr_SNDATA to write epochs in order of MJD
   This sorting preserves time-order when combining
   data from different telescopes.

   Aug 20, 2007: modify for new index notation where
                 epoch runs over epochs and filters.
  
   Jun 19, 2009: remove "isn" arg

  *****/

  int NEPOCH;
  int epoch, epoch_tmp, rank;
  int EPMIN, EPMIN_tmp;
  int ISRANKED[MXEPOCH];
 
  float mjd, mjd_tmp;

  // --------------- BEGIN ------------------
 
  NEPOCH = SNDATA.NEWMJD;

  // init
  for ( epoch=1; epoch <= NEPOCH; epoch++ ) {
    SNDATA.UNSORTED_EPOCH[epoch] = NULLINT;
    ISRANKED[epoch] = 0;
  }

  // fill  array

  for ( epoch = 1; epoch <= NEPOCH; epoch++ ) {

    EPMIN = SNDATA.EPOCH_RANGE_NEWMJD[epoch][0] ;  

    mjd = SNDATA.MJD[EPMIN] ;

    // now find sorted "rank" of this unsorted "epoch".

    rank = 1;

    for ( epoch_tmp=1; epoch_tmp<=NEPOCH; epoch_tmp++ ) {
      EPMIN_tmp = SNDATA.EPOCH_RANGE_NEWMJD[epoch_tmp][0] ;  
      mjd_tmp  = SNDATA.MJD[EPMIN_tmp] ;
      if ( mjd >  mjd_tmp ) rank++ ;
      if ( mjd == mjd_tmp && epoch < epoch_tmp ) rank++ ;
    }

    SNDATA.UNSORTED_EPOCH[rank] = epoch;

  }  // end of epoch loop


  // check that everything was filled.
  for ( epoch=1; epoch <= NEPOCH; epoch++ ) { 

    sprintf(c1err,"UNSORTED_EPOCH[epoch %d] = %d",
	      epoch, SNDATA.UNSORTED_EPOCH[epoch] );

    if ( SNDATA.CID == -5 )  printf("%s \n", c1err);
   
    if ( SNDATA.UNSORTED_EPOCH[epoch] == NULLINT )
      errmsg(SEV_FATAL, 0, "sort_epochs", c1err, BLANK_STRING );

  }


  return SUCCESS;

} // end of sort_epochs_bymjd



/*  xxxxxxxxxxxxx mark delete Jan 23 2018 xxxxxxxxxxxxxxx
int IDTELESCOPE ( char *telescope ) {
  // returns integer ID of *telescope 
  int ID;
  //  char fnam[] = " IDTELESCOPE" ;
  // ------------- BEGIN ------------
  ID = IDTEL_SDSS ; 
  if ( strcmp(telescope,"sdss") == 0 ) ID = IDTEL_SDSS ;
  if ( strcmp(telescope,"SDSS") == 0 ) ID = IDTEL_SDSS ;
  if ( strcmp(telescope,"mdm24m") == 0 ) ID = IDTEL_MDM ;
  if ( strcmp(telescope,"uh88") == 0 ) ID = IDTEL_UH88 ;
  if ( strcmp(telescope,"UH88") == 0 ) ID = IDTEL_UH88 ;
  //  if ( ID == NULLINT ) {
    //    sprintf(c1err,"Telescope '%s' is not recognized", telescope );
    //    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
  //  }
  return ID ;
} // end of IDTELESCOPE
xxxxxxxxxxxx end delete xxxxxxxxxxxxxx */


// ********************************************
int rd_SNDATA ( void ) {

  /******
   Read SN data text file into SNDATA[isn] structure
   for this "isn".

   Store epoch only if search run is in SN.LIST;
   i.e., purge bad runs.


  Dec 29 2017: use open_TEXTgz to read gzipped files.

  ********/

  char fnam[] = "rd_SNDATA" ;

  char 
    inFile[100], c_get[80], line_passband[100]
    ,c_filt[2], varname[40], filtlist[MXFILTINDX], *ptrtok
    ;

  int 
    cid, ifilt_tmp[20], epoch, EPMIN, EPMAX, eptmp, NEWMJD
    ,NFILT_DEF, MINFILT_DEF, NFILT_NEWMJD, NTMP
    ;
  
  float fluxmax    = 1.0E7 ; 
  float fluxerrmax = 1.0E6 ;

  int   *iptr;  // pointer to integer fitler-band data
  float *fptr;  //            float

  FILE *fp;


  //------------ BEGIN ------------

  // store filename in local variables

  sprintf ( inFile, "%s", SNDATA.SNFILE_INPUT );

  fp = fopen(inFile, "rt"); 
  if ( fp == NULL ) {
      sprintf(c1err,"Cannot open: %s", inFile );
      errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
      fclose(fp);  return ERROR;
    }

  printf(" READ  %s", inFile );
  fflush(stdout);

  // ---------------------------------------
  // start parsing
  
  epoch   = 0;  // epoch * filters
  NEWMJD  = 0;  // just new MJDs
  EPMIN = EPMAX = NFILT_DEF = MINFILT_DEF = NFILT_NEWMJD = 0 ;

  while( fscanf(fp,"%s", c_get ) != EOF ) {


    if ( strcmp(c_get,"SDSS-SN:")==0 ) {
      readint ( fp, 1, &cid );  // read CID
      SNDATA.CID = cid ;          
      epoch   = 0;
      NEWMJD  = 0;
      printf(" (cid=%6d) " , cid );
    }

    if ( strcmp(c_get,"RA:")==0 && NEWMJD == 0 )
      readdouble ( fp, 1, &SNDATA.RA ) ;

    if ( strcmp(c_get,"DEC:")==0 && NEWMJD == 0)
      readdouble ( fp, 1, &SNDATA.DEC ) ;

    if ( strcmp(c_get,"FAKE:")==0 && NEWMJD == 0)
      readint ( fp, 1, &SNDATA.FAKE ) ;

    if ( strcmp(c_get,"MAGTYPE:")==0 && NEWMJD == 0)
      readchar ( fp, SNDATA.MAGTYPE ) ;
    if ( strcmp(c_get,"MAGREF:")==0 && NEWMJD == 0)
      readchar ( fp, SNDATA.MAGREF ) ;

    if ( strcmp(c_get,"FILTERS:")==0 && NEWMJD == 0) {
      readchar ( fp, filtlist ) ;
      NFILT_DEF = PARSE_FILTLIST(filtlist, SNDATA_FILTER.MAP );
      MINFILT_DEF = SNDATA_FILTER.MAP[0] ;
      SNDATA_FILTER.NDEF = NFILT_DEF ;
    }

    if ( strcmp(c_get,"SEARCH_TYPE:")==0 && NEWMJD == 0)
      readint ( fp, 1, &SNDATA.SEARCH_TYPE ) ;
    if ( strcmp(c_get,"SEARCH_PEAKMJD:")==0 && NEWMJD == 0)
      readfloat ( fp, 1, &SNDATA.SEARCH_PEAKMJD ) ;

    if ( strcmp(c_get,"REDSHIFT_FINAL:")==0 && NEWMJD == 0) {
      read_redshift ( fp, &SNDATA.REDSHIFT_FINAL, 
		      &SNDATA.REDSHIFT_FINAL_ERR ) ;

      /* xxxxx mark delete Jan 5 2018 xxxxxxxx
      fluxmax = snfluxmax(SNDATA.REDSHIFT_FINAL);	 
      if ( fluxmax < 5000. ) fluxmax = 5000. ;

      fluxerrmax = 0.1*fluxmax ;
      if ( fluxerrmax < 3000. ) fluxerrmax = 3000. ;
      xxxxxx */
    }
        
    if ( strcmp(c_get,"MWEBV:") == 0 && NEWMJD == 0 )
      { readfloat ( fp, 1, &SNDATA.MWEBV ) ; }

    if ( strcmp(c_get,"NEPOCH_PRESN:") == 0  &&  NEWMJD == 0 ) {
      iptr = &SNDATA.NPRESN[MINFILT_DEF] ;
      readint ( fp, NFILT_DEF, iptr );
    }

    // ----
    if ( strcmp(c_get,"HOSTGAL_SB_FLUX:") == 0  &&  NEWMJD == 0 ) {
      fptr = &SNDATA.HOSTGAL_SB_FLUX[MINFILT_DEF] ;
      readfloat ( fp, NFILT_DEF, fptr );
      SNDATA.HOSTGAL_USEMASK |= 4 ;
    }
    if ( strcmp(c_get,"HOSTGAL_SB_FLUXCAL:") == 0  &&  NEWMJD == 0 ) {
      fptr = &SNDATA.HOSTGAL_SB_FLUX[MINFILT_DEF] ;
      readfloat ( fp, NFILT_DEF, fptr );
      SNDATA.HOSTGAL_USEMASK |= 4 ;
    }
    if ( strcmp(c_get,"HOSTGAL_SB_FLUXERR:") == 0  &&  NEWMJD == 0 ) {
      fptr = &SNDATA.HOSTGAL_SB_FLUXERR[MINFILT_DEF] ;
      readfloat ( fp, NFILT_DEF, fptr );
      SNDATA.HOSTGAL_USEMASK |= 8 ;
    }
    if ( strcmp(c_get,"HOSTGAL_SB_FLUXCAL_ERR:") == 0  &&  NEWMJD == 0 ) {
      fptr = &SNDATA.HOSTGAL_SB_FLUXERR[MINFILT_DEF] ;
      readfloat ( fp, NFILT_DEF, fptr );
      SNDATA.HOSTGAL_USEMASK |= 8 ;
    }


    // check for sim stuff (only if fake > 0 )
    if ( SNDATA.FAKE > 0 ) {

      if ( strcmp(c_get,"SIM_COMMENT:")==0 ) 
	{ readchar ( fp, SNDATA.SIM_COMMENT ) ; }

      if ( strcmp(c_get,"SIM_REDSHIFT_HELIO:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_REDSHIFT_HELIO ) ; }

      if ( strcmp(c_get,"SIM_REDSHIFT_CMB:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_REDSHIFT_CMB ) ; }

      if ( strcmp(c_get,"SIM_REDSHIFT_HOST:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_REDSHIFT_HOST ) ; }

      if ( strcmp(c_get,"SIM_REDSHIFT_FLAG:")==0 ) 
	{ readint ( fp, 1, &SNDATA.SIM_REDSHIFT_FLAG ) ; }

      if ( strcmp(c_get,"SIM_VPEC:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_VPEC ) ; }

      if ( strcmp(c_get,"SIM_DLMU:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_DLMU ) ; }

      if ( strcmp(c_get,"SIM_LENSDMU:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_LENSDMU ) ; }

      if ( strcmp(c_get,"SIM_RA:")==0 ) 
	readfloat ( fp, 1, &SNDATA.SIM_RA ) ;
      if ( strcmp(c_get,"SIM_DEC:")==0 ) 
	readfloat ( fp, 1, &SNDATA.SIM_DEC ) ;
      if ( strcmp(c_get,"SIM_PEAKMJD:")==0 ) 
	readfloat ( fp, 1, &SNDATA.SIM_PEAKMJD ) ;

      if ( strcmp(c_get,"SIM_MWEBV:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_MWEBV ) ; }

      if ( strcmp(c_get,"SIMOPT_MWCOLORLAW:")==0 ) 
	{ readint ( fp, 1, &SNDATA.SIMOPT_MWCOLORLAW ) ; }
      if ( strcmp(c_get,"SIMOPT_MWEBV:")==0 ) 
	{ readint ( fp, 1, &SNDATA.SIMOPT_MWEBV ) ; }

      if ( strcmp(c_get,"SIM_AVTAU:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_AVTAU ) ; }
      if ( strcmp(c_get,"SIM_AV:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_AV ) ; }
      if ( strcmp(c_get,"SIM_RV:")==0 ) 
	{ readfloat ( fp, 1, &SNDATA.SIM_RV ) ; }

      if ( strcmp(c_get,"SIM_GALFRAC:") == 0 ) {
	fptr = &SNDATA.SIM_GALFRAC[MINFILT_DEF] ;
	readfloat ( fp, NFILT_DEF, fptr );
      }

      if ( strcmp(c_get,"SIM_PEAKMAG:")==0 ) {
	fptr = &SNDATA.SIM_PEAKMAG[MINFILT_DEF] ;
	readfloat ( fp, NFILT_DEF, fptr );
      }

      if ( strcmp(c_get,"SIM_EXPOSURE:")==0 ) {
	fptr = &SNDATA.SIM_EXPOSURE_TIME[MINFILT_DEF] ;
	readfloat ( fp, NFILT_DEF, fptr ); 
      }

      if ( strcmp(c_get,"SIM_STRETCH:")==0 ) 
	readfloat ( fp, 1, &SNDATA.SIM_STRETCH ) ;
      if ( strcmp(c_get,"SIM_DELTA:")==0 ) 
	readfloat ( fp, 1, &SNDATA.SIM_DELTA ) ;
      if ( strcmp(c_get,"SIM_DM15:")==0 ) 
	readfloat ( fp, 1, &SNDATA.SIM_DM15 ) ;

      if ( strcmp(c_get,"SIM_NON1a:")==0 ) 
	readint ( fp, 1, &SNDATA.SIM_TEMPLATE_INDEX ) ;
      if ( strcmp(c_get,"SIM_NONIa:")==0 )     // Mar 2013
	readint ( fp, 1, &SNDATA.SIM_TEMPLATE_INDEX ) ;
      if ( strcmp(c_get,"SIM_NONIA:")==0 )    // Mar 2013
	readint ( fp, 1, &SNDATA.SIM_TEMPLATE_INDEX ) ;
      if ( strcmp(c_get,"SIM_TEMPLATE_INDEX:")==0 )    // 7.31.2018
	readint ( fp, 1, &SNDATA.SIM_TEMPLATE_INDEX ) ;

    }  // end of fake > 0 if-block


	// get global epoch info

    if ( strcmp(c_get,"EPOCH:")==0 ) {
      //	  readint ( fp, 1, &NEWMJD );  // read EPOCH number


	  // increment epoch instead of reading it ...
	  // we may have to exclude some epochs so use local index

      EPMIN = epoch + 1; 
      NEWMJD++ ;
      NFILT_NEWMJD = 0;
      if ( NEWMJD >= MXEPOCH ) {
	sprintf(c1err,"MEWMJD %d  exceeds MXEPOCH=%d", NEWMJD, MXEPOCH );
	sprintf(c2err,"Check file: %s", inFile);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      
      SNDATA.NEWMJD         = NEWMJD;
      SNDATA.EPOCH_RANGE_NEWMJD[NEWMJD][0] = EPMIN;

    }

    /* xxxxxxxx mark delete Jan 23 2018 xxxxxxxxxxxxx
    if ( strcmp(c_get,"TELESCOPE:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"TELESCOPE:");
      readchar ( fp, SNDATA.TELESCOPE[EPMIN] ) ;
      SNDATA.IDTEL[EPMIN] = IDTELESCOPE( SNDATA.TELESCOPE[EPMIN] ) ;
    }
    xxxxxxxxxxxxxxxxxxx */


    if ( strcmp(c_get,"SEARCH_RUN:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"SEARCH_RUN:");
      readint ( fp, 1, &SNDATA.SEARCH_RUN[EPMIN] ) ;
    }

    if ( strcmp(c_get,"TEMPLATE_RUN:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"TEMPLATE_RUN:");
      readint ( fp, 1, &SNDATA.TEMPLATE_RUN[EPMIN] ) ;
    }
    if ( strcmp(c_get,"QMASK:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"QMASK:");
      readint ( fp, 1, &SNDATA.QMASK[EPMIN] ) ;
    }


    if ( strcmp(c_get,"MJD:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"MJD:");
      readdouble ( fp, 1, &SNDATA.MJD[EPMIN] ) ;
    }

    /* xxxxxxxxxxxxxxxxx mark delete xxxxxxxxxx
    if ( strcmp(c_get,"COLUMN:")==0 ) {  // obsolete legacy name for SDSS
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"COLUMN:");
      readint ( fp, 1, &SNDATA.IDCCD[EPMIN] ) ;
    }

    if ( strcmp(c_get,"IDCCD:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"COLUMN:");
      readint ( fp, 1, &SNDATA.IDCCD[EPMIN] ) ;
    }
    xxxxxxxxxx end delete xxxxxxxxxxxx */

    if ( strcmp(c_get,"STRIPE:")==0 ) {  // obsolete legacy name for SDSS
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"STRIPE:");
      readchar ( fp, SNDATA.FIELDNAME[EPMIN] ) ;
    }
    if ( strcmp(c_get,"FIELD:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"STRIPE:");
      readchar ( fp, SNDATA.FIELDNAME[EPMIN] ) ;
    }


    /* xxxxxxxxx mark delete June 24 2019 xxxxxxxxxxx
    if ( strcmp(c_get,"CLOUDCAM_SIG:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"CLOUDCAM_SIG:");
      readfloat ( fp, 1, &SNDATA.CLOUDCAM_SIG[EPMIN] ) ;
    }
    if ( strcmp(c_get,"CLOUDCAM_AVG:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"CLOUDCAM_AVG:");
      readfloat ( fp, 1, &SNDATA.CLOUDCAM_AVG[EPMIN] ) ;
    }


    if ( strcmp(c_get,"MOONPHASE:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"MOONPHASE:");
      readfloat ( fp, 1, &SNDATA.MOONPHASE[EPMIN] ) ;
    }
    if ( strcmp(c_get,"MOONDIST:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"MOONDIST:");
      readfloat ( fp, 1, &SNDATA.MOONDIST[EPMIN] ) ;
    }

    if ( strcmp(c_get,"AIRMASS:")==0 ) {
      if ( NEWMJD == 0 ) parse_err(inFile,NEWMJD,"AIRMASS:");
      readfloat ( fp, 1, &SNDATA.AIRMASS[EPMIN] ) ;
    }
    xxxxxxxxxxxxxxxxxxxxx */

    // ------------------------------------
    // FILTER-DEPENDENT information


    if ( strcmp(c_get,"PASSBAND:")==0  || strcmp(c_get,"BAND")==0 ) {

      if ( fgets(line_passband, 60, fp) == NULL ) { continue; }

      //      printf("\n\n  LINE_PASSBAND = '%s' \n", line_passband);

      ptrtok = strtok(line_passband," ") ; // split string

      // skip blank spaces (NULL) and new-line (10)

      while ( ptrtok != NULL && ptrtok[0] != 10 ) {
	sprintf(c_filt, "%c", ptrtok[0] );

	NTMP = PARSE_FILTLIST(c_filt, ifilt_tmp );

	//   printf("  FOUND filter %s = %d \n", c_filt, ifilt_tmp[0] );
	epoch++;	EPMAX = epoch ;
	NFILT_NEWMJD++ ;
	SNDATA.FILTINDX[epoch] = ifilt_tmp[0];
	SNDATA.NEPOCH          = epoch ;

	SNDATA.IDTEL[epoch]       =  SNDATA.IDTEL[EPMIN] ;
	SNDATA.MJD[epoch]         =  SNDATA.MJD[EPMIN] ;
	// xxx mark delete   SNDATA.IDCCD[epoch]  =  SNDATA.IDCCD[EPMIN] ;

	ptrtok = strtok(NULL, " ");
      }

      SNDATA.EPOCH_RANGE_NEWMJD[NEWMJD][1] = EPMAX ;

    }


    // now read filter-band info using 'rd_filtband' utility.

    if ( strcmp(c_get,"SEARCH_FIELD:")==0 ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"SEARCH_FIELD:");
      iptr = &SNDATA.SEARCH_FIELD[EPMIN] ;
      readint ( fp, NFILT_NEWMJD, iptr );
    }

    if ( strcmp(c_get,"TEMPLATE_FIELD:")==0 ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"TEMPLATE_FIELD:");
      iptr = &SNDATA.TEMPLATE_FIELD[EPMIN];
      readint ( fp, NFILT_NEWMJD, iptr );      
    }

    if ( strcmp(c_get,"PHOTOMETRYFLAG:")==0 ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"PHOTOMETRYFLAG:");
      iptr = &SNDATA.PHOTFLAG[EPMIN];
      readint ( fp, NFILT_NEWMJD, iptr );      
    }


    if ( strcmp(c_get,"GAIN:")==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"GAIN:");
      fptr = &SNDATA.GAIN[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );      
    }

    if ( strcmp(c_get,"RDNOISE:")==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"RDNOISE:");
      fptr = &SNDATA.READNOISE[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
    }

    if ( strcmp(c_get,"XPIXEL:")==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"XPIXEL:");
      fptr = &SNDATA.XPIX[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
    }
    if ( strcmp(c_get,"YPIXEL:")==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"YPIXEL:");
      fptr = &SNDATA.YPIX[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
    }


	// conditions
    if ( strcmp(c_get,"SKY_SIG:")==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"SKY_SIG:");
      fptr = &SNDATA.SKY_SIG[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
    }


    if ( strcmp(c_get,"PSF_SIG1:")==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"PSF_SIG1:");
      fptr = &SNDATA.PSF_SIG1[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
    }

    if ( strcmp(c_get,"PSF_SIG2:")==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"PSF_SIG2:");
      fptr = &SNDATA.PSF_SIG2[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
    }

    if ( strcmp(c_get,"PSF_RATIO:")==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,"PSF_RATIO:");
      fptr = &SNDATA.PSF_RATIO[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
    }





    /* xxxxxxxxxxx mark delete Jun 24 2019 xxxxxxxxxxxx
    // read fluxes
    sprintf(varname,"FLUX:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.FLUX[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -1.0E5, fluxmax );
      readchar(fp, FLUXUNIT );  // Oct 6, 2008
    }
    sprintf(varname,"FLUX_ERRTOT:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.FLUX_ERRTOT[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, fluxerrmax );
    }

    xxxxxxxxxxxxxx */

    // read CALIBRATED fluxes

    sprintf(varname,"FLUXCAL:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.FLUXCAL[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -1.0E4, fluxmax );
    }



    sprintf(varname,"%s","FLUXCAL_ERRTOT:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.FLUXCAL_ERRTOT[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, fluxerrmax );
    }


	// read magnitudes

    sprintf(varname,"%s", "MAG:" );
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.MAG[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -20.0, 1.0E4 );
    }
    sprintf(varname,"MAG_ERRPLUS:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.MAG_ERRPLUS[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, 2.0E2 );
    }

    sprintf(varname,"MAG_ERRMINUS:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.MAG_ERRMINUS[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, 2.0E2 );
    }


	// read zero points

    sprintf(varname,"ZEROPT:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.ZEROPT[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, 50.0 );
    }

    sprintf(varname,"ZEROPT_ERR:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.ZEROPT_ERR[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, 20.0 );
    }

    sprintf(varname,"ZEROPT_SIG:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.ZEROPT_SIG[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, 20.0 );
    }


    sprintf(varname,"SKYSUB_ERR:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.SKYSUB_ERR[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, 100000.0 );
    }

    sprintf(varname,"GALSUB_ERR:");
    if ( strcmp(c_get,varname)==0  ) {
      if ( NFILT_NEWMJD == 0 ) parse_err(inFile,NEWMJD,varname);
      fptr = &SNDATA.GALSUB_ERR[EPMIN];
      readfloat ( fp, NFILT_NEWMJD, fptr );
      checkval_F(varname, NFILT_NEWMJD, fptr, -10.0, 50000.0 );
    }

	// check for sim stuff (only if fake > 0 )
	if ( SNDATA.FAKE > 0 ) {

	  if ( strcmp(c_get,"SIM_MAG:")==0 ) 
	    readfloat ( fp, 1, &SNDATA.SIMEPOCH_MAG[EPMIN] ) ;

	  if ( strcmp(c_get,"SIM_MAGSMEAR:")==0 ) 
	    readfloat ( fp, 1, &SNDATA.SIMEPOCH_MAGSMEAR[EPMIN] ) ;

	  if ( strcmp(c_get,"SIM_TREST:")==0 ) 
	    readfloat ( fp, 1, &SNDATA.SIMEPOCH_TREST[EPMIN] ) ;

	  if ( strcmp(c_get,"SIM_TOBS:")==0 ) 
	    readfloat ( fp, 1, &SNDATA.SIMEPOCH_TOBS[EPMIN] ) ;

	  if ( strcmp(c_get,"SIM_WARPCOL_SYMBOL:")==0 ) {
	    for ( eptmp=EPMIN; eptmp <= EPMAX; eptmp++ )
	      readchar(fp,SNDATA.SIMEPOCH_WARPCOLNAM[eptmp] );
	  }

	  if ( strcmp(c_get,"SIM_WARPCOL_VALUE:")==0 ) {
	    fptr = &SNDATA.SIMEPOCH_WARPCOLVAL[EPMIN];
	    readfloat ( fp, NFILT_NEWMJD, fptr );
	  }

	  if ( strcmp(c_get,"SIM_KCOR_SYMBOL:")==0 ) {
	    for ( eptmp=EPMIN; eptmp <= EPMAX; eptmp++ )
	      readchar(fp,SNDATA.SIMEPOCH_KCORNAM[eptmp]  );
	  }

	  if ( strcmp(c_get,"SIM_KCOR_VALUE:")==0 ) {
	    fptr = &SNDATA.SIMEPOCH_KCORVAL[EPMIN];
	    readfloat ( fp, NFILT_NEWMJD, fptr );
	  }

	  if ( strcmp(c_get,"SIM_AVWARP:")==0 ) {
	    fptr = &SNDATA.SIMEPOCH_AVWARP[EPMIN];
	    readfloat ( fp, NFILT_NEWMJD, fptr );
	  }

	}

  }  // end of c_get fscan loop



  printf("  Done. \n");

  fclose(fp);

  return SUCCESS;

}  // end of rd_SNDATA



// *************************************************************
void read_redshift(FILE *fp, float *redshift, float *redshift_err ) {
  // read redshift and error separated by "+-"
  char cdum[10];
  // --------------- BEGIN ----------
  readfloat(fp, 1, redshift );
  readchar(fp, cdum);                  // read "+-" symbol
  readfloat(fp, 1, redshift_err );
}


// *****************************************
int PARSE_FILTLIST (char *filtlist_string, int *filtlist_array ) {

  // June 27 2008
  // parse input string 'filtlist_string' and return integer-array
  // "filtlist_array" of absolute filter indices.
  // Function arg is the number of filters.

  char fnam[] = "PARSE_FILTLIST";

  int LENLIST, NF_USER ;
  int ifilt_user, ifilt_list, ifilt_match;
  char cfilt_user[2], cfilt_list[2];

  //------------- BEGIN ----------


  LENLIST = strlen(FILTERSTRING);
  NF_USER = strlen(filtlist_string);

  if ( NF_USER >= LENLIST ) {
    sprintf(c1err,"%d filters is too many (> %d)", NF_USER, LENLIST);
    sprintf(c2err,"Check filter list  '%s' ", filtlist_string);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // loop over all user-defined filters to make sure
  // that they are all defined in FILTERSTRING.

  for ( ifilt_user = 0 ; ifilt_user < NF_USER; ifilt_user++ ) {

    sprintf(cfilt_user, "%c", *(filtlist_string+ifilt_user) ) ;

    ifilt_match = -9 ;

    for ( ifilt_list = 0 ; ifilt_list < LENLIST; ifilt_list++ ) {
      sprintf(cfilt_list, "%c", FILTERSTRING[ifilt_list] ) ;
      if ( strcmp(cfilt_user,cfilt_list) == 0 && ifilt_match < 0 ) {
	ifilt_match = ifilt_list;
	*(filtlist_array+ifilt_user) = ifilt_match ;
      }
    }

    if ( ifilt_match < 0 ) {
      sprintf(c1err,"User-defined filter '%s' is not defined.", cfilt_user);
      sprintf(c2err,"Check filter list  %s' ", filtlist_string);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  }


  return NF_USER ;

}  // end of function PARSE_FILTLIST



/* xxxxxxxx mark delete Jan 5 2018 xxxxxxxxxxx
float snfluxmax ( float z ) {

  // returns max possible SN flux for redsfhit =z .
  // Used for idiot checks on flux.

  double ztmp8,  magtmp8, fluxmax8 ;
  double H0 = 65.0 / ( 1.0E6 * PC_km) ;
  double W0 = -1.0 ;
  double OM =  0.3 ;
  double OE =  0.7;


  ztmp8   = (double)z ;
  //  magtmp8 = dLmag( H0, OM, OE, W0, ztmp8 ) - 19.4 ;
  fluxmax8 = 1.0E12/powf((double)10.0,0.4*magtmp8);

  fluxmax8 *= 2.0;  // Apr 16, 2008

  //  printf("\n xxxx z=%f  => magtmp=%f  fluxmax=%f \n", ztmp8, magtmp8, fluxmax);

  return  (float)fluxmax8 ;

} // end  of snfluxmax
xxxxxxxxxxxx */

// ************************************************
void checkArrayBound(int i, int MIN, int MAX, 
		     char *varName, char *comment, char *funName) {

  // One-line call to check array bound.
  // Abort if index 'i' is outside index range MIN to MAX.
  // *varName, *comment and *fun are used to construct 
  // error message. Note that *fun is the name of the calling function.
  //

  char fnam[] = "checkArrayBound" ;

  // ------------------ BEGIN ---------------

  if ( (i >= MIN) && (i <= MAX) ) { return ; }

  sprintf(c1err,"%s = %d is outside valid range %d - %d (fun=%s) \n", 
	  varName ,i, MIN, MAX, funName);
  sprintf(c2err,"%s", comment);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  

} // end of checkArrayBound

void  checkArrayBound_(int *i, int *MIN, int *MAX, 
		       char *varName, char *comment, char *funName) {
  checkArrayBound(*i,*MIN,*MAX, varName, comment, funName) ;
}


// ******************************************************
void read_VARNAMES_KEYS(FILE *fp, int MXVAR, int NVAR_SKIP, char *callFun, 
			int *NVAR, int *NKEY,
			int *UNIQUE, char **VARNAMES ) {

  // Mar 2019
  // Read file (*fp) and return information about "VARNAMES: <varList>"
  //
  // Inputs
  //   fp      : file pointer to read
  //    MXVAR  : max number of variables to return after VARNAMES keys.
  //    NVAR_SKIP : number of variables to skip at end of list
  //  *callFun : name of calling function; for error message only
  //
  // Output:
  //   *NVAR     : total number of variables after all VARNAMES keys
  //   *NKEY     : total number of VARNAMES keys found
  //   *UNIQUE   : for each variable, 1=> unique, 0=> duplicate from previous
  //  **VARNAMES : list of all variables (0 to *NVAR-1)

  int  NVAR_LOCAL = 0 ;
  int  NKEY_LOCAL = 0 ;
  int  FOUND_VARNAMES, IVAR, ivar, ivar2, NVAR_TMP ;
  char c_get[60], LINE[100] ;
  char fnam[] = "read_VARNAMES_KEYS" ;

  // -------------- BEGIN ------------

  while( (fscanf(fp, "%s", c_get )) != EOF) {
    FOUND_VARNAMES = ( strcmp(c_get,"VARNAMES:")==0 ) ;
    if ( FOUND_VARNAMES ) {
      NKEY_LOCAL++ ;
      fgets(LINE, 100, fp ); // scoop up varnames
      NVAR_TMP  = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
      NVAR_TMP -= NVAR_SKIP ;
      for ( ivar=0; ivar < NVAR_TMP; ivar++ ) {
	IVAR = ivar+NVAR_LOCAL ;
	if ( IVAR < MXVAR ) { get_PARSE_WORD(0,ivar,VARNAMES[IVAR]); }
      }
      NVAR_LOCAL += NVAR_TMP ;
    } // end FOUND_VARNAMES
  } // end while    



  if ( NVAR_LOCAL > MXVAR ) {
    sprintf(c1err,"NVAR=%d exceeds MXVAR=%d", NVAR_LOCAL, MXVAR);
    sprintf(c2err,"called by %s, VARNAMES[0]=%s", callFun, VARNAMES[0] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // ------------------------------------
  // load output args
  *NVAR = NVAR_LOCAL ;
  *NKEY = NKEY_LOCAL ;

  // check which VARNAMES are unique
  char *NAME, *NAME2 ;
  for(ivar=0; ivar < NVAR_LOCAL; ivar++ ) {
    UNIQUE[ivar] = 1;
    NAME = VARNAMES[ivar] ;

    // loop over previous VARNAMES to check for duplicate
    for(ivar2=0; ivar2 < ivar; ivar2++ ) {
      NAME2 = VARNAMES[ivar2];
      if ( strcmp(NAME,NAME2) == 0 ) { UNIQUE[ivar]=0; } // duplicate
    }
  }

  
  return ;

} // end read_VARNAMES_KEYS

// **********************
void check_uniform_bins(int NBIN,double *VAL_ARRAY, char *comment_forAbort) {

  // July 2016: abort on non-uniform bin in *VAL array

  int i, j ;
  double VAL, VAL_LAST, DIF, DIF_LAST ;
  char fnam[] = "check_uniform_bins" ;

  // ------------ BEGIN ------------
  
  VAL_LAST = DIF_LAST = 0.0 ;
  for(i=0; i < NBIN; i++ ) {

    VAL = VAL_ARRAY[i];
    DIF = VAL - VAL_LAST ;

    if ( i > 1 && DIF != DIF_LAST ) {
      printf("\n %s PRE-ABORT DUMP for %s: \n", fnam, comment_forAbort);
      for(j=i-2; j < i+3; j++ ) {
	if ( j<0 || j >= NBIN ) { continue ; }
	printf("\t VAL_ARRAY[%d] = %f  (DIF=%f) \n", 
	       j, VAL_ARRAY[j], VAL_ARRAY[j]-VAL_ARRAY[j-1] );
      }
      sprintf(c1err,"Detected non-uniform grid (%s)", comment_forAbort);
      sprintf(c2err,"DIF=%f but DIF_LAST=%f", DIF, DIF_LAST);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    VAL_LAST = VAL;
    DIF_LAST = DIF ;
  }
  
  return ;

} // end check_uniform_bins

// ******************************************************
void  check_magUndefined(double mag, char *varName, char *callFun) {

  // Created Jun 23 2016
  char fnam[] = "check_magUndefined" ;

  if ( mag == MAG_UNDEFINED ) {
    sprintf(c1err,"Undefined %s = %f", varName, mag);
    sprintf(c2err,"Found in function %s", callFun );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );   
  }
}  // end check_magUndefined


// ******************************************************
void checkval_I(char *varname, int nval, int *iptr, int imin, int imax){

  // check that all iptr values are between imin and imax;
  // if not then abort. 
  // Jan 11 2017: isnan(ival) -> isnan ( float)ival )

  int i;
  int  ival ;
  char fnam[] = "checkval_I" ;

  // ----------------- BEGIN -----------------

  for ( i = 0; i < nval;  i++ ) {
    ival = *(iptr+i);
    
    if ( isnan( (float)ival ) ) {
      sprintf(c1err,"%s(item %d) = nan", varname, i );
      sprintf(c2err,"Expected  %i < %s < %i", imin, varname, imax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( ival < imin || ival > imax  ) {
      sprintf(c1err,"%s(item %d) = %d fails bounds-check", varname, i, ival );
      sprintf(c2err,"Expected  %d < %s < %d", imin, varname, imax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

}  // end of checkval_I


void checkval_i__(char *varname, int *nval, int *iptr, int *imin, int *imax){
  checkval_I(varname, *nval, iptr, *imin, *imax) ;
}


// ******************************************************
void checkval_F(char *varname, int nval, float *fptr, float fmin, float fmax){

  // check that all fptr values are between fmin and fmax;
  // if not then abort. 
  // Nov 2013: check for nan and rename function checkval -> checkval_F
  // May 27 2014: isnanf -> isnan (recommended by S. Rodney)

  int i;
  float val ;
  char fnam[] = "checkval_F" ;

  // ----------------- BEGIN -----------------

  for ( i = 0; i < nval;  i++ ) {
    val = *(fptr+i);
    
    if ( isnan(val) ) {
      sprintf(c1err,"%s(item %d) = nan", varname, i );
      sprintf(c2err,"Expected  %f < %s < %f", fmin, varname, fmax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( val < fmin || val > fmax  ) {
      sprintf(c1err,"%s(item %d) = %f fails bounds-check", varname, i, val );
      sprintf(c2err,"Expected  %f < %s < %f", fmin, varname, fmax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

}  // end of checkval_F

// ******************************************************
void checkval_D(char *varname, int nval, 
		double *dptr, double dmin, double dmax){

  // check that all fptr values are between fmin and fmax;
  // if not then abort. 
  // Nov 2013: check for nan and rename function checkval -> checkval_F

  int i;
  double val ;
  char fnam[] = "checkval_D" ;

  // ----------------- BEGIN -----------------

  for ( i = 0; i < nval;  i++ ) {
    val = *(dptr+i);
    
    if ( isnan(val) ) {
      sprintf(c1err,"%s(item %d) = nan", varname, i );
      sprintf(c2err,"Expected  %f < %s < %f", dmin, varname, dmax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( val < dmin || val > dmax  ) {
      sprintf(c1err,"%s(item %d) = %f fails bounds-check", varname, i, val );
      sprintf(c2err,"Expected  %f < %s < %f", dmin, varname, dmax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

}  // end of checkval_D

// ********************************
int landolt_ini__(int *opt, float *mag, float *kshift) {  
  // fortran wrapper
  int istat;
  istat = Landolt_ini(*opt, mag, kshift);
  return istat;
}

// ********************************
int Landolt_ini(
		int opt                // (I) BD17 or Vega
		,float *mag            // (I) primary mag  to store
		,float *kshift         // (I) shift k's for systematics
		) {

  // Read Landolt transformation parameters from file
  // and store primary mags in global var.
  // Landolt UBVRI,BX <=> 012345
  //
  // Jan 2, 2009: allow k_i values up to 1.0 instead of 0.1
  //              to allow for filter-adjustment tests that
  //              have larger color-transformations

  int ifilt, k ;

  float kval, kerr,  magtmp ;

  char 
    fnam[] = "Landolt_ini" 
    ,c_get[40]
    ,c_tmp[60]
    ,c_k[6]
    ,kfile[40]
    ,kfile_full[120]
    ;

  FILE *fp;

  // --------- BEGIN -------------


  print_banner("INIT  BESSELL <=> LANDOLT  TRANSFORMATIONS" );


  // init color terms to crazy value.

  for ( k=0; k <= 4 ; k++ ) {
    LANDOLT_COLOR_VALUE[k] = 0.0 ;
    LANDOLT_COLOR_ERROR[k] = 0.0 ;
  }

  // store mag in global array

  printf("   UBVRI,BX offsets: ");
  for ( ifilt=0; ifilt < NFILT_LANDOLT; ifilt++ ) {
    magtmp = *(mag + ifilt) ;
    LANDOLT_MAGPRIMARY[ifilt] = (double)magtmp ;
    printf(" %7.3f", magtmp );
  }
  printf("\n\n");

  if ( opt == 0 ) 
    goto PRINT_COLOR_TERMS ;
  else if ( opt < 4 ) 
    sprintf(kfile, "LANDOLT_COLOR_TERMS_BD17.DAT" );
  else
    sprintf(kfile, "LANDOLT_COLOR_TERMS_VEGA.DAT" );


  // read color terms from file
  // First try opening file in user dir ;
  // if not there than open official file in $SNDATA_ROOT.

  sprintf(c_tmp,"Ready to read Landolt color terms from");

  if ( (fp = fopen(kfile, "rt")) != NULL ) {
    printf("   %s user k-file: \n   %s \n", c_tmp, kfile );
    goto READ_KFILE ;
  }

  sprintf( PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );
  sprintf(kfile_full,"%s/filters/Landolt/%s", 
	  PATH_SNDATA_ROOT, kfile );

  if ( (fp = fopen(kfile_full, "rt"))==NULL ) {
    sprintf(c1err,"Cannot open: %s", kfile_full );
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
    fclose(fp); 
  }

  printf("   %s default k-file: \n   %s \n", c_tmp, kfile_full );

 READ_KFILE:

  while( fscanf(fp,"%s", c_get )!= EOF ) {

    for ( k=0; k <= 4; k++ ) {
      sprintf(c_k, "k%d:", k );
      if ( strcmp(c_get,c_k) == 0 ) {

	readfloat(fp, 1, &kval );
	LANDOLT_COLOR_VALUE[k] = (double)kval + *(kshift+k);

	readchar(fp, c_tmp);  // skip over '+-' symbol

	readfloat(fp, 1, &kerr );
	LANDOLT_COLOR_ERROR[k] = (double)kerr ;

      }
    }
  } // end of fscanf loop


  fclose(fp);


  // print color terms k0 - k4

 PRINT_COLOR_TERMS:

  for ( k=0; k <=4 ; k++ ) {
    kval = LANDOLT_COLOR_VALUE[k] ;
    kerr = LANDOLT_COLOR_ERROR[k] ;

    printf("\t   k%d = %6.3f +- %6.3f \n", k, kval, kerr ) ;

    if ( fabsf(kval) > 1.0 ) {
      sprintf(c1err,"k%d has invalid value above", k);
      errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
    }
  }

  printf("\t   U(reported) -> [UX - BX + B]_(synthetic) \n");


  return SUCCESS;

}  // end 




/**********************************************
  SALT-II color correction formula
**********************************************/
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
  double 
    arg
    ,lr, lr2, lr3
    ,numerator 
    ,denominator 
    ,CLAM
    ;

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
  // Mar 3, 2011: fix dumb bug causing crazy value at lam = LAM_MIN

  // define local args for *colorPar inputs
  double 
    LAM_B, LAM_V
    ,LAM_MIN, LAM_MAX
    ,XN, params[10]
    ;

  int nparams  ;

  double constant = log(10.0)/2.5 ;
  double alpha    = 1 ;
  double val      = 0 ;

  double rl, rlmin, rlmax, tmp ;
  int i;

  char fnam[] = "SALT2colorlaw1" ;

  // --------------- BEGIN ------------

  sprintf(c2err,"Check colorlaw parameters");

  // strip of colorPar values into local variables.
  LAM_B   = *(colorPar+0);   // mean lambda of B filter
  LAM_V   = *(colorPar+1);   // mean labmda of V filter
  LAM_MIN = *(colorPar+2);   // special extrap function below this
  LAM_MAX = *(colorPar+3);   // idem for upper wavelength
  XN      = *(colorPar+4);   // Number of parameters for function

  nparams = (int)XN ;        // number of parameters to follow
  for ( i=0; i < nparams; i++ ) {
    tmp       = *(colorPar+5+i); 
    params[i] = tmp;
    if ( fabs(tmp) > 10.0 ) {
      sprintf(c1err, "insane params[%d] = %f", i, tmp );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  // --------------------------------------
  // make a few sanity checks on passed parameters

  if ( LAM_B < 4000 || LAM_B > 4500 ) {
    sprintf(c1err, "insane LAM_B = %6.0f", LAM_B );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( LAM_V < 5000 || LAM_V >6000 ) {
    sprintf(c1err, "insane LAM_V = %6.0f", LAM_V );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  if ( LAM_MIN < 1000 || LAM_MIN > 6000 ) {
    sprintf(c1err, "insane LAM_MIN = %6.0f", LAM_MIN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  if ( LAM_MAX < 6000 || LAM_MAX > 16000 ) {
    sprintf(c1err, "insane LAM_MAX = %6.0f", LAM_MAX );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // ------------------------------

  for(i=0; i < nparams; i++)
    { alpha = alpha - params[i]; }

  // compute reduced wavelengths
  rl    = (lambda  - LAM_B) / ( LAM_V - LAM_B );
  rlmin = (LAM_MIN - LAM_B) / ( LAM_V - LAM_B );
  rlmax = (LAM_MAX - LAM_B) / ( LAM_V - LAM_B );

  if(lambda >= LAM_MIN  && lambda <= LAM_MAX ) {
    val = SALT2colorfun_pol(rl,nparams,params,alpha);

  }else{

    if(lambda < LAM_MIN ) {
      double Pmin  = SALT2colorfun_pol(rlmin,nparams,params,alpha);
      double dPmin = SALT2colorfun_dpol(rlmin,nparams,params,alpha);
      val =  Pmin + dPmin* (rl-rlmin);
    }else{
      double Pmax  = SALT2colorfun_pol(rlmax,nparams,params,alpha);
      double dPmax = SALT2colorfun_dpol(rlmax,nparams,params,alpha);
      val =  Pmax + dPmax* (rl-rlmax);

      /*
      printf(" xxx lambda=%6.0f rl=%5.2f rlmax=%5.2f alpha=%6.2f Pmax=%6.2f dPmax=%6.2f nparam=%d val=%f \n",
	     lambda, rl, rlmax, alpha, Pmax, dPmax, nparams, val );
      */
    }
  }



  /*
  if ( fabs(lambda-2271.) < 2.0 && abs(c+0.4) < 0.1 ) {
    printf(" xxx %s : lam=%6.1f c=%5.2f -> const=%f val=%f (LAM_MIN=%5.0f) \n",
	   fnam, lambda, c, constant, val, LAM_MIN);
    fflush(stdout);
  }
  */

  return  exp(c*constant*val);

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



// ********************************
// fortran wrapper
int landolt_convert__(int *opt, double *mag_in, double *mag_out) {
  int istat ;
  istat = Landolt_convert(*opt, mag_in, mag_out);
  return istat ;
}

// ********************************
int Landolt_convert(int opt, double *mag_in, double *mag_out) {

  /*****
   Created Mar 11, 2007 by R.Kessler

   opt = +1 : *mag_in  = synthetic Bessell UBVRI,BX
              *mag_out = reported  Landolt UBVRI 

   opt = -1 : *mag_in  = reported  Landolt UBVRI, BX-B
              *mag_out = synthetic Bessell UBVRI, BX

         (but note that reported U is (UX - BX + B)_synth


  ******/

  int ifilt;
  double k0, k1, k2, k3, k4 ;

  int off_U  = 0 ;
  int off_B  = 1 ;
  int off_V  = 2 ;
  int off_R  = 3 ;
  int off_I  = 4 ;
  int off_BX = 5 ;

  double del, delref;
  double DEL_V, DEL_BV, DEL_UBX, DEL_VR, DEL_RI;
  double Vout, Utmp, DEL_BXB ;

  // ------------ BEGIN ----------------

  // init *mag_out
  for ( ifilt=0; ifilt < NFILT_LANDOLT ; ifilt++ ) {
    *(mag_out+ifilt) = -99.0 ; 
  }


  k0 = LANDOLT_COLOR_VALUE[0];
  k1 = LANDOLT_COLOR_VALUE[1];
  k2 = LANDOLT_COLOR_VALUE[2];
  k3 = LANDOLT_COLOR_VALUE[3];
  k4 = LANDOLT_COLOR_VALUE[4];

  // apply magdif array

  if ( opt > 0 ) {  // convert Bessell -> Landolt


    del    = *(mag_in+off_B)       - *(mag_in+off_V);
    delref = LANDOLT_MAGPRIMARY[off_B] - LANDOLT_MAGPRIMARY[off_V] ;
    DEL_V  = k0*(del-delref);
    DEL_BV = k1 * (del-delref);

    del     = *(mag_in+off_U)       - *(mag_in+off_BX);
    delref  = LANDOLT_MAGPRIMARY[off_U] - LANDOLT_MAGPRIMARY[off_BX] ;
    DEL_UBX = k2 * (del-delref);

    del     = *(mag_in+off_V)       - *(mag_in+off_R);
    delref  = LANDOLT_MAGPRIMARY[off_V] - LANDOLT_MAGPRIMARY[off_R] ;
    DEL_VR  = k3 * (del-delref);

    del     = *(mag_in+off_R)       - *(mag_in+off_I);
    delref  = LANDOLT_MAGPRIMARY[off_R] - LANDOLT_MAGPRIMARY[off_I] ;
    DEL_RI  = k4 * (del-delref);

    *(mag_out + off_V) = *(mag_in + off_V) + DEL_V ;
    *(mag_out + off_B) = *(mag_in + off_B) + DEL_V + DEL_BV ;
    *(mag_out + off_R) = *(mag_in + off_R) + DEL_V - DEL_VR ;
    *(mag_out + off_I) = *(mag_in + off_I) + DEL_V - DEL_VR - DEL_RI ;

    *(mag_out + off_U) = *(mag_in+off_U) - *(mag_in+off_BX) + *(mag_in+off_B)
     + DEL_V + DEL_BV + DEL_UBX ;

  } 

  else if ( opt < 0 ) {  // convert Landolt -> Bessell

    del    = *(mag_in+off_B)       - *(mag_in+off_V) ;
    delref = LANDOLT_MAGPRIMARY[off_B] - LANDOLT_MAGPRIMARY[off_V] ;
    DEL_BV = (del + k1*delref) / ( 1. + k1 ) ;  // (B-V)_Bess
    DEL_V  = k0 * (delref - DEL_BV );  // V_Bess - V_Land

    del    = *(mag_in+off_V)       - *(mag_in+off_R) ;
    delref = LANDOLT_MAGPRIMARY[off_V] - LANDOLT_MAGPRIMARY[off_R] ;
    DEL_VR = (del + k3*delref) / ( 1. + k3 ) ;  // (V-R)_Bess

    del    = *(mag_in+off_R)       - *(mag_in+off_I) ;
    delref = LANDOLT_MAGPRIMARY[off_R] - LANDOLT_MAGPRIMARY[off_I] ;
    DEL_RI = (del + k4*delref) / ( 1. + k4 ) ;  // (R-I)_Bess

    del    = *(mag_in+off_U)       - *(mag_in+off_B) ;
    delref = LANDOLT_MAGPRIMARY[off_U] - LANDOLT_MAGPRIMARY[off_B] ;
    DEL_UBX = (del + k2*delref) / ( 1. + k2 );   // (UX-BX)_Bess


    Vout = *(mag_in  + off_V) + DEL_V ;
    *(mag_out+off_V) = Vout ;
    *(mag_out+off_B) = Vout + DEL_BV ;
    *(mag_out+off_R) = Vout - DEL_VR ;
    *(mag_out+off_I) = *(mag_out+off_R) - DEL_RI ;

    Utmp = DEL_UBX + *(mag_out+off_B) ;  // reported U = UX-BX+B

    // to get synthetic U, add synthetic BX-B

    DEL_BXB = *(mag_in+off_BX) ;
    *(mag_out+off_U) = Utmp + DEL_BXB ;

    // BX = (BX-B)_in + B_out
    *(mag_out+off_BX) = DEL_BXB + *(mag_out+off_B);

  }

  return SUCCESS ; // add Aug 7 2014 to avoid compile warning.

}  // end of function Landolt_convert



// *************************************
float edgedist ( float X, float Y, int NXPIX, int NYPIX ) {

  float XX, YY, DX, DY, DMIN;

  XX = (float)NXPIX - X;
  YY = (float)NYPIX - Y;      

  // take min by brute force since fminf does not work ???
  if ( X < XX ) DX = X; else DX = XX;
  if ( Y < YY ) DY = Y; else DY = YY;

  if ( DX < DY ) DMIN = DX; else DMIN = DY;

  return DMIN ;

} // end of edgedist


// ==================================
FILE *open_TEXTgz(char *FILENAME, const char *mode, int *GZIPFLAG ) {

  // Dec 1 2017:
  // Shell to call fopen.
  // If mode is read, then check for gzipped file.
  // If both gzip and unzip file exist, ABORT.
  // Return GZIPFLAG=1 if gzip file is opened
  //
  // Mar 6 2019: replace zcat with 'gunzip -c' so that it works on Mac.

  FILE *fp ;
  struct stat statbuf ;
  int istat_gzip, istat_unzip, LEN;
  char gzipFile[MXPATHLEN], unzipFile[MXPATHLEN];
  char cmd_zcat[MXPATHLEN] ;
  char fnam[]=  "open_TEXTgz" ;

  // -------------- BEGIN ------------

  *GZIPFLAG = 0 ;

  // for reading, check for gz file
  if ( strstr(mode,"r") != NULL ) {

    if ( strstr(FILENAME,".gz") != NULL )  { 
      // .gz file name passed as argument
      LEN = strlen(FILENAME);
      sprintf(gzipFile,   "%s", FILENAME); 
      sprintf(unzipFile,  "%s", FILENAME);
      unzipFile[LEN-3] = 0 ; // remove .gz extension
    } 
    else {
      // check for .gz file if not given
      sprintf(gzipFile,  "%s.gz", FILENAME);
      sprintf(unzipFile, "%s",    FILENAME);
    } 

    istat_gzip  = stat(gzipFile,  &statbuf) ;
    istat_unzip = stat(unzipFile, &statbuf) ;

    //    printf(" xxx ------------------------------- \n");
    //    printf(" xxx istat=%3d for '%s' \n", istat_gzip,  gzipFile);
    //    printf(" xxx istat=%3d for '%s' \n", istat_unzip, unzipFile);

    if ( istat_gzip==0 && istat_unzip==0 ) {
      printf("\n PRE-ABORT DUMP \n");
      printf("  Found %s \n", gzipFile );
      printf("  Found %s \n", unzipFile );
      sprintf(c1err, "Found gzipped and unzipped file."); 
      sprintf(c2err, "Cannot open both files.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    if ( istat_gzip == 0 ) {
      // xxx mark delete      sprintf(cmd_zcat, "zcat %s", gzipFile);
      sprintf(cmd_zcat, "gunzip -c %s", gzipFile);
      fp = popen(cmd_zcat,"r");
      *GZIPFLAG = 1 ;
      return(fp);
    }
  }

  // if we get here, do regular text open
  fp = fopen(FILENAME,mode);
  return(fp);

} // end open_TEXTgz


// =====================================
void snana_rewind(FILE *fp, char *FILENAME, int GZIPFLAG) {

  int gzipFlag ;
  char fnam[] = "snana_rewind" ;

  // --------------- BEGIN ----------------

  if ( GZIPFLAG == 0 ) {
    rewind(fp);
  }
  else {
    // cannot rewind popen for gz file; must close and re-open
    int istat = pclose(fp);

    if ( istat == -1 ) {
      sprintf(c1err,"pclose Error = -1 for file=");
      sprintf(c2err,"%s", FILENAME);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    fp = open_TEXTgz(FILENAME, "rt", &gzipFlag);
  }

} // end snana_rewind


// *************************************************
FILE *snana_openTextFile (int vboseFlag, char *SNPATH, char *fileName, 
			  char *fullName, int *gzipFlag ) {

  /* ----------------------------------------------
    Shell to open text file for reading.
    First search local directory; if file not found
    then search $SNPATH

    This allows user to easily over-ride default
    with private version.

   Function returns file pointer and gzipFile.

   Dec 29 2017: use open_TEXTgz to allow reading gzipped files.
   Jan 11 2018: add gzipFile output arg
   Mar 20 2019: padd vboseFlag to print comment to stdout
  ----------------------------------------------- */

#define TEXTMODE_read "rt"

  int LDMP = (vboseFlag>0) ;
  FILE *fp ;
  //  char fnam[] = "snana_openTextFile" ;

  // --------------- BEGIN ----------------

  // First try current working directory

  sprintf(fullName, "%s", fileName );

  //  printf("xxx %s : fileName = '%s' \n", fnam, fullName); // DDDDDDDDDD

  //  fp = fopen(fullName, "rt");
  fp = open_TEXTgz(fullName,TEXTMODE_read, gzipFlag );
  if ( fp != NULL ) {       
    if ( LDMP )  { printf("\t Opened : %s \n", fullName ); }
    return fp;
  } 


   // if we get here, try official location
   sprintf(fullName, "%s/%s", SNPATH,  fileName );
   //   fp = fopen(fullName, "rt") ;
   fp = open_TEXTgz(fullName,TEXTMODE_read, gzipFlag );

   if ( LDMP && fp != NULL ) { printf("\t Opened : %s \n", fullName ); }

   // return pointer regardless of status
   return fp;

}  // end of snana_openTextFile



// ***************************
int INTFILTER ( char *cfilt ) {

  // returns absolute filter index  for string *cfilt
  // Oct 29 2019: use last char if cfilt to work with arbitrary string

  int len = strlen(cfilt);
  int ifilt, itmp;
  char ctmp[2], cfilt1[2];
  //---------- BEGIN ----------------

  sprintf(cfilt1, "%c", cfilt[len-1]);
  ifilt = 0 ;
  for ( itmp=0; itmp < MXFILTINDX; itmp++ ) {
    sprintf(ctmp,"%c", FILTERSTRING[itmp] );
    if (strcmp(ctmp,cfilt1) == 0 ) { ifilt = itmp; return ifilt ; }
  }

  return ifilt;
} // end of INTFILTER



// *****************************************************
int wr_filtband_int ( 
		     FILE *fp        // write to this  *fp
		     ,char *keyword  //keyword before integers
		     ,int NINT       // number of integers to write
		     ,int *iptr      // pointer to first integer
		     ,char *comment  // comment after integers
		     ,int opt        // 0=> header format; 1=> BAND format
		     ) {

  int i;

  // ------------- BEGIN -----------

  if ( opt == 1 ) {
    fprintf(fp,"  %16s  ", keyword);

    for ( i = 0; i < NINT ; i++ ) 
      fprintf(fp, "%8d ", *(iptr+i) );
  }
  else {
    fprintf(fp,"%s  ", keyword);

    for ( i = 0; i < NINT ; i++ ) 
      fprintf(fp, "%4d ", *(iptr+i) );
  }

  fprintf(fp," %s \n" , comment ) ;

  return(0); // added Aug 7 2014 to avoid compile warning

} // end of wr_filtband_int


// *****************************************************
int wr_filtband_float(
		      FILE *fp 
		      ,char *keyword  // keyword
		      ,int NFLOAT     // number of floats to write
		      ,float *fptr    // point to 1st float
		      ,char *comment  // comment after floats
		      ,int idec       // number of digits after decimal
     ) {

  int i;

  // ---------------- BEGIN --------------


  if ( idec > 0 ) 
    fprintf(fp,"  %16s  ", keyword  ) ;  // for epoch
  else
    fprintf(fp,"%20s", keyword  ) ;  // for header

  for ( i=0; i<NFLOAT; i++ ) {
    if ( idec == 2 ) 
      fprintf(fp,"%8.2f ", *(fptr+i) );
    else if ( idec == 3 ) 
      fprintf(fp,"%8.3f ", *(fptr+i) );
    else if ( idec == 4 ) 
      fprintf(fp,"%8.4f ", *(fptr+i) );
    else if ( idec == 5 ) 
      fprintf(fp,"%8.5f ", *(fptr+i) );

    else if ( idec < 0  ) 
      fprintf(fp,"%6.1f ", *(fptr+i) );
  }

  fprintf(fp," %s \n" , comment ) ;

  return(0); // add Aug 7 2014 to avoid compile warnings.

} // end of wr_filtband_float


// ***********************************
void check_argv(void) {

  // make sure that there are no unused command-line args

  int NBAD, i ;

  // ----------- BEGIN ---------

  NBAD = 0;
  for ( i = 1; i < NARGV_LIST; i++ ) {
    if ( USE_ARGV_LIST[i] == 0 ) {
      NBAD++ ;
      if ( NBAD == 1 ) printf("  CHECK_ARGV ERRORS: \n" );

      printf(" \t ERROR: detected un-used command-line arg: %s \n", 
	     ARGV_LIST[i]);
    }
  }

  if ( NBAD > 0 )  madend(1);

} // end check_argv

// *******************************************************
void parse_err ( char *inFile, int NEWMJD, char *keyword ) {

  // print standard error message for parsing.
  // errmsg called with ABORT flag !

  char fnam[] = "parse_err" ;

  sprintf(c1err,"Problem parsing %s", inFile );
  sprintf(c2err,"Check NEWMJD-Epoch %d, keyword '%s' ", NEWMJD, keyword);

  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
}

// *********************************************************
void  legacyKey_abort(char *callFun, char *legacyKey, char *newKey) {

  if ( strlen(newKey) > 0 ) {
    sprintf(c1err,"Legacy key '%s' is no longer valid.", legacyKey);
    sprintf(c2err,"Use %s key instead .", newKey );
    errmsg(SEV_FATAL, 0, callFun, c1err, c2err );
  }
  else {
    sprintf(c1err,"Legacy key '%s' is no longer valid.", legacyKey);
    sprintf(c2err,"Remove key from file ."  );
    errmsg(SEV_FATAL, 0, callFun, c1err, c2err );
  }

} // end legacyKey_abort

// *********************************************************
void missingKey_ABORT(char *key, char *file, char *callFun) {
  // created Jan 2014
  char fnam[] = "missingKey_ABORT";
  sprintf(c1err,"%s could not find KEY='%s'", callFun, key);
  sprintf(c2err,"in file = %s", file);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );

} // end of missingKey_ABORT

void tabs_ABORT(int NTAB, char *fileName, char *callFun) {
  char fnam[] = "tabs_ABORT";
  if ( NTAB > 0 ) {
    sprintf(c1err,"%s found %d invalid tabs in file ", 
	    callFun, NTAB);
    sprintf(c2err,"%s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
} // end tabs_ABORT

// ************************************************
void errmsg(
       int isev            /* (I) severity flag */
      ,int iprompt         /* (I) 1=>prompt user to continue */
      ,char *fnam          /* (I) name of function calling ktrigerr */
      ,char *msg1          /* (I) message to print           */
      ,char *msg2          /* (I) 2nd msg to print ("" => no 2nd msg) */
          ) 
/*************************************
   Created Jan 30, 2004  by R.S. Kessler.
     [copied from ktrig.c]

   ABSTRACT:
      Print error message(s), and count no. of errors.
      Abort program if isevere == SEV_ABORT
        
*************************************/
{
   char c_severe[12];  /* char string for serverity */ 
   char cmsg[200];

        /* ------- begin function execution ------- */

   switch ( isev )
     {
        case SEV_INFO  : sprintf(c_severe,"INFO");     break;
        case SEV_WARN  : sprintf(c_severe,"WARNING");  break;
        case SEV_ERROR : sprintf(c_severe,"ERROR");    break;
        case SEV_FATAL : sprintf(c_severe,"FATAL");    break;
     }

 /* print SEVERITY, FUNCTION name, and 1st message. */
   sprintf(cmsg, "%s[%s]: \n\t %s", c_severe, fnam, msg1 );
   printf("\n %s \n", cmsg );

  /* write error message to file if global flag is set. */


 /* print 2nd message if non-null */

   if( strlen(msg2) > 0 ) { printf("\t %s", msg2 ); }

   printf("\n");

   if ( isev == SEV_FATAL ) { madend(1); }

   fflush(stdout);

}   /* end of function "errmsg" */




// ************************************************
void madend(int flag) {

   char cmsg[40] = { "ABORT program on Fatal Error." };

   printf("\n");
   printf("\n");
   printf("\n   `|```````|`    ");
   printf("\n   <| o\\ /o |>    ");
   printf("\n    | ' ; ' |     ");
   printf("\n    |  ___  |     %s ", cmsg);
   printf("\n    | |' '| |     ");
   printf("\n    | `---' |     ");
   printf("\n    \\_______/    ");
   printf("\n");

   printf("\n");   
   fflush(stdout);

   //   if ( flag == 1 ) { exit(1); }
   if ( flag == 1 ) { exit(EXIT_ERRCODE); }

}    //  end of "madend"  



// ************************************************
void happyend(void) {

   fflush(stdout);
   printf("\n Program stopping gracefully. Bye. \n");
   exit(0);

}


// ************************************************
void readint(FILE *fp, int nint, int *list)   
/*****
   Created 7-feb-95 by R.Kessler

   ABSTRACT: 
     Read "nint" integers from configuration file, 
     and put in array "list." 

  Apr 2019: abort on reading string
*****/
{
    char c_get[80];
    int  i, itmp, istat, scanStat ;
    char fnam[] = "readint";

    for ( i=0; i<nint; i++) {
      scanStat = fscanf(fp,"%s",c_get);         // read next string 
      itmp     = NOINT  ;
      istat    = sscanf ( c_get, "%d", &itmp );

      if ( itmp == NOINT ) {
	sprintf(c1err,"Could not read int from string='%s' ; ", c_get);
	sprintf(c2err,"reading item %d of %d items.", i+1, nint);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      list[i] = itmp; // xxxx atoi(c_get);       // convert to int   
     }

}  // end of function "readint" 


// ************************************************
void readlong(FILE *fp, int nint, long long *list)   
/*****
   Created Feb 2015 by R.Kessler

   ABSTRACT: 
     Read "nint" long long integers from configuration file, 
     and put in array "list." 
*****/
{
    int  i, scanStat ;

    for ( i=0; i<nint; i++)
      {  scanStat = fscanf(fp,"%lld", &list[i] );  }

}  // end of function "readlong" 




// ************************************************
void readfloat(FILE *fp, int nint, float *list)   
/************* 
   Created 22-JAN-96 by A. Roodman
   From READINT

   ABSTRACT: 
     Read "nint" floats from configuration file, 
     and put in array "list." 

   April 2019: abort reading string by mistake.

***************/
{
    char c_get[80];
    int  i, scanStat, fstat ;
    float ftmp;
    char fnam[] = "readfloat";

    for ( i=0; i<nint; i++) {
      scanStat = fscanf(fp,"%s",c_get);         // read next string 
      ftmp     = NOFLOAT  ;
      fstat    = sscanf ( c_get, "%f", &ftmp );
      if ( ftmp == NOFLOAT ) {
	sprintf(c1err,"Could not read float from string='%s' ; ", c_get);
	sprintf(c2err,"reading item %d of %d items.", i+1, nint);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      list[i] = ftmp;
    }

}  // end of function "readfloat"


// ************************************************
void readdouble(FILE *fp, int nint, double *list)   
/************* 
   Feb 2009 R.Kessler
   ABSTRACT: 
     Read "nint" doubles from configuration file, 
     and put in array "list." 

  Apr 2019: abort reading string by mistake.

***************/
{
    char c_get[80];
    int  i,scanStat, dstat ;
    double dtmp;
    char fnam[] = "readdouble";

    for ( i=0; i<nint; i++) {
      scanStat = fscanf(fp,"%s",c_get);         // read next string 
      dtmp     = NODOUBLE ;
      dstat    = sscanf ( c_get, "%le", &dtmp );

      if ( dtmp == NODOUBLE ) {
	sprintf(c1err,"Could not read double from string='%s' ; ", c_get);
	sprintf(c2err,"reading item %d of %d items.", i+1, nint);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

      list[i] = dtmp;
    }

}  // end of function "readdouble"



// ************************************************
void readchar(FILE *fp, char *clist)    
/*****
   Created 7-feb-95 by R.Kessler
   ABSTRACT: 
     Read next char. string from configuration file, 
     and put in the array "clist". 

******/
{
  int scanStat ;
  char fnam[] = "readchar";
  scanStat = fscanf(fp,"%s",clist);   // read next string 
  if ( clist[0] == '\t' ) {
    sprintf(c1err,"Invalid tab read in string '%s'", clist);
    sprintf(c2err,"Check file being read.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

}  // end of function "readchar"



// ************************************************
void print_banner (const char *banner ) {

  printf("\n ********************************************"
	 "********************** \n");
  printf("   %s  \n" , banner );
  fflush(stdout);

}
void debugexit(char *string) {
  printf("\n xxx DEBUG EXIT: %s \n", string);
  fflush(stdout);
  exit(1);
}


// ************************************************
double zhelio_zcmb_translator (double z_input, double RA, double DEC, 
			       char *coordSys, int OPT ) {

  /**********************

   Created Dec 2013 by R.Kessler
   General redshift-translator function to between helio and cmb frame.
   [replaces Z2CMB that translated only in 1 direction]

   OPT > 0  -> z_input = z_helio, return z_out = zcmb
   OPT < 0  -> z_input = z_cmb,   return z_out = zhelio
   RA,DEC   = sky coordinates
   coordSys = coordinate system; e.g., 'eq' or 'gal' or 'J2000'


   l = longitude = RA (deg)
   b = lattitue  = DEC  (deg)

   Use exact formuala,
   
    1 + z_cmb = ( 1 + z_helio ) / ( 1 - V0 . nhat/c )

   where V0 is the CMB velocity vector, 
   nhat is the unit vector between us and the SN,
   and c = 3E5 km/s.

   Note that the NED-calculator used in JRK07 is 
   an approximation, z_cmb = z_helio + V0.nhat/c,
   that is OK when z_helio << 1.

 ****************/

  double 
     ra_gal, dec_gal
    ,ss, ccc, c1, c2, c3, vdotn, z_out
    ;

  char fnam[] = "zhelio_zcmb_translator" ;

  // --------------- BEGIN ------------

  // on negative redshift, just return input redshift with
  // no calculation. Allows flags such as -9 to be unperturbed.
  if ( z_input < 1.0E-10 ) { return z_input ; }

  if ( strcmp(coordSys,"eq"   ) == 0 || 
       strcmp(coordSys,"J2000") == 0 ) {

    // input and output in degrees
    slaEqgal( RA, DEC, &ra_gal, &dec_gal ) ;
  }
  else if ( strcmp(coordSys,"gal") == 0 ) {
    ra_gal  = RA ;
    dec_gal = DEC ;
  }
  else {
    sprintf(c1err,"Invalid coordSys = '%s' ", coordSys );
    sprintf(c2err,"OPT=%d z_in=%f RA=%f DEC=%f", OPT, z_input, RA, DEC);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // get projection

  ss  = sin(RADIAN*dec_gal) * sin(RADIAN*CMBapex_b);
  c1  = cos(RADIAN*dec_gal) ;
  c2  = cos(RADIAN*CMBapex_b) ;
  c3  = cos(RADIAN*(ra_gal-CMBapex_l));
  ccc = c1 * c2 * c3 ;
  vdotn = CMBapex_v * ( ss + ccc ) / LIGHT_km ;

  z_out = -9.0 ;

  if ( OPT > 0 ) {
    z_out  = ( 1. + z_input ) / ( 1. - vdotn ) - 1. ; 
  }
  else if ( OPT < 0 )  {
    z_out  = ( 1. + z_input) * ( 1. - vdotn ) - 1.0 ;  
  }
  else if ( OPT == 0 ) {
    sprintf(c1err,"Invalid OPT=0" );
    sprintf(c2err,"z_input=%f  RA=%f  DEC=%f", z_input,RA,DEC);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
   
  return(z_out) ;


} // end of zhelio_zcmb_translator 


double zhelio_zcmb_translator__ (double *z_input, double *RA, double *DEC,
                                 char *coordSys, int *OPT ) {
  return zhelio_zcmb_translator(*z_input, *RA, *DEC, coordSys, *OPT) ;
}

// **********************************
double Z2CMB(
	     double z_helio   // (I) heliocentric redshift
	     ,double RA       // (I) RA, deg
	     ,double DEC      // (I) DEC, deg
	     ,char *coordSys  // (I) 'eq' or 'gal' or 'J2000'
	     ) {

  /**********

 
   Feb 28, 2011 (R.Kessler)
   General function to translate z_helio into z_cmb (return arg).
   Allows input coords in either J2000 or equatorial.

   l = longitude = RA (deg)
   b = lattitue  = DEC (deg)

   Use exact formuala,
   
    1 + z_cmb = ( 1 + z_helio ) / ( 1 - V0 . nhat/c )

   where V0 is the CMB velocity vector, 
   nhat is the unit vector between us and the SN,
   and c = 3E5 km/s.

   Note that the NED-calculator used in JRK07 is 
   an approximation, z_cmb = z_helio + V0.nhat/c,
   that is OK when z_helio << 1.

   Installed in sntools.c on Oct 29 2013 by RK : NOT TESTED YET !!!
 
  Nov 18 2013: tested ! Note that David's slaEqgal has all args in 
               degrees instead of radians.

  Dec 6 2013: for negative redshift, do nothing and return z_helio

  *************/

  double 
     ra_gal, dec_gal
    ,ss, ccc, c1, c2, c3, vdotn, z_cmb
    ;

  // define location and velocity of CMB dipole
  double  l_CMBapex  = 264.14 ;   // deg (RA galactic coords !!!)
  double  b_CMBapex  = 48.26 ;    // deg (DEC)
  double  v_CMBapex  = 371.0 ;    // km/sec

  char fnam[] = "Z2CMB";

  // --------------- BEGIN ------------

  if ( z_helio < 0 ) { return z_helio ; }

  if ( strcmp(coordSys,"eq"   ) == 0 || 
       strcmp(coordSys,"J2000") == 0 ) {

    // input and output in degrees
    slaEqgal( RA, DEC, &ra_gal, &dec_gal ) ;
  }
  else if ( strcmp(coordSys,"gal") == 0 ) {
    ra_gal  = RA ;
    dec_gal = DEC;
  }
  else {
    sprintf(c1err,"coordSys = '%s' is invalid", coordSys );
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
  }

  // get projection

  ss  = sin(RADIAN*dec_gal) * sin(RADIAN*b_CMBapex);

  c1  = cos(RADIAN*dec_gal) ;
  c2  = cos(RADIAN*b_CMBapex) ;
  c3  = cos(RADIAN*(ra_gal-l_CMBapex));
  ccc = c1 * c2 * c3 ;

  vdotn = v_CMBapex * ( ss + ccc ) / LIGHT_km ;

  z_cmb    = ( 1. + z_helio ) / ( 1. - vdotn ) - 1. ; // exact

  // z_cmb = z_helio + vdotn;  // approx xxxxx
   
  return z_cmb ;

}   //  end of Z2CMB

double z2cmb_(double *z_helio, double *RA, double *DEC, char *coordSys) {
  return Z2CMB(*z_helio, *RA ,*DEC, coordSys);
}


// Altered from the fortran SLALIB by David Cinabro, June 2006.
// Translates equatorial coordinats (RA,DEC) to galactic
// longitude and latitude.  All in degrees and double precision.
// All the subroutines needed are included below.
// Usage: 
//    slaEqgal ( double RA, double DEC, double *GalLat, double *GalLong );

void slaEqgal ( double dr, double dd, double *dl, double *db )
/*
**  - - - - - - - - -
**   s l a E q g a l
**  - - - - - - - - -
**
**  Transformation from J2000.0 equatorial coordinates to
**  IAU 1958 Galactic coordinates.
**
**  (double precision)
**
**  Given:
**     dr,dd       double       J2000.0 RA,Dec
**
**  Returned:
**     *dl,*db     double       Galactic longitude and latitude l2,b2
**
**  (all arguments were radians, but translation from and to degrees done below)
**
**  Called:
**     slaDcs2c, slaDmxv, slaDcc2s, slaDranrm, slaDrange
**
**  Note:
**     The equatorial coordinates are J2000.0.  Use the routine
**     slaEg50 if conversion from B1950.0 'FK4' coordinates is
**     required.
**
**  Reference:
**     Blaauw et al, Mon.Not.R.astron.Soc.,121,123 (1960)
**
**  Last revision:   21 September 1998
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double v1[3], v2[3];
   double drr, ddr;
   double DPI = 3.1415926535897932384626433832795028841971693993751;

/*
**  l2,b2 system of Galactic coordinates
**
**  p = 192.25       RA of Galactic north pole (mean B1950.0)
**  q =  62.6        inclination of Galactic to mean B1950.0 equator
**  r =  33          longitude of ascending node
**
**  p,q,r are degrees
**
**  Equatorial to Galactic rotation matrix (J2000.0), obtained by
**  applying the standard FK4 to FK5 transformation, for zero proper
**  motion in FK5, to the columns of the B1950 equatorial to
**  Galactic rotation matrix:
*/
   static double rmat[3][3];

   rmat[0][0] = -0.054875539726;
   rmat[0][1] = -0.873437108010;
   rmat[0][2] = -0.483834985808;
   rmat[1][0] =  0.494109453312;
   rmat[1][1] = -0.444829589425;
   rmat[1][2] =  0.746982251810;
   rmat[2][0] = -0.867666135858;
   rmat[2][1] = -0.198076386122;
   rmat[2][2] =  0.455983795705;

   // Translate to radians
   drr = dr*DPI/180.0;
   ddr = dd*DPI/180.0;

/* Spherical to Cartesian */
   slaDcs2c ( drr, ddr, v1 );

/* Equatorial to Galactic */
   slaDmxv ( rmat, v1, v2 );

/* Cartesian to spherical */
   slaDcc2s ( v2, dl, db );

/* Express in conventional ranges */
   *dl = slaDranrm ( *dl );
   *db = slaDrange ( *db );
   // Translate back to degrees
   *dl = *dl*180.0/DPI;
   *db = *db*180.0/DPI;
}

void slaDcs2c ( double a, double b, double v[3] )
/*
**  - - - - - - - - -
**   s l a D c s 2 c
**  - - - - - - - - -
**
**  Spherical coordinates to direction cosines (double precision)
**
**  Given:
**     a,b       double      spherical coordinates in radians
**                           (RA,Dec), (long,lat) etc
**
**  Returned:
**     v         double[3]   x,y,z unit vector
**
**  The spherical coordinates are longitude (+ve anticlockwise looking
**  from the +ve latitude pole) and latitude.  The Cartesian coordinates
**  are right handed, with the x axis at zero longitude and latitude,
**  and the z axis at the +ve latitude pole.
**
**  Last revision:   22 July 2004
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double cosb;

   cosb = cos ( b );
   v[0] = cos ( a ) * cosb;
   v[1] = sin ( a ) * cosb;
   v[2] = sin ( b );
}

void slaDmxv ( double dm[3][3], double va[3], double vb[3] )
/*
**  - - - - - - - -
**   s l a D m x v
**  - - - - - - - -
**
**  Performs the 3-d forward unitary transformation:
**     vector vb = matrix dm * vector va
**
**  (double precision)
**
**  Given:
**     dm       double[3][3]    matrix
**     va       double[3]       vector
**
**  Returned:
**     vb       double[3]       result vector
**
**  Note:  va and vb may be the same array.
**
**  Last revision:   22 July 2004
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int i, j;
   double w, vw[3];


/* Matrix dm * vector va -> vector vw. */
   for ( j = 0; j < 3; j++ ) {
      w = 0.0;
      for ( i = 0; i < 3; i++ ) {
         w += dm[j][i] * va[i];
      }
      vw[j] = w;
   }

/* Vector vw -> vector vb. */
   for ( j = 0; j < 3; j++ ) {
      vb[j] = vw[j];
   }
}

void slaDcc2s ( double v[3], double *a, double *b )
/*
**  - - - - - - - - -
**   s l a D c c 2 s
**  - - - - - - - - -
**
**  Cartesian to spherical coordinates.
**
**  (double precision)
**
**  Given:
**     v       double[3]   x,y,z vector
**
**  Returned:
**     *a,*b   double      spherical coordinates in radians
**
**  The spherical coordinates are longitude (+ve anticlockwise looking
**  from the +ve latitude pole) and latitude.  The Cartesian coordinates
**  are right handed, with the x axis at zero longitude and latitude,
**  and the z axis at the +ve latitude pole.
**
**  If v is null, zero a and b are returned.  At either pole, zero a is
**  returned.
**
**  Last revision:   22 July 2004
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double x, y, z, r;

   x = v[0];
   y = v[1];
   z = v[2];
   r = sqrt ( x * x + y * y );

   *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
   *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
}

double slaDranrm ( double angle )
/*
**  - - - - - - - - - -
**   s l a D r a n r m
**  - - - - - - - - - -
**
**  Normalize angle into range 0-2 pi.
**
**  (double precision)
**
**  Given:
**     angle     double      the angle in radians
**
**  The result is angle expressed in the range 0-2 pi (double).
**
**  Defined in slamac.h:  D2PI, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double w;
   double D2PI = 6.2831853071795864769252867665590057683943387987502;
/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                                        :(A)+(B)*floor(-(A)/(B))):(A))

   w = dmod ( angle, D2PI );
   return ( w >= 0.0 ) ? w : w + D2PI;
}

double slaDrange ( double angle )
/*
**  - - - - - - - - - -
**   s l a D r a n g e
**  - - - - - - - - - -
**
**  Normalize angle into range +/- pi.
**
**  (double precision)
**
**  Given:
**     angle     double      the angle in radians
**
**  The result is angle expressed in the range +/- pi.
**
**  Defined in slamac.h:  DPI, D2PI, dmod
**
**  Last revision:   22 July 2004
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
  double w;
  double DPI = 3.1415926535897932384626433832795028841971693993751;
  double D2PI = 6.2831853071795864769252867665590057683943387987502;
/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                                        :(A)+(B)*floor(-(A)/(B))):(A))
/* dsign(A,B) - magnitude of A with sign of B (double) */
#define dsign(A,B) ((B)<0.0?-(A):(A))

  w = dmod ( angle, D2PI );
  return ( fabs ( w ) < DPI ) ? w : w - dsign ( D2PI, angle );
}

