// Created Oct 2014 by R.Kessler
//
// GET_MULTISEASON function to determine seasons and to
// return reduced chi2 per season, where the chi2 is based
// on a model with zero flux. The goal is that large reduced
// chi2 corresponds to transient activity, while chi2red ~ 1
// corresponds to no activity. Initial motivation is to find
// AGN, but this function is not linked to any particular
// physical model and can thus be used for any mult-season
// sources.
//


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "multiseason.h" 

// =============================================================

void INIT_MULTISEASON(float *parList) {

  // input parList is
  //   OPTMASK:  option mask (currently ignored)
  //   TGAP   :  gap in days defines season
  //   NREJECT_OUTLIER  : number of outliers to reject per season

  char  fnam[] = "INIT_MULTISEASON" ;

  // ------------ BEGIN ---------

  INPUTS.OPTMASK          = (int)parList[0] ;
  INPUTS.TGAP             = (float)parList[1] ;
  INPUTS.NREJECT_OUTLIER  = (int)parList[2] ;

  if ( INPUTS.OPTMASK == 0 ) { return ; }
  if ( INPUTS.TGAP > 1.0E8 ) { return ; }

  sprintf(BANNER,"%s: init multi-season analysis.", fnam );
  print_banner(BANNER);

  if ( INPUTS.OPTMASK == 1 ) 
    { printf("\t Flux-activity reference is zero. \n"); }
  else if ( INPUTS.OPTMASK == 2 ) 
    { printf("\t Flux-activity reference is AVG-Flux per season \n"); }
  else {
    sprintf(c1err,"Invalid OPTMASK=%d", INPUTS.OPTMASK);
    sprintf(c2err,"Valid MULTISEASON_OPTMASK is 1 or 2");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
    
  fflush(stdout);

  
} // end of INIT_MULTISEASON

void init_multiseason__(float *parList) { INIT_MULTISEASON(parList); }


// ========================================================
void GET_MULTISEASON(char *CCID, int NOBS, double *MJD_LIST, 
		     double *FLUX_LIST, double *FLUXERR_LIST,
		     int *REJECT, int *NSEASON, double *CHI2RED,      // (O)
		     double *MJDMIN, double *MJDMAX, double *AVGFLUX) // (O)
{     

  // Inputs:
  //   CCID : character SNID (for error messages only)
  //   NOBS : number of observations
  //   MJD_LIST  : list of NOBS MJDs
  //   FLUX_LIST : list of NOBS fluxes
  //   FLUXERR_LIST: list of NOBS flux-errors 
  //          (ignore obs with negative FLUXERR)
  //
  // Outputs:
  //   REJECT:   list of rejected epochs
  //   NSEASON:  number of seasons found (1-MXSEASON-1)
  //   CHI2RED: reduced chi2 per season
  //   MJDMIN:  min MJD for each season
  //   MJDMAX:  max MJD for each season
  //   AVGFLUX: avg flux per season (chi2 is w.r.t. this flux)
  //
  // May 1 2019: 
  //   + convert float to double
  //   + add new option to compare flux to season-average 
  //     instead of zero.

  int   OPT_FLUXREF_ZERO = ( INPUTS.OPTMASK == 1 ) ;
  int   OPT_FLUXREF_AVG  = ( INPUTS.OPTMASK == 2 ) ;

  int   obs, Nobs[MXSEASON], Ndof[MXSEASON] ;
  int   *INDX_SORT, *INDX_SORT2,  *ISEASON;
  int   Nseason, iSeason, i, isort, Ntmp ;
  double MJD_LAST, XN, MJD, SNR, CHI2, WGT, SUMCHI2[MXSEASON] ;
  double *CHI2_LIST, *chi2_list, FLUX, ERR, FDIF, CHI2TMP ;
  double SUMFLUX[MXSEASON], SUMWGT[MXSEASON] ;
  int   *OBS_LIST,  *obs_list ;
  int   LDMP = 0;
  char fnam[] = "GET_MULTISEASON" ;

  // --------------------- BEGIN --------------------

  if ( INPUTS.TGAP > 1.0E8 ) { return ; }

  // init output args.
  *NSEASON = 0 ;

  for(i=0; i < MXSEASON; i++ ) {
    CHI2RED[i] = -9.0 ;
    MJDMIN[i]  = -9.0 ;
    MJDMAX[i]  = -9.0 ;    
    
    SUMWGT[i]  = 0.0 ;
    SUMFLUX[i] = 0.0 ;
    SUMCHI2[i] = 0.0 ;   // local  var
    Ndof[i]    = 0 ;     // local var
    Nobs[i]    = 0 ;     // number of obs per season
  }

  // ---------------------------------------------
  // malloc local arrays
  int MEMI = NOBS * sizeof(int);
  int MEMD = NOBS * sizeof(double);

  INDX_SORT  = (int   *) malloc( MEMI ) ;
  INDX_SORT2 = (int   *) malloc( MEMI ) ;
  ISEASON    = (int   *) malloc( MEMI ) ;
  CHI2_LIST  = (double*) malloc( MEMD );
  chi2_list  = (double*) malloc( MEMD );

  OBS_LIST   = (int*) malloc( MEMI ) ;
  obs_list   = (int*) malloc( MEMI ) ;

  // -----------------------------------------
  // first find season(s) with sorted MJD

  sortDouble(NOBS, MJD_LIST, +1, INDX_SORT);

  MJD_LAST = -1.0E8 ;  iSeason  = -1 ;

  for(obs=0; obs < NOBS; obs++ ) {

    isort = INDX_SORT[obs];
    CHI2_LIST[isort]   = 0.0 ;
    REJECT[isort]      = 1 ; // init all to rejected
    OBS_LIST[isort]    = -9;

    if ( FLUXERR_LIST[isort] <= 0.0 ) { continue ; }

    SNR   = FLUX_LIST[isort] / FLUXERR_LIST[isort] ;
    CHI2  = SNR * SNR ;
    CHI2_LIST[isort] = CHI2; 
    OBS_LIST[isort]  = obs ;

    MJD = MJD_LIST[isort] ;
    if ( MJD > MJD_LAST + INPUTS.TGAP ) {
      iSeason++ ;      
      MJDMIN[iSeason] = MJD ;
    }
    MJD_LAST = MJD ;
    ISEASON[isort]  = iSeason;
    MJDMAX[iSeason] = MJD ;
    
    // compute weighted avg
    ERR   = FLUXERR_LIST[isort] ;
    WGT   = 1.0/(ERR*ERR);
    FLUX  = FLUX_LIST[isort];
    SUMFLUX[iSeason] += (WGT*FLUX);
    SUMWGT[iSeason]  += WGT ;
    Nobs[iSeason]++ ;

  } // end of obs

  Nseason = iSeason + 1 ;
  *NSEASON = Nseason ;  // load output arg.

  if ( LDMP ) { 
    printf(" xxx -------------------------------------------- \n" ); 
    CHI2TMP = 0.0 ;
  }

  // - - - - - - - - - -
  // compute avg flux per season, and reset Nobs
  for(i=0; i < Nseason; i++ ) {
    if ( SUMWGT[i] > 0 ) { AVGFLUX[i] = SUMFLUX[i]/SUMWGT[i]; }

    if ( LDMP ) 
      { printf(" xxx AVGFLUX[%d] = %.1f \n",i,AVGFLUX[i]); fflush(stdout);}
    
    Nobs[i] = 0;
  }

  // ------------------------------------------------------------
  // compute reduced chi2 for each season, with outlier rejection

  
  for(iSeason=0; iSeason < Nseason; iSeason++ ) {

    /* xxxxxxxx mark delete May 1 2019 xxxxxxxxxxxxx
    FSUM = WSUM = 0.0 ;
    for(obs=0; obs < NOBS; obs++ ) {
      isort = INDX_SORT[obs] ;
      if ( iSeason == ISEASON[isort] ) {
	W     = 1.0/(FLUXERR_LIST[isort]*FLUXERR_LIST[isort]);
	FSUM += W * FLUX_LIST[isort] ;	WSUM += W ;
      }
    }
    AVGFLUX[iSeason] = FSUM / WSUM ;
    xxxxxxxxxxxx */


    // build chi2_list for this season
    for(obs=0; obs < NOBS; obs++ ) {
      isort = INDX_SORT[obs] ;
      if ( iSeason == ISEASON[isort] ) {
	Nobs[iSeason]++ ;  
	Ntmp = Nobs[iSeason] ;

	MJD  = MJD_LIST[isort] ;
	FLUX = FLUX_LIST[isort] ;
	ERR  = FLUXERR_LIST[isort] ;

	if ( OPT_FLUXREF_ZERO ) 
	  { FDIF = FLUX - 0.0 ; }
	else if ( OPT_FLUXREF_AVG ) 
	  { FDIF = FLUX - AVGFLUX[iSeason] ; }

	CHI2              = (FDIF*FDIF)/(ERR*ERR);
	chi2_list[Ntmp-1] = CHI2 ;
	obs_list[Ntmp-1]  = OBS_LIST[isort] ;

	if ( LDMP && iSeason == 1 ) {
	  printf(" xxx MJD=%9.3f  F=%8.1f +- %5.1f  CHI2=%5.1f \n",
		 MJD, FLUX, ERR, CHI2); fflush(stdout);
	  CHI2TMP += CHI2;
	}

      }
    } // obs

    // - - - - - - - - - - - 
    // sort chi2_list (smallest to largest)
    Ntmp = Nobs[iSeason] ;
    sortDouble(Ntmp, chi2_list, +1, INDX_SORT2);

    // Now sum chi2 excluding last NREJECT_OUTLIERS
    Ntmp = Nobs[iSeason] - INPUTS.NREJECT_OUTLIER ;
    for(obs=0; obs < Ntmp; obs++ ) {
      isort = INDX_SORT2[obs] ;
      CHI2  = chi2_list[isort];
      SUMCHI2[iSeason] += CHI2 ;
      Ndof[iSeason]++ ;

      REJECT[obs_list[isort]]  = 0 ; // turn off reject flag for used epochs
    }

    XN = (double)Ndof[iSeason] ;
    if ( XN > 0.001 ) 
      { CHI2RED[iSeason] = SUMCHI2[iSeason]/ (double)Ndof[iSeason] ; }

  } // iSeason



  if ( LDMP ) {
    printf(" xxx CHI2TMP/Nobs(ALL) = %.1f/%d = %.1f \n",
	   CHI2TMP, Nobs[1], CHI2TMP/(double)Nobs[1] );
    fflush(stdout);
  }


  // ----------------------------------
  if ( LDMP ) {
    //    printf(" xxx ----------------------------------------- \n");
    printf(" xxx %s : SNID = %s  with  NOBS=%d \n", fnam, CCID, NOBS);
    for (i=0; i < Nseason; i++ ) {
      printf(" xxx Season %2d : MJD = %.3f to %.3f  "
	     "chi2red = %.1f/%d = %6.2f\n",
	     i, MJDMIN[i], MJDMAX[i], SUMCHI2[i], Ndof[i], CHI2RED[i] );
    }
    fflush(stdout);
  }

  // ------- cleanup ------
  free(INDX_SORT);
  free(INDX_SORT2);
  free(ISEASON);
  free(CHI2_LIST);
  free(chi2_list);
  free(OBS_LIST);
  free(obs_list);

} // end of GET_MULTISEASON


void get_multiseason__(char *CCID, int *NOBS, double *MJD,
		       double *FLUX, double *FLUXERR,
                       int *REJECT, int *NSEASON, double *CHI2RED,
                       double *MJDMIN, double *MJDMAX, double *AVGFLUX ) {

  GET_MULTISEASON(CCID, *NOBS, MJD, FLUX, FLUXERR,     // (I)
		  REJECT, NSEASON, CHI2RED, MJDMIN, MJDMAX, AVGFLUX); // (O)

}  // end of get_multiseason
