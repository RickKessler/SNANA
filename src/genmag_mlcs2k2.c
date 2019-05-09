/**********************************************************************/
/* Generate rest-frame mags for an SN using the MLCS template         */
/* description                                                        */
/*                                                                    */
/* G. Miknaitis, May 31, 2006                                         */
/* R. Kessler                                                         */
/* J. Bernstein (JPB) April 2009                                      */
/*                                                                    */
/* init_genmag_mlcs2k2():  load template data into array              */
/*                                                                    */
/*   inputs:                                                          */ 
/*   outputs:                                                         */  
/*      *template_mags   = 2D array of mag vs. time for 9 bands       */
/*      *template_dates  = 1D array of template dates,                */
/*                                  relative to Bmax                  */
/*                                                                    */
/* genmag_mlcs2k2():    return photometry in a single passband,       */
/*                      for a set of input dates and MLCS delta       */
/*                                                                    */
/*   inputs:                                                          */
/*      ifilt            = filter number, 012345679<-> UBVRIYJHK      */
/*      rest_dates        = array of dates for which we want template */
/*                         photometry, relative to Bmax               */
/*      nobs             = number of obs dates                        */
/*      delta            = MLCS delta parameter                       */
/*   outputs:                                                         */
/*      mags             = template mags at requested dates and filt  */
/*                                                                    */
/*                                                                    */
/* NOTES:                                                             */
/*  - genmag_mlcs should abort if given a crazy delta, but doesn't.   */
/*                                                                    */
/*  - The format for the template data file is:		              */
/*    							              */
/*    day since B-band maximum				              */
/*    template						              */
/*    offset							      */
/*    pvalue							      */
/*    qvalue							      */
/*    							              */
/*    from which template mags for a given delta can be computed:     */
/*    							              */
/*       mag = template + offset + pvalue*delta + qvalue*delta*delta  */
/*                                                                    */
/*  - Code assumes that each template data file has the same          */
/*    range and number of days, and that the lines are time-ordered   */
/*                                                                    */
/*                                                                    */

/************ HISTORY *********************

 Feb 22 2017: in init_genmag_mlcs2k2, allow input VERSION to 
              include arbitrary path.

 Dec 29 2017: use open_TEXTgz() to allow reading gzipped text files.
      
*******************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sntools.h"

#include "genmag_mlcs2k2.h"

// -------------------------------------------------------------
//     define mangled routines with underscore for fortran
// -------------------------------------------------------------


int init_genmag_mlcs2k2__(char *version, char *covFile, float *scale_covar,
		       char *filtlist ) {
  int istat;
  istat = init_genmag_mlcs2k2 ( version, covFile, *scale_covar, filtlist) ;
  return istat;
}


 
int genmag_mlcs2k2__(
		int *ifilt           // (I) filter index 012345678 => UBVRIYJHK
		,double *delta       // (I) luminosity parameter
		,int *nobs           // (I) number of epochs
		,double *rest_dates  //  (I) list of Trest(days)
		,double *rest_magval   //  (O) rest-frame mags values
		,double *rest_magerr   //  (O) rest-frame mag errors
		)
{
  int istat;
  istat = genmag_mlcs2k2(*ifilt,*delta,*nobs, rest_dates,   // input
		      rest_magval, rest_magerr);         // output

  return istat;

}

int gencovar_mlcs2k2__ (
		  int *matsize         // (I) row-len (or col-len)
		  ,int *ifilt          // (I) list of 'matsize' filter indices
		  ,double *rest_epoch  // (I) list of 'matsize' rest days
		  ,double *covar       // (O) covariance matrix
		  ) {
  int istat ;
  istat = gencovar_mlcs2k2 ( *matsize, ifilt, rest_epoch, covar ) ;
  return istat;
}


void get_lamrange_mlcs2k2__(double *lammin, double *lammax) {
  get_LAMRANGE_mlcs2k2(lammin,lammax);
}

// ******************************************
int init_genmag_mlcs2k2(
		     char *VERSION   // (I) version or PATH/VERSION
		     ,char *covFile  // (I) name of cov file
		     ,float scale_covar  // (I) scale cov (i.e, sqerr)
		     ,char *filtlist     // (O) char-string of model filters
		     ) {
  
  char tempFilename[200];
  char fnam[] = "init_genmag_mlcs2k2" ;
  char BANNER[80], cfilt[2], version[60];

  int ifilt, nday, jdum, GZIPFLAG ;
  FILE *file_ptr ;
  double tmp_day, tmp_offset, tmp_pvalue, tmp_qvalue;
  double T0, tmpmin, tmpmax ;

  // ---------- BEGIN -------

  if ( scale_covar <= 0.0 ) scale_covar = 1.0;  // temp until bug is fixed



  sprintf(BANNER, "%s : Initialize mlcs2k2 Light Curve Generator", fnam );
  print_banner ( BANNER );

  extract_MODELNAME(VERSION,
		    PATHMODEL_MLCS2k2, version); // returned

  if ( getenv(PRIVATE_MODELPATH_NAME) != NULL ) {
    // path set by ENV
    sprintf( PATHMODEL_MLCS2k2, "%s/%s", 
	     getenv(PRIVATE_MODELPATH_NAME), version );    
  }
  else if ( strlen(PATHMODEL_MLCS2k2) > 0 ) {
    // do nothing for user-defined path/model
  }
  else {
    // default location under $SNDATA_ROOT
    sprintf(PATHMODEL_MLCS2k2, "%s/models/mlcs2k2/%s", 
	    getenv("SNDATA_ROOT"), version );
  }  


  printf("\t Read mlcs2k2 model params from PATHMODEL_MLCS= \n\t   %s\n", 
	 PATHMODEL_MLCS2k2);


  set_LAMRANGE_mlcs2k2();

  filtlist[0] = 0 ;

  sprintf(FILTSTRING_MLCS2k2, "UBVRIYJHK" );

  // Loop over filters and load up template data array 
  NFILT_MLCS2k2 = 0 ;

  for (ifilt=0; ifilt < MAXFILT_MLCS2k2 ; ifilt++){

    nday=0;  /* Initialize day counter */

    /* Get filename for template data for this filt */
    sprintf(tempFilename,"%s/vectors_%c.dat",
	    PATHMODEL_MLCS2k2, FILTSTRING_MLCS2k2[ifilt] ) ;

    // Open file for reading 
    // Abort on error only if optical filter is missing.
    // IR filters are allowed to be missing.

    file_ptr = open_TEXTgz(tempFilename, "rt", &GZIPFLAG ) ;
    // xxx mark delete if (( file_ptr = fopen(tempFilename, "rt")) == NULL ){
      
    if ( file_ptr == NULL ) {
      if ( ifilt < NFILT_MLCS2k2_REQUIRED ) {
	printf("\n Missing required mlcs2k2 template: \n");
	printf("    %s\n", tempFilename );
	printf(" ***** ABORT ***** \n");
	exit(1);
	//	error++ ;
      }

    } 
    else {

      while (fscanf(file_ptr, "%lf %lf %lf %lf", 
		    &tmp_day, &tmp_offset,
		    &tmp_pvalue, &tmp_qvalue)>0){

	TEMPLATE_DAYS_MLCS2k2[nday]               = tmp_day;
        TEMPLATE_MAGS_MLCS2k2[nday][ifilt].MAGOFF = tmp_offset;
        TEMPLATE_MAGS_MLCS2k2[nday][ifilt].D1     = tmp_pvalue ;
        TEMPLATE_MAGS_MLCS2k2[nday][ifilt].D2     = tmp_qvalue ;
	nday++;
	NDAY_MLCS2k2 = nday ;

      }  // end of fscanf


      NFILT_MLCS2k2++ ;
      sprintf(cfilt,"%c", FILTSTRING_MLCS2k2[ifilt] );
      sprintf(filtlist,"%s%s", filtlist, cfilt );

      if ( strcmp(cfilt,"Y") == 0 ) LAMRANGE_MLCS2k2[1] = MAXLAM_MLCS2k2_Y ;
      if ( strcmp(cfilt,"J") == 0 ) LAMRANGE_MLCS2k2[1] = MAXLAM_MLCS2k2_J ;
      if ( strcmp(cfilt,"H") == 0 ) LAMRANGE_MLCS2k2[1] = MAXLAM_MLCS2k2_H ;
      if ( strcmp(cfilt,"K") == 0 ) LAMRANGE_MLCS2k2[1] = MAXLAM_MLCS2k2_K ;

      // make sure that each filter has exactly the same epoch range

      if ( NFILT_MLCS2k2 == 1 ) {
	// store min/max Trest for template
	TMINDAY_MLCS2k2 = TEMPLATE_DAYS_MLCS2k2[0] ;  
	TMAXDAY_MLCS2k2 = TEMPLATE_DAYS_MLCS2k2[NDAY_MLCS2k2-1] ; 
      }
      else {
	tmpmin = TEMPLATE_DAYS_MLCS2k2[0] ; 
	tmpmax = TEMPLATE_DAYS_MLCS2k2[NDAY_MLCS2k2-1] ; 
	if ( tmpmin != TMINDAY_MLCS2k2 || tmpmax != TMAXDAY_MLCS2k2 ) {
	  printf("\n\n\n FATAL ERROR in %s (filtlist=%s) : \n", 
		 fnam, filtlist );
	  printf(" Found different epoch-range in  \n");
	  printf("    %s .\n", tempFilename );
	  printf(" All epoch ranges must be the same for mlcs2k2 vectors.\n");
	  printf(" ***** ABORT ***** \n");
	  exit(1);
	}
      }

      fclose(file_ptr);

    }  // end of file-open if-block

  }  // end of ifilt loop

  /* Keep track of number of points */


  printf("\n");

  printf("\t RESTLAMBDA-range: %d to %d A \n",
	 (int)LAMRANGE_MLCS2k2[0], (int)LAMRANGE_MLCS2k2[1] );

  printf("\t TREST-range: %5.1f  to  %5.1f  days (NDAY=%d) \n", 
	 TMINDAY_MLCS2k2, TMAXDAY_MLCS2k2, NDAY_MLCS2k2 );

  // Mar 2010: special flag to fix buggy vectors for SNchallenge
  if ( strcmp(version,"mlcs2k2.SNchallenge") == 0 ) {
    QWGT_FLAG = 1;
    MOFF_MLCS2k2 = +0.16; // correct rest mags to H0=70
    printf("\t Suppress Q for Trest > %3.0f and DELTA > %3.1f (MOFF=%4.2f mag) \n",
	   QWGT_TRESTMIN, QWGT_DELTAMIN, MOFF_MLCS2k2 );
  }
  else 
    { QWGT_FLAG = 0; MOFF_MLCS2k2 = 0.0 ; }


  // read covariance matrix.
  rd_mlcs2k2_cov(covFile,scale_covar);

  // screen dump info about covariance matrix.
  dmp_mlcs2k2_cov();


  // get IDAYPEAK_MLCS2k2 = day-index at peakmag
  T0 = 0.0 ;
  JDATES ( T0, &IDAYPEAK_MLCS2k2, &jdum ) ; 


  // misc. inits
  USE_PEAKERR_ONLY = 0 ; 

  return +1 ;

}	   // end of init_genmag_mlcs2k2



// ***************************************************
int genmag_mlcs2k2(
		int ifilt        // (I) filter index 012345678 => UBVRIYJHK
		,double delta          // (I) luminosity parameter
		,int nobs              // (I) number of epochs
		,double *rest_dates    // (I) list of Trest(days)
		,double *rest_magval   // (O) rest-frame mag values
		,double *rest_magerr   // (O) rest-frame mag errors
		){

  double 
     w1, w2
    ,tempmag1, tempmag2
    ,temperr1, temperr2
    ,Trest
    ,sqdelta
    ,M1, P1, Q1
    ,M2, P2, Q2
    ;

  int i, j1, j2, idum;

  // ---------- BEGIN function ------------

  sqdelta = delta * delta ;

  for(i=0; i < nobs; i++) {
    
    Trest = rest_dates[i];

    QWGT_MLCS2k2 = get_QWGT_mlcs2k2(Trest, delta);

    // check of Trest is within valid template range.

    if ( Trest > TMAXDAY_MLCS2k2 ) {  
      j2 = NDAY_MLCS2k2-1 ;

      get_MPQ_mlcs2k2(j2,ifilt, &M2, &P2, &Q2 );
      rest_magval[i] = M2 + P2*delta + Q2*sqdelta 
	+ 0.03 * (Trest - TMAXDAY_MLCS2k2); // add .03 mag per day past template

      rest_magerr[i] = 1.0 ; 
      continue ; 
    }

    // If Trest < defined min, then extrapolate smoothly between
    // TMINDAY_MLCS2k2 and Tdummag
    
    if ( Trest < TMINDAY_MLCS2k2 ) {
      mag_extrap_mlcs2k2(ifilt,delta,Trest, &rest_magval[i], &rest_magerr[i] );
      continue ;
    }

    // get j1,j2 date-indices that bracket Trest

    JDATES ( Trest, &j1, &j2 ) ; 

    // Now do simple linear interpolation by taking a weighted average

    w1 = fabs(TEMPLATE_DAYS_MLCS2k2[j2] - Trest );
    w2 = fabs(TEMPLATE_DAYS_MLCS2k2[j1] - Trest );
    
    get_MPQ_mlcs2k2(j1,ifilt, &M1, &P1, &Q1 );
    get_MPQ_mlcs2k2(j2,ifilt, &M2, &P2, &Q2 );

    tempmag1 = M1 + P1*delta + Q1*sqdelta ;
    tempmag2 = M2 + P2*delta + Q2*sqdelta ;

    rest_magval[i] = 
      (w1*tempmag1 + w2*tempmag2) / (w1 + w2);

    // now interploate the mag-error

    if ( USE_PEAKERR_ONLY == 0 ) {
      temperr1 = TEMPLATE_MAGS_MLCS2k2[j1][ifilt].MAGERR;
      temperr2 = TEMPLATE_MAGS_MLCS2k2[j2][ifilt].MAGERR;
    } 
    else {
      temperr1 = TEMPLATE_MAGS_MLCS2k2[IDAYPEAK_MLCS2k2][ifilt].MAGERR;
      temperr2 = TEMPLATE_MAGS_MLCS2k2[IDAYPEAK_MLCS2k2][ifilt].MAGERR;
    }

    rest_magerr[i] = 
      (w1*temperr1 + w2*temperr2) / (w1 + w2);


       idum = 1 ;

       if ( ifilt == -999 
	    && fabsf(Trest-44.50) < .01 
	    && fabs(delta+.0726)<.001 ) {

	 printf(" -------------------------------------------------- \n");
	 printf("Trest=%f temp_date1=%6.1f temp_date2=%6.1f %d,%d\n", 
		Trest,TEMPLATE_DAYS_MLCS2k2[j1],
		TEMPLATE_DAYS_MLCS2k2[j2],j1,j2);
   
	 printf("rest_mag=%f temp_mag[1,2]=%lf,%lf wgt[1,2]=%f,%f\n", 
		rest_magval[i],tempmag1,tempmag2,w1,w2); 
       }


   
  }  // end of observation loop

  return(SUCCESS);

}  // end of genmag_mlcs2k2


// ************************************************
void JDATES ( double Trest, int *j1, int *j2 ) {

  // return integer indices j1 and j2 such that
  // TEMPLATE_DAYS_MLCS2k2[j1,j2] bracket Trest

  int j;

  // ------------- BEGIN -----------

  *j1 = 0 ;  // init
  *j2 = 1 ;

  for(j=0; j < NDAY_MLCS2k2; j++ ) {
    if ( TEMPLATE_DAYS_MLCS2k2[j] > Trest ){
      *j1 = j-1 ;
      *j2 = j ;
      break ;
    }
  }

  if ( *j1 < 0 ) { *j1=0; *j1=1 ; }

}  // end of JDATES

// ***************************************
void mag_extrap_mlcs2k2(int ifilt, double delta, double Trest,
		     double *mag, double *magerr ) {

  // returns extrapolated *mag and *magerr when 
  // Trest is below min-defined range.
  // Extrapolation uses days defined by iepoch=0,1
  // (days -10 and -9)

  // extrapolate with exponential fall-off in magnitude
  // such that slope at Tminrest is the slope between
  // Tminrest and Tminrest + 1 day.


  double 
    arg
    ,sqdelta
    ,magmindef0,  errmindef0
    ,magmindef1,  errmindef1
    ,dmag_dTrest, derr_dTrest
    ;

  // ---------------- BEGIN ------------

  *mag    = 99.0;
  *magerr = 5.0;
  if ( Trest < -20.0 ) return ;

  sqdelta = delta * delta;

  // mag at min-defined Trest ... 

  magmindef0 
    = TEMPLATE_MAGS_MLCS2k2[0][ifilt].MAGOFF 
    + TEMPLATE_MAGS_MLCS2k2[0][ifilt].D1 * delta
    + TEMPLATE_MAGS_MLCS2k2[0][ifilt].D2 * sqdelta ;

  // and mag 1 day after min Trest

  magmindef1 
    = TEMPLATE_MAGS_MLCS2k2[1][ifilt].MAGOFF 
    + TEMPLATE_MAGS_MLCS2k2[1][ifilt].D1 * delta
    + TEMPLATE_MAGS_MLCS2k2[1][ifilt].D2 * sqdelta ;

  dmag_dTrest = (magmindef1 - magmindef0) ;
  arg         = (Trest - TMINDAY_MLCS2k2) * dmag_dTrest / magmindef0;
  *mag        = magmindef0 * exp(arg);

  // repeat for error.

  errmindef0 = TEMPLATE_MAGS_MLCS2k2[0][ifilt].MAGERR;
  errmindef1 = TEMPLATE_MAGS_MLCS2k2[1][ifilt].MAGERR;
  derr_dTrest = (errmindef1 - errmindef0) ;
  arg         = (Trest - TMINDAY_MLCS2k2) * derr_dTrest / errmindef0;
  *magerr     = errmindef0 * exp(arg);

}

// **************************************
void rd_mlcs2k2_cov(char *covFile, float scale_covar) {

  /***
      read covariance matrix from text file.
      Peakmag errors (with x4.1 scale_covar) are:

      Read MLCS covariance for U : mag-err at peak = 0.123 
      Read MLCS covariance for B : mag-err at peak = 0.102 
      Read MLCS covariance for V : mag-err at peak = 0.063 
      Read MLCS covariance for R : mag-err at peak = 0.059 
      Read MLCS covariance for I : mag-err at peak = 0.099 
      Read MLCS covariance for Y : mag-err at peak = XXXX 
      Read MLCS covariance for J : mag-err at peak = XXXX 
      Read MLCS covariance for H : mag-err at peak = XXXX
      Read MLCS covariance for K : mag-err at peak = XXXX 

    To save space, elements with zero covariance can be suppressed
    from the text file.

  **/


  char covFile_loc[100];
  char fnam[] = "rd_mlcs2k2_cov" ;
  FILE *fp ;

  char  cfilt1[2], cfilt2[2] ;
  float epoch1, epoch2, cov12, magerr ;
  int iep1, iep2, ifilt1, ifilt2, nline ;
  int ifilt1_last, GZIPFLAG ;

  // ------------ BEGIN -----------

  // if argument is null, then use default covariance matrix.

  if ( strlen(covFile) ==  0 ) 
    { sprintf(covFile_loc,"%s/smatrix.dat", PATHMODEL_MLCS2k2 ); }
  else
    { sprintf(covFile_loc,"%s", covFile ); }


  // init covariance matrix to zeros

  for (iep1=0; iep1 < MAXDAYS_MLCS2k2; iep1++ ) {
    for (iep2=0; iep2 < MAXDAYS_MLCS2k2; iep2++ ) {
      for (ifilt1=0; ifilt1 < MAXFILT_MLCS2k2; ifilt1++ ) {
	for (ifilt2=0; ifilt2 < MAXFILT_MLCS2k2; ifilt2++ ) {

	  MLCS2k2_COVAR[ifilt1][iep1][ifilt2][iep2] = 0.0 ;

	}
      }
    }
  }


  // Open file for reading 
  fp = open_TEXTgz(covFile_loc, "rt", &GZIPFLAG );
  if ( fp  == NULL ) {
    printf("\n ERROR in %s \n", fnam );
    printf("    Could not open '%s' \n ", covFile_loc );
    printf("    ***** ABORT ***** \n");
    exit(1);
  }
  
  printf("\n   Read  cov matrix from: \n\t %s \n", covFile_loc);
  printf("   Scale cov matrix by %f \n", scale_covar );


  nline = 0; 
  NONZERO_COVAR = 0;
  ifilt1_last = -1;

  while (fscanf(fp, "%s %f %s %f %f ", 
		cfilt1, &epoch1, cfilt2, &epoch2, &cov12 ) != EOF ) {

    cov12 *= scale_covar ;  // scale covariance matrix

    // get indices for filter and epoch.
    // Note that both indices start at zero.

    ifilt1 = ifilt_mlcs2k2(cfilt1);  // note ifilt starts at zero
    ifilt2 = ifilt_mlcs2k2(cfilt2);

    iep1 = (int)(epoch1 - TMINDAY_MLCS2k2 + 0.001 ) ; // index starts at zero
    iep2 = (int)(epoch2 - TMINDAY_MLCS2k2 + 0.001 ) ;


    // check for bad index

    if ( iep1 < 0 || iep1 >= MAXDAYS_MLCS2k2 ) {
      printf("\t ERROR: iep1=%d at epoch1=%f \n", iep1, epoch1 );
      printf("\t ***** ABORT ***** \n"); exit(1);
    }
    if ( iep2 < 0 || iep2 >= MAXDAYS_MLCS2k2 ) {
      printf("\t ERROR: iep2=%d at epoch2=%f \n", iep2, epoch2 );
      printf("\t ***** ABORT ***** \n"); exit(1);
    }
    if ( ifilt1 < 0 || ifilt1 >= NFILT_MLCS2k2 ) {
      printf("\t ERROR: ifilt1=%d at cfilt1=%s (NFILT_MLCS2k2=%d)\n", 
	     ifilt1, cfilt1, NFILT_MLCS2k2 );
      printf("\t ***** ABORT ***** \n"); exit(1);
    }
    if ( ifilt2 < 0 || ifilt2 >= NFILT_MLCS2k2 ) {
      printf("\t ERROR: ifilt2=%d at cfilt2=%s (NFILT_MLCS2k2=%d)\n",  
	     ifilt2, cfilt2, NFILT_MLCS2k2 );
      printf("\t ***** ABORT ***** \n"); exit(1);
    }

    // print comment when new filter is found.

    if ( ifilt1 != ifilt1_last ) {
      printf("\t Read MLCS2k2 covariance for %s : ", cfilt1 );
      fflush(stdout);
    }


    // print mag err at peak
    if ( epoch1 == 0.0 && epoch2 == 0.0 && ifilt1 == ifilt2 ) {
      magerr = sqrtf(cov12);
      printf("mag-err at peak = %5.3f \n", magerr );
   }

    nline++ ; 
    if ( cov12 != 0.0 ) NONZERO_COVAR++ ;

    ifilt1_last = ifilt1 ;

   /***
    printf(" %3d : %s(%d) at %6.1f(%2d)  %s(%d) at %6.1f(%2d) => cov12 = %f \n"
	   ,nline
	   ,cfilt1, ifilt1, epoch1, iep1
	   ,cfilt2, ifilt2, epoch2, iep2
	   ,cov12 );
    **/

    MLCS2k2_COVAR[ifilt1][iep1][ifilt2][iep2] = cov12 ;

    // fill errors (diagonal covariant elements only)

    if ( iep1 == iep2 && ifilt1 == ifilt2 ) 
      TEMPLATE_MAGS_MLCS2k2[iep1][ifilt1].MAGERR = sqrtf(cov12);

  }


  fclose(fp);

  printf("    Filled %d covariance elements (%d non-zero) \n", 
	 nline, NONZERO_COVAR );

  fflush(stdout);


} // end of rd_mlcs2k2_cov


// **********************************
int ifilt_mlcs2k2(char *cfilt) {

  // return filter index (0-8) from character name (UBVRIYJHK)

  int ifilt, ifilt_save ;

  ifilt_save = 0;

  for ( ifilt=0; ifilt < NFILT_MLCS2k2; ifilt++ ) {
    if ( cfilt[0] == FILTSTRING_MLCS2k2[ifilt] ) {
      ifilt_save = ifilt;
    }
  }

  return ifilt_save ;

} // end of ifilt_mlcs2k2


// *********************************
void dmp_mlcs2k2_cov(void) {

  // screen-dump normalized cov matrix at t=0
  // to show filter-correlations at peak.

  int ifilt_row, ifilt_col, iep0 ;
  double t0, cov, sig_row, sig_col, rho ;

  // ----------- BEGIN -----------

  if ( NONZERO_COVAR <= 0 ) return ;

  printf("\n   RHO(MLCS2k2) = COV[i,j]/sig(i)*sig(j) at Trest=0 : \n");

  // get epoch bin at peak
  t0  = (double)0.0 ;
  iep0 = (int)(t0 - TMINDAY_MLCS2k2 + 0.001 ) ;

  // make table header

  printf("    ");
  for ( ifilt_row=0; ifilt_row < NFILT_MLCS2k2 ; ifilt_row++ )
    printf("       %c", FILTSTRING_MLCS2k2[ifilt_row] );


  printf("\n       ");
  for ( ifilt_row=0; ifilt_row < NFILT_MLCS2k2 ; ifilt_row++ )
    printf("--------" );

  printf("\n");

  for ( ifilt_row=0; ifilt_row < NFILT_MLCS2k2 ; ifilt_row++ ) {

    printf("    %c : ", FILTSTRING_MLCS2k2[ifilt_row] );

    for ( ifilt_col=0; ifilt_col < NFILT_MLCS2k2 ; ifilt_col++ ) {
      cov     = MLCS2k2_COVAR[ifilt_row][iep0][ifilt_col][iep0] ;
      sig_row = sqrt( MLCS2k2_COVAR[ifilt_row][iep0][ifilt_row][iep0] );
      sig_col = sqrt( MLCS2k2_COVAR[ifilt_col][iep0][ifilt_col][iep0] );

      if ( sig_row != 0.0 && sig_col != 0.0 ) 
	rho     = cov / ( sig_row * sig_col );
      else
	rho = 0.0;

      printf("%7.4f ", rho);
    }
    printf("\n");
  }

  printf("       ");
  for ( ifilt_row=0; ifilt_row < NFILT_MLCS2k2 ; ifilt_row++ )
    printf("--------" );

  printf("\n");
  fflush(stdout);


}  // end of dmp_mlcs2k2_cov



// *********************************
void dmp_mlcs2k2_magerr(void) {

  // screen-dump UBVRI magerrs at each epoch.
  // format is paw-readable:
  //    Trest  U-err  B-err  V-err  R-err   I-err
  // Use the genmag_mlcs2k2 function as part of the test.

  int ifilt, iep; 
  double Trest, mag, err; 
  double delta = 0.0;

  // ----------- BEGIN -----------


  // make table header
  printf("\n               MLCS2k2 MAGERR for \n");

  printf(" EPOCH  ");
  for ( ifilt=0; ifilt < NFILT_MLCS2k2 ; ifilt++ )
    printf("       %c", FILTSTRING_MLCS2k2[ifilt] );

  printf("\n    ");
  for ( ifilt=0; ifilt < NFILT_MLCS2k2 ; ifilt++ )
    printf("--------" );

  printf("\n");


  for (iep= -20; iep < 90; iep++ ) {

    Trest = (double)iep ;
    printf(" %8.2f ", Trest );

    for ( ifilt=0; ifilt < NFILT_MLCS2k2 ; ifilt++ ) {

      // fetch mag and error (ignore mag)

      genmag_mlcs2k2(ifilt, delta, 1, &Trest, &mag, &err );
      printf(" %6.3f " , err );

    }  // end of ifilt loop
    printf("\n");

  }   // end of epoch loop


  printf("\n      ");
  for ( ifilt=0; ifilt < NFILT_MLCS2k2 ; ifilt++ )
    printf("--------" );

  printf("\n");


  fflush(stdout);


}  // end of dmp_mlcs2k2_magerr


// ***********************************
int gencovar_mlcs2k2 (
		  int MATSIZE          // (I) matrix size (row or col)
		  ,int *ifilt          // (I) list of filter indices
		  ,double *rest_epoch  // (I) list of rest-epochs (days)
		  ,double *covar       // (O) covariance matrix
		  ) {

  /******

  Created Aug 10, 2006 by R.Kessler
  Returns covariance matrix for specified filters and epochs.
  Example:
     matsize    = 6
     ifilt      = 2 2 2 3 3 3  (B B B V V V)
     rest_epoch = 3 6 7 3 6 7 

  is for rest epochs 3,6,7 days in filters B and V.

  Then the returned 6x6 matrix is returned as a linear array:

  (1st row)
     *covar[00] = C11(B,epoch1 : B,epoch1)
     *covar[01] = C12(B,epoch1 : B,epoch2)
     *covar[02] = C13(B,epoch1 : B,epoch3)
     *covar[03] = C14(B,epoch1 : V,epoch1)
     *covar[04] = C15(B,epoch1 : V,epoch2)
     *covar[05] = C16(B,epoch1 : V,epoch3)

 (2nd row)
     *covar[06] = C21(B,epoch2 : B,epoch1)
     *covar[07] = C22(B,epoch2 : B,epoch2)
      ...

  ****/

  int 
     icovar
    ,irow, icol
    ,ifilt_row, ifilt_col
    ,iep1_row, iep1_col
    ,iep2_row, iep2_col
    ,j1, j2
    ;

  double 
    ep_row
    ,ep_col
    ,wgt[4][4]
    ,cov[4][4]
    ,dif1_row
    ,dif2_row
    ,dif1_col
    ,dif2_col
    ,covsum
    ,wsum
    ;

  // --------------- BEGIN -------------


  icovar = 0;

  for ( irow=0; irow < MATSIZE; irow++ ) {
    
    ep_row    = *(rest_epoch + irow);  // epoch,  days
    if ( ep_row < TMINDAY_MLCS2k2 )  ep_row = TMINDAY_MLCS2k2;
    if ( ep_row > TMAXDAY_MLCS2k2 )  ep_row = TMAXDAY_MLCS2k2 - 0.01 ;

    ifilt_row = *(ifilt + irow) ;     // filter index

    // get epoch indices that bracket ep_row
    JDATES ( ep_row, &iep1_row, &iep2_row ) ; 

    dif1_row = fabs ( TEMPLATE_DAYS_MLCS2k2[iep1_row] - ep_row ) ;
    dif2_row = fabs ( TEMPLATE_DAYS_MLCS2k2[iep2_row] - ep_row ) ;

    for ( icol=0; icol < MATSIZE; icol++ ) {

      ep_col    = *(rest_epoch + icol);  // epoch,  days
      if ( ep_col < TMINDAY_MLCS2k2 )  ep_col = TMINDAY_MLCS2k2;
      if ( ep_col > TMAXDAY_MLCS2k2 )  ep_col = TMAXDAY_MLCS2k2 - 0.01;

      ifilt_col = *(ifilt + icol);     // filter index

      // get epoch indices that bracket ep_col
      JDATES ( ep_col, &iep1_col, &iep2_col ) ; 

      dif1_col = fabs ( TEMPLATE_DAYS_MLCS2k2[iep1_col] - ep_col ) ;
      dif2_col = fabs ( TEMPLATE_DAYS_MLCS2k2[iep2_col] - ep_col ) ;

      // now we have ep_row,ep_col in the middle of the square grid
      // containing this point. Take average of 4 corners wgted
      // by distance from ep_row,ep_col to corner

      wgt[1][1] = dif2_row * dif2_col ;
      wgt[1][2] = dif2_row * dif1_col ;
      wgt[2][1] = dif1_row * dif2_col ;
      wgt[2][2] = dif1_row * dif1_col ;


      cov[1][1] = MLCS2k2_COVAR[ifilt_row][iep1_row][ifilt_col][iep1_col] ;
      cov[1][2] = MLCS2k2_COVAR[ifilt_row][iep1_row][ifilt_col][iep2_col] ;
      cov[2][1] = MLCS2k2_COVAR[ifilt_row][iep2_row][ifilt_col][iep1_col] ;
      cov[2][2] = MLCS2k2_COVAR[ifilt_row][iep2_row][ifilt_col][iep2_col] ;

      wsum   = 0.0 ;
      covsum = 0.0 ;
      for ( j1=1; j1 <= 2; j1++ ) {
	for ( j2=1; j2 <= 2; j2++ ) {
	  covsum += wgt[j1][j2] * cov[j1][j2] ;
	  wsum   += wgt[j1][j2] ;
	}
      }

      if ( irow == -6 && icol == -4 ) {
	printf(" ----------- row,col = %d %d ------------ \n", irow, icol );

	printf(" ifilt_row,col = %d %d \n", ifilt_row, ifilt_col );

	printf(" ep_row,col = %f %f \n", ep_row, ep_col );
	printf(" dif_row = %f,%f    dif_col = %f,%f\n", 
	       dif1_row, dif2_row, dif1_col, dif2_col ) ;
   
	printf(" iep_row = %d,%d   iep_col=%d,%d \n", 
	       iep1_row, iep2_row, iep1_col, iep2_col );

	printf(" wgt = %10.6f %10.6f %10.6f %10.6f \n", 
	       wgt[1][1] , wgt[1][2] , wgt[2][1] ,  wgt[2][2] );

	printf(" COV = %10.6f %10.6f %10.6f %10.6f \n", 
	       cov[1][1] , cov[1][2] , cov[2][1] ,  cov[2][2] );
	     
	printf(" wsum=%f  covsum=%f  <COV>=%f \n",  
	       wsum, covsum, covsum/wsum);

      }

      *(covar + icovar) = covsum / wsum ;

      icovar++ ;      

    }  // end of icol
  }  // end of irow


  return 1;

} // end of gennmag_covar


// ==============================================
// misc. functions to return info.

int mlcs2k2_Tmin(void) {  return (int)TMINDAY_MLCS2k2 ; }
int mlcs2k2_Tmax(void) {  return (int)TMAXDAY_MLCS2k2 ; }

// ==============================================
void set_LAMRANGE_mlcs2k2(void) {
  // Set rest-frame lambda range (internal use only)
  char tempFilename[200];
  char restlamFile[40] = "RESTLAMBDA_RANGE.DAT" ;
  FILE *fp ;

  // set default lambda range for UBVRI
  LAMRANGE_MLCS2k2[0] = MINLAM_MLCS2k2_U ;  // start with hard-coded values
  LAMRANGE_MLCS2k2[1] = MAXLAM_MLCS2k2_I ;

  // if RESTLAMBDA_RANGE.DAT exists, use that range instead

  sprintf(tempFilename,"%s/%s", PATHMODEL_MLCS2k2, restlamFile );
  if (( fp = fopen(tempFilename, "rt")) != NULL ) {
    printf("\t Read RESTLAMBDA-range from \n\t   %s/%s\n", 
	   "$PATHMODEL_MLCS2k2", restlamFile );
    readdouble(fp, 2, LAMRANGE_MLCS2k2 );
    fclose(fp);
  }
}

// ==============================================
void get_MPQ_mlcs2k2(int j, int ifilt, 
		     double *M, double *P, double *Q) {  // return args
  
  // Created March 12, 2010 by R.Kessler
  // inputs: j = day index; ifilt = filter index
  // Returns M, P, Q model parameters

    *M = TEMPLATE_MAGS_MLCS2k2[j][ifilt].MAGOFF + MOFF_MLCS2k2 ;
    *P = TEMPLATE_MAGS_MLCS2k2[j][ifilt].D1;
    *Q = TEMPLATE_MAGS_MLCS2k2[j][ifilt].D2 * QWGT_MLCS2k2 ;

} // end of get_MPQ_mlcs2k2 

// =================================================		  
double get_QWGT_mlcs2k2(double Trest, double delta) {

  // Created March 12, 2010 by R.Kessler

  double QTMP, dif, sig, arg;

  // -------------- BEGIN -----------

  QTMP = 1.0;

  if ( QWGT_FLAG == 0 )        return QTMP ;
  if ( Trest < QWGT_TRESTMIN ) return QTMP ;
  if ( delta < QWGT_DELTAMIN ) return QTMP ;

  // Use half-Gaussian to suppress Q
  dif = (Trest - QWGT_TRESTMIN);
  sig = QWGT_TRESTSIG ;
  arg = 0.5 * (dif*dif)/(sig*sig) ;
  QTMP = exp(-arg);

  //  printf(" xxxx QTMP(Trest=%6.2f) = %f \n", Trest, QTMP );

  return QTMP; 


} // end of QWGT_mlcs2k2

// ==============================================
void get_LAMRANGE_mlcs2k2(double *lammin, double *lammax) {
  // this function is used by external programs
  *lammin = LAMRANGE_MLCS2k2[0];
  *lammax = LAMRANGE_MLCS2k2[1];
}

