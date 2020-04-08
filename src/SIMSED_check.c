/**************************************************
  Created July 18 2018 by R.Kessler

  Read SIMSED model and run diagnostic checks.
  + Report FATAL errors such as inconsistent binning,
  + Report WARNINGS such as extreme mag/day changes,
     and Flux(last Day) not falling enough w.r.t., peak flux.
  + report memory requirememts
  + report read time
 
  Example:
    SIMSED_check.exe <inpDir>
    SIMSED_check.exe <inpDir>  -lamrange 1000 11000

 ***********/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "MWgaldust.h"

#include "genmag_SIMSED.h"
#include "genmag_SEDtools.h"


// ========= GLOBAL DECLARATIONS ============

// give warnings for:
#define WARNVAL_MAGDAY_SLOPE1    2.0  // warning of mag/day 
#define WARNVAL_MAGDAY_SLOPE2    0.5  // warning of d2mag/d2day (2nd deriv.)
#define WARNVAL_MAGDIF_LASTPEAK  5.0  //  mag(last) - mag(peak) < 5  

// parameters used to set warnings
//#define WARNPAR_DAYDIF      1.0  // require at least one day slope warnings
#define WARNPAR_PEAKFRAC    0.05 // this fract of peak flux for magslope test  

#define BOXFILT_WIDTH_APPROX 1000.0  // aim for this width
#define MXBOXFILT 50   // max number of box filters to define

// NERR counts once per SED so that the grand summary
// tells how many SEDs had this error

int NWARNTMP_MAGDAY_SLOPE1;    // number of warnings for current SED
int NWARNTMP_MAGDAY_SLOPE2;    
int NWARNTMP_MAGDIF_LASTPEAK ;

int NWARNTOT_MAGDAY_SLOPE1;    // total number of warnings
int NWARNTOT_MAGDAY_SLOPE2; 
int NWARNTOT_MAGDIF_LASTPEAK ;

int NSEDWARN_MAGDAY_SLOPE1;    // number of SED files with warning
int NSEDWARN_MAGDAY_SLOPE2;
int NSEDWARN_MAGDIF_LASTPEAK ;

struct {
  char   INPDIR_SIMSED[MXPATHLEN];
  double LAMRANGE[2];
} INPUTS;


int     NBOXFILT ;
double  LAMRANGE_BOXFILT[MXBOXFILT][2];
double  LAMBIN_BOXFILT;

struct {
  int    IDAY_PEAK;
  double FLUX_FIRST, FLUX_LAST, FLUX_PEAK;
  double FLUX[MXBIN_DAYSED_SEDMODEL] ;
} BOXFILT[MXBOXFILT] ;


// -----------
void  parse_args(int argc, char **argv);
void  SIMSED_check_DRIVER(int ised);
void  define_boxFilters(void);
void  load_FLUX_BOXFILT(void);
void  SIMSED_check_slope(int ised, int ifilt);
void  SIMSED_check_lastDay(int ised, int ifilt);
void  SIMSED_WARNSTRING(int ised, int ifilt, char *string) ;
void  SIMSED_warning_summary(void);

// ====================================
int main(int argc, char **argv) {

  int OPTMASK = OPTMASK_SIMSED_TESTMODE;
  char fnam[] = "main" ;
  
  // ------------ BEGIN ------------

  set_EXIT_ERRCODE(EXIT_ERRCODE_UNKNOWN);

  sprintf(BANNER,"Begin execution of SIMSED_check" );
  print_banner(BANNER);

  // parse user arguments
  parse_args(argc, argv );

  init_genmag_SIMSED(INPUTS.INPDIR_SIMSED, "", "", "", OPTMASK);

  NSEDWARN_MAGDAY_SLOPE1    = 0 ;
  NSEDWARN_MAGDAY_SLOPE2    = 0 ;
  NSEDWARN_MAGDIF_LASTPEAK = 0 ;
  NWARNTOT_MAGDAY_SLOPE1   = 0 ;
  NWARNTOT_MAGDAY_SLOPE2   = 0 ;
  NWARNTOT_MAGDIF_LASTPEAK = 0 ;

  int ised, NDAY, NLAM ;
  char *sedFile, SEDFILE[MXPATHLEN] ;
  char comment[60] ;
  for(ised=1; ised <= SEDMODEL.NSURFACE; ised++ ) {
    sedFile = SEDMODEL.FILENAME[ised];
    //    printf(" Check %s \n", sedFile);

    sprintf(SEDFILE, "%s/%s", INPUTS.INPDIR_SIMSED, sedFile);
    sprintf(comment,"(%d of %d)", ised, SEDMODEL.NSURFACE);
    read_SIMSED_flux(SEDFILE,comment);
    
    if ( TEMP_SEDMODEL.LAMMIN < INPUTS.LAMRANGE[0] ) 
      { TEMP_SEDMODEL.LAMMIN  = INPUTS.LAMRANGE[0] ; }

    if ( TEMP_SEDMODEL.LAMMAX > INPUTS.LAMRANGE[1] ) 
      { TEMP_SEDMODEL.LAMMAX  = INPUTS.LAMRANGE[1] ; }

    printf("\t %.1f < DAY < %.1f     %.1f < LAM < %.1f \n",
	   TEMP_SEDMODEL.DAYMIN, TEMP_SEDMODEL.DAYMAX,
	   TEMP_SEDMODEL.LAMMIN, TEMP_SEDMODEL.LAMMAX );

    SIMSED_check_DRIVER(ised);
  } 

  SIMSED_warning_summary();
  printf("\n Done. \n"); fflush(stdout);

} // end main


// ******************************************
void parse_args(int argc, char **argv) {

  int i, ipar,  ep, ilast, iuse ;
  char fnam[] = "parse_args" ;

  // ---------- BEGIN --------

  NARGV_LIST = argc;
  for( i=0; i < NARGV_LIST; i++ ) {
    USE_ARGV_LIST[i] = 0 ;
    sprintf(ARGV_LIST[i], "%s", argv[i] );
  }
  
  i=1; ilast=i ;
  sscanf(ARGV_LIST[i], "%s", INPUTS.INPDIR_SIMSED );
  USE_ARGV_LIST[i] = 1;
  i++; ilast=i ;

  INPUTS.LAMRANGE[0] =  500.0 ;
  INPUTS.LAMRANGE[1] = 25000.0 ;

  while ( i < NARGV_LIST ) {

    
    if ( strcmp(ARGV_LIST[i],"-lamrange") == 0 ) {
      i++;  sscanf(ARGV_LIST[i], "%le", &INPUTS.LAMRANGE[0] );
      i++;  sscanf(ARGV_LIST[i], "%le", &INPUTS.LAMRANGE[1] );
    }


    if ( i > ilast ) {
      for ( iuse = ilast; iuse <= i; iuse++ )
        { USE_ARGV_LIST[iuse] = 1; }
    }
    i++ ;  ilast=i;

  } 
  
  check_argv();

  ENVreplace("init",fnam,1);
  ENVreplace(INPUTS.INPDIR_SIMSED,fnam,1);

  printf("   INPDIR_SIMSED: %s \n", INPUTS.INPDIR_SIMSED );
  printf("   Test LAMRANGE: %.1f to %.1f A \n",
	 INPUTS.LAMRANGE[0], INPUTS.LAMRANGE[1] );
  printf("\n\n");

  return ;

} // end parse_args


// =====================================
void SIMSED_check_DRIVER(int ised) {

  int iday, ifilt ;
  char fnam[] = "SIMSED_check_DRIVER" ;

  // --------- BEGIN ---------

  NWARNTMP_MAGDAY_SLOPE1    = 0 ;
  NWARNTMP_MAGDAY_SLOPE2    = 0 ;
  NWARNTMP_MAGDIF_LASTPEAK  = 0 ;

  // define box filters so that each filter width is approx 1000 A
  define_boxFilters();

  // load flux for each day and each box filter
  load_FLUX_BOXFILT();


  for(ifilt = 0; ifilt < NBOXFILT; ifilt++  ) {
    SIMSED_check_lastDay(ised,ifilt);
    SIMSED_check_slope(ised, ifilt); //  check mag/day
  }

  if ( NWARNTMP_MAGDAY_SLOPE1   > 0 ) { NSEDWARN_MAGDAY_SLOPE1++ ; }
  if ( NWARNTMP_MAGDAY_SLOPE2   > 0 ) { NSEDWARN_MAGDAY_SLOPE2++ ; }
  if ( NWARNTMP_MAGDIF_LASTPEAK > 0 ) { NSEDWARN_MAGDIF_LASTPEAK++ ; }


  return ;
} // end SIMSED_check_DRIVER


// =====================================
void define_boxFilters(void) {

  int ifilt;
  double xtmp, xi, LAMRANGE, LAMBIN ;
  char fnam[] = "define_boxFilters" ;

  // ----------- BEGIN ----------

  // start by definine box filters
  LAMRANGE = (TEMP_SEDMODEL.LAMMAX - TEMP_SEDMODEL.LAMMIN) ;
  xtmp     = LAMRANGE/BOXFILT_WIDTH_APPROX ; ;
  NBOXFILT = (int)(xtmp+0.5);
  LAMBIN   = LAMRANGE/(double)NBOXFILT ;
  LAMBIN_BOXFILT = LAMBIN ; // set global

  if ( NBOXFILT >= MXBOXFILT ) {
    sprintf(c1err,"NBOXFILT=%d exceeds boundf MXBOXFILT=%d",
	    NBOXFILT, MXBOXFILT );
    sprintf(c2err,"Check LAM range, or update MXBOXFILT");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  //  printf(" xxx NBOXFILT=%d , LAMBIN=%.2f \n",  NBOXFILT,LAMBIN);

  for(ifilt=0; ifilt < NBOXFILT; ifilt++ ) {
    xi = (double)ifilt ;
    LAMRANGE_BOXFILT[ifilt][0] = TEMP_SEDMODEL.LAMMIN + LAMBIN*xi;
    LAMRANGE_BOXFILT[ifilt][1] = LAMRANGE_BOXFILT[ifilt][0] + LAMBIN ;
  }


  return ;

} // end define_boxFilters



// =====================================
void  load_FLUX_BOXFILT(void) {

  // compute and store FLUXes in BOXFILT struct.

  int NDAY = TEMP_SEDMODEL.NDAY; 
  int NLAM = TEMP_SEDMODEL.NLAM; 
  double LAMMIN = TEMP_SEDMODEL.LAMMIN ;
  double LAMBIN = LAMBIN_BOXFILT ;

  double FLUXTMP, FLUXMAX[MXBOXFILT], LAM ;
  int iday, ifilt, ilam, jflux ;
  char fnam[] = "load_FLUX_BOXFILT" ;

  // ---------- BEGIN -----------

  for(ifilt=0; ifilt < NBOXFILT; ifilt++ )  { 
    BOXFILT[ifilt].IDAY_PEAK  = -999 ;
    BOXFILT[ifilt].FLUX_FIRST = -9999.0 ;
    BOXFILT[ifilt].FLUX_LAST  = -9999.0 ;
    BOXFILT[ifilt].FLUX_PEAK  = -9999.0 ;
    for(iday = 0 ; iday < NDAY; iday++ ) { BOXFILT[ifilt].FLUX[iday] = 0.0 ; }
  }
  
  //  double DAY_PEAK;
  //  double FLUX_FIRST, FLUX_LAST, FLUX_PEAK;
  //  double FLUX[MXBIN_DAYSED_SEDMODEL] ;

  for(iday = 0 ; iday < NDAY; iday++ ) {

    // compute flux[ifilt] for this day
    for(ilam=0 ; ilam < NLAM; ilam++ ) {    
      jflux   = NLAM*iday + ilam;
      FLUXTMP = TEMP_SEDMODEL.FLUX[jflux];
      FLUXTMP *= SEDMODEL.FLUXSCALE ;
      // figure out which of the box filters
      LAM = TEMP_SEDMODEL.LAM[ilam];
      ifilt = (int)( (LAM-LAMMIN)/LAMBIN ) ;
      BOXFILT[ifilt].FLUX[iday] += FLUXTMP ;
    } // end ilam
      

    // keep track of first, last and peak flux
    for(ifilt = 0; ifilt < NBOXFILT; ifilt++  ) {
      FLUXTMP = BOXFILT[ifilt].FLUX[iday] ;
      
      if ( iday==0 ) 
	{ BOXFILT[ifilt].FLUX_FIRST = FLUXTMP ; }
      
      if ( iday==NDAY-1 ) 
	{ BOXFILT[ifilt].FLUX_LAST   = FLUXTMP ; }
      
      if ( FLUXTMP > BOXFILT[ifilt].FLUX_PEAK ) {
	BOXFILT[ifilt].FLUX_PEAK = FLUXTMP ;
	BOXFILT[ifilt].IDAY_PEAK = iday ;
      }
    }

  }   // iday

  int LDMP=0;

  if ( LDMP ) {
    printf(" xxx \n");
    printf(" xxx                  FLUX      FLUX      FLUX        \n");
    printf(" xxx    BOXFILT       FIRST     PEAK      LAST    IDAY_PEAK \n");
    printf(" xxx  ----------------------------------------------------- \n");
    for(ifilt = 0; ifilt < NBOXFILT; ifilt++  ) {
      printf(" xxx %7.1f-%7.1f  %9.3le %9.3le %9.3le  %3d \n",
	     LAMRANGE_BOXFILT[ifilt][0], LAMRANGE_BOXFILT[ifilt][1], 
	     BOXFILT[ifilt].FLUX_FIRST,	BOXFILT[ifilt].FLUX_PEAK,
	     BOXFILT[ifilt].FLUX_LAST,  BOXFILT[ifilt].IDAY_PEAK );
      fflush(stdout);
    }
    debugexit(fnam);
  } // end LDMP

  return ;

} // end load_FLUX_BOXFILT


// =====================================
void SIMSED_check_lastDay(int ised, int ifilt) {

  // check mag(last) - mag(peak) > WARN_MAGDIF_LASTPEAK

  double FLUX_PEAK = BOXFILT[ifilt].FLUX_PEAK ;
  double FLUX_LAST = BOXFILT[ifilt].FLUX_LAST ;
  double MAGDIF;
  char warnString[80];
  char fnam[] = "SIMSED_check_lastDay" ;

  // ------------ BEGIN --------

  if ( FLUX_LAST <= 0.0 ) { FLUX_LAST = 1.0E-30; }
  
  MAGDIF = -2.5*log10(FLUX_LAST/FLUX_PEAK);
  if ( MAGDIF < WARNVAL_MAGDIF_LASTPEAK ) {
    SIMSED_WARNSTRING(ised,ifilt, warnString);    
    printf("    %s: mag(Last)-mag(peak)=%.1f < %.1f \n",
	   warnString, MAGDIF, WARNVAL_MAGDIF_LASTPEAK );
    fflush(stdout);
    NWARNTMP_MAGDIF_LASTPEAK++; 
    NWARNTOT_MAGDIF_LASTPEAK++; 
  }

  return ;

} // end  SIMSED_check_lastDay

// ==============================================
void SIMSED_check_slope(int ised, int ifilt) {
  
  // check mag/day slopes

  int NDAY = TEMP_SEDMODEL.NDAY; 
  int JMINLAM  = (int)LAMRANGE_BOXFILT[ifilt][0] ;
  int JMAXLAM  = (int)LAMRANGE_BOXFILT[ifilt][1] ;
  double FPEAK =  BOXFILT[ifilt].FLUX_PEAK ;

  int    IDAY, iday, SKIP ;
  double FLUX[3], DAY[3], MAGDIF01, MAGDIF12;
  double DAYDIF01, DAYDIF12, DAYAVG01, DAYAVG12 ;
  double MAGSLOPE01, MAGSLOPE12, D2MAGD2DAY ;
  char warnString[80] ;
  char fnam[] = "SIMSED_check_slope";

  // ------------- BEGIN ---------


  for(IDAY=0; IDAY < NDAY-2; IDAY++ ) {

    SKIP = 0 ;
    for(iday=0; iday <=2; iday++ ) {
      FLUX[iday] = BOXFILT[ifilt].FLUX[IDAY+iday];
      DAY[iday]  = TEMP_SEDMODEL.DAY[IDAY+iday] ;      
      if ( FLUX[iday] == 0.0 ) { SKIP=1; }
      if ( FLUX[0]/FPEAK < WARNPAR_PEAKFRAC ) { SKIP=1; }
      if ( FLUX[1]/FPEAK < WARNPAR_PEAKFRAC ) { SKIP=1; }
    }
    if ( SKIP ) { continue ; }

    MAGDIF01   = -2.5*log10(FLUX[1]/FLUX[0]);
    MAGDIF12   = -2.5*log10(FLUX[2]/FLUX[1]);
    DAYDIF01   =  DAY[1] - DAY[0];
    DAYDIF12   =  DAY[2] - DAY[1];
    DAYAVG01   = (DAY[1] + DAY[0])/2.0;
    DAYAVG12   = (DAY[2] + DAY[1])/2.0;
    MAGSLOPE01 = MAGDIF01/DAYDIF01 ;
    MAGSLOPE12 = MAGDIF12/DAYDIF12 ;
    D2MAGD2DAY = fabs(MAGSLOPE12 - MAGSLOPE01) / (DAYAVG12-DAYAVG01);

    // skip rising part
    if ( MAGSLOPE01 < 0.0 && MAGSLOPE12 < 0.0 ) { continue ; }

    // xxxxxxxx
    if ( fabs(DAY[0] - 999303.3)<.2 ) {
      printf(" xxx -------------------------------- \n");
      printf(" xxx DUMP ised=%d ifilt=%d \n", ised, ifilt);
      printf(" xxx DAY[0,1,2] = %.1f, %.1f, %.1f \n",
	     DAY[0], DAY[1], DAY[2]);
      printf(" xxx FLUX[0,1,2] = %le, %le, %le \n",
	     FLUX[0], FLUX[1], FLUX[2]);
      printf(" xxx DAYAVG01=%.1f  DAYAVG12=%.1f \n", 
	     DAYAVG01, DAYAVG12);
      printf(" xxx SLOPE01=%.1f  SLOPE12=%.1f \n", 
	     MAGSLOPE01, MAGSLOPE12);
      printf(" xxx d^2Mag/d^2Day = %.3f \n", D2MAGD2DAY );
    }
    // xxxxxx

    if ( fabs(D2MAGD2DAY) > WARNVAL_MAGDAY_SLOPE2 ) {
      SIMSED_WARNSTRING(ised,ifilt, warnString);
      printf("  %s: d^2 = %.2f/%.2f = %.2f (DAY=%.2f, PeakFrac=%.2f) \n", 
	     warnString, MAGSLOPE12-MAGSLOPE01, DAYAVG12-DAYAVG01,
	     D2MAGD2DAY, DAY[0], FLUX[0]/FPEAK );
      NWARNTMP_MAGDAY_SLOPE2++ ;
      NWARNTOT_MAGDAY_SLOPE2++ ;
    }

  }

  return;

} // end SIMSED_check_slope


/* xxxx
// ==============================================
void SIMSED_check_slope(int ised, int ifilt) {
  
  // check mag/day slopes

  int NDAY = TEMP_SEDMODEL.NDAY; 
  int JMINLAM  = (int)LAMRANGE_BOXFILT[ifilt][0] ;
  int JMAXLAM  = (int)LAMRANGE_BOXFILT[ifilt][1] ;
  double FPEAK =  BOXFILT[ifilt].FLUX_PEAK ;

  int iday ;
  double FLAST, FNOW, DAYNOW, DAYLAST, DAYDIF, MAGDIF;
  double MAGSLOPE1, MAGSLOPE2, MAGSLOPE1_LAST ;
  char warnString[80] ;
  char fnam[] = "SIMSED_check_slope";

  // ------------- BEGIN ---------

  FLAST = BOXFILT[ifilt].FLUX[0];  
  DAYLAST  = TEMP_SEDMODEL.DAY[0];

  for(iday=1; iday < NDAY; iday++ ) {

    DAYNOW   = TEMP_SEDMODEL.DAY[iday] ;
    DAYDIF   = DAYNOW - DAYLAST ;
    if ( DAYDIF < WARNPAR_DAYDIF ) { continue ; }

    FNOW  = BOXFILT[ifilt].FLUX[iday];
    if ( FLAST/FPEAK < WARNPAR_PEAKFRAC ) { goto STORE_LAST ; }
    if ( FNOW/FPEAK  < WARNPAR_PEAKFRAC ) { goto STORE_LAST ; }

    //    if ( FLAST == 0 || FNOW == 0 ) { continue ; }

    MAGDIF   = -2.5*log10(FNOW/FLAST);
    MAGSLOPE1 = MAGDIF/DAYDIF ;

    if ( fabs(MAGSLOPE1) > WARNVAL_MAGDAY_SLOPE1 ) {
      SIMSED_WARNSTRING(ised,ifilt, warnString);
      printf("    %s: MAG/DAY = %.2f/%.2f = %.2f (DAY=%.2f) \n", 
	     warnString, MAGDIF, DAYDIF, MAGSLOPE1, DAYNOW );
      NWARNTMP_MAGDAY_SLOPE1++ ;
      NWARNTOT_MAGDAY_SLOPE1++ ;
    }

  STORE_LAST:
    FLAST = FNOW ;   DAYLAST = DAYNOW ;
  }

  return;

} // end SIMSED_check_slope
xxx */

void SIMSED_WARNSTRING(int ised, int ifilt, char *string) {
  //  return string = "WARNING[ised=%d, %d<LAM<%d]"
  int JMINLAM = (int)LAMRANGE_BOXFILT[ifilt][0] ;
  int JMAXLAM = (int)LAMRANGE_BOXFILT[ifilt][1] ;
  sprintf(string, "WARNING[ised=%d, %d<LAM<%d]", 
	  ised, JMINLAM, JMAXLAM );  
  return ;
} // end SIMSED_WARNSTRING 


void  SIMSED_warning_summary(void) {

  int NSED = SEDMODEL.NSURFACE;
  printf("\n\t\t WARNING SUMMARY \n");

  printf(" Number of SED files with: \n"); 

  printf("\t MAG/DAY warning: %d of %d (NwarnTOT=%d)\n",
	 NSEDWARN_MAGDAY_SLOPE1, NSED, NWARNTOT_MAGDAY_SLOPE1 );

  printf("\t d^2MAG/d^2DAY warning: %d of %d (NwarnTOT=%d)\n",
	 NSEDWARN_MAGDAY_SLOPE2, NSED, NWARNTOT_MAGDAY_SLOPE2 );

  printf("\t MAG(LAST)-MAG(PEAK) warning: %d of %d (NwarnTOT=%d) \n",
	 NSEDWARN_MAGDIF_LASTPEAK, NSED, NWARNTOT_MAGDIF_LASTPEAK  );
  
} // end SIMSED_warning_summary
