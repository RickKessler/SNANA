/*****************************************


  Created Mar 12, 2009 by R.Kessler

  Simulate (actually, calculate) the filter-response
  as determined by a calibration system that scans
  with a particular wavelength step and gaussian sigma.
  Then connect the dots from calib-measurement and 
  compute synthetic mag. Finally, compare synthetic
  mag from true filter response to that from calibration.
  
  Output is a text file defined by the user,
  and the columns are explained in <outfile>_README

  USAGE: 
   filtercal_sim.exe  <input file>



  HISTORY
  ~~~~~~~~~~~~~


*****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"

// =======================================
//   declare function prototypes
// =======================================

void  parse_args(int argc, char **argv);
void  read_input(void);
void  proc_input(void);

void  read_filterTrans(int ifile);
void  read_sed(int ifile);

void filtercal_init(void);
void filtercal_end(void);
void filtercal_rebin(int NBIN, double *ptrlam, double *ptrfun, 
		     double *ptrfun_rebin);
void filtercal_lampspec(void);  // init LAMP spectrum

void open_output(void);

double SYNMAG( double z, double *ptrlam, double *ptrtrans, double *ptrflux );

void MAGEVAL(int ised, int ifilt, int iz, int isig, int istp, int ioff);

void FILTCAL_CURVE (
		    int ifilt, 
		    double LAM_STEP, double LAM_SIGMA, double LAM_START,
		    int *NBIN_FILTCAL, double *FILTCAL_LAMBDA, 
		    double *FILTCAL_TRANS, double *LAM_AVG ); 

double FILTCAL_POINT( int ifilt, double LAM_CEN, double LAM_SIGMA ) ;

void   DUMP_FILTER(int ifilt, int isig, int istp,  double *TRANS_REBIN );

// =======================================
//   declare global variables
// =======================================

#define MXFILTCAL    20    // max # defined filters
#define MXSEDCAL     20    // max # defined SEDs
#define MXLAM_FILT   1500  // max # input lambda bins for filters
#define MXLAM_SED    4000  // max # input lambda bins for SED
#define MXLAMINT    20000  // max # lambda bins for integration
#define MXBIN         100  // max # bins for Z, LAMSTEP, LAMSIGMA

// xxx mark #define FNU_AB  3.631E-20 // flat Fnu for AB, erg/cm^2*s*Hz

double GLOBAL_LAMBDA_MIN; //  global lambda min,max among filters;
double GLOBAL_LAMBDA_MAX; //  used to limit SED range that is read


FILE *fpout ;
double  MAGDIF_MAX ;


struct INPUTS {
  char inputFile[200] ;
  char outFile[200] ;
  char readmeFile[200];

  int NFILT ;
  char filterFile[MXFILTCAL][200] ;
  char filterName[MXFILTCAL][8] ;   // short-hand symbol for filter

  int NSED ;
  char sedFile[MXSEDCAL][200] ;
  char sedName[MXSEDCAL][20] ;  // short-hand name

  // lambda quantities in Angstroms

  int LAMCAL_SIGMA_RANGE[2] ;
  int LAMCAL_SIGMA_BINSIZE ;
  int NBIN_LAMCAL_SIGMA ;

  int LAMCAL_STEP_RANGE[2] ;
  int LAMCAL_STEP_BINSIZE ;
  int NBIN_LAMCAL_STEP ;

  float LAMCAL_OFFSETFRAC_RANGE[2] ;
  float LAMCAL_OFFSETFRAC_BINSIZE ;
  int   NBIN_LAMCAL_OFFSETFRAC ;

  // xxx obsolete  int LAMCAL_MEAN_OFFSET ; // offset from <lambda> of filter

  int LAMINT_BINSIZE ;   // integration bin-size (A)


  // redshift range & binning
  float REDSHIFT_RANGE[2] ;
  float REDSHIFT_BINSIZE ;
  int NBIN_REDSHIFT ;
  float REDSHIFT_MAX; // needed to know how much SED to read

  // filter-dump options
  int  NDUMP_FILTER;
  char DUMP_FILTER_NAME[MXBIN][8];
  int  DUMP_FILTER_LAMSIG[MXBIN];
  int  DUMP_FILTER_LAMSTEP[MXBIN];

  // convert range and binsize into array for instant access
  double ARRAY_LAMSIGMA[MXBIN];
  double ARRAY_LAMSTEP[MXBIN];
  double ARRAY_LAMOFFSETFRAC[MXBIN];
  double ARRAY_REDSHIFT[MXBIN];

} INPUTS ;



// define lambda vs. bin (defined by LAMINT_BINSISE)  
// used for integration. 
// Note that this is common for all filters & SEDs.

double LAMBDA_INTEG[MXLAMINT]; 
double *ptrlam_integ;  // pointer to above
double LAMP_SPECTRUM[MXLAMINT]; 

double MAGREF[MXFILTCAL][MXSEDCAL][MXBIN] ;
double BADMAG;

struct FILTER {
  int    NBIN_LAM ; 
  double LAMBDA_MIN;
  double LAMBDA_MAX;
  double LAMBDA_MEAN;

  double LAMBDA[MXLAM_FILT];
  double TRANS[MXLAM_FILT];

  double TRANS_INTEG[MXLAMINT]; // re-binned for integration

} FILTER[MXFILTCAL] ;



struct SED {
  int    NBIN_LAM ; 
  double LAMBDA[MXLAM_SED];
  double FLUX[MXLAM_SED];

  double FLUX_INTEG[MXLAMINT];

} SED[MXSEDCAL] ;




// ******************************************
int main(int argc, char **argv) {

  int i, ifilt, ised, iz, ioff, isig, istp, NBIN ;
  //  char fnam[] = "main" ;

  double *ptrlam, *ptrfun, *ptrfun_integ ;
  double *ptrflux, *ptrtrans ;
  double z, mag;

  // ---------------- BEGIN ---------------

  sprintf(BANNER,"Begin execution of filtercal_sim.exe " );
  print_banner(BANNER);

  filtercal_init();
  init_SNPATH();

  parse_args(argc, argv );

  read_input();

  proc_input(); // process input

  filtercal_lampspec();

  // --------------------------------------

  // read filter transmissions
  printf("\n");
  for ( i = 1; i <= INPUTS.NFILT ; i++ )
    read_filterTrans(i);

  fflush(stdout);

  // add safety factors onto global lambda ranges used
  // to read SED

  GLOBAL_LAMBDA_MIN -= 2000. ;
  GLOBAL_LAMBDA_MAX += 5000. ;

  // read SEDs
  printf("\n");
  for ( i = 1; i <= INPUTS.NSED ; i++ )
    read_sed(i);


  fflush(stdout);

  // ----------------------------------------------------
  // re-bin everything to use same absolute 
  // lambda bins based on LAMINT_BINSIZE

  printf("\n   Rebin Filter Transmissions");
  for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {
    NBIN          = FILTER[ifilt].NBIN_LAM ;
    ptrlam        = &FILTER[ifilt].LAMBDA[1];
    ptrfun        = &FILTER[ifilt].TRANS[1];
    ptrfun_integ  = FILTER[ifilt].TRANS_INTEG;
    filtercal_rebin(NBIN,ptrlam, ptrfun, ptrfun_integ);
  }
  printf(" and SED Flux. \n");
  for ( ised=1; ised <= INPUTS.NSED; ised++ ) {
    NBIN          = SED[ised].NBIN_LAM ;
    ptrlam        = &SED[ised].LAMBDA[1];
    ptrfun        = &SED[ised].FLUX[1];
    ptrfun_integ  = SED[ised].FLUX_INTEG;
    filtercal_rebin(NBIN, ptrlam, ptrfun, ptrfun_integ);
  }



  // -------------------------------------------------------------------
  printf("\n  Compute reference synthetic AB mag for each SED, Filter, Z \n");

  BADMAG = 99.0;
  for ( ised=1; ised <= INPUTS.NSED; ised++ ) {
    for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {
      for ( iz=1; iz <= INPUTS.NBIN_REDSHIFT; iz++ ) {

	z         = INPUTS.ARRAY_REDSHIFT[iz];
	ptrflux   = SED[ised].FLUX_INTEG ;
	ptrtrans  = FILTER[ifilt].TRANS_INTEG ;
	mag       = SYNMAG( z, ptrlam_integ, ptrtrans, ptrflux );
	MAGREF[ifilt][ised][iz] = mag;
	printf("\t Reference mag(%s-%s, z=%4.2f) = %7.4f \n",
	       INPUTS.sedName[ised], INPUTS.filterName[ifilt], z, mag);
      }
    }
  }


  // change BADMAG after MAGREF is filled so that
  // mag - magref = 99. or -99 if either mag is bad

  BADMAG += 99. ;

  // -------------------------------------
  // open output text file & write header
  open_output();

  printf("\n");

  // -----------------------------------------------------------
  // BIG LOOP over each filter, sed & filter-calib properties

  for ( ised=1; ised <= INPUTS.NSED; ised++ ) {
    for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {

      printf("   Process filtercal for %s-%s  Z = ",
	     INPUTS.sedName[ised], INPUTS.filterName[ifilt] );

      for ( iz=1; iz <= INPUTS.NBIN_REDSHIFT; iz++ ) {

	printf(" %4.2f", INPUTS.ARRAY_REDSHIFT[iz] ); 
	fflush(stdout);

	for ( isig=1; isig <= INPUTS.NBIN_LAMCAL_SIGMA; isig++ ) {
	  for ( istp=1; istp <= INPUTS.NBIN_LAMCAL_STEP;  istp++ ) {
	    for ( ioff=1; ioff <= INPUTS.NBIN_LAMCAL_OFFSETFRAC+1; ioff++ ) {

	      MAGEVAL(ised,ifilt,iz,isig,istp,ioff);

	    }// LAM OFFSET-FRAC
	  }  // LAM STEP
	}   // LAM SIGMA
      }     // REDSHIFT

        printf("\n");

    }       // filter
  }         // SED


  filtercal_end();

  return(0);

} // end of main.


// ******************************************
void parse_args(int argc, char **argv) {
  // ---------- BEGIN --------
  sprintf( INPUTS.inputFile, "%s", argv[1] );
} // end of parse_args



// ****************************
void filtercal_init(void) {

  // init user INPUTS

  INPUTS.NFILT = 0;
  INPUTS.NSED  = 0;

  INPUTS.REDSHIFT_RANGE[0] = 0.0 ;
  INPUTS.REDSHIFT_RANGE[1] = 0.0 ;
  INPUTS.REDSHIFT_BINSIZE  = 0.0 ;

  INPUTS.LAMCAL_SIGMA_RANGE[0] = 0;
  INPUTS.LAMCAL_SIGMA_RANGE[1] = 0;
  INPUTS.LAMCAL_SIGMA_BINSIZE  = 0;

  INPUTS.LAMCAL_STEP_RANGE[0] = 0;
  INPUTS.LAMCAL_STEP_RANGE[1] = 0;
  INPUTS.LAMCAL_STEP_BINSIZE  = 0;


  INPUTS.LAMCAL_OFFSETFRAC_RANGE[0] = 0;
  INPUTS.LAMCAL_OFFSETFRAC_RANGE[1] = 0;
  INPUTS.LAMCAL_OFFSETFRAC_BINSIZE  = 0;

  INPUTS.LAMINT_BINSIZE = 0;

  INPUTS.NDUMP_FILTER = 0;

  sprintf(INPUTS.inputFile, "NULL" );
  sprintf(INPUTS.outFile,   "filtercal_sim.out" );

  // init other global variables.

  ptrlam_integ = LAMBDA_INTEG; // useful local pointer

  GLOBAL_LAMBDA_MIN = +999999. ;
  GLOBAL_LAMBDA_MAX = -999999. ;  

} // 



// ****************************
void filtercal_end(void) {

  fclose(fpout);

  printf("\n GRACEFULL FINISH.  \n");
  printf(" Results written to   %s \n", INPUTS.outFile );
  printf(" Results explained in %s \n", INPUTS.readmeFile );

} // end of filtercal_end

// ****************************
void  read_input(void) {

  FILE *fp;
  char *ptrFile ;
  char fnam[] = "read_input";

  char c_get[60];
  int N;

  // -------- BEGIN --------


  // open user input file.

  ptrFile = INPUTS.inputFile;
  if ( ( fp = fopen(ptrFile, "rt") ) == NULL ) {
    sprintf(c1err, "Cannot find input file:");
    sprintf(c2err,"'%s'", ptrFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n Read input instructions from : %s \n", ptrFile);
  sprintf(c2err,"Check inputs file.");

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"FILTER_TRANS_FILE:") == 0 ) {
      INPUTS.NFILT++ ;
      N = INPUTS.NFILT ;
      if ( N >= MXFILTCAL ) {
	sprintf(c1err,"%d filter files exceeds array bound.", N);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      readchar(fp, INPUTS.filterFile[N] );
      readchar(fp, INPUTS.filterName[N] );
      printf("\t Add filter file : %s (%s) \n"
	     , INPUTS.filterFile[N], INPUTS.filterName[N] ); 
    }

    if ( strcmp(c_get,"SED_FILE:") == 0 ) {
      INPUTS.NSED++ ;
      N = INPUTS.NSED ;
      if ( N >= MXSEDCAL ) {
	sprintf(c1err,"%d SED files exceeds array bound.", N);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      readchar(fp, INPUTS.sedFile[N] );
      readchar(fp, INPUTS.sedName[N] );
      printf("\t Add SED file : %s (%s) \n"
	     , INPUTS.sedFile[N], INPUTS.sedName[N] ); 

    }

    if ( strcmp(c_get,"OUTFILE:") == 0 )
      readchar(fp, INPUTS.outFile );

    if ( strcmp(c_get,"LAMCAL_SIGMA_RANGE:") == 0 )
      readint ( fp, 2, INPUTS.LAMCAL_SIGMA_RANGE );

    if ( strcmp(c_get,"LAMCAL_SIGMA_BINSIZE:") == 0 )
      readint ( fp, 1, &INPUTS.LAMCAL_SIGMA_BINSIZE );

    if ( strcmp(c_get,"LAMCAL_STEP_RANGE:") == 0 )
      readint ( fp, 2, INPUTS.LAMCAL_STEP_RANGE );

    if ( strcmp(c_get,"LAMCAL_STEP_BINSIZE:") == 0 )
      readint ( fp, 1, &INPUTS.LAMCAL_STEP_BINSIZE );

    if ( strcmp(c_get,"LAMCAL_OFFSETFRAC_RANGE:") == 0 )
      readfloat ( fp, 2, INPUTS.LAMCAL_OFFSETFRAC_RANGE );

    if ( strcmp(c_get,"LAMCAL_OFFSETFRAC_BINSIZE:") == 0 )
      readfloat ( fp, 1, &INPUTS.LAMCAL_OFFSETFRAC_BINSIZE );

    if ( strcmp(c_get,"LAMINT_BINSIZE:") == 0 )
      readint ( fp, 1, &INPUTS.LAMINT_BINSIZE );

    if ( strcmp(c_get,"REDSHIFT_RANGE:") == 0 )
      readfloat ( fp, 2, INPUTS.REDSHIFT_RANGE );

    if ( strcmp(c_get,"REDSHIFT_BINSIZE:") == 0 )
      readfloat ( fp, 1, &INPUTS.REDSHIFT_BINSIZE );

    // --
    if ( strcmp(c_get,"DUMP_FILTER:") == 0 ) {
      INPUTS.NDUMP_FILTER++ ;
      N = INPUTS.NDUMP_FILTER ;
      readchar ( fp, INPUTS.DUMP_FILTER_NAME[N] ) ;
      readint  ( fp, 1, &INPUTS.DUMP_FILTER_LAMSIG[N] );
      readint  ( fp, 1, &INPUTS.DUMP_FILTER_LAMSTEP[N] );
    }


  } // end of fscanf

  fclose(fp);

} // end of read_input


// ****************************
void  proc_input(void) {


  int 
    IDIF, NB, ilam, i
    ,ILAM, IBIN, IMIN, IMAX, IVAL
    ;

  double DDIF, DBIN, DMIN, DMAX, di ;

  char fnam[] = "proc_input";

  // ----------- BEGIN -----------


  IMIN = INPUTS.LAMCAL_SIGMA_RANGE[0];
  IMAX = INPUTS.LAMCAL_SIGMA_RANGE[1];
  IDIF = IMAX - IMIN ;
  IBIN = INPUTS.LAMCAL_SIGMA_BINSIZE ;
  if ( IBIN > 0 ) NB   = IDIF / IBIN + 1;
  else            NB   = 1;
  INPUTS.NBIN_LAMCAL_SIGMA = NB;
  if ( NB >= MXBIN ) {
    sprintf(c1err,"NBIN_LAMCAL_SIGMA=%d exceeds array bound of %d", NB,MXBIN);
    sprintf(c2err,"Check LAMCAL_SIGMA_XXX parameters in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  for ( i=1; i <= NB; i++ ) {
    IVAL = IMIN + (i-1)*IBIN ;
    INPUTS.ARRAY_LAMSIGMA[i] = (double)IVAL ;
  }



  IMIN = INPUTS.LAMCAL_STEP_RANGE[0];
  IMAX = INPUTS.LAMCAL_STEP_RANGE[1];
  IDIF = IMAX - IMIN ;
  IBIN = INPUTS.LAMCAL_STEP_BINSIZE ;
  if ( IBIN > 0 )   NB   = IDIF / IBIN + 1 ;
  else              NB   = 1 ;
  INPUTS.NBIN_LAMCAL_STEP = NB  ;
  if ( NB >= MXBIN ) {
    sprintf(c1err,"NBIN_LAMCAL_STEP=%d exceeds array bound of %d", NB,MXBIN);
    sprintf(c2err,"Check LAMCAL_STEP_XXX parameters in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  for ( i=1; i <= NB; i++ ) {
    IVAL = IMIN + (i-1)*IBIN ;
    INPUTS.ARRAY_LAMSTEP[i] = (double)IVAL ;
  }



  DMIN = INPUTS.REDSHIFT_RANGE[0];
  DMAX = INPUTS.REDSHIFT_RANGE[1]; 
  DDIF = DMAX - DMIN ;
  DBIN = INPUTS.REDSHIFT_BINSIZE  ;
  if ( DBIN > 0.0 ) NB = (int)( (DDIF+0.000001) / DBIN ) + 1;
  else              NB = 1 ;
  INPUTS.NBIN_REDSHIFT = NB ;
  if ( NB >= MXBIN ) {
    sprintf(c1err,"NBIN_REDSHIFT=%d exceeds array bound of %d", NB,MXBIN);
    sprintf(c2err,"Check REDSHIFT_XXX parameters in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  for ( i=1; i <= NB; i++ ) {
    di = (double)(i-1) ;
    INPUTS.ARRAY_REDSHIFT[i] = DMIN + di * DBIN;
    INPUTS.REDSHIFT_MAX = INPUTS.ARRAY_REDSHIFT[i] ;
  }



  DMIN = INPUTS.LAMCAL_OFFSETFRAC_RANGE[0];
  DMAX = INPUTS.LAMCAL_OFFSETFRAC_RANGE[1]; 
  DDIF = DMAX - DMIN ;
  DBIN = INPUTS.LAMCAL_OFFSETFRAC_BINSIZE  ;
  if ( DBIN > 0.0 ) NB = (int)( (DDIF+0.000001) / DBIN ) + 1;
  else              NB = 1 ;
  INPUTS.NBIN_LAMCAL_OFFSETFRAC = NB ;
  if ( NB >= MXBIN ) {
    sprintf(c1err,"NBIN_LAMCAL_OFFSETFRAC=%d exceeds array bound of %d", 
	    NB,MXBIN);
    sprintf(c2err,"Check LAMCAL_OFFSETFRAC_XXX parameters in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  for ( i=1; i <= NB; i++ ) {
    di = (double)(i-1) ;
    INPUTS.ARRAY_LAMOFFSETFRAC[i] = DMIN + di * DBIN;
  }

  // add one more offset bin which is really the max magdif
  // over all offsets.

  INPUTS.ARRAY_LAMOFFSETFRAC[NB+1] = 9999. ;


  // ---------------------------------------------
  // setup global lambda bins for integration

  ILAM = 0 ;
  for ( ilam=0; ilam < MXLAMINT; ilam++ ) {
    LAMBDA_INTEG[ilam] = (double)ILAM ;
    ILAM += INPUTS.LAMINT_BINSIZE ;
  }

} //

// ****************************
void read_filterTrans(int ifile) {

  FILE *fp;
  char tmpFile[200];
  double lambda, trans;
  double sum1, sum2;
  int NBIN, ITMP;
  char fnam[] = "read_filterTrans" ;

  // ----------- BEGIN ------------

  // first try to read from current directory;
  // if not there, then read from $SNDATA_ROOT/filters

  sprintf(tmpFile, "%s", INPUTS.filterFile[ifile] );
  if ( ( fp = fopen(tmpFile, "rt") ) != NULL ) 
    goto PRLIB ;

  sprintf(tmpFile, "%s/filters/%s", 
	  PATH_SNDATA_ROOT, INPUTS.filterFile[ifile] );
  if ( ( fp = fopen(tmpFile, "rt") ) == NULL ) {
    sprintf(c1err,"Could not find filter-transmission file");
    sprintf(c2err,"%s", INPUTS.filterFile[ifile] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


 PRLIB:
  printf("   Read filter-trans file: %s \n", tmpFile );

  FILTER[ifile].LAMBDA_MIN  =  999999. ;
  FILTER[ifile].LAMBDA_MAX  = -999999. ;
  FILTER[ifile].NBIN_LAM = 0;

  sum1 = sum2 = 0.0 ;

  while( (fscanf(fp, "%le %le", &lambda, &trans)) != EOF) {

    FILTER[ifile].NBIN_LAM++ ;
    NBIN = FILTER[ifile].NBIN_LAM ;

    if ( NBIN >= MXLAM_FILT ) {
      sprintf(c1err,"%d lambda bins exceeds array bound.", NBIN);
      sprintf(c2err,"Check filter file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    FILTER[ifile].LAMBDA[NBIN] = lambda ;
    FILTER[ifile].TRANS[NBIN]  = trans ;

    sum1 += trans * lambda ;
    sum2 += trans ;

    if ( lambda < FILTER[ifile].LAMBDA_MIN ) 
      FILTER[ifile].LAMBDA_MIN = lambda ;

    if ( lambda > FILTER[ifile].LAMBDA_MAX ) 
      FILTER[ifile].LAMBDA_MAX = lambda ;

  } // end of fscanf 


  // round off MEAN to nearest Angstrom
  ITMP =  (int)(sum1/sum2 + 0.5 ) ;
  FILTER[ifile].LAMBDA_MEAN  = (double)ITMP ;


  printf("\t %4d lambda bins: %5.0f - %5.0f A, <lambda> = %5.0f A \n"
	 , FILTER[ifile].NBIN_LAM
	 , FILTER[ifile].LAMBDA_MIN
	 , FILTER[ifile].LAMBDA_MAX
	 , FILTER[ifile].LAMBDA_MEAN
	 );


  // set global lambda range for reading SEDs,
  // and leave 200 A safety margin

  if ( FILTER[ifile].LAMBDA_MIN < GLOBAL_LAMBDA_MIN ) 
    GLOBAL_LAMBDA_MIN = FILTER[ifile].LAMBDA_MIN  ;

  if ( FILTER[ifile].LAMBDA_MAX > GLOBAL_LAMBDA_MAX ) 
    GLOBAL_LAMBDA_MAX = FILTER[ifile].LAMBDA_MAX  ;


  fclose(fp);

} // end of read_filterTrans


// ****************************
void read_sed(int ifile) {

  FILE *fp;
  char tmpFile[200];
  double lambda, flux;
  double LAM_MIN, LAM_MAX;
  int NBIN;
  char fnam[] = "read_sed" ;

  // ----------- BEGIN ------------

  sprintf(tmpFile, "%s", INPUTS.sedFile[ifile] );

  if ( ( fp = fopen(tmpFile, "rt") ) == NULL ) {
    sprintf(c1err,"Could not find SED file");
    sprintf(c2err,"%s", INPUTS.sedFile[ifile] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  printf("   Read SED file: %s \n", tmpFile );
  SED[ifile].NBIN_LAM = 0;

  LAM_MIN =  999999. ;
  LAM_MAX = -999999. ;


  while( (fscanf(fp, "%le %le", &lambda, &flux)) != EOF) {

    // make max-lam cut, but not min-lam cut to allow for red-shifting
    // if ( lambda < GLOBAL_LAMBDA_MIN ) continue ;
    if ( lambda > GLOBAL_LAMBDA_MAX ) continue ;
    if ( lambda == 0.0              ) continue ;

    if ( lambda > LAM_MAX ) LAM_MAX = lambda ;
    if ( lambda < LAM_MIN ) LAM_MIN = lambda ;

    //    printf(" xxx read SED lambda = %f and flux=%le \n", lambda, flux );

    SED[ifile].NBIN_LAM++ ;
    NBIN = SED[ifile].NBIN_LAM ;

    if ( NBIN >= MXLAM_SED ) {
      sprintf(c1err,"%d lambda bins exceeds array bound.", NBIN);
      sprintf(c2err,"Check SED file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    SED[ifile].LAMBDA[NBIN] = lambda ;
    SED[ifile].FLUX[NBIN]   = flux   ;

  } // end of fscanf 


  printf("\t %4d flux-lambda bins: %5.0f - %5.0f A \n"
	 , SED[ifile].NBIN_LAM, LAM_MIN, LAM_MAX  );


  fclose(fp);

} // end of read_sed


// ****************************
void open_output(void) {

  char *cout, *crdme;

  FILE *fptmp;
  int ifilt, ised;

  // ---------------- BEGIN ---------------

  // open standard output file
  cout  = INPUTS.outFile ;
  fpout = fopen(cout, "wt");

  // open readme file and write relevant info

  crdme = INPUTS.readmeFile ;
  sprintf(crdme,"%s_README", cout );
  fptmp = fopen(crdme, "wt");

  fprintf(fptmp,"Colum definitions: \n");
  fprintf(fptmp,"Colum 1: filter index (see below) \n");
  fprintf(fptmp,"Colum 2: SED index (see below) \n");

  fprintf(fptmp,"Colum 3: redshift \n");
  fprintf(fptmp,"\t (%4.2f - %4.2f, step=%4.2f, %d grid-pts) \n"
	  ,INPUTS.REDSHIFT_RANGE[0]
	  ,INPUTS.REDSHIFT_RANGE[1]
	  ,INPUTS.REDSHIFT_BINSIZE
	  ,INPUTS.NBIN_REDSHIFT
	  );

  fprintf(fptmp,"Colum 4: sigma-width (A) of calib beam \n");
  fprintf(fptmp,"\t (%d - %d A, step=%d A, %d grid-pts) \n"
	  ,INPUTS.LAMCAL_SIGMA_RANGE[0]
	  ,INPUTS.LAMCAL_SIGMA_RANGE[1]
	  ,INPUTS.LAMCAL_SIGMA_BINSIZE
	  ,INPUTS.NBIN_LAMCAL_SIGMA
	  );

  fprintf(fptmp,"Colum 5: step size (A) of calib beam \n");
  fprintf(fptmp,"\t (%d - %d A, step=%d A, %d grid-pts) \n"
	  ,INPUTS.LAMCAL_STEP_RANGE[0]
	  ,INPUTS.LAMCAL_STEP_RANGE[1]
	  ,INPUTS.LAMCAL_STEP_BINSIZE
	  ,INPUTS.NBIN_LAMCAL_STEP
	  );

  fprintf(fptmp,"Colum 6: LAMBDA offset from <lambda> of filter (A) \n");
  fprintf(fptmp,"\t (%4.2f - %4.2f * LAMSTEP,  %d grid-pts) \n"
	  ,INPUTS.LAMCAL_OFFSETFRAC_RANGE[0]
	  ,INPUTS.LAMCAL_OFFSETFRAC_RANGE[1]
	  ,INPUTS.NBIN_LAMCAL_OFFSETFRAC
	  );

  fprintf(fptmp,"Colum 7: average wavelength of measured filter \n");

  fprintf(fptmp,"Colum 8: synthetic mag-shift (measured-true filter) \n" );

  fprintf(fptmp,"\n" );

  for ( ifilt=1; ifilt <= INPUTS.NFILT ; ifilt++ ) {
    fprintf(fptmp,"\t Filter Index = %2.2d => %s   (<lambda> = %6.0f) \n"
	    ,ifilt
	    , INPUTS.filterFile[ifilt] 
	    , FILTER[ifilt].LAMBDA_MEAN
	    );
  }


  for ( ised=1; ised <= INPUTS.NSED ; ised++ ) {
    fprintf(fptmp,"\t SED Index = %2.2d => %s \n", 
	    ised, INPUTS.sedFile[ised] );
  }

  fclose(fptmp);


} // open_output



// ****************************************************************
void filtercal_lampspec(void) {

  // init Lamp spectrum used for calibration

  int ilam;

  // just a flat spectrum until we know more.
  for ( ilam = 0; ilam < MXLAMINT; ilam++ ) {
    LAMP_SPECTRUM[ilam] = 1.0 ;
  }


} // end of filtercal_lampspec

// ****************************************************************
void filtercal_rebin(
		     int NBIN          // (I) number of input bins
		     ,double *ptrlam    // (I) input lambda bins
		     ,double *ptrfun   // (I) input contents (trans or SED)
		     ,double *ptrfun_rebin // (O) rebined contents
		     ) {

  /****
    Start with input array ptrfun vs. ptrlam, and rebin
    using existing LAMBDA_INTEG array.  Output array is
    therefore ptrfun_rebin vs. LAMBDA_INTEG. When all filter
    transmissions and SED are rebined to the same 
    lambda-bins, the integrations and filter-flux convolutions
    will be much easier.

  ***/

  double 
    *ptrlam_rebin
    ,LAM_MIN
    ,LAM_MAX
    ,lam_rebin
    ,lam, fun
    ,lamdif
    ,lamdif_min
    ;

  // use double precision for interpolation function
  double
    lam8, fun8
    ,interp8_lam[4]
    ,interp8_fun[4]
    ;

  int 
    i
    ,ilam
    ,ilam_rebin
    ,ilam_near
    ,ilam_0
    ,ILAM_REBIN_MIN
    ,ILAM_REBIN_MAX
    ;

  int OPT_INTERP = 1 ; // 1,2 => linear, quadratic interp

  // -------- BEGINN ------------
 
  ptrlam_rebin = LAMBDA_INTEG; // useful local pointer

  // get lambda range of input, and translate
  // to re-binned lambda range
  LAM_MIN = *ptrlam ;
  LAM_MAX = *(ptrlam+NBIN-1);
  ILAM_REBIN_MIN = (int)LAM_MIN / INPUTS.LAMINT_BINSIZE - 1 ;
  ILAM_REBIN_MAX = (int)LAM_MAX / INPUTS.LAMINT_BINSIZE + 1 ;

  // init rebined output to zero for all bins

  for ( ilam_rebin = 0; ilam_rebin < MXLAMINT; ilam_rebin++ ) 
    *(ptrfun_rebin + ilam_rebin) = 0.0 ;


  // loop over rebin lambdas that overlap the input contents

  for ( ilam_rebin = ILAM_REBIN_MIN; 
	ilam_rebin < ILAM_REBIN_MAX; ilam_rebin++ ) {

    lam_rebin = *(ptrlam_rebin + ilam_rebin);

    ilam_near = -999;
    lamdif_min = 99999. ;

    // find nearest input lambda bin

    for ( ilam=0; ilam < NBIN; ilam++ ) {
      lam    = *(ptrlam + ilam);
      lamdif = fabsf(lam - lam_rebin) ;
      if ( lamdif < lamdif_min ) 
	{ lamdif_min = lamdif; ilam_near = ilam ; }
    } // ilam

    // load interp arrays centered in ilam_near

    if ( ilam_near == 0 ) 
      ilam_0 = 0 ;
    else if ( ilam_near == NBIN-1 ) 
      ilam_0 = ilam_near - 2 ;
    else
      ilam_0 = ilam_near - 1 ;

    for ( i=1; i <= 3; i++ ) {
      lam   = *(ptrlam + ilam_0 + i-1) ;
      fun   = *(ptrfun + ilam_0 + i-1) ;
      interp8_lam[i]  = (double)lam ;
      interp8_fun[i]  = (double)fun ;
    }
    lam8  = (double)lam_rebin ;

    // interpolate function to re-binned lambda value
    fun8 = interp_1DFUN(OPT_INTERP, lam8, 
			3, &interp8_lam[1], &interp8_fun[1], "rebin" );

    if ( fun8 < 0. )  fun8 = 0.0 ;  // don't allow negative transmission

    // load output array
    *(ptrfun_rebin + ilam_rebin) = fun8;

    // debug-dump 
    if ( ilam_rebin > 1500 && ilam_rebin < -1505 ) {
      printf(" ---------------------------------------------------- \n");
      printf("\n xxx lam_rebin=%6.1f  interpfun=%le \n", lam_rebin, fun8 );
      printf("\n xxx lam[1,2,3]=%6.1f %6.1f %6.1f \n",
	     interp8_lam[1], interp8_lam[2], interp8_lam[3] );
      printf("\n xxx fun[1,2,3]=%le %le %le \n",
	     interp8_fun[1], interp8_fun[2], interp8_fun[3] );
    } 


  }  // ilam_rebin

} // end of filtercal_rebin


// ****************************************************************
double SYNMAG( double z, double *ptrlam, double *ptrtrans, double *ptrflux ) {

  /****
   Returns synthetic AB magnitude for filter transmission (ptrtrans)
   and flux (ptrflux) at redshift 'z'.
   The flux is given in rest-frame in erg/s/cm^2/Hz,
   and transmission in counts:


                   | int [ Flux(lam) * Trans(lam) * lam * dlam ]    |
  mag = -2.5*log10 |----------------------------------------------- |
                   | int [  Trans(lam) * dlam/lam ]                 |


  *****/

  double 
    dlam
    ,lam, lamz
    ,trans
    ,flux_E
    ,flux_count
    ,sumFT
    ,sumTrans
    ,magoff, mag
    ,z1
    ,tmp
    ;

  int ilam, ilamz, NLAM, Lzeroflux  ;

  // ---------- BEGIN -------

  sumFT = sumTrans = 0.0 ;
  dlam = (double)INPUTS.LAMINT_BINSIZE ;
  NLAM = 0 ;
  z1 = 1. + z;

  Lzeroflux = 0;

  for ( ilam=0; ilam < MXLAMINT; ilam++ ) {

    lam     = *(ptrlam   + ilam);
    trans   = *(ptrtrans + ilam);  // trans per count

    if ( trans <= 0.0 ) continue ;

    // get blue rest-frame flux that we see in redshifted (observer) filter
    lamz  = lam / z1;
    tmp   = lamz / dlam + 0.5 ;
    ilamz = (int)tmp ;

    flux_E     = *(ptrflux  + ilamz);  // erg/s/cm^2/A
    flux_count = flux_E * lam/hc;     // per s per cm^2 per A 

    //    printf(" xxxx lam=%6.0f  flux_E = %le ilamz=%d \n", lam, flux_E, ilamz );

    if ( flux_E == 0.0 ) Lzeroflux = 1;

    sumFT    += trans * flux_E * lam  ;
    sumTrans += trans / lam ;

    NLAM++ ;
  }

  sumFT    *= dlam ;
  sumTrans *= dlam ;

  magoff =  2.5 * log10( FNU_AB * LIGHT_A ) ;

  // compute mag if filter trans is non-zero and if
  // SED flux is non-zero in every bin.

  if ( sumTrans > 0.0  && Lzeroflux == 0 ) 
    mag = -2.5 * log10(sumFT/sumTrans) + magoff ;
  else {
    mag = BADMAG ;
  }

  return mag;

} // end of SYNMAG

// ****************************************************************
void MAGEVAL(int ised, int ifilt, int iz, int isig, int istp, int ioff) {

  /****
   Use input indices to extract parameters for measuring
   filter response and measuring synthetic mag.
  ***/


  int       NBIN_FILTCAL    ;

  double 
    Z
    ,OFFSETFRAC
    ,LAM_OFF
    ,LAM_START
    ,LAM_SIGMA
    ,LAM_STEP
    ,LAM_AVG
    ,mag, magdif, absmagdif, magref
    ,FILTCAL_LAMBDA[MXLAM_FILT]
    ,FILTCAL_TRANS[MXLAM_FILT]
    ,FILTCAL_REBIN[MXLAMINT]
    , *ptrlam
    , *ptrtrans
    , *ptrflux 
    ;

  // --------- BEGIN ---------


  Z          = INPUTS.ARRAY_REDSHIFT[iz];
  OFFSETFRAC = INPUTS.ARRAY_LAMOFFSETFRAC[ioff];
  LAM_SIGMA  = INPUTS.ARRAY_LAMSIGMA[isig];
  LAM_STEP   = INPUTS.ARRAY_LAMSTEP[istp];
  LAM_OFF    = OFFSETFRAC * LAM_STEP ;
  LAM_START  = FILTER[ifilt].LAMBDA_MEAN + LAM_OFF ;

  // check special ioff bin that means to store max magdif
  if ( ioff == 1 )  MAGDIF_MAX = 0.0 ;

  if (  ioff == INPUTS.NBIN_LAMCAL_OFFSETFRAC+1 ) {
    LAM_OFF  = OFFSETFRAC ;
    magdif   = MAGDIF_MAX ;
    goto UPDATE_OUTFILE ;
  }

  FILTCAL_CURVE(ifilt, LAM_STEP, LAM_SIGMA, LAM_START           // inputs
		,&NBIN_FILTCAL    // (O) Number of filter bin
		,FILTCAL_LAMBDA   // (O) array of lambdas
		,FILTCAL_TRANS    // (O) array of transmission
		,&LAM_AVG );      // (O) average lambda

  /***
  for ( i=0; i < NBIN_FILTCAL;  i++ ) {
    printf(" xxxxx  FILTER_CURVE point %2d LAM=%6.0f  TRANS=%6.3f \n",
	   i, FILTCAL_LAMBDA[i], FILTCAL_TRANS[i] );
  }
  ***/

  // rebin measured curve to prepare for integration

  filtercal_rebin(NBIN_FILTCAL, FILTCAL_LAMBDA, FILTCAL_TRANS,  // inputs 
		  FILTCAL_REBIN );  // output

  // use measured curve to compute synthetic mag for this SED

  ptrlam   = LAMBDA_INTEG ;
  ptrtrans = FILTCAL_REBIN ;
  ptrflux  = SED[ised].FLUX_INTEG ;
 
  mag    = SYNMAG( Z, ptrlam, ptrtrans, ptrflux );
  magref = MAGREF[ifilt][ised][iz] ;
  magdif = mag - magref ;

  absmagdif = fabs(magdif) ;
  if ( absmagdif > MAGDIF_MAX )  MAGDIF_MAX = absmagdif ;

 UPDATE_OUTFILE:

  fprintf(fpout,"%2d  %2d  %4.2f  %4.0f  %4.0f  %4.0f %6.0f %7.4f \n"
	  ,ifilt, ised, Z, LAM_SIGMA, LAM_STEP, LAM_OFF, LAM_AVG, magdif );


  // check for filter dump
  if ( iz == 1 && ised == 1 && ioff == 1 ) 
    DUMP_FILTER ( ifilt, isig, istp, ptrtrans );


} // end of MAGEVAL


// ************************************************
void  FILTCAL_CURVE(
		    int ifilt            // (I) filter index
		    , double LAM_STEP     // (I) calib step size (A)
		    , double LAM_SIGMA    // (I) calib beam width (A)
		    , double LAM_START    // (I) start at this lambda
		    , int *NBIN_FILTCAL      // (O) number of measures
		    , double *FILTCAL_LAMBDA  // (O) lambda at each measure 
		    , double *FILTCAL_TRANS   // (O) measure trans
		    , double *LAM_AVG         // (O) <lambda>
		    ) {


  /****
   Perform filter-trans measurement based on input lambda-step size
   and beam-width.  Output array of lambdas and measured transmissions.
   Note that LAM_START is near the middle of the filter;
   this defines the offset for the lambda-scan.

  *****/

  double 
    LAM 
    ,LAM_DIF, LAM_ADD
    ,LAM_FIRST, LAM_LAST
    ,TRANS
    ,SUM1, SUM2
    ;

  int N ;

  // -------------------- BEGIN ---------------

  // LAM_ADD accounts for beam size
  LAM_ADD = 3. * LAM_SIGMA ;

  // now figure exact start/end lambdas
  LAM_DIF = LAM_START - FILTER[ifilt].LAMBDA_MIN + LAM_ADD ;
  N = (int)(LAM_DIF/LAM_STEP) + 2 ;
  LAM_FIRST = LAM_START - (double)N * LAM_STEP ;

  LAM_DIF = FILTER[ifilt].LAMBDA_MAX - LAM_START + LAM_ADD ;
  N = (int)(LAM_DIF/LAM_STEP) + 2 ;
  LAM_LAST = LAM_START + (double)N * LAM_STEP ;

  /***
  printf(" xxxx LAM_OFF=%f  LAM_CEN=%f  LAM_DIF=%f LAM_START=%f \n",
	 LAM_OFF,  LAM_CEN,   LAM_DIF,  LAM_START );
  **/

  SUM1 = SUM2 = 0.0 ;
  LAM = LAM_FIRST ;
  N   = 0;
  while ( LAM < LAM_LAST ) {

    LAM += LAM_STEP ;
    TRANS = FILTCAL_POINT(ifilt, LAM, LAM_SIGMA );

    SUM1 += TRANS * LAM ;
    SUM2 += TRANS ;

    // load output array
    *(FILTCAL_LAMBDA + N)  = LAM   ;
    *(FILTCAL_TRANS  + N)  = TRANS ;

    N++;
    //    printf(" xxxx LAM = %f   TRANS = %f \n", LAM, TRANS );

  }

  *NBIN_FILTCAL = N ;

  if ( SUM2 > 0.0 ) 
    *LAM_AVG = SUM1/SUM2;
  else
    *LAM_AVG = 0.0 ;

} // end of FILTCAL_CURVE


// ************************************************
double  FILTCAL_POINT(
		    int ifilt            // (I) filter index
		    ,double LAM_CEN      // (I) central lambda
		    ,double LAM_SIGMA    // (I) Guassian sigma of calib beam
		    ) {

  double 
    XNSIG = 3.0
    ,LAM
    ,LAM_MIN
    ,LAM_MAX  // integration range
    ,LAM_DIF
    ,sumLT   // sum product of lamp * transmission
    ,sumL    // sum of lamp
    ,TRANS_CALIB
    ,TRANS_EXACT
    ,LAMP_SPEC
    ,LAMP_FLUXCAL
    ,SQSIG, ARG
    ;

  int ilam, ILAM_MIN, ILAM_MAX;
  int LDMP;

  //  char fnam[] = "FILTCAL_POINT" ;

  // ---------- BEGIN ------------

  LAM_DIF   = fabs(LAM_CEN-8500.);
  if ( LAM_SIGMA == -99. && LAM_DIF < -5. ) LDMP = 1 ;  else  LDMP = 0 ;

  LAM_MIN = LAM_CEN - XNSIG * LAM_SIGMA ;
  LAM_MAX = LAM_CEN + XNSIG * LAM_SIGMA ;

  ILAM_MIN = (int)LAM_MIN / INPUTS.LAMINT_BINSIZE - 1 ;
  ILAM_MAX = (int)LAM_MAX / INPUTS.LAMINT_BINSIZE + 1 ;

  sumLT = sumL = 0.0;
  SQSIG = LAM_SIGMA * LAM_SIGMA ;
  if ( SQSIG == 0.0 ) SQSIG = 1.0;

  for ( ilam=ILAM_MIN; ilam <= ILAM_MAX; ilam++ ) {
    LAM         = LAMBDA_INTEG[ilam];
    TRANS_EXACT = FILTER[ifilt].TRANS_INTEG[ilam];
    LAMP_SPEC   = LAMP_SPECTRUM[ilam]; 

    LAM_DIF   = (LAM - LAM_CEN);
    ARG       = 0.5*LAM_DIF*LAM_DIF/SQSIG ;
    LAMP_FLUXCAL  = LAMP_SPEC * exp(-ARG);

    if ( LDMP > 0 ) 
      printf("\t LAMPFLUX(lam=%6.0f) = %7.4f   LAMSPEC=%6.3f \n",
	     LAM, LAMP_FLUXCAL, LAMP_SPEC );

    sumLT += LAMP_FLUXCAL * TRANS_EXACT ;
    sumL  += LAMP_FLUXCAL ;
  }


  if ( sumL <= 0.0 ) 
    TRANS_CALIB = 0.0 ;
  else 
    TRANS_CALIB = sumLT / sumL ;


  //  printf(" zzzz SIG=%4.0f LAM=%5.0f INT-range=%5.0f-%5.0f  sumL=%6.2f trans=%5.3f\n",
  //	 LAM_SIGMA, LAM_CEN, LAM_MIN, LAM_MAX, sumL, TRANS_CALIB );


  return TRANS_CALIB ;

} // end of FILTCAL_POINT


// ***********************************************
void  DUMP_FILTER(int ifilt, int isig, int istp, double *TRANS ) {

  int 
    ISIGMA, ISTEP
    ,isigma_tmp
    ,istep_tmp
    ,i, iname
    ,DODUMP
    ,ilam
    ,LTMP1, LTMP2
    ;

  double 
    lam
    ,trans_exact
    ,trans_calib
  ;

  char 
    dmpFile_ref[200]
    ,dmpFile_cal[200]
    ,*ptrName
    ,*ptrName_tmp
    ;

  FILE *fpdmp_ref, *fpdmp_cal;

  // ---------- BEGIN ------------

  ptrName = INPUTS.filterName[ifilt] ;
  ISIGMA  = (int)INPUTS.ARRAY_LAMSIGMA[isig]; // value, not index
  ISTEP   = (int)INPUTS.ARRAY_LAMSTEP[istp];  // value, not index

  // check for match

  DODUMP = 0;
  for ( i=1; i <= INPUTS.NDUMP_FILTER; i++ ) {
    ptrName_tmp = INPUTS.DUMP_FILTER_NAME[i];
    isigma_tmp  = INPUTS.DUMP_FILTER_LAMSIG[i];
    istep_tmp   = INPUTS.DUMP_FILTER_LAMSTEP[i];

    iname = strcmp(ptrName_tmp,ptrName);

    if ( iname == 0 && isigma_tmp == ISIGMA && istep_tmp == ISTEP ) {
      DODUMP = 1;
    }
  }

  if ( DODUMP == 0 ) return;

  // --------------
  // if we get here, do the dump

  sprintf(dmpFile_ref, "filterTrans-%s-SIG0-STEP0.dat",  ptrName );

  sprintf(dmpFile_cal, "filterTrans-%s-SIG%d-STEP%d.dat", 
	  ptrName, ISIGMA, ISTEP );

  printf("\n\t Filter dump: '%s' \n", dmpFile_cal );

  fpdmp_ref = fopen(dmpFile_ref, "wt");
  fpdmp_cal = fopen(dmpFile_cal, "wt");


  for ( ilam=0; ilam < MXLAMINT; ilam++ ) {

    lam = LAMBDA_INTEG[ilam];
    trans_exact = FILTER[ifilt].TRANS_INTEG[ilam];
    trans_calib = *(TRANS+ilam);
    LTMP1 = LTMP2 = 0 ;

    // check if either transmission is non-zero
    if ( trans_exact > 0.0 || trans_calib > 0.0 ) LTMP1 = 1;


    // check lambda range for exact filter trans
    if ( lam >= FILTER[ifilt].LAMBDA_MIN &&
	 lam <= FILTER[ifilt].LAMBDA_MAX ) LTMP2 = 1;

    if ( LTMP1 > 0 || LTMP2 > 0 ) {
      fprintf(fpdmp_ref,"%7.1f  %7.4f\n", lam, trans_exact);
      fprintf(fpdmp_cal,"%7.1f  %7.4f\n", lam, trans_calib);
    }

  }  // end of ilam loop

  fclose(fpdmp_ref);
  fclose(fpdmp_cal);

} // end of DUMP_FILTER
