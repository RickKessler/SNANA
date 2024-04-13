 /*****************************************

  Created June, 2010 by R.Kessler

  Simulate (actually, calculate) magnitude changes as a function
  of changes in the atmospheric extinction as determined  by
  MODTRAN calculations.
  
  Reference (stellar) mags (REFSED_FILE keyword) are used 
  to compute color corrections which are then applied to 
  SN and GALAXY mags as a function of redshift.

  USAGE: 
   atmosphere.exe  <input file>

  Beware that MODTRAN files have wavelength units of nm,
  but all other wavelength units are in Angstroms.

  Also beware that MODTRAN names are hard-coded below
  in the array MODTRAN_TYPE. At some point an upgrade
  is needed to handle arbitray MODTRAN parameter names.


  HISTOGRAMS

  hid               description
  imod           transmission of 'imod' component (imod=1 => total)
  20 + ised      SED flux vs. wavelength for ised
  100            mean wavelength vs. filter index
  100 + ifilt    transmission vs. wavelength
  200 + ifilt    transmission * extinction

  HOFF = 10000 * ichange (ichange = atmosphere-change id)

  HOFF + 1       [1-Trans(ichange)]/[ 1-Trans(default)]

  HOFF + 10 + ifilt     raw Dm vs. SED index (redshift=0)
  HOFF + 20 + ifilt     raw Dm vs. SED index (redshift=0)
  HOFF + 100*ised + 90  raw Dm vs. filt index (redshift=0)

  HOFF + 100*ised + ifilt        raw Dm vs. redshift 
  HOFF + 100*ised + ifilt + 20   cor Dm vs. redshift 



  HISTORY
  ~~~~~~~~~~~~~

  July 9, 2011: add new key MODTRAN_REPLACE to replace any compoment
                with a text-data file.

  Feb 23, 2013: 
    - mask out mkplots_hbook using undefined HBOOK flag.
      To get plots, call functions in sntools_output.c 
      and remove hbook calls in this file.
        
 May 4 2013: 
    - init a few arrays in get_SYNMAG_raw to fix a few valgrig bugs
    - fix bug looping over NREF in  fitColorcor().
    - use SNHIST_INIT() and SNHIST_FILL to work for both hbook and root.
    - rename mkplots_hbook() -> mkplots() 

 Apr 26 2014:
   replace HFILE_OPEN call with wrapper TABLEFILE_OPEN call.
    --> compiles & links, but not tested.

 Apr 2024: fix few compilier bugs

*****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_output.h"

// =======================================
//   declare function prototypes
// =======================================

void  parse_args(int argc, char **argv);
void  read_input(void);
void  read_sedinput(FILE *fp, int NSED, int sedType );
void  read_atmos_change(FILE *fp, int ichange);
void  proc_input(void);
int   INDEX_MODTRAN_TYPE(char *type);

void  read_filterTrans(int ifile);
void  read_sed(int ifile);
void  read_modtran(void);
void  replace_modtran(int ireplace);

void atmoshpere_init(void);
void atmosphere_end(void);

int char2ifilt(char *cfilt); // return integer filter index

void lamrebin(int NBIN, double *ptrlam, double *ptrfun, 
	      double *ptrfun_rebin);

void   rebin_all(void);
void   AtmosCHANGE( int ichange, double *AtmTRANS_CHANGE ); 

double SYNMAG( double z,  double *ptrLam, double *ptrTrans
	       ,double *ptrFlux_rest, double *ptrAtmRatio );

void loopShell(char *copt);

void get_SYNMAG_raw(int ised, int ichange, int ifilt, int iz) ;
void get_SYNMAG_colorcor(int ised, int ichange, int ifilt, int iz) ;

void GET_FILTER_TRANS ( int ifilt, int imod, int ifrac, double *ftrans ) ; 

void get_colorcor(int ichange, int ifilt);
void fitColorcor(int ichange, int NREF, double *mref, double *mxt, 
		 double *color, double *color_slope, double *color_off );

void filter_extinct(void);
void mkplots(void); 

// =======================================
//   declare global variables
// =======================================

#define MXFILTCAL       10    // max # defined filters
#define MXSEDCAL        20    // max # defined SEDs
#define MXLAM_FILT    1500  // max # input lambda bins for filters
#define MXLAM_SED     4000  // max # input lambda bins for SED
#define MXLAM_MODTRAN 3000  // max # input lambda bins for MODTRAN file
#define MXLAMINT     20000  // max # lambda bins for integration
#define MXZBIN        102   // max # bins for Z
#define FNU_AB  3.631E-20 // flat Fnu for AB, erg/cm^2*s*Hz
#define MXCHANGES  10       // max number of atmos changes to define

#define SEDTYPE_REF  1
#define SEDTYPE_SN   2
#define SEDTYPE_GAL  3
char SEDTYPE_STRING[4][8] = { "NULL", "REF", "SN", "GAL" } ;

#define OPT_COLORCOR_RAW  1  // use raw (apparent) color for color cor
#define OPT_COLORCOR_TRUE 2  // use true color for color cor

#define OPT_FILTCOR_DEFAULT 1  // use nominal atmosTrans for all changes
#define OPT_FILTCOR_EXACT   2  // use exact atmosTrans for all changes

#define NMOD_COL        11   // number of MODTRAN extinction columns
#define IMOD_TOTAL       1
#define IMOD_H2Otrans    2
#define IMOD_UMIX        3
#define IMOD_03          4
#define IMOD_TRACE       5
#define IMOD_N2          6 
#define IMOD_H2Ocont     7
#define IMOD_MOLEC       8
#define IMOD_AEROtrans   9
#define IMOD_HNO3       10
#define IMOD_AEROabtrns 11

#define OPT_INTERP  1    // 1,2 => linear, quadratic interp

char MODTRAN_TYPE[NMOD_COL+1][20] = {
  "NULL", "TOTAL", "H2Otrans", "UMIX", "O3", "TRACE", "N2",
  "H2Ocont", "MOLEC", "AEROtrans", "HNO3", "AEROabtrns" } ;


double GLOBAL_LAMBDA_MIN; //  global lambda min,max among filters;
double GLOBAL_LAMBDA_MAX; //  used to limit SED range that is read

double ONEARRAY[MXLAMINT+1] ;

struct {
  char inputFile[200] ;
  char hbook_outFile[200] ;
  char readmeFile[200];

  int  NFILT ;
  char filterDir[200] ; // optional filter path
  char filterFile[MXFILTCAL][200] ;
  char filterName[MXFILTCAL][8] ;      // short-hand symbol for filter
  char filterColorcorString[MXFILTCAL][4];
  int  ifilterColorcor[MXFILTCAL][2];  // ifilt's used for color cor

  int NSED ;
  char sedDir[200];
  char sedFile[MXSEDCAL][200] ;
  char sedName[MXSEDCAL][20] ;  // short-hand name
  int  sedType[MXSEDCAL];       // ref, SN or galaxy
 
  char ModtranFile[200] ; // extinction file from MODTRAN
  char ModtranReplace[NMOD_COL+1][2][200] ; // replace component with file
  int  NREPLACE_Modtran ;

  int FILTER_MODTRAN_FLAG ;

  int NCHANGES ; // number of changes to try
  float ATMOS_CHANGE[MXCHANGES+1][NMOD_COL+1];

  // lambda quantities in Angstroms
  int LAMINT_BINSIZE ;   // integration bin-size (A)

  // redshift range & binning
  float REDSHIFT_RANGE[2] ;
  float REDSHIFT_BINSIZE ;
  int   NBIN_REDSHIFT ;
  float REDSHIFT_MAX; // needed to know how much SED to read
  double ARRAY_REDSHIFT[MXZBIN]; // calculate array

  int OPT_COLORCOR; // see OPT_COLORCOR_XXX parameters
  int OPT_FILTCOR;  // which atmosTrans to use

} INPUTS ;



// define lambda vs. bin (defined by LAMINT_BINSISE)  
// used for integration. 
// Note that this is common for all filters & SEDs.

double LAMBDA_REBIN[MXLAMINT+1]; 
double *ptrlam_rebin;  // pointer to above

double BADMAG;

struct {
  int    NBIN_LAM ; 
  double LAMBDA_MIN;
  double LAMBDA_MAX;
  double LAMBDA_MEAN;

  double LAMBDA[MXLAM_FILT];
  double TRANS[MXLAM_FILT];

  double TRANS_REBIN[MXLAMINT+1];   // input filter trans
  double TRANSTOT_REBIN[MXLAMINT+1]; // filter trans *= atmTrans

  double color_slope[MXCHANGES] ;
  double color_offset[MXCHANGES] ;

} FILTER[MXFILTCAL] ;



struct {
  int    NBIN_LAM ; 
  double LAMBDA[MXLAM_SED];
  double FLUX[MXLAM_SED];
  double FLUX_REBIN[MXLAMINT+1] ;

} SED[MXSEDCAL] ;


struct {
  int    NBIN_LAM;
  double LAMBDA[MXLAM_MODTRAN];

  // define atm trans vs. component (IMOD_TOTAL = net trans)
  double AtmTRANS[NMOD_COL+1][MXLAM_MODTRAN];  // from MODTRAN file
  double AtmTRANS_REBIN[NMOD_COL+1][MXLAMINT+1]; // rebinned

  // total atmos trans vs. user-requested changes
  double AtmTRANS_CHANGE[MXCHANGES][MXLAMINT+1] ;
} MODTRAN ;


// nominal mags with/without extinction.
// 'avecXT' is that reference for the changes below
double SYNMAG_sansXT[MXFILTCAL][MXSEDCAL][MXZBIN] ; 
double SYNMAG_avecXT[MXFILTCAL][MXSEDCAL][MXZBIN] ; // filter includes atmos

// define mag changes from atmosph using nominal filter trans;
// effectively the SED changes based on the atmos changes
double SYNMAG_change_raw[MXCHANGES+1][MXFILTCAL][MXSEDCAL][MXZBIN] ; 
double SYNMAG_change_colorcor[MXCHANGES+1][MXFILTCAL][MXSEDCAL][MXZBIN] ; 

int FLAGTMP;

// ******************************************
int main(int argc, char **argv) {

  int i, ifilt, ichange;
  double *ptrchange ;
  //  char fnam[] = "main" ;
  // ---------------- BEGIN ---------------

  sprintf(BANNER,"Begin execution of atmosphere.exe " );
  print_banner(BANNER);

  atmoshpere_init();
  init_SNPATH();

  parse_args(argc, argv );

  read_input();

  proc_input(); // process input

  // --------------------------------------

  // read filter transmissions
  printf("\n");
  for ( i = 1; i <= INPUTS.NFILT ; i++ )
    { read_filterTrans(i); }


  GLOBAL_LAMBDA_MIN -= 2. * INPUTS.LAMINT_BINSIZE ;
  GLOBAL_LAMBDA_MAX += 2. * INPUTS.LAMINT_BINSIZE ;
  printf("\n\t GLOBAL_LAMBDA_[MIN,MAX] = %6.1f,%6.1f \n",
	 GLOBAL_LAMBDA_MIN, GLOBAL_LAMBDA_MAX );
  
  fflush(stdout);

  // read SEDs
  printf("\n   Read SEDs ... \n");
  for ( i = 1; i <= INPUTS.NSED ; i++ )
    { read_sed(i); }


  fflush(stdout);

  read_modtran();

  for ( i=1; i <= INPUTS.NREPLACE_Modtran; i++ ) 
    { replace_modtran(i); }

  // ----------------------------------------------------
  // re-bin everything to use same absolute 
  // lambda bins based on LAMINT_BINSIZE

  rebin_all();

  // ------ apply nominal extinction to filter transmissions
  filter_extinct();
  fflush(stdout);


  // compute total atmos trans for each user change
  for ( ichange = 1; ichange <= INPUTS.NCHANGES; ichange++ ) {     
    ptrchange = MODTRAN.AtmTRANS_CHANGE[ichange];
    AtmosCHANGE( ichange, ptrchange );    // returns ptrchange
    fflush(stdout);
  }

  // -------------------------------------------------------------------

  BADMAG = 99.0;

  printf("\n  Compute RAW synthetic mags ... \n");
  loopShell("raw");  // get apparent 'raw' mags


  // get color corrections for each 'ichange'
  for ( ichange=1; ichange <= INPUTS.NCHANGES ; ichange++ ) {
    for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {
      get_colorcor(ichange,ifilt);
    }
  }

  printf("\n  Compute Color-Corrected synthetic mags ... \n");
  loopShell("colorcor");  // get apparent 'raw' mags


  mkplots();

  atmosphere_end();

  return(0);

} // end of main.


// ******************************************
void parse_args(int argc, char **argv) {
  // ---------- BEGIN --------
  sprintf( INPUTS.inputFile, "%s", argv[1] );
} // end of parse_args



// ****************************
void atmoshpere_init(void) {

  int i, imod;

  // init user INPUTS

  INPUTS.NFILT = 0;
  INPUTS.NSED  = 0;

  INPUTS.filterDir[0] = 0 ;
  INPUTS.sedDir[0]    = 0 ; 

  INPUTS.REDSHIFT_RANGE[0] = 0.0 ;
  INPUTS.REDSHIFT_RANGE[1] = 0.0 ;
  INPUTS.REDSHIFT_BINSIZE  = 0.0 ;

  INPUTS.LAMINT_BINSIZE = 0;

  INPUTS.OPT_COLORCOR = -1 ;
  INPUTS.OPT_FILTCOR  = -1 ;

  INPUTS.NREPLACE_Modtran = 0;

  sprintf(INPUTS.inputFile,      "NULL" );
  sprintf(INPUTS.hbook_outFile,  "atmosphere.his" );

  ptrlam_rebin = LAMBDA_REBIN; // useful local pointer

  INPUTS.FILTER_MODTRAN_FLAG = 0 ; // default: Atmos XT not part of filterTrans

  INPUTS.NCHANGES = 0 ;
  for ( i=0; i<= MXCHANGES; i++ ) {
    for ( imod=0 ; imod <= NMOD_COL; imod++ ) {
      INPUTS.ATMOS_CHANGE[i][imod] = 1.0 ; // default is no change
    }
  }

  FLAGTMP = 0;

  GLOBAL_LAMBDA_MIN =  9999999. ;
  GLOBAL_LAMBDA_MAX = -9999999. ;

}   // atmosphere_init



// ****************************
void atmosphere_end(void) {


  printf(" GRACEFULL FINISH.  \n");
  //  printf(" Results explained in %s \n", INPUTS.readmeFile );

} // end of atmosphere_end

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

    if ( strcmp(c_get,"FILTER_DIR:") == 0 )
      { readchar(fp, INPUTS.filterDir ); }

    if ( strcmp(c_get,"FILTER_TRANS_FILE:") == 0 ) {
      INPUTS.NFILT++ ;
      N = INPUTS.NFILT ;
      if ( N >= MXFILTCAL ) {
	sprintf(c1err,"%d filter files exceeds array bound.", N);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      readchar(fp, INPUTS.filterFile[N] );
      readchar(fp, INPUTS.filterName[N] );
      readchar(fp, INPUTS.filterColorcorString[N] );

      printf("\t Add filter file : %s (%s) \n"
	     , INPUTS.filterFile[N], INPUTS.filterName[N] ); 
    }

    // ----------
    if ( strcmp(c_get,"SED_DIR:") == 0 )
      {  readchar(fp, INPUTS.sedDir ); }

    if ( strcmp(c_get,"REFSED_FILE:") == 0 ) {
      INPUTS.NSED++ ;
      read_sedinput(fp, INPUTS.NSED, SEDTYPE_REF );
    }
    if ( strcmp(c_get,"SNSED_FILE:") == 0 ) {
      INPUTS.NSED++ ;
      read_sedinput(fp, INPUTS.NSED, SEDTYPE_SN );
    }
    if ( strcmp(c_get,"GALSED_FILE:") == 0 ) {
      INPUTS.NSED++ ;
      read_sedinput(fp, INPUTS.NSED, SEDTYPE_GAL );
    }
    // --------


    if ( strcmp(c_get,"HBOOK_OUTFILE:") == 0 )
      readchar(fp, INPUTS.hbook_outFile );

    if ( strcmp(c_get,"FILTER_MODTRAN_FLAG:") == 0 )
      readint(fp, 1, &INPUTS.FILTER_MODTRAN_FLAG );

    if ( strcmp(c_get,"REDSHIFT_RANGE:") == 0 )
      readfloat ( fp, 2, INPUTS.REDSHIFT_RANGE );

    if ( strcmp(c_get,"REDSHIFT_BINSIZE:") == 0 )
      readfloat ( fp, 1, &INPUTS.REDSHIFT_BINSIZE );

    if ( strcmp(c_get,"LAMINT_BINSIZE:") == 0 )
      readint ( fp, 1, &INPUTS.LAMINT_BINSIZE );

    if ( strcmp(c_get,"OPT_COLORCOR:") == 0 ) {
      readint ( fp, 1, &INPUTS.OPT_COLORCOR );
      printf("\t OPT_COLORCOR: %d \n", INPUTS.OPT_COLORCOR );
    }

    if ( strcmp(c_get,"OPT_FILTCOR:") == 0 ) {
      readint ( fp, 1, &INPUTS.OPT_FILTCOR );
      printf("\t OPT_FILTCOR: %d \n", INPUTS.OPT_FILTCOR );
    }

    // ---- atmos stuff

    if ( strcmp(c_get,"MODTRAN_FILE:") == 0 )
      {  readchar(fp, INPUTS.ModtranFile ); }

    if ( strcmp(c_get,"MODTRAN_REPLACE:") == 0 ) {
      INPUTS.NREPLACE_Modtran++ ;
      N = INPUTS.NREPLACE_Modtran ;
      readchar(fp, INPUTS.ModtranReplace[N][0] ); // component
      readchar(fp, INPUTS.ModtranReplace[N][1] ); // file with Trans
    }

    if ( strcmp(c_get,"@ATMOS_CHANGE") == 0 ) {
      INPUTS.NCHANGES++ ;
      read_atmos_change(fp, INPUTS.NCHANGES );
    }

  } // end of fscanf

  fclose(fp);


} // end of read_input


// ********************************************
void read_sedinput(FILE *fp, int NSED, int sedType ) {

  char fnam[] = "read_sedinput";


  if ( NSED >= MXSEDCAL ) {
    sprintf(c1err,"%d SED files exceeds array bound.", NSED);
    sprintf(c2err,"Check inputs file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  readchar(fp, INPUTS.sedFile[NSED] );
  readchar(fp, INPUTS.sedName[NSED] );
  INPUTS.sedType[NSED] = sedType;

  printf("\t Add %3s-SED(%2d) file : %s (%s) \n"
	 ,SEDTYPE_STRING[sedType], NSED
	 ,INPUTS.sedFile[NSED]
	 ,INPUTS.sedName[NSED]   ); 

}

// *****************************************
void read_atmos_change(FILE *fp, int ichange) {

  int   imod ;
  float xfac;
  char  ctype[40];
  //  char  fnam[] = "read_atmos_change" ;

  // -------- BEGIN ----------

 RDMODTYPE:
  readchar(fp, ctype );
  if ( strcmp(ctype,"@DONE") == 0 ) { return ; }
  readfloat(fp, 1, &xfac );
  imod = INDEX_MODTRAN_TYPE(ctype); // translate string to integer index

  INPUTS.ATMOS_CHANGE[ichange][imod] = xfac ;

  printf("\t ATMOS_CHANGE[%d] : %5.2f x %s \n", ichange, xfac, ctype );

  goto RDMODTYPE;

} // end of read_atmos_change


// ****************************
int INDEX_MODTRAN_TYPE(char *type) {

  int imod ;
  char fnam[] = "INDEX_MODTRAN_TYPE" ;

  // -------- BEGIN -----------

  for ( imod=1; imod <= NMOD_COL; imod++ ) {
    if ( strcmp(MODTRAN_TYPE[imod],type) == 0 ) { 
      return imod ; 
    }
  }


  // if we get here, then abort on error

  sprintf(c1err,"Invalid MODTRAN TYPE: '%s'", type );
  sprintf(c2err,"Check  char MODTRAN_TYPE  for valid strings");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 

  return(-9);

} // end of INDEX_MODTRAN_TYPE


// ****************************
void  proc_input(void) {

  int  NB, ilam, i, ILAM, ifilt, if0, if1    ;
  double DDIF, DBIN, DMIN, DMAX, di ;
  char string[4], cf0[2], cf1[2];
  char fnam[] = "proc_input";

  // ----------- BEGIN -----------

  printf(" \n");

  DMIN = INPUTS.REDSHIFT_RANGE[0];
  DMAX = INPUTS.REDSHIFT_RANGE[1]; 
  DDIF = DMAX - DMIN ;
  DBIN = INPUTS.REDSHIFT_BINSIZE  ;
  if ( DBIN > 0.0 ) NB = (int)( (DDIF+0.000001) / DBIN ) + 1;
  else              NB = 1 ;
  INPUTS.NBIN_REDSHIFT = NB ;
  if ( NB >= MXZBIN ) {
    sprintf(c1err,"NBIN_REDSHIFT=%d exceeds array bound of %d", NB,MXZBIN);
    sprintf(c2err,"Check REDSHIFT_XXX parameters in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  for ( i=1; i <= NB; i++ ) {
    di = (double)(i-1) ;
    INPUTS.ARRAY_REDSHIFT[i] = DMIN + di * DBIN;
    INPUTS.REDSHIFT_MAX = INPUTS.ARRAY_REDSHIFT[i] ;
  }

  // ---------------------------------------------
  // setup global lambda bins for integration

  ILAM = 0 ;
  for ( ilam=0; ilam < MXLAMINT; ilam++ ) {
    LAMBDA_REBIN[ilam] = (double)ILAM ;
    ILAM += INPUTS.LAMINT_BINSIZE ;
    ONEARRAY[ilam] = 1.0 ;
  }


  // decode filterColorcorString into integer ifilt indices
  // => used for color correction


  for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {
    sprintf(string, "%s", INPUTS.filterColorcorString[ifilt] );
    sprintf(cf0, "%c", string[0] );
    sprintf(cf1, "%c", string[2] );

    if0 = char2ifilt(cf0);
    if1 = char2ifilt(cf1);

    INPUTS.ifilterColorcor[ifilt][0] = if0 ;
    INPUTS.ifilterColorcor[ifilt][1] = if1 ;

    /*
    printf(" xxx %s colorcor = %s =>  '%s' and '%s' => IFILT = %d - %d \n",
	   INPUTS.filterName[ifilt], string, cf0, cf1, if0, if1 );
    */

  }


} // proc_input


// **********************************
int char2ifilt ( char *cfilt ) {

  // Return integer ifilt index for char string cfilt

  int ifilt;
  char fnam[] = "char2ifilt" ;


  for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {
    if ( strcmp(INPUTS.filterName[ifilt],cfilt) == 0 ) {
      return ifilt;
    }
  }

  sprintf(c1err,"Could not find filter index for '%s'", cfilt ) ;
  errmsg(SEV_FATAL, 0, fnam, c1err, ""); 

  return(-9);

} // end of char2ifilt

// ****************************
void read_filterTrans(int ifile) {

  FILE *fp;
  double lambda, trans;
  double sum1, sum2;
  int NBIN, ITMP, gzipFlag ;

  char 
    PATH_SNANA_FILTERS[200]
    ,filterFile[200], *fptr
    ,tmpFile[200]
    ,fnam[] = "read_filterTrans" 
    ;

  // ----------- BEGIN ------------

  // first try to read from current directory;
  // if not there, then read from $SNDATA_ROOT/filters

  fptr = INPUTS.filterFile[ifile] ;
  if ( strlen(INPUTS.filterDir) == 0 ) 
    { sprintf(filterFile,"%s", fptr ); }
  else
    { sprintf(filterFile,"%s/%s", INPUTS.filterDir, fptr ); }


  sprintf(PATH_SNANA_FILTERS,"%s/filters", PATH_SNDATA_ROOT );
  fp = snana_openTextFile(0, PATH_SNANA_FILTERS, filterFile, 
			  tmpFile, &gzipFlag );

  if ( fp == NULL ) {
    sprintf(c1err,"Could not find filter-transmission file");
    sprintf(c2err,"%s", INPUTS.filterFile[ifile] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


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
      { FILTER[ifile].LAMBDA_MIN = lambda ; }

    if ( lambda > FILTER[ifile].LAMBDA_MAX ) 
      { FILTER[ifile].LAMBDA_MAX = lambda ; }

  } // end of fscanf 


  // round off MEAN to nearest Angstrom

  if ( sum2 <= 0.0 ) {
    sprintf(c1err,"trans-sum(sum2) = 0 for filter ifile = %d", ifile);
    sprintf(c2err,"File: %s", fptr);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    
  }

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
    { GLOBAL_LAMBDA_MIN = FILTER[ifile].LAMBDA_MIN  ; }

  if ( FILTER[ifile].LAMBDA_MAX > GLOBAL_LAMBDA_MAX ) 
    { GLOBAL_LAMBDA_MAX = FILTER[ifile].LAMBDA_MAX  ; }


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

  sprintf(tmpFile, "%s/%s", INPUTS.sedDir, INPUTS.sedFile[ifile] );

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


  printf("\t %4d flux-lambda bins: %5.0f - %5.0f A   (%s) \n"
	 , SED[ifile].NBIN_LAM, LAM_MIN, LAM_MAX, INPUTS.sedName[ifile]  );


  fclose(fp);

} // end of read_sed


// ****************************************************************
void read_modtran(void) {

  FILE *fp ;

  int icol, ilam ;
  double lam, trans ;

  char *ptrFile;
  char c_get[40];
  char fnam[] = "read_modtran";

  // ------------ BEGIN ------------

  ptrFile = INPUTS.ModtranFile ;

  if ( ( fp = fopen(ptrFile, "rt") ) == NULL ) {
    sprintf(c1err,"Could not find MODTRAN file");
    sprintf(c2err,"%s", ptrFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  printf("\n Read atmospheric extinction parameters from: \n");
  printf("  %s \n", ptrFile );

  ilam = 0;
  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"LAM:") == 0 ) {

      ilam++ ;

      readdouble ( fp, 1, &lam ) ;
      MODTRAN.LAMBDA[ilam] = lam * 10. ; // convert nm to A

      for ( icol = 1; icol <= NMOD_COL; icol++ ) {
	readdouble ( fp, 1, &trans );
	MODTRAN.AtmTRANS[icol][ilam] = trans ;

	if ( icol == -2 && trans < 1.0 ) 
	  printf("\t xxx LAM=%f TRANS(2)=%f \n", lam, trans );
 
      } // icol
    }

  } 

  if ( ilam >= MXLAM_MODTRAN ) {
    sprintf(c1err,"%d LAM entries exceeds MODTRAN array bound", ilam );
    sprintf(c2err,"of MXLAM_MODTRAN =  %d ", MXLAM_MODTRAN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  MODTRAN.NBIN_LAM = ilam;
  printf("\t Read %d lambda bins. \n\n", MODTRAN.NBIN_LAM );

  fclose(fp);

} // end of read_modtran



// ****************************************************************
void replace_modtran(int ireplace ) {

  int imod, Nrow, ilam ;

  double 
    lam
    ,Trans_orig
    ,Trans_replace
    ,TMP_LAM[MXLAM_MODTRAN] 
    ,TMP_TRANS[MXLAM_MODTRAN]
    ;

  char
    component[60]
    ,transFile[200]
    ,fnam[] = "replace_modtran" 
    ;

  // -------------- BEGIN --------------

  sprintf(component, "%s", INPUTS.ModtranReplace[ireplace][0] );
  sprintf(transFile, "%s", INPUTS.ModtranReplace[ireplace][1] );

  imod = INDEX_MODTRAN_TYPE(component);

  printf(" Replace modtran component %s(%d) with %s \n", 
	 component, imod, transFile );

  rd2columnFile(transFile, MXLAM_MODTRAN, 
		&Nrow, 	TMP_LAM, TMP_TRANS, 0 ); // returned

  // convert lambda from nm to A
  for ( ilam=0 ; ilam < Nrow; ilam++ ) {
    TMP_LAM[ilam] *= 10.0 ;
  }

  // now translate the TMP_TRANS into the lambda grid of
  // the original MODTRAN grid.

  for ( ilam=1; ilam <= MODTRAN.NBIN_LAM ; ilam++ ) {
    lam   = MODTRAN.LAMBDA[ilam] ;
    Trans_replace = 
      interp_1DFUN(OPT_INTERP, lam, Nrow, TMP_LAM, TMP_TRANS, fnam);
    Trans_orig = MODTRAN.AtmTRANS[imod][ilam] ;
    MODTRAN.AtmTRANS[imod][ilam] = Trans_replace ;

    // update total transmission by multiplying by Trans_replace/Trans_orig
    MODTRAN.AtmTRANS[1][ilam] *= Trans_replace/Trans_orig ;
  }

  //  debugexit("replace modtrans");

} // end of  replace_modtran

// ****************************************************************
void loopShell(char *copt) {

  int ised, ichange, ifilt, iz, NZ;

  char fnam[] = "loopShell";

  // ------------ BEGIN -----------

  for ( ichange=1; ichange <= INPUTS.NCHANGES; ichange++ ) {
    for ( ised=1; ised <= INPUTS.NSED; ised++ ) {
 
      if ( INPUTS.sedType[ised] == SEDTYPE_REF ) 
	NZ = 1;
      else
	NZ = INPUTS.NBIN_REDSHIFT ;

      if ( ichange == 1 ) {
	printf("\t Process SED %s  (NZBIN=%d) \n", INPUTS.sedName[ised], NZ );
	fflush(stdout) ;
      }

      for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {
	for ( iz=1; iz <= NZ ; iz++ ) {

	  if ( strcmp(copt,"raw") == 0 ) 
	    get_SYNMAG_raw(ised,ichange,ifilt,iz); 
	  else if ( strcmp(copt,"colorcor") == 0 ) 
	    get_SYNMAG_colorcor(ised,ichange,ifilt,iz); 
	  else {
	    sprintf(c1err,"Invalid option '%s' ", copt);
	    errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
	  }

	} // iz
      } // ifilt

    } // ised
  } // ichange

} // end of loopShell

// ****************************************************************
void get_SYNMAG_raw(int ised, int ichange, int ifilt, int iz) {

  double 
    *ptrFlux, *ptrTrans
    ,z, mag, magxt, trans
    ,TRnominal, TRchange,  lam
    ,TR_ratio
    ,atmRatio[MXLAMINT]
    ,tmpTrans[MXLAMINT]
    ;

  int ilam ;

  char fnam[] = "get_SYNMAG_raw" ;

  // ------------ BEGIN ------------

  // first evaluate mags with/without atmos extinction
  // just to see the net effect.

  z         = INPUTS.ARRAY_REDSHIFT[iz];
  ptrFlux   = SED[ised].FLUX_REBIN ;  // rest-frame SED

  // mag with no atmos extinction
  ptrTrans  = FILTER[ifilt].TRANS_REBIN ;
  mag       = SYNMAG( z,  ptrlam_rebin, ptrTrans, ptrFlux, ONEARRAY );

  // mag with nominal atmos extinction
  ptrTrans  = FILTER[ifilt].TRANSTOT_REBIN ;
  magxt     = SYNMAG( z, ptrlam_rebin, ptrTrans, ptrFlux, ONEARRAY );

  SYNMAG_sansXT[ifilt][ised][iz] = mag;    // without atmos extinction
  SYNMAG_avecXT[ifilt][ised][iz] = magxt;  // with atmos extinction

  /*
  printf("\t mag(%s-%s, z=%4.2f) = %7.4f  -> %7.4f (with AtmosXtinc)\n",
	 INPUTS.sedName[ised], INPUTS.filterName[ifilt], z, mag, magxt );
  */

  // --------------

  for ( ilam=0; ilam < MXLAMINT; ilam++ ) {     
    lam = LAMBDA_REBIN[ilam] ;

    tmpTrans[ilam] =  0.0 ; // init entire array (May 2013)
    atmRatio[ilam] =  0.0 ;

    if ( lam < GLOBAL_LAMBDA_MIN ) continue ;
    if ( lam > GLOBAL_LAMBDA_MAX ) continue ;

    trans     = FILTER[ifilt].TRANSTOT_REBIN[ilam] ;
    TRnominal = MODTRAN.AtmTRANS_REBIN[IMOD_TOTAL][ilam];
    TRchange  = MODTRAN.AtmTRANS_CHANGE[ichange][ilam]; 

    if ( TRnominal > 0.0 ) 
      { TR_ratio  = TRchange/TRnominal ; }
    else
      { TR_ratio = 0.0 ; }

    if ( INPUTS.OPT_FILTCOR == OPT_FILTCOR_DEFAULT ) {
      atmRatio[ilam] = TR_ratio ; // modify atmoTrans for SED
      tmpTrans[ilam] = trans ;    // use default atmos Trans for filter
    }
    else if ( INPUTS.OPT_FILTCOR == OPT_FILTCOR_EXACT ) {
      atmRatio[ilam] = ONEARRAY[ilam];    // do NOT modify SED
      tmpTrans[ilam] = trans * TR_ratio ; // correct filter trans
    }
    else {
      sprintf(c1err,"Invalid OPT_FILTCOR = %d", INPUTS.OPT_FILTCOR );
      errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
    }
  }   // ilam 
  
  magxt = SYNMAG( z, ptrlam_rebin, tmpTrans, ptrFlux, atmRatio ) ;
  SYNMAG_change_raw[ichange][ifilt][ised][iz] = magxt ;



} // end of get_SYNMAG_raw


// ****************************************************************
void get_SYNMAG_colorcor(int ised, int ichange, int ifilt, int iz) {

  // Apply color correction for the above set of indices

  int if0, if1, OPT;

  double 
    slope, off
    ,mag_raw, mag0=0.0, mag1=0.0, color
    ,mag_cor
    ;

  // ---------- BEGIN -----------

  slope = FILTER[ifilt].color_slope[ichange] ;
  off   = FILTER[ifilt].color_offset[ichange] ;

  if0 = INPUTS.ifilterColorcor[ifilt][0];
  if1 = INPUTS.ifilterColorcor[ifilt][1];

  OPT = INPUTS.OPT_COLORCOR;

  mag_raw   = SYNMAG_change_raw[ichange][ifilt][ised][iz];

  if ( OPT == OPT_COLORCOR_RAW  ){
    mag0  = SYNMAG_change_raw[ichange][if0][ised][iz];
    mag1  = SYNMAG_change_raw[ichange][if1][ised][iz];
  }
  else if ( OPT == OPT_COLORCOR_TRUE  ){
    mag0 = SYNMAG_avecXT[if0][ised][iz] ; // mag with nominal extintion
    mag1 = SYNMAG_avecXT[if1][ised][iz] ;
  }

  color = mag1 - mag0 ;

  mag_cor = mag_raw + slope * ( color - off ) ;

  SYNMAG_change_colorcor[ichange][ifilt][ised][iz] = mag_cor ;

} // end of get_SYNMAG_colorcor


// ****************************************************************
double SYNMAG( double z,  double *ptrLam, double *ptrTrans, 
	       double *ptrFlux_rest, double *ptrAtmRatio ) {

  /****
   Returns synthetic AB magnitude for filter transmission (ptrtrans)
   and rest-frame flux (ptrflux_rest) at redshift 'z'.
   The flux is given in rest-frame in erg/s/cm^2/Hz,
   and transmission in counts:


                   | int [ Flux(lam) * Trans(lam) * lam * dlam ]    |
  mag = -2.5*log10 |----------------------------------------------- |
                   | int [  Trans(lam) * dlam/lam ]                 |


  May 4 2013: check Lzeroflux BEFORE multiplying atmos trans
              because filter leakage can extend the lambda range
              to where atmosTrans = 0.

  *****/

  double 
    dlam
    ,lam, lamz
    ,trans
    ,atmRatio
    ,flux_E
    ,sumFT
    ,sumTrans
    ,magoff, mag
    ,z1
    ,tmp
    ;

  int ilam, ilamz, Lzeroflux  ;

  // ---------- BEGIN -------

  sumFT = sumTrans = 0.0 ;
  dlam = (double)INPUTS.LAMINT_BINSIZE ;
  z1 = 1. + z ;

  Lzeroflux = 0;


  for ( ilam=0; ilam < MXLAMINT; ilam++ ) {

    lam      = *(ptrLam   + ilam);
    trans    = *(ptrTrans + ilam);     // trans per count
    atmRatio = *(ptrAtmRatio + ilam);  // ratio of atmosTrans to nominal
  

    if ( trans <= 0.0 ) { continue ; }

    // get blue rest-frame flux that we see in redshifted (observer) filter
    lamz  = lam / z1;
    tmp   = lamz / dlam + 0.5 ;
    ilamz = (int)tmp ;

    flux_E  = *(ptrFlux_rest  + ilamz);  // erg/s/cm^2/A
    if ( flux_E == 0.0 )  { Lzeroflux = 1; }
    flux_E *= atmRatio ;

    sumFT    += ( trans * flux_E * lam ) ;
    sumTrans += ( trans / lam ) ;

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
void filter_extinct(void) {

  int ifilt, ilam ;

  double xt, trans ;
  // ----------- BEGIN -----------

  printf("\n");

  if ( INPUTS.FILTER_MODTRAN_FLAG > 0 ) 
    { printf("\t Atmospheric extinction already applied to filter trans \n"); }
  else
    { printf("\t FilterTrans *= (Atmospheric Extinction) \n" ); }

  for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {

    for ( ilam=0; ilam < MXLAMINT; ilam++ ) {

      if ( INPUTS.FILTER_MODTRAN_FLAG > 0 ) 
	{ xt = 1.0; }
      else
	{ xt = MODTRAN.AtmTRANS_REBIN[IMOD_TOTAL][ilam] ; }

      trans = FILTER[ifilt].TRANS_REBIN[ilam] ;
      FILTER[ifilt].TRANSTOT_REBIN[ilam] = trans * xt ; 

    } // ilam
  } // ifilt

} // end of filter_extinct

// ****************************************************************
void rebin_all ( void ) {

  int ifilt, ised, NBIN, imod;
  double *ptrlam, *ptrfun, *ptrfun_rebin ;

  // ------------ BEGIN ----------

  printf("\n");

  printf("\t Rebin Filter Transmissions \n ");
  for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {
    NBIN          = FILTER[ifilt].NBIN_LAM ;
    ptrlam        = &FILTER[ifilt].LAMBDA[1];
    ptrfun        = &FILTER[ifilt].TRANS[1];
    ptrfun_rebin  = FILTER[ifilt].TRANS_REBIN;  
    lamrebin(NBIN,ptrlam, ptrfun, ptrfun_rebin);
  }

  printf("\t Rebin SED Flux  \n");
  for ( ised=1; ised <= INPUTS.NSED; ised++ ) {
    NBIN          = SED[ised].NBIN_LAM ;
    ptrlam        = &SED[ised].LAMBDA[1];
    ptrfun        = &SED[ised].FLUX[1];
    ptrfun_rebin  = SED[ised].FLUX_REBIN;
    lamrebin(NBIN, ptrlam, ptrfun, ptrfun_rebin);
  }

  printf("\t Rebin MODTRAN Extinction  \n");
  for ( imod = 1; imod <= NMOD_COL; imod++ ) {
    NBIN          = MODTRAN.NBIN_LAM ;
    ptrlam        = &MODTRAN.LAMBDA[1] ;
    ptrfun        = &MODTRAN.AtmTRANS[imod][1] ;
    ptrfun_rebin  = MODTRAN.AtmTRANS_REBIN[imod] ;
    lamrebin(NBIN, ptrlam, ptrfun, ptrfun_rebin);
  }


} // rebin_all


// ****************************************************************
void lamrebin(
	      int NBIN          // (I) number of input bins
	      ,double *ptrlam    // (I) input lambda bins
	      ,double *ptrfun   // (I) input contents (trans or SED)
	      ,double *ptrfun_rebin // (O) rebined contents
	      ) {

  /****
    Start with input array ptrfun vs. ptrlam, and rebin
    using existing LAMBDA_REBIN array.  Output array is
    therefore ptrfun_rebin vs. LAMBDA_REBIN. When all filter
    transmissions and SED are rebined to the same 
    lambda-bins, the integrations and filter-flux convolutions
    will be much easier.

   Jul 9, 2011: pass entire input array to interp_1DFUN()
                and get rid of local logic to isolate nearest bin.

  ***/

  double 
    *ptrlam_rebin
    ,LAM_MIN
    ,LAM_MAX
    ,lam_rebin
    ;

  // use double precision for interpolation function
  double
    fun8
    ,interp8_lam[4]
    ,interp8_fun[4]
    ;

  int ilam_rebin, ILAM_REBIN_MIN, ILAM_REBIN_MAX    ;

  // -------- BEGINN ------------
 
  ptrlam_rebin = LAMBDA_REBIN; // useful local pointer

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
    if ( lam_rebin < LAM_MIN ) { continue ; }
    if ( lam_rebin > LAM_MAX ) { continue ; }

    fun8 = interp_1DFUN(OPT_INTERP, lam_rebin, NBIN, ptrlam, ptrfun,
			"lamrebin" );

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

} // end of lamrebin


// *************************************
void AtmosCHANGE(int ichange, double *TRANS ) {

  // for input indices 'ised' and  'ichange',
  // returns modified *sedchange.
  // For 'ichange' corresponding to user-requested ATMOS_CHANGE,
  // return total atmos trans *TRANS.

  double 
     TRtot_default
    ,TRmod_default
    ,TRtot_change
    ,TRmod_change 
    ,xtmp, xchange
    ,one = 1.0 
    ;

  int imod, ilam ;

  //  char fnam[] = "AtmosCHANGE" ;

  // ------------ BEGIN -----------

  for ( ilam=0; ilam < MXLAMINT; ilam++ ) {

    TRtot_default = MODTRAN.AtmTRANS_REBIN[IMOD_TOTAL][ilam] ;
    *(TRANS+ilam) = TRtot_default ;

    if ( TRtot_default <= 0.0 ) continue ;

    // apply all changed components to get changed transmission.

    TRtot_change = TRtot_default ;

    // Tchange = 1 - (1-Tdefault)*scale  for each component
    // Tchange_TOTAL = product(Tchange)

    for ( imod = 1; imod <= NMOD_COL; imod++ ) {
      if ( imod == IMOD_TOTAL ) continue ;

      TRmod_default   = MODTRAN.AtmTRANS_REBIN[imod][ilam] ;
      xchange         = (double)INPUTS.ATMOS_CHANGE[ichange][imod]  ;

      // make change only when requested
      if ( xchange != one  && TRmod_default > 0.0 ) {
	xtmp      = one - TRmod_default ; // extinction
	xtmp     *= xchange ;
	if ( xtmp >= one ) { xtmp = one ; }
	TRmod_change = one  - xtmp ;

	TRtot_change *= ( TRmod_change / TRmod_default ) ; 
      }

    } // imod

    *(TRANS+ilam) = TRtot_change ;

  } // ilam


} // end of AtmosCHANGE


// ******************************************
void get_colorcor(int ichange, int ifilt) {

  // get color correction for this filter and 'change'

  int ised, if0, if1, iz, NREF, itype;
  int OPT ;

  double 
    color_slope, color_off
    ,magdif
    ,magref[MXSEDCAL]
    ,magxt[MXSEDCAL] 
    ,mag0[MXSEDCAL]
    ,mag1[MXSEDCAL]
    ,color[MXSEDCAL]
    ;

  int LDMP = 0 ;
  FILE *fpdmp;
  char c_color[12], dmpFile[200];
  char fnam[] = "get_colorcor";

  // ---------- BEGIN ----------

  // get filter indices used for color correction
  if0 = INPUTS.ifilterColorcor[ifilt][0];
  if1 = INPUTS.ifilterColorcor[ifilt][1];
  iz  = 1 ; //  only redshift=0 (1st index)  for regs

 
  if ( LDMP > 0 ) {

    sprintf(dmpFile,"dumpColorcor-%s-change%d.dat", 
	    INPUTS.filterName[ifilt], ichange );
    fpdmp = fopen(dmpFile, "wt");

    /*
    printf(" -------------------------------------------- \n");
    printf("  DMP: Get %s colorCor for Filter = %s and ichange=%d \n"
	   ,INPUTS.filterColorcorString[ifilt]
	   ,INPUTS.filterName[ifilt], ichange );
    */
  }



  OPT = INPUTS.OPT_COLORCOR ;
  NREF = 0;

  // strip off apparent mag dimming for each SED

  for ( ised=1; ised <= INPUTS.NSED; ised++ ) {

    itype = INPUTS.sedType[ised];
    if ( itype != SEDTYPE_REF ) { continue ; }

    NREF++;

    // mag with nominal atmos extinction
    magref[NREF] = SYNMAG_avecXT[ifilt][ised][iz] ;

    // mag with change in atmos extinction
    magxt[NREF]  = SYNMAG_change_raw[ichange][ifilt][ised][iz];

    // get  mags needed for color correction

    if ( OPT == OPT_COLORCOR_RAW ) {
      mag0[NREF]  = SYNMAG_change_raw[ichange][if0][ised][iz];
      mag1[NREF]  = SYNMAG_change_raw[ichange][if1][ised][iz];
      sprintf(c_color,"raw");
    }
    else if ( OPT == OPT_COLORCOR_TRUE ) {
      mag0[NREF]  = SYNMAG_avecXT[if0][ised][iz] ;
      mag1[NREF]  = SYNMAG_avecXT[if1][ised][iz] ;
      sprintf(c_color,"true");
    }
    else {
      sprintf(c1err,"Invalid  OPT_COLORCOR = %d", OPT );
      errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
    }

    color[NREF]  = mag1[NREF]  - mag0[NREF] ;

    if ( LDMP > 0 ) {
      magdif = magref[NREF] - magxt[NREF] ;
      fprintf(fpdmp,"%f  %f  %f  %f \n", 
	      color[NREF], magxt[NREF], magref[NREF], magdif );
      fflush(fpdmp);

      /**  xxxxxxxxxxxx
      sprintf(ctmp,"%s", INPUTS.sedName[ised] );
      printf("\t DMP %6s: Mref-mxt=%9.5f   xtcolor=%8.4f \n"
	     ,ctmp, magdif, color[NREF]  ); 
	     xxxxxxxxxxxxx */
      
    }
  }

  // mref = mxt + b * [ (m1-m0) - color_off ]


  fitColorcor(ichange, NREF, &magref[1], &magxt[1], &color[1], 
	      &color_slope, &color_off);  // outputs
  

  FILTER[ifilt].color_slope[ichange]  = color_slope ;
  FILTER[ifilt].color_offset[ichange] = color_off ;

  printf("  Atmos Change %2d : %s_cor = %s_raw + %11.3le [ (%s-%s)_%s - %7.3f ] \n"
	 ,ichange
	 ,INPUTS.filterName[ifilt]
	 ,INPUTS.filterName[ifilt]
	 ,color_slope
	 ,INPUTS.filterName[if1]
	 ,INPUTS.filterName[if0]
	 ,c_color
	 ,color_off
	 );

  fflush(stdout) ;


} // end of get_colorcor


// *************************************
void fitColorcor(int ichange, int NREF, 
		 double *mref, double *mxt, double *color,  // inputs
		 double *color_slope, double *color_off     // outputs
		 ) {

  // Solve for the color slope & offset for a set of
  // Mref, Mxt and color,
  //
  // mref[i] = mxt[i] + color_slope * [ color[i] - color_off ]
  //

  int iref, irefmin, irefmax;
  double 
    c, colormin, colormax
    ,dm, dc
    ,MDIF[MXSEDCAL]
    ,slope, offmin, offmax
    ;
  
  // --------------- BEGIN ------------

  *color_slope =  0.0 ; // init outputs
  *color_off   =  0.0 ;

  // quick & dirty; use extreme two points and draw line

  irefmin = irefmax = -9 ;
  colormin =  999999.9 ;
  colormax = -999999.9 ;

  //  printf("  000000 %d ----------------------- \n", ichange); // DDDDDDDDDDDD
  for ( iref = 0 ; iref < NREF; iref++ ) {  // NREF-1 ==> NREF
    MDIF[iref] = mref[iref] - mxt[iref] ;
    c =  color[iref];

    if ( fabs(c) < 1.0E-8 ) {  continue ; }
    // printf("  0000000 %f  %f  \n", c, MDIF[iref] ); // DDDDDDDD
    if  ( c < colormin ) { colormin = c ; irefmin = iref ; }
    if  ( c > colormax ) { colormax = c ; irefmax = iref ; }
  }

  if ( irefmin < 0 ) { return ; }
  if ( irefmax < 0 ) { return ; }

  dm = MDIF[irefmax]  - MDIF[irefmin];
  dc = color[irefmax] - color[irefmin] ;

  if ( dc == 0.0 ) { return ; }


  slope = dm/dc ;
  if ( fabs(slope) > 1.0E-6 ) {
    offmin = color[irefmin] - MDIF[irefmin]/slope;
    offmax = color[irefmax] - MDIF[irefmax]/slope;
  }
  else {
    offmin = offmax = 0.0 ;
  }

  *color_slope = slope ;
  *color_off   = (offmin+offmax)/2.0 ;

} // end of fitColorcor


// *************************************
void  mkplots(void) {

  // Write results to hbook or root file.
  // Apr 2013: need to switch to using sntools_output
  // May 04, 2013: use OPEN_HFILE and CLOSE_HFILE

  char chis[80];
  char *ptrHis;  
  double mag, magref;

  int 
    hid, hid2, hoff_change
    ,hidraw, hidcor
    ,ised, ifilt, imod, ichange
    ,nblam, ilam, nbz, iz
    ,nb, nbf
    ,NSED
    ;

  /*
  float minlam4, maxlam4, lam4, w4, magdif_raw4, magdif_cor4;
  float zmin4, zmax4, z4, zhalf4;
  float t4, xt4, flux4, xfilt4, xsed4 ;
  float xmin4, xmax4, xmod4 ;
  float fmin4, fmax4, smin4, smax4 ;
  */

  double minlam8, maxlam8, lam8, magdif_raw8, magdif_cor8;
  double zmin8, zmax8, z8, zhalf8;
  double t8, xt8, flux8, xfilt8, xsed8 ;
  double xmin8, xmax8, xmod8 ;
  double fmin8, fmax8, smin8, smax8 ;

  // ------------- BEGIN -------------

  ptrHis = INPUTS.hbook_outFile ;

  TABLEFILE_INIT();
  TABLEFILE_OPEN(ptrHis,"new");

  // hard-wire lambda binning to cover UV+optical

  minlam8 = 2000. ;
  maxlam8 = 12000. ;
  nblam   = (int)( (maxlam8-minlam8)/INPUTS.LAMINT_BINSIZE ) ;
  NSED = INPUTS.NSED ;

  // start with atmosphere transmisson vs. lambda by component
  for ( imod=1; imod <= NMOD_COL; imod++ ) {
    hid = imod ;
    sprintf(chis,"Atmos Trans from %s", MODTRAN_TYPE[imod] );
    SNHIST_INIT(1, hid, chis, &nblam, &minlam8, &maxlam8 );

    for ( ilam=0 ; ilam <= MXLAMINT; ilam++ ) {
      lam8 = (double)LAMBDA_REBIN[ilam];

      if ( lam8 >= minlam8 && lam8 <= maxlam8 ) {
	t8   = MODTRAN.AtmTRANS_REBIN[imod][ilam] ;
	SNHIST_FILL(1, hid, &lam8, t8 );
      }

    } // ilam
  } // imod



  // now the SEDs
  for ( ised = 1; ised <= NSED ; ised++ ) {
    hid = 20 + ised;
    sprintf(chis,"%s SED flux vs. wavelength", INPUTS.sedName[ised] );
    SNHIST_INIT(1, hid, chis, &nblam, &minlam8, &maxlam8 );

    for ( ilam=0 ; ilam <= MXLAMINT; ilam++ ) {
      lam8   = (double)LAMBDA_REBIN[ilam];
      if ( lam8 < minlam8 ) continue ;
      if ( lam8 > maxlam8 ) continue ;
      flux8  = (double)SED[ised].FLUX_REBIN[ilam];
      SNHIST_FILL(1, hid, &lam8, flux8 );
    }
  }


  hid = 100 ;
  sprintf(chis, "mean wavelength vs. IFILTER index");
  nbf = INPUTS.NFILT ;
  fmin8 = 0.5; fmax8 = fmin8 + (double)nbf ;
  SNHIST_INIT(1, hid, chis, &nbf, &fmin8, &fmax8 );

  // now the filter transmissions with/without extinction
  for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {

    xfilt8  = (double)ifilt;
    lam8    = (double)FILTER[ifilt].LAMBDA_MEAN  ;
    hid = 100;   SNHIST_FILL(1, hid, &xfilt8, lam8 );

    hid = 100 + ifilt ;
    sprintf(chis,"%s band transmission from file", 
	    INPUTS.filterName[ifilt] );
    SNHIST_INIT(1, hid, chis, &nblam, &minlam8, &maxlam8 );

    hid = 200 + ifilt ;
    sprintf(chis,"%s band transmission * extinction", 
	    INPUTS.filterName[ifilt] );
    SNHIST_INIT(1, hid, chis, &nblam, &minlam8, &maxlam8 );

    for ( ilam=0 ; ilam <= MXLAMINT; ilam++ ) {
      lam8 = (double)LAMBDA_REBIN[ilam];
      t8   = (double)FILTER[ifilt].TRANS_REBIN[ilam] ;
      xt8  = (double)FILTER[ifilt].TRANSTOT_REBIN[ilam] ;

      if ( lam8 >= minlam8 && lam8 <= maxlam8 ) {
	hid = 100 + ifilt; SNHIST_FILL(1, hid, &lam8, t8  );  
	hid = 200 + ifilt; SNHIST_FILL(1, hid, &lam8, xt8 );
      }
    }


  } // ifilt



  //  plot mag-diff  vs. redshift, where mag-dif = mag(XT) - mag(no XT)

  zhalf8 = INPUTS.REDSHIFT_BINSIZE/2.0 ;
  nbz    = INPUTS.NBIN_REDSHIFT ;
  zmin8  = INPUTS.REDSHIFT_RANGE[0] - zhalf8 ;
  zmax8  = INPUTS.REDSHIFT_RANGE[1] + zhalf8 ;

  if ( zmax8 == 0.0 && zmin8 == 0.0 ) {
    zmin8 = -0.01; zmax8 = 0.01 ;
  }

  for ( ichange = 1; ichange <= INPUTS.NCHANGES; ichange++ ) {

    // plot atm trans from change 

    hoff_change = 10000*ichange ;
    
    hid = hoff_change  ;
    sprintf(chis,"AtmosTrans(total) for ATMOS-CHANGE(%d)", ichange);
    SNHIST_INIT(1, hid, chis, &nblam, &minlam8, &maxlam8 );

    for ( ilam=0 ; ilam <= MXLAMINT; ilam++ ) {
      lam8 = (double)LAMBDA_REBIN[ilam];
      if ( lam8 < minlam8 ) continue ;
      if ( lam8 > maxlam8 ) continue ;
      xt8  = (double)MODTRAN.AtmTRANS_CHANGE[ichange][ilam];
      SNHIST_FILL(1, hid, &lam8, xt8  );  
    } // ilam

    // now book atm change by component
    hid = hoff_change  + 1;
    nb = NMOD_COL ;
    xmin8 = 0.5; xmax8 = xmin8 + (double)(nb);
    sprintf(chis,"Atmos XT(%d)/XT(default) vs. component", ichange);
    SNHIST_INIT(1, hid, chis, &nb, &xmin8, &xmax8 );

    for ( imod=1; imod <= NMOD_COL; imod++ ) {
      xmod8 = (double)imod ;
      xt8 = INPUTS.ATMOS_CHANGE[ichange][imod];
      if ( imod != IMOD_TOTAL ) { SNHIST_FILL(1, hid, &xmod8, xt8  ); }
    }

    smin8 = 0.5; smax8 = smin8 + (double)NSED ;
    for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {
      hid = hoff_change  + 10 + ifilt ;
      sprintf(chis,"raw [D]m(%s) vs. SED index for ichange=%d (redshift=0)", 
	      INPUTS.filterName[ifilt], ichange );
      SNHIST_INIT(1, hid, chis, &NSED, &smin8, &smax8 );

      hid = hoff_change  + 20 + ifilt ;
      sprintf(chis,"cor [D]m(%s) vs. SED index for ichange=%d (redshift=0)", 
	      INPUTS.filterName[ifilt], ichange );
      SNHIST_INIT(1, hid, chis, &NSED, &smin8, &smax8 );
    }

    for ( ised=1; ised <= NSED; ised++ ) {
      xsed8 = (double)ised ;

      hid2 = hoff_change  + 100*ised + 90 ;
      sprintf(chis,"raw [D]m vs. ifilt index (z=0) for %s (AtmosChange=%d)"
	      ,INPUTS.sedName[ised],ichange );
      SNHIST_INIT(1, hid2, chis, &nbf, &fmin8, &fmax8 );

      for ( ifilt=1; ifilt <= INPUTS.NFILT; ifilt++ ) {

	hidraw = hoff_change  + 100*ised + ifilt ;
	sprintf(chis,"raw [D]m?%s!  vs. redshift for %s (AtmosChange=%d)"
		,INPUTS.filterName[ifilt]
		,INPUTS.sedName[ised]
		,ichange
		);
	SNHIST_INIT(1, hidraw, chis, &nbz, &zmin8, &zmax8 );

	hidcor = hidraw + 20;
	sprintf(chis,"cor [D]m?%s!  vs. redshift for %s (AtmosChange=%d)"
		,INPUTS.filterName[ifilt]
		,INPUTS.sedName[ised]
		,ichange
		);
	SNHIST_INIT(1, hidcor, chis, &nbz, &zmin8, &zmax8 );

	xfilt8 = (double)ifilt ;

	for ( iz=1; iz <= INPUTS.NBIN_REDSHIFT; iz++ ) {

	  if ( INPUTS.sedType[ised] == SEDTYPE_REF && iz > 1 ) continue ;

	  z8     = INPUTS.ARRAY_REDSHIFT[iz];
	  magref = SYNMAG_avecXT[ifilt][ised][iz] ;

	  mag         = SYNMAG_change_raw[ichange][ifilt][ised][iz];
	  magdif_raw8 = (double)(mag - magref) ;
	  mag         = SYNMAG_change_colorcor[ichange][ifilt][ised][iz];
	  magdif_cor8 = (double)(mag - magref) ;

	  SNHIST_FILL(1, hidraw, &z8, magdif_raw8  ); 
	  SNHIST_FILL(1, hidcor, &z8, magdif_cor8  ); 
	  
	  if ( iz == 1 ) {
	    hid = hoff_change  + 100*ised + 90 ;
	    SNHIST_FILL(1, hid, &xfilt8, magdif_raw8  ); 

	    hid = hoff_change  + 10 + ifilt ;
	    SNHIST_FILL(1, hid, &xsed8, magdif_raw8  ); 

	    hid = hoff_change  + 20 + ifilt ;
	    SNHIST_FILL(1, hid, &xsed8, magdif_cor8  ); 
	  }

	} // iz
      }  // ifilt
    } // ised
  } // ichange


  TABLEFILE_CLOSE(ptrHis);

  printf("\n Results written to hbook file:  %s \n", 
	 INPUTS.hbook_outFile );
  fflush(stdout);


} // end of   mkplots

