/**************************************************
  Created Jun 1 2018 by R.Kessler

  Rebin SED files in SIMSED directory to reduce file size
  and memory usage. Input SIMSED dir must already be in 
  SNANA format. Rebinning is performed by discarding bins;
  so that information is removed but not altered.
  Since lambda bins must be uniform, a lambda-cut option is provided.
  No need for DAY cut because DAY bins can be increased vs. DAY.

  Example:
    SIMSED_rebin.exe <inpDir>  -rebin_day 2   
      --> keep every other day bin

    SIMSED_rebin.exe <inpDir>  -rebin_day 2  -rebin_lam 2
      --> keep every other day bin and every other lambda bin

    SIMSED_rebin.exe <inpDir>  '-rebin_day(40:400)' 5
      --> keep every 5th day bin for 40 < day < 400

    SIMSED_rebin.exe <inpDir>  '-rebin_day(40:100)' 2  '-rebin_day(100:400)' 5
      --> keep every 2nd day bin for 40 < day < 100, and then
          keep every 4th day bin for 100<day < 400.

    SIMSED_rebin.exe <inpDir>  -lamwin 2000 12000
      --> only keep wavelengths 2000 to 12000 A

  BEWARE NOTES:
    + day rebinning can be a function of day,
      by LAM rebinning must be uniform.
    + must use single quotes around arguments with ()

 ***********/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
//#include "sntools_cosmology.h"
#include "MWgaldust.h"

#include "genmag_SIMSED.h"
#include "genmag_SEDtools.h"


// ========= GLOBAL DECLARATIONS ============

#define MXREBIN_DAY 4
#define LAMMAX_REBIN 1.0E5

#define	MXBIN_DAY_REBIN 2000
#define	MXBIN_LAM_REBIN 5000
#define MXBIN_SED_REBIN MXBIN_DAY_REBIN*MXBIN_LAM_REBIN

struct {
  char   INPDIR_SIMSED[MXPATHLEN];
  char   OUTDIR_SIMSED[MXPATHLEN];  // determined from INPDIR

  int    NREBIN_DAY;
  int    REBIN_DAY[MXREBIN_DAY];
  double DAYRANGE_REBIN[MXREBIN_DAY][2];

  double CUTWIN_LAM[2];
  int REBIN_LAM;

} INPUTS;


char SEDINFO_FILE_INP[MXPATHLEN];
char SEDINFO_FILE_OUT[MXPATHLEN];

FILE  *FP_SEDINFO_INP;
FILE  *FP_SEDINFO_OUT;

SEDMODEL_FLUX_DEF  TEMP_SEDMODEL_INP ;
SEDMODEL_FLUX_DEF  TEMP_SEDMODEL_OUT ;

void  parse_args(int argc, char **argv);
void  parse_REBIN_DAY(char *string, int rebin);
void  print_rebin_info(FILE *fp);
void  make_OUTDIR_SIMSED(void);
void  SIMSED_DRIVER(void);
void  write_rebinned_SED(char *sed_fileName);
void  rebin_SED(void);

// ====================================
int main(int argc, char **argv) {

  char fnam[] = "main" ;
  // ------------ BEGIN ------------

  set_EXIT_ERRCODE(EXIT_ERRCODE_UNKNOWN);

  sprintf(BANNER,"Begin execution of SIMSED_rebin" );
  print_banner(BANNER);

  // parse user arguments
  parse_args(argc, argv );

  sprintf(SEDINFO_FILE_INP, "%s/SED.INFO", INPUTS.INPDIR_SIMSED);
  sprintf(SEDINFO_FILE_OUT, "%s/SED.INFO", INPUTS.OUTDIR_SIMSED);

  // create output directory and open new SED.INFO file
  make_OUTDIR_SIMSED();

  // allocate TEMP array to read SED 
  malloc_SEDFLUX_SEDMODEL(&TEMP_SEDMODEL_INP,
			  MXBIN_DAY_REBIN, MXBIN_LAM_REBIN, MXBIN_SED_REBIN);

  // malloc rebinned array
  malloc_SEDFLUX_SEDMODEL(&TEMP_SEDMODEL_OUT,
			  MXBIN_DAY_REBIN, MXBIN_LAM_REBIN, MXBIN_SED_REBIN);


  SIMSED_DRIVER();

  fclose(FP_SEDINFO_INP);
  fclose(FP_SEDINFO_OUT);

} // end main


// ******************************************
void parse_args(int argc, char **argv) {

  int i, ipar,  ep, ilast, iuse, rebin_tmp ;
  char fnam[] = "parse_args" ;


  // ---------- BEGIN --------

  NARGV_LIST = argc;
  for( i=0; i < NARGV_LIST; i++ ) {
    USE_ARGV_LIST[i] = 0 ;
    ARGV_LIST[i] = (char*) malloc( MXPATHLEN*sizeof(char) );
    sprintf(ARGV_LIST[i], "%s", argv[i] );
  }
  
  INPUTS.REBIN_LAM = 1;
  INPUTS.NREBIN_DAY = 0 ;
  INPUTS.CUTWIN_LAM[0] = 0.0 ;
  INPUTS.CUTWIN_LAM[1] = LAMMAX_REBIN ;

  i=1; ilast=i ;
  sscanf(ARGV_LIST[i], "%s", INPUTS.INPDIR_SIMSED );
  sprintf(INPUTS.OUTDIR_SIMSED, "%s_rebin", INPUTS.INPDIR_SIMSED );
  USE_ARGV_LIST[i] = 1;
  i++; ilast=i ;

  while ( i < NARGV_LIST ) {
   
    if ( strcmp(ARGV_LIST[i],"-rebin_lam") == 0 ) {
      i++;  sscanf(ARGV_LIST[i], "%d", &INPUTS.REBIN_LAM );
    }

    if ( strcmp(ARGV_LIST[i],"-cutwin_lam") == 0 ) {
      i++;  sscanf(ARGV_LIST[i], "%le", &INPUTS.CUTWIN_LAM[0] );
      i++;  sscanf(ARGV_LIST[i], "%le", &INPUTS.CUTWIN_LAM[1] );
    }

    if ( strstr(ARGV_LIST[i],"-rebin_day") != NULL ) {
      i++ ; sscanf(ARGV_LIST[i], "%d",  &rebin_tmp );
      parse_REBIN_DAY(ARGV_LIST[i-1], rebin_tmp);
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
  printf("   OUTDIR_SIMSED: %s \n", INPUTS.OUTDIR_SIMSED );

  print_rebin_info(stdout);

  printf("\n\n");

  return ;

} // end parse_args


// ******************************************
void  parse_REBIN_DAY(char *string, int rebin) {

  // parse string = -rebin_day([daymin]:[daymax])
  // and store in INPUTS arrays.

  int    irebin,  NSPLIT, j ;
  double DAYRANGE[2];
  char stringDayRange[100], *ptrSPLIT[2], stringSPLIT[2][40];
  char colon[] = ":";
  char fnam[] = "parse_REBIN_DAY";

  // --------- BEGIN ---------

  extractStringOpt(string, stringDayRange);

  /*
  printf(" xxx string         = '%s' \n", string);
  printf(" xxx stringDayRange = '%s' \n", stringDayRange );
  */

  if ( strlen(stringDayRange) == 0 ) {
    // just one rebin_day for all days
    INPUTS.NREBIN_DAY=1; 
    irebin=0; 
    INPUTS.REBIN_DAY[irebin] = rebin;
    INPUTS.DAYRANGE_REBIN[irebin][0] = -20000.0 ;
    INPUTS.DAYRANGE_REBIN[irebin][1] = +20000.0 ;
  }
  else {
    // check specif day range for this rebin factor
    ptrSPLIT[0] = stringSPLIT[0];
    ptrSPLIT[1] = stringSPLIT[1];
    splitString(stringDayRange, colon, fnam, 3,
		&NSPLIT, ptrSPLIT);  // <== returned
    if ( NSPLIT != 2 ) {
      sprintf(c1err,"Expected 2 values in '%s'", stringDayRange);
      sprintf(c2err,"but found %d values.", NSPLIT);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    irebin = INPUTS.NREBIN_DAY ;
    INPUTS.REBIN_DAY[irebin] = rebin;
    for(j=0; j< 2; j++ ) {
      sscanf(stringSPLIT[j], "%le", &DAYRANGE[j] );
      INPUTS.DAYRANGE_REBIN[irebin][j] = DAYRANGE[j] ;
    }

    INPUTS.NREBIN_DAY++ ;

  }

  return ;

} // end parse_REBIN_DAY


// ***************************************
void  print_rebin_info(FILE *fp) {

  int i;

  if ( INPUTS.CUTWIN_LAM[0] > 1.0 ) {
    fprintf(fp, "#   CUTWIN_LAM:    %.0f to %.0f \n",
	    INPUTS.CUTWIN_LAM[0], INPUTS.CUTWIN_LAM[1] );
  }

  if ( INPUTS.REBIN_LAM > 1 ) {
    fprintf(fp, "#   REBIN_LAM:     %d \n", INPUTS.REBIN_LAM );
  }

  for(i=0; i < INPUTS.NREBIN_DAY; i++ ) {
    fprintf(fp, "#   REBIN_DAY:     %d for %.1f < DAY < %.1f \n"
	   , INPUTS.REBIN_DAY[i]
	   , INPUTS.DAYRANGE_REBIN[i][0]
	   , INPUTS.DAYRANGE_REBIN[i][1] );
  }

  fprintf(fp,"#\n");
  fflush(fp);

  return;

} // end print_rebin_info


// *************************************
void  make_OUTDIR_SIMSED(void) {

  // Create output SIMSED dir, and open SED.INFO file.
  FILE *fp;
  char cmd_make[400], cmd_rm[400] ;
  // ----------- BEGIN ---------

  fp = fopen(SEDINFO_FILE_OUT,"rt");
  if ( fp ) {
    sprintf(cmd_rm,"rm -r %s/", INPUTS.OUTDIR_SIMSED);
    if ( strlen(INPUTS.OUTDIR_SIMSED) > 5 ) 
      { system(cmd_rm); }
  }

  sprintf(cmd_make, "mkdir -m g+wr %s", INPUTS.OUTDIR_SIMSED);
  system(cmd_make);
  printf("  mkdir %s \n", INPUTS.OUTDIR_SIMSED);

  FP_SEDINFO_OUT = fopen(SEDINFO_FILE_OUT,"wt");
  print_rebin_info(FP_SEDINFO_OUT) ;

  return ;

} // end void  make_OUTDIR_SIMSED



// ****************************
void  SIMSED_DRIVER(void) {

#define MXCHAR_LINE 200

  int NLINE=0, NSED=0 ;
  char LINE[MXCHAR_LINE];
  char *ptrtok, KEY[40], sed_fileName[100], SED_FILENAME[200];
  char sedComment[] = "";
  char fnam[] = "SIMSED_DRIVER" ;

  int    SEDMODEL_OPTMASK = 2;  // allow non-uniform day bins
  double Trange_SIMSED[2] = { -1000.0, 1000.0 } ;
  double Lrange_SIMSED[2] = {  100.0,  LAMMAX_REBIN } ;

  // ---------- BEGIN -------

  FP_SEDINFO_INP = fopen(SEDINFO_FILE_INP, "rt");
  if ( !FP_SEDINFO_INP ) {
    sprintf(c1err,"Could not open input SED.INFO file:");
    sprintf(c2err,"%s", SEDINFO_FILE_INP);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  printf(" Begin reading input SED.INFO file: \n");

  while ( fgets(LINE, MXCHAR_LINE, FP_SEDINFO_INP) != NULL ) {
    
    // copy original LINE to new SED.INFO file
    fprintf(FP_SEDINFO_OUT, "%s", LINE);  NLINE++;

    // extract first two words of string; looking for "SED:" and file name.
    ptrtok = strtok(LINE," ");
    if ( ptrtok != NULL ) { sprintf(KEY, "%s", ptrtok); }
    if ( strcmp(KEY,"SED:") != 0 ) { continue; }

    ptrtok = strtok(NULL," " );
    if ( ptrtok != NULL ) { sprintf(sed_fileName, "%s", ptrtok); }
    
    sprintf(SED_FILENAME, "%s/%s", INPUTS.INPDIR_SIMSED, sed_fileName);
    rd_sedFlux(SED_FILENAME, sedComment
	     ,Trange_SIMSED, Lrange_SIMSED
	     ,MXBIN_DAY_REBIN, MXBIN_LAM_REBIN
	     ,SEDMODEL_OPTMASK
	     ,&TEMP_SEDMODEL_INP.NDAY, 
	       TEMP_SEDMODEL_INP.DAY, &TEMP_SEDMODEL_INP.DAYSTEP
	     ,&TEMP_SEDMODEL_INP.NLAM, 
	       TEMP_SEDMODEL_INP.LAM, &TEMP_SEDMODEL_INP.LAMSTEP
	     ,TEMP_SEDMODEL_INP.FLUX,  TEMP_SEDMODEL_INP.FLUXERR );  


    // do the dirty work of rebinng and modify TEMP_SEDMODEL
    rebin_SED();

    write_rebinned_SED(sed_fileName);

    NSED++ ;
  }

  printf("  Done readig %d lines and %d SEDs from SED.INFO file.\n", 
	 NLINE, NSED );

  return ;

} // end void  SIMSED_DRIVER


// **********************************************
void  rebin_SED(void) {

  // Do the rebinning here and modify global TEMP_SEDMODEL
  //

  int NDAY_INP = TEMP_SEDMODEL_INP.NDAY;
  int NLAM_INP = TEMP_SEDMODEL_INP.NLAM;

  int *KEEP_DAY, *KEEP_LAM ;
  int NDAY_OUT, NLAM_OUT;
  int iday, ilam, iday_out, ilam_out, jflux, nskip, irebin, REBIN;
  double DAY, LAM, FLUX, DAYMIN, DAYMAX ;
  char fnam[] = "rebin_SED" ;

  // ---------------- BEGIN -------------

  // allocate flag per DAY and LAM bin indicatig KEEP or DISCARD

  KEEP_DAY = (int*) malloc ( NDAY_INP * sizeof(int) );
  KEEP_LAM = (int*) malloc ( NLAM_INP * sizeof(int) );

  NDAY_OUT = NLAM_OUT = 0 ;

  // select DAY bins to keep
  nskip = 0 ;
  for(iday=0; iday < NDAY_INP; iday++ ) {
      KEEP_DAY[iday] = 0;
      DAY   = TEMP_SEDMODEL_INP.DAY[iday];

      // find rebin factor based on DAY
      REBIN=1;  // default is keep every DAY
      for(irebin=0; irebin < INPUTS.NREBIN_DAY; irebin++ ) {
	DAYMIN = INPUTS.DAYRANGE_REBIN[irebin][0] ;
	DAYMAX = INPUTS.DAYRANGE_REBIN[irebin][1] ;
	if ( DAY >= DAYMIN && DAY <= DAYMAX ) 
	  { REBIN = INPUTS.REBIN_DAY[irebin] ; }
      }

      nskip++ ;
      if ( nskip != REBIN )   
	{ continue ; }
      else                               
	{ nskip = 0 ; }

      KEEP_DAY[iday] = 1;  NDAY_OUT++ ;

  }
  TEMP_SEDMODEL_OUT.NDAY = NDAY_OUT ;


  // select lam bins to keep
  nskip = 0 ;
  for (ilam=0; ilam < NLAM_INP; ilam++ ) {
    KEEP_LAM[ilam] = 0;
    LAM   = TEMP_SEDMODEL_INP.LAM[ilam];
    if ( LAM < INPUTS.CUTWIN_LAM[0] ) { continue ; }
    if ( LAM > INPUTS.CUTWIN_LAM[1] ) { continue ; }

    nskip++ ;
    if ( nskip != INPUTS.REBIN_LAM )   { continue ; }
    else                               { nskip = 0 ; }

    KEEP_LAM[ilam] = 1;  NLAM_OUT++ ;
  }
  TEMP_SEDMODEL_OUT.NLAM = NLAM_OUT ;

  // ------------------------------------------
  // transfer SEDMODEL_INP to SEDMODEL_OUT

  iday_out = ilam_out = 0 ;
  for(iday=0; iday < NDAY_INP; iday++ ) {

    ilam_out = 0 ;
    if ( KEEP_DAY[iday] == 0 ) { continue ; }

    for (ilam=0; ilam < NLAM_INP; ilam++ ) {
      if ( KEEP_LAM[ilam] == 0 ) { continue ; }

      jflux = NLAM_INP*iday + ilam;
      DAY   = TEMP_SEDMODEL_INP.DAY[iday] ;
      LAM   = TEMP_SEDMODEL_INP.LAM[ilam] ;
      FLUX  = TEMP_SEDMODEL_INP.FLUX[jflux];
      
      jflux = NLAM_OUT*iday_out + ilam_out ;
      TEMP_SEDMODEL_OUT.DAY[iday_out] = DAY ;
      TEMP_SEDMODEL_OUT.LAM[ilam_out] = LAM ;
      TEMP_SEDMODEL_OUT.FLUX[jflux]   = FLUX ;
      
      ilam_out++ ;
    } // end ilam
    iday_out++ ;

  } // end iday

  
  free(KEEP_DAY); free(KEEP_LAM);

  return;

} // end rebin_SED

// **********************************************
void  write_rebinned_SED(char *sed_fileName) {

  // write TEMP_SEDMODEL_OUT array to sed_fileName in output directory.

  int NDAY = TEMP_SEDMODEL_OUT.NDAY;
  int NLAM = TEMP_SEDMODEL_OUT.NLAM;
  
  int  iday, ilam, jflux;
  char SED_FILENAME[MXPATHLEN];
  double DAY, LAM, FLUX;
  FILE *fp;

  char fnam[] = "write_rebinned_SED" ;
  // ---------- BEGIN --------


  // construct full SED file name, including output path
  sprintf(SED_FILENAME, "%s/%s", INPUTS.OUTDIR_SIMSED, sed_fileName);
  
  fp = fopen(SED_FILENAME,"wt");
  if ( !fp ) {
    sprintf(c1err,"Could not open output SED file:");
    sprintf(c2err,"%s", SED_FILENAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  for(iday=0; iday < NDAY; iday++ ) {
    for (ilam=0; ilam < NLAM; ilam++ ) {
      jflux = NLAM*iday + ilam;
      DAY   = TEMP_SEDMODEL_OUT.DAY[iday];
      LAM   = TEMP_SEDMODEL_OUT.LAM[ilam];
      FLUX  = TEMP_SEDMODEL_OUT.FLUX[jflux];
      fprintf(fp,"%7.3f  %8.3f  %10.4E \n", DAY, LAM, FLUX );
    }
  }

  fclose(fp);

  return ;

} // end void  write_rebinned_SED
