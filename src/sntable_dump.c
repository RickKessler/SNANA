/********************************************
  sntable_dump.c

  Created Apr 2013 by R.Kessler  
  [based on example root code by B.Dilday in 
      /home/s1/rkessler/snana_debug/rootTreeDump/src/dumpTree.cc

  The purpose of this program is to extract table info into a text 
  file for snana users who do not want to use root for their analysis.
  The table is assumed to be an hbook-ntuple for .hbook extension, 
  or a root-tree for a .root extension.  An optional argument 
  HBOOK or ROOT can be used to explicitly define which format.

  This code is intended to be used by perl scripts in $SNANA_DIR/util
  but it can also run in stand-alone mode.

  USAGE:
    sntable_dump.exe <tableFile>  
       (print list of tables)

    sntable_dump.exe <tableFile>  <tableName> 
       (print list of variables, then quit)

  Options:
    sntable_dump.exe <tableFile>  <tableName>  -V <varNames>
    sntable_dump.exe <tableFile>  <tableName>  -v <varNames>   

    sntable_dump.exe <tableFile>  <tableName>  -O  <outFile>
    sntable_dump.exe <tableFile>  <tableName>  -o  <outFile>
    sntable_dump.exe <tableFile>  <tableName>  -o  stdout  ! screen dump
                                              
    sntable_dump.exe <tableFile>  <tableName>  NOHEADER
       (multi-column table with no header info and no 'SN:' key)

    sntable_dump.exe <tableFile>  <tableName>  -o  <outFile> -format csv
        (csv format)

    sntable_dump.exe <tableFile>  <tableName>  NEVT
        (print number of events in table)

    sntable_dump.exe <tableFile>  <tableName>  -outlier_fit 3 4
    sntable_dump.exe <tableFile>  <tableName>  -outlier_sim 3 4
       (print 3-4 sigma outliers: '-outlier' for fit-data,
          '--outlier_sim' for data-sim outliers)

   sntable_dump.exe <tableFile>  <tableName>  OBS
     (dump all observations, intended for fluxerrmap-analysis)
  
   If only the tableFile and tableName are given, then a list of
   each variable is given.  If an optional list of <varNames>
   is given, then the values are dumped (for each SN) into a 
   self-documented (fitres-formatted) text file.  If the <outFile>
   option is not given, then a default filename 'sntable_dump.fitres'
   is created.

   For hbook file, tableName is the integer id; 
   for root file, tableName is the string name of the table.

   Example hbook commands:
   sntable_dump.exe snana.hbook  7100  -v CCID RA DECL TYPE
   sntable_dump.exe snfit.hbook  7788  -v CCID RA DECL TYPE -o snfit.fitres

   Example root commands:
   sntable_dump.exe snana.root  SNANA  -v CCID RA DECL TYPE
   sntable_dump.exe snfit.root  FITRES -v CCID RA DECL TYPE -o snfit.fitres



                            HISTORY

  Apr 26, 2014: replace HBOOK-specific and ROOT-specific functions
                with wrappers in sntools_output.c[h].

  Aug 2 2014: new command-line option '-OUTLIER <n1> <n2>' to 
              print fit-outliers between n1 and n2-sigma.

  Aug 4,2014  allow -v arguments to append -outlier variables.

  Oct 26 2014: use refactored system (no changes, just re-compile)

  Jul 20 2016: allow  "-o stdout" to dump fitres output to screen 
               instead of to a file. Done for SNANA_tester script.

  Jul 21 2017: optional argumen "-format IGNORE"
               New function write_IGNORE_FILE();

  Aug 4 2017:
     + for SIMLIB table, fisrt colum is ROW, not CID.
     + check for both CCID and ROW

  Aug 17 2017:
    + checking "CCID ROW" causes problem with APPEND_FITRlES,
      so remove ROW-check for now as quick fix.

 Feb 22 2018: new OBS argument to dump all observations

 FEb 24 2019:
  +  allow CID as first varname in -v list.
  +  new option  "-format csv"  

 Nov 6 2019: 
    + add FLUXERRCALC_SIM to outlier table
    + increase MXVAR_DUMP 20 -> 40
    + abort trap if NVAR > MXVAR_DUMP.

 Dec 28 2019: for OBS option, include FLUXCAL_MODEL

 May 13 2020: replace --v with -v, etc ...

 Aug 08 2020: new NEVT option to print number of events

 Mar 14 2021: for outlier output table, write string BAND before IFILTOBS.
              See OPT arg to load_DUMPLINE in sntools_output_hbook[root].c

 Jun 22 2021: replace 200 -> MXPATHLEN for input file names.
                [fixes failure found by Dillon]

********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "sntools.h"
#include "sntools_output.h"

// ----------------------------

char msgerr1[80], msgerr2[80];


#define MXVAR_DUMP 40 

struct INPUTS {
  char TABLE_FILE[MXPATHLEN];
  char TABLE_ID[60];
  bool IS_SNANA, IS_FITRES;
  char VARNAMES[MXVAR_DUMP][60] ;
  char VARLIST[MXVAR_DUMP*60];     // space-separate list
  char OUTFILE_FITRES[MXPATHLEN] ;
  char OUTFILE_IGNORE[MXPATHLEN] ;   // if --format IGNORE
  char FORMAT_OUTFILE[40];   
  int  ISFORMAT_CSV ;

  int  ADD_HEADER ;  // add fitres header of NVAR > 0
  int  NVAR ;        // number of VARNAMES variables to dump
  int  IVAR_NPT;     // ivar index for NOBS or NPTFIT

  float OUTLIER_NSIGMA[2] ;   // Nsigma range for OUTLIER option
  char  OUTLIER_VARNAME_CHI2FLUX[40];
  char  OUTLIER_VARNAME_FLUXCAL[40]; // FLUXCAL_MODEL or FLUXCAL_SIM

  int   OPT_OBS ;
  bool  SNTABLE_NEVT ;

  char  COMMENT_FLUXREF[40];
} INPUTS ;


void  parse_args(int argc, char **argv) ;
void  open_fitresFile(void);
void  set_outlier_varnames(void); 
void  write_IGNORE_FILE(void) ;
void  write_headerInfo(FILE *FP) ;
bool  keyMatch_dash(char *arg, char *key_base);

FILE *FP_OUTFILE ;
char LINEKEY_DUMP[40];  // 'SN:' or  ''
char SEPKEY_DUMP[2];   // ',' or ' '

// ==================================
int main(int argc, char **argv) {

  char *TFILE, *TID, *TLIST[MXVAR_DUMP] ;
  int   NVAR, ivar, NDUMP, DO_IGNORE, IVAR_NPT ;

  // ------------- BEGIN -----------

  printf(" Begin sntable_dump.exe \n"); fflush(stdout);

  parse_args(argc,argv);

  DO_IGNORE = ( strcmp(INPUTS.FORMAT_OUTFILE,"IGNORE") == 0 ) ;

  set_FILTERSTRING(FILTERSTRING);
  set_EXIT_ERRCODE(EXIT_ERRCODE_sntable_dump);

  set_outlier_varnames(); // only if option is set

  // open fitres file if VARNAMES is specified.
  open_fitresFile();

  // set local variables 
  TFILE    = INPUTS.TABLE_FILE ;
  TID      = INPUTS.TABLE_ID ;
  NVAR     = INPUTS.NVAR ;
  IVAR_NPT = INPUTS.IVAR_NPT ;

  for(ivar=0; ivar<NVAR; ivar++ )   
    {  TLIST[ivar] = INPUTS.VARNAMES[ivar] ;     }


  TABLEFILE_INIT();

  if ( strlen(TID) == 0 ) { 
    // no table id/name -> print list of ntuples/trees, then abort.
    SNTABLE_LIST(TFILE); 
    goto DONE ;
  }


  if ( INPUTS.SNTABLE_NEVT ) {
    int NEVT = SNTABLE_NEVT(TFILE,TID);  // Aug 2020
    printf(" NEVT:  %d \n", NEVT);   // script-parsable output
  }
  else if ( NVAR == 0 ) {
    // no input variables --> list variable names
    SNTABLE_DUMP_VARNAMES(TFILE,TID);  
  }
  else if ( INPUTS.OUTLIER_NSIGMA[0] >= 0.0 ) {
    // dump fit-outliers for each epoch to ascii/fitres file
    // Must set SNTABLE_LIST = 'FITRES+RESIDUALS'
    NDUMP = SNTABLE_DUMP_OUTLIERS(TFILE, TID, NVAR, TLIST, IVAR_NPT, 
				  INPUTS.OUTLIER_NSIGMA, FP_OUTFILE,
				  LINEKEY_DUMP, SEPKEY_DUMP );

    if ( DO_IGNORE ) { write_IGNORE_FILE(); }
  } 
  else {
    // dump variable values to ascii/fitres file.
    NDUMP = SNTABLE_DUMP_VALUES(TFILE, TID, NVAR, TLIST, IVAR_NPT,
				FP_OUTFILE, LINEKEY_DUMP, SEPKEY_DUMP );  
  }

  if ( NVAR > 0 ) {
    printf("\n OUTFILE_DUMP: %s \n",  
	   INPUTS.OUTFILE_FITRES );
    
    if ( DO_IGNORE )
      { printf(" OUTFILE_IGNORE: %s \n", INPUTS.OUTFILE_IGNORE); }

    fflush(stdout);
  }

 DONE:
  return(0);  // return value for sntable_dump.pl script

} // end of main



// ============================
void parse_args(int NARG, char **argv) {

  // Feb 2 2018: parse OBS argument

  int i, NVAR, IFLAG_VARNAMES ;
  char fnam[] = "parse_args" ;

  // -------------- BEGIN -------------

  INPUTS.NVAR         =  0 ;
  INPUTS.ADD_HEADER   =  1 ;
  INPUTS.OUTLIER_NSIGMA[0]  = -9.0 ;
  INPUTS.OUTLIER_NSIGMA[1]  = -9.0 ;
  INPUTS.OPT_OBS      = 0 ;
  INPUTS.SNTABLE_NEVT = false;

  sprintf(LINEKEY_DUMP, "SN:");
  SEPKEY_DUMP[0] = 0;
  
  INPUTS.TABLE_FILE[0]     = 0;
  INPUTS.TABLE_ID[0]       = 0;
  INPUTS.IS_SNANA = INPUTS.IS_FITRES = false ;

  INPUTS.FORMAT_OUTFILE[0] = 0;
  INPUTS.VARLIST[0]        = 0;
  sprintf(INPUTS.OUTFILE_FITRES, "sntable_dump.fitres" );
  sprintf(INPUTS.OUTFILE_IGNORE, "sntable_dump.ignore" );
  sprintf(INPUTS.OUTLIER_VARNAME_CHI2FLUX, "NULL_CHI2FLUX" );

  INPUTS.ISFORMAT_CSV = 0 ;

  IFLAG_VARNAMES   =  0 ; 

  if ( NARG <= 1 ) {
    sprintf(msgerr1,"At least one argument required.") ;
    sprintf(msgerr2,"USAGE: sntable_dump.exe <fileName>  <tableName>" );
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }

  sprintf(INPUTS.TABLE_FILE, "%s", argv[1] );  // required

  if ( NARG > 2 )  { 
    sprintf(INPUTS.TABLE_ID, "%s", argv[2] ); 

    if (strcmp(INPUTS.TABLE_ID,TABLENAME_SNANA) == 0 ) 
      { INPUTS.IS_SNANA = true; }

    if (atoi(INPUTS.TABLE_ID) == TABLEID_SNANA )
      { INPUTS.IS_SNANA = true; }

    if (strcmp(INPUTS.TABLE_ID,TABLENAME_FITRES) == 0 ) 
      { INPUTS.IS_FITRES = true; }

    if (atoi(INPUTS.TABLE_ID) == TABLEID_FITRES )
      { INPUTS.IS_FITRES = true; }
  }

  // ----

  for(i=3; i<NARG; i++ ) {

    if ( strcmp(argv[i],"NOHEADER" ) == 0 )  { 
      IFLAG_VARNAMES    = 0;
      INPUTS.ADD_HEADER = 0  ;   
      LINEKEY_DUMP[0] = 0;
    }


    if ( strcmp_ignoreCase(argv[i],"nevt" ) == 0 ) {
      INPUTS.SNTABLE_NEVT = true ;  // Aug 2020
    }

    // xxx mark    if ( strcmp_ignoreCase(argv[i],"-outlier" ) == 0 ) {
    if ( keyMatch_dash(argv[i],"outlier") || 
	 keyMatch_dash(argv[i],"outlier_fit") ) {
      sscanf(argv[i+1], "%f", &INPUTS.OUTLIER_NSIGMA[0] );
      sscanf(argv[i+2], "%f", &INPUTS.OUTLIER_NSIGMA[1] );
      sprintf(INPUTS.OUTLIER_VARNAME_CHI2FLUX, "CHI2FLUX" );
      sprintf(INPUTS.OUTLIER_VARNAME_FLUXCAL,  "FLUXCAL_MODEL" );
      sprintf(INPUTS.COMMENT_FLUXREF,          "fitFlux" );
    }

    // xxx mark if ( strcmp_ignoreCase(argv[i],"-outlier_sim" ) == 0 ) {
    if ( keyMatch_dash(argv[i],"outlier_sim") ) {
      sscanf(argv[i+1], "%f", &INPUTS.OUTLIER_NSIGMA[0] );
      sscanf(argv[i+2], "%f", &INPUTS.OUTLIER_NSIGMA[1] );
      sprintf(INPUTS.OUTLIER_VARNAME_CHI2FLUX, "CHI2FLUX_SIM" );
      sprintf(INPUTS.OUTLIER_VARNAME_FLUXCAL,  "FLUXCAL_SIM" );
      sprintf(INPUTS.COMMENT_FLUXREF,          "simFlux" );
    }

    if ( strcmp_ignoreCase(argv[i],"obs" ) == 0 ) {
      INPUTS.OPT_OBS      = 1; 
      INPUTS.OUTLIER_NSIGMA[0] = 0.0 ;
      INPUTS.OUTLIER_NSIGMA[1] = 1.0E8 ;
      sprintf(INPUTS.OUTLIER_VARNAME_CHI2FLUX, "CHI2FLUX_SIM" ) ;
      sprintf(INPUTS.OUTLIER_VARNAME_FLUXCAL,  "FLUXCAL_SIM"  ) ;
      sprintf(INPUTS.COMMENT_FLUXREF,          "simFlux"      ) ;
    }

    if ( strcmp(argv[i],"-O" ) == 0  || strcmp(argv[i],"-o")== 0 )  { 
      IFLAG_VARNAMES = 0;
      sprintf(INPUTS.OUTFILE_FITRES,"%s", argv[i+1])  ; 
    }
    
    if ( strcmp_ignoreCase(argv[i],"-format" ) == 0 ) {
      sscanf(argv[i+1], "%s", INPUTS.FORMAT_OUTFILE );
      IFLAG_VARNAMES = 0;
      if ( strcmp(INPUTS.FORMAT_OUTFILE,"csv") == 0 ||
	   strcmp(INPUTS.FORMAT_OUTFILE,"CSV") == 0 )  { 
	INPUTS.ISFORMAT_CSV = 1; 
	LINEKEY_DUMP[0] = 0 ;
	SEPKEY_DUMP[0]  = 0 ;
      }
    }


    if ( IFLAG_VARNAMES ) {
      NVAR = INPUTS.NVAR ;
      sprintf(INPUTS.VARNAMES[NVAR],"%s", argv[i] );
      // buggy if ( NVAR==0 ) { sprintf(INPUTS.VARNAMES[NVAR],"CCID ROW"); }
      if ( NVAR==0  && strcmp(argv[i],"CID")!=0 ) 
	{ sprintf(INPUTS.VARNAMES[NVAR],"CCID"); } // Aug 17 2017
      sprintf(INPUTS.VARLIST,"%s %s", INPUTS.VARLIST, argv[i] );
      INPUTS.NVAR++ ; 
    }

    if ( strcmp(argv[i],"-V" ) == 0 )  { IFLAG_VARNAMES = 1; }
    if ( strcmp(argv[i],"-v" ) == 0 )  { IFLAG_VARNAMES = 1; }
    
  }


  // print summary of inputs

  printf(" Table File : %s \n", INPUTS.TABLE_FILE );

  if ( strlen(INPUTS.TABLE_ID) > 0  ) 
    {  printf(" Table ID   : %s \n", INPUTS.TABLE_ID ); }

  for(i=0; i < INPUTS.NVAR; i++ ) 
    { printf(" Table VARNAME(%d) : %s \n", i, INPUTS.VARNAMES[i] );  }


  fflush(stdout);
  return ;

} // end of parse_args


// ==================================
bool keyMatch_dash(char *arg, char *key_base) {

  // Created march 2021
  // for input *arg, check if it matches
  //    key_base
  //   -key_base
  //  --key_base

  int i;
  char key[60];
  char dash_list[3][4] = { "", "-", "--" };

  // ---------- BEGIN -----------

  for(i=0; i < 3; i++ ) {
    sprintf(key, "%s%s", dash_list[i], key_base);
    if ( strcmp_ignoreCase(arg,key) == 0 ) { return true; }
  }

  return false;

} // end keyMatch


// ==================================
void  set_outlier_varnames(void) {

  // Created Aug 2 2014 by R.Kessler
  // for the OUTLIER option, set a hard-wired list of variables
  //
  // Feb 22 2018: add more variables for OPTOBS
  // Nov 01 2019: add FLUXERRCALC_SIM for 'OBS' option
  // Nov 30 2020: check IS_SNANA and IS_FITRES

  int  ivar, NVAR = 0, NVAR_APPEND = 0 ;
  char VARNAMES_APPEND[MXVAR_DUMP][60];
  char fnam[] = "set_outlier_varnames";

  // ----------- BEGIN -------

  if ( INPUTS.OUTLIER_NSIGMA[0] < -0.0001 ) { return ; }

  // if --v option gave additional variables, then store
  // them here and append them below. The user variables
  // are stored in INPUTS.VARNAMES, but below we need to
  // over-write this array with a hard-wired 'outlier' list.

  NVAR_APPEND = INPUTS.NVAR ;
  for(ivar=0; ivar < NVAR_APPEND; ivar++ ) {
    sprintf(VARNAMES_APPEND[ivar], "%s", INPUTS.VARNAMES[ivar] );
  }

  // First three elements correspond to optional IGNORE format  
  sprintf(INPUTS.VARNAMES[NVAR],"%s", "CCID" );   
  NVAR++ ;

  // - - - - -
  sprintf(INPUTS.VARNAMES[NVAR],"%s", "FIELD" );   
  NVAR++ ;
  
  sprintf(INPUTS.VARNAMES[NVAR],"%s", "zHEL" );   
  NVAR++ ;

  if ( INPUTS.OPT_OBS ) {
    sprintf(INPUTS.VARNAMES[NVAR],"%s", "zHD"     );    NVAR++ ;
    if ( INPUTS.IS_FITRES ) 
      { sprintf(INPUTS.VARNAMES[NVAR],"%s", "FITPROB" );    NVAR++ ; }
  }

  // append extra user variables here
  for(ivar=0; ivar < NVAR_APPEND; ivar++ ) {
    sprintf(INPUTS.VARNAMES[NVAR], "%s", VARNAMES_APPEND[ivar] );
    NVAR++ ;
  }

  /* xxxxxxxxxx
  if ( INPUTS.IS_SNANA )
    { sprintf(INPUTS.VARNAMES[NVAR],"%s", VARNAME_CUTFLAG_SNANA ); NVAR++ ; }
  xxxxxx */

  if ( INPUTS.IS_SNANA ) 
    { sprintf(INPUTS.VARNAMES[NVAR],"%s", "NOBS" );  }    
  else 
    { sprintf(INPUTS.VARNAMES[NVAR],"%s", "NPTFIT" );  }
  INPUTS.IVAR_NPT = NVAR;
  NVAR++ ;

  // epoch/vector quantities are after NPTFIT

  sprintf(INPUTS.VARNAMES[NVAR],"%s", "MJD" );   
  NVAR++ ;


  sprintf(INPUTS.VARNAMES[NVAR],"%s", "IFILTOBS" );   
  NVAR++ ;

  if ( INPUTS.IS_FITRES ) {
    sprintf(INPUTS.VARNAMES[NVAR],"%s", "REJECT" );  // 1--> excluded from fit
    NVAR++ ;

    sprintf(INPUTS.VARNAMES[NVAR],"%s", "TOBS" );   
    NVAR++ ;
  }

  sprintf(INPUTS.VARNAMES[NVAR],"%s", "PSF" );   
  NVAR++ ;

  sprintf(INPUTS.VARNAMES[NVAR],"%s", "ZP" );   
  NVAR++ ;

  sprintf(INPUTS.VARNAMES[NVAR],"%s", "GAIN" );   // Feb 27 2018
  NVAR++ ;

  //  sprintf(INPUTS.VARNAMES[NVAR],"%s", "SKY_SIG" );  
  sprintf(INPUTS.VARNAMES[NVAR],"%s", "SKYSIG" );   
  NVAR++ ;

  sprintf(INPUTS.VARNAMES[NVAR],"%s", "SKYSIG_T" );   
  NVAR++ ;

  sprintf(INPUTS.VARNAMES[NVAR],"%s", "FLUXCAL_DATA" );   
  NVAR++ ;

  // output either FLUXCAL_MODEL or FLUXCAL_SIM
  sprintf(INPUTS.VARNAMES[NVAR],"%s", INPUTS.OUTLIER_VARNAME_FLUXCAL );   
  NVAR++ ;

  sprintf(INPUTS.VARNAMES[NVAR], "FLUXCAL_DATA_ERR" );  
  NVAR++ ;
    
  sprintf(INPUTS.VARNAMES[NVAR],"%s", INPUTS.OUTLIER_VARNAME_CHI2FLUX );   
  NVAR++ ;

  if ( INPUTS.OPT_OBS  &&  INPUTS.IS_FITRES ) {
    sprintf(INPUTS.VARNAMES[NVAR], "FLUXCAL_MODEL" );    NVAR++ ; 
    sprintf(INPUTS.VARNAMES[NVAR], "SBFLUXCAL" );         NVAR++ ;
    sprintf(INPUTS.VARNAMES[NVAR], "ERRTEST"   );         NVAR++ ;
    sprintf(INPUTS.VARNAMES[NVAR], "FLUXERRCALC_SIM" );   NVAR++ ;
  }

  if ( NVAR >= MXVAR_DUMP ) {
    sprintf(msgerr1,"NVAR=%d exceeds bound of MXVAR_DUMP=%d",
	    NVAR, MXVAR_DUMP );
    sprintf(msgerr2,"Check MXVAR_DUMP in sntable_dump.c");
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
    
  }

  INPUTS.NVAR = NVAR ;

  return ;

} // end of set_outlier_varnames


// ==================================
void  open_fitresFile(void) {

  // open ascii file to write.
  // Use global file pointer FP_OUTFILE.
  // Aug 4 2017: for SIMLIB table, fisrt colum is ROW, not CID.

  int i ;
  char VARNAME_CID[] = "CID";
  char VARNAME_ROW[] = "ROW";
  char SEP[4] = " " ;
  char *ptrVar;
  //  char fnam[] = "open_fitresFile" ;
  // -------------- BEGIN --------------

  FP_OUTFILE = NULL ;

  if ( INPUTS.NVAR == 0 ) { return ; }

  if ( strcmp(INPUTS.OUTFILE_FITRES,"stdout") == 0 ) 
    { FP_OUTFILE = stdout ; }
  else 
    { FP_OUTFILE = fopen(INPUTS.OUTFILE_FITRES, "wt") ; }

  if ( INPUTS.ADD_HEADER == 0 ) { return ; }

  write_headerInfo(FP_OUTFILE); 

#ifdef TEXTFILE_NVAR
  if ( INPUTS.ISFORMAT_CSV == 0 ) 
    { fprintf(FP_OUTFILE,"NVAR: %d\n", INPUTS.NVAR ); } 
#endif

  if ( INPUTS.ISFORMAT_CSV ) 
    { sprintf(SEP,", "); }
  else
    { fprintf(FP_OUTFILE,"VARNAMES:  ");  }


  for ( i=0; i < INPUTS.NVAR; i++ )  { 
    
    ptrVar = INPUTS.VARNAMES[i] ;  // default

    // for text file, only use CID as varname (CCID for hbook & root)
    if ( i==0 ) {
      if ( strcmp(INPUTS.VARNAMES[i],"CCID") == 0 ) 
	{ ptrVar = VARNAME_CID ; }  // replace CCID with CID in text VARNAMES
      else if ( strcmp(INPUTS.TABLE_ID,"SIMLIB") == 0 ) 
	{ ptrVar = VARNAME_ROW ; }  // use ROW for 1st SIMLIB column 
    }

    // Mar 14 2021: insert band string before IFILTOBS
    if ( ISTABLEVAR_IFILT(ptrVar) )
      { fprintf(FP_OUTFILE,"%s%s", "BAND", SEP );  }

    if ( i == INPUTS.NVAR-1 ) { sprintf(SEP," ") ; }
    fprintf(FP_OUTFILE,"%s%s", ptrVar, SEP ); 
  }


  fprintf(FP_OUTFILE,"\n");

  return ;

} // end of  open_fitresFile


// ==================================
void write_headerInfo(FILE *FP) {

  if ( INPUTS.ISFORMAT_CSV ) { return ; }

  fprintf(FP,"# Variables are extracted from \n");
  fprintf(FP,"#   FILE:  %s \n", INPUTS.TABLE_FILE );
  fprintf(FP,"#   TABLE: %s \n", INPUTS.TABLE_ID   );
  fprintf(FP,"# \n" );

  float SIG0 = INPUTS.OUTLIER_NSIGMA[0] ;
  float SIG1 = INPUTS.OUTLIER_NSIGMA[1] ;

  if ( SIG0 >= 0.0 && INPUTS.OPT_OBS==0 ) {
    fprintf(FP,"# Select epoch fit-outliers defined by \n" );
    fprintf(FP,"#    %.1f < |(dataFlux - %s)/err| < %.1f \n",
	    SIG0, INPUTS.COMMENT_FLUXREF, SIG1 );
    fprintf(FP,"# --> \n" ) ;
    fprintf(FP,"#    %.1f < %s < %.1f \n",
	    SIG0*SIG0, INPUTS.OUTLIER_VARNAME_CHI2FLUX, SIG1*SIG1 );
    fprintf(FP,"# \n" );
  }

  //  fprintf(FP,"# REJECT=1 --> excluded from fit.\n" );
  //  fprintf(FP,"# \n" );

} // end write_headerInfo

// ==================================
void write_IGNORE_FILE(void) {

  // Created July 22 2017
  //
  // Read newly-created OUTLIER file (created by this program), 
  // and create new file in INGORE format for SNANA:
  //    IGNORE:  <SNID>  <MJD>  <BAND>  # misc
  //
  // This file should be copied to
  //   $SNDATA_ROOT/lcmerge/[VERSION]/[VERSION].IGNORE
  //

  FILE *FP ;
  int  OPTMASK=3, NROW, irow, ISTAT, IFILTOBS, REJECT ;
  double DVAL, MJD, PSF, TOBS, CHI2 ;
  char VARLIST[200], CCID[32], CCID_LAST[32], CVAL[32], BAND[2];
  char *VARNAME_CHI2 = INPUTS.OUTLIER_VARNAME_CHI2FLUX ;
  char TBLNAME[] = "OUTLIERS" ; 
  //  char fnam[] = "write_IGNORE_FILE" ;

  // ---------------- BEGIN ---------------

  sprintf(BANNER,"Write IGNORE format to %s", INPUTS.OUTFILE_IGNORE);
  print_banner(BANNER);

  // read newly-created table
  sprintf(VARLIST,"CID,MJD,IFILTOBS,REJECT,TOBS,PSF,%s", VARNAME_CHI2);
  NROW = 
    SNTABLE_AUTOSTORE_INIT(INPUTS.OUTFILE_FITRES, TBLNAME, VARLIST, OPTMASK );


  // ------ open IGNORE file ------
  FP = fopen(INPUTS.OUTFILE_IGNORE, "wt") ;
  CCID_LAST[0] = 0 ;

  write_headerInfo(FP); 

  for(irow=0; irow < NROW; irow++ ) {

    sprintf(CCID, "%s", SNTABLE_AUTOSTORE[0].CCID[irow] );

    // add blank space between CIDs to visually identify
    // CIDs with multiple outliers.
    if ( strcmp(CCID,CCID_LAST) != 0 ) 
      { fprintf(FP,"\n"); }

    SNTABLE_AUTOSTORE_READ(CCID, "MJD",      &ISTAT, &DVAL, CVAL ); 
    MJD = DVAL ;

    SNTABLE_AUTOSTORE_READ(CCID, "IFILTOBS", &ISTAT, &DVAL, CVAL ); 
    IFILTOBS = (int)DVAL ;

    SNTABLE_AUTOSTORE_READ(CCID, "REJECT", &ISTAT, &DVAL, CVAL ); 
    REJECT = (int)DVAL ;

    SNTABLE_AUTOSTORE_READ(CCID, "PSF", &ISTAT, &DVAL, CVAL ); 
    PSF = DVAL ;

    SNTABLE_AUTOSTORE_READ(CCID, "TOBS", &ISTAT, &DVAL, CVAL ); 
    TOBS = DVAL ;

    SNTABLE_AUTOSTORE_READ(CCID, VARNAME_CHI2, &ISTAT, &DVAL, CVAL ); 
    CHI2 = DVAL ;

    sprintf(BAND, "%c", FILTERSTRING[IFILTOBS] );
    fprintf(FP, "IGNORE: %s %.3f %s    "
	    "# PSF=%.2f  Tobs=%6.1f  CHI2=%5.1f  REJ=%d\n", 
	    CCID, MJD, BAND, PSF, TOBS, CHI2, REJECT ) ;  

    sprintf(CCID_LAST,"%s", CCID);
  }

  fclose(FP);

  return ;

} // end write_IGNORE_FILE
