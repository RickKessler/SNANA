/*********************************
  Jun 23, 2009 by R.Kessler
 
  ------------
  
  Combine multiple fitres contents into a single table in 
  TEXT/ROOT/HBOOK format. Designed to study correlations
  between different fitres outputs for the same set of SN.

 Usage:
  >  combine_fitres.exe <fitres1> <fitres2> ...

   where each command-line argument is the name of a
   fitres file to combine into the ntuple.
   Output files are
      combine_fitres.hbook  (hbook column-wise ntuple)
      combine_fitres.root   (root tree)
      combine_fitres.txt    (combined fitres text-file)

 Options below accept python like args, but with no dashs, 1 dash
 or 2 dashes. I.e., outprefix, -outprefix, --outprefix are all 
 accepted.

  >  combine_fitres.exe <fitres1> <fitres2> --outprefix <outprefix>
    produces output files <outprefix>.hbook and  <outprefix>.text

  >  combine_fitres.exe <fitres1> <fitres2> --outfile_text <outfile>
       (text outpout name is <outfile>)

  >  combine_fitres.exe <fitres1> <fitres2> --outfile_text <outfile.gz>
       (produce gzipped outfile.gz)

  >  combine_fitres.exe <fitres1>   R      ! .ROOT extension 
  >  combine_fitres.exe <fitres1>  -R      ! same with python syntax
  >  combine_fitres.exe <fitres1>   r      ! .root extension
       (create root file instead of hbook file)

  >  combine_fitres.exe <fitres1>  H      ! .HBOOK extension
  >  combine_fitres.exe <fitres1>  h      ! .hbook extension
       (create default hbook file ; same as with no H or h arg)

  >  combine_fitres.exe <fitres1>  R H 
  >  combine_fitres.exe <fitres1>  H R
      (create both hbook and root files)

  >  combine_fitres.exe <fitres1>  --mxrow 50000

  >  combine_fitres.exe <fitres1>  --varnames zHD,c,x1,SIM_DLMAG
          [select only these varnames (along with CID)]

  >  combine_fitres.exe  <fitres1> --zcut <zmin> <zmax>
         [cut on zHD]

  >  combine_fitres.exe <fitres1> <fitres2> ..  T   ! [outPrefix].TEXT
  >  combine_fitres.exe <fitres1> <fitres2> .. -T   ! same in python syntax
  >  combine_fitres.exe <fitres1> <fitres2> ..  t   ! [outPrefix].text
      (create only text output; disable default hbook output)

  >  combine_fitres.exe  <fitres1> <fitres2> ... --nullval_float -12345
         (override default nullval_float = -888)

 WARNINGS/NOTES:
 * If fitres files contain different SN, then first
   fitres file determines the list of SN; extra SN
   in the 2nd, 3rd ... fitres file are ignored.


                HISTORY
 
  Mar 31 2016: MXFILE -> 20 ( was 10)

  Dec 04 2016: allow VERSION column

  Feb 27 2017: add check for CCID being string or integer;
               CIDint -> isn if CCID is a string.

  Mar 8 2017: 
    make VARNAME_1ONLY list of SNANA variables to include only once. 
    Also check for first SNANA file (IFILE_FIRST_SNANA) in case 
    first file was not made by SNANA.

  Apr 29 2019: option T or t for text-only (no hbook)
  May 03 2019: refactor indices to start at 0 instead of 1

  Jun 06 2019:
    +  new argument -outfile_text <outFile>
    +  if outFile has .gz extension, gzip it.

  Sep 19 2019
   + naive trick to speed up CID-matching, but works only
     if both FITRES files have same list of CIDs. 
   + MXSTRLEN -> MXCHAR_VARNAME (instead of 28)

 Oct 7 2019:
   use hash table for cid matching using uthash.h from
      http://troydhanson.github.io/uthash/index.html
   It's way way faster for large (e.g., biasCor) files.
   To revert back to slow matching method,
      combine_fitres.exe <argList>  -matchflag 1
   Beware of significant refactor.

  Jan 16 2020: 
    + implement -mxrow (was read, but not implemented)
    + new input -zcut <zmin> <zmax>
  Jan 23 2020
    + abort if NEVT_APPROX > MXSN
    + USEDCID -> bool instead of short int.
    + MXSN -> 5M (was 2M)

 Apr 27 2020:
   + init strings to NULL (see INIVAL_COMBINE_STR)

 Jun 7 2020: 
   + match variables with ignore-case so that zCMB and ZCMB are
     written out as zCMB and ZCMB_2.

 Sep 15 2020: 
   + new input -nullval_float 0 (override default -888). Initial use is
     for Pippin to set classifier PROB_CC=0 when classifier returns
     no value.
   + allow 1 or 2 dashes in front of input args to allow pythonic structure.

 Oct 27 2020: fix a few warnings from -Wall flag 

 Nov 18 2020: allow comma-sep list of space-sep list of fitres files.
              See new function parse_FFILE(arg).

 Dec 09 2020: new -varnames arg to select subset of variables to write out.

 Jan 4 2021: add print_stats() for NEVT, NEVT_MISSING, NEVT_COMMON

 Jun 08 2021: switch to using match_cid_hash  utility in sntools.c
                (matchflag=5).
   
 Apr 27 2022: 
   + print full command to stdout.
   + new parsing util keyarg_match to allow args with no dashes, 
     1 dash, or 2 dashes ... to allow python-like inputs.
     All previous args should still work; this update allows
     more variations of inputs.

******************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "sntools.h"
#include "sntools_output.h" 
#include "uthash.h" 

// ================================
// Function Declarations
// ================================


void  PARSE_ARGV(int argc, char **argv);
bool  keyarg_match(char *arg, char *keyname_base);
void  parse_FFILE(char *arg);

void  init_misc(void);
void  INIT_TABLEVAR(void);
void  ADD_FITRES(int ifile);
int   match_CID_orig(int ifile, int isn2);
int   match_CID_hash_local(int ifile, int isn2);
void  ADD_FITRES_VARLIST(int ifile, int isn, int isn2);

int   NMATCH_VARNAME(char *ctag , int ntlist ) ;
int   maxlen_varString(char *varName);
int   SKIP_VARNAME(int file, int ivar) ;

void  print_stats(void);

void  WRITE_SNTABLE(void); // output table in selected format
void  outFile_text_override(char *outFile, int *GZIPFLAG);
void  ADD_SNTABLE_COMMENTS(void) ;

void  fitres_malloc_flt(int ifile, int NVAR, int MAXLEN); 
void  fitres_malloc_str(int ifile, int NVAR, int MAXLEN); 
void  freeVar_TMP(int ifile, int NVARTOT, int NVARSTR, int MAXLEN); 

void malloc_NMATCH_PER_EVT(int N);

// ================================
// Global variables
// ================================


#define MXFFILE  20       // max number of fitres files to combine
#define MXSN     10000000   // max SN to read per fitres file
#define MXVAR_PERFILE  50  // max number of NTUP variables per file
#define MXVAR_TOT  MXVAR_TABLE     // max number of combined NTUP variables

#define DEFAULT_NULLVAL_FLOAT  -888.0 // default float value
#define DEFAULT_NULLVAL_STRING  "NULL" // default string value

#define MXSTRLEN       MXCHAR_VARNAME  // changed from 28 (Sep 20 2019)
#define MXSTRLEN_BAND      4
#define MXSTRLEN_CID      20
#define MXSTRLEN_FIELD    20
#define MXSTRLEN_VERSION  32 // photometry versoin, added Dec 2016

#define IVARSTR_CCID  0   // CCID index for CVAR_XXX arrays

// logicals to control which output files to create
int CREATEFILE_HBOOK, CREATEFILE_ROOT, CREATEFILE_TEXT ;

int NVARALL_FILE[MXFFILE]; // Number of variables per file
int NVARSTR_FILE[MXFFILE]; // number of variables that are string

// cast for each variable; negative -> do not write
// See ICAST_FITRES[D,F,I,C] parameters in sntools.h
int ICAST_FITRES_COMBINE[MXVAR_TOT]; 

int NWRITE_SNTABLE ;

// for SNtable
#define  TABLEID_COMBINE   TABLEID_FITRES
#define  TABLENAME_COMBINE TABLENAME_FITRES

#define MATCHFLAG_ORIG       1 // original slow CID-matching
#define MATCHFLAG_HASH_LOCAL 3 // use hash table recommended by Sam
#define MATCHFLAG_HASH_UTIL  5 // generical utility in sntools (Jun 2021)

// inputs
struct INPUTS {
  int  NFFILE ;   // number of input fitres files to combine
  char FFILE[MXFFILE][2*MXPATHLEN] ;
  char OUTPREFIX_COMBINE[MXPATHLEN] ;
  char OUTFILE_TEXT[MXPATHLEN] ;
  int  MXROW_READ ;
  int  MATCHFLAG ;
  double CUTWIN_zHD[2];
  int    DOzCUT;
  float  NULLVAL_FLOAT ; // Sep 2020

  int  NVARNAMES_KEEP ;
  char **VARNAMES_KEEP; // select and save only these varnames
  char VARLIST_KEEP[MXPATHLEN];
} INPUTS ;


#define  NVARNAME_1ONLY 4  // NVAR to include only once
char VARNAME_1ONLY[NVARNAME_1ONLY][20] = 
  { "FIELD", "CIDint", "VERSION", "IDSURVEY" } ;
int IFILE_FIRST_SNANA; // first SNANA input file

char  SUFFIX_HBOOK[8] =   "HBOOK" ;
char  suffix_hbook[8] =   "hbook" ;
char  SUFFIX_ROOT[8]  =   "ROOT"  ;
char  suffix_root[8]  =   "root"  ;
char  SUFFIX_TEXT[8]  =   "TEXT"  ;
char  suffix_text[8]  =   "text" ;

char *ptrSuffix_hbook = suffix_hbook ; // default is lower case, unless H
char *ptrSuffix_root  = suffix_root ;
char *ptrSuffix_text  = suffix_text ;
bool USEDCID[MXSN];

int IVARSTR_STORE[MXVAR_TOT] ; // keep track of string vars
int IVAR_zHD;

int NLIST_FIRST_FITRES ;   // number of SN in 1st fitres file
int NLIST2_FITRES ;        // idem for 2nd, 3rd, etc ...

//  table variables
int   NVAR_WRITE_COMBINED ;  // total number of variables in combined file
int   NVARALL_FITRES ;       // idem, but including multiple CIDs
int   NVARSTR_FITRES ;

int   NVARALL_FITRES_LAST ;    
int   NVARSTR_FITRES_LAST ;

struct FITRES_VALUES {

  char  ***STR_ALL ;  // string values for all files [ivar][isn]
  char  ***STR_TMP ;  // string values for current file

  float  **FLT_ALL ;  // float values for all files  [ivar][isn]
  float  **FLT_TMP ;  // idem for 1 file

} FITRES_VALUES ;

char  VARNAME_COMBINE[MXVAR_TOT][MXSTRLEN]; // name of each combined variable


struct ROW_VALUES {
  float FLT[MXVAR_TOT] ;  
  char  STR[MXVAR_TOT][MXSTRLEN] ;
  char  CCID[MXSTRLEN_CID];     // added Dec 8 2014
  int   CIDint ;
} TABLEROW_VALUES ;       // values for 1 row; used for SNTABLE


struct hash_table {
  int id;                    /* key */
  char name[10];
  UT_hash_handle hh;         /* makes this structure hashable */
} ;

struct hash_table *users = NULL; 

char *NMATCH_PER_EVT ; // just 1 byte each to save memory
int  NEVT_COMMON;   // number of events common to all FITRES files
int  NEVT_READ[MXFFILE];    // NEVT read from each file
int  NEVT_MISSING[MXFFILE]; // NEVT missing w.r.t. 1st FITRES file

// =========================================
int main(int argc, char **argv) {

  int ifile ;
  //  char fnam[] = "main" ;

  // ----------------- BEGIN --------

  sprintf(BANNER,"Begin execution of combine_fitres.exe \n" );
  print_banner(BANNER);

  set_EXIT_ERRCODE(EXIT_ERRCODE_combine_fitres);

  init_misc();

  PARSE_ARGV(argc,argv);

  INIT_TABLEVAR();

  print_banner("Begin Reading Fitres Files.\n");

  TABLEFILE_INIT(); // Oct 27 2014

  for ( ifile = 0; ifile < INPUTS.NFFILE; ifile++ ) {
    ADD_FITRES(ifile);
  }

  if ( INPUTS.MATCHFLAG == MATCHFLAG_HASH_LOCAL ) 
    { match_CID_hash_local(-1,0); } // remove hash table

  if ( INPUTS.MATCHFLAG == MATCHFLAG_HASH_UTIL ) 
    { match_cid_hash("", -1,0); } // remove hash table

  // ---------------

  WRITE_SNTABLE() ;

  printf("   Done writing %d events. \n", NWRITE_SNTABLE );
  fflush(stdout);

  print_stats();

  return(0);

} // end of main


// ===============================
void  init_misc(void) {

  int i;
  // set default output to hbook if both root and hbook are compiled
  // (specified by order of inits below)
#ifdef USE_ROOT  
  CREATEFILE_ROOT = 1 ;  CREATEFILE_HBOOK = 0 ;
#endif

#ifdef USE_HBOOK
  CREATEFILE_ROOT = 0 ;  CREATEFILE_HBOOK = 1 ;
#endif

#ifdef USE_TEXT
  CREATEFILE_TEXT = 1;
#endif

  NEVT_COMMON = 0 ;
  for(i=0; i < MXFFILE; i++ )  
    { NEVT_MISSING[i] = NEVT_READ[i] = 0; }

} // end init_misc

// ===============================
void  PARSE_ARGV(int argc, char **argv) {

  int  i, NFFILE=0 ;
  char fnam[] = "PARSE_ARGV" ;

  //----------- BEGIN --------------

  NARGV_LIST = argc ; // set global 

  INPUTS.NFFILE       = 0;
  INPUTS.MXROW_READ   = 1000000000 ;

  INPUTS.MATCHFLAG    = MATCHFLAG_HASH_UTIL ;  // Jun 2021  

  INPUTS.OUTFILE_TEXT[0]  = 0 ;
  INPUTS.CUTWIN_zHD[0] = -9.0 ;  
  INPUTS.CUTWIN_zHD[1] = +9.0 ; 
  INPUTS.DOzCUT = 0 ;
  sprintf(INPUTS.OUTPREFIX_COMBINE, "combine_fitres" );
  INPUTS.NULLVAL_FLOAT   =  DEFAULT_NULLVAL_FLOAT ;
  INPUTS.NVARNAMES_KEEP  = 0 ;
  INPUTS.VARLIST_KEEP[0] = 0 ;

  
  printf("\n Full command: ");
  for ( i = 0; i < NARGV_LIST ; i++ ) {  printf("%s ", argv[i]);   }
  printf("\n\n"); fflush(stdout);

  // - - - - -
  for ( i = 1; i < NARGV_LIST ; i++ ) {
    
    // check for optional args

    if ( keyarg_match(argv[i],"outprefix"))  {
      i++ ; sprintf(INPUTS.OUTPREFIX_COMBINE,"%s", argv[i]);
      continue ;
    }

    /* xxx mark delete xxxxx
    if ( strcmp(argv[i],"-outPrefix") == 0 ) { // allow Fermi-spell
      i++ ; sprintf(INPUTS.OUTPREFIX_COMBINE,"%s", argv[i]);
      continue ;
    }
    xxxxxxx */

    if ( keyarg_match(argv[i],"outfile_text"))  {
      i++ ; sprintf(INPUTS.OUTFILE_TEXT,"%s", argv[i]);
      continue ;
    }

    if ( keyarg_match(argv[i],"mxrow"))  {
      i++ ; sscanf(argv[i], "%d", &INPUTS.MXROW_READ);
      continue ;
    }

    if ( keyarg_match(argv[i],"varnames"))  {
      i++ ; sscanf(argv[i], "%s", INPUTS.VARLIST_KEEP);
      parse_commaSepList("VARNAMES_KEEP", INPUTS.VARLIST_KEEP, MXVAR_TOT,
			 40, &INPUTS.NVARNAMES_KEEP, &INPUTS.VARNAMES_KEEP);
      continue ;
    }

    if ( keyarg_match(argv[i],"zcut"))  {
      i++ ; sscanf(argv[i], "%le", &INPUTS.CUTWIN_zHD[0] );
      i++ ; sscanf(argv[i], "%le", &INPUTS.CUTWIN_zHD[1] );
      INPUTS.DOzCUT = 1; 
      continue ;
    }

    if ( keyarg_match(argv[i],"nullval_float"))  {
      i++ ; sscanf(argv[i], "%f", &INPUTS.NULLVAL_FLOAT);
      continue ;
    }

    if ( keyarg_match(argv[i],"r") )  {
      CREATEFILE_ROOT = 1;
      if (CREATEFILE_HBOOK==1) {CREATEFILE_HBOOK=0;}  // root on, hbook off
      if ( strstr(argv[i],"R")!=NULL ) { ptrSuffix_root = SUFFIX_ROOT ; }
      continue ;
    }

    if ( keyarg_match(argv[i],"h") )  {
      CREATEFILE_HBOOK = 2 ;  // turn hbook back on and leave it on
      if ( strstr(argv[i],"H") !=NULL ) { ptrSuffix_hbook = SUFFIX_HBOOK ; }
      continue ;
    }

    if ( keyarg_match(argv[i],"t") )  {
      CREATEFILE_HBOOK = 0 ; 
      CREATEFILE_ROOT  = 0 ; 
      if ( strstr(argv[i],"T") !=NULL ) { ptrSuffix_text = SUFFIX_TEXT ; }
      continue ;
    }

    if ( keyarg_match(argv[i],"matchflag") )  {
      i++ ; sscanf(argv[i], "%d", &INPUTS.MATCHFLAG);
      continue ;
    }

    // parse FITRES file(s) and add to INPUTS.FFILE list
    parse_FFILE(argv[i]);

  } // end loop over arg list
  
  
  if ( INPUTS.MATCHFLAG == MATCHFLAG_ORIG ) {
    printf("   CID-match method: brute force loop over each file.\n");
  }
  else if ( INPUTS.MATCHFLAG == MATCHFLAG_HASH_LOCAL ) {
    printf("   CID-match method: hash table with local util.\n");
  }
  else if ( INPUTS.MATCHFLAG == MATCHFLAG_HASH_UTIL ) {
    printf("   CID-match method: hash table with sntools util.\n");
  }

  if ( INPUTS.NFFILE <= 0 ) {
    sprintf(c1err, "Bad args. Must give fitres file(s)");
    sprintf(c2err, "  combine_fitres.exe <fitresFile List> ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // idiot checks for file type

#ifndef USE_ROOT
  if ( CREATEFILE_ROOT ) {
    sprintf(c1err, "Cannot create output ROOT file because");
    sprintf(c2err, "'#define ROOT' is not set in sntools_output.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
#endif


#ifndef USE_HBOOK
  if ( CREATEFILE_HBOOK ) {
    sprintf(c1err, "Cannot create output HBOOK file because");
    sprintf(c2err, "'#define HBOOK' is not set in sntools_output.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
#endif

  return;

} // end of PARSE_ARGV

// ===============================================
bool keyarg_match(char *arg, char *keyname_base) {

  // return true if
  // arg == keyname_base  OR
  // arg == -keyname_base OR
  // arg == --keyname_base
  //
  // This allows python-like args, and also allows
  // SNANA style without dashes.

#define NKEY_CHECK 3
  bool match = false;
  int  k;
  char key_list[NKEY_CHECK][80];
  char fnam[] = "keyarg_match";

  // ------------- BEGIN -------------

  sprintf(key_list[0],"%s", keyname_base);
  sprintf(key_list[1],"-%s", keyname_base);
  sprintf(key_list[2],"--%s", keyname_base);

  for(k=0; k < NKEY_CHECK; k++ ) {
    if ( strcmp_ignoreCase(arg,key_list[k]) == 0 ) { match = true; }
  }

  return match;

} // end keyarg_match

// ====================================
void parse_FFILE(char *arg) {

  // Created Nov 18 2020
  // add arg to INPUTS.FFILE list ... if arg is comma-sep list,
  // split each item and add each item to list. This allows
  // both space-separated and comma-separated file list.

  int ifile, nfile_add=0, NFFILE = INPUTS.NFFILE;
  char **file_list;
  char fnam[] = "parse_FFILE" ;
  
  // ----------- BEGIN -------------

  parse_commaSepList("FITRES_FILE_LIST", arg, 10, MXPATHLEN,
		     &nfile_add, &file_list );

  for(ifile=0; ifile < nfile_add; ifile++ ) {
    sprintf( INPUTS.FFILE[NFFILE], "%s", file_list[ifile] );
    printf("  Will combine fitres file: %s \n", 
	   INPUTS.FFILE[NFFILE] );
    NFFILE++ ;
    free(file_list[ifile]);
  }

  INPUTS.NFFILE = NFFILE ;
  free(file_list);

} // end parse_FFILE

// ==========================
void INIT_TABLEVAR(void) {

  NVAR_WRITE_COMBINED = 0 ;

  NVARALL_FITRES = 0 ;
  NVARSTR_FITRES = 0 ;

  NLIST_FIRST_FITRES = NLIST2_FITRES = 0 ;

  IFILE_FIRST_SNANA = -9;
  IVAR_zHD = -9;

} // end INIT_TABLEVAR


// ==============================
void ADD_FITRES(int ifile) {

  // Jan 2011: set CTAG[1] =  VARNAMES_FITRES[0] so that either
  //           CID or LIBID can be the first element.
  //
  // Dec 21, 2013: avoid too-small LEN_APPROX on 2nd file
  // Oct 27, 2014: use refactored SNTABLE_XXX functions.
  //               Beware that READTABLE_POINTERS uses C-like index from 0.
  //
  // Dec 1 2017: 
  //   + use new function SNTABLE_NEVT_APPROX_TEXT to check for .gz file.
  //
  // May 2 2019: 
  //   + remove redundant call to TABLEFILE_CLOSE
  // Sep 19 2019: abort if column 0 is not a valid column
  // Sep 24: slightly improve matching method so that it is very 
  //         fast when both files have exactly the same CIDs.
  //

  int  ivar, ivarstr, j, isn, isn2, NMATCH2 ;
  int  NVARALL, NVARSTR, NVAR, NTAG_DEJA, NLIST, ICAST ;
  int  index=-9, REPEATCID, NEVT_APPROX, IFILETYPE, iappend ;

  char 
    *VARNAME, VARNAME_F[MXCHAR_VARNAME], VARNAME_C[MXCHAR_VARNAME]
    ,*ptr_CTAG, ccid[60]
    ,fnam[] = "ADD_FITRES"
    ;

  // ----------- BEGIN -----------

  printf("\n --------------------- %s -------------------------\n", fnam );

  // open file & read header; use generic table name SNTABLE
  // since for ascii files the table name is not used.
  IFILETYPE = TABLEFILE_OPEN( INPUTS.FFILE[ifile], "read text" );
  NVARALL   = SNTABLE_READPREP(IFILETYPE,"SNTABLE");
  NVARALL_FILE[ifile] = NVARALL;

  // check if this is an SNANA file; mark first SNANA file
  if ( IFILE_FIRST_SNANA < 0 ) {
    for ( ivar=0; ivar < NVARALL; ivar++ ) {    
      VARNAME = READTABLE_POINTERS.VARNAME[ivar] ;
      for(j=0; j < NVARNAME_1ONLY; j++ ) {
	if ( strcmp(VARNAME,VARNAME_1ONLY[j])==0 ) 
	  {  IFILE_FIRST_SNANA = ifile; }	
      }
    }
    if ( IFILE_FIRST_SNANA >= 0 ) {
      printf("\t First SNANA ifile = %d (%s) \n",
	     IFILE_FIRST_SNANA, INPUTS.FFILE[ifile] );
    }
  }

  // --------------------------------------------------
  //                 MEMORY ALLOC
  //
  // count how many string-variables: needed to allocate memory
  NVARSTR = 0 ;
  for ( ivar=0; ivar < NVARALL; ivar++ ) {    
    if ( SKIP_VARNAME(ifile, ivar) ) { continue ; }
    if ( READTABLE_POINTERS.ICAST_STORE[ivar] == ICAST_C ) 
      { NVARSTR++ ; }
  }

  // get number of SN for memory allocation
  NEVT_APPROX = SNTABLE_NEVT(INPUTS.FFILE[ifile],"TABLE");

  if ( NEVT_APPROX >= MXSN-1 ) { 
    sprintf(c1err,"NEVT_APPROX=%d exceeds MXSN=%d", NEVT_APPROX, MXSN);
    sprintf(c2err,"Probably need to increase MXSN");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    //     NEVT_APPROX = MXSN-1 ; 
  }

  // if 2nd file is much smaller, avoid too small malloc
  if ( ifile > 0 && NEVT_APPROX < NLIST_FIRST_FITRES ) 
    { NEVT_APPROX = NLIST_FIRST_FITRES + 2 ;  }


  // allocate new memory each time
  fitres_malloc_flt(ifile, NVARALL, NEVT_APPROX);
  fitres_malloc_str(ifile, NVARSTR, NEVT_APPROX);

  // ----------------------------------------------------

  ivarstr = 0 ;
  NVAR    = NVARALL_FITRES ;

  NVARALL_FITRES_LAST = NVARALL_FITRES ;
  NVARSTR_FITRES_LAST = NVARSTR_FITRES ;

  for ( ivar=0; ivar < NVARALL; ivar++ ) {

    // don't duplicate FIELD (Dec 8 2014)
    if ( SKIP_VARNAME(ifile,ivar) ) { continue ; }

    // get VARNAME and ICAST
    VARNAME = READTABLE_POINTERS.VARNAME[ivar] ;
    ICAST   = READTABLE_POINTERS.ICAST_STORE[ivar] ;

    
    if ( ifile==0 && strcmp(VARNAME,"zHD") == 0 ) { IVAR_zHD = ivar; }

    // Sep 19 2019: make sure first column is CID
    if ( ivar == IVARSTR_CCID ) {
      if ( ICAST_for_textVar(VARNAME) != ICAST_C ) {
	sprintf(c1err,"Unrecognized first column: %s", VARNAME);
	sprintf(c2err,"Check %s", INPUTS.FFILE[ifile] );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );       
      }
    }
 
    if ( ICAST == ICAST_C ) {
      NVARSTR_FITRES++ ;   // summed over fitres files
      sprintf(VARNAME_C, "%s:C*20", VARNAME) ;  //specify char ptr
      index = SNTABLE_READPREP_VARDEF(VARNAME_C, 
			     &FITRES_VALUES.STR_TMP[ivarstr][0], MXSN, 1 ) ;
      ivarstr++ ;          // sum for each separate fitres file
    }
    else {     
      sprintf(VARNAME_F, "%s:F", VARNAME) ;  // specify float pointer
      index = SNTABLE_READPREP_VARDEF(VARNAME_F, 
			     &FITRES_VALUES.FLT_TMP[ivar][0], MXSN, 1 ) ; 
    }


    if ( NVAR >= MXVAR_TOT ) {
      sprintf(c1err,"NVARALL_COMBINE=%d exceeds arround bound of", NVAR);
      sprintf(c2err,"MXVAR_TOT = %d ", MXVAR_TOT ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }      
    ptr_CTAG = VARNAME_COMBINE[NVAR] ;
    sprintf(ptr_CTAG, "%s", VARNAME );

    // use CID (1st var) in 1st file only
    REPEATCID = ( ifile > 0 && ivar == IVARSTR_CCID ) ;
    if ( REPEATCID ) { 
      ICAST_FITRES_COMBINE[NVAR] = -1 ;  // don't write repeated CIDs
      NVARALL_FITRES++ ;     NVAR = NVARALL_FITRES ;
      continue ;
    } 
    else { 
      ICAST_FITRES_COMBINE[NVAR] = ICAST ;
      NVAR_WRITE_COMBINED++ ; 
    } 
    
    // check if this tag name already exists;
    // if so, append '_$ifile' to end of tagname

    NTAG_DEJA = NMATCH_VARNAME(ptr_CTAG, NVAR) ;
    if ( NTAG_DEJA > 0 ) {
      iappend = ifile+1; // first iappend is 2
      printf("\t ADD_FITRES WARNING: VARNAME=%s already exists: ", ptr_CTAG);
      printf("append %d \n", iappend );
      sprintf( VARNAME_COMBINE[NVAR], "%s_%d", ptr_CTAG, iappend );
    }

    NVARALL_FITRES++ ;     NVAR = NVARALL_FITRES ;

    fflush(stdout);

  } // end of ivar loop


  // read everything and close file.
  NLIST = SNTABLE_READ_EXEC();
  NEVT_READ[ifile] = NLIST;

  if ( ifile == 0 ) {
    NLIST_FIRST_FITRES  = NLIST ; 
    malloc_NMATCH_PER_EVT(NLIST);
  }

  // always fill 2nd NLIST2
  NLIST2_FITRES = NLIST ;

  printf("\n"); fflush(stdout);

  // ===============================================

  if ( ifile > 0 ) 
    { printf("\t begin CID matching ... \n"); fflush(stdout); }

  NMATCH2 = 0 ;

  for(isn2=0; isn2 < NLIST2_FITRES; isn2++ ) {
    
    if ( INPUTS.MATCHFLAG == MATCHFLAG_HASH_LOCAL ) {
      isn = match_CID_hash_local(ifile,isn2);   // isn is for ifile=0
    }
    else if ( INPUTS.MATCHFLAG == MATCHFLAG_HASH_UTIL ) {
      sprintf(ccid, "%s", FITRES_VALUES.STR_TMP[IVARSTR_CCID][isn2]);
      isn = match_cid_hash(ccid,ifile,isn2);   // isn is for ifile=0
    }
    else 
      { isn = match_CID_orig(ifile,isn2);  }

    if ( isn >= 0 ) { 
      ADD_FITRES_VARLIST(ifile,isn,isn2);   
      NMATCH2++ ;
      NMATCH_PER_EVT[isn]++ ; 
    }

  }

  NEVT_MISSING[ifile] = NLIST_FIRST_FITRES - NMATCH2 ;

  fflush(stdout);

  /*
  if ( ifile==0 && IVAR_zHD < 0 ) {
    sprintf(c1err,"Cannot apply cut on zHD (%.2f to %.2f).",
	    INPUTS.CUTWIN_zHD[0], INPUTS.CUTWIN_zHD[1] );
    sprintf(c2err,"Could not find zHD column.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  */

  // free temp arrays
  freeVar_TMP(ifile, NVARALL, NVARSTR, NEVT_APPROX);

} // end of ADD_FITRES

// =====================================
int match_CID_orig(int ifile, int isn2) {

  // Brute force CID matching by searching both files.
  // isn2 is current SN index; functions returns isn of 
  // first fitres file such that cid(isn2) == cid(isn).
  //
  // The logic here optimizes the speed of finding a cid match
  // between two lists ... CPU becomes noticable when N_SN 
  // is order 10^5.

  char ccid[MXSTRLEN], ccid2[MXSTRLEN] ;
  int isn;
  //  char fnam[] = "match_CID_orig" ;

  // ----------- BEGIN ------------

  if ( ifile == 0 ) { return(isn2); }

  sprintf(ccid2, "%s", FITRES_VALUES.STR_TMP[IVARSTR_CCID][isn2] );   

  for ( isn = 0; isn < NLIST_FIRST_FITRES ; isn++ ) {

    // checking logical is faster than checking string
    if ( USEDCID[isn] ) { continue ; }

    sprintf(ccid,"%s", FITRES_VALUES.STR_ALL[IVARSTR_CCID][isn] );
    if ( strcmp(ccid,ccid2) == 0 ) {
      USEDCID[isn] = true ;
      return(isn);
    }
  }
  
  return(-9); 

  //  if ( ifile>0 )  { printf("\t (NLOOP_TOT = %lld)\n", NLOOP_TOT); }

} // end if match_CID_orig

// =====================================
int match_CID_hash_local(int ifile, int isn2) {

  // Created Oct 7 2019
  // Use hash table to speed matching.
  // Inputs:
  //   ifile = file index 
  //   isn2  = current SN index
  // 
  // Function returns isn index of ifile=0.

  int   isn;
  char  ccid[MXSTRLEN_CID];
  struct hash_table *s, *tmp;
  //  char fnam[] = "match_CID_hash_local" ;

  // ----------- BEGIN ------------

  if ( ifile < 0 ) {
    /* free the hash table contents */
    HASH_ITER(hh, users, s, tmp) {
      HASH_DEL(users, s);
      free(s);
    }
    return(-1);
  }


  if ( ifile == 0 ) {
    // create hash table
    sprintf(ccid, "%s", FITRES_VALUES.STR_TMP[IVARSTR_CCID][isn2]);
    s     = malloc(sizeof(struct hash_table));
    s->id = isn2;
    strcpy(s->name, ccid);
    HASH_ADD_STR( users, name, s );      
    return(isn2) ;
  }

  // if we get here, match each CID to ifile=0
  
  sprintf(ccid, "%s", FITRES_VALUES.STR_TMP[IVARSTR_CCID][isn2] );
  HASH_FIND_STR( users, ccid, s);
  if ( s ) {
    isn = s->id;
    //    printf("\t xxx CID=%8s matches to %d in ifile=0\n",  ccid, isn);
    return(isn);	     
  }
  else {
    //    printf("\t xxx CID=%s does not match.\n", ccid);
    return(-9);
  }
  

  return(-9);

} // end if match_CID_hash_local


// =============================
void ADD_FITRES_VARLIST(int ifile, int isn, int isn2) {

  int  NVARALL      = NVARALL_FILE[ifile];
  int  LTMP = 0 ;
  int  MXUPDATE = 50;
  int  ivarstr, IVARSTR, IVARTOT, ivar, ICAST, TMPMOD ;
  char ccid2[MXSTRLEN_CID];
  char fnam[] = "ADD_FITRES_VARLIST" ;

  // --------- BEGIN ------------

  // check for screen-update

  TMPMOD = (isn2 % 10000);
  if ( TMPMOD == 0.0               )    { LTMP = 1; }
  if ( isn2 <= MXUPDATE            )    { LTMP = 1; }
  if ( isn2 == NLIST_FIRST_FITRES-1)    { LTMP = 1; } 
  if ( isn2 == NLIST_FIRST_FITRES  )    { LTMP = 1; }

  if ( ifile <= 1 && LTMP  ) {
    sprintf(ccid2, "%s", FITRES_VALUES.STR_TMP[IVARSTR_CCID][isn2] ); 
    printf("\t %s = '%12s'  ->  isn = %6d   \n", 
	   VARNAME_COMBINE[0], ccid2, isn2 );  fflush(stdout);
  }
    
  if ( isn2 == MXUPDATE+1 ) { printf("\t ... \n"); }

  ivarstr = 0;
  IVARTOT = NVARALL_FITRES_LAST ;

  for ( ivar=0; ivar < NVARALL; ivar++ ) {

    if ( SKIP_VARNAME(ifile, ivar) ) { continue ; }

    // VARNAME = READTABLE_POINTERS.VARNAME[ivar] ;
    ICAST   = READTABLE_POINTERS.ICAST_STORE[ivar] ;

    if ( ICAST != ICAST_C )  {   // not a string
      
      if (  isnan(FITRES_VALUES.FLT_TMP[ivar][isn2]) !=0 ) {
	sprintf(c1err,"isnan for FLT_TMP[ivar=%d][isn=%d] = %f", 
		ivar, isn2, FITRES_VALUES.FLT_TMP[ivar][isn2]  );
	sprintf(c2err,"varname = %s", VARNAME_COMBINE[ivar] );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      
      FITRES_VALUES.FLT_ALL[IVARTOT][isn] = 
	FITRES_VALUES.FLT_TMP[ivar][isn2] ; 
    }
    else {
      IVARSTR = NVARSTR_FITRES_LAST + ivarstr ;
      IVARSTR_STORE[IVARTOT] = IVARSTR ;
      
      sprintf(FITRES_VALUES.STR_ALL[IVARSTR][isn],"%s", 
	      FITRES_VALUES.STR_TMP[ivarstr][isn2] ); 
      
      ivarstr++ ;
    }
    
    IVARTOT++ ;
    
  } // ivar
  
} // end ADD_FITRES_VARLIST

 
// =====================================
int SKIP_VARNAME(int ifile, int ivar) {

  // Dec 8 2014
  // Return 1 if this variable should be ignored.
  // Dec 2020: check KEEP_VARNAMES 

  char *VARNAME = READTABLE_POINTERS.VARNAME[ivar] ;
  bool KEEP = false;
  int  j, k ;

  // -------- BEGIN ------------

  // check option to keep only use-selected list of varnames
  for ( k=0; k < INPUTS.NVARNAMES_KEEP; k++ ) {
    if ( strcmp(VARNAME,INPUTS.VARNAMES_KEEP[k]) == 0 ) { KEEP=true;}
  }
  if ( ivar>0 && INPUTS.NVARNAMES_KEEP > 0 && !KEEP ) { return(1); }

  // - - - - - -
  // no SNANA file yet
  if ( IFILE_FIRST_SNANA < 0 ) { return(0); } 

  // never skip 1st SNANA file.
  if ( ifile == IFILE_FIRST_SNANA ) { return(0) ; }  

  for(j=0; j < NVARNAME_1ONLY; j++ ) {
    if ( strcmp(VARNAME,VARNAME_1ONLY[j]) == 0 ) { return(1) ; }
  }

  return(0);

} // end of SKIP_VARNAME


// ====================================
void freeVar_TMP(int ifile, int NVARTOT, int NVARSTR, int MAXLEN) {

  // Aug 2013
  // Free _TMP arrays so that they can be re-allocated
  // with a different number of variables and SN.

  int ivar, isn ;

  for ( ivar=0; ivar < NVARTOT; ivar++ ) {

    if ( ivar < NVARSTR ) {

      for ( isn=0; isn < MAXLEN; isn++ ) 
	{ free(FITRES_VALUES.STR_TMP[ivar][isn]) ; }

      free( FITRES_VALUES.STR_TMP[ivar]  ) ;
    }

    free( FITRES_VALUES.FLT_TMP[ivar] ) ;

  }  // ivar

  free(FITRES_VALUES.FLT_TMP);
  free(FITRES_VALUES.STR_TMP);

} // end of freeVar_TMP

// =====================================
void  fitres_malloc_flt(int ifile, int NVAR, int MAXLEN) {

  // allocate memory for FITRES_VALUES.TMP to read current 
  // fitres file, and extend memory for FITRES_VALUES.FLT_ALL 
  // to store all fitres values.
  // NVAR is the number of variables to read from this fitres file.
  // MAXLEN is an estimate of the max array length to allocate.

  int ivar, isn, IVAR_ALL, NTOT, MEMF ;
  //  char fnam[] = "fitres_malloc_flt" ;

  // ---------- BEGIN ------------

  printf("\t Allocate fitres memory for NVAR=%d and MAXLEN=%d (ifile=%d)\n",
	 NVAR, MAXLEN, ifile );

  // redo malloc on TMP arrays for each fitres file
  MEMF      = (NVAR+1) * sizeof(float*) ;
  FITRES_VALUES.FLT_TMP = (float **)malloc(MEMF) ;

  // -----------------------------------
  NTOT = NVARALL_FITRES + NVAR + 1 ;
  MEMF = sizeof(float*) * NTOT ;

  if (ifile == 0 ) 
    { FITRES_VALUES.FLT_ALL = (float**)malloc(MEMF); }
  else  
    { FITRES_VALUES.FLT_ALL = (float**)realloc(FITRES_VALUES.FLT_ALL,MEMF); }

  
  for ( ivar=0; ivar < NVAR; ivar++ ) {

    MEMF = sizeof(float  ) * MAXLEN ;
    IVAR_ALL = NVARALL_FITRES + ivar ;

    FITRES_VALUES.FLT_TMP[ivar]     = (float  *)malloc(MEMF);
    FITRES_VALUES.FLT_ALL[IVAR_ALL] = (float  *)malloc(MEMF);    

    for ( isn=0; isn < MAXLEN; isn++ ) {
      USEDCID[isn] = false ;
      FITRES_VALUES.FLT_TMP[ivar][isn]     = INPUTS.NULLVAL_FLOAT ;
      FITRES_VALUES.FLT_ALL[IVAR_ALL][isn] = INPUTS.NULLVAL_FLOAT ;
    }
  }


} // end of fitres_malloc_flt




// =====================================
void  fitres_malloc_str(int ifile, int NVAR, int MAXLEN) {

  // allocate memory for string-variables
  // NVAR is the number of string variables in this fitres file.
  // Note that NVAR >= 1 because the CID string must always
  // be there.
  //
  // Apr 27 2020: init STR_ALL and STR_TMP to 'NULL'

  //  char fnam[] = "fitres_malloc_str" ;
  int ivar, IVAR_ALL, isn, MEMC, NTOT ;
  
  // ---------- BEGIN ------------

  printf("\t Allocate string memory for NVAR=%d and MAXLEN=%d \n",
	 NVAR, MAXLEN );

  // re-do malloc on TMP arrays for each fitres file

  MEMC                    = (NVAR+1) * sizeof(char**) ;
  FITRES_VALUES.STR_TMP   = (char***)malloc(MEMC) ; 

  // -----------------------------------
  NTOT = NVARSTR_FITRES + NVAR + 1 ;
  MEMC = NTOT * sizeof(char**);

  if (ifile == 0 ) 
    { FITRES_VALUES.STR_ALL = (char***)malloc(MEMC) ; }
  else 
    { FITRES_VALUES.STR_ALL = (char***)realloc(FITRES_VALUES.STR_ALL,MEMC); }

  
  for ( ivar=0; ivar < NVAR; ivar++ ) {

    IVAR_ALL = NVARSTR_FITRES + ivar ;

    // allocate SN-dimension
    MEMC = sizeof(char*) * MAXLEN ;
    FITRES_VALUES.STR_TMP[ivar]      = (char**)malloc(MEMC);
    FITRES_VALUES.STR_ALL[IVAR_ALL]  = (char**)malloc(MEMC);    

    // allocate chars for each string argument
    MEMC = MXSTRLEN * sizeof(char)  ;
    for ( isn=0; isn < MAXLEN; isn++ ) {
      FITRES_VALUES.STR_TMP[ivar][isn]     = (char*)malloc(MEMC);
      FITRES_VALUES.STR_ALL[IVAR_ALL][isn] = (char*)malloc(MEMC);
      sprintf(FITRES_VALUES.STR_TMP[ivar][isn],"%s"    , DEFAULT_NULLVAL_STRING) ;
      sprintf(FITRES_VALUES.STR_ALL[IVAR_ALL][isn],"%s", DEFAULT_NULLVAL_STRING) ;
   
    } // isn
  }  // ivar


} // end of fitres_malloc_str

// ===================================
void malloc_NMATCH_PER_EVT(int N) {
  int i;
  NMATCH_PER_EVT = (char*) malloc( N * sizeof(char) );
  for(i=0; i < N; i++ ) { NMATCH_PER_EVT[i] = 0 ; }
} 

// ===================================
int NMATCH_VARNAME(char *ctag , int ntlist ) {

  // Returns number of *ctag matches in first "ntlist"
  // elements of CTAG array
  // Jun 7 2020: match with ignoreCase

  int i, NMATCH, OVP ;

  NMATCH = 0;

  for ( i=0; i < ntlist; i++ ) {
    OVP    = strcmp_ignoreCase(ctag,VARNAME_COMBINE[i]) ;
    if ( OVP == 0 )  { NMATCH++ ; }
  }

  return NMATCH;

} // end of NMATCH_VARNAME


// =========================================
void  print_stats(void) {

  // Created Jan 4 2021
  // print summary as follows

  //   ifile  NEVT_READ  NEVT_MISSING 
  //
  // and then
  //   NEVT_COMMON: 
  //

  int NFFILE = INPUTS.NFFILE ;
  int ifile, isn ;
  char key_grep[] = "SUMMARY:" ;

  // ------------ BEGIN ------------

  sprintf(BANNER,"FILE_STATISTICS %s", key_grep );
  print_banner(BANNER);

  printf("%s \n", key_grep);
  printf("%s   ifile   NEVT_READ   NEVT_MISSING \n", key_grep);
  printf("%s  ------------------------------------\n", key_grep);
  for(ifile=0; ifile < NFFILE; ifile++ ) {
    printf("%s   %2d    %8d   %8d \n",
	   key_grep, ifile,  NEVT_READ[ifile], NEVT_MISSING[ifile] );
    fflush(stdout);
  }
  printf("%s  ------------------------------------\n", key_grep);


  if ( NFFILE > 1 ) {
    for(isn=0; isn < NLIST_FIRST_FITRES; isn++ ) {
      if ( NMATCH_PER_EVT[isn] == NFFILE ) { NEVT_COMMON++; }
    }
    printf("%s NEVT_COMMON: %d  (%d missing in at least one file)\n\n", 
	   key_grep, NEVT_COMMON, NEVT_READ[0]-NEVT_COMMON );
    fflush(stdout);
  }

  return;

} // end print_stats

// =========================================
void WRITE_SNTABLE(void) {

  // Feb 2, 2013
  // create table.
  // by pre-processor flag.
  // This routine replaces NTUP_SHELL and WR_CWNT().
  //
  // Jan  7, 2014: for hbook, call remove_string_termination(...)
  // Apr 26, 2014: return gracefully if neither hbook & root are defined.
  // Feb 26, 2017: if CIDint already exists, don't write out another one.

  char 
    OUTFILE[6][200] 
    ,tableVar[60]
    ,*ptrSTR
    ,BLOCKVAR[]   = "VAR"     
    ,openOpt[40], CCIDint[40]
    ;

  double zHD;
  int GZIPFLAG = 0 ;
  int ivar, ivarstr, isn, ICAST, CIDint ;
  int IFILETYPE, NOUT, out, SKIP ;

  char  fnam[] = "WRITE_SNTABLE" ;
  // --------------- BEGIN ------------

  NOUT = 0 ;
  OUTFILE[NOUT][0] = 0 ;
  NWRITE_SNTABLE = 0 ;


#ifdef USE_HBOOK
  if (CREATEFILE_HBOOK)  { 
    sprintf(OUTFILE[NOUT], "%s.%s", 
	    INPUTS.OUTPREFIX_COMBINE, ptrSuffix_hbook );  
    sprintf(openOpt,"%s new", ptrSuffix_hbook);
    IFILETYPE = TABLEFILE_OPEN(OUTFILE[NOUT],openOpt);
    NOUT++ ;
  }
#endif


#ifdef USE_ROOT
  if ( CREATEFILE_ROOT )  { 
    sprintf(OUTFILE[NOUT], "%s.%s", 
	    INPUTS.OUTPREFIX_COMBINE, ptrSuffix_root ); 
    sprintf(openOpt,"%s new", ptrSuffix_root);
    IFILETYPE = TABLEFILE_OPEN(OUTFILE[NOUT],openOpt);
    NOUT++ ;
  }
#endif


#ifdef USE_TEXT
  if ( CREATEFILE_TEXT )  { 
    sprintf(OUTFILE[NOUT], "%s.%s", 
	    INPUTS.OUTPREFIX_COMBINE, ptrSuffix_text ); 
    outFile_text_override(OUTFILE[NOUT],&GZIPFLAG); 
    sprintf(openOpt,"%s new", ptrSuffix_text);
    //    IFILETYPE = TABLEFILE_OPEN(OUTPREFIX_COMBINE,openOpt);
    IFILETYPE = TABLEFILE_OPEN(OUTFILE[NOUT],openOpt);
    NOUT++ ;
  }
#endif

  printf("\n   Create combined SNTable with %d variables \n", 
	 NVAR_WRITE_COMBINED );

  SNTABLE_CREATE(TABLEID_COMBINE, TABLENAME_COMBINE, "KEY" ) ;

  ADD_SNTABLE_COMMENTS() ;

  // check of CIDint is already there (Feb 2017)
  int CIDint_EXISTS = 0 ;
  for ( ivar=0; ivar < NVARALL_FITRES; ivar++ ) {
    if ( strcmp(VARNAME_COMBINE[ivar],"CIDint")==0 ) { CIDint_EXISTS=1; }
  }


  for ( ivar=0; ivar < NVARALL_FITRES; ivar++ ) {

    ICAST = ICAST_FITRES_COMBINE[ivar] ;

    if ( ICAST < 0 ) {  continue ; }

    if ( ICAST == ICAST_C ) {
      int MXLEN ;
      MXLEN = maxlen_varString(VARNAME_COMBINE[ivar]);

      // include char CCID & integer CID representation for CID
      if ( ivar == 0 ) {

	sprintf(tableVar, "CCID:C*%d ", MXLEN );  
	SNTABLE_ADDCOL(TABLEID_COMBINE, BLOCKVAR, 
		       TABLEROW_VALUES.CCID, tableVar, 1 ) ;

	if ( CIDint_EXISTS == 0 ) {
	  sprintf(tableVar, "CIDint:I" );  
	  SNTABLE_ADDCOL(TABLEID_COMBINE, BLOCKVAR, 
			 &TABLEROW_VALUES.CIDint, tableVar, 1 ) ;
	}
      }
      else {       
	ivarstr = IVARSTR_STORE[ivar];
	sprintf(tableVar,"%s:C*%d", VARNAME_COMBINE[ivar], MXLEN );  
	SNTABLE_ADDCOL(TABLEID_COMBINE, BLOCKVAR, 
	TABLEROW_VALUES.STR[ivarstr], tableVar, 1 ) ;
      }

    }
    else {
      sprintf(tableVar,"%s:F", VARNAME_COMBINE[ivar] );  
      SNTABLE_ADDCOL(TABLEID_COMBINE, BLOCKVAR, 
		     &TABLEROW_VALUES.FLT[ivar], tableVar, 1) ;
    }

  }    // ivar


  // ------------------------------
  printf("   Fill combined table with %d rows ... \n", NLIST_FIRST_FITRES );
  fflush(stdout);

  if ( INPUTS.DOzCUT ) {
    printf("\n ONLY WRITE EVENTS with %.3f < zHD < %.3f \n\n",
	   INPUTS.CUTWIN_zHD[0], INPUTS.CUTWIN_zHD[1] );
    fflush(stdout);
  }

  for ( isn = 0; isn < NLIST_FIRST_FITRES ; isn++ ) {

    TABLEROW_VALUES.CIDint = -999;
    sprintf(TABLEROW_VALUES.CCID,"xxx ");

    for ( ivar=0; ivar < NVARALL_FITRES; ivar++ ) {

      ICAST = ICAST_FITRES_COMBINE[ivar] ;
      if ( ICAST < 0 ) { continue ; }

      if ( ICAST == ICAST_C )  { 
	ivarstr = IVARSTR_STORE[ivar];
	ptrSTR  = TABLEROW_VALUES.STR[ivarstr] ;

	sprintf(ptrSTR,"%s", FITRES_VALUES.STR_ALL[ivarstr][isn]);

	if ( ivar == 0 ) {
	  sprintf(TABLEROW_VALUES.CCID, "%s ", ptrSTR); // Dec 8 2014
	  CIDint = atoi(ptrSTR); 
	  sprintf(CCIDint,"%d", CIDint);
	  if ( strcmp(CCIDint,ptrSTR) == 0 ) 
	    { TABLEROW_VALUES.CIDint = CIDint; } // CCID is already int
	  else { 
	    // CCID is a string, so integer CIDint increments
	    TABLEROW_VALUES.CIDint = isn ; 
	  }
	}

	if ( CREATEFILE_HBOOK ) 
	  { remove_string_termination(ptrSTR, MXSTRLEN); }

      }
      else
	{ TABLEROW_VALUES.FLT[ivar] = FITRES_VALUES.FLT_ALL[ivar][isn];  }
      
      /*
      float VAL = TABLEROW_VALUES.FLT[ivar];
      if ( ivar == 54 && VAL < 0.0 ) {
	printf(" xxx %s: CID=%d  %s = %le \n", 
	       fnam, CIDint,  VARNAME_COMBINE[ivar], 
	       TABLEROW_VALUES.FLT[ivar] );
	fflush(stdout);
	}  */


    } // ivar

    SKIP = 0 ;
    if ( INPUTS.DOzCUT ) {
      zHD = TABLEROW_VALUES.FLT[IVAR_zHD];
      if ( zHD < INPUTS.CUTWIN_zHD[0] ) { SKIP = 1; }
      if ( zHD > INPUTS.CUTWIN_zHD[1] ) { SKIP = 1; }
    }
    if ( SKIP ) { continue ; }

    NWRITE_SNTABLE++ ;
    SNTABLE_FILL(TABLEID_COMBINE);

    // Jan 2020: stop if -mxrow
    if ( isn >= INPUTS.MXROW_READ-1 ) {
      printf("\n\t STOP AFTER WRITING %d ROWS. \n\n", isn);
      fflush(stdout);  goto DONE_FILL ;
    }
    
  } // isn

 DONE_FILL:

  // close it
  for(out=0; out < NOUT; out++ ) 
    { TABLEFILE_CLOSE(OUTFILE[out]);  }

  // check gzip option
  if ( GZIPFLAG )  { 
    char cmd[400];
    sprintf(cmd,"gzip %s", INPUTS.OUTFILE_TEXT);
    system(cmd); 
  }

  return;

} // end of WRITE_SNTABLE


// ====================================================
void outFile_text_override(char *outFile, int *GZIPFLAG) {

  // June 2019
  // if user-input OUTFILE_TEXT is set, then set 
  //  outFile = OUTFILE_TEXT.
  //
  // If there is a .gz extension on OUTFILE_TEXT,
  // remove .gz extension and set *GZIPFLAG=1.

  char *OUTFILE_TEXT = INPUTS.OUTFILE_TEXT ;
  int LEN = strlen(OUTFILE_TEXT) ;

  // -------------- BEGIN -------------
  if ( LEN <2 ) { return; }

  // check for .gz extension

  if ( strstr(OUTFILE_TEXT,".gz") != NULL )
    { OUTFILE_TEXT[LEN-3] = 0;  *GZIPFLAG=1; }

  sprintf(outFile,"%s", OUTFILE_TEXT);


  return ;

} // end outFile_text_override

// ==================================
void  ADD_SNTABLE_COMMENTS(void) {

  // Call STORE_TABLEFILE_COMMENT(comment)
  // May 2 2019: [120] -> [MXPATHLEN]

  int ifile ;
  char comment[MXPATHLEN];

  // ------------ BEGIN ------------
  sprintf(comment,"Created by combine_fitres.exe");
  STORE_TABLEFILE_COMMENT(comment) ;

  sprintf(comment,"Number of combined files: %d", INPUTS.NFFILE);
  STORE_TABLEFILE_COMMENT(comment) ;

  for(ifile=0; ifile < INPUTS.NFFILE; ifile++ ) {
    sprintf(comment,"   + %s", INPUTS.FFILE[ifile] );
    STORE_TABLEFILE_COMMENT(comment) ;
  }

  return;

} // end of WRITE_TABLE_COMMENTS

// ======================================
int  maxlen_varString(char *varName) {

  // Jan 2014
  // return max string length for variable nameed *varName.
  // Dec 4 2016: add VERSION

  int MXLEN ;
  MXLEN = MXSTRLEN ;

  if ( strcmp(varName,"CID") == 0 ) 
    { MXLEN = MXSTRLEN_CID; }

  if ( strcmp(varName,"CCID") == 0 ) 
    { MXLEN = MXSTRLEN_CID; }

  if ( strcmp(varName,"STARID") == 0 ) 
    { MXLEN = MXSTRLEN_CID; }

  if ( strcmp(varName,"ROW") == 0 ) 
    { MXLEN = MXSTRLEN_CID; }

  if ( strcmp(varName,"BAND") == 0 ) 
    { MXLEN = MXSTRLEN_BAND; }
  
  if ( strcmp(varName,"FIELD") == 0 ) 
    { MXLEN = MXSTRLEN_FIELD; }

  if ( strcmp(varName,"IAUC_NAME") == 0 ) 
    { MXLEN = MXSTRLEN_CID; }

  if ( strcmp(varName,"VERSION") == 0 ) 
    { MXLEN = MXSTRLEN_VERSION; }


  return MXLEN ;

} // end of void  maxlen_varString

