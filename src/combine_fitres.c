/*********************************
  Jun 23, 2009 by R.Kessler
 
  ------------
  Re-written in C to use utilities that allow
  arbitrary-length fitres file. Current fortran
  version is memory limited because it reads the
  entire file as char strings.
  
  Combine multiple fitres contents into a
  single ntuple. Designed to study correlations
  between different fitres outputs for the same
  set of SN.


 Usage:
  >  combine_fitres.exe <fitres1> <fitres2> ...

   where each command-line argument is the name of a
   fitres file to combine into the ntuple.
   Output files are
      combine_fitres.hbook  (hbook column-wise ntuple)
      combine_fitres.root   (root tree)
      combine_fitres.txt    (combined fitres text-file)

 Option:
  >  combine_fitres.exe <fitres1> <fitres2> --outprefix <outprefix>
    produces output files <outprefix>.hbook and  <outprefix>.txt

  >  combine_fitres.exe <fitres1>  R      ! .ROOT extension 
  >  combine_fitres.exe <fitres1>  r      ! .root extension
       (create root file instead of hbook file)

  >  combine_fitres.exe <fitres1>  H      ! .HBOOK extension
  >  combine_fitres.exe <fitres1>  h      ! .hbook extension
       (create default hbook file ; same as with no H or h arg)

  >  combine_fitres.exe <fitres1>  R H 
  >  combine_fitres.exe <fitres1>  H R
      (create both hbook and root files)

  >  combine_fitres.exe <fitres1>  -mxrow 50000


 WARNINGS/NOTES:
 * If fitres files contain different SN, then first
   fitres file determines the list of SN; extra SN
   in the 2nd, 3rd ... fitres file are ignored.


                HISTORY

 Feb 02, 2013
   - remove legacy Row-wise ntuple option

   - use new make_SNTABLE function that uses sntools_output.c;
     NTUP_SHELL and WR_CWNTUP are commented out and marked
     for deletion.

 Fev 24, 2013: END_HFILE  -> CLOSE_HFILE and add ROOTFILE calls.

 Apr 14, 2013: fix bug setting ccid in BREAK_CCIDSTRING(int ifile) 

 Apr 23, 2013: for repeat variable append '_2' instead of '2'.
               Idem for '_3' instead of '3', etc ...

 May 1, 2013: 
    - fix bug calling OPEN_ROOTFILE (was missing IERR arg)
    - new optional args 'H' and 'R' to turn on HBOOK, ROOT optoin.
      No table format -> default hbook format.

 Aug 28 2013: major overhaul to handle strings: CID, FIELD, etc ...

 Dec 21 2013: fix bug setting LEN_APPROX when 2nd file is much
              smaller than 1st file.

 Jan 6 2014: in WRITE_SNTABLE, pad string variables (CID,GALID ...)
	     with blanks instead of '\0' to pass as fortran arg.
             Allows paw-logic such as cid='[cid]' that previously
             did not work.
        
 Mar 3, 2014: 
    - MXSTRLEN_CID-> 20 (was 12)
    - for ntuple/tree, add CIDint = atoi(CCID) to have integer represenation
       --> allows making CID cuts and plots.

 Apr 26 2014:
   replace file-type specific functions (HFILE_OPEN, ROOTFILE_OPEN)
   with generic wrapper functions, TABLEFILE_xxx.
   Only call wrapper functions in sntools_output.c, 
   and NOT from sntools_output_[type].c

   Change .his -> .hbook (same in split_and_fit.pl)

 Sep 13 2014:  H,R -> upper case extension .HBOOK/.ROOT.

 Sep 14, 2014: replace access to SNTOOLS_FITRES struct with
               call to SNTABLE_GET_VARINFO(); because all the
               fitres-read code is now in sntools_output_text.c 

  Oct 17 2014: add extra call to SNTABLE_GET_VARINFO to fix bad bug.

  Oct 27, 2014: used refactored file-reading. No more WRITE_FITRES
                since WRITE_SNTABLE writes all table formats in one shot.

  Dec 8 2014: 
     in ADD_FITRES, skip duplicate FIELD. Logic is a bit clumsy,
     so may need to be refactored. See new function SKIP_VARNAME().

     On CID, add CCID in addition to existing CIDint.
     --> CCID must exist in the hbook/table file.

  Feb 15 2015: MXSN -> 2 million (was 1 million)
 
  Mar 31 2016: MXFILE -> 20 ( was 10)

  Dec 04 2016: allow VERSION column

  Feb 27 2017: add check for CCID being string or integer;
               CIDint -> isn if CCID is a string.

  Mar 8 2017: 
    make VARNAME_1ONLY list of SNANA variables to include only once. 
    Also check for first SNANA file (IFILE_FIRST_SNANA) in case 
    first file was not made by SNANA.

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

// ================================
// Function Declarations
// ================================


void  PARSE_ARGV(int argc, char **argv);
void  INIT_TABLEVAR(void);
void  ADD_FITRES(int ifile);
int   NMATCH_VARNAME(char *ctag , int ntlist ) ;
int   maxlen_varString(char *varName);
int   SKIP_VARNAME(int file, int ivar) ;

void  WRITE_SNTABLE(void); // output table in selected format
void  ADD_SNTABLE_COMMENTS(void) ;

void  fitres_malloc_flt(int ifile, int NVAR, int MAXLEN); 
void  fitres_malloc_str(int ifile, int NVAR, int MAXLEN); 
void  freeVar_TMP(int ifile, int NVARTOT, int NVARSTR, int MAXLEN); 

// ================================
// Global variables
// ================================


#define MXFFILE  20       // max number of fitres files to combine
#define MXSN     2000000   // max SN to read per fitres file
#define MXVAR_PERFILE  50  // max number of NTUP variables per file
#define MXVAR_TOT  MXVAR_TABLE     // max number of combined NTUP variables
#define INIVAL_COMBINE  -888.0
#define MXSTRLEN         28
#define MXSTRLEN_BAND      4
#define MXSTRLEN_CID      20
#define MXSTRLEN_FIELD    20
#define MXSTRLEN_VERSION  32 // photometry versoin, added Dec 2016

#define IVARSTR_CCID  1   // CCID index for CVAR_XXX arrays


int USEFILE_HBOOK, USEFILE_ROOT, USEFILE_TEXT ;

// cast for each variable; negative -> do not write
// See ICAST_FITRES[D,F,I,C] parameters in sntools.h
int ICAST_FITRES_COMBINE[MXVAR_TOT]; 


// for SNtable
#define  TABLEID_COMBINE   TABLEID_FITRES
#define  TABLENAME_COMBINE TABLENAME_FITRES

// files
int  NFFILE_INPUT ;   // number of input fitres files to combine
char FFILE_INPUT[MXFFILE][2*MXPATHLEN] ;
char OUTPREFIX_COMBINE[MXPATHLEN] = "combine_fitres" ;
int  MXROW_READ ;

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

short int USEDCID[MXSN];

int IVARSTR_STORE[MXVAR_TOT] ; // keep track of string vars

int NLIST_FIRST_FITRES ;   // number of SN in 1st fitres file
int NLIST2_FITRES ;        // idem for 2nd, 3rd, etc ...

//  table variables
int   NVAR_WRITE_COMBINED ;  // total number of variables in combined file
int   NVARALL_FITRES ;       // idem, but including multiple CIDs
int   NVARSTR_FITRES ;

struct FITRES_VALUES {

  char  ***STR_ALL ;  // string values for all files [ivar][isn]
  char  ***STR_TMP ;  // string values for current file

  float  **FLT_ALL ;  // values for all files  [ivar][isn]
  float  **FLT_TMP ;  // idem for 1 file

} FITRES_VALUES ;

char  VARNAME_COMBINE[MXVAR_TOT][MXSTRLEN]; // name of each combined variable


struct ROW_VALUES {
  float FLT[MXVAR_TOT] ;  
  char  STR[MXVAR_TOT][MXSTRLEN] ;
  char  CCID[MXSTRLEN_CID];     // added Dec 8 2014
  int   CIDint ;
} TABLEROW_VALUES ;       // values for 1 row; used for SNTABLE


// =========================================
int main(int argc, char **argv) {

  int i, ifile ;
  char fnam[] = "main" ;

  // ----------------- BEGIN --------

  sprintf(BANNER,"Begin execution of combine_fitres.exe \n" );
  print_banner(BANNER);

  set_EXIT_ERRCODE(EXIT_ERRCODE_combine_fitres);

  // set default output to hbook if both root and hbook are compiled
  // (specified by order of inits below)
#ifdef USE_ROOT  
  USEFILE_ROOT = 1 ;  USEFILE_HBOOK = 0 ;
#endif

#ifdef USE_HBOOK
  USEFILE_ROOT = 0 ;  USEFILE_HBOOK = 1 ;
#endif

#ifdef USE_TEXT
  USEFILE_TEXT = 1;
#endif

  PARSE_ARGV(argc,argv);

  INIT_TABLEVAR();

  print_banner("Begin Reading Fitres Files.\n");

  TABLEFILE_INIT(); // Oct 27 2014

  for ( ifile = 1; ifile <= NFFILE_INPUT; ifile++ ) {
    ADD_FITRES(ifile);
  }

  // ---------------

  WRITE_SNTABLE() ;

  printf("   Done. \n");
  fflush(stdout);

} // end of main



// ===============================
void  PARSE_ARGV(int argc, char **argv) {

  int  i;
  char ctmp[80];
  char fnam[] = "PARSE_ARGV" ;

  //----------- BEGIN --------------

  NARGV_LIST = argc ; // set global 

  NFFILE_INPUT = 0;
  MXROW_READ   = 1000000000 ;

  for ( i = 1; i < NARGV_LIST ; i++ ) {
    
    // check for optional args

    if ( strcmp(argv[i],"--outprefix") == 0 ) {
      i++ ; sprintf(OUTPREFIX_COMBINE,"%s", argv[i]);
      continue ;
    }
    if ( strcmp(argv[i],"-outprefix") == 0 ) {
      i++ ; sprintf(OUTPREFIX_COMBINE,"%s", argv[i]);
      continue ;
    }

    if ( strcmp(argv[i],"-mxrow") == 0 ) {
      i++ ; sscanf(argv[i], "%d", &MXROW_READ);
      continue ;
    }

    if ( strcmp_ignoreCase(argv[i],"r") == 0 ) { 
      USEFILE_ROOT = 1;  
      if ( USEFILE_HBOOK == 1 ) { USEFILE_HBOOK = 0 ; }  // root on, hbook off
      if ( strcmp(argv[i],"R")==0 ) { ptrSuffix_root = SUFFIX_ROOT ; }
      continue ;
    }

    if ( strcmp_ignoreCase(argv[i],"h") == 0 ) { 
      USEFILE_HBOOK = 2 ;  // turn hbook back on and leave it on
      if ( strcmp(argv[i],"H")==0) { ptrSuffix_hbook = SUFFIX_HBOOK ; }
      continue ;
    }

    NFFILE_INPUT++ ;
    sprintf( FFILE_INPUT[NFFILE_INPUT], "%s", argv[i] );
    printf("  Will combine fitres file: %s \n", FFILE_INPUT[NFFILE_INPUT] );
  }

  if ( NFFILE_INPUT <= 0 ) {
    sprintf(c1err, "Bad args. Must give fitres file(s)");
    sprintf(c2err, "  combine_fitres.exe <fitresFile List> ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // idiot checks for file type

#ifndef USE_ROOT
  if ( USEFILE_ROOT ) {
    sprintf(c1err, "Cannot create output ROOT file because");
    sprintf(c2err, "'#define ROOT' is not set in sntools_output.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
#endif


#ifndef USE_HBOOK
  if ( USEFILE_HBOOK ) {
    sprintf(c1err, "Cannot create output HBOOK file because");
    sprintf(c2err, "'#define HBOOK' is not set in sntools_output.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
#endif


} // end of PARSE_ARGV



// ==========================
void INIT_TABLEVAR(void) {

  int isn, ivar;

  NVAR_WRITE_COMBINED = 0 ;

  NVARALL_FITRES = 0 ;
  NVARSTR_FITRES = 0 ;

  NLIST_FIRST_FITRES = NLIST2_FITRES = 0 ;

  IFILE_FIRST_SNANA = -9;

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
  // Dec 1 2017: use new function SNTABLE_NEVT_APPROX_TEXT to check for .gz file.
  //
  int 
    ivar, IVARTOT, IVARSTR, ivartot, ivarstr, j
    ,isn, isn2,  ISN, ICID, ICAST, LTMP, IGNORE
    ,NVARALL_FILE, NVARSTR_FILE, NVAR
    ,NTAG_DEJA, NLIST
    ,MXUPDATE = 50
    ,index, REPEATCID
    ,NVARALL_LAST
    ,NVARSTR_LAST
    ,istat, isize, NEVT_APPROX, IFILETYPE
    ;

  struct stat statbuf;

  char 
    *VARNAME, VARNAME_F[MXCHAR_VARNAME], VARNAME_C[MXCHAR_VARNAME]
    ,*ptr_CTAG
    ,ccid[MXSTRLEN], ccid2[MXSTRLEN]
    ,fnam[] = "ADD_FITRES"
    ;

  // ----------- BEGIN -----------


  printf("\n --------------------- %s -------------------------\n", fnam );

  NVARALL_LAST = NVARALL_FITRES ;
  NVARSTR_LAST = NVARSTR_FITRES ;

  // open file & read header; use generic table name SNTABLE
  // since for ascii files the table name is not used.
  IFILETYPE = TABLEFILE_OPEN( FFILE_INPUT[ifile], "read text" );
  NVARALL_FILE = SNTABLE_READPREP(IFILETYPE,"SNTABLE");

  // check if this is an SNANA file; mark first SNANA file
  if ( IFILE_FIRST_SNANA < 0 ) {
    for ( ivar=1; ivar <= NVARALL_FILE; ivar++ ) {    
      VARNAME = READTABLE_POINTERS.VARNAME[ivar-1] ;
      for(j=0; j < NVARNAME_1ONLY; j++ ) {
	if ( strcmp(VARNAME,VARNAME_1ONLY[j])==0 ) 
	  {  IFILE_FIRST_SNANA = ifile; }	
      }
    }
    if ( IFILE_FIRST_SNANA >= 0 ) {
      printf("\t First SNANA ifile = %d (%s) \n",
	     IFILE_FIRST_SNANA, FFILE_INPUT[ifile] );
    }
  }

  // --------------------------------------------------
  //                 MEMORY ALLOC
  //
  // count how many string-variables: needed to allocate memory
  NVARSTR_FILE = 0 ;
  for ( ivar=1; ivar <= NVARALL_FILE; ivar++ ) {    
    if ( SKIP_VARNAME(ifile, ivar) ) { continue ; }
    if ( READTABLE_POINTERS.ICAST_STORE[ivar-1] == ICAST_C ) 
      { NVARSTR_FILE++ ; }
  }

  // get approx number of SN for memory allocation

  NEVT_APPROX = SNTABLE_NEVT_APPROX_TEXT(FFILE_INPUT[ifile], NVARALL_FILE);

  if ( NEVT_APPROX >= MXSN-1 ) { NEVT_APPROX = MXSN-1 ; }


  // Dec 2013: if 2nd file is much smaller, avoid too small malloc
  if ( ifile > 1 && NEVT_APPROX < NLIST_FIRST_FITRES ) {
    NEVT_APPROX = NLIST_FIRST_FITRES + 2 ;
  }


  // allocate new memory each time
  fitres_malloc_flt(ifile, NVARALL_FILE, NEVT_APPROX);
  fitres_malloc_str(ifile, NVARSTR_FILE, NEVT_APPROX);

  // ----------------------------------------------------

  ivarstr = 0 ;

  for ( ivar=1; ivar <= NVARALL_FILE; ivar++ ) {

    // don't duplicate FIELD (Dec 8 2014)
    if ( SKIP_VARNAME(ifile, ivar) ) { continue ; }

    // get VARNAME and ICAST
    VARNAME = READTABLE_POINTERS.VARNAME[ivar-1] ;
    ICAST   = READTABLE_POINTERS.ICAST_STORE[ivar-1] ;

    NVARALL_FITRES++ ;
    NVAR = NVARALL_FITRES ;

    if ( ICAST == ICAST_C ) {
      NVARSTR_FITRES++ ;   // summed over fitres files
      ivarstr++ ;          // sum for each separate fitres file
      sprintf(VARNAME_C, "%s:C*20", VARNAME) ;  //specify char ptr
      index = SNTABLE_READPREP_VARDEF(VARNAME_C, 
			     &FITRES_VALUES.STR_TMP[ivarstr][1], MXSN, 1 ) ;

    }
    else {     
      sprintf(VARNAME_F, "%s:F", VARNAME) ;  // specify float pointer
      index = SNTABLE_READPREP_VARDEF(VARNAME_F, 
			     &FITRES_VALUES.FLT_TMP[ivar][1], MXSN, 1 ) ; 
    }


    if ( NVAR >= MXVAR_TOT ) {
      sprintf(c1err,"NVARALL_COMBINE=%d exceeds arround bound of", NVAR);
      sprintf(c2err,"MXVAR_TOT = %d ", MXVAR_TOT ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }      
    ptr_CTAG = VARNAME_COMBINE[NVAR] ;
    sprintf(ptr_CTAG, "%s", VARNAME );

    // use CID (1st var) in 1st file only
    REPEATCID = ( ifile > 1 && ivar == 1 ) ;
    if ( REPEATCID ) { 
      ICAST_FITRES_COMBINE[NVAR] = -1 ;  // don't write repeated CIDs
      continue ;
    } 
    else { 
      ICAST_FITRES_COMBINE[NVAR] = ICAST ;
      NVAR_WRITE_COMBINED++ ; 
    } 
    
    // check if this tag name already exists;
    // if so, append '_$ifile' to end of tagname

    NTAG_DEJA = NMATCH_VARNAME(ptr_CTAG, NVAR-1) ;
    if ( NTAG_DEJA > 0 ) {
      printf("\t ADD_FITRES WARNING: VARNAME=%s already exists: ", ptr_CTAG);
      printf("append %d \n", ifile );
      sprintf( VARNAME_COMBINE[NVAR], "%s_%d", ptr_CTAG, ifile );
    }

    fflush(stdout);

  } // end of ivar loop


  // read everything.
  NLIST = SNTABLE_READ_EXEC();
  TABLEFILE_CLOSE(FFILE_INPUT[ifile]);

  if ( ifile == 1 ) { NLIST_FIRST_FITRES  = NLIST ; }

  // always fill 2nd NLIST2
  NLIST2_FITRES = NLIST ;

  printf("\n"); fflush(stdout);

  // ===============================================

  // isn2 is current SN index; ISN is SN index from  1st file
  // The logic here optimizes the speed of finding a cid match
  // between two lists ... CPU becomes noticable when N_SN 
  // is order 10^5.


  for ( isn2 = 1; isn2 <= NLIST2_FITRES; isn2++ ) {

    sprintf(ccid2, "%s", FITRES_VALUES.STR_TMP[IVARSTR_CCID][isn2] );   

    // first find "ISN" match in first list; used to load NTUP array

    if ( ifile == 1 ) { ISN = isn2; goto FOUND_ISN ; }

    ISN = -9;
    for ( isn = 1; isn <= NLIST_FIRST_FITRES ; isn++ ) {

      // checking logical is faster than checking string
      if ( USEDCID[isn] ) { continue ; }

      sprintf(ccid,"%s", FITRES_VALUES.STR_ALL[IVARSTR_CCID][isn] );

      if ( strcmp(ccid,ccid2) == 0 ) {
	ISN          = isn ;
	USEDCID[isn] = 1;
	goto FOUND_ISN ;
      }
    }

    if ( ISN <= 0 ) { continue ; }

  FOUND_ISN:

    LTMP = 0 ;

    // check for screen-update

    float TMPMOD ;
    TMPMOD = fmodf( (float)(ISN), 10000. ) ;
    if ( TMPMOD == 0.0        )    { LTMP = 1; }
    if ( ISN <= MXUPDATE      )    { LTMP = 1; }
    if ( ISN == NLIST_FIRST_FITRES  )    { LTMP = 1; }

    if ( ifile == 1 && LTMP == 1 )
      { printf("\t %s = '%12s'  ->  ISN = %6d \n", 
	       VARNAME_COMBINE[1], ccid2, ISN ); }
    
    if ( ISN == MXUPDATE+1 ) { printf("\t ... \n"); }

    ivarstr = 0;
    IVARTOT = NVARALL_LAST ;

    for ( ivar=1; ivar <= NVARALL_FILE; ivar++ ) {

      // Dec 2014
      if ( SKIP_VARNAME(ifile, ivar) ) { continue ; }

      VARNAME = READTABLE_POINTERS.VARNAME[ivar-1] ;
      ICAST   = READTABLE_POINTERS.ICAST_STORE[ivar-1] ;      
      IVARTOT++ ;

      if ( ICAST != ICAST_C )  {   // not a string

	if (  isnan(FITRES_VALUES.FLT_TMP[ivar][isn2]) !=0 ) {
	  sprintf(c1err,"isnan for FLT_TMP[ivar=%d][isn=%d]", ivar, isn2 );
	  sprintf(c2err,"varname = %s", VARNAME_COMBINE[ivar] );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	}
	
	FITRES_VALUES.FLT_ALL[IVARTOT][ISN] = 
	  FITRES_VALUES.FLT_TMP[ivar][isn2] ; 
      }
      else {
	ivarstr++ ;
	IVARSTR = NVARSTR_LAST + ivarstr ;
	IVARSTR_STORE[IVARTOT] = IVARSTR ;

	sprintf(FITRES_VALUES.STR_ALL[IVARSTR][ISN],"%s", 
		FITRES_VALUES.STR_TMP[ivarstr][isn2] ); 
      }
 
    } // ivar


  } // end of isn2 loop

  fflush(stdout);

  // free temp arrays
  freeVar_TMP(ifile, NVARALL_FILE, NVARSTR_FILE, NEVT_APPROX);

} // end of ADD_FITRES


// =====================================
int SKIP_VARNAME(int ifile, int ivar) {

  // Dec 8 2014
  // Return 1 if this variable should be ignored.
  
  char *VARNAME ;
  int  j;

  // no SNANA file yet
  if ( IFILE_FIRST_SNANA < 0 ) { return(0); } 

  // don't skip 1st SNANA file.
  if ( ifile == IFILE_FIRST_SNANA ) { return(0) ; }  

  VARNAME = READTABLE_POINTERS.VARNAME[ivar-1] ;
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

  for ( ivar=1; ivar <= NVARTOT; ivar++ ) {

    if ( ivar <= NVARSTR ) {

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

  int ivar, isn, IVAR_ALL, NTOT, MEMF, MEMD ;
  char fnam[] = "fitres_malloc_flt" ;

  // ---------- BEGIN ------------

  printf("\t Allocate fitres memory for NVAR=%d and MAXLEN=%d (ifile=%d)\n",
	 NVAR, MAXLEN, ifile );

  // redo malloc on TMP arrays for each fitres file
  MEMF      = (NVAR+1) * sizeof(float*) ;
  MEMD      = (NVAR+1) * sizeof(double*) ;
  FITRES_VALUES.FLT_TMP = (float **)malloc(MEMF) ;

  // -----------------------------------
  NTOT = NVARALL_FITRES + NVAR + 1 ;
  MEMF = sizeof(float*) * NTOT ;

  if (ifile == 1 ) 
    { FITRES_VALUES.FLT_ALL = (float**)malloc(MEMF); }
  else  
    { FITRES_VALUES.FLT_ALL = (float**)realloc(FITRES_VALUES.FLT_ALL,MEMF); }

  
  for ( ivar=1; ivar <= NVAR; ivar++ ) {

    MEMF = sizeof(float  ) * MAXLEN ;
    MEMD = sizeof(double ) * MAXLEN ;
    IVAR_ALL = NVARALL_FITRES + ivar ;

    FITRES_VALUES.FLT_TMP[ivar]     = (float  *)malloc(MEMF);
    FITRES_VALUES.FLT_ALL[IVAR_ALL] = (float  *)malloc(MEMF);    

    for ( isn=0; isn < MAXLEN; isn++ ) {
      USEDCID[isn] = 0 ;
      FITRES_VALUES.FLT_TMP[ivar][isn]     = INIVAL_COMBINE ;
      FITRES_VALUES.FLT_ALL[IVAR_ALL][isn] = INIVAL_COMBINE ;
    }
  }


} // end of fitres_malloc_flt




// =====================================
void  fitres_malloc_str(int ifile, int NVAR, int MAXLEN) {

  // allocate memory for string-variables
  // NVAR is the number of string variables in this fitres file.
  // Note that NVAR >= 1 because the CID string must always
  // be there.

  char fnam[] = "fitres_malloc_str" ;
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


  if (ifile == 1 ) 
    { FITRES_VALUES.STR_ALL = (char***)malloc(MEMC) ; }
  else 
    { FITRES_VALUES.STR_ALL = (char***)realloc(FITRES_VALUES.STR_ALL,MEMC); }

  
  for ( ivar=1; ivar <= NVAR; ivar++ ) {

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
    } // isn
  }  // ivar


} // end of fitres_malloc_str


// ===================================
int NMATCH_VARNAME(char *ctag , int ntlist ) {

  // Returns number of *ctag matches in first "ntlist"
  // elements of CTAG array
  

  int i, NMATCH, OVP ;

  NMATCH = 0;

  for ( i=1; i <= ntlist; i++ ) {
    OVP    = strcmp(ctag,VARNAME_COMBINE[i]) ;
    if ( OVP == 0 )  { NMATCH++ ; }
  }

  return NMATCH;

} // end of NMATCH_VARNAME



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
    ,BLOCKID[]    = "ID" 
    ,BLOCKVAR[]   = "VAR"     
    ,fnam[] = "WRITE_SNTABLE" 
    ,openOpt[40], CCIDint[40]
    ;

  int ivar, ivarstr, isn, IERR, ICAST, LENC, ic, CIDint ;
  int IFILETYPE, NOUT ;

  // --------------- BEGIN ------------

  IERR = -9 ;
  NOUT = 0 ;
  sprintf(OUTFILE[NOUT],"");

#ifdef USE_HBOOK
  if (USEFILE_HBOOK)  { 
    sprintf(OUTFILE[NOUT], "%s.%s", OUTPREFIX_COMBINE, ptrSuffix_hbook );  
    sprintf(openOpt,"%s new", ptrSuffix_hbook);
    IFILETYPE = TABLEFILE_OPEN(OUTFILE[NOUT],openOpt);
    NOUT++ ;
  }
#endif


#ifdef USE_ROOT
  if ( USEFILE_ROOT )  { 
    sprintf(OUTFILE[NOUT], "%s.%s", OUTPREFIX_COMBINE, ptrSuffix_root ); 
    sprintf(openOpt,"%s new", ptrSuffix_root);
    IFILETYPE = TABLEFILE_OPEN(OUTFILE[NOUT],openOpt);
    NOUT++ ;
  }
#endif


#ifdef USE_TEXT
  if ( USEFILE_TEXT )  { 
    sprintf(OUTFILE[NOUT], "%s.%s", OUTPREFIX_COMBINE, ptrSuffix_text ); 
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
  int CIDint_EXISTS=0 ;
  for ( ivar=1; ivar <= NVARALL_FITRES; ivar++ ) {
    if ( strcmp(VARNAME_COMBINE[ivar],"CIDint")==0 ) { CIDint_EXISTS=1; }
  }


  for ( ivar=1; ivar <= NVARALL_FITRES; ivar++ ) {

    ICAST = ICAST_FITRES_COMBINE[ivar] ;
    if ( ICAST < 0 ) {  continue ; }

    if ( ICAST == ICAST_C ) {
      int MXLEN ;
      MXLEN = maxlen_varString(VARNAME_COMBINE[ivar]);

      // include char CCID & integer CID representation for CID
      if ( ivar == 1 ) {

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

  for ( isn = 1; isn <= NLIST_FIRST_FITRES ; isn++ ) {

    TABLEROW_VALUES.CIDint = -999;
    sprintf(TABLEROW_VALUES.CCID,"xxx ");

    for ( ivar=1; ivar <= NVARALL_FITRES; ivar++ ) {

      ICAST = ICAST_FITRES_COMBINE[ivar] ;
      if ( ICAST < 0 ) { continue ; }

      if ( ICAST == ICAST_C )  { 
	ivarstr = IVARSTR_STORE[ivar];
	ptrSTR  = TABLEROW_VALUES.STR[ivarstr] ;

	sprintf(ptrSTR,"%s", FITRES_VALUES.STR_ALL[ivarstr][isn]);

	// xyz
	if ( ivar == 1 ) {
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

	if ( USEFILE_HBOOK ) 
	  { remove_string_termination(ptrSTR, MXSTRLEN); }

      }
      else
	{ TABLEROW_VALUES.FLT[ivar] = FITRES_VALUES.FLT_ALL[ivar][isn];  }
      
    } // ivar

    SNTABLE_FILL(TABLEID_COMBINE);

  } // isn

  // close it
  int out;
  for(out=0; out < NOUT; out++ ) 
    { TABLEFILE_CLOSE(OUTFILE[out]);  }

} // end of WRITE_SNTABLE


// ==================================
void  ADD_SNTABLE_COMMENTS(void) {

  // Call STORE_TABLEFILE_COMMENT(comment)

  int ifile ;
  char comment[120];

  // ------------ BEGIN ------------
  sprintf(comment,"Created by combine_fitres.exe");
  STORE_TABLEFILE_COMMENT(comment) ;

  sprintf(comment,"Combined files:", "%s", NFFILE_INPUT);
  STORE_TABLEFILE_COMMENT(comment) ;

  for(ifile=1; ifile <= NFFILE_INPUT; ifile++ ) {
    sprintf(comment,"\t + %s", FFILE_INPUT[ifile] );
    STORE_TABLEFILE_COMMENT(comment) ;
  }

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

