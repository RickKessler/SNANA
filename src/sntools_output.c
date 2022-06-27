/************************************************
 Created Feb 2013, R.Kessler


 Interface tools for snana output
  SNTABLE:  store table structure(s), e.g., hbook ntuple, root tree, etc ..
  SNHIST:   1D and 2D histograms
  SNLCPAK:  pack SN ligthcurve and best-fit model(s) for plotting
  SPECPAK:  pack spectra (4/2019)

 Function tree, where MAIN refers to main program
 calling functions in this file.

 MAIN  (see SNFILE_OUT_INIT in snana.car)
   --> TABLEFILE_INIT()     ! required
   --> SNLCPAK_INIT(...)    ! required for SNLCPAK
   --> TABLEFILE_OPEN(...)  ! required
   --> SNLCPAK_SURVEY()     ! required for SNLCPAK

       The above routines are called only from snana.car


  The following CCID loop is written specific to each program
  (snana.exe, snlc_fit.exe, psnid.exe)

  foreach CCID {
     SNLCPAK_NFIT(NFIT) ;  // only for multiple fits per SN
     MAKEDIR_OUTPUT(CCID);

     (fill STRING arrays and [DATA] arrays)

     SNLCPAK_DISPLAYTEXT(STRING1) ;  // text to display avove plots
     SNLCPAK_DISPLAYTEXT(STRING2) ;
     SNLCPAK_DISPLAYTEXT(STRING3) ;
     
     SNLCPAK_DATA(CCID, [DATA], SNLCPAK_EPFLAG_FLUXDATA );  // data array
     SNLCPAK_DATA(CCID, [DATA], SNLCPAK_EPFLAG_REJECT );  // exluded from fit
     SNLCPAK_DATA(CCID, [DATA], SNLCPAK_EPFLAG_CHI2   );  // data-fit chi2
     SNLCPAK_DATA(CCID, [DATA], SNLCPAK_EPFLAG_FITFUN );  // smooth fit-fun

     SNLCPAK_FILL();  // store [DATA] in output file(s).

     CDTOPDIR_OUTPUT();
  }

 To make an optional table (ntuple-like object)
     SNTABLE_CREATE(ID, NAME);
     SNTABLE_ADDCOL(ID, BLOCK, ptrVar_1, VARNAMES_1 );
     SNTABLE_ADDCOL(ID, BLOCK, ptrVar_2, VARNAMES_2 );
     ...
     SNTABLE_ADDCOL(ID, BLOCK, ptrVar_N, VARNAMES_N );
     SNTABLE_FILL(ID); 

  To make histogram,
     SNHIST_INIT(...)
     SNHIST_FILL(...)


  HISTORY

 Mar 28 2016:
   + check NVAR_ADDCOL_TOT < MXVAR_TABLE; abort otherwise

 Sep 12 2016: allow reading column into 1 or 2 different arrays.
              See .NPTR[ivar]. Need by SALT2mu to read some info
              into redundant CUTWIN array.

 Jan 27 2020: abort of no extension for root or hbook file;
              see TABLEFILE_OPEN.

 May 04 2020: add MARZ output option.

 May 30 2020: include sndata.h and remove a few redundant define statements
               in sntools_outout.h

************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <sys/stat.h>

// #include "sntools.h"
#include "sndata.h"
#include "sntools_output.h"

// include the package-specific code(s) here.
#ifdef USE_HBOOK
#include "sntools_output_hbook.c"
#endif

#ifdef USE_ROOT
#include "sntools_output_root.c"
#endif

#ifdef USE_TEXT
#include "sntools_output_text.c"
#endif

#ifdef USE_MARZ
#include "sntools_output_marz.c"
#endif


// ===============================================
void SNTABLE_DEBUG_DUMP(char *fnam, int idump) {
  printf(" xxx DEBUG_DUMP(%s , %d): SNLCPAK_USE[HBOOK/ROOT]=%d/%d \n",
	 fnam, idump, SNLCPAK_USE_HBOOK, SNLCPAK_USE_ROOT);
  fflush(stdout);
} // end of   void SNTABLE_DEBUG_DUMP
void sntable_debug_dump__(char *fnam, int *idump) 
{ SNTABLE_DEBUG_DUMP(fnam,*idump); }
 

// ===============================================
void TABLEFILE_INIT(void) {

  // Created April 2014
  // general init of global flags.  
  // This function must be called once before any
  // other function below. If you forget, then the
  // next function call will abort.

  int o,t ;
  char *s ;
  char U[] = "UNKNOWN" ;

  // -------------- BEGIN --------------

  if ( CALLED_TABLEFILE_INIT == 740 ) { return ; } // Jan 11, 2017
  CALLED_TABLEFILE_INIT = 740 ; // anything but 0
  NOPEN_TABLEFILE = 0 ;
  NFILE_AUTOSTORE = 0 ;
  NREAD_AUTOSTORE = 0 ;
  // -------------------------------

  // ------ misc inits -------
  OUTLIER_INFO.USEFLAG = 0 ;
  NLINE_TABLECOMMENT = 0 ;

  for(o=0; o < MXOPENFLAG; o++ ) {
    sprintf(STRING_TABLEFILE_OPENFLAG[o],"%s", U) ; 
    for(t=0; t < MXTABLEFILETYPE; t++ ) {

      if ( o == 0 ) {
	FIRSTCALL_TABLEFILE_OPEN[t] = 1; // init to true
	sprintf(STRING_TABLEFILE_TYPE[t], "%s", U) ;
	sprintf(STRING_IDTABLE_SNANA[t],  "%s", U) ;
	sprintf(STRING_IDTABLE_FITRES[t], "%s", U) ;
	sprintf(STRING_IDTABLE_OUTLIER[t], "%s", U) ;
      } 

      NAME_TABLEFILE[o][t][0] = 0 ;
      USE_TABLEFILE[o][t] = 0 ;
    }
  }

  get_SNANA_VERSION(SNANA_VERSION); // Dec 10 2017

  ADDCOL_VARLIST_LAST[0] = 0 ;

  s = STRING_TABLEFILE_TYPE[IFILETYPE_NULL]  ;  sprintf(s,"NULL");
  s = STRING_TABLEFILE_TYPE[IFILETYPE_HBOOK] ;  sprintf(s,"HBOOK");
  s = STRING_TABLEFILE_TYPE[IFILETYPE_ROOT]  ;  sprintf(s,"ROOT");
  s = STRING_TABLEFILE_TYPE[IFILETYPE_TEXT]  ;  sprintf(s,"TEXT");

  s = STRING_TABLEFILE_OPENFLAG[OPENFLAG_NULL]  ;  sprintf(s,"NULL");
  s = STRING_TABLEFILE_OPENFLAG[OPENFLAG_NEW]   ;  sprintf(s,"NEW" );
  s = STRING_TABLEFILE_OPENFLAG[OPENFLAG_READ]  ;  sprintf(s,"READ");

  s = STRING_IDTABLE_SNANA[IFILETYPE_HBOOK] ; sprintf(s,"7100");
  s = STRING_IDTABLE_SNANA[IFILETYPE_ROOT]  ; sprintf(s,"SNANA");
  s = STRING_IDTABLE_SNANA[IFILETYPE_TEXT]  ; sprintf(s,"SNANA");

  s = STRING_IDTABLE_FITRES[IFILETYPE_HBOOK] ; sprintf(s,"7788"  );
  s = STRING_IDTABLE_FITRES[IFILETYPE_ROOT]  ; sprintf(s,"FITRES");
  s = STRING_IDTABLE_FITRES[IFILETYPE_TEXT]  ; sprintf(s,"FITRES");

  s = STRING_IDTABLE_OUTLIER[IFILETYPE_HBOOK] ; sprintf(s,"7800"  );
  s = STRING_IDTABLE_OUTLIER[IFILETYPE_ROOT]  ; sprintf(s,"OUTLIER");
  s = STRING_IDTABLE_OUTLIER[IFILETYPE_TEXT]  ; sprintf(s,"OUTLIER");

  // useful string for cast manipulations
  sprintf(CCAST_TABLEVAR," CISF---D-------L--" );

  // misc. for SNLCPAK & SPECPAK

  SNLCPAK_USE_HBOOK = SNLCPAK_USE_ROOT = SNLCPAK_USE_TEXT = 0 ;  
  SPECPAK_USE_HBOOK = SPECPAK_USE_ROOT = SPECPAK_USE_TEXT = 0 ;
  SPECPAK_USE_MARZ  = 0;

  sprintf(SNLCPAK_OUTPUT.SURVEY,             "NULL" );
  sprintf(SNLCPAK_OUTPUT.VERSION_PHOTOMETRY, "NULL" );
  sprintf(SNLCPAK_OUTPUT.VERSION_SNANA,      "NULL" );
  sprintf(SNLCPAK_OUTPUT.SURVEY_FILTERS,     "NULL" );
  SNLCPAK_OUTPUT.TEXT_FORMAT[0] = 0;

#ifdef USE_TEXT
  FILEPREFIX_TEXT[0] = 0 ;
#endif

} // end of TABLEFILE_INIT

void tablefile_init__(void) {  TABLEFILE_INIT();  }

// =====================================
int get_TABLEFILE_TYPE(char *FILENAME) {

#ifdef USE_HBOOK
  if ( ISFILE_HBOOK(FILENAME) ) { return IFILETYPE_HBOOK; }
#endif

#ifdef USE_ROOT
  if ( ISFILE_ROOT (FILENAME) ) { return IFILETYPE_ROOT ; }
#endif

#ifdef USE_ROOT
  if ( ISFILE_TEXT (FILENAME) ) { return IFILETYPE_TEXT ; }
#endif

#ifdef USE_MARZ
  if ( ISFILE_MARZ (FILENAME) ) { return IFILETYPE_MARZ ; }
#endif

  // if we get here we did not find a valid type so return null value
  return -9;

} // end of TABLEFILE_TYPE

// ==================================================
int TABLEFILE_OPEN(char *FILENAME, char *STRINGOPT) {

  // Created Apr 26 2014 by R.Kessler
  // Determine file type (hbook, root, ...) based on suffix
  // of file name, then open it. 
  //
  // STRINGOPT is a list of option-keys:
  // -  new   -> open and create new file for writing
  // -  read  -> open and existing file for readonly
  // -  q     -> quiet mode; don't print anything to screen
  // -  root or hbook or text -> use explicit file type; ignore suffix.
  //
  // and note that all string-options are case-insensitive so
  // that new or NEW will work, q or Q, etc ...
  //
  // STRINGOPT Examples:
  //   STRINGOPT = 'new'       ! open new file, use suffix for file-type
  //   STRINGOPT = 'READ q'    ! open existing file, quiet mode
  //   STRINGOPT = 'NEW ROOT'  ! open new root file, ignore suffix.
  //
  // Finally, at any given time 1 NEW and 1 OLD file can be
  // opened for each file type so that reading and writing
  // can be done simultaneously.
  // Also, 1 hbook and 1 root file can be opened for output;
  // but cannot open 2 output-hbook files or 2 output-root files.
  //
  // Oct 14 2014: call new function OPEN_TEXTFILE(...) for read-mode
  // Jul 13 2020: declare *ENV and *FMT (used if HBOOK is NOT defined)

  int  OPEN_FLAG, TYPE_FLAG, OPT_Q, USE_CURRENT, IERR ;
  char *ptrtok, local_STRINGOPT[80], ctmp[20], *FMT, ENV[200] ;
  char fnam[] = "TABLEFILE_OPEN" ;

  // ---------------------- BEGIN ---------------------

  TABLEFILE_INIT_VERIFY(fnam, FILENAME) ;

  OPEN_FLAG = TYPE_FLAG = OPT_Q = IERR = 0 ;
  
  // define allowed keys in STRINGOPT
  char key_new[]   = "new" ;  // --> write
  char key_read[]  = "read" ;
  char key_q[]     = "q" ;

  // define extensions
  char key_root[]  = "root";
  char key_hbook[] = "hbook" ;
  char key_text[]  = "text" ;

  sprintf(local_STRINGOPT,"%s", STRINGOPT);
  ptrtok = strtok(local_STRINGOPT," "); // split string

  while ( ptrtok != NULL ) {
    sprintf(ctmp, "%s", ptrtok );

    if ( strcmp_ignoreCase(ctmp,key_new)  == 0 ) 
      { OPEN_FLAG  = OPENFLAG_NEW ; }

    else if ( strcmp_ignoreCase(ctmp,key_read) == 0 ) 
      { OPEN_FLAG = OPENFLAG_READ ;  }

    else if ( strcmp_ignoreCase(ctmp,key_q)    == 0 ) 
      { OPT_Q    = 1; }  // quiet mode

    else if ( strcmp_ignoreCase(ctmp,key_root) == 0 ) 
      { TYPE_FLAG = IFILETYPE_ROOT ; }

    else if ( strcmp_ignoreCase(ctmp,key_hbook) == 0 ) 
      { TYPE_FLAG = IFILETYPE_HBOOK ; }

    else if ( strcmp_ignoreCase(ctmp,key_text) == 0 ) 
      { TYPE_FLAG = IFILETYPE_TEXT ; }

    else {
      sprintf(MSGERR1,"Invalid option '%s'", ctmp);
      sprintf(MSGERR2,"in STRINGOPT = '%s' ", STRINGOPT);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }

    ptrtok = strtok(NULL," " );
  }


  
  // make sure table file has some kind of extension (Jan 2020)
  if ( TYPE_FLAG==IFILETYPE_ROOT || TYPE_FLAG==IFILETYPE_HBOOK ||
       TYPE_FLAG==IFILETYPE_MARZ ) {
    if ( strchr(FILENAME,'.') == NULL ) {
      sprintf(MSGERR1,"Missing extension for FILENAME = '%s' ", FILENAME);
      sprintf(MSGERR2,"Add valid extension: "
	      "e.g., '.ROOT', '.HBOOK', '.TEXT', '.fits' " );
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }
  }


  // --------------------------------------------------
  // if file-type is not passed explicitly, then  analyze suffix.


  if ( TYPE_FLAG == 0 ) {   
#ifdef USE_HBOOK
    if ( ISFILE_HBOOK(FILENAME) )  { TYPE_FLAG = IFILETYPE_HBOOK ; }
    if ( TYPE_FLAG > 0 ) { goto ISFILE_DONE ; }
#endif

#ifdef USE_ROOT
    if ( ISFILE_ROOT(FILENAME) )  {  TYPE_FLAG = IFILETYPE_ROOT ; }
    if ( TYPE_FLAG > 0 ) { goto ISFILE_DONE ; }
#endif

#ifdef USE_TEXT
    // Oct 2014: this works for read-mode only because for writing
    // FILENAME is a prefix.
    if ( ISFILE_TEXT(FILENAME) )  {  TYPE_FLAG = IFILETYPE_TEXT ; }
    if ( TYPE_FLAG > 0 ) { goto ISFILE_DONE ; }
#endif

#ifdef USE_MARZ
    if ( ISFILE_MARZ(FILENAME) )  {  TYPE_FLAG = IFILETYPE_MARZ ; }
    if ( TYPE_FLAG > 0 ) { goto ISFILE_DONE ; }
#endif

  } // end of OPT_TYPE==0 check

  
 ISFILE_DONE:

  if( OPEN_FLAG > 0 && TYPE_FLAG > 0 ) 
    { USE_CURRENT = USE_TABLEFILE[OPEN_FLAG][TYPE_FLAG];  }
  else
    { USE_CURRENT = 0 ; }

  /*
  printf(" xxx %s: OPEN_FLAG=%d  TYPE_FLAG=%d  (%s) \n",
	 fnam, OPEN_FLAG, TYPE_FLAG, FILENAME);
  */

  // ---------------------------------------
  // check for fatal errors -> abort

  if ( OPEN_FLAG == 0 ) {
    print_preAbort_banner(fnam);
    printf("  file=\n %s\n", FILENAME);
    sprintf(MSGERR1,"Must specify 'new' or 'read' i STRINGOPT.");
    sprintf(MSGERR2,"Check STRINGTOP ='%s' argument.", STRINGOPT);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  // abort if we still do not know the file type.
  if ( TYPE_FLAG == 0 ) {   
    print_preAbort_banner(fnam);
    printf("  file=\n %s\n", FILENAME);
    sprintf(MSGERR1,"Unknown table-file type");
    sprintf(MSGERR2,"STRINGOPT = '%s' ", STRINGOPT);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);     
  }

  if ( USE_CURRENT ) {
    print_preAbort_banner(fnam);
    printf("   file=\n %s\n", FILENAME);
    sprintf(MSGERR1,"%s %s file alread used.",
	    STRING_TABLEFILE_OPENFLAG[OPEN_FLAG],
	    STRING_TABLEFILE_TYPE[TYPE_FLAG] );
    sprintf(MSGERR2,"STRINGOPT='%s'", STRINGOPT);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  // -----------------------------------------------------------
  // Determine stringOpt to pass to type-specific open function.

  char stringOpt[20] ;
  stringOpt[0] = 0 ;
  if ( OPEN_FLAG == OPENFLAG_NEW ) 
    { sprintf(stringOpt,"%sN", stringOpt); }

  if ( OPT_Q  ) 
    { sprintf(stringOpt,"%sQ", stringOpt); }

  // open file based on its type.

 
#ifdef USE_HBOOK
  if ( TYPE_FLAG == IFILETYPE_HBOOK ) {
    int  LUN;
    char TOPDIR[20] ;    

    LUN = 40 + 10*OPENFLAG_NEW + TYPE_FLAG; 

    if ( OPEN_FLAG == OPENFLAG_NEW ) 
      { sprintf(TOPDIR,"%s", HBOOK_TOPDIRNAM_NEW ); }
    else
      { sprintf(TOPDIR,"%s", HBOOK_TOPDIRNAM_READ ); }
          
    OPEN_HFILE(FILENAME, LUN, stringOpt, TOPDIR, &IERR);
    NOPEN_TABLEFILE++ ;
  }
#endif

#ifndef USE_HBOOK
  if ( TYPE_FLAG == IFILETYPE_HBOOK ) {
    FMT = STRING_TABLEFILE_TYPE[TYPE_FLAG]; // "HBOOK"
    sprintf(ENV,"CERN_DIR");                // required ENV for make
    TABLEFILE_notCompiled_ABORT(FILENAME, FMT, ENV);  
  }
#endif


#ifdef USE_ROOT
  if ( TYPE_FLAG == IFILETYPE_ROOT ) {
    OPEN_ROOTFILE(FILENAME, stringOpt, &IERR);
    NOPEN_TABLEFILE++ ;
  }
#endif

#ifndef USE_ROOT
  if ( TYPE_FLAG == IFILETYPE_ROOT ) {
    FMT = STRING_TABLEFILE_TYPE[TYPE_FLAG];  // "ROOT"
    sprintf(ENV,"ROOT_DIR");                 // required ENV for make
    TABLEFILE_notCompiled_ABORT(FILENAME, FMT, ENV);  
  }
#endif

#ifdef USE_TEXT
  if ( TYPE_FLAG == IFILETYPE_TEXT ) {

    if ( OPEN_FLAG == OPENFLAG_NEW ) {
      char *PREFIX = FILENAME ;
      INIT_TEXTFILES(PREFIX) ; 
      NOPEN_TABLEFILE++ ;
      IERR = 0 ;
    }
    else {
      char mode[] = "rt" ;
      OPEN_TEXTFILE(FILENAME,mode); // read mode
    }
  }
#endif

#ifdef USE_MARZ
  if ( TYPE_FLAG == IFILETYPE_MARZ ) {
    OPEN_MARZFILE(FILENAME, &IERR);
    NOPEN_TABLEFILE++ ;
  }
#endif
  
  // store USE-flag and filename
  sprintf(NAME_TABLEFILE[OPEN_FLAG][TYPE_FLAG], "%s", FILENAME);
  USE_TABLEFILE[OPEN_FLAG][TYPE_FLAG] = 1; 

  if ( IERR != 0 ) 
    { return -9 ; }  // generaic error code
  else
    { return TYPE_FLAG ; }

} // end of TABLEFILE_OPEN

int  tablefile_open__(char *FILENAME, char *STRINGOPT) {
  return TABLEFILE_OPEN(FILENAME,STRINGOPT);
}
		      

// ===================================
void TABLEFILE_CLOSE(char *FILENAME) {

  // Apr 2014, R.Kessler
  // Close file based on its type (book or root) and how it was 
  // opened (new or read).  Note that input FILENAME here must 
  // match the name used to open the file (otherwise ABORT).
  //
  // Apr 2015: add CLOSE_TEXTFILE call.

  int o, t, USE, NMATCH, OPEN_FLAG, TYPE_FLAG ;
  char fnam[] = "TABLEFILE_CLOSE" ;

  // --------------- BEGIN ---------------

  // get open & type info for FILENAME by string-matching to 
  // store list of opened files.

  USE = NMATCH = OPEN_FLAG = TYPE_FLAG = 0 ;
  for(o=0; o < MXOPENFLAG; o++ ) {
    for(t=0; t < MXTABLEFILETYPE; t++ ) {
      if ( strcmp(NAME_TABLEFILE[o][t],FILENAME) == 0 ) {
	USE     = USE_TABLEFILE[o][t] ;
	NMATCH += 1;  
	OPEN_FLAG = o; TYPE_FLAG = t ;
      }
    }
  }

  
  if ( IGNOREFILE(FILENAME)  ) {
    TYPE_FLAG = IFILETYPE_TEXT ;
  }
  else {
    if ( NMATCH != 1 ) {    
      sprintf(MSGERR1,"Invalid NMATCH = %d  (USE=%d) for ", 
	      NMATCH,USE );
      sprintf(MSGERR2,"%s", FILENAME);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);     
    }
    
    if ( USE == 0 ) {
      sprintf(MSGERR1,"Something weird: found MATCH to unused file=");
      sprintf(MSGERR2,"%s", FILENAME);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);     
    }
  }

#ifdef USE_HBOOK
  if(TYPE_FLAG == IFILETYPE_HBOOK)  
    { CLOSE_HFILE(FILENAME,OPEN_FLAG);  }
#endif

#ifdef USE_ROOT
  if(TYPE_FLAG == IFILETYPE_ROOT ) 
    { CLOSE_ROOTFILE(FILENAME,OPEN_FLAG); }
#endif

#ifdef USE_TEXT
  if(TYPE_FLAG == IFILETYPE_TEXT ) 
    { CLOSE_TEXTFILE(); }
#endif

  
#ifdef USE_MARZ
  if(TYPE_FLAG == IFILETYPE_MARZ ) 
    { CLOSE_MARZFILE(FILENAME); }
#endif
 

  // ----------------------------------------------
  // reset info for this IO-flag and file-type;
  // e.g.., allows opening another read-only file.

  NAME_TABLEFILE[OPEN_FLAG][TYPE_FLAG][0] = 0 ;
  USE_TABLEFILE[OPEN_FLAG][TYPE_FLAG] = 0;

} // end of TABLEFILE_CLOSE


void tablefile_close__(char *FILENAME) { TABLEFILE_CLOSE(FILENAME); }


// =============================================================
void STORE_TABLEFILE_COMMENT(char *comment) {

  // Created Oct 23 2014
  // store global comment to write at top of table
  // (not part of table columns or rows)

  int N = NLINE_TABLECOMMENT ;
  sprintf(LINE_TABLECOMMENT[N], "%s", comment);
  NLINE_TABLECOMMENT++ ;

} // end of STORE_TABLEFILE_COMMENT

void store_tablefile_comment__(char *comment) 
{ STORE_TABLEFILE_COMMENT(comment); }

// =============================================================
void SNTABLE_CREATE(int IDTABLE, char *NAME, char *TEXT_FORMAT) {

  // Create/Initialize new table.
  // *IDTABLE and *NAME are integer and character udentifiers.
  // HBOOK must have an integer id, while ROOT uses a char ID.
  // TEXT_FORMAT used only for TEXT: 'key', 'csv', 'col'
  
  int  ID, USE ; 
  char fnam[] = "SNTABLE_CREATE" ;

  // -------- BEGIN ----------

  if ( NOPEN_TABLEFILE == 0 ) {
    char comment[60];
    sprintf(comment, "Cannot init IDTABLE=%d(%s)", IDTABLE, NAME);
    TABLEFILE_notOpen_ABORT(fnam, comment );
  }

  ID = IDTABLE ; // local ID has a real address

  NVAR_ADDCOL_TOT = 0; // Mar 28 2016

#ifdef USE_MARZ
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_MARZ] ;
  if ( USE && IDTABLE == TABLEID_MARZ ) 
    { SNTABLE_CREATE_MARZ(IDTABLE,NAME);  return; } 
#endif
  

#ifdef USE_HBOOK
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_HBOOK] ; 
  if ( USE ) { SNTABLE_CREATE_HBOOK(IDTABLE,NAME);  }
#endif

  if ( IDTABLE == 1600 ) { return ; } // skip special hbook-only table

#ifdef USE_ROOT
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_ROOT] ; 
  if ( USE ) { SNTABLE_CREATE_ROOT(IDTABLE,NAME);  }
#endif


#ifdef USE_TEXT
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_TEXT] ;
  if ( USE ) { SNTABLE_CREATE_TEXT(IDTABLE,NAME,TEXT_FORMAT);  } 
#endif

 
  fflush(stdout);

} // end of SNTABLE_CREATE

void sntable_create__(int *ID, char *NAME, char *TEXT_FORMAT) {
  SNTABLE_CREATE(*ID,NAME,TEXT_FORMAT)  ;
}


// ======================================
void SNTABLE_FILL(int IDTABLE ) {

  // Called for each even after all of the pointers are filled
  // with their values, this function fills the table with one 
  // more row.

  int USE ;
  //  char fnam[] = "SNTABLE_FILL" ;
  // -------- BEGIN ----------


#ifdef USE_HBOOK
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_HBOOK] ; 
  if ( USE ) { SNTABLE_FILL_HBOOK(IDTABLE);  }
#endif

  if ( IDTABLE == 1600 ) { return ; } // skip special hbook-only table
       
#ifdef USE_ROOT
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_ROOT] ; 
  if ( USE ) { SNTABLE_FILL_ROOT(IDTABLE); }
#endif


#ifdef USE_TEXT
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_TEXT] ; 
  if ( USE ) { SNTABLE_FILL_TEXT(IDTABLE); }
#endif


} // end of SNTABLE_FILL

void sntable_fill__(int *ID) {  SNTABLE_FILL(*ID);  }



// ======================================================
void SNTABLE_ADDCOL(int IDTABLE, char *BLOCK, void* PTRVAR, 
		    char *VARLIST, int USE4TEXT ) {

  // define new column(s) in table.
  // Inputs:
  // - IDTABLE :  integer specifier previously passed to SNTABLE_CREATE
  // - BLOCK   :  name of variable block (used for hbook, but not root)
  // - PTRVAR  :  pointer to variable that will be filled each event
  // - VARLIST :  comma-separated or space-separated list of variables 
  //              with cast;
  //              e.g.,  'RA:D, DEC:D',   'ITYPE:I',  'FLUX:F'
  //              Note that 'RA:D', and 'RA' are both double
  //              For list, must all have same cast.
  //
  // - USE4TEXT : logical flag for TEXT format since TEXT table keeps
  //              a subset. Ignored for ROOT and HBOOK.
  //
  //  To load a vector as a column element (e.g., fit resids per epoch),
  //  First call this function with
  //       VARLIST = 'NFIT[$maxsize]:I'
  //  and then with
  //       VARLIST = 'EPOCH_VAR1(NFIT):F'
  //       VARLIST = 'EPOCH_VAR2(NFIT):D'
  //       VARLIST = 'EPOCH_VAR1(NFIT):I'
  // etc ...
  // Note that $maxsize is an actual number giving the max size
  // of the vector. However, 'NFIT' is a string giving the name
  // of the array-size for each SN.

  int USE ;
  char fnam[] = "SNTABLE_ADDCOL" ;
  SNTABLE_ADDCOL_VARDEF ADDCOL_VARDEF ;

  // ---------------- BEGIN --------------
  
  parse_ADDCOL_VARLIST(VARLIST, &ADDCOL_VARDEF); // return ADDCOL_VARDEF

  //Mar 28 2016: add protection against too many variables
  NVAR_ADDCOL_TOT += ADDCOL_VARDEF.NVAR;
  if ( NVAR_ADDCOL_TOT >= MXVAR_TABLE ) {
    sprintf(MSGERR1,"NVAR_ADDCOL_TOT=%d exceeds MXVAR_TABLE=%d",
	    NVAR_ADDCOL_TOT, MXVAR_TABLE );
    sprintf(MSGERR2,"IDTABLE=%d", IDTABLE);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);     
  }

  
#ifdef USE_HBOOK
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_HBOOK] ; 
  if ( USE ) 
    { SNTABLE_ADDCOL_HBOOK(IDTABLE, BLOCK, PTRVAR, &ADDCOL_VARDEF); }
#endif

  if ( IDTABLE == 1600 ) { return ; } // skip special hbook-only table

#ifdef USE_ROOT
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_ROOT] ; 
  if ( USE ) 
    { SNTABLE_ADDCOL_ROOT(IDTABLE, PTRVAR, &ADDCOL_VARDEF); }
#endif


#ifdef USE_TEXT
  USE = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_TEXT] ; 
  if ( USE && USE4TEXT ) 
    { SNTABLE_ADDCOL_TEXT(IDTABLE, PTRVAR, &ADDCOL_VARDEF); }
#endif


} // end of SNTABLE_ADDCOL

void sntable_addcol__(int *ID, char *BLOCK, void* PTRVAR, 
		      char *VARLIST, int *USE4TEXT) {
  SNTABLE_ADDCOL(*ID, BLOCK, PTRVAR, VARLIST, *USE4TEXT );
}

// =====================================
void parse_ADDCOL_VARLIST(char *VARLIST,
			  SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF) {

  // Created May 2014
  // parse input VARLIST and fill output struct  *ADDCOL_VARDEF.
  // Each  VARLIST element is separated by a blank space or comma.
  //
  // Allow variable declarations in *VARLIST
  //
  // *  RA          -> implicit double
  // *  RA:D        -> declared double
  // *  FLUX:F      -> declared float
  // *  ITYPE:I     -> declared int
  // *  FIELD:C*20  -> 20 char 
  // *  NFIT[600]   -> vector array size, max size = 600
  // *  FLUX(NFIT)  -> vector of size NFIT
  //
  // and lists such as
  // *  'RA:D,DEC:D'
  // *  'SNRMAX1:F SNRMAX2:F SNRMAX3:F'
  //

  int NVAR, ICAST, ICAST_FIRST, VECFLAG, ISIZE, LDMP ;
  char local_VARLIST[MXCHAR_VARLIST];
  char *ptrtok, tmpVarName_withCast[80], tmpVarName[80] ;

  char fnam[] = "parse_ADDCOL_VARLIST" ;

  // ------------ BEGIN -------------

  LDMP = 0 ; // ( strlen(VARLIST) > 100) ;

  // store original [unparsed] VARLiST
  sprintf(ADDCOL_VARDEF->VARLIST_ORIG,"%s", VARLIST);

  NVAR = ICAST_FIRST = 0;
  sprintf(local_VARLIST,"%s", VARLIST);
  ptrtok = strtok(local_VARLIST,", "); // split string

  if ( LDMP ) {
    printf(" xxx ---------------------- \n" ) ;
    printf(" xxx VARLIST = '%s' \n", VARLIST);
    fflush(stdout);
  }

  while ( ptrtok != NULL ) {

    sprintf(tmpVarName_withCast, "%s", ptrtok );

    parse_TABLEVAR(tmpVarName_withCast,   // input
		   tmpVarName, &ICAST, &VECFLAG, &ISIZE); // returned

    if ( NVAR == 0 ) { ICAST_FIRST = ICAST ; }

    if ( ICAST != ICAST_FIRST ) {
      print_preAbort_banner(fnam);
      printf(" Invalid VARLIST = '%s' (ICAST_FIRST=%d  ICAST=%d) \n", 
	     VARLIST, ICAST_FIRST, ICAST );

      sprintf(MSGERR1,"Cannot change CAST in VARLIST.");
      sprintf(MSGERR2,"Invalid VARLIST printed above.");
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);         
    }

    if( NVAR < MXVAR_ADDCOL ) {
      // load global struct to be used by any SNTABLE implemenation
      sprintf(ADDCOL_VARDEF->VARNAME[NVAR],"%s", tmpVarName);
      sprintf(ADDCOL_VARDEF->CCAST[NVAR],  "%c", CCAST_TABLEVAR[ICAST] );
      ADDCOL_VARDEF->ICAST[NVAR] = ICAST ;
      ADDCOL_VARDEF->ISIZE[NVAR] = ISIZE ;
      ADDCOL_VARDEF->VECTOR_FLAG[NVAR] = VECFLAG ;
    }

    if ( LDMP ) {
      printf(" xxx VARNAME[%d]='%s'  ICAST=%d  ISIZE=%d  VECFLAG=%d\n",
	     NVAR, tmpVarName, ICAST, ISIZE, VECFLAG);
      fflush(stdout);
    }

    NVAR++ ;

    if ( NVAR > MXVAR_ADDCOL ) {
      sprintf(MSGERR1,"NVAR  exceeds bound of MXVAR_ADDCOL=%d for",
	      MXVAR_ADDCOL);
      sprintf(MSGERR2,"VARLIST = %s", VARLIST);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);   
    }

    ptrtok = strtok(NULL,", " );
  } // end while


  ADDCOL_VARDEF->NVAR = NVAR ;
  
  if ( NVAR <=0 ) {
    sprintf(MSGERR1, "NVAR=0 for VARLIST = '%s' ", 
	    ADDCOL_VARDEF->VARLIST_ORIG );
    sprintf(MSGERR2, "Last valid VARLIST = '%s'", ADDCOL_VARLIST_LAST);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  sprintf(ADDCOL_VARLIST_LAST,"%s", VARLIST); // for abort message only

} // end of parse_ADDCOL_VARLIST

// ===============================
void parse_TABLEVAR(char *varName_with_cast, char *varName, 
		    int *ICAST, int *VECFLAG, int *ISIZE) {

  // Created Apr 24 2014 by R.Kessler
  //
  // Parse string containing single var-name plus cast, 
  // and possible vector-information. Can use for writing or reading.
  //
  // Input : *varName_with_cast
  // Ouptut:
  //    *varName = variable name without cast or [] or ().
  //    *ICAST   = integer cast
  //    *VECFLAG = 0=scalar, 1=vec size, 2=vector
  //    *ISIZE   = char size (for char) or  vector length
  //
  //
  // Input *varName_with_cast can have optional casting suffix
  //  ':F', '/F', ':I', '/I', etc ...
  // If there is no suffix-cast, assume double precision.
  //  
  // Return function value ICAST and function argument
  // *varName with cast-suffix removed.
  //
  // Examples:
  //  *varName_with_cast = 'REDSHIFT' 
  //    -> varName = 'REDSHIFT'  and ICAST = ICAST_D
  //
  //  *varName_with_cast = 'REDSHIFT:F' 
  //    -> varName = 'REDSHIFT'  and ICAST = ICAST_F
  //
  //  *varName_with_cast = 'REDSHIFT:R'   [fortran real]
  //    -> varName = 'REDSHIFT'  and ICAST = ICAST_F
  //
  //  *varName_with_cast = 'ITYPE:/I' 
  //    -> varName = 'ITYPE'  and ICAST = ICAST_I
  //
  //  *varName_with_cast = 'ITYPE:/S' 
  //    -> varName = 'ITYPE'  and ICAST = ICAST_S
  //
  // May 7 2014:  for fortran specifier, allow
  //     ':R' and '/R' -> float
  //     ':8' --> double
  // 
  // Jan 8 2017: if char cast has no *20, then set isize=MXCHAR_CCID
  //             to avoid parsing bug.
  //
  // Jun 07 2019: add S for short int
  //

  int  icast, icast_tmp, vecFlag, isize ;
  int  i, ncp, ibr0, ibr1 ;
  char varName2[60], c1[2], cBRACKET[20] ;
  //  char fnam[] = "parse_TABLEVAR" ;

  // --------------- BEGIN ------------------

  varName[0] = 0 ;
  varName2[0] = 0 ;

  icast = vecFlag = isize = 0 ; // init outputs

  // in case of options, get varname2 = VARNAME without :
  ncp = strlen(varName_with_cast)  ;
  int ENDVAR = 0 ;
  ibr0 = ibr1 = 9999 ;  cBRACKET[0] = 0 ;

  for(i=0; i < ncp; i++ ) { 
    sprintf(c1, "%c", varName_with_cast[i] );
    if ( *c1 == ':' ) { ENDVAR = 1; }
    if ( *c1 == '/' ) { ENDVAR = 1; } // added Jun 9 2013 (bug-fix)
    if ( *c1 == '[' ) { ENDVAR = 1; ibr0  = i; }
    if ( *c1 == ']' ) { ENDVAR = 1; ibr1  = i; }

    // construct varName2 = name without cast or vector info
    if ( ENDVAR == 0 ) { strcat(varName2,c1); }

    // construct string-contents of []; should be integer array size
    if ( i > ibr0 && i < ibr1 )  { strcat(cBRACKET,c1); }

  } // end of i loop

  // check for /F or :F  or /I or :I  at end of string.

  if ( strstr(varName_with_cast,"/F") != NULL ) 
    { icast = ICAST_F ;  }
  else if ( strstr(varName_with_cast,":F") != NULL ) 
    { icast = ICAST_F ;  }

  else if ( strstr(varName_with_cast,"/R") != NULL ) 
    { icast = ICAST_F ;  }
  else if ( strstr(varName_with_cast,":R") != NULL ) 
    { icast = ICAST_F ;  }

  else if ( strstr(varName_with_cast,"/I") != NULL ) 
    { icast = ICAST_I ;  }
  else if ( strstr(varName_with_cast,":I") != NULL ) 
    { icast = ICAST_I ;  }

  else if ( strstr(varName_with_cast,":S") != NULL ) // short int
    { icast = ICAST_S ;  }

  else if ( strstr(varName_with_cast,"/L") != NULL ) 
    { icast = ICAST_L ; }
  else if ( strstr(varName_with_cast,":L") != NULL ) 
    { icast = ICAST_L ; }

  else if ( strstr(varName_with_cast,"/D") != NULL ) 
    { icast = ICAST_D ;  }
  else if ( strstr(varName_with_cast,":D") != NULL ) 
    { icast = ICAST_D ;  }

  else if ( strstr(varName_with_cast,"/8") != NULL ) 
    { icast = ICAST_D ;  }
  else if ( strstr(varName_with_cast,":8") != NULL ) 
    { icast = ICAST_D ;  }

  else if ( strstr(varName_with_cast,"/C") != NULL ) 
    { icast = ICAST_C ; }
  else if ( strstr(varName_with_cast,":C") != NULL ) 
    { icast = ICAST_C ; }


  if ( icast == ICAST_C ) {
    // get size of char array from :C*20 -> isize=20
    int lenvar = strlen(varName2);
    if ( lenvar+3 < ncp ) 
      { sscanf( &varName_with_cast[lenvar+3], "%d", &isize ); }
    else
      { isize = MXCHAR_CCID ; } 
      
  }

  // check for vector options :
  // [] --> this variable is an array size, return isize = inside []
  // () --> this variable is a vector ; set vecFlag and return it with ().
  if ( ibr0 < 9999 ) {
    vecFlag = 1 ;  
    icast   = ICAST_I ; // vector size must be integer
    sscanf(cBRACKET,"%d", &isize);
  }

  if ( strstr(varName_with_cast,"(") != NULL )  
    { vecFlag = 2; }
  
  // July 2017: if icast<0, check for char castings
  if ( icast <= 0 ) {
    icast_tmp = ICAST_for_textVar(varName_with_cast) ;
    // printf(" xxx icast_tmp=%d for '%s' \n", icast_tmp, varName_with_cast);
    if ( icast_tmp == ICAST_C ) { icast = ICAST_C; }
  }

  // ------------
  if ( icast > 0 ) 
    { sprintf(varName, "%s", varName2); }
  else { 
    // default is double
    sprintf(varName, "%s", varName_with_cast) ; 
    icast = ICAST_D ;  
  } 


  // load output args
  *ICAST   = icast ;
  *VECFLAG = vecFlag ;
  *ISIZE   = isize ;

  return ;

} // end of parse_TABLEVAR


// ==========================================
int SNTABLE_READPREP(int IFILETYPE, char *TABLENAME) {

  // Created Oct 14 2014
  // read list of variables and fill struct READTABLE_POINTERS
  // with NVARTOT and VARLIST array.
  // Functions returns number of variables in table.
  //
  // Inputs:
  //   IFILETYPE  : specifiies hbook, root, text ...
  //   TABLENAME  : name of table to read
  //
  // For hbook, translate table name text into ntuple id.

  char fnam[] = "SNTABLE_READPREP" ;
  int  NVAR, ISOPEN_DEJA, ivar ;
  char msg[200], TBNAME_LOCAL[40] ;
  // --------------- BEGIN --------------

  sprintf(msg,"read VARLIST from table = %s", TABLENAME);
  TABLEFILE_INIT_VERIFY(fnam, msg) ;

  // init pointers
  READTABLE_POINTERS.NVAR_TOT  = 0 ;
  READTABLE_POINTERS.NVAR_READ = 0 ;
  READTABLE_POINTERS.IFILETYPE = IFILETYPE_NULL ;
  READTABLE_POINTERS.NROW      = 0 ;
  READTABLE_POINTERS.FP_DUMP   = NULL ;
  for(ivar=0; ivar < MXVAR_TABLE; ivar++ ) {
    READTABLE_POINTERS.NPTR[ivar] = 0 ;
    sprintf(READTABLE_POINTERS.VARNAME[ivar], "unknown");
    READTABLE_POINTERS.ICAST_READ[ivar]     = -9 ;
    READTABLE_POINTERS.ICAST_STORE[ivar]    = -9 ;
    READTABLE_POINTERS.PTRINDEX[ivar]       = -9 ;
  }

  sprintf(TBNAME_LOCAL, "%s", TABLENAME);

  NVAR         = -9;

  ISOPEN_DEJA = USE_TABLEFILE[OPENFLAG_READ][IFILETYPE] ;
  if ( ISOPEN_DEJA == 0 ) {
    sprintf(MSGERR1,"Cannot read VARLIST from table=%s", TABLENAME);
    sprintf(MSGERR2,"because %s file is not open.",
	    STRING_TABLEFILE_TYPE[IFILETYPE] );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }


#ifdef USE_HBOOK
  int ID ;
  if ( IFILETYPE == IFILETYPE_HBOOK ) {
    ID = TABLEID_HBOOK(TABLENAME);
    sprintf(TBNAME_LOCAL, "%d", ID); //hbook table must be int ID
    NVAR = SNTABLE_READPREP_HBOOK(ID);
  }
#endif


#ifdef USE_ROOT
  if ( IFILETYPE == IFILETYPE_ROOT ) {
    NVAR = SNTABLE_READPREP_ROOT(TABLENAME); 
  }
#endif

#ifdef USE_TEXT
  if ( IFILETYPE == IFILETYPE_TEXT ) {
    NVAR = SNTABLE_READPREP_TEXT(); 
  }
#endif

  // store file type and name of table
  READTABLE_POINTERS.NVAR_TOT  = NVAR ;
  READTABLE_POINTERS.IFILETYPE = IFILETYPE ;
  sprintf(READTABLE_POINTERS.TABLENAME, "%s", TBNAME_LOCAL );

  // ---------------------------------------------  
  printf("   Read %d table varNames from %s table=%s \n",
	 NVAR, STRING_TABLEFILE_TYPE[IFILETYPE], TABLENAME );
  fflush(stdout);
  // -----------------

  if ( NVAR <= 0 ) {
    sprintf(MSGERR1,"NVAR=%d for IFILETYPE=%d ??? ", NVAR, IFILETYPE);
    sprintf(MSGERR2,"Table = %s ", TBNAME_LOCAL );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  return(NVAR);

} // end of SNTABLE_READPREP



// ========================================================
int SNTABLE_READPREP_VARDEF(char *VARLIST, void *ptr, 
			    int mxlen, int optMask) {

  // Oct 2014
  // Define and store pointer to read entire table column.
  // Check each element of *VARLIST to allow for ambiguous names;
  // e.g., VARLIST = 'z Z zcmb'.
  // If multiple names are given, only one is allowed to exist,
  // otherwise code aborts.
  //
  // This function is just a shell to examine VARLIST, find the
  // defined VARNAME, and call sntable_readprep_vardef1() with the 
  // single valid variable name.
  //
  //
  // Inputs:
  //   - ptr : pointer to store
  //   - mxlen = max length of input ptr
  //   - optMask: bit0 -> print for each var, bit1 -> abort on missing var
  //
  // Function returns absolute IVAR index.

  int  istat, ISTAT,  NVAR_TOT, NVAR_FOUND, FLAG_VBOSE, FLAG_ABORT ;
  char VARLIST_LOCAL[MXCHAR_VARLIST], VARLIST_FOUND[MXCHAR_VARLIST];
  char VARNAME_withCast[MXCHAR_VARNAME]; 
  char VARNAME_noCast[MXCHAR_VARNAME]; 
  char *ptrtok ;
  char fnam[] = "SNTABLE_READPREP_VARDEF" ;

  // ---------------- BEGIN -------------------

  NVAR_TOT = NVAR_FOUND = 0 ;
  ISTAT = -9 ;

  sprintf(VARLIST_LOCAL, "%s", VARLIST);
  VARLIST_FOUND[0] = 0 ;

  FLAG_VBOSE = (optMask & 1); // flag to print comment for each variable
  FLAG_ABORT = (optMask & 2); // flag to abort on missing variable

  ptrtok = strtok(VARLIST_LOCAL," ");
  while ( ptrtok != NULL ) {
    sprintf(VARNAME_withCast,"%s", ptrtok);
    istat=sntable_readprep_vardef1(VARNAME_withCast, ptr, mxlen, FLAG_VBOSE,
				     VARNAME_noCast );

    NVAR_TOT++ ;
    if ( istat >= 0 ) { 
      ISTAT = istat ;  
      sprintf(VARLIST_FOUND, "%s %s", VARLIST_FOUND, VARNAME_withCast);
      NVAR_FOUND++ ; 
    }
    ptrtok = strtok(NULL," " );
  }

  // error checking.

  // check flag to abort on missing variable
  if ( FLAG_ABORT && NVAR_FOUND == 0 ) {
    sprintf(MSGERR1,"Could not find variable among VARLIST='%s'",
	    VARLIST);
    sprintf(MSGERR2,"Check variables in table = %s",
	    READTABLE_POINTERS.TABLENAME );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
  }

  return ISTAT ;

} // end of SNTABLE_READPREP_VARDEF

// ===========================================================
int sntable_readprep_vardef1(char *varName_withCast, void *ptr, 
			     int mxlen, int vboseflag, 
			     char *varName_noCast) {

  // Oct 2014
  // do the dirty work described in SNTABLE_READPREP_VARDEF.
  // Find place of *varname and store *ptr for reading
  // Returns index of *varname if *varname exists;
  // returns -1 otherwise. 
  //
  // Inputs:
  //  varName_withCast : name of variable including cast: e.g. z:D
  //  ptr              : pointer to array to fill
  //  mxlen            : max len of ptr array to protect during read;
  //                     mxlen=0 --> do NOT store ptr (for text dump)
  //  vboseflag        : bit0-> print info to stdout if non-zero
  //                     bit1-> abort on missing variable
  // 
  // Output:
  //   varName_noCast 
  //
  // Default cast is double. Options are
  //   [VARNAME]/D or [VARNAME]:D  -> double
  //   [VARNAME]/F or [VARNAME]:F  -> float
  //   [VARNAME]/I or [VARNAME]:I  -> int
  //   [VARNAME]/C or [VARNAME]:C  -> char
  //
  // Sep 12 2016: allow up to 2 pointers for each column.

  int ivar, i, NVAR_TOT, NVAR_READ, VECFLAG, ISIZE, NPTR, MATCH, LEN ;
  int ICAST_STORE ;
  char varName[MXCHAR_VARNAME*2], *varTmp ;
  char fnam[] = "sntable_readprep_vardef1" ;

  // ---------------- BEGIN ---------------

  NVAR_TOT  = READTABLE_POINTERS.NVAR_TOT ;
  NVAR_READ = READTABLE_POINTERS.NVAR_READ ;


  if ( NVAR_TOT <= 0 ) {
    sprintf(MSGERR1,"Cannot find ' %s' because NVAR_TOT=%d",
	    varName_withCast, NVAR_TOT );
    MSGERR2[0] = 0 ;
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
  }

  // parse to get ICAST and 'varName' without cast-suffix.
  // Ignore VECFLAG and ISIZE args.
  // Note that ICAST is the cast of the array to be filled
  // on read-back, and it is NOT the cast of the stored column data.

  parse_TABLEVAR(varName_withCast,                    // (I)
		 varName,  &ICAST_STORE, &VECFLAG, &ISIZE);  // (O)

  LEN = strlen(varName);
  if ( LEN > MXCHAR_VARNAME ) {
    print_preAbort_banner(fnam);
    printf("\t varName = '%s' \n", varName);
    sprintf(MSGERR1,"len(varName) = %d exceeds bound of MXCHAR_VARNAME=%d .",
	    LEN, MXCHAR_VARNAME);
    sprintf(MSGERR2,"Try shorter name.");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  sprintf(varName_noCast, "%s", varName); // load output arg.

  // store max length of array to make check when reading later.
  READTABLE_POINTERS.MXLEN = mxlen ;


  // abort if this varname has already been stored by user
  // more than once

  for ( i=0; i < NVAR_READ; i++ ) {
    ivar = READTABLE_POINTERS.PTRINDEX[i] ;

    varTmp = READTABLE_POINTERS.VARNAME[ivar] ;
    NPTR   = READTABLE_POINTERS.NPTR[ivar];
    MATCH  = ( strcmp(varName,varTmp) == 0 );

    if ( MATCH && NPTR > 1 && ivar > 0 ) {
      sprintf(MSGERR1, "'%s' already stored (ivar=%d).", 
	      varName, ivar);
      sprintf(MSGERR2, "Cannot store the same variable %d times", 
	      NPTR+1);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }
  }

  ivar = -777 ;

  // search VARNAMES list
  for ( ivar=0; ivar < NVAR_TOT; ivar++ ) {

    if ( strcmp(varName,READTABLE_POINTERS.VARNAME[ivar]) != 0 ) 
      { continue ; }
    
    i = READTABLE_POINTERS.NVAR_READ ;
    READTABLE_POINTERS.PTRINDEX[i] = ivar ;
    READTABLE_POINTERS.NVAR_READ++ ;    
    READTABLE_POINTERS.ICAST_STORE[ivar]  = ICAST_STORE;
    READTABLE_POINTERS.NPTR[ivar]++ ;  // Sep 2016 
    NPTR = READTABLE_POINTERS.NPTR[ivar] ;

    if ( mxlen == 0 ) { goto  PTRSTORE_DONE ; }

    if ( ICAST_STORE == ICAST_D ) 
      { READTABLE_POINTERS.PTRVAL_D[NPTR-1][ivar] = (double*)ptr ; } 

    else if ( ICAST_STORE == ICAST_F )
      { READTABLE_POINTERS.PTRVAL_F[NPTR-1][ivar] = (float*)ptr ; }  

    else if ( ICAST_STORE == ICAST_I ) 
      { READTABLE_POINTERS.PTRVAL_I[NPTR-1][ivar] = (int*)ptr ; }

    else if ( ICAST_STORE == ICAST_S ) // short int 
      { READTABLE_POINTERS.PTRVAL_S[NPTR-1][ivar] = (short int*)ptr ; }

    else if ( ICAST_STORE == ICAST_L ) 
      { READTABLE_POINTERS.PTRVAL_L[NPTR-1][ivar] = (long long int*)ptr ; }

    else if ( ICAST_STORE == ICAST_C ) 
      { READTABLE_POINTERS.PTRVAL_C[NPTR-1][ivar] = (char**)ptr ; }

    else {
      sprintf(MSGERR1,"Unknown ICAST = %d", ICAST_STORE );
      sprintf(MSGERR2,"See ICAST_[D,F,I] parameters in sntools_output.h");
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
    }
    
  PTRSTORE_DONE:
    if ( (vboseflag & 1) > 0  ) {
      //      sprintf(ctmp,"");
      //      if ( NPTR > -1 ) { sprintf(ctmp, " <== fill %d arrays", NPTR); }
      printf("\t Prepare read-variable '%s' (ivar=%d, CAST=%c) \n", 
	     varName, ivar, CCAST_TABLEVAR[ICAST_STORE]  );
      fflush(stdout);
    }

    return ivar ;
      
  } // end of ivar


  if  ( (vboseflag & 2) > 0 ) {
    sprintf(MSGERR1,"Could not find required variable '%s' in table.", 
	    varName);
    sprintf(MSGERR2, "Check VARNAMES list (VBOSEFLAG=%d)", vboseflag );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    return (-1);
  }
  else 
    { return(-1) ; }


}  // end of  sntable_readprep_vardef1


// ===============================================================
void load_READTABLE_POINTER(int IROW, int IVAR, double DVAL, char *CVAL) {

  // Oct 26 2014
  // load user-defined array for
  // Inputs:
  //   IROW:  row or event number (starts at 0)
  //   IVAR:  sprase variable index
  //   DVAL:  value (double,float or int)
  //   CVAL:  char val (only if char)
  //

  char fnam[] = "load_READTABLE_POINTER" ;
  char *VARNAME ;
  int  IVAR_TOT, ICAST ;
  // ---------------- BEGIN ----------

  // translate sparse IVAR index into absolute index
  IVAR_TOT =  READTABLE_POINTERS.PTRINDEX[IVAR] ;
  ICAST    =  READTABLE_POINTERS.ICAST_STORE[IVAR_TOT] ;
  VARNAME  =  READTABLE_POINTERS.VARNAME[IVAR_TOT] ; // for err-msg only

  // avoid over-writing user-define array bound
  if ( IROW >= READTABLE_POINTERS.MXLEN ) {
    sprintf(MSGERR1, "IROW=%d exceeds user-defined array bound for", IROW) ;
    sprintf(MSGERR2, "VARNAME='%s'  (IVAR_TOT=%d  ICAST=%d)",
	    VARNAME, IVAR_TOT, ICAST);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    return ;
  }

  int nptr;
  for(nptr=0; nptr<READTABLE_POINTERS.NPTR[IVAR_TOT]; nptr++ ) {
    if ( ICAST == ICAST_D ) 
      { READTABLE_POINTERS.PTRVAL_D[nptr][IVAR_TOT][IROW] = (double)DVAL ; }

    else if ( ICAST == ICAST_L ) 
      { READTABLE_POINTERS.PTRVAL_L[nptr][IVAR_TOT][IROW] =
	  (long long int)DVAL;}
    
    else if ( ICAST == ICAST_F ) 
      { READTABLE_POINTERS.PTRVAL_F[nptr][IVAR_TOT][IROW] = (float)DVAL ; }
    
    else if ( ICAST == ICAST_I ) 
      { READTABLE_POINTERS.PTRVAL_I[nptr][IVAR_TOT][IROW] = (int)DVAL ; }

    else if ( ICAST == ICAST_S )  // short int
      { READTABLE_POINTERS.PTRVAL_S[nptr][IVAR_TOT][IROW] = (int)DVAL ; }
    
    else if ( ICAST == ICAST_C ) 
      { sprintf(READTABLE_POINTERS.PTRVAL_C[nptr][IVAR_TOT][IROW],"%s",CVAL);}
  }

} // end of load_READTABLE_POINTER


// ============================================
void load_DUMPLINE(int OPT, char *LINE, double DVAL) {

  // Oct 26 2014
  // update input *LINE with DVAL.
  // If DVAL is an int, write int format.
  // *LINE is intended for a dump to ascii file.
  // If OPT=1, DVAL is IFILTOBS and write char band before IFILTOBS

  // Mar 11 2019: use 'long long' instead of int.
  // Mar 04 2021: add OPT arg

  long long int LVAL = (long long int)DVAL;
  char STRVAL[40], BAND[2];
  int  IFILTOBS;
  // -------------- BEGIN ----------------
  if ( OPT == 1 ) {
    IFILTOBS = (int)DVAL ;
    sprintf(BAND,"%c", FILTERSTRING[IFILTOBS]);
    strcat(LINE," ");      strcat(LINE,BAND);
  }

  if ( (DVAL - (double)LVAL) == 0.0 ) 
    { sprintf(STRVAL," %lld", LVAL ); }
  else
    { sprintf(STRVAL," %.4f", DVAL ); }

  strcat(LINE,STRVAL);

  return ;

} // end of load_DUMPLINE

// ============================================
void load_DUMPLINE_STR(char *LINE, char *STRING) {
  strcat(LINE," "); strcat(LINE,STRING);
}

// ========================================
int  get_ICAST_READTBLE_POINTER(char *varName) {

  // Created Jan 2017
  // Return integet ICAST for input *varName
  // This ICAST is what is stored in the file,
  // not what is stored in memoery after readback.

  int ivar, ICAST ;
  char *tmpVar;
  // ----------- BEGIN ---------

  ICAST = -9 ;
  for( ivar=1; ivar <= READTABLE_POINTERS.NVAR_TOT; ivar++ ) {    
    tmpVar = READTABLE_POINTERS.VARNAME[ivar];  
    if ( strcmp(tmpVar,varName) == 0 ) 
      { ICAST = READTABLE_POINTERS.ICAST_READ[ivar];  }
  }

  return(ICAST);

} // end get_ICAST_READTBLE_POINTER

// ============================================
int SNTABLE_READ_EXEC(void) {

  // Oct 2014: 
  // execute table-read and fill arrays defined by
  // previous calls to SNTABLE_READPREP_VARDEF.
  // Function returns number of table rows read.

  int IFILETYPE = READTABLE_POINTERS.IFILETYPE ;
  int NROW ;
  char fnam[] = "SNTABLE_READ_EXEC" ;

  // --------------- BEGIN ---------------

  NROW = -777 ;

#ifdef USE_HBOOK
  if ( IFILETYPE == IFILETYPE_HBOOK ) {
    NROW = SNTABLE_READ_EXEC_HBOOK();
  }
#endif


#ifdef USE_ROOT
  if ( IFILETYPE == IFILETYPE_ROOT ) {
    NROW = SNTABLE_READ_EXEC_ROOT();
  }
#endif

#ifdef USE_TEXT
  if ( IFILETYPE == IFILETYPE_TEXT ) {
    NROW = SNTABLE_READ_EXEC_TEXT();
  }
#endif
  
  // sanity check
  if ( NROW == -777 ) {
    sprintf(MSGERR1,"Unknown file type -> cannot exec table-read for");
    sprintf(MSGERR2,"table = '%s'", READTABLE_POINTERS.TABLENAME );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);   
  }


  printf("\t ==> Finished reading %d of %d variables with %d entries\n",
         READTABLE_POINTERS.NVAR_READ, 
         READTABLE_POINTERS.NVAR_TOT,
	 NROW );
  fflush(stdout);

  READTABLE_POINTERS.NROW = NROW ; // store in global

  return NROW ;

} // end of  SNTABLE_READ_EXEC

// =====================================
void SNTABLE_LIST(char *FILENAME) {

  // list tables inside *FILENAME.

  char msg[] = "none." ;
  char fnam[] = "SNTABLE_LIST" ;
  int  FOUND_FILE = 0 ;

  TABLEFILE_INIT_VERIFY(fnam, msg) ;

#ifdef USE_HBOOK
  if ( ISFILE_HBOOK(FILENAME) )  { 
      SNTABLE_LIST_HBOOK(FILENAME) ; 
      FOUND_FILE = 1; 
  }
#endif

#ifdef USE_ROOT
  if ( ISFILE_ROOT(FILENAME) )  { 
    SNTABLE_LIST_ROOT(FILENAME); 
    FOUND_FILE = 1; 
  }
#endif

  // if we get here, abort.
  if( FOUND_FILE==0 ) { TABLEFILE_noFile_ABORT(fnam,FILENAME); }
  
} // end of SNTABLE_LIST


// ===================================================
void TABLEFILE_noFile_ABORT(char *FUNNAME, char *FILENAME) {
  sprintf(MSGERR1,"Unable to determine file type for TABLEFILE=");
  sprintf(MSGERR2,"%s", FILENAME);
  errmsg(SEV_FATAL, 0, FUNNAME, MSGERR1, MSGERR2);   
}

// =========================================
void TABLEFILE_notOpen_ABORT(char *FUNNAME, char *comment) {

  // Sep 2014
  // abort because no table files are open.
  char fnam[] = "TABLEFILE_notOpen_ABORT" ;

  print_preAbort_banner(fnam);
  printf("\t Check &SNLCINP namelist variables \n");
  printf("\t HFILE_OUT, ROOTFILE_OUT, TEXTFILE_PREFIX \n" );

  sprintf(MSGERR1,"%s", comment);
  sprintf(MSGERR2,"because no table-files are open.");
  errmsg(SEV_FATAL, 0, FUNNAME, MSGERR1, MSGERR2);     
  
} // end of TABLEFILE_notOpen_ABORT

// =========================================
void TABLEFILE_notCompiled_ABORT(char*FILENAME, char*FORMAT, char *ENV) {

  // Apr 2015
  // Leave Abort message because *FORMAT is not compiled in code,
  // and therefore *FILENAME cannot be opened.
  // Also leave message that *ENV is required for *FORMAT 
  // to be compiled.

  char fnam[] = "TABLEFILE_notCompiled_ABORT" ;
  
  print_preAbort_banner(fnam);
  printf("   Cannot open %s out-file \n\t '%s' \n", 
	 FORMAT, FILENAME);
  printf("   because code is not compiled with %s.\n", 
	 FORMAT);
  printf("   To compile with %s, make sure ENV $%s is set.\n",
	 FORMAT, ENV);

  sprintf(MSGERR1,"Cannot open %s file above", FORMAT);
  sprintf(MSGERR2,"because code is compiled without %s", FORMAT);
  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);

}  // end of TABLEFILE_notCompiled_ABORT

// ===================================================
void TABLEFILE_INIT_VERIFY(char *FUNNAME, char *MSG) {

  if ( CALLED_TABLEFILE_INIT != 740 ) {
    sprintf(MSGERR1,"Must call TABLEFILE_INIT ! ");
    sprintf(MSGERR2,"USER MSG: %s", MSG);
    errmsg(SEV_FATAL, 0, FUNNAME, MSGERR1, MSGERR2); 
  }
}

// ===========================================================
void SNTABLE_DUMP_VARNAMES(char *FILENAME, char *TABLENAME) {

  char fnam[] = "SNTABLE_DUMP_VARNAMES" ;
  char msg[80];
  int  FOUND_FILE = 0 ;

  // --------------- BEGIN --------------

  sprintf(msg,"Dump Table = '%s'", TABLENAME);
  TABLEFILE_INIT_VERIFY(fnam, msg) ;

#ifdef USE_HBOOK
  if ( ISFILE_HBOOK(FILENAME) ) {
    int NTID = atoi(TABLENAME);
    SNTABLE_DUMP_VARNAMES_HBOOK(FILENAME,NTID) ;
    FOUND_FILE = 1; 
  }
#endif

#ifdef USE_ROOT
  if ( ISFILE_ROOT(FILENAME) )  { 
    SNTABLE_DUMP_VARNAMES_ROOT(FILENAME,TABLENAME); 
    FOUND_FILE = 1; 
  }
#endif


  if ( FOUND_FILE == 0 )  { TABLEFILE_noFile_ABORT(fnam,FILENAME); }

} // end of SNTABLE_DUMP_VARNAMES


// ============================================================
int SNTABLE_DUMP_VALUES(char *FILENAME, char *TABLENAME, 
			int NVAR, char **VARLIST, int IVAR_NPT, 
			FILE *FP_OUTFILE, 
			char *LINEKEY_DUMP, char *SEPKEY_DUMP ) {

  // create text/fitres file with list if variables
  // from input VARLIST. Function returns number of
  // rows dumped.
  // 
  // Inputs:
  //  FILENAME     name of table file to read (hbook or root)
  //  TABLENAME    name of table to read
  //  NVAR:        number of variables to read from table
  //  **VARLIST    list of variable names in table
  //  FP_OUTFILE   pointer to ascii file to write to
  //                 (note that file must be opened by calling function)
  //  LINEKEY_DUMP: key at start of each row.
  // 
  // Dec 8 2014: change VBOSE from 1 to 3 to that it aborts on missing var.
  //
  // Jul 22 2017: if LINEKEY == "IGNORE:" then write out char BAND
  // Oct 31 2019: give better error message if NPTFIT is missing.

  int  NREAD = 0 ;
  char msg[80] ;
  char fnam[] = "SNTABLE_DUMP_VALUES" ;

  // --------------- BEGIN -------------

  sprintf(msg,"Dump Table = '%s'", TABLENAME);
  TABLEFILE_INIT_VERIFY(fnam, msg) ;


  int ivar, ICAST, NVAR_TOT, NC=0, IFILETYPE, MXLEN=0 ;
  int    VBOSE = 3 ; // 1-->print each var; 2--> abort in missing var
  double DDUMMY ;
  char   CDUMMY[80], *ptrVar, varName_NPT[20];
  char stringOpt[] = "read";

  
  IFILETYPE = TABLEFILE_OPEN(FILENAME, stringOpt); // open file to read
  NVAR_TOT  = SNTABLE_READPREP(IFILETYPE,TABLENAME); // read varList

  // Oct 2019: 
  //   Abort if OUTLIER flag is set but NPTFIT isn't defined in table.
  //   Code below aborts anyway if NPTFIT is missing, but generic
  //   'missing variable' error message is cryptic. Here the error 
  //   message explains what to do and how to set SNTABLE_LIST.
  
  if ( OUTLIER_INFO.USEFLAG ) {
    sprintf(varName_NPT, "%s", VARLIST[IVAR_NPT] ) ; // Nov 30 2020
    ivar = IVAR_READTABLE_POINTER(varName_NPT) ;
    if ( ivar < 0 ) {
      sprintf(MSGERR1,"Cannot extract obs/outiers "
	      "because %s is missing.", varName_NPT);
      sprintf(MSGERR2,"Try SNTABLE_LIST = 'FITRES+RESIDUALS' ");
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }
  }


  // store each variable to read and dump
  for(ivar=0; ivar < NVAR; ivar++ ) {   
    ptrVar = VARLIST[ivar] ;

    ICAST = ICAST_for_textVar(ptrVar);

    if ( ICAST == ICAST_C ) {
      SNTABLE_READPREP_VARDEF(ptrVar, CDUMMY, MXLEN, VBOSE); 
      NC++ ;
    }    
    else {
      SNTABLE_READPREP_VARDEF(ptrVar, &DDUMMY, MXLEN, VBOSE); 
    }
  }


  // store misc info.
  READTABLE_POINTERS.FP_DUMP   = FP_OUTFILE ;
  sprintf(READTABLE_POINTERS.LINEKEY_DUMP,"%s", LINEKEY_DUMP);
  sprintf(READTABLE_POINTERS.SEPKEY_DUMP, "%s", SEPKEY_DUMP);

  // do the read & write
  NREAD = SNTABLE_READ_EXEC();

  // close file that was read.
  TABLEFILE_CLOSE(FILENAME);

  return NREAD ;

} // end of SNTABLE_DUMP_VALUES


// =========================================
int  SNTABLE_DUMP_OUTLIERS(char *FILENAME, char *TABLENAME, 
			   int NVAR, char **VARLIST, int IVAR_NPT, 
			   float *OUTLIER_NSIGMA, 
			   FILE *FP_OUTFILE, char *LINEKEY, char *SEPKEY ) {

  // Created Aug 2014
  // create text/fitres file with list if variables
  // from input VARLIST.
  // Function returns number of outliers dumped.
  //
  // Oct 26 2014: use refactored system.
  // Jul 22 2017: add LINEKEY argument.
  // Feb    2018: check opton to dump everything
  // Nov 30 2020: pass IVAR_NPT to get column name NOBS or NPTFIT
  
  int  NDUMP, ivar, indx_store ;
  bool match ;
  char msg[80], *varName, simVar[80] ;
  char fnam[] = "SNTABLE_DUMP_OUTLIERS" ;

  // --------------- BEGIN -------------

  sprintf(msg,"Dump Outliers from Table = '%s'", TABLENAME);
  TABLEFILE_INIT_VERIFY(fnam, msg) ;

  // set global OUTLIER_INFO struct

  float Nsig0 = OUTLIER_NSIGMA[0] ;
  float Nsig1 = OUTLIER_NSIGMA[1] ;

  OUTLIER_INFO.USEFLAG = 1 ; // set global flag
  if ( Nsig0 == 0.0 && Nsig1 > 0.99E8 ) 
    { OUTLIER_INFO.USEFLAG = 2; } // flag to dump all OBS

  OUTLIER_INFO.CUTWIN_NSIGMA[0] = Nsig0 ;
  OUTLIER_INFO.CUTWIN_NSIGMA[1] = Nsig1 ;

  OUTLIER_INFO.CUTWIN_CHI2FLUX[0] = Nsig0*Nsig0 ;
  OUTLIER_INFO.CUTWIN_CHI2FLUX[1] = Nsig1*Nsig1 ;

  varName = OUTLIER_INFO.VARNAME[INDX_OUTLIER_NPTFIT] ;
  sprintf(varName, "%s", VARLIST[IVAR_NPT] ); // Nov 30 2020

  varName = OUTLIER_INFO.VARNAME[INDX_OUTLIER_IFILT] ;
  sprintf(varName, "%s", OUTLIER_VARNAME_IFILT) ;

  varName = OUTLIER_INFO.VARNAME[INDX_OUTLIER_CHI2] ;
  sprintf(varName, "%s", OUTLIER_VARNAME_CHI2) ;

  for(indx_store=0; indx_store < NVAR_OUTLIER_DECODE; indx_store++ ) 
    { OUTLIER_INFO.IVAR[indx_store] = -9 ; }

  for(ivar=0; ivar < NVAR; ivar++ ) {
    for(indx_store=0; indx_store < NVAR_OUTLIER_DECODE; indx_store++ ) {
      varName = OUTLIER_INFO.VARNAME[indx_store] ;
      match = false;

      if ( strcmp(VARLIST[ivar],varName) == 0 ) 
	{  OUTLIER_INFO.IVAR[indx_store] = ivar;  match = true;  }

      // try CHI2FLUX_SIM
      sprintf(simVar,"%s_SIM", OUTLIER_INFO.VARNAME[indx_store]);
      if ( strcmp(VARLIST[ivar],simVar) == 0 ) 
	{ OUTLIER_INFO.IVAR[indx_store] = ivar;  match = true; }

      /* xxx
      if ( match ) {
	printf(" xxx %s: varname = %s  indx_store=%d -> ivar=%d\n", 
	       fnam, varName, indx_store, ivar ); 
      }
      xxxx*/

    }
  } // end ivar loop


  NDUMP = 
    SNTABLE_DUMP_VALUES(FILENAME,TABLENAME, NVAR, VARLIST, IVAR_NPT,
			FP_OUTFILE, LINEKEY, SEPKEY );

  // print summary-stats to stdout
  if ( OUTLIER_INFO.USEFLAG == 1 ) { SNTABLE_SUMMARY_OUTLIERS(); }

  return NDUMP ;

} // end of SNTABLE_DUMP_OUTLIERS


// ================================
bool ISTABLEVAR_IFILT(char *VARNAME) {
  bool ISVAR = ( strcmp(VARNAME,OUTLIER_VARNAME_IFILT) == 0 )  ;
  return ISVAR;
} 

// ============================================================
void SNTABLE_SUMMARY_OUTLIERS(void) {

  // print outlier summary for each band, and for grand total.
  // print to stdout.

  float frac = 0.0 ;
  int   N1, N0, ISBAND, IFILT ;
  char  txt[20];

  printf("\n");
  
  for(IFILT=1; IFILT <= MXFILTINDX; IFILT++ ) {
    frac = 0.0 ;

    if ( IFILT < MXFILTINDX ) {
      N0 = OUTLIER_INFO.NEP_TOT[IFILT] ;
      N1 = OUTLIER_INFO.NEP_SELECT[IFILT] ; 
      sprintf(txt, "%c-band", FILTERSTRING[IFILT] );
      ISBAND = 1;
    }
    else {  
      // total over all filters
      N0 = OUTLIER_INFO.NEP_TOT[0] ;
      N1 = OUTLIER_INFO.NEP_SELECT[0] ; 
      sprintf(txt, "Total " );
      ISBAND = 0 ;
    }

    // bail for filter with no epochs, but always print
    // something for the total, even if there are no epochs.
    if ( N0 == 0 && ISBAND ) { continue ; }
    if ( N0 >  0 ) { frac = (float)N1 / (float)N0 ; }

    printf("  %s Outlier fraction : %5d/%6d = %.5f ",  txt, N1, N0, frac);
    if ( ISBAND )  { printf("  (IFILTOBS=%2d) ", IFILT ); }
    printf("\n");
    fflush(stdout);
  } // end IFILT loop

} // end of  SNTABLE_SUMMARY_OUTLIERS

void SNTABLE_AUTOSTORE_RESET(void) {
  NFILE_AUTOSTORE = 0 ;
}

void sntable_autostore_reset__(void)
{ SNTABLE_AUTOSTORE_RESET(); }

// =====================================================
int SNTABLE_AUTOSTORE_INIT(char *fileName, char *tableName, 
			   char *varList, int optMask) {

  // Created Oct 2014
  // One-call init to replace 3 functions:
  // SNTABLE_READPREP, SNTABLE_READPREP_VARDEF, SNTABLE_READ_EXEC.
  // *varList is a comma-separated list of variables.
  // Note that pointers are NOT passed, thus variables
  // are all auto-stored internally as double.
  //
  // Retreive values with 
  //   SNTABLE_AUTOSTORE_READ(CCID, varName, *istat, &VAL_D, VAL_C);
  //
  // This function is useful if you do NOT want to bother allocating 
  // your own memory, or if you need a fortran interface.
  // 
  // Inputs:
  //  *fileName  : name of file, any format (root, hbook or ascii)
  //  *tableName : name of table to read
  //  *varList   : comma-separated list of variables to read/store
  //               varList = 'ALL' --> read everything.
  //   optMask   : mask of options (was vboseflag)
  //      1=print, 2=abort if no varname matches, 4=append next file
  //
  // Output:
  //   Function returns number of table entries/rows.
  //
  // Nov 29 2016: avoid filling array when indx=-9 --> no variable.
  //              Free un-used storage memory.
  //
  // Dec 7 2016: allow ROW column name in addition to CID or CCID
  //
  // Jan 6 2017: (optMask & 4) -> append next file.
  //             Allow multiple files to be read & stored.
  //             Check for UNIQUE varName.
  //
  // Jan 11 2017: call TABLEFILE_CLOSE at end.
  //
  // May 11 2017: call IVAR_READTABLE_POINTER so that non-existent
  //              variables are NOT stored in autoStore arrays.
  //
  // July 11 2017: pass 'CID CCID' instead of 'CCID CID' so that
  //               last element (CCID) has priority if both are defined.
  //               Needed to work with ML_APPLY in NN pipeline.
  //
  // Oct 14 2020: 
  //   + fix ABORT feature if no variable name matches
  //   + use catVarList_with_comma util

  bool APPEND_FLAG, ABORT_FLAG;
  int  IFILETYPE, NF, ICAST, UNIQUE ;
  int  NVAR_USR, ivar, NROW, i, indx ;
  char *ptrtok, *tmpVar, varName_withCast[MXCHAR_VARNAME];
  char *varList_table, *varList_table_ptrtok;
  char varName[MXCHAR_VARNAME] ;
  char readOpt[] = "read";
  char blankFile[] = "" ;
  void *ptrStore;
  char fnam[] = "SNTABLE_AUTOSTORE_INIT" ;

  // -------------- BEGIN --------------

  // do global init if not already done.
  if ( CALLED_TABLEFILE_INIT != 740 ) { TABLEFILE_INIT();  }
  printf("   %s\n", fnam); fflush(stdout);

  // May 2022:
  // after reading autostore, reset NFILE for different autostore usage
  if ( NREAD_AUTOSTORE > 0 ) { NFILE_AUTOSTORE = NREAD_AUTOSTORE = 0; } 

  // check option to store multiple files
  ABORT_FLAG  = ( optMask & 2 ) ;
  APPEND_FLAG = ( optMask & 4 ) ; // append more variables
  if ( APPEND_FLAG || NFILE_AUTOSTORE==0 ) { NFILE_AUTOSTORE++ ; }
	

  varList_table        = (char*) malloc( MXCHAR_VARLIST * sizeof(char) );
  varList_table_ptrtok = (char*) malloc( MXCHAR_VARLIST * sizeof(char) );
  varList_table[0] = 0;

  NF = NFILE_AUTOSTORE-1;  // file index
  SNTABLE_AUTOSTORE[NF].NVAR = 0 ;
  SNTABLE_AUTOSTORE[NF].NROW = 0 ;
  SNTABLE_AUTOSTORE[NF].IFILETYPE = -9;

  // open file and return file type(root,hbook,text)
  IFILETYPE = TABLEFILE_OPEN(fileName,readOpt) ;
  SNTABLE_READPREP(IFILETYPE,tableName );

  // make comma-sep list of variables in table header;
  // skip CID by starting at ivar=1
  for( ivar = 1 ; ivar < READTABLE_POINTERS.NVAR_TOT; ivar++ ) {    
    tmpVar = READTABLE_POINTERS.VARNAME[ivar];  
    catVarList_with_comma(varList_table, tmpVar); 
  }

  if ( strcmp(varList,"ALL") == 0 ) 
    { sprintf(varList_table_ptrtok, "%s", varList_table); }
  else
    { sprintf(varList_table_ptrtok, "%s", varList); }

  // printf("\n xxx varList = '%s' \n\n", varList_local); fflush(stdout);

  // get list of variables
  ptrtok = strtok(varList_table_ptrtok,",");
  NVAR_USR = 0 ;
  while ( ptrtok != NULL ) {
    sprintf(varName, "%s",   ptrtok );
    ptrtok = strtok(NULL,"," );	

    if ( strcmp(varName,"CID" )  == 0 ) { continue ; } // avoid duplicate CCID
    if ( strcmp(varName,"CIDint")== 0 ) { continue ; } 
    if ( strcmp(varName,"CCID" ) == 0 ) { continue ; } // avoid duplicate CCID
    if ( strcmp(varName,"ROW"  ) == 0 ) { continue ; } // May 1 2017
    if ( strcmp(varName,"SNID" ) == 0 ) { continue ; } // Mar 13 2021
    if ( strcmp(varName,"GALID") == 0 ) { continue ; } // May 2021

    if ( IVAR_READTABLE_POINTER(varName) < 0 ) { continue ; }

    if ( EXIST_VARNAME_AUTOSTORE(varName) ) {
      sprintf(MSGERR1,"varName='%s' already defined in AUTOSTORE", varName);
      sprintf(MSGERR2,"Check previous AUTOSTORE declarations.");
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }

    sprintf(SNTABLE_AUTOSTORE[NF].VARNAME[NVAR_USR],"%s", varName);
    NVAR_USR++ ;
  }


  // abort if no variables are found 
  if ( ABORT_FLAG && NVAR_USR == 0 ) {
    print_preAbort_banner(fnam);

    printf("\n Requested variables to find:\n\t%s\n", varList);   
    printf("\n Available VARNAMES in file: \n\t%s\n", varList_table);

    sprintf(MSGERR1,"Found no varName matches in");
    sprintf(MSGERR2,"%s", fileName);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  // get number of rows needed to allocate memory
  NROW = SNTABLE_NEVT(blankFile,tableName);

  // - - - - - -
  // load scalars into global
  SNTABLE_AUTOSTORE[NF].NVAR      = NVAR_USR ; 
  SNTABLE_AUTOSTORE[NF].NROW      = NROW ;
  SNTABLE_AUTOSTORE[NF].IFILETYPE = IFILETYPE ; 

  printf("\t Malloc %d autoStore arrays for %d rows (IFILE=%d).\n", 
	 NVAR_USR, NROW, NF ); fflush(stdout);

  //  maloc CCID and pointers.
  SNTABLE_AUTOSTORE_malloc(0,NF,-9); 

  // init each variable with auto-generated memory
  // Tack on CID since user will fetch values based on CID.
  sprintf(varName_withCast,"CID:C  CCID:C  ROW:C  SNID:C  GALID:C");
  ivar = SNTABLE_READPREP_VARDEF(varName_withCast, 
				 SNTABLE_AUTOSTORE[NF].CCID, NROW, 1);

  if ( ivar < 0 ) {
    sprintf(MSGERR1,"Could not find required '%s'", varName_withCast);
    sprintf(MSGERR2,"Check %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  for(ivar=0; ivar < NVAR_USR; ivar++ ) {

    SNTABLE_AUTOSTORE[NF].EXIST[ivar] = 0 ;

    sprintf(varName,"%s", SNTABLE_AUTOSTORE[NF].VARNAME[ivar]) ;
    ICAST = get_ICAST_READTBLE_POINTER(varName); // cast in file

    // if varName is a duplicate from previous file, then append
    // varName -> varName_[NF+1] so that each stored varName is unique.
    // Allows combining files with same varNames.
    UNIQUE = UNIQUE_AUTOSTORE_VARNAME(NF,varName);
    if ( UNIQUE == 0 ) {
      sprintf(SNTABLE_AUTOSTORE[NF].VARNAME[ivar], "%s_%d", varName, NF+1);
    }

    if ( ICAST < 0 ) { continue ; }   

    if ( ICAST == ICAST_C )  {  
      sprintf(varName_withCast,"%s:C", varName ); 
      SNTABLE_AUTOSTORE_malloc(ICAST_C,NF,ivar); 
      ptrStore = SNTABLE_AUTOSTORE[NF].CVAL[ivar] ; 
    }
    else  {
      sprintf(varName_withCast,"%s:D", varName ); 
      SNTABLE_AUTOSTORE_malloc(ICAST_D,NF,ivar); 
      ptrStore = SNTABLE_AUTOSTORE[NF].DVAL[ivar] ;  
    }

    indx = SNTABLE_READPREP_VARDEF(varName_withCast, ptrStore, NROW, 1);

    if ( indx >= 0 ) { 
      SNTABLE_AUTOSTORE[NF].IVARMAP[indx] = ivar ; 
      SNTABLE_AUTOSTORE[NF].EXIST[ivar]   = 1; 
      SNTABLE_AUTOSTORE[NF].ICAST_READ[ivar]   = ICAST ; 
    }
    else {
      sprintf(MSGERR1,"Invalid indx=%d for %s", indx, varName_withCast);
      sprintf(MSGERR2,"returned from SNTABLE_READPREP_VARDEF" );
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }
    
  } // end ivar loop


  // read table and store values.
  NROW = SNTABLE_READ_EXEC();
  SNTABLE_AUTOSTORE[NF].NROW = NROW ;

  // ------------------------------------------
  // store string length for each CCID for faster lookup
  char *ptrCCID ;
  for(i=0; i < NROW; i++ )  { 
    ptrCCID = SNTABLE_AUTOSTORE[NF].CCID[i] ;
    SNTABLE_AUTOSTORE[NF].LENCCID[i] = strlen(ptrCCID);
  } 

  // init LASTREAD quantities
  LASTREAD_AUTOSTORE.IFILE = -9;
  LASTREAD_AUTOSTORE.IROW  = -9;
  sprintf(LASTREAD_AUTOSTORE.CCID,"XXX");

  free(varList_table); free(varList_table_ptrtok);

  return NROW ;

} // SNTABLE_AUTOSTORE_INIT


// fortran/mangled function
int sntable_autostore_init__(char *fileName, char *tableName, 
			     char *varList, int *optMask ) {
  return SNTABLE_AUTOSTORE_INIT(fileName,tableName,varList,*optMask);
}


// ===================================================
void  SNTABLE_AUTOSTORE_malloc(int OPT, int IFILE, int IVAR) {

  // Created Jan 2017
  // OPT =  0  --> mallic CCID, and create NVAR pointers
  // OPT = +1  --> malloc char for each row
  // OPT = +8  --> malloc double for each row
  // OPT = -1  --> free char mem 
  // OPT = -8  --> free double mem 
  //
  // IVAR used only to free memory for undefined variable.
  //

  int MEMC, MEMI, MEMD, NROW, NVAR_USR, i ;
  char fnam[] = "SNTABLE_AUTOSTORE_malloc" ;

  // -------------- BEGIN -------------

  MEMC = sizeof(char);
  MEMI = sizeof(int);
  MEMD = sizeof(double);
  NROW = SNTABLE_AUTOSTORE[IFILE].NROW ;
  NVAR_USR =  SNTABLE_AUTOSTORE[IFILE].NVAR ;

  if ( OPT == 0 ) {

    // setup CCID for each row
    SNTABLE_AUTOSTORE[IFILE].CCID = (char**)malloc( NROW*sizeof(char*) );
    for(i=0; i < NROW; i++ ) { 
      SNTABLE_AUTOSTORE[IFILE].CCID[i] =  (char*)malloc(MEMC*MXCHAR_CCID);  
    }

    SNTABLE_AUTOSTORE[IFILE].LENCCID = (int*)malloc( NROW*MEMI ) ;

    // set up char pointer for each variable (but not over rows)
    SNTABLE_AUTOSTORE[IFILE].CVAL  = (char***)malloc( NVAR_USR*sizeof(char**));

    // set up double pointer for each variable
    SNTABLE_AUTOSTORE[IFILE].DVAL = 
      (double**) malloc( NVAR_USR*sizeof(double*));

  }

  else if ( OPT == ICAST_C) {

    SNTABLE_AUTOSTORE[IFILE].CVAL[IVAR] = (char**)malloc(NROW*sizeof(char*)); 
    for(i=0; i < NROW; i++ )  { 
      SNTABLE_AUTOSTORE[IFILE].CVAL[IVAR][i] = 
	(char*)malloc(MEMC*MXCHAR_CCID);
    }
    
  }
  else if ( OPT == ICAST_D ) {
    SNTABLE_AUTOSTORE[IFILE].DVAL[IVAR] = (double*)malloc(NROW*MEMD); 

  }
  else if ( OPT == -ICAST_C ) {	  
    for(i=0; i<NROW; i++ ) { free(SNTABLE_AUTOSTORE[IFILE].CVAL[IVAR][i]); }
    free(SNTABLE_AUTOSTORE[IFILE].CVAL[IVAR]);
  }
  else if ( OPT == -ICAST_D ) {
    free( SNTABLE_AUTOSTORE[IFILE].DVAL[IVAR] );
  }

  if ( OPT < 0 ) {
    char *ptrVar = SNTABLE_AUTOSTORE[IFILE].VARNAME[IVAR] ;
    printf("\t %s: free mem for IFILE=%d IVAR=%d(%s) \n",
	   fnam, IFILE, IVAR, ptrVar ); fflush(stdout);
  }

  return ;

} // end SNTABLE_AUTOSTORE_malloc

// ============================================
int EXIST_VARNAME_AUTOSTORE(char *varName) {

  // Created Nov 29 2016
  // Returns 1 if varName exists.
  // 
  // Jan 6, 2017: check all files. 
  // Oct 27 2020: if varName == LIST, list all varNames and return 0
  // Dec 13 2021: refactor to use IVAR_VARNAME_AUTOSTORE
  //
  int ivar;
  // ------- BEGIN ---------

  ivar = IVAR_VARNAME_AUTOSTORE(varName);

  if ( ivar >=  0 ) 
    { return 1; }  // true
  else
    { return(0); } // false
  
} // end   int EXIST_VARNAME_AUTOSTORE

int exist_varname_autostore__(char *varName) {
  return EXIST_VARNAME_AUTOSTORE(varName);
}


// ============================================
int IVAR_VARNAME_AUTOSTORE(char *varName) {

  // Created Dec 13 2021
  // Returns ivar [0:NVAR-1] if varName exists; else return -9
  // 

  int ivar, ifile, ivar_tot, NVAR_USR ;
  char *varName_autostore;
  bool PRINT_LIST = ( strcmp(varName,"LIST") == 0 ) ;
  // ------- BEGIN ---------

  ifile = ivar_tot = 0;
  for(ifile=0; ifile < NFILE_AUTOSTORE; ifile++ ) {
    NVAR_USR = SNTABLE_AUTOSTORE[ifile].NVAR ;
    for(ivar=0; ivar < NVAR_USR; ivar++ ) {
      varName_autostore = SNTABLE_AUTOSTORE[ifile].VARNAME[ivar];
      if ( PRINT_LIST ) {
	printf("\t VARNAME[ifile=%d,ivar=%2.2d] = %s\n", 
	       ifile, ivar, varName_autostore); fflush(stdout);
      }
      if ( strcmp(varName_autostore,varName)==0 ) {
	if ( SNTABLE_AUTOSTORE[ifile].EXIST[ivar] ) { return(ivar_tot) ; }
      }
      ivar_tot++ ;
    }
  }

  // varName does NOT exist
  return(-9);
  
} // end   int IVAR_VARNAME_AUTOSTORE


// ===========================================
int IVAR_READTABLE_POINTER(char *varName) {

  // Created May 11 2017
  // For input varName, returns IVAR index in READTABLE_POINTERS.
  // If no match, returns -9.

  int ivar, IVAR = - 9;
  char *tmp_varName;
  //  char fnam[] = "IVAR_READTABLE_POINTER";

  // ------------- BEGIN ------------

  for( ivar = 1 ; ivar < READTABLE_POINTERS.NVAR_TOT; ivar++ ) {    
    tmp_varName = READTABLE_POINTERS.VARNAME[ivar];  
    if ( strcmp(varName,tmp_varName) == 0 ) { IVAR = ivar; }
  }

  return(IVAR);

} // end IVAR_READTABLE_POINTER

// ===========================================
int  UNIQUE_AUTOSTORE_VARNAME(int IFILE, char *VARNAME) {

  // Created Jan 2017
  // Returns true (=1) of this VARNAME is unique.
  // Returns false (=0) of VARNAME exists in previous ifile.

  int ifile, ivar, NVAR_USR ;
  char *tmpVarName ;
  //  char fnam[] =  "UNIQUE_AUTOSTORE_VARNAME" ;

  // ------------- BEGIN ------------

  // on first file, all varNames are unique by definition
  if ( IFILE == 0 ) { return(1); }

  for(ifile=0; ifile < IFILE; ifile++ ) {
    NVAR_USR = SNTABLE_AUTOSTORE[ifile].NVAR ;
    for(ivar=0; ivar < NVAR_USR; ivar++ ) {
      tmpVarName = SNTABLE_AUTOSTORE[ifile].VARNAME[ivar] ;
      if ( strcmp(tmpVarName,VARNAME) == 0 ) { return(0); }
    }
  }

  return(1);

} // end UNIQUE_AUTOSTORE_VARNAME

// =====================================
void SNTABLE_AUTOSTORE_READ(char *CCID, char *VARNAME, int *ISTAT,
			    double *DVAL, char *CVAL ) {

  // Created Oct 2013 by R.Kessler
  //
  // For input CCID and variable name *varName,
  // If *VARNAME is double, return *DVAL
  // If *VARNAME is char,   return *CVAL
  //
  // *ISTAT =  0  if CCID is found;
  // *ISTAT = -1  if CCID is NOT found 
  // *ISTAT = -2  if VARNAME is NOT found
  //
  //
  // Jan 6, 2017: 
  //   + check each autoStore file for varName.
  //   + refactor to return both double (DVAL) and char (CVAL), 
  //     but only one is set according to the cast of VARNAME.
  //
  // Oct 20 2020:
  //   + do NOT abort on missing varname; instead, return ISTAT = -2

  int IVAR_READ, IFILE_READ, ivar, i ;
  int NVAR_USR, NROW_TOT, ICAST, LENCCID, IROW ;
  char *tmpCCID, *tmpVar;
  
  char fnam[] = "SNTABLE_AUTOSTORE_READ" ;

  // ------------- BEGIN --------------

  *ISTAT = -1 ;       // default is that CCID is not found.
  IVAR_READ = IFILE_READ = IROW = -9 ;
  NREAD_AUTOSTORE++ ;

  // search file and variable indices
  for(i=0; i < NFILE_AUTOSTORE; i++ ) {
    NVAR_USR = SNTABLE_AUTOSTORE[i].NVAR ;
    for(ivar=0; ivar < NVAR_USR; ivar++ ) {
      tmpVar = SNTABLE_AUTOSTORE[i].VARNAME[ivar] ;
      if ( strcmp(tmpVar,VARNAME) == 0 ) 
	{ IVAR_READ = ivar ; IFILE_READ = i ;  goto FIND_CCID; }   
    }
  }


  if ( IVAR_READ < 0  || IFILE_READ < 0 ) { *ISTAT = -2; return ; }

 FIND_CCID:

  ICAST    = SNTABLE_AUTOSTORE[IFILE_READ].ICAST_READ[IVAR_READ] ;    
  NROW_TOT = SNTABLE_AUTOSTORE[IFILE_READ].NROW ; // total number of rows read
  LENCCID  = strlen(CCID);

  // if IFILE and CCID are the same as last time, 
  // skip slow check of all CCIDs

  bool IS_SAME_FILE = ( IFILE_READ == LASTREAD_AUTOSTORE.IFILE );
  bool IS_SAME_CCID = ( strcmp(CCID,LASTREAD_AUTOSTORE.CCID)==0);
  if (IS_SAME_FILE && IS_SAME_CCID ) 
    { IROW =  LASTREAD_AUTOSTORE.IROW ;  goto SET_OUTVAL; }

  // do slow loop over each row and do CCID string match each row.
  for(i=0; i < NROW_TOT; i++ ) {
    tmpCCID = SNTABLE_AUTOSTORE[IFILE_READ].CCID[i] ;
    if ( tmpCCID[0] != CCID[0] ) { continue ; } // quick check first char
    if ( strcmp(tmpCCID,CCID) == 0 )  { IROW = i ;  goto SET_OUTVAL ; }
  }

  if ( IROW < 0 ) { return ; } // could not find CCID

 SET_OUTVAL:
  *ISTAT = 0 ;
  if ( ICAST == ICAST_C ) {  
    sprintf(CVAL,"%s",SNTABLE_AUTOSTORE[IFILE_READ].CVAL[IVAR_READ][IROW]) ; 
  }
  else  { 
    // return double for D,F,I
    *DVAL = SNTABLE_AUTOSTORE[IFILE_READ].DVAL[IVAR_READ][IROW]; 
  }

  LASTREAD_AUTOSTORE.IFILE = IFILE_READ ;
  LASTREAD_AUTOSTORE.IROW  = IROW;
  sprintf(LASTREAD_AUTOSTORE.CCID,"%s", CCID );

  // if we get here, return 'not found' value.
  return ;

} // end of SNTABLE_AUTOSTORE_READ

// fortran/mangle function
void sntable_autostore_read__(char *CCID, char *varName, int *ISTAT, 
			      double *DVAL, char *CVAL ) {
  SNTABLE_AUTOSTORE_READ(CCID,varName,ISTAT,DVAL,CVAL);
}


void fetch_autostore_ccid(int ifile, int isn, char *ccid) {
  // Created Jan 4 2021
  // If autostore CID names are not known, pass indices here to
  // get ccid.
  sprintf(ccid, "%s", SNTABLE_AUTOSTORE[ifile].CCID[isn] ) ;

} // end fetch_autostore_ccid

void fetch_autostore_ccid__(int *ifile, int *isn, char *ccid) 
{ fetch_autostore_ccid(*ifile, *isn, ccid); }



// ========================================
int SNTABLE_NEVT(char *FILENAME, char *TABLENAME) {

  // Created April 2014
  // return number of rows in *TABLENAME (ntuple ID or tree).
  // If input FILENAME is non-null, open/read Nrow/close.
  // If input FILENAME == "", then assume that file is 
  // alread open and just read Nrow.
  // Note that file-type logic is 
  //  [already open (USE_TABLEFILE)] || [ or FILENAME has type ]
  //
  // If there is no file, return NROW=0 (do NOT abort).
  //
  // Oct 13 2014: add call to SNTABLE_NEVT_TEXT(FILENAME)
  //

  int  ID ;
  int  ISOPEN_DEJA, ISTYPE_HBOOK, ISTYPE_ROOT, ISTYPE_TEXT ;
  int  NEVT = 0 ;
  char msg[80] ;
  char fnam[] = "SNTABLE_NEVT";

  // -------------- BEGIN -------------

  sprintf(msg,"search NROW from table = %s", TABLENAME);
  TABLEFILE_INIT_VERIFY(fnam, msg) ;

  ISTYPE_HBOOK = ISTYPE_ROOT = ISTYPE_TEXT = 0 ;

#ifdef USE_HBOOK
  ISOPEN_DEJA = USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_HBOOK] ;
  ISTYPE_HBOOK = ISFILE_HBOOK(FILENAME) ;
  if ( ISOPEN_DEJA || ISTYPE_HBOOK ) {
    ID   = TABLEID_HBOOK(TABLENAME) ;
    NEVT = SNTABLE_NEVT_HBOOK(FILENAME,ID) ; 
    return(NEVT);
  }
#endif


#ifdef USE_ROOT
  ISOPEN_DEJA = USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_ROOT] ;
  ISTYPE_ROOT = ISFILE_ROOT(FILENAME) ;
  if ( ISOPEN_DEJA || ISTYPE_ROOT )  { 
    NEVT = SNTABLE_NEVT_ROOT(FILENAME,TABLENAME); 
    return(NEVT);
  }
#endif

#ifdef USE_TEXT
  ISOPEN_DEJA = USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_TEXT] ;
  ISTYPE_TEXT = ISFILE_TEXT(FILENAME); 
  if ( ISOPEN_DEJA || ISTYPE_TEXT ) {
    NEVT = SNTABLE_NEVT_TEXT(FILENAME); 
    return(NEVT);
  }
#endif

  return NEVT ;

} // end of SNTABLE_NEVT

int  sntable_nevt__(char *FILENAME, char *NAME) {
  return SNTABLE_NEVT(FILENAME,NAME);
}


// ============================================================
// ==================== SNHIST FUNCTIONS ======================
// ============================================================


void SNHIST_INIT(int NDIM, int ID, char *TITLE, 
		 int *NBIN, double *XMIN, double *XMAX ) {

  char fnam[] = "SNHIST_INIT" ;

  if ( NOPEN_TABLEFILE == 0 ) {
    char comment[100];
    sprintf(comment, "Cannot init HISTID=%d (%s)", ID, TITLE );
    TABLEFILE_notOpen_ABORT(fnam, comment );
  }
  
#ifdef USE_HBOOK
  if ( SNLCPAK_USE_HBOOK ) 
    { SNHIST_INIT_HBOOK( NDIM, ID, TITLE, NBIN, XMIN, XMAX ) ; }
#endif


#ifdef USE_ROOT
  if ( SNLCPAK_USE_ROOT ) 
    { SNHIST_INIT_ROOT( NDIM, ID, TITLE, NBIN, XMIN, XMAX ) ; }
#endif

} // end of SNHIST_INIT


void snhist_init__(int *NDIM, int *ID, char *TITLE, 
		   int *NBIN, double *XMIN, double *XMAX ) {
  SNHIST_INIT(*NDIM, *ID, TITLE, NBIN, XMIN, XMAX);
}

void SNHIST_FILL(int NDIM, int ID, double *VALUE, double WGT) {
  
  // ---------------- BEGIN ----------
  
#ifdef USE_HBOOK
  if ( SNLCPAK_USE_HBOOK ) 
    { SNHIST_FILL_HBOOK(NDIM,ID,VALUE,WGT); }
#endif

#ifdef USE_ROOT
  if ( SNLCPAK_USE_ROOT ) 
    { SNHIST_FILL_ROOT(NDIM,ID,VALUE,WGT); }
#endif

}  // end of SNHIST_FILL

void snhist_fill__(int *NDIM, int *ID, double *VALUE, double *WGT) {
  SNHIST_FILL(*NDIM, *ID, VALUE, *WGT);
}




// ====================
void SNHIST_RDBINS(int NDIM, int ID, char *CTIT, int *NBIN,
		   double *XMIN, double *XMAX) {

  // Jun 2013
  // Read & return histogram binning info and text-title (CTIT).  
  // NBIN, XMIN and XMAX are 1d or 2d arrays depending 
  // on dimension of histogram.
  //
  // Input: ID = integer histogram identifier
  // Outputs:
  //   CTIT = title
  //   NBIN = Number of bins for each dimension
  //   XMIN = min value for each dimension
  //   XMAX = max value for each dimension
  //
  
  char fnam[] = "SNHIST_RDBINS" ;
  int  ISOPEN_DEJA ;
  // ------------- BEGIN -----------

  if ( NDIM < 1 || NDIM > 2 ) {
    sprintf(MSGERR1,"Invalid NDIM=%d for ID=%d", NDIM, ID );
    sprintf(MSGERR2,"NDIM = 1 or 2" ); 
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

#ifdef USE_HBOOK
  ISOPEN_DEJA = USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_HBOOK] ;
  if ( ISOPEN_DEJA ) { 
    SNHIST_RDBINS_HBOOK(NDIM,ID,CTIT,NBIN,XMIN,XMAX); 
    return ;
  }
#endif

#ifdef USE_ROOT
  ISOPEN_DEJA = USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_ROOT] ;
  if ( ISOPEN_DEJA )  { 
    SNHIST_RDBINS_ROOT(NDIM,ID,CTIT,NBIN,XMIN,XMAX); 
    return ;
  }
#endif


  // if we get here, abort.
  sprintf(MSGERR1,"Could not find %dD histogram, HID=%d", NDIM, ID);
  sprintf(MSGERR2,"Check histograms with sntable_dump.pl" ); 
  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 

} // end of SNHIST_RDBINS

void snhist_rdbins__(int *NDIM, int *ID, char *CTIT, int *NBIN,
		     double *XMIN, double *XMAX) {
  SNHIST_RDBINS(*NDIM, *ID, CTIT, NBIN, XMIN, XMAX);
}



void SNHIST_RDCONT(int NDIM, int ID, int NRDBIN, double *CONTENTS) {

  // read histogram contents
  // Inputs: 
  //   NDIM = 1 or 2 ==> dimenstion of histogram
  //   ID   = histogram id
  //   NRDBIN = number of bins
  //
  // Ouptut:
  //   *CONTENTS = histogram contents
  //

  char fnam[] = "SNHIST_RDCONT" ;
  int  ISOPEN_DEJA ;
  // ------------- BEGIN -----------

  if ( NDIM < 1 || NDIM > 2 ) {
    sprintf(MSGERR1,"Invalid NDIM=%d for ID=%d", NDIM, ID );
    sprintf(MSGERR2,"NDIM = 1 or 2" ); 
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

#ifdef USE_HBOOK
  ISOPEN_DEJA = USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_HBOOK] ;
  if ( ISOPEN_DEJA ) { 
    SNHIST_RDCONT_HBOOK(NDIM,ID,NRDBIN,CONTENTS); 
    return ;
  }
#endif

#ifdef USE_ROOT
  ISOPEN_DEJA = USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_ROOT] ;
  if ( ISOPEN_DEJA ) { 
    SNHIST_RDCONT_ROOT(NDIM,ID,NRDBIN,CONTENTS); 
    return ;
  }
#endif

  // if we get here, abort.
  sprintf(MSGERR1,"Could not find %dD histogram, HID=%d", NDIM, ID);
  sprintf(MSGERR2,"Check histograms with sntable_dump.pl" ); 
  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 

} // end of SNHIST_RDCONT

void snhist_rdcont__(int *NDIM, int *ID, int *NRDBIN, double *CONTENTS) {
  SNHIST_RDCONT(*NDIM, *ID, *NRDBIN, CONTENTS) ;
}


// =======================================================
// ================== SUBDIR FUNCTIONS ===================
// =======================================================


// =======================================================
void MAKEDIR_OUTPUT( char *CCID, int CID ) {

  // create directory for this SN.
  // (I) CCID : character representation
  // (I) CID  : integer representation (needed for HBOOK)
  //
  // if CID >= 0 then create subdir SN[CCID]
  // if CID < 0  then CCID is just a generic dirName 
  //

  //  char fnam[] = "MAKEDIR_OUTPUT" ;

  // --------------- BEGIN --------------

  if ( NOPEN_TABLEFILE == 0 ) { return ; }

#ifdef USE_HBOOK
  if ( SNLCPAK_USE_HBOOK ) { MAKEDIR_HBOOK(CCID,CID); }
#endif

#ifdef USE_ROOT
  if ( SNLCPAK_USE_ROOT ) { MAKEDIR_ROOT(CCID,CID); }
#endif

#ifdef USE_TEXT
#endif

  fflush(stdout);

}  // end of MAKEDIR_OUTPUT

void makedir_output__(char *CCID, int *CID ) {
  MAKEDIR_OUTPUT(CCID, *CID);
}


// ===================================================
void CDTOPDIR_OUTPUT(void) {

  // cd top topdir if in OPEN-NEW mode; do nothing for READ mode.

  int  NEW ;
  //  char fnam[] = "CDTOPDIR_OUTPUT" ;

  if ( NOPEN_TABLEFILE == 0 ) { return ; }

#ifdef USE_HBOOK
  NEW = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_HBOOK] ;
  if ( NEW && SNLCPAK_USE_HBOOK ) { CDTOPDIR_HBOOK();  }
#endif

#ifdef USE_ROOT
  NEW = USE_TABLEFILE[OPENFLAG_NEW][IFILETYPE_ROOT] ;
  if ( NEW && SNLCPAK_USE_ROOT )  { CDTOPDIR_ROOT();  }
#endif
  
} // end of CDTOPDIR_OUTPUT

void cdtopdir_output__ (void) {
  CDTOPDIR_OUTPUT();
}


// ==========================================
//    SNLCPAK functions
// ==========================================

void SNLCPAK_INIT(char *SURVEY, char *VERSION_PHOT, char *VERSION_SNANA,
		  char *SURVEY_FILTERS, int SIMFLAG, 
		  int NFIT_PER_SN, char *TEXT_FORMAT) {

  // One-time global init called before opening output file.
  // NFIT = Number of fits per SN
  //
  // Sep 7 2014: add TEXT_FORMAT (used only if TEXT option is selected)
  // Dec 19 2019: pass SIMFLAG argument.

  char BANNER[100];
  char fnam[] = "SNLCPAK_INIT" ;

  // ---------------- BEGIN ------------

  sprintf(BANNER,"%s for SURVEY=%s  (%d fits per SN)",
	  fnam, SURVEY, NFIT_PER_SN);

  print_banner(BANNER) ;

  SNLCPAK_USE_HBOOK = 0 ; 
  SNLCPAK_USE_ROOT  = 0 ; 
  SNLCPAK_USE_TEXT  = 0 ; // Sep 7 2014
  SNLCPAK_USE_HDF5  = 0 ; 
  NCALL_SNLCPAK_FILL = 0 ;

  set_FILTERSTRING(FILTERSTRING);

  // store info in global
  sprintf(SNLCPAK_OUTPUT.SURVEY,             "%s", SURVEY         );
  sprintf(SNLCPAK_OUTPUT.VERSION_PHOTOMETRY, "%s", VERSION_PHOT   );
  sprintf(SNLCPAK_OUTPUT.VERSION_SNANA,      "%s", VERSION_SNANA  );
  sprintf(SNLCPAK_OUTPUT.SURVEY_FILTERS,     "%s", SURVEY_FILTERS );
  sprintf(SNLCPAK_OUTPUT.TEXT_FORMAT,        "%s", TEXT_FORMAT    );

  SNLCPAK_OUTPUT.NFILTDEF_SURVEY = strlen(SNLCPAK_OUTPUT.SURVEY_FILTERS) ;
  SNLCPAK_OUTPUT.SIMFLAG = SIMFLAG ;

  // store max possible fits per SN (can be fewer)
  SNLCPAK_OUTPUT.NFIT_PER_SN = NFIT_PER_SN ;
  if ( NFIT_PER_SN <= 0 ) {
    sprintf(MSGERR1,"Invalid NFIT_PER_SN = %d", NFIT_PER_SN);
    sprintf(MSGERR2,"NFIT_PER_SN must be > 0") ;
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  // clear arrays to allow for SNLCPAK_XXX calls.
  SNLCPAK_CLEAR_PLOT();     

  // ---------------------------------------------
  // set flag that init has been done.
  SNLCPAK_OUTPUT.INITDONE = SNLCPAK_INITDONE ;


}  // end of SNLCPAK_INIT

void snlcpak_init__(char *SURVEY, char *VER_PHOT, char *VER_SNANA,
		    char *SURVEY_FILTERS, int *SIMFLAG, 
		    int *NFIT_PER_SN, char *TEXTFMT) {
  SNLCPAK_INIT(SURVEY, VER_PHOT, VER_SNANA, SURVEY_FILTERS, 
	       *SIMFLAG, *NFIT_PER_SN, TEXTFMT );
}


// ======================================================
void SNLCPAK_NFIT_PER_SN(int NFIT_PER_SN) {
  // override value passed to SNLCPAK_INIT
  //  char fnam[] = "SNLCPAK_NFIT_PER_SN" ;
  printf("\t SNLCPAK_NFIT_PER_SN:  NFIT_PER_SN -> %d\n", NFIT_PER_SN);
  SNLCPAK_OUTPUT.NFIT_PER_SN = NFIT_PER_SN ;
}

void snlcpak_nfit_per_sn__(int *NFIT_PER_SN) {
  SNLCPAK_NFIT_PER_SN(*NFIT_PER_SN) ;
}

// ======================================================
void SNLCPAK_SURVEY(void) {

  // one-time global init after output file is opened.
  // Store global info such as name of survey, filters, etc ...

  //  char fnam[] = "SNLCPAK_SURVEY" ;

#ifdef USE_HBOOK
  if ( SNLCPAK_USE_HBOOK ) { SNLCPAK_SURVEY_HBOOK(); }
#endif

#ifdef USE_ROOT
  if ( SNLCPAK_USE_ROOT ) { SNLCPAK_SURVEY_ROOT(); }
#endif


} // end of SNLCPAK_SURVEY

void snlcpak_survey__(void) {
  SNLCPAK_SURVEY();
}


void SNLCPAK_CLEAR_PLOT(void) {

  // Called after filling each subdir;
  // clear SNLCPAK counters to prepare for next SUBDIR.
  
  int  iflag, ifilt, NOBS, i ;
  //  char fnam[] = "SNLCPAK_CLEAR_PLOT" ;

  // init counters
  SNLCPAK_OUTPUT.MXTEXT   = MXTEXT_SNLCPAK ;
  SNLCPAK_OUTPUT.NTEXT    = 0 ;
  SNLCPAK_OUTPUT.NOBS_MAX = 0 ;

  for ( iflag=0; iflag < MXFLAG_SNLCPAK_EPOCH; iflag++ )  {   
    NOBS = SNLCPAK_OUTPUT.NOBS[iflag] ;
    //    if ( NOBS == 0 ) { continue ; }

    SNLCPAK_OUTPUT.NOBS[iflag] = 0 ;   
    free(SNLCPAK_OUTPUT.MJD[iflag]);
    free(SNLCPAK_OUTPUT.TOBS[iflag]);
    free(SNLCPAK_OUTPUT.EPDATA[iflag]);
    free(SNLCPAK_OUTPUT.EPDATA_ERR[iflag]);
    free(SNLCPAK_OUTPUT.IFILTOBS[iflag]);
    free(SNLCPAK_OUTPUT.IFILT[iflag]);

    for ( ifilt=0; ifilt<MXFILTINDX; ifilt++ )  { 
      SNLCPAK_OUTPUT.NOBS_FILT[iflag][ifilt] = 0; 
    }
  }


  for ( iflag=0; iflag < MXFLAG_SNLCPAK_FILTER; iflag++ )  {   
    for ( ifilt=0; ifilt<MXFILTINDX; ifilt++ )  { 
      SNLCPAK_OUTPUT.FILTDATA[iflag][ifilt] = 0.0 ; 
      SNLCPAK_OUTPUT.FILTDATA_ERR[iflag][ifilt] = 0.0 ; 
    }
  }

  // set CCID to blank
  SNLCPAK_OUTPUT.CCID[0] = 0 ;

  // reset plot-filter info
  SNLCPAK_OUTPUT.NFILTDEF_PLOT = 0;
  SNLCPAK_OUTPUT.FILTLIST_PLOT[0] = 0 ;

  // rest dislay text
  for ( i=0; i < MXTEXT_SNLCPAK; i++ ) 
    { sprintf( SNLCPAK_OUTPUT.DISPLAYTEXT[i],"NULL" ) ; }

} // end of SNLCPAK_CLEAR_PLOT

void snlcpak_clear_plot__(void) {
  SNLCPAK_CLEAR_PLOT();
}

void SNLCPAK_CLEAR_SN(void) {
  // May 2014: !!! OBSOLETE !!!
  SNLCPAK_OUTPUT.NLCPAK      = 0 ; // incremental counter 
  SNLCPAK_OUTPUT.NLCPAK_TOT  = 1 ; // expected total number of fits/subdirs
  printf(" xxx SNLCPAK_CLEAR_SN: NLCPAK_TOT -> 1 \n" ); fflush(stdout);
}


void SNLCPAK_DISPLAYTEXT(char *CCID, char *DISPLAYTEXT) {

  // Store text strings to display on plot.
  // NOTHING to do with text output in TEXTFILE_OPEN.

  int N;
  char fnam[] = "SNLCPAK_DISPLAYTEXT" ;

  // ------------- BEGIN -------------

  SNLCPAK_CHECK(CCID,fnam);

  SNLCPAK_OUTPUT.NTEXT++ ;
  N = SNLCPAK_OUTPUT.NTEXT ;

  if ( N >= MXTEXT_SNLCPAK ) {
    sprintf(MSGERR1,"NTEXT exceeds bound of MXTEXT_SNLCPAK=%d",
	    MXTEXT_SNLCPAK);
    sprintf(MSGERR2,"CCID=%s, TEXT='%s'", CCID, DISPLAYTEXT);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  sprintf( SNLCPAK_OUTPUT.DISPLAYTEXT[N],"%s", DISPLAYTEXT) ;

} // end of SNLCPAK_TEXT

void snlcpak_displaytext__(char *CCID, char *DISPLAYTEXT) {
  SNLCPAK_DISPLAYTEXT(CCID,DISPLAYTEXT);
}

void SNLCPAK_NFIT(int NFIT) {
  // store number of light curve fits (NFIT) for this SN;
  // note that NFIT <= NFIT_PER_SN
  //  char fnam[] = "SNLCPAK_NFIT" ;
  SNLCPAK_OUTPUT.NLCPAK_TOT = NFIT;
  SNLCPAK_OUTPUT.NLCPAK     = 0 ; // reset counter (5/01/2014)
} 
void snlcpak_nfit__(int *NFIT) {
  SNLCPAK_NFIT(*NFIT);
}

void SNLCPAK_DATA(char *CCID, int NOBS, double *MJD, double *TOBS, 
		  double *DATA, double *DATA_ERR, int *IFILTOBS, int FLAG) {

  // store NOBS data values (epochs) passed from main program.
  // FLAG indicates data type:  DATA, REJECT-FLAG, FIT-CHI2, FIT-FUNCTION.
  // The first 3  must have NOBS and TOBS corresponding to the data;
  // the last one (FIT-FUNCTION) has different NOBS and TOBS, usually
  // in very fine TOBS-bins to show a smooth best-fit model function.
  //
  // May 20 2016: new FLAG for simulated flux

  int i, NFILT, IFILT, ifiltobs, flag, MEMD, MEMI ;
  char cfilt[2], comment[80];
  char fnam[] = "SNLCPAK_DATA";

  // -------- BEGIN ---------

  if ( NOPEN_TABLEFILE == 0 ) {
    sprintf(comment, "Cannot store light curve");
    TABLEFILE_notOpen_ABORT(fnam, comment );
  }

  // increment plot-counter for FLUX-DATA since this is the
  // only required FLAG
  if ( FLAG  == SNLCPAK_EPFLAG_FLUXDATA )  
    { SNLCPAK_OUTPUT.NLCPAK++ ; }

  sprintf(comment,"%s(FLAG=%d)", fnam, FLAG);
  SNLCPAK_CHECK(CCID,comment);

  if ( NOBS <= 0  ) {
    sprintf(MSGERR1,"Invalid NOBS = %d for CCID=%s and FLAG=%d . ", 
	    NOBS, CCID, FLAG);
    sprintf(MSGERR2,"NOBS must be > 0. " );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }


  // keep track of max NOBS among all FLAG values
  if ( NOBS >  SNLCPAK_OUTPUT.NOBS_MAX ) 
    {  SNLCPAK_OUTPUT.NOBS_MAX = NOBS ; }


  if ( FLAG < MXFLAG_SNLCPAK_EPOCH ) {

    // allocate memory
    MEMD = (NOBS+10) * sizeof(double);
    MEMI = (NOBS+10) * sizeof(int) ;  
    SNLCPAK_OUTPUT.TOBS[FLAG]            = (double*)malloc(MEMD) ;
    SNLCPAK_OUTPUT.MJD[FLAG]             = (double*)malloc(MEMD) ;
    SNLCPAK_OUTPUT.EPDATA[FLAG]          = (double*)malloc(MEMD) ;
    SNLCPAK_OUTPUT.EPDATA_ERR[FLAG]      = (double*)malloc(MEMD) ;
    SNLCPAK_OUTPUT.IFILTOBS[FLAG]        = (int*)malloc(MEMI) ;
    SNLCPAK_OUTPUT.IFILT[FLAG]           = (int*)malloc(MEMI) ;

    SNLCPAK_OUTPUT.NOBS[FLAG] = NOBS ;
    for ( i=0; i < NOBS; i++ ) {
      SNLCPAK_OUTPUT.MJD[FLAG][i]         = MJD[i] ;
      SNLCPAK_OUTPUT.TOBS[FLAG][i]        = TOBS[i] ;
      SNLCPAK_OUTPUT.EPDATA[FLAG][i]      = DATA[i] ;
      SNLCPAK_OUTPUT.EPDATA_ERR[FLAG][i]  = DATA_ERR[i] ;
      SNLCPAK_OUTPUT.IFILTOBS[FLAG][i]    = IFILTOBS[i] ;
    }
  }
  else {
    NFILT = NOBS ;  // one value per filter
    for (i=0; i < NFILT; i++ ) {
      ifiltobs = IFILTOBS[i] ;
      sprintf(cfilt,"%c", FILTERSTRING[ifiltobs] );
      IFILT = strcspn( SNLCPAK_OUTPUT.SURVEY_FILTERS, cfilt) ;
      flag  = FLAG - 100; // storage index

      if ( flag < 0 || flag >= MXFLAG_SNLCPAK_FILTER ) {
	sprintf(MSGERR1,"Invalid flag=%d for FLAG=%d", flag, FLAG);
	sprintf(MSGERR2,"IFILT=%d  IFILTOBS=%d(%s)", IFILT, ifiltobs, cfilt);
	errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
      }

      SNLCPAK_OUTPUT.FILTDATA[flag][IFILT]      = DATA[i] ; 
      SNLCPAK_OUTPUT.FILTDATA_ERR[flag][IFILT]  = DATA_ERR[i] ;

    }    // i loop
  }


}  // end of SNLCPAK_DATA

void snlcpak_data__(char *CCID, int *NOBS, double *MJD, double *TOBS, 
		    double *DATA, double *DATA_ERR, int *IFILTOBS,int *FLAG) {
  SNLCPAK_DATA(CCID, *NOBS, MJD, TOBS, DATA, DATA_ERR, IFILTOBS, *FLAG);
}

// ========================================
void SNLCPAK_FILL(char *CCID) {
 
  char fnam[] = "SNLCPAK_FILL" ;

  // ----------- BEGIN ------------

  SNLCPAK_CHECK(CCID,fnam);

  // do some global stuff needed for any output format.
  SNLCPAK_FILL_PREP();
  NCALL_SNLCPAK_FILL++ ;  

#ifdef USE_HBOOK
  if ( SNLCPAK_USE_HBOOK  ) { SNLCPAK_FILL_HBOOK(); }
#endif


#ifdef USE_ROOT
  if ( SNLCPAK_USE_ROOT  ) { SNLCPAK_FILL_ROOT(); }
#endif


#ifdef USE_TEXT
  if ( SNLCPAK_USE_TEXT  ) { SNLCPAK_FILL_TEXT(); }
#endif

 
  SNLCPAK_CLEAR_PLOT();

}  // end  SNLCPAK_FILL

void snlcpak_fill__(char *CCID) {
  SNLCPAK_FILL(CCID);
}

// ====================================
void SNLCPAK_FILL_PREP() {

  // compute quantities useful for for the _FILL_ routines.

  int  IFILT, FLAG, obs, NFILT, USE, NOBS, IFILTOBS, IFMIN, IFMAX ;
  char cfilt[2], *SURVEY_FILTERS ;
  //  char fnam[] = "SNLCPAK_FILL_PREP" ;

  // ------------- BEGIN --------------

  for ( FLAG=0; FLAG < MXFLAG_SNLCPAK_EPOCH; FLAG++ ) {
    for ( IFILT=0; IFILT < MXFILTINDX; IFILT++ )  { 
      SNLCPAK_OUTPUT.USEFILT[FLAG][IFILT] = 0 ; 
    }
  }


  // use required data array to check for which filters to plot.
  // Filter order is determined by SURVEY_FILTERS.

  SURVEY_FILTERS = SNLCPAK_OUTPUT.SURVEY_FILTERS ;

  for ( FLAG=0; FLAG < MXFLAG_SNLCPAK_EPOCH; FLAG++ ) {
    NOBS      = SNLCPAK_OUTPUT.NOBS[FLAG] ;
    for(obs=0; obs<NOBS; obs++ ) {   

      IFILTOBS = SNLCPAK_OUTPUT.IFILTOBS[FLAG][obs] ;
      sprintf(cfilt,"%c", FILTERSTRING[IFILTOBS] );
      IFILT = strcspn( SURVEY_FILTERS, cfilt);
      SNLCPAK_OUTPUT.USEFILT[FLAG][IFILT] = 1 ;  
      
      SNLCPAK_OUTPUT.IFILT[FLAG][obs] = IFILT ;  // sparse IFILT for this obs.
      SNLCPAK_OUTPUT.NOBS_FILT[FLAG][IFILT]++ ; //  NOBS for this filter

    }  // obs loop
  }  // FLAG loop


  // count how many defined filters for DATA, and set FILTLIST_PLOT
  NFILT  = 0 ;
  IFMIN  = IFMAX = -9 ;
  FLAG   = SNLCPAK_EPFLAG_FLUXDATA ;
  for ( IFILT=0; IFILT < MXFILTINDX; IFILT++ ) {
    USE = SNLCPAK_OUTPUT.USEFILT[FLAG][IFILT] ;
    sprintf(cfilt,"%c", SURVEY_FILTERS[IFILT] );
    if ( USE ) { 
      NFILT++ ; 
      sprintf(SNLCPAK_OUTPUT.FILTLIST_PLOT,"%s%s",
	      SNLCPAK_OUTPUT.FILTLIST_PLOT, cfilt) ;

      if ( IFMIN < 0 ) { IFMIN = IFILT ; }
      IFMAX = IFILT ;
    }
  }

  SNLCPAK_OUTPUT.NFILTDEF_PLOT = NFILT ;
  SNLCPAK_OUTPUT.IFILT_MIN     = IFMIN ;
  SNLCPAK_OUTPUT.IFILT_MAX     = IFMAX ;

} // end of SNLCPAK_FILL_PREP



// ===============================================
char *replace_str(char *st, const char *orig, const char *repl) {
  static char buffer[4096];
  char *ch;
  if (!(ch = strstr(st, orig)))
    return st;
  strncpy(buffer, st, ch-st);  
  buffer[ch-st] = 0;
  sprintf(buffer+(ch-st), "%s%s", repl, ch+strlen(orig));
  return buffer;
}


// ===========================================
void SNLCPAK_CHECK(char *CCID, char *comment) {

  // either store CCID in global array,
  // or check the CCID matches what is in global array.
  // This CHECK prevents switching CIDs before calling SNLCPAK_FILL.

  char *ptrCCID ;
  char fnam[] = "SNLCPAK_CHECK" ;

  // ------------ BEGIN -----------

  // make sure that SNLCPAK_INIT was called
  if ( SNLCPAK_OUTPUT.INITDONE != SNLCPAK_INITDONE ) {
    sprintf(MSGERR1,"%s called before SNLCPAK_INIT.", comment);
    sprintf(MSGERR2,"Must call SNLCPAK_INIT first.");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }


  ptrCCID = SNLCPAK_OUTPUT.CCID ;

  // store CCID or make sure that CCID is correct
  if ( strlen(ptrCCID) == 0 ) {
    // define CCID
    sprintf(ptrCCID,"%s", CCID);    
  }
  else if ( strcmp(CCID,ptrCCID) != 0 ) {
    // abort on error
    sprintf(MSGERR1,"Invalid CCID = '%s' passed to %s", CCID, comment);
    sprintf(MSGERR2,"Expected CCID = '%s'", ptrCCID) ;
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

} // end of SNLCPAK_CHECK



// ==========================================
//    SPECPAK functions (April 2019)
//    Only for TEXT output to feed external plotter. Ignore ROOT,HBOOK.
// ==========================================


void SPECPAK_INIT(char *SURVEY, char *VERSION_PHOT, char *TEXT_FORMAT) {

  // April 2019
  // One-time global init called before opening output file.


  char BANNER[100];
  char fnam[] = "SPECPAK_INIT" ;

  // ---------------- BEGIN ------------

  sprintf(BANNER,"%s for SURVEY=%s ", fnam, SURVEY);

  print_banner(BANNER) ;

  SPECPAK_USE_HBOOK = 0 ; 
  SPECPAK_USE_ROOT  = 0 ; 
  SPECPAK_USE_TEXT  = 0 ; // Sep 7 2014
  SPECPAK_USE_MARZ  = 0 ;
  NCALL_SPECPAK_FILL = 0 ;

  // store info in global
  sprintf(SPECPAK_OUTPUT.SURVEY,             "%s", SURVEY         );
  sprintf(SPECPAK_OUTPUT.VERSION_PHOTOMETRY, "%s", VERSION_PHOT   );
  sprintf(SPECPAK_OUTPUT.TEXT_FORMAT,        "%s", TEXT_FORMAT    );

  // clear arrays
  SPECPAK_CLEAR_PLOT();     

  // ---------------------------------------------
  // set flag that init has been done.
  SPECPAK_OUTPUT.INITDONE = SPECPAK_INITDONE ;


}  // end of SPECPAK_INIT

void specpak_init__(char *SURVEY, char *VER_PHOT, char *TEXTFMT) {
  SPECPAK_INIT(SURVEY, VER_PHOT, TEXTFMT );
}


// =========================================
void SPECPAK_CLEAR_PLOT(void) {

  // May 13 2019:
  //  + fix bug: NCALL_SPECPAK_FILL>10 --> >0 so that
  //    allocated memory is freed. Fixes valgrind errors.

  //  char fnam[] = "SPECPAK_CLEAR_PLOT" ;

  //  printf(" xxx %s  NCALL=%d \n", fnam, NCALL_SPECPAK_FILL ) ;

  SPECPAK_OUTPUT.NSPEC = 0 ;
  SPECPAK_OUTPUT.NLAMBIN_TOT = 0 ;

  if ( NCALL_SPECPAK_FILL > 0 ) {
    free(SPECPAK_OUTPUT.ID    ) ;
    free(SPECPAK_OUTPUT.LAMMIN) ;
    free(SPECPAK_OUTPUT.LAMMAX) ;
    free(SPECPAK_OUTPUT.FLAM) ;
    free(SPECPAK_OUTPUT.FLAMERR) ;
  }


} // end SPECPAK_CLEAR_PLOT

void specpak_clear_plot__(void) {  SPECPAK_CLEAR_PLOT();  }


// =======================================================
void SPECPAK_DATA(char *CCID, int ID, double MJD, double Tobs, double Texpose,
		  int NLAMBIN, double *LAMMIN, double *LAMMAX, 
		  double *FLAM, double *FLAMERR)
{

  // pack/store multiple spectra for one event.
  // Note that MJD = Tobs = -9 for host galaxy.

  int MEMD, MEMI, i, ii, NSPEC, NLAMTOT, ILAM_START ;
  char fnam[] = "SPECPAK_DATA" ;

  // ------------ BEGIN ------------

  if ( NLAMBIN <= 0  ) {
    sprintf(MSGERR1,"Invalid NLAMBIN = %d for CCID=%s and ID=%d . ", 
	    NLAMBIN, CCID, ID);
    sprintf(MSGERR2,"NLAMBIN must be > 0. " );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }

  NSPEC      = SPECPAK_OUTPUT.NSPEC;
  ILAM_START = SPECPAK_OUTPUT.NLAMBIN_TOT;
  SPECPAK_OUTPUT.NLAMBIN_TOT += NLAMBIN ;
  NLAMTOT    = SPECPAK_OUTPUT.NLAMBIN_TOT ;

  sprintf(SPECPAK_OUTPUT.CCID, "%s", CCID );
  SPECPAK_OUTPUT.ID_LIST[NSPEC]       = ID ;
  SPECPAK_OUTPUT.NLAMBIN_LIST[NSPEC]  = NLAMBIN ;
  SPECPAK_OUTPUT.MJD_LIST[NSPEC]      = MJD;
  SPECPAK_OUTPUT.TOBS_LIST[NSPEC]     = Tobs ;
  SPECPAK_OUTPUT.TEXPOSE_LIST[NSPEC]  = Texpose ;

  MEMD    = NLAMTOT * sizeof(double);
  MEMI    = NLAMTOT * sizeof(int) ;  

  // malloc wave-dependent storage
  if ( NSPEC == 0 ) {
    SPECPAK_OUTPUT.ID      = (int   *)malloc(MEMI);
    SPECPAK_OUTPUT.LAMMIN  = (double*)malloc(MEMD); // array of lammin in each bin
    SPECPAK_OUTPUT.LAMMAX  = (double*)malloc(MEMD); // array of lammax in each bin
    SPECPAK_OUTPUT.FLAM    = (double*)malloc(MEMD); // array of flam in each bin
    SPECPAK_OUTPUT.FLAMERR = (double*)malloc(MEMD); 
  }
  else {
    // realloc for more wave bins
    SPECPAK_OUTPUT.ID      = (int   *)realloc(SPECPAK_OUTPUT.ID,     MEMI );
    SPECPAK_OUTPUT.LAMMIN  = (double*)realloc(SPECPAK_OUTPUT.LAMMIN, MEMD );
    SPECPAK_OUTPUT.LAMMAX  = (double*)realloc(SPECPAK_OUTPUT.LAMMAX, MEMD );
    SPECPAK_OUTPUT.FLAM    = (double*)realloc(SPECPAK_OUTPUT.FLAM,   MEMD );
    SPECPAK_OUTPUT.FLAMERR = (double*)realloc(SPECPAK_OUTPUT.FLAMERR,MEMD );
  }

  // load SPECPAK_OUTPUT
  for(i=0; i < NLAMBIN; i++ ) {
    ii = i + ILAM_START ;
    SPECPAK_OUTPUT.ID[ii]      = ID ;
    SPECPAK_OUTPUT.LAMMIN[ii]  = LAMMIN[i];
    SPECPAK_OUTPUT.LAMMAX[ii]  = LAMMAX[i];
    SPECPAK_OUTPUT.FLAM[ii]    = FLAM[i];
    SPECPAK_OUTPUT.FLAMERR[ii] = FLAMERR[i];
  }

  SPECPAK_OUTPUT.NSPEC++ ;
  if ( SPECPAK_OUTPUT.NSPEC >= MXSPEC_SPECPAK ) {
    sprintf(MSGERR1, "NSPEC exceeds MXSPEC_SPECPAK=%d bound.", 
	    MXSPEC_SPECPAK);
    sprintf(MSGERR2, "Reduce NSPEC or increase MXSPEC_SPECPAK"); 
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  return ;

} // end SPECPAK_DATA


void specpak_data__(char *CCID, int *ID, double *MJD,double *Tobs,
		    double *Texpose,int *NLAMBIN,double *LAMMIN,double *LAMMAX,
		    double *FLAM, double *FLAMERR) {
  SPECPAK_DATA(CCID, *ID, *MJD, *Tobs, *Texpose,
	       *NLAMBIN, LAMMIN, LAMMAX, FLAM, FLAMERR);
}


// ========================================
void SPECPAK_FILL(char *CCID) {
 
  //  char fnam[] = "SPECPAK_FILL" ;

  // ----------- BEGIN ------------

  NCALL_SPECPAK_FILL++ ; 

  // do nothing for HBOOK or ROOT ... only for TEXT.

#ifdef USE_TEXT
  if ( SPECPAK_USE_TEXT  ) { SPECPAK_FILL_TEXT(); }
#endif

#ifdef USE_MARZ
  if ( SPECPAK_USE_MARZ  ) { SPECPAK_FILL_MARZ(); }
#endif


  SPECPAK_CLEAR_PLOT(); 

}  // end  SPECPAK_FILL

void specpak_fill__(char *CCID) {  SPECPAK_FILL(CCID); }



/*  =======================================
      ISFILE_XXX function

   These functions are usable regardless of whether hbook,root
   are defined, so at least we can figure out the file type.

 =========================================== */


// ==============================
int ISFILE_HBOOK(char *fileName) {

  // returns true if suffix corresponds to hbook.

#define NSUFFIX_HBOOK 6
  int   isuf ;
  char  SUFFIX_HBOOK[NSUFFIX_HBOOK][8] = 
    { 
      ".his" ,  ".HIS", 
      ".tup",   ".TUP", 
      ".hbook", ".HBOOK"  
    } ;
  
  for ( isuf=0; isuf<NSUFFIX_HBOOK; isuf++ ) {
    if ( strstr(fileName, SUFFIX_HBOOK[isuf] ) != NULL )  { return 1 ; }
  }
  
  return 0;

} // end of ISFILE_HBOOK

// ==============================
int ISFILE_ROOT(char *fileName) {

  // returns true if suffix corresponds to a root file.

#define NSUFFIX_ROOT 2
  int   isuf ;
  char  SUFFIX_ROOT[NSUFFIX_ROOT][8] = 
    { ".root" , ".ROOT"  } ;

  for ( isuf=0; isuf<NSUFFIX_ROOT; isuf++ ) {
    if ( strstr(fileName, SUFFIX_ROOT[isuf] ) != NULL )  { return 1 ; }
  }
  
  return 0;

} // end of ISFILE_ROOT

// ==============================              
int ISFILE_TEXT(char *fileName) {

  // Created Oct 2014
  // returns 1 (true) if suffix corresponds to ascii/text file,
  // or if VARNAMES key is found in header.
  // Returs 0 (false) otherwise.
  //
  // Note that the suffix method must be used for files
  // that don't yet exist; i.e., to be written.
  // The header-key method works only for reading.
  //
  // May 04 2020: return false for fits or FITS file extension.
  // Jan 22 2021: check HOSTLIB extension
  // Oct 04 2021: add M0DIF

#define NSUFFIX_TEXT 20
  int   isuf ;
  char  SUFFIX_TEXT_LIST[NSUFFIX_TEXT][10] = 
    { 
      ".text" ,   ".TEXT",
      ".txt" ,    ".TXT",
      ".fitres",  ".FITRES",
      ".m0dif",   ".M0DIF",
      ".snana",   ".SNANA",    // added Apr 17 2021 (for merged SNANA table)
      ".dat",     ".DAT",
      ".out",     ".OUT",      // added Feb 7 2017
      ".table",   ".TABLE",    // idem
      ".hostlib", ".HOSTLIB",  // added Jan 22 2021
      ".dump",    ".DUMP"      // added Apr 16 2021
    } ;

  char fnam[] = "ISFILE_TEXT" ;

  // ---------- BEGIN ---------

  // ----------------------
  // if fileName is blank then return 0 since we don't know what it is.
  if ( strlen(fileName) == 0 ) { return 0; }

  // ------ try each suffix --------
  for ( isuf=0; isuf < NSUFFIX_TEXT;  isuf++ ) {
    if ( strstr(fileName, SUFFIX_TEXT_LIST[isuf] ) != NULL )  
      { return 1 ; }
  }


  // May 4 2020: bail on fits extention
  if ( strstr(fileName,".fits") != NULL ) { return 0; }
  if ( strstr(fileName,".FITS") != NULL ) { return 0; }

  // if we get here, the extension is not known so open file
  // and check for for recognizable header key VARNAMES

  FILE *fp ;
  char ctmp[60];
  int NWD=0 ;

  if ( (fp = fopen(fileName,"rt")) == NULL ) {
    sprintf(MSGERR1, "Could not open TEXT FILE = ");
    sprintf(MSGERR2, "'%s'", fileName);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  while( (fscanf(fp, "%s", ctmp)) != EOF && NWD < 1000 ) {   
    if ( strcmp(ctmp,"VARNAMES:") == 0 ) { return 1; }
    NWD++ ;
  }

  return 0;

} // end of ISFILE_TEXT

// ==============================              
int ISFILE_MARZ(char *fileName) {

  // Created May 2 2020 by R.Kessler
  // To be a MARZ file, filename must
  //  * have fits extension
  //  * have marz or MARZ in the name.

  bool ISFITS = false ;
  bool ISMARZ = false ;
  // --------------- BEGIN ---------------

  // if fileName is blank then return 0 since we don't know what it is.
  if ( strlen(fileName) == 0 ) { return 0; }

  // ------ must have .fits suffix -------------

  if ( strstr(fileName, ".fits" ) != NULL )  { ISFITS = true ; }
  if ( !ISFITS ) { return 0; }

  // ----- must have marz or MARZ somewhere in file name ------
  if ( strstr(fileName, "marz" ) != NULL )  { ISMARZ = true ; }
  if ( strstr(fileName, "MARZ" ) != NULL )  { ISMARZ = true ; }
  if ( !ISMARZ ) { return 0; }

  // if we get here, it's a valid MARZ file
  return 1 ;
  
} // end ISFILE_MARZ
