// **********************************************
// Created Feb 23 2013 by R. Kessler
//
// functions to write and read ascii tables.
//
// - - - - - - - - - - - 
// HISTORY
// July 22 2013: include two new columns: CID and MJD
//
// Oct 15, 2013: in SNLCPAK_FILL_TEXT(void) fix bug closing each text file.
// Jan 22, 2014: write optional rest-frame filter and flux, and simRestFlux.
//
// Oct 27 2014: major overhaul for refactor
//
// Jan 10 2015: fix bug in SNTABLE_FILL_TEXT; add padded blank after CVAL
//              just before trim call.
//
// Apr 26 2015: 
//   For LCPLOT table, all SN written to one  [REPFIX].LCPLOT.TEXT file
//   instead of one file per SN.
//   New function CLOSE_TEXTFILE() to close LCPLOT file(s).
//
// Nov 7 2015: TEXTFILE_PREFIX is always a prefix, even if there is a dot.
//             See SNTABLE_CREATE_TEXT
//
// Nov 19 2015: undo Nov 7 change; if there's a dot, assume it's the 
//              whole file-name.
//
// Jul 17 2016:
//   + remove IDTABLE_OVERRIDE
//   + write forst 20 chars of .ptr_C into CVAL
//
// Dec 4 2016: allow char VERSION column (i.e., photometry) 
//
// Jan 14 2017: move ISFILE_TEXT to sntools_output.c
//
// Dec 2 2017: 
//  + minor overhaul to process .gz files. See *open_TEXTgz.
//    Beware that rewind() does not work on gzipped file.
//
// Dec 20 2017:
//   SNTABLE_NEVT_TEXT() reads lines and is much faster
//   READ_EXEC uses only fgets
//
// Apr 17 2019: in SNTABLE_NEVT_TEXT, rewind -> snana_rewind.
//
// Aug 25 2019: MXCHAR_LINE -> 2500 (was 2000)
//
// Oct 30 2019: begin adding optional ERRCALC column, but got
//    stuck passing flag to INIT_TEXTFILES. Also realized that
//    for fluxError maps we need more variables (PSF,SKY ..)
//    so better to use EPOCH table from ROOT or HBOOK file.
//    Bottom line is no change.
//
// Dec 19 2019: for sims, write SIM_FLUXCAL column to LCPLOT file.
//
// Jun 03 2020: if MARZ option is set, remove [z=xxx] from CCID
//              in TEXT SPECPLOT table.
//
// Jan 4 2021: MXCHAR_LINE -> 3200 (was 2500)
// **********************************************

char FILEPREFIX_TEXT[100];

#define MXTABLE_TEXT 10 
#define MXVAR_TEXT   MXVAR_TABLE
#define MXEPVAR_TEXT 50   // for light curve epoch
#define MXCHAR_LINE  3200 

#define OPT_FORMAT_KEY   1
#define OPT_FORMAT_CSV   2
#define OPT_FORMAT_COL   3
#define OPT_FORMAT_NONE  4  // no ascii output

#define SUFFIX_TEXT           "TEXT"   
#define SUFFIX_LCPLOT_TEXT    "LCPLOT.TEXT" 
#define SUFFIX_LCLIST_TEXT    "LCLIST.TEXT" 
#define SUFFIX_SPECPLOT_TEXT  "SPECPLOT.TEXT" 
#define SUFFIX_SPECLIST_TEXT  "SPECLIST.TEXT" 
#define TEXTMODE_rt           "rt"
#define TEXTMODE_wt           "wt"

#define MSKOPT_PARSE_WORDS_STRING 2 // must match same param in sntools.h

FILE *PTRFILE_TEXT ;                   // generic ascii file pointer
char FILENAME_TEXT[MXCHAR_FILENAME];   // name of opened text file
int  GZIPFLAG_TEXT;                    // gzipped or not

FILE *PTRFILE_LCLIST ;  // list-file for LC ascii files
FILE *PTRFILE_LCPLOT ;      
FILE *PTRFILE_SPECLIST ; 
FILE *PTRFILE_SPECPLOT ;    

int  NVAR_LCPLOT, NVAR_LCPLOT_REQ, NVAR_SPECPLOT ;
char VARNAME_SNLC[MXEPVAR_TEXT][40] ; 
char VARDEF_SNLC[MXEPVAR_TEXT][80] ; // short definition


// structure for writing table.
struct TABLEINFO_TEXT {
  int    NTABLE ;
  int    IDTABLE[MXTABLE_TEXT] ;
  char   TBNAME[MXTABLE_TEXT][40] ;

  FILE  *FP[MXTABLE_TEXT] ;
  char   FILENAME[MXTABLE_TEXT][MXCHAR_FILENAME] ;
  char   FORMAT[MXTABLE_TEXT][12] ;
  int    OPT_FORMAT[MXTABLE_TEXT] ;

  int    NVAR[MXTABLE_TEXT];
  char  *VARLIST[MXTABLE_TEXT];
  char **VARNAME[MXTABLE_TEXT];

  int    NFILL[MXTABLE_TEXT] ; // number of SNTABLE_FILL_TEXT calls.

  int        *ICAST[MXTABLE_TEXT] ;
  double    **ptr_D[MXTABLE_TEXT] ;
  float     **ptr_F[MXTABLE_TEXT] ;
  int       **ptr_I[MXTABLE_TEXT] ;
  short int **ptr_S[MXTABLE_TEXT] ;
  long long int  **ptr_L[MXTABLE_TEXT] ;
  char   **ptr_C[MXTABLE_TEXT] ;

} TABLEINFO_TEXT ;


// -----------------------------

#ifdef __cplusplus
extern"C" {          
#endif

  void INIT_TEXTFILES(char *PREFIX) ;  // analog of OPEN_ROOTFILE[HFILE]

  void SNTABLE_CREATE_TEXT(int IDTABLE, char *TBNAME, char *TEXT_FORMAT);

  void SNTABLE_ADDCOL_TEXT(int IDTABLE, void *PTRVAR, 
			   SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF) ;

  void SNTABLE_WRITE_HEADER_TEXT(int ITAB) ;

  void SNTABLE_FILL_TEXT(int IDTABLE);

  int  ITABLE_TEXT(int IDTABLE, char *FUNNAM, int OPT_ABORT) ;

  void formatFloat_TEXT(char *VARNAME, double VAL, char *CVAL); // return CVAL

  void OPEN_TEXTFILE(char *FILENAME, char *mode) ;
  void CLOSE_TEXTFILE(void);

  // fill SNLC plots
  void OPEN_TEXTFILE_LCLIST(char *PREFIX) ;  
  void open_textfile_lclist__(char *PREFIX) ;   
  void SNLCPAK_FILL_TEXT(void) ;
  void snlcpak_textLine(FILE *fp, int FLAG, int obs, int ifilt, int OUTFLAG);
  void SNLCPAK_WRITE_HEADER_TEXT(FILE *fp) ;

  // fill SPEC
  void OPEN_TEXTFILE_SPECLIST(char *PREFIX) ;  
  void open_textfile_speclist__(char *PREFIX) ;   
  void SPECPAK_FILL_TEXT(void) ;
  void specpak_textLine(FILE *fp, char *CCID, int ilam, int OUTFLAG);
  void SPECPAK_WRITE_HEADER_TEXT(void) ;

  void get_sepchar(int OPT_FORMAT, char *comment, char *sep) ;
  int  get_OPT_FORMAT(char *FORMAT, char *comment) ;

  void set_FILENAME_OVERRIDE_TEXT(int IDTABLE, char *fileName);
  void set_filename_override_text__(int *IDTABLE, char *fileName);

  int  SNTABLE_NEVT_TEXT(char *FILENAME);
  int  SNTABLE_NEVT_APPROX_TEXT(char *FILENAME, int NVAR);

  int  SNTABLE_READPREP_TEXT(void);
  int  SNTABLE_READ_EXEC_TEXT(void);
  void SNTABLE_CLOSE_TEXT(void) ;

  int validRowKey_TEXT(char *string) ;
  // xxx mvoed to sntools_output.h  int ICAST_for_textVar(char *varName) ;

  int count_varnames_TEXT();
  int get_varname_TEXT(int ivar, char *varName );

  // misc. sntools functions
  void  readint(FILE *fp, int nint, int *list) ;
  void  readchar(FILE *fp, char *clist) ;
  FILE *open_TEXTgz(char *FILENAME, const char *mode, int *GZIPFLAG);
  int   store_PARSE_WORDS(int OPT, char *FILENAME);
  void  get_PARSE_WORD(int langFlag, int iwd, char *word);
  void  trim_blank_spaces(char *string);
  void  debugexit(char *string);
  void  snana_rewind(FILE *fp, char *FILENAME, int GZIPFLAG);
#ifdef __cplusplus
}
#endif


// =============================================
//
//   BEGIN FUNCTIONS 
//
// =============================================


void INIT_TEXTFILES(char *PREFIX) {

  // May 2014
  // store prefix, and misc inits.
  // Cannot open table file(s) here because we don't know
  // which one(s). Table File(s) opened at SNTABLE_CREATE stage.
  // If lightcurve text-format has been specified, then call
  // function to open list-file.

  char fnam[] = "INIT_TEXTFILES" ;
  char *FMT, comment[200];

  sprintf(FILEPREFIX_TEXT,"%s", PREFIX);

  TABLEINFO_TEXT.NTABLE = 0 ;
  
  NVAR_LCPLOT = NVAR_LCPLOT_REQ = NVAR_SPECPLOT = 0 ; // for LCPLOT

  // check for lightcurve dump in ascii formats
  SNLCPAK_OUTPUT.OPT_TEXT_FORMAT = -9 ; // init no ascii light curves
  FMT = SNLCPAK_OUTPUT.TEXT_FORMAT ;
  if ( strlen(FMT) > 0 && strcmp(FMT,"none") != 0 ) {  
    sprintf(comment,"called from %s with PREFIX=%s", fnam, PREFIX);
    SNLCPAK_OUTPUT.OPT_TEXT_FORMAT = get_OPT_FORMAT(FMT,comment);
    OPEN_TEXTFILE_LCLIST(PREFIX) ; 
  }


  SPECPAK_OUTPUT.OPT_TEXT_FORMAT = -9 ; // init no ascii light curves
  FMT = SPECPAK_OUTPUT.TEXT_FORMAT ;
  if ( strlen(FMT) > 0 && strcmp(FMT,"none") != 0 ) {  
    sprintf(comment,"called from %s with PREFIX=%s", fnam, PREFIX);
    SPECPAK_OUTPUT.OPT_TEXT_FORMAT = get_OPT_FORMAT(FMT,comment);
    OPEN_TEXTFILE_SPECLIST(PREFIX);
  }


} // end of INIT_TEXTFILES


// ============================================
void SNTABLE_CREATE_TEXT(int IDTABLE, char *TBNAME, char *TEXT_FORMAT) {

  // open file and store pointer to file.
  // File name is [PREFIX].[NAME].TEXT.
  // Store TEXT_FORMAT to be used when writing.

  int  NTAB, ivar, OPT_FORMAT, GZIPFLAG ;
  char FILENAME[MXCHAR_FILENAME], comment[100];
  char fnam[] = "SNTABLE_CREATE_TEXT" ;

  // --------- BEGIN ---------

  NTAB = TABLEINFO_TEXT.NTABLE ;

  // increment number of tables
  TABLEINFO_TEXT.NTABLE++ ;

  printf("  %s: init %s TEXT-table in %s format. \n", 
	 fnam, TBNAME, TEXT_FORMAT); fflush(stdout);
 
  // construct file name from stored prefix and table name.

  sprintf(FILENAME, "%s.%s.%s", FILEPREFIX_TEXT,TBNAME,SUFFIX_TEXT);


  if ( strchr(FILEPREFIX_TEXT,'.') == NULL ) {   
    // no dot --> add on table and suffix to file name
    sprintf(FILENAME, "%s.%s.%s", FILEPREFIX_TEXT,TBNAME,SUFFIX_TEXT);
  }
  else {
    // if there is a dot in the prefix, assume it's the full filename
    sprintf(FILENAME, "%s", FILEPREFIX_TEXT);
  }


  // store other info
  TABLEINFO_TEXT.IDTABLE[NTAB] = IDTABLE ;
  TABLEINFO_TEXT.NFILL[NTAB]   = 0 ;
  sprintf(TABLEINFO_TEXT.TBNAME[NTAB],   "%s", TBNAME);
  sprintf(TABLEINFO_TEXT.FILENAME[NTAB], "%s", FILENAME);
  sprintf(TABLEINFO_TEXT.FORMAT[NTAB],   "%s", TEXT_FORMAT);  

  sprintf(comment,"Called from %s with IDTABLE=%d (%s)", 
	  fnam, IDTABLE, TBNAME);
  OPT_FORMAT = get_OPT_FORMAT(TEXT_FORMAT,comment);
  TABLEINFO_TEXT.OPT_FORMAT[NTAB] = OPT_FORMAT ; 

  // if format = 'none', return here after defining this text-table
  // so that it can be skipped in the ADDCOL and FILL functions.
  if ( OPT_FORMAT == OPT_FORMAT_NONE ) { return ; }

  // init NVAR and VARNAMES list
  TABLEINFO_TEXT.NVAR[NTAB] = 0 ;

  // - - - - - - - - - - - -
   // open text file  
  TABLEINFO_TEXT.FP[NTAB] = open_TEXTgz(FILENAME,TEXTMODE_wt, &GZIPFLAG) ;
  if ( !TABLEINFO_TEXT.FP[NTAB] ) {
    sprintf(MSGERR1, "Could not open TEXT FILE = ");
    sprintf(MSGERR2, "%s", FILENAME);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }


  TABLEINFO_TEXT.VARLIST[NTAB] = (char*)malloc( MXCHAR_LINE*sizeof(char));
  TABLEINFO_TEXT.VARLIST[NTAB][0] = 0 ;

  TABLEINFO_TEXT.VARNAME[NTAB] = (char**)malloc(MXVAR_TABLE * sizeof(char*));
  for(ivar=0; ivar < MXVAR_TABLE; ivar++ ) {
    TABLEINFO_TEXT.VARNAME[NTAB][ivar] = (char*)malloc(40 * sizeof(char));
  }


  // allocate memory for pointers
  int MX = MXVAR_TABLE ;
  TABLEINFO_TEXT.ICAST[NTAB] = (int*)     malloc ( MX * sizeof(int )    );  
  TABLEINFO_TEXT.ptr_I[NTAB] = (int**)    malloc ( MX * sizeof(int*)    );
  TABLEINFO_TEXT.ptr_S[NTAB] = 
    (short int**) malloc ( MX * sizeof(short int*)    );
  TABLEINFO_TEXT.ptr_L[NTAB] = 
    (long long int**) malloc( MX*sizeof(long long int*) );
  TABLEINFO_TEXT.ptr_F[NTAB] = (float**)  malloc ( MX * sizeof(float*)  );
  TABLEINFO_TEXT.ptr_D[NTAB] = (double**) malloc ( MX * sizeof(double*) );  
  TABLEINFO_TEXT.ptr_C[NTAB] = (char**)   malloc ( MX * sizeof(char*)   );


  int IVAR;
  for(IVAR=0; IVAR < MXVAR_TABLE; IVAR++ )
    { TABLEINFO_TEXT.ICAST[NTAB][IVAR] = -9 ; }


} // end of SNTABLE_CREATE_TEXT



// =============================================
void SNTABLE_ADDCOL_TEXT(int IDTABLE, void *PTRVAR, 
			 SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF) {

  // store column info name and pointer.
  // Note that only a subset of all variables is kept;
  // The HBOOK & ROOT tables keep everything.
  //
  // Apr 12 2019: abort of VARLIST > MXCHAR_LINE

  int  IVAR, ivar, ITAB, ICAST, LENV ;
  char VARNAME[MXCHAR_VARNAME], *varList, *TBNAME ;
  char fnam[] = "SNTABLE_ADDCOL_TEXT" ;

  // ------------- BEGIN --------------


  // find sparse index ITAB for this IDTABLE
  ITAB   = ITABLE_TEXT(IDTABLE,fnam, 0 );
  TBNAME = TABLEINFO_TEXT.TBNAME[ITAB] ;

  if ( ITAB < 0 ) {
    print_preAbort_banner(fnam);
    printf("   VARLIST_ORIG = '%s' \n", ADDCOL_VARDEF->VARLIST_ORIG);
    printf("   NVAR = %d \n", ADDCOL_VARDEF->NVAR );
    sprintf(MSGERR1, "Could not find text table with ");
    sprintf(MSGERR2, "IDTABLE = %d  TBNAME=%s",   IDTABLE, TBNAME ); 
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }


  int OPT_FORMAT =   TABLEINFO_TEXT.OPT_FORMAT[ITAB] ;
  if ( OPT_FORMAT == OPT_FORMAT_NONE ) { return ; }

  for(ivar=0 ; ivar < ADDCOL_VARDEF->NVAR; ivar++ ) {

    TABLEINFO_TEXT.NVAR[ITAB]++ ;
    IVAR =  TABLEINFO_TEXT.NVAR[ITAB] - 1 ; // global counter
    if ( IVAR >= MXVAR_TEXT ) { continue ; }

    // construct VARNAMES list for header
    sprintf(VARNAME, "%s", ADDCOL_VARDEF->VARNAME[ivar] );    

    // replace CCID -> CID in text file
    if ( strcmp(VARNAME,"CCID") == 0 )  { sprintf(VARNAME,"CID") ; }

    varList = TABLEINFO_TEXT.VARLIST[ITAB] ;

    // check varlist bound
    LENV = strlen(varList) + strlen(VARNAME) + 1 ;
    if ( LENV >= MXCHAR_LINE ) {
      print_preAbort_banner(fnam);
      printf(" varList = '%s' \n", varList);
      printf(" VARNAME to add: '%s' \n", VARNAME);
      sprintf(MSGERR1, "len(VARLIST)=%d exceeds bound of %d",
	      LENV, MXCHAR_LINE );
      sprintf(MSGERR2, "after adding VARNAME='%s'", VARNAME );
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
    }
   
    strcat(varList," "); strcat(varList,VARNAME);
    sprintf(TABLEINFO_TEXT.VARNAME[ITAB][IVAR],"%s", VARNAME);

    // store pointer based on cast.

    ICAST = ADDCOL_VARDEF->ICAST[ivar] ;
    TABLEINFO_TEXT.ICAST[ITAB][IVAR] = ICAST ;    

    if ( ICAST == ICAST_D ) 
      { TABLEINFO_TEXT.ptr_D[ITAB][IVAR] = (double*)PTRVAR + ivar ; }
    else if ( ICAST == ICAST_F ) 
      { TABLEINFO_TEXT.ptr_F[ITAB][IVAR] = (float*)PTRVAR + ivar ; }
    else if ( ICAST == ICAST_I ) 
      { TABLEINFO_TEXT.ptr_I[ITAB][IVAR] = (int*)PTRVAR + ivar ; }    
    else if ( ICAST == ICAST_S ) 
      { TABLEINFO_TEXT.ptr_S[ITAB][IVAR] = (short int*)PTRVAR + ivar ; }    
    else if ( ICAST == ICAST_L ) 
      { TABLEINFO_TEXT.ptr_L[ITAB][IVAR] = (long long int*)PTRVAR+ivar;}    
    else if ( ICAST == ICAST_C )
      { TABLEINFO_TEXT.ptr_C[ITAB][IVAR] = (char*)PTRVAR ; }
    
  } // end of ivar


} // end of SNTABLE_ADDCOL_TEXT


// =============================================
int ITABLE_TEXT(int IDTABLE, char *FUNNAM, int OPT_ABORT ) {

  // return sparse table index for IDTABLE
  // On error, print calling function FUNNAM.
 
  int i, ITAB;
  char fnam[] = "ITABLE_TEXT" ;

  ITAB = -9 ;
  for(i=0; i < TABLEINFO_TEXT.NTABLE ; i++ ) {   
    if ( IDTABLE == TABLEINFO_TEXT.IDTABLE[i] ) { ITAB = i ; }
  }

  if ( ITAB < 0 && OPT_ABORT > 0 ) {
    sprintf(MSGERR1, "%s could not find text table with ", FUNNAM );
    sprintf(MSGERR2, "IDTABLE = %d", IDTABLE ); 
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }

  return ITAB ;

} // end of ITABLE_TEXT


// ==================================================
void formatFloat_TEXT(char *VARNAME, double VAL, char *VALSTRING) {

  // return VALSTRING formatted according to the input value VAL
  // Feb 5 2015: write .6f for RA and DEC

  char c1[4];
  long long int IVAL8 ;
  double ABSVAL = fabs(VAL) ;

  // -------------- BEGIN ------------

  sprintf(c1, "%c", VARNAME[0]); // 1st char
  IVAL8 = (long long int)VAL ;

  /* xxxx
  if ( strstr(VARNAME,"HOST") != NULL ) {
    printf(" xxx VAL=%le  IVAL8=%lld  diff = %le \n",
	   VAL, IVAL8, ((double)IVAL8-VAL) );
  }
  xxxxx */

  if ( strstr(VARNAME,"MJD") != NULL ) {
    sprintf(VALSTRING, "%.3f", VAL); // some kind of MJD
  }

  else if ( strstr(VARNAME,"RA") != NULL && ABSVAL<400 ) {
    sprintf(VALSTRING, "%.6f", VAL);    
  }
  else if ( strstr(VARNAME,"DEC") != NULL && ABSVAL<400 ) {
    sprintf(VALSTRING, "%.6f", VAL);    
  }
  else if ( strcmp(c1,"z") == 0 ) {
    // probably a redshift
    sprintf(VALSTRING, "%.5f", VAL);    
  }
  else if ( VAL > 1.0E10 ) {
    sprintf(VALSTRING, "%.4le", VAL);  // Aug 2019
  }
  else if ( (VAL - IVAL8) == 0.0 ) {
    // it's really an integer
    sprintf(VALSTRING, "%lld", IVAL8);  
  }
  else if ( VAL > 0.1 ) {
    sprintf(VALSTRING, "%.5f", VAL);  // .4f -> .5f Nov 10 2018
  }
  else {
    // small number -> expon.
    sprintf(VALSTRING, "%.5e", VAL);    
  }

  
  /*
  printf(" xxx '%s' = %f --> '%s' \n",	 VARNAME, VAL, VALSTRING);
  fflush(stdout); // xxxxxx
  */

  return ;

} // end of formatFloat_TEXT

// ==================================================
void SNTABLE_FILL_TEXT(int IDTABLE) {

  // write ROW to text file for this IDTABLE.
  // On NFILL=0, also write header info based on format.
  //
  // Jul 17 2016: write first 20 chars of ptr_C into CVAL
  //              to avoid crazy-long strings.
  //
  // Jun 24 2017: ROW[1000] -> ROW[2000]

  int ITAB, NFILL, NVAR, IVAR, ICAST, OPT_FORMAT ;

  FILE *FP ;
  char ROW[MXCHAR_LINE], CVAL[80], *FORMAT, *VARNAME, sep[4], comment[200] ;
  char fnam[] = "SNTABLE_FILL_TEXT" ;

  // ------------- BEGIN ------------

  // --------------------------------------------
  // find sparse index ITAB for this IDTABLE
  ITAB = ITABLE_TEXT(IDTABLE,fnam,1); // abort if no table.
  
  FP         = TABLEINFO_TEXT.FP[ITAB]  ;
  NFILL      = TABLEINFO_TEXT.NFILL[ITAB];
  OPT_FORMAT = TABLEINFO_TEXT.OPT_FORMAT[ITAB];
  FORMAT     = TABLEINFO_TEXT.FORMAT[ITAB];

  if ( OPT_FORMAT == OPT_FORMAT_NONE ) { return ; }

  // write header on first FILL since this is the only way
  // to know that all columns are defined.
  if ( NFILL == 0 ) { SNTABLE_WRITE_HEADER_TEXT(ITAB) ; }
    
  ROW[0] = 0 ;
  CVAL[0] = 0 ;

  sprintf(comment,"called by %s for IDTABLE=%d and FORMAT='%s'",
	  fnam, IDTABLE, FORMAT);
  get_sepchar(OPT_FORMAT, comment, sep) ; // return sep

  long long int VAL_L;
  double VAL_D ;
  float  VAL_F ;
  int    VAL_I ;
  short int VAL_S;

  NVAR = TABLEINFO_TEXT.NVAR[ITAB] ; 

  for ( IVAR=0; IVAR < NVAR ; IVAR++ ) {
    ICAST   = TABLEINFO_TEXT.ICAST[ITAB][IVAR] ;
    VARNAME = TABLEINFO_TEXT.VARNAME[ITAB][IVAR] ;
      
    if (ICAST < 0 ) {
      sprintf(MSGERR1, "Undefined ICAST for IVAR=%d  IDTABLE=%d", 
	      IVAR, IDTABLE);
      sprintf(MSGERR2, "Something is messed up."); 
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
    }

    if ( ICAST == ICAST_D ) {
      VAL_D = *TABLEINFO_TEXT.ptr_D[ITAB][IVAR] ;
      formatFloat_TEXT(VARNAME, VAL_D, CVAL); // return CVAL
    }
    else if ( ICAST == ICAST_F ) {
      VAL_F = *TABLEINFO_TEXT.ptr_F[ITAB][IVAR] ;
      formatFloat_TEXT(VARNAME, (double)VAL_F, CVAL); // return CVAL
    }
    else if ( ICAST == ICAST_I ) {
      VAL_I = *TABLEINFO_TEXT.ptr_I[ITAB][IVAR] ;
      sprintf(CVAL, "%d", VAL_I );
    }
    else if ( ICAST == ICAST_S ) {
      VAL_S = *TABLEINFO_TEXT.ptr_S[ITAB][IVAR] ;
      sprintf(CVAL, "%d", VAL_S );
    }
    else if ( ICAST == ICAST_L ) {
      VAL_L = *TABLEINFO_TEXT.ptr_L[ITAB][IVAR] ;
      sprintf(CVAL, "%lld", VAL_L );
    }
    else if ( ICAST == ICAST_C ) {
      // leave extra blank space at the end of CVAL
      // to ensure that the trim function works 
      sprintf(CVAL,"%.*s ", MXCHAR_CCID, TABLEINFO_TEXT.ptr_C[ITAB][IVAR] );
      trim_blank_spaces(CVAL) ;
    }

    strcat(ROW,CVAL);  
    if ( IVAR < NVAR-1 ) { strcat(ROW,sep); }

  } // end of IVAR loop


  if ( OPT_FORMAT == OPT_FORMAT_KEY ) {
    fprintf(FP, "SN: %s\n", ROW);
  }
  else if ( OPT_FORMAT == OPT_FORMAT_CSV ) {
    fprintf(FP, "%s\n", ROW);
  }
  else if ( OPT_FORMAT == OPT_FORMAT_COL ) {
    fprintf(FP, "%s\n", ROW);
  }

  fflush(FP);

  // increment number of FILL calls
  TABLEINFO_TEXT.NFILL[ITAB]++ ;

} // end of SNTABLE_FILL_TEXT

// ===============================================================
int get_OPT_FORMAT(char *FORMAT, char *comment) {

  // get integer index for text *FORMAT.
  // *comment is for abort message.

  char fnam[] = "get_OPT_FORMAT" ;

  if ( strcmp_ignoreCase(FORMAT,(char*)"key") == 0 ) 
    { return OPT_FORMAT_KEY ; }
  else if ( strcmp_ignoreCase(FORMAT,(char*)"csv") == 0 ) 
    { return OPT_FORMAT_CSV ; }
  else if ( strcmp_ignoreCase(FORMAT,(char*)"col") == 0 ) 
    { return OPT_FORMAT_COL ; }
  else if ( strcmp_ignoreCase(FORMAT,(char*)"none") == 0 ) 
    { return OPT_FORMAT_NONE ; }
  else {
    sprintf(MSGERR1, "Invalid TEXT_FORMAT = '%s'", FORMAT);
    sprintf(MSGERR2, "%s", comment);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  return(-9); // can never get here

} // end of get_OPT_FORMAT

// ===============================================================
void get_sepchar(int OPT_FORMAT, char *comment, char *sep) {

  // Created Sep 2014
  // for input OPT_FORMAT, return *sep string to write
  // between values: either comma or blank string.
  // *comment is used for error message.

  char fnam[] = "getformat_sepchar" ;
  //  char *sep   = (char*) malloc( sizeof(char) * 4) ;

  if ( OPT_FORMAT == OPT_FORMAT_KEY ) 
    { sprintf(sep," "); }
  else if ( OPT_FORMAT == OPT_FORMAT_CSV ) 
    { sprintf(sep,","); }  
  else if ( OPT_FORMAT == OPT_FORMAT_COL ) 
    { sprintf(sep," "); }
  else {
    sprintf(MSGERR1, "Unknown OPT_FORMAT=%d", OPT_FORMAT);
    sprintf(MSGERR2, "%s", comment);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }

} // end of get_sepchar


// =============================================
void SNTABLE_WRITE_HEADER_TEXT(int ITAB) {

  // Write optional header info at top of file.
  // See FORMAT string below.

  int   NVAR, IVAR, OPT_FORMAT, IDTABLE ;
  char *FORMAT, *VARLIST, *TBNAME;
  char fnam[] = "SNTABLE_WRITE_HEADER_TEXT" ;
  FILE *FP ;
  // ------------- BEGIN --------------

  OPT_FORMAT = TABLEINFO_TEXT.OPT_FORMAT[ITAB] ;
  FORMAT     = TABLEINFO_TEXT.FORMAT[ITAB] ;
  NVAR       = TABLEINFO_TEXT.NVAR[ITAB] ;
  VARLIST    = TABLEINFO_TEXT.VARLIST[ITAB] ;
  TBNAME     = TABLEINFO_TEXT.TBNAME[ITAB] ;
  IDTABLE    = TABLEINFO_TEXT.IDTABLE[ITAB] ;
  FP         = TABLEINFO_TEXT.FP[ITAB] ; 


  if ( NVAR >= MXVAR_TEXT ) {
    sprintf(MSGERR1, "NVAR=%d exceeds bound of MXVAR_TEXT=%d",
	    NVAR, MXVAR_TEXT);
    sprintf(MSGERR2, "for IDTABLE=%d  TableName='%s' ",
	    TABLEINFO_TEXT.IDTABLE[ITAB], TABLEINFO_TEXT.TBNAME[ITAB] ) ;
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  // for KEY format, write comments at top
  if ( OPT_FORMAT == OPT_FORMAT_KEY ) {

    // Oct 23 2014: print comment lines
    int iline ;
    for(iline=0; iline < NLINE_TABLECOMMENT; iline++ ) {
      fprintf(FP, "# %s \n", LINE_TABLECOMMENT[iline] );
    }
    fflush(FP);
    
    fprintf(FP, "# SNANA_VERSION : %s \n", SNANA_VERSION) ;
    fprintf(FP, "# TABLE NAME: %s \n", TBNAME);
    fprintf(FP, "# \n" );
#ifdef TEXTFILE_NVAR
    fprintf(FP, "NVAR: %d \n", NVAR ); 
#endif
    fprintf(FP, "VARNAMES: %s \n", VARLIST );
    fprintf(FP, "#\n" );
    fflush(FP);	
  }

  else if ( OPT_FORMAT == OPT_FORMAT_CSV ) {

    char CSVLIST[1000], *varName;
    for(IVAR=0; IVAR < NVAR ; IVAR++ ) {
      varName = TABLEINFO_TEXT.VARNAME[ITAB][IVAR] ;

      if ( IVAR == 0 ) 
	{ sprintf(CSVLIST, "%s", varName ); }
      else
	{ sprintf(CSVLIST, "%s,%s", CSVLIST, varName ); }
    }

    fprintf(FP,"%s\n", CSVLIST); 
  }
  else if ( OPT_FORMAT == OPT_FORMAT_COL ) {
    // do nothing
  }
  else {
    sprintf(MSGERR1, "Unknown OPT_FORMAT=%d (FORMAT='%s')",
	    OPT_FORMAT, FORMAT);
    sprintf(MSGERR2, "writing header for IDTABLE=%d   TBNAME='%s'", 
	    IDTABLE, TBNAME ); 
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }


} // end of SNTABLE_WRITE_HEADER_TEXT


void OPEN_TEXTFILE(char *FILENAME, char *mode) {

  // Created Oct 2014
  // e.g., mode = 'rt' or 'wt'

  char fnam[] = "OPEN_TEXTFILE" ;

  PTRFILE_TEXT = open_TEXTgz(FILENAME,mode, &GZIPFLAG_TEXT);

  sprintf(FILENAME_TEXT, "%s", FILENAME);  // Dec 2 2017

  if ( !PTRFILE_TEXT ) {
    sprintf(MSGERR1, "Could not open text-file in mode='%s' : ", mode);
    sprintf(MSGERR2, "%s", FILENAME);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  char msg[MXCHAR_FILENAME+20];
  sprintf(msg,"%s: %s", fnam, FILENAME);
  printf("   %s \n", msg);


} // end of OPEN_TEXTFILE


void CLOSE_TEXTFILE(void) {

  // Created Apr 2015
  // Close open TEXT file(s) that were created as output.
  // Beware: does not close TEXT files opened for reading.
  
  int itab;
  char *FNAM;
  FILE *FP;
  //  char fnam[] = "CLOSE_TEXTFILE" ;

  // ------------- BEGIN -------------

  for(itab=0; itab < TABLEINFO_TEXT.NTABLE; itab++ ) {
      FP   = TABLEINFO_TEXT.FP[itab] ;
      FNAM = TABLEINFO_TEXT.FILENAME[itab] ;
      if ( FP != NULL ) {
	printf("   Close %s \n", FNAM); fflush(stdout);
	fclose(FP);
      }
  }

  if ( NVAR_LCPLOT ) {
    printf("   Close TEXT LCPLOT files \n"); fflush(stdout);
    fclose(PTRFILE_LCLIST);
    fclose(PTRFILE_LCPLOT);
  }

  if ( NVAR_SPECPLOT ) {
    printf("   Close TEXT SPECPLOT files \n"); fflush(stdout);
    fclose(PTRFILE_SPECLIST);
    fclose(PTRFILE_SPECPLOT);
  }

} // end of CLOSE_TEXTFILE

// =============================================
void OPEN_TEXTFILE_LCLIST(char *PREFIX) {
  
  // Sep 2014: (changed from OPEN_TEXTFILE_LEGACY)
  // Open ascii LIST file for each SN light curves dumped in text format,
  // and define variable names and definitions in string arrays.
  // Used for SNTABLE_LIST = 'LCPLOT(text:[fmt])' and also for
  // legacy option LDMP_SNLFUX=T .
  //
  // Apr 26 2015: open LCPLOT file to store all SN LC in one file.
  //              See lcplotFile[]
  //

  int GZIPFLAG ;
  char listFile[MXCHAR_FILENAME];
  char lcplotFile[MXCHAR_FILENAME];
  char fnam[] = "OPEN_TEXTFILE_LCLIST" ;

  // ------------- BEGIN -----------

  sprintf(listFile,   "%s.%s", PREFIX, SUFFIX_LCLIST_TEXT );
  sprintf(lcplotFile, "%s.%s", PREFIX, SUFFIX_LCPLOT_TEXT );

  PTRFILE_LCLIST = open_TEXTgz(listFile,TEXTMODE_wt, &GZIPFLAG);
  if ( !PTRFILE_LCLIST ) {
    sprintf(MSGERR1, "Could not open ascii LC list-file = ");
    sprintf(MSGERR2, "%s", listFile);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  //  PTRFILE_LCPLOT = fopen(lcplotFile,"wt");
  PTRFILE_LCPLOT = open_TEXTgz(lcplotFile,TEXTMODE_wt, &GZIPFLAG );
  if ( !PTRFILE_LCPLOT ) {
    sprintf(MSGERR1, "Could not open ascii LCPLOT file = ");
    sprintf(MSGERR2, "%s", lcplotFile);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }


  // -------------------------------------------
  // define columns, including optional rest-frame variables

  int NVAR = 0 ;

  sprintf(VARNAME_SNLC[NVAR], "CID" );           
  sprintf(VARDEF_SNLC[NVAR],  "Candidate ID");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "MJD" );  
  sprintf(VARDEF_SNLC[NVAR],  "time of observation");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "Tobs" );          
  sprintf(VARDEF_SNLC[NVAR],  "MJD-PKMJD");
  NVAR++ ;


  if ( SNLCPAK_OUTPUT.SIMFLAG ) {  // Dec 2019
    sprintf(VARNAME_SNLC[NVAR], "SIM_FLUXCAL" );  
    sprintf(VARDEF_SNLC[NVAR],  "true flux"   );
    NVAR++ ;
  }

  sprintf(VARNAME_SNLC[NVAR], "FLUXCAL" );  
  sprintf(VARDEF_SNLC[NVAR],  "calibrated flux");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "FLUXCAL_ERR" );   
  sprintf(VARDEF_SNLC[NVAR],  "error on calibrated flux");
  NVAR++ ;


  sprintf(VARNAME_SNLC[NVAR], "DATAFLAG" );      
  sprintf(VARDEF_SNLC[NVAR],  "1=>data, 0=>fit, -1=>data rejected in fit");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "BAND" );          
  sprintf(VARDEF_SNLC[NVAR],  "1 char band/filter");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "CHI2" );          
  sprintf(VARDEF_SNLC[NVAR],  "data-fit chi2");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "IFIT" );     
  sprintf(VARDEF_SNLC[NVAR],  "fit number out of %d fits per SN", 
	  SNLCPAK_OUTPUT.NFIT_PER_SN );
  NVAR++ ;

  NVAR_LCPLOT_REQ = NVAR ; // required

  // - - - - - -

  /*
  // optional if ERRCALC is requested
  sprintf(VARNAME_SNLC[NVAR], "ERRCALC" );     
  sprintf(VARDEF_SNLC[NVAR],  "optional calculated error");
  NVAR++ ;
  */

  // below are optional if using rest-frame model
  sprintf(VARNAME_SNLC[NVAR], "BAND_REST" );     
  sprintf(VARDEF_SNLC[NVAR],  "optional rest-frame band");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "FLUX_REST" );   
  sprintf(VARDEF_SNLC[NVAR],  "optional rest-frame flux");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "KCOR" );   
  sprintf(VARDEF_SNLC[NVAR],  "optional kcor");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "AVWARP" );   
  sprintf(VARDEF_SNLC[NVAR],  "optional AVwarp");
  NVAR++ ;

  sprintf(VARNAME_SNLC[NVAR], "SIM_FLUX_REST" ); 
  sprintf(VARDEF_SNLC[NVAR],  "optional rest-frame SIM flux");
  NVAR++ ;

  NVAR_LCPLOT = NVAR;

  // ----------------------------------
  fprintf(PTRFILE_LCLIST,
	  "# column defintion in each lightcurve file: \n");

  int ivar;
  for(ivar=0; ivar < NVAR_LCPLOT ; ivar ++ ) {
    fprintf(PTRFILE_LCLIST,
	    "# col %2d -> %-14s :  %s \n", 
	    ivar+1, VARNAME_SNLC[ivar], VARDEF_SNLC[ivar] );
  }

#ifdef TEXTFILE_NVAR
  fprintf(PTRFILE_LCLIST, "\nNVAR: 2\n");
#endif
  fprintf(PTRFILE_LCLIST, "VARNAMES: CID IFIT\n");
  fflush(PTRFILE_LCLIST);

  SNLCPAK_USE_TEXT = true ;
  sprintf( SNLCPAK_OUTPUT.TEXTFILE_PREFIX, "%s", PREFIX );


  char banner[200];
  sprintf(banner,"%s: open text list-file with PREFIX=%s", 
	  fnam, PREFIX);
  print_banner(banner);


} // end of OPEN_TEXTFILE_LCLIST

void open_textfile_lclist__(char *PREFIX) {
  OPEN_TEXTFILE_LCLIST(PREFIX);
}

// ===========================================
int SNTABLE_NEVT_TEXT(char *FILENAME) {

  // Created Oct 13 2014
  // return number of table rows in ascii file.
  // Note that the table name is not required as an argument
  // because there can only be 1 table per file.
  // If FILENAME is blank, then assume file has already been opened.
  //
  // Dec 20 2017: huge speed-up reading lines instead of words.
  //
  // Apr 17 2019: rewind -> snana_rewind

  int NROW, LENF, GZIPFLAG ;
  FILE *fp ;
  char  LINE[MXCHAR_LINE], *ptrtok ;
  char fnam[] = "SNTABLE_NEVT_TEXT" ;

  // ------------ BEGIN --------------

  NROW = 0 ;
  LENF = strlen(FILENAME) ;

  if ( LENF > 0 ) 
    { fp = open_TEXTgz(FILENAME,TEXTMODE_rt, &GZIPFLAG); }
  else
    { fp = PTRFILE_TEXT ; } // already opened


  if ( !fp ) {
    sprintf(MSGERR1, "Could not open ascii table file: ");
    sprintf(MSGERR2, "'%s' ", FILENAME);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  while ( fgets(LINE,MXCHAR_LINE, fp) != NULL ) {
    ptrtok = strtok(LINE, " ");
    if ( validRowKey_TEXT(ptrtok) ) { NROW++ ; }
  }


  //  printf(" xxx %s: LENF=%d FILENAME=%s \n", fnam, LENF, FILENAME);

  if ( LENF > 0 )  { 
    fclose(fp); 
  }
  else {
    snana_rewind(fp, FILENAME_TEXT, GZIPFLAG_TEXT);
    //    rewind(fp); 
  } 

  return NROW ;

} // end of SNTABLE_NEVT_TEXT

// ========================================================
int SNTABLE_NEVT_APPROX_TEXT(char *FILENAME, int NVAR) {

  // Dec 1 2017
  // [code moved from combine_fitres.c]
  // estimate number of rows in text file based on its size;
  // avoids reading entire file.

  int istat, isize, NEVT_APPROX, gzFlag=0;
  char gzFile[MXCHAR_FILENAME];
  struct stat statbuf ;

  // ----------- BEGIN ---------

  if ( strstr(FILENAME,".gz") != NULL ) { gzFlag=1; }
  istat = stat(FILENAME, &statbuf );

  // check for gz file
  if ( istat != 0 ) {
    sprintf(gzFile,"%s.gz", FILENAME);
    istat = stat(gzFile, &statbuf );
    gzFlag=1;
  }
  

  isize = statbuf.st_size ;
  NEVT_APPROX = isize / (4*NVAR) ; // approx number of rows

  if ( gzFlag ) { NEVT_APPROX *= 2.5 ; }

  /* 
  printf(" xxx NEVT_APPROX=%d for file='%s' (isize=%d)\n", 
	 NEVT_APPROX, FILENAME, isize );
  */

  return(NEVT_APPROX);

}  // end SNTABLE_NEVT_APPROX_TEXT

// ==========================================
int  SNTABLE_READPREP_TEXT(void) {

  int NVAR, ivar, ISTAT, FOUNDKEY, NRD, GZIPFLAG ; 
  FILE *FP ;
  char ctmp[MXCHAR_FILENAME], *VARNAME, *VARLIST;
  char fnam[]=  "SNTABLE_READPREP_TEXT" ;

  // ---------- BEGIN -------------

  FP    = PTRFILE_TEXT ;
  sprintf(ctmp,"BLANK"); 
  ISTAT = 999;
  NVAR = FOUNDKEY = NRD = 0 ;

  while ( validRowKey_TEXT(ctmp) == 0 && ISTAT != EOF ) {

    NRD++ ;
    
    ISTAT = fscanf(FP, "%s", ctmp) ;
    if ( ISTAT == EOF && FOUNDKEY == 0 ) {
      sprintf(MSGERR1,"%s", 
	      "End of file reached without finding any HEADER keys.");
      sprintf(MSGERR2,"%s", "Check text file.");
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
    }

#ifndef TEXTFILE_NVAR
    // no NVAR key, so read entire VARNAMES line and parse VARLIST
    if ( strcmp(ctmp,"VARNAMES:") == 0 ) {
      FOUNDKEY = 1;
      VARLIST = (char*)malloc( MXCHAR_LINE * sizeof(char) );
      fgets(VARLIST, MXCHAR_LINE, FP );
      NVAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,VARLIST);
      free(VARLIST);
      for ( ivar=0; ivar < NVAR; ivar++ ) {
	VARNAME = READTABLE_POINTERS.VARNAME[ivar] ;
	get_PARSE_WORD(0,ivar,VARNAME);
	// printf(" xxx VARNAME[%2d] = '%s' \n", ivar, VARNAME);
	READTABLE_POINTERS.ICAST_STORE[ivar] = ICAST_for_textVar(VARNAME);
	READTABLE_POINTERS.ICAST_READ[ivar]  = ICAST_for_textVar(VARNAME);
      }
    }
#endif



#ifdef TEXTFILE_NVAR
    if ( strcmp(ctmp,"NVAR:") == 0 ) 
      { readint( FP, 1, &NVAR ); FOUNDKEY = 1 ; }

    if ( strcmp(ctmp,"VARNAMES:") == 0 ) {
      for ( ivar=0; ivar < NVAR; ivar++ ) {
	VARNAME = READTABLE_POINTERS.VARNAME[ivar] ;
	readchar(FP, VARNAME );
	READTABLE_POINTERS.ICAST_STORE[ivar] = ICAST_for_textVar(VARNAME);
	READTABLE_POINTERS.ICAST_READ[ivar]  = ICAST_for_textVar(VARNAME);
      } 
    }  // not row key
#endif

    // abort if a valid KEY cannot be found
    if ( FOUNDKEY == 0 && NRD >= 10000 ) {
      sprintf(MSGERR1,"Cannot find valid HEADER key after NRD=%d items.",NRD);
      sprintf(MSGERR2,"%s", "Check text file.");
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
    }

  } // end while


  // reset NVARTOT_FITRES to exclude CID
  READTABLE_POINTERS.NVAR_TOT = NVAR ;

  // rewind, but do not close file.
  // xxx mark delete  rewind(FP);

  // close file and re-open because rewind does not work on
  // gzipped files 
  fclose(PTRFILE_TEXT);
  PTRFILE_TEXT = open_TEXTgz(FILENAME_TEXT, TEXTMODE_rt, &GZIPFLAG);

  return NVAR;

} // end of SNTABLE_READPREP_TEXT

// ==============================================
int SNTABLE_READ_EXEC_TEXT(void) {

  // Oct 2014:  
  // Execute ascii read over all rows and fill
  // pointers passed previously to SNTABLE_READPREP_VARDEF.
  // Functions returns number of rows read. 
  //
  // July 29 2016: abort on NVAR key with different value.
  // Dec  20 2017: use fgets to reduce read-time 
  //
  int NROW = 0 ;
  int i, ivar, isn, ICAST, nptr;

  char ctmp[MXCHAR_FILENAME], LINE[MXCHAR_LINE], *ptrtok, cvar[100];
  char KEYNAME_ID[40];
  long double DVAR[MXVAR_TABLE];
  char        CVAR[MXVAR_TABLE][60];
  
  char fnam[]    = "SNTABLE_READ_EXEC_TEXT" ;
  int  NVAR_TOT  = READTABLE_POINTERS.NVAR_TOT ;  // all variables
  int  NVAR_READ = READTABLE_POINTERS.NVAR_READ ; // subset to read
  FILE *FP       = PTRFILE_TEXT ; 
  
  // ------------ BEGIN -----------  
   
  // get key name of ID varname such as CID, GALID, etc.
  sprintf(KEYNAME_ID,"%s", READTABLE_POINTERS.VARNAME[0] ); 

  while ( fgets(LINE, MXCHAR_LINE, FP ) != NULL ) {

    // check first word in the line
    ptrtok = strtok(LINE, " ");
    sprintf(ctmp, "%s", ptrtok);
    if ( ctmp[0] == '#' ) { continue ; }  // skip comment lines
	
#ifdef TEXTFILE_NVAR
    // for catenated TEXT files, NVAR key can appear
    // multiple times; ABORT if any NVAR is different
    // to avoid mistakes
    if ( strcmp(ctmp,"NVAR:") == 0 ) {
      ptrtok = strtok(NULL," " );
      sscanf(ptrtok, "%d", &NVAR_TMP ) ;
      NKEY_NVAR++ ;

      if ( NVAR_TMP != NVAR_TOT ) {
	sprintf(MSGERR1,"Invalid 'NVAR: %d' (NKEY=%d) -> column change.",
		NVAR_TMP, NKEY_NVAR);
	sprintf(MSGERR2,"Expect NVAR=%d.", NVAR_TOT);
	errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
      }
    }
#endif

    if ( validRowKey_TEXT(ctmp) == 0 ) { continue ; }

    // if we get here, we have a valid ROW key so read rest of row.

    NROW++ ;   

    ptrtok = strtok(NULL," " ); ivar=0 ;
    //  xxx mark delete while ( ptrtok != NULL && ptrtok != '\0' && ivar < NVAR_TOT) { 
    while ( ptrtok != NULL && ivar < NVAR_TOT) { 

      // Dec 20 2017: extract only variables on READ-list
      if ( READTABLE_POINTERS.NPTR[ivar] > 0 ) {      
	sprintf(cvar,"%s",   ptrtok );      
	sscanf(cvar, "%Lf",  &DVAR[ivar] );
	sscanf(cvar, "%s",   CVAR[ivar] );
      }
      ptrtok = strtok(NULL," " );
      ivar++ ;
    }
            
    if ( fmodf( (float)(NROW), 100000. ) == 0 && NROW > 0 )  { 
      printf("\t Reading table row %d  (%s=%s) \n", 
	     NROW, KEYNAME_ID, CVAR[0] );  fflush(stdout);
    }


    // set user arrays via pointer
    for ( i = 0; i < NVAR_READ; i++ ) {
      
      ivar  = READTABLE_POINTERS.PTRINDEX[i] ; // starts at 1
      ICAST = READTABLE_POINTERS.ICAST_STORE[ivar] ;     
      
	if ( ivar < 0 || ivar >= MXVAR_TABLE ) {
	  sprintf(MSGERR1,"Invalid PTRINDEX[%d] = %d", i, ivar );
	  sprintf(MSGERR2,"PTRINDEX must be %d to %d", 0, MXVAR_TABLE-1 );
	  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
	}
    
	isn = NROW - 1;  // isn is a C-like index

	for(nptr=0; nptr<READTABLE_POINTERS.NPTR[ivar]; nptr++ ) {

	  if ( ICAST == ICAST_D )  { 
	    READTABLE_POINTERS.PTRVAL_D[nptr][ivar][isn] = 
	      (double)DVAR[ivar] ; 
	  }	  
	  else if ( ICAST == ICAST_F )  { 
	    READTABLE_POINTERS.PTRVAL_F[nptr][ivar][isn] = 
	      (float)DVAR[ivar] ; 
	  }	  
	  else if ( ICAST == ICAST_I )  { 
	    READTABLE_POINTERS.PTRVAL_I[nptr][ivar][isn] = 
	      (int)DVAR[ivar] ; 
	  }	  
	  else if ( ICAST == ICAST_S )  { 
	    READTABLE_POINTERS.PTRVAL_S[nptr][ivar][isn] = 
	      (short int)DVAR[ivar] ; 
	  }	  
	  else if ( ICAST == ICAST_L )  { 
	    READTABLE_POINTERS.PTRVAL_L[nptr][ivar][isn] = 
	      (long long int)DVAR[ivar] ; 
	  }	  
	  else if ( ICAST == ICAST_C )  { 
	    sprintf(READTABLE_POINTERS.PTRVAL_C[nptr][ivar][isn],"%s",
		    CVAR[ivar]); 
	  }	  
	  else {
	    sprintf(MSGERR1,"Unknown ICAST=%d  var[%d]=%s  nptr=%d", 
		    ICAST, ivar, READTABLE_POINTERS.VARNAME[ivar], nptr );
	    sprintf(MSGERR2,"See ICAST_  parameters in sntools_output.h");
	    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
	  }
	} // end nptr loop
	
      } // end of i loop      

  } // end fscanf 
  
  fclose(FP);

  // May 2 2019: reset flags to allow opening another file.
  NAME_TABLEFILE[OPENFLAG_READ][IFILETYPE_TEXT][0] = 0 ;
  USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_TEXT]     = 0;

  return(NROW) ;

} // end of SNTABLE_READ_EXEC_TEXT

void SNTABLE_CLOSE_TEXT(void) {
  // May 2020
  // can call this function after SNTABLE_NEVT
  // so that there is no need to read entire file.
  NAME_TABLEFILE[OPENFLAG_READ][IFILETYPE_TEXT][0] = 0 ;
  USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_TEXT]     = 0;
} 

// =========================================
int validRowKey_TEXT(char *string) {

  // Aug 28 2013
  //  Returns 1 if string is a valid row key in a 
  //  text-formatted table file.
  //  Return 0 otherwise.

  if ( strcmp(string,"SN:")    == 0 ) { return 1; }
  if ( strcmp(string,"ROW:")   == 0 ) { return 1; }
  if ( strcmp(string,"GAL:")   == 0 ) { return 1; }
  if ( strcmp(string,"OBS:")   == 0 ) { return 1; } // May 30 2020 (SPECTRA)
  //  if ( strcmp(string,"LIBID:") == 0 ) { return 1; }

  return 0 ;
} // validRowKey_TEXT

// ====================================
int ICAST_for_textVar(char *varName) {

  // Oct 2014
  // Since ascii files do not include a cast, this function
  // has hard-coded string values to return ICAST_C for
  // string variables such as CCID and FIELD. All other 
  // variables are returned as float (ICAST_F).
  //
  // Note that users can explicitly cast a string by appending 
  // :C to the variables passed to SNTABLE_READPREP_VARDEF(..)
  // This ICAST function allows sloppy read-usage when the 
  // explicit cast is left out of SNTABLE_READPREP_VARDEF(..).
  //
  // Nov 02 2014: use strcmp_ignoreCase instead of strcmp.
  // Feb 17 2016: remove GALID from ICAST_C list
  // Dec 04 2016: add VERSION to ICAST_C list
  // Aug 04 2017: allow for "CCID ROW" using strstr
  // Apr 29 2019: allow PARNAME

  // 1st-key identifiers  
  if ( strcmp_ignoreCase(varName,(char*)"CID"  )  == 0 ) 
    { return ICAST_C; }
  if ( strcmp_ignoreCase(varName,(char*)"SNID" )  == 0 ) 
    { return ICAST_C; }
  if ( strcmp_ignoreCase(varName,(char*)"CCID" )  == 0 ) 
    { return ICAST_C; }
  if ( strcmp_ignoreCase(varName,(char*)"GALID")  == 0 ) 
    { return ICAST_C; }
  if ( strcmp_ignoreCase(varName,(char*)"ROW"  )  == 0 ) 
    { return ICAST_C; }
  if ( strcmp_ignoreCase(varName,(char*)"STARID") == 0 ) 
    { return ICAST_C; }

  // allow for things like "CCID ROW"
  if ( strstr(varName,"CCID") != NULL ) { return ICAST_C; }
  if ( strstr(varName,"ROW" ) != NULL ) { return ICAST_C; }

  // misc string-keys
  if ( strcmp_ignoreCase(varName,(char*)"FIELD")  == 0 ) 
    { return ICAST_C; }

  if ( strcmp_ignoreCase(varName,(char*)"BAND" )  == 0 ) 
    { return ICAST_C; }

  if ( strcmp_ignoreCase(varName,(char*)"IAUC_NAME" ) == 0 ) 
    { return ICAST_C;}

  if ( strcmp_ignoreCase(varName,(char*)"CATALOG" )   == 0 ) 
    { return ICAST_C;}

  if ( strcmp_ignoreCase(varName,(char*)"VERSION" )   == 0 ) 
    { return ICAST_C;}

  if ( strcmp_ignoreCase(varName,(char*)"PARNAME" )   == 0 )  // 4.29.2019
    { return ICAST_C;}

  // if not a string, return float cast which really means
  // that it's not a char.

  return ICAST_F ;

} // end of ICAST_for_textVar


// ==========================================
void SNLCPAK_FILL_TEXT(void) {

  // Oct 15 2013: fix aweful bug; move close(ptr) inside ifilt loop.
  // Jan 22 2014: write optional rest-frame info (filter and flux)
  // Sep 07 2014: open 1 file per SN rather than per SN/filter

  int ifilt, USE, obs, FLAG, OUTFLAG, i ;
  int NOBS_DATA, NOBS_FITFUN, NOBS, REJECT, ISDATA, ISFIT ;
  int FLAGLIST[2] ;
  int NFIT, IFIT ;

  char *PREFIX, *CCID, BAND[2] ;
  char  fnam[] = "SNLCPAK_FILL_TEXT" ;

  // --------------- BEGIN ------------

  PREFIX = SNLCPAK_OUTPUT.TEXTFILE_PREFIX ;
  CCID   = SNLCPAK_OUTPUT.CCID ;
 
  NOBS_DATA   = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FLUXDATA];
  NOBS_FITFUN = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FITFUN];

  NFIT   = SNLCPAK_OUTPUT.NFIT_PER_SN ;
  IFIT   = SNLCPAK_OUTPUT.NLCPAK     ;

  // update list file
  fprintf(PTRFILE_LCLIST, "SN: %s  %d \n", CCID, IFIT );
  fflush(PTRFILE_LCLIST);


  // Aug 4 2015: write header on first call to make sure that
  //             all variables are defined.
  if ( NCALL_SNLCPAK_FILL == 1 ) 
    { SNLCPAK_WRITE_HEADER_TEXT(PTRFILE_LCPLOT);  }

  for(ifilt=0; ifilt < SNLCPAK_OUTPUT.NFILTDEF_SURVEY; ifilt++ ) {

    USE = SNLCPAK_OUTPUT.USEFILT[SNLCPAK_EPFLAG_FLUXDATA][ifilt] ;
    sprintf(BAND,"%c", SNLCPAK_OUTPUT.SURVEY_FILTERS[ifilt] );

    if ( USE == 0 ) { continue ; }

    FLAGLIST[0] = SNLCPAK_EPFLAG_FLUXDATA ;
    FLAGLIST[1] = SNLCPAK_EPFLAG_FITFUN ;

    for(i=0; i < 2; i++ ) {
      FLAG = FLAGLIST[i] ;
      NOBS = SNLCPAK_OUTPUT.NOBS[FLAG];

      ISDATA = ( FLAG == SNLCPAK_EPFLAG_FLUXDATA  ) ;
      ISFIT  = ( FLAG == SNLCPAK_EPFLAG_FITFUN    ) ;
      OUTFLAG = -9 ;  REJECT=0;

      for(obs=0; obs<NOBS; obs++ )  {  
    
	if ( ISDATA ) { OUTFLAG = 1; }
	if ( ISFIT  ) { OUTFLAG = 0; }
	
	if ( ISDATA ) {
	  REJECT = (int)SNLCPAK_OUTPUT.EPDATA[SNLCPAK_EPFLAG_REJECT][obs] ;
	  if ( REJECT ) { OUTFLAG = -1; }
	}
		
	if ( SNLCPAK_OUTPUT.IFILT[FLAG][obs] == ifilt ) { 
	  if ( OUTFLAG == -9 ) {
	    sprintf(MSGERR1,"Invalid OUTFLAG=%d for obs=%d", OUTFLAG, obs);
	    sprintf(MSGERR2,"ISDATA=%d  ISFIT=%d  REJECT=%d",
		    ISDATA, ISFIT, REJECT);
	    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
	  }
	  snlcpak_textLine(PTRFILE_LCPLOT,FLAG,obs,ifilt,OUTFLAG);  
	}

      } // obs   
    }  // FLAG loop

  }  // ifilt


}  // end of SNLCPAK_FILL_TEXT

// ==============================
void SNLCPAK_WRITE_HEADER_TEXT(FILE *fp) {
  
  // Created Sep 7 2014 by R.Kessler
  // write header info to LC ascii file, based on user-select format
  // in SNTABLE_LIST = 'LCPLOT(text:[fmt])' ; see OPT_FORMAT below.
  
  char fnam[] = "SNLCPAK_WRITE_HEADER_TEXT" ;
  int  OPT_FORMAT = SNLCPAK_OUTPUT.OPT_TEXT_FORMAT ;
  int  ivar, NVAR, NOBS_REST, NOBS_SIMREST ;

  // ----------- BEGIN ------------

  NOBS_REST    = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FLUXREST] ;
  NOBS_SIMREST = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_SIMFLUXREST] ;

  NVAR = NVAR_LCPLOT_REQ ;
  if ( NOBS_REST    > 0 ) { NVAR +=4 ; }
  if ( NOBS_SIMREST > 0 ) { NVAR +=1 ; }

  if ( OPT_FORMAT == OPT_FORMAT_COL ) {
    return ; // do nothing 
  }
  else if ( OPT_FORMAT == OPT_FORMAT_CSV ) {
    for(ivar=0; ivar < NVAR; ivar++ ) {
      fprintf(fp, "%s", VARNAME_SNLC[ivar] );
      if ( ivar < NVAR-1 ) { fprintf(fp,", "); }
    }
  }
  else if ( OPT_FORMAT == OPT_FORMAT_KEY ) {
#ifdef TEXTFILE_NVAR
    fprintf(fp, "NVAR:  %d\n", NVAR);
#endif
    fprintf(fp, "VARNAMES: ");
    for(ivar=0; ivar < NVAR; ivar++ ) 
      { fprintf(fp, "%s ", VARNAME_SNLC[ivar] ); }
  }
  else {
    sprintf(MSGERR1, "Invalid OPT_FORMAT = %d", OPT_FORMAT);
    sprintf(MSGERR2, "Check &SNLCINP variable SNTABLE_LIST");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  fprintf(fp, "\n");   fflush(fp);

} // end of   SNLCPAK_WRITE_HEADER_TEXT

// ==============================
void snlcpak_textLine(FILE *fp, int FLAG, int obs, int ifilt, int OUTFLAG) {

  // Inputs
  //  FLAG    : indicates data or fit, and used as index
  //  obs     : epoch/obs index
  //  ifilt   : sparse filter index
  //  OUTFLAG : 1->data, 0->best-fit, -1->rejected data
  //
  // Dec 19 2019: write SIM_FLUXCAL for sim
  // Jan 23 2020: write KCOR and AVWARP (for rest-frame model only)
  //
  int ISFIT       = ( FLAG == SNLCPAK_EPFLAG_FITFUN    ) ;
  int ISDATA      = ( FLAG == SNLCPAK_EPFLAG_FLUXDATA  ) ;
  int NOBS_FITFUN = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FITFUN];

  int IFILT_REST, NOBS_REST, NOBS_SIMREST, flag ;
  int OPT_FORMAT ;
  double chi2,  FLUX_REST, KCOR, AVWARP, SIM_FLUXCAL ;
  char BAND[2], BAND_REST[2], sep[4], comment[200], LINE[400], CVAL[80];

  char fnam[] = "snlcpak_textLine" ;

  // ------------ BEGIN --------------


  if ( NOBS_FITFUN == 0 || ISFIT ) 
    { chi2 = 0.0 ; }
  else
    { chi2 = SNLCPAK_OUTPUT.EPDATA[SNLCPAK_EPFLAG_CHI2][obs] ; }

  sprintf(BAND,"%c", SNLCPAK_OUTPUT.SURVEY_FILTERS[ifilt] );


  // get value-separate string for specified format.
  OPT_FORMAT = SNLCPAK_OUTPUT.OPT_TEXT_FORMAT ;
  sprintf(comment,"Called from %s", fnam);
  get_sepchar(OPT_FORMAT,comment, sep);


  if ( OPT_FORMAT == OPT_FORMAT_KEY ) 
    {  sprintf(LINE,"OBS: ");  }
  else
    {  LINE[0] = 0 ;  }

  strcat(LINE,SNLCPAK_OUTPUT.CCID);

  
  sprintf(CVAL,"%s %.3f",  sep, SNLCPAK_OUTPUT.MJD[FLAG][obs] );
  strcat(LINE,CVAL);

  sprintf(CVAL,"%s %7.2f", sep, SNLCPAK_OUTPUT.TOBS[FLAG][obs] );
  strcat(LINE,CVAL);

  
  if ( SNLCPAK_OUTPUT.SIMFLAG ) {
    int EPFLAG_SIM = SNLCPAK_EPFLAG_FLUXSIM ;
    SIM_FLUXCAL = 0.0 ;
    if ( ISDATA ) { SIM_FLUXCAL = SNLCPAK_OUTPUT.EPDATA[EPFLAG_SIM][obs]; }
    sprintf(CVAL,"%s %11.4le", sep, SIM_FLUXCAL );
    strcat(LINE,CVAL);    
  }


  sprintf(CVAL,"%s %11.4le", sep, SNLCPAK_OUTPUT.EPDATA[FLAG][obs] );
  strcat(LINE,CVAL);

  sprintf(CVAL,"%s %11.4le", sep, SNLCPAK_OUTPUT.EPDATA_ERR[FLAG][obs] );
  strcat(LINE,CVAL);

  sprintf(CVAL,"%s %2d",   sep, OUTFLAG );
  strcat(LINE,CVAL);

  sprintf(CVAL,"%s %s",    sep, BAND );
  strcat(LINE,CVAL);

  sprintf(CVAL,"%s %.2f",  sep, chi2 );
  strcat(LINE,CVAL);

  sprintf(CVAL,"%s %d",    sep, SNLCPAK_OUTPUT.NLCPAK );
  strcat(LINE,CVAL);

  /*
  fprintf(fp,"%-8s %.3f  %7.2f   %11.4le  %11.4le  %2d  %s  %.2f %d"
	  ,SNLCPAK_OUTPUT.CCID 
	  ,SNLCPAK_OUTPUT.MJD[FLAG][obs]
	  ,SNLCPAK_OUTPUT.TOBS[FLAG][obs]
	  ,SNLCPAK_OUTPUT.EPDATA[FLAG][obs]		
	  ,SNLCPAK_OUTPUT.EPDATA_ERR[FLAG][obs]		
	  ,OUTFLAG              // DATAFLAG		
	  ,BAND		
	  ,chi2
	  ,SNLCPAK_OUTPUT.NLCPAK   	  ); */

  // check optoin to write rest-frame flux for both data and sim

  NOBS_REST    = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FLUXREST] ;
  NOBS_SIMREST = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_SIMFLUXREST] ;

  if ( NOBS_REST > 0 ) {  // Jan 2014: rest-frame info

    // default is no rest-frame info
    FLUX_REST  = -999. ;
    KCOR       = -99. ;
    AVWARP     = -99. ;
    IFILT_REST = -9;
    sprintf(BAND_REST, "!" );
    
    if ( OUTFLAG == 1 )   {  
      flag        = SNLCPAK_EPFLAG_FLUXREST ;
      IFILT_REST  = SNLCPAK_OUTPUT.IFILTOBS[flag][obs]; 
    }

    if ( OUTFLAG == 1 && IFILT_REST > 0 ) {
      // rest info for each data point used in fit
      FLUX_REST  = SNLCPAK_OUTPUT.EPDATA[SNLCPAK_EPFLAG_FLUXREST][obs];
      KCOR   = SNLCPAK_OUTPUT.EPDATA[SNLCPAK_EPFLAG_KCOR][obs];
      AVWARP = SNLCPAK_OUTPUT.EPDATA[SNLCPAK_EPFLAG_AVWARP][obs];
      sprintf(BAND_REST, "%c", FILTERSTRING[IFILT_REST] );
    }

    //    sprintf(CVAL,"%s  %s%s %.5le",  sep, BAND_REST, sep, FLUX_REST);    
    sprintf(CVAL,"%s  %s%s %.5le%s %.3f%s %.3f",  
	    sep, BAND_REST, sep, FLUX_REST, sep,KCOR, sep, AVWARP);    
    strcat(LINE,CVAL);

  }  // NOBS_REST


  // and now check the sim rest-flux
  if ( NOBS_SIMREST > 0 ) { 
    flag   = SNLCPAK_EPFLAG_SIMFLUXREST ;
    
    if ( OUTFLAG == 1 ) {
      // rest info for each data point used in fit
      FLUX_REST  = SNLCPAK_OUTPUT.EPDATA[flag][obs];
      sprintf(CVAL,"%s %.5le", sep, FLUX_REST);
      strcat(LINE,CVAL);
      //      fprintf(fp," %.5le ", FLUX_REST);
    }
    else {
      // no info for fit-function
      FLUX_REST  = -999. ;
      sprintf(CVAL,"%s %.0f", sep, FLUX_REST);
      strcat(LINE,CVAL);
      //      fprintf(fp," %.0f ", FLUX_REST);
    } 
  }  // NOBS_REST

  fprintf(fp,"%s\n", LINE );
  fflush(fp);

  return ;

}  // end of snlcpak_textLine


// ==============================================
void OPEN_TEXTFILE_SPECLIST(char *PREFIX) {

  int  GZIPFLAG;
  char specFile[MXCHAR_FILENAME] ;
  char fnam[] = "OPEN_TEXTFILE_SPECLIST" ;

  // ------------- BEGIN -------------

  sprintf(specFile, "%s.%s", PREFIX, SUFFIX_SPECLIST_TEXT );
  PTRFILE_SPECLIST = open_TEXTgz(specFile,TEXTMODE_wt, &GZIPFLAG );
  if ( !PTRFILE_SPECLIST ) {
    sprintf(MSGERR1, "Could not open ascii SPECLIST file = ");
    sprintf(MSGERR2, "%s", specFile);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  sprintf(specFile, "%s.%s", PREFIX, SUFFIX_SPECPLOT_TEXT );
  PTRFILE_SPECPLOT = open_TEXTgz(specFile,TEXTMODE_wt, &GZIPFLAG );
  if ( !PTRFILE_SPECPLOT ) {
    sprintf(MSGERR1, "Could not open ascii SPECPLOT file = ");
    sprintf(MSGERR2, "%s", specFile);    
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  SPECPAK_USE_TEXT = true ;

  return ;

} // end OPEN_TEXTFILE_SPECLIST

void open_textfile_speclist__(char *PREFIX) {OPEN_TEXTFILE_SPECLIST(PREFIX);}

// ==============================
void SPECPAK_WRITE_HEADER_TEXT(void) {

  // April 2019
  // write VARNAMES header line to both SPECLIST and SPECPLOT file.
  // Beware that VARNAMES lists are hard-coded here.

  FILE *FP_LIST      = PTRFILE_SPECLIST ;
  char VARLIST_KEY[] = "CID ID NLAMBIN MJD TOBS TEXPOSE";
  char VARLIST_CSV[] = "CID,ID,NLAMBIN,MJD,TOBS,TEXPOSE";

  FILE *FP_PLOT      =  PTRFILE_SPECPLOT ;
  char VARPLOT_KEY[] = "CID ID LAMMIN LAMMAX FLAM FLAMERR";
  char VARPLOT_CSV[] = "CID,ID,LAMMIN,LAMMAX,FLAM,FLAMERR";

  int  OPT_FORMAT = SPECPAK_OUTPUT.OPT_TEXT_FORMAT ;

  char fnam[] = "SPECPAK_WRITE_HEADER_TEXT" ;

  // ----------- BEGIN ------------

  NVAR_SPECPLOT = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,VARPLOT_KEY);

  if ( OPT_FORMAT == OPT_FORMAT_COL ) {
    return ; // do nothing 
  }
  else if ( OPT_FORMAT == OPT_FORMAT_CSV ) {
    fprintf(FP_LIST, "%s", VARLIST_CSV );
    fprintf(FP_PLOT, "%s", VARPLOT_CSV );
  }
  else if ( OPT_FORMAT == OPT_FORMAT_KEY ) {
    fprintf(FP_LIST, "VARNAMES: %s", VARLIST_KEY);
    fprintf(FP_PLOT, "VARNAMES: %s", VARPLOT_KEY);
  }
  else {
    sprintf(MSGERR1, "Invalid OPT_FORMAT = %d", OPT_FORMAT);
    sprintf(MSGERR2, "Check &SNLCINP variable SNTABLE_LIST");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  fprintf(FP_LIST, "\n");   fflush(FP_LIST);
  fprintf(FP_PLOT, "\n");   fflush(FP_PLOT);

} // end of   SNLCPAK_WRITE_HEADER_TEXT

// ==========================================
void SPECPAK_FILL_TEXT(void) {

  // Apr 2019
  // Write spectra data
  // If CID is of the form nnn[z=xxx], remove stuff in brackets
  // since [z=xxx] is only for marz table.

  int   NLAMBIN_TOT = SPECPAK_OUTPUT.NLAMBIN_TOT;
  int   NSPEC       = SPECPAK_OUTPUT.NSPEC;
  int   ispec, ilam, LENCCID;
  char  CCID[MXCHAR_CCID], *ctmp;
  //  char  fnam[] = "SPECPAK_FILL_TEXT" ;

  // --------------- BEGIN ------------

  if ( NCALL_SPECPAK_FILL == 1 ) 
    { SPECPAK_WRITE_HEADER_TEXT();  }
 
  sprintf(CCID,"%s", SPECPAK_OUTPUT.CCID);
  if ( SPECPAK_USE_MARZ ) {
    ctmp = strchr(CCID,'[') ;
    if ( ctmp != NULL ) {  LENCCID = (int)(ctmp-CCID); CCID[LENCCID] = 0 ; }
  }

  // fill SPEC-LIST table vs. ID
  for(ispec=0; ispec < NSPEC; ispec++ ) {

    fprintf(PTRFILE_SPECLIST,"OBS: %8s  %2d  %4d  %.3f %6.1f  %.1f \n",
	    CCID,
	    SPECPAK_OUTPUT.ID_LIST[ispec],
	    SPECPAK_OUTPUT.NLAMBIN_LIST[ispec],
	    SPECPAK_OUTPUT.MJD_LIST[ispec],
	    SPECPAK_OUTPUT.TOBS_LIST[ispec],
	    SPECPAK_OUTPUT.TEXPOSE_LIST[ispec] );	    
    fflush(PTRFILE_SPECLIST);
  }

  // fill SPEC-PLOT table vs. wavelength
  for(ilam=0; ilam < NLAMBIN_TOT; ilam++ )
    { specpak_textLine(PTRFILE_SPECPLOT, CCID, ilam, 0); }


}  // end of SNLCPAK_FILL_TEXT

// ==============================
void specpak_textLine(FILE *fp, char *CCID, int ilam, int OUTFLAG ) {

  // Inputs
  //  *fp      : file pointer to write to
  //  ilam     : wave bin
  //  OUTFLAG  : not used

  int    ID       = SPECPAK_OUTPUT.ID[ilam];
  double LAMMIN   = SPECPAK_OUTPUT.LAMMIN[ilam];
  double LAMMAX   = SPECPAK_OUTPUT.LAMMAX[ilam];
  double FLAM     = SPECPAK_OUTPUT.FLAM[ilam];
  double FLAMERR  = SPECPAK_OUTPUT.FLAMERR[ilam];
  int  OPT_FORMAT = SPECPAK_OUTPUT.OPT_TEXT_FORMAT ;
  char LINE[400]  = "" ;
  char sep[2]     = "" ;
  char CVAL[80];
  char fnam[]     = "specpak_textLine" ;

  // ------------ BEGIN --------------

  get_sepchar(OPT_FORMAT, fnam, sep); // return sep (blank or comma)

  if ( OPT_FORMAT == OPT_FORMAT_KEY ) 
    {  sprintf(LINE,"OBS: ");  }
  else
    {  LINE[0] = 0 ;  }

  sprintf(CVAL,"%8s",             CCID    ) ; strcat(LINE,CVAL);
  sprintf(CVAL,"%s %d",      sep, ID      ) ; strcat(LINE,CVAL);
  sprintf(CVAL,"%s %7.1f",   sep, LAMMIN  ) ; strcat(LINE,CVAL);
  sprintf(CVAL,"%s %7.1f",   sep, LAMMAX  ) ; strcat(LINE,CVAL);
  sprintf(CVAL,"%s %10.3le", sep, FLAM    ) ; strcat(LINE,CVAL);
  sprintf(CVAL,"%s %10.3le", sep, FLAMERR ) ; strcat(LINE,CVAL);

  fprintf(fp,"%s\n", LINE );
  fflush(fp);

}  // end of specpak_textLine

