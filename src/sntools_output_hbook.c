/***************************************************
 sntools_output_hbook.c
 Created Feb 23 2013 by R.Kessler
 Contains all hbook-specific code for snana.
 Included in sntools_output.c if HBOOK flag is set.


 
                 HISTORY

 Apr 01 2013: in SNTABLE_ADDCOL_HBOOK for ntuples, allow
               c-like ':D' and ':F'

 Apr 02, 2013: HLIM -> 2 million (was 1 million) to handle
               larger MXLC_PLOT.

 Apr 15 2013: in SNTABLE_ADDCOL_HBOOK(), add (char*) cast in front of malloc.
              Fixed c++ compile bug found by STK 

 Apr 20 2013: Add SNTABLE_DUMP_* functions.
              Use hbnamc for char.

 Apr 25 2013: add SNTABLE_LIST_HBOOK() -> just abort.

 Jun 17, 2013
   - new global HBOOK_STATUS for NEW/OLD;
     affects CLOSE_HBOOK function

   - new functions SNHIST_RD* 

 Feb 8 2014: use tolower() function to make case-insensitive
             comparisons in get_cwnt_ivar().

 Apr 24 2014:
  New utility sntable_fetch_HBOOK(...) called by
  SNTABLE_DUMP_VALUES_HBOOK  and SNTABLE_READ_VALUES_HBOOK
  to dump or to read/pass values .

    
 May 11 2014
   SNTABLE_LOAD_HBOOK -> SNTABLE_ADDCOL_HBOOK  and major overhaul.

 Dec 2014: reduce hbook_mem from 1.5M to 1.0M

 Jan 2017:
   +  add function get_TABLEID_HBOOK()
   +  move ISFILE_HBOOK to stools_output.c

 Aug 11 2017:
   + convert HBOOK_CWNT_READROW.VAL_[CAST] into pointers
     that are allocated for each cast. Reduces program
     memory by 38 MB, and allocated memory is way smaller
     than previous static [38MB] memory.
     See new function  malloc_READROW_HBOOK().

 May 4 2019: 
   + HBOOK_MEM -> 500,000 (was 1 million). Needed for Git & Midway.

***************************************************/


#ifdef __cplusplus
extern"C" {
#endif

  void OPEN_HFILE(char *FILENAME, int LUN, 
		  char *COPT, char *TOPDIR, int *IERR) ;
  void CLOSE_HFILE(char *FILENAME, int OPENFLAG);

  void SNTABLE_CREATE_HBOOK(int IDTABLE, char *name) ;
  void SNTABLE_FILL_HBOOK(int IDTABLE ) ;
  void SNTABLE_ADDCOL_HBOOK(int IDTABLE, char *BLOCK, void* PTRVAR,
			    SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF) ;

  int  SNTABLE_NEVT_HBOOK(char *FILENAME, int NTID) ;
  void SNTABLE_LIST_HBOOK(char *file) ;
  void SNTABLE_DUMP_VARNAMES_HBOOK(char *fileName, int NTID);

  int SNTABLE_READPREP_HBOOK(int NTID);
  int SNTABLE_READ_EXEC_HBOOK(void);

  void malloc_READROW_HBOOK(int OPT) ;
  void sntable_pushRowOut_hbook(int EVT, int OPT, int IFIT);
  void outlier_dump_fromHbook(int irow) ;

  void removeParenth(char *string) ;

  void SNHIST_INIT_HBOOK(int NDIM, int ID, char *TITLE, 
			 int *NBIN, double *XMIN, double *XMAX ) ;
  void SNHIST_FILL_HBOOK(int NDIM, int ID, double *VALUE, double WGT) ;

  void SNHIST_RDBINS_HBOOK(int NDIM, int ID, char *CTIT, int *NBIN,
			   double *XMIN, double *XMAX);
  void SNHIST_RDCONT_HBOOK(int NDIM, int ID, int NBIN, double *CONTENTS);

  void MAKEDIR_HBOOK(char *CCID, int CID) ;
  void CDTOPDIR_HBOOK(void);
  
  void SNLCPAK_SURVEY_HBOOK(void);

  void modify_text_4paw(char *TEXT) ;
  void SNLCPAK_FILL_HBOOK(void) ;

  void get_cwnt_info (char *fileName, int NTID, int OPT) ;
  int  get_cwnt_ivar (char *varName) ;
  int  parse_cwnt_vardef(char *VARCAST_NAME, char *VARNAME, char *CAST ) ;
  int  TABLEID_HBOOK(char *tableName); // Jan 2017

  // define HBOOK functions
  void hlimit_(int *hlim);
  void hropen_(int *LUN, char *TOPDIR, char *FILENAME, char *COPT_CAP, 
	       int *LREC, int *IERR, int LEN1, int LEN2, int LNE3);
  void hcdir_(char *dir, char *COPT, int LEN1, int LEN2);
  void hmdir_( char *subDirName, char *COPT, int LEN1, int LEN2) ;
  
  void hrout_(int *HID, int *ICYC, char *HOPT, int LEN1 );
  void hrend_(char *TOPDIR, int LEN1 ) ;
  
  void hbook1_(int *IDLOC, char *TITLE, int *NBIN, float *XMIN, float *XMAX, 
	       float *w, int LENTIT ); 
  
  void hbook2_(int *IDLOC, char *TITLE, 
	       int *NXBIN, float *XMIN, float *XMAX, 
	       int *NYBIN, float *YMIN, float *YMAX, 
	       float *w, int LENTIT ); 

  void hf1_ ( int *IDLOC, float *VAL, float *WGT) ;
  void hf2_ ( int *IDLOC, float *XVAL, float *YVAL, float *WGT) ;

  void hbnt_(int *ID, char *name, char *blank, int LEN1, int LEN2);
  void hfnt_(int *ID);
  void hbname_(int *ID, char *BLOCK, void *PTRVAR, char *VARLIST, 
	       int LEN1, int LEN2);

  void hbnamc_(int *ID, char *BLOCK, char *PTRVAR, char *VARLIST, 
	       int LENBLK, int LENPTR, int LENVAR);

  void hpak_ ( int *HID, float *DATA     ) ;
  void hpake_( int *HID, float *DATA_ERR ) ;
  
  void hrin_(int *NTID, int *ICYC, int *IOFF);
  void hgiven_( int *NTID, char *CHTIT, int *NVAR, 
		char **CTAG, float *RLOW, float *RHI, int LEN1, int LEN2);

  void hgive_(int *HID, char *CTIT, int *NBX, float *X0, float *X1, 
	      int *NBY, float *Y0, float *Y1, int *NWT, int *LOC, int LENTIT);

  void hunpak_(int *HID, float *farray, char *choice, int *NUM, int LENCHOICE);

  void hnoent_(int *NTID, int *NEVENT);

  void hntvdef_(int *NTID, int *IVAR, char *CHTAG, char *BLOCK, int *ITYPE,
	       int LENTAG, int LENBLOCK );

  void hgnt_(int *NTID, int *IDEVT, int *IERR);

  void hmerge_wrapper_(int *OPT, char *fileName, int LENF);

  // fortran wrapper to close LUN (in snana.car)
  void forclose_(int *LUN);



#ifdef __cplusplus
}
#endif


// globals
#define HBOOK_MEM  500000
struct pawC  { float hmem[HBOOK_MEM];  } pawc_ ;
//struct Quest { int   iquest[100]; } quest_ ;

const char *HBOOK_TOPDIRNAM_NEW  = "snana" ;
const char *HBOOK_TOPDIRNAM_READ = "read"  ;

char HBOOK_TOPDIR_NEWFILE[20] ;  // store topdir only for NEW file
char HBOOK_TOPDIR_READFILE[20] ; // store topdir for existing hbook file

int LUNHBOOK_NEWFILE ;
int LUNHBOOK_READFILE ;

// define structure for CWNT info read from file 
// (for SNTABLE_DUMP_*)
#define MXVAR_CWNT MXVAR_TABLE  // max number of CWNT variables
#define CBUFF_CWNT 32           // max size of char tag
struct  HBOOK_CWNT_INFO {
  int NVAR ;

  // ivar indics are global over all variables
  char *VARNAME[MXVAR_CWNT] ;  // final list of variable names
  char *BLOCK[MXVAR_CWNT] ;    // block name for each variable
  char *CAST[MXVAR_CWNT] ;
  int  ICAST[MXVAR_CWNT] ;     // integer cast of each variable, (D,F,I,C)  

  int NEVENT ; // number of rows/entries

} HBOOK_CWNT_INFO ;


struct HBOOK_CWNT_READROW {
  // address for hbname: 
  // Aug 2 2014: scalar -> vector[MXEPOCH]

  char   CCID[40] ;

  // define sparse arrays to map table column into sparse ivar for each cast
  int NVARCAST[MXCAST];    // Nvar for each cast
  int IVARCAST_MAP[MXVAR_CWNT][MXCAST];  // sparse ivar vs. [IVAR(table)][ICAST]

  // all ivar indices are sparse over the dump variables
  int    ICAST[MXVAR_CWNT] ;
  double **VAL_D ; // to be malloced vs. sparse NVARCAST and epoch.
  float  **VAL_F ;
  int    **VAL_I ;
  char  ***VAL_C ;

} HBOOK_CWNT_READROW ;


#define MXL_HBOOK 10 // max number of long long ints to store
struct HBOOK_L {
  int    NSTORE;         // number of long long int stored as R*8
  double    DUMMYSTORE_D[MXL_HBOOK]; // store in CWNT instead of long int
  long long int *DUMMYSTORE_L[MXL_HBOOK]; 
} HBOOK_L ;


// ======================================================
//
//               OPEN/CLOSE FUNCTIONS
//
// ======================================================


// --------------------------------------------------
void OPEN_HFILE(char *FILENAME, int LUN, 
		char *COPT, char *TOPDIR, int *RETERR) {

  // Feb 2013: wrapper to initialize file for hbook.
  //
  // Open hbook file, and call hilimit on first call.
  // (I) FILENAME = name of file.
  // (I) LUN      = log unit number
  // (I) COPT     = 'N' for new file, ' ' for existing file
  // (I) TOPDIR   = name of top directory
  // (O) IERR     : returns status, 0=OK
  //
  // Apr 25 2014: if 'Q' is part of COPT, then this is QUIET mode
  //              --> print nothing to screen.
  //

  int 
    hlim 
    ,LREC = 1024
    ,IERR, LENOPT, i
    ,LVBOSE, LNEW, FIRST
    ;

  char COPT_CAP[8], c1[2];

  // ----------------- BEGIN --------------

  // include 'P' in option to allow capital letters in FILENAME
  // If Q or q is specified, set QUIET flag and remove from 
  // COPT_CAP that is passed.
  LVBOSE = 1 ;  LNEW=0 ;
  LENOPT = strlen(COPT);
  sprintf(COPT_CAP, "P") ;  
  for(i=0; i < LENOPT; i++ ) {
    sprintf(c1,"%c", COPT[i] );
    if ( strcmp(c1,"Q") == 0 ) { LVBOSE = 0; continue ; }
    if ( strcmp(c1,"q") == 0 ) { LVBOSE = 0; continue ; }
    if ( strcmp(c1,"N") == 0 ) { LNEW   = 1;  }
    sprintf(COPT_CAP, "%s%s ", COPT_CAP, c1 ) ;  
  }
  sprintf(COPT_CAP, "%s" , COPT_CAP) ;  

  // - - - - - - - - - - - - - - -
  

  if ( LVBOSE ) 
    { print_banner("OPEN_HFILE: Initialize HBOOK file" ) ; }


  FIRST = FIRSTCALL_TABLEFILE_OPEN[IFILETYPE_HBOOK] ;
  if ( FIRST ) {
    hlim = HBOOK_MEM ;
    if ( LVBOSE ) 
      { printf("\t call hlimit(%d) ... \n", hlim);  fflush(stdout); }
    hlimit_(&hlim);
    FIRSTCALL_TABLEFILE_OPEN[IFILETYPE_HBOOK] = 0;
  }

  if ( LNEW ) { SNLCPAK_USE_HBOOK = 1 ; }

  hropen_(&LUN, TOPDIR, FILENAME, COPT_CAP, &LREC, &IERR,
	  strlen(TOPDIR), strlen(FILENAME), strlen(COPT_CAP) );

  *RETERR = IERR ; // load output function arg

  if ( IERR == 0 ) {

    char tmpDir[100] = "  " ;
    char R[4]        = "R "; 
    char cstat[8];

    // read full path of current hbook directory
    hcdir_(tmpDir, R, strlen(TOPDIR)+3, strlen(R) ) ;

    if ( strchr(COPT,'N') != NULL ) { 
      sprintf(cstat,"NEW");
      sprintf(HBOOK_TOPDIR_NEWFILE,"%s", tmpDir);
      LUNHBOOK_NEWFILE = LUN ;
    }
    else { 
      sprintf(cstat,"OLD"); 
      sprintf(HBOOK_TOPDIR_READFILE,"%s", tmpDir);
      LUNHBOOK_READFILE = LUN ;
    }

    if ( LVBOSE ) {
      printf("   Opened %s HBOOK file\n\t %s \n", cstat, FILENAME);
      printf("\t Current directory: '%s' (IERR=%d  LREC=%d) \n", 
	     tmpDir, IERR, LREC );
      fflush(stdout);
    }

  }  // end if IERR==0 block


  HBOOK_L.NSTORE = 0;

} // end of OPEN_HFILE



// ===================================================
void CLOSE_HFILE(char *FILENAME, int OPENFLAG) {   
  
  // Input FILENAME is used only for print statement.
  // Input OPENFLAG (NEW or READ) determines how file is closed.

  int  HID    = 0 ;
  int  ICYC   = 1 ;
  int  LUN ;
  char HOPT[] = "T ";
  char TOPDIR[40];

  if ( OPENFLAG == OPENFLAG_NEW )   {  
    hrout_(&HID, &ICYC, HOPT, strlen(HOPT) ); 
    sprintf(TOPDIR, "%s", HBOOK_TOPDIR_NEWFILE);
    LUN = LUNHBOOK_NEWFILE ;
  }
  else {
    sprintf(TOPDIR, "%s", HBOOK_TOPDIR_READFILE);
    LUN = LUNHBOOK_READFILE ;
  }

  hrend_(TOPDIR, strlen(TOPDIR) ) ;

  // 4/27/2014: optional fortran wrapper to close LUN as 
  // recommended in hbook manual. This is not really needed,
  // but to include it requires updating $OBJ_OUTPUT
  // defined in the Makefile.
  //  forclose_(&LUN) ; 

  printf("   Close hbook file %s \n", FILENAME ); 
  fflush(stdout);

} // end of CLOSE_HFILE


// ======================================================
//
//               SNTABLE FUNCTIONS
//
// ======================================================

// --------------------------------------------------
void SNTABLE_CREATE_HBOOK(int IDTABLE, char *name) {
  int ID; 
  char blank[] = " " ;
  char fnam[] = "SNTABLE_CREATE_HBOOK" ;
  // -------- BEGIN ----------
  ID = IDTABLE ; // local ID has a real address
  printf("\n  %s: Book %s colunn-wise ntuple (NTID=%d) \n",
	 fnam, name, ID) ;
  hbnt_(&ID, name, blank, strlen(name), strlen(blank) ); 
  fflush(stdout);
}


// --------------------------------------------------
void SNTABLE_FILL_HBOOK(int IDTABLE ) {


  // Mar 2019: crazy hack to store long int as REAL*8
  int i;
  for(i=0; i < HBOOK_L.NSTORE; i++ ) 
    { HBOOK_L.DUMMYSTORE_D[i] = *HBOOK_L.DUMMYSTORE_L[i]; }


  // fill CWNT
  int ID = IDTABLE ;  hfnt_(&ID); 
}


// --------------------------------------------------
void SNTABLE_ADDCOL_HBOOK(int IDTABLE, char *BLOCK, void* PTRVAR, 
			  SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF) {
  
  // Add VARDEF column(s) to hbook ntuple. 
  //
  // May 10 2014: major change to use input struct ADDCOL_VARDEF

  int ID, NVAR, LENLIST, LENV, LENB, ivar, ICAST, IVEC, ISIZE, NLL ;
  char *VARLIST_FINAL, *VARLIST_ORIG, VARNAME[60], CCAST[4], tmpVar[60] ;
  char fnam[] = "SNTABLE_ADDCOL_HBOOK" ;

  // ---------------- BEGIN --------------
  ID   = IDTABLE ;

  LENLIST       = strlen(ADDCOL_VARDEF->VARLIST_ORIG);
  VARLIST_FINAL = (char*)malloc(sizeof(char) * (LENLIST+100) );
  VARLIST_ORIG  = (char*)malloc(sizeof(char) * (LENLIST+100) );


  // make comma-separated list
  NVAR = ADDCOL_VARDEF->NVAR ; 
  ICAST = IVEC = -9 ;

  for ( ivar=0; ivar < NVAR ; ivar++ ) {
    ICAST = ADDCOL_VARDEF->ICAST[ivar] ;
    IVEC  = ADDCOL_VARDEF->VECTOR_FLAG[ivar] ;
    ISIZE = ADDCOL_VARDEF->ISIZE[ivar] ;
    sprintf(CCAST,   "%s", ADDCOL_VARDEF->CCAST[ivar]) ;
    sprintf(VARNAME, "%s", ADDCOL_VARDEF->VARNAME[ivar]) ;

    // replace C-like cast with fortran cast for HBNAME
    if( strcmp(CCAST,"D")  == 0 ) { sprintf(CCAST,"R*8" ) ; }
    if( strcmp(CCAST,"F")  == 0 ) { sprintf(CCAST,"R"   ) ; }
    if( strcmp(CCAST,"I")  == 0 ) { sprintf(CCAST,"I"   ) ; }
    if( strcmp(CCAST,"L")  == 0 ) { sprintf(CCAST,"R*8" ) ; } //long int->R*8 

    // construct comma-separated list for HBNAME
    if ( IVEC == 1 ) {
      sprintf(tmpVar, "%s[0,%d]:%s", VARNAME, ISIZE, CCAST); 
    }
    else if ( ICAST == ICAST_C ) {
      sprintf ( tmpVar, "%s:%s*%d", VARNAME, CCAST, ISIZE ); 
    }
    else {
      sprintf(tmpVar, "%s:%s", VARNAME, CCAST ); 
    }

    if ( ivar == 0 ) 
      { sprintf(VARLIST_FINAL,"%s", tmpVar ); }
    else {
      strcat(VARLIST_FINAL,","); 
      strcat(VARLIST_FINAL,tmpVar); 
    }


  } // end ivar loop


  // ---------------------------------

  LENB = strlen(BLOCK) ;
  LENV = strlen(VARLIST_FINAL) ;

  if ( IVEC == 999 ) {
    printf(" xxx ---------------------------y--- \n");
    printf(" xxx VARLIST(ICAST=%d) = '%s' ->  '%s'  \n", 
	   ICAST, ADDCOL_VARDEF->VARLIST_ORIG, VARLIST_FINAL);
    fflush(stdout);
  }

  if ( ICAST == ICAST_C ) { 
    hbnamc_(&ID, BLOCK, (char*)PTRVAR, VARLIST_FINAL, 
	    LENB, sizeof(char*), LENV);  
  }
  else if ( ICAST == ICAST_L ) {
    // store double DUMMYSTORE instead of long int
    NLL = HBOOK_L.NSTORE ;
    HBOOK_L.DUMMYSTORE_L[NLL] = (long long int*)PTRVAR;
    hbname_(&ID, BLOCK, &HBOOK_L.DUMMYSTORE_D[NLL],VARLIST_FINAL,LENB,LENV);
    HBOOK_L.NSTORE++ ;
    if ( HBOOK_L.NSTORE >= MXL_HBOOK ) {
      sprintf(MSGERR1, "Storing too many long long ints");
      sprintf(MSGERR2, "NSTORE(L)=%d exceeds MXL_HBOOK=%d",
	      HBOOK_L.NSTORE, MXL_HBOOK );
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }
  }
  else  { 
    // D,F,I
    hbname_(&ID, BLOCK, PTRVAR, VARLIST_FINAL, LENB, LENV);  
  }

  free(VARLIST_ORIG);
  free(VARLIST_FINAL);

} // end of SNTABLE_ADDCOL_HBOOK


// ======================================
void SNTABLE_LIST_HBOOK(char *file) {

  // Apr 25 2013
  // just abort if we get here.
  // Analogous root function will list trees.

  char *s[MXTABLEFILETYPE] ;
  char fnam[] = "SNTABLE_LIST_HBOOK" ;
 
  // --------------- BEGIN ---------

  s[1] =  STRING_IDTABLE_SNANA[IFILETYPE_HBOOK] ;
  s[2] =  STRING_IDTABLE_FITRES[IFILETYPE_HBOOK] ;
  s[3] =  STRING_IDTABLE_OUTLIER[IFILETYPE_HBOOK] ; // Mar 2021

  sprintf(MSGERR1, "Cannot list ntuples.");
  sprintf(MSGERR2, "Try %s or %s or %s", s[1], s[2], s[3] );
  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 

  return ;

} // end of SNTABLE_LIST_HBOOK


// ====================================
int SNTABLE_NEVT_HBOOK(char *FILENAME, int NTID) {
  // Apr 2014
  // return Number of events in this ntuple
  get_cwnt_info(FILENAME, NTID, -1); // open/read/close
  return  HBOOK_CWNT_INFO.NEVENT ;

} // end of SNTABLE_NEVT_HBOOK


// ===============================================
int SNTABLE_READPREP_HBOOK(int NTID) {

  // Created Oct 14 2014
  // Read list of variables in table and store in 
  // struct READTABLE_POINTERS

  int NVAR, ivar ;
  char noFile[] = "" ;
  char tmpName[60];
  //  char fnam[] = "SNTABLE_READPREP_HBOOK" ;

  // --------------- BEGIN ---------------

  // read variable list and load into HBOOK_CWNT_INFO struct.
  get_cwnt_info(noFile, NTID, +1 ); // file already open, leave open

  NVAR = HBOOK_CWNT_INFO.NVAR ;
  
  // transfer info to more global READTABLE_POINTERS struct.
  for ( ivar=0; ivar < NVAR; ivar++ )  { 

    // remove brackets [] or ()  from varName
    sprintf(tmpName,"%s", HBOOK_CWNT_INFO.VARNAME[ivar]  );
    removeParenth(tmpName); // Aug 2 2014

    sprintf(READTABLE_POINTERS.VARNAME[ivar], "%s", tmpName);
    READTABLE_POINTERS.ICAST_READ[ivar] = HBOOK_CWNT_INFO.ICAST[ivar] ;
  }
  
  return NVAR ;

} // end of SNTABLE_READPREP_HBOOK


// ======================================
int SNTABLE_READ_EXEC_HBOOK(void) {

  // Created Oct 2014:
  // Driver function to read hbook ntuple rows and do one
  // of the following:
  // * fill user-defined pointers from calls to SNTABLE_READPREP_VARDEF.
  // * write each row to ascii file
  // * write outlier epochs to ascii file.
  //

  int NEVT = 0 ;
  int NLL  = 0 ;
  int IVAR_L[10];
  char fnam[] = "SNTABLE_READ_EXEC_HBOOK" ;

  int  IVAR_READ, NVAR_READ, IVAR_TOT, IVAR_CWNT, ep ;
  int  ICAST_CWNT, OPT_READ, ivarcast ;
  int  LENB, LENC, LENV;
  char *varName, *blkName, *ptrC[MXEPOCH] ;
  void *ptrTmp ;

  // hbname args
  char CNULL[]         = " " ;
  int  INULL           = 0   ;
  char VARNAME_CLEAR[] = "$CLEAR" ;
  char VARSET[80] ;

  // ------------ BEGIN -----------

  // check OPT input to see if we are dumping to text,
  // or storing to ARRARY_OUT.

  OPT_READ = 0 ;

  if ( READTABLE_POINTERS.FP_DUMP ) { 
    OPT_READ = OPT_SNTABLE_READ_forDUMP ;
    printf("\t Dump SNTABLE values for each event. \n");
  }
  else if ( READTABLE_POINTERS.MXLEN > 1 ) { 
    OPT_READ = OPT_SNTABLE_READ_forARRAY ;
    printf("\t Read/store values for each CWNT row. \n" );
  } 
  else if ( OUTLIER_INFO.USEFLAG ) {
    OPT_READ = OPT_SNTABLE_READ_forOUTLIERS ; // check OUTLIER option 
  }
  else {
    sprintf(MSGERR1,"Cannot figure out what to do");
    sprintf(MSGERR2,"Check calls to SNTABLE_READ_XXX");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);     
  }

  int NTID = atoi(READTABLE_POINTERS.TABLENAME);

  // define globals to read back variables
  
  hbname_(&NTID, CNULL,  &INULL,  VARNAME_CLEAR, 
	  strlen(CNULL), strlen(VARNAME_CLEAR) );

  // for each variable to read or dump, find IVAR_CWNT so that
  // we can get the BLOCK and CAST needed to read it back

  NVAR_READ = READTABLE_POINTERS.NVAR_READ ;
  ptrTmp =  fnam; // to avoid compile warning with -Wall

  // ------------------------------------------------------
  // determine number of variabes for each cast,
  // then allocate memory for each cast.
  malloc_READROW_HBOOK(+1);

  // ------------------------------------------------------

  for ( IVAR_READ=0; IVAR_READ < NVAR_READ; IVAR_READ++ ) {
    
    IVAR_TOT     = READTABLE_POINTERS.PTRINDEX[IVAR_READ] ;
    varName      = READTABLE_POINTERS.VARNAME[IVAR_TOT] ;
    IVAR_CWNT    = get_cwnt_ivar(varName);

    ICAST_CWNT   = HBOOK_CWNT_INFO.ICAST[IVAR_CWNT] ;
    blkName      = HBOOK_CWNT_INFO.BLOCK[IVAR_CWNT] ;
    sprintf(VARSET, "$SET:%s", varName);
        
    HBOOK_CWNT_READROW.ICAST[IVAR_READ] = ICAST_CWNT ; 

    if ( ICAST_CWNT == ICAST_D ) { 
      ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_READ][ICAST_D];
      ptrTmp   = HBOOK_CWNT_READROW.VAL_D[ivarcast]; 

      // check block-name to see if this is really a 'long long int'
      if(strstr(blkName,"_L")!= NULL) { IVAR_L[NLL]=IVAR_READ; NLL++; }

    }

    else if ( ICAST_CWNT == ICAST_F ) { 
      ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_READ][ICAST_F];
      ptrTmp = HBOOK_CWNT_READROW.VAL_F[ivarcast]; 
    }

    else if ( ICAST_CWNT == ICAST_I ) {
      ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_READ][ICAST_I];
      ptrTmp   = HBOOK_CWNT_READROW.VAL_I[ivarcast]; 
    }

    else if ( ICAST_CWNT == ICAST_C ) {
      //  MEMC = CBUFF_CWNT * sizeof(char);
      //  HBOOK_CWNT_READROW.VAL_C[IVAR_READ] = (char*)malloc(MEMC);
      ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_READ][ICAST_C];
      for(ep=0; ep < MXEPOCH; ep++ ) { 
	ptrC[ep]  = HBOOK_CWNT_READROW.VAL_C[ivarcast][ep]; 
	sprintf( ptrC[ep]," ");
      } 
    }

    else {
      sprintf(MSGERR1,"Invalid ICAST=%d for IVAR_READ=%d (%s)",
	      ICAST_CWNT, IVAR_READ, varName);
      sprintf(MSGERR2,"Check valid ICAST_D/F/I/C " );
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);       
    }

    LENB    = strlen(blkName);
    LENV    = strlen(VARSET);
    LENC    = CBUFF_CWNT ;   // for hbnamc only
    if ( ICAST_CWNT == ICAST_C )  
      { hbnamc_(&NTID, blkName, *ptrC,  VARSET, LENB,LENC,LENV);  }    
    else
      { hbname_(&NTID, blkName, ptrTmp, VARSET, LENB,LENV); }
        
  }  // IVAR_READ


  
  // loop over each event (ROW) and extract CWNT values
  int   EVT, irow, IERR ;
  NEVT = HBOOK_CWNT_INFO.NEVENT ; 

  for ( EVT=1; EVT <= NEVT; EVT++ ) {  // note fortran-line index

    hgnt_(&NTID, &EVT, &IERR);  // read table-row here 

    irow = EVT - 1; // C-like row index

    // fill pointers or write to dump file; pass C-like index
    if ( OUTLIER_INFO.USEFLAG ) 
      { outlier_dump_fromHbook(irow); }
    else 
      { sntable_pushRowOut_hbook(irow,OPT_READ,0); }

  } // end evt  loop

  return NEVT ;


} // end of SNTABLE_READ_EXEC_HBOOK


// =========================================================
void malloc_READROW_HBOOK(int OPT) {

  // Malloc  HBOOK_CWNT_READROW.VAL_[cast] arrays.
  // OPT > 0 --> malloc
  // OPT < 0 --> free

  int icast, ICAST_CWNT, NVARCAST, MEM, NEP, ep ;
  int ivar, IVAR_READ, IVAR_TOT, IVAR_CWNT ;
  int NVAR_READ = READTABLE_POINTERS.NVAR_READ ;
  char *varName ;
  //  char fnam[] = "malloc_READROW_HBOOK" ;

  // ----------- BEGIN ------------

  for(icast=0; icast < MXCAST; icast++ ) {
    HBOOK_CWNT_READROW.NVARCAST[icast] = 0 ;
    for ( IVAR_READ=0; IVAR_READ < NVAR_READ; IVAR_READ++ ) 
      { HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_READ][icast] = -9 ; }
  }

  for ( IVAR_READ=0; IVAR_READ < NVAR_READ; IVAR_READ++ ) {
    IVAR_TOT     = READTABLE_POINTERS.PTRINDEX[IVAR_READ] ;
    varName      = READTABLE_POINTERS.VARNAME[IVAR_TOT] ;
    IVAR_CWNT    = get_cwnt_ivar(varName);
    ICAST_CWNT   = HBOOK_CWNT_INFO.ICAST[IVAR_CWNT] ;

    HBOOK_CWNT_READROW.NVARCAST[ICAST_CWNT]++ ;
    NVARCAST = HBOOK_CWNT_READROW.NVARCAST[ICAST_CWNT];
    HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_READ][ICAST_CWNT] = NVARCAST-1;
  }

  printf("  Malloc %d READROWs : %d(D) + %d(F) + %d(I) + %d(LL) + %d(C) \n"
	 ,NVAR_READ
	 ,HBOOK_CWNT_READROW.NVARCAST[ICAST_D]
	 ,HBOOK_CWNT_READROW.NVARCAST[ICAST_F]
	 ,HBOOK_CWNT_READROW.NVARCAST[ICAST_I]
	 ,HBOOK_CWNT_READROW.NVARCAST[ICAST_L]
	 ,HBOOK_CWNT_READROW.NVARCAST[ICAST_C]  
	 );  fflush(stdout);

  // allocate READROW memory
  NEP = MXEPOCH; // check laster for NPT

  NVARCAST = HBOOK_CWNT_READROW.NVARCAST[ICAST_D] + 1 ;
  HBOOK_CWNT_READROW.VAL_D = (double**)malloc( NVARCAST * sizeof(double*) );
  MEM = NEP * sizeof(double) ;
  for(ivar=0; ivar < NVARCAST; ivar++ ) 
    { HBOOK_CWNT_READROW.VAL_D[ivar] = (double*)malloc(MEM); }

  NVARCAST = HBOOK_CWNT_READROW.NVARCAST[ICAST_F] + 1 ;
  HBOOK_CWNT_READROW.VAL_F = (float**)malloc( NVARCAST * sizeof(float*) );
  MEM = NEP * sizeof(float) ;
  for(ivar=0; ivar < NVARCAST; ivar++ ) 
    { HBOOK_CWNT_READROW.VAL_F[ivar] = (float*)malloc(MEM); }

  NVARCAST = HBOOK_CWNT_READROW.NVARCAST[ICAST_I] + 1 ;
  HBOOK_CWNT_READROW.VAL_I = (int**)malloc( NVARCAST * sizeof(int*) );
  MEM = NEP * sizeof(int) ;
  for(ivar=0; ivar < NVARCAST; ivar++ ) 
    { HBOOK_CWNT_READROW.VAL_I[ivar] = (int*)malloc(MEM);    }

  NVARCAST = HBOOK_CWNT_READROW.NVARCAST[ICAST_C] + 1 ;
  HBOOK_CWNT_READROW.VAL_C = (char***)malloc( NVARCAST * sizeof(char**) );
  MEM = NEP * sizeof(char*) ;
  for(ivar=0; ivar < NVARCAST; ivar++ )  { 
    HBOOK_CWNT_READROW.VAL_C[ivar] = (char**)malloc(MEM); 
    for(ep=0; ep < NEP; ep++ ) {
      HBOOK_CWNT_READROW.VAL_C[ivar][ep] = (char*)malloc(32*sizeof(char) ); 
    }
  }

 
  //  printf("\n xxx DEBUG ABORT xxxx \n");  exit(666); // xxx REMOVE

  return ;

} // end malloc_READROW_HBOOK

// =========================================================
void outlier_dump_fromHbook(int irow) {

  // Created Oct 26 2014
  // hgnt has already been called; loop over each epoch
  // and dump to ascii file. 
  // Nov 30 2020: fix nasty bug setting NPTFIT

  float CHI2, CHI2MIN, CHI2MAX; ;
  int IVAR_NPT, IVAR_IFILT, IVAR_CHI2, IVAR_CUTFLAG, IFILT ;
  int NPT, ep, ivar, ivarcast ;
  int LDMP = (irow == -888 ) ;
  char fnam[] = "outlier_dump_fromHbook" ;

  // ------------- BEGIN ------------

  // strip off global OUTLIER info into local variables
  CHI2MIN     = OUTLIER_INFO.CUTWIN_CHI2FLUX[0] ;
  CHI2MAX     = OUTLIER_INFO.CUTWIN_CHI2FLUX[1] ;
  IVAR_NPT    = OUTLIER_INFO.IVAR[INDX_OUTLIER_NPTFIT];
  IVAR_IFILT  = OUTLIER_INFO.IVAR[INDX_OUTLIER_IFILT];
  IVAR_CHI2   = OUTLIER_INFO.IVAR[INDX_OUTLIER_CHI2];

  ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_NPT][ICAST_I];
  NPT  = HBOOK_CWNT_READROW.VAL_I[ivarcast][0]; 

  if ( LDMP ) 
    {  printf(" xxx %s: irow=%3d  IVAR_NPT=%3d of %3d\n", 
	      fnam, irow, IVAR_NPT, NPT ); fflush(stdout); 
    }

  // loop over each fitted epoch and check if it's an outlier
  for(ep=0; ep < NPT; ep++ ) {

    ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_CHI2][ICAST_F];
    CHI2  = HBOOK_CWNT_READROW.VAL_F[ivarcast][ep];

    ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_IFILT][ICAST_I];
    IFILT = HBOOK_CWNT_READROW.VAL_I[ivarcast][ep];
    
    OUTLIER_INFO.NEP_TOT[0]++ ;      // total summed over filters
    OUTLIER_INFO.NEP_TOT[IFILT]++ ; // total for this band only
    
    if ( CHI2 >= CHI2MIN && CHI2 <= CHI2MAX ) {

      // copy scalar quanties into each vector element
      for( ivar=0; ivar <= IVAR_NPT; ivar++ ) {

	ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[ivar][ICAST_D];
	if ( ivarcast >=0 ) {
	  HBOOK_CWNT_READROW.VAL_D[ivarcast][ep] = 
	    HBOOK_CWNT_READROW.VAL_D[ivarcast][0] ;
	}

	ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[ivar][ICAST_F];
	if ( ivarcast >=0 ) {
	  HBOOK_CWNT_READROW.VAL_F[ivarcast][ep] = 
	    HBOOK_CWNT_READROW.VAL_F[ivarcast][0] ;
	}

	ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[ivar][ICAST_I];
	if ( ivarcast >=0 ) {
	  HBOOK_CWNT_READROW.VAL_I[ivarcast][ep] = 
	    HBOOK_CWNT_READROW.VAL_I[ivarcast][0] ;
	}

	ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[ivar][ICAST_C];
	if ( ivarcast >=0 ) {
	  sprintf(HBOOK_CWNT_READROW.VAL_C[ivar][ep], "%s",
		  HBOOK_CWNT_READROW.VAL_C[ivar][0]);
	}

      }
      
      //  sntable_dump_update(FP_DUMP, LINEKEY_SN, ep, FARRAY_OUT);  
      sntable_pushRowOut_hbook(irow,OPT_SNTABLE_READ_forDUMP,ep); 
      
      OUTLIER_INFO.NEP_SELECT[0]++ ;
      OUTLIER_INFO.NEP_SELECT[IFILT]++ ;
    }
  } // end ep loop

  return ;

} // end of outlier_dump_fromHbook

// =========================================================
void sntable_pushRowOut_hbook(int IROW, int OPT_READ, int IFIT) {

  // Oct 26 2014:
  // Use contents of ntuple row (HBOOK_CWNT_READROW struct)
  // and either fill user-defined pointers to arrays
  // or write to user-defined ascii file.

  // Input: 
  //  IROW       : row number (starts at 0)
  //  OPT_READ   : see #define OPT_SNTABLE_XXX above.
  //  IFIT       : vector index (pass 0 for scaler)

  // Feb 24 2019:  check SEPKEY

  //  char fnam[] =  "sntable_pushRowOut_hbook" ;

  int NVAR_DUMP = READTABLE_POINTERS.NVAR_READ ;

  int IVAR_DUMP, IVAR_TOT, ICAST_READ, LDUMP, LREAD, ivarcast, ADDSEPKEY ;
  int OPT, NVAR_WR=0;

  int    VAL_I ;
  double VAL_D ;
  float  VAL_F ;
  char  *ptrC  ;

  char   *SEPKEY = READTABLE_POINTERS.SEPKEY_DUMP ;
  char   LINE[MXCHAR_VARLIST], *VARNAME, band[2] ;
  char   BLANK[] = " " ;
  char   fnam[] = "pushRowOut" ;

  // ------------ BEGIN ---------
  
  LDUMP = ( OPT_READ == OPT_SNTABLE_READ_forDUMP  );
  LREAD = ( OPT_READ == OPT_SNTABLE_READ_forARRAY );

  sprintf(LINE,"%s ", READTABLE_POINTERS.LINEKEY_DUMP); 
  ADDSEPKEY = (strlen(SEPKEY) > 0 ) ;

  for ( IVAR_DUMP = 0; IVAR_DUMP < NVAR_DUMP; IVAR_DUMP++ ) {

    IVAR_TOT   = READTABLE_POINTERS.PTRINDEX[IVAR_DUMP] ;
    VARNAME    = READTABLE_POINTERS.VARNAME[IVAR_TOT] ;
    ICAST_READ = HBOOK_CWNT_READROW.ICAST[IVAR_DUMP] ;
    VAL_D      = -9999. ;

    OPT=0;  if ( ISTABLEVAR_IFILT(VARNAME) ) { OPT=1; }

    if ( ICAST_READ == ICAST_I ) { 
      ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_DUMP][ICAST_I];
      VAL_I  = HBOOK_CWNT_READROW.VAL_I[ivarcast][IFIT] ;
      VAL_D   = (double)VAL_I ; 

      if ( LDUMP )  { load_DUMPLINE(OPT, LINE, VAL_D); }
      if ( LREAD )  { load_READTABLE_POINTER(IROW,IVAR_TOT,VAL_D,BLANK); }
    }
    
    else if ( ICAST_READ == ICAST_F )  { 
      ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_DUMP][ICAST_F];
      VAL_F = HBOOK_CWNT_READROW.VAL_F[ivarcast][IFIT]  ;
      VAL_D = (double)VAL_F ;   
      if ( LDUMP )  { load_DUMPLINE(OPT,LINE, VAL_D); }
      if ( LREAD )  { load_READTABLE_POINTER(IROW,IVAR_TOT,VAL_D,BLANK); }
    }
    
    else if ( ICAST_READ == ICAST_D )  { 
      ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_DUMP][ICAST_D];
      VAL_D = HBOOK_CWNT_READROW.VAL_D[ivarcast][IFIT] ;
      if ( LDUMP )  { load_DUMPLINE(OPT, LINE, VAL_D); }
      if ( LREAD )  { load_READTABLE_POINTER(IROW,IVAR_TOT,VAL_D,BLANK); }
    }
    
    else if ( ICAST_READ == ICAST_C ) {
      ivarcast = HBOOK_CWNT_READROW.IVARCAST_MAP[IVAR_DUMP][ICAST_C];
      ptrC   =  HBOOK_CWNT_READROW.VAL_C[ivarcast][IFIT]; 
      trim_blank_spaces(ptrC);  
      load_DUMPLINE_STR(LINE,ptrC);
      if ( LREAD )  { load_READTABLE_POINTER(IROW,IVAR_TOT,VAL_D,ptrC); }
    }

    /*
    printf(" xxx IVAR=%d of %d: LDUMP=%d ICAST=%d \n", 
	   IVAR, NVAR, LDUMP, ICAST_READ); fflush(stdout);
    */

    // Check option to write sep-string between variables;
    // in particular, a comma for csv format.
    if ( LDUMP && ICAST_READ > 0 && ADDSEPKEY ) { 
      NVAR_WR++ ;
      if ( NVAR_WR < NVAR_DUMP ) { strcat(LINE,SEPKEY); }
    }

  } // end of IVAR


  // ------------------------------------------
  // update ascii file for dump option
  if ( LDUMP )  { 
    FILE *FP = READTABLE_POINTERS.FP_DUMP ;
    fprintf(FP, "%s\n", LINE);  
    fflush(FP); 
  }

  return ;

} // end of  void sntable_pushRowOut_hbook

// =========================================================
void SNTABLE_DUMP_VARNAMES_HBOOK(char *fileName, int NTID) {

  // Screen-dump list of all variable names in CWNT.
  // Input args are the hbook *fileName and ntuple id (NTID).
  // If fileName is blank, assume that file is already open.

  int ivar, NVAR ;
  //  char fnam[] = "SNTABLE_DUMP_VARNAMES_HBOOK" ;

  // -------------- BEGIN ---------------

  get_cwnt_info(fileName, NTID, -1 ); // open/read/close

  NVAR = HBOOK_CWNT_INFO.NVAR ;
  for ( ivar=0; ivar < NVAR; ivar++ )  { 
    printf("%-24s  %s  BLOCK=%s   IVAR=%d\n"
	   , HBOOK_CWNT_INFO.VARNAME[ivar] 
	   , HBOOK_CWNT_INFO.CAST[ivar] 
	   , HBOOK_CWNT_INFO.BLOCK[ivar] 
	   , ivar
	   );   
    fflush(stdout); 
  }
  

  printf("\n Done listing %d variables for NTID=%d \n", NVAR, NTID);
  printf(" Number of table rows: %d\n", HBOOK_CWNT_INFO.NEVENT);
  fflush(stdout) ; 


  return ;

} // end of SNTABLE_DUMP_VARNAMES_HBOOK



// ==============================================
void get_cwnt_info(char *fileName, int NTID, int OPT) {

  // Apr 2013
  // open hbook file in READ-mode and hrin NTID so that it 
  // is ready for access. For each variable store name, 
  // cast and CHBLOK
  //
  // (I) fileName = file name
  // (I) NTID     = ntuple id
  // (I) OPT      = 1-> leave file open;  -1-> close if fileName != ""
  //                   (if fileName=='', always leave open)
  //
  
  // OPEN_HFILE args
  int IERR ;
  int LUN = 20 ;
  char copt[8]  = "Q" ;  // quiet mode opening file

  // hrin args
  int ICYC = 999 ;
  int IOFF = 0 ;

  // misc args
  int  NVAR, MEMC, LENF, ITYPE, IVAR_CWNT, ICAST  ;
  char VARCAST_NAME[CBUFF_CWNT]; // name of variable with cast, such as RA:R*8
  char BLOCK[CBUFF_CWNT];
  char *ptrVar, *ptrBlk, *ptrCast ;

  char fnam[] = "get_cwnt_info" ;

  // --------- BEGIN ----------

  if ( NTID <= 0 ) {
    sprintf(MSGERR1,"Invalid ntuple id  NTID = %d .", NTID);
    sprintf(MSGERR2,
	    "Check that ntuple TABLE_ID is integer instead of string.");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }


  // open hbook file if fileName is given; 
  // otherwise assume it's already open

  LENF = strlen(fileName) ;

  if ( LENF > 0 ) {

    sprintf(HBOOK_TOPDIR_READFILE, "%s", HBOOK_TOPDIRNAM_READ );
    OPEN_HFILE(fileName, LUN, copt, HBOOK_TOPDIR_READFILE, &IERR) ;
  
    if ( IERR != 0 ) {
      sprintf(MSGERR1, "Could not open hbook file");
      sprintf(MSGERR2, "%s", fileName);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }
  }


  // read in the CWNT; note that ICYC is large to read
  // last version created by merging program.
  hrin_(&NTID, &ICYC, &IOFF);
  hnoent_(&NTID, &HBOOK_CWNT_INFO.NEVENT);
  if ( HBOOK_CWNT_INFO.NEVENT <= 0 ) {
    sprintf(MSGERR1, "No CWNT entries for NTID=%d", NTID);
    sprintf(MSGERR2, "Check hbook file '%s'", fileName);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }
  
  // ------------------------------------------
  // loop over each variable to get block & cast

  
  ITYPE = 1 ;
  NVAR  = 0 ;
  MEMC  = CBUFF_CWNT * sizeof(char) ;

  while ( ITYPE != 0 ) {

    IVAR_CWNT = NVAR + 1 ;
    sprintf(VARCAST_NAME," ");    sprintf(BLOCK," ");

    hntvdef_(&NTID, &IVAR_CWNT, VARCAST_NAME, BLOCK, &ITYPE, 
	     CBUFF_CWNT, CBUFF_CWNT );

    if ( ITYPE > 0 ) {

      HBOOK_CWNT_INFO.VARNAME[NVAR] =  (char*)malloc(MEMC);
      HBOOK_CWNT_INFO.BLOCK[NVAR]   =  (char*)malloc(MEMC);
      HBOOK_CWNT_INFO.CAST[NVAR]    =  (char*)malloc(MEMC) ;

      ptrVar   = HBOOK_CWNT_INFO.VARNAME[NVAR] ;
      ptrBlk   = HBOOK_CWNT_INFO.BLOCK[NVAR] ;
      ptrCast  = HBOOK_CWNT_INFO.CAST[NVAR] ;

      // break up VARCAST_NAME info ptrVar and ptrCast; return ICAST index
      ICAST =  parse_cwnt_vardef(VARCAST_NAME, 
				 ptrVar, ptrCast);  // return args

      HBOOK_CWNT_INFO.ICAST[NVAR]    = ICAST ;  // store in global

      // strip off extra spaces in BLOCK name
      trim_blank_spaces(BLOCK);
      sprintf(ptrBlk,"%s", BLOCK);


      /*
      printf(" xxx %3d : %'-24.24s' -> ICAST=%d  BLOCK='%s'\n", 
	     NVAR, VARCAST_NAME, ICAST, ptrBlk );
      */

      NVAR++ ;
    }
  }

  HBOOK_CWNT_INFO.NVAR = NVAR ;

  // close file if OPT is negative, and also if fileName is non-null.
  if ( OPT < 0 && LENF > 0 ) 
    { CLOSE_HFILE(fileName,OPENFLAG_READ); }

  return ;

} // end of get_cwnt_info


// =============================================
int parse_cwnt_vardef(char *VARCAST_NAME, char *VARNAME, char *CAST ) {

  // parse input VARCAST_NAME which has the form  VARNAME:CAST
  // such as RA:R*8  or ZCMB:R*4 or  TYPE:I*4.
  //
  // Returns 
  //  * VARNAME = string before semicolon
  //  * CAST    = string after semicolon
  //  * ICAST   = function arg = integer code for CAST = ICAST_D/F/I/C
  //

  int  icast, i, FOUND_COLON, DONE ;
  char ctmp[2], CAST1[2];

  char fnam[] = "parse_cwnt_vardef" ;
 
  // ---------- BEGIN ---------

  icast = -9 ;
  VARNAME[0] = 0 ;
  CAST[0] = 0 ;

  FOUND_COLON = DONE = 0; 
  sprintf(ctmp,"X");

  i=0;

  while ( DONE == 0 ) {
    sprintf(ctmp, "%c", VARCAST_NAME[i] );
    i++ ;

    if( strcmp(ctmp,":") == 0 ) { FOUND_COLON = 1 ; continue ; }
    if( strcmp(ctmp," ") == 0 ) { DONE        = 1 ; continue ; }

    if ( FOUND_COLON ) 
      { sprintf(CAST,"%s%s", CAST, ctmp); }
    else
      { sprintf(VARNAME,"%s%s", VARNAME, ctmp); }


    if ( i > 50 ) {
      sprintf(MSGERR1,"i=%d -> too many characters for", i);
      sprintf(MSGERR2,"VARCAST_NAME = '%s'", VARCAST_NAME);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
    }
  }

  sprintf(CAST1,"%c", CAST[0] );
  if ( strcmp(CAST, "R*8") == 0 ) { icast = ICAST_D ; }
  if ( strcmp(CAST, "R*4") == 0 ) { icast = ICAST_F ; }
  if ( strcmp(CAST, "I*4") == 0 ) { icast = ICAST_I ; }
  if ( strcmp(CAST, "L*8") == 0 ) { icast = ICAST_L; } 
  if ( strcmp(CAST1,"C"  ) == 0 ) { icast = ICAST_C ; }
  

  /*  
  printf(" xxx %-20s -> \n", VARCAST_NAME );
  printf("     VARNAME = '%s'  CAST = '%s' \n", VARNAME, CAST);
  fflush(stdout);
  */

  if ( icast < 0 ) {
    sprintf(MSGERR1,"Could not determine CWNT cast for CAST = '%s'", CAST);
    sprintf(MSGERR2,"VARCAST_NAME = '%s'", VARCAST_NAME);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  return icast ;

} // end of parse_cwnt_vardef


// ==========================
int get_cwnt_ivar(char *varName) {

  // return absolute CWNT variable index for *varName.
  // ABORT if *varName is not found.
  //
  // Feb 2014: make check case-insensitive by converting
  //            everything to lower case before checking.
  //
  // Aug 2014: exclude () or [] from VARNAME to allow array
  //           See new function call to removeParenth().
  //

  int  NVAR, ivar, ivar_cwnt, LEN1, LEN2;
  char tmpName[60] ;
  char fnam[] = "get_cwnt_ivar" ;

  // --------------------- BEGIN -----------------

  ivar_cwnt = -9 ;


  LEN1 = strlen(varName);
  NVAR = HBOOK_CWNT_INFO.NVAR ;

  for(ivar=0; ivar < NVAR ; ivar++ ) {

    sprintf(tmpName, "%s", HBOOK_CWNT_INFO.VARNAME[ivar]);
    removeParenth(tmpName); // Aug 2 2014

    LEN2  = strlen(tmpName);
    if ( LEN1 != LEN2 ) { continue ; }

    if ( strcmp_ignoreCase(tmpName,varName) == 0 ) 
      { ivar_cwnt = ivar; return ivar_cwnt ; }

  } // ivar
  
  
  // if we get here then ABORT on error
  sprintf(MSGERR1,"Could not find CWNT variable '%s'", varName);
  sprintf(MSGERR2,"Check VARNAMES argument(s)" );
  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 

  return ivar_cwnt;

} // end of get_cwnt_ivar


// =================================
void removeParenth(char *string) {

  // Created Aug 2 2014
  // if input string is 'bla(NFIT)'   then return 'bla'.
  // if input string is 'NFIT[0,600]' then return 'NFIT'.
  // Input string is overwritten

  char newString[200], c1[2];
  int  LEN, i ;

  sprintf(newString,"%s", "");
  LEN = strlen(string);

  for(i=0; i < LEN; i++ ) {
    sprintf(c1, "%c", string[i] );
    if ( *c1 == '[' ) { goto COPYSTRING ; }
    if ( *c1 == '(' ) { goto COPYSTRING ; }
    strcat(newString,c1);
  }

 COPYSTRING:
  sprintf(string, "%s", newString);

} // end of removeParenth


// ======================================================
//
//               SNHIST FUNCTIONS
//
// ======================================================


// --------------------------------------------------
void SNHIST_INIT_HBOOK(int NDIM, int ID, char *TITLE, 
			int *NBIN, double *XMIN, double *XMAX ) {

  int idim ;
  int IDLOC = ID ;
  float XMIN4[10], XMAX4[10] ;
  float w4 = 0.0 ;
  char pawTITLE[100];

  // ------------- BEGIN ---------------

  sprintf(pawTITLE,"%s", TITLE);
  modify_text_4paw(pawTITLE) ; // make substitutions for hbook chars

  // generate float list of limits in case it's needed.
  for(idim=0; idim<NDIM; idim++ )  {  
    XMIN4[idim] = (float)XMIN[idim];
    XMAX4[idim] = (float)XMAX[idim];
  }
  
  if ( NDIM == 1 ) {
    hbook1_(&IDLOC, pawTITLE, &NBIN[0], &XMIN4[0], &XMAX4[0], &w4, 
	    strlen(pawTITLE)  );
  }
  else if ( NDIM == 2 ) {     
    hbook2_(&IDLOC, pawTITLE
	    ,&NBIN[0], &XMIN4[0], &XMAX4[0]
	    ,&NBIN[1], &XMIN4[1], &XMAX4[1], &w4, strlen(pawTITLE) ); 
  }


}  // end of  SNHIST_INIT_HBOOK


// --------------------------------------------------
void SNHIST_FILL_HBOOK(int NDIM, int ID, double *VALUE, double WGT) {
  
  int idim, IDLOC;
  float W4, VAL4[10] ;
  
  // ---------------- BEGIN ----------
  
  IDLOC = ID ;
  W4    = (float)WGT ;
  for ( idim=0; idim<NDIM; idim++ ) { VAL4[idim] = (float)VALUE[idim]; }

  if ( NDIM == 1 ) 
    { hf1_ ( &IDLOC, &VAL4[0], &W4) ; }  // pass float instead of double
  else if ( NDIM == 2 ) 
    { hf2_ ( &IDLOC, &VAL4[0], &VAL4[1], &W4) ; }

}  // end of SNHIST_FILL_HBOOK



// =======================================================
void SNHIST_RDBINS_HBOOK(int NDIM, int ID, char *CTIT, int *NBIN,
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

  int HID  = ID  ;
  int ICYC = 999 ;
  int IOFF = 0   ;

  int   NBX, NBY, NWT, LOC ;
  int   LENTIT = 80 ;
  float X0, X1, Y0, Y1 ;

  // -------------- BEGIN -------------

  NBIN[0] = 0 ; // init output number of bins

  hrin_(&HID, &ICYC, &IOFF);  
  hgive_(&HID, CTIT, &NBX,&X0,&X1, &NBY,&Y0,&Y1, &NWT, &LOC, LENTIT);
	  
  
  // load output function args
  if ( NBX > 0 ) 
    {  NBIN[0] = NBX ;  XMIN[0] = X0 ;  XMAX[0] = X1 ; } // x-axis

  if ( NBY > 0 ) 
    {  NBIN[1] = NBY ;  XMIN[1] = Y0 ;  XMAX[1] = Y1 ; } // y-axis

}  // end of SNHIST_RDBINS_HBOOK



void SNHIST_RDCONT_HBOOK(int NDIM, int ID, int NRDBIN, double *CONTENTS) {

  // June 2013
  // return CONTENTS of histogram ID
  // Inputs:
  //    ID = histo id
  //    NRDBINS = number of bins to read (see SNHIST_RDBINS)
  // Output:
  //   CONTENTS[0 - NRDBINS-1]
  
  int HID = ID ;
  int NUM = 1  ; 
  int ibin ;
  float *farray ;
  char  choice[] = "HIST" ;
  
  farray = (float*)malloc(NRDBIN * sizeof(float) );
  hunpak_(&HID, farray, choice, &NUM, strlen(choice) );

  for (ibin=0; ibin < NRDBIN; ibin++ ) 
    {  CONTENTS[ibin] = (double)farray[ibin] ; }


  free(farray);

} // end of SNHIST_RDCONT_HBOOK


// ======================================================
//
//               DIRECTORY FUNCTIONS
//
// ======================================================


// --------------------------------------------------
void MAKEDIR_HBOOK(char *CCID, int CID) {

  // Oct 24 2017:  CID< 1000000 -> CID <= 999999
  //       to avoid install script from replacing 1000000 with 500000
  //                   
  int  NFIT, IFIT ;
  char sdir[60];

  char topDirName[20] = "       ";
  char subDirName[60] ;
  char pawDirName[60] ;
  char R[4] = "R " ;
  char S[4] = "S " ;
  char fnam[] = "MAKEDIR_HBOOK" ;

  // ------------------- BEGIN --------------

  if ( CID < 0 ) 
    {  sprintf(sdir,"%s", CCID); }  // CCID is just a dirName
  else if ( CID <= 999999 )         
    {  sprintf(sdir,"SN%6.6d", CID); } // SN[CID] dir
  else
    {  sprintf(sdir,"SN%8.8d", CID); } // SN[CID] dir with longer name

 
  //  add _[IFIT] suffix to subdir if there are multiple fits per SN
  if ( CID >= 0 ) {
    NFIT   = SNLCPAK_OUTPUT.NFIT_PER_SN  ;
    IFIT   = SNLCPAK_OUTPUT.NLCPAK  ; // starts at 0 since it sums later
    if ( NFIT > 1 ) { sprintf(sdir,"%s_%2.2d", sdir, IFIT); }
  }

  hcdir_(topDirName, R, strlen(topDirName), strlen(R)); // return topDir
  sprintf(pawDirName,"//PAWC/%s", sdir);
  sprintf(subDirName,"%s/%s", topDirName, sdir);

  hmdir_ ( subDirName, S, strlen(subDirName), strlen(S) ) ;
  hmdir_ ( pawDirName, S, strlen(pawDirName), strlen(S) ) ;

  printf("\t %s created  %s \n", fnam, subDirName );
  fflush(stdout);
  
} // end of MAKEDIR_HBOOK


// --------------------------------------------------
void CDTOPDIR_HBOOK(void) {

  char fnam[] = "CDTOPDIR_HBOOK" ;
  char *topDir ;
  char PAWC[]  = "//PAWC" ;
  char blank[] = " ";

  // -------------- BEGIN ----------
  topDir = HBOOK_TOPDIR_NEWFILE ; 
  hcdir_( PAWC,  blank, strlen(PAWC), strlen(blank) ) ;
  hcdir_( topDir, blank, strlen(topDir), strlen(blank)  ) ;
  printf("\t %s: return to %s \n", fnam, topDir ) ;
  fflush(stdout) ;

} // end of CDTOPDIR_HBOOK



// ======================================================
//
//               SNLCPAK FUNCTIONS
//
// ======================================================


// --------------------------------------------------
void SNLCPAK_SURVEY_HBOOK(void) {

  // Aug 29 2017: extend NFIT from 4 to 10

  //  char fnam[] = "SNLCPAK_SURVEY_HBOOK" ;
  // ------------------------------------------------
  // book info histograms in top directory
  int    NBIN, HID ;
  float  XMIN, XMAX, X ;
  float  WGT0 = 0.0 ;
  float  WGT1 = 1.0 ;
  char   CTIT[200];

  NBIN = 11 ;  XMIN = -0.5; XMAX = 10.5; 

  HID = 12 ;
  sprintf(CTIT,"%s  %s  %s"
	  ,SNLCPAK_OUTPUT.SURVEY
	  ,SNLCPAK_OUTPUT.VERSION_PHOTOMETRY
	  ,SNLCPAK_OUTPUT.SURVEY_FILTERS      );

  hbook1_(&HID, CTIT, &NBIN, &XMIN, &XMAX, &WGT0, strlen(CTIT) ) ;
  X   = (float)SNLCPAK_OUTPUT.NFIT_PER_SN ;  
  hf1_(&HID, &X, &WGT1);


  // store list of all filters in histogram title so that 
  // plot-macro can access cfilt = FILTERSTRING[ifiltobs] 
  HID = 14 ;
  hbook1_(&HID, FILTERSTRING, &NBIN,&XMIN,&XMAX,&WGT0, strlen(FILTERSTRING));
  X = 0.0 ;
  hf1_(&HID, &X, &WGT1);

  return ;

} // end of SNLCPAK_SURVEY_HBOOK


// --------------------------------------------------
void SNLCPAK_FILL_HBOOK(void) {

  // Mar 24 2013: fix bug booking Tobs for data: NOBS instead of NBIN
  // 
  // Jan 26, 2014: add FLAG_FLUXREST and FLAG_SIMFLUXREST
  //               Also add abort-trap if any FLAG is not defined.
  //               
  // May 24 2016: define HOFF_EPOCH and CTAG_EPOCH for 
  //              SNLCPAK_EPFLAG_FLUXSIM, even though nothing is
  //              plotted for true simFlux.
  //
  // Oct 3 2018: NBIN -> NOBS_FITFUN in Tobs plot for FITFUN (bugfix)
  //
  // Mar 8 2020: define HOFF_EPOCH and CTAG_EPOCH for kcor and AVwarp.
  //
  int  
    NOBS, FLAG, obs, firstObs, NOBS_FITFUN, NOBS_CHI2
    ,IFILT, USE, OVFIT, iplot, NFIT
    ,NFILTDEF_SURVEY, itext
    ;   

  char 
    cfilt[2]
    ,*SURVEY_FILTERS
    ,CTAG_EPOCH[MXFLAG_SNLCPAK_EPOCH][40]
    ,CTAG_FILTER[MXFLAG_SNLCPAK_FILTER][40]
    ,CTIT[80]
    ,*CCID
    ,HTEXT[100]
    ,fnam[] = "SNLCPAK_FILL_HBOOK" 
    ;

  int 
    HID, NBIN
    ,HOFF_EPOCH[MXFLAG_SNLCPAK_EPOCH] 
    ,HID_FILTER[MXFLAG_SNLCPAK_FILTER] 
    ;

  float  XMIN_F, XMAX_F, X_F ;
  float  WGT_F, WGT0_F = 0.0 ;
  double  X, EX ;
  double *ptrTOBS, *ptrDATA, *ptrERR ;

  float 
     *DATA4
    ,*DATA_ERR4
    ,NULL4[MXFLAG_SNLCPAK_EPOCH] 
    ;

  int HID_FILTLIST     = 1 ;     // ctit = filter list
  int HID_NFIT         = 2 ;
  int HID_CCID         = 3 ;
  int HOFF_TEXT        = 50 ;
  int HID_DATA_TOBS    = 1000 ;
  int HID_FITFUN_TOBS  = 3000 ;  // TOBS list

  // -------------- BEGIN ---------------

  // allocate DATA4 and DATA_ERR4 arrays for HPAK
  NOBS      = SNLCPAK_OUTPUT.NOBS_MAX ;
  DATA4     = (float*)malloc(NOBS * sizeof(float) ) ;
  DATA_ERR4 = (float*)malloc(NOBS * sizeof(float) ) ;

  // set dummy values to make sure everything is set below.
  for ( FLAG=0; FLAG < MXFLAG_SNLCPAK_EPOCH ; FLAG++ ) {
    HOFF_EPOCH[FLAG]    =  -9 ;
    sprintf(CTAG_EPOCH[FLAG],  "UNDEFINED" );
  }

  // hard-wire hbook offsets and tag-strings
  HOFF_EPOCH[SNLCPAK_EPFLAG_FLUXDATA]    =  100 ;
  HOFF_EPOCH[SNLCPAK_EPFLAG_REJECT]      = 1300 ;
  HOFF_EPOCH[SNLCPAK_EPFLAG_CHI2]        =  400 ;
  HOFF_EPOCH[SNLCPAK_EPFLAG_KCOR]        = 2000 ; // Mar 2020
  HOFF_EPOCH[SNLCPAK_EPFLAG_AVWARP]      = 2400 ; // Mar 2020
  HOFF_EPOCH[SNLCPAK_EPFLAG_FITFUN]      = 3200 ;
  HOFF_EPOCH[SNLCPAK_EPFLAG_FLUXSIM]     = 3900 ;  
  HOFF_EPOCH[SNLCPAK_EPFLAG_FLUXREST]    = 4000 ;  
  HOFF_EPOCH[SNLCPAK_EPFLAG_SIMFLUXREST] = 5000 ;  

  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_FLUXDATA],    "Data Flux"      );
  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_REJECT],      "Reject Flag"    );
  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_CHI2],        "[h]^2!"         );
  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_KCOR],        "kcor"           );
  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_AVWARP],      "AVWARP"         );
  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_FITFUN],      "Fit-Model Flux" );
  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_FLUXSIM],     "Sim Flux"       );
  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_FLUXREST],    "Flux-rest"      );
  sprintf( CTAG_EPOCH[SNLCPAK_EPFLAG_SIMFLUXREST], "Sim Flux-rest"  );

  NULL4[SNLCPAK_EPFLAG_FLUXDATA]    = -99999. ;
  NULL4[SNLCPAK_EPFLAG_REJECT]      =  0.0 ;
  NULL4[SNLCPAK_EPFLAG_CHI2]        =  0.0 ;
  NULL4[SNLCPAK_EPFLAG_FITFUN]      =  0.0 ;
  NULL4[SNLCPAK_EPFLAG_FLUXREST]    = -99999. ;
  NULL4[SNLCPAK_EPFLAG_SIMFLUXREST] = -99999. ;

  // -----------
  HID_FILTER[SNLCPAK_BANDFLAG_NDOF   - 100]  = 4900 ;
  HID_FILTER[SNLCPAK_BANDFLAG_CHI2   - 100]  = 4901 ;
  HID_FILTER[SNLCPAK_BANDFLAG_PKFLUX - 100]  = 4910 ;
  HID_FILTER[SNLCPAK_BANDFLAG_PKMJD  - 100]  = 4920 ;

  sprintf( CTAG_FILTER[SNLCPAK_BANDFLAG_NDOF   -100], "Ndof"    );
  sprintf( CTAG_FILTER[SNLCPAK_BANDFLAG_CHI2   -100], "[h]^2!"  );
  sprintf( CTAG_FILTER[SNLCPAK_BANDFLAG_PKFLUX -100], "Peak Flux"    );
  sprintf( CTAG_FILTER[SNLCPAK_BANDFLAG_PKMJD  -100], "Peak MJD[-]MJDOFF" );


  SURVEY_FILTERS = SNLCPAK_OUTPUT.SURVEY_FILTERS ;
  CCID = SNLCPAK_OUTPUT.CCID ;
  
  // --------------------------------------------------
  // book info histo hid=1 with ctit = FILTLIST_PLOT and 
  // containing IFILT_MIN-MAX. Note we add 1 for fortran-like index
  HID = HID_FILTLIST ;
  sprintf( CTIT, "%s", SNLCPAK_OUTPUT.FILTLIST_PLOT );
  NBIN=2;  XMIN_F=0.5;  XMAX_F=2.5;
  hbook1_(&HID, CTIT, &NBIN, &XMIN_F, &XMAX_F, &WGT0_F, strlen(CTIT));
  X_F = 1.0 ;  WGT_F = (float)(SNLCPAK_OUTPUT.IFILT_MIN + 1) ;
  hf1_(&HID, &X_F, &WGT_F);
  X_F = 2.0 ;  WGT_F = (float)(SNLCPAK_OUTPUT.IFILT_MAX + 1) ;
  hf1_(&HID, &X_F, &WGT_F);

  NFIT = SNLCPAK_OUTPUT.NLCPAK_TOT ;
  if ( NFIT > 1 ) {
    HID = HID_NFIT ;
    sprintf( CTIT, "NFIT = %d for this SN", NFIT );
    NBIN=5;  XMIN_F=0.5;  XMAX_F=5.5;
    hbook1_(&HID, CTIT, &NBIN, &XMIN_F, &XMAX_F, &WGT0_F, strlen(CTIT));
    X_F = 2.0 ;  WGT_F = (double)NFIT ;
    hf1_(&HID, &X_F, &WGT_F);
  }


  HID = HID_CCID ;
  sprintf( CTIT, "%s", CCID );
  NBIN=1;  XMIN_F=0.0;  XMAX_F=1.0 ;
  hbook1_(&HID, CTIT, &NBIN, &XMIN_F, &XMAX_F, &WGT0_F, strlen(CTIT));

  // book histos with text strings. HOFF_TEXT is filled with NTEXT.
  // Mandatory string is constructed from SURVEY-FILTERS and VERSION.
  // Include only the 1st 15 filters, then add +... to indicate more.
  char cFILT_TMP[60];

  if ( strlen(SURVEY_FILTERS) < 15 )
    { sprintf(cFILT_TMP, "%s", SURVEY_FILTERS) ; }
  else
    { sprintf(cFILT_TMP, "%.*s...", 15, SURVEY_FILTERS); }


  NBIN=1;  XMIN_F=0.;  XMAX_F=1.0 ;
  sprintf(SNLCPAK_OUTPUT.DISPLAYTEXT[0],"%s-%s   SNID=%s (%s)" 
	  ,SNLCPAK_OUTPUT.SURVEY
	  ,cFILT_TMP
	  ,SNLCPAK_OUTPUT.CCID
	  ,SNLCPAK_OUTPUT.VERSION_PHOTOMETRY
	  );
  

  for(itext=0; itext <= SNLCPAK_OUTPUT.NTEXT; itext++ ) {
    HID     = HOFF_TEXT + itext ;
    sprintf(HTEXT,"%s", SNLCPAK_OUTPUT.DISPLAYTEXT[itext] ) ;
    modify_text_4paw(HTEXT) ; // make substitutions for hbook chars
    hbook1_(&HID, HTEXT, &NBIN, &XMIN_F, &XMAX_F, &WGT0_F, strlen(HTEXT));
    X_F = 0.5 ; WGT_F = (float)SNLCPAK_OUTPUT.NTEXT; 
    hf1_(&HID, &X_F, &WGT_F);
  }


  // --------------------------------------------------
  // book/fill 1D histogram for each defined FLAG and filter.

  NOBS_FITFUN = 0;
  NFILTDEF_SURVEY = SNLCPAK_OUTPUT.NFILTDEF_SURVEY ; 

  for ( FLAG = 0 ; FLAG < MXFLAG_SNLCPAK_EPOCH; FLAG++ ) {

    // skip FLUX_REST since there is nothing to plot it. (Jan 2014)
    if ( FLAG == SNLCPAK_EPFLAG_FLUXREST    ) { continue ; }
    if ( FLAG == SNLCPAK_EPFLAG_SIMFLUXREST ) { continue ; }

    NOBS = SNLCPAK_OUTPUT.NOBS[FLAG] ;
    if ( NOBS <= 0 ) { continue ; }

    if ( HOFF_EPOCH[FLAG] < 0 ) {
      sprintf(MSGERR1,"HOFF_EPOCH and CTAG_EPOCH are not defined for");
      sprintf(MSGERR2,"SNLCPAK_EPFLAG = %d", FLAG );
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }

    ptrDATA  = SNLCPAK_OUTPUT.EPDATA[FLAG];
    ptrERR   = SNLCPAK_OUTPUT.EPDATA_ERR[FLAG] ;

    for (IFILT=0; IFILT < NFILTDEF_SURVEY; IFILT++ ) {
      USE = SNLCPAK_OUTPUT.USEFILT[FLAG][IFILT] ;

      if ( USE ) {

	// for FITFUN, plot only for current filter since we know
	// that all filters have exactly the same TOBS range.
	if ( FLAG == SNLCPAK_EPFLAG_FITFUN ) {
	  NOBS_FITFUN = SNLCPAK_OUTPUT.NOBS_FILT[FLAG][IFILT];
	  NBIN        = NOBS_FITFUN ;
	}
	else
	  { NBIN      = NOBS ; }

	sprintf(cfilt,"%c", SURVEY_FILTERS[IFILT] );
	sprintf(CTIT,"%s %s vs. epoch", cfilt, CTAG_EPOCH[FLAG] );
	HID  = HOFF_EPOCH[FLAG] + IFILT + 1 ; // fortran-like index
	XMIN_F = 0.5 ;  XMAX_F = XMIN_F + (float)NBIN ;    
	hbook1_(&HID, CTIT, &NBIN, &XMIN_F, &XMAX_F, &WGT0_F, strlen(CTIT));

	firstObs = -9 ;
	for(obs=0; obs < NOBS; obs++ ) {
	  if ( SNLCPAK_OUTPUT.IFILT[FLAG][obs] == IFILT ) {
	    if ( firstObs < 0 ) { firstObs = obs; }
	    DATA4[obs]     = (float)ptrDATA[obs] ;
	    DATA_ERR4[obs] = (float)ptrERR[obs] ;
	  }
	  else {
	    DATA4[obs]    = NULL4[FLAG] ;
	    DATA_ERR4[obs] = 0.0 ;
	  }
	}  // obs loop


	if ( FLAG != SNLCPAK_EPFLAG_FITFUN ) { firstObs=0; }
	hpak_ ( &HID, &DATA4[firstObs] );
	hpake_( &HID, &DATA_ERR4[firstObs] ) ;


      } // USE if-block
	
    }  // IFILT loop    
  }  // FLAG loop



  // book Tobs for DATA
  NOBS = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FLUXDATA] ;
  HID = HID_DATA_TOBS ;
  XMIN_F = 0.5 ;  XMAX_F = XMIN_F + (float)NOBS ;
  sprintf(CTIT,"T?obs! vs. epoch for Data" );
  hbook1_(&HID, CTIT, &NOBS, &XMIN_F, &XMAX_F, &WGT0_F, strlen(CTIT));
  ptrTOBS = SNLCPAK_OUTPUT.TOBS[SNLCPAK_EPFLAG_FLUXDATA] ;
  for(obs=0; obs<NOBS; obs++ )  { DATA4[obs] = ptrTOBS[obs];  }
  hpak_ ( &HID, DATA4 );


  // ----------------------------------------
  // book info for  FITFUN

  NOBS_CHI2 = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_CHI2] ;
  OVFIT     = ( NOBS_FITFUN > 0 || NOBS_CHI2 > 0 ) ;  

  if ( OVFIT ) {

    FLAG = SNLCPAK_EPFLAG_FITFUN ;
    // list of Tobs in fine bins (smooth function)
    HID  = HID_FITFUN_TOBS ;
    XMIN_F = 0.5 ;  XMAX_F = XMIN_F + (float)NOBS_FITFUN ;    
    sprintf(CTIT,"T?obs! vs. epoch for FITFUN (fine bins)" );
    hbook1_(&HID, CTIT, &NOBS_FITFUN, &XMIN_F, &XMAX_F, &WGT0_F, strlen(CTIT));
    ptrTOBS = SNLCPAK_OUTPUT.TOBS[FLAG] ;
    for(obs=0; obs<NOBS_FITFUN; obs++ )  
      { DATA4[obs] = (float)ptrTOBS[obs];  }
    hpak_ ( &HID, DATA4 );
  } 

  
  // plot filter-dependent quantities
  for ( iplot=0; iplot < MXFLAG_SNLCPAK_FILTER; iplot++ ) {
    FLAG = 100 + iplot ;

    // plot NDOF and PKFLUX only for fit-model overlay ;
    // always book PKMJD
    if ( OVFIT == 0 && FLAG == SNLCPAK_BANDFLAG_NDOF   ) { continue ; }
    if ( OVFIT == 0 && FLAG == SNLCPAK_BANDFLAG_CHI2   ) { continue ; }
    if ( OVFIT == 0 && FLAG == SNLCPAK_BANDFLAG_PKFLUX ) { continue ; }

    HID  = HID_FILTER[iplot];   
    NBIN = NFILTDEF_SURVEY ;
    XMIN_F = 0.5;   XMAX_F = XMIN_F + (float)NBIN ;
    sprintf( CTIT, "%s  vs. filter", CTAG_FILTER[iplot] );
    hbook1_(&HID, CTIT, &NBIN, &XMIN_F, &XMAX_F, &WGT0_F, strlen(CTIT));

    for(IFILT=0; IFILT < NFILTDEF_SURVEY; IFILT++ ) {   
      X     = SNLCPAK_OUTPUT.FILTDATA[iplot][IFILT] ; 
      EX    = SNLCPAK_OUTPUT.FILTDATA_ERR[iplot][IFILT] ; 
      if ( FLAG == SNLCPAK_BANDFLAG_PKMJD && X > 20000. ) 
	{ X -= (double)53000.0; }      

      DATA4[IFILT]     = (float)X ;
      DATA_ERR4[IFILT] = (float)EX ;
    }

    hpak_  ( &HID, DATA4     );
    hpake_ ( &HID, DATA_ERR4 );
    
  } // iplot loop
  

  free(DATA4);
  free(DATA_ERR4);

  return ;

} // end of SNLCPAK_FILL_HBOOK()




// --------------------------------------------------------
void modify_text_4paw(char *HTEXT) {

  // to get nice symbols in paw replace the normal text substrings
  // with special paw symbols,
  //
  //   +-       -->   "A#
  //   mu,MU    -->   [m]
  //   ZPHOT    -->   z?ph!
  //   E-4       -->  10^-4!
  //   chi2     -->   [h]^2!
  //   <        -->   "L#
  //   etc ...
  //
  // If first char is #, then do nothing since it's a comment
  //
  // ----------- BEGIN -----------

  char HTEXT_ORIG[200], c0[2]; 

  sprintf(c0, "%c", HTEXT[0] ) ;
  if ( strcmp(c0,"#") == 0 ) { return ; }

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "mu", "[m]") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "MU", "[m]") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "x0", "x?0!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "x1", "x?1!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "DELTA", "[D]") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "ZPHOT", "z?ph!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "E-2", "10^[-]2!") ) ;
  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "E-3", "10^[-]3!") ) ;
  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "E-4", "10^[-]4!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "chi2", "[h]^2!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "AV", "A?V!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "RV", "R?V!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "BETA", "[b]") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "+-", "\"A#") ) ;
  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "+-", "\"A#") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "+_", "\"A#") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "<", "\"L#") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "OM", "[W]?M!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "OL", "[W]?[L]!") ) ;

  sprintf(HTEXT_ORIG,"%s", HTEXT);
  sprintf(HTEXT, "%s", replace_str(HTEXT_ORIG, "Or", "[W]?r!") ) ;

} // end of modify_text_4paw


int TABLEID_HBOOK(char *tableName) {

  // Created Jan 2017
  // Return integer ntuple ID for input string *tableName

  int ID ;

  if ( strcmp(tableName,TABLENAME_FITRES)==0 )
    { ID = TABLEID_FITRES; }
  else if ( strcmp(tableName,TABLENAME_SNANA)==0 )
    { ID = TABLEID_SNANA; }
  else
    { ID = atoi(tableName); }

  return(ID);

} // end TABLEID_HBOOK


