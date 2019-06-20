/***************************************************
   sntools_output_root.c
   Created Feb 23 2013 by R.Kessler
   Contains root-specific code for snana output.


             PROGRESS:
      ~~~~~~~~~~~~~~~~~~~~~~~~~~

  Mar 03 2013: root TTree working for SNANA and FITRES.
               These are the analogs of CWNTup IDs 7100 & 7788.

  Apr 2013: Add SNTABLE_DUMP_* functions.
            Called by new program sntable_dump.c .

  Apr 2013: modify SNTABLE_ADDCOL_ROOT to load one variable at a time
            instead of loading a string with variable names separated
            by a colon. This request is from ANL group.

  Apr 6 2014: 
    Fix bug that did not allow histograms to be updated for each SN.
    Problem was that NH_ROOT was reset to zero for each SN subdir
    to allow duplicate HID in each subdir. 
    Instead, fix SNHIST_INIT_ROOT to allow duplicate HIDs,
    but now there is no protection against duplicates in same sdir.
 
    In SNHIST_INIT_ROOT, abort if returned IHIST < 0.

    Add char optoin *COPT in OPEN_ROOTFILE to specify old or new file.

    Finally write code for SNHIST_RDBINS_ROOT and SNHIST_RDCONT_ROOT;
    needed by nearnbr_maxFoM.c


 May 11 2014
   SNTABLE_LOAD_ROOT -> SNTABLE_ADDCOL_ROOT  and major overhaul.

 Oct 2 2014: define MXEP_PER_FILT=400 (doubles old hard-wired value)

 Oct 26 2014: refactor table-reading.
 
 Feb 5 2016: write MJD into SNLCPAK

 Jan 14 2017: remove ISFILE_ROOT; move it to sntools_output.c

 Jun 20 2019: fix dump-output to include comma for csv format.

***************************************************/

#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "TKey.h" 
#include "TTreePlayer.h"
#include "TLeaf.h"
#include "TSQLRow.h" 
#include "TSQLResult.h" 

#include "TSystem.h"
#include "TFileMerger.h"

// ================================
//             Globals
// ================================

#define VBOSE_FLAG_ROOT  0   // 0=quiet   1=extra printf

#define MXEP_PER_FILT 400   // allow this many epochs per filter(data+FITFUN)

TFile  *TFILE_ROOT_NEW  ; // NEW file pointer
TFile  *TFILE_ROOT_READ ; // READ-mode file pointer

#define MXTABLE_ROOT  10   // max number of TTrees
int     NTABLE_ROOT ;
int     SNTABLE_ROOT_MAP[MXTABLE_ROOT];
TTree  *SNTable_ROOT[MXTABLE_ROOT] ;


#define MXHIST_ROOT 200
int   NH_ROOT[3];               // index = NDIM (1,2 => 1D,2D)
int   HROOTMAP_ID[3][MXHIST_ROOT] ;     // HID vs. sparse index
int   HROOTMAP_DIR[3][MXHIST_ROOT] ;    // dir index vs. 
char  HROOT_DIRNAMES[MXHIST_ROOT][40] ; // store dirname for each histo


int ROOT_OPENFLAG ;

TH1F *H1_ROOT[MXHIST_ROOT] ;
TH2F *H2_ROOT[MXHIST_ROOT] ;


#define MONSN_DIRNAME  "MONSN"  // contains subdir(s) for each SN with plots
TDirectory *TD_MONSN, *TD_MONSN_SUBDIR, *TD_TEMP ;

// define compact structre to fill for tree.

struct SNLCPAK_TREEDATA {
  int   MXEP  ; // max epochs to allocate
  int   NDATA, NFITFUN, NTOT ;

  char TREENAMES[200]; // space-separated list of SNLCPAK tree names

  // pointers to externally-generated text strings to display avove plots
  char *DISPLAYTEXT[MXFIT_PER_SN];

  // quantities vs. epoch
  float  *TOBS_F, *FLUX_F, *FLUXERR_F, *CHISQ_F ;
  float  *SIMFLUX_F ; // May 2016
  double *MJD_D ;  // added Feb 2016
  int   *IFILTOBS, *IREJECT, *IFLAGDATA;
  char  **CFILTOBS ;

  // quantities vs. filter index  (0 to NSURVEY_FILTERS-1)
  int     NFILT ;
  double *BAND_PKMJD, *BAND_PKMJDERR ;
  float  *BAND_PKFLUX_F,  *BAND_CHI2_F ;
  int    *BAND_NDOF ;

} SNLCPAK_TREEDATA ;


TTree *SNLCPAK_TREE[MXFIT_PER_SN] ;


// define info for reading back a root tree
struct TREE_INFO_READ {
  //  TFile       *file ;
  TTree       *tree ;
  TObjArray   *obj;
  TIterator   *iter;
  TLeaf       *leaf;

  TTreePlayer *player;
  TSQLResult  *SQL_Result;
  TSQLRow     *SQL_Row;

  char VARNAME[MXVAR_TABLE][60];
  int  ICAST[MXVAR_TABLE];

} TREE_INFO_READ ;



// ================================
// function declarations
// ================================

#ifdef __cplusplus
extern"C" {
#endif

  void OPEN_ROOTFILE(char *FILENAME, char *COPT, int *IERR);
  void CLOSE_ROOTFILE(char *FILENAME, int OPENFLAG);

  void SNTABLE_CREATE_ROOT(int IDTABLE, char *name);

  void SNTABLE_FILL_ROOT(int IDTABLE); 

  void SNTABLE_ADDCOL_ROOT(int IDTABLE, void *PTRVAR, 
			   SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF); 

  int  ITABLE_ROOT ( int IDTABLE ) ;

  void SNTABLE_LIST_ROOT(char *fileName); // list all trees

  void SNTABLE_DUMP_VARNAMES_ROOT(char *fileName, char *treeName);

  void SNTABLE_LIST_ROOT(char *fileName) ;
  int  SNTABLE_NEVT_ROOT(char *FILENAME, char *treeName) ;

  int  SNTABLE_READPREP_ROOT(char *TABLENAME);
  int  SNTABLE_READ_EXEC_ROOT(void);
  int  sntable_read_exec_root(int ievt_min, int ievt_max);

  int select_outlier_root(double *DARRAY);

  void SNLCPAK_SURVEY_ROOT(void);
  void SNLCPAK_makeTree_ROOT(int ipak) ;
  void SNLCPAK_MALLOC_TREEDATA(void);
  void SNLCPAK_Branch_DISPLAYTEXT(int ipak);
  void SNLCPAK_Branch_lightCurve(int ipak);    // light curve info
  void SNLCPAK_Branch_filter(int ipak);    // filter info
  void SNLCPAK_FILL_TREEDATA(void) ;
  void SNLCPAK_FILL_ROOT(void) ;

  void SNHIST_INIT_ROOT( int NDIM, int ID, char *TITLE, 
			 int *NBIN, double *XMIN, double *XMAX ) ; 
  void SNHIST_FILL_ROOT(int NDIM, int ID, double *VALUE, double WGT); 
  int  IHIST_ROOT ( int NDIM, int IDHIST ) ;
  int  ITABLE_ROOT ( int IDTABLE ) ;
  void SNHIST_RDBINS_ROOT(int ID, char *CTIT, int *NBIN,
			   double *XMIN, double *XMAX);
  void SNHIST_RDCONT_ROOT(int ID, int NRDBIN, double *CONTENTS);

  void MAKEDIR_ROOT(char *CCID);
  void CDTOPDIR_ROOT(void);

  void modify_text_4root(char *Text) ;

  void get_tree(char *fileName, char *treeName);
  int  prepTreeString(char *treeName, char *string) ;

  void MERGE_ROOT(int NFILE, char **INFILES, char *OUTFILE);

#ifdef __cplusplus
}
#endif


// ==========================================
//
//          OPEN/CLOSE FUNCITONS
//
// ==========================================

// ------------------------------------------------
void OPEN_ROOTFILE(char *FILENAME, char *COPT, int *IERR) {

  // Feb 2013: wrapper to initialize file for root.

  int LVBOSE = 1 ;
  char cstat[12], banner[80] ;

  // check for quiet option 
  if ( strchr(COPT,'Q') != NULL )  { LVBOSE = 0 ; }
  if ( strchr(COPT,'q') != NULL )  { LVBOSE = 0 ; }

  if ( strchr(COPT,'N') != NULL )  { 
    sprintf(cstat,"RECREATE"); //   new file; clobber old one 
    ROOT_OPENFLAG    = OPENFLAG_NEW ;
    TFILE_ROOT_NEW   = new TFile(FILENAME,cstat);
    SNLCPAK_USE_ROOT = 1 ;  

    NH_ROOT[1]  = 0 ;  // 1D histos
    NH_ROOT[2]  = 0 ;  // 2D histos
    NTABLE_ROOT = 0;   // tables (Trees)
  }
  else  { 
    sprintf(cstat,"READ");  // read only
    ROOT_OPENFLAG   = OPENFLAG_READ ;
    TFILE_ROOT_READ = new TFile(FILENAME,cstat);
  }


  if ( LVBOSE ) {
    sprintf(banner,"OPEN_ROOTFILE: use  %s  mode:", cstat ) ;
    print_banner(banner);
    printf("   Opened %s\n", FILENAME); fflush(stdout);
  }


  FIRSTCALL_TABLEFILE_OPEN[IFILETYPE_ROOT] = 0;

  *IERR = 0 ; // return SUCCESS flag


} // end of OPEN_ROOTFILE



// ===================================================
void CLOSE_ROOTFILE(char *FILENAME, int OPENFLAG) {

  if ( ROOT_OPENFLAG == OPENFLAG_NEW ) {
    printf("   Write root file %s \n", FILENAME); 
    TFILE_ROOT_NEW->Write() ;

    printf("   Close root file %s \n", FILENAME); 
    TFILE_ROOT_NEW->Close() ;
  }
  else {
    printf("   Close root file %s \n", FILENAME ); 
    TFILE_ROOT_READ->Close() ;
  }

  fflush(stdout);

} // end of CLOSE_ROOTFILE


// ==========================================
//
//          SNTABLE FUNCTIONS
//
// ==========================================


// ------------------------------------------------
void SNTABLE_CREATE_ROOT(int IDTABLE, char *name) {

  // initialize table for root by defining TTree.
  // This is an ntuple-like structure.

  char title[100];
  int  i ;
  char fnam[] = "SNTABLE_CREATE_ROOT" ;

  // --------------- BEGIN ---------------

  if ( IDTABLE == 1600 ) { return; } // this one is only for hbook

  NTABLE_ROOT++ ; // keep track of each table.

  if ( NTABLE_ROOT >= MXTABLE_ROOT ) {
    sprintf(MSGERR1,"NTABLE_ROOT = %d exceeds bound.", NTABLE_ROOT);
    sprintf(MSGERR2,"See   #define MXTABLE_ROOT");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }


  printf("\n  %s: Book %s tree \n",  fnam, name) ; fflush(stdout);

  i = NTABLE_ROOT - 1;
  sprintf(title,"%s table, ID = %d", name, IDTABLE);  
  SNTable_ROOT[i]     = new TTree(name,title) ;
  SNTABLE_ROOT_MAP[i] = IDTABLE ;
 
} // end of  SNTABLE_CREATE_ROOT" 


// ===============================================================
void SNTABLE_ADDCOL_ROOT(int IDTABLE, void *PTRVAR, 
			 SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF) {

  // Created Feb 2013 by R.Kessler
  // Add tree column(s) based on variables in input struct ADDCOL_VARDEF.
  

  int ITABLE, ivar, ICAST, IVEC, ISIZE, LDMP ;

  char 
    CCAST[4]
    ,VARNAME[80]
    ,BRNAME[80]
    ,treeVarName[80]
    ,varDum[80]
    ,fnam[] = "SNTABLE_ADDCOL_ROOT" 
      ;

  long long int *ptr_L;
  double *ptr_D ;
  float  *ptr_F ;
  int    *ptr_I ;

  // -------------- BEGIN ------------------

  if ( IDTABLE == 1600 ) { return; } // this one is only for hbook

  ITABLE = ITABLE_ROOT(IDTABLE) ; // sparse table index


  LDMP = 0 ;
  if ( LDMP ) {
    printf(" xxx ------------------------------ \n");
    printf(" xxx VARLIST_ORIG = '%s' \n",
	   ADDCOL_VARDEF->VARLIST_ORIG); fflush(stdout);
  }

  // load table
  for(ivar=0; ivar < ADDCOL_VARDEF->NVAR; ivar++ ) {

    ICAST   = ADDCOL_VARDEF->ICAST[ivar] ;
    IVEC    = ADDCOL_VARDEF->VECTOR_FLAG[ivar] ;
    ISIZE   = ADDCOL_VARDEF->ISIZE[ivar] ;
    sprintf(CCAST,  "%s", ADDCOL_VARDEF->CCAST[ivar]   ) ;
    sprintf(VARNAME,"%s", ADDCOL_VARDEF->VARNAME[ivar] ) ;
    sprintf(BRNAME, "%s", ADDCOL_VARDEF->VARNAME[ivar]    ) ;

    if ( IVEC == 1 ) {
      // define integer  that sets vector size
      sprintf(treeVarName,"%s/I", VARNAME );
    }
    else if ( IVEC == 2 ) {
      // replace () with [] ; treeVarName has the form 'MJD(NFIT)'
      sprintf(varDum, "%s/%s", VARNAME, CCAST); 
      sprintf(treeVarName, "%s", replace_str(varDum, "(", "[") ) ;

      sprintf(varDum, "%s", treeVarName); 
      sprintf(treeVarName, "%s", replace_str(varDum, ")", "]") ) ;
    }
    else {
      sprintf(treeVarName,"%s/%s", VARNAME, CCAST ) ;
    }

    if ( LDMP ) {
      printf(" xxx treeVarName[%d] = '%s' \n", ivar, treeVarName);
      fflush(stdout);
    }

    if ( ICAST == ICAST_D ) { 
      ptr_D = (double*)PTRVAR ;
      SNTable_ROOT[ITABLE]->Branch( BRNAME, (ptr_D+ivar), treeVarName );
    }
    else if ( ICAST == ICAST_L ) { 
      ptr_L = (long long int*)PTRVAR ;
      SNTable_ROOT[ITABLE]->Branch( BRNAME, (ptr_L+ivar), treeVarName );
    }
    else if ( ICAST == ICAST_F ) {
      ptr_F = (float*)PTRVAR ; 
      SNTable_ROOT[ITABLE]->Branch( BRNAME, (ptr_F+ivar), treeVarName );
    }
    else if ( ICAST == ICAST_I ) { 
      ptr_I = (int*)PTRVAR ; 
      SNTable_ROOT[ITABLE]->Branch( BRNAME, (ptr_I+ivar), treeVarName );
    }
    else if ( ICAST == ICAST_C )  { 
      SNTable_ROOT[ITABLE]->Branch( BRNAME, PTRVAR, treeVarName );
    }
    else {
      sprintf(MSGERR1,"Invalid ICAST=%d for VARNAME=", ICAST );
      sprintf(MSGERR2,"%s", VARNAME);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }
  }

}  // end of SNTABLE_ADDCOL_ROOT


// ------------------------------------------------
void SNTABLE_FILL_ROOT(int IDTABLE) {

  //  char fnam[] = "SNTABLE_FILL_ROOT" ;

  if ( IDTABLE == 1600 ) { return; } // this one is only for hbook
  int ITABLE = ITABLE_ROOT(IDTABLE); // sparse table index
  SNTable_ROOT[ITABLE]->Fill();

} // end of SNTABLE_FILL


int ITABLE_ROOT ( int IDTABLE ) {
  // return sparse index for IDTABLE

  int i ;
  char fnam[] = "ITABLE_ROOT" ;

  for ( i = 0; i < NTABLE_ROOT ; i++ ) {
    if ( IDTABLE == SNTABLE_ROOT_MAP[i] ) { return i ; }
  }

  // if we get here, then abort
  sprintf(MSGERR1,"Could not find sparse index for IDTABLE=%d", IDTABLE);
  sprintf(MSGERR2,"NTABLE_ROOT = %d", NTABLE_ROOT);
  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 

  return(-9); // should never get here

} // end of ITABLE_ROOT


// ===================================
void SNTABLE_LIST_ROOT(char *fileName) {

  // List all trees in *fileName

  TFile *f ;
  TKey  *key;
  //  char  fnam[] = "SNTABLE_LIST_ROOT" ;

  // --------- BEGIN -------------

  printf("\n List trees in %s \n", fileName);  fflush(stdout);

  f = new TFile(fileName,"READ");
  f->ls();  //Lists all objects in file
  
  // f->GetListOfKeys()->Print() ;
  TIter next(f->GetListOfKeys());

  while ((key=(TKey*)next())) {

    printf("key: %-14s points to object of class %-20s .\n"
	   ,key->GetName()
           ,key->GetClassName()  );

    /*
    if(strstr("TTree",key->GetClassName()) != NULL) 
      { printf(" Hello.\n"); } */
  }


} // end of SNTABLE_LIST_ROOT


// =================================================
int SNTABLE_NEVT_ROOT(char *FILENAME, char *treeName) {

  // Apri 2014: 
  // return number of rows in tree.
  // If input FILENAME is blank (LENF=0), assume file is already open.

  int Nrow, LENF, IERR;
  char blank[] = " " ;
  //  char fnam[] = "SNTABLE_NEVT_ROOT" ;
  // ------------ BEGIN --------------

  Nrow = 0 ;
  LENF = strlen(FILENAME);

  if ( LENF > 0 ) { 
    OPEN_ROOTFILE(FILENAME, blank, &IERR) ; 
    if ( IERR != 0 ) { return Nrow ; }
  }

  get_tree(FILENAME, treeName);

  char scanString[MXCHAR_VARLIST];
  scanString[0] = 0 ;
  TREE_INFO_READ.SQL_Result = TREE_INFO_READ.tree->Query(scanString, "");
  Nrow   = TREE_INFO_READ.SQL_Result->GetRowCount();

  if ( LENF > 0 ) { CLOSE_ROOTFILE(FILENAME,OPENFLAG_READ); }

  return Nrow ;

} // end of SNTABLE_NEVT_ROOT


// =================================================
int  SNTABLE_READPREP_ROOT(char *treeName) {

  // Jan 30 2016: add return(NVAR). Bug found by -O1 optimization.


  int NVAR, ICAST ;
  char leafName[80], leafType[80];
  char noFile[] = "" ; // because file should already be opened
  //  char fnam[] = "SNTABLE_READPREP_ROOT" ;

  // ------------ BEGIN -------------

  get_tree(noFile, treeName);

  NVAR = 0;

  while( (TREE_INFO_READ.leaf =(TLeaf *)TREE_INFO_READ.iter->Next()) ) { 
   
    sprintf(leafName, "%s", TREE_INFO_READ.leaf->GetName());      // varName
    sprintf(leafType, "%s", TREE_INFO_READ.leaf->GetTypeName());  // cast

    // determine cast
    ICAST = 0 ;
    if ( strcmp(leafType,"Float_t") == 0 ) 
      { ICAST = ICAST_F ; }
    else if ( strcmp(leafType,"Double_t") == 0 ) 
      { ICAST = ICAST_D ; }
    else if ( strcmp(leafType,"Long_t") == 0 ) 
      { ICAST = ICAST_L ; }
    else if ( strcmp(leafType,"Long64_t") == 0 ) 
      { ICAST = ICAST_L ; }
    else if ( strcmp(leafType,"Int_t") == 0 ) 
      { ICAST = ICAST_I ; }
    else if ( strcmp(leafType,"Char_t") == 0 ) 
      { ICAST = ICAST_C ; } 

    sprintf(TREE_INFO_READ.VARNAME[NVAR], "%s", leafName);
    TREE_INFO_READ.ICAST[NVAR] = ICAST ;
    READTABLE_POINTERS.ICAST_READ[NVAR] = ICAST; // Jan 2017
    sprintf(READTABLE_POINTERS.VARNAME[NVAR], "%s", leafName);

    NVAR++ ;
  }

  READTABLE_POINTERS.NVAR_TOT = NVAR ;

  return(NVAR);

} // end of   SNTABLE_READPREP_ROOT

// ==============================================   
int SNTABLE_READ_EXEC_ROOT(void) {

  int  NROW_TOT, Nrow, irow, IROW_MIN, IROW_MAX ;
  int  NROW_READ_ITER = 20000 ;
  char BLANKFILE[] = "";
  char fnam[] = "SNTABLE_READ_EXEC_ROOT" ;
  
  // ----------- BEGIN ----------


  NROW_TOT = SNTABLE_NEVT_ROOT(BLANKFILE, READTABLE_POINTERS.TABLENAME) ;
  Nrow     = 0 ;

  for(irow=0; irow < NROW_TOT; irow+= NROW_READ_ITER ) {
    IROW_MIN = irow;
    IROW_MAX = irow + NROW_READ_ITER-1;
    if ( IROW_MAX >= NROW_TOT-1 ) { IROW_MAX = NROW_TOT-1; }
    Nrow += sntable_read_exec_root(IROW_MIN,IROW_MAX);

    if ( fmodf( (float)(Nrow), 100000. ) == 0 && Nrow > 0 )  {
      printf("\t Finished reading table row %d of %d \n", 
	     Nrow, NROW_TOT ); fflush(stdout);
    }

  } // end irow loop

  return(NROW_TOT) ;

} // end SNTABLE_READ_EXEC_ROOT

// ========================================
int sntable_read_exec_root(int IROW_MIN, int IROW_MAX) {

  // Oct 2014:                                                                
  // Execute read over all tree rows and fill
  // pointers passed previously to SNTABLE_READPREP_VARDEF.
  // Functions returns number of rows read.                                   
  //
  // Mar 3 2017: 
  //   refactored to read limited irow range and delete tree.
  //   Goal is to limit instantaneous memory usage.
  //  
  // Mar 11 2019: query args -> long long (intead of just long)
  //
  // Jun 20 2019: include SEPKEY for csv output format.
  //
  
  int NVAR_READ_TOT = READTABLE_POINTERS.NVAR_READ ;
  int FIRST = ( IROW_MIN == 0 ) ;

  int  LREAD, LDUMP, LOUTLIER, DO_DUMP ;
  int  LEN_SCANSTRING, IVAR_READ, IVAR_TOT ;
  char *scanString, *tmpString, *ptrVar, *treeName, *optString ;
  char BLANKFILE[] = "";
  char   *SEPKEY = READTABLE_POINTERS.SEPKEY_DUMP ;
  char fnam[] = "sntable_read_exec_root";
  FILE *FP_DUMP ;
  // ------------ BEGIN -----------                                           

  /*
  printf(" xxx %s: IROW_MIN/MAX = %d / %d  \n",
	 fnam, IROW_MIN, IROW_MAX ); fflush(stdout);
  */
	 
  treeName      = READTABLE_POINTERS.TABLENAME ;
  FP_DUMP       = READTABLE_POINTERS.FP_DUMP ;

  LREAD = LDUMP = LOUTLIER = 0 ;

  if ( OUTLIER_INFO.USEFLAG ==1 ) {
    LOUTLIER = LDUMP = 1;
    if(FIRST) { printf("\t Dump fit-outliers for each event. \n"); }
  }
  else if ( OUTLIER_INFO.USEFLAG == 2 ) {
    LOUTLIER = LDUMP = 1;
    if(FIRST) { printf("\t Dump all observations for each event. \n"); }
  }
  else if ( FP_DUMP ) {
    LDUMP = 1;
    if(FIRST) { printf("\t Dump SNTABLE values for each event. \n"); }
  }
  else if ( READTABLE_POINTERS.MXLEN > 1 ) {
    LREAD = 1;
    if (FIRST) { printf("\t Read/store table values for each tree row.\n"); }
  }
  else {
    sprintf(MSGERR1,"Cannot figure out what to do.");
    sprintf(MSGERR2,"Check calls to SNTABLE_READ_XXX");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);   
  }

  fflush(stdout);

  // --------------------------------------------
  // read tree from file that is already open

  get_tree(BLANKFILE, treeName); 

  // allocate local memory
  LEN_SCANSTRING = MXCHAR_VARLIST ;
  scanString     = (char*) malloc( LEN_SCANSTRING * sizeof(char) );
  tmpString      = (char*) malloc( LEN_SCANSTRING * sizeof(char) );
  optString      = (char*) malloc( 100 * sizeof(char) );
  scanString[0]  = 0 ;
  tmpString[0]   = 0 ;
  sprintf(optString,"precision=16");

  // construct colon-separated "scan" string and make sure 
  // that each variable exists.

  for(IVAR_READ=0; IVAR_READ < NVAR_READ_TOT ; IVAR_READ++ ) {

    IVAR_TOT = READTABLE_POINTERS.PTRINDEX[IVAR_READ] ;
    ptrVar   = READTABLE_POINTERS.VARNAME[IVAR_TOT] ;

    if ( IVAR_READ == 0  )
      { sprintf(scanString,"%s", ptrVar); }
    else
      { sprintf(scanString,"%s:%s", scanString, ptrVar); }

    // make sure VARNAME exists.
    if ( TREE_INFO_READ.tree->GetLeaf(ptrVar) == 0 ) {
      sprintf(MSGERR1,
	      "VARNAME = '%s' does not exist in tree = '%s'",ptrVar, treeName);
      sprintf(MSGERR2,
	      "Check valid list by removing '--VARNAMES <varNames>' arg,");
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }
  } // IVAR_READ


  //  TREE_INFO_READ.tree->Scan(scanString) ;

  int    Nfield, Nrow, irow, i, istat, ICAST ;
  char   LINE[MXCHAR_VARLIST] ;
  char   valStore[MXVAR_TABLE][40];
  double DVAL, DARRAY[MXVAR_TABLE] ;

  long long int IROW_START = IROW_MIN;
  long long int NROW_READ  = IROW_MAX - IROW_MIN + 1;
  
  //  TREE_INFO_READ.SQL_Result = TREE_INFO_READ.tree->Query(scanString, "");

  
  // this read is the slow part.
  // 1st arg is varList, 2nd arg is selection, 3rd arg is option
  TREE_INFO_READ.SQL_Result = 
    TREE_INFO_READ.tree->Query(scanString, "", optString, 
			       NROW_READ, IROW_START);

  Nfield = TREE_INFO_READ.SQL_Result->GetFieldCount();
  Nrow   = TREE_INFO_READ.SQL_Result->GetRowCount();
  irow   = 0 ;
  LINE[0] = 0;

  while( (TREE_INFO_READ.SQL_Row =
	  (TSQLRow *)TREE_INFO_READ.SQL_Result->Next()) ) {

    sprintf(LINE, "%s", READTABLE_POINTERS.LINEKEY_DUMP ) ; 

    for (i=0; i < Nfield; i++) {  // ivar_read loop

      sprintf(tmpString,"%s ", TREE_INFO_READ.SQL_Row->GetField(i)); 

      if ( strlen(tmpString) > 1 )  { 
	istat = prepTreeString(treeName, tmpString ) ; 
	sprintf(valStore[i], "%s", tmpString); 
      }
      else  { 
	// for string values we get a blank after the 1st epoch (??)
	// so paste 1st epoch value back here. This is a very strange
	// feature ??
	sprintf(tmpString,"%s", valStore[i] );
	istat = 0 ; 	
      }

      IVAR_TOT = READTABLE_POINTERS.PTRINDEX[i] ;
      ICAST    = TREE_INFO_READ.ICAST[IVAR_TOT] ;
      DVAL     = -99999. ;

      // load local DVAL 
      if ( ICAST != ICAST_C ) {	sscanf( tmpString, "%le", &DVAL );  }

      DARRAY[i] = DVAL ;

      if( LREAD )   { 
	load_READTABLE_POINTER(irow+IROW_MIN, i, DVAL, tmpString); 
      }
      else if ( LDUMP )  { 
	if ( ICAST == ICAST_C ) 
	  { sprintf(LINE,"%s %s", LINE, tmpString); }
	else
	  { load_DUMPLINE(LINE,DVAL); }        // .xyz
      }

      if ( i < Nfield-1 ) { strcat(LINE,SEPKEY); }  // June 20 2019

    } // end of i loop over Nfield


    DO_DUMP = 0;
    if ( LDUMP    ) { DO_DUMP = 1; }
    if ( LOUTLIER ) { DO_DUMP = select_outlier_root(DARRAY) ; }
    
    if ( DO_DUMP )  { 
      fprintf(FP_DUMP,"%s\n", LINE ); 
      fflush(FP_DUMP); 
    }

    irow++ ;
    
  } // end while loop

  TREE_INFO_READ.tree->Delete("");           // Mar 3 2017
  TREE_INFO_READ.SQL_Result->Delete("") ;

  free(scanString);
  free(tmpString);

  return(Nrow) ;

} // end of  sntable_read_exec_root



// =================================================
void SNTABLE_DUMP_VARNAMES_ROOT(char *fileName, char *treeName) {

  // Apr 20 2013: utility to dump name of each variable in 
  //              file *fileName and tree *treeName

  //  char fnam[] = "SNTABLE_DUMP_VARNAMES_ROOT" ;

  char leafName[100];
  char leafType[100];
  int NVAR_LIST ;

  // -------------- BEGIN -------------

  get_tree(fileName, treeName);

  NVAR_LIST = 0;
  printf("\n List of SNANA tree-variables: \n" );

  while( (TREE_INFO_READ.leaf = (TLeaf *)TREE_INFO_READ.iter->Next()) ) {
    //    printf("%s \n", tbr->GetName());
    sprintf(leafName, "%s", TREE_INFO_READ.leaf->GetName());
    sprintf(leafType, "%s", TREE_INFO_READ.leaf->GetTypeName());
    printf("%-30s     (%s) \n", leafName, leafType);
    fflush(stdout);
    NVAR_LIST++ ;
  }

  printf("\n Done listing %d variables in %s tree. \n",
	 NVAR_LIST, treeName );
  fflush(stdout);
 

} // end of SNTABLE_DUMP_VARNAMES_ROOT


// =============================================================
int select_outlier_root(double *DARRAY) {

  // Created Oct 2014
  // Input DARRAY is an array of variable values for current row.
  // Examine values to make CHI2FLUX cut, and to increment 
  // statistics vs. IFILTOBS.
  //
  // Function returns 1 of OUTLIER cut is satisfied; 0 otherwise.
  //

  int   IVAR_CHI2, IVAR_IFILT, IFILT ;
  float CHI2, CHI2MIN, CHI2MAX ;
  char fnam[] = "select_outlier_root" ;


  // -------------- BEGIN ----------

  CHI2MIN     = OUTLIER_INFO.CUTWIN_CHI2FLUX[0] ;
  CHI2MAX     = OUTLIER_INFO.CUTWIN_CHI2FLUX[1] ;
  IVAR_IFILT  = OUTLIER_INFO.IVAR[INDX_OUTLIER_IFILT];
  IVAR_CHI2   = OUTLIER_INFO.IVAR[INDX_OUTLIER_CHI2]; 

  IFILT = (int)DARRAY[IVAR_IFILT] ;
  CHI2  = (float)DARRAY[IVAR_CHI2];

  // sanity checks
  if ( IFILT < 0  ) {
    sprintf(MSGERR1,"Could not find %s = %d", 
	    OUTLIER_VARNAME_IFILT, IFILT );
    sprintf(MSGERR2,"Check variable list");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  if ( CHI2 < 0.  ) {
    sprintf(MSGERR1,"Could not find %s = %f .", 
	    OUTLIER_VARNAME_CHI2, CHI2 );
    sprintf(MSGERR2,"Check variable list");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }


  // apply CHI2 cut
  OUTLIER_INFO.NEP_TOT[0]++ ;
  OUTLIER_INFO.NEP_TOT[IFILT]++ ;

  if ( CHI2 >= CHI2MIN && CHI2 <= CHI2MAX ) {
    OUTLIER_INFO.NEP_SELECT[0]++ ;
    OUTLIER_INFO.NEP_SELECT[IFILT]++ ;
    return 1 ;
  }
  else {
    return 0;
  }

} // end of  select_outlier_root


// =====================================
void get_tree(char *fileName, char *treeName) {

  // open root file and get tree
  // if fileName is blank, then file is already opened.

  char fnam[] = "get_tree" ;

  // ---------------- BEGIN -------------

  //  TREE_INFO_READ.file     = new TFile(fileName, "READ");
  if ( strlen(fileName) > 0 )  { 
    TFILE_ROOT_READ    = new TFile(fileName, "READ");

    if ( TFILE_ROOT_READ->IsZombie()  ) {
      sprintf(MSGERR1,"Root file = '%s' does not exist or is not valid.", 
	      fileName);
      sprintf(MSGERR2,"Try another root file.");
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
    }
  }

  TREE_INFO_READ.tree      = (TTree *)TFILE_ROOT_READ->Get(treeName);
  if ( TREE_INFO_READ.tree == 0 ) {
    sprintf(MSGERR1,"Tree = '%s' does not exist.", treeName);
    sprintf(MSGERR2,"Try another tree name.");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  TREE_INFO_READ.obj      = TREE_INFO_READ.tree->GetListOfLeaves();
  TREE_INFO_READ.iter     = TREE_INFO_READ.obj->MakeIterator();


} // end of get_tree


// ==============================================
int prepTreeString(char *treeName, char *string) {

  // if treeName is one of the INFO_TABLES (GLOBAL or SNLCPAK_*), 
  // then do NOT modify string, and return integer flag 1.
  // For any other treeName, remove spaces from string
  // and return integer flag 0.

#define NTABLE_INFO 2
  int istat, i, itab, LENNAM ; 
  char c1[2];  
  char INFO_TABLES[NTABLE_INFO][20] = { "GLOBAL", "SNLCPAK" } ;
  char TREETMP[100];

  istat = 0; // default is SN-dependent tree

  // get tree7 = first 7 chars
  for (itab=0; itab < NTABLE_INFO; itab++ ) {
    LENNAM = strlen(INFO_TABLES[itab]);
    TREETMP[0] = 0 ;
    
    if ( (int)strlen(treeName) >= LENNAM ) {
      for(i=0; i < LENNAM; i++ ) {
	sprintf(c1, "%c", treeName[i] );
	sprintf(TREETMP, "%s%s", TREETMP, c1);
      }
    }     

    if ( strcmp(TREETMP,INFO_TABLES[itab] ) == 0 )   { istat = 1 ; }
  }  // end itab


  if( istat == 0 && strlen(string)>0 ) 
    { trim_blank_spaces(string); }
 
  return istat ;

} // end of prepTreeString

// ==========================================
//
//          SNHIST FUNCTIONS
//
// ==========================================


void SNHIST_INIT_ROOT( int NDIM, int ID, char *TITLE, 
		       int *NBIN, double *XMIN, double *XMAX ) {

  // Created Mar 3, 2013 by R.Kessler
  //
  // Initialize 1D or 2D histogram
  // NDIM = 1 or 2 (dimension)
  // ID   = integer ID
  // TITLE = charString title
  // NBIN[idim] = Number of bins for each dimension
  // XMIN[idim] = min axis value for each dim.
  // XMAX[idim] = max axis value for each dim.
  //

  int N, indx;
  char 
    HNAM[20]
    ,fnam[] = "SNHIST_INIT_ROOT" 
    ;

  // ---------------- BEGIN ----------------


  // if this histogram is already used, then create a new one
  // with the same name. Allows the same histo name in each subdir.
  // Warning: no protection for duplicates in same directory !
  indx = IHIST_ROOT ( NDIM, ID ) ;
  if ( indx >= 0 ) { goto BOOKIT ; }

  NH_ROOT[NDIM]++ ;      N = NH_ROOT[NDIM] ;

  if ( N >= MXHIST_ROOT ) {
    sprintf(MSGERR1,"Number of %dD histos exceeds bound of MXHIST_ROOT=%d", 
	    NDIM, MXHIST_ROOT);
    sprintf(MSGERR2,"ID=%d  TITLE='%s'", ID, TITLE);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  indx = N-1 ;  
  HROOTMAP_ID[NDIM][indx]  = ID ;
  HROOTMAP_DIR[NDIM][indx] = 444 ; // ???

  //  sprintf(HROOT_DIRNAMES[next],"%s", "BLA" );


 BOOKIT:
  sprintf(HNAM,"H%5.5d", ID);

  //  printf(" xxx init '%s' TITLE='%s' \n", HNAM, TITLE); fflush(stdout);

  if ( NDIM == 1 ) {
    H1_ROOT[indx] = new TH1F( HNAM, TITLE, 
			      NBIN[0], (float)XMIN[0], (float)XMAX[0] ); 
  }
  else if ( NDIM == 2 ) {
    H2_ROOT[indx] = new TH2F( HNAM, TITLE
			      ,NBIN[0], (float)XMIN[0], (float)XMAX[0]
			      ,NBIN[1], (float)XMIN[1], (float)XMAX[1]
			      );  
  }
  else {
    sprintf(MSGERR1,"Invalid NDIM = %d  for ID=%d", NDIM, ID); 
    sprintf(MSGERR2,"TITLE  = '%s'", TITLE);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

} // end of SNHIST_INIT_ROOT


// ------------------------------------------------
void SNHIST_FILL_ROOT(int NDIM, int ID, double *VALUE, double WGT) {

  // Created Mar 3, 2013 by R.Kessler

  int IHIST, idim ;
  float VAL_F[2], WGT_F;
  char fnam[] = "SNHIST_FILL_ROOT" ;

  // --------- BEGIN ---------

  IHIST = IHIST_ROOT(NDIM, ID); // sparse index

  if ( IHIST < 0 ) {
    
    printf("\n\t\t PRE-ABORT DUMP: \n\n" );
    TFILE_ROOT_NEW->pwd();  
    TFILE_ROOT_NEW->ls();  
    
    sprintf(MSGERR1,"Could not find HID=%d  (IHIST=%d)", ID, IHIST);
    sprintf(MSGERR2,"NH_ROOT[NDIM=%d] = %d", NDIM, NH_ROOT[NDIM] );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }


  WGT_F = (float)WGT ;
  for (idim=0; idim < NDIM; idim++ ) 
    { VAL_F[idim] = (float)VALUE[idim] ;  }

  if ( NDIM == 1 ) 
    { H1_ROOT[IHIST]->Fill(VAL_F[0], WGT_F); }
  else if ( NDIM == 2 ) 
    { H2_ROOT[IHIST]->Fill(VAL_F[0], VAL_F[1], WGT_F); }

} // end of SNHIST_FILL_ROOT


// ------------------------------------------------
int IHIST_ROOT ( int NDIM, int IDHIST ) {

  // return sparse index for IDHIST

  int i ;
  //  char fnam[] = "IHIST_ROOT" ;

  for ( i = 0; i < NH_ROOT[NDIM] ; i++ ) {
    if ( IDHIST == HROOTMAP_ID[NDIM][i] ) { return i ; }
  }

  return -9 ;


} // end of IHIST_ROOT



// =======================================================
void SNHIST_RDBINS_ROOT(int NDIM, int ID, char *CTIT, int *NBIN,
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

  double XBIN[2] ;
  char HNAM[12] ;
  char fnam[] = "SNHIST_RDBINS_ROOT" ;
  int  LDMP = 0 ;
  
  // -------------- BEGIN -------------

  sprintf(HNAM,"H%5.5d", ID);

  if (LDMP) { 
    printf(" xxx -------------------------------------- \n");
    printf(" xxx %s  for  HNAM='%s'  (NDIM=%d)\n", fnam, HNAM, NDIM); 
  }

  XBIN[0] = XBIN[1] = -999.0 ;

  if ( NDIM == 1 ) { 
    TH1F *H = (TH1F*)TFILE_ROOT_READ->Get(HNAM) ; 
    NBIN[0]  = H->GetXaxis()->GetNbins();
    XBIN[0]  = (double)H->GetXaxis()->GetBinWidth(1);
    XMIN[0]  = (double)H->GetXaxis()->GetBinLowEdge(1);
    sprintf(CTIT, "%s", H->GetTitle() );
  }
  else  { 
    TH2F *H = (TH2F*)TFILE_ROOT_READ->Get(HNAM) ; 
    NBIN[0]  = H->GetXaxis()->GetNbins();
    XBIN[0]  = (double)H->GetXaxis()->GetBinWidth(1);
    XMIN[0]  = (double)H->GetXaxis()->GetBinLowEdge(1);
    
    NBIN[1]  = H->GetYaxis()->GetNbins();
    XBIN[1]  = (double)H->GetYaxis()->GetBinWidth(1);
    XMIN[1]  = (double)H->GetYaxis()->GetBinLowEdge(1);

    sprintf(CTIT, "%s", H->GetTitle() );
  }
 
  if ( NDIM >= 1 ) {
    XMAX[0]  = XMIN[0] + ( XBIN[0] * (double)NBIN[0] ) ;
    if ( LDMP ) { 
      printf(" xxx X-axis: NBIN=%d  RANGE=%.3f to %.3f  BIN=%.3f \n", 
	     NBIN[0], XMIN[0], XMAX[0], XBIN[0] ); 
      printf(" xxx CTIT = '%s' \n", CTIT);
    }
  }

  if ( NDIM == 2 ) {
    XMAX[1]  = XMIN[1] + ( XBIN[1] * (double)NBIN[1] ) ;
    if ( LDMP) { 
      printf(" xxx Y-axis: NBIN=%d  RANGE=%.3f to %.3f  BIN=%.3f \n", 
	     NBIN[1], XMIN[1], XMAX[1], XBIN[1] ); 
    }
  }



  /*
  sprintf(MSGERR1,"Cannot read root histo-bins for ID=%d", ID);
  sprintf(MSGERR2,"  " ); 
  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  */

}  // end of SNHIST_RDBINS_ROOT


// ========================================================
void SNHIST_RDCONT_ROOT(int NDIM, int ID, int NRDBIN, double *CONTENTS) {

  int ibin, NBX, NBY, ix, iy ;
  char HNAM[40];
  //  char fnam[] = "SNHIST_RDCONT_ROOT" ;

  // --------------- BEIGN ------------------

  sprintf(HNAM,"H%5.5d", ID);

  if ( NDIM == 1 ) { 
    TH1F *H = (TH1F*)TFILE_ROOT_READ->Get(HNAM) ; 
    for(ibin=1; ibin <= NRDBIN; ibin++ ) 
      { CONTENTS[ibin-1] = (double)H->GetBinContent(ibin);  }
  }


  if ( NDIM == 2 ) { 
    TH2F *H = (TH2F*)TFILE_ROOT_READ->Get(HNAM) ; 

    ibin = 0;
    NBX  = H->GetXaxis()->GetNbins();
    NBY  = H->GetYaxis()->GetNbins();

    // fill CONTENTS the same way that HBOOK does for 2D.
    for (iy=1; iy <= NBY; iy++ ) {
      for(ix=1; ix <= NBX; ix++ ) {
	CONTENTS[ibin] = (double)H->GetBinContent(ix,iy); 
	ibin++ ;
      }
    }
  }



} // end of SNHIST_RDCONT_ROOT

// ==========================================
//
//          DIRECTORY FUNCTIONS
//
// ==========================================

void MAKEDIR_ROOT(char *CCID, int CID) {
  // Created Mar 3, 2013 by R.Kessler
  // if CID >= 0 then append 'SN' prefix
  // if CID <  0 then CCID is an arbitrary  dirName

  char fnam[] = "MAKEDIR_ROOT" ;
  char dirName[40];
  int IFIT, NFIT ;


  if ( CID >= 0 ) 
    {  sprintf(dirName, "SN%s", CCID);   }
  else
    {  sprintf(dirName, "%s", CCID); }


  //  add _FIT[IFIT] suffix to subdir if there are multiple fits per SN
  if ( CID >= 0 ) {
    NFIT   = SNLCPAK_OUTPUT.NFIT_PER_SN  ;
    IFIT   = SNLCPAK_OUTPUT.NLCPAK  ; // starts at 0 since it sums later
    if ( NFIT > 1 ) { sprintf(dirName,"%s_FIT%d", dirName, IFIT); }
  }



  if ( CID >= 0 ) {
    // create subdir under /MONSN dir.
    TD_MONSN->cd(); 
    TD_MONSN_SUBDIR = TD_MONSN->mkdir(dirName);  
    TD_MONSN_SUBDIR->cd();
    printf("\t %s : Created %s/%s\n", fnam, MONSN_DIRNAME, dirName);
  }
  else {
    // create subdir under top directory
    TD_TEMP = TFILE_ROOT_NEW->mkdir(dirName) ;
    TD_TEMP->cd();
    printf("\t %s : Created %s\n", fnam, dirName);
  }

  fflush(stdout);

} // end of  MAKEDIR_ROOT


void CDTOPDIR_ROOT(void) {
  char fnam[] = "CDTOPDIR_ROOT" ;
  // Created Mar 3, 2013 by R.Kessler
  printf("\t %s : return to top dir.\n", fnam ); fflush(stdout);
  TFILE_ROOT_NEW->cd();  
}

// ==========================================
//
//          SNLCPAK FUNCTIONS
//
// ==========================================


void SNLCPAK_SURVEY_ROOT(void) {

  // one-time init
  // * store global information 
  // * init Tree(s) for plotting

  int ipak, NLCPAK ;

  // allocate TREEDATA memory
  SNLCPAK_MALLOC_TREEDATA();


  // set NLCPAK = number of plots to make per SN.
  // Note that NFIT_PER_SN = NLCPAK if fits are done.
  // But to sotre/plot raw data with no fit, 
  // NFIT_PER_SN=0 and NLCPAK=1.
  NLCPAK = SNLCPAK_OUTPUT.NFIT_PER_SN ;
  if ( NLCPAK == 0 ) { NLCPAK = 1; }


  // create and fill global tree with global info
  // (mostly char strings)
  TTree *Tree_global ;
  Tree_global = new TTree("GLOBAL", "GLOBAL INFO") ;

  Tree_global->Branch( "SURVEY", SNLCPAK_OUTPUT.SURVEY, 
		       "SURVEY/C"  ) ;

  Tree_global->Branch( "SURVEY_FILTERS", SNLCPAK_OUTPUT.SURVEY_FILTERS,
		       "SURVEY_FILTES/C"  ) ;

  Tree_global->Branch( "FILTERSTRING", FILTERSTRING, 
		       "FILTERSTRING/C"  ) ;

  Tree_global->Branch( "VERSION_PHOTOMETRY", 
		       SNLCPAK_OUTPUT.VERSION_PHOTOMETRY,
		       "VERSION_PHOTOMETRY/C"  ) ;

  Tree_global->Branch( "VERSION_SNANA", 
		       SNLCPAK_OUTPUT.VERSION_SNANA,
		       "VERSION_SNANA/C"  ) ;

  Tree_global->Branch( "NLCPAK", &NLCPAK, "NLCPAK/I" );

  Tree_global->Branch( "SNLCPAK_TREENAMES", 
		       SNLCPAK_TREEDATA.TREENAMES,
		       "SNLCPAK_TREENAMES/C"  ) ;

  // -------------------------------
  // create tree for each set of plots/fit

  SNLCPAK_TREEDATA.TREENAMES[0] = 0 ;
  for ( ipak=0; ipak < NLCPAK; ipak++ ) 
    { SNLCPAK_makeTree_ROOT(ipak) ; }


  // fill the global tree after the makeTree calls so that
  // we have the list of TREENAMES
  Tree_global->Fill();


  // create MONSN subdir to contain subdir for each SN with plots
  TD_MONSN = TFILE_ROOT_NEW->mkdir(MONSN_DIRNAME);

} // end of SNLCPAK_SURVEY_ROOT


// ------------------------------------------------
void SNLCPAK_makeTree_ROOT(int ipak) {

  // Create tree to store light curve data and best-fit.
  // For each fit (ipak), create a separate tree.

  int NFIT ;
  char  TREE_NAME[40], TREE_TITLE[100] ;
  //  char  fnam[] = "SNLCPAK_makeTree_ROOT"  ;

  // --------------- BEGIN ---------------

  NFIT = SNLCPAK_OUTPUT.NFIT_PER_SN ; 


  if ( NFIT > 0  ) 
    {  sprintf(TREE_TITLE,"light curves and best-fit function"); }
  else
    {  sprintf(TREE_TITLE,"light curves" ); }


  // get name of tree based on number of expected fits per SN
  if ( NFIT <= 1 ) 
    { sprintf(TREE_NAME,"SNLCPAK"); }
  else
    { sprintf(TREE_NAME,"SNLCPAK_FIT%d", ipak ); }

  // increment space-separated list of tree names
  strcat(SNLCPAK_TREEDATA.TREENAMES,TREE_NAME);
  strcat(SNLCPAK_TREEDATA.TREENAMES,"  ");

  printf("\t Create LightCurve Tree %s \n", TREE_NAME ) ;
  fflush(stdout);

  // define new TTree
  SNLCPAK_TREE[ipak] = new TTree(TREE_NAME, TREE_TITLE) ;

  SNLCPAK_TREE[ipak]->Branch( "CCID", SNLCPAK_OUTPUT.CCID, "CCID/C"  ) ;

  // init variable-length arrays for Tree-Branches
  SNLCPAK_Branch_DISPLAYTEXT(ipak);  // display strings
  SNLCPAK_Branch_lightCurve(ipak);  // light curve info
  SNLCPAK_Branch_filter(ipak);      // filter info: MJD, NDOF, CHI2




} // end of SNLCPAK_makeTree_ROOT


// ------------------------------------------------
void SNLCPAK_FILL_TREEDATA(void) {

  // fill memory SNLCPAK_TREEDATA struct for root.
  // There are NO root calls here.
  //
  // Apr 17, 2013: fix aweful index bug setting ep2 for FITFUN
  // May 20, 2016: fill SIMFLUX_F

  int FLAG, ep, ep2, NTOT, NDATA, NFITFUN, IFILTOBS, NSIM ;
  double dval;
  char *cfiltobs;
  char fnam[] = "SNLCPAK_FILL_TREEDATA" ;

  // ----------- BEGIN ----------


  NDATA   = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FLUXDATA] ;  
  NSIM    = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FLUXSIM] ;  
  NFITFUN = SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FITFUN] ;
  NTOT    = NDATA + NFITFUN ;

  SNLCPAK_TREEDATA.NDATA    = NDATA ;
  SNLCPAK_TREEDATA.NFITFUN  = NFITFUN ;
  SNLCPAK_TREEDATA.NTOT     = NTOT ;
  SNLCPAK_TREEDATA.NFILT    = SNLCPAK_OUTPUT.NFILTDEF_SURVEY ; 


  if ( NTOT >= SNLCPAK_TREEDATA.MXEP ) {
    sprintf(MSGERR1,"NEPTOT=%d(DATA) + %d(FITFUN) = %d",
	    NDATA, NFITFUN, NTOT);
    sprintf(MSGERR2,"exceeds memory alloc. of MXEP=%d",
	    SNLCPAK_TREEDATA.MXEP ) ;
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2); 
  }

  // now fill relevant structure

  FLAG = SNLCPAK_EPFLAG_FLUXDATA ;
  for(ep=0; ep < NDATA; ep++ ) {

    dval          = SNLCPAK_OUTPUT.TOBS[FLAG][ep] ;
    SNLCPAK_TREEDATA.TOBS_F[ep]    = (float)dval ;

    dval          = SNLCPAK_OUTPUT.MJD[FLAG][ep] ;
    SNLCPAK_TREEDATA.MJD_D[ep]    = dval ;

    IFILTOBS      = SNLCPAK_OUTPUT.IFILTOBS[FLAG][ep] ;
    SNLCPAK_TREEDATA.IFILTOBS[ep]  = IFILTOBS ;

    cfiltobs = SNLCPAK_TREEDATA.CFILTOBS[ep] ;
    sprintf(cfiltobs,"%c", FILTERSTRING[IFILTOBS]) ; 

    dval          = SNLCPAK_OUTPUT.EPDATA[FLAG][ep] ;
    SNLCPAK_TREEDATA.FLUX_F[ep]    = (float)dval ;

    dval          = SNLCPAK_OUTPUT.EPDATA_ERR[FLAG][ep] ;
    SNLCPAK_TREEDATA.FLUXERR_F[ep] = (float)dval ;

    if ( NSIM > 0 ) 
      { dval = SNLCPAK_OUTPUT.EPDATA[SNLCPAK_EPFLAG_FLUXSIM][ep] ; }
    else
      { dval = 0.0 ; }	
    SNLCPAK_TREEDATA.SIMFLUX_F[ep]    = (float)dval ;    

    dval          = SNLCPAK_OUTPUT.EPDATA[SNLCPAK_EPFLAG_REJECT][ep] ;
    SNLCPAK_TREEDATA.IREJECT[ep]   = (int  )dval ;

    if ( NFITFUN > 0 ) 
      { dval = SNLCPAK_OUTPUT.EPDATA[SNLCPAK_EPFLAG_CHI2][ep] ; }
    else
      { dval = 0.0 ; }
    SNLCPAK_TREEDATA.CHISQ_F[ep]   = (float)dval ;

    SNLCPAK_TREEDATA.IFLAGDATA[ep] = 1 ;

  }


  // include the optional best-fit function

  FLAG = SNLCPAK_EPFLAG_FITFUN ;
  for(ep2=0; ep2 < NFITFUN; ep2++ ) {

    //    ep = ep2 + NDATA-1;  // aweful bug .....
    ep = ep2 + NDATA ;         // fixed Apr 17 2013

    dval          = SNLCPAK_OUTPUT.TOBS[FLAG][ep2] ;
    SNLCPAK_TREEDATA.TOBS_F[ep]    = (float)dval ;

    dval          = SNLCPAK_OUTPUT.MJD[FLAG][ep2] ;
    SNLCPAK_TREEDATA.MJD_D[ep]    = dval ;

    IFILTOBS      = SNLCPAK_OUTPUT.IFILTOBS[FLAG][ep2] ;
    SNLCPAK_TREEDATA.IFILTOBS[ep]  = IFILTOBS ;

    cfiltobs = SNLCPAK_TREEDATA.CFILTOBS[ep2] ;
    sprintf(cfiltobs,"%c", FILTERSTRING[IFILTOBS]) ; 
    
    dval          = SNLCPAK_OUTPUT.EPDATA[FLAG][ep2] ;
    SNLCPAK_TREEDATA.FLUX_F[ep]    = (float)dval ;

    dval          = SNLCPAK_OUTPUT.EPDATA_ERR[FLAG][ep2] ;
    SNLCPAK_TREEDATA.FLUXERR_F[ep] = (float)dval ;

    SNLCPAK_TREEDATA.IREJECT[ep]   = 0   ;
    SNLCPAK_TREEDATA.CHISQ_F[ep]   = 0.0 ;
    SNLCPAK_TREEDATA.IFLAGDATA[ep] = 0   ;
  }


  // -------

  // load filter part of TREEDATA structure

  int IFILT, indx ;
  for ( IFILT=0; IFILT < SNLCPAK_TREEDATA.NFILT; IFILT++ ) {

    FLAG = SNLCPAK_BANDFLAG_PKMJD ; indx=FLAG-100 ;
    SNLCPAK_TREEDATA.BAND_PKMJD[IFILT] = 
      SNLCPAK_OUTPUT.FILTDATA[indx][IFILT];
    SNLCPAK_TREEDATA.BAND_PKMJDERR[IFILT] = 
      SNLCPAK_OUTPUT.FILTDATA_ERR[indx][IFILT];

    FLAG = SNLCPAK_BANDFLAG_PKFLUX ; indx=FLAG-100 ;
    SNLCPAK_TREEDATA.BAND_PKFLUX_F[IFILT] = 
      (float)SNLCPAK_OUTPUT.FILTDATA[indx][IFILT] ;

    FLAG = SNLCPAK_BANDFLAG_NDOF   ; indx=FLAG-100 ;
    SNLCPAK_TREEDATA.BAND_NDOF[IFILT] = 
      (int)SNLCPAK_OUTPUT.FILTDATA[indx][IFILT];

    FLAG = SNLCPAK_BANDFLAG_CHI2   ; indx=FLAG-100 ;
    SNLCPAK_TREEDATA.BAND_CHI2_F[IFILT] = 
      (float)SNLCPAK_OUTPUT.FILTDATA[indx][IFILT];
  }



} // end of SNLCPAK_FILL_TREEDATA

// ------------------------------------------------
void SNLCPAK_FILL_ROOT(void) {

  // Store light curve data with (optional) best-fit function.
  // Include REJECT flags and chi2 per obs.
  
  int NTOT, NFILT, ipak ;
  //  char fnam[] = "SNLCPAK_FILL_ROOT"  ;

  // -------------- BEGIN ---------------

  // load array sizes.
  NTOT   = 
    SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FLUXDATA]  +
    SNLCPAK_OUTPUT.NOBS[SNLCPAK_EPFLAG_FITFUN] ;
  NFILT    = SNLCPAK_OUTPUT.NFILTDEF_SURVEY ; 


  // fill TREEDATA 
  SNLCPAK_FILL_TREEDATA();

  // fill the table with root call
  ipak = SNLCPAK_OUTPUT.NLCPAK - 1 ;

  if ( VBOSE_FLAG_ROOT ) 
    { printf("\t    Fill Tree for ipak = %d \n", ipak); fflush(stdout);  }

  SNLCPAK_TREE[ipak]->Fill();   

} // end of SNLCPAK_FILL_ROOT



// --------------------------------------
void  SNLCPAK_MALLOC_TREEDATA(void) {

  int MEM_F, MEM_I, MEM_D, MEM_C ;
  int NFILT, NEP, i ;
  // ------------- BEGIN ---------------

  NFILT  = SNLCPAK_OUTPUT.NFILTDEF_SURVEY ; 
  NEP    = NFILT * MXEP_PER_FILT ;      // data + FITFUN

  SNLCPAK_TREEDATA.MXEP = NEP ; // store to check for mem leaks

  if ( VBOSE_FLAG_ROOT ) {
    printf("\t Allocate SNLCPAK.TREEDATA:  NEP=%d, NFILT=%d\n",
	   NEP, NFILT ) ;
    fflush(stdout);
  }

  // -------
  MEM_F = NEP * sizeof(float);
  MEM_D = NEP * sizeof(double);
  MEM_I = NEP * sizeof(int);  
  MEM_C = NEP * sizeof(char*); 

  SNLCPAK_TREEDATA.TOBS_F    = (float*)malloc( MEM_F ) ;
  SNLCPAK_TREEDATA.MJD_D     = (double*)malloc( MEM_D ) ;
  SNLCPAK_TREEDATA.IFILTOBS  = (int  *)malloc( MEM_I ) ;
  SNLCPAK_TREEDATA.FLUX_F    = (float*)malloc( MEM_F ) ;
  SNLCPAK_TREEDATA.FLUXERR_F = (float*)malloc( MEM_F ) ;
  SNLCPAK_TREEDATA.SIMFLUX_F = (float*)malloc( MEM_F ) ;
  SNLCPAK_TREEDATA.IREJECT   = (int  *)malloc( MEM_I ) ;
  SNLCPAK_TREEDATA.CHISQ_F   = (float*)malloc( MEM_F ) ;
  SNLCPAK_TREEDATA.IFLAGDATA = (int  *)malloc( MEM_I ) ;

  SNLCPAK_TREEDATA.CFILTOBS   = (char**)malloc( MEM_C ) ;
  for(i=0; i < NEP; i++ ) 
    { SNLCPAK_TREEDATA.CFILTOBS[i]  = (char*)malloc( 2*sizeof(char) ) ; }

  // -------
  MEM_D = NFILT * sizeof(double);
  MEM_F = NFILT * sizeof(float);
  MEM_I = NFILT * sizeof(int);  
  SNLCPAK_TREEDATA.BAND_PKMJD     = (double * )malloc( MEM_D );
  SNLCPAK_TREEDATA.BAND_PKMJDERR  = (double * )malloc( MEM_D );
  SNLCPAK_TREEDATA.BAND_PKFLUX_F  = (float  * )malloc( MEM_F );
  SNLCPAK_TREEDATA.BAND_CHI2_F    = (float  * )malloc( MEM_F );
  SNLCPAK_TREEDATA.BAND_NDOF      = (int    * )malloc( MEM_I );
    

}  // end of SNLCPAK_MALLOC_TREEDATA



// --------------------------------------
void  SNLCPAK_Branch_DISPLAYTEXT(int ipak) {

  char *ptrTEXT ;
  int  MXTEXT, i ;
  //  char fnam[] = "SNLCPAK_Branch_DISPLAYTEXT" ;

  // ----------- BEGIN ---------

  MXTEXT = SNLCPAK_OUTPUT.MXTEXT ;


  if ( VBOSE_FLAG_ROOT ) {
    printf("\t Init Tree Branch for DISPLAYTEXT (MXTEXT=%d)\n", MXTEXT);
    fflush(stdout);
  }

  // set pointers to display text


  // define branch for all MXTEXT strings
  
  SNLCPAK_TREE[ipak]->Branch( "NTEXT",  &SNLCPAK_OUTPUT.NTEXT,      
			     "NTEXT/I"  ) ;
  
  
  for ( i=0; i<MXTEXT; i++ ) 
    { SNLCPAK_TREEDATA.DISPLAYTEXT[i] = SNLCPAK_OUTPUT.DISPLAYTEXT[i] ; }

  /*  
  SNLCPAK_TREE[ipak]->Branch ( "DISPLAYTEXT", &SNLCPAK_TREEDATA.DISPLAYTEXT,
			      "DISPLAYTEXT[NTEXT]/C" ) ;
  */  


    
  ptrTEXT = SNLCPAK_OUTPUT.DISPLAYTEXT[0] ;
  SNLCPAK_TREE[ipak]->Branch ( "DISPLAYTEXT0", ptrTEXT, "DISPLAYTEXT0/C") ;

  ptrTEXT = SNLCPAK_OUTPUT.DISPLAYTEXT[1] ;
  SNLCPAK_TREE[ipak]->Branch ( "DISPLAYTEXT1", ptrTEXT, "DISPLAYTEXT1/C" ) ;

  ptrTEXT = SNLCPAK_OUTPUT.DISPLAYTEXT[2] ;
  SNLCPAK_TREE[ipak]->Branch ( "DISPLAYTEXT2", ptrTEXT, "DISPLAYTEXT2/C" ) ;

  ptrTEXT = SNLCPAK_OUTPUT.DISPLAYTEXT[3] ;
  SNLCPAK_TREE[ipak]->Branch ( "DISPLAYTEXT3", ptrTEXT, "DISPLAYTEXT3/C" ) ;
  

  /*  xxx this loop does not even compile ??? xxxx
  for ( i=0; i<MXTEXT; i++ ) {
    ptrTEXT = SNLCPAK_OUTPUT.DISPLAYTEXT[i] ;
    sprintf(VARNAME, "DISPLAYTEXT%d", i );
    sprintf(LEAF,    "%s/C", VARNAME ) ;
    SNLCPAK_TREE[ipak]->Branch ( VARNAME, ptrTEXT, *LEAF ) ;
  }
  xxx  */


} // end of SNLCPAK_Branch_DISPLAYTEXT

// --------------------------------------
void  SNLCPAK_Branch_lightCurve(int ipak) {

  // init light curve info for SNLCPAK_TREEDATA,
  // and fill SNLCPAK_TREEDATA structure needed for ->Fill.

  //  char fnam[] = "SNLCPAK_Branch_lightCurve" ;

  // ----------- BEGIN ---------

  if ( VBOSE_FLAG_ROOT ) {
    printf("\t Init Branch for Lightcurve\n" );
    fflush(stdout);
  }


  // Define tree columns:
  // NDATA Tobs  cfiltobs ifiltobs flux fluxErr reject chi2  DATAFLAG

  // start with scalars with size of arrays.
  SNLCPAK_TREE[ipak]->Branch( "NDATA",      &SNLCPAK_TREEDATA.NDATA,      
			      "NDATA/I"   ) ;

  SNLCPAK_TREE[ipak]->Branch( "NFITFUN",    &SNLCPAK_TREEDATA.NFITFUN, 
			      "NFITFUN/I" ) ;

  SNLCPAK_TREE[ipak]->Branch( "NTOT",       &SNLCPAK_TREEDATA.NTOT,  
			      "NTOT/I"    ) ;  

  // now book the arrays.
  SNLCPAK_TREE[ipak]->Branch( "TOBS",        SNLCPAK_TREEDATA.TOBS_F,     
			      "TOBS[NTOT]/F" );  

  SNLCPAK_TREE[ipak]->Branch( "MJD",        SNLCPAK_TREEDATA.MJD_D,     
			      "MJD[NTOT]/D" );  

  SNLCPAK_TREE[ipak]->Branch( "IFILTOBS",    SNLCPAK_TREEDATA.IFILTOBS, 
			      "IFILTOBS[NTOT]/I" );  

  SNLCPAK_TREE[ipak]->Branch( "CFILTOBS",   SNLCPAK_TREEDATA.CFILTOBS[0], 
			      "CFILTOBS[NTOT]/C" );  

  SNLCPAK_TREE[ipak]->Branch( "FLUXCAL",     SNLCPAK_TREEDATA.FLUX_F, 
			      "FLUXCAL[NTOT]/F" );

  SNLCPAK_TREE[ipak]->Branch( "FLUXCAL_ERR", SNLCPAK_TREEDATA.FLUXERR_F, 
			      "FLUXCAL_ERR[NTOT]/F" );

  SNLCPAK_TREE[ipak]->Branch( "FLUXCAL_SIM",     SNLCPAK_TREEDATA.SIMFLUX_F, 
			      "FLUXCAL_SIM[NTOT]/F" );

  SNLCPAK_TREE[ipak]->Branch( "REJECT",      SNLCPAK_TREEDATA.IREJECT,  
			      "REJECT[NTOT]/I" ); 

  SNLCPAK_TREE[ipak]->Branch( "FITCHI2",     SNLCPAK_TREEDATA.CHISQ_F, 
			      "FITCHI2[NTOT]/F" );

  SNLCPAK_TREE[ipak]->Branch( "IFLAGDATA",   SNLCPAK_TREEDATA.IFLAGDATA,  
			      "IFLAGDATA[NTOT]/I" ); 


} // SNLCPAK_Branch_lightCurve

// --------------------------------------
void  SNLCPAK_Branch_filter(int ipak) {

    // filter info: MJD, NDOF, CHI2 
  //  char fnam[] = "SNLCPAK_Branch_filter" ;

  // ---------- BEGIN ----------

  if ( VBOSE_FLAG_ROOT ) {
    printf("\t Init Branch for Filter\n" );
    fflush(stdout);
  }

  SNLCPAK_TREE[ipak]->Branch( "NFILTDEF_SURVEY",  &SNLCPAK_TREEDATA.NFILT,  
			 "NFILTDEF_SURVEY/I"    ) ;  

  SNLCPAK_TREE[ipak]->Branch( "BAND_PEAKMJD",  SNLCPAK_TREEDATA.BAND_PKMJD,
			 "BAND_PEAKMJD[NFILTDEF_SURVEY]/D" );  

  SNLCPAK_TREE[ipak]->Branch("BAND_PEAKMJD_ERR",SNLCPAK_TREEDATA.BAND_PKMJDERR,
			 "BAND_PEAKMJD_ERR[NFILTDEF_SURVEY]/D" );  
  
  SNLCPAK_TREE[ipak]->Branch( "BAND_PEAKFLUX",  SNLCPAK_TREEDATA.BAND_PKFLUX_F,
			 "BAND_PEAKFLUX[NFILTDEF_SURVEY]/F" );  

  SNLCPAK_TREE[ipak]->Branch( "BAND_NDOF",  SNLCPAK_TREEDATA.BAND_NDOF,
			 "BAND_NDOF[NFILTDEF_SURVEY]/I" );  

  SNLCPAK_TREE[ipak]->Branch( "BAND_CHI2",  SNLCPAK_TREEDATA.BAND_CHI2_F,
			 "BAND_CHI2[NFILTDEF_SURVEY]/F" );  


}  // end of SNLCPAK_Branch_filter


// --------------------------------------
void modify_text_4root(char *Text) {

  char TEXT_ORIG[200]; 
  //  char fnam[] = "modify_text_4root" ;

  // ------------ BEGIN --------------

  // below are dummies; need to use TLatex.

  sprintf(TEXT_ORIG,"%s", Text);
  sprintf(Text, "%s", replace_str(TEXT_ORIG, "MU", "Mu") ) ;

  sprintf(TEXT_ORIG,"%s", Text);
  sprintf(Text, "%s", replace_str(TEXT_ORIG, "+-", "+_") ) ;

} // end of modify_text_4root





// ========================================================
void MERGE_ROOT(int NFILE, char **INFILES, char *OUTFILE) {

  // Extract code pieces from
  // http://root.cern.ch/svn/root/trunk/main/src/hadd.cxx


  int ifile ;
  char fnam[] = "MERGE_ROOT" ;

  // ------------- BEGIN --------------

  printf("\n"); fflush(stdout);

  gSystem->Load("libTreePlayer");

  TFileMerger merger(kFALSE,kFALSE);
  //  merger.SetMsgPrefix(fnam);
  merger.SetPrintLevel(1);
 
  for(ifile=0; ifile < NFILE; ifile++ ) 
    { merger.AddFile( INFILES[ifile] ) ; }


  merger.OutputFile(OUTFILE,kFALSE,1) ;

  Bool_t status = merger.Merge();

  printf("\t %s status = %d \n", fnam, status);
  fflush(stdout);

} // end of MERGE_ROOT

