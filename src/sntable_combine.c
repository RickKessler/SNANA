/* ====================================

  Created Jan 2017

  Combine multiple tables with same set of CIDs, but different 
  variables. Default output table has same format as first input 
  table. User can specify -outFile <outFile> with any format as
  specified by the suffix (HBOOK, ROOT, TEXT, FITRES). Default 
  tableName to combine is "FITRES", but any tableName can be 
  specified with -tableName argument.

  If a variable name is repeated after the first input file,
  _[NFILE] is appended. For example, if 2nd file has duplicate 
  variable names "c x1" (e.g., re-running SALT2 with different 
  options), then the duplicate varNames are c_2 and x1_2 in the 
  output "combined" tabe file.

  This program includes the functionality of "combine_fitres.exe".
  This older program works only on input TEXT files, while this
  new program (sntable_combine) works on input files in any valid
  SNANA format: HBOOK, ROOT, TEXT.
  

  USAGE:
  sntable_combine.exe -inFile FITOPT000.ROOT FITOPT000_EXTRAS.FITRES
     [create default outFile COMBINE_FITRES.ROOT]

  sntable_combine.exe -inFile F0.ROOT F1.FITRES F2.FITRES F3.FITRES
     [create default outFile COMBINE_FITRES.ROOT]
    
  sntable_combine.exe -inFile FITOPT000.FITRES FITOPT000_EXTRAS.FITRES
     [create default outFile COMBINE_FITRES.TEXT ]

  sntable_combine.exe -inFile FITOPT000.FITRES FITOPT000_EXTRAS.FITRES \
                      -outFile COMBINE_FITOPT000.ROOT
     [outFile format is ROOT]

  sntable_combine.exe -inFile FITOPT000.ROOT FITOPT000_EXTRAS.FITRES \
                      -tableName SNANA
     [create default outFile COMBINE_SNANA.ROOT ]


    HISTORY
  ~~~~~~~~~~~~~

 Jun 16 2019:  call  TABLEFILE_CLOSE(f), EXCEPT for TEXT format

===================================== */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

/* xxx
// Mar 22 2017: define larger HBOOK memory for this program 
#define USER_HBOOK_MEM  // flag to NOT define HBOOK_MEM in sntools_output
#define HBOOK_MEM  1000000
struct pawC  { float hmem[HBOOK_MEM];  } pawc_ ;
xxx */

#include "sntools.h"
#include "sntools_output.h" 

// -------- global variables ------------

#define MXFILE_COMBINE MXFILE_AUTOSTORE 

int NROW_STDOUT_UPDATE ;
time_t t_start, t_read, t_end ;

struct {
  int  NFILE_COMBINE ;
  char INFILE_LIST[MXFILE_COMBINE][MXCHAR_FILENAME];
  char OUTFILE[MXCHAR_FILENAME];

  char TABLENAME[40];
  int  TABLEID;

  // computed from inputs
  int IFILETYPE[MXFILE_COMBINE] ;
} INPUTS ;


struct {
  int NVAR_TOT; // NVAR summed over all input files 
  int NROW;     // NROW to write  = NROW in first file
  int IFILETYPE ;
  char   CCID[MXCHAR_CCID];
  int    *VAL_I;  // store per ivar
  float  *VAL_F ;
  double *VAL_D;
  char   **VAL_C ;
} OUTPUT ;

// ================================
// Function Declarations
// ================================

void  PARSE_ARGV(int argc, char **argv);
void  testRead(void);
void  malloc_OUTPUT_STORAGE(void);
void  openFile_combine(void);
void  sntable_combine_init(void);
void  sntable_combine_fill(int irow);
void  sntable_combine_summary(void);

// =========================================
int main(int argc, char **argv) {

  int iFile, optMask, irow, IFILETYPE  ;
  char *f ;
  // ------------ BEGIN ------------

  t_start = time(NULL);

  PARSE_ARGV(argc,argv);

  optMask = 5; // 1=verbose; 4=append each file
  OUTPUT.NVAR_TOT = 0 ;
  for(iFile=0; iFile < INPUTS.NFILE_COMBINE; iFile++ ) {
    f = INPUTS.INFILE_LIST[iFile] ;
    printf("\n@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@\n");
    fflush(stdout);
    SNTABLE_AUTOSTORE_INIT(f, INPUTS.TABLENAME, "ALL", optMask);
    IFILETYPE = SNTABLE_AUTOSTORE[iFile].IFILETYPE ;
    INPUTS.IFILETYPE[iFile] = IFILETYPE ;
    if(IFILETYPE != IFILETYPE_TEXT) { TABLEFILE_CLOSE(f); } // June 16 2019
    OUTPUT.NVAR_TOT += SNTABLE_AUTOSTORE[iFile].NVAR ;
  }

  OUTPUT.NROW = SNTABLE_AUTOSTORE[0].NROW ;
  printf("\n Prepare output table with %d variables & %d rows. \n", 
	 OUTPUT.NVAR_TOT, OUTPUT.NROW );
  fflush(stdout);
  

  // set frequency of screen update.
  NROW_STDOUT_UPDATE = 200; 
  if ( OUTPUT.NROW > 5000  ) { NROW_STDOUT_UPDATE = 1000; }
  if ( OUTPUT.NROW > 50000 ) { NROW_STDOUT_UPDATE = 5000; }

  t_read = time(NULL);

  // malloc globals to use for output storage
  malloc_OUTPUT_STORAGE();

  // -------------------------------------
  // init output table
  openFile_combine();

  sntable_combine_init();

  printf("\n Begin filling combined %s table . . . .\n", INPUTS.TABLENAME); 
  fflush(stdout);

  for(irow=0; irow < OUTPUT.NROW; irow++ ) 
    { sntable_combine_fill(irow); }

  TABLEFILE_CLOSE(INPUTS.OUTFILE);

  // check total process time
  sntable_combine_summary();

  //  testRead();

  return(0);

} // end main

// ===============================
void  PARSE_ARGV(int argc, char **argv) {

  int  NARGV, i, NF, INFILE_FLAG ;

  // ------------- BEGIN ---------------

  NARGV = argc ;
  INPUTS.NFILE_COMBINE = NF = 0 ;
  INFILE_FLAG=0;
  INPUTS.OUTFILE[0] = 0;
  sprintf(INPUTS.TABLENAME, "%s", TABLENAME_FITRES);

  for ( i = 0; i < NARGV ; i++ ) {
        
    if ( strcmp(argv[i],"-outFile" ) == 0 ) {
      i++ ; sprintf(INPUTS.OUTFILE,"%s", argv[i]);
      INFILE_FLAG=0;
      continue ;
    }

    if ( strcmp(argv[i],"-tableName" ) == 0 ) {
      i++ ; sprintf(INPUTS.TABLENAME,"%s", argv[i]);
      INFILE_FLAG=0;
      continue ;
    }

    if ( strcmp(argv[i],"-inFile" ) == 0 ) 
      { INFILE_FLAG=1; continue; }

    if ( INFILE_FLAG ) {
      sprintf(INPUTS.INFILE_LIST[NF], "%s", argv[i] );
      NF++ ;
      continue ;
    }
  }

  INPUTS.NFILE_COMBINE = NF ;

  printf("\n Inputs TABLE files to combine: \n");
  for(i=0 ; i < NF; i++ ) 
    { printf("\t %s \n", INPUTS.INFILE_LIST[i] ); }

  printf("\n Output file with combined table: %s \n", INPUTS.OUTFILE );
  fflush(stdout);

  return ;

} // end PARSE_ARGV


// ==============================
void  malloc_OUTPUT_STORAGE(void) {

  int ivar;
  int MEMD  = sizeof(double);
  int MEMF  = sizeof(float);
  int MEMI  = sizeof(int);
  int MEMC1 = sizeof(char*);
  int MEMC0 = sizeof(char);

  int NVAR = OUTPUT.NVAR_TOT ;
  // ------------- BEGIN -----------

  OUTPUT.VAL_D = (double*)malloc( MEMD  * NVAR );
  OUTPUT.VAL_F = (float *)malloc( MEMF  * NVAR );
  OUTPUT.VAL_I = (int   *)malloc( MEMI  * NVAR );
  OUTPUT.VAL_C = (char **)malloc( MEMC1 * NVAR );

  for(ivar=0; ivar < OUTPUT.NVAR_TOT; ivar++ ) {
    OUTPUT.VAL_C[ivar] = (char*)malloc( MEMC0 * MXCHAR_CCID );
  }

  return ;

} // end malloc_OUTPUT_STORE


// ============================================
void openFile_combine(void) {

  // if user does NOT give outFile argument, then default
  // output filename is 
  //       COMBINE_[tableName].[SUFFIX]
  // where SUFFIX = ROOT or HBOOK or TEXT or FITRES
  //

  int  IFILETYPE ;
  char openOpt[20], SUFFIX[20] ;
  //  char fnam[] = "openFile_combine" ;

  // ----------- BEGIN -------------

  if ( strlen(INPUTS.OUTFILE) == 0 ) {
    IFILETYPE = INPUTS.IFILETYPE[0] ; // use fileType of 1st file
    sprintf(SUFFIX, "%s", STRING_TABLEFILE_TYPE[IFILETYPE]); 
    sprintf(INPUTS.OUTFILE, "COMBINE_%s.%s", INPUTS.TABLENAME, SUFFIX);
  }

  sprintf(openOpt,"new");
  OUTPUT.IFILETYPE = TABLEFILE_OPEN(INPUTS.OUTFILE,openOpt);
  CDTOPDIR_OUTPUT(1);


  return ;
} // end openFile_combine

// ============================================
void  sntable_combine_init(void) {

  // Create FITRES table and init each variable.

  int USE4TEXT=1;
  int TABLEID ;
  int ifile, ivar, NVAR_TOT, ICAST, NFILE, NVAR ;
  char CCAST[12], VARNAME[60], FORMAT[20] ;
  void *PTRVAR = NULL;

  char BLOCK[] = "COMBINE" ;
  //  char fnam[]  = "sntable_combine_init" ;

  // ----------- BEGIN -----------
  
  FORMAT[0] = 0 ;
  if ( OUTPUT.IFILETYPE == IFILETYPE_TEXT ) { sprintf(FORMAT,"key") ; }

  if ( strcmp(INPUTS.TABLENAME,TABLENAME_FITRES) == 0 ) 
    { INPUTS.TABLEID = TABLEID_FITRES ;  }
  else
    { INPUTS.TABLEID = TABLEID_SNANA ;  }

  SNTABLE_CREATE(INPUTS.TABLEID, INPUTS.TABLENAME, FORMAT );

  NFILE = NFILE_AUTOSTORE;
  NVAR_TOT  = 0 ;
  TABLEID = INPUTS.TABLEID ;

  sprintf(VARNAME, "%s:C*%d", "CCID", MXCHAR_CCID );
  SNTABLE_ADDCOL(TABLEID, BLOCK, OUTPUT.CCID, VARNAME, USE4TEXT );

  for(ifile=0; ifile < NFILE; ifile++ ) {
    NVAR = SNTABLE_AUTOSTORE[ifile].NVAR ;
    for(ivar=0; ivar < NVAR; ivar++ ) {

      ICAST = SNTABLE_AUTOSTORE[ifile].ICAST_READ[ivar];
      if ( ICAST < 0 ) { continue ; }
      sprintf(CCAST, "%c", CCAST_TABLEVAR[ICAST] );
      if ( ICAST == ICAST_C ) { sprintf(CCAST,"%s*%d", CCAST, MXCHAR_CCID); }

      sprintf(VARNAME, "%s:%s", 
	      SNTABLE_AUTOSTORE[ifile].VARNAME[ivar], CCAST );

      if ( ICAST == ICAST_D ) 
	{ PTRVAR = &OUTPUT.VAL_D[NVAR_TOT];  }
      else if ( ICAST == ICAST_F ) 
	{ PTRVAR =& OUTPUT.VAL_F[NVAR_TOT]; }
      else if ( ICAST == ICAST_I ) 
	{ PTRVAR = &OUTPUT.VAL_I[NVAR_TOT]; }
      else if ( ICAST == ICAST_C ) 
	{ PTRVAR = OUTPUT.VAL_C[NVAR_TOT];  }

      SNTABLE_ADDCOL(TABLEID, BLOCK, PTRVAR, VARNAME, USE4TEXT );
      fflush(stdout);

      NVAR_TOT++ ;
    }
  }

  return ;

} //  end sntable_combine_init


// ============================================
void  sntable_combine_fill(int irow) {

  int ISTAT, iFile, NVAR, NVAR_TOT, ivar ,ICAST ;
  double DVAL;
  char CCID[MXCHAR_CCID], CVAL[40], *VARNAME ;
  //  char fnam[] = "sntable_combine_fill" ;

  // ------------- BEGIN ------------
  
  // get CCID from first file
  sprintf(CCID, "%s", SNTABLE_AUTOSTORE[0].CCID[irow] );
  sprintf(OUTPUT.CCID, "%s", CCID); // load global for output table
  NVAR_TOT = 0 ;

  if ( (irow%NROW_STDOUT_UPDATE) == 0 ) {
    printf("\t Update table row %7d of %7d  (CID=%s)\n", 
	   irow, OUTPUT.NROW, CCID );
    fflush(stdout);
  }

  for(iFile=0; iFile < NFILE_AUTOSTORE ; iFile++ ) {
    NVAR = SNTABLE_AUTOSTORE[iFile].NVAR ;

    for(ivar=0; ivar < NVAR; ivar++ ) {
      VARNAME = SNTABLE_AUTOSTORE[iFile].VARNAME[ivar]; 
      ICAST   = SNTABLE_AUTOSTORE[iFile].ICAST_READ[ivar];

      DVAL = -3333.0 ; sprintf(CVAL,"NULL_COMBINE");
      SNTABLE_AUTOSTORE_READ(CCID, VARNAME, &ISTAT, &DVAL, CVAL );

      if ( ICAST == ICAST_C ) 
	{ sprintf(OUTPUT.VAL_C[NVAR_TOT], "%s ", CVAL) ;  }
      else if ( ICAST == ICAST_D ) 
	{ OUTPUT.VAL_D[NVAR_TOT] = DVAL ; }
      else if ( ICAST == ICAST_F ) 
	{ OUTPUT.VAL_F[NVAR_TOT] = (float)DVAL ; }
      else if ( ICAST == ICAST_I ) 
	{ OUTPUT.VAL_I[NVAR_TOT] = (int)DVAL ; }

      NVAR_TOT++ ;
    }
  }

  SNTABLE_FILL(INPUTS.TABLEID);

  return ;
} // end sntable_combine_fill


// =================================================
void sntable_combine_summary(void) {

  double t_min1, t_min2, t_sec1, t_sec2, t_rate;

  // get wall-time info
  t_end   = time(NULL);
  t_sec1  = (t_read - t_start) ; // total read time
  t_sec2  = (t_end  - t_read ) ; // total write time
  t_min1  = t_sec1 / 60.0 ;   // convert to minutes
  t_min2  = t_sec2 / 60.0 ;   // convert to minutes

  t_rate  = 0.0 ;
  if ( t_sec2 > 0.0  ) { t_rate = (double)OUTPUT.NROW/t_sec2 ; }

  printf("\n SUMMARY: \n");
  printf("   Time to read %d input tables: %.2f minutes. \n", 
	 NFILE_AUTOSTORE, t_min1 );
  printf("   Time to combine %d tables: %.2f minutes (%.0f rows/sec). \n",
	 NFILE_AUTOSTORE, t_min2, t_rate );
  printf("   Combined %s table is in %s \n",
	 INPUTS.TABLENAME, INPUTS.OUTFILE);
  printf("\n");
  fflush(stdout);

} // end sntable_combine_summary


// =================================================
void testRead(void) {

  // test-read
#define MXLIST_TEST 10
  char CCID[MXLIST_TEST][20], varName[MXLIST_TEST][60];
  char CVAL[40];
  int istat, NLIST, i ;
  double DVAL ;

  // --------- BEGIN -----------

  printf("\n #### TEST READ AUTOSTORE ##### \n" );
  NLIST=0;

  sprintf(CCID[NLIST],"8771181");  sprintf(varName[NLIST], "NNTEST_PROB_IA");
  NLIST++ ;

  sprintf(CCID[NLIST],"8771181");  sprintf(varName[NLIST], "x1");
  NLIST++ ;

  sprintf(CCID[NLIST],"8771181");  sprintf(varName[NLIST], "SNRMAX");
  NLIST++ ;

  sprintf(CCID[NLIST],"8771181");  sprintf(varName[NLIST], "FIELD");
  NLIST++ ;

  for(i=0; i < NLIST; i++ ) {
    DVAL=0; sprintf(CVAL,"NULL");
    SNTABLE_AUTOSTORE_READ(CCID[i], varName[i], &istat, &DVAL, CVAL );
    printf("   CCID = %10.10s   %16.16s = %.4f(%s)  (istat=%d)\n", 
	   CCID[i], varName[i], DVAL, CVAL, istat);
    fflush(stdout);
  }

  return ;

} // end testRead
