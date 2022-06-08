/* =================================================
  Created Jan 9 2017
  Apply NN to data sample. 
   Output is FITRES_formatted file with list of
      CID ITYPE, NCELL, PROB_IA
  for each event in data sample


  USAGE: 
   nearnbr_apply.exe \
     -inFile_data   <inFile_data>
     -inFile_MLpar  <inFile_MLpar>  !formatted output of nearnbr_maxFoM.exe
     -outFile       <file>
    
     -varName_prob <varName>  ! default = NN_PROB_IA
     -nchop         <nchop>   ! default = 15
     -nproc         <nproc>   ! Num events to process; default=0=all


  Feb 7 2017: 
    + remove -inFile_sim and read this from inFile_MLpar.
    + read VARNAME_TRUE from inFile_MLpar (no longer hard-wired)

 Mar 2 2017:
    + fix NN_PROB_IA calculation in nearnbr_apply_exec();
      now works even if ITYPE_BEST<0.

 Mar 9 2020
   in read_NNpar, replace zPHOT with zHD.

 ==================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <sys/stat.h>

#include "sntools_nearnbr.h"
#include "sntools.h"
#include "sntools_output.h"


struct {
  char inFile_data[MXCHAR_FILENAME];
  char inFile_MLpar[MXCHAR_FILENAME];
  char outFile[MXCHAR_FILENAME];
  char varName_prob[60];   // name of output colum with SNIa prob
  
  int  NCHOP, NPROC;
  int  WFALSE;
} INPUTS ;


#define VARNAME_WFALSE  "WFALSE"
#define SEPMAX_PREFIX   "SEPMAX_" 
#define OUT_TABLEID     TABLEID_FITRES
#define OUT_TABLENAME   TABLENAME_FITRES
#define ITYPE_BEST_1A   1  

#define MXVAR_SEPMAX 10
int   NVAR_SEPMAX ;
char  VARNAME_SEPMAX[MXVAR_SEPMAX][20];
float VALUE_SEPMAX[MXVAR_SEPMAX];
char  VARLIST_SEPMAX[1000];

char  VARNAME_TRUETYPE[MXCHAR_VARNAME];
char  TRAIN_FILENAME[MXCHAR_FILENAME] ;
float TRAIN_SCALE_NON1A;   // June 2019  
int   TRUETYPE_SNIa;

int IFILE_DATA, IFILE_SIM ;
char msgerr1[80], msgerr2[80];

int NPROC_TOT;

struct {
  char  CCID[MXCHAR_CCID];
  int   CIDint ;
  int   ITYPE_BEST, NCELL ;
  float PROB_IA ;
} NNRESULTS ;

time_t tinit_start, tinit_end, tloop_start, tloop_end;

// =======================================
void  parse_args(int argc, char **argv) ;
void  read_NNpar(void);
void  read_data(void);
void  nearnbr_apply_init(void);
void  nearnbr_apply_exec(int ievt);

void  open_outFile(void);
void  time_summary(void);

// ==================================
int main(int argc, char **argv) {

  int i,i1;
  // --------------- BEGIN --------------------

  tinit_start = time(NULL) ;
  parse_args(argc,argv);

  read_NNpar(); // read sepmax variable names and values
  read_data();  // read data
  nearnbr_apply_init(); // read large sim
  
  tinit_end = time(NULL) ;

  open_outFile();

  printf("\n Apply NN to each data event . . . \n");

  tloop_start = time(NULL) ;
  int NROW = SNTABLE_AUTOSTORE[IFILE_DATA].NROW;
  if ( INPUTS.NPROC > 0 && INPUTS.NPROC < NROW ) 
    { NROW = INPUTS.NPROC; }

  NPROC_TOT = 0 ;
  for(i=0; i < NROW; i++ )  { 
    nearnbr_apply_exec(i);  
    NPROC_TOT++ ;
    // check screen update
    i1 = i+1;
    if ( (i1%1000)==0 || i1==NROW )
      { printf("\t Process event %6d of %d \n", i1,NROW); fflush(stdout); }
  }

  tloop_end = time(NULL) ;

  TABLEFILE_CLOSE(INPUTS.outFile);
  printf("\n Done with NN on %d events. \n", NROW);

  time_summary();

  return(0) ;

} // end main


// ====================
void parse_args(int NARG, char **argv) {

  int i ;
  char fnam[] = "parse_args" ;

  // ------------ BEGIN ----------

  // print full command (Sep 9 2015)
  printf("   Full command: ");
  for ( i = 0; i < NARG ; i++ ) { printf("%s ", argv[i]); }
  printf("\n\n"); fflush(stdout);

  INPUTS.inFile_data[0]  = 0 ;
  INPUTS.inFile_MLpar[0] = 0 ;
  sprintf(INPUTS.varName_prob,"NN_PROB_IA");
  sprintf(INPUTS.outFile, "out_NNresults.text" );
  INPUTS.WFALSE = 1 ;
  INPUTS.NCHOP  = 15 ;
  INPUTS.NPROC  = 0 ;  // 0 --> all

  if ( NARG < 2 ) {
    sprintf(msgerr1,"Must give 3 input files as arguments:");
    sprintf(msgerr2,"-inFile_data   -inFile_MLpar" );
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }

  for(i=1; i < NARG; i++ ) {

    if ( strcmp_ignoreCase(argv[i],"-inFile_data") == 0 ) 
      { sscanf(argv[i+1], "%s", INPUTS.inFile_data) ; }

    if ( strcmp_ignoreCase(argv[i],"-inFile_MLpar") == 0 ) 
      { sscanf(argv[i+1], "%s", INPUTS.inFile_MLpar) ; }

    if ( strcmp_ignoreCase(argv[i],"-outFile") == 0 ) 
      { sscanf(argv[i+1], "%s", INPUTS.outFile) ; }

    if ( strcmp_ignoreCase(argv[i],"-varName_prob") == 0 ) 
      { sscanf(argv[i+1], "%s", INPUTS.varName_prob) ; }

    if ( strcmp_ignoreCase(argv[i],"-nchop") == 0 ) 
      { sscanf(argv[i+1], "%d", &INPUTS.NCHOP) ; }

    if ( strcmp_ignoreCase(argv[i],"-nproc") == 0 ) 
      { sscanf(argv[i+1], "%d", &INPUTS.NPROC) ; }

  } 

  return ;

} // end parse_args


// ========================================
void  read_NNpar(void) {

  int  irow, Wfalse, optMask, ISTAT, Found_WFALSE, ivar, LEN_SEPMAX ;
  double DVAL;  char CVAL[40];
  char *f = INPUTS.inFile_MLpar ;
  char CROW[4], sepmaxVar[20], *PTRVAR, tmpVar[20], tmpVarList[20] ;
  char tableName[] = "NNsepmax" ;
  char fnam[]      = "read_NNpar" ;

  // ----------- BEGIN -----------

  optMask = 1; 
  SNTABLE_AUTOSTORE_INIT(f, tableName, "ALL", optMask);
  LEN_SEPMAX = strlen(SEPMAX_PREFIX);

  // examine VARNAMES to find SEPMAX_[varName]
  NVAR_SEPMAX = 0;
  VARLIST_SEPMAX[0] = 0 ;

  for(ivar=0; ivar < SNTABLE_AUTOSTORE[0].NVAR; ivar++ ) {
    PTRVAR = SNTABLE_AUTOSTORE[0].VARNAME[ivar] ;
    strncpy(sepmaxVar, PTRVAR, LEN_SEPMAX);   sepmaxVar[LEN_SEPMAX] = 0 ;

    if ( strcmp(sepmaxVar,SEPMAX_PREFIX) == 0 ) {

      sprintf(tmpVar, "%s", &PTRVAR[LEN_SEPMAX]);
      sprintf(VARNAME_SEPMAX[NVAR_SEPMAX], "%s", tmpVar );
      NVAR_SEPMAX++ ;      

      // variable name in read-list is same, except for z->zHD
      sprintf(tmpVarList, "%s", tmpVar);
      if ( strcmp(tmpVar,"z"    )==0 ) {  sprintf(tmpVarList,"zHD");  }
      if ( strcmp(tmpVar,"zPHOT")==0 ) {  sprintf(tmpVarList,"zHD");  }

      if ( NVAR_SEPMAX==1 ) 
	{ sprintf(VARLIST_SEPMAX,"%s", tmpVarList ); }
      else
	{ strcat(VARLIST_SEPMAX,","); strcat(VARLIST_SEPMAX,tmpVarList); }
    }
  }

  // ---------------
  Found_WFALSE = 0 ;

  // look for line with correct Wfalse and store sepmax values
  for(irow=0; irow < SNTABLE_AUTOSTORE[0].NROW; irow++ ) {
    sprintf(CROW, "%s", SNTABLE_AUTOSTORE[0].CCID[irow] );

    SNTABLE_AUTOSTORE_READ(CROW, VARNAME_WFALSE, &ISTAT, &DVAL, CVAL);
    Wfalse = (int)DVAL;

    if ( ISTAT != 0 ) {
      sprintf(msgerr1,"Could not find required '%s'  column.",VARNAME_WFALSE);
      sprintf(msgerr2,"Check %s", f);
      errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
    }

    if ( Wfalse != INPUTS.WFALSE ) { continue ; }

    Found_WFALSE = 1 ;

    for(ivar=0; ivar < NVAR_SEPMAX; ivar++ ) {
      sprintf(sepmaxVar, "%s%s", SEPMAX_PREFIX, VARNAME_SEPMAX[ivar] );
      SNTABLE_AUTOSTORE_READ(CROW, sepmaxVar, &ISTAT, &DVAL, CVAL);
      VALUE_SEPMAX[ivar] = (float)DVAL ;
      printf("\t Found SEPMAX Variable  %3s = %.4f \n", 
	     VARNAME_SEPMAX[ivar], VALUE_SEPMAX[ivar] );
    }
  }

  printf("\t VARLIST = '%s' \n", VARLIST_SEPMAX );

  if ( Found_WFALSE == 0 ) {
    sprintf(msgerr1,"Could not find WFALSE=%d", INPUTS.WFALSE);
    sprintf(msgerr2,"Check %s", f);
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );    
  }



  // re-open NNpar file and read VARNAME_TRUETYPE and the
  // name of the training file.
  char c_get[MXCHAR_FILENAME];
  FILE *fp = fopen(f,"rt");
  
  TRAIN_SCALE_NON1A = 1.0; // default value if not in file.
  TRUETYPE_SNIa     = 1 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"VARNAME_TRUE:") ==  0 ) 
      { readchar(fp, VARNAME_TRUETYPE ); }

    if ( strcmp(c_get,"TRAIN_FILENAME:") == 0 ) 
      { readchar(fp, TRAIN_FILENAME ); ENVreplace(TRAIN_FILENAME,fnam,1); }

    if ( strcmp(c_get,"TRAIN_SCALE_NON1A:") == 0 ) 
      { readfloat(fp, 1, &TRAIN_SCALE_NON1A ); }

    if ( strcmp(c_get,"TRUETYPE_SNIa:") == 0 ) 
      { readint(fp, 1, &TRUETYPE_SNIa ); }

  }
  fclose(fp);
 

  return ;

} // end read_NNpar

// ===============================================
void read_data(void) {

  char tableName[] = "FITRES" ;
  char *varList    = VARLIST_SEPMAX ;
  int optMask = 1+4 ; // 1=print, 4=next file
  //  char fnam[]      = "read_data";
  // --------------- BEGIN --------------

  SNTABLE_AUTOSTORE_INIT(INPUTS.inFile_data, tableName, varList, optMask);
  IFILE_DATA = NFILE_AUTOSTORE-1 ;

  return ;

} // end read_data


// ============================
void  nearnbr_apply_init(void) {

  int ivar ;
  char  *VARNAME ;
  double SEPMAX[8] ;
  // ------------- BEGIN -----------

  NEARNBR_INIT();
  NEARNBR_SET_TRAINFILE(TRAIN_FILENAME,  TRAIN_SCALE_NON1A) ;    
  NEARNBR_SET_TRUETYPE(VARNAME_TRUETYPE, TRUETYPE_SNIa    ) ; 

  for(ivar=0; ivar < NVAR_SEPMAX; ivar++ ) {
    VARNAME  = VARNAME_SEPMAX[ivar] ;
    SEPMAX[0]   = (double)VALUE_SEPMAX[ivar] ;  // min
    SEPMAX[1]   = (double)VALUE_SEPMAX[ivar] ;  // = max
    SEPMAX[2]   = (double)0.0 ; // binSize=0 --> apply (do NOT train)
    NEARNBR_SET_SEPMAX(VARNAME,SEPMAX) ; 
  }

  NEARNBR_CELLMAP_INIT(INPUTS.NCHOP); 

  NEARNBR_INIT2(1) ;

} // end nearnbr_apply_init


// ==================================
void nearnbr_apply_exec(int ievt) {

  int ivar ;
  double DVAL; char CCID[20], *VARNAME ;
  //  char fnam[] = "nearnbr_apply_exec" ;

  // --------------- BEGIN -----------------

  sprintf(CCID, "%s", SNTABLE_AUTOSTORE[IFILE_DATA].CCID[ievt]);

  for(ivar=0; ivar < NVAR_SEPMAX; ivar++ ) {
    VARNAME  = VARNAME_SEPMAX[ivar] ;

    // xxx slow SNTABLE_AUTOSTORE_READ(CCID, VARNAME, &ISTAT, &DVAL, CVAL);
    DVAL = SNTABLE_AUTOSTORE[IFILE_DATA].DVAL[ivar][ievt];

    NEARNBR_LOADVAL(CCID,VARNAME,DVAL) ;
  }

  NEARNBR_APPLY(CCID);

  int ITYPE_BEST, NTYPE, ITYPE_LIST[10], NCELL_TRAIN_LIST[10] ;
  NEARNBR_GETRESULTS(CCID, &ITYPE_BEST, &NTYPE, ITYPE_LIST, NCELL_TRAIN_LIST );

  // -----------------------------------
  // convert to NN_PROB_Ia and update outFile . . . .

  // ACCOUNT FOR SCALE_NON1A !!!
  int NCELL_TOT=0, NCELL_1A=0, i;
  for(i=0; i < NTYPE; i++ ) {
    NCELL_TOT += NCELL_TRAIN_LIST[i];

    if ( ITYPE_LIST[i] == ITYPE_BEST_1A ) 
      {  NCELL_1A += NCELL_TRAIN_LIST[i]; }
    
  }

  
  sprintf(NNRESULTS.CCID, "%s", CCID);
  NNRESULTS.CIDint     = ievt ;
  NNRESULTS.ITYPE_BEST = ITYPE_BEST ;
  NNRESULTS.NCELL      = NCELL_TOT ;
  if ( NCELL_TOT > 0 )
    { NNRESULTS.PROB_IA  = (float)NCELL_1A / (float)NCELL_TOT; }
  else
    { NNRESULTS.PROB_IA  = 0.0 ; }

  SNTABLE_FILL(OUT_TABLEID);

  return ;

} // end nearnbr_apply_exec


// ==============================
void  open_outFile(void) {

  char BLOCK[] = "NEARNBR" ;
  char textFormat[] = "key" ;
  char varDef[60];
  int USE4TEXT = 1; 
  //  char fnam[]  = "open_outFile";

  // ------------ BEGIN ----------

  TABLEFILE_OPEN(INPUTS.outFile, "new");
  CDTOPDIR_OUTPUT();

  SNTABLE_CREATE(OUT_TABLEID, OUT_TABLENAME, textFormat );

  // define columns to fill & write

  SNTABLE_ADDCOL(OUT_TABLEID, BLOCK, NNRESULTS.CCID,
		 "CCID:C*20",  USE4TEXT );

  SNTABLE_ADDCOL(OUT_TABLEID, BLOCK, &NNRESULTS.CIDint,
		 "CIDint:I", 0 );

  SNTABLE_ADDCOL(OUT_TABLEID, BLOCK, &NNRESULTS.ITYPE_BEST,
		 "NN_ITYPE:I",  USE4TEXT );

  SNTABLE_ADDCOL(OUT_TABLEID, BLOCK, &NNRESULTS.NCELL,
		 "NN_NCELL:I",  USE4TEXT );

  sprintf(varDef, "%s:F", INPUTS.varName_prob);
  SNTABLE_ADDCOL(OUT_TABLEID, BLOCK, &NNRESULTS.PROB_IA,
		 varDef, USE4TEXT );

  return;

} // end open_outFile



void  time_summary(void) {

  //  int NROW = SNTABLE_AUTOSTORE[IFILE_DATA].NROW;

  double tinit = tinit_end-tinit_start ;
  double tloop = tloop_end - tloop_start ;
  double rate  = (double)NPROC_TOT / (tloop+1.0E-9) ;

  printf("\n Summary of processing time: \n");
  printf("\t Init  time: %.2f seconds \n", tinit );
  printf("\t Apply time: %.2f minutes  (%.2f/sec)\n", 
	 tloop/60.0, rate );
 
  fflush(stdout);

} // end time_summary
