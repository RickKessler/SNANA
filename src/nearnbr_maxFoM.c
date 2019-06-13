/***********************************************************
  Created April 2014 by R.Kessler


  Stand-alone program to analyze output based on sntools_nearnbr.c[h]
  and do the NN training. Writes results to stdout.

  HID 800  variable name with true type (added Jan 2017)
  HID 801  Number of merged files
  HID 802  TrueType vs. sparse type-index (NBIN = number of possible types_
  HID 803  SEPMAX bins
  HID 810 to 810 + NTYPE-1 :  tabulated results for each TrueType
  HID 840  tot number of training events (for info message)

 Usage:
    nearnbr_maxFoM.exe <inpFile>
        inpFile = hbook or root file output from fitting program

 Optional Usage
    nearnbr_maxFoM.exe <inpFile>  --truetype <type> -outFile <outFile>
       [print results for only one truetype]
  
    nearnbr_maxFoM.exe <inpFile>  --truetype <type>  VBOSE
       [print verbose mode]

    nearnbr_maxFoM.exe <inpFile>  --truetype <type>  -wfalse_pipeline <wfalse>
       [used by NEARNBR_pipeline to make greppable summary]



                   HISTORY
           ~~~~~~~~~~~~~~~~~~~~~~~~~

  May 6 2014: new command-line options 'VBOSE' and 'VARDEF'

  Sep 9 2015: change PIPELINE arg to -wfalse_pipeline <wfalse>
              to select which wfalse.

              Print full command.

  Jun 21 2016: fix bug computing Eff in getFoM().

  Jan 09 2017: read VARNAME_TRUETYPE form histogram title so that
               we don't have to pass it later to apply training.

  Apr 11 2019: read HID 840 with total number of training events

***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "sntools.h"
#include "sntools_output.h"
#include "sntools_nearnbr.h"

//#define MXTRUETYPE      100 
#define MXVAR           10

#define HID_Nmerge      801
#define HID_TypeList    802
#define HID_SEPMAXbins  803
#define HID_UntrainedPurity  805
#define HIDOFF_SEPMAXresults 810
#define HID_NEVTOT           840
#define HID_TRAIN_FILENAME   850
#define HID_varName_trueType 860

#define N_WFALSE        3
int WFALSE_LIST[N_WFALSE] = { 1,3,5  } ;

struct INPUTS {
  char NEARNBR_INFILE[200] ;  // name of hbook or root file
  char NEARNBR_OUTFILE[200];  // formatted outFile (Jan 2017)

  int  IFILETYPE ;            // specifies hbook or root or ...
  int  TRUETYPE ;       // process only this TRUETYPE (default is all)

  int FLAG_VBOSE ;     // flag for more verbose output
  int FLAG_VARDEF ;    // flag to write &NNINP lines 
  int WFALSE_PIPELINE ;  // when called by NEARNBR_pipeline.pl
  
  // below are inputs read from INFILE histograms
  int Nmerge ;
  int NTrueType ;
  int TrueType[MXTRUETYPE];
  double UntrainedPurity[MXTRUETYPE] ;

  int    NVAR ;
  int    NBIN_SEPMAX ;
  char   **VARNAME ;
  double **SEPMAX ;  // sepmax values vs. IVAR:ISEP
  double *SEPMAX_RANGE[2] ; // min,max sepmax vs. IVAR

  int  ***NTRAIN ; // vs. ispType(True) : ispType(train) : isep 
  int  **NFAIL ;   // vs. ispType(True) : isep

  int *NSIMTYPE ; // number of simulated types vs. true type (6/21/2016)
  
  int NEVTOT;  // total number of training events (HID 840)

} INPUTS ;

FILE *FP_OUT;
int   IROW_OUT;
char  NNINP_VARDEF[MXTRUETYPE][N_WFALSE][200] ;
char  VARNAME_TRUETYPE[40] ;
char  TRAIN_FILENAME[MXCHAR_FILENAME];
float TRAIN_NON1A_SCALE;
int   TRUETYPE_SNIa ;

char msgerr1[80], msgerr2[80];

// define logical flag to TRAIN on each TRUETYPE
int DOTRAIN_TRUETYPE[MXTRUETYPE] ;


struct {
  double  PURITY_noNN ;
  double  EFF[N_WFALSE];
  double  PPURITY[N_WFALSE]; // pseudo-purity
  double  FOM[N_WFALSE];
  char    WARN_BOUNDARY[N_WFALSE][40]; // asterisk(s) for boundary warning
} NNFOM[MXTRUETYPE] ;

// =====================
void  parse_args(int argc, char **argv) ;
void  open_inFile(void);
void  open_outFile(void);

void  RDNN_VARNAME_TRUETYPE(void); 
void  RDNN_Nmerge(void);
void  RDNN_TypeList(void);
void  RDNN_SEPMAXbins(void);
void  RDNN_NTRAIN(int iTypeTrue) ; 
void  RDNN_UntrainedPurity(void);
void  RDNN_TRAIN_FILENAME(void);
void  RDNN_NEVTOT(void);

void  malloc_train(void) ;
void  dump_NTRAIN(int iType_true, int iType_train, int isep ) ;

void   optimizeNN(int iTypeTrue, int ifalse);
double getFoM(int isep, int iTypeTrue, double Wfalse, 
	      double *eff, double *purity) ;

void dumpLine_UntrainedPurity(int iTypeTrue);
void dump_NNINP_VARDEF(int iTypeTrue) ;
void dump_forPipeline(int iTypeTrue) ;

// ==================================
int main(int argc, char **argv) {

  parse_args(argc,argv);

  // open input table file
  open_inFile();

  // read histograms (from hbook or root)
  RDNN_VARNAME_TRUETYPE();  // Jan 2017
  RDNN_Nmerge();
  RDNN_TypeList();
  RDNN_SEPMAXbins();
  RDNN_UntrainedPurity();
  RDNN_TRAIN_FILENAME();
  RDNN_NEVTOT();           // read NEVTOT (Apr 11 2019)

  // allocate memory for training
  malloc_train();

  // read training results (from hbook or root)
  int iTypeTrue ;
  for(iTypeTrue=0; iTypeTrue < INPUTS.NTrueType; iTypeTrue++ ) {  
    RDNN_NTRAIN(iTypeTrue) ; 
  }

  
  // ----------------------------------
  //  pseudo-purity = Ntrue/( Ntrue + Wfalse*Nfalse)
  //
  
  optimizeNN(-9,0) ; // make table header

  open_outFile();

  int ifalse ;
  for(iTypeTrue=0 ; iTypeTrue < INPUTS.NTrueType; iTypeTrue++ ) {

    if ( DOTRAIN_TRUETYPE[iTypeTrue] == 0 ) { continue ; }

    // print untrained purity for reference
    dumpLine_UntrainedPurity(iTypeTrue);

    for(ifalse = 0; ifalse < N_WFALSE; ifalse ++ ) 
      {  optimizeNN(iTypeTrue,ifalse);  }

    printf("\t\t\t\t\t +/- => hi/low SEPMAX boundary\n");
    printf("\n"); fflush(stdout);

    if ( INPUTS.FLAG_VARDEF   ) { dump_NNINP_VARDEF(iTypeTrue); }

    if ( INPUTS.WFALSE_PIPELINE > 0 ) { dump_forPipeline(iTypeTrue); }
  }
  
  
  // xxx obsolete  TABLEFILE_CLOSE(INPUTS.NEARNBR_INFILE);
  fclose(FP_OUT);
  printf(" Formatted text output: %s \n", INPUTS.NEARNBR_OUTFILE);

  return(0);

} // end of main


// ====================
void parse_args(int NARG, char **argv) {

  int i ;
  char fnam[] = "parse_args" ;

  // ---------------- BEGIN ------------- 
  if ( NARG < 2 ) {
    sprintf(msgerr1,"Must give FITRES-TABLE input file as arg,");
    msgerr2[0] = 0 ;
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }

  // print full command (Sep 9 2015)
  printf("   Full command: ");
  for ( i = 0; i < NARG ; i++ ) { printf("%s ", argv[i]); }
  printf("\n\n"); fflush(stdout);

  INPUTS.TRUETYPE      = -99;
  INPUTS.FLAG_VBOSE    = 0 ;
  INPUTS.FLAG_VARDEF   = 1 ;
  INPUTS.WFALSE_PIPELINE = 0 ;
  INPUTS.NEARNBR_OUTFILE[0] = 0 ;
  sprintf(INPUTS.NEARNBR_INFILE,"%s", argv[1]) ; 

  for(i=2; i < NARG; i++ ) {

    if ( strcmp_ignoreCase(argv[i],"--truetype") == 0 ) 
      { sscanf(argv[i+1], "%d", &INPUTS.TRUETYPE) ; }

    if ( strcmp_ignoreCase(argv[i],"-truetype") == 0 ) 
      { sscanf(argv[i+1], "%d", &INPUTS.TRUETYPE) ; }

    if ( strcmp_ignoreCase(argv[i],"truetype") == 0 ) 
      { sscanf(argv[i+1], "%d", &INPUTS.TRUETYPE) ; }

    if ( strcmp_ignoreCase(argv[i],"-outFile") == 0 ) 
      { sscanf(argv[i+1], "%s", INPUTS.NEARNBR_OUTFILE) ; }

    if ( strcmp_ignoreCase(argv[i],"VBOSE") == 0 ) 
      { INPUTS.FLAG_VBOSE=1; }

    if ( strcmp_ignoreCase(argv[i],"NEARNBR_SEPMAX_VARDEF") == 0 ) 
      { INPUTS.FLAG_VARDEF = 1; }

    if ( strcmp_ignoreCase(argv[i],"-wfalse_pipeline") == 0 ) 
      { sscanf(argv[i+1], "%d", &INPUTS.WFALSE_PIPELINE) ; }
  }

  printf("\n# ================================================ \n");
  printf(" INFILE = '%s' \n", INPUTS.NEARNBR_INFILE );
  fflush(stdout);

} // end of parse_args



// ==================================
void open_inFile(void) {
  char *inFile, copt[40] ;
  char fnam[] = "open_inFile" ;
  // --------------- BEGIN -----------

  TABLEFILE_INIT();

  inFile   = INPUTS.NEARNBR_INFILE ;

  if ( INPUTS.FLAG_VBOSE  ) 
    {  sprintf(copt,"read"); }
  else
    {  sprintf(copt,"read q"); }

  
  INPUTS.IFILETYPE = TABLEFILE_OPEN(inFile,copt);

  if ( INPUTS.IFILETYPE < 0 ) {
    sprintf(msgerr1,"Could not open table file: ");
    sprintf(msgerr2,"%s", inFile);
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }

} // end of open_inFile


// ==========================================
void  open_outFile(void) {

  int NVAR = INPUTS.NVAR  + 7 ;
  int ivar ;
  //  char fnam[] = "open_outFile" ;

  // ------------- BEGIN --------------

  if ( strlen(INPUTS.NEARNBR_OUTFILE) == 0 ) 
    { sprintf(INPUTS.NEARNBR_OUTFILE,"nearnbr_maxFoM.out"); }

  FP_OUT = fopen(INPUTS.NEARNBR_OUTFILE,"wt");
  IROW_OUT = 0 ;

  fprintf(FP_OUT, "MLAPPLY_CODENAME:  nearnbr_apply.exe\n\n" ) ;

  fprintf(FP_OUT, "TRAIN_FILENAME: %s \n\n", TRAIN_FILENAME );
  fprintf(FP_OUT, "NEVTOT_TRAIN:   %d \n",   INPUTS.NEVTOT );
  fprintf(FP_OUT, "NBIN_SEPMAX:    %d \n",   INPUTS.NBIN_SEPMAX );
  fprintf(FP_OUT, "VARNAME_TRUE:   %s \n",   VARNAME_TRUETYPE ); 
  fprintf(FP_OUT, "TRUETYPE_SNIa:  %d \n",   TRUETYPE_SNIa ); 
  fprintf(FP_OUT, "TRAIN_NON1A_SCALE: %.2f \n",   TRAIN_NON1A_SCALE ); 

  // print min,max training range for each variable
  for(ivar=0; ivar < INPUTS.NVAR; ivar++ ) {
    fprintf(FP_OUT,"TrainRange(%2s): %8.4f to %8.4f \n"
	    ,INPUTS.VARNAME[ivar]
	    ,INPUTS.SEPMAX_RANGE[0][ivar] 
	    ,INPUTS.SEPMAX_RANGE[1][ivar] );
  }


  fprintf(FP_OUT, "\n" ) ;
  fprintf(FP_OUT, "NVAR: %d\n", NVAR);
  fprintf(FP_OUT, "VARNAMES: ROW TRUETYPE WFALSE  EFF PURITY FOM     ISEP  ");
  for(ivar=0; ivar < INPUTS.NVAR; ivar++ ) 
    { fprintf(FP_OUT,"SEPMAX_%s ", INPUTS.VARNAME[ivar]  );  }
  fprintf(FP_OUT, "\n");
  fflush(FP_OUT);

  return ;

} // end open_outFile

// ==========================================
void  RDNN_VARNAME_TRUETYPE(void) {

  // Feb 07 2017
  // Read name of variable with trueType (from hist title)
  int HID = HID_varName_trueType ;  
  int NB; double XMIN, XMAX, X;

  sprintf(VARNAME_TRUETYPE, "unknown");
  SNHIST_RDBINS(1, HID, VARNAME_TRUETYPE, &NB, &XMIN, &XMAX);
  trim_blank_spaces(VARNAME_TRUETYPE);
  //  printf(" xxx VARNAME_TRUETYPE = '%s' \n", VARNAME_TRUETYPE);

  // Jun 12 2019: read true SNIa type from content
  SNHIST_RDCONT(1, HID, NB, &X);
  TRUETYPE_SNIa = (int)X;

  return ;

} // end RDNN_VARNAME_TRUETYPE

// ==========================================
void RDNN_NEVTOT(void) {
  int HID = HID_NEVTOT ;
  int NBRD;
  double X, XMIN, XMAX ;
  char CTIT[200];
  // ---------- BEGIN -----------
  SNHIST_RDBINS(1, HID,CTIT, &NBRD, &XMIN, &XMAX);
  SNHIST_RDCONT(1, HID, NBRD, &X);
  INPUTS.NEVTOT = (int)X ;
  return ;
} // end RDNN_NEVTOT

// ==========================================
void RDNN_TRAIN_FILENAME(void) {

  // Created Feb 7 2017
  // Apr 11 2019: fix unitialized bug by setting strTmp[i][0]=0.

  char strTmp[NSPLIT_TITLE][MXCHAR_FILENAME];
  int HID = HID_TRAIN_FILENAME ;  
  int NB,i ; double XMIN, XMAX, X;

  TRAIN_FILENAME[0] = 0;

  for(i=0; i < NSPLIT_TITLE ; i++ ) {
    strTmp[i][0] = 0 ;
    SNHIST_RDBINS(1, HID+i, strTmp[i], &NB, &XMIN, &XMAX);
    trim_blank_spaces(strTmp[i]);
    //    printf(" xxx %d : strTmp = '%s'  (len=%d)\n", 
    //	   i, strTmp[i], strlen(strTmp[i]) ); fflush(stdout);

    if ( strlen(strTmp[i]) > 0 ) { strcat(TRAIN_FILENAME, strTmp[i] ); }
  }


  // Jun 12 2019: read NONIA_SCALE from y-axis content
  SNHIST_RDCONT(1, HID_TRAIN_FILENAME, NB, &X);
  TRAIN_NON1A_SCALE = (float)X;

  return ;
} // end RDNN_TRAIN_FILENAME

// ==========================================
void  RDNN_Nmerge(void) {

  int HID = HID_Nmerge ;

  char CTIT[80];
  int  NB, Nmerge ;
  double XMIN, XMAX, X;
  char fnam[] = "RDNN_Nmerge" ;

  // --------- BEGIN -----------

  SNHIST_RDBINS(1, HID,CTIT, &NB, &XMIN, &XMAX);
  SNHIST_RDCONT(1, HID, NB,  &X);
  Nmerge = (int)X ;

  if ( INPUTS.FLAG_VBOSE ) 
    {  printf(" Read HID %d : \n\t Nmerge = %d \n", HID, Nmerge ); }

  if ( Nmerge <= 0 || Nmerge > 9999 ) {
    printf("\t PRE-ABORT DUMP: \n");
    printf("   NB=%d  XMIN=%f  XMAX=%f\n", NB, XMIN, XMAX);
    printf("   CTIT = '%s' \n", CTIT);
    
    sprintf(msgerr1,"Invalid Nmerge=%d for HID=%d", Nmerge, HID);
    sprintf(msgerr2,"Nmerge must be at least 1." );
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }

  INPUTS.Nmerge = Nmerge ;  // store in global

  fflush(stdout);
  
  return ;

} // end of RDNN_Nmerge

// ==========================================
void  RDNN_TypeList(void) {

  // read hid 802 and fill INPUTS.NTrueType and .TrueType array.

  int  HID = HID_TypeList ;
  char   CTIT[80];
  int    NB, i, TrueType, DOFLAG ;
  double XMIN, XMAX, X[MXTRUETYPE];  

  SNHIST_RDBINS(1, HID,CTIT, &NB, &XMIN, &XMAX);
  SNHIST_RDCONT(1, HID, NB, X);

  if ( INPUTS.FLAG_VBOSE ) {  printf(" Read HID %d : \n", HID); }

  INPUTS.NTrueType = NB ;
  DOFLAG = 1 ; // default is to train on all TRUETYPEs

  for(i=0; i < NB; i++ ) {
    TrueType           = (int)X[i]/INPUTS.Nmerge ;
    INPUTS.TrueType[i] = TrueType ;

    if ( INPUTS.FLAG_VBOSE ) 
      { printf("\t TrueType[%2d] = %d \n", i, TrueType ); }

    if ( INPUTS.TRUETYPE >= 0 ) 
      { DOFLAG = ( INPUTS.TRUETYPE == TrueType) ; }

    DOTRAIN_TRUETYPE[i] = DOFLAG ;
  }

  fflush(stdout);

} // end of RDNN_TypeList

// ==========================================
void  RDNN_SEPMAXbins(void) {

  // read SEPMAX bins for each variable.

  int HID = HID_SEPMAXbins ;
  char   CTIT[80], *ptrtok ;
  int    NB[2], NXY, ivar, isep, j, NVAR, NBSEP ;
  double XMIN[2], XMAX[2], *CONTENTS, XNorm, SEPVAL ;

  SNHIST_RDBINS(2, HID, CTIT, NB, XMIN, XMAX);
  XNorm = (double)INPUTS.Nmerge ;

  NBSEP              = NB[0] ; // local var
  NVAR               = NB[1] ;
  INPUTS.NBIN_SEPMAX = NB[0] ; // global var

  INPUTS.NVAR        = NB[1] ;
  
  if ( INPUTS.FLAG_VBOSE ) {
    printf(" Read HID %d : \n", HID);
    printf("\t NBIN_SEPMAX = %d     NVAR=%d\n",  NBSEP, NVAR );
  }

  // allocat memory for variable names
  INPUTS.VARNAME = (char  **)malloc( sizeof(char  *) * NVAR ) ;
  INPUTS.SEPMAX  = (double**)malloc( sizeof(double*) * NVAR ) ;
  INPUTS.SEPMAX_RANGE[0]  = (double*)malloc( sizeof(double) * NVAR );
  INPUTS.SEPMAX_RANGE[1]  = (double*)malloc( sizeof(double) * NVAR );

  for(ivar=0; ivar < INPUTS.NVAR; ivar++ ) {
    INPUTS.VARNAME[ivar] = (char  *)malloc( sizeof(char)   * 60 ); 
    INPUTS.SEPMAX[ivar]  = (double*)malloc( sizeof(double) * NBSEP);
  }

  ptrtok = strtok(CTIT," ");
  
  for(ivar=0; ivar < NVAR; ivar++ ) {
    ptrtok = strtok(NULL," " );
    sprintf( INPUTS.VARNAME[ivar], "%s", ptrtok);

    if ( INPUTS.FLAG_VBOSE ) 
      { printf("\t VAR[%d] = '%s' \n", ivar, INPUTS.VARNAME[ivar] ); }
  }

  // read contents

  CONTENTS = (double*) malloc( sizeof(double) * NBSEP * NVAR ) ;
  NXY = NB[0] * NB[1] ;
  SNHIST_RDCONT(2, HID, NXY, CONTENTS);

  for(ivar=0; ivar < NVAR; ivar++ ) {
    for( isep=0; isep < NBSEP ; isep++ ) {
      j = ivar*NBSEP + isep;
      SEPVAL = CONTENTS[j]/XNorm ;
      INPUTS.SEPMAX[ivar][isep] = SEPVAL ;
    }
  }
  free(CONTENTS);

  char cline_first[200], cline_last[200];

  sprintf(cline_first, "\t SEPMAX(1st  bin): ");
  for(ivar=0; ivar < NVAR; ivar++ ) { 
    SEPVAL = INPUTS.SEPMAX[ivar][0] ;
    INPUTS.SEPMAX_RANGE[0][ivar] = SEPVAL ;
    sprintf(cline_first, "%s %.4f " , cline_first, SEPVAL );
  }


  sprintf(cline_last, "\t SEPMAX(last bin): " ); 
  for(ivar=0; ivar < NVAR; ivar++ ) { 
    SEPVAL = INPUTS.SEPMAX[ivar][NBSEP-1] ;
    INPUTS.SEPMAX_RANGE[1][ivar] = SEPVAL ;
    sprintf(cline_last, "%s %.4f " , cline_last, SEPVAL ); 
  }

  if ( INPUTS.FLAG_VBOSE ) { 
    printf("%s \n", cline_first );
    printf("%s \n", cline_last  );
    fflush(stdout);
  }


} // end of RDNN_SEPMAXbins


// ========================
void  malloc_train(void) {

  int NTYPE, NBSEP, itype, itype2, MEM, MTOT ;

  NTYPE = INPUTS.NTrueType ;
  NBSEP = INPUTS.NBIN_SEPMAX ;

  MTOT = 0 ;

  INPUTS.NTRAIN    = (int***) malloc ( sizeof(int**) * NTYPE ) ; 
  INPUTS.NFAIL     = (int **) malloc ( sizeof(int *) * NTYPE ) ;
  INPUTS.NSIMTYPE  = (int  *) malloc ( sizeof(int  ) * NTYPE ) ;

  for(itype = 0; itype < NTYPE; itype++ ) {  // true type    

    INPUTS.NSIMTYPE[itype]      = 0 ;
    
    MEM = sizeof(int*) * NTYPE ;      MTOT += MEM ;  
    INPUTS.NTRAIN[itype] = (int**) malloc ( MEM ) ; 

    MEM = sizeof(int ) * NBSEP ;    MTOT += MEM ; 
    INPUTS.NFAIL[itype]  = (int *) malloc ( MEM ) ; 

    for(itype2=0; itype2 < NTYPE; itype2++ ) {  // train type
      MEM = sizeof(int) * NBSEP ;      MTOT += MEM ;
      INPUTS.NTRAIN[itype][itype2] = (int*) malloc ( MEM );      
    }
  }

  double dmem = (double)MTOT/1.0E6 ;

  if ( INPUTS.FLAG_VBOSE ) {
    printf("\n Allocate %.3f MB memory for NTRAIN and NFAIL \n\n" , dmem) ;
    fflush(stdout);
  }

  return ;

} // end of malloc_train


// ========================================
void  RDNN_UntrainedPurity(void) {

  // Apr 15 2014
  // read/store untrained purities for reference;
  // these purities are not used for any calculations.

  int HID, NTrueType, i, NB ;
  double XMIN, XMAX, XNorm, tmpPurity[MXTRUETYPE] ;
  char CTIT[100];
  // ----------------- BEGIN ------------

  HID       =  HID_UntrainedPurity ;
  NTrueType =  INPUTS.NTrueType ;  
  XNorm     =  (double)INPUTS.Nmerge ;  

  for ( i=0; i < NTrueType; i++ )  { INPUTS.UntrainedPurity[i] = -9.0 ; }

  SNHIST_RDBINS(1, HID, CTIT, &NB, &XMIN, &XMAX ) ;
  SNHIST_RDCONT(1, HID, NTrueType,  tmpPurity   ) ;

  for ( i=0; i < NTrueType; i++ ) 
    { INPUTS.UntrainedPurity[i] = tmpPurity[i]/XNorm ; }

  return ;
  
} // end of   RDNN_UntrainedPurity

// ========================================
void RDNN_NTRAIN(int iTypeTrue) {

  // ispType is the sparse index for the True Type.
  // Vertical axis of HID is the trained sparse type 
  // that has the most training-SN with SQDIST<1.

  // ----------- BEGIN ----------
  
  int    HID = HIDOFF_SEPMAXresults + iTypeTrue ;
  int    NB[2], NBTOT, NTYPE, NBSEP, iTypeTrain, isep, j, NCONTENTS ;
  double XMIN[2], XMAX[2], *CONTENTS ;
  char   CTIT[80] ;

  NTYPE = INPUTS.NTrueType ;
  NBSEP = INPUTS.NBIN_SEPMAX ;
  NBTOT = (NTYPE+1) * NBSEP ;

  if ( INPUTS.FLAG_VBOSE ) 
    {  printf(" Read HID %d : \n", HID); }

  // we know the binning, but read it anyway to define hist.
  SNHIST_RDBINS(2, HID, CTIT, NB, XMIN, XMAX);

  // read contents. Note that itype goes from -1 to NTYPE-1,
  // where -1 indicates a failure to find a train-type.
  CONTENTS = (double*) malloc( sizeof(double) * NBTOT ) ;
  SNHIST_RDCONT(2, HID, NBTOT, CONTENTS);

  for(iTypeTrain=-1; iTypeTrain < NTYPE; iTypeTrain++ ) { 
    for(isep=0; isep < NBSEP; isep++ ) {
      
      j = (iTypeTrain+1)*NBSEP + isep;
      NCONTENTS = (int) CONTENTS[j];
      
      if ( iTypeTrain >= 0 ) 
	{ INPUTS.NTRAIN[iTypeTrue][iTypeTrain][isep] = NCONTENTS ; }
      else 
	{ INPUTS.NFAIL[iTypeTrue][isep] = NCONTENTS;  }

      if ( isep==0 ) { INPUTS.NSIMTYPE[iTypeTrue] += NCONTENTS ; }
      
    } // isep
  } // iTypeTrain
  

  //  isep = 0 ;  dump_NTRAIN(iTypeTrue, 0, isep);    

  printf(" NSIM[trueType=%d] = %6d  \n",
	 iTypeTrue,  INPUTS.NSIMTYPE[iTypeTrue] ) ;
  fflush(stdout);
  
  return ;
  
} // end of RDNN_NTRAIN

// ===============
void dump_NTRAIN(int iType_true, int iType_train, int isep ) {

  int N;

  if ( iType_train >= 0 ) 
    { N = INPUTS.NTRAIN[iType_true][iType_train][isep] ; }
  else
    { N = INPUTS.NFAIL[iType_true][isep] ; }

  
  printf("\t XXX DUMP NTRAIN[%2d][%2d][%4d] = %d \n",
	 iType_true, iType_train, isep, N);
  
  return ;

} // end of dump_NTRAIN



// ==============================================
void  optimizeNN(int iTypeTrue, int ifalse) {

  // April 2014
  // For input iTypeTrue and Wfalse, evaluate optimal (maximum)
  //         FoM = eff x purity 
  // among the SEPMAX bins, and print results to screen.
  // Purity = Ntrue/(Ntrue + Wfalse*Nfalse);
  // Eff    = Ntag(trueIa)/[ Ntrue(all) ]
  //
  // Wfalse < 0 is a flag to
  //  - compute/print purity before NN (where EFF=100%)
  //  - print table header
  //
  // WARNING: this function calls getFoM and prints results ...
  //          should probably untangle this spagetti.
  //

  int NBSEP, isep, isep_SAVE, ivar ;
  double FoM, Eff, Pur, FoM_SAVE, Eff_SAVE, Pur_SAVE, Wfalse ;

  char cWARN[4]; // warning if SEPMAX is on a boundary
  char dashLine[] =
    " --------------------------------------------------"
    "--------------------------" ;


  // ------------------ BEGIN ------------------

  Wfalse = (double)WFALSE_LIST[ifalse] ;

  if ( iTypeTrue == -9 ) {

    // print table header and return

    printf("\n");

    printf(" True               Pseudo-        "); 
    for(ivar=0; ivar <= INPUTS.NVAR; ivar++ )
      { printf("%10s", "Optimized" ); }
    printf("\n");

    printf(" Type Wfalse  Eff x Purity = FoM    isep ");    
    for(ivar=0; ivar < INPUTS.NVAR; ivar++ )
      {  printf("%6s-sep ", INPUTS.VARNAME[ivar] ); }
    printf("\n");

    printf(" %s \n", dashLine );
    fflush(stdout);
    return ;
  }


  // ----------------------------------------
  // determine optimzed SEPMAX values here
  NBSEP = INPUTS.NBIN_SEPMAX ;
  FoM_SAVE = Eff_SAVE = Pur_SAVE = 0.0 ;
  isep_SAVE = -9 ;
  for(isep=0; isep < NBSEP; isep++ ) {
    FoM = getFoM(isep, iTypeTrue, Wfalse, &Eff, &Pur); // return Eff,Pur
    if ( FoM > FoM_SAVE ) {
      FoM_SAVE  = FoM ;
      Eff_SAVE  = Eff ;
      Pur_SAVE  = Pur ;
      isep_SAVE = isep ;
    }
  }

  // --------------------------------------------
  // print optimized results to line in table.

  printf(" %4d  %.1f   %.3f x %.3f = %.3f   %5d ",
	 INPUTS.TrueType[iTypeTrue],
	 Wfalse, Eff_SAVE, Pur_SAVE, FoM_SAVE, isep_SAVE );

  // update formatted text file (Jan 2017)
  IROW_OUT++ ;
  fprintf(FP_OUT,"ROW:      %2d      %d       %d   ", 
	  IROW_OUT, INPUTS.TrueType[iTypeTrue], (int)Wfalse );
  fprintf(FP_OUT,"%5.3f %5.3f %5.3f  %5d ",
	  Eff_SAVE, Pur_SAVE, FoM_SAVE, isep_SAVE );

  // xyz
  for(ivar=0; ivar < INPUTS.NVAR; ivar++ ) 
    { fprintf(FP_OUT," %7.4f ", INPUTS.SEPMAX[ivar][isep_SAVE] );  }

  fprintf(FP_OUT,"\n"); fflush(FP_OUT);

  // ------- store results in global struct ----------
  NNFOM[iTypeTrue].PURITY_noNN = INPUTS.UntrainedPurity[iTypeTrue] ;
  NNFOM[iTypeTrue].EFF[ifalse]     = Eff_SAVE ;
  NNFOM[iTypeTrue].PPURITY[ifalse] = Pur_SAVE ;
  NNFOM[iTypeTrue].FOM[ifalse]     = FoM_SAVE ; 
  NNFOM[iTypeTrue].WARN_BOUNDARY[ifalse][0] = 0 ; // init boundary warning

  double VAL, VALMIN, VALMAX ;
  char   VARDEF_LINE[200], VARDEF_STRING[100] ;
  
  for(ivar=0; ivar < INPUTS.NVAR; ivar++ )  { 
    VAL    =  INPUTS.SEPMAX[ivar][isep_SAVE] ;
    VALMIN =  INPUTS.SEPMAX_RANGE[0][ivar] ;
    VALMAX =  INPUTS.SEPMAX_RANGE[1][ivar] ;

    if ( VAL == VALMIN  ) { // VAL is on a low SEPMAX boundary
      sprintf(cWARN,"-");
      strcat(NNFOM[iTypeTrue].WARN_BOUNDARY[ifalse],INPUTS.VARNAME[ivar]);
      strcat(NNFOM[iTypeTrue].WARN_BOUNDARY[ifalse],"- ");
    } 
    else if ( VAL == VALMAX ) { // VAL is on a high SEPMAX boundary
      sprintf(cWARN,"+");
      strcat(NNFOM[iTypeTrue].WARN_BOUNDARY[ifalse],INPUTS.VARNAME[ivar]);
      strcat(NNFOM[iTypeTrue].WARN_BOUNDARY[ifalse],"+ ");
    } 
    else
      { sprintf(cWARN," "); }

    printf("%9.4f%s", VAL, cWARN );


    // -- construct strinf for &NNINP
    sprintf(VARDEF_STRING,"%s %.3f", INPUTS.VARNAME[ivar], VAL);
    if ( ivar == 0 ) {
      sprintf(VARDEF_LINE,"%s",  VARDEF_STRING);      
    }
    else {
      sprintf(VARDEF_LINE,"%s %s",  VARDEF_LINE, VARDEF_STRING);
    }

  } // end of ivar
 
  printf("\n");

  sprintf(NNINP_VARDEF[iTypeTrue][ifalse],
	  "NEARNBR_SEPMAX_VARDEF = '%s' ", VARDEF_LINE );

  fflush(stdout) ;

  return ;
  
} // end of optimizeNN


// ====================================================
double  getFoM(int isep, int iTypeTrue, double Wfalse, 
	       double *Eff, double *Pur) {


  // Function returns FoM = Eff x Purity at this isep bin.
  // iTypeTrue is the sparse index for the true type.
  // *Eff and *Pur are the output efficiency and purity.
  //
  // Jun 21 2016: fix bug computing Eff ... need separate
  //              iType loops for Eff and Purity.
  
  int    NTRAIN, NTRUE_TOT, NTRUE_TYPE, NFAIL, iType ;
  double dNTRAIN, dNtrue, dNfalse, FoM, Eff_local, Pur_local ;
  int    NTYPE = INPUTS.NTrueType ;
  //  char  fnam[] = "getFoM" ;
  
  // ------------- BEGIN -----------

  Eff_local = Pur_local = 0.0 ;

  // ------------------------

  // compute efficiency
  NFAIL       = INPUTS.NFAIL[iTypeTrue][isep] ; 
  NTRUE_TOT   = NFAIL ;   // include training failures
  NTRUE_TYPE  = 0 ;
  for(iType = 0; iType < NTYPE; iType++ ) {    
    NTRAIN     = INPUTS.NTRAIN[iTypeTrue][iType][isep] ;
    NTRUE_TOT += NTRAIN ;
    if ( iTypeTrue == iType )  { NTRUE_TYPE += NTRAIN ; } 
  }
  if ( NTRUE_TOT > 0 )
    { Eff_local = (double)NTRUE_TYPE / (double)NTRUE_TOT ; }


  // ---------------------------------------
  // compute  purity  
  dNtrue = dNfalse = 0.0  ;
  for(iType = 0; iType < NTYPE; iType++ ) {    
    // For NTRAIN, note that iType loops over the true-type index,
    // and iTypeTrue --> train-type index.
    NTRAIN    = INPUTS.NTRAIN[iType][iTypeTrue][isep] ;  
    dNTRAIN = (double)NTRAIN ;

    if ( iTypeTrue == iType ) 
      { dNtrue  = dNTRAIN ; }  // true types identified correctly
    else  
      { dNfalse += dNTRAIN; }  // false NN types

    /*
    // xxxxxxxxxxxxxxx
    if ( isep == 1979 && iTypeTrue == 1 && Wfalse == 1.0 ) {
      printf(" iType=%d(%2d)  NTRUE=%5d(NFAIL=%3d)  NTRAIN=%5d  "
	     "Ntrue/Nfalse=%5d/%d \n",
	     iType,INPUTS.TrueType[iType],
	     NTRUE,NFAIL,  NTRAIN, (int)dNtrue, (int)dNfalse );
      fflush(stdout);
    }
    // xxxxxxxxx
    */

  }  // iType
  if ( dNtrue > 0.0 || dNfalse > 0.0 ) 
    { Pur_local = dNtrue / ( dNtrue + (Wfalse * dNfalse) ) ;  }



  // ------------------------------
  FoM = Eff_local * Pur_local ;

  // load output variables
  *Eff = Eff_local ;
  *Pur = Pur_local ;

  return(FoM) ;

} // end of getFoM


// ==============================================
void  dumpLine_UntrainedPurity(int iTypeTrue) {
  
  int ivar;
  float PUR = INPUTS.UntrainedPurity[iTypeTrue] ;
  //  char fnam[] = "dumpLine_UntrainedPurity" ;

  printf(" %4d  ---   %5.3f x %5.3f = %5.3f  <------ Untrained Purity  \n",
	 INPUTS.TrueType[iTypeTrue], 1.0, PUR, PUR );

  IROW_OUT++ ;
  fprintf(FP_OUT,"ROW:      %2d      %d       %d   ",
	  IROW_OUT, INPUTS.TrueType[iTypeTrue], 0 );
  fprintf(FP_OUT,"%5.3f %5.3f %5.3f  %5d  ",
	  1.0, PUR, PUR, 0 );

  for(ivar=0; ivar < INPUTS.NVAR; ivar++ ) 
    { fprintf(FP_OUT, "%7.4f  ", 0.0 );  }

  fprintf(FP_OUT,"\n");

} // end


// =======================================
void dump_NNINP_VARDEF(int iTypeTrue) {
  int ifalse ;

  printf(" &NNINP Implementation for snlc_fit.exe :\n" ) ;

  for(ifalse = 0; ifalse < N_WFALSE; ifalse ++ ) {
    printf("    Wfalse=%d :  %s \n", 
	   WFALSE_LIST[ifalse],
	   NNINP_VARDEF[iTypeTrue][ifalse] );    
  }
  fflush(stdout);
  return ;

}  // end of dump_NNINP_VARDEF


// ======================================
void dump_forPipeline(int iTypeTrue) {

  char fnam[] = "dump_forPipeline" ;
  char key[] = "SUMMARY_for_PIPELINE" ;
  int ifalse, IFALSE;

  // ---------

  // determine which ifalse index to use
  IFALSE = -1;
  for(ifalse=0 ; ifalse < N_WFALSE ; ifalse++ ) {
    if ( INPUTS.WFALSE_PIPELINE == WFALSE_LIST[ifalse] ) 
      { IFALSE = ifalse; }
  }

  if ( IFALSE < 0 ) {
    sprintf(msgerr1,"Invalid WFALSE_PIPELINE = %d", 
	    INPUTS.WFALSE_PIPELINE );
    sprintf(msgerr2,"Check WFALSE_LIST");
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }

  /*
  printf("\nSUMMARY_for_PIPELINE:  TrueType= %d  Pur_noNN= %.3f   "
	 "FOM=  %.3f x %.3f = %.3f   WARN = '%s' \n\n"
	 ,INPUTS.TrueType[iTypeTrue]
	 ,NNFOM[iTypeTrue].PURITY_noNN
	 ,NNFOM[iTypeTrue].EFF[ifalse]
	 ,NNFOM[iTypeTrue].PPURITY[ifalse]
	 ,NNFOM[iTypeTrue].FOM[ifalse]
	 ); */

  printf("\n");

  printf("%s[TrueType]:     %d \n", 
	 key, INPUTS.TrueType[iTypeTrue] );

  printf("%s[WFALSE]:  %d \n", 
	 key, INPUTS.WFALSE_PIPELINE );

  printf("%s[Purity_noNN]:  %.3f \n", 
	 key, NNFOM[iTypeTrue].PURITY_noNN );

  printf("%s[FOM]:          %.3f x %.3f = %.3f \n"
	 ,key
	 ,NNFOM[iTypeTrue].EFF[IFALSE]
	 ,NNFOM[iTypeTrue].PPURITY[IFALSE]
	 ,NNFOM[iTypeTrue].FOM[IFALSE] );
  printf("%s[WARNINGS]:     '%s' \n", 
	 key,NNFOM[iTypeTrue].WARN_BOUNDARY[IFALSE]);

  printf("\n");
  return ;

} // end of dump_forPipeline
