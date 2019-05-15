/***********************************************************

 Created Apr 2013 by R.Kessler
 Toolkit for Nearest Neighbor (NN) method.
 Can be used by any fitting program.
 See manual Sec 8 for usage with snana.

           Functions to call by external program

  EXTERNAL_INIT
  -> NVAR = SNANA_NEARNBR_RDINPUT()  ! read inputs from &NNINP namelist
                                     ! returns NVAR>0 if option is set

            -> NEARNBR_INIT()
            -> NEARNBR_SET_TRAINPATH(trainPath)  // optional path
            -> NEARNBR_SET_TRAINFILE(trainFile)  // called for each file
            -> NEARNBR_SET_TRUETYPE(VARNAME_TRUETYPE)
            -> NEARNBR_SET_SEPMAX(VARNAME,*SEPMAX) // call for each var
            -> NEARNBR_INIT2(ISPLIT)

  EXTERNAL_FIT (after fit)
            -> NEARNBR_LOADVAL(CCID,VARNAME,VALUE) // call for each var
            -> NEARNBR_APPLY(CCID)
            -> NEARNBR_GETRESULTS(CCID, ...)  // return results


                     HISTORY

  Oct 26 2014: switch to refactored table-read system so that all
               table formats are treated the same and no more special 
               code for reading ascii vs. hbook/root.

  Feb 2016: re-do NEARNBR_GETRESULTS( . .. ) to return NCELL for each
            type, instead of only returning FRAC_CELL for best ITYPE.
            Allows calling function to determine separate probabilities
            for Ibc and II.

  Jun 21 2016: refactor to compute NN_PTRAIN and NN_PROB, where
               the latter includes NN cut. No change in ITYPE_BEST.

  July 6, 2018: 
    + fix inputs to get_1DINDEX to account for previous refactor
      where indices go from 0 to N-1 instead of 1-N.

  Apr 11 2019: fill HID = HOFF+40 with total number of training events.

**********************************************/

#include <stdio.h> 
#include <string.h>
#include <stdlib.h>   
#include <unistd.h>
#include <math.h>  

#include "sntools.h" 

#include "sntools_nearnbr.h"
#include "sntools_output.h" 


// ============================================
void NEARNBR_INIT(void) {

  char fnam[] = "NEARNBR_INIT" ;

  sprintf(BANNER,"%s: Init Nearest Neighbor method", fnam );
  print_banner(BANNER);
  fflush(stdout);

  NEARNBR_INPUTS.TRAINFILE_PATH[0] = 0 ;
  NEARNBR_INPUTS.VARNAME_TRUETYPE[0] = 0 ;

  NEARNBR_INPUTS.NVAR        = 0 ;
  NEARNBR_INPUTS.NTRAINFILE  = 0 ;
  NEARNBR_TRAINLIB.NTRUETYPE = 0 ;
  NEARNBR_INPUTS.FILLHIST    = 0; 

  NEARNBR_INPUTS.TRAIN_ODDEVEN = 0;

  // hard-wire PROB cuts for now .. maybe later define via &NNINP.
  NEARNBR_INPUTS.CUTPROB     = 0.5 ; 
  NEARNBR_INPUTS.NSIGMA_PROB = 1.0 ;

  NN_TRAINFLAG = 0 ;
  NN_APPLYFLAG = 0 ;

  NEARNBR_CELLMAP.DOFLAG=0;
  NEARNBR_CELLMAP.NCHOP_PER_VAR = 0;

  nearnbr_reset();


}
void nearnbr_init__(void) { NEARNBR_INIT(); }


// ============================================
void NEARNBR_SET_TRAINPATH(char *path) {
  sprintf(NEARNBR_INPUTS.TRAINFILE_PATH,"%s", path);
  printf("\t Train-file path: %s\n", path);
  fflush(stdout);
}
void nearnbr_set_trainpath__(char *path) { NEARNBR_SET_TRAINPATH(path); }


// ============================================
void NEARNBR_SET_TRAINFILE(char *file) {
  int N ;
  N =  NEARNBR_INPUTS.NTRAINFILE ;
  sprintf(NEARNBR_INPUTS.TRAINFILE_LIST[N], "%s", file);
  NEARNBR_INPUTS.NTRAINFILE++ ;
  printf("\t Include train-file: %s\n", file);
  fflush(stdout);
}
void nearnbr_set_trainfile__(char *file) { NEARNBR_SET_TRAINFILE(file); }




// ============================================
void NEARNBR_SET_TRUETYPE(char *varName) {
  printf("\t Train-file variable with true type: %s\n", varName);
  sprintf(NEARNBR_INPUTS.VARNAME_TRUETYPE, "%s", varName);
} 
void nearnbr_set_truetype__(char *varName) 
{ NEARNBR_SET_TRUETYPE(varName); }



void NEARNBR_SET_ODDEVEN(void) 
{ NEARNBR_INPUTS.TRAIN_ODDEVEN = 1; }
void nearnbr_set_oddeven__(void) { NEARNBR_SET_ODDEVEN(); }



// ============================================
void NEARNBR_SET_SEPMAX(char *varName, double *SEPMAX) {
  int N, DOTRAIN ;

  N = NEARNBR_INPUTS.NVAR ;
  sprintf(NEARNBR_INPUTS.VARNAMES[N], "%s", varName) ;

  NEARNBR_INPUTS.SEPMAX_MIN[N] = SEPMAX[0] ; 
  NEARNBR_INPUTS.SEPMAX_MAX[N] = SEPMAX[1] ; 
  NEARNBR_INPUTS.SEPMAX_BIN[N] = SEPMAX[2] ; 
  DOTRAIN = (SEPMAX[2] > 0.0) ; // non-zero bin size 
  
  // set global TRAINFLAG if any bin-size is > 0
  if ( DOTRAIN ) 
    { NN_TRAINFLAG = 1; }
  else
    { NN_APPLYFLAG = 1; }

  NEARNBR_INPUTS.NVAR++ ;
}
void nearnbr_set_sepmax__(char *varName, double *SEPMAX) {
  NEARNBR_SET_SEPMAX(varName,SEPMAX);
}


// ==============================================
void NEARNBR_CELLMAP_INIT(int NCHOP_PER_VAR ) {

  // Created Jan 2017
  // Init cells for faster lookup of NN.
  // kdtree would be better, but this algorithm is simpler.
  // Basic idea is to define NVAR-dimensional cells on a
  // regular grid, and for each cell store all training events 
  // landing in the cell. Also include training events lying
  // with SEPMAX of any cell edge, to ensure that all NN
  // are included for any event inside the cell.
  //

  int  MEMI = sizeof(int);
  int  MEMF = sizeof(float);
  int  NCHOP_LOCAL, NVAR ;
  int  NTRAIN_TOT = NEARNBR_TRAINLIB.NTOT; 

  int  NCHOP_ARRAY[MXVAR_NEARNBR], ICELL_ARRAY[MXVAR_NEARNBR];
  int  DOFLAG = NEARNBR_CELLMAP.DOFLAG ;
  char fnam[] = "NEARNBR_CELLMAP_INIT" ;

  // ---------- BEGIN ----------

  if ( DOFLAG==0 && NCHOP_PER_VAR>1 ) {
    NEARNBR_CELLMAP.DOFLAG = 1 ;
    NEARNBR_CELLMAP.NCHOP_PER_VAR = NCHOP_PER_VAR ;
    return ;
  }

  if ( DOFLAG == 0 ) { return ; }

  NCHOP_LOCAL = NEARNBR_CELLMAP.NCHOP_PER_VAR ;
  sprintf(BANNER,"%s: Chop Cells for Faster Lookup ", fnam );
  print_banner(BANNER);

  NVAR = NEARNBR_INPUTS.NVAR ;
  NEARNBR_CELLMAP.NVAR = NVAR ;

  // - - - - - - - - - - - -
  int itrain, ivar, ICELL_1D, icell, NLIST, NCELL_TOT=1 ;
  float VAL, VAL_MIN, VAL_MAX, VAL_BIN, RSQ, DIF, SQSEPMAX ;
  float VAL_ARRAY[MXVAR_NEARNBR], SEPMAX[MXVAR_NEARNBR];
  float VAL_CENTMP[MXVAR_NEARNBR], VAL_CENREF[MXVAR_NEARNBR];
  char *VARNAME ;

  for ( ivar=0; ivar < NEARNBR_INPUTS.NVAR ; ivar++ ) {
    NEARNBR_CELLMAP.VAL_MIN[ivar] = +1.0E8;
    NEARNBR_CELLMAP.VAL_MAX[ivar] = -1.0E8;
    NCHOP_ARRAY[ivar] = NCHOP_LOCAL ;
    NCELL_TOT *= NCHOP_LOCAL ;
  }
  NEARNBR_CELLMAP.NCELL_TOT = NCELL_TOT ;

  // ----- malloc arrays to hold training list per cell -----
  NEARNBR_CELLMAP.NLIST = 
    (int*)malloc( NCELL_TOT * MEMI );
  NEARNBR_CELLMAP.ITRAIN_LIST = 
    (int**)malloc( NCELL_TOT * sizeof(int*) );

  NEARNBR_CELLMAP.MEMTOT = 0 ;
  NEARNBR_CELLMAP.MEMTOT = 0 ;
  for(icell=0; icell < NCELL_TOT; icell++ ) {  
    NEARNBR_CELLMAP.NLIST[icell] = 0;     
    NEARNBR_CELLMAP.ITRAIN_LIST[icell] = 
      (int*) malloc(MEMI*BUFFSIZE_CELLMAP_NEARNBR); 
    NEARNBR_CELLMAP.MEMTOT += (MEMI*BUFFSIZE_CELLMAP_NEARNBR);
  }


  // get min & max value for each training variable
  for(itrain=0; itrain < NTRAIN_TOT ; itrain++ ) { 
    for ( ivar=0; ivar < NEARNBR_INPUTS.NVAR ; ivar++ ) {
      VAL = NEARNBR_TRAINLIB.FITRES_VALUES[ivar][itrain] ;
      if ( VAL < NEARNBR_CELLMAP.VAL_MIN[ivar] )
	{ NEARNBR_CELLMAP.VAL_MIN[ivar] = VAL; }
      if ( VAL > NEARNBR_CELLMAP.VAL_MAX[ivar] )
	{ NEARNBR_CELLMAP.VAL_MAX[ivar] = VAL; }
    }
  }

  printf("\n  Cell-Range of each training variable: \n");
  for ( ivar=0; ivar < NEARNBR_INPUTS.NVAR ; ivar++ ) {    
    VARNAME = NEARNBR_INPUTS.VARNAMES[ivar] ;
    VAL_MIN = NEARNBR_CELLMAP.VAL_MIN[ivar] ;
    VAL_MAX = NEARNBR_CELLMAP.VAL_MAX[ivar] ;
    VAL_BIN = (VAL_MAX-VAL_MIN)/ (float)NCHOP_LOCAL ;
    NEARNBR_CELLMAP.VAL_BIN[ivar] = VAL_BIN ;
    SEPMAX[ivar] = sqrtf(NEARNBR_LIST_SQSEPMAX[ivar][0]);

    printf("   %4s : %8.4f to %8.4f  (cell size=%5.3f, SEPMAX=%.3f) \n",
	   VARNAME, VAL_MIN, VAL_MAX, VAL_BIN, SEPMAX[ivar]);
  }
  printf("  Total number of lookup cells: %d^%d = %d \n",
	 NCHOP_LOCAL, NVAR, NCELL_TOT);

  fflush(stdout);

  
  // init map for multi-dim indices -> 1D index
  init_1DINDEX(ID1D_CELLMAP_NEARNBR, NVAR, NCHOP_ARRAY);



  // loop over training and assign each cell
  for(itrain=0; itrain < NTRAIN_TOT ; itrain++ ) { 

    for ( ivar=0; ivar < NVAR ; ivar++ ) 
      { VAL_ARRAY[ivar] = NEARNBR_TRAINLIB.FITRES_VALUES[ivar][itrain] ;  }

    getInfo_CELLMAP(+1, VAL_ARRAY,  &ICELL_1D );   // return ICELL_1D
    getInfo_CELLMAP(-1, VAL_CENREF, &ICELL_1D );   // return VAL_CENREF

    // check distance to all other cells
    for(icell=0; icell < NCELL_TOT; icell++ ) {

      // pass icell and return binCenter VAL_CENTMP
      getInfo_CELLMAP(-1, VAL_CENTMP, &icell ); 

      RSQ = 0.0 ;
      for ( ivar=0; ivar < NVAR ; ivar++ ) {
	if ( icell == ICELL_1D ) { continue ; }
	VAL_BIN  = NEARNBR_CELLMAP.VAL_BIN[ivar] ;
	DIF      = fabsf( VAL_CENREF[ivar] - VAL_CENTMP[ivar] ) ;
	if ( DIF >= 0.99*VAL_BIN ) 
	  { DIF -= VAL_BIN ; } // compare distance between closest edges
	if ( fabsf(DIF) > SEPMAX[ivar] ) 
	  { ivar=NVAR; RSQ=999.; continue; }
	SQSEPMAX = NEARNBR_LIST_SQSEPMAX[ivar][0] ;
	RSQ     += ( (DIF*DIF)/SQSEPMAX ) ;
      }

      if ( RSQ < 1.0 ) {
	realloc_NEARNBR_CELLMAP(icell); 
	NLIST = NEARNBR_CELLMAP.NLIST[icell] ;	
	NEARNBR_CELLMAP.ITRAIN_LIST[icell][NLIST] = itrain ;
	NEARNBR_CELLMAP.NLIST[icell]++ ;
      }

    }  // end icell loop

  } // itrain


  int  NCELL_USE=0, NLIST_MAX=0 ;
  for(icell=0; icell < NCELL_TOT; icell++ ) {
    NLIST = NEARNBR_CELLMAP.NLIST[icell] ;
    if ( NLIST == 0 ) { continue ; }
    if ( NLIST > NLIST_MAX ) { NLIST_MAX = NLIST; }
    NCELL_USE++ ;
    //    printf("\t  icell = %3d --> NSTORE = %d \n", icell, NLIST );
  }

  printf("  Total number of used cells:  %d of %d. \n", 
	 NCELL_USE, NCELL_TOT);
  
  float FMAX = (float)NLIST_MAX / (float)NTRAIN_TOT ;
  printf("  Max list size for 1 cell: %d (%.3f of total train sample) \n",
	 NLIST_MAX, FMAX);

  float XMEM = (float)NEARNBR_CELLMAP.MEMTOT/1.0E6 ;
  printf("  Total memory storage for cells: %.3f Mb \n", XMEM);

  fflush(stdout);

  return ;

} // end NEARNBR_CELLMAP_INIT

// ========================================
void getInfo_CELLMAP(int OPT, float *VAL_ARRAY, int *ICELL_1D ) {

  // Inputs:
  //  OPT > 0 --> input VAL_ARRAY and output ICELL_1D
  //  OPT < 0 --> input ICELL_1D and output VAL_ARRAY at bin center

  int NVAR =   NEARNBR_CELLMAP.NVAR ;
  int ivar, icell, ICELL_ARRAY[MXVAR_NEARNBR] ;
  float VAL_MIN, VAL_MAX, VAL_BIN, VAL, xi ;
  char fnam[] = "getInfo_CELLMAP" ;

  // --------------- BEGIN --------------

  if ( OPT > 0 ) {
    // VAL_ARRAY is input; return 1D index ICELL_1D
    for(ivar=0; ivar < NVAR; ivar++ ) {
      VAL_MIN = NEARNBR_CELLMAP.VAL_MIN[ivar] ;
      VAL_MAX = NEARNBR_CELLMAP.VAL_MAX[ivar] ;
      VAL_BIN = NEARNBR_CELLMAP.VAL_BIN[ivar] ;
      VAL     = VAL_ARRAY[ivar] ;
      
      if ( VAL <= VAL_MIN ) { VAL = VAL_MIN + (VAL_MAX-VAL_MIN)*1.0E-6; }
      if ( VAL >= VAL_MAX ) { VAL = VAL_MAX - (VAL_MAX-VAL_MIN)*1.0E-6; }
      icell  = (int)( (VAL-VAL_MIN)/VAL_BIN ) ;
      ICELL_ARRAY[ivar] = icell  ; 
    } // ivar

    *ICELL_1D = get_1DINDEX(ID1D_CELLMAP_NEARNBR, NVAR, ICELL_ARRAY);
    
  }
  else {
    // 1D index ICELL_1D is input, return bin center VAL_ARRAY
    int NCHOP = NEARNBR_CELLMAP.NCHOP_PER_VAR ;
    int IDIGIT[MXVAR_NEARNBR];
    int J1, J2 ;
    for(ivar=0; ivar < NVAR; ivar++ ) { IDIGIT[ivar] = 0 ; }
    J1 = *ICELL_1D ;  J2=0; ivar=NVAR;
    while ( J1 != 0 ) {
      J2 = (J1%NCHOP);	J1 = (J1/NCHOP);  ivar-- ;
      IDIGIT[NVAR-ivar-1] = (float)J2 ; 
    }
    for(ivar=0; ivar < NVAR; ivar++ ) { 
      VAL_MIN = NEARNBR_CELLMAP.VAL_MIN[ivar] ;
      VAL_BIN = NEARNBR_CELLMAP.VAL_BIN[ivar] ;
      xi      = (float)IDIGIT[ivar] ;
      VAL_ARRAY[ivar] = VAL_MIN + VAL_BIN * (xi+0.5) ;
    }
    /*
    printf(" xxx IDIGIT = %d %d %d \n",
    IDIGIT[0], IDIGIT[1], IDIGIT[2] ); */
  }

  return ;
  

} // end getInfo_CELLMAP


// ============================
void realloc_NEARNBR_CELLMAP(int ICELL_1D) {

  float XNLIST, XNADD;
  int  MEM_REALLOC, NLIST;
  int  MEMI = sizeof(int);
  char fnam[] = "realloc_NEARNBR_CELLMAP";

  // ------------ BEGIN ------------

    NLIST  = NEARNBR_CELLMAP.NLIST[ICELL_1D] ;  
    XNLIST = (float)NLIST;  
    XNADD  = (float)BUFFSIZE_CELLMAP_NEARNBR ;
    
    // check if more memory is needed    
    if ( fmodf(XNLIST,XNADD) == 0.0 && NLIST>3 ) {
      // printf("\t Add %d to ITRAIN_LIST[%3d] = %d\n",NBUF_ADD,ICELL_1D,NTMP);
      MEM_REALLOC = MEMI * ( NLIST + BUFFSIZE_CELLMAP_NEARNBR);

      NEARNBR_CELLMAP.ITRAIN_LIST[ICELL_1D] = 
	(int*)realloc(NEARNBR_CELLMAP.ITRAIN_LIST[ICELL_1D],MEM_REALLOC ); 
      NEARNBR_CELLMAP.MEMTOT += (MEMI*BUFFSIZE_CELLMAP_NEARNBR);
    }

  return ;

} // end realloc_NEARNBR_CELLMAP

// ==============================================
void NEARNBR_INIT2(int ISPLIT) {

  // 2nd init after all input has been passed to NEARNBR_INPUTS
  // If ISPLIT=1 then book all plots.
  // If ISPLIT>1 then book only the 2D analysis plot

  int  ifile, ivar ;
  char fnam[] = "NEARNBR_INIT2" ;

  // ------------ BEGIN ---------------

  // abort if required input is missing
  nearnbr_check_inputs();

  // read each file and malloc/store array for each variable.
  NEARNBR_TRAINLIB.NTOT     = 0 ;
  NEARNBR_TRAINLIB.NTOT_USE = 0 ;

  for ( ifile=0; ifile < NEARNBR_INPUTS.NTRAINFILE; ifile++ ) 
    { nearnbr_read_trainLib(ifile); }

  float FMEM = 1.0E-6 * (float)NEARNBR_TRAINLIB.MEMTOT ;
  printf("\t TOTAL TRAINING SET MEMORY: %6.3f MB for %d EVENTS\n", 
	 FMEM, NEARNBR_TRAINLIB.NTOT );
  fflush(stdout) ;

  // alloate memory to store SQSEPMAX for each trainlib entry
  int MEMF = sizeof(float) * NEARNBR_TRAINLIB.NTOT ;
  for ( ivar=0; ivar < NEARNBR_INPUTS.NVAR ; ivar++ ) 
    { NEARNBR_STORE.SQSEP[ivar] = (float*)malloc(MEMF); }


  // create linear array of variable sets ; set NN_TRAINFLAG
  nearnbr_init_SEPMAX();

  // examine TRUETYPE array and store info for each possible TRUETYPE
  nearnbr_storeTrueTypes();

  // init SUBSET to be entire training lib
  nearnbr_init_SUBSET() ;

  // init speedup for APPLY mode
  if ( NN_APPLYFLAG ) { NEARNBR_CELLMAP_INIT(0) ; }

  // create histograms for training mode (i.e., multiple SEPMAX bins),
  if( NN_TRAINFLAG ) { nearnbr_makeHist(ISPLIT); }

  char mode[12] = "UNKNONW" ;
  if ( NN_TRAINFLAG )  { sprintf(mode,"TRAINING"); }
  if ( NN_APPLYFLAG )  { sprintf(mode,"ANALYSIS"); }
  printf("\n  Finished Nearest Neighbor (NN) init --> %s mode. \n\n",
	 mode); fflush(stdout);

}  // end of NEARNBR_INIT2

void nearnbr_init2__(int *ISPLIT) { NEARNBR_INIT2(*ISPLIT); }


// ====================================================
void nearnbr_check_inputs(void) {

  char *cptr ;
  char fnam[] = "nearnbr_check_inputs" ;

  // ------------- BEGIN -----------------

  if ( NEARNBR_INPUTS.NVAR <= 0 ) {
    sprintf(c1err,"No NN variables specified ??") ;
    sprintf(c2err,"Check &SNLCINP nml variable NEARNBR_SEPMAX_LIST");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  if ( NEARNBR_INPUTS.NTRAINFILE <= 0 ) {
    sprintf(c1err,"No NN train files specified ??") ;
    sprintf(c2err,"Check &SNLCINP nml variable NEARNBR_TRAINFILE_LIST");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }



  // require true-type  variable for training option
  cptr = NEARNBR_INPUTS.VARNAME_TRUETYPE ;
  if ( strlen(cptr) == 0 ) {
    sprintf(c1err,"Did not specify variable name for true type.");
    sprintf(c2err,"Check &SNLCINP nml variable NEARNBR_TRUETYPE_VARNAME");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


} // nearnbr_check_inputs 


// ======================================
void nearnbr_read_trainLib(int ifile) {

  // Jun 14 2016: read only zHD (z no longer valid)
  // Apr 06 2019: read CID to use with ODDEVEN option.

  int OPT_ABORT_noVAR = 2;
  int NROW, irow, ivar, NVAR, NVAR_TOT, MEMF, MEMI, LFIRST, IFILETYPE ;
  int   *IPTR ;
  float *FPTR ;
  char  INFILE_FULL[MXCHAR_FILENAME] ;
  char  VARARG[80], TABLENAME[40], *ptrVar;
  char  fnam[] = "nearnbr_read_trainLib" ;

  // --------------- BEGIN ---------------

  // returns INFILE_FULL = full name of trainin file
  nearnbr_TRAIN_FILENAME(ifile,INFILE_FULL); 

  // ---------------------------------
  // open file based on its type, and read number of rows.

  NROW = 0 ;    IFILETYPE = -9 ;

  IFILETYPE = TABLEFILE_OPEN(INFILE_FULL, "read "); // q-> quiet mode 
  sprintf(TABLENAME, "%s", STRING_IDTABLE_FITRES[IFILETYPE] );

  NROW     = SNTABLE_NEVT(INFILE_FULL, TABLENAME); 
  NVAR_TOT = SNTABLE_READPREP(IFILETYPE,TABLENAME);

  // leave table file open
  
  if ( NROW == 0 ) {
    printf("\n PRE-ABORT DUMP: \n INFILE_FULL = '%s' \n", INFILE_FULL);	   
    sprintf(c1err,"Could not find training table");
    sprintf(c2err,"Bad table or non-existing input file (above)");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // -----------------------------------------
  // allocate memory based on NROW
  NEARNBR_TRAINLIB.MEMTOT = 0 ;
  LFIRST = (NEARNBR_TRAINLIB.NTOT == 0) ;
  MEMF   = (NEARNBR_TRAINLIB.NTOT + NROW) * sizeof(float) ; 
  MEMI   = (NEARNBR_TRAINLIB.NTOT + NROW) * sizeof(int) ; 

  // ------------------------------------------
  // Apr 5 2019: allocate CID and TRUETYPE
  if ( LFIRST )  { 
    NEARNBR_TRAINLIB.CID_VALUES = (int*)malloc(MEMI); 
    NEARNBR_TRAINLIB.TRUETYPE   = (int*)malloc(MEMI); 
  }
  else {
    NEARNBR_TRAINLIB.CID_VALUES = 
      (int*)realloc(NEARNBR_TRAINLIB.CID_VALUES, MEMI); 
    NEARNBR_TRAINLIB.TRUETYPE = 
      (int*)realloc(NEARNBR_TRAINLIB.TRUETYPE, MEMI); 
  }
  NEARNBR_TRAINLIB.MEMTOT += (2*MEMI) ;

  
  sprintf(VARARG, "CID:I" );
  IPTR   = &NEARNBR_TRAINLIB.CID_VALUES[NEARNBR_TRAINLIB.NTOT] ;
  SNTABLE_READPREP_VARDEF( VARARG, IPTR, NROW, OPT_ABORT_noVAR);

  // loop over training-set variables
  // Last element is integer TRUETYPE that is treated differently.

  NVAR = NEARNBR_INPUTS.NVAR;
  for(ivar=0; ivar <= NVAR; ivar++ ) {

    if ( ivar < NVAR )
      { ptrVar = NEARNBR_INPUTS.VARNAMES[ivar] ; }   
    else  
      { ptrVar = NEARNBR_INPUTS.VARNAME_TRUETYPE ;  }

    if ( LFIRST )  { 
      NEARNBR_TRAINLIB.FITRES_VALUES[ivar] = 
	(float*)malloc(MEMF); 
    }
    else {
      NEARNBR_TRAINLIB.FITRES_VALUES[ivar] = 
	(float*)realloc(NEARNBR_TRAINLIB.FITRES_VALUES[ivar], MEMF); 
    }
    NEARNBR_TRAINLIB.MEMTOT += MEMF ;

    FPTR = &NEARNBR_TRAINLIB.FITRES_VALUES[ivar][NEARNBR_TRAINLIB.NTOT] ;

    if ( strcmp(ptrVar,"z") == 0 ) {
      sprintf(VARARG,"zHD:F") ;
    }
    else {
      sprintf(VARARG,"%s/F", ptrVar);
    }

    SNTABLE_READPREP_VARDEF( VARARG, FPTR, NROW, OPT_ABORT_noVAR);
    
  } // end of ivar

 

  // -------------------------------
  // read table file and increment total number of entries 

  NROW  = SNTABLE_READ_EXEC();
  // xxx obolete  TABLEFILE_CLOSE(INFILE_FULL) ; // close TABLE file
  CDTOPDIR_OUTPUT();


  // Apr 6 2019: check option to use only ODD CID .xyz
  int NROW_ODD;
  if ( NEARNBR_INPUTS.TRAIN_ODDEVEN )  
    { NROW_ODD = nearnbr_storeODD_trainLIB(NROW);  NROW=NROW_ODD; }


  NEARNBR_TRAINLIB.NTOT  +=  NROW;

  // -----------------------------------
  // transfer TRUETYPE from float to integer array;
  // make sure first to allocate integer array.

  for(irow = 0; irow < NEARNBR_TRAINLIB.NTOT; irow++ ) {
    NEARNBR_TRAINLIB.TRUETYPE[irow] = 
      (int)NEARNBR_TRAINLIB.FITRES_VALUES[NVAR][irow] ;
  }

  return ;

} // end of nearnbr_read_trainLib


// =============================================
int nearnbr_storeODD_trainLIB(int NROWTOT_NEW) {

  // Apr 2019
  // Store only ODD CID for trainLib read from FITRES file.

  int NVAR = NEARNBR_INPUTS.NVAR;
  int NROWSTORE_ODD = 0 ;
  int CID, IROW, irow, ivar ;
  double VALUE;
  char fnam[] = "nearnbr_storeODD_trainLIB";

  // ------------ BEGIN ---------------
  
  for(irow = 0; irow < NROWTOT_NEW; irow++ ) {
    IROW = NEARNBR_TRAINLIB.NTOT + irow;
    CID  = NEARNBR_TRAINLIB.CID_VALUES[IROW];
    if ( CID%2  == 1 ) {  // select ODD            
      for(ivar=0; ivar <= NVAR; ivar++ ) {
	VALUE = NEARNBR_TRAINLIB.FITRES_VALUES[ivar][IROW];
	NEARNBR_TRAINLIB.FITRES_VALUES[ivar][NROWSTORE_ODD] = VALUE;	 
      }
      NROWSTORE_ODD++;
    } // end ODD if-block
  } // end irow loop

  printf(" %s: store %d odd CID rows from %d total rows\n",
	 fnam, NROWSTORE_ODD, NROWTOT_NEW);
  fflush(stdout);

  return(NROWSTORE_ODD);

} // end nearnbr_storeODD_trainLIB

// =============================================
void nearnbr_TRAIN_FILENAME(int ifile, char *TRAIN_FILENAME) {

  char *INFILE, *PATH;
  char fnam[] = "nearnbr_TRAIN_FILENAME" ;

  // --------- BEGIN ----------

  INFILE = NEARNBR_INPUTS.TRAINFILE_LIST[ifile] ;
  PATH   = NEARNBR_INPUTS.TRAINFILE_PATH ;

  if ( strlen(PATH) > 0 ) 
    {  sprintf(TRAIN_FILENAME, "%s/%s", PATH, INFILE ); }
  else
    {  sprintf(TRAIN_FILENAME, "%s", INFILE ); }

  return ;

} // end nearnbr_TRAIN_FILENAME

// =============================================
void nearnbr_apply_trainLib(void) {

  // Created June 21 2016 by R.Kessler
  //
  // determine training prob (P_TRAIN) for each event in the 
  // training sample. Needed to compute P_BAYES, including
  // the NN requirement, otherwise PTRAIN(Ia) is too high.
  // The reason is that within each NN cell, the NN cut
  // removes some of the CC, and thus P_BAYES > P_TRAIN.
  // P_BAYES is not used for NN selection, but is used 
  // in fitting programs requiring a Bayesian type-probability

  int  NROW      = NEARNBR_TRAINLIB.NTOT ;
  int  NVAR      = NEARNBR_INPUTS.NVAR;
  int  NTRUETYPE = NEARNBR_TRAINLIB.NTRUETYPE ;
  int  MEMF      = NROW * sizeof(float);
  char fnam[]    = "nearnbr_apply_trainLib" ;
  
  int  irow, ivar, itype, ITYPE_BEST, NTYPE;
  int  ITYPE_LIST[NTRUETYPE_MAX];
  int  NCELL_TRAIN_LIST[NTRUETYPE_MAX] ;
  double d_val, P_TRAIN, XNTOT ;
  char *varNameList[MXVAR_NEARNBR], CCID[12] ;

  // ------------ BEGIN --------------

  NN_APPLYFLAG = 2 ;  // set to "internal apply mode"

  printf("\n  %s: compute TRAIN-PROB for each event in training sample.\n", 
	 fnam);
  printf("\t\t (needed to compute proper P_BAYES) \n");
  fflush(stdout);

  // allocate memory to store P_TRAIN for each possible type.
  NEARNBR_TRAINLIB.P_TRAIN = (float**)malloc(NTRUETYPE * sizeof(float*) ); 
  for(itype=0; itype < NTRUETYPE; itype++ ) 
    { NEARNBR_TRAINLIB.P_TRAIN[itype] =  (float*)malloc(MEMF);   }

  // store local list of varNames before trainLib loop
  for(ivar = 0; ivar <= NVAR; ivar++ ) {
    if ( ivar < NVAR )
      { varNameList[ivar] = NEARNBR_INPUTS.VARNAMES[ivar] ; }   
    else  
      { varNameList[ivar] = NEARNBR_INPUTS.VARNAME_TRUETYPE ;  }
  }

  // ----------------------------------
  for(irow=0; irow < NROW; irow++ ) {      

    // didn't read trainLib CID, so just use row number
    sprintf(CCID,"TRAIN%8.8d", irow);
    
    // load values used for NN (e.g, z, x1, c )
    for(ivar = 0; ivar <= NVAR; ivar++ )  { 
      d_val = (double)NEARNBR_TRAINLIB.FITRES_VALUES[ivar][irow] ;  
      NEARNBR_LOADVAL(CCID, varNameList[ivar], d_val ) ;
    }

    NEARNBR_APPLY(CCID) ;

    NEARNBR_GETRESULTS(CCID, &ITYPE_BEST, &NTYPE,
		       ITYPE_LIST, NCELL_TRAIN_LIST ) ;

    XNTOT = 0.0 ;
    for(itype=0; itype < NTYPE; itype++ ) 
      { XNTOT += (double)NCELL_TRAIN_LIST[itype] ;  }

    for(itype=0; itype < NTYPE; itype++ ) {
      NEARNBR_TRAINLIB.P_TRAIN[itype][irow] = 0.0 ;
      if ( ITYPE_LIST[itype] >= 0 )  {
	P_TRAIN = (double)NCELL_TRAIN_LIST[itype]/XNTOT ;
	NEARNBR_TRAINLIB.P_TRAIN[itype][irow] = (float)P_TRAIN ;
      }
    }
    
  } // end irow loop over trainLib rows


  NN_APPLYFLAG = 1 ;  // back to external apply mode

  return ;

} // end nearnbr_apply_trainLib


// =============================================
void nearnbr_init_SUBSET(void) {

  // default subset for analysis is entire trainlib.

  int itrain ;
  int NTRAIN = NEARNBR_TRAINLIB.NTOT; 

  NEARNBR_TRAINLIB.ITRAIN = (int  *)malloc(sizeof(int) * NTRAIN); 
  NEARNBR_TRAINLIB.NSUBSET = NTRAIN ;
  for(itrain=0; itrain < NTRAIN ; itrain++ ) 
    { NEARNBR_TRAINLIB.ITRAIN[itrain] = itrain ;  }

} // end of nearnbr_init_SUBSET

// ==================================
void  nearnbr_init_SEPMAX(void) {

  int  NVAR, IVAR, NBIN[MXVAR_NEARNBR],  NTOT ;
  float SEPMAX_MIN, SEPMAX_MAX, SEPMAX_BIN, SEPMAX, tmpDif ;
  char *VARNAME ;
  char fnam[] = "nearnbr_init_SEPMAX" ;

  // ---------- BEGIN --------------

  // determine how many bins per dimension from input
  NTOT = 1;
  NVAR = NEARNBR_INPUTS.NVAR ;

  for(IVAR=0; IVAR < NVAR ; IVAR++ ) {  
    SEPMAX_MIN = NEARNBR_INPUTS.SEPMAX_MIN[IVAR] ;
    SEPMAX_MAX = NEARNBR_INPUTS.SEPMAX_MAX[IVAR] ;
    SEPMAX_BIN = NEARNBR_INPUTS.SEPMAX_BIN[IVAR] ;
    VARNAME    = NEARNBR_INPUTS.VARNAMES[IVAR] ;
    
    if ( SEPMAX_MAX == 0.0 ) {
      SEPMAX_MAX = SEPMAX_MIN ;
      NEARNBR_INPUTS.SEPMAX_MAX[IVAR] = SEPMAX_MAX ;
    }

    if ( SEPMAX_MIN <= 1.0E-9 ) {
      sprintf(c1err,"Invalid SEPMAX_MIN = %f", SEPMAX_MIN);
      sprintf(c2err,"for NEARNBR var = '%s'  (range=%f to %f)", 
	      VARNAME, SEPMAX_MIN, SEPMAX_MAX );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }


    if ( SEPMAX_MIN == SEPMAX_MAX ) 
      { NBIN[IVAR] = 1 ; }
    else { 
      if ( SEPMAX_BIN <= 1.0E-9 ) {
	sprintf(c1err,"Invalid SEPMAX_BIN = %f", SEPMAX_BIN);
	sprintf(c2err,"for NEARNBR var = '%s'  (range=%f to %f)", 
		VARNAME, SEPMAX_MIN, SEPMAX_MAX );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      tmpDif = (SEPMAX_MAX - SEPMAX_MIN) + 0.00001 ;
      NBIN[IVAR] = ceil(tmpDif/SEPMAX_BIN) ;
    }

    NTOT *= NBIN[IVAR];

    printf("\t NBIN_SEPMAX(%-12s) = %4d  (%6.3f to %6.3f / %6.3f)\n", 
	   VARNAME, NBIN[IVAR], SEPMAX_MIN, SEPMAX_MAX , SEPMAX_BIN);
    
  }  // end of IVAR


  NBINTOT_SEPMAX_NEARNBR = NTOT ;
  if ( NTOT > MXBIN_SEPMAX_NEARNBR ) {
    sprintf(c1err,"NBINTOT_SEPMAX_NEARNBR = %d exceeds bound.", NTOT );
    sprintf(c2err,"MXBIN_SEPMAX_NEARNBR = %d", MXBIN_SEPMAX_NEARNBR ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // ------------------
  // allocate memory
  int MEMTOT, MEMF, MEMI, NTYPE ;

  MEMF = sizeof(float)*NTOT;
  MEMI = sizeof(int)  *NTOT;
  MEMTOT = 0;

  NTYPE =   NEARNBR_TRAINLIB.NTRUETYPE ;


  NEARNBR_LIST_SQSEPMAX = (float**)malloc(sizeof(float*)*NVAR ); 
  NEARNBR_LIST_SEPMAX   = (float**)malloc(sizeof(float*)*NVAR ); 
  NEARNBR_LIST_NTYPE    = (int  **)malloc(sizeof(int  *)*NVAR ); 
  for(IVAR=0; IVAR < NVAR ; IVAR++ ) {  
    NEARNBR_LIST_SQSEPMAX[IVAR] = (float*)malloc(MEMF); 
    NEARNBR_LIST_SEPMAX[IVAR]   = (float*)malloc(MEMF); 
    NEARNBR_LIST_NTYPE[IVAR]    = (int  *)malloc(MEMI); 
    MEMTOT += (MEMF+MEMI) ;
  }


  if ( NTOT > 1 )  { 
    char   MEM_UNIT[4];
    double XMEM ;
    if ( MEMTOT < 1000000 ) 
      { XMEM = (double)MEMTOT * 1.0E-3 ;   sprintf(MEM_UNIT,"kB");  }
    else 
      { XMEM = (double)MEMTOT * 1.0E-6 ;   sprintf(MEM_UNIT,"MB");  }

    printf("\t Allocate %.2f %s memory for %d NEARNBR-SEPMAX bins. \n", 
	   XMEM, MEM_UNIT, NTOT); 
    fflush(stdout);
  }


  // define each SEPMAX variable for each NVAR-dim IBIN
  int IBTOT, NBTMP, IBIN_LIST[MXVAR_NEARNBR];
  int LDMP = 0 ;

  for(IBTOT = 0; IBTOT < NTOT; IBTOT++ ) {

    if ( LDMP ) { printf(" xxx IBTOT=%5d: SEPMAX= ", IBTOT); }
    NBTMP = 1; 
    for(IVAR=0; IVAR < NVAR ; IVAR++ ) {  

      IBIN_LIST[IVAR] = (IBTOT/NBTMP) % NBIN[IVAR] ;
      SEPMAX_MIN = NEARNBR_INPUTS.SEPMAX_MIN[IVAR] ;
      SEPMAX_BIN = NEARNBR_INPUTS.SEPMAX_BIN[IVAR] ;  
      SEPMAX     = SEPMAX_MIN + SEPMAX_BIN*(float)IBIN_LIST[IVAR] ;
    
      NEARNBR_LIST_SQSEPMAX[IVAR][IBTOT] = SEPMAX*SEPMAX ;
      NEARNBR_LIST_SEPMAX[IVAR][IBTOT]   = SEPMAX ;


      NBTMP *= NBIN[IVAR];
      if ( LDMP ) { printf(" %6.3f(%3d) ", SEPMAX, IBIN_LIST[IVAR] ); }
    }

    if ( LDMP ) { printf("\n"); } 

  }  // NTOT loop

  fflush(stdout);

} // end of nearnbr_init_SEPMAX



// ===================================
void NEARNBR_LOADVAL(char *ccid, char *varName, double d_val ) {

  // store value d_val for variable varName
  // Oct 7 2016: padd ccid for error message.

  int  IVAR, ivar, NVAR, NTYPE, ITYPE, ITYPE_TMP, isp, ISP ;
  char *tmpVar ;
  char fnam[] = "NEARNBR_LOADVAL" ;

  // ---------------- BEGIN ----------------

  // first check for true type
  if ( strcmp(varName, NEARNBR_INPUTS.VARNAME_TRUETYPE) == 0 ) {
    ITYPE = (int)d_val ;
    NEARNBR_STORE.TRUETYPE_LOAD = ITYPE ;


    // get sparse ITYPE
    ISP = -9 ;
    NTYPE = NEARNBR_TRAINLIB.NTRUETYPE ;
    for ( isp=0 ; isp < NTYPE ; isp++ ) {
      ITYPE_TMP = NEARNBR_TRAINLIB.TRUETYPE_LIST[isp] ;
      if ( ITYPE == ITYPE_TMP ) { ISP = isp ; }
    }

    if ( ISP < 0 ) {
      printf("\n PRE-ABORT DUMP: \n");
      for ( isp=0 ; isp < NTYPE ; isp++ ) {
	ITYPE_TMP = NEARNBR_TRAINLIB.TRUETYPE_LIST[isp] ;
	printf("\t Available %s = %d \n", varName, ITYPE_TMP); 
      }

      sprintf(c1err,"Could not find sparse TYPE/index");
      sprintf(c2err,"VARNAME='%s'  ITYPE=%d  CID=%s", 
	      varName,   ITYPE, ccid);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    NEARNBR_STORE.TRUETYPE_SPARSE = ISP ;
    return ;
  }

  
  // --------------------------------------
  NEARNBR_STORE.NVAL_LOAD++ ;

  // find IVAR that matches varName
  NVAR = NEARNBR_INPUTS.NVAR ; 
  IVAR = -9 ;
  for ( ivar=0; ivar < NVAR ; ivar++ ) {
    tmpVar = NEARNBR_INPUTS.VARNAMES[ivar] ;
    if ( strcmp(varName,tmpVar) == 0 ) 
      { IVAR = ivar ; goto FOUND_IVAR ; }
  }

  if ( IVAR < 0 ) {

    printf("\n PRE-ABORT DUMP: \n");
    for ( ivar=0; ivar < NVAR ; ivar++ ) {
      printf("\t Valid NN varName: '%s' \n", 
	     NEARNBR_INPUTS.VARNAMES[ivar]);
    }

    sprintf(c1err,"Could not find NN variable '%s' = %f (CID=%s)", 
	    varName, d_val, ccid );
    sprintf(c2err,"Check NN variables.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

 FOUND_IVAR:
  NEARNBR_STORE.VALUE_LOAD[IVAR] = (float)d_val ;
  NEARNBR_STORE.NLOAD[IVAR]++ ;

 
  return ;

} // end of NEARNBR_LOADVAL

void nearnbr_loadval__(char *ccid, char *varName, double *d_val ) {
  NEARNBR_LOADVAL(ccid, varName, *d_val ) ;
}


// ========================================
void nearnbr_reset(void) {

  int IVAR ;

  for ( IVAR=0; IVAR < MXVAR_NEARNBR; IVAR++ ) {  
    NEARNBR_STORE.VALUE_LOAD[IVAR] = -9.0 ; 
    NEARNBR_STORE.NLOAD[IVAR] =  0   ; 
  }

  NEARNBR_STORE.NVAL_LOAD = 0 ; // reset

} // end of nearnbr_reset


// ========================================
void nearnbr_storeTrueTypes(void) {

  // Examine every TRUE TYPE in the training set and
  // store sparse list of each possible type.

  int NTYPE, i, TRUETYPE, MAP ;
  int TRUETYPE_LIST[MXTRUETYPE] ; 
  char fnam[] = "nearnbr_storeTrueTypes" ;

  // ------------- BEGIN -----------
  
  // initialize TRUETYPE storage arrays and counters
  NEARNBR_TRAINLIB.NTRUETYPE = NTYPE = 0;
  for(i=0; i< NTRUETYPE_MAX; i++ ) {
    NEARNBR_TRAINLIB.TRUETYPE_LIST[i] = -9;
    NEARNBR_TRAINLIB.NSN[i]           =  0;
  }

  for(i=0; i< MXTRUETYPE;    i++ ) 
    {NEARNBR_TRAINLIB.TRUETYPE_MAP[i]  = -9;}
  

  // loop over each entry in TRAILIB and store TRUETYPE
  
  for(i=0; i< NEARNBR_TRAINLIB.NTOT; i++ ) {
    TRUETYPE = NEARNBR_TRAINLIB.TRUETYPE[i] ;

    if ( TRUETYPE < 0 ) { continue ; }

    if ( TRUETYPE >= MXTRUETYPE ) {
      sprintf(c1err,"TRUETYPE = %d exceeds bound.", TRUETYPE   ) ;
      sprintf(c2err,"Check bound MXTRUETYPE = %d", MXTRUETYPE ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    MAP = NEARNBR_TRAINLIB.TRUETYPE_MAP[TRUETYPE] ;
    if ( MAP < 0 ) { 
      TRUETYPE_LIST[NTYPE] = TRUETYPE ;
      NTYPE++ ;
      NEARNBR_TRAINLIB.TRUETYPE_MAP[TRUETYPE]  = NTYPE ; // sparse inex
    }

  }  // i


  // sort list of types to be in increasing order since
  // they are random in the TRAINLIB.
  int isort, INDEX_SORT[MXTRUETYPE];
  int ORDER = +1 ;
  sortInt(NTYPE, TRUETYPE_LIST, ORDER, INDEX_SORT);
  for(i=0; i < NTYPE; i++ ) {
    isort = INDEX_SORT[i];
    TRUETYPE = TRUETYPE_LIST[isort] ;
    NEARNBR_TRAINLIB.TRUETYPE_LIST[i]        = TRUETYPE ;
    NEARNBR_TRAINLIB.TRUETYPE_MAP[TRUETYPE]  = i ; // update sparse map
  }

  NEARNBR_TRAINLIB.NTRUETYPE = NTYPE ;

  // summarize list of TRUE TYPES

  printf("\t Found %d TRUE-TYPEs in TRAINLIB: %s = ", 
	 NTYPE, NEARNBR_INPUTS.VARNAME_TRUETYPE );

  for ( i=0; i<NTYPE; i++ ) 
    { printf("%d ", NEARNBR_TRAINLIB.TRUETYPE_LIST[i] ); }

  printf("\n"); fflush(stdout);


  // ------------------------------------------------------
  // compute/store/print untrainined purity for each type
  // These un-trained purities are for visual reference,
  // but are not used in any calculations.

  double P, x1, x0 ;

  for(i=0; i< NEARNBR_TRAINLIB.NTOT; i++ ) {
    TRUETYPE = NEARNBR_TRAINLIB.TRUETYPE[i] ;
    if ( TRUETYPE < 0 ) { continue ; }

    // increment total number per TRUETYPE
    MAP = NEARNBR_TRAINLIB.TRUETYPE_MAP[TRUETYPE] ;
    NEARNBR_TRAINLIB.NSN[MAP]++ ;
    NEARNBR_TRAINLIB.NTOT_USE++ ;
  }

  x0 = (double)NEARNBR_TRAINLIB.NTOT_USE ; 
  for(i=0; i < NTYPE; i++ ) {
    
    x1 = (double)NEARNBR_TRAINLIB.NSN[i] ;
    P  = x1/x0 ;
    NEARNBR_TRAINLIB.UNTRAINED_PURITY[i] = P ;
    printf("\t Untrained Purity (TrueType=%3d) = %6.3f \n",
	   NEARNBR_TRAINLIB.TRUETYPE_LIST[i], P);
    fflush(stdout);
  }


} // end of nearnbr_storeTrueTypes


// ===========================
void nearnbr_makeHist(int ISPLIT) {

  // Feb 7 2017: 
  //  + include new histogram whose title contains
  //    file name with training sample. 

  char   TITLE[MXCHAR_FILENAME], ctmp[100];
  double xmin[2], xmax[2], xval[2], w ;
  int    ivar, NVAR, ibin, NTYPE, jtype, ID, NBIN[2];
  
  // --------------- BEGIN ------------

  if ( NN_APPLYFLAG  ) { return ; }

  NVAR  = NEARNBR_INPUTS.NVAR ;
  NTYPE = NEARNBR_TRAINLIB.NTRUETYPE ;

  printf("\t Create NEARNBR training histograms (ISPLIT=%d) \n", ISPLIT); 
  fflush(stdout);

  if ( ISPLIT > 1 ) { goto H2DANAL ; }

  // start with number of merged histograms
  ID = HIDOFF_NEARNBR + 1 ;
  sprintf(TITLE, "Number of merged files");
  xmin[0] = 0.5 ;  xmax[0] = 1.5 ;  NBIN[0] = 1 ;  xval[0] = w = 1.0 ;
  SNHIST_INIT(1, ID, TITLE, NBIN, xmin, xmax );
  SNHIST_FILL(1, ID, xval, w);


  // plot TRUETYPE map between sparse & absolute indices
  ID = HIDOFF_NEARNBR + 2 ;
  sprintf(TITLE,"TRUETYPE vs. SparseTRUETYPE");
  xmin[0] = -0.5 ;  xmax[0] = xmin[0] + (double)NTYPE ;  NBIN[0] = NTYPE ;
  SNHIST_INIT(1, ID, TITLE, NBIN, xmin, xmax );
  for ( jtype=0 ; jtype < NTYPE ; jtype++ ) {
    xval[0] = (double)jtype ;
    w       = (double)NEARNBR_TRAINLIB.TRUETYPE_LIST[jtype] ;
    SNHIST_FILL(1, ID, xval, w);
  }

  // book and fill SEPMAX values for each variable and isep bin
  ID = HIDOFF_NEARNBR + 3 ;
  sprintf(TITLE,"# ");
  for(ivar=0; ivar < NVAR ; ivar++ ) { 
    sprintf(ctmp,"%s  ", NEARNBR_INPUTS.VARNAMES[ivar] ); 
    strcat(TITLE,ctmp);
  } 
  strcat(TITLE," vs. SEPMAX bin");
  xmin[0] = -0.5 ;  xmax[0] = xmin[0] + (int)NBINTOT_SEPMAX_NEARNBR ;
  xmin[1] = -0.5 ;  xmax[1] = xmin[1] + (int)NVAR ;
  NBIN[0] = NBINTOT_SEPMAX_NEARNBR ;  NBIN[1] = NVAR;
  SNHIST_INIT(2, ID, TITLE, NBIN, xmin, xmax );

  for(ivar=0; ivar < NVAR ; ivar++ ) { 
    for(ibin=0; ibin < NBINTOT_SEPMAX_NEARNBR ; ibin++ ) {

      xval[0] = (double)ibin ;
      xval[1] = (double)ivar ;
      w       = NEARNBR_LIST_SEPMAX[ivar][ibin] ;
      SNHIST_FILL(2, ID, xval, w );

    } // NBINTOT
  }  // NVAR


  // - - - - -
  // plot untrained Purity vs. sparse type
  ID = HIDOFF_NEARNBR + 5 ;
  sprintf(TITLE,"Untrained Purity vs. SparseTRUETYPE");
  xmin[0] = -0.5 ;  xmax[0] = xmin[0] + (double)NTYPE ;  NBIN[0] = NTYPE ;
  SNHIST_INIT(1, ID, TITLE, NBIN, xmin, xmax );
  for ( jtype=0 ; jtype < NTYPE ; jtype++ ) {
    xval[0] = (double)jtype ;
    w       = (double)NEARNBR_TRAINLIB.UNTRAINED_PURITY[jtype] ;
    SNHIST_FILL(1, ID, xval, w);
  }


  // ----------------------------------------
  // now book plots that are filled in the NEARNBR_APPLY. For each
  // SNTYPE make 2D plot of TrainSet-TRUETYPE(train) vs. SEPBIN
  // Assume that the list of TYPEs in the training set is the
  // same as in the data set under analysis.
  // Note that the Train-type (vertical axis) is a sparse TYPE index,
  // and -1 means that there is not valid TYPE.

 H2DANAL:

  for ( jtype=0 ; jtype < NTYPE ; jtype++ ) {
    ID = HIDOFF_NEARNBR + 10 + jtype ;
    sprintf(TITLE,"TrainSet-SparseTYPE vs. SEPMAX BIN for SNTYPE=%d",
	    NEARNBR_TRAINLIB.TRUETYPE_LIST[jtype] ) ;
    xmin[0] = -0.5 ;  xmax[0] = xmin[0] + (int)NBINTOT_SEPMAX_NEARNBR ;
    xmin[1] = -1.5 ;  xmax[1] = xmin[1] + (int)(NTYPE+1) ;
    NBIN[0] = NBINTOT_SEPMAX_NEARNBR ;  NBIN[1] = NTYPE+1;
    SNHIST_INIT(2, ID, TITLE, NBIN, xmin, xmax );
  }


  // Apr 11 2019:
  // store total number of training events to read/print later.
  // This is just for information, not needed for training.
  xmin[0]=0.0; xmax[0]=1.0; NBIN[0]=1;
  ID = HIDOFF_NEARNBR + 40 ;  
  sprintf(TITLE, "Number of training events");
  SNHIST_INIT(1, ID, TITLE, NBIN, xmin, xmax );  

  // -------------------------------
  // store fileName with training (SIM) sample (Feb 7, 2017)
  // Due to stupid limit on hbook title length (char 80),
  // write three histograms to allow up to 240 chars.
  
  xmin[0]=0.0; xmax[0]=1.0; NBIN[0]=1;

  char strTmp[NSPLIT_TITLE][MXCHAR_TITLE+1]; 
  int i, LEN_TITLE, LONGNAME ;
  ID = HIDOFF_NEARNBR + 50 ;    
  nearnbr_TRAIN_FILENAME(0,TITLE);  // return TITLE = trainFile name
  LEN_TITLE = strlen(TITLE);

  // init each split histo-title to \0

  for(i=0; i < NSPLIT_TITLE; i++ ) {
    strTmp[i][0] = 0 ; 
    LONGNAME  = ( LEN_TITLE > (i+1)*MXCHAR_TITLE ) ;

    if ( LONGNAME ) { 
      strncpy(strTmp[i], &TITLE[MXCHAR_TITLE*i] , MXCHAR_TITLE); 
      strTmp[i][MXCHAR_TITLE] = 0; 
    }
    else if ( (LEN_TITLE - MXCHAR_TITLE*i) > 0 ) {
      sprintf(strTmp[i],"%s", &TITLE[MXCHAR_TITLE*i] ); 
    }

    // printf(" xxx strTmp[%d] = '%s' \n", i, strTmp[i] ); fflush(stdout);

  } // end i loop


  for(i=0; i < NSPLIT_TITLE ; i++ ) 
    { SNHIST_INIT(1, ID+i, strTmp[i], NBIN, xmin, xmax );  }


  // -----------------------------------------------------
  // finally, write name of variable which stores TRUETYPE
  ID = HIDOFF_NEARNBR + 60 ;  
  sprintf(TITLE, "%s", NEARNBR_INPUTS.VARNAME_TRUETYPE ); 
  SNHIST_INIT(1, ID, TITLE, NBIN, xmin, xmax );  

  // ---------
  NEARNBR_INPUTS.FILLHIST = 1;

  return ;

} // end of nearnbr_makeHist


// =====================================
void NEARNBR_APPLY(char *CCID) {

  // Main analysis driver to apply nearnbr.
  // Called by external program.
  // Jan 12 2017: few speed-up tricks for APPLY mode; see NN_APPLYFLAG

  int  ivar, itrain, LDMP, isep, i, NTRAIN_SUBSET ;
  int NTYPE      = NEARNBR_TRAINLIB.NTRUETYPE ;
  char fnam[] = "NEARNBR_APPLY" ;

  // ----------- BEGIN --------------

  NEARNBR_RESULTS_TRAIN.ITYPE     = -9 ;  // default is no evaluation
  NEARNBR_RESULTS_TRAIN.NCELL_TOT =  0 ;
  for(i=0; i < NEARNBR_TRAINLIB.NTRUETYPE ; i++ ) 
    { NEARNBR_RESULTS_TRAIN.NCELL[i]  = 0; }


  // make a few sanity checks.
  nearnbr_preAnal_verify();

  float   SQDIST ;
  int     TRUETYPE, TYPE_CUTPROB, isparse_TYPE, isubset, NNTOT ;
  int     NCUTDIST_TRAIN[MXTRUETYPE], NCUTDIST_FINAL[MXTRUETYPE] ;

  // for training only: get training subset inside largest SEPMAX sphere
  if ( NN_TRAINFLAG ) 
    { nearnbr_fill_SUBSET_TRAIN(CCID); }
  else
    { nearnbr_fill_SUBSET_APPLY(CCID); }


  NTRAIN_SUBSET = NEARNBR_TRAINLIB.NSUBSET ;

  // -----------------------------------
  // Two nested loops: 1) SEPMAX, 2) trainlib entries 
  // Count NTYPE for each SQSEP combo

  for(isep=0; isep < NBINTOT_SEPMAX_NEARNBR ; isep++ ) {

    for(i=0; i<NTYPE; i++ )  {
      NCUTDIST_TRAIN[i] = NCUTDIST_FINAL[i] = 0 ;
      NEARNBR_RESULTS_TRAIN.NCELL[i]  = 0 ;
      NEARNBR_RESULTS_FINAL.NCELL[i]  = 0 ;
    } 
    
    // - - - - - - - - - - - - - - - - - - - - - 
    NNTOT = 0 ;

    for(isubset=0; isubset < NTRAIN_SUBSET; isubset++ ) {
      itrain = NEARNBR_TRAINLIB.ITRAIN[isubset] ; 

      TRUETYPE  = NEARNBR_TRAINLIB.TRUETYPE[itrain] ;
      if ( TRUETYPE < 0 ) { continue ; }
      SQDIST    = nearnbr_SQDIST(isep,itrain) ;

      if ( SQDIST < 1.0 ) { 
	i         = NEARNBR_TRAINLIB.TRUETYPE_MAP[TRUETYPE]; // sparse index 
	if( i<0 || i >= NTYPE ) {
	  sprintf(c1err,"Invalid sparse index i=%d for TRUETYPE=%d", 
		  i, TRUETYPE);
	  sprintf(c2err, "isep=%d  itrain=%d  SQDIST=%f", 
		  isep, itrain, SQDIST);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
	}
	NCUTDIST_TRAIN[i]++ ; 
	NNTOT++ ;
	NEARNBR_RESULTS_TRAIN.NCELL[i]  = NCUTDIST_TRAIN[i] ;
	NEARNBR_RESULTS_FINAL.NCELL[i]  = NCUTDIST_FINAL[i] ; 
      } // SQDIST<1

    } // itrain
    // - - - - - - - - - - - - - - - - - - - - - 

    // analyze to get Type;
    // function returns TYPE_CUTPROB and isparse_TYPE
    isparse_TYPE = nearnbr_whichType(NTYPE, NCUTDIST_TRAIN, &TYPE_CUTPROB);
   
    // check to fill histogram(s)
    nearnbr_fillHist(isep,isparse_TYPE) ;
    
  }  // isep

  // ==============================================  

  // if just one sepmax bin, then evaluate and store nominal result.
  if ( NBINTOT_SEPMAX_NEARNBR == 1 ) {
    NEARNBR_RESULTS_TRAIN.ITYPE     = TYPE_CUTPROB ;
  } // 1 SEPMAX bin


  // -------------------------------------------
  // reset variables and counters for next SN
  nearnbr_reset();

  if ( NN_TRAINFLAG ) {  printf(" Done.\n"); fflush(stdout); }

  return ;

} // end of  NEARNBR_APPLY


void nearnbr_apply__( char *CCID) {  NEARNBR_APPLY(CCID) ; }

// =====================================================
float  nearnbr_SQDIST(int isep,int itrain) {

  // return NN distance

  int ivar ;
  int NVAR   = NEARNBR_INPUTS.NVAR ;
  float SQSEP, SQSEPMAX, SQDIST ;
  //  char fnam[] = "nearnbr_SQDIST" ;

  // --------------- BEGIN ---------------

  SQDIST = 0.0 ;
  for(ivar=0; ivar < NVAR; ivar++ ) {
    SQSEP    = NEARNBR_STORE.SQSEP[ivar][itrain] ;
    SQSEPMAX = NEARNBR_LIST_SQSEPMAX[ivar][isep] ;
    if ( SQSEP > SQSEPMAX) { return(99.0); } 
    SQDIST  += (SQSEP/SQSEPMAX);
  } // ivar

  return SQDIST ;

} // end of nearnbr_SQDIST

// ============================================
void nearnbr_preAnal_verify(void) {

  int ivar ;
  char fnam[] = "nearnbr_preAnal_verify" ;

  // -------------- BEGIN -------------

  // make sure that number of loaded variables matches expectation
  if ( NEARNBR_STORE.NVAL_LOAD != NEARNBR_INPUTS.NVAR ) {

    printf("\n PRE-ABORT DUMP: \n");
    for ( ivar=0; ivar < NEARNBR_INPUTS.NVAR ; ivar++ ) {
      printf("\t ivar=%d : %-16s  = %f (loaded %d times) \n"
	     ,ivar
	     ,NEARNBR_INPUTS.VARNAMES[ivar]
	     ,NEARNBR_STORE.VALUE_LOAD[ivar]
	     ,NEARNBR_STORE.NLOAD[ivar]
	     );
    }
    
    sprintf(c1err,"%d NN variables loaded by nearbnr_loadVal",
	    NEARNBR_STORE.NVAL_LOAD );
    sprintf(c2err,"but expected NVAR=%d", NEARNBR_INPUTS.NVAR );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

} // end of nearnbr_preAnal_verify


// ====================================
int nearnbr_whichType(int NTYPE, int *NCUTDIST,  int *TYPE_CUTPROB ) {

  // Determine TRUETYPE_CUTPROB.
  // Function args:
  //  (I) NTYPE             = total number of possible types
  //  (I) NCUTDIST[itype]   = Number with (R < 1) vs. itype
  //  (I) ISEPMAX           = SEPMAX bin (for plotting only)
  //  (O) TYPE_CUTPROB      = type passing PROB cut (-9 of no type found)
  //  (O) Function arg returns sparse TYPE index (-1 if no type)
  //
  //

  int   ISPARSE, ISPARSE_PROBMAX, NTOT, i  ;
  float PROB, PROBMAX, VAR_PROB, SIG_PROB, PROB4CUT ;
  float XN, XNTOT, XNTOT_CUBE ;

  //  char fnam[]  = "nearnbr_whichType" ;

  // ----------------- BEGIN ------------------

  *TYPE_CUTPROB = -9 ;  ISPARSE = ISPARSE_PROBMAX = -1 ;
  PROBMAX = -9.0 ;
  
  NTOT = 0 ;
  for(i=0; i < NTYPE; i++ ) { NTOT += NCUTDIST[i] ; }

  if ( NTOT == 0 ) { return ISPARSE ; }

  XNTOT      = (float)NTOT ;
  XNTOT_CUBE = XNTOT * XNTOT * XNTOT ;

  // find which type passes the CUTPROB cut
  for(i=0; i< NTYPE; i++ ) { 

    XN   = (float)NCUTDIST[i] ;
    PROB = XN / XNTOT ;

    if ( PROB > PROBMAX ) { PROBMAX=PROB; ISPARSE_PROBMAX=i; }
    
    VAR_PROB = XN*(XNTOT-XN) / XNTOT_CUBE ;

    // do not allow VAR_PROB = 1; at least 1 event counts toward error
    if ( VAR_PROB == 0.0 ) { VAR_PROB = 1.0/(XNTOT*XNTOT) ; }

    SIG_PROB = sqrtf(VAR_PROB) ;
    PROB4CUT = PROB - (SIG_PROB * NEARNBR_INPUTS.NSIGMA_PROB) ;

    if ( PROB4CUT  > NEARNBR_INPUTS.CUTPROB ) {
      *TYPE_CUTPROB = NEARNBR_TRAINLIB.TRUETYPE_LIST[i] ;
      ISPARSE = i ;
    }
  }

  // if nothing passes cut, at least return index of type
  // with max prob (Mar 2 2017)
  if ( ISPARSE < 0 ) { ISPARSE = ISPARSE_PROBMAX;  }
  
  return(ISPARSE) ;

} // end of nearnbr_whichType

// ===================================================
void nearnbr_fill_SUBSET_TRAIN(char *CCID) {

  // for training, store subset of training sample
  // that lies within largest SEPMAX sphere.
  // --> speeds up the training.
  //
  // Jan 29 2017: restore filling NEARNBR_STORE.SQSEP that got
  //              left out in the refactor earlier this month.
  //

  int NVAR       = NEARNBR_INPUTS.NVAR ;
  int isep, NSUBSET, itrain, ivar, TRUETYPE, NTRAIN ;
  double SQDIST, f_subset, VAL_DATA, VAL_TRAIN, SEP, SQSEP ;
  //  char   fnam[] = "nearnbr_fill_SUBSET_TRAIN" ;

  // ------------ BEGIN -----------

  isep    = NBINTOT_SEPMAX_NEARNBR - 1 ;
  NTRAIN  = NEARNBR_TRAINLIB.NTOT; 
  NSUBSET = 0 ;

  for ( itrain=0; itrain < NTRAIN; itrain++ ) {

    TRUETYPE  = NEARNBR_TRAINLIB.TRUETYPE[itrain] ;
    if ( TRUETYPE < 0 ) { continue ; }

    for(ivar=0; ivar < NVAR; ivar++ ) {
      VAL_DATA  = NEARNBR_STORE.VALUE_LOAD[ivar] ;
      VAL_TRAIN = NEARNBR_TRAINLIB.FITRES_VALUES[ivar][itrain] ;
      SEP       = VAL_DATA - VAL_TRAIN ;
      SQSEP     = SEP*SEP ;
      NEARNBR_STORE.SQSEP[ivar][itrain] = SQSEP ;
    } // end ivar

    SQDIST    = nearnbr_SQDIST(isep,itrain);      
    if ( SQDIST < 1.0 ) { 
      NEARNBR_TRAINLIB.ITRAIN[NSUBSET] = itrain;
      NSUBSET++ ;
    }
  }  // end of itrain loop


  NEARNBR_TRAINLIB.NSUBSET = NSUBSET ;

  f_subset = (double)NSUBSET/(double)NTRAIN ;
  printf("  Do NN train on CID=%s with TrueType=%d (subsetFrac=%.3f) ",
	 CCID, NEARNBR_STORE.TRUETYPE_LOAD, f_subset );
  fflush(stdout); 

  return ;

} // end of nearnbr_fill_SUBSET_TRAIN


// ===================================================
void nearnbr_fill_SUBSET_APPLY(char *CCID) {

  int isep = 0;  
  int NTRAIN_TOT = NEARNBR_TRAINLIB.NTOT; 
  int NVAR       = NEARNBR_INPUTS.NVAR ;
  int itrain, NVAR_NEAR, ivar, LDMP,  NTRAIN_SUBSET ;
  float VAL_DATA, VAL_TRAIN, SEP, SQSEP, SQRATIO, SQSEPMAX ;
  float VAL_ARRAY[MXVAR_NEARNBR];
  char fnam[] = "nearnbr_fill_SUBSET_APPLY" ;

  // ------------- BEGIN -------------

  // check option to use pre-defined cubes defined in CELMAP

  int ICELL_1D, NLIST, ilist ;
  if ( NEARNBR_CELLMAP.DOFLAG ) {
    
    getInfo_CELLMAP(1, NEARNBR_STORE.VALUE_LOAD  ,  // (I)
		    &ICELL_1D )  ;                  // (O)

    NLIST = NEARNBR_CELLMAP.NLIST[ICELL_1D] ;	
    for(ilist=0; ilist < NLIST; ilist++ ) {
      itrain = NEARNBR_CELLMAP.ITRAIN_LIST[ICELL_1D][ilist] ;
      NEARNBR_TRAINLIB.ITRAIN[ilist] = itrain ;

      if ( itrain >= NTRAIN_TOT || itrain < 0 ) {
	sprintf(c1err,"itrain=%d exceeds NTRAIN_TOT=%d", itrain, NTRAIN_TOT);
	sprintf(c2err,"Check CELLMAP functions" );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

      for(ivar=0; ivar < NVAR; ivar++ ) {
	VAL_DATA  = NEARNBR_STORE.VALUE_LOAD[ivar] ;
	VAL_TRAIN = NEARNBR_TRAINLIB.FITRES_VALUES[ivar][itrain] ;
	SEP       = VAL_DATA - VAL_TRAIN ;
	SQSEP     = SEP*SEP ;
	NEARNBR_STORE.SQSEP[ivar][itrain] = SQSEP ;
      } // end ivar

    } // end ilist/itrain

    NEARNBR_TRAINLIB.NSUBSET = NLIST ;

    return ;
  }

  // slower method; loop over every training event
  // and keep subset withing SEPMAX cube.

  NTRAIN_SUBSET = 0 ;
  for ( itrain=0; itrain < NTRAIN_TOT; itrain++ ) {
    NVAR_NEAR = 0 ;
    for(ivar=0; ivar < NVAR; ivar++ ) {
      VAL_DATA  = NEARNBR_STORE.VALUE_LOAD[ivar] ;
      VAL_TRAIN = NEARNBR_TRAINLIB.FITRES_VALUES[ivar][itrain] ;
      SEP       = VAL_DATA - VAL_TRAIN ;
      SQSEP     = SEP*SEP ;
      NEARNBR_STORE.SQSEP[ivar][itrain] = SQSEP ;

      SQSEPMAX = NEARNBR_LIST_SQSEPMAX[ivar][isep] ;
      SQRATIO = SQSEP / SQSEPMAX ;
      if ( SQRATIO > 1.0001 ) 
	{ ivar=NVAR ; } // end ivar loop
      else
	{ NVAR_NEAR++ ; } // increment number of near-training variables
    
      LDMP = (itrain < -5 ) ;
      if ( LDMP ) {
	printf(" xxx %3d: Delta(%s) = %6.3f - %6.3f = %f \n",
	       itrain, NEARNBR_INPUTS.VARNAMES[ivar], 
	       VAL_DATA, VAL_TRAIN, SEP);  fflush(stdout);
      }

    }  // ivar

    // Jan 2017: update sparse subset for which each variable
    //           satisfies SEP < SEPMAX constraint
    if ( NVAR_NEAR == NVAR  ) {
      NEARNBR_TRAINLIB.ITRAIN[NTRAIN_SUBSET] = itrain ;
      NTRAIN_SUBSET++ ;
      NEARNBR_TRAINLIB.NSUBSET = NTRAIN_SUBSET ;
    }

  } // itrain


  return ;

} // end nearnbr_fill_SUBSET_APPLY

// ===================================================
void nearnbr_fillHist(int ISEPMAX, int ISPARSE_TYPE) {

  // ----------- histo-filling ----------
  double  xval[2],   w1 = 1.0 ;
  int     HID ;

  if ( NEARNBR_INPUTS.FILLHIST == 0 ) { return ; }

  HID     = HIDOFF_NEARNBR + 10 + NEARNBR_STORE.TRUETYPE_SPARSE ;
  xval[0] = (double)ISEPMAX ;
  xval[1] = (double)ISPARSE_TYPE ;   // points to trained type
  SNHIST_FILL(2, HID, xval, w1 ) ;

  // Apr 11 2019: load total number of training events
  if ( ISEPMAX == 0 ) {
    HID     = HIDOFF_NEARNBR + 40;
    xval[0] = (double)0.5 ;
    SNHIST_FILL(1, HID, xval, w1 ) ;
  }

} // end of nearnbr_fillHist



// ============================================================
void NEARNBR_GETRESULTS(char *CCID, int *ITYPE_BEST, 
			int *NTYPE, int *ITYPE_LIST,
			int *NCELL_TRAIN_LIST )
{

  // Apr 7 2014
  // Return results of NN analysis ONLY if there is one SEPMAX bin;
  // with >1 SEPMAX bin, need to analyze the histograms with
  // 'nearnbr_maxFoM.exe' to find the optimal sepmax coefficients.
  //
  // Return 
  //   ITYPE_BEST = best integer type (-1 => no type, -9 => not evaluated)
  //   NTYPE      = total number of true types
  //   ITYPE_LIST = list of true types
  //   NCELL_TRAIN_LIST = list of NCELL for each true type, training cuts
  //   NCELL_FINAL_LIST = idem with NN cuts
  //
  // If number of SEPMAX bins > 1, then return ITYPE = -9.
  //
  // Feb 2016: re-write to return NCELL for each type.
  // Jun 2016: refactor to return NCELL_FINAL_LIST.
  //           NCELL_TRAIN_LIST = old NCELL_LIST.
  //

  int  i, LDMP  ;
  //  char fnam[] = "NEARNBR_GETRESULTS" ;

  // ----------- BEGIN --------------

  LDMP = 0 ; // ( NN_APPLYFLAG == 1 && NEARNBR_RESULTS_FINAL.NCELL[0]>0 ) ;

  if ( LDMP ) { printf(" xxx ---------------------------- \n"); }

  *ITYPE_BEST  = NEARNBR_RESULTS_TRAIN.ITYPE ;
  *NTYPE       = NEARNBR_TRAINLIB.NTRUETYPE ;

  for(i=0; i < *NTYPE; i++ ) {
    ITYPE_LIST[i]         = NEARNBR_TRAINLIB.TRUETYPE_LIST[i] ;
    NCELL_TRAIN_LIST[i]   = NEARNBR_RESULTS_TRAIN.NCELL[i] ;
  }


} // end of NEARNBR_GETRESULTS


void nearnbr_getresults__(char *CCID, int *ITYPE, int *NTYPE, 
			  int *ITYPE_LIST,
			  int *NCELL_TRAIN_LIST) {   
  NEARNBR_GETRESULTS(CCID, ITYPE, NTYPE,
		     ITYPE_LIST, NCELL_TRAIN_LIST );
}
