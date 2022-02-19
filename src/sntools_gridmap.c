/*****************************************
  July 2021: extracted from sntools.c[h]

  Utilities to read and store multi-dimensional map in GRIDMAP struct,
  and use map for multi-D interpolation.

  Note: these utils are NOT related to those in sntools_modelgrid_gen.c[h]
        and sntools_modelgrid_read.c[h]

 *****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "sntools.h"


// ********************************************************
void read_GRIDMAP(FILE *fp, char *MAPNAME, char *KEY_ROW, char *KEY_STOP, 
		  int IDMAP, int NDIM, int NFUN, int OPT_EXTRAP, int MXROW,
                  char *callFun, GRIDMAP *GRIDMAP_LOAD ) {

  // Mar 2019
  // Utility to read mutil-D map from file and call
  // init_inter_GRIDMAP to create & store GRIDMAP_LOAD.
  // Beware that uniform map-bins are strictly enforce;
  // ABORTs on non-uniform bin.
  //
  // Inputs:
  //   *fp        : already-opened file to read
  //  MAPNAME     : name of map
  //  KEY_ROW     : NVAR columns follows this row-key
  //  KEY_STOP    : stop reading when this key is reached;
  //              " default is to stop reading on blank line.
  //  IDMAP       : integer ID of GRIDMAP_LOAD
  //  NDIM        : number of dimensions of map
  //  NFUN        : number of functions of map
  //  OPT_EXTRAP  : flag for extrapolation outside map range
  //                1-> extrap, 0->return error, -1->abort outside range
  //  MXROW       : abort if NROW > MXROW 
  //  callFun     : name of calling function (for error messages)
  //
  // Output:
  //    GRIDMAP_LOAD
  //
  //
  // Apr 12 2019: abort if 10 or more rows read without valid key
  // Jun 12 2020: pass MAPNAME as input arg.
  // Mar 27 2021: fix MSKOPT to allow comments with commas
  // May 13 2021: require grid uniformity at 1E-6 instead of E-3

  int   READ_NEXTLINE = 1 ;
  int   NROW_READ     = 0 ;
  int   NROW_UPDATE   = 500000; // update reading after this many
  int   NVARTOT = NDIM + NFUN;
  char  *VARLIST = GRIDMAP_LOAD->VARLIST ; // for comment only
  double DUMVAL = -999.0 ;

  int   MEMD   = sizeof(double);
  int   MEMVAR = NVARTOT  * sizeof(double*);
  int   MEMROW = MXROW    * MEMD ;

  int  NBADBIN = 0 ;
  int  NLINE   = 0 ;
  int  MSKOPT =
    MSKOPT_PARSE_WORDS_STRING        + 
    MSKOPT_PARSE_WORDS_IGNORECOMMENT +  // ignore comments on valid rows
    MSKOPT_PARSE_WORDS_IGNORECOMMA ;    // allow commas in comments

  double **TMPMAP2D ;  // [0:NVARTOT-1][MXROW-1]
  double *TMPVAL, *TMPVAL_LAST, *DIFVAL_LAST, DDIF, DIF;

  int   ivar, NWD, ISKEY_ROW, EXTRA_WORD_OK ;
  int   LDIF1, LDIF2, ivar2, NROW_SKIP=0 ;
  char  LINE[200], word[40] ;
  char fnam[] = "read_GRIDMAP" ;
 
  // ----------- BEGIN -------------

  // create generic MAPNAME using row key and IDMAP

  // allocate arrays to monitor uniform binning.
  TMPVAL      = (double*) malloc(NVARTOT * MEMD );
  TMPVAL_LAST = (double*) malloc(NVARTOT * MEMD );
  DIFVAL_LAST = (double*) malloc(NVARTOT * MEMD );
  for(ivar=0; ivar<NVARTOT; ivar++) {
    TMPVAL[ivar] = DUMVAL;
    TMPVAL_LAST[ivar] = DUMVAL;
    DIFVAL_LAST[ivar] = DUMVAL;
  }

  // alloate temp 2D array to read map
  TMPMAP2D = (double**) malloc(MEMVAR);
  for(ivar=0; ivar<NVARTOT; ivar++) {TMPMAP2D[ivar]=(double*)malloc(MEMROW);} 


  while ( READ_NEXTLINE ) {
    LINE[0] = 0 ;
    fgets(LINE,200,fp);  NLINE++ ;
    NWD = store_PARSE_WORDS(MSKOPT,LINE);

    // abort if we read too many lines without finding any valid row keys
    if ( NLINE > 20 && NROW_READ==0 ) {
      sprintf(c1err,"Found no '%s' keys after reading %d lines.",
	      KEY_ROW, NLINE);
      sprintf(c2err,"NDIM=%d, NFUN=%d, callFun=%s", 
	      NDIM, NFUN, callFun );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    // Skip blank line.
    // However, stop reading only after reading at least one valid row;
    // this allows blank line(s) between VARNAMES and first map row.
    if ( NROW_READ > 0  && NWD == 0 )  { READ_NEXTLINE=0; }
    if ( NWD == 0 ) { continue ; }

    get_PARSE_WORD(0,0,word);  

    ISKEY_ROW = 0 ;
    if ( strcmp(word,KEY_ROW) ==0 ) { ISKEY_ROW = 1; }
    if ( strcmp(word,KEY_STOP)==0 ) { READ_NEXTLINE=0; continue; }

    if ( ISKEY_ROW ) {
      
      if ( NROW_READ > 0 && (NROW_READ % NROW_UPDATE)==0 ) {
	printf("\t Reading %s row %8d \n", MAPNAME, NROW_READ );
      }

      NROW_SKIP = 0 ;
      // allow comment string on same line as grid data
      EXTRA_WORD_OK = 1 ;
      if ( NWD-1 > NVARTOT ) {
	get_PARSE_WORD(0,NVARTOT+1,word);
	EXTRA_WORD_OK = ( word[0] == '#' ) ;
      }
      //  printf(" xxx extra word = '%s'  OK=%d \n",word, EXTRA_WORD_OK);

      if ( (NWD-1 < NVARTOT) || (!EXTRA_WORD_OK) ) {
	print_preAbort_banner(fnam);
	printf("  MAPNAME = '%s' \n", MAPNAME);
	printf("  LINE = '%s' \n", LINE);
	sprintf(c1err,"Expected NVARTOT=%d words after '%s' key,",
		NVARTOT, KEY_ROW);
	sprintf(c2err,"but found %d. (NDIM=%d, NFUN=%d, ROW=%d, callFun=%s)", 
		NWD-1, NDIM, NFUN, NROW_READ, callFun );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }

      for(ivar=0; ivar < NVARTOT; ivar++ ) {
	get_PARSE_WORD(0,1+ivar,word);
	sscanf ( word, "%le", &TMPVAL[ivar] );
	TMPMAP2D[ivar][NROW_READ] = TMPVAL[ivar];

	// check for uniform binning
	DIF = TMPVAL[ivar] - TMPVAL_LAST[ivar];
	if ( DIF > 0.0  && ivar < NDIM && TMPVAL_LAST[ivar]!=DUMVAL ) { 
	  DDIF  = DIF - DIFVAL_LAST[ivar] ;
	  LDIF1 = ( fabs(DDIF/DIF) > 1.0E-6 ) ; 
	  LDIF2 = ( DIFVAL_LAST[ivar] > 0.0 ) ;
	  if ( LDIF1 && LDIF2 ) {
	    NBADBIN++ ;
	    printf(" ERROR: non-uniform bin at '%s'=", VARLIST );
	    for(ivar2=0; ivar2 < NVARTOT; ivar2++ ) 
	      { printf("%.3f ", TMPVAL[ivar2] ); }
	    printf(" (row %d)\n", NROW_READ );	    fflush(stdout);
	  }
	  DIFVAL_LAST[ivar] = DIF; 
	} // end DIF>0
	// end uniform bin check

	TMPVAL_LAST[ivar] = TMPVAL[ivar];
      } // end ivar loop

      NROW_READ++ ;
      if ( NROW_READ >= MXROW ) {
	sprintf(c1err,"NROW_READ=%d exceeds MXROW=%d", NROW_READ, MXROW);
	sprintf(c2err,"NDIM=%d  NFUN=%d  callFun=%s", 
		NDIM, NFUN, callFun );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }

    } // end KEY_ROW
    else {
      // 4.2019: abort if too many rows have invalid key
      NROW_SKIP++ ;
      if ( NROW_SKIP >= 10 ) { 
	print_preAbort_banner(fnam);  
	printf("   Last line read: %s\n", LINE);
	sprintf(c1err,"Read %d rows without valid row-key, "
		"stop-key, or blank line.", NROW_SKIP );
	sprintf(c2err,"KEY_ROW='%s'  KEY_STOP='%s'  callFun=%s", 
		KEY_ROW, KEY_STOP, callFun );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }
    }

  } // end while


  // -------------------------------------------------
  // ABORT on non-uniform bins
  if ( NBADBIN > 0 ) {
    sprintf(c1err,"%d non-uniform bin errors detected", NBADBIN);
    sprintf(c2err,"Check %s map. ", KEY_ROW );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  // ----------------
  printf("    Load GRIDMAP-%3.3d '%s(%s)'  NROW=%d \n",
	 IDMAP, MAPNAME, VARLIST, NROW_READ); fflush(stdout);

  init_interp_GRIDMAP(IDMAP, MAPNAME, NROW_READ, NDIM, NFUN, OPT_EXTRAP,
		      TMPMAP2D, 
		      &TMPMAP2D[NDIM],
		      GRIDMAP_LOAD  );       // <== returned

  // free temp memory
  for(ivar=0; ivar < NVARTOT; ivar++ )  { free(TMPMAP2D[ivar]); }
  free(TMPMAP2D);
  free(TMPVAL); free(TMPVAL_LAST); free(DIFVAL_LAST);
  return ;

} // end read_GRIDMAP


// ==============================================================
void malloc_GRIDMAP(int OPT, GRIDMAP *gridmap, int NFUN, int NDIM, int MAPSIZE) {

  // Created May 26 2021
  // OPT > 0 -> malloc
  // OPT < 0 -> free

  int ifun;
  int I4  = sizeof(int) ;
  int I8  = sizeof(double) ;
  int I8p = sizeof(double*) ;
  char string[12];
  char fnam[] = "malloc_GRIDMAP";

  // ---------------- BEGIN ----------- 

  if ( OPT > 0 ) {
    sprintf(string,"allocate");
    gridmap->NBIN      = (int     *)malloc(I4*NDIM+I4);
    gridmap->VALMIN    = (double  *)malloc(I8*NDIM+I8);
    gridmap->VALMAX    = (double  *)malloc(I8*NDIM+I8);
    gridmap->VALBIN    = (double  *)malloc(I8*NDIM+I8);
    gridmap->RANGE     = (double  *)malloc(I8*NDIM+I8);
    gridmap->FUNMIN    = (double  *)malloc(I8*NFUN);
    gridmap->FUNMAX    = (double  *)malloc(I8*NFUN);
    gridmap->INVMAP    = (int     *)malloc(I4*MAPSIZE+I4);

    gridmap->FUNVAL    = (double **)malloc(I8p*NFUN);
    for(ifun=0; ifun < NFUN; ifun++ ) 
      {  gridmap->FUNVAL[ifun] = (double *)malloc(I8*MAPSIZE);  }
  }
  else {
    sprintf(string,"free GRIDMAP %d ", gridmap->ID );

    free(gridmap->NBIN);
    free(gridmap->VALMIN);
    free(gridmap->VALMAX);
    free(gridmap->VALBIN);
    free(gridmap->RANGE);
    free(gridmap->FUNMIN);
    free(gridmap->FUNMAX);
    free(gridmap->INVMAP);

    for(ifun=0; ifun < NFUN; ifun++ ) { free(gridmap->FUNVAL[ifun]); }
    free(gridmap->FUNVAL);
  }

  printf("\t %s: %s\n", fnam, string);
  fflush(stdout);

  return ;

} // end malloc_GRIDMAP

// ==============================================================
void init_interp_GRIDMAP(int ID, char *MAPNAME, int MAPSIZE, 
			 int NDIM, int NFUN, int OPT_EXTRAP,
			 double **GRIDMAP_INPUT, double **GRIDFUN_INPUT,
			 GRIDMAP *gridmap ) {

  // Created July 2011 by R.Kessler
  // Return struct *gridmap to assist in mult-dimensional interp .
  // This struct contains all the binning info for each dimension.
  //
  // Arguments
  // (I) ID        reference id:  grep IDGRIDMAP_ sntools.h | grep define
  // (I) MAPNAME   human-readable name for error message
  // (I) MAPSIZE   total number of bins in gridmap
  // (I) NDIM      number of dimensions = number of variables
  // (I) NFUN      Number of functions on same GRID
  // (I) OPT_EXTRAP  1=>extrap, 0=>return error, -1=>ABORT
  // (I) **GRIDMAP_INPUT[idim][i=0 to MAPSIZE-1] = grid vals; e.g., abscissa 
  // (I) **GRIDFUN_INPUT[ifun][i=0 to MAPSIZE-1] = function values
  // (O) *gridmap  structure to return 
  //         (to be passed later to interp_GRIDMAP)
  //
  // Feb 12 2018: 
  //   + malloc gridmap here, instead of externally
  //   + refactor so that all local indices are 0 to N-1 (not 1-N)
  //
  // Mar 13 2018:
  //   + fix bug malloc-ing FUNVAL : I8p -> I8p * NFUN
  //
  // May 26 2021: move malloc calls into malloc_GRIDMAP()
  //

  int idim, ifun, i, NBIN, igrid_tmp, igrid_1d[100] ;
  double VAL, VALMIN, VALMAX, VALBIN, LASTVAL, RANGE, DIF ;
  double FUNVAL, RANGE_CHECK, RATIO ;
  char fnam[] = "init_interp_GRIDMAP" ;

  // --------- BEGIN ------------

  malloc_GRIDMAP(+1, gridmap, NFUN, NDIM, MAPSIZE) ;

  VALBIN = 0.0 ;
  for ( idim=0; idim < NDIM; idim++ ) {

    NBIN    =  0 ;
    VALMIN  = +1.0E12 ;
    VALMAX  = -1.0E12 ;
    LASTVAL = -999. ;

    for ( i=0; i < MAPSIZE; i++ ) {
      VAL = GRIDMAP_INPUT[idim][i] ;

      if ( VAL > VALMAX  ) { VALMAX = VAL ; }
      if ( VAL < VALMIN  ) { VALMIN = VAL ; }
      if ( VAL > LASTVAL && LASTVAL != -999. ) 
	{ VALBIN = VAL - LASTVAL ; }

      LASTVAL = VAL ;
    } 

    RANGE = VALMAX - VALMIN ;
    if ( RANGE > 1.0E-9 ) 
      { NBIN = (int)( (RANGE+0.001*VALBIN) / VALBIN ) + 1; }
    else
      { NBIN = 1; }

    // load output struct
    gridmap->ID           = ID ;
    gridmap->NDIM         = NDIM ;
    gridmap->NBIN[idim]   = NBIN ;
    gridmap->VALMIN[idim] = VALMIN ;
    gridmap->VALMAX[idim] = VALMAX ;
    gridmap->VALBIN[idim] = VALBIN ;
    gridmap->RANGE[idim]  = RANGE ;
    gridmap->NFUN         = NFUN ;
    gridmap->NROW         = MAPSIZE ;
    gridmap->OPT_EXTRAP   = OPT_EXTRAP ;

    // make sure that VALBIN x integer = RANGE
    RANGE_CHECK = (double)(NBIN-1) * VALBIN;
    RATIO       = RANGE_CHECK/RANGE ;

    if ( fabs(RATIO-1.0) > 1.0E-4 ) {
      print_preAbort_banner(fnam);
      printf("\t VALMAX - VALMIN  = %le (%le to %le)\n", 
	     RANGE, VALMIN, VALMAX );
      printf("\t (NBIN-1)*BINSIZE = %le (%d x %le) \n",
	     RANGE_CHECK, NBIN-1, VALBIN);
      printf("\t Ratio-1 = %le \n", RATIO-1. );

      sprintf(c1err,"Non-uniform binning for idim=%d", idim);
      sprintf(c2err,"Check map = '%s' ", MAPNAME );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

  } // idim

  // now store function value at each node along with
  // mapping between mult-D indices and 1D index
  // Load gridmap->FUNVAL  and gridmap->IGRIDMAP

  init_1DINDEX(ID, NDIM, &gridmap->NBIN[0] ) ;

  for ( ifun=0; ifun < NFUN; ifun++ ) { 
    gridmap->FUNMIN[ifun] = +999999.0 ;
    gridmap->FUNMAX[ifun] = -999999.0 ;
  }


  for ( i=0; i < MAPSIZE; i++ )  {  

    for ( ifun=0; ifun < NFUN; ifun++ ) { 
      FUNVAL = GRIDFUN_INPUT[ifun][i] ;
      gridmap->FUNVAL[ifun][i] = FUNVAL ;
      if(FUNVAL < gridmap->FUNMIN[ifun]) { gridmap->FUNMIN[ifun] = FUNVAL; }
      if(FUNVAL > gridmap->FUNMAX[ifun]) { gridmap->FUNMAX[ifun] = FUNVAL; }
    }

      for ( idim=0; idim < NDIM; idim++ ) {      
	VAL    = GRIDMAP_INPUT[idim][i] ;
	VALMIN = gridmap->VALMIN[idim] ;
	VALBIN = gridmap->VALBIN[idim] ;
	DIF    = VAL - VALMIN ;
	if ( VALBIN < 1.0E-9 ) 
	  { igrid_1d[idim] = 0 ; }
	else
	  { igrid_1d[idim] = (int)((DIF+1.0E-9)/VALBIN); }  //  + 1 ; }

      }

      igrid_tmp = get_1DINDEX(ID, NDIM, &igrid_1d[0] ) ; 

      if ( igrid_tmp < 0 || igrid_tmp >= MAPSIZE ) {
	print_preAbort_banner(fnam);
	printf("   MAPNAME=%s: \n", MAPNAME );
	for ( idim=0; idim < NDIM; idim++ ) {  
	  VAL    = GRIDMAP_INPUT[idim][i] ;
	  VALMIN = gridmap->VALMIN[idim] ;
	  VALBIN = gridmap->VALBIN[idim] ;
	  printf("   idim=%d : VAL=%10.3f  BIN=%10.3f  MIN=%10.3f  "
		 "igrid_1d=%d \n",
		 idim, VAL, VALBIN, VALMIN, igrid_1d[idim] );
	}
	printf("\t Probably have non-uniform binning.\n");
	       
	sprintf(c1err,"Invalid igrid_tmp=%d (ID=%d, MAPSIZE=%d)", 
		igrid_tmp, ID, MAPSIZE );
	sprintf(c2err,"original NDIM=%d  igrid=%d ", NDIM, i ) ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }

      gridmap->INVMAP[igrid_tmp] = i ;

  } // end loop over MAPSIZE
  
  return ;

} // end of init_interp_GRIDMAP


int interp_GRIDMAP(GRIDMAP *gridmap, double *data, double *interpFun ) {

  // Created Jul 3, 2011
  // Do multi-dimensional interpolation.
  // (I) *gridmap is returned from init_GRIDMAP
  // (I) *data is the multi-dimensional data point to interpolate
  // (O) *interpFun is the interpolated function value, or array
  //                of function values for multiple functions
  //
  // Function returns  0  if *data is withing the grid;
  // Function returns -1  if *data is outside the grid.
  //
  // Note that init_interp_GRIDMAP must be called first
  // to initialize the gridmap structure.
  //
  // Jul 25, 2011; remove   WGT_SUM <= 0 test since it can be
  //               zero if the function is zero nearby.
  //
  // Aug 28, 2011: 
  //     fix bug for when TMPVAL is within "EPSILON" of TMPMAX
  //
  // Mar 14 2019: 
  //  + check OPT_EXTRAP option
  //  + return SUCCESS or ERROR instead of hard-coded values.
  //
  // Mar 15 2020: allow numerical glitches in TMPMIN and TMPMAX

  int 
    ivar, ifun, NFUN, NVAR, ID, igrid, MSK, NBIN, OPT_EXTRAP
    ,NCORNERS, icorner, igrid_tmp, igrid_1D, g
    ,igrid_cell[100], igrid_var[100], IGRID_VAR[100]
    ;

  double  
    WGT_SUM[100], CORNER_WGTSUM, CORNER_WGT
    ,TMPVAL, TMPDIF, TMPMIN, TMPMAX, TMPBIN, TMPRANGE, xgrid, XNBIN
    ,GRIDFRAC[100], FUNVAL[100]
    ;

  double EPSILON = 1.0E-8 ;

  int  LDMP=0 ;
  bool outside_bound, too_lo, too_hi ;
  char fnam[] = "interp_GRIDMAP" ;

  // ---------- BEGIN ------------

  ID   = gridmap->ID ;
  NVAR = gridmap->NDIM ;
  NFUN = gridmap->NFUN ;
  OPT_EXTRAP = gridmap->OPT_EXTRAP ;

  // Mar 27 2021: check trivial case with NDIM=1 and 1 bin
  if ( NFUN==1 && NVAR == 1 && gridmap->NBIN[0]==1 ) {
    ifun = igrid=0;
    igrid_tmp    = get_1DINDEX( ID, NVAR, &igrid);
    igrid_1D     = gridmap->INVMAP[igrid_tmp] ;
    interpFun[0] = gridmap->FUNVAL[ifun][igrid_1D];
    return(SUCCESS);
  }

  for  ( ifun=0; ifun < NFUN; ifun++ )   {  
    interpFun[ifun] = 0.0 ; 
    WGT_SUM[ifun] = 0.0 ;
  }
  CORNER_WGTSUM = 0.0 ;


  if ( LDMP ) 
    { printf(" xxxxx ------------- START DUMP ----------------- \n"); }

  // get central index and grid-frac in each dimension
  for ( ivar=0; ivar < NVAR; ivar++ ) {
    TMPVAL   = data[ivar] ;
    TMPMIN   = gridmap->VALMIN[ivar] ;
    TMPMAX   = gridmap->VALMAX[ivar] ;
    TMPBIN   = gridmap->VALBIN[ivar] ;
    TMPRANGE = TMPMAX - TMPMIN ;

    // Mar 15 2020: allow numerical glitches
    TMPMAX += (1.0E-14*TMPRANGE);
    TMPMIN -= (1.0E-14*TMPRANGE);

    too_lo        = ( TMPVAL < TMPMIN ) ;
    too_hi        = ( TMPVAL > TMPMAX ) ;
    outside_bound = ( too_lo || too_hi );

    if ( outside_bound ) {
      // check extrap option
      if ( OPT_EXTRAP > 0 ) {
	if ( too_lo ) { TMPVAL = TMPMIN + (TMPRANGE*1.0E-12); }
	if ( too_hi ) { TMPVAL = TMPMAX - (TMPRANGE*1.0E-12); }
      }
      else if ( OPT_EXTRAP < 0 ) {
	// ??
      }
      else 
	{ return(ERROR); }
	
    } // end outside_bound


    TMPDIF  = TMPVAL - TMPMIN ;
    if ( TMPBIN == 0.0 )
      { XNBIN = 0.0 ; igrid = 0; }
    else if ( (TMPMAX - TMPVAL)/TMPRANGE < EPSILON  )  { 
      XNBIN = (TMPDIF - TMPRANGE*EPSILON)/TMPBIN ;
      igrid = (int)(XNBIN) ; 
    }
    else {
      XNBIN = (TMPDIF + TMPRANGE*EPSILON ) / TMPBIN ;
      igrid = (int)(XNBIN); //  + 1; 
    }

    xgrid   = (double)igrid ;

    // store relative cell location: 0-1
    if ( TMPBIN > 0.0 ) 
      {  GRIDFRAC[ivar]  = TMPDIF/TMPBIN - xgrid ; }
    else
      {  GRIDFRAC[ivar]  = 1.0 ; }

    IGRID_VAR[ivar] = igrid ; // store central bin  for each var

    if ( LDMP ) {
      printf(" xxxx VAL=%f  BIN=%f  XNBIN=%f  igrid=%2d GRIDFRAC=%le \n",
	     TMPVAL, TMPBIN, XNBIN, igrid, GRIDFRAC[ivar] );
      fflush(stdout);
    }

  } // ivar


  // determine the grid points at the corners of the
  // NVAR-dimentional cell containing *galpar.
  // Then take weighted average of WGTMAP at each corner.

  double XN ;
  XN = (double)NVAR ;

  NCORNERS = (int)pow(2.0,XN);
  for ( icorner=0; icorner < NCORNERS; icorner++ ) {
   
    if ( LDMP ) {
      printf(" xxx --------- Next icorner = %d / %d ------------ \n", 
	     icorner, NCORNERS ); 
      fflush(stdout);
    }

    CORNER_WGT = 1.0 ;
    for ( ivar=0; ivar < NVAR; ivar++ ) {
      //      MSK = 1 << (ivar-1);
      MSK = 1 << (ivar);
      igrid_cell[ivar] = (icorner & MSK)/MSK ; // 0 or 1 only
      igrid_var[ivar]  = IGRID_VAR[ivar] + igrid_cell[ivar];

      if ( LDMP  ) {
	printf("\t xxxxxx ivar=%d : cell=%d  igrid_var=%d \n",
	       ivar, igrid_cell[ivar], igrid_var[ivar] );
	fflush(stdout);
      }

      // make sure that igrid_var is valid
      NBIN = gridmap->NBIN[ivar] ;   g = igrid_var[ivar];
      if ( g < 0 || g >= NBIN  ) {
	TMPVAL  = data[ivar] ;
	TMPMIN  = gridmap->VALMIN[ivar] ;
	TMPMAX  = gridmap->VALMAX[ivar] ;

	print_preAbort_banner(fnam);
	printf("\t ID=%d  NBIN=%d  cell=%d icorner=%d\n",
	       ID, NBIN, igrid_cell[ivar], icorner);
	printf("\t grep IDGRIDMAP $SNANA_DIR/src/sntools_gridmap.h");
	sprintf(c1err, "Invalid igrid_var[ivar=%d]=%d", ivar, g);
	sprintf(c2err, "VAL=%f  VALMIN/MAX = %f / %f", 
		TMPVAL, TMPMIN, TMPMAX);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      
      if ( igrid_cell[ivar] ) 
	{ CORNER_WGT *= GRIDFRAC[ivar] ; }
      else
	{ CORNER_WGT *= (1.0 - GRIDFRAC[ivar]) ; }
    }

    CORNER_WGTSUM += CORNER_WGT ;

    // translate 1d indices for each variable into absolute lookup index
    igrid_tmp    = get_1DINDEX( ID, NVAR, &igrid_var[0]);
    igrid_1D     = gridmap->INVMAP[igrid_tmp] ;        

    for  ( ifun=0; ifun < NFUN; ifun++ )   {  
      FUNVAL[ifun]       = gridmap->FUNVAL[ifun][igrid_1D];
      WGT_SUM[ifun]     += (CORNER_WGT * FUNVAL[ifun]) ;
    }
    
    if ( LDMP ) {
      printf(" xxx CORNER_[WGT,SUM](%d)=%6.4f,%6.4f  FUNVAL=%f  WGT_SUM=%f\n",
	     icorner, CORNER_WGT, CORNER_WGTSUM, FUNVAL[0], WGT_SUM[0] );
      fflush(stdout);
    }
    
  } // corner

  if ( CORNER_WGTSUM <= 0.0 ) {
    sprintf(c1err,"Could not compute CORNER_WGT for gridmap ID=%d", 
	    gridmap->ID );
    sprintf(c2err,"%s", "data = ") ;
    for ( ivar=0; ivar < NVAR; ivar++ ) 
      { sprintf(c2err,"%s %f", c2err, *(data+ivar) ) ; }
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  for  ( ifun=0; ifun < NFUN; ifun++ )   {  
    *(interpFun+ifun) = WGT_SUM[ifun]   / CORNER_WGTSUM ; 
  }


  if ( LDMP ) {
    printf("  xxxx data=%6.2f , %6.2f  interpFun=%f (WGT_SUM=%f,%f)\n",
	   data[0], data[1], interpFun[0], WGT_SUM[0], CORNER_WGTSUM );
    printf("  xxxx DUMP DONE. \n");
    fflush(stdout);
    
  }

  return(SUCCESS) ;

} // end of interp_GRIDMAP


// ================================================
int  get_1DINDEX(int ID, int NDIM, int *indx ) {

  // Created April 2011 (initial use for hostlib weight-map)
  //
  // Return 1d index for NDIM-dimensional grid.
  // *indx is an array of indices for each dimension.
  // Each *indx value 0 to N-1 
  //
  // Note: must call init_1DINDEX(ID ...) once per ID
  // before calling this function.
  // 
  // If the number of grid-points in each dimension is
  // N1, N2 ... N_NDIM, then the returned index is an
  // integer from 0 to N1*N2* ... N_NDIM-1.
  // For example, for a 3 dimensional array [4][5][4],
  // the returned index is from 0 to 4*5*4-1 = 80-1 = 79
  //
  // Feb 25, 2013: ABORT if *indx exceeds NPT
  // Feb 12, 2018: indx is 0 to N-1 (no longer 1-N)

  char fnam[] = "get_1DINDEX" ;
  int i, offset, INDEX_1D, index_1d, NPT ;

  //------------ BEGIN -------------

  //  printf(" xxxx %s called with ID = %d \n", fnam, ID) ;

  if ( NPT_PERDIM_1DINDEX[ID][0] == 0 ) {
    sprintf(c1err,"ID=%d  is not defined.", ID );
    sprintf(c2err,"%s", "Must first call init_1DINDEX()");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  INDEX_1D = 0;

  for ( i=0; i < NDIM; i++ ) {   
    offset   = OFFSET_1DINDEX[ID][i];
    index_1d =  indx[i] ;       // index in this dimension
    INDEX_1D += (index_1d ) * offset ; // global 1D index

    /*
    printf(" xxx %s: i=%d index_1d=%3d  INDEX_1D=%6d\n",
	   fnam, i, index_1d, INDEX_1D); fflush(stdout);
    */

    // make sure that index does not exceed NPT
    NPT =    NPT_PERDIM_1DINDEX[ID][i] ;
    if ( index_1d >= NPT ) {
      sprintf(c1err,"index_1d=%d exceeds NPT=%d (ID=%d)", 
	      index_1d, NPT, ID );
      sprintf(c2err,"for idim = %d of %d", i, NDIM) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

  }
  return INDEX_1D ;

} // end of get_1DINDEX


void init_1DINDEX(int ID, int NDIM, int *NPT_PERDIM ) {

  // Apr 2011
  // init offsets needed to quickly compute 1d index
  // for multi-dimensional array or grid.
  //
  // ID        = reference ID for this mapping
  // NDIM      = number of dimensions
  // *NPT_PERDIM  = max number of elements in each dimension 
  //
  //
  // Feb 25 2013: store NPT_PERDIM
  // Feb 12 2018: refactor with indices starting at zero
  //

  int LDMP = 0 ;
  int i, NPT, NPT_LAST, OFFSET_LAST, OFFSET ;
  char fnam[] = "init_1DINDEX" ;

  // --------- BEGIN ----------

  if ( ID < 1 || ID >= MXMAP_1DINDEX ) {
    sprintf(c1err,"Invalid ID=%d", ID );
    sprintf(c2err,"Valid ID range is %d to %d", 1, MXMAP_1DINDEX-1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( NDIM < 0 || NDIM >= MXDIM_1DINDEX ) {
    sprintf(c1err,"Invalid NDIM=%d", NDIM );
    sprintf(c2err,"Valid NDIM range is %d to %d", 1, MXDIM_1DINDEX-1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  for ( i=0; i < NDIM; i++ ) {

    NPT_PERDIM_1DINDEX[ID][i]  = NPT_PERDIM[i]; 
    NPT                        = NPT_PERDIM[i]; 
    OFFSET_1DINDEX[ID][i] = OFFSET  = 0 ;

    if ( i > 0 ) {
      NPT_LAST    = NPT_PERDIM[i - 1];  
      OFFSET_LAST = OFFSET_1DINDEX[ID][i-1] ; 
      OFFSET      = OFFSET_LAST * NPT_LAST ;
      OFFSET_1DINDEX[ID][i] = OFFSET;
    }
    else {
      OFFSET_1DINDEX[ID][i] = OFFSET  = 1 ;
    }
      
    
    if ( LDMP ) {
      printf(" xxxx OFFSET_1DINDEX[ID=%d][ivar=%2d] = %7d   "
	     " NPT_PERDIM=%d\n",
	     ID, i, OFFSET, NPT_PERDIM[i] );
    }

  } // end NDIM


} // end of init_1DINDEX

void clear_1DINDEX(int ID) {
  //  printf("  Clear 1DINDEX for ID=%d \n", ID  );
  OFFSET_1DINDEX[ID][0] = 0;
  NPT_PERDIM_1DINDEX[ID][0] = 0 ;
}


// mangled functions for fortran
void clear_1dindex__(int *ID){
  clear_1DINDEX(*ID);
}
void init_1dindex__(int *ID, int *NDIM, int *NPT_PERDIM ) {
  init_1DINDEX(*ID, *NDIM, NPT_PERDIM);
}
int get_1dindex__(int *ID, int *NDIM, int *indx ) {
  return  get_1DINDEX(*ID, *NDIM, indx );
}
  
