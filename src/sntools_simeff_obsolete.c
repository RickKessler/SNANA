/****************************************

  Created July 2021 [extracted from sntools.c]
  Obsolete MLCS-SIMEFF utilities used in SDSS-II cosmology analysis
   (https://arxiv.org/abs/0908.4274)
  These functions are kept in SNANA to enable replicating these
  old MLCS fits and sim, but this MLCS-based methodology is obsolete,
  and hasn't been used since 2009.

 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "sntools.h"
#include "sntools_simeff_obsolete.h"

// =================================================
double get_SIMEFFMAP(int OPTMASK, int NVAR, double *GRIDVALS) {

  // Created July 2011 by R.Kessler
  // return EFFICIENCY from interpolating map.
  // Must call init_SIMEFFMAP() before calling this function.
  // GRIDVALS are the actaul values for each variable
  // on the grid.
  // OPTMASK bit 1 (LSB) => return EFF/EFFMAX (default returns EFF)

  int  ivar, istat, iflag;
  int LDMP = 0 ;

  double 
    EFF
    ,tmpval
    ,TMPVAL
    ,TMPVAL_LIST[MXGENVAR_SIMEFFMAP]
    ,MINVAL_4LOG10 = 1.0E-12
    ,arg
    ;
  char fnam[] =  "get_SIMEFFMAP";

  //  --------------- BEGIN -------------

  if ( NVAR != SIMEFFMAP.NGENVAR || SIMEFFMAP.NGENVAR == 0 ) {
    sprintf(c1err,"Passed NVAR=%d", NVAR);
    sprintf(c2err,"but NGENVAR = %d", SIMEFFMAP.NGENVAR );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // check for LIN, LOG and INV scales

  if ( LDMP ) {
    printf(" 0. xxxxxxxx --------------------------------- \n");
  }

  TMPVAL = 0.0 ;
  for ( ivar=1; ivar <= SIMEFFMAP.NGENVAR; ivar++ ) {
    iflag  = SIMEFFMAP.IFLAGSCALE[ivar];
    tmpval = GRIDVALS[ivar-1];

    if ( iflag == 1 )        // LINEAR
      { TMPVAL = tmpval ; }  
    else if ( iflag == 10 ) {     // LOG10
      arg = tmpval ;
      if ( tmpval < MINVAL_4LOG10 )  { arg = MINVAL_4LOG10 ; }
      TMPVAL = log10(arg) ; 
    }  
    else if ( iflag == -1 )     // INVERSE
      { TMPVAL = 1./tmpval ; }  
    else {
      sprintf(c1err,"Invalid IFLAGSCALE = %d for %s = %f",
	      iflag, SIMEFFMAP.VARNAME[ivar], tmpval );
      sprintf(c2err,"%s", "Something is really screwed up.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // if tmpval is  outside the min/max bound,
    // then pull it to the edge.

    if ( TMPVAL < SIMEFFMAP.VARMIN[ivar] ) {
      TMPVAL = SIMEFFMAP.VARMIN[ivar] ; 
    }

    if ( TMPVAL > SIMEFFMAP.VARMAX[ivar] )  { 
      TMPVAL = SIMEFFMAP.VARMAX[ivar] ;  //  RANGE * 1.0E-6 ; 
    }

    TMPVAL_LIST[ivar-1] = TMPVAL ;

    if ( LDMP ) {
      printf(" 1. xxxxxxx ivar=%d : iflag=%d  VAL(orig,interp) = %le ,%le \n",
	     ivar, iflag, tmpval, TMPVAL );
      fflush(stdout); 
    }

  }

  // do the interpolation
  
  istat = interp_GRIDMAP(&SIMEFF_GRIDMAP, TMPVAL_LIST, &EFF ); // return EFF

  if ( istat != SUCCESS ) {
    print_preAbort_banner(fnam);    
    for ( ivar=1; ivar <= SIMEFFMAP.NGENVAR; ivar++ ) {
      printf("\t %-12s = %f\n", 
	     SIMEFFMAP.VARNAME[ivar], *(GRIDVALS+ivar-1) );
    }

    sprintf(c1err,"interp_GRIDMAP returned istat=%d", istat);
    sprintf(c2err,"%s", "--> could not interplate");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }



  if ( OPTMASK & 1 ) { EFF /= SIMEFFMAP.EFFMAX ; }


  return EFF;

} // end of get_SIMEFFMAP


// =================================================
int init_SIMEFFMAP(char *file, char *varnamesList) {

  // Created July 2011 by R.Kessler
  // initialize effic. map stored in input *file.
  // Return list of blank-seprated variable names in *varnamesList.
  // For example, if the variables are redshift (Z), extinction (AV)
  // and mlcs DELTA, then  varnamesList = "Z AV DELTA" .
  //

  int NFUN=1;
  int NGENVAR, NVAR_READ, NEFF_READ, NBINTOT, ivar, gzipFlag ;
  double TMPVAL[MXGENVAR_SIMEFFMAP];
  double EFF;

  FILE *fp;
  char 
    c_get[200]
    ,PATH_SIMEFF[MXPATHLEN]
    ,FILE[MXPATHLEN]
    ,*cptr
    ,fnam[] = "init_SIMEFFMAP" 
    ;

  // -------------- BEGIN --------------

  varnamesList[0] = 0 ;      // init output
  SIMEFFMAP.NGENVAR = 0 ;       // init global 
  SIMEFFMAP.EFFMAX  = 0.0 ;

  // construct official PATH_SIMEFF
  sprintf(PATH_SIMEFF,"%s/models/simeff", getenv("SNDATA_ROOT") );

  // check user's directory, then check PATH_SIMEFF
  fp = snana_openTextFile(1,PATH_SIMEFF, file, FILE, &gzipFlag );
  
  // abort if there is no SIMEFF file.
  if ( fp == NULL ) {
    sprintf(c1err,"%s", "Could not open SIMEFF file");
    sprintf(c2err,"%s", file );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  sprintf(BANNER,"%s", "init_SIMEFF: read efficiency map from ");
  print_banner(BANNER);
  printf("\t %s\n", FILE );

  NVAR_READ =  NEFF_READ =  NGENVAR   = 0 ;
  NBINTOT   = 1 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"NGENVAR:") == 0 )  {
      readint(fp, 1, &NGENVAR);
      SIMEFFMAP.NGENVAR = NGENVAR ;
      if ( NGENVAR >= MXGENVAR_SIMEFFMAP ) {
	sprintf(c1err,"NGENVAR=%d exceeds bound of", NGENVAR);
	sprintf(c2err,"MXGENVAR_SIMEFFMAP = %d", MXGENVAR_SIMEFFMAP );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    } // NGENVAR: key


    if ( strcmp(c_get,"GENVAR:") == 0 )  {
      NVAR_READ++ ;
      readchar(fp, SIMEFFMAP.VARSCALE[NVAR_READ] );
      readchar(fp, SIMEFFMAP.VARNAME[NVAR_READ] );
      readint (fp,   1, &SIMEFFMAP.NBIN[NVAR_READ] );
      readdouble(fp, 1, &SIMEFFMAP.VARMIN[NVAR_READ] );
      readdouble(fp, 1, &SIMEFFMAP.VARMAX[NVAR_READ] );

      NBINTOT *= SIMEFFMAP.NBIN[NVAR_READ] ;
      SIMEFFMAP.NBINTOT = NBINTOT;
      sprintf(varnamesList,"%s %s ", 
	      varnamesList, SIMEFFMAP.VARNAME[NVAR_READ] );

      cptr = SIMEFFMAP.VARSCALE[NVAR_READ];
      if ( strcmp(cptr,"LIN") == 0 ) 
	{ SIMEFFMAP.IFLAGSCALE[NVAR_READ] = 1; }
      else if ( strcmp(cptr,"LOG") == 0 ) 
	{ SIMEFFMAP.IFLAGSCALE[NVAR_READ] = 10; } 
      else if ( strcmp(cptr,"INV") == 0 ) 
	{ SIMEFFMAP.IFLAGSCALE[NVAR_READ] = -1; } 
      else {
	sprintf(c1err,"Invalid scale = '%s'", cptr);
	sprintf(c2err,"%s", "valid scales are LIN, LOG and  INV");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }	       

      printf("\t Found MAP variable %s %-8s : %2d bins from %7.3f to %7.3f\n"
	     ,SIMEFFMAP.VARSCALE[NVAR_READ]
	     ,SIMEFFMAP.VARNAME[NVAR_READ]
	     ,SIMEFFMAP.NBIN[NVAR_READ]
	     ,SIMEFFMAP.VARMIN[NVAR_READ]
	     ,SIMEFFMAP.VARMAX[NVAR_READ]
	     );

      // allocate temp memory when all of the GENVAR keys are read
      if ( NVAR_READ == NGENVAR ) { malloc_SIMEFFMAP(+1); }

    } // GENVAR: key


    if ( strcmp(c_get,"EFF:") == 0 )  {

      NEFF_READ++ ;
      readdouble(fp, NGENVAR+1, TMPVAL ) ;

      for ( ivar=0; ivar < NGENVAR; ivar++ ) {
	SIMEFFMAP.TMPVAL[ivar][NEFF_READ-1] = TMPVAL[ivar];
      }

      EFF = TMPVAL[NGENVAR];

      SIMEFFMAP.TMPEFF[NEFF_READ-1] = EFF ;
      if ( EFF > SIMEFFMAP.EFFMAX )
	{ SIMEFFMAP.EFFMAX  = EFF ; }

    }


  } // fscanf


  //  printf("\t EFF(MAX) = %6.4f \n", SIMEFFMAP.EFFMAX);

  // init multi-dimensional interpolation
  init_interp_GRIDMAP(IDGRIDMAP_SIMEFFMAP, "SIMEFF",
		      NBINTOT, NGENVAR, NFUN, 0,
		      SIMEFFMAP.TMPVAL, &SIMEFFMAP.TMPEFF,
		      &SIMEFF_GRIDMAP ); // <== return this struct

  // free temp memory
  malloc_SIMEFFMAP(-1);

  // make sure we read what was expected.
  if ( NVAR_READ != NGENVAR ) {
    sprintf(c1err, "Expected NGENVAR = %d", NGENVAR );
    sprintf(c2err, "but found %d  'GENVAR:' keys.", NVAR_READ);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( NEFF_READ != NBINTOT ) {
    sprintf(c1err, "Expected %d 'EFF:' keys", NBINTOT );
    sprintf(c2err, "but found %d  keys.", NEFF_READ);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return NGENVAR ;

} // end of init_SIMEFFMAP

// mangled routines for fortrans
int init_simeffmap__(char *file, char *varnamesList) {
  return init_SIMEFFMAP(file, varnamesList);
}
double get_simeffmap__(int *OPTMASK, int *NVAR, double *GRIDVALS) {
  return get_SIMEFFMAP(*OPTMASK, *NVAR, GRIDVALS);
}


// ===========================================
void malloc_SIMEFFMAP(int flag) {

  // Created July 2011 by R.Kessler
  // flag > 0 => allocate
  // flag < 0 => free 

  int  NVAR, NBINTOT, MEM, ivar, I8, I8p, I4 ;
  double XMEM;
  //  char fnam[] = "malloc_SIMEFFMAP" ;

  // --------------- BEGIN ---------------

  NVAR    = SIMEFFMAP.NGENVAR; 
  NBINTOT = SIMEFFMAP.NBINTOT ;

  MEM = (NVAR+1)*(NBINTOT) * 8;
  XMEM = 1.0E-6 * (double)MEM ;

  I8p = sizeof(double*);
  I8  = sizeof(double );
  I4  = sizeof(int);

  if ( flag > 0 ) {
    printf("\t Allocate %6.3f MB of temp SIMEFFMAP memory (NBIN=%d). \n", 
	   XMEM, NBINTOT );

    SIMEFFMAP.TMPEFF = (double  *)malloc(NBINTOT*I8); // EFF memory

    SIMEFFMAP.TMPVAL = (double **)malloc(NVAR*I8p);   // memory for each var
    for ( ivar=0; ivar < NVAR; ivar++ ) {
      SIMEFFMAP.TMPVAL[ivar] = (double *)malloc(NBINTOT*I8);
    }


  }
  else if ( flag < 0 ) {
    printf("\t Free temp SIMEFFMAP memory. \n" );
    free(SIMEFFMAP.TMPEFF);
    for ( ivar=0; ivar < NVAR; ivar++ ) {
      free(SIMEFFMAP.TMPVAL[ivar]);
    }
    free(SIMEFFMAP.TMPVAL) ;
  }

  return ;

} // end of malloc_SIMEFFMAP
