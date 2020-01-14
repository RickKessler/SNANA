/**********************************************

  July 20 2019 R.Kessler
  Complete overhaul: only reads GRID in MAG space; no more code.
  Read GRID made by D.Jones, based on Burns 2018.

********************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fitsio.h"
#include "sntools.h"
#include "sntools_grid.h" 
#include "sntools_genSmear.h" // Aug 30 2019
#include "genmag_snoopy.h" 


// ============ MANGLED FORTRAN FUNCTIONS =============

int init_genmag_snoopy__( char *version, int *optmask, char *filtlist) {
  int istat;
  istat = init_genmag_snoopy ( version, *optmask, filtlist) ;
  return istat;
}


int genmag_snoopy__(int *ifilt, double *stretch, int *nobs, 
		  double *Trest, double *rest_mags, double *rest_magerrs ) {
  int istat;
  istat = genmag_snoopy(*ifilt, *stretch, *nobs, 
			Trest, rest_mags, rest_magerrs ) ;
  return istat;
}

void get_lamrange_snoopy__(double *lammin, double *lammax) {
  get_LAMRANGE_snoopy(lammin,lammax);
}


// =========== CODE ================

/* Initiallize reading rset-frame GRID for SNOOPY model
 *    Inputs:  version  (where templates.dat and templates/ folder reside)
 *             optmask   bitmask of options (unused)
 *
 *    Outputs: filtlist:  "BVugriYJH"
 */
int init_genmag_snoopy( char *VERSION, int optmask, char *filtlist) {
  /* Initialization function for the SNooPy templates */
  int  opt_fits_read ;
  char msg[100], tmpFile[400], version[60];
  char fnam[] = "init_genmag_snoopy";

  // ----------- BEGIN -------------

  // hard-wire a few strings
  sprintf(SNOOPYNAME,"snoopy");
  sprintf(FILTLIST_SNOOPY, "BVugriYJH" );
  sprintf(filtlist, "%s", FILTLIST_SNOOPY);

  sprintf(msg,"%s for %s", fnam, FILTLIST_SNOOPY );
  print_banner(msg);

  // construct SNOOPY_MODELPATH
  extract_MODELNAME(VERSION,
		    SNOOPY_MODELPATH, version); // returned


  if ( strlen(SNOOPY_MODELPATH) == 0 ) {
    // default location under $SNDATA_ROOT
    sprintf( SNOOPY_MODELPATH, "%s/models/snoopy/%s", 
             getenv("SNDATA_ROOT"), version );
  }


  // read snoopy.info file to get name of GRID file, etc ...
  parse_snoopy_modelinfo(SNOOPY_MODELPATH);

  sprintf(tmpFile,"%s/%s", SNOOPY_MODELPATH, SNOOPY_MODELINFO.GRIDFILE ) ;

  opt_fits_read = 1;  // bit1 => verbose 
  fits_read_SNGRID(opt_fits_read, tmpFile, &SNGRID_SNOOPY );  
  
  if ( strcmp(filtlist,SNGRID_SNOOPY.FILTERS) != 0 ) {
    sprintf(c1err,"Expected snoopy filters = '%s'.", 
	    filtlist );
    sprintf(c2err,"But found filters '%s' in snoopy grid.", 
	    SNGRID_SNOOPY.FILTERS );
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err);
  }
  
  return(SUCCESS);

} // end of init_genmag_snoopy



/* Generate SNooPy templates and return them
 *
 *  Inputs:   
 *     ifilt       (1 <= ifilt <= 9)  the filter you want.
 *     stretch      The shape/stretch parameter
 *     nobs         Number of rest-frame epochs
 *     Trest_list   rest-frame epochs (0 at peak) at which to
 *                    generate mags.
 *
 *  Outputs:  
 *    mag_list     rest-frame magnitudes 
 *    magerr_list  mag errors
*/
int genmag_snoopy(int ifilt, double stretch, int nobs, 
		  double *Trest_list, double *mag_list, double *magerr_list ) {

   /* These are the default window function paramters.  */

  int  o ;
  char fnam[] = "genmag_snoopy" ;

   // ----------- BEGIN -------------
  
  for(o=0; o < nobs; o++ ) { mag_list[o] = 0.0;  magerr_list[o] = 1.0; }

   /* check inputs */
   if (ifilt < IFILT_SNOOPY_MIN || ifilt > IFILT_SNOOPY_MAX ) {
     sprintf(c1err,"Invalid ifilt=%d", ifilt );
     sprintf(c2err,"Valid ifilt range is %d - %d ", 
	     IFILT_SNOOPY_MIN,  IFILT_SNOOPY_MAX ) ;
     errmsg(SEV_FATAL, 0, fnam, c1err,c2err);
   }


   gridinterp_snoopy(ifilt, stretch, nobs, Trest_list, 
		     mag_list, magerr_list );    // <== returned

   // Aug 30 2019 notes: if smear model in simulation (see istat_genSmear function)
   //
   // get_genSmear(double Trest_list, int NFILTER, double *Lam,
   //		     double *magSmear) 
   //
   //  mag_list[i] += magSmear[i] for each band i
   //
   // NOT used for LC fitting.

   return(SUCCESS);

} // end of genmag_snoopy





// ************************************************
void parse_snoopy_modelinfo(char *modelPath ) {

  char infoFile[MXPATHLEN], c_get[80]; 
  FILE *fp;
  char fnam[] = "parse_snoopy_modelinfo" ;

  // ------- BEGIN ----------

  sprintf(infoFile,"%s/%s.info", modelPath, SNOOPYNAME );

  if (! (fp = fopen(infoFile, "rt"))) { 
    sprintf(c1err,"Could not find info file:" ); 
    sprintf(c2err," '%s'", infoFile ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
   
   sprintf(SNOOPY_MODELINFO.GRIDFILE,"%s", "BLANK");
   SNOOPY_MODELINFO.OPT_GRIDFILE = 0;

   while( fscanf(fp, "%s", c_get) != EOF) {

     if ( strcmp(c_get,"GRIDFILE:") == 0 ) {
       SNOOPY_MODELINFO.OPT_GRIDFILE = 1 ;
       fscanf(fp, "%s", SNOOPY_MODELINFO.GRIDFILE );
     }


   }
   fclose(fp);

} // end of parse_snoopy_modelinfo


void get_LAMRANGE_snoopy(double *lammin, double *lammax) {
  *lammin = MINLAM_SNOOPY ;
  *lammax = MAXLAM_SNOOPY ;
} // end of get_LAMRANGE_snoopy


// *****************************************************************
void gridinterp_snoopy(int ifilt, double shape, 
		       int nobs, double *Trest_list,
		       double *mag_list, double *magerr_list ) {


  // Interpolate mag and magerr on snoopy grid.
  // Note that fits_read_SNGRID() must be called before
  // calling this function.

  int  i, ioff, im, index_shape, index_Trest ;
  int  ILC, IPTRLC_OFF[2], IPTRLC, NBIN_TREST, NBIN_SHAPE ;    

  double  
    Trest, ratio_Trest, ratio_shape
    ,gridFlux2[2][2]     // [Trest][shape]
    ,gridFlux1[2]
    ,gridErr2[2][2]     // [Trest][shape]
    ,gridErr1[2]
    ,dif
    ,TREST_MIN,  SHAPE_MIN 
    ,TREST_MAX,  SHAPE_MAX
    ,TREST_BIN,  SHAPE_BIN
    ;

  short I2TMP, I2ERR ;

  char fnam[] = "genmag_snoopy";

  // ----------- BEGIN -------------

  NBIN_TREST    = SNGRID_SNOOPY.NBIN[IPAR_GRIDGEN_TREST] ;
  NBIN_SHAPE    = SNGRID_SNOOPY.NBIN[IPAR_GRIDGEN_SHAPEPAR] ;

  TREST_BIN  = SNGRID_SNOOPY.BINSIZE[IPAR_GRIDGEN_TREST] ;
  SHAPE_BIN  = SNGRID_SNOOPY.BINSIZE[IPAR_GRIDGEN_SHAPEPAR] ;

  // get shape grid-index

  index_shape = INDEX_GRIDGEN(IPAR_GRIDGEN_SHAPEPAR, shape, &SNGRID_SNOOPY );
  if ( index_shape >= NBIN_SHAPE ) 
    { index_shape-- ; }

  ILC           = index_shape; 
  IPTRLC_OFF[0] = SNGRID_SNOOPY.PTR_GRIDGEN_LC[ILC+0] ; 
  IPTRLC_OFF[1] = SNGRID_SNOOPY.PTR_GRIDGEN_LC[ILC+1] ;

  SHAPE_MIN   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_SHAPEPAR][index_shape] ;
  SHAPE_MAX   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_SHAPEPAR][index_shape+1] ;
  ratio_shape = (shape -  SHAPE_MIN)/SHAPE_BIN ;

  // make sure that 1st word is BEGIN-LC marker
  I2TMP = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC_OFF[0]];
  if ( I2TMP != MARK_GRIDGEN_LCBEGIN ){
    sprintf(c1err,"First I*2 word of ILC=%d is %d .", ILC, I2TMP );
    sprintf(c2err,"But expected %d", MARK_GRIDGEN_LCBEGIN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // make sure that 2nd word is first 8 bits of ILC
  I2TMP = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC_OFF[0]+1];
  if ( I2TMP != ( ILC & 127 ) ){
    sprintf(c1err,"2nd I*2=%d  for ILC=%d, PTRLC_OFF=%d .", 
	    I2TMP, ILC, IPTRLC_OFF[0]+1  );
    sprintf(c2err,"But expected ILC&127 = %d", (ILC & 127) );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  ioff       = (ifilt-1)*NBIN_TREST + NPADWD_LCBEGIN-1 ;

  for ( i=0; i < nobs; i++ ) {

    Trest       = Trest_list[i];
    index_Trest = INDEX_GRIDGEN(IPAR_GRIDGEN_TREST,Trest, &SNGRID_SNOOPY );
    if ( index_Trest >= NBIN_TREST ) 
      { index_Trest-- ; }

    TREST_MIN   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_TREST][index_Trest] ;
    TREST_MAX   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_TREST][index_Trest+1] ;
    ratio_Trest = (Trest - TREST_MIN)/TREST_BIN ;

    /*
    if ( ratio_Trest < 0.0 || ratio_Trest > 1.000000001 ) {
      sprintf(c1err, "Bad ratio_Trest = %f at Trest(%s)=%6.2f  dm15=%6.2f", 
	      ratio_Trest, SNOOPY_TEMPLATE[ifilt].FILTERNAME, Trest, dm15 );
      sprintf(c2err,"index_Trest=%d  TREST(MIN,MAX,BIN)=%6.2f,%6.2f,%6.2f", 
	      index_Trest, TREST_MIN, TREST_MAX, TREST_BIN );
      MADABORT(c1err,c2err);
    }
    */

    // get gridFlux2 at the 4 corners bounding Trest,dm15
    for ( im=0; im <=1; im++ ) {    // dm15 edges
      IPTRLC      = IPTRLC_OFF[im] + ioff + index_Trest;
      I2TMP       = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC] ;
      I2ERR       = SNGRID_SNOOPY.I2GRIDGEN_LCERR[IPTRLC] ;
      gridFlux2[0][im]  = (double)I2TMP / GRIDGEN_I2LCPACK ;
      gridErr2[0][im]   = (double)I2ERR / GRIDGEN_I2LCPACK ;

      IPTRLC      = IPTRLC_OFF[im] + ioff + (index_Trest+1) ;
      I2TMP       = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC] ;
      I2ERR       = SNGRID_SNOOPY.I2GRIDGEN_LCERR[IPTRLC] ;
      gridFlux2[1][im]  = (double)I2TMP / GRIDGEN_I2LCPACK ;
      gridErr2[1][im]   = (double)I2ERR / GRIDGEN_I2LCPACK ;

      dif = gridFlux2[1][im] - gridFlux2[0][im] ;
      gridFlux1[im] = gridFlux2[0][im] + (dif * ratio_Trest);

      dif = gridErr2[1][im] - gridErr2[0][im] ;
      gridErr1[im] = gridErr2[0][im] + (dif * ratio_Trest);
    }

    dif = gridFlux1[1] - gridFlux1[0] ;
    mag_list[i]  = gridFlux1[0] + (dif * ratio_shape);

    dif = gridErr1[1] - gridErr1[0] ;
    magerr_list[i]  = gridErr1[0] + (dif * ratio_shape);

  } // i-loop over observations
  

} // end of gridinterp_snoopy

