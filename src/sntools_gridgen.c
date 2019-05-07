/****************************************
 Created Oct 2010 by R.Kessler
 
 Tools to generate GRID of templates with snlc_sim.exe
 (see manual for more details).

 The light curve structure for each SN is
 
   pad1  = MARK_GRIDGEN_LCBEGIN
   pad2  = first 8 bits of ILC (absolute lc index)
   i2mag(filt1,ep1) = mag * MAGPACK_GRIDGEN 
   i2mag(filt1,ep2) 
   i2mag(filt1,ep3)
   ...
   i2mag(filt1,Nep)
   i2mag(filt2,ep1) 
   ... 
   i2mag(filt2,Nep)
   ...
   ...
   i2mag(Nfilt,Nep)
   pad3 = MARG_GRIDGEN_LCEND
   pad4 = MARG_GRIDGEN_LCEND

 After reading each LC with your external program, it is 
 strongly recommended to verify the pad words to make sure 
 that your fits reader is not lost. 

 As of this time (Oct 2010) there are no SNANA utilites
 to read the FITS files produced by these functions.
 
 With one logz bin and the snoopy model, write relative flux
 instead of mags.


     HISTORY


*********************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#include "sntools.h"      // general snana stuff
#include "fitsio.h"
#include "sntools_grid.h"
#include "snlc_sim.h"
#include "genmag_SIMSED.h"

// =============================
//     GRID  functions 
// =============================

// ************************************
void init_GRIDsource(int opt) {

  // Created Oct 16, 2010 by R.Kessler
  // Initialize for "GRID" option to generate light curves
  // on a fixed grid and then dump them out.
  //
  // opt=0,1 => two different init stages;
  //            latter stage is for filter-LamAvg and non1a
  //

  // ----------- BEGIN -----------

  if (opt == 0 ) {
    init0_GRIDsource();
    return ;
  }

  if (opt == 1 ) {
    init1_GRIDsource();
    return ;
  }



} // end of init_GRIDsource


// ******************************************
void init0_GRIDsource(void) {

  // Jun 2013
  // Separated from init_GRIDsource, and fixed to work with SIMSED model.
  //
  // Apr 28 2019:
  //  Fix to work with SIMSED model using INDEX param.

  float GENRANGE_LOCAL[NPAR_GRIDGEN+1][2];
  float PMIN, PMAX;

  int 
    IPAR, IPAR2, i,  NFILT, NBIN, NINDEX
    ,indx[NPAR_GRIDGEN+1]
    ,ilcoff, ilc, ilc_last, itmp
    ,iz, ish, ic1, ic2
    ,Nz, Nsh, Nc1, Nc2, Nf, Nep
    ,NPAR_STORE, iptroff, MEM
    ;

  char valname[20] ;
  char *NAME ;
  char fnam[] = "init0_GRIDsource";

  // --------------- BEGIN ----------

  NROW_WRITE_TOT = 0;

  // transer NBIN from GRIDGEN_INPUTS struct to SNGRID struct
  for ( IPAR = 1; IPAR <= NPAR_GRIDGEN; IPAR++ ) {
    SNGRID_WRITE.NBIN[IPAR] = GRIDGEN_INPUTS.NBIN[IPAR];
  }

  SIMLIB_INIT_DRIVER();

  SNGRID_WRITE.IVERSION = IVERSION_GRID_WRITE ; // Apr 2013

  sprintf(SNGRID_WRITE.SURVEY, "%s", GENLC.SURVEY_NAME );

  sprintf(SNGRID_WRITE.MODEL,  "%s", INPUTS.MODELNAME );

  sprintf(BANNER,"init_GRIDsource for SURVEY=%s  and  GENMODEL=%s : \n", 
	  SNGRID_WRITE.SURVEY, SNGRID_WRITE.MODEL );
  print_banner(BANNER);


  // load gen-ranges into local array GENRANGE_LOCAL

  IPAR = IPAR_GRIDGEN_LOGZ;
  GENRANGE_LOCAL[IPAR][0] = log10f(INPUTS.GENRANGE_REDSHIFT[0]);
  GENRANGE_LOCAL[IPAR][1] = log10f(INPUTS.GENRANGE_REDSHIFT[1]);
  sprintf(SNGRID_WRITE.NAME[IPAR], "%s", "LOGZ" );


  IPAR = IPAR_GRIDGEN_SHAPEPAR;

  if ( INDEX_GENMODEL == MODEL_NON1ASED ) {

    if ( INPUTS.NON1ASED.INDEX[1] > 0 ) 
      { NINDEX = INPUTS.NON1ASED.NINDEX ; }
    else  { 
      NINDEX = count_NON1A_LIST(INPUTS.NON1ASED.PATH);   // 5.2019
      printf("\t ALL-NON1A option -> count %d NON1A keys \n", NINDEX);
    }
    GENRANGE_LOCAL[IPAR][0] = 1.0 ;  
    GENRANGE_LOCAL[IPAR][1] = (float)NINDEX ;
    SNGRID_WRITE.NBIN[IPAR] = NINDEX ;
    sprintf(SNGRID_WRITE.NAME[IPAR], "%s", "NONIA" );  // sparse index
  }
  else if ( INDEX_GENMODEL == MODEL_SIMSED ) {
    itmp = INPUTS.IPAR_SIMSED_SHAPE ; 
    if ( itmp > 0 ) {
      // for SALT2 grid in Mosher 2014 
      PMIN = INPUTS.GENGAUSS_SIMSED[itmp].RANGE[0] ;
      PMAX = INPUTS.GENGAUSS_SIMSED[itmp].RANGE[1] ;
      NAME = INPUTS.PARNAME_SIMSED[itmp] ;	
    }
    else {      
      NBIN  = count_SIMSED_INFO(INPUTS.GENMODEL);
      //  printf(" xxx %s: NBIN -> %d \n", fnam, NBIN);
      PMIN  = 1.0 ;
      PMAX  = (float)NBIN;
      NAME  = INPUTS.PARNAME_SIMSED[0] ;  
    }

    GENRANGE_LOCAL[IPAR][0]  = PMIN ;
    GENRANGE_LOCAL[IPAR][1]  = PMAX ;
    SNGRID_WRITE.NBIN[IPAR]  = NBIN ;
    sprintf(SNGRID_WRITE.NAME[IPAR],  "%s", NAME );
  }
  else {
    // SNIa: mlcs, snoopy, ...
    GENRANGE_LOCAL[IPAR][0] = INPUTS.GENGAUSS_SHAPEPAR.RANGE[0] ;
    GENRANGE_LOCAL[IPAR][1] = INPUTS.GENGAUSS_SHAPEPAR.RANGE[1] ;
    sprintf(SNGRID_WRITE.NAME[IPAR], "%s", GENLC.SHAPEPAR_NAME);
  }



  IPAR  = IPAR_GRIDGEN_COLORPAR ;  // AV or color
  IPAR2 = IPAR_GRIDGEN_COLORLAW ;  // RV or beta
  if  ( INDEX_GENMODEL == MODEL_SALT2 ) {
    GENRANGE_LOCAL[IPAR][0]  = INPUTS.GENGAUSS_SALT2c.RANGE[0];
    GENRANGE_LOCAL[IPAR][1]  = INPUTS.GENGAUSS_SALT2c.RANGE[1];
    GENRANGE_LOCAL[IPAR2][0] = INPUTS.GENGAUSS_SALT2BETA.RANGE[0] ;
    GENRANGE_LOCAL[IPAR2][1] = INPUTS.GENGAUSS_SALT2BETA.RANGE[1] ;

    sprintf(SNGRID_WRITE.NAME[IPAR],  "%s", "c" );
    sprintf(SNGRID_WRITE.NAME[IPAR2], "%s", "BETA" );
  }
  else if ( INDEX_GENMODEL == MODEL_SIMSED ) {
    itmp = INPUTS.IPAR_SIMSED_COLOR ;
    if ( itmp > 0 ) {
      GENRANGE_LOCAL[IPAR][0]  = INPUTS.GENGAUSS_SIMSED[itmp].RANGE[0];
      GENRANGE_LOCAL[IPAR][1]  = INPUTS.GENGAUSS_SIMSED[itmp].RANGE[1];
      sprintf(SNGRID_WRITE.NAME[IPAR], "%s", INPUTS.PARNAME_SIMSED[itmp] );
    }
    else {
      GENRANGE_LOCAL[IPAR][0]  = 0.0 ;
      GENRANGE_LOCAL[IPAR][1]  = 0.0 ;
      sprintf(SNGRID_WRITE.NAME[IPAR],  "noColor" );
    }
     
    // color law is optional to generate GRID
    itmp = INPUTS.IPAR_SIMSED_CL ;
    if ( itmp > 0 ) {
      GENRANGE_LOCAL[IPAR2][0]  = INPUTS.GENGAUSS_SIMSED[itmp].RANGE[0];
      GENRANGE_LOCAL[IPAR2][1]  = INPUTS.GENGAUSS_SIMSED[itmp].RANGE[1];
      sprintf(SNGRID_WRITE.NAME[IPAR2],  "%s", INPUTS.PARNAME_SIMSED[itmp] );
    }
    else {
      GENRANGE_LOCAL[IPAR2][0]  = 0.0 ;
      GENRANGE_LOCAL[IPAR2][1]  = 0.0 ;
      sprintf(SNGRID_WRITE.NAME[IPAR2],  "noCL" );
    }

  }

  else {
    // mlcs, snoopy ...
    GENRANGE_LOCAL[IPAR][0]  = INPUTS.GENRANGE_AV[0];
    GENRANGE_LOCAL[IPAR][1]  = INPUTS.GENRANGE_AV[1];
    GENRANGE_LOCAL[IPAR2][0] = INPUTS.GENGAUSS_RV.RANGE[0] ;
    GENRANGE_LOCAL[IPAR2][1] = INPUTS.GENGAUSS_RV.RANGE[1] ;

    sprintf(SNGRID_WRITE.NAME[IPAR],  "%s", "AV" );
    sprintf(SNGRID_WRITE.NAME[IPAR2], "%s", "RV" );
  }


  NFILT = INPUTS.NFILTDEF_OBS;
  IPAR  = IPAR_GRIDGEN_FILTER;
  GENRANGE_LOCAL[IPAR][0] = 0.0 ;
  GENRANGE_LOCAL[IPAR][1] = (float)(NFILT-1);
  SNGRID_WRITE.NBIN[IPAR] = NFILT ;
  sprintf(SNGRID_WRITE.NAME[IPAR], "%s", "FILTER" );
  //  printf(" xxx NFILT = %d \n", NFILT); // xxxx remove this

  IPAR = IPAR_GRIDGEN_TREST;
  GENRANGE_LOCAL[IPAR][0] = INPUTS.GENRANGE_TREST[0];
  GENRANGE_LOCAL[IPAR][1] = INPUTS.GENRANGE_TREST[1];
  sprintf(SNGRID_WRITE.NAME[IPAR], "%s", "TREST" );


  Nz    = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_LOGZ] ;
  Nsh   = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_SHAPEPAR] ;
  Nc1   = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_COLORPAR] ;
  Nc2   = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_COLORLAW] ;
  Nf    = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_FILTER] ;
  Nep   = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_TREST] ;


  // make sure that snlc_sim  structures are big enough.
  if ( Nep*Nf >= MXEPSIM ) {
    sprintf(c1err,"Nep x Nfilt = %d x %d = %d exceeds bound of MXEPSIM=%d", 
	    Nep, Nf, Nep*Nf, MXEPSIM );
    sprintf(c2err,"Reduce Nep*Nfilt or increase MXEPSIM" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  SNGRID_WRITE.NGRIDGEN_LC    = Nz * Nsh * Nc1 * Nc2  ;
  INPUTS.NGEN_LC = SNGRID_WRITE.NGRIDGEN_LC ;
  INPUTS.NGEN    = SNGRID_WRITE.NGRIDGEN_LC ;
  set_screen_update(SNGRID_WRITE.NGRIDGEN_LC);

  NPAR_STORE = NPAR_GRIDGEN - 2 ; // store all but FILTER & TREST


  OPT_SNOOPY_FLUXPACK =  ( Nz == 1 &&  INDEX_GENMODEL == MODEL_SNOOPY ) ;
  if ( OPT_SNOOPY_FLUXPACK ) {
    GRIDGEN_I2LCPACK = FLUXPACK_GRIDGEN ;
    sprintf(valname,"%s", "relFlux");
  }
  else {
    GRIDGEN_I2LCPACK = MAGPACK_GRIDGEN;
    sprintf(valname,"%s", "mag");
  }
    
  // compute bin-size and grid value at each bin

  for ( IPAR = 1; IPAR <= NPAR_GRIDGEN; IPAR++ ) {

    NBIN = SNGRID_WRITE.NBIN[IPAR];
    NAME = SNGRID_WRITE.NAME[IPAR] ;
    if ( NBIN >= MXGRIDGEN ) {
      sprintf(c1err,"NBIN=%d exceeds bound for %s", NBIN, NAME);
      sprintf(c2err,"Check MXGRIDGEN = %d", MXGRIDGEN );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
						
    fillbins(1           // (I) binning option
	     ,NAME       // (I) 
	     ,NBIN       // (I)
	     ,GENRANGE_LOCAL[IPAR]          // (I)
	     ,&SNGRID_WRITE.BINSIZE[IPAR]   // (O)
	     ,&SNGRID_WRITE.VALUE[IPAR][1]  // (O)
	     );


    SNGRID_WRITE.VALMIN[IPAR] = SNGRID_WRITE.VALUE[IPAR][1] ;
    SNGRID_WRITE.VALMAX[IPAR] = SNGRID_WRITE.VALUE[IPAR][NBIN] ;

    SNGRID_WRITE.ILCOFF[IPAR] = 1;

    for ( IPAR2 = 1; IPAR2 <= IPAR-1 ; IPAR2++ ) {
      SNGRID_WRITE.ILCOFF[IPAR] *= SNGRID_WRITE.NBIN[IPAR2] ;
    }
      
    printf("\t ILCOFF(%s,IPAR=%d) = %d \n", 
	   NAME, IPAR, SNGRID_WRITE.ILCOFF[IPAR] ) ; fflush(stdout);

    if ( IPAR <= NPAR_STORE ) {
      MEM = 4*(SNGRID_WRITE.NGRIDGEN_LC+1) ;
      SNGRID_WRITE.PTR_VALUE[IPAR] = (int *)malloc(MEM) ;
    }

  }  // IPAR  loop

  
  SNGRID_WRITE.NGRIDGEN_PER_LC = NFILT * Nep ;
  SNGRID_WRITE.NWD_I2GRIDGEN   = 
    SNGRID_WRITE.NGRIDGEN_PER_LC + NPADWD_LCBEGIN + NPADWD_LCEND ;

  
  SNGRID_WRITE.I2GRIDGEN_LCMAG = 
    (short *)malloc(4*SNGRID_WRITE.NWD_I2GRIDGEN+4);
  SNGRID_WRITE.I2GRIDGEN_LCERR = 
    (short *)malloc(4*SNGRID_WRITE.NWD_I2GRIDGEN+4);
  SNGRID_WRITE.PTR_GRIDGEN_LC  = 
    (int   *)malloc(4*SNGRID_WRITE.NGRIDGEN_LC+8);

  // total size include 2 2-bytes words (MAG + ERR),
  // hence the factor of 4 for each LC point.
  SNGRID_WRITE.SIZEOF_GRIDGEN = 
    4 * SNGRID_WRITE.NWD_I2GRIDGEN * SNGRID_WRITE.NGRIDGEN_LC ;

  // load global comments
  for ( i=0 ; i < NCOMMENT_GRIDGEN; i++) 
    { sprintf(COMMENT_GRIDGEN[0],"%s", "" ); }


  i = -1;
  i++; sprintf(COMMENT_GRIDGEN[i],"SURVEY=%s  GENMODEL=%s", 	  
	  SNGRID_WRITE.SURVEY, SNGRID_WRITE.MODEL );

  i++; sprintf(COMMENT_GRIDGEN[i],"SNANA_DIR = %s", PATH_SNANA_DIR);	  

  i++; sprintf(COMMENT_GRIDGEN[i],
	       "Generate %d light curves, each with %s and %d epochs."
	       ,SNGRID_WRITE.NGRIDGEN_LC, INPUTS.GENFILTERS
	       ,SNGRID_WRITE.NBIN[IPAR_GRIDGEN_TREST] );
  printf("\t %s\n", COMMENT_GRIDGEN[i] );
  
  i++; sprintf(COMMENT_GRIDGEN[i],
	       "GRID storage (mag+err, excluding header): %7.3f MB.",
	       1.0E-6 * (float)SNGRID_WRITE.SIZEOF_GRIDGEN );
  printf("\t %s\n", COMMENT_GRIDGEN[i] );


  i++; sprintf(COMMENT_GRIDGEN[i],
	       "All %s and errors are multiplied by %d for I*2 storage.",
	       valname, (int)GRIDGEN_I2LCPACK );

  // for each set of IPAR indices, store them as a function of
  // absolute light curve index 'ilc'.   

  ilc_last = 0 ;

  for ( ish = 1; ish <= Nsh;  ish++ ) {   
    for ( ic2 = 1; ic2 <= Nc2;  ic2++ ) {
      for ( ic1 = 1; ic1 <= Nc1;  ic1++ ) {
	for ( iz = 1; iz <= Nz;  iz++ ) {

	  indx[IPAR_GRIDGEN_LOGZ]     = iz ;   // redshift
	  indx[IPAR_GRIDGEN_COLORPAR] = ic1 ;  // color
	  indx[IPAR_GRIDGEN_COLORLAW] = ic2 ;  // RV/beta
	  indx[IPAR_GRIDGEN_SHAPEPAR] = ish ;  // shape or NON1A index

	  ilc = 1 ; 

	  for ( IPAR = 1; IPAR <= NPAR_STORE; IPAR++ ) {
	    ilcoff = SNGRID_WRITE.ILCOFF[IPAR] ;
	    ilc   += ( ilcoff * (indx[IPAR]-1) ) ;

	    if ( ilc > SNGRID_WRITE.NGRIDGEN_LC || ilc < 1 ) {
	      sprintf(c1err,"Invalid ilc=%d ", ilc );
	      sprintf(c2err,"at iz=%d ishape=%d ic1=%d ic2=%d",
		      iz, ish, ic1, ic2 );
	      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	    }
	  } // end of IPAR

	  SNGRID_WRITE.PTR_GRIDGEN_LC[ilc] = 
	    1 + SNGRID_WRITE.NWD_I2GRIDGEN * (ilc-1);

	  for ( IPAR = 1; IPAR <= NPAR_STORE; IPAR++ ) {
	    iptroff = (IPAR-1) * SNGRID_WRITE.NGRIDGEN_LC ;
	    SNGRID_WRITE.PTR_VALUE[IPAR][ilc] = indx[IPAR] ;
	  }

	  if ( ilc != ilc_last+1 ) {
	    sprintf(c1err,"ilc=%d but ilc_last=%d", ilc, ilc_last);
	    sprintf(c2err,"iz=%2d ic1=%2d ic2=%d ish=%d \n",
		 iz, ic1, ic2, ish );
	    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 	    
	  }

	  ilc_last = ilc ;

	  /*
	  printf(" xxx ilc(iz=%2d ic1=%2d ic2=%d ish=%d  ) = %4d \n",
	  iz, ic1, ic2, ish, ilc ); */

	} // ilum
      } // iz
    } // ic1
  } // ic2

} // end of init0_GRIDsource(void) {


// ******************************************
void init1_GRIDsource(void) {

  int   ifilt, ifilt_obs, NFILT, NBIN, i, i2 ;
  float LAM; 
  char fnam[] = "init1_GRIDsource" ;

  // -------------- BEGIN --------------

  NFILT = INPUTS.NFILTDEF_OBS;
  for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
    ifilt_obs  = INPUTS.IFILTMAP_OBS[ifilt];
    LAM        = INPUTS.LAMAVG_OBS[ifilt_obs] ;
    SNGRID_WRITE.FILTER_LAMAVG[ifilt] = LAM ; 
  }
  
  
  // transfer NON1A info here
  // xxx mark delete  if ( SNTYPE_GRIDGEN() == SNTYPE_GRIDGEN_NONIa ) {
  if (INDEX_GENMODEL == MODEL_NON1ASED  )  {
    NBIN = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_SHAPEPAR] ;
    for ( i = 1; i <= NBIN; i++ ) { 
      i2  = INPUTS.NON1ASED.INDEX[i] ; // i = sparse index; i2=absolute
      SNGRID_WRITE.NON1A_INDEX[i] = i2 ;

      sprintf(SNGRID_WRITE.NON1A_NAME[i], "%s", 
	      INPUTS.NON1ASED.LIST_NAME[i2] );

      sprintf(SNGRID_WRITE.NON1A_CTYPE[i], "%s", 
	      INPUTS.NON1ASED.LIST_TYPE[i2] );

      SNGRID_WRITE.NON1A_WGT[i]        = INPUTS.NON1ASED.WGT[i] ;
      SNGRID_WRITE.NON1A_MAGOFF[i]     = INPUTS.NON1ASED.MAGOFF[i] ;
      SNGRID_WRITE.NON1A_MAGSMEAR[i]   = INPUTS.NON1ASED.MAGSMEAR[i] ;
      SNGRID_WRITE.NON1A_ITYPE_USER[i] = INPUTS.NON1ASED.SNTAG[i] ;  
    }

  }


  if (INDEX_GENMODEL == MODEL_SIMSED  ) {
    //    SNGRID_WRITE.NBIN[IPAR_GRIDGEN_SHAPEPAR] = SEDMODEL.NSURFACE ;
  }


  return ;

} // end of init1_GRIDsource(void) {


// ******************************************
void gen_GRIDevent(int ilc) {

  // called from snlc_sim.
  // Jun 2013; fixed to work with SIMSED model

  int indx, IPAR, ifilt, NFILT, ifilt_obs ;
  int iep, NEP, NEPTOT, NEWMJD, itmp    ;

  double cpar, claw, Trest, Tobs, z1, logz, MU8, ARG8 ;

  char fnam[] = "gen_GRIDevent";

  // --------- BEGIN ----------

  /* xxxxx
  IPAR = IPAR_GRIDGEN_LOGZ ;
  for(iz=0; iz <= SNGRID_WRITE.NBIN[IPAR]; iz++ ) {
    printf(" xxx logz(iz=%3d) = %.5f\n", iz, SNGRID_WRITE.VALUE[IPAR][iz] ); 
    fflush(stdout);
  }
  xxxxxx */


  GENLC.CID  = ilc ; 

  IPAR = IPAR_GRIDGEN_LOGZ ; 
  indx = SNGRID_WRITE.PTR_VALUE[IPAR][ilc] ;  // iz index
  logz = SNGRID_WRITE.VALUE[IPAR][indx] ;   
  GENLC.REDSHIFT_CMB = pow(10.0,logz);

  gen_distanceMag(GENLC.REDSHIFT_CMB, GENLC.REDSHIFT_CMB,
		  &GENLC.DLMU, &GENLC.LENSDMU );

  z1 = 1.0 + GENLC.REDSHIFT_CMB ;
  SNGRID_WRITE.CURRENT_VALUE[IPAR] = logz ;
  GENLC.REDSHIFT_HELIO = GENLC.REDSHIFT_CMB ; // no coords to translate

  IPAR = IPAR_GRIDGEN_SHAPEPAR ; 
  indx = SNGRID_WRITE.PTR_VALUE[IPAR][ilc] ;   // ilum index
  GENLC.SHAPEPAR = SNGRID_WRITE.VALUE[IPAR][indx] ;
  SNGRID_WRITE.CURRENT_VALUE[IPAR] = GENLC.SHAPEPAR ;
  GENLC.DELTA = GENLC.DM15 
    = GENLC.STRETCH = GENLC.SALT2x1 = GENLC.SHAPEPAR;


  IPAR = IPAR_GRIDGEN_COLORPAR ;               // color/AV
  indx = SNGRID_WRITE.PTR_VALUE[IPAR][ilc] ;
  cpar = SNGRID_WRITE.VALUE[IPAR][indx] ;
  SNGRID_WRITE.CURRENT_VALUE[IPAR] = cpar ;

  IPAR = IPAR_GRIDGEN_COLORLAW ;              // Beta/RV
  indx = SNGRID_WRITE.PTR_VALUE[IPAR][ilc] ;
  claw = SNGRID_WRITE.VALUE[IPAR][indx] ;
  SNGRID_WRITE.CURRENT_VALUE[IPAR] = claw ;

  if ( INDEX_GENMODEL == MODEL_SALT2 ) {
    GENLC.SALT2c     = cpar ;
    GENLC.SALT2beta  = claw ;
    GENLC.SALT2alpha = INPUTS.GENGAUSS_SALT2ALPHA.PEAK ;
    MU8 = GENLC.DLMU ;
    GENLC.SALT2x0 = SALT2x0calc(GENLC.SALT2alpha, GENLC.SALT2beta, 
				GENLC.SALT2x1, GENLC.SALT2c, MU8 );
  }
  else if ( INDEX_GENMODEL == MODEL_SIMSED ) {  

    GENLC.SIMSED_PARVAL[0]  = GENLC.SHAPEPAR ; 

    itmp = INPUTS.IPAR_SIMSED_SHAPE ;    
    if ( itmp >= 0 ) { GENLC.SIMSED_PARVAL[itmp]  = GENLC.SHAPEPAR ; }

    itmp = INPUTS.IPAR_SIMSED_COLOR ;
    if ( itmp >= 0 ) { GENLC.SIMSED_PARVAL[itmp]  = cpar ; }

    itmp = INPUTS.IPAR_SIMSED_CL ;    
    if ( itmp >= 0 ) { GENLC.SIMSED_PARVAL[itmp]  = claw ; }

    // set flux scale for distance modulus.
    MU8  =  GENLC.DLMU ;
    ARG8 = -0.4 * MU8 ;
    GENLC.SALT2x0 = pow(TEN, ARG8 );
  }
  else {
    GENLC.AV = cpar ;
    GENLC.RV = claw ;
  }


  // now fill Trest array that SIMLIB_read() would normally fill.

  GENLC.MWEBV        =  0.0;
  GENLC.MWEBV_SMEAR  =  0.0;  // fixed Mar 14 2016
  GENLC.PEAKMJD      =  0.0 ;
  GENLC.NFILTDEF_OBS = INPUTS.NFILTDEF_OBS ;
  SIMLIB_OBS_GEN.PIXSIZE[0]  = -9.0 ;

  NEP    = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_TREST] ; // NEP per filter
  NEPTOT = NEWMJD = 0;
  NFILT  = INPUTS.NFILTDEF_OBS ;

  for ( iep = 1; iep <= NEP; iep++ ) {

    NEWMJD++; 
    GENLC.NEWMJD = NEWMJD;

    GENLC.EPOCH_RANGE_NEWMJD[NEWMJD][0] = NEPTOT + 1 ;
    GENLC.EPOCH_RANGE_NEWMJD[NEWMJD][1] = NEPTOT + NFILT ;

    Trest  = SNGRID_WRITE.VALUE[IPAR_GRIDGEN_TREST][iep] ;
    Tobs   = Trest*z1;

    for ( ifilt = 0; ifilt < NFILT; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];  

      NEPTOT++ ;
      GENLC.NEPOCH = NEPTOT;

      if ( NEPTOT >= MXEPSIM) {
	sprintf(c1err,"NEPTOT = %d exceeds bound of MXEPSIM", NEPTOT);
	sprintf(c2err,"Occured at ifilt_obs=%d and iep=%d (Trest=%6.2f)",
		ifilt_obs, iep, Trest );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      GENFILT.Trest8[ifilt_obs][NEPTOT] = Trest ;
      GENLC.MJD[NEPTOT]         = 53000.0 + Tobs ;
      GENLC.epoch8_obs[NEPTOT]  = Tobs  ;
      GENLC.epoch8_rest[NEPTOT] = Trest ;
      GENLC.IFILT_OBS[NEPTOT]   = ifilt_obs;

      /*
      printf("\t xxx ifilt=%d  ifilt_obs=%2d (%c) Tobs = %6.1f \n",
	     ifilt, ifilt_obs, FILTERSTRING[ifilt_obs], Tobs );
      */

    } // ifilt
  } // iep

  //   debugexit("gen_GRIDevent");

  return ;

} // end of gen_GRIDevent

// ****************************************
void wr_GRIDfile(int OPT, char *GRIDGEN_FILE) {

  // OPT=1: open  file and write header info
  // OPT=2: append with next SN
  // OPT=3: end of job checks (Jun 2012)

  int istat;
  char fnam[] = "wr_GRIDfile" ;
  char clobberFile[200];

  // ----------- BEGIN ----------

  if ( OPT == 1 ) {

    // open file and write header
    load_EXTNAME_GRIDGEN();

    if ( strcmp(GRIDGEN_INPUTS.FORMAT,"TEXT") == 0 ) 
      { OPT_GRIDGEN_FORMAT = OPT_GRIDGEN_FORMAT_TEXT ; } // human-readable
    else if ( strcmp(GRIDGEN_INPUTS.FORMAT,"FITS") == 0 ) 
      { OPT_GRIDGEN_FORMAT = OPT_GRIDGEN_FORMAT_FITS ; } // for psnid 
    else {
      sprintf(c1err, "Invalid GRID_FORMAT option: '%s'", 
	      GRIDGEN_INPUTS.FORMAT);
      sprintf(c2err, "%s", "Check your sim-input file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
    }


    // open GRIDGEN file for output.

    //    sprintf(GRIDGEN_FILE , "%s", GRIDfile); // transfter to global

    // check TEXT vs. FITS format.
    if ( OPT_GRIDGEN_FORMAT == OPT_GRIDGEN_FORMAT_TEXT ) {
      if ( (fp_GRIDGEN_TEXT = fopen(GRIDGEN_FILE, "wt")) == NULL ) {  
	sprintf(c1err,"%s", "Could not open GRIDGEN output file");
	sprintf(c2err,"%s", GRIDGEN_FILE );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      printf("\n  Opened GRIDGEN text file: \n\t %s \n", GRIDGEN_FILE );
      wrhead_GRIDfile_text();
    }
    else if ( OPT_GRIDGEN_FORMAT == OPT_GRIDGEN_FORMAT_FITS ) {

      sprintf(clobberFile, "!%s", GRIDGEN_FILE );
      istat = 0;
      fits_create_file(&fp_GRIDGEN_FITS, clobberFile, &istat) ;
      check_fitserror("fits_create_file", istat) ;
      printf("\n  Opened GRIDGEN fits file (istat=%d): \n\t %s \n", 
	     istat, GRIDGEN_FILE );
      wrhead_GRIDfile_fits();
    }

  }
  else if ( OPT == 2 ) {
    
    //    append with next SN
    if ( OPT_GRIDGEN_FORMAT == OPT_GRIDGEN_FORMAT_TEXT ) 
      { append_GRIDfile_text(); }
    else if ( OPT_GRIDGEN_FORMAT == OPT_GRIDGEN_FORMAT_FITS ) 
      { append_GRIDfile_fits(); }

  }
  else if ( OPT == 3 ) {
    end_GRIDfile();
  }

} // end of wr_GRIDfile


// ****************************************
void wrhead_GRIDfile_fits(void) {

  // Write header info and GRID tables 
  //
  // Dec 13, 2010: 
  //  - move supplemental nonIa table to end so that ordering
  //    of other tables is unchanges
  //  - DELTA-GRID and x1-GRID -> LUMI-GRID  (generic name)
  //  - AV-GRID and c-GRID     -> COLOR-GRID (generic name)
  //  - RV-GRID and BETA-GRID  -> RV/BETA-GRID 
  //
  // Feb 21, 2011: 
  //  fix dumb bug by creating I2LCMAG table last so that it's
  //  not over-written by the NONIa-INFO table.
  //  Add another PAR_NAME column to SNPAR_INFO table ;
  //  now we have both the TBL_NAME and PAR_NAME.
  //
  // Apr 2 2013: call get_GRIDKEY and write KEY 
  //             Write LAMAVG as 2nd column in FILTER table.
  //
  //
  // Apr 25 2013: fix aweful casting bug that shows up on
  //               64 bit machines, but not on 32 bit ?!?!?!
  //                  NBIN: TSHORT -> TUSHORT
  //                  Define ILCOFF as unsigned int instead of unsigned long.
  //
  // Aug 12 2016: 
  //   write TYPE as C16 instead of C8 in case PEC1A-xxx is too long
  // ----------------

  long
     NAXIS = 1
    ,NAXES = 0
    ,FPIX  = 1
    ,ZERO  = 0
    ,NROW  = 0
    ;

  int
    IVERSION
    ,colnum, firstrow, firstelem, nrow, nbin, istat 
    ,ipar, i
    ,ncol = 1
    ;
  
  char 
     blank[8] = " " 
    ,stringF4[4]  = "1E"   // float 
    //    ,stringU2[4]  = "1U"   // 16-bit unsigned short (TUSHORT)
    ,stringI2[4]  = "1I"   // 16-bit signed short
    ,stringU4[4]  = "1V"   // 32-bit unsigned int
    ,stringC8[4]  = "8A"   // ascii 
    ,stringC16[4] = "16A"   // ascii 
    ,stringC20[4] = "20A"    // ascii
    ,stringC120[8] = "120A"  // long ascii 
    ,extname[20]  = " "
    ,ctmp[100]
    ,*ttype[20]
    ,*tform[20]
    ,*tunit[20]
    ;
  
  char *tparcol[] = 
    { "TBL_NAME" , "PAR_NAME", "NBIN" , "VALMIN" , "VALMAX" , "BINSIZE" , "ILCOFF" } ;

  char *tparmag[] = { "I2MAG" , "I2MAGERR" } ;
  char *tptr[]    = { "PTR_LCMAG" } ;

  // for NONIa-info table
  char  *tparnon1a[] = { "GRID-INDEX" , "SNANA-INDEX", "TYPE", "NAME", 
			 "WGT", "MAGOFF" , "MAGSMEAR", "ITYPE_USER", } ;
  short  I2NONIA_GRID[MXNON1A_TYPE];
  short  I2NONIA_SNANA[MXNON1A_TYPE];
  char  *TBLPAR_NON1A_CTYPE[MXNON1A_TYPE] ;
  char  *TBLPAR_NON1A_NAME[MXNON1A_TYPE] ;
  float  TBLPAR_NON1A_WGT[MXNON1A_TYPE] ;
  float  TBLPAR_NON1A_MAGOFF[MXNON1A_TYPE] ;
  float  TBLPAR_NON1A_MAGSMEAR[MXNON1A_TYPE] ;
  short  I2NONIA_ITYPE_USER[MXNON1A_TYPE]; // added Jan 2017

  float dumarray[2], ftmpval, *fptr ;

  // define arrays for par-info table
  char             *TBLPAR_EXTNAME[NPAR_GRIDGEN+1] ;
  char             *TBLPAR_PHYSNAME[NPAR_GRIDGEN+1] ;
  unsigned short    TBLPAR_NBIN[NPAR_GRIDGEN+1] ;
  float             TBLPAR_VALMIN[NPAR_GRIDGEN+1] ;
  float             TBLPAR_VALMAX[NPAR_GRIDGEN+1] ;
  float             TBLPAR_BINSIZE[NPAR_GRIDGEN+1] ;
  unsigned int      TBLPAR_ILCOFF[NPAR_GRIDGEN+1] ;

  //  char     fnam[]  = "wrhead_GRIDfile_fits"  ;

  // -------------- BEGIN --------------

  IVERSION = SNGRID_WRITE.IVERSION ;

  istat = 0;
  for ( ipar=0; ipar < 20; ipar++ )  tunit[ipar] = blank;

  // create mandatory primary image (length=0)
  fits_create_img(fp_GRIDGEN_FITS, FLOAT_IMG, NAXIS, &NAXES, &istat) ;
  check_fitserror("Create zero-len primary fits image", istat) ;

  // -------------------------------
  printf("\t create GRIDGEN header keys. \n");


  fits_update_key(fp_GRIDGEN_FITS, TINT, "IVERSION", 
		  &IVERSION, "Internal GRID version", &istat );

  get_GRIDKEY();  // unique grid identifier
  fits_update_key(fp_GRIDGEN_FITS, TSTRING, "GRIDKEY", 
		  SNGRID_WRITE.UNIQUE_KEY, "Unique GRID key (time in UTC)", &istat );

  fits_update_key(fp_GRIDGEN_FITS, TSTRING, "SURVEY", 
		  SNGRID_WRITE.SURVEY, "Simulated Survey", &istat );

  fits_update_key(fp_GRIDGEN_FITS, TSTRING, "GENMODEL", 
		  SNGRID_WRITE.MODEL, "Simulated model", &istat );

  fits_update_key(fp_GRIDGEN_FITS, TSTRING, "FILTERS", 
		  INPUTS.GENFILTERS, "Simulated filters", &istat );

  sprintf(ctmp,"Total number of %s light curves", INPUTS.GENFILTERS );
  fits_update_key(fp_GRIDGEN_FITS, TINT, "NTOT_LC", 
		  &SNGRID_WRITE.NGRIDGEN_LC, ctmp, &istat );

  ftmpval = 1.0E-6 * (float)SNGRID_WRITE.SIZEOF_GRIDGEN ;
  fits_update_key(fp_GRIDGEN_FITS, TFLOAT, "SIZEI2LC", 
		  &ftmpval, "Memory (MB) to store I2LCMAG GRID", &istat );

  fits_update_key(fp_GRIDGEN_FITS, TSTRING, "SNANA", 
		  PATH_SNANA_DIR, "SNANA path", &istat );

  // ----------
  if ( OPT_SNOOPY_FLUXPACK ) {
    sprintf(ctmp,"%s", "relFLux = I2Flux / I2LCPACK" );
  }
  else {
    sprintf(ctmp,"%s", "mag = I2mag / I2LCPACK" );
  }
  fits_update_key(fp_GRIDGEN_FITS, TFLOAT, "I2LCPACK", 
		  &GRIDGEN_I2LCPACK, ctmp, &istat);

  fits_update_key(fp_GRIDGEN_FITS, TFLOAT, "MAGPACK", 
		  &GRIDGEN_I2LCPACK, "same as I2LCPACK", &istat);
  // ---------

  ftmpval = NULLFLOAT ;
  fits_update_key(fp_GRIDGEN_FITS, TFLOAT, "NULLMAG", 
		  &ftmpval, "value for undefined model mag", &istat);

  ftmpval = INPUTS.OMEGA_MATTER;
  fits_update_key(fp_GRIDGEN_FITS, TFLOAT, "Omega_M", 
		  &ftmpval, "Simulated Dark Matter content.", &istat);

  ftmpval = INPUTS.OMEGA_LAMBDA;
  fits_update_key(fp_GRIDGEN_FITS, TFLOAT, "Omega_DE", 
		  &ftmpval, "Simulated Dark Energy content.", &istat);

  ftmpval = INPUTS.W0_LAMBDA;
  fits_update_key(fp_GRIDGEN_FITS, TFLOAT, "w0", 
		  &ftmpval, "Simulated w for DE.", &istat);

  ftmpval = INPUTS.H0;
  fits_update_key(fp_GRIDGEN_FITS, TFLOAT, "H0", 
		  &ftmpval, "Simulated Hubble param at z=0.", &istat);

  // write zero-length image 
  fits_write_img(fp_GRIDGEN_FITS, TFLOAT, FPIX, ZERO, dumarray, &istat);
  check_fitserror("Write zero-len primary fits image", istat) ;

  // ------------------------------------------------

  printf("\t create GRIDGEN info-table for each parameter : \n" );

  NROW     =  0 ;
  ncol     =  7 ; 

  tform[0] = stringC20 ; // name of parameter
  tform[1] = stringC8  ; // name of parameter
  tform[2] = stringI2  ; // NBIN
  tform[3] = stringF4  ; // MIN value
  tform[4] = stringF4  ; // MAX value
  tform[5] = stringF4  ; // binsize
  tform[6] = stringU4  ; // ILCOFF

  sprintf(extname, "%s", EXTNAME_SNPAR_INFO );
  fits_create_tbl(fp_GRIDGEN_FITS, BINARY_TBL, NROW, ncol,
		  tparcol, tform, tunit, extname, &istat );
  sprintf(BANNER,"fits_create_tbl for %s", extname );

  check_fitserror(BANNER, istat) ;


  // fill table values
  for ( ipar = 1; ipar <= NPAR_GRIDGEN; ipar++ ) {

    // EXTNAME is the generic table name such as SHAPEPAR
    // SNGRID_WRITE.NAME is the specific name such as DELTA or DM15

    TBLPAR_EXTNAME[ipar]   = EXTNAME_GRIDGEN[ipar]      ; // table name
    TBLPAR_PHYSNAME[ipar]  = SNGRID_WRITE.NAME[ipar]    ; // physical name
    TBLPAR_NBIN[ipar]      = SNGRID_WRITE.NBIN[ipar]    ;
    TBLPAR_VALMIN[ipar]    = SNGRID_WRITE.VALMIN[ipar]  ;
    TBLPAR_VALMAX[ipar]    = SNGRID_WRITE.VALMAX[ipar]  ;
    TBLPAR_BINSIZE[ipar]   = SNGRID_WRITE.BINSIZE[ipar] ;
    TBLPAR_ILCOFF[ipar]    = SNGRID_WRITE.ILCOFF[ipar] ;
  }


  firstrow=1;  firstelem=1;  nrow = NPAR_GRIDGEN;

  colnum = 1;
  fits_write_col(fp_GRIDGEN_FITS, TSTRING, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_EXTNAME[1], &istat);  
  sprintf(BANNER,"fits_write_col for grid-par column: %s", "TBL-EXTNAME" );
  check_fitserror(BANNER, istat) ;

  colnum = 2 ;
  fits_write_col(fp_GRIDGEN_FITS, TSTRING, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_PHYSNAME[1], &istat); 
  sprintf(BANNER,"fits_write_col for grid-par column: %s", "PHYS-NAME" );
  check_fitserror(BANNER, istat) ;

  colnum = 3;
  fits_write_col(fp_GRIDGEN_FITS, TUSHORT, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_NBIN[1], &istat);
  sprintf(BANNER,"fits_write_col for grid-par column: %s", "NBIN" );
  check_fitserror(BANNER, istat) ;

  colnum = 4;
  fits_write_col(fp_GRIDGEN_FITS, TFLOAT, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_VALMIN[1], &istat);
  sprintf(BANNER,"fits_write_col for grid-par column: %s", "VALMIN" );
  check_fitserror(BANNER, istat) ;

  colnum = 5;
  fits_write_col(fp_GRIDGEN_FITS, TFLOAT, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_VALMAX[1], &istat);
  sprintf(BANNER,"fits_write_col for grid-par column: %s", "VALMAX" );
  check_fitserror(BANNER, istat) ;
 
  colnum = 6;
  fits_write_col(fp_GRIDGEN_FITS, TFLOAT, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_BINSIZE[1], &istat);
  sprintf(BANNER,"fits_write_col for grid-par column: %s", "BINSIZE" );
  check_fitserror(BANNER, istat) ;

  colnum = 7;
  fits_write_col(fp_GRIDGEN_FITS, TUINT, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_ILCOFF[1], &istat);
  sprintf(BANNER,"fits_write_col for grid-par column: %s", "ILCOFF" );
  check_fitserror(BANNER, istat) ;


  // ------------------------------------------------
  char NAME_LAMAVG[] = "LAMAVG"; 

  printf("\t create GRIDGEN table of values for : " );
  NROW  =  0;  tform[0] =  stringF4 ;

  for ( ipar = 1; ipar <= NPAR_GRIDGEN; ipar++ ) {
    nbin     =  SNGRID_WRITE.NBIN[ipar] ;
    fptr     = &SNGRID_WRITE.VALUE[ipar][1];
    ttype[0] =  SNGRID_WRITE.NAME[ipar] ;
    printf("%s ", ttype[0] );
    sprintf(extname, "%s", EXTNAME_GRIDGEN[ipar] );

    ncol = 1;
    if ( ipar == IPAR_GRIDGEN_FILTER ) 
      { ncol = 2 ;  ttype[1] = NAME_LAMAVG ; tform[1] = stringF4 ; }
    
    fits_create_tbl(fp_GRIDGEN_FITS, BINARY_TBL, NROW, ncol,
		    ttype, tform, tunit, extname, &istat );
    sprintf(BANNER,"fits_create_tbl for %s", ttype[0] );
    check_fitserror(BANNER, istat) ;

    colnum=1;  firstrow=1;  firstelem=1;  nrow = nbin;
    fits_write_col(fp_GRIDGEN_FITS, TFLOAT, colnum, firstrow, firstelem, nrow,
		   fptr, &istat);
    sprintf(BANNER,"fits_write_tbl for %s", ttype[0] );
    check_fitserror(BANNER, istat) ;

    
    // for filter table, include LAMAVG
    if ( ipar == IPAR_GRIDGEN_FILTER ) {
      colnum = 2 ; 
      fptr   = &SNGRID_WRITE.FILTER_LAMAVG[0] ;
      fits_write_col(fp_GRIDGEN_FITS, TFLOAT, colnum, firstrow, firstelem, 
		     nrow, fptr, &istat);
      sprintf(BANNER,"fits_write_tbl for %s", ttype[1] );
      check_fitserror(BANNER, istat) ;
    }
    

  }   // end of ipar loop
  printf("\n" );


  // -------------------------------
  // write NONIa-INFO  table

  if ( SNTYPE_GRIDGEN() == SNTYPE_GRIDGEN_NONIa ) {
    printf("\t create GRIDGEN table of nonIa types : \n" );
    fflush(stdout);

    if ( IVERSION <=2 ) 
      { ncol = 4; }
    else if ( IVERSION == 3 ) 
      { ncol = 7 ; }
    else if ( IVERSION >= 4 ) 
      { ncol = 8 ; }

    NROW = 0;    ipar = IPAR_GRIDGEN_SHAPEPAR ; 
    firstrow=1;  firstelem=1;  nrow = SNGRID_WRITE.NBIN[ipar] ;

    tform[0] =  stringI2  ;  // GRID-index (same as shapepar)
    tform[1] =  stringI2  ;  // nonIa SNANA index (same as shapepar)
    tform[2] =  stringC16 ;  // type (II, Ib, ...) 

    tform[3] =  stringC120 ;  // allow long name with path 
    tform[4] =  stringF4  ; // WGT    (Aug 30 2013)
    tform[5] =  stringF4  ; // MAGOFF (Aug 30 2013)
    tform[6] =  stringF4  ; // MAGSMEAR (Aug 30 2013)
    tform[7] =  stringI2  ;  // ITYPE_USER (Jan 2017)

    fptr     = &SNGRID_WRITE.VALUE[ipar][1];
    for ( i=1; i <= nrow; i++ ) {
      I2NONIA_GRID[i]          = (int)SNGRID_WRITE.VALUE[ipar][i];
      I2NONIA_SNANA[i]         = SNGRID_WRITE.NON1A_INDEX[i] ;
      TBLPAR_NON1A_CTYPE[i]    = SNGRID_WRITE.NON1A_CTYPE[i] ;
      TBLPAR_NON1A_NAME[i]     = SNGRID_WRITE.NON1A_NAME[i] ;
      TBLPAR_NON1A_WGT[i]      = SNGRID_WRITE.NON1A_WGT[i] ;
      TBLPAR_NON1A_MAGOFF[i]   = SNGRID_WRITE.NON1A_MAGOFF[i] ;
      TBLPAR_NON1A_MAGSMEAR[i] = SNGRID_WRITE.NON1A_MAGSMEAR[i] ;
      I2NONIA_ITYPE_USER[i]    = SNGRID_WRITE.NON1A_ITYPE_USER[i] ;
    }

    sprintf(extname, "%s", EXTNAME_NONIa_INFO );
    fits_create_tbl(fp_GRIDGEN_FITS, BINARY_TBL, NROW, ncol,
		    tparnon1a, tform, tunit, extname, &istat );
    sprintf(BANNER,"fits_create_tbl for %s", extname );
    check_fitserror(BANNER, istat) ;

    colnum = 1;
    fits_write_col(fp_GRIDGEN_FITS, TSHORT, colnum, firstrow, firstelem, nrow,
		 &I2NONIA_GRID[1], &istat);  
    sprintf(BANNER,"fits_write_col for  %s", "NONIA GRID-INDEX" );
    check_fitserror(BANNER, istat) ;

    colnum = 2;
    fits_write_col(fp_GRIDGEN_FITS, TSHORT, colnum, firstrow, firstelem, nrow,
		 &I2NONIA_SNANA[1], &istat);  
    sprintf(BANNER,"fits_write_col for  %s", "NONIA SNANA-INDEX" );
    check_fitserror(BANNER, istat) ;

    colnum = 3;
    fits_write_col(fp_GRIDGEN_FITS, TSTRING, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_NON1A_CTYPE[1], &istat);  
    sprintf(BANNER,"fits_write_col for  %s", "NONIA-TYPE" );
    check_fitserror(BANNER, istat) ;

    colnum = 4;
    fits_write_col(fp_GRIDGEN_FITS, TSTRING, colnum, firstrow, firstelem, nrow,
		 &TBLPAR_NON1A_NAME[1], &istat);  
    sprintf(BANNER,"fits_write_col for %s", "NONIA-NAME" );
    check_fitserror(BANNER, istat) ;

    // add WGT and MAGOFF, Aug 30 2013
    if ( SNGRID_WRITE.IVERSION >= 3 ) {
      colnum = 5;
      fits_write_col(fp_GRIDGEN_FITS, TFLOAT, colnum, firstrow, firstelem,
		     nrow, &TBLPAR_NON1A_WGT[1], &istat);  
      sprintf(BANNER,"fits_write_col for %s", "NONIA-WGT" );
      check_fitserror(BANNER, istat) ;

      colnum = 6;
      fits_write_col(fp_GRIDGEN_FITS, TFLOAT, colnum, firstrow, firstelem,
		     nrow, &TBLPAR_NON1A_MAGOFF[1], &istat);  
      sprintf(BANNER,"fits_write_col for %s", "NONIA-MAGOFF" );
      check_fitserror(BANNER, istat) ;

      colnum = 7 ;
      fits_write_col(fp_GRIDGEN_FITS, TFLOAT, colnum, firstrow, firstelem,
		     nrow, &TBLPAR_NON1A_MAGSMEAR[1], &istat);  
      sprintf(BANNER,"fits_write_col for %s", "NONIA-MAGSMEAR" );
      check_fitserror(BANNER, istat) ;

      colnum = 8 ; // added Jan 2017
      fits_write_col(fp_GRIDGEN_FITS, TSHORT, colnum, firstrow, firstelem,
		     nrow, &I2NONIA_ITYPE_USER[1], &istat );
      sprintf(BANNER,"fits_write_col for %s", "NONIA-ITYPE_USER" );
      check_fitserror(BANNER,istat) ; 
    }

  } // end of nonIa if-block

  // -------------------------------
  // write pointer for each ILC index


  printf("\t create GRIDGEN table of pointers to %d light curves \n", 
	 SNGRID_WRITE.NGRIDGEN_LC );
  fflush(stdout);

  ncol = 1 ;  NROW = 0 ;  tform[0] = stringU4 ; 

  sprintf(extname, "%s", EXTNAME_PTRI2LCMAG );
  fits_create_tbl(fp_GRIDGEN_FITS, BINARY_TBL, NROW, ncol,
		  tptr, tform, tunit, extname, &istat );
  sprintf(BANNER,"fits_create_tbl for %s", "PTR_I2LCMAG" );
  check_fitserror(BANNER, istat) ;

  colnum=1;  firstrow=1;  firstelem=1;  nrow = SNGRID_WRITE.NGRIDGEN_LC ;
  fits_write_col(fp_GRIDGEN_FITS, TUINT, colnum, firstrow, firstelem, nrow,
		 &SNGRID_WRITE.PTR_GRIDGEN_LC[1], &istat);
  sprintf(BANNER,"fits_write_tbl for %s", extname );
  check_fitserror(BANNER, istat) ;


  // ---------------------------------------------------
  // finally, create the LC table with mag and magerr
  // This must be the last table created since it is filled later.

  printf("\t create GRIDGEN table for model MAGS and ERRORS \n" );

  NROW     =  0 ; 
  ncol     =  2 ;
  tform[0] = stringI2 ; 
  tform[1] = stringI2 ;
  sprintf(extname, "%s", EXTNAME_I2LCMAG );

  fits_create_tbl(fp_GRIDGEN_FITS, BINARY_TBL, NROW, ncol,
		  tparmag, tform, tunit, extname, &istat );
  sprintf(BANNER,"fits_create_tbl for %s", "I*2 MAGS & ERRORS" );
  check_fitserror(BANNER, istat) ;


} // end of wrhead_GRIDfile_fits



// ********************************
void get_GRIDKEY(void) {

  // xxxx  char vers_snana[20], vers_photom[200], ctkey[100] ;
  char ctkey[100] ;
  time_t tkey ;
  struct tm * ptm;

  /* xxxxxx mark delete Dec 10 2017 xxxxxxxxxx
  get_snana_versions__(vers_snana, vers_photom, 20, 200);
  xxx */

  time(&tkey);
  ptm = gmtime ( &tkey );

  
  sprintf(ctkey,"%4.4d-%2.2d-%2.2d_%2.2d:%2.2d:%2.2d", 
	  1900+ptm->tm_year, ptm->tm_mon, ptm->tm_mday ,
	  ptm->tm_hour, ptm->tm_min, ptm->tm_sec );
	  
  sprintf(SNGRID_WRITE.UNIQUE_KEY, "%s--%s--%s--%s", 
	  SNANA_VERSION_CURRENT, SNGRID_WRITE.SURVEY, 
	  SNGRID_WRITE.MODEL, ctkey );

  
} // end of get_GRIDKEY

// ****************************************
void wrhead_GRIDfile_text(void) {

  int IPAR, NBIN, i  ;
  char  *name;
  float VAL ;

  // ---------- BEGIN ---------


  for ( i = 0; i < NCOMMENT_GRIDGEN; i++ )
    { fprintf(fp_GRIDGEN_TEXT, "\t %s\n", COMMENT_GRIDGEN[i] ); }


  for ( IPAR = 1; IPAR <= NPAR_GRIDGEN; IPAR++ ) {

    NBIN = SNGRID_WRITE.NBIN[IPAR] ;
    name = SNGRID_WRITE.NAME[IPAR] ;

    fprintf(fp_GRIDGEN_TEXT, " %2d GRID-BINS for %8s : ", 
	    NBIN, name );

    for ( i = 1; i <= NBIN; i++ ) {
      VAL = SNGRID_WRITE.VALUE[IPAR][i] ;
      fprintf(fp_GRIDGEN_TEXT, "%7.3f ", VAL); 
    }
    fprintf(fp_GRIDGEN_TEXT, "\n");

    // for LOGZ, also print list of Z
    if ( strcmp(name,"LOGZ") == 0 ) {
      fprintf(fp_GRIDGEN_TEXT, "\t\t      => Z : " );
      for ( i = 1; i <= NBIN; i++ ) {
	VAL = powf(10.0,SNGRID_WRITE.VALUE[IPAR][i]);
	fprintf(fp_GRIDGEN_TEXT, "%6.4f ", VAL);
      }
      fprintf(fp_GRIDGEN_TEXT, "\n");
    }

    // for NON1A, print list of NON1a
    if ( strcmp(name,"NONIA") == 0 || strcmp(name,"nonIa") == 0	) {

      // print absolute index vs. sparse index
      fprintf(fp_GRIDGEN_TEXT, "\t    => NON1A INDEX : " );
      for ( i = 1; i <= NBIN; i++ )
	{ fprintf(fp_GRIDGEN_TEXT, "%6d ", SNGRID_WRITE.NON1A_INDEX[i] ); }
      fprintf(fp_GRIDGEN_TEXT, "\n");

      // print type vs. sparse index
      fprintf(fp_GRIDGEN_TEXT, "\t    => NON1A TYPES : " );
      for ( i = 1; i <= NBIN; i++ )
	{ fprintf(fp_GRIDGEN_TEXT, "%6s ", SNGRID_WRITE.NON1A_CTYPE[i] ); }
      fprintf(fp_GRIDGEN_TEXT, "\n");

      // print full name vs. sparse index
      fprintf(fp_GRIDGEN_TEXT, "\t    => NON1A NAMES : " );
      for ( i = 1; i <= NBIN; i++ )
	{ fprintf(fp_GRIDGEN_TEXT, "%6s ", SNGRID_WRITE.NON1A_NAME[i] ); }
      fprintf(fp_GRIDGEN_TEXT, "\n");
    }

  } // IPAR

    fprintf(fp_GRIDGEN_TEXT, "\n");

} // end of wrhead_GRIDfile_text


// ****************************************
void update_GRIDarrays(void) {

  // Oct 2010 R.Kessler
  // Update GRIDGEN arrays that get written out in append_GRIDfile().
  //

  int NFILT, ifilt, ifilt_obs, NEP, ep, jindx, NTMP ; 
  int I2MAGVAL, I2MAGERR, JINDX,  NEPFILT[MXFILTINDX] ;
  long int ilc;
  double magval, magerr, Trest ;
  char cfilt[2];
  char fnam[] = "update_GRIDarrays";

  //------- BEGIN ----------


  ilc    = GENLC.CID;
  NFILT  = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_FILTER] ; 
  NEP    = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_TREST] ; // NEP per filter

  // start with LC-begin markers
  SNGRID_WRITE.I2GRIDGEN_LCMAG[0] = MARK_GRIDGEN_LCBEGIN;
  SNGRID_WRITE.I2GRIDGEN_LCMAG[1] = ( ilc & 127 );
  SNGRID_WRITE.I2GRIDGEN_LCERR[0] = MARK_GRIDGEN_LCBEGIN;
  SNGRID_WRITE.I2GRIDGEN_LCERR[1] = ( ilc & 127 );

  // tranfer light curve mags to I2GRIDGEN_LCMAG array

  for ( ifilt = 0; ifilt < NFILT ; ifilt++ )
    { NEPFILT[ifilt] = 0; }

  for ( ep = 1; ep <= GENLC.NEPOCH; ep++ ) {   
    ifilt_obs = GENLC.IFILT_OBS[ep] ;
    ifilt     = GENLC.IFILTINVMAP_OBS[ifilt_obs];
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );

    magval    = GENLC.genmag8_obs[ep] ;
    magerr    = GENLC.generr8_obs[ep] ;
    magerr   *= INPUTS.GENMODEL_ERRSCALE ;
    Trest     = GENLC.epoch8_obs[ep]/(1. + GENLC.REDSHIFT_CMB) ;

    /*
    if ( fabs(Trest) < 0.1 && ifilt_obs == 5 ) {
      printf(" xxx %s: write magval(z,T=%.1f,ep=%d) = %.3f  z=%.3f\n", 
	     fnam, Trest, ep, magval, GENLC.REDSHIFT_CMB ); fflush(stdout);
      fflush(stdout);
    }
    */

    // for snoopy  model and just 1 logz-bin, use rest-frame
    // mag & error which is really just the relative flux (0-1)
    if ( OPT_SNOOPY_FLUXPACK ) {
      magval = GENLC.genmag8_rest[ep] ;
      magerr = GENLC.generr8_rest[ep] ;
    }
    
    if ( magval > MAXMAG_GRIDGEN ) { magval = MAXMAG_GRIDGEN ; }
    if ( magerr > MAXMAG_GRIDGEN ) { magerr = MAXMAG_GRIDGEN ; }

    // convert mag and error for I*2 format

    magval *= GRIDGEN_I2LCPACK;  if ( magval > 0.0) { magval += 0.5 ; }
    magerr *= GRIDGEN_I2LCPACK;

    I2MAGVAL  = (int)(magval) ;
    I2MAGERR  = (int)(magerr) ;

    NEPFILT[ifilt]++ ;
    NTMP  = NEPFILT[ifilt] ; 
    jindx = ifilt*NEP + NTMP-1 ;

    // make sure indices are withing bounds
    if ( NTMP > NEP ) {
      sprintf(c1err,"NEPFILT[ifilt=%d(%s)]=%d exceeds %d bins/filter", 
	      ifilt, cfilt, NTMP, NEP );
      errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
    }
    if ( jindx >= SNGRID_WRITE.NGRIDGEN_PER_LC || jindx < 0 ) {
      sprintf(c1err,"jindx = %d exceeds  size of I2GRIDGEN_LCMAG", 
	      jindx );
      sprintf(c2err, "NEPFILT[ifilt=%d(%s)]=%d  NEP=%d", 
	      ifilt, cfilt, NTMP, NEP) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    // load I*2 array for FITS-storage
    JINDX = jindx + NPADWD_LCBEGIN;
    SNGRID_WRITE.I2GRIDGEN_LCMAG[JINDX] = I2MAGVAL;
    SNGRID_WRITE.I2GRIDGEN_LCERR[JINDX] = I2MAGERR;

  } // ep


  // tack on  end-of-lc markers
  JINDX = NPADWD_LCBEGIN + SNGRID_WRITE.NGRIDGEN_PER_LC - 1 ;
  JINDX++ ; 
  SNGRID_WRITE.I2GRIDGEN_LCMAG[JINDX] = MARK_GRIDGEN_LCEND ;
  SNGRID_WRITE.I2GRIDGEN_LCERR[JINDX] = MARK_GRIDGEN_LCEND ;

  JINDX++ ; 
  SNGRID_WRITE.I2GRIDGEN_LCMAG[JINDX] = MARK_GRIDGEN_LCEND ;
  SNGRID_WRITE.I2GRIDGEN_LCERR[JINDX] = MARK_GRIDGEN_LCEND ;

  // check that NEPFILT[ifilt] = NEP for all filters.

  for ( ifilt = 0; ifilt < NFILT; ifilt++ ) {

    NTMP = NEPFILT[ifilt];
    if ( NTMP != NEP ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];  
      sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );

      sprintf(c1err,"Invalid NEPFILT[ifilt=%d(%s)] = %d", 
	      ifilt, cfilt, NTMP );
      sprintf(c2err,"Should be equal to NEP = %d ", NEP );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  }

} // end of update_GRIDarrays


// ******************************
void append_GRIDfile_fits(void) {

  int ilc, colnum, firstrow, firstelem, nrow, istat ;
  //  char fnam[] = "append_GRIDfile_fits";

  // -------------- BEGIN --------------

  istat     = 0;
  ilc       = GENLC.CID;
  firstelem = 1 ;  
  firstrow  = SNGRID_WRITE.PTR_GRIDGEN_LC[ilc];
  nrow      = SNGRID_WRITE.NWD_I2GRIDGEN ;

  NROW_WRITE_TOT += nrow ;

  colnum    = 1 ;  
  fits_write_col(fp_GRIDGEN_FITS, TSHORT, colnum, firstrow, firstelem, nrow,
		 SNGRID_WRITE.I2GRIDGEN_LCMAG, &istat);
  sprintf(BANNER,"fits_write_tbl for I2MAG and ilc = %d", ilc );
  check_fitserror(BANNER, istat) ;


  colnum    = 2 ;  
  fits_write_col(fp_GRIDGEN_FITS, TSHORT, colnum, firstrow, firstelem, nrow,
		 SNGRID_WRITE.I2GRIDGEN_LCERR, &istat);
  sprintf(BANNER,"fits_write_tbl for I2ERR and ilc = %d", ilc );
  check_fitserror(BANNER, istat) ;


}  // append_GRIDfile_fits

// ******************************
void append_GRIDfile_text(void) {

  // Oct 16,  2010 R.Kessler
  // write SN grid-lightcurve in human-readable format ...
  // Just for visual inspection and not meant to be parsed 
  // by other programs.

  int 
    IPAR, NFILT, ifilt, ifilt_obs, NEP, ep, jindx
    ,i2mag, i2err, iwd
    ;

  double Tobs, Trest ;
  char cfilt[2];
  //  char fnam[] = "append_GRIDfile_text";

  //------- BEGIN ----------

  fprintf(fp_GRIDGEN_TEXT,
	  "# ============================================== \n"  );
  fprintf(fp_GRIDGEN_TEXT,"SNLC = %5d : " , GENLC.CID );

  // start by printing gen-parameters

  for ( IPAR=1; IPAR <= 4; IPAR++ ) {
    fprintf(fp_GRIDGEN_TEXT,"%s=%6.3f  "
	    ,SNGRID_WRITE.NAME[IPAR]
	    ,SNGRID_WRITE.CURRENT_VALUE[IPAR]
	    );
  }
  fprintf(fp_GRIDGEN_TEXT, "  (z=%6.4f)\n", GENLC.REDSHIFT_CMB );


  // write begin pad words
  fprintf(fp_GRIDGEN_TEXT,"%s",  "LCBEGIN PAD WORDS: ");
  for ( iwd=0; iwd < NPADWD_LCBEGIN; iwd++ ) {
    fprintf(fp_GRIDGEN_TEXT, "%5.5d  ", SNGRID_WRITE.I2GRIDGEN_LCMAG[iwd] );
  }
  fprintf(fp_GRIDGEN_TEXT,"%s",  "\n");


  NFILT  = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_FILTER] ; 
  NEP    = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_TREST] ; // NEP per filter

  // print Tobs value before mags
  fprintf(fp_GRIDGEN_TEXT,"   %5s: ", " Tobs" );
  for ( ep = 1; ep <= NEP; ep++ ) {
    NROW_WRITE_TOT++ ;
    jindx  = ifilt*NEP + (ep-1) ;
    Trest  = SNGRID_WRITE.VALUE[IPAR_GRIDGEN_TREST][ep];
    Tobs   = Trest * ( 1. + GENLC.REDSHIFT_HELIO );
    fprintf(fp_GRIDGEN_TEXT,"%5.1f ", Tobs );
  }
  fprintf(fp_GRIDGEN_TEXT,"\n" );


  for ( ifilt = 0; ifilt < NFILT; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];  
    sprintf(cfilt, "%c", FILTERSTRING[ifilt_obs] );

    // dump  mags
    fprintf(fp_GRIDGEN_TEXT,"   %s-mag: ", cfilt);
    for ( ep = 1; ep <= NEP; ep++ ) {
      jindx = ifilt*NEP + (ep-1) + NPADWD_LCBEGIN ;
      i2mag = SNGRID_WRITE.I2GRIDGEN_LCMAG[jindx]  ;
      fprintf(fp_GRIDGEN_TEXT,"%5d ", i2mag );
    }
    fprintf(fp_GRIDGEN_TEXT,"\n" );

    // dump mag errors

    fprintf(fp_GRIDGEN_TEXT,"   %s-err: ", cfilt);
    for ( ep = 1; ep <= NEP; ep++ ) {
      jindx  = ifilt*NEP + (ep-1)  + NPADWD_LCBEGIN ;
      i2err  = SNGRID_WRITE.I2GRIDGEN_LCERR[jindx]  ;
      fprintf(fp_GRIDGEN_TEXT,"%5d ", i2err );
    }
    fprintf(fp_GRIDGEN_TEXT,"\n" );

  } // ifilt


  // write end pad words
  fprintf(fp_GRIDGEN_TEXT,"%s",  "LCEND PAD WORDS: ");
  for ( iwd=0; iwd < NPADWD_LCBEGIN; iwd++ ) {
    jindx = NPADWD_LCBEGIN + SNGRID_WRITE.NGRIDGEN_PER_LC + iwd ;
    fprintf(fp_GRIDGEN_TEXT, "%5.5d  ", 
	    SNGRID_WRITE.I2GRIDGEN_LCMAG[jindx] );
  }
  fprintf(fp_GRIDGEN_TEXT,"%s",  "\n");


} // end of append_GRIDfile_text



// ============================================
void end_GRIDfile(void) {


  // Created Jun 9, 2012
  // Close files and 
  // make sure last light curve pointers is as expected.

  int istat;
  char fnam[] = "end_GRIDfile";

  // --------------- BEGIN --------

  sprintf(BANNER,"%s", fnam);
  print_banner(BANNER);


  if ( OPT_GRIDGEN_FORMAT == OPT_GRIDGEN_FORMAT_TEXT ) 
    { fclose(fp_GRIDGEN_TEXT); }
  else if ( OPT_GRIDGEN_FORMAT == OPT_GRIDGEN_FORMAT_FITS ) {
    istat = 0;
    fits_close_file(fp_GRIDGEN_FITS, &istat);
    check_fitserror("close fits file", istat) ;
  }

  // check light curve pointer from last update

  int ilc, NFILT, NEP, NPAD, IPTR_CALC;
 
  ilc = GENLC.CID ;

  printf("\t Total number of generated points: %d \n", NROW_WRITE_TOT);

  printf("\t Actual PTR_GRIDGEN_LC[last ilc=%8d] = %d \n", 
	 ilc, SNGRID_WRITE.PTR_GRIDGEN_LC[ilc]) ;

  NPAD   = NPADWD_LCBEGIN + NPADWD_LCEND ;
  NFILT  = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_FILTER] ; 
  NEP    = SNGRID_WRITE.NBIN[IPAR_GRIDGEN_TREST] ; // NEP per filter
  IPTR_CALC = 1 + ((NFILT*NEP)+NPAD)*(ilc-1);
  printf("\t Expected PTR_GRIDGEN_LC[last ilc]        = %d \n", IPTR_CALC);

  if ( IPTR_CALC != SNGRID_WRITE.PTR_GRIDGEN_LC[ilc] ) {
    sprintf(c1err,"Expected PTR_GRIDGEN_LC(last ilc) does not match.");
    sprintf(c1err,"See Actual vs. Expected above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
	  
  }

  fflush(stdout);

} // end of end_GRIDfile

// ===================================
int SNTYPE_GRIDGEN(void) {

  int ISMODEL_NON1ASED = (INDEX_GENMODEL == MODEL_NON1ASED);
  int ISMODEL_SIMSED   = (INDEX_GENMODEL == MODEL_SIMSED  );

  // Returns 1 for SNIa and 2 for non-Ia

  if ( ISMODEL_NON1ASED || ISMODEL_SIMSED ) {
    return SNTYPE_GRIDGEN_NONIa ;
  }
  else {
    return SNTYPE_GRIDGEN_Ia ;  // SALT2, snoopy, mlcs ...
  }


} // end of SNTYPE_GRIDGEN



// ************************************
void load_EXTNAME_GRIDGEN(void) {

  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_LOGZ],     "%s", EXTNAME_LOGZ);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_COLORPAR], "%s", EXTNAME_COLORPAR);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_COLORLAW], "%s", EXTNAME_COLORLAW);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_SHAPEPAR], "%s", EXTNAME_SHAPEPAR);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_FILTER],   "%s", EXTNAME_FILTER);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_TREST],    "%s", EXTNAME_TREST);

} // end  of load_EXTNAME_GRIDGEN(void)



// ========= END ===============
