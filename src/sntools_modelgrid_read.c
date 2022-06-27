/***********
 Created Oct 2010 by R.Kessler
 
 Tools to read MODEL GRID for simulation, psnid, etc ...
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


***************/

#include "sntools.h"      // general snana stuff
#include "fitsio.h"
#include "sntools_modelgrid.h"

#define ITYPE_PEC1A  11
#define ITYPE_MODEL  12

// ************************************
void load_EXTNAME_GRIDREAD(int IVERSION) {

  // Duplicate load_EXTNAME_GRIDGEN with different name to avoid
  // linking conflicts between  ifdef MODELGRID_GEN and MODELGRID_READ.
  // This is clearly quite clumsy ?!?!?!
  // Simulation needs GRIDGEN version while fitter needs the
  // READ version.
  //
  // Apr 4 2013: add IVERSION argument for back-compatibility

  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_LOGZ],     "%s", EXTNAME_LOGZ);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_COLORPAR], "%s", EXTNAME_COLORPAR);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_COLORLAW], "%s", EXTNAME_COLORLAW);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_SHAPEPAR], "%s", EXTNAME_SHAPEPAR);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_FILTER],   "%s", EXTNAME_FILTER);
  sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_TREST],    "%s", EXTNAME_TREST);

  if ( IVERSION < 2 ) {
    sprintf(EXTNAME_GRIDGEN[IPAR_GRIDGEN_SHAPEPAR], "%s", "LUMI-GRID" );
  }

} // end  of load_EXTNAME_GRIDREAD(void)


// *************************************************
void fits_read_SNGRID(int OPTMASK, char *sngridFile, 
		      SNGRID_DEF *SNGRID) {

  // Nov 5, 2010 R.Kessler
  // General utility to read fits file created by snlc_sim.exe
  // using the GRID option. This utility essentially fills
  // the same global arrays that were used to create the grid file.
  // 
  // OPTMASK bits (lsb=1)
  // bit1 : verbose mode
  //
  //
  // Aug 12 2016: check for PEC1A

  fitsfile *fp_SNGRID;

  int VERBOSE, IVERSION, ITYPE, istat, anynul, irow, colnum;
  int ipar, ifilt, nwdtmp, extver = 0 ;

  short  I2TMP[NPAR_GRIDGEN+1];
  int    I4TMP[NPAR_GRIDGEN+1];

  int   I4NULVAL  = -9 ;
  short I2NULVAL  = -9 ;
  float FNULVAL   = -9.0 ;
  char  CNULVAL[] = "NULL" ;

  long int NROW, FIRSTELEM, FIRSTROW ;

  char 
    fnam[] = "fits_read_SNGRID" 
    , keyname[60]
    , comment[100]
    , msg1[100]
    , *parname
    , *tparnam[NPAR_GRIDGEN+1]
    , *nonIa_type[MXGRIDGEN]
    , *nonIa_name[MXGRIDGEN]
    , extname[100]
    ;

  float     *ptrFloat ;

  // ---------------- BEGIN -------------------

  VERBOSE = (OPTMASK & 1 );

  printf("\n  Begin %s for \n  %s \n", fnam, sngridFile ) ;

  // open fits file for reading
  istat = 0;
  fits_open_file(&fp_SNGRID, sngridFile, READONLY, &istat );
  check_fitserror("Open FITS file", istat);


  // read internal version for back-compatibility
  sprintf(keyname, "%s", "IVERSION" );
  fits_read_key(fp_SNGRID, TINT, keyname, 
		&IVERSION, comment, &istat );
  if ( istat != 0 ) { IVERSION = 1 ; }
  SNGRID->IVERSION = IVERSION ;
  //  check_fitserror("read IVERSION key", istat);
  if ( VERBOSE ) {
    printf("    read fits keyword %s = %d \n", keyname, IVERSION );
    fflush(stdout);
  }


  // read optional UNIQUE_KEY (no call to check_fitserror)
  SNGRID->UNIQUE_KEY[0] = 0 ;
  sprintf(keyname, "%s", "GRIDKEY" );
  fits_read_key(fp_SNGRID, TSTRING, keyname, 
		&SNGRID->UNIQUE_KEY, comment, &istat );
  if ( istat == 0  )  { 
    printf("    read fits keyword %s = %s \n", keyname, SNGRID->UNIQUE_KEY ); 
    fflush(stdout);
  }

  istat = 0; // reset in case istat != 0 for UNIQUE_KEY

  // read name of survey
  sprintf(keyname, "%s", "SURVEY" );
  fits_read_key(fp_SNGRID, TSTRING, keyname, 
		&SNGRID->SURVEY, comment, &istat );
  check_fitserror("read SURVEY key", istat);
  if ( VERBOSE )  { 
    printf("    read fits keyword %s = %s \n", keyname, SNGRID->SURVEY ); 
    fflush(stdout);
  }


  // read name of SN model
  sprintf(keyname, "%s", "GENMODEL" );
  fits_read_key(fp_SNGRID, TSTRING, keyname, 
		&SNGRID->MODEL, comment, &istat );
  check_fitserror("read GENMODEL key", istat);
  if ( VERBOSE )    { 
    printf("    read fits keyword %s = %s \n", keyname, SNGRID->MODEL ); 
    fflush(stdout);
  }


  // read I2LCPACK to convert I*2 value back into relative flux
  sprintf(keyname, "%s", "I2LCPACK" );
  fits_read_key(fp_SNGRID, TFLOAT, keyname, 
		&GRIDGEN_I2LCPACK, comment, &istat );
  check_fitserror("read I2LCPACK key", istat);
  if ( VERBOSE ) {
    printf("    read fits keyword %s = %f \n", keyname, GRIDGEN_I2LCPACK ); 
    fflush(stdout);
  }

  
  // read list of filters
  sprintf(keyname, "%s", "FILTERS" );
  fits_read_key(fp_SNGRID, TSTRING, keyname, 
		&SNGRID->FILTERS, comment, &istat );
  check_fitserror("read FILTERS key", istat);
  if ( VERBOSE ) { 
    printf("    read fits keyword %s = %s \n", keyname, SNGRID->FILTERS ); 
    fflush(stdout);
  }
 

  // ===================================================
  //             Read SNPAR_INFO table
  // ===================================================


  sprintf(extname, "%s", EXTNAME_SNPAR_INFO );
  fits_movnam_hdu(fp_SNGRID, ANY_HDU, extname, extver, &istat );
  sprintf(msg1,"move to %s table", extname ) ;
  check_fitserror(msg1,  istat);

  // set read-params that are fixed for SNPAR_INFO table.
  NROW = NPAR_GRIDGEN;  FIRSTELEM = 1; FIRSTROW = 1;  

  for ( ipar=1; ipar <= NPAR_GRIDGEN; ipar++ ) 
    { tparnam[ipar] = SNGRID->NAME[ipar] ; }

  // read column of physical par-names
  colnum = 2;
  fits_read_col_str(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		    CNULVAL, &tparnam[1], &anynul, &istat );
  sprintf(msg1,"read PHYS-NAME column from %s table.", extname);
  check_fitserror(msg1, istat);

  // read NBIN column
  colnum = 3 ; 
  fits_read_col_sht(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		    I2NULVAL, &I2TMP[1], &anynul, &istat );

  sprintf(msg1,"read NBIN column from %s table.", extname );
  check_fitserror(msg1, istat);
  for ( ipar=1; ipar <= NPAR_GRIDGEN; ipar++ ) {
    SNGRID->NBIN[ipar]    = I2TMP[ipar] ;
  }

  // read ILCOFF column 
  // change 'sht' to 'int'
  colnum = 7 ; 
  fits_read_col_int(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		    I4NULVAL, &I4TMP[1], &anynul, &istat );
  sprintf(msg1,"read ILCOFF column from %s table.", extname );
  check_fitserror(msg1,istat);
  for ( ipar=1; ipar <= NPAR_GRIDGEN; ipar++ ) {
    SNGRID->ILCOFF[ipar]    = I4TMP[ipar] ;
  }


  // ================================================
  //    read bin values for each variable
  // ================================================


  load_EXTNAME_GRIDREAD(IVERSION) ;

  for ( ipar = 1; ipar <= NPAR_GRIDGEN ; ipar++ ) {

    parname = SNGRID->NAME[ipar] ;

    sprintf(extname, "%s", EXTNAME_GRIDGEN[ipar] );
    fits_movnam_hdu(fp_SNGRID, ANY_HDU, extname, extver, &istat );
    sprintf(msg1,"move to %s table (%s)", extname, parname );
    check_fitserror(msg1, istat);

    NROW    = SNGRID->NBIN[ipar] ;
    colnum  = 1;

    ptrFloat  = &SNGRID->VALUE[ipar][1] ;
    fits_read_col_flt(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		      FNULVAL, ptrFloat, &anynul, &istat );
    sprintf(msg1,"read %s bin-values from table", parname);
    check_fitserror(msg1, istat);

    SNGRID->VALMIN[ipar]  = SNGRID->VALUE[ipar][1] ;
    SNGRID->VALMAX[ipar]  = SNGRID->VALUE[ipar][NROW] ;

    // compute bin size
    float dif = SNGRID->VALUE[ipar][NROW] - SNGRID->VALUE[ipar][1] ;
    SNGRID->BINSIZE[ipar] = dif/(float)(NROW-1);
    //      SNGRID->VALUE[ipar][2] - SNGRID->VALUE[ipar][1] ;

    
    // read 2nd LAMAVG column for FILTER (Apr 2013)
    // Do NOT check for error in case of older GRID files.
    if ( ipar == IPAR_GRIDGEN_FILTER ) {

      if ( IVERSION >= 2 ) {
	colnum    = 2 ;
	ptrFloat  = &SNGRID->FILTER_LAMAVG[0] ;
	fits_read_col_flt(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
			  FNULVAL, ptrFloat, &anynul, &istat );
	sprintf(msg1,"read %s LAMAVG from table", parname);
	check_fitserror(msg1, istat);
      }
      else {
	for(ifilt=0 ; ifilt <= NROW; ifilt++ ) 
	  { SNGRID->FILTER_LAMAVG[ifilt] = -999. ; }
      }

      // Mar 2016: load IFILTOBS
      char band[2];
      for(ifilt=0 ; ifilt < NROW; ifilt++ )  {
	sprintf(band, "%c", SNGRID->FILTERS[ifilt] );
	SNGRID->IFILTOBS[ifilt] = INTFILTER(band); 
      }
    }
    
    if ( VERBOSE ) {
      printf("    %10.10s has %3d bins from %10.4f to %10.4f  (ILCOFF=%d)\n"
	     ,SNGRID->NAME[ipar]
	     ,SNGRID->NBIN[ipar]
	     ,SNGRID->VALMIN [ipar]
	     ,SNGRID->VALMAX[ipar]
	     ,SNGRID->ILCOFF[ipar]
	     );
    }

  } // end of ipar loop for table-binning


  // ===============================
  // read NONIa-info table 

  SNGRID->FRAC_PEC1A = 0.0 ; // init user-settable param (Aug 2016)
  ipar = IPAR_GRIDGEN_SHAPEPAR ;
  if ( strcmp(SNGRID->NAME[ipar],"NONIA" )==0 ) {
    sprintf(extname, "%s", EXTNAME_NONIa_INFO );
    fits_movnam_hdu(fp_SNGRID, ANY_HDU, extname, extver, &istat );
    sprintf(msg1,"move to %s table", extname );
    check_fitserror(msg1, istat);

    // set pointers
    NROW    = SNGRID->NBIN[ipar] ;
    for ( irow=1; irow <= NROW; irow++ )  { 
      nonIa_type[irow]    = SNGRID->NON1A_CTYPE[irow] ; 
      nonIa_name[irow]    = SNGRID->NON1A_NAME[irow] ; 
    }

    colnum = 2 ; 
    fits_read_col_int(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		      I4NULVAL, &SNGRID->NON1A_INDEX[1], &anynul, &istat );
    sprintf(msg1,"read %s table, NONIA-INDEX", extname );
    check_fitserror(msg1, istat);

    colnum = 3 ; 
    fits_read_col_str(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		      CNULVAL, &nonIa_type[1], &anynul, &istat );
    sprintf(msg1,"read NON1A_TYPE column from %s table.", extname);
    check_fitserror(msg1, istat);

    // July 26, 2012
    // classify each nonIa type as an  integer:
    // 1(Ib,Ic,Ibc ...)  or  2 (II, IIP, IIL ...) or 3(PEC1A)
    // This is to help group similar types into general I or II classes.
    // Feb 13 2017: ITYPE_AUTO(PEC1A) -> 3 (not 1)
    for ( irow=1; irow <= NROW; irow++ )  { 
      ITYPE = 
	get_NON1A_ITYPE_SNGRID( SNGRID->NON1A_CTYPE[irow] );

      SNGRID->ISPEC1A[irow]     = 0 ; // default
      
      if ( ITYPE == ITYPE_PEC1A )  { 
	SNGRID->NON1A_ITYPE_AUTO[irow] = 3 ; // Feb 2017 
	SNGRID->ISPEC1A[irow]     = 1 ;      // Aug 2016
      }
      else if ( ITYPE == ITYPE_MODEL ) {
	SNGRID->NON1A_ITYPE_AUTO[irow] = 4 ;
      }
      else  { 
	SNGRID->NON1A_ITYPE_AUTO[irow] = ITYPE ;  // Ia,II,Ibc
      } 
    }  // end irow

    colnum = 4 ; 
    fits_read_col_str(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		      CNULVAL, &nonIa_name[1], &anynul, &istat );
    sprintf(msg1,"read NON1A_NAME column from %s table.", extname);
    check_fitserror(msg1, istat);

    // read WGT and MAGOFF ( Aug 30 2013)
    if ( IVERSION >= 3 ) {
      colnum = 5 ; 
      fits_read_col_flt(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
			FNULVAL, &SNGRID->NON1A_WGT[1], &anynul, &istat );
      sprintf(msg1,"read NON1A_WGT column from %s table.", extname);
      check_fitserror(msg1, istat);

      colnum = 6 ; 
      fits_read_col_flt(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
			FNULVAL, &SNGRID->NON1A_MAGOFF[1], &anynul, &istat );
      sprintf(msg1,"read NON1A_MAGOFF column from %s table.", extname);
      check_fitserror(msg1, istat);

      colnum = 7 ; 
      fits_read_col_flt(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
			FNULVAL, &SNGRID->NON1A_MAGSMEAR[1], &anynul, &istat );
      sprintf(msg1,"read NON1A_MAGSMEAR column from %s table.", extname);
      check_fitserror(msg1, istat);
    }

    // Jan 2017: read  ITYPE_USER
    if ( SNGRID->IVERSION >= 4 ) {
      colnum = 8 ;
      fits_read_col_int(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
			I4NULVAL, &SNGRID->NON1A_ITYPE_USER[1], 
			&anynul, &istat );      
      sprintf(msg1,"read ITYPE_USER column from %s table.", extname );
      check_fitserror(msg1, istat);
    }

    /* 
    for ( irow=1; irow <= NROW; irow++ )  { 
      printf("\t xxxx NON1A row=%2d  INDEX=%3d  TYPE=%s  NAME=%s\n"
	     ,irow
	     ,SNGRID->NON1A_INDEX[irow]
	     ,SNGRID->NON1A_TYPE[irow]
	     ,SNGRID->NON1A_NAME[irow] );
    }
    */

  }

  // ===============================
  // read PTR_GRIDGEN_LC pointers

  sprintf(extname, "%s", EXTNAME_PTRI2LCMAG );
  fits_movnam_hdu(fp_SNGRID, ANY_HDU, extname, extver, &istat );
  sprintf(msg1,"move to %s table", extname );
  check_fitserror(msg1, istat);

  // read size of table
  sprintf(keyname, "%s", "NAXIS2" );
  fits_read_key(fp_SNGRID, TINT, keyname, &nwdtmp, comment, &istat );
  check_fitserror("read NAXIS2 key from PTR_I2LCMAG table", istat);

  // read pointer table
  SNGRID->NGRIDGEN_LC        = nwdtmp ;
  SNGRID->PTR_GRIDGEN_LC     = (int*)malloc( sizeof(int) * (nwdtmp+1) );
  NROW = nwdtmp ; colnum = 1; 
  fits_read_col_int(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		    I4NULVAL, &SNGRID->PTR_GRIDGEN_LC[1], &anynul, &istat );
  sprintf(msg1,"read %s table", extname );
  check_fitserror(msg1, istat);

  if ( VERBOSE ) {
    printf("    Read %d I2LCMAG pointers: %d %d %d .... %d \n", nwdtmp
	   ,SNGRID->PTR_GRIDGEN_LC[1]
	   ,SNGRID->PTR_GRIDGEN_LC[2]
	   ,SNGRID->PTR_GRIDGEN_LC[3]
	   ,SNGRID->PTR_GRIDGEN_LC[nwdtmp] );
  }

  // ============================================
  // finally read the I*2 mag & error arrays

  sprintf(extname, "%s", EXTNAME_I2LCMAG );
  fits_movnam_hdu(fp_SNGRID, ANY_HDU, extname, extver, &istat );

  sprintf(msg1,"move to %s table", extname );
  check_fitserror(msg1, istat);

  // read size of table
  sprintf(keyname, "%s", "NAXIS2" );
  fits_read_key(fp_SNGRID, TINT, keyname, &nwdtmp, comment, &istat );
  check_fitserror("read NAXIS2 key from I2LCMAG table", istat);
  NROW = nwdtmp ;

  SNGRID->SIZEOF_GRIDGEN = 4 * nwdtmp ; // July 9, 2012

  if ( VERBOSE ) {
    printf("\t read %d  I2LCMAG  values for GRID \n", nwdtmp );
    printf("\t read %d  I2LCERR  values for GRID \n", nwdtmp );

  }
  fflush(stdout) ;
  SNGRID->I2GRIDGEN_LCMAG = (short*)malloc( sizeof(short) * (nwdtmp+1) ) ;
  SNGRID->I2GRIDGEN_LCERR = (short*)malloc( sizeof(short) * (nwdtmp+1) ) ;
 
  colnum = 1 ; 
  fits_read_col_sht(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		    I2NULVAL, &SNGRID->I2GRIDGEN_LCMAG[1], &anynul, &istat );
  sprintf(msg1,"read %s table", extname );
  check_fitserror(msg1, istat);


  colnum = 2 ; 
  fits_read_col_sht(fp_SNGRID, colnum, FIRSTROW, FIRSTELEM, NROW, 
		    I2NULVAL, &SNGRID->I2GRIDGEN_LCERR[1], &anynul, &istat );
  check_fitserror("read I2LCERR table", istat);

  // close grid file
  fits_close_file(fp_SNGRID, &istat);
  check_fitserror("Close FITS file", istat) ;


  //  debugexit("fits read"); // xxxxxxxxx

} // end of read_GRIDfile_fits

// ==================================
void check_fitserror(char *comment, int status) {

  char fnam[] = "check_fitserror" ;

  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/

  if (status) {
    fits_report_error(stderr, status); /* print error report */
    errmsg(SEV_FATAL, 0, fnam, comment, "Check cfitsio routines." ); 
  }
  return;
}



// ********************************************
int INDEX_GRIDGEN(int ipar, double parval, SNGRID_DEF *SNGRID) {

  // return index of ipar-parameter with value parval.
  double valmin, valmax, valbin, dif, ratio ;
  int indx ;
  // ---------- BEGIN
  indx = -9;
  valmin = (double)SNGRID->VALMIN[ipar] ;
  valmax = (double)SNGRID->VALMAX[ipar] ;
  valbin = (double)SNGRID->BINSIZE[ipar] ;

  if ( parval <= valmin ) 
    { indx = 1 ; }
  else if ( parval >= valmax ) 
    { indx = SNGRID->NBIN[ipar] ; }
  else  { 
    dif   = parval - valmin ;
    ratio = dif / valbin ;
    indx  = 1 + (int)ratio  ;  
  }

  return(indx) ;

} // end of INDEX_GRIDGEN


// ===============================================
int  get_NON1A_ITYPE_SNGRID(char *NONIA_CTYPE ) {

  // Aug 12 2016: check for PEC1A
  // Aug 27 2017: check MODEL1[2,3,4]
  
  int L0, L1, ISPEC1A, ISMODEL ;
  char c0[2], c1[2];
  char fnam[] = "get_NON1A_ITYPE_SNGRID" ;

  // ------------------ BEGIN ----------------

  sprintf( c0,"%c", NONIA_CTYPE[0] ); 
  sprintf( c1,"%c", NONIA_CTYPE[1] ); 
  ISPEC1A = ( strstr(NONIA_CTYPE,"PEC1A") != NULL ) ;
  ISMODEL = ( strstr(NONIA_CTYPE,"MODEL") != NULL ) ;

  L0 = L1 = 0;
  if ( strcmp(c0,"I") == 0 ) { L0 = 1; }
  if ( strcmp(c1,"I") == 0 ) { L1 = 1; }

  if ( L0 && L1 ) 
    { return 2; }  // II[P,n,L ..]
  else if ( L0 ) 
    { return 1 ; } // I[b,c ...]
  else if ( ISPEC1A ) 
    { return ITYPE_PEC1A ; }  // pec Ia (Aug 2016)
  else if ( ISMODEL ) 
    { return ITYPE_MODEL ; }  // MODEL1[2,3,3] Aug 28 2017
  else {
    sprintf(c1err,"Invalid NONIA_TYPE = '%s'", NONIA_CTYPE);
    sprintf(c2err,"Must begin with 'I', 'II', or 'PEC1A' ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return(-9); // should never get here

} // end of get_NON1A_ITYPE_SNGRID


// - - - - - - - - - - - - - - - - -
//   Utilities after reading
// - - - - - - - - - - - - - - - - -



void renorm_wgts_SNGRID(SNGRID_DEF *SNGRID ) {

  // Created Mar 2016 by R.Kessler 
  // Re-normalize wgts so that sum = 1.0      
  //
  // Aug 15 2016: refactor to handle FRAC_PEC1A

  int i, ISPEC1A ;
  double SUM, SUMALL[2], FRAC[2]  ;
  double FRAC_PEC1A  = (double)SNGRID->FRAC_PEC1A ;
  double FRAC_NON1A  = 1.0 - FRAC_PEC1A ;
  int    NTMPL       = SNGRID->NBIN[IPAR_GRIDGEN_SHAPEPAR];
  //  char fnam[] = "renorm_wgts_SNGRID" ;

  // ----------- BEGIN -------------

  SUMALL[0] = SUMALL[1] = 0.0 ; // 0->NON1A, 1->PEC1A
  FRAC[0] = FRAC_NON1A ;
  FRAC[1] = FRAC_PEC1A ;

  for ( i=0; i <= NTMPL ; i++ ) {

    SNGRID->NON1A_WGTSUM[i] = 0.0 ; // init

    if ( i == 0 ) { continue ; }

    // ignore excluded templates.          
    if ( SNGRID->NON1A_ITYPE_AUTO[i] < 0 ) { continue; }

    // check for PEC1A (Aug 2016)
    ISPEC1A = SNGRID->ISPEC1A[i] ;

    // increment sum separately for NON1A & PEC1A
    SUMALL[ISPEC1A] += SNGRID->NON1A_WGT[i];
  }

  if ( SUMALL[0] == 0.0 ) { return ; }

  // apply renorm.
  SUM = 0.0 ;
  for ( i=1; i <= NTMPL ; i++ )  {
    ISPEC1A = SNGRID->ISPEC1A[i] ; 
    SNGRID->NON1A_WGT[i]   /= SUMALL[ISPEC1A] ; 
    SNGRID->NON1A_WGT[i]   *= FRAC[ISPEC1A];
    SUM                    += SNGRID->NON1A_WGT[i] ;
    SNGRID->NON1A_WGTSUM[i] = SUM ;   // incremental sum
  }

  return ;

} // end of renorm_wgts_SNGRID


// ==================================================
void sort_by_PEC1A_SNGRID(SNGRID_DEF *SNGRID ) {

  // Created Aug 15 2016
  // Sort the NON1A elements so that all of the PEC1A
  // are at the end of the list. This allows the RANDOM 
  // weight to determine NON1A or PEC1A.

  SNGRID_DEF *SNGRID_TMP ;
  int    NTMPL       = SNGRID->NBIN[IPAR_GRIDGEN_SHAPEPAR];
  int    i, itmp, inew, IOFF_PEC1A, ISPEC1A, NNON1A, NPEC1A=0 ;
  //  char fnam[] = "sort_by_PEC1A_SNGRID" ;

  // ----------- BEGIN --------------

  // count how many are PEC1A
  for(i=1; i<=NTMPL; i++ ) { if ( SNGRID->ISPEC1A[i] ) { NPEC1A++ ; }  }

  if ( NPEC1A == 0  ) { return ; }

  // create local copy of SNGRID->NON1A_xxx
  SNGRID_TMP = (SNGRID_DEF*) malloc ( sizeof(SNGRID_DEF) ) ;
  for(i=1; i<=NTMPL; i++ ) { copy_NON1A_SNGRID(i,i, SNGRID, SNGRID_TMP ); }

  // now transfer NON1A and PEC1A from SNGRID_TMP back to SNGRID
  IOFF_PEC1A = NTMPL - NPEC1A  ; // fortran-like indices start at 1
  NNON1A = NPEC1A = 0 ;
  for(itmp=1; itmp<=NTMPL; itmp++ ) { 

    ISPEC1A = SNGRID_TMP->ISPEC1A[itmp] ;
    
    if ( ISPEC1A == 0 ) 
      { NNON1A++; inew = NNON1A ;    }
    else  
      { NPEC1A++ ; inew = IOFF_PEC1A+NPEC1A ;  }

    copy_NON1A_SNGRID(itmp,inew, SNGRID_TMP, SNGRID );  
  }

  free(SNGRID_TMP);
  return ;

} //end sort_by_PEC1A_SNGRID


// ==================================================
void copy_NON1A_SNGRID(int IROW1, int IROW2, 
		       SNGRID_DEF *SNGRID1, SNGRID_DEF *SNGRID2 ) {

  // Created Aug 15 2016
  // For NON1A row IROW, copy all 
  //  SNGRID2->NON1A_xxxx[IROW2] = SNGRID1->NON1A_xxx[IROW1]

  //  char fnam[] = "copy_NON1A_SNGRID" ;

  // ----------- BEGIN --------------

  SNGRID2->NON1A_INDEX[IROW2]       = SNGRID1->NON1A_INDEX[IROW1];
  SNGRID2->NON1A_ITYPE_AUTO[IROW2]  = SNGRID1->NON1A_ITYPE_AUTO[IROW1];
  SNGRID2->NON1A_ITYPE_USER[IROW2]  = SNGRID1->NON1A_ITYPE_USER[IROW1];
  SNGRID2->NON1A_WGT[IROW2]      = SNGRID1->NON1A_WGT[IROW1];
  SNGRID2->NON1A_WGTSUM[IROW2]   = SNGRID1->NON1A_WGTSUM[IROW1];
  SNGRID2->NON1A_MAGOFF[IROW2]   = SNGRID1->NON1A_MAGOFF[IROW1];
  SNGRID2->NON1A_MAGSMEAR[IROW2] = SNGRID1->NON1A_MAGSMEAR[IROW1];
  SNGRID2->ISPEC1A[IROW2]        = SNGRID1->ISPEC1A[IROW1];

  sprintf(SNGRID2->NON1A_NAME[IROW2],  "%s", SNGRID1->NON1A_NAME[IROW1] );
  sprintf(SNGRID2->NON1A_CTYPE[IROW2], "%s", SNGRID1->NON1A_CTYPE[IROW1] );
  return ;

} // end copy_NON1A_SNGRID


// ==================================================
void dump_SNGRID(SNGRID_DEF *SNGRID ) {

  int i, NTMPL, IPAR, LEN ;
  char *NAME, dumpName[40] ;
  //  char fnam[] = "dump_SNGRID" ;

  // --------------- BEGIN ---------------

  printf("\n NONIA Grid Summary: \n") ;


  printf("          Sim-  \n");
  printf("   sparse input                                   ITYPE      \n");
  printf("   index  index     Name             Type       USER(AUTO)  WGT Msmear \n");
  printf("   ------------------------------------------------------------------ \n");
  IPAR     = IPAR_GRIDGEN_SHAPEPAR ;
  NTMPL    = SNGRID->NBIN[IPAR];

  for ( i=1; i <= NTMPL ; i++ ) {

    // if name is more than 20 chars then print last 20 chars
    NAME = SNGRID->NON1A_NAME[i] ;
    LEN  = strlen(NAME);
    // printf(" xxx LEN=%d for NAME='%s' \n", LEN, NAME); fflush(stdout);

    if ( LEN >= 20 ) 
      { sprintf(dumpName, "%s", &NAME[LEN-20] ); }
    else
      { sprintf(dumpName, "%s", NAME ); }
    
      printf("   %3d  %4d    %-20.20s  %-12.12s  %2d(%d)  %.3f %.2f\n"
	     ,i
	     ,SNGRID->NON1A_INDEX[i]
	     ,dumpName                     // SNGRID->NON1A_NAME[i]
	     ,SNGRID->NON1A_CTYPE[i]
	     ,SNGRID->NON1A_ITYPE_USER[i]
	     ,SNGRID->NON1A_ITYPE_AUTO[i]
	     ,SNGRID->NON1A_WGT[i]
	     ,SNGRID->NON1A_MAGSMEAR[i]
	     );
      
  }

  printf("\n");

  return ;
}  // end dump_SNGRID



// ========= END ===============
