/***********************
  Created May 2020 by R. Kessler

  Write host-galaxy spectra to .fits file that is specifically
  formatted so that it can be input to the online MARZ program 
  (http://samreay.github.io/Marz) to determine redshifts.

  This utility is invoked with interactive snana.exe program and
    &SNLCINP
      MARZFILE_OUT = '***marz***.fits'

  The output file name must have a .fits extension, and 
  either 'marz' or 'MARZ' must be part of the fileName.
  MARZFILE_OUT will not work with fitting codes (snlc_fit, psnid),
  and will not work in batch mode.

  Since images and tables must be written sequentially to fits files,
  each host spectrum is stored in memory (MARZ_SPECDATA struct) 
  and everything is written in the CLOSE_MARZFILE process.
  MXSPEC_MARZ limits the number of stored spectra to avoid 
  using excessive memory.
 
 ***********************/

#include "fitsio.h"

// ======= GLOBALS =========

#define MXSPEC_MARZ 1000   // stop writing after this many to conserve memory
#define MARZTABLE_FIBRES    "FIBRES"
#define MARZTABLE_FLAM      "FLAM"
#define MARZTABLE_VARIANCE  "VARIANCE"
#define MARZTABLE_WAVE      "WAVELENGTH"
#define MXCHAR_CCID_MARZ    32   // allows [CID](z=x.xxx)
#define FLAMSCALE_MARZ 1.0E19 // scale FLAM to avoid float overflow

fitsfile *FP_MARZ ;

struct MARZ_SPECDATA {  
  int    NBIN_WAVE; // number of wave bins 
  float  *WAVE ;    // wavelength in each bin
  int    NSPEC;     // number of stored spectra
  char  **NAME;     // name of each SN vs. ispec
  float **FLAM ;    // FLAM vs. wave bin and ispec
  float **VARFLAM ; // variance vs. wave bin and ispec
} MARZ_SPECDATA ;

// ======== FUNCTIONS ==========

#ifdef __cplusplus
extern"C" {
#endif
  void OPEN_MARZFILE(char *FILENAME, int *ERR);
  void WRITE_TABLES_MARZ(char *FILENAME);
  void CLOSE_MARZFILE(char *FILENAME);
  void SNTABLE_CREATE_MARZ(int IDTABLE, char *NAME);
  void SPECPAK_FILL_MARZ(void);
  void marzfitsio_errorCheck(char *comment, int status);

#ifdef __cplusplus
}
#endif



// ========================================  
void marzfitsio_errorCheck(char *comment, int status) {
  char fnam[] = "marzfitsio_errorCheck" ;
  char comment2[] = "Check cfitsio routines.";
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/
  if (status) {
    fits_report_error(stderr, status); /* print error report */
    errmsg(SEV_FATAL, 0, fnam, comment, comment2);
  }
  return;
} 

// ======================================
void OPEN_MARZFILE(char *FILENAME, int *ERR) {

  int istat = 0 ;
  char msg[200], clobberFile[200];
  char fnam[] = "OPEN_MARZFILE";
  // ------------ BEGIN ------------

  SPECPAK_USE_MARZ = true ;
  *ERR = 0;

  sprintf(clobberFile, "!%s", FILENAME);
  sprintf(msg,"%s %s", fnam, FILENAME);
  print_banner(msg);
  fits_create_file(&FP_MARZ, clobberFile, &istat );
  sprintf(msg,"Open %s", FILENAME );
  marzfitsio_errorCheck(msg, istat);

  return;

} // end OPEN_MARZFILE


// ======================================
void CLOSE_MARZFILE(char *FILENAME) {

  int istat = 0;
  char msg[200];
  char fnam[] = "CLOSE_MARZFILE" ;
  // -------------- BEGIN -------------

  WRITE_TABLES_MARZ(FILENAME);

  printf("   Close marz file %s \n", FILENAME);
  fits_close_file(FP_MARZ, &istat);
  sprintf(msg, "Close FITS file %s", FILENAME );
  marzfitsio_errorCheck(msg, istat);

  fflush(stdout);

  return ;

} // end CLOSE_MARZFILE


// ===================================
void WRITE_TABLES_MARZ(char *FILENAME) {

  int NSPEC     = MARZ_SPECDATA.NSPEC ;
  int NBIN_WAVE = MARZ_SPECDATA.NBIN_WAVE ;
  int istat     = 0;
  long NROW     = (long)NSPEC ;
  long fpixel=1, N2D = (long)(NSPEC*NBIN_WAVE) ;
  long NAXIS=2, NAXES[2] = { (long)NBIN_WAVE, (long)NSPEC } ;
  int  ncol, icol, iwave, ispec, colnum, firstelem, firstrow, i2 ;

  int  MXCOL   = 10 ;
  int  MEMC1   = MXCOL*sizeof(char*);
  int  MEMC0   = MXCHAR_CCID_MARZ*sizeof(char);
  
  char extName[40], **varNames, **formats, **units;
  char tblName[40], msg[100] ;
  char fnam[] = "WRITE_TABLES_MARZ" ;

  // ------------ BEGIN ------------

  printf("   Write marz file %s \n", FILENAME);
  printf("      (NSPEC=%d, %d WAVE bins from %.1f to %.1f A) \n", 
	 NSPEC, NBIN_WAVE, MARZ_SPECDATA.WAVE[0], MARZ_SPECDATA.WAVE[NBIN_WAVE-1] ) ;

  varNames = (char**) malloc(MEMC1);
  formats  = (char**) malloc(MEMC1);
  units    = (char**) malloc(MEMC1);
  for(icol=0; icol < MXCOL; icol++ ) {
    varNames[icol] = (char*) malloc(MEMC0);
    formats[icol]  = (char*) malloc(MEMC0);
    units[icol]    = (char*) malloc(MEMC0);
    sprintf(units[icol]," ");
  }

  // - - - - - - -
  // Create NBIN_WAVE x NSPEC 1D array for image storage
  int MEMF = NBIN_WAVE * NSPEC * sizeof(float);
  float *ARRAY1D;  ARRAY1D = (float*) malloc(MEMF);

  // - - - - - - - - - - - 
  // make primary image for FLAM for each spectrum
  fits_create_img(FP_MARZ, FLOAT_IMG, NAXIS, NAXES, &istat) ;
  sprintf(msg,"Create FLAM primary image") ;
  marzfitsio_errorCheck(msg, istat) ;

  i2 = 0;
  for(ispec=0; ispec < NSPEC; ispec++ ) {
    for(iwave=0; iwave < NBIN_WAVE; iwave++ ) {
      ARRAY1D[i2] = MARZ_SPECDATA.FLAM[ispec][iwave];      i2++;
    }
  }
  fits_write_img(FP_MARZ, TFLOAT, fpixel, N2D, ARRAY1D, &istat );
  sprintf(msg,"Write FLAM primary image") ;
  marzfitsio_errorCheck(msg, istat) ;

  // - - - - - - - - - - - - 
  // load some header keys

  // name of survey                                                               
  fits_update_key(FP_MARZ, TSTRING, "SURVEY",
                  SPECPAK_OUTPUT.SURVEY, "Survey", &istat );

  // name of survey                                                               
  fits_update_key(FP_MARZ, TSTRING, "VERSION_PHOTOMETRY",
                  SPECPAK_OUTPUT.VERSION_PHOTOMETRY, "Phot version", &istat );

  // flux scale
  float scale = (float)FLAMSCALE_MARZ;
  fits_update_key(FP_MARZ, TFLOAT, "FLAM_SCALE",
                  &scale, "FLAM here scaled by this factor", &istat );

  // - - - - - - - - - - - 
  // make image for VARFLAM for each spectrum
  fits_create_img(FP_MARZ, FLOAT_IMG, NAXIS, NAXES, &istat) ;
  sprintf(msg,"Create VARIANCE image") ;
  marzfitsio_errorCheck(msg, istat) ;

  i2 = 0;
  for(ispec=0; ispec < NSPEC; ispec++ ) {
    for(iwave=0; iwave < NBIN_WAVE; iwave++ ) {
      ARRAY1D[i2] = MARZ_SPECDATA.VARFLAM[ispec][iwave];      i2++;
    }
  }
  fits_write_img(FP_MARZ, TFLOAT, fpixel, N2D, ARRAY1D, &istat );
  sprintf(msg,"Write VARIANCE image") ;
  marzfitsio_errorCheck(msg, istat) ;

  sprintf(extName,"VARIANCE");
  fits_update_key(FP_MARZ, TSTRING, "EXTNAME", extName, "image name", &istat );


  // - - - - - - - - - - - 
  // make image for WAVELENGTH for each spectrum
  long NAXIS_WAVE = 1;
  long NAXES_WAVE = (long)NBIN_WAVE;
  fits_create_img(FP_MARZ, FLOAT_IMG, NAXIS_WAVE, &NAXES_WAVE, &istat) ;
  sprintf(msg,"Create WAVELENGTH image") ;
  marzfitsio_errorCheck(msg, istat) ;

  fits_write_img(FP_MARZ, TFLOAT, fpixel, NAXES_WAVE, MARZ_SPECDATA.WAVE, &istat );
  sprintf(msg,"Write WAVELENGTH image") ;
  marzfitsio_errorCheck(msg, istat) ;

  sprintf(extName,"WAVELENGTH");
  fits_update_key(FP_MARZ, TSTRING, "EXTNAME", extName, "image name", &istat );

  // -----------------------
  // write table with SN NAMEs
  ncol = 1;  
  sprintf(tblName,"%s", MARZTABLE_FIBRES);
  sprintf(varNames[0],"NAME");
  sprintf(formats[0], "%dA", MXCHAR_CCID_MARZ);
  fits_create_tbl(FP_MARZ, BINARY_TBL, NROW, ncol,
		  varNames, formats, units, tblName,  &istat );
  sprintf(msg, "Create %s table", tblName );
  marzfitsio_errorCheck(msg, istat);

  firstrow = firstelem = 1; colnum=1;  
  fits_write_col(FP_MARZ, TSTRING, colnum, firstrow, firstelem, NSPEC,
		 MARZ_SPECDATA.NAME, &istat);
  sprintf(msg, "Write %s table", tblName );
  marzfitsio_errorCheck(msg, istat);



  return ;

} // end WRITE_TABLES_MARZ


// =====================================
void SNTABLE_CREATE_MARZ(int IDTABLE, char *NAME) {

  char fnam[] = "SNTABLE_CREATE_MARZ" ;

  // -------------- BEGIN --------------

  printf("  %s: prep host-spectra storage for MARZ \n\n", fnam );
  fflush(stdout);

  MARZ_SPECDATA.NSPEC     = 0;
  MARZ_SPECDATA.NBIN_WAVE = 0;

  // malloc EVT dimension here; for each event, malloc 
  // WAVE dimension in SPECPAK_FILL_MARZ
  int MEMF = MXSPEC_MARZ * sizeof(float*);
  int MEMC = MXSPEC_MARZ * sizeof(char*);
  MARZ_SPECDATA.FLAM    = (float**)malloc(MEMF);
  MARZ_SPECDATA.VARFLAM = (float**)malloc(MEMF);
  MARZ_SPECDATA.NAME    = (char**)malloc(MEMC);

  return ;

} // end SNTABLE_CREATE_MARZ

// ====================================
void SPECPAK_FILL_MARZ(void) {

  // Store all spectra in MARZ_SPECDATA structure.
  // Everything is written at the end in CLOSE_MARZFILE
  // becauyse FITS tables must be written sequentially.

  int NSPEC     = MARZ_SPECDATA.NSPEC;
  int NBIN_WAVE = SPECPAK_OUTPUT.NLAMBIN_LIST[0] ;
  int MEMF      = NBIN_WAVE * sizeof(float);
  int MEMC      = MXCHAR_CCID_MARZ*sizeof(char);
  double *PTR_FLAM    = SPECPAK_OUTPUT.FLAM ;
  double *PTR_FLAMERR = SPECPAK_OUTPUT.FLAMERR ;
  char *CCID          = SPECPAK_OUTPUT.CCID;

  double  FLAM, FLAMERR, VARFLAM, LAMMIN, LAMMAX, LAMAVG ;
  int iwave;
  char fnam[] = "SPECPAK_FILL_MARZ";

  // ------------- BEGIN --------------

  // avoid using too much memory
  if ( NSPEC >= MXSPEC_MARZ ) { return; }

  // allocate memory for this spectrum
  MARZ_SPECDATA.FLAM[NSPEC]    = (float*)malloc(MEMF);
  MARZ_SPECDATA.VARFLAM[NSPEC] = (float*)malloc(MEMF);
  MARZ_SPECDATA.NAME[NSPEC]    = (char *)malloc(MEMC);
  if ( NSPEC == 0 ) { 
    MARZ_SPECDATA.NBIN_WAVE = NBIN_WAVE ;
    MARZ_SPECDATA.WAVE      = (float*)malloc(MEMF); 
  }

  /*  
  printf(" xxx %s: NBIN_WAVE=%d  CCID=%s  ID=%d\n", 
	 fnam, NBIN_WAVE, CCID, SPECPAK_OUTPUT.ID_LIST[0]);
  */

  fflush(stdout);

  sprintf(MARZ_SPECDATA.NAME[NSPEC],"%s", CCID);

  for(iwave = 0; iwave < NBIN_WAVE; iwave++ ) {
    FLAM     = FLAMSCALE_MARZ * PTR_FLAM[iwave];
    FLAMERR  = FLAMSCALE_MARZ * PTR_FLAMERR[iwave];
    VARFLAM  = FLAMERR*FLAMERR;

    MARZ_SPECDATA.FLAM[NSPEC][iwave]    = (float)FLAM;
    MARZ_SPECDATA.VARFLAM[NSPEC][iwave] = (float)VARFLAM;

    // store wavelength grid on first spectrum
    if ( NSPEC == 0 ) {
      LAMMIN = SPECPAK_OUTPUT.LAMMIN[iwave];
      LAMMAX = SPECPAK_OUTPUT.LAMMAX[iwave];
      LAMAVG = 0.5*(LAMMIN + LAMMAX);
      MARZ_SPECDATA.WAVE[iwave] = LAMAVG ;
    }
  }

  // abort if wavelength binning changes
  if ( MARZ_SPECDATA.NBIN_WAVE != SPECPAK_OUTPUT.NLAMBIN_LIST[0] ) {
    sprintf(MSGERR1,"Stored NBIN_WAVE=%d for first event ...", 
	    MARZ_SPECDATA.NBIN_WAVE);
    sprintf(MSGERR2,"but NBIN(CID=%s) = %d. Cannot change binning.",
	    CCID, SPECPAK_OUTPUT.NLAMBIN_LIST[0] );      
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  MARZ_SPECDATA.NSPEC++ ;

  return;

} // end SPECPAK_FILL_MARZ

