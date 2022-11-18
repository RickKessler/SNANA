/**************************************************
  Created Oct 2019
  Translate fortran kcor-read utilities into C

  Test/debug with 
    snlc_sim.exe <inFile> USE_KCOR_REFACTOR 1

  November 2022:
  The original "kcor" notation referred to the kcor.exe program that
  produced a table of K-corrections for MLCS. Later, filter and calib
  info was added to the k-cor table storage, and with the invention
  of SALT2, kcor tables no longer contained K-corrections.

  Here the generic "kcor" notation has been replaced with "calib",
  except when referring to explicit K-cor information.


  NOT READY !!

***************************************************/

#include "fitsio.h"
#include "sntools.h"
#include "sntools_data.h"
#include "sntools_dataformat_fits.h"
#include "sntools_calib.h" 
#include "sntools_spectrograph.h"
#include "MWgaldust.h"


// ======================================
void READ_CALIB_DRIVER(char *kcorFile, char *FILTERS_SURVEY,
		      double *MAGREST_SHIFT_PRIMARY,
                      double *MAGOBS_SHIFT_PRIMARY ) {

  int ifilt;
  char BANNER[100];
  char fnam[] = "READ_CALIB_DRIVER" ;

  // ----------- BEGIN ---------------

  sprintf(BANNER,"%s: read calib/filters/kcor", fnam );
  print_banner(BANNER);

  // - - - - - - - 
  // store passed info in global struct
  sprintf(CALIB_INFO.FILENAME, "%s", kcorFile) ;
  sprintf(CALIB_INFO.FILTERS_SURVEY, "%s", FILTERS_SURVEY);
  CALIB_INFO.NFILTDEF_SURVEY = strlen(FILTERS_SURVEY);
  for(ifilt=0; ifilt < MXFILT_CALIB; ifilt++ ) {
    CALIB_INFO.MAGREST_SHIFT_PRIMARY[ifilt] = MAGREST_SHIFT_PRIMARY[ifilt];
    CALIB_INFO.MAGOBS_SHIFT_PRIMARY[ifilt]  = MAGOBS_SHIFT_PRIMARY[ifilt];
  }

  // - - - - - - 
  read_calib_init();

  read_calib_open();

  read_calib_head();

  read_calib_zpoff();

  read_calib_snsed();

  read_kcor_tables();

  read_calib_mags();

  read_calib_filters();

  // pass dump flags
  int DO_DUMP = KCOR_VERBOSE_FLAG ;
  if ( DO_DUMP ) {
    char BLANK[] = "";
    addFilter_kcor(777, BLANK, &CALIB_INFO.FILTERCAL_REST);
    addFilter_kcor(777, BLANK, &CALIB_INFO.FILTERCAL_OBS );
    printf("\n\n");
  }

  read_calib_primarysed();

  print_calib_summary();


  int istat = 0 ;
  fits_close_file(CALIB_INFO.FP, &istat); 

  printf("\n %s: Done \n", fnam);
  fflush(stdout);

  return ;

} // end READ_CALIB_DRIVER

void read_calib_driver__(char *kcorFile, char *FILTERS_SURVEY,
                        double *MAGREST_SHIFT_PRIMARY,
                        double *MAGOBS_SHIFT_PRIMARY ) {
  READ_CALIB_DRIVER(kcorFile, FILTERS_SURVEY, 
		   MAGREST_SHIFT_PRIMARY, MAGOBS_SHIFT_PRIMARY );
} // end read_calib_driver__

// ===============================
void read_calib_init(void) {

  int i, i2;
  char *kcorFile = CALIB_INFO.FILENAME;
  char BLANK[] = "";
  char fnam[] = "read_calib_init" ;

  // ------------ BEGIN ------------

  if ( IGNOREFILE(kcorFile) ) {
    sprintf(c1err,"Must specifiy kcor/calib file with KCOR_FILE key");
    sprintf(c2err,"KCOR_FILE contains filter trans, primary ref, AB off, etc");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  CALIB_INFO.NCALL_READ++ ;
  CALIB_INFO.NKCOR             = 0 ;
  CALIB_INFO.NKCOR_STORE       = 0 ;
  CALIB_INFO.NFILTDEF          = 0 ;
  CALIB_INFO.NPRIMARY          = 0 ;
  CALIB_INFO.NFLAMSHIFT        = 0 ;
  CALIB_INFO.RVMW              = RV_MWDUST ;
  CALIB_INFO.OPT_MWCOLORLAW    = OPT_MWCOLORLAW_ODON94 ;
  CALIB_INFO.STANDALONE        = false ;
  CALIB_INFO.MASK_EXIST_BXFILT = 0 ;

  for(i=0; i < MXFILT_CALIB; i++ ) {
    CALIB_INFO.IFILTDEF[i]   = -9 ;
    CALIB_INFO.ISLAMSHIFT[i] = false;
    CALIB_INFO.MASK_FRAME_FILTER[i] = 0 ;
    CALIB_INFO.IS_SURVEY_FILTER[i]  = false ;

    for(i2=0; i2 < MXFILT_CALIB; i2++ ) 
      { CALIB_INFO.EXIST_KCOR[i][i2] = false ; }
  }

  for(i=0; i < MXTABLE_KCOR; i++ ) {
    CALIB_INFO.IFILTMAP_KCOR[OPT_FRAME_REST][i] = -9; 
    CALIB_INFO.IFILTMAP_KCOR[OPT_FRAME_OBS][i]  = -9;
  }

  KCOR_VERBOSE_FLAG = 0 ;

  addFilter_kcor(0, BLANK, &CALIB_INFO.FILTERCAL_REST); // zero map
  addFilter_kcor(0, BLANK, &CALIB_INFO.FILTERCAL_OBS ); // zero map
  CALIB_INFO.FILTERCAL_REST.OPT_FRAME = OPT_FRAME_REST ;
  CALIB_INFO.FILTERCAL_OBS.OPT_FRAME  = OPT_FRAME_OBS  ;

  IFILTDEF_BESS_BX = INTFILTER("X");
 
  return ;

} // end read_calib_init


// ===============================
void  read_calib_open(void) {

  int istat=0;
  char kcorFile[MXPATHLEN] ;
  char fnam[] = "read_calib_open" ;

  // -------- BEGIN --------

  sprintf(kcorFile, "%s", CALIB_INFO.FILENAME);
  fits_open_file(&CALIB_INFO.FP, kcorFile, READONLY, &istat);

  // if kcorFile doesn't exist, try reading from SNDATA_ROOT
  if ( istat != 0 ) {
    istat = 0;
    sprintf(kcorFile,"%s/kcor/%s", 
	    getenv("SNDATA_ROOT"), CALIB_INFO.FILENAME);
    fits_open_file(&CALIB_INFO.FP, kcorFile, READONLY, &istat);
  }

  if ( istat != 0 ) {					       
    sprintf(c1err,"Could not open %s", CALIB_INFO.FILENAME );  
    sprintf(c2err,"Check local dir and $SNDATA_ROOT/kcor");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("  Opened %s\n", kcorFile); fflush(stdout);

  return ;

} // end read_calib_open


// =============================
void read_calib_head(void) {

  // Read FITS header keys

  fitsfile *FP = CALIB_INFO.FP;
  int istat = 0, i, IFILTDEF, len ;
  char KEYWORD[40], comment[100], tmpStr[100], cfilt[2] ;
  int  NUMPRINT =  14;
  int  MEMC     = 80*sizeof(char); // for malloc
  //  int  memc     =  8*sizeof(char); // for malloc

  int *IPTR; double *DPTR; char *SPTR;
  char fnam[] = "read_calib_head" ;

  // --------- BEGIN ---------- 
  printf("   %s \n", fnam); fflush(stdout);

  sprintf(KEYWORD,"VERSION");  IPTR = &CALIB_INFO.VERSION;
  fits_read_key(FP, TINT, KEYWORD, IPTR, comment, &istat);
  snfitsio_errorCheck("can't read VERSION", istat);
  printf("\t\t Read %-*s  = %d  (kcor.exe version) \n", 
	 NUMPRINT, KEYWORD, *IPTR);

  sprintf(KEYWORD,"NPRIM");  IPTR = &CALIB_INFO.NPRIMARY;
  fits_read_key(FP, TINT, KEYWORD, IPTR, comment, &istat);
  snfitsio_errorCheck("can't read NPRIM", istat);
  printf("\t\t Read %-*s  = %d  primary refs \n", 
	 NUMPRINT, KEYWORD, *IPTR );

  // read name of each primary
  for(i=0; i < CALIB_INFO.NPRIMARY; i++ ) {
    sprintf(KEYWORD,"PRIM%3.3d", i+1); 
    fits_read_key(FP, TSTRING, KEYWORD, tmpStr, comment, &istat);
    // extract first word 
    SPTR = strtok(tmpStr," "); 
    CALIB_INFO.PRIMARY_NAME[i] = (char*)malloc(MEMC);
    sprintf(CALIB_INFO.PRIMARY_NAME[i],"%s",SPTR);

    sprintf(c1err,"can't read %s", KEYWORD);
    snfitsio_errorCheck(c1err, istat);
    printf("\t\t Read %-*s  = %s \n", NUMPRINT, KEYWORD, SPTR );
  }

  // read NFILTERS
  sprintf(KEYWORD,"NFILTERS");  IPTR = &CALIB_INFO.NFILTDEF;
  fits_read_key(FP, TINT, KEYWORD, IPTR, comment, &istat);
  snfitsio_errorCheck("can't read NFILTERS", istat);
  printf("\t\t Read %-*s  = %d  filters \n", NUMPRINT, KEYWORD, *IPTR );
  
  // read name of each filter
  for(i=0; i < CALIB_INFO.NFILTDEF; i++ ) {
    sprintf(KEYWORD,"FILT%3.3d", i+1);
    CALIB_INFO.FILTER_NAME[i] = (char*)malloc(MEMC);
    CALIB_INFO.SURVEY_NAME[i] = (char*)malloc(MEMC);
    SPTR=CALIB_INFO.FILTER_NAME[i];
    fits_read_key(FP, TSTRING, KEYWORD, SPTR, comment, &istat);

    sprintf(c1err,"can't read %s", KEYWORD);
    snfitsio_errorCheck(c1err, istat);
    //    printf("\t\t Read %-*s  = '%s'  comment=%s\n", 
    //	   NUMPRINT, KEYWORD, SPTR, comment );

    if ( SPTR[0] == '*' ) { CALIB_INFO.ISLAMSHIFT[i] = true; }

    // store absolute filter index 
    len = strlen(CALIB_INFO.FILTER_NAME[i]);
    sprintf(cfilt, "%c", CALIB_INFO.FILTER_NAME[i][len-1] ) ;
    IFILTDEF = INTFILTER(cfilt) ;
    CALIB_INFO.IFILTDEF[i] = IFILTDEF ;
	    
    // mark survey filters
    if ( strchr(CALIB_INFO.FILTERS_SURVEY,cfilt[0]) != NULL )  
      { CALIB_INFO.IS_SURVEY_FILTER[IFILTDEF] = true;  }

    // Nov 2022 
    // comment has the form "name: SURVEY=[SURVEY]"
    // so strip off SURVEY name per filter and store it.
    char *eq    = strchr(comment, '=');
    int  ind_eq = (int)(eq - comment);
    sprintf(CALIB_INFO.SURVEY_NAME[i], "%s", &comment[ind_eq+1]);

  }  // end NFILTDEF loop

  int NFLAMSHIFT = CALIB_INFO.NFLAMSHIFT ;
  int NFILTDEF   = CALIB_INFO.NFILTDEF ;

  // - - - - - - 

  // check lam-shifted filters
  if ( NFLAMSHIFT > 0 ) {
    if ( NFLAMSHIFT == NFILTDEF/2 ) {
      CALIB_INFO.ISLAMSHIFT[NFILTDEF] = true;
      NFILTDEF = CALIB_INFO.NFILTDEF  = NFLAMSHIFT;
    }
    else {
      sprintf(c1err,"NFLAMSHIFT=%d != (NFILTDEF/2 =%d/2",
	      NFLAMSHIFT, NFILTDEF);
      sprintf(c2err,"Check LAM-shifted bands in kcor file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }
  }

  // - - - - - - - - 
  // read optional Galactic RV and color law used to warp spectra for k-corr.
  // Note that there is no abort on missing key.

  sprintf(KEYWORD,"RV");  DPTR = &CALIB_INFO.RVMW ;
  fits_read_key(FP, TDOUBLE, KEYWORD, DPTR, comment, &istat);
  if ( istat == 0 ) 
    { printf("\t\t Read %-*s  = %4.2f \n",   NUMPRINT, KEYWORD, *DPTR ); }
  istat = 0 ;

  sprintf(KEYWORD,"OPT_MWCOLORLAW");  IPTR = &CALIB_INFO.OPT_MWCOLORLAW ;
  fits_read_key(FP, TINT, KEYWORD, IPTR, comment, &istat);
  if ( istat == 0 ) 
    { printf("\t\t Read %-*s  = %d \n",   NUMPRINT, KEYWORD, *IPTR ); }
  istat = 0 ;
  
  // read number of KCOR tables
  sprintf(KEYWORD,"NKCOR");  IPTR = &CALIB_INFO.NKCOR ;
  fits_read_key(FP, TINT, KEYWORD, IPTR, comment, &istat);
  snfitsio_errorCheck("can't read NKCOR", istat);
  printf("\t\t Read %-*s  = %d  K-COR tables \n", NUMPRINT, KEYWORD, *IPTR );
  
  int NKCOR = CALIB_INFO.NKCOR;
  if ( NKCOR > MXTABLE_KCOR ) {
    sprintf(c1err,"NKCOR=%d exceeds bound of MXTABLE_KCOR=%d",
	    NKCOR, MXTABLE_KCOR );
    sprintf(c2err,"Check NKCOR in FITS header");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
  }

  // kcor strings that define rest- and obs-frame bands
  for(i=0; i < NKCOR; i++ ) {
    sprintf(KEYWORD,"KCOR%3.3d", i+1);
    CALIB_INFO.STRING_KCORLINE[i] = (char*)malloc(MEMC);
    SPTR = CALIB_INFO.STRING_KCORLINE[i];
    fits_read_key(FP, TSTRING, KEYWORD, SPTR, comment, &istat);

    sprintf(c1err,"can't read %s", KEYWORD);
    snfitsio_errorCheck(c1err, istat);
    // printf("\t Read %-*s = %s \n", NUMPRINT, KEYWORD, SPTR );
  }
  

  // - - - - - - - 
  // read bin info 

  read_kcor_binInfo("wavelength", "L"  , MXLAMBIN_FILT_CALIB,
		    &CALIB_INFO.BININFO_LAM ); 
  read_kcor_binInfo("epoch",      "T"  , MXTBIN_KCOR,
		    &CALIB_INFO.BININFO_T ); 
  read_kcor_binInfo("redshift",   "Z"  , MXZBIN_KCOR,
		    &CALIB_INFO.BININFO_z ); 
  read_kcor_binInfo("AV",         "AV" , MXAVBIN_KCOR,
		    &CALIB_INFO.BININFO_AV ); 

  CALIB_INFO.zRANGE_LOOKUP[0] = CALIB_INFO.BININFO_z.RANGE[0] ;
  CALIB_INFO.zRANGE_LOOKUP[1] = CALIB_INFO.BININFO_z.RANGE[1] ;

  return ;

} // end read_calib_head


// ==================================================
void read_kcor_binInfo(char *VARNAME, char *VARSYM, int MXBIN,
		       KCOR_BININFO_DEF *BININFO) {

  // Read the following from header and load BININFO struct with
  //   Number of bins  :  NB[VARSYM]
  //   binsize         :  [VARSYM]BIN
  //   min val         :  [VARSYM]MIN
  //   max val         :  [VARSYM]MAX
  //
  // Abort if number of bins exceeds MXBIN

  fitsfile *FP = CALIB_INFO.FP;
  int istat = 0, i, NBIN ;
  char KEYWORD[20], comment[100];
  int *IPTR; double *DPTR;
  char fnam[] = "read_kcor_binInfo" ;

  // ----------- BEGIN -----------

  sprintf(BININFO->VARNAME, "%s", VARNAME);

  sprintf(KEYWORD,"NB%s", VARSYM);  IPTR = &BININFO->NBIN;
  fits_read_key(FP, TINT, KEYWORD, IPTR, comment, &istat);
  sprintf(c1err,"can't read %s", KEYWORD);
  snfitsio_errorCheck(c1err, istat);

  NBIN = BININFO->NBIN ;

  sprintf(KEYWORD,"%sBIN", VARSYM);  DPTR = &BININFO->BINSIZE;
  fits_read_key(FP, TDOUBLE, KEYWORD, DPTR, comment, &istat);
  sprintf(c1err,"can't read %s", KEYWORD);
  snfitsio_errorCheck(c1err, istat);

  sprintf(KEYWORD,"%sMIN", VARSYM);  DPTR = &BININFO->RANGE[0];
  fits_read_key(FP, TDOUBLE, KEYWORD, DPTR, comment, &istat);
  sprintf(c1err,"can't read %s", KEYWORD);
  snfitsio_errorCheck(c1err, istat);

  sprintf(KEYWORD,"%sMAX", VARSYM);  DPTR = &BININFO->RANGE[1];
  fits_read_key(FP, TDOUBLE, KEYWORD, DPTR, comment, &istat);
  sprintf(c1err,"can't read %s", KEYWORD);
  snfitsio_errorCheck(c1err, istat);

  printf("\t\t Read %4d %-10s bins (%.2f to %.2f)\n",
	 NBIN, VARNAME, BININFO->RANGE[0], BININFO->RANGE[1]);
  fflush(stdout);


  if ( BININFO->NBIN >= MXBIN ) {
    sprintf(c1err, "NBIN(%s) = %d exceeds bounf of %d",
	    VARNAME, NBIN, MXBIN );
    sprintf(c2err, "Check FITS HEADER keys: NB%s, %sBIN, %sMIN, %sMAX",
	    VARSYM, VARSYM, VARSYM, VARSYM );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
  }

  // store value in each bin
  double xi, xmin=BININFO->RANGE[0], xbin=BININFO->BINSIZE ;
  BININFO->GRIDVAL = (double*) malloc(NBIN * sizeof(double) );
  for(i=0; i < NBIN; i++ ) {
    xi = (double)i ;
    BININFO->GRIDVAL[i] = xmin + ( xi * xbin ) ;
  }

  return ;

} // end read_kcor_binInfo

// =============================
void read_calib_zpoff(void) {

  // read filter-dependent info:
  // PrimaryName  PrimaryMag  ZPTOFF(prim)  ZPTOFF(filt)

  fitsfile *FP = CALIB_INFO.FP ;
  int  NFILTDEF = CALIB_INFO.NFILTDEF ;
  int  NPRIMARY = CALIB_INFO.NPRIMARY ; 

  int  hdutype, istat=0, ifilt, anynul, iprim, iprim_store ;
  char **NAME_PRIM, *tmpName ;
  char fnam[] = "read_calib_zpoff" ;

  //  int ICOL_FILTER_NAME        = 1 ;
  int ICOL_PRIMARY_NAME       = 2 ;
  int ICOL_PRIMARY_MAG        = 3 ;
  int ICOL_PRIMARY_ZPOFF_SYN  = 4 ;
  int ICOL_PRIMARY_ZPOFF_FILE = 5 ;

  long FIRSTROW, NROW, FIRSTELEM = 1;
  
  // --------- BEGIN ----------

  printf("   %s  \n", fnam); fflush(stdout);

  fits_movrel_hdu(FP, 1, &hdutype, &istat);
  snfitsio_errorCheck("Cannot move to ZPOFF table", istat);

  NAME_PRIM = (char**) malloc(2*sizeof(char*) );
  NAME_PRIM[0] = (char*) malloc(60*sizeof(char));
  NAME_PRIM[1] = (char*) malloc(60*sizeof(char));
    
  for(ifilt=0; ifilt < NFILTDEF; ifilt++ ) {
           
    FIRSTROW = ifilt + 1;  NROW=1;
    fits_read_col_str(FP, ICOL_PRIMARY_NAME, FIRSTROW, FIRSTELEM, NROW,
		      NULL_A, NAME_PRIM, &anynul, &istat )  ;      

    // find IPRIM index
    iprim_store = -9;
    for(iprim=0; iprim < NPRIMARY; iprim++ ) {
      tmpName = CALIB_INFO.PRIMARY_NAME[iprim] ;
      if ( strcmp(NAME_PRIM[0],tmpName) == 0 ) { iprim_store = iprim; }      
    }

    if( iprim_store < 0 ) {
      sprintf(c1err,"Unrecognized PRIMARY_NAME = %s", NAME_PRIM[0] );
      sprintf(c2err,"Check ZPT table");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }

    CALIB_INFO.PRIMARY_INDX[ifilt] = iprim_store ;   
  } // end ifilt


  FIRSTROW = 1;  NROW = NFILTDEF ;

  // read primary mag, ZPOFF for primary, and photometry offsets
  // passed from optional ZPOFF.DAT file in $SNDATA_ROOT/filters.

  fits_read_col_dbl(FP, ICOL_PRIMARY_MAG, FIRSTROW, FIRSTELEM, NROW,
		    NULL_1D, CALIB_INFO.PRIMARY_MAG, &anynul, &istat )  ;      
  snfitsio_errorCheck("Read PRIMARY_MAG", istat);

  fits_read_col_dbl(FP, ICOL_PRIMARY_ZPOFF_SYN, FIRSTROW, FIRSTELEM, NROW,
		    NULL_1D, CALIB_INFO.PRIMARY_ZPOFF_SYN, &anynul, &istat );
  snfitsio_errorCheck("Read PRIMARY_ZPOFF", istat);

  
  // read optional ZPOFF from ZPOFF.DAT file in filter subDir.
  // This is typoically a post-publication hack to get mags
  // back on the desired system.
  fits_read_col_dbl(FP, ICOL_PRIMARY_ZPOFF_FILE, FIRSTROW, FIRSTELEM, NROW,
		    NULL_1D, CALIB_INFO.PRIMARY_ZPOFF_FILE, &anynul, &istat );  
  snfitsio_errorCheck("Read PRIMARY_ZPOFF_FILE", istat);

  if ( KCOR_VERBOSE_FLAG  ) {
    printf("\n");
    printf(" xxx  %s DUMP: \n", fnam);
    printf(" xxx                   Prim.   Prim.   Prim.    Filter \n");
    printf(" xxx Filter            name    Mag     ZPOFF    ZPOFF \n");
    printf(" xxx ----------------------------------------------------- \n");
    for(ifilt=0; ifilt < NFILTDEF; ifilt++ ) {

      iprim = CALIB_INFO.PRIMARY_INDX[ifilt];

      printf(" xxx %-14s %6s(%d)  %6.3f  %7.4f  %7.4f \n",
	     CALIB_INFO.FILTER_NAME[ifilt],
	     CALIB_INFO.PRIMARY_NAME[iprim], CALIB_INFO.PRIMARY_INDX[ifilt],
	     CALIB_INFO.PRIMARY_MAG[ifilt],
	     CALIB_INFO.PRIMARY_ZPOFF_SYN[ifilt],
	     CALIB_INFO.PRIMARY_ZPOFF_FILE[ifilt] );
	     
    }
    printf("\n");
  } // end verbose

  return ;

} // end read_calib_zpoff

// =============================
void read_calib_snsed(void) {

  fitsfile *FP  = CALIB_INFO.FP ;
  int  NBL      = CALIB_INFO.BININFO_LAM.NBIN;
  int  NBT      = CALIB_INFO.BININFO_T.NBIN;
  long NROW     = NBL * NBT ;
  int  MEMF     = NROW*sizeof(float);
  long FIRSTROW = 1, FIRSTELEM=1 ;
  int  ICOL=1, istat=0, hdutype, anynul ;
  char fnam[] = "read_calib_snsed" ;
  // --------- BEGIN ----------
  printf("   %s \n", fnam); fflush(stdout);

  CALIB_INFO.FLUX_SNSED_F = (float*) malloc(MEMF);

  fits_movrel_hdu(FP, 1, &hdutype, &istat);
  snfitsio_errorCheck("Cannot move to SNSED table", istat);

  fits_read_col_flt(FP, ICOL, FIRSTROW, FIRSTELEM, NROW,
		    NULL_1E, CALIB_INFO.FLUX_SNSED_F, &anynul, &istat )  ;      
  snfitsio_errorCheck("Read FLUX_SNSED", istat);

  return ;

} // end read_calib_snsed

// =============================
void read_kcor_tables(void) {

  // Examine K_xy to flag which filters are OBS and REST.
  // If there are no K-cor tables, then any SURVEY filter
  // is defined as an OBS filter. 

  fitsfile *FP      = CALIB_INFO.FP ;
  int NKCOR         = CALIB_INFO.NKCOR;
  int NFILTDEF_KCOR = CALIB_INFO.NFILTDEF;

  int NKCOR_STORE = 0 ;
  int k, i, i2, icol, hdutype, istat=0, anynul, ifilt_rest, ifilt_obs, len ;
  int IFILTDEF, ifilt, LBX0, LBX1;
  char *STRING_KLINE, *STRING_KSYM, cfilt_rest[40], cfilt_obs[40];
  char cband_rest[2], cband_obs[2], *FILTER_NAME ;
  char fnam[] = "read_kcor_tables" ;

  // --------- BEGIN ----------

  printf("   %s \n", fnam); fflush(stdout);

  fits_movrel_hdu(FP, 1, &hdutype, &istat);
  snfitsio_errorCheck("Cannot move to KCOR table", istat);

  if ( NKCOR == 0 ) { return; }

  for(k=0; k < NKCOR; k++ ) {

    CALIB_INFO.STRING_KCORSYM[k] = (char*)malloc(8*sizeof(char) ); 
    STRING_KLINE = CALIB_INFO.STRING_KCORLINE[k] ;
    STRING_KSYM  = CALIB_INFO.STRING_KCORSYM[k] ; // e.g., K_xy

    parse_KCOR_STRING(STRING_KLINE, STRING_KSYM, cfilt_rest, cfilt_obs);
    ifilt_rest = INTFILTER(cfilt_rest);
    ifilt_obs  = INTFILTER(cfilt_obs);

    // if this IFILT_OBS is not a SURVEY filter, then bail.
    // i.e., ignore K_xy that include extra unused filters.

    len = strlen(cfilt_rest); sprintf(cband_rest,"%c", cfilt_rest[len-1]);
    len = strlen(cfilt_obs ); sprintf(cband_obs, "%c", cfilt_obs[len-1]);
    
    if ( strchr(CALIB_INFO.FILTERS_SURVEY,cband_obs[0]) == NULL )  
      { continue; }

    CALIB_INFO.EXIST_KCOR[ifilt_rest][ifilt_obs] = true ;
    
    // set mask for rest frame and/or observer frame
    for(ifilt=0; ifilt < NFILTDEF_KCOR; ifilt++ ) {
      IFILTDEF    = CALIB_INFO.IFILTDEF[ifilt];
      FILTER_NAME = CALIB_INFO.FILTER_NAME[ifilt];

      if ( strcmp(FILTER_NAME,cfilt_rest) == 0 )
	{ CALIB_INFO.MASK_FRAME_FILTER[ifilt] |= MASK_FRAME_REST ; }

      if ( strcmp(FILTER_NAME,cfilt_obs) == 0 )
	{ CALIB_INFO.MASK_FRAME_FILTER[ifilt] |= MASK_FRAME_OBS ; }

    } // end ifilt

    LBX0 = ISBXFILT_KCOR(cfilt_rest);
    if ( LBX0 ) { CALIB_INFO.MASK_EXIST_BXFILT |= MASK_FRAME_REST; }

    LBX1 = ISBXFILT_KCOR(cfilt_obs);
    if ( LBX1 ) { CALIB_INFO.MASK_EXIST_BXFILT |= MASK_FRAME_OBS; }


    /* can't remember purpose of STANDALONE mode ... fix later 
         IF ( RDKCOR_STANDALONE .and. 
     &         IFILTDEF_INVMAP_SURVEY(ifilt_obs) .LE. 0 ) THEN
            NFILTDEF_SURVEY = NFILTDEF_SURVEY + 1
            IFILTDEF_MAP_SURVEY(NFILTDEF_SURVEY) = IFILT_OBS
            IFILTDEF_INVMAP_SURVEY(ifilt_obs)    = NFILTDEF_SURVEY   
         ENDIF
    */

    // abort on undefined filter

    if ( ifilt_rest <= 0 || ifilt_obs <= 0 ) {
      sprintf(c1err,"Undefined %s: IFILTDEF(REST,OBS)=%d,%d",
	      STRING_KSYM, ifilt_rest, ifilt_obs);
      sprintf(c2err,"Check filters in Kcor file.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // define new filter only if not already defined.
    addFilter_kcor(ifilt_rest, cfilt_rest, &CALIB_INFO.FILTERCAL_REST);
    addFilter_kcor(ifilt_obs,  cfilt_obs,  &CALIB_INFO.FILTERCAL_OBS );

    // ??? if ( EXIST_BXFILT_OBS .and. RDKCOR_STANDALONE ) { continue; }   

    CALIB_INFO.IFILTMAP_KCOR[OPT_FRAME_REST][NKCOR_STORE] = ifilt_rest; 
    CALIB_INFO.IFILTMAP_KCOR[OPT_FRAME_OBS][NKCOR_STORE]  = ifilt_obs ;
    CALIB_INFO.k_index[NKCOR_STORE] = k;
    NKCOR_STORE++; 
    CALIB_INFO.NKCOR_STORE = NKCOR_STORE ;

  } // end k loop over KCOR tables


  // sanity check on number of rest,obs filters
  int NFILT_REST = CALIB_INFO.FILTERCAL_REST.NFILTDEF ;
  int NFILT_OBS  = CALIB_INFO.FILTERCAL_OBS.NFILTDEF ;
  if ( NFILT_REST >= MXFILT_REST_CALIB ) {
    sprintf(c1err,"NFILT_REST = %d exceeds bound of %d.",
	    NFILT_REST, MXFILT_REST_CALIB);
    sprintf(c2err,"Check kcor input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  if ( NFILT_OBS >= MXFILT_OBS_CALIB ) {
    sprintf(c1err,"NFILT_OBS = %d exceeds bound of %d.",
	    NFILT_OBS, MXFILT_OBS_CALIB);
    sprintf(c2err,"Check kcor input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // ------------------------------------------------------
  // BX check:
  // loop over defined filters (from header) ; if undefined rest-frame 
  // X filter exists, then add it to the rest-frame list. This allows 
  // Landolt option in fitting program without having to explicitly
  // define a K-correction wit the X filter.
  // Beware to set BX before INIT_KCOR_INDICES !!!


  for(ifilt=0; ifilt < NFILTDEF_KCOR ; ifilt++ ) {
    FILTER_NAME  = CALIB_INFO.FILTER_NAME[ifilt];
    IFILTDEF = CALIB_INFO.IFILTDEF[ifilt];
    if( ISBXFILT_KCOR(FILTER_NAME) ) {  // ensure 'BX', not BLABLA-X 
      CALIB_INFO.MASK_EXIST_BXFILT        |= MASK_FRAME_REST ; 
      CALIB_INFO.MASK_FRAME_FILTER[ifilt] |= MASK_FRAME_REST ;
      addFilter_kcor(IFILTDEF, FILTER_NAME, &CALIB_INFO.FILTERCAL_REST);
    }
  }


  // init multi-dimensional array to store KCOR tables
  init_kcor_indices();


  /* 
c =============================
c if user has specified any MAGOBS_SHIFT_PRIMARY, MAGOBS_SHIFT_ZP,
c (and same for REST) for a non-existant filter, then ABORT ...
c 
c   WARNING: this check stays in snana.car because it's
c            based on user input to &SNLCINP

      CALL  RDKCOR_CHECK_MAGSHIFTS
  */

  // read the actual KCOR table(s)
  long FIRSTROW = 1, FIRSTELEM=1 ;
  int NBINTOT = CALIB_INFO.MAPINFO_KCOR.NBINTOT;
  int NBT     = CALIB_INFO.MAPINFO_KCOR.NBIN[KDIM_T];
  int NBz     = CALIB_INFO.MAPINFO_KCOR.NBIN[KDIM_z];
  int NBAV    = CALIB_INFO.MAPINFO_KCOR.NBIN[KDIM_AV];
  int NROW    = NBT * NBz * NBAV;

  int ifilto, ifiltr, IBKCOR[NKDIM_KCOR], IBIN_FIRST, IBIN_LAST;;
  double KCOR_SHIFT;

  int MEMF    = NBINTOT * sizeof(float);
  CALIB_INFO.KCORTABLE1D_F = (float*) malloc(MEMF);
  for(k=0; k < NBINTOT; k++ ) { CALIB_INFO.KCORTABLE1D_F[k] = 9999.9; }

  // loop over kcor tables
  for(i=0; i < NKCOR_STORE; i++ ) {
    ifilt_obs  = CALIB_INFO.IFILTMAP_KCOR[OPT_FRAME_OBS][i];
    ifilt_rest = CALIB_INFO.IFILTMAP_KCOR[OPT_FRAME_REST][i];
    k          = CALIB_INFO.k_index[i]; // original index in FITS file
    icol       = k + 4; // skip T,z,AV columns
    STRING_KLINE = CALIB_INFO.STRING_KCORLINE[k] ; 
    STRING_KSYM  = CALIB_INFO.STRING_KCORSYM[k] ; 
    
    //     get sparse filter indices
    ifilto = CALIB_INFO.FILTERCAL_OBS.IFILTDEF_INV[ifilt_obs]; 
    ifiltr = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifilt_rest];

    if ( ifilto < 0 ) {
      sprintf(c1err, "Unknown sparse index for ifilt_obs=%d", ifilt_obs);
      sprintf(c2err, "Check '%s' ", STRING_KSYM);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    if ( ifiltr < 0 ) {
      sprintf(c1err,"Unknown sparse index for ifilt_rest=%d", ifilt_rest);
      sprintf(c2err,"Check '%s' ", STRING_KSYM);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    IBKCOR[KDIM_IFILTr] = ifiltr ;
    IBKCOR[KDIM_IFILTo] = ifilto ;	  
    IBKCOR[KDIM_T]      = 0 ;     
    IBKCOR[KDIM_z]      = 0 ;
    IBKCOR[KDIM_AV]     = 0 ;    
    IBIN_FIRST  = get_1DINDEX(IDMAP_KCOR_TABLE, NKDIM_KCOR, IBKCOR);
    IBIN_LAST   = IBIN_FIRST + NROW - 1 ;

    fits_read_col_flt(FP, icol, FIRSTROW, FIRSTELEM, NROW,
		      NULL_1E, &CALIB_INFO.KCORTABLE1D_F[IBIN_FIRST], 
		      &anynul, &istat )  ;      
    snfitsio_errorCheck("Read KCOR TABLE", istat);

    // apply user-defined parimary mag shifts  (e.g., systematic tests)
    KCOR_SHIFT = 
      CALIB_INFO.MAGOBS_SHIFT_PRIMARY[ifilt_obs] -
      CALIB_INFO.MAGREST_SHIFT_PRIMARY[ifilt_rest] ;

    for(i2=IBIN_FIRST; i2 <= IBIN_LAST; i2++ ) 
      { CALIB_INFO.KCORTABLE1D_F[i2] += (float)KCOR_SHIFT;  }

  } // end i-loop to NKCOR_STORE

  return ;

} // end read_kcor_tables

// ===========================
int ISBXFILT_KCOR(char *cfilt) {
  // return true if BX is part of filter name
  if ( strstr(cfilt,"BX") != NULL ) { return(1); }
  return(0);
} // end ISBXFILT_KCOR

// ===========================================
void addFilter_kcor(int ifiltdef, char *FILTER_NAME, FILTERCAL_DEF *MAP){

  // ifiltdev = 0     --> zero map, return
  // ifiltdef = 777   --> dump map
  // ifiltdef = 1 - N --> load map
  //
  // Inputs:
  //    ifiltdef : absolute filter index
  //    NAME     : full name of filter
  //    MAP      : structure to malloc/load

  int OPT_FRAME = MAP->OPT_FRAME; // indicates REST or OBS
  int ifilt, NF;
  char cfilt1[2], *SURVEY_NAME ;
  char fnam[] = "addFilter_kcor" ;

  // --------- BEGIN -----------

  if ( ifiltdef == 0 ) {
    // zero map, then return
    MAP->NFILTDEF = 0;
    MAP->FILTERSTRING[0] =  0 ; 
    MAP->NFILT_DUPLICATE =  0 ;
    for(ifilt=0; ifilt < MXFILT_CALIB; ifilt++ ) {
      MAP->IFILTDEF[ifilt]       = -9 ; 
      MAP->IFILTDEF_INV[ifilt]   = -9 ;
      MAP->FILTER_NAME[ifilt]    = (char*)malloc(40*sizeof(char) ) ;
      MAP->FILTER_NAME[ifilt][0] = 0;
      MAP->SURVEY_NAME[ifilt]    = (char*)malloc(80*sizeof(char) ) ;
      MAP->SURVEY_NAME[ifilt][0] = 0;
      MAP->NDEFINE[ifilt] = 0 ;
      MAP->PRIMARY_MAG[ifilt]   = 99.0 ;
      MAP->PRIMARY_ZPOFF_SYN[ifilt]  =  0.0 ;  // required
      MAP->PRIMARY_ZPOFF_FILE[ifilt] =  0.0 ;  // optional
      MAP->PRIMARY_KINDX[ifilt] = -9 ;
      MAP->NBIN_LAM[ifilt]  =  0 ;
    }
    return ;
  }


  if ( ifiltdef == 777 ) {
    // dump map, then return
    int IFILTDEF ;
    NF = MAP->NFILTDEF;
    printf("\n");
    printf("\t xxx %s dump: \n", fnam);
    printf("\t xxx FILTERSTRING = '%s' \n", MAP->FILTERSTRING);
    for(ifilt=0; ifilt < NF; ifilt++ ) {
      IFILTDEF = MAP->IFILTDEF[ifilt];
      sprintf(cfilt1, "%c", FILTERSTRING[IFILTDEF] );

      /*
      printf("\t xxx IFILTDEF[%2d,%s] = %2d  (%s)  PRIM_MAG=%.3f NBL=%d\n",
	     ifilt, cfilt1, IFILTDEF, MAP->FILTER_NAME[ifilt],
	     MAP->PRIMARY_MAG[ifilt], MAP->NBIN_LAM[ifilt]  ); 
	     fflush(stdout); */
    }
    return ;
  }

  // return if this filter is already defined
  if ( MAP->IFILTDEF_INV[ifiltdef] >= 0 ) { return; }

  NF = MAP->NFILTDEF ;
  sprintf(cfilt1, "%c", FILTERSTRING[ifiltdef] );

  MAP->IFILTDEF_INV[ifiltdef] = NF;
  MAP->IFILTDEF[NF]           = ifiltdef;
  strcat(MAP->FILTERSTRING,cfilt1);
  sprintf(MAP->FILTER_NAME[NF], "%s", FILTER_NAME);


  // find FILTER_NAME in CALIB_INFO and store assocated SURVEY_NAME
  for (ifilt=0; ifilt < CALIB_INFO.NFILTDEF; ifilt++ ) {
    if ( strcmp(FILTER_NAME,CALIB_INFO.FILTER_NAME[ifilt]) == 0 ) {
      SURVEY_NAME = CALIB_INFO.SURVEY_NAME[ifilt] ;
      sprintf(MAP->SURVEY_NAME[NF],"%s", SURVEY_NAME);
    }
  }
  
  MAP->NFILTDEF++ ;


  // find original filter index from header to get primary mag & zpoff
  int k, kfilt=-9, NFILTDEF_KCOR = CALIB_INFO.NFILTDEF ;
  int  MASK_FRAME;
  bool MATCH_FRAME, MATCH_FILTER;
  char *NAME;
  for(k=0; k < NFILTDEF_KCOR; k++ ) {
    MATCH_FILTER = ( CALIB_INFO.IFILTDEF[k] == ifiltdef );
    if ( !MATCH_FILTER ) { continue; }
    
    MASK_FRAME   = CALIB_INFO.MASK_FRAME_FILTER[k] ; // rest or obs
    MATCH_FRAME  = ( MASK_FRAME & (1<<OPT_FRAME) ) > 0 ; 

    /* xxx
    printf(" xxx %s: k=%2d ifiltdef=%2d(%s) " 
	   "MASK_FRAME=%2d OPT_FR=%d  MATCH(FRAME,FILT)=%d,%d \n",
	   fnam, k, ifiltdef, cfilt1, MASK_FRAME, OPT_FRAME,
	   MATCH_FRAME, MATCH_FILTER); xxx */
    if ( MATCH_FILTER && MATCH_FRAME )  { kfilt = k; }
  }

  if ( kfilt < 0 ) {
    sprintf(c1err,"Could not find kfilt for ifiltdef=%d (%s)", 
	    ifiltdef, FILTER_NAME);
    sprintf(c2err,"Probably a code bug.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // store index of primary to read later
  MAP->PRIMARY_KINDX[NF] = CALIB_INFO.PRIMARY_INDX[kfilt] ;

  double *ptr_SHIFT;
  if ( OPT_FRAME == OPT_FRAME_REST ) 
    { ptr_SHIFT = CALIB_INFO.MAGREST_SHIFT_PRIMARY; }
  else
    { ptr_SHIFT = CALIB_INFO.MAGOBS_SHIFT_PRIMARY; }

  MAP->PRIMARY_MAG[NF]  = 
    CALIB_INFO.PRIMARY_MAG[kfilt] + ptr_SHIFT[ifiltdef] ;

  MAP->PRIMARY_ZPOFF_SYN[NF] = 
    CALIB_INFO.PRIMARY_ZPOFF_SYN[kfilt] + ptr_SHIFT[ifiltdef] ;

  MAP->PRIMARY_ZPOFF_FILE[NF] = 
    CALIB_INFO.PRIMARY_ZPOFF_FILE[kfilt];

  return ;

} // end  addFilter_kcor


// ==============================
void parse_KCOR_STRING(char *STRING, 
		       char *strKcor, char *cfilt_rest, char *cfilt_obs) {

  // Parse input STRING of the form:
  //   Kcor K_XY for rest [filterName]X to obs [filterName]Y
  //
  // and return
  //   strKcor    = 'K_XY'
  //   cfilt_rest = [filterName]X
  //   cfilt_obs  = [filterName]Y
  //

  int NWD, iwd, len ;
  char word[40], cband_rest[2], cband_obs[2];
  char fnam[] = "parse_KCOR_STRING" ;

  // ---------- BEGIN ----------

  strKcor[0] = cfilt_rest[0] = cfilt_obs[0] = 0;

  NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING, STRING);

  for(iwd=0; iwd < NWD-1; iwd++ ) {

    get_PARSE_WORD(0, iwd, word );

    if( word[0] == 'K' && word[1] == '_' ) 
      { sprintf(strKcor, "%s", word); }

    if ( strcmp(word,"rest") == 0 ) 
      { get_PARSE_WORD(0, iwd+1, cfilt_rest); }

    if ( strcmp(word,"obs") == 0 ) 
      { get_PARSE_WORD(0, iwd+1, cfilt_obs); }
  } // end iwd loop

  // - - - - - -
  // sanity tests

  int NERR = 0;
  if ( strlen(strKcor) == 0 ) {
    printf(" ERROR: couldn't find required K_XY in KCOR_STRING\n");
    NERR++ ;  fflush(stdout);
  }

  if ( strlen(cfilt_rest) == 0 ) {
    printf(" ERROR: couldn't find rest-frame filter in KCOR_STRING\n");
    NERR++ ;  fflush(stdout);
  }

  if ( strlen(cfilt_obs) == 0 ) {
    printf(" ERROR: couldn't find obs-frame filter in KCOR_STRING\n");
    NERR++ ;  fflush(stdout);
  }

  cband_rest[0] = strKcor[2]; // rest-frame band char
  cband_obs[0]  = strKcor[3]; // obs-frame band char

  len = strlen(cfilt_rest);
  if ( cband_rest[0] != cfilt_rest[len-1] ) {
    printf(" ERROR: rest band='%c' is not compatible with rest filter='%s'\n",
	   cband_rest[0], cfilt_rest);
    NERR++ ;  fflush(stdout);
  }

  len = strlen(cfilt_obs);
  if ( cband_obs[0] != cfilt_obs[len-1] ) {
    printf(" ERROR: obs band='%c' is not compatible with obs filter='%s'\n",
	   cband_obs[0], cfilt_obs );
    NERR++ ;  fflush(stdout);
  }

  if ( NERR > 0 ) {
    sprintf(c1err,"Problem parsing KCOR STRING: ");
    sprintf(c2err,"%s", STRING);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return ;
} // end parse_KCOR_STRING


// ======================================
void init_kcor_indices(void) {

  // Init index maps for KCOR-related tables:
  //   KCOR, AVWARP, LCMAP and MWXT.
  // These tables are multi-dimensional, but we use a
  // 1-d table allocation to save memory. This routine
  // prepares index-mappings from the multi-dimensional
  // space to the 1-d space.

  int  NBIN_T          = CALIB_INFO.BININFO_T.NBIN ;
  int  NBIN_z          = CALIB_INFO.BININFO_z.NBIN ;
  int  NBIN_AV         = CALIB_INFO.BININFO_AV.NBIN ;
  int  NFILTDEF_REST   = CALIB_INFO.FILTERCAL_REST.NFILTDEF ;
  int  NFILTDEF_OBS    = CALIB_INFO.FILTERCAL_OBS.NFILTDEF ;
    char fnam[] = "init_kcor_indices";

  // -------------- BEGIN ----------------

  if ( CALIB_INFO.NKCOR_STORE == 0 ) { return; }
 
  sprintf(CALIB_INFO.MAPINFO_KCOR.NAME, "KCOR");
  CALIB_INFO.MAPINFO_KCOR.IDMAP             = IDMAP_KCOR_TABLE ;
  CALIB_INFO.MAPINFO_KCOR.NDIM              = NKDIM_KCOR;
  CALIB_INFO.MAPINFO_KCOR.NBIN[KDIM_T]      = NBIN_T ;
  CALIB_INFO.MAPINFO_KCOR.NBIN[KDIM_z]      = NBIN_z ;
  CALIB_INFO.MAPINFO_KCOR.NBIN[KDIM_AV]     = NBIN_AV ;
  CALIB_INFO.MAPINFO_KCOR.NBIN[KDIM_IFILTr] = NFILTDEF_REST ;
  CALIB_INFO.MAPINFO_KCOR.NBIN[KDIM_IFILTo] = NFILTDEF_OBS ;
  get_MAPINFO_KCOR("NBINTOT", &CALIB_INFO.MAPINFO_KCOR);

  sprintf(CALIB_INFO.MAPINFO_AVWARP.NAME, "AVWARP");
  CALIB_INFO.MAPINFO_AVWARP.IDMAP             = IDMAP_KCOR_AVWARP ;
  CALIB_INFO.MAPINFO_AVWARP.NDIM              = N4DIM_KCOR;
  CALIB_INFO.MAPINFO_AVWARP.NBIN[0]           = NBIN_T ;
  CALIB_INFO.MAPINFO_AVWARP.NBIN[1]           = MXCBIN_AVWARP ;
  CALIB_INFO.MAPINFO_AVWARP.NBIN[2]           = NFILTDEF_REST ;
  CALIB_INFO.MAPINFO_AVWARP.NBIN[3]           = NFILTDEF_REST ;
  get_MAPINFO_KCOR("NBINTOT", &CALIB_INFO.MAPINFO_AVWARP);

  sprintf(CALIB_INFO.MAPINFO_LCMAG.NAME, "LCMAG");
  CALIB_INFO.MAPINFO_LCMAG.IDMAP             = IDMAP_KCOR_LCMAG ;
  CALIB_INFO.MAPINFO_LCMAG.NDIM              = N4DIM_KCOR;
  CALIB_INFO.MAPINFO_LCMAG.NBIN[0]           = NBIN_T ;
  CALIB_INFO.MAPINFO_LCMAG.NBIN[1]           = NBIN_z ;
  CALIB_INFO.MAPINFO_LCMAG.NBIN[2]           = NBIN_AV ;
  CALIB_INFO.MAPINFO_LCMAG.NBIN[3]           = NFILTDEF_REST ;
  get_MAPINFO_KCOR("NBINTOT", &CALIB_INFO.MAPINFO_LCMAG);
  
  sprintf(CALIB_INFO.MAPINFO_MWXT.NAME, "MWXT");
  CALIB_INFO.MAPINFO_MWXT.IDMAP             = IDMAP_KCOR_MWXT ;
  CALIB_INFO.MAPINFO_MWXT.NDIM              = N4DIM_KCOR;
  CALIB_INFO.MAPINFO_MWXT.NBIN[0]           = NBIN_T ;
  CALIB_INFO.MAPINFO_MWXT.NBIN[1]           = NBIN_z ;
  CALIB_INFO.MAPINFO_MWXT.NBIN[2]           = NBIN_AV ;
  CALIB_INFO.MAPINFO_MWXT.NBIN[3]           = NFILTDEF_OBS ;
  get_MAPINFO_KCOR("NBINTOT", &CALIB_INFO.MAPINFO_MWXT);
  
  // clear map IDs since RDKCOR can be called multiple times.
  clear_1DINDEX(IDMAP_KCOR_TABLE);
  clear_1DINDEX(IDMAP_KCOR_AVWARP);
  clear_1DINDEX(IDMAP_KCOR_LCMAG);
  clear_1DINDEX(IDMAP_KCOR_MWXT);

  // init multi-dimensional index maps
  init_1DINDEX(IDMAP_KCOR_TABLE,  NKDIM_KCOR, CALIB_INFO.MAPINFO_KCOR.NBIN);
  init_1DINDEX(IDMAP_KCOR_AVWARP, N4DIM_KCOR, CALIB_INFO.MAPINFO_AVWARP.NBIN);
  init_1DINDEX(IDMAP_KCOR_LCMAG,  N4DIM_KCOR, CALIB_INFO.MAPINFO_LCMAG.NBIN);
  init_1DINDEX(IDMAP_KCOR_MWXT,   N4DIM_KCOR, CALIB_INFO.MAPINFO_MWXT.NBIN);

  return ;

} // end init_kcor_indices

// ===============================================
void get_MAPINFO_KCOR(char *what, KCOR_MAPINFO_DEF *MAPINFO) {
  int i, NDIM, NBIN ;
  char string_NBIN[40];
  // ------------- BEGIN ------------
  if ( strcmp(what,"NBINTOT") == 0 ) {
    NDIM = MAPINFO->NDIM;
    MAPINFO->NBINTOT = 1;
    string_NBIN[0] = 0;
    for(i=0; i < NDIM; i++ ) { 
      NBIN = MAPINFO->NBIN[i] ;
      MAPINFO->NBINTOT *= NBIN ; 
      if ( i == 0 ) 
	{ sprintf(string_NBIN,"%d", NBIN); }
      else
	{ sprintf(string_NBIN,"%s x %d", string_NBIN, NBIN); }
    }
    
    printf("\t\t NBINMAP(%-6s) = %s = %d \n",
	   MAPINFO->NAME, string_NBIN, MAPINFO->NBINTOT );
    fflush(stdout);
  }

  return ;

} // end get_MAPINFO_KCOR


// =============================
void read_calib_mags(void) {


  // Read LCMAG table for each rest-filter,
  // and read MWXT-slope for each obs-filter.
  //The columns for this table are
  //	 1-3: T,Z,AV
  //	 4  : 3+NFILTDEF_RDKCOR      :  MAGOBS(T,Z,AV)
  //	 Next NFILTDEF_RDKCOR bins  :  MWXTSLP(T,Z,AV)
  //
  // Store info only for filters that are used in a K-cor.

  fitsfile *FP         = CALIB_INFO.FP ;
  int  NBIN_T          = CALIB_INFO.BININFO_T.NBIN ;
  int  NBIN_z          = CALIB_INFO.BININFO_z.NBIN ;
  int  NBIN_AV         = CALIB_INFO.BININFO_AV.NBIN ;
  int  NFILTDEF_KCOR   = CALIB_INFO.NFILTDEF;
  int  NFILTDEF_REST   = CALIB_INFO.FILTERCAL_REST.NFILTDEF ;
  int  NFILTDEF_OBS    = CALIB_INFO.FILTERCAL_OBS.NFILTDEF ;

  int  NBINTOT_LCMAG   = NBIN_T * NBIN_z * NBIN_AV * NFILTDEF_REST;
  int  MEMF_LCMAG      = NBINTOT_LCMAG * sizeof(float);
  CALIB_INFO.LCMAG_TABLE1D_F = (float*)malloc(MEMF_LCMAG);

  int  NBINTOT_MWXT    = NBIN_T * NBIN_z * NBIN_AV * NFILTDEF_OBS;
  int  MEMF_MWXT       = NBINTOT_MWXT * sizeof(float);
  CALIB_INFO.MWXT_TABLE1D_F = (float*)malloc(MEMF_MWXT);

  int istat=0, hdutype, anynul, ifilt, ifiltr, ifilto, IFILTDEF ;
  int MASK, ISREST, ISOBS, ICOL_LCMAG, ICOL_MWXT;
  int IBLCMAG[N4DIM_KCOR], IBMWXT[N4DIM_KCOR];
  int IBIN_FIRST, IBIN_LAST, ibin, LBX;
  long long FIRSTROW=1, FIRSTELEM=1, NROW;
  char *CFILT ;
  char fnam[] = "read_calib_mags" ;

  // --------- BEGIN ----------
  printf("   %s \n", fnam); fflush(stdout);

  fits_movrel_hdu(FP, 1, &hdutype, &istat);
  snfitsio_errorCheck("Cannot move to MAG table", istat);

  if ( CALIB_INFO.NKCOR_STORE == 0 ) { return; }

  NROW = NBIN_T * NBIN_z * NBIN_AV;
  for(ibin=0; ibin < N4DIM_KCOR; ibin++ ) 
    { IBLCMAG[ibin] = IBMWXT[ibin] = 0 ;  }

  for (ifilt=0; ifilt < NFILTDEF_KCOR; ifilt++ ) {

    MASK   = CALIB_INFO.MASK_FRAME_FILTER[ifilt] ;
    ISREST = ( MASK & MASK_FRAME_REST);
    ISOBS  = ( MASK & MASK_FRAME_OBS);

    /*
    printf("\t xxx %s: ifilt=%d MASK=%d for '%s' \n",
	   fnam, ifilt, MASK, CALIB_INFO.FILTER_NAME[ifilt] );
    */

    if ( ! ( ISREST || ISOBS ) )  { continue ; }

    IFILTDEF = CALIB_INFO.IFILTDEF[ifilt] ;
    CFILT    = CALIB_INFO.FILTER_NAME[ifilt];
    ICOL_LCMAG   = 4 + ifilt;  // skip T,z,AV columns
    ICOL_MWXT    = ICOL_LCMAG + NFILTDEF_KCOR ;

    if ( ISREST ) {
      ifiltr = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[IFILTDEF];
      IBLCMAG[KDIM_IFILTr] = ifiltr ;
      IBIN_FIRST = get_1DINDEX(IDMAP_KCOR_LCMAG, N4DIM_KCOR, IBLCMAG) ;
      IBIN_LAST  = IBIN_FIRST + NROW - 1;

      //      IBIN_FIRST = -66;
      if ( IBIN_FIRST < 0  || IBIN_LAST >= NBINTOT_LCMAG ) {
	sprintf(c1err,"Invalid IBIN_FIRST,LAST(LCMAG) = %d,%d", 
		IBIN_FIRST, IBIN_LAST );
	sprintf(c2err,"for REST-filter = %s (ifiltr=%d, IFILTDEF=%d)",
		CFILT, ifiltr, IFILTDEF);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      fits_read_col_flt(FP, ICOL_LCMAG, FIRSTROW, FIRSTELEM, NROW,
			NULL_1E, &CALIB_INFO.LCMAG_TABLE1D_F[IBIN_FIRST], 
			&anynul, &istat );
      sprintf(c1err,"read LCMAG(%s)", CFILT);
      snfitsio_errorCheck(c1err, istat);

      // apply user mag-shifts
      for(ibin=IBIN_FIRST; ibin<=IBIN_LAST; ibin++ ) {
	CALIB_INFO.LCMAG_TABLE1D_F[ibin] += 
	  ( CALIB_INFO.MAGREST_SHIFT_PRIMARY[IFILTDEF] - 19.6);	  
      }
    } // end ISREST
    
    // now get MWXT cor for obs-frame filters.
    LBX = ISBXFILT_KCOR(CFILT) ;

    if ( ISOBS  && !LBX ) {
      ifilto = CALIB_INFO.FILTERCAL_OBS.IFILTDEF_INV[IFILTDEF];
      IBMWXT[KDIM_IFILTr] = ifilto ; // note index here is 3, not 4
      IBIN_FIRST = get_1DINDEX(IDMAP_KCOR_MWXT, N4DIM_KCOR, IBMWXT) ;
      IBIN_LAST  = IBIN_FIRST + NROW - 1;
      
      //      IBIN_FIRST = -66;
      if ( IBIN_FIRST < 0 || IBIN_LAST > NBINTOT_MWXT ) {
	sprintf(c1err,"Invalid IBIN_FIRST,LAST(MWXT) = %d,%d", 
		IBIN_FIRST, IBIN_LAST );
	sprintf(c2err,"for OBS-filter = %s (ifilto=%d, IFILTDEF=%d)",
		CFILT, ifilto, IFILTDEF);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      fits_read_col_flt(FP, ICOL_MWXT, FIRSTROW, FIRSTELEM, NROW,
			NULL_1E, &CALIB_INFO.MWXT_TABLE1D_F[IBIN_FIRST],
			&anynul, &istat );
      sprintf(c1err,"read MWXT-slope(%s)", CFILT);
      snfitsio_errorCheck(c1err, istat);

    } // end ISOBS    

  } // end ifilt loop


  /*xxxxxx dump entire table to compare with original fortran
  for(ibin=0; ibin < NBINTOT_LCMAG; ibin++ ) {
    printf(" CHECK LCMAG_TABLE1D[%6d] = %8.4f\n", 
	   ibin+1, CALIB_INFO.LCMAG_TABLE1D_F[ibin]); fflush(stdout);
  }

  for(ibin=0; ibin < NBINTOT_MWXT; ibin++ ) {
    printf(" CHECK MWXT_TABLE1D[%6d] = %8.4f\n", 
	   ibin+1, CALIB_INFO.MWXT_TABLE1D_F[ibin]); fflush(stdout);
  }
  xxxxxxxx */

  return ;

} // end read_calib_mags

// =============================
void read_calib_filters(void) {

  fitsfile *FP          = CALIB_INFO.FP ;
  int    NFILTDEF_KCOR  = CALIB_INFO.NFILTDEF;
  int    NBL            = CALIB_INFO.BININFO_LAM.NBIN;
  int    MEMF           = NBL * sizeof(float);
  float *ARRAY_LAM      = (float*)malloc(MEMF);
  float *ARRAY_TRANS    = (float*)malloc(MEMF);

  //  int  NFILTDEF_REST   = CALIB_INFO.FILTERCAL_REST.NFILTDEF ;
  //  int  NFILTDEF_OBS    = CALIB_INFO.FILTERCAL_OBS.NFILTDEF ;

  int istat=0, hdutype, anynul, ICOL, NMATCH_OBS, ifilt ;
  int MASK, IFILTDEF, IFILT_REST, IFILT_OBS ;
  long long FIRSTROW=1, FIRSTELEM=1 ;
  char FILTERLIST_READ[MXFILTINDX], *FILTER_NAME, FILTER_BAND[2];
  char *SURVEY_NAME ;
  char FRAME_REST[] = "REST" ;
  char FRAME_OBS[]  = "OBS" ;
  char fnam[] = "read_calib_filters" ;

  // --------- BEGIN ----------

  printf("   %s \n", fnam); fflush(stdout);

  fits_movrel_hdu(FP, 1, &hdutype, &istat);
  snfitsio_errorCheck("Cannot move to FILTERS table", istat);

  // read array of wavelength bins 
  ICOL=1 ;
  fits_read_col_flt(FP, ICOL, FIRSTROW, FIRSTELEM, NBL,
		    NULL_1E, ARRAY_LAM,	&anynul, &istat );
  sprintf(c1err,"read LAM array" );
  snfitsio_errorCheck(c1err,istat);

  // - - - - -

  NMATCH_OBS = 0 ;
  FILTERLIST_READ[0] = 0 ;

  for(ifilt=0; ifilt < NFILTDEF_KCOR; ifilt++ ) {

    ICOL        = 2 + ifilt;
    IFILTDEF    = CALIB_INFO.IFILTDEF[ifilt] ;
    FILTER_NAME = CALIB_INFO.FILTER_NAME[ifilt] ;
    SURVEY_NAME = CALIB_INFO.SURVEY_NAME[ifilt] ;
    sprintf(FILTER_BAND, "%c", FILTERSTRING[IFILTDEF] );
    strcat(FILTERLIST_READ,FILTER_BAND);

    fits_read_col_flt(FP, ICOL, FIRSTROW, FIRSTELEM, NBL,
		      NULL_1E, ARRAY_TRANS, &anynul, &istat );
    sprintf(c1err,"read %s filter trans", FILTER_NAME );
    snfitsio_errorCheck(c1err,istat);

    // match filter name to get absolute filter indices IFILT_REST & IFILT_OBS
    filter_match_kcor(FILTER_NAME, &IFILT_REST, &IFILT_OBS);

    /*
    printf(" xxx C: %s -> IFILT[REST,OBS] = %d, %d \n",
	   FILTER_NAME, IFILT_REST, IFILT_OBS ); fflush(stdout);
    */

    if ( IFILT_REST > 0 ) {
      CALIB_INFO.MASK_FRAME_FILTER[ifilt] = MASK_FRAME_REST;
      check_duplicate_filter(FRAME_REST, IFILT_REST, FILTER_NAME );
      loadFilterTrans_kcor(IFILT_REST, NBL, ARRAY_LAM, ARRAY_TRANS,
			   &CALIB_INFO.FILTERCAL_REST );	     
    } // end IFILT_REST

    if ( IFILT_OBS > 0 ) {
      NMATCH_OBS++ ;
      CALIB_INFO.MASK_FRAME_FILTER[ifilt] = MASK_FRAME_OBS;
      addFilter_kcor(IFILT_OBS, FILTER_NAME, &CALIB_INFO.FILTERCAL_OBS) ;
      check_duplicate_filter(FRAME_OBS, IFILT_OBS, FILTER_NAME );
      loadFilterTrans_kcor(IFILT_OBS, NBL, ARRAY_LAM, ARRAY_TRANS,
			   &CALIB_INFO.FILTERCAL_OBS );	

      // SHIFT_FILTTRANS function obsolete since lam shifts are in fit code
      // ?? FILTOBS_ZPOFF_SNPHOT(ifilt_obs) = ZPOFF_SNPHOT_RDKCOR(ifilt)
    }

  } // end ifilt loop


  if ( NMATCH_OBS == 0 ) {
    print_preAbort_banner(fnam);
    printf("\t Obs filters in kcor file: '%s' \n", 
	   CALIB_INFO.FILTERCAL_OBS.FILTERSTRING );
    printf("\t SURVEY_FILTERS: '%s' \n", CALIB_INFO.FILTERS_SURVEY);

    sprintf(c1err, "Observer filters do not match any SURVEY_FILTERS.");
    sprintf(c2err, "see PRE-ABORT dump above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // abort if any duplicate filters were found
  check_duplicate_filter(0, -1, 0 );

  return ;
} // end read_calib_filters

// =============================
void filter_match_kcor(char *NAME, int *IFILT_REST, int *IFILT_OBS) {

  // For input filter NAME, returns IFILT_REST and IFILT_OBS
  // Usually only one of the IFILT_[REST,OBS] indices is valid,
  // but note that both can be used. Also note that NAME is the
  // full filter-name, not just the last character.
  //

  int  NFILTDEF_REST   = CALIB_INFO.FILTERCAL_REST.NFILTDEF ;
  int  NFILTDEF_OBS    = CALIB_INFO.FILTERCAL_OBS.NFILTDEF ;
  int  IFILTDEF        = INTFILTER(NAME);
  bool IS_SURVEY_FILTER = CALIB_INFO.IS_SURVEY_FILTER[IFILTDEF];
  int  ifilt ;
  char *NAME_REST, *NAME_OBS ;
  char fnam[] = "filter_match_kcor" ;

  // ---------- BEGIN -----------

  *IFILT_REST = *IFILT_OBS = -9;

  // check if rest-frame filter
  for(ifilt=0; ifilt < NFILTDEF_REST; ifilt++ ) {
    NAME_REST = CALIB_INFO.FILTERCAL_REST.FILTER_NAME[ifilt];
    if ( strcmp(NAME,NAME_REST) == 0 )  { *IFILT_REST = IFILTDEF ; }
  }


  //  printf(" xxx %s: %-10.10s  IFILTDEF=%2d  IS_SURVEY_FILTER=%d \n",
  //	 fnam, NAME, IFILTDEF, IS_SURVEY_FILTER );

  // continue only if this obs-frame filter is a survey filter
  if ( !IS_SURVEY_FILTER ) { return; }

  // if not a rest-frame filter, then it MUST  be an obs-frame filter.  
  // (e..g, explicit or from spectrograph)
  if ( *IFILT_REST < 0 ) {  *IFILT_OBS = IFILTDEF; return ;    }


  // check if obs-frame filter
  for(ifilt=0; ifilt < NFILTDEF_OBS; ifilt++ ) {
    NAME_OBS  = CALIB_INFO.FILTERCAL_OBS.FILTER_NAME[ifilt] ;
    if ( strcmp(NAME_OBS,NAME) == 0 ) { *IFILT_OBS = IFILTDEF ; }
  }

  return ;

} // end filter_match_kcor


// ==============================================
void check_duplicate_filter(char *FRAME, int IFILTDEF, char *FILTER_NAME ) {

  // IFILTDEF > 0 -> give warning if dupliate
  // IFILTDEF < 0 -> abort with summary of duplicates

  int  NDEFINE=0 ;
  char fnam[] = "check_duplicate_filter" ;

  // --------- BEGIN ------------


  if ( IFILTDEF > 0 ) {
    if ( strcmp(FRAME,"REST") == 0 )  { 
      CALIB_INFO.FILTERCAL_REST.NDEFINE[IFILTDEF]++ ;
      NDEFINE = CALIB_INFO.FILTERCAL_REST.NDEFINE[IFILTDEF];  
      if ( NDEFINE == 2 ) { CALIB_INFO.FILTERCAL_REST.NFILT_DUPLICATE++; }

    }
    else if ( strcmp(FRAME,"OBS") == 0 ) { 
      CALIB_INFO.FILTERCAL_OBS.NDEFINE[IFILTDEF]++ ;
      NDEFINE = CALIB_INFO.FILTERCAL_OBS.NDEFINE[IFILTDEF]; 
      if ( NDEFINE == 2 ) { CALIB_INFO.FILTERCAL_OBS.NFILT_DUPLICATE++; }
    }
  
    // give warning on duplicate, but do not abort (yet).
    if ( NDEFINE > 1 ) {
      sprintf(c1err,"NDEFINE=%d -> duplicate %s filter '%s' (%d) ",
	      NDEFINE, FRAME, FILTER_NAME, IFILTDEF);
      sprintf(c2err,"Check kcor-input");
      errmsg(SEV_WARN, 0, fnam, c1err, c2err); 
    }
  }   
  else if ( IFILTDEF < 0 ) {
    // check final summary  
    int NDUP_REST = CALIB_INFO.FILTERCAL_REST.NFILT_DUPLICATE ;
    int NDUP_OBS  = CALIB_INFO.FILTERCAL_OBS.NFILT_DUPLICATE ;
    if ( NDUP_REST > 0 || NDUP_OBS > 0 ) {
      sprintf(c1err,"%d/%d duplicate REST/OBS filters", 
	      NDUP_REST, NDUP_OBS);
      sprintf(c2err,"Check duplicate warnings above.");
      errmsg(SEV_WARN, 0, fnam, c1err, c2err); 
    }
  }


  return ;

} // end check_duplicate_filter


// ======================================================
void loadFilterTrans_kcor(int IFILTDEF, int NBL, 
			  float *ARRAY_LAM, float *ARRAY_TRANS,
			  FILTERCAL_DEF *MAP) {

  //
  // store filter trans info in MAP structure.
  // Note that storage is float, but calculations (rms, mean, ...)
  // are done with double precision.
  //
  // Inputs:
  //   IFILTDEF     : absolute filter index
  //   NBL          : number of lambda bins
  //   ARRAY_LAM    : lambda array to store
  //   ARRAY_TRANS  : transmmission array to store
  //

  int OPT_FRAME = MAP->OPT_FRAME ;
  int MEMF  = NBL * sizeof(float);
  int ilam, ifilt, ilam_min, ilam_max ;
  double LAM, TRANS, MEAN, SQRMS, LAMMIN=1.1E5, LAMMAX=0.0 ;
  double TMAX=0.0, SUM0=0.0, SUM1=0.0, SUM2=0.0 ;
  char fnam[] = "loadFilterTrans_kcor" ;

  // ---------------- BEGIN ---------------

  if ( OPT_FRAME == OPT_FRAME_REST ) 
    { ifilt = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[IFILTDEF]; }
  else
    { ifilt = CALIB_INFO.FILTERCAL_OBS.IFILTDEF_INV[IFILTDEF]; }

  MAP->LAM[ifilt]      = (float*)malloc(MEMF);
  MAP->TRANS[ifilt]    = (float*)malloc(MEMF);
  
  // determine min and max ilam that contains Trans>0 
  ilam_min = 9999999;  ilam_max = -9;
  for(ilam=0; ilam < NBL; ilam++ ) {
    TRANS = (double)ARRAY_TRANS[ilam];
    if ( TRANS > 1.0E-6  ) {
      if ( ilam_min > 99999 ) { ilam_min = ilam; }
      ilam_max = ilam;
    }
  }
  // add extra edge bin with zero transmission
  if ( ilam_min > 0     ) { ilam_min-- ; }
  if ( ilam_max < NBL-1 ) { ilam_max++ ; }

  // - - - - - 
  int NBL_STORE = 0 ;
  for(ilam=ilam_min; ilam <= ilam_max; ilam++ ) {
    LAM   = (double)ARRAY_LAM[ilam];
    TRANS = (double)ARRAY_TRANS[ilam];

    SUM0 += TRANS;
    SUM1 += (TRANS * LAM);
    SUM2 += (TRANS * LAM * LAM);

    MAP->LAM[ifilt][NBL_STORE]   = (float)LAM ;
    MAP->TRANS[ifilt][NBL_STORE] = (float)TRANS ;

    NBL_STORE++ ;
  } // end ilam

  MEAN  = SUM1/SUM0;
  SQRMS = SUM2/SUM0 - MEAN*MEAN;

  // load extra info about transmission function
  MAP->NBIN_LAM[ifilt]  = NBL_STORE  ;
  MAP->TRANS_MAX[ifilt] = TMAX;    // max trans
  MAP->LAMMEAN[ifilt]   = MEAN ;   // mean wavelength
  MAP->LAMRMS[ifilt]    = sqrt(SQRMS) ; // RMS wavelength
  MAP->LAMRANGE[ifilt][0] = MAP->LAM[ifilt][0] ;
  MAP->LAMRANGE[ifilt][1] = MAP->LAM[ifilt][NBL_STORE-1] ;

  MAP->LAMRANGE_KCOR[ifilt][0] = -9.0 ; // compute later for rest-frame
  MAP->LAMRANGE_KCOR[ifilt][1] = -9.0 ;
  
  /* xxx
  printf(" xxx C: IFILTDEF=%2d  PRIMARY(MAG,ZPOFF) = %.3f, %.3f  (NBL=%d)\n",
	 IFILTDEF, MAP->PRIMARY_MAG[ifilt], 
	 MAP->PRIMARY_ZPOFF_SYN[ifilt], NBL );
  xxx */

  fflush(stdout);
  
  return ;

} // end loadFilterTrans_kcor

// =====================================
void read_calib_primarysed(void) {

  fitsfile *FP          = CALIB_INFO.FP ;
  int  NFILTDEF_OBS     = CALIB_INFO.FILTERCAL_OBS.NFILTDEF ;
  int  NBL              = CALIB_INFO.BININFO_LAM.NBIN; // from SED
  int istat=0, hdutype, anynul, ifilt;
  int KINDX=-9, KINDX_FIRST=-9, KINDX_2ND=-9, NERR_PRIM=0 ;
  char *NAME ;
  char fnam[] = "read_calib_primarysed" ;

  // --------- BEGIN ----------

  printf("   %s \n", fnam); fflush(stdout);

  fits_movrel_hdu(FP, 1, &hdutype, &istat);
  snfitsio_errorCheck("Cannot move to PRIMARYSED table", istat);

  for(ifilt=0; ifilt < NFILTDEF_OBS; ifilt++ ) {
    KINDX = CALIB_INFO.FILTERCAL_OBS.PRIMARY_KINDX[ifilt];
    if ( KINDX_FIRST < 0 ) { KINDX_FIRST = KINDX; }
    if ( KINDX != KINDX_FIRST ) { NERR_PRIM++ ; KINDX_2ND=KINDX; } 
  }

  if ( KINDX < 0 ) {
    sprintf(c1err,"Could not find primary reference");
    sprintf(c2err,"Something is messed up.");
    errmsg(SEV_WARN, 0, fnam, c1err, c2err); 
  }

  if ( NERR_PRIM > 0 ) {
    print_preAbort_banner(fnam);
    printf("   Found Primary %s\n", CALIB_INFO.PRIMARY_NAME[KINDX_FIRST] );
    printf("   Found Primary %s\n", CALIB_INFO.PRIMARY_NAME[KINDX_2ND] );
    sprintf(c1err, "More than one PRIMARY ref not allowed");
    sprintf(c2err, "Check kcor file and above list of primary refs.");
    errmsg(SEV_WARN, 0, fnam, c1err, c2err); 
  }

  
  NAME = CALIB_INFO.PRIMARY_NAME[KINDX];
  // xxx mark  printf("\t\t Primary Reference: %s\n", NAME);

  int       MEMF = NBL * sizeof(float); 
  int       ICOL_LAM, ICOL_FLUX ;
  long long FIRSTROW=1, FIRSTELEM=1;
  float  *ptr_f;

  // read lambda array
  CALIB_INFO.FILTERCAL_OBS.NBIN_LAM_PRIMARY = NBL;
  CALIB_INFO.FILTERCAL_OBS.PRIMARY_LAM  = (float*) malloc(MEMF);
  CALIB_INFO.FILTERCAL_OBS.PRIMARY_FLUX = (float*) malloc(MEMF);

  // read wavelength array (should be same as array for SN SED)
  ICOL_LAM=1;    ptr_f =  CALIB_INFO.FILTERCAL_OBS.PRIMARY_LAM ;
  fits_read_col_flt(FP, ICOL_LAM, FIRSTROW, FIRSTELEM, NBL,
		    NULL_1E, ptr_f, &anynul, &istat );
  sprintf(c1err,"read lam array for primary = '%s'", NAME);
  snfitsio_errorCheck(c1err, istat);

  // read primary flux array
  ICOL_FLUX= ICOL_LAM + 1 + KINDX;  
  ptr_f =  CALIB_INFO.FILTERCAL_OBS.PRIMARY_FLUX ;
  fits_read_col_flt(FP, ICOL_FLUX, FIRSTROW, FIRSTELEM, NBL,
		    NULL_1E, ptr_f, &anynul, &istat );
  sprintf(c1err,"read flux array for primary = '%s'", NAME);
  snfitsio_errorCheck(c1err, istat);

  return ;

} // end read_calib_primarysed

// ====================================
void print_calib_summary(void)  {

  int  NFILTDEF_REST   = CALIB_INFO.FILTERCAL_REST.NFILTDEF ;
  int  NFILTDEF_OBS    = CALIB_INFO.FILTERCAL_OBS.NFILTDEF ;
  int  KINDX           = CALIB_INFO.FILTERCAL_OBS.PRIMARY_KINDX[0];

  int  ifilt, ifiltdef, k;
  double lamavg, lamrms, *lamrange, prim_mag, prim_zpoff;
  char *NAME;
  char dashLine[]  = 
    "----------------------------------------------------------------------" ;
  
  char fnam[] = "print_calib_summary" ;

  // ------------- BEGIN ------------

  printf("\n   %s \n", fnam); fflush(stdout);
  // xxx  print_banner(fnam);
  
  printf("\t  Primary Spectrum: %s \n", CALIB_INFO.PRIMARY_NAME[KINDX]);


  printf("\n  FILTER SUMMARY (RDKCOR) : \n");
  printf("   internal                   LAM  LAM   rest-lam      primary\n");
  printf("  index name                  AVG  RMS   kcor range    mag   ZP\n");
  printf(" %s \n", dashLine);

  // print filters in the same order as in the KCOR file,
  // rather than in the order they appear in the FILTDEF_STRING .

  /* xxx ?  
  char fmt[] = 
    "  %3d %-20.20s %6d %4d "
    "%5d-%5d "
    "%7.3f %7.3f\n" ;
  xxx */

  
  for(ifilt=0; ifilt < NFILTDEF_OBS; ifilt++ ) {
    ifiltdef = CALIB_INFO.FILTERCAL_OBS.IFILTDEF[ifilt]; 
    NAME     = CALIB_INFO.FILTERCAL_OBS.FILTER_NAME[ifilt] ;
    lamavg   = CALIB_INFO.FILTERCAL_OBS.LAMMEAN[ifilt]; 
    lamrms   = CALIB_INFO.FILTERCAL_OBS.LAMRMS[ifilt];     
    prim_mag   = CALIB_INFO.FILTERCAL_OBS.PRIMARY_MAG[ifilt];
    prim_zpoff = CALIB_INFO.FILTERCAL_OBS.PRIMARY_ZPOFF_SYN[ifilt];

    printf("  %3d %-20.20s %6d %4d %5d-%5d %7.3f %7.3f\n",
	   ifiltdef, NAME, (int)(lamavg+0.5), (int)(lamrms+0.5),
	   0,0,  prim_mag, prim_zpoff);
    fflush(stdout);
  }


  // - - - - 
  if ( NFILTDEF_REST > 0 )     { printf("\n"); }


  for(ifilt=0; ifilt < NFILTDEF_REST; ifilt++ ) {
    ifiltdef = CALIB_INFO.FILTERCAL_REST.IFILTDEF[ifilt]; 
    NAME     = CALIB_INFO.FILTERCAL_REST.FILTER_NAME[ifilt] ;
    lamavg   = CALIB_INFO.FILTERCAL_REST.LAMMEAN[ifilt]; 
    lamrms   = CALIB_INFO.FILTERCAL_REST.LAMRMS[ifilt];     

    prim_mag   = CALIB_INFO.FILTERCAL_REST.PRIMARY_MAG[ifilt];
    prim_zpoff = CALIB_INFO.FILTERCAL_REST.PRIMARY_ZPOFF_SYN[ifilt];

    set_lamrest_range_KCOR(ifiltdef);
    lamrange = CALIB_INFO.FILTERCAL_REST.LAMRANGE_KCOR[ifilt];     

    printf("  %3d %-20.20s %6d %4d %5d-%5d %7.3f %7.3f\n",
	   ifiltdef, NAME, (int)(lamavg+0.5), (int)(lamrms+0.5),
	   (int)(lamrange[0]+0.5), (int)(lamrange[1]+0.5),  
	   prim_mag, prim_zpoff);
    fflush(stdout);
  }


  printf("\t NFILTDEF[SURVEY,REST] = %d %d \n",
	 NFILTDEF_OBS, NFILTDEF_REST);


  // - - - - - KCORs - - - - -
  int NKCOR = CALIB_INFO.NKCOR_STORE ;
  int ifiltdef_o, ifiltdef_r, ifilt_o, ifilt_r;
  char *name_o, *name_r;
  if ( NKCOR == 0 ) { return; }
  printf("\n KCOR SUMMARY (%d tables)\n", NKCOR);
  for(k=0; k < NKCOR; k++ ) {
    ifiltdef_o = CALIB_INFO.IFILTMAP_KCOR[OPT_FRAME_OBS][k] ;
    ifiltdef_r = CALIB_INFO.IFILTMAP_KCOR[OPT_FRAME_REST][k] ;
    ifilt_o    = CALIB_INFO.FILTERCAL_OBS.IFILTDEF_INV[ifiltdef_o];
    ifilt_r    = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifiltdef_r];
    name_o     = CALIB_INFO.FILTERCAL_OBS.FILTER_NAME[ifilt_o];
    name_r     = CALIB_INFO.FILTERCAL_REST.FILTER_NAME[ifilt_r];
    printf("     %3d Found %s table    (%d -> %d : %s -> %s)\n",
	   k, CALIB_INFO.STRING_KCORSYM[k], ifiltdef_r, ifiltdef_o,
	   name_r, name_o);
  }

  fflush(stdout);
  return;

} // end print_calib_summary


// =====================================
void set_lamrest_range_KCOR(int ifiltdef) {

  // Nov 2022
  // Translate fortran subroutine SET_LAMREST_RANGE here,
  // but leave out user-override feature based on SNLCINP
  // input OVERRIDE_RESTLAM_BOUNDARY ... may add this later.
  //
  // Load LAMRANGE_KCOR[ifilt] for each rest-frame filter;
  // used to pick rest-frame band for KCOR.

  int  ifilt            = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifiltdef];
  int  NFILTDEF_REST    = CALIB_INFO.FILTERCAL_REST.NFILTDEF ;
  char *FILTERSTRING    = CALIB_INFO.FILTERCAL_REST.FILTERSTRING ;
  double *LAMRANGE_KCOR = CALIB_INFO.FILTERCAL_REST.LAMRANGE_KCOR[ifilt];
  double *LAMRANGE      = CALIB_INFO.FILTERCAL_REST.LAMRANGE[ifilt]; 
  double lamavg         = CALIB_INFO.FILTERCAL_REST.LAMMEAN[ifilt];

  int  ifilt_tmp, ifiltdef_tmp, ifiltdef_near[2], j ;
  double lamavg_tmp, lamdif, lamdif_near[2];
  char fnam[] = "set_lamrest_range_KCOR" ;
  int LDMP = ( ifiltdef == -12 ) ; // 11=U, 12=B

  // ------------ BEGIN -------------

  // if rest-frame bands are Bessell UBVRI, use hard-wired boundaries
  // for mlcs2k2.
  if ( strncmp(FILTERSTRING, "UBVRI", 5) == 0 ) {
    set_lamrest_range_UBVRI(ifiltdef);
    return;
  }

  // find nearest rest-frame filter
  lamdif_near[0] = 99999. ;
  lamdif_near[1] = 99999. ;
  ifiltdef_near[0] = -9;
  ifiltdef_near[1] = -9;

  for ( ifilt_tmp=0; ifilt_tmp < NFILTDEF_REST; ifilt_tmp++ ) {
    ifiltdef_tmp = CALIB_INFO.FILTERCAL_REST.IFILTDEF[ifilt_tmp];
    if ( ifiltdef_tmp == ifiltdef ) { continue; }
    if ( ifiltdef_tmp == IFILTDEF_BESS_BX ) { continue; }

    lamavg_tmp = CALIB_INFO.FILTERCAL_REST.LAMMEAN[ifilt_tmp];
    lamdif     = fabs(lamavg_tmp - lamavg);

    if ( lamavg_tmp < lamavg && lamdif < lamdif_near[0] ) {
      lamdif_near[0]   = lamdif;
      ifiltdef_near[0] = ifiltdef_tmp ;
    }

    if ( lamavg_tmp > lamavg && lamdif < lamdif_near[1] ) {
      lamdif_near[1]   = lamdif;
      ifiltdef_near[1] = ifiltdef_tmp;
    }

  } // end ifilt_tmp loop

  // - - - -
  // set extreme bound if there is no nbr filter
  if ( ifiltdef_near[0] < 0 ) {
    LAMRANGE_KCOR[0] = LAMRANGE[0] / 2.0 ;
  }
  if ( ifiltdef_near[1] < 0 ) {
    LAMRANGE_KCOR[1] = LAMRANGE[1] + 1000.0 ;
  }

  // - - - - - - -
  // if voisin filter has smaller LAMAVG,
  // then look for where transmissions are the same
  int iminmax, ifiltdef_nbr, ifilt_nbr, nbin_ovp, ilam, ilam_nbr;
  double lamavg_nbr, lamovp, lamovp_avg;
  double transmax, transmax_nbr, transdif_min=99999.9 ;
  double trans, trans_nbr, transdif, trans_ovp, lam, lam_nbr ;
  int nbin, nbin_nbr;
  bool VALID;

  for(iminmax=0; iminmax < 2; iminmax++ ) {

    ifiltdef_nbr = ifiltdef_near[iminmax];
    if ( ifiltdef_nbr < 0 ) { continue; }

    //    printf(" xxx %s: iminmax=%d  filtdef=%d  ifiltdef_nbr=%d \n",
    //	   fnam, iminmax ,ifiltdef, ifiltdef_nbr);

    ifilt_nbr   = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifiltdef_nbr];
    if ( ifilt_nbr < 0 ) { continue; }

    lamavg_nbr  = CALIB_INFO.FILTERCAL_REST.LAMMEAN[ifilt_nbr];
    if ( lamavg_nbr == lamavg ) { continue; }

    lamovp       = LAMRANGE[iminmax];
    transmax     = CALIB_INFO.FILTERCAL_REST.TRANS_MAX[ifilt];
    transmax_nbr = CALIB_INFO.FILTERCAL_REST.TRANS_MAX[ifilt_nbr];

    nbin     = CALIB_INFO.FILTERCAL_REST.NBIN_LAM[ifilt];
    nbin_nbr = CALIB_INFO.FILTERCAL_REST.NBIN_LAM[ifilt_nbr];
    
    transdif_min = 9999.0 ;
    trans_ovp    = -9.0 ;
    nbin_ovp     = 0 ;

    for(ilam=0; ilam < nbin; ilam++ ) {
      trans = CALIB_INFO.FILTERCAL_REST.TRANS[ifilt][ilam];
      trans /= transmax; 
      if ( trans < 1.0E-10 ) { continue; }

      for(ilam_nbr=0; ilam_nbr < nbin_nbr; ilam_nbr++ ) {
	trans_nbr = CALIB_INFO.FILTERCAL_REST.TRANS[ifilt_nbr][ilam_nbr];
	trans_nbr /= transmax_nbr ;
	if ( trans_nbr < 1.0E-10 ) { continue; }

	lam     = CALIB_INFO.FILTERCAL_REST.LAM[ifilt][ilam];
	lam_nbr = CALIB_INFO.FILTERCAL_REST.LAM[ifilt_nbr][ilam_nbr];
	if ( fabs(lam-lam_nbr) > 1.0 ) { continue; }

	transdif = fabs(trans - trans_nbr);       
	if ( transdif < transdif_min ) {
	  transdif_min = transdif ;
	  trans_ovp    = trans_nbr;
	  lamovp       = lam;
	  nbin_ovp++ ;

	  if ( LDMP ) {
	    printf(" zzz transdif = %f - %f = %f  lam=%.1f/%1.f "
		   "(ilam=%d/%d)\n",
		   trans, trans_nbr, transdif, lam, lam_nbr, ilam, ilam_nbr);
	  }

	}

      } // end ilam_nbr
    } // end ilam 

    if ( LDMP ) {
      printf(" xxx ------------- iminmax=%d ------------- \n", iminmax);
      printf(" xxx %s: ifiltdef=%d  transdif_min=%.3f at trans=%.5f (%d)\n",
	     fnam, ifiltdef, transdif_min, trans_ovp, nbin_ovp);
      printf(" xxx %s: lamavg=%.1f  lamavg_nbr=%.1f  lamovp=%.1f\n",
	     fnam, lamavg, lamavg_nbr, lamovp );
    }

    // if LAMOVP is unphysical, take average of filter and nbr
    lamovp_avg = 0.5 * (lamavg + lamavg_nbr);

    VALID = lamovp >= lamavg_nbr  && lamovp <= lamavg ;
    if ( iminmax==0 && !VALID ) { lamovp = lamovp_avg;  }

    VALID = lamovp <= lamavg_nbr  && lamovp >= lamavg ;
    if ( iminmax == 1 && !VALID ) { lamovp = lamovp_avg; }

    // if the wavelength ranges do not overlap, take average
    if ( nbin_ovp == 0 ) { lamovp = lamovp_avg; }

    // finally, load wavelength edge in global
    LAMRANGE_KCOR[iminmax] = lamovp;

  } // end iminmax

  // - - - - -

  return;

} // end set_lamrest_range_KCOR

void set_lamrest_range_UBVRI(int ifiltdef) {

  int  ifilt = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifiltdef];
  double *LAMRANGE_KCOR = CALIB_INFO.FILTERCAL_REST.LAMRANGE_KCOR[ifilt];
  char fnam[] = "set_lamrest_range_UBVRI";

  // ----------- BEGIN ----------

  
  if ( ifiltdef == INTFILTER("U") ) {
    LAMRANGE_KCOR[0] = 1000.0 ;
    LAMRANGE_KCOR[1] = 3900.0 ;
  }
  else if ( ifiltdef == INTFILTER("B") ) {
    LAMRANGE_KCOR[0] = 3900.0 ;
    LAMRANGE_KCOR[1] = 4850.0 ;
  }
  else if ( ifiltdef == INTFILTER("X") ) {
    LAMRANGE_KCOR[0] = 3615.0 ;
    LAMRANGE_KCOR[1] = 5595.0 ;
  }
  else if ( ifiltdef == INTFILTER("V") ) {
    LAMRANGE_KCOR[0] = 4850.0 ;
    LAMRANGE_KCOR[1] = 5850.0 ;
  }
  else if ( ifiltdef == INTFILTER("R") ) {
    LAMRANGE_KCOR[0] = 5850.0 ;
    LAMRANGE_KCOR[1] = 7050.0 ;
  }
  else if ( ifiltdef == INTFILTER("I") ) {
    LAMRANGE_KCOR[0] = 7050.0 ;

    int  ifiltdef_Y = INTFILTER("Y");
    int  ifilt_Y    = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifiltdef_Y];
    if ( ifilt_Y > 0 ) {
      LAMRANGE_KCOR[1] = CALIB_INFO.FILTERCAL_REST.LAMRANGE[ifilt_Y][1];
    }
    else {
      LAMRANGE_KCOR[1] = CALIB_INFO.FILTERCAL_REST.LAMRANGE[ifilt][1];
    }
  }

  return;

} // end set_lamrest_range_UBVRI


// =====================================================
void get_calib_primary(char *primary_name, int *nblam,
                      double *lam, double *flux) {

  // Created Nov 2022
  // Return information about primary reference.
 
  int  ilam, kindx;
  char fnam[] = "get_calib_primary" ;
  FILTERCAL_DEF *FILTERCAL_OBS = &CALIB_INFO.FILTERCAL_OBS ;
  int  NBIN_LAM = FILTERCAL_OBS->NBIN_LAM_PRIMARY ;

  // ----------- BEGIN ----------


  // make sure to get primary for OBS filters (not rest-frame filters)
  kindx = FILTERCAL_OBS->PRIMARY_KINDX[0];
  sprintf(primary_name, "%s", CALIB_INFO.PRIMARY_NAME[kindx]);

  *nblam = NBIN_LAM ;

  for(ilam=0; ilam < NBIN_LAM; ilam++ ) {
    lam[ilam]  = (double)FILTERCAL_OBS->PRIMARY_LAM[ilam];
    flux[ilam] = (double)FILTERCAL_OBS->PRIMARY_FLUX[ilam];
  } // end ilam loop

  return;
} // end get_calib_primary

void get_calib_primary__(char *primary_name, int *NBLAM,
                        double *lam, double *flux) {
  get_calib_primary(primary_name, NBLAM, lam, flux);
} 



void get_KCOR_FILTERCAL(int OPT_FRAME, char *fnam, FILTERCAL_DEF *MAP) {

  // !!! xxx works locally, but returned MAP is corrupt xxx !!!
  //
  // for input OPT_FRAME, return pointer to FILTERCAL_DEF struct.
  // Note that fnam argument is used for abort message.
  if ( OPT_FRAME == OPT_FRAME_OBS ) { 
    MAP = &CALIB_INFO.FILTERCAL_OBS ;
  }
  else if ( OPT_FRAME == OPT_FRAME_REST ) {
    MAP = &CALIB_INFO.FILTERCAL_REST ; 
  }
  else {
    sprintf(c1err,"Invalid OPT_FRAME = %d", OPT_FRAME);
    sprintf(c2err,"Must be either OPT_FRAME_REST=%d or OPT_FRAME_OBS=%d",
	    OPT_FRAME_REST, OPT_FRAME_OBS);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

} // end get_KCOR_FILTERCAL


// =====================================================
void get_calib_filterTrans(int OPT_FRAME, int ifiltdef, char *surveyName,
                          char *filterName, double *magprim,
			  int *nblam, double *lam,
                          double *transSN, double *transREF) {

  // Created Nov 2022
  // Return information about filter "ifiltdef" for
  // OPT_FRAME = 0(REST) or 1(OBS)
  
  FILTERCAL_DEF *FILTERCAL;
  int ifilt, ilam ;
  char fnam[] = "get_calib_filterTrans" ;

  // ------------- BEGIN --------------
  
  // xxx returns corrupt FILTERCAL??? get_KCOR_FILTERCAL(OPT_FRAME, fnam, &FILTERCAL );

  if ( OPT_FRAME == OPT_FRAME_OBS ) { 
    FILTERCAL = &CALIB_INFO.FILTERCAL_OBS ;
  }
  else if ( OPT_FRAME == OPT_FRAME_REST ) {
    FILTERCAL = &CALIB_INFO.FILTERCAL_REST ; 
  }
  else {
    sprintf(c1err,"Invalid OPT_FRAME = %d", OPT_FRAME);
    sprintf(c2err,"Must be either OPT_FRAME_REST=%d or OPT_FRAME_OBS=%d",
	    OPT_FRAME_REST, OPT_FRAME_OBS);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  ifilt      = FILTERCAL->IFILTDEF_INV[ifiltdef];    

  sprintf(surveyName,"%s", FILTERCAL->SURVEY_NAME[ifilt] );
  sprintf(filterName,"%s", FILTERCAL->FILTER_NAME[ifilt] );

  *magprim   = FILTERCAL->PRIMARY_MAG[ifilt];

  /*
  printf(" xxx %s: survey=%s  filter=%s  nblam=%d magprim=%.3f\n",
	 fnam, surveyName, filterName, *nblam, *magprim); fflush(stdout);
  */

  *nblam     = FILTERCAL->NBIN_LAM[ifilt] ;
  for(ilam=0; ilam < *nblam; ilam++ ) {
    lam[ilam]      = FILTERCAL->LAM[ifilt][ilam];
    transSN[ilam]  = FILTERCAL->TRANS[ifilt][ilam];
    transREF[ilam] = FILTERCAL->TRANS[ifilt][ilam];
  }

  // TO-DO: check filter-update for SNLS that has filterTrans per SN
  //.xyz

  return;

} // end of get_calib_filterTrans

void get_calib_filtertrans__(int *OPT_FRAME, int *ifilt_obs, char *surveyName,
                            char *filterName,
                            double *magprim, int *nblam, double *lam,
                            double *transSN, double *transREF) {
  get_calib_filterTrans(*OPT_FRAME, *ifilt_obs, surveyName, filterName,
		       magprim, nblam, lam, transSN, transREF);
}

// ===========================================================
void get_calib_filtlam_stats(int OPT_FRAME, int ifiltdef,
			    double *lamavg, double *lamrms,
			    double *lammin, double *lammax) {

  int  ifilt;
  FILTERCAL_DEF *FILTERCAL;
  char fnam[] = "get_calib_filtlam_stats" ;

  // ----------- BEGIN ----------

  if ( OPT_FRAME == OPT_FRAME_OBS ) 
    { FILTERCAL = &CALIB_INFO.FILTERCAL_OBS; }
  else if ( OPT_FRAME == OPT_FRAME_REST ) 
    { FILTERCAL = &CALIB_INFO.FILTERCAL_REST; }

  ifilt      = FILTERCAL->IFILTDEF_INV[ifiltdef];    

  *lamavg = FILTERCAL->LAMMEAN[ifilt];
  *lamrms = FILTERCAL->LAMRMS[ifilt];
  *lammin = FILTERCAL->LAMRANGE[ifilt][0];
  *lammax = FILTERCAL->LAMRANGE[ifilt][1];

  return;
} // end get_calib_filtlam_stats

void get_calib_filtlam_stats__(int *OPT_FRAME, int *ifiltdef,
			      double *lamavg, double *lamrms,
			      double *lammin, double *lammax) {
  get_calib_filtlam_stats(*OPT_FRAME, *ifiltdef,
			 lamavg, lamrms, lammin, lammax);

}

// ==========================================================
void get_calib_filtindex_map(int OPT_FRAME, int *NFILTDEF, int *IFILTDEF_MAP,
			   int *IFILTDEF_INVMAP) {

  // return list of IFILTDEF filter indices vs. sparse index (ifilt)
  // and inverse map.

  int  ifilt, ifiltdef ;
  FILTERCAL_DEF *FILTERCAL;
  char fnam[] = "get_calib_filtindex_map" ;

  // ---------- BEGIN ------------

  if ( OPT_FRAME == OPT_FRAME_OBS ) 
    { FILTERCAL = &CALIB_INFO.FILTERCAL_OBS; }
  else if ( OPT_FRAME == OPT_FRAME_REST ) 
    { FILTERCAL = &CALIB_INFO.FILTERCAL_REST; }

  //  ifilt      = FILTERCAL->IFILTDEF_INV[ifiltdef];    

  *NFILTDEF = FILTERCAL->NFILTDEF ;
  for(ifilt=0; ifilt < *NFILTDEF; ifilt++ ) {
    ifiltdef  = FILTERCAL->IFILTDEF[ifilt];  
    IFILTDEF_MAP[ifilt]       = ifiltdef ;
    IFILTDEF_INVMAP[ifiltdef] = ifilt;
  }

  return;

} // end get_calib_filtindex_map

void get_calib_filtindex_map__(int *OPT_FRAME, int *NFILTDEF, int *IFILTDEF_MAP,
			      int *IFILTDEF_INVMAP) {
  get_calib_filtindex_map(*OPT_FRAME, NFILTDEF, IFILTDEF_MAP, IFILTDEF_INVMAP);
} 


double get_calib_zpoff_file(int ifiltdef) {
  // Return ZPOFF read from ZPOFF.DAT file in filter directory
  int    ifilt = CALIB_INFO.FILTERCAL_OBS.IFILTDEF_INV[ifiltdef];
  double zpoff = CALIB_INFO.PRIMARY_ZPOFF_FILE[ifilt] ;
  return zpoff;
} // end get_calib_zpoff_file

double get_calib_zpoff_file__(int *ifiltdef)
{ return get_calib_zpoff_file(*ifiltdef); }

// === END ===
