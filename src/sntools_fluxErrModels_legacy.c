// **************************************************
// **************************************************
// **************************************************
//
//  Created Feb 12 2018 (SNANA v10_58j)
//
//  Catenate sntools_fluxErrModels.c  sntools_imageNoiseModels.c
//  into 'legacy' models, to eventually be replaced by more
//  general set of maps. Current code is a mess because different
//  pieces are in different places.
//
// **************************************************
// **************************************************
// **************************************************

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_fluxErrModels_legacy.h"




// ===================================
//    sntools_fluxErrModels.c
// ===================================


/**************************************

  Oct 2015 R.Kessler

  Detailed model of flux error correction after simulation
  computes error based on ZP, PSF, SKYSIG. These routines
  are based on comparing sim-errors to data or to fakes
  passed through the same image pipeline as data.

  To ensure simlib syncing, the name of the simlib file
  is passed as the primary 'option' key. An additional 
  integer OPT flag is also passed to allow for systematic 
  variations.

  Oct 25 2017: update from Dillon Brout to work with DES_SMP

*********************************************************/



void init_fluxErrModel_legacy(char *SIMLIB_FILENAME, int OPTFLAG) {

  char fnam[] = "init_fluxErrModel_legacy" ;

  // -------------- BEGIN ---------------

  // must have OPTFLAG > 0 to continue
  if ( OPTFLAG <= 0 ) { return ; }

  sprintf(BANNER,"%s for %s and OPTFLAG=%d \n",
	  fnam, SIMLIB_FILENAME, OPTFLAG);

  print_banner(BANNER);

  // store arguments in global
  sprintf(FLUXERRMODEL.SIMLIB_FILENAME,"%s", SIMLIB_FILENAME);
  FLUXERRMODEL.OPTFLAG = OPTFLAG ;

}  // end of init_fluxErrModels_legacy



// =======================================================
double scale_fluxErrModel_legacy(char *BAND, char *FIELD, double MJD,
				 double ZP, double SKYSIG, double PSF) {

  // return scale to apply to simulated error.
  // In snlc_sim,   FLUXERR *= scale_fluxErrModel

  char fnam[] = "scale_fluxErrModel_legacy";
  char *ptrFile = FLUXERRMODEL.SIMLIB_FILENAME ;
  double scale ;

  // -------------- BEGIN ---------------

  // store inputs in global
  sprintf(FLUXERRMODEL.BAND,  "%s", BAND);
  sprintf(FLUXERRMODEL.FIELD, "%s", FIELD);
  FLUXERRMODEL.IBAND   = INTFILTER(FLUXERRMODEL.BAND);
  FLUXERRMODEL.MJD     = MJD ;  // days
  FLUXERRMODEL.ZP      = ZP ;
  FLUXERRMODEL.SKYSIG  = SKYSIG ;  // ADU per pixel
  FLUXERRMODEL.PSF     = PSF ;     // sigma, pixels

  scale = 0.0 ;

  if ( strstr(ptrFile,"DES_DIFFIMG") != NULL ) {
    scale =  scale_fluxErr_DES_DIFFIMG();
  }
  else if ( strstr(ptrFile,"DES_SMP") != NULL ) {
    scale = scale_fluxErr_DES_SMP();
  }
  else {    
    sprintf(c1err,"Unrecognized SIMLIB_FILENAME");
    errmsg(SEV_FATAL, 0, fnam, c1err, ptrFile ); 
  }

  return scale ;

} // end of scale_fluxErrModel_legacy


// ==========================================
double scale_fluxErr_DES_DIFFIMG(void) {

  // Nov 2015
  // Scale snana-sim errors for DES-DIffImg photometry,
  // a.k.a "force photo".

  int ISDEEP;
  int  IBAND  = FLUXERRMODEL.IBAND;
  char *FIELD = FLUXERRMODEL.FIELD; 
  double PSF  = FLUXERRMODEL.PSF ; // sigma, pixels 
  double SCALE ;

  // PAW > test#errfudge year=1+Y2 deep=0
  // PAW > test#errfudge year=1+Y2 deep=1

  double scale_psfCor[2][6][5] =  
    { // [isdeep][ifiltobs][polyCoeff]
      {
	{0.0,      0.0,       0.0,     0.0,       0.0}, // unused
	{1.0,      0.0,       0.0,       0,         0},          // u, shallow
	{4.53181, -4.33347,  1.9626,   -0.385538,  0.0280264},  // g
	{3.70768, -3.44778,  1.61242,  -0.326832,  0.0245651},  // r
	{2.31471, -1.67548,  0.773386, -0.1527,    0.0112755},  // i
	{1.60241, -0.829311, 0.417534, -0.0892309, 0.00715445}, // z
      },
      //
      {
	{0.0,     0.0,        0, 0, 0},   // unused
	{1.0,     0.0,        0, 0, 0},   // u, deep
	{1.08856, 0.0290076,  0, 0, 0},   // g
	{1.06293, 0.020909,   0, 0, 0},   // r
	{1.03998, 0.0442453,  0, 0, 0},   // i
	{1.17307, 0.00649704, 0, 0, 0}    // z
      }
    } ;

  int i;
  double ppow[5] ;
    
  // --------------- BEGIN ---------------

  // determine DEEP or shallow field
  ISDEEP = 0 ;
  if ( strcmp(FIELD,"X3") == 0 ) { ISDEEP = 1; }
  if ( strcmp(FIELD,"C3") == 0 ) { ISDEEP = 1; }

  ppow[0]  =  1.0;
  ppow[1]  =  PSF;
  ppow[2]  =  PSF*PSF ;
  ppow[3]  =  PSF*PSF*PSF ;
  ppow[4]  =  ppow[2] * ppow[2] ;   

  SCALE = 0.0 ;
  for(i = 0; i < 5; i++ ) {
    SCALE += ppow[i] * scale_psfCor[ISDEEP][IBAND][i];
  }

  // Mar 7 2016: scale shallow field errors by another 10% ;
  //             Reason is not understood.
  //  if ( ISDEEP == 0 ) {  SCALE *= 1.1 ;  }


  /*
  printf(" xxx BAND=%s(%d) FIELD=%s ISDEEP=%d   PSF=%.3f  SCALE=%.2f \n",
	 BAND, IBAND, FIELD, ISDEEP, PSF, SCALE); // xxxxxx
  */

  return SCALE ;

} // end scale_fluxErr_DES_DIFFIMG


// ==========================================
double scale_fluxErr_DES_SMP(void) {

  // Oct 2017 - Dillon Brout
  // Scale snana-sim errors for DES-SMP photometry,
  

  int ISDEEP;
  int  IBAND  = FLUXERRMODEL.IBAND;
  char *FIELD = FLUXERRMODEL.FIELD;
  double PSF  = FLUXERRMODEL.PSF ; // sigma, pixels
  double SCALE ;

  
  // define scale on 1/ERROR
  double scale_psfCor[2][6][5] =
    { // [isdeep][ifiltobs][polyCoeff]
      {
	{0.0,      0.0,       0.0,     0.0,       0.0}, // unused
	{1.0,      0.0,       0.0,       0,         0},          // u, shallow
	{0.0607565362076, -0.544855960539, 1.64207510476, -0.523081799288,  0},  // g
	{0.034007152959, -0.279534495415, 0.798301326773, 0.309383698661,  0},  // r
	{-0.0731411268374, 0.535407222377, -1.19447182689, 1.86355683887, 0},  // i
	{0.0458880745923, -0.312013858808, 0.74115602632, 0.472181059646, 0}, // z
      },
      //
      {
	{0.0,     0.0,        0, 0, 0},   // unused
	{1.0,     0.0,        0, 0, 0},   // u, deep
	{-0.0546075051192, 0.240486543405, 0.136561715417, 0.223413716625, 0},   // g
	{0.149078822923, -1.20998524892, 3.25822786077, -1.76253852481, 0},   // r
	{0.047048994048, -0.360533841825, 0.991150247057, 0.205524045574, 0},   // i
	{-0.19138434084, 1.35164735143, -2.79634068599, 2.94513446907, 0}    // z
      }
    } ;
  
  int i;
  double ppow[5] ;
    
  // --------------- BEGIN ---------------
  
  // determine DEEP or shallow field
  ISDEEP = 0 ;
  if ( strcmp(FIELD,"X3") == 0 ) { ISDEEP = 1; }
  if ( strcmp(FIELD,"C3") == 0 ) { ISDEEP = 1; }
  
  ppow[0]  =  PSF*PSF*PSF ;
  ppow[1]  =  PSF*PSF;
  ppow[2]  =  PSF ;
  ppow[3]  =  1.0 ;
  ppow[4]  =  ppow[2] * ppow[2] ;
  
  SCALE = 0.0 ;
  for(i = 0; i < 5; i++ ) {
    SCALE += ppow[i] * scale_psfCor[ISDEEP][IBAND][i];
  }


  return (1./SCALE) ;
    
} // end scale_fluxErr_DES_SMP



// ======================================
//     sntools_imageNoiseModels.c
// ======================================

/**********************************


  Aug 2014

  Models of anomalous noise in the data
  (i.e., NOT photon statistics and NOT intrinsic scatter).


 **********************************/


// ********************************
void INIT_NOISEMODEL_HOST_LEGACY(char *HOSTNOISE_FILE) {

  // Created Aug 2014 by R.Kessler
  // Read HOSTNOISE_FILE (if defined) and setup quick lookup maps.
  // This is anomalous subtraction noise in pre-explosion epochs
  // that is correlated with the observed galaxy mag.
  // This noise it is NOT the photon shot-noise from 
  // the HOSTLIB routines in sntools_host.c[h].
  //
  // Feb 21 2018: initialize each new map to have FIELDLIST=ALL

  int  IMAP, IBIN, LIBID, gzipFlag ;
  FILE *fp ;
  char  msg[100];
  char *ptrFile = HOSTNOISE_FILE; 
  char fnam[]   = "INIT_NOISEMODEL_HOST_LEGACY" ;
  char c_get[60];
  char PATH_SIMLIB[MXPATHLEN] ;

  // ------------- BEGIN ----------------

  NOISEMODEL_NAME[0]   = 0 ;
  NMAP_NOISEMODEL_HOST = 0 ;
  NPAR_NOISEMODEL      = 0;

  if ( IGNOREFILE(ptrFile)  ) { return ; }

  // xxx  sprintf(NOISEMODEL_FILE, "%s", HOSTNOISE_FILE); // store in global

  sprintf(msg,"%s: init model for anomalous image noise from host."
	  ,fnam);
  print_banner(msg);

  // use utility to check local dir and path.        
  sprintf(PATH_SIMLIB, "%s/simlib", PATH_SNDATA_ROOT);
  fp = snana_openTextFile(1,PATH_SIMLIB, ptrFile,
			  NOISEMODEL_FILE, &gzipFlag ); // returned (and fill global)

  if ( !fp ) {
    abort_openTextFile("HOSTNOISE_FILE", PATH_SIMLIB, ptrFile, fnam);

    /* xxxxxxxx mark delete Feb 1 2020 xxxxxxxx
    sprintf(c1err, "could not open HOSTNOISE_FILE:");
    sprintf(c2err, "%s", ptrFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    xxxxxxxx */
  }
  else {
    //    printf("\t Read noise model from: \n\t %s\n", NOISEMODEL_FILE );
  }

  // ----------------

  IMAP = -9;

  while( (fscanf(fp, "%s", c_get)) != EOF) {
    
    if ( strcmp(c_get,"NOISEMODEL_NAME:") == 0 )  { 
      readchar(fp, NOISEMODEL_NAME);  
      NPAR_NOISEMODEL = 1;
    }

    if ( strcmp(c_get,"LIBID:") == 0 ) {
      // start next map of noise params vs. HOSTMAG
      NMAP_NOISEMODEL_HOST++ ;  
      if ( NMAP_NOISEMODEL_HOST >= MXMAP_NOISEMODEL_HOST ) {
	sprintf(c1err, "NMAP_NOISEMODEL_HOST = %d", 
		NMAP_NOISEMODEL_HOST);
	sprintf(c2err, "in LIBID after LIBID=%d", 
		NOISEMODEL_HOST[IMAP].LIBID );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      IMAP = NMAP_NOISEMODEL_HOST-1 ; 
      readint(fp, 1, &LIBID);
      NOISEMODEL_HOST[IMAP].LIBID  = LIBID ;
      NOISEMODEL_HOST[IMAP].NBIN_HOSTMAG = 0 ;
      sprintf(NOISEMODEL_HOST[IMAP].FIELDLIST,"ALL");
    }

    if ( strcmp(c_get,"BAND:") == 0 ) {
      readchar(fp, NOISEMODEL_HOST[IMAP].BANDLIST );
    }

    if ( strcmp(c_get,"FIELD:") == 0 ) {
      readchar(fp, NOISEMODEL_HOST[IMAP].FIELDLIST );
    }

    if ( strcmp(c_get,"HOSTNOISE:") == 0 ) {
      IBIN = NOISEMODEL_HOST[IMAP].NBIN_HOSTMAG ;  // hostmag bin index
      readdouble(fp, 1,      
		 &NOISEMODEL_HOST[IMAP].HOSTMAG[IBIN] );
      readdouble(fp, NPAR_NOISEMODEL, 
		 NOISEMODEL_HOST[IMAP].PARAM[IBIN] );

      NOISEMODEL_HOST[IMAP].NBIN_HOSTMAG++ ;
      
      IBIN = NOISEMODEL_HOST[IMAP].NBIN_HOSTMAG ;  // hostmag bin index
      if ( IBIN >= MXBIN_HOSTMAG) {
	sprintf(c1err, "IBIN=%d for LIBID=%d", IBIN, LIBID);
	sprintf(c2err, "Check %s", ptrFile);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
    }

  } // end while

  // ------------------

  for(IMAP=0; IMAP < NMAP_NOISEMODEL_HOST ; IMAP++ ) {

    printf("\t IMAP=%2d: LIBID=%d   BANDLIST=%s   FIELDLIST=%s \n",
	   IMAP,
	   NOISEMODEL_HOST[IMAP].LIBID,
	   NOISEMODEL_HOST[IMAP].BANDLIST,
	   NOISEMODEL_HOST[IMAP].FIELDLIST );
    fflush(stdout);
  }


  return ;

} // end of INIT_NOISEMODEL_HOST_LEGACY


void  init_noisemodel_host_legacy__(char *HOSTNOISE_FILE) 
{ INIT_NOISEMODEL_HOST_LEGACY(HOSTNOISE_FILE); }


void gen_noisemodel_host_legacy__(char *BAND,  char *FIELD, int *GALID, 
				  double *GALMAG, double *SBMAG, 
				  double *SNSEP, double *noisePar) {
  GEN_NOISEMODEL_HOST_LEGACY(BAND,FIELD, *GALID, *GALMAG, *SBMAG, 
			     *SNSEP, noisePar);
}


// ====================================================
void GEN_NOISEMODEL_HOST_LEGACY(char *BAND,  char *FIELD, 
			 int GALID, double GALMAG, double SBMAG, double SNSEP,
			 double *noisePar ) {

  // return anomalous host-subtraction noise in FLUXCAL units.
  // Inputs:
  //  BAND :   char representation, e.g., 'r'
  //  FIELD:   name of field
  //  GALID:   hostlib identifier (for debug statements only)
  //  GALMAG:  host galaxy mag in this band (includes total flux)
  //  SBMAG:   surface brightness mag in this band, in 1 sq arcsec.
  //  SNSEP:   SN-host sep, arcsec

  int    LDMP = 0 ;
  int    imap, IMAP, NMATCH, LIBID, NBIN, ibin, IPAR, NPAR ;
  char   BANDLIST[MXPATHLEN], FIELDLIST[MXPATHLEN];
  double MAG_LOCAL ;
  double GRID_HOSTMAG[MXBIN_HOSTMAG], GRID_NOISE[MXBIN_HOSTMAG] ;

  char   fnam[] = "GEN_NOISEMODEL_HOST_LEGACY" ;

  // ----------------- BEGIN --------------------

  NPAR = 2;
  noisePar[0] = 0.0 ;  // FLUXCAL noise per pixe;
  noisePar[1] = 0.0 ;  // ERRSCALE to multiply skynoise

  // figure out which LIBID based on BAND and FIELD

  NMATCH = 0 ;  LIBID = IMAP = -9;
  for(imap=0; imap < NMAP_NOISEMODEL_HOST ; imap++ ) {
    
    sprintf(BANDLIST,  "%s", NOISEMODEL_HOST[imap].BANDLIST  ) ;
    sprintf(FIELDLIST, "%s", NOISEMODEL_HOST[imap].FIELDLIST ) ;

    if ( strstr(BANDLIST, BAND)  == NULL ) { continue ; }
    if ( strstr(FIELDLIST,FIELD) == NULL ) { continue ; }

    NMATCH++ ;
    IMAP  = imap ;
    LIBID = NOISEMODEL_HOST[IMAP].LIBID ;
  }
  
  if ( NMATCH != 1 ) {
    sprintf(c1err, "Invalid NMATCH=%d for BAND=%s  FIELD=%s",
	    NMATCH, BAND, FIELD);
    sprintf(c2err, "Check file: %s", NOISEMODEL_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }

  if ( IMAP < 0 ) {
    sprintf(c1err, "Invalid IMAP=%d (BAND=%s FIELD=%s", 
	    IMAP, BAND, FIELD );
    sprintf(c2err, "Check file: %s", NOISEMODEL_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }

  if ( LDMP ) {
    printf(" xxx ---------------------------------------------- \n");
    printf(" xxx %s: FIELD=%s-%s   MAG(GAL->SB)=%.2f->%.2f  GALID=%d\n", 
	   fnam, FIELD, BAND, GALMAG, SBMAG, GALID );
    fflush(stdout);
  }

  // interpolate mag to find noise, FLUXCAL per pixel
  NBIN = NOISEMODEL_HOST[IMAP].NBIN_HOSTMAG ;
  for(ibin=0; ibin < NBIN; ibin++ ) {
    GRID_HOSTMAG[ibin] = NOISEMODEL_HOST[IMAP].HOSTMAG[ibin] ;
    GRID_NOISE[ibin]   = NOISEMODEL_HOST[IMAP].PARAM[ibin][0] ;

    //    GRID_NOISE[ibin] += 1.0 ; // idiot test xxxxxxxxxxxxx

    if ( LDMP ) {
      printf(" xxx bin=%2d  HOSTMAG=%f  NOISE=%f \n",
	     ibin, GRID_HOSTMAG[ibin], GRID_NOISE[ibin] );
      fflush(stdout);
    }
  } // end of ibin
  

  // -------------------------------------------------------
  // if GALMAG is outside the grid, use the edge grid value

  if ( strcmp(NOISEMODEL_NAME,"SB_ERRSCALE") == 0 ) {
    MAG_LOCAL = SBMAG ;
    IPAR = 1; // index of noisePar to fill
  }
  else {
    MAG_LOCAL = 99. ;
    sprintf(c1err, "Unknown NOISEMODEL_NAME = '%s'", NOISEMODEL_NAME);
    sprintf(c2err, "Check file: %s", NOISEMODEL_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }


  int ibin_first = 0;
  int ibin_last  = NBIN - 1;
  int OPT_INTERP = 1 ; // linear
  double TMP, TMP0, TMP1, MAG0, MAG1, SLOPE ;

  if ( MAG_LOCAL < GRID_HOSTMAG[ibin_first] )  {  
    // extrapolate at bright end
    TMP0 = GRID_NOISE[ibin_first] ; 
    TMP1 = GRID_NOISE[ibin_first+1] ;
    MAG0 = GRID_HOSTMAG[ibin_first] ; 
    MAG1 = GRID_HOSTMAG[ibin_first+1] ; 
    SLOPE = (TMP1-TMP0)/(MAG1-MAG0);
    TMP   = GRID_NOISE[ibin_first] + (MAG_LOCAL-MAG0)*SLOPE ;
  }

  else if ( MAG_LOCAL > GRID_HOSTMAG[ibin_last] )  {  
    // at dim end, use value in last bin (i.e., no extrapolation)
    TMP = GRID_NOISE[ibin_last] ;  
  }
  else {
    // if we get here, then interpolate GRID_NOISE
    TMP = interp_1DFUN(OPT_INTERP, MAG_LOCAL, NBIN, 
		       GRID_HOSTMAG, GRID_NOISE, fnam );
  }

  
  // do not let ERRSCALE fall below 1
  if ( IPAR == 1 && TMP < 1.0 ) { TMP = 1.0 ; }

  noisePar[IPAR] = TMP ;

  if ( LDMP ) {
    printf(" xxx found LIBID=%d  NoisePar=%f \n",  LIBID, TMP );
    fflush(stdout);
  }

  return ;

} // end of GEN_NOISEMODEL_HOST_LEGACY
