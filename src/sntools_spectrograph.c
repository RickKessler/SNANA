/**************************************************
 Created July 2016 by R.Kessler

 utility to read info needed to simulation spectra and/or IFU.
 Used by kcor and snlc_sim programs.

 Strategy:
   User provides table of wavelength bins, and in each bin gives
   SNR1 & SNR2 for the two MAGREF values. This code then solves for 
   effective zeroPoint (Zpe) and skyNoise (includes readout), 
   which can be applied to any MAG.  Each pair of SNR values can 
   be given for several different exposure times,
   allowing interpolation between exposure times. This feature
   allows chaning the exposure time for SN at different redshifts.

  To use this function, kcor input file should contain a key
     SPECTROGRAPH:  <fileName>

  where <fileName> contains

  INSTRUMENT:  <name>
  MAGREF_LIST:  <magref1>  <magref2>
  TEXPOSE_LIST: <t_1>  <t_2> ... <t_n>
 
  SPECBIN: <minL> <maxL>  <sigL> SNR1(t_1) SNR2(t_1) . . SNR1(t_n) SNR2(t_n)
  SPECBIN: <minL> <maxL>  <sigL> SNR1(t_1) SNR2(t_1) . . SNR1(t_n) SNR2(t_n)
  SPECBIN: <minL> <maxL>  <sigL> SNR1(t_1) SNR2(t_1) . . SNR1(t_n) SNR2(t_n)

  etc ...

  Note that minLam,maxLam are specified for each bin so that
  non-uniform bins are allowed.

           HISTORY 
       ~~~~~~~~~~~~~~

  May 06 2020:
    + in getSNR_spectrograph, return SNR=0 if ZP is undefined.

  Aug 20 2021:
    + in getSNR_spectrograph, interpolate ZP vs. log10(Texpose) instead
      of ZP vs. Texpose. Works much better with sparse Texpose grid.
      [issue found by comparing SNR against D.Rubin]

*********************************************************/

#include "fitsio.h"
#include "sntools.h"
#include "sntools_spectrograph.h"
#include "sntools_dataformat_fits.h"

// =======================================
void init_spectrograph(char *inFile, char *stringOpt ) {

  //  char fnam[] = "init_spectrograph" ;

  // ------------- BEGIN ---------------

  // check option(s) 
  parse_spectrograph_options(stringOpt);

  // read wavelenth bins and list of SNR1,SNR2 from ascii file
  read_spectrograph_text(inFile);

  // solve for ZP and SIGSKY for each SNR[1,2] pair
  solve_spectrograph(); 

  printf("\n"); fflush(stdout);

} // end init_spectrograph

// ===============================================
void parse_spectrograph_options(char *stringOpt) {

#define MXOPT_SPEC 10
  int  opt, NOPT, NTMP ; 
  char *ptrOpt[MXOPT_SPEC], s[MXOPT_SPEC][40];
  char *ptr2[MXOPT_SPEC], s2[2][40];
  char comma[]=",", equal[]="=" ;
  //   char fnam[] = "parse_spectrograph_options" ;

  // ------------- BEGIN --------------

  INPUTS_SPECTRO.NREBIN_LAM=1 ;

  if ( IGNOREFILE(stringOpt) ) { return ; }

  // split by comma-separated values
  for(opt=0; opt < MXOPT_SPEC; opt++ )  { ptrOpt[opt] = s[opt] ; }

  splitString(stringOpt, comma, MXOPT_SPEC,           // inputs               
              &NOPT, ptrOpt );                        // outputs             

  // split each option by equal sign:  bla=val
  ptr2[0] = s2[0];
  ptr2[1] = s2[1];
  for(opt=0; opt < NOPT; opt++ ) {
    splitString(s[opt], equal, MXOPT_SPEC, &NTMP, ptr2);

    if ( strcmp(s2[0],"rebin")==0 ) {
      sscanf(s2[1] , "%d", &INPUTS_SPECTRO.NREBIN_LAM ); 
      printf("\t Spectrograph option: rebin wavelength by %d\n",
	     INPUTS_SPECTRO.NREBIN_LAM) ;
    }
  }

  return ;

} // end parse_spectrograph_options

// ===============================================
void read_spectrograph_text(char *inFile) {

  FILE *fp ;
  int NROW_FILE, ikey, NKEY_FOUND, DONE_MALLOC ;
  int NBL, NBT, NRDCOL, t ;
  char c_get[60], tmpLine[200] ;
  char fnam[] = "read_spectrograph_text" ;

#define NKEY_REQ_SPECTROGRAPH 3
#define IKEY_INSTRUMENT 0
#define IKEY_MAGREF     1
#define IKEY_TEXPOSE    2

  const char KEYREQ_LIST[NKEY_REQ_SPECTROGRAPH][40] =  { 
    "INSTRUMENT:", 
    "MAGREF_LIST:", 
    "TEXPOSE_LIST:"  
  } ;
  int KEYFLAG_FOUND[NKEY_REQ_SPECTROGRAPH];

  // ---------------- BEGIN -------------

  printf("\n %s: \n", fnam); fflush(stdout);

  // get number of rows in file to get malloc size.
  NROW_FILE = nrow_read(inFile,fnam);

  fp = fopen(inFile,"rt");
  if ( !fp ) {
    sprintf(c1err,"Could not find SPECTROGRAPH table file:");
    sprintf(c2err,"%s", inFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  sprintf(INPUTS_SPECTRO.INFILE_NAME,      "inFile"  );
  sprintf(INPUTS_SPECTRO.INSTRUMENT_NAME,  "UNKNOWN" );
  INPUTS_SPECTRO.MAGREF_LIST[0] = 99. ;
  INPUTS_SPECTRO.MAGREF_LIST[1] = 99. ;
  INPUTS_SPECTRO.NBIN_TEXPOSE = 0 ;
  INPUTS_SPECTRO.NBIN_LAM     = 0 ;
  INPUTS_SPECTRO.NBIN_LAM_noREBIN = 0 ;  
  INPUTS_SPECTRO.LAM_MIN     = INPUTS_SPECTRO.LAM_MAX     = -9.0 ;
  INPUTS_SPECTRO.TEXPOSE_MIN = INPUTS_SPECTRO.TEXPOSE_MAX = -9.0 ;

  INPUTS_SPECTRO.SNR_POISSON_RATIO_ABORT_vsMAGREF  = 5.0 ;
  INPUTS_SPECTRO.SNR_POISSON_RATIO_ABORT_vsTEXPOSE = 1.2 ;
  INPUTS_SPECTRO.MAGSNR_TOLERANCE_ABORT            = 0.001 ;

  NERR_SNR_SPECTROGRAPH = 0 ;
  NERR_BADSNR_SPECTROGRAPH = 0 ;

  NKEY_FOUND = DONE_MALLOC = 0 ;
  for(ikey=0; ikey < NKEY_REQ_SPECTROGRAPH; ikey++ ) 
    { KEYFLAG_FOUND[ikey] = 0; }

  NBL = NBT = 0 ;

  printf("    Open %s \n", inFile ); fflush(stdout);

  reset_VALUES_SPECBIN(); 

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    // if comment key is found, read remainder of line into dummy string  
    // so that anything after comment key is ignored (even a valid key)  
    if ( c_get[0] == '#' || c_get[0] == '!' || c_get[0] == '%' )
      { fgets(tmpLine, 80, fp) ; continue ; }

    // check all default keys
    for(ikey=0; ikey < NKEY_REQ_SPECTROGRAPH; ikey++ ) {
      if ( strcmp(c_get, KEYREQ_LIST[ikey] ) == 0 ) { 
	NKEY_FOUND++ ;  KEYFLAG_FOUND[ikey] = 1;

	switch(ikey) {
	case IKEY_INSTRUMENT :
	  readchar(fp, INPUTS_SPECTRO.INSTRUMENT_NAME) ; break ;
	case IKEY_MAGREF :
	  readdouble(fp, 2, INPUTS_SPECTRO.MAGREF_LIST ) ; break ;
	case IKEY_TEXPOSE :
	  NBT = read_TEXPOSE_LIST(fp); break ;
	default:
	  sprintf(c1err,"Unknow KEY[%d] = %s", ikey, KEYREQ_LIST[ikey] );
	  sprintf(c2err,"Check code and file=%s", inFile );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
	  break ;
	}

      } // end string match
    } // end ikey

    // - - - - - - - - - - - - - - - - - - 
    if( NKEY_FOUND == NKEY_REQ_SPECTROGRAPH && DONE_MALLOC==0 ) {
      // malloc array vs. lambda and TEXPOSE
      malloc_spectrograph(+1, NROW_FILE, NBT);

      NRDCOL = 2 + 2*NBT ; // number of columns to read below
      DONE_MALLOC = 1;
    }

    // check optional key(s)
    if ( strcmp(c_get,"SNR_POISSON_RATIO_ABORT_vsMAGREF:") == 0 ) {
      readdouble(fp, 1, &INPUTS_SPECTRO.SNR_POISSON_RATIO_ABORT_vsMAGREF);
    }

    if ( strcmp(c_get,"SNR_POISSON_RATIO_ABORT_vsTEXPOSE:") == 0 ) {
      readdouble(fp, 1, &INPUTS_SPECTRO.SNR_POISSON_RATIO_ABORT_vsTEXPOSE);
    }
    
    if ( strcmp(c_get,"MAGSNR_TOLERANCE_ABORT:") == 0 ) {
      readdouble(fp, 1, &INPUTS_SPECTRO.MAGSNR_TOLERANCE_ABORT );
    }

    // - - - - - - - - - -

    if ( strcmp(c_get,"SPECBIN:") == 0 ) {

      if ( NKEY_FOUND < NKEY_REQ_SPECTROGRAPH ) { 
	print_preAbort_banner(fnam);
	for(ikey=0; ikey < NKEY_REQ_SPECTROGRAPH; ikey++ ) {
	  printf("   Required header key:  '%s'   (FOUND=%d) \n", 
		 KEYREQ_LIST[ikey], KEYFLAG_FOUND[ikey] );
	}
	sprintf(c1err,"Found SPECBIN key before required header keys.");
	sprintf(c2err,"Read %d of %d required keys.",
		NKEY_FOUND, NKEY_REQ_SPECTROGRAPH );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
      } // end ABORT

      NBL = read_SPECBIN_spectrograph(fp);

    } // end SPECBIN line

  }  // end while
  

  if ( NERR_SNR_SPECTROGRAPH > 0 ) {
    sprintf(c1err,"Found %d errors for which", NERR_SNR_SPECTROGRAPH);
    sprintf(c2err,"SNR(Texpose) is NOT monotically increasing");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( NERR_BADSNR_SPECTROGRAPH > 0 ) {
    sprintf(c1err,"Found %d SNR<=0 errors.", NERR_BADSNR_SPECTROGRAPH);
    sprintf(c2err,"Check spectrograph table.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // ---------------------------------------

  printf("    Read %d LAMBDA bins from %.0f to %.0f A \n",
	 NBL, INPUTS_SPECTRO.LAM_MIN, INPUTS_SPECTRO.LAM_MAX ); 

  printf("    Read %d TEXPOSE values: ", NBT ); 
  for(t=0; t < NBT; t++ ) 
    { printf("%d ", (int)INPUTS_SPECTRO.TEXPOSE_LIST[t]); }
  printf(" sec \n");

  printf("    SNR_POISSON_RATIO_ABORT(MAGREF,TEXPOSE) = %.2f, %.2f \n",
	 INPUTS_SPECTRO.SNR_POISSON_RATIO_ABORT_vsMAGREF,
	 INPUTS_SPECTRO.SNR_POISSON_RATIO_ABORT_vsTEXPOSE) ;
 
  printf("    Abort if [(m1-m0) - 2.5*log10(SNR0/SNR1)] < %f \n",
	 INPUTS_SPECTRO.MAGSNR_TOLERANCE_ABORT );
 

  // -----
  if ( NBL >= MXLAM_SPECTROGRAPH ) {
    sprintf(c1err,"NBIN_LAM=%d exceeds MXLAM_SPECTROGRAPH=%d",
	    NBL, MXLAM_SPECTROGRAPH );
    sprintf(c2err,"Check spectrograph file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  fflush(stdout);
  fclose(fp);

  return ;

} // end read_spectrograph_text

// ========================================================
int read_SPECBIN_spectrograph(FILE *fp) {

  // read values after SPECBIN key

  int NBT    = INPUTS_SPECTRO.NBIN_TEXPOSE ;
  int NBL    = INPUTS_SPECTRO.NBIN_LAM ;
  int NRDCOL = NCOL_noSNR + 2*NBT ;  // LAMMIN,LAMMAX,LAMRES + SNR values
  int NREBIN_LAM = INPUTS_SPECTRO.NREBIN_LAM ;
  int jtmp, t, NBMOD ;
  double XTMP[MXVALUES_SPECBIN] ;
  double SNR_OLD, SNR_NEW ;
  char fnam[] = "read_SPECBIN_spectrograph" ;

  // ------------ BEGIN -------------

  readdouble(fp, NRDCOL, XTMP );
  if ( NREBIN_LAM == 1 ) { goto PARSE_XTMP ; }

  // ----------------------
  // check rebinning
  INPUTS_SPECTRO.NBIN_LAM_noREBIN++ ;
  NBMOD = INPUTS_SPECTRO.NBIN_LAM_noREBIN%NREBIN_LAM ;
  if ( NBMOD == 1 )  { VALUES_SPECBIN[0] = XTMP[0] ; }
  VALUES_SPECBIN[1] = XTMP[1] ;

  // always increment SNR in quadrature
  for(jtmp=NCOL_noSNR; jtmp < NRDCOL; jtmp++ ) {
    SNR_OLD = VALUES_SPECBIN[jtmp] ;
    SNR_NEW = XTMP[jtmp];
    VALUES_SPECBIN[jtmp] = sqrt(SNR_OLD*SNR_OLD + SNR_NEW*SNR_NEW);
  }
  if ( NBMOD != 0 ) { return(NBL); }

  // -----------------------------------------
  // transfer global VALUES_SPECBIN back to local XTMP
  for(jtmp=0; jtmp < NRDCOL; jtmp++ )
    { XTMP[jtmp] = VALUES_SPECBIN[jtmp] ; }

 PARSE_XTMP:

  if ( XTMP[1] < XTMP[0] ) {
    sprintf(c1err,"LAMMAX=%f < LAMMIN=%f ???", XTMP[1], XTMP[0] );
    sprintf(c2err,"Check LAMBDA binning in SPECTROGRAPH file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }
 

  INPUTS_SPECTRO.LAMMIN_LIST[NBL] = XTMP[0] ;
  INPUTS_SPECTRO.LAMMAX_LIST[NBL] = XTMP[1] ;
  INPUTS_SPECTRO.LAMSIGMA_LIST[NBL] = XTMP[2] ;
  INPUTS_SPECTRO.LAMAVG_LIST[NBL] = ( XTMP[0] + XTMP[1] ) / 2.0 ;
  INPUTS_SPECTRO.LAMBIN_LIST[NBL] = ( XTMP[1] - XTMP[0] ) ;

  jtmp=NCOL_noSNR ;
  for(t=0; t < NBT; t++ ) {
    INPUTS_SPECTRO.SNR0[NBL][t]  = XTMP[jtmp] ;  jtmp++ ;
    INPUTS_SPECTRO.SNR1[NBL][t]  = XTMP[jtmp] ;  jtmp++ ;
    check_SNR_SPECTROGRAPH(NBL,t);
  }
  NBL++ ;  
  INPUTS_SPECTRO.NBIN_LAM = NBL ;


  // keep track of global min/max
  if ( NBL == 1 ) { INPUTS_SPECTRO.LAM_MIN = XTMP[0]; }
  INPUTS_SPECTRO.LAM_MAX = XTMP[1];
  
  reset_VALUES_SPECBIN(); 

  return(NBL) ;

} // end read_SPECBIN_spectrograph


// ========================================================
void check_SNR_SPECTROGRAPH(int l, int t) {

  // May 21 2020: Refactor and update
  //
  // ABORT if
  //   SNR(t) is not increasing
  //   SNR1/SNR0 ratio is outside tolerance 
  //   SNR(t)/SNR(t-1) is outside tolerance
  //
  // SNR-Tolerance ratios are user-input in spectrograph file:
  //   SNR_POISSON_RATIO_ABORT_vsMAGREF:  <val>
  //   SNR_POISSON_RATIO_ABORT_vsTEXPOSE: <val>
  //
  // Inputs:
  //   t = TEXPOSURE index
  //   l = lambda index
  //
  //
  // BEWARE SNR INDEX:
  //   In code we have SNR0,SNR1, but manual has SNR1,SNR2 ..
  //   so print error messages to match manual notation.
  //

  double LAM     = INPUTS_SPECTRO.LAMAVG_LIST[l];
  double TEXPOSE = INPUTS_SPECTRO.TEXPOSE_LIST[t] ;
  double MAGREF0 = INPUTS_SPECTRO.MAGREF_LIST[0];
  double MAGREF1 = INPUTS_SPECTRO.MAGREF_LIST[1];
  double SNR0    = INPUTS_SPECTRO.SNR0[l][t];
  double SNR1    = INPUTS_SPECTRO.SNR1[l][t];
  double LOGSNR_RATIO_TOL;

  double arg, FLUXREF_RATIO, SNR_RATIO_POISSON, SNR_RATIO, LOGSNR_RATIO ;
  double *ptrSNR;
  double *ptrTexpose = INPUTS_SPECTRO.TEXPOSE_LIST;
  int    iSNR ;
  char fnam[]    = "check_SNR_SPECTROGRAPH" ;
  //  int  LDMP = 0;
  // ---------------BEGIN ----------

    sprintf(c2err,"LAM=%.1f  TEXPOSE=%.2f  (l=%d, t=%d)", 
	    LAM, TEXPOSE, l, t );

  // check for negative SNR that are obviously bad
  if ( SNR0 < 1.0E-12 || SNR1 < 1.0E-12 ) { 
    sprintf(c1err,"Invalid/negative SNR[0,1] = %le, %le ", SNR0, SNR1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // check that SNR1/SNR0 ratio is with factor of TOL of Poisson expectations
  arg               = 0.4*(MAGREF0-MAGREF1);
  FLUXREF_RATIO     = pow(TEN,arg);  // fluxref1/fluxref0
  SNR_RATIO_POISSON = sqrt(FLUXREF_RATIO); // expected SNR ratio from Poisson
  SNR_RATIO         = SNR1/SNR0;  
  LOGSNR_RATIO = log10(SNR_RATIO/SNR_RATIO_POISSON);

  LOGSNR_RATIO_TOL = log10(INPUTS_SPECTRO.SNR_POISSON_RATIO_ABORT_vsMAGREF);
  if ( fabs(LOGSNR_RATIO) > LOGSNR_RATIO_TOL ) {
    print_preAbort_banner(fnam);
    printf("\t MAGREF[0,1]   = %.3f, %.3f \n", MAGREF0, MAGREF1);
    printf("\t FLUXREF_RATIO = %le \n", FLUXREF_RATIO);
    printf("\t SNR1/SNR0[real,Poisson] = %f, %f \n", 
	   SNR_RATIO, SNR_RATIO_POISSON);
    printf("\t LOGSNR_RATIO = %f (TOL=%.3f) \n",
	   LOGSNR_RATIO, LOGSNR_RATIO_TOL) ;
    sprintf(c1err,"Possible crazy SNR1/SNR0 ratio (way off from Poisson)");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  // - - - - - - - - - 
  if ( t == 0 ) { return ; }

  LOGSNR_RATIO_TOL = log10(INPUTS_SPECTRO.SNR_POISSON_RATIO_ABORT_vsTEXPOSE);

  // check that SNR is monotonically increasig
  for(iSNR=0; iSNR < 2; iSNR++ ) {
    if ( iSNR == 0 ) { ptrSNR = INPUTS_SPECTRO.SNR0[l] ; }
    else             { ptrSNR = INPUTS_SPECTRO.SNR1[l] ; }

    SNR_RATIO = ptrSNR[t] / ptrSNR[t-1] ;
    SNR_RATIO_POISSON = sqrt(ptrTexpose[t]/ptrTexpose[t-1]); // prediction
    LOGSNR_RATIO = log10(SNR_RATIO/SNR_RATIO_POISSON);

    // abort if SNR is way off of Poisson prediction
    if ( fabs(LOGSNR_RATIO) > LOGSNR_RATIO_TOL ) {
      print_preAbort_banner(fnam);      
      printf("\t SNR%d(Texpose=%.2f) = %f (t=%d)\n", 
	     iSNR+1, ptrTexpose[t-1], ptrSNR[t-1], t-1);
      printf("\t SNR%d(Texpose=%.2f) = %f (t=%d)\n", 
	     iSNR+1, ptrTexpose[t],  ptrSNR[t], t );
      printf("\t SNR1/SNR0[real,Poisson] = %f, %f \n", 
	     SNR_RATIO, SNR_RATIO_POISSON);
      printf("\t LOGSNR_RATIO = %f (TOL=%.3f) \n",
	     LOGSNR_RATIO, LOGSNR_RATIO_TOL) ;
      sprintf(c1err,"Crazy SNR%d vs. Texpose", iSNR+1);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // abort if SNR is not increasing
    if ( ptrSNR[t] <= ptrSNR[t-1] ) {
      print_preAbort_banner(fnam);      
      printf("\t SNR%d(Texpose=%.2f) = %f (t=%d)\n", 
	     iSNR+1, ptrTexpose[t-1], ptrSNR[t-1], t-1);
      printf("\t SNR%d(Texpose=%.2f) = %f (t=%d)\n", 
	     iSNR+1, ptrTexpose[t],  ptrSNR[t], t );
      sprintf(c1err,"SNR%d is not increasing with Texpose", iSNR+1);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  }



  return ;
;
} // end check_SNR_SPECTROGRAPH

// ========================================================
void reset_VALUES_SPECBIN(void) {

  int jtmp;
  for(jtmp=0; jtmp < MXVALUES_SPECBIN; jtmp++) 
    { VALUES_SPECBIN[jtmp] = 0.0; } 

} // end reset_VALUES_SPECBIN

// ========================================================
int read_TEXPOSE_LIST(FILE *fp) {

  // read all exposure times on this line; 
  // stop reading at end-of-line or at comment string

#define MXCHAR_TEXPOSE_LIST 400  // max chars to read on line
  int NBT=0 ;
  double T, TLAST;
  char tmpLine[MXCHAR_TEXPOSE_LIST], *ptrtok, cval[40] ;
  char fnam[] = "read_TEXPOSE_LIST";

  // ---------- BEGIN --------------

  fgets(tmpLine, MXCHAR_TEXPOSE_LIST, fp) ; 
  ptrtok = strtok(tmpLine," ") ; // split string   

  while ( ptrtok != NULL  ) {
    sprintf(cval, "%s", ptrtok );
    
    //    printf(" xxx cval[NBT=%d] = '%s' \n", NBT, cval );

    if ( cval[0] == '#'  ) { goto DONE ; }
    if ( cval[0] == '%'  ) { goto DONE ; }
    if ( cval[0] == '!'  ) { goto DONE ; }
    if ( cval[0] == '\r' ) { goto DONE ; } // <CR>

    if ( NBT < MXTEXPOSE_SPECTROGRAPH ) 
      { sscanf(cval , "%le", &INPUTS_SPECTRO.TEXPOSE_LIST[NBT] ); }

    // make sure each TEXPOSE is increasing 
    if ( NBT>0 ) {
      T     = INPUTS_SPECTRO.TEXPOSE_LIST[NBT] ;
      TLAST = INPUTS_SPECTRO.TEXPOSE_LIST[NBT-1] ;
      if ( T < TLAST ) {
	sprintf(c1err,"TEXPOSE_LIST must be in increasing order.");
	sprintf(c2err,"TEXPOSE_LIST[%d,%d] = %.2f , %.2f ",
		NBT-1,NBT, TLAST,T);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }

    NBT++ ;     
    ptrtok = strtok(NULL, " ");
  }

 DONE:
  INPUTS_SPECTRO.NBIN_TEXPOSE = NBT ;
  INPUTS_SPECTRO.TEXPOSE_MIN  = INPUTS_SPECTRO.TEXPOSE_LIST[0] ;
  INPUTS_SPECTRO.TEXPOSE_MAX  = INPUTS_SPECTRO.TEXPOSE_LIST[NBT-1] ;

  if ( NBT >= MXTEXPOSE_SPECTROGRAPH ){ 
    sprintf(c1err,"Found %d TEXPOSE_LIST values", NBT);
    sprintf(c2err,"but MXTEXPOSE_SPECTROGRAPH=%d", MXTEXPOSE_SPECTROGRAPH);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  return(NBT);

} // end read_TEXPOSE_LIST



// ========================================================
void malloc_spectrograph(int OPT, int NBIN_LAM, int NBIN_TEXPOSE) {

  // OPT = +1 --> malloc
  // OPT = -1 --> free
  //
  // May 2021; allocate MXLAM_EXTEND more in wavelength

  int  MXLAM_EXTEND = MXLAM_SPECTROGRAPH_EXTEND;

  int  NBIN_LAM_LOCAL  = NBIN_LAM + MXLAM_EXTEND;
  int  MEML0 = NBIN_LAM_LOCAL     * sizeof(double);
  int  MEML1 = NBIN_LAM_LOCAL     * sizeof(double*);
  int  MEMT  = NBIN_TEXPOSE       * sizeof(double);
  int  MEMB  = NBIN_LAM_LOCAL     * sizeof(bool);
  int  ilam;
  //  char fnam[] = "malloc_spectrograph" ;

  // ------------ BEGIN ---------------

  if ( OPT > 0 ) {
    INPUTS_SPECTRO.LAMMIN_LIST   = (double*) malloc(MEML0);
    INPUTS_SPECTRO.LAMMAX_LIST   = (double*) malloc(MEML0);
    INPUTS_SPECTRO.LAMAVG_LIST   = (double*) malloc(MEML0);
    INPUTS_SPECTRO.LAMSIGMA_LIST = (double*) malloc(MEML0);
    INPUTS_SPECTRO.LAMBIN_LIST   = (double*) malloc(MEML0);

    INPUTS_SPECTRO.SNR0      = (double**) malloc(MEML1);
    INPUTS_SPECTRO.SNR1      = (double**) malloc(MEML1);
    INPUTS_SPECTRO.ZP        = (double**) malloc(MEML1);
    INPUTS_SPECTRO.SQSIGSKY  = (double**) malloc(MEML1);

    INPUTS_SPECTRO.ISLAM_EXTEND_LIST = (bool*) malloc(MEMB);

    for(ilam=0; ilam < NBIN_LAM_LOCAL; ilam++ ) {
      INPUTS_SPECTRO.SNR0[ilam]      = (double*) malloc(MEMT);
      INPUTS_SPECTRO.SNR1[ilam]      = (double*) malloc(MEMT);
      INPUTS_SPECTRO.ZP[ilam]        = (double*) malloc(MEMT);
      INPUTS_SPECTRO.SQSIGSKY[ilam]  = (double*) malloc(MEMT);

      INPUTS_SPECTRO.LAMMIN_LIST[ilam]   = -999.0 ;
      INPUTS_SPECTRO.LAMMAX_LIST[ilam]   = -99.0 ;
      INPUTS_SPECTRO.LAMAVG_LIST[ilam]   = -9.0 ;
      INPUTS_SPECTRO.LAMSIGMA_LIST[ilam] =  0.0 ;
      INPUTS_SPECTRO.LAMBIN_LIST[ilam]   =  0.0 ;
    }
  }
  else {
    free(INPUTS_SPECTRO.LAMMIN_LIST);
    free(INPUTS_SPECTRO.LAMMAX_LIST);
    free(INPUTS_SPECTRO.LAMAVG_LIST);
    free(INPUTS_SPECTRO.LAMSIGMA_LIST);
    free(INPUTS_SPECTRO.LAMBIN_LIST);
    free(INPUTS_SPECTRO.SNR0);
    free(INPUTS_SPECTRO.SNR1);
    free(INPUTS_SPECTRO.ZP);
    free(INPUTS_SPECTRO.SQSIGSKY);
    for(ilam=0; ilam < NBIN_LAM_LOCAL; ilam++ ) {
      free( INPUTS_SPECTRO.SNR0[ilam]  );
      free( INPUTS_SPECTRO.SNR1[ilam]  );
      free( INPUTS_SPECTRO.ZP[ilam] );
      free( INPUTS_SPECTRO.SQSIGSKY[ilam] );
    }
  }


  return ;

} // end malloc_spectrograph


// =================================
void  solve_spectrograph(void) {

  // in each lambda bin, solve for Zpe and SQSIGSKY
  //
  // July 8 2020: 
  //  + HACK SNR1 if m1-m0 is slightly less than 2.5log10(SNR0/SNR1)

  int l,t, iref, ITexpose ;
  int NBL = INPUTS_SPECTRO.NBIN_LAM ;
  int NBT = INPUTS_SPECTRO.NBIN_TEXPOSE ;

  double MAGREF[2], POWMAG[2], SQPOWMAG[2], ARG, SNR[2], SNR1_ORIG ;
  double TOP, BOT, ZP, SQSIGSKY, F[2], DUM0, DUM1, LAMMIN, LAMMAX, LAMAVG ;
  double SNR_check[2], check[2], MAGREF_DIF, MAGSNR_DIF, magCheck ;
  int    LDMP_SNRFIX = 1;
  char  msg[100];
  char fnam[] = "solve_spectrograph" ;

  // ------------- BEGIN ---------------

  sprintf(msg,"%s: solve for ZP and SQSIG for each LAMBDA and Texpose", fnam);
  print_banner(msg);

  for(iref=0; iref < 2; iref++ ) {
    MAGREF[iref] = INPUTS_SPECTRO.MAGREF_LIST[iref];
    ARG  = -0.4 * MAGREF[iref] ; 
    POWMAG[iref] = pow(TEN,ARG); 
    SQPOWMAG[iref] = POWMAG[iref] * POWMAG[iref] ;
  }

  
  for(l=0; l < NBL; l++ ) {

    LAMMIN = INPUTS_SPECTRO.LAMMIN_LIST[l] ;
    LAMMAX = INPUTS_SPECTRO.LAMMAX_LIST[l] ;
    LAMAVG = 0.5 * ( LAMMIN + LAMMAX );

    for(t=0; t < NBT; t++ ) {

      ITexpose = (int)INPUTS_SPECTRO.TEXPOSE_LIST[t]; // for error msg only
      SNR[0] = INPUTS_SPECTRO.SNR0[l][t] ;
      SNR[1] = INPUTS_SPECTRO.SNR1[l][t] ;
      TOP    = POWMAG[0] - POWMAG[1] ;

      // sanity check:
      if ( SNR[0] <= 1.0E-9 || SNR[1] < 1.0E-9 ) {
        sprintf(c1err,"Invalid SNR[0,1] = %le, %le", SNR[0], SNR[1]);
	sprintf(c2err,"Check <LAM> = %.1f, Texpose=%d", LAMAVG, ITexpose);
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }
      
      MAGREF_DIF = MAGREF[1] - MAGREF[0] ;
      MAGSNR_DIF = 2.5*log10(SNR[0]/SNR[1]);
      magCheck   = (MAGREF_DIF - MAGSNR_DIF);
      if ( magCheck < INPUTS_SPECTRO.MAGSNR_TOLERANCE_ABORT ) {
        sprintf(c1err,"failed solution check at LAM=%.1f, Texpose=%d",
		LAMAVG, ITexpose );
	sprintf(c2err,"mref1-mref0 = %f < 2.5log10(SNR0/SNR1) = %f",
		MAGREF_DIF, MAGSNR_DIF);
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }

      if ( magCheck < 0.0 ) {  
	// allow a little slop and increase SNR[1] a tiny bit
	// so that a ZP,SIG solution is possible
	SNR1_ORIG = SNR[1];
	ARG    = -0.4 * MAGREF_DIF ;
	SNR[1] = 1.005 * pow(TEN,ARG) * SNR[0] ;
	if ( LDMP_SNRFIX ) {
	  printf(" HACK-SNR1: LAM=%.1f  SNR0=%f  SNR1=%f -> %f \n",
		 fnam, LAMAVG, SNR[0], SNR1_ORIG, SNR[1] );
	  fflush(stdout);
	}
      }

      DUM0 = POWMAG[0]/SNR[0];
      DUM1 = POWMAG[1]/SNR[1];
      BOT  = DUM0*DUM0 - DUM1*DUM1 ;
	
      if ( TOP <= 0.0 || BOT <= 0.0 ) {
	print_preAbort_banner(fnam);
	printf("\t BOT = %le  and  TOP = %le\n", BOT, TOP);
	printf("\t BOT = (%le)^2 - (%le)^2 \n", DUM0, DUM1);
	printf("\t TOP = %le - %le \n", POWMAG[0], POWMAG[1]);
	printf("\t SNR[0]=%le  SNR[1]=%le \n", SNR[0], SNR[1] );
	sprintf(c1err,"Cannot solve ZP for LAM=%.1f to %.1f,  and t=%d sec",
		LAMMIN, LAMMAX, ITexpose );
	sprintf(c2err,"Check SPECTROGRAPH") ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }
      
      ZP     = 2.5*log10(TOP/BOT) ;  // photo-electrons
      F[0]   = pow(TEN, -0.4*(MAGREF[0]-ZP) );
      F[1]   = pow(TEN, -0.4*(MAGREF[1]-ZP) );
      SQSIGSKY = (F[0]/SNR[0])*(F[0]/SNR[0]) - F[0] ;
      
      // store in global 
      INPUTS_SPECTRO.ZP[l][t]        = ZP ;
      INPUTS_SPECTRO.SQSIGSKY[l][t]  = SQSIGSKY ;

      // make sure original SNR can be computed from ZP and SQSIGSKY
      for(iref=0; iref<2; iref++ ) {  
	SNR_check[iref] = F[iref] / sqrt( SQSIGSKY + F[iref] ) ;
	check[iref] = fabs(SNR[iref]/SNR_check[iref]-1.0) ;
      }

      if ( check[0] > 0.001  ||  check[1] > 0.001 ) {
	print_preAbort_banner(fnam);
	printf("   SNR0(input/check) = %f/%f = %f \n",
	       SNR[0], SNR_check[0], SNR[0]/SNR_check[0] );
	printf("   SNR1(input/check) = %f/%f = %f \n",
	       SNR[1], SNR_check[1], SNR[1]/SNR_check[1] );
	printf("   Lambda bin: %f to %f \n", LAMMIN, LAMMAX);
	printf("   F0=%le   F1=%le  \n", F[0], F[1] );
	sprintf(c1err,"Problem computing ZP and SQSIGSKY" );
	sprintf(c2err,"ZP=%f  SQSKYSIG=%f", ZP, SQSIGSKY );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
      }

    }  // end t bins
  }  // end l bins

  return ;

} // end solve_spectrograph



// ================================================
void get_FILTERtrans_spectrograph(double *LAMFILT_MIN, double *LAMFILT_MAX, 
				  int MXTRANS, int *NBIN_LAMFILT,
				  double *LAMFILT_ARRAY, 
				  double *TRANSFILT_ARRAY ) 
{  
  // Return IFU transmission efficiency (TRANS_ARRAY) vs. 
  // wavelength (LAM_ARRAY), in the wavelength region LMIN to LMAX.
  // If number of wavelength bins (NBIN_TRANS) exceeds MXTRANS --> ABORT.
  // 
  // In each wavelength bin, Eff ~ Flux ~ 10^[-0.4*(MAGREF-ZP)
  //
  // and normlize max trans = 1.0.
  //
  // If there are multuple T_EXPOSE, then use first one since
  // relative trans should be the same regardless of T_EXPOSE.
  //
  // Inputs:
  //  LAMFILT_MIN[MAX] : min,max wavelength of filter slice
  //  MXTRANS    : max size of output array; abort on overflow
  //
  // Outputs:
  //   LAMFILT_MIN[LMAX]  : extended lambda range so that Trans->0 at edges.
  //   NBIN_TRANS  : number of filter-transmission bins
  //   LAM_ARRAY   : wavelength array
  //   TRANS_ARRAY : transmission array
  //
  // Oct 27 2016: protect interpolation in 1st half of 1st bin
  //              and last half of last bin
  //
  // Feb 01 2017: 
  //  + fix malloc and bugs related to mixing up NBL_TRANS [NB]
  //    and NBL_SPECTRO.
  //    
  int l, t, MEMD ;
  double *FLUX_TMP, *ZP_TMP, LCEN, F, ZP, ARG, FLUX_MAX ;
  double LAMRANGE, LAMSTEP, LAMSTEP_APPROX, xlam ;
  int NBL_SPECTRO = INPUTS_SPECTRO.NBIN_LAM ;
  int NBL_TRANS ;
  int OPT_INTERP=1;

  char fnam[] = "get_FILTERtrans_spectrograph" ;

  // -------------- BEGIN -------------

  // define finer lambda bins for intergration, with a minimum of 10 bins.
  LAMSTEP_APPROX = 5.0 ;
  LAMRANGE       = *LAMFILT_MAX - *LAMFILT_MIN ;
  NBL_TRANS      = (int)(LAMRANGE/LAMSTEP_APPROX);   
  if(NBL_TRANS<10) {NBL_TRANS=10;}
  LAMSTEP        = LAMRANGE/(double)NBL_TRANS;

  // add one bin on the low and high side where transmission=0
  *LAMFILT_MIN -= LAMSTEP ;  NBL_TRANS++ ;
  *LAMFILT_MAX += LAMSTEP ;  NBL_TRANS++ ;

  if ( NBL_TRANS >= MXTRANS ) {
    sprintf(c1err,"NBL_TRANS=%d exceeds bound MXTRANS=%d", 
	    NBL_TRANS, MXTRANS);
    sprintf(c2err,"filter lambda range %.1f to %.1f, lamstep=%.3f", 
	    *LAMFILT_MIN, *LAMFILT_MAX, LAMSTEP ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // allocate local memory
  MEMD     = sizeof(double) ;
  FLUX_TMP = (double*) malloc ( MEMD * NBL_TRANS   ) ;
  ZP_TMP   = (double*) malloc ( MEMD * NBL_SPECTRO ) ;
  FLUX_MAX = 0.0 ;

  // load local ZP_TMP[ilam] array to use for interpolation
  t=0; // always use first T_EXPOSE bin
  for(l=0; l < NBL_SPECTRO ; l++ ) { ZP_TMP[l] = INPUTS_SPECTRO.ZP[l][t];  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for(l=0; l < NBL_TRANS; l++ ) {

    FLUX_TMP[l] = 0.0 ;

    xlam = (double)l ;
    LCEN = *LAMFILT_MIN + LAMSTEP*(xlam+0.5); // center of bin
    F    = 0.0 ;

    // interpolate ZP
    if ( l==0 || l == NBL_TRANS-1 ) { 
      F  = 0.0 ; // first & last bin has FLUX_TRANS=0
      ZP = 0.0 ;
    }
    else if ( LCEN <= INPUTS_SPECTRO.LAMAVG_LIST[0] ) {
      ZP = ZP_TMP[0];
    }
    else if ( LCEN >= INPUTS_SPECTRO.LAMAVG_LIST[NBL_SPECTRO-1] ) {
      ZP = ZP_TMP[NBL_SPECTRO-1];
    }
    else {
      ZP = interp_1DFUN (OPT_INTERP, LCEN, NBL_SPECTRO, 
			 INPUTS_SPECTRO.LAMAVG_LIST, ZP_TMP, fnam );      
    }

    if ( ZP > 0.001 ) {
      ARG = -0.4*(INPUTS_SPECTRO.MAGREF_LIST[0] - ZP);
      F   = pow(TEN,ARG) ;
      if ( F > FLUX_MAX ) { FLUX_MAX = F ; }
    }

    FLUX_TMP[l]      = F ; 
    LAMFILT_ARRAY[l] = LCEN ;

  } // end lambda loop
  
  // --------------------
  if ( FLUX_MAX < 1.0E-9 ) {
    sprintf(c1err,"Synthetic FLUX_MAX=%f ", FLUX_MAX);	   
    sprintf(c2err,"LAMFILT_MIN/MAX=%.1f/%.2f  NBL_TRANS=%d", 
	    *LAMFILT_MIN, *LAMFILT_MAX, NBL_TRANS) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  *NBIN_LAMFILT = NBL_TRANS ;
  for(l=0; l < NBL_TRANS ; l++ ) 
    { TRANSFILT_ARRAY[l] =  FLUX_TMP[l] / FLUX_MAX ; }


  free(FLUX_TMP);  free(ZP_TMP);
  return ;

} // end get_FILTERtrans_spectrograph


// =====================================
void read_spectrograph_fits(char *inFile) {

  // Read spectrograph info from FITS-formatted *inFile
  // that was created by kcor.exe.
  // The read-back includes each wavelength bin range,
  // and its ZP and SQSIGSKY. The original SNR1,SNR2 
  // are not stored.
  //
  // if SPECTROGRAPH key is in header, read corresponding table.
  // if SPECTROGRAPH key is not in header, return with no message.
  //
  // Oct 14 2016: read LAMSIGMA_LIST
  // Sep 19 2018: fill INPUTS_SPECTRO.ISFIX_LAMBIN (used for output FORMAT)
  // May 06 2020: default format is to format LAMMIN & LAMMAX 
  //  instead of just LAMCEN
  // Aug 20 2021: store LOGTEXPOSE_LIST

  int istat, hdutype, extver, icol, anynul ;
  fitsfile *fp ;

  float  tmpVal_f  ;
  int    NBL, NBT, l, t ;
  double L0, L1;

  char keyName[40], comment[80], TBLname[40], INFILE[MXPATHLEN] ;
  char fnam[] = "read_spectrograph_fits" ;

  // --------------- BEGIN -----------------

  SPECTROGRAPH_USEFLAG = 0;

  // open fits file
  istat = 0 ;
  sprintf(INFILE, "%s", inFile);
  fits_open_file(&fp, INFILE, READONLY, &istat );
  
  // if inFile does not exist, try official area
  if ( istat != 0 ) {
    sprintf(INFILE, "%s/kcor/%s", PATH_SNDATA_ROOT, inFile);
    istat = 0 ;
    fits_open_file(&fp, INFILE, READONLY, &istat );
  }

  sprintf(c1err,"Open %s", INFILE );
  snfitsio_errorCheck(c1err, istat);

  // -------------------------------------------------
  sprintf(TBLname, "%s", FITSTABLE_NAME_SPECTROGRAPH );

  // -------------------------------------------------
  // check for SPECTROGRAPH keys in header.
  sprintf(comment,"Read spectrograph from kcor file");

  sprintf(keyName, "%s", "SPECTROGRAPH_INSTRUMENT" );
  fits_read_key(fp, TSTRING, keyName, &INPUTS_SPECTRO.INSTRUMENT_NAME, 
		comment, &istat );
  if( istat != 0 ) { goto FCLOSE ; }


  sprintf(keyName, "%s", "SPECTROGRAPH_FILTERLIST" );
  fits_read_key(fp, TSTRING, keyName, &INPUTS_SPECTRO.SYN_FILTERLIST_BAND, 
		comment, &istat );

  printf("\n Read spectrograph instrument '%s' \n", 
	 INPUTS_SPECTRO.INSTRUMENT_NAME );
  fflush(stdout);

  SPECTROGRAPH_USEFLAG = 1 ; // set global flag that spectrograph is defined.

  extver = istat = 0 ;  hdutype = BINARY_TBL ;
  fits_movnam_hdu( fp, hdutype, TBLname, extver, &istat);
  sprintf(c1err,"movnam to %s table (hdutype=%d)",  TBLname, hdutype ) ;
  snfitsio_errorCheck(c1err, istat);

  // - - - - - -
  // read binning
  sprintf(keyName, "%s", "NBL" );
  fits_read_key(fp, TINT, keyName, &NBL, comment, &istat );
  sprintf(c1err,"read number of lambda bins");
  snfitsio_errorCheck(c1err, istat);
  printf("   Found %d wavelength bins \n", NBL);
  
  sprintf(keyName, "%s", "NBT" );
  fits_read_key(fp, TINT, keyName, &NBT, comment, &istat );
  sprintf(c1err,"read number of TEXPOSE bins");
  snfitsio_errorCheck(c1err, istat);
  printf("   Found %d TEXPOSE bins \n", NBT );

  fflush(stdout);
  INPUTS_SPECTRO.NBIN_LAM     = NBL ;
  INPUTS_SPECTRO.NBIN_TEXPOSE = NBT ;
  
  // now read each exposure time:
  printf("\t TEXPOSE(seconds) = ");
  for(t=0; t < NBT; t++ ) {
    sprintf(keyName, "TEXPOSE%2.2d", t );
    fits_read_key(fp, TFLOAT, keyName, &tmpVal_f, comment, &istat );
    sprintf(c1err,"read %s", keyName);
    snfitsio_errorCheck(c1err, istat);

    printf("%d ", (int)tmpVal_f );
    INPUTS_SPECTRO.TEXPOSE_LIST[t] = tmpVal_f ;
    INPUTS_SPECTRO.LOGTEXPOSE_LIST[t] = log10(tmpVal_f) ; // Aug 20 2021
  }
  printf("\n"); fflush(stdout);

  // set global min & max
  INPUTS_SPECTRO.TEXPOSE_MIN = INPUTS_SPECTRO.TEXPOSE_LIST[0] ;
  INPUTS_SPECTRO.TEXPOSE_MAX = INPUTS_SPECTRO.TEXPOSE_LIST[NBT-1] ;

  // malloc arrays before reading
  malloc_spectrograph(1,NBL,NBT);

  // ---------------------------
  long NROW, FIRSTELEM, FIRSTROW ;
  FIRSTELEM = FIRSTROW = 1 ;  NROW=NBL ;  anynul=istat=0 ;
  
  // read lambda range for each wavelength bin
  icol = 1 ;
  fits_read_col_dbl(fp, icol, FIRSTROW, FIRSTELEM, NROW, NULL_1D, 
		    INPUTS_SPECTRO.LAMMIN_LIST, &anynul, &istat );
  sprintf(c1err,"read LAMMIN_LIST column" );
  snfitsio_errorCheck(c1err, istat);

  icol = 2 ;
  fits_read_col_dbl(fp, icol, FIRSTROW, FIRSTELEM, NROW, NULL_1D, 
		    INPUTS_SPECTRO.LAMMAX_LIST, &anynul, &istat );
  sprintf(c1err,"read LAMMAX_LIST column" );
  snfitsio_errorCheck(c1err, istat);

  printf("   Wavelength range stored: %.2f to %.2f A \n",
	 INPUTS_SPECTRO.LAMMIN_LIST[0], INPUTS_SPECTRO.LAMMAX_LIST[NBL-1]);

  icol = 3 ;
  fits_read_col_dbl(fp, icol, FIRSTROW, FIRSTELEM, NROW, NULL_1D, 
		    INPUTS_SPECTRO.LAMSIGMA_LIST, &anynul, &istat );
  sprintf(c1err,"read LAMSIGMA_LIST column" );
  snfitsio_errorCheck(c1err, istat);

  // compute LAMAVG & LAMBIN
  double LBIN, LASTBIN=0.0 ;
  INPUTS_SPECTRO.FORMAT_MASK = 2; // default: write LAMMIN & LAMMAX
  for(l=0; l <NBL; l++ ) {
    L0 = INPUTS_SPECTRO.LAMMIN_LIST[l] ;
    L1 = INPUTS_SPECTRO.LAMMAX_LIST[l] ;
    LBIN = L1-L0 ;
    INPUTS_SPECTRO.LAMAVG_LIST[l] = ( L0 + L1 ) / 2.0 ;
    INPUTS_SPECTRO.LAMBIN_LIST[l] = LBIN;

    INPUTS_SPECTRO.ISLAM_EXTEND_LIST[l] = false ;

    LASTBIN=LBIN; 
  }


  // set global min & max
  INPUTS_SPECTRO.LAM_MIN = INPUTS_SPECTRO.LAMMIN_LIST[0] ;
  INPUTS_SPECTRO.LAM_MAX = INPUTS_SPECTRO.LAMMAX_LIST[NBL-1] ;

  // read ZP and SQSIGSKY for each Texpose.
  // They are stored in fits file as float, but array is double,
  // so use ZP_f and SQ_f as intermediate array.

  int NBL_MALLOC = NBL + MXLAM_SPECTROGRAPH_EXTEND ;
  float *ZP_f = (float*)malloc( NBL_MALLOC * sizeof(float) ) ;
  float *SQ_f = (float*)malloc( NBL_MALLOC * sizeof(float) ) ;
  double  LAMMIN_ZP=1.0E9, LAMMAX_ZP=0.0, LAM, ZP, SQ ;

  for(t=0; t < NBT; t++ ) {
    icol++ ;
    fits_read_col_flt(fp, icol, FIRSTROW, FIRSTELEM, NROW, NULL_1E, 
		      ZP_f, &anynul, &istat );
    sprintf(c1err,"read ZP  column" );
    snfitsio_errorCheck(c1err, istat);

    icol++ ;
    fits_read_col_flt(fp, icol, FIRSTROW, FIRSTELEM, NROW, NULL_1E, 
		      SQ_f, &anynul, &istat );
    sprintf(c1err,"read SQSIGSKY column" );
    snfitsio_errorCheck(c1err, istat);

    for(l=0; l <NBL; l++ ) {
      ZP = (double)ZP_f[l] ;
      SQ = (double)SQ_f[l] ;
      LAM = INPUTS_SPECTRO.LAMAVG_LIST[l];
      INPUTS_SPECTRO.ZP[l][t]       = ZP ;
      INPUTS_SPECTRO.SQSIGSKY[l][t] = SQ ;

      if ( ZP != 0.0 ) {
	if ( LAMMIN_ZP > 0.9E9 ) { LAMMIN_ZP = LAM; }
	LAMMAX_ZP = LAM;
      }

    } // end l loop over lambda
  } // end t loop over Texpose

  /* xxx
  printf("   Wavelength range with valid ZP: %.1f to %.1f \n",
  LAMMIN_ZP, LAMMAX_ZP); xxx */

  free(ZP_f);  free(SQ_f);

  // ---------------------------------------------------
  // ---------------------------------------------------
  // read 2nd table of lambda range vs. SYN_FILTER_SPECTROGRAPH
  // ---------------------------------------------------
  // ---------------------------------------------------

  float LAMMIN_f[MXFILTINDX], LAMMAX_f[MXFILTINDX];
  int ifilt ;
  char *cName[MXFILTINDX] ;

  sprintf(TBLname, "SYN_FILTER_SPECTROGRAPH" );

  extver = istat = 0 ;  hdutype = BINARY_TBL ;
  fits_movnam_hdu( fp, hdutype, TBLname, extver, &istat);
  sprintf(c1err,"movnam to %s table (hdutype=%d)",  TBLname, hdutype ) ;
  snfitsio_errorCheck(c1err, istat);

  FIRSTELEM = FIRSTROW = 1 ;  anynul=istat=0 ;
  NROW = strlen(INPUTS_SPECTRO.SYN_FILTERLIST_BAND);
  for(ifilt=0; ifilt<NROW; ifilt++ )  { 
    cName[ifilt] = INPUTS_SPECTRO.SYN_FILTERLIST_NAME[ifilt] ; 
    cName[ifilt][0] = 0 ;
    INPUTS_SPECTRO.SYN_FILTERLIST_LAMMIN[ifilt] = -9.0 ;
    INPUTS_SPECTRO.SYN_FILTERLIST_LAMMAX[ifilt] = -9.0 ;
  }


  icol = 1 ;
  fits_read_col_str(fp, icol, FIRSTROW, FIRSTELEM, NROW, NULL_A, 
		    &cName[0], &anynul, &istat );
  sprintf(c1err,"read SYN_FILTERLIST_NAME column" );
  snfitsio_errorCheck(c1err, istat);

  icol = 2 ;
  fits_read_col_flt(fp, icol, FIRSTROW, FIRSTELEM, NROW, NULL_1E, 
		    LAMMIN_f, &anynul, &istat );
  sprintf(c1err,"read SYN_FILTERLIST_LAMMIN  column " );
  snfitsio_errorCheck(c1err, istat);

  icol = 3 ;
  fits_read_col_flt(fp, icol, FIRSTROW, FIRSTELEM, NROW, NULL_1E, 
		    LAMMAX_f, &anynul, &istat );
  sprintf(c1err,"read SYN_FILTERLIST_LAMMAX  column" );
  snfitsio_errorCheck(c1err, istat);

  
  // note this is a sparse "ifilt" over SYN_FILTERLIST,
  // and not over all kcor filters.
  for(ifilt=0 ; ifilt < NROW; ifilt++ ) {
    INPUTS_SPECTRO.SYN_FILTERLIST_LAMMIN[ifilt] = LAMMIN_f[ifilt] ;
    INPUTS_SPECTRO.SYN_FILTERLIST_LAMMAX[ifilt] = LAMMAX_f[ifilt] ;

    /*
    printf(" xxx '%s' : LAMRANGE = %.1f to %.1f \n"
	   ,INPUTS_SPECTRO.SYN_FILTERLIST_NAME[ifilt]
	   ,INPUTS_SPECTRO.SYN_FILTERLIST_LAMMIN[ifilt]
	   ,INPUTS_SPECTRO.SYN_FILTERLIST_LAMMAX[ifilt] );  */
  }

  // ------------------------------------------
  // close fits file
 FCLOSE:
  istat = 0 ;
  fits_close_file(fp, &istat);
  sprintf(c1err, "Close Spectrograph FITS file"  );
  snfitsio_errorCheck(c1err, istat);


  return ;

} // end read_spectrograph_fits


// ====================================================
void extend_spectrograph_lambins(void) {

  // May 2021
  // extend spectrograph wavelength bins blue and red
  // to handle line smearing. The extra bins are used only
  // to smear flux into the nominal bins; the extra bins
  // are never recorded in the output spectra

  int    MXLAM_EXTEND = MXLAM_SPECTROGRAPH_EXTEND ;
  int    NLAM_ORIG    = INPUTS_SPECTRO.NBIN_LAM;
  int t, NBT          = INPUTS_SPECTRO.NBIN_TEXPOSE ;  
  double NSIG_EXTEND  = 2.5 ;

  double *LAMMIN_LIST = INPUTS_SPECTRO.LAMMIN_LIST ;
  double *LAMMAX_LIST = INPUTS_SPECTRO.LAMMAX_LIST ;
  double *LAMAVG_LIST = INPUTS_SPECTRO.LAMAVG_LIST ;
  double *LAMBIN_LIST = INPUTS_SPECTRO.LAMBIN_LIST ;
  double *LAMSIG_LIST = INPUTS_SPECTRO.LAMSIGMA_LIST ;
  double **ZP_LIST    = INPUTS_SPECTRO.ZP ;
  double **SQSIG_LIST = INPUTS_SPECTRO.SQSIGSKY ;
  //  double **SNR0_LIST  = INPUTS_SPECTRO.SNR0 ;
  //  double **SNR1_LIST  = INPUTS_SPECTRO.SNR1 ;

  int    firstlast, ilam, NLAM_TMP ;
  int    NLAM_EXTEND[2], NLAM_EXTEND_TOT=0 ; 
  double LAMSIG, LAMMIN, LAMMAX, LAMBIN ;
  bool   ISLAM_EXTEND ;
  char TEXT_FIRSTLAST[2][8] = { "FIRST", "LAST" } ;
  char TEXT_COLOR[2][8]     = { "BLUE",  "RED " } ;
  char fnam[] = "extend_spectrograph_lambins" ;

  // ---------- BEGIN ----------


  for ( firstlast=0; firstlast <2; firstlast++ ) {
    
    if ( firstlast == 0 ) {
      ilam = 0 ;
    }
    else {
      ilam = NLAM_ORIG-2; // avoid last bin in case it's too small
    }

    LAMMIN = LAMMIN_LIST[ilam];
    LAMMAX = LAMMAX_LIST[ilam];
    LAMBIN = LAMMAX - LAMMIN ;
    LAMSIG = LAMSIG_LIST[ilam];
    NLAM_TMP = (int)((LAMSIG*NSIG_EXTEND)/LAMBIN); 
    NLAM_TMP++ ;

    NLAM_EXTEND_TOT += NLAM_TMP;
    NLAM_EXTEND[firstlast] = NLAM_TMP ;

    printf("\t Extend SPECTROGRAPH %s-edge by %d bins for wave-resolution\n",
	   TEXT_COLOR[firstlast], NLAM_TMP );
  }
  
  int NLAM_NEW = NLAM_ORIG + NLAM_EXTEND_TOT ;
  printf("\t Extend NBLAM = %d -> %d \n", NLAM_ORIG, NLAM_NEW);

  if ( NLAM_EXTEND_TOT >= MXLAM_EXTEND ) {
    sprintf(c1err,"NLAM_EXTEND = %d exceeds MXLAM_SPECTROGRAPH_EXTEND=%d",
	    NLAM_EXTEND_TOT, MXLAM_EXTEND);
    sprintf(c2err,"Check lam-resolution, or inrease bound") ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  
  // re-write arrays in reverse order so that arrays can be
  // over-written without creating a copy
  int ilam_copy, ilam_new ;
  int NLAM_BLUE = NLAM_EXTEND[0];  // number of extended blue bins
  int ILAM_MARK[3] = { 0, NLAM_BLUE, NLAM_ORIG+NLAM_BLUE } ;
  double xtmp;

  //  printf(" xxx %s: NLAM = %d -> %d\n", fnam, NLAM_ORIG, NLAM_NEW);

  for ( ilam_new=NLAM_NEW-1; ilam_new >=0; ilam_new-- ) {

    ISLAM_EXTEND = false;

    if ( ilam_new >= ILAM_MARK[2] ) {
      // extended lambda bins on RED side
      ISLAM_EXTEND = true ;
      ilam_copy = NLAM_ORIG - 2;
      xtmp   = (double)(ilam_new - ILAM_MARK[2] + 1);
      LAMBIN = LAMBIN_LIST[ilam_copy];
      LAMMIN = LAMMIN_LIST[NLAM_ORIG-1] + xtmp*LAMBIN ;
      LAMMAX = LAMMIN + LAMBIN ;
    }
    else if ( ilam_new >= ILAM_MARK[1] ) {
      // original lambda bins
      ilam_copy = ilam_new - NLAM_BLUE ;
      LAMMIN = LAMMIN_LIST[ilam_copy];
      LAMMAX = LAMMAX_LIST[ilam_copy];
    }
    else {
      // extended lambda bins on BLUE side
      ISLAM_EXTEND = true ;
      ilam_copy = NLAM_BLUE;
      xtmp   = (double)(NLAM_BLUE-ilam_new);
      LAMBIN = LAMBIN_LIST[ilam_copy];
      LAMMIN = LAMMIN_LIST[ilam_copy] - xtmp*LAMBIN ;
      LAMMAX = LAMMIN + LAMBIN ;

      //   printf(" xxx %s: ilam_new/cp=%d/%d  xtmp=%.1f LAMMIN=%.1f \n",
      //   fnam, ilam_new, ilam_copy, xtmp, LAMMIN);
    }

    LAMMIN_LIST[ilam_new] = LAMMIN ;
    LAMMAX_LIST[ilam_new] = LAMMAX ;
    LAMAVG_LIST[ilam_new] = (LAMMIN + LAMMAX)/2.0;

    LAMSIG_LIST[ilam_new]  =  LAMSIG_LIST[ilam_copy];
    LAMBIN_LIST[ilam_new]  =  LAMBIN_LIST[ilam_copy];

    INPUTS_SPECTRO.ISLAM_EXTEND_LIST[ilam_new] = ISLAM_EXTEND;

    for ( t=0; t < NBT; t++ ) {
      ZP_LIST[ilam_new][t]    = ZP_LIST[ilam_copy][t];
      SQSIG_LIST[ilam_new][t] = SQSIG_LIST[ilam_copy][t];
    }


  } // end ilam

  // - - - -
  INPUTS_SPECTRO.NBIN_LAM = NLAM_NEW ;

  dump_INPUTS_SPECTRO(4,"+extended wave bins");

  return ;

} // end extend_spectrograph_lambins


// ====================================================
void dump_INPUTS_SPECTRO(int nbin_dump, char *comment) {

  // Created May 2021
  // dump spectrograph table for first and last "nbin_dump" lambda bins.

  int l, NBL = INPUTS_SPECTRO.NBIN_LAM ;  
  int t, NBT = INPUTS_SPECTRO.NBIN_TEXPOSE ;  
  double LAMAVG, LAMBIN, LAMSIG, ZP, SIGSKY, VARSKY, SNR0, SNR1;
  bool   ISLAM_EXTEND ;

  // ---------- BEGIN -----------

  t=0;
  printf("\n DUMP SPECTROGRAPH TABLE: %s\n", comment);
  printf("  lamBin   LAMAVG  LAMBIN  LAMSIG   ZP[%d] SIGSKY[%d] "
	 "Extended\n", t, t);

  for(l=0; l < NBL; l++ ) {
    if ( l >= nbin_dump && l < NBL-nbin_dump ) { continue; }

    LAMAVG = INPUTS_SPECTRO.LAMAVG_LIST[l];
    LAMBIN = INPUTS_SPECTRO.LAMBIN_LIST[l];
    LAMSIG = INPUTS_SPECTRO.LAMSIGMA_LIST[l];
    ISLAM_EXTEND = INPUTS_SPECTRO.ISLAM_EXTEND_LIST[l] ;

    ZP     = INPUTS_SPECTRO.ZP[l][t] ;
    VARSKY = INPUTS_SPECTRO.SQSIGSKY[l][t] ;

    if ( VARSKY > 0.0 ) { SIGSKY = sqrt(VARSKY); }
    else                { SIGSKY = VARSKY; }

    printf(" %6d  %9.2f  %6.2f  %4.1f  " 
	   "%6.1f   %7.2f    %d\n",
	   l, LAMAVG, LAMBIN, LAMSIG, 
	   ZP, SIGSKY, ISLAM_EXTEND );
  }

  fflush(stdout);

  return;

} // end dump_INPUTS_SPECTRO

// ====================================================
double getSNR_spectrograph(int ILAM, double TEXPOSE_S, double TEXPOSE_T,
			   bool ALLOW_TEXTRAP, double GENMAG, double *ERRFRAC_T) {

  // Return SNR for inputs
  //  + SPECTROGRAPH wavelength bin (ILAM)
  //  + search exposure time (TEXPOSE_S)
  //  + template exposure time (TEXPOSE_T)
  //  + ALLOW_TEXTRAP=T -> allow extrapolating SNR ~ sqrt(TEXPOSE); else abort
  //  + magnitude in wavelength bin (GENMAG)
  //
  // *ERRFRAC_T  = sigma_template/FluxErrTot
  //  is the fraction of error associated with template noise.
  //  Used externally to generate correlated template noise.
  //
  // If Texpose is outside valid range, abort.
  // SQSIGSKY is returned as well.
  //
  // Feb  2 2017: fix awful bug and scale template noise to search-zp
  // May  6 2020: return SNR=0 if ZP=-9 (undefined)
  // May 22 2020: return SNR=0 if variance < 0 (see SQ_SUM)
  // May 27 2020: pass & implement new option ALLOW_TEXTRAP 
  // Jun 04 2020: init *ERRFRAC_T=0
  // Aug 20 2021: minor refac to interpolate ZP vs. log(Texpose)
  //               intead of ZP vs. Texpose

  int OPT_INTERP=1;
  int NBT       = INPUTS_SPECTRO.NBIN_TEXPOSE ;
  double Tmin   = INPUTS_SPECTRO.TEXPOSE_LIST[0] ;
  double Tmax   = INPUTS_SPECTRO.TEXPOSE_LIST[NBT-1] ;
  double TEXPOSE_S_local = TEXPOSE_S ;
  double LOGTEXPOSE_S, LOGTEXPOSE_T;
  //  double TEXPOSE_T_local = TEXPOSE_T ;
  double SNR, ZP_S, ZP_T, arg, SQ_S, SQ_T, SQ_SUM, Flux, FluxErr ;
  bool   DO_TEXTRAP = false;
  char fnam[] = "getSNR_spectrograph" ;
  char errmsg_ZP_S[] = "getSNR_spectrograph(ZP_S)";
  char errmsg_ZP_T[] = "getSNR_spectrograph(ZP_T)";
  char errmsg_SQ_S[] = "getSNR_spectrograph(SQ_S)";
  char errmsg_SQ_T[] = "getSNR_spectrograph(SQ_T)";
  int  LDMP = (ILAM < -3);
  // int  LDMP = ( fabs(INPUTS_SPECTRO.LAMAVG_LIST[ILAM]-8000.) < 2.0 );
  bool REFAC_ZP = true ;

  // -------------- BEGIN --------------

  SNR = SQ_S = SQ_T = ZP_S = ZP_T = 0.0 ;
  *ERRFRAC_T = 0.0 ;

  if ( ALLOW_TEXTRAP ) {
    if ( TEXPOSE_S < Tmin ) 
      { TEXPOSE_S_local = Tmin + 0.00001 ; DO_TEXTRAP = true;}
    if ( TEXPOSE_S > Tmax ) 
      { TEXPOSE_S_local = Tmax - 0.00001 ; DO_TEXTRAP = true ; }
  }
  else if ( TEXPOSE_S < Tmin  || TEXPOSE_S > Tmax ) {
    sprintf(c1err,"Invalid TEXPOSE_S = %f", TEXPOSE_S );
    sprintf(c2err,"Valid TEXPOSE_S range: %.2f to %.2f \n", Tmin, Tmax);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  LOGTEXPOSE_S = log10(TEXPOSE_S_local); // Aug 20 2021

  // May 2020: if ZP is undefined in this ILAM bin, return SNR=0
  if ( INPUTS_SPECTRO.ZP[ILAM][0] < 0.0 ) { return(SNR); }

  // interpolate ZP(Texpose) and SQSIG(Texpose)

  if ( REFAC_ZP ) {
    ZP_S = interp_1DFUN (OPT_INTERP, LOGTEXPOSE_S, NBT, 
			 INPUTS_SPECTRO.LOGTEXPOSE_LIST,
			 INPUTS_SPECTRO.ZP[ILAM], errmsg_ZP_S );
  }
  else {
    // legacy 
    ZP_S = interp_1DFUN (OPT_INTERP, TEXPOSE_S_local, NBT, 
			 INPUTS_SPECTRO.TEXPOSE_LIST,
			 INPUTS_SPECTRO.ZP[ILAM], errmsg_ZP_S );
  }

  SQ_S = interp_1DFUN (OPT_INTERP, TEXPOSE_S_local, NBT, 
		       INPUTS_SPECTRO.TEXPOSE_LIST,
		       INPUTS_SPECTRO.SQSIGSKY[ILAM], errmsg_SQ_S );
  
  if ( TEXPOSE_T > 0.01 ) {
    LOGTEXPOSE_T = log10(TEXPOSE_T);
    if ( REFAC_ZP ) {
      ZP_T = interp_1DFUN (OPT_INTERP, LOGTEXPOSE_T, NBT, 
			   INPUTS_SPECTRO.LOGTEXPOSE_LIST,
			   INPUTS_SPECTRO.ZP[ILAM], errmsg_ZP_T );

    }
    else {
      // legacy
      ZP_T = interp_1DFUN (OPT_INTERP, TEXPOSE_T, NBT, 
			   INPUTS_SPECTRO.TEXPOSE_LIST,
			   INPUTS_SPECTRO.ZP[ILAM], errmsg_ZP_T );
    }

    SQ_T = interp_1DFUN (OPT_INTERP, TEXPOSE_T, NBT, 
			 INPUTS_SPECTRO.TEXPOSE_LIST,
			 INPUTS_SPECTRO.SQSIGSKY[ILAM], errmsg_SQ_T );

    // Feb 2 2017: 
    //  scale template noise to search-image exposure. Nominaly
    //  the scaling would be TEXPOSE_S/TEXPOSE_T, but the noise
    //  has a more complex function; hence use the ZP to scale
    //  the template noise.  
    //   SQNOISE_SCALE = FLUXSCALE^2 = 10^( 0.8 * ZPdif )
    SQ_T *= pow( TEN, 0.8*(ZP_S-ZP_T) ) ;
     
  }


  arg     = -0.4*(GENMAG-ZP_S);
  Flux    = pow(TEN,arg) ;      // in p.e.

  SQ_SUM  = (SQ_S + SQ_T + Flux);
  if ( SQ_SUM >= 0.0 ) 
    {  FluxErr = sqrt(SQ_SUM);  SNR = Flux/FluxErr ;  }
  else
    { FluxErr = -9.0 ; }

  // check extrapolation beyond defined range of TEXPOSE (May 27 2020)
  if ( DO_TEXTRAP )
    { SNR *= sqrt(TEXPOSE_S / TEXPOSE_S_local); }

  if ( SQ_T >= 0.0 )
    {  *ERRFRAC_T = sqrt(SQ_T)/FluxErr ; } 
  
  if ( LDMP ) {
    print_preAbort_banner(fnam);
    printf(" xxx ILAM=%d LAM=%f \n", ILAM, INPUTS_SPECTRO.LAMAVG_LIST[ILAM] );
    printf(" xxx TEXPOSE_[S,T] = %f , %f \n", TEXPOSE_S, TEXPOSE_T );
    printf(" xxx GENMAG = %f \n", GENMAG);
    printf(" xxx SQ[S,T] = %le , %le    Flux=%le \n", SQ_S, SQ_T, Flux);    
    printf(" xxx SQ_SUM(SQ_S+SQ_T+Flux) = %le \n", SQ_SUM );

    int i;
    for(i=0;  i< NBT; i++ ) {
      printf(" xxx ZP-interpFun: i=%d: log10(Texpose)=%6.3f, ZP=%.3f \n", i,
             INPUTS_SPECTRO.LOGTEXPOSE_LIST[i],
             INPUTS_SPECTRO.ZP[ILAM][i] );
    }

    printf(" xxx ZP[S,T] interp values = %.3f , %.3f  \n", ZP_S, ZP_T );    

    printf(" xxx SNR = %f / %f = %f \n", Flux, FluxErr,SNR);
    printf(" xxx \n");


    fflush(stdout);
  }

  return(SNR);

} // end getSNR_spectrograph
