
/*
 * MWgaldust_av is an interface to schlegel's dust_getval for
 * use inside of snana.  All the inputs are set to default values
 * and the calling sequence sequence is internal rather than via the
 * command line.  The standard input is RA,DEC in degrees which are translated
 * to galactic coordinates which is what dust_getval expects.
 * Below is the original documentation of dust_getval, and all the
 * needed routines from Schlegel are combined in one file below this
 * simple interface subroutine.
 * D.Cinabro 6 June 2007
 *
 * Usage: void MWgaldust(double RA,double DEC,double *avgal, double galebmv)
 * RA,DEC are in degrees, avgal is the extnction in Sloan u,g,r,i,z filters,
 * and galebmv is the Milky Way galaxy E(B-V) redening from the Schlegal map
 */

/******************************************************************************/
/*
 * NAME:
 *   dust_getval
 *
 * PURPOSE:
 *   Read values from BH files or our dust maps.
 *
 *   Either the coordinates "gall" and "galb" must be set, or these coordinates
 *   must exist in the file "infile".  Output is written to standard output
 *   or the file "outfile".
 *
 * CALLING SEQUENCE:
 *   dust_getval gall galb map=map ipath=ipath interp=interp noloop=noloop \
 *    infile=infile outfile=outfile verbose=verbose
 *
 * OPTIONAL INPUTS:
 *   gall:       Galactic longitude(s) in degrees
 *   galb:       Galactic latitude(s) in degrees
 *   map:        Set to one of the following (default is "Ebv"):
 *               I100: 100-micron map in MJy/Sr
 *               X   : X-map, temperature-correction factor
 *               T   : Temperature map in degrees Kelvin for n=2 emissivity
 *               Ebv : E(B-V) in magnitudes
 *               mask: Mask values
 *   infile:     If set, then read "gall" and "galb" from this file
 *   outfile:    If set, then write results to this file
 *   interp:     Set this flag to "y" to return a linearly interpolated value
 *               from the 4 nearest pixels.
 *               This is disabled if map=='mask'.
 *   noloop:     Set this flag to "y" to read entire image into memory
 *               rather than reading pixel values for one point at a time.
 *               This is a faster option for reading a large number of values,
 *               but requires reading up to a 64 MB image at a time into
 *               memory.  (Actually, the smallest possible sub-image is read.)
 *   verbose:    Set this flag to "y" for verbose output, printing pixel
 *               coordinates and map values
 *   ipath:      Path name for dust maps; default to path set by the
 *               environment variable $DUST_DIR/maps, or to the current
 *               directory.
 *
 * EXAMPLES:
 *   Read the reddening value E(B-V) at Galactic (l,b)=(12,+34.5),
 *   interpolating from the nearest 4 pixels, and output to the screen:
 *   % dust_getval 12 34.5 interp=y
 *
 *   Read the temperature map at positions listed in the file "dave.in",
 *   interpolating from the nearest 4 pixels, and output to file "dave.out".
 *   The path name for the temperature maps is "/u/schlegel/".
 *   % dust_getval map=T ipath=/u/schlegel/ interp=y \
 *     infile=dave.in outfile=dave.out 
 *
 * DATA FILES FOR SFD MAPS:
 *   SFD_dust_4096_ngp.fits
 *   SFD_dust_4096_sgp.fits
 *   SFD_i100_4096_ngp.fits
 *   SFD_i100_4096_sgp.fits
 *   SFD_mask_4096_ngp.fits
 *   SFD_mask_4096_sgp.fits
 *   SFD_temp_ngp.fits
 *   SFD_temp_sgp.fits
 *   SFD_xmap_ngp.fits
 *   SFD_xmap_sgp.fits
 *
 * DATA FILES FOR BH MAPS:
 *   hinorth.dat
 *   hisouth.dat
 *   rednorth.dat
 *   redsouth.dat
 *
 * REVISION HISTORY:
 *   Written by D. Schlegel, 19 Jan 1998, Durham
 *   5-AUG-1998 Modified by DJS to read a default path from an environment
 *              variable $DUST_DIR/map.
 *
 *  Jan 28, 2007 D.Cinabro : add MWEBV = E(B-V) as additional output arg
 *
 *
 *  Oct 7, 2007 R.Kessler move contents of $DUST_DIR/maps to
 *              $SNDATA_ROOT/MWDUST/ and change local path 
 *
 * Mar 02, 2012: uchar    pLabel_temp[8] -> pLabel_temp[9] 
 *                  (bug found by S.Rodney)
 *
 * Feb 22, 2013: RK - update to compile with c++
 *
 * Sep 21 2013 RK - move GAL-related functions from sntools.c to here
 * 
 * Oct 29 2013 RK - move slalib routines to sntools.c
 *
 * Jan 28 2020 RK - abort if WAVE>12000 and using Fitz99 color law
 * 
 * Oct 9 2021 DB and DS - update Fitz/Odonell ratio and extend WAVE to 15000
 */
/**************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h> 
#include <ctype.h>

#include "MWgaldust.h"
#include "sntools.h"
//#include "sntools_cosmology.h"


// #####################################################
//
//   snana-GALextinct functions (moved from sntools)
//
// #####################################################



// *****************************
// mangled routines for fortran
double galextinct_(double *RV, double *AV, double *WAVE, int *OPT ) {
  return GALextinct(*RV, *AV, *WAVE, *OPT);
}
void text_mwoption__(char *nameOpt, int  *OPT, char *TEXT) {
  text_MWoption(nameOpt,*OPT,TEXT);
}
void modify_mwebv_sfd__(int *OPT, double *RA, double *DECL, 
			double *MWEBV, double *MWEBV_ERR) {
  modify_MWEBV_SFD(*OPT, *RA, *DECL, MWEBV, MWEBV_ERR) ;
}

// **********************************************
void text_MWoption(char *nameOpt, int OPT, char *TEXT) {

  // Created Sep 19 2013
  // Return corresponding *TEXT description of 
  // integer option OPT for 
  // *nameOpt = "MWCOLORLAW" or "COLORLAW" or "MWEBV" or "EBV"
  // ABORT on invalid OPT.

  char fnam[] = "text_MWoption" ;

  // ------------------ BEGIN ------------------

  sprintf(TEXT,"NULL");


  // ----------------------------------------
  if ( strcmp(nameOpt,"MWCOLORLAW") == 0  || 
       strcmp(nameOpt,"COLORLAW"  ) == 0 ) {

    if ( OPT == OPT_MWCOLORLAW_OFF ) 
      { sprintf(TEXT,"No Extinction");  }

    else if ( OPT == OPT_MWCOLORLAW_CCM89 ) 
      { sprintf(TEXT,"CCM89");  }

    else if ( OPT == OPT_MWCOLORLAW_ODON94 ) 
      { sprintf(TEXT,"CCM89+ODonell94");  }  

    else if ( OPT == OPT_MWCOLORLAW_FITZ99_APPROX ) 
      { sprintf(TEXT,"Fitzpatrick99 (approx fit to F99/ODonnel94)");  }

    else if ( OPT == OPT_MWCOLORLAW_FITZ99_EXACT ) 
      { sprintf(TEXT,"Fitzpatrick99 (exact cubic spline)");  }
    
    else if ( OPT == OPT_MWCOLORLAW_FITZ04 ) 
      { sprintf(TEXT,"Fitzpatrick04 (exact cubic spline)");  }

    else if ( OPT == OPT_MWCOLORLAW_GORD16 ) 
      { sprintf(TEXT,"Gordon16 (exact cubic spline)");  }

    else if ( OPT == OPT_MWCOLORLAW_GORD23 ) 
      { sprintf(TEXT,"Gordon23");  }

    else {
      sprintf(c1err,"Invalid OPT_MWCOLORLAW = %d", OPT);
      sprintf(c2err,"Check OPT_MWCOLORAW_* in sntools.h");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  }

  // ----------------------------------------
  else if ( strcmp(nameOpt,"MWEBV")== 0 || 
	    strcmp(nameOpt,"EBV"  )== 0 ) {

    if ( OPT == OPT_MWEBV_OFF ) 
      { sprintf(TEXT,"No Extinction");  }

    else if ( OPT == OPT_MWEBV_FILE ) 
      { sprintf(TEXT,"FILE value (SIMLIB or data header)");  }

    else if ( OPT == OPT_MWEBV_SFD98 ) 
      { sprintf(TEXT,"SFD98");  }

    else if ( OPT == OPT_MWEBV_Sch11_PS2013 ) 
      { sprintf(TEXT,"Schlafly11+PS2013: 0.86*MWEBV(SFD98)" );  }

    else {
      sprintf(c1err,"Invalid OPT_MWEBV = %d", OPT);
      sprintf(c2err,"Check OPT_MWEBV_* in sntools.h");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  }

  // ----------------------------------------
  else {
    sprintf(c1err,"Invalid nameOpt = %s", nameOpt );
    sprintf(c2err,"Valid nameOpt are COLORLAW and EBV");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
  }
  
} // end of text_MWoption

// **************************************
void modify_MWEBV_SFD(int OPT, double RA, double DECL, 
		      double *MWEBV, double *MWEBV_ERR) {

  // Created Sep 21 2013 by R.Kessler
  // According to input integer option OPT, modify MWEBV & MWEBV_ERR. 
  // Note that input MWEBV arguments are the FILE values, 
  // and may be changed !
  // Input RA and DEC may or may not be used depending on the option.

  double MWEBV_INP, MWEBV_OUT=0.0, MWEBV_ERR_OUT=0.0, MWEBV_SFD98=0.0 ;
  double dumXT[10] ;
  //  char fnam[] = "modify_MWEBV_SFD" ;

  // ----------- BEGIN -----------

  MWEBV_INP = *MWEBV ;
  MWEBV_OUT = -999.0 ;

  // check trival option with no Galactic extiction
  if ( OPT == OPT_MWEBV_OFF )  { 
    MWEBV_OUT = MWEBV_ERR_OUT = 0.0 ; 
    goto LOAD_OUTPUT ; 
  }  

  // ----------------------------------------------------
  // always compute MWEBV_SFD98 since many options need this.
  if ( OPT >= OPT_MWEBV_SFD98 )  
    {  MWgaldust(RA, DECL, dumXT, &MWEBV_SFD98 );  }
  else 
    { MWEBV_SFD98 = -999. ; }

  // ----------------------------------------------------

  if ( OPT == OPT_MWEBV_FILE )  {  
      // already defined extinction from file -> use it 
      MWEBV_OUT     = *MWEBV ;      // do nothing
      MWEBV_ERR_OUT = *MWEBV_ERR ;    
  }    

  else if ( OPT == OPT_MWEBV_SFD98 )  { 
    // force SFD98 regardless of input/FILE value
    MWEBV_OUT     = MWEBV_SFD98 ;
    MWEBV_ERR_OUT = MWEBV_SFD98/6.0 ;
  } 

  else if ( OPT == OPT_MWEBV_Sch11_PS2013 ) {

    MWEBV_OUT = 0.86 * MWEBV_SFD98 ;  // Apr 13 2018: Dan suggests this

    //    MWEBV_ERR_OUT = 0.1 * MWEBV_OUT ; // error is 10% of EBV
    MWEBV_ERR_OUT = 0.05 * MWEBV_OUT ; // Apr 13 2018

  }

  // load output
 LOAD_OUTPUT:
  *MWEBV     = MWEBV_OUT ;
  *MWEBV_ERR = MWEBV_ERR_OUT ;


} // end of modify_MWEBV_SFD


// **********************************************
double GALextinct(double RV, double AV, double WAVE, int OPT) {

/*** 
  
  Input : 
    AV   = V band (defined to be at 5495 Angstroms) extinction
    RV   = assumed A(V)/E(B-V) (e.g., = 3.1 in the LMC)
    WAVE = wavelength in angstroms

    OPT=89 => use original CCM89 :
              Cardelli, Clayton, & Mathis (1989) extinction law.

    OPT=94 => use update from O'Donell

    OPT=-99 => use Fitzpatrick 99 (PASP 111, 63) as implemented by 
              D.Scolnic with polynomial fit to ratio of F'99/O'94 vs. lambda.
              O'94 is the opt=94 option and F'99 was computed from 
              http://idlastro.gsfc.nasa.gov/ftp/pro/astro/fm_unred.pro
              Deprecated to OPT=-99 from OPT=99 September 25 2024.
              Only reliable for RV=3.1.

    OPT=99 => use Fitzpatrick 99 (PASP 111, 63) as implemented by S. Thorp.
                This version directly evaluates the cubic spline in inverse
                wavelength, as defined by the fm_unred.pro code. Consistent
                with extinction.py by K. Barbary, and BAYESN F99 implementation.
                Promoted to OPT=99 September 25 2024.

    OPT=204 => use Fitzpatrick 04 (ASP Conf. Ser. 309, 33) as implemented by S. Thorp.
                This uses the same curve as Fitzpatrick 99, but with updated
                behaviour in the IR. Checked against implementation by
                Gordon 24 (JOSS 9, 7023): github.com/karllark/dust_extinction.

    OPT=216 => use Gordon et al. 16 (ApJ 826, 104) as implemented by S. Thorp.
                This has an extra free parameter, FA. The final curve is a
                mixture of Fitzpatrick 99 and Gordon et al. 03 (ApJ 594, 279), 
                where the latter is SMC bar-like dust with RV=2.74 and no
                UV bump. FA=1 reverts to Fitzpatrick 99. FA=0 gives
                Gordon et al. 03. Effective RV = 1/[FA/RV + (1-FA)/2.74]
                Tested against Gordon 24 implementation.
                Currently has FA=0.0 hardcoded => Gordon et al. 03.

  Returns magnitudes of extinction.

 Nov 1, 2006: Add option to use new/old NIR coefficients
              (copied from Jha's MLCS code)

;     c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085,    $ ;Original
;                 0.01979, -0.77530,  0.32999 ]               ;coefficients
;     c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434,    $ ;from CCM89
;                -0.62251,  5.30260, -2.09002 ]

      c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137,    $    ;New coefficients
                 -1.718,   -0.827,    1.647, -0.505 ]        ;from O'Donnell
      c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985,    $    ;(1994)
                 11.102,    5.491,  -10.805,  3.347 ]

  Sep 18 2013 RK 
    - add opt=99 option to use Fitzpatrick 99 update.
    - reduce/remove pow calls to save CPU
    - rename CCMextinct -> GALextinct

  Aug 4 2019 RK
   + fix subtle bug by returning XT=0 only if AV=0, and not if AV<1E-9.
     Recall that negative AV are used for warping spectra in kcor.c.
     This bug caused all AV<0 to have same mag as AV=0.

  Sep 19 2024 ST
   + add an exact Fitzpatrick 99 implementation with opt=9999.

  Sep 25 2024 ST
   + Exact F'99 spline implementation promoted to opt=99
   - Old F'99 based on F'99/O'94 ratio deprecated to opt=-99

  Oct 19 2024 ST
   + Begun adding more dust laws
 ***/

  int i, DO94  ;
  double XT, x, y, a, b, fa, fb, xpow, xx, xx2, xx3 ;
  double y2, y3, y4, y5, y6, y7, y8 ;
  double FA ;
  char fnam[] = "GALextinct" ;

  // ------------------- BEGIN --------------

  XT = 0.0 ;
  FA = 0.0 ; // hardcoded for now

  if ( AV == 0.0  )  {  return XT ; }

  // -----------------------------------------
  // if seleting exact Fitz99-like option,
  // bypass everything else and call S. Thorp's functions
  
  if ( OPT == OPT_MWCOLORLAW_FITZ99_EXACT || OPT == OPT_MWCOLORLAW_FITZ04 )  {
    XT = GALextinct_Fitz99_exact(RV, AV, WAVE, OPT);
    return XT ;
  }
  else if ( OPT == OPT_MWCOLORLAW_GORD16 ) {
      double XTA, XTB;
      XTA = GALextinct_Fitz99_exact(RV, AV, WAVE, 99);
      XTB = GALextinct_Fitz99_exact(2.74, AV, WAVE, 203);
      return FA*XTA + (1-FA)*XTB ;
    }
  
  // -----------------------------------------
  DO94 = (OPT == OPT_MWCOLORLAW_ODON94 ||
	  OPT == OPT_MWCOLORLAW_FITZ99_APPROX ) ;

  x = 10000./WAVE;    // inverse wavelength in microns
  y = x - 1.82;

  if (x >= 0.3 && x < 1.1) {           // IR
    xpow = pow(x,1.61) ;
    a =  0.574 * xpow ;
    b = -0.527 * xpow ;
  } 
  else if (x >= 1.1 && x < 3.3) {    // Optical/NIR

    y2 = y  * y ;
    y3 = y2 * y ;
    y4 = y2 * y2 ;
    y5 = y3 * y2;
    y6 = y3 * y3;
    y7 = y4 * y3;
    y8 = y4 * y4;

    if ( DO94 ) {
    a = 1. + 0.104*y - 0.609*y2 + 0.701*y3 + 1.137*y4
      - 1.718*y5 - 0.827*y6 + 1.647*y7 -0.505*y8 ;

    b = 1.952*y + 2.908*y2 -3.989*y3 - 7.985*y4
      + 11.102*y5 + 5.491*y6 - 10.805*y7 + 3.347*y8;
    }
    else {
    a = 1. + 0.17699*y - 0.50447*y2 - 0.02427*y3
      + 0.72085*y4 + 0.01979*y5 - 0.77530*y6 + 0.32999*y7 ;

    b = 1.41338*y + 2.28305*y2 + 1.07233*y3 - 5.38434*y4
      - 0.62251*y5 + 5.30260*y6 - 2.09002*y7 ;
    }

  } 
  else if (x >= 3.3 && x < 8.0 ) {    // UV
    if (x >= 5.9) {
      xx  = x - 5.9 ;
      xx2 = xx  * xx ;
      xx3 = xx2 * xx ;

      fa = -0.04473*xx2 - 0.009779*xx3 ;
      fb =  0.21300*xx2 + 0.120700*xx3 ;

    } else {
      fa = fb = 0.0;
    }

    xx = x - 4.67 ; xx2 = (xx*xx);
    a =  1.752 - 0.316*x - 0.104/(xx2 + 0.341) + fa;

    xx = x - 4.62 ; xx2 = (xx*xx);
    b = -3.090 + 1.825*x + 1.206/(xx2 + 0.263) + fb;
  } 
  else if (x >= 8.0 && x <= 10.0) {  // Far-UV
    xx  = x - 8.0  ;
    xx2 = xx  * xx ;
    xx3 = xx2 * xx ; 

    a = -1.073 - 0.628*xx + 0.137*xx2 - 0.070*xx3 ;
    b = 13.670 + 4.257*xx - 0.420*xx2 + 0.374*xx3 ;
  } else {
    a = b = 0.0;
  }

  XT = AV*(a + b/RV);

  // Sep 18 2013 RK/DS - Check option for Fitzptrack 99

#define NPOLY_FITZ99 11 //Dillon and Dan upped to 10, Oct 9 2021
  if ( OPT == OPT_MWCOLORLAW_FITZ99_APPROX ) {  

    double XTcor, wpow[NPOLY_FITZ99] ;    
    double F99_over_O94[NPOLY_FITZ99] = {  // Dillon and Dan, Oct 9 2021
      8.55929205e-02,  1.91547833e+00, -1.65101945e+00,  7.50611119e-01,
      -2.00041118e-01,  3.30155576e-02, -3.46344458e-03,  2.30741420e-04,
      -9.43018242e-06,  2.14917977e-07, -2.08276810e-09
    };

    if ( WAVE > WAVEMAX_FITZ99  ) {
      sprintf(c1err,"Invalid WAVE=%.1f A for Fitzpatrick 99 color law.",
	      WAVE );
      sprintf(c2err,"Avoid NIR (>%.1f), or update Fitz99 in NIR",
	      WAVEMAX_FITZ99 );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // compute powers of wavelength without using slow 'pow' function
    wpow[0]  = 1.0 ;
    wpow[1]  = WAVE/1000. ;
    wpow[2]  = wpow[1] * wpow[1] ;
    wpow[3]  = wpow[1] * wpow[2] ;
    wpow[4]  = wpow[2] * wpow[2] ;
    wpow[5]  = wpow[3] * wpow[2] ;
    wpow[6]  = wpow[3] * wpow[3] ;
    wpow[7]  = wpow[4] * wpow[3] ;
    wpow[8]  = wpow[4] * wpow[4] ;
    wpow[9]  = wpow[5] * wpow[4] ;
    wpow[10]  = wpow[5] * wpow[5] ;

    XTcor = 0.0 ;
    for(i=0; i < NPOLY_FITZ99; i++ ) 
      {  XTcor += (wpow[i] * F99_over_O94[i]) ; }
    
    XT *= XTcor ;
  }

  return XT ;

}  // end of GALextinct


// ============= EXACT F99 EXTINCTION LAW ==============
double GALextinct_Fitz99_exact(double RV, double AV, double WAVE, int OPT) {
/*** 
  Created by S. Thorp, Sep 19 2024

  Default Fitzpatrick (1999) implementation since Sep 25 2024

  Also used to compute Fitzpatrick (2004), Gordon et al. (2003),
  and Gordon et al. (2016) laws.

  Input : 
    AV   = V band (defined to be at 5495 Angstroms) extinction
    RV   = assumed A(V)/E(B-V) (e.g., = 3.1 in the LMC)
    WAVE = wavelength in angstroms
    OPT  = Option from (99, 203, 204, 216)
Returns :
    XT = magnitudes of extinction
***/

    char fnam[] = "GALextinct_Fitz99_exact" ;

    if ( OPT == OPT_MWCOLORLAW_GORD03 && RV != 2.74 ) {
      sprintf(c1err,"Requested OPT=%d and RV=%.2f", OPT, RV);
      sprintf(c2err,"Gordon et al. 03 only valid for RV=2.74");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    //number of knots
    int Nk;
    // constants
    double x02, gamma2, c1, c2, c3, c4, c5;
    // target wavelength in inverse microns
    double x = 10000.0/WAVE;
    // spline result
    double y;

    // constants
    c2 = -0.824 + 4.717/RV;
    c5 = 5.90;
    if ( OPT == OPT_MWCOLORLAW_FITZ99_EXACT ) {
        x02 = 21.123216; // 4.596*4.596
        gamma2 = 0.9801; // 0.99*0.99
        c1 = 2.03 - 3.007*c2;
        c3 = 3.23;
        c4 = 0.41;
        Nk = 9;
    } 
    else if ( OPT == OPT_MWCOLORLAW_FITZ04 ) {
        x02 = 21.086464; // 4.592*4.592
        gamma2 = 0.850084; // 0.922*0.922
        c1 = 2.18 - 2.91*c2;
        c3 = 2.991;
        c4 = 0.319;
        Nk = 10;
    }
    else {
      sprintf(c1err,"Requested OPT=%d", OPT);
      sprintf(c2err,"Only 99, 204 are implemented!");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    if (WAVE <= 2700.0) { //FM90 curve in UV
        //extra terms
        double x2, y2, d, k;
        x2 = x*x;
        y = x2 - x02;
        d = x2 / (y*y + x2*gamma2);
        k = c1 + c2*x + c3*d;
        if (x >= c5) {
            y = x - c5;
            y2 = y * y;
            k += c4 * (0.5392*y2 + 0.05644*y2*y);
        }
        return AV * (1.0 + k/RV); 
    } else { //spline for optical/IR
        // extra constants
        double d1, d2, x82, x92, x8x02, x9x02;
        // powers of RV
        double RV2, RV3, RV4;
        // terms in the cubic spline equation
        double a, b, c, d, deltax, deltax2;

        // spline knot locations in inverse microns
        double xF[Nk];
        if ( OPT == OPT_MWCOLORLAW_FITZ04 ) {
            // Use FM07 knots for Fitzpatrick (2004) curve
            xF[0] = 0.0;
            xF[1] = 0.5;
            xF[2] = 0.75;
            xF[3] = 1.0;
        } else {
            xF[0] = 0.0;
            xF[1] = 1.0/2.65;
            xF[2] = 1.0/1.22;
        }
        xF[Nk-6] = 1.0/0.60;
        xF[Nk-5] = 1.0/0.547;
        xF[Nk-4] = 1.0/0.467; 
        xF[Nk-3] = 1.0/0.411;
        xF[Nk-2] = 1.0/0.270;
        xF[Nk-1] = 1.0/0.260;
        // spline knot values
        double yF[Nk];

        // y''/y matrix
        const double F04_KINVD[8][10] = {
                { 8.363001444, -29.445021665, 26.136103994, -5.423174767, 0.894742338,
                    -0.579592918, 0.067504416, -0.014025779, 0.001743637, -0.001280700 },
                { -2.178008666, 32.670129992, -60.816623962, 32.539048603, -5.368454030,
                    3.477557505, -0.405026495, 0.084154674, -0.010461821, 0.007684200 },
                { 0.349033220, -5.235498303, 25.130391856, -28.733019647, 20.579073782,
                    -13.330637103, 1.552601566, -0.322592918, 0.040103646, -0.029456099 },
                { -0.143088106, 2.146321587, -10.302343618, 17.313660802, -41.079282639,
                    35.355167970, -4.117769370, 0.855572522, -0.106361845, 0.078122696 },
                { 0.026683504, -0.400252565, 1.921212312, -3.228704024, 50.566388411,
                    -77.512222582, 35.824720884, -7.443507408, 0.925351342, -0.679669874 },
                { -0.007102699, 0.106540481, -0.511394311, 0.859426550, -13.459919650,
                    36.916413476, -45.296062632, 22.122269570, -2.750164769, 2.019993983 },
                { 0.000811556, -0.012173341, 0.058432036, -0.098198282, 1.537933619,
                    -4.218078179, 13.229095656, -13.261965375, 10.411053066, -7.646910755 },
                { -0.000364872, 0.005473077, -0.026270768, 0.044149485, -0.691447712,
                    1.896428085, -5.947739106, 7.633399832, -21.255428514, 18.341800494 }
        };
        const double F99_KINVD[7][9] = {
                { 10.250658602, -20.740327299, 11.788924891, -3.715854641, 2.664586253, 
                    -0.310340815, 0.064481288, -0.008016093, 0.005887814 },
                { -2.044608805, 10.254039610, -13.024855859, 13.772048671, -9.875739260, 
                    1.150214211, -0.238986593, 0.029709995, -0.021821969 },
                { 0.871617871, -4.371302789, 9.117884190, -31.624036073, 28.674517818, 
                    -3.339682936, 0.693905048, -0.086263899, 0.063360769 },
                { -0.162541946, 0.815173814, -1.700330722, 48.803145329, -76.266394662, 
                    35.679620964, -7.413359167, 0.921603416, -0.676917025 },
                { 0.043265924, -0.216985519, 0.452599357, -12.990574083, 36.584795098, 
                    -45.257439481, 22.114244617, -2.749167134, 2.019261221 },
                { -0.004943575, 0.024792817, -0.051714110, 1.484306083, -4.180187385, 
                    13.224682565, -13.261048442, 10.410939076, -7.646827029 },
                { 0.002222608, -0.011146734, 0.023250420, -0.667337025, 1.879392562,
                    -5.945755002, 7.632987583, -21.255377265, 18.341762851 }
        };

        // counters
        int i, q;

        // constants
        x82 = xF[Nk-2]*xF[Nk-2];
        x92 = xF[Nk-1]*xF[Nk-1];
        x8x02 = x82 - x02;
        x9x02 = x92 - x02;
        d1 = x82 / (x8x02*x8x02 + gamma2*x82);
        d2 = x92 / (x9x02*x9x02 + gamma2*x92);
        // powers of RV
        RV2 = RV*RV;
        RV3 = RV2*RV;
        RV4 = RV2*RV2;

        // RV-dependent spline knot values
        // polynomial coeffs match FM_UNRED.pro and extinction.py
        // NOTE: the optical coefficients differ from Gordon 24 implementation
        double yFNIR;
        if ( OPT == OPT_MWCOLORLAW_FITZ04 ) {
            yFNIR = (0.63*RV -0.84);
            yF[0] = -RV;
            yF[1] = yFNIR*pow(xF[1], 1.84) - RV;
            yF[2] = yFNIR*pow(xF[2], 1.84) - RV;
            yF[3] = yFNIR*pow(xF[3], 1.84) - RV;
        }
        else {
            yF[0] = -RV;
            yF[1] = -0.914616129*RV; // 0.26469*(RV/3.1) - RV
            yF[2] = -0.7325*RV; // 0.82925*(RV/3.1) - RV
        }
        yF[Nk-6] = -0.422809 + 0.00270*RV +  2.13572e-04*RV2;
        yF[Nk-5] = -5.13540e-02 + 0.00216*RV - 7.35778e-05*RV2;
        yF[Nk-4] =  7.00127e-01 + 0.00184*RV - 3.32598e-05*RV2;
        yF[Nk-3] =  1.19456 + 0.01707*RV - 5.46959e-03*RV2 +  
            7.97809e-04*RV3 - 4.45636e-05*RV4;
        yF[Nk-2] = c1 + c2*xF[Nk-2] + c3*d1;
        yF[Nk-1] = c1 + c2*xF[Nk-1] + c3*d2;


        // find index in knot list
        q = 0;
        while (q < Nk-1) {
            if (x < xF[q+1]) { 
                break; 
            } else {
                q++;
            }
        }

        // evaluate the spline
        deltax = xF[q+1] - xF[q];
        deltax2 = deltax * deltax;
        a = (xF[q+1] - x) / deltax;
        b = 1.0 - a;
        c = (a*a*a - a) * deltax2 / 6.0;
        d = (b*b*b - b) * deltax2 / 6.0;
        y = a*yF[q] + b*yF[q+1];
        // compute 2nd derivatives
        double d2yq = 0;
        double d2yq1 = 0;
        for (i=0; i<Nk; i++) {
            if (q > 0) {
                if ( OPT == OPT_MWCOLORLAW_FITZ04 ) {
                    d2yq  += F04_KINVD[q-1][i] * yF[i];
                } else {
                    d2yq  += F99_KINVD[q-1][i] * yF[i];
                }
            }
            if (q < Nk-1) {
                if ( OPT == OPT_MWCOLORLAW_FITZ04 ) {
                    d2yq1 += F04_KINVD[q][i] * yF[i];
                } else {
                    d2yq1 += F99_KINVD[q][i] * yF[i];
                }
            }
        }
        y += c*d2yq + d*d2yq1;

        return AV*(1.0 + y/RV);
    }

} // end of F99exact


// ========== FUNCTION TO RETURN EBV(SFD) =================
void MWgaldust(
	       double RA          // (I) RA
	       ,double DEC        // (I) DEC
	       ,double *galxtinct // (O) u,g,r,i,z extinction
	       ,double *galebmv   // (O) E(B-V)
	       )

{
   int      imap;
   int      qInterp;
   int      qVerbose;
   int      qNoloop;
   int      nGal;
   double    tmpl;
   double    tmpb;
   float *  pGall = NULL;
   float *  pGalb = NULL;
   float *  pMapval;
   double dustval;

   /* Declarations for keyword input values */
   char  *  pMapName;
   char  *  pIPath     = NULL ;
   char     pDefMap[]  = "Ebv" ;
   char     pDefPath[200];

   /* Declarations for data file names */
   char     pFileN[MAX_FILE_NAME_LEN];
   char     pFileS[MAX_FILE_NAME_LEN];
   struct   mapParms {
      char *   pName;
      char *   pFile1;
      char *   pFile2;
   } ppMapAll[] = {
     { "Ebv" , "SFD_dust_4096_ngp.fits", "SFD_dust_4096_sgp.fits" },
     { "I100", "SFD_i100_4096_ngp.fits", "SFD_i100_4096_sgp.fits" },
     { "X"   , "SFD_xmap_ngp.fits"     , "SFD_xmap_sgp.fits"      },
     { "T"   , "SFD_temp_ngp.fits"     , "SFD_temp_sgp.fits"      },
     { "mask", "SFD_mask_4096_ngp.fits", "SFD_mask_4096_sgp.fits" }
   };

   double RV[5];

   const int nmap = sizeof(ppMapAll) / sizeof(ppMapAll[0]);

   /* Set defaults */

   // xxx old   sprintf(pDefPath, "%s/maps", getenv("DUST_DIR") );
   sprintf(pDefPath, "%s/MWDUST", getenv("SNDATA_ROOT") );
   pIPath = pDefPath ;

   pMapName = pDefMap;
   qInterp = 1; /* interpolation */
   qVerbose = 0; /* not verbose */
   qNoloop = 0; /* do not read entire image into memory */

   RV[0] = 5.155; // u
   RV[1] = 3.793; // g
   RV[2] = 2.751; // r
   RV[3] = 2.086; // i
   RV[4] = 1.479; // z

   // Translate from RA and DEC to galactic

   slaEqgal(RA,DEC,&tmpl,&tmpb);
   nGal = 1;
   pGall = ccvector_build_(nGal);
   pGalb = ccvector_build_(nGal);
   pGall[0] = (float)tmpl;
   pGalb[0] = (float)tmpb;

   /* Determine the file names to use */
   for (imap=0; imap < nmap; imap++) {
      if (strcmp(pMapName,ppMapAll[imap].pName) == 0) {
	         sprintf(pFileN, "%s/%s", pIPath, ppMapAll[imap].pFile1);
	         sprintf(pFileS, "%s/%s", pIPath, ppMapAll[imap].pFile2);
      }
   }

   /* Read values from FITS files in Lambert projection */
   pMapval = lambert_getval(pFileN, pFileS, nGal, pGall, pGalb,
    qInterp, qNoloop, qVerbose);

   dustval = (double) *pMapval;
   galxtinct[0] = RV[0]*dustval;
   galxtinct[1] = RV[1]*dustval;
   galxtinct[2] = RV[2]*dustval;
   galxtinct[3] = RV[3]*dustval;
   galxtinct[4] = RV[4]*dustval;
   *galebmv      = dustval;
   return;

}



char Label_lam_nsgp[]  = "LAM_NSGP";
char Label_lam_scal[]  = "LAM_SCAL";

uchar *  label_lam_nsgp  = (uchar*)Label_lam_nsgp;
uchar *  label_lam_scal  = (uchar*)Label_lam_scal;

/******************************************************************************/
/* Fortran wrapper for reading Lambert FITS files */
void DECLARE(fort_lambert_getval)
  (char  *  pFileN,
   char  *  pFileS,
   void  *  pNGal,
   float *  pGall,
   float *  pGalb,
   void  *  pQInterp,
   void  *  pQNoloop,
   void  *  pQVerbose,
   float *  pOutput)
{
   int      iChar;
   int      qInterp;
   int      qNoloop;
   int      qVerbose;
   long     iGal;
   long     nGal;
   float *  pTemp;

   /* Truncate the Fortran-passed strings with a null,
    * in case they are padded with spaces */
   for (iChar=0; iChar < 80; iChar++)
    if (pFileN[iChar] == ' ') pFileN[iChar] = '\0';
   for (iChar=0; iChar < 80; iChar++)
    if (pFileS[iChar] == ' ') pFileS[iChar] = '\0';

   /* Select the 4-byte words passed by a Fortran call */
   if (sizeof(short) == 4) {
      nGal = *((short *)pNGal);
      qInterp = *((short *)pQInterp);
      qNoloop = *((short *)pQNoloop);
      qVerbose = *((short *)pQVerbose);
   } else if (sizeof(int) == 4) {
      nGal = *((int *)pNGal);
      qInterp = *((int *)pQInterp);
      qNoloop = *((int *)pQNoloop);
      qVerbose = *((int *)pQVerbose);
   } else if (sizeof(long) == 4) {
      nGal = *((long *)pNGal);
      qInterp = *((long *)pQInterp);
      qNoloop = *((long *)pQNoloop);
      qVerbose = *((long *)pQVerbose);
   }

   pTemp = lambert_getval(pFileN, pFileS, nGal, pGall, pGalb,
    qInterp, qNoloop, qVerbose);

   /* Copy results into Fortran-passed location for "pOutput",
    * assuming that memory has already been allocated */
   for (iGal=0; iGal < nGal; iGal++) pOutput[iGal] = pTemp[iGal];
}

/******************************************************************************/
/* Read one value at a time from NGP+SGP polar projections.
 * Set qInterp=1 to interpolate, or =0 otherwise.
 * Set qVerbose=1 to for verbose output, or =0 otherwise.
 */
float * lambert_getval
  (char  *  pFileN,
   char  *  pFileS,
   long     nGal,
   float *  pGall,
   float *  pGalb,
   int      qInterp,
   int      qNoloop,
   int      qVerbose)
{
   int      iloop;
   int      iGal;
   int      ii;
   int      jj;
   int   *  pNS; /* 0 for NGP, 1 for SGP */
   int      nIndx;
   int   *  pIndx;
   int   *  pXPix;
   int   *  pYPix;
   int      xPix;
   int      yPix;
   int      xsize;
   DSIZE    pStart[2];
   DSIZE    pEnd[2];
   DSIZE    nSubimg;
   float *  pSubimg;
   float    dx;
   float    dy;
   float    xr;
   float    yr;
   float    pWeight[4];
   float    mapval;
   float *  pOutput;
   float *  pDX = NULL;
   float *  pDY = NULL;

   /* Variables for FITS files */
   int      qRead;
   int      numAxes;
   DSIZE *  pNaxis;
   char  *  pFileIn = NULL;
   HSIZE    nHead;
   uchar *  pHead;

   /* Allocate output data array */
   pNS = ccivector_build_(nGal);
   pOutput = ccvector_build_(nGal);

   /* Decide if each point should be read from the NGP or SGP projection */
   for (iGal=0; iGal < nGal; iGal++)
    pNS[iGal] = (pGalb[iGal] >= 0.0) ? 0 : 1; /* ==0 for NGP, ==1 for SGP */

   if (qNoloop == 0) {  /* LOOP THROUGH ONE POINT AT A TIME */

      /* Loop through first NGP then SGP */
      for (iloop=0; iloop < 2; iloop++) {
         qRead = 0;

         /* Loop through each data point */
         for (iGal=0; iGal < nGal; iGal++) {
            if (pNS[iGal] == iloop) {

               /* Read FITS header for this projection if not yet read */
               if (qRead == 0) {
                  if (iloop == 0) pFileIn = pFileN; else pFileIn = pFileS;
                  fits_read_file_fits_header_only_(pFileIn, &nHead, &pHead);
                  qRead = 1;
               }

               if (qInterp == 0) {  /* NEAREST PIXELS */

                  /* Determine the nearest pixel coordinates */
                  lambert_lb2pix(pGall[iGal], pGalb[iGal], nHead, pHead,
                   &xPix, &yPix);

                  pStart[0] = xPix;
                  pStart[1] = yPix;

                  /* Read one pixel value from data file */
                  fits_read_point_(pFileIn, nHead, pHead, pStart, &mapval);
                  pOutput[iGal] = mapval;

                  if (qVerbose != 0)
                   printf("%8.3f %7.3f %1d %8d %8d %12.5e\n",
                   pGall[iGal], pGalb[iGal], iloop, xPix, yPix, mapval);

               } else {  /* INTERPOLATE */

                  fits_compute_axes_(&nHead, &pHead, &numAxes, &pNaxis);

                  /* Determine the fractional pixel coordinates */
                  lambert_lb2fpix(pGall[iGal], pGalb[iGal], nHead, pHead,
                   &xr, &yr);
/* The following 4 lines introduced an erroneous 1/2-pixel shift
   (DJS 18-Mar-1999).
                  xPix = (int)(xr-0.5);
                  yPix = (int)(yr-0.5);
                  dx = xPix - xr + 1.5;
                  dy = yPix - yr + 1.5;
 */
                  xPix = (int)(xr);
                  yPix = (int)(yr);
                  dx = xPix - xr + 1.0;
                  dy = yPix - yr + 1.0;

                  /* Force pixel values to fall within the image boundaries */
                  if (xPix < 0) { xPix = 0; dx = 1.0; }
                  if (yPix < 0) { yPix = 0; dy = 1.0; }
                  if (xPix >= pNaxis[0]-1) { xPix = pNaxis[0]-2; dx = 0.0; }
                  if (yPix >= pNaxis[1]-1) { yPix = pNaxis[1]-2; dy = 0.0; }

                  pStart[0] = xPix;
                  pStart[1] = yPix;
                  pEnd[0] = xPix + 1;
                  pEnd[1] = yPix + 1;

                  /* Create array of weights */
                  pWeight[0] =    dx  *    dy  ;
                  pWeight[1] = (1-dx) *    dy  ;
                  pWeight[2] =    dx  * (1-dy) ;
                  pWeight[3] = (1-dx) * (1-dy) ;

                  /* Read 2x2 array from data file */
                  fits_read_subimg_(pFileIn, nHead, pHead, pStart, pEnd,
                   &nSubimg, &pSubimg);

                  pOutput[iGal] = 0.0;
                  for (jj=0; jj < 4; jj++)
                   pOutput[iGal] += pWeight[jj] * pSubimg[jj];

                  fits_free_axes_(&numAxes, &pNaxis);
                  ccfree_((void **)&pSubimg);

                  if (qVerbose != 0)
                   printf("%8.3f %7.3f %1d %9.3f %9.3f %12.5e\n",
                   pGall[iGal], pGalb[iGal], iloop, xr, yr, pOutput[iGal]);

               }  /* -- END NEAREST PIXEL OR INTERPOLATE -- */
            }
         }
      }

   } else {  /* READ FULL IMAGE */

      pIndx = ccivector_build_(nGal);
      pXPix = ccivector_build_(nGal);
      pYPix = ccivector_build_(nGal);
      if (qInterp != 0) {
         pDX = ccvector_build_(nGal);
         pDY = ccvector_build_(nGal);
      }

      /* Loop through first NGP then SGP */
      for (iloop=0; iloop < 2; iloop++) {

         /* Determine the indices of data points in this hemisphere */
         nIndx = 0;
         for (iGal=0; iGal < nGal; iGal++) {
            if (pNS[iGal] == iloop) {
               pIndx[nIndx] = iGal;
               nIndx++;
            }
         }

         /* Do not continue if no data points in this hemisphere */
         if (nIndx > 0) {

            /* Read FITS header for this projection */
            if (iloop == 0) pFileIn = pFileN; else pFileIn = pFileS;
            fits_read_file_fits_header_only_(pFileIn, &nHead, &pHead);

            if (qInterp == 0) {  /* NEAREST PIXELS */

               /* Determine the nearest pixel coordinates */
               for (ii=0; ii < nIndx; ii++) {
                  lambert_lb2pix(pGall[pIndx[ii]], pGalb[pIndx[ii]],
                   nHead, pHead, &pXPix[ii], &pYPix[ii]);
               }

               pStart[0] = ivector_minimum(nIndx, pXPix);
               pEnd[0] = ivector_maximum(nIndx, pXPix);
               pStart[1] = ivector_minimum(nIndx, pYPix);
               pEnd[1] = ivector_maximum(nIndx, pYPix);

               /* Read smallest subimage containing all points in this hemi */
               fits_read_subimg_(pFileIn, nHead, pHead, pStart, pEnd,
                &nSubimg, &pSubimg);
               xsize = pEnd[0] - pStart[0] + 1;

               /* Determine data values */
               for (ii=0; ii < nIndx; ii++) {
                  pOutput[pIndx[ii]] = pSubimg[ pXPix[ii]-pStart[0] +
                   (pYPix[ii]-pStart[1]) * xsize ];

               }

               ccfree_((void **)&pSubimg);

            } else {  /* INTERPOLATE */

               fits_compute_axes_(&nHead, &pHead, &numAxes, &pNaxis);

               /* Determine the fractional pixel coordinates */
               for (ii=0; ii < nIndx; ii++) {
                  lambert_lb2fpix(pGall[pIndx[ii]], pGalb[pIndx[ii]],
                   nHead, pHead, &xr, &yr);
/* The following 4 lines introduced an erroneous 1/2-pixel shift
   (DJS 03-Mar-2004).
                  pXPix[ii] = (int)(xr-0.5);
                  pYPix[ii] = (int)(yr-0.5);
                  pDX[ii] = pXPix[ii] - xr + 1.5;
                  pDY[ii] = pYPix[ii] - yr + 1.5;
*/
                  pXPix[ii] = (int)(xr);
                  pYPix[ii] = (int)(yr);
                  pDX[ii] = pXPix[ii] - xr + 1.0;
                  pDY[ii] = pYPix[ii] - yr + 1.0;

                  /* Force pixel values to fall within the image boundaries */
                  if (pXPix[ii] < 0) { pXPix[ii] = 0; pDX[ii] = 1.0; }
                  if (pYPix[ii] < 0) { pYPix[ii] = 0; pDY[ii] = 1.0; }
                  if (pXPix[ii] >= pNaxis[0]-1)
                   { pXPix[ii] = pNaxis[0]-2; pDX[ii] = 0.0; }
                  if (pYPix[ii] >= pNaxis[1]-1)
                   { pYPix[ii] = pNaxis[1]-2; pDY[ii] = 0.0; }

               }

               pStart[0] = ivector_minimum(nIndx, pXPix);
               pEnd[0] = ivector_maximum(nIndx, pXPix) + 1;
               pStart[1] = ivector_minimum(nIndx, pYPix);
               pEnd[1] = ivector_maximum(nIndx, pYPix) + 1;

               /* Read smallest subimage containing all points in this hemi */
               fits_read_subimg_(pFileIn, nHead, pHead, pStart, pEnd,
                &nSubimg, &pSubimg);
               xsize = pEnd[0] - pStart[0] + 1;

               /* Determine data values */
               for (ii=0; ii < nIndx; ii++) {
                  /* Create array of weights */
                  pWeight[0] =    pDX[ii]  *    pDY[ii]  ;
                  pWeight[1] = (1-pDX[ii]) *    pDY[ii]  ;
                  pWeight[2] =    pDX[ii]  * (1-pDY[ii]) ;
                  pWeight[3] = (1-pDX[ii]) * (1-pDY[ii]) ;

                  pOutput[pIndx[ii]] =
                    pWeight[0] * pSubimg[
                     pXPix[ii]-pStart[0]   + (pYPix[ii]-pStart[1]  )*xsize ]
                   +pWeight[1] * pSubimg[
                     pXPix[ii]-pStart[0]+1 + (pYPix[ii]-pStart[1]  )*xsize ]
                   +pWeight[2] * pSubimg[
                     pXPix[ii]-pStart[0]   + (pYPix[ii]-pStart[1]+1)*xsize ]
                   +pWeight[3] * pSubimg[
                     pXPix[ii]-pStart[0]+1 + (pYPix[ii]-pStart[1]+1)*xsize ] ;

               }

               fits_free_axes_(&numAxes, &pNaxis);
               ccfree_((void **)&pSubimg);

            }  /* -- END NEAREST PIXEL OR INTERPOLATE -- */
         }

      }

      ccivector_free_(pIndx);
      ccivector_free_(pXPix);
      ccivector_free_(pYPix);
      if (qInterp != 0) {
         ccvector_free_(pDX);
         ccvector_free_(pDY);
      }
   }

   /* Free the memory allocated for the FITS header 
      (Moved outside previous brace by Chris Stoughton 19-Jan-1999) */
   fits_dispose_array_(&pHead);

   /* Deallocate output data array */
   ccivector_free_(pNS);

   return pOutput;
}

/******************************************************************************/
/* Transform from galactic (l,b) coordinates to fractional (x,y) pixel location.
 * Latitude runs clockwise from X-axis for NGP, counterclockwise for SGP.
 * This function returns the ZERO-INDEXED pixel position.
 * Updated 04-Mar-1999 to allow ZEA coordinate convention for the same
 * projection.
 */
void lambert_lb2fpix
  (float    gall,   /* Galactic longitude */
   float    galb,   /* Galactic latitude */
   HSIZE    nHead,
   uchar *  pHead,
   float *  pX,     /* X position in pixels from the center */
   float *  pY)     /* Y position in pixels from the center */
{
   int      q1;
   int      q2;
   int      nsgp;
   float    scale;
   float    crval1;
   float    crval2;
   float    crpix1;
   float    crpix2;
   float    cdelt1;
   float    cdelt2;
   float    cd1_1;
   float    cd1_2;
   float    cd2_1;
   float    cd2_2;
   float    lonpole;
   float    xr;
   float    yr;
   float    theta;
   float    phi;
   float    Rtheta;
   float    denom;
   static double dradeg = 180 / 3.1415926534;
   char  *  pCtype1;
   char  *  pCtype2;

   fits_get_card_string_(&pCtype1, label_ctype1, &nHead, &pHead);
   fits_get_card_string_(&pCtype2, label_ctype2, &nHead, &pHead);
   fits_get_card_rval_(&crval1, label_crval1, &nHead, &pHead);
   fits_get_card_rval_(&crval2, label_crval2, &nHead, &pHead);
   fits_get_card_rval_(&crpix1, label_crpix1, &nHead, &pHead);
   fits_get_card_rval_(&crpix2, label_crpix2, &nHead, &pHead);

   if (strcmp(pCtype1, "LAMBERT--X")  == 0 &&
       strcmp(pCtype2, "LAMBERT--Y")  == 0) {

      fits_get_card_ival_(&nsgp, label_lam_nsgp, &nHead, &pHead);
      fits_get_card_rval_(&scale, label_lam_scal, &nHead, &pHead);

      lambert_lb2xy(gall, galb, nsgp, scale, &xr, &yr);
      *pX = xr + crpix1 - crval1 - 1.0;
      *pY = yr + crpix2 - crval2 - 1.0;

   } else if (strcmp(pCtype1, "GLON-ZEA")  == 0 &&
              strcmp(pCtype2, "GLAT-ZEA")  == 0) { 

      q1 = fits_get_card_rval_(&cdelt1, label_cdelt1, &nHead, &pHead);
      q2 = fits_get_card_rval_(&cdelt2, label_cdelt2, &nHead, &pHead);
      if (q1 == TRUE_MWDUST && q2 == TRUE_MWDUST) {
          cd1_1 = cdelt1;
          cd1_2 = 0.0;
          cd2_1 = 0.0;
          cd2_2 = cdelt2;
       } else {
         fits_get_card_rval_(&cd1_1, label_cd1_1, &nHead, &pHead);
         fits_get_card_rval_(&cd1_2, label_cd1_2, &nHead, &pHead);
         fits_get_card_rval_(&cd2_1, label_cd2_1, &nHead, &pHead);
         fits_get_card_rval_(&cd2_2, label_cd2_2, &nHead, &pHead);
      }
      q1 = fits_get_card_rval_(&lonpole, label_lonpole, &nHead, &pHead);
      if (q1 == FALSE_MWDUST) lonpole = 180.0; /* default value */

      /* ROTATION */
      /* Equn (4) - degenerate case */
      if (crval2 > 89.9999) {
         theta = galb;
         phi = gall + 180.0 + lonpole - crval1;
      } else if (crval2 < -89.9999) {
         theta = -galb;
         phi = lonpole + crval1 - gall;
      } else {
         printf("ERROR: Unsupported projection!!!\n");
         /* Assume it's an NGP projection ... */
         theta = galb;
         phi = gall + 180.0 + lonpole - crval1;
      }

      /* Put phi in the range [0,360) degrees */
      phi = phi - 360.0 * floor(phi/360.0);

      /* FORWARD MAP PROJECTION */
      /* Equn (26) */
      Rtheta = 2.0 * dradeg * sin((0.5 / dradeg) * (90.0 - theta));

      /* Equns (10), (11) */
      xr = Rtheta * sin(phi / dradeg);
      yr = - Rtheta * cos(phi / dradeg);

      /* SCALE FROM PHYSICAL UNITS */
      /* Equn (3) after inverting the matrix */
      denom = cd1_1 * cd2_2 - cd1_2 * cd2_1;
      *pX = (cd2_2 * xr - cd1_2 * yr) / denom + (crpix1 - 1.0);
      *pY = (cd1_1 * yr - cd2_1 * xr) / denom + (crpix2 - 1.0);

   } else {

      *pX = -99.0;
      *pY = -99.0;

   }

   ccfree_((void **)&pCtype1);
   ccfree_((void **)&pCtype2);
}

/******************************************************************************/
/* Transform from galactic (l,b) coordinates to (ix,iy) pixel location.
 * Latitude runs clockwise from X-axis for NGP, counterclockwise for SGP.
 * This function returns the ZERO-INDEXED pixel position.
 * 
 */
void lambert_lb2pix
  (float    gall,   /* Galactic longitude */
   float    galb,   /* Galactic latitude */
   HSIZE    nHead,
   uchar *  pHead,
   int   *  pIX,    /* X position in pixels from the center */
   int   *  pIY)    /* Y position in pixels from the center */
{
   int      naxis1;
   int      naxis2;
   float    xr;
   float    yr;

   fits_get_card_ival_(&naxis1, label_naxis1, &nHead, &pHead);
   fits_get_card_ival_(&naxis2, label_naxis2, &nHead, &pHead);

   lambert_lb2fpix(gall, galb, nHead, pHead, &xr, &yr);
   *pIX = (int)floor(xr + 0.5);
   *pIY = (int)floor(yr + 0.5);

   /* Force bounds to be valid at edge, for ex at l=0,b=0 */
   //printf("NAXES %d %d\n", naxis1,naxis2);
   if (*pIX >= naxis1) *pIX = naxis1 - 1;
   if (*pIY >= naxis2) *pIY = naxis2 - 1;
}

/******************************************************************************/
/* Transform from galactic (l,b) coordinates to (x,y) coordinates from origin.
 * Latitude runs clockwise from X-axis for NGP, counterclockwise for SGP.
 */
void lambert_lb2xy
  (float    gall,   /* Galactic longitude */
   float    galb,   /* Galactic latitude */
   int      nsgp,   /* +1 for NGP projection, -1 for SGP */
   float    scale,  /* Radius of b=0 to b=90 degrees in pixels */
   float *  pX,     /* X position in pixels from the center */
   float *  pY)     /* Y position in pixels from the center */
{
   double   rho;
   static double dradeg = 180 / 3.1415926534;

   rho = sqrt(1. - nsgp * sin(galb/dradeg));
/* The following two lines were modified by Hans Schwengeler (17-Mar-1999)
   to get this to work on a Tur64 Unix 4.0E (DEC Alpha).  It appears that
   float and double on not the same on this machine.
   *pX = rho * cos(gall/dradeg) * scale;
   *pY = -nsgp * rho * sin(gall/dradeg) * scale;
*/
   *pX = rho * cos((float)((double)gall/dradeg)) * scale;
   *pY = -nsgp * rho * sin((float)((double)gall/dradeg)) * scale;
}
/******************************************************************************/
/* Find the min value of the nData elements of an integer array pData[].
 */
int ivector_minimum
  (int      nData,
   int   *  pData)
{
   int      i;
   int      vmin;

   vmin = pData[0];
   for (i=1; i<nData; i++) if (pData[i] < vmin) vmin=pData[i];

   return vmin;
}

/******************************************************************************/
/* Find the max value of the nData elements of an integer array pData[].
 */
int ivector_maximum
  (int      nData,
   int   *  pData)
{
   int      i;
   int      vmax;

   vmax = pData[0];
   for (i=1; i<nData; i++) if (pData[i] > vmax) vmax=pData[i];

   return vmax;
}


/**********************************************************************/
/*
 * Read an ASCII file as a 2-dimensional array of floating point numbers.
 * The number of columns is determined by the number of entries on the
 * first non-comment line.  Missing values are set to zero.
 * Comment lines (preceded with a hash mark) are ignored.
 * The returned matrix is in COLUMN-MAJOR order, and as such is
 * addressed as (*ppData)[iCol*NRows+iRow].
 * This is the Fortran storage scheme, but it is useful in addressing
 * a column as a vector that is contiguous in memory.
 * If the data array is empty, it should be passed with the value ppData = NULL.
 *
 * Return IO_GOOD if the file exists, and IO_FALSE_MWDUST otherwise.
 */
int asciifile_read_colmajor
  (char     pFileName[],
   int      numColsMax,
   int   *  pNRows,
   int   *  pNCols,
   float ** ppData)
{
   int      iCol;
   int      iRow;
   int      qExist;
   MEMSZ    memSize;
   float *  pNewData;

   qExist = asciifile_read_rowmajor(pFileName, numColsMax,
    pNRows, pNCols, ppData);

   if (qExist == IO_GOOD) {
      /* Create a new array of the same dimensions */
      memSize = sizeof(float) * (*pNRows)*(*pNCols);
      ccalloc_(&memSize, (void **)&pNewData);

      /* Copy the data into this array in column-major order */
      for (iCol=0; iCol < (*pNCols); iCol++) {
      for (iRow=0; iRow < (*pNRows); iRow ++) {
         pNewData[iCol*(*pNRows)+iRow] = (*ppData)[iRow*(*pNCols)+iCol];
      } }

      /* Toss out the old array */
      ccfree_((void **)ppData);
      *ppData = pNewData;
   }

   return qExist;
}

/******************************************************************************/
/*
 * Read an ASCII file as a 2-dimensional array of floating point numbers.
 * The number of columns is determined by the number of entries on the
 * first non-comment line.  Missing values are set to zero.
 * Comment lines (preceded with a hash mark) are ignored.
 * The returned matrix is in row-major order, and as such is
 * addressed as (*ppData)[iRow*NCols+iCol].
 * If the data array is empty, it should be passed with the value ppData = NULL.
 *
 * Return IO_GOOD if the file exists, and IO_BAD otherwise.
 */
int asciifile_read_rowmajor
  (char     pFileName[],
   int      numColsMax,
   int   *  pNRows,
   int   *  pNCols,
   float ** ppData)
{
   int      fileNum;
   int      iCol;
   int      nValues;
   int      qExist;
   const int numAddRows = 10;
   MEMSZ    memSize;
   MEMSZ    newMemSize;
   char  *  iq;
   float *  pValues;
   float *  pData;
   char     pPrivR[] = "r\0";

   *pNCols = 0;
   *pNRows = 0;

   qExist = inoutput_open_file(&fileNum, pFileName, pPrivR);

   if (qExist == IO_GOOD) {

      /* Allocate a starting block of memory for the data array */
      /* Start with enough memory for numAddRows rows */
      memSize = numAddRows * sizeof(float) * numColsMax;
      /* SHOULD BE ABLE TO USE REALLOC BELOW !!!??? */
      ccalloc_(&memSize, (void **)ppData);
      pData = *ppData;

      /* Allocate the temporary memory for each line of values */
      pValues = ccvector_build_(numColsMax);

      /* Read the first line, which determines the # of cols for all lines */
      iq = asciifile_read_line(fileNum, numColsMax, pNCols, pData);

      /* Read the remaining lines if a first line was read successfully */
      if (iq != NULL) {
         *pNRows=1;

         while ((iq = asciifile_read_line(fileNum, numColsMax,
          &nValues, pValues)) != NULL) {

            /* Allocate more memory for the data array if necessary */
            /* Always keep enough memory for at least one more row */
            newMemSize = sizeof(float) * (*pNRows + 1) * (*pNCols);
            if (newMemSize > memSize) {
               newMemSize = newMemSize + sizeof(float)*(numAddRows);
               ccalloc_resize_(&memSize, &newMemSize, (void **)ppData);
               pData = *ppData;
               memSize = newMemSize;
            }
   
            /* Case where the line contained fewer values than allowed */
            if (nValues < *pNCols) {
               for (iCol=0; iCol<nValues; iCol++)
                  pData[iCol+(*pNCols)*(*pNRows)] = pValues[iCol];
               for (iCol=nValues; iCol < *pNCols; iCol++)
                  pData[iCol+(*pNCols)*(*pNRows)] = 0;

            /* Case where line contained as many or more values than allowed */
            } else {
               for (iCol=0; iCol < *pNCols; iCol++)
                  pData[iCol+(*pNCols)*(*pNRows)] = pValues[iCol];
            }

            (*pNRows)++;
         }
      }

      inoutput_close_file(fileNum);
      ccvector_free_(pValues);
   }

   return qExist;
}

/******************************************************************************/
/*
 * Read all values from the next line of the input file.
 * Return NULL if end of file is reached.
 * Comment lines (preceded with a hash mark) are ignored.
 */
char * asciifile_read_line
  (int      filenum,
   int      numColsMax,
   int   *  pNValues,
   float *  pValues)
{
   char  *  pRetval;
   const    char whitespace[] = " \t\n";
   char     pTemp[MAX_FILE_LINE_LEN];
   char  *  pToken;

   /* Read next line from open file into temporary string */
   /* Ignore comment lines that are preceded with a hash mark */
   while (( (pRetval = fgets(pTemp, MAX_FILE_LINE_LEN, pFILEfits[filenum]))
    != NULL) && (pTemp[0] == '#') );

   if (pRetval != NULL) {
      /* Read one token at a time as a floating point value */
      *pNValues = 0;
      pToken = pTemp;
      while ((pToken = strtok(pToken, whitespace)) != NULL &&
             (*pNValues) < numColsMax ) {
         sscanf(pToken, "%f", &pValues[*pNValues]);
         pToken = NULL;
         (*pNValues)++;
      }
   }

   return pRetval;
}

/******************************************************************************/
/*
 * Subroutines to read and write FITS format files.
 * Note all variables are passed as pointers, so the routines can be called
 * by either C or Fortran programs.
 * Remember to omit the final underscore for calls from Fortran,
 * so one says 'call fits_add_card(...)' or 'i=fits_add_card(...)' in Fortran,
 * but 'i=fits_add_card_(...)' in C.
 *
 * D Schlegel -- Berkeley -- ANSI C
 * Mar 1992  DJS  Created
 * Dec 1993  DJS  Major revisions to allow dynamic memory allocations.
 */
/******************************************************************************/

#define min(a,b) ( ((a) < (b)) ? (a) : (b) )
#define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#define TRUE_MWDUST  1
#define FALSE_MWDUST 0

char Datum_zero[]    = "\0\0\0\0";
char Label_airmass[] = "AIRMASS ";
char Label_bitpix[]  = "BITPIX  ";
char Label_blank[]   = "BLANK   ";
char Label_bscale[]  = "BSCALE  ";
char Label_bzero[]   = "BZERO   ";
char Label_ctype1[]  = "CTYPE1  ";
char Label_ctype2[]  = "CTYPE2  ";
char Label_cdelt1[]  = "CDELT1  ";
char Label_cdelt2[]  = "CDELT2  ";
char Label_cd1_1[]   = "CD1_1   ";
char Label_cd1_2[]   = "CD1_2   ";
char Label_cd2_1[]   = "CD2_1   ";
char Label_cd2_2[]   = "CD2_2   ";
char Label_latpole[] = "LATPOLE ";
char Label_lonpole[] = "LONPOLE ";
char Label_crpix1[]  = "CRPIX1  ";
char Label_crpix2[]  = "CRPIX2  ";
char Label_crval1[]  = "CRVAL1  ";
char Label_crval2[]  = "CRVAL2  ";
char Label_date_obs[]= "DATE-OBS";
char Label_dec[]     = "DEC     ";
char Label_empty[]   = "        ";
char Label_end[]     = "END     ";
char Label_exposure[]= "EXPOSURE";
char Label_extend[]  = "EXTEND  ";
char Label_filtband[]= "FILTBAND";
char Label_filter[]  = "FILTER  ";
char Label_ha[]      = "HA      ";
char Label_instrume[]= "INSTRUME";
char Label_lamord[]  = "LAMORD  ";
char Label_loss[]    = "LOSS    ";
char Label_naxis[]   = "NAXIS   ";
char Label_naxis1[]  = "NAXIS1  ";
char Label_naxis2[]  = "NAXIS2  ";
char Label_object[]  = "OBJECT  ";
char Label_observer[]= "OBSERVER";
char Label_pa[]      = "PA      ";
char Label_platescl[]= "PLATESCL";
char Label_ra[]      = "RA      ";
char Label_rnoise[]  = "RNOISE  ";
char Label_rota[]    = "ROTA    ";
char Label_seeing[]  = "SEEING  ";
char Label_skyrms[]  = "SKYRMS  ";
char Label_skyval[]  = "SKYVAL  ";
char Label_slitwidt[]= "SLITWIDT";
char Label_st[]      = "ST      ";
char Label_telescop[]= "TELESCOP";
char Label_time[]    = "TIME    ";
char Label_tub[]     = "TUB     ";
char Label_ut[]      = "UT      ";
char Label_vhelio[]  = "VHELIO  ";
char Label_vminusi[] = "VMINUSI ";
char Card_simple[] =
   "SIMPLE  =                    T          "\
   "                                        ";
char Card_empty[] =
   "                                        "\
   "                                        ";
char Card_null[] =
   "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"\
   "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"\
   "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"\
   "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
char Card_end[] =
   "END                                     "\
   "                                        ";
char Text_T[] = "T";
char Text_F[] = "F";

uchar *  datum_zero    = (uchar*)Datum_zero;
uchar *  label_airmass = (uchar*)Label_airmass;
uchar *  label_bitpix  = (uchar*)Label_bitpix;
uchar *  label_blank   = (uchar*)Label_blank;
uchar *  label_bscale  = (uchar*)Label_bscale;
uchar *  label_bzero   = (uchar*)Label_bzero;
uchar *  label_ctype1  = (uchar*)Label_ctype1;
uchar *  label_ctype2  = (uchar*)Label_ctype2;
uchar *  label_cdelt1  = (uchar*)Label_cdelt1;
uchar *  label_cdelt2  = (uchar*)Label_cdelt2;
uchar *  label_cd1_1   = (uchar*)Label_cd1_1;
uchar *  label_cd1_2   = (uchar*)Label_cd1_2;
uchar *  label_cd2_1   = (uchar*)Label_cd2_1;
uchar *  label_cd2_2   = (uchar*)Label_cd2_2;
uchar *  label_latpole = (uchar*)Label_latpole;
uchar *  label_lonpole = (uchar*)Label_lonpole;
uchar *  label_crpix1  = (uchar*)Label_crpix1;
uchar *  label_crpix2  = (uchar*)Label_crpix2;
uchar *  label_crval1  = (uchar*)Label_crval1;
uchar *  label_crval2  = (uchar*)Label_crval2;
uchar *  label_date_obs= (uchar*)Label_date_obs;
uchar *  label_dec     = (uchar*)Label_dec;
uchar *  label_empty   = (uchar*)Label_empty;
uchar *  label_end     = (uchar*)Label_end;
uchar *  label_exposure= (uchar*)Label_exposure;
uchar *  label_extend  = (uchar*)Label_extend;
uchar *  label_filtband= (uchar*)Label_filtband;
uchar *  label_filter  = (uchar*)Label_filter;
uchar *  label_ha      = (uchar*)Label_ha;
uchar *  label_instrume= (uchar*)Label_instrume;
uchar *  label_lamord  = (uchar*)Label_lamord;
uchar *  label_loss    = (uchar*)Label_loss;
uchar *  label_naxis   = (uchar*)Label_naxis;
uchar *  label_naxis1  = (uchar*)Label_naxis1;
uchar *  label_naxis2  = (uchar*)Label_naxis2;
uchar *  label_object  = (uchar*)Label_object;
uchar *  label_observer= (uchar*)Label_observer;
uchar *  label_pa      = (uchar*)Label_pa;
uchar *  label_platescl= (uchar*)Label_platescl;
uchar *  label_ra      = (uchar*)Label_ra;
uchar *  label_rnoise  = (uchar*)Label_rnoise;
uchar *  label_rota    = (uchar*)Label_rota;
uchar *  label_seeing  = (uchar*)Label_seeing;
uchar *  label_skyrms  = (uchar*)Label_skyrms;
uchar *  label_skyval  = (uchar*)Label_skyval;
uchar *  label_slitwidt= (uchar*)Label_slitwidt;
uchar *  label_st      = (uchar*)Label_st;
uchar *  label_telescop= (uchar*)Label_telescop;
uchar *  label_time    = (uchar*)Label_time;
uchar *  label_tub     = (uchar*)Label_tub;
uchar *  label_ut      = (uchar*)Label_ut;
uchar *  label_vhelio  = (uchar*)Label_vhelio;
uchar *  label_vminusi = (uchar*)Label_vminusi;
uchar *  card_simple   = (uchar*)Card_simple;
uchar *  card_empty    = (uchar*)Card_empty;
uchar *  card_null     = (uchar*)Card_null;
uchar *  card_end      = (uchar*)Card_end;
uchar *  text_T        = (uchar*)Text_T;
uchar *  text_F        = (uchar*)Text_F;

/******************************************************************************/
/*
 * Read in FITS format data.  Assume the header is a multiple of
 * 2880-byte blocks, with the last block containing an END card.
 * Any data that follows is ignored.
 */
void fits_read_file_fits_header_only_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      fileNum;
   char     pPrivR[] = "r\0";

   inoutput_open_file(&fileNum, pFileName, pPrivR);

   /* Read header */
   fits_read_fits_header_(&fileNum, pNHead, ppHead);

   inoutput_close_file(fileNum);
}

/******************************************************************************/
/*
 * Read in header cards that are in ASCII format, for example as output
 * by IRAF with carraige returns after lines that are not even 80 characters.
 * The data is read into an array in FITS format.
 *
 * Return IO_GOOD if the file exists, and IO_BAD otherwise.
 */
int fits_read_file_ascii_header_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      cLen;
   int      fileNum;
   int      maxLen = 80;
   int      qExist;
   char     pPrivR[] = "r\0";
   uchar    pCard[80];

   qExist = inoutput_open_file(&fileNum, pFileName, pPrivR);

   if (qExist == IO_GOOD) {
      /* Read the header into memory until the end of file */
      *pNHead = 0;
      while (fgets((char *)pCard, maxLen, pFILEfits[fileNum]) != NULL) {
         /* Replace the /0 and remainder of card with blanks */
         for (cLen=strlen((const char *)pCard); cLen < 80; cLen++)
          pCard[cLen]=' ';

         fits_add_card_(pCard, pNHead, ppHead);
      }

      /* If no END card, then add one */
      if (fits_find_card_(label_end, pNHead, ppHead) == *pNHead) {
         fits_add_card_(card_end, pNHead, ppHead);
      }

      /* Close the file */
      inoutput_close_file(fileNum);
   }

   return qExist;
}

/******************************************************************************/
/*
 * Read in FITS format data.  Assume the header is a multiple of
 * 2880-byte blocks, with the last block containing an END card.
 * The data follows as either real values or as integer values
 * that should be scaled by the BZERO and BSCALE values.  The data
 * format is determined by the BITPIX card in the header.
 * The data is rescaled to 32-bit reals.
 * Also, the BITPIX card in the header is changed to -32.
 * Memory is dynamically allocated for the header and data arrays.
 *
 * Returned value is 0 unless the FITS file contains fewer data points
 * than indicated in the header, in which case the difference is returned.
 */
DSIZE fits_read_file_fits_r4_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   float ** ppData)
{
   int      bitpix;
   DSIZE    retval;

   retval = fits_read_file_fits_noscale_(pFileName, pNHead,
    ppHead, pNData, &bitpix, (uchar **)ppData);

   /* Convert data to real*4 if not already */
   fits_data_to_r4_(pNHead, ppHead, pNData, (uchar **)ppData);

   return retval;
}

/******************************************************************************/
/*
 * Read in FITS format data.  Assume the header is a multiple of
 * 2880-byte blocks, with the last block containing an END card.
 * The data follows as either real values or as integer values
 * that should be scaled by the BZERO and BSCALE values.  The data
 * format is determined by the BITPIX card in the header.
 * The data is rescaled to 16-bit integers.
 * Also, the BITPIX card in the header is changed to 16.
 * Memory is dynamically allocated for the header and data arrays.
 *
 * Returned value is 0 unless the FITS file contains fewer data points
 * than indicated in the header, in which case the difference is returned.
 */
DSIZE fits_read_file_fits_i2_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   short int   ** ppData)
{
   int      bitpix;
   DSIZE    retval;

   retval = fits_read_file_fits_noscale_(pFileName, pNHead,
    ppHead, pNData, &bitpix, (uchar **)ppData);

   /* Convert data to integer*4 if not already */
   fits_data_to_i2_(pNHead, ppHead, pNData, (uchar **)ppData);

   return retval;
}

/******************************************************************************/
/*
 * Read subimage from a FITS format data file, indexed from pStart to pEnd
 * in each dimension.
 *
 * The header is assumed to already be read using
 * fits_read_file_fits_header_only_(), to avoid reading it upon every
 * call to this routine.  The axis dimensions and BITPIX are read from
 * the header that is passed.  The dimensions of pLoc must agree with
 * the dimensions specified by NAXIS in this header.
 *
 * The data values are rescaled to 32-bit reals.
 *
 * Returned value is 0 unless the FITS file contains fewer data points
 * than requested, in which case the difference is returned.
 */

DSIZE fits_read_subimg_
  (char     pFileName[],
   HSIZE    nHead,
   uchar *  pHead,
   DSIZE *  pStart,
   DSIZE *  pEnd,
   DSIZE *  pNVal,
   float ** ppVal)
{
   int      bitpix;
   DSIZE    iloc;
   DSIZE    nExpect;
   int      size;
   MEMSZ    memSize;
   int      iAxis;
   int      numAxes;
   DSIZE *  pNaxis;
   float    bscale;
   float    bzero;
   uchar *  pData;

   int      fileNum;
   char     pPrivR[] = "r\0";

   inoutput_open_file(&fileNum, pFileName, pPrivR);

   /* Skip header */
   fits_skip_header_(&fileNum);

   /* From the given header, read BITPIX and PNAXIS */
   fits_get_card_ival_(&bitpix, label_bitpix, &nHead, &pHead);
   fits_compute_axes_(&nHead, &pHead, &numAxes, &pNaxis);

   /* Allocate memory for output */
   nExpect = 1;
   for (iAxis=0; iAxis < numAxes; iAxis++)
    nExpect *= (pEnd[iAxis] - pStart[iAxis] + 1);
   size = fits_size_from_bitpix_(&bitpix);
   memSize = size * nExpect;
   ccalloc_(&memSize, (void **)&pData);

   *pNVal = 0;
   fits_read_subimg1(numAxes, pNaxis, pStart, pEnd, fileNum, bitpix,
    pNVal, pData);
#ifdef LITTLE_ENDIAN
   fits_byteswap(bitpix, *pNVal, pData);
#endif

   /* Convert data to real*4 if not already */
   if (bitpix == -32) {
      *ppVal = (float *)pData;
   } else {
      /* Get the scaling parameters from the header */
      if (fits_get_card_rval_(&bscale, (uchar *)Label_bscale, &nHead, &pHead)
       == FALSE_MWDUST) {
         bscale = 1.0;  /* Default value for BSCALE */
      }
      if (fits_get_card_rval_(&bzero , (uchar *)Label_bzero , &nHead, &pHead)
       == FALSE_MWDUST) {
         bzero = 0.0;  /* Default value for BZERO */
      }
 
      memSize = sizeof(float) * nExpect;
      ccalloc_(&memSize, (void **)ppVal);
      for (iloc=0; iloc < *pNVal; iloc++)
       (*ppVal)[iloc] = fits_get_rval_(&iloc, &bitpix, &bscale, &bzero, &pData);
   }
 
   inoutput_close_file(fileNum);

   /* Plug a memory leak - Chris Stoughton 19-Jan-1999 */
   fits_free_axes_(&numAxes, &pNaxis);

   return (nExpect - (*pNVal));
}

void fits_read_subimg1
  (int      nel,
   DSIZE *  pNaxis,
   DSIZE *  pStart,
   DSIZE *  pEnd,
   int      fileNum,
   int      bitpix,
   DSIZE *  pNVal,
   uchar *  pData)
{
   int      iloop;
   int      ii;
   int      ipos;
   int      size;
   DSIZE    nskip;
   DSIZE    nread;
   FILE  *  pFILEin;

   pFILEin = pFILEfits[fileNum];
   size = fits_size_from_bitpix_(&bitpix);

   /* Skip "nskip" points */
   nskip = pStart[nel-1];
   for (ii=0; ii < nel-1; ii++) nskip = nskip * pNaxis[ii];
   ipos = ftell(pFILEin);
   fseek(pFILEin, (ipos + size*nskip), 0);

   if (nel > 1) {
      for (iloop=0; iloop < pEnd[nel-1]-pStart[nel-1]+1; iloop++)
       fits_read_subimg1(nel-1, pNaxis, pStart, pEnd, fileNum, bitpix,
        pNVal, pData);
   } else {
      nread = pEnd[0]-pStart[0]+1;

      /* Read in "nread" points */
      *pNVal += (int)fread(&pData[(*pNVal)*size], size, nread, pFILEin);
   }

   /* Skip "nskip" points */
   nskip = pNaxis[nel-1] - pEnd[nel-1] - 1;
   for (ii=0; ii < nel-1; ii++) nskip = nskip * pNaxis[ii];
   ipos = ftell(pFILEin);
   fseek(pFILEin, (ipos + size*nskip), 0);
}

/******************************************************************************/
/*
 * Read in one element from a FITS format data file, indexed by the
 * values in pLoc.
 *
 * The header is assumed to already be read using
 * fits_read_file_fits_header_only_(), to avoid reading it upon every
 * call to this routine.  The axis dimensions and BITPIX are read from
 * the header that is passed.  The dimensions of pLoc must agree with
 * the dimensions specified by NAXIS in this header.
 *
 * The data value is rescaled to a 32-bit real.
 *
 * Returned value is 0 unless the FITS file contains fewer data points
 * than requested (1), in which case the difference (1) is returned.
 */

DSIZE fits_read_point_
  (char     pFileName[],
   HSIZE    nHead,
   uchar *  pHead,
   DSIZE *  pLoc,
   float *  pValue)
{
   int      bitpix;
   DSIZE    iloc;
   int      nmult;
   int      size;
   MEMSZ    memSize;
   int      iAxis;
   int      numAxes;
   DSIZE *  pNaxis;
   float    bscale;
   float    bzero;
   uchar *  pData;
   DSIZE    retval;

   int      fileNum;
   int      ipos;
   char     pPrivR[] = "r\0";
   FILE  *  pFILEin;

   inoutput_open_file(&fileNum, pFileName, pPrivR);

   /* Skip header */
   fits_skip_header_(&fileNum);

   /* From the given header, read BITPIX and PNAXIS */
   fits_get_card_ival_(&bitpix, label_bitpix, &nHead, &pHead);
   fits_compute_axes_(&nHead, &pHead, &numAxes, &pNaxis);

   /* Find the 1-dimensional index for the data point requested */
   iloc = 0;
   nmult = 1;
   for (iAxis=0; iAxis < numAxes; iAxis++) {
      iloc = iloc + pLoc[iAxis] * nmult;
      nmult = nmult * pNaxis[iAxis];
   }

   /* Read one element from the data file */
   pFILEin = pFILEfits[fileNum];
   ipos = ftell(pFILEin);
   size = fits_size_from_bitpix_(&bitpix);
   memSize = size;
   ccalloc_(&memSize, (void **)&pData);
   fseek(pFILEin, (ipos + size*iloc), 0);
   retval = 1 - (int)fread(pData, size, 1, pFILEin);
#ifdef LITTLE_ENDIAN
   fits_byteswap(bitpix, 1, pData);
#endif

   /* Convert data to real*4 if not already */
   if (bitpix == -32) {
      *pValue = *( (float *)pData );
   } else {
      /* Get the scaling parameters from the header */
      if (fits_get_card_rval_(&bscale, (uchar *)Label_bscale, &nHead, &pHead)
       == FALSE_MWDUST) {
         bscale = 1.0;  /* Default value for BSCALE */
      }
      if (fits_get_card_rval_(&bzero , (uchar *)Label_bzero , &nHead, &pHead)
       == FALSE_MWDUST) {
         bzero = 0.0;  /* Default value for BZERO */
      }

      iloc = 0;
      *pValue = fits_get_rval_(&iloc, &bitpix, &bscale, &bzero, &pData);
   }

   inoutput_close_file(fileNum);

   /* Plug a memory leak - D. Schlegel 06-Feb-1999 */
   fits_free_axes_(&numAxes, &pNaxis);

   return retval;
}

/******************************************************************************/
/*
 * Read in FITS format data.  Assume the header is a multiple of
 * 2880-byte blocks, with the last block containing an END card.
 * The data follows as either real values or as integer values
 * that should be scaled by the BZERO and BSCALE values.  The data
 * format is determined by the BITPIX card in the header.
 * Memory is dynamically allocated for the header and data arrays.
 *
 * Returned value is 0 unless the FITS file contains fewer data points
 * than indicated in the header, in which case the difference is returned.
 */
DSIZE fits_read_file_fits_noscale_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   int   *  pBitpix,
   uchar ** ppData)
{
   int      fileNum;
   DSIZE    retval;
   char     pPrivR[] = "r\0";

   inoutput_open_file(&fileNum, pFileName, pPrivR);

   /* Read header */
   fits_read_fits_header_(&fileNum, pNHead, ppHead);

   /* From the header, read BITPIX and determine the number of data points */
   *pNData = fits_compute_ndata_(pNHead, ppHead);
   fits_get_card_ival_(pBitpix, label_bitpix, pNHead, ppHead);

   /* Read data */
   retval = fits_read_fits_data_(&fileNum, pBitpix, pNData, ppData);

   inoutput_close_file(fileNum);
   return retval;
}

/******************************************************************************/
/*
 * Read in EXTENDED FITS format data.  Assume the header is a multiple of
 * 2880-byte blocks, with the last block containing an END card.
 * In addition, allow for an extended header file.  (But always reads
 * exactly one additional header.)
 * The data follows as either real values or as integer values
 * that should be scaled by the BZERO and BSCALE values.  The data
 * format is determined by the BITPIX card in the header.
 * Memory is dynamically allocated for the header and data arrays.
 *
 * Returned value is 0 unless the FITS file contains fewer data points
 * than indicated in the header, in which case the difference is returned.
 */
DSIZE fits_read_file_xfits_noscale_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   HSIZE *  pNXhead,
   uchar ** ppXHead,
   DSIZE *  pNData,
   int   *  pBitpix,
   uchar ** ppData)
{

   int      fileNum;
   HSIZE    iCard;
   HSIZE *  pTempN;
   DSIZE    retval;
   uchar    pExtend[40];
   uchar *  pTempHead;
   char     pPrivR[] = "r\0";

   inoutput_open_file(&fileNum, pFileName, pPrivR);

   /* Read header */
   fits_read_fits_header_(&fileNum, pNHead, ppHead);

   /* Read extended header(if it exists) */
   pTempN = pNHead;
   pTempHead = *ppHead;
   iCard = fits_find_card_(label_extend, pNHead, ppHead);
   if (iCard < *pNHead)
     { sscanf( (const char*)&(*ppHead)[iCard*80+10], "%s", pExtend); }

   if (strcmp((const char*)pExtend, (const char *)text_T) == 0) {
      fits_read_fits_header_(&fileNum, pNXhead, ppXHead);
      pTempN = pNXhead;
      pTempHead = *ppXHead;
   }

   /* From the header, read BITPIX and determine the number of data points */
   /* (from the extended header if it exists) */
   *pNData = fits_compute_ndata_(pNHead, ppHead);
   fits_get_card_ival_(pBitpix, label_bitpix, pNXhead, ppXHead);

   /* Read data */
   retval = fits_read_fits_data_(&fileNum, pBitpix, pNData, ppData);

   inoutput_close_file(fileNum);
   return retval;
}

/******************************************************************************/
/*
 * Write FITS format data.  Write the header as a multiple of
 * 2880-byte blocks, with the last block containing an END card.
 * The data follows is written as 32-bit reals.  If necessary, the
 * data is converted to that format first.
 *
 * Returned value is 0 unless not all of the data points are written,
 * in which case the difference is returned.  (Does not work!!!???  Values
 * returned from fwrite() are bizarre!)
 */

DSIZE fits_write_file_fits_r4_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   float ** ppData)
{
   int      bitpix = -32;
   DSIZE    retval;

   fits_data_to_r4_(pNHead, ppHead, pNData, (uchar **)ppData);
   retval = fits_write_file_fits_noscale_(pFileName, pNHead,
    ppHead, pNData, &bitpix, (uchar **)ppData);
   return retval;
}

/******************************************************************************/
/*
 * Write FITS format data.  Write the header as a multiple of
 * 2880-byte blocks, with the last block containing an END card.
 * The data follows is written as 16-bit integers.  If necessary, the
 * data is converted to that format first.
 *
 * Returned value is 0 unless not all of the data points are written,
 * in which case the difference is returned.  (Does not work!!!???  Values
 * returned from fwrite() are bizarre!)
 */

DSIZE fits_write_file_fits_i2_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   short int ** ppData)
{
   int      bitpix = 16;
   DSIZE    retval;

   fits_data_to_i2_(pNHead, ppHead, pNData, (uchar **)ppData);
   retval = fits_write_file_fits_noscale_(pFileName, pNHead,
    ppHead, pNData, &bitpix, (uchar **)ppData);
   return retval;
}

/******************************************************************************/
/*
 * Write FITS format data.  Write the header as a multiple of
 * 2880-byte blocks, with the last block containing an END card.
 * The data follows as either real values or as integer values
 * that should be scaled by the BZERO and BSCALE values.  The data
 * format is determined by the BITPIX card in the header.
 *
 * Returned value is 0 unless not all of the data points are written,
 * in which case the difference is returned.  (Does not work!!!???  Values
 * returned from fwrite() are bizarre!)
 */

DSIZE fits_write_file_fits_noscale_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   int   *  pBitpix,
   uchar ** ppData)
{
   int      fileNum;
   DSIZE    retval;
   char     pPrivW[] = "w\0";

   inoutput_open_file(&fileNum,pFileName, pPrivW);

   /* Write header */
   fits_write_fits_header_(&fileNum, pNHead, ppHead);

   /* Write data */
   retval = fits_write_fits_data_(&fileNum, pBitpix, pNData, ppData);

   inoutput_close_file(fileNum);
   return retval;
}

/******************************************************************************/
/*
 * Read data blocks from an open FITS file.
 * One contiguous area of memory is dynamically allocated.
 *
 * Returned value is 0 unless the FITS file contains fewer data points
 * than indicated in the header, in which case the difference is returned.
 */
DSIZE fits_read_fits_data_
  (int   *  pFilenum,
   int   *  pBitpix,
   DSIZE *  pNData,
   uchar ** ppData)
{
   int      size;
   DSIZE    retval;

   /* Allocate the minimum number of 2880-byte blocks for the data */
   fits_create_fits_data_(pBitpix, pNData, ppData);

   /* Read the data until the number of data points or until the end
      of file is reached. */
   size = fits_size_from_bitpix_(pBitpix);
   retval = *pNData - (int)fread(*ppData, size, *pNData, pFILEfits[*pFilenum]);
#ifdef LITTLE_ENDIAN
   fits_byteswap(*pBitpix, *pNData, *ppData);
#endif

   return retval;
}

/******************************************************************************/
/*
 * Write data blocks to an open FITS file.
 *
 * Returned value is 0 unless not all of the data points are written,
 * in which case the difference is returned.  (Does not work!  Values
 * returned from fwrite() are bizarre!)
 */
DSIZE fits_write_fits_data_
  (int   *  pFilenum,
   int   *  pBitpix,
   DSIZE *  pNData,
   uchar ** ppData)
{
   int      i;
   int      j;
   int      k;
   int      size;
   int      retval;

   /* Write the number of data points indicated */
   size = fits_size_from_bitpix_(pBitpix);
#ifdef LITTLE_ENDIAN
   fits_byteswap(*pBitpix, *pNData, *ppData);
#endif
   retval = *pNData - (int)fwrite(*ppData, size, *pNData, pFILEfits[*pFilenum]);
#ifdef LITTLE_ENDIAN
   fits_byteswap(*pBitpix, *pNData, *ppData);
#endif
 
   /* Write some zeros such that the data takes up an integral number
      of 2880 byte blocks */
   j = (ftell(pFILEfits[*pFilenum]) % 2880)/size ;
   if (j != 0) {
      k = 1;
      for (i=j; i<(2880/size); i++)
       fwrite(datum_zero, size, k, pFILEfits[*pFilenum]);
   }

   return retval;
}

/******************************************************************************/
/*
 * Read header blocks from an open FITS file.
 * Memory for new blocks are dynamically allocated when needed.
 */
void fits_read_fits_header_
  (int   *  pFilenum,
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   HSIZE    jCard;
   uchar    pCard[80];

   /* Read the header into memory until the END card */
   *pNHead = 0;
   while (fits_get_next_card_(pFilenum, pCard)) {
      /* Only include this card if it is not blank */
      if (strncmp((const char *)pCard, (const char *)card_empty, 80) != 0) {
         fits_add_card_(pCard, pNHead, ppHead);
      }
   }
   fits_add_card_(card_end, pNHead, ppHead);
 
   /* Finish reading to the end of the last header block (the one w/END) */
   /* ignoring, and in effect deleting, any header cards after the END card */
   jCard = (ftell(pFILEfits[*pFilenum]) % 2880)/80 ;
   if (jCard != 0) {
      for (iCard=jCard; iCard<=35; iCard++) {
         fits_get_next_card_(pFilenum, pCard);
      }
   }

   /* Delete all cards where the label is blank */
   fits_purge_blank_cards_(pNHead, ppHead);

   /* Add missing cards to the FITS header */
   fits_add_required_cards_(pNHead, ppHead);
}

/******************************************************************************/
/*
 * Skip header blocks from an open FITS file.
 * This is a modified version of fits_read_fits_header_().
 */
void fits_skip_header_
  (int   *  pFilenum)
{
   HSIZE    iCard;
   HSIZE    jCard;
   uchar    pCard[80];

   /* Read the header into memory until the END card */
   while (fits_get_next_card_(pFilenum, pCard));
 
   /* Finish reading to the end of the last header block (the one w/END) */
   jCard = (ftell(pFILEfits[*pFilenum]) % 2880)/80 ;
   if (jCard != 0) {
      for (iCard=jCard; iCard<=35; iCard++) {
         fits_get_next_card_(pFilenum, pCard);
      }
   }
}

/******************************************************************************/
/*
 * Add any cards to the header that are required by the FITS definition
 * but are missing.
 */
void fits_add_required_cards_
  (HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      iAxis;
   int      numAxes;
   int      naxis;
   int      naxisX;
#if 0
   int      crpixX;
   float    crvalX;
   float    cdeltX;
#endif
   DSIZE *  pNaxis;
   uchar    pLabel_temp[9]; /* Must be long enough for 8 chars + NULL */

   if (fits_get_card_ival_(&naxis, label_naxis, pNHead, ppHead) == FALSE_MWDUST) {
      naxis = 0; /* default to no data axes */
      fits_change_card_ival_(&naxis, label_naxis, pNHead, ppHead);
   }

   fits_compute_axes_(pNHead, ppHead, &numAxes, &pNaxis);

   for (iAxis=0; iAxis < numAxes; iAxis++) {
      /* For each axis, be sure that a NAXISx, CRPIXx, CRVALx and CDELTx
       * card exists.  If one does not exist, then create it.
       * Create the labels for each axis for which to look as pLabel_temp.
       * Be certain to pad with spaces so that a NULL is not written.
       */

      sprintf((char *)pLabel_temp, "NAXIS%-3d", iAxis+1);
      if (fits_get_card_ival_(&naxisX, pLabel_temp, pNHead, ppHead) == FALSE_MWDUST) {
         naxisX = 1; /* default to 1 */
         fits_change_card_ival_(&naxisX, pLabel_temp, pNHead, ppHead);
         printf("Adding a card %s\n", pLabel_temp);
      }

#if 0
      sprintf(pLabel_temp, "CRPIX%-3d  ", iAxis+1);
      if (fits_get_card_ival_(&crpixX, pLabel_temp, pNHead, ppHead) == FALSE_MWDUST) {
         crpixX = 1; /* default to start numbering at the first pixel */
         fits_change_card_ival_(&crpixX, pLabel_temp, pNHead, ppHead);
         printf("Adding a card %s\n", pLabel_temp);
      }

      sprintf(pLabel_temp, "CRVAL%-3d  ", iAxis+1);
      if (fits_get_card_rval_(&crvalX, pLabel_temp, pNHead, ppHead) == FALSE_MWDUST) {
         crvalX = 0.0; /* default to the first pixel value to be zero */
         fits_change_card_rval_(&crvalX, pLabel_temp, pNHead, ppHead);
         printf("Adding a card %s\n", pLabel_temp);
      }

      sprintf(pLabel_temp, "CDELT%-3d  ", iAxis+1);
      if (fits_get_card_rval_(&cdeltX, pLabel_temp, pNHead, ppHead) == FALSE_MWDUST) {
         cdeltX = 1.0; /* default to spacing each pixel by a value of 1 */
         fits_change_card_rval_(&cdeltX, pLabel_temp, pNHead, ppHead);
         printf("Adding a card %s\n", pLabel_temp);
      }
#endif
   }

   /* Plug a memory leak - Chris Stoughton 19-Jan-1999 */
   fits_free_axes_(&numAxes, &pNaxis);
}

/******************************************************************************/
/*
 * Write header blocks to an open FITS file.
 */
void fits_write_fits_header_
  (int   *  pFilenum,
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      iCard;
   int      jCard;
   uchar *  pHead = *ppHead;

   /* Write the number of header cards indicated */
   for (iCard=0; iCard < *pNHead; iCard++) {
      fits_put_next_card_(pFilenum, &pHead[iCard*80]);
   }
 
   /* Write some more blank cards such that the header takes up an
      integral number of 2880 byte blocks */
   jCard = (ftell(pFILEfits[*pFilenum]) % 2880)/80 ;
   if (jCard != 0) {
      for (iCard=jCard; iCard <= 35; iCard++) {
         fits_put_next_card_(pFilenum, card_empty);
      }
   }
}

/******************************************************************************/
/*
 * Create a FITS header that only contains an END card.
 * Memory for new blocks are dynamically allocated when needed.
 */
void fits_create_fits_header_
  (HSIZE *  pNHead,
   uchar ** ppHead)
{
   /* First dispose of any memory already allocated by ppHead. */
   fits_dispose_array_(ppHead);

   /* Create a header with only a SIMPLE and END card.
    * Note that an entire 2880 byte block will be created
    * by the call to fits_add_card_().
    */
   *pNHead = 0;
   fits_add_card_(card_end, pNHead, ppHead);
   fits_add_card_(card_simple, pNHead, ppHead);
}

/******************************************************************************/
/*
 * Copy a FITS header into a newly created array.  Dynamically allocate
 * memory for the new header.
 */
void fits_duplicate_fits_header_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   uchar ** ppHeadCopy)
{
   MEMSZ    memSize;
 
   /* Allocate the minimum number of 2880-byte blocks for the header */
   memSize = ((int)((80*(*pNHead)-1)/2880) + 1) * 2880;
   ccalloc_(&memSize, (void **)ppHeadCopy);

   /* Copy all of the header bytes verbatim */
   memmove((void *)(*ppHeadCopy), (const void *)(*ppHead), memSize);
}

/******************************************************************************/
/*
 * Copy a FITS data array of real*4 into a newly created array.
 */
void fits_duplicate_fits_data_r4_
  (DSIZE *  pNData,
   float ** ppData,
   float ** ppDataCopy)
{
   int      bitpix = -32;

   fits_duplicate_fits_data_(&bitpix, pNData, (uchar **)ppData,
    (uchar **)ppDataCopy);
}

/******************************************************************************/
/*
 * Copy a FITS data array into a newly created array.
 */
void fits_duplicate_fits_data_
  (int   *  pBitpix,
   DSIZE *  pNData,
   uchar ** ppData,
   uchar ** ppDataCopy)
{
   int      size;
   MEMSZ    memSize;
 
   /* Allocate the minimum number of 2880-byte blocks for the header */
   size = fits_size_from_bitpix_(pBitpix);
   memSize = ((int)((size*(*pNData)-1)/2880) + 1) * 2880;
   ccalloc_(&memSize, (void **)ppDataCopy);

   /* Copy all of the data bytes verbatim */
   memmove((void *)(*ppDataCopy), (const void *)(*ppData), memSize);
}

/******************************************************************************/
/*
 * Create a FITS data array of real*4 with the number of elements specified.
 * Memory for new blocks are dynamically allocated when needed.
 */
void fits_create_fits_data_r4_
  (DSIZE *  pNData,
   float ** ppData)
{
   int      bitpix = -32;
   fits_create_fits_data_(&bitpix, pNData, (uchar **)ppData);
}

/******************************************************************************/
/*
 * Create a FITS data array with the number of elements specified.
 * Memory for new blocks are dynamically allocated when needed.
 */
void fits_create_fits_data_
  (int   *  pBitpix,
   DSIZE *  pNData,
   uchar ** ppData)
{
   int      size;
   MEMSZ    memSize;

   /* Allocate the minimum number of 2880-byte blocks for the data */
   size = fits_size_from_bitpix_(pBitpix);
   memSize = ((int)((size*(*pNData)-1)/2880) + 1) * 2880;
   ccalloc_(&memSize, (void **)ppData);
}

/******************************************************************************/
/*
 * Free the memory allocated for the header and data arrays.
 * Return TRUE_MWDUST if both arrays existed and were freed, and FALSE_MWDUST otherwise.
 */
int fits_dispose_header_and_data_
  (uchar ** ppHead,
   uchar ** ppData)
{
   int      retval;

   retval = fits_dispose_array_(ppHead) && fits_dispose_array_(ppData);
   return retval;
}

/******************************************************************************/
/*
 * Free the memory allocated for a FITS header or data array.
 * Return TRUE_MWDUST if the array existed and was freed, and FALSE_MWDUST otherwise.
 */
int fits_dispose_array_
  (uchar ** ppHeadOrData)
{
   int      retval;

   retval = FALSE_MWDUST;
   if (*ppHeadOrData != NULL) {
      ccfree_((void **)ppHeadOrData);
      retval = TRUE_MWDUST;
   }
   return retval;
}

/******************************************************************************/
/*
 * Compute the total number of data points.
 * This information is determined from the header cards NAXIS and NAXISx.
 */
DSIZE fits_compute_ndata_
  (HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      numAxes;
   DSIZE    iAxis;
   DSIZE *  pNaxis;
   DSIZE    nData;

   fits_compute_axes_(pNHead, ppHead, &numAxes, &pNaxis);
   if (numAxes == 0)
      nData = 0;
   else {
      nData = 1;
      for (iAxis=0; iAxis < numAxes; iAxis++) nData *= pNaxis[iAxis];
   }

   /* Plug a memory leak - D. Schlegel 06-Feb-1999 */
   fits_free_axes_(&numAxes, &pNaxis);

   return nData;
}

/******************************************************************************/
/*
 * Compute the number of axes and the dimension of each axis.
 * This information is determined from the header cards NAXIS and NAXISx.
 */
void fits_compute_axes_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   int   *  pNumAxes,
   DSIZE ** ppNaxis)
{
   int      iAxis;
   int      ival;
   DSIZE *  pNaxis;
   MEMSZ    memSize;
   uchar    pLabel_temp[9];

   fits_get_card_ival_(pNumAxes, label_naxis, pNHead, ppHead);
   if (*pNumAxes > 0) {
      memSize = (*pNumAxes) * sizeof(DSIZE);
      ccalloc_(&memSize, (void **)ppNaxis);
      pNaxis = *ppNaxis;
      for (iAxis=0; iAxis < *pNumAxes; iAxis++) {
         /* Create the label for this axis for which to look.
          * Be certain to pad with spaces so that a NULL is not written.
          */
         sprintf((char *)pLabel_temp, "NAXIS%d  ", iAxis+1);
         fits_get_card_ival_(&ival, pLabel_temp, pNHead, ppHead);
         pNaxis[iAxis] = ival;
      }
   }
}
       
/******************************************************************************/
/*
 * Free memory for axes dimensions allocated by "fits_compute_axes_".
 */
void fits_free_axes_
  (int   *  pNumAxes,
   DSIZE ** ppNaxis)
{
   if (*pNumAxes > 0) {
      ccfree_((void **)ppNaxis);
   }
}

/******************************************************************************/
/*
 * Compute the wavelength for a given pixel number using Vista coefficients.
 * This must be preceded with a call to fits_compute_vista_poly_coeffs_
 * to find the polynomial coefficients from the header cards.
 * The first element of pCoeff is a central pixel number for the fit
 * and the remaining LAMORD elements are the polynomial coefficients.
 * The wavelength of pixel number iPix (zero-indexed):
 *   wavelength(iPix) = SUM{j=1,nCoeff-1} Coeff(j) * [iPix - Coeff(0)]**(j-1)
 */
float compute_vista_wavelength_
  (DSIZE *  pPixelNumber,
   int   *  pNCoeff,
   float ** ppCoeff)
{
   int      iCoeff;
   DSIZE    centralPixelNumber;
   float    wavelength;

   centralPixelNumber = (DSIZE)(*ppCoeff)[0];
   wavelength = 0.0;
   for (iCoeff=1; iCoeff < *pNCoeff; iCoeff++) {
      wavelength += (*ppCoeff)[iCoeff]
       * pow(*pPixelNumber - centralPixelNumber, (float)(iCoeff-1));
   }
   return wavelength;
}

/******************************************************************************/
/*
 * Compute the number of coefficients for a polynomial wavelength fit
 * and the values of those coefficients.
 * This information is determined from the header cards LAMORD and LPOLYx.
 * Set nCoeff=LAMORD+1, and an array ppCoeff is created that has LAMORD+1
 * elements.  The first element is a central pixel number for the fit
 * and the remaining LAMORD elements are the polynomial coefficients.
 * The wavelength of pixel number iPix (zero-indexed):
 *   wavelength(iPix) = SUM{j=1,nCoeff-1} Coeff(j) * [iPix - Coeff(0)]**(j-1)
 * The coefficients are stored 4 on a line, so that the LPOLY0 card
 * contains up to the first 4 coefficients, and LPOLY1 up to the next 4, etc.
 */
void fits_compute_vista_poly_coeffs_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   int   *  pNCoeff,
   float ** ppCoeff)
{
   int      iCoeff;
   int      iLpolyNum;
   int      nLpolyNum;
   MEMSZ    memSize;
   uchar    pLabel_temp[9]; /* Must be long enough for 8 chars + NULL */
   char  *  pStringVal;
   char  *  pCardLoc;
   const char pCharSpace[] = " \'";

   fits_get_card_ival_(pNCoeff,label_lamord,pNHead,ppHead);
   if (*pNCoeff > 0) {
      (*pNCoeff)++;
      memSize = (*pNCoeff) * sizeof(float);
      ccalloc_(&memSize, (void **)ppCoeff);
      nLpolyNum = (*pNCoeff+3) / 4;
      iCoeff = 0;
      for (iLpolyNum=0; iLpolyNum < nLpolyNum; iLpolyNum++) {
         /* Create the label for this coefficient for which to look.
          * Be certain to pad with spaces so that a NULL is not written.
          */
         sprintf((char *)pLabel_temp, "LPOLY%-3d  ", iLpolyNum);
         fits_get_card_string_(&pStringVal, pLabel_temp, pNHead, ppHead);
         pCardLoc = pStringVal;
         for (iCoeff=iLpolyNum*4; iCoeff < min(iLpolyNum*4+4, *pNCoeff);
          iCoeff++) {
            sscanf(strtok(pCardLoc,pCharSpace), "%f", &(*ppCoeff)[iCoeff]);
            pCardLoc=NULL;
         }
      }
   }
}
       
/******************************************************************************/
/* 
 * Convert a data array to real*4 data, if it is not already.
 * A new array is created for the data, and the old array is discarded.
 * Change the BITPIX card in the header to -32 to indicate the data is real*4.
 */
void fits_data_to_r4_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   uchar ** ppData)
{
   int      bitpix;
   int      newBitpix;
   DSIZE    iData;
   HSIZE    iCard;
   MEMSZ    memSize;
   float    bscale;
   float    bzero;
   float    blankval;
   float    newBlankval;
   float *  pNewData;

   fits_get_card_ival_(&bitpix, label_bitpix, pNHead, ppHead);

   /* Convert data to real*4 if not already */
   if (bitpix != -32) {

      /* Get the scaling parameters from the header */
      if (fits_get_card_rval_(&bscale, (uchar *)Label_bscale, pNHead, ppHead)
       == FALSE_MWDUST) {
         bscale = 1.0;  /* Default value for BSCALE */
      }
      if (fits_get_card_rval_(&bzero , (uchar *)Label_bzero , pNHead, ppHead)
       == FALSE_MWDUST) {
         bzero = 0.0;  /* Default value for BZERO */
      }

      /* Allocate the minimum number of 2880-byte blocks for the data */
      memSize = ((int)((4*(*pNData)-1)/2880) + 1) * 2880;
      ccalloc_(&memSize, (void **)&pNewData);

      /* Convert the data and write to the new array */
      /* Note that nothing is done to rescale BLANK values properly */
      for (iData=0; iData < *pNData; iData++) {
         pNewData[iData] =
          fits_get_rval_(&iData, &bitpix, &bscale, &bzero, ppData);
      }

      /* Free the memory from the old array, and change the ppData pointer
         to point to the new array */
      ccfree_((void **)ppData);
      *ppData = (uchar *)pNewData;

      /* Change the BITPIX card to -32, indicating the data is real*4 */
      newBitpix = -32;
      fits_change_card_ival_(&newBitpix, label_bitpix, pNHead, ppHead);

      /* Delete the BSCALE and BZERO cards which are no longer used */
      fits_delete_card_(label_bscale, pNHead, ppHead);
      fits_delete_card_(label_bzero , pNHead, ppHead);

      /* Rescale the BLANK card if it exists */
      if ((iCard = fits_find_card_(label_blank, pNHead, ppHead)) != *pNHead) {
         fits_get_card_rval_(&blankval, label_blank, pNHead, ppHead);
         if      (bitpix ==  8) newBlankval = blankval * bscale + bzero;
         else if (bitpix == 16) newBlankval = blankval * bscale + bzero;
         else if (bitpix == 32) newBlankval = blankval * bscale + bzero;
         else if (bitpix == -8) newBlankval = blankval;
         else if (bitpix ==-32) newBlankval = blankval;
         else if (bitpix ==-64) newBlankval = blankval;
         fits_change_card_rval_(&newBlankval, label_blank, pNHead, ppHead);
      }

   }
}

/******************************************************************************/
/* 
 * Convert a data array to integer*2 data, if it is not already.
 * A new array is created for the data, and the old array is discarded.
 * Change the BITPIX card in the header to 16 to indicate the data is integer*2.
 */
void fits_data_to_i2_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   uchar ** ppData)
{
   int      bitpix;
   int      newBitpix;
   DSIZE    iData;
   HSIZE    iCard;
   MEMSZ    memSize;
   float    bscale;
   float    bzero;
   float    blankval;
   float    newBlankval;
   short int *  pNewData;

   fits_get_card_ival_(&bitpix, label_bitpix, pNHead, ppHead);

   /* Convert data to integer*2 if not already */
   if (bitpix != 16) {

      /* Get the scaling parameters from the header */
      if (fits_get_card_rval_(&bscale, (uchar *)Label_bscale, pNHead, ppHead)
       == FALSE_MWDUST) {
         bscale = 1.0;  /* Default value for BSCALE */
      }
      if (fits_get_card_rval_(&bzero , (uchar *)Label_bzero , pNHead, ppHead)
       == FALSE_MWDUST) {
         bzero = 0.0;  /* Default value for BZERO */
      }

      /* Allocate the minimum number of 2880-byte blocks for the data */
      memSize = ((int)((2*(*pNData)-1)/2880) + 1) * 2880;
      ccalloc_(&memSize, (void **)&pNewData);

      /* Convert the data and write to the new array */
      /* Note that nothing is done to rescale BLANK values properly */
      for (iData=0; iData < *pNData; iData++) {
         pNewData[iData] =
          fits_get_ival_(&iData, &bitpix, &bscale, &bzero, ppData);
      }

      /* Free the memory from the old array, and change the ppData pointer
         to point to the new array */
      ccfree_((void **)ppData);
      *ppData = (uchar *)pNewData;

      /* Change the BITPIX card to 16, indicating the data is integer*2 */
      newBitpix = 16;
      fits_change_card_ival_(&newBitpix, label_bitpix, pNHead, ppHead);

      /* Delete the BSCALE and BZERO cards which are no longer used */
      fits_delete_card_(label_bscale, pNHead, ppHead);
      fits_delete_card_(label_bzero , pNHead, ppHead);

      /* Rescale the BLANK card if it exists */
      if ((iCard = fits_find_card_(label_blank, pNHead, ppHead)) != *pNHead) {
         fits_get_card_rval_(&blankval, label_blank, pNHead, ppHead);
         if      (bitpix ==  8) newBlankval = blankval * bscale + bzero;
         else if (bitpix == 16) newBlankval = blankval * bscale + bzero;
         else if (bitpix == 32) newBlankval = blankval * bscale + bzero;
         else if (bitpix == -8) newBlankval = blankval;
         else if (bitpix ==-32) newBlankval = blankval;
         else if (bitpix ==-64) newBlankval = blankval;
         fits_change_card_rval_(&newBlankval, label_blank, pNHead, ppHead);
      }

   }
}

/******************************************************************************/
/*
 * Add a card immediately before the END card, or as the next card
 * (if no blank or END card), whichever comes first.  Return the card
 * number of the added card.
 * Memory is dynamically allocated if necessary by adding another 2880-byte
 * block.
 */
HSIZE fits_add_card_
  (uchar    pCard[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    numCardEnd;
   MEMSZ    memSize;
   uchar    pCardTemp[80];
   uchar *  pNewHeader;

   fits_string_to_card_(pCard, pCardTemp);

   numCardEnd=fits_find_card_(card_end, pNHead, ppHead);

   /* Test to see if more memory is needed for the header */
   if ((*pNHead)%36 == 0) {
      /* Copy header to new location and change appropriate pointers */
      memSize = (36+(*pNHead)) * 80;
      ccalloc_(&memSize, (void **)&pNewHeader);
      if (*pNHead > 0) {
         memmove(pNewHeader, *ppHead, (*pNHead)*80);
         ccfree_((void **)ppHead);
      }
      *ppHead = pNewHeader;
      numCardEnd += (pNewHeader - *ppHead);
   }

   if ((*pNHead > 0) && (numCardEnd<*pNHead) ) {
      /* Copy the end card forward 80 bytes in memory */
      memmove(&(*ppHead)[(numCardEnd+1)*80], &(*ppHead)[numCardEnd*80], 80);
      /* Add the new card where the END card had been */
      memmove(&(*ppHead)[numCardEnd*80], pCardTemp, 80);
      (*pNHead)++;
      return numCardEnd;
   }
   else {
      /* There is no end card, so simply add the new card at end of header */
      memmove(&(*ppHead)[(*pNHead)*80], pCardTemp, 80);
      return (*pNHead)++;
   }
}

/******************************************************************************/
/*
 * Add a card in the first card with a blank label or immediately before
 * the END card, or as the next card (if no blank or END card), whichever
 * comes first.  Return the card number of the added card.
 * Memory is dynamically allocated if necessary.
 */
HSIZE fits_add_cardblank_
  (uchar    pCard[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    numCardEmpty;
   HSIZE    numCardEnd;
   MEMSZ    memSize;
   uchar *  pHead = *ppHead;
   uchar *  pNewHeader;

   numCardEmpty = fits_find_card_(card_empty, pNHead, ppHead);
   numCardEnd   = fits_find_card_(card_end  , pNHead, ppHead);

   /* First case finds a blank card before the end card that is overwritten  */
   if ((*pNHead > 0) && (numCardEmpty < numCardEnd)) {
      memmove(&pHead[numCardEmpty*80], pCard, 80);
      return numCardEmpty;
   }
   else {
      /* Test to see if more memory is needed for the header */
      if ((*pNHead)%36 == 0) {
         /* Copy header to new location and change appropriate pointers */
         memSize = (36+(*pNHead)) * 80;
         ccalloc_(&memSize, (void **)&pNewHeader);
         memmove(pNewHeader, pHead, (*pNHead)*80);
         ccfree_((void **)&pHead);
         pHead = pNewHeader;
         numCardEmpty += (pNewHeader-pHead);
         numCardEnd   += (pNewHeader-pHead);
      }
      if ((*pNHead > 0) && (numCardEnd < *pNHead) ) {
         /* Copy the end card forward 80 bytes in memory */
         memmove(&pHead[(numCardEnd+1)*80], &pHead[numCardEnd*80], 80);
         /* Add the new card where the END card had been */
         memmove(&pHead[numCardEnd*80], pCard, 80);
         (*pNHead)++;
         return numCardEnd;
      } else {
         /* There is no end card, so simply add the new card at end of header */
         memmove(&pHead[(*pNHead)*80], pCard, 80);
         return (*pNHead)++;
      }
   }
}

/******************************************************************************/
/*
 * Create a new card with the given label and integer value.
 * Note that this can create multiple cards with the same label.
 * Return the card number of the new card.
 */
HSIZE fits_add_card_ival_
  (int   *  pIval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[81]; /* Last character is for null from sprintf() */

   sprintf((char *)pCardTemp, "%-8.8s= %20d", pLabel, *pIval);
   iCard = fits_add_card_(pCardTemp, pNHead, ppHead);

   return iCard;
}

/******************************************************************************/
/*
 * Create a new card with the given label and real value.
 * Note that this can create multiple cards with the same label.
 * Return the card number of the new card.
 */
HSIZE fits_add_card_rval_
  (float *  pRval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[81]; /* Last character is for null from sprintf() */

   sprintf((char *)pCardTemp, "%-8.8s= %20.7e", pLabel, *pRval);
   iCard = fits_add_card_(pCardTemp, pNHead, ppHead);

   return iCard;
}

/******************************************************************************/
/*
 * Create a new card with the given label and string value.
 * Note that this can create multiple cards with the same label.
 * Return the card number of the new card.
 */
HSIZE fits_add_card_string_
  (char  *  pStringVal,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[81]; /* Last character is for null from sprintf() */

   /* !!!??? NOTE: A QUOTE SHOULD BE WRITTEN AS 2 SINGLE QUOTES */
   sprintf((char *)pCardTemp, "%-8.8s= '%-1.68s'", pLabel, pStringVal);
   iCard = fits_add_card_(pCardTemp, pNHead, ppHead);

   return iCard;
}

/******************************************************************************/
/*
 * Create a new COMMENT card with the given string.
 * Note that this can create multiple cards with the same label.
 * Return the card number of the new card.
 */
HSIZE fits_add_card_comment_
  (char  *  pStringVal,
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[81]; /* Last character is for null from sprintf() */

   /* !!!??? NOTE: A QUOTE SHOULD BE WRITTEN AS 2 SINGLE QUOTES */
   sprintf((char *)pCardTemp, "COMMENT %-1.72s", pStringVal);
   iCard = fits_add_card_(pCardTemp, pNHead, ppHead);

   return iCard;
}

/******************************************************************************/
/*
 * Create a new HISTORY card with the given string.
 * Note that this can create multiple cards with the same label.
 * Return the card number of the new card.
 */
HSIZE fits_add_card_history_
  (char  *  pStringVal,
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[81]; /* Last character is for null from sprintf() */

   /* !!!??? NOTE: A QUOTE SHOULD BE WRITTEN AS 2 SINGLE QUOTES */
   sprintf((char *)pCardTemp, "HISTORY %-1.72s", pStringVal);
   iCard = fits_add_card_(pCardTemp, pNHead, ppHead);

   return iCard;
}

/******************************************************************************/
/*
 * Delete all cards where the label is blank.
 * Return the number of cards that were discarded.
 */
HSIZE fits_purge_blank_cards_
  (HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    numDelete;

   numDelete = 0;
   while (fits_delete_card_(label_empty, pNHead, ppHead) != *pNHead) {
      numDelete++;
   }

   return numDelete;
}

/******************************************************************************/
/*
 * Delete the first card that matches the given label, or do nothing if no
 * matches are found.  Return the card number of the deleted card,
 * or return nHead (out of range) if no match was found.
 */
HSIZE fits_delete_card_
  (uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   HSIZE    jCard;
   uchar *  pHead = *ppHead;

   iCard = fits_find_card_(pLabel, pNHead, ppHead);
   if (iCard < *pNHead) {
      (*pNHead)--;
      for (jCard=iCard; jCard <* pNHead; jCard++) {
         memmove(&pHead[jCard*80], &pHead[(jCard+1)*80], 80);
      }
      memmove(&pHead[jCard*80], card_empty, 80);
   }
   return iCard;
}

/******************************************************************************/
/*
 * Return the card number of the 1st header card with the label passed,
 * or return nHead (out of range) if no match was found.
 */
HSIZE fits_find_card_
  (uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar *  pHead;

   if (*pNHead == 0) iCard=0;
   else {
      pHead = *ppHead;
      for (iCard=0;
	   (iCard<*pNHead) && (strncmp((const char*)pLabel, (const char*)&pHead[iCard*80],8)!=0); iCard++);
   }
   return iCard;
}

/******************************************************************************/
/* Swap the integer values in the cards that match the passed labels.
 */
void fits_swap_cards_ival_
  (uchar *  pLabel1,
   uchar *  pLabel2,
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      ival1;
   int      ival2;
 
   fits_get_card_ival_(&ival1, pLabel1, pNHead, ppHead);
   fits_get_card_ival_(&ival2, pLabel2, pNHead, ppHead);
   fits_change_card_ival_(&ival2, pLabel1, pNHead, ppHead);
   fits_change_card_ival_(&ival1, pLabel2, pNHead, ppHead);
}
 
/******************************************************************************/
/* Swap the integer values in the cards that match the passed labels.
 */
void fits_swap_cards_rval_
  (uchar *  pLabel1,
   uchar *  pLabel2,
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   float    rval1;
   float    rval2;
 
   fits_get_card_rval_(&rval1, pLabel1, pNHead, ppHead);
   fits_get_card_rval_(&rval2, pLabel2, pNHead, ppHead);
   fits_change_card_rval_(&rval2, pLabel1, pNHead, ppHead);
   fits_change_card_rval_(&rval1, pLabel2, pNHead, ppHead);
}

/******************************************************************************/
/*
 * Find the 1st header card whose label matches the label passed,
 * and return the integer value of the argument after the label.
 * Return TRUE_MWDUST if there is a match, and FALSE_MWDUST if there is none.
 */
int fits_get_card_ival_
  (int   *  pIval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   HSIZE    iret;
   uchar *  pHead = *ppHead;
   char     pTemp[21];

   for (iCard=0;
	(iCard<*pNHead) && (strncmp((const char*)pLabel, (const char*)&pHead[iCard*80], 8)!=0); iCard++);
   if (iCard < *pNHead) {
#if 0
     sscanf(&pHead[iCard*80+10], "%20d", pIval);
#endif
     memmove(pTemp, &pHead[iCard*80+10], 20);
     pTemp[20] = '\0';
     sscanf(pTemp, "%d", pIval);
     iret = TRUE_MWDUST;
   }
   else {
     iret = FALSE_MWDUST;
   }
   return iret;
}

/******************************************************************************/
/*
 * Find the 1st header card whose label matches the label passed,
 * and return the real (float) value of the argument after the label.
 * Return TRUE_MWDUST if there is a match, and FALSE_MWDUST if there is none.
 */
int fits_get_card_rval_
  (float *  pRval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      iCard;
   int      iret;
   uchar *  pHead = *ppHead;
   char     pTemp[21];

   
   for (iCard=0; (iCard<*pNHead) && (strncmp((const char*)pLabel, (const char*)&pHead[iCard*80], 8)!=0);
    iCard++);
   if (iCard < *pNHead) {
#if 0
     sscanf(&pHead[iCard*80+10], "%20f", pRval);
#endif
     memmove(pTemp, &pHead[iCard*80+10], 20);
     pTemp[20] = '\0';
     sscanf(pTemp, "%f", pRval);
     iret = TRUE_MWDUST;
   }
   else {
     iret = FALSE_MWDUST;
   }
   return iret;
}

#if 0
/******************************************************************************/
/*
 * Return TRUE_MWDUST if there is a match, and FALSE_MWDUST if there is none.
 */
int fits_get_julian_date_
  (float *  pJulianDate,
   uchar    pLabelDate[],
   uchar    pLabelTime[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      iret;
   int      month;
   int      date;
   int      year;
   float    time;

   if (iret=fits_get_card_date_(month,date,year,pLabelDate,pNHead,ppHead)
    == TRUE_MWDUST) {
      *pJulianDate=...
      if (fits_get_card_time_(&time,pLabelTime,pNHead,ppHead) == TRUE_MWDUST) {
         *pJulianDate+=...
      }
   } else {
      *pJulianDate=0.0;
   }
   return iret;
}
#endif

/******************************************************************************/
/*
 * Find the 1st header card whose label matches the label passed,
 * and return the date as three integers month, date and year.
 * Return TRUE_MWDUST if there is a match, and FALSE_MWDUST if there is none.
 */
int fits_get_card_date_
  (int   *  pMonth,
   int   *  pDate,
   int   *  pYear,
   uchar    pLabelDate[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      iret;
   char  *  pStringVal;

   iret = fits_get_card_string_(&pStringVal, pLabelDate, pNHead, ppHead);
   if (iret == TRUE_MWDUST) {
      sscanf(pStringVal, "%d/%d/%d", pMonth, pDate, pYear);
      if (*pYear < 1900) *pYear += 1900;
      /* Free the memory used for the string value of this card */
      ccfree_((void **)&pStringVal);
   }
   return iret;
}

/******************************************************************************/
/*
 * Find the 1st header card whose label matches the label passed,
 * and return the time (TIME, RA, DEC or HA, for example) converted
 * to a real value.  Typically, this is used with TIME, RA or HA
 * to return a value in hours, or it is used with DEC to return a
 * value in degrees.
 * Return TRUE_MWDUST if there is a match, and FALSE_MWDUST if there is none.
 */
int fits_get_card_time_
  (float *  pTime,
   uchar    pLabelTime[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      iret;
   int      timeHour;
   int      timeMin;
   float    timeSec;
   char  *  pStringVal;

   iret = fits_get_card_string_(&pStringVal, pLabelTime, pNHead, ppHead);
   if (iret == TRUE_MWDUST) {
      sscanf(pStringVal, "%d:%d:%f", &timeHour, &timeMin, &timeSec);
      *pTime=abs(timeHour) + timeMin/60.0 + timeSec/3600.0;
      /* Make the returned value negative if a minus sign is in the string */
      if (strchr(pStringVal, '-') != NULL) *pTime=-(*pTime);
      /* Free the memory used for the string value of this card */
      ccfree_((void **)&pStringVal);
   } else {
      *pTime = 0.0;
   }
   return iret;
}

/******************************************************************************/
/*
 * Find the 1st header card whose label matches the label passed,
 * and return a pointer to the string argument after the label.
 * Memory is dynamically allocated for the string argument.
 * Return TRUE_MWDUST if there is a match, and FALSE_MWDUST if there is none.
 * If there is not match, then create and return the string pStringUnknown.
 */
int fits_get_card_string_
  (char  ** ppStringVal,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   int      iChar;
   int      iret;
   HSIZE    iCard;
   MEMSZ    memSize;
   uchar *  pHead = *ppHead;
   char  *  pTemp;
   char     pStringUnknown[] = "?";

   memSize = 70;
   ccalloc_(&memSize, (void **)&pTemp);
   for (iCard=0;
    (iCard<*pNHead) && (strncmp((const char*)pLabel, (const char*)&pHead[iCard*80], 8)!=0); iCard++);
   if (iCard < *pNHead) {
   /* It must start with a single quote in column 11 (1-indexed) if not blank.
      Otherwise, an empty string is returned, which is O.K. */
     iChar = 11;
     /* Copy characters from column 12 until column 80 or another single
        quote is reached. */
     /* !!!??? NOTE: TWO SINGLE QUOTES SHOULD BE READ IN AS A QUOTE */
     if (pHead[iCard*80+10]=='\'') {
       while (iChar<80 && (pTemp[iChar-11]=pHead[iCard*80+iChar]) != '\'')
        iChar++;
     }

     pTemp[iChar-11]='\0';  /* Pad with a NULL at the end of the string */
     /* Remove trailing blanks; leading blanks are significant */
     iChar = strlen(pTemp);
     while (iChar>0 && pTemp[--iChar]==' ') pTemp[iChar]='\0';

     iret = TRUE_MWDUST;
   }
   else {
     strcpy(pTemp, pStringUnknown);

     iret = FALSE_MWDUST;
   }

   *ppStringVal=pTemp;
   return iret;
}

/******************************************************************************/
/*
 * Change the 1st card that matches the passed label, or add a card if there
 * is not a match.  Return the card number of the changed or added card.
 */
HSIZE fits_change_card_
  (uchar    pCard[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[80];
   uchar *  pHead = *ppHead;

   fits_string_to_card_(pCard, pCardTemp);

   iCard = fits_find_card_(pCardTemp, pNHead, ppHead);
   if (iCard < *pNHead) {
      memmove(&pHead[iCard*80], pCardTemp, 80);
   } else {
      iCard = fits_add_card_(pCardTemp, pNHead, ppHead);
   }

   return iCard;
}

/******************************************************************************/
/*
 * Find the 1st header card whose label matches the label passed,
 * and change the integer value of the argument after the label.
 * If no card exists, then create one.  Return the card number of
 * the changed or added card.
 */
HSIZE fits_change_card_ival_
  (int   *  pIval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[81]; /* Last character is for null from sprintf() */

   sprintf( (char*)pCardTemp, "%-8.8s= %20d", pLabel, *pIval);
   iCard = fits_change_card_(pCardTemp, pNHead, ppHead);

   return iCard;
}

/******************************************************************************/
/*
 * Find the 1st header card whose label matches the label passed,
 * and change the real (float) value of the argument after the label.
 * If no card exists, then create one.  Return the card number of
 * the changed or added card.
 */
HSIZE fits_change_card_rval_
  (float *  pRval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[81]; /* Last character is for null from sprintf() */

   sprintf( (char*)pCardTemp, "%-8.8s= %20.7e", pLabel, *pRval);
   iCard = fits_change_card_(pCardTemp, pNHead, ppHead);

   return iCard;
}

/******************************************************************************/
/*
 * Find the 1st header card whose label matches the label passed,
 * and change the string value of the argument after the label.
 * If no card exists, then create one.  Return the card number of
 * the changed or added card.
 */
HSIZE fits_change_card_string_
  (char  *  pStringVal,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead)
{
   HSIZE    iCard;
   uchar    pCardTemp[81]; /* Last character is for null from sprintf() */

   /* !!!??? NOTE: A QUOTE SHOULD BE WRITTEN AS 2 SINGLE QUOTES */
   sprintf( (char*)pCardTemp, "%-8.8s= '%-1.68s'", pLabel, pStringVal);
   iCard = fits_change_card_(pCardTemp, pNHead, ppHead);

   return iCard;
}

/******************************************************************************/
/* Convert a character string to a FITS-complient 80-character card.
 * The string is copied until either a NULL or CR is reached or the 80th
 * character is reached.  The remainder of the card is padded with spaces.
 * In addition, the label part of the card (the first 8 characters)
 * are converted to upper case.
 *
 * Note that pCard[] must be dimensioned to at least the necessary 80
 * characters.
 */
void fits_string_to_card_
  (uchar    pString[],
   uchar    pCard[])
{
   int      iChar;
   int      qNull;

   /* Copy the raw string into the card array */
   memmove(pCard, pString, 80);

   /* Search for a NULL or CR in the card, and replace that character and
    * all following characters with a space.
    */
   qNull = FALSE_MWDUST;
   iChar = 0;
   while (iChar < 80) {
      if (pCard[iChar] == '\0' || pCard[iChar] == '\n') qNull = TRUE_MWDUST;
      if (qNull == TRUE_MWDUST) pCard[iChar] = ' ';
      iChar++;
   }

   /* Convert the label (the first 8 characters) to upper case) */
   for (iChar=0; iChar < 8; iChar++) {
      pCard[iChar] = toupper(pCard[iChar]);
   }
}

/******************************************************************************/
/*
 * Return the (float) value of the data array indexed by the iloc'th elements,
 * taking care to use the proper data format as specified by bitpix.
 * Several unconventional values for bitpix are supported: 32, 8, -8.
 * For a 2-dimensional array, set iloc=x+y*naxis1.
 */
float fits_get_rval_
  (DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBscale,
   float *  pBzero,
   uchar ** ppData)
{
   float    rval;
   uchar     * pIdata8  = (uchar     *)(*ppData);
   short int * pIdata16 = (short int *)(*ppData);
   long  int * pIdata32 = (long  int *)(*ppData);
   float     * pRdata32 = (float     *)(*ppData);
   double    * pRdata64 = (double    *)(*ppData);

   if      (*pBitpix ==-32) rval = pRdata32[*pIloc];
   else if (*pBitpix == 16) rval = pIdata16[*pIloc] * (*pBscale) + (*pBzero);
   else if (*pBitpix == 32) rval = pIdata32[*pIloc] * (*pBscale) + (*pBzero);
   else if (*pBitpix ==-64) rval = pRdata64[*pIloc];
   else if (*pBitpix ==  8) rval = pIdata8 [*pIloc] * (*pBscale) + (*pBzero);
   else if (*pBitpix == -8) rval = pIdata8 [*pIloc];
   else                     rval = 0.0; /* Invalid BITPIX! */
   return rval;
}

/******************************************************************************/
/*
 * Return the (int) value of the data array indexed by the iloc'th elements,
 * taking care to use the proper data format as specified by bitpix.
 * Several unconventional values for bitpix are supported: 32, 8, -8.
 * For a 2-dimensional array, set iloc=x+y*naxis1.
 */
int fits_get_ival_
  (DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBscale,
   float *  pBzero,
   uchar ** ppData)
{
   int      ival;
   float    rval;
   uchar     * pIdata8  = (uchar     *)(*ppData);
   short int * pIdata16 = (short int *)(*ppData);
   long  int * pIdata32 = (long  int *)(*ppData);
   float     * pRdata32 = (float     *)(*ppData);
   double    * pRdata64 = (double    *)(*ppData);

   if      (*pBitpix ==-32) rval = pRdata32[*pIloc];
   else if (*pBitpix == 16) rval = pIdata16[*pIloc] * (*pBscale) + (*pBzero);
   else if (*pBitpix == 32) rval = pIdata32[*pIloc] * (*pBscale) + (*pBzero);
   else if (*pBitpix ==-64) rval = pRdata64[*pIloc];
   else if (*pBitpix ==  8) rval = pIdata8 [*pIloc] * (*pBscale) + (*pBzero);
   else if (*pBitpix == -8) rval = pIdata8 [*pIloc];
   else                     rval = 0.0; /* Invalid BITPIX! */

   /* Round to the nearest integer */
   if (rval >= 0.0) {
     ival = (int)(rval + 0.5);
   } else {
     ival = (int)(rval - 0.5);
   }

   return ival;
}

/******************************************************************************/
/*
 * Put a (float) value into location iloc of the data array, taking care to
 * use the proper data format as specified by bitpix.  For a 2-dimensional
 * array, set iloc=x+y*naxis1.
 * Several unconventional values for bitpix are supported: 32, 8, -8.
 * Note: Is the rounding done properly!!!???
 */
void fits_put_rval_
  (float *  pRval,
   DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBscale,
   float *  pBzero,
   uchar ** ppData)
{
   uchar     * pIdata8  = (uchar     *)(*ppData);
   short int * pIdata16 = (short int *)(*ppData);
   long  int * pIdata32 = (long  int *)(*ppData);
   float     * pRdata32 = (float     *)(*ppData);
   double    * pRdata64 = (double    *)(*ppData);

   if      (*pBitpix ==-32) 
     { pRdata32[*pIloc] = *pRval; }
   else if (*pBitpix == 16) 
     { pIdata16[*pIloc] = (short)( (*pRval - *pBzero) / (*pBscale) ); }
   else if (*pBitpix == 32) 
     { pIdata32[*pIloc] = (long)((*pRval - *pBzero) / (*pBscale)); }
   else if (*pBitpix ==-64) 
     { pRdata64[*pIloc] = *pRval; }
   else if (*pBitpix ==  8) 
     { pIdata8 [*pIloc] = (uchar)((*pRval - *pBzero) / (*pBscale)); }
   else if (*pBitpix == -8) 
     { pIdata8 [*pIloc] = (uchar)*pRval; }
}

/******************************************************************************/
/*
 * Ask whether a particular pixel position in the data array is equal
 * to the value specified by the BLANK card.  This test is performed WITHOUT
 * first rescaling the data.  Pass the blank value as
 * a real variable, even if it was originally integer.  For a 2-dimensional
 * array, set iloc=x+y*naxis1.  Return TRUE_MWDUST (iq!=0) if the pixel is
 * BLANK, or FALSE_MWDUST (iq==0) if it is not.
 * The value blankval must be found first and passed to this routine.
 * Several unconventional values for bitpix are supported: 32, 8, -8.
 */
int fits_qblankval_
  (DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBlankval,
   uchar ** ppData)
{
   int      iq;
   uchar     * pIdata8  = (uchar     *)(*ppData);
   short int * pIdata16 = (short int *)(*ppData);
   long  int * pIdata32 = (long  int *)(*ppData);
   float     * pRdata32 = (float     *)(*ppData);
   double    * pRdata64 = (double    *)(*ppData);

   if      (*pBitpix ==-32) iq = ( pRdata32[*pIloc] == (*pBlankval) );
   else if (*pBitpix == 16) iq = ( pIdata16[*pIloc] == (*pBlankval) );
   else if (*pBitpix == 32) iq = ( pIdata32[*pIloc] == (*pBlankval) );
   else if (*pBitpix ==-64) iq = ( pRdata64[*pIloc] == (*pBlankval) );
   else if (*pBitpix ==  8) iq = ( pIdata8 [*pIloc] == (*pBlankval) );
   else if (*pBitpix == -8) iq = ( pIdata8 [*pIloc] == (*pBlankval) );
   else                     iq = FALSE_MWDUST; /* Invalid BITPIX! */

   return iq;
}

/******************************************************************************/
/*
 * Replace a data element by a BLANK value as determined by blankval.
 * The value is assigned WITHOUT rescaling.
 * The value blankval must be found first and passed to this routine.
 * Several unconventional values for bitpix are supported: 32, 8, -8.
 */
void fits_put_blankval_
  (DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBlankval,
   uchar ** ppData)
{
   uchar     * pIdata8  = (uchar     *)(*ppData);
   short int * pIdata16 = (short int *)(*ppData);
   long  int * pIdata32 = (long  int *)(*ppData);
   float     * pRdata32 = (float     *)(*ppData);
   double    * pRdata64 = (double    *)(*ppData);

   if      (*pBitpix ==-32) pRdata32[*pIloc] = *pBlankval;
   else if (*pBitpix == 16) pIdata16[*pIloc] = (short)*pBlankval;
   else if (*pBitpix == 32) pIdata32[*pIloc] = (long)*pBlankval;
   else if (*pBitpix ==-64) pRdata64[*pIloc] = *pBlankval;
   else if (*pBitpix ==  8) pIdata8 [*pIloc] = (uchar)*pBlankval;
   else if (*pBitpix == -8) pIdata8 [*pIloc] = (uchar)*pBlankval;
}

/******************************************************************************/
/*
 * Replace any nulls in a card by spaces.
 */
void fits_purge_nulls
  (uchar    pCard[])
{
   int iChar;

   for (iChar=0; iChar < 80; iChar++) {
      if (pCard[iChar] == '\0') pCard[iChar] = ' ';
   }
}

/******************************************************************************/
/*
 * Read the next 80-character card from the specified device.
 * Return 0 if the END card is reached.
 */
int fits_get_next_card_
  (int   *  pFilenum,
   uchar    pCard[])
{
   int      iChar;

   for (iChar=0; iChar < 80; iChar++) {
      pCard[iChar] = fgetc(pFILEfits[*pFilenum]);
   }
   return strncmp((const char *)card_end, (const char *)pCard, 8);
}

/******************************************************************************/
/*
 * Write passed card to open file.  Return FALSE_MWDUST for a write error.
 */
int fits_put_next_card_
  (int   *  pFilenum,
   uchar    pCard[])
{
   int      iChar;
   int      retval;

   retval = TRUE_MWDUST;
   for (iChar=0; iChar < 80; iChar++) {
      if (fputc(pCard[iChar], pFILEfits[*pFilenum]) == EOF) retval = FALSE_MWDUST;
   }
   return retval;
}

/******************************************************************************/
/*
 * Determine the size of an individual datum based upon the FITS definitions
 * of the BITPIX card.
 */
int fits_size_from_bitpix_
  (int *pBitpix)
{
   int size;

   if      (*pBitpix ==   8) size = 1;
   else if (*pBitpix ==  16) size = 2;
   else if (*pBitpix ==  32) size = 4;
   else if (*pBitpix ==  64) size = 8;
   else if (*pBitpix == -16) size = 2;
   else if (*pBitpix == -32) size = 4;
   else if (*pBitpix == -64) size = 8;
   else                      size = 0; /* Bitpix undefined! */

   return size;
}

/******************************************************************************/
/*
 * For data of arbitrary dimensions, shift the pixels along the "*pSAxis"
 * axis by "*pShift" pixels, wrapping data around the image boundaries.
 * For example, if the middle dimension of a 3-dimen array is shifted:
 *   new[i,j+shift,k] = old[i,j,k]
 * The data can be of any data type (i.e., any BITPIX).
 */
void fits_pixshift_wrap_
  (int   *  pSAxis,
   DSIZE *  pShift,
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   uchar ** ppData)
{
   int      size;
   DSIZE    posShift;
   int      iAxis;
   int      numAxes;
   DSIZE *  pNaxis;
   DSIZE    dimBig;
   DSIZE    dimSml;
   DSIZE    indxBig;
   DSIZE    indxSml;
   DSIZE    offset;
   DSIZE    iloc;
   int      bitpix;
   MEMSZ    memSize;
   DSIZE    nVector;
   DSIZE    iVector;
   uchar *  pVector;

   fits_compute_axes_(pNHead, ppHead, &numAxes, &pNaxis);
   nVector = pNaxis[*pSAxis];
   posShift = *pShift;
   while (posShift < 0) posShift += nVector; /* Must be positive value */

   /* Allocate an array equal in size to one vector in the *pSAxis dimension */
   fits_get_card_ival_(&bitpix, label_bitpix, pNHead, ppHead);
   size = fits_size_from_bitpix_(&bitpix);
   memSize = size * nVector;
   ccalloc_(&memSize, (void **)&pVector);

   /* Compute the number of larger and smaller indices */
   dimBig = 1;
   for (iAxis=0; iAxis < *pSAxis; iAxis++) dimBig *= pNaxis[iAxis];
   dimSml = 1;
   for (iAxis=*pSAxis+1; iAxis < numAxes; iAxis++) dimSml *= pNaxis[iAxis];

   /* Loop through each of the larger and smaller indices */
   for (indxBig=0; indxBig < dimBig; indxBig++) {
   for (indxSml=0; indxSml < dimSml; indxSml++) {
      offset = indxBig * nVector * dimSml + indxSml;

      /* Copy vector into temporary vector */
      for (iVector=0; iVector < nVector; iVector++) {
         iloc = offset + iVector * dimSml;
         memmove(&pVector[iVector*size], &(*ppData)[iloc*size], size);
      }

      /* Copy the shifted vector back into the main data array */
      for (iVector=0; iVector < nVector; iVector++) {
         /* Use the MOD operator below to wrap the dimensions */
         iloc = offset + ((iVector+(posShift)) % nVector) * dimSml;
         memmove(&(*ppData)[iloc*size], &pVector[iVector*size], size);
      }
   } }

   /* Free memory */
   ccfree_((void **)&pVector);

   /* Plug a memory leak - D. Schlegel 06-Feb-1999 */
   fits_free_axes_(&numAxes, &pNaxis);
}

/******************************************************************************/
/*
 * For data of 2 dimensions, transpose the data by setting
 * pData[i][j] = pData[j][i].
 * A new data array is created and the old one is destroyed.
 * Also, swap the appropriate header cards.
 */
void fits_transpose_data_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   uchar ** ppData)
{
   int      bitpix;
   int      size;
   int      iByte;
   int      numAxes;
   DSIZE *  pNaxis;
   DSIZE    nData;
   DSIZE    iRow;
   DSIZE    iCol;
   DSIZE    ilocOld;
   DSIZE    ilocNew;
   MEMSZ    memSize;
   uchar *  pNewData;

   fits_compute_axes_(pNHead, ppHead, &numAxes, &pNaxis);
   if (numAxes == 2) {
      /* Allocate an array equal in size to the data array */
      nData = fits_compute_ndata_(pNHead, ppHead);
      fits_get_card_ival_(&bitpix, label_bitpix, pNHead, ppHead);
      size = fits_size_from_bitpix_(&bitpix);
      memSize = size * nData;
      ccalloc_(&memSize, (void **)&pNewData);

      /* Copy the data into the new data array, transposing the first 2 axes */
      for (iRow=0; iRow < pNaxis[1]; iRow++) {
         for (iCol=0; iCol < pNaxis[0]; iCol++) {
            ilocOld = size * (iRow*pNaxis[0] + iCol);
            ilocNew = size * (iCol*pNaxis[1] + iRow);
            /* For each data element, copy the proper number of bytes */
            for (iByte=0; iByte < size; iByte++) {
               pNewData[ilocNew+iByte] = (*ppData)[ilocOld+iByte];
            }
         }
      }

      /* Discard the old data array and return the new one */
      ccfree_((void **)ppData);
      *ppData = pNewData;

      /* Switch the values in the header of NAXIS1 and NAXIS2,
       * and the cards used to label the pixel numbers on those axes.
       */
      fits_swap_cards_ival_(label_naxis1, label_naxis2, pNHead, ppHead);
      fits_swap_cards_rval_(label_crpix1, label_crpix2, pNHead, ppHead);
      fits_swap_cards_rval_(label_crval1, label_crval2, pNHead, ppHead);
      fits_swap_cards_rval_(label_cdelt1, label_cdelt2, pNHead, ppHead);
   }

   /* Plug a memory leak - D. Schlegel 06-Feb-1999 */
   fits_free_axes_(&numAxes, &pNaxis);
}

/******************************************************************************/
/*
 * Average several rows (or columns) of a 2-dimensional array of floats.
 * If (*iq)==0 then average rows; if (*iq)==1 then average columns.
 * Note that *pNaxis1=number of columns and *pNaxis2=number of rows.
 * Memory is dynamically allocated for the output vector.
 */
void fits_ave_rows_r4_
  (int   *  iq,
   DSIZE *  pRowStart,
   DSIZE *  pNumRowAve,
   DSIZE *  pNaxis1,
   DSIZE *  pNaxis2,
   float ** ppData,
   float ** ppOut)
{
   DSIZE   iCol;
   DSIZE   iRow;
   DSIZE   rowStart;
   DSIZE   rowEnd;
   MEMSZ   memSize;
   float   weight;
   float * pData = *ppData;
   float * pOut;

   if (*iq == 0) {
      memSize = sizeof(float) * (*pNaxis1);
      ccalloc_(&memSize, (void **)ppOut);
      pOut = *ppOut;
      rowStart = max(0, *pRowStart);
      rowEnd = min(*pRowStart + *pNumRowAve, *pNaxis2);
      weight = (rowEnd + 1 - rowStart);
      for (iCol=0; iCol < *pNaxis1; iCol++) {
         pOut[iCol] = 0.0;
         for (iRow=rowStart; iRow <= rowEnd; iRow++) {
            pOut[iCol] += pData[iRow*(*pNaxis1) + iCol];
         }
         pOut[iCol] /= weight;
      }
   } else if (*iq == 1) {
      memSize = sizeof(float) * (*pNaxis2);
      ccalloc_(&memSize, (void **)ppOut);
      pOut = *ppOut;
      rowStart = max(0, *pRowStart);
      rowEnd = min(*pRowStart + *pNumRowAve, *pNaxis1);
      weight = (rowEnd + 1 - rowStart);
      for (iRow=0; iRow < *pNaxis2; iRow++) {
         pOut[iRow] = 0.0;
         for (iCol=rowStart; iCol <= rowEnd; iCol++) {
            pOut[iRow] += pData[iRow*(*pNaxis1) + iCol];
         }
         pOut[iRow] /= weight;
      }
   }

}

/******************************************************************************/
/*
 * Average several rows (or columns) of a 2-dimensional array of floats
 * with their standard deviations.  For each combined set of points:
 *    obj_ave = SUM_i {obj_i / sig_i^2} / SUM_i {1 / sig_i^2}
 *    sig_ave = 1 / SUM_i {1 / sig_i^2}
 * If (*iq)==0 then average rows; if (*iq)==1 then average columns.
 * Note that *pNaxis1=number of columns and *pNaxis2=number of rows.
 * Memory is dynamically allocated for the output vector.
 */
void fits_ave_obj_and_sigma_rows_r4_
  (int   *  iq,
   DSIZE *  pRowStart,
   DSIZE *  pNumRowAve,
   DSIZE *  pNaxis1,
   DSIZE *  pNaxis2,
   float ** ppObjData,
   float ** ppSigData,
   float ** ppObjOut,
   float ** ppSigOut)
{
   DSIZE   iCol;
   DSIZE   iRow;
   DSIZE   rowStart;
   DSIZE   rowEnd;
   DSIZE   iloc;
   MEMSZ   memSize;
   float   weight;
   float   oneOverSumVar;
   float * pObjData = *ppObjData;
   float * pSigData = *ppSigData;
   float * pObjOut;
   float * pSigOut;

   if (*iq == 0) {
      memSize = sizeof(float) * (*pNaxis1);
      ccalloc_(&memSize, (void **)ppObjOut);
      ccalloc_(&memSize, (void **)ppSigOut);
      pObjOut = *ppObjOut;
      pSigOut = *ppSigOut;
      rowStart = max(0, (*pRowStart));
      rowEnd = min((*pRowStart) + (*pNumRowAve) - 1, (*pNaxis2) - 1);
      for (iCol=0; iCol < *pNaxis1; iCol++) {
         pObjOut[iCol] = 0.0;
         oneOverSumVar = 0.0;
         for (iRow=rowStart; iRow <= rowEnd; iRow++) {
            iloc = iRow*(*pNaxis1) + iCol;
            weight = 1.0 / (pSigData[iloc] * pSigData[iloc]);
            pObjOut[iCol] += pObjData[iloc] * weight;
            oneOverSumVar += weight;
         }
         pObjOut[iCol] /= oneOverSumVar;
         pSigOut[iCol] = 1.0 / sqrt(oneOverSumVar);
      }

   } else if (*iq == 1) {
      memSize = sizeof(float) * (*pNaxis2);
      ccalloc_(&memSize, (void **)ppObjOut);
      ccalloc_(&memSize, (void **)ppSigOut);
      pObjOut = *ppObjOut;
      pSigOut = *ppSigOut;
      rowStart = max(0, (*pRowStart));
      rowEnd = min((*pRowStart) + (*pNumRowAve) - 1, (*pNaxis1) - 1);
      for (iRow=0; iRow < *pNaxis2; iRow++) {
         pObjOut[iRow] = 0.0;
         oneOverSumVar = 0.0;
         for (iCol=rowStart; iCol <= rowEnd; iCol++) {
            iloc = iRow*(*pNaxis1) + iCol;
            weight = 1.0 / (pSigData[iloc] * pSigData[iloc]);
            pObjOut[iRow] += pObjData[iloc] * weight;
            oneOverSumVar += weight;
         }
         pObjOut[iRow] /= oneOverSumVar;
         pSigOut[iRow] = 1.0 / sqrt(oneOverSumVar);
      }
   }

}

/******************************************************************************/
/*
 * Swap bytes between big-endian and little-endian.
 */
void fits_byteswap
  (int      bitpix,
   DSIZE    nData,
   uchar *  pData)
{
   int      ibits;
   DSIZE    idata;

   ibits = abs(bitpix);
   if (ibits == 16) {
      for (idata=0; idata < nData; idata++) {
         fits_bswap2( &pData[2*idata  ], &pData[2*idata+1] );
      }
   } else if (ibits == 32) {
      for (idata=0; idata < nData; idata++) {
         fits_bswap2( &pData[4*idata  ], &pData[4*idata+3] );
         fits_bswap2( &pData[4*idata+1], &pData[4*idata+2] );
      }
   } else if (ibits == 64) {
      for (idata=0; idata < nData; idata++) {
         fits_bswap2( &pData[8*idata  ], &pData[8*idata+7] );
         fits_bswap2( &pData[8*idata+1], &pData[8*idata+6] );
         fits_bswap2( &pData[8*idata+2], &pData[8*idata+5] );
         fits_bswap2( &pData[8*idata+3], &pData[8*idata+4] );
      }
   }

}

void fits_bswap2
  (uchar *  pc1,
   uchar *  pc2)
{
   uchar    ct;
   ct = *pc1;
   *pc1 = *pc2;
   *pc2 = ct;
}


#ifdef OLD_SUNCC

/******************************************************************************/
/*
 * Copy one section of memory (a string) to another, even if they overlap.
 * The advertised C library routine by this name does not actually exist
 * in old SunOS.
 */
void memmove
  (void  *  s,
   const void  *  ct,
   MEMSZ    n)
{
   MEMSZ    i;
   char  *  ps = (char *)s;
   const char  *  pct = (const char *)ct;
 
   /* Do nothing if ps == pct */
   if (ps > pct) for (i=0; i < n; i++) *(ps+n-i-1) = *(pct+n-i-1);
   else if (ps < pct) for (i=0; i < n; i++) *(ps+i) = *(pct+i);
}

#endif

/******************************************************************************/
/*
 * Change the size of memory for data, and return the new address as *ppData.
 * Copy contents of memory in common.
 */
void ccalloc_resize_
  (MEMSZ *  pOldMemSize,
   MEMSZ *  pNewMemSize,
   void  ** ppData)
{
   void  *  pNewData;

   if (*pNewMemSize > *pOldMemSize) {
      ccalloc_(pNewMemSize,&pNewData);
      memmove((void *)pNewData,(void *)(*ppData),*pOldMemSize);
      ccfree_(ppData);
      *ppData = pNewData;
   } else if (*pNewMemSize < *pOldMemSize) {
      ccalloc_(pNewMemSize,&pNewData);
      memmove((void *)pNewData,(void *)(*ppData),*pNewMemSize);
      ccfree_(ppData);
      *ppData = pNewData;
   }
}

/******************************************************************************/
/*
 * Allocate *pMemSize bytes of data.  The starting memory location is *ppData.
 * If the array has previously been allocated, then resize it.
 */
void ccrealloc_
  (MEMSZ *  pMemSize,
   void  ** ppData)
{
   if (*ppData == NULL) {
      *ppData = (void *)malloc((size_t)(*pMemSize));
   } else {
      *ppData = (void *)realloc(*ppData,(size_t)(*pMemSize));
   }
}

/******************************************************************************/
/*
 * Allocate *pMemSize bytes of data.  The starting memory location is *ppData.
 * Also zero all of the data byes.
 */
void ccalloc_init
  (MEMSZ *  pMemSize,
   void  ** ppData)
{
   size_t   nobj = 1;
   *ppData = (void *)calloc(nobj, (size_t)(*pMemSize));
}

/******************************************************************************/
/*
 * Allocate *pMemSize bytes of data.  The starting memory location is *ppData.
 */
void ccalloc_
  (MEMSZ *  pMemSize,
   void  ** ppData)
{
   *ppData = (void *)malloc((size_t)(*pMemSize));
}

/******************************************************************************/
/*
 * Free the memory block that starts at address *ppData.
 */
void ccfree_
  (void  ** ppData)
{
   free(*ppData);
   *ppData = NULL;
}

/******************************************************************************/
float * ccvector_build_
  (MEMSZ    n)
{
   float * pVector = (float *)malloc((size_t)(sizeof(float) * n));
   return pVector;
}

/******************************************************************************/
double * ccdvector_build_
  (MEMSZ    n)
{
   double * pVector = (double *)malloc((size_t)(sizeof(double) * n));
   return pVector;
}

/******************************************************************************/
int * ccivector_build_
  (MEMSZ    n)
{
   int * pVector = (int *)malloc((size_t)(sizeof(int) * n));
   return pVector;
}
 
/******************************************************************************/
float ** ccpvector_build_
  (MEMSZ    n)
{
   float ** ppVector = (float **)malloc((size_t)(sizeof(float *) * n));
   return ppVector;
}
 
/******************************************************************************/
/* Build a vector of pointers to arrays of type (float **) */
float *** ccppvector_build_
  (MEMSZ    n)
{
   float *** pppVector = (float ***)malloc((size_t)(sizeof(float **) * n));
   return pppVector;
}
 
/******************************************************************************/
float * ccvector_rebuild_
  (MEMSZ    n,
   float *  pOldVector)
{
   float * pVector;

   if (pOldVector == NULL) {
      pVector = (float *)malloc((size_t)(sizeof(float) * n));
   } else {
      pVector = (float *)realloc(pOldVector,(size_t)(sizeof(float) * n));
   }

   return pVector;
}

/******************************************************************************/
double * ccdvector_rebuild_
  (MEMSZ    n,
   double *  pOldVector)
{
   double * pVector;

   if (pOldVector == NULL) {
      pVector = (double *)malloc((size_t)(sizeof(double) * n));
   } else {
      pVector = (double *)realloc(pOldVector,(size_t)(sizeof(double) * n));
   }

   return pVector;
}

/******************************************************************************/
int * ccivector_rebuild_
  (MEMSZ    n,
   int   *  pOldVector)
{
   int   *  pVector;

   if (pOldVector == NULL) {
      pVector = (int *)malloc((size_t)(sizeof(int) * n));
   } else {
      pVector = (int *)realloc(pOldVector,(size_t)(sizeof(int) * n));
   }

   return pVector;
}

/******************************************************************************/
float ** ccpvector_rebuild_
  (MEMSZ    n,
   float ** ppOldVector)
{
   float ** ppVector;

   if (ppOldVector == NULL) {
      ppVector = (float **)malloc((size_t)(sizeof(float *) * n));
   } else {
      ppVector = (float **)realloc(ppOldVector,(size_t)(sizeof(float *) * n));
   }

   return ppVector;
}

/******************************************************************************/
/* Build a vector of pointers to arrays of type (float **) */
float *** ccppvector_rebuild_
  (MEMSZ    n,
   float *** pppOldVector)
{
   float *** pppVector;
 
   if (pppOldVector == NULL) {
      pppVector = (float ***)malloc((size_t)(sizeof(float **) * n));
   } else {
      pppVector = (float ***)realloc(pppOldVector,(size_t)(sizeof(float **) * n)
);
   }
 
   return pppVector;
}

/******************************************************************************/
void ccvector_free_
  (float *  pVector)
{
   free((void *)pVector);
}

/******************************************************************************/
void ccdvector_free_
  (double *  pVector)
{
   free((void *)pVector);
}

/******************************************************************************/
void ccivector_free_
  (int   *  pVector)
{
   free((void *)pVector);
}

/******************************************************************************/
void ccpvector_free_
  (float ** ppVector)
{
   free((void *)ppVector);
}

/******************************************************************************/
void ccppvector_free_
  (float *** pppVector)
{
   free((void *)pppVector);
}

/******************************************************************************/
/* Build an nRow x nCol matrix, in pointer-style.
 * Allocate one contiguous array of floats with  nRow*nCol elements.
 * Then create a set of nRow pointers, each of which points to the next nCol
 *  elements.
 */
float ** ccarray_build_
  (MEMSZ    nRow,
   MEMSZ    nCol)
{
   MEMSZ    iRow;
 
   float *  pSpace  = (float *)malloc(sizeof(float ) * nRow * nCol);
   float ** ppArray = (float**)malloc(sizeof(float*) * nRow);

   for (iRow = 0; iRow < nRow; iRow++) {
      /* Quantity (iRow*nCol) scales by sizeof(float) */
      ppArray[iRow] = pSpace + (iRow * nCol);
   }

   return ppArray;
}

/******************************************************************************/
/* Build an nRow x nCol matrix, in pointer-style.
 * Allocate one contiguous array of doubles with  nRow*nCol elements.
 * Then create a set of nRow pointers, each of which points to the next nCol
 *  elements.
 */
double ** ccdarray_build_
  (MEMSZ    nRow,
   MEMSZ    nCol)
{
   MEMSZ    iRow;
 
   double *  pSpace  = (double *)malloc(sizeof(double ) * nRow * nCol);
   double ** ppArray = (double**)malloc(sizeof(double*) * nRow);

   for (iRow = 0; iRow < nRow; iRow++) {
      /* Quantity (iRow*nCol) scales by sizeof(double) */
      ppArray[iRow] = pSpace + (iRow * nCol);
   }

   return ppArray;
}

/******************************************************************************/
/* Build an nRow x nCol matrix, in pointer-style.
 * Allocate one contiguous array of ints with  nRow*nCol elements.
 * Then create a set of nRow pointers, each of which points to the next nCol
 *  elements.
 */
int ** cciarray_build_
  (MEMSZ    nRow,
   MEMSZ    nCol)
{
   MEMSZ    iRow;
 
   int *  pSpace  = (int *)malloc(sizeof(int ) * nRow * nCol);
   int ** ppArray = (int**)malloc(sizeof(int*) * nRow);

   for (iRow = 0; iRow < nRow; iRow++) {
      /* Quantity (iRow*nCol) scales by sizeof(int) */
      ppArray[iRow] = pSpace + (iRow * nCol);
   }

   return ppArray;
}

/******************************************************************************/
/* Build an nRow x nCol matrix, in pointer-style.
 * Allocate one contiguous array of floats with  nRow*nCol elements.
 * Then create a set of nRow pointers, each of which points to the next nCol
 *  elements.
 */
float ** ccarray_rebuild_
  (MEMSZ    nRow,
   MEMSZ    nCol,
   float ** ppOldArray)
{
   MEMSZ    iRow;
 
   float *  pSpace;
   float ** ppArray;

   if (ppOldArray == NULL) {
      pSpace  = (float *)malloc(sizeof(float ) * nRow * nCol);
      ppArray = (float**)malloc(sizeof(float*) * nRow);
   } else {
      ppArray = (float**)realloc(ppOldArray, sizeof(float*) * nRow);
      pSpace  = (float *)realloc(ppOldArray[0],sizeof(float ) * nRow * nCol);
   }

   for (iRow = 0; iRow < nRow; iRow++) {
      /* Quantity (iRow*nCol) scales by sizeof(float) */
      ppArray[iRow] = pSpace + (iRow * nCol);
   }

   return ppArray;
}

/******************************************************************************/
/* Build an nRow x nCol matrix, in pointer-style.
 * Allocate one contiguous array of doubles with  nRow*nCol elements.
 * Then create a set of nRow pointers, each of which points to the next nCol
 *  elements.
 */
double ** ccdarray_rebuild_
  (MEMSZ    nRow,
   MEMSZ    nCol,
   double ** ppOldArray)
{
   MEMSZ    iRow;
 
   double *  pSpace;
   double ** ppArray;

   if (ppOldArray == NULL) {
      pSpace  = (double *)malloc(sizeof(double ) * nRow * nCol);
      ppArray = (double**)malloc(sizeof(double*) * nRow);
   } else {
      ppArray = (double**)realloc(ppOldArray, sizeof(double*) * nRow);
      pSpace  = (double *)realloc(ppOldArray[0],sizeof(double ) * nRow * nCol);
   }

   for (iRow = 0; iRow < nRow; iRow++) {
      /* Quantity (iRow*nCol) scales by sizeof(double) */
      ppArray[iRow] = pSpace + (iRow * nCol);
   }

   return ppArray;
}

/******************************************************************************/
/* Build an nRow x nCol matrix, in pointer-style.
 * Allocate one contiguous array of ints with  nRow*nCol elements.
 * Then create a set of nRow pointers, each of which points to the next nCol
 *  elements.
 */
int ** cciarray_rebuild_
  (MEMSZ    nRow,
   MEMSZ    nCol,
   int   ** ppOldArray)
{
   MEMSZ    iRow;
 
   int   *  pSpace;
   int   ** ppArray;

   if (ppOldArray == NULL) {
      pSpace  = (int *)malloc(sizeof(int ) * nRow * nCol);
      ppArray = (int**)malloc(sizeof(int*) * nRow);
   } else {
      ppArray = (int**)realloc(ppOldArray, sizeof(int*) * nRow);
      pSpace  = (int *)realloc(ppOldArray[0],sizeof(int ) * nRow * nCol);
   }

   for (iRow = 0; iRow < nRow; iRow++) {
      /* Quantity (iRow*nCol) scales by sizeof(int) */
      ppArray[iRow] = pSpace + (iRow * nCol);
   }

   return ppArray;
}

/******************************************************************************/
/* Free all memory allocated for an nRow x nCol matrix, as set up
 * with the routine ccarray_build_().
 */
void ccarray_free_
  (float ** ppArray,
   MEMSZ    nRow)
{
   /* Memory has been allocated in one contiguous chunk, so only free
    * ppArray[0] rather than every ppArray[i] for all i.
    */
   free((void *)ppArray[0]);
   free((void *)ppArray);
}

/******************************************************************************/
/* Free all memory allocated for an nRow x nCol matrix, as set up
 * with the routine ccdarray_build_().
 */
void ccdarray_free_
  (double ** ppArray,
   MEMSZ    nRow)
{
   /* Memory has been allocated in one contiguous chunk, so only free
    * ppArray[0] rather than every ppArray[i] for all i.
    */
   free((void *)ppArray[0]);
   free((void *)ppArray);
}

/******************************************************************************/
/* Free all memory allocated for an nRow x nCol matrix, as set up
 * with the routine cciarray_build_().
 */
void cciarray_free_
  (int   ** ppArray,
   MEMSZ    nRow)
{
   /* Memory has been allocated in one contiguous chunk, so only free
    * ppArray[0] rather than every ppArray[i] for all i.
    */
   free((void *)ppArray[0]);
   free((void *)ppArray);
}

/******************************************************************************/
/* Set all elements of an array equal to zero.
 * The array is assumed to be pointer-style, as allocated by ccarray_build_.
 */
void ccarray_zero_
  (float ** ppArray,
   MEMSZ    nRow,
   MEMSZ    nCol)
{
   MEMSZ    iRow;
   MEMSZ    iCol;

   for (iRow=0; iRow < nRow; iRow++) {
      for (iCol=0; iCol < nCol; iCol++) {
         ppArray[iRow][iCol] = 0.0;
      }
   }
}

/******************************************************************************/
/* Set all elements of a vector equal to zero.
 */
void ccvector_zero_
  (float *  pVector,
   MEMSZ    n)
{
   MEMSZ    i;
   for (i=0; i < n; i++) pVector[i] = 0.0;
}

/******************************************************************************/
/* Set all elements of a vector equal to zero.
 */
void ccdvector_zero_
  (double *  pVector,
   MEMSZ    n)
{
   MEMSZ    i;
   for (i=0; i < n; i++) pVector[i] = 0.0;
}

/******************************************************************************/
/* Set all elements of a vector equal to zero.
 */
void ccivector_zero_
  (int    *  pVector,
   MEMSZ    n)
{
   MEMSZ    i;
   for (i=0; i < n; i++) pVector[i] = 0;
}



 



/* Initialize file pointers to NULL -- this is done automatically since
 * this is an external variable.
 * Files are numbered 0 to IO_FOPEN_MAX-1.
 */

FILE  *  pFILEfits[IO_FOPEN_MAX];

/******************************************************************************/
/* Return IO_GOOD if a file exists, and IO_BAD otherwise.
 */
int inoutput_file_exist
  (char  *  pFileName)
{
   int      retval;

   if (access(pFileName, F_OK) == 0) {
      retval = IO_GOOD;
   } else {
      retval = IO_BAD;
   }
   return retval;
}

/******************************************************************************/
/* Return the index number of the first unused (NULL) file pointer.
 * If there are no free file pointers, then return IO_FOPEN_MAX.
 */
int inoutput_free_file_pointer_()
{
   int retval = 0;
 
   while(retval <= IO_FOPEN_MAX && pFILEfits[retval] != NULL) retval++;
   return retval;
}

/******************************************************************************/
/* 
 * Open a file for reading or writing.
 * A file number for the file pointer is returned in *pFilenum.
 *
 * Return IO_GOOD if the file was successfully opened, IO_BAD otherwise.
 * Failure can be due to either lack of free file pointers, or
 * the file not existing if pPriv=="r" for reading.
 */
int inoutput_open_file
  (int   *  pFilenum,
   char     pFileName[],
   char     pPriv[])
{
   int      iChar;
   int      retval;
   char     tempName[IO_FORTRAN_FL];

   if ((*pFilenum = inoutput_free_file_pointer_()) == IO_FOPEN_MAX) {
      printf("ERROR: Too many open files\n");
      retval = IO_BAD;
   } else {
      /* Truncate the Fortran-passed file name with a null,
       * in case it is padded with spaces */
      iChar = IO_FORTRAN_FL;
      strncpy(tempName, pFileName, iChar);
      for (iChar=0; iChar < IO_FORTRAN_FL; iChar++) {
         if (tempName[iChar] == ' ') tempName[iChar] = '\0';
      }

      retval = IO_GOOD;
      if (pPriv[0] == 'r') {
         if (inoutput_file_exist(tempName) == IO_GOOD) {
            if ((pFILEfits[*pFilenum] = fopen(tempName, pPriv)) == NULL) {
               printf("ERROR: Error opening file: %s\n", tempName);
               retval = IO_BAD;
            }
         } else {
            printf("ERROR: File does not exist: %s\n", tempName);
            retval = IO_BAD;
         }
      } else {
         if ((pFILEfits[*pFilenum] = fopen(tempName, pPriv)) == NULL) {
            printf("ERROR: Error opening file: %s\n", tempName);
            retval = IO_BAD;
         }
      }
   }

   return retval;
}

/******************************************************************************/
/*
 * Close a file.  Free the file pointer so it can be reused.
 * Return IO_BAD if an error is encountered; otherwise return IO_GOOD.
 */
int inoutput_close_file
  (int      filenum)
{
   int      retval;

   if (fclose(pFILEfits[filenum]) == EOF) {
      retval = IO_BAD;
   } else {
      retval = IO_GOOD;
   }
   pFILEfits[filenum] = NULL;

   return retval;
}

