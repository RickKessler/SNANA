// =======================================
//      SNANA control flags
//
// Sep 19 2024: S. Thorp
//    define OPT_MWCOLORLAW_FITZ99_EXACT = 9999 and leave original
//    OPT_MWCOLORLAW_FITZ99_APPROX = 99. After more testing, the plan is
//    make OPT_MWCOLORLAW_FITZ99_EXACT=99 the new default and allow
//    OPT_MWCOLORLAW_FITZ99_APPROX = -99 to revert back to old approximation.
// Sep 25 2024: S. Thorp, R. Kessler
//    define OPT_MWCOLORLAW_FITZ99_EXACT = 99
//    define OPT_MWCOLORLAW_FITZ99_APPROX = -99
//  Oct 19 2024: S. Thorp
//    define OPT_MWCOLORLAW_GORD03 = 203
//           OPT_MWCOLORLAW_FITZ04 = 204
//           OPT_MWCOLORLAW_GORD16 = 216
//           OPT_MWCOLORLAW_GORD23 = 223
// =======================================


#define OPT_MWCOLORLAW_OFF      0  // No Extinction applied.
#define OPT_MWCOLORLAW_CCM89   89  // Clayton,Cardelli,Matheson, 1989
#define OPT_MWCOLORLAW_FM90    90  // Fitzpatrick & Massa, 1990
#define OPT_MWCOLORLAW_ODON94  94  // O'Donnel 1994 update
#define OPT_MWCOLORLAW_FITZ99_APPROX  -99   // approx Fitzpatrick 1999 (D.Scolnic, 2013)
#define OPT_MWCOLORLAW_FITZ99_EXACT   99 // exact Fitzpatrick 1999 (S.Thorp, 2024)
#define OPT_MWCOLORLAW_GORD03  203 // Gordon et al. 2003 (S. Thorp, 2024)
#define OPT_MWCOLORLAW_FITZ04  204 // Fitzpatrick 2004 (S.Thorp, 2024)
#define OPT_MWCOLORLAW_GORD16  216 // Gordon et al. 2016 (S.Thorp, 2024)
#define OPT_MWCOLORLAW_GORD23  223 // Gordon et al. 2023 (S.Thorp, 2024)

#define OPT_MWEBV_OFF            0  // no extinction
#define OPT_MWEBV_FILE           1  // FILE value (simlib or data header)
#define OPT_MWEBV_SFD98          2  // use SFD98 value
#define OPT_MWEBV_Sch11_PS2013   3  // PS1-2013 implementation of Schlafly 2011

#define WAVEMAX_FITZ99 25000.0  // Oct 2021 Dillon and Dan switched from 12000


// =======================================
//      SNANA-interface functons
// =======================================

void MWgaldust(double RA,double DEC, double *avgal, double *EBV );

// functions moved from sntools.c (Sep 2013)
double GALextinct (double  RV, double  AV, double  WAVE, int  OPT);
double galextinct_(double *RV, double *AV, double *WAVE, int *OPT);
double GALextinct_Fitz99_exact(double RV, double AV, double WAVE, int OPT);
double GALextinct_FM90(double x, double c1, double c2, double c3, double c4,
                        double c5, double x02, double g2);
double GALextinct_Gord23(double RV, double AV, double WAVE);

// xxx mark double F99exact(double RV, double AV, double WAVE);

void   text_MWoption(  char *what, int  OPT, char *TEXT) ; // return TEXT
void   text_mwoption__(char *what, int *OPT, char *TEXT) ; 

void   modify_MWEBV_SFD  (int OPT, double RA, double DECL,    // (I)
			  double *MWEBV, double *MWEBV_ERR) ; // (I->O)
void   modify_mwebv_sfd__(int *OPT, double *RA, double *DECL,
			  double *MWEBV, double *MWEBV_ERR) ;

// =======================================
#ifndef __INCinterface_h
#define __INCinterface_h

 
/*
 * Deal with Fortran-C calling conventions.  Most unix Fortran
 * compilers add an extra trailing _ to fortran names, and so the
 * default values of FORTRAN_PREPEND and FORTRAN_APPEND are '' and '_'.
 */

/* Deal with Fortran--C name differences; the default is to add an _ */
#if defined(__SUN)
#  define FORTRAN_APPEND __
#endif
#if defined(__LINUX)
#  define FORTRAN_APPEND __
#endif
#if !defined(FORTRAN_APPEND)
#  define FORTRAN_APPEND _
#endif

#define _CONCATENATE(P,F) P ## F
#define CONCATENATE(P,F) _CONCATENATE(P,F)

#if defined(FORTRAN_PREPEND)
#   if defined(FORTRAN_APPEND)
#      define DECLARE(F) \
                CONCATENATE(FORTRAN_PREPEND,CONCATENATE(F,FORTRAN_APPEND))
#   else
#      define DECLARE(F) CONCATENATE(FORTRAN_PREPEND,F)
#   endif
#else
#   if defined(FORTRAN_APPEND)
#      define DECLARE(F) CONCATENATE(F,FORTRAN_APPEND)
#   else
#      define DECLARE(F) F
#   endif
#endif

#endif /* __INCinterface_h */


#ifndef __INCsubs_asciifile_h
#define __INCsubs_asciifile_h
 
int asciifile_read_colmajor
  (char     pFileName[],
   int      numColsMax,
   int   *  pNRows,
   int   *  pNCols,
   float ** ppData);
int asciifile_read_rowmajor
  (char     pFileName[],
   int      numColsMax,
   int   *  pNRows,
   int   *  pNCols,
   float ** ppData);
char * asciifile_read_line
  (int      filenum,
   int      numColsMax,
   int   *  pNValues,
   float *  values);

#endif /* __INCsubs_asciifile_h */

#ifndef __INCsubs_common_math_h
#define __INCsubs_common_math_h

#define min(a,b) ( ((a) < (b)) ? (a) : (b) )
#define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SIGN(a)  ( ((a) >= 0.0) ? (1.0) : (-1.0) )
#define TRUE_MWDUST  1
#define FALSE_MWDUST 0

#if 0
static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#endif

float FMIN
  (float    a,
   float    b);
float FMAX
  (float    a,
   float    b);
void vector_copy
  (int      nData,
   float *  pDataIn,
   float *  pDataOut);
float vector_minimum
  (int      nData,
   float *  pData);
float vector_maximum
  (int      nData,
   float *  pData);
float vector_mean
  (int      nData,
   float *  pData);
int ivector_sum
  (int      nData,
   int   *  pData);
float vector_sum
  (int      nData,
   float *  pData);
float vector_dispersion_about_zero
  (int      nData,
   float *  pData);
void vector_mean_and_dispersion
  (int      nData,
   float *  pData,
   float *  pMean,
   float *  pDispersion);
void vector_mean_and_dispersion_with_errs
  (int      nData,
   float *  pData,
   float *  pErr,
   float *  pMean,
   float *  pDispersion);
void vector_moments
  (int      nData,
   float *  pData,
   float *  pMean,
   float *  pVar,
   float *  pSkew,
   float *  pKurt);
int vector_closest_value
  (float    x,
   int      nData,
   float *  pData);
int vector_furthest_value
  (float    x,
   int      nData,
   float *  pData);
int vector_lower_value
  (float    x,
   int      nData,
   float *  pData);
int vector_higher_value
  (float    x,
   int      nData,
   float *  pData);
float vector_interpolated_index
  (float    x,
   int      nData,
   float *  pData);
float quadrature_sum
  (float    term1,
   float    term2);
void subtract_med_filter_vector
  (int      filtSize,
   int      nData,
   float *  pData);
void boxcar_filter_vector
  (int      filtSize,
   int      nData,
   float *  pData);
void median_filter_vector
  (int      filtSize,
   int      nData,
   float *  pData);
float vector_median
  (int      nData,
   float *  pData);
float vector_selip
  (int      kData,
   int      nData,
   float *  pData);
void ivector_sort_and_unique
  (int   *  pNData,
   int   ** ppData);
int vector_subtract_poly_fit_with_rej
  (int      nData,
   float *  pDataX,
   float *  pDataY,
   float    sigRejLo,
   float    sigRejHi,
   int      growRej,
   int      nRejIter,
   int      nCoeff);
int vector_poly_fit_with_rej
  (int      nData,
   float *  pDataX,
   float *  pDataY,
   float    sigRejLo,
   float    sigRejHi,
   int      growRej,
   int      nRejIter,
   int      nCoeff,
   float *  pCoeff,
   float *  pChi2);
float vector_percentile_value_with_rej
  (int      nData,
   float *  pDataVal,
   float *  pDataErr,
   float    fSelip,
   float    sigRejLo,
   float    sigRejHi,
   int      nRejIter);
void sort_structure_on_float
  (int      nData,
   int      structLen,
   float *  pFloatAddress,
   unsigned char *  pStructAddress);
int fit_exponential
  (int      nData,
   float *  pXaxis,
   float *  pYaxis,
   float *  pSigma,
   int      nCCoeff,
   float *  pCCoeff,
   int   *  pDOF,
   float *  pChi2);
void func_wderiv_exponential
  (float    x,
   float    a[],
   float *  y,
   float    dyda[],
   int      na);
 
#endif /* __INCsubs_common_math_h */

#ifndef __INCsubs_common_string_h
#define __INCsubs_common_string_h

#include <stdio.h> /* Include for definition of FILE */

/* removed Sep 9 2010
int getline
  (char  *  pLine,
   int      max);
*/

float sexagesimal_to_float
  (char  *  pTimeString);
void trim_string
  (char  *  pName,
   int      nChar);
char * get_next_line
  (char  *  pNextLine,
   int      maxLen,
   FILE  *  pFILE);
void string_to_integer_list
  (char  *  pString,
   int   *  pNEntry,
   int   ** ppEntry);

#endif /* __INCsubs_common_string_h */

#ifndef __INCsubs_fits_h
#define __INCsubs_fits_h

typedef unsigned char uchar;
typedef long int HSIZE;
typedef long int DSIZE;

void fits_read_file_fits_header_only_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead);
int fits_read_file_ascii_header_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead);
DSIZE fits_read_file_fits_r4_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   float ** ppData);
DSIZE fits_read_file_fits_i2_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   short int ** ppData);
DSIZE fits_read_subimg_
  (char     pFileName[],
   HSIZE    nHead,
   uchar *  pHead,
   DSIZE *  pStart,
   DSIZE *  pEnd,
   DSIZE *  pNVal,
   float ** ppVal);
void fits_read_subimg1
  (int      nel,
   DSIZE *  pNaxis,
   DSIZE *  pStart,
   DSIZE *  pEnd,
   int      fileNum,
   int      bitpix,
   DSIZE *  pNVal,
   uchar *  pData);
DSIZE fits_read_point_
  (char     pFileName[],
   HSIZE    nHead,
   uchar *  pHead,
   DSIZE *  pLoc,
   float *  pValue);
DSIZE fits_read_file_fits_noscale_
  (char     pFileName[],
   DSIZE *  pNHead,
   uchar ** ppHead,
   HSIZE *  pNData,
   int   *  pBitpix,
   uchar ** ppData);
DSIZE fits_read_file_xfits_noscale_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   HSIZE *  pNxhead,
   uchar ** ppXhead,
   DSIZE *  pNData,
   int   *  pBitpix,
   uchar ** data);
DSIZE fits_write_file_fits_r4_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   float ** ppData);
DSIZE fits_write_file_fits_i2_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   short int ** ppData);
DSIZE fits_write_file_fits_noscale_
  (char     pFileName[],
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   int   *  pBitpix,
   uchar ** ppData);
DSIZE fits_read_fits_data_
  (int   *  pFilenum,
   int   *  pBitpix,
   DSIZE *  pNData,
   uchar ** ppData);
DSIZE fits_write_fits_data_
  (int   *  pFilenum,
   int   *  pBitpix,
   DSIZE *  pNData,
   uchar ** ppData);
void fits_read_fits_header_
  (int   *  pFilenum,
   HSIZE *  pNHead,
   uchar ** ppHead);
void fits_skip_header_
  (int   *  pFilenum);
void fits_add_required_cards_
  (HSIZE *  pNHead,
   uchar ** ppHead);
void fits_write_fits_header_
  (int   *  pFilenum,
   HSIZE *  pNHead,
   uchar ** ppHead);
void fits_create_fits_header_
  (HSIZE *  pNHead,
   uchar ** ppHead);
void fits_duplicate_fits_header_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   uchar ** ppHeadCopy);
void fits_duplicate_fits_data_r4_
  (DSIZE *  pNData,
   float ** ppData,
   float ** ppDataCopy);
void fits_duplicate_fits_data_
  (int   *  pBitpix,
   DSIZE *  pNData,
   uchar ** ppData,
   uchar ** ppDataCopy);
void fits_create_fits_data_r4_
  (DSIZE *  pNData,
   float ** ppData);
void fits_create_fits_data_
  (int   *  pBitpix,
   DSIZE *  pNData,
   uchar ** ppData);
int fits_dispose_header_and_data_
  (uchar ** ppHead,
   uchar ** ppData);
int fits_dispose_array_
  (uchar ** ppHeadOrData);
DSIZE fits_compute_ndata_
  (HSIZE *  pNHead,
   uchar ** ppHead);
void fits_compute_axes_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   int   *  pNumAxes,
   DSIZE ** ppNaxis);
void fits_free_axes_
  (int   *  pNumAxes,
   DSIZE ** ppNaxis);
float compute_vista_wavelength_
  (DSIZE *  pPixelNumber,
   int   *  pNCoeff,
   float ** ppCoeff);
void fits_compute_vista_poly_coeffs_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   int   *  pNCoeff,
   float ** ppCoeff);
void fits_data_to_r4_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   uchar ** ppData);
void fits_data_to_i2_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   uchar ** ppData);
HSIZE fits_add_card_
  (uchar    pCard[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_add_cardblank_
  (uchar    pCard[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_add_card_ival_
  (int   *  pIval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_add_card_rval_
  (float *  pRval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_add_card_string_
  (char  *  pStringVal,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_add_card_comment_
  (char  *  pStringVal,
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_add_card_history_
  (char  *  pStringVal,
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_purge_blank_cards_
  (HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_delete_card_
  (uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_find_card_
  (uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
void fits_swap_cards_ival_
  (uchar *  pLabel1,
   uchar *  pLabel2,
   HSIZE *  pNHead,
   uchar ** ppHead);
void fits_swap_cards_rval_
  (uchar *  pLabel1,
   uchar *  pLabel2,
   HSIZE *  pNHead,
   uchar ** ppHead);
int fits_get_card_ival_
  (int   *  pIval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
int fits_get_card_rval_
  (float *  pRval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
int fits_get_card_date_
  (int   *  pMonth,
   int   *  pDate,
   int   *  pYear,
   uchar    pLabelDate[],
   HSIZE *  pNHead,
   uchar ** ppHead);
int fits_get_card_time_
  (float *  pTime,
   uchar    pLabelTime[],
   HSIZE *  pNHead,
   uchar ** ppHead);
int fits_get_card_string_
  (char  ** ppStringVal,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_change_card_
  (uchar    pCard[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_change_card_ival_
  (int   *  pIval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
HSIZE fits_change_card_rval_
  (float *  pRval,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHhead);
HSIZE fits_change_card_string_
  (char  *  pStringVal,
   uchar    pLabel[],
   HSIZE *  pNHead,
   uchar ** ppHead);
void fits_string_to_card_
  (uchar    pString[],
   uchar    pCard[]);
float fits_get_rval_
  (DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBscale,
   float *  pBzero,
   uchar ** ppDdata);
int fits_get_ival_
  (DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBscale,
   float *  pBzero,
   uchar ** ppDdata);
void fits_put_rval_
  (float *  pRval,
   DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBscale,
   float *  pBzero,
   uchar ** ppData);
int fits_qblankval_
  (DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBlankval,
   uchar ** ppDdata);
void fits_put_blankval_
  (DSIZE *  pIloc,
   int   *  pBitpix,
   float *  pBlankval,
   uchar ** ppData);
void fits_purge_nulls
  (uchar    pCard[]);
int fits_get_next_card_
  (int   *  pFilenum,
   uchar    pCard[]);
int fits_put_next_card_
  (int   *  pFilenum,
   uchar    pCard[]);
int fits_size_from_bitpix_
  (int   *  pBitpix);
void fits_pixshift_wrap_
  (int   *  pSAxis,
   DSIZE *  pShift,
   HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   uchar ** ppData);
void fits_transpose_data_
  (HSIZE *  pNHead,
   uchar ** ppHead,
   DSIZE *  pNData,
   uchar ** ppData);
void fits_ave_rows_r4_
  (int   *  iq,
   DSIZE *  pRowStart,
   DSIZE *  pNumRowAve,
   DSIZE *  pNaxis1,
   DSIZE *  pNaxis2,
   float ** ppData,
   float ** ppOut);
void fits_ave_obj_and_sigma_rows_r4_
  (int   *  iq,
   DSIZE *  pRowStart,
   DSIZE *  pNumRowAve,
   DSIZE *  pNaxis1,
   DSIZE *  pNaxis2,
   float ** ppObjData,
   float ** ppSigData,
   float ** ppObjOut,
   float ** ppSigOut);
void fits_byteswap
  (int      bitpix,
   DSIZE    nData,
   uchar *  pData);
void fits_bswap2
  (uchar *  pc1,
   uchar *  pc2);

extern uchar * datum_zero;
extern uchar * label_airmass;
extern uchar * label_bitpix;
extern uchar * label_blank;
extern uchar * label_bscale;
extern uchar * label_bzero;
extern uchar * label_ctype1;
extern uchar * label_ctype2;
extern uchar * label_cdelt1;
extern uchar * label_cdelt2;
extern uchar * label_cd1_1;
extern uchar * label_cd1_2;
extern uchar * label_cd2_1;
extern uchar * label_cd2_2;
extern uchar * label_latpole;
extern uchar * label_lonpole;
extern uchar * label_crpix1;
extern uchar * label_crpix2;
extern uchar * label_crval1;
extern uchar * label_crval2;
extern uchar * label_date_obs;
extern uchar * label_dec;
extern uchar * label_empty;
extern uchar * label_end;
extern uchar * label_exposure;
extern uchar * label_extend;
extern uchar * label_filtband;
extern uchar * label_filter;
extern uchar * label_ha;
extern uchar * label_instrume;
extern uchar * label_lamord;
extern uchar * label_loss;
extern uchar * label_naxis;
extern uchar * label_naxis1;
extern uchar * label_naxis2;
extern uchar * label_object;
extern uchar * label_observer;
extern uchar * label_pa;
extern uchar * label_platescl;
extern uchar * label_ra;
extern uchar * label_rnoise;
extern uchar * label_rota;
extern uchar * label_seeing;
extern uchar * label_skyrms;
extern uchar * label_skyval;
extern uchar * label_slitwidt;
extern uchar * label_st;
extern uchar * label_telescop;
extern uchar * label_time;
extern uchar * label_tub;
extern uchar * label_ut;
extern uchar * label_vhelio;
extern uchar * label_vminusi;
extern uchar * card_simple;
extern uchar * card_empty;
extern uchar * card_null;
extern uchar * card_end;
extern uchar * text_T;
extern uchar * text_F;

#endif /* __INCsubs_fits_h */

#ifndef __INCsubs_inoutput_h
#define __INCsubs_inoutput_h

#include <stdio.h> /* Include for definition of FILE */

#define MAX_FILE_LINE_LEN 500 /* Maximum line length for data files */
#define MAX_FILE_NAME_LEN  200
#define IO_FOPEN_MAX       20  // Files must be numbered 0 to IO_FOPEN_MAX-1
#define IO_FORTRAN_FL     200  // Max length of file name from a Fortran call
#define IO_GOOD             1
#define IO_BAD              0

extern FILE  *  pFILEfits[];

int inoutput_file_exist
  (char  *  pFileName);
int inoutput_free_file_pointer_();
int inoutput_open_file
  (int   *  pFilenum,
   char     pFileName[],
   char     pPriv[]);
int inoutput_close_file
  (int      filenum);

#endif /* __INCsubs_inoutput_h */

#ifndef __INCsubs_lambert_h
#define __INCsubs_lambert_h

void DECLARE(fort_lambert_getval)
  (char  *  pFileN,
   char  *  pFileS,
   void  *  pNGal,
   float *  pGall,
   float *  pGalb,
   void  *  pQInterp,
   void  *  pQNoloop,
   void  *  pQVerbose,
   float *  pOutput);
float * lambert_getval
  (char  *  pFileN,
   char  *  pFileS,
   long     nGal,
   float *  pGall,
   float *  pGalb,
   int      qInterp,
   int      qNoloop,
   int      qVerbose);
void lambert_lb2fpix
  (float    gall,   /* Galactic longitude */
   float    galb,   /* Galactic latitude */
   HSIZE    nHead,
   uchar *  pHead,
   float *  pX,     /* X position in pixels from the center */
   float *  pY);    /* Y position in pixels from the center */
void lambert_lb2pix
  (float    gall,   /* Galactic longitude */
   float    galb,   /* Galactic latitude */
   HSIZE    nHead,
   uchar *  pHead,
   int   *  pIX,    /* X position in pixels from the center */
   int   *  pIY);   /* Y position in pixels from the center */
void lambert_lb2xy
  (float    gall,   /* Galactic longitude */
   float    galb,   /* Galactic latitude */
   int      nsgp,   /* +1 for NGP projection, -1 for SGP */
   float    scale,  /* Radius of b=0 to b=90 degrees in pixels */
   float *  pX,     /* X position in pixels from the center */
   float *  pY);    /* Y position in pixels from the center */
int ivector_minimum
  (int      nData,
   int   *  pData);
int ivector_maximum
  (int      nData,
   int   *  pData);


extern uchar * label_lam_nsgp;
extern uchar * label_lam_scal;

#endif /* __INCsubs_lambert_h */

#ifndef __INCsubs_memory_h
#define __INCsubs_memory_h

typedef size_t MEMSZ;

#ifdef OLD_SUNCC
void memmove
  (void  *  s,
   const void  *  ct,
   MEMSZ    n);
#endif
void ccalloc_resize_
  (MEMSZ *  pOldMemSize,
   MEMSZ *  pNewMemSize,
   void  ** ppData);
void ccrealloc_
  (MEMSZ *  pMemSize,
   void  ** ppData);
void ccalloc_init
  (MEMSZ *  pMemSize,
   void  ** ppData);
void ccalloc_
  (MEMSZ *  pMemSize,
   void  ** ppData);
void ccfree_
  (void  ** ppData);
float * ccvector_build_
  (MEMSZ    n);
double * ccdvector_build_
  (MEMSZ    n);
int * ccivector_build_
  (MEMSZ    n);
float ** ccpvector_build_
  (MEMSZ    n);
float *** ccppvector_build_
  (MEMSZ    n);
float * ccvector_rebuild_
  (MEMSZ    n,
   float *  pOldVector);
double * ccdvector_rebuild_
  (MEMSZ    n,
   double *  pOldVector);
int * ccivector_rebuild_
  (MEMSZ    n,
   int   *  pOldVector);
float ** ccpvector_rebuild_
  (MEMSZ    n,
   float ** ppOldVector);
float *** ccppvector_rebuild_
  (MEMSZ    n,
   float *** pppOldVector);
void ccvector_free_
  (float *  pVector);
void ccdvector_free_
  (double *  pVector);
void ccivector_free_
  (int   *  pVector);
void ccpvector_free_
  (float ** ppVector);
void ccppvector_free_
  (float *** pppVector);
float ** ccarray_build_
  (MEMSZ    nRow,
   MEMSZ    nCol);
double ** ccdarray_build_
  (MEMSZ    nRow,
   MEMSZ    nCol);
int ** cciarray_build_
  (MEMSZ    nRow,
   MEMSZ    nCol);
float ** ccarray_rebuild_
  (MEMSZ    nRow,
   MEMSZ    nCol,
   float ** ppOldArray);
double ** ccdarray_rebuild_
  (MEMSZ    nRow,
   MEMSZ    nCol,
   double ** ppOldArray);
int ** cciarray_rebuild_
  (MEMSZ    nRow,
   MEMSZ    nCol,
   int   ** ppOldArray);
void ccarray_free_
  (float ** ppArray,
   MEMSZ    nRow);
void ccdarray_free_
  (double ** ppArray,
   MEMSZ    nRow);
void cciarray_free_
  (int   ** ppArray,
   MEMSZ    nRow);
void ccarray_zero_
  (float ** ppArray,
   MEMSZ    nRow,
   MEMSZ    nCol);
void ccvector_zero_
  (float *  pVector,
   MEMSZ    n);
void ccdvector_zero_
  (double *  pVector,
   MEMSZ    n);
void ccivector_zero_
  (int    *  pVector,
   MEMSZ    n);

#endif /* __INCsubs_memory_h */




// ======
