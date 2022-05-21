/* =================================================
  June, 2011  R.Kessler
  Separate snhost.h from the snhost.c file

  Feb 5 2015: int GALID -> long long GALID
       
  Aug 18 2015
    New struct HOSTLIB_ZPHOTEFF

  Sep 28 2015:  increase array bounds for WFIRST studies:
     NMAGPSF_HOSTLIB -> 8 (was 5)
     MXVAR_HOSTLIB   -> 200 (was 100)

 Feb 21 2016:  MINLOGZ_HOSTLIB -> -2.3 (was -2.0)

 May 22 2017: NMAGPSF_HOSTLIB->9 (was 8)

 Dec 30 2017: MALLOCSIZE_HOSTLIB -> 40000 (was 10000)

 Jan 15 2018:  MINLOGZ_HOSTLIB -> -2.523 (was -2.3)

 May 10 2018:  extend max z from 3.16 to 4.0
    NZPTR_HOSTLIB -> 320 (was 280)
    MAXLOGZ_HOSTLIB -> 0.6 (was 0.3)

 Feb 4 2019: add DLR and d_DLR

 Sep 19 2019:
   +  MINLOGZ_HOSTLIB -> -3.0 (was -2.523) for Dan/H0
   +  NZPTR_HOSTLIB -> HOSTLIB.NZPTR is a variable (not constant param) 
       computed from MAX/MIN LOGZ (no longer hard-wired)

 Apr 7 2020: MXWGT_HOSTLIB -> 50,000 (was 5,000)
 May 23 2020: add VPEC & VPEC_ERR (optional)
 Jun 02 2021: 
   + MXWGT_HOSTLIB -> 500,000 for GHOST (was 50k)
   + MXCHAR_LINE_HOSTLIB -> 800 (was 600)
   + define MXCHAR_LINE_APPEND

 Jul 1 2021
   + MXWGT_HOSTLIB -> 20 million (was 500,000) for GHOST [tested by AG]

 Feb 9 2022: MXVAR_HOSTLIB -> 300 (was 200) to allow for up to
             100 zPHOT quantiles

 May 5 2022: MXCHAR_LINE_HOSTLIB->900

==================================================== */

#define HOSTLIB_MSKOPT_USE           1 // internally set if HOSTLIB_FILE
#define HOSTLIB_MSKOPT_GALMAG        2 // get host mag & noise in SN aperture
#define HOSTLIB_MSKOPT_SNMAGSHIFT    4 // adust SN mag from wgtmap
#define HOSTLIB_MSKOPT_SN2GAL_RADEC  8 // transfer SN coords to galaxy coords
#define HOSTLIB_MSKOPT_SN2GAL_Z     16 // transfer SN Z to ZTRUE
#define HOSTLIB_MSKOPT_USEONCE      32 // use host galaxy only once
#define HOSTLIB_MSKOPT_USESNPAR     64 // use SN color & shape from hostlib
#define HOSTLIB_MSKOPT_SWAPZPHOT   128 // swap ZTRUE with ZPHOT
#define HOSTLIB_MSKOPT_USEVPEC     512 // use VPEC from hostlib 

#define HOSTLIB_MSKOPT_VERBOSE     256 // print extra info during init
#define HOSTLIB_MSKOPT_DUMP       1024 // screen-dump for each host 
#define HOSTLIB_MSKOPT_DUMPROW    2048 // DUMP 1 row per host for parsing

#define HOSTLIB_MSKOPT_APPEND     4096  // append columns from file
#define HOSTLIB_MSKOPT_PLUSMAGS   8192  // compute & write host mags from host spectra
#define HOSTLIB_MSKOPT_PLUSNBR   16384  // append list of nbr to HOSTLIB
#define HOSTLIB_MSKOPT_ZPHOT_QGAUSS 32768  // write Gauss quantiles for zPHOT

#define HOSTLIB_1DINDEX_ID 10    // ID for 1DINDEX transformations

#define MXCHAR_LINE_HOSTLIB 900  // max number of chars per HOSTLIB line
#define MXCHAR_LINE_APPEND  500  // max number of appended chars per line
#define MXVAR_HOSTLIB       400  // max number of variables (NVAR:) in HOSTLIB
#define MXVAR_WGTMAP_HOSTLIB 10  // max no. weight-map variables
#define MXROW_WGTMAP      25000000  // 20 million, Alex Gagliano 09/2021
#define MXROW_HOSTLIB     10000000  // 10 million, Alex Gagliano 09/2021
#define MXCHECK_WGTMAP     1000  // max no. galaxies to check wgt map
#define MALLOCSIZE_HOSTLIB 40000 // incremental size of internal HOSTLIB array
#define MXCOMMENT_HOSTLIB  40    // max number of lines for README file
#define MXGauss2dTable     200   // max length of Gauss2d table
#define NVAR_Gauss2d       3     // Number of variables in Gauss2d table
#define MXBIN_ZPHOTEFF     100   // 

#define NSERSIC_TABLE        50    // number of integral tables
#define SERSIC_INDEX_MIN   0.15
//#define SERSIC_INDEX_MAX   8.00    // increase from 5 (6/24/2015)
#define SERSIC_INDEX_MAX  11.00    // increase from 8 (Mar 2 2022)
#define MXSERSIC_HOSTLIB      9    // max number of summed profiles per host
#define NBIN_RADIUS_SERSIC  200    // Number of R/Re bins to store integrals
#define MAXRADIUS_SERSIC   100.0   // max R/Re value for integ table
#define MINRADIUS_SERSIC  1.0E-4   // min R/Re value for integ table

#define NRBIN_GALMAG        100    // No. of radius bins for Galmag 
#define NTHBIN_GALMAG        36    // No. theta bins for galmag
#define MXBIN_SERSIC_bn     2000   // max bins in Sersic_bn file
 
// hard wire logarithmic z-bins
#define DZPTR_HOSTLIB      0.01   // logz-binning for Z-pointers
#define MINLOGZ_HOSTLIB   -3.00    // zmin = 0.001
#define MAXLOGZ_HOSTLIB    0.61    // zmax = 4.07
#define LOGZRANGE_HOSTLIB  MAXLOGZ_HOSTLIB-MINLOGZ_HOSTLIB
#define ZMIN_HOSTLIB       pow(10.0,MINLOGZ_HOSTLIB)
#define ZMAX_STAR          0.001   // give warnings for ZTRUE < ZMAX_STAR

#define NMAGPSF_HOSTLIB    9    // number of aperture mags vs. PSF to compute
#define DEG_ARCSEC    1./3600.  // 1 arcsec in deg.
#define DEBUG_WGTFLUX2    0.0    // fix 2nd WGTFLUX if non-zero

#define MXUSE_SAMEGAL 50     // max number of times to re-use hostGal
                             // with MINDAYSEP_SAMEGAL option

// define required keys in the HOSTLIB
#define HOSTLIB_VARNAME_GALID     "GALID"  // required 
#define HOSTLIB_VARNAME_ZTRUE     "ZTRUE"  // required

// define optional keys
#define HOSTLIB_VARNAME_TRUE_MATCH   "TRUE"
#define HOSTLIB_VARNAME_ZPHOT        "ZPHOT"
#define HOSTLIB_VARNAME_ZPHOT_ERR    "ZPHOT_ERR" 
#define HOSTLIB_VARNAME_VPEC         "VPEC"         
#define HOSTLIB_VARNAME_VPEC_ERR     "VPEC_ERR"     
#define HOSTLIB_VARNAME_LOGMASS_TRUE "LOGMASS_TRUE"  // log10(Mgal/Msolar)
#define HOSTLIB_VARNAME_LOGMASS_ERR  "LOGMASS_ERR" 
#define HOSTLIB_VARNAME_LOGMASS_OBS  "LOGMASS_OBS"
#define HOSTLIB_VARNAME_RA           "RA"  
#define HOSTLIB_VARNAME_DEC          "DEC" 
#define HOSTLIB_VARNAME_RA_HOST      "RA_HOST"  
#define HOSTLIB_VARNAME_DEC_HOST     "DEC_HOST" 
#define HOSTLIB_VARNAME_RA_GAL       "RA_GAL"  
#define HOSTLIB_VARNAME_DEC_GAL      "DEC_GAL" 
#define HOSTLIB_VARNAME_ANGLE        "a_rot"    // rotation angle
#define HOSTLIB_VARNAME_FIELD        "FIELD" 
#define HOSTLIB_VARNAME_NBR_LIST     "NBR_LIST" // Nov 2019
#define HOSTLIB_VARNAME_ELLIPTICITY  "ellipticity" // Sept 2021 Alex Gagliano
#define HOSTLIB_VARNAME_GALID2       "GALID2"
#define HOSTLIB_VARNAME_SQRADIUS     "sqradius"
#define HOSTLIB_SUFFIX_MAGOBS        "_obs"     // key = [filt]$SUFFIX
#define HOSTLIB_SUFFIX_MAGOBS_ERR    "_obs_err"     // key = [filt]$SUFFIX
#define HOSTLIB_PREFIX_ZPHOT_Q       PREFIX_ZPHOT_Q // see sndata.h
#define HOSTLIB_VARNAME_A_DLR        "a_DLR" // use this to measure DLR
#define HOSTLIB_VARNAME_B_DLR        "b_DLR"


// for SNMAGSHIFT, allow hostlib param instead of wgtmap.
// To save storage memory, SNMAGSHIFT is stored as 2 byte short int 
// with I2 = SNMAGSHIFT * I2SNMAGSHIFT_HOSTLIB 
#define HOSTLIB_VARNAME_SNMAGSHIFT  "SNMAGSHIFT" 
#define I2MAGSCALE_HOSTLIB 20000.0

// define ascii files with tables needed for numerical computation
#define FILENAME_Gauss2d    "$SNDATA_ROOT/simlib/Gauss2dIntegrals.dat" 
#define FILENAME_Sersic_bn  "$SNDATA_ROOT/simlib/Sersic_bn.dat" 

int NCALL_GEN_SNHOST_DRIVER ;
char PATH_DEFAULT_HOSTLIB[2*MXPATHLEN]; // e.g., $SNDATA_ROOT/simlib

#define MXTMPWORD_HOSTLIB MXVAR_HOSTLIB  // xxx 100
char *TMPWORD_HOSTLIB[MXTMPWORD_HOSTLIB]; // used for splitString

int OPTMASK_OPENFILE_HOSTLIB ;

//for developers only
bool REFAC_HOSTLIB;

/* xxx mark delete Apr 24 2022 
// define generic host properties e.g. LOGMASS, LOGsSFR, COLOR...
#define HOSTGAL_PROPERTY_BASENAME_LOGMASS  "LOGMASS"
#define HOSTGAL_PROPERTY_BASENAME_LOGSFR   "LOGSFR"
#define HOSTGAL_PROPERTY_BASENAME_LOGsSFR  "LOGsSFR"
#define HOSTGAL_PROPERTY_BASENAME_COLOR    "COLOR"

#define HOSTGAL_PROPERTY_NAME_LIST HOSTGAL_PROPERTY_BASENAME_LOGMASS " " HOSTGAL_PROPERTY_BASENAME_LOGSFR " " HOSTGAL_PROPERTY_BASENAME_LOGsSFR " " HOSTGAL_PROPERTY_BASENAME_COLOR
xxxxxx end mark xxxxx */

// Mar 16 2022: beware that -9 for sSFR is valid, so hostless sSFR is -99;
// the other hostless values are -9 as before.
#define HOSTLESS_PROPERTY_VALUE_LIST (float[]){ -9.0, -9.0, -99.0, -9.0 }

int N_HOSTGAL_PROPERTY;


typedef struct {
  double VAL_TRUE, VAL_OBS, VAL_ERR;
} HOSTGAL_PROPERTY_VALUE_DEF;

typedef struct {
  int  IVAR_TRUE, IVAR_OBS, IVAR_ERR;
  char BASENAME[100];
  double SCALE_ERR;
} HOSTGAL_PROPERTY_IVAR_DEF;


struct HOSTLIB_DEF {
  char FILENAME[MXPATHLEN] ; // full file name of HOSTLIB
  int  GZIPFLAG;

  int  NGAL_READ  ; // total number of galaxies read
  int  NGAL_STORE ; // number of GAL entries read and stored
  
  int  NVAR_REQUIRED ;
  int  NVAR_OPTIONAL ;
  int  NVAR_ALL    ; // total no, variables in HOSTLIB
  int  NVAR_STORE  ; // NVAR wgtmap + required + optional
  int  NERR_NAN;     // number of NaN read (Apr 2021)

  char VARNAME_REQUIRED[MXVAR_HOSTLIB][40];
  char VARNAME_OPTIONAL[MXVAR_HOSTLIB][40];
  char VARNAME_ALL[MXVAR_HOSTLIB][40];  // all names
  char VARNAME_ORIG[MXVAR_HOSTLIB][40]; // all original names
  char VARNAME_STORE[MXVAR_HOSTLIB][40];  
  
  int  IVAR_ALL[MXVAR_HOSTLIB];  // [sparse store index] = ALL-ivar

  double VALMIN[MXVAR_HOSTLIB];  // vs. IVAR_STORE
  double VALMAX[MXVAR_HOSTLIB]; 

  int  NVAR_SNPAR  ; // subset of SN parameters (e.g., shape, color ...)
  char VARSTRING_SNPAR[100] ;     // string of optional SN par names 
  int  IS_SNPAR_OPTIONAL[MXVAR_HOSTLIB];   // flag subset that are SN par
  int  FOUND_SNPAR_OPTIONAL[MXVAR_HOSTLIB] ;
  int  IS_SNPAR_STORE[MXVAR_HOSTLIB];   // flag subset that are SN par

  // define pointers used to malloc memory with MALLOCSIZE_HOSTLIB
  double *VALUE_ZSORTED[MXVAR_HOSTLIB];  // sorted by redshift
  double *VALUE_UNSORTED[MXVAR_HOSTLIB]; // same order as in HOSTLIB
  int    *LIBINDEX_UNSORT;    // map between z-sorted and unsorted (w/cuts)
  int    *LIBINDEX_ZSORT;     // inverse map 
  int     SORTFLAG ; // 1=> sorted

  int   *INDEX_FIELD; // corresponds to user-inputs HOSTLIB_FIELDMATCH
  char **FIELD_UNSORTED ;
  char **FIELD_ZSORTED ;

  char **NBR_UNSORTED ; // read from NBR_LIST column, Nov 11 2019
  char **NBR_ZSORTED ;

  int *LIBINDEX_READ; // map between read index (no cuts) and unsorted

  int MALLOCSIZE_D, MALLOCSIZE_I, MALLOCSIZE_Cp ;
  int NGAL_STORE_MALLOC ;

  // pointers to stored variables
  int IVAR_GALID ;
  int IVAR_TRUE_MATCH ;  // optional column: 1->use for true match
  int IVAR_ZTRUE  ;
  int IVAR_ZPHOT ;
  int IVAR_ZPHOT_ERR  ;
  int IVAR_ZPHOT_Q0; // index of first ZPHOT_Q (not necessarily 0th quantile)
  int NZPHOT_Q;
  int IVAR_VPEC ;
  int IVAR_VPEC_ERR  ;
  int IVAR_LOGMASS_TRUE ; // legacy 
  int IVAR_LOGMASS_ERR ; // legacy 
  int IVAR_LOGMASS_OBS ;  // legacy
  HOSTGAL_PROPERTY_IVAR_DEF *HOSTGAL_PROPERTY_IVAR ;
  int IVAR_RA ;
  int IVAR_DEC ; 
  int IVAR_ANGLE ;  // rot angle of a-axis w.r.t. RA
  int IVAR_FIELD ;                // optional FIELD key
  int IVAR_NBR_LIST;              // NBR_LIST column added by +HOSTNBR arg
  int IGAL_NBR_LIST;              // AG 08/2021
  int IVAR_GALID2;                // AG 09/2021
  int IVAR_ELLIPTICITY;
  int IVAR_SQRADIUS;
  int IVAR_a[MXSERSIC_HOSTLIB];   // semi-major  half-light
  int IVAR_b[MXSERSIC_HOSTLIB];   // semi-minor 
  int IVAR_w[MXSERSIC_HOSTLIB];   // weight
  int IVAR_n[MXSERSIC_HOSTLIB];   // Sersic index
  int IVAR_a_DLR;   // to measure DLR: e.g. a_IMAGE from sextractor
  int IVAR_b_DLR;   // to measure DLR
  int IVAR_MAGOBS[MXFILTINDX] ;     // pointer to oberver-mags
  int IVAR_MAGOBS_ERR[MXFILTINDX] ; // pointer to obs-mag errs (Aug 6 2021)
  int IVAR_WGTMAP[MXVAR_HOSTLIB] ;  // wgtmap-ivar vs [ivar_STORE]
  int IVAR_STORE[MXVAR_HOSTLIB]  ;  // store-ivar vs [ivarmap]
  int NFILT_MAGOBS;  // NFILT with host mag info read

  char filterList[MXFILTINDX]; // filter list for gal-mag
  char VARNAME_ZPHOT_Q[MXBIN_ZPHOTEFF][12];
  int  PERCENTILE_ZPHOT_Q[MXBIN_ZPHOTEFF]; // list of percentiles
  double SIGMA_QGAUSS[MXBIN_ZPHOTEFF];   // for forced Gauss quantiles

  // redshift information
  double ZMIN,ZMAX ;         // helio
  double ZGAPMAX ;           // max z-gap in library
  double ZGAPAVG ;           // avg z-gap in library
  double Z_ATGAPMAX[2];  // redshift at max ZGAP (to find big holes)
  int    NSTAR;          // number of entries with ZTRUE < ZMAX_STAR

  // vpec info (to print)
  double FIX_VPEC_ERR; // optional read from header
  double VPEC_RMS, VPEC_AVG, VPEC_MIN, VPEC_MAX; // info to print

  int   NZPTR;
  int  *IZPTR;            // pointers to nearest z-bin with .01 bin-size
  int   MINiz, MAXiz ;    // min,max valid iz arg for IZPTR

  int NLINE_COMMENT ;
  char COMMENT[MXCOMMENT_HOSTLIB][120] ; // comment lines for README file.

  // PSF-aperture info
  double Aperture_Radius[NMAGPSF_HOSTLIB+1]; // integ. radius for each PSF
  double Aperture_PSFSIG[NMAGPSF_HOSTLIB+1]; // fixed list of sigma-PSFs (asec)
  double Aperture_Rmax ;   // max radius (arsec) to integrage gal flux
  double Aperture_Rbin ;   // radial integration bins size (arcsec)
  double Aperture_THbin ;  // azim. integration binsize (radians)

  // pre-computed 2d-Gaussian integrals
  int    NGauss2d ;        // Number of Gauss2d entries
  int    NBIN_Gauss2dRadius ;
  int    NBIN_Gauss2dSigma ;
  double Gauss2dTable[NVAR_Gauss2d+1][MXGauss2dTable] ;
  double Gauss2dRadius[3]; // binsize, min, max
  double Gauss2dSigma[3];  // binsize, min, max

  // pre-computed cos and sin to speed gal-flux integration
  double Aperture_cosTH[NTHBIN_GALMAG+1] ;
  double Aperture_sinTH[NTHBIN_GALMAG+1] ;

  long long IGAL_FORCE; // set if HOSTLIB_GALID_FORCE is set

} HOSTLIB ;


#define MXCHAR_NBR_LIST 200 // Apr 25 2022 -> 200 (was 100)
#define MXNBR_LIST       50

struct {
  // optional command-line inputs
  double SEPNBR_MAX;     // optional command-line input (default=10 arcsec)
  int    NNBR_WRITE_MAX;  // idem for how many NBRs to write (default=10)

  // internal arrays for +HOSTNBR command-line option
  int    NNBR_MAX; // actual max of NNBR
  double *SKY_SORTED_DEC, *SKY_SORTED_RA ; 
  int    *SKY_SORTED_IGAL_zsort;
  int    *SKY_SORTED_IGAL_DECsort;
  long long GALID_atNNBR_MAX;  

  // internal diagnostics for stdout dump
  int NGAL_PER_NNBR[100] ; // store histogram of NNBR distribution
  int NGAL_TRUNCATE;       // NGAL cliiped by NNBR_WRITE_MAX or MXCHAR

} HOSTLIB_NBR_WRITE ;


struct {
  double ZWIN[2], RAWIN[2], DECWIN[2];
} HOSTLIB_CUTS;


struct SAMEHOST_DEF {
  int REUSE_FLAG ;          // 1-> re-use host
  unsigned short  *NUSE ;     // number of times each host is used.

  // define array to store all PEAKMJDs for each host; allows re-using
  // host after NDAYDIF_SAMEGAL. Note 2-byte integers to save memory
  unsigned short **PEAKDAY_STORE ; // add PEAKMJD_STORE_OFFSET to get PEAKMJD
  int PEAKMJD_STORE_OFFSET ;       // min generated PEAKMJD
} SAMEHOST ;

// Sersic quantities to define galaxy profile
// these are all defined during init
struct SERSIC_PROFILE_DEF {
  int  NPROF ;    // number of defined Sersic/profile components  
  char VARNAME_a[MXSERSIC_HOSTLIB][12];     // name of major axis; i.e, a1
  char VARNAME_b[MXSERSIC_HOSTLIB][12];     // name of minor axis; i.e, b1
  char VARNAME_w[MXSERSIC_HOSTLIB][12];     // name of weight
  char VARNAME_n[MXSERSIC_HOSTLIB][12];     // name of index

  int  IVAR_a[MXSERSIC_HOSTLIB];        // pointer to 'a' values
  int  IVAR_b[MXSERSIC_HOSTLIB];        // 
  int  IVAR_w[MXSERSIC_HOSTLIB];        // 
  int  IVAR_n[MXSERSIC_HOSTLIB];        // 

  double  FIXn[MXSERSIC_HOSTLIB];

  int    NFIX; // number of Sersic profiles with fixed index
  double FIX_VALUE[MXSERSIC_HOSTLIB];     // list of fixed Sersic indices
  char   FIX_NAME[MXSERSIC_HOSTLIB][12];  // list of fixed Sersic indices

} SERSIC_PROFILE ;


// Sersic integral tables to speed calculations
struct SERSIC_TABLE_DEF {
  int     TABLEMEMORY ;  // total memory of integral tables (bytes)
  double  INVINDEX_MIN ; // min 1/n 
  double  INVINDEX_MAX ; // max 1/n
  double  INVINDEX_BIN ; // binsize of 1/n

  double  inv_n[NSERSIC_TABLE+1] ; // list of 1/n
  double  n[NSERSIC_TABLE+1] ;     // i.e, n=4 for deVauc, n=1 for exponent ...
  double  bn[NSERSIC_TABLE+1] ;   // calculated b_n values

  double *INTEG_CUM[NSERSIC_TABLE+1];   // cumulative integral per Sersic index
  double  INTEG_SUM[NSERSIC_TABLE+1];   // store total integral for each index

  int  BIN_HALFINTEGRAL[NSERSIC_TABLE+1]; // bin where integral is total/2

  int  NBIN_reduced ;  // number of reduced R/Re bins
  double *reduced_logR ;     // list of R/Re upper-interal limit
  double  reduced_logRmax ;  // max R/Re in table
  double  reduced_logRmin ;  // max R/Re in table
  double  reduced_logRbin ;  // bin size  of R/Re

  // define table-grid of b_n vs. n read from ascii file FILENAME_SERSIC_BN
  int    Ngrid_bn ;
  double grid_n[MXBIN_SERSIC_bn] ;
  double grid_bn[MXBIN_SERSIC_bn] ;

} SERSIC_TABLE ;


struct HOSTLIB_ZPHOTEFF_DEF {
  int    NZBIN ;
  double ZTRUE[MXBIN_ZPHOTEFF] ;
  double EFF[MXBIN_ZPHOTEFF] ;
} HOSTLIB_ZPHOTEFF ; 


// define the weight  map .. 
struct HOSTLIB_WGTMAP_DEF {

  // GRID value vs. [ivar][igrid]

  // parameters describing each WGT
  char  VARNAME[MXVAR_HOSTLIB][40]; 
  bool  READSTAT ;    // T => wgtmap has been read
  bool  USE_SALT2GAMMA_GRID;
  bool  FOUNDVAR_SNMAGSHIFT;  // T -> SNMAGSHIFT column was found

  double MEMTOT_MB;  // memory (MB) allocated for wgtmap

  // params for SN variables that are NOT in the hostlib (Mar 2020)
  int   IDMAP_INDEX_SNVAR;       // for get_1DINDEX
  int   N_SNVAR ;                // number of SN properties
  bool  IS_SNVAR[MXVAR_WGTMAP_HOSTLIB]; // true -> SN property NOT in HOSTLIB
  int   ISPARSE_SNVAR[MXVAR_WGTMAP_HOSTLIB]; //sparse list point to ivar_WGTMAP
  int   INVSPARSE_SNVAR[MXVAR_HOSTLIB]; // HOSTLIB ivar -> sparse WGTMAP ivar

  char  VARNAME_SNVAR[MXVAR_WGTMAP_HOSTLIB][40];
  int   NB1D_SNVAR[MXVAR_WGTMAP_HOSTLIB];      // numer of bins per SNvar
  int   NBTOT_SNVAR;                         // product of NB1D
  int  *IBIN1D_SNVAR[MXVAR_WGTMAP_HOSTLIB];  // global 1D bin -> bin per var
  double *VALGRID_SNVAR[MXVAR_WGTMAP_HOSTLIB];   // value vs. global 1D bin

  // SNVAR variables update each event
  double *ptrVal_SNVAR[MXVAR_WGTMAP_HOSTLIB];  // value(s) for each event
  int     ibin_SNVAR ;                         // SNVAR bin per event

  // weigt storage for each galaxy
  double  WGTMAX ; // max weight for entire  wgtmap
  double *WGTSUM ;      // cumulative sum of weights over entire HOSTLIB
  //double *WGT ;         // wgt for each hostlib entry
  short int *I2SNMAGSHIFT ;  // SN mag shift for each hostlib entry
  double **WGTSUM_SNVAR;              // idem
  short int **I2SNMAGSHIFT_SNVAR ;    // idem

  // define  arrays to store list of GALIDs to check wgtmap interpolation
  int      NCHECKLIST ;
  int       CHECKLIST_IGAL[MXCHECK_WGTMAP] ;  // sparse pointer
  long long CHECKLIST_GALID[MXCHECK_WGTMAP] ; // absolute GALID
  double    CHECKLIST_ZTRUE[MXCHECK_WGTMAP] ;
  double    CHECKLIST_WGT[MXCHECK_WGTMAP] ;
  double    CHECKLIST_SNMAG[MXCHECK_WGTMAP] ;

  struct  GRIDMAP  GRIDMAP ;       // all WGTMAP vars

  int OPT_EXTRAP; // 1 ==> pull out-of-range values to edge of grid

} HOSTLIB_WGTMAP ;



typedef struct { 
  int     NPROF; // number of Sersic profiles 
  // Sersic profiles for this host
  double  INDEX ; // Sersic 'n[JPROF]'  for selected JPROF term
  double  a[MXSERSIC_HOSTLIB]  ;
  double  b[MXSERSIC_HOSTLIB]  ;
  double  n[MXSERSIC_HOSTLIB]  ;
  double  w[MXSERSIC_HOSTLIB]  ;
  double  wsum[MXSERSIC_HOSTLIB] ;
  double  bn[MXSERSIC_HOSTLIB] ;

  double a_rot; // rot angle (deg) w.r.t. RA
} SERSIC_DEF ; // created Nov 2019


// SNHOSTGAL below contains information about TRUE host;
// here we define information for each nerby galaxy that
// is sorted by DDRL
typedef struct {
  long long GALID ;
  double ZPHOT, ZPHOT_ERR ;     // photoZ of host
  double ZSPEC, ZSPEC_ERR ;     // ZTRUE
  double RA, DEC, SNSEP, DLR, DDLR ;
  double LOGMASS_TRUE, LOGMASS_ERR, LOGMASS_OBS ; // legacy
  HOSTGAL_PROPERTY_VALUE_DEF *HOSTGAL_PROPERTY_VALUE ;
  double MAG[MXFILTINDX]; 
  double MAG_ERR[MXFILTINDX];
  double ZPHOT_Q[MXBIN_ZPHOT_Q];
  bool   TRUE_MATCH ;
  // Added for LSST but maybe of more general utility
  // Alex Gagliano 09/2021
  long long GALID2 ; // Second ID e.g., from external catalog
  double SQRADIUS; // Ixx + Iyy
  double ELLIPTICITY;

  int GALID_UNIQUE; // see input HOSTLIB_GALID_UNIQUE_OFFSET (Oct 2021)

} SNHOSTGAL_DDLR_SORT_DEF ;

SNHOSTGAL_DDLR_SORT_DEF SNHOSTGAL_DDLR_SORT[MXNBR_LIST] ;


// define structure to hold information for one event ...
// gets over-written for each generated SN
struct SNHOSTGAL {

  int   IGAL  ;     // sequential sparse galaxy index 
  int   IGAL_SELECT_RANGE[2] ; // range to select random IGAL

  long long GALID ;   // Galaxy ID from library
  int  IMATCH_TRUE ; // true index for SNHOSTGAL_DDLR_SORT

  // redshift info
  double ZGEN  ;     // saved ZSN passed to driver
  double ZTRUE ;     // host galaxy redshift 
  double ZDIF ;      // zSN(orig) - zGAL, Nov 2015
  double ZPHOT, ZPHOT_ERR ;     // photoZ of host
  double ZPHOT_Q[MXBIN_ZPHOT_Q];
  double ZSPEC, ZSPEC_ERR ;     // = zSN or z of wrong host
  double VPEC,  VPEC_ERR  ;     // peculiar velocity
  double PEAKMJD ;

  int    NNBR;    // number of nearby galaxies
  int    IGAL_NBR_LIST[MXNBR_LIST];   // IGAL list of neighbors
  double DDLR_NBR_LIST[MXNBR_LIST];   // DDLR per NBR
  double SNSEP_NBR_LIST[MXNBR_LIST];

  SERSIC_DEF SERSIC ; // Nov 2019

  // coordinate info

  double reduced_R ;     // reduced R/Re (randomly chosen)
  double phi ;           // randomly chosen azim. angle, radians  

  double a_SNGALSEP_ASEC ;  // angle-coord along major axis
  double b_SNGALSEP_ASEC ;  // idem for minor axis

  double RA_GAL_DEG ;       // Galaxy sky coord from library (DEG)
  double DEC_GAL_DEG ;  
  double RA_SN_DEG ;            // SN sky coord (DEG)
  double DEC_SN_DEG ;   
  double RA_SNGALSEP_ASEC ;     // SN-galaxy sep in RA, arcsec
  double DEC_SNGALSEP_ASEC ;    // idem in DEC
  double SNSEP ;        // SN-gal sep, arcsec
  double DLR ;          // directional light radius
  double DDLR;          // SNSEP/DLR (following Gupta 2016)

  // aperture-mag info
  double SB_MAG[MXFILTINDX] ;  // surface brightness mag in 1 sq-arcsec
  double SB_FLUXCAL[MXFILTINDX] ;

  // xxx delete Jan 31 2020  double GALMAG_TOT[MXFILTINDX];  
  double GALMAG[MXFILTINDX][NMAGPSF_HOSTLIB+1] ; // mag per PSF bin
  double GALFRAC[NMAGPSF_HOSTLIB+1]; // true gal light-frac in each aperture
  double GALFRAC_SBRADIUS[NMAGPSF_HOSTLIB+1]; // gal light-frac in SB radius
  double WGTMAP_SNMAGSHIFT ;        // SN mag shift from wgtmap
  double WGTMAP_WGT ;               // selection weight

  // parameters used to interpolate selection WGT
  double WGTMAP_VALUES[MXVAR_HOSTLIB]; 

  // random numbers (0-1)
  double FlatRan1_GALID ;
  double FlatRan1_radius[2] ;
  double FlatRan1_phi ;  // random 0-1, not 0 to 2*PI

} SNHOSTGAL ;


// store arbitrary extra variables to write out to data file
struct SNTABLEVAR_DEF {
  int    NOUT ;  // number of variables to write out
  int    IVAR_STORE[MXVAR_HOSTLIB] ;
  char   NAME[MXVAR_HOSTLIB][40];
  int    USED_IN_WGTMAP[MXVAR_HOSTLIB]; // to avoid printing same var twice

  // values updated for each event
  // add second dimension -AG 08/2021
  double VALUE[MXVAR_HOSTLIB][MXNBR_LIST] ;   
} HOSTLIB_OUTVAR_EXTRA ;


// Jun 27 2019; define structed used to determine host spectrum
#define MXSPECBASIS_HOSTLIB 20  // max number of spectral templates
#define MXBIN_SPECBASIS  20000  // max number of wave bins
#define PREFIX_SPECBASIS "SPECBASIS"    // e.g., specbasis00, specbasis01, ...
#define PREFIX_SPECBASIS_HOSTLIB "COEFF_"      // extra prefix for HOSTLIB

#define PREFIX_SPECDATA "SPECDATA"    // e.g., specdata00, specdata01, ...
#define VARNAME_SPECDATA_HOSTLIB "IDSPECDATA"   // extra prefix for HOSTLIB

#define ITABLE_SPECBASIS 0
#define ITABLE_SPECDATA  1

struct {
  int  ITABLE ;       // either SPECBASIS or SPECDATA (Feb 23 2021)
  char TABLENAME[12] ;   // "BASIS" or "DATA"

  int  NSPECBASIS ; 
  int  NSPECDATA, IDSPECDATA  ; 

  int  NBIN_WAVE;  // number of wavelength bins
  int  ICOL_WAVE;  // table column with wavelength
  int  ICOL_SPECTABLE[MXSPECBASIS_HOSTLIB]; // colum for each template
  int  NUM_SPECBASIS[MXSPECBASIS_HOSTLIB];  // number for each template
  char VARNAME_SPECBASIS[MXSPECBASIS_HOSTLIB][28];

  int  IVAR_HOSTLIB[MXSPECBASIS_HOSTLIB]; // identified HOSTLIB ivar with coeff

  double  FLAM_SCALE, FLAM_SCALE_POWZ1 ;
  double *WAVE_CEN, *WAVE_MIN, *WAVE_MAX, *WAVE_BINSIZE ; // rest-frame
  double *FLAM_BASIS[MXSPECBASIS_HOSTLIB];
  
  int NWARN_INTEG_HOSTMAG[MXFILTINDX];

  double *FLAM_EVT; // updated each event.

} HOSTSPEC ;


typedef struct {
  char VARNAMES_APPEND[MXFILTINDX*5]; // allow "obs_X " per band

  int  NLINE_COMMENT ;
  char *COMMENT[MXFILTINDX];
  char FILENAME_SUFFIX[40];
  int  NLINE_APPEND;
  char **LINE_APPEND ;
} HOSTLIB_APPEND_DEF ;

time_t TIME_INIT_HOSTLIB[2];

// =====================================
//
// prototype declarations (moved from sntools_host.c, Feb 2013)
//
// =====================================

void   INIT_HOSTLIB(void);  // one-time init
void   init_SNHOSTGAL(void);  // init each event
void   GEN_SNHOST_DRIVER(double ZGEN_HELIO, double PEAKMJD);
void   GEN_SNHOST_GALID(double ZGEN);
void   GEN_SNHOST_POS(int IGAL);
void   SIMLIB_SNHOST_POS(int IGAL, SERSIC_DEF *SERSIC, int DEBUG_MODE);
void   GEN_SNHOST_ANGLE(double a, double b, double *ANGLE);
void   GEN_SNHOST_NBR(int IGAL);
void   GEN_SNHOST_DDLR(int i_nbr);
void   SORT_SNHOST_byDDLR(void);
void   reset_SNHOSTGAL_DDLR_SORT(int MAXNBR);

void   TRANSFER_SNHOST_REDSHIFT(int IGAL);
void   GEN_SNHOST_GALMAG(int IGAL);
void   GEN_SNHOST_ZPHOT(int IGAL);
double GEN_SNHOST_ZPHOT_QUANTILE(int IGAL, int q);
void   GEN_SNHOST_VPEC(int IGAL);
void   GEN_SNHOST_LOGMASS(void); // Feb 2020
void   GEN_SNHOST_PROPERTY(int ivar_property); 
int    USEHOST_GALID(int IGAL) ;
void   FREEHOST_GALID(int IGAL) ;
void   checkAbort_noHOSTLIB(void) ;
void   checkAbort_HOSTLIB(void) ;

void   STORE_SNHOST_MISC(int IGAL, int ibin_SNVAR);
double modelPar_from_SNHOST(double parVal_orig, char *parName);
void   DUMPROW_SNHOST(void) ;
void   DUMP_SNHOST(void);
void   initvar_HOSTLIB(void);
void   init_OPTIONAL_HOSTVAR(void) ;
void   init_OPTIONAL_HOSTVAR_PROPERTY(char *basename, int *NVAR_PROPERTY) ;

void   init_REQUIRED_HOSTVAR(void) ;
int    load_VARNAME_STORE(char *varName) ;
void   open_HOSTLIB(FILE **fp);
void   close_HOSTLIB(FILE *fp);

void   init_HOSTLIB_WGTMAP(int OPT_INIT, int IGAL_START, int IGAL_END);
void   read_HOSTLIB_WGTMAP(void);
void   parse_HOSTLIB_WGTMAP(FILE *fp, char *string);
int    read_VARNAMES_WGTMAP(char *VARLIST);
bool   checkSNvar_HOSTLIB_WGTMAP(char *varName);
void   runCheck_HOSTLIB_WGTMAP(void);
void   malloc_HOSTLIB_WGTMAP(void); 
void   malloc_HOSTGAL_PROPERTY(void);
int    getindex_HOSTGAL_PROPERTY(char *PROPERTY);

// xxxvoid   malloc_SNVAR_HOSTLIB_WGTMAP(int NGAL, int NBTOT, double ***PTR);
void   prep_SNVAR_HOSTLIB_WGTMAP(void);
void   getVal_SNVAR_HOSTLIB_WGTMAP(int ibin, double *VAL_WGTMAP); // init
int    getBin_SNVAR_HOSTLIB_WGTMAP(void); // for each event

void   parse_Sersic_n_fixed(FILE *fp, char *string); 
void   read_head_HOSTLIB(FILE *fp);
bool   match_varname_HOSTLIB(char *varName0, char *varName1);
void   checkAlternateVarNames_HOSTLIB(char *varName) ;
void   read_gal_HOSTLIB(FILE *fp);
void   read_galRow_HOSTLIB(FILE *fp, int nval, double *values, 
			   char *field, char *nbr_list  );
int    passCuts_HOSTLIB(double *xval);
void   summary_snpar_HOSTLIB(void) ;
void   malloc_HOSTLIB(int NGAL_STORE, int NGAL_READ);
void   sortz_HOSTLIB(void);
void   zptr_HOSTLIB(void);
void   init_HOSTLIB_ZPHOTEFF(void);
void   init_HOSTLIB_ZPHOT_QUANTILE(void);
void   init_GALMAG_HOSTLIB(void);
void   init_Gauss2d_Overlap(void);
void   init_SAMEHOST(void);
void   init_Sersic_VARNAMES(void);
void   init_Sersic_HOSTLIB(void);
void   init_Sersic_integrals(int j);
void   read_Sersic_bn(void);
void   Sersic_names(int j, char *a, char *b, char *w, char *n);
void   get_Sersic_info(int IGAL, SERSIC_DEF *SERSIC) ;
void   test_Sersic_interp(void);
double get_Sersic_bn(double n);
void   init_OUTVAR_HOSTLIB(void) ;
void   LOAD_OUTVAR_HOSTLIB(int IGAL) ;
void   append_HOSTLIB_STOREPAR(void);
bool   QstringMatch(char *varName0, char *varName1);

void   readme_HOSTLIB(void);
void   check_duplicate_GALID(void);
int    IVAR_HOSTLIB(char *varname, int ABORTFLAG);
int    IVAR_HOSTLIB_PREFIX(char *prefix, int ABORTFLAG);
bool   ISCHAR_HOSTLIB(int IVAR);

long long get_GALID_HOSTLIB(int igal);
double get_ZTRUE_HOSTLIB(int igal);
double get_VALUE_HOSTLIB(int ivar, int igal);
double get_GALFLUX_HOSTLIB(double a, double b);

double interp_GALMAG_HOSTLIB(int ifilt_obs, double PSF ); 
double Gauss2d_Overlap(double offset, double sig);
void   magkey_HOSTLIB(int  ifilt_obs, char *key, char *key_err);
void   set_usebit_HOSTLIB_MSKOPT(int MSKOPT);

void setbit_HOSTLIB_MSKOPT(int MSKOPT) ; // added Jan 2017

void GEN_SNHOST_ZPHOT_from_CALC(double ZGEN, double *ZPHOT, double *ZPHOT_ERR);
void zphoterr_asym(double ZTRUE, double ZPHOTERR, 
		   GENGAUSS_ASYM_DEF *asymGaussPar );

void GEN_SNHOST_ZPHOT_from_HOSTLIB(int INBR, double ZGEN, 
				   double *ZPHOT, double *ZPHOT_ERR); 
double snmagshift_salt2gamma_HOSTLIB(int GALID);

void   set_GALID_UNIQUE(int i);

// SPECBASIS functions
void   read_specTable_HOSTLIB(void);
void   match_specTable_HOSTVAR(void);
void   checkVarName_specbasis(char *varName);
int    ICOL_SPECTABLE(char *varname, int ABORTFLAG) ;
void   genSpec_HOSTLIB(double zhel, double MWEBV, int DUMPFLAG,
		       double *GENFLUX_LIST, double *GENMAG_LIST);

// fetch_HOSTPAR function for GENMODEL (e.g., BYOSED)
int fetch_HOSTPAR_GENMODEL(int OPT, char *NAMES_HOSTPAR, double *VAL_HOSTPAR);

void   rewrite_HOSTLIB(HOSTLIB_APPEND_DEF *HOSTLIB_APPEND);
void   malloc_HOSTLIB_APPEND(int NGAL, HOSTLIB_APPEND_DEF *HOSTLIB_APPEND);
void   addComment_HOSTLIB_APPEND(char *COMMENT,
				 HOSTLIB_APPEND_DEF *HOSTLIB_APPEND);
void   rewrite_HOSTLIB_plusNbr(void) ;
void   get_LINE_APPEND_HOSTLIB_plusNbr(int igal_unsort, char *LINE_APPEND);
void   rewrite_HOSTLIB_plusMags(void);
void   monitor_HOSTLIB_plusNbr(int OPT, HOSTLIB_APPEND_DEF *HOSTLIB_APPEND);

void   rewrite_HOSTLIB_plusAppend(char *append_file);

double integmag_hostSpec(int IFILT_OBS, double z, int DUMPFLAG);

// END
