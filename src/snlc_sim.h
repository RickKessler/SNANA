/*******************************************
  snlc_sim.h :  Created Nov 2006 by R.Kessler

  HISTORY


 Feb 23 2016: add INPUTS.GENGAUPEAK_AV
 Apr 14 2016: SIMLIB_MXGEN_LIBID-> 1000 (was 10,000)
 Nov 27 2016: MXREAD_SIMLIB->12,000 (was 6000)
 Mar 28 2017: add GENLC.IDSURVEY
 Apr 20 2017: MXPAR_ZVAR -> 100 (was 50)
 Aug 08 2017: remove all FILTSUM variables and functions
 Aug 10 2017: MXEPSIM_PERFILT -> 500 (was 200)
 Aug 14 2017: MXCID_SIM->199 million (was 99 million)
 Sep 25 2017: MXCID_SIM->299 million
 Jan 25 2018:
     + MXFIELD_OVP_SIMLIB -> 10 (was 5); thanks Zoheyr !
     + MXREAD_SIMLIB      -> 100000 (was 12000)

 Jul 30 2018: define input_file_include2
 Jan 06 2020: genmag8 -> genmag, same for epoch8 & peakmag8
 Mar 20 2020: Dillon: MXPAR_ZVAR -> 150 (was 100)

 Nov 05 2020: MXEPSIM_PERFILT -> 1000 (was 500) for ZTF sims.
 Feb 05 2021:
   + remove definition of MXFIELD_OVP_SIMLIB
   + replace MXFIELD_OVP_SIMLIB with MXFIELD_OVP from sndata.h

 Apr 27 2021: MXSIMGEN_DUMP -> 600 (was 400)

 Jun 02 2021:
    + MXEPSIM_PERFILT -> 2000 (was 1000) for 0.1 day cadence
    + MXOBS_SIMLIB -> 15,000 (was 10k) for Roman synthetic bands

 Jun 25 2021: MXINPUT_FILE_SIM -> 4 (was 3)
 Jan 28 2022: MXEPSIM -> 15k (was 10k)

********************************************/


// ************ GLOBAL VARIABLES *************

#define  MXINPUT_FILE_SIM   4    // 1 input file + 3 includes
#define  MXCID_SIM  299999999   // max sim CID and max number of SN
#define  MXEPSIM_PERFILT  2500       //
#define  MXEPSIM       15000  // really big for sntools_grid
#define  MXLAMSIM      4000   // mx number of lambda bins
#define  MXCUTWIN       20
#define  MXCUTWIN_SNRMAX 5    // mx number of SNRMAX cuts
#define  MXCUTWIN_PEAKMJD_BYFIELD 10
#define  MXNON1A_TYPE 1000     // max number of non1a types/indices
#define  MXNON1A_KEY  10      // max number of non1a keys
#define  MXCHAR_FIELDNAME 20
#define  MXSIMGEN_DUMP 600    // max number of variables to dump
#define  TABLEID_DUMP  7100   // for SNTABLE functions
#define  MXZRAN        10     // max randoms to store for z-smearing
#define  MXPAR_SIMSED  30     // max number of SIMSED params

#define  MXREAD_SIMLIB 100000  // max number of SIMLIB observations/entries
#define  MXOBS_SIMLIB  15000    // max number of observ. per simlib
#define  MXOBS_SPECTROGRAPH 50 // max number of spectra per event

#define  MXGENSKIP_PEAKMJD_SIMLIB  10
#define  MXSEASON_SIMLIB  20      // max number of seasons
#define  MXFLUXERR_COR_SIMLIB 100  // max number of FLUXERR_COR keys in header
#define  TGAP_SEASON_SIMLIB 90.0  // gap (days) to define new season
#define  GENRANGE_TOBS_PAD  0.1   // Tobs padding for MJD-selection
#define  SIMLIB_ID_REWIND -7   // rewind flag
#define  ISOURCE_PEAKMJD_RANDOM 1 // PEAKMJD is randomly generated
#define  ISOURCE_PEAKMJD_SIMLIB 2 // PEAKMJD is read from SIMLIB header

#define  OPTLINE_SIMLIB_S             1   // is a SIMLIB line with 'S:'
#define  OPTLINE_SIMLIB_T             2   // obsolete
#define  OPTLINE_SIMLIB_SPECTROGRAPH  3

#define  INDEX_RATEMODEL_POWERLAW  1
#define  INDEX_RATEMODEL_POWERLAW2 2
#define  INDEX_RATEMODEL_AB       3
#define  INDEX_RATEMODEL_ZPOLY    4  // cubic polynomial in redshift
#define  INDEX_RATEMODEL_FLAT     5  //
#define  INDEX_RATEMODEL_CCS15    6  // CANDELES/Strolger 2015 CC rate
#define  INDEX_RATEMODEL_MD14     7  // Madau & Dickey 2014
#define  INDEX_RATEMODEL_PISN     8  //
#define  INDEX_RATEMODEL_TDE      9  // Tidal disrupt event (was EPM)
#define  INDEX_RATEMODEL_BPOLY    10  // relRate vs. GalLat
#define  INDEX_RATEMODEL_COSBPOLY 11  // relRate vs. cos(GalLat)
#define  MXPOLY_GALRATE            5   // 5th order poly for BPOLY or COSB

#define  MXRATEPAR_ZRANGE  9    // max number of z-ranges for rate model
#define  RATEMODELNAME_CCS15       "CC_S15"
#define  RATEMODELNAME_MD14        "MD14"
#define  RATEMODELNAME_PISN        "PISN_PLK12"
#define  RATEMODELNAME_TDE         "TDE"

#define  REDSHIFT_FLAG_NONE      0
#define  REDSHIFT_FLAG_SNSPEC    1
#define  REDSHIFT_FLAG_HOSTSPEC  2
#define  REDSHIFT_FLAG_HOSTPHOT  3
#define  REDSHIFT_FLAG_WRONGHOST 4
#define  STRING_REDSHIFT_FLAG (char*[5]){ "NONE", "SNSPEC", "HOSTSPEC", "HOSTPHOT", "WRONGHOST" }

// define write flags
#define WRMASK_VBOSE     1   // legacy verbose output
#define WRMASK_TEXT      2   // terse/text output: MJD FLUX etc ...
#define WRMASK_MODEL     4   // dump obs model mag vs. Tobs
#define WRMASK_BLINDTEST 8   // suppress SIM_XXX and other info
#define WRMASK_CIDRAN   16   // use random CID (1-MXCID)
#define WRMASK_FITS     32   // write to fits file instead of ascii
#define WRMASK_COMPACT  64   // suppress non-essential PHOT output
#define WRMASK_FILTERS  256  // write filterTrans files (Aug 2016)

// xxx #define KEYSOURCE_FILE 1
// xxx #define KEYSOURCE_ARG  2
#define FLAG_NWD_ZERO 100 // flag that override word is a key with no arg

#define IFLAG_GENSMEAR_FILT 1 // intrinsic smear at central LAMBDA of filter
#define IFLAG_GENSMEAR_LAM  2 // intrinsic smear vs. wavelength
int IFLAG_GENSMEAR ;


#define SMEARMASK_HOSTGAL_PHOT   1 // add Poisson shot noise
#define SMEARMASK_HOSTGAL_IMAGE  2 // anomalouse image noise in HOSTNOISE_FILE

int WRFLAG_VBOSE     ;
int WRFLAG_TEXT      ;
int WRFLAG_MODEL     ;
int WRFLAG_BLINDTEST ;
int WRFLAG_CIDRAN    ;
int WRFLAG_FITS      ;
int WRFLAG_FILTERS   ; // Aug 2016
int WRFLAG_COMPACT   ; // Jan 2018

#define SIMLIB_PSF_PIXEL_SIGMA   "PIXEL_SIGMA"        // default
#define SIMLIB_PSF_ARCSEC_FWHM   "ARCSEC_FWHM"        // option
#define SIMLIB_PSF_NEA_PIXEL     "NEA_PIXEL"          // option
#define SIMLIB_PSF_NEA_ARCSECSQ  "NEA_ARCSECSQ"       // option

#define SIMLIB_SKYSIG_SQPIX    "ADU_PER_SQPIXEL"    // default
#define SIMLIB_SKYSIG_SQASEC   "ADU_PER_SQARCSEC"   // option

//#define SIMLIB_MXGEN_LIBID 10000  // stop gen on this many for one LIBID
#define SIMLIB_MXGEN_LIBID 1000
#define SIMLIB_MSKOPT_REPEAT_UNTIL_ACCEPT      2 // force each LIBID to accept
#define SIMLIB_MSKOPT_QUIT_NOREWIND            4 // quit after one pass
#define SIMLIB_MSKOPT_RANDOM_TEMPLATENOISE     8 // random template noise
#define SIMLIB_MSKOPT_IGNORE_TEMPLATENOISE    16 // ignore template coherence
#define SIMLIB_MSKOPT_IGNORE_FLUXERR_COR      32 // ignore FLUXERR_COR map
#define SIMLIB_MSKOPT_ENTIRE_SEASON          128 // keep entire SIMLIB season
#define SIMLIB_MSKOPT_ENTIRE_SURVEY          256 // keep entire SIMLIB survey

#define METHOD_TYPE_SPEC 1    // spec id
#define METHOD_TYPE_PHOT 2    // phot id
#define OFFSET_TYPE_PHOT 100   // add this to specType to get photoType

FILE  *fp_SIMLIB ;

struct {
  time_t t_start, t_end, t_end_init, t_update_last ;
  int    NGENTOT_LAST ;
} TIMERS ;

// define auxillary files produced with data files.
typedef struct {

  // required outputs
  FILE *FP_LIST;    char  LIST[MXPATHLEN] ;
  FILE *FP_README;  char  README[MXPATHLEN] ;
  FILE *FP_IGNORE;  char  IGNORE[MXPATHLEN] ;
  FILE *FP_DUMP;    char  DUMP[MXPATHLEN] ;
  FILE *FP_YAML;    char  YAML[MXPATHLEN] ;  // Aug 10 2020, for batch mode only
  char PATH_FILTERS[MXPATHLEN]; // directory instead of file

  // optional outputs (just filename, not pointer)
  char  ZVAR[MXPATHLEN] ;   // optional zvariation file
  char  GRIDGEN[MXPATHLEN]; // optional GRID-output

  // char string to write each line to memory, and then
  // one ASCI write per line instead of per value.
  char *OUTLINE ;

} SIMFILE_AUX_DEF ;


// Mar 2016: create typedefs for NON1A
typedef struct {

  int IFLAG_GEN;

  // NON1ASED inputs from sim-input file
  char  PATH[MXPATHLEN];                // user-defined path to NON1A seds
  char  LISTFILE[MXPATHLEN];           // $PATH/NON1A.LIST
  int   NKEY ;                         // number of user NON1A keys
  char  KEYLIST[MXNON1A_KEY][40];      // list of user NON1A keys

  int   INDEX[MXNON1A_TYPE];      // NON1ASED index list vs. sparse index
  float WGT[MXNON1A_TYPE];       // relative wgt for simulation
  float MAGOFF[MXNON1A_TYPE];    // magoff
  float MAGSMEAR[MXNON1A_TYPE][2];  // mag-smear (bifurc Gaussian)
  int   SNTAG[MXNON1A_TYPE];     // user index-tag
  float KEYVAL[MXNON1A_TYPE][MXNON1A_KEY]; // same as above
  int   ISPEC1A[MXNON1A_TYPE];  // 1 --> PEC1A instead of CC (Aug 2016)
  int   NNON1A ;                // number of NON1A SED templates
  int   NPEC1A ;                // number of PEC1A SED templates
  int   STOP ;                  //  flag to stop reading NON1a

  // computed from NON1ASED inputs
  int   NINDEX;                   // number of NON1A indices to simulate
  int   NGEN[MXNON1A_TYPE] ;      // NGEN per NON1ASED index
  int   MXGEN[MXNON1A_TYPE] ;     // MAX NGEN per NON1ASED index
  char  SED_FILE[MXNON1A_TYPE][200] ;   // SED  file per NON1ASED index

  // copied from INPUTS
  int   NGENTOT ;           // copied from INPUTS.NGEN
  int   CIDOFF ;            // copied from INPUTS.CIDOFF

  // read from NON1A.LIST:
  int   LIST_INDEX[MXNON1A_TYPE];     // sparse list of valid indices
  char  LIST_TYPE[MXNON1A_TYPE][20];  // 'II', 'Ib', '91bg', etc ...
  int   LIST_ISPEC1A[MXNON1A_TYPE];   // 1 --> PEC1A (Aug 2016)
  char  LIST_NAME[MXNON1A_TYPE][MXPATHLEN];  // name of entry
  int   NLIST ;
  int    INDEXVALID[MXNON1A_TYPE];   // index is argument
  double FLUXSCALE[MXNON1A_TYPE];       // 5.2019
  double RESTLAMBDA_RANGE[2];           // 5.2019

} INPUTS_NON1ASED_DEF ;


// parameters to define random systematic shifts at start of sim;
// thus changing RANSEED alone results in random systematic shifts
// for all parameters without user specifying each random shift externally.
typedef struct {
  int   USE ;

  int  IDSURVEY;    // override true IDSURVEY for survey-dependent randoms

  int  RANSEED_SYST; // RANSEED for picking syst variations (if >0)
  int  RANSEED_GEN;  // ranseed for AFTER picking syst params (if >0)

  float SIGSHIFT_ZP[MXFILTINDX];
  float SIGSHIFT_LAMFILT[MXFILTINDX]; // filterTrans shifts, Ang
  float SIGSCALE_FLUXERR;    // scale true & measured errors by 1+Gran*SIG
  float SIGSCALE_FLUXERR2;   // scale measured errors, not true errors
  float SIGSCALE_MWEBV;      // scale Galactic extinction by 1+Gran*SIG
  float SIGSHIFT_MWRV;       // shift RV
  float SIGSHIFT_REDSHIFT;   // shift measured redshift PA 2020
  float SIGSHIFT_zPHOT_HOST; // shift in host photo-z, May 2022
  char GENMODEL_WILDCARD[MXPATHLEN]; // choose between wildcard models PA 2020
  char GENPDF_FILE_WILDCARD[MXPATHLEN]; // choose between wildcard GENPDF_FILEs

  float SIGSHIFT_OMEGA_MATTER ;
  float SIGSHIFT_W0 ;
  float RANGESHIFT_OMEGA_MATTER[2] ;
  float RANGESHIFT_W0[2] ;

} INPUTS_RANSYSTPAR_DEF ;


typedef struct  {
  int   IFLAG_GEN ;               // indicates RANDOM or GRID
  int   ISPARSE;                  // current sparse index
  int   NGENWR[MXNON1A_TYPE];     // actual # written per index
  int   NGENTOT[MXNON1A_TYPE];    // actual # generated per index
  int   CIDRANGE[MXNON1A_TYPE][2]; // CID range for each non1a index
  char  TYPE[MXNON1A_TYPE][40];    // char def such as II, Ib, Ic, ..
  float FRAC_PEC1A ;              // computed N(pec1a)/NTOT(NON1A+Pec1a)
} GENLC_NON1ASED_DEF ;


typedef struct {
  char   NAME[40] ;           // filled internally

  double DNDZ_ZEXP_REWGT;     // re-wgt dN/dz by z^ZEXP_REWGT
  double DNDZ_ZPOLY_REWGT_LEGACY[4]; // xxx legacy var
  GENPOLY_DEF DNDZ_ZPOLY_REWGT ;     // poly(z) to reweight rate-vs-z

  double DNDZ_SCALE[2] ;      // scale DNDZ for Ia and NON1A (4/19/2017)
  double DNDZ_ALLSCALE ;      // scale all SN models (Ia, SIMSED, etc ... )
  double DNDB_SCALE ;         // scale rate for galactic models (LCLIB)

  // rate model pars: A&B or R0&beta, zpoly ...
  int    NMODEL_ZRANGE;     // number of rate models glued together
  double MODEL_PARLIST[MXRATEPAR_ZRANGE+1][20];
  double MODEL_ZRANGE[MXRATEPAR_ZRANGE+1][2];  // z-range of rate model
  GENPOLY_DEF MODEL_ZPOLY;  // some models are polynom functions of z
  GENPOLY_DEF MODEL_BPOLY;  // rate vs. Galact b or cosb (LCLIB model)
  int    INDEX_MODEL ;  //

  double RATEMAX;  // used with DNDB (Galactic models)

  // max redshift wgt for generation in genz_hubble()
  double  ZGENWGT_MAX, ZGENMIN_STORE, ZGENMAX_STORE ;

  // predicted SN count
  double SEASON_COUNT ;     // nominal SN count per season
  double SEASON_FRAC ;      // fracion among RATEPARs (SN,PEC1A)

} RATEPAR_DEF ;


#define SPECTROGRAPH_OPTMASK_LAMSIGMAx0  1  // no lam smearing
#define SPECTROGRAPH_OPTMASK_LAMSIGMAx2  2  // double LAMSMEAR
#define SPECTROGRAPH_OPTMASK_LAMSPIKE    4  // flux only in center bin
#define SPECTROGRAPH_OPTMASK_SNRx100     8  // multiply SNR x 100
#define SPECTROGRAPH_OPTMASK_noTEMPLATE 16  // no template noise
#define SPECTROGRAPH_OPTMASK_onlyTNOISE 32  // only template noise
#define SPECTROGRAPH_OPTMASK_TEXTRAP    64  // extrap TEXPOSE outside range
#define SPECTROGRAPH_OPTMASK_NOSPEC   2048  // skip spectra
#define SPECTROGRAPH_OPTMASK_FITS_REFAC 1024 // temp: refactor fits format
#define SPECTROGRAPH_OPTMASK_FITS_LEGACY 4096 // later: preserve legacy format
#define SPECTROGRAPH_OPTMASK_noNOISE 32768  // internal only: turn off noise

typedef struct {
  int    DOFLAG_SPEC ; // logical flag for spectra
  int    OPTMASK;    // sim-input key SPECTROGRAPH_OPTMASK: <MASK>
  double NLAMSIGMA ; // how far out to smear flux in lambda bins

  // below are for tests & debugging, based on OPTMASK
  int    ILAM_SPIKE ;     // set by OPTMASK
  double SCALE_LAMSIGMA ; // set by OPTMASK
  double SCALE_SNR ;      // set by OPTMASK
  double SCALE_TEXPOSE ;  // from user input SPECTROGRAPH_SCALE_TEXPOSE

} SPECTROGRAPH_OPTIONS_DEF;


#define MXPEREVT_TAKE_SPECTRUM MXSPECTRA
int     NPEREVT_TAKE_SPECTRUM ;
typedef struct {
  float   EPOCH_RANGE[4];    // Trest or TOBS range, or MJD range
  char    FIELD[80];         // restrict spectra to particular field(s)

  GENPOLY_DEF GENLAMPOLY_WARP ;   // calibration warp as poly fun of wavelength
  GENPOLY_DEF GENZPOLY_TEXPOSE ;  // TEXPOSE =poly fun of z
  GENPOLY_DEF GENZPOLY_SNR ;      // SNR = poly fun of z
  float   SNR_LAMRANGE[2];   // lam-range to define SNR
  char    EPOCH_FRAME[8];    // either 'REST' or 'OBS' or 'HOST'
  int   OPT_FRAME_EPOCH  ; // epoch is GENFRAME_REST or GENFRAME_OBS or MJD
  int   OPT_FRAME_LAMBDA ; // for SNR opt below, LAMREST or LAMOBS

  // OPT_TEXPOSE = 1(TEXPOSE), or 2(SNR)
  int   OPT_TEXPOSE ;

} TAKE_SPECTRUM_DEF;

#define NTYPE_FLUXNOISE      6
#define TYPE_FLUXNOISE_S     0  // from search image -> feeds pipe effic.
#define TYPE_FLUXNOISE_SZ    1  // Search plus zero point
#define TYPE_FLUXNOISE_Z     2  // from zero point
#define TYPE_FLUXNOISE_T     3  // from template
#define TYPE_FLUXNOISE_F     4  // from error fudge
#define TYPE_FLUXNOISE_SUM   5  // quadrature sum S+Z+T+F

typedef struct {

  // all SQSIG are in p.e.
  double SQSIG_CALC_TRUE[NTYPE_FLUXNOISE];   // actual scatter
  double SQSIG_CALC_DATA ;
  double SQSIG_FUDGE_TRUE[NTYPE_FLUXNOISE];
  double SQSIG_FUDGE_DATA ;
  double SQSIG_FINAL_TRUE[NTYPE_FLUXNOISE];
  double SQSIG_FINAL_DATA, SIG_FINAL_DATA ;

  double SNR_CALC_S;     // for legacy fluxerr scale option
  double SNR_CALC_ST;    // feeds pipe detection effic
  double SNR_CALC_SZT ;  // feeds FLUXERRMODEL_FILE

  double SQSIG_SRC;       // source noise = Npe
  double SQSIG_SKY;       // sky+CCD noise
  double SQSIG_ZP ;        // ZP error
  double SQSIG_TSRC;      // template source noise
  double SQSIG_TSKY ;     // corrsky_ccd  noise from template
  double SQSIG_HOST_PHOT ;  // galaxy shot noise
  double SQSIG_HOST_IMAGE ; // anomalous noise from HOSTNOISE_FILE->obsolete
  double SQSIG_RAN ;        // error shift due to Poisson noise

  // SNR for monitor
  double SNR_CALC_MON ;   // calc SNR for MAGMONITOR_SNR input
  double SNR_FINAL_MON ;   // calc SNR for MAGMONITOR_SNR input

  // misc.
  double Npe_over_FLUXCAL, NADU_over_Npe, NEA, GALMAG_NEA ;
  char BAND[2];
  int  IFILT_OBS, INDEX_REDCOV ;

  // store actual flux-shift for each term.
  double FLUX_SHIFT_TRUE[NTYPE_FLUXNOISE];

} FLUXNOISE_DEF ;


typedef struct {
  int    NSUM ;
  double RHO_EVT, RHO_SUM, RHO_AVG ;
} MONITOR_REDCOV_FLUXNOISE_DEF ;


// Nov 2021: define struct for interpolating host photo-z resolution vs z
//   Sim input key is HOSTLIB_GENZPHOT_FUDGEMAP: <STRING>
typedef struct {
  char STRING[200]; // e.g., 'z0:RMS0,z1:RMS1,z2:RMS2,etc...'
  int    NzBIN; 
  double z_LIST[100]  ;
  double RMS_LIST[100] ;
} HOSTLIB_GENZPHOT_FUDGEMAP_DEF ;

// -------------------------------------
// define user INPUTS

struct INPUTS {

  int USE_KCOR_REFACTOR; //1-> run both legacy and new; 2-> new only
  int USE_KCOR_LEGACY;   //use legacy fortran code to read & apply

  bool DASHBOARD_DUMPFLAG ;  // dump all input maps and libraries
  bool KEYNAME_DUMPFLAG;     // dump input key names and quit (broken!!)
  bool README_DUMPFLAG;      // dump readme and stop (Feb 2022)

  // input file list includes nominal, plus up to few INCLUDE files
  char INPUT_FILE_LIST[MXINPUT_FILE_SIM][MXPATHLEN]; // input file names
  char TIME_START[100];  // override start time for batch uniformity 

  int  NREAD_INPUT_FILE;  // number of input files read: 1,2 or 3

  int  NWORDLIST ;      // number of words read from input file
  char **WORDLIST ;     // list of words read from input file

  int  TRACE_MAIN;            // debug to trace progress through main loop
  int  DEBUG_FLAG ;           // arbitrary debug usage

  bool RESTORE_DES3YR;          // restore DES3YR bugs
  bool RESTORE_HOSTLIB_BUGS ;   // set if DEBUG_FLAG==3 .or. RESTORE_DES3YR
  bool RESTORE_FLUXERR_BUGS ;   // set if DEBUG_FLAG==3 .or. idem
  bool RESTORE_WRONG_VPEC   ;   // restore incorrect VPEC sign convention

  char SIMLIB_FILE[MXPATHLEN];  // read conditions from simlib file
  char SIMLIB_OPENFILE[MXPATHLEN];  // name of opened files (internal)
  int  SIMLIB_GZIPFLAG ;            // gzip flag (needed to rewind)

  char SIMLIB_SURVEY[40];     // override name of SURVEY in simlib file
  char SIMLIB_FIELDLIST[200]; // default=ALL, or, e.g., C1+C2+C3
  int  SIMLIB_FIELDSKIP_FLAG ; // INTERNAL: 1->count skipped fields for NGENTOT
  STRING_DICT_DEF DICT_SIMLIB_FIELDLIST_PRESCALE;   // SIMLIB ps per FIELD
  STRING_DICT_DEF DICT_SPECTRUM_FIELDLIST_PRESCALE; // spectrum ps per FIELD

  int  SIMLIB_IDSTART;      // start at this LIBID (default=1)
  int  SIMLIB_MAXRANSTART;  // start at random LIBID among this many
  int  SIMLIB_IDLOCK;    // >0 => use same LIBID;
                         // e.g., 1 => always 1st used LIBID
  int  SIMLIB_MINOBS;    // min NOBS to use (default=5)
  int  SIMLIB_NREPEAT ;  // repeat each ID this many times (for less reading)
  int  SIMLIB_MXREPEAT ; // used only with BPOLY Galactic rate model
  double SIMLIB_MINSEASON ; // min season length (days); default=0

  int    SIMLIB_NSKIPMJD ;  // number of MJDs to skip to change cadence
  int    SIMLIB_IDSKIP[MXREAD_SIMLIB]; // list of SIMLIB IDs to skip
  int    NSKIP_SIMLIB ;       // number of SIMLIB_IDSKIP values read

  int    SIMLIB_DUMP;  // dump this simlib id, then quit (0=all)
  float  SIMLIB_CADENCEFOM_ANGSEP; // controls calc of cadence FoM
  double SIMLIB_CADENCEFOM_PARLIST[10] ; // optional *parList for SNcadenceFoM

  int  USE_SIMLIB_GENOPT ;    // use some optional gen-keys in simlib header
  int  USE_SIMLIB_REDSHIFT ;  // 1 => use redshift in LIB (if it's there)
  int  USE_SIMLIB_DISTANCE ;  // 1 => use distance in LIB (if it's there)
  int  USE_SIMLIB_PEAKMJD ;   // idem for optional PEAKMJD
  int  USE_SIMLIB_MAGOBS ;    // use MAGOBS column instead of SN model
  int  USE_SIMLIB_SPECTRA;    // use TAKE_SPECTRUM keys in SIMLIB header
  int  USE_SIMLIB_SALT2 ;     // use SALT2c and SALT2x1 from SIMLIB header
  int  SIMLIB_MSKOPT ;        // special SIMLIB options (see manaul)

  // ---- end simlib inputs -----

  int  CIDOFF ;             // user CID offset
  int *CIDRAN_LIST ;        // for internal use to create random CID list
  int  CIDRAN_MAX ;         // max CID (used for random CID only)
  int  CIDRAN_MIN ;         // min CID (used for random CID only)

  int  JOBID;       // command-line only (for batch) to compute SIMLIB_IDSTART
  int  NJOBTOT;     // idem, for submit_batch_jobs.py

  int  HOSTLIB_USE ;            // 1=> used; 0 => not used, 2=>rewrite HOSTLIB
  char HOSTLIB_PLUS_COMMAND[60];        //e.g., +HOSTMAGS, +HOSTNBR, +HOSTAPPEND
  char HOSTLIB_FILE[MXPATHLEN]; // lib of Ztrue, Zphot, Zerr ...
  char HOSTLIB_APPEND_FILE[MXPATHLEN];  // argument of +APPEND
  char HOSTLIB_WGTMAP_FILE[MXPATHLEN];  // optional wgtmap override
  char HOSTLIB_ZPHOTEFF_FILE[MXPATHLEN];  // optional EFF(zphot) vs. ZTRUE
  char HOSTLIB_SPECBASIS_FILE[MXPATHLEN]; // spec basis vec for host spec
  char HOSTLIB_SPECDATA_FILE[MXPATHLEN]; // spec data for host spec
  int  HOSTLIB_MSKOPT ;         // user bitmask of options
  int  HOSTLIB_MAXREAD ;        // max entries to read (def= infinite)
  int  HOSTLIB_GALID_NULL ;     // value for no galaxy; default is -9
  int  HOSTLIB_GALID_PRIORITY[2] ;  // preferentially select this GALID range
  int  HOSTLIB_GALID_UNIQUE;         // flag to force unique galid
  int  HOSTLIB_MINDAYSEP_SAMEGAL ;    // min DAYs before re-using host gal
  float  HOSTLIB_MNINTFLUX_SNPOS; // gen SNPOS greater than this flux-fraction (.00)
  float  HOSTLIB_MXINTFLUX_SNPOS; // gen SNPOS within this flux-fraction (.99)
  float  HOSTLIB_GENRANGE_NSIGZ[2];  // allowed range of (Zphot-Z)/Zerr
  float  HOSTLIB_MAXDDLR ;             // keep hosts with DDLR < MAXDDLR
  float  HOSTLIB_GENZPHOT_FUDGEPAR[5]; // analytic ZPHOT & ZPHOTERR   

  HOSTLIB_GENZPHOT_FUDGEMAP_DEF HOSTLIB_GENZPHOT_FUDGEMAP;
  float  HOSTLIB_GENZPHOT_OUTLIER[2];  // range for FLAT outlier distribution.
  float  HOSTLIB_GENZPHOT_BIAS[5];     // poly(z) bias on ZPHOT (Mar 28 2018)
  int    USE_HOSTLIB_GENZPHOT;         // T if any of above are used

  double HOSTLIB_GENRANGE_RA[2];
  double HOSTLIB_GENRANGE_DEC[2];
  double HOSTLIB_SBRADIUS ; // arcsec, determine SB using this radius

  double HOSTLIB_DZTOL[3] ; // define zSN-zGAL tol vs z
  GENPOLY_DEF HOSTLIB_GENPOLY_DZTOL; // zSN-zGAL tol vs zPOLY

  double HOSTLIB_SCALE_LOGMASS_ERR ; // default is 1.0
  char   HOSTLIB_SCALE_PROPERTY_ERR[200] ; // e.g. '0.8(LOGMASS),0.0(LOGSFR),0.1(LOGsSFR)' , default is 1.0 for every host prop  
  double HOSTLIB_SCALE_SERSIC_SIZE ; // default is 1.0
  char   HOSTLIB_STOREPAR_LIST[2*MXPATHLEN]; // (I) comma-sep list

  // debug options
  int    HOSTLIB_GALID_FORCE ;    // force this GALID
  double HOSTLIB_ABMAG_FORCE ;    // for ABmag on galmag and gal spec
  double HOSTLIB_FIXRAN_RADIUS ;  // fix random number of radius
  double HOSTLIB_FIXRAN_PHI ;     // fix random number for phi
  double HOSTLIB_FIXSERSIC[4];    // fix sersic a,b,n,a_rot

  char FLUXERRMODEL_FILE[MXPATHLEN];   // input err-scale map(s)
  char FLUXERRMAP_IGNORE_DATAERR[100]; // list of MAPNAMES to ignore in data error
  int  FLUXERRMODEL_OPTMASK ;
  char FLUXERRMODEL_REDCOV[200];  // overwrite REDCOR key in _FILE

  // define anomalous subtraction noise in separate file to be
  // used in both the simulation and in snana to inflate errors.
  char HOSTNOISE_FILE[MXPATHLEN];

  int  NGEN_LC;            // number of LCs to generate & write
  int  NGENTOT_LC;         // number of LCs to generate
  float NGEN_SCALE ;       // scale NGEN_LC or NGENTOT_LC (May 2013)
  float NGEN_SCALE_NON1A;  // scale only NONIA rate
  float NGEN_SEASON ;      // Number of seasons to generate;
                           // season defined by redshift, MJDrange and dOmega
  int  NGEN ;               // one of the above
  int  NGEN_SCREEN_UPDATE ; // screen-update after this many

  int  INIT_ONLY;           // 1=quit after init_RateModel
                            // 2=quit before generation

  int  NSUBSAMPLE_MARK ;    // mark independent sub-samples with index

  char PATH_SNDATA_SIM[MXPATHLEN];
  char GENVERSION[200];      // output version
  char GENPREFIX[200];       // filename prefix (default=GENVERSION)
  char GENSOURCE[20];       // 'RANDOM'  or  'DATA-VERSION'
  char GENMODEL[MXPATHLEN] ; // source model name, with optional path
  char MODELPATH[MXPATHLEN]; // path to model (formerly GENLC.MODELPATH)
  char MODELNAME[100];       // stripped from GENMODEL

  char GENPDF_FILE[MXPATHLEN];   // PDF for color, stretch, etc ...
  char GENPDF_IGNORE[MXPATHLEN];
  int  GENPDF_OPTMASK;           // bit-mask of options
  char GENPDF_FLAT[100];         // force flat PDF; e.g., SALT2x1,SALT2c,RV

  char GENMODEL_EXTRAP_LATETIME[MXPATHLEN];
  char GENSNXT[20] ;        // SN hostgal extinction: CCM89 or SJPAR
  int  GENMODEL_MSKOPT;     // bit-mask of model options
  char GENMODEL_ARGLIST[400] ;
  int  GENMAG_SMEAR_MSKOPT;   // bit-mask of GENSMEAR options
  unsigned int ISEED;         // random seed
  unsigned int ISEED_ORIG;    // for readme output
  int          NSTREAM_RAN;   // number of independent random streams

  int    RANLIST_START_GENSMEAR;  // to pick different genSmear randoms

  double OMEGA_MATTER;   // used to select random Z and SN magnitudes
  double OMEGA_LAMBDA;
  double w0_LAMBDA;
  double wa_LAMBDA;    // w = w0 + wa*(1-a)
  double H0;           // km/s per MPc
  double MUSHIFT;      // coherent MU shift at all redshifts (Oct 2020)
  char   HzFUN_FILE[MXPATHLEN];  // 2 column file with zCMB H(z,theory)
  HzFUN_INFO_DEF HzFUN_INFO;     // store cosmo theory info here.

  double GENRANGE_RA[2];        // RA range (deg) to generate
  double GENRANGE_DEC[2];       // idem for DEC
  float SOLID_ANGLE;           // non-zero => overwrite default calc.
  float MXRADIUS_RANDOM_SHIFT; // random coord shift within MXRADIUS, deg

  double GENRANGE_REDSHIFT[2];  // generated zCMB range
  double GENSIGMA_REDSHIFT;     // smear reported redshift
  double GENBIAS_REDSHIFT ;     // measurement + local bubble
  float  GENSIGMA_VPEC ;        // Gaussian sigma on Vpec, km/sec (default=0)
  double VEL_CMBAPEX ;          // CMB dipole vel; =0 --> zhelio = zcmb
  float  VPEC_ERR;              // error for vpec in data file

  char   WRONGHOST_FILE[MXPATHLEN];

  RATEPAR_DEF RATEPAR ;
  RATEPAR_DEF RATEPAR_PEC1A ; // only for PEC1A in NON1A input

  char  STRING_FUDGE_SNRMAX[20];  // '[SNRMAX]'  or  '[BAND]=[SNRMAX]'
  float FUDGE_SNRMAX;      // adjust EXPOSURE_TIME to force SNR at peak
  int   IFILTOBS_FUDGE_SNRMAX; // -1=all, or get scale from this IFILTOBS
  int   OPT_FUDGE_SNRMAX ;     // 1=adjust EXPOSURE_TIME; 2=adjust sigSKY only

  double GENRANGE_MJD[2];         // range of MJD: allows rigid end
  double GENRANGE_PEAKMJD[2];     // range of PEAKMJD to generate
  double MJD_EXPLODE ;          // define explosion time for NON1A or SIMSED
  double GENRANGE_PEAKMAG[2] ;  // OR among filters (Mar 2016)
  float  GENRANGE_TREST[2];     // relative to peak, days
  float  GENRANGE_TOBS[2];      // for GRID option
  float  GENSIGMA_PEAKMJD;      // option to estimate PEAKMJD for data file

  int    OPT_SETPKMJD;          // opton to estimate PEAKMJD for data file
  float  MJDWIN_SETPKMJD;
  float  SNRCUT_SETPKMJD;
  float  NEWMJD_DIF;            // min MJD dif to count NEW MJD

  int  MAGMONITOR_SNR ;   // compute SNR for this mag -> monitor

  GENGAUSS_ASYM_DEF GENGAUSS_SHAPEPAR ;    // MEAN, SIGMA, RANGE
  GENGAUSS_ASYM_DEF GENGAUSS_DELTA ;
  GENGAUSS_ASYM_DEF GENGAUSS_DM15 ;
  GENGAUSS_ASYM_DEF GENGAUSS_STRETCH ;
  GENGAUSS_ASYM_DEF GENGAUSS_RV ;

  GEN_EXP_HALFGAUSS_DEF  GENPROFILE_AV ;
  GEN_EXP_HALFGAUSS_DEF  GENPROFILE_EBV_HOST ;

  SPECTROGRAPH_OPTIONS_DEF  SPECTROGRAPH_OPTIONS ;
  TAKE_SPECTRUM_DEF         TAKE_SPECTRUM[MXPEREVT_TAKE_SPECTRUM] ;
  float                     TAKE_SPECTRUM_TEMPLATE_TEXPOSE_SCALE ;
  int                       TAKE_SPECTRUM_DUMPCID;
  float                     TAKE_SPECTRUM_HOSTFRAC;
  float                     TAKE_SPECTRUM_HOSTSNFRAC;

  char                      WARP_SPECTRUM_STRING[200];
  int                       NWARP_TAKE_SPECTRUM ; // set internally
  int                       NHOST_TAKE_SPECTRUM ; // set internally

  char  STRETCH_TEMPLATE_FILE[200];

  /* xxxxx mark delete Mar 17 2022 xxxxxxxxxx
     [use struct in sntools_genExpHalfGauss.h]
  double GENRANGE_AV[2];      // host extinction range
  double GENEXPTAU_AV ;       // define exponential distribution of AV
  double GENGAUSIG_AV ;       // AV-sigma of Gaussian core
  double GENGAUPEAK_AV;       // location of Gauss peak (degfault=0)
  double GENRATIO_AV0;        // Expon/Gauss ratio at AV0
  xxxxxxxxxxx end mark xxxxxxx */

  int    DOGEN_AV ;
  int    OPT_SNXT ;  // option for hostgal extinction
  bool   DOGEN_SHAPE, DOGEN_COLOR ; // generate with function or GENPDF

  int    WV07_GENAV_FLAG;   //  1=> use published ESSENCE-WV07 distrib.
  double WV07_REWGT_EXPAV;   //  re-weight exp component (default=1)

  GENGAUSS_ASYM_DEF GENGAUSS_SALT2c ;    // MEAN, SIGMA, RANGE
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2x1 ;
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2ALPHA ;
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2BETA ;
  GENPOLY_DEF SALT2BETA_cPOLY;
  double BIASCOR_SALT2GAMMA_GRID[2]; // gamma range for BBC-biasCor sample
  char   GENPOP_ASYMGAUSS_FILE[MXPATHLEN];   // file with population params
  char   GENPOP_ASYMGAUSS_MODEL[40];         // e.g., G10_HIZ, G10_LOWZ, etc .

  float  GENALPHA_SALT2 ; // legacy variable: same as GENMEAN_SALT2ALPHA
  float  GENBETA_SALT2 ;  // legacy variable: same as GENMEAN_SALT2BETA

  GENGAUSS_ASYM_DEF GENGAUSS_RISETIME_SHIFT ;
  GENGAUSS_ASYM_DEF GENGAUSS_FALLTIME_SHIFT ;

  double  FIXMAG[2] ; // for observer-frame FIXMAG model
  int     GENFRAME_FIXMAG;  // GENFRAME_OBS or GENFRAME_REST

  // SIMSED parameters & ranges
  int   USE_BINARY_SIMSED;  // 1 => use binary files fof faster I/O
  char  PATH_BINARY_SIMSED[MXPATHLEN]; // location of binaries (default = ./)

  int   NPAR_SIMSED_PARAM;     // continuous interp params
  int   NPAR_SIMSED_GRIDONLY;  // GRIDONLY parameters (random)
  int   NPAR_SIMSED_param;     // baggage parameters
  int   IPAR_SIMSED_SHAPE ;    // for SIMSED_SHAPEPAR key
  int   IPAR_SIMSED_COLOR ;    // for SIMSED_COLORPAR key
  int   IPAR_SIMSED_CL ;       // for SIMSED_COLORLAW key

  int   NPAR_SIMSED;           // sum of above
  int   NPAR_SIMSED_MODEL;     // total NPAR in model
  char  KEYWORD_SIMSED[MXPAR_SIMSED][20];  // input keyword per param
  int   GENFLAG_SIMSED[MXPAR_SIMSED];      // bitmask corresponding to above
  char  PARNAME_SIMSED[MXPAR_SIMSED][40];
  int   OPTMASK_SIMSED; // argument for genmag_SIMSED (Dec 2018)

  int   NINDEX_SUBSET_SIMSED_GRIDONLY ;
  int   INDEX_SUBSET_SIMSED_GRIDONLY[MXPAR_SIMSED];

  // inputs for SIMSED covariances among parameters in SED.INFO file
  int   NPAIR_SIMSED_COV ;                 // number of input COV pairs
  float COVPAIR_SIMSED_COV[MXPAR_SIMSED];  // COV among two variables
  int   IPARPAIR_SIMSED_COV[MXPAR_SIMSED][2] ; // IPAR indices of two corr vars

  CHOLESKY_DECOMP_DEF SIMSED_DECOMP;
  int    IPARLIST_SIMSED_COV[MXPAR_SIMSED];  // vs. IROW
  int    IROWLIST_SIMSED_COV[MXPAR_SIMSED];  // vs. IPAR
  int    NROW_SIMSED_COV;

  GENGAUSS_ASYM_DEF GENGAUSS_SIMSED[MXPAR_SIMSED];

  int    MWEBV_FLAG;             // Milky way extinction flag (0,1)
  double MWEBV_SIGRATIO;         // ERR(MWEBV)/MWEBV (=0.16 in SJ code)
  double MWEBV_SIG ;             // fixed error on MWEBV
  double MWEBV_SHIFT;
  double MWEBV_SCALE;            // multiplicative scale (Jun 2017)
  double GENRANGE_MWEBV[2] ;     // useful for effic. maps.
  double RV_MWCOLORLAW ;          // Galactic RV (default=3.1)
  int    OPT_MWCOLORLAW ;        // 89(CCM89), 94(Odonnel), 99(Fitzpat)
  int    OPT_MWEBV ;             // option to modify MWEBV_SFD
  int    APPLYFLAG_MWEBV;        // flag to apply MWEBV analysis corrections
  char   STR_MWCOLORLAW[60] ;    // char-string for comments
  char   STR_MWEBV[60] ;         // char-string for comments

  float GENMAG_OFF_GLOBAL ;            // add same offset to all bands
  float GENMAG_OFF_NON1A ;             // idem, but for NON1A only
  float GENMAG_OFF_ZP[MXFILTINDX];     // add ZP off to gen-mags (calib error)
  float GENMAG_OFF_MODEL[MXFILTINDX];  // add user-model offset to mags

  float TMPOFF_ZP[MXFILTINDX];         // vs. sparse ifilt
  float TMPOFF_MODEL[MXFILTINDX];
  float TMPOFF_LAMFILT[MXFILTINDX];    // vs. sparse index

  // use float for MJD_TEMPLATE to use convenient parsing utility.
  // Precision doesn't matter for this.
  int   USE_MJD_TEMPLATE ;     // logical flag
  float MJD_TEMPLATE;          // Flux(search) -= Flux(MJD_template)
  float MJD_TEMPLATE_FILTER[MXFILTINDX];  // idem, but filter

  int   NFILT_SMEAR;
  int   IFILT_SMEAR[MXFILTINDX];
  float SIGMACLIP_MAGSMEAR[2];  // sigma clipping for mag-smearing

  float  GENMAG_SMEAR_FILTER[MXFILTINDX]; // smear by filter
  float  GENMAG_SMEAR[2];             // intrinsic mag-smear (asymm Gauss)
  char   GENMAG_SMEAR_MODELNAME[100]; // name of specific smear-model
  char   GENMAG_SMEAR_MODELARG[MXPATHLEN];  // optional arg after colon
  float  GENMAG_SMEAR_ADDPHASECOR[2]; // magSmear and coh time
  char   GENMAG_SMEAR_SCALE[100];
  // GENMAG_SMEAR_SCALE(c): 0.98,.1,-0.02,0.3,-5.0
  // GENMAG_SMEAR_SCALE(x1): 0.98,.1,-0.02

  int    NPAR_GENSMEAR_USRFUN ;
  double GENMAG_SMEAR_USRFUN[100];   // intrinsic smear with function

  double GENSMEAR_RANGauss_FIX ;   // if >=0 then set Gauss randoms to this
  double GENSMEAR_RANFlat_FIX ;    // if >=0 then set Flat randoms to this

  char   STRONGLENS_FILE[MXPATHLEN] ;
  char   WEAKLENS_PROBMAP_FILE[MXPATHLEN];
  float  WEAKLENS_DMUSCALE;            // scale width of DMU profile
  float  WEAKLENS_DSIGMADZ ;           // symmetric Gaussian model

  float GENMODEL_ERRSCALE ;    // scale model-errors for intrinsic color-smear
  float GENMODEL_ERRSCALE_CORRELATION; // correlation with GENMAG_SMEAR
  int   GENMODEL_ERRSCALE_OPT ;  // 1=> peak MLCS error; 2=> nominal MLCS err

  int   DO_MODELSMEAR;  // flag to do some kind of model smearing.

  char  GENFILTERS[MXFILTINDX];        // 'gri', 'grizY', etc ...
  int   NFILTDEF_OBS;
  int   IFILTMAP_OBS[MXFILTINDX];     // converts ifilt to ifilt_obs
  float LAMAVG_OBS[MXFILTINDX];
  float LAMRMS_OBS[MXFILTINDX];
  float LAMMIN_OBS[MXFILTINDX];
  float LAMMAX_OBS[MXFILTINDX];

  int   GENRANGE_DMPEVENT[2];  //dump events in this evt # range
  float GENRANGE_DMPTREST[2];  //dump rest-mags for this Trest-range
  int   USEFLAG_DMPTREST ;     // logical flag for above

  int  FORMAT_MASK ;         ;  // 1=verbose, 2=text, 32=FITS ...
  int  WRITE_MASK ;          ;  // computed from FORMAT_MASK
  int  WRFLAG_MODELPAR;    // write model pars to data files (e.g,SIMSED,LCLIB)
  int  WRFLAG_YAML_FILE ;  // write YAML file (Aug 12 2020)

  int   SMEARFLAG_FLUX ;        // 0,1 => off,on for photo-stat smearing
  int   SMEARFLAG_ZEROPT ;      // 0,1 => off,on for zeropt smearing
  int   SMEARFLAG_HOSTGAL;      // host-gal logical flag

  int   KCORFLAG_COLOR;          //   1=> closest color, 2=>Jha color
  float EXPOSURE_TIME;           // global exposure time (1=default)
  float EXPOSURE_TIME_FILTER[MXFILTINDX]; // exposure per filter
  int   EXPOSURE_TIME_MSKOPT ;   // bits 1,2,3 => scale ZPT, SKYSIG,READNOISE

  char KCOR_FILE[MXPATHLEN];        // name of kcor Lookup file

  // define fudges on seeing conditions
  float FORCEVAL_PSF ;         // force PSF value if > 0
  float FUDGESCALE_PSF ;       // scale PSF
  float FUDGESCALE_NOISE_SKY ;  // scale SKY noise
  float FUDGESCALE_NOISE_READ ; // scale CCD/readout noise
  float FUDGESCALE_NOISE_TEMPLATE; // scale template noise
  float FUDGESHIFT_ZPT ;    ;  // shift zero point
  float FUDGESHIFT_ZPT_FILTER[MXFILTINDX]; // ZP shift per filter
  float FUDGESHIFT_LAM;
  float FUDGESHIFT_LAM_FILTER[MXFILTINDX];    // lam shift per filter
  float FUDGESCALE_FLUXERR  ;  // global fudge on true error
  float FUDGESCALE_FLUXERR2 ;  // global fudge on measured error only
  float FUDGESCALE_FLUXERR_FILTER[MXFILTINDX];  // true & measured errors
  float FUDGESCALE_FLUXERR2_FILTER[MXFILTINDX]; // measured error only

  float FUDGE_MAGERR ;      ;  // global mag error fudge added to all errors
  float FUDGE_MAGERR_FILTER[MXFILTINDX]; // idem for each band
  float FUDGE_ZPTERR ;         // global fudge for ZPTERR in SIMLIB file
  float FUDGE_ZPTERR_FILTER[MXFILTINDX] ;

  int   FUDGEOPT_FLUXERR;      // option passed to model of flux-err fudge

  int GENPERFECT;   // 1 => perfect lightcurves with x1000 photostats
  int APPLY_SEARCHEFF_OPT ;      // bit 0,1,2 => trigger, spec, zhost

  // define Gaussian sigmas to assign random syst error(s) during init stage
  INPUTS_RANSYSTPAR_DEF RANSYSTPAR ;

  // define snana-style cut windows
  int APPLY_CUTWIN_OPT ;      // 0= ignore cuts; 1=> apply cuts;
                              // 3=> apply cuts to data, not to SIMGEN-DUMP

  int   NCUTWIN_TOT;
  float EPCUTWIN_LAMREST[2];    // lambda-requirement on all epochs
  float EPCUTWIN_SNRMIN[2];     // min SNR needed for each epoch
  float CUTWIN_TRESTMIN[2];
  float CUTWIN_TRESTMAX[2];
  float CUTWIN_TGAPMAX[2];     // rest-frame maxgap
  float CUTWIN_T0GAPMAX[2];    // rest-frame maxgap near peak
  float CUTWIN_REDSHIFT_TRUE[2];    // cut on true ZCMB
  float CUTWIN_REDSHIFT_FINAL[2];    // cut on "best" measured redshift
  float CUTWIN_HOST_ZPHOT[2];       // cut on host photoz (Aug 19 2018)

  int   NCUTWIN_SNRMAX ;                      // number of SNRMAX cuts
  int   OVERRIDE_CUTWIN_SNRMAX ;              // flag to override SNRMAX cuts
  float CUTWIN_SNRMAX[MXCUTWIN_SNRMAX][2];    // SNRMAX window
  // list of obs filters to require SNRMAX
  char  CUTWIN_SNRMAX_FILTERS[MXCUTWIN_SNRMAX][MXFILTINDX];
  int   CUTWIN_SNRMAX_NFILT[MXCUTWIN_SNRMAX] ;      // Nfilt passing SNRMAX
  char  CUTWIN_SNRMAX_LOGIC[MXCUTWIN_SNRMAX][8];    // specify "N", AND or OR
  float CUTWIN_SNRMAX_TREST[MXCUTWIN_SNRMAX][2];    // rest-frame epoch window

  float CUTWIN_NEPOCH[2];  // [0]=>NOBS ; [1] -> SNRMIN
  int   CUTWIN_NOBSDIF[2]; // NOBS window after MJDDIF cut
  float CUTWIN_MJDDIF[2];  // count NOBS only if passing MJDDIF test
  float CUTWIN_MWEBV[2];
  float CUTWIN_PEAKMAG[2] ;         // OR among filters (Sep 2012)
  float CUTWIN_PEAKMAG_ALL[2] ;     // AND among filters (Aug 2018)

  int   NCUTWIN_PEAKMAG_BYFIELD;
  float CUTWIN_PEAKMAG_BYFIELD[MXCUTWIN_PEAKMJD_BYFIELD][2] ; // idem,by field
  char  CUTWIN_BYFIELDLIST[MXCUTWIN_PEAKMJD_BYFIELD][40]; // for above

  // define CUTWIN_EPOCHS_SNRMIN: 5 20  iz  # SNR>5 for < 20 days, i or z
  float CUTWIN_EPOCHS_SNRMIN;         // Trange applied to SNR>SNRMIN
  float CUTWIN_EPOCHS_TRANGE[2];      // helps pick fast transients
  char  CUTWIN_EPOCHS_FILTERS[MXFILTINDX]; // filters to use
  int   CUTWIN_EPOCHS_NFILT;  // number of filters for TIME_ABOVE cut
  int   CUTWIN_EPOCHS_IFILTLIST[MXFILTINDX];

  // May 2018: define cuts on number of saturated epochs, by filter
  int   NCUTWIN_SATURATE;
  int   CUTWIN_SATURATE_NOBS[MXCUTWIN_SNRMAX][2];
  char  CUTWIN_SATURATE_FILTERS[MXCUTWIN_SNRMAX][MXFILTINDX];

  int   NPE_PIXEL_SATURATE ; // override SIMLIB header value
  int   PHOTFLAG_SATURATE ;
  int   PHOTFLAG_SNRMAX ;
  int   PHOTFLAG_NEARPEAK ;

  // May 2018: define cuts on numbef of UN-saturated epopchs, by filter
  int   NCUTWIN_NOSATURATE;
  int   CUTWIN_NOSATURATE_NOBS[MXCUTWIN_SNRMAX][2];
  char  CUTWIN_NOSATURATE_FILTERS[MXCUTWIN_SNRMAX][MXFILTINDX];

  float EFFERR_STOPGEN ;  // abort when ERREFF is this small per 1000

  int   SNTYPE_Ia_SPEC ;     // SNTYPE for spectropsic subset
  int   SNTYPE_Ia_PHOT ;     // SNTYPE for photo-id subset (-> no spectrum)

  int   GENTYPE_SPEC;   // idem for any model; overrides any other type
  int   GENTYPE_PHOT;   // idem for photo-id subset

  // - - - - - -
  INPUTS_NON1ASED_DEF NON1ASED ;
  // - - - - - -

  int   NON1A_MODELFLAG ;     // specifies NON1ASED or NON1AGRID

  // separate method using MAG-GRID
  char NON1AGRID_FILE[MXPATHLEN];  // FITS file with mag vs. Epoch & z.

  int  CLEARPROMPT; // 1 => prompt before removing old version
  int  REQUIRE_DOCANA ;  // 1 => require DOCUMENTATION keys in maps

  int  NVAR_SIMGEN_DUMP;  // number of SIMGEN variables to write to fitres file
  char VARNAME_SIMGEN_DUMP[MXSIMGEN_DUMP][40] ; // var-names
  bool IS_SIMSED_SIMGEN_DUMP[MXSIMGEN_DUMP];
  int  IFLAG_SIMGEN_DUMPALL ;  // 1 -> dump every generated SN
  int  PRESCALE_SIMGEN_DUMP ;  // prescale on writing to SIMGEN_DUMP file

  // inputs for intrinsic scatter matrix (July 27, 2011)
  int    NCOVMAT_SCATTER ;           // number of non-zero elements
  double COVMAT_SCATTER[3][3] ;
  double COVMAT_CHOLESKY_SCATTER[3][3];
  double COVMAT_SCATTER_SQRT[3][3] ;// Do later
  double COVMAT_SCATTER_REDUCED[3][3] ;

  // silly option to interpolate model on epoch grid
  // (to mimic fake overlays on DES images)
  float TGRIDSTEP_MODEL_INTERP;

  char NONLINEARITY_FILE[MXPATHLEN];

  // stuff for LCLIB model
  char LCLIB_FILE[MXPATHLEN];
  char LCLIB_TEMPLATE_EPOCHS[100];

} INPUTS ;


// define GENLC structure
struct GENLC {

  char SURVEY_NAME[40];    // name of survey in SIMLIB
  char SUBSURVEY_NAME[40]; // subsurvey; e.g,, CFA3 is subsurvey for LOWZ
  int  SUBSURVEY_ID;       // IDSURVEY for SUBSURVEY
  int  IDSURVEY ;
  char primary[40];              // name of primary (AB, VEGA, BD17 ...)

  int  SIMLIB_USEFILT_ENTRY[MXFILTINDX];   // 1=> filter used for this entry

  int  SDSS_SIM ;        // 1= SDSS; 0= non-SDSS (logical)
  int  SIMLIB_ID;        // LIB ID from simlib
  int  SIMLIB_IDMAX;     // max ID in libraray
  //  int  SIMLIB_NWRAP;     // number of library wrap-arounds
  int  SIMLIB_IDLOCK ;   // LIBID that is locked for all gen_event calls
  int  NGEN_SIMLIB_ID ;  // Nuse on same SIMLIB ID until event is accepted

  double SIMLIB_FLUXERR_ADDPAR[MXFILTINDX] ; // fudge-error to add in quad.

  RATEPAR_DEF *RATEPAR ; // selects RATEPAR or RATEPAR_PEC1A

  double RA, DEC ;          // generated position
  double random_shift_RA, random_shift_DEC;     // random coord shift
  double random_shift_RADIUS, random_shift_PHI;
  double GLON, GLAT;        // for LCLIB-galactic models
  double REDSHIFT_HELIO ;   // true Helio redshift of SN
  double REDSHIFT_CMB   ;   // true CMB   redshift of SN
  double REDSHIFT_HOST  ;   // true Helio redshift of host
  int    REDSHIFT_FLAG  ;   // indicates source of redshift

  double DLMU;               // true distMod = 5.0 * log10(DL/10pc),
  double LENSDMU;            // weak lensing DMU (Apr 2017)
  double SL_MAGSHIFT;        // magshift from strong lens magnification

  double PEAKMJD;
  double MJD_RANGE[2];       // MJD range to accept in SIMLIB
  double MJD_RANGE_LC[2];    // MJD range of lightcurve withing GENRANGE_TREST
  int    ISOURCE_PEAKMJD;    // either RANDOM or read from SIMLIB
  double MJD_EXPLODE ;       // explosion time for NON1A or SIMSED
  double DTSEASON_PEAK;      // |PEAKMJD-MJD_seasonEdge|

  double REDSHIFT_RAN[MXZRAN];    // store randoms for redshift smear
  double REDSHIFT_CMB_SMEAR ;     // smeared [measured] CMB redshift
  double REDSHIFT_HELIO_SMEAR ;   // smeared [measured] helio redshift
  double REDSHIFT_SMEAR_ERR ;     // reported error on above
  double VPEC ;                   // true radial peculiar velocity, km/sec
  double VPEC_SMEAR ;             // measured VPEC (Jan 2018)
  double REDSHIFT_MAX_SNR5 ;      // zmax with SNR > 5 (for monitor only)

  int    CORRECT_HOSTMATCH ;  // 1=correct match, 0=wrong host

  double PEAKMJD_SMEAR  ; // smeared PEAKMJD to simulate survey fit error
  double PEAKMJD_RANGauss ;

  double RESTLAM_MODEL[2]; // rest-frame wavelength range of model

  double  NOSHAPE ;     // undefined, such as nonIa or FIXMAG
  double  FIXMAG ;      // for FIXMAG model only
  double  STRETCH;
  double  DELTA;        // for MLCS2k2 model
  double  DM15;         // for DM15 model
  double  SHAPEPAR ;
  double *ptr_SHAPEPAR ;

  double SALT2x0 ;
  double SALT2x1 ;
  double SALT2c ;
  double SALT2mB;     // peak B-band mag
  double SALT2alpha ;
  double SALT2beta ;
  double SALT2gammaDM ;       // mag shift from gamma/host correlation

  double GENMAG_OFF_GLOBAL ;  // INPUTS.GENMAG_OFF_GLOBAL + z-dependence

  char  SNTYPE_NAME[80];   // 1a, 1b, 1c, II, etc ...
  char  SNTEMPLATE[80];

  char  DISTANCE_NAME[40];    // DLMAG, mB, ...
  char  COLORPAR_NAME[40];    // c, AV ...
  char  SHAPEPAR_NAME[40] ;   // name of SHAPEPAR: DM15, DELTA, etc ...
  char  SHAPEPAR_GENNAME[40]; // GENPEAK_[name], GENSIGMA_[name], etc ...
  char  COLORPAR2_NAME[40] ;  // RV, beta ...
  char  SHAPEPAR2_NAME[40] ;  // alpha ...

  int    SIMSED_IPARMAP[MXPAR_SIMSED];     // IPAR list in COV mat
  double SIMSED_PARVAL[MXPAR_SIMSED];    // params for SIMSED model

  double AVTAU;       // exponential tau for host extinction
  double AVSIG;       // Gauss sigma
  double AV0RATIO ;   // Guass/expon ratio at AV=0


  GEN_EXP_HALFGAUSS_DEF GENPROFILE_AV; // 3/2020: added by djb for z-dep
  GEN_EXP_HALFGAUSS_DEF GENPROFILE_EBV_HOST; // 3/2020: added by djb for z-dep

  double AV ;         // host extinction param for CCM89 law
  double RV ;         // host extinction param for CCM89 law

  double MWEBV ;             // E(B-V) from Shlagel maps
  double MWEBV_SMEAR;        // smeared E(B-V) applied to SN mag.
  double MWEBV_ERR ;         // assigned error (for fitting unc.)
  double FLUXCOR_MWEBV_MAP[MXFILTINDX]; // analysis flux-cor per filter
  double MAGCOR_MWEBV_MAP[MXFILTINDX];  // analysis mag-cor per filter
  double MAGCOR_MWEBV_TRUE[MXFILTINDX];
  double RA_LAST, DEC_LAST;  // re-compute MWEBV only if RA,DEC change

  int NFILTDEF_OBS;     // number of user-defined observer filters
  int IFILTMAP_OBS[MXFILTINDX];  // obs-filter index vs. sparse index
  int IFILTINVMAP_OBS[MXFILTINDX]; // sparse index vs. obs-filter index

  int IFILTMAP_REST1[MXFILTINDX]; // nearest rest-filter index vs. sparse index
  int IFILTMAP_REST2[MXFILTINDX]; // 2nd near rest-filt index vs. sparse index
  int IFILTMAP_REST3[MXFILTINDX]; // 3rd nearest

  int  NFILTDEF_SPECTROGRAPH ;
  int  IFILTDEF_SPECTROGRAPH[MXFILTINDX];     // vs. sparse synth spectro bands
  int  IFILTINV_SPECTROGRAPH[MXFILTINDX];     // vs. absolute ifilt_obs
  char FILTERLIST_SPECTROGRAPH[MXFILTINDX] ;  // idem
  int  IFLAG_SYNFILT_SPECTROGRAPH[MXFILTINDX] ;  // vs. ifilt_obs

  double LAMDIF_REST1[MXFILTINDX];
  double LAMDIF_REST2[MXFILTINDX];
  double LAMDIF_REST3[MXFILTINDX];

  int NFILTDEF_SIMLIB;
  int IFILTMAP_SIMLIB[MXFILTINDX];  // simlib-filter index vs. sparse index
  int IFILTINVMAP_SIMLIB[MXFILTINDX];  // inverse map

  int  NFILTDEF_REST;    // for model
  int  IFILTMAP_REST[MXFILTINDX] ;
  int  IFILTINVMAP_REST[MXFILTINDX];
  char FILTLIST_REST[MXFILTINDX] ;

  int DOFILT[MXFILTINDX] ; // logical to generate ifilt_obs

  // define max allowed filter indices. For example, if we simulate
  // SDSS gri only (ifilt=2,3,4), then MINDEF,MAXDEF=12345 for ugriz.
  // MINDEF_REST is Bessell-U for MLCS2k2, even if U-band is not used.
  // these variables are used as offsets.


  int   CID ;           // internal data CID or 40000 + ilc
  int   CIDOFF ;       // CID offset depends on MJD range (random only)
  int   CIDRAN ;       // use this random CID (if INPUTS.CIDRAN > 0)
  //  int   YEAR ;         // survey year simulated
  int   SUBSAMPLE_INDEX ; // only if NSUBSAMPLE_MARK > 0 (June 2017)
  int   NEPOCH;        // includes model-epoch at T=0 and epoch with fluxerr<0

  int   ACCEPTFLAG ;       // True of event is accepted
  int   ACCEPTFLAG_FORCE ; // true if same LIBID is read too many times
  int   ACCEPTFLAG_LAST ;

  double  MJD[MXEPSIM];
  int     INDEX_MJDSORT[MXEPSIM];
  int     IFILT_OBS[MXEPSIM];    // observer filter vs. epoch
  int     IEPOCH_NEARPEAK ;    // epoch-index closest to peak (for trigger)
  double  DTPEAK_MIN ;        // t-Tpeak for nearest obs.

  // variables used for cuts
  int   NOBS_MJDDIF ;   // NOBS passing MJD dif cut
  int   NOBS_SNR ;      // NOBS passing SNR cut
  int   NOBS ;          // real obs only
  int   NOBS_FILTER[MXFILTINDX]; // NOBS per band
  int   NOBS_UNDEFINED;   // suppresed NOBS because model is undefined
  int   NOBS_SATURATE[2]; // 0=unsatured, 1=saturated
  int   NOBS_SATURATE_FILTER[2][MXFILTINDX]; // idem, filter-dependent

  float TRESTMIN ;
  float TRESTMAX ;
  float TGAPMAX  ;
  float T0GAPMAX ;
  float TIME_ABOVE_SNRMIN ;

  double SNRMAX_GLOBAL ;
  double SNRMAX_FILT[MXFILTINDX];   // vs. filter
  double SNRMAX_SORTED[MXFILTINDX]; // [1]=>max among all filt, [2]=>2nd max,
  double SNRMAX_TRUE[MXFILTINDX];   // trueFlux/error

  double FLUXADU_at_SNRMAX[MXFILTINDX] ;   // needed for FUDGE2_SNRMAX option
  double FLUXpe_at_SNRMAX[MXFILTINDX] ;   // needed for FUDGE2_SNRMAX option

  int NEWMJD ;
  int EPOCH_RANGE_NEWMJD[MXEPSIM][2];

  char  TELESCOPE[MXEPSIM][40];   // name of telescope
  char  FIELDNAME[MXEPSIM][MXCHAR_FIELDNAME] ;  // name of field for each obs
  int   IFLAG_GENSOURCE ;     // specifies GENSOURCE

  // define scale factors that depend on user-input EXPOSURE_TIME
  double SCALE_NOISE[MXFILTINDX] ;       // scale noise in pe = sqrt(above)
  double SHIFT_ZPTSIMLIB[MXFILTINDX] ;   // zpt shift for EXPOSURE time

  // - - - - - - - - - - - -

  // flux and magnitudes (float-> double, Jan 2014)
  double flux[MXEPSIM] ;            // flux in ADU
  double fluxerr_data[MXEPSIM];     // reported error in data file
  double template_err[MXEPSIM];     // correlated template error
  double trueSNR[MXEPSIM];          // true/generated SNR
  int    npe_above_sat[MXEPSIM];    // nphotoelectrons above saturation

  FLUXNOISE_DEF *FLUXNOISE;    // Dec 27 2019 - refactor for noise cov.

  MONITOR_REDCOV_FLUXNOISE_DEF MONITOR_REDCOV_FLUXNOISE[MXFILTINDX][NTYPE_FLUXNOISE];

  // xxxx -----------------------------------------------------
  // xxxxx legacy arrays to remove after GENFLUX_DRIVER refactor
  // noise contributions (in photoelectrons)
  double NOISE_SN[MXEPSIM] ;
  double NOISE_SKY[MXEPSIM] ;
  double NOISE_TEMPLATE[MXEPSIM] ; //correlated noise from TEMPLATE_SKY[CCD]SIG
  double NOISE_CCD[MXEPSIM];
  double NOISE_HOSTGAL_PHOT[MXEPSIM] ;  // galaxy shot noise
  double NOISE_HOSTGAL_IMAGE[MXEPSIM] ; // anomalous noise from HOSTNOISE_FILE
  double NOISE_TOTAL[MXEPSIM] ;

  double NOISE_AREA[MXFILTINDX]; // effective area for sky noise
  double NOISE_PSF[MXFILTINDX];  // effective PSF = sqrt[ AREA/(4*PI) ]

  double SNR_CALC[MXEPSIM] ;    // used for trigger effic (Aug 24 2014)
  double SNR_MON[MXEPSIM];      // calculated SNR for MAGMONITOR_SNR input
  //xxx ----------------------- end legacy ------------

  // Gaussian randoms for broadband measurement noise
  double RANGauss_NOISE_SEARCH[MXEPSIM];   // search noise, per epoch
  double RANGauss_NOISE_FUDGE[MXEPSIM];   //  fudged search noise, per epoch
  double RANGauss_NOISE_ZP[MXEPSIM];       // ZP noise, per epoch (Feb 2018)
  double RANGauss_NOISE_TEMPLATE[MXFIELD_OVP][MXFILTINDX];  // template noise per filter and field-overlap

  // GENSMEAR refers to intrinsic scatter models
  double  MAGSMEAR_COH[2];              // coherent part of scatter
  double  GENSMEAR_RANGauss_FILTER[MXFILTINDX+1]  ;  // filter smear

  double  SPECEFF_RAN[MXFILTINDX+1]  ;
  double  magsmear8[MXEPSIM];        // actual intrinsic mag-smear

  double  NON1AGRID_RANFlat;    // to pick random template id
  double  NON1AGRID_RANGauss ;  // to pick random magSmear

  double  COVMAT_SCATTER_GRAN[3];   // Gaussian randoms for intrinsic scatter
  double  COVMAT_SCATTER[3];        // scattered values
  char    COVMAT_SCATTER_NAME[3][40]; // name of each scatter term

  // - - - - -
  double  epoch_obs_range[2];   // min and max epoch, all bands
  double  epoch_obs[MXEPSIM];       // observer epoch = MJD - PEAKMJD
  double  epoch_rest[MXEPSIM];      // rest epoch relative to peak, days
  double  genmag_obs[MXEPSIM] ;     // generated obs  magnitude
  double  generr_obs[MXEPSIM] ;     // obs mag err from model
  double  peakmag_obs[MXFILTINDX] ;
  double  genmag_obs_template[MXFILTINDX]; // for LCLIB model (stars)

  double  genmag_rest[MXEPSIM] ;    // idem, rest frame
  double  generr_rest[MXEPSIM] ;    // idem, rest frame
  double  peakmag_rest[MXFILTINDX] ;

  double  genmag_rest2[MXEPSIM] ;    // 2nd nearest rest-frame mag.
  double  generr_rest2[MXEPSIM] ;    // 2nd nearest rest-frame mag.
  double  peakmag_rest2[MXFILTINDX] ;

  double  genmag_rest3[MXEPSIM] ;    // 2nd nearest rest-frame mag.
  double  generr_rest3[MXEPSIM] ;    // 2nd nearest rest-frame mag.
  double  peakmag_rest3[MXFILTINDX] ;

  int     NWIDTH_SIMGEN_DUMP;
  double  WIDTH[MXFILTINDX];  // generated LC width per band (for monitor)

  double AVwarp8[MXEPSIM];
  int    ifilt_AVwarp[MXEPSIM][2];

  double kcorval8[MXEPSIM];    // store kcor for archive
  char   kcornam[MXEPSIM][8] ; // kcor name
  double warpcolval8[MXEPSIM];    // store rest color used for warping
  char   warpcolnam[MXEPSIM][8] ;  // color used for warping

  double kcor_RVMW ;
  int    kcor_OPT_MWCOLORLAW ;

  float mag[MXEPSIM] ;       // observed mag
  float mag_err[MXEPSIM] ;   // error onabove

  // Jan 19 2020:
  //   ISPEAK    -> OBSFLAG_PEAK
  //   ISOBS     -> OBSFLAG_GEN
  //   USE_EPOCH -> OBSFLAG_WRITE
  //   ISTEMPLATE -> OBSFLAG_TEMPLATE
  bool OBSFLAG_GEN[MXEPSIM];     // flag to generate mags & flux
  bool OBSFLAG_WRITE[MXEPSIM];   // write these to data file
  bool OBSFLAG_PEAK[MXEPSIM] ;   // extra epochs with peak mags
  bool OBSFLAG_TEMPLATE[MXEPSIM]; // extra epochs that are templates

  int IEPOCH_PEAK[MXFILTINDX] ; // artificial peak epoch vs.ifilt_obs
  int IEPOCH_SNRMAX_GLOBAL;       // global epoch with SNRMAX (Jun 2018)
  int IEPOCH_SNRMAX[MXFILTINDX] ; // idem vs. ifilt_obs
  int IEPOCH_TEMPLATE[MXFILTINDX]; // identifies template epochs

  float COVAR[MXFILT_COVAR][MXEPCOV][MXFILT_COVAR][MXEPCOV] ;  // cov matrix
  float PEAKMAGERR_MODEL[MXFILTINDX];

  int     SEARCHEFF_MASK;     // search eff mask
  double  SEARCHEFF_SPEC;     // EFF(spec-confirmed)
  double  SEARCHEFF_zHOST;    // zHOST efficiency (only of not spec-confirmed)
  double  MJD_TRIGGER;           // min MJD when search trigger is satisfied.
  double  MJD_DETECT_FIRST;   // mjd of first detection
  double  MJD_DETECT_LAST;    // mjd of last detection

  int     CUTBIT_MASK;           // snana cutbit mask
  float   GENEFF, GENEFFERR;     // generation efficiency & error

  double RISETIME_SHIFT;  // shift in rise-time, days
  double FALLTIME_SHIFT;  // shift in 15-day fall-time, days

  int   SIMTYPE;              // always set; SNTYPE not always set
  int   SNTYPE;               // user index-tag (see 'SNTYPE' non1a key)
  int   TEMPLATE_INDEX ;    // template index for NONA1SED, SIMSED, LCLIB

  GENLC_NON1ASED_DEF NON1ASED ;  // Mar 24 2016

  // spectro and phot tags
  int METHOD_TYPE ;    // SPEC or PHOT
  int NTYPE_SPEC;      // number of spectroscopic tags
  int NTYPE_SPEC_CUTS; // # spec-tags after cuts
  int NTYPE_PHOT;      // number of photometric events
  int NTYPE_PHOT_CUTS; // idem after cuts
  int NTYPE_PHOT_WRONGHOST;  // idem with wrong host (4.2019)
  float FRAC_PHOT_WRONGHOST; // true fraction with wrong host

  // misc.
  int   STOPGEN_FLAG;
  int   FUDGE_SNRMAX_FLAG ;  // 0 or 1 or 2

} GENLC ;


// strong lens structure (July 2019)
struct GENSL {
  int INIT_FLAG ;
  int REPEAT_FLAG;     // T => repeat image
  int NIMG;            // number of images to process
  int IMGNUM;          // image-num being processed
  int IDLENS;
  int BLEND_FLAG;
  double zSN, zLENS;
  double PEAKMJD_noSL;    // undelayed PEAKMJD
  double RA_noSL, DEC_noSL;
  double MJDMIN, MJDMAX;  // used for SIMLIB read
  double *TDELAY_LIST, *XIMG_LIST, *YIMG_LIST;
  double *MAGNIF_LIST, *MAGSHIFT_LIST ;
} GENSL ;


// temp structure used by NEPFILT_GENLC
// Jan 28 2022: define pointers to be malloced (save memory)
struct GENFILT {
  
  double *Trest[MXFILTINDX]; 
  double *Tobs[MXFILTINDX]; 
  double *genmag_obs[MXFILTINDX];  
  double *genmag_smear[MXFILTINDX];

  double *genmag_rest[MXFILTINDX]; 
  double *genmag_rest2[MXFILTINDX];
  double *genmag_rest3[MXFILTINDX];

  double *generr_obs[MXFILTINDX]; 
  double *generr_rest[MXFILTINDX];
  double *generr_rest2[MXFILTINDX];
  double *generr_rest3[MXFILTINDX];
} GENFILT ;


int NGENLC_TOT ;             // actual number of generated LC
int NGENLC_WRITE ;           // number written
int NGENLC_TOT_SUBSURVEY[MXIDSURVEY];
int NGENLC_WRITE_SUBSURVEY[MXIDSURVEY];

int NGENSPEC_TOT;            // total number of generated spectra
int NGENSPEC_WRITE ;         // number of spectra written
int NGENFLUX_DRIVER;         // number of calls to GENFLUX_DRIVER

struct NGEN_REJECT {
  int GENRANGE, GENMAG;
  int GENPAR_SELECT_FILE ;
  int HOSTLIB ;    // subset of GENRANGE with no valid HOST
  int SEARCHEFF;
  int CUTWIN ;
  int NEPOCH ;   // counts NEPOCH < NEPOCH_MIN
} NGEN_REJECT ;


// valid Z-range with defined rest-frame model for each obs-filter
// (for README comment only)
double ZVALID_FILTER[2][MXFILTINDX] ;
int    NSKIP_FILTER[MXFILTINDX]; // number of times obs-filter is skipped
int    NGEN_ALLSKIP ;

  // define SIMLIB_DUMP struct vs. SIMLIB  entry
int NREAD_SIMLIB ;
typedef struct {
  double MJDMIN;
  double MJDMAX;

  double NEPFILT[MXFILTINDX]; // 0=> all epochs
  double GAPMAX[MXFILTINDX];  // max gap in days, [0] => all filters
  double GAPAVG[MXFILTINDX];  // avg gap in days, [0] => all filters
  double GAIN[MXFILTINDX];      // CCD GAIN
  double ZPT[MXFILTINDX];       // avg ZPT
  double PSF[MXFILTINDX];       // avg PSF
  double SKYSIG_ADU[MXFILTINDX];    // avg sky noise/pix
  double SKYSIG_pe[MXFILTINDX];    // avg sky noise/pix
  double SKYMAG[MXFILTINDX];    // sky brightness (mag/arcsec^2)
  double M5SIG[MXFILTINDX];     // 5-sigma limiting mag (calculated quantity)
  double FOM[MXFILTINDX];       // figure of merit
} SIMLIB_DUMP_DEF ;

SIMLIB_DUMP_DEF  SIMLIB_DUMP_AVG1   ;  // avg for one SIMLIB entry
SIMLIB_DUMP_DEF  SIMLIB_DUMP_AVGALL ;  // average over all SIMLIB entries
SIMLIB_DUMP_DEF  SIMLIB_DUMP_NAVGALL;  // how many used for avg in AVGALL

char SIMLIB_DUMPFILE_SEQ[MXPATHLEN]; // TEXT: one row per SEQUENCE
char SIMLIB_DUMPFILE_OBS[MXPATHLEN]; // TEXT: one row per OBS
char SIMLIB_DUMPFILE_ROOT[MXPATHLEN];  // covert to root

#define SIMLIB_DUMPMASK_SEQ   1  // one-row summary for each SEQUENCE (LIBID)
#define SIMLIB_DUMPMASK_OBS   2  // one-row summary for each OBS

#define TABLEID_SIMLIB_DUMP    7788
#define TABLENAME_SIMLIB_DUMP  "SIMLIB"


struct SIMLIB_GLOBAL_HEADER {

  // stuff read from global SIMLIB header
  char SURVEY_NAME[60];
  char SUBSURVEY_LIST[MXPATHLEN];
  char FILTERS[MXFILTINDX];  // global list of all filters
  char FIELD[60];            // Nov 2021
  char TELESCOPE[60];
  char PSF_UNIT[40] ;
  bool NEA_PSF_UNIT;
  char SKYSIG_UNIT[40];
  char USERNAME[40];
  int  NLIBID, NLIBID_VALID ;
  double PIXSIZE, SOLID_ANGLE ;
  int  NPE_PIXEL_SATURATE;    // Jan 3, 2018
  int  PHOTFLAG_SATURATE ;
  int  PHOTFLAG_SNRMAX ;
  int  PHOTFLAG_NEARPEAK ;

  int    NGENSKIP_PEAKMJD ; // number of PEAKMJD ranges to skip
  double GENSKIP_PEAKMJD[MXGENSKIP_PEAKMJD_SIMLIB][2] ;

  // flux-error correction maps
  char   FLUXERR_ADD_FILTERS[MXFILTINDX];
  double FLUXERR_ADD_VALUES[MXFILTINDX]; // values to add in quadrature

  // LEGACY: error scales vs. filter and LOG10(SNR)
  int    NFLUXERR_COR ;
  char   FLUXERR_COR_FILTERS[MXFLUXERR_COR_SIMLIB][MXFILTINDX];
  float  FLUXERR_COR_LOG10SNR[MXFLUXERR_COR_SIMLIB];
  float  FLUXERR_COR_SCALE[MXFLUXERR_COR_SIMLIB][MXFILTINDX];

} SIMLIB_GLOBAL_HEADER ;


struct SIMLIB_HEADER {
  // header info for each LIBID entry

  // required
  char   LIBNAME[20];        // LIB[LIBID] (for SNTABLE)
  char   SUBSURVEY_NAME[40]; // optional sub-survey (e..g, CFA3)
  int    NOBS, LIBID, NWRAP ;
  int    NOBS_APPEND ;  // these obs are not MJD-sorted (Jan 2018)
  int    NOBS_SIM_MAGOBS; // NOBS with SIM_MAGOBS<99
  double RA, DEC ;

  // optional stuff
  double MWEBV, PIXSIZE ;
  int    FAKEID, CCDNUM ;
  long long GALID;

  // these header keys can be changed anywhere in the simlib entry
  char TELESCOPE[60] ; // July 2016
  char FIELD[60], FIELDLIST_OVP[MXFIELD_OVP][MXCHAR_FIELDNAME];
  int  NFIELD_OVP ;

  // optional GENRANGES to re-generate
  int    REGEN_FLAG ;
  GENGAUSS_ASYM_DEF GENGAUSS_PEAKMJD ;
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2x1 ;
  GENGAUSS_ASYM_DEF GENGAUSS_SALT2c ;
  double GENRANGE_REDSHIFT[2] ;
  double CUTWIN_REDSHIFT[2] ;

  // season stuff computed from sorted MJDs (Apr 2017)
  int    NSEASON ;
  int    USE_SEASON[MXSEASON_SIMLIB]; // logical flag for each used season
  double MJDRANGE_SEASON[MXSEASON_SIMLIB][2];
  double TLEN_SEASON[MXSEASON_SIMLIB] ;
  double MJDRANGE_SURVEY[2];

  int NREPEAT; // repeat counter

  // counters to help flag reason for aborting
  int NFOUND_TOT  ; // total number of SIMLIBs read
  int NFOUND_MJD  ; // idem with valid MJD
  int NFOUND_RA   ; // idem with valid RA
  int NFOUND_DEC  ; // idem with valid DEC
  int NFOUND_FIELD ;
  int NFOUND_GENCUTS; // May 30 2020

} SIMLIB_HEADER ;



// Define SIMLIB  struct for reading
// Allocate for writing/reading entire array.
typedef struct  {
  int     NOBS;       // everything, including SPECTROGRAPH and TAKE_SPECTRUM
  int     NOBS_READ ; // orginal NOBS read from cadence (never changes)
  int     NOBS_SPECTROGRAPH ;
  int     NOBS_TAKE_SPECTUM ;

  int     OPTLINE[MXOBS_SIMLIB];
  int     IFILT_OBS[MXOBS_SIMLIB];    // absolute filter index

  char    *PTR_BAND[MXOBS_SIMLIB];
  char    BAND[MXOBS_SIMLIB][4];

  int     IDEXPT[MXOBS_SIMLIB];
  int     NEXPOSE[MXOBS_SIMLIB];  // Jan 2018 (for saturation calc)
  double  MJD[MXOBS_SIMLIB];
  double  CCDGAIN[MXOBS_SIMLIB];
  double  READNOISE[MXOBS_SIMLIB];
  double  SKYSIG[MXOBS_SIMLIB];
  double  PSFSIG1[MXOBS_SIMLIB];
  double  PSFSIG2[MXOBS_SIMLIB];
  double  PSFRATIO[MXOBS_SIMLIB];
  double  NEA[MXOBS_SIMLIB];
  double  ZPTADU[MXOBS_SIMLIB];    // ZPT in ADU for entire exposure
  double  ZPTERR[MXOBS_SIMLIB];    // ZPT error
  double  MAG[MXOBS_SIMLIB];       // optional mag
  double  PIXSIZE[MXOBS_SIMLIB] ;   // Nov 26, 2011

  char    *PTR_FIELDNAME[MXOBS_SIMLIB];
  char    FIELDNAME[MXOBS_SIMLIB][MXCHAR_FIELDNAME];
  char    TELESCOPE[MXOBS_SIMLIB][40];
  int     APPEND_PHOTFLAG[MXOBS_SIMLIB];  // Jan 201

  double  TEMPLATE_SKYSIG[MXOBS_SIMLIB] ;
  double  TEMPLATE_READNOISE[MXOBS_SIMLIB] ;
  double  TEMPLATE_ZPT[MXOBS_SIMLIB] ;

  int ISTORE_RAW[MXOBS_SIMLIB] ;
  int ISEASON[MXOBS_SIMLIB] ;  // bracketed by 90+ day gaps

  // spectrograph info from SIMLIB or TAKE_SPECTRUM key
  int    OBSLIST_SPECTROGRAPH[MXOBS_SPECTROGRAPH];   // sparse list of SPECTROGRAPH obs
  int    OBSLIST_TAKE_SPECTRUM[MXOBS_SPECTROGRAPH];  // sparse list of TAKE_SPECTRUM

  int    IFILT_SPECTROGRAPH[MXOBS_SIMLIB] ;  // synthetic filter index
  int    INDX_TAKE_SPECTRUM[MXOBS_SIMLIB] ;  // 0 to NPEREVT_TAKE_SPECTRUM-1
  double TEXPOSE_SPECTROGRAPH[MXOBS_SIMLIB]; // exposure time

} SIMLIB_OBS_DEF ;


SIMLIB_OBS_DEF SIMLIB_OBS_RAW ;  // read from simlib
SIMLIB_OBS_DEF SIMLIB_OBS_GEN ;  // used to generate SN



// Jan 6 2016 - define contiguous temp arrays used to sort SIMLIB by MJD.
struct {
  int     NMJD ;
  int     INDEX_SORT[MXOBS_SIMLIB] ;
  double  MJD[MXOBS_SIMLIB] ;
  double  MJD_LAST ;
  double  MJDOFF ;
} SIMLIB_LIST_forSORT ;


#define ALLFIELDS  "ALLFIELDS"
struct SIMLIB_TEMPLATE {

  int    USEFLAG ;      // logical to use correlated template noise
  int    NFIELD_OVP ;   // number of overlapping fields; resets each SIMLIB
  char   FIELDNAME[MXFIELD_OVP][MXCHAR_FIELDNAME];

  double SKYSIG[MXFIELD_OVP][MXFILTINDX] ;
  double READNOISE[MXFIELD_OVP][MXFILTINDX] ;
  double ZPT[MXFIELD_OVP][MXFILTINDX] ;

  double TEXPOSE_SPECTROGRAPH ;    // for SPECTROGRAPH

} SIMLIB_TEMPLATE ;


// LEGACY FLUXERR_COR map structure (Dec 2011); COR <-> correction
struct SIMLIB_FLUXERR_COR {
  int     USE ;
  double **LOG10SNR ;  // snana-computed SNR (filter, SNR-bin)
  double **SCALE ;     // error correction scale (filter, SNR-bin)
  int    *MAPSIZE ;    // size of map per filter
} SIMLIB_FLUXERR_COR ;


int GENFRAME_OPT;        // one of below, based on model option
#define GENFRAME_REST 1  // => generate in rest frame; then boost
#define GENFRAME_OBS  2  // => generate directly in obs frame
#define GENFRAME_HOST 3  // => used by TAKE_SPECTRUM on host instead of SN
#define GENFRAME_MJD  4  // => used by TAKE_SPECTRUM for MJD option

int INDEX_GENMODEL         ;  // index for model
int SUBINDEX_GENMODEL      ;  // sub-index; e.g, separate SALT2, SALT3
int INDEX_GENMODEL_TWEAK   ;  // index for model tweak
int LGEN_SNIA              ;  // =1  => generating a IA model
bool IS_PySEDMODEL         ;  // python SED model (BYOSED, SNEMO)

#define MODEL_TWEAK_SHOCK 1 // early Shock (Kasen 2006)

#define  IFLAG_GENRANDOM   1
#define  IFLAG_GENGRID     4

#define OPT_SNXT_CCM89  1  // use exact CCM89 model to apply host extinc
#define OPT_SNXT_SJPAR  2  // use Saurabh alpha,beta,zeta paramitrization

#define CUTBIT_TRESTMAX     0   // (1)
#define CUTBIT_TRESTMIN     1   // (2)
#define CUTBIT_SNRMAX       2   // (4)   max SNR cut
#define CUTBIT_TGAPMAX      3   // (8)   max Tgap
#define CUTBIT_T0GAPMAX     4   // (16)  idem near peak
#define CUTBIT_NOBS_SNR     5   // (32)  NOBS with special SNR cut
#define CUTBIT_NOBS_MJDDIF  6   // (64)  NOBS with MJDDIF cut
#define CUTBIT_MWEBV        7   // (128) galactic extinction
#define CUTBIT_REDSHIFT     8   // (256) redshift
#define CUTBIT_PEAKMAG      9   // (512) peak mag
#define CUTBIT_TIME_ABOVE  10   // (1024) time above SNRMIN
#define CUTBIT_SATURATE    11   // (2048) saturation cuts

#define ALLBIT_CUTMASK    4095   // 2^(maxbit+1)-1

// define strings to contain info about simulated volume & time
int  NLINE_RATE_INFO;
char LINE_RATE_INFO[MXEPSIM][80];


// define up to four allowed names for each model (e.g., 'mlcs2k2', 'mlcs')
#define MXNAME_PER_MODEL 4
char GENMODEL_NAME[MXMODEL_INDEX][MXNAME_PER_MODEL][20];

struct GENPERFECT {
  int NVAR;
  char  parnam[20][40];   // name modified INPUTS var
  float parval[20][2];    // orig[0] and modified[1] par value
  int   partype[20];          // 1=int ; 2=float
} GENPERFECT ;


int NEVT_SIMGEN_DUMP ; // NEVT written to SIMGEN_DUMP file
int NVAR_SIMGEN_DUMP ; // total define SIMGEN variables
                 // note that INPUTS.NVAR_SIMGEN_DUMP is how many user var
int NVAR_SIMGEN_DUMP_GENONLY; // variables for generation only
int INDEX_SIMGEN_DUMP[MXSIMGEN_DUMP]; // gives strucdt index vs. [user ivar]

struct SIMGEN_DUMP {
  float      *PTRVAL4 ;
  double     *PTRVAL8 ;
  int        *PTRINT4 ;
  long long  *PTRINT8 ;  // added Feb 2015 to handle galaxy ID
  char       *PTRCHAR ;  // allows FIELD (7.30.2014)

  char  VARNAME[20] ;
} SIMGEN_DUMP[MXSIMGEN_DUMP] ;

struct SIMGEN_DUMMY {
  float      VAL4 ;
  double     VAL8 ;
  int        IVAL4 ;
  long long  IVAL8 ; // Feb 2015
  char   CVAL[60];
} SIMGEN_DUMMY ;



int  NAVWARP_OVERFLOW[MXFILTINDX];
char WARNING_AVWARP_OVERFLOW[100];


// arrays used inside init_genSEDMODEL (Jan 2014)
struct genSEDMODEL {
  double *lam ;
  double *primaryFlux ;
  double *TransSN, *TransREF ;
} genSEDMODEL ;


// store list of genpar cut-windows from GENPAR_SELECT_FILE
struct {
  int    NSTORE ;
  double *ZCMBmin,  *ZCMBmax, *PKMJDmin, *PKMJDmax ;
  double *RAmin, *RAmax, *DECmin, *DECmax ;
  char   **SNID,  **FIELDNAME ;
} GENPAR_SELECT ;

// ================================================
//  Z-dependent variations of sim-parameters

#define MXPAR_ZVAR     150   // max number of ZVAR variables
#define MXZBIN_ZVAR    100   // max number of zins to define Z-variation
#define POLYORDER_ZVAR   3  // order of polynomial with ZPOLY option
#define FLAG_ZPOLY_ZVAR  1
#define FLAG_ZBINS_ZVAR  2

int  USE_ZVAR_FILE ; // logical flag if ZVARIATION_FILE is set
int  NPAR_ZVAR_TOT;  // total number of defined sim-params to choose from
int  NPAR_ZVAR_USR;  // user-number of sim-params specified with z-dependence
char INPUT_ZVARIATION_FILE[MXPATHLEN] ;   // optional file with input info
struct {
  int     FLAG  ;         // POLYnomial or ZBINS
  char    PARNAME[40];    // name of sim parameter to vary with z
  GENPOLY_DEF POLY;

  // optional variation in z-bins to allow any functional form
  int     NZBIN;
  double  ZBIN_VALUE[MXZBIN_ZVAR];       // z-bins
  double  ZBIN_PARSHIFT[MXZBIN_ZVAR];    // shift at each z-bin

} INPUT_ZVARIATION[MXPAR_ZVAR];

// define valid SIM parameters to apply Z-dependence.
char PARDEF_ZVAR[MXPAR_ZVAR+1][40] ;


// ------------------------------------------
//  Declare functions
// ------------------------------------------

void   init_commandLine_simargs(int argc, char **argv);
int    LUPDGEN(int N);

void   SIMLIB_INIT_DRIVER(void);
void   SIMLIB_initGlobalHeader(void);
void   SIMLIB_readGlobalHeader_TEXT(void);
void   SIMLIB_prepGlobalHeader(void);
void   SIMLIB_prep_fluxerrScale_LEGACY(void);
void   SIMLIB_findStart(void);

void   SIMLIB_READ_DRIVER(void);
void   SIMLIB_readNextCadence_TEXT(void);
void   SIMLIB_addCadence_SPECTROGRAPH(void) ;
void   SIMLIB_prepCadence(int REPEAT_CADENCE);
void   SIMLIB_prepMJD_forSORT(int ISTORE);
void   SIMLIB_sortbyMJD(void);
void   SIMLIB_randomize_skyCoords(void);

void   init_SIMLIB_HEADER(void);
int    keep_SIMLIB_HEADER(void);
int    keep_SIMLIB_OBS(int isort);

int    SIMLIB_read_templateNoise(char *field, char *whatNoise, char **wdlist);
void   SIMLIB_TAKE_SPECTRUM(void) ;

void   set_TIMERS(int flag);

int    SKIP_SIMLIB_FIELD(char *field);
int    USE_SAME_SIMLIB_ID(int IFLAG) ;
void   set_SIMLIB_NREPEAT(void);

void   store_SIMLIB_SEASONS(void);
void   set_SIMLIB_MJDrange(int OPT, double *MJDrange);
void   remove_short_SIMLIB_SEASON(void);

void   store_SIMLIB_SPECTROGRAPH(int ifilt, double *VAL_STORE, int ISTORE);
void   store_GENSPEC(double *VAL_STORE);

void   get_SPECTROGRAPH_ZPTPSFSKY(int OBSRAW, int ifilt,
				  double TEXPOSE_S, double TEXPOSE_T,
				  double *ZPT, double *PSFSIG,
				  double *SKYSIG_S, double *SKYSIG_T);

void   ABORT_SIMLIB_FILTER(int isort);
//void   ABORT_SIMLIB_FILTER(int OPTLINE, double MJD, char *cfilt) ;

int    IFIELD_OVP_SIMLIB(int OPT, char *FIELD) ;
void   GENFILTERS_CHECK(void);


double get_SIMLIB_fluxerrScale_LEGACY(int ifiltobs, double SNR ) ;

void   get_SIMLIB_SCALES( int ifilt_obs, double *SHIFT_ZPT,
			  double *SCALE_SKYSIG, double *SCALE_SKYSIG_T,
			  double *SCALE_READNOISE ) ;

double SIMLIB_angsep_min(int NSTORE, double RA, double DEC,
			 double *RA_STORE, double *DEC_STORE);
int    parse_SIMLIB_ZPT(char *cZPT, double *ZPT,
			char *cfiltList, int *ifiltList) ;
// xxx mark void   parse_SIMLIB_GENRANGES(FILE *fp_SIMLIB, char *KEY) ;
void   parse_SIMLIB_GENRANGES(char **WDLIST) ;
void   parse_SIMLIB_IDplusNEXPOSE(char *inString, int *IDEXPT, int *NEXPOSE) ;

int    regen_SIMLIB_GENRANGES(void); // regenerate after reading SIMLIB header
int    check_SIMLIB_GENRANGE(double *GENRANGE_ORIG, double *GENRANGE_NEW);

void   ENDSIMLIB_check(void);

void   get_user_input(void);       // top function to get user inputs
void   set_user_defaults(void);    // set INPUTS.xxx defaults
void   set_user_defaults_SPECTROGRAPH(void);
void   set_user_defaults_RANSYSTPAR(void);
void   set_GENMODEL_NAME(void);
int    read_input_file(char *inFile, int keySource);    
int    parse_input_key_driver(char **WORDLIST, int keySource); // Jul 20 2020

void   parse_input_GENPOP_ASYMGAUSS(void);
void   parse_input_GENZPHOT_OUTLIER(char *string);
void   parse_input_GENZPHOT_FUDGEMAP(char *string);
void   parse_input_FIXMAG(char *string);
int    parse_input_RANSYSTPAR(char **WORDS, int keySource );
int    parse_input_HOSTLIB(char **WORDS, int keySource );
int    parse_input_SIMLIB(char **WORDS, int keySource );
int    parse_input_GENMODEL_ARGLIST(char **WORDS, int keySource );
int    parse_input_GENMODEL(char **WORDS, int keySource );
int    parse_input_NON1ASED(char **WORDS, int keySource );
void   parse_GENMAG_SMEAR_MODELNAME(void);
int  parse_input_KEY_PLUS_FILTER(char **WORDS, int keySource, char *KEYCHECK,
				 float *VALUE_GLOBAL,float *VALUE_FILTERLIST);
int    parse_input_SOLID_ANGLE(char **WORDS, int keySource);
// xxx mark void   parse_input_FIELDLIST(void);

int    parse_input_RATEPAR(char **WORDS, int keySource, char *WHAT,
			   RATEPAR_DEF *RATEPAR );
int    parse_input_ZVARIATION(char **WORDS, int keySource);
int    parse_input_SIMGEN_DUMP(char **WORDS, int keySource);
int    parse_input_SIMSED(char **WORDS, int keySource);
int    parse_input_SIMSED_PARAM(char **WORDS);
int    parse_input_SIMSED_COV(char **WORDS, int keySource );
void   parse_input_SIMSED_SUBSET(char *PARNAME, char *STRINGOPT);

bool   keyContains_SIMSED_PARAM(char *KEYNAME);

int    parse_input_LCLIB(char **WORDS, int keySource );
int    parse_input_CUTWIN(char **WORDS, int keySource );
int    parse_input_GRIDGEN(char **WORDS, int keySource);
int    parse_input_LENS(char **WORDS, int keySource );
int    parse_input_TAKE_SPECTRUM(char **WORDS, int keySource, FILE *fp );
int    parse_input_SPECTRUM(char **WORDS, int keySource);
void   expand_TAKE_SPECTRUM_MJD(float *MJD_RANGE);
int    parse_input_GENMAG_SMEAR_SCALE(char **WORDS, int keySource );
int    parse_input_GENMAG_SMEARPAR_OVERRIDE(char **WORDS, int keySource );

// xxx to do ...
int    parse_input_MWEBV(char **WORDS, int keySource );
int    parse_input_GENMAG_OFF(char **WORDS, int keySource );
int    parse_input_GENMAG_SMEAR(char **WORDS, int keySource );
// xxxx

void   parse_input_OBSOLETE(char **WORDS, int keySource );
bool   valid_DNDZ_KEY(char *WHAT, int keySource, char *KEYNAME );


void   checkVal_GENGAUSS(char *varName, double *val, char *fromFun ) ;

void   sim_input_override(void) ;  // parse command-line overrides
void   prep_user_input(void);      // prepare user input for sim.
void   prep_user_cosmology(void);
void   prep_user_CUTWIN(void);
void   prep_user_SIMSED(void);
void   prep_dustFlags(void);
void   prep_GENPDF_FLAT(void);

void   prep_genmag_offsets(void) ;
void   prep_RANSYSTPAR(void); 
void   prep_RANSYSTPAR_LEGACY(void); 
void   pick_RANSYSTFILE_WILDCARD(char *wildcard, char *keyName, char *randomFile);
void   genmag_offsets(void) ;
void   prioritize_genPDF_ASYMGAUSS(void);
void   compute_lightCurveWidths(void);

void   prep_simpath(void);
int    get_NON1A_MODELFLAG(char *GENVERSION) ;

void   genperfect_override(void);
void   gen_event_driver(int ilc);    // generate RA, DEC, Z, PEAKMJD, etc.
void   gen_event_reject(int *ILC, SIMFILE_AUX_DEF *SIMFILE_AUX,
			char *REJECT_STAGE );
void   gen_event_stronglens(int ilc, int istage);
void   gen_random_coord_shift(void);
void   gen_filtmap(int ilc);  // generate filter-maps
void   gen_modelPar(int ilc, int OPT_FRAME);
void   gen_modelPar_SALT2(int OPT_FRAME);
void   gen_modelPar_SIMSED(int OPT_FRAME);
double pick_gridval_SIMSED(int ipar, int ipar_model);
void   gen_modelPar_dust(int OPT_FRAME);
double gen_MWEBV(double RA, double DEC);       // generate MWEBV
void   override_modelPar_from_SNHOST(void) ;

void pick_NON1ASED(int ilc,
                   INPUTS_NON1ASED_DEF *INP_NON1ASED,
                   GENLC_NON1ASED_DEF *GEN_NON1ASED) ;
// ----------

void   genran_modelSmear(void); // gen randoms for mag-smearing per filter.
double gen_redshift_cmb(void);   // gen ran z from dN/dz
double gen_redshift_helio(void);   // translate zCMB -> zHEL + vPEC
double gen_redshift_helio_OBSOLETE(double zcmb, double RA, double DEC, double vpec) ;
void   gen_redshift_LCLIB(void);

double gen_peakmjd(void);
double gen_peakmjd_smear(void);
void   gen_zsmear(double zerr);
void   genshift_risefalltimes(void);

double gen_dLmag (double zCMB, double zHEL );
void   gen_distanceMag(double zCMB, double zHEL, double *MU, double *lensDMU);

double genz_hubble(double zmin, double zmax, RATEPAR_DEF *RATEPAR );

void   init_RATEPAR ( RATEPAR_DEF *RATEPAR ) ;
void   set_RATEPAR(int ilc, INPUTS_NON1ASED_DEF *INP_NON1ASED ) ;
double SNrate_model(double z, RATEPAR_DEF *RATEPAR ) ;
double SNcount_model(double zmin, double zmax, RATEPAR_DEF *RATEPAR ) ;
double GALrate_model(double b, double l, RATEPAR_DEF *RATEPAR ) ;

double genz_wgt(double z, RATEPAR_DEF *RATEPAR ) ;

void init_simFiles(SIMFILE_AUX_DEF *SIMFILE_AUX);
void update_simFiles(SIMFILE_AUX_DEF *SIMFILE_AUX);
void end_simFiles(SIMFILE_AUX_DEF *SIMFILE_AUX);

void update_accept_counters(void);

void    simEnd(SIMFILE_AUX_DEF *SIMFILE_AUX);
double  gen_AV(void);          // generate AV from model

double  GENAV_WV07(void);
double  gen_RV(void);          // generate RV from model
void    gen_conditions(void);  // generate conditions for each field

//int    gen_PEAKMAG_SPEC_TRIGGER(void); // call GENMAG_DRIVER for peak only
int    gen_TRIGGER_PEAKMAG_SPEC(void); // call GENMAG_DRIVER for peak only
int    gen_TRIGGER_zHOST(void);        // evaluate zHOST trigger early

void   GENMAG_DRIVER(void);    // driver to generate true mags
void   GENFLUX_DRIVER(void);   // driver to generate observed fluxes
void   GENFLUX_DRIVER_LEGACY(void);   // driver to generate observed fluxes
void   set_GENFLUX_FLAGS(int ep);
void   gen_fluxNoise_randoms(void);
void   gen_fluxNoise_calc(int ep, int vbose, FLUXNOISE_DEF *FLUXNOISE);
void   gen_fluxNoise_fudge_diag(int ep, int vbose, FLUXNOISE_DEF *FLUXNOISE);
void   gen_fluxNoise_fudge_cov(int icov);
void   gen_fluxNoise_apply(int ep, int vbose, FLUXNOISE_DEF *FLUXNOISE);
void   dumpLine_fluxNoise(char *fnam, int ep, FLUXNOISE_DEF *FLUXNOISE);
void   dumpEpoch_fluxNoise_apply(char *fnam, int ep, FLUXNOISE_DEF *FLUXNOISE);
void   dumpCovMat_fluxNoise(int icov, int NOBS, double *COV);
void   monitorCov_fluxNoise(void);
void   check_crazyFlux(int ep, FLUXNOISE_DEF *FLUXNOISE);

void   GENSPEC_DRIVER(void);    // driver to generate all spectra for event
void   GENSPEC_MJD_ORDER(int *imjd_order); // order to generate spectra
bool   GENSPEC_PRESCALE_REJECT_SN(void) ;
bool   DO_GENSPEC(int imjd);
void   GENSPEC_INIT(int opt, int imjd);  // init arrays
void   GENSPEC_OBSFLUX_INIT(int imjd, int ILAM_MIN, int ILAM_MAX) ;
void   GENSPEC_TRUE(int imjd);  // generate true MAGs and FLUXes
void   GENSPEC_HOST_CONTAMINATION(int imjd);
void   GENSPEC_TEXPOSE_TAKE_SPECTRUM(int imjd);
double GENSPEC_SMEAR_LEGACY(int imjd, double LAMMIN, double LAMMAX );
double GENSPEC_SMEAR(int imjd, double LAMMIN, double LAMMAX );
double GENSPEC_OBSFLUX_RANSMEAR(int imjd, double OBSFLUXERR, double ERRFRAC_T,
				double *GAURAN_T) ;
void   GENSPEC_FLAM(int imjd);
void   GENSPEC_LAMSMEAR(int imjd, int ilam, double GenFlux );

void   GENSPEC_LAMSMEAR_LEGACY(int imjd, int ilam, double GenFlux,
			       double GenFluxErr, double GenFluxErr_T );

void   GENSPEC_LAMOBS_RANGE(int INDX, double *LAMOBS_RANGE);
double GENSPEC_PICKMJD(int OPT, int INDX, double z,
		       double *TOBS, double *TREST );
void   GENSPEC_FUDGES(int imjd);
void GENSPEC_VERIFY_PEAKMAG(int ifilt_obs, double *GENFLUX_LIST) ;

void   genmodel(int ifilt_obs, int inear);   // generate model-mags
void   genmodelSmear(int NEPFILT, int ifilt_obs, int ifilt_rest,
		     double z, double *epoch, double *genmag, double *generr);

double genmodelSmear_interp(int ifilt_interp, int iep);
double genmodel_Tshift(double T, double z);

void   init_simvar(void);        // one-time init of counters, etc ..
void   init_genmodel(void);      // init above
void   init_genSpec(void);        // one-time init for SPECTROGRAPH
void   init_genSEDMODEL(void); // generic init for SEDMODEL
void   init_kcor_legacy(char *kcorFile);
void   init_kcor_refactor(void);
void   init_covar_mlcs2k2(void);    // init GENLC.COVAR array
void   init_zvariation(void);      // z-dependent sim parameters
void   init_hostNoise(void) ;
void   update_PARDEF_ZVAR(char *parName); // update PARDEF_ZVAR
void   init_CIDRAN(void);
void   sort_CIDRAN(void);
void   init_modelSmear(void);
void   dump_modelSmearSigma(void);
void   init_genSmear_filters(void);
double genSmear_ERRSCALE(double *generr, int iep, int NEPFILT);
void   check_model_default(int indx_model);
void   checkpar_SIMSED(void);
double get_zvariation(double z, char *parname);
GENGAUSS_ASYM_DEF  get_zvariation_GENGAUSS(double z, char *parName,
					   GENGAUSS_ASYM_DEF *GENGAUSS_LOCAL );
void   cp_zvariation(char *outFile_zvar);
void   genmag_boost(void);
void   genmag_MWXT_fromKcor(void);   // apply MW extinct for rest-frame models

void   LOAD_SEARCHEFF_DATA(void);
void   LOAD_SEARCHEFF_DATA_LEGACY(void);

void   gen_spectype(void);

void   setz_unconfirmed(void);
double zHEL_WRONGHOST(void);


int    gen_cutwin(void);
int    gen_cutwin_PEAKMAG(int OPT, int ifilt_obs);
int    geneff_calc(void);
void   magdim_calc(void);
void   screen_update(void);
void   set_screen_update(int NGEN);

int  setEpochGrid( double TMIN, double TMAX, double *TGRID);
void interpEpochGrid(int NEP_LC, double *TList_LC, int NGRID,
                     double *TList, double *magList, double *magerrList );

void   init_RANDOMsource(void);    // init stuff for randoms
void   init_GENLC(void);
int    fudge_SNR(void);

int    NEPFILT_GENLC(int opt, int ifilt_obs);
void   malloc_GENFILT(void);

void   dmp_event(int ilc);
void   dmp_trace_main(char *string, int ilc);
void   snlc_to_SNDATA(int FLAG) ;
void   hostgal_to_SNDATA(int FLAG, int ifilt_obs);
void   MWEBVfluxCor_to_SNDATA(int epoch) ;

void   sprintf_GENGAUSS(char *string, char *name,
			GENGAUSS_ASYM_DEF *genGauss);

void   append_SNPHOT_TEXT(void);  // photometry
void   append_SNSPEC_TEXT(void);  // spectra

void   SORT_GENMJD(void);
void   init_RateModel(void);
void   init_DNDZ_Rate(void) ; // extraGalactic rate vs. redshift
void   init_DNDB_Rate(void) ; // Galactic Rate vs. l & b

int  GENRANGE_CUT(void);
int  GENMAG_CUT(void);

void DUMP_GENMAG_DRIVER(void);

void DASHBOARD_DRIVER(void);
void SIMLIB_DUMP_DRIVER(void);
void SIMLIB_DUMP_openTable(int LDMP_MJD_TEXT,int LDMP_ROOT) ;
void SIMLIB_DUMP_makeTable(int LDMP_MJD_TEXT,int LDMP_ROOT) ;
void zero_SIMLIB_DUMP(SIMLIB_DUMP_DEF *SIMLIB_DUMP) ;
void update_SIMLIB_DUMP_AVGALL(int OPT);

void MJDGAP(int N, double *MJDLIST,  double MJDGAP_IGNORE,
	    double *GAPMAX, double *GAPAVG ) ;

void test_fortran(void);
void test_igm(void);
void test_ran(void);
void test_PARSE_WORDS(void);
void test_zcmb_dLmag_invert(void);

void wr_HOSTLIB_info(void);    // write hostgal info
void wr_SIMGEN_FITLERS(char *path);
void wr_SIMGEN_DUMP(int OPT_DUMP, SIMFILE_AUX_DEF *SIMFILE_AUX);
void wr_SIMGEN_YAML(SIMFILE_AUX_DEF *SIMFILE_AUX);
void rewrite_HOSTLIB_DRIVER(void);

int  MATCH_INDEX_SIMGEN_DUMP(char *varName ) ;
void PREP_SIMGEN_DUMP(int OPT_DUMP );
void PREP_SIMGEN_DUMP_TAKE_SPECTRUM(void);
int  doReject_SIMGEN_DUMP(char *rejectStage) ;

// functions for instrinsic scatter (from R.Biswas, Jul 27, 2011)
void ZERO_COVMAT_SCATTER(void);
void INIT_COVMAT_SCATTER( void );
void PARSE_COVMAT_SCATTER(char * CKEY, char *CVAL);
int  GEN_COVMAT_SCATTER ( double * scatter, double *randoms );

void INIT_FUDGE_SNRMAX(void);

extern void init_snvar__(int *IERR);

extern void rdkcor_(char *kcorFile, int *IERR, int len);

extern double get_maglc8__(int *ifilt, double *t8, double *z8, double *av8);

extern double get_mwxt8__(int *ifilt, double *t8, double *z8, double *av8,
			  double *mwebv8, double *RV, int *OPT_COLORLAW );

extern double get_kcor8__(int *ifilt_rest, int *ifilt_obs,
			  double *t8, double *z8, double *av8 );

extern double get_avwarp8__(double *t8, double *z8,
			    double *mag8_1, double *mag8_2,
			    int *ifilt_1, int *ifilt_2, int *istat );

extern void get_filtlam__(int *opt_frame, int *ifilt,
			  float *lamavg, float *lamrms,
			  float *lammin, float *lammax ) ;

extern void get_kcor_info__(int *NKCOR, double *RV, int *OPT_MWCOLORLAW );

extern void   rdxtpar_(char *xtDir, int len);
extern void   init_xthost__( int *opt );
extern double get_snxt8__(int *opt, int *ifilt,
			 double *T8, double *AV8, double *RV8);

extern int nearest_ifilt_rest__( int *opt, int *ifilt_obs, int *rank, float *z, float *lamdif_min );


extern int filtindx_(char *cfilt, int len);

extern int get_filtmap__ ( char *copt, float *filtmap, int len );


extern void get_filttrans__(int *maskFrame, int *ifilt,
			    char *survey_name, char *filter_name,
			    double *magPrim, int *NLAM, double *lam,
			    double *TransSN, double *TransREF,
			    int len1, int len2);

extern void set_survey__ ( char *name, int *NFILTDEF, int *IFILTDEF,
			   float *LAMSHIFT, int len  );

extern int  get_idsurvey__(char *survey, int len );

extern double kcorfun8_ ( int *ifilt_obs, int *ifilt_rest,
			  double *mag_rest, double *lamdif,
			  double *Trest, double *Z, double *AVwarp ) ;

extern void  get_primary__(char *primary, int *NLAM,
			   double *lam, double *primFlux, int len);

// -----------------------------
//   genmag_xxx functions
// -----------------------------

int gen_smearMag  ( int epoch, int VBOSE );
int npe_above_saturation ( int epoch, double flux_pe);

int init_genmag_mlcs2k2(char *version, char *covFile, float scale_covar,
		     char *filtlist );

int genmag_mlcs2k2(int ifilt, double delta, int nobs,
                double *rest_dates, double *rest_magval, double *rest_magerr);

int gencovar_mlcs2k2(int matsize, int *ifilt, double *rest_epoch,
                  double *covar );


int init_genmag_snoopy( char *modelPath, int optmask, char *filtlist);

int genmag_snoopy(int ifilt, double dm15, int nobs, double *rest_dates,
		  double *rest_mags, double *rest_magerrs );

int init_genmag_stretch0 ( char *templateFilename, char *filtersystem);

int init_genmag_stretch (
                         double H0_ini
                         ,double OMEGAM_ini
                         ,double OMEGAL_ini
                         ,double W0_ini
                         ,double RV_ini
                         ,char  *templateFilename
                         ,char  *filtersystem
                         );


int init_genmag_SALT2(char *model_version, char *model_extrap_latetime,
		      int OPTMASK);

void genmag_SALT2(int OPTMASK, int ifilt,
		  double *parList_SN, double *parList_HOST, double mwebv,
		  double z, double z_forErr, int nobs, double *Tobs,
		  double *magobs_list, double *magerr_list );

double SALT2x0calc(double alpha, double beta, double x1, double c,
		   double dlmag);
double SALT2mBcalc(double x0);

int init_genmag_SIMSED(char *version, char *PATH_BINARY,
		       char *SURVEY, char *kcorFile, int OPTMASK );

void genmag_SIMSED(int OPTMASK, int ifilt, double x0,
		   int NLUMIPAR, int *iflagpar, int *iparmap, double *lumipar,
		   double RV_host, double AV_host,
		   double mwebv, double z, int nobs, double *Tobs,
		   double *magobs_list, double *magerr_list, int *index_sed );

// ------------------------------------
// generic functions for SEDMODELs

int init_filter_SEDMODEL(int ifilt, char *filter_name, char *survey_name,
			 double magprimary, int NLAM, double *lam,
			 double *transSN, double *transREF, double lamshift);

int init_primary_SEDMODEL(char *refname, int NLAM,
			  double *LAMLIST, double *FLUXLIST );

int IFILTSTAT_SEDMODEL(int ifilt_obs, double z) ;

double gridval_SIMSED(int ipar,  int ibin);
double nearest_gridval_SIMSED (int ipar, double lumipar );

// ========== END OF FILE ============
