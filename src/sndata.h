/*****************************************************
 Created Mar 7, 2006 by R.Kessler


 Define SNDATA  structures to be 
 used by programs that generate SNDATA files
  (including simulation)

  Jan 31 2019: MXTYPE -> 1000 (was 200)
  Aug 13 2019: MXPATHLEN -> 300 (was 200)
  Apr 05 2021: MXSPECTRA -> 200 (was 40)
  May 27 2021: MXSPECTRA -> 300 (was 200)
  Jun 04 2021: MXEPOCH -> 5000 (was 2000)
  Apr 24 2022: define HOSTGAL_PROPERTY_xxx [moved from sntools_host.h]
  May 19 2023: MXSPECTRA -> 500 (was 300)
  May 24 2023: MXSPECTRA -> 1000 (was 500)
  May 24 2023: MXDOCLINE -> 2000 (was 1000) to handle long list of NON1A keys
  Aug 05 2023: MXEPOCH -> 15000 (was 5000) for Roman+LSST catalog sim
  Sep 30 2023: MXSPECTRA -> 2000 (was 1000); for Roman+LSST SED TRUE
  Jul 26 2024: define TRANSIENT_NAME and move MXCHAR_CCID here (from sntools_output.h)
  Nov 20 2024: defie MXCHAR_FIELDNAME and MXCHAR_FILTNAME

*****************************************************/

#define MXEPOCH  30000   // max Nepoch for data & simlib

#define MXEPCOV  112     // max epochs to store in covariance matrix

#define MXFIELD_OVP  10  // max number of overlap fields (Feb 2021)
#define MXFILT_COVAR  9  // max number of filters per obs.
#define MXFILTINDX  100  // max filter index
#define MXIDSURVEY 200   // max number of SURVEYS in SURVEY.DEF file
#define MXSPECTRA  2000   // max number of spectra in data files

#define MXDOCLINE 2000    // max number of lines in README.DOC file
#define MXTYPE    1000    // max TYPE id in data base
#define MXBRIGHT  20     // max number of bright times (for MJD ranges)
#define MXVAR_PRIVATE 40 // max number of private variables
#define MXHOSTGAL      2 // max number of matched hosts to write out
#define MXHOSTGAL_PROPERTY 10 // max number of host properites;e.g. logmass
#define MXVAR_HOSTGAL 100 // max number of host params to write out Alex Gagliano 09/2021
#define MXBIN_ZPHOT_Q 101 // max number of quantile percent bins (0,1,2 ...100)
#define MXIMG_STRONGLENS 8  // max number of strong lens images per lens
#define ZEROPOINT_FLUXCAL_DEFAULT 27.5

#define WRITE_MASK_LCMERGE       2  // idem to write lcmerge data files.
#define WRITE_MASK_SIM_SNANA     4  // idem to write SNANA-SIM  
#define WRITE_MASK_SIM_MAGOBS    8  // write data-like with SIM_MAGOBS only 
#define WRITE_MASK_SIM_SNRMON   16  // write SNR(MAGMONITOR)
#define WRITE_MASK_SIM_MODELPAR 32  // write model par for SIMSED, LCLIB
#define WRITE_MASK_SPECTRA     128  // write spectra (Oct 14 2021)
#define WRITE_MASK_ATMOS       256  // write RA,DEC,AIRMASS per obs, for atmos cor
#define WRITE_MASK_SIM_TEMPLATEMAG 512 // write template mags for recurring events
#define WRITE_MASK_SED_TRUE   1024     // write true SED instead of spectra (sim only)
#define WRITE_MASK_COMPACT      64  // suppress non-essential PHOT output
#define WRITE_MASK_COMPACT_noFLUXCAL  4096 // internally set if SMEARFLAG_FLUX=0

#define OPT_ZPTSIG_TRUN  1   // option to use ZPTSIG from template
#define OPT_ZPTSIG_SRUN  2   // idem for search run


#define FAKEFLAG_DATA    0   // SNDATA[isn].FAKE = 0 for data
#define FAKEFLAG_FAKES   1   //                  = 1 for survey fakes
#define FAKEFLAG_LCSIM   2   //                  = 2 for lcsim
#define FAKEFLAG_LCSIM_BLINDTEST 3  // blind-test sim

// strings for DATATYPE key in FITS header (Feb 2021)
#define DATATYPE_DATA        "DATA"       // real data
#define DATATYPE_SIM_SNANA   "SIM_SNANA"  // SNANA sim
#define DATATYPE_SIM_MAGOBS  "SIM_MAGOBS" // e.g., fakes

#define RV_MWDUST 3.1               // RV for MilkyWay

// define generic host properties e.g. LOGMASS, LOGsSFR, COLOR...                                
#define HOSTGAL_PROPERTY_BASENAME_LOGMASS  "LOGMASS"
#define HOSTGAL_PROPERTY_BASENAME_LOGSFR   "LOGSFR"
#define HOSTGAL_PROPERTY_BASENAME_LOGsSFR  "LOGsSFR"
#define HOSTGAL_PROPERTY_BASENAME_COLOR    "COLOR"

#define HOSTGAL_PROPERTY_NAME_LIST HOSTGAL_PROPERTY_BASENAME_LOGMASS " " HOSTGAL_PROPERTY_BASENAME_LOGSFR " " HOSTGAL_PROPERTY_BASENAME_LOGsSFR " " HOSTGAL_PROPERTY_BASENAME_COLOR

// May 2023 define PySEDMODEL names here (instead of in genmag_PySEDMODEL.h)    
//   so that both sim and analysis codes have access and the analysis           
//   (LCFIT) code need not link to python   
#define MODEL_NAME_BYOSED   "BYOSED"
#define MODEL_NAME_SNEMO    "SNEMO"
#define MODEL_NAME_PYBAYESN "PYBAYESN" // don't confuse with C code BAYESN
#define MODEL_NAME_AGN      "AGN"
#define NCHOICE_PySEDMODEL  4
char    PySEDMODEL_CHOICE_LIST[NCHOICE_PySEDMODEL][20] ;


// Mar 2024:  define mask bits for new MASK_REDSHIFT_SOURCE variables to make
//  clerar what is/isn't available. Immediate need is to choose between
//  using quantiles or zspec, and not be fooled by crazy zphot with tiny error.
#define MASK_REDSHIFT_SOURCE_ZHOST_SPEC      1
#define MASK_REDSHIFT_SOURCE_ZSN_SPEC        2
#define MASK_REDSHIFT_SOURCE_ZHOST_PHOT     4
#define MASK_REDSHIFT_SOURCE_ZHOST_QUANTILE  8

// ------------------------------------------

//  disk pointers defined in init_SNDATA

#define MXCHAR_CCID       20  // should be same as MXCHAR_CCID in snana.car
#define MXCHAR_FIELDNAME  20
#define MXCHAR_FILTNAME   20

#define MXPATHLEN            300  // max length of path of full file-name
#define MXLEN_VERSION         72  // max length of VERSION name
#define MXLEN_VERSION_PREFIX  52  // max len of prefix in data or sim version

#define PREFIX_ZPHOT_Q  "ZPHOT_Q" // for zphot quantiles
#define STRING_NZPHOT_Q "NZPHOT_Q"

char PATH_SNDATA_ROOT[MXPATHLEN];        // top dir for SN data
char PATH_SNDATA_PHOTOMETRY[MXPATHLEN];
char PATH_SNDATA_LCMERGE[MXPATHLEN];
char PATH_SNDATA_SIM[MXPATHLEN];
char PATH_SNANA_DIR[MXPATHLEN];
char PATH_USER_INPUT[MXPATHLEN]; // Jan 31 2020
   
struct SURVEY_INFO_DEF {
  int  NSURVEYDEF ;  // number of surveys in SURVEY.DEF file
  char SURVEYDEF_LIST[MXIDSURVEY][40];  // SURVEY-string vs. IDSURVEY
  int  SURVEYFLAG[MXIDSURVEY]; // status of use in survey or field group.
  double  geoLAT[MXIDSURVEY], geoLON[MXIDSURVEY]; // instrument location on Earth
  double  geoALT[MXIDSURVEY];                   // telescope altitude, meters
  double  sin_geoLAT[MXIDSURVEY], cos_geoLAT[MXIDSURVEY];
  double  pressure_atmos[MXIDSURVEY], temperature_atmos[MXIDSURVEY];
  double  pwv_atmos[MXIDSURVEY];
} SURVEY_INFO ;


struct SNDATA_FILTER {
  int  NDEF;
  int  MAP[MXFILTINDX] ;   // gives ifilt_obs for sparse index
  char LIST[MXFILTINDX] ;  // filter list; i.e, ugriz
} SNDATA_FILTER ;  // observer filters


char FLUXUNIT[8] ; // default flux unit is raw ADU

// define structure with info about SN photometry version

struct VERSION
{
  int   N_SNLC;                            // Number of SN in SNDATA struct
  char  NAME[2*MXLEN_VERSION];             // name of version
  char  PREFIX[2*MXLEN_VERSION_PREFIX];    // filename prefix

  int  N_SNFILE;           // number of SN files ( >= N_SNLC)
  char SNLIST_FILE[MXPATHLEN];   // full name of SN list file
  char README_FILE[MXPATHLEN];   // full name of README file
  char IGNORE_FILE[MXPATHLEN];   // added Jun 2012

  int  NLINE_README;         // number of lines of README doc
  int  NLINE_README_INIT;    // subset for initialization (used by sim)
  char README_DOC[MXDOCLINE][2*MXPATHLEN]; 

  int NEPOCH_TOT;          // total number of epochs in all SN

  // xxx mark  int GENFRAME_SIM ;   // sim only: 1=rest+boost(=> KCOR), 2=obs

} VERSION_INFO ;



// define main SNDATA data structure

struct SNDATA {

  // global stuff
  char SNANA_VERSION[20] ; // Feb 2021
  char DATATYPE[20];       // e.g., DATA, SIM_SNANA ...

  // name of SURVEY and SUBSURVEY
  char SURVEY_NAME[40];       // SDSS, SNLS, LSST, etc ...
  char SUBSURVEY_NAME[40];    // e.g., LOWZ_ALL(CFA3) --> CFA3 is subsurvey
  char SUBSURVEY_LIST[MXPATHLEN] ; // optional list in global simlib header
  int  SUBSURVEY_FLAG ;
  
  bool  WRFLAG_BLINDTEST ;  
  bool  WRFLAG_PHOTPROB ;
  bool  WRFLAG_SKYSIG_T ;
  bool  WRFLAG_ATMOS ;         // include RA,DEC,AIRMASS per obs (May 2023)
  bool  WRFLAG_SPECTRA ;       // write spectra

  int   APPLYFLAG_MWEBV;           // T=> correct FLUXCAL
  int   MASK_FLUXCOR;     // indicates SNANA fudges applied to flux[err]
  char  VARNAME_SNRMON[40];

  int   SIMLIB_MSKOPT ;                // mask of options (Dec 2015)
  char  SIMLIB_FILE[MXPATHLEN];        // name of simlib file
  char  CALIB_FILE[MXPATHLEN];
  char  SPEC_FILE[MXPATHLEN];

  int   PHOTFLAG_DETECT; // mask in PHOTFLAG for detection

  // ---- SN-dependent stuff -------
  char SNFILE_INPUT[MXPATHLEN];

  char snfile_output[MXPATHLEN];   // name of data file (no path)
  char SNFILE_OUTPUT[MXPATHLEN];   // full name of data file

  char AUXHEADER_FILE[MXPATHLEN];  // extra info from file to dump into header
  char AUXHEADER_LINES[MXFILTINDX][100]; // extra lines to dump into header
  int  NLINES_AUXHEADER ;

  // start with info read from input flux file

  char  CCID[MXCHAR_CCID];   // char name of SNID
  int   CID ;             // candidate id  
  int   FAKE ;            // 1=FAKE, 0=DATA
  int   NEPOCH;           // total NEPOCH including peak and unused filters
  int   NOBS ;            // total Num of observations (<= NEPOCH)

  // list of observations to store; they pass select_MJD_SNDATA func
  int   NOBS_STORE;
  int   OBS_STORE_LIST[MXEPOCH];

  int   SUBSAMPLE_INDEX ; // if user-input NSUBSAMPLE_MARK > 0

  // Note that for non-SDSS surveys, NEPOCH_NEWMJD = NEPOCH and
  // EPOCH_NEWMJD is just an incremental integer list from 1:NEPOCH.
  // For SDSS, NEPOCH_NEWMJD = NEPOCH/5 and EPOCH_NEWMJD skips by 5.

  int  NEWMJD;                         // # epochs with NEW MJD 
  int  EPOCH_RANGE_NEWMJD[MXEPOCH][2];   // epoch-range for each NEW MJD

  float  AIRMASS[MXEPOCH];
  float  RA[MXEPOCH], DEC[MXEPOCH];   // coord per epoch, degrees
  float  dRA[MXEPOCH], dDEC[MXEPOCH]; // shift w.r.t. band-avg, arcsec
  float  COORDRES[MXEPOCH];           // for both RA & DEC, arcsec
  double RA_AVG, DEC_AVG;  // wgted-average among all RA/DEC (need double precision)
  double RA_AVG_BAND[MXFILTINDX];  // coord-AVG among each band
  double DEC_AVG_BAND[MXFILTINDX];

  float PIXSIZE;                 // pixel size, arcsec
  int   NXPIX, NYPIX;

  int   CCDNUM[MXEPOCH] ; // CCD number or sensor id
  int   IMGNUM[MXEPOCH] ; // 10.13.2021 image number (e.g., EXPNUM, VISIT_ID)

  bool   OBSFLAG_WRITE[MXEPOCH];
  double MJD[MXEPOCH];            // MJD for each epoch

  char  MAGTYPE[20];   // LOG10 or ASINH
  char  MAGREF[20];    // VEGA or AB

  int  SNTYPE;                     // user-defined integer type 
  int  IDTEL[MXEPOCH];             // integer telescope id

  int   FILTINDX[MXEPOCH];        // integer filter indx
  char *FILTNAME[MXEPOCH];        // full name of filter
  char  FILTCHAR_1D[MXEPOCH*2];   // for fortran interface

  int   SEARCH_RUN[MXEPOCH] ;
  int   TEMPLATE_RUN[MXEPOCH] ;
  int   SEARCH_FIELD[MXEPOCH] ;
  int   TEMPLATE_FIELD[MXEPOCH] ;

  float XPIX[MXEPOCH] ;      // X-pixel location
  float YPIX[MXEPOCH] ;      // Y-pixel location
  float EDGEDIST[MXEPOCH] ;      // distance from edge, pixels

  float READNOISE[MXEPOCH] ;   // read noise in e-
  float GAIN[MXEPOCH] ;        // e/ADU
  float FLUX[MXEPOCH] ;        // flux in ADU (per filter/epoch)
  float FLUX_ERRTOT[MXEPOCH] ;  
  float FLUX_ERRTEMPLATE[MXEPOCH] ;  // template contribution to error
  int   NPE_ABOVE_SAT[MXEPOCH];      // Npe above saturation

  // now define stuff computed or found by this program

  char  NAME_IAUC[MXCHAR_CCID];       // public name in IAUC
  char  NAME_TRANSIENT[MXCHAR_CCID];  // internal string name (July 2024)
  char *FIELDNAME[MXEPOCH] ;    // survey field (generalize SDSS STRIPE)
  char  FIELDNAME_1D[MXEPOCH*20] ;  // for fortran interface

  float FLUXCAL[MXEPOCH] ;         // calibrated flux for fitter
  float FLUXCAL_ERRTOT[MXEPOCH] ;  
  float FLUXCAL_ERRTEMPLATE[MXEPOCH] ;  // correlated template error

  float MAG[MXEPOCH] ;            // magnitude (per filter/epoch)
  float MAG_ERR[MXEPOCH] ;    // magnitude error (per filter/epoch)
  float ZEROPT[MXEPOCH] ;         // zero point for template 
  float ZEROPT_ERR[MXEPOCH] ;     // zero point error on mean
  float ZEROPT_SIG[MXEPOCH] ;     // zero point sigma
  float TEXPOSE[MXEPOCH];         // exposure time (sec)
  float SKYSUB_ERR[MXEPOCH] ;     // sky-subtraction error
  float GALSUB_ERR[MXEPOCH] ;     // gal-subtraction error
  int   NPRESN[MXFILTINDX] ;      // number of pre-SN epochs

  int   OPT_ZEROPT_SIG;          // determines trun or srun zptsig
  int    PHOTFLAG[MXEPOCH] ;     // photometry flags (0 => OK)
  float  PHOTPROB[MXEPOCH];      // fit-prob or FoM per epoch
 
  float SKY_SIG[MXEPOCH] ;     // sky noise (ADU/pix)
  float SKY_SIG_T[MXEPOCH] ;   // template sky noise (ADU/pix)
  float PSF_SIG1[MXEPOCH] ;       // PSF sigma of inner gaussian
  float PSF_SIG2[MXEPOCH] ;       // PSF isgma of outer
  float PSF_RATIO[MXEPOCH] ;    // PSF
  float PSF_NEA[MXEPOCH];        // write NEA if > 0 (Feb 28 2021)
  bool  NEA_PSF_UNIT; 
  float MWEBV ;                    // MilyWay Galactic E(B-V)
  float MWEBV_ERR;                 // error on  above
  float MWXT_MAG[MXFILTINDX];      // MWXT mag at mean wave of band

  int     HOSTGAL_USEMASK ;  // bits 0,1,2,3 --> MAGOBS, MAGOBSERR, SB, SBERR
  int     HOSTGAL_NMATCH[2] ; // NMATCH and NMATCH2 (tight/loose DLR cut)
  int     HOSTGAL_FLAG[MXHOSTGAL];    // indicate problems (May 2021)
  long long HOSTGAL_OBJID[MXHOSTGAL] ;            // up to 4 host matches
  float   HOSTGAL_SB_FLUXCAL[MXFILTINDX];   // surface bright (FLUXCAL/asec)
  float   HOSTGAL_SB_FLUXCALERR[MXFILTINDX];  // error on above
  float   HOSTGAL_SB_MAG[MXFILTINDX];      // SB mag in 1 sq-arcsec
  float   HOSTGAL_MAG[MXHOSTGAL][MXFILTINDX];         // host mag
  float   HOSTGAL_MAGERR[MXHOSTGAL][MXFILTINDX];   
  float   HOSTGAL_SNSEP[MXHOSTGAL] ;
  float   HOSTGAL_DDLR[MXHOSTGAL]   ;                 // SNSEP/DLR
  double  HOSTGAL_RA[MXHOSTGAL];
  double  HOSTGAL_DEC[MXHOSTGAL];
  float   HOSTGAL_CONFUSION ;         // note: does NOT depend on each host
  float   HOSTGAL_PHOTOZ[MXHOSTGAL] ;
  float   HOSTGAL_PHOTOZ_ERR[MXHOSTGAL] ;
  float   HOSTGAL_SPECZ[MXHOSTGAL] ;
  float   HOSTGAL_SPECZ_ERR[MXHOSTGAL] ;

  float   HOSTGAL_ZPHOT_Q[MXHOSTGAL][MXBIN_ZPHOT_Q] ;  // redshifts
  int     HOSTGAL_PERCENTILE_ZPHOT_Q[MXBIN_ZPHOT_Q] ;  // percentiles
  int     HOSTGAL_NZPHOT_Q ;

  float  *PTR_HOSTGAL_PROPERTY_TRUE[MXHOSTGAL_PROPERTY];
  float  *PTR_HOSTGAL_PROPERTY_OBS[MXHOSTGAL_PROPERTY];
  float  *PTR_HOSTGAL_PROPERTY_ERR[MXHOSTGAL_PROPERTY];
  float   HOSTGAL_LOGMASS_TRUE[MXHOSTGAL] ;
  float   HOSTGAL_LOGMASS_OBS[MXHOSTGAL] ;  
  float   HOSTGAL_LOGMASS_ERR[MXHOSTGAL] ;
  float   HOSTGAL_LOGSFR_TRUE[MXHOSTGAL] ;
  float   HOSTGAL_LOGSFR_OBS[MXHOSTGAL] ;
  float   HOSTGAL_LOGSFR_ERR[MXHOSTGAL] ;
  float   HOSTGAL_LOGsSFR_TRUE[MXHOSTGAL] ;  
  float   HOSTGAL_LOGsSFR_OBS[MXHOSTGAL] ;                 
  float   HOSTGAL_LOGsSFR_ERR[MXHOSTGAL] ;
  float   HOSTGAL_COLOR_TRUE[MXHOSTGAL] ;
  float   HOSTGAL_COLOR_OBS[MXHOSTGAL] ;
  float   HOSTGAL_COLOR_ERR[MXHOSTGAL] ;

  long long HOSTGAL_OBJID2[MXHOSTGAL] ;
  long long HOSTGAL_OBJID_UNIQUE[MXHOSTGAL] ;
  float   HOSTGAL_ELLIPTICITY[MXHOSTGAL] ;
  float   HOSTGAL_SQRADIUS[MXHOSTGAL] ;
  int     HOSTGAL_NFILT_MAGOBS ;      // NFILT with magobs (fixed number)

  float REDSHIFT_HELIO;         // final (best) redshift, Helio frame
  float REDSHIFT_HELIO_ERR;     // final (best) redshift, Helio frame
  float REDSHIFT_FINAL;         // idem, CMB frame
  float REDSHIFT_FINAL_ERR;     // error on above
  float VPEC, VPEC_ERR;         // Jan 2018
  float LENSDMU, LENSDMU_ERR;   // Aug 2024
  int   REDSHIFT_QUALITYFLAG;   // quality flag: survey dependent meaning
  int   MASK_REDSHIFT_SOURCE;   // Mar 2024: determines source(s) of redshift

  // info obtained during survey (in SQL). 
  // Note that LC fit => Masao's fitter

  float SEARCH_PEAKMJD ;     // approx MJD at g-band peak, from LC fit
  int   SEARCH_TYPE ;        // type from search (i.e., 120=confirmed Ia)

  float MJD_TRIGGER ;      //  MJD when trigger is satisfied (Apr 2017)
  float MJD_DETECT_FIRST ; // mjd of 1st detection
  float MJD_DETECT_LAST;

  // declare generation quantities for simulation (fake flag = 2)

  char SIM_MODEL_NAME[60];   // model name
  int  SIM_MODEL_INDEX;      //integer id for model or class
  int  SIM_TEMPLATE_INDEX ;  // template index for NON1ASED, SIMSED, LCLIB ...
  char SIM_COMMENT[200]; 
  // xxx mark  int  SIM_TYPE_INDEX;        // same as SNTYPE (if set).
  int  SIM_GENTYPE;          // same as SNTYPE (if set).
  char SIM_TYPE_NAME[60];    // Ia, Ib, II, etc ...

  int   SIM_LIBID;
  int   SIM_NGEN_LIBID;
  int   SIM_SEARCHEFF_MASK;    // bits 1,2 => passes software,human eff.
  float SIM_SEARCHEFF_SPEC ;   // spec-search efficiency
  float SIM_SEARCHEFF_zHOST;   // EFF(zHOST) when not spec-confirmed
  float SIM_REDSHIFT_HELIO ;   // for SN
  float SIM_REDSHIFT_CMB   ;   // for SN
  float SIM_REDSHIFT       ; // legacy variable, same as z_CMB
  float SIM_REDSHIFT_HOST  ; // zhelio of host (Jan 2016)
  int   SIM_REDSHIFT_FLAG  ; // indicates source of redshift (4.19.2019)
  char  SIM_REDSHIFT_COMMENT[40];
  float SIM_VPEC ;           // peculiar velocity, km/sec
  float SIM_DLMU ;
  float SIM_LENSDMU ;  
  float SIM_MUSHIFT;           //  May 2025 user-define MU shift (random or fixed)
  float SIM_RA, SIM_DEC ;      //  simulated RA and DEC
  float SIM_PEAKMJD ;          //  peak MJD in g-band
  float SIM_MJD_EXPLODE;       //  MJD of explosion (see OPTMASK_T0SHIFT_EXPLODE)
  float SIM_AVTAU ;            //  
  float SIM_AV, SIM_RV ;       //  host extinction parameters
  float SIM_MWRV ;             //  MilkyWay RV
  float SIM_MWEBV;             //  simulated Milky Way extinction
  int   SIM_NOBS_UNDEFINED ;   //  NOBS for which GENMODEL is not defined
  int   SIMOPT_MWCOLORLAW ;    //  integer option for color law
  int   SIMOPT_MWEBV ;         //  integer option to modify MWEBV_SFD
  int   SIMOPT_FLUXERR ;       //  sim-option to fludge flux errors
  float SIM_MAGSMEAR_COH ;     // coherent part of intrinsic smear

  float SIM_GALFRAC[MXFILTINDX];     //  effective galaxy mag under SN
  float SIM_PEAKMAG[MXFILTINDX] ;    //  peak mag in each filter
  float SIM_TEMPLATEMAG[MXFILTINDX];  //  for LCLIB model only
  float SIM_EXPOSURE_TIME[MXFILTINDX] ;   // relative exposure time
        
  // - - - - - HOSTLIB properties - - - -
  char  HOSTLIB_FILE[MXPATHLEN];       // name of hostlib file (Feb 2014)
  int   SIM_HOSTLIB_MSKOPT ;          // non-zero => simulate HOSTLIB
  int   NPAR_SIM_HOSTLIB;             // number of host params
  char  SIM_HOSTLIB_KEYWORD[MXVAR_HOSTGAL][60];   // keyword for ascii
  char  SIM_HOSTLIB_PARNAME[MXVAR_HOSTGAL][40];   // name of host params to store
  float SIM_HOSTLIB_PARVAL[MXVAR_HOSTGAL][MXHOSTGAL];   // host param values per neighbo

  float SIM_ANGSEP_GROUPID;  // ang-sep between SIMLIB and HOSTLIB coord (GROUPID only)

  long long SIM_HOSTLIB_GALID ; // true HOST GALID -> OBJID
  //  float     SIM_HOSTLIB_DDLR  ; // true DDLR

  // - - - - -
  float SIM_RISETIME_SHIFT;    // rise time shift relative to model
  float SIM_FALLTIME_SHIFT;
  float SIM_TRESTMIN ;
  float SIM_TRESTMAX ;

  float SIM_WGT_POPULATION; // product of population wgts (color x stretch x AV)

  // luminosity params
  float SIM_STRETCH ;  // for history lesson
  float SIM_DELTA ;    // for mlcs
  float SIM_DM15 ;     // for Snoopy
  float SIM_THETA ;    // for BayeSN
  float SIM_SALT2alpha ;
  float SIM_SALT2beta ;
  float SIM_SALT2gammaDM;    //  mag shift from host-SN correlations
  float SIM_SALT2x0;
  float SIM_SALT2x1;
  float SIM_SALT2c;
  float SIM_SALT2mB;

  int   SIMFLAG_COVMAT_SCATTER ; // 0 or 1 => off or ON
  float SIM_COVMAT_SCATTER[3];
  char  SIM_COVMAT_SCATTER_NAME[3][40];

  int   NPAR_SIMSED;
  char  SIMSED_KEYWORD[100][80]; // SIMSED_[PARAM,param,GRIDONLY ...]
  char  SIMSED_PARNAME[100][40]; // SIMSED parameter names
  float SIMSED_PARVAL[100];      // SIMSED parameter values

  // Python SED models: BYOSED, SNEMO ...
  int   NPAR_PySEDMODEL;
  char  PySEDMODEL_KEYWORD[100][80]; // keyword for data files
  char  PySEDMODEL_PARNAME[100][40]; // parameter names
  float PySEDMODEL_PARVAL[100];      // parameter values
  char  PySEDMODEL_NAME[60];         // e.g., BYOSED, SNEMO
  
  int   NPAR_LCLIB;
  char  LCLIB_KEYWORD[100][40]; // LCLIB keyword for data file
  char  LCLIB_PARNAME[100][40]; // LCLIB parameter names
  float LCLIB_PARVAL[100];      // LCLIB parameter values

  // strong lens info (July 20 2019)
  int    SIM_SL_FLAG;    // strong lens flag
  int    SIM_SL_IDLENS;  // lens ID in library
  long long SIM_SL_GALID; // hostlib GALID matched to lens properties
  int    SIM_SL_NIMG;    // number of strong lens images
  int    SIM_SL_IMGNUM;  // image num [0 to NIMG-1]
  double SIM_SL_MINSEP ; // separation from nearest image (Jul 2022)
  double SIM_SL_zLENS, SIM_SL_TDELAY, SIM_SL_MAGSHIFT;
  double SIM_SL_LOGMASS, SIM_SL_LOGMASS_ERR;
  double SIM_SL_XIMG, SIM_SL_YIMG ;

  // simulation quantities for each epoch
  float SIMEPOCH_TREST[MXEPOCH] ;       // Trest - Tpeak, days
  float SIMEPOCH_TOBS[MXEPOCH] ;        // Tobs - Tpeak, days
  float SIMEPOCH_MAG[MXEPOCH] ;         // generated obs mag
  float SIMEPOCH_MAGSMEAR[MXEPOCH] ;    // intrinsic mag-smear
  float SIMEPOCH_AVWARP[MXEPOCH] ;      // AVwarp parameter
  float SIMEPOCH_KCORVAL[MXEPOCH] ;     // KCOR value
  char  SIMEPOCH_KCORNAM[MXEPOCH][8] ;  // name of KCOR
  float SIMEPOCH_WARPCOLVAL[MXEPOCH] ;     // warp color value
  char  SIMEPOCH_WARPCOLNAM[MXEPOCH][8] ;  // warp color name (i.e, B-V)
  float SIMEPOCH_SNRMON[MXEPOCH];          // SNR of monitor mag

  float SIMEPOCH_DCR[MXEPOCH];        // arcsec
  float SIMEPOCH_DCR_dRA[MXEPOCH];    // arcsec
  float SIMEPOCH_DCR_dDEC[MXEPOCH];   // arcsec
  float SIMEPOCH_DCR_dMAG[MXEPOCH];   // mag shift
  int   MAGMONITOR_SNR;            // transferred from INPUTS.MAGMONITOR

  // private variables (Nov 24, 2012)
  int     NVAR_PRIVATE ;
  char    PRIVATE_VARNAME[MXVAR_PRIVATE][60];
  char    PRIVATE_KEYWORD[MXVAR_PRIVATE][60];
  double  PRIVATE_VALUE[MXVAR_PRIVATE] ;

  // biasCor mask that is internally set by sim; 
  // designed to inform analysis codes
  // +1=alphaGrid, +2=betaGrid ...
  int  SIM_BIASCOR_MASK; 

} SNDATA ;


// end of file.
