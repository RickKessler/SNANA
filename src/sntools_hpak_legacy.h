/***********************************************************
  sntools_hpak_legacy.h

  Legacy SNLCPAK (light curve pack) and SPECPAK (spectra pack)
  global declarations, structs, and function prototypes.

  Nov 09 2026: REFACTOR ONLY -- moved verbatim from sntools_output.h
               in preparation for removing SNLCPAK/SPECPAK.
               No functional changes.

  This file is included from sntools_output.h (so that other codes
  that include sntools_output.h, e.g. psnid_BEST.c, still see the
  prototypes) and from sntools_output.cpp.
  It must be included after sndata.h (MXFILTINDX, MXCHAR_CCID, 
  MXSPECTRA).

  Note that the globals below are intentionally defined (not extern)
  in this header, exactly as they were in sntools_output.h.
***********************************************************/

#ifndef SNTOOLS_HPAK_LEGACY_H
#define SNTOOLS_HPAK_LEGACY_H

#define MXSPEC_SPECPAK MXSPECTRA

// -------------------------------------------------
// SNLCPAK global declarations for light curves

#define SNLCPAK_INITDONE  12345
#define SPECPAK_INITDONE  98765
#define MXFIT_PER_SN      10

bool SNLCPAK_USE_ROOT  ;
bool SNLCPAK_USE_HDF5  ;
bool SNLCPAK_USE_TEXT  ;
bool SNLCPAK_USE_NPZ   ;
int NCALL_SNLCPAK_FILL ; // number of calls to SNLCPAK_FILL

bool SPECPAK_USE_ROOT  ;
bool SPECPAK_USE_HDF5  ;
bool SPECPAK_USE_TEXT  ;
bool SPECPAK_USE_NPZ   ;
bool SPECPAK_USE_MARZ  ;
int NCALL_SPECPAK_FILL ; // number of calls to SPECPAK_FILL

// define flags for epoch info
#define SNLCPAK_EPFLAG_FLUXDATA    1  // data flux 
#define SNLCPAK_EPFLAG_FLUXMODEL   11 // fit-model flux 
#define SNLCPAK_EPFLAG_REJECT      2  // for data only
#define SNLCPAK_EPFLAG_CHI2        3  // for fitted data only
#define SNLCPAK_EPFLAG_FITFUN      4
#define SNLCPAK_EPFLAG_FLUXSIM     5  // simulated flux (May 20 2016)
#define SNLCPAK_EPFLAG_FLUXREST    6  // k-corrected rest-frame fluxes
#define SNLCPAK_EPFLAG_KCOR        7  // k-cor (Jan 2020)
#define SNLCPAK_EPFLAG_AVWARP      8  // AVwarp (Jan 2020)
#define SNLCPAK_EPFLAG_SIMFLUXREST 9  // simulated rest-flux
#define MXFLAG_SNLCPAK_EPOCH      10  

// define flags for filter info (one entry per filter)
#define SNLCPAK_BANDFLAG_NDOF      100
#define SNLCPAK_BANDFLAG_PKFLUX    101
#define SNLCPAK_BANDFLAG_PKMJD     102
#define SNLCPAK_BANDFLAG_CHI2      103
#define MXFLAG_SNLCPAK_FILTER      4

#define MXTEXT_SNLCPAK 5  // max number of displayText strings


// define structure containing information for 1 light curve.
struct  SNLCPAK_OUTPUT {
  int NCALL ;
  int NFIT_PER_SN ;  // number of fits per SN

  // define user-passed global variables set just once from SNLCPAK_INIT
  char SURVEY[60] ;
  char VERSION_PHOTOMETRY[200] ;
  char VERSION_SNANA[40];
  char SURVEY_FILTERS[MXFILTINDX] ;

  char TEXT_FORMAT[20]; 
  int  OPT_TEXT_FORMAT ;  // integer index for above (-9 -> no text output)

  // set user-passed variables that change with each SN
  char CCID[40];

  int  NTEXT ;         // number of text strings
  int  MXTEXT ;        // to init text in tree
  char DISPLAYTEXT[MXTEXT_SNLCPAK][80] ;  // display strings

  int      NOBS_MAX ;      // max of all NOBS (for malloc)
  int      NOBS[MXFLAG_SNLCPAK_EPOCH] ;  // NOBS for each SNLCPAK_FLAG
  double  *MJD[MXFLAG_SNLCPAK_EPOCH] ;   // pointer to dynam allocate
  double  *TOBS[MXFLAG_SNLCPAK_EPOCH] ;  // pointer to dynam allocate
  double  *EPDATA[MXFLAG_SNLCPAK_EPOCH] ;
  double  *EPDATA_ERR[MXFLAG_SNLCPAK_EPOCH] ;
  int     *IFILTOBS[MXFLAG_SNLCPAK_EPOCH] ;

  double  FILTDATA[MXFLAG_SNLCPAK_FILTER][MXFILTINDX];  
  double  FILTDATA_ERR[MXFLAG_SNLCPAK_FILTER][MXFILTINDX];  

  int  NLCPAK_TOT ;  // expected number of LC to PACK 
  int  NLCPAK ;      // incremental counter
  int  SIMFLAG ;     // =1 for simulated fluxes (May 2016)

  // internally defined variable
  int  INITDONE  ;
  int  NFILTDEF_SURVEY ;
  int  NFILTDEF_PLOT ;             // Number of filters to plot  
  char FILTLIST_PLOT[MXFILTINDX] ; // char-list of plotted filters
  int  IFILT_MIN, IFILT_MAX;       // min,max sparse filter index


  int  USEFILT[MXFLAG_SNLCPAK_EPOCH][MXFILTINDX];   // vs. sparse IFILT indx
  int  NOBS_FILT[MXFLAG_SNLCPAK_EPOCH][MXFILTINDX] ;  // idem per filter
  int  *IFILT[MXFLAG_SNLCPAK_EPOCH] ;               // sparse IFILT

  char  TEXTFILE_PREFIX[200];

} SNLCPAK_OUTPUT ;


// Apr 2019: create struture for multiple spectra per event (SN & HOST)
struct SPECPAK_OUTPUT {
  int  NSPEC;  // per event
  char SURVEY[60] ;
  char VERSION_PHOTOMETRY[200] ;
  int  INITDONE ;

  char TEXT_FORMAT[20]; 
  int  OPT_TEXT_FORMAT ;  // integer index for above (-9 -> no text output)

  // set user-passed variables that change with each SN
  char    CCID[MXCHAR_CCID];
  int     NLAMBIN_TOT;    // sum of lam bins for each event
  int     ID_LIST[MXSPEC_SPECPAK];
  int     NLAMBIN_LIST[MXSPEC_SPECPAK];
  double  MJD_LIST[MXSPEC_SPECPAK];
  double  TOBS_LIST[MXSPEC_SPECPAK]; 
  double  TEXPOSE_LIST[MXSPEC_SPECPAK];

  int    *ID ;
  double *LAMMIN, *LAMMAX, *FLAM, *FLAMERR ;

} SPECPAK_OUTPUT ;


#ifdef __cplusplus
extern "C" {
#endif

  // SNLCPAK functions

  void SNLCPAK_INIT(char *SURVEY, char *VERSION_PHOT, char *VERSION_SNANA,
		    char *SURVEY_FILTERS, int SIMFLAG, 
		    int NFIT, char *TEXT_FORMAT );

  void snlcpak_init__(char *SURVEY, char *VERSION_PHOT, char *VERSION_SNANA,
		      char *SURVEY_FILTERS, int *SIMFLAG, 
		      int *NFIT, char *TEXTFMT ); 

  void SNLCPAK_NFIT_PER_SN  (int  NFIT_PER_SN);
  void snlcpak_nfit_per_sn__(int *NFIT_PER_SN);

  void SNLCPAK_SURVEY(void);  // called after output file is opened
  void snlcpak_survey__(void);

  void SNLCPAK_DISPLAYTEXT(char *CCID, char *DISPLAYTEXT) ;
  void snlcpak_displaytext__(char *CCID, char *DISPLAYTEXT) ;

  void SNLCPAK_NFIT(int NFIT) ;
  void snlcpak_nfit__(int *NFIT) ;
  
  void SNLCPAK_DATA(char *CCID, int NOBS, double *MJD, double *TOBS, 
		    double *DATA, double *DATA_ERR, int *IFILTOBS, int FLAG);
  void snlcpak_data__(char *CCID, int *NOBS, double *MJD, double *TOBS, 
		      double *DATA, double *DATA_ERR, int *IFILTOBS,int *FLAG);
  
  void SNLCPAK_FILL(char *CCID) ;
  void snlcpak_fill__(char *CCID) ;
  
  void SNLCPAK_CLEAR_SN(void); 
  void SNLCPAK_CLEAR_PLOT(void); 
  void snlcpak_clear_plot__(void);
  
  // internal wrappers that do NOT need mangled fortran fun.
  void SNLCPAK_FILL_PREP(void);
  void SNLCPAK_CHECK(char *CCID, char *callFun);

  // SPECPAK functions (Apr 2019)
  void SPECPAK_INIT(char *SURVEY, char *VERSION_PHOT, char *TEXT_FORMAT );
  void specpak_init__(char *SURVEY, char *VERSION_PHOT, char *TEXTFMT ); 
  void SPECPAK_CLEAR_PLOT(void); 
  void specpak_clear_plot__(void);

  void SPECPAK_DATA(char *CCID, int IDSPEC, double MJD, double Tobs, 
		    double Texpose,int NLAMBIN,double *LAMMIN, double *LAMMAX,
		    double *FLAM,double *FLAMERR);
  void specpak_data__(char *CCID,int *IDSPEC, double *MJD, double *Tobs,
		      double *Texpose, 
		      int *NLAMBIN, double *LAMMIN, double *LAMMAX,
		      double *FLAM,double *FLAMERR);

  void SPECPAK_FILL(char *CCID) ;
  void specpak_fill__(char *CCID) ;

#ifdef __cplusplus
}
#endif

#endif  // SNTOOLS_HPAK_LEGACY_H

// END
