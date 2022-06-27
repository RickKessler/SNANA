/* ==============================================================
 Created July 2012 by M. Sako and R.Kesser

 Functions to initialize and execute fitting
 for the psnid 'BEST 'method, where
     BEST = "Bayesian Evidence from Supernova Templates"

 and is described in
 Sako et al, ApJ, 738, (2011). 
 http://adsabs.harvard.edu/abs/2011ApJ...738..162S

 Calling functions are in psnid.cra.

?!?!?!??!?!?!?!??!?!?!??!?!?!?!??!?!?!??!?!?!?!??!?!?!??!?!?!?!??!?!?!?
  ISSUES found by RK Feb 28 2013;

  code does NOT abort when IIN is a template but is not defined in 
  PSNID_CHARTYPE_LISTS. Why not abort (Search Q1) ?

  chir < 2 cut is hard-wired, but should be a namelist parameter.
  There is also a cut P_BAYES>0.9 for something. All cuts should
  be nml CUT_XXX parameters with reasonable defaults. (see Q2).

  Some parts of the code are fragile and won't work if another
  global type is added such as pec-Ia or Ix. Explicit boolean
  logic and else-if blocks need to be replaced by functions
  that loop over all non-Ia types. See, for example, LIA_Ibc,LIA_II

?!?!?!??!?!?!?!??!?!?!??!?!?!?!??!?!?!??!?!?!?!??!?!?!??!?!?!?!??!?!?!?

 --------------------------------
 Integration notes for Masao:

   - to get started:
     cp $SNANA_DIR/src/psnid_BEST.c .
     cp $SNANA_DIR/src/psnid_tools.c .
     cp $SNANA_DIR/src/psnid.cra .
     cp /home/s1/rkessler/snana_debug/psnid/Makefile .
     make
   - fix NDOF to only count used epochs with valid model-mags
     -> this is a problem for a grid search; the best fit almost
        always favors a model with no overlap with data (MS).

   - psnid_fit.c -> psnid_BEST.c
     -> DONE (MS)

   - config parameters now read from namelist file.
     see instructions at top of psnid.cra, and private
     keys at the end of
        /home/s1/rkessler/snana_debug/psnid/PSNID_SDSS.nml

   - use only the filter for which
        PSNID_INPUTS.USEFILT[ifilt_obs] = 1
     so that filter choice can be controled via namelist
     input FILTLIST_FIT = 'griz'
     -> DONE (MS)

   - BEST_PAR, BEST_TYPE, BEST_ITYPE are bad names because
     'BEST' does not really mean the acronym for your method.
     Please use a different name such as 'FINAL' or 'FIT'
     --> DONE (MS)

   - Move additional 'general-use' functions into psnid_tools.c[h], 
     where 'general-use' means that it works on any method.
     Your MCMC utility probably goes into psnid_tools.c.

   - please start replacing your 'psnid' prefix with 'psnid_best'.
     The psnid prefix is for general tools in psnid_tools.c.

   - The following params make your code rather fragile:
       #define PSNID_ITYPE_SNIA     0
       #define PSNID_ITYPE_SNIBC    1
       #define PSNID_ITYPE_SNII     2
     It would be better to get the list directly from SNGRID.
     A user-namelist string can define a subset, or a list to
     ignore.

   - change ISEED1, ISEED2 to a more general name to avoid
     possible conflict with internal variables elsewhere.

   - Suggestion for comment string:
       FIT-FUN=[TYPE] for z-prior=[zprior]  P(Ia,Ibc,II) = PIa,PIbc,PII     
     so that we see all Bayesian probs on each plot.

   - write results to output text-file in 
         PSNID_BEST_INIT_OUTFILE
         PSNID_BEST_UPDATE_OUTFILE


   - to make light curve plots,
       mkfitplots.pl  --h <hisFileName>  --tmin -20 --tmax 80

     HISTORY
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Dec 02 2012  RK :  
  * fix compile bug found by Ishida (in psnid_best_set_cid)
  * write generic comment-info at top of fitres file;
     see PSNID_BEST_INIT_OUTFILE(char *outFile)
   
  Dec 19, 2012 RK
     init SIMSTRING_NAMES[200] = ""  (caught by valgrind)

  Jan 4, 2013 RK
     replace ISTAT with ERRFLAG in PSNID_BEST_DOFIT().
     Define example ERRFLAG_FAIL and set ERRFLAG = ERRFLAG_FAIL,
     but more explicit ERRFLAG_XXX values should be defined.
     This ERRFLAG return code goes into the snana ntuple
     (with LTUP_SNANA=T).

  Jan 7 ,2013 RK
    Fix 64 bit bugs in 
      void psnid_best_model_alloc()
      void psnid_best_model_free()
    by replacing '1' arguments with long ONE8;

  Jan 21, 2013
    bug gix by M.Sako to implemtn MCMC_NSTEP > 0.

  Jan 23, 2013 RK - print 'Skip' message with ERRFLAG.

  Jan 27, 2013 RK - new function psnid_strTypeMatch(ITYPE,CTYPE,0)
                    and char PSNID_CHARTYPE_LISTS[][].

  Feb 8, 2013 RK (for v10_24b)
     Add new plotting interface:  SNLCPLOT_PSNID_BEST()
     For this version, both the old and new interfaces work.
     Next version the old interface will be removed.

     move psnid_best_pogson2flux* functions to psnid_tools.c
     since these functions are of general use.
     They are re-named 'psnid_pogson2flux* functions (without 'best').

  Feb 21, 2013 RK - fix to compile with c++ .

  Feb 27, 2013 RK - 
     fix NONIA DISPLAYSTRING to show sim-index and name instead of 
     Masao's internal index. Still need to fix LUMIPAR_Ibc and LUMIPAR_II
     in the output-fitres file because they show the local internal index
     that cannot be traced to the original template.
     See how INDX_SIM is computed.

  Feb 28, 2013 RK
    - add IIN to PSNID_CHARTYPE_LISTS.

  Mar 04, 2013 RK 
    - fix string-len args to getsim_string().
      Define new parameter MXLEN_SIMSTRING 400.

    - add and load CHI2_FILT array in SNLCPLOT_PSNID_BEST().

  Mar 7 2013 MS - fix NaN-output bugs.

  Mar 10, 2013 RK - 
     add new module  MONPLOT_PSNID_BEST() with examples of how
     to add arbitrary monitor plots. This function is meant
     to be modified, while SNLCPLOT_PSNID_BEST() is a standard 
     function that should be left alone.
             
  Mar 27, 2013 RK -
     for SIM, tack on SIM_ITYPE in fitres file. Get char-type
     from SNANA function GET_SIMNAME_TYPE, and then convert to
     PSNID-integer ITYPE. This translation is needed because
     the snana code that reads the char-type (Ia,Ib,II ...)
     cannot determine the integer type based on PSNID convention.


     PSNID_BEST_RESULTS.FINAL_PAR -> PSNID_BEST_RESULTS.GRID_FINAL_PAR -
     and define 
     PSNID_BEST_RESULTS.FINAL_PARVAL & PARERR


  Mar 30 2013: MS - tweak PSNID_MODEL_MAGERR[m][type_count][j][k+1]  

  Mar 31, 2013 RK
    - remove PSNID_BEST_INIT_OUTFILE & PSNID_BEST_UDPATE_OUTFILE
    - add    PSNID_BEST_INIT_OUTPUT  & PSNID_BEST_UDPATE_OUTPUT
        (includes both SNTABLE for hbook/root & text-output)

    - new function psnid_best_define_VARNAMES(..) defines variables
      for text-ouptut, SNTABLE, and NN. Thus all variables are defined
      in only one place.


  Apr 2, 2013:
     Central filter wavelengths are in SNGRID_PSNID[itype].FILTER_LAMAVG.
     -> use these an CCMextinct to compute Galactic extinction for model.

     psnid_best_nearnbr() does basic nearest-nbr analysis with hard-wired
     coefficients.

  Apr 4, 2013: naming convention: lumi -> shape and LUMI -> SHAPE

  Apr 8 2013
    - add Galactic extinction in psnid_best_mwxt (RK)
    - few model-error tweaks (M.S.)

  Apr 21 2013
    - major bug fix.  PSNID_PEAK_GUESS was not initialized correctly.

  Apr 23 2013
    - added LIGHTCURVE_QUALITY parameter to output


  Jul 2 2013:
    - major updates/bug fix (MS)
    - in psnid_best_store_finalPar(), translate internal NONIA_INDEX
      to absolute index.

  Jul 19, 2013: add MJD argument to SNLCPAK

  Jul 26, 2013: Dillon Brout, added psnid nml cuts for PBayes and Fitprob
    - Default Fitprobs cuts are 0.01
    - Default PBayes_Ia cut is 0.0
    - Default PBayes_Ibc and II cuts are 0.9
    - This is essentially what was hardcoded before, but with a few added
    conditions.

 Aug 22 2013: RK - fix bugs for plotting.

 Sep 7 2013 RK :
  - use new wrapper "psnid_best_modelErr" to compute model errors.

  - new  function psnid_best_ratePrior() to return relative rate.
    Default rates are all 1 unless namelist OPT_RATEPRIOR = 1.
     [also see changes in psnid.car]

  - little cleanup for speed; runs about x2 faster.

  - dLmag: replace hard-wired args with PSNID_INPUTS.[H0,OMAT,OLAM,W0]
     passed from &SNLCINP INPUTS H0_REF, OMAT_REF, OLAM_REF, W0_REF.


 Nov 1 2013: new option to include peakmags in ascii & hbook/root files.
      See &PSNIDINP variable FILTLIST_PEAKMAG_STORE = 'ri' .
      This is NOT automatic to avoid too many variables for
      surveys like JPAS that would results in 3 x 56 = 168 extra
      variables.  Default is '', or no peak mags.

 Feb 5 2014 MS : new option to reject outlier light curve points.
  - new keywords are CHISQMIN_OUTLIER and NREJECT_OUTLIER, default to
      CHISQMIN_OUTLIER = 999.0
      NREJECT_OUTLIER  = 0
  - see psnid_best_flag_outlier function.

 Mar 7 2014 MS : major bug fix in MCMC mode
  - fixed problem psnid_best_run_mcmc; redshift prior was not
    used properly.
  - also added new model error.

 Mar 21 2014 MS : few more bux fixes in MCMC routine
  - mu step was not initialized properly.
  - changed center of mu range from LDCM to best-fit value from
    grid search.
  - better memory management in MCMC routine.

  May 20 2014 RK
    PSNID_BEST_INIT_OUTPUT ->
    PSNID_BEST_INIT_SNTABLE &  PSNID_BEST_INIT_OUTFILE_LEGACY
    where the latter is only for the legacy outfile defined
    by &PSNIDINP variable FITRES_DMPFILE. The nominal output
    is now controled by &SNLCINP string SNTABLE_LIST.

  Sep 10 2015 RK
    Fix bug in  PSNID_BEST_GET_FITFUN() calling hunt.
    See trest1 argument to hunt() call.

  Feb 2017 (RK) 
    +  add new PEC1A type (4th class)
    +  lots of refactoring 
    +  mark lots of 'fragile alerts' that need to be cleaned up
       (e..g, explicit if-block over each NON1A type, instead of 
         looping or using indexed array)
 
  Aug 2017 (RK)
    + updates to allow using NON1A templates WITHOUT SNIa templates.
    + output FITRES table includes only types for which templates
      are specified.

  Jan 6 2018: add extra z-arg in dLmag. Should be zCMB and zHEL,
              but here it doesn't matter which z.

  July 25 2018: define PSNID_MINOBS, and abort if too few obs.

  May 20 2019: PSNID_NONIA_MXTYPES = 100 ->  1000 (for lots of KN models)

  Feb 26 2020: + print out warning message when skipping candidate
                 due to user-specified OPT_ZPRIOR

  May 22 2020: + added error model of H20 templates (S13_H20)
  Oct 23 2020: use HzFUN_INFO and new sntools_cosmology.c[h]

  May 23 2022 RK 
    + add Ic-BL for V19 templates
    + abort if a type is undefined

 ================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_sf_gamma.h>

#include "fitsio.h"
#include "sntools.h"              // snana stuff
#include "sntools_cosmology.h"    // dLmag
#include "sntools_modelgrid.h"
#include "sntools_output.h"
#include "sntools_nearnbr.h"
#include "psnid_tools.h"


// Delcare functions
// Upper case name are called by fortran/snana program psnid
// Lower case names are internal.

void PSNID_BEST_INIT(void);
void psnid_best_init__(void) {
  PSNID_BEST_INIT();
}


void PSNID_BEST_INIT_SNTABLE(int OPT, char *TEXTFORMAT, int LSIM) ;
void psnid_best_init_sntable__(int *OPT, char *TEXTFORMAT, int *LSIM) 
{ PSNID_BEST_INIT_SNTABLE(*OPT,TEXTFORMAT,*LSIM); }


void PSNID_BEST_TABLEFILE_COMMENTS(void) ; // store comments for table

void PSNID_BEST_UPDATE_OUTPUT(char *CCID) ;
void psnid_best_update_output__(char *CCID) {
  PSNID_BEST_UPDATE_OUTPUT(CCID) ;
}


int PSNID_BEST_DOFIT(char *CCID, int NOBS, int *IFILTOBS, double *MJD,
		     double *FLUXDATA, double *FLUXERR, double *FLUXSIM,
		     double *REDSHIFT, double *REDSHIFT_ERR, 
		     double MWEBV, double MWEBVERR, int SIM_NON1A_INDEX);
int psnid_best_dofit__(char *CCID, int *NOBS, int *IFILTOBS, double *MJD,
		       double *FLUXDATA, double *FLUXERR, double *FLUXSIM,
		       double *REDSHIFT, double *REDSHIFT_ERR, 
		       double *MWEBV, double *MWEBVERR, int *SIM_NON1A_INDEX) 
{
  return PSNID_BEST_DOFIT(CCID, *NOBS, IFILTOBS,MJD,FLUXDATA,FLUXERR,FLUXSIM,
			  REDSHIFT, REDSHIFT_ERR, *MWEBV, *MWEBVERR,
			  *SIM_NON1A_INDEX );
}


void PSNID_BEST_SUMMARY(void);
void psnid_best_summary__(void) {
  PSNID_BEST_SUMMARY();
}


extern void getsim_string__(int *MODEL_INDEX, int *NVAR, 
			    char *STRING, char *fromWho, int len1, int len2);

extern void get_simname_type__(char *name, int len);

extern void init_table_snanavar__(int *ID_TABLE, char *BLOCK, int *IFLAG, 
				  int len ) ;

extern void init_table_simvar__(int *ID_TABLE, char *BLOCK, int len);

extern int snana_nearnbr_rdinput__(void);  // Apr 16 2013 - RK

/************************************************************************/
/*************  BEGIN: Masao's routines and global variables ************/

#define PSNID_NITER          3

#define PSNID_ZPRIOR_FLAT    0
#define PSNID_ZPRIOR_SNLC    1
#define PSNID_ZPRIOR_HOST    2
#define PSNID_NZPRIOR        3

#define PSNID_ITYPE_SNIA     0
#define PSNID_ITYPE_SNIBC    1
#define PSNID_ITYPE_SNII     2
#define PSNID_ITYPE_PEC1A    3  // Feb 13 2017
#define PSNID_ITYPE_MODEL1   4  // Aug 28 2017
#define PSNID_ITYPE_MODEL2   5  // Aug 28 2017
#define PSNID_ITYPE_MODEL3   6  // Aug 28 2017
#define PSNID_ITYPE_MODEL4   7  // Aug 28 2017
#define PSNID_NTYPES         8  //


#define MXLEN_SIMSTRING   400 // RK Mar 2013

#define PSNID_TABLE_ID    7788     // RK - for SNTABLE
#define PSNID_TABLE_NAME "FITRES"  // idem

#define PSNID_BEST_NEARNBR_NVAR 3 // 3 variables for NN analysis
#define IGNORE_for_NEARNBR -9

#define H0_PSNID 72.0/(1.e6*PC_km) // Sep 2013 RK

#define PSNID_MINOBS 2   // July 25 2018

// ---------------------------------------------------
// Jan 27, 2013 RK - define char-lists here in the header
//                   instead of hard-wired in the code
// Feb 13 2017 RK - add IIb to list (for Jones 2017 paper)
// May 23 2022 RK - add Ic-BL to handle V19 templates
char PSNID_CHARTYPE_LISTS[PSNID_NTYPES][80] =   
  { 
    "Ia" ,                          // itype = PSNID_ITYPE_SNIA 
    "Ib  Ic   Ibc Ic-BL" ,          // itype = PSNID_ITYPE_SNIBC
    "II  IIP  IIL  IIn  IIN IIb",   // itype = PSNID_ITYPE_SNII
    "PEC1A-Iabg",                   // itype = PSNID_ITYPE_PEC1A
    "MODEL1" ,     // generic MODEL 1
    "MODEL2" ,
    "MODEL3" ,
    "MODEL4"       // generic MODEL 4
  } ;

char PSNID_ITYPE_STRING[PSNID_NTYPES][8] =
  { "Ia", "Ibc", "II", "PEC1A", "MODEL1", "MODEL2", "MODEL3", "MODEL4" } ; 

// ---------------------------------------------------

// PSNID fit parameter index
#define PSNID_PARAM_LOGZ     0
#define PSNID_PARAM_SHAPEPAR 1
#define PSNID_PARAM_COLORPAR 2
#define PSNID_PARAM_COLORLAW 3
#define PSNID_PARAM_TMAX     4
#define PSNID_PARAM_DMU      5
#define PSNID_NPARAM         6

#define PSNID_MCMC_NLIMITS   5

char PSNID_TYPE_NAME[PSNID_NTYPES][10];

#define PSNID_NONIA_MXTYPES   1000  //max number of templates per class
int PSNID_NONIA_ABSINDEX[PSNID_NTYPES][PSNID_NONIA_MXTYPES];

int PSNID_PARAM_MAX_INDEX[PSNID_NPARAM];

int PSNID_NFILTER, PSNID_MAXND ;
int PSNID_MAXNZ, PSNID_MAXNL, PSNID_MAXNA, PSNID_MAXNU;
double PSNID_ZMIN, PSNID_ZSTEP;
double PSNID_LMIN, PSNID_LSTEP;
double PSNID_AMIN, PSNID_ASTEP;
double PSNID_UMIN, PSNID_USTEP, PSNID_UMAX ;
double PSNID_TBIN ;
int PSNID_MUFINE;
int PSNID_MAXNL_NONIA;
int PSNID_NGRID[PSNID_NTYPES];

double PSNID_FITPROB_CUTLIST[PSNID_NTYPES]; // Feb 2017 RK
double PSNID_PBAYES_CUTLIST[PSNID_NTYPES];  // Feb 2017 RK

int PSNID_ZRBN[PSNID_NITER], PSNID_ZWID[PSNID_NITER],
  PSNID_ARBN[PSNID_NITER], PSNID_AWID[PSNID_NITER],
  PSNID_URBN[PSNID_NITER], PSNID_UWID[PSNID_NITER];
double PSNID_TSTART[PSNID_NITER], PSNID_TSTOP[PSNID_NITER], PSNID_TSTEP[PSNID_NITER];

int PSNID_FITDMU, PSNID_FITDMU_CC;
int PSNID_USE_AV_PRIOR, PSNID_USE_DM_PRIOR, PSNID_USE_Z_PRIOR;
int PSNID_THIS_TYPE;
double PSNID_PEAK_START, PSNID_PEAK_STOP, PSNID_PEAK_GUESS, PSNID_PEAKMJD;
double PSNID_FIRST_MJD;

// Grid models
double ***PSNID_MODEL_EPOCH, ****PSNID_MODEL_MAG, ****PSNID_MODEL_MAGERR,
  ****PSNID_MODEL_EXTINCT, ****PSNID_MODEL_MWEXTINCT;


double PSNID_BIGN=1.e30, PSNID_SMALLN=1.e-30, PSNID_SPECZSIG=20.0;
double PSNID_BASE_COLOR[PSNID_NTYPES];


struct SIMVAR_PSNID {
  int N, MODEL, ITYPE ;    
} SIMVAR_PSNID ;



/*************************************************************************/
/********************         MCMC variables     *************************/
int MCMC_RUN, MCMC_NSTEP, MCMC_NBURN;
double PSNID_BEST_MCMC_DELTA_Z, PSNID_BEST_MCMC_DELTA_DM,
  PSNID_BEST_MCMC_DELTA_AV, PSNID_BEST_MCMC_DELTA_TMAX,
  PSNID_BEST_MCMC_DELTA_DMU;
double PSNID_BEST_MCMC_DELTA_Z_DEFAULT, PSNID_BEST_MCMC_DELTA_DM_DEFAULT,
  PSNID_BEST_MCMC_DELTA_AV_DEFAULT, PSNID_BEST_MCMC_DELTA_TMAX_DEFAULT,
  PSNID_BEST_MCMC_DELTA_DMU_DEFAULT;
double MCMC_REDSHIFT;
#define PSNID_BEST_ISEED1      -9283
#define PSNID_BEST_ISEED2      -8134
long psnid_idum1, psnid_idum2;


// hard code default MCMC step sizes
#define PSNID_BEST_MCMC_DELTA_Z_DEFAULT       0.0080
#define PSNID_BEST_MCMC_DELTA_DM_DEFAULT      0.0200
#define PSNID_BEST_MCMC_DELTA_AV_DEFAULT      0.0100
#define PSNID_BEST_MCMC_DELTA_TMAX_DEFAULT    0.0300
#define PSNID_BEST_MCMC_DELTA_DMU_DEFAULT     0.0300

// RK define generic error flag as example; 
// MS should define these properly
#define ERRFLAG_FAIL 1 // generic undefined failure


long ONE8 = 1 ;

///////////////////////////////////////////////////////////////////////////
/////      Storing all relevant PSNID results in this structure.      /////
///////////////////////////////////////////////////////////////////////////
typedef struct {

  // RK Mar 28 2013: below is either GRID-search of MCMC values
  //                 so that final values have unique address.
  double  FINAL_PARVAL[PSNID_NZPRIOR][PSNID_NTYPES][PSNID_NPARAM];
  double  FINAL_PARERR[PSNID_NZPRIOR][PSNID_NTYPES][PSNID_NPARAM];
  double  FINAL_PEAKMAG[PSNID_NTYPES][MXFILTINDX] ;
  int     FINAL_NONIA_INDX_SPARSE[PSNID_NTYPES]; // SHAPEPAR for NONIA

  char CID[64]; // RK Aug 2013: [60] -> [64] for double-alignment

  int NZPRIOR;
  double ZPRIOR[PSNID_NZPRIOR];
  double ZPRIOR_ERR[PSNID_NZPRIOR];
  int ZPRIOR_DO[PSNID_NZPRIOR];

  // Grid-search results
  int     NGOOD[PSNID_NZPRIOR][PSNID_NTYPES];
  int     MINCHISQ_IND[PSNID_NZPRIOR][PSNID_NTYPES][PSNID_NPARAM];
  double  TOBSMIN[PSNID_NZPRIOR][PSNID_NTYPES];
  double  TOBSMAX[PSNID_NZPRIOR][PSNID_NTYPES];
  double  MINCHISQ[PSNID_NZPRIOR][PSNID_NTYPES];
  double  FITPROB[PSNID_NZPRIOR][PSNID_NTYPES];
  double  GRID_FINAL_PAR[PSNID_NZPRIOR][PSNID_NTYPES][PSNID_NPARAM];

  double  PBAYESIAN[PSNID_NZPRIOR][PSNID_NTYPES];
  int     LIGHTCURVE_QUALITY[PSNID_NZPRIOR][PSNID_NTYPES];
  int     FINAL_ITYPE[PSNID_NZPRIOR];
  char    FINAL_TYPE[PSNID_NZPRIOR][10];

  // MCMC parameters (median and +/- 1,2 sigma limits)
  double MCMC_FINAL_PAR[PSNID_NZPRIOR][PSNID_NTYPES][PSNID_NPARAM][PSNID_MCMC_NLIMITS];


} PSNID_BEST_RESULTS_DEF;

PSNID_BEST_RESULTS_DEF PSNID_BEST_RESULTS;


// Oct 2013 (RK); define lc-residual structure for each fit 
  RESIDS_PSNID_DOFIT_DEF    RESIDS_PSNID_DOFIT[PSNID_NTYPES] ;
F_RESIDS_PSNID_DOFIT_DEF  F_RESIDS_PSNID_DOFIT ; // for best-type only

///////////////////////////////////////////////////////////////////////////


void psnid_best_define_TableVARNAMES(int DO_ADDCOL, int DO_NN);
void psnid_best_define_TableRESIDS();
void psnid_best_malloc_resids(void);

void psnid_best_reset_results();
void psnid_best_set_cid(char *cid);
void psnid_best_setup_searchgrid();

void psnid_best_model_alloc();
void psnid_best_model_free();

void psnid_best_split_nonia_types(int *types, int optdebug);
void psnid_best_set_grid_limits(int typeindex);
void psnid_best_set_grid_values(int typeindex, int *nonia_types);
void psnid_best_dump_model_lc(int shape_index, int z_index);


void psnid_best_set_zprior(double *zin, double *zinerr);
void psnid_best_set_zprior_onez(double *zin, double *zinerr);
void psnid_best_get_z_grid(double *grid);
void psnid_best_get_l_grid(double *grid);
void psnid_best_get_c_grid(double *grid);
void psnid_best_get_u_grid(double *grid);

void psnid_best_flag_epochs(int nobs, int *data_filt, double *data_mjd,
			    double *data_fluxcal,
			    double *data_fluxcalerr,
			    int *useobs);

double psnid_best_ratePrior(int itype, int indx_template, double z);

void psnid_best_check_for_avprior();
void psnid_best_check_for_dmu();
double psnid_best_avprior1(int itype, double av);

void psnid_best_print_zprior_info(int i);
void psnid_best_store_results(char *CCID, int itype, int z,
			      int ***ind, double **evidence, int optdump);
void psnid_best_store_PBayes(char *CCID,int z,double **evidence,int optdump);
int  psnid_bestType_cuts(int z);
void psnid_best_store_finalPar(void) ;
void psnid_best_store_fitResids(int itype, int **obsflag);

void psnid_best_dump_results();
void psnid_best_calc_fitpar(int itype, int ***ind,
			    double *redshift, double *shapepar, double *colorpar,
			    double *colorlaw, double *pkmjd, double *dmu);

void psnid_best_peak_mjd_guess(int nobs, int *data_filt, double *data_mjd,
			       double *data_fluxcal,
			       double *data_fluxcalerr);
void psnid_best_check_lc_quality(int itype, int z, int ***ind,
				 int nepoch, double *mjd, double *flux,
				 double *fluxerr);
void psnid_best_check_ind_bounds(int ipass, int itype, int ***ind);


void psnid_best_nearnbr(char *CCID); 

/* core psnid functions */
void psnid_best_grid_compare(int itype, int zpind, int nobs,
			     int *data_filt, double *data_mjd,
			     double *data_fluxcal, double *data_fluxcalerr,
			     int *useobs,
			     int ***ind, double **evidence);
void psnid_best_calc_chisq(int nobs, int *useobs,
			   int *data_filt, double *data_mjd,
			   double *data_fluxcal, double *data_fluxcalerr,
			   double *model_mjd,
			   double **model_mag, double **model_magerr,
			   int *ngood, double *chisq, int optDebug,
			   int doReject);

int psnid_best_flag_outlier(int itype, int zpind, int nobs,
			    int *data_filt, double *data_mjd,
			    double *data_fluxcal, double *data_fluxcalerr,
			    int *useobs, int ***ind);

int psnid_strTypeMatch(int itype, char *string, int optDebug ); // RK


// Sep 2013 RK - wrappers for model mag-error.
double psnid_best_modelErr(int TYPEINDX, double Trest, double modelErr_grid);
double psnid_best_modelErr_S11(int TYPEINDX, double Trest) ; 
double psnid_best_modelErr_S13(int TYPEINDX, double Trest) ; 
double psnid_best_modelErr_S14(int TYPEINDX, double Trest) ; 


/**********************************/
/***       MCMC functions       ***/
//void psnid_best_run_mcmc(char *CCID, int itype, int nepoch,
//			 double *data_mjd,
//			 double **data_mag, double **data_magerr,
//			 double ebv,
//			 int zfix, double zprior, double zprior_err);

void psnid_best_init_mcmc();
void psnid_best_run_mcmc(char *CCID, int itype, int nobs,
			 int *data_filt, double *data_mjd,
			 double *data_fluxcal, double *data_fluxcalerr,
			 int *useobs,
			 double mwebv, double mwebverr,
			 int zfix, double zprior, double zprior_err);

void psnid_best_mcmc_zref(int itype, double *zref);

void psnid_best_lc_interp(int itype,
			  double z, double shape,
			  double color, double tmax, double dmu,
			  double ***day,
			  double ****mags, double ****magserr,
			  double ****extinct, double ****mwextinct,
			  double *out_day, double **out_mags,
			  double **out_magserr);

void psnid_best_z_l_weights(double z, double shape, double **w);
void psnid_best_histstats(int ngrid, double *xgrid, int *ygrid,
			  double *moments);
void psnid_best_param_limits(int ngrid, double *xgrid, int *ygrid,
			     double *params);

double psnid_best_mwxtmag(int ifilt) ;

/***************************************/


/* routines to get fit parameters and best-fit models for plotting*/

void SNLCPLOT_PSNID_BEST(int iplot) ;    // standard LC+fit plots
void snlcplot_psnid_best__(int *iplot);  // (RK Feb 8 2013)

void MONPLOT_PSNID_BEST(int iplot) ;    // arbitrary monitor plots
void monplot_psnid_best__(int *iplot);  // (RK Mar 10 2013)


int PSNID_BEST_GET_NFIT(void);
int psnid_best_get_nfit__(void) {
  return PSNID_BEST_GET_NFIT();
}

void PSNID_BEST_GET_FITPAR(char *CCID, int ind,
			   double *redshift, double *shapepar, double *colorpar,
			   double *colorlaw, double *pkmjd, double *dmu,
			   double *chi2, int *ndof,
			   char *plotlabel, char *funLabel);
void psnid_best_get_fitpar__(char *CCID, int *ind,
			     double *redshift, double *shapepar, double *colorpar,
			     double *colorlaw, double *pkmjd, double *dmu,
			     double *chi2, int *ndof,
			     char *plotlabel, char *funLabel) {
  PSNID_BEST_GET_FITPAR(CCID, *ind, redshift, shapepar, colorpar, colorlaw,
			pkmjd, dmu, chi2, ndof, plotlabel, funLabel);
}


void PSNID_BEST_GET_FITFUN(char *CCID, int ind,
			   int ifilt, int nepoch, double *epoch,
			   double *mag, double *magerr);
void psnid_best_get_fitfun__(char *CCID, int *ind,
			     int *ifilt, int *nepoch, double *epoch,
			     double *mag, double *magerr) {
  PSNID_BEST_GET_FITFUN(CCID, *ind, *ifilt, *nepoch, epoch, mag, magerr );
}

int PSNID_BEST_GET_BESTTYPE(void);
int psnid_best_get_besttype__(void) {
  return PSNID_BEST_GET_BESTTYPE();
}

void PSNID_BEST_GET_TYPENAME(int itype, char *name) ;
void psnid_best_get_typename__(int *itype, char *name) {
  PSNID_BEST_GET_TYPENAME(*itype,name);
}

void psnid_best_get_plotind(int ind, int *this_z, int *this_t);
int  doPlot_psnid(int iplot) ;
int  doplot_psnid__(int *iplot) ;

  

/*************    END: Masao's routines and global variables ************/
/************************************************************************/




/************************************************************************/
int PSNID_BEST_DOFIT(char *CCID, int NOBS, int *IFILT, 
		     double *MJD, double *FLUXDATA, double *FLUXERR, 
		     double *FLUXSIM,
		     double *REDSHIFT, double *REDSHIFT_ERR, 
		     double MWEBV, double MWEBVERR, int SIM_NON1A_INDEX)
/*
  PSNID main engine and SNANA/FORTRAN entry point to run the PSNID fit
  for this light curve.

   Inputs:
   CCID               = name of candidate ID (string)
   NOBS               = total number of observations (filter-epochs)
   IFILT[iobs]        = filter index for each filter-epoch  (0..N-1)
   MJD[iobs]          = MJD for each obs
   FLUXDATA,FLUXERR[iobs] = calibrated flux and error (SNANA "fluxcal" units)
   FLUXSIM            = true sim flux (added Jan 2020)
   REDSHIFT           = redshift prior
   REDSHIFT_ERR       = redshift prior error
   MWEBV,MWEBVERR     = Milky Way E(B-V)

   Returns ERRFLAG = 0 if fit converges; 
   returns ERRFLAG > 0 otherwise -> skip it.
  


  Feb 9, 2013 RK -   
    IFILTOBS -> IFILT (variable name) since this is a 
    sparse index that goes from 0 to NFILT-1.
    "IFILTOBS" should be used only for absolute filter indices.

  Jan 2020: pass FLUXSIM

 */
/************************************************************************/
{

  int  i, j, k, z, ERRFLAG, *nonia_types, *useobs, **obsflag;
  int ***minchisq_ind;
  int nOutlier=0;
  double **evidence;
  char fnam[] = "PSNID_BEST_DOFIT" ;

  // ---------- BEGIN -------------

  ERRFLAG = ERRFLAG_FAIL ;  // default is not OK

  printf("\t Begin PSNID fit on CID = %s \n", CCID);
  fflush(stdout);


  if ( NOBS < PSNID_MINOBS ) {
    sprintf(c1err,"NOBS=%d < MINOBS(%d) for CID=%s", 
	    NOBS, PSNID_MINOBS, CCID);
    sprintf(c2err,"Check light curve");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }

  // Feb 8 2013: RK store light curve info for use in other functions
  //             such as dumping or plotting.
  psnid_store_data(CCID, NOBS, IFILT, MJD,
		   FLUXDATA, FLUXERR, FLUXSIM,
		   REDSHIFT, REDSHIFT_ERR,
		   MWEBV, MWEBVERR, SIM_NON1A_INDEX );
  
  /****************************************************************
    Structure:

    - Loop over redshift priors:
      - Loop over types (Ia/Ibc/II):
        - Set up model grid:
        - Grid search:
          - Input  = data and grid of models
          - Output = indices of best-fit parameters and chi-squared
        - Calculate Bayesian evidence
      - Calculate Bayesian probabilities
      - Optionally run MCMC
      - Dump output


    psnid_best_reset_results
    psnid_best_set_cid
    psnid_best_setup_searchgrid
    psnid_best_set_zprior_onez

    Loop over redshifts (but currently only one OPT_ZPRIOR)
      Loop over SN types

        psnid_best_set_grid_limits
        psnid_best_model_alloc
        psnid_best_set_grid_values

        psnid_best_grid_compare
            psnid_best_calc_chisq
        psnid_best_store_results
        psnid_best_run_mcmc

      psnid_best_store_PBayes
    psnid_best_dump_results

  ****************************************************************/



  /****************************************************/
  /*****     Some diagnostic print statements     *****/
  //    psnid_dumpInput_data(CCID, NOBS, IFILT, MJD,
  //    		       FLUX, FLUXERR,
  //    		       REDSHIFT, REDSHIFT_ERR,
  //    		       MWEBV, MWEBVERR);
  // This prints out useful information about the grid
  //  psnid_dumpBins_SNGRID();
  /****************************************************/


  /*******************************************************************/
  /***  Preliminary setup -- read these parameters from name-list file  ***/
  /********************************************************************/

  psnid_best_reset_results();
  psnid_best_set_cid(CCID);


  // split non-Ia types into Ibc and II
  PSNID_MAXNL_NONIA  = 
    SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_SHAPEPAR];
  nonia_types        = ivector(0, PSNID_MAXNL_NONIA);
  psnid_best_split_nonia_types(nonia_types, 0 ); 

  // local variables that keep track of best-fit model
  minchisq_ind = i3tensor(0,PSNID_NITER+1, 0,PSNID_NTYPES, 0,PSNID_NPARAM);
  evidence     = dmatrix(0,PSNID_NZPRIOR, 0,PSNID_NTYPES);

  for (i=0; i<=PSNID_NITER; i++) {
    for (j=0; j<PSNID_NTYPES; j++) {
      for (k=0; k<PSNID_NPARAM; k++) {
	minchisq_ind[i][j][k] = 0;
      }
    }
  }

  useobs  = ivector(0,NOBS);
  obsflag = imatrix(0,PSNID_NTYPES, 0,NOBS);
  for (i=0; i<PSNID_NTYPES; i++) {
    for (j=0; j<=NOBS; j++) {  obsflag[i][j] = 1;  }
  }

  // redshift priors (0=flat, [1,2=other priors])
  psnid_best_set_zprior_onez(REDSHIFT, REDSHIFT_ERR);
  if ( PSNID_BEST_RESULTS.ZPRIOR_DO[0] == 0 )
  {
    printf("\t WARNING:  skipping CID = %s for OPT_ZPRIOR %d\n\n",
	   CCID, PSNID_INPUTS.OPT_ZPRIOR);
    return ERRFLAG ;
  }

  /****************************************/
  /***     Main part of grid search     ***/
  /****************************************/

  // Loop over:
  // 1) redshift priors (0=FLAT, 1=SNLC, 2=ZHOST; see PSNID_ZPRIOR_XXXX)
  // 2) SN tpes         (0=Ia, 1=Ibc, 2=II; see PSNID_ITYPE_SNXX)

  z  =  0 ; // only one z prior

  
  // loop over SN tpes (0=Ia, 1=Ibc, 2=II; see PSNID_ITYPE_SNXX)
  for (i=0; i<PSNID_NTYPES; i++) {

    evidence[z][i] = 0.0 ;

    // skip this type if there are no templates
    if ( PSNID_NGRID[i] == 0 ) { continue ; }

    // set limits, allocate memory, and fill in global model arrays
    psnid_best_set_grid_limits(i);
    psnid_best_model_alloc();
    psnid_best_set_grid_values(i, nonia_types);

    // flag early and late epochs
    for (j=0; j<=NOBS; j++) { useobs[j] = 1 ; }  // use all points by default
    psnid_best_flag_epochs(NOBS, IFILT, MJD, FLUXDATA, FLUXERR, useobs);

    //////////////////////////////////////////////////////////////////
    /////                      GRID SEARCH                       /////
    psnid_best_grid_compare(i, z, NOBS, IFILT, MJD, FLUXDATA, FLUXERR,
			    useobs, minchisq_ind, evidence);

    // Outlier Rejection
    if ( PSNID_INPUTS.NREJECT_OUTLIER > 0 ) {
      // Flag chi2 outlier points using best-fit model of each SN type.
      //  useobs = 1  (good)
      //    "    = 0  (pre/post max)
      //    "    = 2  (chi2 outlier)
      nOutlier = psnid_best_flag_outlier(i, z, NOBS, IFILT,
					 MJD, FLUXDATA, FLUXERR,
					 useobs, minchisq_ind);
      // then refit if any point is flagged
      if ( nOutlier > 0 ) {
	psnid_best_grid_compare(i, z, NOBS, IFILT, MJD, FLUXDATA, FLUXERR,
				useobs, minchisq_ind, evidence);
      }
    }
    // end of outlier rejection
    for (j=0; j<=NOBS; j++) {
      obsflag[i][j] = useobs[j];
    }

    psnid_best_check_lc_quality(i, z, minchisq_ind,
    				NOBS, MJD, FLUXDATA, FLUXERR);
    psnid_best_store_results(CCID, i, z, minchisq_ind, evidence, 0);
    //////////////////////////////////////////////////////////////////
    

    //////////////////////////////////////////////////////////////////
    /////                    MCMC (SN Ia only)                   /////
    if (i==PSNID_ITYPE_SNIA && PSNID_INPUTS.MCMC_NSTEP>0) {
      
      psnid_best_init_mcmc();
      
      psnid_best_run_mcmc(CCID, i, NOBS, IFILT, MJD, FLUXDATA, FLUXERR,
			  useobs, MWEBV, MWEBVERR, z,
			  PSNID_BEST_RESULTS.ZPRIOR[z],
			  PSNID_BEST_RESULTS.ZPRIOR_ERR[z]);
    } // end MCMC block
    //////////////////////////////////////////////////////////////////


    psnid_best_model_free();

  }  // end i loop over PSNID_NTYPES


  psnid_best_store_PBayes(CCID, z, evidence, 0);
  

  ERRFLAG = 0;  // valid redshift prior


  // print results to screen
  psnid_best_dump_results();

  // RK - load FINAL_PARVAL and FINAL_PARERR arrays.
  psnid_best_store_finalPar();

  // store fit resids after ITYPE_BEST and finalPar are determined 
  F_RESIDS_PSNID_DOFIT.NFIT  = 0 ;
  for (i=0; i<PSNID_NTYPES; i++) {
    if ( PSNID_NGRID[i] > 0 )
      { psnid_best_store_fitResids(i, obsflag); }
  }

  // memory cleanup
  free_ivector(nonia_types, 0, PSNID_MAXNL_NONIA);
  free_i3tensor(minchisq_ind, 0,PSNID_NITER+1, 0,PSNID_NTYPES, 0,PSNID_NPARAM);
  free_dmatrix(evidence, 0,PSNID_NZPRIOR, 0,PSNID_NTYPES);

  free_imatrix(obsflag, 0,PSNID_NTYPES, 0,NOBS);
  free_ivector(useobs, 0,NOBS);

  // RK - add message if fit is skipped
  if ( ERRFLAG != 0 ) {
    printf("\t Skip PSNID fit on CID = %s : ERRFLAG=%d\n", CCID, ERRFLAG);
    fflush(stdout);
  }

  // check option for nearest nbr analysis
  if ( ERRFLAG == 0 ) {
    psnid_best_nearnbr(CCID);  // check for nearest nbr analysis
  }


  printf("\n");  fflush(stdout);
  return   ERRFLAG ;

}
// end of PSNID_BEST_DOFIT



// ===================================
void psnid_best_store_finalPar(void) {

  // Created Mar 28 2013 by RK
  // load PSNID_BEST_RESULTS.FINAL_PARVAL and FINAL_PARERR so that 
  // the output arrays have a well-define address.
  // The final params are either from GRID-search or from MCMC.
  // And errors are filled only for Ia using MCMC.
  //
  // Jul 2 2013 RK - 
  // for NONIA_INDEX (shape par), translate internal index to absolute
  //
  // Aug 22 2013: fix aweful bug setting SHAPEPAR for NONIA
  //
  // Nov 1 2013 RK - store PEAKMAGS

  // --------------- BEGIN -------------

  int t, z, ipar, ISIA, ifilt_obs, USE, IFILT ;
  int ONE = 1 ;
  double dif, PARVAL, TOBS, MAG, MAGERR ;

  z=0;
  for ( t=0; t < PSNID_NTYPES; t++ ) {  // loop over itypes

    ISIA = ( t == PSNID_ITYPE_SNIA ) ;

    for ( ipar=0; ipar < PSNID_NPARAM; ipar++ ) {


    if ( MCMC_RUN && ISIA ) {
      PARVAL = PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][t][ipar][2] ;
      PSNID_BEST_RESULTS.FINAL_PARVAL[z][t][ipar] = PARVAL ;
	
      dif = 	
	PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][t][ipar][3] -
	PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][t][ipar][1] ;
      PSNID_BEST_RESULTS.FINAL_PARERR[z][t][ipar] = dif/2.0 ;
    }
    else {
      // grid search
      PARVAL = PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][t][ipar] ;
      PSNID_BEST_RESULTS.FINAL_PARVAL[z][t][ipar] = PARVAL ;	
      PSNID_BEST_RESULTS.FINAL_PARERR[z][t][ipar] = -999. ; // no error
    }

    // for NONIA_INDEX (shape par), translate internal index to absolute
    if ( ipar == PSNID_PARAM_SHAPEPAR && ISIA == 0 ) {
      int ISHAPE, INDX_SPARSE, INDX_SIM;
      ISHAPE      = (int)PARVAL ;

      // store SIM-INDX instead of local internal index so that
      // the FINAL NONIA index (in plots and fitres file) is tracable
      INDX_SPARSE = PSNID_NONIA_ABSINDEX[t][ISHAPE]; 
      INDX_SIM    = SNGRID_PSNID[2].NON1A_INDEX[INDX_SPARSE] ;
      PSNID_BEST_RESULTS.FINAL_PARVAL[z][t][ipar] = (double)INDX_SIM ;

      // store internal sparse index
      PSNID_BEST_RESULTS.FINAL_NONIA_INDX_SPARSE[t] = ISHAPE ;
    }

    } // end loop over ipar

    // now the peakmags.
    TOBS = (double)0.0 ; // peak
    for(ifilt_obs=0; ifilt_obs < MXFILTINDX; ifilt_obs++ ) {
      PSNID_BEST_RESULTS.FINAL_PEAKMAG[t][ifilt_obs] = 999.0 ;
      USE = PSNID_INPUTS.USEFILT_PEAKMAG_STORE[ifilt_obs] ;
      if ( USE == 0 ) { continue ; }
      
      IFILT = PSNID_INPUTS.IFILTLIST_INV[ifilt_obs] ; // sparse filt index
      PSNID_BEST_GET_FITFUN(PSNID_BEST_RESULTS.CID, t,
			  IFILT, ONE,          // (I) ifilt, Nep
			  &TOBS       ,        // (I) Tobs
			  &MAG, &MAGERR ) ;    // (O) best-fit mag and error
      PSNID_BEST_RESULTS.FINAL_PEAKMAG[t][ifilt_obs] = MAG ;
    }

  } // end loop over types


} // end of psnid_best_store_finalPar




/************************************************************************/
/************************************************************************/
/**************                                          ****************/
/**************                MS Functions              ****************/
/**************                                          ****************/
/************************************************************************/
/************************************************************************/


/**********************************************************************/
void psnid_best_split_nonia_types(int *types, int optdebug)
/**********************************************************************/
/***

    Sets types[0..PSNID_MAXNL_NONIA] to:
          = 1 for Ibc
          = 2 for II

    Set optdebug to print out diagnostics.


  Jan 27, 2013 RK
   replace hard-wired string comparisons with 
   logical function psnid_strTypeMatch(ITYPE,CTYPE,0).
   List of string-types are now declared in the header above
   in char PSNID_CHARTYPE_LISTS[][].

  Feb 13 2017 RK
   Check for PEC1A. Still using fragile logic.

  Aug 28 2017: check for MODEL1, MODEL2, etc ...

  May 2022 RK - abort if any type is not defined; see NERR
 ***/
  
/**********************************************************************/
{
  int i, nibc=0, nii=0, npec1a=0, nignore=0, ITYPE, NERR=0 ;
  int nmodel[PSNID_NTYPES], nmodel_sum=0 ;
  char *CTYPE ;
  char fnam[] = "psnid_best_split_nonia_types";

  // ------------ BEGIN -----------------

  for(i=0; i<PSNID_NTYPES; i++ ) { nmodel[i]=0; }
  
  for (i=1; i <= PSNID_MAXNL_NONIA; i++) {

    types[i-1] = -9 ;

    // define local vars (RK Jan 2013)
    ITYPE = SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NON1A_ITYPE_AUTO[i] ;
    CTYPE = SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NON1A_CTYPE[i] ;
    
    if ( optdebug )
      { printf(" xxx --- '%s'(%d) ----- \n", CTYPE,ITYPE); }
    
    if ( ITYPE > 0 ) {

      // ------------------------
      // Note to Masao from RK: 
      // would be more robust to loop over types here
      // in case more types are added later.
      // ------------------------

      if ( psnid_strTypeMatch(PSNID_ITYPE_SNIBC,CTYPE,0) ) {
	// IBC
	types[i-1] = PSNID_ITYPE_SNIBC;
	nibc = nibc + 1;
	PSNID_NONIA_ABSINDEX[PSNID_ITYPE_SNIBC][nibc] = i;

	if (optdebug == 1) {
	  printf("\t xxx index = %2d is a Ibc   type=%d     itype=%d\n",
		 i, types[i-1], ITYPE);
	}

      } 
      else if ( psnid_strTypeMatch(PSNID_ITYPE_SNII,CTYPE,0) ) {
	// II
	types[i-1] = PSNID_ITYPE_SNII;
	nii = nii + 1;
	PSNID_NONIA_ABSINDEX[PSNID_ITYPE_SNII][nii] = i;

	if (optdebug == 1) {
	  printf("\t xxx index = %2d is a II    type=%d     itype=%d\n",
		 i, types[i-1], ITYPE);
	}
      }
      else if ( psnid_strTypeMatch(PSNID_ITYPE_PEC1A,CTYPE,0) ) {
	// PEC1A (added Feb 2017 by RK)  fragile
	types[i-1] = PSNID_ITYPE_PEC1A;
	npec1a = npec1a + 1;
	PSNID_NONIA_ABSINDEX[PSNID_ITYPE_PEC1A][npec1a] = i;

	if (optdebug == 1) {
	  printf("\t xxx index = %2d is a PEC1A  type=%d     itype=%d\n",
		 i, types[i-1], ITYPE);
	}      
      }
      else if ( psnid_strTypeMatch(PSNID_ITYPE_MODEL1,CTYPE,0) ) {
	// MODEL1 (added Aug 2017 by RK)  fragile
	types[i-1] = PSNID_ITYPE_MODEL1;
	nmodel[1]++ ;   nmodel_sum++ ;
	PSNID_NONIA_ABSINDEX[PSNID_ITYPE_MODEL1][nmodel[1]] = i;	
      }
      else if ( psnid_strTypeMatch(PSNID_ITYPE_MODEL2,CTYPE,0) ) {
	// MODEL2 (added Aug 2017 by RK)  fragile
	types[i-1] = PSNID_ITYPE_MODEL2;
	nmodel[2]++ ; nmodel_sum++ ;
	PSNID_NONIA_ABSINDEX[PSNID_ITYPE_MODEL2][nmodel[2]] = i;	
      }
      else if ( psnid_strTypeMatch(PSNID_ITYPE_MODEL3,CTYPE,0) ) {
	// MODEL3 (added Aug 2017 by RK)  fragile
	types[i-1] = PSNID_ITYPE_MODEL3;
	nmodel[3]++ ; nmodel_sum++ ;
	PSNID_NONIA_ABSINDEX[PSNID_ITYPE_MODEL3][nmodel[3]] = i;	
      }
      else if ( psnid_strTypeMatch(PSNID_ITYPE_MODEL4,CTYPE,0) ) {
	// MODEL4 (added Aug 2017 by RK)  fragile
	types[i-1] = PSNID_ITYPE_MODEL4;
	nmodel[4]++ ;  nmodel_sum++ ;
	PSNID_NONIA_ABSINDEX[PSNID_ITYPE_MODEL4][nmodel[4]] = i;	
      }
    }
    
    else {  // ITYPE <= 0
      types[i-1] = ITYPE;
      nignore = nignore + 1;
      if (optdebug == 1) {
	printf("\t xxx index = %2d is ignored            itype=%d\n",
	       i, ITYPE);  }
    }

    if ( types[i-1] < 0 )  {
      printf(" ERROR: unknown type for '%s'(%d) \n", CTYPE,ITYPE);
      NERR++ ;
    }

  }

  if (NERR > 0 || optdebug == 1) {
    printf("\t There are %2d Ibc   templates\n", nibc);
    printf("\t There are %2d II    templates\n", nii);
    printf("\t There are %2d PEC1A templates\n", npec1a);

    for(i=1; i <= 4; i++ ) {
      printf("\t There are %2d MODEL%d templates\n", nmodel[i],i);
    }
    
    printf("\t There are %2d  templates ignored\n\n", nignore);
  }

  if ( NERR > 0 ) {
    sprintf(c1err,"%d templates with unknown types", NERR);
    sprintf(c2err,"Check ERROR messages above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }

  PSNID_NGRID[PSNID_ITYPE_SNIA]  =
    USEFLAG_TEMPLATES_PSNID[TYPEINDX_SNIA_PSNID];
  PSNID_NGRID[PSNID_ITYPE_SNIBC] = nibc ;
  PSNID_NGRID[PSNID_ITYPE_SNII]  = nii ;
  PSNID_NGRID[PSNID_ITYPE_PEC1A] = npec1a ;
  for(i=1; i<=4; i++) { PSNID_NGRID[PSNID_ITYPE_MODEL1+i-1]=nmodel[i]; }
  
  if (nibc + nii + + npec1a + nmodel_sum + nignore != PSNID_MAXNL_NONIA) {
    // Feb 28 2013  Q1: why not abort here ?
    printf("\n");
    printf(" %s WARNING: number of non-Ia types don't add up!\n", fnam);
    printf("\t nibc=%d nii=%d nignore=%d  PSNID_MAXNL_NONIA=%d \n",
	   nibc, nii, nignore, PSNID_MAXNL_NONIA) ;
    fflush(stdout);
  }

  return;
}
// end of psnid_best_split_nonia_types


/***************************************************************************/
int psnid_strTypeMatch(int itype, char *string, int optDebug ) {

  // Ceated Jan 27, 2013 by R.Kessler
  // Returns 1(TRUE) of input *string matches integer itype.
  // Returns 0(FALSE) if not.
  // optdebug =1 => print result to screen for debug.
  
  
  int  ISTAT ;
  char cList[80], *ptrtok, ctmp[20] ;
  char fnam[] = "psnid_strTypeMatch" ;

  // --------------- BEGIN -------------

  ISTAT = 0;  // default = FALSE
  sprintf(cList,"%s", PSNID_CHARTYPE_LISTS[itype] ) ;
  
  // check each sub-string separated by blank space
  ptrtok = strtok(cList," ");

  while ( ptrtok != NULL ) {
    sprintf(ctmp,"%s", ptrtok);
    if ( strcmp(string,ctmp) == 0 ) { ISTAT=1; goto DONE ; }
    ptrtok = strtok(NULL," " );

  } // ptrtok

 DONE:

  if ( optDebug ) {
    printf("\n");
    printf(" XXX ---------------------------------------------- \n");
    printf(" XXX %s DUMP \n", fnam);
    printf(" XXX itype = %d  string = '%s' \n", itype, string);
    printf(" XXX ==>  cList='%s'  ISTAT=%d \n", 
	   PSNID_CHARTYPE_LISTS[itype], ISTAT);
    printf("\n");
    fflush(stdout);
  }  

  return ISTAT ;

} // end of psnid_strTypeMatch

/***************************************************************************/
void psnid_best_set_grid_limits(int typeindex)
/*
  Sets following global variables:

  PSNID_THIS_TYPE   = 1 (Ia) or 2 (non-Ia)
  PSNID_NFILTER
  PSNID_MAXND

  PSNID_MAXNZ
  PSNID_ZMIN
  PSNID_ZSTEP

  PSNID_MAXNL
  PSNID_LMIN
  PSNID_LSTEP

  Feb 13 2017: check for PEC1A, but still using fragile logic.
  Aug 28 2017: include MODEL1[234]

 */
/******************************************************************************/
{
  int this_type=0, tmp_maxnl ;
  char fnam[] = "psnid_best_set_grid_limits" ;

  // ---------------- BEGIN ---------------

  // grid models should distinguish SN II from SN Ibc!
  if (typeindex == PSNID_ITYPE_SNIA) {
    this_type       = TYPEINDX_SNIA_PSNID;
    PSNID_THIS_TYPE = TYPEINDX_SNIA_PSNID;;
  }
  else if (typeindex >  PSNID_ITYPE_SNIA && typeindex <= PSNID_NTYPES ) {
    this_type       = TYPEINDX_NONIA_PSNID ;
    PSNID_THIS_TYPE = TYPEINDX_NONIA_PSNID ;
  }
  else {
    sprintf(c1err,"Unknown typeindex = %d", typeindex );
    sprintf(c2err,"Check grid-read summary.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }
  

  PSNID_NFILTER = SNGRID_PSNID[this_type].NBIN[IPAR_GRIDGEN_FILTER];
  PSNID_MAXND   = SNGRID_PSNID[this_type].NBIN[IPAR_GRIDGEN_TREST];
  PSNID_TBIN    = SNGRID_PSNID[this_type].BINSIZE[IPAR_GRIDGEN_TREST];
  
  // redshift bins
  PSNID_MAXNZ  = SNGRID_PSNID[this_type].NBIN[IPAR_GRIDGEN_LOGZ];
  PSNID_ZMIN   = pow(10,SNGRID_PSNID[this_type].VALMIN[IPAR_GRIDGEN_LOGZ]);
  PSNID_ZSTEP  = SNGRID_PSNID[this_type].BINSIZE[IPAR_GRIDGEN_LOGZ];

  // shape bins
  tmp_maxnl     = SNGRID_PSNID[this_type].NBIN[IPAR_GRIDGEN_SHAPEPAR];
  if (typeindex == PSNID_ITYPE_SNIA) {
    PSNID_MAXNL = tmp_maxnl;
  }
  else {    
    PSNID_MAXNL = PSNID_NGRID[typeindex] ;
  }

  
  PSNID_LMIN    = SNGRID_PSNID[this_type].VALMIN[IPAR_GRIDGEN_SHAPEPAR];
  PSNID_LSTEP   = SNGRID_PSNID[this_type].BINSIZE[IPAR_GRIDGEN_SHAPEPAR];
  //  printf("typeindex = %d  this_type = %d  PSNID_MAXNL = %d\n",
  //	 typeindex, this_type, PSNID_MAXNL);


  PSNID_PARAM_MAX_INDEX[PSNID_PARAM_LOGZ]     = PSNID_MAXNZ;
  PSNID_PARAM_MAX_INDEX[PSNID_PARAM_SHAPEPAR] = PSNID_MAXNL;
  PSNID_PARAM_MAX_INDEX[PSNID_PARAM_COLORPAR] = PSNID_MAXNA;
  PSNID_PARAM_MAX_INDEX[PSNID_PARAM_COLORLAW] = 1;
  PSNID_PARAM_MAX_INDEX[PSNID_PARAM_TMAX]     = PSNID_MAXND;
  PSNID_PARAM_MAX_INDEX[PSNID_PARAM_DMU]      = PSNID_MAXNU;


  // Aug 29 2017: if Trest bin size is smaller than last-iter bin-size,
  //              reset last-iter bin-size
  PSNID_TBIN = SNGRID_PSNID[this_type].BINSIZE[IPAR_GRIDGEN_TREST];
  if ( PSNID_TBIN < PSNID_TSTEP[2] )
    { PSNID_TSTEP[2]  =  PSNID_TBIN ; }
  
  return;
}
// end of psnid_best_set_grid_limits


/******************************************************************************/
void psnid_best_set_grid_values(int typeindex, int *nonia_types)
/******************************************************************************/
{

  // Sep  7 2013 RK - use wrapper "psnid_best_modelErr" to get model error
  // Oct 23 2013 RK - set redshift = 10**logz outside epoch loop ; 
  //        see logz variable.
  //
  // Oct 10 2014 RK - allocate memory for test1,mag1,magerr1 and for 1->2
  // Feb 13 2017 RK - check for PEC1A

  int i,j,k,m, mfilt, this_type=0, type_count=0;
  int ind1c1, ind2c1;
  int tmp_maxnl=-9 , this_filt;
  double trest1[MXEP_PSNID], mag1[MXEP_PSNID], magerr1[MXEP_PSNID];
  double trest2[MXEP_PSNID], mag2[MXEP_PSNID], magerr2[MXEP_PSNID];
  double color1, color2, redshift, logz, dmdc;
  int nepoch1, nepoch2, optDump=0;
  char fnam[] = "psnid_best_set_grid_values" ;

  // ------------ BEGIN -------------------

  // grid models should distinguish SN II from SN Ibc!
  if (typeindex == PSNID_ITYPE_SNIA) {
    this_type = TYPEINDX_SNIA_PSNID;
    tmp_maxnl = PSNID_MAXNL;
  } 
  else if ( typeindex > PSNID_ITYPE_SNIA && typeindex <= PSNID_NTYPES ) {
    this_type = TYPEINDX_NONIA_PSNID;
    tmp_maxnl = PSNID_MAXNL_NONIA;
  }
  else {   
    sprintf(c1err,"Invalid typeindex = %d", typeindex);
    sprintf(c2err,"         ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }


  // get models of first and last color bins
  // mags from first color bin
  ind1c1 = 1;
  ind2c1 = SNGRID_PSNID[this_type].NBIN[IPAR_GRIDGEN_COLORPAR];

  color1 = SNGRID_PSNID[this_type].VALUE[IPAR_GRIDGEN_COLORPAR][ind1c1];
  color2 = SNGRID_PSNID[this_type].VALUE[IPAR_GRIDGEN_COLORPAR][ind2c1];

  // color of base mags (global)
  PSNID_BASE_COLOR[typeindex] = color1;
  //  printf("  color1=%8.3f\n", color1);
  //  printf("  color2=%8.3f\n", color2);


  // shape, z, epoch, filter for 0 color
  for (i=1; i<=tmp_maxnl; i++) {                // shape/lumin

    if ( this_type == TYPEINDX_SNIA_PSNID ||
	(this_type == TYPEINDX_NONIA_PSNID && nonia_types[i-1] == typeindex)) {

      type_count++;
      //      printf("\t %d %d %d %d %d\n",
      //	     i,this_type,nonia_types[i], typeindex, type_count);

      for (j=1; j<=PSNID_MAXNZ; j++) {          // redshift

	logz     = SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_LOGZ][j] ;
	redshift = pow(10., logz);
	    
	for (m=1; m<=PSNID_NFILTER; m++) {      // internal filter index
	  mfilt = m-1 ;

	  // absolute filter index
	  this_filt = PSNID_INPUTS.IFILTLIST[m-1];
	  if (PSNID_INPUTS.USEFILT[this_filt] == 1) {  // use this filter (from nml)
	    // Get first and last color bins
  	    //                         z      c1 c2 lu  f
	    get_template_lc(this_type, j, ind1c1, 1, i, this_filt, optDump,
			    &nepoch1, trest1, mag1, magerr1);
	    get_template_lc(this_type, j, ind2c1, 1, i, this_filt, optDump,
			    &nepoch2, trest2, mag2, magerr2);

	    for (k=0; k<PSNID_MAXND; k++) {       // epoch

	      PSNID_MODEL_EPOCH[type_count][j][k+1]   
		= trest1[k]*(1.+redshift);

	      PSNID_MODEL_MAG[m][type_count][j][k+1]  
		= mag1[k]  + psnid_best_mwxtmag(mfilt);

	      // RK Sep 7 2013 - call wrapper for model mag-error
	      PSNID_MODEL_MAGERR[m][type_count][j][k+1]  = 
		psnid_best_modelErr( this_type, trest1[k], magerr1[k] );  

	      dmdc = (mag1[k]-mag2[k])/(color2-color1) ;
	      PSNID_MODEL_EXTINCT[m][type_count][j][k+1]    = dmdc ;
	      PSNID_MODEL_MWEXTINCT[m][type_count][j][k+1]  = dmdc ;

	    }
	  }
	}
      }

    }
  }

  return;
}
// end of psnid_best_set_grid_values

// ================================
double psnid_best_modelErr(int TYPEINDX, double Trest, double modelErr_grid){

  // Created Sep 9 2013 by R.Kessler
  // return mag model-error based on user option xxx

  // Inputs:
  //  - TYPEINDX = IA or CC
  //  - Trest    = rest-frame epoch relative to peak brightness (days)
  //  - modelErr_grid = model error read from grid (for GRID option)

  int    ISIA ;
  double MAGERR=0.0 ;
  char   *modelName ;
  char   fnam[] = "psnid_best_modelErr" ;

  // --------------- BEGIN -----------------

  // strip off string-name of mag-error model into local pointer
  modelName = PSNID_INPUTS.MODELNAME_MAGERR ; 
  
  ISIA = (TYPEINDX == TYPEINDX_SNIA_PSNID) ;

  if ( strcmp(modelName,"GRID") == 0 ) 
    {  MAGERR = modelErr_grid ; }

  else if ( strcmp(modelName,"S14") == 0 ) 
    {  MAGERR = psnid_best_modelErr_S14(TYPEINDX,Trest) ; }

  else if ( strcmp(modelName,"S13") == 0 ) 
    {  MAGERR = psnid_best_modelErr_S13(TYPEINDX,Trest) ; }

  else if ( strcmp(modelName,"S11") == 0 ) 
    {  MAGERR = psnid_best_modelErr_S11(TYPEINDX,Trest) ; }

  else if ( strcmp(modelName,"S13_H20") == 0 ) {
    if ( ISIA ) {
      MAGERR = psnid_best_modelErr_S13(TYPEINDX,Trest);
    } else {
      MAGERR = modelErr_grid;
    }
  }

  else 
    {
      sprintf(c1err,"Invalid  MODELNAME_MAGERR = '%s' ", modelName );
      sprintf(c2err,"Check &PSNIDINP namelist.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }


  // optional fudge on template errors
  ISIA = (TYPEINDX == TYPEINDX_SNIA_PSNID) ;
  if ( ISIA ) {
    MAGERR *= PSNID_INPUTS.TEMPLERR_SCALE_SNIA ;
  } else {
    MAGERR *= PSNID_INPUTS.TEMPLERR_SCALE_NONIA ;
  }

  return  MAGERR ;

} // end of psnid_best_modelErr


// ==============================================
double psnid_best_modelErr_S14(int TYPEINDX, double Trest) {

  // Feb 17 2014 RK
  // exploring new model errors
  // model mag-err for this TYPEINDEX and Trest.

  int    ISIA ;
  double MAGERR, absT ;

  MAGERR = 0.0 ;
  absT   = fabs(Trest);

  ISIA = (TYPEINDX == TYPEINDX_SNIA_PSNID) ;

  if ( ISIA ) {
    if (Trest < -10.0) 
      {  MAGERR  =  0.08*absT/10.0;  }

    else if ( absT < 20.0 ) 
      {  MAGERR  =  0.06 + 0.04*absT/20.0;  } 

    else if ( absT >= 20.0) 
      {  MAGERR =  0.10 + 0.18*(absT-20.0)/60.0;  } 
  }
  else 
    {  MAGERR =  0.06 + 0.08*absT/60.0;  } // CC 


  //  MAGERR = 0.06;


  return MAGERR ;

} // end of psnid_best_modelErr_S14


// ==============================================
double psnid_best_modelErr_S13(int TYPEINDX, double Trest) {

  // Sep 7 2013 RK 
  // wrapper to return S13 (Sako et al 2013 SDSS data release) 
  // model mag-err for this TYPEINDEX and Trest.

  int    ISIA ;
  double MAGERR, absT ;

  MAGERR = 0.0 ;
  absT   = fabs(Trest);

  ISIA = (TYPEINDX == TYPEINDX_SNIA_PSNID) ;

  if ( ISIA ) {
    if (Trest < -10.0) 
      {  MAGERR  = 0.3; } 

    else if ( absT < 20.0 ) 
      {  MAGERR  =  0.06 + 0.04*absT/20.0;  } 

    else if ( absT >= 20.0) 
      {  MAGERR =  0.10 + 0.18*(absT-20.0)/60.0;  } 
  }
  else 
    {  MAGERR =  0.06 + 0.08*absT/60.0;  } // CC 



  return MAGERR ;

} // end of psnid_best_modelErr_S13


// ==============================================
double psnid_best_modelErr_S11(int TYPEINDX, double Trest) {

  // Sep 7 2013 RK 
  // wrapper to return S11 (Sako et al 2011; Eqs 6-7)
  // model mag-err for this TYPEINDEX and Trest.

  int    ISIA ;
  double MAGERR, absT ;

  MAGERR = 0.0 ;
  absT   = fabs(Trest);

  ISIA = (TYPEINDX == TYPEINDX_SNIA_PSNID) ;

  if ( ISIA ) {


    if ( absT < 20.0 ) 
      {  MAGERR  =  0.08 + 0.04*absT/20.0;  }    // Eq. 6
    else 
      {  MAGERR =  0.12 + 0.08*(absT-20.0)/60.0;  }   // Eq. 6
  }
  else 
    {  MAGERR =  0.08 + 0.08*absT/60.0;  } // CC : Eq. 7


  return MAGERR ;

} // end of psnid_best_modelErr_S11


/******************************************************************************/
double psnid_best_mwxtmag(int ifilt)
/******************************************************************************/
{
  // Apr 8 2013
  // return mag of Galactic extinction for ifilt = 0 - N-1

  double AV, MWEBV, XTMW ;
  MWEBV = DATA_PSNID_DOFIT.MWEBV ;
  AV    = RV_MWDUST * MWEBV ;
  XTMW  = AV * PSNID_INPUTS.XTMW_at_AV1[ifilt] ;
  return  XTMW ;
}
// end of psnid_best_mwxtmag


/******************************************************************************/
void psnid_best_dump_model_lc(int shape_index, int z_index)
/******************************************************************************/
{

  int i,j,k,m;
  double mJy, mJyerr, fluxcal, fluxcalerr;

  i = shape_index;       // shape/luminosity bin
  j = z_index;           // redshift bin

  //
  printf("\n XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
  printf(" XXXXXXXXXX  BEGIN: light curve model dump  XXXXXXXXXX\n");
  printf("\t shape_idnex = %d  z_index = %d\n", i, j);

  for (k=1; k<=PSNID_MAXND; k++) {        // epoch
    printf("%2d  %12.5f  ", k, PSNID_MODEL_EPOCH[i][j][k]);
    for (m=1; m<=PSNID_NFILTER; m++) {    // filter

      psnid_pogson2mJy(PSNID_MODEL_MAG[m][i][j][k],
		       PSNID_MODEL_MAGERR[m][i][j][k], &mJy, &mJyerr);
      psnid_pogson2fluxcal(PSNID_MODEL_MAG[m][i][j][k],
			   PSNID_MODEL_MAGERR[m][i][j][k], &fluxcal, &fluxcalerr);

      printf("   %8.3f %7.3f",
	     PSNID_MODEL_MAG[m][i][j][k],
	     PSNID_MODEL_MAGERR[m][i][j][k]);
      //      printf("   %8.3f %7.3f", PSNID_MODEL_EXTINCT[m][i][j][k],
      //	     PSNID_MODEL_MWEXTINCT[m][i][j][k]);
      //      printf("   %9.3f %9.3f", mJy, mJyerr);
      //      printf("   %10.2f %10.2f", fluxcal, fluxcalerr);
    }
    printf("\n");
  }
  printf(" XXXXXXXXXX   END: light curve model dump  XXXXXXXXXX\n");
  printf(" XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");

  return;
}
// end of psnid_best_dump_model_lc





/**********************************************************************/
/**********************************************************************/
/***********************                        ***********************/
/***********************  CORE PSNID FUNCTIONS  ***********************/
/***********************                        ***********************/
/**********************************************************************/
/**********************************************************************/




/**********************************************************************/
void psnid_best_model_alloc()
/**********************************************************************/
{
  
  //                            shape        redshift          epoch
  PSNID_MODEL_EPOCH     = d3tensor(ONE8,PSNID_MAXNL, ONE8,PSNID_MAXNZ, 
				   ONE8,PSNID_MAXND);
  //                                  filter     shape       redshift          epoch
  PSNID_MODEL_MAG       = d4tensor(ONE8,PSNID_NFILTER, ONE8,PSNID_MAXNL, 
				   ONE8,PSNID_MAXNZ, ONE8,PSNID_MAXND);

  PSNID_MODEL_MAGERR    = d4tensor(ONE8,PSNID_NFILTER, ONE8,PSNID_MAXNL, 
				   ONE8,PSNID_MAXNZ, ONE8,PSNID_MAXND);

  PSNID_MODEL_EXTINCT   = d4tensor(ONE8,PSNID_NFILTER, ONE8,PSNID_MAXNL, 
				   ONE8,PSNID_MAXNZ, ONE8,PSNID_MAXND);

  PSNID_MODEL_MWEXTINCT = d4tensor(ONE8,PSNID_NFILTER, ONE8,PSNID_MAXNL, 
				   ONE8,PSNID_MAXNZ, ONE8,PSNID_MAXND);


  return;
}
// end of psnid_best_model_alloc


/**********************************************************************/
void psnid_best_model_free()
/**********************************************************************/
{

  free_d3tensor(PSNID_MODEL_EPOCH,  ONE8, PSNID_MAXNL, 
		ONE8,PSNID_MAXNZ, ONE8,PSNID_MAXND);

  free_d4tensor(PSNID_MODEL_MAG, ONE8,PSNID_NFILTER, ONE8,
		PSNID_MAXNL, ONE8,PSNID_MAXNZ, ONE8,PSNID_MAXND);

  free_d4tensor(PSNID_MODEL_MAGERR,  ONE8,PSNID_NFILTER, ONE8,
		PSNID_MAXNL, ONE8, PSNID_MAXNZ, ONE8, PSNID_MAXND);

  free_d4tensor(PSNID_MODEL_EXTINCT, ONE8, PSNID_NFILTER, ONE8,
		PSNID_MAXNL, ONE8, PSNID_MAXNZ, ONE8, PSNID_MAXND);

  free_d4tensor(PSNID_MODEL_MWEXTINCT, ONE8,PSNID_NFILTER, ONE8,
		PSNID_MAXNL, ONE8, PSNID_MAXNZ, ONE8,PSNID_MAXND);


  return;
}
// end of psnid_best_model_free


/**********************************************************************/
void psnid_best_setup_searchgrid()
/**********************************************************************/
{

  int N ;
  double DIF;
  //  char fnam[] = "psnid_best_setup_searchgrid" ;
  
  // --------------- BEGIN ---------------
  printf("\t Setup search grid. \n"); fflush(stdout);

  PSNID_ZRBN[0] =  6;
  PSNID_ZRBN[1] =  3;
  PSNID_ZRBN[2] =  1;
  PSNID_ZWID[0] = 10;
  PSNID_ZWID[1] = 15;
  PSNID_ZWID[2] = 30;

  // color/extinction

  N   = PSNID_INPUTS.NCOLOR ;
  DIF = (PSNID_INPUTS.COLOR_MAX - PSNID_INPUTS.COLOR_MIN) ;
  PSNID_AMIN  = PSNID_INPUTS.COLOR_MIN;
  PSNID_MAXNA = PSNID_INPUTS.NCOLOR;
  if ( N > 0 ) 
    { PSNID_ASTEP = (DIF)/(double)(N-1.0); }
  else
    { PSNID_ASTEP = 0.0 ; }


  PSNID_ARBN[0] =  4;
  PSNID_ARBN[1] =  2;
  PSNID_ARBN[2] =  1;
  PSNID_AWID[0] =  3;
  PSNID_AWID[1] =  5;
  PSNID_AWID[2] = 10;

  // standalone PSNID
  PSNID_ARBN[0] =  6;
  PSNID_ARBN[1] =  3;
  PSNID_ARBN[2] =  1;
  PSNID_AWID[0] =  2;
  PSNID_AWID[1] =  2;
  PSNID_AWID[2] =  5;
  /*
    PSNID_ARBN[0] = 10;
    PSNID_ARBN[1] =  6;
    PSNID_ARBN[2] =  1;
    PSNID_AWID[0] =  6;
    PSNID_AWID[1] =  4;
    PSNID_AWID[2] = 10;
  */


  // mu
  PSNID_URBN[0] =  5;
  PSNID_URBN[1] =  3;
  PSNID_URBN[2] =  1;
  PSNID_UWID[0] =  5;
  PSNID_UWID[1] =  5;
  PSNID_UWID[2] = 10;


  // Tmax (passed from namelist, Sep 15 2017)
  int iter; 
  for(iter = 0 ; iter < MXITER_PSNID; iter++ ) {
    PSNID_TSTART[iter] = PSNID_INPUTS.TMAX_START[iter]; 
    PSNID_TSTOP[iter]  = PSNID_INPUTS.TMAX_STOP[iter]; 
    PSNID_TSTEP[iter]  = PSNID_INPUTS.TMAX_STEP[iter]; 
  }


  return;

}
// end of psnid_best_setup_searchgrid


/**********************************************************************/
void psnid_best_grid_compare(int itype, int zpind, int nobs,
			     int *data_filt, double *data_mjd,
			     double *data_fluxcal, double *data_fluxcalerr,
			     int *useobs,
			     int ***ind, double **evidence)
/***
  input:

    int itype                 = 0, 1, 2 for Ia, Ibc, II

    (These are also global variables NOBS, IFILTOBS, MJD, FLUX, FLUXERR)
    int nobs                  = number of filter-epochs
    int *data_filt            = ivector(0,NOBS);
    double *data_mjd          = dvector(0,NOBS);
    double *data_fluxcal      = dvector(0,NOBS);
    double *data_fluxcalerr   = dvector(0,NOBS);

    ind zpind                 = index for redshift priors

  global:

    PSNID_BEST_RESULTS.ZPRIOR           = redshift prior
    PSNID_BEST_RESULTS.ZPRIOR_ERR       = redshift prior uncertainty

    double ***PSNID_MODEL_EPOCH    = d3tensor(1,MAXNL, 1,MAXNZ, 1,MAXND);
    double ****PSNID_MODEL_MAG     = d4tensor(1,NFILTER, 1,MAXNM, 1,MAXNZ, 1,MAXND);
    double ****PSNID_MODEL_MAGERR  = d4tensor(1,NFILTER, 1,MAXNM, 1,MAXNZ, 1,MAXND);
    double ****PSNID_MODEL_EXTINCT = d4tensor(1,NFILTER, 1,MAXNM, 1,MAXNZ, 1,MAXND);


  output;

    int ***ind         = i3tensor(0,NITER+1, 0,NTYPES ,0,PSNID_NPARAM);
                       = indices of best-fit model (see PSNID_PARAM_XXXXX)
    double **evidence  = dmatrix(0,PSNID_NZPRIOR, 0,PSNID_NTYPES);
                       = Bayesian evidence for each z prior and type

***/
/**********************************************************************/
{
  int i, ipass, ngood=0, thisngood=1, z, d, a, f, t, u, jtmp, this_z=1;
  int minz=1, maxz=1, zstep=1, mind=1, maxd=1, dstep=1,
    mina=1, maxa=1, astep=1, mini=1, maxi=1, minu=1, maxu=1, ustep=1;
  int npeak;
  double peak_guess, chisq, chisqlo, chisqloz, chisqlozu, chisq_z, wgt;
  double istep = 1.0;
  double *c_grid, *z_grid, *u_grid, ushift;
  double *fit_epoch, **fit_mag, **fit_magerr;
  double  pav=0.0;
  double ZPRIOR, ZPRIOR_ERR, DZ, ZSIG ;
  int    DOPRIOR_ZSPEC, DOPRIOR_ZPHOT, optDebug, NON1A_INDEX, isp;
  int    AWID, ZWID, ZRBN, ARBN, UWID, URBN, indTmp ;

  //  char fnam[] = "psnid_best_grid_compare" ;

  // --------------------- BEGIN ------------------

  chisqlo = PSNID_BIGN;

  fit_epoch  = dvector(1,PSNID_MAXND);
  fit_mag    = dmatrix(1,PSNID_NFILTER, 1,PSNID_MAXND);
  fit_magerr = dmatrix(1,PSNID_NFILTER, 1,PSNID_MAXND);


  c_grid   = dvector(1,PSNID_MAXNA);
  z_grid   = dvector(1,PSNID_MAXNZ);
  u_grid   = dvector(1,PSNID_MAXNU);
  psnid_best_get_c_grid(c_grid);
  psnid_best_get_z_grid(z_grid);
  psnid_best_get_u_grid(u_grid);


  hunt(z_grid, PSNID_MAXNZ, PSNID_BEST_RESULTS.ZPRIOR[zpind], &this_z);
  // z_grid[this_z] <= PSNID_BEST_RESULTS.ZPRIOR[zpind] < z_grid[this_z+1]
  //  if (PSNID_BEST_RESULTS.ZPRIOR[zpind] - z_grid[this_z] >
  //      z_grid[this_z+1] - PSNID_BEST_RESULTS.ZPRIOR[zpind]) this_z += 1;
  if (this_z < 1)           { this_z = 1; }
  if (this_z > PSNID_MAXNZ) { this_z = PSNID_MAXNZ; }


  // local variables and logical for ZPRIOR (RK)
  ZPRIOR      = PSNID_BEST_RESULTS.ZPRIOR[zpind] ;
  ZPRIOR_ERR  = PSNID_BEST_RESULTS.ZPRIOR_ERR[zpind] ;
  DOPRIOR_ZSPEC = DOPRIOR_ZPHOT = 0;
  if ( ZPRIOR > 0.0 ) {
    ZSIG           = ZPRIOR/ZPRIOR_ERR ;
    DOPRIOR_ZSPEC  = ZSIG > PSNID_SPECZSIG ;
    DOPRIOR_ZPHOT  = ZSIG < PSNID_SPECZSIG ;
  }

  int AVOID_SIMCHEAT = 
    ( PSNID_INPUTS.OPT_SIMCHEAT==0 && PSNID_INPUTS.LSIM  && itype > 0 ) ;

  /* xxxxxxxxxxx
  printf(" xxx %s: AVOID=%d LSIM=%d  itype=%d  SIM_NON1A_INDEX=%d\n",
	 fnam, AVOID_SIMCHEAT, PSNID_INPUTS.LSIM, itype,
	 DATA_PSNID_DOFIT.SIM_NON1A_INDEX ) ; fflush(stdout);
	 xxxxxxxxxxx */


  // do all shape/lumin, CC templates each pass
  mind  = 1;
  maxd  = PSNID_MAXNL;
  dstep = 1;

  int NITER = PSNID_NITER ;

  for (ipass = 0; ipass <= NITER; ipass++) { 

    optDebug = 0 ;  // RK
    chisqlo = PSNID_BIGN;

    if (ipass == 0) {        // coarse grid

      // redshift
      if ( DOPRIOR_ZSPEC ) {  // spectro-z
	minz  = this_z;
	maxz  = this_z;
	zstep = 1;
      } else {
	minz  = 1;
	maxz  = PSNID_MAXNZ;
	zstep = PSNID_ZRBN[ipass];
      }

      // color
      mina  = 1;
      maxa  = PSNID_MAXNA;
      astep = PSNID_ARBN[ipass];

      // peak date
      //      PSNID_FIRST_MJD = data_mjd[0];
      psnid_best_peak_mjd_guess(nobs, data_filt, 
				data_mjd, data_fluxcal, data_fluxcalerr);
      PSNID_FIRST_MJD = PSNID_PEAK_GUESS;
      mini  = 1;
      istep = PSNID_TSTEP[ipass];  // step size in days
      PSNID_PEAK_START = (float) (int) PSNID_FIRST_MJD + PSNID_TSTART[ipass];
      PSNID_PEAK_STOP  = (float) (int) PSNID_FIRST_MJD + PSNID_TSTOP[ipass];
      //      PSNID_PEAK_START =  PSNID_FIRST_MJD + PSNID_TSTART[ipass];
      //      PSNID_PEAK_STOP  =  PSNID_FIRST_MJD + PSNID_TSTOP[ipass];
      
      npeak            = (int)((PSNID_PEAK_STOP - PSNID_PEAK_START)/istep);
      if (npeak<1) npeak=1;
      maxi  = npeak;
      //    printf("SMJD = %f  EMJD = %f  npeak = %d\n",PSNID_PEAK_START, PSNID_PEAK_STOP, npeak);

      // dmu
      minu  = 1;
      maxu  = PSNID_MAXNU;
      ustep = PSNID_URBN[ipass];

    } else if (ipass > 0 && ipass < PSNID_NITER) {

      // finer grid
      // redshift

      if (DOPRIOR_ZSPEC ) {  // spectro-z
	minz  = this_z;
	maxz  = this_z;
	zstep = 1;
      } else {
	indTmp = ind[ipass-1][itype][PSNID_PARAM_LOGZ] ; 
	ZWID   = PSNID_ZWID[ipass-1] ;    ZRBN = PSNID_ZRBN[ipass-1] ;
	jtmp   = indTmp - (ZWID*ZRBN);	
	minz   = (jtmp < 1          ) ? 1           : jtmp ;
	jtmp   = indTmp + (ZWID*ZRBN);	
	maxz   = (jtmp > PSNID_MAXNZ) ? PSNID_MAXNZ : jtmp ;
	zstep = PSNID_ZRBN[ipass-1];
      }
      if (minz > maxz) {
	minz  = maxz;
	zstep = 1;
      }

      indTmp = ind[ipass-1][itype][PSNID_PARAM_COLORPAR] ;
      AWID   = PSNID_AWID[ipass-1] ;   ARBN = PSNID_ARBN[ipass-1] ;
      jtmp   = indTmp - (AWID*ARBN);	
      mina   = (jtmp < 1          ) ? 1           : jtmp ;     
      jtmp   = indTmp + (AWID*ARBN);	
      maxa   = (jtmp > PSNID_MAXNA) ? PSNID_MAXNA : jtmp ;

      //    astep = PSNID_ARBN[ipass];
      astep = PSNID_ARBN[ipass-1];
      if (mina > maxa) {
	mina  = maxa;
	astep = 1;
      }

      // peak date
      mini  = 1;
      istep = PSNID_TSTEP[ipass];  // step size in days
      PSNID_PEAK_START = (float) (int) PSNID_PEAK_GUESS + PSNID_TSTART[ipass];
      PSNID_PEAK_STOP  = (float) (int) PSNID_PEAK_GUESS + PSNID_TSTOP[ipass];
      npeak         = (int) ((PSNID_PEAK_STOP - PSNID_PEAK_START)/istep);
      if (npeak<1) npeak=1;
      maxi  = npeak;

      // dmu
      indTmp = ind[ipass-1][itype][PSNID_PARAM_DMU] ;
      UWID   = PSNID_UWID[ipass-1] ;   URBN = PSNID_URBN[ipass-1] ;
      jtmp   = indTmp - (UWID*URBN);	
      minu   = (jtmp < 1          ) ? 1           : jtmp ;     
      jtmp   = indTmp + (UWID*URBN);	
      maxu   = (jtmp > PSNID_MAXNU) ? PSNID_MAXNU : jtmp ;
      
      //    ustep = PSNID_URBN[ipass];
      ustep = PSNID_URBN[ipass-1];
      if (minu > maxu) {
	minu  = maxu;
	ustep = 1;
      }

    } else if (ipass == PSNID_NITER) {

      ///////////////////////////////
      ///////  Bayesian mode  ///////
      ///////////////////////////////
      // redshift
      if ( DOPRIOR_ZSPEC ) {  // spectro-z
	minz  = this_z;
	maxz  = this_z;
	zstep = 1;
      } else {
	indTmp = ind[ipass-1][itype][PSNID_PARAM_LOGZ] ; 
	ZWID   = PSNID_ZWID[ipass-1] ;    ZRBN = PSNID_ZRBN[ipass-1] ;
	jtmp   = indTmp - (ZWID*ZRBN);	
	minz   = (jtmp < 1          ) ? 1           : jtmp ;
	jtmp   = indTmp + (ZWID*ZRBN);	
	maxz   = (jtmp > PSNID_MAXNZ) ? PSNID_MAXNZ : jtmp ;
	zstep = PSNID_ZRBN[ipass-1];
	//	minz  = 1;
	//	maxz  = PSNID_MAXNZ;
	//	zstep = 1;
      }
      if (minz > maxz) {
	minz  = maxz;
	zstep = 1;
      }
      
      indTmp = ind[ipass-1][itype][PSNID_PARAM_COLORPAR] ;
      AWID   = PSNID_AWID[ipass-1] ;  ARBN = PSNID_ARBN[ipass-1] ;
      jtmp   = indTmp - (AWID*ARBN);	
      mina   = (jtmp < 1          ) ? 1           : jtmp ;     
      jtmp   = indTmp + (AWID*ARBN);	
      maxa   = (jtmp > PSNID_MAXNA) ? PSNID_MAXNA : jtmp ; 
      astep  = PSNID_ARBN[ipass-1];
      if (mina > maxa) {
	mina  = maxa;
	astep = 1;
      }

      // peak date
      mini  = 1;
      istep = PSNID_TSTEP[ipass-1];  // step size in days
      PSNID_PEAK_START = (float) (int) PSNID_PEAK_GUESS + PSNID_TSTART[ipass-1];
      PSNID_PEAK_STOP  = (float) (int) PSNID_PEAK_GUESS + PSNID_TSTOP[ipass-1];
      npeak         = (int)((PSNID_PEAK_STOP - PSNID_PEAK_START)/istep);
      if (npeak<1) npeak=1;
      maxi  = npeak;

      // dmu
      indTmp = ind[ipass-1][itype][PSNID_PARAM_DMU] ;
      UWID   = PSNID_UWID[ipass-1] ;
      URBN   = PSNID_URBN[ipass-1] ;
      jtmp   = indTmp - (UWID*URBN);	
      minu   = (jtmp < 1          ) ? 1           : jtmp ;     
      jtmp   = indTmp + (UWID*URBN);	
      maxu   = (jtmp > PSNID_MAXNU) ? PSNID_MAXNU : jtmp ;
      ustep = PSNID_URBN[ipass-1] ; 
      if (minu > maxu) {
	minu  = maxu;
	ustep = 1;
      }

      ///////////////////////////////

    } else if (ipass > PSNID_NITER) {

      // Aug 17 2014 RK - use this to set optDebug for the chi2 function.
      minz = maxz =  ind[ipass-1][itype][PSNID_PARAM_LOGZ] ;
      mind = maxd =  ind[ipass-1][itype][PSNID_PARAM_SHAPEPAR] ;
      mina = maxa =  ind[ipass-1][itype][PSNID_PARAM_COLORPAR] ;
      mini = maxi =  ind[ipass-1][itype][PSNID_PARAM_TMAX] ;
      minu = maxu =  ind[ipass-1][itype][PSNID_PARAM_DMU] ;
      zstep = astep = ustep = istep = dstep = 1;
    }

    evidence[zpind][itype] = 0.0;

    // dmu grid
    if (itype != PSNID_ITYPE_SNIA && PSNID_FITDMU_CC == 0) {
      minu  = 1 ;
      maxu  = 1 ;
    }
    if (PSNID_FITDMU == 0) {
      minu  = 1;
      maxu  = PSNID_MAXNU;    // 1
      ustep = PSNID_MAXNU;    // 1
    }

    /*
    // use this to force full dmu grid
    minu  = 1;
    maxu  = PSNID_MAXNU;
    ustep = 1;
    */

    /*
    printf("\t itype = %d, ipass = %d"
	   "   minz,maxz,zstep = %d %d %d"
	   "   minu,maxu,ustep = %d %d %d"
	   "   mind,maxd,dstep = %d %d %d"
	   "   mina,maxa,astep = %d %d %d"
	   "   mini,maxi,istep = %d %d %7.1f %10.2f\n",
	   itype, ipass,
	   minz, maxz, zstep,
	   minu, maxu, ustep,
	   mind, maxd, dstep,
	   mina, maxa, astep,
	   mini, maxi, istep, PSNID_PEAK_START);
    fflush(stdout);
    */

    /**********************************************************/
    /*****  compare data with grid of light curve models  *****/
    /**********************************************************/
    for (z = minz; z <= maxz; z = z + zstep) {        // redshift
      chisqloz = PSNID_BIGN;

      // Compute photo-z redshift prior contribution to chi2
      // Do chi2-calc here before u,d,a,i loops below (RK)
      if ( DOPRIOR_ZPHOT ) {
	DZ      = ZPRIOR - z_grid[z] ;  
	ZSIG    = DZ/ZPRIOR_ERR ;
	chisq_z = PSNID_INPUTS.WGT_ZPRIOR * (ZSIG*ZSIG) ;
      }
      else  { 
	chisq_z = 0.0 ; 
      }


      for (u = minu; u <= maxu; u = u + ustep) {      // dmu
	chisqlozu = PSNID_BIGN;
	ushift = u_grid[u];
	// for Ia opton, set ushift=0 (9/15/2017)

	for (d = mind; d <= maxd; d = d + dstep) {    // shapepar

	  if ( AVOID_SIMCHEAT ) {
	    isp         = PSNID_NONIA_ABSINDEX[itype][d]; 
	    NON1A_INDEX = 
	      SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NON1A_INDEX[isp];

	    if ( NON1A_INDEX == DATA_PSNID_DOFIT.SIM_NON1A_INDEX ) 
	      { continue ; }
	  } // end AVOID_SIMCHEAT


	  for (a = mina; a <= maxa; a = a + astep) {  // colorpar

	    for (i = mini; i <= maxi; i++) {           // peak MJD (RK fix?)
	      //   for (i = mini; i < maxi; i++) {  // peak MJD (bug?)

	      // shift model along time axis
	      peak_guess   = PSNID_PEAK_START + i*istep;
	      for (t = 1; t <= PSNID_MAXND; t++) {
		fit_epoch[t] = PSNID_MODEL_EPOCH[d][z][t] + peak_guess;
		for (f = 1; f <= PSNID_NFILTER; f++) {
		  // apply extinction and dmu
		  fit_mag[f][t]    = PSNID_MODEL_MAG[f][d][z][t] -
		    (c_grid[a]-PSNID_BASE_COLOR[itype])*PSNID_MODEL_EXTINCT[f][d][z][t] + ushift;
		  fit_magerr[f][t] = PSNID_MODEL_MAGERR[f][d][z][t];
		}
	      }
	      
	      optDebug = (ipass > PSNID_NITER && itype == 0 ) ; 

	      // calculate chi-squared between data and model
	      chisq = 0.0;
	      ngood = 0;
	      //	      printf("\t z,d,a,i = %d %d %d %d  peak_mjd = %10.2f\n",
	      //		     z,d,a,i, peak_guess);
	      psnid_best_calc_chisq(nobs, useobs, data_filt,
				    data_mjd, data_fluxcal, data_fluxcalerr,
				    fit_epoch, fit_mag, fit_magerr,
				    &ngood, &chisq, optDebug, 0);

	      //////////////////////////
	      //  PRIORS

	      // photo-z(host)
	      if ( DOPRIOR_ZPHOT ) { chisq += chisq_z ; } 

	      // AV
	      if (PSNID_USE_AV_PRIOR == 1) {
		pav = psnid_best_avprior1(itype, c_grid[a]);
		chisq += -2.0*PSNID_INPUTS.AV_PRIOR_STR*log(pav);
	      }


	      // make sure the fit doesn't favor a model with all 99.99
	      if (ngood == 0) {
		ngood = 1;
		chisq = 9999.99;
	      }


	      ///////////////////////
	      // Bayesian evidence //
	      if (ipass == PSNID_NITER) {
		if (chisq < 10000.) {
		  wgt = psnid_best_ratePrior(itype,d,z_grid[z]); // Sep 6 2013 - RK
		  evidence[zpind][itype] += wgt * exp(-(chisq)/2.);
		}
	      }
	      ///////////////////////

	      //	    if (itype == 0 && ipass == 2) {
	      //	      printf("zdaui = %d %d %d %d %d\n",z,d,a,u,i);
	      //	    }

	      // if fit is better, replace chisq and indices
	      if (chisq < chisqlo && ipass <= PSNID_NITER ) {  
		chisqlo   = chisq;
		// xxx		dmulo     = ushift;
		thisngood = ngood;
		ind[ipass][itype][PSNID_PARAM_LOGZ]     = (z > 0) ? z : 1;
		ind[ipass][itype][PSNID_PARAM_SHAPEPAR] = (d > 0) ? d : 1;
		ind[ipass][itype][PSNID_PARAM_COLORPAR] = (a > 0) ? a : 1;
		ind[ipass][itype][PSNID_PARAM_TMAX]     = (i > 0) ? i : 1;
		ind[ipass][itype][PSNID_PARAM_DMU]      = (u > 0) ? u : 1;
	      }

	      if (chisq < chisqloz) {    // min chi-sq for this z
		chisqloz = chisq;
	      }
	      if (chisq < chisqlozu) {   // min chi-sq for this z and dmu
		chisqlozu = chisq;
	      }

	    }
	  }
	}
      }
    }
    /**********************************************************/

    // zero index is not allowed
    psnid_best_check_ind_bounds(ipass, itype, ind);


    /*
      printf("     z (min,max,step) = %3d %3d %3d   %3d\n",minz, maxz, zstep, (int) (maxz-minz+1)/zstep);
      printf("   dmu (min,max,step) = %3d %3d %3d   %3d\n",minu, maxu, ustep, (int) (maxu-minu+1)/ustep);
      printf("  dm15 (min,max,step) = %3d %3d %3d   %3d\n",mind, maxd, dstep, (int) (maxd-mind+1)/dstep);
      printf("    av (min,max,step) = %3d %3d %3d   %3d\n",mina, maxa, astep, (int) (maxa-mina+1)/astep);
      printf("  tmax (min,max,step) = %3d %3d %3d   %3d\n",mini, maxi, 1, (int) maxi-mini+1);
      printf("                Total =           %7d\n",
      (int) (maxz-minz+1)/zstep * (int) (maxu-minu+1)/ustep * (int) (maxd-mind+1)/dstep *
      (int) (maxa-mina+1)/astep * (int) maxi-mini+1);
    

    
    printf(" xxx ---------------------------- \n");
    printf(" xxx ipass = %d  itype = %d \n", ipass, itype);
    printf(" xxx zind = %3d, dind = %3d, aind = %3d, iind = %3d, uind = %3d\n",
	   ind[ipass][itype][PSNID_PARAM_LOGZ],
	   ind[ipass][itype][PSNID_PARAM_SHAPEPAR],
	   ind[ipass][itype][PSNID_PARAM_COLORPAR],
	   ind[ipass][itype][PSNID_PARAM_TMAX],
	   ind[ipass][itype][PSNID_PARAM_DMU]);
	   fflush(stdout);
    
    */

    /*
      if (thisngood <= 4 + PSNID_FITDMU) {
      // reduced chi-sq not defined; output chi-sq instead
      minchisq[itype] = chisqloz;
      } else {
      minchisq[itype] = chisqlo/(thisngood - (4 + PSNID_FITDMU));
      }
    */

    PSNID_BEST_RESULTS.MINCHISQ[zpind][itype] = chisqlo;
    //  if (ipass == PSNID_NITER)
    //    printf("COMMENT: chisqlo = %8.2f;   ngood = %3d\n", chisqlo, ngood);

    if ( ipass <= PSNID_NITER ) {
      PSNID_PEAK_GUESS = PSNID_PEAK_START + 
	ind[ipass][itype][PSNID_PARAM_TMAX]*istep;
    }

    if (ipass == PSNID_NITER) 
      { PSNID_PEAKMJD = PSNID_PEAK_GUESS; }


    /********************/
    /* begin diagnostic */
    /*
      if (itype == 0) {
      //  if (itype < 4) {
      z = ind[ipass][itype][PSNID_PARAM_LOGZ];
      d = ind[ipass][itype][PSNID_PARAM_SHAPEPAR];
      a = ind[ipass][itype][PSNID_PARAM_COLORPAR];
      i = ind[ipass][itype][PSNID_PARAM_TMAX];


      printf("\n");
      printf(" zgrid   =   %d %d %d\n",minz,maxz,zstep);
      printf(" ugrid   =   %d %d %d\n",minu,maxu,ustep);
      printf(" dgrid   =   %d %d %d\n",mind,maxd,dstep);
      printf(" agrid   =   %d %d %d\n",mina,maxa,astep);

      printf("\n");
      printf(" ipass = %d\n", ipass);
      printf("ind[%d] = %d %d %d %d   (z d a i)\n", ipass,
      ind[ipass][itype][PSNID_PARAM_LOGZ],ind[ipass][itype][PSNID_PARAM_SHAPEPAR],
      ind[ipass][itype][PSNID_PARAM_COLORPAR],ind[ipass][itype][PSNID_PARAM_TMAX]);
      printf("     z = %f\n", z_grid[ind[ipass][itype][PSNID_PARAM_LOGZ]]);
      printf("     d = %d\n", ind[ipass][itype][PSNID_PARAM_SHAPEPAR]);
      printf("    av = %f\n", c_grid[ind[ipass][itype][PSNID_PARAM_COLORPAR]]);
      printf("  Tmax = %f\n", PSNID_PEAK_GUESS);
	 
      printf("thisngood = %d\n", thisngood);
      printf("minchisq  = %f\n", PSNID_BEST_RESULTS.MINCHISQ[zpind][itype]);
      printf("     dmu  = %f\n", dmulo);
      }


      if (ipass == PSNID_NITER) {

      int this_filt;
      this_filt = 4;

      for (i=0; i<nobs; i++) {
      if (data_filt[i] == this_filt) {  // r-band
      printf(" data  %1d  %12.4f   %10.2f %10.2f\n",
      data_filt[i], data_mjd[i],
      data_fluxcal[i], data_fluxcalerr[i]
      );
      }
      }

      double this_mag, this_magerr, model_flux, model_fluxe;
      for (i=1; i<=PSNID_MAXND; i++) {
      this_mag    = PSNID_MODEL_MAG[this_filt-1][d][z][i] -
      c_grid[a]*PSNID_MODEL_EXTINCT[this_filt-1][d][z][i];
      this_magerr = PSNID_MODEL_MAGERR[this_filt-1][d][z][i];
      psnid_pogson2fluxcal(this_mag, this_magerr,
      &model_flux, &model_fluxe);
      printf("model  %1d  %12.4f   %8.3f %8.3f     %10.2f %10.2f\n",
      this_filt, PSNID_PEAK_GUESS+PSNID_MODEL_EPOCH[d][z][i],
      this_mag, this_magerr, model_flux, model_fluxe
      );
      }
      }  
    */
    /* end diagnostic */
    /******************/

    //  printf("%12.2f\n",PSNID_PEAK_GUESS);
    //  printf("dmulo = %8.2f\n",dmulo);

  }  // loop over passes


  PSNID_BEST_RESULTS.NGOOD[zpind][itype] = thisngood;
  //  printf("\t thisngood = %d\n", thisngood);
  //  fflush(stdout);

  free_dvector(fit_epoch, 1,PSNID_MAXND);
  free_dmatrix(fit_mag, 1,PSNID_NFILTER, 1,PSNID_MAXND);
  free_dmatrix(fit_magerr, 1,PSNID_NFILTER, 1,PSNID_MAXND);
  free_dvector(c_grid, 1,PSNID_MAXNA);
  free_dvector(z_grid, 1,PSNID_MAXNZ);
  free_dvector(u_grid, 1,PSNID_MAXNU);


  return;
}
// end of psnid_best_grid_compare



/**********************************************************************/
void psnid_best_calc_chisq(int nobs, int *useobs,
			   int *data_filt, double *data_mjd,
			   double *data_fluxcal, double *data_fluxcalerr,
			   double *model_mjd,
			   double **model_mag, double **model_magerr,
			   int *ngood, double *chisq, int optDebug,
			   int doReject)
/*

  Input:
    int nobs                 =  number of filter epochs       = NOBS
    double *data_filt        =  vector of filter ID           = IFILTOBS = 0..
    double *data_mjd         =  vector of data MJD            = MJD
    double *data_fluxcal     =  vector of data fluxcal        = FLUX
    double *data_fluxcalerr  =  vector of data fluxcal error  = FLUXERR
    double *model_mjd        =  vector of model MJD
    double **model_mag       =  array of model mag [filter][epoch]
    double **model_magerr    =  array of model mag error [filter][epoch]

  Output:
    int ngood                =  number of good filter epochs
    double chisq             =  chi-squared



  Sep 2013 RK : cleanup for speed. 
    -  Replace pow(X,2.0) with X*X
    -  Compute tFrac only once inside loop

  Feb 5 2014 MS: support for outlier rejection
    - replaced "continue" statements

 */

/**********************************************************************/
{
  int t, f, this_t=0, count=0, this_filt;
  int    USEMJD, NOBS_FILT[MXFILTINDX] ;
  double this_mag=30.0, this_magerr=0.0;
  double chisq_tmp=0.0, CHISQ_FILT[MXFILTINDX];
  double data_flux, data_fluxe, model_flux, model_fluxe, mjd;
  double MJDDIF_MODEL, MJDDIF_DATA, tFrac, MAGDIF, ERRDIF ;
  double SQFERR, FDIF ;
  int *chisq_ind;
  double *chisq_obs;


  chisq_obs = dvector(1,nobs);
  chisq_ind = ivector(1,nobs);
  for (t=1; t<=nobs; t++) {
    chisq_obs[t] = 0.0;
    chisq_ind[t] = 0;
  }

  *chisq = 0.0;

  if ( optDebug ) {
    for(this_filt=0; this_filt<MXFILTINDX; this_filt++ ) {
      NOBS_FILT[this_filt]  = 0 ; 
      CHISQ_FILT[this_filt] = 0.; 
    }
  }

  if (optDebug > 0) { printf(" CCCC: --------------- \n"); fflush(stdout);  }

  for (t = 0; t < nobs; t++) {  // filter-epochs

    f = data_filt[t]+1;   // data_filt = 0 .. NFILTER-1, so add 1
    this_filt = PSNID_INPUTS.IFILTLIST[f-1];

    // check if this filter is used
    //    if ( PSNID_INPUTS.USEFILT[this_filt] == 0 ) { continue ; }
    //    if ( useobs[t]                       != 1 ) { continue ; }

    if ( PSNID_INPUTS.USEFILT[this_filt] == 1 && useobs[t] == 1 ) {

      model_flux  = model_fluxe = 0.0 ;
    
      // check if MJD is within range
      mjd    = data_mjd[t] ;
      USEMJD = ( mjd >= model_mjd[1]  &&  mjd < model_mjd[PSNID_MAXND] ) ;
    
      if ( USEMJD ) {

	/********************************************
	  data point is within model time range
	********************************************/

	// First get model mags

	hunt(model_mjd, PSNID_MAXND, data_mjd[t], &this_t);

	// ----
	if (model_mag[f][this_t]    < PSNID_GOODMAG_HI &&
	    model_mag[f][this_t]    > PSNID_GOODMAG_LO &&
	    model_mag[f][this_t+1]  < PSNID_GOODMAG_HI &&
	    model_mag[f][this_t+1]  > PSNID_GOODMAG_LO &&
	    model_magerr[f][this_t] < PSNID_GOODMAGERR_HI &&
	    model_magerr[f][this_t] > PSNID_GOODMAGERR_LO &&
	    model_magerr[f][this_t+1] < PSNID_GOODMAGERR_HI &&
	    model_magerr[f][this_t+1] > PSNID_GOODMAGERR_LO) {
	  // model mag is within valid range (defined in psnid_tools.h)
	
	  // RK - compute tFrac just once for both MAG and MAGERR 
	  MJDDIF_DATA  = data_mjd[t] - model_mjd[this_t] ;
	  MJDDIF_MODEL = model_mjd[this_t+1] - model_mjd[this_t] ;
	  tFrac        = MJDDIF_DATA / MJDDIF_MODEL ;
	  MAGDIF       = model_mag[f][this_t+1]    - model_mag[f][this_t] ;
	  ERRDIF       = model_magerr[f][this_t+1] - model_magerr[f][this_t] ;

	  this_mag    = model_mag[f][this_t]    + MAGDIF * tFrac ;
	  this_magerr = model_magerr[f][this_t] + ERRDIF * tFrac ;

	  psnid_pogson2fluxcal(this_mag, this_magerr,
			       &model_flux, &model_fluxe); 
      
	} else {
	
	  // invalid model mag
	  model_flux  = model_fluxe = 0.0 ;  
	}
	// ----

      
      } else {    // USEMJD
       
	// invalid MJD
	model_flux  = model_fluxe = 0.0 ;
      
      }  // end of USEMJD block


      // sum chi2 outside USEMJD if-block
      if (data_fluxcalerr[t] > 0.0 ) {
	
	data_flux  = data_fluxcal[t];
	data_fluxe = data_fluxcalerr[t];
	
	FDIF      = model_flux - data_flux ;
	SQFERR    = model_fluxe*model_fluxe + data_fluxe*data_fluxe ;
	chisq_tmp = (FDIF * FDIF) / SQFERR ;      

	chisq_obs[t+1] = chisq_tmp;
	*chisq        += chisq_tmp ; // increment output chisq
      
	count++;

	if ( optDebug ) { 
	  NOBS_FILT[this_filt]++ ; 
	  CHISQ_FILT[this_filt]+= chisq_tmp ; 
	}
      
      }  // end of fluxerr > 0 
    

      if (optDebug > 0) {
	printf(" CCCC: filter,mjd,flux,fluxe,model,modele,chisq = \n");
	printf(" CCCC:    "
	       "%d %10.2f  %7.1f %6.1f  %7.1f  %6.1f  %8.3f\n",
	       this_filt, data_mjd[t],
	       data_fluxcal[t], data_fluxcalerr[t], 
	       model_flux, model_fluxe, chisq_tmp);  fflush(stdout);    
	fflush(stdout);
      }

    }
  }   // end of nobs loop


  // --------------------

  if (optDebug > 0 )  {
    printf(" CCCC:  chisq(tot) = %le\n", *chisq); fflush(stdout);  
    int ifilt;
    for(ifilt=0; ifilt < PSNID_NFILTER; ifilt++ ) {
      this_filt = PSNID_INPUTS.IFILTLIST[ifilt];  // absolute filter index
      printf(" CCCC: filter=%d : NOBS = %3d  chi2 = %f \n",
	     this_filt, NOBS_FILT[this_filt], CHISQ_FILT[this_filt] );
      printf(" CCCC: ------------------------------------ \n");
      fflush(stdout);
    }
  }

  *ngood += count ;


  if ( doReject == 1 ) {
    // reject up to PSNID_INPUTS.NREJECT_OUTLIER points if
    // chisq_obs >= PSNID_INPUTS.CHISQMIN_OUTLIER
    if ( PSNID_INPUTS.NREJECT_OUTLIER > 0 ) {
      indexx(nobs, chisq_obs, chisq_ind);
      for (t = nobs-PSNID_INPUTS.NREJECT_OUTLIER+1; t <= nobs; t++) {
	// NOTE chisq_obs(1,nobs) and useobj(0,nobs)
	if (chisq_obs[chisq_ind[t]] >= PSNID_INPUTS.CHISQMIN_OUTLIER) {
	  useobs[chisq_ind[t]-1] = 2;
	}
      }

      /*
      // debug print statements
      for (t = 1; t <= nobs; t++) {
	f = data_filt[chisq_ind[t]-1]+1;   // data_filt = 0 .. NFILTER-1, so add 1
	this_filt = PSNID_INPUTS.IFILTLIST[f-1];
	printf("\t t = %3d  f = %2d  this_filt = %2d  chisq_ind = %3d  "
	       "useobs = %2d  chisq_obs = %8.3f  "
	       "USEFILT = %d  "
	       "MJD = %f  FLUX = %f  FLUXERR = %f  "
	       "\n",
	       t, f, this_filt, chisq_ind[t],
	       useobs[chisq_ind[t]-1], chisq_obs[chisq_ind[t]],
	       PSNID_INPUTS.USEFILT[this_filt],
	       data_mjd[chisq_ind[t]-1], data_fluxcal[chisq_ind[t]-1],
	       data_fluxcalerr[chisq_ind[t]-1]);
	fflush(stdout);
      }
      */

    }
  }

  free_dvector(chisq_obs, 1,nobs);
  free_ivector(chisq_ind, 1,nobs);


  return;
}
// end of psnid_best_calc_chisq


/**********************************************************************/
int psnid_best_flag_outlier(int itype, int zpind, int nobs,
			    int *data_filt, double *data_mjd,
			    double *data_fluxcal, double *data_fluxcalerr,
			    int *useobs, int ***ind)
/**********************************************************************/
{
  int z, d, a, u;
  int t, f;
  int ngood=0, optDebug=0, nout=0;
  double *fit_epoch, **fit_mag, **fit_magerr;
  double *c_grid, *u_grid;
  double ushift, peak_guess, chisq;


  fit_epoch  = dvector(1,PSNID_MAXND);
  fit_mag    = dmatrix(1,PSNID_NFILTER, 1,PSNID_MAXND);
  fit_magerr = dmatrix(1,PSNID_NFILTER, 1,PSNID_MAXND);

  z = ind[PSNID_NITER-1][itype][PSNID_PARAM_LOGZ];
  d = ind[PSNID_NITER-1][itype][PSNID_PARAM_SHAPEPAR];
  a = ind[PSNID_NITER-1][itype][PSNID_PARAM_COLORPAR];
  //  i = ind[PSNID_NITER-1][itype][PSNID_PARAM_TMAX];
  u = ind[PSNID_NITER-1][itype][PSNID_PARAM_DMU];

  c_grid   = dvector(1,PSNID_MAXNA);
  u_grid   = dvector(1,PSNID_MAXNU);
  psnid_best_get_c_grid(c_grid);
  psnid_best_get_u_grid(u_grid);

  ushift = u_grid[u];

  // shift model along time axis
  //  peak_guess   = PSNID_PEAK_START + i*istep;
  peak_guess   = PSNID_PEAKMJD;
  for (t = 1; t <= PSNID_MAXND; t++) {
    fit_epoch[t] = PSNID_MODEL_EPOCH[d][z][t] + peak_guess;
    for (f = 1; f <= PSNID_NFILTER; f++) {
      // apply extinction and dmu
      fit_mag[f][t]    = PSNID_MODEL_MAG[f][d][z][t] -
	(c_grid[a]-PSNID_BASE_COLOR[itype])*PSNID_MODEL_EXTINCT[f][d][z][t] + ushift;
      fit_magerr[f][t] = PSNID_MODEL_MAGERR[f][d][z][t];
    }
  }

  // calculate chi-squared between data and model
  chisq = 0.0;
  ngood = 0;
  //	      printf("\t z,d,a,i = %d %d %d %d  peak_mjd = %10.2f\n",
  //		     z,d,a,i, peak_guess);
  psnid_best_calc_chisq(nobs, useobs, data_filt,
			data_mjd, data_fluxcal, data_fluxcalerr,
			fit_epoch, fit_mag, fit_magerr,
			&ngood, &chisq, optDebug, 1 );

  for (t = 0; t < nobs; t++) {
    if (useobs[t] == 2) {
      nout++;
    }
  }

  free_dvector(fit_epoch,  1,PSNID_MAXND);
  free_dmatrix(fit_mag,    1,PSNID_NFILTER, 1,PSNID_MAXND);
  free_dmatrix(fit_magerr, 1,PSNID_NFILTER, 1,PSNID_MAXND);

  free_dvector(u_grid, 1,PSNID_MAXNU);
  free_dvector(c_grid, 1,PSNID_MAXNA);


  return nout;
}
// end of psnid_best_flag_outlier


/**********************************************************************/
void psnid_best_peak_mjd_guess(int nobs, int *data_filt, double *data_mjd,
			       double *data_fluxcal,
			       double *data_fluxcalerr)
/**********************************************************************/
{
  int i;
  //  double mjdsum, snrsum, maxsnr=0.0;
  double maxsnr=0.0;


  /*
  mjdsum = 0.0;
  snrsum = 0.0;
  for (i=0; i<nobs; i++) {
    mjdsum += data_mjd[i]*(data_fluxcal[i]/data_fluxcalerr[i]);
    snrsum += (data_fluxcal[i]/data_fluxcalerr[i]);
  }
  printf("\t MJD_GUESS = %10.2f\n", mjdsum/snrsum);
  PSNID_PEAK_GUESS = mjdsum/snrsum;
  */

  maxsnr=0.0;
  PSNID_PEAK_GUESS = data_mjd[0];
  for (i=0; i<nobs; i++) {
    //    printf("data,err = %8.2f %8.2f   %8.2f\n",
    //	   data_fluxcal[i],data_fluxcalerr[i],data_fluxcal[i]/data_fluxcalerr[i]);
    if (data_fluxcal[i]/data_fluxcalerr[i]>maxsnr) {
      PSNID_PEAK_GUESS = data_mjd[i];
      maxsnr = data_fluxcal[i]/data_fluxcalerr[i];
    }
  }
  //  printf("\t MJD_GUESS = %10.2f\n", PSNID_PEAK_GUESS);


  return;
}

// =============================================		 
double psnid_best_ratePrior(int itype, int indx_template, double z) {
  // Inputs:
  //   itype         = 0 for SNIA, >0 for CC
  //   indx_template = SNCC sparse template index
  //   z             = redshift (not the prior index)
  //

  double rate ;
  double a,b,zmax, z1, z1pow, frac_subtype ;
  int    INDX_SPARSE ;

  // -------------- BEGIN -------------
  
  rate = 1.0 ;

  if ( PSNID_INPUTS.OPT_RATEPRIOR == 0 ) { return rate ; }

  if ( itype == PSNID_ITYPE_SNIA ) {
    a    = PSNID_INPUTS.ZRATEPRIOR_SNIA[0];
    b    = PSNID_INPUTS.ZRATEPRIOR_SNIA[1];
    zmax = PSNID_INPUTS.ZRATEPRIOR_SNIA[2];
    INDX_SPARSE = -9 ; // undefined for Ia
    frac_subtype = 1.0 ;
  }
  else {
    a    = PSNID_INPUTS.ZRATEPRIOR_NONIA[0];
    b    = PSNID_INPUTS.ZRATEPRIOR_NONIA[1];
    zmax = PSNID_INPUTS.ZRATEPRIOR_NONIA[2];

    INDX_SPARSE  = PSNID_NONIA_ABSINDEX[itype][indx_template]; 

    // NONIA_WGT are normalized to sum to 1, so take fraction
    frac_subtype = SNGRID_PSNID[2].NON1A_WGT[INDX_SPARSE] ;
  }

  if ( z < zmax ) 
    {  z1 = 1.0 + z ; }
  else
    {  z1 = 1.0 + zmax ; }

  z1pow = pow(z1,b);  // note: lookup might be faster
  rate  = a * z1pow * frac_subtype;

  /* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  if ( itype > 0 ) {
    printf(" xxx itype=%d INDX_[templ,SPARSE]=%d,%d  frac=%f  rate=%.3le\n",
	   itype, indx_template, INDX_SPARSE, frac_subtype, rate);
    fflush(stdout);
  }
  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

  return rate ;
  
} //end of psnid_best_ratePrior

/**********************************************************************/
void psnid_best_check_lc_quality(int itype, int z, int ***ind,
				 int nepoch, double *mjd, double *flux,
				 double *fluxerr)
/**********************************************************************/
{

  // Jan 4 2014 RK - also compute and store TOBSMIN/MAX
  //                 Replace epoch[i] with scalar Tobs since we
  //                 don't need to store each Tobs

  int i;
  int t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, lcq;
  double zdelay;
  double redshift, shapepar, colorpar, colorlaw, pkmjd, dmu;
  double Tobs, Tmin, Tmax ;

  //  epoch = dvector(1,nepoch);

  // determine best-fit parameters from indices
  psnid_best_calc_fitpar(itype, ind,
			 &redshift, &shapepar, &colorpar, &colorlaw,
			 &pkmjd, &dmu);

  zdelay = 1. + redshift;
  //  zdelay = 1. + ZLC;

  Tmin = 99999. ; Tmax = -99999.;

  for (i = 1; i <= nepoch; i++) {

    // should we bail if mjd[i] < 10 ? RK - Jan 2014

    //xxxx    epoch[i] = mjd[i] - pkmjd; 
    Tobs =  mjd[i] - pkmjd;

    if (Tobs <  -15.0*zdelay                       ) t1 =  4;
    if (Tobs >= -15.0*zdelay && Tobs <  -5.0*zdelay) t2 =  8;
    if (Tobs >=  -5.0*zdelay && Tobs <=  5.0*zdelay) t3 = 32;
    if (Tobs >    5.0*zdelay && Tobs <= 15.0*zdelay) t4 = 16;
    if (Tobs >   15.0*zdelay && Tobs <= 30.0*zdelay) t5 =  2;
    if (Tobs >   30.0*zdelay                       ) t6 =  1;
    
    if ( mjd[i] > 10.0 ) {
      if ( Tobs < Tmin ) { Tmin = Tobs ; } // RK, Jan 2014
      if ( Tobs > Tmax ) { Tmax = Tobs ; }
    }

    //    printf("%12.2f  %+6.2f  %12.2f\n", mjd[i], epoch[i], PEAK_GUESS);
  }
  //  printf("%d %d %d %d\n", t1, t2, t3, t4);


  lcq = t1 + t2 + t3 + t4 + t5 + t6;
  PSNID_BEST_RESULTS.LIGHTCURVE_QUALITY[z][itype] = lcq;
  //  printf("lcq = %2d\n", lcq);


  // load global TOBSMIN/MAX (RK - Jan 2014)
  PSNID_BEST_RESULTS.TOBSMIN[z][itype] = Tmin ;
  PSNID_BEST_RESULTS.TOBSMAX[z][itype] = Tmax ;

  return;
}
// end of check_lc_quality


/**********************************************************************/
void psnid_best_reset_results()
/**********************************************************************/
/***
    Set all variables in PSNID_BEST_RESULTS structure to default values.
    Called for each event
 ***/
{
  int i, j, k, l, ifilt ;


  // ------- BEGIN -------------
  for(i=0; i < PSNID_NTYPES; i++ )  {
    sprintf(PSNID_TYPE_NAME[i], "%s", PSNID_ITYPE_STRING[i] );
    PSNID_BEST_RESULTS.FINAL_NONIA_INDX_SPARSE[i] = -9 ;
    
    for(ifilt=0; ifilt< MXFILTINDX; ifilt++ ) {
      PSNID_BEST_RESULTS.FINAL_PEAKMAG[i][ifilt] = -9.0 ;
    }
  }
  
  sprintf(PSNID_BEST_RESULTS.CID, "UNKNOWN");

  for (i=0; i<PSNID_NZPRIOR; i++) {

    PSNID_BEST_RESULTS.ZPRIOR[i]      = -9.000;
    PSNID_BEST_RESULTS.ZPRIOR_ERR[i]  = -9.000;
    PSNID_BEST_RESULTS.ZPRIOR_DO[i]   = -9;
    PSNID_BEST_RESULTS.FINAL_ITYPE[i] = -9;

    sprintf(PSNID_BEST_RESULTS.FINAL_TYPE[i],  "UNKNOWN");

    for (j=0; j<PSNID_NTYPES; j++) {

      PSNID_BEST_RESULTS.NGOOD[i][j]      = 0;
      PSNID_BEST_RESULTS.MINCHISQ[i][j]   = 9999.999;
      PSNID_BEST_RESULTS.TOBSMIN[i][j]   =  9999.999;
      PSNID_BEST_RESULTS.TOBSMAX[i][j]   = -9999.999;
      PSNID_BEST_RESULTS.FITPROB[i][j]    = -9.000;
      PSNID_BEST_RESULTS.PBAYESIAN[i][j]  = -9.000;
      PSNID_BEST_RESULTS.LIGHTCURVE_QUALITY[i][j]  = -9;

      for (k=0; k<PSNID_NPARAM; k++) {

	PSNID_BEST_RESULTS.FINAL_PARVAL[i][j][k]   = -9.0;
	PSNID_BEST_RESULTS.FINAL_PARERR[i][j][k]   = -9.0;
		
	PSNID_BEST_RESULTS.MINCHISQ_IND[i][j][k]   = 0;
	PSNID_BEST_RESULTS.GRID_FINAL_PAR[i][j][k]  = -9.000;

	for (l=0; l<PSNID_MCMC_NLIMITS; l++) {
	  PSNID_BEST_RESULTS.MCMC_FINAL_PAR[i][j][k][l]  = -9.000;
	}

      }
    }
  }


  return;
}
// end of psnid_best_reset_results


/**********************************************************************/
void psnid_best_set_cid(char *cid)
/**********************************************************************/
{

  sprintf(PSNID_BEST_RESULTS.CID, "%s",  cid); //fix compile bug, RK Dec 2 2012

}
// end of psnid_best_set_cid


/**********************************************************************/
void psnid_best_set_zprior(double *zin, double *zinerr)
/**********************************************************************/
/***

               ?!?!?!? OBSOLETE ?!?!?

    zin[0], zinerr[0] = SNLC_REDSHIFT, SNLC_REDSHIFT_ERR
    zin[1], zinerr[1] = SNHOST_PHOTOZ, SNHOST_PHOTOZ_ERR
    zin[2], zinerr[2] = -9.0000, -9.0000 (for now)

    See psnid.cra.

    Fills in:   PSNID_BEST_RESULTS.ZPRIOR
                PSNID_BEST_RESULTS.ZPRIOR_ERR
                PSNID_BEST_RESULTS.ZPRIOR_DO

 ***/
{
  int i, n=1;

  // flat
  PSNID_BEST_RESULTS.ZPRIOR[0]      = -1.00;
  PSNID_BEST_RESULTS.ZPRIOR_ERR[0]  = -1.00;
  PSNID_BEST_RESULTS.ZPRIOR_DO[0]   =  1;

  for (i=1; i<PSNID_NZPRIOR; i++) {
    PSNID_BEST_RESULTS.ZPRIOR[i]     = zin[i-1];
    PSNID_BEST_RESULTS.ZPRIOR_ERR[i] = zinerr[i-1];
  }

  for (i=1; i<PSNID_NZPRIOR; i++) {    // skip flat
    if (PSNID_BEST_RESULTS.ZPRIOR[i]     > 0.0 && 
	PSNID_BEST_RESULTS.ZPRIOR_ERR[i] > 0.0 ) {
      PSNID_BEST_RESULTS.ZPRIOR_DO[i] = 1 ;
      n++;
    } else {
      PSNID_BEST_RESULTS.ZPRIOR_DO[i] = 0 ;
    }
  }

  PSNID_BEST_RESULTS.NZPRIOR = n;


  return;
}
// end of psnid_best_set_zprior


/**********************************************************************/
void psnid_best_set_zprior_onez(double *zin, double *zinerr)
/**********************************************************************/
/***

    zin[0], zinerr[0] = SNLC_REDSHIFT, SNLC_REDSHIFT_ERR
    zin[1], zinerr[1] = SNHOST_PHOTOZ, SNHOST_PHOTOZ_ERR
    zin[2], zinerr[2] = -9.0000, -9.0000 (for now)

    See psnid.cra.

    Fills in:   PSNID_BEST_RESULTS.ZPRIOR
                PSNID_BEST_RESULTS.ZPRIOR_ERR
                PSNID_BEST_RESULTS.ZPRIOR_DO

    Set up to run PSNID for only a single redshift prior.

    OPT_ZPRIOR = PSNID_ZPRIOR_FLAT   --> flat
    OPT_ZPRIOR = PSNID_ZPRIOR_SNLC   --> spectro
    OPT_ZPRIOR = PSNID_ZPRIOR_HOST   --> host photo-z

    For spectro and host photo-z priors, redshift error range must satisfy:

    PSNID_INPUTS.CUTWIN_ZERR[0] <= zinerr <= PSNID_INPUTS.CUTWIN_ZERR[1]

 ***/
{
  int i; 

  for (i=0; i<PSNID_NZPRIOR; i++) 
    { PSNID_BEST_RESULTS.ZPRIOR_DO[i] = 0 ; } // init (RK)

  if (PSNID_INPUTS.OPT_ZPRIOR == PSNID_ZPRIOR_FLAT) {
    PSNID_BEST_RESULTS.ZPRIOR[0]      = -1.00;
    PSNID_BEST_RESULTS.ZPRIOR_ERR[0]  = -1.00;
    PSNID_BEST_RESULTS.ZPRIOR_DO[0]   =  1 ;
  } else if (PSNID_INPUTS.OPT_ZPRIOR == PSNID_ZPRIOR_SNLC) {
    if (zinerr[0] >= PSNID_INPUTS.CUTWIN_ZERR[0] &&
	zinerr[0] <= PSNID_INPUTS.CUTWIN_ZERR[1]) {
      PSNID_BEST_RESULTS.ZPRIOR[0]      = zin[0];
      PSNID_BEST_RESULTS.ZPRIOR_ERR[0]  = zinerr[0];
      PSNID_BEST_RESULTS.ZPRIOR_DO[0]   = 1 ;
    }
  } else if (PSNID_INPUTS.OPT_ZPRIOR == PSNID_ZPRIOR_HOST) {
    if (zinerr[1] >= PSNID_INPUTS.CUTWIN_ZERR[0] &&
	zinerr[1] <= PSNID_INPUTS.CUTWIN_ZERR[1]) {
      PSNID_BEST_RESULTS.ZPRIOR[0]      = zin[1];
      PSNID_BEST_RESULTS.ZPRIOR_ERR[0]  = zinerr[1];
      PSNID_BEST_RESULTS.ZPRIOR_DO[0]   = 1 ;
    }
  }

  PSNID_BEST_RESULTS.NZPRIOR = 1;

  
  /*
  for (i=0; i<PSNID_NZPRIOR; i++) {
    printf("\t xxx ZPRIOR[z,err,do] = %f %f %d \n",
	   PSNID_BEST_RESULTS.ZPRIOR[i],
	   PSNID_BEST_RESULTS.ZPRIOR_ERR[i],
	   PSNID_BEST_RESULTS.ZPRIOR_DO[i]);
  }
  fflush(stdout); 
  */

  return ;
}
// end of psnid_best_set_zprior_onez


/**********************************************************************/
void psnid_best_flag_epochs(int nobs, int *data_filt, double *data_mjd,
			    double *data_fluxcal,
			    double *data_fluxcalerr,
			    int *useobs)
/**********************************************************************/
{
  int i, t, f, this_filt;
  double *peak_guess, *maxsnr, snr;
  double sum_maxsnr=0.0, sum_peak_guess=0.0, final_peak_guess=0.0;


  maxsnr     = dvector(0, PSNID_NFILTER);
  peak_guess = dvector(0, PSNID_NFILTER);

  //  for (i=0; i<PSNID_NFILTER; i++) {  // bug line
  for (i=0; i<=PSNID_NFILTER; i++) {  // RK:  <= instead of <
    maxsnr[i]     = -9.0;
    peak_guess[i] = -9.0;
  }

  for (t=0; t<nobs; t++) {

    f = data_filt[t]+1;   // data_filt = 0 .. NFILTER-1, so add 1
    //    printf("\t iobs = %d   f = %d\n", t, f);
    this_filt = PSNID_INPUTS.IFILTLIST[f-1];

    if (PSNID_INPUTS.USEFILT[this_filt] == 1) {  // use this filter (from nml)
      snr = data_fluxcal[t] / data_fluxcalerr[t] ;
      if ( snr > maxsnr[data_filt[t]] ) {
      	peak_guess[data_filt[t]] = data_mjd[t];
      	maxsnr[data_filt[t]]     = snr ;
      }

    }
  }

  for (i=0; i<PSNID_NFILTER; i++) {
    //    printf("\t maxsnr=%8.3f peak_guess=%8.3f\n", maxsnr[i],peak_guess[i]);
    if (maxsnr[i]>0.0 && peak_guess[i]>0.0) {
      sum_maxsnr     += maxsnr[i];
      sum_peak_guess += peak_guess[i]*maxsnr[i];
    }
  }
  if (sum_maxsnr>0.0) {
    final_peak_guess = sum_peak_guess/sum_maxsnr;
  }
  //  printf("\t final_peak_guess=%10.3f\n", final_peak_guess);


  // NOTE: put these later into NML
  for (t=0; t<nobs; t++) {

    f = data_filt[t]+1;   // data_filt = 0 .. NFILTER-1, so add 1
    this_filt = PSNID_INPUTS.IFILTLIST[f-1];

    if (PSNID_INPUTS.USEFILT[this_filt] == 1) {

      if (data_mjd[t] > final_peak_guess - 30.0 &&
	  data_mjd[t] < final_peak_guess + 80.0) {
	useobs[t] = 1;
      } else {
	useobs[t] = 0;
      }

    } else {
      useobs[t] = 0;
    }
    //       printf("\t   t=%3d  MJD=%10.3f  useobs=%d\n",
    //    	   t, data_mjd[t], useobs[t]);
  }
  

  free_dvector(maxsnr, 0, PSNID_NFILTER);
  free_dvector(peak_guess, 0, PSNID_NFILTER);


  return;
}
// end of psnid_best_flag_epochs


/**********************************************************************/
void psnid_best_calc_fitpar(int itype, int ***ind,
			    double *redshift, double *shapepar, double *colorpar,
			    double *colorlaw, double *pkmjd, double *dmu)
/**********************************************************************/
{

  // redshift (GRID)
  *redshift = pow(10.,SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_LOGZ][ind[PSNID_NITER][itype][PSNID_PARAM_LOGZ]]);

  // shapepar (GRID)
  *shapepar = SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_SHAPEPAR][ind[PSNID_NITER][itype][PSNID_PARAM_SHAPEPAR]];

  // colorpar (PSNID)
  *colorpar = PSNID_AMIN + PSNID_ASTEP*(ind[PSNID_NITER][itype][PSNID_PARAM_COLORPAR]-1.);

  // colorlaw (GRID, not fit)
  *colorlaw = SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_COLORLAW][1];

  // Tmax (PSNID)
  *pkmjd = PSNID_PEAK_GUESS;

  // dmu (PSNID, optional)
  *dmu = PSNID_UMIN + PSNID_USTEP*(ind[PSNID_NITER][itype][PSNID_PARAM_DMU]-1.);

  return;
}
// end of psnid_best_calc_fitpar


/**********************************************************************/
void psnid_best_check_ind_bounds(int ipass, int itype, int ***ind)
/**********************************************************************/
/***
    Make sure that the best fit indices are within allowed bounds.
 ***/
{
  int i;


  for (i = 0; i < PSNID_NPARAM; i++) {

    if ( ind[ipass][itype][i] < 1 ) {
      //	printf("\t zero index encountered  i = %d  ind = %d\n",
      //	       i, ind[ipass][itype][i]);
      //	fflush(stdout);
      ind[ipass][itype][i] = 1;
    }

    if ( ind[ipass][itype][i] > PSNID_PARAM_MAX_INDEX[i] ) {
      //      printf("\t max index encountered  i = %d  ind = %d\n",
      //	     i, ind[ipass][itype][i]);
      //      fflush(stdout);
      ind[ipass][itype][i] = PSNID_PARAM_MAX_INDEX[i];
    }

  }

  return;
}

// ***************************************
void psnid_best_store_fitResids(int itype, int **obsflag) {

  // Created Oct 2013
  // Store light curve fit-residuals to be used for
  // plotting.  Code moved from SNLCPLOT_PSNID_BEST().
  // Note that store_data() has already been called.
  //
  // Load RESIDS_PSNID_DOFIT[itype].
  //
  // Feb 13, 2014 RK - set REJECT=1 when MAG_ERR<0. Still need to
  //                   set REJECT flag for rejected CHISQ outliers.
  //
  // Aug 17 2014: RK - require .not.REJECT to increment NDOF_FILT
  //                   and CHI2_FILT. Fixes chis/Ndof lables on plots.
  //
  // -------------- BEGIN -----------

  int    NOBS, IFILT, IFILTOBS, obs, IS_BESTFIT, IPAR, REJECT ;
  char  *CCID ;
  double PKMJD, TOBS, TREST, MJD, PULL, CHI2, MAG, MAG_ERR;
  double FITFLUX, FITFLUX_ERR, DATAFLUX, DATAFLUX_ERR ;
  double DIF, SQDIF, SQERR, Redshift ;
  int    izprior = 0 ;
  int    ONE = 1;

  char CFILT[2], CTYPE[20] ;
  //  char fnam[] = "psnid_best_store_fitResids";
  
  // ==========================================================
  // set local pointers; start with those filled in psnid_store_data()
  CCID   = DATA_PSNID_DOFIT.CCID ;
  NOBS   = DATA_PSNID_DOFIT.NOBS ;

  IPAR   = PSNID_PARAM_TMAX ;
  PKMJD  = PSNID_BEST_RESULTS.FINAL_PARVAL[izprior][itype][IPAR];

  IPAR   = PSNID_PARAM_LOGZ ;
  Redshift = PSNID_BEST_RESULTS.FINAL_PARVAL[izprior][itype][IPAR];

  // init filter-dependent variables
  for(IFILT=0; IFILT < MXFILTINDX; IFILT++ ) {
    RESIDS_PSNID_DOFIT[itype].USE_FILT[IFILT]   = 0 ;  
    RESIDS_PSNID_DOFIT[itype].NDOF_FILT[IFILT]  = 0.0 ;  
    RESIDS_PSNID_DOFIT[itype].CHI2_FILT[IFILT]  = 0.0 ;
    RESIDS_PSNID_DOFIT[itype].PKMJD_FILT[IFILT] = 0.0 ; 
  }

  IS_BESTFIT = (itype == PSNID_BEST_RESULTS.FINAL_ITYPE[0] ) ;

  for(obs=0; obs < NOBS; obs++ )  { 
    MJD           = DATA_PSNID_DOFIT.MJD[obs] ;
    TOBS          = MJD - PKMJD ;    
    TREST         = TOBS/(1+Redshift) ;
    IFILT         = DATA_PSNID_DOFIT.IFILT[obs] ;
    IFILTOBS      = DATA_PSNID_DOFIT.IFILTOBS[obs] ;
    sprintf(CFILT,"%c", FILTERSTRING[IFILTOBS] );  // for comments
    DATAFLUX      = DATA_PSNID_DOFIT.FLUXDATA[obs] ;
    DATAFLUX_ERR  = DATA_PSNID_DOFIT.FLUXERR[obs] ;
    REJECT        = 0 ;

    PSNID_BEST_GET_FITFUN(CCID, itype, 
			  IFILT, ONE,          // (I) ifilt, Nepoch
			  &TOBS,               // (I) Tobs
			  &MAG, &MAG_ERR ) ;   // (O) best-fit mag and error

    // convert fitted MAG and MAG_ERR into FLUX and FLUX_ERR
    psnid_pogson2fluxcal(MAG, MAG_ERR, &FITFLUX, &FITFLUX_ERR);

    if ( MAG_ERR > 90.0 ) {  // model is not defined
      FITFLUX = FITFLUX_ERR = 0.0 ; 
      //  REJECT  = 1 ; // Aug 17,2014: not rejected in fit !
    }

    if ( obsflag[itype][obs] != 1 ) {
       REJECT  = 1 ;
       if ( obsflag[itype][obs] == 2 ) {
	 sprintf(CTYPE,"%s", PSNID_ITYPE_STRING[itype]);
	 printf("\t  Rejected %s-%s outlier at MJD=%.3f (%s fit) \n",
		CCID, CFILT, MJD, CTYPE ); fflush(stdout);
       }
    } // end of obsflag check

    DIF   = DATAFLUX - FITFLUX ; 
    SQDIF = DIF*DIF ;
    SQERR = FITFLUX_ERR*FITFLUX_ERR + DATAFLUX_ERR*DATAFLUX_ERR ;
    if ( SQERR > 0.0 )  
      { CHI2 = SQDIF / SQERR ;  PULL = DIF/sqrt(SQERR); }
    else                
      { CHI2 = -9.0 ; PULL = -999.; }

    // filter-dependent variables
    RESIDS_PSNID_DOFIT[itype].USE_FILT[IFILT]   = 1 ;   // set logical flag
    RESIDS_PSNID_DOFIT[itype].PKMJD_FILT[IFILT] = PKMJD ;

    // Aug 17 2014: require not.REJECT to increment NDOF & CHI2
    if ( REJECT == 0 && CHI2 > -0.001 ) {  
      RESIDS_PSNID_DOFIT[itype].NDOF_FILT[IFILT] += 1.0 ; // double, not int
      RESIDS_PSNID_DOFIT[itype].CHI2_FILT[IFILT] += CHI2 ; 
    }

    // epoch-dependent variables
    RESIDS_PSNID_DOFIT[itype].TOBS[obs]          = TOBS ;
    RESIDS_PSNID_DOFIT[itype].REJECT[obs]        = (double)REJECT ;
    RESIDS_PSNID_DOFIT[itype].CHI2[obs]          = CHI2 ;
    RESIDS_PSNID_DOFIT[itype].FITFLUX[obs]       = FITFLUX ;
    RESIDS_PSNID_DOFIT[itype].FITFLUX_ERR[obs]   = FITFLUX_ERR ;

    if ( IS_BESTFIT ) {
      // note there is no itype index here
      F_RESIDS_PSNID_DOFIT.NFIT               = NOBS ;
      F_RESIDS_PSNID_DOFIT.MJD[obs]           = MJD ;
      F_RESIDS_PSNID_DOFIT.TOBS[obs]          = TOBS ;
      F_RESIDS_PSNID_DOFIT.TREST[obs]         = TREST ;
      F_RESIDS_PSNID_DOFIT.CHI2[obs]          = CHI2 ;
      F_RESIDS_PSNID_DOFIT.PULL[obs]          = PULL ;
      F_RESIDS_PSNID_DOFIT.IFILTOBS[obs]      = IFILTOBS ;
      F_RESIDS_PSNID_DOFIT.FITFLUX[obs]       = FITFLUX ;
      F_RESIDS_PSNID_DOFIT.FITFLUX_ERR[obs]   = FITFLUX_ERR ;
      F_RESIDS_PSNID_DOFIT.DATAFLUX[obs]      = DATAFLUX ;
      F_RESIDS_PSNID_DOFIT.DATAFLUX_ERR[obs]  = DATAFLUX_ERR ;
    }

  } // end of obs loop

  // PSNID_INPUTS.NFILT ;

  // -----------------------------------
  // count how many filters are used.
  int USE, NFILT_USE ;
  NFILT_USE = 0 ;
  for(IFILT=0; IFILT < PSNID_INPUTS.NFILT; IFILT++ ) { 
    USE =  RESIDS_PSNID_DOFIT[itype].USE_FILT[IFILT] ;
    if ( USE ) { NFILT_USE++ ; }
  }
  RESIDS_PSNID_DOFIT[itype].NFILT_USE = NFILT_USE ;

} // end of psnid_best_store_fitResids


/*************************************************************/
void SNLCPLOT_PSNID_BEST(int iplot) {
/**************************************************************/

  // Created Feb 9 2013 by R.Kessler
  // Call SNLCPAK_XXX functions for plotting.
  // This replaces the old interface.
  //
  // Mar  4, 2013: add CHI2_FILT[ifilt]
  // Aug 22, 2013 .GRID_FINAL_PAR(grid minimum) -> 
  //              .FINAL_PARVAL(minimum or MCMC)
  //
  // Oct 23, 2013 RK - major overhaul using RESIDS_PSNID_DOFIT struct.
  //                  
  // Jan 08 2020: pass simulated FLUXSIM (true flux) to SNLCPAK.
  //

  int  
    itype, izprior, NOBS, USE
    , NFILT, IFILT, IFILTOBS, NFILT_USE, LDMP
    , ONE = 1
    ,*ptrIFILT, *ptrIFILTOBS, *IFILTLIST 
    ;

  double 
    PKMJD, TOBS, MJD, MAG, MAG_ERR, FITFLUX, FITFLUX_ERR
    ,*ptrTOBS, *ptrMJD, *ptrFLUXDATA, *ptrFLUXSIM, *ptrFLUXERR
    ,*ptrREJ, *ptrCHI2, *ptrDUMERR0
    ;

  char  *CCID ;
  char   fnam[]  = "SNLCPLOT_PSNID_BEST" ;

  // ------------ BEGIN -------------

    
  // get CCID from stored data structure
  CCID     = DATA_PSNID_DOFIT.CCID ;

  // get izprior and itype for this 'iplot'
  psnid_best_get_plotind(iplot, &izprior, &itype);

  if ( PSNID_NGRID[itype] == 0 ) { return ; }
  
  // total number of filters
  NFILT     = PSNID_INPUTS.NFILT ;

  // ========================================================
  // set display-text strings using SNLCPAK_DISPLAYTEXT
  // 1. FUN = Ia(best-fit: yes/no)   z-prior= none or Zerr<Zerr)
  // 2. COL = xx   RV = xx     SHAPE = xx +- xx  
  // 3. DIST = xx   z = xx  chi2/dof = xx/xx

  int  itype_best, OPT_ZPRIOR ;
  char DISPLAYTEXT[80], text_zprior[60], text_bestfit[20];
  double z, zerr;

  itype_best = PSNID_BEST_RESULTS.FINAL_ITYPE[0];
  OPT_ZPRIOR = PSNID_INPUTS.OPT_ZPRIOR ;

  if ( OPT_ZPRIOR == PSNID_ZPRIOR_FLAT ) 
    { sprintf(text_zprior,"none"); }
  else if ( OPT_ZPRIOR == PSNID_ZPRIOR_SNLC ) { 
    z    = DATA_PSNID_DOFIT.REDSHIFT[0] ;
    zerr = DATA_PSNID_DOFIT.REDSHIFT_ERR[0] ;
    sprintf(text_zprior,"%6.4f +- %6.4f", z, zerr );
  }
  else if ( OPT_ZPRIOR == PSNID_ZPRIOR_HOST ) { 
    z    = DATA_PSNID_DOFIT.REDSHIFT[1] ;
    zerr = DATA_PSNID_DOFIT.REDSHIFT_ERR[1] ;
    sprintf(text_zprior,"%6.4f +- %6.4f", z, zerr );
  }
  else {
    sprintf(c1err,"Unknown OPT_ZPRIOR = %d ", OPT_ZPRIOR );
    sprintf(c2err,"Check  PSNID_ZPRIOR_XXX parameters.") ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }

  LDMP = (itype == -8080 ) ;  // RK - OCt 23 2013
  if ( LDMP ) {
    printf(" xxx -------------------------------- \n");
    printf(" xxx CCID=%s  itype = %d  itype_best=%d   izprior=%d\n", 
	   CCID, itype, itype_best, izprior); fflush(stdout);
  }

  fflush(stdout);

  if ( itype == itype_best ) 
    { sprintf(text_bestfit,"YES"); }
  else
    { sprintf(text_bestfit,"NO"); }


  sprintf(DISPLAYTEXT,"FUN = %s (BestFit = %s)     z-prior = %s",
	  PSNID_TYPE_NAME[itype], text_bestfit, text_zprior );
  SNLCPAK_DISPLAYTEXT(CCID, DISPLAYTEXT);


  // next string: color, colorlaw, shape
  // NOTE: sntools_grid.c need to read name of column: SHAPE -> DELTA, etc ..
  int  this_type ;  // which GRID file, 1 or 2
  int  ISHAPE, INDX_SIM, INDX_SPARSE, LEN ;
  char text_shape[40], parName[40], *NAME  ;
  double SHAPE ;

  if   ( itype == PSNID_ITYPE_SNIA ) { this_type = 1; }
  else                               { this_type = 2 ; }


  if( itype == PSNID_ITYPE_SNIA )  { 
    SHAPE = 
      PSNID_BEST_RESULTS.FINAL_PARVAL[izprior][itype][PSNID_PARAM_SHAPEPAR] ;
    sprintf(parName,"%s", 
	    SNGRID_PSNID[this_type].NAME[IPAR_GRIDGEN_SHAPEPAR] ) ;
    sprintf(text_shape, "%s = %6.3f", parName, SHAPE ); 
  }
  else { 
    ISHAPE = PSNID_BEST_RESULTS.FINAL_NONIA_INDX_SPARSE[itype];
    INDX_SPARSE = PSNID_NONIA_ABSINDEX[itype][ISHAPE]; 

    // get INDX_SIM from sim-input file to identify original template
    INDX_SIM    = SNGRID_PSNID[this_type].NON1A_INDEX[INDX_SPARSE] ;

    // get name of original template from sim-input file.
    NAME = SNGRID_PSNID[this_type].NON1A_NAME[INDX_SPARSE] ;
    LEN  = strlen(NAME);

    if ( LEN >= 20 ) 
      { sprintf(parName,"%s", &NAME[LEN-20] ) ; }
    else
      { sprintf(parName,"%s", NAME) ; }
      
    sprintf(text_shape, "NONIA-INDEX=%d(%s)", INDX_SIM, parName);


    /*
    printf(" xxx NONIA: itype=%d ILUM=%d  SIM_INDEX=%d (%s)\n",
    itype, ILUM, INDX_SIM, parName );  */
  }
      

  sprintf(DISPLAYTEXT,"%s=%6.3f     %s=%5.2f     %s "
	  ,SNGRID_PSNID[this_type].NAME[IPAR_GRIDGEN_COLORPAR]
	  ,PSNID_BEST_RESULTS.FINAL_PARVAL[izprior][itype][PSNID_PARAM_COLORPAR]
	  ,SNGRID_PSNID[this_type].NAME[IPAR_GRIDGEN_COLORLAW]
	  ,PSNID_BEST_RESULTS.FINAL_PARVAL[izprior][itype][PSNID_PARAM_COLORLAW]
	  ,text_shape
	  );
  SNLCPAK_DISPLAYTEXT(CCID, DISPLAYTEXT);


  // next string: MU, redshift, chi2/dof
  sprintf(DISPLAYTEXT,"dMU = %6.3f     z=%.4f     chi2/dof = %.1f/%d"
	  ,PSNID_BEST_RESULTS.FINAL_PARVAL[izprior][itype][PSNID_PARAM_DMU]
	  ,PSNID_BEST_RESULTS.FINAL_PARVAL[izprior][itype][PSNID_PARAM_LOGZ]
	  ,PSNID_BEST_RESULTS.MINCHISQ[izprior][itype]
	  ,PSNID_BEST_RESULTS.NGOOD[izprior][itype]
	  );
  SNLCPAK_DISPLAYTEXT(CCID, DISPLAYTEXT);

  // -------------------------------------------------
  // get pointers to light curve data and fit resids

  CCID          = DATA_PSNID_DOFIT.CCID ;
  NOBS          = DATA_PSNID_DOFIT.NOBS ;
  ptrMJD        = DATA_PSNID_DOFIT.MJD ;
  ptrFLUXDATA   = DATA_PSNID_DOFIT.FLUXDATA ;
  ptrFLUXERR    = DATA_PSNID_DOFIT.FLUXERR ;
  ptrFLUXSIM    = DATA_PSNID_DOFIT.FLUXSIM ;
  ptrIFILTOBS   = DATA_PSNID_DOFIT.IFILTOBS ;
  ptrDUMERR0    = DATA_PSNID_DOFIT.DUMERR0 ;

  ptrTOBS     = RESIDS_PSNID_DOFIT[itype].TOBS ;
  ptrREJ      = RESIDS_PSNID_DOFIT[itype].REJECT ;
  ptrCHI2     = RESIDS_PSNID_DOFIT[itype].CHI2 ;

  SNLCPAK_DATA(CCID, NOBS, ptrMJD,ptrTOBS, ptrFLUXDATA, ptrFLUXERR, ptrIFILTOBS,
	       SNLCPAK_EPFLAG_FLUXDATA  ) ;

  SNLCPAK_DATA(CCID, NOBS, ptrMJD, ptrTOBS, ptrREJ,  ptrDUMERR0, ptrIFILTOBS,
	       SNLCPAK_EPFLAG_REJECT  ) ;

  SNLCPAK_DATA(CCID, NOBS, ptrMJD, ptrTOBS, ptrCHI2, ptrDUMERR0, ptrIFILTOBS,
	       SNLCPAK_EPFLAG_CHI2  ) ;

  if ( PSNID_INPUTS.LSIM ) {
    // pass true sim fluxes, Jan 2020
    SNLCPAK_DATA(CCID, NOBS, ptrMJD,ptrTOBS,ptrFLUXSIM,ptrFLUXERR,ptrIFILTOBS,
	       SNLCPAK_EPFLAG_FLUXSIM  ) ;
  }

  // pak filter-dependent quantities
  IFILTLIST = PSNID_INPUTS.IFILTLIST ;

   
  SNLCPAK_DATA(CCID, NFILT, ptrMJD, ptrTOBS, 
	       RESIDS_PSNID_DOFIT[itype].NDOF_FILT, ptrDUMERR0, 
	       IFILTLIST, SNLCPAK_BANDFLAG_NDOF) ;

  SNLCPAK_DATA(CCID, NFILT, ptrMJD, ptrTOBS, 
	       RESIDS_PSNID_DOFIT[itype].CHI2_FILT, ptrDUMERR0, 
	       IFILTLIST, SNLCPAK_BANDFLAG_CHI2) ;

  SNLCPAK_DATA(CCID, NFILT, ptrMJD, ptrTOBS, 
	       RESIDS_PSNID_DOFIT[itype].PKMJD_FILT, ptrDUMERR0, 
	       IFILTLIST, SNLCPAK_BANDFLAG_PKMJD) ;

  // -------------------------------------------
  // pak best-fit function in 2 day bins
  // to have a smooth curve to overlay on the data.
  
  int     NBT, i, MEM_D, MEM_I ;
  double  TMIN = -60.0, TMAX = 140.0, TBIN = 2.0, xi ;

  // Aug 2017: use template grid if its Trest GRID is less than 1 day
  if ( PSNID_TBIN < 1.0 )
    { TBIN = PSNID_TBIN ;  TMIN/=2.0;  TMAX/=2.0; }
  
  NFILT_USE = RESIDS_PSNID_DOFIT[itype].NFILT_USE ;
  PKMJD = PSNID_BEST_RESULTS.FINAL_PARVAL[izprior][itype][PSNID_PARAM_TMAX];

  NBT = (int)((TMAX - TMIN + 1.0E-6)/TBIN) ;
  MEM_D = NBT * NFILT_USE * sizeof(double);
  MEM_I = NBT * NFILT_USE * sizeof(int);
  

  // allocate
  FITFUN_PSNID_DOFIT.MJD       = (double*)malloc(MEM_D) ;
  FITFUN_PSNID_DOFIT.TOBS      = (double*)malloc(MEM_D) ;
  FITFUN_PSNID_DOFIT.FLUX      = (double*)malloc(MEM_D) ;
  FITFUN_PSNID_DOFIT.FLUX_ERR  = (double*)malloc(MEM_D) ;
  FITFUN_PSNID_DOFIT.IFILT     = (int   *)malloc(MEM_I) ;
  FITFUN_PSNID_DOFIT.IFILTOBS  = (int   *)malloc(MEM_I) ;

  // set local pointers
  ptrMJD      = FITFUN_PSNID_DOFIT.MJD ;
  ptrTOBS     = FITFUN_PSNID_DOFIT.TOBS ;
  ptrFLUXDATA = FITFUN_PSNID_DOFIT.FLUX ;
  ptrFLUXERR  = FITFUN_PSNID_DOFIT.FLUX_ERR ;
  ptrIFILT    = FITFUN_PSNID_DOFIT.IFILT ;
  ptrIFILTOBS = FITFUN_PSNID_DOFIT.IFILTOBS ;

  NOBS = 0;
  for ( IFILT=0; IFILT < NFILT; IFILT++ ) {
    
    USE = RESIDS_PSNID_DOFIT[itype].USE_FILT[IFILT] ;
    if ( USE == 0 ) { continue ; }

    IFILTOBS  = PSNID_INPUTS.IFILTLIST[IFILT]; // absolute index

    for(i=0; i < NBT; i++ ) {
      xi = (double)i + 0.5 ;
      TOBS = TMIN + (TBIN * xi) ;
      MJD  = TOBS + PKMJD ;
      
      // call generic GET_FITFUN for this TOBS
      PSNID_BEST_GET_FITFUN(CCID, iplot, 
			  IFILT, ONE,          // (I) ifilt, Nepo
			  &TOBS       ,        // (I) Tobs
			  &MAG, &MAG_ERR ) ;   // (O) best-fit mag and error
      
      // convert MAG and MAG_ERR into FLUX and FLUX_ERR
      // Note: this function should be promoted to psnid_tools.c
      psnid_pogson2fluxcal(MAG, MAG_ERR, &FITFLUX, &FITFLUX_ERR);
    
      // load arrays
      ptrMJD[NOBS]        = MJD ;
      ptrTOBS[NOBS]       = TOBS ;
      ptrFLUXDATA[NOBS]   = FITFLUX ;
      ptrFLUXERR[NOBS]    = FITFLUX_ERR ;
      ptrIFILTOBS[NOBS]   = IFILTOBS ;
      ptrIFILT[NOBS]      = IFILT ;
      NOBS++ ;
      
    }  // i (Tobs loop)
  } // IFILT
    

  SNLCPAK_DATA(CCID, NOBS, ptrMJD, ptrTOBS, ptrFLUXDATA, ptrFLUXERR, 
	       ptrIFILTOBS, SNLCPAK_EPFLAG_FITFUN  ) ;

  // free
  free(FITFUN_PSNID_DOFIT.MJD) ;
  free(FITFUN_PSNID_DOFIT.TOBS) ;
  free(FITFUN_PSNID_DOFIT.FLUX) ;
  free(FITFUN_PSNID_DOFIT.FLUX_ERR);
  free(FITFUN_PSNID_DOFIT.IFILT);
  free(FITFUN_PSNID_DOFIT.IFILTOBS);

  // -------------------------------------------
  // fill output
  SNLCPAK_FILL(CCID) ;

  return ;
  
} // end of  SNLCPLOT_PSNID_BEST

void snlcplot_psnid_best__(int *iplot) {
  SNLCPLOT_PSNID_BEST(*iplot) ;
}


/*************************************************************/
void MONPLOT_PSNID_BEST(int iplot) {

  // create arbitrary 1D and 2D plots.
  // These plots are stored in SN-specifid subdirs named
  // SN[CCID]_FIT[iplot]

  int    NDIM, NB[2], HID, i0, i1 ;
  double XMIN[2], XMAX[2], X[2], WGT ;
  char   CHTIT[80];
  //  char   fnam[] = "MONPLOT_PSNID_BEST" ;

  // ----------- BEGIN ----------

  NDIM = 1 ;  // 1D hist
  HID  = 801 ; // pick arbitrary number
  sprintf(CHTIT,"Example 1D histogram");
  NB[0] = 10; XMIN[0]=0.0; XMAX[0]=10.0 ;
  SNHIST_INIT(NDIM, HID, CHTIT, NB, XMIN, XMAX ) ;

  // now fill histogram
  for(i0=0; i0<NB[0]; i0++ ) {
    X[0] = (double)i0 + 0.5 ; WGT = X[0];
    SNHIST_FILL(NDIM, HID, X, WGT );
  }


  // now make and fill a 2D plot
  NDIM = 2 ;  // 2D hist
  HID  = 802 ; // pick arbitrary number
  sprintf(CHTIT,"Example 2D histogram");
  NB[0] = 10; XMIN[0]=0.0; XMAX[0]=10.0 ;
  NB[1] = 20; XMIN[1]=0.0; XMAX[1]=20.0 ;
  SNHIST_INIT(NDIM, HID, CHTIT, NB, XMIN, XMAX ) ;

  // now fill histogram
  for(i0=0; i0<NB[0]; i0++ ) {
    for(i1=0; i1<NB[1]; i1++ ) {
      X[0] = (double)i0  + 0.5 ;  // horizont axis
      X[1] = (double)i1  + 0.5 ;  // vertical axis
      WGT = X[0] * X[1];
      SNHIST_FILL(NDIM, HID, X, WGT );
    }
  }

  return ;
  
} // end of MONPLOT_PSNID_BEST

void monplot_psnid_best__(int *iplot) {
  MONPLOT_PSNID_BEST(*iplot) ;
}


/**************************************************************/

/**********************************************************************/
int PSNID_BEST_GET_NFIT(void)
/**********************************************************************/
{
  int n;
  //  n = PSNID_NTYPES*PSNID_BEST_RESULTS.NZPRIOR ;
  n = PSNID_NTYPES ;  // RK - Mar 9 2013
  return(n);
}
// end of  PSNID_BEST_GET_NFIT


/**********************************************************************/
void psnid_best_get_plotind(int ind, int *this_z, int *this_t)
/**********************************************************************/
{
  int i, z, count = 0;
  int DO ;

  *this_z = -9 ;
  *this_t = -9 ;

  for (z=0; z<PSNID_NZPRIOR; z++) {
    DO = PSNID_BEST_RESULTS.ZPRIOR_DO[z] ;

    if (DO == 1) {
      for (i=0; i<PSNID_NTYPES; i++) {
	if (ind == count) {
	  *this_z = z;
	  *this_t = i;
	}
	count++;
      }
    }
    
  }

  return;
}
// end of psnid_best_get_plotind


// ***********************************
int doPlot_psnid(int iplot) {

  // return 1 to make this plot; 0 to skip plot.

  int izprior, itype;
  psnid_best_get_plotind(iplot, &izprior, &itype);
	 
  if ( PSNID_NGRID[itype] > 0 )
    { return 1 ; }
  else
    { return 0 ; }
  
} // end doPlot_psnid

int doplot_psnid__(int *iplot) { return doPlot_psnid(*iplot); }

/**********************************************************************/
void PSNID_BEST_GET_FITPAR(char *CCID, int ind,
		      double *redshift, double *shapepar, double *colorpar,
		      double *colorlaw, double *pkmjd, double *dmu,
		      double *chi2, int *ndof,
		      char *plotlabel, char *funLabel)
/*
  Returns PSNID best-fit parameters for a given SN type and redshift prior.

  Input:

    CCID    = candidate ID (not yet used for anything)
    itype   = SN type (see PSNID_ITYPE_SNIA/IBC/II)
    izprior = redshift prior (see PSNID_ZPRIOR_FLAT/SNLC/HOST)

  Output:

    redshift = redshift
    shapepar = x1,Delta or  CC index
    colorpar = c or AV
    colorlaw = RV (not fit)
    pkmjd    = peak MJD
    dmu      = delta mu relative to templates

 */
/**********************************************************************/
{
  int itype, izprior;

  psnid_best_get_plotind(ind, &izprior, &itype);

  *redshift = PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_LOGZ];
  *shapepar = PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_SHAPEPAR];
  *colorpar = PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_COLORPAR];
  *colorlaw = PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_COLORLAW];
  *pkmjd    = PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_TMAX];
  *dmu      = PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_DMU];

  *chi2     = PSNID_BEST_RESULTS.MINCHISQ[izprior][itype];
  *ndof     = PSNID_BEST_RESULTS.NGOOD[izprior][itype];


  // RK also return  errors on fit params
  *(redshift+1) = 0.001 ;
  *(shapepar+1) = 0.001 ;
  *(colorpar+1) = 0.001 ;
  *(colorlaw+1) = 0.001 ;
  *(pkmjd+1)    = 0.001 ;
  *(dmu+1)      = 0.001 ;


  // RK - slightly informative comment for debugging.
  // RK - note to Masao: need text strings for izprior.

  sprintf(plotlabel,"CCID=%s  ind=%d  izprior=%d", 
	  CCID, ind,  izprior);
  sprintf(funLabel, "%s", PSNID_TYPE_NAME[itype]);

  /*
  printf("\t z = %8.4f  shapepar = %5.2f  colorpar = %5.2f"
	 "  colorlaw = %5.2f  pkmjd = %10.2f  dmu = %5.2f"
	 "  chi2 = %10.2e  ndof = %3d  label = %s\n",
	 *redshift, *shapepar, *colorpar,
	 *colorlaw, *pkmjd, *dmu,
	 *chi2, *ndof,
	 plotlabel);
  */

  return;
}
// end of PSNID_BEST_GET_FITPAR



/**********************************************************************/
void PSNID_BEST_GET_FITFUN(char *CCID, int ind,
			   int ifilt, int nepoch, double *epoch,
			   double *mag, double *magerr)
/*
  Returns PSNID best-fit model for a given SN type, redshift prior, 
  epoch array, and filter.  Linear interpolation between grid points.

  Input:

    CCID    = candidate ID (not yet used for anything)
    itype   = SN type (see PSNID_ITYPE_SNIA/IBC/II)
    ifilt   = filter index
    nepoch  = number of epochs to return fit-function
    *epoch  = MJD (>10000) vector or relative to peak (TOBS)

  Output:
    mag     = magnitude
    magerr  = magnitude error


 Sep 9 2013: define dtDT and dcDC to avoid redundant computations.

 Oct 23, 2013 RK - fix aweful bug: first argument to  psnid_best_modelErr
                   is "this_type" instead of 't'. Affected chi2 output
                   in LDMP_SNFLUX and also for plotting light curves.
                   Does NOT affect fitres/ascii dump of psnid results.

 Aug 17, 2014 RK - fix awful bug computing redshift with wrong
                   array. Was messing up light curve plots, 
                   but not the fits.

 Sep 10 2015 RK - fix bug calling hunt() function; need to pass 
                  fortran-like array starting at 1 instead of 0.
                  See new array trest1_forHunt that starts at 1.

 Oct 10 2015 RK - fix memory leak by adding free(trest1_forHunt) 
 
 Feb 18 2017 RK - fix another mem leak about returning BEFORE mallocs.

 **********************************************************************/
{
  int i;
  int z,l,t,f, this_t=0, this_type=0, this_filt=0, this_l=0;
  double color, ushift, this_epoch=0.0;
  double magatmp, magbtmp, dtDT, dcDC, gridErr;
  int ind1c1, ind2c1;  
  double *trest1, *mag1, *magerr1, *trest1_forHunt;
  double *trest2, *mag2, *magerr2;
  double color1, color2, redshift, TMAX;
  int nepoch1=0, nepoch2=0, optDump=0;
  int itype, izprior;
  //  char fnam[] = "PSNID_BEST_GET_FITFUN";

  // ------------------- BEGIN -------------------

    // init output arrays
  for(i=0; i < nepoch ; i++ )   { mag[i] = magerr[i] = 99.0 ; }

  
  psnid_best_get_plotind(ind, &izprior, &itype); // return izprior and itype

  // bail if there are no templates for this type (Feb 2017)
  // Make sure to bail BEFORE malloc.
  if ( PSNID_NGRID[itype] == 0 ) { return ; } 

  // allocate temp memory
  int MEM = sizeof(double) * MXEP_PSNID ;
  trest1    = (double*) malloc(MEM);
  trest2    = (double*) malloc(MEM);
  mag1      = (double*) malloc(MEM);
  mag2      = (double*) malloc(MEM);
  magerr1   = (double*) malloc(MEM);
  magerr2   = (double*) malloc(MEM);
  trest1_forHunt = (double*) malloc(MEM);


  // init local arrays.
  for (i=0; i < MXEP_PSNID ; i++) {
    trest1[i]  = -9.9;
    mag1[i]    = -9.9;
    magerr1[i] = -9.9;
    trest2[i]  = -9.9;
    mag2[i]    = -9.9;
    magerr2[i] = -9.9;
    trest1_forHunt[i] = -9.9 ;
  }


  if (itype == PSNID_ITYPE_SNIA) {
    this_type = TYPEINDX_SNIA_PSNID;
  }
  else if ( itype > PSNID_ITYPE_SNIA && itype <= PSNID_NTYPES ) {
    this_type = TYPEINDX_NONIA_PSNID;
  }


  f = ifilt;  // filter index (0..N-1)
  this_filt = PSNID_INPUTS.IFILTLIST[f];  // absolute filter index

  l = PSNID_BEST_RESULTS.MINCHISQ_IND[izprior][itype][PSNID_PARAM_SHAPEPAR];
  z = PSNID_BEST_RESULTS.MINCHISQ_IND[izprior][itype][PSNID_PARAM_LOGZ];
  if (itype == PSNID_ITYPE_SNIA) {
    this_l = l;
  }
  else if ( itype > PSNID_ITYPE_SNIA && itype <= PSNID_NTYPES ) {
    this_l = PSNID_NONIA_ABSINDEX[itype][l];
  }
  //  printf("  l=%d\n", l);

  redshift = 
    PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_LOGZ];

  color  = 
    PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_COLORPAR];

  ushift = 
    PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_DMU];

  TMAX = 
    PSNID_BEST_RESULTS.GRID_FINAL_PAR[izprior][itype][PSNID_PARAM_TMAX];

  // get models of first and last color bins
  // mags from first color bin
  ind1c1 = 1;
  ind2c1 = SNGRID_PSNID[this_type].NBIN[IPAR_GRIDGEN_COLORPAR];
  
  color1 = SNGRID_PSNID[this_type].VALUE[IPAR_GRIDGEN_COLORPAR][ind1c1];
  color2 = SNGRID_PSNID[this_type].VALUE[IPAR_GRIDGEN_COLORPAR][ind2c1];
  //  printf("  color =%8.3f\n", color);
  //  printf("  color1=%8.3f\n", color1);
  //  printf("  color2=%8.3f\n", color2);

  /*
  printf("\t this_type=%d  z=%d  ind1c1=%d  this_l=%d  this_filt=%d  "
	 "l=%d f=%d\n",
	 this_type, z, ind1c1, this_l, this_filt, l, f);
  */

  // bail if filter not used
  if (PSNID_INPUTS.USEFILT[this_filt] != 1) { goto FREE_MEM ; }

  //                         z      c1 c2      lu          f
  get_template_lc(this_type, z, ind1c1, 1, this_l, this_filt, optDump,
		  &nepoch1, trest1, mag1, magerr1);
  get_template_lc(this_type, z, ind2c1, 1, this_l, this_filt, optDump,
		  &nepoch2, trest2, mag2, magerr2);
  
  if (nepoch1 != nepoch2) { goto FREE_MEM ; } // sanity check
   

  // Sep 10 2015: create fortran-like array for hunt function
  for(i=0; i < nepoch1; i++ ) { trest1_forHunt[i+1] = trest1[i];  }

  // ------------------------------
  for (i=0; i < nepoch; i++) {	
	
    if (epoch[i] > 10000.) {
      this_epoch = (epoch[i]-TMAX)/(1.+redshift) ;
    } else {
      this_epoch = epoch[i]/(1.+redshift);
    }
    

    //    hunt(trest1, nepoch1, this_epoch, &this_t); // return this_t
    hunt(trest1_forHunt, nepoch1, this_epoch, &this_t); // return this_t
    
    if (this_t < 1 || this_t > nepoch1-2) {     
      mag[i]    = 99.00;
      magerr[i] = 99.00;
    } 
    else {
      
      t = this_t;
	  
      dcDC = (color-color1) / (color2-color1) ;
      magatmp = mag1[t]   - dcDC * (mag1[t]   - mag2[t]  ) + ushift;
      magbtmp = mag1[t+1] - dcDC * (mag1[t+1] - mag2[t+1]) + ushift;
      
      dtDT = (this_epoch - trest1[t]) / (trest1[t+1] - trest1[t]) ;

      if ( magbtmp > 31.0 )
	{ mag[i] = magbtmp; } 
      else
	{ mag[i]  = magatmp  + ( magbtmp   - magatmp    ) * dtDT ; }
      
      mag[i] += psnid_best_mwxtmag(ifilt);

      gridErr = magerr1[t] + ( magerr1[t+1] - magerr1[t] ) * dtDT ;
      
      // RK 9.07.2013
      magerr[i] = psnid_best_modelErr( this_type, this_epoch, gridErr ); 
       
      // xxxxxx
      if ( nepoch == 1  && ifilt == -3  &&
	   fabs(epoch[0]-9.2) < 0.5 ) {  
	printf(" xxx ------ call modelErr (t=%d) ---------- \n", t);
	printf(" xxx nepoch1=%d  nepoch2=%d \n", nepoch1, nepoch2 );
	printf(" xxx this_type=%d  TMAX=%.3f  this_ep=%.2f  \n",
	       this_type, TMAX, this_epoch);
	printf(" xxx magatmp=%.2f  magbtmp=%.2f    mwxtmag=%.2f \n",
	       magatmp, magbtmp, psnid_best_mwxtmag(ifilt) );
	printf(" xxx mag1[t+1]=%.3f  mag2[t+1]=%.3f \n",
	       mag1[t+1], mag2[t+1] );	  
	printf(" xxx dtDT=%.3f  dcDC=%.3f  ushift=%.3f \n", dtDT, dcDC, ushift );
	printf(" xxx mag=%.3f  gridErr = %.3f  magerr=%.3f \n",
	       mag[i], gridErr, magerr[i] );
	fflush(stdout);
      }
      // xxxxxxxx
      
    }
  }  // end i loop over nepoch
  
  // free temp memory
 FREE_MEM:
  free(trest1) ;
  free(trest2) ;
  free(mag1);
  free(mag2);
  free(magerr1);
  free(magerr2);
  free(trest1_forHunt);

  return;

}   // end of psnid_best_get_fitfun



// ***********************************************
int psnid_bestType_cuts(int z) {

  // Created Feb 13 2017 by R.Kessler
  // Apply cuts and return best type
  // [move fragile code from psnid_best_store_PBayes and 
  //   replace with loops & arrays]
  //
  // Inputs: z = zPrior index

  int ITYPE_BEST_PBAYES=-9, ITYPE_BEST_CUTS = -9 ;
  int i, npt, LDMP=0 ;
  double chi2red, chi2, PBAYES, FITPROB, PBAYES_MAX=-1.0;
  //  char fnam[] = "psnid_bestType_cuts" ;

  // 816 895 are now IA, but were UNKNOWN before

  // --------------- BEGIN -----------------

  if ( LDMP ) {  printf("\n xxx ----------------------------- \n" ); }

  for(i=0; i < PSNID_NTYPES; i++ ) {
    if ( PSNID_NGRID[i] == 0 ) { continue ; }
    PBAYES = PSNID_BEST_RESULTS.PBAYESIAN[z][i] ;
    if ( PBAYES > PBAYES_MAX ) 
      { PBAYES_MAX = PBAYES; ITYPE_BEST_PBAYES=i; }

    if ( LDMP ) { printf(" xxx PBAYES[%d] = %f \n", i, PBAYES); }
  }

  if ( LDMP ) {
    printf(" xxx ITYPE_BEST_PBAYES = %d  PBAYES_MAX=%f\n", 
	   ITYPE_BEST_PBAYES, PBAYES_MAX ); fflush(stdout);
  }

  if ( ITYPE_BEST_PBAYES < 0 ) { return(ITYPE_BEST_PBAYES); }

  // Apply user cuts 

  PBAYES  = PBAYES_MAX ;
  FITPROB = PSNID_BEST_RESULTS.FITPROB[z][ITYPE_BEST_PBAYES];
  npt     = PSNID_BEST_RESULTS.NGOOD[z][ITYPE_BEST_PBAYES];
  chi2    = PSNID_BEST_RESULTS.MINCHISQ[z][ITYPE_BEST_PBAYES];
  chi2red = 9999.999 ;
  if ( npt > 4 ) { chi2red = chi2/(double)(npt-4); }


  if ( LDMP ) {
    printf(" xxx FITPROB = %f (FITPROB_CUT=%f)\n", 
	   FITPROB,PSNID_FITPROB_CUTLIST[ITYPE_BEST_PBAYES] ) ;
    printf(" xxx chi2red = %f / (%d-4) = %f \n",
	   chi2, npt, chi2red); 
    fflush(stdout);
  }

  if ( PBAYES < PSNID_PBAYES_CUTLIST[ITYPE_BEST_PBAYES] ) 
    { return(ITYPE_BEST_CUTS); }

  if ( FITPROB < PSNID_FITPROB_CUTLIST[ITYPE_BEST_PBAYES] ) 
    { return(ITYPE_BEST_CUTS); }

  if ( (ITYPE_BEST_PBAYES == PSNID_ITYPE_SNIA) && chi2red > 20.0 ) 
    { return(ITYPE_BEST_CUTS); }


  // if we get here, all cuts pass so return best ITYPE
  ITYPE_BEST_CUTS = ITYPE_BEST_PBAYES ;
  return(ITYPE_BEST_CUTS); 

} // end psnid_bestType_cuts

/**********************************************************************/
int PSNID_BEST_GET_BESTTYPE(void)
/**********************************************************************/
{
  // return already-determined best type to external program.
  int type;
  type = PSNID_BEST_RESULTS.FINAL_ITYPE[0];
  return type;
}


/**********************************************************************/
void PSNID_BEST_GET_TYPENAME(int itype, char *name) 
/**********************************************************************/
{
  // return outut *name associated with input itype.
  sprintf(name,"%s", PSNID_TYPE_NAME[itype]);
}

/**********************************************************************/
void psnid_best_print_zprior_info(int i)
/**********************************************************************/
{

  printf("\n\t %-10s  ", PSNID_BEST_RESULTS.CID);

  if (i>0) {
    printf("using z = %8.5f +/- %8.5f\n",
	   PSNID_BEST_RESULTS.ZPRIOR[i], PSNID_BEST_RESULTS.ZPRIOR_ERR[i]);
  } else {
    printf("using flat redshift prior\n");
  }

  return;
}
// end of psnid_best_print_zprior_info


/**********************************************************************/
void psnid_best_store_results(char *CCID, int itype, int z,
			      int ***ind, double **evidence, int optdump)
/**********************************************************************/
{
  int i ;
  double redshift, shapepar, colorpar, colorlaw, pkmjd, dmu;

  // ------------- BEGIN ----------------
  
  if ( PSNID_NGRID[itype] == 0 ) { return ; }  // Aug 2017
  
  // determine best-fit parameters from indices
  psnid_best_calc_fitpar(itype, ind,
			 &redshift, &shapepar, &colorpar, &colorlaw,
			 &pkmjd, &dmu);


  // fit probability
  PSNID_BEST_RESULTS.FITPROB[z][itype] = 
    PROB_Chi2Ndof(PSNID_BEST_RESULTS.MINCHISQ[z][itype], 
		  PSNID_BEST_RESULTS.NGOOD[z][itype]);

  // indices of best-fit model
  for (i=0; i<PSNID_NPARAM; i++) {
    PSNID_BEST_RESULTS.MINCHISQ_IND[z][itype][i] = ind[PSNID_NITER][itype][i];
  }

  PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_LOGZ]     = redshift;
  PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_SHAPEPAR] = shapepar;
  PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_COLORPAR] = colorpar;
  PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_COLORLAW] = colorlaw;
  PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_TMAX]     = pkmjd;
  PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_DMU]      = dmu;


  if (optdump == 1) {

    printf("\t %-10s", CCID);

    // redshift
    printf("  z[%03d] = %8.5f",
	   ind[PSNID_NITER][itype][PSNID_PARAM_LOGZ],
	   PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_LOGZ]);

    // shapepar
    printf("   %5s[%02d] = %5.2f",
	   SNGRID_PSNID[PSNID_THIS_TYPE].NAME[IPAR_GRIDGEN_SHAPEPAR],
	   ind[PSNID_NITER][itype][PSNID_PARAM_SHAPEPAR],
	   PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_SHAPEPAR]);

    // colorpar
    printf("   %2s[%02d] = %5.2f",
	   SNGRID_PSNID[PSNID_THIS_TYPE].NAME[IPAR_GRIDGEN_COLORPAR],
	   ind[PSNID_NITER][itype][PSNID_PARAM_COLORPAR],
	   PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_COLORPAR]);

    // colorlaw
    printf("   %4s[%02d] = %5.2f",
	   SNGRID_PSNID[PSNID_THIS_TYPE].NAME[IPAR_GRIDGEN_COLORLAW],
	   1,
	   PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_COLORLAW]);

    // pkmjd
    printf("   Tmax[%02d] = %10.2f",
	   ind[PSNID_NITER][itype][PSNID_PARAM_TMAX],
	   PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_TMAX]);

    // dmu
    printf("   dmu[%02d] = %5.2f",
	   ind[PSNID_NITER][itype][PSNID_PARAM_DMU],
	   PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][itype][PSNID_PARAM_DMU]);

    // chisq, npt, Pfit, PBayes, lcq
    printf("   chi2 = %10.3e  E = %10.3e  npt = %3d  lcq = %2d\n",
	   PSNID_BEST_RESULTS.MINCHISQ[z][itype],
	   evidence[z][itype], PSNID_BEST_RESULTS.NGOOD[z][itype],
	   PSNID_BEST_RESULTS.LIGHTCURVE_QUALITY[z][itype]);
  }


  return;
}   // end of psnid_best_store_results


/**********************************************************************/
void psnid_best_store_PBayes(char *CCID, int z, double **evidence, int optdump)
/**********************************************************************/
{

  // Oct 24 2013 RK - change hard-wired chir cut from 2 to 20
  //                  to allow loser FITPROB_IA_CUT. No need for CHIR cut.

  int i ;
  double sum_evidence=0.0;

  // -------------- BEGIN --------------

  // init PBAYES
  for (i=0; i<PSNID_NTYPES; i++) 
    {  PSNID_BEST_RESULTS.PBAYESIAN[z][i] = -1.0; }

  for (i=0; i<PSNID_NTYPES; i++) {
    if ( PSNID_NGRID[i] > 0 ) { sum_evidence += evidence[z][i]; }
  }

  if (sum_evidence > 0.0) {
    for (i=0; i<PSNID_NTYPES; i++) {
      if ( PSNID_NGRID[i] > 0 ) {
	PSNID_BEST_RESULTS.PBAYESIAN[z][i] = 
	  evidence[z][i]/sum_evidence;
      }
    }
  } 
 
  // get best itype with cuts using function (RK, Feb 2017)
  int itype_best = psnid_bestType_cuts(z);
  
  if ( itype_best >= 0 ) {
    sprintf(PSNID_BEST_RESULTS.FINAL_TYPE[z],"%s", 
	    PSNID_TYPE_NAME[itype_best]);
    PSNID_BEST_RESULTS.FINAL_ITYPE[z] = itype_best ;
  }
  else {
    sprintf(PSNID_BEST_RESULTS.FINAL_TYPE[z],"UNKNOWN");
    PSNID_BEST_RESULTS.FINAL_ITYPE[z] = -9;
  }


  return ;


}  // end of psnid_best_store_PBayes


/**********************************************************************/
void psnid_best_dump_results()
/**********************************************************************/
{
  int z, t, this_type=0;
  double Z, ZERR;
  char   *CID, cprior[60];

  // Nov  1 2013 - RK : for negative zerr, print '-> flat' message
  // Aug 27 2017 - RK : print only types used in fit.

  // -------------- BEGIN --------------
  
  CID = PSNID_BEST_RESULTS.CID ;

  for (z=0; z<PSNID_NZPRIOR; z++) {

    if (PSNID_BEST_RESULTS.ZPRIOR_DO[z] == 1) {
      Z    = PSNID_BEST_RESULTS.ZPRIOR[z] ;
      ZERR = PSNID_BEST_RESULTS.ZPRIOR_ERR[z] ;

      
      if ( ZERR < 0.0 ) 
	{ sprintf(cprior,"-> flat") ; }
      else              
	{ cprior[0] = 0  ; }

      printf("\t %-10s  z_prior = %8.5f +/- %8.5f  %s \n",
	     CID, Z, ZERR, cprior );

      for (t=0; t<PSNID_NTYPES; t++) {

	if( PSNID_NGRID[t] == 0 ) { continue ; }
	
	if (t==PSNID_ITYPE_SNIA) {
	  this_type = TYPEINDX_SNIA_PSNID;
	}
	else if ( t>PSNID_ITYPE_SNIA && t <= PSNID_NTYPES ) {
	  this_type = TYPEINDX_NONIA_PSNID;
	}

	printf("\t %-8s  %3s", CID, PSNID_TYPE_NAME[t]);

	printf("  z[%03d] = %8.5f",
	       PSNID_BEST_RESULTS.MINCHISQ_IND[z][t][PSNID_PARAM_LOGZ],
	       PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][t][PSNID_PARAM_LOGZ]);

	printf("  %5s[%02d] = %5.2f",
	       SNGRID_PSNID[this_type].NAME[IPAR_GRIDGEN_SHAPEPAR],
	       PSNID_BEST_RESULTS.MINCHISQ_IND[z][t][PSNID_PARAM_SHAPEPAR],
	       PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][t][PSNID_PARAM_SHAPEPAR]);

	printf("  %2s[%02d] = %5.2f",
	       SNGRID_PSNID[this_type].NAME[IPAR_GRIDGEN_COLORPAR],
	       PSNID_BEST_RESULTS.MINCHISQ_IND[z][t][PSNID_PARAM_COLORPAR],
	       PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][t][PSNID_PARAM_COLORPAR]);

	printf("  %4s[%02d] = %5.2f",
	       SNGRID_PSNID[this_type].NAME[IPAR_GRIDGEN_COLORLAW],
	       1,
	       PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][t][PSNID_PARAM_COLORLAW]);

	printf("  Tmax[%02d] = %10.2f",
	       PSNID_BEST_RESULTS.MINCHISQ_IND[z][t][PSNID_PARAM_TMAX],
	       PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][t][PSNID_PARAM_TMAX]);

	printf("  dmu[%02d] = %5.2f",
	       PSNID_BEST_RESULTS.MINCHISQ_IND[z][t][PSNID_PARAM_DMU],
	       PSNID_BEST_RESULTS.GRID_FINAL_PAR[z][t][PSNID_PARAM_DMU]);

	printf("  chi2 = %9.3e  npt = %3d  prob = %9.2e  PBayes = %7.4f  lcq = %2d",
	       PSNID_BEST_RESULTS.MINCHISQ[z][t],
	       PSNID_BEST_RESULTS.NGOOD[z][t],
	       PSNID_BEST_RESULTS.FITPROB[z][t],
	       PSNID_BEST_RESULTS.PBAYESIAN[z][t],
	       PSNID_BEST_RESULTS.LIGHTCURVE_QUALITY[z][t]);

	printf("\n");

      }

      if (MCMC_RUN==1) {  // SN Ia only

	this_type = TYPEINDX_SNIA_PSNID;

	//	printf("\t %-10s  MCMC results\n", PSNID_BEST_RESULTS.CID);
	printf("\t %-10s  MCMC     z =  %9.4f ( %+9.4f %+9.4f )\n",
	       CID,
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_LOGZ][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_LOGZ][3]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_LOGZ][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_LOGZ][1]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_LOGZ][2]);
	printf("\t %-10s  MCMC %5s =   %8.3f ( %+8.3f  %+8.3f  )\n",
	       CID,
	       SNGRID_PSNID[this_type].NAME[IPAR_GRIDGEN_SHAPEPAR],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_SHAPEPAR][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_SHAPEPAR][3]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_SHAPEPAR][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_SHAPEPAR][1]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_SHAPEPAR][2]);
	printf("\t %-10s  MCMC %5s =   %8.3f ( %+8.3f  %+8.3f  )\n",
	       PSNID_BEST_RESULTS.CID,
	       SNGRID_PSNID[this_type].NAME[IPAR_GRIDGEN_COLORPAR],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_COLORPAR][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_COLORPAR][3]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_COLORPAR][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_COLORPAR][1]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_COLORPAR][2]);
	printf("\t %-10s  MCMC    mu =   %8.3f ( %+8.3f  %+8.3f  )\n",
	       PSNID_BEST_RESULTS.CID,
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_DMU][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_DMU][3]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_DMU][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_DMU][1]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_DMU][2]);
	printf("\t %-10s  MCMC  Tmax = %10.2f ( %+7.2f   %+7.2f   )\n",
	       PSNID_BEST_RESULTS.CID,
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_TMAX][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_TMAX][3]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_TMAX][2],
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_TMAX][1]-
	       PSNID_BEST_RESULTS.MCMC_FINAL_PAR[z][PSNID_ITYPE_SNIA][PSNID_PARAM_TMAX][2]);

      }

    }
  }


  z = 0 ;
  printf("\t %-10s  best type = %s (%d) \n"
	 ,PSNID_BEST_RESULTS.CID
	 ,PSNID_BEST_RESULTS.FINAL_TYPE[z]
	 ,PSNID_BEST_RESULTS.FINAL_ITYPE[z]
	 );  
  fflush(stdout);

  return;
}
// end of psnid_best_dump_results


/**********************************************************************/
void psnid_best_check_for_avprior()
/**********************************************************************/
/***
    Determine whether to use AV prior based on input AV_TAU, AV_SMEAR
    Default is AV_TAU   = -9.9
               AV_SMEAR = -9.9   (see psnid.car)
	       --> PSNID_USE_AV_PRIOR = 0
 ***/
{

  char cAV[] = "AV-extinction prior"; // RK: Mar 31 2013

  if (PSNID_INPUTS.AV_TAU>0.0 && PSNID_INPUTS.AV_SMEAR>0.0) {
    PSNID_USE_AV_PRIOR = 1;
    printf("\t %s is ON \n", cAV);  fflush(stdout);
  } else {
    PSNID_USE_AV_PRIOR = 0;
    printf("\t %s is OFF \n", cAV);  fflush(stdout);
  }
  
  return;
}


/**********************************************************************/
void psnid_best_check_for_dmu()
/**********************************************************************/
/***

 ***/

{

  char cDMU[] = "Fit for dMU:" ;

  if (PSNID_INPUTS.NDMU>1) {
    PSNID_FITDMU    = 1;
    PSNID_FITDMU_CC = 1;
    PSNID_UMIN  = PSNID_INPUTS.DMU_MIN;
    PSNID_USTEP = (PSNID_INPUTS.DMU_MAX-PSNID_INPUTS.DMU_MIN)/(PSNID_INPUTS.NDMU-1.0);
    PSNID_MAXNU = PSNID_INPUTS.NDMU;
    printf("\t %s YES \n" , cDMU ); fflush(stdout);
  } else {
    PSNID_FITDMU    = 0;
    PSNID_FITDMU_CC = 0;
    PSNID_UMIN      = 0.0;
    PSNID_USTEP     = 0.0;
    PSNID_MAXNU     = 1;
    printf("\t %s NO \n" , cDMU ); fflush(stdout);
  }


  return;
}



/**********************************************************************/
double psnid_best_avprior1(int itype, double av)
/**********************************************************************/
{
  double pav=0.0;

  if (itype == 0) {         // Ia
    if (av >= 0.0) {
      pav = exp(-av/PSNID_INPUTS.AV_TAU);
    } else {
      pav = exp(av/PSNID_INPUTS.AV_SMEAR);
    }
  } else if (itype == 1) {  // Ibc
    pav = 1.0;
  } else if (itype == 2) {  // II
    pav = 1.0;
  } else {                  // something else
    pav = 1.0;
  }

  return pav;
}
// end of psnid_best_avprior1


/**********************************************************************/
void psnid_best_get_c_grid(double *grid)
/**********************************************************************/
{
  int c;

  for (c = 1; c <= PSNID_MAXNA; c++)
    grid[c] = PSNID_AMIN + PSNID_ASTEP*(c-1);

  return;

}
// end of psnid_best_get_c_grid


/**********************************************************************/
void psnid_best_get_z_grid(double *grid)
/**********************************************************************/
{
  int z;

  for (z = 1; z <= PSNID_MAXNZ; z++)
    grid[z] = pow(10.,SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_LOGZ][z]);


  return;

}
// end of psnid_best_get_z_grid


/**********************************************************************/
void psnid_best_get_l_grid(double *grid)
/**********************************************************************/
{
  int d;

  for (d = 1; d <= PSNID_MAXNL; d++)
    grid[d] = SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_SHAPEPAR][d];

  return;

}
// end of psnid_best_get_l_grid


/**********************************************************************/
void psnid_best_get_u_grid(double *grid)
/**********************************************************************/
{
  int d;

  for (d = 1; d <= PSNID_MAXNU; d++)
    grid[d] = PSNID_UMIN + PSNID_USTEP*(d-1);

  return;

}
// end of psnid_best_get_u_grid


/**********************************************************************/
void psnid_best_init_mcmc()
/**********************************************************************/
{

  // this function should be called if PSNID_INPUTS.MCMC_NSTEP>0 and
  // type=PSNID_ITYPE_SNIA
  MCMC_RUN   = 1;

  // no point in running MCMC with too small number of steps
  if (PSNID_INPUTS.MCMC_NSTEP < 12000) {
    MCMC_NSTEP = 12000;
    MCMC_NBURN =  2000;
  } else {
    MCMC_NSTEP =  PSNID_INPUTS.MCMC_NSTEP;
    MCMC_NBURN =  MCMC_NSTEP/10;
  }

  return;
}


/**********************************************************************/
void psnid_best_run_mcmc(char *CCID, int itype, int nobs,
			 int *data_filt, double *data_mjd,
			 double *data_fluxcal, double *data_fluxcalerr,
			 int *useobs,
			 double mwebv, double mwebverr,
			 int zfix, double zprior, double zprior_err)
/**********************************************************************/
/*
  MCMC Main Engine

  Input:
    char *CCID  =
    int itype   =
    int nobs    =
    int *data_filt   =
    double *data_mjd = 
   

  Outline:
  - open output histogram files
  - read model light curves
  - set limits and step sizes
  - create histogram bins
  - MAIN PART OF MCMC; loop over
    - step in parameter space
    - interpolate between light curve model points
    - calculate chisq
    - decide to take step or not
  - histogram stats


  read_lc_all_mw
  get_mcmc_start
  lc_interp
  calc_chisq

 */
{
  int i, j, ngood;
  int take_step=0, *no_and_yes;
  double *this_model_day, **this_model_mags, **this_model_magserr;
  double *z_grid, *dm_grid;
  double zmin, zmax, dmmin, dmmax;
  double new_z, old_z, new_dm, old_dm, new_av, old_av,
    new_tmax, old_tmax, new_dmu, old_dmu;
  double new_chisq, old_chisq, rannum;
  double alpha, dmratio=1.0, avratio=1.0, zpratio=1.0,
    zratio=1.0, dmuratio=1.0;

  // z, dm, av, tmax, dmu, mu histograms
  int this_z=0, mcmc_nz_grid, *mcmc_z_hist,
    this_dm=0, mcmc_ndm_grid, *mcmc_dm_hist,
    this_av=0, mcmc_nav_grid, *mcmc_av_hist,
    this_tmax=0, mcmc_ntmax_grid, *mcmc_tmax_hist,
    this_dmu=0, mcmc_ndmu_grid, *mcmc_dmu_hist,
    this_mu=0, mcmc_nmu_grid, *mcmc_mu_hist;
  double *mcmc_z_grid, *mcmc_dm_grid, *mcmc_av_grid, *mcmc_tmax_grid,
    *mcmc_dmu_grid, *mcmc_mu_grid,
    mcmc_z_bin, mcmc_dm_bin, mcmc_av_bin, mcmc_tmax_bin,
    mcmc_dmu_bin, mcmc_mu_bin;

  // z vs mu histogram
  int mcmc_nz2_grid, mcmc_nmu2_grid, mcmc_nav2_grid,
    mcmc_z2_rebin, mcmc_mu2_rebin, mcmc_av2_rebin,
    **mcmc_zmu_hist, **mcmc_zav_hist;
  double *mcmc_z2_grid, *mcmc_mu2_grid, *mcmc_av2_grid,
    mcmc_z2_bin, mcmc_mu2_bin, mcmc_av2_bin;
  int mu_tmpi;
  double mu_first;

  double mu;
  double zref;
  double *z_moments, *dm_moments, *av_moments, *tmax_moments,
    *dmu_moments, *mu_moments;
  double *z_limits, *dm_limits, *av_limits, *tmax_limits,
    *dmu_limits, *mu_limits;

  //  char fnam[] = "psnid_best_run_mcmc";
  
  no_and_yes = ivector(0,2);
  no_and_yes[0] = 0;
  no_and_yes[1] = 0;


  // output MCMC chains and histograms
  //  zlab         = (char *)malloc(2);

  psnid_idum1=PSNID_BEST_ISEED1;
  psnid_idum2=PSNID_BEST_ISEED2;


  // read template and av files

  this_model_day     = dvector(1,PSNID_MAXND);
  this_model_mags    = dmatrix(1,PSNID_NFILTER, 1,PSNID_MAXND);
  this_model_magserr = dmatrix(1,PSNID_NFILTER, 1,PSNID_MAXND);



  z_grid  = dvector(1,PSNID_MAXNZ);
  dm_grid = dvector(1,PSNID_MAXNL);
  psnid_best_get_z_grid(z_grid);
  psnid_best_get_l_grid(dm_grid);

  zmin  = z_grid[1];
  zmax  = z_grid[PSNID_MAXNZ];
  dmmin = dm_grid[1];
  dmmax = dm_grid[PSNID_MAXNL];


  //  diagnostic
  /*
  printf("PSNID_MAXNL = %d\n", PSNID_MAXNL);
  for(i=1; i<=PSNID_MAXNL; i++) {
    printf("  i = %2d, THISND = %d\n", i-1, THISND[i-1]);
  }
  printf("  zmin = %8.5f\n", zmin);
  printf("  zmax = %8.5f\n", zmax);
  printf(" dmmin = %8.5f\n", dmmin);
  printf(" dmmax = %8.5f\n", dmmax);
  */


  /*  use best-fit as initial guess  */
  old_z   = PSNID_BEST_RESULTS.GRID_FINAL_PAR[zfix][itype][PSNID_PARAM_LOGZ];
  old_dm  = PSNID_BEST_RESULTS.GRID_FINAL_PAR[zfix][itype][PSNID_PARAM_SHAPEPAR];
  old_av   = PSNID_BEST_RESULTS.GRID_FINAL_PAR[zfix][itype][PSNID_PARAM_COLORPAR];
  old_tmax = PSNID_BEST_RESULTS.GRID_FINAL_PAR[zfix][itype][PSNID_PARAM_TMAX];
  old_dmu  = PSNID_BEST_RESULTS.GRID_FINAL_PAR[zfix][itype][PSNID_PARAM_DMU];
  /*
  printf("%f %f %f %f %f\n", old_z, old_dm, old_av, old_tmax, old_dmu);
  */


  // zref = reference redshift for determiming MCMC step sizes
  psnid_best_mcmc_zref(itype, &zref);
  //  printf("zref = %8.4f\n", zref);
  /*  scale step sizes according to photo-z  */
  PSNID_BEST_MCMC_DELTA_Z    = PSNID_BEST_MCMC_DELTA_Z_DEFAULT*(old_z/zref);
  PSNID_BEST_MCMC_DELTA_DM   = PSNID_BEST_MCMC_DELTA_DM_DEFAULT*(old_z/zref);
  PSNID_BEST_MCMC_DELTA_AV   = PSNID_BEST_MCMC_DELTA_AV_DEFAULT*(old_z/zref);
  PSNID_BEST_MCMC_DELTA_TMAX = PSNID_BEST_MCMC_DELTA_TMAX_DEFAULT*(old_z/zref);
  PSNID_BEST_MCMC_DELTA_DMU  = PSNID_BEST_MCMC_DELTA_DMU_DEFAULT*(old_z/zref);

  // non-flat z-prior
  //  if (zfix > 0) {
  if (zprior > 0.0) {
    old_z = zprior;
    PSNID_BEST_MCMC_DELTA_Z = PSNID_BEST_MCMC_DELTA_Z_DEFAULT*(zprior_err/0.05);
  }
  //  fprintf(outmcmc_ptr, "%8.5f %8.5f\n", zprior, zprior_err);
	  
  // don't step in dm for non-Ia
  if (itype != 0) PSNID_BEST_MCMC_DELTA_DM   = 0.000;
  //  if (itype == 2) PSNID_BEST_MCMC_DELTA_TMAX = PSNID_BEST_MCMC_DELTA_TMAX*2.0;


  /*  output histograms  */
  // z
  /*
  if (strcmp(SURVEY,"sdss") == 0) {
    mcmc_nz_grid = 351;
    mcmc_z_bin   = 0.002;
  } else if (strcmp(SURVEY,"des") == 0) {
    mcmc_nz_grid = 501;
    mcmc_z_bin   = 0.004;
  } else if (strcmp(SURVEY,"desY") == 0) {
    mcmc_nz_grid = 351;
    mcmc_z_bin   = 0.002;
  } else if (strcmp(SURVEY,"lsst") == 0) {
    mcmc_nz_grid = 501;
    mcmc_z_bin   = 0.004;
  } else if (strcmp(SURVEY,"ps1") == 0) {
    mcmc_nz_grid = 501;
    mcmc_z_bin   = 0.004;
  } else if (strcmp(SURVEY,"hstsn") == 0) {
    mcmc_nz_grid = 501;
    mcmc_z_bin   = 0.004;
  } else if (strcmp(SURVEY,"jpas56") == 0) {
    mcmc_nz_grid = 351;
    mcmc_z_bin   = 0.002;
  } else {
    mcmc_nz_grid = 351;
    mcmc_z_bin   = 0.002;
  }
  */
  mcmc_nz_grid = 501;
  //  mcmc_z_bin   = 0.002;
  mcmc_z_bin   = (pow(10.,SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_LOGZ][PSNID_MAXNZ])-
		  pow(10.,SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_LOGZ][1]))/(mcmc_nz_grid-1.0);
  mcmc_z_grid  = dvector(1,mcmc_nz_grid);
  mcmc_z_hist  = ivector(1,mcmc_nz_grid);
  for (i=1; i<=mcmc_nz_grid; i++) {
    mcmc_z_grid[i] = pow(10.,SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_LOGZ][1]) + (i-1)*mcmc_z_bin;
    mcmc_z_hist[i] = 0;
    //    printf("\t %3d %8.4f\n", i,mcmc_z_grid[i]);
  }
  // dm
  mcmc_ndm_grid = 651;
  //  mcmc_dm_bin   = 0.002;
  mcmc_dm_bin   = (SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_SHAPEPAR][PSNID_MAXNL]-
		   SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_SHAPEPAR][1])/(mcmc_ndm_grid-1.0);
  mcmc_dm_grid  = dvector(1,mcmc_ndm_grid);
  mcmc_dm_hist  = ivector(1,mcmc_ndm_grid);
  for (i=1; i<=mcmc_ndm_grid; i++) {
    mcmc_dm_grid[i] = SNGRID_PSNID[PSNID_THIS_TYPE].VALUE[IPAR_GRIDGEN_SHAPEPAR][1] + (i-1)*mcmc_dm_bin;
    mcmc_dm_hist[i] = 0;
  }
  // av
  mcmc_nav_grid = 1201;
  //  mcmc_av_bin   = 0.005;
  mcmc_av_bin   = (PSNID_ASTEP*(PSNID_MAXNA-1.0))/(mcmc_nav_grid-1.0);
  mcmc_av_grid  = dvector(1,mcmc_nav_grid);
  mcmc_av_hist  = ivector(1,mcmc_nav_grid);
  for (i=1; i<=mcmc_nav_grid; i++) {
    mcmc_av_grid[i] = PSNID_AMIN + (i-1)*mcmc_av_bin;
    mcmc_av_hist[i] = 0;
  }
  // tmax
  mcmc_ntmax_grid = 1001;
  mcmc_tmax_bin   = 0.02;
  mcmc_tmax_grid  = dvector(1,mcmc_ntmax_grid);
  mcmc_tmax_hist  = ivector(1,mcmc_ntmax_grid);
  for (i=1; i<=mcmc_ntmax_grid; i++) {
    mcmc_tmax_grid[i] = ((int) old_tmax - 10.0) + (i-1)*mcmc_tmax_bin;
    mcmc_tmax_hist[i] = 0;
  }
  // dmu
  mcmc_ndmu_grid = 1001;
  //  mcmc_dmu_bin   = 0.005;
  mcmc_dmu_bin   = (PSNID_USTEP*(PSNID_MAXNU-1.0))/(mcmc_ndmu_grid-1.0);
  mcmc_dmu_grid  = dvector(1,mcmc_ndmu_grid);
  mcmc_dmu_hist  = ivector(1,mcmc_ndmu_grid);
  for (i=1; i<=mcmc_ndmu_grid; i++) {
    mcmc_dmu_grid[i] = PSNID_UMIN + (i-1)*mcmc_dmu_bin;
    mcmc_dmu_hist[i] = 0;
  }
  // mu
  MCMC_REDSHIFT = old_z;
  //  mu       = dl(2);
  mu    = dLmag(old_z, old_z, &PSNID_INPUTS.HzFUN_INFO);
  mu = mu + old_dmu;
  mcmc_nmu_grid = 1001;
  //  mcmc_mu_bin   = 0.005;
  mcmc_mu_bin   = (PSNID_USTEP*(PSNID_MAXNU-1.0))/(mcmc_ndmu_grid-1.0);
  mcmc_mu_grid  = dvector(1,mcmc_nmu_grid);
  mcmc_mu_hist  = ivector(1,mcmc_nmu_grid);
  // let first mu bin to be XX.X0
  mu_tmpi  = (int) ((mu+0.05)*10);
  mu_first = 0.1* (double) mu_tmpi;
  for (i=1; i<=mcmc_nmu_grid; i++) {
    mcmc_mu_grid[i] = (mu_first + PSNID_UMIN) + (i-1)*mcmc_mu_bin;
    mcmc_mu_hist[i] = 0;
  }
  /*
  printf("\t old_z = %9.4f   mu = %9.4f  mu_tmpi = %d   mu_first = %9.4f "
	 "\t mcmc_mu_grid[1] = %9.4f  mcmc_mu_grid[last] = %9.4f\n",
 	 old_z, mu, mu_tmpi, mu_first,
	 mcmc_mu_grid[1], mcmc_mu_grid[1001]);
  fflush(stdout);
  */

  //////////////////////
  // z vs mu contours //
  mcmc_z2_rebin  = 2;
  mcmc_mu2_rebin = 4;
  mcmc_nz2_grid  = ((mcmc_nz_grid-1)/mcmc_z2_rebin)+1;
  mcmc_nmu2_grid = ((mcmc_nmu_grid-1)/mcmc_mu2_rebin)+1;
  mcmc_z2_bin    = mcmc_z_bin*mcmc_z2_rebin;
  mcmc_mu2_bin   = mcmc_mu_bin*mcmc_mu2_rebin;
  mcmc_z2_grid   = dvector(1,mcmc_nz2_grid);   // xgrid
  mcmc_mu2_grid  = dvector(1,mcmc_nmu2_grid);  // ygrid
  mcmc_zmu_hist  = imatrix(1,mcmc_nz2_grid, 1,mcmc_nmu2_grid);  // hist
  for (i=1; i<=mcmc_nz2_grid; i++) {
    mcmc_z2_grid[i] = (i-1)*mcmc_z2_bin;
    for (j=1; j<=mcmc_nmu2_grid; j++) {
      mcmc_zmu_hist[i][j] = 0;
    }
  }
  for (i=1; i<=mcmc_nmu2_grid; i++) {
    mcmc_mu2_grid[i] = (mu_first - 2.5) + (i-1)*mcmc_mu2_bin;
  }
  //////////////////////

  //////////////////////
  // z vs av contours //
  //  mcmc_z2_rebin  = 2;
  mcmc_av2_rebin = 4;
  //  mcmc_nz2_grid  = ((mcmc_nz_grid-1)/mcmc_z2_rebin)+1;
  mcmc_nav2_grid = ((mcmc_nav_grid-1)/mcmc_av2_rebin)+1;
  //  mcmc_z2_bin    = mcmc_z_bin*mcmc_z2_rebin;
  mcmc_av2_bin   = mcmc_av_bin*mcmc_av2_rebin;
  //  mcmc_z2_grid   = dvector(1,mcmc_nz2_grid);   // xgrid
  mcmc_av2_grid  = dvector(1,mcmc_nav2_grid);  // ygrid
  mcmc_zav_hist  = imatrix(1,mcmc_nz2_grid, 1,mcmc_nav2_grid);  // hist
  for (i=1; i<=mcmc_nz2_grid; i++) {
    //    mcmc_z2_grid[i] = (i-1)*mcmc_z2_bin;
    for (j=1; j<=mcmc_nav2_grid; j++) {
      mcmc_zav_hist[i][j] = 0;
    }
  }
  for (i=1; i<=mcmc_nav2_grid; i++) {
    mcmc_av2_grid[i] = -1.00 + (i-1)*mcmc_av2_bin;
  }
  //////////////////////


  //////////////////////////////
  ////  main part of MCMC   ////

  old_chisq = PSNID_BIGN;

  for(i = 1; i <= MCMC_NSTEP; i++) {

    if (i == 1) {
      new_z    = old_z;
      new_dm   = old_dm;
      new_av   = old_av;
      new_tmax = old_tmax;
      new_dmu  = old_dmu;
    } else {
      new_z    = old_z    +    PSNID_BEST_MCMC_DELTA_Z*gasdev(&psnid_idum2);
      new_dm   = old_dm   +   PSNID_BEST_MCMC_DELTA_DM*gasdev(&psnid_idum2);
      new_av   = old_av   +   PSNID_BEST_MCMC_DELTA_AV*gasdev(&psnid_idum2);
      new_tmax = old_tmax + PSNID_BEST_MCMC_DELTA_TMAX*gasdev(&psnid_idum2);
      new_dmu  = old_dmu  +  PSNID_BEST_MCMC_DELTA_DMU*gasdev(&psnid_idum2);
    }

    /*
    // DEBUG
    if (i <= 10) {
      printf(" i = %d  ", i);
      printf("old z dm av tmax dmu = %9.5f %8.3f %8.3f %10.2f %8.3f   ",
	     old_z, old_dm, old_av, old_tmax, old_dmu);
      printf("new z dm av tmax dmu = %9.5f %8.3f %8.3f %10.2f %8.3f\n",
	     new_z, new_dm, new_av, new_tmax, new_dmu);
      printf("\t MCMC i = %d, old_chisq = %8.3f  new_chisq = %8.3f\n",
	     i, old_chisq, new_chisq);
    }
    // DEBUG
    */

    // peg at min/max values
    if (new_z < zmin) new_z = zmin;
    if (new_z > zmax) new_z = zmax;
    if (new_dm < dmmin) new_dm = dmmin;
    if (new_dm > dmmax) new_dm = dmmax;

    // compute light curve from grid of models
    //  - interpolate for (z,dm) values

    psnid_best_lc_interp(itype,
			 new_z, new_dm, new_av, new_tmax, new_dmu,
			 PSNID_MODEL_EPOCH,
			 PSNID_MODEL_MAG, PSNID_MODEL_MAGERR,
			 PSNID_MODEL_EXTINCT, PSNID_MODEL_MWEXTINCT,
			 this_model_day, this_model_mags,
			 this_model_magserr);

    /*
    // DEBUG
    if (i==1) {
      for (j=1; j<PSNID_MAXND; j++) {
	printf("\t epoch = %10.2f  mag = %8.4f  mage = %8.4f\n",
	       this_model_day[j],
	       this_model_mags[3][j],
	       this_model_magserr[3][j]);
      }
    }
    // DEBUG
    */

    ngood     = 0;
    new_chisq = 0.0;

    psnid_best_calc_chisq(nobs, useobs,
			  data_filt, data_mjd, data_fluxcal, data_fluxcalerr,
			  this_model_day, this_model_mags, this_model_magserr,
			  &ngood, &new_chisq, 0, 0);


    /////   dm15 prior   /////
    /*
    if (PSNID_USE_DM_PRIOR == 1 && itype == 0) {
      //      dmratio = dmprior(itype, -1, new_dm)/dmprior(itype, -1, old_dm);
      dmratio = 1.0;
    } else {
      dmratio = 1.0;
    }
    */
    //    printf("%12.5e  %12.5e  %12.5e  %12.5e  %12.5e    ",
    //	   dmratio, new_dm, old_dm,
    //	   dmprior(itype, -1, new_dm), dmprior(itype, -1, old_dm));

    /////   av prior   /////
    if (PSNID_USE_AV_PRIOR == 1) {
      avratio = psnid_best_avprior1(itype, new_av)/
	psnid_best_avprior1(itype, old_av);
    } else {
      avratio = 1.0;
    }
    //    printf("%12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n",
    //	   avratio, new_av, old_av,
    //	   avprior1(itype, new_av), avprior1(itype, old_av));

    /////    z distribution prior   /////
    /*
    if (PSNID_USE_Z_PRIOR == 1) {
      if (PSNID_Z_PRIOR_TYPE == 1) {
	zpratio = pow(new_z,PSNID_Z_PRIOR_SLOPE)/pow(old_z,PSNID_Z_PRIOR_SLOPE);
      }
    } else {
      zpratio = 1.0;
    }
    */

    /////    zprior    /////
    //    if (zfix > 0) {
    if (zprior > 0.0) {
      zratio = exp(-pow((zprior-new_z)/zprior_err,2)/2.0)/
	exp(-pow((zprior-old_z)/zprior_err,2)/2.0);
    } else {
      zratio = 1.0;
    }

    //////  dmu prior for CC SNe  /////
    /*
    if (MCMC_CC_DMU_PRIOR == 1 && itype > 0) {
      dmuratio = exp(-pow(new_dmu/MCMC_CC_DMU_SIGMA,2)/2.0)/
	exp(-pow(old_dmu/MCMC_CC_DMU_SIGMA,2)/2.0);
    } else {
      dmuratio = 1.0;
    }
    */

    alpha = zratio*zpratio*avratio*dmratio*dmuratio*
      exp((old_chisq - new_chisq)/2.);

    if (alpha >= 1.0) {
      take_step = 1;
    } else {
      rannum = ran2(&psnid_idum1);
      take_step = (rannum < alpha) ? 1 : 0;
    }

    /*
    fprintf(outmcmc_ptr,
	    "%3d   old %12.5e  -->  new %12.5e    %1d ",
	    ngood, old_chisq, new_chisq, take_step);
    */

    if (take_step == 1) {

      MCMC_REDSHIFT = new_z;
      //      mu       = dl(2);  // distance modulus

      mu  = dLmag(old_z, old_z, &PSNID_INPUTS.HzFUN_INFO);

      //      printf("%8.4f  %13.5e\n",MCMC_REDSHIFT, mu);

      // output new parameter values
      /*
      if (i <= MCMC_NOUT) {
	fprintf(outmcmc_ptr,
		"1  %9d   %12.5e  %12.5e  %12.5e  %10.4f  %11.4e  %12.5e\n",
		i, new_z, new_dm, new_av, new_tmax, new_dmu, mu+new_dmu);
      }
      */

      old_z     = new_z;
      old_dm    = new_dm;
      old_av    = new_av;
      old_tmax  = new_tmax;
      old_dmu   = new_dmu;
      old_chisq = new_chisq;

    } else if (take_step == 0) {

      MCMC_REDSHIFT = old_z;

      mu  = dLmag(old_z, old_z, &PSNID_INPUTS.HzFUN_INFO);

      // output old parameter values
      /*
      if (i <= MCMC_NOUT) {
	fprintf(outmcmc_ptr,
		"0  %9d   %12.5e  %12.5e  %12.5e  %10.4f  %11.4e  %12.5e\n",
		i, old_z, old_dm, old_av, old_tmax, old_dmu, mu+old_dmu);
      }
      */

    }

    no_and_yes[take_step] += 1;

    if (i > MCMC_NBURN) {

      // histogram z, dm, av, tmax, dmu, mu
      hunt(mcmc_z_grid, mcmc_nz_grid, old_z, &this_z);
      if (this_z >= 1 && this_z <= mcmc_nz_grid) mcmc_z_hist[this_z] += 1;
      hunt(mcmc_dm_grid, mcmc_ndm_grid, old_dm, &this_dm);
      if (this_dm >= 1 && this_dm <= mcmc_ndm_grid) mcmc_dm_hist[this_dm] += 1;
      hunt(mcmc_av_grid, mcmc_nav_grid, old_av, &this_av);
      if (this_av >= 1 && this_av <= mcmc_nav_grid) mcmc_av_hist[this_av] += 1;
      hunt(mcmc_tmax_grid, mcmc_ntmax_grid, old_tmax, &this_tmax);
      if (this_tmax >= 1 && this_tmax <= mcmc_ntmax_grid) mcmc_tmax_hist[this_tmax] += 1;
      hunt(mcmc_dmu_grid, mcmc_ndmu_grid, old_dmu, &this_dmu);
      if (this_dmu >= 1 && this_dmu <= mcmc_ndmu_grid) mcmc_dmu_hist[this_dmu] += 1;
      hunt(mcmc_mu_grid, mcmc_nmu_grid, mu+old_dmu, &this_mu);
      if (this_mu >= 1 && this_mu <= mcmc_nmu_grid) mcmc_mu_hist[this_mu] += 1;

      // histogram z vs mu
      /*
      hunt(mcmc_z2_grid, mcmc_nz2_grid, old_z, &this_z2);
      hunt(mcmc_mu2_grid, mcmc_nmu2_grid, mu+old_dmu, &this_mu2);
      if (this_z2 >= 1 && this_z2 <= mcmc_nz2_grid &&
	  this_mu2 >= 1 && this_mu2 <= mcmc_nmu2_grid) {
	mcmc_zmu_hist[this_z2][this_mu2] += 1;
      }
      */
      // histogram z vs mu
      /*
      hunt(mcmc_av2_grid, mcmc_nav2_grid, old_av, &this_av2);
      if (this_z2 >= 1 && this_z2 <= mcmc_nz2_grid &&
	  this_av2 >= 1 && this_av2 <= mcmc_nav2_grid) {
	mcmc_zav_hist[this_z2][this_av2] += 1;
      }
      */
    }

  }
  //////////////////////////////



  //////////////////////////////////////////////////////////
  // calculate first 4 moments of distributions
  // NOTE:
  //   mean     = moments[1]
  //   sigma    = sqrt(moments[2])
  //   skewness = moments[3]/pow(moments[2],1.5)
  //   kurtosis = moments[4]/pow(moments[2],2)

  z_moments    = dvector(1,4);
  dm_moments   = dvector(1,4);
  av_moments   = dvector(1,4);
  tmax_moments = dvector(1,4);
  dmu_moments  = dvector(1,4);
  mu_moments   = dvector(1,4);
  psnid_best_histstats(mcmc_nz_grid, mcmc_z_grid, mcmc_z_hist, z_moments);
  psnid_best_histstats(mcmc_ndm_grid, mcmc_dm_grid, mcmc_dm_hist, dm_moments);
  psnid_best_histstats(mcmc_nav_grid, mcmc_av_grid, mcmc_av_hist, av_moments);
  psnid_best_histstats(mcmc_ntmax_grid, mcmc_tmax_grid, mcmc_tmax_hist, tmax_moments);
  psnid_best_histstats(mcmc_ndmu_grid, mcmc_dmu_grid, mcmc_dmu_hist, dmu_moments);
  psnid_best_histstats(mcmc_nmu_grid, mcmc_mu_grid, mcmc_mu_hist, mu_moments);

  //////////////////////////////////////////////////////////


  /*
  printf("\nz   %13.5e  %13.5e  %13.5e  %13.5e\n",
	 z_moments[1], z_moments[2], z_moments[3], z_moments[4]);
  printf("dm  %13.5e  %13.5e  %13.5e  %13.5e\n",
	 dm_moments[1], dm_moments[2], dm_moments[3], dm_moments[4]);
  printf("AV  %13.5e  %13.5e  %13.5e  %13.5e\n",
	 av_moments[1], av_moments[2], av_moments[3], av_moments[4]);
  printf("Tm  %13.5e  %13.5e  %13.5e  %13.5e\n",
	 tmax_moments[1], tmax_moments[2], tmax_moments[3], tmax_moments[4]);
  printf("dmu %13.5e  %13.5e  %13.5e  %13.5e\n",
	 dmu_moments[1], dmu_moments[2], dmu_moments[3], dmu_moments[4]);
  printf("mu  %13.5e  %13.5e  %13.5e  %13.5e\n",
	 mu_moments[1], mu_moments[2], mu_moments[3], mu_moments[4]);
  */


  /*
  printf(" %9d/%9d = %8.3f percent accept\n",
	 no_and_yes[1], MCMC_NSTEP,
	 (double) 100.*no_and_yes[1]/MCMC_NSTEP);
  */

  //////////////////////////////////////////////////////////
  //  calculate 50%, +/- 1 and 2 sigma limits

  z_limits    = dvector(1,5);
  dm_limits   = dvector(1,5);
  av_limits   = dvector(1,5);
  tmax_limits = dvector(1,5);
  dmu_limits  = dvector(1,5);
  mu_limits   = dvector(1,5);

  psnid_best_param_limits(mcmc_nz_grid, mcmc_z_grid, mcmc_z_hist, z_limits);
  psnid_best_param_limits(mcmc_ndm_grid, mcmc_dm_grid, mcmc_dm_hist, dm_limits);
  psnid_best_param_limits(mcmc_nav_grid, mcmc_av_grid, mcmc_av_hist, av_limits);
  psnid_best_param_limits(mcmc_ntmax_grid, mcmc_tmax_grid, mcmc_tmax_hist, tmax_limits);
  psnid_best_param_limits(mcmc_ndmu_grid, mcmc_dmu_grid, mcmc_dmu_hist, dmu_limits);
  psnid_best_param_limits(mcmc_nmu_grid, mcmc_mu_grid, mcmc_mu_hist, mu_limits);


  // store results in global variable
  for (i=0; i<5; i++) {
    PSNID_BEST_RESULTS.MCMC_FINAL_PAR[zfix][itype][PSNID_PARAM_LOGZ][i]     = z_limits[i+1];
    PSNID_BEST_RESULTS.MCMC_FINAL_PAR[zfix][itype][PSNID_PARAM_SHAPEPAR][i] = dm_limits[i+1];
    PSNID_BEST_RESULTS.MCMC_FINAL_PAR[zfix][itype][PSNID_PARAM_COLORPAR][i] = av_limits[i+1];
    PSNID_BEST_RESULTS.MCMC_FINAL_PAR[zfix][itype][PSNID_PARAM_TMAX][i]     = tmax_limits[i+1];
    // MCMC dmu stores mu!
    PSNID_BEST_RESULTS.MCMC_FINAL_PAR[zfix][itype][PSNID_PARAM_DMU][i]      = mu_limits[i+1];
  }

  /*
  printf("%s       z = %8.4f  (%+9.4f %+9.4f )\n", outpref,
  	 z_limits[3],z_limits[4]-z_limits[3],z_limits[2]-z_limits[3]);
  printf("%s      mu = %8.4f  (%+9.4f %+9.4f )\n", outpref,
	 mu_limits[3],mu_limits[4]-mu_limits[3],mu_limits[2]-mu_limits[3]);
  printf("%s      av = %8.4f  (%+9.4f %+9.4f )\n", outpref,
	 av_limits[3],av_limits[4]-av_limits[3],av_limits[2]-av_limits[3]);
  printf("%s      dm = %8.4f  (%+9.4f %+9.4f )\n", outpref,
	 dm_limits[3],dm_limits[4]-dm_limits[3],dm_limits[2]-dm_limits[3]);
  printf("%s    tmax = %8.2f  (%+9.2f %+9.2f )\n", outpref,
  	 tmax_limits[3],tmax_limits[4]-tmax_limits[3],tmax_limits[2]-tmax_limits[3]);
  */

  free_ivector(no_and_yes,  0,2);

  free_dvector(z_moments,   1,4);
  free_dvector(dm_moments,  1,4);
  free_dvector(av_moments,  1,4);
  free_dvector(tmax_moments,1,4);
  free_dvector(dmu_moments, 1,4);
  free_dvector(mu_moments,  1,4);

  free_dvector(z_limits,    1,5);
  free_dvector(dm_limits,   1,5);
  free_dvector(av_limits,   1,5);
  free_dvector(tmax_limits, 1,5);
  free_dvector(dmu_limits,  1,5);
  free_dvector(mu_limits,   1,5);

  free_dvector(mcmc_z_grid, 1,mcmc_nz_grid);
  free_ivector(mcmc_z_hist, 1,mcmc_nz_grid);
  free_dvector(mcmc_dm_grid, 1,mcmc_ndm_grid);
  free_ivector(mcmc_dm_hist, 1,mcmc_ndm_grid);
  free_dvector(mcmc_av_grid, 1,mcmc_nav_grid);
  free_ivector(mcmc_av_hist, 1,mcmc_nav_grid);
  free_dvector(mcmc_tmax_grid, 1,mcmc_ntmax_grid);
  free_ivector(mcmc_tmax_hist, 1,mcmc_ntmax_grid);
  free_dvector(mcmc_dmu_grid, 1,mcmc_ndmu_grid);
  free_ivector(mcmc_dmu_hist, 1,mcmc_ndmu_grid);
  free_dvector(mcmc_mu_grid, 1,mcmc_nmu_grid);
  free_ivector(mcmc_mu_hist, 1,mcmc_nmu_grid);

  free_dvector(mcmc_z2_grid,  1,mcmc_nz2_grid);
  free_dvector(mcmc_mu2_grid, 1,mcmc_nmu2_grid);
  free_imatrix(mcmc_zmu_hist, 1,mcmc_nz2_grid, 1,mcmc_nmu2_grid);

  free_dvector(mcmc_av2_grid, 1,mcmc_nav2_grid);
  free_imatrix(mcmc_zav_hist, 1,mcmc_nz2_grid, 1,mcmc_nav2_grid);

  //////////////////////////////////////////////////////////




  free_dvector(this_model_day, 1,PSNID_MAXND);
  free_dmatrix(this_model_mags, 1,PSNID_NFILTER,1,PSNID_MAXND);
  free_dmatrix(this_model_magserr, 1,PSNID_NFILTER,1,PSNID_MAXND);

  free_dvector(z_grid, 1,PSNID_MAXNZ);
  free_dvector(dm_grid, 1,PSNID_MAXNL);


  return;
}
// end of psnid_best_run_mcmc


/**********************************************************************/
void psnid_best_mcmc_zref(int itype, double *zref)
/**********************************************************************/
{

  if (itype == 0) {  // SN Ia
    *zref = 0.40*pow(10,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMAX[IPAR_GRIDGEN_LOGZ]);
  } else {
    *zref = 0.10*pow(10,SNGRID_PSNID[TYPEINDX_SNIA_PSNID].VALMAX[IPAR_GRIDGEN_LOGZ]);
  }


  return;
}
// end of mcmc_zref


/**********************************************************************/
void psnid_best_histstats(int ngrid, double *xgrid, int *ygrid,
			  double *moments)
/**********************************************************************/
{
  int i, ysum;
  double xysum, xave, x2ysum, x3ysum, x4ysum;


  for (i=1; i<=4; i++) moments[i] = 0.0;

  // calculate mean
  ysum  = 0 ;
  xysum = 0.0;
  for (i = 1; i < ngrid; i++) {
    xave   = (xgrid[i] + xgrid[i+1])/2.;
    ysum  += ygrid[i];
    xysum += xave*ygrid[i];
  }

  if (ysum > 0.0 || ysum < 0.0) {
    moments[1] = xysum/((double) ysum);

    // and higher moments
    x2ysum = 0.0;
    x3ysum = 0.0;
    x4ysum = 0.0;
    for (i = 1; i < ngrid; i++) {
      xave   = (xgrid[i] + xgrid[i+1])/2.;
      //    ysum  += ygrid[i];
      x2ysum += pow(xave - moments[1],2)*ygrid[i];
      x3ysum += pow(xave - moments[1],3)*ygrid[i];
      x4ysum += pow(xave - moments[1],4)*ygrid[i];
    }
    moments[2] = x2ysum/((double) ysum);
    moments[3] = x3ysum/((double) ysum);
    moments[4] = x4ysum/((double) ysum);

  } else {

    for (i=1; i<=4; i++) moments[i] = PSNID_BIGN;

  }


  return;
}
// end of psnid_best_histstats


/**********************************************************************/
void psnid_best_param_limits(int ngrid, double *xgrid, int *ygrid,
			     double *params)
/**********************************************************************/
{
  int i, this_p=0;
  double *cumdist, *pval;

  pval    = dvector(1,5);
  /*
  pval[1] = 0.500000;    // median
  pval[2] = 0.158657;    // -1 sigma
  pval[3] = 0.841343;    // +1 sigma
  pval[4] = 0.022750;    // -2 sigma
  pval[5] = 0.977250;    // +2 sigma
  */

  pval[1] = 0.022750;    // -2 sigma
  pval[2] = 0.158657;    // -1 sigma
  pval[3] = 0.500000;    // median
  pval[4] = 0.841343;    // +1 sigma
  pval[5] = 0.977250;    // +2 sigma

  // first calculate cumulative distribution
  cumdist = dvector(1, ngrid);
  for (i=1; i<=ngrid; i++) cumdist[i] = 0.0;
  // xxx  nmcmc   = (double) (MCMC_NSTEP - MCMC_NBURN);

  cumdist[1] = (double) ygrid[1];
  for (i=2; i<=ngrid; i++) {
    cumdist[i] = cumdist[i-1] + ygrid[i];
    //    cumdist[i] = cumdist[i-1] + ygrid[i]/nmcmc;
    //    printf("%f %f\n", xgrid[i], cumdist[i]);
  }
  for (i=1; i<=ngrid; i++) cumdist[i] = cumdist[i]/cumdist[ngrid];

  for (i=1; i<=5; i++) {
    // cumdist[this_p] <= pval[i] < cumdist[this_p+1]
    hunt(cumdist, ngrid, pval[i], &this_p);
    if (this_p < ngrid) {
      if (this_p < 1) {
	params[i] = xgrid[1];
      } else {
	params[i] = 0.5*(xgrid[this_p]+xgrid[this_p+1]) + 
	  (xgrid[this_p+1]-xgrid[this_p])*(pval[i]-cumdist[this_p])/
			 		  (cumdist[this_p+1]-cumdist[this_p]);
      }
    } else {
      params[i] = xgrid[this_p];
    }
    //    printf("i = %d    %f\n", this_p, params[i]);
  }
  //  printf("nmcmc = %f   total = %f\n",nmcmc,cumdist[ngrid]);
  //  printf("\n");

  free_dvector(pval,    1,5);
  free_dvector(cumdist, 1,ngrid);

  return;
}
// end of psnid_best_param_limits


/**********************************************************************/
void psnid_best_lc_interp(int itype, 
			  double z, double shape,
			  double color, double tmax, double dmu,
			  double ***day,
			  double ****mags, double ****magserr,
			  double ****extinct, double ****mwextinct,
			  double *out_day, double **out_mags,
			  double **out_magserr)
/**********************************************************************/
/*
  Interpolates light curves between REDSHIFT and SHAPEPAR grid points.

  input:

    itype     = 0,1,2 for Ia,Ibc,II
    z         = REDSHIFT
    shape     = SHAPEPAR
    color     = COLORPAR
    tmax      = Tmax
    dmu       = delta mu
    day       = grid model epoch
    mags      = grid model mags
    magserr   = grid model mag error
    extinct   = grid model host extinction
    mwextinct = grid model MW extinction

  output:

    out_day     =
    out_mags    =
    out_magserr =

 */
{
  int i, j, k, l, this_l=0, this_z=0;
  double fracz, fracl, **weights;
  double *z_grid, *l_grid;


  // REDSHIFT
  z_grid = dvector(1, PSNID_MAXNZ);
  psnid_best_get_z_grid(z_grid);
  hunt(z_grid, PSNID_MAXNZ, z, &this_z);
  // z_grid[this_z] <= z < z_grid[this_z+1]
  if (this_z < 1)   this_z = 1;
  if (this_z > PSNID_MAXNZ) this_z = PSNID_MAXNZ;

  // SHAPEPAR
  l_grid = dvector(1, PSNID_MAXNL);
  psnid_best_get_l_grid(l_grid);
  hunt(l_grid, PSNID_MAXNL, shape, &this_l);
  // l_grid[this_l] <= shape < l_grid[this_l+1]
  if (this_l < 1)  this_l = 1;
  if (this_l > PSNID_MAXNL) this_l = PSNID_MAXNL;

  fracz  = (z - z_grid[this_z])/(z_grid[this_z+1] - z_grid[this_z]);
  fracl = (shape - l_grid[this_l])/(l_grid[this_l+1] - l_grid[this_l]);

  // don't interpolate dm for non-Ia
  if (itype != 0) fracl = 0.0;

  weights = dmatrix(0,2, 0,2);
  psnid_best_z_l_weights(fracz, fracl, weights);


  /*
  // DEBUG
  printf("%13.5e  %8.5f %8.5f %8.5f   %3d       ",
  	 fracz, z_grid[this_z], z, z_grid[this_z+1], this_z);
  printf("%13.5e  %8.5f %8.5f %8.5f   %3d       ",
  	 fracl, l_grid[this_l], shape, l_grid[this_l+1], this_l);
  printf("%13.5e  %13.5e  %13.5e  %13.5e\n",
	 weights[0][0],weights[0][1],weights[1][0],weights[1][1]);
  // DEBUG
  */

  // lc is weighted average of 4 (REDSHIFT,SHAPEPAR) grid points
  for (i = 1; i <= PSNID_MAXND; i++) {
    out_day[i] = 0.0;
    for (k=0; k<2; k++) {    // z
      for (l=0; l<2; l++) {  // shape
	out_day[i] += weights[k][l]*day[this_l+l][this_z+k][i];
      }
    }
    out_day[i] += tmax;
    for (j = 1; j <= PSNID_NFILTER; j++) {
      out_mags[j][i]    = 0.0;
      out_magserr[j][i] = 0.0;
      for (k=0; k<2; k++) {
	for (l=0; l<2; l++) {
	  out_mags[j][i] +=
	    weights[k][l]*(mags[j][this_l+l][this_z+k][i] -
			   (color-PSNID_BASE_COLOR[itype])*extinct[j][this_l+l][this_z+k][i]);
	  out_magserr[j][i] +=
	    weights[k][l]*(magserr[j][this_l+l][this_z+k][i]);
	}
      }
      out_mags[j][i] += dmu;
    }
  }



  free_dvector(z_grid, 1, PSNID_MAXNZ);
  free_dvector(l_grid, 1, PSNID_MAXNL);
  free_dmatrix(weights, 0,2, 0,2);


  return;
}
// end of psnid_best_lc_interp


/**********************************************************************/
void psnid_best_z_l_weights(double z, double shape, double **w)
/**********************************************************************/
{

  // w[z][shape]

  w[0][0] = (1.-shape)*(1.-z);
  w[0][1] = shape*(1.-z);
  w[1][0] = (1.-shape)*z;
  w[1][1] = shape*z;

  /*
  printf("%13.5e  %13.5e    %13.5e  %13.5e  %13.5e  %13.5e\n",
	 z, shape, w[0][0],w[0][1],w[1][0],w[1][1]);
  */

  return;
}
// end of psnid_best_z_l_weights


/************************************************************************/
/************************************************************************/
/*************                                           ****************/
/*************                RK Functions               ****************/
/*************                                           ****************/
/************************************************************************/
/************************************************************************/



// ================================================
void PSNID_BEST_INIT(void) {

  // Mar 31 2013
  // One-time initializations before fitting and before data are read.

  int NVAR ;
  char fnam[] = "PSNID_BEST_INIT" ;

  // ------------- BEGIN -----------

  print_banner(fnam);

  psnid_best_setup_searchgrid();
  psnid_best_check_for_avprior();
  psnid_best_check_for_dmu();
  psnid_best_reset_results(); // set PSNID_TYPE_NAME

  NVAR = snana_nearnbr_rdinput__();

  if ( NVAR>0 && NVAR != PSNID_BEST_NEARNBR_NVAR ) {
    sprintf(c1err,"&SNLCINP NEARNBR_XXX specified NVAR=%d", NVAR);
    sprintf(c2err,"but PSNID_BEST requires NVAR=%d", PSNID_BEST_NEARNBR_NVAR);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  PSNID_NEARNBR.NVAR = NVAR ;

  // Feb 2017 RK - setup cut-array to fix fragile logic elsewhere
  PSNID_FITPROB_CUTLIST[PSNID_ITYPE_SNIA]   = PSNID_INPUTS.FITPROB_IA_CUT;
  PSNID_FITPROB_CUTLIST[PSNID_ITYPE_SNIBC]  = PSNID_INPUTS.FITPROB_IBC_CUT;
  PSNID_FITPROB_CUTLIST[PSNID_ITYPE_SNII]   = PSNID_INPUTS.FITPROB_II_CUT;
  PSNID_FITPROB_CUTLIST[PSNID_ITYPE_PEC1A]  = PSNID_INPUTS.FITPROB_PEC1A_CUT;
  PSNID_FITPROB_CUTLIST[PSNID_ITYPE_MODEL1] = PSNID_INPUTS.FITPROB_MODEL1_CUT;
  PSNID_FITPROB_CUTLIST[PSNID_ITYPE_MODEL2] = PSNID_INPUTS.FITPROB_MODEL2_CUT;
  PSNID_FITPROB_CUTLIST[PSNID_ITYPE_MODEL3] = PSNID_INPUTS.FITPROB_MODEL3_CUT;
  PSNID_FITPROB_CUTLIST[PSNID_ITYPE_MODEL4] = PSNID_INPUTS.FITPROB_MODEL4_CUT;

  PSNID_PBAYES_CUTLIST[PSNID_ITYPE_SNIA]   = PSNID_INPUTS.PBAYES_IA_CUT;
  PSNID_PBAYES_CUTLIST[PSNID_ITYPE_SNIBC]  = PSNID_INPUTS.PBAYES_IBC_CUT;
  PSNID_PBAYES_CUTLIST[PSNID_ITYPE_SNII]   = PSNID_INPUTS.PBAYES_II_CUT;
  PSNID_PBAYES_CUTLIST[PSNID_ITYPE_PEC1A]  = PSNID_INPUTS.PBAYES_PEC1A_CUT;
  PSNID_PBAYES_CUTLIST[PSNID_ITYPE_MODEL1] = PSNID_INPUTS.PBAYES_MODEL1_CUT;
  PSNID_PBAYES_CUTLIST[PSNID_ITYPE_MODEL2] = PSNID_INPUTS.PBAYES_MODEL2_CUT;
  PSNID_PBAYES_CUTLIST[PSNID_ITYPE_MODEL3] = PSNID_INPUTS.PBAYES_MODEL3_CUT;
  PSNID_PBAYES_CUTLIST[PSNID_ITYPE_MODEL4] = PSNID_INPUTS.PBAYES_MODEL4_CUT;

} // end of   PSNID_BEST_INIT


// *****************************************************
void PSNID_BEST_INIT_SNTABLE(int OPT, char *TEXTFORMAT, int LSIM) {

  // Created May 19 2014
  // Init output TABLE based on &SNLCINP namelist string SNTABLE_LIST.
  //
  // Inputs:
  //   OPT=1 --> do FITRES table
  //   OPT=2 --> include epoch residuals
  //   TEXTFORMAT --> write [TEXTFILE_PREFIX].TEXT.FITRES in this format.
  //         (TEXTFORMAT='' --> no text file)
  //   LSIM=1 --> is a simulation
  //
  // May 20 2019: abort if   PSNID_MAXNL_NONIA exceeds bound.

  char fnam[] = "PSNID_BEST_INIT_SNTABLE" ;

  // ------------ BEGIN --------------

  // store flags in globals
  PSNID_INPUTS.WRSTAT_TABLE   = OPT ;
  PSNID_INPUTS.LSIM           = LSIM ;
  PSNID_FITRES.NVAR           = 0;

  // allocate memory to store best fit residuals
  psnid_best_malloc_resids();

  if ( OPT == 0 ) { return ; }
  // ----------------------------------

  int  ID_TABLE     = PSNID_TABLE_ID;
  int  IFLAG        = 1 ;   // use regular SNANA variables
  char NAME_TABLE[] = PSNID_TABLE_NAME ; // "FITRES" ;
  char NAME_SNANABLOCK[] = "SNANAVAR" ;
  char NAME_SIMBLOCK[]   = "SIMVAR" ;

  PSNID_BEST_TABLEFILE_COMMENTS() ; // store comments for table

  SNTABLE_CREATE(ID_TABLE, NAME_TABLE, TEXTFORMAT ) ;

  // include standard SNANA-analysis variables
  init_table_snanavar__(&ID_TABLE, NAME_SNANABLOCK, &IFLAG, 
			strlen(NAME_SNANABLOCK) );

  // include SNANA-SIM variables
  if ( LSIM ) {
    init_table_simvar__(&ID_TABLE, NAME_SIMBLOCK, strlen(NAME_SIMBLOCK) );
  }

  // Aug 28 2017 KLUGE:
  // call function to set PSNID_NGRID[itype] so that table-init
  // knows which types to include.
  int nonia_types[PSNID_NONIA_MXTYPES] ;
  PSNID_MAXNL_NONIA  = 
    SNGRID_PSNID[TYPEINDX_NONIA_PSNID].NBIN[IPAR_GRIDGEN_SHAPEPAR];

  if ( PSNID_MAXNL_NONIA > PSNID_NONIA_MXTYPES ) {
    sprintf(c1err,"PSNID_MAXNL_NONIA=%d exceeds bound.", PSNID_MAXNL_NONIA);
    sprintf(c2err,"Check bound PSNID_NONIA_MXTYPES=%d",  PSNID_NONIA_MXTYPES);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  psnid_best_split_nonia_types(nonia_types, 0 );
  
  int DO_ADDCOL = 1; // add columns to FITRES table
  int DO_NN     = 0 ; // hard-wire until passed via &PSNIDNML
  psnid_best_define_TableVARNAMES(DO_ADDCOL,DO_NN);

  // check option to init residuals for table
  if ( OPT == 2 ) { psnid_best_define_TableRESIDS() ; }

  printf("\n"); fflush(stdout);
  return ;

} // end of void PSNID_BEST_INIT_SNTABLE


// =============================================
void  PSNID_BEST_TABLEFILE_COMMENTS(void) {

  // Oct 23 2014
  // define comments to write at top of ascii file ...
  // maybe later will aslo go into global header
  // for hbook or root file.

  int itype ;
  char comment[100], *KEY;
  // --------------- BEGIN ------------

  for ( itype=0; itype < PSNID_NTYPES; itype++ ) {   
    sprintf(comment,"ITYPE = %2d --> '%s' ", itype, PSNID_TYPE_NAME[itype]);
    STORE_TABLEFILE_COMMENT(comment);
  }

  
  // write unique keys for templates
  KEY = SNGRID_PSNID[TYPEINDX_SNIA_PSNID].UNIQUE_KEY ;
  if ( strlen(KEY) > 1 )  { 
    sprintf(comment,"GRIDKEY(TEMPLATES_SNIA)  = %s", KEY); 
    STORE_TABLEFILE_COMMENT(comment);
  }
	  
  KEY = SNGRID_PSNID[TYPEINDX_NONIA_PSNID].UNIQUE_KEY  ;
  if ( strlen(KEY) > 1 )  { 
    sprintf(comment, "GRIDKEY(TEMPLATES_NONIA) = %s", KEY ); 
    STORE_TABLEFILE_COMMENT(comment);
  }
	  

  sprintf(comment, "FILTLIST_FIT:  %s", PSNID_INPUTS.CFILTLIST );
  STORE_TABLEFILE_COMMENT(comment);

  sprintf(comment, "OPT_ZPRIOR: %d", PSNID_INPUTS.OPT_ZPRIOR );
  STORE_TABLEFILE_COMMENT(comment);

  sprintf(comment, "COLOR_BINS: %d  %6.3f  %6.3f ", 
	  PSNID_INPUTS.NCOLOR,  
	  PSNID_INPUTS.COLOR_MIN, PSNID_INPUTS.COLOR_MAX );
	  STORE_TABLEFILE_COMMENT(comment);

  sprintf(comment, "DMU_BINS: %d  %6.3f  %6.3f ", 
	  PSNID_INPUTS.NDMU,  PSNID_INPUTS.DMU_MIN, PSNID_INPUTS.DMU_MAX);
  STORE_TABLEFILE_COMMENT(comment) ;

  sprintf(comment, "MCMC_NSTEP: %d ", PSNID_INPUTS.MCMC_NSTEP );
  STORE_TABLEFILE_COMMENT(comment);


} // end of PSNID_BEST_TABLEFILE_COMMENTS



// ************************************
void  psnid_best_malloc_resids(void) {

  // Created Oct 2013 by RK.
  // Allocate memory for fit residuals vs. epoch.
  // Used for light curve plots and for fitres table/ntuple.

  int itype, MEMD, MEMI, MEMD2, MEMI2, MEMF ;
  
  MEMD  = MXEP_PSNID    * sizeof(double);
  MEMI  = MXEP_PSNID    * sizeof(int);
  MEMF  = MXEP_PSNID    * sizeof(float);

  MEMD2 = MXFILTINDX * sizeof(double);
  MEMI2 = MXFILTINDX * sizeof(int);


  for(itype=0; itype < PSNID_NTYPES; itype++ ) {
    RESIDS_PSNID_DOFIT[itype].TOBS          = (double*)malloc( MEMD ) ;
    RESIDS_PSNID_DOFIT[itype].REJECT        = (double*)malloc( MEMD ) ;
    RESIDS_PSNID_DOFIT[itype].CHI2          = (double*)malloc( MEMD ) ;    
    RESIDS_PSNID_DOFIT[itype].FITFLUX       = (double*)malloc( MEMD ) ;
    RESIDS_PSNID_DOFIT[itype].FITFLUX_ERR   = (double*)malloc( MEMD ) ;
  
    RESIDS_PSNID_DOFIT[itype].USE_FILT   = (int   *)malloc( MEMI2 ) ;
    RESIDS_PSNID_DOFIT[itype].NDOF_FILT  = (double*)malloc( MEMD2 ) ;
    RESIDS_PSNID_DOFIT[itype].CHI2_FILT  = (double*)malloc( MEMD2 ) ;
    RESIDS_PSNID_DOFIT[itype].PKMJD_FILT = (double*)malloc( MEMD2 ) ;    
  }
    
  // allocate table memory for best-fit only,
  // using 4-byte (F->float) instead of 8-byte. 
  F_RESIDS_PSNID_DOFIT.MJD           = (double*)malloc( MEMD ) ;
  F_RESIDS_PSNID_DOFIT.TOBS          = (float*)malloc( MEMF ) ;
  F_RESIDS_PSNID_DOFIT.TREST         = (float*)malloc( MEMF ) ;
  F_RESIDS_PSNID_DOFIT.CHI2          = (float*)malloc( MEMF ) ; 
  F_RESIDS_PSNID_DOFIT.PULL          = (float*)malloc( MEMF ) ; 
  F_RESIDS_PSNID_DOFIT.IFILTOBS      = (int  *)malloc( MEMI ) ; 
  F_RESIDS_PSNID_DOFIT.FITFLUX       = (float*)malloc( MEMF ) ;
  F_RESIDS_PSNID_DOFIT.FITFLUX_ERR   = (float*)malloc( MEMF ) ;
  F_RESIDS_PSNID_DOFIT.DATAFLUX      = (float*)malloc( MEMF ) ;
  F_RESIDS_PSNID_DOFIT.DATAFLUX_ERR  = (float*)malloc( MEMF ) ;
  
} // end of   psnid_best_malloc_resids



// *****************************************************
void psnid_best_define_TableVARNAMES(int DO_ADDCOL, int DO_NN ) {

  // Created Mar 2013  by RK
  // Define output variables for the following options,
  //
  //  DO_ADDCOL -> load columns for SNTABLE
  //  DO_NN    -> for nearest neighbor (sets USE4NN array)
  //
  // Here there is just one place to define a variable that
  // serves any of three functions.
  //
  // Jun 18, 2013: fix aweful bug: set NONIA_INDEX to SHAPEPAR
  //               instead of TMAX
  //

  int   ivar, ID, t, z, ISIA, IPAR, *USE4NN ;
  int   ifilt_obs, USE  ;

  int    *IPTR ;
  double *DPTR ;

  char 
    *ptrVARNAME
    ,TABLENAME[60]
    ,*CTYPE, cfilt[2]
    ,PREFIX[40]
    ,PREFIX_VARNAME[PSNID_NPARAM][40]   // list of variable name prefixes
    ,CBLOCK[] = "PSNIDVAR"              // for SNTABLE_ADDCOL only
    ,fnam[]   = "psnid_best_define_TableVARNAMES" 
    ;

  // --------------- BEGIN --------------

  ID   = PSNID_TABLE_ID ;

  printf("     Create BLOCK = %s for TABLE ID = %d \n", CBLOCK, ID);
  fflush(stdout);

  if ( DO_NN ) {
    for(ivar=0; ivar < PSNID_MXVAR_FITRES; ivar++ ) 
      { PSNID_FITRES.USE4NN[ivar] = IGNORE_for_NEARNBR ;  }
  }

  ivar = 0 ;

  // now the variables for both SNTABLE and legacy FITRES_DMPFILE

  ivar++ ;  ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
  sprintf(ptrVARNAME,"ZPRIOR");  
  //  psnid_dump_VARNAME(LDMP,ptrVARNAME);
  if ( DO_ADDCOL) {
    sprintf(TABLENAME,"%s:I", ptrVARNAME);
    IPTR = &PSNID_INPUTS.OPT_ZPRIOR ;
    SNTABLE_ADDCOL(ID, CBLOCK, IPTR, TABLENAME, 0 );
  }


  ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
  sprintf(ptrVARNAME,"ITYPE_BEST");   
  //  psnid_dump_VARNAME(LDMP,ptrVARNAME);
  if ( DO_ADDCOL) {
    sprintf(TABLENAME,"%s:I", ptrVARNAME);
    IPTR = &PSNID_BEST_RESULTS.FINAL_ITYPE[0] ;
    SNTABLE_ADDCOL(ID, CBLOCK, IPTR, TABLENAME, 1); 
  }

  // always book SIM_ITYPE; for data it's just set to -9.
  ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
  sprintf(ptrVARNAME,"SIM_ITYPE");   
  //  psnid_dump_VARNAME(LDMP,ptrVARNAME);
  if ( DO_NN ) { PSNID_FITRES.USE4NN[ivar] = 99; }
  if ( DO_ADDCOL) {
    sprintf(TABLENAME,"%s:I", ptrVARNAME);
    IPTR = &SIMVAR_PSNID.ITYPE ;
    SNTABLE_ADDCOL(ID, CBLOCK, IPTR, TABLENAME, 1);
  }

  

  z = 0; // sparse prior index is fixed

  // -----------------------------------------
  // hard-wire prefix for each PSNID_PARAM_XXX
  IPAR = 0;
  IPAR++ ; sprintf ( PREFIX_VARNAME[PSNID_PARAM_LOGZ],     "Z" );
  IPAR++ ; sprintf ( PREFIX_VARNAME[PSNID_PARAM_SHAPEPAR], "SHAPEPAR" );
  IPAR++ ; sprintf ( PREFIX_VARNAME[PSNID_PARAM_COLORPAR], "COLORPAR" );
  IPAR++ ; sprintf ( PREFIX_VARNAME[PSNID_PARAM_COLORLAW], "COLORLAW" );
  IPAR++ ; sprintf ( PREFIX_VARNAME[PSNID_PARAM_TMAX],     "TMAX" );
  IPAR++ ; sprintf ( PREFIX_VARNAME[PSNID_PARAM_DMU],      "DMU" );
  if ( IPAR != PSNID_NPARAM ) {
    sprintf(c1err,"PREFIX_VARNAME set for %d params", IPAR);
    sprintf(c2err,"but should have been PSNID_NPARAM = %d times.", 
	    PSNID_NPARAM);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  // ------ loop over types (Ia, Ibc, II ...)

  for ( t=0 ; t < PSNID_NTYPES ; t++ ) {

    if ( PSNID_NGRID[t] == 0 ) { continue ; } // Aug 28 2017
    
    CTYPE  = PSNID_ITYPE_STRING[t] ;
    ISIA   = (t == PSNID_ITYPE_SNIA ) ;

    for ( IPAR = 0; IPAR < PSNID_NPARAM ; IPAR++ ) {
      ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;

      sprintf(PREFIX,"%s", PREFIX_VARNAME[IPAR] );

      // check for special prefix for NONIA index
      if ( IPAR == PSNID_PARAM_SHAPEPAR  &&  ISIA == 0 )
	{ sprintf(PREFIX, "NONIA_INDEX" ); }


      sprintf(ptrVARNAME,"%s_%s", PREFIX, CTYPE);   
      //      psnid_dump_VARNAME(LDMP,ptrVARNAME);

      if ( DO_NN && ISIA ) {
	USE4NN  = &PSNID_FITRES.USE4NN[ivar] ;
	if ( IPAR == PSNID_PARAM_LOGZ     )    { *USE4NN = IPAR; }
	if ( IPAR == PSNID_PARAM_COLORPAR )    { *USE4NN = IPAR; }
	if ( IPAR == PSNID_PARAM_SHAPEPAR )    { *USE4NN = IPAR; }
      }

      if ( DO_ADDCOL) {
	sprintf(TABLENAME,"%s:D", ptrVARNAME);
	DPTR  = &PSNID_BEST_RESULTS.FINAL_PARVAL[z][t][IPAR] ;
	SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 1 );
      }

    } // end of IPAR loop


    // ---

    // add Tobsmin and Tobsmax (RK, Jan 2014)
    ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
    sprintf(ptrVARNAME,"TOBSMIN_%s", CTYPE);   
    //    psnid_dump_VARNAME(LDMP,ptrVARNAME);
    if ( DO_ADDCOL) {
      sprintf(TABLENAME,"%s:D", ptrVARNAME);
      DPTR  = &PSNID_BEST_RESULTS.TOBSMIN[z][t] ;
      SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 1 );
    }

    ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
    sprintf(ptrVARNAME,"TOBSMAX_%s", CTYPE);   
    //    psnid_dump_VARNAME(LDMP,ptrVARNAME);
    if ( DO_ADDCOL) {
      sprintf(TABLENAME,"%s:D", ptrVARNAME);
      DPTR  = &PSNID_BEST_RESULTS.TOBSMAX[z][t] ;
      SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 1 );
    }

    // CHI2
    ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
    sprintf(ptrVARNAME,"CHI2_%s", CTYPE);   
    //    psnid_dump_VARNAME(LDMP,ptrVARNAME);
    if ( DO_ADDCOL) {
      sprintf(TABLENAME,"%s:D", ptrVARNAME);
      DPTR  = &PSNID_BEST_RESULTS.MINCHISQ[z][t] ;
      SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 1 );
    }

    // NPT
    ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
    sprintf(ptrVARNAME,"NPT_%s", CTYPE);   
    //    psnid_dump_VARNAME(LDMP,ptrVARNAME);
    if ( DO_ADDCOL) {
      sprintf(TABLENAME,"%s:I", ptrVARNAME);
      IPTR  = &PSNID_BEST_RESULTS.NGOOD[z][t] ;
      SNTABLE_ADDCOL(ID, CBLOCK, IPTR, TABLENAME, 1 );    
    }

    // FITPROB
    ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
    sprintf(ptrVARNAME,"FITPROB_%s", CTYPE);   
    //    psnid_dump_VARNAME(LDMP,ptrVARNAME);
    if ( DO_ADDCOL) {
      sprintf(TABLENAME,"%s:D", ptrVARNAME);
      DPTR  = &PSNID_BEST_RESULTS.FITPROB[z][t] ;
      SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 1 );
    }

    // PBAYES
    ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
    sprintf(ptrVARNAME,"PBAYES_%s", CTYPE);   
    //    psnid_dump_VARNAME(LDMP,ptrVARNAME);
    if ( DO_ADDCOL) {
      sprintf(TABLENAME,"%s:D", ptrVARNAME);
      DPTR  = &PSNID_BEST_RESULTS.PBAYESIAN[z][t] ;
      SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 1 );    
    }

    // LIGHTCURVE_QUALITY
    ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
    sprintf(ptrVARNAME,"LCQ_%s", CTYPE);   
    //    psnid_dump_VARNAME(LDMP,ptrVARNAME);
    if ( DO_ADDCOL) {
      sprintf(TABLENAME,"%s:I", ptrVARNAME);
      IPTR  = &PSNID_BEST_RESULTS.LIGHTCURVE_QUALITY[z][t] ;
      SNTABLE_ADDCOL(ID, CBLOCK, IPTR, TABLENAME, 1 );    
    }

    // include errors for Ia only
    if ( t == PSNID_ITYPE_SNIA ) {
      for ( IPAR = 0; IPAR < PSNID_NPARAM ; IPAR++ ) {
	ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
	sprintf(ptrVARNAME,"%se_%s", PREFIX_VARNAME[IPAR], CTYPE);   
	//	psnid_dump_VARNAME(LDMP,ptrVARNAME);
	if ( DO_ADDCOL) {
	  sprintf(TABLENAME,"%s:D", ptrVARNAME);
	  DPTR  = &PSNID_BEST_RESULTS.FINAL_PARERR[z][t][IPAR] ;
	  SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 1 );   
	}
      } // end of IPAR loop
    } // end of Ia if-block


    // check option to store peak mags (Nov 1 2013 - RK)

    for(ifilt_obs=0; ifilt_obs < MXFILTINDX; ifilt_obs++ ) {
      USE = PSNID_INPUTS.USEFILT_PEAKMAG_STORE[ifilt_obs] ;
      if ( USE == 0 ) { continue ; }
	
      sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );
      ivar++ ; ptrVARNAME =  PSNID_FITRES.VARNAMES[ivar] ;
      sprintf(ptrVARNAME,"PEAKMAG_%s_%s", cfilt, CTYPE);
      //      psnid_dump_VARNAME(LDMP,ptrVARNAME);       
      
      if ( DO_ADDCOL) {
	sprintf(TABLENAME,"%s:D", ptrVARNAME);
	DPTR  = &PSNID_BEST_RESULTS.FINAL_PEAKMAG[t][ifilt_obs];
	SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 1 );
      }
	
    } // end ifilt_obs


    // make sure we don't exceed array bound.
    if ( ivar >= PSNID_MXVAR_FITRES ) {
      sprintf(c1err,"NVAR(fitres) = %d exceeds bound of", ivar );
      sprintf(c2err,"PSNID_MXVAR_FITRES = %d", PSNID_MXVAR_FITRES ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  } // end loop over types

  
  PSNID_FITRES.NVAR = ivar ; // store as global

  return ;

} // end of psnid_best_define_TableVARNAMES


// *********************************************
void psnid_best_define_TableRESIDS(void) {

  // Created Oct 2013
  // init light curve fit residuals for best-fit type

  int ID, *IPTR ;
  float  *FPTR ;
  double *DPTR ;
  char CBLOCK[] = "RESIDS" ;
  char TABLENAME[60];

  // ---------- BEGIN --------------

  ID   = PSNID_TABLE_ID ;

  printf("     Create BLOCK = %s for TABLE ID = %d \n", CBLOCK, ID);
  fflush(stdout);

  //  sprintf(TABLENAME,"NFIT[0,600]:I") ;
  sprintf(TABLENAME,"NFIT[600]:I") ;  // May 11 2014: use C-like declaration
  IPTR  = &F_RESIDS_PSNID_DOFIT.NFIT ;
  SNTABLE_ADDCOL(ID, CBLOCK, IPTR, TABLENAME, 0 );
 
  sprintf(TABLENAME,"IFILTOBS(NFIT):I") ;
  IPTR  = F_RESIDS_PSNID_DOFIT.IFILTOBS ;
  SNTABLE_ADDCOL(ID, CBLOCK, IPTR, TABLENAME, 0 ); 

  sprintf(TABLENAME,"TOBS(NFIT):F");
  FPTR  = F_RESIDS_PSNID_DOFIT.TOBS ;
  SNTABLE_ADDCOL(ID, CBLOCK, FPTR, TABLENAME, 0 );

  sprintf(TABLENAME,"TREST(NFIT):F");
  FPTR  = F_RESIDS_PSNID_DOFIT.TREST ;
  SNTABLE_ADDCOL(ID, CBLOCK, FPTR, TABLENAME, 0);

  sprintf(TABLENAME,"MJD(NFIT):D");
  DPTR  = F_RESIDS_PSNID_DOFIT.MJD ;
  SNTABLE_ADDCOL(ID, CBLOCK, DPTR, TABLENAME, 0 );
  
  sprintf(TABLENAME,"DATAFLUX(NFIT):F");
  FPTR  = F_RESIDS_PSNID_DOFIT.DATAFLUX ;
  SNTABLE_ADDCOL(ID, CBLOCK, FPTR, TABLENAME, 0 );

  sprintf(TABLENAME,"DATAFLUX_ERR(NFIT):F");
  FPTR  = F_RESIDS_PSNID_DOFIT.DATAFLUX_ERR ;
  SNTABLE_ADDCOL(ID, CBLOCK, FPTR, TABLENAME, 0 );

  sprintf(TABLENAME,"FITFLUX(NFIT):F");
  FPTR  = F_RESIDS_PSNID_DOFIT.FITFLUX ;
  SNTABLE_ADDCOL(ID, CBLOCK, FPTR, TABLENAME, 0 );

  sprintf(TABLENAME,"FITFLUX_ERR(NFIT):F");
  FPTR  = F_RESIDS_PSNID_DOFIT.FITFLUX_ERR ;
  SNTABLE_ADDCOL(ID, CBLOCK, FPTR, TABLENAME, 0 );

  sprintf(TABLENAME,"CHI2(NFIT):F");
  FPTR  = F_RESIDS_PSNID_DOFIT.CHI2 ;
  SNTABLE_ADDCOL(ID, CBLOCK, FPTR, TABLENAME, 0 );

  sprintf(TABLENAME,"FITPULL(NFIT):F");
  FPTR  = F_RESIDS_PSNID_DOFIT.PULL ;
  SNTABLE_ADDCOL(ID, CBLOCK, FPTR, TABLENAME, 0 );

  return ;

} // end of  psnid_best_define_TableRESIDS


// *****************************************************
void PSNID_BEST_UPDATE_OUTPUT(char *CCID) {

  int itype ;
  //  char fnam[] = "PSNID_BEST_UPDATE_OUTPUT" ;
  
  // ---------------- BEGIN --------------

  SIMVAR_PSNID.ITYPE = -9 ; // init to unknown ITYPE

  // for simulation determine SIM_ITYPE using PSNID indices 
  // (not snana indices)
  if ( PSNID_INPUTS.LSIM ) {

    char SIMNAME_TYPE[12];
    get_simname_type__(SIMNAME_TYPE, 12);

    for ( itype=0; itype < PSNID_NTYPES; itype++ ) { 
      if ( psnid_strTypeMatch(itype,SIMNAME_TYPE, 0) == 1 )  
	{ SIMVAR_PSNID.ITYPE=itype ; }
    } 
  }  // end of LSIM block

  
  if ( PSNID_INPUTS.WRSTAT_TABLE   ) 
    { SNTABLE_FILL(PSNID_TABLE_ID); }


}  // end of PSNID_BEST_UPDATE_OUTPUT



// ===========================================
void PSNID_BEST_SUMMARY(void) {

  // optional summary.

} // end of PSNID_BEST_SUMMARY



// ============================================
void psnid_best_nearnbr(char *CCID) {

  double DVAL ;
  int z, t, IPAR, ivar, USE4NN ;
  char *VARNAME_PSNID;
  //  char fnam[] = "psnid_best_nearnbr" ;

  // --------------- BEGIN ----------

  if ( PSNID_NEARNBR.NVAR <= 0 ) { return ; }

  z = 0 ; t = PSNID_ITYPE_SNIA ; 

  // search for variables with USE4NN flag

  for ( ivar=0; ivar <= PSNID_FITRES.NVAR; ivar++ ) {
    USE4NN = PSNID_FITRES.USE4NN[ivar] ;
    //    if ( USE4NN < 0  || USE4NN == 99 ) { continue ; }
    if ( USE4NN < 0 ) { continue ; }
    IPAR    = USE4NN ;
    DVAL    = PSNID_BEST_RESULTS.FINAL_PARVAL[z][t][IPAR] ;
    VARNAME_PSNID  = PSNID_FITRES.VARNAMES[ivar] ;
    NEARNBR_LOADVAL(CCID, VARNAME_PSNID, DVAL ) ; // see sntools_nearnbr.c
  } 

  NEARNBR_APPLY(CCID); // do the NN analysis

} // end of  psnid_best_nearnbr
