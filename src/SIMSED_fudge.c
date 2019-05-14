/**************************************
 Created July 2, 2010 by R.Kessler

 The two main types of SIMSED fudges are
 1. Apply arbitrary color law to a set of explosion model SEDs.
    The SEDs must be in SNANA format and reside in 
     $SNDATA_ROOT/models/SIMSED

 2. Apply arbitrary mag-smearing model.

 To apply this program to the Kasen/Flash models, a conversion
 program such as Benedikt's SIMSED_prep must run first
 to generate the SEDs in SNANA format.

 Beware that the program's  SIMSED output can be very large
 when the color law is added (option 1 above).
 For example, if the initial SIMSED set is a function of
 two variables (such as MNi56 and COSANGLE), the output
 SIMSED set will be a function of three variables, where 
 the 3rd variable is the user choice of RV, SALT2 color, 
 etc ... If 10 'color' bins are requested, the output 
 SIMSED set will be x10 larger !  You can use input keys
 EPOCH_RANGE and LAMBDA_RANGE to trim the SEDs.

 In the new SIMSED version (defined by user input keyword
 SIMSED_VERSION_OUTPUT), the new SED.INFO file is udpated 
 to include the extra color variable.

 To use PSNID (Sako et al. 2008, Appendix C ) bolometric
 luminosity normalization, set user input keyword
 PSNID_NORM: 1.0 . Default setting is 0.0 (off). 

 To print "debug" version of color & stretch smeared
 SIMSED templates, set user input keyword
 DEBUG_SIMSED: 1.0 . Default setting is 0.0 (off).
 template file format is then:
 <day><lam><flux>(<fluxerr>)<origflux> 
      <polyfudge><cfudge><sfudge><shapecor>
      (<peak_lumin><psnid_norm>).
 Columns in parentheses are output only if the appropriate
 flags are set. 

 To add a "NEW" color law option, modify this code as follows: 
   - increment NCHOICE_COLORFUN
   - define INDEX_[NEW]
   - define new COLORFUN_LIST  NAME, NPAR and PARNAMES 
     in init_colorfuns().

 The stretch function converts the user shape parameter to
 stretch, for use in SED/light curve stretching.  

 To add a "NEW" STRETCHFUN option, modify this code as follows:
   - increment NCHOICE_STRETCHFUN
   - define INDEX_STRETCHFUN_[NEW]
   - define new STRETCHFUN_LIST NAME, NPAR and PARNAMES
     in init_stretchfuns().

 For a given shape parameter value, the shapeFudge function 
 alters SED flux to create the shape-brightness relation. 

 To add a "NEW" shapeFudge option, modify this code as follows:
   - increment NCHOICE_SHFUDGEFUN
   - define INDEX_SHAPE_[NEW]
   - define new SHFUDGEFUN_LIST NAME, NPAR and PARNAMES
     in init_shFudgefuns().

 To add a new MAG-SMEAR model 'XXX', modify this code as follows:

   - define new  element of struct SMEARMODEL_ID
   - define new SMEARMODEL values in init_SMEAROPT(char *modelName)
   - defgin new function SMEARMAG_XXX(lam,day)
   - call new function from sedFudge_smear()



    HISTORY
  ~~~~~~~~~~

 Aug 02, 2010:
  - write color fudge (flux-ratio) to  colorFudge.dat
  - remove C_OFF parameters for SALT2colorlaw0
  - change colorlaw names: 
      SALT2_Guy07 -> SALT2colorlaw0
      SALT2_Guy10 -> SALT2colorlaw1


 Aug 28, 2010:
  update to allow for flux-error column
  See FLUX_ERRFLAG.

 Sep 7, 2010:  pass fudged sed to sedSmooth to allow smoothing algorithm.
               New function sedWrite to write final smoothed sed.

 Sep 13, 2010 : call sedSmooth after reading SED but before applying fudge.
                No point is smoothing after the analytic fudge.
                Also made new 'sedRead(ised)' function.

 Sep 15, 2010: install Rahul's sedSmooth function.
               See OPT_SMOOTH_TOPHAT and OPT_SMOOTH_GAUSS

 Sep 16, 2010: sedFudge.dat renamed back to colorFudge.dat
               Add new options to fudge shape parameter ...
               analogous to color fudge. For now only
               SALT2 x1 is treated and only mags are effected;
               there is no correction (yet?) to spectral surface.

 Sep 22, 2010: translate SHAPEPAR into stretch and apply stretch-fudge
               as well as the shape-mag fudge. See new function
               'stretchFlux()'.

               New utility 'update_SEDINFO_file()' to isolate this stuff
               from sedFudge().

               New user-inut CORRECT_DM15: <parname> 

 Sep 26, 2010: compute stretch relative to lambda-dependent Tpeak.
               Treak = TPEAK_LAMPOLY(lambda) is a polynomial function 
               of lambda based on user  INPUTS.LAMPOLY_TPEAK.

 Feb 19, 2011: fix  colorlaw1 option to work in colorFudge().

 Mar 04, 2011: print colorlaw parameters in SED.INFO
               so that JLM can extract BETA

 Mar 11, 2011: new function check_for_S2x1 to check for x1 in the
               input SIMSED; if there then scale output SED fluxes
                by 10**[.4*x1*S2ALPHA]

 May 6, 2011: MXCBIN = 100 (was 40)

 Oct 02, 2011: add new MAG-SMEAR options; see input key 
                 SMEAROPT: <modelName> <model param values ...>

 Nov 18, 2011: insert sedSynMag to compute mags and colors to check
               color smearing  for the smearing options.
               Just a dummy calc for now to get the structure in place.

 Dec 20, 2011:  add SIMSED.KRW09  SMEAROPT.

 Dec 28, 2011: declare NPAR_SEDINFO_OUPUT to fix dumb bug.

 April 02, 2012: better abort trap on invalid day in interp_SEDFLUX_KRW09()

 May 17, 2012: fix bug related to CORRECT_PEAKMB and CORRECT_DM15.
               See new function check_for_CORRECT(void)

 Jun 7, 2012: add G10 and C11 options for smear models.
              Use functions in genSmear_models.c

 Jun 28, 2012:  New 'Blondin8' option for KRW09 models to read the 8-summed
                models. Also writes out MB and DM15.

 Aug 31, 2012: in init_SMEARMODEL() function, 
                print smearopt summary with each parameter.

 Nov 27, 2012: fix a few compilation errors found by Ishida
               (related to printing error messages)

 Jun 17, 2013: Add hooks for user-input lcstretch and 
               shapeFudge functions. Add optional PSNID 
               bolometric luminosity normalization and
               shapeFudge sedWrite debug routines. 
               grep "FLAG_PSNID_NORM" to see bol norm code. 
               grep "FLAG_DEBUG_SIMSED" to see debug code. 
               (JLM)

 Aug 25 2013 RK -
    replace local Btransfun() with call to filterTrans_BessB() 
    in genmag_SEDtools.c. Same function, but it is not a utility.

 Jan 21, 2014 RK - inside init_SMEARMODEL() call get_NRAN_genSmear

 Mar 11 2014: for hard-wired Bessell filter location, Bessell90/ -> 
              Bessell90/Bessell90_K09/

 Mar 13 2014: add VCR model, and new option to write SIMSED with 
              extra dimension using SMEARPAR_GRID info.              

 Mar 16, 2014: remove PSNID_NORM stuff since S11DM15 is a real
               model in the simulation.

 Mar 22 2014: new input option "SALT2_VERSION_INPUT: SALT2"
              to convert into SIMSED with using the intermediate
              paw macro simsed.kumac . Created 2D grid in both
              x1 and color.

 Aug 2017: remove SEDMODEL.FLUX_ERRFLAG and replace with SEDMODEL.OPTMASK
           NOT TESTED.

**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_genSmear.h"
#include "MWgaldust.h"

#include "genmag_SIMSED.h"
#include "genmag_SEDtools.h"
#include "genmag_SALT2.h" 

// =========== DECLARE FUNCTIONS ==================
void  parse_args(int argc, char **argv);
void  init_colorfuns(void);
void  init_shFudgefuns(void);
void  init_stretchfuns(void);
void  read_input(void);
void  summarize_input(void);
void  init_colorpar(void);
void  init_shapepar(void);
int   init_SMEAROPT(char *modelName); 
void  init_SALT2_fudge(char *version) ;
void  set_SIMSED_PATH(void);
void  mkdir_SIMSED_VERSION_OUTPUT(void);
void  wr_colorFudge(void);
void  open_SEDINFO_file(void);
void  wrhead_SEDINFO_color(void);
void  wrhead_SEDINFO_smear(void);

void  update_SEDINFO_file(FILE *fp, char *sed_outFile, int ised, 
			  double cval, double sval,
			  double MB0cor, double DM15cor );

void  sedRead(int ised);
void  sedSmooth_driver(int ised);
void  sedFudge_color(int ised);
void  sedFudge_smear(int ised, int ismear);
void  load_TEMP_SEDFUDGE_SMEAR(int iday, int ilam, int ismear);
void  sedWrite(FILE *fp);

void  sedSynMag(int ised, int ifilt); 

void   init_SMEARMODEL(void) ;
void   load_SMEARPAR_GRID(char *parName, int indx,
			  double parMin, double parMax, double parBin );
double SMEARMAG_COHERENT(double lam, double day);
double SMEARMAG_SINX(double lam, double day);
double SMEARMAG_KRW09(double lam, double day);
void   init_SMEARMAG_KRW09(void);
void   read_SEDINFO_KRW09(void) ;
void   read_SED_KRW09(int imodel, int icos);
void   malloc_KRW09(int imodel) ;
void   smooth_SED_KRW09(int imodel, int icos) ;
double interp_SEDFLUX_KRW09(int imodel, int icos, double lam, double day) ;

//void   mktup_KRW09(void);
void   fluxavg_KRW09(void);


void  print_colorfun_options(void);
void  print_stretchfun_options(void);
void  print_shfudgefun_options(void);

void  init_REMOVE(void);
void  check_for_S2x1(void);
void  check_for_CORRECT(void);

double TPEAK_LAMPOLY(double lam);
double FLUXFUDGE_LAMPOLY(double lam);

double Btransfun(double lam);
double colorFudge(double lambda, double colorVal);
double shapeFudge(double lambda, double shapeVal);
double lcstretch(double x1); // return light curve stretch from shape par
double stretchFlux(double str, int iday, int ilam);

void sedSmooth_nbravg(int NLAM, double *lambda,
		      double *FLUX_IN, double *FLUXERR_IN,
		      double *FLUX_SMOOTH, double *FLUXERR_SMOOTH );

void sedSmooth_dump(int ised, double Trest, int NLAM, double *lam, 
		    double *flux_in,  double *fluxerr_in,
		    double *flux_out, double *fluxerr_out );

void make_param_ntuple(void);

void get_intBins(int imodel, int icos, double lam, double day, 
		 int *ilam, int *iday, char *abortComment) ;

int ALLMODELS_VALID(double DAY, double LAM) ;

// Rahul's smoothing functions

int sedSmooth( int opt, double smoothpar, int Nlam, double *lambda,
		double *flux_in, double *fluxerr_in, 
		double *flux_smooth, double *fluxerr_smooth );

double *setfilter(int l, int opt);
double *setbox(int l, int j, int nBins, double *fl);
void printsedSmooth_usage();


// =========== DECLARE GLOBALS =================

#define MXPARAM   20   // max no. of parameters for color function (RV, etc ..)
#define MXCBIN   100   // max no. of color bins
#define MXREMOVE  50   // max no SED bins to remove
#define MXLAMPOLY  4

#define MXLAM_SMOOTH  MXBIN_LAMSED_SEDMODEL 
#define MXLAM_FILT    MXBIN_LAMSED_SEDMODEL/10  // for synthetic mag calc
#define OPT_NOSMOOTH       0
#define OPT_SMOOTH_TOPHAT  1  // tophat filtetr
#define OPT_SMOOTH_GAUSS   2  // Guass filter
#define OPT_SMOOTH_NBRAVG  9  // average nearest neighbor (silly test)
#define MXDUMP_SMOOTH  50   // max no. of SEDs to DUMP smoothing results
#define S2ALPHA  0.11 
#define S2BETA   3.20
#define MBOFF    10.635 // see mBoff_SALT2 in genmag_SALT2.c

// Log of peak luminosity (erg/s) of standard SN Ia                                  
#define LUMIN_STANDARD 43.2638

int NPAR_SEDINFO_OUTPUT; // Number of params written to SED.INFO file

// char SNDATA_ROOT[200];

// define structure with info about all smear models

// #define MXFILE_SMEARMODEL 100 // if SMEARMODEL uses files
#define MXPAR_SMEARMODEL 10   // max number of params per smearmodel
#define MXDEF_SMEARMODEL 20   // max number of smear-models
#define MAGSMEAR_MIN    -3.0  // don't allow magSmear less than this
#define MAGSMEAR_MAX    +3.0  // don't allow magSmear above this
#define MAGSMEAR_ABORT  10.0

#define MXBIN_SMEARPAR 500   // for smearpar 

int     NDEF_SMEARMODEL ;    // number of defined smear models
double  GAURANLIST_SMEARMODEL[MXDEF_SMEARMODEL] ;
double  FLATRANLIST_SMEARMODEL[MXDEF_SMEARMODEL] ;
double  SYNMAGDIF_SMEARMODEL[MXPAR_SMEARMODEL] ;

char INFO_FILENAME2[] = "SED_PARAM.DAT"; // fitres-format for ntuple

struct SMEARMODEL_DEF {
  int  NPAR ;    
  int  NGAURAN ;    // number of Gaussian random numbers per SED sequence
  int  NFLATRAN ;   // number of flat-randoms [0-1] per sequence
  char NAME[80];
  char COMMENT[80];

  int  NFILTMAG ;   // number of synthetic mags to compute
  char FILT[20] ;
  
} SMEARMODEL_DEF[MXDEF_SMEARMODEL] ;


// define physical [optional] smearing parameter; e.g., v_Si for VCR model
// Note that number of SEDs written is multiplied by NBIN_SMEARPAR
struct SMEARPAR_GRID {
  int    NBIN ;
  double VALUE[MXBIN_SMEARPAR];
  char   NAME[40];
} SMEARPAR_GRID ;


struct SMEARMODEL_ID {
  int NOSMEAR ;
  int COHERENT ;
  int SINX;
  int KRW09 ;
  int G10 ;
  int C11 ;
  int VCR ;
  int USE_genSmear ; // T => use model in genSmear_models.c
} SMEARMODEL_ID ;


// specialized structures for specific models.

#define MXMODEL_KRW09 10
#define MXCOS_KRW09   40
struct KRW09 {
  int     NMODELS ; // number of models defined by DC and IGNIT
  double  DC[MXMODEL_KRW09];
  double  IGNIT[MXMODEL_KRW09] ; // vs. NMODELS


  int     SUPER8FLAG; // -> use all SEDs in KRW09_Blondin8

  int     NPAR ;
  int     IPAR_COSANGLE, IPAR_DC, IPAR_IGNIT, IPAR_DM15, IPAR_MB ;  
  int     N_COSANGLE ;
  char    SEDFILE_LIST[MXMODEL_KRW09][MXCOS_KRW09][MXPATHLEN] ;
  double  COSANGLE_LIST[MXCOS_KRW09] ;
  double  COSANGLE_MIN, COSANGLE_MAX;
  double  DM15_LIST[MXMODEL_KRW09][MXCOS_KRW09] ;
  double  MB_LIST[MXMODEL_KRW09][MXCOS_KRW09] ;
  int     ICOS_NEAR90DEG[2] ; // just below and above 0
  int     FLUX_ERRFLAG ;

  char    DIR[MXPATHLEN];

  int     NDAY[MXMODEL_KRW09][MXCOS_KRW09] ;
  int     NLAM[MXMODEL_KRW09][MXCOS_KRW09] ;
  double  DAYSTEP[MXMODEL_KRW09][MXCOS_KRW09] ;
  double  LAMSTEP[MXMODEL_KRW09][MXCOS_KRW09] ;
  double  **DAY[MXMODEL_KRW09] ;
  double  **LAM[MXMODEL_KRW09] ;
  double  **FLUX[MXMODEL_KRW09] ;
  double  **FLUXERR[MXMODEL_KRW09] ;

  int MXBIN_FLUX; // max number of flux bins (for malloc)

  double SEDRANGE_DAY[2];
  double SEDRANGE_LAM[2];

  
  // temp values used by SMEARMAG_KRW09
  int    TEMP_ICOS_NEAR[2] ;
  double TEMP_RANFLAT ;

  // current random values of COSANGLE and IMODEL
  double  COSANGLE_RAN, DM15_RAN, MB_RAN ;
  int     IMODEL_RAN ;

} KRW09 ;

// ------------------
// user INPUTS structure
struct INPUTS {
  char inputFile[200];

  int  MODELTYPE ;  
  char SALT2_VERSION_INPUT[200];
  char SIMSED_VERSION_INPUT[200];
  char SIMSED_VERSION_OUTPUT[200];

  char   SHAPEPAR_NAME[40] ;  // x1, DELTA,stretch ...
  double SHAPEPAR_RANGE[2] ;
  double SHAPEPAR_BINSIZE ;
  int    NBIN_SHAPEPAR ;
  double SHAPEPAR_VALUE[MXCBIN]; // shape value in each bin
  double SHAPEFUN_PARVAL[MXPARAM];

  char   COLORPAR_NAME[40] ;  // AV, c, ...
  double COLORPAR_RANGE[2] ;
  double COLORPAR_BINSIZE ;
  int    NBIN_COLORPAR;
  double COLORPAR_VALUE[MXCBIN]; // color value in each bin

  char   COLORFUN_NAME[40];    // SALT2colorlaw0, CCM89, etc ...
  double COLORFUN_PARVAL[MXPARAM];
  int    NPAR_COLORFUN;

  char   SHFUDGEFUN_NAME[40];  // DEFAULT, SAKO08, etc..
  double SHFUDGEFUN_PARVAL[MXPARAM];
  int    NPAR_SHFUDGEFUN;

  char   STRETCHFUN_NAME[40];  // DEFAULT, SAKO08, etc..
  double STRETCHFUN_PARVAL[MXPARAM];
  int    NPAR_STRETCHFUN;

  int    NPAR_SHAPEFUN; // duplicate NPAR_STRETCHFUN

  double SIMSED_TREST_RANGE[2] ; // range to read SED (days)
  double SIMSED_LAM_RANGE[2] ;   // range to read SED (A)

  char   CORRECT_PEAKMB[40];  // correct original peak MB
  char   CORRECT_DM15[40];    // correct original DM15
  int    IPAR_PEAKMB;         // parameter index of peak MB
  int    IPAR_DM15;           // parameter index of DM15

  double LAMPOLY_TPEAK[MXLAMPOLY];   // Tpeak = polyfun(lambda) for stretch
  double LAMPOLY_FLUX[MXLAMPOLY];    // 3rd order poly fudge on flux

  // option to remove SED(s) based on parameter name & value
  int    NSEDPAR_REMOVE;      // number of SEDs to remove
  char   SEDPARNAME_REMOVE[MXREMOVE][20];
  double SEDPARVAL_REMOVE[MXREMOVE];
  int    IPAR_REMOVE[MXREMOVE];

  int    OPT_SMOOTH; // 0 = no smoothing
  double LAMSCALE_SMOOTH ;
  char   COMMENT_SMOOTH[100];

  int    NDUMP_SMOOTH;  // number of SEDs to dump 
  char   DUMP_SMOOTH_SEDNAME[MXDUMP_SMOOTH][200];
  double DUMP_SMOOTH_TREST[MXDUMP_SMOOTH];

  double DUMP_LAMBDA; 

  // for mag-smearing only
  int     INDEX_SMEARMODEL ;
  char    SMEARMODEL_PARSTRING[MXPAR_SMEARMODEL][60] ;
  double  SMEARMODEL_PARVAL[MXPAR_SMEARMODEL] ;

  // alternate SIMSED modelpath 
  char SIMSED_MODELPATH[200];

} INPUTS ;

char SIMSED_PATHMODEL_INPUT[200];
char SIMSED_PATHMODEL_OUTPUT[200];


#define NCHOICE_COLORFUN      4
#define INDEX_SALT2colorlaw0  1
#define INDEX_SALT2colorlaw1  2
#define INDEX_CCM89           3   // Cardelli, Clayton Mathis 1989
#define INDEX_CCM94           4   // O'Donnel update to CCM89 law
#define NULLVAL -9999999.

#define NCHOICE_STRETCHFUN              3
#define INDEX_STRETCHFUN_NONE           1  // no conversion to stretch
#define INDEX_STRETCHFUN_DEFAULT        2  // SALT2x1 -> stretch
#define INDEX_STRETCHFUN_SAKO08         3  // Sako et. al, 2008 Appx C

#define NCHOICE_SHFUDGEFUN           2
#define INDEX_SHFUDGE_DEFAULT          1
#define INDEX_SHFUDGE_SAKO08           2  // Sako et. al, 2008 Appx C

int INDEX_COLORFUN ;
int INDEX_STRETCHFUN ; 
int INDEX_SHFUDGEFUN ;
int NSED_REMOVED[MXREMOVE];  // actual number of SED bins removed 
int SKIPSED[MXSEDMODEL] ;    // logical flag to skip SED
int IPAR_S2x1 ;

double FLAG_DEBUG_SIMSED ; // if FLAG == 0, print all fudge factors

struct COLORFUN_LIST {
  char NAME[40] ;
  int  NPAR ;
  char PARNAMES[MXPARAM+1][40];
};

struct COLORFUN_LIST COLORFUN_LIST[NCHOICE_COLORFUN+1];
struct COLORFUN_LIST SHFUDGEFUN_LIST[NCHOICE_SHFUDGEFUN+1];
struct COLORFUN_LIST STRETCHFUN_LIST[NCHOICE_STRETCHFUN+1];

FILE *fp_sedinfo ;
FILE *fp_sedinfo2 ;

int NSED_OUTPUT;

struct TEMP_SEDSMOOTH {
  double FLUX[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL];
  double FLUXERR[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL];
} TEMP_SEDSMOOTH ;



struct TEMP_SEDFUDGE {
  double FLUX[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL];
  double FLUXERR[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL];

  double ORIG_FLUX[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL]; // FLAG_DEBUG_SIMSED
  double POLYFUDGE[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL]; // FLAG_DEBUG_SIMSED
  double CFUDGE[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL]; // FLAG_DEBUG_SIMSED
  double SFUDGE[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL]; // FLAG_DEBUG_SIMSED
  double SHAPECOR[MXBIN_DAYSED_SEDMODEL][MXBIN_LAMSED_SEDMODEL]; // FLAG_DEBUG_SIMSED

} TEMP_SEDFUDGE ;


char colorFudge_datFile[60] = "colorFudge.dat" ;
char shapeFudge_datFile[60] = "shapeFudge.dat" ;


char CTAGNTUP[20][8], *ptr_CTAGNTUP[20] ;

// Rahul's stuff for smoothing.

#define pi 3.1415926535
double  BINSIZE_LAMBDA ;
int     NCALL_SMOOTH;


// ====================================
int main(int argc, char **argv) {

  char fnam[] = "main" ;
  int ised, ismear;
  // ------------ BEGIN ------------

  set_EXIT_ERRCODE(EXIT_ERRCODE_UNKNOWN);

  sprintf(PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );

  sprintf(BANNER,"Begin execution of SIMSED_fudge" );
  print_banner(BANNER);

  parse_args(argc, argv );

  set_FILTERSTRING(FILTERSTRING); // Mar 13 2014 - RK

  init_colorfuns();
  init_stretchfuns();
  init_shFudgefuns();

  read_input();

  set_SIMSED_PATH();

  if ( strcmp(INPUTS.SIMSED_VERSION_INPUT,"NULL") != 0 ) {
    read_SIMSED_INFO(SIMSED_PATHMODEL_INPUT);
  }
  else  { 
    init_SALT2_fudge(INPUTS.SALT2_VERSION_INPUT); 
  }


  if ( INDEX_COLORFUN  > 0 ) { check_for_S2x1();  }

  check_for_CORRECT(); // May 17, 2012

  init_REMOVE(); // make sure REMOVEd variables exist.

  mkdir_SIMSED_VERSION_OUTPUT();

  // write color fudge vs. lambda to text file
  if ( INDEX_COLORFUN  > 0 ) {  wr_colorFudge();   }

  init_SMEARMODEL();

  open_SEDINFO_file();


  printf("\n");
  NSED_OUTPUT = NCALL_SMOOTH = 0 ;

  // begin master loop over input SEDs
  for ( ised = 1 ; ised <= SEDMODEL.NSURFACE ; ised++ ) {

    if ( SKIPSED[ised] > 0 ) { continue ; }

    sedRead(ised);            // read SED from  file
    sedSmooth_driver(ised);   // smooth SED

    if ( INDEX_COLORFUN  > 0 ) { 
      sedFudge_color(ised);  // color fudge
    }

    if ( INPUTS.INDEX_SMEARMODEL > 0 ) {      
      int NBIN = SMEARPAR_GRID.NBIN ;
      if ( NBIN == 0 )
	{  sedFudge_smear(ised,-9); }
      else {
	printf("  %s grid-index: ", SMEARPAR_GRID.NAME );
	for(ismear=0; ismear < NBIN; ismear++ )  { 
	  printf("%d ", ismear);  fflush(stdout);
	  sedFudge_smear(ised,ismear); 	  
	}
	printf("\n"); fflush(stdout);
      }

    } // SMEARMODEL

  } // ised

  fclose(fp_sedinfo);

  if ( INPUTS.INDEX_SMEARMODEL > 0 ) {
    fclose(fp_sedinfo2);
    make_param_ntuple();
  }

  printf("\n");
  printf(" Finished writing %d new SEDs to \n", NSED_OUTPUT );
  printf("   %s\n", SIMSED_PATHMODEL_OUTPUT );
  printf(" PROGRAM ENDING GRACEFULLY \n");

} // end of main


// ******************************************
void parse_args(int argc, char **argv) {
  char fnam[] = "parse_args" ;
  // ---------- BEGIN --------
  sprintf( INPUTS.inputFile, "%s", argv[1] );
} // end of parse_args




// **************************************
void init_colorfuns(void) {

  // Initialize allowed color function names,
  // number of parameters, and parameter names.
  // Must be called BEFORE read_input().

  char fnam[] = "init_colorfuns" ;

  int INDEX, ipar, NPAR ;

  // ------------ BEGIN --------------

  // init COLORFUN structure
  for ( INDEX=0; INDEX <= NCHOICE_COLORFUN; INDEX++ ) {
    COLORFUN_LIST[INDEX].NPAR = 0 ;
    for ( ipar=0; ipar <= MXPARAM; ipar++ ) {
      sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"UNDEFINED");
    }
  }

  // =========

  INDEX = INDEX_SALT2colorlaw0 ;
  sprintf(COLORFUN_LIST[INDEX].NAME,"SALT2colorlaw0");
  COLORFUN_LIST[INDEX].NPAR = 5 ;
  ipar = 0 ;
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"LAM_B");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"LAM_V");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"COR0");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"COR1");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"BETA");


  INDEX = INDEX_SALT2colorlaw1 ;
  sprintf(COLORFUN_LIST[INDEX].NAME,"SALT2colorlaw1");
  COLORFUN_LIST[INDEX].NPAR = 10 ;
  ipar = 0;
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"LAM_B");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"LAM_V");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"MINLAM_INTERP");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"MAXLAM_INTERP");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"NPARAMS");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"PAR0");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"PAR1");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"PAR2");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"PAR3");
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"BETA");


  INDEX = INDEX_CCM89 ;
  sprintf(COLORFUN_LIST[INDEX].NAME,"CCM89");
  COLORFUN_LIST[INDEX].NPAR = 1 ;
  ipar = 0;
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"RV");


  INDEX = INDEX_CCM94 ;
  sprintf(COLORFUN_LIST[INDEX].NAME,"CCM94");
  COLORFUN_LIST[INDEX].NPAR = 1 ;
  ipar = 0;
  ipar++ ; sprintf(COLORFUN_LIST[INDEX].PARNAMES[ipar],"RV");

  // check that NPAR is at least 1, but less than MXPARAM


  for ( INDEX=1; INDEX <= NCHOICE_COLORFUN; INDEX++ ) {
    NPAR = COLORFUN_LIST[INDEX].NPAR ;
    if ( NPAR < 1 || NPAR >= MXPARAM ) {
      sprintf(c1err,"Invalid NPAR=%d for INDEX_COLORFUN=%d",
	      NPAR, INDEX );
      sprintf(c2err,"Check COLORFUN option %s", 
	      COLORFUN_LIST[INDEX].NAME );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    } 
  }   // INDEX



} // end of init_colorfuns

// **************************************
void init_stretchfuns(void) {

  // Initialize allowed stretch function names,
  // number of parameters, and parameter names.
  // Must be called BEFORE read_input().

  char fnam[] = "init_stretchfuns" ;

  int INDEX, ipar, NPAR ;

  // ------------ BEGIN --------------


  // init stretchfun structure
  for ( INDEX=0; INDEX <= NCHOICE_STRETCHFUN; INDEX++ ) {
    STRETCHFUN_LIST[INDEX].NPAR = 0 ;
    for ( ipar=0; ipar <= MXPARAM; ipar++ ) {
      sprintf(STRETCHFUN_LIST[INDEX].PARNAMES[ipar],"UNDEFINED");
    }
  }

  // =========

  INDEX = INDEX_STRETCHFUN_NONE ;
  sprintf(STRETCHFUN_LIST[INDEX].NAME,"NONE");
  STRETCHFUN_LIST[INDEX].NPAR = 1 ; //these are hardcoded!
  ipar = 0;
  ipar++ ; sprintf(STRETCHFUN_LIST[INDEX].PARNAMES[ipar],"HARDCODED");

  INDEX = INDEX_STRETCHFUN_DEFAULT ;
  sprintf(STRETCHFUN_LIST[INDEX].NAME,"DEFAULT");
  STRETCHFUN_LIST[INDEX].NPAR = 1 ; //these are hardcoded!
  ipar = 0;
  ipar++ ; sprintf(STRETCHFUN_LIST[INDEX].PARNAMES[ipar],"HARDCODED");


  INDEX = INDEX_STRETCHFUN_SAKO08 ;
  sprintf(STRETCHFUN_LIST[INDEX].NAME,"SAKO08");
  STRETCHFUN_LIST[INDEX].NPAR = 4 ;
  ipar = 0;
  ipar++ ; sprintf(STRETCHFUN_LIST[INDEX].PARNAMES[ipar],"PAR0");
  ipar++ ; sprintf(STRETCHFUN_LIST[INDEX].PARNAMES[ipar],"PAR1");
  ipar++ ; sprintf(STRETCHFUN_LIST[INDEX].PARNAMES[ipar],"PAR2");
  ipar++ ; sprintf(STRETCHFUN_LIST[INDEX].PARNAMES[ipar],"PAR3");

  // check that NPAR is at least 1, but less than MXPARAM

  for ( INDEX=1; INDEX <= NCHOICE_STRETCHFUN; INDEX++ ) {
    NPAR = STRETCHFUN_LIST[INDEX].NPAR ;
    if ( NPAR < 1 || NPAR >= MXPARAM ) {
      sprintf(c1err,"Invalid NPAR=%d for INDEX_STRETCHFUN=%d",
	      NPAR, INDEX );
      sprintf(c2err,"Check STRETCHFUN option %s", 
	      STRETCHFUN_LIST[INDEX].NAME );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    } 
  }   // INDEX



} // end of init_stretchfuns

// **************************************
void init_shFudgefuns(void) {

  // Initialize allowed shFudge function names,
  // number of parameters, and parameter names.
  // Must be called BEFORE read_input().

  char fnam[] = "init_shFudgefuns" ;

  int INDEX, ipar, NPAR ;

  // ------------ BEGIN --------------

  // init SHFUDGE structure
  for ( INDEX=0; INDEX <= NCHOICE_SHFUDGEFUN; INDEX++ ) {
    SHFUDGEFUN_LIST[INDEX].NPAR = 0 ;
    for ( ipar=0; ipar <= MXPARAM; ipar++ ) {
      sprintf(SHFUDGEFUN_LIST[INDEX].PARNAMES[ipar],"UNDEFINED");
    }
  }

  // =========

  INDEX = INDEX_SHFUDGE_DEFAULT ;
  sprintf(SHFUDGEFUN_LIST[INDEX].NAME,"DEFAULT");
  SHFUDGEFUN_LIST[INDEX].NPAR = 1 ;
  ipar = 0 ;
  ipar++ ; sprintf(SHFUDGEFUN_LIST[INDEX].PARNAMES[ipar],"COR0");

  INDEX = INDEX_SHFUDGE_SAKO08 ;
  sprintf(SHFUDGEFUN_LIST[INDEX].NAME,"SAKO08");
  SHFUDGEFUN_LIST[INDEX].NPAR = 4 ;
  ipar = 0;
  ipar++ ; sprintf(SHFUDGEFUN_LIST[INDEX].PARNAMES[ipar],"UNITY_VAL");
  ipar++ ; sprintf(SHFUDGEFUN_LIST[INDEX].PARNAMES[ipar],"a_PAR0");
  ipar++ ; sprintf(SHFUDGEFUN_LIST[INDEX].PARNAMES[ipar],"a_PAR1");
  ipar++ ; sprintf(SHFUDGEFUN_LIST[INDEX].PARNAMES[ipar],"b");

  // check that NPAR is at least 1, but less than MXPARAM

  for ( INDEX=1; INDEX <= NCHOICE_SHFUDGEFUN; INDEX++ ) {
    NPAR = SHFUDGEFUN_LIST[INDEX].NPAR ;
    if ( NPAR < 1 || NPAR >= MXPARAM ) {
      sprintf(c1err,"Invalid NPAR=%d for INDEX_SHFUDGEFUN=%d",
	      NPAR, INDEX );
      sprintf(c2err,"Check SHFUDGEFUN option %s", 
	      SHFUDGEFUN_LIST[INDEX].NAME );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    } 
  }   // INDEX



} // end of init_shFudgefuns



// **************************************
void read_input(void) {

  FILE *fp;

  char *ptrFile, *cptr ;
  char fnam[] = "read_input";

  char c_get[60], name[40] ;
  int N, i, ipar, NPAR, indx ;
  double tmp2[2] ;

  // -------- BEGIN --------


  // open user input file.

  ptrFile = INPUTS.inputFile;
  if ( ( fp = fopen(ptrFile, "rt") ) == NULL ) {
    sprintf(c1err, "Cannot find input file:");
    sprintf(c2err,"'%s'", ptrFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n Read input instructions from : %s \n", ptrFile);
  sprintf(c2err,"Check inputs file.");

  // set defaults  

  sprintf(INPUTS.SALT2_VERSION_INPUT,   "NULL");
  sprintf(INPUTS.SIMSED_VERSION_INPUT,  "NULL");
  sprintf(INPUTS.SIMSED_VERSION_OUTPUT, "NULL");

  FLAG_DEBUG_SIMSED = 0.0;

  sprintf(INPUTS.SHAPEPAR_NAME,"NULL");
  INPUTS.SHAPEPAR_RANGE[0]   = -9.0 ;
  INPUTS.SHAPEPAR_RANGE[1]   = -9.0 ;
  INPUTS.SHAPEPAR_BINSIZE    =  0.0 ;
  INPUTS.NBIN_SHAPEPAR       =  0   ;

  sprintf(INPUTS.COLORPAR_NAME,"NULL");
  INPUTS.COLORPAR_RANGE[0]   = -9.0 ;
  INPUTS.COLORPAR_RANGE[1]   = -9.0 ;
  INPUTS.COLORPAR_BINSIZE    =  0.0 ;
  INPUTS.NBIN_COLORPAR       =  0   ;

  INPUTS.SIMSED_TREST_RANGE[0] = -25.0 ;    // days
  INPUTS.SIMSED_TREST_RANGE[1] = +50.0 ;
  INPUTS.SIMSED_LAM_RANGE[0]   =  2000.0 ;  // Angstromgs
  INPUTS.SIMSED_LAM_RANGE[1]   = 10000.0 ;

  sprintf(INPUTS.COLORFUN_NAME,"NULL");
  INPUTS.NPAR_COLORFUN   = 0 ;
  INDEX_COLORFUN         = 0 ;

  sprintf(INPUTS.SHFUDGEFUN_NAME,"DEFAULT");
  INPUTS.NPAR_SHFUDGEFUN   = 1 ;
  INDEX_SHFUDGEFUN         = 1 ;

  sprintf(INPUTS.STRETCHFUN_NAME,"DEFAULT");
  INPUTS.NPAR_STRETCHFUN   = 1 ;
  INDEX_STRETCHFUN         = INDEX_STRETCHFUN_DEFAULT ;

  INPUTS.NPAR_SHAPEFUN   = 0 ; // redundant with SHFUDGE_FUN, 
  // val is set to 1 if user sets a SHAPEPAR_NAME


  for ( ipar=0; ipar < MXPARAM; ipar++ )  { 
    INPUTS.COLORFUN_PARVAL[ipar]   = NULLVAL ; 
    INPUTS.STRETCHFUN_PARVAL[ipar] = NULLVAL ;
    INPUTS.SHFUDGEFUN_PARVAL[ipar] = NULLVAL ; 
    INPUTS.SHAPEFUN_PARVAL[ipar]   = NULLVAL ; //redundant SHFUDGEFUN_PARVAL!
  }

  for ( i=0 ; i < MXLAMPOLY; i++ ) {
    INPUTS.LAMPOLY_TPEAK[i] = 0.0;
    INPUTS.LAMPOLY_FLUX[i]  = 0.0;
  }
  INPUTS.LAMPOLY_FLUX[0] = 1.0 ;  // default is no fudge

  sprintf(INPUTS.CORRECT_PEAKMB,"NULL");
  INPUTS.IPAR_PEAKMB = 0;
  INPUTS.IPAR_DM15   = 0;
  INPUTS.DUMP_LAMBDA = 0.0;

  INPUTS.NSEDPAR_REMOVE = 0;

  INPUTS.OPT_SMOOTH = OPT_NOSMOOTH ;
  INPUTS.NDUMP_SMOOTH = 0;
  INPUTS.LAMSCALE_SMOOTH = 50.0; // 50 A default

  INPUTS.INDEX_SMEARMODEL = -9 ;

  sprintf(INPUTS.SIMSED_MODELPATH,"NULL");

  // read input file

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"SALT2_VERSION_INPUT:") == 0 )  { 
      readchar(fp, INPUTS.SALT2_VERSION_INPUT ); 
      INDEX_STRETCHFUN = INDEX_STRETCHFUN_NONE ;
      INDEX_SHFUDGEFUN = INDEX_SHFUDGE_DEFAULT ;
      sprintf(INPUTS.STRETCHFUN_NAME,"NONE");
      sprintf(INPUTS.SHFUDGEFUN_NAME,"NONE");
    }

    if ( strcmp(c_get,"SIMSED_VERSION_INPUT:") == 0 )
      { readchar(fp, INPUTS.SIMSED_VERSION_INPUT ); }

    if ( strcmp(c_get,"SIMSED_VERSION_OUTPUT:") == 0 )
      { readchar(fp, INPUTS.SIMSED_VERSION_OUTPUT ); }

    if ( strcmp(c_get,"LAMPOLY_FLUX:") == 0 )
      { readdouble(fp, 4, INPUTS.LAMPOLY_FLUX ); }

    if ( strcmp(c_get,"DEBUG_SIMSED:") == 0 )
      { readdouble(fp, 1, &FLAG_DEBUG_SIMSED); }

    // -----

    if ( strcmp(c_get,"SHAPEPAR_NAME:") == 0 ) {
      readchar(fp, INPUTS.SHAPEPAR_NAME );
      INPUTS.NPAR_SHAPEFUN = 1 ;
    }

    if ( strcmp(c_get,"SHAPEPAR_RANGE:") == 0 )
      { readdouble(fp, 2, INPUTS.SHAPEPAR_RANGE ); }

    if ( strcmp(c_get,"SHAPEPAR_BINSIZE:") == 0 )
      { readdouble(fp, 1, &INPUTS.SHAPEPAR_BINSIZE ); }

    if ( strcmp(c_get,"LAMPOLY_TPEAK:") == 0 )
      { readdouble(fp, 3, INPUTS.LAMPOLY_TPEAK ); }

    // -------

    if ( strcmp(c_get,"COLORPAR_NAME:") == 0 )
      { readchar(fp, INPUTS.COLORPAR_NAME ); }

    if ( strcmp(c_get,"COLORPAR_RANGE:") == 0 )
      { readdouble(fp, 2, INPUTS.COLORPAR_RANGE ); }

    if ( strcmp(c_get,"COLORPAR_BINSIZE:") == 0 )
      { readdouble(fp, 1, &INPUTS.COLORPAR_BINSIZE ); }

    // --

    if ( strcmp(c_get,"SHFUDGEFUN_NAME:") == 0 ) {
      readchar(fp, name );
      sprintf(INPUTS.SHFUDGEFUN_NAME, "%s", name );   
      
      for ( i=1; i <= NCHOICE_SHFUDGEFUN; i++ ) {
	if ( strcmp(name,SHFUDGEFUN_LIST[i].NAME) == 0 ) {
	  INPUTS.NPAR_SHFUDGEFUN = SHFUDGEFUN_LIST[i].NPAR ;      
	  INDEX_SHFUDGEFUN       = i ;
	  }
      } // NCHOICE_SHFUDGEFUN
      
      if ( INPUTS.NPAR_SHFUDGEFUN == 0 ) {
	print_shfudgefun_options();
	sprintf(c1err,"Invalid SHFUDGEFUN_NAME: '%s' ", name );
	errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
      }
    } // SHFUDGEFUN_NAME

    if ( strcmp(c_get,"SHAPEFUN_PAR:") == 0 ) {
      readdouble(fp, INPUTS.NPAR_SHAPEFUN, &INPUTS.SHAPEFUN_PARVAL[1] );
      INPUTS.SHFUDGEFUN_PARVAL[1] = INPUTS.SHAPEFUN_PARVAL[1];
    }

    if ( strcmp(c_get,"SHFUDGEFUN_PAR:") == 0 ) {
      if ( INPUTS.NPAR_SHFUDGEFUN == 0 ) {
	sprintf(c1err,"Cannot  read SHFUDGEFUN_PAR because" );
	sprintf(c2err,"SHFUDGEFUN_NAME has  not been specified.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
      readdouble(fp, INPUTS.NPAR_SHFUDGEFUN, &INPUTS.SHFUDGEFUN_PARVAL[1] );
    }

    // --

    if ( strcmp(c_get,"STRETCHFUN_NAME:") == 0 ) {
      readchar(fp, name );
      sprintf(INPUTS.STRETCHFUN_NAME, "%s", name );   
      
      for ( i=1; i <= NCHOICE_STRETCHFUN; i++ ) {
	if ( strcmp(name,STRETCHFUN_LIST[i].NAME) == 0 ) {
	  INPUTS.NPAR_STRETCHFUN = STRETCHFUN_LIST[i].NPAR ;      
	  INDEX_STRETCHFUN       = i ;
	  }
      } // NCHOICE_STRETCHFUN
      
      if ( INPUTS.NPAR_STRETCHFUN == 0 ) {
	print_stretchfun_options();
	sprintf(c1err,"Invalid STRETCHFUN_NAME: '%s' ", name );
	errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
      }

    } // STRETCHFUN_NAME

    if ( strcmp(c_get,"STRETCHFUN_PAR:") == 0 ) {
      if ( INPUTS.NPAR_STRETCHFUN == 0 ) {
	sprintf(c1err,"Cannot  read STRETCHFUN_PAR because" );
	sprintf(c2err,"STRETCHFUN_NAME has  not been specified.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      } 
      if ( INDEX_STRETCHFUN == INDEX_STRETCHFUN_DEFAULT ) {
	sprintf(c1err, "No input STRETCHFUN_PAR allowed for %s", 
		INPUTS.STRETCHFUN_NAME);
	sprintf(c2err, "All parameters are hardwired.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      readdouble(fp, INPUTS.NPAR_STRETCHFUN, &INPUTS.STRETCHFUN_PARVAL[1] );
    } // STRETCHFUN_PAR
 
    // -------

    if ( strcmp(c_get,"COLORFUN_NAME:") == 0 ) {
      readchar(fp, name );
      sprintf(INPUTS.COLORFUN_NAME, "%s", name );

      for ( i=1; i <= NCHOICE_COLORFUN; i++ ) {
	if ( strcmp(name,COLORFUN_LIST[i].NAME) == 0 ) {
	  INPUTS.NPAR_COLORFUN = COLORFUN_LIST[i].NPAR ;      
	  INDEX_COLORFUN       = i ;
	  }
      } // NCHOICE_COLORFUN

      if ( INPUTS.NPAR_COLORFUN == 0 ) {
	print_colorfun_options();
	sprintf(c1err,"Invalid COLORFUN_NAME: '%s' ", name );
	errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
      }
    } // COLORFUN_NAME


    if ( strcmp(c_get,"COLORFUN_PAR:") == 0 ) {
      if ( INPUTS.NPAR_COLORFUN == 0 ) {
	sprintf(c1err,"Cannot  read COLORFUN_PAR because" );
	sprintf(c2err,"COLORFUN_NAME has  not been specified.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
      readdouble(fp, INPUTS.NPAR_COLORFUN, &INPUTS.COLORFUN_PARVAL[1] );

    }

    // -----

    if ( strcmp(c_get,"LAMBDA_RANGE:") == 0 ) {
      readdouble(fp, 2, tmp2 );

      if ( tmp2[1] < tmp2[0] || tmp2[0] <= 0 || tmp2[1] > 1E6 ) {
	sprintf(c1err,"Invalid LAMBDA_RANGE: %6.0f  %6.0f ", 
		tmp2[0],  tmp2[1]);
	sprintf(c2err,"Check input file '%s'", INPUTS.inputFile);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      INPUTS.SIMSED_LAM_RANGE[0] = tmp2[0] ;
      INPUTS.SIMSED_LAM_RANGE[1] = tmp2[1] ;
    }


    if ( strcmp(c_get,"EPOCH_RANGE:") == 0 ) {
      readdouble(fp, 2, tmp2 );
      if ( tmp2[1] < tmp2[0] || tmp2[0] > 50 || tmp2[1] < -20 ) {
	sprintf(c1err,"Invalid EPOCH_RANGE: %6.1f  %6.1f ", 
		tmp2[0],  tmp2[1]);
	sprintf(c2err,"Check input file '%s'", INPUTS.inputFile);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
      INPUTS.SIMSED_TREST_RANGE[0] = tmp2[0] ;
      INPUTS.SIMSED_TREST_RANGE[1] = tmp2[1] ;
    }

    if ( strcmp(c_get,"DUMP_LAMBDA:") == 0 )
      readdouble(fp, 1, &INPUTS.DUMP_LAMBDA );

    if ( strcmp(c_get,"CORRECT_PEAKMB:") == 0 ) {
      readchar(fp, INPUTS.CORRECT_PEAKMB );
      INPUTS.IPAR_PEAKMB = 99; // temp flag
    }
    if ( strcmp(c_get,"CORRECT_DM15:") == 0 ) {
      readchar(fp, INPUTS.CORRECT_DM15 );
      INPUTS.IPAR_DM15 = 99; // temp flag
    }


    //check for SED bins to remove
    if ( strcmp(c_get,"SED_REMOVE:") == 0 ) {
      readchar(fp, name );
      readint(fp, 1, &N);
      for ( i=1; i <= N; i++ ) {
	INPUTS.NSEDPAR_REMOVE++; ipar = INPUTS.NSEDPAR_REMOVE; 
	readdouble(fp, 1, &INPUTS.SEDPARVAL_REMOVE[ipar] );
	sprintf(INPUTS.SEDPARNAME_REMOVE[ipar], "%s", name );
      }
    }

    // check smooth options

    if ( strcmp(c_get,"OPT_SMOOTH:") == 0 )
      {  readint(fp, 1, &INPUTS.OPT_SMOOTH); }

    if ( strcmp(c_get,"LAMSCALE_SMOOTH:") == 0 )
      { readdouble(fp, 1, &INPUTS.LAMSCALE_SMOOTH); }


    if ( strcmp(c_get,"DUMP_SMOOTH:") == 0 ) {
      INPUTS.NDUMP_SMOOTH++ ;
      N = INPUTS.NDUMP_SMOOTH ;
      readchar(fp, INPUTS.DUMP_SMOOTH_SEDNAME[N] );
      readdouble(fp, 1, &INPUTS.DUMP_SMOOTH_TREST[N] );
    }


    // check for SMEARMODEL options
    if ( strcmp(c_get,"SMEAROPT:") == 0 ) {

      readchar(fp, name ) ;

      // fetch index to SMEARMODEL struct
      indx = init_SMEAROPT(name); 
      INPUTS.INDEX_SMEARMODEL = indx ;

      NPAR = SMEARMODEL_DEF[indx].NPAR ;

      // read each parameters as a string; then covert to double

      for ( ipar=1; ipar <= NPAR; ipar++ ) {
	cptr = INPUTS.SMEARMODEL_PARSTRING[ipar] ;
	readchar(fp, cptr );
	INPUTS.SMEARMODEL_PARVAL[ipar] = atof(cptr) ;
      }

    }


    if ( strcmp(c_get,"SIMSED_MODELPATH:") == 0 ) {
      readchar(fp, INPUTS.SIMSED_MODELPATH );
    }

    if ( strcmp(c_get,"END:") == 0 ) goto CLOSEINPUT ;

  }  // fscanf

  // ============

 CLOSEINPUT:
  fclose(fp);

  init_colorpar();
  init_shapepar();

  summarize_input();

} // end of read_input


// ====================================
void summarize_input(void) {
  
  // summarize a few things
  char *inVER, *ptrName, *ptrTxt ;
  int i; 
  int LSIMSED = ( strcmp(INPUTS.SIMSED_VERSION_INPUT,"NULL") != 0 ) ;
  double lam, parval, stmp ;

  printf("\n");

  if ( LSIMSED ) 
    { inVER = INPUTS.SIMSED_VERSION_INPUT ; }
  else
    { inVER = INPUTS.SALT2_VERSION_INPUT ; }

  printf(" Input  version: %s \n\n", inVER);
  printf(" Output version: %s \n", INPUTS.SIMSED_VERSION_OUTPUT );  

  if ( INPUTS.NBIN_COLORPAR > 1 ) 
    printf(" Output version includes '%s' with %d bins from %6.3f to %6.3f \n"
	   ,INPUTS.COLORPAR_NAME
	   ,INPUTS.NBIN_COLORPAR
	   ,INPUTS.COLORPAR_RANGE[0]
	   ,INPUTS.COLORPAR_RANGE[1]
	   );

  if ( INPUTS.NBIN_SHAPEPAR > 1 ) {
    printf(" Output version includes '%s' with %d bins from %6.3f to %6.3f \n"
	   ,INPUTS.SHAPEPAR_NAME
	   ,INPUTS.NBIN_SHAPEPAR
	   ,INPUTS.SHAPEPAR_RANGE[0]
	   ,INPUTS.SHAPEPAR_RANGE[1]
	   );

    for  ( i=1; i <= INPUTS.NBIN_SHAPEPAR; i++ ) {
      stmp = INPUTS.SHAPEPAR_VALUE[i];
      printf("\t stretch(%s = %6.3f) = %6.3f \n", 
	     INPUTS.SHAPEPAR_NAME, stmp, lcstretch(stmp) ); 
    }


    printf("\n STRETCHFUN NAME: %s \n", INPUTS.STRETCHFUN_NAME );
    if ( INDEX_STRETCHFUN > INDEX_STRETCHFUN_DEFAULT ) {
      for ( i=1; i <= INPUTS.NPAR_STRETCHFUN; i++ ) {
	parval  = INPUTS.STRETCHFUN_PARVAL[i];
	ptrName = STRETCHFUN_LIST[INDEX_STRETCHFUN].PARNAMES[i] ;
	printf(" STRETCHFUN PARAM %s = %6.3f \n", ptrName, parval);
      }
    }

    printf("\n SHFUDGEFUN NAME: %s \n", INPUTS.SHFUDGEFUN_NAME );
    if ( INDEX_SHFUDGEFUN > 1 ) {
      for ( i=1; i <= INPUTS.NPAR_SHFUDGEFUN; i++ ) {
	parval  = INPUTS.SHFUDGEFUN_PARVAL[i];
	ptrName = SHFUDGEFUN_LIST[INDEX_SHFUDGEFUN].PARNAMES[i] ;
	printf(" SHFUDGEFUN PARAM %s = %6.3e \n", ptrName, parval);
      }
    }

    printf("\n  Tpeak = %4.2f + %9.3le*(LAMBDA - %7.1f) \n"
	   ,INPUTS.LAMPOLY_TPEAK[0]
	   ,INPUTS.LAMPOLY_TPEAK[1]
	   ,INPUTS.LAMPOLY_TPEAK[2]
	   );
    lam = 3000.0;
    printf("  Tpeak(U2:%6.0f) = %5.2f \n", lam, TPEAK_LAMPOLY(lam) ) ;

    lam = 3560.0;
    printf("  Tpeak(U :%6.0f) = %5.2f \n", lam, TPEAK_LAMPOLY(lam) ) ;

    lam = 4380.0;
    printf("  Tpeak(B :%6.0f) = %5.2f \n", lam, TPEAK_LAMPOLY(lam) ) ;

    lam = 5490.0;
    printf("  Tpeak(V :%6.0f) = %5.2f \n", lam, TPEAK_LAMPOLY(lam) ) ;

    lam = 6520.0;
    printf("  Tpeak(R :%6.0f) = %5.2f \n", lam, TPEAK_LAMPOLY(lam) ) ;

    lam = 8030.0;
    printf("  Tpeak(I :%6.0f) = %5.2f \n", lam, TPEAK_LAMPOLY(lam) ) ;

    printf("\n");
  }

  printf(" Read/write SIMSEDs with  %6.1f < EPOCH  < %6.1f days \n",
	 INPUTS.SIMSED_TREST_RANGE[0], INPUTS.SIMSED_TREST_RANGE[1] );
  printf(" Read/write SIMSEDs with  %6.0f < LAMBDA < %6.0f days \n",
	 INPUTS.SIMSED_LAM_RANGE[0], INPUTS.SIMSED_LAM_RANGE[1] );	 

  printf("\n");

  printf(" COLORFUN NAME: %s \n", INPUTS.COLORFUN_NAME );
  for ( i=1; i <= INPUTS.NPAR_COLORFUN; i++ ) {
    parval  = INPUTS.COLORFUN_PARVAL[i];
    ptrName = COLORFUN_LIST[INDEX_COLORFUN].PARNAMES[i] ;
    printf(" COLORFUN PARAM %s = %6.3f \n", ptrName, parval);
  }

  // check smoothing options
  ptrTxt = INPUTS.COMMENT_SMOOTH ;

  sprintf(ptrTxt, "OPT_SMOOTH=%d: ", INPUTS.OPT_SMOOTH );
  
  if ( INPUTS.OPT_SMOOTH == OPT_NOSMOOTH )
    { sprintf(ptrTxt, "%s No spectral smoothing.", ptrTxt ) ; }
  else if ( INPUTS.OPT_SMOOTH == OPT_SMOOTH_TOPHAT )
    { sprintf(ptrTxt, "%s Rebin with TOP-HAT filter", ptrTxt ); }
  else if ( INPUTS.OPT_SMOOTH == OPT_SMOOTH_GAUSS )
    { sprintf(ptrTxt, "%s Rebin with Gaussian filter", ptrTxt); }
  else if ( INPUTS.OPT_SMOOTH == OPT_SMOOTH_NBRAVG )
    { sprintf(ptrTxt, "%s Average nearest neighbor.", ptrTxt); }
 

  if ( INPUTS.OPT_SMOOTH > 0 ) {
    sprintf(ptrTxt,"%s (%d Angstrom smooth scale)", 
	    ptrTxt, (int)INPUTS.LAMSCALE_SMOOTH);
  }
  printf(" %s\n", ptrTxt );


  // LAMPOLY fudge
  printf("\n Lambda-Polynomial Fudge: \n");
  for ( lam=2000.0 ; lam <= 8000; lam+= 1000.0 ) {
    printf("\t LAMPOLY_FLUX  Fudge(lam=%6.0f) = %6.3f \n", 
	  lam, FLUXFUDGE_LAMPOLY(lam) ) ;

  }

  for ( i=1; i <= INPUTS.NDUMP_SMOOTH; i++ ) {
    printf("   DUMP_SMOOTH: %s at Trest=%6.2f \n"
	   ,INPUTS.DUMP_SMOOTH_SEDNAME[i]
	   ,INPUTS.DUMP_SMOOTH_TREST[i] );
  }


  // -------------
  printf("\n");
  fflush(stdout);


} // end of read_input


// *************************************
void init_colorpar(void) {

  double dif, bin, XN, cmin, ctmp ;
  int i, N;
  char fnam[] = "init_colorpar" ;

  // calculate a few quantities
  bin = INPUTS.COLORPAR_BINSIZE ;
  dif = INPUTS.COLORPAR_RANGE[1] - INPUTS.COLORPAR_RANGE[0] ;
  if ( bin > 0.0 ) 
    XN  = dif / bin + 1.0 + 1.0E-8;
  else
    XN = 1.0 ;

  N = (int)(XN);
  //  printf(" xxxx N=%d  XN=%f \n", N, XN );

  if ( fabs((double)N - XN) > .00001 ) {
    sprintf(c1err,"COLOR_BINSIZE=%f does not match COLOR_RANGE.", bin);
    sprintf(c2err,"N=%d  XN=%f ", N, XN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  INPUTS.NBIN_COLORPAR = N ;

  if ( N >= MXCBIN ) {
    sprintf(c1err,"NBIN_COLORPAR=%d exceeds bound of %d", N, MXCBIN );
    sprintf(c2err,"Either increase MXCBIN or reduce no. color bins");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // load colorpar values
  cmin = INPUTS.COLORPAR_RANGE[0];
  for  ( i=1; i <= N; i++ ) {
    ctmp = cmin + bin * (double)(i-1) ;
    INPUTS.COLORPAR_VALUE[i] = ctmp ;
  }


} // end of init_colorpar



// *************************************
void init_shapepar(void) {


  double dif, bin, XN, smin, stmp ;
  int i, N;
  char fnam[] = "init_shapepar" ;

  // calculate a few quantities
  bin = INPUTS.SHAPEPAR_BINSIZE ;
  dif = INPUTS.SHAPEPAR_RANGE[1] - INPUTS.SHAPEPAR_RANGE[0] ;
  if ( bin > 0.0 ) 
    XN  = dif / bin + 1.0 + 1.0E-8;
  else
    XN = 1.0 ;

  N = (int)(XN);
  //  printf(" xxxx N=%d  XN=%f \n", N, XN );

  if ( fabs((double)N - XN) > .00001 ) {
    sprintf(c1err,"SHAPE_BINSIZE=%f does not match SHAPE_RANGE.", bin);
    sprintf(c2err,"N=%d  XN=%f ", N, XN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  INPUTS.NBIN_SHAPEPAR = N ;

  if ( N >= MXCBIN ) {
    sprintf(c1err,"NBIN_SHAPEPAR=%d exceeds bound of %d", N, MXCBIN );
    sprintf(c2err,"Either increase MXCBIN or reduce no. shape bins");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // load shapepar values
  smin = INPUTS.SHAPEPAR_RANGE[0];
  for  ( i=1; i <= N; i++ ) {
    stmp = smin + bin * (double)(i-1) ;
    INPUTS.SHAPEPAR_VALUE[i] = stmp ;
  }

} // end of init_shapepar


// *************************************
void print_colorfun_options(void) {

  int i;

  printf("\n COLORFUN_NAME options: \n" );

  for ( i = 1; i <= NCHOICE_COLORFUN; i++ ) {
    printf("\t %s  (%d params) \n", 
	   COLORFUN_LIST[i].NAME, COLORFUN_LIST[i].NPAR );
  }

  fflush(stdout);
} 

// *************************************
void print_stretchfun_options(void) {

  int i;

  printf("\n STRETCHFUN_NAME options: \n" );

  for ( i = 1; i <= NCHOICE_STRETCHFUN; i++ ) {
    printf("\t %s  (%d params) \n", 
	   STRETCHFUN_LIST[i].NAME, STRETCHFUN_LIST[i].NPAR );
  }

  fflush(stdout);
} 

// *************************************
void print_shfudgefun_options(void) {

  int i;

  printf("\n SHFUDGE_NAME options: \n" );

  for ( i = 1; i <= NCHOICE_SHFUDGEFUN; i++ ) {
    printf("\t %s  (%d params) \n", 
	   SHFUDGEFUN_LIST[i].NAME, SHFUDGEFUN_LIST[i].NPAR );
  }

  fflush(stdout);
} 


// ****************************************************
void init_SALT2_fudge(char *version) {

  int  ISED, ised ;
  char tmpFile[200], sedcomment[60] ;
  char fnam[] = "init_SALT2_fudge" ;

  // ------------- BEGIN ---------------

  printf(" Read tempaltes from %s \n", version); fflush(stdout);

  // ==========================================
  // construct path to SALT2 surfaces
  
  if ( getenv(PRIVATE_MODELPATH_NAME) == NULL ) {
    sprintf( SALT2_MODELPATH, "%s/models/SALT2/%s", 
	     getenv("SNDATA_ROOT"), version );
  }
  else {
    sprintf( SALT2_MODELPATH, "%s/%s", 
	     getenv(PRIVATE_MODELPATH_NAME), version );
  }
  
  SEDMODEL.NSURFACE = 2 ;
  SEDMODEL.OPTMASK  = 0 ;

  // get template name for each sed surface

  sprintf(SIMSED_PATHMODEL_INPUT, "%s", SALT2_MODELPATH);

  // store each template so that the M1 surface can be
  // used to stretch-correct the M0 surface.

  for ( ISED = 1 ; ISED <= SEDMODEL.NSURFACE ; ISED++ ) {

    ised = ISED-1 ; // indices for SALT2
    sprintf(SEDMODEL.FILENAME[ISED], "salt2_template_%d.dat", ised);

    sprintf(tmpFile, "%s/salt2_template_%d.dat", SALT2_MODELPATH, ised );    
    sprintf(sedcomment,"SALT2-%d", ised);

    rd_sedFlux(tmpFile, sedcomment
	       ,INPUTS.SIMSED_TREST_RANGE, INPUTS.SIMSED_LAM_RANGE
	       ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL, 0
	       ,&TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY, &TEMP_SEDMODEL.DAYSTEP
	       ,&TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM, &TEMP_SEDMODEL.LAMSTEP
	       ,TEMP_SEDMODEL.FLUX,  TEMP_SEDMODEL.FLUXERR );

    // make sure that DAY and LAM binning is identical for each surface

    check_sedflux_bins(ised, "DAY", 
       TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY[0], TEMP_SEDMODEL.DAYSTEP);
    check_sedflux_bins(ised, "LAM", 
       TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM[0], TEMP_SEDMODEL.LAMSTEP);

    // transfer TEMP_SEDMODEL to permanent storage
    fill_SALT2_TABLE_SED(ised);

  } //  end loop over SED templates

  SEDMODEL.NSURFACE   = 1;
  SEDMODEL.FLUXSCALE  = 0.78 ;
  SEDMODEL.MAGERR_FIX = 0.10 ;

} // end of init_SALT2_fudge


// ****************************************************
int init_SMEAROPT(char *modelName) {

  // Created Oct 3, 2011 by R.Kessler
  // initialize hard-wired values for each model.
  // Then match input modelName and return appropriate index
  // to quickly access information about modelName.
  // 

  int NDEF, i, indx; 
  char *cptr, cfilt[20];
  char fnam[] = "nit_SMEAROPT" ; 

  // fill all possible quantities for requested model *name

  NDEF = 0;
  SMEARMODEL_ID.USE_genSmear = 0;
  SMEARPAR_GRID.NBIN = 0 ;   // default

  // ----------------
  NDEF++ ; 
  sprintf(SMEARMODEL_DEF[NDEF].NAME,    "%s", "NOSMEAR");
  sprintf(SMEARMODEL_DEF[NDEF].COMMENT, "%s", 
	  "No Smearing; re-write same SEDs");
  SMEARMODEL_DEF[NDEF].NPAR     = 0 ;
  SMEARMODEL_DEF[NDEF].NGAURAN  = 0 ;
  SMEARMODEL_DEF[NDEF].NFLATRAN = 0 ;
  SMEARMODEL_ID.NOSMEAR         = NDEF ;

  // ----------------
  NDEF++ ; 
  sprintf(SMEARMODEL_DEF[NDEF].NAME,    "%s", "COHERENT");
  sprintf(SMEARMODEL_DEF[NDEF].COMMENT, "%s", 
	  "mag-shift = Gran*P1");
  SMEARMODEL_DEF[NDEF].NPAR     = 1 ;
  SMEARMODEL_DEF[NDEF].NGAURAN  = 1 ;
  SMEARMODEL_DEF[NDEF].NFLATRAN = 0 ;
  SMEARMODEL_ID.COHERENT        = NDEF ;

  // ----------------
  NDEF++ ; 
  sprintf(SMEARMODEL_DEF[NDEF].NAME, "%s", "SIN(LAMBDA)");
  sprintf(SMEARMODEL_DEF[NDEF].COMMENT, "%s", 
	  "mag-shift = Gran*P1*sin[2PI*(lam-P2)/P3] ");
  SMEARMODEL_DEF[NDEF].NPAR = 3 ; // Amp, lambda(node),  Del-lambda(2PI)
  SMEARMODEL_DEF[NDEF].NGAURAN  = 1 ;
  SMEARMODEL_DEF[NDEF].NFLATRAN = 0 ;
  SMEARMODEL_ID.SINX            = NDEF ;

  // ----------------
  NDEF++ ; 
  sprintf(SMEARMODEL_DEF[NDEF].NAME, "%s", "SIMSED.KRW09");
  sprintf(SMEARMODEL_DEF[NDEF].COMMENT, "%s", 
	  "KRW09 cosangle : P1=sub-model, P2=ampl. P3=smoothLam" );
  SMEARMODEL_DEF[NDEF].NPAR     = 3 ; // sub-model  amplitude-scale
  SMEARMODEL_DEF[NDEF].NGAURAN  = 0 ;
  SMEARMODEL_DEF[NDEF].NFLATRAN = 2 ; // IMODEL and COSANGLE
  SMEARMODEL_ID.KRW09           = NDEF ;

  
  NDEF++ ; 
  sprintf(SMEARMODEL_DEF[NDEF].NAME, "%s", "G10");
  sprintf(SMEARMODEL_DEF[NDEF].COMMENT, "%s", 
	  "Guy10 : no free params" );
  SMEARMODEL_DEF[NDEF].NPAR     = 0  ; 
  SMEARMODEL_ID.G10           = NDEF ;


  NDEF++ ; 
  sprintf(SMEARMODEL_DEF[NDEF].NAME, "%s", "C11");
  sprintf(SMEARMODEL_DEF[NDEF].COMMENT, "%s", 
	  "Chotard 2011 thesis : P1 = 0,1,2 -> far UV option" );
  SMEARMODEL_DEF[NDEF].NPAR     = 1  ; // sub-model  amplitude-scale
  SMEARMODEL_ID.C11             = NDEF ;


  NDEF++ ; 
  sprintf(SMEARMODEL_DEF[NDEF].NAME, "%s", "VCR");
  sprintf(SMEARMODEL_DEF[NDEF].COMMENT, "%s", 
	  "Velocity-Color-Relation (MFK14): no free params" );
  SMEARMODEL_DEF[NDEF].NPAR   = 0  ; 
  SMEARMODEL_ID.VCR           = NDEF ;


  // ===============================
  // load global
  NDEF_SMEARMODEL = NDEF;

  // now match input modelName and also summary all available models
  // so there is always a reminder.

  indx = -9;
  for ( i = 1; i <= NDEF; i++ ) {

    cptr = SMEARMODEL_DEF[i].NAME ;
    printf(" Available SMEAROPT:  %-22s + %d parameters \n",
	   cptr, SMEARMODEL_DEF[i].NPAR );

    if ( strcmp(modelName,cptr) == 0 )  { indx = i ; }
  }
  
  // abort if model is not defined
  if ( indx < 0 ) {
    sprintf(c1err,"Unknown SMEARMODEL '%s'", modelName );
    sprintf(c2err,"See list of valid SMEARMODELs above." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // define filters for synthetic mags
  sprintf(cfilt, "UBVRI" );
  SMEARMODEL_DEF[indx].NFILTMAG = strlen(cfilt) ;
  sprintf(SMEARMODEL_DEF[indx].FILT, "%s", cfilt );

  // init random number generator
  int ISEED = 12345 ;
  srandom(ISEED);
 

  // return index to this SMEARMODEL
  return indx ;

} // end of init_SMEAROPT



// *************************************
void set_SIMSED_PATH(void) {


  if ( strcmp(INPUTS.SALT2_VERSION_INPUT,"NULL") != 0 ) {
    int Nc, Ns ;
    double DIF, BIN ;

    // for SALT2 -> SIMSED option, set output SIMSED dirName 
    // to combine number of color bins and number of x1 bins.
    DIF = INPUTS.SHAPEPAR_RANGE[1] - INPUTS.SHAPEPAR_RANGE[0] ;
    BIN = INPUTS.SHAPEPAR_BINSIZE ;
    Nc  = (int)((DIF+0.0001)/BIN) + 1;

    DIF = INPUTS.COLORPAR_RANGE[1] - INPUTS.COLORPAR_RANGE[0] ;
    BIN = INPUTS.COLORPAR_BINSIZE ;
    Ns  = (int)((DIF+0.0001)/BIN) + 1;

    sprintf(INPUTS.SIMSED_VERSION_OUTPUT,"SIMSED.%s_%2.2dc+%2.2dx1",
	    INPUTS.SALT2_VERSION_INPUT, Nc, Ns);	    
  }

  // ------------
  if ( strcmp(INPUTS.SIMSED_MODELPATH,"NULL") == 0 ) {
    // default path
    sprintf( SIMSED_PATHMODEL_INPUT, "%s/models/SIMSED/%s", 
	     PATH_SNDATA_ROOT, INPUTS.SIMSED_VERSION_INPUT );
    
    sprintf( SIMSED_PATHMODEL_OUTPUT, "%s/models/SIMSED/%s", 
	     PATH_SNDATA_ROOT, INPUTS.SIMSED_VERSION_OUTPUT );

  }
  else {
    // use-specified path
    sprintf( SIMSED_PATHMODEL_INPUT, "%s/%s", 
	     INPUTS.SIMSED_MODELPATH, INPUTS.SIMSED_VERSION_INPUT );

    sprintf( SIMSED_PATHMODEL_OUTPUT, "%s/%s", 
	     INPUTS.SIMSED_MODELPATH, INPUTS.SIMSED_VERSION_OUTPUT );
  }


} // end of set_SIMSED_PATH


// ***********************************
void  mkdir_SIMSED_VERSION_OUTPUT(void) {

  // use system call to create output SIMSED version

  char 
    cmd[400]
    ,rmpath[400]
    ,infoFile[200]
    ,newpath[200]
    ,fnam[] = "mkdir_SIMSED_VERSION_OUTPUT"
    ;

  int isys;
  FILE *fp;

  // ------------ BEGIN ----------

  sprintf(newpath,"%s", SIMSED_PATHMODEL_OUTPUT ) ;
  sprintf(infoFile,"%s/%s", newpath, INFO_SIMSED_FILENAME );

  // if SED.INFO file already exists, then replace it

  if ( ( fp = fopen(infoFile, "rt") ) != NULL ) {
    sprintf(c1err,"Replace SIMSED directory");
    sprintf(rmpath,"rm -rf %s ; ", newpath );
  }
  else {
    sprintf(c1err,"Create SIMSED directory");
    sprintf(rmpath,"" );
  }

  printf("\n %s : \n   %s \n", c1err, newpath );

  // construct command to create new SIMSED version

  sprintf(cmd,"%s mkdir -m g+wr %s", rmpath, newpath );

  printf("\n%s \n\n", cmd);
  isys = system(cmd);

  fflush(stdout);

} // end of  mkdir_SIMSED_VERSION_OUTPUT


// *************************************
void wr_colorFudge(void) {

  char tmpFile[200];
  char fnam[] = "wr_colorFudge" ;

  double c = 1.0;
  double lam, fudge ;

  FILE *fp;

  // ------------ BEGIN -------------

  sprintf(tmpFile,"%s/%s", SIMSED_PATHMODEL_OUTPUT, colorFudge_datFile );
  fp = fopen(tmpFile, "wt") ;

  for ( lam=2000.0 ; lam <= 9000.0; lam+=100.0 ) {
    fudge  = colorFudge(lam,c);
    fprintf(fp,"%8.1f  %7.5f \n", lam, fudge );
  }
  fclose(fp);

} // end of wr_colorFudge


// *************************************
void open_SEDINFO_file(void) {

  // Open new SED.INFO file and write header info
  // Note that extra color variable is added to VARNAMES list

  char
    infoFile[MXPATHLEN]
    ,infoFile2[MXPATHLEN]
    ,fnam[] = "open_SEDINFO_file"
    ;

  int ipar, indx, i, NRAN ;

  // ----- BEGIN -----

  sprintf(infoFile,"%s/%s", SIMSED_PATHMODEL_OUTPUT, INFO_SIMSED_FILENAME );
  fp_sedinfo = fopen(infoFile, "wt") ;

  printf(" Open %s file \n", INFO_SIMSED_FILENAME);

  fprintf(fp_sedinfo,
	  "# ------------------------------------------------------ \n");
  fprintf(fp_sedinfo,"# Created by program $SNANA_DIR/bin/SIMSED_fudge.exe \n");

  fprintf(fp_sedinfo,"# $SNANA_DIR = %s \n", getenv("SNANA_DIR") );
  fprintf(fp_sedinfo,"# $HOST      = %s \n", getenv("HOST") );

  fprintf(fp_sedinfo,"# \n");

  fprintf(fp_sedinfo,"# SEDs modified from version %s \n", 
	  INPUTS.SIMSED_VERSION_INPUT);

  fprintf(fp_sedinfo,"# \n");
  fprintf(fp_sedinfo,"# %s \n", INPUTS.COMMENT_SMOOTH );
  fprintf(fp_sedinfo,"# \n");

  if ( INDEX_COLORFUN > 0          ) { wrhead_SEDINFO_color();  }

  if ( INPUTS.INDEX_SMEARMODEL > 0 ) { wrhead_SEDINFO_smear(); }

  fprintf(fp_sedinfo,
	  "# ------------------------------------------------------ \n");
  fprintf(fp_sedinfo,"\n");

  // now re-write the header info, tacking on the extra color param

  if ( SEDMODEL.FLUXSCALE < 1000.0 ) 
    fprintf(fp_sedinfo,"FLUX_SCALE: %6.3f \n", SEDMODEL.FLUXSCALE );
  else
    fprintf(fp_sedinfo,"FLUX_SCALE: %9.3le \n", SEDMODEL.FLUXSCALE );

  fprintf(fp_sedinfo,"MAG_ERR:    %6.3f \n", SEDMODEL.MAGERR_FIX );
  fprintf(fp_sedinfo,"RESTLAMBDA_RANGE: %5.0f %6.0f \n", 
	  INPUTS.SIMSED_LAM_RANGE[0], INPUTS.SIMSED_LAM_RANGE[1] );
  fprintf(fp_sedinfo,"FLUX_ERRFLAG: %d \n", 
	  (SEDMODEL.OPTMASK & OPTMASK_FLUXERR_SEDMODEL)  );
  fprintf(fp_sedinfo,"\n");

  fprintf(fp_sedinfo,"NPAR:  %d \n", NPAR_SEDINFO_OUTPUT );
  fprintf(fp_sedinfo,"PARNAMES:  " );

  //open 2nd SED-info file in fitres-format to make ntuple
  if ( INPUTS.INDEX_SMEARMODEL > 0 ) {
    int NVAR ;
    NVAR = NPAR_SEDINFO_OUTPUT ;     NVAR++ ;  // add CID
    sprintf(infoFile2,"%s/%s", SIMSED_PATHMODEL_OUTPUT, INFO_FILENAME2 );
    fp_sedinfo2 = fopen(infoFile2, "wt") ;
    fprintf(fp_sedinfo2,"NVAR:  %d \n", NVAR );
    fprintf(fp_sedinfo2,"VARNAMES: CID   " );
  }

  for ( ipar = 1; ipar <= NPAR_SEDINFO_OUTPUT; ipar++ ) {
    fprintf(fp_sedinfo, "%s ", SEDMODEL.PARNAMES[ipar] );

    if ( INPUTS.INDEX_SMEARMODEL > 0 ) 
      {  fprintf(fp_sedinfo2,"%s ", SEDMODEL.PARNAMES[ipar] ); }

  } // ipar

  fprintf(fp_sedinfo,  "\n\n");

  if ( INPUTS.INDEX_SMEARMODEL > 0 ) {  fprintf(fp_sedinfo2, "\n\n"); }

  fflush(fp_sedinfo);

  //  fclose(fp_sedinfo);


} // end of open_SEDINFO_file



// *******************************************
void wrhead_SEDINFO_color(void) {

  char fnam[] = "wrhead_SEDINFO_color" ;

  int NPAR, ipar, i;
  double parval;
  char *cptr;

  // ------------ BEGIN -------------


  NPAR = SEDMODEL.NPAR;
  if ( INPUTS.NBIN_COLORPAR > 1 ) {
    NPAR++ ;
    sprintf(SEDMODEL.PARNAMES[NPAR],"%s", INPUTS.COLORPAR_NAME );
  }
  if ( INPUTS.NBIN_SHAPEPAR > 1 ) {
    NPAR++ ;
    sprintf(SEDMODEL.PARNAMES[NPAR],"%s", INPUTS.SHAPEPAR_NAME );
  }

  // Note that SEDMODEL.NPAR is for the input model;
  // NPAR_SEDINFO_OUTPUT is for the output model.

  NPAR_SEDINFO_OUTPUT = NPAR ;

  if ( INPUTS.NBIN_COLORPAR > 1 ) {
    fprintf(fp_sedinfo,
	    "# This version includes '%s' with %d bins from %6.3f to %6.3f \n"
	    ,INPUTS.COLORPAR_NAME
	    ,INPUTS.NBIN_COLORPAR
	    ,INPUTS.COLORPAR_RANGE[0]
	    ,INPUTS.COLORPAR_RANGE[1]
	    );
  }

  if ( INPUTS.NBIN_SHAPEPAR > 1 ) {
    fprintf(fp_sedinfo,
	    "# This version includes '%s' with %d bins from %6.3f to %6.3f \n"
	    ,INPUTS.SHAPEPAR_NAME
	    ,INPUTS.NBIN_SHAPEPAR
	    ,INPUTS.SHAPEPAR_RANGE[0]
	    ,INPUTS.SHAPEPAR_RANGE[1]
	    );
  }


  if ( INDEX_STRETCHFUN > INDEX_STRETCHFUN_DEFAULT ) {
    fprintf(fp_sedinfo, "# \n# This version uses stretch function %s \n", 
	    INPUTS.STRETCHFUN_NAME );
      for ( i=1; i <= INPUTS.NPAR_STRETCHFUN; i++ ) {
	parval  = INPUTS.STRETCHFUN_PARVAL[i];
	cptr = STRETCHFUN_LIST[INDEX_STRETCHFUN].PARNAMES[i] ;
	fprintf(fp_sedinfo,"#     STRETCHFUN PARAM %s = %6.3e \n", cptr, parval);
      }
  }


  if ( IPAR_S2x1 > 0 ) {
    cptr = SEDMODEL.PARNAMES[IPAR_S2x1] ;    

    fprintf(fp_sedinfo, "# This version is scaled by 10^[0.4 * %s * %4.2f] \n", 
	   cptr, S2ALPHA);

  } else if (INDEX_SHFUDGEFUN == 1) {
    cptr = INPUTS.SHAPEPAR_NAME;

    fprintf(fp_sedinfo, "# This version is scaled by 10^[0.4 * %s * %4.2f] \n", 
	    cptr, INPUTS.SHAPEFUN_PARVAL[1]);

  } else if (INDEX_SHFUDGEFUN == 2) {
    cptr = INPUTS.SHAPEPAR_NAME;

    fprintf(fp_sedinfo, "# \n# This version is scaled by 10^[-0.4 * (a * %s + b * %s * %s))\n",
	    "X", "X", "X");
    fprintf(fp_sedinfo, "#    HERE, X = %s - UNITY_VAL \n", cptr);
    fprintf(fp_sedinfo, "#    HERE, a(lambda) = c_0 - c_1 * lambda [Angstroms] \n");
      for ( i=1; i <= INPUTS.NPAR_SHFUDGEFUN; i++ ) {
	parval  = INPUTS.SHFUDGEFUN_PARVAL[i];
	cptr = SHFUDGEFUN_LIST[INDEX_SHFUDGEFUN].PARNAMES[i] ;
	fprintf(fp_sedinfo,"#     SHFUDGEFUN PARAM %s = %6.3e \n", cptr, parval);
      }

  }

  if ( INPUTS.IPAR_PEAKMB > 0 ) {
    fprintf(fp_sedinfo,"# %s  is corrected. \n", INPUTS.CORRECT_PEAKMB );
  }

  if ( INPUTS.IPAR_DM15 > 0 ) {
    fprintf(fp_sedinfo,"# %s  is corrected. \n", INPUTS.CORRECT_DM15 );
  }


  fprintf(fp_sedinfo,"# \n");
  fprintf(fp_sedinfo,"# Flux(fudged color law with %s=1)/Flux  vs. wavelength \n",
	  INPUTS.COLORPAR_NAME );
  fprintf(fp_sedinfo,"# is in %s \n",  colorFudge_datFile );


  // print the color-fun parameters
  fprintf(fp_sedinfo,"# %s Parameters: \n", 
	  COLORFUN_LIST[INDEX_COLORFUN].NAME );
  for ( ipar= 1; ipar <= INPUTS.NPAR_COLORFUN; ipar++ ) {
    fprintf(fp_sedinfo,"#  %15s:  %14.6le  \n" 
	    ,COLORFUN_LIST[INDEX_COLORFUN].PARNAMES[ipar]
	    ,INPUTS.COLORFUN_PARVAL[ipar] );
  }


  // check for MB abd DM15 parameters to correct
  ipar = INPUTS.IPAR_PEAKMB ; 
  if ( ipar == 99 ) {
    sprintf(c1err,"Could not find '%s' variable to correct", 
	    INPUTS.CORRECT_PEAKMB);
    sprintf(c2err,"peak B-band mag. Check CORRECT_PEAKMB key.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  if ( ipar > 0 ) {
    printf(" Will correct output '%s' (ipar=%d) \n", 
	   INPUTS.CORRECT_PEAKMB, ipar );
  }


  ipar = INPUTS.IPAR_DM15 ; 
  if ( ipar == 99 ) {
    sprintf(c1err,"Could not find '%s' variable to correct", 
	    INPUTS.CORRECT_DM15 );
    sprintf(c2err,"DM15. Check CORRECT_DM15 key.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  if ( ipar > 0 ) {
    printf(" Will correct output '%s' (ipar=%d) \n", 
	   INPUTS.CORRECT_DM15, ipar );
  }



} // end of wrhead_SEDINFO_color



// *******************************************
void wrhead_SEDINFO_smear(void) {

  int indx, ipar, i, NRAN, NPAR, NPAR_SMEAR;
  int NFILT, ifilt ;
  char fnam[] = "wrhead_SEDINFO_smear" ;
  char *cptr;

  // ------------ BEGIN -------------

  indx = INPUTS.INDEX_SMEARMODEL ;

  NPAR =  SEDMODEL.NPAR ;  // starting NPAR

  // tack on random number(s) to PARNAMES list
  NRAN = SMEARMODEL_DEF[indx].NGAURAN ;
  for ( i = 1; i <= NRAN; i++ ) {
    NPAR++ ;    SEDMODEL.NPAR = NPAR ;
    cptr  = SEDMODEL.PARNAMES[NPAR] ;
    sprintf(cptr, "GAURAN%d", i );
  }


  NRAN = SMEARMODEL_DEF[indx].NFLATRAN ;
  for ( i = 1; i <= NRAN; i++ ) {
    NPAR++ ;    SEDMODEL.NPAR = NPAR ;
    cptr  = SEDMODEL.PARNAMES[NPAR] ;
    sprintf(cptr, "FLATRAN%d", i );
  }

  // tack on synthetic mag-difs
  NFILT = SMEARMODEL_DEF[indx].NFILTMAG ;
  for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
    NPAR++ ;  SEDMODEL.NPAR = NPAR ;
    cptr  = SEDMODEL.PARNAMES[NPAR] ;
    sprintf(cptr,"MAGDIF_%c", SMEARMODEL_DEF[indx].FILT[ifilt] );
  }

  // Jun 29, 2012: 
  // tack on DM15 and MB for KRW09 models (if they are defined)
  if ( KRW09.N_COSANGLE > 0 ) {

    NPAR++ ;  SEDMODEL.NPAR = NPAR ;
    cptr  = SEDMODEL.PARNAMES[NPAR] ;
    sprintf(cptr,"COSANG" );

    NPAR++ ;  SEDMODEL.NPAR = NPAR ;
    cptr  = SEDMODEL.PARNAMES[NPAR] ;
    sprintf(cptr,"DM15" );
  
    NPAR++ ;  SEDMODEL.NPAR = NPAR ;
    cptr  = SEDMODEL.PARNAMES[NPAR] ;
    sprintf(cptr,"MB" );
  }


  // load global
  NPAR_SEDINFO_OUTPUT = NPAR;


  // -------------------------------------------------
  // tack on optional SMEARPAR 

  if ( SMEARPAR_GRID.NBIN > 0 ) {
    NPAR++ ;    NPAR_SEDINFO_OUTPUT = NPAR ;
    sprintf(SEDMODEL.PARNAMES[NPAR],"%s", SMEARPAR_GRID.NAME);
  }

  

  // ---------------------------------
  // write comments to SED.INFO header
  fprintf(fp_sedinfo,"# SMEAROPT: %s \n", 
	  SMEARMODEL_DEF[indx].NAME );
  fprintf(fp_sedinfo,"#\t %s \n", 
	  SMEARMODEL_DEF[indx].COMMENT );


  // write each SMEARMODEL parameter value
  NPAR_SMEAR = SMEARMODEL_DEF[indx].NPAR ;
  for ( ipar=1; ipar <= NPAR_SMEAR; ipar++ ) {
    fprintf(fp_sedinfo,"#\t P%d = %s \n",
	    ipar, INPUTS.SMEARMODEL_PARSTRING[ipar] );
  }
  

} // end of wrhead_SEDINFO_smear


// **************************************
void check_for_CORRECT(void) {

  // if CORRECT_XXX key is set then load INPUTS.IPAR_PEAKMB/DM15

  int ipar ;
  char fnam[] = "check_for_CORRECT" ;

  // ------------  BEGIN -------------

  for ( ipar = 1; ipar <= SEDMODEL.NPAR; ipar++ ) {

    // check for peakMB variable
    if ( strcmp(SEDMODEL.PARNAMES[ipar],INPUTS.CORRECT_PEAKMB) == 0 ) 
      { INPUTS.IPAR_PEAKMB = ipar ; }

    // check for DM15 variable
    if ( strcmp(SEDMODEL.PARNAMES[ipar],INPUTS.CORRECT_DM15) == 0 ) 
      { INPUTS.IPAR_DM15 = ipar ; }

  } // ipar

} // end of check_for_CORRECT

// *******************************************
void check_for_S2x1(void) {

  // Mar 11, 2011 R.Kessler
  // if any input parameter is x1 or S2x1, then store
  // IPAR_S2x1 = ipar and use this parameter later
  // to normalize the output SEDs.


  int  ipar;
  char *cptr;

  char fnam[] = "check_for_S2x1";

  // -------------- BEGIN -------------------

  IPAR_S2x1 = -9 ;

  for ( ipar = 1 ; ipar <= SEDMODEL.NPAR ; ipar++ ) {
    cptr = SEDMODEL.PARNAMES[ipar];
    if ( strcmp(cptr,"S2x1")    == 0 ) { IPAR_S2x1 = ipar ; }
    if ( strcmp(cptr,"x1")      == 0 ) { IPAR_S2x1 = ipar ; }
    if ( strcmp(cptr,"SALT2x1") == 0 ) { IPAR_S2x1 = ipar ; }
  }

  if ( IPAR_S2x1 > 0 ) {
    INPUTS.SHAPEFUN_PARVAL[1] = S2ALPHA ;
    cptr = SEDMODEL.PARNAMES[IPAR_S2x1];    
    printf("\t Found '%s' input parameter => \n", cptr );
    printf("\t will scale output SED flux by 10^[0.4 * %s * %4.2f] \n", 
	   cptr, S2ALPHA);
  }


} // end of  check_for_S2x1

// ********************************
void init_REMOVE(void) {

  // for each SED bin to remove, store IPAR index
  // to speed up the lookup for later.
  // Also screen dump each SED bin to remove,
  // and abort on any undefined parameter names.
  // Finally, set SKIPSED[ised] logicals.

  int i, ipar, ipar_remove, ised, Ntmp ;
  double parval_user, parval_sed, dif ;

  char *ptrname;
  char fnam[] = "init_REMOVE" ;

  // ---------- BEGIN --------

  if ( INPUTS.NSEDPAR_REMOVE <= 0 ) return ;

  printf("\n");

  for ( i = 1; i <= INPUTS.NSEDPAR_REMOVE; i++ ) {
    NSED_REMOVED[i] = 0; // init counter
    ptrname = INPUTS.SEDPARNAME_REMOVE[i];

    ipar_remove = 0;
    for ( ipar=1; ipar <= SEDMODEL.NPAR; ipar++ ) {
      if ( strcmp(ptrname,SEDMODEL.PARNAMES[ipar]) == 0 ) {
	ipar_remove = ipar ;
      }
    } // ipar

    if ( ipar_remove > 0 ) { 
      INPUTS.IPAR_REMOVE[i] = ipar_remove ;
      printf("\t remove SED bins with %s = %f  (ipar=%d) \n"
	     ,INPUTS.SEDPARNAME_REMOVE[i]
	     ,INPUTS.SEDPARVAL_REMOVE[i]
	     ,ipar_remove ) ;
    }
    else {
      sprintf(c1err,"Cannot remove non-existing SIMSED parameter '%s'", 
	      ptrname);
      sprintf(c2err,"Check SED_REMOVE keyword(s).");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  } // i loop over REMOVE


  // -----------------
  // Now set SKIPSED logicals

  // first init all SKIPSED to false
  for ( i=0; i < MXSEDMODEL; i++ )  { SKIPSED[i] = 0; }

  for ( ised=1; ised <= SEDMODEL.NSURFACE; ised++ ) {

    for ( i=1; i <= INPUTS.NSEDPAR_REMOVE; i++ ) {
      ipar_remove = INPUTS.IPAR_REMOVE[i] ;
      parval_sed  = SEDMODEL.PARVAL[ised][ipar_remove];
      parval_user = INPUTS.SEDPARVAL_REMOVE[i];
      dif         = fabs(parval_user-parval_sed);

      if ( ised == -2 ) 
	printf(" xxx ipar_remove=%d(%s) parval[sed,user]=%6.3f,%6.3f \n",
	       ipar_remove, INPUTS.SEDPARNAME_REMOVE[i],
	       parval_sed, parval_user );

      if ( dif < 1.0E-6 ) {
	SKIPSED[ised] = 1;    
	NSED_REMOVED[i]++ ;
      }

    } // i
  } // ised

  // finally check that each remove-bin actually did something
  for ( i=1; i <= INPUTS.NSEDPAR_REMOVE; i++ ) {
    if ( NSED_REMOVED[i] <= 0 ) {
      sprintf(c1err,"Did not find any SED_REMOVE bins for");
      sprintf(c2err,"%s = %f", 
	      INPUTS.SEDPARNAME_REMOVE[i], INPUTS.SEDPARVAL_REMOVE[i] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

} // end of init_REMOVE



// ********************************
void sedRead(int ised) {

  char
    fnam[] = "sedRead"
    ,sed_inFile_full[200]  
    ,sedcomment[100] 
    ;

  int N;

  // ------------ BEGIN ------------

  sprintf(sed_inFile_full, "%s/%s", 
	  SIMSED_PATHMODEL_INPUT, SEDMODEL.FILENAME[ised] );

  sprintf(sedcomment,"(ised=%d/%d)", ised, SEDMODEL.NSURFACE );
  rd_sedFlux(sed_inFile_full, sedcomment
	     ,INPUTS.SIMSED_TREST_RANGE, INPUTS.SIMSED_LAM_RANGE
	     ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL
	     ,SEDMODEL.OPTMASK
             ,&TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY, &TEMP_SEDMODEL.DAYSTEP
             ,&TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM, &TEMP_SEDMODEL.LAMSTEP
             ,TEMP_SEDMODEL.FLUX, TEMP_SEDMODEL.FLUXERR );
  
  N = TEMP_SEDMODEL.NDAY ;
  TEMP_SEDMODEL.MINDAY  = TEMP_SEDMODEL.DAY[0] ;
  TEMP_SEDMODEL.MAXDAY  = TEMP_SEDMODEL.DAY[N-1] ;

  N = TEMP_SEDMODEL.NLAM ;
  TEMP_SEDMODEL.MINLAM  = TEMP_SEDMODEL.LAM[0] ;
  TEMP_SEDMODEL.MAXLAM  = TEMP_SEDMODEL.LAM[N-1] ;


} // end of sedRead



// ********************************
void sedFudge_color(int ised) {

  // Apply color fudge to this SED,  write out new set of SEDs,
  // and update SED.INFO file

  char 
    fnam[] = "sedFudge_color"
    ,sed_outFile[200]
    ,sed_outFile_full[200]
    ,ctmp[80], ctmperr[80]
    ,*ptrsed, snam[20]
    ,outPrefix[100]
    ;

  int 
    icbin, isbin, NcBIN, NsBIN
    ,ipar, iday, ilam, jflux, iday1, ilam1
    ,IDAY_PEAK
    ,IDAY_T15, DOSALT2

    ;
  double 
    cval, sval, stretch, day, lam
    ,flux, fluxerr, x1, Flux_M0, Flux_M1
    ,flux_new, fluxerr_new
    ,cfudge, sfudge, shapecor, polyfudge,  fudge, DIF
    ,daymin, daymin15
    ,Btrans
    ,FBsum0,  FBsum0_new,  FBcor0,  MBcor0, MB0
    ,FBsum15, FBsum15_new, FBcor15, MBcor15
    ,DM15cor
    ;


  FILE *fp_out;

  // -------- BEGIN -----------


  // find IDAY_PEAK closest to peak
  daymin = daymin15 = 99999.0; 
  IDAY_PEAK = IDAY_T15 = -9; 
  for ( iday = 0; iday < TEMP_SEDMODEL.NDAY; iday++ ) {
    day = TEMP_SEDMODEL.DAY[iday];
    if ( fabs(day) < daymin ) {
      IDAY_PEAK = iday ;
      daymin = fabs(day);
    }
    if ( fabs(day-15.0) < daymin15 ) {
      IDAY_T15 = iday ;
      daymin15 = fabs(day-15.0);
    }

  }
  if ( IDAY_PEAK <= 0 ) {
    sprintf(c1err,"Could not find peak day for ised=%d", ised);
    errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
  }

  /*
  printf(" xxxxx ised=%d  IDAY_[PEAK,T15] = %d  %d \n", 
	 ised, IDAY_PEAK, IDAY_T15 );
  */


  DOSALT2 = ( strcmp(INPUTS.SALT2_VERSION_INPUT,"NULL") != 0 );

  ptrsed = SEDMODEL.FILENAME[ised] ;

  if ( DOSALT2 ) { 
    sprintf(outPrefix,"%s", INPUTS.SALT2_VERSION_INPUT); 
    sprintf(snam,"x1");
  }
  else { 
    sprintf(outPrefix,"%s", ptrsed); 
    sprintf(snam,"stretch");
  }

  // loop over color bins
  NcBIN = INPUTS.NBIN_COLORPAR; 
  NsBIN = INPUTS.NBIN_SHAPEPAR; 


  for ( isbin = 1; isbin <= NsBIN; isbin++ ) {

    sval    = INPUTS.SHAPEPAR_VALUE[isbin] ;
    stretch = lcstretch(sval); // conver sval into lc stretch

    if ( NsBIN > 1 ) {
      printf("\t Process %s = %6.3f (%s=%6.3f) \n",
	     INPUTS.SHAPEPAR_NAME, sval, snam, stretch );
      fflush(stdout);
    }

    for ( icbin = 1; icbin <= NcBIN; icbin++ ) {

      cval  = INPUTS.COLORPAR_VALUE[icbin] ;
      
      // create name of new SED file 

      if ( NsBIN <= 1 ) {
	sprintf(sed_outFile, "%s_c%2.2d", outPrefix, icbin );
      }
      else {
	sprintf(sed_outFile, "%s_c%2.2ds%2.2d", outPrefix, icbin, isbin );
      }

      sprintf(sed_outFile_full, "%s/%s", 
	      SIMSED_PATHMODEL_OUTPUT, sed_outFile );

      // open output SED file
      fp_out = fopen(sed_outFile_full, "wt") ;

      NSED_OUTPUT++ ;

      FBsum0  = FBsum0_new  = 0.0;
      FBsum15 = FBsum15_new = 0.0 ;

      for ( iday = 0; iday < TEMP_SEDMODEL.NDAY; iday++ ) {
	
	for ( ilam = 0; ilam < TEMP_SEDMODEL.NLAM; ilam++ ) {

	  day      = TEMP_SEDMODEL.DAY[iday];
	  lam      = TEMP_SEDMODEL.LAM[ilam];
	  flux     = TEMP_SEDSMOOTH.FLUX[iday][ilam] ;
	  fluxerr  = TEMP_SEDSMOOTH.FLUXERR[iday][ilam] ;

	  // get color fudge
	  cfudge  = colorFudge(lam,cval); 

	  // get lambda-polynominal fduge
	  polyfudge = FLUXFUDGE_LAMPOLY(lam);

	  // get shape fudge
	  if ( DOSALT2 ) {
	    Flux_M0   = flux ;
	    Flux_M1   = SALT2_TABLE.SEDFLUX[1][iday][ilam] ;
	    x1        = sval ;
	    flux_new  = Flux_M0 + (x1 * Flux_M1) ;
	    sfudge    = shapeFudge(lam,x1); // correct flux for Phillips
	    shapecor  = flux_new/flux ;
	    polyfudge = 1.0 ;
	  }
	  else if ( INPUTS.NBIN_SHAPEPAR > 1 ) {
	    sfudge   = shapeFudge(lam,sval); 
	    shapecor = stretchFlux(stretch,iday,ilam)/flux ;
	  }
	  else if ( IPAR_S2x1 > 0 ) {
	    x1       = SEDMODEL.PARVAL[ised][IPAR_S2x1];
	    sfudge   = shapeFudge(lam,x1); 
	    shapecor = 1.0 ;
	  }
	  else 
	    { sfudge = shapecor = 1.0 ; }

	  // multiply all fudges into grand fudge
	  fudge   = polyfudge * cfudge * sfudge * shapecor ;

	  /*	  
	  if ( ilam == 100 && ( iday == 20 || iday == 35 ) )
	    printf("\t xxxx day=%5.1f  sfudge = %f  shapecor = %f \n", 
		   day, sfudge, shapecor );
	  */

	  // check for dump option (1st sed at peak only)
	  if ( ised == 1 && 
	       fabs(day+.5)<.6 && 
	       fabs(lam-INPUTS.DUMP_LAMBDA) < 5. ) {
	    printf(" DUMP: day=%5.1f  lam=%6.0f   %s=%6.3f  colorFudge=%8.3f \n", 
		   day, lam, INPUTS.COLORPAR_NAME, cval, fudge);
	  }

	  // apply grand fudge
	  flux_new    = fudge * flux ;
	  fluxerr_new = fudge * fluxerr ;	  

	  // save components of fudges (including extra for debug)
	  TEMP_SEDFUDGE.FLUX[iday][ilam]      = flux_new ;
	  TEMP_SEDFUDGE.FLUXERR[iday][ilam]   = fluxerr_new ;
	  TEMP_SEDFUDGE.ORIG_FLUX[iday][ilam] = flux ;
	  TEMP_SEDFUDGE.POLYFUDGE[iday][ilam] = polyfudge ;
	  TEMP_SEDFUDGE.CFUDGE[iday][ilam]    = cfudge ;
	  TEMP_SEDFUDGE.SFUDGE[iday][ilam]    = sfudge ;
	  TEMP_SEDFUDGE.SHAPECOR[iday][ilam]  = shapecor ;

	  if ( INPUTS.IPAR_PEAKMB > 0 &&  lam > 3500. && lam < 6000. ) {

	    //	    Btrans      = Btransfun(lam); // xxx removed Aug 25 2013
	    Btrans      = filterTrans_BessB(lam);

	    if ( iday == IDAY_PEAK ) {
	      FBsum0     += Btrans * flux ;  // energy-flux sum => no lambda
	      FBsum0_new += Btrans * flux_new ;
	    }
	    if ( iday == IDAY_T15 ) {
	      FBsum15     += Btrans * flux ; 
	      FBsum15_new += Btrans * flux_new ;
	    }
	  }
      
	} // ilam
	

      } // iday
      
   
      sedWrite(fp_out);   // write TEMP_SEDFUDGE to file
      fclose (fp_out);    // close file

      if ( INPUTS.IPAR_PEAKMB > 0 ) {
	FBcor0 = FBsum0_new / FBsum0 ;
	MBcor0 = -2.5*log10(FBcor0);

	FBcor15 = FBsum15_new / FBsum15 ;
	MBcor15 = -2.5*log10(FBcor15);

	DM15cor = (MBcor15 - MBcor0) ;
 
	/*
	printf(" xxxx ised=%2d  %s=%6.3f  MBcor(0,15) = %6.3f,%6.3f  DM15cor=%6.3f\n", 
	       ised, INPUTS.COLORPAR_NAME, cval, MBcor0, MBcor15, DM15cor ) ;
	*/

      }

      // update new info file
      update_SEDINFO_file(fp_sedinfo, sed_outFile, ised,
			  cval, sval, MBcor0, DM15cor );

    } // NBIN loop over color values
  } // NBIN loop over shape values

  fprintf(fp_sedinfo,"\n" );

  fflush(stdout);

} // end of sedFudge_color


// ***********************************
void sedFudge_smear(int ised, int ismear) {

  // October 2011.
  // Apply smearing to sed based on user SMEAROPT name and param-values.
  // Inputs:
  //   * ised refers to input sed that is read in
  //   * ismear is optional smearpar index 
  //
  // Jan 31, 2012: don't let magshift outside MAGSMEAR_MIN/MAX
  // Mar 16 2014: add ismear arg.

  char 
    sed_outFile[MXPATHLEN]
    ,sed_outFile_full[MXPATHLEN]
    ,smearSuffix[60]
    ,fnam[] = "sedFudge_smear" 
    ;

  FILE *fp_out ;
  int indx, iday, ilam, iran, NGAURAN, NFLATRAN,  NFILT, ifilt;

  // ------------ BEGIN -------------

  indx = INPUTS.INDEX_SMEARMODEL;

  // create name of output sed file.

  if ( ismear < 0 ) 
    {  sprintf(smearSuffix,"");  }
  else  {  
    sprintf(smearSuffix,"_%s%2.2d", SMEARPAR_GRID.NAME, ismear+1);  
  }


  sprintf(sed_outFile, "%s%s", SEDMODEL.FILENAME[ised], smearSuffix );


  sprintf(sed_outFile_full, "%s/%s", 
	      SIMSED_PATHMODEL_OUTPUT, sed_outFile );

      // open output SED file
  fp_out = fopen(sed_outFile_full, "wt") ;
  NSED_OUTPUT++ ;


  // generate Gaussian random numbers for this SED 

  init_RANLIST();

  NGAURAN = SMEARMODEL_DEF[indx].NGAURAN ;
  for ( iran=1; iran <= NGAURAN; iran++ ) 
    { GAURANLIST_SMEARMODEL[iran] = GaussRan(1);  }

  NFLATRAN = SMEARMODEL_DEF[indx].NFLATRAN ;
  for ( iran=1; iran <= NFLATRAN; iran++ ) 
    { FLATRANLIST_SMEARMODEL[iran] = FlatRan1(1);  }


  // if using one of the genSmear_model.c models, 
  // set randoms for this SED.
  if ( SMEARMODEL_ID.USE_genSmear ) {
    SETRANGauss_genSmear(NGAURAN,  &GAURANLIST_SMEARMODEL[1]  );
    SETRANFlat_genSmear (NFLATRAN, &FLATRANLIST_SMEARMODEL[1] );
  }


  for ( iday = 0; iday < TEMP_SEDMODEL.NDAY; iday++ ) {
    for ( ilam = 0; ilam < TEMP_SEDMODEL.NLAM; ilam++ ) {
      load_TEMP_SEDFUDGE_SMEAR(iday,ilam,ismear);
    } 
  }  

  // write TEMP_SEDFUDGE to file
  sedWrite(fp_out);  fclose(fp_out);


  NFILT = SMEARMODEL_DEF[indx].NFILTMAG ;
  for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
    sedSynMag(ised,ifilt);  // compute synthetic mag at peak
  }

  // update new info file
  update_SEDINFO_file(fp_sedinfo, sed_outFile, ised,
		      NULLVAL, NULLVAL, NULLVAL, NULLVAL );


  update_SEDINFO_file(fp_sedinfo2, "noSED", ised,
		      NULLVAL, NULLVAL, NULLVAL, NULLVAL );

  fflush(stdout);

} // end of void sedFudge_smear



// =========================================
void load_TEMP_SEDFUDGE_SMEAR(int iday, int ilam, int ismear) {

  // Mar 2014
  // Wrapper to compute magShift from intrinsic scatter model for this
  // input iday and ilam; add intrinsic scatter "magshift" to flux 
  // and load TEMP_SEDFUDGE array.
  // If optional 'ismear' is set, use this index to pick
  // physical SMEARPAR value; e.g., v_Si for VCR model.

  int ONE=1 ;
  int indx, i ;
  double day, lam, flux, fluxerr, magshift, fscale ;
  char fnam[] = "load_TEMP_SEDFUDGE_SMEAR" ;

  // -------------- BEGIN ------------

  indx    = INPUTS.INDEX_SMEARMODEL;
  day     = TEMP_SEDMODEL.DAY[iday];
  lam     = TEMP_SEDMODEL.LAM[ilam];
  flux    = TEMP_SEDSMOOTH.FLUX[iday][ilam] ;
  fluxerr = TEMP_SEDSMOOTH.FLUXERR[iday][ilam] ;
  
  if ( indx == SMEARMODEL_ID.NOSMEAR )
    { magshift = 0.0 ; } 
  else if ( indx == SMEARMODEL_ID.COHERENT ) 
    { magshift = SMEARMAG_COHERENT(lam,day); }
  
  else if ( indx == SMEARMODEL_ID.SINX ) 
    { magshift = SMEARMAG_SINX(lam,day); }
  
  else if ( indx == SMEARMODEL_ID.KRW09 ) 
    { magshift = SMEARMAG_KRW09(lam,day); }
  
  else if ( indx == SMEARMODEL_ID.G10 ) 
    {  get_genSmear_SALT2(day, ONE, &lam, &magshift);   }
  
  else if ( indx == SMEARMODEL_ID.C11 ) 
    { get_genSmear_Chotard11(day, ONE, &lam, &magshift); }
  
  else if ( indx == SMEARMODEL_ID.VCR )  {
    GENSMEAR_VCR.USE        = 2 ;  // flag to use fixed value below
    GENSMEAR_VCR.GENRAN_VSI = SMEARPAR_GRID.VALUE[ismear] ;
    get_genSmear_VCR(day, ONE, &lam, &magshift);
  }
  else {
    sprintf(c1err,"No function for SMEARMODEL = '%s'",
	    SMEARMODEL_DEF[indx].NAME);
    sprintf(c2err,"Either add function of change SMEARMODEL");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 		
  }
  
  // abort on crazy mag shift that give large fluxes
  if ( -magshift > MAGSMEAR_ABORT ) {
    printf("\n");
    printf(" ============================================= \n");
    printf(" PRE-ABORT DUMP of randoms: \n");
    
    int NGAURAN = SMEARMODEL_DEF[indx].NGAURAN ;
    int NFLATRAN = SMEARMODEL_DEF[indx].NFLATRAN ;
    for ( i=1; i <= NGAURAN; i++ ) 
      {  printf("\t GAURAN(%d) = %f \n", i, GAURANLIST_SMEARMODEL[i] ); }
    
    for ( i=1; i <= NFLATRAN; i++ ) 
      {  printf("\t FLATRAN(%d) = %f \n", i, FLATRANLIST_SMEARMODEL[i] ); }
    
    sprintf(c1err,"Crazy magshift = %f > MAGSMEAR_ABORT=%5.2f",
	    magshift,  MAGSMEAR_ABORT );
    sprintf(c2err,"at lam=%7.1f  day=%5.2f", lam, day);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  if ( magshift > MAGSMEAR_MAX ) { magshift = MAGSMEAR_MAX ; }
  if ( magshift < MAGSMEAR_MIN ) { magshift = MAGSMEAR_MIN ; }
  
  fscale = 1.0/pow(10.0,(0.4*magshift));
  TEMP_SEDFUDGE.FLUX[iday][ilam]    = fscale*flux ;
  TEMP_SEDFUDGE.FLUXERR[iday][ilam] = fscale*fluxerr ;
  
} // end of  load_TEMP_SEDFUDGE_SMEAR

// =========================================
void  sedSynMag(int ised, int ifilt) {

  // Compute synethetic mag at peak for this ifilt.
  // Fill SYNMAGDIF_SMEARMODEL[ifilt]  = mag(fudge) - mag(orig)

  int indx, NLAM_FILT, NLAM_SED, ilam, IDAY_PEAK, jflux ;
  int OPT_INTERP = 1; 
  char cfilt[2];
  char transFile[MXPATHLEN];
  char fnam[] = "sedSynMag";

  double 
    lam, trans, flux_orig, flux_fudge
    ,FSUM_ORIG, FSUM_FUDGE
    ,LAMMAX_FILT
    ,LAMMIN_FILT
    ,LAM_FILT[MXLAM_FILT]
    ,TRANS_FILT[MXLAM_FILT]
    ,DAYSTEP, DIF
    ,LAMSTEP
    ;


  // -------------- BEGIN -------------

  // get character for this filter
  indx = INPUTS.INDEX_SMEARMODEL;
  sprintf(cfilt,"%c", SMEARMODEL_DEF[indx].FILT[ifilt] );

  LAMSTEP    = TEMP_SEDMODEL.LAMSTEP ;   
  DAYSTEP    = TEMP_SEDMODEL.DAYSTEP ;   

  // set filter trans-file based on filter
  if ( strcmp(cfilt,"U") == 0 ||
       strcmp(cfilt,"B") == 0 ||
       strcmp(cfilt,"V") == 0 ||
       strcmp(cfilt,"R") == 0 ||
       strcmp(cfilt,"I") == 0 ) 
    {
    sprintf(transFile,"%s/filters/Bessell90/Bessell90_K09/Bessell90_%s.dat", 
	    PATH_SNDATA_ROOT, cfilt) ;
  }
  else {
    // abort
  }

  // determine day-index at peak
  DIF        = - TEMP_SEDMODEL.DAY[0] ;   
  IDAY_PEAK  = (int)(DIF/DAYSTEP) ;


  // read filter trans vs. wavelength
  rd2columnFile(transFile, MXLAM_FILT, &NLAM_FILT, LAM_FILT, TRANS_FILT);

  LAMMIN_FILT = LAM_FILT[0];
  LAMMAX_FILT = LAM_FILT[NLAM_FILT-1];

  NLAM_SED = TEMP_SEDMODEL.NLAM ;
  FSUM_ORIG = FSUM_FUDGE = 0.0 ;

  for ( ilam=0 ; ilam < NLAM_SED; ilam++ ) {
    lam = TEMP_SEDMODEL.LAM[ilam] ;

    if ( lam < LAMMIN_FILT ) { continue ; }
    if ( lam > LAMMAX_FILT ) { continue ; }

    jflux      = NLAM_SED*IDAY_PEAK + ilam ;
    flux_orig  = TEMP_SEDMODEL.FLUX[jflux] ;
    flux_fudge = TEMP_SEDFUDGE.FLUX[IDAY_PEAK][ilam] ;

    // get transmission from linear interp
    trans = interp_1DFUN(OPT_INTERP, lam, 
			 NLAM_FILT, LAM_FILT, TRANS_FILT, fnam);

    FSUM_ORIG  += ( flux_orig  * trans );
    FSUM_FUDGE += ( flux_fudge * trans );
  }


  SYNMAGDIF_SMEARMODEL[ifilt] = -2.5*log10(FSUM_FUDGE/FSUM_ORIG) ;

  
} // end of sedSynMag


// **********************************
void update_SEDINFO_file(FILE *fp, char *sed_outFile, 
			 int ised, double cval, double sval,
			 double MB0cor, double DM15cor ) {

  int ipar, iran, indx;
  double MB0, DM15;

  char fnam[] = "update_SEDINFO_file";

  // ---------------- BEGIN ----------------

  if ( strcmp(sed_outFile,"noSED") == 0 ) 
    { fprintf(fp,"SN: %4.4d  ", ised );  }
  else 
    { fprintf(fp,"SED: %s  ",  sed_outFile );  }


  for( ipar=1; ipar <= SEDMODEL.NPAR; ipar++ ) {

    if ( ipar == INPUTS.IPAR_PEAKMB ) {
      MB0 = SEDMODEL.PARVAL[ised][ipar] + MB0cor ;
      fprintf(fp, "%6.3f ", MB0 ) ;
      continue ;
    }
    if ( ipar == INPUTS.IPAR_DM15 ) {
      DM15 = SEDMODEL.PARVAL[ised][ipar] + DM15cor ;
      fprintf(fp, "%6.3f ", DM15 ) ;
      continue ;
    }

    fprintf(fp, "%s ", SEDMODEL.PARVAL_STRING[ised][ipar] );

  }


  if ( INPUTS.NBIN_COLORPAR >1 ) 
    { fprintf(fp,"%6.3f ", cval ); } // tack on color param

  if ( INPUTS.NBIN_SHAPEPAR > 1 ) 
    { fprintf(fp,"%6.3f ", sval ); } // tack on shape param

  // tack on randoms for SMEARING
  indx = INPUTS.INDEX_SMEARMODEL ; 
  if ( indx > 0 ) {
    for ( iran = 1; iran <= SMEARMODEL_DEF[indx].NGAURAN; iran++ ) 
      {  fprintf(fp,"%6.3f ", GAURANLIST_SMEARMODEL[iran] ); }

    for ( iran = 1; iran <= SMEARMODEL_DEF[indx].NFLATRAN; iran++ ) 
      {  fprintf(fp,"%6.3f ", FLATRANLIST_SMEARMODEL[iran] ); }
  }


  // tack on synthetic mags
  int NFILT, ifilt ;
  double mag;
  NFILT = SMEARMODEL_DEF[indx].NFILTMAG ;
  for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
    mag = SYNMAGDIF_SMEARMODEL[ifilt] ;
    fprintf(fp,"%6.3f ", mag );    
  }

  // tack on KRW09 variables

  if ( KRW09.N_COSANGLE > 0 ) {

    fprintf(fp,"%6.3f %6.3f %8.3f ", 
	    KRW09.COSANGLE_RAN, KRW09.DM15_RAN, KRW09.MB_RAN );
  }


  // tack on special/extra params
  if ( INPUTS.INDEX_SMEARMODEL == SMEARMODEL_ID.VCR ) 
    { fprintf(fp," %.2f ", GENSMEAR_VCR.GENRAN_VSI );  }


  // always tack on <CR>
  fprintf(fp, " \n");

  fflush(fp);

} // end of update_SEDINFO_file


// **************************************************
double SMEARMAG_COHERENT(double lam, double day) {

  double magshift, ran, P1 ;
  char fnam[] = "SMEARMAG_COHERENT" ;

  // -------------- BEGIN -----------------
  ran = GAURANLIST_SMEARMODEL[1] ;
  P1  = INPUTS.SMEARMODEL_PARVAL[1];
  magshift = ran*P1;
  return magshift ;
} // end of SMEARMAG_COHERENT

// **************************************************
double SMEARMAG_SINX(double lam, double day) {

  // Sinusoidal mag-smearing vs. wavelength.
  // 

  double magshift, ran, P1, P2, P3, arg, RAD ;
  char fnam[] = "SMEARMAG_SINX" ;

  // -------------- BEGIN -----------------

  ran = GAURANLIST_SMEARMODEL[1] ;
  P1  = INPUTS.SMEARMODEL_PARVAL[1];
  P2  = INPUTS.SMEARMODEL_PARVAL[2];
  P3  = INPUTS.SMEARMODEL_PARVAL[3];

  arg      = TWOPI*(lam-P2)/P3 ;
  magshift = ran * P1 * sin(arg);

  return magshift ;

} // end of SMEARMAG_SINX

// **************************************************
double SMEARMAG_KRW09(double lam, double day) {
  
  // Dec 2011 RK
  // Function to smear spectra using KRW09 model.
  // Flux-ratio shift is FLUX(cosangle)/FLUX(cosangle=0).
  //
  // Jan 31, 2012:  require FLUX>0 and FLUX0 >0 before taking
  //                log10(RATIO)

  int 
    NLAM, NCOS, icos, i, IMODEL
    ,icos_near[2], icos_ref[2]
    ,ilam, iday, jflux
    ,LDMP_JPM_COMPARE = 0
    ,LDMP1, LDMP2, LLAM
    ,LDMP
    ,IMODEL_AVG = 0 
    ,ICOS_AVG   = 0 
    ;

  double 
    ran, magshift, cosangle, RANGE, A, x
    ,COSANGLE_GRID, COSANGLE_RAN, DIF, DIFMIN, dif
    ,FLUXREF[2], FLUX0, FLUXCOS[2], FLUX, RATIO_AVG, RATIO_0
    ,FLUXAVG
    ,DM15SED[2], MBSED[2]
    ;

  char fnam[] = "SMEARMAG_KRW09";

  // -------------- BEGIN -----------

  magshift = 0.0 ;

  // select random model (DC,IGNIT values)
  ran    = FLATRANLIST_SMEARMODEL[1] ;
  x      = ran * ((double)KRW09.NMODELS - 1.0E-9) ;
  IMODEL = (int)(x) + 1 ;
  if ( IMODEL < 1  ||  IMODEL > KRW09.NMODELS ) {
    sprintf(c1err,"Invalid IMODEL=%d (ran=%le  x=%le)", 
	    IMODEL, ran, x );
    sprintf(c2err,"Valid IMODELS are 1 - %d", 
	    KRW09.NMODELS);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  KRW09.IMODEL_RAN = IMODEL ;

  // pick random cos(theta)
  ran          = FLATRANLIST_SMEARMODEL[2] ;
  RANGE        = KRW09.COSANGLE_MAX - KRW09.COSANGLE_MIN ;
  COSANGLE_RAN = KRW09.COSANGLE_MIN + (ran * RANGE) ;
  KRW09.COSANGLE_RAN = COSANGLE_RAN ;

  // find closest COSANGLE

  if ( ran != KRW09.TEMP_RANFLAT ) {
    NCOS = KRW09.N_COSANGLE ;
    DIFMIN = 9999. ;  icos_near[0] = icos_near[1] = -9;
    for ( icos=1; icos <= NCOS; icos++  ) {
      COSANGLE_GRID = KRW09.COSANGLE_LIST[icos] ;
      DIF = fabs(COSANGLE_GRID - COSANGLE_RAN);
      if ( DIF < DIFMIN )  { DIFMIN = DIF ;   icos_near[0] = icos; }
    } // icos

    if ( icos_near[0] <= 0 ) {
      sprintf(c1err,"Could not find nearest COSANGLE bin");
      sprintf(c2err,"for ran=%6.4f  COSANGLE=%6.3f", ran, COSANGLE_RAN);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    COSANGLE_GRID = KRW09.COSANGLE_LIST[icos_near[0]] ;
    if ( COSANGLE_RAN < COSANGLE_GRID && icos_near[0] > 1 ) 
      { icos_near[0]-- ; }

    icos_near[1] = icos_near[0] + 1;
    
    KRW09.TEMP_ICOS_NEAR[0] = icos_near[0] ;
    KRW09.TEMP_ICOS_NEAR[1] = icos_near[1] ;
  }


  KRW09.TEMP_RANFLAT = ran ; 

  for ( i=0; i <=1 ; i++ ) {
    icos_ref[i]   = KRW09.ICOS_NEAR90DEG[i] ;  
    icos_near[i]  = KRW09.TEMP_ICOS_NEAR[i] ;

    // near cosang=0
    FLUXREF[i] = interp_SEDFLUX_KRW09(IMODEL, icos_ref[i],  lam, day); 
    FLUXCOS[i] = interp_SEDFLUX_KRW09(IMODEL, icos_near[i], lam, day);

    DM15SED[i] = KRW09.DM15_LIST[IMODEL][icos_near[i]] ;
    MBSED[i]   = KRW09.MB_LIST[IMODEL][icos_near[i]] ;
  }


  // interpolate in cos-angle space to get flux at 90 deg.
  dif = 0.0 - KRW09.COSANGLE_LIST[icos_ref[0]];
  DIF = 
    KRW09.COSANGLE_LIST[icos_ref[1]] - 
    KRW09.COSANGLE_LIST[icos_ref[0]] ;
  FLUX0   = FLUXREF[0] + (dif/DIF)*(FLUXREF[1]-FLUXREF[0]);


  // get average flux for this lam and day
  FLUXAVG = interp_SEDFLUX_KRW09(IMODEL_AVG, ICOS_AVG, lam, day); 
  
  // interpolate FLUX at COSANGLE
  dif = COSANGLE_RAN - KRW09.COSANGLE_LIST[icos_near[0]];
  DIF = 
    KRW09.COSANGLE_LIST[icos_near[1]] - 
    KRW09.COSANGLE_LIST[icos_near[0]];
  FLUX     = FLUXCOS[0] + (dif/DIF)*(FLUXCOS[1]-FLUXCOS[0]);

  // also interpolate DM15 and MB (Jun 29, 2012)
  KRW09.DM15_RAN = DM15SED[0] + (dif/DIF) * ( DM15SED[1] - DM15SED[0] );
  KRW09.MB_RAN   = MBSED[0]   + (dif/DIF) * ( MBSED[1]   - MBSED[0]   );

  // convert FLUX/FLUXREF to a change in mag
  RATIO_0 = RATIO_AVG = 1.0 ;  // default is no KRW09 info
  if ( FLUXAVG > 0.0 && FLUX > 0.0 )  { RATIO_AVG = FLUX / FLUXAVG ; }
  if ( FLUX0   > 0.0 && FLUX > 0.0 )  { RATIO_0   = FLUX / FLUX0   ; }

  A        = INPUTS.SMEARMODEL_PARVAL[2] ; // user ampl of perturbation
  magshift = -2.5*log10(RATIO_AVG) ;

  /*
  if ( fabs(day-0.1)  < 0.5  && 
       fabs(lam-5000) < 20.0 && 
       FLUX0 > 0 && FLUXAVG > 0  ) {
    printf(" 000000  %d  %6.3f  %6.2f  %7.4f  %7.4f  %f\n", 
	   IMODEL, COSANGLE_RAN, day, lam, RATIO_AVG, RATIO_0 ) ;
	   }  */ 
  

  LDMP = ( lam == -2160.0 && day == 20.0 ) ;
  if ( LDMP ) {
    printf("\n xxx ------------ %s DUMP ------------------ \n", fnam);
    printf(" xxx lam = %f  day=%f \n", lam, day);
    printf(" xxx FLUXREF[0,1] = %le , %le \n", FLUXREF[0], FLUXREF[1] );
    printf(" xxx dif=%f  DIF=%f \n", dif, DIF);
    printf(" xxx FLUX=%le  FLUX0=%le \n", FLUX, FLUX0 );
    printf(" xxx RATIO_AVG=%le  A=%f \n", RATIO_AVG, A);
    printf(" xxx magshift = %f \n", magshift );
  }

  // check dump-option to compare with Marriner's stand-alone code
  if ( LDMP_JPM_COMPARE ) {
    LLAM = fabs(lam-4000.0) < 5.0 ;

    LDMP1 = ( LLAM && fabs(day-2.0) < .2  );
    if ( LDMP1 ) 
      { printf(" 111111  %6.3f  %6.3f \n", COSANGLE_RAN, RATIO_AVG); }

    LDMP2 = ( LLAM && COSANGLE_RAN > 0.90 );
    if ( LDMP2 ) 
      { printf(" 222222  %6.3f  %6.3f \n", day, RATIO_AVG); }

  }

  return  magshift ;

} // end of SMEARMAG_KRW09


// **********************************
double interp_SEDFLUX_KRW09(int imodel, int icos, double lam, double day) {

  // interp flux to the input lam, day.
  // icos & imodel just fixes the model (COSANGLE, DC,IGNIT values)

  double 
    LAMSTEP, DAYSTEP
    ,LAMGRID, DAYGRID
    ,frac_DAY, frac_LAM
    ,CORNER_FLUX[2][2]
    ,CORNER_WGT[2][2]
    ,CORNER_WGTSUM
    ,FSUM, FLUX
    ;

  int ilam, iday, NLAM, NDAY, iiday, iilam, jflux ;
  int ERRFLAG ;

  char fnam[] = "interp_SEDFLUX_KRW09";

  // ------------- BEGIN -----------


  LAMSTEP = KRW09.LAMSTEP[imodel][icos];   
  NLAM    = KRW09.NLAM[imodel][icos] ;
  DAYSTEP = KRW09.DAYSTEP[imodel][icos];   
  NDAY    = KRW09.NDAY[imodel][icos] ;

  // get ilam and iday bins
  get_intBins(imodel, icos, lam, day, &ilam, &iday, fnam) ; 

  // get flux at each corner of the DAY vs LAMBDA square
  CORNER_WGTSUM = 0.0 ;
  FSUM          = 0.0 ;
  ERRFLAG = 0;

  for ( iilam=0; iilam <=1; iilam++ ) {
    for ( iiday=0; iiday <=1; iiday++ ) {

      jflux   = NLAM*(iday+iiday) + (ilam+iilam) ;

      DAYGRID = KRW09.DAY[imodel][icos][iday+iiday] ;
      LAMGRID = KRW09.LAM[imodel][icos][ilam+iilam] ;

      // determine relative distance (fraction) from this grid point;
      // 0 => right on the grid point and 1 is max dist away
      frac_DAY = fabs(day - DAYGRID) / DAYSTEP ;
      frac_LAM = fabs(lam - LAMGRID) / LAMSTEP ;

      if ( frac_DAY > 1.0 || frac_DAY < 0 ) {
	ERRFLAG = 1;
	sprintf(c1err,"invalid frac_DAY = %f (should be 0-1)", frac_DAY);
	sprintf(c2err,"for day(inp,grid)=%6.2f,%6.2f iday=%d iiday=%d",
        	day, DAYGRID, iday, iiday);
	errmsg(SEV_WARN, 0, fnam, c1err, c2err );
      }
      if ( frac_LAM > 1.0 || frac_LAM < 0 ) {
	ERRFLAG = 1 ;
	sprintf(c1err,"invalid frac_LAM = %f (should be 0-1)", frac_LAM);
	sprintf(c2err,"for lam(inp,grid)=%6.1f,%6.1f ilam=%d iilam=%d",
		lam, LAMGRID, ilam, iilam);
	errmsg(SEV_WARN, 0, fnam, c1err, c2err );
      }

      if ( ERRFLAG ) {
	sprintf(c1err,"More ABORT info:" );
	sprintf(c2err,"imodel=%d icos=%d  LAM0=%7.1f  DAY0=%6.2f"
		, imodel, icos
		, KRW09.LAM[imodel][icos][0]
		, KRW09.DAY[imodel][icos][0] );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );	       
      }

      CORNER_FLUX[iiday][iilam] = KRW09.FLUX[imodel][icos][jflux];

      CORNER_WGT[iiday][iilam] = 
	pow( (1.0 - frac_DAY), 2.0) + 
	pow( (1.0 - frac_LAM), 2.0) ;

      CORNER_WGTSUM += CORNER_WGT[iiday][iilam];
      FSUM += (CORNER_FLUX[iiday][iilam] * CORNER_WGT[iiday][iilam]);

    } // iiday
  } // iilam


  FLUX = FSUM / CORNER_WGTSUM ;

  return FLUX ;

} // end of interp_SEDFLUX_KRW09


// ********************************************************
void get_intBins(int imodel, int icos, double lam, double day,
		 int *ilam, int *iday, char *abortComment) {

  // for input lam and day, return integer bins *ilam and *iday
  // *abortComment is for ABORT message.

  double LAMSTEP, DAYSTEP, DIF_LAM, DIF_DAY ;
  int    NLAM, NDAY ;
  
  char fnam[] = "get_intBins";

  // ---------------- BEGIN ---------------


  LAMSTEP = KRW09.LAMSTEP[imodel][icos];   
  DAYSTEP = KRW09.DAYSTEP[imodel][icos];   
  NLAM    = KRW09.NLAM[imodel][icos] ;
  NDAY    = KRW09.NDAY[imodel][icos] ;

  /*
  printf("# ----------------------------------------- \n");
  printf(" YOOOOO %s inputs are imod=%d icos=%d  lam=%f day=%f \n",
	 fnam, imodel, icos, lam, day); fflush(stdout);
  printf(" YOOOOO LAMSTEP=%f  DAYSTEP=%f  LAM0=%f  DAY0=%f \n"
	 ,LAMSTEP, DAYSTEP
	 , KRW09.LAM[imodel][icos][0]
	 , KRW09.DAY[imodel][icos][0] );
  */

  DIF_LAM     = lam - KRW09.LAM[imodel][icos][0] ;
  *ilam        = (int)(DIF_LAM/LAMSTEP) ;

  DIF_DAY     = day - KRW09.DAY[imodel][icos][0] ;
  *iday        = (int)(DIF_DAY/DAYSTEP) ;

  // error checking 
  if ( DIF_DAY < 0.0 ) {
    sprintf(c1err,"day=%6.2f is before EPMIN=%6.2f (icos=%d,imodel=%d)",
	    day, KRW09.DAY[imodel][icos][0], icos, imodel );
    sprintf(c2err,"Try restricting EPOCH_RANGE key in input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  int ERRFLAG = 0;
  int ibin, NBIN ;
  char ctmp[8];

  if ( *ilam < 0 || *ilam >= NLAM ) 
    { ERRFLAG = 1; ibin= *ilam ; NBIN=NLAM ; sprintf(ctmp,"ilam"); }
  if ( *iday < 0 || *iday >= NDAY ) 
    { ERRFLAG = 1; ibin= *iday ; NBIN=NDAY ; sprintf(ctmp,"iday"); }

  if ( ERRFLAG >  0 ) {  
    sprintf(c1err,"%s=%d is out of range (0 to %d), '%s'",  
	    ctmp, ibin, NBIN, abortComment);
    sprintf(c2err,"for lam = %8.2f  day=%7.2f  icos=%d  imodel=%d", 
	    lam, day, icos, imodel);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    
  } // ERRFLAG

} // end of get_intBins


// ***********************************
void init_SMEARMODEL(void) {

  double SMEAR_SCALE = 1.0 ;
  int indx, NGAURAN, NFLATRAN, OPT,  NPAR, ipar ;
  double sigcoh = -9.0 ;
  char fnam[] = "init_SMEARMODEL" ;

  // ------------------- BEGIN -------------------

  // always init few things 
  KRW09.N_COSANGLE = 0;

  // get smear-model index
  indx = INPUTS.INDEX_SMEARMODEL ; 

  // check for special init functions
  if ( indx == SMEARMODEL_ID.KRW09 )
    { init_SMEARMAG_KRW09(); }
  else if ( indx == SMEARMODEL_ID.G10 ) {
    char smearFile[200] = "" ;
    init_genSmear_FLAGS(SMEAR_SCALE); // internal inits
    init_genSmear_SALT2("SALT2.Guy10", "G10", sigcoh) ; 
    get_NRAN_genSmear(&NGAURAN, &NFLATRAN); // Jan 2014, RK
    SMEARMODEL_DEF[indx].NGAURAN  = NGAURAN ; 
    SMEARMODEL_ID.USE_genSmear = 1 ;
  }
  else if ( indx == SMEARMODEL_ID.C11 ) {
    init_genSmear_FLAGS(SMEAR_SCALE); // internal inits
    OPT     = (int)INPUTS.SMEARMODEL_PARVAL[1] ;
    init_genSmear_Chotard11(OPT) ;
    get_NRAN_genSmear(&NGAURAN, &NFLATRAN); // Jan 2014, RK
    SMEARMODEL_DEF[indx].NGAURAN  = NGAURAN ; 
    SMEARMODEL_ID.USE_genSmear = 1 ;
  }
  else if ( indx == SMEARMODEL_ID.VCR ) {
    init_genSmear_FLAGS(SMEAR_SCALE); // internal inits
    init_genSmear_VCR("VCR", MODEL_SALT2) ;
    get_NRAN_genSmear(&NGAURAN, &NFLATRAN); 
    SMEARMODEL_DEF[indx].NGAURAN   = NGAURAN ; 
    SMEARMODEL_DEF[indx].NFLATRAN  = NFLATRAN ; 
    SMEARMODEL_ID.USE_genSmear = 1 ;

    // load VSI grid -> extra SIMSED dimension
    double VSI, VSI_MIN, VSI_MAX, VSI_BINSIZE = 0.2 ;
    VSI_MIN = GENSMEAR_VCR.VSI[0] ;
    VSI_MAX = GENSMEAR_VCR.VSI[GENSMEAR_VCR.NBIN_VSI-1] ;
    load_SMEARPAR_GRID("VSI", indx, VSI_MIN, VSI_MAX, VSI_BINSIZE );
  }

  // print summary (Aug 31, 2012)
  NPAR =  SMEARMODEL_DEF[indx].NPAR ;
  printf("\n ====== SMEAR-MODEL SUMMARY ============ \n");
  printf(" Selected SMEAROPT: %s \n",  SMEARMODEL_DEF[indx].NAME ) ;
  printf("\t (%s) \n", SMEARMODEL_DEF[indx].COMMENT );
  for ( ipar=1; ipar <= NPAR; ipar++ ) {
    printf("\t User SmearParam(%2.2d) = '%s' (atof = %f)\n"
	   ,ipar
	   ,INPUTS.SMEARMODEL_PARSTRING[ipar] 
	   ,INPUTS.SMEARMODEL_PARVAL[ipar] 
	   );    
  }
  printf(" Number of Gaussian  randoms : %d \n", 
	 SMEARMODEL_DEF[indx].NGAURAN );
  printf(" Number of flat(0-1) randoms : %d \n", 
	 SMEARMODEL_DEF[indx].NFLATRAN );
  printf("\n ==================================================== \n\n") ;


  return ;

} // end of init_SMEARMODEL


// ======================================
void load_SMEARPAR_GRID(char *parName, int indx,
			double parMin, double parMax, double parBin ) {

  // Mar 2014
  // load SMEARMODEL_DEF[indx].SMEARPAR_GRID array.

  int    NBIN ;
  double parVal ;
  char   fnam[] = "load_SMEARPAR_GRID" ;

  // ---------- BEGIN ------------------

  sprintf(SMEARPAR_GRID.NAME,"%s", parName);

  NBIN = 0 ;
  for(parVal=parMin; parVal <= parMax; parVal += parBin ) {
    if ( NBIN < MXBIN_SMEARPAR ) 
      { SMEARPAR_GRID.VALUE[NBIN] = parVal ; }
    NBIN++ ;
  }

  if ( NBIN >= MXBIN_SMEARPAR ) {
    sprintf(c1err, "NBIN_SMEARPAR=%d exceeds bound of MXBIN_SMEARPAR=%d", 
	    NBIN, MXBIN_SMEARPAR);
    sprintf(c2err, "%s  range=%f to %f with %f binsize",
	    parName, parMin, parMax, parBin);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  SMEARPAR_GRID.NBIN = NBIN ;
  printf("\n\t !!! %s SMEARPAR-Grid for SIMSED: %.2f to %.2f (%d bins) !!!\n"
	 ,parName
	 ,SMEARPAR_GRID.VALUE[0]
	 ,SMEARPAR_GRID.VALUE[NBIN-1], NBIN );
  
} // end of load_SMEARPAR_GRID



#ifdef CERNLIB
// ======================================
void mktup_KRW09(void) {

  // XXXXXX use sntools_output XXXXXXXXXX
  // make monitor ntuple (lam,day,magshift) for smear model.

  int  
    hlim, ivar, hid
    ,NVAR   = 7
    ,NTID   = 7788
    ,NBUFF  = 1024*NVAR
    ,LREC   = 1024
    ,LUNHIS = 20
    ,ICYC   = 1
    ,IERR
    ,ilam, iday, indx, iran
    ,IMODEL, ICOS
    ;

  float XVAR[20];

  char 
    ntfile[200]
    ,ntname[20]   = "SMEAR"
    ,topdir[20]   = "snana"
    ,tmpdir[20]   = "          "
    ,copt[8]      = "PN " // allow capital letters, new file
    ,blank[4]     = "  "
    ,R[4]         = "R " 
    ,fnam[]     = "mktup_KRW09" 
    ;

#define NRANTUP 5
  double 
    lam8, day8, mag8, MINLAM8, MAXLAM8
    ,RANLIST[NRANTUP] = {0.001, 0.25, 0.500, 0.75, 0.999 }
  ;

  // ------------- BEGIN ------------

  sprintf(ntfile,"%s/SMEAROPT_MON.TUP", SIMSED_PATHMODEL_OUTPUT );
  // sprintf(ntfile,"SMEAROPT_MON.TUP" );

  printf(" Create monitor ntuple: %s \n", ntfile);

  sprintf ( CTAGNTUP[0], "IMODEL"   ); // specifies DC/IGNIT values
  sprintf ( CTAGNTUP[1], "COSANG"   );
  sprintf ( CTAGNTUP[2], "DM15"     );
  sprintf ( CTAGNTUP[3], "MB"       );
  sprintf ( CTAGNTUP[4], "LAM"      );
  sprintf ( CTAGNTUP[5], "DAY"      );
  sprintf ( CTAGNTUP[6], "MSMEAR"   );

  for ( ivar=0; ivar < NVAR; ivar++ ) { 
    ptr_CTAGNTUP[ivar] = CTAGNTUP[ivar]; 
    
    /*
    printf(" xxx load CTAGNTUP[%d] = '%s'   ptr = (%s) \n", 
	   ivar, CTAGNTUP[ivar], ptr_CTAGNTUP[ivar] );
    fflush(stdout); 
    */
  }

  hlim = 2000000 ;
  hlimit_(&hlim);
  hropen_(&LUNHIS, topdir, ntfile, copt, &LREC, &IERR,
	  strlen(topdir), strlen(ntfile), strlen(copt) );

  hcdir_(tmpdir, R, strlen(tmpdir), strlen(R) );

  hbookn_(&NTID, ntname, &NVAR, topdir, &NBUFF, ptr_CTAGNTUP[0],
	  strlen(ntname), strlen(topdir), 8);

  indx = INPUTS.INDEX_SMEARMODEL;
  MINLAM8 = INPUTS.SIMSED_LAM_RANGE[0] ;
  MAXLAM8 = INPUTS.SIMSED_LAM_RANGE[1] ;

  for ( iran=0; iran < NRANTUP; iran++ ) {

    FLATRANLIST_SMEARMODEL[1] = RANLIST[iran] ; // for ran IMODEL
    FLATRANLIST_SMEARMODEL[2] = RANLIST[iran] ; // for ran cosangle

    for ( day8=-10.0; day8 <= 10.0; day8+=10.0 ) {
      for ( lam8 = MINLAM8; lam8 < MAXLAM8; lam8+=100. ) {
	
	mag8    =  SMEARMAG_KRW09(lam8,day8); 
	IMODEL  = KRW09.IMODEL_RAN ;
	ICOS    =  KRW09.TEMP_ICOS_NEAR[0] ;

	XVAR[0] = (float)IMODEL ;
	XVAR[1] = (float)KRW09.COSANGLE_RAN ;
	XVAR[2] = (float)KRW09.DM15_LIST[IMODEL][ICOS] ;
	XVAR[3] = (float)KRW09.MB_LIST[IMODEL][ICOS] ;
	XVAR[4] = (float)lam8 ;
	XVAR[5] = (float)day8 ;
	XVAR[6] = (float)mag8 ;
	hfn_( &NTID, XVAR ); 
	  
      } // ilam
    } // iday
  } // iran


  hid = NTID ;
  hrout_(&hid, &ICYC, blank, strlen(blank) );
  hrend_(topdir, strlen(topdir) );

} // end of mktup_KRW09
#endif

// **************************************************
void init_SMEARMAG_KRW09(void) {

  int indx, NCOS, icos, imodel ;
  char 
     *name
    ,*modelString
    ,*par1
    ,fnam[] = "init_SMEARMAG_KRW09" 
    ;

  FILE *fp;

  // --------------- BEGIN ----------------

  indx = SMEARMODEL_ID.KRW09 ;
  name = SMEARMODEL_DEF[indx].NAME ;
  par1 = INPUTS.SMEARMODEL_PARSTRING[1] ;

  printf("\n");
  printf(" ##################################################### \n");
  printf(" Init %s (%s) \n", name, par1);


  if ( strcmp(par1,"Blondin8") == 0 ) {
    sprintf(KRW09.DIR,"%s/models/SIMSED/%s_%s", PATH_SNDATA_ROOT, name, par1);
    KRW09.SUPER8FLAG = 1 ;
  }
  else {
    sprintf(KRW09.DIR,"%s/models/SIMSED/%s", PATH_SNDATA_ROOT, name);
    KRW09.SUPER8FLAG = 0 ;
  }

  printf("    Read model from : %s \n", KRW09.DIR );

  KRW09.SEDRANGE_DAY[0] = INPUTS.SIMSED_TREST_RANGE[0] - 3.0 ;
  KRW09.SEDRANGE_DAY[1] = INPUTS.SIMSED_TREST_RANGE[1] + 3.0 ;

  KRW09.SEDRANGE_LAM[0] = INPUTS.SIMSED_LAM_RANGE[0] - 40.0 ;
  KRW09.SEDRANGE_LAM[1] = INPUTS.SIMSED_LAM_RANGE[1] + 40.0 ;

  // read SED.INFO file
  read_SEDINFO_KRW09();

  NCOS = KRW09.N_COSANGLE ;

  // allocate memory and read SEDs
  for ( imodel=1; imodel <= KRW09.NMODELS; imodel++ ) {
    for ( icos=1; icos <= NCOS; icos++ ) 
      { read_SED_KRW09(imodel,icos);  }
  }

  // smooth the SEDs since they can be noisy
  printf("\n    Smooth the %s  SEDs ... \n", name );
  for ( imodel=1; imodel <= KRW09.NMODELS; imodel++ ) {
    for ( icos=1; icos <= NCOS; icos++ ) 
      {  smooth_SED_KRW09(imodel,icos); }
  }


  // fill FLUXAVG array with average spectral flux
  fluxavg_KRW09();

  //  mktup_KRW09();

  printf(" Done initializing %s \n\n", name);


} // end of  init_SMEARMAG_KRW09


// ----------------------
void read_SEDINFO_KRW09(void) {
  
  // read SED.INFO file for KRW09 and store list of sed files 
  // containing *modelString along with it's COSANGLE value.
  // This SED.INFO file is different than open_SEDINFO_FILE().


  char 
    SEDINFO_FILE[MXPATHLEN]
    ,PARNAME_COSANGLE[] = "COSANGLE"
    ,PARNAME_DC[]       = "DC"
    ,PARNAME_IGNIT[]    = "IGNIT"
    ,PARNAME_DM15[]     = "DM15SED"
    ,PARNAME_MB[]       = "MBSED"
    ,parName[60]
    ,sedFile[MXPATHLEN]
    ,c_get[60]
    ,*name
    ,*modelString
    ,*ptrTmp
    ,fnam[] = "read_SEDINFO_KRW09" 
    ;

  int NPAR, ipar, indx, NCOS, NMODELS ;
  int IPAR_COSANGLE, IPAR_DC, IPAR_IGNIT, IPAR_DM15, IPAR_MB, USESED ;
  double parTmp[20], COSANGLE, COSANGLE_POS, COSANGLE_NEG ; 
  double DC, IGNIT, DC_LAST, IGNIT_LAST ;
  double DM15, MB ;

  FILE *fp ;

  // ------------ BEGIN ------------

  KRW09.NPAR = KRW09.N_COSANGLE = -9;
  KRW09.IPAR_COSANGLE = KRW09.IPAR_DC = KRW09.IPAR_IGNIT ;
  KRW09.IPAR_DM15 = KRW09.IPAR_MB = 0;
  KRW09.FLUX_ERRFLAG = 0 ;

  indx = SMEARMODEL_ID.KRW09 ;
  name = SMEARMODEL_DEF[indx].NAME ;
  modelString = INPUTS.SMEARMODEL_PARSTRING[1] ;


  sprintf(SEDINFO_FILE,"%s/SED.INFO", KRW09.DIR );
  if ( ( fp = fopen(SEDINFO_FILE, "rt") ) == NULL ) {
    sprintf(c1err, "Cannot open SED-INFO file:");
    sprintf(c2err,"'%s'", SEDINFO_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  NPAR = IPAR_COSANGLE = -9 ;
  IPAR_DM15 = IPAR_MB = -9 ;
  NCOS = 0;
  NMODELS  = 0;

  KRW09.COSANGLE_MIN = +999. ;
  KRW09.COSANGLE_MAX = -999. ;
  KRW09.ICOS_NEAR90DEG[0] = -9;
  KRW09.ICOS_NEAR90DEG[1] = -9;

  COSANGLE_POS =  999.0 ;
  COSANGLE_NEG = -999.0 ;

  DC_LAST = IGNIT_LAST = -9.0 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"NPAR:") == 0 ) 
      { readint(fp, 1, &NPAR); }

    if ( strcmp(c_get,"PARNAMES:") == 0 ) {
      for ( ipar=1; ipar <= NPAR; ipar++ ) {

	readchar(fp, parName);

	if ( strcmp(parName, PARNAME_COSANGLE) == 0 ) 
	  { IPAR_COSANGLE = ipar ; }	  
	if ( strcmp(parName, PARNAME_DC) == 0 ) 
	  { IPAR_DC = ipar ; }	  
	if ( strcmp(parName, PARNAME_IGNIT) == 0 ) 
	  { IPAR_IGNIT = ipar ; }	  
	if ( strcmp(parName, PARNAME_DM15) == 0 ) 
	  { IPAR_DM15 = ipar ; }	  
	if ( strcmp(parName, PARNAME_MB) == 0 ) 
	  { IPAR_MB = ipar ; }	  
      }
    }

    if ( strcmp(c_get,"FLUX_ERRFLAG:") == 0 ) 
      { readint(fp, 1, &KRW09.FLUX_ERRFLAG); }

    // Now store SEDs that contain the string given by
    // the first user-input parameter.
    
    if ( strcmp(c_get,"SED:") == 0 ) {
      readchar(fp, sedFile);

      USESED = 0;
      ptrTmp = strstr(sedFile,modelString);

      if ( ptrTmp != NULL    ) { USESED = 1; } // select from KRW09
      if ( KRW09.SUPER8FLAG  ) { USESED = 1; } // use all in Blondin's 8

      if ( USESED ) {
	readdouble(fp, NPAR, &parTmp[1] );
      
	COSANGLE = parTmp[IPAR_COSANGLE];
	DC       = parTmp[IPAR_DC];
	IGNIT    = parTmp[IPAR_IGNIT];
	DM15     = parTmp[IPAR_DM15];
	MB       = parTmp[IPAR_MB];

	if ( DC != DC_LAST || IGNIT != IGNIT_LAST ) { 
	  NMODELS++ ; 
	  NCOS = 0 ;
 	  if ( NMODELS >= MXMODEL_KRW09 ) {
	    sprintf(c1err,"NMODELS=%d exceeds bound.", NMODELS);
	    sprintf(c2err,"See  #definne MXMODEL_KRW09");
	    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
	  }
	  KRW09.DC[NMODELS]    = DC ;
	  KRW09.IGNIT[NMODELS] = IGNIT ;
	}
	DC_LAST    = DC ;
	IGNIT_LAST = IGNIT ;
	NCOS++ ;
	
	sprintf(KRW09.SEDFILE_LIST[NMODELS][NCOS],"%s", sedFile);
	KRW09.DM15_LIST[NMODELS][NCOS] = DM15 ;
	KRW09.MB_LIST[NMODELS][NCOS]   = MB ;

	if ( NMODELS == 1 ) {

	  KRW09.COSANGLE_LIST[NCOS] = COSANGLE;
	  
	  if ( COSANGLE < KRW09.COSANGLE_MIN ) 
	    { KRW09.COSANGLE_MIN = COSANGLE; }
	  
	  if ( COSANGLE > KRW09.COSANGLE_MAX ) 
	    { KRW09.COSANGLE_MAX = COSANGLE; }

	  if ( COSANGLE >= 0.0 && COSANGLE < COSANGLE_POS ) {
	    KRW09.ICOS_NEAR90DEG[1] = NCOS ;
	    COSANGLE_POS = COSANGLE ;
	  }
	
	  if ( COSANGLE <= 0.0 && COSANGLE > COSANGLE_NEG ) {
	    KRW09.ICOS_NEAR90DEG[0] = NCOS ;
	    COSANGLE_NEG = COSANGLE ;
	  }
	}
	/*
	printf("\t Select sedFile = %s  :  COSANGLE = %f \n",
	       sedFile, COSANGLE );
	*/
      }

    }  // end of SED: key

  }  // end of fscanf


  KRW09.NPAR          = NPAR ;
  KRW09.N_COSANGLE    = NCOS;
  KRW09.IPAR_COSANGLE = IPAR_COSANGLE ;
  KRW09.IPAR_IGNIT    = IPAR_IGNIT ;
  KRW09.IPAR_DC       = IPAR_DC ;
  KRW09.IPAR_DM15     = IPAR_DM15 ;
  KRW09.IPAR_MB       = IPAR_MB ;
  KRW09.NMODELS       = NMODELS ;

  fclose(fp);

  printf("    %s : NPAR=%d  IPAR_COSANGLE=%d   N_COSANGLE=%d \n",
	 name, NPAR, IPAR_COSANGLE, NCOS);
  printf("    %s : ICOS(90deg)= %d - %d \n",
	 name, KRW09.ICOS_NEAR90DEG[0], KRW09.ICOS_NEAR90DEG[1] );
  printf("    %s : Found %d models defined by DC,IGNIT \n", 
	 name, NMODELS) ;


  sprintf(c2err,"Check SED.INFO file in %s", name);
  if ( IPAR_COSANGLE <= 0 ) {
    sprintf(c1err,"Could not find required PARNAME: '%s'", 
	    PARNAME_COSANGLE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( IPAR_DC <= 0 ) {
    sprintf(c1err,"Could not find required PARNAME: '%s'", 
	    PARNAME_DC);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( IPAR_IGNIT <= 0 ) {
    sprintf(c1err,"Could not find required PARNAME: '%s'", 
	    PARNAME_IGNIT );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( IPAR_DM15 <= 0 ) {
    sprintf(c1err,"Could not find required PARNAME: '%s'", 
	    PARNAME_DM15 );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( IPAR_MB <= 0 ) {
    sprintf(c1err,"Could not find required PARNAME: '%s'", 
	    PARNAME_MB );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  //  debugexit(fnam); // DDDDDDDDDDDDDDD

} // end of read_SEDINFO_KRW09



// *****************************************
void malloc_KRW09(int imodel) {

  long int MEMTOT ;
  int MEM, NCOS, i, MXBIN_FLUX ; 
  double DMEMTOT ;
  char ctmp[60];
  char fnam[] = "malloc_KRW09" ;

  // ------------ BEGIN ------------

  MXBIN_FLUX = MXBIN_SED_SEDMODEL/5  ;
  KRW09.MXBIN_FLUX = MXBIN_FLUX ;

  NCOS = KRW09.N_COSANGLE ;

  MEM = 0 ;  MEMTOT=0;

  MEM       = (NCOS+1) * sizeof(double *);
  KRW09.DAY[imodel] = (double**)malloc(MEM);
  KRW09.LAM[imodel] = (double**)malloc(MEM);
  MEMTOT   += (2*MEM) ;
    
  MEM = MXBIN_DAYSED_SEDMODEL * sizeof(double); 
  for ( i=0; i <= NCOS; i++ ) 
    { KRW09.DAY[imodel][i] = (double*)malloc(MEM);    MEMTOT += MEM ; }
    
  MEM = MXBIN_LAMSED_SEDMODEL * sizeof(double); 
  for ( i=0; i <= NCOS; i++ ) 
    { KRW09.LAM[imodel][i] = (double*)malloc(MEM);    MEMTOT += MEM ; }
    

  // Now allocate for FLUX and FLUXERR arrays.
  MEM           = (NCOS+1) * sizeof(double*);
  KRW09.FLUX[imodel]    = (double**)malloc(MEM);
  KRW09.FLUXERR[imodel] = (double**)malloc(MEM);
    
  MEM = MXBIN_FLUX * sizeof(double); 
  for ( i=0; i <= NCOS; i++ )  { 
    KRW09.FLUX[imodel][i]    = (double*)malloc(MEM) ;  MEMTOT += MEM ; 
    KRW09.FLUXERR[imodel][i] = (double*)malloc(MEM) ;  MEMTOT += MEM ; 
  }

  DMEMTOT = (double)MEMTOT ;
  printf("\n");
  
  if ( imodel == 0 ) 
    { sprintf(ctmp,"AVG"); }
  else
    { sprintf(ctmp,"DC=%d IGNIT=%d"
	      ,(int)KRW09.DC[imodel] 
	      ,(int)KRW09.IGNIT[imodel]  );
    }

  int indx;
  char *name ;
  indx = SMEARMODEL_ID.KRW09 ;
  name = SMEARMODEL_DEF[indx].NAME ;
  printf("  Allocated %6.2f MB for %s SEDs  (model %d: %s)\n", 
	 DMEMTOT/1.0E6, name, imodel, ctmp);
  fflush(stdout);

} // end of malloc_KRW09


// *****************************************
void  read_SED_KRW09(int imodel, int icos) {

  int  i, NCOS, indx, NBIN_FLUX, NDAY, NLAM ;

  char 
    *name
    ,sedComment[60]
    ,sedFile[MXPATHLEN]
    ,fnam[] = "read_SED_KRW09" 
    ;

  // --------- BEGIN ---------------

  indx = SMEARMODEL_ID.KRW09 ;
  name = SMEARMODEL_DEF[indx].NAME ;

  NCOS = KRW09.N_COSANGLE ;

  // on first SED allocate memory for DAY and LAM arrays.
  

  if ( icos == 1 ) { malloc_KRW09(imodel);  }

  // Now read each SED and store it


  sprintf(sedComment,"(icos=%d/%d)", icos, NCOS);
  sprintf(sedFile,"%s/%s", KRW09.DIR, KRW09.SEDFILE_LIST[imodel][icos] );


  KRW09.NDAY[imodel][icos] = -9;
  KRW09.NLAM[imodel][icos] = -9;

  rd_sedFlux(sedFile, sedComment
	     ,KRW09.SEDRANGE_DAY, KRW09.SEDRANGE_LAM
	     ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL
	     ,KRW09.FLUX_ERRFLAG
	     ,&KRW09.NDAY[imodel][icos]
	     ,KRW09.DAY[imodel][icos]
	     ,&KRW09.DAYSTEP[imodel][icos]
	     ,&KRW09.NLAM[imodel][icos]
	     ,KRW09.LAM[imodel][icos]
	     ,&KRW09.LAMSTEP[imodel][icos]
	     ,KRW09.FLUX[imodel][icos]
	     ,KRW09.FLUXERR[imodel][icos] 
	     ) ; 
  

  NBIN_FLUX = KRW09.NDAY[imodel][icos] * KRW09.NLAM[imodel][icos] ;
  if ( NBIN_FLUX >= KRW09.MXBIN_FLUX ) {
    sprintf(c1err,"NDAY(%d) x NLAM(%d) = %d exceeds bound.",
	    KRW09.NDAY[imodel][icos], KRW09.NLAM[imodel][icos], NBIN_FLUX);
    sprintf(c2err,"Bound is MXBIN_FLUX = %d", KRW09.MXBIN_FLUX);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  
  NDAY = KRW09.NDAY[imodel][icos] ;
  NLAM = KRW09.NLAM[imodel][icos] ;
  printf("\t %d DAY bins (%6.2f to %6.2f)  %d LAM bins (%6.0f to %6.0f)\n"
	 ,NDAY
	 ,KRW09.DAY[imodel][icos][0]
	 ,KRW09.DAY[imodel][icos][NDAY-1]
	 ,NLAM
	 ,KRW09.LAM[imodel][icos][0]
	 ,KRW09.LAM[imodel][icos][NLAM-1]
	 );
  fflush(stdout); 
  
    
} // end of read_SED_KRW09


// *****************************************
void  smooth_SED_KRW09(int imodel, int icos) {

  int 
    iday, ilam, jflux 
    ,NDAY, NLAM, OPT_SMOOTH, LDMP
    ;

  double 
    flux, fluxerr, LAMSMOOTH
    ,TEMPFLUX_IN[MXLAM_SMOOTH]
    ,TEMPFLUXERR_IN[MXLAM_SMOOTH]
    ;

  char
    *name
    ,fnam[] = "smooth_SED_KRW09" 
    ;

  // ------------ BEGIN --------------

  NDAY = KRW09.NDAY[imodel][icos] ;
  NLAM = KRW09.NLAM[imodel][icos] ;

  OPT_SMOOTH = OPT_SMOOTH_GAUSS ;
  LAMSMOOTH  = INPUTS.SMEARMODEL_PARVAL[3] ; // smoothing, A

  for ( iday = 0; iday < NDAY; iday++ ) {

    for ( ilam = 0; ilam < NLAM; ilam++ ) {
      jflux   = NLAM*iday + ilam ;
      flux    = KRW09.FLUX[imodel][icos][jflux];
      fluxerr = KRW09.FLUXERR[imodel][icos][jflux];
      TEMPFLUX_IN[ilam]    = flux;
      TEMPFLUXERR_IN[ilam] = fluxerr;
    } // end of ilam

    
    sedSmooth( OPT_SMOOTH, LAMSMOOTH
	      ,NLAM, KRW09.LAM[imodel][icos]
	      ,TEMPFLUX_IN, TEMPFLUXERR_IN
	      ,TEMP_SEDSMOOTH.FLUX[iday]     // output
	      ,TEMP_SEDSMOOTH.FLUXERR[iday]  // output
	      );

    // copy the smoothed SED back into the original array;
    // i.e, overwrite original [unsmoothed] KRW09 SED
    for ( ilam = 0; ilam < NLAM; ilam++ ) {
      jflux   = NLAM*iday + ilam ;

      LDMP = ( iday == -4 && icos == 12 && ilam >= 300  && ilam < 500 ) ; 
      if ( LDMP  ) {
	printf(" 222222  %7.2f  %le  %le \n"
	       , KRW09.LAM[imodel][icos][ilam]    // lambda (A)
	       , KRW09.FLUX[imodel][icos][jflux]  // unsmoothed flux
	       , TEMP_SEDSMOOTH.FLUX[iday][ilam] ); // smoothed flux
	fflush(stdout);
      }

      KRW09.FLUX[imodel][icos][jflux]    = TEMP_SEDSMOOTH.FLUX[iday][ilam];
      KRW09.FLUXERR[imodel][icos][jflux] = TEMP_SEDSMOOTH.FLUXERR[iday][ilam];
    } // end of ilam
    

  } // iday

}  // end of smooth_SED_KRW09


// **************************************************
void  fluxavg_KRW09(void) {

  // Created Jun 23, 2012 by R.Kessler
  // compute grand average flux vs. DAY and LAM.
  // Beware that each SED can have different binning
  // for DAY and LAM, so define specific bins for FLUXAVG.

  int 
    jflux_avg
    ,NCOS, NMODELS
    ,icos, imodel
    ,iday_avg, ilam_avg, NDAY_AVG, NLAM_AVG
    ;

  int IMODEL_AVG = 0;
  double FLUX, DAYSTEP, LAMSTEP ;
  double LAM_AVG, DAY_AVG, DAY, LAM;

  char fnam[] = "fluxavg_KRW09" ;

  // -------------- BEGIN ---------------

  NMODELS = KRW09.NMODELS; 
  NCOS    = KRW09.N_COSANGLE ;

  printf("\t Get Spectral flux averaged over %d MODELS and %d COS bins\n",
	 NMODELS, NCOS);

  // construct IMODEL=0 that is the average of all models.
  // The min/max values of the AVG model are chosen to be 
  // valid for all of the imodel>0 models.
  // Use 1 day bins and 20 A bins for AVG model

  malloc_KRW09(IMODEL_AVG);
  KRW09.DAYSTEP[IMODEL_AVG][0] =  1.0 ;  
  KRW09.LAMSTEP[IMODEL_AVG][0] = 20.0 ;  
  

  DAYSTEP = KRW09.DAYSTEP[IMODEL_AVG][0] ;
  LAMSTEP = KRW09.LAMSTEP[IMODEL_AVG][0] ;

  NDAY_AVG = NLAM_AVG = 0 ;

  DAY = 1.0 ;
  for ( LAM=1000.0; LAM < 12000.0; LAM+=LAMSTEP ) {

    if ( ALLMODELS_VALID(DAY,LAM) == 1 ) { 
      NLAM_AVG++ ;
      KRW09.NLAM[IMODEL_AVG][0] = NLAM_AVG;
      KRW09.LAM[IMODEL_AVG][0][NLAM_AVG-1] = LAM ;
    }
  }

  LAM = 5000.0 ;
  for ( DAY = -20.9; DAY < 25.1; DAY+= DAYSTEP ) {
    if ( ALLMODELS_VALID(DAY,LAM) == 1 ) { 
      NDAY_AVG++ ;
      KRW09.NDAY[IMODEL_AVG][0] = NDAY_AVG;
      KRW09.DAY[IMODEL_AVG][0][NDAY_AVG-1] = DAY ;
    }
  }

  if ( NDAY_AVG < 2 || NLAM_AVG < 2 ) {
    sprintf(c1err,"Invalid: NDAY=%d and NLAM=%d", NDAY_AVG, NLAM_AVG);
    sprintf(c2err,"for AVGMODEL");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\t AVGMODEL DAY-RANGE: %6.1f to %6.1f  (%d bins)\n"
	 ,KRW09.DAY[IMODEL_AVG][0][0]
	 ,KRW09.DAY[IMODEL_AVG][0][NDAY_AVG-1], NDAY_AVG );
  printf("\t AVGMODEL LAM-RANGE: %6.0f to %6.0f  (%d bins)\n"
	 ,KRW09.LAM[IMODEL_AVG][0][0]
	 ,KRW09.LAM[IMODEL_AVG][0][NLAM_AVG-1], NLAM_AVG );
  fflush(stdout);



  for ( iday_avg = 0; iday_avg < NDAY_AVG; iday_avg++ ) {
    for ( ilam_avg = 0; ilam_avg < NLAM_AVG; ilam_avg++ ) {

      jflux_avg   = NLAM_AVG*iday_avg + ilam_avg ;

      KRW09.FLUX[IMODEL_AVG][0][jflux_avg] = 0.0 ;

      LAM_AVG = KRW09.LAM[IMODEL_AVG][0][ilam_avg] ;
      DAY_AVG = KRW09.DAY[IMODEL_AVG][0][iday_avg] ;

      // now interpolate each model and increment AVG model
      for ( imodel=1; imodel <= NMODELS; imodel++ ) {
	for ( icos=1; icos <= NCOS; icos++ ) {
	  FLUX = interp_SEDFLUX_KRW09(imodel, icos,  LAM_AVG, DAY_AVG); 
	  KRW09.FLUX[IMODEL_AVG][0][jflux_avg] += FLUX ;

	} // icos
      } // imodel

    } // ilam_avg
  } // day_avg


  // divide sum by NMODELS*NCOS to get grand average
  double XN = (double)(NCOS*NMODELS);
  for ( iday_avg = 0; iday_avg < NDAY_AVG; iday_avg++ ) {
    for ( ilam_avg = 0; ilam_avg < NLAM_AVG; ilam_avg++ ) {
      jflux_avg   = NLAM_AVG*iday_avg + ilam_avg ;
      KRW09.FLUX[IMODEL_AVG][0][jflux_avg] /= XN ;
    } // ilam
  } // day


} // end of   fluxavg_KRW09

// ******************************************
int ALLMODELS_VALID(double DAY, double LAM) {

  // Returns TRUE (1) if the input DAY and LAM
  // are valid for all models and COSANGLEs.

  int imodel, icos, NDAY, NLAM, LDMP,  DAYMAX ;
  char fnam[] = "ALLMODELS_VALID" ;

  LDMP = ( fabs(DAY-21.5) < 0.001 ) ;

  for(imodel=1; imodel <= KRW09.NMODELS; imodel++ ) {
    for ( icos=1; icos <= KRW09.N_COSANGLE; icos++ ) {
      NDAY = KRW09.NDAY[imodel][icos];
      NLAM = KRW09.NLAM[imodel][icos];

      /*
      DAYMAX = KRW09.DAY[imodel][icos][NDAY-1] ;
      if ( LDMP && DAYMAX < 23.0 ) {
	
	printf(" xxxx day23: imodel=%d icos=%2d  DAY(min,max)=%7.2f,%7.2f\n"
	       ,imodel, icos
	       ,KRW09.DAY[imodel][icos][0]
	       ,KRW09.DAY[imodel][icos][NDAY-1]   ) ;
	fflush(stdout);
      }      */

      if ( DAY < KRW09.DAY[imodel][icos][0]      ) { return 0; }
      if ( DAY > KRW09.DAY[imodel][icos][NDAY-1] ) { return 0; }
      if ( LAM < KRW09.LAM[imodel][icos][0]      ) { return 0; }
      if ( LAM > KRW09.LAM[imodel][icos][NLAM-1] ) { return 0; }
    }
  }

  //  if we get here, return true.
  return 1; 

} // end of ALLMODELS_VALID

// **************************************************
double shapeFudge(double lambda, double shapeVal ) {

  // Return shape fudge for this wavelength "lambda"
  // and this "shapeVal" (x1, Delta  ...)
  // based on the INDEX_SHFUDGEFUN value selected by user

  double arg, alpha, x0;
  double unity, c0, c1, a, b, x, sq ;
  char fnam[] = "shapeFudge";
  double *ptrpar ;

  // ----------------- BEGIN -----------

    x0 = 1.0 ;
    if ( INPUTS.NBIN_SHAPEPAR <=  1 && IPAR_S2x1 <= 0 ) {
      return x0 ;
    }

    if ( INDEX_SHFUDGEFUN == INDEX_SHFUDGE_DEFAULT ) {

      alpha = INPUTS.SHAPEFUN_PARVAL[1] ; 
      arg   = 0.4 * (alpha * shapeVal ) ;
      x0    = pow(TEN,arg);
      
    } else if ( INDEX_SHFUDGEFUN == INDEX_SHFUDGE_SAKO08 ) {

      // see equation C2 in Sako et al. 2008

      ptrpar = &INPUTS.SHFUDGEFUN_PARVAL[1] ;

      unity = *(ptrpar+0);
      c0    = *(ptrpar+1);
      c1    = *(ptrpar+2);
      b     = *(ptrpar+3);

      x   = shapeVal - unity;
      sq  = x*x;
      a   = c0 - c1*lambda ;
      arg = -0.4 * (a*x + b*sq);
      x0  = pow(10, arg);
    }

  return x0 ;

} // end of shapeFudge


// *********************************************
double stretchFlux(double str, int iday, int ilam ) {

  // return Flux at stretched value  F[(t-t0)/stretch]

  int    iday2, idaystr, NDAY, LDMP ;
  double 
     F, F0, F1, Fstr, logF, logF0, logF1
    ,day, day0, day2, day3, daystr, dayfrac, tmp
    ,DAYSTEP, DAYMIN, DAYMAX
    ,lam, Tpeak
    ; 

  // ------------- BEGIN --------------

  day    = TEMP_SEDMODEL.DAY[iday];
  lam    = TEMP_SEDMODEL.LAM[ilam];
  Tpeak  = TPEAK_LAMPOLY(lam);


  // get stretched day relative to Tpeak
  daystr = (day-Tpeak)/str + Tpeak ;

  // will have to interpolate since daystr does not fall on a grid point

  DAYSTEP  = TEMP_SEDMODEL.DAYSTEP ;
  DAYMIN   = TEMP_SEDMODEL.MINDAY ;
  DAYMAX   = TEMP_SEDMODEL.MAXDAY - TEMP_SEDMODEL.DAYSTEP - 0.001 ;
  NDAY     = TEMP_SEDMODEL.NDAY ;

  if ( daystr < DAYMIN ) {
    idaystr = 1 ;
  }
  else if ( daystr >= DAYMAX ) {
    idaystr = NDAY - 2 ;
  }
  else {

    for ( iday2 = 0; iday2 < NDAY-1 ; iday2++ ) {
      day2 = TEMP_SEDMODEL.DAY[iday2];
      day3 = TEMP_SEDMODEL.DAY[iday2+1];
      if ( daystr >= day2 && daystr < day3 ) {
	idaystr = iday2;
      }
    }

  }

  day0 = TEMP_SEDMODEL.DAY[idaystr];

  // interpolate in log space to avoid negative fluxes.
  F0 = TEMP_SEDSMOOTH.FLUX[idaystr][ilam] ;
  F1 = TEMP_SEDSMOOTH.FLUX[idaystr+1][ilam] ;
  logF0   = log10(F0);
  logF1   = log10(F1);
  dayfrac = (daystr-day0)/DAYSTEP ;
  logF    = logF0 + (logF1-logF0) * dayfrac ; 

  Fstr  = pow(10.0,logF);

  LDMP = ( iday == -998 && ilam == 101 && str > 1.1 );
  if ( LDMP ) {
    printf("\n");
    printf(" xxx DUMP FLUX INTERP inside stretchFlux \n");
    printf(" xxx DAY[MIN,MAX,STEP] = %5.1f , %5.1f , %5.1f \n",
	   DAYMIN, DAYMAX, DAYSTEP );
    printf(" xxx str = %6.3f  iday=%d  ilam=%d \n", str, iday, ilam );
    printf(" xxx day=%5.1f -> daystr = %6.2f (day0=%5.1f, dayfrac=%5.3f) \n", 
	   day, daystr, day0, dayfrac );
    printf(" xxx F0=%9.3le  F1=%9.3le  Fstr = %9.3le \n", F0, F1, Fstr );
    debugexit("flux interp");
  }

  /*
  // dummy correction to test code xxxx
  tmp  = 1.0 + (str-1.0) * fabs(day/10.0) ;
  Fstr = F * tmp ;
  */

  return Fstr ;
}


// *******************************
double TPEAK_LAMPOLY(double lam) {
  double LAMDIF, SQLAMDIF, Tpeak;

  LAMDIF    = lam - INPUTS.LAMPOLY_TPEAK[2];
  SQLAMDIF  = LAMDIF * LAMDIF;
  Tpeak     = INPUTS.LAMPOLY_TPEAK[0]
            + INPUTS.LAMPOLY_TPEAK[1]*SQLAMDIF;

  return Tpeak ;
} //  end of TPEAK_LAMPOLY


// *******************************
double FLUXFUDGE_LAMPOLY(double lam) {
  // return fduge multiplier from lambda polynomial
  int i;
  double sum, tmp, xi ; 

  sum = 0.0 ;
  for ( i = 0; i <= 3; i++ ) {
    xi  = (double)i;
    tmp = INPUTS.LAMPOLY_FLUX[i] * pow(lam,xi);
    sum += tmp ;

  }
  return sum ;

} //  end of FLUXFUDGE_LAMPOLY

// **************************************************
double colorFudge(double lambda, double colorVal) {

  // Return color fudge for this wavelength "lambda"
  // and this "colorVal" (SALT2c, AV  ...)

  double *ptrpar, XTMAG, arg, x0, beta, colorlaw ;
  char fnam[] = "colorFudge";

  // ---------- BEGIN -------------

  ptrpar = &INPUTS.COLORFUN_PARVAL[1] ;


  if ( INDEX_COLORFUN == INDEX_SALT2colorlaw0 ) {
    colorlaw = SALT2colorlaw0 ( lambda, colorVal, ptrpar );

    beta = *(ptrpar+4); //   5'th element is SALT2 beta
    arg  = -0.4 * (beta * colorVal ) ;
    x0   = pow(TEN,arg);
    return  (x0 * colorlaw);

    // fluxFudge = colorlaw * 10**[-0.4*beta*c]

  } 
  else if ( INDEX_COLORFUN ==  INDEX_SALT2colorlaw1 ) {

    colorlaw = SALT2colorlaw1 ( lambda, colorVal, ptrpar );

    beta = ptrpar[9] ; //   10'th element is SALT2 beta
    arg  = -0.4 * (beta * colorVal ) ;
    x0   = pow(TEN,arg);
    return  (x0 * colorlaw);

  }
  else if ( INDEX_COLORFUN ==  INDEX_CCM89 ) {
    XTMAG = GALextinct( *ptrpar, colorVal, lambda, 89 ) ;
    arg   = -0.4*XTMAG ;
    return pow(TEN,arg);
  }
  else if ( INDEX_COLORFUN ==  INDEX_CCM94 ) {
    XTMAG = GALextinct( *ptrpar, colorVal, lambda, 94 );
    arg   = -0.4*XTMAG ;
    return pow(TEN,arg);
  }
  else {
    sprintf(c1err,"Invalid INDEX_COLORFUN = %d", INDEX_COLORFUN);
    errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
  }


} // end of colorFudge



// ************************************
double Btransfun(double lam) {

  // Aug 25 2013: moved to genmag_SEDtools.c
  //
  // return Bessell B transmission at this lambda.
  // Transmission is energy fraction, not count fraction.
  // just hard-wire the B-transmission vs. lambda.
  // to avoid more I/O

#define NBIN_Btrans 21
  double Btrans_LAMBDA[NBIN_Btrans] = { 3600., 3700., 3800., 3900., 4000., 4100., 4200., 4300., 4400., 4500., 4600., 4700., 4800., 4900., 5000., 5100., 5200., 5300., 5400., 5500., 5600. } ;

  double Btrans_VALUE[NBIN_Btrans]  = { 0.0,   0.03,  0.134, 0.567, 0.920, 0.978, 1.00, 0.978, 0.935, 0.853, 0.740, 0.640, 0.536, 0.424, 0.325, 0.235, 0.150, 0.095, 0.043, 0.009, 0.00 } ;


  double trans, lam0, lam1, frac, tr0,tr1 ;
  int ilam;

  char fnam[] = "Btransfun";

  //  ----------- BEGIN -----------

  trans = 0.0 ;

  if ( lam < Btrans_LAMBDA[1] )             return trans ;
  if ( lam > Btrans_LAMBDA[NBIN_Btrans-1] ) return trans ;

  for ( ilam=0; ilam < NBIN_Btrans-1 ;  ilam++ ) {
    lam0 = Btrans_LAMBDA[ilam] ;
    lam1 = Btrans_LAMBDA[ilam+1] ;

    tr0  = Btrans_VALUE[ilam];
    tr1  = Btrans_VALUE[ilam+1];

    if ( lam >= lam0 && lam <= lam1 ) {
      frac = (lam - lam0)/(lam1-lam0);
      trans = tr0 + frac*(tr1-tr0);
      return trans;
    }

  } // ilam

  sprintf(c1err,"Could not evaluate Btrans at lam=%f", lam);
  errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 

} // end of Btransfun



// ********************************
void sedSmooth_driver(int ised) {

  // Sep 7, 2010
  // start with TEMP_SEDMODEL array and create TEMP_SEDSMOOTH array
  // This function is called after reading input SED, but before
  // applying the analytic fudge.

  int iday, ilam, jflux, NDAY, NLAM, NDUMP, i ;

  double 
    day, lam, flux, fluxerr, Trest
    ,smoothArgs[20]
    ,Trest_dump[MXDUMP_SMOOTH]
    ,TEMPFLUX_IN[MXLAM_SMOOTH]
    ,TEMPFLUXERR_IN[MXLAM_SMOOTH]
    ,TEMPFLUX_OUT[MXLAM_SMOOTH]
    ,TEMPFLUXERR_OUT[MXLAM_SMOOTH]
    ;

  char *ptrSmooth, *ptrSed;

  // ---------- BEGIN -----------
  
  NDAY = TEMP_SEDMODEL.NDAY ;
  NLAM = TEMP_SEDMODEL.NLAM ;

  // check if this SED name is scheduled for a dump

  ptrSed = SEDMODEL.FILENAME[ised] ;
  NDUMP  = 0;
  for  ( i = 1; i <= INPUTS.NDUMP_SMOOTH; i++ ) {
    ptrSmooth = INPUTS.DUMP_SMOOTH_SEDNAME[i] ;
    if ( strcmp(ptrSed,ptrSmooth) == 0 ) { 
      NDUMP++ ;
      Trest_dump[NDUMP] = INPUTS.DUMP_SMOOTH_TREST[i] ;
    }
  }


  for ( iday = 0; iday < NDAY; iday++ ) {

    // first load flux and  fluxerr arrays

    for ( ilam = 0; ilam < NLAM; ilam++ ) {
      day     = TEMP_SEDMODEL.DAY[iday];
      lam     = TEMP_SEDMODEL.LAM[ilam];
      jflux   = TEMP_SEDMODEL.NLAM*iday + ilam ;

      flux    = TEMP_SEDMODEL.FLUX[jflux];
      fluxerr = TEMP_SEDMODEL.FLUXERR[jflux];

      TEMPFLUX_IN[ilam]    = flux;
      TEMPFLUXERR_IN[ilam] = fluxerr;

      TEMP_SEDSMOOTH.FLUX[iday][ilam]     = TEMPFLUX_IN[ilam] ;
      TEMP_SEDSMOOTH.FLUXERR[iday][ilam]  = TEMPFLUXERR_IN[ilam] ;

    } // ilam


    // now check smoothing option
    // The default NOSMOOTH option has already been set,
    // so only check for non-trivial  options here.


    if ( INPUTS.OPT_SMOOTH == OPT_SMOOTH_NBRAVG ) {

      sedSmooth_nbravg(NLAM, TEMP_SEDMODEL.LAM,
		       TEMPFLUX_IN, TEMPFLUXERR_IN,
		       TEMP_SEDSMOOTH.FLUX[iday],    // output
		       TEMP_SEDSMOOTH.FLUXERR[iday]  // output
		       );
    }
    else if ( INPUTS.OPT_SMOOTH == OPT_SMOOTH_TOPHAT ||
	      INPUTS.OPT_SMOOTH == OPT_SMOOTH_GAUSS ) {

      smoothArgs[0] = INPUTS.LAMSCALE_SMOOTH ;
      sedSmooth(INPUTS.OPT_SMOOTH, smoothArgs[0],
		NLAM, TEMP_SEDMODEL.LAM,
		TEMPFLUX_IN, TEMPFLUXERR_IN,
		TEMP_SEDSMOOTH.FLUX[iday],    // output
		TEMP_SEDSMOOTH.FLUXERR[iday]  // output
		);

    }

    for ( i = 1; i <= NDUMP; i++ ) {
      Trest = Trest_dump[i] ;
      if ( fabs(day-Trest) < 0.5 ) {
	sedSmooth_dump(ised, Trest, NLAM, 
		       TEMP_SEDMODEL.LAM,
		       TEMPFLUX_IN, 
		       TEMPFLUXERR_IN,
		       TEMP_SEDSMOOTH.FLUX[iday],
		       TEMP_SEDSMOOTH.FLUXERR[iday]
		       );
      }

    } // NDUMP loop


  } // iday

} // end of sedSmooth_driver



// **********************************
void sedSmooth_nbravg(int NLAM, double *lambda,
		      double *FLUX_IN, double *FLUXERR_IN,
		      double *FLUX_SMOOTH, double *FLUXERR_SMOOTH ) {

  // Sep 2010, R.Kessler
  // very trivial algorithm of just averaging each bin with
  // its four neighboring bins.

  int ilam;
  double flux, fluxerr, f1, f2, f3, f4, f5;

  // ------------ BEGIN -------------

  for ( ilam = 0; ilam < NLAM; ilam++ ) {

    *(FLUX_SMOOTH+ilam)    = *(FLUX_IN+ilam);
    *(FLUXERR_SMOOTH+ilam) = *(FLUXERR_IN+ilam);

    if ( ilam < 2        ) continue ;
    if ( ilam > NLAM - 2 ) continue ;

    f1 = *(FLUX_IN+ilam-2) ;
    f2 = *(FLUX_IN+ilam-1) ;
    f3 = *(FLUX_IN+ilam+0) ;
    f4 = *(FLUX_IN+ilam+1) ;
    f5 = *(FLUX_IN+ilam+2) ;

    flux    = (f1 + f2 + f3 + f4 + f5)/5.0 ;
    fluxerr = *(FLUXERR_IN+ilam) ;

    *(FLUX_SMOOTH+ilam)    = flux ;
    *(FLUXERR_SMOOTH+ilam) = fluxerr ;

  } // ilam

} // end of sedSmooth_nbravg


// ***********************************
void sedSmooth_dump(int ised, double Trest, int NLAM, double *lam, 
		    double *flux_in,  double *fluxerr_in,
		    double *flux_smooth, double *fluxerr_smooth
		    ) {

  // dump with following format:
  // LAMBDA  FLUX_IN  FLUXERR_IN   FLUX_SMOOTH  FLUXERR_SMOOTH
  //

  int ilam ;
  char dumpFile[200];

  FILE *fp;

  // ----------- BEGIN -------------

  sprintf(dumpFile,"%s_day%d.DUMP",
	  SEDMODEL.FILENAME[ised], (int)Trest );

  printf(" DUMP_SMOOTH -> %s \n", dumpFile );

  fp = fopen(dumpFile, "wt");

  for ( ilam = 0; ilam < NLAM; ilam++ ) {

    fprintf(fp, "%7.1f  %9.3le %9.3le   %9.3le %9.3le \n",
	    *(lam         + ilam),
	    *(flux_in     + ilam),
	    *(fluxerr_in  + ilam),
	    *(flux_smooth    + ilam),
	    *(fluxerr_smooth + ilam)
	    );

  } // lam
    
  fclose(fp);

} // end of sedSmooth_dump

// ********************************
void sedWrite(FILE *fp_out) {

  // Sep 7, 2010
  // Write out TEMP_SEDFUDGE array to SED file.


  int iday, ilam;
  double day, lam, flux, fluxerr ;

  char ctmp[80], ctmperr[80];
  char debug1[300], debug2[200]; // FLAG_DEBUG_SIMSED

  // ---------- BEGIN -----------

  for ( iday = 0; iday < TEMP_SEDMODEL.NDAY; iday++ ) {

    for ( ilam = 0; ilam < TEMP_SEDMODEL.NLAM; ilam++ ) {
      day     = TEMP_SEDMODEL.DAY[iday];
      lam     = TEMP_SEDMODEL.LAM[ilam];

      flux    = TEMP_SEDFUDGE.FLUX[iday][ilam] ; 
      fluxerr = TEMP_SEDFUDGE.FLUXERR[iday][ilam] ; 
      

      if ( floor(day) == day ) 
	// write day as integer to save disk space
	{ sprintf(ctmp,"%.0f %.0f %.4le ", day, lam, flux ); }
      else
	{ sprintf(ctmp,"%7.3f %.0f %.4le ", day, lam, flux ); }
      
      if ( (SEDMODEL.OPTMASK & OPTMASK_FLUXERR_SEDMODEL) > 0 ) 
	{ sprintf(ctmperr,"%le", fluxerr );  }
      else
	{ sprintf(ctmperr,""); }

      if ( FLAG_DEBUG_SIMSED == 1 ) {
	  sprintf(debug1, "%le   %le   %le   %le   %le ", 
		  TEMP_SEDFUDGE.ORIG_FLUX[iday][ilam], 
		  TEMP_SEDFUDGE.POLYFUDGE[iday][ilam], 
		  TEMP_SEDFUDGE.CFUDGE[iday][ilam],
		  TEMP_SEDFUDGE.SFUDGE[iday][ilam],
		  TEMP_SEDFUDGE.SHAPECOR[iday][ilam]);


	} else {
	sprintf(debug1, ""); 
	sprintf(debug2, ""); 
      } // FLAG_DEBUG_SIMSED

      fprintf(fp_out,"%s %s %s \n", ctmp, ctmperr, debug1, debug2 );      
    }
  }

} // end of sedWrite


double lcstretch(double x1) {

  // returns light curve stretch as a function of the input shape parameter
  // function used depends on selected INDEX_STRETCHFUN

  double s, x, xsq, xcube ;
  double c0, c1, c2, c3, tau ;
  double *ptrpar ;

  // ---------- BEGIN -------------

  if ( INDEX_STRETCHFUN == INDEX_STRETCHFUN_NONE ) { return x1 ; }

  if ( INDEX_STRETCHFUN == INDEX_STRETCHFUN_DEFAULT ) {

    x    = x1 ;
    xsq  = x*x;  xcube = xsq*x;

    // From Guy 2007, convert x1 into stretch
    s = 0.98 + .091*x + .003*xsq - 7.5E-4*xcube ;
    // Nov 8, 2010: fudge so that s=1 at x1=0.
    s += 0.02 ; 

    // from fitting stretch model to SALT2 sim
    // s = 1.0 + 0.10258*x + 0.16449E-02*xsq ;

  } 
  else if ( INDEX_STRETCHFUN == INDEX_STRETCHFUN_SAKO08 ) {

    // see equation C4 of Sako et al. 2008

    double dm15 = x1 ; // RK - Aug 25 2013 
    ptrpar = &INPUTS.STRETCHFUN_PARVAL[1] ;

    c0 = *(ptrpar+0);
    c1 = *(ptrpar+1);
    c2 = *(ptrpar+2);
    c3 = *(ptrpar+3);

    x    = dm15 ;
    xsq  = x*x;  xcube = xsq*x;

    tau = c0 + c1*x + c2*xsq + c3*xcube ;
    s = (15.0/tau);

    fflush(stdout);

  }
  return s;
}  // end of lcstretch



// *******************************
void make_param_ntuple(void) {

  char cmd[600] ;
  char mvcmd[200];
  char ntupFile[] = "SED_PARAM.TUP" ;
  // ------------ BEGIN -------------

  sprintf(mvcmd,"mv combine_fitres.his %s", ntupFile );

  sprintf(cmd,"cd %s ; combine_fitres.exe %s > DMP ; %s ; rm DMP",
	  SIMSED_PATHMODEL_OUTPUT, INFO_FILENAME2, mvcmd);

  printf("\n Create SED-param ntuple : %s \n", ntupFile );
  fflush(stdout);
  int isys = system(cmd);
  
} // end of make_param_ntuple


// ===========================================================
//
// Rahul's smoothing code
//
// ===========================================================


/************************************
Created Sep 2010 by Rahul Biswass
Installed by RK

sedSmooth produces a smoothed spectrum and corresponding errors 
by smoothing with a filter. 

  INPUT:
  opt 	   =1 uses a tophat filter, opt=2 uses a Gaussian filter,
            any other opt echoes the usage. 

  smoothpar  = pointer to smoothing parameters:
             +0 => approx smoothing width in Angstroms

  Nlam 	     = Number of lambda bins in the original spectrum
  lambda     = array (0,Nlam-1) wavelength in each bin
  flux_in    = array (0,Nlam-1)  values of flux	
  fluxerr_in = array (0,Nlam-1)  values of flux	

  OUTPUT:
  flux_smooth    = array (0,Nlam -1) values of smoothed flux
  fluxerr_smooth = array (0,Nlam -1) values of errors on the smoothed flux



  Installation modifications (RK):
  - compute BINSIZE_LAMBDA instead of assuming 10 A.
  - use global parameters OPT_SMOOTH_TOPHAT[GAUSS] 
  - replace 'len' with NBIN_SMOOTH
  - put smoothpar as 2nd argument instead of last argument
    (so that all input args are before output args)

******************************/

int sedSmooth(int opt, double smoothpar, int Nlam, double *lambda,
	      double *flux_in, double *fluxerr_in,
	      double *flux_smooth, double *fluxerr_smooth	      
	      ) {

  double *filter , *fluxbox, *fluxboxerr;
  int i,j, NBIN_SMOOTH ;
  double   widthTobin;

  // -------------- BEGIN --------------

  NCALL_SMOOTH++ ;

  // compute binsize
  BINSIZE_LAMBDA  = *(lambda+1) - *(lambda+0);

  widthTobin = floor(smoothpar/BINSIZE_LAMBDA);
  NBIN_SMOOTH = (int)(2.*widthTobin + 1.); //made to be odd

  //printf("The number of bins smoothed over is %d\n" ,len);
  //Check that averaging length < total length of array

  if (NBIN_SMOOTH > Nlam) {
    printf("Should average over lengths less than total length");
    exit(EXIT_FAILURE);
  }	


filter = setfilter(NBIN_SMOOTH-1,opt);

for (j=0;j<Nlam;j++){
	fluxbox = setbox(NBIN_SMOOTH,j,Nlam,flux_in);
	flux_smooth[j] =0.;
	for (i=0;i<NBIN_SMOOTH;i++){
		flux_smooth[j]+=fluxbox[i]*filter[i];
	}
}
for (j=0;j<Nlam;j++){
	fluxboxerr =setbox(NBIN_SMOOTH,j,Nlam,fluxerr_in);
	fluxerr_smooth[j] =0.;
	for (i=0;i<NBIN_SMOOTH;i++){
		fluxerr_smooth[j]+=\
                         fluxboxerr[i]*fluxboxerr[i]*filter[i]*filter[i];
         }
         fluxerr_smooth[j] = sqrt(fluxerr_smooth[j]);
}
return 0;
}
/*
setfilter sets the filter of size determined by l and the option opt)
	opt =1 : tophat
	opt =2 : Gaussian
	opt = anything else: print usage and return
*/
double *setfilter(int l, int opt){
	static double filt[MXLAM_SMOOTH];
	double sum,t1,t2,width;
	int i;
	int k; 

	k = l/2; 
	width = k*BINSIZE_LAMBDA ;
	
	sum =0.;
	if (opt == OPT_SMOOTH_TOPHAT ){
	  for (i=0;i<l;i++){
	    filt[i] =1.0;
	    sum+=filt[i];
	  }
	  if ( NCALL_SMOOTH == 1 ) 
	    printf("\t ==> Smoothing with a tophat filter of %d bins (%d Ang) \n",  
		   l+1, (l+1)*(int)BINSIZE_LAMBDA );

	}
	else if (opt == OPT_SMOOTH_GAUSS ){
	  sum =0.;
	  for (i=0;i<l;i++){
	    t1 = pow((double)(k-i)*10.,2.);
	    t2 = pow(width,2.);
	    filt[i] = exp(-t1/(2.0*t2));
	    filt[i] = filt[i]/sqrt(2.*pi*t2);
	    sum+=filt[i];
	  }
	  if ( NCALL_SMOOTH == 1 ) 
	    printf("\t ==> Smoothing with Gaussian filter (width=%5.0f Ang) \n", width);
	}
	else {
	printsedSmooth_usage();
	}
	for (i=0;i<l;i++){
		filt[i] = filt[i]/sum;
	}
	return filt;

}

/*
setbox sets the temporary array around flux[j] or fluxerr[j]
taking into account padding choices used. 

Usage:  l = number of bins being used
	j = center bin value
	nBins = Tot number of bins in the main dataset (used to decide
		when to pad the arrays)
 	fl	= pointer to the array containing the variable (flux, or errors)
*/
double *setbox(int l, int j, int nBins, double *fl){
	static double filt[MXLAM_SMOOTH];
	int i, temp;


	for (i =0;i<l;i++){
		temp = j- (l-1)/2 +i;
		if (temp < 0) filt[i] = fl[0];//pad at beginning
		else if (temp>nBins -1) filt[i] =fl[nBins-1];//pad at end
		else filt[i] = fl[temp];//no padding required
	}
	return filt;
}
/*
printsedSmooth_usage is called when the opt is not defined as an algorith 
and prints out the usage of sedSmooth
*/
void printsedSmooth_usage(){
  printf("sedSmooth produces a smoothed spectrum and corresponding errors\n");
  printf("	by smoothing with a filter.\n ");
  printf("	INPUT\n");
  printf("	opt 	   =1 uses a tophat filter, opt=2 uses a Gaussian filter, any \n");
  printf("			other opt echoes the usage. \n");
  printf("	Nlam 	   = Number of observations in the original spectrum\n");
  printf("	lambda     = array (0,Nlam-1) wavelength taken to be at binnings of \n");
  printf("			10. A^o (change by changing binsize)\n");
  printf("	flux_in    = array (0,Nlam-1)  values of flux	\n");
  printf("	fluxerr_in = array (0,Nlam-1)  values of flux	\n");
  printf("	flux_smooth= array (0,Nlam -1) values of smoothed flux\n");
  printf("	fluxerr_smooth= array (0,Nlam -1) values of errors on the smoothed flux\n");
  printf("	smoothpar  = pointer to a width that you would like in angstron (This gets\n"); 
  printf("		translated number of bins used)\n");
  exit(EXIT_FAILURE);
}

