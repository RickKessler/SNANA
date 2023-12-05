/**************************************

 Created July 7, 2010 by R.Kessler

 utility to extract single spectrum (SED) from a SIMSED version
 for a specified set of SED parameters. The output SED is 
 interpolated based on the parameters defining the SED grid.

 Syntax:

 SIMSED_extractSpec.exe  --v <version> --t <Trest> --plist <param list>


 where 
   <version>     is the SIMSED version in $SNDATA_ROOT/models/SIMSED
                    or in $SNANA_MODELPATH
   <Trest>       is the epoch in days 
   <param list>  is a list of parameters describing the SIMSED version
                 (give NPAR float values without names)

 optional args
   --tlist <Trest-file>    file with list of Trest values and
                           optional output spec filenames, 
                           optional obs. frame lambda range, and
                           optional snr per lambin input. 

   --lam 2000 9200         rest frame lambda range: default is 1000 to 20000

   --z <z>                 transform to redshift = z: default is rest-frame

   --mu <mu>               use this distance modulus for redshift z

   --snr <snr> <lambin>    use this to set an ideal obs frame  snr per lambin
                           for flux and fluxerr output

   --seed <iseed>          Random  seed for flux-fluctuations from <snr>

   SNANA                   use this to get SNANA-style output format. Default 
                           output is ASCII ( plain 2 or 3 column, no headers).

   --T0 <T0>               use this to get observer date in SNANA-style 
                           headers. must include with SNANA format. 

   --dump <lambda>        dump flux for this wavelength (RK Dec 2011)

   LDUMP                   use this to have debugging lines dumped to screen

   OPT_ZEROFLUX           use ZFLUXVAL = 1e-20 for computing fluxerror in cases where flux is zero

 The extracted spectrum is written to a text file
 with format "wavelength  flux" and optional "dflux"
 if either SEDMODEL.FLUX_ERRFLAG or snr are non-zero. 

 Trest-file:
  <Trest1> <specFile1> <laminit1> <lamfinal1><snrlambin1><snr2000_1><snr3000_1><...><snr10000_1>
  <Trest2> <specFile2> <laminit2> <lamfinal2><snrlambin2><snr2000_2><snr3000_2><...><snr10000_2>
   ...
 where the spectrum is dumped into a text file called <specFile#>.
 !! Trest-file laminit and lamfinal values are in rest-frame !!
 !! Trest-file SNR values are in observer frame !!
 If <specFile#> exists, then the spectrum is appended so that
 header info can be prepared before calling this program.
 If <specFile#> does not exist then it is created. If the optional
 <specFile#> file names are not given, then default file-names
 SIMSED_[Trest].SPEC are used.

 If <laminit#> and <lamfinal#> exist, they are transformed to observer frame
 and used to mask the output spectrum data. If these values do not exist, the 
 default values of laminit and lamfinal are used instead. 

 If <snrlambin><snr values> exist, the appropriate SNR will be created at each
 wavelength, and linearly interpolated in between. A negative number should be
 used to indicate no value ( if for instance, the spectrum doesn't start until 3000A).
 If no values are given for particular wavelength slots, the nearest non-zero SNR
 value will be used. 

 Processing time to read and interpolate 4 SEDs 
 (2 parameter interp) is about 0.7 sec for a set
 of one dozen TREST values in a list. 


     HISTORY
  ~~~~~~~~~~~~

 Feb 08 2011 (JLM) - adapted getSpec_approx(void) to return exact day. 
                     still only finds nearest grid parameters.

 Feb 11, 2011 (RK) - include SEDMODEL.FLUXSCALE factor

 Feb 14, 2011 (RK) - add hook for --z <z> option

 Feb 14, 2011 (JLM) - fleshed out obsSpec(void), included --snr <snr> option

 Feb 19, 2011 RK - new '--tlist <file>' argument to give list of Trest values.
                   specFile argument is now ignored; output filename is
                   SIMSED_<Trest>.SPEC

 Apr 3, 2011 (RK) -
    - replace interp8() with interp_1DFUN(...)
    - print SNANA_DIR near top
    - finally do SED-parameter interpolation; see interp_prep().

 April 11, 2011 RK - allow optional spec filename to follow each
                     Trest value in the --tlist file.

 Jun 8, 2011 RK - use SNANA_MODELPATH if it is defined.
                  Allow SNDATA_ROOT to be undefined as long
                  as SNANA_MODELPATH is defined.
                  ABORT if neither SNDATA_ROOT nor SNANA_MODELPATH exist.

 Jul 17, 2011 JLM - new SNANA argument option allows SNANA-formatted output style.
                    new --T0 <t0> option allows peak date input. This option
                    is required with SNANA flag. 

 Aug 23, 2011 JLM - allow optional laminit and lamfinal values to follow
                    each Trest filename pair in the --tlist file. 

 Aug 24, 2011 RK - fix parse_args() to abort on any undefined argument. 

 Aug 25, 2011 JLM - 
    Replaced --snr <SNR> with --snr <SNR> <SNRLAMBIN>
    Upgraded fudgeSNR to appropriately scale SNR with SED binsize. 
    Upgraded fudgeSNR to randomize flux in keeping with fluxerror. 
    Added --seed <seed> to allow random seed to be input by user. 

  RK fix a few bugs
   - INPUTS.specLAMinit[200] -> [MXEPOCH]
   - INPUTS.specLAMfinal[200] -> [MXEPOCH]
   - remove NARGV_LIST declaration inside parse_args since NARGV_LIST is global
   - in fudgeSNR() set FLUXERR to measured instead of IDEAL
   - remove --snrlamrange arg and require two args after '--snr'

  JM fix one more bug
   - Altered write_header() to accurately reflect number of spec lines output

Sep 22, 2011 JLM - 
   Allow optional SNRLAMBIN, SNR values to follow Trest,laminit,lamfinal 
   input in the --tlist file. SNR values should go from 2000 - 10000 Angstroms, 
   in 1000 Angstrom steps. 

Nov 21, 2011 RK:
     if a parameter value is -999 then ignore it.
     Allows skipping baggage parameters.
     See new logical array INPUTS.USEPAR[ipar]

Nov 21, 2011 JLM:
     add galaxy inputs --gfile and --gfrac for
     passing in a galaxy spectrum and percent
     galaxy contamination respectively. 
     contamination is defined at peak, between 5000-6000 Angstroms. 

 Dec 10, 2011 RK: 
     add --dump <lamdump> option: needed for SNANA_tester

 Dec 12, 2011 JLM:
     add OPT_ZEROFLUX option (with ZFLUXVAL) to compute fluxerrs when flux is 0.0
     add LDUMP option to dump debugging lines 
     fixed bug in fudgeSNR() allowing negative fluxerrors when flux is negative
           

 Jan 2, 2012 RK
   specLAMinit -> specLAMinit_rest and specLAMinit_obs
   and same for specLAMfinal. Fixes bug in which the rest-frame
   range in the --tlist file was used on the obs. lambda.

   In parse_epochs(), abort if Snrlambin <= 0

 Jan 4, 2012 JLM
    created new SPEC_LAMBIN value which is redshifted along with 
    SPECLAM, SPECFLUX values, and used in fudgeSNR to get correct
    (observer frame) wavelength bin size. 

 Jan 10, 2012 JLM
     altered read_simsed to read sufficient range of epochs when galaxy 
     contamination is desired 

 Feb 21, 2013 RK - update to work with c++ compiler.

 May 23, 2013 JLM - defined MXCHAR_LINE to set max EPOCHFILE line length

 Aug 2017 SEDMODEL.FLUX_ERRFLAG replaced with SEDMODEL.OPTMASK
          See new local global FLUX_ERRFLAG.
             BEWARE:  NOT TESTED

 Jun 4 2020: call init_random_seed(...)

**************************************/

#include "genmag_SIMSED.c"


// =========== DECLARE FUNCTIONS ==================

void  init_inputs(void);
void  parse_args(int argc, char **argv);
void  parse_epochs(char *epFile);
void  write_banner(int ep);
void  write_header(int ep);
void  write_footer(int ep);
void  prep_gal(void);
void  add_gal(void);
void  interp_prep(void);
void  read_simsed(void);
void  getSpec(int ep);
void  init_SIMSED_PATH(void);
void  obsSpec(void);
void  fudgeSNR(int ep);
void init_getSNR(int *wsnr, int *wlam, double lam, int ep);
double getSNR(int *wsnr, int *wlam, double lam, int ep);
void  wrSpec(int ep);

int ISED_MATCH(double *parval);

// =========== DECLARE GLOBALS =================

#define NXSNR 9 // no. of SNR inputs (2000A - 10000A)
#define MXPARAM 20    // max no. of parameters
#define MXCORNERS 100 // max no. corners stored for interpolation
#define ZFLUXVAL 1e-20 // use to compute fluxerror when flux is zero
#define MXCHAR_LINE 600 // max length of EPOCHFILE lines

#define MXEPOCH_SIMSED 2000

struct INPUTS {
  char   SIMSED_VERSION[200];

  int    NEPOCH ;
  int    ISEED ;
  double EPOCH[MXEPOCH_SIMSED];
  double MINEPOCH, MAXEPOCH ;
  char   EPOCHFILE[200];  // optional file with list of epochs
  char   char_EPOCH[MXEPOCH_SIMSED][40] ;  // used to create specFile name

  char   GALFILE[200]; // optional galaxy file

  int    USEPAR[MXPARAM];
  double PARAM_LIST[MXPARAM];

  double LAM_RANGE[2];

  double specLAMinit_rest[MXEPOCH_SIMSED];
  double specLAMfinal_rest[MXEPOCH_SIMSED];
  double specLAMinit_obs[MXEPOCH_SIMSED];
  double specLAMfinal_obs[MXEPOCH_SIMSED];

  double specSNR[MXEPOCH_SIMSED][NXSNR]; 
  double specSNRLAMBIN[MXEPOCH_SIMSED];

  char   specFile[MXEPOCH_SIMSED][200];
  char   TYPE[20];
  char   ZFLUX[20];
  double REDSHIFT, DLMAG, SNR, SNRLAMBIN;
  double PEAKMJD;
  double GALFRAC;

  double LAMDUMP ;

  int    DEBUG ;

} INPUTS ;

// char SIMSED_PATHMODEL[200];
char specFile[200]; // name of output spec file

double SNRLAM[NXSNR];
  
// define spectrum
int NBIN_LAM;
double SPEC_LAMSTEP;
double SPECLAM[MXBIN_LAMSED_SEDMODEL];
double SPECFLUX[MXBIN_LAMSED_SEDMODEL];
double SPECFLUXERR[MXBIN_LAMSED_SEDMODEL];

// define galaxy
int NBIN_GALLAM;
double GALLAM[MXBIN_LAMSED_SEDMODEL];
double GALFLUX[MXBIN_LAMSED_SEDMODEL];
double GALNORM;

struct CORNER {
  int NSTORE;

  double WGT[MXCORNERS] ;
  double WGTSUM ;
  int    ISED[MXCORNERS] ;

  double *FLUX[MXCORNERS] ;
  double *FLUXERR[MXCORNERS] ;

} CORNER ;

int FLUX_ERRFLAG ;

// ====================================
int main(int argc, char **argv) {

  char fnam[] = "main" ;
  int ep, i;
  
  // ------------ BEGIN ------------

  set_EXIT_ERRCODE(EXIT_ERRCODE_UNKNOWN);

  sprintf(BANNER,"Begin execution of SIMSED_extractSpec.exe" );
  print_banner(BANNER);

  SNRLAM[0] = 2000.;
  for ( i=1; i<NXSNR; i++){
    SNRLAM[i]= SNRLAM[i-1] + 1000;
  }

  printf("   SNANA_DIR = %s \n", getenv("SNANA_DIR") );
  init_inputs();
  parse_args(argc, argv );

  printf("   Init randoms, SEED = %d \n", INPUTS.ISEED);

  init_random_seed(INPUTS.ISEED, 1);

  // prepare the corners and SEDs needed to interpolate
  // the SEDs in SED-parameter space
  interp_prep();

  // read in the SED surface(s) for this set of --plist parameters
  read_simsed();

  // if GALFRAC larger than zero, read & normalize GAL FLUX
  prep_gal();

  //  debugexit("parse_epochs"); // xxxxxxxxxxx

  for ( ep = 1; ep <= INPUTS.NEPOCH; ep++ ) {

    // get rest-frame spectrum 
    getSpec(ep);

    // if GALFRAC gt 0, add gal contamination
    add_gal();

    // (optional) transform to observer-frame
    obsSpec();

    // if SNR larger than zero, fudge flux error
    fudgeSNR(ep);

    // write Hbanner if INPUTS.TYPE is SNANA and ep is 1
    write_banner(ep);

    // write epoch banner if INPUTS.TYPE is SNANA 
    write_header(ep);

    // write spectrum to text file
    wrSpec(ep);

    // write footer to text file
    write_footer(ep);
  }

  printf(" PROGRAM ENDING GRACEFULLY \n");

} // end of main



// **********************
void init_inputs(void) {

  int ep;
  int iSNR;
  char fnam[] = "init_inputs";

  INPUTS.ISEED = 12345;
  INPUTS.NEPOCH = 0;
  INPUTS.MINEPOCH = +999.0 ;
  INPUTS.MAXEPOCH = -999.0 ;

  sprintf(INPUTS.GALFILE, "%s", "NULL");

  for ( ep=0; ep < MXEPOCH_SIMSED; ep++ ) {
    INPUTS.EPOCH[ep] = -999. ; 
    sprintf(INPUTS.specFile[ep], "%s", "NULL" );
  }

  INPUTS.LAM_RANGE[0]   =  1000.0 ;
  INPUTS.LAM_RANGE[1]   = 20000.0 ;
  INPUTS.REDSHIFT =  0.0 ;
  INPUTS.DLMAG    =  0.0 ;
  INPUTS.SNR      = -101.0; 
  INPUTS.SNRLAMBIN = 100.0; // Angstroms
  INPUTS.PEAKMJD  = -1.0;
  INPUTS.GALFRAC = 0.0;
  sprintf(INPUTS.TYPE,"%s", "ASCII"); // ASCII is default

  INPUTS.DEBUG = 0 ;

  for ( ep=0; ep < MXEPOCH_SIMSED; ep++ ) {
    INPUTS.specLAMinit_rest[ep]  = INPUTS.LAM_RANGE[0];
    INPUTS.specLAMfinal_rest[ep] = INPUTS.LAM_RANGE[1];
    for ( iSNR=0; iSNR < NXSNR; iSNR++) {
      INPUTS.specSNR[ep][iSNR] = INPUTS.SNR;
    }
    INPUTS.specSNRLAMBIN[ep] = INPUTS.SNRLAMBIN; 
  }

  INPUTS.LAMDUMP = -999. ;

} // end of init_inputs


// ** START NEW PARSE ARGS **
// ******************************************
void parse_args(int argc, char **argv) {

  int i, ipar,  ep, ilast, iuse ;
  char fnam[] = "parse_args" ;
  double Trest, parval, z1 ;

  // ---------- BEGIN --------


  NARGV_LIST = argc;

  for( i=0; i < NARGV_LIST; i++ ) {
    USE_ARGV_LIST[i] = 0 ;
    ARGV_LIST[i] = (char*) malloc( MXPATHLEN*sizeof(char) );
    sprintf(ARGV_LIST[i], "%s", argv[i] );
  }

  i=1; ilast=1 ;

  SEDMODEL.NPAR = 0;
  
  while ( i < NARGV_LIST ) {
   
    if ( strcmp(ARGV_LIST[i],"--v") == 0 ) {
      i++;  sscanf(ARGV_LIST[i], "%s", INPUTS.SIMSED_VERSION );
      init_SIMSED_PATH();
      read_SIMSED_INFO(SIMSED_PATHMODEL);
      FLUX_ERRFLAG = ( SEDMODEL.OPTMASK & OPTMASK_FLUXERR_SEDMODEL) ;
    }

    if ( strcmp(ARGV_LIST[i],"--t") == 0 ) {
      INPUTS.NEPOCH++ ;
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.EPOCH[1] );
    }

    if ( strcmp(ARGV_LIST[i],"--tlist") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%s",  INPUTS.EPOCHFILE) ;
      parse_epochs(INPUTS.EPOCHFILE);
    }

    if ( strcmp(ARGV_LIST[i],"--plist") == 0 ) {
      if ( SEDMODEL.NPAR == 0 ) {
        sprintf(c1err,"Cannot read plist before SIMSED version is set.");
        sprintf(c2err,"Check syntax at top of SIMSED_extractSpec.c");
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

      for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {
        i++ ; sscanf(ARGV_LIST[i], "%le",  &parval ); 
        INPUTS.PARAM_LIST[ipar] = parval; 
	if ( parval == -999. )
	  { INPUTS.USEPAR[ipar] = 0 ; }
	else
	  { INPUTS.USEPAR[ipar] = 1 ; }
      }
    } // --plist


    // now check optionalargs

    if ( strcmp(ARGV_LIST[i],"--gfile") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%s",  INPUTS.GALFILE) ;
    }

    if ( strcmp(ARGV_LIST[i],"--gfrac") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.GALFRAC );
    }

    if ( strcmp(ARGV_LIST[i],"--iseed") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%d", &INPUTS.ISEED );
    }
    if ( strcmp(ARGV_LIST[i],"--seed") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%d", &INPUTS.ISEED );
    }

    if ( strcmp(ARGV_LIST[i],"--lam") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.LAM_RANGE[0] );
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.LAM_RANGE[1] );
    }

    if ( strcmp(ARGV_LIST[i],"--z") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.REDSHIFT );
    }
    if ( strcmp(ARGV_LIST[i],"--redshift") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.REDSHIFT );
    }

    if ( strcmp(ARGV_LIST[i],"--mu") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.DLMAG );
    }

    if ( strcmp(ARGV_LIST[i],"--dm") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.DLMAG );
    }

    if ( strcmp(ARGV_LIST[i],"--snr") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.SNR );
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.SNRLAMBIN );
    }
    
    if ( strcmp(ARGV_LIST[i],"SNANA") == 0 ) {
      sprintf(INPUTS.TYPE, "%s", "SNANA");
      USE_ARGV_LIST[i] = 1; 
    }

    if ( strcmp(ARGV_LIST[i],"--T0") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.PEAKMJD );
    }

    if ( strcmp(ARGV_LIST[i],"--dump") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i], "%le",  &INPUTS.LAMDUMP ) ;
    }

    if ( strcmp(ARGV_LIST[i],"LDUMP") == 0 ) {
      INPUTS.DEBUG = 1 ;
      USE_ARGV_LIST[i] = 1; 
    }

    if ( strcmp(ARGV_LIST[i],"OPT_ZEROFLUX") == 0 ) {
      sprintf(INPUTS.ZFLUX, "%s", "ZFLUX");
      USE_ARGV_LIST[i] = 1; 
    }

    // -------------------
   
    if ( i > ilast ) {
      for ( iuse = ilast; iuse <= i; iuse++ )
        { USE_ARGV_LIST[iuse] = 1; }
    }
    i++ ;  ilast=i;

  } // end while


  // make sure that every argument is defined.

  check_argv();


  // -------------

  if ( strcmp(INPUTS.TYPE,"SNANA") == 0 && INPUTS.PEAKMJD == -1 ) {
    sprintf(c1err,"Must include option --T0 <peakmjd> with SNANA argument.");
    sprintf(c2err,"Fix input parameters ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( INPUTS.GALFRAC != 0.0 && strcmp(INPUTS.GALFILE, "NULL") == 0 ) {
    sprintf(c1err,"To use --gfrac must also use option --gfile <galaxy file>.");
    sprintf(c2err,"Fix input parameters ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  printf("\n");
  printf(" Extract spectrum from '%s' for \n", INPUTS.SIMSED_VERSION );

  z1 = 1.0 + INPUTS.REDSHIFT ;
  for ( ep=1; ep <= INPUTS.NEPOCH; ep++ ) {
    Trest = INPUTS.EPOCH[ep] ;
    printf("\t Trest  \t = %8.3f days \n", Trest );

    if ( Trest < INPUTS.MINEPOCH )
      { INPUTS.MINEPOCH = Trest ; }
    if ( Trest > INPUTS.MAXEPOCH )
      { INPUTS.MAXEPOCH = Trest ; }

    INPUTS.specLAMinit_obs[ep]  = z1 * INPUTS.specLAMinit_rest[ep];
    INPUTS.specLAMfinal_obs[ep] = z1 * INPUTS.specLAMfinal_rest[ep];

  }
  
  if ( INPUTS.DEBUG  ) {
    printf("\t xxx DUMP: MIN,MAX epoch = %f, %f \n",
	   INPUTS.MINEPOCH, INPUTS.MAXEPOCH );
  }
 

  for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {
    printf("\t %s   \t = %f ",
           SEDMODEL.PARNAMES[ipar], INPUTS.PARAM_LIST[ipar] );
    if ( INPUTS.USEPAR[ipar] == 0 ) { printf(" ==> IGNORE"); }
    printf("\n");
  }

  printf("\t Rest-frame wavelength \t = %6.0f to %6.0f \n",
         INPUTS.LAM_RANGE[0], INPUTS.LAM_RANGE[1] );
  printf("\t FLUX_SCALE \t = %10.3le \n",
         SEDMODEL.FLUXSCALE );
  printf("\t FLUX_ERRFLAG\t = %d \n", FLUX_ERRFLAG); 

  if ( INPUTS.REDSHIFT > 0.0 ) {
    printf("\t Transform to z=%5.3f and distMod=%6.3f \n",
           INPUTS.REDSHIFT, INPUTS.DLMAG );
  }

  if ( INPUTS.SNR > -100.0 ) {
    printf("\t Fudge flux & error with S/N = %5.2f per %5.1f A \n",
           INPUTS.SNR, INPUTS.SNRLAMBIN );
  }

  if ( INPUTS.GALFRAC > 0.0 ) {
    printf("\t Add %5.1f percent galaxy contamination using template %s \n", 
	   INPUTS.GALFRAC * 100, INPUTS.GALFILE );
  }

  if ( strcmp(INPUTS.ZFLUX, "ZFLUX") == 0 ) {
    printf("\t Will use flux = %e to compute fluxerror if flux == 0.0 \n", 
      ZFLUXVAL);
  }

  printf("\t Ouput spectra in %s format. \n\n", INPUTS.TYPE);


} // end of parse_args
// ** END NEW PARSE ARGS **



// ********************************
void parse_epochs(char *epFile) {
  
  // read Trest and optional spectrum filename on each line of *epFile.

  // format of Trest-file:
  //<Trest1> <specFile1> <laminit1> <lamfinal1><snrlambin1><snr2000_1><snr3000_1><...><snr10000_1> 
  //<Trest2> <specFile2> <laminit2> <lamfinal2><snrlambin2><snr2000_2><snr3000_2><...><snr10000_2>
  // laminit, lamfinal are assumed to be in rest frame. 
  // snr information is taken to be in observer frame. 
  
  FILE *fp;
  char 
    *ptrtok
    ,line[MXCHAR_LINE]
    ,tmpline[MXCHAR_LINE]
    ,strTrest[100], specFile[200]
    ,strLaminit[100], strLamfinal[100] 
    ,strSNR[100], strSNRLAMBIN[100]
    ,fnam[] = "parse_epochs";
    ;

  int N, iSNR;
  double 
    Trest
    ,Laminit, Lamfinal
    ,Snr, Snrlambin
    ,z1
    ;

  // ------------- BEGIN --------------

  fp = fopen(epFile, "rt") ;

  if ( !fp ) {
    sprintf(c1err,"Could not open epoch file: %s", epFile);
    sprintf(c2err,"Check --tlist argument");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  while ( fgets (line, MXCHAR_LINE, fp) !=NULL  ) {
    // skip blank lines
    if ( strlen(line) == 0 ) { continue ; }

    // break of tmpline into blank-separated strings
    sprintf(tmpline,"%s", line);
    ptrtok = strtok(tmpline," ");

    // 1st string
    sscanf ( ptrtok, "%s", strTrest );
    sscanf(strTrest, "%le" , &Trest ) ;

    // get 2nd string 
    ptrtok = strtok(NULL, " ");
    if ( ptrtok != NULL ) 
      { sscanf ( ptrtok, "%s", specFile ); }
    else
      { sprintf(specFile, "%s", "NULL" ); }

    // get 3rd string
    ptrtok = strtok(NULL, " ");
    if ( ptrtok != NULL )
      { sscanf (ptrtok, "%s", strLaminit);
        sscanf(strLaminit, "%le", &Laminit); }
    else
      { Laminit = INPUTS.LAM_RANGE[0]; }
    
    // get 4th string
    ptrtok = strtok(NULL, " ");
    if ( ptrtok != NULL ) { 
      sscanf (ptrtok, "%s", strLamfinal); 
      sscanf(strLamfinal, "%le", &Lamfinal ) ; 
    }
    else
      { Lamfinal = INPUTS.LAM_RANGE[1]; }
    
    // get 5th string
    ptrtok = strtok(NULL, " ");
    if ( ptrtok != NULL ) { 
      sscanf (ptrtok, "%s", strSNRLAMBIN); 
      sscanf(strSNRLAMBIN, "%le", &Snrlambin); 

      // abort on invalid Snrlambin (RK, Jan 2 2012)
      if ( Snrlambin <= .0001 ) {
	sprintf(c1err,"Invalid Snrlambin='%f'", Snrlambin );
	sprintf(c2err,"Check 5th argument of file '%s' " , epFile );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
    }
    else
      { Snrlambin = INPUTS.SNRLAMBIN; }


    if ( INPUTS.DEBUG ) {
      printf(" xxx DUMP: Trest = %7.3f  specFile = '%s' \n", 
	     Trest, specFile );
      printf(" xxx DUMP: LamRest(init,final) = %7.3f, %7.3f\n", 
	     Laminit, Lamfinal  );
      printf(" xxx DUMP: SNRLAMBIN=%7.3f \n", Snrlambin );
    }

    // store Trest, specFile, Laminit, and Lamfinal in global INPUTS structure.
    INPUTS.NEPOCH++ ;
    N = INPUTS.NEPOCH ;
    sprintf(INPUTS.specFile[N],   "%s", specFile) ;
    sprintf(INPUTS.char_EPOCH[N], "%s", strTrest );
    INPUTS.EPOCH[N] = Trest ;
    INPUTS.specLAMinit_rest[N]  = Laminit;
    INPUTS.specLAMfinal_rest[N] = Lamfinal;

    INPUTS.specSNRLAMBIN[N] = Snrlambin; 

    // get 6th - 14th strings
    if ( strcmp(INPUTS.TYPE,"SNANA") == 0 ) {printf(" xxx SNR vals: ");}
    for (iSNR=0; iSNR < NXSNR; iSNR++){
      ptrtok = strtok(NULL, " ");
      if ( ptrtok != NULL )
	{ sscanf (ptrtok, "%s", strSNR); 
	  sscanf(strSNR, "%le", &Snr); }
      else
	{ Snr = INPUTS.SNR; }

      INPUTS.specSNR[N][iSNR] = Snr;
      if ( strcmp(INPUTS.TYPE,"SNANA") == 0 ) {printf("%d %7.3f ",iSNR,  Snr );}
    }
    printf("\n");
  }

  fclose(fp);

  //  debugexit("parse epochs"); // xxxxx

} // end of parse_epochs

// ********************************
void write_banner(int ep) {
  
  // read Trest and optional spectrum filename on each line op *epFile.
  
  FILE *fp;
  char 
    *ptrtok
    ,line[200]
    ,tmpline[200]
    ,strTrest[100]
    ,specFile[200]
    ,ctmp[80]
    ,fnam[] = "write_banner" ;
    ;

  int 
    istat
    ,N
    ,ERRFLAG
    ;
  double Trest ;

  // ------------- BEGIN --------------


  if ( strcmp(INPUTS.TYPE,"SNANA") != 0 || ep > 1 ) { return; }

  if ( strcmp(INPUTS.specFile[ep],"NULL") != 0 ) {
    sprintf(specFile, "%s", INPUTS.specFile[ep] );
    // if file exists, then append; else create new file.
    istat = access(specFile, F_OK);
  }
  else {
    istat = -1;
    if ( INPUTS.NEPOCH == 1 ) {
      sprintf(specFile, "%s", "SIMSED.SPEC" );
    }
    else {
      sprintf(specFile,"SIMSED_%s.SPEC", INPUTS.char_EPOCH[ep] );
    }
  }

  if ( istat == 0 )  { 
    fp = fopen(specFile, "at") ;  // append
    sprintf(ctmp,"%s" , "Append" );
  }
  else {
    fp = fopen(specFile, "wt") ;    // new file
    sprintf(ctmp,"%s" , "Write" );
  }


  printf(" %s Hbanner to : %s \n", ctmp, specFile );

  ERRFLAG = (FLUX_ERRFLAG  || INPUTS.SNR > 1.0E-9 ) ;

  fprintf(fp, "#============\n");
  fprintf(fp, "# \n");
  fprintf(fp, "# Simulated spectra below.  \n");
  fprintf(fp, "# \n");
  fprintf(fp, "#============\n\n");

  fprintf(fp, "NSPEC: %d \n\n", INPUTS.NEPOCH );
  
  if ( ERRFLAG == 1 ) {
    fprintf(fp, "NVAR_SPEC: 3\n");
    fprintf(fp, "VARNAMES_SPEC: WAVELENGTH FLUX FLUXERR\n\n");
  } else {
    fprintf(fp, "NVAR_SPEC: 2\n");
    fprintf(fp, "VARNAMES_SPEC: WAVELENGTH FLUX \n\n");
  }

  fclose(fp);


} // end of write_banner

// ********************************
void write_header(int ep) {
  
  // If INPUTS.TYPE is SNANA, write banner for each epoch
  
  FILE *fp;
  char 
    *ptrtok
    ,line[200]
    ,tmpline[200]
    ,strTrest[100]
    ,specFile[200]
    ,ctmp[80]
    ;

  int 
    istat
    ,ilam
    ,Nlam
    ;
  double Tobs, z, Laminit, Lamfinal ;

  // ------------- BEGIN --------------

  if ( strcmp(INPUTS.TYPE,"SNANA") != 0 ) { return; }

  if ( strcmp(INPUTS.specFile[ep],"NULL") != 0 ) {
    sprintf(specFile, "%s", INPUTS.specFile[ep] );
    // if file exists, then append; else create new file.
    istat = access(specFile, F_OK);
  }
  else {
    istat = -1;
    if ( INPUTS.NEPOCH == 1 ) {
      sprintf(specFile, "%s", "SIMSED.SPEC" );
    }
    else {
      sprintf(specFile,"SIMSED_%s.SPEC", INPUTS.char_EPOCH[ep] );
    }
  }

  // count number of output lambdas for use in header
  Laminit  = INPUTS.specLAMinit_obs[ep];
  Lamfinal = INPUTS.specLAMfinal_obs[ep];
  Nlam = 0; 
  for ( ilam = 0; ilam < NBIN_LAM; ilam++) {
    if ( SPECLAM[ilam] <= Lamfinal && SPECLAM[ilam] >= Laminit ) {
      Nlam++;
    }
  }

  if ( istat == 0 )  { 
    fp = fopen(specFile, "at") ;  // append
    sprintf(ctmp,"%s" , "Append" );
  }
  else {
    fp = fopen(specFile, "wt") ;    // new file
    sprintf(ctmp,"%s" , "Write" );
  }

  z = INPUTS.REDSHIFT;
  Tobs = INPUTS.EPOCH[ep]*(1+z)+INPUTS.PEAKMJD;

  fprintf(fp, "#--------------\n\n");
  fprintf(fp, "SPECNUM: %d  \n", ep );
  fprintf(fp, "NSPECBIN_WAVELENGTH: %d \n", Nlam);
  fprintf(fp, "SPECOBSMJD: %0.3f  \n", Tobs );
  fprintf(fp, "SPECUNITS: STD (erg/s/A/cm^2)  \n\n");

  fclose(fp);

} // end of write_header

// ******************************
void init_SIMSED_PATH(void) {

  char 
    *VERSION 
    ,SNDATA_ROOT[200]
    ,PRIVATE_PATH[200] 
    ;

  int EXIST_PRIVATE_PATH, EXIST_SNDATA_ROOT ;
  char fnam[] = "init_SIMSED_PATH" ;

  // ------------ BEGIN --------

  VERSION = INPUTS.SIMSED_VERSION ;

  if ( getenv("SNDATA_ROOT") != NULL ) {
    sprintf(SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );
    EXIST_SNDATA_ROOT = 1;
  }
  else 
    { EXIST_SNDATA_ROOT = 0 ; }

  if ( getenv(PRIVATE_MODELPATH_NAME) != NULL ) {
    sprintf(PRIVATE_PATH,  "%s", getenv(PRIVATE_MODELPATH_NAME) );
    EXIST_PRIVATE_PATH = 1 ;
  }
  else 
    { EXIST_PRIVATE_PATH = 0 ; }


  if ( EXIST_SNDATA_ROOT == 0 && EXIST_PRIVATE_PATH == 0 ) {
    sprintf(c1err,"Neither $SNDATA_ROOT or $%s are defined.", 
	    PRIVATE_MODELPATH_NAME );
    sprintf(c2err,"At least one of these must be defined with setenv.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  if ( EXIST_PRIVATE_PATH ) {
    sprintf( SIMSED_PATHMODEL, "%s/%s", PRIVATE_PATH, VERSION );
    printf("   %s = %s \n", PRIVATE_MODELPATH_NAME, PRIVATE_PATH  );
  }
  else {
    sprintf( SIMSED_PATHMODEL, "%s/models/SIMSED/%s", 
             SNDATA_ROOT, VERSION );
    printf("   %s: SNDATA_ROOT = %s \n", fnam, SNDATA_ROOT );
    fflush(stdout);
  }
  
  printf(" xxx %s: done here. \n", fnam);
  fflush(stdout);

} // end of init_SIMSED_PATH


// ============================================
void interp_prep(void) {

  char fnam[] = "interp_prep" ;

  int 
    ipar, istat, I0SED, I1SED, ISED, RDFLAG
    ,icorner
    ,NCORNERS
    ,NSTORE
    ,MSK, J01
    ,LDMP_CORNER = 1
    ;

  double 
    parval, p0, p1, WGT
    ,NEAREST_GRIDVAL[MXPARAM][2]
    ,GRIDFRAC[MXPARAM]
    ,PARVAL_CORNER[MXPARAM]
    ;

  char ctmp[80];

  // ------------- BEGIN --------------

  NCORNERS = 1;

  for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {

    if ( INPUTS.USEPAR[ipar] == 0 ) { continue ; }

    parval = INPUTS.PARAM_LIST[ipar];

    NCORNERS *= 2 ;

    istat  = get_SEDMODEL_INDICES( ipar, parval, &I0SED, &I1SED ); 

    if ( I1SED == I0SED ) {
      parval += 1.0E-8 ;
      istat  = get_SEDMODEL_INDICES( ipar, parval, &I0SED, &I1SED ); 
    }

    p0 = SEDMODEL.PARVAL[I0SED][ipar] ;
    p1 = SEDMODEL.PARVAL[I1SED][ipar] ;
    NEAREST_GRIDVAL[ipar][0] = p0 ;
    NEAREST_GRIDVAL[ipar][1] = p1 ;

    GRIDFRAC[ipar] = (parval - p0)/(p1-p0);

    /*
    printf(" xxxx %12s = %7.3f => bounded by %6.3f and %6.3f (istat=%d)\n"
	   ,SEDMODEL.PARNAMES[ipar], parval
	   ,NEAREST_GRIDVAL[ipar][0]
	   ,NEAREST_GRIDVAL[ipar][1] 
	   ,istat
	   ) ;
    */
  }


  CORNER.WGTSUM = 0.0 ;
  NSTORE        = 0 ;
  int NPAR_USE  = 0 ;

  for ( icorner = 0; icorner < NCORNERS; icorner++ ) {

    if ( LDMP_CORNER ) 
      { sprintf(ctmp, "CORNER %3.3d (", icorner ); }

    WGT = 1.0 ;
    NPAR_USE = 0;
    for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {
      if ( INPUTS.USEPAR[ipar] == 0 ) { continue ; }
      NPAR_USE++;

      MSK = 1 << (NPAR_USE-1);
      J01 = (icorner & MSK)/MSK ; // 0 or 1 only
      PARVAL_CORNER[ipar] = NEAREST_GRIDVAL[ipar][J01] ;

      if ( J01 ) 
        { WGT *= GRIDFRAC[ipar] ; }
      else
        { WGT *= (1.0 - GRIDFRAC[ipar]) ; }

      if ( LDMP_CORNER ) { sprintf(ctmp,"%s%d ", ctmp, J01); }
    } // ipar


    ISED   = ISED_MATCH(PARVAL_CORNER);
    RDFLAG = (WGT > 0.00001) ;

    if ( RDFLAG ) {
      CORNER.ISED[NSTORE]   = ISED ;
      CORNER.WGT[NSTORE]    = WGT ;
      CORNER.WGTSUM        += WGT ;
      NSTORE++;
    }
    
    if ( LDMP_CORNER && RDFLAG ) { 
      sprintf(ctmp, "%s)", ctmp ) ;
      printf("\t %s ISED=%4d  WGT=%10.6f (%s) \n", 
	     ctmp, ISED, WGT, SEDMODEL.FILENAME[ISED] ) ;
    }

  } // icorner
  
  CORNER.NSTORE  = NSTORE ;

  if ( NSTORE >= MXCORNERS ) {
    sprintf(c1err,"CORNER.NSTORE = %d exceeds bound of ", NSTORE);
    sprintf(c2err,"MXCORNERS = %d", MXCORNERS );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( fabs(CORNER.WGTSUM-1.0) > 0.001 ) {
    sprintf(c1err, "CORNER.WGTSUM = %f", CORNER.WGTSUM );
    sprintf(c2err, "%s", "Check individual weights.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // debugexit(fnam);

} // end of interp_prep

// ============================================
void add_gal(void) {

  // Add galaxy flux to existing SED flux

  double interp_lam; 
  double interp_galflux;
  int ilam;
  char fnam[] = "add_gal" ;

  if ( INPUTS.GALFRAC <= 0.0 ) {
    return;
  }

  //Calculate normalization
  for ( ilam=0; ilam < NBIN_LAM; ilam++ ) {
    interp_lam = SPECLAM[ilam];
      interp_galflux = 
	interp_1DFUN(2, interp_lam, NBIN_GALLAM, GALLAM, GALFLUX, fnam );
      SPECFLUX[ilam] = SPECFLUX[ilam] + GALNORM*interp_galflux; 
  }
  
 
} // end of add_gal


// ============================================
void prep_gal(void) {

  // Read in galaxy SED information ( lambda, flux)
  // corresponding to --gfile. 
  // Store results in globals GALLAM, GALFLUX

  double gal_sum, sed_sum;
  double interp_lam; 
  double interp_galflux;
  int epoch = -1;
  int ilam;
  char fnam[] = "prep_gal" ;

  if ( INPUTS.GALFRAC <= 0.0 ) {
    return;
  }

  // Get galaxy SED
  rd2columnFile(INPUTS.GALFILE, MXBIN_LAMSED_SEDMODEL, &NBIN_GALLAM , GALLAM , GALFLUX, 0 );
  printf(" Read %d rows of galaxy SED \n\n", NBIN_GALLAM );

  // Get peak SN spectrum
  // use negative epoch to force peak 
  getSpec(epoch); 

  //Calculate normalization
  gal_sum = 0.0;
  sed_sum = 0.0;
  for ( ilam=0; ilam < NBIN_LAM; ilam++ ) {
    interp_lam = SPECLAM[ilam];
    if ( interp_lam >= 5000. && interp_lam <= 6000. ) {
      interp_galflux = 
	interp_1DFUN(2, interp_lam, NBIN_GALLAM, GALLAM, GALFLUX, fnam );
      gal_sum = gal_sum + interp_galflux;
      sed_sum = sed_sum + SPECFLUX[ilam]; 
    }
  }
  GALNORM = INPUTS.GALFRAC * ( sed_sum / gal_sum );
  if ( INPUTS.DEBUG ) { printf(" GALNORM is %e \n", GALNORM);}
  
 
} // end of prep_gal

// ============================================
void read_simsed(void) {

  // Read simsed spectral surface (lambda vs. Trest) corresponding
  // to --plist.  Store results in TEMP_SEDMODEL struct.
  //
  // This function needs modification to propery interpolate
  // from SIMSED parameters

  int ised, icorner, nflux_nan ;

  double TREST_RANGE[2] ;

  char 
    sedFile_full[200]
    ,sedcomment[100]
    ,fnam[] = "read_simsed" 
    ;

  // ------------- BEGIN --------------

  TREST_RANGE[0] =  INPUTS.MINEPOCH - 2.0 ;
  TREST_RANGE[1] =  INPUTS.MAXEPOCH + 2.0 ;

  // Need to ensure that spectral surface extends to 0.0 ep,
  // in case of galaxy contamination
  if ( INPUTS.GALFRAC > 0.0 ) {
    if ( INPUTS.MINEPOCH > 0.0 ) TREST_RANGE[0] = -2.0;
    if ( INPUTS.MAXEPOCH < 0.0 ) TREST_RANGE[1] = 2.0;
  }
 

  for ( icorner = 0; icorner < CORNER.NSTORE ; icorner++ ) {

    ised = CORNER.ISED[icorner];

    sprintf(sedFile_full,"%s/%s", 
	    SIMSED_PATHMODEL, SEDMODEL.FILENAME[ised] );

    sprintf(sedcomment,"(ised=%d/%d)", ised, SEDMODEL.NSURFACE );

    // allocate memory for this SED
    CORNER.FLUX[icorner]    = (double *)malloc(8*MXBIN_SED_SEDMODEL);
    CORNER.FLUXERR[icorner] = (double *)malloc(8*MXBIN_SED_SEDMODEL);

    malloc_SEDFLUX_SEDMODEL(&TEMP_SEDMODEL,0,0,0);

    rd_sedFlux(sedFile_full, sedcomment
	       ,TREST_RANGE, INPUTS.LAM_RANGE
	       ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL
	       ,SEDMODEL.OPTMASK
	       ,&TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY, &TEMP_SEDMODEL.DAYSTEP
	       ,&TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM, &TEMP_SEDMODEL.LAMSTEP
	       ,CORNER.FLUX[icorner], CORNER.FLUXERR[icorner]
	       ,&nflux_nan);
  }

  
  //  debugexit(fnam);


} // end of read_simsed



// =====================================
int ISED_MATCH(double *parval) {

  double sumdif, sumdif_min, dif ;
  int ised_match, ised, ipar ;

  char fnam[] = "ISED_MATCH" ;

  // ------------- BEGIN -------------

  // now find SED index that matches NEAREST_GRIDVAL
  sumdif_min = 999999. ;
  ised_match = 0 ;

  for ( ised = 1; ised <= SEDMODEL.NSURFACE; ised++ ) {

    sumdif = 0.0 ;
    for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {
      if ( INPUTS.USEPAR[ipar] == 0 ) { continue ; }
      dif = parval[ipar] - SEDMODEL.PARVAL[ised][ipar] ;
      sumdif += fabs(dif);
    }

    if ( sumdif < sumdif_min ) {
      sumdif_min = sumdif;
      ised_match = ised;
    }

  } // ised

  if ( ised_match < 1 || ised_match > SEDMODEL.NSURFACE ) {
    sprintf(c1err,"Invalid  ised_match = %d ", ised_match );
    errmsg(SEV_FATAL, 0, fnam, c1err, "" );       
  }


  return ised_match ;

} // end of ISED_MATCH

// ****************************
void getSpec(int ep) {

  // get spectrum corresponding to rest-frame epoch "Trest".
  //
  // Fill the following global arrays that depend on the
  // lambda 'ilam' bin
  //    SPECLAM[ilam] 
  //    SPECFLUX[ilam]
  //    SPECFLUXERR[ilam]
  //
  // Fill the global SPEC_LAMSTEP variable (assumed constant)
  //
  // Apr 3, 2011 (RK)
  //  replace interp8() with interp_1DFUN() and fix bug passing
  //  arguments ... not sure how it worked before ??
  //  Previous args to interp8 should have been &a_day[1] instead
  //  of a_day, etc ...
  //
  //  Finally do the SED-parameter interpolation using the
  //  corners defined in interp_prep().
  //
  // -----------------------

  int 
     iday_match, iday_mid, iday
    ,ilam, jflux, jjflux, k, kk, icorner
    ,RDFLAG, ISED
    ;

  double 
    Trest
    ,tdif_min, dif
    ,flux, fluxerr, lam
    ,a_day[3], a_flux[3], a_fluxerr[3]
    ,Fscale, Ftmp, WGT
    ;

  char fnam[] = "getSpec" ;
  int OPT_INTERP  = 1 ; // 1=linear, 2=quadratic
  int NBIN_INTERP = 3 ;

  // ------------- BEGIN ------------

  if ( ep < 0 ) { 
    Trest = 0.0;
  } else {
    Trest = INPUTS.EPOCH[ep] ; // strip off Trest
  }

  // now find closest day
  tdif_min = 99999. ;
  iday_match = -9;
  for ( iday=0; iday < TEMP_SEDMODEL.NDAY; iday++ ) {
    dif = fabs(Trest - TEMP_SEDMODEL.DAY[iday]);
    if ( dif < tdif_min ) {
      iday_match = iday ;
      tdif_min   = dif ;
    }

  } // iday

  //check location of iday_match within iday range
  if ( iday_match == 0 ) {
    iday_mid = iday_match + 1;
  } else if ( iday_match == TEMP_SEDMODEL.NDAY - 1 ) {
    iday_mid = iday_match - 1;
  } else {
    iday_mid = iday_match; 
  }// iday_mid


  // prepare for interpolation to input rest frame epoch
  NBIN_LAM = TEMP_SEDMODEL.NLAM; 

  // same a_day for each lambda, so set up once, then dummy check
  a_day[0] = TEMP_SEDMODEL.DAY[iday_mid-1];
  a_day[1] = TEMP_SEDMODEL.DAY[iday_mid];
  a_day[2] = TEMP_SEDMODEL.DAY[iday_mid+1];

  // dummy check
  if ( Trest < a_day[0] || Trest > a_day[2] ){
    sprintf(c1err,"EPOCH=%3.1f is outside expected range of", Trest );
    sprintf(c2err,"%3.1f - %3.1f epochs", a_day[0], a_day[2] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }


    // finally, lambda loop to interpolate flux, fluxerror at proper day, 
    // load into global arrays

    Fscale = SEDMODEL.FLUXSCALE ;
    for ( ilam=0; ilam < NBIN_LAM; ilam++ ) {
      lam      = TEMP_SEDMODEL.LAM[ilam] ;
    
      jflux = NBIN_LAM*(iday_mid) + ilam; 

      for ( k=0; k <= 2; k++ ) 
	{ a_flux[k] = a_fluxerr[k] = 0.0 ; }

      // interpolate in SED-parameter space
      for ( icorner=0; icorner < CORNER.NSTORE; icorner++ ) {

	WGT    = CORNER.WGT[icorner] / CORNER.WGTSUM ;
	ISED   = CORNER.ISED[icorner] ;

	for ( k=0; k <= 2; k++ ) {
	  kk       =  k - 1 ;  // -1, 0, +1
	  jjflux   = jflux + kk*NBIN_LAM ; 
	  flux     = *(CORNER.FLUX[icorner]    + jjflux  );
	  fluxerr  = *(CORNER.FLUXERR[icorner] + jjflux  );
	  a_flux[k]    +=  WGT * flux ;
	  a_fluxerr[k] +=  WGT * fluxerr ;
	}
      }

      // now interpolate the epoch
      flux    = interp_1DFUN (OPT_INTERP, Trest, 
			      NBIN_INTERP, a_day, a_flux,    fnam );
      fluxerr = interp_1DFUN (OPT_INTERP, Trest, 
			      NBIN_INTERP, a_day, a_fluxerr, fnam );

      SPECLAM[ilam]     = lam ;
      SPECFLUX[ilam]    = Fscale * flux ;
      SPECFLUXERR[ilam] = Fscale * fluxerr ;

    }// end of interpolation and global array load  


   SPEC_LAMSTEP = TEMP_SEDMODEL.LAMSTEP;


} // end of getSpec


// ************************************
void obsSpec(void) {

  // Created Feb 2011 by R.Kessler
  // Transform spectrum to observer-frame using inputs
  // INPUTS.REDSHIFT (if > 0) and distance modulus INPUTS.DLMAG.
  // Note that SPECLAM, SPECFLUX and SPECFLUXERR are 
  // all input & output arrays.  SPEC_LAMSTEP is also transformed
  // to observer frame. 

  double 
    LAM
    ,SPECLAM_REST[MXBIN_LAMSED_SEDMODEL]
    ,SPECFLUX_REST[MXBIN_LAMSED_SEDMODEL]
    ,SPECFLUXERR_REST[MXBIN_LAMSED_SEDMODEL]
    ,z, mu, snr
    ,z1, zfactor, DLratio, SQDLratio
    ;

  int ilam ;

  // ---------- BEGIN ------------

  if ( INPUTS.REDSHIFT <= 0.0 ) { return ; }

  // store rest-frame spectrum & error in local array

  for ( ilam = 0; ilam <  NBIN_LAM; ilam++ ) {    
    SPECLAM_REST[ilam]     = SPECLAM[ilam] ;
    SPECFLUX_REST[ilam]    = SPECFLUX[ilam] ;
    SPECFLUXERR_REST[ilam] = SPECFLUXERR[ilam] ;
  }

  z   = INPUTS.REDSHIFT ;
  mu  = INPUTS.DLMAG ;
  snr = INPUTS.SNR ;

  z1      = 1.0 + z;
  zfactor = 1. / z1 ;
  SQDLratio = pow(10.,-0.4*mu);
  
  // apply distance and cosmological redshift, move to output arrays
  for ( ilam = 0; ilam < NBIN_LAM; ilam++ ) {
    SPECLAM[ilam]     = SPECLAM_REST[ilam]*z1 ;
    SPECFLUX[ilam]    = SPECFLUX_REST[ilam]   * SQDLratio * zfactor;
    SPECFLUXERR[ilam] = SPECFLUXERR_REST[ilam]* SQDLratio * zfactor; 
  }

  SPEC_LAMSTEP = SPEC_LAMSTEP*z1; 

} // end of obsSpec


// ************************************
void fudgeSNR(int ep) {

  // Created Feb 2011 by J. Mosher
  // Include fudged fluxerror using input
  // INPUTS.SNR ( if > -101 ) . 
  // NOTE that SPECFLUX and SPECFLUXERR are 
  // input, output arrays. 

  double 
    snr, snrScale, snrlambin, snrSpec 
    ,lam_binsize
    ,ranGauss, ranFlux
    ,F, FERR, LAM, Z1
    ; 

  int ilam, iSNR, whereLAM, whereSNR;

  char fnam[] = "fudgeSNR" ;

  // ---------- BEGIN ------------

  whereSNR = 0;
  whereLAM = 0;
  LAM = SPECLAM[0];

  init_getSNR(&whereSNR,&whereLAM,LAM, ep);
  //printf ("LAM : %f , whereSNR gt -99.9 %d, whereLAM fits %d \n", LAM, whereSNR, whereLAM);
  snr         = INPUTS.specSNR[ep][whereSNR]; //check if user input

  if ( snr <= -100.0 ) { return ; }

  snrlambin   = INPUTS.specSNRLAMBIN[ep];  // bin size to define SNR
  lam_binsize = SPEC_LAMSTEP;     // bin size of input spectrum
  fill_RANLISTs();

  // scale SNR to binsize of spectrum
  snrScale = sqrt(lam_binsize/snrlambin);

  // if ideal snr is desired, re-calculate SPECFLUXERR

  for (ilam  = 0; ilam < NBIN_LAM; ilam++ ){
    LAM      = SPECLAM[ilam];
    snr      = getSNR(&whereSNR, &whereLAM, LAM, ep);
    snrSpec  = snr * snrScale ;    // spectrum SNR in spectrum bins
    F        = SPECFLUX[ilam] ;
    if ( F == 0.0 ) {
      if (strcmp(INPUTS.ZFLUX,"ZFLUX") == 0.0) {
	F = ZFLUXVAL;
      }
    }
    FERR     = fabs(F/snrSpec) ; //FERR always positive
    ranFlux  = SPECFLUX[ilam] + FERR * getRan_Gauss(1) ;

    SPECFLUX[ilam]    = ranFlux;                // measured flux
    SPECFLUXERR[ilam] = FERR * sqrt(fabs(ranFlux/F)); // measured error

    // DDDDDDDDDDDDD
    if ( fabs(LAM-5000) < 3.0 ) {
      printf(" xxxx LAM_OBS=%8.2f  FERR=%le   snr=%6.2f  ZScale=%le \n",
	     LAM, FERR, snr, snrScale );
    }
    // DDDDDDDDDDDDD

  }


} // end of fudgeSNR


// ************************************
void wrSpec(int ep) {

  FILE *fp;
  int ilam, ERRFLAG, SNANAFLAG, istat ;
  double Laminit, Lamfinal, LAMDIF ;
  char 
    ctmp[80]
    ,ctmperr[80] 
    ,ctmpheader[80]
    ,line[80]
    ,specFile[200]
    ,fnam[] = "wrSpec"
    ;

  // -------------- BEGIN ---------------

  sprintf(specFile, "%s" , "NULL");

  if ( strcmp(INPUTS.specFile[ep],"NULL") != 0 ) {
    sprintf(specFile, "%s", INPUTS.specFile[ep] );
    // if file exists, then append; else create new file.
    istat = access(specFile, F_OK);
  }
  else {
    istat = -1;
    if ( INPUTS.NEPOCH == 1 ) {
      sprintf(specFile, "%s", "SIMSED.SPEC" );
    }
    else {
      sprintf(specFile,"SIMSED_%s.SPEC", INPUTS.char_EPOCH[ep] );
    }
  }



  Laminit  = INPUTS.specLAMinit_obs[ep];
  Lamfinal = INPUTS.specLAMfinal_obs[ep];

  if ( istat == 0 )  { 
    fp = fopen(specFile, "at") ;  // append
    sprintf(ctmp,"%s" , "Append" );
  }
  else {
    fp = fopen(specFile, "wt") ;    // new file
    sprintf(ctmp,"%s" , "Write" );
  }


  printf(" %s extracted spectrum to : %s \n", ctmp, specFile );

  ERRFLAG = (FLUX_ERRFLAG  || INPUTS.SNR > 1.0E-9 ) ;
  SNANAFLAG = (strcmp(INPUTS.TYPE,"SNANA") == 0) ;


  for ( ilam = 0; ilam <  NBIN_LAM; ilam++ ) {

    if ( SPECLAM[ilam] <= Lamfinal && SPECLAM[ilam] >= Laminit ) {
      
      sprintf(ctmp," %8.2f  %le  ", SPECLAM[ilam], SPECFLUX[ilam] );
      sprintf(ctmperr," ");
      sprintf(ctmpheader," ");
      
      if ( ERRFLAG ) {
	sprintf(ctmperr,"%le", SPECFLUXERR[ilam] );
      }
      
      if ( SNANAFLAG ) {
	sprintf(ctmpheader,"SPEC%02d:", ep); 
      }
      

      sprintf(line,"%s %s %s", ctmpheader, ctmp, ctmperr );
      fprintf(fp,"%s\n", line );

      // check dump option
      LAMDIF = (SPECLAM[ilam] - INPUTS.LAMDUMP);
      if ( fabs(LAMDIF) < 2.0 ) {
	printf("%s : %s \n", specFile, line);
      }

    }
  }

    if ( SNANAFLAG )  {
       fprintf(fp,"END_SPEC%02d:\n\n", ep);
    }


  fclose(fp);

} // end of wrSpec

// ********************************
void write_footer(int ep) {

  // read Trest and optional spectrum filename on each line op *epFile.

  FILE *fp;
  char
    *ptrtok
    ,line[200]
    ,tmpline[200]
    ,strTrest[100]
    ,specFile[200]
    ,ctmp[80]
    ;

  int
    istat
    ,N
    ,ERRFLAG
    ;
  double Trest ;

  // ------------- BEGIN --------------


  if ( strcmp(INPUTS.TYPE,"SNANA") != 0 || ep < INPUTS.NEPOCH ) { return; }

  if ( strcmp(INPUTS.specFile[ep],"NULL") != 0 ) {
    sprintf(specFile, "%s", INPUTS.specFile[ep] );
    // if file exists, then append; else create new file.
    istat = access(specFile, F_OK);
  }
  else {
    istat = -1;
    if ( INPUTS.NEPOCH == 1 ) {
	sprintf(specFile, "%s", "SIMSED.SPEC" );
      }
      else {
	sprintf(specFile,"SIMSED_%s.SPEC", INPUTS.char_EPOCH[ep] );
      }
    }

    if ( istat == 0 ) {
      fp = fopen(specFile, "at") ;  // append
      sprintf(ctmp,"%s" , "Append" );
    }
    else {
      fp = fopen(specFile, "wt") ;    // new file
      sprintf(ctmp,"%s" , "Write" );
    }


    printf(" %s footer to : %s \n", ctmp, specFile );

    ERRFLAG = (FLUX_ERRFLAG  || INPUTS.SNR > 1.0E-9 ) ;

    fprintf(fp, "#END:\n");

fclose(fp);

} // end of write_footer

// ********************************
void init_getSNR(int *wsnr, int *wlam, double lam, int ep) {

  // given an inital (obs) wavelength lam and epoch ep,
  // this subroutine finds the appropriate first SNR wavelength
  // index wlam and first non-zero SNR value index wsnr

  int iSNR, i;

  iSNR = 0;
  while(lam > SNRLAM[*wlam+1] && iSNR < NXSNR-1 ) {
    *wlam = *wlam + 1;
    iSNR = iSNR + 1;
  }

  iSNR = 0;
  while(INPUTS.specSNR[ep][iSNR] < 0 && iSNR < NXSNR-1) {
    *wsnr = *wsnr + 1;
    iSNR = iSNR + 1;
  }

} // end of init_getSNR

// ********************************
double getSNR(int *wsnr, int *wlam, double lam, int ep) {

  // for a given (obs) wavelength lam and epoch ep,
  // this subroutine uses the snr vector index *wsnr and 
  // wavelength vector index *wlam to find the appropriate
  // (obs) snr to be used. If index *wlam  is outside range of valid
  // user input snrs, the first (or last) valid snr is returned. 
  // Otherwise, snr is linearly interpolated to the input lam. 

  // JLM_ARGS: please ABORT if lam is invalid 

  int    iSNR, i;
  double slope, intercept;
  double snr, DIF_SNR, DIF_LAM ;

  // check *wlam 
  if (lam > SNRLAM[*wlam+1]  && *wlam < NXSNR-1 ) {
    *wlam = *wlam + 1;
  }

  // check *wsnr
  if (*wsnr < *wlam && INPUTS.specSNR[ep][*wsnr+1] > 0 && *wsnr < NXSNR-1){
    *wsnr = *wsnr + 1;
  }

  if (*wlam < *wsnr) {  
    snr = INPUTS.specSNR[ep][*wsnr]; 
  } 
  else if (*wlam == *wsnr && INPUTS.specSNR[ep][*wsnr+1] > 0) { 

    DIF_SNR = INPUTS.specSNR[ep][*wsnr+1] - INPUTS.specSNR[ep][*wsnr] ;
    DIF_LAM = SNRLAM[*wsnr+1]-SNRLAM[*wsnr] ;
    slope   = DIF_SNR / DIF_LAM ;
    intercept = INPUTS.specSNR[ep][*wsnr];
    snr = slope*(lam-SNRLAM[*wsnr]) + intercept;
  } 
  else { 
    snr = INPUTS.specSNR[ep][*wsnr]; 
  }

  return snr;

} // end of getSNR



