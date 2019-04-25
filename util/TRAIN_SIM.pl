#!/usr/bin/perl
#
# --------------
# revamp of SALT2train_prep.pl -JLM
# tips from SALT2train_prep-comments.README (RK)
#
# Usage:
#
# ./SALT2train_prep.pl  <input file>
#
# This code was written to make simulated LC + SPECTRUM training sets. 
# Given a light curve simulation input file, which must includes the 
#  necessary SIMSED arguments, script will generate LC+SPECTRUM 
# data files in $SNDATA_ROOT/SIM. 
#
# Typical input file syntax is as follows:
# 
# XX ORIGINAL "IDEAL" SPECTRUM DISTRIBUTION
# LCSIM_INPUT: SIMSED_SALT2_6.input GENVERSION TESTpcafit7 GENMAG_SMEAR .1 NGEN_LC 100
# SPECSIM_ARGUMENTS: IDEAL_EPOCHRANGE -15 40 IDEAL_CADENCE 10
# END:
# 
# The LCSIM_INPUT: keyword should be followed by a simulation inputfile, a GENVERSION name
# (here TESTpcafit7 ), and any alterations to the simulation input file. These use the 
# snana simulation input file keywords. The last line, beginning with keyword
# SPECSIM_ARGUMENTS: is how alterations to spec sim defaults are made. 
#
# Spectrum simulation [spec sim] defaults are written into this script. 
# To see them, 
# > grep DEFAULT SALT2train_sim.pl .
#
# If you wish to change spec sim defaults, input file keywords are as follows:
# IDEAL_EPOCHRANGE -15 50 [rest frame]
# IDEAL_SNR_OBS 20 100 [20 is SNR, 100 is binsize for SNR] 
# IDEAL_CADENCE 5 [ first spectrum epoch chosen randomly w/in cadence of first lc ep ]
# IDEAL_LAMINIT_REST 2000 [ rest frame ]
# IDEAL_LAMFINAL_REST 9200 [ rest frame ]
#  
# ***THESE SHOULD BE ENTERED IN THE INPUT FILE, BEHIND THE KEYWORD
# SPECSIM_ARGUMENTS: *** 
# Multiple SPECSIM_ARGUMENTS: lines may be included for any one
# LCSIM: keyword call.
# 
# To run a spectrum-library based simulation, provide the path to the library file on the
# SPECSIM_ARGUMENTS: line. Here's an example:
# SPECSIM_ARGUMENTS: SPECGEN_LIB /home/s1/jmosher/RICK-simprep/REALSPEC_SIM/SALT2-trainingmap/temp.out .
#
#
#
# 
###################
### For easy file navigation: sections are as follows:
#
### DECLARATIONS
### OUTPUT FILE FORMATS 
### MAIN CODE
### subroutines (listed individually)
#
##### TO DO LIST ##
# to do list:
#
#   * tarball training set
#   * use SIMSED epoch range as default
#   * keyword for alternate training set location path
#
#
##### HISTORY #####
#
# Feb 28 2011: functional - can make sims. see to do list. 
#
# Aug 18 2011: updated init_INPUTKEYS to allow multiple pipeline END: calls
#              JLM
#
# Aug 23 2011: updated runLib_SPECSIM to pass each library spectrum's 
#              FIRSTWAVE, LASTWAVE data to SIMSED_extractSpec.exe .
#              JLM
#
# Aug 24 2011: updated input KEYWORDS to allow IDEAL_SNR_OBS <SNR><binsize> input option
#              default values are 20 100[Angstroms]. 
#
# Aug 25 2011: RK  for SIMSED_extractSpec args, remove --dlam arg and combine
#                  snr and snrlamrange args into --snr <snr> <snrlambin>
#
# Aug 25 2011: JLM fixed bug in IDEAL sim SIMSED_extractSpec input file creation.
#                  File now conforms to <epoch><filename><laminit><lamfinal> convention. 
#                  IDEAL_SNR_OBS input is now two parameters IDEAL_SNR_OBS <snr> <snrlambin>
#
# Sep 14 2011: JLM added get_RANSEED subroutine to use same RANSEED for LC, initial spec gen. 
#                  added make_SUMMARY subroutine to dump sim summary for each GENVERSION. 
#                  fixed some small bugs in distribution-based spec default param allocation. 
#                  defined $SNDATA_ROOT and $SNANA_DIR in main
#                  wrote wrapper for SIMSED_extractSpec call
#                  
# Oct 17 2011: JLM added option to choose SNR epoch by epoch for simulated spectra. 
#                  Writes extra lines to SIMSED_extractSpec input file 
#                  New format is <epoch><filename><laminit><lamfinal><snrlambin><snr2000><...><snr10000>.
#                  This change only affects LIBRARY-based simulations. 
#
# Oct 18 2011: JLM fixed bug that was allowing build up of junk files.
#            
# Oct 20 2011: JLM replaced Perl FORMAT calls with sprintf calls. 
#
# Nov 29 2011: JLM library now matches only by redshift ( rather than by z,c,x1 ).
#                  added necessary components for galaxy cantamination;
#                  these use SPECSIM_ARGUMENTS: keywords GALFRAC and GALFILE.
#                  updated SIMSED model param_extract to account for dummy parameters
#                  (i.e. Hsiao SIMSED models ). 
#            
# Dec 12 2011: JLM updated call_SIMSEDextractSpec()  to use OPT_ZEROFLUX option
#                     
# Dec 15 2011: JLM checkparam_SPECSIM() now searches SIM .README to determine real SIMSED params  
#
# Jan 4 2011:  JLM Updated wavelength, SNR variable names to include frame.
#                  Input keywords have been amended as follows. 
#                  IDEAL_LAMINIT -> IDEAL_LAMINIT_REST
#                  IDEAL_LAMFINAL -> IDEAL_LAMFINAL_REST
#                  IDEAL_SNR -> IDEAL_SNR_OBS
#                  (although old keywords should also still work). 
#                   
# May 11 2012: JLM Added crazy input SIMSED parameter check
#                  Code aborts if param vals are > 100 or < -100
#                  SIMSED parameter parse compatible with SNANAv9.95
#                  
#
# Jun 01 2013: JLM Made small changes to SIMSED_extractSpec cleanup to allow TRAIN_SIM.pl to be 
#                  run in batch mode. 
#
###################
# --------------

use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use sntools;

use File::Basename;
use File::Copy;
use POSIX;
use Cwd;
use strict;

sub make_banner($);
sub make_subbanner($);
sub jobinfo_dump(@);
sub check_fileexist(@);
sub check_goodexit(@);
sub var_assign(@);		   
sub readCOLUMN(@);		   
sub readCOLUMN_CUT(@);		   
sub readCOLUMN_LIBCUT(@);		   
sub readCOLUMN_EQCUT(@);		   
#
sub get_inFile(@);
sub parse_inFile_SALT2train($);
sub parse_SPECDIST(@);			    
sub init_INPUTKEYS();
sub parse_LCSIMLINE($);
sub init_RUNJOBS($);
sub init_LCSIM($);
sub run_LCSIM();
#
sub init_SPECSIM($);
sub run_SPECSIM();
sub prep_SPECLOC($);
sub get_RANSEED($);
sub increment_RANSEED();
sub undefHASH;
sub readLIB($);
sub readSIM($);
sub match_LIB_SIM(@);
sub runLib_SPECSIM();
sub runIdeal_SPECSIM();
sub split_param(@);
sub getparam_SPECSIM($); ### reads LCSIM dump to get SIMSED parameters
sub checkparam_SPECSIM(); ### reads SED.INFO to evaluate SIMSED parameters
sub genspec_GAMMA_SPECSIM(@);
sub phase_generator($);
sub call_SIMSEDextractSpec(@);
#
sub initvar_SPECSIM();	### sets all dist parvals to zero. needed?
sub fill_EPOCHvar_SPECSIM(@); ### allocates dist parvals to hashes
sub read_EPOCHvar_SPECSIM($); ### reads dist parvals from hashes
sub fill_PHASEvar_SPECSIM(@); ### allocates dist parvals to hashes
sub read_PHASEvar_SPECSIM($); ### reads dist parvals from hashes
sub fillvar_SPECARG(@); ### parses fill SPECARG line into hashes
sub readvar_SPECARG($);	### reads SPECARGS from hashes	    
#
sub gaussian_rand(@);
sub gaussian_assign(@);
sub poisson_rand($);
sub poisson_assign(@);
sub exp_rand($);
sub exp_assign(@);
sub gamma_rand(@);
sub gamma_assign(@);
#
sub make_SUMMARY;
### DECLARATIONS ###

my $LOCAL_DIR; ### current working directory
my $SNDATA_ROOT;
my $SNANA_DIR;
my $USER_inFile; ### User input file name
my $ABORT; ### Use for successful run check
###
my $RUN_LCSIM; ### if T, run lc simulation
my $RUN_SPECSIM; ### if T, run spec simulation
my $RUN_NORMTEST; ### if T, run normalization check
###
#PARSE VARIABLES
my $OPTABORT = 0;
my $OPTWARN = 1;
my $OPTQUIET = 2;
my $NARG1 = 1;
my $NARG2 = 2;
my $NARG7 = 7;
my $NARGALL = 99;
my @TMP; ### consistent with rest of snana sntools.pl usage
my $KEY; 
###
my %SALT2MAGSYS_HASH; ###
my %LCSIM_HASH; ###
my %SIMMINZ_HASH;
my %SIMMAXZ_HASH;
my %SPECSIM_HASH; ###
my %GENLIB_HASH;
my %GALFRAC_HASH;
my %GALFILE_HASH;
# HASHES FOR DISTRIBUTIONS
my %EPOCHDIST_HASH;
my %PHASEDIST_HASH;
my %EPOCH_DPAR1_HASH;
my %EPOCH_DPAR2_HASH;
my %EPOCH_DPAR3_HASH;
my %EPOCH_DPAR4_HASH;
my %EPOCH_DPAR5_HASH;
my %EPOCH_DPAR6_HASH;
my %PHASE_DPAR1_HASH;
my %PHASE_DPAR2_HASH;
my %PHASE_DPAR3_HASH;
my %PHASE_DPAR4_HASH;
my %PHASE_DPAR5_HASH;
my %PHASE_DPAR6_HASH;
# HASHES FOR LIBRARY
my %LIB_c_HASH;
my %LIB_x1_HASH;
my %LIB_z_HASH;
my %LIB_nspec_HASH;
my %LIB_ispec_HASH;
my %LIB_t0_HASH;
my %LIB_usedval_HASH;
my %SIM_z_HASH;
my %SIM_S2x1_HASH;
my %SIM_S2c_HASH;
my %SIM_matchval_HASH;
my %SIM_bestmatch_HASH;
my %LIB_specmjd_HASH;
my %LIB_firstwave_HASH;
my %LIB_lastwave_HASH;
my %LIB_SNR2000_HASH;
my %LIB_SNR9000_HASH;
my %LIB_SNR10000_HASH;
my %LIB_SNR3000_HASH;
my %LIB_SNR4000_HASH;
my %LIB_SNR5000_HASH;
my %LIB_SNR6000_HASH;
my %LIB_SNR7000_HASH;
my %LIB_SNR8000_HASH;
my %LIB_SNRNBIN_HASH;
my %LIB_BINWAVE_HASH;
#
my %DLAM_HASH;
my %EPOCHRANGE_LOW_HASH;
my %EPOCHRANGE_HIGH_HASH;
my %SNR_HASH;
my %SNRLAMBIN_HASH;
my %CADENCE_HASH;
my %LAMINIT_HASH;
my %LAMFINAL_HASH;
###
my $LCSIM_INPUT; ### User sim input file name
my $LCSIM_CALL; ### Full set of command line vars for lc_sim.exe
my $SIM_ZMIN;	    
my $SIM_ZMAX;	    
my $INPUT_FILE_INCLUDE; ### sim SIMSED INCLUDE file name
my $JOBNAME_SNLCSIM; ### path to snlc_sim.exe binary
my $SNLC_SIM_LOG; ### log file to use for sim
my $GENVERSION; ### sim input keyword value
###
my $LCDUMP_FILE; ### lc sim file from which to grab SNe cids
my $SPECDUMP_FILE;
my $SIMREADME; ### lc sim README file
my $SPECSIM_dir; ### dir in which to put spec sims
my $UNPACK_dir; ### dir in which to unpack SALT2 formatted sims
###
my $SPECGEN_LIB; ###(DEFAULT NONE) path to spectrum sim gen library
my @IDEAL_SPECARGS; 
my @IDEAL_EPOCHRANGE; ###(DEFAULT -15 50 ) SIMSED rfepoch range
my @IDEAL_SNR_OBS; ### (DEFAULT 20 100 ) SIMSED SNR and SNR binrange 
my $IDEAL_CADENCE; ### (DEFAULT  5 ) SIMSED cadence
my $IDEAL_LAMINIT_REST; ###  (DEFAULT 2000 ) initial SIMSED rest frame wavelength
my $IDEAL_LAMFINAL_REST; ### (DEFAULT 9200 ) final SIMSED rest frame wavelength
my $IDEAL_DLAM; ### ( DEFAULT 0.0 ) SIMSED observer frame dlam to use		 
my $GENMODEL; ### SIMSED model to use
my $RANSEED; ## LCSIM SEED TO PASS TO SPECSIM
my ($GALFRAC, $GALFILE); ##INFO ON GALAXY CONTAMINATION TO PASS TO SPECSIM	    
my (@USEPARVALS, @PARNAMES);
my @SN_SIM_PARAMS; ### HOLDING LCSIM PARAMS	     
my $JOBNAME_SIMSEDextract; ### SED extract script to use
my $OUTPUT_FORMAT; ### extraction type ( SNANA or ASCII )
my $TRAINLIST_file; ### subroutine generates rfphase list of SEDS to make
###		
# DECL FOR SPEC DISTRIB GEN
my $EPOCH_DISTRIBUTION;
my $PHASE_DISTRIBUTION;
my @NEPOCH_LIST;
my @VAR_ARRAY;
my @PHASEVAR_ARRAY;
my @EPOCHVAR_ARRAY;
my $DIST_PAR1;
my $DIST_PAR2;
my $DIST_PAR3;
my $DIST_PAR4;
my $DIST_PAR5;
my $DIST_PAR6;
my $CADENCE;
# put these into file 
my $DEFAULT_FRAC1 = 0.85;
my $DEFAULT_FRAC0 = 0.10;
my $DEFAULT_FRAC2 = 0.05;	     
my $DEFAULT_POISSON = 1;
my $DEFAULT_EXPLAM = 0.026;
my $DEFAULT_EXPOFFSET = -35.63;
my $DEFAULT_GAUSSMEAN = 4;
my $DEFAULT_GAUSSSTDEV = 5;
my $DEFAULT_GAUSSCADENCE = 1;
my $DEFAULT_GK1 = 25;
my $DEFAULT_GK2 = 15;
my $DEFAULT_GTHETA1 = 2.;
my $DEFAULT_GTHETA2 = 2.;
my $RANSEED_OFFSET = 833; ## CHANGES SEED BY FIXED AMOUNT EACH SPECTRUM GEN
# available types: poisson, gamma, gaussian
### 
my $printprefix; ### all print names are for perl FORMAT calls
my $printphase;
my $printdate;
my $printz;
my $printCID;
my $printtype;
my $printfile;
my $printlaminit_rest;
my $printlamfinal_rest;
my $printspecfile;
my ($printSNR1, $printSNR2, $printSNR3, $printSNR4);
my ($printSNR5, $printSNR6, $printSNR7, $printSNR8, $printSNR9);
my $printSNRLAMBIN; 


### OUTPUT FILE FORMATS ###

format TRAINLIST = 
@* @* @*
$printCID $printtype $printfile
. 


format SPECHEAD =
@FORMAT STD
"@"
@Date @####.##
"@", $printdate
@Redshift @.####
"@", $printz
@SN @>>>>>>>>
"@", $printCID
@Source SIMU
"@"
@WAVE : wavelength, observer frame, in A
"#"
@FLUX : erg/s/cm2/A
"#"
@ERR  : statistical uncertainty on flux
"#"
@end
"#"
.

###  MAIN CODE ###

$LOCAL_DIR = getcwd; 
$SNDATA_ROOT = $ENV{'SNDATA_ROOT'}; 
$SNANA_DIR = $ENV{'SNANA_DIR'};

$USER_inFile = get_inFile(@ARGV);

init_RUNJOBS($USER_inFile);

parse_inFile_SALT2train($USER_inFile);


if ($RUN_LCSIM eq 'T' ) {
    init_LCSIM($USER_inFile);
    run_LCSIM();
}

if ($RUN_SPECSIM eq 'T' ) {
    init_SPECSIM($USER_inFile);
    run_SPECSIM();
}

make_SUMMARY(); 

print " FINISHED SUCCESSFULLY \n ";

###  END  ###
# ----------

sub get_inFile(@)
{
    
    my $NARG = scalar(@_);
    my $tempFile; 
    
    ## check that user gave input filename
    if ( $NARG < 1 )  {
    print " Must give input filename as argument. \n" ;
    die " ***** ABORT ***** \n" ;
}
    
    $tempFile = $ARGV[0];

    check_fileexist("input", $tempFile);

    return($tempFile); 	

}

    
sub init_RUNJOBS($)
{
    my ($inFile) = shift;

    ### CHECK INPUT FILE FOR USER RUNJOBS
    make_banner("CHECKING FOR USER-SELECTED RUNJOBS") ; 
    print "USER RUNJOBS inputs are: \n";

    $KEY = "RUN_LCSIM:" ; ## OPTIONAL
    @TMP = sntools::parse_line($inFile, $NARG1, $KEY, $OPTQUIET) ;
    $RUN_LCSIM = @TMP[0];
    unless ( $RUN_LCSIM eq '' ) { print "\t RUN_LCSIM \t $RUN_LCSIM \n"; }

    $KEY = "RUN_SPECSIM:" ; ## OPTIONAL
    @TMP = sntools::parse_line($inFile, $NARG1, $KEY, $OPTQUIET) ;
    $RUN_SPECSIM = @TMP[0];
    unless ( $RUN_SPECSIM eq '' ) { print "\t RUN_SPECSIM \t $RUN_SPECSIM \n"; }

    $KEY = "RUN_NORMTEST:" ; ## OPTIONAL
    @TMP = sntools::parse_line($inFile, $NARG1, $KEY, $OPTQUIET) ;
    $RUN_NORMTEST = @TMP[0]; 
    unless ( $RUN_NORMTEST eq '' ) { print "\t RUN_NORMTEST \t $RUN_NORMTEST \n"; }

    $KEY = "OUTPUT_FORMAT:" ; ## OPTIONAL
    @TMP = sntools::parse_line($inFile, $NARG1, $KEY, $OPTQUIET) ;
    $OUTPUT_FORMAT = @TMP[0];
    if ( $OUTPUT_FORMAT eq '' ) {
	$OUTPUT_FORMAT = "SNANA";
    }
    print "\nOUTPUT_FORMAT is $OUTPUT_FORMAT \n\n";


    ### LOAD DEFAULT RUNJOBS
    if ( $RUN_LCSIM eq '' ) { 
	$RUN_LCSIM = "T";
	print "USING DEFAULT RUNJOB RUN_LCSIM $RUN_LCSIM\n";
	}
    if ( $RUN_SPECSIM eq '' ) { 
	$RUN_SPECSIM = "T";
	print "USING DEFAULT RUNJOB RUN_SPECSIM $RUN_SPECSIM\n";
	}
    if ( $RUN_NORMTEST eq '' ) { 
	$RUN_NORMTEST = "F";
	print "USING DEFAULT RUNJOB RUN_NORMTEST $RUN_NORMTEST\n";
	}

}

sub init_LCSIM($)
{

### This subroutine called if RUN_LCSIM = T
### Subroutine checks for sim.input file
###            establishes needed variables
###            by parsing input, sim_input files

    my ($inFile) = shift;

    ## extract lcsim job name from user input file

    make_banner("PREPARING FOR LCSIM");

    $KEY = "JOBNAME_SNLCSIM:" ; ## OPTIONAL
    @TMP = sntools::parse_line($inFile, $NARG1, $KEY, $OPTQUIET) ;
    $JOBNAME_SNLCSIM = @TMP[0];
    unless ( $JOBNAME_SNLCSIM eq '' ) {
	print "\t JOBNAME_SNLCSIM \t $JOBNAME_SNLCSIM \n"; }

    ## load default information
    if ( $JOBNAME_SNLCSIM eq '' ) {
	$JOBNAME_SNLCSIM = $SNANA_DIR . "/bin" . "/snlc_sim.exe" ;
    }

}	    

sub run_LCSIM()
{

    my @output; #use to dump any unix returns from job call
    my $ABORT_die; #use if ABORT found in log file
    my $GEN_SALT2INSTRUMENT;
    my $SPECSYM_i;

    make_banner("STARTING LIGHT CURVE SIMULATION");

    ### FOR EACH APPROPRIATE LCSIM ENTRY IN HASH, LCSIM JOB IS RUN ###

    ## For each instrument, magsys pair
	while ( ($GENVERSION, $LCSIM_CALL) = each(%LCSIM_HASH)) {
		## Run LCSIM job
                unless ( $LCSIM_CALL eq "NULL" ) {

		    make_banner("RUNNING LCSIM FOR GENVERSION: $GENVERSION");

		    $SNLC_SIM_LOG = $LOCAL_DIR . '/' . 'snlcsim' . '_' . $GENVERSION . '.log';

		    jobinfo_dump($JOBNAME_SNLCSIM, $LCSIM_CALL, $SNLC_SIM_LOG);

		    @output = `$JOBNAME_SNLCSIM $LCSIM_CALL >& $SNLC_SIM_LOG `;

		    ### check log file for ABORT - if ABORT, then die with error

		    print "\t\t CHECKING LCSIM OUTPUT \n\n";

		    $KEY = "ABORT" ; ## OPTIONAL

		    check_fileexist("log file", $SNLC_SIM_LOG);

		    check_goodexit("SNLCSIM", $SNLC_SIM_LOG, $KEY);

		    print "\t\t FINISHED WITH LIGHT CURVE SIM FOR $GENVERSION \n\n"; 
   
		}
            }
}
    

sub run_SPECSIM()
{

    make_banner("RUNNING SPECTRUM SIMULATIONS");

    ### FOR EACH APPROPRIATE SPECSIM ENTRY IN HASH, SPECSIM JOB IS RUN ###

	while ( ($GENVERSION, $LCSIM_CALL) = each(%LCSIM_HASH)) {

                    make_banner("RUNNING SPECSIM FOR GENVERSION: $GENVERSION ");

		    ### getting JOBNAME information
		    init_SPECSIM($USER_inFile);

		    ### getting distribution information
		    $EPOCH_DISTRIBUTION = $EPOCHDIST_HASH{$GENVERSION};
		    $PHASE_DISTRIBUTION = $PHASEDIST_HASH{$GENVERSION};
		    $SPECGEN_LIB = $GENLIB_HASH{$GENVERSION};

		    ### getting galaxy contamination information
		    $GALFRAC = $GALFRAC_HASH{$GENVERSION};
		    $GALFILE = $GALFILE_HASH{$GENVERSION}; 

		    ### getting specarg information

		    ( $IDEAL_EPOCHRANGE[0], $IDEAL_EPOCHRANGE[1], $IDEAL_CADENCE, $IDEAL_SNR_OBS[0], $IDEAL_SNR_OBS[1], $IDEAL_LAMINIT_REST, 
		      $IDEAL_LAMFINAL_REST, $IDEAL_DLAM ) = readvar_SPECARG($GENVERSION);

		    ### finding LCSIM directory, LCSIM README, and parsing for JOB keys
		    prep_SPECLOC($USER_inFile);

		    ### checking SIMSED parameters for GENMODEL
		    checkparam_SPECSIM();

		    ### get LCSIM ran seed to pass to SIMSED_extractSpec
		    $RANSEED = get_RANSEED($SIMREADME);

		    if ( $SPECGEN_LIB eq 'NONE' || $SPECGEN_LIB eq ''){
			runIdeal_SPECSIM();
		    } else {
			$SIM_ZMIN = $SIMMINZ_HASH{$GENVERSION};
			$SIM_ZMAX = $SIMMAXZ_HASH{$GENVERSION};
			runLib_SPECSIM();
		    }

                    print "\t\t FINISHED WITH SPECTRUM SIM FOR $GENVERSION \n\n";


	}            
} 

sub make_SUMMARY()
{

    my $lcsimdir; 
    my $listfile;
    my $datfile;
    my @LISTline;
    my $nlc;
    my $nspec;
    my $mean_spec;

    make_banner("SUMMARY OF SIMULATIONS");

    while ( ($GENVERSION, $LCSIM_CALL) = each(%LCSIM_HASH)) {

	$lcsimdir= $SNDATA_ROOT . "/SIM/" . $GENVERSION . "/";
	$listfile = $lcsimdir . $GENVERSION . ".LIST";
	check_fileexist("GENVERSION list file", $listfile);

	$nlc = 0;
	$nspec = 0;

	open(LISTFILE, $listfile) || die "ABORT: LCSIM LIST FILE not found. \n ";		
	while(<LISTFILE>){
	    chomp;
	    (@LISTline) = split(' ', $_ );
	
	    $datfile = $lcsimdir . $LISTline[0];
	    check_fileexist("LCSIM data file", $datfile);
	    $nlc = $nlc + 1;
	    
	    $KEY = "NSPEC:"; ### Must have this key to unpack
	    @TMP = sntools::parse_line($datfile, $NARG1, $KEY, $OPTQUIET) ;
	    if ( $TMP[0] eq '') {
		$nspec = $nspec + 0;
	    } else {
		$nspec = $nspec + $TMP[0];
	    }
	}

	$mean_spec = $nspec/$nlc; 
	print "\t GENVERSION   : $GENVERSION \n"; 
	print "\t Light Curves : $nlc \n";
	print "\t Spectra      : $nspec \n";
	print "\t Mean Spec/LC : $mean_spec \n\n";
    }

}


sub get_RANSEED($)
{
    my ($inFile) = shift;
    my $seed;

    $KEY = "SEED:"; ### Must have this key to unpack
    @TMP = sntools::parse_line($inFile, $NARG1, $KEY, $OPTABORT) ;
    $seed = @TMP[0];
    print "\t RANSEED \t $seed \n";
    
    return($seed);
}


sub prep_SPECLOC($)
{

    my ($inFile) = shift;
    my $lcsimdir;
    my @tmp2; ### use for parsing @TMP in NARG2 cases

    print "\nMaking required directory locations, checking for required files. \n";
    
       $lcsimdir= $SNDATA_ROOT . "/SIM/" . $GENVERSION . "/";
       $LCDUMP_FILE = $lcsimdir . $GENVERSION . ".DUMP";
       $SIMREADME = $lcsimdir . $GENVERSION . ".README"; 

       print "\t LCSIM GENVERSION \t $GENVERSION \n";

       ### MUST have access to dump file, readme file

       check_fileexist("LCSIM dump", $LCDUMP_FILE);
       print "\t LCSIM dump file located. \n";

       check_fileexist("LCSIM readme", $SIMREADME);
       print "\t LCSIM readme file located. \n";

    print "\nChecking for more input keys. \n";

    $KEY = "JOBNAME_SIMSEDextract:" ; ## OPTIONAL
    @TMP = sntools::parse_line($inFile, $NARG1, $KEY, $OPTQUIET) ;
    $JOBNAME_SIMSEDextract = @TMP[0];
    if ( $JOBNAME_SIMSEDextract eq '' ) {
	$JOBNAME_SIMSEDextract = $SNANA_DIR . "/bin" . "/SIMSED_extractSpec.exe" ;
	## $JOBNAME_SIMSEDextract = "/home/s1/rkessler/snana/bin/SIMSED_extractSpec.exe"; 
    }

    print "\nUSING required keys: \n\n";

    $KEY = "BRIEF_DESCRIPTION:"; ### Must have this key to unpack
    @TMP = sntools::parse_line($SIMREADME, $NARG7, $KEY, $OPTABORT) ;
    @tmp2 = split(' ', @TMP[0]);
    $GENMODEL = @tmp2[6];
    print "\t GENMODEL \t $GENMODEL \n";

    print "\nDone finding required files and keys. \n"; 
    
    make_banner("SPECSIM PREPARATION SUCCESSFUL");

} ### end prep_SPECLOC()

sub init_SPECSIM($)
{

    my ($inFile) = shift;

    print "\nChecking for JOBNAME_SIMSEDextract input key. \n";

    $KEY = "JOBNAME_SIMSEDextract:" ; ## OPTIONAL
    @TMP = sntools::parse_line($inFile, $NARG1, $KEY, $OPTQUIET) ;
    $JOBNAME_SIMSEDextract = @TMP[0];
    if ( $JOBNAME_SIMSEDextract eq '' ) {
	$JOBNAME_SIMSEDextract = $SNANA_DIR . "/bin" . "/SIMSED_extractSpec.exe" ;
	## $JOBNAME_SIMSEDextract = "/home/s1/rkessler/snana/bin/SIMSED_extractSpec.exe"; 
    }

    make_banner("SPECSIM INITIALIZATION SUCCESSFUL");

} ### end init_SPECSIM()

sub readLIB($)
{

    my $LIB_FILE = shift;
    my $nLIB;
    my $FILE_ID; 
    my $FILE_z; 
    ### ALL LIB HASHES ARE GLOBAL VARIABLES

    ### read LIB_FILE file into LIB hashes
    ### readCOLUMN_LIBCUT includes ISPEC <= 1, Z_HELIO between SIM_ZMIN, SIM_ZMAX. 
    readCOLUMN_LIBCUT($LIB_FILE, "SNID", "Z_HELIO", "ISPEC", 1, \%LIB_z_HASH);
    readCOLUMN_LIBCUT($LIB_FILE, "SNID", "t0", "ISPEC", 1, \%LIB_t0_HASH);
    readCOLUMN_LIBCUT($LIB_FILE, "SNID", "c", "ISPEC", 1, \%LIB_c_HASH);
    readCOLUMN_LIBCUT($LIB_FILE, "SNID", "x1", "ISPEC", 1, \%LIB_x1_HASH);
    readCOLUMN_LIBCUT($LIB_FILE, "SNID", "NSPEC", "ISPEC", 1, \%LIB_nspec_HASH);
    readCOLUMN_LIBCUT($LIB_FILE, "SNID", "ISPEC", "ISPEC", 1, \%LIB_ispec_HASH);
    
    ### debug check hash contents
    while ( ($FILE_ID, $FILE_z) = each(%LIB_z_HASH)) {
	print "\t\t  SNID:  $FILE_ID \n";
	print "\t\t  Z:     $FILE_z \n ";
	print "\t\t  t0:    $LIB_t0_HASH{$FILE_ID} \n";
	print "\t\t  c:     $LIB_c_HASH{$FILE_ID} \n";
	print "\t\t  x1:    $LIB_x1_HASH{$FILE_ID} \n";
	print "\t\t  NSPEC: $LIB_nspec_HASH{$FILE_ID} \n\n";
    }

    ### count keys in LIB hash
    $nLIB = keys %LIB_z_HASH;
    return $nLIB; 

}

sub readSIM($)
{
    my $FILE = shift;
    my $nSIM;

    ### read LCDUMP FILE into SIM hashes
    readCOLUMN($FILE, "CID", "Z", \%SIM_z_HASH);
    readCOLUMN($FILE, "CID", "S2x1", \%SIM_S2x1_HASH);
    readCOLUMN($FILE, "CID", "S2c", \%SIM_S2c_HASH);

    ### count keys in SIM hash
    $nSIM = keys %SIM_z_HASH;

} ### end sub readSIM($)

sub match_LIB_SIM(@)
{

    my ( $nLIB, $nSIM ) = @_;
    my $lastloop_iSIM;
    my $percent;
    my $FILE_ID;
    my $FILE_z;
    my $iSIM;
    my $SIM_ID;
    my $SIM_z;
    my $deltac;
    my $deltax1;
    my $deltaz;
    my $SIM_matchval;
    my @SIM_ID_sorted;
    my $tempkey;
    my $tempval;
	

    ### undefine, then initialize $LIB_usedval_HASH
    while ( ($FILE_ID, $FILE_z) = each(%LIB_z_HASH)) {
	$LIB_usedval_HASH{$FILE_ID} = 0.0 
    }

    $lastloop_iSIM = int ($nSIM / $nLIB ) * $nLIB - 1 ;
    $percent = ( $nSIM % $nLIB ) / ( 1.0 * $nLIB );
    $iSIM = 0;
    print "MATCH ONLY CONSIDERING Z VALUES\n"; 
    while ( ($SIM_ID, $SIM_z) = each(%SIM_z_HASH)) {
        ## calculate matchval for each SNID - CID pair
	while( ($FILE_ID, $FILE_z) = each(%LIB_z_HASH)) {
	    $deltac = abs( $SIM_S2c_HASH{$SIM_ID} - $LIB_c_HASH{$FILE_ID} );
	    $deltaz = abs( $SIM_z_HASH{$SIM_ID} - $LIB_z_HASH{$FILE_ID} );
	    $deltax1 = abs( $SIM_S2x1_HASH{$SIM_ID} - $LIB_x1_HASH{$FILE_ID} );
	    ## $SIM_matchval_HASH{$FILE_ID} = $deltac+$deltax1+$deltaz+$LIB_usedval_HASH{$FILE_ID};
            $SIM_matchval_HASH{$FILE_ID} = $deltaz+$LIB_usedval_HASH{$FILE_ID};
	}
        ## sort matchvals low to high to find best match
        ## increment used value for best match
        @SIM_ID_sorted =
            sort { $SIM_matchval_HASH{$a} cmp $SIM_matchval_HASH{$b} } keys %SIM_matchval_HASH;
	$SIM_bestmatch_HASH{$SIM_ID} = $SIM_ID_sorted[0];
	$LIB_usedval_HASH{$SIM_ID_sorted[0]} = $LIB_usedval_HASH{$SIM_ID_sorted[0]} + 10.0;
	print "\t\tTHE BEST MATCH SNID for $SIM_ID is $SIM_ID_sorted[0] \n" ;
        ## use rand to keep subset of SIM SNe for last partial SIM loop
	if ( $iSIM == $lastloop_iSIM ) {
	    while ( ($tempkey, $tempval) = each %LIB_usedval_HASH ) {
		if ( rand() >= $percent ) {
		    $LIB_usedval_HASH{$tempkey} = $LIB_usedval_HASH{$tempkey} + 10.0;
		}
	    }
	}
	$iSIM = $iSIM + 1;
    }
    print "\t Done matching LIB, SIM SNe \n";
} ### end sub match_LIB_SIM(@)

sub undefHASH()
{

    undef %LIB_c_HASH;
    undef %LIB_x1_HASH;
    undef %LIB_z_HASH;
    undef %LIB_nspec_HASH;
    undef %LIB_ispec_HASH;
    undef %LIB_t0_HASH;
    undef %LIB_usedval_HASH;
    undef %SIM_z_HASH;
    undef %SIM_S2x1_HASH;
    undef %SIM_S2c_HASH;
    undef %SIM_matchval_HASH;
    undef %SIM_bestmatch_HASH;
    undef %LIB_specmjd_HASH;
    undef %LIB_firstwave_HASH;
    undef %LIB_lastwave_HASH;
    undef %LIB_SNR6000_HASH;
    undef %LIB_SNR7000_HASH;
    undef %LIB_SNR8000_HASH;
    undef %LIB_SNR2000_HASH;
    undef %LIB_SNR3000_HASH;
    undef %LIB_SNR4000_HASH;
    undef %LIB_SNR5000_HASH;
    undef %LIB_SNR9000_HASH;
    undef %LIB_SNR10000_HASH;
    undef %LIB_SNRNBIN_HASH;
    undef %LIB_BINWAVE_HASH;

}

sub runLib_SPECSIM()
{
    my @DUMPline;
    my @LIBline;
    my $lcsim_dir;
    my $filename;
    my $SEDSIM_lcfile; ## check for training lcfile
    my $LCSIM_datfile; ## use lc sim dat file to get params
    my $trainlist_file; ##
    my $die_missing;
    my $LIB_count;
    my $SIM_count;
    my $nspec;
    my $check_spec; ### debuggin
    my $specmjd;
    my $ispec;
    my $iSIM;
    my $CID;
    my $SNID;
    my $SNID_t0;
    my $SNID_z;
    my $SNID_z1;
    my $CID_t0;
    my $CID_temp1;
    my $CID_temp2;
    my $CID_temp3;
    my $CID_z;
    my $CID_dlmu;
    my @CID_params;
    my $phase_file;
    my $log_file; 

    make_banner("BEGINNING < $GENVERSION > LIB SPECTRUM SIMULATION") ;
    print (" \t SIMULATION ZMIN: $SIM_ZMIN, SIMULATION ZMAX: $SIM_ZMAX\n" );

    ### cleaning HASHES
    undefHASH;

    ### read LIB FILE into LIB hashes
    $LIB_count = readLIB($SPECGEN_LIB);
    print "\t\t READ IN $LIB_count UNIQUE LIB SNe\n";

    ### read LCDUMP FILE into SIM hashes 
    $SIM_count = readSIM($LCDUMP_FILE);
    print "\t\t READ IN $SIM_count UNIQUE SIM SNe\n\n";    

    ### assign LIB SNe to SIM SNe 
    match_LIB_SIM($LIB_count, $SIM_count);

    ### make spectra
    $iSIM = 0;
    while ( ($CID, $SNID) = each(%SIM_bestmatch_HASH)) {
	undef %LIB_firstwave_HASH;
	undef %LIB_lastwave_HASH;
	undef %LIB_specmjd_HASH;
	undef %LIB_SNR2000_HASH;
	undef %LIB_SNR3000_HASH;
	undef %LIB_SNR4000_HASH;
	undef %LIB_SNR5000_HASH;
	undef %LIB_SNR6000_HASH;
	undef %LIB_SNR7000_HASH;
	undef %LIB_SNR8000_HASH;
	undef %LIB_SNR9000_HASH;
	undef %LIB_SNR10000_HASH;
	undef %LIB_SNRNBIN_HASH;
	undef %LIB_BINWAVE_HASH;

	### find number, mjds of spectra for SNID match SN. 
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SPECMJD", "SNID", $SNID, \%LIB_specmjd_HASH);	
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "FIRSTWAVE", "SNID", $SNID, \%LIB_firstwave_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "LASTWAVE", "SNID", $SNID, \%LIB_lastwave_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR2000", "SNID", $SNID, \%LIB_SNR2000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR3000", "SNID", $SNID, \%LIB_SNR3000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR4000", "SNID", $SNID, \%LIB_SNR4000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR5000", "SNID", $SNID, \%LIB_SNR5000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR6000", "SNID", $SNID, \%LIB_SNR6000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR7000", "SNID", $SNID, \%LIB_SNR7000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR8000", "SNID", $SNID, \%LIB_SNR8000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR9000", "SNID", $SNID, \%LIB_SNR9000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNR10000", "SNID", $SNID, \%LIB_SNR10000_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "SNRNBIN", "SNID", $SNID, \%LIB_SNRNBIN_HASH);
	readCOLUMN_EQCUT($SPECGEN_LIB, "ISPEC", "BINWAVE", "SNID", $SNID, \%LIB_BINWAVE_HASH);
	
	### find number of spectra expected for SNID match SN
	$nspec = $LIB_nspec_HASH{$SNID};
	$check_spec = keys %LIB_specmjd_HASH;
        print "debug: $CID $SNID genspec NPHASE is $nspec, check_spec NPHASE is $check_spec \n";
	unless ( $nspec <= $check_spec ) { die "ABORT: number of spectra, SNID LIB mjds don't match \n"; }

	make_banner("Generating Spectra for SN $CID");
	print "\t Creating $nspec $SNID match spectra for SIM SN $CID\n";

	if ( $nspec > 0 ) {
	    ### make file addresses to be used
	    $phase_file = $GENVERSION . "_" . $CID . "_phaselist.temp"; 
	    $log_file = $GENVERSION . "_" . $CID . "_phaselist.log"; 
	    $lcsim_dir= $SNDATA_ROOT . "/SIM/" . $GENVERSION . "/";
	    $LCSIM_datfile = $lcsim_dir . $GENVERSION . "_SN0" . $CID . ".DAT";
	    $SEDSIM_lcfile = $LCSIM_datfile;
	    check_fileexist("SEDSIM lc", $SEDSIM_lcfile);
	    
	    ### get parameters from dat file and hashes
	    @SN_SIM_PARAMS = ();
	    @SN_SIM_PARAMS = getparam_SPECSIM($LCSIM_datfile);
	    $SNID_t0 = $LIB_t0_HASH{$SNID};
	    $SNID_z = $LIB_z_HASH{$SNID};
	    $SNID_z1 = 1. + $SNID_z;
	    
	    ### make phase file ( required by SIMSED_extractSpec.exe ) 
	    if ( -e $phase_file ) { unlink($phase_file);}
	    open(PHASEOUT, ">>$phase_file") || die "Could not open temporary phase file.\n";
	    $printtype = 'SPEC';
	    while (($ispec, $specmjd) = each(%LIB_specmjd_HASH)){
                ## print ("debug: $SNID, $ispec , $specmjd, $SNID_t0, $SNID_z, $CID_z, $CID_t0 \n");
		## print ("debug: $IDEAL_EPOCHRANGE[0], $IDEAL_EPOCHRANGE[1]\n"); 
		$printphase = ($specmjd-$SNID_t0)/$SNID_z1;
		$printdate = ($printphase*(1+$CID_z)) + $CID_t0;
		$printfile = $SEDSIM_lcfile;
		$printlaminit_rest = $LIB_firstwave_HASH{$ispec}/$SNID_z1;
		$printlamfinal_rest = $LIB_lastwave_HASH{$ispec}/$SNID_z1;
		$printSNR1 = $LIB_SNR2000_HASH{$ispec}; 
		$printSNR2 = $LIB_SNR3000_HASH{$ispec}; 
		$printSNR3 = $LIB_SNR4000_HASH{$ispec}; 
		$printSNR4 = $LIB_SNR5000_HASH{$ispec}; 
		$printSNR5 = $LIB_SNR6000_HASH{$ispec}; 
		$printSNR6 = $LIB_SNR7000_HASH{$ispec}; 
		$printSNR7 = $LIB_SNR8000_HASH{$ispec}; 
		$printSNR8 = $LIB_SNR9000_HASH{$ispec}; 
		$printSNR9 = $LIB_SNR10000_HASH{$ispec}; 
		$printSNRLAMBIN = $LIB_BINWAVE_HASH{$ispec};

		if ( $printphase >= $IDEAL_EPOCHRANGE[0] && $printphase <= $IDEAL_EPOCHRANGE[1] ) {
		    if ( $OUTPUT_FORMAT eq "ASCII") {
			open(SPECHEAD, ">>$printfile");
			write(SPECHEAD);
			close(SPECHEAD);
		    }
		    print PHASEOUT sprintf("%13.3f %20s %13.3f %13.3f %13.3f ", 
					   $printphase, $printfile, $printlaminit_rest, $printlamfinal_rest, $printSNRLAMBIN);
		    print PHASEOUT sprintf("%13.3f %13.3f %13.3f %13.3f %13.3f ", 
					   $printSNR1, $printSNR2, $printSNR3, $printSNR4, $printSNR5);
		    print PHASEOUT sprintf("%13.3f %13.3f %13.3f %13.3f \n", 
					   $printSNR6, $printSNR7, $printSNR8, $printSNR9);

		} else {
		}
	    }
	    close(PHASEOUT);
	    unless ( -e $phase_file ) { die "missing phase file \n"; }
	    
	    call_SIMSEDextractSpec($phase_file, $log_file);

	    increment_RANSEED();
	    
	    print "\Cleaning up local directory\n";
	    #unlink(<*.SPEC>);
	    unlink $phase_file;
	    unlink $log_file;
	} else {
	    print "\t No spectra to generate. Skipping to next SIM SN \n";
	}



    }
    
    

    make_banner("DONE WITH < $GENVERSION > LIB SPECTRUM SIMULATION") ;

} 

sub readCOLUMN_LIBCUT(@)
{

    ### THIS SPECIAL VERSION OF readCOLUMN_CUT()
    ### ALSO CUTS BY SIM_ZMIN, SIM_ZMAX. 
    ### REQUIRES LIB COLUMN HEADER Z_HELIO

    ### ASSUMING FILE SELF-DOCUMENTING
    ### FIRST LINE is format: NVAR=<##>
    ### SECOND LINE is format: VARNAMES: <name1> <name2> .. <namen>
    ### REMAINING LINES are: SN: <col1> <col2> .. <coln>

    my ( $FILE, $KEY, $VAR, $CUT, $CUTVAL, $HASHREF ) = @_;
    my @readline;
    my $nline;
    my $iline;
    my $line_flag;
    my $var_column;
    my $key_column;
    my $cut_column;
    my $z_column;

    ### print "debug In CUT - CUT is $CUT, CUTVAL is $CUTVAL \n" ;

    open(FILE, $FILE) || die "ABORT: FILE $FILE not found. \n ";

    while(<FILE>){
	chomp;
        (@readline) = split(' ', $_ );
        if ( $readline[0] =~ /^NVAR:/ ) {
            $nline = $readline[1];
	}   elsif ( $readline[0] =~ /^VARNAMES:/ ) {
            $line_flag = 0;
            $iline = 0;
            while ( $line_flag < 4 and $iline <= $nline ) {
                ### print "debug : In VAR loop: $iline, $readline[$iline] \n";
                if ($readline[$iline] eq $VAR) {
                    ### print "$VAR, $readline[$iline], $var_column \n";
                    $line_flag = $line_flag+1;
                    $var_column = $iline;
                }
                if ($readline[$iline] eq $KEY) {
                    $line_flag = $line_flag+1;
                    $key_column = $iline
                    }
                if ($readline[$iline] eq $CUT){
                    $line_flag = $line_flag+1;
                    $cut_column = $iline;
                }
		if ($readline[$iline] eq "Z_HELIO"){
		    $line_flag = $line_flag+1;
		    $z_column = $iline;
		}
                $iline = $iline + 1;
            }
            unless ( $line_flag == 4 ) {
                die "ABORT: COLUMN HEADER $VAR , $KEY, or Z_HELIO NOT FOUND in FILE $FILE.\n";
            }
        }   elsif ( $readline[0] =~ /^SN:/ ) {
            if ( $readline[$cut_column] <= $CUTVAL and 
		 $readline[$z_column] >= $SIM_ZMIN and 
		 $readline[$z_column] <= $SIM_ZMAX ) {
                $ { $HASHREF }{ $readline[$key_column] } =  $readline[$var_column];
            }
        }
    }
    close(FILE);
} ## end sub readCOLUMN_LIBCUT

sub readCOLUMN_EQCUT(@)
{
    ### ASSUMING FILE SELF-DOCUMENTING
    ### FIRST LINE is format: NVAR=<##>
    ### SECOND LINE is format: VARNAMES: <name1> <name2> .. <namen>
    ### REMAINING LINES are: SN: <col1> <col2> .. <coln>

    my ( $FILE, $KEY, $VAR, $CUT, $CUTVAL, $HASHREF ) = @_;
    my @readline;
    my $nline;
    my $iline;
    my $line_flag;
    my $var_column;
    my $key_column;
    my $cut_column;

    open(FILE, $FILE) || die "ABORT: FILE $FILE not found. \n ";

    while(<FILE>){
	chomp;
        (@readline) = split(' ', $_ );
        if ( $readline[0] =~ /^NVAR:/ ) {
            $nline = $readline[1];
	}   elsif ( $readline[0] =~ /^VARNAMES:/ ) {
            $line_flag = 0;
            $iline = 0;
            while ( $line_flag < 3 and $iline <= $nline ) {
                ### print "debug : In VAR loop: $iline, $readline[$iline] \n";
                if ($readline[$iline] eq $VAR) {
                    ### print "$VAR, $readline[$iline], $var_column \n";
                    $line_flag = $line_flag+1;
                    $var_column = $iline;
                }
                if ($readline[$iline] eq $KEY) {
                    $line_flag = $line_flag+1;
                    $key_column = $iline
                    }
                if ($readline[$iline] eq $CUT){
                    $line_flag = $line_flag+1;
                    $cut_column = $iline;
                }
                $iline = $iline + 1;
            }
            unless ( $line_flag == 3 ) {
                die "ABORT: COLUMN HEADER $VAR or $KEY NOT FOUND in FILE $FILE.\n";
            }
        }   elsif ( $readline[0] =~ /^SN:/ ) {
            if ( $readline[$cut_column] eq $CUTVAL ){
                $ { $HASHREF }{ $readline[$key_column] } =  $readline[$var_column];
            }
        }
    }
    close(FILE);
} ## end sub readCOLUMN_EQCUT

sub readCOLUMN_CUT(@)
{
    ### ASSUMING FILE SELF-DOCUMENTING
    ### FIRST LINE is format: NVAR=<##>
    ### SECOND LINE is format: VARNAMES: <name1> <name2> .. <namen>
    ### REMAINING LINES are: SN: <col1> <col2> .. <coln>

    my ( $FILE, $KEY, $VAR, $CUT, $CUTVAL, $HASHREF ) = @_;
    my @readline;
    my $nline;
    my $iline;
    my $line_flag;
    my $var_column;
    my $key_column;
    my $cut_column;

    ### print "debug In CUT - CUT is $CUT, CUTVAL is $CUTVAL \n" ;

    open(FILE, $FILE) || die "ABORT: FILE $FILE not found. \n ";

    while(<FILE>){
        chomp;
        (@readline) = split(' ', $_ );
        if ( $readline[0] =~ /^NVAR:/ ) {
            $nline = $readline[1];
        }   elsif ( $readline[0] =~ /^VARNAMES:/ ) {
            $line_flag = 0;
            $iline = 0;
            while ( $line_flag < 3 and $iline <= $nline ) {
		### print "debug : In VAR loop: $iline, $readline[$iline] \n"; 
                if ($readline[$iline] eq $VAR) {
                    ### print "$VAR, $readline[$iline], $var_column \n";
                    $line_flag = $line_flag+1;
                    $var_column = $iline;
                }
                if ($readline[$iline] eq $KEY) {
                    $line_flag = $line_flag+1;
                    $key_column = $iline
		    }
		if ($readline[$iline] eq $CUT){
		    $line_flag = $line_flag+1;
		    $cut_column = $iline;
		}
                $iline = $iline + 1;
            }
            unless ( $line_flag == 3 ) {
                die "ABORT: COLUMN HEADER $VAR or $KEY NOT FOUND in FILE $FILE.\n";
            }
        }   elsif ( $readline[0] =~ /^SN:/ ) {
	    if ( $readline[$cut_column] <= $CUTVAL ){ 
		$ { $HASHREF }{ $readline[$key_column] } =  $readline[$var_column];
	    }
        }
    }
    close(FILE);
} ## end sub readCOLUMN_CUT


sub readCOLUMN(@)
{
    ### ASSUMING FILE SELF-DOCUMENTING
    ### FIRST LINE is format: NVAR=<##>
    ### SECOND LINE is format: VARNAMES: <name1> <name2> .. <namen>
    ### REMAINING LINES are: SN: <col1> <col2> .. <coln>

    my ( $FILE, $KEY, $VAR, $HASHREF ) = @_;
    my @readline;
    my $nline;
    my $iline;
    my $line_flag;
    my $var_column; 
    my $key_column;

    open(FILE, $FILE) || die "ABORT: FILE $FILE not found. \n ";

    while(<FILE>){
	chomp;
	(@readline) = split(' ', $_ );
	if ( $readline[0] =~ /^NVAR:/ ) {
	    $nline = $readline[1];
	}   elsif ( $readline[0] =~ /^VARNAMES:/ ) {
	    $line_flag = 0;
	    $iline = 0;
	    while ( $line_flag < 2 and $iline <= $nline ) {
		if ($readline[$iline] eq $VAR) {
		    ### print "$VAR, $readline[$iline], $var_column \n"; 
		    $line_flag = $line_flag+1;
		    $var_column = $iline;
		}
		if ($readline[$iline] eq $KEY) {
		    $line_flag = $line_flag+1;
		    $key_column = $iline 
		}
		$iline = $iline + 1;
	    }
	    unless ( $line_flag == 2 ) { 
		die "ABORT: COLUMN HEADER $VAR or $KEY NOT FOUND in FILE $FILE.\n"; 
	    }  
	}   elsif ( $readline[0] =~ /^SN:/ ) {
	    $ { $HASHREF }{ $readline[$key_column] } =  $readline[$var_column];
	}
    }
    close(FILE);
} ## end sub_readCOLUMN


sub runIdeal_SPECSIM()
{
    my @DUMPline;
    my $CID;
    my $NSNE;
    my $iSNE;
    my $lcsim_dir;
    my $filename; 
    my $SEDSIM_lcfile; ## check for training lcfile
    my $LCSIM_datfile; ## use lc sim dat file to get params
    my $trainlist_file; ## 
    my $die_missing;
    my @output;

make_banner("BEGINNING < $GENVERSION > IDEAL SPECTRUM SIMULATION") ;

print "\nReading SN IDs from LCSIM DUMP file. \n";
    print "\t Dump file: $LCDUMP_FILE \n"; 

    ### READ DUMP FILE ONCE TO GET NSNE ###
    $NSNE = 0;
    open(DUMPFILE, $LCDUMP_FILE) || die "ABORT: LCSIM DUMP FILE not found. \n ";
    while(<DUMPFILE>){
        chomp;
        (@DUMPline) = split(' ', $_ );
        if ( $DUMPline[0] =~ /^SN:/ ) {
	    $NSNE ++ ;
	}
    }
    close(DUMPFILE);

    ###   GET LIST OF NEPOCHS FOR NSNE    ###	
    ### (@NEPOCH_LIST made in subroutine) ###
    get_nepochlist($NSNE, $PHASE_DISTRIBUTION);

    ### GO THROUGH DUMP FILE AGAIN TO GENERATE SPECTRA ###
    $iSNE = 1;
    open(DUMPFILE, $LCDUMP_FILE) || die "ABORT: LCSIM DUMP FILE not found. \n ";
    while(<DUMPFILE>){
        chomp;
        (@DUMPline) = split(' ', $_ );
        if ( $DUMPline[0] =~ /^SN:/ ) {

                ## for each SN, make file names based on CID
           	$CID = $DUMPline[1];
		$printCID = $CID; 

		make_banner("Generating Spectra for SN $CID");
		
                $lcsim_dir= $SNDATA_ROOT . "/SIM/" . $GENVERSION . "/";
                $LCSIM_datfile = $lcsim_dir . $GENVERSION . "_SN0" . $CID . ".DAT";
		$SEDSIM_lcfile = $LCSIM_datfile;

		check_fileexist("SEDSIM lc", $SEDSIM_lcfile);

               ## get parameters from dat file
		@SN_SIM_PARAMS = ();
                @SN_SIM_PARAMS = getparam_SPECSIM($LCSIM_datfile);
 
	       ## generate spectra for this LC

		genspec_GAMMA_SPECSIM($NEPOCH_LIST[$iSNE], $SEDSIM_lcfile);
		$iSNE = $iSNE +1 ; 

		    #### GAMMA version genspec_SPECSIM(@SN_SIM_PARAMS) incorporates FLAT option
		    ###  may want to change name back to genspec_SPECSIM eventually


    
	    }
    }
        close(DUMPFILE);

    make_banner("DONE WITH < $GENVERSION > IDEAL SPECTRUM SIMULATION") ;


} 

sub checkparam_SPECSIM()
{

### read SIMSED SED.INFO file to get binary vector of USEPARAMS

    my $sedinfofile;
    my ($npar, $ipar);


    $sedinfofile = $SNDATA_ROOT . "/models/SIMSED/" . $GENMODEL . "/SED.INFO"; 
    check_fileexist("SED.INFO file", $sedinfofile);

    print "\n Starting SED.INFO file $sedinfofile parameter check.\n";

    $KEY = "NPAR:" ;
    @TMP = sntools::parse_line($sedinfofile, $NARG1, $KEY, $OPTABORT) ;
    $npar = @TMP[0];

    $KEY = "PARNAMES:" ;
    @TMP = sntools::parse_line($sedinfofile, $NARGALL, $KEY, $OPTABORT) ;
    @PARNAMES = split(' ', @TMP[0]);

    unless ( @PARNAMES == $npar ) { die "Number of parameters and parameter names conflict\n"; }

### read SIM .README file to distinguish between real and dummy parameters
### $SIMREADME variable was set in prep_SPECLOC()

    $ipar = 0;
    while ( $ipar < $npar ) {
        @TMP = `grep Peak  $SIMREADME | grep $PARNAMES[$ipar] `;
	printf("DEBUG: double grep result is @TMP \n");
	if ( @TMP[0] ne '' ) { $USEPARVALS[$ipar] = 1; }
	else { $USEPARVALS[$ipar] = 0; }
	printf("$PARNAMES[$ipar]: useval is $USEPARVALS[$ipar]\n");
        $ipar++;
    }

}


sub split_param(@)
{
## pass in a line from DAT file beginning with SIMSED_ 
## determine if line is parameter value, and if so, 
## return which parameter and parameter value

    my ($linezero, $lineone) = @_;
    my @name_split;
    my ($par_val, $par_location); 
    my $ival; 

    if ( $linezero eq "SIMSED_NPAR:" ) {

	$par_val = -999.;
	$par_location = -9;
    } else { 
	$par_val = $lineone; 
	@name_split = split( /\Q(\E/ , $linezero );
	$name_split[1] =~ s#\Q)\E:##;

       for ( $ival = 0; $ival < @PARNAMES; $ival++ ){
	   if ( $name_split[1] eq $PARNAMES[$ival] ) {
	       $par_location = $ival;
	   }
       }

       unless ( $USEPARVALS[$par_location] == 1 ) { $par_val = -999. };

    }
			 

    return ($par_location, $par_val);

}

sub getparam_SPECSIM($) 
{

    my ( $paramfile ) = shift;
    my ( @line);
    my ( @debug_params);
    my ( @params, $ival);
    my ( $sim_z);
    my ( $sim_dlmu);
    my ( $peakmjd);
    my ( $sim_trestmin_round );
    my ( $sim_trestmax_round );
    my ( $first_obs_date);
    my ( $first_obs_counter) = 0;

print "\nStarting SN simulation parameter read.\n"; 

    ### initializing param array
    for ( $ival = 0; $ival < @PARNAMES; $ival++ ){
	$params[$ival] = 0.0 ;
    }

    open (PARAMFILE, $paramfile) || die "Couldn't open PARAMFILE $paramfile \n";
    while (<PARAMFILE>){
        chomp;
        @line = split(' ', $_);
        if ( $line[0] =~ /^SIM_DLMU/) { $sim_dlmu = $line[1]; }

	### multiple SIMSED param keywords exist - include all possibilities here
        if ( $line[0] =~ /^SIMSED_/) { 
	    ### split_param returns param_location and param_value
            ### param_location is negative if line was SIMSED_NPAR: 
	    @debug_params = split_param($line[0], $line[1]);	    
	    printf("DEBUG split_param: $line[0], $line[1], @debug_params, @PARNAMES\n");
	    if ( $debug_params[0] >= 0 ) { 
		$params[$debug_params[0]] = $debug_params[1] ; }
	}
    
        
        if ( $line[0] =~ /^SIM_REDSHIFT/) { $sim_z = $line[1]; }
        if ( $line[0] =~ /^SIM_TRESTMIN/) { $sim_trestmin_round = floor($line[1]); }
        if ( $line[0] =~ /^SIM_TRESTMAX/) { $sim_trestmax_round = ceil($line[1]); }
        if ( $line[0] =~ /^SIM_PEAKMJD/ ) { $peakmjd = sprintf("%.3f", $line[1]); }
        if ( ($line[0] =~ /^OBS/) && ($first_obs_counter == 0 )) { $first_obs_date = $line[1]; $first_obs_counter = 1; }

    }
    close(PARAMFILE);

print "\nDone reading SN simulation parameters.\n";
print "\n params are @params \n"; 

return ($sim_z, $sim_trestmin_round, $sim_trestmax_round, $peakmjd, $sim_dlmu, $first_obs_date, @params );
}


sub parse_inFile_SALT2train($)
{
    use Switch;

    my $inFile_SALT2train = shift;
    my ( @line ) ; ## use to hold lines in file
    my $trash_value; 
    my $narg;
    my $iarg;
    my @tmp2;
    my $STYPE;
    my $DTYPE;
    my $NVAR;
    my $check_flag;


    ### ENSURE inFILE existence HAS ALREADY BEEN CHECKED

    ### Certain Format Assumptions are being made
    ### 1) LCSIM_INPUT should be first non-job keyword input
    ### 2) Each LCSIM_INPUT line may be followed by the 
    ###    optional spectrum simulation keyword lines, in
    ###    any order. Multiple SPECSIM_ARGUMENTS lines may
    ###    be included. 
    ###
    ###         -- SPECGEN_LIB 
    ###                   OR 
    ###         -- NPHASE_DISTRIBUTION:
    ###         -- EPOCH_DISTRIBUTION:
    ###
    ###                   AND
    ###         -- SPECSIM_ARGUMENTS:

    open (INFILE_SALT2TRAIN, $inFile_SALT2train) or  die "ERROR: Cannot open input file $inFile_SALT2train \n";
    make_banner("Parsing input file: $inFile_SALT2train");
    

    while ( <INFILE_SALT2TRAIN> ){
        chomp;
        @line = split(' ', $_ );
        switch ($line[0]){
	    case ("END:")           {
		if ( $GENVERSION ne '' ) { write_GENVERSIONhash() };
		init_INPUTKEYS;
	    }
            case("LCSIM_INPUT:")    {
		parse_LCSIMLINE(\@line);
                if ( $GENVERSION eq '' ) {die "Must have GENVERSION as part of LCSIM_INPUT \n";}
            }
	    case("SPECSIM_ARGUMENTS:"){
		shift(@line);
		push( @IDEAL_SPECARGS, @line ); 
	    }
            else {
		if ( $line[1] eq "LCSIM_INPUT:" ) {
		    $trash_value = shift(@line); 
		    parse_LCSIMLINE(\@line);
		    if ( $GENVERSION eq '' ) {die "Must have GENVERSION as part of LCSIM_INPUT \n";}
		}
	    }

	}
    }

    close(INFILE_SALT2TRAIN);

    make_banner("SIMULATIONS TO BE RUN ARE:");

    if ($RUN_LCSIM == "T" ) { 
	print "The following LCSIM calls will be run: \n";
	while ( ($GENVERSION, $LCSIM_CALL) = each(%LCSIM_HASH)) {
	    unless ( $LCSIM_CALL eq "NULL" ) {
		print "\t\t snlc_sim.exe $LCSIM_CALL \n";
	    }
	}
    }
    if ($RUN_SPECSIM == "T" ) {
	print "\nThe following SPECSIM jobs will be run: \n";
		print "\t\t *********************************\n";
	while ( ($GENVERSION, $LCSIM_CALL) = each(%LCSIM_HASH)) {
		print "\t\t  GENVERSION: $GENVERSION \n";
		print "\t\t  GENLIB:     $GENLIB_HASH{$GENVERSION} \n ";
		print "\t\t  EPOCHDIST:  $EPOCHDIST_HASH{$GENVERSION} ";
		    @VAR_ARRAY = read_EPOCHvar_SPECSIM($GENVERSION);
		    print "\t\t  EPOCH DIST parameters: @VAR_ARRAY \n";
		print "\t\t  PHASEDIST:  $PHASEDIST_HASH{$GENVERSION} ";
		    @VAR_ARRAY = read_PHASEvar_SPECSIM($GENVERSION);
		    print "\t\t  PHASE DIST parameters: @VAR_ARRAY \n";
		@VAR_ARRAY = readvar_SPECARG($GENVERSION);
		print "\t\t  SPECARG parameters: $VAR_ARRAY[0]";
		print " $VAR_ARRAY[1] $VAR_ARRAY[2] $VAR_ARRAY[3]";
		print " $VAR_ARRAY[4] $VAR_ARRAY[5] $VAR_ARRAY[6] \n";
		print "\t\t *********************************\n";
	    }
    }

    

}




sub make_banner($)
{

    my $banner_string = shift;
    my $SYMLINE = "\n _____________________________ \n\n" ;

    print "$SYMLINE";
    print "$banner_string";
    print "$SYMLINE";

} # End of make_banner

sub make_subbanner($)
{
    my $banner_string = shift;
    my $SYMLINE = "\n\t\t _____________________________ \n\n" ;

    print "$SYMLINE";
    print "\t\t $banner_string";
    print "$SYMLINE";

}

sub jobinfo_dump(@)
{ 
    my ($job, $input, $log) = @_;
    
    # if no job, log, or input is required, call with "none" in the relevant place

    print "JOB CALLED AS: \n";
    print "\t <JOB NAME> <JOB INPUT>  >& <JOB LOG> \n\n";
    print "\t USING \n";
    unless ($job eq 'none') {print "\t JOB NAME:  $job \n";}
    unless ($job eq 'none') {print "\t JOB INPUT: $input \n";}
    unless ($job eq 'none') {print "\t JOB LOG:   $log \n";}
}


sub check_fileexist(@)
{
    my ($key, $file) = @_;
    my $die_string; 

    $die_string = "\n\nABORT: Cannot find $key file $file \n\n";
    die $die_string unless -e $file; 

} # End of check_fileexist

sub check_goodexit(@)
{
    ## ideally here, "key" is the name of the job whose log is being checked
    my $JOBname;
    my $logfile;
    my $ABORTkey;
    my $ABORT; 
    my $die_string;

    ($JOBname, $logfile, $ABORTkey) = @_;
    if ( $ABORTkey eq '' ) { $ABORTkey = "ABORT"; }
    
    @TMP = sntools::parse_line($logfile, 99, $ABORTkey, $OPTQUIET) ;
    $ABORT = @TMP; ### use number of entries in @TMP to test ABORT
    $die_string = "$JOBname job failed. ABORT simprep. \n\n";
    if ( $ABORT > 0 ) { 
	print "$ABORTkey string found in log $logfile .\n"; 
	die $die_string ;}

    
} # End of check_goodexit


sub parse_LCSIMLINE($) 
{

    my $narg;
    my $iarg;
    my @tmp2;
 
    my ( @lcsimline ) = @{$_[0]};

    $narg = @lcsimline; ## pull out number of elements in line vector
    $LCSIM_INPUT = @lcsimline[1];

    if ( $narg == 2 ) {
	$LCSIM_CALL = @lcsimline[1];

        ### parse LCSIM_INPUT file for GENVERSION ###
	check_fileexist("LCSIM_INPUT", $LCSIM_INPUT);
	$KEY = "GENVERSION:" ;
	@TMP = sntools::parse_line($LCSIM_INPUT, $NARG1, $KEY, $OPTABORT) ;
	$GENVERSION = @TMP[0];
        ### done parsing ###
	if ( $GENVERSION eq '' ) {die "GENVERSION key not found in file $LCSIM_INPUT \n";}

    } elsif ($narg > 2 ) {
	$LCSIM_CALL = join(' ', @lcsimline[1...$narg-1]);
	$iarg = 0;

        ### get GENVERSION from LCSIM_INPUT vector ###
	while ($iarg < $narg ) {
	    if ( @lcsimline[$iarg] eq 'GENVERSION' ) { $GENVERSION = @lcsimline[$iarg+1] }
	    $iarg = $iarg + 1;
	}

    }

    ### parse LCSIM_CALL for zmin, zmax ###
    ### if they're not there, parse LCSIM_INPUT file for zmin, zmax ###
    $iarg = 0;
    $SIM_ZMIN = ''; 
    $SIM_ZMAX = ''; 
    while ( $iarg < $narg ){
	if ( @lcsimline[$iarg] eq 'GENRANGE_REDSHIFT' ) {
	    $SIM_ZMIN = @lcsimline[$iarg+1];
	    $SIM_ZMAX = @lcsimline[$iarg+2];
	    $iarg = $iarg + 2;
	} else {
	    $iarg = $iarg + 1;
	}
    }
    if ( $SIM_ZMIN eq '' ) {
        $KEY = "GENRANGE_REDSHIFT:"; ### Must have this key to unpack
	@TMP = sntools::parse_line($LCSIM_INPUT, $NARG7, $KEY, $OPTABORT) ;
	@tmp2 = split(' ', @TMP[0]);
	$SIM_ZMIN = @tmp2[0];
	$SIM_ZMAX = @tmp2[1];
	
    }

}

sub poisson_rand($) {

    my $lambda = shift;
    my ($L);
    my $k = 0;
    my $p = 1;

    $L = exp(-1*$lambda);
    while ( $p > $L ){
	$k = $k + 1;
	$p = $p * rand();
    }
    
    return ($k-1);
}

sub exp_rand($) {

    my $lambda = shift;
    my $T;

    $T = -1 * log(1-rand()) / $lambda ;
    
    print "in EXP_RAND: T is $T \n";
    return $T;
}

sub gamma_rand(@) {
    
    # this is really erlang function 
    # see http://heather.cs.ucdavis.edu/~matloff/156/PLN/RandNumGen.pdf 

    my ( $k, $theta ) = @_;
    my $lambda = 1.0/($theta*1.0) ;
    my $sum = 0;
    my $i = 0;

    for ( $i <= $k ) {
	$sum = $sum + exp_rand($lambda);
	$i = $i + 1;
    }
	    
    return ($sum/$lambda);

}

sub gaussian_rand(@) {
    
    my ($mean, $sigma) = @_;
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $sigma * $u1 * $w + $mean;
    $g1 = $sigma * $u2 * $w + $mean;
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}

sub genspec_GAMMA_SPECSIM(@)
{
    my ( $nphase, $SEDSIM_lcfile) = @_; 
    my ($z, $trestmin, $trestmax, $peakmjd, $sim_dlmu, $firstobsdate, @params ) = @SN_SIM_PARAMS;
    my $phase_file = "phase_list.temp" ; ## trash
    my @phase_line;
    my $log_file = 'log.temp' ;
    my $temp_header_file;
    my $temp_readspec_file;
    my $temp_writespec_file;
    my $temp_filename1;
    my $temp_filename2;
    my $temp_filename;
    my @phaselist;
    my $nphaselist;
    my $ilist;
    my $iphase;
    my $phaseflag;
    my $absval;

    print "\tSetting up for spectrum generation.\n";

    ($nphase, @phaselist ) = phase_generator($nphase);
    print "\tGrabbed first $nphase, @phaselist \n";

    if ( -e $phase_file ) { unlink($phase_file);}
    open(PHASEOUT, ">>$phase_file") || die "Could not open temporary phase file.\n";

    $printtype = 'SPEC';

    $iphase=0;
    while ( $iphase < $nphase ){
        $printphase = $phaselist[$iphase];
        $printdate = ($iphase*(1+$z))+$peakmjd;
	$temp_writespec_file = $SEDSIM_lcfile; 
        $printfile = $temp_writespec_file;
	$printlaminit_rest = $IDEAL_LAMINIT_REST;
	$printlamfinal_rest = $IDEAL_LAMFINAL_REST;
	$printSNR1 = $printSNR2 = $printSNR3 = $printSNR4 = $IDEAL_SNR_OBS[0];
	$printSNR5 = $printSNR6 = $printSNR7 = $printSNR8 = $printSNR9 = $IDEAL_SNR_OBS[0];
	$printSNRLAMBIN = $IDEAL_SNR_OBS[1];
	if ( $OUTPUT_FORMAT eq "ASCII" ) {
	    open(SPECHEAD, ">>$temp_writespec_file");
	    write(SPECHEAD);
	    close(SPECHEAD);
	}
                    print PHASEOUT sprintf("%13.3f %20s %13.3f %13.3f %13.3f ",
                                           $printphase, $printfile, $printlaminit_rest, $printlamfinal_rest, $printSNRLAMBIN);
                    print PHASEOUT sprintf("%13.3f %13.3f %13.3f %13.3f %13.3f ",
                                           $printSNR1, $printSNR2, $printSNR3, $printSNR4, $printSNR5);
                    print PHASEOUT sprintf("%13.3f %13.3f %13.3f %13.3f \n",
                                           $printSNR6, $printSNR7, $printSNR8, $printSNR9);

        $iphase = $iphase + 1;
    }
    close(PHASEOUT);
    unless ( -e $phase_file ) { die "missing phase file \n"; }

    call_SIMSEDextractSpec($phase_file, $log_file);

    increment_RANSEED();

    print "\Cleaning up local directory\n";
    unlink(<*.SPEC>);
    unlink(<*.temp>);

}


sub call_SIMSEDextractSpec(@)
{
    ## uses globals $GENMODEL, $IDEAL_LAMINIT_REST, $IDEAL_LAMFINAL_REST, @IDEAL_SNR_OBS, $RANSEED
    ##              $GALFRAC, $GALFILE, 
    ##              $JOBNAME_SIMSEDextract, $OUTPUT_FORMAT, @SN_SIM_PARAMS

    my ($phaselist_file, $log_file) = @_; 
    my ($z, $trestmin, $trestmax, $peakmjd, $sim_dlmu, $firstobsdate, @params ) = @SN_SIM_PARAMS;
    my $arglist1;
    my $arglist2;
    my $arglist3;
    my $arglist4;
    my $arglist5;
    my @output;
    my $ABORT_die;
    my ($nparams, $iparams, $ngoodpar);

    $nparams = @params;
    $ngoodpar = 0;
    
    for ($iparams = 1; $iparams <= $nparams; $iparams++){
	unless ( $params[$iparams] < -100 || $params[$iparams] > 100) {
	    $ngoodpar = $ngoodpar + 1;
	}
    }    
    
    if ($ngoodpar < 2) {
	    $ABORT_die = "\nABORT: Too many crazy SIMSED parameters. \n";
	    die $ABORT_die;
	}

    $arglist1 = " --v  $GENMODEL --tlist $phaselist_file ";
    $arglist2 = " --plist @params --lam $IDEAL_LAMINIT_REST $IDEAL_LAMFINAL_REST ";
    $arglist3 = " --z $z --dm $sim_dlmu --snr $IDEAL_SNR_OBS[0] $IDEAL_SNR_OBS[1] ";
    $arglist4 = " $OUTPUT_FORMAT --T0 $peakmjd --seed $RANSEED  OPT_ZEROFLUX ";    
    $arglist5 = " --gfrac $GALFRAC --gfile $GALFILE ";

    print "Calling SIMSED_extractSpec.exe \n" ;
    print "SNSIMPARAMS are @SN_SIM_PARAMS\n"; 
    print "params are @params\n";
    print "Using:    \n" ;
    print "\t JOBNAME_SIMSEDextract \t $JOBNAME_SIMSEDextract \n";
    print "\t GENMODEL \t--v $GENMODEL \n";
    print "\t tlist    \t--tlist $phaselist_file \n";
    print "\t params  \t\t $arglist2 $arglist3 \n";
    print "\t logfile  \t\t $log_file \n";
    print "\t output type \t $OUTPUT_FORMAT \n";
    print "\t peakmjd \t\t--T0 $peakmjd \n";
    print "\t seed \t\t--seed $RANSEED \n";
    print "\t galfrac \t\t--gfrac $GALFRAC \n";
    print "\t galfile \t\t--gfile $GALFILE \n";

    @output = `$JOBNAME_SIMSEDextract $arglist1 $arglist2 $arglist3 $arglist4 $arglist5 > $log_file`;

    print "\t $JOBNAME_SIMSEDextract $arglist1 $arglist2 $arglist3 $arglist4 $arglist5 \n";

    print "\tChecking for successful run\n";
    if ($? == -1 ) {
        $ABORT_die = "\nABORT: SIMSEDextract routine failed. \n";
        die $ABORT_die;
    }

    $KEY = "ABORT" ; ## OPTIONAL
    check_goodexit("SIMSEDextract", $log_file, $KEY);

}

sub increment_RANSEED()
{
    ## increment $RANSEED 
    $RANSEED = $RANSEED + $RANSEED_OFFSET; 
}

sub phase_generator($) {

    my $nphase = shift;

    my $random;
    my $nlist;
    my @phaselist;
    my $printphase;
    my $nphaselist;
    my $ilist;
    my $check_flag; 
    my $phaseflag;
    my $absval;
    my $d1PAR;
    my $d2PAR; 
    my $d3PAR; 
    my $d4PAR; 
    my $d5PAR;
    my $dist_cadence;

	### dist #    d1   #    d2     #   d3      #   d4     #   d5      #
	###      #         #           #           #          #           #
	### exp  #  lambda # t offset  # (cadence) #          #           #
	### pois #  lambda # (cadence) #           #          #           #
        ### gaus #   mean  #   stdev   # (cadence) #          #           #
	### gam1 #    k    #   theta   # (cadence) #          #           #
	### gam2 #    k1   #   theta1  #    k2     #  theta2  # (cadence) #
	### SNLS #  frac1  #  frac0    #  frac2    #          #           #

    print " IN PHASE_GENERATOR - nphase is $nphase \n";
	$d1PAR = $PHASE_DPAR1_HASH{$GENVERSION};
	$d2PAR = $PHASE_DPAR2_HASH{$GENVERSION};
	$d3PAR = $PHASE_DPAR3_HASH{$GENVERSION};
	$d4PAR = $PHASE_DPAR4_HASH{$GENVERSION};
	$d5PAR = $PHASE_DPAR5_HASH{$GENVERSION};

		if ($PHASE_DISTRIBUTION eq 'EXP') { $dist_cadence = $d3PAR ; }
		elsif ($PHASE_DISTRIBUTION eq 'POISSON' ) { $dist_cadence = $d2PAR; }
		elsif ($PHASE_DISTRIBUTION eq 'GAMMA' ) { $dist_cadence = $d3PAR; }
		elsif ($PHASE_DISTRIBUTION eq 'GAUSSIAN' ) { $dist_cadence = $d3PAR; }
                else { 
		    $dist_cadence = ( $IDEAL_EPOCHRANGE[1] - $IDEAL_EPOCHRANGE[0])/$nphase; 
		}

    print "\t Generating phase list. NPHASE is $nphase. dist_cadence is $dist_cadence\n";

    if ( $PHASE_DISTRIBUTION eq 'FLAT' or $PHASE_DISTRIBUTION eq '') {
	$printphase = rand($dist_cadence) + $IDEAL_EPOCHRANGE[0]; 
	while ( $printphase < @IDEAL_EPOCHRANGE[1] ){
	    push(@phaselist, $printphase);
	    $printphase = $printphase + $dist_cadence;
	}
    } else {
	@phaselist = ();
	$nlist = @phaselist;

	while ( $nlist < $nphase ){
	    $check_flag = 0;
	    ### repeat loop until check_flag equals 1
	    while ( $check_flag < 1 ) {
		### get a phase value
		if ($PHASE_DISTRIBUTION eq 'EXP') { 
		    ## print "debug: calling EXP rand \n"; 
		    $printphase = exp_rand($d1PAR)*exp($d1PAR*$d2PAR);
		    ## print "debug: printphase is $printphase \n"; 
		}
		elsif ($PHASE_DISTRIBUTION eq 'POISSON' ) { $printphase = poisson_rand($d1PAR); }
		elsif ($PHASE_DISTRIBUTION eq 'GAMMA' ) { $printphase = gamma_rand($d1PAR, $d2PAR); }
		elsif ($PHASE_DISTRIBUTION eq 'GAUSSIAN' ) { $printphase = gaussian_rand($d1PAR, $d2PAR); }
		## print "debug: nphase: $nphase, nlist: $nlist, phase: $printphase, @phaselist \n"; 
		### calculate phaseflag value for current printphase and phaselist values. 
		$ilist = 0;
		$phaseflag = 0;
		while ( $ilist < $nlist ) {
		    $absval = abs($printphase-$phaselist[$ilist]);
		    if ( $absval >= $dist_cadence ) {
			$phaseflag = $phaseflag + 1;
		    }
		    $ilist = $ilist + 1;
		}
		### check that phase value meets all requirements
		## print "debug: $IDEAL_EPOCHRANGE[0], $printphase, $IDEAL_EPOCHRANGE[1], $phaseflag, $nlist.  \n"; 
		if ( $printphase >= $IDEAL_EPOCHRANGE[0] &&
		     $printphase <= $IDEAL_EPOCHRANGE[1] &&
		     $phaseflag == $nlist ) {
		    	    push(@phaselist, $printphase);
			    $check_flag  = 1; 
			    $nlist = @phaselist; 
			    ## print "debug: good phase is $printphase \n";
			    ## print "debug: new list is @phaselist \n";
			}
	    }
	    $nlist = @phaselist;
	}
    }
    return ( $nphase, @phaselist );
}


sub var_assign(@) 
{
    my ($default, $input) = @_;
    my $output_val;

    if ($input eq '') {
	$output_val = $default;
    } else {
	$output_val = $input;
    }

    return $output_val;
}

sub parse_SPECDIST(@)
{

     my ($STYPE, $DTYPE, $nvar,  @distpar) = @_;
     my $check_flag;
     
     if ( $STYPE eq 'PHASE' ) {
	 $PHASE_DISTRIBUTION = $DTYPE;
     } elsif ( $STYPE eq 'EPOCH' ) {
	 $EPOCH_DISTRIBUTION = $DTYPE;
     } else {
	 die "ABORT: SPEC_DISTRIBUTION argument #1 invalid\n";
     }
		
     if ($DTYPE eq 'POISSON') {
	 $check_flag = poisson_assign(@distpar);
     } elsif ($DTYPE eq 'GAMMA') {
	 $check_flag = gamma_assign(@distpar);
     } elsif ($DTYPE eq 'GAMMA2') {
	 $check_flag = gamma2_assign(@distpar);
     } elsif ($DTYPE eq 'GAUSSIAN') {
	 $check_flag = gaussian_assign(@distpar);
     } elsif ($DTYPE eq 'EXP') {
	 $check_flag = exp_assign(@distpar);
     } elsif ($DTYPE eq 'SNLS') {
	 $check_flag = snls_assign(@distpar);
	 ## print "debug $GENVERSION $DTYPE var_array @VAR_ARRAY \n"; 
     } else {
	 die "ABORT: SPEC_DISTRIBUTION argument #2 invalid\n";
     }

     if ( $check_flag == 1 ) {
	 die "ABORT: wrong number of variables for spec distribution\n";
     }

     if ( $STYPE eq 'PHASE' ) {
	 ## print "debug $GENVERSION var_array @VAR_ARRAY \n";
	  @PHASEVAR_ARRAY = @VAR_ARRAY;
	 ## print "debug $GENVERSION phasevar_array @PHASEVAR_ARRAY \n";
     } elsif ( $STYPE eq 'EPOCH' ) {
	 ## print "debug $GENVERSION var_array @VAR_ARRAY \n";
	  @EPOCHVAR_ARRAY = @VAR_ARRAY;
	 ## print "debug $GENVERSION epochvar_array @EPOCHVAR_ARRAY \n";
      }
      
}


sub fillvar_SPECARG(@)
{
    my ($GENVERSION, @inputs) = @_;
    my $iarg = 0;
    my $skip_arg;
    my $die_flag;
    my $narg = @inputs;

    while ( $iarg  < $narg ) {
	if ( $inputs[$iarg] eq "IDEAL_EPOCHRANGE" ) {
	    $IDEAL_EPOCHRANGE[0] = $inputs[$iarg+1];
	    $IDEAL_EPOCHRANGE[1] = $inputs[$iarg+2];
	    $iarg = $iarg+3;
	} elsif ( $inputs[$iarg] eq "IDEAL_CADENCE" ) {
	    $IDEAL_CADENCE = $inputs[$iarg+1];
	    $iarg = $iarg+2;
	} elsif ( $inputs[$iarg] eq "IDEAL_SNR_OBS") {
	    $IDEAL_SNR_OBS[0] = $inputs[$iarg+1];
	    $IDEAL_SNR_OBS[1] = $inputs[$iarg+2];
	    $iarg = $iarg+3;
	} elsif ( $inputs[$iarg] eq "IDEAL_LAMINIT_REST" || $inputs[$iarg] eq "IDEAL_LAMINIT") {
	    $IDEAL_LAMINIT_REST = $inputs[$iarg+1];
	    $iarg = $iarg+2;
	} elsif ( $inputs[$iarg] eq "IDEAL_LAMFINAL_REST" || $inputs[$iarg] eq "IDEAL_LAMFINAL") {
	    $IDEAL_LAMFINAL_REST = $inputs[$iarg+1];
	    $iarg = $iarg+2;
	} elsif ( $inputs[$iarg] eq "IDEAL_DLAM") {
	    $IDEAL_DLAM = $inputs[$iarg+1];
	    $iarg = $iarg+2;
	} elsif ( $inputs[$iarg] eq "SPEC_DISTRIBUTION") {
	    $skip_arg = 3 + $inputs[$iarg+3];
	    parse_SPECDIST(@inputs[$iarg+1..$iarg+$skip_arg]);
	    $iarg = $iarg + $skip_arg;	
	} elsif ( $inputs[$iarg] eq "SPECGEN_LIB") {
	    $SPECGEN_LIB = $inputs[$iarg+1];
	    print "SPECGEN_LIB: $SPECGEN_LIB\n";
	    $iarg = $iarg+2;
	} elsif ( $inputs[$iarg] eq "GALFRAC"){
	    $GALFRAC = $inputs[$iarg+1];
	    $iarg = $iarg + 1;
	} elsif ( $inputs[$iarg] eq "GALFILE"){
	    $GALFILE = $inputs[$iarg+1]; 
	    $iarg = $iarg + 1;
	} else {
	    $iarg = $iarg + 1;	
	}
    }

    if ( $GALFRAC eq '' ) {
	$GALFRAC = 0.0;
	$GALFILE = "NULL"; 
    } else {
    }

    if ( $IDEAL_SNR_OBS[0] eq '' ) {
	$IDEAL_SNR_OBS[0] = 20;
	$IDEAL_SNR_OBS[1] = 100;
    } else {
    }

    if ( ( $IDEAL_CADENCE eq '' ) && ( $PHASE_DISTRIBUTION eq '' ) ) {
	$IDEAL_CADENCE = 5;
    } else {
    }

    if ( $IDEAL_DLAM eq '' ) {
	$IDEAL_DLAM = 0.0;
	### print "\t USING DEFAULT IDEAL_DLAM \t\t $IDEAL_DLAM \n";
    } else {
	### print "\t IDEAL_DLAM \t $IDEAL_DLAM \n";
    }

    if ( $IDEAL_LAMINIT_REST eq '' ) {
	$IDEAL_LAMINIT_REST = 2000;
    } else {
    }

    if ( $IDEAL_LAMFINAL_REST eq '' ) {
	$IDEAL_LAMFINAL_REST = 9200;
    } else {
    }


	    $EPOCHRANGE_LOW_HASH{$GENVERSION} = $IDEAL_EPOCHRANGE[0];
	    $EPOCHRANGE_HIGH_HASH{$GENVERSION} = $IDEAL_EPOCHRANGE[1];
	    $CADENCE_HASH{$GENVERSION} = $IDEAL_CADENCE;
	    $SNR_HASH{$GENVERSION} = $IDEAL_SNR_OBS[0];
            $SNRLAMBIN_HASH{$GENVERSION} = $IDEAL_SNR_OBS[1];
	    $LAMINIT_HASH{$GENVERSION} = $IDEAL_LAMINIT_REST;
	    $LAMFINAL_HASH{$GENVERSION} = $IDEAL_LAMFINAL_REST;
	    $DLAM_HASH{$GENVERSION} = $IDEAL_DLAM;
            $GENLIB_HASH{$GENVERSION} = $SPECGEN_LIB;
            $GALFRAC_HASH{$GENVERSION} = $GALFRAC; 
            $GALFILE_HASH{$GENVERSION} = $GALFILE; 
    


}


sub fill_EPOCHvar_SPECSIM(@) 
{
    my ($GENVERSION, $DPAR1, $DPAR2, $DPAR3, $DPAR4, $DPAR5, $DPAR6) = @_;

    $EPOCH_DPAR1_HASH{$GENVERSION} = $DPAR1;
    $EPOCH_DPAR2_HASH{$GENVERSION} = $DPAR2;
    $EPOCH_DPAR3_HASH{$GENVERSION} = $DPAR3;
    $EPOCH_DPAR4_HASH{$GENVERSION} = $DPAR4;
    $EPOCH_DPAR5_HASH{$GENVERSION} = $DPAR5;
    $EPOCH_DPAR6_HASH{$GENVERSION} = $DPAR6;
 
}

sub read_EPOCHvar_SPECSIM($)
{
    my $GENVERSION = shift;

    return ( $EPOCH_DPAR1_HASH{$GENVERSION}, $EPOCH_DPAR2_HASH{$GENVERSION},
	     $EPOCH_DPAR3_HASH{$GENVERSION}, $EPOCH_DPAR4_HASH{$GENVERSION},
	     $EPOCH_DPAR5_HASH{$GENVERSION}, $EPOCH_DPAR6_HASH{$GENVERSION} );
}

sub fill_PHASEvar_SPECSIM(@) 
{
    my ($GENVERSION, $DPAR1, $DPAR2, $DPAR3, $DPAR4, $DPAR5, $DPAR6) = @_;

    $PHASE_DPAR1_HASH{$GENVERSION} = $DPAR1;
    $PHASE_DPAR2_HASH{$GENVERSION} = $DPAR2;
    $PHASE_DPAR3_HASH{$GENVERSION} = $DPAR3;
    $PHASE_DPAR4_HASH{$GENVERSION} = $DPAR4;
    $PHASE_DPAR5_HASH{$GENVERSION} = $DPAR5;
    $PHASE_DPAR6_HASH{$GENVERSION} = $DPAR6;
 
}

sub read_PHASEvar_SPECSIM($)
{
    my $GENVERSION = shift;

    return ( $PHASE_DPAR1_HASH{$GENVERSION}, $PHASE_DPAR2_HASH{$GENVERSION},
	     $PHASE_DPAR3_HASH{$GENVERSION}, $PHASE_DPAR4_HASH{$GENVERSION},
	     $PHASE_DPAR5_HASH{$GENVERSION}, $PHASE_DPAR6_HASH{$GENVERSION} );
}

sub readvar_SPECARG($)
{
    my $GENVERSION = shift;

    return ( $EPOCHRANGE_LOW_HASH{$GENVERSION}, $EPOCHRANGE_HIGH_HASH{$GENVERSION}, 
	     $CADENCE_HASH{$GENVERSION}, $SNR_HASH{$GENVERSION}, $SNRLAMBIN_HASH{$GENVERSION},
	     $LAMINIT_HASH{$GENVERSION}, $LAMFINAL_HASH{$GENVERSION}, $DLAM_HASH{$GENVERSION} );
}


sub initvar_SPECSIM()
{
    ## $GAUSS_MEAN = $GAUSS_STDEV = $CADENCE = '';
    ## $DIST_PAR1 = $GK1 = $GK2 = $GTHETA1 = $GTHETA2 = '';

}

sub init_INPUTKEYS()
{
    $GENVERSION = ''; 
    $SPECGEN_LIB = '';
    $EPOCH_DISTRIBUTION = '';
    $PHASE_DISTRIBUTION = '';
    @EPOCHVAR_ARRAY = '';
    @PHASEVAR_ARRAY = '';
    @VAR_ARRAY = '';
    @IDEAL_SPECARGS = '';

}

sub poisson_assign(@)
{
    
    my (@temp) = @_;
    my ($NVAR) = @temp;
    my $POISSON_LAM;
    my $CADENCE;

    if ( $NVAR >= 1 || $NVAR == 0) {
	$POISSON_LAM = var_assign($DEFAULT_POISSON, @temp[0]);
	$CADENCE = var_assign($DEFAULT_GAUSSCADENCE, $temp[1]);
	@VAR_ARRAY = ( $POISSON_LAM, $CADENCE );
	print "\t POISSON VAR: @VAR_ARRAY \n";
	return 0;
    } else {
	return 1;
    }
}

sub gamma_assign(@)
{

    my (@temp) = @_;
    my ($NVAR) = @temp;
    my $GK1;
    my $GTHETA1;

    if ( $NVAR == 2 || $NVAR == 0 ) {
	$GK1 = var_assign($DEFAULT_GK1, $temp[0]);
	$GTHETA1 = var_assign($DEFAULT_GTHETA1, $temp[1]);
	print "\t GAMMA VARS: @VAR_ARRAY \n";
	@VAR_ARRAY = ($GK1, $GTHETA1);
	return 0;
    } else {
	return 1;
    }
}

sub gamma2_assign(@)
{

    my (@temp) = @_;
    my ($NVAR) = @temp;
    my $GK1;
    my $GTHETA1;
    my $GK2;
    my $GTHETA2; 


    if ( $NVAR == 4 || $NVAR == 0 ) {
	$GK1 = var_assign($DEFAULT_GK1, $temp[0]);
	$GTHETA1 = var_assign($DEFAULT_GTHETA1, $temp[1]);
	$GK2 = var_assign($DEFAULT_GK2, $temp[2]);
	$GTHETA2 = var_assign($DEFAULT_GTHETA2, $temp[3]);
	@VAR_ARRAY = ($GK1, $GTHETA1,  $GK2, $GTHETA2 );
	print "\t GAMMA2 VARS: @VAR_ARRAY \n";
	return 0;
    } else {
	return 1;
    }
}

sub gaussian_assign(@)
{

    my (@temp) = @_;
    my ($NVAR) = @temp;
    my $GAUSS_MEAN;
    my $GAUSS_STDEV;
    my $CADENCE; ## OPTIONAL. DEF iS 1.


    if ( $NVAR >= 2 || $NVAR == 0 ) {
	$GAUSS_MEAN = var_assign($DEFAULT_GAUSSMEAN, $temp[0]);
	$GAUSS_STDEV = var_assign($DEFAULT_GAUSSSTDEV, $temp[1]);
	$CADENCE = var_assign($DEFAULT_GAUSSCADENCE, $temp[2]);
	@VAR_ARRAY = ( $GAUSS_MEAN, $GAUSS_STDEV, $CADENCE );
	### print "\t GAUSSIAN VARS: @VAR_ARRAY \n";
	return 0;
    } else {
	return 1;
    }
}
	
sub exp_assign(@)
{

    my (@temp) = @_;
    my ($NVAR) = @temp;
    my $EXP_LAMBDA; 
    my $EXP_DELTAt; ### DEFAULT is -35. 

    if ( $NVAR >= 2 || $NVAR == 0) {
	$EXP_LAMBDA = var_assign($DEFAULT_EXPLAM, $temp[0]);
	$EXP_DELTAt = var_assign($DEFAULT_EXPOFFSET, $temp[1]);
	$CADENCE = var_assign($DEFAULT_GAUSSCADENCE, $temp[2]);
	@VAR_ARRAY = ( $EXP_LAMBDA, $EXP_DELTAt, $CADENCE );
	### print "\t EXP VAR: @VAR_ARRAY \n";
	return 0;
    } else {
	return 1;
    }
}

sub snls_assign(@)
{

    my (@temp) = @_;
    my ($NVAR) = @temp;
    my $FRAC1; # default 0.85
    my $FRAC0; # default 0.10
    my $FRAC2; # default 0.05

    if ( $NVAR == 3 || $NVAR == 0) {
	$FRAC1 = var_assign($DEFAULT_FRAC1, $temp[0]);
	$FRAC0 = var_assign($DEFAULT_FRAC0, $temp[1]);
	$FRAC2 = var_assign($DEFAULT_FRAC2, $temp[2]);
	@VAR_ARRAY = ( $FRAC1, $FRAC0, $FRAC2 );
	### print "\t SNLS VAR: @VAR_ARRAY \n";
	return 0;
    } else {
	return 1;
    }
}

sub write_GENVERSIONhash()
{
    ### Make LCSIM_HASH, SPECSIM_HASH Entries Now ###
    ### Use Genversion as key - these should be unique ###
    $LCSIM_HASH{$GENVERSION} = $LCSIM_CALL;
    $SIMMINZ_HASH{$GENVERSION} = $SIM_ZMIN;
    $SIMMAXZ_HASH{$GENVERSION} = $SIM_ZMAX;
    fillvar_SPECARG($GENVERSION, @IDEAL_SPECARGS);
    $EPOCHDIST_HASH{$GENVERSION} = $EPOCH_DISTRIBUTION;
    $PHASEDIST_HASH{$GENVERSION} = $PHASE_DISTRIBUTION;
    ## print "debug write_GENVERSIONhash $GENVERSION @EPOCHVAR_ARRAY and @PHASEVAR_ARRAY\n";
    fill_EPOCHvar_SPECSIM($GENVERSION, @EPOCHVAR_ARRAY);    
    fill_PHASEvar_SPECSIM($GENVERSION, @PHASEVAR_ARRAY);    

}

sub get_nepochlist(@)
{
    my ($nSNE, $dist) = @_;
    my $iSNE;
    my $flat_nphase; 
    my $d1PAR; 
    my $d2PAR;
    my $d3PAR;
    my $d4PAR;
    my $d5PAR;
    my $printepochs;

	### dist #    d1   #    d2     #   d3      #   d4     #   d5      #
	###      #         #           #           #          #           #
	### exp  #  lambda # t offset  # (cadence) #          #           #
	### pois #  lambda # (cadence) #           #          #           #
	### gam1 #    k    #   theta   # (cadence) #          #           #
	### gam2 #    k1   #   theta1  #    k2     #  theta2  # (cadence) #
	### SNLS #  frac1  #  frac0    #  frac2    #          #           #
        ### gaus #   mean  #   stdev   # (cadence) #          #           #


	$d1PAR = $EPOCH_DPAR1_HASH{$GENVERSION};
	$d2PAR = $EPOCH_DPAR2_HASH{$GENVERSION};
	$d3PAR = $EPOCH_DPAR3_HASH{$GENVERSION};
	$d4PAR = $EPOCH_DPAR4_HASH{$GENVERSION};
	$d5PAR = $EPOCH_DPAR5_HASH{$GENVERSION};

    @NEPOCH_LIST = '';

        if ( $EPOCH_DISTRIBUTION eq 'EXP' ) {
	    for ( $iSNE=0; $iSNE<$nSNE; $iSNE++){ 
	        $printepochs = exp_rand($d1PAR);
		## print "EXP NEPOCHS: $printepochs \n"; 
		push(@NEPOCH_LIST, $printepochs); 
	    }
	} elsif ( $EPOCH_DISTRIBUTION eq 'POISSON' ) {
	    for ( $iSNE=0; $iSNE<$nSNE; $iSNE++){ 
	        $printepochs = poisson_rand($d1PAR);
		## print "POISSON NEPOCHS: $printepochs \n"; 
		push(@NEPOCH_LIST, $printepochs); 
	    }
	} elsif ( $EPOCH_DISTRIBUTION eq 'GAUSSIAN' ) {
	    for ( $iSNE=0; $iSNE<$nSNE; $iSNE++){ 
	        $printepochs = gaussian_rand($d1PAR, $d2PAR);
		## print "GAUSS NEPOCHS: $printepochs \n"; 
		push(@NEPOCH_LIST, $printepochs); 
	    }
	} elsif ( $EPOCH_DISTRIBUTION eq 'GAMMA' ) {
	    for ( $iSNE=0; $iSNE<$nSNE; $iSNE++){ 
	        $printepochs = gamma_rand($d1PAR, $d2PAR);
		## print "GAMMA NEPOCHS: $printepochs \n"; 
		push(@NEPOCH_LIST, $printepochs); 
	    }
	} elsif ( $EPOCH_DISTRIBUTION eq 'SNLS' ) { 
	    for ( $iSNE=0; $iSNE<$nSNE; $iSNE++){ 
	        $printepochs = SNLS_rand($d1PAR, $d2PAR, $d3PAR);
		## print "SNLS NEPOCHS: $printepochs \n"; 
		push(@NEPOCH_LIST, $printepochs); 
		## print "SNLS pars: d1 $d1PAR d2 $d2PAR d3 $d3PAR \n";
		print "get_nepochlist: @NEPOCH_LIST \n" ;
	    }
	} elsif ( $EPOCH_DISTRIBUTION eq 'GAMMA2' ) { 
	    for ( $iSNE=0; $iSNE<$nSNE; $iSNE++){ 
	        $printepochs = gamma2_rand($d1PAR, $d2PAR, $d3PAR, $d4PAR);
		## print "GAMMA2 NEPOCHS: $printepochs \n"; 
		push(@NEPOCH_LIST, $printepochs); 
	    }
	} elsif ( $EPOCH_DISTRIBUTION eq 'FLAT' or $EPOCH_DISTRIBUTION eq '' ) {
	    $flat_nphase = ( $IDEAL_EPOCHRANGE[1] - $IDEAL_EPOCHRANGE[0] ) / $IDEAL_CADENCE;
	    for ( $iSNE=0; $iSNE<$nSNE; $iSNE++){ 
		push(@NEPOCH_LIST, $flat_nphase) ;
	    }
	}
	

}

sub SNLS_rand(@)
{
    my ( $frac1, $frac2, $frac3) = @_;
    my $random; 
    my $nphase;

	$random = rand();
	if ( $random < $frac1 ) {
	    $nphase = 1;
	} elsif ( $random >= $frac1 && $random < ($frac1+$frac2) ) {
	    $nphase = 0;
	} else {
	    $nphase = 2;
	}


    return $nphase;
}

sub gamma2_rand(@)
{
 
    my ( $GK1, $GTHETA1, $GK2, $GTHETA2 ) = @_;
    my $random;
    my $nphase;
    
    	$random = rand();
	if ( $random < (1.5/8.0) ){
	    $nphase = int(gamma_rand($GK1, $GTHETA1));
	} elsif ($random <(2.5/8.0) && $random >= (1.5/8.0)){
	    $nphase = int(gamma_rand($GK2, $GTHETA2));
	} else {
	    $nphase = 0;
	}

    return $nphase; 

}
