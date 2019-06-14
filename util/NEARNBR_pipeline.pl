#!/usr/bin/perl
# 
# Created Dec 2014 by R.Kessler
# Run nearest neighbor (NN) pipeline:
#
#  * two simjobs using sim_SNmix
#  * split_and_fit job on REF
#  * split_and_fit job on TRAIN for each SEPMAX_VARDEF grid
#  * nearnbr_FoM on each fit-output
#  * final split_and_fit on data using NN-trained parameters
#
# There are three kinds of data that are analyzed in the final stage:
#
#  1. REALDATA    defined by VERSION_REALDATA key
#
#  2. SIMDATA     defined by SIMGEN_MASTER_DATA key 
#                   (SN models can differ from training models) 
#
#  3. VALIDATA    internally generated with same models used in training;
#                 intended for validation check.
#
#
# Usage:
#   NEARNBR_pipeline.pl <inputFile>
#
#   NEARNBR_pipeline.pl <inputFile>  -RESTART <istage>
#   NEARNBR_pipeline.pl <inputFile>  -END     <istage>
#   NEARNBR_pipeline.pl <inputFile>  -ONLY   <istage>
#      (see $ISTAGE_XXX params below)
#   NEARNBR_pipeline.pl <inputFile>  -RESTART <istage> NOLAUNCH
#   NEARNBR_pipeline.pl <inputFile>  NOSUBMIT
#      (create NNSTAGE input file, then abort)
#
#
# NGEN_UNIT in SIMGEN_MASTER file defines stats for REF and TRAINING.
# BIASCOR sample is x8 bigger
# DATA sample is x10 smaller
#
# -----------------
# inputFile keys:
# 
#  TOPDIR_OUTPUT:         <topdir location for all output>
#  INFILE_SIMDATA_MASTER:  <sim_SNmix master file for SIMDATA>
#  INFILE_SIMTRAIN_MASTER: <sim_SNmi master template for training & BIASCOR>
#  INFILE_FIT_MASTER:      <file for snlc_fit.exe>
#  WFALSE:                <Wfalse used for pseudo-purity: default=1>
#  NEARNBR_SEPMAX_VARDEF: <sepmax grid>
#  NEARNBR_SEPMAX_VARDEF: <another sepmax grid (optional)>
#  VERSION_REALDATA:      <photometry data version>

#  
# NGEN_UNIT in SIMGEN-master file is for training.
# Here is how to specifiy NGEN_UNIT for data and biasCor:
#  SCALE_NGEN_UNIT_VALIDATA: NGEN_UNIT(VALIDATA) = SCALE x NGEN_UNIT(TRAIN)
#  SCALE_NGEN_UNIT_BIASCOR:  NGEN_UNIT(BIASCOR)  = SCALE x NGEN_UNIT(TRAIN)
#     or
#  NGEN_UNIT_VALIDATA:    <absolute NGEN_UNIT for vali-data>
#  NGEN_UNIT_BIASCOR:    <absolute NGEN_UNIT for biasCor>
#    (note that absolute NGEN_UNIT overrides SCALE_NGEN_UNIT_XXX)
#
# Beware that LCFIT stages may contain symbolic links, so to copy use
# rsync -rl <outDir>  instead of scp 
#
#
#          HISTORY
#
# Sep 9 2015: add WFALSE input key (default=1)
#
# Feb 29 2016: parse updated SALT2mu file that now has NSNFIT key
#
# Apr 12 2016: new optional key  "SYMLINKS_FLAG: 0" to turn off
#              the creation of NNFIT_DATA and NNFIT_SIM
#
# Jun 14 2016: allow VERSION_REALDATA to be data or sim.
# Oct 06 2016: allow VERSION_REALDATA to be in PRIVATE_DATA_PATH
#
# Jan 24 2017: 
#    + set DO_BBC if varname_PIa is defined in SALT2mu_infile.
#      for pure photo-z analysis, BBC fit doesn't work, so don't bother. 
#      Note that stages 5 & 6 are skipped without BBC option.
#
# July 11 2017: check RANSEED and RANSEED_REPEAT
# July 24 2017: check TOPDIR_OUTPUT for ENV
# July 06 2018: check PATH_SNDATA_SIM.LIST file
# Jan  23 2019: parsing sigint and scalePCC is optional, not requred.
#
# Feb  05 2019: 
#  + refactor and break SIMGEN stage into 3 stages: DATA,TRAIN,BIASCOR
#  + SIM-DATA is now a separate stage; no longer take data as 
#     pre-scaled sub-sample from training. This change needed
#     to enable different SN models for DATA and TRAINING.
#  + new keys to specify NGEN_UNIT for data and biasCor
#      (see instructions above)
#
# Feb 18 2019: 
#  major overhaul to include SIMDATA with SN models that can differ
#  from SN models used in training.
#
# Apr 11 2019: --truetype --> -truetype (just one dash)
#
# Jun 14 2019: 
#   + read new optional key SIMTRAIN_SCALE_NON1A to enhance CC
#     without increasing SNIa size. See GENOPT_GLOBAL appended.
#     NNtrain C-code has been updated to account for SCALE_NON1A > 1.
#
# -------------------------------------

use List::Util qw(first);
use IO::Handle ;
use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# hard-wired definitions

my $ISTAGE_SIMGEN_SIMDATA  =  1 ; # gen sim-data w. arbitrary models
my $ISTAGE_SIMGEN_VALIDATA =  2 ; # gen validation data w. same model as train
my $ISTAGE_SIMGEN_REF      =  3 ; # gen sim-ref for training
my $ISTAGE_SIMGEN_TRAIN    =  4 ; # gen sim-training set
my $ISTAGE_SIMGEN_BIASCOR  =  5 ; # generate sim biasCor
my $ISTAGE_SIMFIT_REF      =  6 ; # fit sim reference
my $ISTAGE_SIMFIT_TRAIN    =  7 ; # train with SEPMAX grid
my $ISTAGE_NNFOM           =  8 ; # run nearnbr_maxFoM to get SEPMAX values
my $ISTAGE_SIMFIT_BIASCOR  =  9 ; # use trained SEPMAX for biasCor 
my $ISTAGE_SIMFIT_CCPRIOR  = 10 ; # use trained SEPMAX for CCprior
my $ISTAGE_FINALFIT        = 11 ; # use trained SEPMAX on data and sim
my $ISTAGE_BBC_SUMMARY     = 12 ; # make BBC summary file
my $NSTAGE                 = 12 ;
my @DOSTAGE ;
my $ISTAGE_GLOBAL ;  # current stage 

my $q  = "\'" ;
my $qq = '"' ;

my $OPT_ABORT = 0 ; # parsing option
my $OPT_WARN  = 1 ; # 1=> leave warning message; 2=>quiet
my $OPT_QUIET = 2 ; # no warnong if key not found

my $CURRENT_DIR = `pwd` ; $CURRENT_DIR  =~ s/\s+$// ;   # trim trailing space
my $MOI = `whoami`      ; 
my $MOI4 = substr($MOI,0,4);
my $GLOBAL_PREFIX  = "NN_${MOI4}" ;

my @SIMGEN_SUFFIX = ( "NULL", "SIMDATA", "VALIDATA", "REF", 
		      "TRAIN", "BIASCOR" );
my $SIMGEN_SUMMARY_FILE = "SIMJOB_SUMMARY.LOG" ;

my $SNDATA_ROOT  = $ENV{'SNDATA_ROOT'};
my $ROOT_DIR     = $ENV{'ROOT_DIR'};

my $dashLine =  
    "---------------------------------------------" .
    "-------------------------------";

my $SYMLINK_REALDATA  = "NNFIT_REALDATA" ;
my $SYMLINK_VALIDATA  = "NNFIT_VALIDATA" ; # was NFIT_SIM before
my $SYMLINK_SIMDATA   = "NNFIT_SIMDATA" ;  # user models

# code/script names used below
my $SCRIPT_FIT   = "split_and_fit.pl" ;
my $SCRIPT_SIM   = "sim_SNmix.pl" ;
my $SCRIPT_WAIT  = "wait_for_files.pl" ;
my $BINARY_NNFOM = "nearnbr_maxFoM.exe" ;

my $SCRIPT_FITARGS = "NOPROMPT 1";  # always these args for split_and_fit

# declare inputs
my ($INFILE_NEARNBR_MASTER);
my ($TOPDIR_OUTPUT, $INFILE_SIMTRAIN_MASTER, $INFILE_SIMDATA_MASTER);
my ($INFILE_FIT_MASTER);
my ($INFILE_SIMGEN_Ia_DATA,  $INFILE_SIMGEN_NONIa_DATA,  $NGEN_UNIT_DATA );
my ($INFILE_SIMGEN_Ia_TRAIN, $INFILE_SIMGEN_NONIa_TRAIN, $NGEN_UNIT_TRAIN );
my ($ISTAGE_START, $ISTAGE_END, $VERSION_REALDATA );
my ($SYMLINKS_FLAG, $DO_LAUNCH, @CONTENTS_SIMTRAIN, @CONTENTS_SIMDATA) ;
my ($SCALE_NGEN_UNIT_VALIDATA, $SCALE_NGEN_UNIT_BIASCOR);
my ($NGEN_UNIT_VALIDATA, $NGEN_UNIT_BIASCOR, $SIMTRAIN_SCALE_NON1A);

my ($INPUT_WFALSE);
my $WFALSE_VALID = "1 3 5";

my (@GENV_LIST_SIMTRAIN, @GENV_LIST_SIMDATA, @GENV_LIST_NN );
my ($GENPREFIX_ORIG_TRAIN, $GENPREFIX_ORIG_DATA );
my (@SIMGEN_GENPREFIX, $GENOPT_GLOBAL_TRAIN );
my (@RANSEED, @SEPMAX_VARDEF_LIST);
my ($NTRAIN, $NFITOPT, $NGENVER_SIMTRAIN, $NGENVER_SIMDATA);
my ($DO_SALT2mu, $DO_BBC );
my ($PRIVATE_DATA_PATH);

# -- - 

my (@CLONE_INFILE_SIMGEN_MASTER );
my (@MSGERR, $SIMGEN_ALLDONE_FILE);
my (@SIMVER_PREFIX, $OUTDIR_REF, @OUTDIR_TRAIN);
my (@BIASCOR_PATH, @CCPRIOR_PATH  );
my ($TRAINFILE_PATH_MAPFILE, @FITOPT_STRING , @NNLOG_LIST );
my ($ARCDIR_INPUTS, $ARCDIR_STDOUT);
my ($TOPDIR_SIMLOGS, $NFITSIM_REF );

# functions
sub setDefaults ;
sub parse_args ;
sub parse_inFile_master ;
sub checkFiles ;
sub prep_outDir ;

sub NNstagePrefix ;
sub run_simgen ;
sub parse_SIMDATA_MASTER ;
sub parse_SIMTRAIN_MASTER ;
sub parse_SNFIT_MASTER ;
sub parse_SIMGEN_RANSEED ;
sub make_SIMGEN_MASTER ;
sub submit_SIMGEN_MASTER ;
sub print_DOSTAGE ;

sub run_lcfit ;
sub run_nearnbr_maxFoM ;
sub update_NNFOM_summary ; 
sub update_BBC_summary ;
sub BBC_summary ;
sub get_SEPMAX_VARDEF ;
sub get_NGEN_UNIT ;
sub print_indexLegend ;

# ============ BEGIN MAIN ===============

$ISTAGE_GLOBAL = -1 ;

&setDefaults();

&parse_args() ;

&parse_inFile_master();

&prep_outDir();

&parse_SIMDATA_MASTER();
&parse_SIMTRAIN_MASTER();
&parse_SNFIT_MASTER();  # parse fit-input file

&checkFiles();

if ( $DO_LAUNCH == 0 ) {
    die "\n# =============================================== \n" .
	"\t USER-ABORT JOB-LAUNCH. \n" .
	"# =============================================== \n"
	;
}

$ISTAGE_GLOBAL = $ISTAGE_SIMGEN_SIMDATA ;
&run_simgen($ISTAGE_GLOBAL);

$ISTAGE_GLOBAL = $ISTAGE_SIMGEN_VALIDATA ;
&run_simgen($ISTAGE_GLOBAL);

$ISTAGE_GLOBAL = $ISTAGE_SIMGEN_REF ;
&run_simgen($ISTAGE_GLOBAL);

$ISTAGE_GLOBAL = $ISTAGE_SIMGEN_TRAIN ;
&run_simgen($ISTAGE_GLOBAL);

if ( $DO_BBC ) {
    $ISTAGE_GLOBAL = $ISTAGE_SIMGEN_BIASCOR ;
    &run_simgen($ISTAGE_GLOBAL);
}
else {
    print "\n SKIP BIASCOR STAGE $ISTAGE_SIMGEN_BIASCOR \n";
}

# fit the REF and TRAIN, but NOT the BiasCor sample
$ISTAGE_GLOBAL = $ISTAGE_SIMFIT_REF ;
&run_lcfit(-1, 0, "NULL");  # only the -1 arg is used

my ($itrain, $iver);

# train each sim using the same NN grid
$ISTAGE_GLOBAL = $ISTAGE_SIMFIT_TRAIN ;
for ($itrain=0; $itrain < $NTRAIN; $itrain++ ) { 
    &run_lcfit($itrain, -9, "TRAIN");   # -9 arg ignored
}


# ------------------------------------------------------------
# read large histograms and search for optimal SEPMAX params
$ISTAGE_GLOBAL = $ISTAGE_NNFOM ;
for ($itrain=0; $itrain < $NTRAIN; $itrain++ ) {  
    &run_nearnbr_maxFoM($itrain);
}


# ------------------------------------------------------------
# Apr 2016: run fits on BIASCOR sample that has only SNIa
if ( $DO_BBC ) {
    $ISTAGE_GLOBAL = $ISTAGE_SIMFIT_BIASCOR ;
    for ($itrain=0; $itrain < $NTRAIN; $itrain++ ) {  
	for($iver=0; $iver < $NGENVER_SIMTRAIN; $iver++ ) {
	    &run_lcfit($itrain, $iver, "BIASCOR"); # all args used
	}
    }
}

# ------------------------------------------------------------
# Apr 2016: run fits on CCPRIOR sample; process only the CC
if ( $DO_BBC ) {
    $ISTAGE_GLOBAL = $ISTAGE_SIMFIT_CCPRIOR ;
    for ($itrain=0; $itrain < $NTRAIN; $itrain++ ) {  
	for($iver=0; $iver < $NGENVER_SIMTRAIN; $iver++ ) {
	    &run_lcfit($itrain, $iver, "CCPRIOR"); # all args used
	}
    }
}


# ------------------------------------------------------------
# run final data+sim fit on each sim since the trained SEPMAX
# values depend on the sim.
$ISTAGE_GLOBAL = $ISTAGE_FINALFIT ;
for ($itrain=0; $itrain < $NTRAIN; $itrain++ ) {  
    for($iver=0; $iver < $NGENVER_SIMTRAIN; $iver++ ) {
	&run_lcfit($itrain, $iver,"FINAL"); # all args used
    }
}


# ------------------------------------------
if ( $DO_SALT2mu && $DO_BBC ) {
    $ISTAGE_GLOBAL = $ISTAGE_BBC_SUMMARY ;
    &update_BBC_summary(-1,-1,-1,-1,-1);  #open file, write table header
    for ($itrain=0; $itrain < $NTRAIN; $itrain++ ) {  
	for($iver=0; $iver < $NGENVER_SIMTRAIN; $iver++ ) {
	    &BBC_summary($itrain, $iver); 
	}
    }
    &update_BBC_summary(-7,-7,-7,-7,-7);  # print legend and close file
}


print "\n ================================================ \n";
print " Done with all NNSTAGEs. \n" ;

# ========== END MAIN ===================


sub setDefaults {

    my ($i);

    $ISTAGE_START = 1 ;
    $ISTAGE_END   = $NSTAGE ;

    $INFILE_SIMTRAIN_MASTER = "" ;
    $INFILE_SIMDATA_MASTER  = "" ;
    $INFILE_FIT_MASTER      = "" ;

# NGEN_UNIT in SIMGEN-master file is for training
    $SCALE_NGEN_UNIT_VALIDATA = 0.10 ;
    $SCALE_NGEN_UNIT_BIASCOR  = 8.0  ;
    $NGEN_UNIT_VALIDATA       = 0.0  ; # non-zero -> override SCALE
    $NGEN_UNIT_BIASCOR        = 0.0  ; # idem

    for($i=0; $i <= $NSTAGE ; $i++ ) 
    { $DOSTAGE[$i] = 0; }

    $INPUT_WFALSE =  1 ;

    $SYMLINKS_FLAG = 1 ;
    
    $DO_SALT2mu = $DO_BBC = 0 ;

    $DO_LAUNCH =  1;

} # end of sub setDefaults 

# =============================================
sub parse_args {

    my ($NARG, $arg, $nextArg, $i, $i1, @USE_ARG ) ;

    $NARG = scalar(@ARGV);
    if ( $NARG < 1 ) {
	$MSGERR[0] = "Must give input filename as argument";
	sntools::FATAL_ERROR(@MSGERR);
    }

    for($i=0; $i<$NARG; $i++ ) { $USE_ARG[$i]=0; }

    $INFILE_NEARNBR_MASTER = $ARGV[0] ;
    $USE_ARG[0] = 1;

    for ( $i = 1; $i < $NARG ; $i++ ) {
	$i1 = $i + 1;

        $arg     = $ARGV[$i] ;
	$nextArg = $ARGV[$i1] ;

	if ( $arg eq "-RESTART" ) {
	    $ISTAGE_START = $nextArg ;
	    $USE_ARG[$i]   = 1;
	    $USE_ARG[$i1]  = 1;
	}

	if ( $arg eq "-END" ) {
	    $ISTAGE_END = $nextArg ;
	    $USE_ARG[$i]   = 1;
	    $USE_ARG[$i1]  = 1;
	}

	if ( $arg eq "-ONLY" ) {
	    $ISTAGE_START = $nextArg ;
	    $ISTAGE_END   = $nextArg ;
	    $USE_ARG[$i]   = 1;
	    $USE_ARG[$i1]  = 1;
	}

	if ( $arg eq "NOLAUNCH" )  { $DO_LAUNCH = 0 ; $USE_ARG[$i]=1; }
	if ( $arg eq "NOSUBMIT" )  { $DO_LAUNCH = 0 ; $USE_ARG[$i]=1; }
    }

    # check for unused args.
    my $NARG_UNKNOWN = 0 ;
    for($i=0; $i<$NARG; $i++ ) { 
	if ( $USE_ARG[$i] == 0 ) {
	    print "\t ERROR: unknown input arg  $ARGV[$i] \n";
	    $NARG_UNKNOWN++ ;
	}
    }

    if ( $NARG_UNKNOWN > 0 ) {
	$MSGERR[0] = "$NARG_UNKNOWN unknown input arg(s)";
	$MSGERR[1] = "Scroll up to see list." ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    # ------------------------------
    # set DOSTAGE flags
    my $i;
    for($i=$ISTAGE_START; $i <= $ISTAGE_END; $i++ ) 
    { $DOSTAGE[$i] = 1 ; }

    print " Process STAGE $ISTAGE_START to $ISTAGE_END \n";
    &print_DOSTAGE("parse_args");

    return ;

}  # end of parse_args

sub print_DOSTAGE {
    my($CALLFUN) = @_ ;

    my($i);
    print " $CALLFUN: DOSTAGE = ";
    for($i=1; $i <= $NSTAGE; $i++ )  { print "$DOSTAGE[$i] " ; }
    print "\n";
    $| = 1; # flush stdout

    return ;

}  # end print_DOSTAGE

# ===================================
sub parse_inFile_master {

    my ($KEY, @CONTENTS_MASTER, @tmp, $i );

    print " Parse $INFILE_NEARNBR_MASTER \n";

    # scoop contents of input file
    @CONTENTS_MASTER = ();
    sntools::loadArray_fromFile($INFILE_NEARNBR_MASTER, \@CONTENTS_MASTER);
    
    $KEY  = "TOPDIR_OUTPUT:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_ABORT, @CONTENTS_MASTER);
    $TOPDIR_OUTPUT = "$tmp[0]" ;
    # check for ENV  
    $TOPDIR_OUTPUT  = qx(echo $TOPDIR_OUTPUT);
    $TOPDIR_OUTPUT  =~ s/\s+$// ;   # trim trailing whitespace

    print "\t TOPDIR_OUTPUT:  $TOPDIR_OUTPUT \n";

    # ---------------
    $KEY  = "INFILE_SIMGEN_MASTER:" ;  # legacy key 
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0 ) { $INFILE_SIMTRAIN_MASTER = "$tmp[0]" ; }

    $KEY  = "INFILE_SIMTRAIN_MASTER:" ;  
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0 ) { $INFILE_SIMTRAIN_MASTER = "$tmp[0]" ;  }

    if ( length($INFILE_SIMTRAIN_MASTER) == 0 ) {
	$MSGERR[0] = "Missing required key: INFILE_SIMTRAIN_MASTER";
	sntools::FATAL_ERROR(@MSGERR);	
    }
    else {
	print "\t INFILE_SIMTRAIN_MASTER:  $INFILE_SIMTRAIN_MASTER \n";
    }
    # ---------------

    $KEY  = "INFILE_SIMDATA_MASTER:" ; 
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0 ) {
	$INFILE_SIMDATA_MASTER = "$tmp[0]" ;
	print "\t INFILE_SIMDATA_MASTER:  $INFILE_SIMDATA_MASTER \n";
    }
    else
    { $DOSTAGE[$ISTAGE_SIMGEN_SIMDATA] = 0; }


    $KEY  = "INFILE_FIT_MASTER:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_ABORT, @CONTENTS_MASTER);
    $INFILE_FIT_MASTER = "$tmp[0]" ;
    print "\t INFILE_FIT_MASTER:  $INFILE_FIT_MASTER \n";

    $KEY  = "WFALSE:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0 ) { 
	$INPUT_WFALSE = $tmp[0] ; 
	if ( index("$WFALSE_VALID",$INPUT_WFALSE) < 0 ) {
	    die "\n ERROR: Invalid WFALSE = $INPUT_WFALSE \n" .
		"\t Valid WFALSE: $WFALSE_VALID \n";
	}
    }
    print "\t WFALSE:  $INPUT_WFALSE \n" ;

    $KEY  = "NEARNBR_SEPMAX_VARDEF:" ;
    @SEPMAX_VARDEF_LIST  = 
	sntools::parse_array($KEY, 49, $OPT_ABORT, @CONTENTS_MASTER);   
    $NTRAIN = scalar(@SEPMAX_VARDEF_LIST);
    for($i=0; $i < $NTRAIN; $i++ ) {
	$SEPMAX_VARDEF_LIST[$i]  =~ s/$q//g ; # remove optional single quotes
	print "\t SEPMAX_VARDEF_LIST[$i] = '$SEPMAX_VARDEF_LIST[$i]' \n";
    }

    # - - - - - - - - - - - 
    $KEY  = "VERSION_DATA:" ;  # legacy key
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp)>0 )  { $VERSION_REALDATA = "$tmp[0]" ; }
    $KEY  = "VERSION_REALDATA:" ;  
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp)>0 )  { $VERSION_REALDATA = "$tmp[0]" ; }

    print "\t VERSION_REALDATA = $VERSION_REALDATA \n";

# - - - - -
# check SCALE_NON1A[NONIA]

    $KEY  = "SIMTRAIN_SCALE_NON1A:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0  ) { $SIMTRAIN_SCALE_NON1A = "$tmp[0]" ; }

    $KEY  = "SIMTRAIN_SCALE_NONIA:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0  ) { $SIMTRAIN_SCALE_NON1A = "$tmp[0]" ; }

# - - - - - - - - - - 
    $KEY  = "SCALE_NGEN_UNIT_VALIDATA:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0  ) { $SCALE_NGEN_UNIT_VALIDATA = "$tmp[0]" ;  }
    $KEY  = "SCALE_NGEN_UNIT_BIASCOR:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0  ) { $SCALE_NGEN_UNIT_BIASCOR = "$tmp[0]" ;  }

    $KEY  = "NGEN_UNIT_VALIDATA:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0  ) { $NGEN_UNIT_VALIDATA = "$tmp[0]" ;  }
    $KEY  = "NGEN_UNIT_BIASCOR:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0  ) { $NGEN_UNIT_BIASCOR = "$tmp[0]" ;  }


    if ( $NGEN_UNIT_VALIDATA > .0001 ) 
    { print "\t NGEN_UNIT(VALIDATA): $NGEN_UNIT_VALIDATA \n"; }
    else
    { print "\t NGEN_UNIT(VALIDATA): " . 
	  "$SCALE_NGEN_UNIT_VALIDATA x NGEN_UNIT(TRAIN) \n"; }

    if ( $NGEN_UNIT_BIASCOR > .0001 ) 
    { print "\t NGEN_UNIT(BIASCOR): $NGEN_UNIT_BIASCOR \n"; }
    else
    { print "\t NGEN_UNIT(BIASCOR): $SCALE_NGEN_UNIT_BIASCOR x NGEN_UNIT(TRAIN) \n"; }

# - - - - - - - - - - - - - -
# option to turn OFF automatic symLinks for duplicate generation  
    $KEY  = "SYMLINKS_FLAG:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_MASTER);
    if ( scalar(@tmp) > 0 ) {
	$SYMLINKS_FLAG = $tmp[0] ;
	print "\t SYMLINKS_FLAG:  $SYMLINKS_FLAG  ". 
	    " (for NNFIT_DATA & NNFIT_SIM)\n";
    }
    # -------------------------------------
    print " Done parsing $INFILE_NEARNBR_MASTER. \n\n";

} # end of parse_inFile_master


# ==========================================
sub checkFiles {

    # make some basic checks on input files before launching
    # jobs; abort on any problem.
    #
    # Jun 14 2016: allow VERSION_REALDATA to be in /SIM area
    # Jul 06 2018: refactor and check PATH_SNDATA_SIM.LIST file

    my (@bla, @KEYLIST, @FILELIST, $NKEY, $ikey, $KEY, $FILE );

    #  ----------- BEGIN --------------
    # make sure input files exist

    unless ( -e $INFILE_SIMTRAIN_MASTER ) {
	$MSGERR[0] = "INFILE_SIMTRAIN_MASTER = $INFILE_SIMTRAIN_MASTER" ;
	$MSGERR[1] = "does not exist.";
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ( length($INFILE_SIMDATA_MASTER) > 0 ) {
	unless ( -e $INFILE_SIMDATA_MASTER ) {
	    $MSGERR[0] = "INFILE_SIMDATA_MASTER = $INFILE_SIMDATA_MASTER" ;
	    $MSGERR[1] = "does not exist.";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

    unless ( -e $INFILE_FIT_MASTER ) {
	$MSGERR[0] = "INFILE_FIT_MASTER = $INFILE_FIT_MASTER" ;
	$MSGERR[1] = "does not exist.";
	sntools::FATAL_ERROR(@MSGERR);
    }

    # -------------------------------------------------------
    # make list of invalid keys/files in fit-namelist file
    
    $NKEY = 0;

    $KEYLIST[$NKEY]  = "NNINP" ;
    $FILELIST[$NKEY] = "$INFILE_FIT_MASTER" ;
    $NKEY++ ;

    $KEYLIST[$NKEY]  = "NEARNBR" ;
    $FILELIST[$NKEY] = "$INFILE_FIT_MASTER" ;
    $NKEY++ ;

    for($ikey=0; $ikey < $NKEY; $ikey++ ) {
	$KEY  = $KEYLIST[$ikey] ;
	$FILE = $FILELIST[$ikey] ; 
	@bla = qx(grep $KEY $FILE);
	if ( scalar(@bla) > 0 ) {
	    $MSGERR[0] = "Must remove '$KEY' from file=$FILE";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

    # -------------------------------------------------------
    # make sure that data version exists, and also allow using
    # simulation as the data version.

    my (@DIRLIST, $PATH_SNDATA_LISTFILE, $DIR, $FOUND );
    @DIRLIST = ( "$SNDATA_ROOT/lcmerge/$VERSION_REALDATA" ) ;
    @DIRLIST = ( @DIRLIST, "$SNDATA_ROOT/SIM/$VERSION_REALDATA" );

    if ( length($PRIVATE_DATA_PATH) > 3 )
    { @DIRLIST = ( @DIRLIST, "$PRIVATE_DATA_PATH/$VERSION_REALDATA" ); }

    # July 2018: check more sim paths
    $PATH_SNDATA_LISTFILE = "$SNDATA_ROOT/SIM/PATH_SNDATA_SIM.LIST";
    if ( -e $PATH_SNDATA_LISTFILE ) {
	my @tmpDirs = `cat $PATH_SNDATA_LISTFILE` ;
	foreach $DIR ( @tmpDirs ) {
	    $DIR  = qx(echo $DIR);  # unpack optional ENV
	    $DIR  =~ s/\s+$// ;     # trim trailing whitespace 
	    @DIRLIST = (@DIRLIST, $DIR);
	}	
    }
    
    $FOUND=0;
    foreach $DIR ( @DIRLIST ) {
	@bla  = qx(ls $DIR  2>/dev/null );
	if ( scalar(@bla) > 0 ) { $FOUND = 1; }
    }

    if ( $FOUND == 0 ) {
	$MSGERR[0] = "VERSION_REALDATA  = '$VERSION_REALDATA'";
	$MSGERR[1] = "does not exist.";
	sntools::FATAL_ERROR(@MSGERR);
    }

    # -------------------------------------------------------
    # make sure that root is defined for the training that
    # requires very large histograms.
    if ( length($ROOT_DIR) < 2 ) {
	my $ss     =  "ISTAGE=$ISTAGE_SIMFIT_TRAIN" ;
	$MSGERR[0] = "ROOT_DIR is not defined, but is required";
	$MSGERR[1] = "for the training stage ($ss)";
	sntools::FATAL_ERROR(@MSGERR);
    }


} # end of checkFiles


# ==========================================
sub NNstagePrefix {
    my $c2     = sprintf("%2.2d", $ISTAGE_GLOBAL);
    my $prefix = "NNSTAGE$c2" ;
    return  $prefix ;
} # end of NNstagePrefix

# ==========================================
sub prep_outDir {

    # create $TOPDIR_OUTPUT.
    # If it already exists, clobber it.
    # If there are no slashes, append pwd to have full path.
    # Ju 11 2019: fix bug avoiding TOPDIR delete.

    if ( length($TOPDIR_OUTPUT) < 2 ) {
	die "\n ERROR: TOPDIR_OUTPUT = '$TOPDIR_OUTPUT' \n";
    }

    my ($jslash);
    $jslash = index($TOPDIR_OUTPUT,'/');
    if ($jslash < 0 ) {
	$TOPDIR_OUTPUT = "$CURRENT_DIR/$TOPDIR_OUTPUT" ;
    }

    # define archive subdirs for input files and stdout files.
    $ARCDIR_INPUTS = "$TOPDIR_OUTPUT/Archive_inputs" ;
    $ARCDIR_STDOUT = "$TOPDIR_OUTPUT/Archive_stdout" ;

    # if skipping first stage, keep current TOPDIR_OUTPUT
    if ( $DOSTAGE[1] == 0 ) { return ; } 
# xxx mark delete  if ( $DOSTAGE[$ISTAGE_SIMGEN_VALIDATA] == 0 ) { return ; }

    if ( -d $TOPDIR_OUTPUT ) { qx(rm -r $TOPDIR_OUTPUT); }

    print " mkdir $TOPDIR_OUTPUT \n" ;
    qx(mkdir $TOPDIR_OUTPUT);
    qx(mkdir $ARCDIR_INPUTS);
    qx(mkdir $ARCDIR_STDOUT);

    qx(cp $INFILE_NEARNBR_MASTER $ARCDIR_INPUTS);
    qx(echo $CURRENT_DIR > $ARCDIR_INPUTS/SUBMIT_DIRNAME.LOG );


} # end of prep_outDir

# ============================
sub run_simgen {
    my($ISTAGE) = @_ ;

    # make copy of SIMGEN_MASTER file depending on stage
    # and run sim_SNmix. Each SIMGEN_MASTER copy has a different
    #
    #  * GENVERSION  -->  NN[1,2]_[whoami]_GENVERSION
    #  * RANSEED
    #  * PREFIX      --> unique SIMLOGS_XXX output
    #
    # Also, store GENVERSION(s) to use later for fitting.
    #

    my $SUFFIX = "$SIMGEN_SUFFIX[$ISTAGE]" ;

    print "\n# ==========================================================\n";
    print "\t  PREPARE & RUN SIMULATIONS for $SUFFIX (ISTAGE=$ISTAGE)\n";

    # ---------- start by parsing SIMGEN_MASTER ------------------

    my $stagePrefix = &NNstagePrefix();
    $TOPDIR_SIMLOGS = "$TOPDIR_OUTPUT/${stagePrefix}_SIMGEN_$SUFFIX" ;

    if ( $DOSTAGE[$ISTAGE] )  { qx(mkdir $TOPDIR_SIMLOGS); }

    # construct commands to modify the SIMGEN_MASTER   
    &make_SIMGEN_MASTER($ISTAGE);    

    # check stage-quit AFTER making SIMGEN_MASTER in order to
    # set things needed for later stages.
    if ( $DOSTAGE[$ISTAGE] == 0 ) { 
	print "\t --> SKIP SIMGEN $SUFFIX \n";
	return ; 
    }


    # submit jobs and wait for DONE stamp.
    &submit_SIMGEN_MASTER($ISTAGE) ;
    
    return ;

} # end of run_simgen


# =============================
sub parse_SIMDATA_MASTER {

    my ($KEY, @tmp);

    unless ( -e $INFILE_SIMDATA_MASTER ) { return; }

    print "\n Parse SIMDATA_MASTER file:  $INFILE_SIMDATA_MASTER \n";

    @CONTENTS_SIMDATA = ();
    sntools::loadArray_fromFile($INFILE_SIMDATA_MASTER,
				\@CONTENTS_SIMDATA);

    
    $KEY  = "GENVERSION:" ;
    @GENV_LIST_SIMDATA = sntools::parse_array($KEY,1, $OPT_ABORT, 
					      @CONTENTS_SIMDATA);
    print "\t GENVERSION_LIST: @GENV_LIST_SIMDATA\n" ;
    $NGENVER_SIMDATA = scalar(@GENV_LIST_SIMDATA);

    $KEY  = "GENPREFIX:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_ABORT, @CONTENTS_SIMDATA);
    $GENPREFIX_ORIG_DATA = "$tmp[0]" ;
    print "\t GENPREFIX:  $GENPREFIX_ORIG_DATA \n";

    # get names of snlc_sim input files in case they are needed
    $INFILE_SIMGEN_Ia_DATA = $INFILE_SIMGEN_NONIa_DATA = "" ;
    
    $KEY = "SIMGEN_INFILE_Ia:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_ABORT, @CONTENTS_SIMDATA);
    $INFILE_SIMGEN_Ia_DATA = $tmp[0] ;

    $KEY = "SIMGEN_INFILE_NONIa:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_SIMDATA);
    if ( scalar(@tmp) > 0 ) { $INFILE_SIMGEN_NONIa_DATA = $tmp[0] ; }

    # get NGEN_UNIT in case it's needed
    $KEY = "NGEN_UNIT:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_ABORT, @CONTENTS_SIMDATA) ;
    $NGEN_UNIT_DATA = $tmp[0] ; 

} # end  parse_SIMDATA_MASTER

# =============================
sub parse_SIMTRAIN_MASTER {

    my ($KEY, @tmp);

    print "\n Parse SIMTRAIN_MASTER file:  $INFILE_SIMTRAIN_MASTER \n";

    @CONTENTS_SIMTRAIN = ();
    sntools::loadArray_fromFile($INFILE_SIMTRAIN_MASTER,
				\@CONTENTS_SIMTRAIN);

    
    $KEY  = "GENVERSION:" ;
    @GENV_LIST_SIMTRAIN = sntools::parse_array($KEY,1, $OPT_ABORT, 
					       @CONTENTS_SIMTRAIN);
    print "\t GENVERSION_LIST: @GENV_LIST_SIMTRAIN \n" ;
    $NGENVER_SIMTRAIN = scalar(@GENV_LIST_SIMTRAIN);

    $KEY  = "GENPREFIX:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_ABORT, @CONTENTS_SIMTRAIN);
    $GENPREFIX_ORIG_TRAIN = "$tmp[0]" ;
    print "\t GENPREFIX:  $GENPREFIX_ORIG_TRAIN \n";

    &parse_SIMGEN_RANSEED() ;

    # get names of snlc_sim input files in case they are needed
    $INFILE_SIMGEN_Ia_TRAIN = $INFILE_SIMGEN_NONIa_TRAIN = "" ;
    
    $KEY = "SIMGEN_INFILE_Ia:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_ABORT, @CONTENTS_SIMTRAIN);
    $INFILE_SIMGEN_Ia_TRAIN = $tmp[0] ;

    $KEY = "SIMGEN_INFILE_NONIa:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_SIMTRAIN);
    if ( scalar(@tmp) > 0 ) { $INFILE_SIMGEN_NONIa_TRAIN = $tmp[0] ; }

    # get NGEN_UNIT in case it's needed
    $KEY = "NGEN_UNIT:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_ABORT, @CONTENTS_SIMTRAIN) ;
    $NGEN_UNIT_TRAIN = $tmp[0] ; 

    # Jun 2019: get current GENOPT_GLOBAL in case we need to append
    $KEY = "GENOPT_GLOBAL:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_SIMTRAIN) ;
    $GENOPT_GLOBAL_TRAIN  = "$tmp[0]" ;
    $GENOPT_GLOBAL_TRAIN  =~ s/\s+$// ;     # trim trailing whitespace 

} # end of parse_SIMTRAIN_MASTER

# =======================================
sub parse_SIMGEN_RANSEED {

    # Created July 11 2017
    # Check both RANSEED and RANSEED_REPEAT

    my ($KEY, @tmp, @wdlist);
    my $FOUND=0;

    $KEY  = "RANSEED:" ;
    @tmp  = sntools::parse_array($KEY,1, $OPT_QUIET, @CONTENTS_SIMTRAIN);
    if ( scalar(@tmp) > 0 ) 
    { $RANSEED[0] = $tmp[0] ;  $FOUND=1;  }

    if ( $FOUND == 0 ) {
	$KEY  = "RANSEED_REPEAT:" ;
	@tmp  = sntools::parse_array($KEY,2, $OPT_QUIET, @CONTENTS_SIMTRAIN);
	if ( scalar(@tmp) > 0 )  { 
	    @wdlist     = split(/\s+/, $tmp[0] ) ;
	    $RANSEED[0] = $wdlist[1] ;  
	    $FOUND=1;  
	}
    }

    if ( $FOUND == 0 ) {
	$MSGERR[0] = "Could not find RANSEED or RANSEED_REPEAT key";
	$MSGERR[1] = "Check simgen-master file" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    $RANSEED[$ISTAGE_SIMGEN_SIMDATA]   = int(2.773*$RANSEED[0]);  
    $RANSEED[$ISTAGE_SIMGEN_VALIDATA]  = int(1.000*$RANSEED[0]);  
    $RANSEED[$ISTAGE_SIMGEN_REF]       = int(2.141*$RANSEED[0]);   
    $RANSEED[$ISTAGE_SIMGEN_TRAIN]     = int(3.072*$RANSEED[0]);   
    $RANSEED[$ISTAGE_SIMGEN_BIASCOR]   = int(4.093*$RANSEED[0]);   

    print "\t RANSEED: @RANSEED \n";

} # end parse_SIMGEN_RANSEED

# =======================================
sub make_SIMGEN_MASTER {
    my ($ISTAGE) = @_ ; 

    # use 'sed' to make several substitutions in the simgen-master file.

    my ($INFILE_MASTER, $CLOSE_INFILE_MASTER);
    my ($VER_PREFIX, $genPrefix, $SUFFIX );
    my ($SEDCMD, @GENV_LIST, $GENV, $GENV_NN, $KEY, $iver, $stagePrefix);
    my ($NEW_NGEN_UNIT, $GENPREFIX_ORIG);

    if( $ISTAGE == $ISTAGE_SIMGEN_SIMDATA ) {
	$INFILE_MASTER  = "$INFILE_SIMDATA_MASTER";
	$GENPREFIX_ORIG = $GENPREFIX_ORIG_DATA ;
	@GENV_LIST      = @GENV_LIST_SIMDATA;
    }
    else {
	$INFILE_MASTER  = "$INFILE_SIMTRAIN_MASTER";
	$GENPREFIX_ORIG = $GENPREFIX_ORIG_TRAIN ;
	@GENV_LIST      = @GENV_LIST_SIMTRAIN;
    }

    # set arrays needed whether or not this stage runs
    $stagePrefix   = &NNstagePrefix();
    $SUFFIX        = $SIMGEN_SUFFIX[$ISTAGE];
    $VER_PREFIX    = "${GLOBAL_PREFIX}_$SUFFIX" ;
    $CLOSE_INFILE_MASTER = "${stagePrefix}_SIMGEN_MASTER_${SUFFIX}.INPUT" ;

    $CLONE_INFILE_SIMGEN_MASTER[$ISTAGE] = "$CLOSE_INFILE_MASTER" ;

    $genPrefix = "${SUFFIX}";
    $SIMGEN_GENPREFIX[$ISTAGE] = "$genPrefix" ;

    $iver = 0 ;
    foreach $GENV (@GENV_LIST) {	
	$GENV_NN                      = "${VER_PREFIX}_${GENV}" ; 
	$GENV_LIST_NN[$ISTAGE][$iver] = "$GENV_NN" ;
	$iver++ ;
    }


    # check to run or skip this stage
    if ( $DOSTAGE[$ISTAGE] == 0 ) { return ; }

    # ---------- start making sed command to modify input file -----------
    $SEDCMD = "sed " ;    
    $SEDCMD = "$SEDCMD -e 's/$RANSEED[0]/$RANSEED[$ISTAGE]/g'" ;

    # check option to change NGEN_UNIT (for VALIDATA and BIASCOR)
    $NEW_NGEN_UNIT = &get_NGEN_UNIT($ISTAGE); 
    if ( $NEW_NGEN_UNIT > 0.001 )  { 
	$SEDCMD = "$SEDCMD -e '/NGEN_UNIT/c NGEN_UNIT: $NEW_NGEN_UNIT '" ;
    }

####    $INFILE_MASTER

    # For biasCor, remove CC and add global options
    if ( $ISTAGE == $ISTAGE_SIMGEN_BIASCOR ) {
	my $sigArg = "$GENOPT_GLOBAL_TRAIN " .
	    "GENSIGMA_SALT2ALPHA 9000 9000  GENSIGMA_SALT2BETA 9000 9000 " ;  
	$SEDCMD = "$SEDCMD -e " . "'/SIMGEN_INFILE_NONIa/d'" ;
	$SEDCMD = "$SEDCMD -e " . "'/GENOPT_GLOBAL/d'" ;
	$SEDCMD = "$SEDCMD -e " . " '\$ i\GENOPT_GLOBAL: $sigArg'" ;
    }

    # Jun 14 2019: for TRAIN, add global option to scale NON1A (.xyz
    if ( $ISTAGE == $ISTAGE_SIMGEN_TRAIN && $SIMTRAIN_SCALE_NON1A != 1.0 ) {
	my $sigArg = "$GENOPT_GLOBAL_TRAIN " .
	    "DNDZ_SCALE_NON1A $SIMTRAIN_SCALE_NON1A " ;
	$SEDCMD = "$SEDCMD -e " . "'/GENOPT_GLOBAL/d'" ;
	$SEDCMD = "$SEDCMD -e " . " '\$ i\GENOPT_GLOBAL: $sigArg'" ;
    }

    $iver = 0 ;
    foreach $GENV (@GENV_LIST) {	
	$GENV_NN = $GENV_LIST_NN[$ISTAGE][$iver];
	$SEDCMD  = "$SEDCMD -e 's/ $GENV/ $GENV_NN/g'" ;
	$iver++ ;
    }

    $KEY = "TOPDIR_SIMLOGS:" ;
    $SEDCMD = "$SEDCMD -e '1i$KEY $TOPDIR_SIMLOGS'" ;
    $SEDCMD = "$SEDCMD -e '/$KEY/d'" ;
    $SEDCMD = "$SEDCMD -e 's/ $GENPREFIX_ORIG/ $genPrefix/g'" ;       
    print " Create $CLOSE_INFILE_MASTER \n";
    qx($SEDCMD $INFILE_MASTER > $CLOSE_INFILE_MASTER);
    
    return ;

} # end of  make_SIMGEN_MASTER 


# ===================================
sub get_NGEN_UNIT {
    my ($ISTAGE) = @_ ;

    my $NEW_NGEN_UNIT = -9.0 ;

    # for data, reduce stats ... later, change model
    if ( $ISTAGE == $ISTAGE_SIMGEN_VALIDATA ) {
	if ( $NGEN_UNIT_VALIDATA > 0.001 ) 
	{ $NEW_NGEN_UNIT = $NGEN_UNIT_VALIDATA; } 
	else 
	{ $NEW_NGEN_UNIT = $NGEN_UNIT_TRAIN * $SCALE_NGEN_UNIT_VALIDATA; }
    }

    if ( $ISTAGE == $ISTAGE_SIMGEN_BIASCOR ) {
	if ( $NGEN_UNIT_BIASCOR > 0.001 ) 
	{ $NEW_NGEN_UNIT = $NGEN_UNIT_BIASCOR; }
	else 
	{ $NEW_NGEN_UNIT = $NGEN_UNIT_TRAIN * $SCALE_NGEN_UNIT_BIASCOR; }
    }

    return($NEW_NGEN_UNIT);

} # end get_NGEN_UNIT

# ===================================
sub submit_SIMGEN_MASTER  {

    my ($ISTAGE) = @_;

    my $GENPREFIX     = $SIMGEN_GENPREFIX[$ISTAGE];
    my $SUFFIX        = $SIMGEN_SUFFIX[$ISTAGE] ;
    my $SIMGEN_MASTER = $CLONE_INFILE_SIMGEN_MASTER[$ISTAGE] ;

    my ($DONELIST, $SIMGEN_STDOUT );

    # ------------- BEGIN ---------------

    print "\n" ;

    $SIMGEN_STDOUT = "$ARCDIR_STDOUT/SIMGEN_${SUFFIX}.STDOUT" ;
    qx(cp $SIMGEN_MASTER $ARCDIR_INPUTS);
    print "   Launch $SCRIPT_SIM  $SIMGEN_MASTER ... \n" ;
    system("$SCRIPT_SIM $SIMGEN_MASTER >& $SIMGEN_STDOUT &" );
    sleep(2);

    # setup the wait

    $DONELIST = "$TOPDIR_SIMLOGS/SIMLOGS_${GENPREFIX}/SIMJOB_ALL.DONE" ;
    $SIMGEN_ALLDONE_FILE = "$TOPDIR_SIMLOGS/SIMGEN_${GENPREFIX}.DONE" ;

    print "\t Wait for $SCRIPT_SIM  jobs to finish ... \n";
    qx($SCRIPT_WAIT  1  $DONELIST  $SIMGEN_ALLDONE_FILE);
    
    print "\t Done with $SCRIPT_SIM  jobs. \n";

    # remove cloned SIMGEN_MASTER files since they have been 
    # copied to TOPDIR_SIMLOGS
    qx(rm $CLONE_INFILE_SIMGEN_MASTER[$ISTAGE]);


    # check for ABORTS in summary file(s)
    my @bla = qx(grep abort $DONELIST);
    if ( scalar(@bla) > 0 ) {
	die "\n FATAL ABORT found in $SCRIPT_SIM jobs; \n\t check $DONELIST\n";
    }

} # end of sub submit_SIMGEN_MASTER


# ================================================
sub parse_SNFIT_MASTER {

    my (@tmp, $KEY, $inFile, $tmpFile, $string);

    $inFile = $INFILE_FIT_MASTER ;
    # read optional FITOPT lines
    $KEY = "FITOPT:" ;
    @tmp =  sntools::parse_line($inFile, 99, $KEY, $OPT_QUIET);
    @FITOPT_STRING = ( "NONE", @tmp );
    $NFITOPT = scalar(@FITOPT_STRING);

    # check if SALT2mu will run, and make sure SALT2mu inFile exists.
    $KEY = "SALT2mu_INFILE:" ;
    @tmp = sntools::parse_line($inFile, 1, $KEY, $OPT_QUIET);
    $DO_SALT2mu = scalar(@tmp);

    # check that SALT2mu inFile exists
    if ( $DO_SALT2mu ) {
	$tmpFile = $tmp[0] ;
	unless ( -e $tmpFile ) {
	    $MSGERR[0] = "SALT2mu_INFILE = $tmpFile" ;
	    $MSGERR[1] = "does not exist.";
	    sntools::FATAL_ERROR(@MSGERR);
	}

	# Jan 2017: check for BBC fit (Beams with BiasCor in SALT2mu fit)
	$KEY = "varname_pIa=" ;
	@tmp = qx(grep $KEY $tmpFile);       
	$DO_BBC = scalar(@tmp) ;
	my $c1 = substr($tmp[0],0,1);
	if ( $c1 eq '#' ) { $DO_BBC=0; }  # check if commented out
	if ( $c1 eq '!' ) { $DO_BBC=0; }
    }

    # read optional PRIVATE_DATA_PATH (Oct 6 2016)
    $KEY = "PRIVATE_DATA_PATH" ;
    $PRIVATE_DATA_PATH = sntools::get_namelistValue($inFile, $KEY);


} # end of parse_SNFIT_MASTER 

# ================================================
sub run_lcfit {

    my ($ITRAIN, $IVER_SIM, $STRINGOPT ) = @_ ;

    # * construct nml file for split_and_fit.
    # * launch it with spit_and_fit.
    # * abort on warnings returned from split_and_fit.
    # * return when all fit jobs are done.
    #
    # ITRAIN is the NN-grid (not SIMTRAIN):
    #   ITRAIN <  0    --> fit reference
    #   ITRAIN >= 0    --> fit sim to train SEPMAX
    #
    # IVER_SIM is the SIMTRAIN sample
    #
    # STRINGOPT = "TRAIN"   --> fit the training on all sims
    # STRINGOPT = "BIASCOR" --> fit sim only for biasCor sample
    # STRINGOPT = "CCPRIOR" --> fit sim only for CCprior sample
    # STRINGOPT = "FINAL"   --> fit REALDATA, SIMDATA, VALIDATA with training
    #
    # IVER_SIM = index for sim-version on FINAL fits only.
    #
    # Jan 20 2016: remove APPEND_FITRES key from REF and TRAIN
    # Apr 07 2016: treat new STRINGOPT = "BIASCOR" and "CCPRIOR"
    # Feb 08 2017: use MLAPPLY mode for split_and_fit
    # Feb 17 2019: include SIMDATA versions
     
    my ($ISIM, $iver, $SEPMAX_VARDEF, @SEPMAX_ARRAY, $PATH );
    my ($SUFFIX, $TXTS2MU, $COMMENT, $IVERMIN, $IVERMAX );
    my ($NMLFILE, $LOGFILE, $DONEFILE, $SEDCMD, $GENV);
    my ($tmpDir, $tmpFile, $OUTDIR, $PREFIX );
    my ($sedRoot, $sedPS );
    my ($line, $fitopt, $NUMOPT, $opt, $cNUM, $textAdd, $NNLOG );
    my ($GENV_REF, $ARG, $SIM_PRESCALE, $VERSION_TMPSIM );
    my ($VERSION_SIMDATA, $VERSION_VALIDATA);

    # ---------------- BEGIN ----------------
    
    my $DOFIT_REF     = ($ITRAIN < 0) ;
    my $DOFIT_TRAIN   = ($ITRAIN >= 0 && $STRINGOPT eq "TRAIN"    );
    my $DOFIT_BIASCOR = ($ITRAIN >= 0 && $STRINGOPT eq "BIASCOR"  );
    my $DOFIT_CCPRIOR = ($ITRAIN >= 0 && $STRINGOPT eq "CCPRIOR"  );
    my $DOFIT_FINAL   = ($ITRAIN >= 0 && $STRINGOPT eq "FINAL"    );

    my $USE_NNINP   = $DOFIT_TRAIN ;  # use SNANA NNINP for training only
    my $DO_MLAPPLY  = ($DOFIT_BIASCOR || $DOFIT_CCPRIOR || $DOFIT_FINAL);
    my $LOCAL_FITARGS = "" ;
    my $stagePrefix = &NNstagePrefix();
    my $TXTS2MU     = "" ;

    
    if ( $DOFIT_REF ) {
        $ISIM       = $ISTAGE_SIMGEN_REF ;
# xxx	$SUFFIX     = "REF_fitSIM" ;
	$SUFFIX     = "SIMFIT_REF" ;
	$PREFIX     = "${stagePrefix}_$SUFFIX" ;
	$OUTDIR     = "$TOPDIR_OUTPUT/${stagePrefix}_$SUFFIX" ;
	$OUTDIR_REF = $OUTDIR;  
	$COMMENT    = "" ;
	$IVERMIN=0;  $IVERMAX = $NGENVER_SIMTRAIN-1 ;
    }
    elsif ( $DOFIT_TRAIN )  {
	$ISIM       = $ISTAGE_SIMGEN_TRAIN ;
#	$SUFFIX     = "TRAINGRID${ITRAIN}_fitSIM" ;
	$SUFFIX     = "SIMFIT_TRAIN_GRID${ITRAIN}" ;
	$PREFIX     = "${stagePrefix}_$SUFFIX" ;
	$OUTDIR     = "$TOPDIR_OUTPUT/${stagePrefix}_$SUFFIX" ;
	$OUTDIR_TRAIN[$ITRAIN] = $OUTDIR; # global
	$SEPMAX_VARDEF = $SEPMAX_VARDEF_LIST[$ITRAIN] ;   
	$COMMENT = "  SEPMAX_VARDEF = $SEPMAX_VARDEF ";
	$IVERMIN=0;  $IVERMAX=$NGENVER_SIMTRAIN-1 ;

	$TRAINFILE_PATH_MAPFILE 
	    = "$TOPDIR_OUTPUT/TRAIN${ITRAIN}_PATH_MAPFILE.TXT";
	$tmpFile = $TRAINFILE_PATH_MAPFILE ;
	if ( -e $tmpFile ) { qx(rm $tmpFile); }
	qx(touch $tmpFile);
    }
    elsif ( $DOFIT_BIASCOR ) {
        $ISIM       = $ISTAGE_SIMGEN_BIASCOR ;
#	$SUFFIX     = "BIASCOR_GRID${ITRAIN}_fitSIM$IVER_SIM" ;
	$SUFFIX     = "SIMFIT_BIASCOR_TRAIN${IVER_SIM}_GRID${ITRAIN}" ;
	$PREFIX     = "${stagePrefix}_$SUFFIX" ;
	$OUTDIR     = "$TOPDIR_OUTPUT/${stagePrefix}_$SUFFIX" ;
	$NNLOG      = $NNLOG_LIST[$ITRAIN][$IVER_SIM][0] ;
	@SEPMAX_ARRAY = &get_SEPMAX_VARDEF($NNLOG,$INPUT_WFALSE) ;
	$SEPMAX_VARDEF = "@SEPMAX_ARRAY";
	$COMMENT    = "\t for SNIa Bias Corrections" ;
	$GENV       = $GENV_LIST_NN[$ISIM][$IVER_SIM] ;
	$BIASCOR_PATH[$ITRAIN][$IVER_SIM] = "$OUTDIR/$GENV" ; # global
	$IVERMIN=$IVER_SIM;  $IVERMAX=$IVER_SIM ; 
	$LOCAL_FITARGS = "TRAILMARK" ; # write OUTDIR in sim-README file
    }
    elsif ( $DOFIT_CCPRIOR ) {
        $ISIM       = $ISTAGE_SIMGEN_REF ;
#	$SUFFIX     = "CCPRIOROR_GRID${ITRAIN}_fitSIM$IVER_SIM" ;
	$SUFFIX     = "SIMFIT_CCPRIOROR_TRAIN${IVER_SIM}_GRID${ITRAIN}" ;
	$PREFIX     = "${stagePrefix}_$SUFFIX" ;
	$OUTDIR     = "$TOPDIR_OUTPUT/${stagePrefix}_$SUFFIX" ;
	$NNLOG      = $NNLOG_LIST[$ITRAIN][$IVER_SIM][0] ;
	@SEPMAX_ARRAY = &get_SEPMAX_VARDEF($NNLOG,$INPUT_WFALSE) ;
	$SEPMAX_VARDEF = "@SEPMAX_ARRAY";
	$COMMENT    = "\t for CC prior" ;
	$GENV       = $GENV_LIST_NN[$ISIM][$IVER_SIM] ;
	$CCPRIOR_PATH[$ITRAIN][$IVER_SIM] = "$OUTDIR/$GENV" ; # global
	$IVERMIN=$IVER_SIM;  $IVERMAX=$IVER_SIM ; 
    }
    elsif ( $DOFIT_FINAL ) {
        $ISIM       = $ISTAGE_SIMGEN_VALIDATA ;
#	$SUFFIX     = "FINALGRID${ITRAIN}_fitDATA+SIM$IVER_SIM" ;
	$SUFFIX     = "DATAFIT_TRAIN${IVER_SIM}_GRID${ITRAIN}";
	$PREFIX     = "${stagePrefix}_$SUFFIX" ;
	$OUTDIR     = "$TOPDIR_OUTPUT/${stagePrefix}_$SUFFIX" ;
	$NNLOG      = $NNLOG_LIST[$ITRAIN][$IVER_SIM][0] ;
	@SEPMAX_ARRAY = &get_SEPMAX_VARDEF($NNLOG,$INPUT_WFALSE) ;
	$SEPMAX_VARDEF = "@SEPMAX_ARRAY";
	$COMMENT    = "" ;
	$TXTS2MU    = "(DO_SALT2mu=$DO_SALT2mu)" ;
	$IVERMIN=$IVER_SIM;  $IVERMAX=$IVER_SIM ; # VALIDATA
    }
    else {
	die "\n FATAL ERROR: UNKNOWN ITRAIN=$ITRAIN\n";
    }

    $NMLFILE  = "${PREFIX}.NML" ;
    $LOGFILE  = "$ARCDIR_STDOUT/${PREFIX}.STDOUT" ;
    $DONEFILE = "$OUTDIR/ALLFITS.DONE" ;

    print "\n ================================================ \n";
    print "  Fit with $NMLFILE  $TXTS2MU \n ";
    print "  run_lcfit args: " .
	"ITRAIN=$ITRAIN, IVER_SIM=$IVER_SIM, " .
	"STRINGOPT=$STRINGOPT (ISIM=$ISIM)\n" ;
    print "$COMMENT \n";

    # ------- check SKIP options -----------
    my $SKIPIT = 0 ;

    if ( $DOFIT_REF    &&  $DOSTAGE[$ISTAGE_SIMFIT_REF]==0 )  
    { $SKIPIT = 1; }

    if ( $DOFIT_TRAIN  &&  $DOSTAGE[$ISTAGE_SIMFIT_TRAIN]==0 ) 
    { $SKIPIT = 1; }

    if ( $DOFIT_BIASCOR  &&  $DOSTAGE[$ISTAGE_SIMFIT_BIASCOR]==0 ) 
    { $SKIPIT = 1; }

    if ( $DOFIT_CCPRIOR  &&  $DOSTAGE[$ISTAGE_SIMFIT_CCPRIOR]==0 ) 
    { $SKIPIT = 1; }

    if ( $DOFIT_FINAL  &&  $DOSTAGE[$ISTAGE_FINALFIT]==0 ) 
    { $SKIPIT = 1; }

    if ( $SKIPIT ) {
	print "\t --> SKIP $SCRIPT_FIT $NMLFILE \n";
	return ;
#	if ( $DOFIT_FINAL ) { return ; }
#	goto  DO_SUMMARY ;
    }


    # -----------------------
    # construct sed command to remove existing OUTDIR and VERSION keys,
    # and write in the same keys with the correct values.
    
    $SEDCMD = "sed" ;

    $SEDCMD = "$SEDCMD " . 
	"-e '1i#### Begin auto-insert from NN-pipeline ---------------'" ;
    $SEDCMD = "$SEDCMD -e " . "'1iOUTDIR:  $OUTDIR'" ;

    # ---------- set 'VERSION:' keys -------------
    
    
    if ( $DOFIT_FINAL ) { 
	$SEDCMD = "$SEDCMD " . "-e '1iVERSION: $VERSION_REALDATA'" ;  
	for ( $iver=0; $iver < $NGENVER_SIMDATA; $iver++ ) {
	    $VERSION_TMPSIM = $GENV_LIST_NN[$ISTAGE_SIMGEN_SIMDATA][$iver];
	    $SEDCMD = "$SEDCMD " . "-e '1iVERSION: $VERSION_TMPSIM'" ; 
	}
    }

    # setup validataion data (VALIDATA)
    for($iver=$IVERMIN; $iver <= $IVERMAX; $iver++ ) {
	
	$GENV   = $GENV_LIST_NN[$ISIM][$iver] ;
	$SEDCMD = "$SEDCMD -e '1iVERSION: $GENV'" ;
	
	# create file that maps training version to path containing
	# reference FITOPT* files.
	if ( $DOFIT_TRAIN  ) {
	    $GENV_REF = $GENV_LIST_NN[$ISTAGE_SIMGEN_REF][$iver] ;
	    $line     = "VERSION:  $GENV  '$OUTDIR_REF/$GENV_REF'";
	    qx(echo "$line" >> $TRAINFILE_PATH_MAPFILE);	    
	}
    }
    $SEDCMD  = "$SEDCMD -e '1i#'" ;


    # For TRAINING, append TRAINING info to each FITOPT line
    if ( $DOFIT_TRAIN  ) {
	for($NUMOPT = 1; $NUMOPT < $NFITOPT; $NUMOPT++ ) {
	    $fitopt = $FITOPT_STRING[$NUMOPT];
	    $cNUM = sprintf("%3.3d", $NUMOPT);
	    $textAdd = "NEARNBR_TRAINFILE_LIST FITOPT${cNUM}.FITRES" ;
	    $SEDCMD  = "$SEDCMD -e '1iFITOPT: $fitopt   $textAdd'" ;
	}
	$SEDCMD  = "$SEDCMD -e '1i#'" ;
	$SEDCMD  = "$SEDCMD -e /FITOPT:/d" ;
    }

    # for training & biasCor, remove HBOOK  and use ROOT to 
    # avoid hbook memory problems. 
    if ( $DOFIT_TRAIN || $DOFIT_BIASCOR || $DOFIT_CCPRIOR ) {
	$sedRoot  = 
	    sntools::sed_nmlInsert($INFILE_FIT_MASTER, 
				   'SNLCINP', 'ROOTFILE_OUT', '"bla.root"');
	$SEDCMD  = "$SEDCMD $sedRoot";
	$SEDCMD  = "$SEDCMD -e /HFILE_OUT/d" ;
    }

    
    # remove some things for REF and TRAIN, but not for FINAL fit
    if ( !$DOFIT_FINAL ) {
	$SEDCMD  = "$SEDCMD -e /SALT2mu_INFILE:/d" ;
    } 
    if ( $DOFIT_REF || $DOFIT_TRAIN )  { 
	$SEDCMD  = "$SEDCMD -e /APPEND_FITRES:/d" ;   # Jan 20 2016
    }
    $SEDCMD = "$SEDCMD -e /SIM_PRESCALE/d" ;  # always remove this key
    
    # ----------------------
    # special options for FINAL fit (on data and SIM)

    if ( $DO_MLAPPLY ) {

	# note that loop skips NUMOPT=0
	for($NUMOPT = 1; $NUMOPT < $NFITOPT; $NUMOPT++ ) {
	    $fitopt = $FITOPT_STRING[$NUMOPT];
	    $cNUM   = sprintf("%3.3d", $NUMOPT);

	    # get trained SEPMAX for this FITOPT > 0
	    $NNLOG = $NNLOG_LIST[$ITRAIN][$IVER_SIM][$NUMOPT] ;
	    @SEPMAX_ARRAY = &get_SEPMAX_VARDEF($NNLOG,$INPUT_WFALSE) ;
	    $SEDCMD  = "$SEDCMD -e ${qq}1iFITOPT: $fitopt ${qq}" ;
	}

	# add MLAPPLY info (Feb 8 2017)
	$tmpDir  = "$OUTDIR_TRAIN[$ITRAIN]" ;
	$GENV    = "$GENV_LIST_NN[$ISTAGE_SIMGEN_TRAIN][$IVER_SIM]" ;
	$tmpFile = "$tmpDir/$GENV/FITOPT\*_NNSEPMAX.TABLE" ;
	$SEDCMD  = "$SEDCMD -e '1i#'" ;
	$SEDCMD  = "$SEDCMD -e ${qq}1iMLAPPLY_INFILE: $tmpFile ${qq}" ;


	if ( $DOFIT_CCPRIOR ) {
	    my $sed_noIa = 
		sntools::sed_nmlInsert($INFILE_FIT_MASTER, 
				       'SNLCINP', 'USESIM_SNIA', 'F' ) ;
	    $SEDCMD  = "$SEDCMD $sed_noIa" ;
	}

	$SEDCMD  = "$SEDCMD -e '1i#'" ;
	$SEDCMD  = "$SEDCMD -e /FITOPT:/d" ;

    }  # end DO_MLAPPLY
    
    # ----------------------
    # misc.
    $SEDCMD = "$SEDCMD -e '1iGZIP_FLAG:  4'" ;
    $SEDCMD = "$SEDCMD -e '1iDONE_STAMP:  $DONEFILE'" ;

    $SEDCMD = "$SEDCMD -e '1i#### End auto-insert ---------------- '" ;

    # remove existing keys from user that are auto-generated by this script.
    $SEDCMD = "$SEDCMD -e /OUTDIR:/d" ;
    $SEDCMD = "$SEDCMD -e /VERSION:/d" ;
    $SEDCMD = "$SEDCMD -e /GZIP_FLAG:/d" ;
    
    # check for SALT2mu
    if ( $DOFIT_FINAL &&  $DO_BBC ) {
	$PATH = $BIASCOR_PATH[$ITRAIN][$IVER_SIM];
	$SEDCMD = "$SEDCMD " . 
	    " -e '/SALT2mu_INFILE:/a\ SALT2mu_BIASCOR_PATH: $PATH' " ;

	$PATH = $CCPRIOR_PATH[$ITRAIN][$IVER_SIM];
	$SEDCMD = "$SEDCMD " . 
	    " -e '/SALT2mu_INFILE:/a\ SALT2mu_CCPRIOR_PATH: $PATH' " ;
    }

    
    # execute sed command to create new nml file.
    qx($SEDCMD $INFILE_FIT_MASTER > $NMLFILE);   

#    if($DOFIT_FINAL) { die "\n\n xxx DIE run_lcfit on FINAL ... \n"; }


    # append &NNINP at end of file (for training only)
    if ( $USE_NNINP ) {
	open  PTRNML , ">> $NMLFILE" ;
	print PTRNML "\n";
	print PTRNML " &NNINP \n";

	if ( $DOFIT_TRAIN ) 
	{ $ARG = $TRAINFILE_PATH_MAPFILE ; }
	else {
	    $GENV_REF = $GENV_LIST_NN[$ISTAGE_SIMGEN_REF][$IVER_SIM] ;
	    $ARG      = "$OUTDIR_REF/$GENV_REF" ;
	}
	print PTRNML "   NEARNBR_TRAINFILE_PATH = \n    '$ARG' \n\n";
	print PTRNML "   NEARNBR_TRAINFILE_LIST = 'FITOPT000.FITRES' \n\n";
	print PTRNML "   NEARNBR_SEPMAX_VARDEF  = \n    '$SEPMAX_VARDEF' \n\n";
	print PTRNML "   NEARNBR_TRUETYPE_VARNAME      = 'SIM_TYPE_INDEX' \n\n";
	print PTRNML "   NEARNBR_TRAINFILE_SCALE_NON1A = $SIMTRAIN_SCALE_NON1A \n\n" ;
	print PTRNML " &END\n";
	close PTRNML ;
    }


    # archive nml-input file just before submitting it
    qx(cp $NMLFILE $ARCDIR_INPUTS);

    # -------------------------------------------------
    # launch it
    my $ALL_FITARGS = "$SCRIPT_FITARGS $LOCAL_FITARGS";
    print "\t Launch $SCRIPT_FIT $NMLFILE  $ALL_FITARGS ... \n" ;
    system("$SCRIPT_FIT $NMLFILE $ALL_FITARGS >& $LOGFILE &" );
    sleep(3);

    # wait for DONE_STAMP 
    print "\t Wait for $SCRIPT_FIT jobs to finish ... \n";
    qx($SCRIPT_WAIT 1  $DONEFILE  $DONEFILE);
    print "\t Done with $SCRIPT_FIT jobs. \n";

    # abort on WARNING in MERGE.LOG [created by split_and_fit]
    my $MERGELOG = "$OUTDIR/MERGE.LOG" ;
    my @bla = qx(grep WARNING $MERGELOG);
    if ( scalar(@bla) > 0 ) {
	die "\n WARNING found in $SCRIPT_FIT job; \n\t check $MERGELOG\n";
    }

    # remove modified nml file since it's already copied into 
    # $OUTDIR by split_and_fit.
  REMOVE_NMLFILE:
    qx(rm $NMLFILE);

    # wait a bit to make sure split_and_fit job is out of queue
    if ( $SKIPIT == 0 ) { sleep(20); }

  
  DO_SUMMARY:

    # for FINAL fit, make generic symbolic links to NNFIT_DATA and NNFIT_SIM
    # to simplify plotting scripts finding directories.
    if ( $DOFIT_FINAL && $SYMLINKS_FLAG ) {
	$VERSION_VALIDATA = $GENV_LIST_NN[$ISTAGE_SIMGEN_VALIDATA][$IVER_SIM];
	unless ( -l "$OUTDIR/$SYMLINK_REALDATA" ) 
	{ qx(cd $OUTDIR ; ln -s $VERSION_REALDATA $SYMLINK_REALDATA ) ; }

	for($iver=0; $iver < $NGENVER_SIMDATA; $iver++ ) {
	    $VERSION_SIMDATA = $GENV_LIST_NN[$ISTAGE_SIMGEN_SIMDATA][$iver];
	    my $itmp=0; my $symLink;
	    unless ( -l "$OUTDIR/$SYMLINK_SIMDATA" ) {
		$symLink = "${SYMLINK_SIMDATA}${itmp}" ;
		qx(cd $OUTDIR ; ln -s $VERSION_SIMDATA $symLink ) ; 	       
	    }
	    $itmp++ ;
	}

	unless ( -l "$OUTDIR/$SYMLINK_VALIDATA" ) 
	{ qx(cd $OUTDIR ; ln -s $VERSION_VALIDATA  $SYMLINK_VALIDATA  ) ; }
    }


    # Jun 10 2016: for BIASCOR, remove large SIMGEN.DAT file 
    if ( $DOFIT_BIASCOR ) {
	$VERSION_TMPSIM = $GENV_LIST_NN[$ISTAGE_SIMGEN_TRAIN][$IVER_SIM]; 
	$tmpFile  = "$OUTDIR/$VERSION_TMPSIM/SIMGEN.DAT" ;
	if ( -e $tmpFile  ) { qx(rm $tmpFile ); }

	$tmpFile  = "$OUTDIR/$VERSION_TMPSIM/SIMGEN.DAT.gz" ;
	if ( -e $tmpFile  ) { qx(rm $tmpFile ); }
    }

    # get average number of SIM events in case we need to pre-scale later
    if ( $DOFIT_REF ) {	
	$NFITSIM_REF = &get_NFIT($OUTDIR_REF);  
	print "\t <NFITSIM_REF> = $NFITSIM_REF \n";
    }

    return ;

} # end of run_lcfit


# =============================================
sub get_NFIT {

    my ($OUTDIR_FIT) = @_ ;

    # examine MERGE.LOG and return average number of 
    # events passing cuts; i.e., after N=

    my ($MERGELOG, @bla, $line, @wdlist, $Nstring, $wd, $jeq);
    my ($NFIT, $NFIT_SUM, $NFIT_AVG, $NMERGE );
    
    $MERGELOG = "$OUTDIR_FIT/MERGE.LOG" ;
    @bla      = qx(grep MERGED $MERGELOG);
    $NFIT_SUM = $NMERGE  = $NFIT_AVG = 0 ;

    foreach $line (@bla) {
	@wdlist   = split(/\s+/,$line ) ;

	# search for 'N='
	$Nstring = "" ;
	foreach $wd (@wdlist) {
	    if ( index($wd,'N=') > 0 ) { $Nstring = $wd ; }
	}
	if ( length($Nstring) == 0 ) {
	    $MSGERR[0] = "Could not find 'N=' in " ;
	    $MSGERR[1] = "line = '$line' " ;
	    $MSGERR[2] = "of $MERGELOG" ;
	    sntools::FATAL_ERROR(@MSGERR);
	}

	# remove ()
	$Nstring  =~ s/\(//g ; 
	$Nstring  =~ s/\)//g ; 
	$jeq = index($Nstring,'=') ;
	
	$NFIT = substr($Nstring,$jeq+1,10);
	$NFIT_SUM += $NFIT ;
	$NMERGE++ ;
    }

    
    if ( $NMERGE > 0 ) {
	$NFIT_AVG = $NFIT_SUM/$NMERGE ;
    }

    return $NFIT_AVG ;

} # end of get_NFIT


# =============================================
sub BBC_summary {

    my ($ITRAIN, $IVER_SIM ) = @_ ;

    # loop over each fit-option and write one-line SALT2mu summary
    my ($iv, $opt);

    my $stagePrefix = sprintf("NNSTAGE%2.2d", $ISTAGE_FINALFIT) ;
#    my $SUFFIX      = "FINALGRID${ITRAIN}_fitDATA+SIM$IVER_SIM" ;
    my $SUFFIX      = "DATAFIT_TRAIN${IVER_SIM}_GRID${ITRAIN}";
    my $PREFIX      = "${stagePrefix}_$SUFFIX" ;
    my $OUTDIR      = "$TOPDIR_OUTPUT/${stagePrefix}_$SUFFIX" ;

    for($opt=0; $opt < $NFITOPT; $opt++ )  {  
	&update_BBC_summary($OUTDIR, $ITRAIN, $IVER_SIM, 
			    $opt, 1, -1); 
	
	for($iv=0; $iv<$NGENVER_SIMDATA; $iv++ ) {
	    &update_BBC_summary($OUTDIR, $ITRAIN, $IVER_SIM, 
				$opt, 1, $iv);	
	}
	&update_BBC_summary($OUTDIR, $ITRAIN, $IVER_SIM, 
			    $opt, 0, -1); 	    
    }


} # end BBC_summary

# =============================================
sub  update_BBC_summary {

    # Inputs:
    # OUTDIR   : argument from split_and_fit job (to find SALT2mu outFiles)
    # ITRAIN   : trainig-grid index
    # IVER_SIM : sim-version index for training
    # IFITOPT  : FITOPT index (refers to FITOPT in fit-namelist)
    # DATAFLAG : 1 for data, 0 for sim
    # ISIMDATA : >=0 for SIMDATA; else REALDATA

    my ($OUTDIR, $ITRAIN, $IVER_SIM, $IFITOPT, 
	$DATAFLAG, $ISIMDATA ) = @_ ;

    my ($FIRST, $LAST, $summaryFile, $stagePrefix);

    if ( $DO_SALT2mu == 0 ) { return ; }

    $FIRST = ( $ITRAIN  == -1 );
    $LAST  = ( $ITRAIN  == -7 );

    # ----------------------------------
    # on first entry, open summary file and make table header


    if ( $FIRST ) {
	$stagePrefix = &NNstagePrefix();
	$summaryFile = "${stagePrefix}_BBC_SUMMARY.README" ;

	print "\n ================================================ \n";
	print "  BBC summary in $summaryFile \n";

	open PTR_SALT2mu , "> $TOPDIR_OUTPUT/$summaryFile" ;
	print "\t Open $summaryFile \n";

	print PTR_SALT2mu
	    " TRAIN-       FIT-\n" ;
	print PTR_SALT2mu
	    " GRID  Sample OPT   NFIT        alpha           beta  " . 
	    "        sigint  scalePCC\n";
	print PTR_SALT2mu "$dashLine \n" ;
	PTR_SALT2mu -> autoflush(1);
	return ;
    }

    # ----------------------------------------------------
    # on last entry, print index legend and close file
    if ( $LAST ) {
	&print_indexLegend(\*PTR_SALT2mu);
	close PTR_SALT2mu ;
	return ;
    }

    # ----------------------------------------------------
    # figure out name of SALT2mu-fitres file that has the
    # goodides (alpha,beta, etc ...)

    my ($VER, $VERSION_SIM, $S2FILE, $F3, @tmp, $KEY, @wdlist );
    my ($NSN, $alpha, $alphaErr, $beta, $betaErr, $sigint, $scalePCC);
    my ($Sample);


    if ( $DATAFLAG && $ISIMDATA<0 )  { 
	$VER    = $VERSION_REALDATA ; 
	$Sample = "REALDATA ";
    }
    elsif ( $DATAFLAG && $ISIMDATA>=0 )  { 
	$VER    = $GENV_LIST_NN[$ISTAGE_SIMGEN_SIMDATA][$ISIMDATA] ; 
	$Sample = "SIMDATA${ISIMDATA} ";	
    }
    else  { 
	$VER    = $GENV_LIST_NN[$ISTAGE_SIMGEN_VALIDATA][$IVER_SIM] ;  
	$Sample = "VALIDATA${IVER_SIM}";
    }
    
    $F3     = sprintf("%3.3d", $IFITOPT);
    $S2FILE = "$OUTDIR/$VER/FITOPT${F3}+SALT2mu.FITRES" ;

    # extract the goodies

    # get number of fit SN
    $KEY = "NSNFIT" ;
    @tmp = sntools::parse_line($S2FILE, 99, $KEY, $OPT_ABORT);
    @wdlist  = split(/\s+/,$tmp[0] ) ;
    $NSN     = $wdlist[1] ;

    $KEY = "alpha0" ;
    @tmp = sntools::parse_line($S2FILE, 99, $KEY, $OPT_ABORT);
    @wdlist  = split(/\s+/,$tmp[0] ) ;
    $alpha   = $wdlist[1];    $alphaErr = $wdlist[3];

    $KEY = "beta0" ;
    @tmp = sntools::parse_line($S2FILE, 99, $KEY, $OPT_ABORT);
    @wdlist  = split(/\s+/,$tmp[0] ) ;
    $beta    = $wdlist[1];    $betaErr = $wdlist[3];

    
    $KEY = "sigint" ;  # Mar 3 2016
    @tmp = sntools::parse_line($S2FILE, 99, $KEY, $OPT_WARN );
    if ( scalar(@tmp) > 0 ) {
	@wdlist  = split(/\s+/,$tmp[0] ) ;    $sigint  = $wdlist[1];
    }
    else {
	$sigint = 0.0 ;
    }


    $scalePCC = 1.0 ;
    if ( $DO_BBC ) {
	$KEY = "scalePCC" ; 
	@tmp = sntools::parse_line($S2FILE, 99, $KEY, $OPT_WARN);
	if ( scalar(@tmp) > 0 ) {
	    @wdlist = split(/\s+/,$tmp[0] ); $scalePCC=$wdlist[1] ;
	}
    }

    # format the numbers
    my $str_NSN      = sprintf("%6d", $NSN);
    my $str_alpha    = sprintf("%6.3f", $alpha);
    my $str_alphaErr = sprintf("%5.3f", $alphaErr);
    my $str_beta     = sprintf("%6.3f", $beta);
    my $str_betaErr  = sprintf("%5.3f", $betaErr);
    my $str_sigint   = sprintf("%5.3f", $sigint);
    my $str_scalePCC = sprintf("%5.2f", $scalePCC);

    print PTR_SALT2mu 
	"  $ITRAIN  $Sample  $IFITOPT  "  .
	"$str_NSN   " . 
	"$str_alpha +- $str_alphaErr  " .
	"$str_beta +- $str_betaErr   " .
	"$str_sigint " .
	"$str_scalePCC " .
	"\n" ;

    if ( $DATAFLAG == 0 ) { print PTR_SALT2mu "\n" ; }

    PTR_SALT2mu -> autoflush(1);

} # end of update_BBC_summary

# =============================================
sub get_SEPMAX_VARDEF {

    my ($NNLOG,$WFALSE) = @_ ;

    # read SEPMAX_VARDEF string from $NNLOG (created by nearnbr_maxFoM.exe)
    # for WFALSE value. Example of string to read:
    #   Wfalse=1 : NEARNBR_SEPMAX_VARDEF = 'zPHOT 0.25 c 0.16 x1 1.1'  
    # give return string of 'zPHOT 0.25 c 0.16 x1 1.1'
    #

    my (@bla, $KEY, @wdlist, $tmp, @SEPMAX_ARRAY );

    $KEY     = "Wfalse=$WFALSE" ;
    @bla     = qx(grep '$KEY' $NNLOG);
    @wdlist  = split(/\s+/,$bla[0] ) ;
    
    $tmp = "@wdlist[5 .. $#wdlist]" ;
    $tmp  =~ s/'//g ;   # remove quotes

    @SEPMAX_ARRAY = split(/\s+/,$tmp ) ;
    return @SEPMAX_ARRAY ;

} # end of get_SEPMAX_VARDEF

# ===================================
sub run_nearnbr_maxFoM {

    # run nearnbr_FoM.exe on each sim-train output for $ITRAIN
    # using the optimized SEPMAX variables.
    # If DOSTAGE is false, process enough to set NNLOG_LIST .
    # Jan 29 2017: add "-outFile $NNOUT" argument.

    my ($ITRAIN) = @_ ;
    my $DO_NNFOM = $DOSTAGE[$ISTAGE_NNFOM] ;

    my ($CMD, $opt, $c3opt, $rootFile, $iver, $GENV, $ARGS, $argW);
    my ($OUTDIR, $NNLOG, $NNOUT, $itrain, $FFF);

    print "\n# ========================================== \n";
    print "  Compute NN-sepmax values for TRAIN$ITRAIN \n";

    # -----------------------------------------------
    # on first itrain, create dummy directory with README file
    # to avoid a missing NNSTAGE## directory that might appear
    # like a failure.
    if ( $ITRAIN == 0  && $DO_NNFOM ) {
	my ($stagePrefix, $tmpDir, $comment, $tmpFile, $NN);
	$stagePrefix = &NNstagePrefix();
        $tmpDir      = "$TOPDIR_OUTPUT/${stagePrefix}_NNFOM" ;
	unless ( -d $tmpDir ) { qx(mkdir $tmpDir); }
	$NN      = "NNSTAGE0$ISTAGE_SIMFIT_TRAIN" ;
	$comment = "See output logs under $NN" ;
	$tmpFile = "$tmpDir/NNFOM_SUMMARY.README" ;
	open  PTR_NNFOM , "> $tmpFile" ;
	print PTR_NNFOM "$comment \n\n";
    }

    # -----------------------------------------------
    # loop over each sim-version and fit-option,
    # and run nearnbr_FoM on the ROOT file created by snlc_sim.exe

    for($iver=0; $iver < $NGENVER_SIMTRAIN; $iver++ ) {

	$GENV    = $GENV_LIST_NN[$ISTAGE_SIMGEN_TRAIN][$iver] ;
	$OUTDIR  = "$OUTDIR_TRAIN[$ITRAIN]/$GENV" ;

	if ( $DO_NNFOM ) 
	{ print "\t --> $GENV \n"; }
	else
	{ print "\t --> SKIP NN-optimazation for $GENV \n"; }

	$argW = "-wfalse_pipeline $INPUT_WFALSE" ;
	for($opt=0; $opt < $NFITOPT; $opt++ ) {
	    $c3opt = sprintf("%3.3d", $opt);
	    $FFF      = "FITOPT${c3opt}" ;
	    $rootFile = "${FFF}.ROOT" ;	   
	    $NNLOG    = "${FFF}_NNSEPMAX.LOG"   ;  # stdout
	    $NNOUT    = "${FFF}_NNSEPMAX.TABLE" ;  # table format
	    $ARGS     = "$rootFile -truetype 1 -outFile $NNOUT  $argW" ;
	    $CMD = "cd $OUTDIR ; $BINARY_NNFOM  $ARGS >& $NNLOG" ;
	    $NNLOG_LIST[$ITRAIN][$iver][$opt] = "$OUTDIR/$NNLOG" ;
	    if ( $DO_NNFOM ) { 
		qx($CMD); 	
		&update_NNFOM_summary($ITRAIN,$iver,$opt);
	    }
	}
    }

    # on last itrain, print index-legend and close summary file.
    if ( $DO_NNFOM && $ITRAIN == $NTRAIN-1 ) 
    {  &update_NNFOM_summary(-9,-9,-9); }

    return ;
    
}  # end of sub run_nearnbr_maxFoM


# ========================
sub update_NNFOM_summary {

    # parse log file and print one-line  update to PTR_NNFOM.

    my ($ITRAIN,$IVER,$IFITOPT) = @_ ;

    my ($NNLOG, @bla, @wdlist, $line);
    my ($Pur_noNN, $EFF,$PUR,$FOM, $WARNINGS  );

    # ----------------------------------
    # make table header on first entry
    if ( $ITRAIN == 0 && $IVER == 0 && $IFITOPT == 0 ) {
	print PTR_NNFOM 
	    " TRAIN-  SIM-  FIT-                                  SEPMAX-bndry\n";
	print PTR_NNFOM 
	    " GRID   TRAIN  OPT  Pur_noNN   EFF  x  PUR = FOM      WARNINGS\n";
	print PTR_NNFOM 
	    " --------------------------------------------------------------- \n";
	PTR_NNFOM -> autoflush(1);
    }

    # ----------------------------------------------------
    # on last entry, print index legend and close file
    if ( $ITRAIN == -9 ) {
	&print_indexLegend(\*PTR_NNFOM);
	close PTR_NNFOM ;
	return ;
    }

    # -------------- STANDARD 1-LINE UPDATE -------------

    # first extract summary info using grep
    $NNLOG = $NNLOG_LIST[$ITRAIN][$IVER][$IFITOPT] ;
    @bla   = qx(grep SUMMARY_for_PIPELINE $NNLOG);
    foreach $line (@bla) {
	@wdlist  = split(/\s+/,$line ) ;
	if ( index($line,'Purity_noNN') > 0 ) {
	    $Pur_noNN = $wdlist[1];  # purity with no NN
	}
	elsif ( index($line,'FOM') > 0 ) {
	    $EFF  = $wdlist[1] ;
	    $PUR  = $wdlist[3] ;
	    $FOM  = $wdlist[5] ;
	}
	elsif ( index($line,'WARNINGS') > 0 ) {
	    $WARNINGS = "@wdlist[1 .. $#wdlist]" ; 
	    $WARNINGS =~ s/$q//g ;        # remove quotes
	    if ( length($WARNINGS) == 0 ) { $WARNINGS = "None" ; }
	}
    }

    print PTR_NNFOM 
	"   $ITRAIN      $IVER     $IFITOPT    " .
	"$Pur_noNN    $EFF x $PUR = $FOM      $WARNINGS\n" ;
    PTR_NNFOM -> autoflush(1);

} # end of sub update_NNFOM_summary 


# ==========================
sub print_indexLegend {
    my ($FH) = @_ ;   # pass file handle as input arg
    
    my ($i, $iver, $opt, $GENV );
    print { $$FH } "\nINDEX LEGEND: \n";

    for($i=0; $i < $NTRAIN; $i++ ) {
	print { $$FH }
	    " TRAINGRID = $i : '$SEPMAX_VARDEF_LIST[$i]' \n";
    }

    print { $$FH } " REALDATA  =     $VERSION_REALDATA \n";

    # print SIMDATA
    for($iver=0; $iver < $NGENVER_SIMDATA; $iver++ ) {
	$GENV    = $GENV_LIST_NN[$ISTAGE_SIMGEN_SIMDATA][$iver] ;
	print { $$FH } " SIMDATA   = $iver : $GENV \n";
    }

    # print SIMTRAIN
    for($iver=0; $iver < $NGENVER_SIMTRAIN; $iver++ ) {
	$GENV    = $GENV_LIST_NN[$ISTAGE_SIMGEN_TRAIN][$iver] ;
	print { $$FH } " SIMTRAIN  = $iver : $GENV \n";
    }

    for($opt=0; $opt < $NFITOPT; $opt++ ) {	    
	print { $$FH } " FITOPT    = $opt : $FITOPT_STRING[$opt] \n";
    }

    print { $$FH } "\n";

    $$FH  -> autoflush(1);
    return ;

} # end of print_indexLegend 
