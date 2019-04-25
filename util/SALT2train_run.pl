#!/usr/bin/perl
#
#__________
# Updated version to work on fulla
#
# Usage:
#
# ./SALT2train_run.pl <input file>
#
# This code was written to run SALT2 training,
# given a training set.
#
#
# The input file has the basic syntax:
#
# CONFIG_FILE: pcafit02_JLMTRAIN_cad10_bothSN30.config
# GENVERSION: JLMTRAIN_UBVRI_t2 JLMTRAIN_SDSS_t2 MODELNAME JLMTRAIN_cad10_bothSN30
#
# Here, the CONFIG_FILE: keyword points to the pcafit configuration file to be used,
# and the GENVERSION: keyword includes all SALT2 trainingsets to be used,
# and a MODELNAME to be used. To run multiple trainings, use multiple sets of these
# two keywords. Each training will result in a "*.runit" file created in the training
# directory, containing all the commands needed to run the training. The code will
# automatically source the runit file.
#
# At the top of the file, additional keywords relating to job management may be included.
# On a system with multiple nodes, such as sdssdp62, multiple trainings can be spawned with
# the keywords below:
#
# NODELIST: des04 des05 des06 des07 des08 des09 des10 sdssdp61 sdssdp60 sdssdp59 sdssdp58
# OUTDIR:   /data/dp62.a/data/analysis/SALT2-training/JLM_TESTS_AUG28t1/STAGE03_TrainRun/workSpace
# DONE_STAMP:  /data/dp62.a/data/analysis/SALT2-training/JLM_TESTS_AUG28t1/RUN03_TrainRun.DONE
#
# Here, the NODELIST: keyword gives the list of nodes available for training.
# The OUTDIR: keyword tells the code where to put the top level training directory. Each
# individual training will have its own subdirectory within OUTDIR: . Finally, a global done
# stamp file will be created at the conclusion of all the jobs. By entering the DONE_STANP:
# keyword, the user can control the name and location of that file.
#
#
# Training defaults are currently written into the script
# with no option for changing. This will be updated in
# future versions.
#
#####
# to-do's
# ______
# 
# rewrite hardcoded() to read defaults from config file
# keywords will have to be reconciled 
#
#
#####
# history
#
# June 9: A DONE_STAMP file may be declared in the user input file, via the keyword DONE_STAMP: . 
# Script touches a file this file when all processes have finished.
# The default DONE_STAMP file is ./ALL.DONE ( where ./ is the directory from which train is run ). 
# This feature should work with both single and multi-train runs.
#
# AUG 31: moved most hard-coded items to sub hardcoded().  
#         printing necessary SALT2training SETENV commands to RUNIT scripts.
#         node and memory information written to head of train.log file. 
#         to avoid overloading nodes, BUSY files will write to SNANA_MODELPATH/.
#         these files will be erased when SALT2train_run.pl train.input CHECKDONE is run. 
#         (i.e. at end of every pipeline RUNIT script)
#         JLM
#
# SEP 08: fixed bug in SALT2.INFO default loading subroutine
#         JLM
#
# SEP 14: rewrote code to generate "runit" scripts for all run modes, 
#         to be able to run in batch mode, if given the keyword input
#         BATCH_INFO: <BATCH_COMMAND> <BATCH_TEMPLATE> <BATCH_NCORE>
#    
#     
# Sep 21, 2011 RK few tweaks to run on fulla, mainly distinguishing
#                 bash and csh
#
# Sep 26, 2011 JLM fixed and tested checkdone bug, reinstated busy files. 
#                  Changed batch memory from 2gb to 4gb. 
# 
# Oct 17, 2011 JLM added generic_salt2_stat as part of post-pcafit processing.
#                  Updated run_ERRORSNAKE: if this option is "T", will run
#                  pcafit a second time with errorsnake. 
#
# Oct 18, 2011 JLM added option to declare output model resolution, via 
#                  pcafit config file keywords @write_lambin and @write_phasebin. 
# 
# Oct 29, 2011 JLM Changed write_pca_uncert argument to have equal lc, spec
#                  ( both set to 1.0 ) . Old values were 2.1 and 0.93. 
#                  Changed run syntax so that RUN_PCAFIT = T or 2 runs training
#                  twice, a second time with errorsnake. 
#                  If only desire to run pcafit once, without error snake, 
#                  set RUN_PCAFIT = 1. 
#
#
# Nov 1, 2011 JLM  Altered script so that @nfit_iteration may be set in the pcafit
#                  configuration file. If set, this key overrides the global RUN_PCAFIT value.
#                  @nfit_iteration 1 runs pcafit once, @nfit_iteration 2 runs it twice.
#
#
# Nov 4, 2011 JLM  Added RUN_PCAFIT = QUICKTEST option. This skips training, but simulates
#                  all other parts of the pipeline.
#
# Nov 16, 2011 JLM Added DRAWCOLORLAW, DRAWTEMPLATES plots to RUN_POSTCOLORLAWFIT routine. 
#                  Will make UV,UBVRI and colorlaw residual plots for each training. 
#
#                  A reference surface and color law are required to make these plots. 
#                  If these are not specified in the input file via the 
#                  TRAINPLOTS_REFSURFACE_FILE: and TRAINPLOTS_REFCOLORLAW_FILE: 
#                  keywords, the script will use the defaults given in the
#                  SALT2training.config file.  
#            
# Dec 5, 2011 JLM  Fixed dumb bugs preventing input file loading of path settings, 
#                   running of salt2stat. 
#
#
# Jan 26, 2012 JLM Upgraded DRAWTEMPLATES usage to allow user to specify phases.
#                  A phase list should be given via TRAINPLOTS_PHASELIST: <list>
#                  keyword. If ENV setting PDFTK_PATH is non-zero, script will
#                  use utility pdftk to join all lc residual files into all_resid.pdf.
#
# Feb 9, 2012 JLM  Upgraded code to allow user to specify salt2_stat training set to be used. 
#                  If SALT2STAT_TESTLIST_FILE: is specified in the input file, 
#                  code will run using that test file. Otherwise, code will run generic_salt2_stat.
#
# Feb 13, 2012 JLM Fixed dumb bug from last upgrade, added salt2-covariance-per-redshift-bin call 
#                  to user-specified salt2-stat job run.        
#          
# Mar 07 2012 JLM  Added user option to specify model parameter SALT2_MAGERRFLOOR
# 
# May 15 2012 JLM  Added user option to specify initial model surfaces
#                  SALT_TEMPLATE0: and SALT_TEMPLATE1: in SALT2train_run.input file
#
# Sep 11 2012 JLM  Altered combine_traininglist syntax to consider LC delimiters
#                  "_" OR "-" for compatibility with hyprid DATA+SIM trains
#
# Sep 29 2012 JLM  Fixed bug in sub run_POSTCOLORLAWFIT($)
#                  Now salt2_color_correction_final.dat is copied to 
#                  salt2_color_correction.dat for proper transfer to model directory.
#
# Sep 30 2012 JLM  Fixed dumb bug in run_POSTCOLORLAWFIT($) (salt2stat job call)
#
# Mar 21 2013 JLM  Upgraded salt2_stat jobcall to use new salt2_stat syntax
#                  Added MAG_OFFSET: ### to SALT2.INFO output. Value of MAG_OFFSET is 
#                  currently set in hardcoded();
#
# May 23 2013 JLM  Changed output file format from perl "FORMAT" to standard printf
#                  for compatibility with UChicago system. 
#
# May 28 2013 JLM  Upgraded SALT2MODEL_DIR cleanup to avoid clobbering sym links
#                  (improves QUICKTEST functionality). 
#
# Jun 01 2013 JLM  Altered post_colorlaw_fit JOBCALL to lengthen wavelength coverage of
#                  cle_final.list / salt2_color_dispersion.dat
#
#####
### For easy file navigation: sections are as follows:
#
### DECLARATIONS
### OUTPUT FILE FORMATS
### MAIN CODE
### subroutines (listed individually)
#
#####
### TO DO LIST
# remove hardcoded filenames from MOVEMODEL
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
sub jobinfo_dump(@);
sub check_fileexist(@);
sub check_goodexit(@);
sub check_badexit(@);
sub promptUser(@);
sub hardcoded(); 

sub get_inFile(@);
sub checkDONE_ARGV(@);		
sub init_RUNJOBS($);
sub check_RUNMETHOD($);
sub check_DONEINFO($);
sub init_PATHS($);
sub check_CONFIGPATHS();
sub check_TRAIN_LIST();
sub genversion_multi(@);
sub setup_ENV($);
sub make_batchhead($); 
sub make_batchsource(@);
sub make_RUN($);	   
sub spawn_RUN();  # RK  is this function obsolete ???
sub submit_RUN($);
sub make_BUSY($);

sub makeSALT2INFO($); 
sub check_grid($);

sub parse_inFile_SALT2train($);
sub parse_configFile_SALT2train($);
sub run_PCAFIT(@);
sub run_ERRORSNAKE(@);
sub prep_PCAFIT_WERRSNAKE($); 	       
sub run_POSTCOLORLAWFIT($);
sub run_MOVEMODEL($);
sub run_QUICKTEST($);		  
sub MOVEMODEL($);
sub check_DONE();
sub if_DONE();

### DECLARATIONS ###

my $USER_inFile;
my $DEFAULT_inFile;
my $CHECKDONE_FLAG;
my $local_dir; 
my $SALT_MAIN;
my $RUN_PCAFIT;
my $RUN_TEST;
my $RUN_ERRORSNAKE;
my $RUN_POSTCOLORLAWFIT;
my $RUN_MOVEMODEL;
###
my $SNDATA_ROOT;
my $PDFTK_PATH;
###
my (@NODELIST, $N_NODES);
my ($BATCH_COMMAND, $BATCH_TEMPLATE, $BATCH_NCORE);   
my ($BATCH_MAXCORE, $BATCH_MEM);
my $OUTDIR;
###
my $N_TRAINS;
my %TRAIN_PATH;
my %TRAIN_CONFIG;
my %TRAIN_GENVERSION;
my %TRAIN_NGENVERSION;
my %TRAIN_TRAINLIST;
my %TRAIN_NODELIST;
my $CONFIG_FILE; 
my ($SALT2TRAININGPATH, $SALT2_DIR, $SALTPATH, $SNANA_MODELPATH, $SHELL_NAME);
my ($TRAINPLOTS_REFSURFACE_FILE, $TRAINPLOTS_REFCOLORLAW_FILE, @TRAINPLOTS_PHASELIST);
my $SALT2STAT_TESTLIST_FILE;
my $ROOTSYS;
my $MODELNAME;
my $GENVERSION;
my $N_GENVERSION;
my $LCSIM_INPUT;
my $TRAIN_LIST;
###
my ($WRITE_PHASEBIN, $WRITE_LAMBIN, $OPT_PCAFIT);
my $SALTTRAIN_DIR;
my $SALT2MODEL_DIR;
my $SALT_TEMPLATE0;
my $SALT_TEMPLATE1;
my $LCRESID_JOBCALL;
my $COMPERRSN_JOBCALL;
my $POSTCOLORFIT_JOBCALL;
my $WRITEPCATEMP_JOBCALL;
my $WRITEPCAUNC_JOBCALL;
my $SALT2STAT_JOBCALL;
my $DRAWTEMPLATE_JOBCALL;
my $DRAWCOLORLAW_JOBCALL;
my $SALT2_ERRMAPINTERPOPT;
my $SALT2_ERRMAPKCOROPT;
my $SALT2_MAGERRFLOOR;
my $SALT2_SEDFLUXINTERPOPT;
my $GRIDLIM_SEDFLUXINTERPOPT; 
my $COLORLAW_VERSION;
my ($COLOR_OFFSET, $SALT2MAG_OFFSET); # set in hardcoded();
###
my $NDONE_CURRENT; ## use in conj with $N_TRAINS
my $DONE_STAMP;
my $SPAWN_TRAINS;
my $DONEFLAG;
my $BUSYFILE; 
my $DONE_LOCAL = "SPAWN.DONE";
my $DONE_GLOBAL = "ALL.DONE"; 
###
format SPAWNINPUT =
RUN_PCAFIT: @*
$RUN_PCAFIT
RUN_ERRORSNAKE: @*
$RUN_ERRORSNAKE
RUN_POSTCOLORLAWFIT: @*
$RUN_POSTCOLORLAWFIT
RUN_MOVEMODEL: @*
$RUN_MOVEMODEL

CONFIG_FILE: @*
$CONFIG_FILE
SALT2TRAININGPATH: @* MODELNAME @*
$SALT2TRAININGPATH $MODELNAME
TRAIN_LIST: @*
$TRAIN_LIST

.


### Main ###


$local_dir = getcwd; 

$SNDATA_ROOT = $ENV{'SNDATA_ROOT'};

$PDFTK_PATH = $ENV{'PDFTK_PATH'};

$USER_inFile = get_inFile(@ARGV);

hardcoded();

$CHECKDONE_FLAG = checkDONE_ARGV(@ARGV);

if ($CHECKDONE_FLAG == 1) {

    make_banner("CHECKING FOR FINISHED JOBS");

    $DONEFLAG = check_DONE();

    if ( $DONEFLAG == 1 ) { if_DONE(); }

} elsif ($CHECKDONE_FLAG == 2){
    
    make_banner("MOVING MODEL");
    
    MOVEMODEL($USER_inFile);

} else {

    init_RUNJOBS($USER_inFile); #which jobs to run
 
    check_RUNMETHOD($USER_inFile); #how to run jobs

    check_DONEINFO($USER_inFile); #where to put global DONE

    parse_inFile_SALT2train($USER_inFile); #all training info read in
    
    run_TRAINING();
    
}



print "FINISHED SUCCESSFULLY \n ";

### SUBROUTINES

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

sub check_DONEINFO($)
{
    my ($inFile) = shift;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my $NARGALL = 99;
    my @tmp; ###
    my $key; ### consistent with rest of snana sntools.pl usage

    my $key = "DONE_STAMP:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $DONE_STAMP = @tmp[0];
    unless ( $DONE_STAMP eq '' ) { 
	print "\t DONE_STAMP \t $DONE_STAMP \n"; }


    if ( $DONE_STAMP eq '') {
        $DONE_STAMP = $local_dir . "/" . $DONE_GLOBAL;
	print "\t USING DEFAULT DONE_STAMP FILE LOCATION $DONE_STAMP\n";
    }

    
}


sub checkDONE_ARGV(@)
{
    my $NARG = scalar(@_);
    my $tempString;
    my $return_flag;

    ## check that user gave CHECKDONE flag
    $tempString = $ARGV[1];
    
    if ( $tempString eq "CHECKDONE") {
	$return_flag = 1;
	
	$SPAWN_TRAINS = $ARGV[2];
	$DONE_STAMP = $ARGV[3];

	if ( $SPAWN_TRAINS eq '' ){
	    $SPAWN_TRAINS = 1;
	    print "\t USING DEFAULT SPAWN_TRAINS FILE LOCATION $SPAWN_TRAINS\n";
	}
	if ( $DONE_STAMP eq '') {
	    $DONE_STAMP = $local_dir . "/" . $DONE_GLOBAL;
	    print "\t USING DEFAULT DONE_STAMP FILE LOCATION $DONE_STAMP\n";
	}
    } elsif ($tempString eq "MOVEMODEL") {
	$return_flag = 2;

	$SALT2MODEL_DIR = $ARGV[2];

    } else {
	$return_flag = 0;
    }
    return($return_flag);
}

sub hardcoded()
{

    make_banner("LOADING HARDCODED DEFAULT SETTINGS");
    
    $BATCH_MAXCORE =    0  ;  # RK
    $BATCH_MEM     = "4gb" ;  # RK
    
    $DEFAULT_inFile = "$SNDATA_ROOT/models/SALT2/SALT2training.config";
    print "\t\t DEFAULT_inFile: $DEFAULT_inFile \n";
    $SALT_TEMPLATE0 = "salt2-2-0/salt2_template_0.dat";
    print "\t\t SALT_TEMPLATE0: $SALT_TEMPLATE0 \n";
    $SALT_TEMPLATE1 = "salt2-2-0/salt2_template_1.dat";
    print "\t\t SALT_TEMPLATE1: $SALT_TEMPLATE1 \n";

    $LCRESID_JOBCALL = " pca_1_opt1_final.list full_weight_1.fits covmat_1_with_constraints.fits ";    
    $COMPERRSN_JOBCALL = " pca_1_opt1_final.list full_weight_1.fits covmat_1_with_constraints.fits ";    
    $POSTCOLORFIT_JOBCALL =  " lcresiduals.list salt2_color_correction.dat -n 3 -l 2000 10000 ";

    $DRAWTEMPLATE_JOBCALL = " salt2_template_0.dat -c salt2_color_correction.dat -x salt2_template_1.dat "; 
    $DRAWCOLORLAW_JOBCALL = " salt2_color_correction.dat -b -o colorlaw.pdf ";

    $WRITEPCATEMP_JOBCALL = "pca_1_opt1_final.list";
    $WRITEPCAUNC_JOBCALL = "pca_1_opt1_final.list full_weight_1.fits 2 1.0 1.0 ";

    $SALT2STAT_JOBCALL = " -s MEGACAM 0.2 1.0 ";

    $GRIDLIM_SEDFLUXINTERPOPT = 10; 
    $SALT2_SEDFLUXINTERPOPT = 2.;
    $SALT2_ERRMAPINTERPOPT = 1.;
    $SALT2_ERRMAPKCOROPT = 1.;
    $SALT2_MAGERRFLOOR = 0.005;

    $COLORLAW_VERSION = 1;
    $COLOR_OFFSET = 0.0;
    $SALT2MAG_OFFSET = 0.27;

}

sub setup_ENV($)
{
    my ($inFile) = shift;
    my $sourceme_file;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my $NARGALL = 99;
    my @tmp; ###
    my @tmp2;
    my $key; ### consistent with rest of snana sntools.pl usage

    make_banner("SETTING ENV VARIABLES") ;    

    $sourceme_file = $DEFAULT_inFile;
    
    ###LOOK FOR SALT2_DIR, ROOTSYS, SALTPATH in input file
    ###if not found, load defaults from DEFAULT_inFile
    $key = "SALT2_DIR:" ; ## OPTIONAL input file key, REQUIRED TO RUN
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    if ( @tmp[0] eq ''){
	@tmp = sntools::parse_line($sourceme_file, $NARG1, $key, $OPTABORT) ;
	$SALT2_DIR = @tmp[0];
    } else {
	$SALT2_DIR = @tmp[0];
    }

    $key = "SALTPATH:" ; ## OPTIONAL input file key, REQUIRED TO RUN
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    if ( @tmp[0] eq ''){
	@tmp = sntools::parse_line($sourceme_file, $NARG1, $key, $OPTABORT) ;
	$SALTPATH = @tmp[0];
    } else {
	$SALTPATH = @tmp[0];
    }

    $key = "ROOTSYS:"; ## OPTIONAL input file key, REQUIRED TO RUN
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    if ( @tmp[0] eq ''){
	@tmp = sntools::parse_line($sourceme_file, $NARG1, $key, $OPTABORT) ;
	$ROOTSYS = @tmp[0];
    } else {
	$ROOTSYS = @tmp[0];
    }

    $key = "TRAINPLOTS_REFSURFACE_FILE:"; ## OPTIONAL input file key, REQUIRED TO RUN
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    
    if ( @tmp[0] eq ''){
	@tmp = sntools::parse_line($sourceme_file, $NARG1, $key, $OPTABORT) ;
	$TRAINPLOTS_REFSURFACE_FILE = @tmp[0];
    } else {
	$TRAINPLOTS_REFSURFACE_FILE = @tmp[0];
    }

    $key = "TRAINPLOTS_PHASELIST:"; ## OPTIONAL input file key
    @tmp = sntools::parse_line($inFile, $NARGALL, $key, $OPTQUIET) ;
    @tmp2 = split(' ', @tmp[0]);

    if ( @tmp2[0] eq ''){
	@tmp = sntools::parse_line($sourceme_file, $NARG1, $key, $OPTABORT) ;
        @tmp2 = split(' ', @tmp[0]);	
	@TRAINPLOTS_PHASELIST = @tmp2;
    } else {
	@TRAINPLOTS_PHASELIST = @tmp2;
    }
    
    $key = "TRAINPLOTS_REFCOLORLAW_FILE:"; ## OPTIONAL input file key, REQUIRED TO RUN
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    if ( @tmp[0] eq ''){
	@tmp = sntools::parse_line($sourceme_file, $NARG1, $key, $OPTABORT) ;
	$TRAINPLOTS_REFCOLORLAW_FILE = @tmp[0];
    } else {
	$TRAINPLOTS_REFCOLORLAW_FILE = @tmp[0];
    }

    $key = "SALT2STAT_TESTLIST_FILE:"; ## OPTIONAL input file key
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    if ( @tmp[0] eq ''){
    } else {
	$SALT2STAT_TESTLIST_FILE = @tmp[0];
	check_fileexist("SALT2STAT_TESTLIST_FILE", $SALT2STAT_TESTLIST_FILE);
    }

    $key = "SALT2_MAGERRFLOOR:"; ## OPTIONAL input file key
        @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    if ( @tmp[0] eq ''){
    } else {
	$SALT2_MAGERRFLOOR = @tmp[0];
    }
    
    ### MUST set SALT2_DIR, SALTPATH, and SALT2TRAININGPATH to run any job
    ### THESE ARE REQUIRED BY PCAFIT
    $ENV{'SALT2_DIR'} = $SALT2_DIR; 
    $ENV{'SALTPATH'} = $SALTPATH; 
    $ENV{'SALT2TRAININGPATH'} = $SALT2TRAININGPATH;

    ### NEED SALT2MODEL_DIR - USE DEFAULT IF NOT PASSED
    $key = 'SALT2MODEL_DIR:' ; ##OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    if ( @tmp[0] eq '' ){
	$key = 'SNANA_MODELPATH:'; ##REQUIRED
	@tmp = sntools::parse_line($sourceme_file, $NARG1, $key, $OPTABORT) ;
	$SNANA_MODELPATH = @tmp[0];
	$SALT2MODEL_DIR = $SNANA_MODELPATH . '/SALT2.' . $MODELNAME;
	print "\t SNANA_MODELPATH is $SNANA_MODELPATH \n"; 
	print "\t sourcemefile is $sourceme_file \n"; 
	print "\t DEFAULT_inFile is $DEFAULT_inFile \n"; 
    } 

    print "\t SALT2MODEL_DIR:\n";
    print "\t\t $SALT2MODEL_DIR \n"; 
    check_fileexist("SALTTRAIN_DIR", $SALTTRAIN_DIR);
    print "\t SALTTRAIN_DIR:\n";
    print "\t\t $SALTTRAIN_DIR \n"; 
    print "\t SALT2_MAGERRFLOOR:\n";
    print "\t\t $SALT2_MAGERRFLOOR \n";
    
    ### IF DIR *NOT* SYMBOLIC LINK THEN CLEAN UP SALT2MODEL_DIR
    if ( ! -l $SALT2MODEL_DIR && -d $SALT2MODEL_DIR ) {
	print "\tCleaning up SALT2MODEL DIR \n";
	unlink(<$SALT2MODEL_DIR/*>);
    } 


    $SHELL_NAME    = sntools::shellName();  # RK
    
}

sub make_BUSY($)
{
    my $runitfile = shift; 

    my $string;
    my $string2;
    my $nodename;
    my @output;

    ### CREATE BUSY FILE
    ###$nodename = `hostname | awk  'BEGIN {FS="fnal"}; {print $1}' `;
    $nodename = `hostname`;
    $BUSYFILE = $SNANA_MODELPATH . "/SALT2." . $MODELNAME . ".BUSY_" . $nodename;
    $string = "BUSYFILE: " . $BUSYFILE;
    $string2 = chomp($string);
    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "touch $BUSYFILE \n";
    print RUNIT "echo '$string'  >BUSY \n\n"; 
    close(RUNIT);

}

sub run_POSTCOLORLAWFIT($)
{

    my ($runitfile) = shift;

    my @output;
    my $log; 
    my $file;
    my $job;
    my $jobcall1;
    my $jobcall2;
    my $jobcall3;
    my $jobcall4;
    my $jobcall;
    my $ENVDEF;
    my $n_phases;
    my $i_phases;
    my $i_dtbands;
    my $n_dtbands = 6;
    my @DRAWTEMP_PREFIX = ("UV", "U", "B", "V", "R", "I");
    my @DRAWTEMP_WAVE = (2500, 3500, 4500, 5500, 6500, 7500);

    make_banner("RUNNING POSTCOLORLAWFIT");

    print "Run write_pca_templates \n\n";
       $job = $SALT2_DIR . "/bin/write_pca_templates";
       $log = "write_pca_template.log"; 
       $jobcall1 = $WRITEPCATEMP_JOBCALL;
    if ($WRITE_PHASEBIN eq '' ){
	$jobcall2 = " "; 
    } else {
	$jobcall2 = " -p $WRITE_PHASEBIN ";
    }
    if ($WRITE_LAMBIN eq '' ){
	$jobcall3 = " "; 
    } else {
	$jobcall3 = " -l $WRITE_LAMBIN ";
    }

    $jobcall = $jobcall1 . $jobcall2 . $jobcall3 ;

    jobinfo_dump($job, $jobcall,  $log);

    ### @output = `$job $jobcall1 >& $log`;
    
    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "$job $jobcall >& $log \n\n";
    close(RUNIT); 


    print "Run lightcurve_residuals \n\n";
       $job = $SALT2_DIR . "/bin/lightcurve_residuals";
       $log = "light_curve_residuals.log";
       $jobcall1 = "$TRAIN_LIST $CONFIG_FILE";
       $jobcall2 = $LCRESID_JOBCALL; #in hardcoded
       $jobcall = $jobcall1 . $jobcall2 ;
       jobinfo_dump($job, $jobcall, $log);

    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "$job $jobcall >& $log \n\n";
    close(RUNIT); 

    print "Run post_color_law_fit \n\n";
       $job = $SALT2_DIR . "/bin/post_color_law_fit";
       $log = "post_color_law_fit.log";
       $jobcall1 = $POSTCOLORFIT_JOBCALL; #in hardcoded
       jobinfo_dump($job, $jobcall1, $log);

    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "$job $jobcall1 >& $log \n\n";
    print RUNIT "cp salt2_color_correction.dat salt2_color_correction_first.dat \n";
    print RUNIT "cp salt2_color_correction_final.dat salt2_color_correction.dat \n";
    close(RUNIT); 

    print "Run residual plots \n\n";
    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    $n_phases = @TRAINPLOTS_PHASELIST;

    $ENVDEF = sntools::setEnvString($SHELL_NAME, "ROOTSYS", $ROOTSYS) ;
    print RUNIT "$ENVDEF \n";  

    $job = $SALT2_DIR . "/bin/drawtemplates";
    $jobcall1 = $TRAINPLOTS_REFSURFACE_FILE;
    $jobcall2 = $DRAWTEMPLATE_JOBCALL; 


    $i_dtbands = 0;
    while ($i_dtbands < $n_dtbands){
        $i_phases = 0;
	while ( $i_phases < $n_phases ){
	    $jobcall3 = "-p $TRAINPLOTS_PHASELIST[$i_phases] -w " .  $DRAWTEMP_WAVE[$i_dtbands] . " 1000 -b -o  ";
	    $jobcall4 = $DRAWTEMP_PREFIX[$i_dtbands] . "_resid_" . $i_phases . ".dat";
	    $jobcall = $jobcall1 . $jobcall2 . $jobcall3 . $jobcall4 ; 
	    print RUNIT "$job $jobcall \n";
	    $i_phases = $i_phases + 1;
	}
	$i_dtbands = $i_dtbands + 1;
    }
    
    if ($PDFTK_PATH ne '' ){
	$job = $PDFTK_PATH . "pdftk";
	$jobcall = " *.pdf cat output all_lcresid.pdf ";
	print RUNIT "$job $jobcall \n";
    }

    $job = $SALT2_DIR . "/bin/drawcolorlaw"; 
    $jobcall1 = $TRAINPLOTS_REFCOLORLAW_FILE; 
    $jobcall2 = $DRAWCOLORLAW_JOBCALL; 
    $jobcall = $jobcall1 . $jobcall2 ;
    print RUNIT "$job $jobcall \n\n";
    
    close(RUNIT); 



    if ($SALT2STAT_TESTLIST_FILE ne ''){
	print "Run salt2_stat with user test file\n\n";
	$job = $SALT2_DIR . "/bin/salt2_stat";
        $jobcall1 = " $SALT2STAT_TESTLIST_FILE $CONFIG_FILE pca_1_opt1_final.list -t $TRAIN_LIST";
	$jobcall3 = " ";
    } else {
	print "Run generic salt2stat \n\n";
       $job = $SALT2_DIR . "/bin/generic_salt2_stat";
       $jobcall1 = " $CONFIG_FILE pca_1_opt1_final.list cle_final.list";
       $jobcall3 = $SALT2STAT_JOBCALL; #in hardcoded
    }
       $jobcall2 = " -m regularization_masks.list -d cle_final.list -a 0.1 -b 3.2 ";
       $log = "salt2stat.log";
       $jobcall = $jobcall1 . $jobcall2 . $jobcall3; 
       jobinfo_dump($job, $jobcall, $log);

    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "$job $jobcall >& $log \n\n";

    if ($SALT2STAT_TESTLIST_FILE ne ''){
	$job = $SALT2_DIR . "/bin/salt2-covariance-per-redshift-bin";
	$log = "s2stat_binned.log";
	$jobcall = "";
	jobinfo_dump($job, $jobcall, $log);
	print RUNIT "$job $jobcall >& $log \n\n";
    }
	    
    close(RUNIT); 


}



sub run_PCAFIT(@)
{
    my ($runitfile, $runnumber) = @_;

    my @output; 
    my $log;
    my $job;
    my $jobcall1;
    my $jobcall2;
    my $jobcall3;
    my $jobcall;

    make_banner("RUNNING PCAFIT"); 

       $job = $SALT2_DIR . "/bin/pcafit";
       $jobcall1 = "-c $CONFIG_FILE -l $TRAIN_LIST ";


    if ( $runnumber == 0 ) {
	print "\t\t First time running pcafit. \n";
	$log = "pcafit.log"; 
	$jobcall2 = "-s 0 $SALTPATH/$SALT_TEMPLATE0 ";
	$jobcall3 = "-s 1 $SALTPATH/$SALT_TEMPLATE1 "; 
    } elsif ( $runnumber == 1 ) {
	print "\t\t Second time running pcafit - use error snake.\n";
	$log = "pcafit_with_error_snake.log";
	$jobcall2 = "-p pca_1_opt1_final_first.list -d ";
	$jobcall3 = " "; 
    }

    $jobcall = $jobcall1 . $jobcall2 . $jobcall3; 
    jobinfo_dump($job, $jobcall, $log);

    ### my ($runitfile) = shift;
    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "$job $jobcall >& $log \n\n";
    close(RUNIT); 

}

sub prep_PCAFIT_WERRSNAKE($)
{

    my ($runitfile) = shift; 

    my $new_config_file = "training_with_error_snake.conf";
    my $addition_to_config1;
    my $addition_to_config2;

    $addition_to_config1 = "\@ERRORSNAKE_COV_MAT model_covmat_for_error_snake_first.fits";
    $addition_to_config2 = "\@ERRORSNAKE_COV_SCALING salt2_lc_dispersion_scaling_first.dat";
    
    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "cp pca_1_opt1_final.list pca_1_opt1_final_first.list \n";
    print RUNIT "cp model_covmat_for_error_snake.fits model_covmat_for_error_snake_first.fits \n";
    print RUNIT "cp salt2_lc_dispersion_scaling.dat salt2_lc_dispersion_scaling_first.dat \n\n";
    print RUNIT "cat $CONFIG_FILE > $new_config_file \n";
    print RUNIT "echo \"$addition_to_config1\" >> $new_config_file \n";
    print RUNIT "echo \"$addition_to_config2\" >> $new_config_file \n";
    close(RUNIT); 

    $CONFIG_FILE = $new_config_file; 

}

sub run_ERRORSNAKE(@)
{
    my ($runitfile, $runnumber) = shift;
    
    my @output; 
    my $log_pcaunc;
    my $log_errsnake;
    my $log; 
    my $job;
    my $jobcall1;
    my $jobcall2;
    my $jobcall3;
    my $jobcall;

    make_banner("RUNNING ERRORSNAKE"); 


    if ( $runnumber == 0 ) {
	print "\t\t First time through. \n";
	$log_pcaunc = "write_pca_uncert.log";
	$log_errsnake = "compute_error_snake.log";
    } elsif ($runnumber == 1){
	print "\t\t Second time through. \n";
	$log_pcaunc = "pca_uncert_with_error_snake.log";
	$log_errsnake = "compute_error_snake_2.log";
    }

    print "Run write_pca_uncertainties \n\n";
       $job = $SALT2_DIR . "/bin/write_pca_uncertainties";
       $log = $log_pcaunc;
       $jobcall1 = $WRITEPCAUNC_JOBCALL;
       jobinfo_dump($job, $jobcall1, $log);

    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "$job $jobcall1 >& $log \n\n";
    close(RUNIT); 

    print "Run compute_error_snake \n\n";
       $job = $SALT2_DIR . "/bin/compute_error_snake";
       $log = $log_errsnake;
       $jobcall1 = "$TRAIN_LIST $CONFIG_FILE"; 
       $jobcall2 = $COMPERRSN_JOBCALL;
       $jobcall = $jobcall1 . $jobcall2 ; 
       jobinfo_dump($job, $jobcall, $log);

    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT "$job $jobcall >& $log \n";
    close(RUNIT); 

}

sub init_RUNJOBS($)
{
    my ($inFile) = shift;
    my $diestring;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my $NARGALL = 99;
    my @tmp; ###
    my $key; ### consistent with rest of snana sntools.pl usage

    ### CHECK INPUT FILE FOR USER RUNJOBS
    make_banner("CHECKING FOR USER-SELECTED RUNJOBS") ;
    print "USER RUNJOBS inputs are: \n";

    my $key = "RUN_PCAFIT:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $RUN_PCAFIT = @tmp[0];
    unless ( $RUN_PCAFIT eq '' ) { 
	print "\t RUN_PCAFIT \t $RUN_PCAFIT \n"; }
    if ( $RUN_PCAFIT eq 'T' ) {
	$RUN_PCAFIT = 2;
    }

    my $key = "RUN_ERRORSNAKE:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $RUN_ERRORSNAKE = @tmp[0];
    unless ( $RUN_ERRORSNAKE eq '' ) { 
	print "\t RUN_ERRORSNAKE \t $RUN_ERRORSNAKE \n"; }

    my $key = "RUN_POSTCOLORLAWFIT:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $RUN_POSTCOLORLAWFIT = @tmp[0];
    unless ( $RUN_POSTCOLORLAWFIT eq '' ) { 
	print "\t RUN_POSTCOLORLAWFIT \t $RUN_POSTCOLORLAWFIT \n"; }

    my $key = "RUN_MOVEMODEL:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $RUN_MOVEMODEL = @tmp[0];
    unless ( $RUN_MOVEMODEL eq '' ) { 
	print "\t RUN_MOVEMODEL \t $RUN_MOVEMODEL \n"; }


    ### LOAD DEFAULT RUNJOBS

    if ( $RUN_PCAFIT eq '' ) {
        $RUN_PCAFIT = "1";
        print "\t USING DEFAULT RUNJOB RUN_PCAFIT $RUN_PCAFIT\n";
        }
    if ( $RUN_ERRORSNAKE eq '' ) {
        $RUN_ERRORSNAKE = "T";
        print "\t USING DEFAULT RUNJOB RUN_ERRORSNAKE $RUN_ERRORSNAKE\n";
        }
    if ( $RUN_POSTCOLORLAWFIT eq '' ) {
        $RUN_POSTCOLORLAWFIT = "T";
        print "\t USING DEFAULT RUNJOB RUN_POSTCOLORLAWFIT $RUN_POSTCOLORLAWFIT\n";
        }
    if ( $RUN_MOVEMODEL eq '' ) {
        $RUN_MOVEMODEL = "T";
        print "\t USING DEFAULT RUNJOB RUN_MOVEMODEL $RUN_MOVEMODEL\n";
        }

    if ($RUN_PCAFIT eq "QUICKTEST") {
	$RUN_ERRORSNAKE = "F";
	$RUN_POSTCOLORLAWFIT = "F";
	$RUN_MOVEMODEL = "F";
	print "\t USING QUICKTEST MODE. ALL RUN_FLAGS SET TO FALSE\n";
    }
    
}

sub check_RUNMETHOD($)
{
    my ($inFile) = shift;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my $NARGALL = 99;
    my @tmp; ###
    my @tmp2;
    my $diestring; 
    my $key; ### consistent with rest of snana sntools.pl usage


    my $key = "NODELIST:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARGALL, $key, $OPTQUIET) ;
    @NODELIST = split(' ', @tmp[0]);
    if ( $NODELIST[0] eq '' ) { 
	$N_NODES = 0; 
    } else {
	print "\t NODELIST \t @tmp \n"; 
	$N_NODES = @NODELIST; 
    }

    my $key = "BATCH_INFO:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARGALL, $key, $OPTQUIET) ;
    if ( scalar(@tmp) > 0 ) {
	@tmp2 = split(' ', @tmp[0]);
	$BATCH_COMMAND  = $tmp2[0] ;
	$BATCH_TEMPLATE = $tmp2[1] ;
	$BATCH_NCORE    = $tmp2[2] ;
    } else {
	$BATCH_NCORE = 0;
    }

    $diestring = "Cannot specify both NODELIST and BATCH_INFO keywords.";
    if ($N_NODES > 0 && $BATCH_NCORE > 0 ) { die $diestring };

    if ( $BATCH_NCORE > 0 ){
	check_fileexist('BATCH_TEMPLATE', $BATCH_TEMPLATE);
    }
    

}

sub check_CONFIGPATHS()
{
    my $missing_die;

    ### CHECK INPUT FILE FOR USER RUNJOBS
    make_banner("CHECKING FOR RUN PATHS AND FILES") ;

    print "USER PATH AND FILE inputs are: \n";
    print "\t CONFIG_FILE \t\t $CONFIG_FILE \n";
    print "\t SALT2TRAININGPATH \t $SALT2TRAININGPATH \n";
    print "\t TRAIN_LIST \t\t $TRAIN_LIST \n";

    ### MUST have a configuration file to run any training job
    check_fileexist('CONFIG_FILE', $CONFIG_FILE);
    $missing_die = "Copy of $CONFIG_FILE to $SALTTRAIN_DIR failed. ABORT. \n"; 
    if ( $SALTTRAIN_DIR ne $local_dir ) { copy($CONFIG_FILE, $SALTTRAIN_DIR) or $missing_die };

    ### COPY training input file to each local directory
    if ( $SALTTRAIN_DIR ne $local_dir ) { copy($USER_inFile, $SALTTRAIN_DIR) or $missing_die };

    ### SALT2TRAININGPATH MUST be a legitimate directory location
    $missing_die = "ABORT: NO TRAINING SAMPLE DIRECTORY  found at $SALT2TRAININGPATH \n";
    check_direxist("SALT2TRAININGPATH", $SALT2TRAININGPATH);
}

sub check_TRAIN_LIST()
{
    my $missing_die;
    my $temp_file; 

    ### TRAIN_LIST MUST EXIST ###
    if ( $TRAIN_LIST ne "DEFAULT" ) {
	check_fileexist('TRAIN_LIST', $TRAIN_LIST);
	$missing_die = "Copy of $TRAIN_LIST to $SALTTRAIN_DIR failed. ABORT. \n"; 
	if ( $SALTTRAIN_DIR ne $local_dir ) { copy($TRAIN_LIST, $SALTTRAIN_DIR) or $missing_die };
    } elsif ( $N_GENVERSION == 1 ) {
	$temp_file = $SALT2TRAININGPATH . "/trainingsample.list";
	check_fileexist('TRAIN_LIST', $temp_file);
	$missing_die = "Copy of $TRAIN_LIST to $SALTTRAIN_DIR failed. ABORT. \n"; 
        copy($temp_file, $SALTTRAIN_DIR) or $missing_die ;
	$TRAIN_LIST = "trainingsample.list";  # RK - another hidden value
    } elsif ( $N_GENVERSION > 1 ) {
        $TRAIN_LIST = genversion_multi($N_GENVERSION, $GENVERSION);
    } else {
	die "Problem with check_CONFIGPATHS TRAIN_LIST checking loop. ABORT.\n";
    }

}

sub genversion_multi(@)
{

    # update such that if $SALTTRAIN_DIR ne $local_dir, 
    # trainlist_file = trainingsample.list ##

    my($n_genvers, $genvers) = @_;
    my $trainlist_filename;
    my $trainlist_file;
    my $i_genvers;
    my @genvers; 
    my $genvers_name;
    my $genvers_list;
    my @firstline;
    my @lastline;
    my (@lcprefix, $TMPSED, $str1, $str2 ) ;
    my @output; ##dump spot for system call
    my %firstcids;
    my %lastcids;
    my @firstsorted;
    my $sortflag = 0;
    my $die_string;

    if ( $SALTTRAIN_DIR eq $local_dir) {
	$trainlist_file = "combine_trainingsample.list";
	$trainlist_filename = $trainlist_file;
    } else {
	$trainlist_file = $SALTTRAIN_DIR . "/trainingsample.list";
	$trainlist_filename = "trainingsample.list"
    }

    print "\n\t MODEL: $MODELNAME , merging traininglist files. \n";
    print "\t Combined trainingsample.list file written to $trainlist_file \n\n";

    if ( -e $trainlist_file ) { unlink($trainlist_file);}

    $i_genvers = 0;
    @genvers = split(' ', $genvers);
    while ( $i_genvers < $n_genvers ){
	print "\t Adding $genvers[$i_genvers] to compiled list \n"; 
	$genvers_name = $genvers[$i_genvers];
	$genvers_list = $SALT2TRAININGPATH . "/" . $genvers_name . "/" .  "trainingsample.list";
	open (LIST, $genvers_list) || die "Unable to open list $genvers_list \n";
	@firstline = split(' ', <LIST>);
	close(LIST);
        @lcprefix = split('_', $firstline[2]);
	if ( $lcprefix[0] =~ m/\./) {
	    @lcprefix = split('-', $firstline[2]);
	    @lcprefix[0] = @lcprefix[0] . "-";
	} else {
	    @lcprefix[0] = @lcprefix[0] . "_";
	}
	@lastline = split(' ', `tail -n 1 $genvers_list`); 
        ### print "$firstline[2], $lcprefix[0], $firstline[0], $lastline[0] \n";
	$firstcids{$genvers_name} = $firstline[0];
	$lastcids{$genvers_name} = $lastline[0];

#        @output = `sed -e 's#'$lcprefix[0]'#'$genvers_name'/'$lcprefix[0]'#' -e 's#spectrum#'$genvers_name'/spectrum#' $genvers_list >> $trainlist_file `;

	# RK above substituion does not work on fulla ;
	#    so let's try another way below:
	$str1 = $lcprefix[0];  # RK
	$str2 = "spectrum" ;   # RK
	$TMPSED = "sed" ;      # RK
	$TMPSED = "$TMPSED " . "-e 's/$str1/$genvers_name\\/$str1/'" ; # RK
	$TMPSED = "$TMPSED " . "-e 's/$str2/$genvers_name\\/$str2/g'" ; # RK
        @output = `$TMPSED $genvers_list >> $trainlist_file ` ;   # RK

	$i_genvers = $i_genvers + 1;
    }

    ### CHECK cid ranges for overlaps
    @firstsorted = sort { $firstcids{$a} <=> $firstcids{$b} } keys %firstcids; 
    ###print "@firstsorted, $firstcids{$firstsorted[0]}, $firstcids{$firstsorted[1]}, $n_genvers \n"; 
    $i_genvers = 0;
    $die_string =  " have overlapping cid ranges. Fix sim offsets and try again. \n\n"  ;
    while ( $i_genvers < $n_genvers - 1 ) {
	$i_genvers = $i_genvers + 1;
	### print "$i_genvers\n"; 
	if ( $firstcids{$firstsorted[$i_genvers]} == $firstcids{$firstsorted[$i_genvers - 1]} ) { $sortflag = 1 } ;
	if ( $lastcids{$firstsorted[$i_genvers]} <= $lastcids{$firstsorted[$i_genvers-1]} )  { $sortflag = 1 } ;
	if ( $firstcids{$firstsorted[$i_genvers]} <= $lastcids{$firstsorted[$i_genvers-1]} )  { $sortflag = 1 } ;
	unless ( $sortflag == 0 )  { die "\nABORT: genversions $firstsorted[$i_genvers], $firstsorted[$i_genvers-1]" . $die_string } ; 
    }



    print "\t Merge complete. \n\n";
    return($trainlist_filename);

}


sub init_PATHS($)
{
    my ($inFile) = shift;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my @tmp; ###
    my $missing_die;
    my $key; ### consistent with rest of snana sntools.pl usage

    ### CHECK INPUT FILE FOR USER RUNJOBS
    make_banner("CHECKING FOR RUN PATHS AND FILES") ;
    print "USER PATH AND FILE inputs are: \n";

    ### MUST have a configuration file to run any training job
    my $key = "CONFIG_FILE:" ; ## REQUIRED
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTABORT) ;
    $CONFIG_FILE = @tmp[0];
    check_fileexist('CONFIG_FILE', $CONFIG_FILE);
    print "\t CONFIG_FILE \t $CONFIG_FILE \n"; 

    ### MUST have a SALT2TRAININGPATH to run any training job
    #
    # Checks for SALT2TRAININGPATH FIRST. 
    # Then checks for GENVERSION OR LCSIM_INPUT
    #
    my $key = "SALT2TRAININGPATH:" ; 
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $SALT2TRAININGPATH = @tmp[0];


    $missing_die = "ABORT: $inFile requires SALT2TRAININGPATH, GENVERSION, or LCSIM_INPUT key. \n";

    if ( $SALT2TRAININGPATH eq '' ) {
       my $key = "GENVERSION:";
       @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
       $GENVERSION = @tmp[0];

       if ( $GENVERSION eq '' ) {
           my $key = "LCSIM_INPUT:" ;
           @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
           $LCSIM_INPUT = @tmp[0];

           if ( $LCSIM_INPUT eq '' ) { die $missing_die; }

           check_fileexist("LCSIM_INPUT", $LCSIM_INPUT);

           my $key = "GENVERSION:" ;
           @tmp = sntools::parse_line($LCSIM_INPUT, $NARG1, $key, $OPTABORT) ;
           $GENVERSION = @tmp[0];
       }

       if ( $GENVERSION eq '' ) { die $missing_die; }

       $SALT2TRAININGPATH = $SNDATA_ROOT . "/SIM/SALT2-training/" . $GENVERSION;

    }

    $missing_die = "ABORT: NO TRAINING SAMPLE DIRECTORY  found at $SALT2TRAININGPATH \n"; 
    check_direxist("SALT2TRAININGPATH", $SALT2TRAININGPATH);
    print "\t SALT2TRAININGPATH \t $SALT2TRAININGPATH \n"; 


    my $key = "TRAIN_LIST:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    if ( @tmp[0] eq '' ) {
	$TRAIN_LIST = $SALT2TRAININGPATH . "/trainingsample.list";
	} else {
	    $TRAIN_LIST = @tmp[0];
	}
    check_fileexist("TRAIN_LIST", $TRAIN_LIST);
    print "\t TRAIN_LIST \t $TRAIN_LIST \n"; 

}

sub make_banner($)
{

    my $banner_string = shift;
    my $SYMLINE = "\n _____________________________ \n\n" ;

    print "$SYMLINE";
    print "$banner_string";
    print "$SYMLINE";

} # End of make_banner


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
    my $die_string1;
    my $die_string2;
    my $die_string;

    $die_string1 = "\n\nABORT: Script can't find file $file ( keyword: $key ). \n";
    $die_string2 = "       Check its location and try again. \n\n";
    $die_string = $die_string1 . $die_string2 ; 
    die $die_string unless -e $file;

} # End of check_fileexist

sub check_direxist(@)
{
    my($key, $dir) = @_;
    my $die_string;

    $die_string = "\n\nABORT: Cannot find $key directory $dir \n\n";
    die $die_string unless -d $dir;

} # End of check_direxist

sub check_badexit(@)
{
    ## ideally here, "key" is the name of the job whose log is being checked
    my $JOBname;
    my $logfile;
    my $ABORTkey;
    my @tmp; # for parse_line
    my $OPTQUIET = 2;
    my $ABORT;
    my $die_string;

    ($JOBname, $logfile, $ABORTkey) = @_;
    if ( $ABORTkey eq '' ) { $ABORTkey = "ABORT"; }

    @tmp = sntools::parse_line($logfile, 99, $ABORTkey, $OPTQUIET) ;
    $ABORT = @tmp; ### use number of entries in @tmp to test ABORT
    $die_string = "$JOBname job failed. ABORT. \n\n";
    if ( $ABORT > 0 ) {
        print "$ABORTkey string found in log $logfile .\n";
	die $die_string ;}

} # End of check_badexit

sub check_goodexit(@)
{
    ## ideally here, "key" is the name of the job whose log is being checked
    my $JOBname;
    my $logfile;
    my $SUCCESSkey;
    my @tmp; # for parse_line
    my $OPTQUIET = 2;
    my $SUCCESS;
    my $die_string;

    ($JOBname, $logfile, $SUCCESSkey) = @_;
    if ( $SUCCESSkey eq '' ) { $SUCCESSkey = "gracefully"; }

    @tmp = sntools::parse_line($logfile, 99, $SUCCESSkey, $OPTQUIET) ;
    $SUCCESS = @tmp; ### use number of entries in @tmp to test SUCCESS
    $die_string = "$JOBname job failed. ABORT.  \n\n";
    if ( $SUCCESS == 0 ) {
        print "$SUCCESSkey string not found in log $logfile .\n";
	die $die_string ;}

} # End of check_goodexit

sub run_QUICKTEST($)
{
    my ($runitfile) = shift;
    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n"; 

    ##if ( -e $SALT2MODEL_DIR && $SALT2MODEL_DIR ne "SALT2.Guy10") { qx(rm -r $SALT2MODEL_DIR); }
    unless ( -e $SALT2MODEL_DIR ) {
	print RUNIT "ln -s $SNANA_MODELPATH/SALT2.Guy10 $SALT2MODEL_DIR \n\n"; }
    close(RUNIT);
}

sub run_MOVEMODEL($)
{
    my ($runitfile) = shift;
    open(RUNIT, ">>$runitfile") or die "Can't open $runitfile \n";
    print RUNIT " $0 $USER_inFile MOVEMODEL $SALT2MODEL_DIR  >& movemodel.log \n\n";
    close(RUNIT); 

}

sub MOVEMODEL($)
{
### Pass code user input file. 
### Auto name model SALT2.GENMODEL 
###
### List of files to be moved is hard-coded. 
### Needs to be updated if files are added to snana implementation.
### Add make_README file function in future
###
my ($inFile) = shift;
my $orig_file;
my $new_file;  
my $nofile_die;
my $fail_message; 
my $file;

make_banner("MOVING MODEL TO SALT2MODEL_DIR $SALT2MODEL_DIR");

       ### GET USER SALT2INFO INPUTS
           init_SALT2INFO($inFile);

       ### CHECK FOR NECESSARY POST PROCESSING FILES
           $file = "salt2_color_correction_final.dat";
           print "\n\tCheck for $file \n\n";
           check_fileexist('COLOR LAW VERS 2', $file);

            print "Generate SNANA SALT2.INFO file. \n\n";
            makeSALT2INFO($file);

            $file = "SALT2.INFO";
            print "\n\tCheck for $file \n\n";
            check_fileexist('SNANA SALT2.INFO', $file);

            $file = "cle_final.list";
            print "\n\tCheck for $file \n\n";
            check_fileexist('SALT2 color dispersion', $file);

            print "Generate SNANA salt2_color_dispersion.dat \n";
            copy($file, "salt2_color_dispersion.dat") or die "salt2_color_dispersion.dat gen failed. \n";


       ### CLEAN UP OR CREATE SALT2MODEL DIR
       if ( -d $SALT2MODEL_DIR ) {
           print "\tCleaning up SALT2MODEL DIR \n";
           unlink(<$SALT2MODEL_DIR/*>);
       } else {
           print "\tCreating SALT2MODEL DIR \n";
           $nofile_die = "ABORT: Unable to create SALT2MODEL DIR. \n";
           mkdir($SALT2MODEL_DIR, 0777) || die $nofile_die;
       }

    $orig_file = 'salt2_color_correction.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\n\tMoved $orig_file\n";

    $orig_file = 'salt2_color_dispersion.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'SALT2.INFO';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_lc_dispersion_scaling.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_lc_relative_covariance_01.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_lc_relative_variance_0.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_lc_relative_variance_1.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_spec_covariance_01.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_spec_variance_0.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_spec_variance_1.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_template_0.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

    $orig_file = 'salt2_template_1.dat';
    $new_file = $SALT2MODEL_DIR . '/' . $orig_file;
    $fail_message = "Move of $orig_file to $new_file failed. ABORT. \n ";
    copy($orig_file, $new_file) or die $fail_message; 
    print "\tMoved $orig_file\n";

make_banner("MODEL MOVED SUCCESSFULLY");

}

sub init_SALT2INFO($)
{
    my ($inFile) = shift;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my @tmp; ###
    my $grid;
    my $key; ### consistent with rest of snana sntools.pl usage

    ### CHECK INPUT FILE FOR USER SALT2INFO PARAMS
    make_banner("CHECKING FOR USER-SELECTED SALT2INFO SETTINGS") ;
    print "SALT2INFO parameters are: \n";

    my $key = "SALT2_ERRMAPINTERPOPT:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    unless ( @tmp[0] eq '' ) { 
	$SALT2_ERRMAPINTERPOPT = @tmp[0];
	print "\t SALT2_ERRMAPINTERPOPT \t $SALT2_ERRMAPINTERPOPT \n"; }

    my $key = "SALT2_ERRMAPKCOROPT:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    unless ( @tmp[0] eq '' ) {
	$SALT2_ERRMAPKCOROPT = @tmp[0];
	print "\t SALT2_ERRMAPKCOROPT \t $SALT2_ERRMAPKCOROPT \n"; }

    my $key = "SALT2_MAGERRFLOOR:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    unless ( @tmp[0] eq '' ) { 
	$SALT2_MAGERRFLOOR = @tmp[0];
	print "\t SALT2_MAGERRFLOOR \t $SALT2_MAGERRFLOOR \n"; }

    ### CHECK GRID SPACING OF TEMPLATE FILES TO SET SALT2_SEDFLUXINTERPOPT
    $grid = check_grid("salt2_template_1.dat");
    if ( $grid < $GRIDLIM_SEDFLUXINTERPOPT ) { $SALT2_SEDFLUXINTERPOPT = 1; }
 
    
}

sub makeSALT2INFO($)
{
### First parse color law file, extracting necessary parts
### THIS SCRIPT ASSUMES COLOR LAW HAS FORMAT AS FOLLOWS:
###
### NPARAM
### PAR1
### PAR2
### ...
### PARN
### 
### MAG_OFFSET: 0.27 ## TO MATCH OVERALL FLUX SCALE
###
### Salt2ExtinctionLaw.version ELVERS
### Salt2ExtinctionLaw.min_lambda MINLAM
### Salt2ExtinctionLaw.max_lambda MAXLAM
### (EOF)
###
### IF FORMAT CHANGES, THIS CODE MUST BE REWRITTEN

### Then write out SALT2.INFO file using information from color law file


    my ( $colorfile ) = shift;
    my ( $salt2infofile ) = "SALT2.INFO" ;
    my ( @line, $iline, $NPARAM, @PARS );
    my ( $ELVERS, $MINLAM, $MAXLAM );
    my ( $ABORT_die);

    ### GET SALT2INFO INFO FROM INPUT FILE

    ### FIGURE OUT NPARAM IN COLORFILE
    open ( COLORFILE, $colorfile ) || die "Couldn't open COLORFILE $colorfile \n";
    $iline = 0;
while(<COLORFILE>) {
    chomp;
    @line = split(' ', $_);
    if ( $iline == 0 ) {
        print "first line of color file is @line \n";
        $NPARAM = $line[0];
        print " Color Law NPARAM is $NPARAM \n";
    } elsif ( $iline > 0 && $iline <= $NPARAM  ) {
        ### GET PARAMS UNTIL HAVE NPARAM OF THEM
        push( @PARS, $line[0]);
    } else {
        ### GET PARAMETERS MINLAM, MAXLAM, ELVERS
        if ( $line[0] =~/^Salt2ExtinctionLaw.min_lambda/ ) { $MINLAM = $line[1]; }
        if ( $line[0] =~/^Salt2ExtinctionLaw.max_lambda/ ) { $MAXLAM = $line[1]; }
        if ( $line[0] =~/^Salt2ExtinctionLaw.version/ ) { $ELVERS = $line[1]; }
    }
    $iline += 1;
}
close(COLORFILE);

print "PARS ARE : @PARS \n";
print "ELVERS: $ELVERS, MINLAM: $MINLAM, MAXLAM: $MAXLAM \n";


    ### WRITE PARAMETERS TO SALT2.INFO FILE. MOVE OLD VERSION IF NECESSARY
if ( -e "SALT2.INFO" ) { move("SALT2.INFO", "origSALT2.INFO")}
print "will write results to output file $salt2infofile \n";

    $ABORT_die = "ABORT: Can't open new SALT2.INFO file for writing \n";
    open(SALT2INFO, ">$salt2infofile") or die $ABORT_die;
    printf SALT2INFO "RESTLAMBDA_RANGE: %d %d \n", $MINLAM, $MAXLAM;
    printf SALT2INFO "COLORLAW_VERSION: %s \n\n", $COLORLAW_VERSION;
    printf SALT2INFO "COLORLAW_PARAMS: %d %d %d ", $MINLAM, $MAXLAM, $NPARAM;
    printf SALT2INFO "%0.6f %0.6f %0.6f %0.6f \n", $PARS[0]*1, $PARS[1]*1, $PARS[2]*1, $PARS[3]*1;
    printf SALT2INFO "COLOR_OFFSET: %0.1f \n\n", $COLOR_OFFSET;
    printf SALT2INFO "MAG_OFFSET: %0.3f \n\n", $SALT2MAG_OFFSET;
    printf SALT2INFO "SEDFLUX_INTERP_OPT: %d \n", $SALT2_SEDFLUXINTERPOPT;
    printf SALT2INFO "ERRMAP_INTERP_OPT: %d \n", $SALT2_ERRMAPINTERPOPT;
    printf SALT2INFO "ERRMAP_KCOR_OPT: %d \n\n", $SALT2_ERRMAPKCOROPT;
    printf SALT2INFO "MAGERR_FLOOR: %0.5f \n", $SALT2_MAGERRFLOOR;
    close(SALT2INFO);

}

sub check_grid($) {


    my $file = shift;

    my @line1;
    my @line2;
    my $grid;

    open (TEMPLATE, $file) || die "Unable to open template $file \n";
    @line1 = split(' ', <TEMPLATE>);
    @line2 = split(' ', <TEMPLATE>);
    $grid = $line2[1]-$line1[1];
	
    close(TEMPLATE);

    return($grid);

}

sub parse_configFile_SALT2train($)
{
    my $config_file = shift;

    my $key;
    my @tmp;
    my $NARG1 = 1;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;

    make_banner("PARSING CONFIG FILE FOR PCAFIT KEYS" );

    my $key = "\@nfit_iteration" ; ## OPTIONAL
    @tmp = sntools::parse_line($config_file, $NARG1, $key, $OPTWARN) ;
    if ( @tmp[0] eq '') { 
	$OPT_PCAFIT = $RUN_PCAFIT;
    } elsif (@tmp[0] ne '' and $RUN_PCAFIT eq "QUICKTEST") {
	$OPT_PCAFIT = $RUN_PCAFIT;
	print "\t OVERRIDING LOCAL nfit_iteration @tmp. \n";
	print "\t USING QUICKTEST mode. \n";
	print "\t NO TRAINING WILL BE RUN. \n";
    } else {
	$OPT_PCAFIT = @tmp[0];
	print "\t OVERRIDING GLOBAL RUN_PCAFIT \n";
	print "\t USING config key nfit_iteration $OPT_PCAFIT \n"; 
    }

    my $key = "\@write_lambin" ; ## OPTIONAL
    @tmp = sntools::parse_line($config_file, $NARG1, $key, $OPTQUIET) ;
    ##print " tmp vector is @tmp \n" ;
    unless ( @tmp[0] eq '' ) { 
	$WRITE_LAMBIN = @tmp[0];
	print "\t WRITE_LAMBIN \t $WRITE_LAMBIN \n"; }

    my $key = "\@write_phasebin" ; ## OPTIONAL
    @tmp = sntools::parse_line($config_file, $NARG1, $key, $OPTQUIET) ;
    ##print " tmp vector is @tmp \n" ;
    unless ( @tmp[0] eq '' ) { 
	$WRITE_PHASEBIN = @tmp[0];
	print "\t WRITE_PHASEBIN \t $WRITE_PHASEBIN \n"; }
    
}

sub parse_inFile_SALT2train($)
{
    use Switch;

    my $inFile_SALT2train = shift;
    my ( @line ) ; ## use to hold lines in file
    my $narg;
    my $iarg;
    my $default_modelname;
    my $die_string1;
    my $die_string2;
    my $die_string;
    my $outdir_okcheck;
    my @path_parse;
    my $key;
    my @tmp;
    my $NARG1 = 1;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $GEN_CONFIGFILE;


    ### ENSURE inFILE existence HAS ALREADY BEEN CHECKED

    ### Certain Format Assumptions are being made
    ### 1) CONFIG_FILE should be the first non-job key listed
    ### 2) USE ONE OF THE FOLLOWING keys per train: 
    ###           GENVERSION, SALT2TRAININGPATH, OR LCSIM_INPUT.
    ###           Only GENVERSION allows multiple values. 
    ###           Script expects SALT2TRAININGPATH, LCSIM_INPUT key lines to have format:
    ###           key, path/input file (optional: , modelname key, modelname).
    ###           More than four values on these lines will trigger ABORT. 
    ### 3) SHOULD ENTER MODELNAME AS LAST ENTRY ON GENVERSION/SALT2TRAININGPATH/LCSIM_INPUT LINE
    ###           preceded by keyword <MODELNAME:>
    ###           If don't enter, will use default. However !!  
    ###           the existence of multiple models with same name will trigger ABORT.  
    ### 4) IF USING USER-SELECTED TRAIN_LIST, TRAIN_LIST KEY MUST DIRECTLY FOLLOW
    ###           GENVERSION/SALT2TRAININGPATH/LCSIM_INPUT KEY 
    ### 5) IF A PARTICULAR MODELNAME HAS NO TRAIN_LIST FILE
    ###           AND N_GENVERSION = 1, DEFAULT WILL BE USED ( pcafit handles this case )
    ###           AND N_GENVERSION > 1, DEFAULT TRAIN_LIST FILE WILL BE GENERATED IN TRAIN DIRECTORY

    open (INFILE_SALT2TRAIN, $inFile_SALT2train) or die "ERROR: Cannot open input file $inFile_SALT2train \n";
    make_banner("Parsing input file: $inFile_SALT2train");

    while ( <INFILE_SALT2TRAIN> ){
        chomp;
        @line = split(' ', $_ );
        switch ($line[0]){
	    case("SALT_TEMPLATE0:"){
		$SALT_TEMPLATE0 = $line[1];
		print "USING INITIAL SALT_TEMPLATE0: $SALT_TEMPLATE0 \n";
	    }
	    case("SALT_TEMPLATE1:"){
		$SALT_TEMPLATE1 = $line[1];
		print "USING INITIAL SALT_TEMPLATE1: $SALT_TEMPLATE1 \n";
	    }
            case("CONFIG_FILE:")     {
                $CONFIG_FILE = $line[1];
            }
            case("GENVERSION:")    {
		$GENVERSION = '';
                $MODELNAME = '';
                $narg = @line; ## pull out number of elements in line vector
		### check for MODELNAME ###
		$default_modelname = "";
		$iarg = 1;
                $N_GENVERSION = 0;
		while ( $iarg < $narg ){
		    if ( $line[$iarg] eq 'MODELNAME:' || $line[$iarg] eq 'MODELNAME' ) { 
			$MODELNAME = $line[$iarg+1];
			print "MODELNAME is $MODELNAME \n";
			$iarg = $iarg + 1;
		    } else {
		    $N_GENVERSION = $N_GENVERSION + 1;
		    if ( $N_GENVERSION > 1 ) { $GENVERSION = "$GENVERSION".' '; }
		    $GENVERSION = $GENVERSION.@line[$iarg];
		    $default_modelname = "${default_modelname}+${line[$iarg]}";
		}
		    $iarg = $iarg + 1;
		}
		if ($N_GENVERSION == 1 ) {
		    $SALT2TRAININGPATH = $SNDATA_ROOT . "/SIM/SALT2-training/" . $GENVERSION;
		}
		if ($N_GENVERSION > 1 ){
		    $SALT2TRAININGPATH = $SNDATA_ROOT . "/SIM/SALT2-training" ;
		}
		if ( $MODELNAME eq '' ) {
		    $default_modelname = substr $default_modelname, 1;
		    $MODELNAME = $default_modelname;
		}
		### TEMP KEEP FOR SYNTAX:    $TRAIN_LIST = genversion_multi($ntmp, @tmp_split); ###
		### NOW HAVE MODELNAME, GENVERSION, N_GENVERSION, and SALT2TRAININGPATH ###
		$die_string = "For input line @line, MODELNAME $MODELNAME already exists. ABORT.\n";
		if ( exists $TRAIN_PATH{$MODELNAME} ) { die $die_string ; }
		$TRAIN_PATH{$MODELNAME} = $SALT2TRAININGPATH;
		$TRAIN_CONFIG{$MODELNAME} = $CONFIG_FILE;
		$TRAIN_GENVERSION{$MODELNAME} = $GENVERSION;
		$TRAIN_NGENVERSION{$MODELNAME} = $N_GENVERSION;
	    }
	    case("LCSIM_INPUT:")  {
		$GENVERSION = '';
		$MODELNAME = '';
		$narg = @line; ## pull out number of elements in line vector
		### expect $narg to be 2-4 ( key, file name, modelname key, modelname ) 
		$die_string = "Unexpected syntax for line: LCSIM_INPUT . ABORT. \n";
		unless ( $narg > 1 && $narg < 5 ) { die $die_string; } 
		### check for MODELNAME ###
                $iarg = 1;
		$N_GENVERSION = 0;
                while ( $iarg < $narg ){
		    if ( $line[$iarg] eq 'MODELNAME:' || $line[$iarg] eq 'MODELNAME' ) { 
                        $MODELNAME = $line[$iarg+1];
			$iarg = $iarg + 1;
                    } else {
			$N_GENVERSION = $N_GENVERSION + 1;
			$LCSIM_INPUT = $line[$iarg];
		    }
		    $iarg = $iarg + 1;
                }
                if ($N_GENVERSION == 1 ) {
                    ### parse LCSIM_INPUT file for GENVERSION ###
                    check_fileexist("LCSIM_INPUT", $LCSIM_INPUT);
                    my $key = "GENVERSION:" ;
                    @tmp = sntools::parse_line($LCSIM_INPUT, $NARG1, $key, $OPTABORT) ;
                    $GENVERSION = @tmp[0];
                    ### done parsing ###
                    if ( $GENVERSION eq '' ) {die "GENVERSION key not found in file $LCSIM_INPUT \n";}
		    $MODELNAME = $GENVERSION;
                    $SALT2TRAININGPATH = $SNDATA_ROOT . "/SIM/SALT2-training/" . $GENVERSION;
		}
		### NOW HAVE MODELNAME, GENVERSION, N_GENVERSION, and SALT2TRAININGPATH ###
		$die_string = "For input line @line, MODELNAME $MODELNAME already exists. ABORT.\n";
		if ( exists $TRAIN_PATH{$MODELNAME} ) {die $die_string ;}
		$TRAIN_PATH{$MODELNAME} = $SALT2TRAININGPATH;
		$TRAIN_CONFIG{$MODELNAME} = $CONFIG_FILE;
		$TRAIN_GENVERSION{$MODELNAME} = $GENVERSION;
                $TRAIN_NGENVERSION{$MODELNAME} = $N_GENVERSION;
	    }
	    case("SALT2TRAININGPATH:")  {
		$GENVERSION = '';
		$MODELNAME = '';
		$narg = @line; ## pull out number of elements in line vector
		### expect $narg to be 2-4 ( key, path, modelname key, modelname ) 
		$die_string = "Unexpected syntax for line: SALT2TRAININGPATH . ABORT. \n";
		unless ( $narg > 1 && $narg < 5 ) { die $die_string; } 
		### check for MODELNAME ###
                $iarg = 1;
		$N_GENVERSION = 0;
                while ( $iarg < $narg ){
		    if ( $line[$iarg] eq 'MODELNAME:' || $line[$iarg] eq 'MODELNAME' ) { 
                        $MODELNAME = $line[$iarg+1];
			$iarg = $iarg + 1;
                    } else {
			$N_GENVERSION = $N_GENVERSION + 1;
			$SALT2TRAININGPATH = $line[$iarg];
		    }
		    $iarg = $iarg + 1;
                }
		if ( $MODELNAME eq '' ) {
		    @path_parse = split('/', $SALT2TRAININGPATH);
		    $MODELNAME = $path_parse[@path_parse-1];
		}
		$GENVERSION = $MODELNAME;
		### NOW HAVE MODELNAME, GENVERSION, N_GENVERSION, and SALT2TRAININGPATH ###
		$die_string = "For input line @line, MODELNAME $MODELNAME already exists. ABORT.\n";
		if ( exists $TRAIN_PATH{$MODELNAME} ) {die $die_string ;}
		$TRAIN_PATH{$MODELNAME} = $SALT2TRAININGPATH;
		$TRAIN_CONFIG{$MODELNAME} = $CONFIG_FILE;
		$TRAIN_GENVERSION{$MODELNAME} = $GENVERSION;
                $TRAIN_NGENVERSION{$MODELNAME} = $N_GENVERSION;
	    }
	    case("TRAIN_LIST:") {
		$TRAIN_LIST = @line[1];
		$die_string = "No MODELNAME currently defined for TRAIN_LIST $TRAIN_LIST. Check input file. ABORT.\n";
		unless ( exists $TRAIN_PATH{$MODELNAME} ) { die $die_string; }
		$TRAIN_TRAINLIST{$MODELNAME} = $TRAIN_LIST; 
		## IF A PARTICULAR MODELNAME HAS NO LIST FILE, AND N_GENVERSION = 1, DEFAULT WILL BE USED ##
		$TRAIN_LIST = '';
	    }
	    else {
	    }
        }
    }
    close(INFILE_SALT2TRAIN);

    make_banner("READING THROUGH HASHES");
    $iarg = 1;

    print "\n\t\tTHE FOLLOWING TRAININGS WILL BE RUN: \n";
    while ( ($MODELNAME, $CONFIG_FILE) = each(%TRAIN_CONFIG)) {
	print "\n\t\t  MODELNAME: $MODELNAME\n";
	print "\t\t  CONFIG_FILE: $TRAIN_CONFIG{$MODELNAME}\n";
	print "\t\t  GENVERSION: $TRAIN_GENVERSION{$MODELNAME}\n";
	print "\t\t  N_GENVERSION: $TRAIN_NGENVERSION{$MODELNAME}\n";
	print "\t\t  TRAIN PATH: $TRAIN_PATH{$MODELNAME}\n";
	if ( exists $TRAIN_TRAINLIST{$MODELNAME} )  {
	    $TRAIN_LIST = $TRAIN_TRAINLIST{$MODELNAME};
	} else {
	    $TRAIN_LIST = "DEFAULT" ;
	    $TRAIN_TRAINLIST{$MODELNAME} = "DEFAULT"; 
	}
	print "\t\t  TRAIN LIST: $TRAIN_TRAINLIST{$MODELNAME}\n";
        unless ( $NODELIST[0]  eq '' ) {
	    $TRAIN_NODELIST{$MODELNAME} = $NODELIST[($iarg%$N_NODES)-1];
	    print "\t\t  NODE: $TRAIN_NODELIST{$MODELNAME}\n";

        }
        $iarg = $iarg + 1;
    }
    $N_TRAINS = $iarg-1; 

    $die_string = "Must give multiple cores/nodes for multiple trainings. ABORT.\n";
    if ( $N_TRAINS > 1 && $N_NODES < 2 && $BATCH_NCORE < 2 ) {die $die_string;}

    if ($N_NODES > 0 ){
	print "N_NODES    = $N_NODES \n";
	print "NODELIST   = @NODELIST \n";
    }

    if ($BATCH_NCORE > 0){
	my $BCMD = "$BATCH_COMMAND $BATCH_TEMPLATE";
	my $B3   = "$BATCH_NCORE cores";
	print "\n\tBATCH_COMMAND: $BCMD ($B3) \n";
    }
    
    my $key = "OUTDIR:" ; ## REQUIRED IF N_TRAINS > 1;
    @tmp = sntools::parse_line($inFile_SALT2train, $NARG1, $key, $OPTQUIET) ;
    $OUTDIR = @tmp[0];
    $die_string1 = "\nWhen running multiple trainings, \n";
    $die_string2 = " absolute path to output directory must be given with input key OUTDIR: . ABORT . \n";
    $die_string = $die_string1 . $die_string2 ;
    if ( $OUTDIR eq '' && $N_TRAINS > 1 ) { 
	die $die_string ; 
     } elsif ($OUTDIR eq '' && $N_TRAINS == 1 ) {
	$OUTDIR = $local_dir;
	print "\n\t\tTraining will be carried out in current directory. \n" ;
    } elsif ($N_TRAINS > 1 ) {
	print "\n\t\t$N_TRAINS trainings will be spawned from directory $OUTDIR.  \n" ;
	print "\n\t\tThis directory will be wiped before proceeding. \n\n" ;
	##$outdir_okcheck = promptUser("Is this ok? Enter Y or N.\n", "N");
	##$die_string = "\n\n\t\tUser selected N. Output dir will not be wiped. ABORT.\n ";
	##unless ( $outdir_okcheck eq "Y" ) { die $die_string ; }
    }


}

sub run_TRAINING() 
{
    my $source_file;
    my $batch_file;
    my $nofile_die;
    my @output;
    my $cmdSpawn;

    make_banner("PREPARING TO RUN TRAINING");

    if ( ($N_NODES > 0 || $BATCH_NCORE > 0) && $N_TRAINS > 1) {
	print "\nPreparing output directory $OUTDIR \n";
	if ( -e $OUTDIR ) {
	    print "\tRemoving existing output directory: $OUTDIR \n";
	    @output = qx(rm -r $OUTDIR);
	} 
	print "\tCreating output directory $OUTDIR \n";
	$nofile_die = "ABORT: Unable to create directory $OUTDIR. \n";
	mkdir($OUTDIR, 0777) || die $nofile_die;
    }

    print "N_TRAINS is $N_TRAINS \n";

    while ( ($MODELNAME, $CONFIG_FILE) = each(%TRAIN_CONFIG)) {
        # reload variables pertaining to train being run
       	$N_GENVERSION = $TRAIN_NGENVERSION{$MODELNAME};
	$GENVERSION = $TRAIN_GENVERSION{$MODELNAME};
	$SALT2TRAININGPATH = $TRAIN_PATH{$MODELNAME};
	$TRAIN_LIST = $TRAIN_TRAINLIST{$MODELNAME};

	parse_configFile_SALT2train($CONFIG_FILE); #particulars on pcafit

	make_banner("PREPARE FOR MODEL: $MODELNAME");
	
	if ( ($N_NODES > 0 || $BATCH_NCORE >0) && $N_TRAINS >1 )  {
	    $SALTTRAIN_DIR = $OUTDIR . '/' . $MODELNAME;
	    if ( -e $SALTTRAIN_DIR ) {
		print "\tRemoving existing directory $SALTTRAIN_DIR \n"; 
		@output = qx(rm -r $SALTTRAIN_DIR);
	    } 
	    print "\tCreating directory $SALTTRAIN_DIR \n"; 
	    $nofile_die = "ABORT: Unable to create directory $SALTTRAIN_DIR. \n";
            mkdir($SALTTRAIN_DIR, 0777) || die $nofile_die;
	    ### $TRAIN_NODE = $TRAIN_NODELIST{$MODELNAME}; ###
	} else {
	    $SALTTRAIN_DIR = $local_dir;
	}
	
	check_CONFIGPATHS();##copies config to trainingdir, checks training path
	check_TRAIN_LIST();##checks training list, merges multiples to trainingdir if nec
        setup_ENV($USER_inFile);##determining env settings and modeldir

	if ($N_TRAINS > 1 && $N_NODES > 0 ) {

	    $source_file = $SALTTRAIN_DIR . "/" . $MODELNAME . ".runit";
	    ## no batch file	    
	    make_RUN($source_file);

	} elsif ($BATCH_NCORE > 0 ){

	    $source_file = $SALTTRAIN_DIR . "/" . $MODELNAME . "_batch";
	    $batch_file  = $SALTTRAIN_DIR . "/" . $MODELNAME . ".batch";
	    
	    make_batchhead($batch_file);
	    make_RUN($batch_file); ### appends RUN commands to end of batchfile
	    make_batchsource($source_file, $batch_file); ### makes _batch file to be sourced
	    
	} else {
	    $source_file = $SALTTRAIN_DIR . "/" .  $MODELNAME .  ".runit";
	    make_RUN($source_file);## writes runit script for training
	}
	
	qx(chmod +x $source_file); # RK
       submit_RUN($source_file); ### submits job

    }

}


sub check_DONE()
{
    my @output;
    my @tmp;
    my $key;
    my $NARG1 = 1;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $file = "BUSY";
    my $busyfile;

    if ( $SPAWN_TRAINS > 1 ) {
	@output = `find ../ -name $DONE_LOCAL | wc | awk '{print $1}'`;
    } else {
	@output = `find ./ -name $DONE_LOCAL | wc | awk '{print $1}'`;
    }
    @tmp = split(' ', @output[0]);
    $NDONE_CURRENT = @tmp[0];

    if ( $NDONE_CURRENT == $SPAWN_TRAINS ) {
	print "\t in check_DONE(): ALL $SPAWN_TRAINS done files have been found.\n";
	$DONEFLAG = 1;
    } else {
	print "\t in check_DONE(): $NDONE_CURRENT of  $SPAWN_TRAINS done files have been found.\n";
	$DONEFLAG = 0;
    }


    my $key = "BUSYFILE:" ; ## OPTIONAL 
    @tmp = sntools::parse_line($file, $NARG1, $key, $OPTQUIET) ;
    $busyfile = @tmp[0];
    unlink($busyfile);
    unlink($file); 
    
    return $DONEFLAG;
     
}

sub if_DONE()
{


    my @output;

    @output = qx(touch $DONE_STAMP);
}

sub make_batchsource(@)
{
    my ($batchsourcefile, $batchfile) = @_;
    my $batch_stem;

    $batch_stem = $MODELNAME . ".batch";
    
    open(RUNIT, ">$batchsourcefile") or die  "Can't open $batchsourcefile \n";
    print RUNIT "cd $SALTTRAIN_DIR \n\n";
    print RUNIT "$BATCH_COMMAND $batch_stem \n";
    close(RUNIT);
    
}

sub make_batchhead($)
{
    my ($batchfile) = shift;
    my ($tmp, $sedcmd);
    my ($batch_stem, $log_stem);

    $batch_stem = $MODELNAME . ".batch";
    $log_stem   = $MODELNAME . ".log";

    # convert BATCH_TEMPLATE into $batchFile by substituting
    # the REPLACE keys in the TEMPLATE

    # prepare sed command to replace keys in batch-template file.
    $sedcmd = "sed" ;
    $sedcmd = "$sedcmd " . "-e 's/REPLACE_NAME/$batch_stem/g'";
    $sedcmd = "$sedcmd " . "-e 's/REPLACE_LOGFILE/$log_stem/g'";
    $sedcmd = "$sedcmd " . "-e 's/REPLACE_MEM/$BATCH_MEM/g'";
    $sedcmd = "$sedcmd " . "-e 's/REPLACE_JOB//g'";

    # replace keys and make new batch file.
    print "debug: $sedcmd $BATCH_TEMPLATE $batchfile \n";
    qx($sedcmd $BATCH_TEMPLATE > $batchfile ; chmod +x $batchfile);  # RK

}


sub make_RUN($)
{

    my $source_file = shift; ## passes in name of source_file
    my $touch_file = $SALTTRAIN_DIR . "/" . $DONE_LOCAL;
    my $q= "";
    my ($cmdSpawn, $ENVDEF); 
    
    print "\t Creating runit file $source_file \n";
    open(RUNIT, ">>$source_file") or die  "Can't open $source_file \n";
    print RUNIT "cd $SALTTRAIN_DIR \n\n";

    
    $ENVDEF = sntools::setEnvString($SHELL_NAME, "SALT2_DIR", $SALT2_DIR) ;
    print RUNIT "$ENVDEF \n";   # RK

    $ENVDEF = sntools::setEnvString($SHELL_NAME, "SALTPATH", $SALTPATH);  
    print RUNIT "$ENVDEF \n";   # RK

    $ENVDEF = sntools::setEnvString($SHELL_NAME, "SALT2TRAININGPATH", 
				    $SALT2TRAININGPATH);
    print RUNIT "$ENVDEF\n";  # RK

    $ENVDEF = sntools::setEnvString($SHELL_NAME, "SNANA_MODELPATH", 
				    $SNANA_MODELPATH);
    print RUNIT "$ENVDEF\n" ;  # RK


    #print RUNIT "setenv SALT2_DIR $SALT2_DIR \n";
    #print RUNIT "setenv SALTPATH $SALTPATH \n"; 
    #print RUNIT "setenv SALT2TRAININGPATH $SALT2TRAININGPATH \n";
    #print RUNIT "setenv SNANA_MODELPATH $SNANA_MODELPATH \n";

    print RUNIT " \n";
    print RUNIT "hostname >&  train.log \n";

    if ( $N_NODES > 0 ) {   # RK
	print RUNIT "cat /proc/meminfo | grep MemFree >>& train.log \n"; 
    }

    print RUNIT " \n"; 
    close(RUNIT); 

    make_BUSY($source_file);

    if ( $OPT_PCAFIT eq "1" || $OPT_PCAFIT eq "2") { 
	run_PCAFIT($source_file, 0);
    }
    if ( $RUN_ERRORSNAKE eq "T" || $OPT_PCAFIT eq "2" ) { 
    	run_ERRORSNAKE($source_file, 0); 
    }
    if ( $OPT_PCAFIT         eq "2"){
	prep_PCAFIT_WERRSNAKE($source_file); 	       
	run_PCAFIT($source_file, 1);
	run_ERRORSNAKE($source_file, 1);
    }
    if ( $RUN_POSTCOLORLAWFIT eq "T") { run_POSTCOLORLAWFIT($source_file); }
    if ( $RUN_MOVEMODEL       eq "T") { run_MOVEMODEL($source_file);}

    if ( $OPT_PCAFIT eq "QUICKTEST") {run_QUICKTEST($source_file);}

    open(RUNIT, ">>$source_file") or die "Can't open $source_file \n";
    print RUNIT " touch $touch_file \n";
    print RUNIT " $0 $USER_inFile CHECKDONE $N_TRAINS $DONE_STAMP >& checkdone.log \n";
    close(RUNIT);

}

sub batch_RUN($)
{
    my $source_file = shift; ### passing in source_file name

    my $q= "";
    my $cmdSpawn; 

    print "\t Launching BATCH job $source_file \n";
    $cmdSpawn = "ssh -x $TRAIN_NODELIST{$MODELNAME} ${q} source $source_file ";
    system("$cmdSpawn &");
    make_banner(" Running $source_file on node $TRAIN_NODELIST{$MODELNAME} .");
}

sub submit_RUN($)
{
    my $source_file = shift; ### passing in source_file name
    my $source_stem;

    my $q= "";
    my $cmdSpawn;


    if ($N_TRAINS > 1 && $N_NODES > 0 ) {
        print "\t Launching SSH job $source_file \n";
        $cmdSpawn = "ssh -x $TRAIN_NODELIST{$MODELNAME} ${q} source $source_file ";

    } elsif ($BATCH_NCORE > 0) {
	$source_stem = $MODELNAME . "_batch";
        $cmdSpawn = "cd $SALTTRAIN_DIR ; ./$source_stem";

    } else {
	$source_stem = $MODELNAME . ".runit";
	$cmdSpawn = "source $source_stem ";
    }

    print ("\t submit_RUN:  $cmdSpawn \n");
    system("$cmdSpawn &");

}



sub ssh_RUN($)
{
    my $source_file = shift; ### passing in source_file name

    my $q= "";
    my $cmdSpawn; 

    print "\t Launching SSH job $source_file \n";
    $cmdSpawn = "ssh -x $TRAIN_NODELIST{$MODELNAME} ${q} source $source_file ";
    system("$cmdSpawn &");
    make_banner(" Running $source_file on node $TRAIN_NODELIST{$MODELNAME} .");
}

sub spawn_RUN()
{
    my $input_file = $SALTTRAIN_DIR . "/train.input"; 
    my $touch_file = $SALTTRAIN_DIR . "/" . $DONE_LOCAL;
    my $source_file = $SALTTRAIN_DIR . "/" . $MODELNAME .  ".runit";
    my $q= "";
    my ($cmdSpawn, $ENVDEF); 
    
    
    print "\t Creating spawn input file $input_file \n";
    open(SPAWNINPUT, ">$input_file") or die "Can't open $input_file \n";
    write(SPAWNINPUT);
    close(SPAWNINPUT);

    print "\t Creating spawn runit file $source_file \n";
    open(RUNIT, ">$source_file") or die  "Can't open $source_file \n";

    $ENVDEF = sntools::setEnvString($SHELL_NAME, "SALT2_DIR", $SALT2_DIR) ;
    print RUNIT "$ENVDEF \n";

    $ENVDEF = sntools::setEnvString($SHELL_NAME, "SALTPATH", $SALTPATH);  
    print RUNIT "$ENVDEF \n"; 

    $ENVDEF = sntools::setEnvString($SHELL_NAME, "SNANA_MODELPATH", 
				    $SNANA_MODELPATH);
    print RUNIT "$ENVDEF\n" ;

#    print RUNIT "setenv SALT2_DIR $SALT2_DIR \n";
#    print RUNIT "setenv SALTPATH $SALTPATH \n"; 
#    print RUNIT "setenv SNANA_MODELPATH $SNANA_MODELPATH \n";

    print RUNIT " \n";
    print RUNIT "cd $SALTTRAIN_DIR \n";
    print RUNIT "hostname >>&  train.log \n";
    print RUNIT "cat /proc/meminfo | grep MemFree >>& train.log \n";
    print RUNIT " \n"; 
    print RUNIT  "$0 train.input >>& train.log \n";
    print RUNIT " touch $touch_file \n"; 
    print RUNIT " $0 train.input CHECKDONE $N_TRAINS $DONE_STAMP >& checkdone.log \n";
    close(RUNIT);

    print "\t Preparing to source ${local_dir}/$source_file \n";
    $cmdSpawn = "ssh -x $TRAIN_NODELIST{$MODELNAME} ${q} source ${local_dir}/$source_file ";
    system("$cmdSpawn &");
    make_banner(" Running $source_file on node $TRAIN_NODELIST{$MODELNAME} .");

}

sub promptUser(@) {

    my ($promptString, $defaultValue) = @_ ;

    if ($defaultValue) {
	print $promptString, "[", $defaultValue, "]: ";
    } else {
	print $promptString, ": ";
    }

    $| = 1;               # force a flush after our print
    $_ = <STDIN>;         # get the input from STDIN 

    chomp; # remove newline from user input

    if ("$defaultValue") {
	return $_ ? $_ : $defaultValue;    # return $_ if it has a value
    } else {
	return $_;
    }
}

