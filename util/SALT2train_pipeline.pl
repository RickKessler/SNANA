#!/usr/bin/perl
#  
# Created May 2011 by R.Kessler
#
# Pipeline script to run end-to-end test for
# the SALT2 training:
#
# 1. Create simulated training sample in SNANA format
#
# 2. Translate into SALT2 format (SALT2train_translate.pl)
#
# 3. Copy any included actual training sets into 
#     $SNDATA_ROOT/SIM/SALT2-training directory
#     for combination with simulated training sets
#
# 4. use SALT2train_run.pl to run the training(s) and to
#    create an SNANA-SALT2 version
#    [for DATA option, start here]
#
# 5. Generate independent simulated data sets such as
#    SDSS, SNLS, Nearby. Use the same SIMSED model
#    as in step 1.
#
# 6. Fit each simulated sample (step 5) with the
#    SALT2-trained model(s) from step 4.
#
# 7. Run SALT2mu on the fitted outputs to get distances.
#
# 8. Run cosmo fitter on the SALT2mu distances to get cosmology.
#
# If DMCORR keywords are set:
#    9.  Generate additional sims with which to perform bias corrections. 
#        Sims based either on the input model, or on the trained models.
#    10. Fit the sims with the trained models.
#    11. Run SALT2mu on the sims to get distances.
#    12. Use the sims to correct the simulated data from steps 5-7
#    13. Run cosmo fitter on the corrected distances to get cosmology
#
#
# 14. Make MU-residual plot for each training test
#     Include the SALT2-model uncertainties.
#
# 15. purge/delete files
#
# 16. archive files
# 
# Standard usage
#  SALT2train_pipeline.pl <inputFile>
#  SALT2train_pipeline.pl <inputFile>  RUN  &
#
# Optional usage
#
#  SALT2train_pipeline.pl <inputFile>  LASTSTAGE $ISTAGE  
#             ($ISTAGE = integer)
#
#  SALT2train_pipeline.pl <inputFile>  RESTART $ISTAGE  
#                ($ISTAGE = integer)
#
#  SALT2train_pipeline.pl <inputFile>  VBOSE
#
#  SALT2train_pipeline.pl <inputFile>  KILL    ! kill jobs on all nodes
#
#  SALT2train_pipeline.pl <inputFile>  QUICKTEST ! replace train with copy
#  SALT2train_pipeline.pl <inputFile>  RUN       ! run QUICKTEST
#
#
#   History
# ~~~~~~~~~~~~~~~~~~~~~~ 
# Aug 15, 2011:  insert "FUDGEALL_ITER1_MAXFRAC = 0.01" into &FITINP
#                (along with MXLC_PLOT=5)
#
# Aug 26, 2011:  SIMDATA GENVERSION is appended with Nsec-since-midnight
#                so that multiple test can run in parallel using same
#                SIMGEN and SIMFIT files for SIMDATA.
#
# Aug 31, 2011   add quick link to find training dirs:
#                  ln -s STAGE03_TrainRun/workSpace/  AAA_PCAFIT_TOPDIR
#
# Sep 01, 2011   Add new $ISTAGE_CLEAN to purge/remove files from 
#                the training. Cleaning can be turned off with
#                "CLEANFLAG: 0"
#
# Sep 9, 2011: read optional key BATCH_INFO.
#              Must have either the BATCH key or NODELIST key (for SSH).
#              The BATCH_INFO keys are written to the train_run input file
#              and to the split_and_fit input file.
#
# Sep 23, 2011: new key PAW_COMMAND: paw [or pawX11] to allow arbitrary
#               alias to call paw.
#               Allow PAWPLOTMACRO or PAW_PLOTMACRO
#
# Sep 27, 2011: add optional archive stage. See new key ARCHIVE_OPT.
#
# Oct 04, 2011: use new -outprefix option for combine_fitres.exe
#               so that parallel jobs don't clobber each other.
#
# Oct 18, 2011: in SALT2mu stage, use h2root to convert SALT2-ntuple
#               into root file.
#
# Oct 20, 2011: abort if any output dir cannot be created.
#
# Nov 4, 2011: lots of changes to handle DATA or SIM
#              See $DATAFLAG and $SIMFLAG.
#              Note that arrays with _DATA or DATA_ indicate
#              that they work for both data and sim.
#              DATAFIT/DATATRAIN is specific to DATA;
#              SIMDATA/SIMTRAIN is specific to the simulation
#
#     Add QUICKTEST option to command-line =>
#     replace training with symbolic link to existing SALT2 model.
#     This options tests the scripting infrastructure.
#
# Nov 10,2011: re-named to SALT2train_pipeline.pl
#
# Nov 17, 2011: call check_split_version() and  abort if
#               any data-version has already been split.
#
# Dec 10, 2011: require GENSPEC keys after each GENPHOT key to allow
#               different spec-sim for each phot sim. Can use [] to repeat.
#
# Dec 29, 2011: fix a few bugs found by Rahul
#   * $VERSION.nml -> FIT_$VERSION.nml to avoid moving original 
#     nml file when prefix is the same as the version name.
#
#   * In STAGE05_DATAFIT, glue GROUP+ VERSION to make directory name
#      (intead of just version as directory name)
#
# Jan 04, 2012
#     Abort if GEN_SNDATA_SIM key is used or if
#     bit2 of FORMAT_MASK is not set. (see sub check_simFile)
#
# Jan 18, 2012: 
#   * in RUN_pipeline set g+rw priv on each $STAGE_PATH.
#   * new sub set_STAGE_PATH to always set $STAGE_PATH and $STAGE_SDIR
#     (previously they were set only for init mode)
#
# Jan 23, 2012: new option "OPT_SPLIT: 0" to NOT split text files;
#               the data reading is slower, but allows multiple
#               pipeline jobs to runs simultaneously.
#
#
# Feb 03, 2012: add "RESTART $ISTAGE" option
# 
# Feb 04, 2012: allow for path in DATATRAIN_LISTFILE 
# Feb 05, 2012: abort if duplicate pipeline job is deteceted (abort_on_duplicate)
# Mar 07, 2012: allow master file to input file key SALT2_MAGERRFLOOR (JLM)
#
# Jun 07, 2012 - JLM
#  Allow master file to input file key SALT_TEMPLATE0, SALT2_TEMPLATE1 (JLM)
#  Lots of changes to add cosmology fitting stage to pipeline.
#  Cosmology fit is off by default. Input file keywords are 
#  COSMO_PROGRAM, COSMO_INFILE, and COSMO_GENOPT.
#  The same cosmology fit is used for all data groups and trains. 
#
#
# Jun 28, 2012:
#  in make_OUTDIRs(), call checkDirExist(..) to make sure that
#  the full path is specified for OUTPUT_TOPDIR. Aborts on the
#  common mistake of setting OUTPUT_TOPDIR to a local subdir path
#  such as OUT.
#
# Aug 28, 2012 - JLM
#  Allow SIMDATA option SIMDATA_EXCLUDE_FITOPT:
#  If set for a SIMDATA_GROUP, no trained model fits
#  of that group will not be performed, only default
#  model fits as determined by SIMDATA nml INFILE. 
#
# Aug 29, 2012 - JLM
#  Also - if SIMDATA_EXCLUDE_FITOPT: key is set,
#  the default use of GLOBAL_GENMODEL will NOT APPLY
#  to that groups' simulations. If a different
#  genmodel is desired, it must be applied in the
#  groups' sim input files or via the groups' 
#  SIMDATA_GENOPT: keyword. 
#
# Sep 04, 2012 - JLM
#  Allow user to skip SPECSIM for a given GENVERSION
#  by EXPLICITLY setting "SIMTRAIN_GENSPEC_INFILE: NONE"   
#
# Sep 07, 2012 - JLM
#  Added extra stage (3) to pipeline to allow real training
#   sets to be combined with simulated data sets.
#   New keywords are:
#   SIMTRAIN_INCDATA_INFILE: <full path to trainlist file>
#   SIMTRAIN_INCDATA_CHANGE: GENVERSION <>  DATAPATH <> 
#   See header of SALT2train_incdata.pl for more info.
# 
# Sep 26, 2012 - JLM
#  Allow master file to input file key SALT2_DIR: to be passed
#  through to SALT2train_run.pl. 
#
# Oct 23, 2012 - JLM
#  Allow master file to input file key SALTPATH: to be passed 
#  through to SALT2train_run.pl. 
#
#
# Feb 8, 2013 RK - MXLC_PLOT -> &SNLCINP instead of &FITINP (v10_24b)
# 
# Mar 25, 2013 JLM - Add optional bias correction stage. 
#                   Requires updated split_and_fit.pl (v10_24i),
#                   and helper scripts SALT2biascorr.pl and
#                   SALT2extend_colordisp.pl. Input keywords are 
#                   described in the parse_DMCORR subroutine header.
#
# Apr 26, 2013: new option LASTSTAGE <ISTAGE> to process up to
#               LASTSTAGE. Hopefully the RESTART option will work
#               to finish, but needs to be checked.
#
# May 03, 2013: added 5th jobid "TRAINLIST" to init_Pipeline SIMFLAG section
#               this allows non-default TRAIN_LIST: keyword to be passed
#               to SALT2train_run.pl input file. If using default TRAIN_LIST,
#               must set SIMTRAIN_TRAINLIST_INFILE: NONE for each SIMTRAIN_MODEL
#               group. 
#
# May 23, 2013: Altered BATCH handling to use two BATCH_INFO input keywords, 
#               BATCH_INFO_TRAIN: and BATCH_INFO_FIT: . This change allows
#               separate memory and node usages to be set for training and fitting. 
#               The old BATCH_INFO: keyword is still operational for 
#               backwords compatibility, but may not be used in combination with
#               either or both of the new keywords.  (JLM)
#
# May 23, 2013: Fixed bug in $ISTAGE_STOP initialization by adding new
#               subroutine set_stage_startstop which will run after the
#               number of stages has been determined by initPipeline(). 
#
# May 25, 2013: Altered subroutine make_DMFIXSIM_scripts to be 
#               QUICKTEST compatible. (JLM)
#       
# May 28, 2013: Altered subroutine make_DMFIXSIM_scripts to use snlc_sim
#               SALT2mu_FILE keyword, direct input of model color dispersion
#               (snpca-2.3.20b and later), and sim_snMix.pl. (JLM)
# 
# Jun 01, 2013: Added stages to check BINARY files prior to batch simulation
#               job submissions. 
#
# Jun 03, 2013: 
#   Altered sim_SNmix call to use number of nodes/cores as simgenversions,
#   added check_parallel_batch sub to ensure more nodes than TRAIN, DATSIM
#   jobs, and cleaned up debug statements. 
#
# Oct 10, 2013 JLM -
#   Minor change to facilitate private versions of bias correction, 
#   cosmology codes. Can use [] to repeat DMCORR_RUNCODEOPT: keywords. 
# 
# Oct 16, 2013: 
#   Added global variable GLOBAL_FIT_RUNOPT to allow command line additions
#   to split_and_fit.pl calls. Variable is currently hardcoded to " NOPROMPT " 
#   in parse_master_infile. 
#    Also added new keywords SIMGEN_RUNOPT: and DMCORR_RUNOPT: to allow 
#    command line options to sim_SNmix.pl jobs.   (JLM)
#
# Mar 14 2014 RK -  fix bug inserting  FUDGEALL_ITER1_MAXFRAC into &FITINP
#
# ---------------------------------------------------------

use IO::Handle;
use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use sntools ;
use strict ;

# declare subs
sub abort_on_duplicate ;
sub set_stage_startstop ;
sub initPipeline ;
sub setDataType ;
sub setDMCORRType ;
sub set_STAGE_PATH ;
sub parse_args ;
sub parse_public_config ;
sub parse_master_infile ;
sub parse_SIMTRAIN ;
sub parse_SIMDATA ;
sub parse_DATATRAIN ;
sub parse_DATAFIT   ;
sub parse_COSMO ;
sub parse_DMCORR ;
sub parse_PAWPLOTMACRO ;
sub check_SSH_or_BATCH ;
sub check_split_version ;
sub check_parallel_batch(@);
sub RUN_pipeline ;
sub make_OUTDIRS ;
sub make_TrainGen_inputs ;
sub make_TRAINGEN_scripts(@); 
sub make_TrainBin_inputs ;
sub make_runTrain_script ;
sub make_fitnmlFiles ;
sub copy_inputFiles ;
sub copy_inputFiles_TRAINGEN ;
sub copy_inputFiles_TRANSLATE ;
sub copy_inputFiles_INCDATA ;
sub copy_inputFiles_TRAINRUN ;
sub copy_inputFiles_SIMGEN(@) ;
sub copy_inputFiles_FIT ;
sub copy_inputFiles_COSMO(@);
sub copy_inputFiles_DMFIXSIM ;
sub copy_inputFiles_DMFIXFIT ;
sub copy_inputFiles_DMFIXRUN ;
sub make_RUNscripts ;
sub make_SIMFIT_scripts ;
sub make_SIMGEN_BINARY_scripts ;
sub make_SIMGEN_scripts ;
sub make_PLOT_scripts ;
sub make_SALT2mu_scripts ;
sub make_COSMO_scripts ;
sub make_CLEAN_script ;
sub make_ARCHIVE_script ;
sub update_TRAIN_INFILE(@);
sub update_TRAIN_CHANGE(@); 
sub scriptFileName(@) ;
sub getGENVERSION(@); 
sub prefixSALT2mu(@);
sub get_SALT2mu_commands(@);
sub get_COSMO_commands(@);
sub get_CATLISTinfo(@);
sub get_FITOPTKEYline(@);
sub wait_for_Done ;
sub check_Done ;
sub I2STAGE_extract ;

# declare globals

my $OPTABORT = 0 ;
my $OPTWARN  = 1 ;
my $OPTQUIET = 2 ;

# start with user-input
my ($MASTER_INFILE, @MASTER_INFILE_CONTENTS, $NLINE_MASTER_INFILE) ;
my ($PUBLIC_CONFIG_DIR, $CONFIG_FILE, $WAIT_SCRIPT);
my ($ISTAGE_STOP, $ISTAGE_START, $INPUT_START, $INPUT_STOP );
my ($SALT2_DIR, $SNANA_MODELPATH, @FILELIST_PURGE );
my ($OPT_PIPELINE, $VBOSEFLAG, @USER_COMMENTS );
my (@NODELIST, $NNODE, $OUTPUT_TOPDIR, $GLOBAL_GENMODEL) ;
my ($GLOBAL_FIT_RUNOPT);
my ($INFILE_SALT2mu, @NODELIST_SALT2mu );
my ($N_PAWPLOTMACRO, @PAWPLOTMACRO_NAME, @PAWPLOTMACRO_ARG, $PAW_COMMAND );
# xxx  new stuff for dual Batch files
my ($BATCH_FLAG, $BATCH_INFO);
my ($BATCH_INFO_TRAIN, $BATCH_TEMPLATE_TRAIN, $BATCH_COMMAND_TRAIN);
my ($BATCH_INFO_FIT, $BATCH_TEMPLATE_FIT, $BATCH_COMMAND_FIT);
my ($BATCH_NODE_TRAIN, $BATCH_NODE_FIT);
my ($OPT_CLEAN, $OPT_ARCHIVE, $OPT_SPLIT, $OPT_COSMO );

# internal globals
my ($JOBID_GENPHOT, $JOBID_GENSPEC, $JOBID_PCAFIT, $JOBID_INCDATA)  ;
my ($JOBID_TRAINLIST);
my ($NJOBID, @JOBID_NAME, $INCLUDE_KEY ) ;
my ($SNANA_DIR, $SNDATA_ROOT, $USER, $SHELL );

my ($SIMFLAG, $DATAFLAG, $DMCORRFLAG);

my ($NSTAGE, $ISTAGE_TRAINGEN, $ISTAGE_TRANSLATE, $ISTAGE_TRAINRUN) ;
my ($ISTAGE_TRAINBINARY, $ISTAGE_DATABINARY);
my ($ISTAGE_SIMGEN, $ISTAGE_SIMFIT, $ISTAGE_SALT2mu, $ISTAGE_COSMO) ;
my ($ISTAGE_DMFIXSIM, $ISTAGE_DMFIXFIT, $ISTAGE_DMFIXPREP);
my ($ISTAGE_DMFIXRUN, $ISTAGE_DMFIXCOSMO);
my ($ISTAGE_PLOTS, $ISTAGE_CLEAN, $ISTAGE_ARCHIVE, $ISTAGE_INCDATA );
my (@STAGE_NAME, @STAGE_SLEEPCHECK, @STAGE_SDIR, @STAGE_PATH);
my (@SCRIPT_PREFIX, @PROGRAM_NAME,  $PROGRAM_COMBINE_FITRES );
my ($INFILE_TRAINGEN, $INFILE_TRAINRUN, $INFILE_INCDATA) ;
my ($INFILE_TRAINBINARY);

my ($NGROUP_DATA, $NSIMGEN, @NVERSION_GROUP ) ;
my (@DATA_GROUPNAME, @SIMDATA_INFILE_GEN, @DATA_INFILE_FIT) ;
my (@DATA_GROUPID, @SIMDATA_GENVERSION); 
my (@DATA_FITNML, @SIMDATA_GENOPT, @SIMDATA_RUNOPT); 
my (@SIMDATA_EXCLUDE_FITOPT) ;
my ($COSMO_GENOPT, $COSMO_INFILE, @CPLIST_FITRES) ;

my ($NGROUP_DMCORR, $NSIMGENDMCORR, @DMCORR_IDEAL) ;
my (@DMCORR_SIMCODE, @DMCORR_RUNCODE ) ;
my (@DMCORR_GROUPNAME, @DMCORR_TYPE,  @DMCORR_INFILE_GEN ) ;
my (@DMCORR_INFILE_FIT, @DMCORR_GENOPT, @DMCORR_EXCLUDE_FITOPT ) ;
my (@DMCORR_GROUPID, @DMCORR_GENVERSION, @DMCORR_SIMTYPE ) ;
my (@DMCORR_FITNML, @NVERSION_GROUP_DMCORR, @DMCORR_RUNOPT  ) ;
my (@DMCORR_MEANC, @DMCORR_MEANX, @DMCORR_SIGLC, @DMCORR_SIGRC ) ;
my (@DMCORR_SIGLX, @DMCORR_SIGRX, @DMCORR_RUNCODEOPT , @CPLIST_CORR) ;
my (@DMCORR_SIMCODEOPT, @DMCORR_SIMALPHA, @DMCORR_SIMBETA) ;

my ($NMODEL_TRAIN, $NFIT_ITER_PCAFIT ) ;
my (@TRAIN_MODELNAME, @TRAIN_CHARINDEX) ;
my (@NINFILE_TRAIN, @NCHANGE_TRAIN, @USE_CHANGE) ;
my (@TRAIN_INFILE, @TRAIN_CHANGE) ;
my (@SIMTRAIN_REPEATFLAG, @SIMTRAIN_GENVERSION) ;
my (@SIMTRAIN_INCDATA_GENVERSION);

my (@DATATRAIN_PATH, @DATATRAIN_LISTFILE) ;
my (@DATAFIT_VERSION );

my (@INFILE_FITNML_LIST, @FITNML_LIST, @VERSION_LIST ) ;
my (@CATLIST_FITRES, @NCATLIST_FITRES );

my (@MSGERR) ;


# hard-wire list of optional keys to copy from MASTER input file
# into the training input file used by SALT2train_run.pl
my @GLOBAL_TRAINRUN_KEYNAMES  
    = ( "TRAINPLOTS_REFSURFACE_FILE:" , 
	"TRAINPLOTS_REFCOLORLAW_FILE:" ,
	"TRAINPLOTS_PHASELIST:" ,
	"SALT2STAT_TESTLIST_FILE:", 
	"SALT2_MAGERRFLOOR:",
	"SALT_TEMPLATE0:",
	"SALT_TEMPLATE1:",
        "SALT2_DIR:",
	"SALTPATH:"
	);
my @GLOBAL_TRAINRUN_KEYVALS ;   # values read from MASTER input file


# --------- BEGIN MAIN -----------

$SHELL=sntools::shellName();

&parse_args();

# determine if training on data or simulations
&setDataType();

# determine if DM bias corrections will be performed
&setDMCORRType();

&initPipeline();

&set_stage_startstop();

&abort_on_duplicate();

# parse master input file
&parse_master_infile();

# parse public training-config file.
&parse_public_config();

&set_STAGE_PATH();


# check command-line options
if ( $OPT_PIPELINE eq "KILL"      ) { sntools::Killjobs(@NODELIST); exit ; }
if ( $OPT_PIPELINE eq "RUN"       ) { &RUN_pipeline();  }
if ( $OPT_PIPELINE eq "RESTART"   ) { &RUN_pipeline();  }
if ( $OPT_PIPELINE eq "LASTSTAGE" ) { &RUN_pipeline();  }

if ( $OPT_PIPELINE eq "QUICKTEST" ) { 
    $NFIT_ITER_PCAFIT = "QUICKTEST" ; $OPT_CLEAN = 0 ;
}

# create output directory for each  stage
&make_OUTDIRS();

# prepare input-file to generate simulated training sample
if ( $SIMFLAG )  { 
    &make_TrainBin_inputs(); 
    &make_TrainGen_inputs(); 
}

# prepare script for SALT2train_run.pl
&make_runTrain_script();  

# create SALT2 fit-namelist file for each DATA sample
print "\n";
my $igroup;
for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) {    
    &make_fitnmlFiles($igroup, $ISTAGE_SIMFIT);
}
# if DMCORRFLAG create SALT2 fit-namelist file for each DMFIX DATA sample
if ( $DMCORRFLAG ) {
    print "\n";
    for ( $igroup = 1; $igroup <= $NGROUP_DMCORR; $igroup++ ) {
	&make_fitnmlFiles($igroup, $ISTAGE_DMFIXFIT);
    }
}

# copy input files for all stages
&copy_inputFiles();

# if ( $DATAFLAG ) { die "\n xxxxx debug abort xxxxx \n" ; }

# prepare RUN#_XXX scripts for each stage
&make_RUNscripts();


# ----------------------------------------------
#  END OF MAIN
# --------------------------------------------


sub initPipeline {

    my ($i,$j,$k);

    if ( $SIMFLAG ) {
	$JOBID_GENPHOT   = 1 ;
	$JOBID_GENSPEC   = 2 ;
	$JOBID_INCDATA   = 3 ;
	$JOBID_PCAFIT    = 4 ;
	$JOBID_TRAINLIST = 5 ;
	$NJOBID          = 5 ;
	@JOBID_NAME = ( "" , "GENPHOT", "GENSPEC" , "INCDATA", "PCAFIT", "TRAINLIST" );
    }
    else {
	$JOBID_PCAFIT = 1 ;
	$NJOBID       = 1 ;
	@JOBID_NAME   = ( "" , "PCAFIT" );
    }

    $SNDATA_ROOT   = $ENV{'SNDATA_ROOT'};
    $SNANA_DIR     = $ENV{'SNANA_DIR'};
    $USER          = `whoami` ;    $USER  =~ s/\s+$// ;

    $OPT_CLEAN    = 1 ;  # flag to purge/remove files at end
    $OPT_ARCHIVE  = 1 ;  # default is to make archive
    $OPT_SPLIT    = 1 ;  # default is to split data (text only)
    $OPT_COSMO    = 0 ;  # default is to skip cosmo fit

    $INCLUDE_KEY  = "INPUT_FILE_INCLUDE:" ;

    $PUBLIC_CONFIG_DIR = "$SNDATA_ROOT/models/SALT2" ;
    $CONFIG_FILE       = "SALT2training.config" ;
    $WAIT_SCRIPT       = "wait_for_files.pl";  # name of utility script

    $N_PAWPLOTMACRO = 0;
    $PAW_COMMAND    = "paw";

    $NNODE      = 0;
    $BATCH_FLAG = 0;


    # SIMDATA // DMFIX parameters
    # k: 1==SIMDATA, 2==DMFIX
    for ( $k=0; $k < 3; $k++ ) {
	for ( $i=0; $i < 200; $i++ ) {
	    $NVERSION_GROUP[$i] = 0 ;
	    for ( $j=0; $j < 20; $j++ )  { $NCATLIST_FITRES[$k][$i][$j] = 0 ; }
	}
    }

    # init pcafit @nfit key to use error snake on 2nd iteration
    $NFIT_ITER_PCAFIT = 2 ; 

    # define stages and name for each stage

    $ISTAGE_TRAINBINARY = 1 ; # check/make BINARY files for SIMSED sims 
    $ISTAGE_TRAINGEN    = 2 ; # simulate training sample(s) in SNANA format
    $ISTAGE_TRANSLATE   = 3 ; # translate into SALT2 format
    $ISTAGE_INCDATA     = 4 ; # make SIM/DATA hybrid data dirs
    $ISTAGE_TRAINRUN    = 5 ; # run the training(s)
    $ISTAGE_DATABINARY  = 6 ; # check/make BINARY files for SIMSED data sims
    $ISTAGE_SIMGEN      = 7 ; # simulate independent test sample(s)
    $ISTAGE_SIMFIT      = 8 ; # fit independent samples with trained model(s)
    $ISTAGE_SALT2mu     = 9 ; # run SALT2mu on fits to get MU-resids
    $ISTAGE_COSMO       = 10 ; # run wfit on MU-resids to get cosmology
    $ISTAGE_DMFIXSIM    = 11 ; # generate train-model-based sims for DM corr
    $ISTAGE_DMFIXFIT    = 12 ; # fit train-model-based sims for DM corr
    $ISTAGE_DMFIXPREP   = 13 ; # prepare files for DM bias correction
    $ISTAGE_DMFIXRUN    = 14 ; # apply DM bias correction
    $ISTAGE_DMFIXCOSMO  = 15 ; # rerun cosmology with corrected DM
    $ISTAGE_PLOTS       = 16 ; # plot MU-resids from SALT2mu
    $ISTAGE_CLEAN       = 17 ; # remove/gzip files
    $ISTAGE_ARCHIVE     = 18 ; # make tarball
    $NSTAGE       =  18;

    $STAGE_NAME[$ISTAGE_TRAINBINARY] = "TrainBin";
    $STAGE_NAME[$ISTAGE_TRAINGEN]    = "TrainGen" ;
    $STAGE_NAME[$ISTAGE_TRANSLATE]   = "Translate" ;
    $STAGE_NAME[$ISTAGE_INCDATA]     = "INCDATA"  ;
    $STAGE_NAME[$ISTAGE_TRAINRUN]    = "TrainRun" ;
    $STAGE_NAME[$ISTAGE_DATABINARY]  = "DataBin";
    $STAGE_NAME[$ISTAGE_SIMGEN]      = "SIMGEN"   ;
    $STAGE_NAME[$ISTAGE_SIMFIT]      = "SIMFIT"   ;
    $STAGE_NAME[$ISTAGE_SALT2mu]     = "SALT2mu"  ;
    $STAGE_NAME[$ISTAGE_COSMO]       = "COSMO"  ;
    $STAGE_NAME[$ISTAGE_DMFIXSIM]    = "DMFIXSIM" ;
    $STAGE_NAME[$ISTAGE_DMFIXFIT]    = "DMFIXFIT" ;
    $STAGE_NAME[$ISTAGE_DMFIXPREP]   = "DMFIXPREP" ;
    $STAGE_NAME[$ISTAGE_DMFIXRUN]    = "DMFIXRUN" ;
    $STAGE_NAME[$ISTAGE_DMFIXCOSMO]  = "DMFIXCOSMO" ;
    $STAGE_NAME[$ISTAGE_PLOTS]       = "PLOTS"    ;
    $STAGE_NAME[$ISTAGE_CLEAN]       = "CLEAN"    ;
    $STAGE_NAME[$ISTAGE_ARCHIVE]     = "ARCHIVE"  ;

    # sleep time between checking when each stage is done.
    $STAGE_SLEEPCHECK[$ISTAGE_TRAINBINARY] = 120 ;
    $STAGE_SLEEPCHECK[$ISTAGE_TRAINGEN]    =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_TRANSLATE]   =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_INCDATA]     =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_TRAINRUN]    = 300 ;
    $STAGE_SLEEPCHECK[$ISTAGE_DATABINARY]  = 120 ;
    $STAGE_SLEEPCHECK[$ISTAGE_SIMGEN]      =  20 ;
    $STAGE_SLEEPCHECK[$ISTAGE_SIMFIT]      = 300 ;
    $STAGE_SLEEPCHECK[$ISTAGE_SALT2mu]     =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_COSMO]       =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_DMFIXSIM]    =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_DMFIXFIT]    = 200 ;
    $STAGE_SLEEPCHECK[$ISTAGE_DMFIXPREP]   =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_DMFIXRUN]    =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_DMFIXCOSMO]  =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_PLOTS]       =  10 ;
    $STAGE_SLEEPCHECK[$ISTAGE_CLEAN]       =  60 ;
    $STAGE_SLEEPCHECK[$ISTAGE_ARCHIVE]     =  10 ;

    # define prefix names of stages that use another script
    $SCRIPT_PREFIX[$ISTAGE_TRAINBINARY]  = "TRAIN_SIM";
    $SCRIPT_PREFIX[$ISTAGE_TRAINGEN]     = "TRAIN_SIM";
    $SCRIPT_PREFIX[$ISTAGE_INCDATA]      = "SALT2train_incdata";
    $SCRIPT_PREFIX[$ISTAGE_TRANSLATE]    = "SALT2train_translate";
    $SCRIPT_PREFIX[$ISTAGE_TRAINRUN]     = "SALT2train_run";
    $SCRIPT_PREFIX[$ISTAGE_SIMFIT]       = "split_and_fit" ;
    $SCRIPT_PREFIX[$ISTAGE_DMFIXSIM]     = "SALT2train_extcolordisp";
    $SCRIPT_PREFIX[$ISTAGE_DMFIXFIT]     = "split_and_fit" ;
    $SCRIPT_PREFIX[$ISTAGE_DMFIXRUN]     = "SALT2train_biascorr" ;

    # define program names used vs. stage
    $PROGRAM_NAME[$ISTAGE_DATABINARY] = "snlc_sim.exe" ;
    $PROGRAM_NAME[$ISTAGE_SIMGEN]     = "sim_SNmix.pl" ;
    $PROGRAM_NAME[$ISTAGE_SALT2mu]    = "SALT2mu.exe" ;
    $PROGRAM_COMBINE_FITRES           = "combine_fitres.exe" ;
    $PROGRAM_NAME[$ISTAGE_COSMO]      = "wfit.exe" ;
    $PROGRAM_NAME[$ISTAGE_DMFIXSIM]   = "sim_SNmix.pl" ;
    $PROGRAM_NAME[$ISTAGE_DMFIXPREP]  = "SALT2mu.exe" ;
    $PROGRAM_NAME[$ISTAGE_DMFIXCOSMO] = "wfit.exe" ; 


    if ( $DATAFLAG ) {
	$STAGE_NAME[$ISTAGE_TRAINBINARY] = "SKIP" ;
	$STAGE_NAME[$ISTAGE_TRAINGEN]    = "SKIP" ;
	$STAGE_NAME[$ISTAGE_TRANSLATE]   = "SKIP" ;
	$STAGE_NAME[$ISTAGE_INCDATA]     = "SKIP" ;
	$STAGE_NAME[$ISTAGE_DATABINARY]  = "SKIP"   ;
	$STAGE_NAME[$ISTAGE_SIMGEN]      = "SKIP"   ;
	$STAGE_NAME[$ISTAGE_SIMFIT]      = "DATAFIT"  ;
	$STAGE_NAME[$ISTAGE_COSMO]       = "SKIP"   ;
	$STAGE_NAME[$ISTAGE_PLOTS]       = "SKIP"   ;
	$STAGE_NAME[$ISTAGE_ARCHIVE]     = "SKIP"   ;

	$STAGE_SLEEPCHECK[$ISTAGE_TRAINBINARY] =  1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_TRAINGEN]    =  1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_TRANSLATE]   =  1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_INCDATA]     =  1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_DATABINARY]  =  1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_SIMGEN]      =  1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_SIMFIT]      =  1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_COSMO]       =  1 ;
    }

    unless ( $DMCORRFLAG ) {
	$STAGE_NAME[$ISTAGE_DMFIXSIM]   = "SKIP" ;
	$STAGE_NAME[$ISTAGE_DMFIXFIT]   = "SKIP" ;
	$STAGE_NAME[$ISTAGE_DMFIXPREP]  = "SKIP" ;
	$STAGE_NAME[$ISTAGE_DMFIXRUN]   = "SKIP" ;
	$STAGE_NAME[$ISTAGE_DMFIXCOSMO] = "SKIP" ;

	$STAGE_SLEEPCHECK[$ISTAGE_DMFIXSIM]   = 1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_DMFIXFIT]   = 1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_DMFIXPREP]  = 1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_DMFIXRUN]   = 1 ;
	$STAGE_SLEEPCHECK[$ISTAGE_DMFIXCOSMO] = 1 ;
    }

}  # end of  initPipeline


# =======================================
sub parse_public_config {

    my ($key, @tmp, $FILE);


    if ( -e $CONFIG_FILE ) {
	$FILE  =  "$CONFIG_FILE" ;
    }
    else {
	$FILE  =  "$PUBLIC_CONFIG_DIR/$CONFIG_FILE" ;
    }

    print " Parse public config file: \n\t $FILE \n";

    $key   =  "SALT2_DIR:" ;
    @tmp   =  sntools::parse_line($FILE, 1, $key, $OPTABORT );
    $SALT2_DIR = $tmp[0] ;

    $key   =  "SNANA_MODELPATH:" ;
    @tmp   =  sntools::parse_line($FILE, 1, $key, $OPTABORT );
    $SNANA_MODELPATH = $tmp[0] ;

    $key   =  "FILELIST_PURGE:" ;
    @tmp   =  sntools::parse_line($FILE, 99, $key, $OPTABORT );
    @FILELIST_PURGE = @tmp ;


} # end of  parse_public_config

# =======================================
sub setDataType {

    # set either $SIMFLAG or $DATAFLAG to TRUE.
    # Check for keys SIMTRAIN_MODEL and  DATATRAIN_MODEL
    
    my ($key, @tmp );

    $key = "SIMTRAIN_MODEL:" ; 
    @tmp   =  sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET );
    $SIMFLAG = scalar(@tmp);

    $key = "DATATRAIN_MODEL:" ; 
    @tmp   =  sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET );
    $DATAFLAG = scalar(@tmp);


    # idiot error checking
    if ( $DATAFLAG == 0 && $SIMFLAG == 0 ) {
	$MSGERR[0] = "Neither SIMTRAIN_MODEL or DATATRAIN_MODEL is given.";
	$MSGERR[1] = "Must specify either SIM or DATA model(s) to train." ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ( $DATAFLAG && $SIMFLAG ) {
	$MSGERR[0] = "Both SIMTRAIN_MODEL and DATATRAIN_MODEL are given.";
	$MSGERR[1] = "Must specify either SIM or DATA model(s) to train." ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ( $DATAFLAG ) { print " ***** Will train on DATA ***** \n";   }
    if ( $SIMFLAG  ) { print " ***** Will train on SIMULATIONS ***** \n"; }
    
    print "\n" ;

} # end of setDataType

# =======================================
sub setDMCORRType {

    # set $DMCORRFLAG to true
    # Check for keys DMCORR_GROUPNAME:
    # print "Checking for DMCORR keywords \n";
    
    my ($key, @tmp );

    $key = "DMCORR_GROUPNAME:" ; 
    @tmp   =  sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET );
    $DMCORRFLAG = scalar(@tmp);

    if ( $DMCORRFLAG ) { print " ***** Will apply DM bias corrections ***** \n";   }
    
    print "\n" ;

} # end of setDMCORRType


# ============================================
sub parse_args {

    # parse command-line argument(s)

    my ($NARG);

    $NARG = scalar(@ARGV);
    $MASTER_INFILE = $ARGV[0];
    $VBOSEFLAG    = 0 ;
    $INPUT_START = -1;
    $INPUT_STOP = -1;

    if ( ! ( -e $MASTER_INFILE ) ) {	       
	$MSGERR[0] = "Master input file '${MASTER_INFILE}'" ;
	$MSGERR[1] = "does not exist.";
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ( $NARG > 1 ) 
    { $OPT_PIPELINE = $ARGV[1] ; }
    else
    { $OPT_PIPELINE = " " ;  return ; }


    # check OPT_PIPELINE

    if ( $OPT_PIPELINE eq "VBOSE"   ) 
    { $VBOSEFLAG = 1 ; }

    elsif ( $OPT_PIPELINE eq "RESTART" ) 
    { $INPUT_START = $ARGV[2] ;  }

    elsif ( $OPT_PIPELINE eq "LASTSTAGE" ) 
    { $INPUT_STOP = $ARGV[2] ;  }

    elsif ( $OPT_PIPELINE eq "QUICKTEST" ) 
    { }
    elsif ( $OPT_PIPELINE eq "RUN" ) 
    { }
    elsif ( $OPT_PIPELINE eq "KILL" ) 
    { }
    else {
	$MSGERR[0] = "Invalid pipeline option: '${OPT_PIPELINE}' ";
	$MSGERR[1] = "Check 2nd argument.";
	sntools::FATAL_ERROR(@MSGERR);	
    }

} # end of parse_args

# ===========================================
sub set_stage_startstop {

  # Created May 23, 2013
  # Run after initPipeline to determine
  # appropriate ISTAGE_START and
  # ISTAGE_STOP values (JLM)

    if ( $INPUT_START < 0 ) {
	$ISTAGE_START = 1;
    } else {
	$ISTAGE_START = $INPUT_START;
    }

    if ( $INPUT_STOP < 0 ) {
	$ISTAGE_STOP = $NSTAGE; 
    } else {
	$ISTAGE_STOP = $INPUT_STOP;
    }

} # end of abort_on_duplicate

# ===========================================
sub abort_on_duplicate {

    # Created Feb 5, 2012
    # Abort if this same job is already running.
    # Avoids interfering jobs if the pipeline hangs
    # and the user re-starts.
    
    my (@ps, $g0, $g1, $g2, $Nps);
    my $Nps_ABORT = 3 ; # ps + current job + duplicate

    if ( $OPT_PIPELINE == "KILL" ) { return ; }

    # define grep-args
    $g0 = $USER ;
    $g1 = "SALT2train_pipeline" ;
    $g2 = $MASTER_INFILE ;
	
    @ps = qx(ps -eaf | grep $g0 | grep $g1 | grep $g2 );

    $Nps = scalar(@ps);
    if  ( $Nps >= $Nps_ABORT  ) {
	$MSGERR[0] = "Detected duplicate ${USER} jobs for";
	$MSGERR[1] = "   $g1 $g2 \n" ;
	$MSGERR[2] = "Either kill duplicate jobs or";
	$MSGERR[3] = "change name of master input file.";
	sntools::FATAL_ERROR(@MSGERR);
    }

#    die "\n xxx Nps = $Nps xxxx\n";

} # &abort_on_duplicate


# =======================================
sub parse_master_infile {

    my ($INFILE, $key, @tmp, $name );

    $INFILE =  $MASTER_INFILE ;

    $key = "COMMENT:" ;
    @USER_COMMENTS = sntools::parse_line($MASTER_INFILE, 99, $key, $OPTQUIET);

    $key = "NODELIST:" ;
    @NODELIST = sntools::parse_line($MASTER_INFILE, 99, $key, $OPTQUIET ) ;
    $NNODE = scalar(@NODELIST);

    $key = "BATCH_INFO:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 3, $key, $OPTQUIET ) ;
    if ( scalar(@tmp) > 0 ) {	
	$BATCH_FLAG += 1;
	$BATCH_INFO = "$tmp[0]" ;
    }

    $key = "BATCH_INFO_TRAIN:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 3, $key, $OPTQUIET ) ;
    if ( scalar(@tmp) > 0 ) {	
	$BATCH_FLAG += 1;
	$BATCH_INFO_TRAIN = "$tmp[0]" ;
    }

    $key = "BATCH_INFO_FIT:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 3, $key, $OPTQUIET ) ;
    if ( scalar(@tmp) > 0 ) {	
	$BATCH_FLAG += 1;
	$BATCH_INFO_FIT = "$tmp[0]" ;
    }

    # make sure we have NODELIST(SSH) or BATCH, but  not both.
    &check_SSH_or_BATCH();

    $key = "OUTPUT_TOPDIR:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 1, $key, $OPTABORT ) ;
    $OUTPUT_TOPDIR = $tmp[0];
    print "  OUTPUT_TOPDIR: $OUTPUT_TOPDIR \n" ;

    if ( $SIMFLAG ) {
	$key = "GLOBAL_GENMODEL:" ;
	@tmp = sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET ) ;
	$GLOBAL_GENMODEL = $tmp[0];
	print "  GLOBAL_GENMODEL: $GLOBAL_GENMODEL \n" ;
    }

    $key = "OPT_CLEAN:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET ) ;
    if ( scalar(@tmp) > 0 ) { $OPT_CLEAN = $tmp[0]; }
    print "  OPT_CLEAN:  $OPT_CLEAN \n" ;

    $key = "OPT_ARCHIVE:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET ) ;
    if ( scalar(@tmp) > 0 ) { $OPT_ARCHIVE = $tmp[0]; }
    print "  OPT_ARCHIVE:  $OPT_ARCHIVE \n" ;

    $key = "OPT_SPLIT:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET ) ;
    if ( scalar(@tmp) > 0 ) { $OPT_SPLIT = $tmp[0]; }
    print "  OPT_SPLIT:  $OPT_SPLIT \n" ;
    

    # -------
    # scoop up the entire input file for linear parsing
    @MASTER_INFILE_CONTENTS  = qx(cat ${MASTER_INFILE}) ;
    $NLINE_MASTER_INFILE     = scalar(@MASTER_INFILE_CONTENTS);

    if ( $SIMFLAG ) {       
	&parse_SIMTRAIN();  # --- parse SIMTRAIN keys ----       
	&parse_SIMDATA();   # --- parse SIMDATA keys ----
    }
    else {
	&parse_DATATRAIN();
	&parse_DATAFIT();
    }

    # parse for COSMO OPTIONS
    &parse_COSMO();

    if ( $DMCORRFLAG ) {
	&parse_DMCORR();  # --- parse DMCORR keys ---
    }

    # parse SALT2mu file
    $key = "INFILE_SALT2mu:" ;
    @tmp = sntools::parse_line($INFILE, 1, $key, $OPTABORT ) ;
    $INFILE_SALT2mu = $tmp[0];
    $name = $PROGRAM_NAME[$ISTAGE_SALT2mu] ;
    print "\n  $name input file (template) : $INFILE_SALT2mu \n" ;

    &parse_PAWPLOTMACRO();

    # read global keys for TRAINRUN (Jan 30 2012)
    my ($NKEY, $ikey);
    $NKEY = scalar(@GLOBAL_TRAINRUN_KEYNAMES);
    for ( $ikey=0; $ikey < $NKEY; $ikey++ ) {
	$key = $GLOBAL_TRAINRUN_KEYNAMES[$ikey] ;
	@tmp = sntools::parse_line($INFILE,99, $key, $OPTQUIET ) ;
	$GLOBAL_TRAINRUN_KEYVALS[$ikey] = "$tmp[0]" ;
    } 

    # Set global FIT_RUNOPT key
    #    to be used with split_and_fit.pl
    # Consider making this a key user input option at a later date.
    # JLM OCT 2013
    $GLOBAL_FIT_RUNOPT = " NOPROMPT ";

    # ---------------
    print "\n";


} # end of parse_master_infile



# ============================
sub parse_DATATRAIN {

    # read master input file linearly and look for DATATRAIN_MODEL keys.

    my ($iline, @wdlist, $tmpLine, $wd0, $wd1, $nwd, $c3 );
    my ($itrain, $KEY, $changes, $model );

    print "\n";
    print "  Parsing DATATRAIN keys ... \n";

    $NMODEL_TRAIN = 0;

    $iline = 0; 
    while ( $iline < $NLINE_MASTER_INFILE ) {

        $tmpLine = "$MASTER_INFILE_CONTENTS[$iline]" ;
        @wdlist  = split(/\s+/,$tmpLine) ;
        $wd0     = $wdlist[0] ;
        $wd1     = $wdlist[1] ;
        $nwd     = scalar(@wdlist);
	
	if ( "$wd0" eq "DATATRAIN_END:" ) { goto DONE_PARSE ; }

	if ( "$wd0" eq "DATATRAIN_MODEL:" ) {
	    $NMODEL_TRAIN++ ;
	    $itrain = $NMODEL_TRAIN ;
	    $TRAIN_MODELNAME[$itrain] = $wd1 ;
	    $c3     = sprintf("%3.3d", $itrain);
	    $TRAIN_CHARINDEX[$itrain] = "$c3" ;

	    # init Number of files and changes for this training
	    $NINFILE_TRAIN[$itrain][$JOBID_PCAFIT]    = 0 ;
	    $NCHANGE_TRAIN[$itrain][$JOBID_PCAFIT]    = 0 ;
	    $TRAIN_CHANGE[$itrain][$JOBID_PCAFIT][1]  = "" ;
	}

	$itrain = $NMODEL_TRAIN ;

	$KEY = "DATATRAIN_PATH:" ;
	if ( "$wd0" eq "$KEY" ) {
	    $DATATRAIN_PATH[$itrain] = $wd1 ;
	    if ( "$wd1" eq "[]" ) {
		$DATATRAIN_PATH[$itrain] = $DATATRAIN_PATH[$itrain-1] ;
	    }
	}

	$KEY = "DATATRAIN_LISTFILE:" ;
	if ( "$wd0" eq "$KEY" ) {	    
	    $DATATRAIN_LISTFILE[$itrain] = $wd1 ;
	    if ( "$wd1" eq "[]" ) {
		$DATATRAIN_LISTFILE[$itrain] = $DATATRAIN_LISTFILE[$itrain-1] ;
	    }
	}

	$KEY = "DATATRAIN_PCAFIT_INFILE:"; # input file	
	if ( "$wd0" eq "$KEY" ) {
	    &update_TRAIN_INFILE($itrain,$JOBID_PCAFIT,$wd1);
	}

	$KEY = "DATATRAIN_PCAFIT_CHANGE:"; # changes to input file
	if ( "$wd0" eq "$KEY" ) {
	    $changes = "@wdlist[1 .. $nwd]";
	    $changes  =~ s/\s+$// ;  # remove trailing blanks
	    &update_TRAIN_CHANGE($itrain,$JOBID_PCAFIT,$changes);
	}

	$iline++ ;
    }    

  DONE_PARSE:

    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {
	$model = $TRAIN_MODELNAME[$itrain] ;
	print "\t  DATATRAIN_MODEL[$itrain] = SALT2.$model \n";

    }

} # end of parse_DATATRAIN


# =========================================================
sub parse_DATAFIT {

    # parse DATAFIT_XXX keys

    my ($iline, $key, $tmpLine, $val, $nwd, @wdlist);
    my ($igroup, $VERSION, $fitnml, $NV, $arg, $GG );
    my (@FITNML_TMP);

    print "\n";
    print "  Parsing DATAFIT keys ... \n";

    $iline = 0;
    while ( $iline < $NLINE_MASTER_INFILE ) {

        $tmpLine = "$MASTER_INFILE_CONTENTS[$iline]" ;
        @wdlist  = split(/\s+/,$tmpLine) ;
        $key     = $wdlist[0] ;
        $val     = $wdlist[1] ;
        $nwd     = scalar(@wdlist);

	if ( "$key" eq "DATAFIT_END:" ) { goto DONE_PARSE ; }

	if ( $key eq "DATAFIT_GROUPNAME:" ) {
	    $NGROUP_DATA++ ;
	    $igroup = $NGROUP_DATA ;
	    $DATA_GROUPNAME[$igroup]  = $val ;
	    $NVERSION_GROUP[$igroup]  = 0 ;
	    $DATA_GROUPID[$igroup]    = "" ;
	    $DATAFIT_VERSION[$igroup] = "" ;
	    $DATA_INFILE_FIT[$igroup] = "" ;
	    $DATA_FITNML[$igroup]     = "" ;
	    print "     Group: $val \n";
	}

	if ( $key eq "DATAFIT_DEF:" ) {
	    $NVERSION_GROUP[$igroup]++ ;
	    $NV       = $NVERSION_GROUP[$igroup] ;
	    $VERSION  = $wdlist[1] ;
	    $fitnml   = $wdlist[2] ;
	    $FITNML_TMP[$NV] = $fitnml ;

	    print "\t VERSION: $VERSION \n";
	    # abort if this DATA version is already split
	    &check_split_version($VERSION);

	    if ( "$fitnml" eq "[]" ) {
		$fitnml = $FITNML_TMP[$NV-1] ;
		$FITNML_TMP[$NV] = $fitnml ;
	    }

	    $arg = $DATA_GROUPID[$igroup] ; 
            $DATA_GROUPID[$igroup] = "$arg" . "$NV " ;

	    $arg = $DATAFIT_VERSION[$igroup] ;
	    $DATAFIT_VERSION[$igroup] = "$arg" . "$VERSION " ;

            $arg  = $DATA_INFILE_FIT[$igroup] ;
            $DATA_INFILE_FIT[$igroup] = "$arg" . "$fitnml " ;

	    $GG   = sprintf("GROUP%2.2d", $igroup); 
            $arg  = $DATA_FITNML[$igroup] ;
            $DATA_FITNML[$igroup] = "$arg" . "${GG}_${VERSION}.nml " ;
	}

	$iline++ ;
    }

  DONE_PARSE:

} # end of parse_DATAFIT

# =========================================================
sub parse_SIMTRAIN {

    # read master input file linearly and look for 
    # SIMTRAIN_MODEL key followed by keys
    #   SIMTRAIN_XXX_INFILE
    #   SIMTRAIN_XXX_CHANGE
    # where XXX = GENPHOT, GENSPEC, INCDATA, PCAFIT, TRAINLIST
    # correspond to the five programs.
    #


    my ($iline, $KEY, $wd0, $wd1, $val, $nwd, $jobid, $JOB, $itrain);
    my ($NTMP,$i,$job, $TMP, @TMP, $tmpLine, @wdlist, $changes, $c3);
    my ($NGEN, $igen, $cgen, $GENV, $model);
    my ($simFile, $simFile_last, $simArgs);
    my ($incdataFile, $incdataArgs);

    print "\n";

    print "  Parsing SIMTRAIN keys ... \n";

    $NMODEL_TRAIN = 0;

    $iline = 0; 
    while ( $iline < $NLINE_MASTER_INFILE ) {

        $tmpLine = "$MASTER_INFILE_CONTENTS[$iline]" ;
        @wdlist  = split(/\s+/,$tmpLine) ;
        $wd0     = $wdlist[0] ;
        $wd1     = $wdlist[1] ;
        $nwd     = scalar(@wdlist);
	
	if ( "$wd0" eq "SIMTRAIN_END:" ) { goto DONE_PARSE ; }

	if ( "$wd0" eq "SIMTRAIN_MODEL:" ) {
	    $NMODEL_TRAIN++ ;
	    $itrain = $NMODEL_TRAIN ;
	    $TRAIN_MODELNAME[$itrain] = $wd1 ;

	    $c3     = sprintf("%3.3d", $itrain);
	    $TRAIN_CHARINDEX[$itrain] = "$c3" ;

	    # init Number of files and changes for this training
	    for ( $jobid = 1; $jobid <= $NJOBID; $jobid++ )  {
		$NINFILE_TRAIN[$itrain][$jobid]    = 0 ;
		$NCHANGE_TRAIN[$itrain][$jobid]    = 0 ;
		$TRAIN_CHANGE[$itrain][$jobid][1]  = "" ;
		$TRAIN_INFILE[$itrain][$jobid][1]  = "" ;
		$SIMTRAIN_REPEATFLAG[$itrain][$jobid] = 0 ;
	    }
	}

	# loop over job-types and check keys in a uniform mannor

	$itrain = $NMODEL_TRAIN ;

	for ( $jobid = 1; $jobid <= $NJOBID; $jobid++ ) {
	    $JOB = $JOBID_NAME[$jobid];

	    $KEY = "SIMTRAIN_${JOB}_INFILE:"; # input file	
	    if ( "$wd0" eq "$KEY" ) {
		&update_TRAIN_INFILE($itrain,$jobid,$wd1);
	    }

	    $KEY = "SIMTRAIN_${JOB}_CHANGE:"; # changes to input file
	    if ( "$wd0" eq "$KEY" ) {
		$changes = "@wdlist[1 .. $nwd]";
		$changes  =~ s/\s+$// ;  # remove trailing blanks
		&update_TRAIN_CHANGE($itrain,$jobid,$changes);
	    }

	    $KEY = "SIMTRAIN_${JOB}_REPEAT:"; # repeat from previous training
	    if ( "$wd0" eq "$KEY" ) {
		$SIMTRAIN_REPEATFLAG[$itrain][$jobid] = 1 ;
		&repeat_SIMTRAIN($itrain,$jobid);
	    }

	}  # jobid loop

	$iline++ ;
    }    

  DONE_PARSE:
    

    # for each GENPHOT extract the GENVERSION by checking both
    # the sim-input file and the GENPHOT args.
    #
    # for each INCDATA extract the GENVERSION by checking the
    # INCDATA args.
    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {

	$NGEN =   $NINFILE_TRAIN[$itrain][$JOBID_GENPHOT] ;
	for ( $igen=1; $igen <= $NGEN; $igen++ ) {
	    $simFile = "$TRAIN_INFILE[$itrain][$JOBID_GENPHOT][$igen]"; 
	    $simArgs = "$TRAIN_CHANGE[$itrain][$JOBID_GENPHOT][$igen]";

	    $GENV    = &getGENVERSION($itrain,$igen,$simArgs,$JOBID_GENPHOT);
	    $SIMTRAIN_GENVERSION[$itrain][$igen] = "$GENV" ;

	    if ( $simFile ne $simFile_last ) { &check_simFile($simFile); }
	    $simFile_last = $simFile ;
	}

	$NGEN = $NINFILE_TRAIN[$itrain][$JOBID_INCDATA] ;
	for ($igen=1; $igen <= $NGEN; $igen++ ){
	    $incdataArgs = "$TRAIN_CHANGE[$itrain][$JOBID_INCDATA][$igen]";
	    $GENV = &getGENVERSION($itrain,$igen,$incdataArgs,$JOBID_INCDATA);
	    $SIMTRAIN_INCDATA_GENVERSION[$itrain][$igen] = "$GENV" ;
	}
    }



  SUMMARY:      # summarize and make sanity checks

    # if batch system, make sure there are more nodes than trains
    if ( $BATCH_FLAG > 0 ) {
	check_parallel_batch($STAGE_NAME[$ISTAGE_TRAINGEN], $NMODEL_TRAIN, $BATCH_NODE_TRAIN);
    }

    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {
	$model = $TRAIN_MODELNAME[$itrain] ;
	print "  SIMTRAIN_MODEL[$itrain] = SALT2.$model \n";

	# make sure we have input files for each job type
	for ( $jobid = 1; $jobid <= $NJOBID; $jobid++ ) {
	    $NTMP =   $NINFILE_TRAIN[$itrain][$jobid] ;
	    if ( $NTMP <= 0 ) {
		$JOB = $JOBID_NAME[$jobid];
		$KEY = "SIMTRAIN_${JOB}_INFILE" ;
		$MSGERR[0] = "No  $KEY  key(s) for this training.";
		$MSGERR[1] = "Check master input file: ${MASTER_INFILE} ";
		sntools::FATAL_ERROR(@MSGERR);
	    }
	}



	if ( $VBOSEFLAG == 0 ) { next ; }

	$NGEN =   $NINFILE_TRAIN[$itrain][$JOBID_GENPHOT] ;

	for ( $igen=1; $igen <= $NGEN; $igen++ ) {	   
	    $cgen = $igen ;
	    
	    $TMP = "$TRAIN_INFILE[$itrain][$JOBID_GENPHOT][$igen]";
	    print "\t GENPHOT[$cgen] input file : $TMP \n";
	    
	    $TMP = "$TRAIN_CHANGE[$itrain][$JOBID_GENPHOT][$igen]";
	    print "\t GENPHOT[$cgen] change args: $TMP \n";
	    
	    $TMP = "$SIMTRAIN_GENVERSION[$itrain][$igen]";
	    print "\t GENPHOT[$cgen] GENVERSION: $TMP \n";
	}  # igen loop
	

	$NGEN  = $NINFILE_TRAIN[$itrain][$JOBID_INCDATA] ;
	
	for ($igen=1; $igen <= $NGEN; $igen++ ) {
	    $cgen = $igen ;

	    $TMP = "$TRAIN_INFILE[$itrain][$JOBID_INCDATA][$igen]";
	    print "\t INCDATA[$cgen] input file : $TMP \n";
	    
	    $TMP = "$TRAIN_CHANGE[$itrain][$JOBID_INCDATA][$igen]";
	    print "\t INCDATA[$cgen] change args: $TMP \n";
	    
	    $TMP = "$SIMTRAIN_INCDATA_GENVERSION[$itrain][$igen]";
	    print "\t INCDATA[$cgen] GENVERSION: $TMP \n";
	} #igen loop

	$TMP = "$TRAIN_INFILE[$itrain][$JOBID_GENSPEC][1]" ;
	print "\t GENSPEC input file : $TMP \n";
	$TMP = "$TRAIN_CHANGE[$itrain][$JOBID_GENSPEC][1]" ;
	print "\t GENSPEC change args: $TMP \n";

	$TMP = "$TRAIN_INFILE[$itrain][$JOBID_PCAFIT][1]" ;
	print "\t PCAFIT input file : $TMP \n";
	$TMP = "$TRAIN_CHANGE[$itrain][$JOBID_PCAFIT][1]" ;
	print "\t PCAFIT change args: $TMP \n";

	$TMP = "$TRAIN_INFILE[$itrain][$JOBID_TRAINLIST][1]" ;
	print "\t TRAIN_LIST input file : $TMP \n";

    }  # itrain

} # end of parse_SIMTRAIN


# ==================================
sub repeat_SIMTRAIN {

    my ($itrain, $jobid) = @_ ;

    my ($NTMP, $JOB, $TMP, $i) ;

    # store INFILE and CHANGE from previous training.

    if ( $itrain <= 1 ) {	
	$JOB = $JOBID_NAME[$jobid];
	$MSGERR[0] = "Cannot repeat SIMTRAIN_$JOB args on first training.";
	sntools::FATAL_ERROR(@MSGERR);
    }
   
    $NTMP =  $NINFILE_TRAIN[$itrain-1][$jobid] ;
    $NINFILE_TRAIN[$itrain][$jobid] = $NTMP ;

    for ( $i=1; $i <= $NTMP; $i++ ) {	   	    
	$TMP = "$TRAIN_INFILE[$itrain-1][$jobid][$i]";
	$TRAIN_INFILE[$itrain][$jobid][$i] = "$TMP"; 

	$TMP = "$TRAIN_CHANGE[$itrain-1][$jobid][$i]";
	$TRAIN_CHANGE[$itrain][$jobid][$i] = "$TMP"; 	    
    }  


} # repeat_SIMTRAIN



# ======================================
sub update_TRAIN_INFILE(@) {
    my ($itrain, $jobid, $inFile) = @_ ;

    my ($NFILE, $JOB, $c3, $inFile2, $model, $NEWFILE) ;

    $NINFILE_TRAIN[$itrain][$jobid]++ ;

    if ( $jobid == $JOBID_GENSPEC ) {
	$NFILE = $NINFILE_TRAIN[$itrain][$JOBID_GENPHOT] ;
    }
    else {
	$NFILE = $NINFILE_TRAIN[$itrain][$jobid] ;
    }

    if ( "$inFile" eq "[]" && $NFILE == 1 ) {
	$NEWFILE = $TRAIN_INFILE[$itrain-1][$jobid][$NFILE] ;
    }
    elsif ( "$inFile" eq "[]" && $NFILE > 1 ) {
	$NEWFILE = $TRAIN_INFILE[$itrain][$jobid][$NFILE-1] ;
    }
    else {
	$NEWFILE = "$inFile" ;
    }

    $TRAIN_INFILE[$itrain][$jobid][$NFILE] = "$NEWFILE";

    # for pcafit, create 2nd filename so that each
    # training has a distinct pcafit file. This is needed
    # because there are no command-line overrides for pcafit.

    if ( $jobid == $JOBID_PCAFIT ) {

	$JOB = $JOBID_NAME[$jobid];

	if ( $NFILE != 1 ) {
	    $MSGERR[0] = "Invalid NINFILE[itrain=$itrain],$JOB] = $NFILE";
	    $MSGERR[1] = "There can be only 1 $JOB-config file per training.";
	    &FATAL_ERROR(@MSGERR);
	}

	$model   = $TRAIN_MODELNAME[$itrain] ;
	$c3      = $TRAIN_CHARINDEX[$itrain];
	$inFile2 = "pcafit${c3}_${model}.config" ;
	$TRAIN_INFILE[$itrain][$jobid][2] = "$inFile2" ;
    }


}  # end of update_TRAIN_INFILE

# ======================================
sub update_TRAIN_CHANGE(@) {

    # Called from parse_SIMTRAIN to update arrays to store change(s) 
    # for this itrain and jobid.
    # if changes = [] then store change from previous training.

    my ($itrain, $jobid, $changes) = @_ ;
    my ($NFILE, $TMP, $CHANGE) ;

    $NCHANGE_TRAIN[$itrain][$jobid]++ ;

    if ( $jobid == $JOBID_GENSPEC ) {
	$NFILE = $NINFILE_TRAIN[$itrain][$JOBID_GENPHOT] ;
    }
    else {
	$NFILE = $NINFILE_TRAIN[$itrain][$jobid] ;
    }

    if ( "$changes" eq "[]" && $NFILE == 1 ) {
	$CHANGE = "$TRAIN_CHANGE[$itrain-1][$jobid][$NFILE]";
    }
    elsif ( "$changes" eq "[]" && $NFILE > 1 ) {
	$CHANGE = "$TRAIN_CHANGE[$itrain][$jobid][$NFILE-1]";
    }
    else {
	$TMP  = "$TRAIN_CHANGE[$itrain][$jobid][$NFILE]" ;
	$CHANGE = "$TMP" . "$changes " ;
    }


    $TRAIN_CHANGE[$itrain][$jobid][$NFILE] = $CHANGE ;

}  # end of update_TRAIN_CHANGE

# ======================================
sub parse_SIMDATA {

    # look for key  "SALT2mu_GROUPNAME: <name>" followed by
    # INFILE_SIMGEN, RUNOPT, GENOPT, and (optional) 
    # SIMDATA_EXCLUDE_FITOPT keys for this group.
    # Note that SIMDATA_XXX is an array over NGROUP_DATA
    # while SIMGEN_XXX is an array over all sim-samples.
    # NSIMGEN >= NGROUP_DATA.
    # To allow multiple jobs with the same SIMDATA-GENVERSION,
    # the GENVERION is appended with the number of seconds
    # since midnight.
    # If SIMDATA_EXCLUDE_FITOPT: is set for a GROUPNAME, 
    # only default model fits will be performed for that group.
    # No trained model fits will be performed. 

    my ($iline, $key, $val, $nwd, $tmpLine, @wdlist, $arg) ;
    my ($name, $igroup, @tmp, $c3n, $nfiles );
    my ($simFile, $fitFile, $GENVERSION );
    my (@T_now, $Nsec, $Nsec5, $GG);

    $NGROUP_DATA  = 0;
    $NSIMGEN         = 0;

    # check if user wants a particular integer id for
    # debugging a previous generation where the 
    # integer id is already specified.
    $key = "SIMDATA_INTID:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET ) ;
    if ( scalar(@tmp) > 0 ) {
	$Nsec = $tmp[0] ;
    }
    else {
	# get number of seconds since midnight to append to GENVERSION
	@T_now = localtime();
	$Nsec  = ($T_now[2] * 3600) + ($T_now[1] * 60) + $T_now[0] ;
    }
    $Nsec5 = sprintf("%5.5d", $Nsec);
    
    print "\n";

    $iline = 0;
    while ( $iline < $NLINE_MASTER_INFILE ) {

        $tmpLine = "$MASTER_INFILE_CONTENTS[$iline]" ;
        @wdlist  = split(/\s+/,$tmpLine) ;
        $key     = $wdlist[0] ;
        $val     = $wdlist[1] ;
        $nwd     = scalar(@wdlist);

	if ( $key eq "SIMDATA_GROUPNAME:" ) {
	    $NGROUP_DATA++ ;
	    $igroup = $NGROUP_DATA ;
	    $DATA_GROUPNAME[$igroup] = $val ;
	    $SIMDATA_EXCLUDE_FITOPT[$igroup] = 0; #initialize to 0
	}

	if ( $key eq "SIMDATA_INFILES:" ) {

	    $simFile = $wdlist[1] ;
	    $fitFile = $wdlist[2] ;

            $arg  = $SIMDATA_INFILE_GEN[$igroup] ;
            $SIMDATA_INFILE_GEN[$igroup] = "$arg" . "$simFile " ;

            $arg  = $DATA_INFILE_FIT[$igroup] ;
            $DATA_INFILE_FIT[$igroup] = "$arg" . "$fitFile " ;

	    # update list of internal gen-indices
	    $NSIMGEN++ ;
            $c3n = sprintf("%3.3d", $NSIMGEN );

	    $arg = $DATA_GROUPID[$igroup] ; 
            $DATA_GROUPID[$igroup] = "$arg" . "$NSIMGEN " ;

	    # --------
	    # update list of gen-versions and fit-nml file for each version
	    @tmp = sntools::parse_line($simFile,1, "GENVERSION:",$OPTABORT );
	    $GENVERSION             = "$tmp[0]-$Nsec5-$c3n" ;
	    $arg = $SIMDATA_GENVERSION[$igroup] ;
	    $SIMDATA_GENVERSION[$igroup] = "$arg" . "$GENVERSION " ;
	    # --------

	    $GG   = sprintf("GROUP%2.2d", $igroup); 
            $arg  = $DATA_FITNML[$igroup] ;
            $DATA_FITNML[$igroup] = "$arg" . "${GG}_${GENVERSION}.nml ";
	    $NVERSION_GROUP[$igroup]++ ;
	}

	if ( $key eq "SIMDATA_EXCLUDE_FITOPT:" ) {
	    $SIMDATA_EXCLUDE_FITOPT[$igroup] = 1;
	}

	if ( $key eq "SIMDATA_GENOPT:" ) {
            $arg  = $SIMDATA_GENOPT[$igroup] ;
            $SIMDATA_GENOPT[$igroup] = "$arg" . "@wdlist[1 .. $nwd] " ;
	}

	if ( $key eq "SIMDATA_RUNOPT:" ) {
            $arg  = $SIMDATA_RUNOPT[$igroup] ;
            $SIMDATA_RUNOPT[$igroup] = "$arg" . "@wdlist[1 .. $nwd] " ;
	}
	
	$iline++ ;
    }  # end of iline

    # if batch requested, check that number of nodes > number of sims

    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) {
	$name = $STAGE_NAME[$ISTAGE_SIMGEN] . " " . $DATA_GROUPNAME[$igroup] ;
	@tmp  = split(/\s+/,$SIMDATA_INFILE_GEN[$igroup]) ;
	$nfiles = @tmp;

	if ( $BATCH_FLAG > 0 ) {
	check_parallel_batch($name, $nfiles, $BATCH_NODE_TRAIN);
	}
    }

    # summarize
    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) {
	$name = $DATA_GROUPNAME[$igroup] ;

	print "  SIMDATA_GROUP '${name}' \n";
	print "\t GEN-INFILES : $SIMDATA_INFILE_GEN[$igroup] \n";
	print "\t FIT-INFILES : $DATA_INFILE_FIT[$igroup] \n";
	print "\t RUN-OPTIONS : $SIMDATA_RUNOPT[$igroup] \n";
	print "\t GEN-OPTIONS : $SIMDATA_GENOPT[$igroup] \n";
	print "\t FIT-OPTIONS : $SIMDATA_EXCLUDE_FITOPT[$igroup] \n";
	print "\t               (1: train model fits off) \n\n";
	# print "\t IDGEN       : $DATA_GROUPID[$igroup] \n";
    }


    printf "\n";

} # end of parse_SIMDATA

# ======================================
sub parse_DMCORR {

    ######
    ### BIAS CORRECTIONS WORK AS FOLLOWS:
    ###  1) simulated SN data are simulated to match cosmology data
    ###  2) simulated SN data are fit in same manner as cosmo data
    ###  3) a comparison between sim and fit distances is used to determine
    ###     bias corrections vs redshift, cosmology data is corrected
    ######
    ###
    ### GROUP KEYS
    ###
    # DMCORR_GROUPNAME: <name to call group>
    # DMCORR_TYPE: <type of bias corr - P11 or R13>
    # DMCORR_COMPNAME: <name of simdata group to correct>
    ###
    ### SIM-RELATED KEYWORDS
    ###
    # SIMDMCORR_INFILES: <simgen input file> <simgen nml file>
    # DMCORR_SIMALPHA: <simalpha value used for simdata>
    # DMCORR_SIMBETA: <simbeta value used for simdata>
    # (opt) DMCORR_SIMINFO: < TRAINS cmean sig_(-)c sig_(+)c x1mean sig_(-)x1 sig_(+)x1 >
    # (opt) DMCORR_GENOPT: <sim generation options>
    # (opt) DMCORR_RUNOPT: <sim code command line options>
    ###
    ### CORRECTION-RELATED KEYWORDS
    ###
    # (opt) DMCORR_RUNCODE: <alternative file to use for bias calculations> 
    #                        SALT2biascorr.pl is default. Set at group level.
    # (opt) DMCORR_RUNCODE_GLOB: <alternative file to use for ALL bias calculations> 
    #                        SALT2biascorr.pl is default. Set for ALL groups.
    # (opt) DMCORR_RUNCODEOPT: <options to use for user input DMCORR_RUNCODE>
    ###
    ### MISCELLANEOUS KEYWORDS
    ###
    # (opt) DMCORR_SIMCODE: <alternative code to extend SALT2 train color dispersion file>
    # (opt) DMCORR_SIMCODEOPT: <options to use for user input DMCORR_SIMCODE>
    # and 
    # (opt) DMCORR_EXCLUDE_FITOPT: <same as SIMGEN_EXCLUDE_FITOPT>
    # keys for this group.
    ###
    ###
    ###
    # Note that SIMDMCORR_XXX is an array over NGROUP_DMCORR
    # while SIMGENDMCORR_XXX is an array over all sim-samples.
    # NSIMGENDMCORR >= NGROUP_DMCORR.
    # To allow multiple jobs with the same SIMDATA-GENVERSION,
    # the GENVERSION is appended with the number of seconds
    # since midnight.
    # If DMCORR_EXCLUDE_FITOPT: is set for a GROUPNAME, 
    # only default model fits will be performed for that group.
    # No trained model fits will be performed. 
    ####

    my ($iline, $key, $val, $nwd, $tmpLine, @wdlist, $arg) ;
    my ($name, $type, $simcomp, $igroup, @tmp, $c3n );
    my ($simFile, $fitFile, $GENVERSION );
    my (@T_now, $Nsec, $Nsec5, $GG);
    my ($tempcomp, $itemp, $matchval);
    my ($glob_runcode);

    $NGROUP_DMCORR  = 0;
    $NSIMGENDMCORR  = 0;

    # get number of seconds since midnight + 100 to append to GENVERSION
      @T_now = localtime();
      $Nsec  = ($T_now[2] * 3600) + ($T_now[1] * 60) + $T_now[0] + 100;
      $Nsec5 = sprintf("%5.5d", $Nsec);
    
    print "\n";

    $iline = 0;
    while ( $iline < $NLINE_MASTER_INFILE ) {

        $tmpLine = "$MASTER_INFILE_CONTENTS[$iline]" ;
        @wdlist  = split(/\s+/,$tmpLine) ;
        $key     = $wdlist[0] ;
        $val     = $wdlist[1] ;
        $nwd     = scalar(@wdlist);

	if ( $key eq "DMCORR_GROUPNAME:" ) {
	    $NGROUP_DMCORR++ ;
	    $igroup = $NGROUP_DMCORR ;
	    $DMCORR_GROUPNAME[$igroup] = $val ;
	    $DMCORR_EXCLUDE_FITOPT[$igroup] = 0; #initialize to 0
	    $DMCORR_SIMTYPE[$igroup] = "AS-IS"; #initialize to AS-IS
	    $DMCORR_IDEAL[$igroup] = -9; #initialize to -9
	    $DMCORR_RUNCODE[$igroup] = "DEFAULT"; #initialize to DEFAULT
	    $DMCORR_SIMCODE[$igroup] = "DEFAULT"; #initialize to DEFAULT
	    $DMCORR_SIMALPHA[$igroup] = 0.097; #initialize to 0.11*0.88
	    $DMCORR_SIMBETA[$igroup] = 3.2; #initialize to 3.2
	}

	if ( $key eq "DMCORR_TYPE:" ) {
	    $DMCORR_TYPE[$igroup] = $val ;
	}

	if ( $key eq "DMCORR_COMPNAME:" ) {
	    $tempcomp = $val ;
	    $matchval = -9 ; #default group number
	    ## FIND igroup number corresponding to tempcomp string
	    for ( $itemp = 1; $itemp <= $NGROUP_DATA; $itemp++ ) {
		$name = $DATA_GROUPNAME[$itemp] ;
		if ( $name eq $tempcomp ) { $matchval = $itemp; }
	    }
	    $DMCORR_IDEAL[$igroup] = $matchval ;
	}
	
	if ( $key eq "DMCORR_SIMALPHA:" ) {
	    $DMCORR_SIMALPHA[$igroup] = $val ;
	}

	if ( $key eq "DMCORR_SIMBETA:" ) {
	    $DMCORR_SIMBETA[$igroup] = $val ;
	}

	if ( $key eq "DMCORR_RUNCODE:" ) {
	    $DMCORR_RUNCODE[$igroup] = $val ;
	}

	if ( $key eq "DMCORR_RUNCODE_GLOB:" ) {
	    $glob_runcode = $val ;
	}

	if ( $key eq "DMCORR_SIMCODE:" ) {
	    $DMCORR_SIMCODE[$igroup] = $val ;
	}

	if ( $key eq "DMCORR_RUNCODEOPT:" ) {
	    if ( "$wdlist[1]" eq "[]" ) {
		$DMCORR_RUNCODEOPT[$igroup] = $DMCORR_RUNCODEOPT[$igroup-1] ;
	    } else {
		$arg  = $DMCORR_RUNCODEOPT[$igroup] ;
		$DMCORR_RUNCODEOPT[$igroup] = "$arg" . "@wdlist[1 .. $nwd] " ;
	    }
	}

	if ( $key eq "DMCORR_SIMCODEOPT:" ) {
            $arg  = $DMCORR_SIMCODEOPT[$igroup] ;
            $DMCORR_SIMCODEOPT[$igroup] = "$arg" . "@wdlist[1 .. $nwd] " ;
	}

	if ( $key eq "DMCORR_SIMINFO:" ) {
	    $DMCORR_SIMTYPE[$igroup] = $wdlist[1] ;
	    $DMCORR_MEANC[$igroup] = $wdlist[2] ;
	    $DMCORR_SIGLC[$igroup] = $wdlist[3] ;
	    $DMCORR_SIGRC[$igroup] = $wdlist[4] ;
	    $DMCORR_MEANX[$igroup] = $wdlist[5] ;
	    $DMCORR_SIGLX[$igroup] = $wdlist[6] ;
	    $DMCORR_SIGRX[$igroup] = $wdlist[7] ;
	}

	if ( $key eq "SIMDMCORR_INFILES:" ) {

	    $simFile = $wdlist[1] ;
	    $fitFile = $wdlist[2] ;

            $arg  = $DMCORR_INFILE_GEN[$igroup] ;
            $DMCORR_INFILE_GEN[$igroup] = "$arg" . "$simFile " ;

            $arg  = $DMCORR_INFILE_FIT[$igroup] ;
            $DMCORR_INFILE_FIT[$igroup] = "$arg" . "$fitFile " ;

	    # update list of internal gen-indices
	    $NSIMGENDMCORR++ ;
            $c3n = sprintf("%3.3d", $NSIMGENDMCORR );

	    $arg = $DMCORR_GROUPID[$igroup] ; 
            $DMCORR_GROUPID[$igroup] = "$arg" . "$NSIMGENDMCORR " ;

	    # --------
	    # update list of gen-versions and fit-nml file for each version
	    @tmp = sntools::parse_line($simFile,1, "GENVERSION:",$OPTABORT );
	    $GENVERSION             = "$tmp[0]-$Nsec5-$c3n" ;
	    $arg = $DMCORR_GENVERSION[$igroup] ;
	    $DMCORR_GENVERSION[$igroup] = "$arg" . "$GENVERSION " ;
	    # --------

	    $GG   = sprintf("GROUP%2.2d", $igroup); 
            $arg  = $DMCORR_FITNML[$igroup] ;
            $DMCORR_FITNML[$igroup] = "$arg" . "${GG}_${GENVERSION}.nml ";
	    $NVERSION_GROUP_DMCORR[$igroup]++ ;
	}

	if ( $key eq "DMCORR_EXCLUDE_FITOPT:" ) {
	    $DMCORR_EXCLUDE_FITOPT[$igroup] = 1;
	}

	if ( $key eq "DMCORR_GENOPT:" ) {
            $arg  = $DMCORR_GENOPT[$igroup] ;
            $DMCORR_GENOPT[$igroup] = "$arg" . "@wdlist[1 .. $nwd] " ;
	}

	if ( $key eq "DMCORR_RUNOPT:" ) {
            $arg  = $DMCORR_GENOPT[$igroup] ;
            $DMCORR_RUNOPT[$igroup] = "$arg" . "@wdlist[1 .. $nwd] " ;
	}
	
	$iline++ ;
    }  # end of iline

    # dummy check
    for ( $igroup = 1; $igroup <= $NGROUP_DMCORR; $igroup++ ) {
	if ( $glob_runcode ne "" ) {
	    if ( $DMCORR_RUNCODE[$igroup] ne "DEFAULT" ) {
		$MSGERR[0] = "CONFLICT IN DMCORR_RUNCODE keys";
		$MSGERR[1] = "Cannot set DMCORR_RUNCODE_GLOB:";
		$MSGERR[2] = " and DMCORR_RUNCODE: simultaneously";
		$MSGERR[3] = "CHECK INPUT FILE";
	    } else {
		$DMCORR_RUNCODE[$igroup] = $glob_runcode;
	    }
	}
    }

    # summarize
    for ( $igroup = 1; $igroup <= $NGROUP_DMCORR; $igroup++ ) {
	$name = $DMCORR_GROUPNAME[$igroup] ;
	$type = $DMCORR_TYPE[$igroup];
	$simcomp = $DATA_GROUPNAME[$DMCORR_IDEAL[$igroup]];
	
	print "  DMCORR_GROUP '${name}' \n";
	print "\t DMCORR-TYPE     : ${type} \n";
	print "\t DMCORR-SIMTYPE  : $DMCORR_SIMTYPE[$igroup] \n";
	print "\t DMCORR-COMP     : '${simcomp}' \n";
	print "\t DMCORR SIMINFO: \n";
	print "\t\t DMCORR-SIMCODE  : $DMCORR_SIMCODE[$igroup] \n";
	print "\t\t DMCORR-RUNOPT   : $DMCORR_RUNOPT[$igroup] \n";
	print "\t\t SIMCODE-OPTIONS : $DMCORR_SIMCODEOPT[$igroup] \n";
	print "\t\t GEN-INFILES     : $DMCORR_INFILE_GEN[$igroup] \n";
	print "\t\t GEN-OPTIONS     : $DMCORR_GENOPT[$igroup] \n";
	print "\t DMCORR FITINFO: \n";
	print "\t\t FIT-INFILES     : $DMCORR_INFILE_FIT[$igroup] \n";
	print "\t\t FIT-OPTIONS     : $DMCORR_EXCLUDE_FITOPT[$igroup] \n";
	print "\t\t                   (1: train model fits off) \n\n";
	print "\t DMCORR BIASCORR INFO: \n";
	print "\t\t DMCORR-RUNCODE  : $DMCORR_RUNCODE[$igroup] \n";
	print "\t\t RUNCODE-OPTIONS : $DMCORR_RUNCODEOPT[$igroup] \n\n";
	# print "\t IDGEN         : $DATA_GROUPID[$igroup] \n";
    }


    printf "\n";

} # end of parse_DMCORR



# ======================================
sub parse_COSMO {

    # default program is wfit.exe
    # look for optional key "COSMO_PROGRAM: <input>"
    # to use alternate program
    #
    # look for optional key "COSMO_GENOPT: <input>"
    # the same options will be used for all data sets
    # options need to be in appropriate cosmo program syntax
    # they will be pasted as is into cosmo program call
    # 
    # look for optional key "COSMO_INFILE: <input>"
    # the same file will be used for all data sets
    # this file is not idiot-checked, use at own risk

    my ($iline, $key, $val, $nwd, $tmpLine, @wdlist, $arg) ;
    my ($name, $igroup, @tmp, $c3n );

    $iline = 0;
    while ( $iline < $NLINE_MASTER_INFILE ) {

        $tmpLine = "$MASTER_INFILE_CONTENTS[$iline]" ;
        @wdlist  = split(/\s+/,$tmpLine) ;
        $key     = $wdlist[0] ;
        $val     = $wdlist[1] ;
        $nwd     = scalar(@wdlist);

	if ( $key eq "COSMO_PROGRAM:" ) {
            $PROGRAM_NAME[$ISTAGE_COSMO] =  $val ;
	    $OPT_COSMO = 1;
	}

	if ( $key eq "COSMO_INFILE:" ) {
	    $COSMO_INFILE = $val ;
	}

	if ( $key eq "COSMO_GENOPT:" ) {
	    $arg = $COSMO_GENOPT;
	    $COSMO_GENOPT = "$arg" . "@wdlist[1 .. $nwd] " ;
	}

        $iline++ ;
    }  # end of iline

    # idiot check for unconfigured cosmo programs
    unless ($PROGRAM_NAME[$ISTAGE_COSMO] =~ /wfit.exe$/) {
	unless ($PROGRAM_NAME[$ISTAGE_COSMO] =~ /sncosmo_mcmc.exe$/ ) {
	    $MSGERR[0] = "Invalid cosmo program $PROGRAM_NAME[$ISTAGE_COSMO].";
	    $MSGERR[1] = "Compatible programs include wfit.exe,";
	    $MSGERR[2] = " sncosmo_mcmc.exe.";
	    $MSGERR[3] = "Rerun with valid cosmo program.";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }
    # idiot check for required input file
    if ($PROGRAM_NAME[$ISTAGE_COSMO] =~ /sncosmo_mcmc.exe$/ ) {
	if ($COSMO_INFILE eq "") {
	    $MSGERR[0] = "Cosmo fits with sncosmo_mcmc.exe";
	    $MSGERR[1] = " require valid COSMO-INFILE inputs.";
	    $MSGERR[2] = "Set COSMO_INFILE: keyword and rerun.";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }
    # idiot check for keyword conflicts
    # 1) COSMO_PROGRAM off, COSMO_GENOPT on
    if ( $OPT_COSMO==0 && $COSMO_GENOPT ne "" ) {
	$MSGERR[0] = "COSMO GENOPT are given but COSMO PROGRAM is missing.";
	$MSGERR[1] = "Must remove COSMO_GENOPT from input";
	$MSGERR[2] = " or set a valid COSMO_PROGRAM value to run pipeline.";
	sntools::FATAL_ERROR(@MSGERR);
    }
    # 2) COSMO_PROGRAM off, COSMO_INFILE on
    if ($OPT_COSMO==0 && $COSMO_INFILE ne "" ) {
	$MSGERR[0] = "COSMO INFILE is defined but OPT_COSMO is 0";
	$MSGERR[1] = "Must remove COSMO INFILE AND OPTIONS from input";
	$MSGERR[2] = " or set OPT_COSMO: 1 to run pipeline.";
	sntools::FATAL_ERROR(@MSGERR);
    }
    # 3) COSMO_PROGRAM wfit.exe, COSMO_INFILE on
    if ( $PROGRAM_NAME[$ISTAGE_COSMO] =~ /wfit.exe$/ && $COSMO_INFILE ne "" ) { 
	$MSGERR[0] = "Cannot run COSMO_PROGRAM wfit.exe with COSMO_INFILE.";
	$MSGERR[1] = "Must change COSMO_PROGRAM to sncosmo_mcmc.exe OR";
	$MSGERR[2] = " remove COSMO_INFILE option.";
	$MSGERR[3] = "For wfit.exe, use optional COSMO_GENOPT keys.";
	sntools::FATAL_ERROR(@MSGERR);
    }
    # 4) COSMO_GENOPT on, COSMO-PROGRAM sncosmo_mcmc
    if ($PROGRAM_NAME[$ISTAGE_COSMO] =~ /sncosmo_mcmc.exe$/ ) {
	if ( $COSMO_GENOPT ne "" ) {
	    $MSGERR[0] = "COSMO PROGRAM incompatible with COSMO-GENOPT.";
	    $MSGERR[1] = "Must remove COSMO-GENOPT and add COSMO_INFILE to run pipeline.";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }




        # summarize
    if ( $OPT_COSMO ){
	print " COSMOLOGY FITTING \n";
	print "\t COSMO-PROGRAM: $PROGRAM_NAME[$ISTAGE_COSMO] \n";
	print "\t COSMO-INFILE: $COSMO_INFILE \n";
	print "\t COSMO-OPTIONS: $COSMO_GENOPT \n";
    } else {
	print " SKIPPING COSMOLOGY FITTING \n";
    }
    printf "\n";


} # end of parse_COSMO


# ==============================
sub parse_PAWPLOTMACRO {

    my ($key, @tmp, @tmp2, @tmpList, $NTMP, @wdlist, $i);
    my ($NAME, $ARG, $ARGSTAT );

    # default paw command should be 'paw', 
    # but some systems may be messed up and need something like pawX11
    $key = "PAW_COMMAND:" ;
    @tmp = sntools::parse_line($MASTER_INFILE, 1, $key, $OPTQUIET ) ;
    if ( scalar(@tmp) > 0 ) { $PAW_COMMAND = $tmp[0] ; };

    # allow either of 2 keys for PLOTMACRO
    $key  = "PAWPLOTMACRO:" ;
    @tmp  = sntools::parse_line($MASTER_INFILE, 2, $key, $OPTQUIET ) ;
    $key  = "PAW_PLOTMACRO:" ;
    @tmp2 = sntools::parse_line($MASTER_INFILE, 2, $key, $OPTQUIET ) ;

    @tmpList = ( @tmp , @tmp2 );
    $NTMP = scalar(@tmpList);    $N_PAWPLOTMACRO = $NTMP ;

    for ( $i=0; $i < $NTMP; $i++ ) {
        @wdlist  = split(/\s+/,$tmpList[$i]) ;
	$NAME    = $wdlist[0];
	$ARG     = $wdlist[1];
	$PAWPLOTMACRO_NAME[$i] = $NAME ;
	$PAWPLOTMACRO_ARG[$i]  = $ARG  ;
	print "  PAW > $NAME $ARG \n";

	# check that ARG is valid
	$ARGSTAT = -1 ;
	if ( $ARG eq "prefix"     ) { $ARGSTAT = 1; }
	if ( $ARG eq "prefixList" ) { $ARGSTAT = 1; }

	if ( $ARGSTAT < 0 ) {
	    $MSGERR[0] = "Invalid PAWPLOTMACRO ARG = '$ARG'" ;
	    $MSGERR[1] = "Valid ARGs : 'prefix' , 'prefixList' ";
	    sntools::FATAL_ERROR(@MSGERR) ;	    
	}
    }


} # end of parse_PAWPLOTMACRO


# ========================================
sub getGENVERSION(@) {

    # return GENVERSION from either the command-line $simArgs
    # or an auto-generated version name.
    # 
    # Dec 14, 2011: if repeatFlag is set then return previous GENVERSION

    my($itrain,$igen,$simArgs,$jobid) = @_;

    my ($iwd, $NWD, $wd, @tmp, $tmpVersion, @wdlist);
    my ($c3gen, $MODEL, $new_simArgs, $repeatFlag, $JOB ) ;

    # loop over the simArgs for a command-line override
    @wdlist     = split(/\s+/,$simArgs) ;
    $NWD        = scalar(@wdlist) ;
    $tmpVersion = "" ;

    $repeatFlag  = 
	$SIMTRAIN_REPEATFLAG[$itrain][$JOBID_GENSPEC] *
	$SIMTRAIN_REPEATFLAG[$itrain][$JOBID_GENPHOT] ;

    if ( $repeatFlag ) {
	return	$SIMTRAIN_GENVERSION[$itrain-1][$igen] ;
    }

    for ( $iwd = 0; $iwd < $NWD ; $iwd++ ) {
	$wd = $wdlist[$iwd] ;
	if ( $wd eq "GENVERSION" )  { $tmpVersion = $wdlist[$iwd+1] ; }
    }

    # if tmpVersion has not been specified by user, then
    # set the GENVERSION to be the name of the training
    # with an index appended.

    if ( $tmpVersion eq '' ) {
	$MODEL      = "$TRAIN_MODELNAME[$itrain]";
	$JOB        = "$JOBID_NAME[$jobid]";
	$c3gen      = sprintf("%3.3d", $igen);
	$tmpVersion = "$MODEL-$JOB-$c3gen" ;
	$new_simArgs = "$simArgs " . "GENVERSION $tmpVersion";
	# update global array to include GENVERSION arg.
	$TRAIN_CHANGE[$itrain][$jobid][$igen] = "$new_simArgs";
    }

    return $tmpVersion ;

} # end of getGENVERSION


# ===========================
sub  check_SSH_or_BATCH {

    # make sure that jobs are submitted via SSH or BATCH,
    # but not both. For BATCH, make sure that BATCH commands
    # are given for both TRAIN and FIT parts.

    my (@wdlist, @TMPNODES, $tmpLine, @tmp );

    # start with lots of error checking

    if ( $NNODE == 0 && $BATCH_FLAG == 0 ) {
	$MSGERR[0] = "Neither SSH or BATCH is defined." ;
	$MSGERR[1] = "Must declare either" ;
	$MSGERR[2] = "  NODELIST: <nodeList>" ;
	$MSGERR[3] = "or declare" ;
	$MSGERR[4] = "  BATCH_INFO_TRAIN: <command> <templateBatchFile>  <Ncore>" ;
	$MSGERR[5] = "  BATCH_INFO_FIT:   <command> <templateBatchFile>  <Ncore>" ;
	sntools::FATAL_ERROR(@MSGERR) ;	    	
    }

    if ( $NNODE && $BATCH_FLAG ) {
	$MSGERR[0] = "Cannot define both  SSH and BATCH." ;
	$MSGERR[1] = "Remove either NODELIST or BATCH_INFO key.";
	sntools::FATAL_ERROR(@MSGERR) ;	   
    }

    if ( $BATCH_FLAG == 3 || ($BATCH_FLAG == 2 && $BATCH_INFO ne "") ) {
	$MSGERR[0] = "WARNING: USING OBSOLETE BATCH_INFO: KEY!!!";
	$MSGERR[1] = "Change input to define BATCH information for both TRAIN and FIT";
	$MSGERR[2] = "by declaring both" ;
	$MSGERR[3] = "  BATCH_INFO_TRAIN: <command> <templateBatchFile>  <Ncore>" ;
	$MSGERR[4] = "  BATCH_INFO_FIT:   <command> <templateBatchFile>  <Ncore>" ;
	$MSGERR[5] = "keywords.";
	sntools::FATAL_ERROR(@MSGERR) ;	   
    }

    if ( $BATCH_FLAG == 1 && $BATCH_INFO eq ""  ) {
	$MSGERR[0] = "Must define BATCH information for both TRAIN and FIT";
	$MSGERR[1] = "by declaring both" ;
	$MSGERR[2] = "  BATCH_INFO_TRAIN: <command> <templateBatchFile>  <Ncore>" ;
	$MSGERR[3] = "  BATCH_INFO_FIT:   <command> <templateBatchFile>  <Ncore>" ;
	$MSGERR[4] = "keywords.";
	sntools::FATAL_ERROR(@MSGERR) ;	   
    }

    if ( $BATCH_FLAG == 1 && $BATCH_INFO ne ""  ) {
	print "\n !!!WARNING: USING OBSOLETE BATCH_INFO: KEY!!! \n";
	print "\n !!!WARNING: MAY CAUSE MEMORY CRASHES OR OTHER PROBLEMS!!! \n";
	print "  COPYING BATCH_INFO to BATCH_INFO_FIT and BATCH_INFO_TRAIN \n";
	$BATCH_INFO_TRAIN = $BATCH_INFO;
	$BATCH_INFO_FIT = $BATCH_INFO;
	##$MSGERR[0] = "WARNING: USING OBSOLETE BATCH_INFO: KEY!!!";
	##$MSGERR[1] = "Change input to define BATCH information for both TRAIN and FIT";
	##$MSGERR[2] = "by declaring both" ;
	##$MSGERR[3] = "  BATCH_INFO_TRAIN: <command> <templateBatchFile>  <Ncore>" ;
	##$MSGERR[4] = "  BATCH_INFO_FIT:   <command> <templateBatchFile>  <Ncore>" ;
	##$MSGERR[5] = "keywords.";
	##sntools::FATAL_ERROR(@MSGERR) ;	   
    }
    
    if ( scalar(@NODELIST) > 0  ) {

	@TMPNODES = @NODELIST ;
	@NODELIST = ();
	foreach $tmpLine ( @TMPNODES ) {
	    @tmp = split(/\s+/,$tmpLine) ;
	    @NODELIST =  ( @NODELIST , @tmp );
	}
	$NNODE    = scalar(@NODELIST);

	print "  Distribute jobs via SSH to \n\t @NODELIST \n" ;
	
     }
    else {
	print "  Distribute training batch jobs '$BATCH_INFO_TRAIN' \n";

	@wdlist = split(/\s+/,$BATCH_INFO_TRAIN ) ;
	$BATCH_COMMAND_TRAIN       = $wdlist[0];
	$BATCH_TEMPLATE_TRAIN = $wdlist[1];
	$BATCH_NODE_TRAIN = $wdlist[2];
	if ( !(-e $BATCH_TEMPLATE_TRAIN ) ) {
	    $MSGERR[0] = "Training Batch-template file $BATCH_TEMPLATE_TRAIN";
	    $MSGERR[1] = "does not exist.";
	    sntools::FATAL_ERROR(@MSGERR) ; 
	}

	print "  Distribute all other batch jobs (inc. fits) '$BATCH_INFO_FIT' \n";
	@wdlist = split(/\s+/,$BATCH_INFO_FIT ) ;
	$BATCH_COMMAND_FIT       = $wdlist[0];
	$BATCH_TEMPLATE_FIT = $wdlist[1];
	$BATCH_NODE_FIT = $wdlist[2];
	if ( !(-e $BATCH_TEMPLATE_FIT ) ) {
	    $MSGERR[0] = "Fitting Batch-template file $BATCH_TEMPLATE_FIT";
	    $MSGERR[1] = "does not exist.";
	    sntools::FATAL_ERROR(@MSGERR) ; 
	}
    }

} # end of check_SSH_or_BATCH 


# ================================
sub check_split_version {

    my ($VERSION) = @_ ;

    # abort if this version has already been split

    my ($fList, @tmp );

    if ( $DATAFLAG ) {
	$fList = "$SNDATA_ROOT/lcmerge/SPLIT*${VERSION}.*" ;
	@tmp   = qx(ls $fList 2>/dev/null );
	if ( scalar(@tmp) > 0 ) {
	    $MSGERR[0] = "Version $VERSION has already been split.";
	    $MSGERR[1] = "Need to remove" ;
	    $MSGERR[2] = "\t $fList";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

} # end of check_split_version

# ===================================
sub check_simFile {

    my ($simFile) = @_ ;
    my ($key, @tmp);

    # check that the GENPHOT simFile has required keys
    # for the training. ALthough this is an SNANA-input file,
    # there must be keys for the SALT2 training as well.

    print "\t Check simFile keys in:  $simFile ... \n";

    $key  = "SALT2_INSTRUMENT:" ;
    @tmp = sntools::parse_line($simFile, 1, $key, $OPTABORT );

    $key  = "SALT2_MAGSYS:" ;
    @tmp = sntools::parse_line($simFile, 1, $key, $OPTABORT );

    # Jan 4, 2012: make sure 'FORMAT_MASK: 2' is set for the
    # training sim-input file because only text format can
    # be used here. For the sim-data set FITS or TEXT is OK,
    # so only check the TRAINING-sim.


    # abort if obsolete key GEN_SNDATA_SIM is found.
    my (@tmp1, @tmp2, $keyfmt, @tmpFile, $incFile, $INCFLAG, $FORMAT_MASK ) ;

    $keyfmt  = "FORMAT_MASK:" ;
    $key     = "GEN_SNDATA_SIM:" ;
    @tmp1 = sntools::parse_line($simFile, 1, $key, $OPTQUIET );

    @tmpFile = sntools::parse_line($simFile, 1, $INCLUDE_KEY, $OPTQUIET);
    if ( scalar(@tmpFile) > 0 )  { 
	$incFile = "@tmpFile[0]" ; 
	$INCFLAG = 1; 
	@tmp2   = sntools::parse_line($incFile, 1, $key, $OPTQUIET );
    }
    else { @tmp2 = (); $INCFLAG = 0 ; }


    if ( scalar(@tmp1) > 0 ) {
	@MSGERR[0] = "Invalid format key '$key' in $simFile" ;
	@MSGERR[1] = "Replace obsolete '$key' with '$keyfmt'";
	sntools::FATAL_ERROR(@MSGERR);
    }
    if ( scalar(@tmp2) > 0 ) {
	@MSGERR[0] = "Invalid format key '$key' in $incFile" ;
	@MSGERR[1] = "Replace obsolete '$key' with '$keyfmt'";
	sntools::FATAL_ERROR(@MSGERR);
    }
    


    # ----------------------------------------------
    # now make sure that FORMAT_MASK is set to 2
    # (ie, don't allow default FORMAT_MASK=32 for FITS format)

    @tmp1 = sntools::parse_line($simFile, 1, $keyfmt, $OPTQUIET );

    if ( $INCFLAG ) 
    { @tmp2 = sntools::parse_line($incFile, 1, $keyfmt, $OPTQUIET ); }
    else 
    { @tmp2 = (); }


    $FORMAT_MASK = 0 ;
    if ( scalar(@tmp1) > 0 ) { $FORMAT_MASK = $tmp1[0] ; }
    if ( scalar(@tmp2) > 0 ) { $FORMAT_MASK = $tmp2[0] ; }


    if ( ($FORMAT_MASK & 2) == 0 ) {
	@MSGERR[0] = "Invalid FORMAT_MASK = $FORMAT_MASK";
	@MSGERR[1] = "Must set 'FORMAT_MASK: 2'  (TEXT format)";
	@MSGERR[2] = "in $simFile (or INCLUDE file)";
	sntools::FATAL_ERROR(@MSGERR);
    }

} # end of check_simFile

# ===================================
sub make_TrainBin_inputs {

  #######
  ### GOAL OF THIS STAGE: 
  ###
  ### FORCE snlc_sim.exe TO HANDLE
  ### BINARY ISSUES (aborts, regenerations)
  ### PRIOR TO PARALLEL BATCH/SSH SUBMISSION
  ### STAGE USES ONE LONG TRAIN_SIM.input FILE,
  ### TURNS OFF SPECSIM, SETS ALL NGEN_LC to 1
  #######

    my ($istage, $prefix, $sedcmd, $itrain, $igen, $NGEN, $jobid) ;
    my ($simFile, $simArgs, $specFile, $specArgs, $MODEL,  $path) ;
    my ($incdataFile, $incdataArgs);
    my ($repeatFlag, $incdataFlag);
    my ($BinArgs);

    $istage = $ISTAGE_TRAINBINARY ;
    $prefix = "$SCRIPT_PREFIX[$istage]" ;
    $path   = "$STAGE_PATH[$istage]";
    print " Prepare input file for script ${prefix}.pl \n";

    $INFILE_TRAINBINARY = "${path}/${prefix}.input" ;

    # append lines
    open  PTR_FILE , ">> $INFILE_TRAINBINARY" ;
    print PTR_FILE "# This input file created by script \n#\t $0 \n";
    print PTR_FILE "\n";
    print PTR_FILE "RUN_LCSIM:    T \n";
    print PTR_FILE "RUN_SPECSIM:  F \n";
    print PTR_FILE "\n";

    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {


	$MODEL   = $TRAIN_MODELNAME[$itrain] ;
	print PTR_FILE "# ----------------------------------------- \n";
	print PTR_FILE "# Args for training model: $MODEL \n";

	$incdataFlag = 0;
	$repeatFlag  = 
	    $SIMTRAIN_REPEATFLAG[$itrain][$JOBID_GENSPEC] *
	    $SIMTRAIN_REPEATFLAG[$itrain][$JOBID_GENPHOT] ;


	$NGEN  = $NINFILE_TRAIN[$itrain][$JOBID_GENPHOT] ;

	for ( $igen=1; $igen <= $NGEN; $igen++ ) {

	    $jobid   = $JOBID_GENPHOT ;  	    
	    $simFile = "$TRAIN_INFILE[$itrain][$jobid][$igen]" ;
	    $simArgs = "$TRAIN_CHANGE[$itrain][$jobid][$igen]" ;
	    $BinArgs = " NGENTOT_LC 0  NGEN_LC 1";

	    $jobid    = $JOBID_GENSPEC ;  
	    $specFile = "$TRAIN_INFILE[$itrain][$jobid][$igen]" ;
	    $specArgs = "$TRAIN_CHANGE[$itrain][$jobid][$igen]" ;

	    $simArgs   =~ s/\s+$// ;  # remove trailing blanks
	    $specArgs  =~ s/\s+$// ;  # remove trailing blanks
	    $BinArgs   =~ s/\s+$// ;  # remove trailing blanks
	    
	    if ( $repeatFlag ) { 
		print PTR_FILE "# Skip repeated $simFile \n";		
	    }
	    else {
		if ( $specArgs eq "" && $specFile eq "" ) {
		    my $ctmp = "SIMTRAIN_GENSPEC_INFILE[CHANGE]";
		    @MSGERR[0] = "Missing required $ctmp key for";
		    @MSGERR[1] = "  SIMTRAIN_MODEL = $MODEL";
		    @MSGERR[2] = "  SIMTRAIN_GENPHOT_SIMFILE = $simFile";
		    @MSGERR[3] = "  SIMTRAIN_GENPHOT_CHANGE  = $simArgs";
		    sntools::FATAL_ERROR(@MSGERR);
		}
		print PTR_FILE "LCSIM_INPUT:  $simFile $simArgs $BinArgs\n";
		unless ( $specFile eq "NONE" ) { ## ALLOW USER TO SKIP SPECSIM
		    print PTR_FILE "SPECSIM_ARGUMENTS: $specFile $specArgs \n";
		}
	    }
	    print PTR_FILE "END: \n\n";
	} # end GENPHOT NGEN loop

	$NGEN  = $NINFILE_TRAIN[$itrain][$JOBID_INCDATA] ;
	
	for ( $igen=1; $igen <= $NGEN; $igen++ ) {

	    $jobid = $JOBID_INCDATA ;
	    $incdataFile = "$TRAIN_INFILE[$itrain][$jobid][$igen]" ;
	    $incdataArgs = "$TRAIN_CHANGE[$itrain][$jobid][$igen]" ;

	    $incdataArgs =~ s/\s+$// ;  # remove trailing blanks

	    unless ( $incdataFile eq "NONE" ) { ## ALLOW USER TO SKIP INCDATA
		$incdataFlag++;
		print PTR_FILE "INCDAT_INPUT: $incdataFile \n";
		print PTR_FILE "INCDAT_ARGUMENTS: $incdataArgs \n";
	    }
	    if ( $incdataFlag > 0) {print PTR_FILE "END: \n\n";}
	} # end INCDATA NGEN loop

    }  # end itrain loop


    print PTR_FILE "\n";
    close PTR_FILE ;   

} # end of make_TrainBin_inputs

# ===================================
sub make_TrainGen_inputs {

    my ($istage, $prefix, $sedcmd, $itrain, $igen, $NGEN, $jobid) ;
    my ($simFile, $simArgs, $specFile, $specArgs, $MODEL,  $path) ;
    my ($incdataFile, $incdataArgs);
    my ($repeatFlag, $incdataFlag);
    my ($c2t);

    $istage = $ISTAGE_TRAINGEN ;
    $prefix = "$SCRIPT_PREFIX[$istage]" ;
    $path   = "$STAGE_PATH[$istage]";
    print " Prepare input files for script ${prefix}.pl \n";

    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {

	$c2t = sprintf("%2.2d", $itrain);
	$INFILE_TRAINGEN = "${path}/${prefix}_${c2t}.input" ;

	# append lines
	open  PTR_FILE , ">> $INFILE_TRAINGEN" ;
	print PTR_FILE "# This input file created by script \n#\t $0 \n";
	print PTR_FILE "\n";
	print PTR_FILE "RUN_LCSIM:    T \n";
	print PTR_FILE "RUN_SPECSIM:  T \n";
	print PTR_FILE "\n";

	$MODEL   = $TRAIN_MODELNAME[$itrain] ;
	print PTR_FILE "# ----------------------------------------- \n";
	print PTR_FILE "# Args for training model: $MODEL \n";

	$incdataFlag = 0;
	$repeatFlag  = 
	    $SIMTRAIN_REPEATFLAG[$itrain][$JOBID_GENSPEC] *
	    $SIMTRAIN_REPEATFLAG[$itrain][$JOBID_GENPHOT] ;


	$NGEN  = $NINFILE_TRAIN[$itrain][$JOBID_GENPHOT] ;

	for ( $igen=1; $igen <= $NGEN; $igen++ ) {

	    $jobid   = $JOBID_GENPHOT ;  	    
	    $simFile = "$TRAIN_INFILE[$itrain][$jobid][$igen]" ;
	    $simArgs = "$TRAIN_CHANGE[$itrain][$jobid][$igen]" ;

	    $jobid    = $JOBID_GENSPEC ;  
	    $specFile = "$TRAIN_INFILE[$itrain][$jobid][$igen]" ;
	    $specArgs = "$TRAIN_CHANGE[$itrain][$jobid][$igen]" ;

	    $simArgs   =~ s/\s+$// ;  # remove trailing blanks
	    $specArgs  =~ s/\s+$// ;  # remove trailing blanks
	    
	    if ( $repeatFlag ) { 
		print PTR_FILE "# Skip repeated $simFile \n";		
	    }
	    else {
		if ( $specArgs eq "" && $specFile eq "" ) {
		    my $ctmp = "SIMTRAIN_GENSPEC_INFILE[CHANGE]";
		    @MSGERR[0] = "Missing required $ctmp key for";
		    @MSGERR[1] = "  SIMTRAIN_MODEL = $MODEL";
		    @MSGERR[2] = "  SIMTRAIN_GENPHOT_SIMFILE = $simFile";
		    @MSGERR[3] = "  SIMTRAIN_GENPHOT_CHANGE  = $simArgs";
		    sntools::FATAL_ERROR(@MSGERR);
		}
		print PTR_FILE "LCSIM_INPUT:  $simFile $simArgs \n";
		unless ( $specFile eq "NONE" ) { ## ALLOW USER TO SKIP SPECSIM
		    print PTR_FILE "SPECSIM_ARGUMENTS: $specFile $specArgs \n";
		}
	    }
	    print PTR_FILE "END: \n\n";
	} # end GENPHOT NGEN loop

	$NGEN  = $NINFILE_TRAIN[$itrain][$JOBID_INCDATA] ;
	
	for ( $igen=1; $igen <= $NGEN; $igen++ ) {

	    $jobid = $JOBID_INCDATA ;
	    $incdataFile = "$TRAIN_INFILE[$itrain][$jobid][$igen]" ;
	    $incdataArgs = "$TRAIN_CHANGE[$itrain][$jobid][$igen]" ;

	    $incdataArgs =~ s/\s+$// ;  # remove trailing blanks

	    unless ( $incdataFile eq "NONE" ) { ## ALLOW USER TO SKIP INCDATA
		$incdataFlag++;
		print PTR_FILE "INCDAT_INPUT: $incdataFile \n";
		print PTR_FILE "INCDAT_ARGUMENTS: $incdataArgs \n";
	    }
	    if ( $incdataFlag > 0) {print PTR_FILE "END: \n\n";}
	} # end INCDATA NGEN loop

    print PTR_FILE "\n";
    close PTR_FILE ;   

    }  # end itrain loop

} # end of make_TrainGen_inputs

# ===================================
sub make_runTrain_script {

    # Create the input file $INFILE_TRAINRUN.
    # This is the input file to $SNANA_DIR/util/SALT2train_run.pl
    #
    # Sep 9, 2011 (RK) : write NODELIST or BATCH_INFO_TRAIN
    # Oct 31,2011 : use argument for RUN_PCAFIT and remove
    #               RUN_ERRORSNAKE and RUN_POSTCOLORLAWFIT keys
    #               Since the default is T in SALT2train_run.pl
    #
    # Nov 4, 2011: update to work for DATA or SIM
    #
    # Jan 30, 2012: paste in the [non-null] GLOBAL_TRAINRUN_KEYNAMES 
    #
    # Feb 04, 2012: allow for path in DATATRAIN_LISTFILE =>
    #               strip off name after last slapsh (/).
    #

    my ($istage, $prefix, $NGEN, $igen, $itrain, $sedcmd);
    my ($NGEN_INCDATA, $incdataFile);
    my ($LISTFILE, $listFile, $islash ) ;
    my ($simFile, $simArgs, $MODEL, $PCAFILE, $GENV, $path, $sdir, $link ) ;
    my ($NKEY, $ikey, $KEY, $VAL );
    my ($TRAINLISTFILE);
    my $OPT_SCRIPT = 0 ;
    
    $istage = $ISTAGE_TRAINRUN ;
    $prefix = "$SCRIPT_PREFIX[$istage]" ;
    $path   = "$STAGE_PATH[$istage]";
    $sdir   = "$STAGE_SDIR[$istage]";
    print " Prepare input file for script ${prefix}.pl \n";

    $INFILE_TRAINRUN = "${path}/${prefix}.input" ;

    # append lines
    open  PTR_FILE , ">> $INFILE_TRAINRUN" ;
    print PTR_FILE "\n# ------------------------------------------ \n";
    print PTR_FILE "# Input file created by script \n#\t $0 \n";
    print PTR_FILE "\n";

    print PTR_FILE "RUN_PCAFIT:  $NFIT_ITER_PCAFIT \n";
    print PTR_FILE "\n";

    if ( $NNODE ) 
    {  print PTR_FILE "NODELIST: @NODELIST  \n"; }
    else
    {  print PTR_FILE "BATCH_INFO: $BATCH_INFO_TRAIN \n"; }

    print PTR_FILE "OUTDIR:   $path/workSpace \n";

    # make quick link from topdir to workSpace
    $link =  "ln -s $sdir/workSpace  AAA_PCAFIT_TOPDIR";
    qx(cd $OUTPUT_TOPDIR ; $link );

    # paste in the GLOBAL_TRAINRUN_KEYNAMES
    print PTR_FILE "\n";
    $NKEY = scalar(@GLOBAL_TRAINRUN_KEYNAMES);
    for ( $ikey=0; $ikey < $NKEY; $ikey++ ) {
	$KEY = "$GLOBAL_TRAINRUN_KEYNAMES[$ikey]" ;
	$VAL = "$GLOBAL_TRAINRUN_KEYVALS[$ikey]" ;
	if ( $VAL ne "") { print  PTR_FILE "$KEY  $VAL \n" ; }
    }
    print PTR_FILE "\n";

    # prepare the DONE_STAMP
    my ($istage,$iscript,$scriptFile,$DONEFILE);

    $istage     = $ISTAGE_TRAINRUN ;  
    $iscript    = 0 ;
    $scriptFile = &scriptFileName($istage,$iscript,0,$OPT_SCRIPT);
    $DONEFILE   = "$OUTPUT_TOPDIR/${scriptFile}.DONE" ;

    print PTR_FILE  "DONE_STAMP: $DONEFILE \n" ;

    print PTR_FILE "\n" ;
    print PTR_FILE "# ============================================== \n";
    print PTR_FILE "\n";

    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {

	&make_pcafitConfig($itrain);

	$MODEL   = $TRAIN_MODELNAME[$itrain] ;
	$NGEN    = $NINFILE_TRAIN[$itrain][$JOBID_GENPHOT] ;
	$PCAFILE = $TRAIN_INFILE[$itrain][$JOBID_PCAFIT][2] ;
	$TRAINLISTFILE = $TRAIN_INFILE[$itrain][$JOBID_TRAINLIST][1] ;

	print PTR_FILE "CONFIG_FILE: $PCAFILE \n" ;

	if ( $SIMFLAG ) {
	    print PTR_FILE "GENVERSION: " ;

	    for ( $igen=1; $igen <= $NGEN; $igen++ ) {
		$GENV  = $SIMTRAIN_GENVERSION[$itrain][$igen] ;
		print PTR_FILE "$GENV ";
	    }
	    
	    # include INCDATA GENVERSIONS in trainrun input file
	    $NGEN_INCDATA    = $NINFILE_TRAIN[$itrain][$JOBID_INCDATA] ;
	    for ( $igen=1; $igen <= $NGEN_INCDATA; $igen++ ) {
		$incdataFile = "$TRAIN_INFILE[$itrain][$JOBID_INCDATA][$igen]" ;
		unless ( $incdataFile eq "NONE" ) {
		    $GENV = $SIMTRAIN_INCDATA_GENVERSION[$itrain][$igen] ;
		    print PTR_FILE "$GENV ";
		}
	    }
	    if ( $NGEN > 1 || $GENV ne $MODEL ) {
		print PTR_FILE "MODELNAME $MODEL \n";
	    }

	    # if not equal to "NONE", include TRAINLIST file in trainrun input file
	    unless ( $TRAINLISTFILE eq "NONE" ) {
		print PTR_FILE "TRAIN_LIST: $TRAINLISTFILE";
	    }
	    
	}
	else {
	    my ($PATH, $MM);
	    $PATH = $DATATRAIN_PATH[$itrain]  ;
	    $MM   = "MODELNAME  $MODEL";
	    print PTR_FILE "SALT2TRAININGPATH: $PATH   $MM\n";

	    $LISTFILE = $DATATRAIN_LISTFILE[$itrain];

	    # strip off name after last slash (/) to allow for path
	    $islash   = rindex($LISTFILE,"/");
	    $listFile = substr($LISTFILE,$islash+1,99); 
	    $listFile =~ s/\s+$// ;  # remove trailing blanks
	    print PTR_FILE "TRAIN_LIST:  $listFile \n" ;
	}



	print PTR_FILE "\n\n";
    }  # end itrain loop

    
    print PTR_FILE "\n";
    close PTR_FILE ;   


}  # end of make_runTrain_script


# ================================
sub  make_pcafitConfig {

    my ($itrain) = @_;

    my ($pcaFile_orig, $pcaFile_final, $ichange, $i, $key, $KEY, $MODEL ) ;
    my ($tmpLine, @wdlist, $iline, $c1, $NERR, $jobid, $NCHANGE) ;
    my ($IGNORE, $KEYMATCH, $TMPCHANGE, @CHANGELIST) ;
    my ($NLINE_UNDO, @LINE_UNDO, @pcaFile_contents,  $NLINE_pcaFile);
    my (@origLine,@newLine,@changeLine);

    $jobid = $JOBID_PCAFIT ;

    $pcaFile_orig  = $TRAIN_INFILE[$itrain][$jobid][1] ;
    $pcaFile_final = $TRAIN_INFILE[$itrain][$jobid][2] ;

    # break up the change-list into keys that start with @
    # so that each key can be replaces separately.
    # The input change-list has then all strung together
    # in the same way as for the sim-phot args.
    # 
    $TMPCHANGE  = "$TRAIN_CHANGE[$itrain][$jobid][1]" ;
    @CHANGELIST   = &PCAFIT_CHANGELIST($TMPCHANGE);
    $NCHANGE = scalar(@CHANGELIST) ;


#    print " ---------------------------------------------- \n";
    print "     Create pcafit-config file : $pcaFile_final \n";

    if ( -e $pcaFile_final) { qx(rm $pcaFile_final); }

    if ( $NCHANGE == 0 ) {
	qx(cp $pcaFile_orig $pcaFile_final);
	return ;
    }
    else {
	open  PTR_PCAFILE , "> $pcaFile_final" ; # open new file

# put the print statements back when we get the new pcafit
# that replaces the top part with keys
	print PTR_PCAFILE "# Original file '${pcaFile_orig}' modified by\n";
	print PTR_PCAFILE "# $0 \n";
	print PTR_PCAFILE "# \n";

	for ( $ichange=0; $ichange < $NCHANGE; $ichange++ )
	{ $USE_CHANGE[$ichange] = 0 ; }		
    }


    $NLINE_UNDO = 0;

    @pcaFile_contents  = qx(cat ${pcaFile_orig}) ;
    $NLINE_pcaFile     = scalar(@pcaFile_contents);

    $iline = 0;
    
    while ( $iline < $NLINE_pcaFile ) {

        $tmpLine  = "$pcaFile_contents[$iline]" ;
        @origLine = split(/\s+/,$tmpLine) ;    
	@newLine  = @origLine ;
	$KEY      = $origLine[0] ;
	$c1       = substr($KEY,0,1) ;

	# if we find a key, check for replacement
	if ( $c1 eq '@' ) { 
	    for ( $ichange=0; $ichange < $NCHANGE; $ichange++ ) {
		if ( $USE_CHANGE[$ichange] ) { next ; }
		$tmpLine     = "$CHANGELIST[$ichange]";
		@changeLine  = split(/\s+/,$tmpLine) ;
		$KEYMATCH    = ( $changeLine[0] eq "$KEY"   ) ;
		$IGNORE      = ( $changeLine[1] eq "IGNORE" ) ;

		if ( $KEYMATCH ) {	       
		    if ( $IGNORE == 0 ) {
			$NLINE_UNDO++ ;
			$LINE_UNDO[$NLINE_UNDO] = "REPLACED_@origLine";
			@newLine = @changeLine ;
		    }
		    $USE_CHANGE[$ichange] = 1 ;
		    goto PRINT_NEWLINE ;
		}
	    } # ichange loop
	} # c1 eq @

      PRINT_NEWLINE:
	print PTR_PCAFILE "@newLine \n";	
	$iline++ ;
    }

    # print UNDO lines at end of config file as a reminder;
    # these lines are ignored by the pcafit parsing code.
    print PTR_PCAFILE "\n\n" ;
    print PTR_PCAFILE "# -------------------------------- \n" ; 
    if ( $NLINE_UNDO ) { 
	print PTR_PCAFILE "# Replaced Keys: $NLINE_UNDO \n"; 
	for ( $i=1; $i <= $NLINE_UNDO ; $i++ ) 
	{ print PTR_PCAFILE "# $LINE_UNDO[$i] \n"; }	
    }
    else {
	print PTR_PCAFILE " Replaced keys: NONE \n"; 
    }

    close PTR_PCAFILE ;

    # abort if any pcafit-change was not used 

    $NERR = 0;
    for ( $ichange=0; $ichange < $NCHANGE; $ichange++ ) {
	if ( $USE_CHANGE[$ichange] == 0 ) {
	    $NERR++ ;
	    $tmpLine     = "$CHANGELIST[$ichange]";
	    @changeLine  = split(/\s+/,$tmpLine) ;
	    $MSGERR[$NERR] = "Could not substitute '@changeLine' \n";
	}
    }
    if ( $NERR ) { 
	$MODEL   = $TRAIN_MODELNAME[$itrain] ;
	$MSGERR[$NERR+1] = "Check SIMTRAIN_PCAFIT_CHANGE keys";
	$MSGERR[$NERR+2] = "for itrain=$itrain  ($MODEL).";
	sntools::FATAL_ERROR(@MSGERR) ;

    }
  
} #  end of  make_pcafitConfig


# ==========================================
sub PCAFIT_CHANGELIST {

    my ($changeString) = @_ ;

    
    # break up CHANGESTRING into strings that each start with @
    my (@changeList, $c1, $ichange, @wdlist, $wd ) ;

    if ( $changeString eq "" ) { return @changeList ; }

    @wdlist    = split(/\s+/,$changeString) ;
    $ichange   = -1 ;

    foreach $wd ( @wdlist  ) {       
	$c1  = substr($wd,0,1) ;

	if ( $c1 eq '@' ) { 
	    $ichange++ ; 
	    $changeList[$ichange] = "$wd";
	}	
	elsif ( $ichange >= 0 ) {
	    $changeList[$ichange] = "$changeList[$ichange] " . "$wd";
	}
	else {
	    $MSGERR[0] = "Invalid SIMTRAIN_PCAFIT_CHANGE key = '$wd'";
	    $MSGERR[1] = "pcafit key must start with '\@'";
	    sntools::FATAL_ERROR(@MSGERR);
	    
	}
    }

    return @changeList ;

} # end of PCAFIT_CHANGELIST 


# ========================================
sub set_STAGE_PATH {

    # Jan 18, 2012
    # for each $istage set  $STAGE_PATH and $STAGE_SDIR

    my ($istage, $name, $c2stage, $sDir, $path);

    for ( $istage = 1; $istage <= $NSTAGE ; $istage++ ) {
	$c2stage  = sprintf("%2.2d", $istage );
	$name     = $STAGE_NAME[$istage];
	$sDir     = "STAGE${c2stage}_${name}" ;
	$path     = "$OUTPUT_TOPDIR/${sDir}";
	$STAGE_PATH[$istage] = "$path" ; # full path 
	$STAGE_SDIR[$istage] = "$sDir" ; # subdir name
    }
}

# ==========================================
sub make_OUTDIRS {

    my ($c2stage, $istage, $name, $sDir, $path, $cmd, $err) ;
    my ($submitDir) ;

    # make output directories.
    $submitDir  = `pwd` ;
    $submitDir  =~ s/\s+$// ;  # remove trailing blanks


    # first make sure that we don't clobber our own working directory
    if ( $OUTPUT_TOPDIR eq $submitDir || $OUTPUT_TOPDIR eq './' ) {
	$MSGERR[0] = "OUTPUT_TOPDIR = $OUTPUT_TOPDIR" ;
	$MSGERR[1] = "Cannot set OUTPUT_TOPDIR to your current dir";
	$MSGERR[2] = "otherwise current dir will get clobbered !";
	$MSGERR[3] = "Set OUTPUT_TOPDIR to something else in ";
	$MSGERR[4] = "$MASTER_INFILE";
	sntools::FATAL_ERROR(@MSGERR);
    }
    
    print "\n Construct output directories in \n    $OUTPUT_TOPDIR \n";
    if ( -d $OUTPUT_TOPDIR ) { qx(rm -r $OUTPUT_TOPDIR) ; }

    for ( $istage = 1; $istage <= $NSTAGE ; $istage++ ) {

	$path = $STAGE_PATH[$istage]  ; # full path 
	$sDir = $STAGE_SDIR[$istage]  ; # subdir name

	$cmd = "mkdir -p $path" ;
	print "\t Create subdir $sDir \n" ;
	$err = system($cmd);

	if ( $err != 0 ) {
	    $MSGERR[0] = "Error=$err creating output directory";
	    $MSGERR[1] = "  $path";
	    $MSGERR[2] = "Check OUTPUT_TOPDIR in $MASTER_INFILE";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

    # Jun 2012: make sure this topdir can be found after
    #           remote login; i.e., that full path is specified
    #           rather than a subdir of current dir.
    sntools::checkDirExist($OUTPUT_TOPDIR, "OUTPUT_TOPDIR");

    print "\n";


} # end of make_OUTDIRS


# ===============================
sub make_fitnmlFiles {

    my($igroup, $istage) = @_;

    # loop over SIMGEN versions for this SALT2mu-group and use 
    # the associated nml-template file to create a separate 
    # 'split_and_fit' file for each GENVERSION.
    # The name of each nml file is $GENVERSION.nml
    #
    # Aug  8, 2011: add GZIP_FLAG key
    # Sep  9, 2011: use NODELIST or BATCH keys
    # Dec 29, 2011: glue GROUP + VERSION to make directory names
    #                ==> avoids conflict if same VERSION is used
    #                in different groups
    # Aug 28, 2012: allow no FITOPT: entries if SIMDATA_EXCLUDE_FITOPT = 1
    #
    # --------------------

    my ($group, $sedcmd, $sedMXLC, $sedITER1, $cmd, @tmp, $c2g);
    my ($KEY0, $KEY1, $KEY2, $KEY3, $KEY4, $KEY5);
    my ($tmpVar, $tmpVal, $nml);
    my ($PATH, $SDIR, $NVER, $iver) ;
    my ($FITDIR, $TEMP_FITDIR, $F3RES) ;
    my ($VERSION,$q,$IFITOPT,$IDTRAIN,$FITMODEL,$FITOPT_KEY);
    my ($FITNML_FILE_ORIG, $FITNML_FILE, $TEMP_FILE) ;
    my ($NTMP);
    my (@local_groupname, @local_nversion_group, @local_genversion);
    my (@local_exclude_fitopt, @local_infile_fitnml_list, @local_fitnml_list );
    my ($nfitopt, $c3, $ISIMOPT, $genversion, $fitoptkeytype);
    my ($whichcatlist, @fitdirlist, $tmpCAT1, $tmpCAT2);

    $nfitopt = $NMODEL_TRAIN; 
    #default, unless SIMDATA & SIMDATA_EXCLUDE_FITOPT[$igroup]==1
    $fitoptkeytype = 0;
    #default, unless $istage == $ISTAGE_DMFIXFIT AND $DMCORR_SIMTYPE[$igroup] != "AS-IS"

    if ( $istage == $ISTAGE_SIMFIT )  {
	$whichcatlist = 1;
	@local_groupname = @DATA_GROUPNAME ; 
	@local_nversion_group = @NVERSION_GROUP ;
	@local_infile_fitnml_list = @DATA_INFILE_FIT ;
	@local_fitnml_list = @DATA_FITNML ;
	if ( $SIMFLAG ) { 
	    @local_genversion = @SIMDATA_GENVERSION ;
	    @local_exclude_fitopt = @SIMDATA_EXCLUDE_FITOPT ;
	    if ( $local_exclude_fitopt[$igroup] ) { $nfitopt = 0; }
	} else {
	    @local_genversion = @DATAFIT_VERSION ;
	}
	    
    }
    if ( $istage == $ISTAGE_DMFIXFIT ) {
	$whichcatlist = 2;
	@local_groupname = @DMCORR_GROUPNAME ; 
	@local_nversion_group = @NVERSION_GROUP_DMCORR ;
	@local_exclude_fitopt = @DMCORR_EXCLUDE_FITOPT ;
	@local_infile_fitnml_list = @DMCORR_INFILE_FIT ;
	@local_fitnml_list = @DMCORR_FITNML ;
	@local_genversion = @DMCORR_GENVERSION ;
	if ( $local_exclude_fitopt[$igroup] ) { $nfitopt = 0; }
	unless ($DMCORR_SIMTYPE[$igroup] eq "AS-IS") {$fitoptkeytype = 1;}
    }

    $group = $local_groupname[$igroup];
    $c2g   = sprintf("%2.2d", $igroup );
    
	print "  Create SALT2-fit namelist file(s) for group '${group}'  \n";

	$NVER = $local_nversion_group[$igroup] ;

	@VERSION_LIST = split(/\s+/,$local_genversion[$igroup]) ;

	@INFILE_FITNML_LIST = split(/\s+/,$local_infile_fitnml_list[$igroup]) ;
	@FITNML_LIST        = split(/\s+/,$local_fitnml_list[$igroup]) ;

	$PATH = $STAGE_PATH[$istage] ; 
	$SDIR = $STAGE_SDIR[$istage] ; 

	for ( $iver = 0; $iver < $NVER; $iver++ ) {

	    $FITNML_FILE_ORIG  = $INFILE_FITNML_LIST[$iver]; 
		$FITNML_FILE       = $FITNML_LIST[$iver]; 
       	        $VERSION           = $VERSION_LIST[$iver];
	    $TEMP_FILE         = "TEMP_${VERSION}.nml" ; 
	    print "     $VERSION : $FITNML_FILE_ORIG -> $FITNML_FILE \n";

	    # create namelist input-file from template and remove existing
	    # split_and_fit keys
	    $sedcmd  = "sed " ;
	    $sedcmd  = "$sedcmd" . "-e '/OUTDIR/d' ";
	    $sedcmd  = "$sedcmd" . "-e '/VERSION/d' ";
	    $sedcmd  = "$sedcmd" . "-e '/NODELIST/d' ";
	    $sedcmd  = "$sedcmd" . "-e '/BATCH_INFO/d' ";
	    $sedcmd  = "$sedcmd" . "-e '/SNANA_MODELPATH/d' ";
	    $sedcmd  = "$sedcmd" . "-e '/OPT_SPLIT/d' ";
	    $sedcmd  = "$sedcmd" . "-e '/SIMTRAINFIT/d' ";

	    if ( -e $TEMP_FILE ) { qx(rm $TEMP_FILE) ; }
	    qx($sedcmd $FITNML_FILE_ORIG > $TEMP_FILE );

	    # put together new OUTPUT dir
	    $TEMP_FITDIR = "GROUP${c2g}_${VERSION}" ;

	    # create new split_and_fit keys
	    $KEY0 = "OUTDIR:     $PATH/${TEMP_FITDIR}" ;

	    if ( $NNODE ) 
	    {  $KEY1 = "NODELIST:   @NODELIST" ; }
	    else   
	    {  $KEY1 = "BATCH_INFO: $BATCH_INFO_FIT" ; }

	    $KEY2 = "VERSION:    $VERSION" ;
	    $KEY3 = "SNANA_MODELPATH:  $SNANA_MODELPATH" ;
	    $KEY4 = "GZIP_FLAG:  1"   ;
	    $KEY5 = "OPT_SPLIT:  $OPT_SPLIT";

	    # insert default &FITINP nml variables if not already defined.
	    $nml      = "SNLCINP" ;  # switch nml, Feb 8 2013
	    $tmpVar   = "MXLC_PLOT" ;
	    $tmpVal   = 5 ; # limit 5 plots per jobs (to limit disk usage)
	    $sedMXLC  = sntools::sed_nmlInsert($FITNML_FILE_ORIG, $nml, $tmpVar,$tmpVal);

	    $nml      = "FITINP" ;
	    $tmpVar   = "FUDGEALL_ITER1_MAXFRAC"; 
	    $tmpVal   = 0.01 ; # add .01*FLUXMAX err in quadrature on 1st iter
	    $sedITER1 = sntools::sed_nmlInsert($FITNML_FILE_ORIG, $nml, $tmpVar,$tmpVal);

	    # now paste in the split_and_fit stuff at the top
	    $q = '"' ;
	    $sedcmd  = "sed " ;
	    $sedcmd  = "$sedcmd" . "-e '1 i\\${KEY0}' " ;
	    $sedcmd  = "$sedcmd" . "-e '1 i\\${KEY1}' " ;
	    $sedcmd  = "$sedcmd" . "-e '1 i\\${KEY3}' " ;
	    $sedcmd  = "$sedcmd" . "-e '1 i\\${KEY4}' " ;
	    $sedcmd  = "$sedcmd" . "-e '1 i\\${KEY5}' " ;

	    if ( $SIMFLAG ) {
		$sedcmd  = "$sedcmd" . "$sedMXLC $sedITER1 ";
	    }

	    # tack on FITMODEL_NAME option for each training
	    $IFITOPT = 0;

	    # make list of VERSIONS as needed
	    
	    ### if $nfitopt == 0 => no versions needed. only fit with default
	    if ($nfitopt == 0){
		$IDTRAIN = 1;
		$sedcmd  = "$sedcmd" . "-e '1 i\\${KEY2}' " ;

		# make symbolic link to split_and_fit's temp dir so that it's 
		# easier to find the directories made by split_and_fit
		$FITDIR      = "${group}_${VERSION}" ;
		qx(cd $PATH ; ln -s $TEMP_FITDIR/${VERSION} $FITDIR ) ;
		$fitdirlist[$IDTRAIN] = $FITDIR; ## keep track of fitdirs for get_CATLISTinfo
		($tmpCAT1, $tmpCAT2) = get_CATLISTinfo(0,1,$whichcatlist, $igroup, $SDIR, $fitdirlist[1]);
		$NTMP = $tmpCAT2;
		$CATLIST_FITRES[$whichcatlist][$igroup][$IDTRAIN][$NTMP] = "$tmpCAT1";
		$NCATLIST_FITRES[$whichcatlist][$igroup][$IDTRAIN] = $tmpCAT2;
	    } 
	    ### do fits with multiple models
	    else {
		# if DMCORR_SIMTYPE == "AS-IS" or if $istage == $ISTAGE_SIMFIT 
		# put in single photometry VERSION 
		if ($fitoptkeytype == 0) {
		    $IDTRAIN=1;
		    $sedcmd  = "$sedcmd" . "-e '1 i\\${KEY2}' " ;

		    # make symbolic link to split_and_fit's temp dir so that it's 
		    # easier to find the directories made by split_and_fit
		    $FITDIR      = "${group}_${VERSION}" ;
		    qx(cd $PATH ; ln -s $TEMP_FITDIR/${VERSION} $FITDIR ) ;
		    $fitdirlist[$IDTRAIN] = $FITDIR; ## keep track of fitdirs for get_CATLISTinfo
		} else {
		    # finish up simtrain fits by making correct symbolic links
		    $ISIMOPT = $NMODEL_TRAIN;
		    for ($IDTRAIN=1; $IDTRAIN <= $ISIMOPT ; $IDTRAIN++ ) {
			$c3 = $TRAIN_CHARINDEX[$IDTRAIN];
			$genversion = "${VERSION}_${c3}";

			# make symbolic link to split_and_fit's temp dir so that it's 
			# easier to find the directories made by split_and_fit
			$FITDIR      = "${group}_${VERSION}_${c3}" ;
			qx(cd $PATH ; ln -s $TEMP_FITDIR/${VERSION}_${c3} $FITDIR ) ;
			$fitdirlist[$IDTRAIN] = $FITDIR;
		    }
		}
		for ( $IDTRAIN=1; $IDTRAIN <= $nfitopt ; $IDTRAIN++ ) {
		    $FITOPT_KEY = get_FITOPTKEYline($fitoptkeytype, $IDTRAIN, $VERSION);
		    $sedcmd      = "$sedcmd" . "-e '1 i\\${FITOPT_KEY}' " ;

		    # keep track of fitres files to be catenated later
		    # just before running  SALT2mu
		    $IFITOPT++ ;
		    if ($fitoptkeytype == 0) {
			# if DMCORR_SIMTYPE == "AS-IS" or if $istage == $ISTAGE_SIMFIT 
			## ALL FITRES FILES IN SAME FITDIR -> replace $fitdirlist[$IDTRAIN] with $fitdirlist[1]
			($tmpCAT1, $tmpCAT2) = get_CATLISTinfo($IFITOPT, $IDTRAIN,  $whichcatlist, $igroup, $SDIR, $fitdirlist[1]);
		    } else {
			## ALL FITRES FILES ARE FITOPT000.FITRES -> replace $IFITOPT with 0
			($tmpCAT1, $tmpCAT2) = get_CATLISTinfo(0, $IDTRAIN, $whichcatlist, $igroup, $SDIR, $fitdirlist[$IDTRAIN]);
		    }
		    $NTMP = $tmpCAT2;
		    $CATLIST_FITRES[$whichcatlist][$igroup][$IDTRAIN][$NTMP] = "$tmpCAT1" ;
		    $NCATLIST_FITRES[$whichcatlist][$igroup][$IDTRAIN] = "$tmpCAT2" ;
		}
	    }

	$cmd = "cat $TEMP_FILE | $sedcmd > $FITNML_FILE ";
#	print " xxx $cmd \n";
	qx($cmd ; rm $TEMP_FILE );

	} # igen

} # end of make_fitnmlFiles



# ================================
sub copy_inputFiles {


    my ($PATH, $istage, $itrain ) ;   
    my ($cdd, $cmdlink);
    my $nfitopt;

    print "\n  Copy input files to: \n";

    # set $nfitopt default
    $nfitopt = $NMODEL_TRAIN;

    # copy master file
    qx(cp $MASTER_INFILE $OUTPUT_TOPDIR);

    # make symbolic link to easy-to-find AAA_MASTER.INPUT
    $cdd     = "cd $OUTPUT_TOPDIR"  ;
    $cmdlink = "ln -s $MASTER_INFILE AAA_MASTER.INPUT" ;
    qx($cdd ; $cmdlink );

    # write user-COMMENT fields into separate file
    &make_readmeFile();


    # --------------------------------------

    if ( $SIMFLAG ) {
	# copy simgen input files (along with INCLUDE files)
	&copy_inputFiles_TRAINGEN();

	# copy translate-input file (same as TRAIN_SIM input file)	
	&copy_inputFiles_TRANSLATE();

	# copy incdata-input file (same as TRAIN_SIM input file)	
	&copy_inputFiles_INCDATA();
    }

    # -----------------------------------------
    # copy training-run file

    &copy_inputFiles_TRAINRUN();


    # -------------------------------------
    if ( $SIMFLAG ) {
	&copy_inputFiles_SIMGEN($ISTAGE_DATABINARY);
	&copy_inputFiles_SIMGEN($ISTAGE_SIMGEN);
    }

    # -------------------------------------

    &copy_inputFiles_FIT();

    # ---------------------------
    $istage = $ISTAGE_SALT2mu ; 
    print "\t $STAGE_SDIR[$istage] \n";
    $PATH = $STAGE_PATH[$istage] ;

    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ )  {
	if ( $SIMFLAG ) {
	    if ($SIMDATA_EXCLUDE_FITOPT[$igroup]) {
		$nfitopt=1;
	    } else {
		$nfitopt = $NMODEL_TRAIN;
	    }
	}
	for ( $itrain = 1; $itrain <= $nfitopt; $itrain++ ) {
	    &copy_SALT2mu_input($itrain,$igroup,$istage); 
	}
    }

    # -------------------------------------
    #
    # Run copy_inputFiles_COSMO even if no input file
    # This subroutine makes the data group subdirectories
    # in addition to copying the input files
    # 
    if ( $OPT_COSMO ) {
	$istage = $ISTAGE_COSMO ;
	print "\t $STAGE_SDIR[$istage] \n";
	$PATH = $STAGE_PATH[$istage];

    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ )  {
	if ( $SIMFLAG ) {
	    if ($SIMDATA_EXCLUDE_FITOPT[$igroup]) {
		$nfitopt=1;
	    } else {
		$nfitopt = $NMODEL_TRAIN;
	    }
	}
	for ( $itrain = 1; $itrain <= $nfitopt; $itrain++ ) {
	    &copy_inputFiles_COSMO($itrain,$igroup,$istage); 
	}
    }
    }

    # -------------------------------------
    #
    # Run these copy subroutines if DMCORRFLAG is 1
    # This subroutine makes the data group subdirectories as needed
    # in addition to copying the input files
    # 

    if ( $DMCORRFLAG ) {
	#copy sim-gen files
	&copy_inputFiles_DMFIXSIM();
	#copy fit files
	&copy_inputFiles_DMFIXFIT();
	#
	#make/copy files for DMFIXPREP, DMFIXRUN, and DMFIXCOSMO
	#
	#$PATH = $STAGE_PATH[$istage] ;
	for ( $igroup = 1; $igroup <= $NGROUP_DMCORR; $igroup++ )  {
	    if ( $SIMFLAG ) {
		if ($DMCORR_EXCLUDE_FITOPT[$igroup]) {
		    $nfitopt=1;
		} else {
		    $nfitopt = $NMODEL_TRAIN;
		}
	    }
	    for ( $itrain = 1; $itrain <= $nfitopt; $itrain++ ) {
		&copy_SALT2mu_input($itrain,$igroup,$ISTAGE_DMFIXPREP); 
		&copy_inputFiles_DMFIXRUN($itrain,$igroup,$ISTAGE_DMFIXRUN);
		&copy_inputFiles_COSMO($itrain,$igroup,$ISTAGE_DMFIXCOSMO);
	    }
	}
	#
    }
	

}   # end of copy_inputFiles


# =================================
sub copy_inputFiles_TRAINGEN {


    ####
    # Despite the implications of its name, 
    # this subroutine copies TRAINGEN input
    # files to all stages needing these files.
    # Currently includes $ISTAGE_TRAINBINARY,
    # $ISTAGE_TRAINGEN, and $ISTAGE_TRANSLATE.
    # To copy files to additional directories,
    # edit subroutine copy_simFile_train(@).
    # (JLM)
    ####

    my ($istage) = $ISTAGE_TRAINGEN;

    my ($PATH, $jobid, $igen, $simFile, $itrain, $NGEN, $igen );
    my ($repeatFlag);

    $PATH   = $STAGE_PATH[$istage] ; 
    print "\t $STAGE_SDIR[$istage] \n";


    $jobid = $JOBID_GENPHOT ;

    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {
	$NGEN =	$NINFILE_TRAIN[$itrain][$jobid] ;
	for ( $igen=1; $igen <= $NGEN; $igen++ ) {

	    $repeatFlag  =  $SIMTRAIN_REPEATFLAG[$itrain][$jobid] ;
	    if ( $repeatFlag ) { next ; }

	    $simFile = "$TRAIN_INFILE[$itrain][$jobid][$igen]" ;
	    if ( ! ( -e "$PATH/$simFile" ) ) {	    
		# check option to substitute GLOBAL_GENMODEL
		&copy_simFile_train($simFile);
	    }
	}
    }

} # end of copy_inputFiles_TRAINGEN

# ===============================
sub copy_inputFiles_TRANSLATE {
    my ($istage, $PATH, $cmd);
    $istage = $ISTAGE_TRANSLATE ; 
    $PATH   = $STAGE_PATH[$istage] ; 
    print "\t $STAGE_SDIR[$istage] \n";
    $cmd = "cp $INFILE_TRAINBINARY $PATH/";
    qx($cmd);
}

# ===============================
sub copy_inputFiles_INCDATA {
    my ($istage, $PATH, $cmd);
    $istage = $ISTAGE_INCDATA ; 
    $PATH   = $STAGE_PATH[$istage] ; 
    print "\t $STAGE_SDIR[$istage] \n";
    $cmd = "cp $INFILE_TRAINBINARY $PATH/";
    qx($cmd);
}


# =================================
sub copy_inputFiles_TRAINRUN {

    my ($istage, $PATH, $pcaFile, $listFile, $itrain );

    # move the TEMP_INFILE_TRAINRUN to outdir location
    $istage = $ISTAGE_TRAINRUN ; 
    $PATH = $STAGE_PATH[$istage] ;   
    print "\t $STAGE_SDIR[$istage] \n";

    # copy training list and pcafit config file for each training
    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {

	if ( $DATAFLAG ) {
	    $listFile = $DATATRAIN_LISTFILE[$itrain];
	    if ( !(-e $listFile) ) {
		$MSGERR[0] = "Cannot find trainling list-file '${listFile}'" ;
		sntools::FATAL_ERROR(@MSGERR);
	    }
	    qx(cp $listFile $PATH);
	}

	$pcaFile = $TRAIN_INFILE[$itrain][$JOBID_PCAFIT][2];
	if ( !(-e $pcaFile) ) {
	    $MSGERR[0] = "Cannot find pcafit-config file '${pcaFile}'" ;
	    sntools::FATAL_ERROR(@MSGERR);
	}
	qx(mv $pcaFile $PATH );
    }
    # copy global config file if it exists.
    if ( -e $CONFIG_FILE ) { qx(cp $CONFIG_FILE $PATH ); }

    # copy BATCH template file for batch mode.
    if ( $BATCH_FLAG ) { qx(cp $BATCH_TEMPLATE_TRAIN $PATH) ; }

} # end of copy_inputFiles_TRAINRUN


# ===============================
sub copy_inputFiles_SIMGEN(@) {

    my ($istage) = (@_);
    # copy SIMGEN files for data set (after training);

    my ($PATH, @simFileList, @incFile, $simFile, $igroup );

    print "\t $STAGE_SDIR[$istage] \n";
    $PATH = $STAGE_PATH[$istage] ; 

    if ( $BATCH_FLAG) { qx(cp $BATCH_TEMPLATE_FIT $PATH) ; }

    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) {
	
        @simFileList  = split(/\s+/,$SIMDATA_INFILE_GEN[$igroup]) ;
	foreach $simFile ( @simFileList ) {
	    qx(cp $simFile $PATH);
	    # check for include file
	    @incFile = sntools::parse_line($simFile,1,$INCLUDE_KEY, $OPTQUIET);
	    if ( scalar(@incFile) > 0 ) { qx(cp $incFile[0] $PATH ); }
	}
    }

} # end of copy_inputFiles_SIMGEN


# ===============================
sub copy_inputFiles_DMFIXSIM {

    # copy SIMGEN files for DMFIXSIM ;
    # this version also makes SUBDIRS for DMCORR groups. 

    my ($istage, $PATH, @simFileList, @incFile, $simFile, $igroup );
    my ($group, $topDir);

    $istage = $ISTAGE_DMFIXSIM ; 
    print "\t $STAGE_SDIR[$istage] \n";
    $PATH = $STAGE_PATH[$istage] ; 
    for ( $igroup = 1; $igroup <= $NGROUP_DMCORR; $igroup++ ) {
	
	$topDir = $PATH;

        @simFileList  = split(/\s+/,$DMCORR_INFILE_GEN[$igroup]) ;
	foreach $simFile ( @simFileList ) {
	    qx(cp $simFile $topDir);
	    # check for include file
	    @incFile = sntools::parse_line($simFile,1,$INCLUDE_KEY, $OPTQUIET);
	    if ( scalar(@incFile) > 0 ) { qx(cp $incFile[0] $topDir ); }
	}
    }
    if ( $BATCH_FLAG) { qx(cp $BATCH_TEMPLATE_FIT $topDir) ; }

} # end of copy_inputFiles_DMFIXSIM



# ==========================
sub copy_inputFiles_FIT {

    # use SIMFIT file as a template to create a fit-namelist file
    # for each GENVERSION.

    my ($istage, $PATH, $fitFile);

    $istage = $ISTAGE_SIMFIT ; 
    print "\t $STAGE_SDIR[$istage] \n";
    $PATH = $STAGE_PATH[$istage] ;
    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) {       
	@FITNML_LIST        = split(/\s+/,$DATA_FITNML[$igroup]) ;
	foreach $fitFile ( @FITNML_LIST ) {
	    if ( -e $fitFile ) {  qx(mv $fitFile $PATH); }
	}
    }
    if ( $BATCH_FLAG ) { qx(cp $BATCH_TEMPLATE_FIT $PATH/) ; }

    
} # end of copy_inputFiles_FIT

# ==========================
sub copy_inputFiles_DMFIXFIT {

    # use SIMFIT file as a template to create a fit-namelist file
    # for each GENVERSION.
    # this version also makes SUBDIRS for DMCORR groups.

    my ($istage, $PATH, $fitFile);
    my ($group, $topDir);

    print "Copying DMFIXFIT files ... \n";
    $istage = $ISTAGE_DMFIXFIT ; 
    print "\t $STAGE_SDIR[$istage] \n";
    $PATH = $STAGE_PATH[$istage] ;
    for ( $igroup = 1; $igroup <= $NGROUP_DMCORR; $igroup++ ) {       

	# make group subdir
	#$group = $DMCORR_GROUPNAME[$igroup];
	#$topDir = "$PATH/$group" ;
	$topDir = $PATH;
	#qx(mkdir $topDir);

	@FITNML_LIST        = split(/\s+/,$DMCORR_FITNML[$igroup]) ;
	foreach $fitFile ( @FITNML_LIST ) {
	    if ( -e $fitFile ) {  
              qx(mv $fitFile $topDir); 
	    }
	}
    }
    if ( $BATCH_FLAG ) { qx(cp $BATCH_TEMPLATE_FIT $topDir/) ; }

    
} # end of copy_inputFiles_DMFIXFIT

# ============================================
sub copy_SALT2mu_input {

    my($itrain,$igroup,$istage) = @_;

    # replace input filename and prefix key in SALT2mu-input file,
    # then move it to the official area.
    # make array of SALT2mu-fitres file names (with path) 
    # for later use by get_COSMO_commands

    my ($topDir, $group) ;
    my ($MODEL, $fresFile_in , $prefix, $inFile, $sedcmd, $cmd) ;
    my ($fresFile_out);

    my (@local_groupname);

    if ($istage == $ISTAGE_SALT2mu) {
	@local_groupname = @DATA_GROUPNAME;
    } elsif ($istage == $ISTAGE_DMFIXPREP) {
	@local_groupname = @DMCORR_GROUPNAME;
    } else {
	$MSGERR[0] = "In subroutine copy_SALT2mu_input ";
	$MSGERR[1] = "Input stage $istage is not recognized ";
	$MSGERR[2] = "Allowed stages are ISTAGE_SALT2mu $ISTAGE_SALT2mu and ";
	$MSGERR[3] = " ISTAGE_DMFIXPREP $ISTAGE_DMFIXPREP ";
	$MSGERR[4] = "Check code and try again. ";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $MODEL       = "$TRAIN_MODELNAME[$itrain]" ;
    $fresFile_in = "${MODEL}.FITRES" ;
    $prefix      = &prefixSALT2mu($itrain, $istage);
    $inFile      = "${prefix}.input" ;   # SALT2mu input file
    $fresFile_out = "${prefix}.fitres" ; # SALT2mu fitres file

    $group    = $local_groupname[$igroup] ;
    $topDir   = "$STAGE_PATH[$istage]/$group" ;
    $CPLIST_FITRES[$istage][$igroup][$itrain] = "${topDir}/${fresFile_out}"; 

    if ( !(-e $INFILE_SALT2mu) ) {
	$MSGERR[0] = "Cannot find SALT2mu file '${INFILE_SALT2mu}'" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    # create subdir for this SALT2mu-group
    if ( $itrain == 1 )  { qx(mkdir $topDir); }

    # copy/modify SALT2mu input file into proper directory.
    # replace file= with correct input-fitres file, and set
    # appropriate 'prefix' for output filenames
    $sedcmd = "sed" ;
    $sedcmd = "$sedcmd " .  "-e '1 i\\file=${fresFile_in}'" ;
    $sedcmd = "$sedcmd " .  "-e '1 i\\prefix=${prefix}'" ;
    $sedcmd = "$sedcmd " .  "-e '/file=/d'" ;
    $cmd    = "$sedcmd $INFILE_SALT2mu > $topDir/$inFile" ;
    qx($cmd) ;

} # end of copy_SALT2mu_input

# ============================================
sub copy_inputFiles_DMFIXRUN(@) {

    my( $itrain, $igroup, $istage) = @_;

    # replace input filename and prefix key in DMFIXRUN-input file,
    # then move it to the official area.

    my ($topDir, $group) ;
    my ($MODEL, $fresFile_in , $prefix, $inFile, $sedcmd, $cmd) ;
    my (@local_groupname);

	@local_groupname = @DMCORR_GROUPNAME;

    $MODEL       = "$TRAIN_MODELNAME[$itrain]" ;
    $prefix      = &prefixSALT2mu($itrain, $ISTAGE_DMFIXPREP);
    $fresFile_in = "${prefix}.FITRES" ;
    $inFile      = "DMFIXRUN.${prefix}.input" ;   

    $group    = $local_groupname[$igroup] ;
    $topDir   = "$STAGE_PATH[$istage]/$group" ;

    # create subdir for this DMFIXRUN-group
    if ( $itrain == 1 )  { qx(mkdir $topDir); }

} # end of copy_inputFiles_DMFIXRUN

# ============================================
sub copy_inputFiles_COSMO(@) {

    my( $itrain, $igroup, $istage ) = @_;

    # replace input filename and prefix key in COSMO-input file,
    # then move it to the official area.

    my ($topDir, $group) ;
    my ($MODEL, $fresFile_in , $prefix, $inFile, $sedcmd, $cmd) ;
    my (@local_groupname);

    if ($istage == $ISTAGE_COSMO ) {
	@local_groupname = @DATA_GROUPNAME;
    } elsif ($istage == $ISTAGE_DMFIXCOSMO) {
	@local_groupname = @DMCORR_GROUPNAME;
    }

    $MODEL       = "$TRAIN_MODELNAME[$itrain]" ;
    $prefix      = &prefixSALT2mu($itrain, $istage);
    $fresFile_in = "${prefix}.FITRES" ;
    $inFile      = "COSMO.${prefix}.input" ;   # COSMO input file

    $group    = $local_groupname[$igroup] ;
    $topDir   = "$STAGE_PATH[$istage]/$group" ;

    if ( $COSMO_INFILE ne "" && !(-e $COSMO_INFILE) ) {
        $MSGERR[0] = "Cannot find COSMO file '${COSMO_INFILE}'" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    # create subdir for this COSMO-group
    if ( $itrain == 1 )  { qx(mkdir $topDir); }

    if ( $COSMO_INFILE ne "" ){
    # CURRENTLY SET UP FOR sncosmo_mcmc.exe INPUT FILE FORMAT
    # copy/modify COSMO input file into proper directory.
    # replace FITRES_FILE: with correct inpu-fitres file
	$sedcmd = "sed" ;
	$sedcmd = "$sedcmd " .  "-e '1 i\\FITRES_FILE: ${fresFile_in}'" ;
	$sedcmd = "$sedcmd " .  "-e '/FITRES\\_FILE:/d'" ;
	$cmd    = "$sedcmd $COSMO_INFILE > $topDir/$inFile" ;
	qx($cmd) ;
    }

} # end of copy_inputFiles_COSMO

# ============================================
sub make_readmeFile {

    my ($igroup, $itrain, $NGEN, $igen, $cdate, $cmodel, $comment) ;
    my ($README_FILE,  $FITMODEL );

    $cdate = `date` ;

    $README_FILE = "${OUTPUT_TOPDIR}/AAA_README.TXT" ;
    open  PTR_README  , ">> $README_FILE" ; 

    print PTR_README "     USER COMMENTS  \n" ; 
    print PTR_README "    ~~~~~~~~~~~~~~~~ \n" ; 


    foreach $comment ( @USER_COMMENTS ) 
    { print PTR_README " $comment \n" ; }

    # ====================================================

    print PTR_README "\n\n" ;

    print PTR_README "     SCRIPT-GENERATED COMMENTS  \n" ; 
    print PTR_README "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n" ; 

    
    print PTR_README " Full name of master script: \n\t $0 \n" ;
    print PTR_README " SNANA_DIR  : $SNANA_DIR \n" ;
    print PTR_README " Start Time : $cdate \n" ;
    print PTR_README 
	"\n Summary of models to train with pcafit: \n" ; 
    print PTR_README
	"   Trained  model                  GENVERSION(s) \n";
    print PTR_README
	"  ----------------------------------------------- \n";

    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {
	$FITMODEL = "SALT2.$TRAIN_MODELNAME[$itrain]" ;
	$cmodel    = sprintf("%-30s", $FITMODEL);


#	$NGEN =	$NINFILE_TRAIN[$itrain][$jobid] ;
	print PTR_README " $cmodel ";
#	for ( $igen=1; $igen <= $NGEN; $igen++ ) {
#	    print PTR_README "$SIMTRAIN_GENVERSION[$itrain][$igen] ";
#	}
	print PTR_README "\n";
    }


    print PTR_README 
	"\n Summary of Simulations to Fit with each training above. \n" ; 

    print PTR_README 
	"    Group                 GENVERSION(s) \n" ;
    print PTR_README 
	"  -------------------------------------- \n" ;
    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) {
        my $group = $DATA_GROUPNAME[$igroup] ;
	my $cgroup = sprintf("%-20s", $group );
	print PTR_README "   $cgroup   @VERSION_LIST \n" ;
    }

    print PTR_README "\n\n End of README file. \n";

    close PTR_README ;

} # end of make_readmeFile

# ============================================
sub copy_simFile_train($simFile) {

    my($simFile) = @_;

    my ($PATH1, $PATH2, $PATH3, $cmd, $cmd1, $cmd2, $cmd3, $q ) ;
    my ($incFile,@tmpFile, $copyFile, @GENMODEL_ORIG );
    my $INCFILE_FLAG = 0;
    my @COPYLIST = ( $simFile ) ;

    $PATH1 = $STAGE_PATH[$ISTAGE_TRAINBINARY] ; 
    $PATH2 = $STAGE_PATH[$ISTAGE_TRAINGEN] ; 
    $PATH3 = $STAGE_PATH[$ISTAGE_TRANSLATE] ; 

    # check for INCLUDE file
    @tmpFile = sntools::parse_line($simFile,1,$INCLUDE_KEY, $OPTQUIET);
    if ( scalar(@tmpFile) > 0 )  { 
	$incFile = "@tmpFile[0]" ; 
	@COPYLIST = ( @COPYLIST , $incFile );
	$INCFILE_FLAG = 1 ;

	if ( !(-e $incFile ) ) {
	    $MSGERR[0] = "Cannot find INCLUDE file '${incFile}'" ;
	    $MSGERR[1] = "for SIMGEN-train input file '${simFile}'" ;
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

    # if there is no global GENMODEL then just copy sim-input file(s)

    if ( $GLOBAL_GENMODEL eq ' ' ) {
	foreach $copyFile ( @COPYLIST )  {
	    $cmd1 = "cp $copyFile $PATH1 ";
	    $cmd2 = "cp $copyFile $PATH2 ";
	    $cmd3 = "cp $copyFile $PATH3 ";
	    qx($cmd1 ; $cmd2 ; $cmd3) ;
	}
	return ;
    }


    # if we get here then substitute GENMODEL: argument
    # with GLOBAL_GENMODEL ... make sure to check both
    # the simFile and the INCLUDE file

    $q = '"' ;

    my ($INSERT1, $INSERT2, $INSERT3, $sedcmd );


    foreach $copyFile ( @COPYLIST )  {
	$cmd = "grep  ${q}GENMODEL:${q}  $copyFile"  ;
	@GENMODEL_ORIG = qx { $cmd };

	if ( scalar(@GENMODEL_ORIG) > 0 ) {

	    my @tmpList  = split(/\s+/,$GENMODEL_ORIG[0]) ;
	    my $tmp1 = "original GENMODEL $tmpList[1]" ;
	    my $tmp2 = "GLOBAL_GENMODEL = $GLOBAL_GENMODEL" ;

	    $INSERT1 = "# Replace $tmp1" ;
	    $INSERT2 = "# with $tmp2" ;
	    $INSERT3 = "GENMODEL: $GLOBAL_GENMODEL";

	    $sedcmd  = "sed" ;
	    $sedcmd  = "$sedcmd " . "-e '1 i\\$INSERT1'" ;
	    $sedcmd  = "$sedcmd " . "-e '1 i\\$INSERT2'" ;
	    $sedcmd  = "$sedcmd " . "-e '1 i\\$INSERT3'" ;
	    $sedcmd  = "$sedcmd " . "-e '/GENMODEL:/d' ";
	    $cmd1    = "$sedcmd $copyFile > $PATH1/$copyFile";
	    $cmd2    = "$sedcmd $copyFile > $PATH2/$copyFile";
	    $cmd3    = "$sedcmd $copyFile > $PATH3/$copyFile";
	}
	else {
	    $cmd1 = "cp $copyFile $PATH1 "; 
	    $cmd2 = "cp $copyFile $PATH2 "; 
	    $cmd3 = "cp $copyFile $PATH3 "; 
	}
	
	# do the copy or sed command
	qx($cmd1 ; $cmd2 ; $cmd3 ) ;

    }

} # end  of copy_simFile_train


# ============================================
sub make_RUNscripts {

    my ($prefix, $lastPrefix, $scriptFile, $istage, $igroup, $iscript) ;
    my (@CMD, $cmdx );
    my $OPT_SCRIPT = 1;
    print "\n  Create stage-scripts: \n ";
    

    # TRAINGEN & TRANSLATE & INCDATA --------------------------------
    if ( $SIMFLAG ) {
	$istage     = $ISTAGE_TRAINBINARY ; $iscript = 0; @CMD = ( );
	$scriptFile = &scriptFileName($istage,$iscript,0,$OPT_SCRIPT );
	$prefix     = "$SCRIPT_PREFIX[$istage]" ;
	$CMD[0]     = "${prefix}.pl ${prefix}.input > ${prefix}.log" ;
	&update_script($scriptFile, 1, @CMD );

	make_TRAINGEN_scripts($ISTAGE_TRAINGEN);

	$lastPrefix = $prefix ;
    
	$istage     = $ISTAGE_TRANSLATE ;  $iscript = 0;  @CMD = ( );
	$scriptFile = &scriptFileName($istage,$iscript, 0, $OPT_SCRIPT ); 
	$prefix = "$SCRIPT_PREFIX[$istage]" ;
	$CMD[0] = "${prefix}.pl ${lastPrefix}.input > ${prefix}.log" ; 
	&update_script($scriptFile, 1, @CMD );

	$istage = $ISTAGE_INCDATA; $iscript =0; @CMD = ( );
	$scriptFile = &scriptFileName($istage,$iscript,0,$OPT_SCRIPT );
	$prefix = "$SCRIPT_PREFIX[$istage]" ;
	$CMD[0] = "${prefix}.pl ${lastPrefix}.input > ${prefix}.log" ;
	&update_script($scriptFile, 1, @CMD );

    }

    # TRAINRUN --------------------------------------
    $istage = $ISTAGE_TRAINRUN ;  $iscript = 0;  @CMD = ( );
    $scriptFile = &scriptFileName($istage,$iscript,0,1);

    $prefix = "$SCRIPT_PREFIX[$istage]" ;
    $CMD[0] = "${prefix}.pl ${prefix}.input > ${prefix}.log" ;
    &update_script($scriptFile, 0, @CMD );


    # SIMGEN -------------------------------------
    if ( $SIMFLAG ) {
	for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) {  
	    &make_SIMGEN_BINARY_scripts($igroup); 
	    &make_SIMGEN_scripts($igroup); 
	}
    }

    # SIMFIT -------------------------------------
    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) 
    { &make_SIMFIT_scripts($igroup); }
    

    # SALT2mu -------------------------------------
    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) 
    { &make_SALT2mu_scripts($igroup, $ISTAGE_SALT2mu); }

    # COSMO ---------------------------------------
    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) 
    { &make_COSMO_scripts($igroup, $ISTAGE_COSMO); }

    # -------------------------------------
    #
    # Run these script subroutines if DMCORRFLAG is 1
    # 

    if ( $DMCORRFLAG ) {
	for ($igroup = 1; $igroup <= $NGROUP_DMCORR; $igroup++ ) { 
	    &make_DMFIXSIM_scripts($igroup); 
	    &make_DMFIXFIT_scripts($igroup);
	    &make_SALT2mu_scripts($igroup, $ISTAGE_DMFIXPREP);
	    &make_DMFIXRUN_scripts($igroup, $ISTAGE_DMFIXRUN);
	    &make_COSMO_scripts($igroup, $ISTAGE_DMFIXCOSMO);
        }
    }



    # PLOTS -------------------------------------
    if ( $SIMFLAG ) {
	for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) 
	{ &make_PLOT_scripts($igroup); }
    }

    # CLEAN -------------------
    &make_CLEAN_script();

    # ARCHIVE -----------------
    if ( $SIMFLAG ) {
	&make_ARCHIVE_script();
    }

    # ----- give +x privs ------
    $cmdx = "chmod +x $OUTPUT_TOPDIR/RUN*" ;
    qx($cmdx);

    # remind user where to run scripts
    print "\n Run scripts from :\n\t $OUTPUT_TOPDIR \n" ;


    if ( $OPT_PIPELINE eq "QUICKTEST" ) {
	print "  $OPT_PIPELINE mode -> replace training with symb. link\n";
	qx(touch $OUTPUT_TOPDIR/QUICKTEST);
    }

} # end of make_RUNscripts


# =================================
sub scriptFileName(@) {

    # create RUN* file and print a "cd $PATH" statement inside
    # opt=0 => just return name, but do NOT create file
    # opt=1 => create file

    my ($istage,$igroup,$iscript,$opt) = @_;

    my ($NAME, $PATH, $c3num, $c2stage, $c2group, $cdd, $tmpFile) ;
    my ($base, $scriptFile );

    $NAME    = $STAGE_NAME[$istage] ;
    $PATH    = $STAGE_PATH[$istage] ;
    $c3num   = sprintf("%3.3d", $iscript );
    $c2stage = sprintf("%2.2d", $istage  );
    $c2group = sprintf("%2.2d", $igroup  );

    $base  = "RUN${c2stage}_${NAME}" ;

    if ( $igroup > 0 ) { 
	if ( $istage == $ISTAGE_SALT2mu || $istage == $ISTAGE_PLOTS || $istage == $ISTAGE_DMFIXPREP ) 
	{ $tmpFile = "${base}-Group${c2group}" ;  }
	elsif ( $istage == $ISTAGE_DMFIXRUN || $istage == $ISTAGE_DMFIXCOSMO || $istage == $ISTAGE_SIMGEN ) 
	{ $tmpFile = "${base}-Group${c2group}" ;  }
	elsif ( $istage == $ISTAGE_DATABINARY || $istage == $ISTAGE_COSMO )
	{ $tmpFile = "${base}-Group${c2group}" ;  }
	elsif ( $istage == $ISTAGE_SIMFIT )
	{ $tmpFile = "${base}-Group${c2group}-${c3num}" ;  }	
	elsif ( $istage == $ISTAGE_DMFIXSIM || $istage == $ISTAGE_DMFIXFIT ) 
	{ $tmpFile = "${base}-Group${c2group}-${c3num}" ;  }	
	else
	{ $tmpFile = "${base}-${c2group}" ;  }    
    }
    else
    { $tmpFile = "${base}" ; }



    $scriptFile = "${tmpFile}" ;

    if ( $opt ) {
	print "\t Create script : '${tmpFile}' \n" ;
	$cdd = "cd $PATH" ;
	qx(cd $OUTPUT_TOPDIR ; echo "$cdd \n" > $scriptFile);
    }

    return $scriptFile ;

} # end of scriptFileName

# =================================
sub SNMixInfileName(@) {

    # create SNMIX* file and print batch/ssh info inside
    # opt=0 => just return name, but do NOT create file
    # opt=1 => create file

    my ($istage,$igroup,$iscript,$nbatch,$opt) = @_;

    my ($NAME, $PATH, $c3num, $c2stage, $c2group, $cdd, $tmpFile) ;
    my ($base, $scriptFile );

    $NAME    = $STAGE_NAME[$istage] ;
    $PATH    = $STAGE_PATH[$istage] ;
    $c3num   = sprintf("%3.3d", $iscript );
    $c2stage = sprintf("%2.2d", $istage  );
    $c2group = sprintf("%2.2d", $igroup  );

    $base  = "SNMIX${c2stage}_${NAME}" ;

    if ( $igroup > 0 ) { 
	if ( $istage == $ISTAGE_SIMGEN )
	{ $tmpFile = "${base}-Group${c2group}-${c3num}.input" ;  }	
	elsif ( $istage == $ISTAGE_DMFIXSIM ) 
	{ $tmpFile = "${base}-Group${c2group}-${c3num}.input" ;  }	
	else
	{ $tmpFile = "${base}-${c2group}.input" ;  }    
    }
    else
    { $tmpFile = "${base}.input" ; }

    $scriptFile = "${tmpFile}" ;

    if ( $opt ) {
	print "\t Create script : '${tmpFile}' \n" ;
	$cdd = "BATCH_INFO: $BATCH_COMMAND_FIT $BATCH_TEMPLATE_FIT $nbatch" ;
	qx(cd $PATH; echo "$cdd \n" > $scriptFile);
    }

    return $scriptFile ;

} # end of SNMixInfileName

# =========================
sub update_snmix {

    my($snmixFile, $PATH, @CMD ) = @_;

    my ($fullName,$i,$NCMD,$touchDone);

    $fullName = "$PATH/$snmixFile" ;
    open  PTR_SCRIPT  , ">> $fullName" ; 

    $NCMD = scalar(@CMD);
    for ( $i=0; $i < $NCMD ; $i++ ) {       
	print PTR_SCRIPT "$CMD[$i]\n" ;
    }

    close PTR_SCRIPT ;

} # end of update_snmix


# =========================
sub update_script {

    my($scriptFile, $DONEFLAG, @CMD ) = @_;

    my ($fullName,$i,$NCMD,$touchDone);

    $fullName = "$OUTPUT_TOPDIR/$scriptFile" ;
    open  PTR_SCRIPT  , ">> $fullName" ; 

    $NCMD = scalar(@CMD);
    for ( $i=0; $i < $NCMD ; $i++ ) {       
	print PTR_SCRIPT "$CMD[$i]\n" ;
    }

    # create a DONE file when the commands finish
    if ( $DONEFLAG ) {
	$touchDone = "touch ${fullName}.DONE" ;
	print PTR_SCRIPT "\n$touchDone \n" ;
    }

    close PTR_SCRIPT ;

} # end of update_script


# ===========================================
sub make_SIMGEN_scripts {

    my($igroup) = @_;

    my ($istage, $igen, $icmd, $cmd, $GROUPID, $OPT_SCRIPT );
    my ($simFile, @simFileList, @GROUPID_LIST, @GENV_LIST );
    my ($RUNOPT, $GENOPT);
    my ($scriptFile, $logFile, $GENVERSION, @CMD, $program );
    my ($snMixinput, $snMixPath, $snMixLog, $idot, $nsim);
    
    $istage   = $ISTAGE_SIMGEN ;  
    $OPT_SCRIPT = 1;

    $program      = $PROGRAM_NAME[$istage];
    @simFileList  = split(/\s+/,$SIMDATA_INFILE_GEN[$igroup]) ;
    @GROUPID_LIST = split(/\s+/,$DATA_GROUPID[$igroup]) ;
    @GENV_LIST    = split(/\s+/,$SIMDATA_GENVERSION[$igroup]) ;
    $RUNOPT       = "$SIMDATA_RUNOPT[$igroup]" ;
    $GENOPT       = "$SIMDATA_GENOPT[$igroup]" ;
    $igen         = 0 ;
    $nsim         = 1 ;

    $scriptFile = &scriptFileName($istage,$igroup,$GROUPID,$OPT_SCRIPT);
    
    foreach $simFile ( @simFileList ) {

	$GROUPID     = $GROUPID_LIST[$igen];
	$GENVERSION  = $GENV_LIST[$igen]; 

	# create SIMGEN script and sim_SNmix input file
        $snMixinput = &SNMixInfileName($istage, $igroup, $GROUPID, $nsim, $OPT_SCRIPT);
        $snMixPath = $STAGE_PATH[$istage];
        $idot   = rindex($snMixinput,".input");
        $snMixLog = substr($snMixinput,0,$idot);

        # write RUNIT script first, then snMix input file
	# RUN snMix
	@CMD = ( ) ; $icmd = 0;
	$CMD[$icmd] = "$program $snMixinput  " ;
	# add snMix SIMRUN-options
	if ( "$RUNOPT" ne '' ) {
	    $CMD[$icmd] = $CMD[$icmd] . $RUNOPT ;
	}
	$icmd++;
	$CMD[$icmd] = " ";
	&update_script($scriptFile, 0, @CMD);

        # DONE WITH RUNIT SCRIPT
        # START snMIX input file

	# add DONE stamp and CIDOFF keyword
	@CMD = ( ); $icmd = 0;
        $CMD[$icmd] = "DONE_STAMP: ${OUTPUT_TOPDIR}/${scriptFile}.DONE " ;
        $icmd++; $CMD[$icmd] = "RESET_CIDOFF: 0 ";
        $icmd++; $CMD[$icmd] = " ";
        &update_snmix($snMixinput, $snMixPath, @CMD);       

	# add GEN keywords
	@CMD = ( ); $icmd = 0;
	$CMD[$icmd]   = "GENVERSION:  ${GENVERSION} ";

	# specify (optional) global GENMODEL
	if ( $GLOBAL_GENMODEL ne '' & $SIMDATA_EXCLUDE_FITOPT[$igroup] == 0 ) {
	    $icmd++ ; $CMD[$icmd] = "GENOPT: GENMODEL ${GLOBAL_GENMODEL} ";
	}
	
	#  include SIMGEN-options
	if ( "$GENOPT" ne '' ) {
	    $icmd++ ; $CMD[$icmd] = "GENOPT: ${GENOPT} " ;
	}

	$icmd++ ; $CMD[$icmd]   = " ";
        # add global input values before starting next script
        $icmd++ ; $CMD[$icmd] = "ENDLIST_GENVERSION ";
	$icmd++ ; $CMD[$icmd]   = " ";
	$icmd++ ; $CMD[$icmd]   = "SIMGEN_INFILE_Ia: $simFile ";
	$icmd++ ; $CMD[$icmd]   = "GENPREFIX: $snMixLog ";
	$icmd++ ; $CMD[$icmd]   = " ";
	# update script with commands
	&update_snmix($snMixinput, $snMixPath, @CMD);
	$igen++ ;  # local loop counter
    }

}  # end of  make_SIMGEN_scripts

# ===========================================
sub make_SIMGEN_BINARY_scripts {

    my($igroup) = @_;

    my ($istage, $igen, $icmd, $cmd, $GROUPID, $OPT_SCRIPT );
    my ($simFile, @simFileList, @GROUPID_LIST, @GENV_LIST, $GENOPT);
    my ($scriptFile, $logFile, $GENVERSION, @CMD, $program );
    my ($c2g);
    
    $istage   = $ISTAGE_DATABINARY ;  
    $OPT_SCRIPT = 1;

    $program      = $PROGRAM_NAME[$istage];
    @simFileList  = split(/\s+/,$SIMDATA_INFILE_GEN[$igroup]) ;
    @GROUPID_LIST = split(/\s+/,$DATA_GROUPID[$igroup]) ;
    @GENV_LIST    = split(/\s+/,$SIMDATA_GENVERSION[$igroup]) ;
    $GENOPT       = "$SIMDATA_GENOPT[$igroup]" ;
    $igen         = 0 ;

    # create SIMGEN script
    $scriptFile = &scriptFileName($istage,$igroup,$GROUPID,$OPT_SCRIPT);
    
    foreach $simFile ( @simFileList ) {

	$GROUPID     = $GROUPID_LIST[$igen];
	$GENVERSION  = $GENV_LIST[$igen]; 
	@CMD = ( ) ;

	$icmd = 0;   $CMD[$icmd] = "$program $simFile \\" ;
	
	#  include SIMGEN-options
	if ( "$GENOPT" ne '' ) {
	    $icmd++ ;   $CMD[$icmd] = "   ${GENOPT} \\" ;
	}

	# specify (optional) global GENMODEL
	if ( $GLOBAL_GENMODEL ne '' & $SIMDATA_EXCLUDE_FITOPT[$igroup] == 0 ) {
	    $icmd++ ; $CMD[$icmd] = "   GENMODEL ${GLOBAL_GENMODEL} \\";
	}

	# specify unique GENVERSION
	$icmd++ ; $CMD[$icmd]   = "   GENVERSION ${GENVERSION} \\";

	# override NGEN_LC for BINARY check/creation
	$icmd++ ; $CMD[$icmd]   = "   NGENTOT_LC 0  NGEN_LC 1 \\";

	# add pipe into log file
	$c2g =  sprintf("%2.2d", $igen);
	$logFile = "${scriptFile}-${igen}.LOG" ;
	$icmd++  ; $CMD[$icmd] = "   >& ${logFile}" ;
	
	$icmd++ ; $CMD[$icmd]   = " ";
	# set o-w priv
	$cmd = "chmod o-w $SNDATA_ROOT/SIM/$GENVERSION" ;
	$icmd++  ; $CMD[$icmd] = "$cmd" ;
	$icmd++  ; $CMD[$icmd] = " " ;

	# update script with commands
	&update_script($scriptFile, 1, @CMD);

	$igen++ ;  # local loop counter
    }

}  # end of  make_SIMGEN_BINARY_scripts

# ===========================================
sub make_DMFIXSIM_scripts {

    my($igroup) = @_;

    my ($istage, $igen, $icmd, $cmd, $GROUPID, $OPT_SCRIPT, $RUNOPT );
    my ($simFile, @simFileList, @GROUPID_LIST, @GENV_LIST, $GENOPT);
    my ($scriptFile, $GENVERSION, @CMD, $program, $usegenversion );
    
    my ($imodel, $nmodel, $model, $modelpath, $genflag, $genmodel, $c3 );
    my ($fitres, $compgroup, $comppath, $script, $modsmear);
    my ($groupname, $meanc, $siglc, $sigrc, $meanx, $siglx, $sigrx);
    my ($prefix, $simcodeopt, $program);

    my ($snMixinput, $snMixPath, $snMixLog, $idot);

    $istage   = $ISTAGE_DMFIXSIM ;  
    $OPT_SCRIPT = 1;

    $groupname    = $DMCORR_GROUPNAME[$igroup];
    $compgroup    = $DATA_GROUPNAME[$DMCORR_IDEAL[$igroup]];
    $program      = $PROGRAM_NAME[$istage];
    @simFileList  = split(/\s+/,$DMCORR_INFILE_GEN[$igroup]) ;
    @GROUPID_LIST = split(/\s+/,$DMCORR_GROUPID[$igroup]) ;
    @GENV_LIST    = split(/\s+/,$DMCORR_GENVERSION[$igroup]) ;
    $RUNOPT       = "$DMCORR_RUNOPT[$igroup]";
    $GENOPT       = "$DMCORR_GENOPT[$igroup]" ;
    $simcodeopt   = "$DMCORR_SIMCODEOPT[$igroup]";
    $igen         = 0 ;

    if ($DMCORR_SIMCODE[$igroup] eq "DEFAULT" )  {
	$script       = "${SCRIPT_PREFIX[$istage]}.pl";
    } else {
	$script       = $DMCORR_SIMCODE[$igroup] ;
    }


    if ( $DMCORR_SIMTYPE[$igroup] eq "AS-IS" ) {
	$nmodel = 1; 
	$genflag = 1;
    } else {
	$nmodel = $NMODEL_TRAIN;
	$genflag = 0;
	$meanc = $DMCORR_MEANC[$igroup];
	$meanx = $DMCORR_MEANX[$igroup];
	$siglx = $DMCORR_SIGLX[$igroup];
	$sigrx = $DMCORR_SIGRX[$igroup];
	$siglc = $DMCORR_SIGLC[$igroup];
	$sigrc = $DMCORR_SIGRC[$igroup];
    }

    foreach $simFile ( @simFileList ) {

       $GROUPID     = $GROUPID_LIST[$igen];
       $GENVERSION  = $GENV_LIST[$igen]; 
	    
       # create SIMGEN script and sim_SNmix input file
       $scriptFile = &scriptFileName($istage,$igroup,$GROUPID,$OPT_SCRIPT);
       $snMixinput = &SNMixInfileName($istage, $igroup, $GROUPID, $nmodel, $OPT_SCRIPT);
       $snMixPath = $STAGE_PATH[$istage];
       $idot   = rindex($snMixinput,".input");
       $snMixLog = substr($snMixinput,0,$idot);

       # write RUNIT script first, then snMix input file
       # RUN snMix
       @CMD = ( ) ; $icmd = 0;
       $CMD[$icmd] = "$program $snMixinput " ;
	# add snMix SIMRUN-options
	if ( "$RUNOPT" ne '' ) {
	    $CMD[$icmd] = $CMD[$icmd] . $RUNOPT ;
	}
       $icmd++;
       $CMD[$icmd] = " ";
       &update_script($scriptFile, 0, @CMD);

       # DONE WITH RUNIT SCRIPT
       # START snMIX input file

       # add DONE stamp and CIDOFF keyword
       @CMD = ( ); $icmd = 0;
       $CMD[$icmd] = "DONE_STAMP: ${OUTPUT_TOPDIR}/${scriptFile}.DONE " ;
       $icmd++; $CMD[$icmd] = "RESET_CIDOFF: 0 ";
       $icmd++; $CMD[$icmd] = " ";
       &update_snmix($snMixinput, $snMixPath, @CMD);       

       # if train-based sims, add SNANA_MODELPATH
       @CMD = ( ); $icmd = 0;
       if ( $genflag == 0 ) { 
	   $CMD[$icmd] = "SNANA_MODELPATH: $SNANA_MODELPATH ";
	   $icmd++;
       }
       $CMD[$icmd] = " ";
       &update_snmix($snMixinput, $snMixPath, @CMD);       

       # loop over nmodel if user requests train-based sims
       for ( $imodel = 1; $imodel <= $nmodel; $imodel ++ ) {
	    @CMD = ( ) ; $icmd = 0; 
	    $usegenversion = $GENVERSION;

	    # for trainbased sims, get model info
	    $c3 = $TRAIN_CHARINDEX[$imodel];
	    $model = $TRAIN_MODELNAME[$imodel];
	    $genmodel = "SALT2.${model}";
	    # for trainbased sims, get comp fitres info
	    $prefix = &prefixSALT2mu($imodel,$ISTAGE_SALT2mu);
	    $comppath = $STAGE_PATH[$ISTAGE_SALT2mu];	    
	    $fitres = "${comppath}/${compgroup}/${prefix}.fitres";
	    # for trainbased sims, get color smear info
            # [JLM] May 25 2013 -- alter modsmear to allow QUICKTEST compatibility
	    #$modsmear = "${STAGE_PATH[$ISTAGE_TRAINRUN]}/workSpace/${model}/salt2_color_dispersion.dat";
	    $modsmear = "${SNANA_MODELPATH}/SALT2.${model}/salt2_color_dispersion.dat";
	    # obsolete w/ updated pcafit $locsmear = "salt2_color_dispersion_${c3}.dat";
	    
      	    # specify (optional) global GENMODEL
	    if ( $DMCORR_EXCLUDE_FITOPT[$igroup] == 0 && $genflag == 1) {
		$CMD[$icmd] = "GENVERSION: $usegenversion ";
		$icmd++; $CMD[$icmd] = "GENOPT: GENMODEL ${GLOBAL_GENMODEL} ";
	    } elsif ( $genflag == 0 ) {
		$usegenversion = "${GENVERSION}_${c3}";
		$CMD[$icmd] = "GENVERSION: $usegenversion ";
		$icmd++ ; $CMD[$icmd] = "GENOPT: GENMODEL ${genmodel} ";
		$icmd++ ; $CMD[$icmd] = "GENOPT: SALT2mu_FILE ${fitres} ";
		$icmd++ ; $CMD[$icmd] = "GENOPT: GENMEAN_SALT2c  ${meanc} GENSIGMA_SALT2c  ${siglc} ${sigrc} ";
		$icmd++ ; $CMD[$icmd] = "GENOPT: GENMEAN_SALT2x1 ${meanx} GENSIGMA_SALT2x1 ${siglx} ${sigrx} ";
		$icmd++ ; $CMD[$icmd] = "GENOPT: GENMAG_SMEAR_MODELNAME G10 ";
		$icmd++ ; $CMD[$icmd] = "GENOPT: MODELSMEAR_FILE ${modsmear} ";
	    } else {
		$CMD[$icmd] = "GENVERSION: $usegenversion ";
	    }
	    #  include SIMGEN-options
	    if ( "$GENOPT" ne '' ) {
		$icmd++ ; $CMD[$icmd] = "GENOPT: ${GENOPT} " ;
	    }

	    $icmd++ ; $CMD[$icmd]   = " ";

	    # update snMix input file with commands
            &update_snmix($snMixinput, $snMixPath, @CMD);       

	} # end of loop over genversions
        
        # add global input values before starting next script
        @CMD = ( ) ; $icmd = 0;
        $CMD[$icmd] = "ENDLIST_GENVERSION ";
        $icmd++ ; $CMD[$icmd] = " ";
        $icmd++ ; $CMD[$icmd] = "SIMGEN_INFILE_Ia: $simFile ";
        $icmd++ ; $CMD[$icmd] = "GENPREFIX: $snMixLog ";
        $icmd++ ; $CMD[$icmd] = " ";
        &update_snmix($snMixinput, $snMixPath, @CMD);       
	$igen++ ;  # local loop counter
   }

}  # end of  make_DMFIXSIM_scripts



# ===========================================
sub make_SIMFIT_scripts {

    # 
    # construct script(s) with  'split_and_fit' command.
    #

    my($igroup) = @_;

    my ($istage, $icmd, $PATH, $igen, $GROUPID, $scriptFile) ;
    my ($fitScript, $DONEKEY, $sedcmd, $nmlFile, $FITRUNOPT) ;
    my (@FITNML_FILE_LIST, @GROUPID_LIST, @CMD ) ;
    my $OPT_SCRIPT=1;

    $istage    = $ISTAGE_SIMFIT ;  
    $PATH      = $STAGE_PATH[$istage] ;
    $fitScript = "$SCRIPT_PREFIX[$istage].pl" ;
    $FITRUNOPT = $GLOBAL_FIT_RUNOPT;

    @FITNML_FILE_LIST  = split(/\s+/,$DATA_FITNML[$igroup]) ;
    @GROUPID_LIST      = split(/\s+/,$DATA_GROUPID[$igroup]) ;
    $igen              = 0 ;

    foreach $nmlFile ( @FITNML_FILE_LIST ) {

	$GROUPID    = $GROUPID_LIST[$igen];

#	print " xxxx GROUP=$GROUPID  nmlFile=$nmlFile ... \n";
	$scriptFile = &scriptFileName($istage,$igroup,$GROUPID,$OPT_SCRIPT);

	@CMD        = ( ) ;
	$icmd       = 0;  
	$CMD[$icmd] = "$fitScript  $nmlFile $FITRUNOPT" ;

	# update nmlFile with DONE_STAMP so that DONE_STAMP appears
	# after all fits are done rather than when the jobs are submitted.
	$DONEKEY = "DONE_STAMP:   $OUTPUT_TOPDIR/${scriptFile}.DONE" ;
	$sedcmd  = "sed -e '1 i\\${DONEKEY}' " ;
	qx(cd $PATH; $sedcmd $nmlFile > DUMDUM ; mv DUMDUM $nmlFile);

	&update_script($scriptFile, 0, @CMD ); # 0 => do NOT touch DONE file
	$igen++ ;  # local loop counter
    }

}  # end of  make_SIMFIT_scripts

# ===========================================
sub make_DMFIXFIT_scripts {

    # 
    # construct script(s) with  'split_and_fit' command.
    #

    my($igroup) = @_;

    my ($istage, $icmd, $PATH, $igen, $GROUPID, $scriptFile) ;
    my ($fitScript, $DONEKEY, $sedcmd, $nmlFile, $FITRUNOPT) ;
    my (@FITNML_FILE_LIST, @GROUPID_LIST, @CMD ) ;
    my $OPT_SCRIPT=1;

    $istage    = $ISTAGE_DMFIXFIT ;  
    $PATH      = $STAGE_PATH[$istage] ;
    $fitScript = "$SCRIPT_PREFIX[$istage].pl" ;
    $FITRUNOPT    = $GLOBAL_FIT_RUNOPT;

    @FITNML_FILE_LIST  = split(/\s+/,$DMCORR_FITNML[$igroup]) ;
    @GROUPID_LIST      = split(/\s+/,$DMCORR_GROUPID[$igroup]) ;
    $igen              = 0 ;

    foreach $nmlFile ( @FITNML_FILE_LIST ) {

	$GROUPID    = $GROUPID_LIST[$igen];

#	print " xxxx GROUP=$GROUPID  nmlFile=$nmlFile ... \n";
	$scriptFile = &scriptFileName($istage,$igroup,$GROUPID,$OPT_SCRIPT);

	@CMD        = ( ) ;
	$icmd       = 0;  
	$CMD[$icmd] = "$fitScript  $nmlFile $FITRUNOPT" ;

	# update nmlFile with DONE_STAMP so that DONE_STAMP appears
	# after all fits are done rather than when the jobs are submitted.
	$DONEKEY = "DONE_STAMP:   $OUTPUT_TOPDIR/${scriptFile}.DONE" ;
	$sedcmd  = "sed -e '1 i\\${DONEKEY}' " ;
	qx(cd $PATH; $sedcmd $nmlFile > DUMDUM ; mv DUMDUM $nmlFile);

	&update_script($scriptFile, 0, @CMD ); # 0 => do NOT touch DONE file
	$igen++ ;  # local loop counter
    }

}  # end of  make_DMFIXFIT_scripts


# ===========================================
sub make_SALT2mu_scripts {
    
    # For this $igroup, construct script to catenate
    # relevant FITRES files and to run SALT2mu.
    # Loop over all training models too.
    #
    # - For SSH/NODELIST, run jobs in parallel on @NODELIST   
    # - For batch-mode, distribute jobs in parallel
    #
    # Use 'wait_for_files.pl' utility to write master doneFile
    # when all of the SALT2mu jobs finish.
    # 
    # ------------

    my($igroup, $istage) = @_;

    my $OPT_SCRIPT = 1 ;
    my $qq         = '"' ;

    my ($mainScriptFile, $subScriptFile, $prefix) ;
    my ($itrain, $group, $PATH, $topDir, $sDir);
    my ($NCMD, $cmd, @CMD, @CMDMAIN, $NMAIN, $inode, $node ) ;
    my $nfitopt;

    my (@local_exclude_fitopt, @local_groupname);

    if ($istage == $ISTAGE_SALT2mu) {
	@local_exclude_fitopt = @SIMDATA_EXCLUDE_FITOPT;
	@local_groupname = @DATA_GROUPNAME;
    } elsif ($istage == $ISTAGE_DMFIXPREP) {
	@local_exclude_fitopt = @DMCORR_EXCLUDE_FITOPT;
	@local_groupname = @DMCORR_GROUPNAME; 
    }


    $nfitopt = $NMODEL_TRAIN; #default, unless SIMFLAG & SIMDATA_EXCLUDE_FITOPT[$igroup]==1
    if ($SIMFLAG) {
	if ($local_exclude_fitopt[$igroup]) {$nfitopt=1;}
    }

    $group    = $local_groupname[$igroup] ;
    $PATH     = "$STAGE_PATH[$istage]" ;
    $topDir   = "$PATH/$group" ;

    $mainScriptFile = &scriptFileName($istage,$igroup,0,$OPT_SCRIPT);

    $NMAIN   = -1;
    $NMAIN++ ; $CMDMAIN[$NMAIN] = "cd $group \n" ;

    # loop over training sets and make subdir for each

    $inode = -1;

#    print "xxx in make_SALT2mu_scripts \n";
#    print "xxx looping over trains. istage is $istage \n\n";

    for ( $itrain = 1; $itrain <= $nfitopt; $itrain++ ) {
	$prefix      = &prefixSALT2mu($itrain, $istage);

	if ( $BATCH_FLAG ) { 
	    $subScriptFile = "${prefix}.batch" ; 
	    $NMAIN++; $CMDMAIN[$NMAIN] = "$BATCH_COMMAND_FIT $subScriptFile" ;
	}
	else { 
	    $inode++ ; 
	    if ( $inode >= $NNODE ) { $inode = 0; }
	    $node = "$NODELIST[$inode]" ;
	    $NODELIST_SALT2mu[$igroup][$itrain] = $node ;
	    $subScriptFile = "${prefix}_$node" ; 
	    $NMAIN++; $CMDMAIN[$NMAIN] = 
		"ssh -x $node $qq cd $topDir ; ./$subScriptFile $qq &" ;
	}

	# construct commands for batch file(s)
	@CMD = &get_SALT2mu_commands($igroup,$itrain,$istage);

#	print " xxx subScript file is $topDir/$subScriptFile \n";
        sntools::checkDirExist($topDir, "can't find $topDir");
	if ( $NNODE > 0 ) {
	    open  PTR_subScript  , ">> $topDir/$subScriptFile" ; 
	    foreach $cmd ( @CMD ) { print PTR_subScript "$cmd \n"; }
	    close PTR_subScript ;
	}
	else {
	    my @REPLACE_KEY = ( "REPLACE_NAME", "REPLACE_LOGFILE", 
				"REPLACE_MEM", "REPLACE_JOB" );
	    my (@REPLACE_STRING, $inF, $outF, $NKEY, $CMDBATCH) ;

	    $CMDBATCH = "";
	    foreach $cmd ( @CMD ) { $CMDBATCH = "$CMDBATCH  $cmd ;" ; }
	    $REPLACE_STRING[0] = "$subScriptFile" ;
	    $REPLACE_STRING[1] = "${subScriptFile}.LOG" ;
	    $REPLACE_STRING[2] = "500mb" ;
	    $REPLACE_STRING[3] = "$CMDBATCH";
	    $inF   = "$BATCH_TEMPLATE_FIT" ;
	    $outF  = "${topDir}/$subScriptFile" ;
	    $NKEY = scalar(@REPLACE_KEY);
	    sntools::sed_strReplace($inF, $outF, $NKEY, 
				    @REPLACE_KEY, @REPLACE_STRING );
	}

	# set +x priv
	qx(chmod +x $topDir/$subScriptFile);

    }  # itrain loop



    # use utility to write global DONE file when all of the
    # individual DONE files exist.
    my  $doneFile   = "$OUTPUT_TOPDIR/${mainScriptFile}.DONE" ;
    my ($waitArgs);
    
    $waitArgs = "$nfitopt  '*.DONE'  $doneFile" ;
    $NMAIN++; $CMDMAIN[$NMAIN] = "\n$WAIT_SCRIPT  $waitArgs" ;

    # -------------------------------------
    # create main script for SALT2mu jobs.
    &update_script($mainScriptFile, 0, @CMDMAIN );

} # end of make_SALT2mu_scripts 

# ===========================================
sub make_TRAINGEN_scripts(@) {
    
    # For this $igroup, construct script to ssh/batch
    # submit individual TRAIN_SIM jobs. 
    #
    # - For SSH/NODELIST, run jobs in parallel on @NODELIST   
    # - For batch-mode, distribute jobs in parallel
    #
    # Use 'wait_for_files.pl' utility to write master doneFile
    # when all of the TRAIN_SIM.pl jobs finish.
    # 
    # ------------

    my($istage) = @_;

    my $OPT_SCRIPT = 1 ;
    my $qq         = '"' ;

    my ($mainScriptFile, $subScriptFile, $inFile, $prefix) ;
    my ($itrain, $group, $PATH, $sDir, $c2t);
    my ($NCMD, $cmd, @CMD, @CMDMAIN, $NMAIN, $inode, $node ) ;
    my ($icmd, $logFile, $doneFile);

    $istage = $ISTAGE_TRAINGEN ; 
    $PATH     = "$STAGE_PATH[$istage]" ;

    $mainScriptFile = &scriptFileName($istage,0,0,$OPT_SCRIPT);

    $NMAIN   = -1;

    # loop over training sets and make subdir for each

    $inode = -1;
    $prefix      = $SCRIPT_PREFIX[$istage] ;

#    print "xxx in make_TRAINGEN_scripts \n";
#    print "xxx looping over trains. istage is $istage \n\n";

    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {
	$c2t =  sprintf("%2.2d", $itrain);
	$inFile = "${prefix}_${c2t}.input";
	$logFile = "${prefix}_${c2t}.log";
	$doneFile = "${prefix}_${c2t}.DONE";

	if ( $BATCH_FLAG ) { 
	    $subScriptFile = "${prefix}_${c2t}.batch" ; 
	    $NMAIN++; $CMDMAIN[$NMAIN] = "$BATCH_COMMAND_FIT $subScriptFile" ;
	}
	else { 
	    $inode++ ; 
	    if ( $inode >= $NNODE ) { $inode = 0; }
	    $node = "$NODELIST[$inode]" ;
	    $NODELIST_SALT2mu[$igroup][$itrain] = $node ;
	    $subScriptFile = "${prefix}_${c2t}_$node" ; 
	    $NMAIN++; $CMDMAIN[$NMAIN] = 
		"ssh -x $node $qq cd $PATH ; ./$subScriptFile $qq &" ;
	}

	# construct commands for batch file(s)
	$CMD[0]     = "${prefix}.pl ${inFile} > ${logFile}" ;
	$CMD[1]     = "touch $doneFile" ;

#	print " xxx subScript file is $topDir/$subScriptFile \n";
        sntools::checkDirExist($PATH, "can't find $PATH");
	if ( $NNODE > 0 ) {
	    open  PTR_subScript  , ">> $PATH/$subScriptFile" ; 
	    foreach $cmd ( @CMD ) { print PTR_subScript "$cmd \n"; }
	    close PTR_subScript ;
	}
	else {
	    my @REPLACE_KEY = ( "REPLACE_NAME", "REPLACE_LOGFILE", 
				"REPLACE_MEM", "REPLACE_JOB" );
	    my (@REPLACE_STRING, $inF, $outF, $NKEY, $CMDBATCH) ;

	    $CMDBATCH = "";
	    foreach $cmd ( @CMD ) { $CMDBATCH = "$CMDBATCH  $cmd ;" ; }
	    $REPLACE_STRING[0] = "$subScriptFile" ;
	    $REPLACE_STRING[1] = "${subScriptFile}.LOG" ;
	    $REPLACE_STRING[2] = "500mb" ;
	    $REPLACE_STRING[3] = "$CMDBATCH";
	    $inF   = "$BATCH_TEMPLATE_TRAIN" ;
	    $outF  = "${PATH}/$subScriptFile" ;
	    $NKEY = scalar(@REPLACE_KEY);
	    sntools::sed_strReplace($inF, $outF, $NKEY, 
				    @REPLACE_KEY, @REPLACE_STRING );
	}

	# set +x priv
	qx(chmod +x $PATH/$subScriptFile);

    }  # itrain loop



    # use utility to write global DONE file when all of the
    # individual DONE files exist.
    my  $ALLdoneFile   = "$OUTPUT_TOPDIR/${mainScriptFile}.DONE" ;
    my ($waitArgs);
    
    $waitArgs = "$NMODEL_TRAIN  '*.DONE'  $ALLdoneFile" ;
    $NMAIN++; $CMDMAIN[$NMAIN] = "\n$WAIT_SCRIPT  $waitArgs" ;

    # -------------------------------------
    # create main script for SALT2mu jobs.
    &update_script($mainScriptFile, 0, @CMDMAIN );

} # end of make_TRAINSIM_scripts

# ===========================================
sub make_COSMO_scripts {

    # For this $igroup, construct script to copy 
    # relevant FITRES files and to run cosmo program.
    # Loop over all training models too.
    #
    # - For SSH/NODELIST, run jobs in parallel on @NODELIST
    # - For batch-mode, distribute jobs in parallel
    #
    # Use 'wait_for_files.pl' utility to write master doneFile
    # when all of the COSMO jobs finish.
    #
    # ------------

    my($igroup, $istage) = @_;


    my $OPT_SCRIPT = 1 ;
    my $qq         = '"' ;

    my ($mainScriptFile, $subScriptFile, $prefix, $COSMOdoneflag) ;
    my ($itrain, $group, $PATH, $topDir, $sDir);
    my ($NCMD, $cmd, @CMD, @CMDMAIN, $NMAIN, $inode, $node ) ;
    my $nfitopt;
    my (@local_exclude_fitopt, @local_groupname);

    if ( $istage == $ISTAGE_COSMO ) {
	@local_exclude_fitopt = @SIMDATA_EXCLUDE_FITOPT;
	@local_groupname = @DATA_GROUPNAME;
    } elsif ( $istage == $ISTAGE_DMFIXCOSMO ) {
	@local_exclude_fitopt = @DMCORR_EXCLUDE_FITOPT;
	@local_groupname = @DMCORR_GROUPNAME;
    }

    $nfitopt = $NMODEL_TRAIN; #default, unless SIMFLAG & SIMDATA_EXCLUDE_FITOPT[$igroup]==1
    if ($SIMFLAG) {
	if ($local_exclude_fitopt[$igroup]) {$nfitopt=1;}
    }

    if ( $OPT_COSMO == 0 ) {
	@CMD[0] = " ";
	$igroup = 0;
	$COSMOdoneflag = 1;
        $mainScriptFile = &scriptFileName($istage,$igroup,0,$OPT_SCRIPT);
        goto UPDATE_SCRIPT ;
    }

    $COSMOdoneflag = 0;
    $group    = $local_groupname[$igroup] ;
    $PATH     = "$STAGE_PATH[$istage]" ;
    $topDir   = "$PATH/$group" ;

    $mainScriptFile = &scriptFileName($istage,$igroup,0,$OPT_SCRIPT);

    $NMAIN   = -1;
    $NMAIN++ ; $CMDMAIN[$NMAIN] = "cd $group \n" ;

    # loop over training sets and make subdir for each

    $inode = -1;

    for ( $itrain = 1; $itrain <= $nfitopt; $itrain++ ) {
        $prefix      = &prefixSALT2mu($itrain, $istage);

        if ( $BATCH_FLAG ) {
            $subScriptFile = "COSMO.${prefix}.batch" ;
            $NMAIN++; $CMDMAIN[$NMAIN] = "$BATCH_COMMAND_FIT $subScriptFile" ;
        }
        else {
            $inode++ ;
            if ( $inode >= $NNODE ) { $inode = 0; }
            $node = "$NODELIST[$inode]" ;
            $NODELIST_SALT2mu[$igroup][$itrain] = $node ;
            $subScriptFile = "COSMO.${prefix}_$node" ;
	    "ssh -x $node $qq cd $topDir ; ./$subScriptFile $qq &" ;
        }

        # construct commands for batch file(s)
	@CMD = &get_COSMO_commands($igroup, $itrain, $istage);

        if ( $NNODE > 0 ) {
            open  PTR_subScript  , ">> $topDir/$subScriptFile" ;
            foreach $cmd ( @CMD ) { print PTR_subScript "$cmd \n"; }
            close PTR_subScript ;
        }
        else {
            my @REPLACE_KEY = ( "REPLACE_NAME", "REPLACE_LOGFILE",
                                "REPLACE_MEM", "REPLACE_JOB" );
            my (@REPLACE_STRING, $inF, $outF, $NKEY, $CMDBATCH) ;

            $CMDBATCH = "";
            foreach $cmd ( @CMD ) { $CMDBATCH = "$CMDBATCH  $cmd ;" ; }
            $REPLACE_STRING[0] = "$subScriptFile" ;
            $REPLACE_STRING[1] = "${subScriptFile}.LOG" ;
            $REPLACE_STRING[2] = "500mb" ;
            $REPLACE_STRING[3] = "$CMDBATCH";
            $inF   = "$BATCH_TEMPLATE_FIT" ;
            $outF  = "${topDir}/$subScriptFile" ;
            $NKEY = scalar(@REPLACE_KEY);
            sntools::sed_strReplace($inF, $outF, $NKEY,
                                    @REPLACE_KEY, @REPLACE_STRING );
        }

        # set +x priv
        qx(chmod +x $topDir/$subScriptFile);

    }  # itrain loop	

    # use utility to write global DONE file when all of the
    # individual DONE files exist.
    my  $doneFile   = "$OUTPUT_TOPDIR/${mainScriptFile}.DONE" ;
    my ($waitArgs);

    $waitArgs = "$nfitopt  '*.DONE'  $doneFile" ;
    $NMAIN++; $CMDMAIN[$NMAIN] = "\n$WAIT_SCRIPT  $waitArgs" ;

    # -------------------------------------
    # finish main script for COSMO jobs.
  UPDATE_SCRIPT:
    &update_script($mainScriptFile, $COSMOdoneflag, @CMDMAIN );

} # end of make_COSMO_scripts

sub make_DMFIXRUN_scripts {

    # For this $igroup, construct script to copy 
    # relevant FITRES files and to run dmfix program.
    # Loop over all training models too.
    #
    # - For SSH/NODELIST, run jobs in parallel on @NODELIST
    # - For batch-mode, distribute jobs in parallel
    #
    # Use 'wait_for_files.pl' utility to write master doneFile
    # when all of the DMFIXRUN jobs finish.
    #
    # ------------

    my($igroup, $istage) = @_;

    my $OPT_SCRIPT = 1 ;
    my $qq         = '"' ;

    my ($mainScriptFile, $subScriptFile, $prefix, $DMFIXRUNdoneflag) ;
    my ($itrain, $group, $PATH, $topDir, $sDir);
    my ($NCMD, $cmd, @CMD, @CMDMAIN, $NMAIN, $inode, $node ) ;
    my ($nfitopt);
    my (@local_exclude_fitopt, @local_groupname);

	@local_exclude_fitopt = @DMCORR_EXCLUDE_FITOPT;
	@local_groupname = @DMCORR_GROUPNAME;

    $nfitopt = $NMODEL_TRAIN; #default, unless SIMFLAG & SIMDATA_EXCLUDE_FITOPT[$igroup]==1
	if ($local_exclude_fitopt[$igroup]) {$nfitopt=1;}

    $DMFIXRUNdoneflag = 0;
    $group    = $local_groupname[$igroup] ;
    $PATH     = "$STAGE_PATH[$istage]" ;
    $topDir   = "$PATH/$group" ;

    $mainScriptFile = &scriptFileName($istage,$igroup,0,$OPT_SCRIPT);

    $NMAIN   = -1;
    $NMAIN++ ; $CMDMAIN[$NMAIN] = "cd $group \n" ;

    # loop over training sets and make subdir for each

    $inode = -1;

    for ( $itrain = 1; $itrain <= $nfitopt; $itrain++ ) {
        $prefix      = &prefixSALT2mu($itrain, $istage);

        if ( $BATCH_FLAG ) {
            $subScriptFile = "DMFIXRUN.${prefix}.batch" ;
            $NMAIN++; $CMDMAIN[$NMAIN] = "$BATCH_COMMAND_FIT $subScriptFile" ;
        }
        else {
            $inode++ ;
            if ( $inode >= $NNODE ) { $inode = 0; }
            $node = "$NODELIST[$inode]" ;
            $NODELIST_SALT2mu[$igroup][$itrain] = $node ;
            $subScriptFile = "DMFIXRUN.${prefix}_$node" ;
	    "ssh -x $node $qq cd $topDir ; ./$subScriptFile $qq &" ;
        }

        # construct commands for batch file(s)
	@CMD = &get_DMFIXRUN_commands($igroup, $itrain, $istage);

        if ( $NNODE > 0 ) {
            open  PTR_subScript  , ">> $topDir/$subScriptFile" ;
            foreach $cmd ( @CMD ) { print PTR_subScript "$cmd \n"; }
            close PTR_subScript ;
        }
        else {
            my @REPLACE_KEY = ( "REPLACE_NAME", "REPLACE_LOGFILE",
                                "REPLACE_MEM", "REPLACE_JOB" );
            my (@REPLACE_STRING, $inF, $outF, $NKEY, $CMDBATCH) ;

            $CMDBATCH = "";
            foreach $cmd ( @CMD ) { $CMDBATCH = "$CMDBATCH  $cmd ;" ; }
            $REPLACE_STRING[0] = "$subScriptFile" ;
            $REPLACE_STRING[1] = "${subScriptFile}.LOG" ;
            $REPLACE_STRING[2] = "500mb" ;
            $REPLACE_STRING[3] = "$CMDBATCH";
            $inF   = "$BATCH_TEMPLATE_FIT" ;
            $outF  = "${topDir}/$subScriptFile" ;
            $NKEY = scalar(@REPLACE_KEY);
            sntools::sed_strReplace($inF, $outF, $NKEY,
                                    @REPLACE_KEY, @REPLACE_STRING );
        }

        # set +x priv
        qx(chmod +x $topDir/$subScriptFile);

    }  # itrain loop	

    # use utility to write global DONE file when all of the
    # individual DONE files exist.
    my  $doneFile   = "$OUTPUT_TOPDIR/${mainScriptFile}.DONE" ;
    my ($waitArgs);

    $waitArgs = "$nfitopt  '*.DONE'  $doneFile" ;
    $NMAIN++; $CMDMAIN[$NMAIN] = "\n$WAIT_SCRIPT  $waitArgs" ;

    # -------------------------------------
    # finish main script for DMFIXRUN jobs.
  UPDATE_SCRIPT:
    &update_script($mainScriptFile, $DMFIXRUNdoneflag, @CMDMAIN );

} # end of make_DMFIXRUN_scripts


# ===========================
sub get_SALT2mu_commands(@) {

    my ($igroup,$itrain,$istage) = @_ ;

    # return string of commands for SALT2mu 

    my ($SALT2mu_cmd, @mktup_cmd, @clean_cmd, $mkroot_cmd );
    my ($program, $tmpArgs );
    my ($NCAT, $NCMD, @CMD, $tmpFile, $i, $prefix, $inFile );
    my ($MODEL, $fresFile_in);
    my ($whichcatlist);

    if ( $istage == $ISTAGE_SALT2mu ) {
	$whichcatlist = 1; 
    } elsif ( $istage == $ISTAGE_DMFIXPREP ) {
	$whichcatlist = 2;
    }

    $prefix      = &prefixSALT2mu($itrain, $istage);
    $inFile      = "${prefix}.input" ;   # SALT2mu input file
    $MODEL       = "$TRAIN_MODELNAME[$itrain]" ;
    $fresFile_in = "${MODEL}.FITRES" ;

    # SALT2mu command
    $program     = $PROGRAM_NAME[$ISTAGE_SALT2mu] ;
    $SALT2mu_cmd = "$program $inFile > ${prefix}.LOG" ;

    # create command to make ntuple from SALT2mu output
    my $h2rootLog = "h2root_${prefix}.log" ;
    my $mktupLog  = "MKTUP_${prefix}.LOG" ; 

    $program      =  $PROGRAM_COMBINE_FITRES ;
    $tmpArgs      = "" ;
    $tmpArgs      = "$tmpArgs  ${prefix}.fitres" ;
    $tmpArgs      = "$tmpArgs  -outprefix ${prefix}" ;
    $mktup_cmd[0] = "$program ${tmpArgs} > ${mktupLog}" ; 
    $mktup_cmd[1] = "mv ${prefix}.tup ${prefix}.TUP" ;
    $mkroot_cmd   = "h2root ${prefix}.TUP >&  $h2rootLog" ;

    # cleanup command
    $clean_cmd[0] = "rm ${prefix}.log" ;
    $clean_cmd[1] = "rm ${mktupLog} " ;
    
	# update commands in mini-script
    
    $NCMD = -1 ;

    # create command to catenate relevant fitres files
    # after the SIMFIT jobs have run	
    $NCAT = $NCATLIST_FITRES[$whichcatlist][$igroup][$itrain] ;
    $tmpFile = "";
    for ( $i = 1 ; $i <= $NCAT ; $i++ ) {
	$tmpFile = "$tmpFile $CATLIST_FITRES[$whichcatlist][$igroup][$itrain][$i]" ;
	##print "xxx CATFITRES: ${tmpFile} \n";
    }

    $NCMD++ ; $CMD[$NCMD] = "cat $tmpFile > $fresFile_in" ;   
    $NCMD++ ; $CMD[$NCMD] = "$SALT2mu_cmd"  ;
    $NCMD++ ; $CMD[$NCMD] = "$mktup_cmd[0]" ;	
    $NCMD++ ; $CMD[$NCMD] = "$mktup_cmd[1]" ; 
    $NCMD++ ; $CMD[$NCMD] = "$mkroot_cmd"   ; 
    $NCMD++ ; $CMD[$NCMD] = "$clean_cmd[0]" ;
    $NCMD++ ; $CMD[$NCMD] = "$clean_cmd[1]" ;
    $NCMD++ ; $CMD[$NCMD] = "touch ${prefix}.DONE" ;
    
    return @CMD ;

} # end of get_SALT2mu_commands

# ===========================

sub get_FITOPTKEYline(@) {

    # $verstype == 0: single VERSION, multiple FITOPTS
    # $verstype == 1: multiple VERSION, each with single FITOPT
    my ( $verstype, $idtrain, $version0) = @_;
    
    my ( $fitmodel, $key, $version, $c3);
    my ( $fitoptline);
	 
    $fitmodel = "SALT2.$TRAIN_MODELNAME[$idtrain]" ;
    $c3 = $TRAIN_CHARINDEX[$idtrain];
    $version = "${version0}_${c3}";
    
    if ( $verstype == 0 ) {
	$fitoptline = "FITOPT: FITMODEL_NAME $fitmodel";
    } else {
	$fitoptline = "VERSION+FITOPT: $version FITMODEL_NAME $fitmodel";
    }
    
    return $fitoptline;

} #end of get_FITOPTKEYline(@)


# ===========================

sub get_CATLISTinfo(@) {

    my ($nfitopt,$idtrain, $whichcatlist, $idgroup, $sdir, $fitdir) = @_ ;
    my ($ntmp, $f3res, $fresfile);

    # return $NCATLIST_FITRES, $CATLIST_FITRES for make_SALT2mu_scripts
    
    $ntmp = $NCATLIST_FITRES[$whichcatlist][$idgroup][$idtrain];
    $ntmp++ ;
    $f3res = sprintf("FITOPT%3.3d.FITRES", $nfitopt);
    $fresfile = "../../$sdir/$fitdir/$f3res" ;
    
    return ($fresfile, $ntmp);
    
} #end of get_CATLISTinfo

# ===========================
sub get_COSMO_commands(@) {

    my ($igroup,$itrain,$istage) = @_ ;

    # return string of commands for COSMO fit

    my ($COSMO_cmd, @mktup_cmd, @clean_cmd, $mkroot_cmd );
    my ($program, $tmpArgs, @temp_CPLIST );
    my ($NCAT, $NCMD, @CMD, $tmpFile, $i, $prefix, $inFile );
    my ($MODEL, $fresFile_in, $cpliststage);

    if ( $istage == $ISTAGE_DMFIXCOSMO ) {
	@temp_CPLIST = @CPLIST_CORR ;
	$cpliststage = 1 ;
    } else {
	@temp_CPLIST = @CPLIST_FITRES ;
	$cpliststage = $ISTAGE_SALT2mu ;
    }

    $prefix = &prefixSALT2mu($itrain, $istage);
    $fresFile_in = "${prefix}.fitres"; # changed for easier post-processing
    $prefix = "COSMO.${prefix}";
    $inFile = "${prefix}.input" ;   # COSMO input file
    $MODEL = "$TRAIN_MODELNAME[$itrain]" ;

    # COSMO command
    $program     = $PROGRAM_NAME[$ISTAGE_COSMO] ;
    if ($program =~ /wfit.exe$/ ) {
	$COSMO_cmd = "$program $fresFile_in $COSMO_GENOPT > ${prefix}.LOG" ;
    } else { # assume is sncosmo_mcmc.exe
	$COSMO_cmd = "$program $inFile > ${prefix}.LOG";
    }

    #update commands in mini-script

    $NCMD = -1 ;

    # create command to copy relevant fitres files
    # after the SALT2mu jobs have run
    $tmpFile = $temp_CPLIST[$cpliststage][$igroup][$itrain];

    $NCMD++ ; $CMD[$NCMD] = "cat $tmpFile > $fresFile_in" ;
    $NCMD++ ; $CMD[$NCMD] = "$COSMO_cmd";
    $NCMD++ ; $CMD[$NCMD] = "touch ${prefix}.DONE";

    return @CMD ;

} # end of get_COSMO_commands

# ===========================
sub get_DMFIXRUN_commands(@) {

    my ($igroup,$itrain,$istage) = @_ ;

    # return string of commands for DMFIXRUN fit

    my ($DMFIXRUN_cmd, @mktup_cmd, @clean_cmd, $mkroot_cmd );
    my ($script, $tmpArgs, $call1, $call2, $call3 );
    my ($NCAT, $NCMD, @CMD, $tmpFile, $i, $prefix, $inFile );
    my ($MODEL, $fresFile1_in, $fresFile2_in, $fresFile3_in, $GENOPT);
    my ($compgroup, $comppath, $compprefix, $compfile);
    my ($group, $topDir, $c3, $simalpha, $simbeta, $grepalpha, $grepbeta);

    $group    = $DMCORR_GROUPNAME[$igroup] ;
    $topDir   = "$STAGE_PATH[$istage]/$group" ;

    $prefix = &prefixSALT2mu($itrain, $ISTAGE_SALT2mu);
    $fresFile2_in = "${prefix}.fitres"; # test fitres file name
    $comppath = $STAGE_PATH[$ISTAGE_SALT2mu];
    $compgroup = $DATA_GROUPNAME[$DMCORR_IDEAL[$igroup]];
    $compfile = "${comppath}/${compgroup}/${prefix}.fitres";

    $prefix = &prefixSALT2mu($itrain, $istage);
    $fresFile1_in = "${prefix}.fitres"; # fix fitres file name

    $prefix = "DMFIXRUN.${prefix}";
    $inFile = "${prefix}.input" ;   # DMFIXRUN input file
    $MODEL = "$TRAIN_MODELNAME[$itrain]" ;
    $c3 = $TRAIN_CHARINDEX[$itrain];
    $GENOPT = "$DMCORR_RUNCODEOPT[$igroup]" ;
    

    $prefix = &prefixSALT2mu($itrain, $ISTAGE_DMFIXCOSMO);
    $CPLIST_CORR[1][$igroup][$itrain] = "${topDir}/${prefix}.fitres";
    $fresFile3_in = "${prefix}.fitres";

    if ( $DMCORR_SIMTYPE[$igroup] eq "AS-IS" ) {
	$grepalpha = sntools::setVarString($SHELL, "inalpha", "echo ${DMCORR_SIMALPHA[$igroup]}");
	$grepbeta = sntools::setVarString($SHELL, "inbeta", "echo ${DMCORR_SIMBETA[$igroup]}");
    } else {
	$grepalpha = sntools::setVarString($SHELL, "inalpha", "grep alpha ${fresFile2_in} | awk '{print \$4}'");
	$grepbeta = sntools::setVarString($SHELL, "inbeta", "grep beta ${fresFile2_in} | awk '{print \$4}'");
	$grepalpha =~ s/`/\\`/g; # need to escape ` for use with sntools::sedreplace function
	$grepalpha =~ s/\$/\\\$/g; # need to escape $ for use with sntools::sedreplace function
	$grepbeta =~ s/`/\\`/g;
	$grepbeta =~ s/\$/\\\$/g; # need to escape $ for use with sntools::sedreplace function
    }

	#print "grepalpha: $grepalpha \n";
	#print "grepbeta: $grepbeta \n";

    $simalpha = "\$inalpha";
    $simbeta = "\$inbeta";

    # DMFIXRUN command
    if ( $DMCORR_RUNCODE[$igroup] eq "DEFAULT" ) {
	$script =  "${SCRIPT_PREFIX[$ISTAGE_DMFIXRUN]}.pl" ;
	$call1 = "-bcorrfile $fresFile1_in -datfile $fresFile2_in ";
	$call2 = "-outfile $fresFile3_in -corrtype $DMCORR_TYPE[$igroup] ";
	$call3 = "-sima ${simalpha} -simb ${simbeta} ";
    } elsif ( $DMCORR_RUNCODE[$igroup] =~ /SALT2train_biascorr.pl$/ ) {
	$script = $DMCORR_RUNCODE[$igroup] ;
	$call1 = "-bcorrfile $fresFile1_in -datfile $fresFile2_in ";
	$call2 = "-outfile $fresFile3_in -corrtype $DMCORR_TYPE[$igroup] ";
	$call3 = "-sima ${simalpha} -simb ${simbeta} ";
    } else {
	$script = $DMCORR_RUNCODE[$igroup] ;
	$call1 = "" ;
	$call2 = "" ;
	$call3 = "" ;
    }

        $DMFIXRUN_cmd = "$script $call1 $call2 $call3 $GENOPT > ${prefix}.LOG";
	$DMFIXRUN_cmd =~ s/\$/\\\$/g; # need to escape $ for use with sntools::sedreplace function


    #update commands in mini-script

    $NCMD = -1 ;

    # create command to copy relevant fitres files
    # after the SALT2mu jobs have run
    $tmpFile = $CPLIST_FITRES[$ISTAGE_DMFIXPREP][$igroup][$itrain];

    $NCMD++ ; $CMD[$NCMD] = "cat $tmpFile > $fresFile1_in" ;
    $NCMD++ ; $CMD[$NCMD] = "cat $compfile > $fresFile2_in" ;
    $NCMD++ ; $CMD[$NCMD] = ${grepalpha} ;
    $NCMD++ ; $CMD[$NCMD] = ${grepbeta} ;
    ##OBSOLETE $NCMD++ ; $CMD[$NCMD] = "echo $fresFile1_in $fresFile2_in $fresFile3_in > $inFile" ;
    $NCMD++ ; $CMD[$NCMD] = "$DMFIXRUN_cmd";
    $NCMD++ ; $CMD[$NCMD] = "touch ${prefix}.DONE";

    return @CMD ;

} # end of get_DMFIXRUN_commands


# ===========================================
sub make_PLOT_scripts {
    
    # for this igroup, construct script to make PLOTs
    # prepare macro(s) with following options:
    # * $prefix for each training
    # * $prefixList to summary all trainings

    my($igroup) = @_;

    
    if ( $N_PAWPLOTMACRO == 0  )  { return ; }
 
    my ($icmd,@CMD,  $itrain, $ntupLink, $itmp );
    my ( $macro, $ipaw, $MACRO_NAME, $MACRO_ARG, $macro_name );
    my ($prefix, $prefixList, $DMP_FILE );
    my ($psFile, $psFile2, $pdfFile, $pscmd, $kfile, $ntupFile);

    my $istage       = $ISTAGE_PLOTS ;
    my $group        = "$DATA_GROUPNAME[$igroup]" ;
    my $topDir       = "$STAGE_PATH[$istage]/$group" ;
    my $ntupDir      = "$STAGE_PATH[$ISTAGE_SALT2mu]/$group" ;
    my $ntup_subDir  = "$STAGE_SDIR[$ISTAGE_SALT2mu]/$group" ;
    my $OPT_SCRIPT   = 1 ;
    my $scriptFile   = &scriptFileName($istage,$igroup,0,$OPT_SCRIPT);
    my $nfitopt;

    $nfitopt = $NMODEL_TRAIN; #default, unless SIMFLAG & SIMDATA_EXCLUDE_FITOPT[$igroup]==1
    if ($SIMFLAG) {
	if ($SIMDATA_EXCLUDE_FITOPT[$igroup]) {$nfitopt=1;}
    } 


    # create subdir for this SALT2mu-group
    qx(mkdir $topDir);

    $icmd = -1 ;
    @CMD = ( ) ;
    $icmd++ ; $CMD[$icmd] = "cd $topDir" ;

    # create symbolic link to each SALT2mu-ntuple
    # and create prefixList

    for ( $itrain=1; $itrain <= $nfitopt ; $itrain++ ) {

	$prefix      = &prefixSALT2mu($itrain, $istage);
	$prefixList  = "$prefixList $prefix";

	# make symbolic link to relative path of SALT2mu ntuple
	# so that the link still works if the topDir is moved.
	$ntupFile    = "${prefix}.TUP" ;
	$ntupLink = "ln -s ../../${ntup_subDir}/${ntupFile}  ${ntupFile}" ;
	
	# update commands in mini-script

	$icmd++ ; $CMD[$icmd] = "$ntupLink" ;
    } # itrain
    $icmd++ ; $CMD[$icmd] = " " ;

    # call prefix-dependent  macros
    for ( $ipaw=0; $ipaw < $N_PAWPLOTMACRO ; $ipaw++ ) {
	$MACRO_NAME =  $PAWPLOTMACRO_NAME[$ipaw] ;
	$MACRO_ARG  =  $PAWPLOTMACRO_ARG[$ipaw] ;
	if ( $MACRO_ARG ne "prefix" ) { next ; }

	for ( $itrain=1; $itrain <= $nfitopt ; $itrain++ ) {
	
	    $prefix      = &prefixSALT2mu($itrain, $istage);
	    $ntupFile    = "${prefix}.TUP" ;
	
	    # create paw-kumac file to make plots using user-defined macro
	    $macro   = "exec $MACRO_NAME ${prefix}" ;
	    $psFile  = "${prefix}.ps" ;
	    $pdfFile = "${prefix}.pdf" ;
	    $kfile   = "${prefix}.kumac" ;
	    qx ( echo  "$macro"     >  ${topDir}/${kfile} );
	    
	    # update commands in mini-script
	    $DMP_FILE = "paw_${prefix}.log"; 
	    $icmd++ ; $CMD[$icmd] = "$PAW_COMMAND -b $kfile > $DMP_FILE";
	    $icmd++ ; $CMD[$icmd] = "rm $DMP_FILE";
	    $icmd++ ; $CMD[$icmd] = "ps2pdf $psFile \n";
	}  # itrain
    } # ipaw


    # call prefixList-dependent macro(s) to give summary
    # of all trainings for this group.

    for ( $ipaw=0; $ipaw < $N_PAWPLOTMACRO ; $ipaw++ ) {

	$MACRO_NAME =  $PAWPLOTMACRO_NAME[$ipaw] ;
	$MACRO_ARG  =  $PAWPLOTMACRO_ARG[$ipaw] ;
	if ( $MACRO_ARG ne "prefixList" ) { next ; }

	$macro   = "exec $MACRO_NAME ${prefixList}" ;
	$kfile   = "${group}.kumac" ;
	qx ( echo  "$macro"     >  ${topDir}/${kfile} );
	$icmd++ ; $CMD[$icmd] = "$PAW_COMMAND -b $kfile > DMP.ALL" ;
	$icmd++ ; $CMD[$icmd] = "rm DMP.ALL " ;

	# convert psFile into pdf file; assume psFile name is
	# the name of the macro
	$itmp        = index($MACRO_NAME,"#");
	$macro_name  = substr($MACRO_NAME,$itmp+1,100);
	$macro_name  =~ s/\s+$// ;  # remove trailing blanks
	$psFile      = "${macro_name}.ps";
	$psFile2     = "${macro_name}_${group}.ps";

	$icmd++ ; $CMD[$icmd] = "mv $psFile $psFile2 \n";  
	$icmd++ ; $CMD[$icmd] = "ps2pdf $psFile2 \n"; 
    }

	# -------------------------
    &update_script($scriptFile, 1, @CMD );

} # end of make_PLOT_scripts {


# ===============================
sub make_CLEAN_script {


    my ($istage, $iscript, $scriptFile, @CMD, $model);
    my ($itrain, $prefix, $icmd, $NPURGE, $purge );
    my $OPT_SCRIPT = 1;

    $istage = $ISTAGE_CLEAN ;  $iscript=0;
    $scriptFile = &scriptFileName($istage, $iscript, 0, $OPT_SCRIPT );

    if ( $OPT_CLEAN == 0 ) { 
	@CMD[0] = " ";
	goto UPDATE_SCRIPT ;
    }

    $NPURGE = scalar(@FILELIST_PURGE);
    $icmd = -1;    
    for ( $itrain = 1; $itrain <= $NMODEL_TRAIN; $itrain++ ) {

	$model = $TRAIN_MODELNAME[$itrain] ;

	$icmd++ ;
	@CMD[$icmd] = "cd $OUTPUT_TOPDIR/AAA_PCAFIT_TOPDIR/$model" ; 

	foreach $purge ( @FILELIST_PURGE ) {
	    $purge  =~ s/\s+$// ;  # remove trailing blanks
	    $icmd++ ;  @CMD[$icmd] = "  rm $purge";
	}

	# gzip remaining files.
	$icmd++ ;  @CMD[$icmd] = "  gzip *.*";	
	$icmd++ ;  @CMD[$icmd] = "  ";	
    }

  UPDATE_SCRIPT:
    &update_script($scriptFile, 1, @CMD );

} # end of  make_CLEAN_script
 
# ===============================
sub make_ARCHIVE_script {

    # Create script to make tarball, and use EXCLUDE_FILE 
    # to exclude large files (fits and his)

    my ($istage, $iscript, $scriptFile, @CMD, $icmd, $PATH, @tmp, $N );
    my ($topDir, $sDir, $tarName, $XARG, $igroup, $igen );
    my (@GENVERSION_LOCAL_LIST, $GENVERSION, $GG);
    my $OPT_SCRIPT = 1;
    my $EXCLUDE_FILE = "EXCLUDE_FROM_TAR.DAT";

    $istage = $ISTAGE_ARCHIVE ;  $iscript=0;
    $scriptFile = &scriptFileName($istage, $iscript, 0, $OPT_SCRIPT );
    $PATH = "$STAGE_PATH[$istage]" ;

    if ( $OPT_ARCHIVE == 0 ) { 
	@CMD[0] = " ";
	goto UPDATE_SCRIPT ;
    }

    # get name of top directory after last forward slash;
    # this topDir name is used to construct the tarball name
    @tmp     = split('/',$OUTPUT_TOPDIR) ;  $N = scalar(@tmp);
    $topDir  = $tmp[$N-1];

    # start with EXCLUDE_FILE

    open  PTR_XFILE , "> $PATH/$EXCLUDE_FILE" ;

    # exclude fits files created by training
    $sDir = "$STAGE_SDIR[$ISTAGE_TRAINRUN]";
    print PTR_XFILE "$topDir/$sDir/workSpace/*/*.fits* \n";

    # exclude SPLIT_JOBS directories from split_and_fit
    $sDir = "$STAGE_SDIR[$ISTAGE_SIMFIT]";
    print PTR_XFILE "$topDir/$sDir/TEMP*/SPLIT_JOBS \n" ;

    # exclude his.gz files in each group
    for ( $igroup = 1; $igroup <= $NGROUP_DATA; $igroup++ ) {    
	@GENVERSION_LOCAL_LIST    = split(/\s+/,$SIMDATA_GENVERSION[$igroup]) ;
	foreach $GENVERSION  ( @GENVERSION_LOCAL_LIST ) {
	    $sDir = "$STAGE_SDIR[$ISTAGE_SIMFIT]" ;
	    $GG   = "TEMP_$GENVERSION/$GENVERSION" ;
	    print PTR_XFILE "$topDir/$sDir/$GG/FITOPT*.his* \n" ;
	}
    }

    close PTR_XFILE ;

    $tarName = "${topDir}.tar" ;

    $XARG = "-X $topDir/$STAGE_SDIR[$ISTAGE_ARCHIVE]/$EXCLUDE_FILE" ;

    # construct tar commands 
    $icmd = -1;
    $icmd++ ; $CMD[$icmd] = "cd ../../";
    $icmd++ ; $CMD[$icmd] = "tar -cf $tarName $XARG $topDir";
    $icmd++ ; $CMD[$icmd] = "gzip $tarName" ;

  UPDATE_SCRIPT:
    &update_script($scriptFile, 1, @CMD );

} # end of  make_ARCHIVE_script

# ===============================
sub prefixSALT2mu(@) {

    my ($itrain, $istage) = @_;
    my ($c3,$MODEL,$prefix,$type);

#    print "xxx in prefixSALT2mu. itrain is $itrain , istage is $istage \n";
    if ($istage == $ISTAGE_SALT2mu || $istage == $ISTAGE_COSMO || $istage == $ISTAGE_PLOTS ){
	$type = "TEST";
    } elsif ($istage == $ISTAGE_DMFIXRUN || $istage == $ISTAGE_DMFIXPREP){
	$type = "FIX";
    } elsif ($istage == $ISTAGE_DMFIXCOSMO){
	$type = "CORR";
    } else {
	$MSGERR[0] = "subroutine prefixSALT2mu doesn't recognize istage value $istage ";
	$MSGERR[1] = "Allowed are SALT2mu $ISTAGE_SALT2mu, COSMO $ISTAGE_COSMO, plots $ISTAGE_PLOTS ";
	$MSGERR[2] = "Allowed are DMFIXRUN $ISTAGE_DMFIXRUN, DMFIXPREP $ISTAGE_DMFIXPREP ";
	$MSGERR[3] = "and DMFIXCOSMO $ISTAGE_DMFIXCOSMO ";
	$MSGERR[4] = "Check inputs and try again. ";
	sntools::FATAL_ERROR(@MSGERR);	    	
    }

    $c3          = sprintf("%3.3d", $itrain);
    $MODEL       = "$TRAIN_MODELNAME[$itrain]" ;
    $prefix      = "SALT2mu${c3}.${MODEL}.${type}";

    return $prefix ;

} # end of prefixSALT2mu

# =====================================
sub KILL_pipeline {

    my $q = "\"";
    my ($inode, $node, $NTMP, $cmdkill, $i ) ;
    my (@NODELIST_SPLIT);

    $NTMP = scalar(@NODELIST);    
    for ( $i=0; $i < $NTMP ; $i++ ) {
	@NODELIST_SPLIT    = split(/\s+/,$NODELIST[$i]) ;
	foreach $node ( @NODELIST_SPLIT ) {	    
	    $cmdkill = "ssh -x $node ${q}kill -KILL -1${q}" ;
	    print "\t $cmdkill \n";
	    system("$cmdkill &");
	}
    }

    die "\n Done killing jobs. \n";

} # end of KILL_pipeline


# ====================================
sub RUN_pipeline {

    # submit ./RUN# scripts in background sequentially.
    # After each script is submitted it waits for the DONE-STAMP.

    my ($iscript, $script, $istage, $cdate);
    my ($STAGELOG, $NSCRIPT, @scriptList, $PATH, $CTMP, @tmp );
    my $dashLine = " ----------------------------------------------------- \n";

    print " Launch RUN-scripts in \n\t $OUTPUT_TOPDIR \n";

    # shorten sleep times for QUICKTEST
    my $tmpFile = "$OUTPUT_TOPDIR/QUICKTEST" ;
    if ( -e $tmpFile ) {
	$STAGE_SLEEPCHECK[$ISTAGE_TRAINRUN] =  10 ;
	$STAGE_SLEEPCHECK[$ISTAGE_SIMFIT]   =  20 ;
	$STAGE_SLEEPCHECK[$ISTAGE_SALT2mu]  =  10 ;
	$STAGE_SLEEPCHECK[$ISTAGE_CLEAN]    =  10 ;
	$OPT_CLEAN = 0 ;
    }

    # Feb 3, 2012:
    # get list of RUN* scripts, but weed out the RUN*DONE files
    # which may exist on a RESTART
    @tmp = qx(cd $OUTPUT_TOPDIR ; ls RUN*);
    @scriptList = ();
    foreach $script ( @tmp ) {
	if  ( index($script,"DONE") > 0 ) { next ; }
        @scriptList = ( @scriptList , $script );
    }
    $NSCRIPT = scalar(@scriptList);


    # open STAGE-LOG as new file or in append mode
    $STAGELOG = "$OUTPUT_TOPDIR/STAGE.LOG" ;
    if ( $ISTAGE_START == 1 ) {
	if ( -e $STAGELOG ) { qx(rm $STAGELOG); }
	open  PTR_STAGELOG  , "> $STAGELOG" ;  # new file 
    }
    else {
	open  PTR_STAGELOG  , ">> $STAGELOG" ;  # append mode for RESTART
	$CTMP = "ISTAGE = $ISTAGE_START ($STAGE_NAME[$ISTAGE_START])"; 
	$cdate  = `date` ;  $cdate  =~ s/\s+$// ; 
	print PTR_STAGELOG "\n";
	print PTR_STAGELOG "# ========================================= \n";
	print PTR_STAGELOG " RESTART at $CTMP \n";
	print PTR_STAGELOG " $cdate \n\n";
    }
    PTR_STAGELOG->autoflush(1); 

    for ( $iscript=0; $iscript < $NSCRIPT ; $iscript++ ) {

	$script  = $scriptList[$iscript] ;
	$script  =~ s/\s+$// ;   # trim trailing whitespace

	# translate name of script into stage number
	$istage  =  &I2STAGE_extract($script);
	$PATH    = $STAGE_PATH[$istage] ;
       
	# make sure previous done files exists and current one does NOT
	&check_Done($script,$istage); 

	# skip stages already done (to allow for RESTART)
	if ( $istage < $ISTAGE_START )  {  next ; }

	# Apr 26 2013: apply option to skip later stages
	if ( $istage > $ISTAGE_STOP  )   { next ; }

	# ---------------------
	# launch this stage

	$cdate  = `date` ;  $cdate  =~ s/\s+$// ; # trim trailing whitespace
	print " exec ./$script  (stage $istage) \n" ;
	print PTR_STAGELOG "$dashLine \n";
	print PTR_STAGELOG " launch ./$script  ($cdate) \n" ;
	system("cd $OUTPUT_TOPDIR ; ./$script & ");
	# -----------------------

	# wait for this stage to finish
	wait_for_Done($script, $istage);

	# set group read+write priv on STAGE dir (Jan 19, 2012)
	print PTR_STAGELOG "\t Set g+rw priv in $PATH \n" ;
	qx(chmod -R g+rw $PATH) ;

    }

    close  PTR_STAGELOG ;

    die "\n Done. \n";

} # end of RUN_pipeline


# ========================================
sub I2STAGE_extract {

    my($runScript) = @_ ;

    # return 2-digit integer $istage from the name of the
    # input runScript. Example: if $runScript = RUN04_XYZ then
    # then 04 is returned

    my $istage  = substr($runScript,3,2);
    
    return $istage ;

} # end of I2STAGE

# ====================
sub check_parallel_batch(@) {

    my ($stagename, $jobs, $nodes) = (@_);

    if ( $jobs > $nodes ) {
	$MSGERR[0] = "Number of jobs in $stagename";
	$MSGERR[1] = "exceeds number of batch nodes";
	$MSGERR[2] = "specified by BATCH_INFO_TRAIN";
	$MSGERR[3] = "Check inputs and try again";
	sntools::FATAL_ERROR(@MSGERR);
    }   # else {
	#print "\n  xxx check_parallel_batch for $stagename \n";
	#print "\n  xxx nodes: $nodes, jobs: $jobs \n";
        #}
       
} # end of check_parallel_batch(@)

# ====================
sub wait_for_Done {

    my ($script,$istage) = @_;

    my ($t_wait, $cdate, $DONESTAMP, $DONEFILE ) ;    

    $t_wait    =  $STAGE_SLEEPCHECK[$istage] ;

    $DONESTAMP = "${script}.DONE" ;
    $DONEFILE  = "$OUTPUT_TOPDIR/${script}.DONE" ;

  CHECKDONE:
    $cdate  = `date` ;  $cdate  =~ s/\s+$// ;   # trim trailing whitespace
    if ( -e $DONEFILE ) { 
	print PTR_STAGELOG "\t Found DONE-STAMP '${DONESTAMP}' ($cdate) \n" ;
	return ;
    }
    else {
	print PTR_STAGELOG "\t Still checking for $DONESTAMP  ($cdate) \n" ;
	sleep($t_wait);
	goto CHECKDONE ;
    }


} # end of wait_for_Done


sub check_Done {

    my ($script,$istage) = @_ ;

    # Feb 03, 2012
    # if $istage < $ISTAGE_START then require DONE file to exist.
    # Otherwise require DONE file to not exist.
    # ABORT if either condition fails.

    my ($DONEFILE, $LDONE, $cstage, $CSTAGE_START, $name, $I2);

    $DONEFILE  = "$OUTPUT_TOPDIR/${script}.DONE" ;
    $LDONE     =  (-e $DONEFILE) ;

    $name        = $STAGE_NAME[$istage] ;
    $I2          = sprintf("%2.2d", $istage);
    $cstage      = "stage = $I2 ($name)";

    $name         = $STAGE_NAME[$ISTAGE_START] ;
    $I2           = sprintf("%2.2d", $ISTAGE_START);
    $CSTAGE_START = "ISTAGE_START = $I2 ($name)" ;

    if ( ($istage < $ISTAGE_START)  && ($LDONE == 0) ) {
	$MSGERR[0] = "DONE-file does NOT exist for $cstage";
	$MSGERR[1] = "Cannot begin at $CSTAGE_START";
	sntools::FATAL_ERROR(@MSGERR);	    	
    }

    if ( ($istage >= $ISTAGE_START)  && $LDONE ) {
	$MSGERR[0] = "DONE-file exists for $cstage";
	$MSGERR[1] = "Cannot begin at $CSTAGE_START";
	sntools::FATAL_ERROR(@MSGERR) ;	    	
    }

} # end of check_Done

