#!/usr/bin/perl
#
# Created Jan 2011 by R.Kessler
#
# Generic script to fit all SN light curves efficiently by 
# splitting the sample into NSPLIT = N_NODE jobs and submiting 
# the jobs in parallel to each node. If no nodes are given, 
# the jobs are run serially in interactive mode. The input 
# file is the usual fit-namelist file, but with extra keys 
# outside the namelists as follows:
#
#  REQUIRED KEYS at TOP of namelist file (outside &SNLCINP):
#
#   OUTDIR:  <output directory for his, fitres, log files>
#
#      BEWARE that current OUTDIR is deleted, then re-created
#             ==> do NOT use existing directory that you want to keep.
#
#  OPTIONAL KEYS at top of namelist file:
#
#   NODELIST: <node1> <node2> ... <nodeN>
#    (no NODELIST => run interactively on current node)
#               or
#   BATCH_INFO:  <command>  <templateFile>  <Ncore>
#     (i.e,, on torque system, <command> = qsub)
#
#   BATCH_MEM:  2000    ! optional override default memory allocation
#
#   NJOB_PER_BATCHFILE: <NJOB>     ! default NJOB=1
#     (NJOB fit jobs per batch submission. Example: 5000 fit jobs,
#      Ncore=50, and NJOB_PER_BATCHFILE=20 --> 5000/20 = 250
#      submitted batch jobs instead of 5000 batch jobs)
#
#   DONE_STAMP:  ALL.DONE"    
#      (global done file to help higher-level scripts)
#
#   SNANA_LOGIN_SETUP:  <whatever command(s) to setup snana>
#     (use this key only if your login does NOT setup snana)
#
#   SALT2mu_INFILE:  <inFile>
#    (run SALT2mu using this input-file as template)
#
#   SALT2mu_SIMVERSION_INPUT:  <version>
#       (use this sim-version for Ia biasCorr and CC prior)
#   SALT2mu_SIMVERSION_INPUT:  2
#       (idem, but use 2nd VERSION in list)
#
#   SALT2mu_BIASCOR_PATH: <path>
#       Pick up sim-FITRES files here for biasCor. e.g., to fit 
#       FITOP002, get bias from FITOPT002.FITRES in SALT2mu_BIASCOR_PATH.
#
#   SALT2mu_CCPRIOR_PATH: <path>
#       Pick up sim-FITRES files here for CCprior
#         [if not given, use BIASCOR_PATH for CCprior]
#       e.g., to fit FITOP003, get CCprior from FITOP003.FITRES in 
#          SALT2mu_CCPRIOR_PATH
#
#
#   DO_FITOPT000:  0        # do NOT process FITOPT000
#
#   DELAY_SUBMIT:   0.5     # .5 sec delay between job-submits
#
#   DELAY_SUBMIT:   0.5,10   # .5 sec delay between job-submits
#                            # and 10 sec between VERSION/FITOPT
#
#   GZIP_FLAG: 1    (zip SPLIT_JOBS)
#   GZIP_FLAG: 3    (zip SPLIT_JOBS and merged hbook files)
#       (gzip output; default=0)
#   GZIP_FLAG: 4    (tar -cf SPLIT_JOBS.tar SPLIT_JOBS, then gzip it)
#
#   H2ROOT_FLAG: 1
#       (convert his files into root files; default=0)
#   H2ROOT_FLAG: 3
#       (idem, but leave log files for debugging)
#
#  APPEND_TABLE_TEXT:  <var1>  <var2> ...
#     (extract variables from hbook or root file and append text table)
#
#  APPEND_FITRES: FITRES  <var1>  <var2> ... 
#     (OBSOLETE: please switch to APPEND_TABLE_TEXT key above)
#
#  MIN_SNANA_VERSION: <version>
#     (abort if using earlier snana version)
#
#  VERSION_AFTERBURNER: <aribtrary user command to run inside each /VERSION subdir>
#
#   PLOTOPT: 
#       run mkfitsplots.pl with defaults
#   PLOTOPT: <options>
#       run mkfitsplots.pl with options
#
#   VERSION: VV1
#   VERSION: VV2
#   VERSION: VV3
#   VERSION: [nickName] VV3
#    [nickName not used, but passed to FITOPT.README for external programs]
#   VERSION: VV4_TEST*   ! use wildcard to select multiple verions
#   ..
#
#   FITOPT: FFF1  
#   FITOPT: FFF2
#   FITOPT: ->FITOPT000  (leave symb link for combining systematics)
#   ...
#
#      or
#   VERSION+FITOPT: VV1   <fit options 1>
#   VERSION+FITOPT: VV2   <fit options 2>
#     (warning: do NOT mix VERSION  and VERSION+FITOPT keys)
#
#   FITRES_COMBINE_FILE:  xyz.fitres  
#     (include this external file in all combine_fitres.exe commands)
#
#   VERSION_FITRES_COMBINE: <any-version-below>
#     (combine each merged fitres file with fitres from this version)
#
#   COMBINE_ALL_FLAG: 1  # default = 0
#     (combine/sum all fitres files into one; useful for large sim jobs)
#
#   TOTAL_WAIT_ABORT:  48  # increase default 24 abort time to 48 hr
#
#   The FITOPT arguments are command-line overrides 
#   used to perform systematic tests and variations.
#   Each VERSION is fit with each FITOPT plus the
#   the default FITOP000 with no extra fit-options
#   (i.e., using exactly what is in the namelist file).
#   If no "FITOPT: <FITOPT>" keys are given, then only 
#   the default FITOP000 is run. If no "VERSION: <VERSION>" 
#   are given, then the version from the VERSION_PHOTOMETRY 
#   argument is processed.
#
# Usage
#  split_and_fit.pl <inputFile>         # submit jobs
#
# Usage OPTIONS:
#  split_and_fit.pl <inputFile> MERGE   # merge root,hbook,text files
#  split_and_fit.pl <inputFile> noMERGE # run jobs without MERGE
#  split_and_fit.pl <inputFile> KILL    # kill all jobs
#  split_and_fit.pl <inputFile> HRECOVER [PREFIX] [FITOPT]  # hmerge recovery
#  split_and_fit.pl <inputFile> NOPROMPT  # do NOT check multiple jobs
#  split_and_fit.pl <inputFile> NOLAUNCH  # prepare, but do NOT launch 
#  split_and_fit.pl <inputFile> NOSUBMIT  # same
#  split_and_fit.pl <inputFile> 1         # 1 second delay between job-submits
#  split_and_fit.pl <inputFile> TRAILMARK # leave note in sim-README (not for data)
#  split_and_fit.pl <inputFile> kicp   # force kicp queue
#  split_and_fit.pl <inputFile> KICP   # idem
#
#  split_and_fit.pl CLEAN            # cleanup all sub-directories
#  split_and_fit.pl CLEAN NOPROMPT   # idem, but no user-authorization
#  split_and_fit.pl CLEANMASK <mask> # see manual for selective cleaning
#                                    # or search CLEANMASK_xxx below
#
#
#  The MERGE call is made automatically by the first call above;
#  use the MERGE option only for debugging.
#  
#  TRAILMARK will write the split_and_fit OUTDIR in the simulated README
#  file so that later split_and_fit jobs finding a symbolic link in a
#  simDir will be able to trace the original split_and_fit OUTDIR. 
#  Enabled only for simulated README (never write to dataDir).
#
#
#              HISTORY
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ---------------------------------------------------------------
#
# Sep 12, 2014: 
#   major overhaul to use SNTABLE_LIST and to remove legacy nml variables
#   after substituting the equivalent SNTABLE_LIST string. Therefore
#   the legacy variables can still be used, but they are replaced
#   internally here.
#
# Oct 28 2014: combine_fitres.txt -> combine_fitre.text 
#              after refactor. Similar change for appending
#              with sntable_dump.pl.
#
# Nov 11 2014: modify to work with data files in lcmerge/ or lcmerge/VERSION.
#
# Dec 13 2014: new function set_OUTDIR() to strip quote and check
#              for sub-dir option if there are no slashes in path.
#
# Dec 26 2014: new GZIP_FLAG option (=4) to make SPLIT_JOBS.tar.gz
# 
# Dec 28 2014: in merge_outFiles, if merge fails try again.
#
# Jan 2 2015:
#   + fix DONE_STAMP to appear after [optional] SALT2mu jobs.
#   + create SALT2mu_JOBS/NJOBS.README with exact count instead
#     of using ls ORIG*
#
# Jan 8 2015: increase BATCH_MEM to 2GB (was 500MB)
# 
# Feb 10 2015: if PRIVATE_DATA_PATH = 'pwd', substitute current dir
#              into nml file. 
#
# Mar 6 2015: add FUDGE_HOSTNOISE_FILE to the copy-list
#
# Apr 6 2015: $DELAY_SUBMIT --> @DELAY_SUBMIT[2]
#             between jobs, and between VERSION/FITOPT
#
# Jul 19 2015: new input key "DO_FITOPT000: 0" to skip FITOPT000
#              and only process the FITOPT options.
#
# Aug 23 2015: if executable is in current path, tack on $CURRENT_DIR
#
# Aug 24 2015: new optional key "MIN_SNANA_VERSION: <version>".
#              Split_and_fit aborts if using lower version.
#
# Aug 31 2015: new key "VERSION_AFTERBURNER: <any user command>"
#
# Feb 19 2016: new key  SALT2mu_SIMVERSION_INPUT  to use simFile as
#              Ia bias corr and CC prior.
#
# Mar 3 2016: make merge-logic wait more robust by making sure that 
#             MERGE-DONE-FILE exists before declaring done.
#             See new MERGE_STATUS[$iver][$ifitopt] = 1 (start merge)
#             or =2 when merging is done.
#
# Apr 7 2016: new key SALT2mu_BIASCOR_PATH
#
# Jun  6 2016: default GZIP_FLAG=4 to create & gzip SPLIT_JOBS.tar
#
# Jun 10 2016: gzip SIMGEN.DAT file (for simulations)
#
# Jun 20 2016: remove checks on legacy keys (LTUP_xxx) since the
#               fitting codes will abort on them.
#
# July 12 2016: allow wild card for simulated versions; e.g., 
#                 VERSION: DESTEST_dataSet-*
#               See new function parse_VERSION().
#
# Aug 29 2016; 
#    + new input key "WAIT_CHECK_DONE: xx" (default=20 sec)
#      to reduce time between checking status
#    + if NSPLIT==1, use 'cp' to copy instead of mergeJob.
#         [to speed up Elise's ABC]
#
# Sep 10 2016: new FITOPT option,
#           FITOPT: ->FITOPT000
#      to just make symbolic link to output of FITOPT000.
#      Useful to synchronize systematics with multiple surveys.
#     
#      Add VPEC_FILE to list of COPY_INPUTFILE_KEYS
# 
# Jan 11 2017: 
#     + add HEADER_OVERRIDE_FILE to COPY_INPUTFILE_KEYS
#
# Jan 24 2017: allow ENV in OUTDIR name.
#
# Jan 30 2017: 
#   + print elapsed time in MERGE.LOG
#   + abort after 24 hr. Can change value with key "TOTAL_WAIT_ABORT: 48"
#
# Feb 2 2017: add MAGCOR_FILE to @COPY_INPUTFILE_KEYS
#
# Feb 8 2017:
#  + remove NN_JOBSUDIR, NN_FLAG and NN_xxx 
#  + replace NN-specific option with more general MLTRAIN and MLAPPLY
#  + SPLIT_JOBS -> SPLIT_JOBS_LCFIT/  SPLIT_JOBS_ML/  SPLIT_JOBS_SALT2mu
#  + some refactoring to handle ML and SALT2mu postProc jobs
#  + new command-line arg TRAILMARK (see above)
#
# Mar 2 2017: define BATCH_MEM_LCFIT[SALT2mu,MLAPPLY]. MLAPPLY needs
#             more memory for the COMBINE-fitres command.
#
# Apr 17 2017: add NONLINEARITY_FILE to  @COPY_INPUTFILE_KEYS
#
# May 15 2017: new feature to allow [nickName] following each FITOPT.
#              nickName is not used here, but passed on to FITOPT.README
#              for external programs.
#
# Jul 13 2017:  
#   + increase memory for MLAPLLY 2GB -> 4GB
#   + increase memory for SALT2mu 1GB -> 2GB
#
# July 17 2017:
#  + if $NBATCH_CREATE > BATCH_NCORE, then add a 'wait_for_files'
#    command to wait for enough DONE files so that number of jobs
#    in the queues does not exceed user-requested BATCH_NCORE
#    See new function batch_delayAdd().
#
# Nov 8 2017: 
#  + new command-line argument 'KICP' or 'kicp' to force using kicp
#    queue (instead of sandyb). Works on Midwawy only.
#
# Nov 16 2017:
#  + Fix DONE-logic bug when NJOB > NCORE and there are FITOPT links.
#    See variable $NDONE_REQ argument to wait_fot_files.
#  + Speed merging in DONE_CHECK_LCFIT() by sleep-waiting only
#    if ifitopt==0
#
# Nov 21 2017:
#  + fix so that VERSION wildcards work with user-defined PATH_SNDATA_SIMDIR
#
# Nov 27 2017
#  + update CLEANFILES_DRIVER to use CLEANMASK to selectively clean.
#    See CLEAN and CLEAMASK arguments in manual, or CLEANMAS_XXX below.
#
# Dec 13 2017:
#  +  in mergeLog_update(), cat in groups of 100 versions to prevent
#     error from 'too-long list'.
#  + In MERGE functions, if NSPLIT==1, use exact string match instead
#    of wild card --> much faster when NJOB >> NCORE.
#
# Dec 20 2017:
#   + new key "NJOB_PER_BATCHFILE: <NJOB>"  (default is 1) 
#     to combine <NJOB> jobs in each batch file submission.
#     Goal is to reduce number of short batch jobs that have lots of overhead.
#   + BATCH file names are now FITSPLIT_[ibatch].BATCH, 
#     and same for BUSY file prefix
#
# Feb 6 2018:
#  new input key  "BATCH_MEM:  <mem>"  to override 1000 MB default 
#
# May 16 2018: if snlc_fit.exe (or psnid) does not exist, 
#              wait 200 sec.  Hopefully works during 'make'.
#
# Sep 11 2018:
#   + if NJOB_PER_BATCHFILE is specified, then increase
#     $WAIT_CHECK_DONE -> 120 (intead of 20 sec)
#
#   + in mergeLog_cat(), trap cat errors to suppress stdout.
#
# Jan 16 2019: make sure CLEAN option return proper exit code
# Feb 18 2019: for SALT2mu option, mv .M0DIF file to ../[VERSION]
#
# Mar 07 2019: 
#   + in parse_GENVERSION(), abort on missing version
#   + fix to work with sims in differnt paths
#   + add time-stamp check (CDATE_XXX) and abort on a rogue job.
#
# Mar 25 2019
#   + fix merging to work with SNANA table; see merge_text_TABLE().
#   + new "APPEND_TABLE_TEXT: <varList>" to replace APPEND_FITRES key.
#     The new key does not need a table-name specifier.
#
# Apr 17 2019:
#   + echo SUCCESS or FAILURE to DONE_STAMP file.
#   + replace FATAL_ERROR with FATAL_ERROR_STAMP to write
#     FAILURE in done-stamp.
#
# Jun 4 2019: "TOTAL_WAIT_ABORT -> 48 hr (was 24 hr)
#
# =======================

use List::Util qw(first);
use File::Path qw(rmtree);
use IO::Handle ;
use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# ----------------------
# declare subs
sub init_stuff ;
sub parse_arg ;
sub check_snana_version ;
sub check_multipleJobs ;
sub legacy_APPEND_FITRES ;
sub check_FITOPT_LINK ;
sub check_waitTime ;
sub parse_nmlFile ;
sub parse_PRIVATE_DATA_PATH ;
sub parse_VERSION ;
sub parse_MLAPPLY ;
sub getVersionList(@) ;
sub namelistLogical ;
sub set_OUTDIR ;
sub set_OUTFLAG ;
sub set_DOTEXT ;
sub get_NSPLIT ;
sub check_BATCH_NCORE ;
sub check_TIME_STAMP ;
sub get_DATADIR(@) ;
sub get_FORMAT(@) ;
sub get_PREFIX_SPLITJOB ;
sub get_PREFIX_BATCHFILE ;
sub get_suffixMerged ;
sub DONE_CHECK_LCFIT ;
sub DONE_CHECK_postProc ;
sub MERGE_DRIVER ;
sub merge_INIT ;

sub merge_outFiles ;
sub merge_outFiles_text ;  # driver to call next three
sub merge_text_TABLE ;
sub merge_text_FITRES ;    # XXX soon be obsolete soon (Mar 25 2019)  
sub merge_text_LCPLOT ;  
sub append_text_TABLE ;
sub append_text_FITRES ;  # XXX soon to be obsolete (Mar 25 2019)
sub mergeLog_update ;
sub mergeLog_cat ;

sub merge_end ;
sub merge_submit ;
sub get_mergeStatus ;

sub get_NSN_MERGED ;
sub copy_SIMGEN_DUMP ;
sub combine_all_text ;
sub combine_subset_text ;
sub combine_all_outFiles ;
sub error_check_LCFIT ;
sub clean_SPLIT_VERSION ; 
sub clean_text_LCPLOT ;  # was clean_DMP_SNFLUX
sub gzip_VERSION ;
sub gzip_SPLIT_JOBS ;
sub write_TRAILMARK ;
sub run_afterBurner ;
sub psnid_FOMsummary ;
sub hrecover ;
sub CLEANFILES_DRIVER ;
sub clean_gzipFITOPT ;

sub make_OUTDIR(@) ;
sub make_OUTDIR_symLink(@) ;
sub open_README ;
sub split_version(@) ;
sub init_SPLITSCRIPT(@) ;
sub split_job_prep ;
sub submitJobs_lcfit ;
sub submitJobs_postProc(@) ; # e.g., ML or SALT2mu

sub fitresNaN(@);
sub createJobs_LCFIT(@) ;
sub createJobs_SALT2mu(@) ;
sub createJobs_MLAPPLY(@) ;
sub create_batchScript(@) ;
sub updmask_CMD_LCFIT_ARRAY ;
sub get_simArg1_SALT2mu ;
sub get_simArg2_SALT2mu ;
sub batch_delayAdd ;
sub LOOPALL(@);

# -----------------
# declare globals

# hard-wired globals

my $SNDATA_ROOT       = $ENV{'SNDATA_ROOT'};
my $SNANA_DIR         = $ENV{'SNANA_DIR'};
my $DEFAULT_DATA_PATH = "$SNDATA_ROOT/lcmerge" ;
my $PATH_SNDATA_SIMDIR_DEFAULT = "$SNDATA_ROOT/SIM" ; # default simDir
my $SIMDIR_LISTFILE   = "$PATH_SNDATA_SIMDIR_DEFAULT/PATH_SNDATA_SIM.LIST";
my (@PATH_SIMDATA_LIST, @SIMDIR_LIST_ALL);

my $SCRIPTNAME       = "split_and_fit.pl" ;
my $SCRIPTNAME_FULL  = "$0";
my $JOBNAME_LCFIT    = "snlc_fit.exe" ;
my $JOBNAME_LCFIT_FULLPATH = `which $JOBNAME_LCFIT` ;

my $PSNID_FLAG       = 0 ;
my $JOBNAME_SPLIT    = "split_version.pl" ;  # used for text format only
my $JOBNAME_SALT2mu  = "SALT2mu.exe";
my $JOBNAME_COMBINE  = "combine_fitres.exe" ;
my $JOBNAME_MLAPPLY  = "UNKNOWN" ;  # read from ML-input file
my $JOBNAME_MKPLOTS  = "mkfitplots.pl"   ;
my $JOBNAME_SNTABLE  = "sntable_dump.pl" ;
my $MAXMERGE         = 120 ;  # abort if trying to merge more than this

my $SPLIT_JOBSUBDIR_LCFIT   = "SPLIT_JOBS_LCFIT";
my $SPLIT_JOBSUBDIR_SALT2mu = "SPLIT_JOBS_SALT2mu";
my $SPLIT_JOBSUBDIR_MLTRAIN = "SPLIT_JOBS_MLTRAIN" ; 
my $SPLIT_JOBSUBDIR_MLAPPLY = "SPLIT_JOBS_MLAPPLY" ; 
my $NJOBS_FILENAME          = "SUBMIT_NJOBS.DAT" ;
my $TIME_STAMP_FILENAME     = "SUBMIT_TIME_STAMP.DAT" ;
 
my $BATCH_TEMPLATE_KICP = '$SBATCH_TEMPLATES/SBATCH_kicp.TEMPLATE' ;

my $OUTINDX_HBOOK     = 0 ; # CERN's HBOOK
my $OUTINDX_ROOT      = 1 ; # CERN's root
my $OUTINDX_TEXT      = 2 ; # text format; must be last for append option
my $NMAX_OUTINDX      = 3 ;

# define quantities that depend on output table format
my @OUTFILE_FORMAT    = ( "HBOOK",  "ROOT" , "TEXT"   ) ;
my @SUFFIX_OUTFILE    = ( "HBOOK" , "ROOT",  "TEXT"   ) ;
my @NMLKEY_OUTFILE    = ( "HFILE_OUT",  "ROOTFILE_OUT", "TEXTFILE_PREFIX" );
my @JOBNAME_MERGE     = ( "merge_hbook.exe",  "merge_root.exe", "internal" ) ;
my @COMBINE_FMTARG    = ( "H",  "R",  "?" );  # arg for combine_fitres

my $ITABLE_SNANA  = 0 ;
my $ITABLE_FITRES = 1 ;
my $ITABLE_LCPLOT = 2 ;
my $MAXTABLE      = 3 ;
my @TABLELIST     = ( "SNANA" , "FITRES" , "LCPLOT" ) ;
my (@DOTEXT_TABLE);

my $SUFFIX_FITRES_TEXT  = "FITRES.TEXT" ;
my $SUFFIX_SNANA_TEXT   = "SNANA.TEXT" ;
my $SUFFIX_LCPLOT_TEXT  = "LCPLOT.TEXT" ; # match def in sntools_output_text.c
my $SUFFIX_LCLIST_TEXT  = "LCLIST.TEXT" ; # idem

my $SUFFIX_FITRES       = "FITRES" ;  # for catenated text
my $SUFFIX_SNANA        = "SNANA"  ;  # for catenated text
my $SUFFIX_LCPLOT       = "LCPLOT" ;  # idem

my $COMBINE_ALL_PREFIX = "COMBINED" ;
my $SUFFIX_LOG         = "LOG" ;
my $SUFFIX_DONE        = "DONE" ;
my $SUFFIX_BUSY        = "BUSY" ;
my $PREFIX_RUNSALT2mu  = "RUN_SALT2mu" ;
my $PREFIX_RUNML       = "RUN_ML" ;

my @COPY_INPUTFILE_KEYS = 
    ( "HFILE_KCOR", "KCOR_FILE", "SIMEFF_FILE", 
      "USERTAGS_FILE", "SNCID_LIST_FILE", "FUDGE_HOSTNOISE_FILE",
      "VPEC_FILE", "HEADER_OVERRIDE_FILE", "MAGCOR_FILE",
      "NONLINEARITY_FILE" );


my $CURRENT_DIR = `pwd` ;  $CURRENT_DIR  =~ s/\s+$// ; 

my $BATCH_MEM_LCFIT    = "1000";     # 1GB
my $BATCH_MEM_SALT2mu  = "2000";     # 2GB
my $BATCH_MEM_MLAPPLY  = "4000";     # 4GB
my $BATCH_MEM = 0 ;                  # user override (Feb 7 2018)

my $MXLC_PLOT_DEFAULT  = 3 ;         # default max plots per SPLIT job

my $qq = '"' ;
my $q  = "\'" ;

my $OPT_ABORT = 0 ; # parsing option
my $OPT_WARN  = 1 ; # 1=> leave warning message; 2=>quiet
my $OPT_QUIET = 2 ; # no warnong if key not found


# ---

my (@MSGERR, $SHELL_NAME, $NSPLIT, $INODE_SALT2mu, $INODE_ML );
my (@OUTDIR_LIST, @SCRIPT_LIST, @script_LIST, $SCRIPT_MASTER_POSTPROC );
my (@NODEMAP, @INODEMAP, @DATA_FORMAT) ;

my ($NFITOPT, $NFITOPT_LINK, $NFITJOB_TOT, $NFITJOB_SUBMIT, $NFITJOB_EXPECT ); 
my ($KILLJOBS_FLAG, $CLEANMASK, $PSNID_SUMLOG );
my ($MERGE_FLAG, $MERGE_LOGDIR, $MERGELOG, $MERGELOG2, $noMERGE_FLAG);
my ($HRECOVER_FLAG, $HRECOVER_PREFIX, $HRECOVER_FITOPT);
my ($PROMPT_FLAG, $LAUNCH_FLAG, $TRAILMARK_FLAG, $DONE_STAMP_FLAG );
my ($NVERSION, $NVERSION_SIM, $NVERSION_DATA,  @DATADIR_LIST);
my (@DONELIST_LCFIT, $NDONE_LCFIT, $NDONE_LINK, @NSN_MERGED, @CMD_LCFIT_ARRAY ) ;
my ($NGRACE, $NGRACE_TOTAL, $MSGWARN, $NABORT, $NABORT_TOTAL );
my ($FIT_README_FILE, @COPY_INPUTFILE_NAMES);
my ($PREFIX_MERGED, $PREFIX_SPLITJOB, $NFITOPT_MERGED);
my (@GZIP_STATUS, @MERGE_STATUS, $DONE_STATUS );
my ($SNTABLE_LIST, $PRIVATE_DATA_PATH );

my ($ALLDONE );
my $IFLAG_WAIT    = 0 ;
my $IFLAG_DOMERGE = 1 ;
my $IFLAG_DONE    = 2 ;
my $WAIT_CHECK_DONE  = 20 ; # wait time (seconds) between checks
my $TOTAL_WAIT_ABORT = 48 ; # abort after 24 hours
my $TOTAL_WAIT       = 0  ; # initialze total wait time

my (@SIM_FLAG_LIST, @SIM_SYMLINKDIR_LIST  ) ;
my $SIM_FLAG_NORMAL   = 1 ; # normal sim
my $SIM_FLAG_SYMLINK  = 2 ; # sim symbolic links to another sim
my $SIM_FLAG_SYMLINK2 = 3 ; # idem, but with OUTDIR in readme

my @OUTFLAG ;

my $GZIP_FLAG_TAR = 4;

# user-inputs
my (@VERSION_LIST, @VERSION_FITOPT_LIST, @PRIVATE_FLAG_LIST );
my (@FITOPT_LIST_ORIG, @FITOPT_LIST,  @FITOPT_LINK_LIST, @IFITOPT_LINK_LIST );
my ($SALT2mu_INFILE, $SALT2mu_SIMVERSION, $OPT_SPLIT, $MXLC_PLOT );
my ($SALT2mu_BIASCOR_PATH, $SALT2mu_CCPRIOR_PATH );
my (@SSH_NODELIST, $SSH_NNODE, $MIN_SNANA_VERSION );
my ($SNANA_LOGIN_SETUP, $SNANA_MODELPATH, $ENVDEF_MODELPATH);
my ($SPLIT_JOBDIR_LCFIT, $SPLIT_JOBDIR_ML, $SPLIT_JOBDIR_SALT2mu, $DONE_STAMP);
my ($GZIP_FLAG, $H2ROOT_FLAG, $ERRCHECK_FLAG );
my ($NFITERR_TOT, $NFITERR_NaN, $DO_FITOPT000 );
my ($VERSION_FITRES_COMBINE, $FITRES_COMBINE_FILE, $COMBINE_ALL_FLAG);
my (@COMBINE_SUBSET_FITOPT, $NJOB_PER_BATCHFILE ); 
my ($APPEND_FITRES_LIST, $APPEND_TABLE_LIST, $VERSION_AFTERBURNER );
my ($OUTDIR, $nmlFile_ORIG, $nmlFile, @CONTENTS_NMLFILE );

# my (@MERGED_TEXTFILE, @MERGED_textFile); # ascii fitres file
my (@MERGED_OUTFILE, @MERGED_outFile, @MERGED_DONEFILE );  # hbook or root
my (@NMISS_OUTFILE );
my ($DEBUG_DIR, $DEBUG_SDIR ); 
my ($BATCH_COMMAND, $BATCH_TEMPLATE, $BATCH_NCORE);
my ($NBATCH_TOTAL, $NBATCH_CREATE_JOB, $NBATCH_CREATE_FILE );
my ($OPT_LCPLOT, $ARG_LCPLOT, @DELAY_SUBMIT );
my ($NEXPECT_PER_MERGE, $NEXPECT_TOTAL );
my ($CURRENT_SNANA_VERSION );
my ($TIME_START, $TIME_END, $CDATE_SUBMIT );


# ML stuff, Feb 2017
my ($MLTRAIN_INFILE, @MLAPPLY_INFILE, $MLFLAG );
my $MLFLAG_TRAIN = 1 ;
my $MLFLAG_APPLY = 2 ;

# ================= BEGIN MAIN ===============

$TIME_START = time ;  # in seconds

my ($iver, $ifitopt, $ifitopt_start, $isplit ) ;

&init_stuff();

&parse_arg();

if ( $KILLJOBS_FLAG == 0 && $MERGE_FLAG == 0 && $PROMPT_FLAG == 1 ) 
{ sntools::check_multipleJobs($SCRIPTNAME) ; } # check for existing job

# check for special hmerge recovery option
if ( $HRECOVER_FLAG && $KILLJOBS_FLAG == 0 )  { &hrecover(); }

if ( $CLEANMASK !=0 ) { &CLEANFILES_DRIVER();  exit(0); }

&parse_nmlFile();

# check to kill jobs after reading list of nodes.
if ( $KILLJOBS_FLAG ) {  sntools::Killjobs(@SSH_NODELIST);  die "\n"; }

for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
    &get_DATADIR($iver);
    &get_FORMAT($iver);  # TEXT or FITS ?
}


&get_NSPLIT();  # move after get_DATADIR (Oct 3, 2011)

&check_BATCH_NCORE();  # sanity checks

# ----------------------------------------

if ( $MERGE_FLAG ) {
    &MERGE_DRIVER();
    print "\n DONE: all jobs have finished and are merged. \n" ;
    exit(0);
} # end of MERGE_FLAG block


# ----------------------------------------


print " -------------------------------------- \n" ;

if ( -d $OUTDIR ) { 
    print "\n Removing old OUTDIR: $OUTDIR ... \n";
    print " Please be patient. \n";
    $| = 1;  # flush stdout
    qx(rm -r $OUTDIR) ; 
}
&make_OUTDIR(0,-1);

&open_README();  # create README file associating ifitopt <-> FITOPT


for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
    &make_OUTDIR(1,$iver) ;
}

print "\n" ;

for  ( $isplit = 1; $isplit <= $NSPLIT; $isplit++ ) {
    &init_SPLITSCRIPT($isplit);
}

$NFITJOB_EXPECT = $NVERSION * $NFITOPT * $NSPLIT ;

# ---------------------------
# append RUN-script for each job

for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
    print "\t Append LCFIT-scripts for $VERSION_LIST[$iver] \n" ;
    for ( $ifitopt=$ifitopt_start ; $ifitopt < $NFITOPT; $ifitopt++ ) {
	for  ( $isplit = 1; $isplit <= $NSPLIT; $isplit++ ) {
	    &createJobs_LCFIT($iver,$ifitopt,$isplit) ;
	} 
    }
}


# ------------------------------
if ( $MLFLAG == $MLFLAG_APPLY ) {
    $INODE_ML = -1 ;
    my $NJOB = 0 ;
    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
	print "\t Create ML-APPLY scripts for $VERSION_LIST[$iver] \n" ;
	for ( $ifitopt=$ifitopt_start ; $ifitopt < $NFITOPT; $ifitopt++ ) {
	    &createJobs_MLAPPLY($iver,$ifitopt) ; # apply ML
	    $NJOB++ ;
	}
    }
    qx(chmod +x $SPLIT_JOBDIR_ML/$PREFIX_RUNML* );
    qx(echo "$NJOB"  > $SPLIT_JOBDIR_ML/$NJOBS_FILENAME );
}

# --------------------------------
if ( $SALT2mu_INFILE ne "" ) {
    $INODE_SALT2mu = -1  ;
    my $NJOB = 0 ;
    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
	print "\t Append SALT2mu-scripts for $VERSION_LIST[$iver] \n" ;
	for ( $ifitopt=$ifitopt_start ; $ifitopt < $NFITOPT; $ifitopt++ ) {
	    &createJobs_SALT2mu($iver,$ifitopt) ; # SALT2mu fit jobs
	    $NJOB++ ;
	}
    }
    qx(chmod +x $SPLIT_JOBDIR_SALT2mu/$PREFIX_RUNSALT2mu* );
    qx(echo "$NJOB"  > $SPLIT_JOBDIR_SALT2mu/$NJOBS_FILENAME );
}


# ------------------------------------
if ( $LAUNCH_FLAG == 0 ) 
{ die "\n\t !!!! ABORT BEFORE LAUNCH !!!! \n"; }

&split_job_prep();

&submitJobs_lcfit();

print "\n See output in $OUTDIR \n";

&merge_submit();

exit(0);

# ===================================
#
#            END OF MAIN
#
# ===================================


# ==========================================
sub init_stuff {

    $MIN_SNANA_VERSION = "v0";

    $SHELL_NAME    = sntools::shellName();

    $NFITJOB_TOT    = 0;
    $NFITJOB_SUBMIT = 0; # = NFITJOB_TOT - (those using LINK option)
    $NFITJOB_EXPECT = 0;
    $NJOB_PER_BATCHFILE = 1 ;

    $MERGE_FLAG    = 0 ; # set when MERGE flag is passed
    $noMERGE_FLAG  = 0;  # flag to NOT launch MERGE script in background
    $PROMPT_FLAG   = 1;  # user prompt; e.g., for multiple jobs
    $LAUNCH_FLAG   = 1;  # flag to actually launch jobs
    $TRAILMARK_FLAG = 0; # default is no TRAILMARK in sim-readme

    $KILLJOBS_FLAG = 0 ;
    $CLEANMASK     = 0 ;
    $DELAY_SUBMIT[0]  = 0. ;
    $DELAY_SUBMIT[1]  = 0. ;

    $HRECOVER_FLAG   = 0 ;
    $HRECOVER_PREFIX =  "" ;
    $HRECOVER_FITOPT = -9 ;

    $NVERSION      = 0 ;
    $NVERSION_SIM  = 0 ;
    $NVERSION_DATA = 0 ;

    $COMBINE_ALL_FLAG = 0;
    @COMBINE_SUBSET_FITOPT = ();

    $GZIP_FLAG     = 4 ; # Jun 6 2016 --> default is gzip SPLIT_JOBS 
    $H2ROOT_FLAG   = 0 ;
    $ERRCHECK_FLAG = 1 ;
    $NFITERR_TOT   = 0 ;
    $NFITERR_NaN   = 0 ;

    $BATCH_TEMPLATE = "" ;
    $NBATCH_CREATE_JOB  = 0 ;
    $NBATCH_CREATE_FILE = 0 ;
    
    $BATCH_NCORE = 0 ;
    $SSH_NNODE   = 0 ;

    $SALT2mu_INFILE = "" ;
    $SALT2mu_SIMVERSION   = "" ; # for Ia bias cor and for CCprior
    $SALT2mu_BIASCOR_PATH = "" ; # alternative for same
    $SALT2mu_CCPRIOR_PATH = "" ; # alternative for same
    
    $MLFLAG = 0;
    $JOBNAME_MLAPPLY = "" ;
    $MLTRAIN_INFILE  = "" ;
    @MLAPPLY_INFILE  = ();
    
    $OPT_SPLIT = 1;  # default is to split version for text file

    $OPT_LCPLOT  = 0 ;
    $ARG_LCPLOT  = '' ;

    $NABORT_TOTAL = 0 ;
    $NGRACE_TOTAL = 0 ;
    $DONE_STATUS  = "SUCCESS" ;

    @VERSION_LIST = () ;
    @VERSION_FITOPT_LIST = () ; # fitopt per version
    
    @FITOPT_LINK_LIST = ();
    $NDONE_LINK       = 0 ;
    
    $DO_FITOPT000  = 1; 
    $ifitopt_start = 0 ;

    my $indx; 
    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) 
	{  $NMISS_OUTFILE[$indx] = 0 ; }

    $VERSION_AFTERBURNER = '' ;

    $CDATE_SUBMIT = "UNKNOWN";

    @SIMDIR_LIST_ALL = ();
    if ( -e $SIMDIR_LISTFILE ) 
    { @SIMDIR_LIST_ALL = `cat $SIMDIR_LISTFILE` ; }

}  # end of init_stuff


# ========================================================
sub parse_arg() {

    # parse command-line arguments.

    my ($NARG, $arg, $argNext, $i ) ;

    $NARG = scalar(@ARGV);
    if ( $NARG < 1 ) {
	$MSGERR[0] = "Must give namelist filename as argument";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $nmlFile_ORIG = $ARGV[0] ;

    if ( $ARGV[0] eq "CLEAN"       )  { $CLEANMASK    = -1 ; }
    if ( $ARGV[0] eq "CLEANMASK"   )  { $CLEANMASK    = $ARGV[1] ; }

    for ( $i = 1; $i < $NARG ; $i++ ) {
        $arg     = $ARGV[$i];
	$argNext = $ARGV[$i+1];

        if ( $arg eq "MERGE"  ) 
	{ $MERGE_FLAG = 1 ; }

        if ( $arg eq "noMERGE"  ) 
	{ $noMERGE_FLAG = 1 ; }

        if ( $arg eq "KILL"   ) 
	{ $KILLJOBS_FLAG = 1 ; }

        if ( $arg eq "1"   ) 
	{ $DELAY_SUBMIT[0] = 1 ; }

        if ( $arg eq "HRECOVER"   )  { 
	    $HRECOVER_FLAG   = 1 ;
	    $HRECOVER_PREFIX = $ARGV[$i+1] ; 
	    $HRECOVER_FITOPT = $ARGV[$i+2] ; 
	}

        if ( $arg eq "NOPROMPT"  )  { $PROMPT_FLAG    = 0 ; }
        if ( $arg eq "NOLAUNCH"  )  { $LAUNCH_FLAG    = 0 ; }
        if ( $arg eq "NOSUBMIT"  )  { $LAUNCH_FLAG    = 0 ; }
        if ( $arg eq "TRAILMARK" )  { $TRAILMARK_FLAG = 1 ; }

	if ( $arg eq "-DATE_SUBMIT" ) { $CDATE_SUBMIT = $argNext ; }

	if ( $arg eq "KICP" || $arg eq "kicp" ) 
	{ $BATCH_TEMPLATE = "$BATCH_TEMPLATE_KICP" ; }
    }

}  # end of parse_arg


# =====================================
sub parse_nmlFile {

    # parse keys outside namelist, and also
    # needed variables inside namelist
    #
    # Dec 4, 2012: strip off optional quotes (' or ") around each VERSION name
    # Apr 15, 2013: check FITOPT: keys for files to copy.
    # Dec 20, 2017: call getVersionList($ver) for VERSION+FITOPT key
    #
    my ($key, @tmp,  @words, $vtmp, $fitopt, $tmp_fitopt );
    my ($NTMP, $FOUND_VERSION_KEY, @tmpVerList, @tmpPathList ) ;
    my ($tmpLine, $line, @TMPVER, $ver, @TMPNODES, $inFile, $idx );

    # ------- BEGIN --------
    
    if ( !(-e $nmlFile_ORIG ) ) {
	@MSGERR[0] = "Cannot find input file: $nmlFile_ORIG ";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $nmlFile = $nmlFile_ORIG ;

    @CONTENTS_NMLFILE = ();
    sntools::loadArray_fromFile($nmlFile_ORIG, \@CONTENTS_NMLFILE);

    # allow modifiying nmlFile here (e..g, removing LTUP keys).
    # NOTHING to change.

    # if nml file has changed (from removing legacy variables),
    # scoop up contents again. Note that comment lines are removed.
    if ( $nmlFile ne $nmlFile_ORIG ) {
	@CONTENTS_NMLFILE = ();
	sntools::loadArray_fromFile($nmlFile, \@CONTENTS_NMLFILE);
    }

    $DONE_STAMP      = ' ' ;
    $DONE_STAMP_FLAG = 0;
    $key   = "DONE_STAMP:" ;
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    $NTMP  = scalar(@tmp) ;
    if ( $NTMP > 0 ) {
	$DONE_STAMP      = "$tmp[0]" ;
	$DONE_STAMP_FLAG = 1 ;
	if ( -e $DONE_STAMP )  { qx(rm $DONE_STAMP) ; }
    }


    # check for min SNANA version (Aug 24 2015)
    &check_snana_version();

    # extract default version from &SNLCINP  namelist
    $key             = "VERSION_PHOTOMETRY" ;
    $VERSION_LIST[0] = sntools::get_namelistValue($nmlFile,$key);
    $NVERSION        = 1;

    # extract list of tables to create (need parse function)    
    $key             = "SNTABLE_LIST" ;
    $SNTABLE_LIST    = sntools::get_namelistValue($nmlFile,$key);
    &set_DOTEXT();

    # ------------------------------------
    # extract misc input files from nmlFile that need to be copied.
    # Check nominal namelist and also check [VERSION+]FITOPT: key

    @COPY_INPUTFILE_NAMES = () ;
    foreach $key ( @COPY_INPUTFILE_KEYS ) {

	print "\t check file key $key \n" ;
	$inFile = sntools::get_namelistValue($nmlFile,$key);
	@COPY_INPUTFILE_NAMES = ( @COPY_INPUTFILE_NAMES , $inFile ) ;

	$idx = 0 ;
	@tmp = sntools::parse_array($key,99, $OPT_QUIET, @CONTENTS_NMLFILE);
	$NTMP = scalar(@tmp) ;
	foreach $tmpLine (@tmp) {
	    @words  = split(/\s+/,$tmpLine) ;
	    $inFile = $words[$1];
	    @COPY_INPUTFILE_NAMES = ( @COPY_INPUTFILE_NAMES , $inFile ) ;
	}              
    }

    # ------------------------------------

    # jobname (default is snlc_fit.exe)
    $key     = "JOBNAME_LCFIT:" ;
    @tmp     = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) {
	$JOBNAME_LCFIT  = $tmp[0] ;
	$JOBNAME_LCFIT  =~ s/\s+$// ;   # trim trailing whitespace

	# check if we have a psnid job
	if ( index($JOBNAME_LCFIT,"psnid") >= 0 ) { $PSNID_FLAG = 1; }

	# check if job is in current dir (Aug 23 2015)
	# If there is a slash, then use full path.
	if ( -e $JOBNAME_LCFIT  && index($JOBNAME_LCFIT,"/")<0 ) {
	    $JOBNAME_LCFIT = "$CURRENT_DIR/$JOBNAME_LCFIT" ;
	}
    }
    $JOBNAME_LCFIT_FULLPATH = `which $JOBNAME_LCFIT` ;
    $JOBNAME_LCFIT_FULLPATH =~ s/\s+$// ;   # trim trailing whitespace

    # -----------------   check ssh nodes --------------
    @SSH_NODELIST = () ;
    $key = "NODELIST:" ;
    @TMPNODES  = sntools::parse_array($key,30,$OPT_QUIET, @CONTENTS_NMLFILE);
    foreach $tmpLine ( @TMPNODES ) {
	@tmp = split(/\s+/,$tmpLine) ;
	@SSH_NODELIST =  ( @SSH_NODELIST , @tmp );
	$SSH_NNODE    = scalar(@SSH_NODELIST);
    }

    # --------- check for batch info --------------
    $key     = "BATCH_INFO:" ;
    @tmp     = sntools::parse_array($key,3,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) {
	@words            = split(/\s+/,$tmp[0]) ;
	$BATCH_COMMAND  = $words[0] ;
	if ( $BATCH_TEMPLATE eq '' ) { $BATCH_TEMPLATE = $words[1] ; }
	$BATCH_NCORE    = $words[2] ;

	# allow for ENV (May 22 2017)
	$BATCH_TEMPLATE = qx(echo $BATCH_TEMPLATE);
	$BATCH_TEMPLATE  =~ s/\s+$// ;   # trim trailing whitespace
    }

    if ( $SSH_NNODE > 0 && $BATCH_NCORE > 0 ) {
	@MSGERR[0] = "Cannot specify both NODELIST and BATCH_INFO keys.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    if ( $SSH_NNODE == 0 && $BATCH_NCORE == 0 ) {
	@MSGERR[0] = "Must specify BATCH_INFO or SSH nodes with:";
	@MSGERR[1] = "NODELIST: <node1> <node2> etc ... " ;
	@MSGERR[2] = "    or " ;
	@MSGERR[3] = "BATCH_INFO: <command> <templateFile> <Ncore>" ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    $key     = "BATCH_MEM:" ;
    @tmp     = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) {
	$BATCH_MEM = $tmp[0] ;
    }

    # - - -  -
    $key = "NJOB_PER_BATCHFILE:" ;
    @tmp     = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) {
	$NJOB_PER_BATCHFILE = $tmp[0] ;
	$WAIT_CHECK_DONE    = 120 ;   # bigger wait here
    }

    # ---------- check for snana 'setup' command --------------
    $SNANA_LOGIN_SETUP = "" ;
    $key = "SNANA_LOGIN_SETUP:" ;
    @tmp = sntools::parse_array($key, 20,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { $SNANA_LOGIN_SETUP = "$tmp[0]" ; }

    $SNANA_MODELPATH  = "" ;
    $ENVDEF_MODELPATH = "" ;
    $key = "SNANA_MODELPATH:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    $SNANA_MODELPATH = "$tmp[0]" ;  # @tmp -> $tmp[0]
    if ( scalar(@tmp) > 0 ) {
	$ENVDEF_MODELPATH = 
	    sntools::setEnvString($SHELL_NAME, 
				  "SNANA_MODELPATH", $SNANA_MODELPATH);
    }

    # --------- check for output directory ----------------
    $key    = "OUTDIR:" ;
    @tmp    = sntools::parse_array($key,1,$OPT_ABORT, @CONTENTS_NMLFILE);
    &set_OUTDIR($tmp[0]);

    $SPLIT_JOBDIR_LCFIT    = "$OUTDIR/$SPLIT_JOBSUBDIR_LCFIT" ;
    $MERGE_LOGDIR          = "$SPLIT_JOBDIR_LCFIT/MERGELOGS" ;
    $SPLIT_JOBDIR_SALT2mu  = "$OUTDIR/$SPLIT_JOBSUBDIR_SALT2mu" ; 

    # ------- check for private data path -------------
    &parse_PRIVATE_DATA_PATH();

    # ---------------- check photometry versions ---------------
    # if they exist, then overwrite the VERSION parsed above.

    $FOUND_VERSION_KEY = &parse_VERSION();

    $key    = "DO_FITOPT000:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { $DO_FITOPT000 = $tmp[0] ;   }

    # Now check for fit-options; tack on the default "NONE".
    $key    = "FITOPT:" ;
    @FITOPT_LIST_ORIG = 
	sntools::parse_array($key,99,$OPT_QUIET, @CONTENTS_NMLFILE);
    for(@FITOPT_LIST_ORIG) { s/"/\'/g ; } # double quote -> single quote 

    @FITOPT_LIST_ORIG = ( '[DEFAULT] NONE' , @FITOPT_LIST_ORIG ) ;
    @FITOPT_LIST = ();

    foreach $fitopt (@FITOPT_LIST_ORIG) {
	$tmp_fitopt  = "$fitopt" ;
	$tmp_fitopt  =~ s/\[[^)]*\]//g ;   # remove optional tagname in []
	$tmp_fitopt  =~ s/^\s+// ;         # remove leading spaces
	@FITOPT_LIST = ( @FITOPT_LIST, "$tmp_fitopt" ) ;
    }
    &check_FITOPT_LINK();

    $NFITOPT     = scalar(@FITOPT_LIST) ;    
   
    # Mar 2013: check for VERSION+FITOPT where a separate 
    #           FITOPT is assigned to each version
    $key    = "VERSION+FITOPT:" ; 
    @TMPVER = sntools::parse_array($key,20,$OPT_QUIET, @CONTENTS_NMLFILE);
    $NTMP   = scalar(@TMPVER) ;
    if ( $NTMP > 0 ) {	

	if (  $FOUND_VERSION_KEY ) {
	    $MSGERR[0] = "Cannot mix VERSION and VERSION+FITOPT keys." ;
	    $MSGERR[1] = "Check $nmlFile" ;
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}

	$NVERSION   = $NTMP ;
	@VERSION_LIST = ();
	# strip off optional quotes around version name; check ' and "
	foreach $line (@TMPVER) {
	    @words     = split(/\s+/,$line) ;
	    $ver       = $words[0];
	    $ver =~ s/$qq//g ;
   	    $ver =~ s/$q//g ;
	    
	    $fitopt = "@words[1 .. $#words]" ;
	    &getVersionList($ver, \@tmpVerList, \@tmpPathList); 
	    
	    @VERSION_LIST        = ( @VERSION_LIST, @tmpVerList );
	    @VERSION_FITOPT_LIST = ( @VERSION_FITOPT_LIST, $fitopt );
	}
    }
    
    # -----------------------------
    $APPEND_TABLE_LIST = "" ;
    $key   = "APPEND_FITRES:" ;  # legacy key
    @tmp   = sntools::parse_array($key, 20, $OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 )  { &legacy_APPEND_FITRES("$tmp[0]");  }

    $key   = "APPEND_TABLE_TEXT:" ; # nominal key as of Mar 25 2019
    @tmp   = sntools::parse_array($key, 20, $OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 )  { $APPEND_TABLE_LIST = "$tmp[0]" ; }

    # - - - - - - 

    $VERSION_AFTERBURNER = "" ;
    $key   = "VERSION_AFTERBURNER:" ;
    @tmp   = sntools::parse_array($key, 20,  $OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 )  { $VERSION_AFTERBURNER = "$tmp[0]" ; }

    $VERSION_FITRES_COMBINE = "" ;
    $key   = "VERSION_FITRES_COMBINE:" ;
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    $VERSION_FITRES_COMBINE = "@tmp" ;

    $FITRES_COMBINE_FILE = "" ;
    $key  = "FITRES_COMBINE_FILE:" ;
    @tmp   = sntools::parse_array($key, 99, $OPT_QUIET, @CONTENTS_NMLFILE);
    $FITRES_COMBINE_FILE = "@tmp" ;

    $COMBINE_ALL_FLAG = 0 ;
    $key  = "COMBINE_ALL_FLAG:" ;
    @tmp  = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    $COMBINE_ALL_FLAG = "$tmp[0]" ;

    # -----------------------------------------
    # read new key, but not yet implemented
    $key  = "COMBINE_SUBSET_FITOPT:" ;
    @tmp  = sntools::parse_array($key, 2, $OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { @COMBINE_SUBSET_FITOPT   = @tmp ;  }
    # -----------------------------------------

    $key  = "OPT_SPLIT:" ;
    @tmp  = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { $OPT_SPLIT = $tmp[0] ; }

    $key    = "GZIP_FLAG:" ;
    @tmp    = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { $GZIP_FLAG = $tmp[0]; }

    $key    = "DELAY_SUBMIT:" ;
    @tmp    = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { 
	@DELAY_SUBMIT =  split(/,/,$tmp[0]); 
	if ( scalar(@DELAY_SUBMIT) == 1 ) { $DELAY_SUBMIT[1] = 0; }
    }

    $key    = "H2ROOT_FLAG:" ;
    @tmp    = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { $H2ROOT_FLAG = $tmp[0]; }

    $key  = "PLOTOPT:" ;
    @tmp  = sntools::parse_array($key, 20, $OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { 
	$OPT_LCPLOT  = 1 ;
	$ARG_LCPLOT  = $tmp[0]; 
	$ARG_LCPLOT  =~ s/\s+$// ;   # trim trailing whitespace
    }


    # -----
    # allow 'SALT2mu' or 'SALT2MU'

    $key = "SALT2MU_INFILE:" ;
    @tmp  = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { $SALT2mu_INFILE = $tmp[0] ; }

    $key = "SALT2mu_INFILE:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { $SALT2mu_INFILE = $tmp[0] ; }

    $key = "SALT2mu_SIMVERSION_INPUT:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { $SALT2mu_SIMVERSION = $tmp[0] ; }

    $key = "SALT2mu_BIASCOR_PATH:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 )  { 
	$SALT2mu_BIASCOR_PATH = $tmp[0] ; 

	if ( length($SALT2mu_CCPRIOR_PATH) == 0 ) 
	{ $SALT2mu_CCPRIOR_PATH = $tmp[0] ;  }
    }

    $key = "SALT2mu_CCPRIOR_PATH:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) {  $SALT2mu_CCPRIOR_PATH = $tmp[0] ;     }

    # check for ML options (Feb 2017)
    &parse_MLAPPLY();
    
    # --------------------------


    # check WAIT_CHECK_DONE
    $key   = "WAIT_CHECK_DONE:" ;   # seconds
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    $NTMP  = scalar(@tmp) ;
    if ( $NTMP > 0 ) { $WAIT_CHECK_DONE = "$tmp[0]" ; }

    $key = "TOTAL_WAIT_ABORT:" ;  # hours
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_NMLFILE);
    $NTMP  = scalar(@tmp) ;
    if ( $NTMP > 0 ) { $TOTAL_WAIT_ABORT = "$tmp[0]" ; }

    # --------------
    # parse namelist quantities needed here 

    $key = "MXLC_PLOT" ;
    $MXLC_PLOT = sntools::get_namelistValue($nmlFile,$key);

    # set OUTFLAG for FITRES, HBOOK, ROOT ...
    &set_OUTFLAG();

    # --------------------------------
    # print summary
    print "\n ------------------------------------------- \n" ;
    print " SNDATA_ROOT = $SNDATA_ROOT \n";

    if ( $SSH_NNODE > 0 ) {
	print " SSH_NNODE       = $SSH_NNODE \n" ;
	print " SSH_NODELIST    = @SSH_NODELIST \n";
    }
    if ( $BATCH_NCORE > 0 ) {
	my $BCMD = "$BATCH_COMMAND $BATCH_TEMPLATE";
	my $B3   = "$BATCH_NCORE cores";
	if ( $BATCH_MEM > 0 ) { $BATCH_MEM_LCFIT = $BATCH_MEM ; }
	print " BATCH_COMMAND: $BCMD ($B3) \n";
	print " NJOB_PER_BATCHFILE: $NJOB_PER_BATCHFILE \n";
	print " BATCH_MEM(LCFIT): $BATCH_MEM MB \n";
    }


    print " SNANA_LOGIN_SETUP = '${SNANA_LOGIN_SETUP}' \n" ;
    print " OUTDIR      = ${OUTDIR} \n";
    print " VERSION     = @VERSION_LIST \n";
    print " VERSION_FITRES_COMBINE = '${VERSION_FITRES_COMBINE}' \n" ;
    print " FITRES_COMBINE_FILE    = '${FITRES_COMBINE_FILE}' \n" ;
    print " APPEND_TABLE_TEXT      = '${APPEND_TABLE_LIST}' \n";
    print " VERSION_AFTERBURNER    = '${VERSION_AFTERBURNER}' \n";

    print " COMBINE_ALL_FLAG = $COMBINE_ALL_FLAG \n" ;
    foreach $tmpLine ( @COMBINE_SUBSET_FITOPT ) {
	print " COMBINE_SUBSET_FITOPT = $tmpLine \n";
    }

    print " SALT2mu_INFILE   = $SALT2mu_INFILE \n" ;
    print " PLOTOPT          = $OPT_LCPLOT ($ARG_LCPLOT) \n" ;
    print " MLFLAG           = $MLFLAG \n" ;
    my ($indx);
    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) { 
	print "  OUTFLAG($OUTFILE_FORMAT[$indx]) = $OUTFLAG[$indx]\n" ;
    }

    foreach $fitopt ( @FITOPT_LIST ) { print " FITOPT    = $fitopt \n"; }    

    print " DELAY(between each batch job) = $DELAY_SUBMIT[0] sec. \n";
    print " DELAY(between VERSION/FITOPT) = $DELAY_SUBMIT[1] sec. \n";
    print " WAIT_CHECK_DONE = $WAIT_CHECK_DONE  seconds. \n";

    if ( $DO_FITOPT000 == 0 ) {
	$ifitopt_start = 1 ;
	print "\n\t NOTE: FITOPT000 will be SKIPPED !!!! \n\n";
    }

    if ( $DONE_STAMP_FLAG ) {  print " DONE_STAMP  = $DONE_STAMP \n"; }

    print " \n";

    return ;
    
} # end of parse_nmlFile


sub parse_PRIVATE_DATA_PATH {

    my ($KEY,$PATH);

    $PRIVATE_DATA_PATH = "" ;

    $KEY    = "PRIVATE_DATA_PATH" ;
    $PATH   = sntools::get_namelistValue($nmlFile, $KEY);

    # echo the PATH in case it's using an ENV
    $PATH   = qx(echo $PATH);
    $PATH  =~ s/\s+$// ;   # trim trailing whitespace

    $PRIVATE_DATA_PATH = "$PATH" ;  # store in global

    if ( length($PATH) > 0 && $PATH ne 'pwd' ) {
	print "  Found PRIVATE_DATA_PATH = $PRIVATE_DATA_PATH \n";
	unless ( -d $PATH ) {
	    $MSGERR[0] = "Cannot find PRIVATE_DATA_PATH = " ;
	    $MSGERR[1] = "$PATH" ;
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}
    }

} # end 


sub parse_VERSION {

    # parse VERSION key(s) in nml File, and check for wild-cards (*)
    # for series of simulated versions.
    # Function returns 1 if "VERSION:" key is found; returns 0 if not.
    #
    # Feb 27 2017; abort if any version has a dot.
    # Jul 18 2017: fix to work with PRIVATE_DATA_PATH
    # Mar 07 2019: abort on missing version.
    # Jun 15 2019: 
    #   + init @tmpVerList & @tmpPathList before call to &getVersionList;
    #     fixes subtle bug with wild card in VERSION key.
    #

    my ($key, $ver, $FOUND_VERSION_KEY, @TMPVER, $NTMP);
    my (@tmpVerList, @tmpPathList );
    my $NOT_FOUND = 0 ;

    $key    = "VERSION:" ; 
    $FOUND_VERSION_KEY = 0 ;
    @TMPVER = sntools::parse_array($key, 1, $OPT_QUIET, @CONTENTS_NMLFILE);
    $NTMP   = scalar(@TMPVER) ;
    if ( $NTMP > 0 ) {	
	@VERSION_LIST = ();
	# strip off optional quotes around version name; check ' and "
	foreach $ver (@TMPVER) {
	    $ver =~ s/$qq//g ;
   	    $ver =~ s/$q//g ;

	    # abort if there is a dot (Feb 27 2017)
	    if ( index($ver,'.') > 0 ) {
		$MSGERR[0] = "Invalid VERSION: $ver  ";
		$MSGERR[1] = "--> no dots allowed for split_and_fit" ;
		sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;	
	    }

            @tmpVerList= ();  @tmpPathList = ();
	    &getVersionList($ver, \@tmpVerList,\@tmpPathList ); 
#	    print " xxx '@tmpVerList'  | @tmpPathList \n";

	    if ( scalar(@tmpVerList) > 0 ) {
		@VERSION_LIST      = ( @VERSION_LIST,      @tmpVerList );
		@PATH_SIMDATA_LIST = ( @PATH_SIMDATA_LIST, @tmpPathList );
		$FOUND_VERSION_KEY = 1 ;
	    }
	    else {
		$NOT_FOUND++ ;
		print " WARNING: could not find VERSION '$ver' \n";
	    }
	}
    }

    if ( $NOT_FOUND > 0 ) {
	$MSGERR[0] = "Could not find $NOT_FOUND VERSION keys.";
	$MSGERR[1] = "Se WARNING messages above.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    $NVERSION = scalar(@VERSION_LIST);
    
    if ( $NVERSION == 0 ) {
	$MSGERR[0] = "  NVERSION = 0 ?!?!?! ";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;	
    }

    my $LDMP = 0 ;
    if ( $LDMP ) {
	foreach $ver (@VERSION_LIST) 
	{ print " xxx found VERSION = '$ver' \n" ; }
    }
    
    return($FOUND_VERSION_KEY) ;
    
} # end parse_VERSION


sub getVersionList(@) {

    # Created Nov 20 2017
    # Without a wild card, just return input $version.
    # If there is a wild card, do 'ls' to get list of 
    # all versions
    #
    # Note that global $PATH_SNDATA_SIMDIR can be changed
    # if versions are found in user-defined PATH_SNDATA_SIM dir.
    #
    # Mar 8 2019: abort if version is found in 2 places.
    #
    my ($version, $verList, $pathList) = @_ ;

    my (@dirList,  @tmpList, $tmpDir );
    my ($i, $ndir, $nver, $lens, $nsimDir);
    my ($NFOUND, $simDir);

    # - - - - - - - - - - - - - - - -
#    @verList = ();  @pathList = ();  # init outputs
    $NFOUND=0;  	$simDir = "NotSim" ;

    # store list of all possible simDirs in local array
    @dirList = @SIMDIR_LIST_ALL ;

    # don't forget default output sim dir and output data dir
    @dirList = ( @dirList, $PATH_SNDATA_SIMDIR_DEFAULT ); 
    $nsimDir = scalar(@dirList); # all possible simDirs

    # tack on possible real-data subDirs
    @dirList = ( @dirList, "$DEFAULT_DATA_PATH" );
    $lens = length($PRIVATE_DATA_PATH);
    if ( $lens > 0 && "$PRIVATE_DATA_PATH" ne "$DEFAULT_DATA_PATH" ) 
    { @dirList = ( @dirList, $PRIVATE_DATA_PATH ); }

    $ndir = scalar(@dirList);

#    print " xxx ---------------- \n";
#    print " xxx look for VERSION '$version' \n";

    for($i=0; $i < $ndir; $i++ ) { 
	$tmpDir   = $dirList[$i] ;
	$tmpDir   = qx(echo $tmpDir);
	$tmpDir   =~ s/\s+$// ;   # trim trailing whitespace
	$dirList[$i]  = $tmpDir ;

	@tmpList = qx(cd $tmpDir; ls -d $version 2>/dev/null ) ;
	if ( scalar(@tmpList) > 0 ) { 
	    @$verList = @tmpList ;  
	    if ( $i < $nsimDir ) { $simDir = "$tmpDir" ; }
	    $NFOUND++ ;
	    print "  Path($version):  $tmpDir \n";
	}
    }  # end i loop over paths to check for data files
        
    # ==============================
    if ( $NFOUND != 1 ) {
	print "\n\n PRE-ABORT DUMP: \n";
	print " List of searched data/sim paths: \n";
	for $tmpDir (@dirList) { print "   $tmpDir \n"; } 
	print "  (PRIVATE_DATA_PATH = $PRIVATE_DATA_PATH) \n";
	$MSGERR[0] = " Found $NFOUND matches for VERSION $version";
	$MSGERR[1] = " Must be 1 and only 1 match.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;	
    }


    # remove trailing spaces, and set path
    $nver = scalar(@$verList);
    for($i=0; $i < $nver; $i++ ) { 
	@$verList[$i]  =~ s/\s+$// ;  
	@$pathList[$i] = "$simDir" ;
    }

    return ;
#    return(\@verList,\@pathList);

} # end getVersionList

sub parse_MLAPPLY {

    # Created Feb 8 2017
    # Parse MLAPPLY_INFILE key and figure out
    # ML-input file for each FITOPT

    my ($key, @tmp, $jstar, $tmpFile, $f0, $ifitopt, $fff );

    $key    = "MLAPPLY_INFILE:" ;
    @tmp = sntools::parse_array($key, 1, $OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) == 0 ) { return ; }

    $tmpFile   = $tmp[0] ;
    $MLFLAG    = $MLFLAG_APPLY ;
    $SPLIT_JOBDIR_ML = "$OUTDIR/$SPLIT_JOBSUBDIR_MLAPPLY" ; 

    # check for wildcar (*) to indicate FITOPT-dependence
    $jstar = index($tmpFile,'*');
    for($ifitopt =0; $ifitopt < $NFITOPT; $ifitopt++ ) {
	if ( $jstar < 0 ) {
	    $MLAPPLY_INFILE[$ifitopt] = "$tmpFile" ;
	}
	else {
	    $fff = sprintf("%3.3d", $ifitopt);
	    $MLAPPLY_INFILE[$ifitopt] = 
		substr($tmpFile,0,$jstar) . "$fff" . 
		substr($tmpFile,$jstar+1,99)
	}
    }  # end ifitopt 
    
    # read job name from first MLAPPLY_INFILE
    $key  = "MLAPPLY_CODENAME:" ;
    $f0   = $MLAPPLY_INFILE[0] ;
    @tmp  = sntools::parse_line($f0, 1, $key, $OPT_QUIET ) ;
    if ( scalar(@tmp) == 0 ) {
	$MSGERR[0] = " Missing 'MLAPPLY_CODENAME:' key in file " ;
	$MSGERR[1] = " specified by  MLAPPLY_INFILE: ";
	$MSGERR[2] = "   $tmpFile ";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;	
    }
    $JOBNAME_MLAPPLY = "$tmp[0]";


    my $LDMP=0 ;
    if ( $LDMP ) {
	print "\n ML-APPLY DUMP: \n";
	print "\t MLFLAG = $MLFLAG \n";
	print "\t SPLIT_JOBDIR_ML = $SPLIT_JOBDIR_ML \n";
	print "\t MLAPPLY_INFILE = '$tmpFile (jstar=$jstar) \n";
	
	for($ifitopt =0; $ifitopt < $NFITOPT; $ifitopt++ ) {
	    print "\t MLAPPLY file for FITOPT=$ifitopt -> " . 
		"$MLAPPLY_INFILE[$ifitopt]\n";
	}

	die "\n\n Done with ML-APPLY DUMP --> ABORT. \n";
    }

    return ;

} # end parse_MLAPPLY

# ====================================
sub check_FITOPT_LINK {

    # check for FITOPT of the form "->FITOPT###"
    # which means to just add a symbolic link without
    # running any fit jobs
    #
    # Load global array  FITOPT_LINK_LIST
    # and set  FITOPT_LIST[$ifitopt] = "NONE" for symbolic links.
    #
    my($FITOPT, $ifitopt, $ifitopt_link, $LINK, $txt);

    $NFITOPT_LINK = 0 ;
    $ifitopt = 0 ;
    foreach $FITOPT (@FITOPT_LIST) {
	$FITOPT_LINK_LIST[$ifitopt] = "" ;
	$IFITOPT_LINK_LIST[$ifitopt] = -1 ;
	if ( substr($FITOPT,0,8) eq "->FITOPT") {
	    $LINK = substr($FITOPT,2,9) ;
	    $ifitopt_link = substr($FITOPT,8,3);
	    $FITOPT_LINK_LIST[$ifitopt] = "$LINK" ;
	    $IFITOPT_LINK_LIST[$ifitopt] = int($ifitopt_link) ;
	    $FITOPT_LIST[$ifitopt]      = "NONE" ;
	    $txt = sprintf("FITTOP%3.3d linked to %s (%d)",
			   $ifitopt, $LINK, $ifitopt_link );
	    printf("  $txt \n");
	    $NFITOPT_LINK++ ;
	}
	$ifitopt++ ;
    }

    return ;

} # end check_FITOPT_LINK


# ====================================
sub namelistLogical {
    my ($file,$var) = @_ ;
    # returns 0 or 1 if nml variable $nmlVar has value of T or .TRUE.   
    my ($tmp);
    $tmp  = sntools::get_namelistValue($file,$var);
    if ( $tmp eq "T"  ||  $tmp eq ".TRUE."  ) { return 1; }

    return 0;

} # end of namelistLogical 

# =========================================
sub set_OUTDIR { 
    my ($outDir) = @_ ;

    # Created Dec 13 2014
    # Set global $OUTDIR based on $outDir that was read from
    # input file. Here we remove quotes (if any).
    # Also, if there are no slashes, then return 
    #    $OUTDIR = $CURRENT_DIR/$outDir 
    # so that $OUTDIR is a subdir below where we are now.

    $outDir =~ s/$qq//g ;  # remove optional quotes
    $outDir =~ s/$q//g ;   # idem
    $OUTDIR = $outDir ;

    # check for ENV
    $OUTDIR = qx(echo $OUTDIR);
    $OUTDIR  =~ s/\s+$// ;   # trim trailing whitespace

    if ( index($outDir,'/') < 0 ) {
	$OUTDIR = "$CURRENT_DIR/$outDir" ;
    }

    return ;

} # end of get_OUTDIR

# =========================================
sub legacy_APPEND_FITRES {

    my($LEGACY_LINE) = @_ ;

    # Mar 25 2019
    # if using legacy APPEND_FITRES key, remove TABLE value,
    # and transfer variable names to APPEND_TABLE_LIST.
    # This allows using the old legacy key.

    my (@wdlist, $VARLIST);
    @wdlist  = split(/\s+/,$LEGACY_LINE) ;
    $VARLIST  = "@wdlist[1 .. $#wdlist]" ;  # leave out table name
    $APPEND_TABLE_LIST = "$VARLIST";

    print "\n";
    print " **********************************************************\n";
    print " WARNING: Processing legacy APPEND_FITRES key ... \n";
    print "          please switch to APPEND_TABLE_TEXT key. \n";
    print " **********************************************************\n";
    print "\n";

    return;

} # end legacy_APPEND_FITRES


# ======================
sub set_OUTFLAG {

    # July 25 2013: parse nml file and set @OUTFLAG array.

    my ($outFile, $key, $indx, $NFORMAT, @bla, $SNANAJOB, $FLAG );

    $NFORMAT = 0 ;
    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) { 
	$OUTFLAG[$indx] = 0 ; 
	$key = "$NMLKEY_OUTFILE[$indx]" ;
	$outFile = sntools::get_namelistValue($nmlFile,$key);
	my $ISBLANK = ( $outFile eq '' || $outFile eq ' ' );
	if ($ISBLANK == 0 ) { $OUTFLAG[$indx] = 1; $NFORMAT++ ; }
    }

    # -------------------------------------------------------
    # FITRES format is treated a little different because
    # it is ALWAYS set except for snana.exe.
    
    @bla = grep { /snana.exe/ } $JOBNAME_LCFIT ;
    if ( scalar(@bla) > 0 ) { $SNANAJOB = 1; } else { $SNANAJOB = 0 ; }
    
    $indx = $OUTINDX_TEXT ;
    $FLAG = $OUTFLAG[$indx] ;
    if ( $FLAG == 0 && $SNANAJOB == 0 ) 
    { $OUTFLAG[$indx] = 1; $NFORMAT++ ; } # set this flag anyway

    if ( $FLAG == 1 && $SNANAJOB == 1 ) 
    { $OUTFLAG[$indx] = 0; $NFORMAT-- ; } # turn off for snana.exe job
    
    # -------------------------------------------------------

    if ( $NFORMAT == 0 ) {
	$MSGERR[0] = "Must specify output table format(s) in &SNLCINP." ;
	$MSGERR[1] = "Set one or more of the following to be non-blank:";
	$MSGERR[2] = "   @NMLKEY_OUTFILE ";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;
    }
	 
    if ( $OUTFLAG[$OUTINDX_ROOT] && $H2ROOT_FLAG ) {
	$MSGERR[0] = "Cannot set both H2ROOT_FLAG and ROOTFILE_OUT" ;
	$MSGERR[1] = "Select only 1 method to produce a root file.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

} # end of set_OUTFLAG

# ========================
sub get_NSPLIT {

    # use NNODE or $BATCH_NCORE to determin $NSPLIT,
    # the argument of split_version. Each VERSNION is 
    # split into $NSPLIT sub-versions for parallel processing.
    
    # Nov  4, 2011: fix bug to prevent NSPLIT=0
    # Feb 01, 2012: set logical OPT_SPLIT=0 if NSPLIT=1
    # Feb 05, 2012: min NSPLIT=1 instead of 2

    my ($NCORE, $XN, $key, @tmp );

    $NSPLIT    = 1 ;

    if ( $SSH_NNODE > 0 ) {  $NSPLIT  = $SSH_NNODE ;  }

    # check for batch system jobs
    $NCORE = $BATCH_NCORE ;
    if ( $NCORE > 0 ) {
	$NBATCH_TOTAL = $NVERSION * ( $NFITOPT-$NFITOPT_LINK ) ;
	$XN     = $NCORE / $NBATCH_TOTAL ;
	$NSPLIT = int($XN);
	if ( $NSPLIT < 2 ) { $NSPLIT = 1; }  # protect agains NSPLIT=0

	$OPT_SPLIT = 0 ;  # never do physical split on batch system

	# mak sure that batch-template script exists
	if ( !(-e $BATCH_TEMPLATE) ) {
	    $MSGERR[0] = "Cannot find BATCH_TEMPLATE file";
	    $MSGERR[1] = "$BATCH_TEMPLATE";
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}	
    }


    # Mar 2013: abort if NSPLIT exceed max number to merge
    if ( $NSPLIT > $MAXMERGE ) {
	$MSGERR[0] = "NSPLIT = $NSPLIT exceeds MAXMERGE=$MAXMERGE";

	if ( $SSH_NNODE )  { 
	    $MSGERR[1] = "Must define <= $MAXMERGE nodes." ; 
	}
	else { 
	    my $MAXCORE = $MAXMERGE * ( $NVERSION * $NFITOPT ) ;
	    $MSGERR[1] = "With NVERSION=$NVERSION and NFITOPT=$NFITOPT" ; 
	    $MSGERR[2] = "NCORE should be <= $MAXCORE ";
	} 
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    return ;
    
} # end of get_NSPLIT


# ============================================
sub check_BATCH_NCORE {

    # Created Dec 20 2017
    # Sanity checks on number of batch cores and NJOB_PER_BATCHFILE
    # (called after get_NSPLIT)
    
    if ( $SSH_NNODE ) { return ; }
    
    if ( $NSPLIT > 1 && $NJOB_PER_BATCHFILE > 1 ) {
	print " WARNING: NJOB_PER_BATCHFILE=$NJOB_PER_BATCHFILE -> 1 " .
	    "because NSPLIT=$NSPLIT>1 \n";
	$NJOB_PER_BATCHFILE = 1;
    }

    # ABORT if NJOB_PER_BATCHFILE is too big, resulting in using
    # fewer cores than reqested by user.
    my $NCORE_NEED = int($NBATCH_TOTAL/$NJOB_PER_BATCHFILE)+1 ;
    if ( $NCORE_NEED < $BATCH_NCORE && $NJOB_PER_BATCHFILE>1 ) {
	$MSGERR[0] = "Bad CORE request." ;
	$MSGERR[1] = "Requested $BATCH_NCORE cores, " . 
	    "but only need $NCORE_NEED cores.";
	$MSGERR[2] = "NJOB_PER_BATCHFILE = $NJOB_PER_BATCHFILE (probably too big)";
	$MSGERR[3] = "NBATCH_TOTAL       = $NBATCH_TOTAL " ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    return ;
    
} # end check_BATCH_NCORE


# ============================================
sub set_DOTEXT{
 
    # examing global SNTABLE_LIST and set global 
    # DOTEXT_TABLE logicals.

    # start with snana defaults.
    $DOTEXT_TABLE[$ITABLE_SNANA]  = 0;
    $DOTEXT_TABLE[$ITABLE_FITRES] = 1 ;
    $DOTEXT_TABLE[$ITABLE_LCPLOT] = 0 ;

    if ( $SNTABLE_LIST eq '' ) { goto DONE; }

    # check for text options
    my (@wdlist, $wd, $itab, $tbname);
    @wdlist     = split(/\s+/,$SNTABLE_LIST ) ;
    foreach $wd (@wdlist) {
	for($itab=0; $itab < $MAXTABLE; $itab++ ) {
	    $tbname = $TABLELIST[$itab] ;
	    if ( index($wd,$tbname) >= 0 ) {
		if ( index($wd,"text") >= 0 ) { $DOTEXT_TABLE[$itab] = 1; }
	    }
	}
    }

  DONE:
    for($itab=0; $itab < $MAXTABLE; $itab++ ) {
	$tbname = $TABLELIST[$itab] ;
	print "\t DOTEXT_TABLE[$tbname] = $DOTEXT_TABLE[$itab] \n";
    }

    $| = 1 ;


}  # end of set_DOTEXT

# ============================================
sub get_DATADIR(@) {

    my ($iver) = @_;

    # set global DATADIR to point to either the data area
    # or the SIM area depending on where the first VERSION resides.
    # Also set global SIM_FLAG = 1 or 0 to indicate data or sim.
    
    # Oct 22, 2012: check for PRIVATE_DATA_PATH option
    # Nov 11, 2014: check lcmerge/ and lcmerge/VERSION
    # Feb 26, 2015: fix logic to work with PRIVATE_DATA_PATH & SIM.

    my ($VERTMP, $SIMDIR, $DATADIR0, $DATADIR1, $KEY, $PATH ) ;

    $PRIVATE_FLAG_LIST[$iver] = 0 ;

    $PATH = "$PRIVATE_DATA_PATH" ;

    # Feb 10 2015: check 'pwd' option to substitute current dir.
    #   Allows using split_and_fit nml on different clusters
    #   with private data in `pwd`.

    if ( $PATH eq 'pwd' ) {
	$PATH = $CURRENT_DIR ;

	my ($sedDir, $SEDCMD) ;
	$sedDir = $CURRENT_DIR ;
	$sedDir =~ s{/}{\\/}g ;  # need \/ for sed command to work

	$SEDCMD = "sed -e 's/pwd/$sedDir/g' " ; 
	$nmlFile = "${nmlFile_ORIG}_MODIFIED" ;
	if ( -e $nmlFile ) { qx(rm $nmlFile); }
	qx($SEDCMD $nmlFile_ORIG > $nmlFile);
    }


    $VERTMP  = "$VERSION_LIST[$iver]";
    $SIMDIR  = "$PATH_SIMDATA_LIST[$iver]/${VERTMP}" ;

    if ( length($PATH) > 0 ) {
	unless ( -d $PATH ) {
	    $MSGERR[0] = "Cannot find PRIVATE_DATA_PATH = " ;
	    $MSGERR[1] = "$PATH" ;
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}

	$PRIVATE_FLAG_LIST[$iver] = 1 ;
	if ( $iver == 0 )  {  print " Use $KEY = \n\t $PATH \n"; }
	$DATADIR0 = "${PATH}" ;
	$DATADIR1 = "${PATH}/${VERTMP}" ;
    }
    else {
	$DATADIR0 = "${SNDATA_ROOT}/lcmerge" ;
	$DATADIR1 = "${SNDATA_ROOT}/lcmerge/${VERTMP}" ;
    }


    if ( $VERTMP eq '' ) {
	$MSGERR[0] = "get_DATADIR():";
	$MSGERR[1] = "VERSION is blank for iver=$iver" ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }


    
    if ( -d $SIMDIR ) {
	$SIM_FLAG_LIST[$iver] = $SIM_FLAG_NORMAL ;
	$DATADIR_LIST[$iver]  = "${SIMDIR}" ;
	$NVERSION_SIM++ ;

	# check fo SYMLINK
	$SIM_SYMLINKDIR_LIST[$iver] = "" ;
	my $lnkFile = "${SIMDIR}/$VERTMP.SYMLINK" ;
	my $lnkVer  = "" ;
	if ( -e $lnkFile ) {
	    $SIM_FLAG_LIST[$iver] = $SIM_FLAG_SYMLINK ;
	    $lnkVer          = qx(cat $lnkFile);
	    $lnkVer          =~ s/\s+$// ; 
	}
	$SIM_SYMLINKDIR_LIST[$iver] = "$lnkVer" ;
    }
    elsif ( -d $DATADIR1 ) {
	$SIM_FLAG_LIST[$iver] = 0 ;
	$DATADIR_LIST[$iver] = "${DATADIR1}" ;
	$NVERSION_DATA++ ;
    }
    elsif ( -d $DATADIR0 ) {
	$SIM_FLAG_LIST[$iver] = 0 ;
	$DATADIR_LIST[$iver] = "${DATADIR0}" ;
	$NVERSION_DATA++ ;
    }

    
    print " SIM_FLAG=$SIM_FLAG_LIST[$iver] for VERSION = $VERTMP \n" ;

} # end of get_DATADIR


# ========================
sub get_FORMAT(@) {

    my ($iver) = @_ ;

    # fill global @DATA_FORMAT = "TEXT" or "FITS"
    # based on extention of first file in LIST file.
    # If extension is ".FITS" then it's in FITS format;
    # else it's assumed to be in TEXT/ASCI format.

    my $fnam = "get_FORMAT" ;
    my ($VERSION, $PATH_DATA, @LIST_FILE, $simFlag );
    my (@bla, $firstFile, $j, $FORMAT );

    $VERSION   = "$VERSION_LIST[$iver]";
    $PATH_DATA = "$DATADIR_LIST[$iver]" ;
    $simFlag   = $SIM_FLAG_LIST[$iver];

    $LIST_FILE[0] = "$PATH_DATA/${VERSION}.LIST" ;
    $LIST_FILE[1] = "$PATH_DATA/${VERSION}/${VERSION}.LIST" ;

    if ( -e $LIST_FILE[0] ) 
    {  @bla = `cat $LIST_FILE[0]` ; }
    elsif ( -e $LIST_FILE[1] ) 
    {  @bla = `cat $LIST_FILE[1]` ; }
    else {
	$MSGERR[0] = "function $fnam cannot find LIST file for";
	$MSGERR[1] = "photometry version = '$VERSION' (iver=$iver)" ;
	$MSGERR[2] = "Tried" ;
	$MSGERR[3] = "  $LIST_FILE[0] and " ;
	$MSGERR[4] = "  $LIST_FILE[1]  " ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    $firstFile = $bla[0];
    $firstFile  =~ s/\s+$// ;   # trim trailing whitespace

    $j = index($firstFile,".FITS") ;

    if ( $j > 0 ) {
	$FORMAT = "FITS" ;
    }
    else {
	$FORMAT = "TEXT" ;
    }

    # load global
    $DATA_FORMAT[$iver] = "$FORMAT" ;

    # print " xxx get format for '$firstFile' => '$FORMAT' \n" ;

} # end of get_FORMAT


# ============================================
sub make_OUTDIR(@) {
    my ($opt,$iver) = @_;

    my ($VERSION, $DATADIR, $SIM_FLAG, $tmpDir, $cmd, $readme) ;
    my ($readmeFile, @readmeList, $miscDir, $inFile );

    # create separate output directory for each $VERSION
    # opt=0 => not a real version, just a subdir
    # opt=1 => real version, so copy README files
    #
    # July 19 2017: to reduce file output, copy FITOPT.README
    #               and input NML to OUTDIR only ... NOT to
    #               each VERSION subdir

    if ( $opt == 0 ) {
	$VERSION = "$SPLIT_JOBSUBDIR_LCFIT" ;
    }
    else {
	$VERSION  = "$VERSION_LIST[$iver]" ;
	$DATADIR  = "$DATADIR_LIST[$iver]" ;
	$SIM_FLAG = "$SIM_FLAG_LIST[$iver]" ;
    }

    $tmpDir = "$OUTDIR/$VERSION" ;
    print " Create fit-output dir: $tmpDir \n";
    if ( $opt ) { $OUTDIR_LIST[$iver] = "$tmpDir" ; }

    
    # -----------------------------------------------
    # Feb 11 2017
    # if DATADIR is a simulation with symbolic link to another
    # sim, create symbolic link for split_and_fit's OUTDIR
    # and avoid wasting CPU with re-fitting.
    if ( $opt &&  ($SIM_FLAG == $SIM_FLAG_SYMLINK)  )  { 
	if ( &make_OUTDIR_symLink($iver) > 0 ) 	{  
	    $SIM_FLAG_LIST[$iver] = $SIM_FLAG_SYMLINK2 ;
	    return;
	}  
    }

    # ----------------------------------------
    qx(mkdir -p $tmpDir);

    # xxx mark delete xxxx qx(cp $nmlFile  $tmpDir );
    if($opt==0) { 
	qx(cp $nmlFile  $tmpDir );  # OUTDIR/SPLIT_JOBS_LCFIT
	qx(cp $nmlFile  $OUTDIR ); 
    }


    if ( $nmlFile ne $nmlFile_ORIG ) 
    { qx(cp ${nmlFile_ORIG}  $tmpDir/${nmlFile_ORIG}_ORIG ); }

    sntools::checkDirExist($OUTDIR,"OUTDIR"); # $tmpComment);

    # copy misc input files if it's here (otherwise SNANA uses public default)
    if ( $opt == 0 ) {
	foreach $inFile ( @COPY_INPUTFILE_NAMES )  { 
	    if ( -e $inFile )  {  qx(cp $inFile $tmpDir ); }  
	} # inFile
    } # opt



    if ( $opt == 1 ) {

	# copy the FITOPT.README file into each version-dir	
	# xxx mark delete July 2017 qx(cp $FIT_README_FILE $tmpDir );

	# copy FITRES_COMBINE_FILE into current  dir if it exists.
	if ( -e $FITRES_COMBINE_FILE ) {
	    qx(cp $FITRES_COMBINE_FILE  $tmpDir );
	}	
    }  # end of opt if-block


    # copy SIM-README file(s) from SIM area only; don't bother 
    # for real data since the data .README files are more stable.

    if ( $opt == 1 &&  ($SIM_FLAG == $SIM_FLAG_NORMAL) ) {

	$readmeFile = "${DATADIR}/${VERSION}.README" ;

	$miscDir    = "${DATADIR}/misc" ;
	if ( -d $miscDir ) {
	    @readmeList = qx(ls ${miscDir}/*.README ) ;
	}
	@readmeList = ( @readmeList , $readmeFile ) ;

    }
	

    # Feb 7 2017: check ML option
    if ( $opt == 0 && $SALT2mu_INFILE ne "" ) {
	my $cp ;  
	$tmpDir = $SPLIT_JOBDIR_SALT2mu ;
	$cp     = "cp $SALT2mu_INFILE $tmpDir" ;
	qx(mkdir -p $tmpDir ; $cp);
	print " Create SALT2mu-output dir: $tmpDir \n";
	if ( $BATCH_NCORE > 0 ) { qx(cp $BATCH_TEMPLATE $tmpDir/); }
    }

    # 10/10/2011: check SALT2mu option
    if ( $opt == 0 && $SALT2mu_INFILE ne "" ) {
	my $cp ;
	$tmpDir = $SPLIT_JOBDIR_SALT2mu ;
	$cp     = "cp $SALT2mu_INFILE $tmpDir" ;
	qx(mkdir -p $tmpDir ; $cp);
	print " Create SALT2mu-output dir: $tmpDir \n";
	if ( $BATCH_NCORE > 0 ) { qx(cp $BATCH_TEMPLATE $tmpDir/); }
    }

    # Feb 2017: check ML-APPLY option
    if ( $opt == 0 &&  ($MLFLAG == $MLFLAG_APPLY) ) {
	qx(mkdir -p $SPLIT_JOBDIR_ML ) ;
	print " Create MLAPPLY output dir: $SPLIT_JOBDIR_ML \n";
    }

    # Sep 2013: make separate MERGE log-dir under SPLIT_JOBS/MERGELOGS
    if ( $opt == 0 ) { qx(mkdir $MERGE_LOGDIR); }

    return ;

} # end of make_OUTDIR


# =========================================
sub make_OUTDIR_symLink(@) { 

    my ($iver)  = @_ ;

    my $VERSION     = "$VERSION_LIST[$iver]" ;
    my $DATADIR     = "$DATADIR_LIST[$iver]" ;
    my $LNKVER      = "$SIM_SYMLINKDIR_LIST[$iver]" ;
    my $SIMPATH     = "$PATH_SIMDATA_LIST[$iver]";
    my $README_FILE = "$SIMPATH/${LNKVER}/${LNKVER}.README" ;

    my $key  = "OUTDIR(split_and_fit):" ;
    my @tmp  = sntools::parse_line($README_FILE, 1, $key, $OPT_QUIET ) ;
    if ( scalar(@tmp) == 0 ) { return(0); }

    my $PREVIOUS_OUTDIR = "$tmp[0]" ;
    qx(cd $OUTDIR; ln -s $PREVIOUS_OUTDIR $VERSION);

#    print "\n"; 
#    print " xxx PREVIOUS_OUTDIR = $PREVIOUS_OUTDIR \n";
#    print " xxx README_FILE     = $README_FILE \n";
#    die " xxx DEBUG ABORT  \n";   

    return 1 ;    

} # end make_OUTDIR_symLink

# =========================================
sub init_SPLITSCRIPT(@) {

    my ($isplit) = @_;
    my ($node, $inode, $SPLITnnn, $script, $SCRIPT);

    $inode   = $isplit - 1;
    if ( $SSH_NNODE ) {
	$node    = $SSH_NODELIST[$inode];
	$node       =~ s/\s+$// ;   # trim trailing whitespace
    }
    elsif ( $BATCH_NCORE > 0 ) {
	$node = "batch" ;
    }
    else {
	$node = "interactive" ;
    }
    

    $SPLITnnn  = sprintf("SPLIT%3.3d",$isplit);

    if ( $BATCH_NCORE > 0 ) {
	# one script with all batch calls
	# xxx mark delete 	$script    = "FITSPLIT_${node}" ;
	$script    = "FITSPLIT_ALL_${node}" ; 
    }
    else {
	# one script per ISPLIT
	$script    = "FIT${SPLITnnn}_${node}" ; 
    }

    $SCRIPT    = "${SPLIT_JOBDIR_LCFIT}/${script}" ; 

    unless ( -e $SCRIPT ) {
	qx(touch $SCRIPT);
	qx(echo ${SNANA_LOGIN_SETUP} >> ${SCRIPT} );
    }

    $SCRIPT_LIST[$isplit] = $SCRIPT ;
    $script_LIST[$isplit] = $script ;
    $INODEMAP[$isplit]    = $inode ;
    $NODEMAP[$isplit]     = $node ;

} # end of init_SPLITSCRIPT

# ===================================
sub open_README {

    # create README file and give FITOPT list
    # Jan 29 2017: include snana version in FITOPT.README

    my ($ifitopt, $c3opt, $FITOPT );

    $FIT_README_FILE = "${OUTDIR}/FITOPT.README" ;
    open READMEFILE , "> $FIT_README_FILE" ;

    print READMEFILE 
	"# USER FIT OPTIONS with SNANA $CURRENT_SNANA_VERSION \n" ;

    for ( $ifitopt=$ifitopt_start; $ifitopt < $NFITOPT; $ifitopt++ ) {
	$c3opt  = sprintf("%3.3d", $ifitopt );
#	$FITOPT = "$FITOPT_LIST[$ifitopt]" ; 
	$FITOPT = "$FITOPT_LIST_ORIG[$ifitopt]" ; 
	print READMEFILE "FITOPT: $c3opt $FITOPT \n";
    }

    if ( -e $FITRES_COMBINE_FILE ) {
	print READMEFILE "FITRES_COMBINE_FILE: $FITRES_COMBINE_FILE \n";
    }

    print READMEFILE "\n";
    close READMEFILE ;

} # end of open_README


# ===================================
sub check_snana_version {

    # Aug 2015
    # if MIN_SNANA_VERSION key is set, then make sure that
    # current snana version >= MIN_SNANA_VERSION; 
    # otherwise abort.
    #
    # Dec 10 2017: get version from sntools.h instead of snana.car
    #
    
    my ($key, @tmp, $snanaFile, @wdlist ) ;

# xxxxxxxxx mark delete Dec 10 2017 xxxxxxxxxx
#    $snanaFile = "$SNANA_DIR/src/snana.car" ;
#    @tmp = qx(grep $qq SNANA_VERSION $qq $snanaFile | grep =);
#    @wdlist     = split(/\s+/,$tmp[0] ) ;
# xxxxxxxxxxxxxxxxxxxxxxxxx

    $snanaFile = "$SNANA_DIR/src/sntools.h" ;
    @tmp = qx(grep $qq SNANA_VERSION_CURRENT $qq $snanaFile | grep define );
    @wdlist     = split(/\s+/,$tmp[0] ) ;

    $CURRENT_SNANA_VERSION = $wdlist[2];
    $CURRENT_SNANA_VERSION =~ s/$qq//g ;  # remove double quotes

    $key  = "MIN_SNANA_VERSION:" ;
    @tmp  = sntools::parse_array($key,1, $OPT_QUIET, @CONTENTS_NMLFILE);
    if ( scalar(@tmp) > 0 ) { 
	$MIN_SNANA_VERSION = "$tmp[0]"; 
	
	if ( $CURRENT_SNANA_VERSION ge $MIN_SNANA_VERSION ) {
	    print " SNANA_VERSION = $CURRENT_SNANA_VERSION " . 
		">= $MIN_SNANA_VERSION \n";
	}
	else {
	    $MSGERR[0] = "Current SNANA version = $CURRENT_SNANA_VERSION";
	    $MSGERR[1] = "but $MIN_SNANA_VERSION or higher is required." ;
	    $MSGERR[2] = "See MIN_SNANA_VERSION key in $nmlFile";
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}
    }

} # end check_snana_version

# ===================================
sub split_version(@) {

    my ($iver) = @_;

    my ($VERSION, $cmdsplit, @tmp, $istat, $FORMAT, $XARG ) ;


    $FORMAT = $DATA_FORMAT[$iver] ;

    if ( $FORMAT eq "FITS" ) { return ; }
    if ( $OPT_SPLIT == 0   ) { return ; }

    $VERSION = $VERSION_LIST[$iver] ;

    $XARG = "" ;
    if ( $PRIVATE_FLAG_LIST[$iver] == 1 ) {
	$XARG = "PRIVATE_DATA_PATH  $DATADIR_LIST[$iver]" ;
    }

    # split version into subset for each node
    $cmdsplit = "$JOBNAME_SPLIT $VERSION $NSPLIT  $XARG";
    print " Split '${VERSION}' into $NSPLIT sub-versions (one per node). \n";
    @tmp = qx($cmdsplit) ;
    @tmp = qx(echo $?);
    $istat = $tmp[0] ;

    # abort if split could not be done.
    if ( $istat != 0 ) {
	$MSGERR[0] = "Unable to split version '$VERSION' " ;
	$MSGERR[1] = "Previous SPLIT*${VERSION}* probably still exists.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

} # end of split_version

# =========================================
sub createJobs_LCFIT(@) {

    # create script for light curve fit.
    # Jul 23 2013: add text-dump args if LDMP_SNFLUX=T
    # Sep 12 2014: add $LEGACY_ARGS_OFF to arg list
    #
    # Apr 06 2015: FIRST logic for header commands, and 
    #              new delay logic
    #
    # Jun 20 2016: remove LEGACY_ARGS_OFF
    # Sep 14 2016: remove optional comments after FITOPT
    #             (anything starting with # or !)

    my ($iver, $ifitopt, $isplit) = @_;

    my ($prefixJob,  $prefixBatch, $SPLITnn, $script, $FITOPT, $FITOPT_VER, $FORMAT);
    my ($FITOPT_LINK, $IFITOPT_LINK, $SIM_FLAG );
    my ($VERSION, $RFILE, $VARG,  $FARG, $SPLARG, $MXARG );
    my ($indx, $OUTFILE, @OUTARG, $DMPARG );
    my ($logFile, $doneFile, $busyFile);
    my ($tmp_prefix, $batchFile, $batchLog, $NTMP );
    my ($CMD_STRING, $FITJOB, $CMD, $FIRST, $KEY, @tmp, $tmp );
    my ($CMD_rmlocf, $CMD_cleanFITRES ) ;
    my ($UPDMASK_CMD, $APPENDFLAG_CMD, $WRITEFLAG_CMD, $icmd); 
    my $LDMP = 0 ;

    if ( $LDMP ) {
	print " xxx ----------------------------------- \n";
	print " xxx Begin createJobs_LCFIT  \n";
	print " xxx iver,ifitopt,isplit = $iver, $ifitopt, $isplit \n";
    }

    # ----------------------------------------------
    # Dec 19 2017: check if batch file is new, or append job    
    $UPDMASK_CMD    = &updmask_CMD_LCFIT_ARRAY();
    $APPENDFLAG_CMD = ( $UPDMASK_CMD & 1 );
    $WRITEFLAG_CMD  = ( $UPDMASK_CMD & 2 );

    if ( $APPENDFLAG_CMD == 0 )  { 
	$icmd = -1 ; 
	@CMD_LCFIT_ARRAY = ();
    }    
    else {
	$icmd = scalar(@CMD_LCFIT_ARRAY) - 1 ;
    }

    # - - - -
    $prefixJob   = &get_PREFIX_SPLITJOB($iver,$ifitopt,$isplit);
    $prefixBatch = &get_PREFIX_BATCHFILE($NBATCH_CREATE_FILE+1);
 
    $SPLITnn  = sprintf("SPLIT%2.2d",$isplit);
    $script   = $SCRIPT_LIST[$isplit] ;
    $FORMAT   = $DATA_FORMAT[$iver];
    
    $tmp_prefix = "" ; 
    if ( $FORMAT eq "TEXT" && $OPT_SPLIT ) 
    {  $tmp_prefix = "${SPLITnn}_" ; }


    $VERSION      = "${tmp_prefix}$VERSION_LIST[$iver]" ;
    $FITOPT_VER   = "$VERSION_FITOPT_LIST[$iver]" ;
    $FITOPT       = "$FITOPT_LIST[$ifitopt]" ;
    $VARG         = "VERSION_PHOTOMETRY $VERSION" ;
    $FITOPT_LINK  = "$FITOPT_LINK_LIST[$ifitopt]" ;
    $IFITOPT_LINK = "$IFITOPT_LINK_LIST[$ifitopt]" ;
    $SIM_FLAG     = "$SIM_FLAG_LIST[$iver]" ;

    # Sep 14 2016: remove optional comments at end of FITOPT
    $FITOPT     =~ s/\#.*//;  # remove # and everything after
    $FITOPT     =~ s/\!.*//;  # remove ! and everything after
    $FITOPT_VER =~ s/\#.*//;
    $FITOPT_VER =~ s/\!.*//;

    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) {

	if ( $indx < $OUTINDX_TEXT ) 
	{ $OUTFILE   = "${prefixJob}.$SUFFIX_OUTFILE[$indx]" ; }
	else
	{ $OUTFILE   = "${prefixJob}" ; }  # just prefix for text
	
	$OUTARG[$indx] = "" ;
	if ( $OUTFLAG[$indx] )  
	{ $OUTARG[$indx] = "$NMLKEY_OUTFILE[$indx]  ${OUTFILE}" ; }
    }


    # Jul 23 2013: check for lightcurve text-dump option
    $DMPARG = "" ;

    $SPLARG = "" ;
    if ( $FORMAT eq "FITS" || $OPT_SPLIT == 0 ) 
    { $SPLARG = "JOBSPLIT $isplit $NSPLIT" ; }
    else
    { $SPLARG = "JOBSPLIT_EXTERNAL $isplit $NSPLIT" ; }
    
    $MXARG = "";
    if ( $MXLC_PLOT eq "" ) 
    { $MXARG = "MXLC_PLOT $MXLC_PLOT_DEFAULT" ; }
    

    $logFile  = "${prefixJob}.${SUFFIX_LOG}" ;
    $doneFile = "${prefixJob}.${SUFFIX_DONE}" ;
    $busyFile = "${prefixBatch}.${SUFFIX_BUSY}" ;

    if ( $LDMP ) {  print " xxx update script = $script \n"; }

    $NFITJOB_TOT++;
    if ( $IFITOPT_LINK < 0 ) { $NFITJOB_SUBMIT++ ; }

    # - - - - - 
    # construct array of commands @CMD_LCFIT_ARRAY
    if ( $SNANA_MODELPATH ne "" ) { 
	$icmd++ ;  @CMD_LCFIT_ARRAY[$icmd] = "$ENVDEF_MODELPATH" ;
    }

    if ( $IFITOPT_LINK >= 0  || ($SIM_FLAG==$SIM_FLAG_SYMLINK2) ) { 
	qx(cd $SPLIT_JOBDIR_LCFIT; echo GRACE > $logFile ; touch $doneFile );
	$NDONE_LINK++ ;
	return ;
    }

    # prepare fit jobs unless it's a sim pointing to a symLink

    $FITJOB = 
	"$JOBNAME_LCFIT " . 
	"$nmlFile $VARG @OUTARG $DMPARG $SPLARG $MXARG " ;
    if ( $FITOPT     ne "NONE" ) { $FITJOB = "$FITJOB  $FITOPT" ;     }
    if ( $FITOPT_VER ne ""     ) { $FITJOB = "$FITJOB  $FITOPT_VER" ; }

# If executable does not exist (e.g., somebody did a 'make')
# then sleep for a few minutes.
    $icmd++ ; @CMD_LCFIT_ARRAY[$icmd] =
	"wait_for_files.pl 1 $JOBNAME_LCFIT_FULLPATH NOSTAMP";
###	"if [ ! -f $JOBNAME_LCFIT_FULLPATH ] ; then sleep 200 ; fi ;" ;

# command to run executable with arguments
    $icmd++ ;  @CMD_LCFIT_ARRAY[$icmd] = "$FITJOB >& $logFile" ;

# create DONE file
    $icmd++ ;  @CMD_LCFIT_ARRAY[$icmd] = "touch  $doneFile" ;

# remove BUSY file created by batch-submit script (Dec 12, 2017)
    if ( $BATCH_NCORE>0 && $WRITEFLAG_CMD ) 
    { $icmd++ ;  @CMD_LCFIT_ARRAY[$icmd] = "rm     $busyFile" ; }
    
    # remove annoying locf messages with 64-bit CERNLIB
    # (does nothing for 32-bit CERNLIB)
    if ( $OUTFLAG[$OUTINDX_HBOOK] ) { 
	$CMD_rmlocf = "remove_locf_messages.pl $logFile QUIET" ;
	$icmd++ ; @CMD_LCFIT_ARRAY[$icmd] = "$CMD_rmlocf" ;
    }

  DONE_CMD:

    # construct CMD_STRING for batch system.
    my $NCMD = scalar(@CMD_LCFIT_ARRAY) ;
    $CMD_STRING = "$CMD_LCFIT_ARRAY[0]" ;
    for($icmd = 1; $icmd < $NCMD; $icmd++ ) 
    { $CMD_STRING = "$CMD_STRING ; $CMD_LCFIT_ARRAY[$icmd]" ;  }
    
    # ----------------------------------------------------
    # check if this is FIRST entry in script    
    $KEY = $SPLIT_JOBDIR_LCFIT ;  # Dec 20 2017
    @tmp = qx(grep $KEY $script );
    $FIRST = ( scalar(@tmp) == 0 ) ;
   
    open PTR_SCRIPT , ">> $script" ;  

    # write one-time command(s) at time of file
    if ( $FIRST ) {
      print PTR_SCRIPT "cd $SPLIT_JOBDIR_LCFIT \n" ;
      if ( $SNANA_MODELPATH ne "" ) 
      { print PTR_SCRIPT "$ENVDEF_MODELPATH \n";}
      print PTR_SCRIPT "\n" ;
    }


    if ( $BATCH_NCORE > 0 ) {

	$NBATCH_CREATE_JOB++ ;  # increment number or executable jobs
      
	if ( $isplit == 1 && $ifitopt == 0 ) {
	    $NTMP = ($NFITOPT-$NFITOPT_LINK) * $NSPLIT ;
	    print "\t Preparing $NTMP batch jobs; please be patient ...\n";
	}

	if ( $WRITEFLAG_CMD ) {
	    $NBATCH_CREATE_FILE++ ; # each batch file --> one batch submission

	    $prefixBatch = &get_PREFIX_BATCHFILE( $NBATCH_CREATE_FILE );
	    $batchFile = "${prefixBatch}.BATCH" ;
	    $batchLog  = "${prefixBatch}.BATCH-LOG" ;
	    
	    print PTR_SCRIPT "\n#Submit batch file number $NBATCH_CREATE_FILE\n";
	    &batch_delayAdd($iver,$ifitopt,$isplit,$NBATCH_CREATE_FILE) ;
	    print PTR_SCRIPT "touch $SPLIT_JOBDIR_LCFIT/$busyFile\n";
	    print PTR_SCRIPT "$BATCH_COMMAND $batchFile \n" ;
	    
	    #  foreach $tmp ( @CMD_LCFIT_ARRAY ) { print " xxx CMD=$tmp \n"; }
	    #  die "\n xxx DEBUG DIE xxx\n" ;
	    
	    sntools::make_batchFile($BATCH_TEMPLATE, $SPLIT_JOBDIR_LCFIT, 
				    $batchFile, $batchLog, $BATCH_MEM_LCFIT, 
				    $CMD_STRING);
	}
    }
    else {
	foreach $CMD ( @CMD_LCFIT_ARRAY ) { print PTR_SCRIPT "$CMD\n\n" ; }
    }
    
    close PTR_SCRIPT ;
    

    if ( $LDMP ) {
	print " xxx \t(iver,ifitopt,isplit=$iver,$ifitopt,$isplit->$SPLITnn)\n";
    }
    

    return ;
    
} # end of createJobs_LCFIT


# ======================================
sub updmask_CMD_LCFIT_ARRAY  {

    # Created Dec 19 2017
    # Set bit 0 to append CMD_LCFIT_ARRAY
    # Set bit 1 to write out BATCH file.
    
    my ($MASK, $NEXT_NBATCH_CREATE, $NMOD);
    my $MASK_APPEND = 1;
    my $MASK_WRITE  = 2;

    # never append jobs for ssh
    if ( $SSH_NNODE ) { return($MASK_WRITE); }

    # never append if 1 job per batch file (this is most common mode)
    if ( $NJOB_PER_BATCHFILE == 1 ) { return($MASK_WRITE); }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # if we get here, combine multuple jobs in one batch file
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    $MASK = 0 ;
    $NEXT_NBATCH_CREATE = $NBATCH_CREATE_JOB + 1;

    if ( $NEXT_NBATCH_CREATE == $NBATCH_TOTAL ) 
    { $MASK |= $MASK_WRITE; }

    $NMOD = ($NEXT_NBATCH_CREATE % $NJOB_PER_BATCHFILE) ;
    if ( $NMOD == 0 )   { $MASK |= $MASK_WRITE; }
    if ( $NMOD != 1 )   { $MASK |= $MASK_APPEND; }

    return($MASK) ;

} # end updmask_CMD_LCFIT_ARRAY

# ======================================
sub batch_delayAdd {
    my ($iver, $ifitopt, $isplit, $ibatch) = @_ ;

    # Created July 17 2017
    # check for delay(s) to add before launching next batch job
    # Write sleep command to PTR_SCRIPT.
    # Logic is tricky; submit-script creates BUSY file, while batch 
    # job removes BUSY file. This fragile code is needed to work
    # in cases where batch queues stall for extended periods rather
    # than launch right away.
    #
    # Dec 12 2017: 
    #   switch from DONE logic to BUSY logic to limit the number 
    #   of files being checked by wait_for_files.
    #
    my ($prefixJob, $prefixBatch, $busyFile, $delay, $j, $versFitopt );
    my ($cmdWait, $cmdBusy, $BUSYLIST, $DONELIST, $NDONE_REQ);
    my $NWAIT=0;

    $prefixJob    = &get_PREFIX_SPLITJOB($iver,$ifitopt,$isplit);

    # if there are more jobs than requested NCORE, then
    # put in a delay to wait for BUSY files.
    if ( $BATCH_NCORE        >  0            &&
	 $NBATCH_CREATE_FILE >  $BATCH_NCORE && 
	 $NBATCH_CREATE_JOB  <  $NBATCH_TOTAL ) {

	$prefixBatch  = &get_PREFIX_BATCHFILE($ibatch);    
	$busyFile     = "${prefixBatch}.${SUFFIX_BUSY}" ;

	$BUSYLIST   = "$SPLIT_JOBDIR_LCFIT/\\*.BUSY";
	$cmdWait    = "wait_for_files.pl -$BATCH_NCORE $BUSYLIST NOSTAMP" ;
	$cmdBusy    = "touch $SPLIT_JOBDIR_LCFIT/$busyFile";

	print PTR_SCRIPT "$cmdWait\n";  
    }

    
    if ( $DELAY_SUBMIT[0] > 0  || $DELAY_SUBMIT[1] > 0 ) { 

	$delay = $DELAY_SUBMIT[0] ;
	# full user-delay between each VERSION and FITOPT
	if ( $isplit==$NSPLIT && $DELAY_SUBMIT[1]>0 ) 
	{ $delay = $DELAY_SUBMIT[1]; }

	print PTR_SCRIPT "sleep $delay \n"; 

	if ( $isplit == $NSPLIT ) { 
	    $j          = index($prefixJob,"SPLIT");
	    $versFitopt = substr($prefixJob,0,$j-1);
	    print PTR_SCRIPT 
		"echo ${qq}Done submitting $versFitopt\n${qq}\n"; 
	}
    }

    return ;

}  # end batch_delayAdd


# =====================================
sub createJobs_MLAPPLY(@) {

    # Created Feb 7 2017
    # Create scripts to 
    #  + run machine-learning (ML) job to classify SN,
    #  + append ML outputs (e.g., Prob_Ia) to original tables.
    #

    my ($iver, $ifitopt) = @_ ;
    
    my $VERSION      = "$VERSION_LIST[$iver]" ;
    my $FITOPT       = "$FITOPT_LIST[$ifitopt]" ;
    my $MERGEDIR     = "$OUTDIR_LIST[$iver]";
    my $PREFIX_ORIG  = sprintf("FITOPT%3.3d", $ifitopt );
    my $PREFIX_OUT   = "MLAPPLY_${VERSION}_${PREFIX_ORIG}" ;
    my $PREFIX_FINAL = "../${VERSION}/${PREFIX_ORIG}+ML" ;

    # get NOMINAL (i.e.,, original) fitres file to apply ML
    my $FRES_NOMINAL = "${PREFIX_ORIG}.${SUFFIX_FITRES}" ;

    my (@CMD, $NCMD, $inFile_data, $inFile_MLpar, $outFile_ML );
    my ($logFile, $doneFile );

    # all MLAPPLY jobs must take standard arguments:
    # 1) input data file, 2) MLpar file, 3) outFile    

    $inFile_data  = "../$VERSION/$FRES_NOMINAL" ;
    $inFile_MLpar = "$MLAPPLY_INFILE[$ifitopt]" ;

    # write outFile in ROOT format because it's x10 faster
    # than writing to TEXT format
#    $outFile_ML   = "${PREFIX_OUT}.${SUFFIX_FITRES}" ;
    $outFile_ML   = "${PREFIX_OUT}.$SUFFIX_OUTFILE[$OUTINDX_ROOT]" ;

    $logFile      = "${PREFIX_OUT}.LOG" ;    # stdout
    $doneFile     = "${PREFIX_OUT}.DONE" ; 
    $NCMD         =  0 ;

    $CMD[$NCMD] = 
	"$JOBNAME_MLAPPLY " .
	"   -inFile_data $inFile_data  "  .
	"   -inFile_MLpar $inFile_MLpar  " .
	"   -outFile $outFile_ML " .
	"   >& $logFile  ";
	;
    $NCMD++ ;
 
    # append ML variables to each output file format
    # that exists (HBOOK,ROOT,TEXT)

    my($FITRES_orig, $FITRES_combine, $SUFFIX, $indx);

    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) {
	if ( $OUTFLAG[$indx] == 0 ) { next ; }

	$SUFFIX         = "$SUFFIX_OUTFILE[$indx]" ;

	# for TEXT, switch suffix to FITRES
	if ( $indx == $OUTINDX_TEXT ) { $SUFFIX = $SUFFIX_FITRES ; }

	$FITRES_orig    = "../${VERSION}/${PREFIX_ORIG}.${SUFFIX}" ;
	$FITRES_combine = "COMBINE_${PREFIX_OUT}.${SUFFIX}" ;
	$logFile        = "COMBINE_${PREFIX_OUT}_${SUFFIX}.LOG" ;
	$CMD[$NCMD] = 
	    "sntable_combine.exe " .
	    "  -inFile $FITRES_orig $outFile_ML " .
	    "  -outFile $FITRES_combine " .
	    "  >& $logFile " ;
	$NCMD++ ;

	# clobber original FITRES table file with appended one
	$CMD[$NCMD] = "mv $FITRES_combine $FITRES_orig " ;
	$NCMD++ ;
    }
    
    $CMD[$NCMD] = "touch $doneFile" ;
    $NCMD++ ;

    &create_batchScript(\$INODE_ML, $PREFIX_RUNML,
			$SPLIT_JOBDIR_ML, $VERSION, $PREFIX_ORIG, 
			$BATCH_MEM_MLAPPLY, @CMD);

    return ;

} # end createJobs_ML


# =====================================
sub createJobs_SALT2mu(@) {
    
    # Create script(s) for SALT2mu jobs.
    # May 2013: major overhaul
    #
    # If FITOPT[nnn].FITRES is the original fitres file
    # from the fitting job, then the scripts here will 
    # create
    #   FITOPT[nnn]+SALT2mu.FITRES
    #   FITOPT[nnn]+SALT2mu.HBOOK   (if hbook format used)
    #   FITOPT[nnn]+SALT2mu.ROOT    (if root  format used)
    #
    # The hbook/root format options are taken from the fitting options
    # based on user input HFILE_OUT and ROOTFILE_OUT.
    # 
    # Nov  6 2013: fix small bug setting $fmtArg to avoid "?" in arg list
    # Sep 13 2014: add SNANA_LOGIN_SETUP to run-script
    # Feb 18 2019: move M0DIF file along with FITRES file
    # ----------

    my ($iver, $ifitopt) = @_ ;

    my ($cdd, $lncmd, $node, $script, $NCMD, @CMD, $cmd, $SIMFILE );
    my ($CMD_lnk, $CMD_S2mu, @CMD_mv, $CMD_clean);
    my ($fArg, $pArg, $simArg, $fmtArg, $Nmove, $indx, $FF);
    my ($icmd, $mv, $outF1, $outF2, $outM1, $outM2 );
    my ($batchFile, $batchLog, $batchJob , $sep );

    # ------------- BEGIN -------------

    my $VERSION      = "$VERSION_LIST[$iver]" ;
    my $FITOPT       = "$FITOPT_LIST[$ifitopt]" ;
    my $MERGEDIR     = "$OUTDIR_LIST[$iver]";
    my $PREFIX_ORIG  = sprintf("FITOPT%3.3d", $ifitopt );
    my $PREFIX_S2mu  = "S2mu_${VERSION}_${PREFIX_ORIG}" ;
    my $PREFIX_FINAL = "../${VERSION}/${PREFIX_ORIG}+SALT2mu" ;

    # get NOMINAL (i.e.,, original) fitres file to fit with SALT2mu
    my $FRES_NOMINAL = "${PREFIX_ORIG}.${SUFFIX_FITRES}" ;

    # define misc log/done files
    my $logFile     = "${PREFIX_S2mu}.${SUFFIX_LOG}" ;
    my $doneFile    = "${PREFIX_S2mu}.${SUFFIX_DONE}" ;
    my $combLogFile  = "COMBINE_${PREFIX_S2mu}.${SUFFIX_LOG}" ;

    # Construct command for  symbolic link to fitres file in MERGEDIR 
    # so that we have a local input FITRES file in the SALT2mu dir.
    my $FRES_LINK  = "ORIG_${VERSION}_${PREFIX_ORIG}.${SUFFIX_FITRES}" ;
    my $CMD_lnk    = "ln -s ../$VERSION/$FRES_NOMINAL  $FRES_LINK" ;

    # construct SALT2mu command
    $fArg = "file=$FRES_LINK" ;
    $pArg = "prefix=${PREFIX_S2mu}" ;

    # check for option simArg for SNIa biasCor and CC prior
    $simArg = 
	&get_simArg1_SALT2mu($iver,$ifitopt) .  # SALT2mu_SIMVERSION key
	&get_simArg2_SALT2mu($iver,$ifitopt) ;  # SALT2mu_BIASCOR_PATH key

    $CMD_S2mu  = "$JOBNAME_SALT2mu $SALT2mu_INFILE $fArg $pArg $simArg" ;

    $Nmove = 0 ;  # number of move-related commands

    # move the final SALT2mu fitres text file, and leave symb. link
    $outF1 = "${PREFIX_S2mu}.fitres" ;
    $outF2 = "${PREFIX_FINAL}.${SUFFIX_FITRES}" ;
    $CMD_mv[$Nmove] = "mv  $outF1  $outF2"   ;   $Nmove++ ;
    $CMD_mv[$Nmove] = "ln -s  $outF2 $outF1" ;   $Nmove++ ;

    # Feb 18 2019: repeat for M0DIF file
    $outM1 = "${PREFIX_S2mu}.M0DIF" ;
    $outM2 = "${PREFIX_FINAL}.M0DIF" ;
    $CMD_mv[$Nmove] = "mv  $outM1  $outM2"   ;   $Nmove++ ;
    $CMD_mv[$Nmove] = "ln -s  $outM2 $outM1" ;   $Nmove++ ;

    # move each his/root file and leave symb. link
    $fmtArg = "" ; # init format specifier for combine_fitres.exe
    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) {
	if ( $OUTFLAG[$indx] ) {
	    if ( $indx == $OUTINDX_TEXT ) { next ; }
	    my $SUFFIX = "$SUFFIX_OUTFILE[$indx]" ;
	    $fmtArg = "$fmtArg $COMBINE_FMTARG[$indx]" ;  # H and/or R ..
	    $outF1  = "${PREFIX_S2mu}.${SUFFIX}" ;  
	    $outF2  = "${PREFIX_FINAL}.${SUFFIX}" ;
	    $CMD_mv[$Nmove] = "mv  $outF1  $outF2"   ; $Nmove++ ;
	    $CMD_mv[$Nmove] = "ln -s $outF2  $outF1" ; $Nmove++ ;
	}
    }  # indx

    # construct command to translate text output into his/root file(s)
    my $OUTARG    = "--outprefix ${PREFIX_S2mu}" ; 
    my $FRES_OUT  = "${PREFIX_S2mu}.fitres" ;
    my $CMD_COMB  = "$JOBNAME_COMBINE $FRES_OUT $fmtArg $OUTARG" ;

    # cleanup 
    $CMD_clean = "rm $combLogFile" ;

    # make list of all commands -> @CMD
    $NCMD = 0 ;
    $NCMD++ ;  $CMD[$NCMD-1] = "$CMD_lnk" ;
    $NCMD++ ;  $CMD[$NCMD-1] = "$CMD_S2mu >& $logFile " ;
    $NCMD++ ;  $CMD[$NCMD-1] = "$CMD_COMB >& $combLogFile" ;
    foreach $mv (@CMD_mv) { $NCMD++ ; $CMD[$NCMD-1] = "$mv"; }
    $NCMD++ ;  $CMD[$NCMD-1] = "$CMD_clean" ;
    $NCMD++ ;  $CMD[$NCMD-1] = "touch $doneFile" ;
    
    # ----------------------------

    &create_batchScript(\$INODE_SALT2mu, $PREFIX_RUNSALT2mu, 
			$SPLIT_JOBDIR_SALT2mu, $VERSION, $PREFIX_ORIG, 
			$BATCH_MEM_SALT2mu, @CMD);

    return ;

} # end of createJobs_SALT2mu


# =========================================
sub get_simArg1_SALT2mu {

    my ($iver,$ifitopt) = @_ ;

    # Created Feb 19 2016 by R.Kessler
    #
    # if SALT2mu_SIMVERSION key is specified, return command-line
    # arguments for SALT2mu of the form
    #   simfile_ccprior=$SIMFILE  simfile_bias=same" ;
    #
    # where $SIMFILE is the fitres filename corresponding
    # to SALT2mu_SIMVERSION. If SALT2mu_SIMVERSION is an integer 'n'
    # then use the n'th VERSION in the list (1st VERSION starts at 1).
    #
    
    my ($simArg, $SIMFILE, $PREFIX, $VERSION);

    $simArg = "" ;  # init to nothing
    if ( length($SALT2mu_SIMVERSION) == 0 ) { return($simArg); }

    if ( (1*$SALT2mu_SIMVERSION) == $SALT2mu_SIMVERSION ) {
	# use integer to look up version
	my $IVER = $SALT2mu_SIMVERSION - 1;
	$VERSION = "$VERSION_LIST[$IVER]" ;
    }
    else {
	# use explicit version
	$VERSION = "${SALT2mu_SIMVERSION}" ;
    }

    # point to actual file instead of symbolic link because
    # the ORIG_ symLink may not exist when job runs.
    $PREFIX   = sprintf("FITOPT%3.3d", $ifitopt );
    $SIMFILE  = "../${VERSION}/${PREFIX}.${SUFFIX_FITRES}" ;
    $simArg   = "simfile_ccprior=$SIMFILE  simfile_bias=same" ;

    return($simArg);

} # end get_simArg1_SALT2mu


# =========================================
sub get_simArg2_SALT2mu {

    # If SALT2mu_BIASCOR_PATH is defined
    
    my ($iver,$ifitopt) = @_ ;
    
    my ($PREFIX, $VERSION, $SIMFILE, $simArg );

    $simArg = "" ;  # init to nothing
    $PREFIX  = sprintf("FITOPT%3.3d", $ifitopt );

    if ( length($SALT2mu_BIASCOR_PATH) == 0 ) {  return $simArg ; }

    $SIMFILE  = "$SALT2mu_BIASCOR_PATH/${PREFIX}.${SUFFIX_FITRES}" ;
    $simArg   = "$simArg simfile_bias=$SIMFILE" ;
    
    if ( $SALT2mu_CCPRIOR_PATH eq $SALT2mu_BIASCOR_PATH ||
	 $SALT2mu_CCPRIOR_PATH eq "same" ) {
	$simArg   = "$simArg simfile_ccprior=same" ;
    }
    elsif ( length($SALT2mu_CCPRIOR_PATH) > 0 ) {
	$SIMFILE  = "$SALT2mu_CCPRIOR_PATH/${PREFIX}.${SUFFIX_FITRES}" ;
	$simArg   = "$simArg simfile_ccprior=$SIMFILE" ;
    }

    return($simArg);
    
} # end get_simArg2_SALT2mu

# ============================
sub create_batchScript(@) {

    my($INODE, $PREFIX, $SPLIT_JOBDIR, $VERSION, 
       $FFF, $BATCH_MEM, @CMD) = @_ ;

    # Created Feb 8 2017
    # create generic batch scripts for any post-LCFIT job
    # (e.g., SALT2mu, ML)
    # This script is NOT intended for LCFIT jobs.
    #
    #   VERSION = data or sim version
    #   FFF     = FITOPTnnn

    my ($script, $batchFile, $batchLog, $batchJob, $cmd, $icmd );
    my ($node, $inode_local, $sep, $NCMD );

    $NCMD = scalar(@CMD);

    # determine name of script
    if ( $SSH_NNODE ) {
	$inode_local = $$INODE;
	$inode_local++ ;
	if ($inode_local >= $SSH_NNODE ) { $inode_local = 0 ; }
	$node = $SSH_NODELIST[$inode_local] ;
	$$INODE = $inode_local ;  # set return arg

	# get name of script for this node
	$script = "$SPLIT_JOBDIR/${PREFIX}_$node";
    }
    else {
	$script = "$SPLIT_JOBDIR/${PREFIX}_BATCH";
    }

    # open script for first time or in append mode.
    if ( -e $script )  { 
	# open in append mode.
	open PTR_SCRIPT , ">> $script" ;  
    }
    else  { 
	# open for first time
	open  PTR_SCRIPT , "> $script" ; 
	print PTR_SCRIPT   "$SNANA_LOGIN_SETUP \n";    # Sep 13 2014
	print PTR_SCRIPT   "cd $SPLIT_JOBDIR \n\n";
    }

    # --------------------------------------------
    
    if ( $BATCH_NCORE > 0 ) {
	$batchFile = "${PREFIX}_${VERSION}_${FFF}.BATCH" ;  
	$batchLog  = "${PREFIX}_${VERSION}_${FFF}.BATCH-LOG" ;
	$batchJob  = "" ; $icmd=0;
	foreach $cmd (@CMD) {
	    $sep = ";" ;
	    $icmd++ ; if ( $icmd == $NCMD ) { $sep = "" ; }
	    $batchJob = "$batchJob $cmd $sep " ;
	}
	print PTR_SCRIPT "$BATCH_COMMAND $batchFile \n";
	sntools::make_batchFile($BATCH_TEMPLATE,$SPLIT_JOBDIR, 
				$batchFile, $batchLog, $BATCH_MEM,$batchJob);
    }
    else {
	foreach $cmd (@CMD) {  print PTR_SCRIPT "$cmd \n"; }
    }

    close PTR_SCRIPT ;

    return ;

} # end create_batchScripts

# =========================================
sub get_PREFIX_SPLITJOB {
  my ($iver,$ifitopt,$isplit) = @_ ;

  my ($FITOPTnnn, $SPLITnnn, $prefix_tmp);

  $FITOPTnnn = sprintf( "FITOPT%3.3d", $ifitopt ) ;

  if ( $isplit >= 0 ) {
      $SPLITnnn   = sprintf( "SPLIT%3.3d" , $isplit  ) ;
      $prefix_tmp = "$VERSION_LIST[$iver]_${FITOPTnnn}_${SPLITnnn}" ;
  }
  else {
      $prefix_tmp = "$VERSION_LIST[$iver]_${FITOPTnnn}" ;
  }

  return $prefix_tmp ;

}  # end of get_PREFIX_SPLITJOB 

# ===============================
sub get_PREFIX_BATCHFILE {
    my ($ifile) = @_ ;
    my $prefix = sprintf("FITSPLIT_%5.5d", $ifile);
    return($prefix);
}

# ======================================
sub split_job_prep {

    my ($cmd, $script, $currentDir, $ftmp);

    # after creating all the FIT* scripts, do some misc. stuff 
    # in the SPLIT_JOB directory

    # give +x permissions
    qx(chmod +x ${SPLIT_JOBDIR_LCFIT}/FIT*);

    # leave file with expected number of
    # JOBS-PER-FITOPT and TOTAL
    # (to compare later with number of DONE files)
    $ftmp = "${SPLIT_JOBDIR_LCFIT}/${NJOBS_FILENAME}" ;
    open  PTR_TMP , "> $ftmp" ;
    print PTR_TMP "$NSPLIT $NFITJOB_TOT\n";
    close PTR_TMP ;


    # Mar 2019: write time stamp so that MERGE process can check later.
    $CDATE_SUBMIT  = `date +%Y%m%d_%H%M_%S` ;
    $CDATE_SUBMIT  =~ s/\s+$// ;   # trim trailing whitespace
    $ftmp = "${SPLIT_JOBDIR_LCFIT}/${TIME_STAMP_FILENAME}" ;
    open  PTR_TMP , "> $ftmp" ;
    print PTR_TMP "$CDATE_SUBMIT\n";
    close PTR_TMP ;

    
} # end of split_job_prep


# ================================
sub submitJobs_lcfit {

    my ($query, $isplit, $node, $cmd,  $NF, $N_usec );
    my ($script, $SCRIPT);

    $NF = $NFITJOB_SUBMIT ;

    if ( $SSH_NNODE ) {
	$query = "Ready to ssh $NF fit-jobs over $SSH_NNODE nodes";
    }
    elsif ( $BATCH_NCORE > 0 ) {
	$query = "Ready to launch $NF fit-jobs in $BATCH_COMMAND batch.";
    } 
    else {
	$query = "Ready to run $NF fit-jobs interactively.";
    }

# --------
    print "\n $query \n";


    # for ASCII format, split the version(s) just before submitting
    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
	&split_version($iver) ;
    }

    my $NSUBMIT = $NSPLIT ;
    if ( $BATCH_NCORE > 0 ) { $NSUBMIT = 1; }

    for  ( $isplit = 1; $isplit <= $NSUBMIT; $isplit++ ) {
	$node     = $NODEMAP[$isplit];
	$SCRIPT   = $SCRIPT_LIST[$isplit] ;
	$script   = $script_LIST[$isplit] ;

	if ( $SSH_NNODE ) {
	    $cmd = "ssh -x $node ${qq} source $SCRIPT ${qq}" ;
	}
	elsif ( $BATCH_NCORE > 0 ) {
	    $cmd = "cd $SPLIT_JOBDIR_LCFIT ; ./$script";
	}  
	else {
	    $cmd = "$SCRIPT" ;
	}
	print " $cmd \n" ; 
	system("$cmd &") ;

    }

}  # end of submitJobs_lcfit

# ========================================
sub submitJobs_postProc(@) {

    my($WHATJOB, $PREFIX, $SPLIT_JOBDIR) = @_ ;

    # Sumbit generic set of postProcessing (after LCFITs)
    # WHATJOBS is a comment string , e.g, "ML" or "SALT2mu"
    #

    my ($cmd, $cdd, $scriptName, $iscript, $node , $ls);
    my (@scriptList, $NSCRIPT, $LEN_PREFIX, $txt);

    $LEN_PREFIX = length($PREFIX) ;
    $cdd        = "cd $SPLIT_JOBDIR";

    if ( $BATCH_NCORE > 0 ) 
    { $ls  = "ls ${PREFIX}_BATCH"; } # only the master script for batch
    else 
    { $ls  = "ls ${PREFIX}*"; }      # all scripts for NODELIST

    @scriptList = qx($cdd; $ls ) ;
    $NSCRIPT    = scalar(@scriptList);

    $txt = "\n Submit $WHATJOB jobs ... \n"; 
    print            "$txt" ;
    print PTR_MERGE2 "$txt" ;

    for ( $iscript = 0; $iscript < $NSCRIPT; $iscript++ ) {
	$scriptName = $scriptList[$iscript];
	$scriptName =~ s/\s+$// ;   # trim trailing whitespace

	if ( $SSH_NNODE ) {
	    $node   = substr($scriptName,$LEN_PREFIX+1,50);
	    $node   =~ s/\s+$// ;   # trim trailing whitespace
	    $cmd    = "ssh -x $node ${qq} $cdd ; source $scriptName ${qq}" ;
	}
	elsif ( $BATCH_NCORE > 0 ) {
	    $cmd = "$cdd ; ./$scriptName";
	}
	else {
	    $cmd = "$scriptName" ;  # interactive
	}
	print " $cmd \n" ;
	print PTR_MERGE2 " $cmd \n" ;
	system("$cmd &") ;
    }

    PTR_MERGE2 -> autoflush(1);

    return ;

} # end submitJobs_postProc

# =============================================
sub DONE_CHECK_postProc {

    my ($WHATJOB, $PREFIX, $SPLIT_JOBDIR) = @_ ;

    my ($NEXPECT, $DONEFILE , @TMPLIST, @DONELIST, $NDONE );

    print PTR_MERGE2 "\n Check on $WHATJOB jobs ... \n";
    $DONEFILE   = "$OUTDIR/${WHATJOB}.$SUFFIX_DONE" ;
    $NEXPECT    = qx(cat $SPLIT_JOBDIR/$NJOBS_FILENAME);

    print PTR_MERGE2 
	"\t waiting for $NEXPECT  $WHATJOB jobs to finish ... \n" ;
    PTR_MERGE2 -> autoflush(1);

  CHECK_AGAIN:

    @DONELIST  = qx(ls $SPLIT_JOBDIR/*.DONE 2>/dev/null ) ;
    $NDONE     = scalar(@DONELIST);

    if ( $NDONE != $NEXPECT ) {
	sleep($WAIT_CHECK_DONE);
	$TOTAL_WAIT += $WAIT_CHECK_DONE ;
	&check_waitTime();
	goto CHECK_AGAIN ;
    }

    print PTR_MERGE2 "\n Finshed $WHATJOB jobs. \n";
    PTR_MERGE2 -> autoflush(1);

    # create local DONE-file for this process
    qx(touch $DONEFILE);

} # end  of DONE_CHECK_postProc


# =======================================
sub DONE_CHECK_LCFIT {

    my ($iver,$ifitopt) = @_ ;

    # Re-write Oct 2012:
    # For this iver and ifitopt, return integer flag
    #   $IFLAG_DONE     if merging already done
    #   $IFLAG_DOMERGE  if all DONE files exist  for FIT jobs
    #   $IFLAG_WAIT     if still processing or still merging
    #
    # July 28, 2013: fix returning IFLAG_DONE to check 
    #                all outFile types: HBOOK, ROOT ...
    #
    # Oct 15 2013: one-line fix to work with either HBOOK or ROOT.
    #               --> must wait for FITRES-merging.
    #
    # Sep 10 2016: check for FITOPT_LINK
    #
    # Nov 17 2017
    #  when fits are all done, sleep(1) only on ifitopt=0 -->
    #  Merge process is much faster.
    # 
    # Dec 13 2017: if NSPLIT=1, search for exact DONE-FILE match to avoid
    #              slow-down with wild card.
    #
    my ($NEXPECT,  @line, @wdlist, $PREFIX, $indx );
    my (@TMPLIST, $cGRACE, $cABORT, $VCOMB );
    my ($VERSION, $SUFFIX, $FF, $SUF, $OUTFILE, $ISDONE_LCFIT ) ;
    my ($DONEFILE, $MERGEDIR, $DONEMERGE, $f0, $f1, $searchString );
    my ($tnow, $tlast, $IFITOPT, $IFITOPT_LINK, $SIM_FLAG );

    # check if merged hbook and/or root file(s) already exist.

    $VERSION  = "$VERSION_LIST[$iver]" ;
    $MERGEDIR = "$OUTDIR_LIST[$iver]";
    $FF       = sprintf("FITOPT%3.3d", $ifitopt );
    $IFITOPT_LINK = $IFITOPT_LINK_LIST[$ifitopt];
    $SIM_FLAG     = $SIM_FLAG_LIST[$iver] ;

    # --------------------------------------------
    # check if all outputs have been merged.
    $DONEMERGE = &get_mergeStatus($iver,$ifitopt);

    if ( $DONEMERGE ) { 
	$MERGE_STATUS[$iver][$ifitopt] = 2 ; # flag that merging ended
	return $IFLAG_DONE ; 
    }

    # -----------------------------------------------------------
    $PREFIX   = &get_PREFIX_SPLITJOB($iver,$ifitopt,-1);

    if ( $NSPLIT == 1 ) 
    { $searchString = "${SPLIT_JOBDIR_LCFIT}/${PREFIX}_SPLIT001.DONE" ; }
    else 
    { $searchString = "${SPLIT_JOBDIR_LCFIT}/${PREFIX}*.DONE" ; }
    @DONELIST_LCFIT  = qx(ls $searchString 2>/dev/null );
    $NDONE_LCFIT     = scalar(@DONELIST_LCFIT);
    $ISDONE_LCFIT    = ( $NDONE_LCFIT == $NEXPECT_PER_MERGE ) ;

    # for LINK option, check if merging is done for IFITOPT_LINK 
    if ( $IFITOPT_LINK >= 0 ) 
    { $ISDONE_LCFIT = &get_mergeStatus($iver,$IFITOPT_LINK); }

    if ( $SIM_FLAG == $SIM_FLAG_SYMLINK2 ) 
    { $ISDONE_LCFIT = 1; }
 
    # ---------------- check status ---------------
    if ( $ISDONE_LCFIT == 0 ) {
	# fitting not done
	return $IFLAG_WAIT ;
    } 
    elsif ( $MERGE_STATUS[$iver][$ifitopt] > 0 ) {
	# fitting done, but still merging (Mar 2016)
	return $IFLAG_WAIT
    }
    else {
	# fits are all done, so start merging.
	# determine number of GRACEFUL finished vs. ABORTs

	if ( $ifitopt==0 ) { sleep(1); }  # Nov 17 2017
	
	if ( $NSPLIT == 1 ) 
	{ $searchString = "${SPLIT_JOBDIR_LCFIT}/${PREFIX}_SPLIT001.LOG" ; }
	else
	{ $searchString = "${SPLIT_JOBDIR_LCFIT}/${PREFIX}*.LOG" ; }

	@TMPLIST  = qx(grep GRACE $searchString );
	$NGRACE   = scalar(@TMPLIST);
	$NGRACE_TOTAL += $NGRACE ;
	
	@TMPLIST  = qx(grep $qq ABORT $qq $searchString);
	$NABORT   = scalar(@TMPLIST);
	$NABORT_TOTAL += $NABORT ;
	
	$cGRACE = sprintf("%3d", $NGRACE );
	$cABORT = sprintf("%3d", $NABORT );

	print PTR_MERGE2 
	    " Found required $NDONE_LCFIT ${PREFIX}*DONE files.\n";
	PTR_MERGE2 -> autoflush(1);
	print PTR_MERGE2 "\t $cGRACE $PREFIX jobs finished gracefully. \n";
	print PTR_MERGE2 "\t $cABORT $PREFIX jobs ABORTED. \n";
	print PTR_MERGE2 "\n" ;	   
	PTR_MERGE2 -> autoflush(1);

	$MERGE_STATUS[$iver][$ifitopt] = 1 ; # flag that merging started

	return $IFLAG_DOMERGE ;		
    }
    
} # end of DONE_CHECK_LCFIT

# =====================================================
sub get_mergeStatus {
    my ($iver,$ifitopt) = @_ ;

    # Created Sep 2016
    # return 1 if merging is done for all outputs (HBOOK,ROOT,TEXT).
    # return 0 otherwise.

    my ($indx, $allDone, $DONEFILE );

    $allDone = 1;  
    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) {
	if ( $OUTFLAG[$indx] == 0 ) { next ; }
	$DONEFILE = $MERGED_DONEFILE[$indx][$iver][$ifitopt] ;
	if ( !( -e $DONEFILE ) ) { $allDone = 0 ; }	
    }

    return $allDone ;

} # end get_mergeStatus

# =====================================================
sub MERGE_DRIVER {

    # Jan 2 2015: touch DONE_STAMP here instead of inside merge_driver
    #             so that STAMP is created after SALT2mu jobs
    #
    # Apr 17, 2019: 
    #   + echo $DONE_STATUS to $DONE_STAMP_FILE (SUCCESS or FAILURE)
    #
    # Apr 29 2019: create DONE_STAMP after gzipping SPLIT_JOBS, 
    #              instead of before (dumb mistake).

    my ($iver, $ifitopt, $IFLAG );

    &merge_INIT();

  CHECK_FOR_DONE:
    $ALLDONE = 1 ;

    my $tnow    = localtime time ;  
    print PTR_MERGE2 " $tnow  ->  Check DONE files. \n" ;
    PTR_MERGE2 -> autoflush(1);

    # merge each version and fit-option
    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {

	$NFITOPT_MERGED = 0 ;

	for ( $ifitopt=$ifitopt_start ; $ifitopt < $NFITOPT; $ifitopt++ ) {

	    if ( $MERGE_STATUS[$iver][$ifitopt] == 2 ) { next; } # already done

	    $IFLAG = &DONE_CHECK_LCFIT($iver,$ifitopt);
	    
	    if ( $IFLAG == $IFLAG_DOMERGE )   { 
		&psnid_FOMsummary($iver,$ifitopt) ; 
		my $indx ;
		for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) {
		    &merge_outFiles($indx,$iver,$ifitopt); 
		}
		&get_NSN_MERGED($iver,$ifitopt) ; 
		&mergeLog_update($iver,$ifitopt,1);

		$ALLDONE = 0 ;   # Mar 3 2016
	    }
	    elsif ( $IFLAG == $IFLAG_DONE )  {  
		# do nothing
	    }  
	    elsif ( $IFLAG == $IFLAG_WAIT ) { 
		$ALLDONE = 0 ; 
	    }

	    if ( $IFLAG != $IFLAG_WAIT ) { $NFITOPT_MERGED++ ; }
	}  # ifitopt
	
	# gzip files for version when all FITOPT are done
	if ( $NFITOPT_MERGED == ($NFITOPT-$ifitopt_start) ) { 
	    &run_afterBurner($iver);  # Aug 31 2015
	    &gzip_VERSION($iver) ; 
	}

    }  # iver

    # -----------------------------------------
    if ( $ALLDONE == 0 ) {
	&check_TIME_STAMP();  # make sure we still have expected time stamp
	sleep($WAIT_CHECK_DONE);
	$TOTAL_WAIT += $WAIT_CHECK_DONE ;
	&check_waitTime();
	goto CHECK_FOR_DONE ;
    }

    # all jobs are done
    for ( $ifitopt=$ifitopt_start ; $ifitopt < $NFITOPT; $ifitopt++ )  {  
	my $indx ;
	for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) 
	{ &combine_all_outFiles($indx,$ifitopt);  }
    }

    &merge_end();

    if ( $MLFLAG == $MLFLAG_APPLY ) {
	my $PREFIX       = $PREFIX_RUNML ;
	my $SPLIT_JOBDIR = $SPLIT_JOBDIR_ML ;
	&submitJobs_postProc("MLAPPLY", $PREFIX, $SPLIT_JOBDIR);
	&DONE_CHECK_postProc("MLAPPLY", $PREFIX, $SPLIT_JOBDIR);
    }


    if ( $SALT2mu_INFILE ne "" ) {
	my $PREFIX       = $PREFIX_RUNSALT2mu ;
	my $SPLIT_JOBDIR = $SPLIT_JOBDIR_SALT2mu ;
	&submitJobs_postProc("SALT2mu", $PREFIX, $SPLIT_JOBDIR);
	&DONE_CHECK_postProc("SALT2mu", $PREFIX, $SPLIT_JOBDIR);
    }

    close PTR_MERGE2 ; # Dec 2014

    # Sep 12 2014: remove modified nml (if Legacy variables were removed)
    if ( index($nmlFile,"MODIFIED") > 0  ) { qx(rm $nmlFile) ; }

    # Feb 9 2017: check option to leave trail-mark in sim-readme
    if ( $TRAILMARK_FLAG ) { &write_TRAILMARK(); }

    # Dec 26 2014: check option to make SPLIT_JOBS.tar.gz
    if ( $GZIP_FLAG == 4 ) { &gzip_SPLIT_JOBS(); }

    # create optional global DONE_STAMP file
    if ( $DONE_STAMP_FLAG ) {  
	qx(touch $DONE_STAMP);  
	qx(echo $DONE_STATUS >> $DONE_STAMP);
    }

    return ;

} # end of MERGE_DRIVER


# ========================================
sub merge_submit {

    my ($CMD_MERGE, $LOG, $argList );

    if ( $noMERGE_FLAG ) { 
	print " BEWARE: MERGE job NOT launched ... \n";
	return ; 
    }

    sleep(1) ;
    print " Launch MERGE job in background ... \n";
    $argList   = "MERGE  -DATE_SUBMIT $CDATE_SUBMIT"  ;
    $CMD_MERGE = "$SCRIPTNAME_FULL $nmlFile  $argList" ;
    if ( $TRAILMARK_FLAG ) { $CMD_MERGE = "$CMD_MERGE TRAILMARK" ; }

    $LOG = "$MERGE_LOGDIR/MERGE2.LOG" ;
    system("$CMD_MERGE > $LOG & " ) ;
    sleep(2);
} 


# =====================================
sub hrecover {

    # Apr 2012
    # recovery option in case hbook files are not properly merged
    # at the end of split_and_fit. Used only for hbook files.
    
    my ($key, @tmp, $JOBNAME, $PREFIX_FULL, $HDIR, $FPREFIX, $suffix ); 

    $FPREFIX      = sprintf("FITOPT%3.3d", $HRECOVER_FITOPT) ;
    $PREFIX_FULL  = "${HRECOVER_PREFIX}_${FPREFIX}" ;

    print "\n Run HBOOK-merge on ${PREFIX_FULL}* \n";

    $key   = "OUTDIR:" ;
    @tmp   = sntools::parse_array($key, 1, $OPT_ABORT, @CONTENTS_NMLFILE);
    $OUTDIR = $tmp[0];
    $HDIR   = "$OUTDIR/SPLIT_JOBS" ;
    print " HDIR = $HDIR \n";

    #if merged hbook exists, then bail.
    $suffix = $SUFFIX_OUTFILE[$OUTINDX_HBOOK] ;
    my $HFILE_MERGED = "$OUTDIR/$HRECOVER_PREFIX/$FPREFIX.${suffix}" ;
    if ( -e $HFILE_MERGED ) {
	print " Merged hbook file already exists: \n";
	print "\t $HFILE_MERGED \n";
	print "\t Bye bye. \n";
	die "\n";
    }


    my ( @HFILE_LIST, @hfile_list, $HFILE, $hfile, $ls, $n, $tmpFile ) ;
    my $prefix_recover = "recover" ;
    my $hfile_recover  = "${prefix_recover}.${suffix}" ;
    my $logFile        = "${prefix_recover}_${PREFIX_FULL}.LOG" ;

    $ls = "${PREFIX_FULL}*.${suffix}" ; 

    # if files are gzipped, unzip them
    @HFILE_LIST = qx(cd $HDIR;  ls ${ls}.gz 2>/dev/null ) ;
    if ( scalar(@HFILE_LIST) > 0 ) {
	qx(cd $HDIR;  gunzip ${PREFIX_FULL}*.${suffix}.gz ) ; 
    }
	
    # -----
    @HFILE_LIST = qx(cd $HDIR;  ls $ls );
    @hfile_list = ();
    $n = 0 ;

    # remove pre-existing files.
    $tmpFile = "$HDIR/$hfile_recover";
    if ( -e $tmpFile) { qx(rm $tmpFile) ; }

    $tmpFile = "$HDIR/$logFile" ;
    if ( -e $tmpFile) { qx(rm $tmpFile) ; }

    foreach $HFILE ( @HFILE_LIST ) {
	$HFILE   =~ s/\s+$// ;   # trim trailing whitespace

	my $gz = index($HFILE,".gz");
	if ( $gz > 0 ) { qx(gunzip $HFILE) ; }

	# make temp symbolic link to lower case name
	# because stupid hmerge does not work on upper case names
	$n++ ;	$hfile = "${prefix_recover}_$n.${suffix}" ;

	print "\t $HFILE -> $hfile \n" ;
	qx(cd $HDIR ; ln -s $HFILE $hfile) ;	
	@hfile_list = ( @hfile_list , $hfile );
    }

    print "\t do hmerge ... \n";

    $JOBNAME = $JOBNAME_MERGE[$OUTINDX_HBOOK] ;
    qx(cd $HDIR ; $JOBNAME @hfile_list  $hfile_recover >& $logFile );
    qx(cd $HDIR ; rm @hfile_list);

    # finally, move recover.his into final location with
    # proper name
    my $mv  = "mv $hfile_recover $HFILE_MERGED";
#    print "mv: $mv \n";
    qx(cd $HDIR ; $mv);

    die "\n Done. \n";
    
} # end of hrecover


# =====================================
sub merge_INIT {

    my (@wdlist, @line, $iver, $ifitopt);

    # create separate subdir for merge-log files

    # setup MERGE.LOG file
    $MERGELOG   = "$OUTDIR/MERGE.LOG" ;
    $MERGELOG2  = "$MERGE_LOGDIR/MERGE2.LOG" ;

    qx(touch $MERGELOG); # create blank MERGE.LOG file

    open  PTR_MERGE2 , "> $MERGELOG2" ;

    print PTR_MERGE2 "Check DONE files in \n $SPLIT_JOBDIR_LCFIT \n\n";
    PTR_MERGE2 -> autoflush(1);

    if ( $PSNID_FLAG ) {
	$PSNID_SUMLOG = "$OUTDIR/PSNID_SUMMARY.LOG" ;
	qx(touch $PSNID_SUMLOG);
    }

    @line       = qx(cat ${SPLIT_JOBDIR_LCFIT}/${NJOBS_FILENAME});
    @wdlist     = split(/\s+/,$line[0]) ;
    $NEXPECT_PER_MERGE    = $wdlist[0] ; 
    $NEXPECT_TOTAL        = $wdlist[1] ; 

    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) { 	
	$GZIP_STATUS[$iver] = 0 ; 

	$OUTDIR_LIST[$iver] = "$OUTDIR/$VERSION_LIST[$iver]" ; # Sep 2016

	for ( $ifitopt=$ifitopt_start ; $ifitopt < $NFITOPT; $ifitopt++ ) 
	{ $MERGE_STATUS[$iver][$ifitopt] = 0 ; }
    }


    # ------------------------------------
    # Sep 10 2013
    # Create MERGE log for each version/fitopt with "notDone" status

    my $doneFlag = 0 ;
    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
	for ( $ifitopt=$ifitopt_start ; $ifitopt < $NFITOPT; $ifitopt++ ) {
	    &mergeLog_update($iver,$ifitopt,0);
	}
    }

} # end of merge_INIT


# =======================================
sub merge_outFiles {
 
    # merge hbook/root files based in $indx argument
    # Sep 14 2014: touch $mergeDone file to avoid conflicts
    #              when merged file exists but is not done merging.
    #
    # Dec 28 2014: if merge fails, try again.
    # Sep 10 2016: check option to make symbolic link to another FITOPT
    #                (see FITOPT_LINK_LIST)
    #
    my ($indx,$iver,$ifitopt) = @_;

    if ( $OUTFLAG[$indx] == 0 ) { return ; }

    my ($VERSION, $DATADIR, $MERGEDIR, $prefix, $cdMDIR, $JOBNAME);
    my ($suf, $outFile, $OUTFILE, $tmpLog, $outList, $logList );
    my ($cmd_h2root, $cmd_mkplots, $cmd_merge, $cmd_rmlocf, $cmd_done);
    my ($mergeLog, $mergeDone, $FLINK, $SIM_FLAG );
    my $LDMP = 0 ;

    $VERSION  = "$VERSION_LIST[$iver]" ;
    $DATADIR  = "$DATADIR_LIST[$iver]" ;
    $MERGEDIR = "$OUTDIR_LIST[$iver]";
    $cdMDIR   = "cd $MERGEDIR";    
    $SIM_FLAG = "$SIM_FLAG_LIST[$iver]" ;

    $PREFIX_SPLITJOB  = &get_PREFIX_SPLITJOB($iver,$ifitopt,-1);
    $PREFIX_MERGED    = sprintf("FITOPT%3.3d", $ifitopt );

    if ( $LDMP ) 
    { printf "xxx ------- MERGE $PREFIX_SPLITJOB -------- \n" ;  $| = 1; }

    # construct name of merged outfile
    $suf = &get_suffixMerged($indx,$ITABLE_FITRES);
    $outFile      = "${PREFIX_MERGED}.${suf}" ;
    $OUTFILE      = "${MERGEDIR}/${outFile}" ;
    $mergeLog     = "$MERGE_LOGDIR/MERGE_${VERSION}_${outFile}.LOG" ; 
    $mergeDone    = "$MERGE_LOGDIR/MERGE_${VERSION}_${outFile}.DONE" ;

    $MERGED_OUTFILE[$indx][$iver][$ifitopt]  = $OUTFILE ;
    $MERGED_outFile[$indx][$iver][$ifitopt]  = $outFile ;
    $MERGED_DONEFILE[$indx][$iver][$ifitopt] = $mergeDone ;

    # ----------------------------
    # check for making symLink instead of merging jobs
    $FLINK  = $FITOPT_LINK_LIST[$ifitopt] ;
    if ( $FLINK ne "" ) {
	my ($SUFFIX, $f0, $f1, $FF );
	$SUFFIX = "$SUFFIX_OUTFILE[$indx]"  ;
	$FF     = sprintf("FITOPT%3.3d", $ifitopt ) ;
	if ( $indx == $OUTINDX_TEXT ) { $SUFFIX = $SUFFIX_FITRES; }
	$f0 = "${FF}.$SUFFIX" ;
	$f1 = "${FLINK}.$SUFFIX" ;
	if ( -e "$MERGEDIR/$f1" ) {
	    qx(cd $MERGEDIR ; ln -s $f1 $f0); 
	    qx(cd $SPLIT_JOBDIR_LCFIT ; touch $mergeDone ) ;
	    $MERGE_STATUS[$iver][$ifitopt] = 2 ; 
	}
	return ;
    }

    # check if this job points symbolic links to another version.
    if ( -l $MERGEDIR ) {
	qx(cd $SPLIT_JOBDIR_LCFIT ; touch $mergeDone ) ;
	$MERGE_STATUS[$iver][$ifitopt] = 2 ; 
	return ; 
    }

    # ---------------------------------------------------
    # FITRES-text output is treated internally instead
    # of using an external merging program.
    if ( $indx == $OUTINDX_TEXT ) {
	&merge_outFiles_text($iver,$ifitopt); 
	return ;
    }

    # --------------------------------------------
    # if we get here than an external merging program is
    # used in which the argument is the list of output files
    # from each split job.
    
    $outList       = "${PREFIX_SPLITJOB}_SPLIT*.${suf}";

    $JOBNAME       = "$JOBNAME_MERGE[$indx]" ;
    $cmd_merge     = "$JOBNAME  $outList  $OUTFILE >& $mergeLog "; 
    $cmd_rmlocf    = "remove_locf_messages.pl $mergeLog QUIET" ;
    $cmd_done      = "touch $mergeLog" ;

    # Aug 29 2016: if just 1 split job, cp is much faster
    if ( $NSPLIT == 1 ) 
    { $cmd_merge = "cp  $outList $OUTFILE >& $mergeLog ";  }

    if ( $LDMP ) { printf "xxx cmd_merge = $cmd_merge \n" ;  $| = 1; }

    # do the hbook/root merge
    qx(cd $SPLIT_JOBDIR_LCFIT ; $cmd_merge ; $cmd_rmlocf ; touch $mergeDone ) ;

    # --------------------------------------------------------
    # Dec 28 2014
    # if OUTFILE does not exist, give it one more try ...
    # mainly for hbook that sometimes fails on 64 bit machines.
    # Make sure to append mergeLog to leave a record.

    unless ( -e $OUTFILE ) {
	$MSGERR[0] = "";
	$MSGERR[1] = "################################################### ";
	$MSGERR[2] = "  ERROR: $JOBNAME failed so try again ... " ;
	$MSGERR[3] = "################################################### ";
	$MSGERR[4] = "";
	my $msg;
	foreach $msg ( @MSGERR)	{ qx(echo "$msg\n" >> $mergeLog); }

	$cmd_merge  = "$JOBNAME  $outList  $OUTFILE >> $mergeLog "; 
	qx(cd $SPLIT_JOBDIR_LCFIT ; $cmd_merge ; $cmd_rmlocf ; touch $mergeDone ) ;
    }

    # check option to convert hbook file into root file with h2root.
    if ( $H2ROOT_FLAG > 0 && $indx == $OUTINDX_HBOOK ) {
	$tmpLog       = "root_${PREFIX_MERGED}.log" ;
	$cmd_h2root   = "h2root ${outFile}" ;
	qx($cdMDIR ; $cmd_h2root >& $tmpLog );
	if ( $H2ROOT_FLAG == 1 ) { qx($cdMDIR ; rm $tmpLog); } 
    }

    # check option to create light curve plots with mkfitplots.pl
    if ( $OPT_LCPLOT  && $OUTFLAG[$OUTINDX_HBOOK] ) {
	$cmd_mkplots  = "$JOBNAME_MKPLOTS --h $outFile $ARG_LCPLOT" ;
	$tmpLog       = "mkfitplots_${PREFIX_MERGED}.log" ;
	qx($cdMDIR ; $cmd_mkplots  >& $tmpLog );
    }

    return ;

} # end of merge_outFiles


# =============================
sub merge_outFiles_text { 
    my ($iver,$ifitopt) = @_ ;

    # Created Sep 12, 2014 : control text-merging

    if ( $DOTEXT_TABLE[$ITABLE_FITRES] ) {	
	# merge/catenate standard FITRES table.
	&merge_text_TABLE($ITABLE_FITRES, $iver, $ifitopt); 
# xxx mark delete	&merge_text_FITRES($iver,$ifitopt); 
	
	# pull variables out of root/hbook and append to text-fitres file
	&append_text_TABLE($ITABLE_FITRES,$iver,$ifitopt);
# xxx	&append_text_FITRES($iver,$ifitopt); # optional
    }

    if ( $DOTEXT_TABLE[$ITABLE_LCPLOT] ) {
	# merge text-LC files
	&merge_text_LCPLOT($iver,$ifitopt); 
    }

# Mar 2019: check for SNANA table too
    if ( $DOTEXT_TABLE[$ITABLE_SNANA] ) {   
	# merge/catenate SNANA table with all events
	&merge_text_TABLE($ITABLE_SNANA, $iver,$ifitopt); 
	&append_text_TABLE($ITABLE_SNANA,$iver,$ifitopt);
    }

} # end of merge_outFiles_text

# =======================================
sub merge_text_LCPLOT {
    ($iver,$ifitopt) = @_ ;

    # Jul 23 2013: 
    # catenate lightcurve text-dumps
    # into a single text file that is analogous to the 
    # .his or .root file.
    # In short, each .FITRES file corresponds to the 
    # .lightcurve file created here.


    my ($VERSION, $MERGEDIR, $cdMDIR, $cdSPLIT);
    my ($PREFIX_SPLITJOB, $PREFIX_MERGED, $catList, $catOut);

    if ( $DOTEXT_TABLE[$ITABLE_LCPLOT] == 0 ) { return ; }

    $VERSION  = "$VERSION_LIST[$iver]" ;
# xxx mark delete    $MERGEDIR = "$OUTDIR/$VERSION";
    $MERGEDIR = "$OUTDIR_LIST[$iver]";
    $cdMDIR   = "cd $MERGEDIR";
    $cdSPLIT  = "cd $SPLIT_JOBDIR_LCFIT" ;

    $PREFIX_SPLITJOB  = &get_PREFIX_SPLITJOB($iver,$ifitopt,-1);
    $PREFIX_MERGED    = sprintf("FITOPT%3.3d", $ifitopt );

    # define list of files to catendate/merge
    $catList = "${PREFIX_SPLITJOB}*.${SUFFIX_LCPLOT_TEXT}" ; 

    # catOut is the final catenated (merged) lightcurve file
    $catOut  = "${MERGEDIR}/${PREFIX_MERGED}.${SUFFIX_LCPLOT}" ;

    print PTR_MERGE2  " cat $catList \n" ;
    print PTR_MERGE2  " --> $catOut \n";
    PTR_MERGE2 -> autoflush(1);

#    print " xxx -------------------------- \n" ;
#    print " xxx LCPLOT catList = $catList \n";
#    print " xxx LCPLOT catOut  = $catOut \n";

    qx($cdSPLIT ; cat $catList > $catOut);


} # end of merge_text_LCPLOT


# =======================================
sub merge_text_TABLE {
    my ($ITABLE,$iver,$ifitopt) = @_;   
   
    # Created Mar 25 2019
    # Merge & atenate SNANA & FITRES text-files.
    # ITABLE is integer table type (SNANA or FITRES)
    #
    # Similar to old merge_text_FITRES, but handles
    # both the SNANA and FITRES tables.

    my ($VERSION, $DATADIR, $MERGEDIR, $cdSPLIT, $cdMDIR);
    my ($suf, $sufCat, $tmpFile1, $dumpFile, @tmp, $NSN );
    my ($cmd1_fres, $cmd2_fres, $searchString );
    my ($sedcmd, $cmd0, $cmd1, $cmd2, $cmd3, $cmd4, $cmd5 );
    my ($f1, $f2, $FRES_FILE, $FITOPT, $VCOMB );
    my ($line_NVAR, $line_VARNAMES, $sedcmd );
    my (@comment_orig, @comment_add, $comment, @tmpList);
    my ($NFILES, $FOUND_TABLEFILES );
    my ($TEXTFILE, $textFile, $TEMPFILE) ;
    my $PND = '#' ;

    # ---------- BEGIN --------------

    $VERSION  = "$VERSION_LIST[$iver]" ;
    $DATADIR  = "$DATADIR_LIST[$iver]" ;
    $MERGEDIR = "$OUTDIR_LIST[$iver]";
    $cdMDIR   = "cd $MERGEDIR";
    $cdSPLIT  = "cd $SPLIT_JOBDIR_LCFIT" ;

    $PREFIX_SPLITJOB  = &get_PREFIX_SPLITJOB($iver,$ifitopt,-1);
    $PREFIX_MERGED    = sprintf("FITOPT%3.3d", $ifitopt );

    print PTR_MERGE2 
	"\n Merge $TABLELIST[$ITABLE] TEXT-table " . 
	"for $VERSION($PREFIX_MERGED) \n";
    PTR_MERGE2 -> autoflush(1);	


    # construct commands to catentate fitres files
    if ( $ITABLE == $ITABLE_FITRES ) {
	# FITRES table
	$suf       = "${SUFFIX_FITRES_TEXT}" ;  # suffix on split-job files
	$sufCat    = "${SUFFIX_FITRES}" ;       # suffix on cat file
	$tmpFile1  = "${PREFIX_SPLITJOB}.${sufCat}" ;
	$textFile  = "$MERGED_outFile[$OUTINDX_TEXT][$iver][$ifitopt]" ;
	$TEXTFILE  = "$MERGED_OUTFILE[$OUTINDX_TEXT][$iver][$ifitopt]" ;
    }    
    else {
	# SNANA table
	$suf       = "${SUFFIX_SNANA_TEXT}" ;  # suffix on split-job files
	$sufCat    = "${SUFFIX_SNANA}" ;       # suffix on cat file
	$tmpFile1  = "${PREFIX_SPLITJOB}.${sufCat}" ;
	$textFile  = "${PREFIX_MERGED}.${sufCat}" ;
	$TEXTFILE  = "$MERGEDIR/$textFile";  
    }


    # check if some fitres files exist
    if ( $ifitopt == 0 ) {
	@tmpList  = qx($cdSPLIT ; ls ${PREFIX_SPLITJOB}*${suf} 2>/dev/null );
	$NFILES   = scalar(@tmpList);
	$FOUND_TABLEFILES = ( $NFILES > 0 ) ;
    } 
    else { 
	$FOUND_TABLEFILES = 1; 
    }


    if ( $FOUND_TABLEFILES == 0 ) {
	$cmd1_fres = "echo WARNING: No $TABLELIST[$ITABLE] " .
	    "table-files found for $PREFIX_SPLITJOB" ;
	$cmd2_fres =  "touch $TEXTFILE" ;
    }
    else {
	if ( $NSPLIT == 1 ) 
	{ $searchString = "${PREFIX_SPLITJOB}_SPLIT001.${suf}" ; }
	else
	{ $searchString = "${PREFIX_SPLITJOB}_SPLIT*${suf}" ; }

	$cmd1_fres = "cat ${searchString} > $tmpFile1 ";
	$cmd2_fres = "mv $tmpFile1 $TEXTFILE" ;
    }

    # check for SIMGEN_DUMP file.
    if ( $ITABLE == $ITABLE_FITRES ) {
	if ( ($SIM_FLAG_LIST[$iver] == $SIM_FLAG_NORMAL) &&  $ifitopt==0 ) 
	{ &copy_SIMGEN_DUMP($iver) ; }
    }

    # do the catenating
    qx($cdSPLIT ; $cmd1_fres ; $cmd2_fres ) ;	    
    
    # Cleanup FITRES file so that NVAR and VARNAMES appear only once.
    # Also add comment at top

    # For FITRES-file comments, replace single quote(s) with 
    # double quote(s) to avoid sed problems.
    $FITOPT = "$FITOPT_LIST[$ifitopt]";
    $FITOPT =~ s/\'/$qq/g ;

    $comment_add[0] = "Catenated from $NSPLIT split-jobs.";
    $comment_add[1] = "SNANA VERSION: $CURRENT_SNANA_VERSION ";   
    $comment_add[2] = "DATA  VERSION: $VERSION " ;
    $comment_add[3] = "FITOPT:  $FITOPT";
    $comment_add[4] = "---------------------------------------- ";

    @tmp           = qx(grep "NVAR: "     $TEXTFILE);
    $line_NVAR     = "$tmp[0]";   
    @tmp           = qx(grep "VARNAMES: " $TEXTFILE);
    $line_VARNAMES = "$tmp[0]" ;

    # get optional job-comments from first split-file;
    # we assume that the comments in each split-job are the same.
    my $firstFile = "$SPLIT_JOBDIR_LCFIT/${PREFIX_SPLITJOB}_SPLIT001.${suf}";
    @comment_orig  = qx(grep '$PND' $firstFile 2>/dev/null ) ;


    $sedcmd = "sed ";
    foreach $comment (@comment_add) { 
	$sedcmd = "$sedcmd " . "-e '1i $PND $comment'" ; 
    } 

    # add in original comments (if they exist)
    foreach $comment (@comment_orig)  { 
	chomp $comment ;
	$sedcmd = "$sedcmd " . "-e '1i $comment'" ; 
    } 

    if ( length($line_NVAR) > 0 ) 
    { $sedcmd = "$sedcmd " . "-e '1i $line_NVAR'" ; }

    if ( length($line_VARNAMES) > 0 ) 
    { $sedcmd = "$sedcmd " . "-e '1i $line_VARNAMES'" ;  }


    $sedcmd = "$sedcmd " . "-e '/NVAR/d'" ;
    $sedcmd = "$sedcmd " . "-e '/VARNAMES/d'" ;
    $sedcmd = "$sedcmd " . "-e '/#/d'"  ;  # remove comments
    $sedcmd = "$sedcmd " . "-e '/^\$/d'" ;  # remove blank lines (Dec 2018) 

    $TEMPFILE = "TEMPXXX_${textFile}" ;
    qx($sedcmd $TEXTFILE > $TEMPFILE ; mv $TEMPFILE $TEXTFILE );


    # -----------------------------------------
    # check option to combine all the FITRES files with a
    # fitres file from a particular version
    # WARNING (5/1/2013: this part needs testing
    $VCOMB = "$VERSION_FITRES_COMBINE" ;
    $FRES_FILE = "${PREFIX_MERGED}.${SUFFIX_FITRES}" ;
    if ( $VCOMB  ne  ''  &&  $ITABLE == $ITABLE_FITRES ) {
	# print PTR_MERGE " Combine each merged FITRES with $VCOMB \n" ;
	$f1   = "$FRES_FILE" ;
	$f2   = "../${VCOMB}/${FRES_FILE}";
	$cmd1 = "$JOBNAME_COMBINE  $f1  $f2" ;
	$cmd2 = "mv $f1 ${f1}_ORIG" ;
	$cmd3 = "mv combine_fitres.text $f1" ;
	$cmd4 = "rm combine_fitres.* " ;
	qx($cdMDIR ; $cmd1 ; $cmd2 ; $cmd3 ; $cmd4  );
    }
	

    # Sep 14 2014: touch done file
    
    if ( $ITABLE == $ITABLE_FITRES ) {
	my $DONEFILE = $MERGED_DONEFILE[$OUTINDX_TEXT][$iver][$ifitopt] ;
	qx(touch $DONEFILE);
    }

    return ;

} # end merge_text_TABLE

# =======================================
sub merge_text_FITRES {

    # !!!!!!!! Mar 25 2019: slated to be OBSOLETE !!!!!!!!
    #
    # catenate fitres text-files.
    # July 2016: require ifitopt==0 to call copy_SIMGEN_DUMP().
    # Dec 13 2016: if NSPLIT==1, use exact string match instead of wildcard.

    my ($iver,$ifitopt) = @_;

    my ($VERSION, $DATADIR, $MERGEDIR, $cdSPLIT, $cdMDIR);
    my ($suf, $sufCat, $tmpFile1, $dumpFile, @tmp, $NSN );
    my ($cmd1_fres, $cmd2_fres, $searchString );
    my ($sedcmd, $cmd0, $cmd1, $cmd2, $cmd3, $cmd4, $cmd5 );
    my ($f1, $f2, $FRES_FILE, $FITOPT, $VCOMB );
    my ($line_NVAR, $line_VARNAMES, $sedcmd );
    my (@comment_orig, @comment_add, $comment, @tmpList);
    my ($NFITRES, $FOUND_FITRES );
    my ($TEXTFILE, $textFile, $TEMPFILE) ;
    my $PND = '#' ;

    # ---------- BEGIN --------------

    $VERSION  = "$VERSION_LIST[$iver]" ;
    $DATADIR  = "$DATADIR_LIST[$iver]" ;
    $MERGEDIR = "$OUTDIR_LIST[$iver]";
    $cdMDIR   = "cd $MERGEDIR";
    $cdSPLIT  = "cd $SPLIT_JOBDIR_LCFIT" ;

    $PREFIX_SPLITJOB  = &get_PREFIX_SPLITJOB($iver,$ifitopt,-1);
    $PREFIX_MERGED    = sprintf("FITOPT%3.3d", $ifitopt );

    # construct commands to catentate fitres files
    $suf       = "${SUFFIX_FITRES_TEXT}" ;  # suffix on split-job files
    $sufCat    = "${SUFFIX_FITRES}" ;       # suffix on cat file
    $tmpFile1  = "${PREFIX_SPLITJOB}.${sufCat}" ;
    $textFile  = "$MERGED_outFile[$OUTINDX_TEXT][$iver][$ifitopt]" ;
    $TEXTFILE  = "$MERGED_OUTFILE[$OUTINDX_TEXT][$iver][$ifitopt]" ;

    # check if some fitres files exist
    if ( $ifitopt == 0 ) {
	@tmpList  = qx($cdSPLIT ; ls ${PREFIX_SPLITJOB}*${suf} 2>/dev/null );
	$NFITRES  = scalar(@tmpList);
	$FOUND_FITRES = ( $NFITRES > 0 ) ;
    }
    else 
    { $FOUND_FITRES = 1; }

    if ( $FOUND_FITRES == 0 ) {
	$cmd1_fres = 
	    "echo WARNING: No fitres files found for $PREFIX_SPLITJOB" ;
	$cmd2_fres =  "touch $TEXTFILE" ;
    }
    else {

	if ( $NSPLIT == 1 ) 
	{ $searchString = "${PREFIX_SPLITJOB}_SPLIT001.${suf}" ; }
	else
	{ $searchString = "${PREFIX_SPLITJOB}_SPLIT*${suf}" ; }

	$cmd1_fres = "cat ${searchString} > $tmpFile1 ";
	$cmd2_fres = "mv $tmpFile1 $TEXTFILE" ;
    }

    # check for SIMGEN_DUMP file.
    if ( ($SIM_FLAG_LIST[$iver] == $SIM_FLAG_NORMAL) &&  $ifitopt==0 ) 
    { &copy_SIMGEN_DUMP($iver) ; }

    # do the catenating
    qx($cdSPLIT ; $cmd1_fres ; $cmd2_fres ) ;	    
    
    # Cleanup FITRES file so that NVAR and VARNAMES appear only once.
    # Also add comment at top

    # For FITRES-file comments, replace single quote(s) with 
    # double quote(s) to avoid sed problems.
    $FITOPT = "$FITOPT_LIST[$ifitopt]";
    $FITOPT =~ s/\'/$qq/g ;

    $comment_add[0] = "Catenated from $NSPLIT split-jobs.";
    $comment_add[1] = "SNANA VERSION: $CURRENT_SNANA_VERSION ";   
    $comment_add[2] = "DATA  VERSION: $VERSION " ;
    $comment_add[3] = "FITOPT:  $FITOPT";
    $comment_add[4] = "---------------------------------------- ";

    @tmp           = qx(grep "NVAR: "     $TEXTFILE);
    $line_NVAR     = "$tmp[0]";   
    @tmp           = qx(grep "VARNAMES: " $TEXTFILE);
    $line_VARNAMES = "$tmp[0]" ;

    # get optional job-comments from first split-file;
    # we assume that the comments in each split-job are the same.
    my $firstFile = "$SPLIT_JOBDIR_LCFIT/${PREFIX_SPLITJOB}_SPLIT001.${suf}";
    @comment_orig  = qx(grep '$PND' $firstFile 2>/dev/null ) ;


    $sedcmd = "sed ";
    foreach $comment (@comment_add) { 
	$sedcmd = "$sedcmd " . "-e '1i $PND $comment'" ; 
    } 

    # add in original comments (if they exist)
    foreach $comment (@comment_orig)  { 
	chomp $comment ;
	$sedcmd = "$sedcmd " . "-e '1i $comment'" ; 
    } 

    if ( length($line_NVAR) > 0 ) 
    { $sedcmd = "$sedcmd " . "-e '1i $line_NVAR'" ; }

    if ( length($line_VARNAMES) > 0 ) 
    { $sedcmd = "$sedcmd " . "-e '1i $line_VARNAMES'" ;  }


    $sedcmd = "$sedcmd " . "-e '/NVAR/d'" ;
    $sedcmd = "$sedcmd " . "-e '/VARNAMES/d'" ;
    $sedcmd = "$sedcmd " . "-e '/#/d'"  ;  # remove comments
    $sedcmd = "$sedcmd " . "-e '/^\$/d'" ;  # remove blank lines (Dec 2018) 

    $TEMPFILE = "TEMPXXX_${textFile}" ;
    qx($sedcmd $TEXTFILE > $TEMPFILE ; mv $TEMPFILE $TEXTFILE );


    # -----------------------------------------
    # check option to combine all the FITRES files with a
    # fitres file from a particular version
    # WARNING (5/1/2013: this part needs testing
    $VCOMB = "$VERSION_FITRES_COMBINE" ;
    $FRES_FILE = "${PREFIX_MERGED}.${SUFFIX_FITRES}" ;
    if ( $VCOMB  ne  '' ) {
	# print PTR_MERGE " Combine each merged FITRES with $VCOMB \n" ;
	$f1   = "$FRES_FILE" ;
	$f2   = "../${VCOMB}/${FRES_FILE}";
	$cmd1 = "$JOBNAME_COMBINE  $f1  $f2" ;
	$cmd2 = "mv $f1 ${f1}_ORIG" ;
	$cmd3 = "mv combine_fitres.text $f1" ;
	$cmd4 = "rm combine_fitres.* " ;
	qx($cdMDIR ; $cmd1 ; $cmd2 ; $cmd3 ; $cmd4  );
    }
	

    # Sep 14 2014: touch done file
  CREATE_DONEFILE:
    my $DONEFILE = $MERGED_DONEFILE[$OUTINDX_TEXT][$iver][$ifitopt] ;
    qx(touch $DONEFILE);

    return ;

} # end of merge_text_FITRES

# ===========================
sub append_text_TABLE {

    # Mar 25 2019
    # 
    # Check option to append text-table file with variables
    # from the hbook/root file. Uses sntable_dump.pl.
    #    

    my ($ITABLE,$iver,$ifitopt) = @_ ;

    my ($indx, @wdlist, $VERSION, $MERGEDIR );
    my ($CMD, $outFile_fit, $OUTFILE_FIT, $dumpFile, $logFile );
    my ($textFile, $TABLENAME, $VARLIST, $ARG_V, $ARG_A, $ARG_O );
    my ($prefix_append, $outFile_append, $SUFFIX);

    # --------------- BEGIN ----------

    # bail if option is NOT set.   
    if ( $APPEND_TABLE_LIST  eq "" ) { return ; }

    $TABLENAME = $TABLELIST[$ITABLE];
    @wdlist   = split(/\s+/,$APPEND_TABLE_LIST) ;
    $textFile = "$MERGED_outFile[$OUTINDX_TEXT][$iver][$ifitopt]" ;
    $VERSION  = "$VERSION_LIST[$iver]" ;
    $MERGEDIR = "$OUTDIR_LIST[$iver]" ;

    if ( $ITABLE == $ITABLE_FITRES ) {
	$textFile = "$MERGED_outFile[$OUTINDX_TEXT][$iver][$ifitopt]" ;
    }
    else {
	# SNANA table
	$textFile  = "${PREFIX_MERGED}.${SUFFIX_SNANA}" ;
    }

    # internal suffix for outfiles produced by sntable_dump
    $SUFFIX = "${PREFIX_MERGED}" ;

    if ( $OUTFLAG[$OUTINDX_HBOOK] ) 
    { $indx =  $OUTINDX_HBOOK ; } 
    elsif ( $OUTFLAG[$OUTINDX_ROOT] ) 
    { $indx =  $OUTINDX_ROOT ; }
    else {
	$MSGERR[0] = "Must define HBOOK and/or ROOT files";
	$MSGERR[1] = "to append TEXT table file.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    # outFile_fit refers to hbook/root output of fit job -> 
    # input to append script
    $outFile_fit    = "$MERGED_outFile[$indx][$iver][$ifitopt]" ;  
    $OUTFILE_FIT    = "$MERGED_OUTFILE[$indx][$iver][$ifitopt]" ;  

    # dumpFile is the intermediate text file containing only the
    # appended variables extracted from outFile_fit

    $dumpFile    = "sntable_dump_${SUFFIX}_${TABLENAME}.text" ; 
    $logFile     = "sntable_dump_${SUFFIX}_${TABLENAME}.log" ;

    $VARLIST  = "@wdlist[0 .. $#wdlist]" ;
    $ARG_V    = "--v $VARLIST" ;
    $ARG_A    = "--a  $textFile" ;
    $ARG_O    = "--o  $dumpFile" ;
    $CMD = "$JOBNAME_SNTABLE $outFile_fit $TABLENAME  $ARG_V  $ARG_A  $ARG_O" ;

    print PTR_MERGE2 "\n APPEND_TABLE_TEXT Command for $VERSION/$SUFFIX: \n";
    print PTR_MERGE2 "$CMD\n";
    PTR_MERGE2 -> autoflush(1);	

    unless ( -e $OUTFILE_FIT ) {
	print PTR_MERGE2 "\t Cannot find '$outFile_fit' for $VERSION ";
	print PTR_MERGE2 " --> skip APPEND_TABLE_TEXT \n";
	PTR_MERGE2 -> autoflush(1);	
	return ;
    }

    qx(cd $MERGEDIR ; $CMD > $logFile);

    # check for abort (June 26 2015)
    my @bla = qx(cd $MERGEDIR ; grep ABORT $logFile);
    if ( scalar(@bla) > 0 ) {
	$MSGERR[0] = "FATAL ERROR: Append-ABORT found in ";
	$MSGERR[1] = "   $MERGEDIR/$logFile";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    } 


    # set outFile_append = name of appended fitres file
    # (or grep 'Write' from logfile to determine name)
    $prefix_append  = "sntable_append_${SUFFIX}" ;
    $outFile_append = "${prefix_append}.text" ;

    # move appended text file to replace original, 
    # and save original with _ORIG suffix.
    qx(cd $MERGEDIR ; mv $textFile        ${textFile}_ORIG  );
    qx(cd $MERGEDIR ; mv $outFile_append  ${textFile}       );

    # cleanup
    qx(cd $MERGEDIR ; rm ${prefix_append}.*  $logFile $dumpFile ) ;
    qx(cd $MERGEDIR ; gzip ${textFile}_ORIG );

    return ;

} # end of append_text_TABLE


# ===========================
sub append_text_FITRES {

    # !!!! Mar 25 2019:LEGACY function soon to be OBSOLETE !!!!
    #
    # Check option to append fitrest-text file with variables
    # from the hbook/root file. Uses sntable_dump.pl.
    #
    # May 9, 2013: give unique sntable file-names to avoid clobber.
    # Feb 8, 2014: if outFile_fit does not exist, leave message in PTR_MERGE2
    #              and return.
    #
    # Sep 13, 2014: when done, remove sntable* files and gzip _ORIG file.
    #

    my ($iver,$ifitopt) = @_ ;

    my ($indx, @wdlist, $VERSION, $MERGEDIR );
    my ($CMD, $outFile_fit, $OUTFILE_FIT, $dumpFile, $logFile );
    my ($textFile, $TABLE_ID, $VARLIST, $ARG_V, $ARG_A, $ARG_O );
    my ($prefix_append, $outFile_append, $SUFFIX);


    # bail if option is NOT set.   
    if ( $APPEND_FITRES_LIST eq "" ) { return ; }

    @wdlist   = split(/\s+/,$APPEND_FITRES_LIST) ;
    $TABLE_ID = $wdlist[0] ;       # first element is the table id
    $textFile = "$MERGED_outFile[$OUTINDX_TEXT][$iver][$ifitopt]" ;
    $VERSION  = "$VERSION_LIST[$iver]" ;
    $MERGEDIR = "$OUTDIR_LIST[$iver]" ;

    # internal suffix for outfiles produced by sntable_dump
    $SUFFIX = "${PREFIX_MERGED}" ;

    # if table id is integer, then it's hbook; else it's root
    if ( (1*$TABLE_ID) eq "$TABLE_ID" ) 
    { $indx =  $OUTINDX_HBOOK ; } 
    else 
    { $indx =  $OUTINDX_ROOT ; }


    # outFile_fit refers to hbook/root output of fit job -> 
    # input to append script
    $outFile_fit    = "$MERGED_outFile[$indx][$iver][$ifitopt]" ;  
    $OUTFILE_FIT    = "$MERGED_OUTFILE[$indx][$iver][$ifitopt]" ;  

    # dumpFile is the intermediate text file containing only the
    # appended variables extracted from outFile_fit
    $dumpFile    = "sntable_dump_${SUFFIX}.fitres" ;
    $logFile     = "sntable_dump_${SUFFIX}.log" ;

    $VARLIST  = "@wdlist[1 .. $#wdlist]" ;
    $ARG_V    = "--v $VARLIST" ;
    $ARG_A    = "--a  $textFile" ;
    $ARG_O    = "--o  $dumpFile" ;
    $CMD = "$JOBNAME_SNTABLE $outFile_fit $TABLE_ID  $ARG_V  $ARG_A  $ARG_O" ;

    print PTR_MERGE2 "\n APPEND_FITRES Command for $VERSION/$SUFFIX: \n";
    print PTR_MERGE2 "$CMD\n";
    PTR_MERGE2 -> autoflush(1);	

    unless ( -e $OUTFILE_FIT ) {
	print PTR_MERGE2 "\t Cannot find '$outFile_fit' for $VERSION ";
	print PTR_MERGE2 " --> skip APPEND_FITRES \n";
	PTR_MERGE2 -> autoflush(1);	
	return ;
    }

    qx(cd $MERGEDIR ; $CMD > $logFile);

    # check for abort (June 26 2015)
    my @bla = qx(cd $MERGEDIR ; grep ABORT $logFile);
    if ( scalar(@bla) > 0 ) {
	$MSGERR[0] = "FATAL ERROR: Append-ABORT found in ";
	$MSGERR[1] = "   $MERGEDIR/$logFile";
	# use die to make sure it gets printed to screen
	die "\n $MSGERR[0]\n $MSGERR[1] \n";
    }


    # set outFile_append = name of appended fitres file
    # (or grep 'Write' from logfile to determine name)
    $prefix_append  = "sntable_append_${SUFFIX}" ;
    $outFile_append = "${prefix_append}.text" ;

    # move appended text file to replace original, 
    # and save original with _ORIG suffix.
    qx(cd $MERGEDIR ; mv $textFile        ${textFile}_ORIG  );
    qx(cd $MERGEDIR ; mv $outFile_append  ${textFile}       );

    # cleanup
    qx(cd $MERGEDIR ; rm ${prefix_append}.*  $logFile $dumpFile ) ;
    qx(cd $MERGEDIR ; gzip ${textFile}_ORIG );

    return ;

} # end of append_text_FITRES


# =========================
sub get_NSN_MERGED {

    my ($iver,$ifitopt) = @_;

    # Created Jul 25 2013
    # get total number of processed SN and fill @NSN_MERGED array.
    # Nominal method is just to grep/count "SN: " lines in the
    # catenated fitres file. For snana jobs that do not have
    # fitres files, grep the log files and add numbers.

    my ( $TEXTFILE, $indx, @tmp, $NSN, $tmpLine, @wdlist);

    $NSN = -999 ;

    $indx = $OUTINDX_TEXT ;
    if ( $OUTFLAG[$indx] ) {
	$TEXTFILE  = "$MERGED_OUTFILE[$indx][$iver][$ifitopt]" ;
	@tmp = qx(grep "SN: " $TEXTFILE);
	$NSN = scalar(@tmp);
    }
    else {
	# grep snana.exe log files
	my $grep_key  = "${qq}Finished process${qq}" ;
	my $iwd_NSN   = 3 ;

	my  $PREFIX_SPLITJOB  = &get_PREFIX_SPLITJOB($iver,$ifitopt,-1);
	my $cdd = "cd $SPLIT_JOBDIR_LCFIT" ;

	@tmp = qx($cdd ; grep $grep_key ${PREFIX_SPLITJOB}_SPLIT*.LOG);
	$NSN = 0 ;
	foreach $tmpLine (@tmp) {
	    @wdlist   = split(/\s+/,$tmpLine) ;
	    $NSN     += $wdlist[$iwd_NSN] ;
	}
    }

    # ------------------------------------------------
    # store NSN in global array for this version/fitopt.
    $NSN_MERGED[$iver][$ifitopt] = $NSN ;

    return ;

} # end of get_NSNTOT


# ======================================
sub  copy_SIMGEN_DUMP { 

    # For simulation, copy SIMGEN_DUMP file to output directory
    # with merged table-files.
    # Dec  8 2014: fix bug defining dumpFile
    # Nov 24 2017: fix to work with gzipped DUMP file

    my ($iver) = @_ ;

    my ($DATADIR, $VERSION, $DUMPFILE, $DUMPFILEgz, $MERGEDIR, $cdMDIR) ;
    my ($tmptup ,$cmd0, $cmd1, $cmd2);

    $VERSION    = "$VERSION_LIST[$iver]" ;
    $DATADIR    = "$DATADIR_LIST[$iver]" ;
    $DUMPFILE   = "${DATADIR}/${VERSION}.DUMP" ;
    $DUMPFILEgz = "${DATADIR}/${VERSION}.DUMP.gz" ;
    $MERGEDIR   = "$OUTDIR_LIST[$iver]" ;
    $cdMDIR     = "cd $MERGEDIR";

    if ( -e $DUMPFILE ) {
	$cmd0     = "cp $DUMPFILE SIMGEN.DAT" ;
	qx($cdMDIR ; $cmd0 ; gzip SIMGEN.DAT );  # Jun 10 2016
    }

    if ( -e $DUMPFILEgz ) {
	$cmd0     = "cp $DUMPFILEgz SIMGEN.DAT.gz" ;
	qx($cdMDIR ; $cmd0 ); 
    }
    
    return ;

} # end of copy_SIMGEN_DUMP

# ======================================
sub mergeLog_update {

    # update MERGE.LOG file and check for errors.
    # If doneFlag  == 0 then write NotDONE status to MERGE-log
    # If doneFlag  == 1 then write final status for merge
    # --> there is always a status for each job.
    #
    # Mar 3 2016: check merge-DONEFILE too

    my ($iver,$ifitopt,$doneFlag) = @_ ;

    my ($CWARN, $indx, $OUTFILE, $DONEFILE, $VERSION, $PREFIX_MERGED );
    my ($EXISTFILE, $C1, $NSN, $CLINK );
    
    $VERSION         = "$VERSION_LIST[$iver]" ;
    $PREFIX_MERGED   = sprintf("FITOPT%3.3d", $ifitopt );

    if ( $doneFlag == 0 ) { goto UPDATE_MERGLOG ; }

    $NSN             = $NSN_MERGED[$iver][$ifitopt] ;

    $CLINK = "" ;
    if ( $IFITOPT_LINK_LIST[$ifitopt] >= 0 )  
    { $CLINK = "->$FITOPT_LINK_LIST[$ifitopt]" ; }

    # make sure that merged hbook/root file(s) actually exist,
    # along with the merged doneFile
    $CWARN = "" ;    
    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) {
	if ( $OUTFLAG[$indx] ) {
	    $OUTFILE   = "$MERGED_OUTFILE[$indx][$iver][$ifitopt]" ;  
	    $DONEFILE  = "$MERGED_DONEFILE[$indx][$iver][$ifitopt]" ;

	    $EXISTFILE = ( (-e $OUTFILE) && (-e $DONEFILE) ) ;
	    $C1        = substr($OUTFILE_FORMAT[$indx],0,1);
	    if ( $EXISTFILE == 0  )   
	    { $CWARN = "${CWARN}-$C1" ; $NMISS_OUTFILE[$indx]++ ; }
	}
    }

    # update MERGE-LOG
    # Create one-line merge-log for each version and fitopt
    # such that catenate MERGE_VERSION* puts them all in the right order
    # regardless of what order they finish.

  UPDATE_MERGLOG:
    my $VTMP = sprintf( "VERSION%3.3d", $iver) ;
    my $MERGELOG_TMP = "$MERGE_LOGDIR/MERGE_${VTMP}_${PREFIX_MERGED}.LOG" ;
    my ($STATUS, $COMMENT);

    if ( -e $MERGELOG_TMP ) { qx(rm $MERGELOG_TMP); }
    open  PTR_MERGE  , "> $MERGELOG_TMP" ;

    if ( $doneFlag ) 
    { $STATUS = "MERGED " ; $COMMENT = "(N=$NSN)  $CWARN"; }
    else
    { $STATUS = "NotDone" ; $COMMENT = "(Running)"; }

    print PTR_MERGE " $STATUS  $VERSION  $PREFIX_MERGED  $COMMENT  $CLINK\n";
    close PTR_MERGE ;

    # Dec 13 2017: grand MERGE summary only on last FITOPT
    # Fragile logic because last FITOPT may not be merged last.
    # Note that mergeLog_cat() is also called later in merge_end
    # in case we miss something here.
    if ( $ifitopt < $NFITOPT-1 ) 
    { return; }

    if ( $doneFlag == 0 && $iver < $NVERSION-1 )
    { return ; }

    # re-make global MERGE.LOG file using catenate so that
    # everything is listed in the correct order.

    &mergeLog_cat();
    
} # end of mergeLog_update

# ======================
sub mergeLog_cat {

    # Dec 20 2017: [code moved from mergeLog_update]
    # create MERGE.LOG with catenate command.
    # Sep 11 2018: append '2>/dev/null' to cat to suppress error messsages
    my ($NCAT, $icat);
    
    qx(rm $MERGELOG; touch $MERGELOG);

    # Dec 12 2017: merge in groups of 100 to avoid "too-long list" error
    $NCAT = int($NVERSION/100);
    for($icat=0;  $icat <= $NCAT; $icat++ ) 
    { qx(cat $MERGE_LOGDIR/MERGE_VERSION${icat}* >> $MERGELOG 2>/dev/null ); }

    return ;
    
}  # end mergeLog_cat 

# ===============================
sub merge_end {

    my ($cn0, $cn1, $CN, $NMISS );

    # re-make MERGE.LOG one last time to be sure it is complete (Dec 20 2017)
    &mergeLog_cat();
    
    # open MERGE LOG  in append mode
    open  PTR_MERGE  , ">> $MERGELOG" ;

    # clean up artificial split-versions
    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
	&clean_SPLIT_VERSION($iver);
    }

    # if text-dump option is on, make grand tar ball of LC text files
    if ( $DOTEXT_TABLE[$ITABLE_LCPLOT] ) {  &clean_text_LCPLOT();  }

    # --------------------------------------
    print PTR_MERGE "\n" ;
    $CN  = sprintf("%4d", $NEXPECT_TOTAL  );
    $cn0 = sprintf("%4d", $NGRACE_TOTAL );
    $cn1 = sprintf("%4d", $NABORT_TOTAL );

    print PTR_MERGE "\t $cn0 of $CN jobs finished with success.\n";   
    print PTR_MERGE "\t $cn1 of $CN jobs ABORTED. \n";
    PTR_MERGE -> autoflush(1);

    if ( $NGRACE_TOTAL != $NEXPECT_TOTAL ) {
	$DONE_STATUS = "FAILURE" ;
	$MSGWARN = "WARNING: some jobs did not finish properly.\n" .
	    "\t (NGRACE=$NGRACE_TOTAL  NEXPECT=$NEXPECT_TOTAL)" ;
	print PTR_MERGE " $MSGWARN \n" ;
	print " $MSGWARN \n" ;
	PTR_MERGE -> autoflush(1);
    }

    
    my ($C1,$indx, $FORMAT);
    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) {
	$NMISS  = $NMISS_OUTFILE[$indx] ;
	$FORMAT = $OUTFILE_FORMAT[$indx];
	$C1     = substr($FORMAT,0,1);
	if ( $NMISS > 0 ) {
	    $DONE_STATUS = "FAILURE" ;
	    print PTR_MERGE 
		"  -$C1 WARNING: $NMISS merged $FORMAT files are missing. \n";
	}
    }

    # compute total elapsed wall time
    my ($t_elapse);
    $TIME_END = time;
    $t_elapse = ($TIME_END - $TIME_START)/3600.0 ; # convert to hr
    $t_elapse = 0.001*int(1000.0*$t_elapse);
    
    print PTR_MERGE "\n" ;
    print PTR_MERGE " Total elapsed wall time: $t_elapse hr. \n" ;
    print PTR_MERGE " Done. \n" ;

    if ( $NFITERR_TOT > 0 ) {
	print PTR_MERGE "\n Found fit-errors: see  /$DEBUG_SDIR  subdir\n";

	print PTR_DEBUG_LIST 
	    " ======================================================\n";
	print PTR_DEBUG_LIST 
	    "\t Found total of $NFITERR_TOT fit-errors. \n";
	print PTR_DEBUG_LIST 
	    "\t See debug commands in DEBUG_COMMANDS.LOG \n";
	print PTR_DEBUG_LIST 
	    " ======================================================\n";
    }

    close PTR_MERGE ;
    close PTR_DEBUG_CMD  ;
    close PTR_DEBUG_LIST  ;

    # give read+write priv to everyone (Jan 18, 2012)
    qx(chmod -R g+rw $OUTDIR);

    print PTR_MERGE2 "\n Finished merging in OUTDIR = \n  $OUTDIR \n" ;
    PTR_MERGE2 -> autoflush(1);

    return ;

} # end of merge_end


# =======================================
sub psnid_FOMsummary {

    # Aug 31, 2013
    # grep out summary information from the log files,
    # sum the stats over each log file, and write grand
    # summary in PTR_PSNID.
    # Note that global PSNID-summary file is created
    # by catenating the [alphabetical] separate psnid-summary
    # files for each version/fitopt.

    my ($iver,$ifitopt) = @_ ;

    my ($cdd, $VERSION, $nnn, $LOGLIST, @bla, $tmpLine, @wdlist, $ctmp );
    my ($wd1, $wd2, $Ntmp, $NTOT, $indx, $strType, $idType, $maxidType) ;
    my (@STRINGTYPE, @IDTYPE, @NTYPE, $TMPLOG, $SUFFIX );

    if ( $PSNID_FLAG == 0 ) { return ; }

    $cdd       = "cd $SPLIT_JOBDIR_LCFIT" ;
    $VERSION   = "$VERSION_LIST[$iver]" ;
    $nnn       = sprintf("%3.3d", $ifitopt);
    $LOGLIST   = "${VERSION}_FITOPT${nnn}_SPLIT*.LOG" ;

    @STRINGTYPE = @IDTYPE = @NTYPE = ();
    $maxidType  = -9 ;
    $NTOT = 0 ;

    # create separate psnid-log, then use catenate to fill
    # the global log
    $SUFFIX  = sprintf("VERSION%3.3d_FITOPT%3.3d", $iver, $ifitopt);
    $TMPLOG  = "$MERGE_LOGDIR/PSNID_SUMMARY_${SUFFIX}.LOG" ;

    open  PTR_PSNID , "> $TMPLOG" ;

    print PTR_PSNID "# =============================================== \n" ;
    print PTR_PSNID "VERSION:    $VERSION   \n";
    print PTR_PSNID "FITOPT${nnn}:  $FITOPT_LIST[$ifitopt] \n";

    PTR_PSNID -> autoflush(1);

    @bla = qx($cdd ; grep PSNID-Type $LOGLIST );
    foreach $tmpLine (@bla) {
	@wdlist   = split(/\s+/,$tmpLine) ;
	$wd1  = $wdlist[1] ;
	$wd2  = $wdlist[2] ;
	$Ntmp = $wdlist[4] ;

	# get string-type from wd1
	$indx    = index($wd1,"=");
	$strType = substr($wd1,$indx+1,12);

	# get integer ID type from wd2: = '(nnn)'
	$ctmp = substr($wd2,1,3);
	$idType = int($ctmp);
	if ( $idType > $maxidType ) { $maxidType = $idType ; }

	$STRINGTYPE[$idType] = $strType ;
	$NTYPE[$idType] += $Ntmp ;
	$NTOT           += $Ntmp ;
#	print PTR_PSNID " strType='$strType'  idType=$idType  N=$Ntmp \n"; 
    }

    # print summary for each type
    my ($pcent, $strLabel);

    for ($idType = -1; $idType <= $maxidType ; $idType++ ) {

	if ( $idType < 0 ) {
	    $strType  = "ALL" ;
	    $strLabel = "ALL" ;
	    $Ntmp     = $NTOT ;
	}
	else {
	    $strType  = $STRINGTYPE[$idType] ;
	    $strLabel = "PSNID-Type=$strType" ;
	    $Ntmp     = $NTYPE[$idType] ;
	}

	if ( length($strType) == 0 ) { next ; }
	
	$strLabel = sprintf("%-20.20s", $strLabel );
	my $Ntmp6 = sprintf("%6d", $Ntmp);
	$pcent = 100.0*$Ntmp / $NTOT ;
	$pcent = sprintf("%6.2f", $pcent);
	print PTR_PSNID "\t N[ $strLabel ] = $Ntmp6  ($pcent %) \n";
    }

    if ( $SIM_FLAG_LIST[$iver] == 0 ) { goto FLUSH ; }

    # -----------------------------------------
    my ($N_SIMGENIa, $N_PSNIDIa_TRUE, $N_PSNIDIa_ALL );
    my (@bla1, @bla2, @bla3, $Nbla, $i, $EFF, $PUR, $FOM );
    my $key1 = "SIMGEN-Ia(ALL)" ;
    my $key2 = "PSNID-Ia(TRUE)" ;
    my $key3 = "PSNID-Ia(ALL)" ;


    print PTR_PSNID "\n    Simulated Figure of Merit: \n" ;

    $N_SIMGENIa = $N_PSNIDIa_TRUE = $N_PSNIDIa_ALL = 0 ;

    
    @bla1 = qx($cdd ; grep "$key1"  $LOGLIST );
    @bla2 = qx($cdd ; grep "$key2"  $LOGLIST );
    @bla3 = qx($cdd ; grep "$key3"  $LOGLIST );
    $Nbla = scalar(@bla1);

    for($i=0; $i < $Nbla ; $i++ ) {
	
	@wdlist = split(/\s+/,$bla1[$i]) ;   $Ntmp  = $wdlist[5];
	$N_SIMGENIa += $Ntmp ;

	@wdlist = split(/\s+/,$bla2[$i]) ;   $Ntmp  = $wdlist[5];
	$N_PSNIDIa_TRUE += $Ntmp ;

	@wdlist = split(/\s+/,$bla3[$i]) ;   $Ntmp  = $wdlist[5];
	$N_PSNIDIa_ALL += $Ntmp ;
    }

    print PTR_PSNID "\t N[ $key1 ] = $N_SIMGENIa \n";
    print PTR_PSNID "\t N[ $key2 ] = $N_PSNIDIa_TRUE \n";
    print PTR_PSNID "\t N[ $key3 ] = $N_PSNIDIa_ALL \n";

    $EFF = $PUR = $FOM = 0.0 ;

    if ( $N_SIMGENIa > 0 ) 
    { $EFF = $N_PSNIDIa_TRUE / $N_SIMGENIa ; }

    if ( $N_PSNIDIa_ALL > 0 ) 
    { $PUR = $N_PSNIDIa_TRUE / $N_PSNIDIa_ALL ; }
    
    $FOM = $EFF * $PUR ;

    $EFF = sprintf("%.3f", $EFF);
    $PUR = sprintf("%.3f", $PUR);
    $FOM = sprintf("%.3f", $FOM);

    print PTR_PSNID "    FOM(Ia) = EFF x PURITY = $EFF x $PUR = $FOM \n";

  FLUSH:
    # print PTR_PSNID "\n";
    PTR_PSNID -> autoflush(1);
    close PTR_PSNID ;

    # catenate all psnid logs into global file (Sep 2013)
    qx(rm $PSNID_SUMLOG) ;
    qx(cd $MERGE_LOGDIR ; cat PSNID_SUMMARY*  >  $PSNID_SUMLOG);

}  # end of psnid_FOMsummary 


# ============================================
sub write_TRAILMARK {

    my($iver, $readmeFile, $DATADIR, $VERSION, $link );

    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
	if ( $SIM_FLAG_LIST[$iver] != $SIM_FLAG_NORMAL ) { next; }
	
	$VERSION    = "$VERSION_LIST[$iver]" ;
	$DATADIR    = "$DATADIR_LIST[$iver]" ;
	$readmeFile = "${DATADIR}/${VERSION}.README" ;
	qx(echo "OUTDIR(split_and_fit): $OUTDIR/$VERSION" >> $readmeFile );
    }

} # end write_TRAILMARK

# ============================================
sub gzip_SPLIT_JOBS {

    # Dec 2014: create SPLIT_JOBS.tar.gz and 'rm -r SPLIT_JOBS'

    my ($CMD, $tarFile, $sdir, @sdirList, $path);

    #bail if there are any warnings in MERGE.LOG
    my @bla = qx(grep WARNING $MERGELOG);
    if ( scalar(@bla) > 0 ) { return ; }
    
    # contruct list of subdirs to gzip
    @sdirList = ( $SPLIT_JOBSUBDIR_LCFIT );
    if ( $SALT2mu_INFILE ne "" ) {
	@sdirList = ( @sdirList, $SPLIT_JOBSUBDIR_SALT2mu );
    }
    if ( $MLFLAG == $MLFLAG_APPLY)  {
	@sdirList = ( @sdirList, $SPLIT_JOBSUBDIR_MLAPPLY );
    }

    foreach $sdir ( @sdirList ) {
	$tarFile = "${sdir}.tar" ;
	qx(cd $OUTDIR; tar -cf $tarFile $sdir ; gzip $tarFile );

	# remove SPLIT_JOBS/ dir only if gzipped tarFile exists
	if ( (-e "$OUTDIR/$tarFile.gz") && length($sdir)>1  ) {
	    use POSIX qw[ close ];
	    POSIX::close( $_ ) for 1 .. 1024; ## Arbitrary upper bound
	    qx(rm $MERGE_LOGDIR/*.*);
	    $path = "$OUTDIR/$sdir" ;
	    rmtree($path);
	}
    }

    return ;

} # end of gzip_SPLIT_JOBS


# ============================================
sub gzip_VERSION {

    my ($iver) = @_ ;

    # Created Dec 2012: 
    #  gzip LOG and his/root files for this version.
    #  Note that gzip includes all FITOP*** for this version.
    #
    # Jun 13 2013: make sure to gzip only $VERSION_FITOPT* and not $VERSION*
    #              in case we have, for example, SDSS and SDSS_bla versions.
    #
    my ($VERSION, $MERGEDIR, $ZIPLIST1, $ZIPLIST3, $indx, $itab );
    my ($suf, @suffixList);

    if ( $GZIP_FLAG          == 0 ) { return ; }  # user flag not set
    if ( $GZIP_FLAG          >  3 ) { return ; }
    if ( $GZIP_STATUS[$iver] == 1 ) { return ; }  # already gzipped

    $VERSION   = "$VERSION_LIST[$iver]" ;
# xxx mark delete    $MERGEDIR  = "$OUTDIR/$VERSION";
    $MERGEDIR  = "$OUTDIR_LIST[$iver]";

    my $dashLine = "- - - - - - - - - - - - - - - - - - - - - - - - - - - " ;
    print PTR_MERGE2 "\n" ;
    print PTR_MERGE2 "\t# $dashLine \n" ;
    print PTR_MERGE2 "\t gzip files for VERSION = '$VERSION' \n";
    print PTR_MERGE2 "\t# $dashLine \n" ;
    print PTR_MERGE2 "\n" ;
    PTR_MERGE2 -> autoflush(1);

    # always gzip SPLIT_JOBS if flag is set
    
    $ZIPLIST1 = "${VERSION}_FITOPT*.LOG" ;  # Jun 13 2013
    $ZIPLIST3 = "" ;
    for($indx=0; $indx < $NMAX_OUTINDX; $indx++ ) {
	if ( $OUTFLAG[$indx] == 0 ) { next ; }

	# get ziplist for SPLIT_JOBS
	$suf      = "$SUFFIX_OUTFILE[$indx]" ;
	$ZIPLIST1 = "$ZIPLIST1 ${VERSION}_FITOPT*.${suf}" ;

	# get ziplist for merged files
	if ( $indx < $OUTINDX_TEXT ) { 	    
	    $suf      = "$SUFFIX_OUTFILE[$indx]" ;
	    $ZIPLIST3 = "$ZIPLIST3 FITOPT*.${suf}" ;
	}
	else {
	    for ( $itab=0; $itab < $MAXTABLE; $itab++ ) {
		if ($DOTEXT_TABLE[$itab] )  { 
		    $suf      = "$TABLELIST[$itab]" ;
		    $ZIPLIST3 = "$ZIPLIST3 FITOPT*.${suf}";
		} 
	    } # end itab

	} # end if-block

    }  # end indx loop
    
    qx(cd $SPLIT_JOBDIR_LCFIT ; gzip ${ZIPLIST1} 2>/dev/null ) ; 

    # check option to gzip merged files
    if ( $GZIP_FLAG == 3 )  { qx(cd $MERGEDIR ; gzip $ZIPLIST3 2>/dev/null); }

    # set GZIP_STATUS flag to done so that we don't gzip again
    $GZIP_STATUS[$iver] = 1 ;

} # end of gzip_VERSION


# ==================================
sub run_afterBurner {
    my ($iver) = @_ ;  # Aug 31 2015

    # run user-specified VERSION_AFTERBURNER commane in version $iver.
    # This command operates on all merged fitres/hbook/root files 
    # under this version.

    my ($VERSION, $MERGEDIR, $CMD);

    $CMD       = "$VERSION_AFTERBURNER";
    $VERSION   = "$VERSION_LIST[$iver]" ;
# xxx mark delete    $MERGEDIR  = "$OUTDIR/$VERSION";
    $MERGEDIR  = "$OUTDIR_LIST[$iver]";

    if ( length($CMD) == 0 ) { return ; }
    qx(cd $MERGEDIR; $CMD);

    print PTR_MERGE2 "\n" ;
    print PTR_MERGE2 "# ================================================== \n" ;
    print PTR_MERGE2 " Run AFTERBURNER command: \n\t $CMD\n" . 
	" for $VERSION \n" ;
    print PTR_MERGE2 "# ================================================== \n" ;
    print PTR_MERGE2 "\n" ;

    PTR_MERGE2 -> autoflush(1);

}   # end run_afterBurner

# ==================================
sub get_suffixMerged {
    
    my ($indx,$itab) = @_ ;

    # return merged suffix for this format-indx and table-index.
    # Only the text suffix depends on the table index.    

    my $SUF ;
    if ( $indx != $OUTINDX_TEXT )
    {  $SUF  = "$SUFFIX_OUTFILE[$indx]" ; }
    else
    {  $SUF  = "$TABLELIST[$itab]" ; }

    return $SUF ;

} # end of get_suffixMerged

# ==============================
sub clean_SPLIT_VERSION  {

    my ($iver) = @_ ;

    # remove SPLIT*_[VERSION] for TEXT-format only !
    # Note that the simulation SPLIT-version is a directory,
    # while the data split-version consists of files.
    # Remove only if in TEXT format that results in an artificial
    # version for each split job.

    my ($FORMAT, $VERSION, $DATADIR, $rm_prefix );

    $FORMAT   = $DATA_FORMAT[$iver];


    # check reasons to bail
    if ($FORMAT ne "TEXT") { return ; }
    if ($OPT_SPLIT == 0  ) { return ; }

    # -----------------------------------
    # now remove the artificial split-versions.
    $VERSION   = "$VERSION_LIST[$iver]" ;
    $DATADIR   = "$DATADIR_LIST[$iver]" ;
    $rm_prefix = "${SNDATA_ROOT}/lcmerge/SPLIT*_${VERSION}" ;

    qx(rm -r ${rm_prefix});   # rm SPLIT-dir



} # end of sub clean_SPLIT_VERSION 


# =========================
sub clean_text_LCPLOT {

    # July 23, 2013
    # We get here at very end if LDMP_SNFLUX = T in &SNLCINP nml.
    # There is one text file per SN and one LIST file per job,
    # so put these all into one big tarball so that the SPLIT_JOBS
    # sub-dir is not such a big mess.
    # Note that files from each VERSION_FITOPT have already
    # been catenated in the ../VERSION subdirs.

    my ( $fList, $cmdTar, $tarFile );

    print PTR_MERGE2 "\n" ;
    print PTR_MERGE2 " Create tarball of  LCPLOT  files.\n" ;
    PTR_MERGE2 -> autoflush(1) ;

    # list allows for optional .gz suffix in case GZIP_FLAG is set.
    $fList   = "*.${SUFFIX_LCPLOT_TEXT}*  *.${SUFFIX_LCLIST_TEXT}*" ;
    $tarFile = "TEXTSNLC.tar " ;

    $cmdTar = "tar -cf $tarFile $fList" ;
    qx(cd $SPLIT_JOBDIR_LCFIT ; $cmdTar ; gzip $tarFile; rm $fList );


} # end of &clean_text_LCPLOT

# ===============================
sub combine_all_text {

    my ($ifitopt) = @_ ;

    # if COMBINE_ALL_FLAG is set then catenate fitres files
    # from all versions into one file (for input $ifitopt) 
    # in a new directory.
    # This option is useful to combine simulations, for example,
    # from different seasons or from different random seeds.
    #

    if ( $COMBINE_ALL_FLAG == 0 ) { return ; }

    my ($FRES_COMBINED, $FFF, $FRES, $COMBINE_DIR, $FRES_COMBINED );

    # create new 'combined' directory 
    $COMBINE_DIR   = "$OUTDIR/${COMBINE_ALL_PREFIX}" ;

    $FFF           = sprintf("FITOPT%3.3d", $ifitopt );
    $FRES          = "${FFF}.FITRES" ;
    $FRES_COMBINED = "$COMBINE_DIR/${COMBINE_ALL_PREFIX}_${FFF}.FITRES" ;

    # use 'cat' to combine all fitres files
    qx(cd $OUTDIR; cat */$FRES > $FRES_COMBINED );

    # -----------------------------------------
    # add comments at top
    my (@comment_add, $comment, $NTOT, @tmp, $sedcmd, $TMPFILE, $j );
    my $PND = '#' ;
    my $KEY = "COMBINE_ALL_COMMENT:" ;

    @tmp = qx(grep "SN: " $FRES_COMBINED);
    $NTOT = scalar(@tmp);
    $j = -1;
    
    $j++ ; $comment_add[$j] = 
	"$KEY COMBINE_ALL_FLAG = $COMBINE_ALL_FLAG";
    $j++ ; $comment_add[$j] = 
	"$KEY Catenated from $NVERSION versions.";
    $j++ ; $comment_add[$j] = 
	"$KEY NSN_COMBINED = $NTOT";
    $j++ ; $comment_add[$j] = 
	"============================================================ " ;
    $j++ ; $comment_add[$j] = 
	"" ;
    
    $sedcmd = "sed ";
    foreach $comment (@comment_add)
	{ $sedcmd = "$sedcmd " . "-e '1i $PND $comment'" ; }

    $TMPFILE = "TMP_COMBINED_${FRES}" ;
    qx(cd $OUTDIR; $sedcmd $FRES_COMBINED > $TMPFILE ; mv $TMPFILE $FRES_COMBINED );


} # end of combine_all_text


sub combine_all_outFiles {
    
    my ($indx,$ifitopt) = @_ ;

    # Jun 2013
    # merge hbook/root files for each version ($iver) as if all
    # versions were a single job. Useful for combining sims with
    # different random seeds.
    # Input $indx indicates outFile format: hbook, ...    

    if ( $COMBINE_ALL_FLAG == 0 ) { return ; }
    if ( $OUTFLAG[$indx]   == 0 ) { return ; }

    my ($COMBINE_DIR, $SUFFIX, $FMT );

    $COMBINE_DIR   = "$OUTDIR/${COMBINE_ALL_PREFIX}" ;
    $SUFFIX        = "$SUFFIX_OUTFILE[$indx]" ;
    $FMT           = "$OUTFILE_FORMAT[$indx]" ;

    # create COMBINE-dir on first fit-option
    if ( $ifitopt == 0 ) {

	if ( !( -d $COMBINE_DIR) )  { qx(mkdir $COMBINE_DIR); }
    
	print PTR_MERGE "\n" ;
	print PTR_MERGE " Catenate all $FMT files for each version in \n" ;
	print PTR_MERGE "   $COMBINE_DIR  \n" ;
	PTR_MERGE -> autoflush(1);
    }

    if ( $indx == $OUTINDX_TEXT ) {
	&combine_all_text($ifitopt); 
	return ;
    }

    # -------- do the grand merge/combine ---------

    my ($iver, $FFF, $OUTFILE, $JOBNAME, $CMD);
    my ($MERGE_LIST, $MERGE_OUTFILE, $MERGE_LOG);
    
    $MERGE_LIST = "" ;
    $FFF        = sprintf("FITOPT%3.3d", $ifitopt );

    for ( $iver=0 ; $iver < $NVERSION; $iver++ ) {
	$OUTFILE    = "$MERGED_OUTFILE[$indx][$iver][$ifitopt]";
	$MERGE_LIST = "$MERGE_LIST $OUTFILE" ;
    }  

    $MERGE_OUTFILE = "${COMBINE_DIR}/COMBINED_${FFF}.${SUFFIX}";
    $MERGE_LOG     = "${COMBINE_DIR}/MERGE_${FFF}.LOG";
    $JOBNAME       = "$JOBNAME_MERGE[$indx]";
    $CMD   = "$JOBNAME $MERGE_LIST $MERGE_OUTFILE >& $MERGE_LOG" ;

    qx($CMD);

    return ;

} # end of combine_all_outFiles


# ========================
sub error_check_LCFIT {

    my ($iver, $ifitopt) = @_ ;

    # check fit-jobs for errors; if errors are detected
    # then create DEBUG_${VERSION} with just those SNe
    # with errors (up to NERRMAX=100) so that it's easier
    # run a debug job. Also prepare fit-namelist pointing
    # to this version and with the correct options.
    #
    # Sep 13, 2014: fix bug, make sure DEBUG_DIR is created only once.

    if ( $ERRCHECK_FLAG == 0 ) { return ; }

    my ($fresFile, $NCID, @CIDLIST_DEBUG, @CIDLIST_NaN );
    my ($DEBUG_FILE);

    if ( $FITOPT_LINK_LIST[$ifitopt] ne "" ) { return ; }

    # create /DEBUG subdir under $OUTDIR
    if ( $iver == 0 &&  $ifitopt == 0 ) {
	$DEBUG_SDIR = "DEBUG";
	$DEBUG_DIR  = "$OUTDIR/$DEBUG_SDIR";

	if ( !(-d $DEBUG_DIR ) ) {
	    qx(mkdir $DEBUG_DIR);

	    $DEBUG_FILE = "$DEBUG_DIR/DEBUG_LIST.LOG" ;
	    open PTR_DEBUG_LIST , "> $DEBUG_FILE" ;
	    
	    $DEBUG_FILE = "$DEBUG_DIR/DEBUG_COMMANDS.LOG" ;
	    open PTR_DEBUG_CMD , "> $DEBUG_FILE" ;
	}
    }

    $fresFile = $MERGED_OUTFILE[$OUTINDX_TEXT][$iver][$ifitopt] ;

    # --------------------
    # check for nan
    @CIDLIST_NaN  = &fitresNaN($fresFile);
    $NFITERR_NaN += scalar(@CIDLIST_NaN);

    # now check for ???

    # --------------------
    # glue together all sources of error.
    @CIDLIST_DEBUG = ( @CIDLIST_NaN );

    # bail if there are no errors.
    $NCID = scalar(@CIDLIST_DEBUG);
    if ( $NCID == 0 ) { return ; }

    $NFITERR_TOT += $NCID ;

} # end of  error_check_LCFIT


# =============================
sub fitresNaN(@) {

    # return CIDLIST for SN entries that have 'nan' somewhere.
    # grep for " nan " instead of "nan" to avoid confusion with 
    # SNANA comment.

    my ($fresFile) = @_ ;
    my (@tmpLines, @CIDLIST, $CID, $NSN, $line, @wdlist );

    @tmpLines = qx(grep -i " nan " $fresFile);
    $NSN      = 0 ;

    # strip of CID from each line that has 'nan'
    foreach $line ( @tmpLines ) {
	@wdlist   = split(/\s+/,$line) ;
	$CID      = "$wdlist[1]" ;
	$CIDLIST[$NSN]  = $CID ;
	print PTR_DEBUG_LIST "Found NaN for CID = $CID  in \n";
	print PTR_DEBUG_LIST "  $fresFile \n\n";
	$NSN++ ;
    }

    return @CIDLIST ;

} # end of fitresNaN


# ===========================================
sub check_TIME_STAMP {

    # Read DATE stamp from file and compare with CDATE_SUBMIT
    # passed as split_and_fit argument. 

    my ($STAMP_FILE, $CDATE_FILE, $STOP_FILE);

    # skip date-stamp check for interactive MERGE option
    if ( $CDATE_SUBMIT eq "UNKNOWN" ) { return ; }

    $STAMP_FILE = "${SPLIT_JOBDIR_LCFIT}/${TIME_STAMP_FILENAME}" ;
    $CDATE_FILE = `cat $STAMP_FILE` ; 
    $CDATE_FILE    =~ s/\s+$// ;   # trim trailing whitespace

    if ( "$CDATE_FILE" eq  "$CDATE_SUBMIT" ) { return; }

    # if we get here, this "split_and_fit MERGE" is most likely
    # rogue/leftover job, so abort and leave message in file.

    $STOP_FILE = "$OUTDIR/STOP-split_and_fit_${CDATE_SUBMIT}.LOG" ;
    open  PTR_LOG  , "> $STOP_FILE" ;  
    print PTR_LOG " Time stamp mis-match:  \n" ;
    print PTR_LOG "   Time of job submit:   '$CDATE_SUBMIT' \n" ;
    print PTR_LOG "   Time in OUTDIR file : '$CDATE_FILE'   \n" ;
    print PTR_LOG "\n";
    print PTR_LOG " STOP process:\n   $0 @ARGV \n" ;
    print PTR_LOG "\n";
    print PTR_LOG " BEWARE: " . 
	"there could be merge-interference in this output ! \n" ;
    print PTR_LOG "\n";

    #close PTR_LOG ;
    
    die "\n STOP lost split_and_fit job ($CDATE_SUBMIT)\n";

#    print " xxx ------------------------------------------------- \n" ;
#    print " xxx CDATE[SUBMIT,FILE] = '$CDATE_SUBMIT' ,  '$CDATE_FILE' \n";

    # .xyz
    return ;

} # end check_TIME_STAMP

# ========================================
sub check_waitTime {

    # Created Jan 30 2017
    # increment wait time ; if too long then abort.

    my ($CDATE, $t_hr);
    $CDATE  = `date +%Y%m%d_%H%M` ;
    $CDATE  =~ s/\s+$// ;   # trim trailing whitespace

    if ( $TOTAL_WAIT > ($TOTAL_WAIT_ABORT*3600.) ) {
	$t_hr = $TOTAL_WAIT/3600.0;
	$t_hr = .01*int(100.*$t_hr) ;
	$MSGERR[0] = "\n";
	$MSGERR[1] = " Total $SCRIPTNAME wait time ($t_hr hr)\n" ;
	$MSGERR[2] = " exceeds limit \n";
	$MSGERR[3] = " --> KILL $SCRIPTNAME background job at $CDATE \n";
	$MSGERR[4] = " BEWARE: batch jobs may still be running on nodes.\n";

	open  PTR_MERGE  , ">> $MERGELOG" ;
	print PTR_MERGE  "@MSGERR";
	PTR_MERGE -> autoflush(1);
#	close PTR_MERGE ;

	print PTR_MERGE2 "@MSGERR \n";

	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    return ;

} # end check_waitTime ;

# ========================================
sub CLEANFILES_DRIVER {

    # Created Jan 2017
    # use 'find' to find all files of the format,
    #     SPLIT_JOBS\*, SIMGEN.DAT\* and ORIG.gz 
    # which were created by split_and_fit.
    # This command drills into all sub-directories 
    # from where "split_and_fit.pl CLEAN" is run.
    #
    # After deleting files, ask user to gzip FITOPT\*.FITRES files.
    #
    # First a list of all files to delete is printed;
    # then use must respond 'yes' to delete them.
    #
    # ----------

    my $CLEANMASK_SPLIT_JOBS = 1 ;  # SPLIT_JOBS_XXX/ and misc
    my $CLEANMASK_SIMGEN     = 2 ;  # SIMGEN.DAT files
    my $CLEANMASK_gzipFITOPT = 4 ;  # FITOPT*FITRES and FITOPT*LCPLOT
    my $CLEANMASK_HBOOKROOT  = 8 ;  # HBOOK and ROOT files

    my $CLEANMASK_DEFAULT =
	$CLEANMASK_SPLIT_JOBS +
	$CLEANMASK_SIMGEN     + 
	$CLEANMASK_gzipFITOPT ;
	
    my ($cmd_find, @tmp, $ftmp, $nline, @wdlist, $response);
    my ($size, $SIZE_TOT );

    my @FLIST_SPLIT_JOBS = ( 
	"FITOPT\*ORIG.gz",    "simLogs" , "DEBUG" ,
	"SPLIT_JOBS",         "SPLIT_JOBS.tar.gz", 
	"SPLIT_JOBS_LCFIT",   "SPLIT_JOBS_LCFIT.tar.gz", 
	"SPLIT_JOBS_MLAPPLY", "SPLIT_JOBS_MLAPPLY.tar.gz", 
	"S2mu\*.text" );

    my @FLIST_SIMGEN = ( "SIMGEN.DAT.gz" );

    my @FLIST_HBOOKROOT = (  "FITOPT\*HBOOK", "FITOPT\*ROOT" );


    if ( $CLEANMASK == -1 ) { $CLEANMASK = $CLEANMASK_DEFAULT ; }
    print " Apply CLEANMASK = $CLEANMASK \n";

    if ( $CLEANMASK & $CLEANMASK_SPLIT_JOBS ) 
    { sntools::clean_fileList($PROMPT_FLAG,@FLIST_SPLIT_JOBS); }

    if ( $CLEANMASK & $CLEANMASK_SIMGEN ) 
    { sntools::clean_fileList($PROMPT_FLAG,@FLIST_SIMGEN); }

    if ( $CLEANMASK & $CLEANMASK_HBOOKROOT ) 
    { sntools::clean_fileList($PROMPT_FLAG,@FLIST_HBOOKROOT); }

    # - - - - - - - -  -
    # ask to compress FITRES files.

    if ( $CLEANMASK & $CLEANMASK_gzipFITOPT ) 
    {  &clean_gzipFITOPT(); }

    print " Done. \n";

} # end CLEANFILES_DRIVER



sub clean_gzipFITOPT {

    my ($cmd_find, $response, @tmp );

    if ( $PROMPT_FLAG == 0 ) {
	$response = 'y' ;
    }
    else {
	print "\n Enter 'y' to gzip all " . 
	    "FITOPT*.FITRES and FITOPT*.LCPLOT files => \n  "; 
	$response = <STDIN> ;     $response =~ s/\s+$// ;
    }

    if ( $response eq "y" ) {
	printf "\t gzipping FITOPT*.FITRES . . . \n";
	$cmd_find = "find . -name FITOPT\*.FITRES -exec gzip {} +";
	@tmp = qx($cmd_find);

	printf "\t gzipping FITOPT*.LCPLOT . . . \n";
	$cmd_find = "find . -name FITOPT\*.LCPLOT -exec gzip {} +";
	@tmp = qx($cmd_find);
    }
    else {
        print " FITOPT Files NOT gzipped.\n" ;
    }

} # end clean_gzipFITOPTs
