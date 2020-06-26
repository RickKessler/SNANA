#!/usr/bin/perl
#
# Created April 2013 by R.Kessler
#
# Run SALT2mu on text/fitres output from split_and_fit.pl.
#
#
# USAGE:
#   SALT2mu_fit.pl  <inputFile>
#   SALT2mu_fit.pl  <inputFile>  NSPLITRAN=<n>  ! divide data into <n> subsets
#   SALT2mu_fit.pl  <inputFile>  NOSUBMIT
#   SALT2mu_fit.pl  <inputFile>  NOPROMPT
#   SALT2mu_fit.pl  <inputFile>  KILL
#   SALT2mu_fit.pl  <inputFile>  SUMMARY  ! make sumary only
#   SALT2mu_fit.pl  <inputFile>  INPDIR+  <DIRNAME1>  INPDIR+  <DIRNAME2>
#   SALT2mu_fit.pl  <inputFile>  CLEAN    ! remove HBOOK, gzip SALT2*FITRES
# This script reads the following additional header keys
# from the SALT2mu input file:
#
#  To process input directories that are independent:
#    JOBNAME:  <optional SALT2mu in user path for debugging>
#    INPDIR: <dir1>           ! ENV are allowed; e.g. $TOPDIR_FITS/systematics
#    INPDIR: <dir2>  bla=xyz  ! fitopt only for this INPDIR
#
#    ....
#    INPDIR: <dirN>
#
#    MUOPT:  bins=20               # apply to all FITOPT by default
#    MUOPT:  bins=30  FITOPT=0     # apply only to FITOPT000
#    MUOPT:  bins=30  FITOPT=0,4   # apply only to FITOPT000 & FITOPT004
#    MUOPT:  [ZDUM] nzbin=20       # add optional nick name in [] for makeCov
#
#    FITOPT000_ONLY: 1   # process only default FITOPT000
#
#    VERSION_EXCLUDE:         <version>  ! exclude this version
#    VERSION_EXCLUDE_STRING:  <version>  ! exclude all versions containing
#                                        ! this string
#    
#    NODELIST:  <node1> <node2> ... <nodeN>
#         or
#    BATCH_INFO:  <command>  <templateFile>  <Ncore>
#
#    DONE_STAMP: <fileName>  ! if not specified, default ALL.DONE is created
#
#    WFITMUDIF_OPT:  <wfit options for fitting M0 vs. z>
#    WFIT_OPT:   <wfit options> MUDIF          ! wfit only MUDIF
#    WFIT_OPT:   <wfit options> MUDIF FITRES   ! wfit both  MUDIF & FITRES
#
#  To run SALT2mu on summed FITRES files from different directories.
#    INPDIR+: <dir1>   ! dir1 created by split_and_fit
#    INPDIR+: <dir2>   ! dir2 created by split_and_fit
#    ....
#    INPDIR+: <dirN>
#
#  To fix one version under a survey with INPDIR, give the full path 
#  to the directory containing the FITRES files. E.g., 
#     INPDIR+:  $SIMFITS/LOWZ/set3
#     INPDIR+:  $SIMFITS/HIZ
#  will loop over each HIZ-survey version and combine with the same 
#  LOWZ/set3 version.
#
# SALT2mu can be run multuple times with different options; e.g., 
#    MUOPT: nzbin=20
#    MUOPT: powzbin=3
#    MUOPT: CUTWIN x1ERR 0 1
#
#    OUTDIR_OVERRIDE: <TOPDIR>  
#    OUTDIR_PREFIX:  SALT2mu    ! prefix for each output directory 
#
#    STRINGMATCH_IGNORE:  LOWZ  SDSS  SNLS
#        (see explanation below)
#
#    STRINGMATCH_IGNORE:  IGNORE  
#       (ignore string matches; requires 1 and only 1 version per INPDIR)
#
#    STRINGMATCH:  <COMBINED>   # combine all versions with COMBINED in name
#    STRINGMATCH:  IGNORE       # same as "STRINGMATCH_IGNORE:  IGNORE"
#
# STRINGMATCH_IGNORE is used to help match which versions to sum with 
# the 'INPDIR+:' keys. For example, if versions LOWZ_SYST1 and SDSS_SYST1 
# are to be summed, then user must specify 'STRINGMATCH_IGNORE: LOWZ SDSS'  
# so that  the "SYST1" parts of the version name can be matched. If
# multiple string matches are found, the script aborts. If no string
# matches are found for a particular version, this version is ignored
# and other versions are checked.
#
# The output directory is determined internally so that there is one
# less [OUTDIR] key to specify. For each INPDIR there is a parallel
# OUTDIR with SALT2mu appended as a prefix. Thus if INPDIR = /BLA/simfits
# then the SALT2mu fits will be in /BLA/SALT2mu_simfits. For summed
# fitres inputs the output subdir names are glued with a plus sign; 
# thus for
#    INPDIR+:  /BLA/LOWZ
#    INPDIR+:  /BLA/SDSS
#    INPDIR+:  /BLA/SNLS
# the SALT2mu fit results are written to  /BLA/SALT2mu_LOWZ+SDSS+SNLS .
# If the OUTDIR_OVERRIDE key is given, this path overrides
# /BLA/SALT2mu_LOWZ+SDSS+SNLS .
#
# This script can operate interactively, ssh-ing into a list of
# nodes specified by NODELIST, or by using a batch system specified
# by the BATCH_INFO keys. The ssh & batch keys work the same way as
# for split_and_fit.pl and sim_SNmix.pl .
#
# NSPLITRAN=<n> argument operates on the native file=<file> argument 
# in the SALT2mu-input file, and ignores the supplmental batch keys 
# INPDIR, INPDIR+,  and MUOPT. Batch keys NODELIST, BATCH_INFO and 
# SNANA_LOGIN_SETUP are used. If NSPLITRAN=<n> appears inside the 
# SALT2mu-input file, this is equivalent to the command-line specifier; 
# command-line NSPLITRAN=<n> overrides value in SALT2mu-input file. 
# The NSPLITRAN mode creates an output sub-directory with a specific 
# name, "OUT_SALT2mu_NSPLITRAN[n]", and can be changed using the 
# "OUTDIR_OVERRIDE: <OUTDIR>"  key. The output file names are hard-wired 
# with prefix=SALT2mu. Note that for <n> split jobs, <n+1> batch jobs
# are submitted because the n+1'th job reads back all of the previous
# output and prepares a summary of averages and RMS in SALT2mu_summary.out.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# History:
#
# Oct 2 2013:
#  Fix  makeSumDir_SALT2mu to grep out FITOPT000 from MERGED.LOG
#  the same way as in makeDir_SALT2mu. 
#
# Apr 26, 2014:  .his -> .hbook
#
# Apr 12, 2016: few fixes to read newer SALT2mu output, and new
#               key VERSION_EXCLUDE to exclude a version.
#
# May 25 2016: hbook -> HBOOK  && root->ROOT  in  @TABLE_SUFFIX_LIST
# Jun 02 2016: new option WFITMUDIF_OPT to run wfit on M0 vs. z
# Jun 07 2016: gzip logs under FITJOBS_xxx
# Jun 20 2016: allow ENV in INPDIR arg.
# Jul 15 2016: allow INPDIR-specific option after each input dir.
#              See usage above
# Jul 22 2016: allow ENV in JOBNAME
# Aug 18 2016: new NOPROMPT command-line option to NOT check for
#               already-running SALT2mu jobs
#
# Sep 13 2016: fix subtle grep bug. grep " FITOPT000 " instead of 
#              grep FITOPT000 since "->FITOPT000" can appear in the
#              MERGE.LOG file.
#
# Sep 26 2016: 
#  + new input OUTDIR_PREFIX (default = "SALT2mu" )
#  + fix bug parsing STRINGMATCH_IGNORE (read up to 99 instead of 9 strings)
#  + abort if INPDIR ends with slash
#
# Dec 4 2016: 
#   + bug-fix so that SNANA_LOGIN_SETUP is applied correctly using NODELIST.
#   + write 2nd summary file that is code-parsable (PTR_SUMDAT)
#      See new function write_SUMMARY_DAT().
#
# Jan 15 2017: 
#   + remove underscore or dash from subDir name; see chop func
#   + for "INPDIR+" option, replace  mkdir with 'mkdir -p'
#
# Mar 12 2017: 
#  + if FITRES files are gzipped, unzip them (See $fgzip).
#    Works on directories cleaned with 'split_and_fit.pl CLEAN'
#    
#  + $OUTDIR is based on longest matching subDir instead of
#    longest matching string --> avoids creating new outdir
#    that includes common prefix of split_and_fit subDirs.
#
#  + common prefixes are removed to construct outDir name:
#    OUTFIT_LOWZ & OUTFIT_SDSS -> OUTDIR = SALT2mu_LOWZ+SDSS
#
# Mar 26 2017: 
#   + add NOSUBMIT option
#   + if no string matches, write diagnostic info to more easily
#     find the version with no string match ; see matchDump().
#
# Apr 23 2017:
#   + gunzipping FITRES files (mar 12 fix) works for INPDIR+.
#     Apply fix in copy_inpFiles to work with "INPDIR:" key.
#
# May 17 2017: 
#   + change FITOPT key -> MUOPT key
#   + new feature to give FITOPT-integer list for each MUOPT;
#     see MUOPT syntax above.
#   + minor refactor to optimally distribute jobs among cores
#
# Jun 25 2017: 
#   + increase memory from 500mb to 1000mb
#   + new keys WFIT_OPT & FITOPT000_ONLY
#
# Jul 13 2017:  increase memory from 1000 -> 2000mb (for large FITRES files)
# Jul 19 2017:  NVERSION_MAX -> 9999 (was 999)
# Aug 09 2017:  new key VERSION_EXCLUDE_STRING
# Aug 25 2017:  write beta1 and its error to summary-dat file
# Sep 04 2017:  fix to work with or without leading space in MERGE.LOG
# Nov 02 2017:  new key STRINGMATCH to explicitly match version strings 
# Nov 08 2017:  command line arg 'KICP' or 'kicp' to force kicp q
#
# Nov 30 2017:  
#   + allow command line args "INPDIR+ <DIRNAM1> INPDIR+ <DIRNAME2> ... "
#   + allow command-line arg  "STRINGMATCH <string>"
#   + write gamma0[_err] to SUMMARY.DAT file
#
# Dec 1 2017: 
#   + if input FITRES files are gzipped, leave them gzipped.
#     For gz files, use zcat. parse_lines function in sntools.pm
#     was modified to use zgrep on .gz files.
#
#   + no longer create HBOOK file by default. Maybe later add option.
#
# Dec 3 2017
#   + add CLEAN option
#
# Dec 20 2017: fix VERSION_STRINGMATCH_IGNORE to insert backslash in
#                front of + in directory name.
#
# Dec 21 2017: in prep_COMMAND(), same prefix with or without MUOPT
# Dec 22 2017: 
#   + add new option to fix one version under INPDIR+
#
# Jan 30 2018:
#   + update so that STRINGMATCH works with INPDIR keys
#     (before it only worked with INPDIR+)
#
# Aug 25 2018: for CLEAN command, include lower-case SALT2mu*hbook
# Jan 16 2019: for CLEAN command, make sure to exit(0)
# Feb 04 2019: abort for interactive mode.
# May 24 2019: new OUTDIR_OVERRIDE key to force output location
# May 30 2019: works with NSPLITRAN=n
# Jun 10 2019: 
#   + count aborts and report SUCCESS or FAIL in ALL.DONE
#   + for INPDIR+ or OUTDIR_OVERRIDE, FITJOBS/ subdir is created
#      under OUTDIR so that everything is one place.
#
# Sep 12 2019: 
#  + add optional DONE_STAMP key to override default ALL.DONE file.
#  + start replacing parse_line with parse_array so that comment
#    lines are ignored.
# Sep 13 2019: add batchName arg to make_batchFile() 
#
# Sep 14 2019:
#  + update to work with FITOPT[nnn].FITRES links that are full paths
#    instead of only file links in same directory. Works now for both.
#  + old FITJOBS subdir renamed to SALT2mu_FITSCRIPTS. 
#  + if STRINGMATCH_IGNORE results in null string, set match-string to 
#    SALT2mu_FITJOBS, which will be the name of the subDir with fit 
#    job logs and output.
#
# Sep 18 2019: fix to work with input filename that includes full path.
# Sep 28 2019:
#   + add quotes around 2nd arg to wait_for_files.pl
#   + new function create_done_file to ensure everything else is done
#     before final DONE file is created.
#   + submit "SALT2mu.exe SUMMARY" job into background and 
#     return unix control. No more "hanging" until it's done,
#     and no more need to pipe stdout to get back control.
#
# Oct 4 2019: if OUTDIR_OVERRIDE has no slash, append LAUNCH_DIR
# Oct 25 2019: OUTDIR_OVERRIDE works for INPDIR too 
# Dec 09 2019: fix bug in makeDir_NSPLITRAN();  check SUMMARY_FLAG
#
# Mar 10 2020: 
#   + new option to ignore string matches if there is only one version;
#     see new function require_one_version().
#
# May 10 2020:
#   To combine input FITRES files, replace unix cat with 
#   SALT2mu.exe cat_only <args> to allow for different columns
#   in each input FITRES file specified by datafile=<fileList>/
#   Later, this SALT2mu.exe call may be replaced by a python
#   wrapper, cat_snana_table.py.
#
# May 27 2020:
#   + in makeDir_NSPLITRAN(), make better abort message when
#     datafile= is missing.
#
# Jun 4 2020:
#   +  replace some FATAL_ERROR calls with FATAL_ERROR_STAMP so that
#      a SUCCESS or FAIL message is passed to Pippin.
#
# Jun 11 2020: in NSPLITRAN_prep_COMMAND(), include wfit commands.
#
# Jun 24 2020: refactor so that NSPLITRAN output goes into separate
#              sub-directories labeled SPLITRAN-[nnnn]/
#
# ------------------------------------------------------

use IO::Handle ;
use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

my $DEBUG_submit_SUMMARY = 1 ;  # internal/temp debug flag

my $NVERSION_MAX = 9999 ;   # nominal = 999; less for debugging

my  $SCRIPTNAME       = "SALT2mu_fit.pl" ;
my  $SCRIPTNAME_FULL  = "$0";

my  $OPT_ABORT = 0 ;
my  $OPT_WARN  = 1 ;
my  $OPT_QUIET = 2 ;
my  $JOBNAME_FIT        =  "SALT2mu.exe";
my  $JOBNAME_MAKETABLE  =  "combine_fitres.exe" ;
my  $JOBNAME_WFIT       =  "wfit.exe" ;
my  $JOBMEMORY          =  "2000mb" ;

my  $OUTPUT_SUFFIX_fitres = "fitres" ;  # output from SALT2mu
my  $OUTPUT_SUFFIX_log    = "log"    ;  # output from SALT2mu
my  $OUTPUT_SUFFIX_mudif  = "M0DIF"  ;

# files/dirs already created by split_and_fit.pl
my  $MERGE_FILENAME    = "MERGE.LOG" ;      
my  $FITOPT_README     = "FITOPT.README";

# work dirs
my  $LAUNCH_DIR     = `pwd`   ; $LAUNCH_DIR  =~  s/\s+$// ;

my  $FITSCRIPTS_PREFIX = "SALT2mu_FITSCRIPTS" ;
my  $FITSCRIPTS_DIR;
my  $FITJOBS_DIR       = "SALT2mu_FITJOBS" ; # default fit output

my @TABLE_FORMAT_LIST    = ( "HBOOK", "ROOT" );
my @TABLE_SUFFIX_LIST    = ( "HBOOK", "ROOT" );
#my @TABLE_SUFFIX_LIST    = ( "hbook", "root" ); # removed May 25 2016
my @COMBINE_FMTARG_LIST  = ( "H"    , "R"     );
my @TABLE_FORMAT_USE ;  # logical flag for each table format 
my $COMBINE_FMTARG ;    # combine_fitres arg:  e.g., "H" ,  "R"   "H R"

my $BATCH_TEMPLATE_KICP = '$SBATCH_TEMPLATES/SBATCH_kicp.TEMPLATE' ;

# ----------------
# global inputs read from SALT2mu input file.
my ($INPUT_FILE, $input_file, @INPDIR_SNFIT_LIST, @INPDIRFIX_SNFIT_LIST );
my (@INPDIR_FITOPT,  @MUOPT_LIST, @MUOPT_LIST_ORIG, @MUOPT_TAGNAME);
my (@FITOPT_LIST,  $SUMFLAG_INPDIR, $FITOPT000_ONLY );
my ($STRINGMATCH_IGNORE,  @STRINGMATCH_IGNORE_LIST );
my ($STRINGMATCH );
my ($VERSION_EXCLUDE,  $VERSION_EXCLUDE_STRING, $PROMPT );
my ($SUBMIT_FLAG, $SUMMARY_FLAG, $NOSUBMIT_FLAG );
my (@SSH_NODELIST, $SSH_NNODE, $SNANA_LOGIN_SETUP );
my ($BATCH_COMMAND, $BATCH_TEMPLATE, $BATCH_NCORE, $DONE_STAMP_FILE );
my ($KILLJOBS_FLAG, $SUMMARY_FLAG, $OUTDIR_PREFIX, $OUTDIR_OVERRIDE );
my ($WFIT_OPT, @WFIT_INPUT, $CLEANFLAG, $NSPLITRAN ) ;

# ----------------
# misc globalas
my (@MSGERR, $N_INPDIR, $NMUOPT_SALT2mu, @NFITOPT_SNFIT);
my (@OUTDIR_SALT2mu_LIST );
my (@VERSION_SNFIT_LIST, @NVERSION_SNFIT, @VERSION_4MATCH_LIST);
my ($NVERSION_4SUM, @IVERSION_4SUM_LIST, @VERSION_4SUM_LIST);
my (@INPDIR_SDIR_LIST, $NTOT_FITRES, $NTOT_JOBS ) ;
my (@NVERSION_FINAL, @VERSION_FINAL_LIST, @SPREFIX_LIST );
my ($NROW_SUMDAT );

my ($NCPU, @NJOB_PER_CPU, $MAXJOB_PER_CPU, $icpu_MAXJOBS);
my (@CMD_PREFIX, @BATCH_NAME, @CMD_FILES);
my (@BATCH_FILES, $NOUTFILE, @CMD_CONTENTS, @NCMDLINE_PER_CPU );
my ($T_START, $T_END, $T_TOT, $NJOB_ABORT, $ALLDONE_FILE );

# ----------------
# functions

sub initStuff ;
sub parse_args ;
sub parse_inpFile ;
sub parse_MUOPT ;	
sub parse_WFIT_OPT ;
sub verify_INPDIR ;
sub makeDir_NSPLITRAN ;
sub makeDirs_SALT2mu ;
sub makeSumDir_SALT2mu ;
sub VERSION_STRINGMATCH_IGNORE ;
sub require_one_version ;
sub GET_IVERMATCH ;
sub subDirName_after_lastSlash ;
sub USE_FITOPT ;
sub copy_inpFiles ;
sub cat_inpFiles ;

sub make_COMMANDS ;
sub write_COMMANDS ;
sub prep_COMMAND ;
sub add_COMMAND ;
sub matchDump ;
sub submit_JOBS ;
sub submit_SUMMARY ;
sub wait_for_done ;
sub create_done_file ;
sub make_SUMMARY ;
sub write_SUMMARY_LOG ;
sub write_SPLITRAN_SUMMARY_LOG ;
sub write_SUMMARY_DAT ;
sub write_SUMMARY_INPDIR ;
sub gzip_logs ;
sub get_TABLE_FORMAT ;
sub CLEANFILES_DRIVER ;
sub clean_gzipSALT2mu ;

# ================== BEGIN MAIN ====================

$T_START = time() ; # gmtime();

&initStuff ;

&parse_args();

if ( $CLEANFLAG ) { &CLEANFILES_DRIVER(); exit(0); }

if ( $KILLJOBS_FLAG == 0  && $PROMPT==1 ) 
{ sntools::check_multipleJobs($SCRIPTNAME) ; } # check for existing job

&parse_inpFile();

if ( $KILLJOBS_FLAG ) {  sntools::Killjobs(@SSH_NODELIST);  die "\n"; }

# xxx if ( $WAIT_FLAG ) { goto WAIT_FOR_DONE; } 

my($IDIR,$IVER) ;

# make sure that each INPDIR exists and was created by split_and_fit
my $INPDIR ;
for($IDIR=0; $IDIR < $N_INPDIR; $IDIR++ ) 
{ &verify_INPDIR($IDIR); }

if ( $NSPLITRAN > 0 ) { 
    &makeDir_NSPLITRAN();
}
elsif ( $SUMFLAG_INPDIR == 0  ) {  
    &makeDirs_SALT2mu() ; 
    for($IDIR=0; $IDIR < $N_INPDIR; $IDIR++ ) {
	for($IVER=0; $IVER < $NVERSION_SNFIT[$IDIR]; $IVER++ ) {	
	    &copy_inpFiles($IDIR,$IVER);   # plain copy of input fitres files
	}
    }
}
else {  
    &makeSumDir_SALT2mu(); 
    my($IFITOPT);
    for($IFITOPT=0 ; $IFITOPT < $NFITOPT_SNFIT[0]; $IFITOPT++ ) {
	for($IVER=0; $IVER < $NVERSION_4SUM ; $IVER++ ) {	
	    &cat_inpFiles($IFITOPT,$IVER);   # catenate input files for sum
	}
    }
}


if ( $NVERSION_MAX < 999 ) { 
    print "\t DEBUG MODE WARNING: only $NVERSION_MAX VERSIONS processed. \n";
}


&make_COMMANDS();  # create all commands per cpu in memory
&write_COMMANDS(); # write them out into files


if ( $SUBMIT_FLAG == 0 && $SUMMARY_FLAG == 0 ) 
{ die "\n NOSUBMIT option --> DO NOT SUBMIT JOBS \n";  }


if ( $SUMMARY_FLAG == 0 ) {
    &submit_JOBS();  # submit all jobs to batch or ssh
    &submit_SUMMARY();  # launch task to wait for DONEs and make summary
    exit(0);
}

# wait for all of the SALT2mu DONE files to appear
&wait_for_done(); 

# make summary file
&make_SUMMARY();

# gzip output
&gzip_logs(); 

# create global DONE file(s) to communicate with higher level pipelines
&create_done_file() ;

exit(0);

# =================== END MAIN =====================


sub initStuff {

    $SUBMIT_FLAG=1;  $NOSUBMIT_FLAG=0;
    $PROMPT=1;  # prompt user if SALT2mu jobs already running
    # init user-inputs
    $SUMFLAG_INPDIR     = 0 ; # default is each INPDIR is independent
    $N_INPDIR           = 0 ;
    $NMUOPT_SALT2mu     = 0 ;
    $STRINGMATCH_IGNORE = "" ;
    $STRINGMATCH        = "" ;

    $BATCH_TEMPLATE    = "" ;
    $BATCH_NCORE       = 0  ;  # for qsub, sbatch, etc ...
    $BATCH_COMMAND     = "" ;
    $SSH_NNODE         = 0  ;  # for ssh
    $SNANA_LOGIN_SETUP = "" ;
    $DONE_STAMP_FILE   = "" ;

    $KILLJOBS_FLAG     = 0  ;
    $SUMMARY_FLAG      = 0 ;  # 0 => make summary after all jobs end

    $VERSION_EXCLUDE = "" ;
    $VERSION_EXCLUDE_STRING = "XQ&H!YI" ; # random garbage, but NOT blank

    # misc. inits
    $NCPU = 1 ;  # at least one for interactive
    @INPDIR_SDIR_LIST   = ();
    @INPDIR_FITOPT      = ();
    $NTOT_FITRES        = 0 ;
    $NVERSION_4SUM  = 0 ;
    
    $WFIT_OPT      = "" ;
    @WFIT_INPUT    = () ;

    $OUTDIR_OVERRIDE = "" ;
    $OUTDIR_PREFIX   = "SALT2mu" ;

    $FITOPT000_ONLY = 0 ;

    $CLEANFLAG = 0 ;

    $NSPLITRAN = -9;

}  # end init_stuff



# ==========================
sub parse_args {

    my ($NARG, $arg, $nextArg, $i, $jslash);
    
    $NARG = scalar(@ARGV);
    if ( $NARG < 1 ) {
	$MSGERR[0] = "Must give SALT2mu input filename as argument";
	sntools::FATAL_ERROR(@MSGERR);
    }

# parse name of input file. 
# INPUT_FILE includes full path; input_file has no path.
#
    $INPUT_FILE = $ARGV[0] ;
    $jslash = rindex($INPUT_FILE,"/");  # location of last slash
    if ( $jslash > 0 ) { 
	$input_file = substr($INPUT_FILE,$jslash+1, 99);
    }
    else {
	$input_file = "$INPUT_FILE" ;
	$INPUT_FILE = "$LAUNCH_DIR/$input_file" ;
    }


    if ( $ARGV[0] eq "CLEAN"  )  { $CLEANFLAG = 1 ; }

    for ( $i = 1; $i < $NARG ; $i++ ) {
	$arg = $ARGV[$i] ;

	if ( $i < $NARG-1 ) { $nextArg = $ARGV[$i+1]; }


	if ( $arg eq "KILL"    ) { $KILLJOBS_FLAG = 1 ; }
	if ( $arg eq "SUMMARY" ) { $SUMMARY_FLAG  = 1 ; $SUBMIT_FLAG=0; }
	if ( $arg eq "NOPROMPT") { $PROMPT        = 0 ; }
	if ( $arg eq "NOSUBMIT") { $SUBMIT_FLAG   = 0 ; $NOSUBMIT_FLAG=1; }

	if ( $arg eq "INPDIR+" ) {
	    $nextArg = qx(echo $nextArg); # allow for ENV
	    $nextArg =~ s/\s+$// ;        # trim trailing whitespace
	    @INPDIR_SNFIT_LIST = ( @INPDIR_SNFIT_LIST, $nextArg);
	    $SUMFLAG_INPDIR = 1;	    
	}

	if ( $arg eq "STRINGMATCH" ) 
	{ $STRINGMATCH = $nextArg ; }

	if ( $arg eq "NSPLITRAN" ) 
	{ $NSPLITRAN = $nextArg ; }

	if ( substr($arg,0,10) eq "NSPLITRAN=" )  { 
	    $NSPLITRAN = substr($arg,10,3) ; 
	    $NSPLITRAN =~ s/\s+$// ;        # trim trailing whitespace
	}

	if ( $arg eq "KICP" || $arg eq "kicp" ) 
	{ $BATCH_TEMPLATE = "$BATCH_TEMPLATE_KICP" ; }
    }


} # end pf parse_args

# ==========================
sub parse_inpFile {

# Sep 12 2019: read contents of INPFILE so that comment lines are ignored

    my ($key, @tmp, $NDIR_SOLO, $idir, $NDIR_SUM, $tmpLine, @words, $i);

    $NDIR_SUM = $NDIR_SOLO = 0 ;
    
    printf " Parse keys from : $INPUT_FILE \n"; 
 
    if ( !(-e $INPUT_FILE ) ) {
	$MSGERR[0] = "Cannot find SALT2mu input file:" ;
	$MSGERR[1] = "'$INPUT_FILE'" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

# Sep 12 2019: read contents of INPFILE so that comment lines are ignored
    my @CONTENTS_INPFILE = ();
    sntools::loadArray_fromFile($INPUT_FILE, \@CONTENTS_INPFILE);


    $key   = "DONE_STAMP:" ;
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) {
        $DONE_STAMP_FILE      = qx(echo "$tmp[0]") ;
	$DONE_STAMP_FILE      =~ s/\s+$// ;   # trim trailing whitespace
        if ( -e $DONE_STAMP_FILE )  { qx(rm $DONE_STAMP_FILE) ; }
    }


    $key = "INPDIR:";
    @tmp   = sntools::parse_array($key,2,$OPT_QUIET, @CONTENTS_INPFILE);
    $NDIR_SOLO = 0;
    foreach $tmpLine (@tmp) {	
	@words  = split(/\s+/,$tmpLine) ;
	$INPDIR_SNFIT_LIST[$NDIR_SOLO] = $words[0] ;
	$INPDIR_FITOPT[$NDIR_SOLO]     = $words[1] ; # optional FITOPT
	$NDIR_SOLO++ ;
    }


    # Nov 29 2017: check INPDIR+ key if NOT already passed via command-line.
    if ( $SUMFLAG_INPDIR == 0 ) {
	$key = "INPDIR+:";
	@tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
	$NDIR_SUM = scalar(@tmp) ;
	if ( $NDIR_SUM > 0 ) 
	{ @INPDIR_SNFIT_LIST = @tmp ; $SUMFLAG_INPDIR = 1; }
    }
    
    $key = "OUTDIR_OVERRIDE:" ;
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { 
	$OUTDIR_OVERRIDE = qx(echo $tmp[0]); # allow for ENV
	$OUTDIR_OVERRIDE =~ s/\s+$// ;   # trim trailing whitespace

	if ( index($OUTDIR_OVERRIDE,'/') < 0 ) 
	{ $OUTDIR_OVERRIDE = "$LAUNCH_DIR/$OUTDIR_OVERRIDE"; }

	print " OUTDIR_OVERRIDE: $OUTDIR_OVERRIDE \n";
    }

    $key = "OUTDIR_PREFIX:" ;
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { $OUTDIR_PREFIX = $tmp[0] ; }

    $key = "JOBNAME:";
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { 
	$JOBNAME_FIT =  qx(echo $tmp[0]); # allow for ENV
	$JOBNAME_FIT =~ s/\s+$// ;        # trim trailing whitespace
    }

    # check for NSPLITRAN key unless already given on command-line override
    if ( $NSPLITRAN < 0 ) {
	$key = "NSPLITRAN=" ;
	@tmp = qx(grep $key $INPUT_FILE);
	if ( scalar(@tmp) > 0 ) { 
	    $NSPLITRAN = substr($tmp[0],10,3); 
	    $NSPLITRAN  =~ s/\s+$// ;   # trim trailing whitespace
	}
    }

    if ( length($STRINGMATCH) == 0 ) {
	$key = "STRINGMATCH:" ;
	@tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
	if ( scalar(@tmp) > 0 ) {
	    $STRINGMATCH = $tmp[0] ; 
	}
    }

    $key = "STRINGMATCH_IGNORE:";
    @tmp   = sntools::parse_array($key,99,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { 
	$STRINGMATCH_IGNORE = $tmp[0] ; 
	@STRINGMATCH_IGNORE_LIST = split(/\s+/,$STRINGMATCH_IGNORE) ;
    }

    $key = "VERSION_EXCLUDE:";
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { $VERSION_EXCLUDE = $tmp[0] ; }

    $key = "VERSION_EXCLUDE_STRING:";
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { $VERSION_EXCLUDE_STRING = $tmp[0] ; }
    
    $key = "FITOPT000_ONLY:" ;
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { $FITOPT000_ONLY = $tmp[0] ; }

    # -------- 
    $key = "MUOPT:";
    @MUOPT_LIST_ORIG = sntools::parse_line($INPUT_FILE, 99, $key, $OPT_QUIET) ;
    @MUOPT_LIST_ORIG = ( "[DEFAULT] NONE", @MUOPT_LIST_ORIG ); # 000 is first
    $NMUOPT_SALT2mu = scalar(@MUOPT_LIST_ORIG);
    if( $NSPLITRAN > 0 ) { $NMUOPT_SALT2mu = 0; }
    if ( $NMUOPT_SALT2mu == 0 ) {
	$key = "FITOPT:";  # check legacy key, May 17 2017
	@MUOPT_LIST_ORIG = 
	    sntools::parse_array($key,99,$OPT_QUIET, @CONTENTS_INPFILE);

	@MUOPT_LIST_ORIG = ( "[DEFAULT] NONE", @MUOPT_LIST_ORIG ); 
	$NMUOPT_SALT2mu = scalar(@MUOPT_LIST_ORIG);
    }

    # check for FITOPT subsets to process (May 2017)
    for($i=0; $i < $NMUOPT_SALT2mu; $i++ ) 
    {	&parse_MUOPT($i);     }

    # ----

    $key = "WFITMUDIF_OPT:";
    @tmp   = sntools::parse_array($key,99,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { 
	$WFIT_OPT   = "$tmp[0]" ;  
	@WFIT_INPUT = ( "M0DIF" ) ; 
    }

    $key = "WFIT_OPT:";
    @tmp   = sntools::parse_array($key,99,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) { 
	$WFIT_OPT   = "$tmp[0]" ;  
	&parse_WFIT_OPT();
    }
    
    # ----------------------------------------
    # parse SSH and BATCH info

    @SSH_NODELIST = () ;
    $key = "NODELIST:" ;
    my @TMPNODES=sntools::parse_array($key,99,$OPT_QUIET, @CONTENTS_INPFILE);
    foreach $tmpLine ( @TMPNODES ) {
	@tmp = split(/\s+/,$tmpLine) ;
	@SSH_NODELIST =  ( @SSH_NODELIST , @tmp );
	$SSH_NNODE    = scalar(@SSH_NODELIST);
	$NCPU = $SSH_NNODE ;
    }

    # batch info
    $key     = "BATCH_INFO:" ;
    @tmp   = sntools::parse_array($key,3,$OPT_QUIET, @CONTENTS_INPFILE);
    if ( scalar(@tmp) > 0 ) {
	@words            = split(/\s+/,$tmp[0]) ;
	$BATCH_COMMAND  = $words[0] ;
	if (length($BATCH_TEMPLATE)==0 ) { $BATCH_TEMPLATE = $words[1] ; }
	$BATCH_NCORE    = $words[2] ;
	$NCPU = $BATCH_NCORE ;

	# allow for ENV (May 22 2017) 
        $BATCH_TEMPLATE = qx(echo $BATCH_TEMPLATE);
        $BATCH_TEMPLATE  =~ s/\s+$// ;   # trim trailing whitespace
    }

    if ( $SSH_NNODE > 0 && $BATCH_NCORE > 0 ) {
	@MSGERR[0] = "Cannot specify both NODELIST and BATCH_INFO keys.";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $SNANA_LOGIN_SETUP = "" ;
    $key = "SNANA_LOGIN_SETUP:" ;
    @tmp = sntools::parse_array($key,99,$OPT_QUIET, @CONTENTS_INPFILE);
    $SNANA_LOGIN_SETUP = "@tmp" ;

    # -----------------------------------------
    # print summary

    my ($tmpDir, @INPDIR_TMP, $iDir) ;
    $N_INPDIR   = scalar(@INPDIR_SNFIT_LIST);
    @INPDIR_TMP = @INPDIR_SNFIT_LIST ;
    $iDir = 0 ;
    foreach $tmpDir ( @INPDIR_TMP ) {

	# check for possible ENV in the directory name
	$tmpDir = qx(echo $tmpDir);
	$tmpDir =~ s/\s+$// ;   # trim trailing whitespace
	@INPDIR_SNFIT_LIST[$iDir] = "$tmpDir" ;
	$iDir++ ;
	
	if ( $SUBMIT_FLAG ) {
	    if ( $SUMFLAG_INPDIR ) 
	    { print " INPDIR+  $tmpDir \n" ; }
	    else
	    { print " INPDIR  $tmpDir \n" ; }
	}

	# abort if slash is at end . . . 
	my $jslash = rindex($tmpDir,"/");  # location of last slash
	my $len    = length($tmpDir);
	if ( $len == $jslash+1 )  {	    
	    $MSGERR[0] = "Cannot end INPDIR with slash:";
	    $MSGERR[1] = "check INPDIR: $tmpDir " ;
	    sntools::FATAL_ERROR(@MSGERR) ;
	}
    }    

    if ( $SUMMARY_FLAG ) { return ; }

    print " JOBNAME            = $JOBNAME_FIT \n" ;
    print " STRINGMATCH        = '$STRINGMATCH \n" ;
    print " STRINGMATCH_IGNORE = '$STRINGMATCH_IGNORE' \n";
    print " N(INPDIR,FITOPT)   = $N_INPDIR , $NMUOPT_SALT2mu \n";
    print " VERSION_EXCLUDE    = $VERSION_EXCLUDE \n";
    print " VERSION_EXCLUDE_STRING = $VERSION_EXCLUDE_STRING \n";
    print " WFIT_OPT           = $WFIT_OPT  (@WFIT_INPUT) \n";
    
    if ( $NSPLITRAN > 0 ) 
    { print " NSPLITRAN       = $NSPLITRAN \n"; }

    if ( $SSH_NNODE  ) 
    { print " $SSH_NNODE SSH NODES: @SSH_NODELIST \n"; }
    elsif ( $BATCH_NCORE ) 
    { print " $BATCH_NCORE batch-cores \n"; }
    else    { 
	$MSGERR[0] = "Interactive mode not supported." ;
	$MSGERR[1] = "SALT2mu input file must specifiy either";
	$MSGERR[2] = "  NODELIST:  <space-sep list of nodes>" ;
	$MSGERR[3] = "     or ";
	$MSGERR[4] = "  BATCH_INFO: [submitCmd] [templateFile] [Ncore]" ;
	sntools::FATAL_ERROR(@MSGERR);
    }


    if ( $NDIR_SUM > 0 && $NDIR_SOLO > 0 ) {
	$MSGERR[0] = "Cannot mix INPDIR and INPDIR+ keys." ;
	$MSGERR[1] = "Check '$INPUT_FILE'" ;
	$MSGERR[2] = "NDIR_SUM=$NDIR_SUM   NDIR_SOLO=$NDIR_SOLO" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    print "\n";

    return ;

} # end of parse_inpFile


# ========================
sub verify_INPDIR {
    my ($idir) = @_ ;

    my ($INPDIR, $TOPDIR, $SDIR );
    
    $INPDIR = $INPDIR_SNFIT_LIST[$idir] ;
    $INPDIRFIX_SNFIT_LIST[$idir] = "" ;
    
    # error checking:
    # - $INPDIR exists
    # - $INPDIR contains required file types
    # - $INPDIR was created by split_and_fit.pl
    # - BEWARE to minimize number of 'ls' commands, otherwise
    #    big jobs will take longer to initialize.
    
    if (! (-d $INPDIR )  ) {
	$MSGERR[0] = "INPDIR does not exist for" ;
	$MSGERR[1] = "$INPDIR" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    # check that INPDIR was created by split_and_fit
    # by checking a few required files
    my @REQUIRED_FILES = ( $MERGE_FILENAME, $FITOPT_README );
    my ($f, $F, $NFILE_MISS) ;
    $NFILE_MISS = 0 ;
    foreach $f ( @REQUIRED_FILES ) {
	$F = "$INPDIR/$f" ;
	if ( !(-e $F)  ) { $NFILE_MISS++ ; }
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - 
    # Dec 2017: check option for fixed version if INPDIR contains
    #           FITOPT*FITRES files instead of MERGE.LOG

    my @EXIST_FITRES   = ();
    if ( $NFILE_MISS > 0 ) {
	{ @EXIST_FITRES = qx(ls $INPDIR/FITOPT00*.FITRES* 2>/dev/null ); }
	if ( scalar(@EXIST_FITRES) > 0 ) {
	    ($TOPDIR,$SDIR) = &subDirName_after_lastSlash($INPDIR);
	    $INPDIR = $INPDIR_SNFIT_LIST[$idir] = $TOPDIR ;
	    $INPDIRFIX_SNFIT_LIST[$idir]        = $SDIR ;
	    $NFILE_MISS =  0;
	    print " Found FIXED INPDIR = '$SDIR' \n" ;
	}
    }

    if ( $NFILE_MISS > 0 ) {
	$MSGERR[0] = "Could not find expected files (@REQUIRED_FILES) from" ;
	$MSGERR[1] = "$INPDIR" ;
	$MSGERR[2] = "This directory is probably NOT from split_and_fit.";
	$MSGERR[3] = "INPDIR must be created by split_and_fit.pl.";
	sntools::FATAL_ERROR(@MSGERR);
    }
        
    # - - - - - - - - - - - - - - - - - - - - - - - 
    # store subdir name after last slash 
    # (used to construct ouptut SALT2mu dirs)
    ($TOPDIR,$SDIR) = &subDirName_after_lastSlash($INPDIR);
    @INPDIR_SDIR_LIST = ( @INPDIR_SDIR_LIST , $SDIR );

    # get number of split_and_fit FITOPTs from FITOPT.README
    # Note that NFITOPT_SNFIT includes 000.

    my ($key, @tmp);
    $F      = "$INPDIR/$FITOPT_README" ;
    $key    = "FITOPT:";
    @tmp    = sntools::parse_line($F, 1, $key, $OPT_ABORT) ;
    $NFITOPT_SNFIT[$idir] = scalar(@tmp); 


} # end of verify_INPDIR 

# ========================
sub parse_WFIT_OPT {

    # Examine $WFIT_OPT string and strip off M0DIF and FITRES.
    # Fill @WFIT_INPUT with M0DIF and/or FITRES

    my ($OPT_ORIG, $OPT_FINAL, @wdlist, $wd );

    $OPT_ORIG = "$WFIT_OPT" ;
    @wdlist  = split(/\s+/,$OPT_ORIG) ;

    foreach $wd (@wdlist) {
	if ( $wd eq "M0DIF" || $wd eq "FITRES" )
	{ @WFIT_INPUT = ( @WFIT_INPUT, "$wd" ); }
	else
	{ $OPT_FINAL = "$OPT_FINAL $wd" }
    }

    $WFIT_OPT = "$OPT_FINAL" ;

} # end parse_WFITOPT

# ========================
sub parse_MUOPT {

    my($imu) = @_ ;

    # Created May 27 2017
    # strip off optional FITOPT=n part of each MUOPT,
    # and store FITOPT-integer list in separate array.
    # Also strip out optional nickname in []
    #
    # MUOPT:  [ZDUM] nzbin=20 FITOPT=0
    # --> MUOPT_LIST -> nzbin=20

    my ($MUOPT_ORIG, $MUOPT, $FITOPT, @wdlist, $wd);
    my ($j0,$j1);

    if( $NSPLITRAN > 0 ) { return ; }

    $MUOPT_ORIG  = "$MUOPT_LIST_ORIG[$imu]" ;
    $MUOPT       = "" ;
    $FITOPT      = "ALL" ;
    $FITOPT_LIST[$imu] = "$FITOPT" ;  # default is all FITOPTs
    
    # -------------------------------------

    @wdlist  = split(/\s+/,$MUOPT_ORIG) ;
    foreach $wd (@wdlist) {
	if ( substr($wd,0,6) eq "FITOPT" ) {
	    $FITOPT = substr($wd,7,99);
	    $FITOPT =~ s/,/ /g ; 
	}
	else {
	    $MUOPT = "$MUOPT $wd" ;
	}
    }

    #extract optional tag name
    $j0 = index($MUOPT,'[');
    $j1 = index($MUOPT,']');
    $MUOPT_TAGNAME[$imu] = "" ;
    if ( $j1 > 0 ) { $MUOPT_TAGNAME[$imu] = substr($MUOPT,$j0+1,$j1-$j0-1); }

    $MUOPT  =~ s/\[[^)]*\]//g ;   # remove optional tagname in []      
    $MUOPT  =~ s/^\s+// ;         # remove leading spaces  
    $MUOPT_LIST[$imu]   = "$MUOPT" ;
    $FITOPT_LIST[$imu]  = "$FITOPT" ;

    if ( $SUBMIT_FLAG ) 
    { print " MUOPT[$imu] --> FITOPT= '$FITOPT' \n";  }

    return ;

}  # end parse_MUOPT


# ======================================
sub subDirName_after_lastSlash {
    my($INPDIR) = @_ ;

    # Created Dec 2017
    # Return name of subDir after last slash.
    # Example:
    #   if INPDIR = /data/college/project/astro
    #   then this function returns 
    #    ( "/data/college/project" , "astro" )
    
    my ($jslash, $TOPDIR, $SDIR);
    $jslash = rindex($INPDIR,"/");  # location of last slash
    $TOPDIR = substr($INPDIR,0, $jslash);
    $SDIR   = substr($INPDIR,$jslash+1, 99);
    return($TOPDIR,$SDIR);
    
} # subDirName_after_lastSlash

# ===================================
sub makeDir_NSPLITRAN {

    my($OUTDIR, @tmp, $datafile, $gzdatafile, $jeq );

    # - - - - - - - - - - - - - - - - - - -
    # prepare OUTDIR
    if ( length($OUTDIR_OVERRIDE)>1 ) 
    { $OUTDIR = "$OUTDIR_OVERRIDE" ; }
    else
    { $OUTDIR = "$LAUNCH_DIR/OUT_SALT2mu_NSPLITRAN${NSPLITRAN}"; }

    $OUTDIR_SALT2mu_LIST[0] = "$OUTDIR" ;

    $FITSCRIPTS_DIR = "${OUTDIR}/$FITSCRIPTS_PREFIX";

    if ( $SUMMARY_FLAG ) { return ; }

    # - - - - - below is run to submit - - - - - -

    if ( -d $OUTDIR ) { qx(rm -r $OUTDIR); }
    qx(mkdir $OUTDIR);

    # copy SALT2mu-input file into OUTDIR
    qx(cp $INPUT_FILE $OUTDIR);

    # copy FITRES file (file-argument)
    # BEWARE that key is either datafile=  or  file=
    @tmp        = qx(grep 'file=' $INPUT_FILE);
    $jeq        = index($tmp[0],'=');  # find location of equal sign
    $datafile   = substr($tmp[0],$jeq+1,500);
    $datafile   = qx(echo $datafile);
    $datafile   =~ s/\s+$// ;   # trim trailing whitespace
    $gzdatafile = "${datafile}.gz" ;

    if ( -e $datafile ) 
    { qx(cp $datafile $OUTDIR); }
    elsif ( -e $gzdatafile ) 
    { qx(cp $gzdatafile $OUTDIR); }
    else {
	$MSGERR[0] = "Unable to find required argument of" ;
	$MSGERR[1] = "datafile=$tmp[0]" ;
	$MSGERR[2] = "Either datafile=  key is missing," ;
	$MSGERR[3] = "or cannot find data file '$tmp[0]' ";
	$MSGERR[4] = "Beware that INPDIR are ignored for NSPLITRAN.";
	    
	sntools::FATAL_ERROR_STAMP($DONE_STAMP_FILE,@MSGERR);
    }

    return ;

} # end makeDir_NSPLITRAN

# ===================================
sub makeDirs_SALT2mu {

    # create separate SALT2mu_$OUTDIR dir for each INPDIR.

    my ($INPDIR, $TOPDIR, $NEWDIR, $jslash, $SDIR);
    my (@tmp, @tmpSorted, @tmpGrep, $tmpLine, @wdlist, $tmpReject, $indx );
    my ($mergeFile, $version, $NVER, $NDIR, $NVER_SKIP ) ;

    # -------------- BEGIN -------------

    $NDIR = $NVER_SKIP = 0 ;
    foreach $INPDIR ( @INPDIR_SNFIT_LIST ) { 

	if ( length($OUTDIR_OVERRIDE) > 1 ) 
	{ $TOPDIR = "$OUTDIR_OVERRIDE" ; }
	else {
	    $jslash = rindex($INPDIR,"/");  # location of last slash
	    $TOPDIR = substr($INPDIR,0, $jslash);
	}

	$SDIR   = $INPDIR_SDIR_LIST[$NDIR] ;	   
	$NEWDIR = "$TOPDIR/${OUTDIR_PREFIX}_$SDIR" ;
	$OUTDIR_SALT2mu_LIST[$NDIR] = $NEWDIR ;

	if ( $SUBMIT_FLAG  ) {
	    if ( -e $NEWDIR ) { qx(rm -r $NEWDIR) ; }
	    if($SUBMIT_FLAG) { print " create $NEWDIR \n" ; }
	    @tmp = qx(mkdir -p $NEWDIR  2>&1 );
	    if ( scalar(@tmp)>0 )  { die "\n @tmp\n ABORT \n"; }
	}
	else {
	    if ( $SUBMIT_FLAG ) { print " define $NEWDIR \n" ; }
	}

	# get list of versions and create subdir for each version
	# Check only FITOPT000 to avoid double-counting versions.
	$mergeFile = "$INPDIR/$MERGE_FILENAME";
	@tmpGrep   = qx(grep MERGED $mergeFile | grep " FITOPT000 ");
	@tmp       = ();
	foreach $tmpLine (@tmpGrep) {
	    $tmpLine  =~ s/^\s+// ;    # remove leading spaces (Sep 4 2017)
	    @wdlist   = split(/\s+/,$tmpLine) ;
	    $version  = $wdlist[1];
	    if ( $version eq "$VERSION_EXCLUDE"             ) { next ; }
	    if ( index($version,$VERSION_EXCLUDE_STRING)>=0 ) { next ; }
	    @tmp = (@tmp , $version);
	}

	@tmpSorted = sort(@tmp) ;  # -> alphabetical order
	$NVER = 0 ;
	foreach $version (@tmpSorted) {
	    if ( $NVER > $NVERSION_MAX ) { $NVER_SKIP++ ; next; }

	    # check for explicit string-match (Nov 2 2017)
	    $tmpReject = 0 ;  $indx=index($version,$STRINGMATCH);
	    if ( length($STRINGMATCH) > 2 && $indx < 0 ) { $tmpReject=1;}
	    if ( $tmpReject ) { next ; }

	    $VERSION_SNFIT_LIST[$NDIR][$NVER] = $version ;
	    $VERSION_FINAL_LIST[$NDIR][$NVER] = $version ;

	    if ( $SUBMIT_FLAG ) {
		qx(mkdir $NEWDIR/$version);
		qx(cp $INPUT_FILE $NEWDIR);
	    }
	    $NVER++ ;
	}
	if ($SUBMIT_FLAG ) { print "\t and $NVER version sub-dirs.\n"; }

	$NVERSION_SNFIT[$NDIR] = $NVER ;
	$NVERSION_FINAL[$NDIR] = $NVER ;
	$NDIR++ ;
    }
    
    print "\n";

    return ;
    
} # end of makeDirs_SALT2mu


# ========================
sub makeSumDir_SALT2mu {

    # creates one output dir in which the fitres files are summed; 
    #
    # find $TOPDIR consisting of maximum overlap of all directory names.
    # Example:
    #    INPDIR1 = /BLA/results/LOWZ
    #    INPDIR2 = /BLA/results/SDSS
    #    INPDIR3 = /BLA/results/SNLS
    #
    #  Then TOPDIR = /BLA/results   and the SALT2mu directory 
    #  created is  /BLA/results/SALT2mu_LOW+SDSS+SNLS
    #
    # Oct 2 2013: see new tmpGrep logic to pick out FITOPT000
    # Dec 22 2017: check INPDIRFIX_SNFIT_LIST for fixed version
    # Mar 29 2019: for single directory, stop string match after
    #              entire string is checked.
    # May 24 2019: allow user-defined OUTDIR_OVERRIDE
    #
    # Mar 10 2020: check require_one_version().
    #
    # -------------------------

    my ( $LEND, $MATCH, $i, $c1ref, $c1, $LENMATCH ) ;
    my ( $TOPDIR, $INPDIR, @tmp, $NVER_SKIP, $tmpReject, $DIRFIX );
    my ( $jslash, $jstart_SDIR, $indx, $OUTDIR);
    
    # -------------------- BEGIN ----------------------


    # - - - - - - - - - 
    # find maximum substring of INPDIR that matches all INPDIRs
    $LEND = length($INPDIR_SNFIT_LIST[0]) ;
    $MATCH = 1;  $i=0 ;

    while ( $MATCH && $i < $LEND ) {
	$c1ref = substr($INPDIR_SNFIT_LIST[0],$i,1);
	foreach $INPDIR (@INPDIR_SNFIT_LIST) {
	    $c1 = substr($INPDIR,$i,1);
	    if ( $c1 ne $c1ref ) { $MATCH = 0 ; $LENMATCH = $i ; }
	}
	$i++ ; 
    }
    
    $TOPDIR = substr($INPDIR_SNFIT_LIST[0], 0, $LENMATCH-1);
    $jstart_SDIR = 0;

    # Mar 12 2017: if $TOPDIR does not exist, go back to previous slash
    if ( !(-d $TOPDIR ) ) {
	$jslash = rindex($TOPDIR,"/");      # location of last slash
	$TOPDIR = substr($TOPDIR,0,$jslash);
	$jstart_SDIR = $LENMATCH - $jslash-1 ;
    }

    # ------------------------------------------------
    # construct SALT2mu subdir name from INPDIR_SNFIT names
    my ($SDIR, $SDIR_SUM, $NDIR );
    $NDIR = 0 ;
    $SDIR_SUM = "$OUTDIR_PREFIX" ; 
    foreach $INPDIR (@INPDIR_SNFIT_LIST) {
	$SDIR = substr($INPDIR_SDIR_LIST[$NDIR],$jstart_SDIR,99);
	if ( $NDIR == 0 ) { $SDIR_SUM = "${SDIR_SUM}_${SDIR}" ; }
	else              { $SDIR_SUM = "${SDIR_SUM}+${SDIR}" ; }
	$NDIR++ ;
    }

    if ( length($OUTDIR_OVERRIDE)>1 ) 
    { $OUTDIR = $OUTDIR_OVERRIDE;   }
    else 
    { $OUTDIR  = "$TOPDIR/$SDIR_SUM" ; }

    # - - - - - - - - - - - - - - - - - - -     

 MAKEDIR:

    $OUTDIR_SALT2mu_LIST[0] = "$OUTDIR" ;
    if ( $SUBMIT_FLAG || $NOSUBMIT_FLAG ) {
	if ( -d $OUTDIR ) { qx(rm -r $OUTDIR) ; }
	print " Create $OUTDIR \n";
	@tmp = qx(mkdir -p $OUTDIR  2>&1 );
	if ( scalar(@tmp) ) { die "\n @tmp \n ABORT "; }	
	qx(cp $INPUT_FILE $OUTDIR);
    }
    else {
	print " Define $OUTDIR \n";
    }


    # --------------------------------------------
    # now the tricky part; find matching subdirs to sum.
    # Ignore subdirs without a match.
    # If there is more than 1 match, abort.   

    print " \n";
    my ($mergeFile, $NVER, $vers, $vers4match, @VTMPLIST, @tmp, @tmpSorted );
    my (@wdlist, $version, @tmpGrep, $tmpLine );

    # start by getting list of versions with STRINGMATCH_IGNORE removed;
    # see VERSION_4MATCH_LIST that are used to match subdirs.

    $NDIR = $NVER_SKIP = 0;
    foreach $INPDIR (@INPDIR_SNFIT_LIST) {

	$mergeFile = "$INPDIR/$MERGE_FILENAME";  
	@tmpGrep   = qx(grep MERGED $mergeFile | grep " FITOPT000 ");
	@tmp       = ();

	if ( scalar(@tmpGrep) == 0 ) {
	    $MSGERR[0] = "Unable to grep 'MERGED' and ' FITOPT000 '";
	    $MSGERR[1] = "from $mergeFile";
	    $MSGERR[2] = "Check for aborts or tabs." ;
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP_FILE,@MSGERR);
	}
	foreach $tmpLine (@tmpGrep) {
	    $tmpLine  =~ s/^\s+// ;    # remove leading spaces  (Sep 4 2017)
	    @wdlist  = split(/\s+/,$tmpLine) ;
	    $version = $wdlist[1];
	    $indx    = index($version,$VERSION_EXCLUDE_STRING) ;
	    if ( $version eq "$VERSION_EXCLUDE"    ) { next ; }
	    if ( $indx >=0                         ) { next ; }

	    # check for explicit string-match (Nov 2 2017)
	    $tmpReject = 0 ;  $indx=index($version,$STRINGMATCH);
	    if ( length($STRINGMATCH) > 2 && $indx < 0 ) { $tmpReject=1;}
	    if ( $tmpReject ) { next ; }

	    @tmp = (@tmp , $version);
	}

	# if one version is required abort on multiple versions
	if ( &require_one_version() == 1 ) {
	    $NVER = scalar(@tmp) ;
	    if ( $NVER != 1 ) {
		$MSGERR[0] = "Found $NVER versions in" ;
		$MSGERR[1] = "$INPDIR/MERGE.LOG ." ;
		$MSGERR[2] = "But ... only 1 version allowed for option to ";
		$MSGERR[3] = "ignore string matches.";
		$i = 4;
		foreach $vers ( @tmp ) 
		{ $MSGERR[$i] = "Found VERSION = $vers" ; $i++; }
		sntools::FATAL_ERROR_STAMP($DONE_STAMP_FILE,@MSGERR);  
	    }
	} 

	@tmpSorted = sort(@tmp) ;  # -> alphabetical order
	$NVER = 0 ;
	$DIRFIX  = "$INPDIRFIX_SNFIT_LIST[$NDIR]" ;

	foreach $vers (@tmpSorted) {
	    if ( $NVER > $NVERSION_MAX ) { $NVER_SKIP++; next; }
	    $vers4match = &VERSION_STRINGMATCH_IGNORE($vers);

            if ( length($STRINGMATCH) > 2 )  { $vers4match = $STRINGMATCH ; }
	    if ( length($DIRFIX)      > 2 )  { $vers = $DIRFIX ; }	   

	    $VERSION_SNFIT_LIST[$NDIR][$NVER]  = $vers ;
#	    print " xxx load VERSION_SNFIT_LIST[$NDIR][$NVER] = $vers\n";
	    $VERSION_4MATCH_LIST[$NDIR][$NVER] = $vers4match ;
#	    print " xxx $vers -> '$vers4match' \n";
	    $NVER++ ;
	}

	$NVERSION_SNFIT[$NDIR] = $NVER ;
	$NDIR++ ;
    }


    if ( $NVER_SKIP > 0 ) {
	$MSGERR[0] = "Skipped $NVER_SKIP versions";
	$MSGERR[1] = "check $NVERSION_MAX = $NVERSION_MAX" ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP_FILE,@MSGERR);
    }

    # -----------

    my ($idir, $iver, $iverMatch, $NVER_REF, $VER_REF);
    my ($NDIR_MATCH );

    $NDIR     = $N_INPDIR ;
    $NVER_REF = $NVERSION_SNFIT[0];
    $NVERSION_4SUM = 0 ;

    for ($iver=0; $iver < $NVER_REF; $iver++ ) {
	$VER_REF    = $VERSION_4MATCH_LIST[0][$iver] ;
	$NDIR_MATCH = 1 ;  # always matches itself

	for($idir=1; $idir < $NDIR; $idir++ ) {
	    if ( &require_one_version() == 1 )
	    { $iverMatch = 0 ; }
	    else
	    { $iverMatch = &GET_IVERMATCH($VER_REF,$idir); }
	    if ( $iverMatch >=0 ) {
		$NDIR_MATCH++ ;
		$IVERSION_4SUM_LIST[$idir][$NVERSION_4SUM] = $iverMatch ;
	    }
	} # idir

	if ( $NDIR_MATCH == $NDIR ) {
	    $IVERSION_4SUM_LIST[0][$NVERSION_4SUM]  = $iver ;
	    $VERSION_4SUM_LIST[$NVERSION_4SUM]      = $VER_REF ;
	    $VERSION_FINAL_LIST[0][$NVERSION_4SUM]  = $VER_REF ;
	    $NVERSION_4SUM++ ;
	    $NVERSION_FINAL[0] = $NVERSION_4SUM ;
#	    print "\t Found match for FITRES-sum: $VER_REF \n"; 

	    if ( $SUBMIT_FLAG || $NOSUBMIT_FLAG ) 
	    { qx(mkdir $OUTDIR/$VER_REF); }
	}
	else { 
#	    print "\t Cound not match $VER_REF (NDIR_MATCH=$NDIR_MATCH)\n"; 
	}

    }  # iver

    print "\n" ;
    return ;

} # end of makeSumDir_SALT2mu 

# ================================
sub GET_IVERMATCH {
    my ($VERSION,$idir) = @_ ;

    # for input string $VERSION, check all versions in INPDIR
    # with index $idir, and look for  a string-match. 
    # Return $iverMatch of matching version.
    # Abort if more than 1 version matches.

    my ($iverMatch, $VMATCH,  $NMATCH, $iver, $NVER);

    $iverMatch = -9 ;
    $NMATCH    =  0 ;
    $NVER      = $NVERSION_SNFIT[$idir] ;  # number of versions for this idir

    for($iver=0; $iver < $NVER; $iver++ ) {
	$VMATCH = "$VERSION_4MATCH_LIST[$idir][$iver]" ;
	if ( $VERSION eq $VERSION_4MATCH_LIST[$idir][$iver]  ) {
	    $iverMatch = $iver ;
	    $NMATCH++ ;
#	    print " xxx iver=$iver  VMATCH = $VMATCH \n" ;
	}
    }

    if ( $NMATCH > 1 ) {
	$MSGERR[0] = "NMATCH = $NMATCH for ";
	$MSGERR[1] = "VERSION=$VERSION and idir=$idir" ;
	$MSGERR[2] = "(idir -> $INPDIR_SNFIT_LIST[$idir])" ;
	$MSGERR[3] = "No more than 1 match allowed.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP_FILE,@MSGERR);
    }

    return $iverMatch ;

} # end of GET_IVERMATCH

# ============================
sub VERSION_STRINGMATCH_IGNORE {
    my ($version) = @_ ;

    my ($version_out, $wd, $WD, $firstChar, $lastChar);
    if ( length($STRINGMATCH_IGNORE) == 0 ) { return $version ; }

    # return $version with ignore-strings removed
    $version_out = $version ;
    
    if ( &require_one_version() == 1 ) { return($version_out); }

    foreach $wd ( @STRINGMATCH_IGNORE_LIST ) {
	$WD  = "$wd" ;
        $WD  =~ s/\+/\\\+/ ;  # add backslash in front of special char + sign 
        $version_out =~ s/$WD// ;
    }

    # Jan 15 2017: remove leading/trailing underscore or dash
    $firstChar = substr($version_out,0,1);
    if ( $firstChar eq "_" || $firstChar eq "-" )
    { $version_out = substr($version_out,1,99); }

    $lastChar  = substr($version_out,-1);
    if ( $lastChar eq "_" || $lastChar eq "-" )
    { chop($version_out); }
    
#    printf " xxx STRINGMATCH: $version -> $version_out \n";
    
    if ( $version_out eq '' ) { $version_out = "$FITJOBS_DIR" ; }

    return $version_out ;

} # end of VERSION_STRINGMATCH_IGNORE {


# ==================================
sub require_one_version {

    # Mar 10 2010
    # check option to ignore string matches, but this
    # works only of there is one version per directory.

    my($istat, $istat0, $istat1);
    $istat0 =  ( "$STRINGMATCH_IGNORE" eq "IGNORE" ) ;
    $istat1 =  ( "$STRINGMATCH"        eq "IGNORE" ) ;
    $istat  =  ( $istat0 || $istat1 ) ;

#    print " xxx require: stat = $istat0, $istat1, $istat \n";

    return($istat);

} # end require_one_version

# ==================================
sub copy_inpFiles {

    my ($idir,$iver) = @_ ;

    # copy FITOPT[nnn].FITRES from each INPDIR
    # to corresponding $OUTDIR.
    #
    # Apr 2017: If FITIOT[nnn].FITRES.gz exists, gunzip it.
    
    my ($VERSION, $INPDIR, $OUTDIR, $TOPDIR );
    my (@cpList, $cpString, @tmp, $f,$F,  $fgzip ) ;

    if ( $SUMMARY_FLAG ) { return ; }

    $INPDIR   = $INPDIR_SNFIT_LIST[$idir] ;
    $VERSION  = $VERSION_SNFIT_LIST[$idir][$iver];

    $OUTDIR   = "$OUTDIR_SALT2mu_LIST[$idir]/$VERSION";

    if ( $iver == 0 ) {
	print "# ---------------------------------------------- \n";
	print " Copy input FITRES files from \n\t $INPDIR \n";
	print "     to  $OUTDIR_SALT2mu_LIST[$idir] \n";
    } 

    print "$iver " ;
    if ( $iver == $NVERSION_SNFIT[$idir] - 1 ) { print "\n"; }

    $TOPDIR   = "$INPDIR/$VERSION" ;
    $cpString = "FITOPT*.FITRES";    

    # get exact cpList of files to copy
    @cpList = () ;
    @tmp = qx(cd $TOPDIR; ls ${cpString}* );
    foreach $f (@tmp) {
	$f     =~  s/\s+$// ;

	if ( index($f,"ORIG") > 0 ) { next; }
	
	@cpList = ( @cpList , $f );
	$NTOT_FITRES++ ;
    }

    # do the copy
    qx(cd $TOPDIR ; cp @cpList  $OUTDIR);

}  # end of copy_inpFiles


# ==================================
sub cat_inpFiles {

    # May 2020
    # catenate FITOPT[nnn].FITRES from each INPDIR
    # to make a summed fitres file.
    #
    # Use cat_snana_table utility which allows different
    # columns in each fitres file.
    #
    # Note that input $IVER is a sparse verson index for the
    # versions to sum.
    #
    # Compared to cat_inpFiles_legacy, this function
    #  + allows different columns
    #  + allows mix of gzipped and unzipped FITRES files

    my ($ifitopt,$iver_sum) = @_ ;

    my ($VERSION_4SUM, $VERSION_SNFIT, $INPDIR, $Ngzip, $Nunzip, $cdout);
    my ($iver, $idir, $nnn, $ffile, $FFILE);
    my ($CAT_DIR, $CATLIST, $CATFILE_OUT, $CMD_CAT, $CMD_GZIP, $CAT_LOG );

    if ( $SUMMARY_FLAG ) { return ; }
    if ( $FITOPT000_ONLY && $ifitopt > 0 ) { return ; }

    $VERSION_4SUM = $VERSION_4SUM_LIST[$iver_sum] ;

    # get comma-sep CATLIST of original VERSION-files to catenate

    $CATLIST  = "" ;
    $nnn      = sprintf("%3.3d", $ifitopt);
    $ffile    = "FITOPT${nnn}.FITRES" ;
    $Ngzip = $Nunzip  = 0 ;
    
    print "# -------------------------------------------------- \n";
    print " Catenate input files for:  $VERSION_4SUM  $ffile\n";
   
    for($idir=0; $idir < $N_INPDIR; $idir++ ) {
	$iver          = $IVERSION_4SUM_LIST[$idir][$iver_sum] ;
	$VERSION_SNFIT = $VERSION_SNFIT_LIST[$idir][$iver];
	$INPDIR        = $INPDIR_SNFIT_LIST[$idir] ;
	
	print "\t Add VERSION[$iver] = $VERSION_SNFIT \n";

	$FFILE  = "${INPDIR}/${VERSION_SNFIT}/$ffile" ;

	if ( -l $FFILE ) {
	    my $symLink = readlink($FFILE);
	    my $jslash = rindex($symLink,"/");  # location of last slash
	    if ( $jslash < 0 ) 
	    { $FFILE      = "${INPDIR}/${VERSION_SNFIT}/$symLink" ; }
	    else
	    { $FFILE = $symLink; }
	}
	if ( -e "$FFILE.gz" ) 	{ $Ngzip++; }
	else                    { $Nunzip++ ; }

	if ( $idir == 0 )    { $CATLIST = "$FFILE" ; }
	else                 { $CATLIST = "$CATLIST,$FFILE" ; }
	
    }  # idir
 

    $CAT_DIR     = "$OUTDIR_SALT2mu_LIST[0]/${VERSION_4SUM}" ;
    $cdout       = "cd $CAT_DIR" ;
    $CATFILE_OUT = "$CAT_DIR/$ffile" ;
    $CAT_LOG     = "TMP_CAT.LOG" ;

    $CMD_CAT  = "SALT2mu.exe cat_only  ";
    $CMD_CAT .= "datafile=$CATLIST ";
    $CMD_CAT .= "append_varname_missing='PROB*' " ;
    $CMD_CAT .= "catfile_out=$CATFILE_OUT ";

    $CMD_GZIP = "" ;
    $CMD_GZIP = "gzip $CATFILE_OUT " ;
    qx($CMD_CAT > $CAT_LOG ; rm $CAT_LOG );
    if ( $Ngzip > 0 ) { qx($CMD_GZIP); }

#    die "\xxx DEBUG DIE xxx\n";
    $NTOT_FITRES++ ;

    return ;
#    print " xxx CATLIST = $CATLIST \n";
#    print " xxx \t --> $CATFILE \n";

} # end of cat_inpFiles


# ==============================
sub make_COMMANDS {

    my ($icpu, $inode, $node, $nnn, $cmdFile, $suffix);
    my ($jdot, $prefix, $SCRIPT, $PREFIX, $NJOB, $SUBDIR );
    my ( $OUTDIR, $NDIR);

    # construct name of FITSCRIPTS directory using prefix if $INPUT_FILE
    $jdot   = index($input_file,".");
    $prefix = substr($input_file,0,$jdot);
    $SUBDIR = "${FITSCRIPTS_PREFIX}_${prefix}" ;
    $FITSCRIPTS_DIR = "$LAUNCH_DIR/$SUBDIR";
    
    if ( length($OUTDIR_OVERRIDE) > 1 ) 
    { $FITSCRIPTS_DIR = "${OUTDIR_OVERRIDE}/${FITSCRIPTS_PREFIX}"; }

    # - - - - - - - - - - - - - - - - - - -

    $NDIR = scalar(@OUTDIR_SALT2mu_LIST) ;
    if ( $NDIR == 1 ) {
	$OUTDIR = $OUTDIR_SALT2mu_LIST[0] ;
	$FITSCRIPTS_DIR = "${OUTDIR}/${FITSCRIPTS_PREFIX}" ; 
    }

    if ( !$SUMMARY_FLAG ) {
	print "\n Create command-scripts in \n\t $FITSCRIPTS_DIR \n";
	if ( -d $FITSCRIPTS_DIR ) { qx(rm -r $FITSCRIPTS_DIR); }
	qx(mkdir $FITSCRIPTS_DIR) ;
    }

    # ----------------------------------------------
    # first create command script for each node/core
    # Store script names in @CMD_FILES

    $inode = 0 ; $MAXJOB_PER_CPU = $icpu_MAXJOBS = -9 ;
    for($icpu = 0; $icpu < $NCPU; $icpu++ ) {
	$nnn               = sprintf("%3.3d", $icpu);
	
	if ( $SSH_NNODE  ) 
	{ $suffix = "CPU${nnn}_$SSH_NODELIST[$icpu]" ;  }
	elsif ( $BATCH_NCORE ) 
	{ $suffix = "CPU${nnn}" ; }  
	else
	{ $suffix = "Interactive" ; }  

	$PREFIX            =  "${FITSCRIPTS_PREFIX}_${suffix}" ;
	$CMD_PREFIX[$icpu] =  "$PREFIX" ;
	$cmdFile           =  "${PREFIX}.CMD" ;
	$CMD_FILES[$icpu]  =  "${cmdFile}" ;

	$BATCH_NAME[$icpu] =  "${input_file}-${suffix}";

	if ( $SUBMIT_FLAG )  { qx(touch $FITSCRIPTS_DIR/$cmdFile ); }
#	print "\t -> $cmdFile \n";

	$NJOB_PER_CPU[$icpu] = 0 ;
	$NCMDLINE_PER_CPU[$icpu] = 0 ;
		
    }  # icpu

    # Loop over OUTDIRs and subdirs.
    # Create SALT2mu job for every FITRES file in every dir.
    # Unitarity check is made at the end.

    if ( $NSPLITRAN > 0 ) {  
	my ($isplit);
	$NTOT_JOBS = $NSPLITRAN + 1;  # add 1 for summary job
	for($isplit=1; $isplit<= $NTOT_JOBS; $isplit++ )
	{ &NSPLITRAN_prep_COMMAND($isplit); }
	goto CHECK_NTOT_JOBS ;
    }

    my ($NDIR, $NVER, $idir, $iver, $ifitopt );
    $NTOT_JOBS = 0 ;    $icpu = 0 ;
    $NDIR = scalar(@OUTDIR_SALT2mu_LIST) ;

    for($idir=0 ; $idir < $NDIR; $idir++ ) {
	$NVER   = $NVERSION_FINAL[$idir] ;
	for ( $iver=0; $iver  < $NVER; $iver++ ) {
	    for($ifitopt=0; $ifitopt < $NFITOPT_SNFIT[$idir]; $ifitopt++ ) {
		&prep_COMMAND($idir,$iver,$ifitopt);
	    }  # ifitopt
	}    # iver
    }  # idir


    # ---------------------------
    # If there are more CPUs then NTOT_JOBS, then reduce NCPU
    # and remove unused CMD files. This is to avoid confusion
    # with extra jobs that don't do anything but leave garbage
    # log/done files.

  CHECK_NTOT_JOBS:
    
    if ( $NTOT_JOBS == 0 && !$SUMMARY_FLAG ) {
	print "\n PRE-ABORT DUMP: \n";
	&matchDump();		    
	$MSGERR[0] = "No SALT2mu jobs found in make_COMMANDS() .";
	$MSGERR[1] = "Probably a stringMatch problem";
	$MSGERR[2] = "NSPLITRAN=$NSPLITRAN  NTOT_JOBS=$NTOT_JOBS" ;
	$MSGERR[3] = "SUBMIT_FLAG=$SUBMIT_FLAG  SUMMARY_FLAG=$SUMMARY_FLAG";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP_FILE,@MSGERR);
    }
    

    if ( $NCPU > $NTOT_JOBS ) {	
	my $NCPU_ORIG = $NCPU ;
	$NCPU = $NTOT_JOBS ;

	if ( $SUMMARY_FLAG )  { return ; }

	print "\n\t";
	print "$NCPU_ORIG CPUs is too many for $NTOT_JOBS jobs: ";
	print "reduce NCPU -> $NTOT_JOBS \n" ;

	# removed unused CMD files.
	for($icpu = $NCPU; $icpu < $NCPU_ORIG; $icpu++ ) {
	    my $cmd_file = "$FITSCRIPTS_DIR/$CMD_FILES[$icpu]" ;
	    if ( -e $cmd_file ) { qx(rm $cmd_file); }
	}
    }

    if ( $SUMMARY_FLAG  )  { return ; }

    # -----------------------------

    for($icpu = 0 ; $icpu < $NCPU; $icpu++ ) {

	$NJOB    = $NJOB_PER_CPU[$icpu] ;
	$PREFIX  = "$CMD_PREFIX[$icpu]" ;

	if ( $BATCH_NCORE == 0 ) {    # batch & interactive
	    # nothing more to do 
	    print "\t Created $cmdFile  with $NJOB jobs. \n";
	}
	else {
	    # batch mode: create BATCH file for each CMD_FILE

	    my ( $script, $batchName, $batchFile, $batchLog, $JOB, $doneFile);

	    $script    = $CMD_FILES[$icpu] ;
# xxx mark delete Oct 13 2019  $batchName = "${PREFIX}" ;
	    $batchName = "$BATCH_NAME[$icpu]";
	    $batchFile = "${PREFIX}.BATCH" ;
	    $batchLog  = "${PREFIX}.LOG" ;
	    $doneFile  = "$FITSCRIPTS_DIR/${PREFIX}.DONE" ;
	    $JOB       = "source $script" ;
	    $BATCH_FILES[$icpu]  =  "${batchFile}" ;

	    print "\t prepare $batchFile  for $BATCH_COMMAND  ($NJOB jobs)\n";
	   
	    sntools::make_batchFile($BATCH_TEMPLATE, $FITSCRIPTS_DIR,
				    $batchName, $batchFile, $batchLog, 
				    $JOBMEMORY, $JOB);
	}
    }  # icpu

    return ;

}  # end of make_COMMANDS


# ========================================
sub NSPLITRAN_prep_COMMAND {
    my ($isplit) = @_ ;

    # May 2019
    # Analog of prep_COMMAND, but here prepare command for
    # each of the individual SPLITRAN jobs using JOBID_SPLITRAN 
    # argument to SALT2mu.exe.
    #
    # Jun 11 2020: include optional wfit commands
    # Jun 24 2020: 
    #  + refactor to make subDir for each splitRan job
    #  + see LEGACY flag to run old-style with all files in OUTDIR
    
    my ($icpu, $OUTDIR, $OUTDIR_SPLITRAN, $CMD, $LOGFILE, $ARGLIST);
    my ($OUT_PREFIX, $OUT_PREFIX_LEGACY, $SUBDIR );
    my $ISJOB_SUMMARY = ( $isplit > $NSPLITRAN );
    my $LEGACY_FLAG = 0 ;  # 1 -> all files in OUTDIR; 0-> subDirs

    # pick CPU on round-robin basis
    $icpu = (($isplit-1) % $NCPU) ;   # 0 to NCPU-1  

    $OUTDIR     = "$OUTDIR_SALT2mu_LIST[0]" ;

    if ( $LEGACY_FLAG ) { 
	$OUT_PREFIX = "OUT_TEST" ;
	$OUT_PREFIX_LEGACY = sprintf("%s-SPLIT%3.3d", $OUT_PREFIX, $isplit);
	$SUBDIR = "./" ;
	$OUTDIR_SPLITRAN = $OUTDIR ;
	$LOGFILE = "${OUT_PREFIX_LEGACY}.LOG" ;
    }
    else {
	$OUT_PREFIX = "SALT2mu_FITOPT000_MUOPT000" ; # 6.24.2020
	$LOGFILE    = "${OUT_PREFIX}.LOG" ;
	$SUBDIR     = sprintf("SPLITRAN-%4.4d", $isplit);
	$OUTDIR_SPLITRAN = "$OUTDIR/$SUBDIR" ; 
	if ( !$SUMMARY_FLAG && !$ISJOB_SUMMARY ) 
	{ qx(mkdir $OUTDIR_SPLITRAN) ;  }
    }


    # put summary job in last CPU so it's likely to finish last
    if ( $ISJOB_SUMMARY ) {
	$icpu    = $icpu_MAXJOBS ; 
	$LOGFILE = "SALT2mu_SPLITRAN_SUMMARY.LOG";
#	$LOGFILE = "${OUT_PREFIX}_summary.log";
    }

    $NJOB_PER_CPU[$icpu]++ ;        # and NJOB for this CPU 

    # keep track of icpu with max number of jobs; to use for summary job.
    if ( $NJOB_PER_CPU[$icpu] >= $MAXJOB_PER_CPU ) {
	$MAXJOB_PER_CPU = $NJOB_PER_CPU[$icpu];
	$icpu_MAXJOBS   = $icpu ;
    }

    $ARGLIST = 
	"NSPLITRAN=$NSPLITRAN  " . 
	"JOBID_SPLITRAN=$isplit  " .  
	"prefix=$OUT_PREFIX" ;
    $CMD     = "$JOBNAME_FIT $INPUT_FILE $ARGLIST \\" ;

    # cd just once
    if ( $NJOB_PER_CPU[$icpu] == 1 ) {
	&add_COMMAND($icpu, "cd $OUTDIR_SPLITRAN" ) ;
    }

    &add_COMMAND($icpu, "" ) ;
    &add_COMMAND($icpu, "# ------------------------------------" ) ;
    if ( $isplit > $NSPLITRAN ) { 
	&add_COMMAND($icpu, "sleep 10"); 
	&add_COMMAND($icpu, "cd $OUTDIR" ) ; # 6.24.2020
    }
    &add_COMMAND($icpu, "$CMD");
    &add_COMMAND($icpu, "    >& $LOGFILE " );

    # Jun 11 2020: check for wfit
    my($tmpInput, $inFile_tmp, $wPREFIX, $CMDwfit );
    foreach $tmpInput ( @WFIT_INPUT ) {  # M0DIF or FITRES

	# replace "SALT2mu" with "wfit" in prefix
	$wPREFIX    = "wfit_${OUT_PREFIX}" ; 	
	$inFile_tmp = "${OUT_PREFIX}.${tmpInput}" ;
	$CMDwfit = 
	    " $JOBNAME_WFIT $inFile_tmp $WFIT_OPT " .
	    "-cospar ${wPREFIX}.COSPAR " .
	    "-resid  ${wPREFIX}.RESID " .
	    "   >&   ${wPREFIX}.LOG " ; 
	&add_COMMAND($icpu, "" ) ;
	&add_COMMAND($icpu, $CMDwfit ) ;
    }

    return ;

} # end NSPLITRAN_prep_COMMAND

# ========================================
sub prep_COMMAND {

    # write SALT2mu commands for input dir, version , fitopt.
    # ifitopt refers to the split_and_fit options.
    # Increment NTOT_JOBS here.
    #
    # Jun 17 2016: remove .text file created by combine_fitres.exe
    # Jul 15 2016: incude INPDIR_FITOPT in OPT_S2MU
    #
    # May 17 2017: 
    #   + call USE_FITOPT to check option to skip
    #   + compute $icpu locally instead of passed as argument
    #   + store commands in memory via add_COMMAND() call
    #     Disk-write is later in write_COMMANDS().
    #
    # Jun 25 2017: pass -label arg to wfit
    # Dec 02 2017: do not create HBOOK table
    # Dec 21 2017: same prefix format with or without MUOPT
    #
    # ---------- BEGIN ----------

    my ($idir,$iver,$ifitopt) = @_ ;

    my ($OUTDIR, $VERSION, $OPT_S2MU, $IOPT_S2MU, $ONEOPT, $MUOPT );
    my ($CMD, $CMDwfit, $FFF, $FPREFIX, $SPREFIX );
    my ($wPREFIX, $tmpInput );
    my ($sss, $LOGFILE, $FARG, $SARG, $icpu );
    my ($out_fitres, $out_FITRES, $JOBNAME );
    my ($argOut, $COMBINE_ARGS, $combLog, $combText, $DO_SNANA_SETUP );

    $OUTDIR  = $OUTDIR_SALT2mu_LIST[$idir] ;
    $VERSION = $VERSION_FINAL_LIST[$idir][$iver] ;

    # name of input fitres file
    $FPREFIX = sprintf("FITOPT%3.3d", $ifitopt);
    $FFF     = "$FPREFIX.FITRES" ;

    if ( $idir==0 && $iver==0 && $ifitopt==0 ) { &get_TABLE_FORMAT(); }

    $ONEOPT = ( $NMUOPT_SALT2mu == 1 ) ;
    for ( $IOPT_S2MU=0; $IOPT_S2MU < $NMUOPT_SALT2mu; $IOPT_S2MU++ ) {

	$MUOPT = "" ;  # for imu=0
	if ( $IOPT_S2MU > 0 ) { $MUOPT = "$MUOPT_LIST[$IOPT_S2MU]" ; }

	$OPT_S2MU = "$INPDIR_FITOPT[$idir] $MUOPT" ;
	$sss = sprintf("%3.3d", $IOPT_S2MU );

	$SPREFIX = "SALT2mu_${FPREFIX}_MUOPT${sss}" ;
	$SPREFIX_LIST[$idir][$iver][$ifitopt][$IOPT_S2MU] = $SPREFIX ;

	# check option to skip this MUOPT/FITOPT (May 2017)
	if ( &USE_FITOPT($IOPT_S2MU,$ifitopt) == 0 ) { next; }

	# pick CPU on round-robin basis
	$icpu = ($NTOT_JOBS % $NCPU) ;  # 0 to NCPU-1
	$NTOT_JOBS++ ;                  # increment total number of jobs
	$NJOB_PER_CPU[$icpu]++ ;        # and NJOB for this CPU

	if ( $SUMMARY_FLAG  ) { next ; }

	$FARG    = "datafile=$FFF";
	$SARG    = "prefix=$SPREFIX" ;
	$LOGFILE = "${SPREFIX}.LOG"  ;
	$JOBNAME = "$JOBNAME_FIT" ;
	$CMD     = "$JOBNAME $INPUT_FILE  $FARG $SARG $OPT_S2MU \\" ;

	$out_fitres = "${SPREFIX}.${OUTPUT_SUFFIX_fitres}" ;
	$out_FITRES = "${SPREFIX}.FITRES";

	# create formatted table from fitres text file
	# Note that out_FITRES is the SALT2mu output and the
	# combine_fitres.exe input.
	    
	&add_COMMAND($icpu, "" ) ;
	&add_COMMAND($icpu, "# ------------------------------------" ) ;
	&add_COMMAND($icpu, "cd $OUTDIR/${VERSION}" ) ;
	&add_COMMAND($icpu, "$CMD");
	&add_COMMAND($icpu, "    >& $LOGFILE " );
	&add_COMMAND($icpu, "mv ${out_fitres} ${out_FITRES}");

	my $MAKETABLE_HBOOK = 0;
	if ( $MAKETABLE_HBOOK ) {
	    $JOBNAME   = "$JOBNAME_MAKETABLE" ;
	    $argOut    = "-outprefix $SPREFIX" ;
	    $COMBINE_ARGS   = "${out_FITRES} $COMBINE_FMTARG $argOut" ;
	    $combLog        = "combine_${SPREFIX}.LOG" ;
	    $combText       = "${SPREFIX}.text" ;  # Jun 17 2016

	    &add_COMMAND($icpu, "$JOBNAME $COMBINE_ARGS > $combLog");
	    &add_COMMAND($icpu, "rm ${combLog} ${combText}") ;
	}

	# Jun 3 2016: check option to run wfit 
	foreach $tmpInput ( @WFIT_INPUT ) {  # M0DIF or FITRES

	    # replace "SALT2mu" with "wfit" in prefix
	    my $inFile_tmp = "${SPREFIX}.${tmpInput}" ;
	    if ( $ONEOPT ) 
	    { $wPREFIX  = "wfit_${tmpInput}_${FPREFIX}" ;  }	    
	    else
	    { $wPREFIX  = "wfit_${tmpInput}_${FPREFIX}_MUOPT${sss}" ;  }

	    $CMDwfit = 
		" $JOBNAME_WFIT $inFile_tmp $WFIT_OPT " .
		"-cospar ${wPREFIX}.COSPAR " .
		"-resid  ${wPREFIX}.RESID " .
		"-label  $MUOPT_TAGNAME[$IOPT_S2MU] " . 
		"   >&   ${wPREFIX}.LOG " ; 
	    &add_COMMAND($icpu, "" ) ;
	    &add_COMMAND($icpu, $CMDwfit ) ;
	}
	
    } # IOPT_S2MU

    return ;

} # end of prep_COMMAND


sub add_COMMAND {
    my ($icpu,$command) = @_ ;
    my $N = $NCMDLINE_PER_CPU[$icpu] ;

    if ( $N==0 ) {
	$CMD_CONTENTS[$icpu][$N] = "$SNANA_LOGIN_SETUP";
	$N++ ;  $NCMDLINE_PER_CPU[$icpu]=$N ;
    }

    $CMD_CONTENTS[$icpu][$N] = "$command" ;
    $NCMDLINE_PER_CPU[$icpu]++ ;

} # end add_COMMAND

# =======================================
sub write_COMMANDS {

    if ( $SUMMARY_FLAG ) { return ; }

    my($icpu, $cmdFile, $i, $CMD, $PREFIX );

    for($icpu = 0 ; $icpu < $NCPU; $icpu++ ) {

	$cmdFile = "$FITSCRIPTS_DIR/$CMD_FILES[$icpu]" ;
	open PTR_CMD , ">> $cmdFile" ;       

	for($i=0; $i < $NCMDLINE_PER_CPU[$icpu]; $i++ ) {
	    $CMD = "$CMD_CONTENTS[$icpu][$i]" ;
	    print PTR_CMD "$CMD\n";
	}

	# create DONE file
	$PREFIX = "$CMD_PREFIX[$icpu]" ;
	print PTR_CMD "\n touch $FITSCRIPTS_DIR/${PREFIX}.DONE \n";

	close PTR_CMD ;
    }

} # end write_COMMANDS

# =======================================
sub USE_FITOPT{

    # Created May 17 2017
    # return 1 of this MUOPT,FITOPT should be processed;
    # return 0 to skip 
    # Based on new user-input option to list  FITOPT integers
    # with each MUOPT.
    # 
    # FITOPT integer list is stored in FITOPT_LIST[$imu]
    #

    my($imu,$ifitopt) = @_ ;

    if ( $FITOPT000_ONLY ) {
	if ( $ifitopt==0 ) { return(1); } else { return(0); }
    }

    my $FITOPT   = "$FITOPT_LIST[$imu]";
    my @IFITOPT  = split(/\s+/,$FITOPT) ;

    if ( $FITOPT eq "ALL" ) { return(1); }

    if ( grep { $_ == $ifitopt } @IFITOPT )
    { return(1); }
    else
    { return(0); }

} # end USE_FITOPT

# ================================
sub matchDump {

    my($idir, $iver, $NVER, $V, $Vmatch);
    
    for($idir=0; $idir < $N_INPDIR; $idir++ ) {
	$NVER = $NVERSION_SNFIT[$idir]; 
	for($iver=0; $iver < $NVER; $iver++ ) {
	    $V      = $VERSION_SNFIT_LIST[$idir][$iver] ;
	    $Vmatch = $VERSION_4MATCH_LIST[$idir][$iver] ;
	    print "     $V -> $Vmatch \n";
	}
    }
    
    #	    $VERSION_SNFIT_LIST[$NDIR][$NVER] = $vers ;
    #	    $VERSION_4MATCH_LIST[$NDIR][$NVER] = $vers4match ;
    
} # matchDump

# ================================
sub get_TABLE_FORMAT {

    my ($INPDIR, $VERSION, $SNFIT_DIR) ;

    # check for his/root files.
    # Fill global logical array  @TABLE_FORMAT_USE
    # and $COMBINE_FMTARG

    my ($suf, @bla, $N, $FMT, $FMTARG);

    $INPDIR  = $INPDIR_SNFIT_LIST[0] ;    # idir=0
    $VERSION = $VERSION_SNFIT_LIST[0][0]; # idir, iver
    $SNFIT_DIR = "$INPDIR/$VERSION" ;

    $N = 0 ;
    $COMBINE_FMTARG = "" ;
    foreach $suf ( @TABLE_SUFFIX_LIST ) {
	@bla    = qx(ls $SNFIT_DIR/FITOPT*.${suf} 2>/dev/null );
	$FMT    = $TABLE_FORMAT_LIST[$N];
	$FMTARG = $COMBINE_FMTARG_LIST[$N] ; 
	$TABLE_FORMAT_USE[$N] = 0 ;

	if ( scalar(@bla) > 0 ) { 
	    $TABLE_FORMAT_USE[$N] = 1 ; 
	    $COMBINE_FMTARG = "$COMBINE_FMTARG $FMTARG" ;
#	    print "  Found input $FMT format -->";
#	    print "  generate $FMT output for SALT2mu ($COMBINE_FMTARG)\n";
	}

	$N++ ;
    }

} # end of get_TABLE_FORMAT 

# ================================
sub submit_JOBS {

    my ($icpu, $script, $node, $cmd );
    my $qq = '"' ;

    if ( $SUMMARY_FLAG  ) { return; }  

    print "\n";

    for($icpu = 0 ; $icpu < $NCPU; $icpu++ ) {

	$script = $CMD_FILES[$icpu] ;	

	if ( $SSH_NNODE  ) {
	    $node   = $SSH_NODELIST[$icpu] ;
	    $cmd = "ssh -x $node ${qq} cd $FITSCRIPTS_DIR ;" . 
		" source $script ${qq}" ;
	    print " Submit SALT2mu jobs to $node . \n";
	    system("$cmd &") ;
	}
	elsif ( $BATCH_NCORE ) {
	    my $batchFile = $BATCH_FILES[$icpu] ;
	    print " Submit SALT2mu batch jobs with $batchFile . \n";
	    qx(cd $FITSCRIPTS_DIR ; $BATCH_COMMAND $batchFile);
	}
	else {
	    # interactive
	    print " Run interactive SALT2mu job $script. \n";
	    system("cd $FITSCRIPTS_DIR ; source $script");
	}
	   
    }  # icpu


    $| = 1;  # auto flush stdout.

    my $OUTDIR ;
    print "\n" ;
    print "# =============================================== \n";
    print " See SALT2mu output in \n" ;
    foreach $OUTDIR ( @OUTDIR_SALT2mu_LIST )  { print "    $OUTDIR \n" ; }
    print " Wait for  $NTOT_JOBS  SALT2mu jobs to finish. \n";
    $| = 1;  # auto flush stdout.

} # end of submit_JOBS


sub submit_SUMMARY {
    my ($CMD, $LOG);

    sleep(2);
    $CMD = "$SCRIPTNAME_FULL $INPUT_FILE SUMMARY NOPROMPT" ;

    # Jun 24 2020: pass NSPLITRAN so that summary job knows how many 
    # CPUs and DONE files to check.
    if ( $NSPLITRAN > 0 )  { $CMD = "$CMD NSPLITRAN=$NSPLITRAN" ; }

    print "\n\t Submit SUMMARY task in background. \n\n";
    $| = 1;  # auto flush stdout.

    system("$CMD  & " ) ;
    sleep(2);

#    $LOG = "$MERGE_LOGDIR/MERGE2.LOG" ;
#    system("$CMD > $LOG & " ) ;

}

# ====================
sub wait_for_done {

    # Called during summary stage to wait for DONE files.
    # Sep 28 2019: $DONESPEC argument -> '$DONESPEC' (in quotes)

    my ($NDONE, $DONESPEC, $CMD_WAIT );

# xxx not needed (Sep 28 2019)  if ( $SUMMARY_FLAG  ) { return; } 

    $NDONE        = $NCPU ;
    $DONESPEC     = "$FITSCRIPTS_DIR/*.DONE" ;
    $ALLDONE_FILE = "$FITSCRIPTS_DIR/ALL.DONE" ;

    $CMD_WAIT = "wait_for_files.pl  $NDONE  '$DONESPEC'  $ALLDONE_FILE" ; 


    system("$CMD_WAIT");

    # June 10 2019: check for ABORTS  
    my ( $cmd, @tmp, $OUTDIR, $msg);
    $NJOB_ABORT = 0 ;
    foreach $OUTDIR ( @OUTDIR_SALT2mu_LIST )  { 
	@tmp = qx(grep ' ABORT ' $OUTDIR/*/SALT2*.LOG );
	$NJOB_ABORT += scalar(@tmp);       
    }

}  # end of wait_for_done

sub create_done_file {

    # Created Sep 29 2019:
    # Create DONE file(s) here to ensure that everything
    # is really finished before DONE file(s) appear.

    my ( $msg );

    if ( $NJOB_ABORT > 0 )  
    { $msg = "FAILED  ($NJOB_ABORT ABORTS)"; }
    else
    { $msg = "SUCCESS" ; }

    print " Final STATUS for DONE file: $msg \n";
    qx(echo '$msg' >> $ALLDONE_FILE);

# Sep 12 2019: optional user-done file from DONE_STAMP key
    if ( $DONE_STAMP_FILE ne "" ) {
	qx(echo '$msg' >> $DONE_STAMP_FILE );
    }

} # end create_done_file

# ==========================
sub make_SUMMARY {

    my ($NDIR, $idir, $NVER, $iver, $ifitopt);

    # create two summary files.
    # 1. LOG file with user-input info such as MUOPT and INPDIR
    #      --> human readable
    # 2. DAT file in fitres-format 
    #      --> parsable 
    #
    # Aug 25 2017: add BETA1[_ERR] to varnames list
    # Nov 30 2017: add GAMMA0[_ERR] to varnames list
    
# xxxx    if ( $NSPLITRAN > 0 ) { return ; } # no summary for SPLITRAN

    my $SUMMARY_LOGFILE = "$FITSCRIPTS_DIR/FITJOBS_SUMMARY.LOG" ;
    my $SUMMARY_DATFILE = "$FITSCRIPTS_DIR/FITJOBS_SUMMARY.DAT" ;

    if ( $NSPLITRAN > 0 ) 
    { write_SPLITRAN_SUMMARY_LOG($SUMMARY_LOGFILE); return; }

    print "\n";
    print " Creating human-readable SUMMARY file \n\t $SUMMARY_LOGFILE \n";
    print " Creating code-parsable  SUMMARY file \n\t $SUMMARY_DATFILE \n";
    print " please wait ... \n" ;

    $NOUTFILE = 0 ;
    open  PTR_SUMLOG , "> $SUMMARY_LOGFILE" ;
    open  PTR_SUMDAT , "> $SUMMARY_DATFILE" ;  $NROW_SUMDAT=0;

    print PTR_SUMLOG  "\t\t SUMMARY of $NTOT_JOBS SALT2mu FITS \n\n";

    print PTR_SUMDAT "VARNAMES: " .
	"ROW VERSION IVERSION FITOPT MUOPT " . # +5
	"NSNFIT CHI2RED_1A " .                 # +2
	"ALPHA ALPHA_ERR BETA BETA_ERR " .     # +4
	"BETA1 BETA1_ERR " .                   # +2  // dbeta/dz
	"GAMMA0 GAMMA0_ERR " .                 # +2  // mag-hostmass
	"SIGINT M0AVG " .                      # +2
	"\n";

#    print PTR_SUMDAT "#  Note: CHI2 = -2log(L) for BCD fit \n";

    $NDIR = scalar(@OUTDIR_SALT2mu_LIST) ;

    print "\n xxx NDIR=$NDIR NVER=$NVERSION_FINAL[0]  NFITOPT=$NFITOPT_SNFIT[0] \n\n";

    for($idir=0 ; $idir < $NDIR; $idir++ ) {
	$NVER   = $NVERSION_FINAL[$idir] ;
	for ( $iver=0; $iver  < $NVER; $iver++ ) {
	    for($ifitopt=0; $ifitopt < $NFITOPT_SNFIT[$idir]; $ifitopt++ ) {
		&write_SUMMARY_LOG($idir,$iver,$ifitopt);
		&write_SUMMARY_DAT($idir,$iver,$ifitopt);

	    }  # ifitopt
	}    # iver
    }  # idir


    $T_END = time();  # gmtime();
    $T_TOT = ($T_END - $T_START)/60.0 ;
    $T_TOT = sprintf("%.1f", $T_TOT );

    print PTR_SUMLOG "\n" ;
    print PTR_SUMLOG " Total wall time: $T_TOT minutes. \n";
    print PTR_SUMLOG "\n" ;
    close PTR_SUMLOG ;

    print "\n";
    print " Found $NOUTFILE of $NTOT_JOBS output fitres files. \n";
    print " Total wall time: $T_TOT minutes. \n";
    print " Done. \n";

    $| = 1;  # auto flush stdout.

} # end of make_summary

# ============================================
sub write_SUMMARY_LOG {
    
    # for this [directory:version:fitopt],
    # write one-row summary to  the human-readable file (PTR_SUMLOG)
    # Loop over MUOPT and write one-row summary for each.

    my ($idir,$iver,$ifitopt) = @_ ;

    my ($INPDIR, $OUTDIR, $VERSION, $SPREFIX, $OUTFILE, $outFile, $IOPT_S2MU );
    my ($cVER, $nnn, $sss, @tmp, @wdlist, $FIRST, $USE, $SKIP);
    my ($chi2, $alpha, $beta, $sigmB, $M0avg );

    $OUTDIR  = $OUTDIR_SALT2mu_LIST[$idir] ;
    $VERSION = $VERSION_FINAL_LIST[$idir][$iver] ;

    $cVER    = sprintf("%-28.28s", $VERSION);
    $nnn     = sprintf("%3.3d", $ifitopt);

    # summarize MUOPT in header for human-readable file
    $FIRST = ( $idir == 0 && $iver == 0 && $ifitopt == 0 ) ;
    if ( $FIRST && $NMUOPT_SALT2mu > 0 ) {
	for ( $IOPT_S2MU=0; $IOPT_S2MU < $NMUOPT_SALT2mu; $IOPT_S2MU++ ) {
	    $sss     = sprintf("%3.3d", $IOPT_S2MU);
	    print PTR_SUMLOG " MUOPT: $sss $MUOPT_LIST_ORIG[$IOPT_S2MU] \n";
	}
	print PTR_SUMLOG "\n";
    }
    
    if ( $iver == 0 && $ifitopt == 0 ) { 
	print PTR_SUMLOG "\n";
	write_SUMMARY_INPDIR($idir);
	print PTR_SUMLOG " OUTDIR: $OUTDIR \n" ; 
	print PTR_SUMLOG "   VERSION  \t\t   FITOPT MUOPT    chi2/dof   " . 
	    "alpha  beta   sigint M0avg\n";
	PTR_SUMLOG -> autoflush(1);
    }

    for ( $IOPT_S2MU=0; $IOPT_S2MU < $NMUOPT_SALT2mu; $IOPT_S2MU++ ) {
	$SPREFIX = $SPREFIX_LIST[$idir][$iver][$ifitopt][$IOPT_S2MU] ;
	$outFile = "${SPREFIX}.FITRES" ;
	$OUTFILE = "${OUTDIR}/${VERSION}/${outFile}" ;
	$USE     = &USE_FITOPT($IOPT_S2MU,$ifitopt);
	$SKIP    = 0 ;
	if ( -e $OUTFILE ) 
	{ $NOUTFILE++ ; } 
	elsif ( $USE == 1 ) { 
	    print PTR_SUMLOG " ${VERSION}/$outFile  not  found ??? \n";
	    PTR_SUMLOG -> autoflush(1);
	    $SKIP = 1;
	}
	elsif ( $USE == 0 ) {
	    print PTR_SUMLOG " ${VERSION}/$outFile -> SKIPPED \n";
	    PTR_SUMLOG -> autoflush(1);
	    $SKIP = 1;
	}

	if ( $SKIP ) { next; }

	$sss     = sprintf("%3.3d", $IOPT_S2MU);

	@tmp  = sntools::parse_line($OUTFILE, 4, "chi2(Ia)/dof", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$chi2    = sprintf("%12.12s", $wdlist[3] );
	
	@tmp    = sntools::parse_line($OUTFILE, 2, "alpha0", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$alpha  = sprintf("%5.3f", $wdlist[1]);

	@tmp     = sntools::parse_line($OUTFILE, 2, "beta0", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$beta    = sprintf("%5.3f", $wdlist[1]);

	@tmp     = sntools::parse_line($OUTFILE, 2, "sigint", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$sigmB   = sprintf("%5.3f", $wdlist[1]);

	@tmp     = sntools::parse_line($OUTFILE, 2, "M0avg", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$M0avg   = sprintf("%7.3f", $wdlist[1]);
	
	print PTR_SUMLOG 
	    " $cVER  $nnn $sss  $chi2  " . 
	    "$alpha  $beta  $sigmB  $M0avg\n";
	PTR_SUMLOG -> autoflush(1);
    }

    return;

}  # end of write_SUMMARY_LOG

# ============================================     
sub write_SPLITRAN_SUMMARY_LOG {
    my ($SUMMARY_LOGFILE) = @_ ;

    # Created Jun 25 2020
    # hack to write SUMMARY.LOG for NSPLITRAN so that Pippin
    # will continue with makeCov and CosmoMC
    # Note that FITOPT000 and MUOPT000 are hard-wired since
    # NSPLITRAN does not work for FITOPTs and MUOPTs.

    open  PTR_SUMLOG , "> $SUMMARY_LOGFILE" ;
    print PTR_SUMLOG  "\t\t SUMMARY of $NTOT_JOBS SALT2mu FITS \n\n" ;
    print PTR_SUMLOG "MUOPT: 000 [DEFAULT] NONE \n";
    print PTR_SUMLOG "OUTDIR: $OUTDIR_OVERRIDE \n";
    print PTR_SUMLOG "\n";
    print PTR_SUMLOG "  VERSION FITOPT MUOPT \n";
    print PTR_SUMLOG "    0001   000    000 \n";

    close PTR_SUMLOG ;

} # 

# ============================================
sub write_SUMMARY_DAT {
    
    # for this [directory:version:fitopt],
    # write one-row summary to  the code-parsable file (PTR_SUMDAT)
    # Loop over MUOPT and write one-row summary for each.
    #
    # Aug 25 2017: include BETA1 and its error.
    # Nov 30 2017: include GAMMA0 and its error
    
    my ($idir,$iver,$ifitopt) = @_ ;

    # ------------ BEGIN ------------

    my ($OUTDIR, $VERSION, $SPREFIX, $OUTFILE, $outFile, $IOPT_S2MU );
    my ($cVER, @tmp, @wdlist, $KEY, $USE );
    my ($alpha, $alpha_err, $beta0, $beta0_err, $sigint, $M0avg );
    my ($beta1, $beta1_err, $gamma0, $gamma0_err);
    my ($NSNFIT, $chi2red);

    $OUTDIR  = $OUTDIR_SALT2mu_LIST[$idir] ;
    $VERSION = $VERSION_FINAL_LIST[$idir][$iver] ;
    if ( length($VERSION) < 28 ) 
    { $cVER = $VERSION ; }
    else 
    { $cVER    = sprintf("%-28.28s", $VERSION); }
    
    for ( $IOPT_S2MU=0; $IOPT_S2MU < $NMUOPT_SALT2mu; $IOPT_S2MU++ ) {
	$SPREFIX = $SPREFIX_LIST[$idir][$iver][$ifitopt][$IOPT_S2MU] ;
	$outFile = "${SPREFIX}.FITRES" ;
	$OUTFILE = "${OUTDIR}/${VERSION}/${outFile}" ;
	$USE     = &USE_FITOPT($IOPT_S2MU,$ifitopt);
	if ( $USE == 0 ) { next; }
	if ( !(-e $OUTFILE) ) { next; }

	@tmp     = sntools::parse_line($OUTFILE, 2, "NSNFIT", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$NSNFIT  = sprintf("%d", $wdlist[1] );

	$KEY     = "chi2(Ia)/dof" ;
	@tmp     = sntools::parse_line($OUTFILE, 4, "$KEY", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$chi2red = sprintf("%.3f", $wdlist[3] );
	
	@tmp       = sntools::parse_line($OUTFILE, 4, "alpha0", $OPT_QUIET) ;
	@wdlist    = split(/\s+/,$tmp[0]) ;
	$alpha     = sprintf("%6.4f", $wdlist[1]);
	$alpha_err = sprintf("%6.4f", $wdlist[3]);

	# - - - - - - - - - - 
	@tmp       = sntools::parse_line($OUTFILE, 4, "beta0", $OPT_QUIET) ;
	@wdlist    = split(/\s+/,$tmp[0]) ;
	$beta0     = sprintf("%5.3f", $wdlist[1]);
	$beta0_err = sprintf("%5.3f", $wdlist[3]);

	@tmp      = sntools::parse_line($OUTFILE, 4, "beta1", $OPT_QUIET) ;
	if ( scalar(@tmp) > 0 ) {
	    @wdlist   = split(/\s+/,$tmp[0]) ;
	    $beta1     = sprintf("%6.3f", $wdlist[1]);
	    $beta1_err = sprintf("%5.3f", $wdlist[3]);
	}
	else {
	    $beta1     = sprintf("%6.3f", 0.0 );
	    $beta1_err = sprintf("%5.3f", 0.0 );
	}

	# - - - - - - - - - - 
	@tmp        = sntools::parse_line($OUTFILE, 4, "gamma0", $OPT_QUIET) ;
	if ( scalar(@tmp) > 0 ) {
	    @wdlist     = split(/\s+/,$tmp[0]) ;
	    $gamma0     = sprintf("%5.3f", $wdlist[1]);
	    $gamma0_err = sprintf("%5.3f", $wdlist[3]);
	}
	else {
	    $gamma0     = sprintf("%6.3f", 0.0 );
	    $gamma0_err = sprintf("%5.3f", 0.0 );
	}

	# - - - - - - - - - -	
	@tmp     = sntools::parse_line($OUTFILE, 2, "sigint", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$sigint  = sprintf("%5.3f", $wdlist[1]);

	@tmp     = sntools::parse_line($OUTFILE, 2, "M0avg", $OPT_QUIET) ;
	@wdlist  = split(/\s+/,$tmp[0]) ;
	$M0avg   = sprintf("%7.3f", $wdlist[1]);
	
	$NROW_SUMDAT++ ;
	print PTR_SUMDAT
	    "ROW: $NROW_SUMDAT $cVER $iver $ifitopt $IOPT_S2MU " .
	    "$NSNFIT $chi2red " . 
	    "$alpha $alpha_err $beta0 $beta0_err $beta1 $beta1_err " .
	    "$gamma0 $gamma0_err " .
	    "$sigint  $M0avg\n";
	PTR_SUMDAT -> autoflush(1);
    }

    return;

}  # end of write_SUMMARY_DAT


# ========================
sub write_SUMMARY_INPDIR {
    my ($idir) = @_ ;

    my ($INPDIR);

    if ( $SUMFLAG_INPDIR == 0 ) {
	# list only the INPDIR corresponding to $idir
	$INPDIR  = $INPDIR_SNFIT_LIST[$idir];
	print PTR_SUMLOG " INPDIR: $INPDIR \n" ; 
    }
    else {
	# list all input dirs
	foreach $INPDIR (@INPDIR_SNFIT_LIST) {
	    print PTR_SUMLOG " INPDIR+: $INPDIR \n" ; 
	}
    }
} # end of  write_SUMMARY_INPDIR 


sub gzip_logs {
    my($tarFile, $prefix);

    if ( $NSPLITRAN > 0 ) { return ; }
    
    $tarFile = "${FITSCRIPTS_PREFIX}.tar" ;
    $prefix  = "${FITSCRIPTS_PREFIX}_CPU" ; # prefix to tar up
    
    qx(cd $FITSCRIPTS_DIR ; tar -cf $tarFile ${prefix}* );
    qx(cd $FITSCRIPTS_DIR ; gzip    $tarFile            );
    sleep(1) ;
    qx(cd $FITSCRIPTS_DIR ; rm      ${prefix}*          );    
} 


# ============================
sub CLEANFILES_DRIVER {

    # Created Dec 3 2017
    # + remove SALT2mu*HBOOK and ROOT files
    # + gzip SALT2*FITRES files.

    # ------------ BEGIN -----------
    my @FLIST_HBOOKROOT = ("SALT2mu\*HBOOK", "SALT2mu\*ROOT", "SALT2mu\*hbook");
    my @FLIST_FITRES    = ("SALT2mu\*FITRES" );

    sntools::clean_fileList($PROMPT,@FLIST_HBOOKROOT);

    &clean_gzipSALT2mu();

    print " Done. \n";
    return ;

} # end CLEANFILES_DRIVER


sub clean_gzipSALT2mu {

    my ($cmd_find, $response, @tmp );

    if ( $PROMPT == 0 ) {
	$response = 'y' ;
    }
    else {
	print "\n Enter 'y' to gzip all " . 
	    "SALT2mu*.FITRES  and SALT2mu*.LOG  files => \n  "; 
	$response = <STDIN> ;     $response =~ s/\s+$// ;
    }

    if ( $response eq "y" ) {
	printf "\t gzipping SALT2mu*.FITRES . . . \n";
	$cmd_find = "find . -name SALT2mu\*.FITRES -exec gzip {} +";
	@tmp = qx($cmd_find);

	printf "\t gzipping SALT2mu*.LOG . . . \n";
	$cmd_find = "find . -name SALT2mu\*.LOG -exec gzip {} +";
	@tmp = qx($cmd_find);
    }
    else {
        print " SALT2mu*FITRES/LOG Files NOT gzipped.\n" ;
    }

} # end clean_gzipSALT2mu
