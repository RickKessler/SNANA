#!/usr/bin/perl
#
# Created Jan 2011 by R.Kessler
# Jan 2019: Major refactor to 
#   + allow multiple SIMGEN_INFILE_NONIa keys
#   + idendity & skip duplicate jobs
#
# General script to simulate a survey with an SN-mix
# of Type Ia (arbitrary set of Ia models) and NONIa SNe, 
# each normalized to the number of seasons specified in
# the input file.  The Ia and NONIa CIDs are randomly 
# mixed so that any subset is a random subset. Before 
# the nominal generation, this script runs a job with N=0
# to determine the appropriate NGENTOT_LC for one season; 
# thus the <inputFile> specifies the number of seasons 
# rather than the number of SNe to generate.
# In addition to specifying an arbitrary mix of SN types,
# a list of systematic variations can also be given
# (see multiple GENVERSION and GENOPT keys in the master 
# input file)
#
#
# Usage  
#  sim_SNmix.pl <SNmixFile>
#  sim_SNmix.pl <SNmixFile> FAST10   (x10  fewer SN)
#  sim_SNmix.pl <SNmixFile> FAST100  (x100 fewer SN)
#  sim_SNmix.pl <SNmixFile> KILL
#  sim_SNmix.pl <SNmixFile> GENSTAT    ! debug GENSTAT
#  sim_SNmix.pl <SNmixFile> KICP       ! force kicp queue
#  sim_SNmix.pl <SNmixFile> NOSUBMIT   ! init, but don't submit
#  sim_SNmix.pl <SNmixFile> NOPROMPT   ! skip prompt if SIMLOGS has no DONE stamp
#                                      !  (make sure batch jobs are killed)
#
#  And for internal use only:
#    sim_SNmix.pl <SNmixFile> -NODEINDX <nodeindx>  -SUFFIX $Nsec5
#
# The <SNmixFile> contains general instructions, Note that <SNmixFile>
# is NOT a valid sim-input file, but rather it includes the names of 
# valid sim-input files for SNIa and SNCC. Example SNmix file is in
#
# $SNDATA_ROOT/sample_input_files/multijob_scripts/SIMGEN_MASTER_DES.INPUT
#
# Generation speed: ~ 10^4 per minute (with 25% written out)
#
# If the SNmix master input file contains
#   NODELIST: <node1> <node2> ...
# then the jobs are submitted in round-robin fashion on each node.
# The jobs are processed by re-calling this same script with the
# extra argument "-NODEINDX <nodeindx> -SUFFIX $Nsec5" for each node.
#
# Misc keys for master-input file:
#
# GENVERSION:  MYNAME_BLABLA
#   GENOPT:         <any command-line override to GENVERSION>
#   GENOPT(SNIA):   <apply option(s) only to SNIA)>
#   GENOPT(SNIa):   <same as above>
#   GENOPT(IA):     <same as above>
#   GENOPT(Ia):     <same as above>
#   GENOPT(NONIA):  <apply option(s) only to NONIA)>
#   GENOPT(NONIa):  <same as above>
#   GENOPT([inFile_subString]):  bla  # apply to jobs with subtring match
#                                     # to arg of SIMGEN_INFILE_[Ia,NONIa]
#   SIMGEN_INFILE_NONIa: <NONIa-override-input-file>
#   SIMGEN_INFILE_Ia:    <idem for SNIa model> 
#   SIMGEN_INFILE_SNIa:  <alternate key for SNIa input file>
#
#  GENVERSION: ANOTHER_SIMGEN_NAME
#    GENOPT:  <another command-line override>
#
# ENDLIST_GENVERSOIN:  (must end list of GENVERSIONs with this key)
#
#   GENPREFIX:      [recommend SURVEY_NAME]
#   FORMAT_MASK:    [MASK]    # can also be GENOPT_GLOBAL arg
#   GENOPT_GLOBAL:  [global command-line arg list for snlc_sim.exe]
# 
#   DONE_STAMP:   /home/bla/BLABLA.DONE  # global "all-done" stamp
#   RESET_CIDOFF: 2   # unique CID in every GENVERSION (if RANCID is set)
#   CIDRAN_MIN:   100000   # e.g., random CIDs > 100000
#                   (to generate sims with unique CID ranges)
#   CONVERT_SIMGEN_DUMP: HBOOK
#   CONVERT_SIMGEN_DUMP: ROOT
#   PATH_SNDATA_SIM:  <path>     # re-direct data file output
#   CLEANUP_FLAG:  0             # leave TMP_[user]_xxx verions (default is on)
#   DOSKIP_DUPLICATE_SIMJOBS: 0  # no duplicate check; brute-force all jobs
#                                # (default is on)
#
#  TOPDIR_OVERRIDE: <dirName> # logs->[TOPDIR_OVERRIDE]/SIMLOGS_[GENPREFIX]
#                             #   (default is current pwd)
#
#  LOGDIR: <dirName>          # write logs to [LOGDIR] or [pwd]/[LOGDIR]
#                             #  (default is [pwd]/SIMLOGS_[GENPREFIX]
#
# - - - - - - - - - - - - - - - - - - - - - - 
#
# HISTORY
#
#
# Dec 01 2017:
#   + in combine_random_simjobs(), remove duplicate NVAR and VARNAME
#     keys in catenated SIMGEN-DUMP file.
#
# Dec 9 2017: in simgen(), adjust location of NGENTOT_LC in argument list.
#
# Dec 18 2017: fix bug in INIT_SIMGEN(), initialize BATCH_TEMPLATE="",
#              not BATCH_TEMPLATE_KICP="". Works now for KICP q.
#
# Jan 11 2018: CIDRAN_MAX *= 1.1 to make finding random CIDs more efficient.
# Jan 13 2018: 
#   + Include "wrote" and Gen stats at end of summary file. 
#     See WROTE_STRING array.
#   + check for INPUT_FILE_INCLUDE among GENOPT, and copy
#      to /misc.
#
# Jan 24 2018: 
#   + gzip FITS files
#   + for RANlSEED_REPEAT, compress when each IVER is done instead of
#      waiting for all jobs to finish
#   + fix to compress RANSEED_CHANGE as well.
#
# Feb 09, 2018:
#   + fix a few issues writing summary file. Add .1 second delay
#     between each batch job and use LOCKFILE_DONE
#
# Apr 20 2018:
#  for get_normalization stage, add "SIMLIB_MAXRANSTART 0" argument
#
# May 10 2018: include NSPEC in summary files.
# 
# May 22 2018: new key FIRST_CIDOFF
# May 29 2018: allow comment string on GENOPT line, following hash (#)
# Jun 09 2018: pass argument "JOBID <JOBID>" to snlc_sim.exe
# Jul 06 2018: read both sim-input and sim-include CONTENTS for parsing.
#
# Jul 22 2018: 
#     + remove FIRST_CIDOFF logic since it doesn't work
#     + Instead, add new master-file input CIDRAN_MIN: <CIDMIN>
# 
# Aug 19 2018: pass "NJOBTOT <NJOBTOT>"
# Jan 16 2019: add exit(0) 
# Jan 21 2019: for combining SIMGEN-DUMP files, check for NVAR
#              in case it isn't there.
#
# Jan 23 2019: major overhaul to allow multiple SIMGEN_NONIa-input files.
#
# Feb 12 2019: 
#    + another major overhauld to specify different set of
#       SIMGEN-input files for each GENVERSION.
#    + abort if any command-line arg is unknown/unused.
#
# Mar 04 2019:
#   in combine_random_simjobs(), ls *SNIa* before *NONIa* for 
#   .LIST file so that fitting code sees Ia params first.
#
# Apri 17 2019: 
#   + fix horrible $m index bug in parse_INFILE_NONIa
#   + refactor to ignore catenating SIMGEN_DUMP files if SIMGEN_DUMP
#       flag is not specified.
#   + write SUCCESS or FAILURE in DONE_STAMP
#   + replace FATAL_ERROR with FATAL_ERROR_STAMP to write   
#     FAILURE in done-stamp.  
#
# Aug 4 2019:
#   sed -i 's/.gz//g' $LIST_FILE --> sed -i 's/\.gz//g' $LIST_FILE
#   to avoid remove _gz from filename. Seems that dot is a special char
#   that needs a backslash.
#
# Sep 7 2019: 
#    replace parse_line with parse_array so that comment lines are ignored.
#
# Sep 13 2019: add batchName arg to make_batchFile()
#
# Sep 30 2019: 
#   + for RANSEED_CHANGE option, fix LIST file to order SNIa before
#     NONIaMODEL, so that fitting code output has correct SIM_xxx
#     names for SNIa model. Previously this ordering worked only
#     for RANSEED_REPEAT, but not for RANSEED_CHANGE.
#   + fix subtle bug checking DONE_STAMP ... only affects
#     warning about missing DONE_STAMP before launching.
#
# Oct 25 2019:
#  Refactor to read things from GENOPT_GLOBAL 
#  (e..g, INCLUDE_FILE, FORMAT_MASK, GENPREFIX, NGENTOT_LC, ...).
#  Abort if H0 key is in master file ... give warning about ZRANGE key.
#
# Oct 29 2019:
#   for RANSEED_CHANGE, turn off option to avoid duplicate jobs,
#   and if RESET_CIDOFF=2, then set RESET_CIDOFF=1 so that 
#   unique CIDs are for each Ia/NONIa set, not for all sets.
#
# Nov 12 2019:
#   ABORT if normalization job fails; see $normLog -->
#   fixes infinite loop bug when QUIT key isn't in $normLog.
#
# Dec 4 2019: new input key  LOGDIR: <dirName>
#
# Jan 22 2020: protect GENOPT_GLOBAL for parentheses in argument.
#
# Feb 19 2020: 
#   + in write_GENSTAT_DRIVER(), increase max BUSYFILE wait time from
#     10s to 60s (for NERSC).
#
# Mar 12 2020: 
#   + add NOPROMPT option to clobber SIMLOGS subDir even if there is no
#     DONE stamp.
#
# Mar 24 2020: fix so that SIMGEN_INFILE_Ia need not be global if
#              there is a SIMGEN_INFILE_Ia for each GENVERSION
#
# Apr 7 2020:
#    leave FAILURE message in DONE stamp if SBATCH file does not exist.
#
# Apr 24 2020:
#   +  abort if PATH_SNDATA_SIM does not exist.
#   +  abort if NPER_SEASON is null string
#
# Apr 28 2020:
#   + abort if SIMGEN_INFILE_[Ia,SNIa,NONIa] is not specified.
#       (github issues from Justin)
#
# May 12 2020: 
#   + increment NABORT if logFile does not exist.
#   + write PATH_SNDATA_SIM in TOTAL_SUMMARY.LOG
#
# Jul 14 2020: 
#   fix $MOI4 (for version prefix) to append zeros if username
#   has less then 4 characters. e.g., user = ab will have
#   prefix = ab00
#     
# Jul 24 2020: no longer copy SIMLIB_FILE to misc/
#
# Aug 15 2020: check GENOPT_GLOBAL for include file.
#
# ---------------------------------------------------------

use strict ;
use File::Copy;
use List::Util qw(first) ;
use IO::Handle ;
use FindBin qw($Bin);
use lib "$Bin" ;
use Time::HiRes qw(usleep);
use sntools ;

# --------------------------------------
# declare subs

sub init_SIMGEN;
sub parse_arg;
sub parse_inFile_master;
sub strip_GENOPT ;
sub syncVerify_files ;
sub syncVerify_GENPAR ;
sub SUBMIT_NODES ;
sub make_logDir ;
sub get_normalization ;
sub get_normalization_model ;
sub get_normalization_SUM ;
sub init_duplicate ;
sub make_symLinks_duplicate ;
sub final_SIMGEN_DUMP_duplicate ;
sub append_normLog ;
sub get_NGEN ;
sub check_node ;
sub PROCESS_SIMJOB;
sub set_timeStamp ;

sub set_version_names ;
sub reset_counters;
sub set_jobFileNames ;
sub create_version ;
sub simgen ;
sub simcopy ;
sub make_AUXFILE_IGNORE ;
sub make_AUXFILE_README ;
sub make_AUXFILE_FILTERS ;
sub make_AUXFILE_LIST ;
sub make_DONEFILE ;
sub printLegend_INFILE ;
sub printLegend_GENOPT ;
sub printSummary ;
sub cleanup ;
sub parse_INFILE_Ia ;
sub parse_INFILE_NONIa ;
sub parse_GENOPT_GLOBAL ;
sub parse_GENOPT ;
sub parse_GENVERSION ;
sub parse_RANSEED ;
sub parse_FORMAT_MASK ;
sub check_FORMAT_MASK ;
sub check_SIMGEN_DUMP ;
sub get_jobidString ;
sub PROCESS_GENOPT_LINE ;
sub read_GENOPT_STRING ;

sub combine_random_simjobs ;
sub move_random_SIMGEN_DUMP ;
sub cat_random_SIMGEN_DUMP ;
sub convert_SIMGEN_DUMP ;
sub compress_output ;
sub check_ALLDONE ;
sub check_VersionDone ;
sub check_COMBINE ;
sub make_TOTAL_SUMMARY_FILE ;

sub write_GENSTAT_DRIVER; # append_GENSTAT
sub write_GENSTAT_open ;
sub write_GENSTAT_job ;
sub write_GENSTAT_sum ;
sub write_GENSTAT_table ; 
sub move_GENSTAT_FILE ;
sub get_FILENAME_GENSTAT;

sub makeTar_SIMLOGS ;
sub store_SIMGEN_INFILE ;

# --------------------------------------
# declare globals

my ($SNDATA_ROOT,$SNANA_DIR,$USER,$HOST,$currentDir);
my ($SNANA_MODELPATH, $ENVDEF_MODELPATH, $SHELL_NAME );
my ($DO_SIMJOB, $NRANJOB_SIMGEN,  $NGENVERSION, $NSIMTYPE);
my (@GENVERSION_NAME, $VERSION_TEMP, $VERSION_FINAL, @BATCH_MEM );
my ($SIMGEN_TEMPDIR, $SIMGEN_FINALDIR, $SIMGEN_MISCDIR, $MISC_SDIR  );
my ($READMEFILE_FINAL, $LISTFILE_FINAL, $IGNOREFILE_FINAL, $DUMPFILE_FINAL);
my ($READMEFILE_TEMP, $DUMPFILE_TEMP, $DONE_STAMP, $DONE_STAMP_FLAG );
my ($LOGDIR, $TOPDIR_SIMLOGS,  $Nsec5 );
my (@NSIM_GEN, @NSIM_WR, @NSIM_SPEC);
my (@VERSION_JOBLIST_FINAL, @VERSION_JOBLIST_TEMP);
my (@STATUS_NORMALIZATION, @TOTAL_STRING );


my (@CONTENTS_INFILE_MASTER, @CONTENTS_INFILE_INCLUDE );
my (@CONTENTS_INFILE_Ia, @CONTENTS_INFILE_NONIa, @CONTENTS_INFILE_SIMGEN);

my ($NARG,$SNMIX_INFILE_MASTER,$SUBMIT_FLAG, $PROMPT_FLAG );
my ($INPUT_FILE_INCLUDE_Ia, $INPUT_FILE_INCLUDE_NONIa, @GENOPT_FILE_INCLUDE);
my ($INPUT_FILE_ZVAR_Ia, $INPUT_FILE_ZVAR_NON1A, @INPUT_FILE_NON1AGRID );
my ($SEARCHEFF_SPEC_FILE, $HOSTNOISE_FILE, $ZPHOTEFF_FILE );
my ($FASTFAC, $SUFFIX_TEMP );
my ($DO_SSH, $DO_BATCH, $SNANA_LOGIN_SETUP);
my ($KILLJOBS_FLAG, $DEBUG_GENSTAT, $INODE_GLOBAL, $PROMPT_FLAG );
my ($DOGEN_SNIa, $DOGEN_NONIa, @NVAR_SIMGEN_DUMP);
my ($CONVERT_SIMGEN_DUMP, $DO_SIMGEN_DUMP);
my ($NCLASS_SIMGEN_DUMP, $NCLASS_SIMGEN_DUMPALL );
my ($OPT_ABORT,$OPT_WARN, $OPT_QUIET, $JOBABORT_FLAG );
my ($MXCID_SIM, $CIDRAN_MAX, $CIDRAN_OFF, @CIDRAN_OFFSET );
my ($CIDOFF, $RESET_CIDOFF, $CIDRAN_MIN );
my (@GENVERSION_NAME, @GENVERSION_GENOPT, $GENPREFIX);
my ($DOSKIP_DUPLICATE_SIMJOBS,@IVERSION_DUPLICATE, @IMODEL_DUPLICATE);
my (@IVERSION_DUPLICATE_INVERSE, @IMODEL_DUPLICATE_INVERSE );
my (@NMODEL_DUPL, $NMODEL_DUPL_TOT); 
my (@DUMPFILE_DUPLICATE, @STRING_COMMENT_DUPLICATE);
my ($LAST_NGEN, $LAST_NSIMTOT_WR);
my ($GENFILTERS,$NFILT, $CLEANUP_FLAG, $GENOPT_GLOBAL);
my ($NGENMODEL_MAX, $mIa );

my (@NGENMODEL, @GENMODEL_NAME, @GENMODEL_CLASS, @SIMGEN_INFILE );
my ($NGENMODEL_GLOBAL, @GENMODEL_NAME_GLOBAL);
my (@GENMODEL_CLASS_GLOBAL, @SIMGEN_INFILE_GLOBAL );

my (@GENMODEL_NGENTOT  );
my ($NSIMTOT_GEN, $NSIMTOT_WR, $NSIMTOT_SPEC);
my (@NSIM_GEN, @NSIM_WR, @NSIM_SPEC);
my ($NSIMARG, @SIMARG_LIST, $SIMARG_ALL);
my ($SIMARG_FORMAT);
my (@SIMARG_RANSEED, $SAMEFLAG_RANSEED );
my ($FORMAT_MASK, $WRFLAG_FITS, $WRFLAG_TEXT, $WRFLAG_CIDRAN);
my ($NGEN_UNIT_VAL, $NGEN_UNIT_NAME );
my (@MSGERR, $NABORT );
my ($PATH_SNDATA_SIM, $USER_PATH_SNDATA_SIM );

my (@GENSTAT_SUMJOB, @GENSTAT_SUMTOT, @GENSTAT_ivLink, @GENSTAT_mLink);


# ssh/batch variables
my ($SSH_NNODE, @SSH_NODELIST, @CPULIST_INDX, $NCPU, $NODEINDX );
my ($BATCH_COMMAND, $BATCH_TEMPLATE, $BATCH_NCORE );
my $BATCH_MEM = "1000" ;  # default 1GB memory, Mar 2 2017

my $JOBNAME_SIM = "snlc_sim.exe" ;
my $JOBNAME_SIM_FULLPATH = `which $JOBNAME_SIM` ;

# define temp prefix for versions that get deleted,
# or should be deleted by user of job fails
my ($MOI4, $PREFIX_TEMP, $SUFFIX_DUMP_TEMP) ;

# xxxxxxx mark delete Jul 2020 xxxxxxxxx
#my $MOI  = `whoami`  ;
#my $MOI4 = substr($MOI,0,4);
#my $PREFIX_TEMP  = "TMP_${MOI4}" ;
#my $SUFFIX_DUMP_TEMP = "DUMP_TEMP" ;
# xxxxxxxxxxxxxxxxxxxxx

my $BATCH_TEMPLATE_KICP = '$SBATCH_TEMPLATES/SBATCH_kicp.TEMPLATE' ;

my $OPT_GENSTAT_TABLE_plusLINKS = 1;
my $OPT_GENSTAT_TABLE_noLINKS   = 2;
my $TABLE_INSERT_MARK           = "TABLE_INSERT_MARK" ;

# =========== BEGIN MAIN ==============


print "\n";
print " ###################################### \n" ;
print " ############### BEGIN ################ \n" ;
print " ###################################### \n" ;
print "\n";

&init_SIMGEN();

&parse_arg();

&parse_inFile_master();

if ( $KILLJOBS_FLAG ) { sntools::Killjobs(@SSH_NODELIST) ;  die "\n"; }

&init_duplicate();

if ( $SUBMIT_FLAG == 0 ) {
    print "\n SUBMIT_FLAG=0 ==> DO NOT SUBMIT JOBS.\n";
    exit(0);
}

&make_logDir();

if ( $DEBUG_GENSTAT ) 
{ &write_GENSTAT_DRIVER(-1,-1); die "\n Done GENSTAT DEBUG\n";}

#die "\n xxx DEBUG DIE xxx \n";

# submit to nodes if NODELIST is given (but not the -NODEINDX argument)
&SUBMIT_NODES();

# -----------------------------------------------------
my ($IVER) ;

if ( $RESET_CIDOFF == 2 ) {
    # option to sum normalization over all GENVERSIONS
    for ( $IVER=1 ; $IVER <= $NGENVERSION; $IVER++ ) {
	&get_normalization($IVER);  
    }   
    print "\n  --> CIDRAN_MAX(sum over GENVERSIONS) = $CIDRAN_MAX \n\n";
}



for ( $IVER=1 ; $IVER <= $NGENVERSION; $IVER++ ) {

    my ( $JOBID,$VERSION, $MOD );
    $VERSION = $GENVERSION_NAME[$IVER] ;
    print "\n ############################################# \n";
    print " Begin Simulations for $VERSION ... \n";

    $STATUS_NORMALIZATION[$IVER] = 0 ;

    # each JOBID is a different RANSEED
    for ( $JOBID=1 ; $JOBID <= $NRANJOB_SIMGEN; $JOBID++ )  {

	&check_node();  # set DO_SIMJOB logical

	if ( $DO_SIMJOB ) {
	    if($RESET_CIDOFF != 2) { &get_normalization($IVER); }
	    &PROCESS_SIMJOB($IVER,$JOBID);
	    if ( &check_VersionDone($IVER) > 0 ) {
		&combine_random_simjobs($IVER); 
	    }
	    &check_ALLDONE();
	}
	else {
	    print "\t ==> SKIP  GENVERSION $VERSION  (JOBID $JOBID) \n";    
	}

    }  # end of JOBID loop

}  # end of iver loop  over GENVERSIONS

print "\n Done. \n";

exit(0);

# ===================================
#
#       END OF MAIN
#
# ===================================


# =====================================
sub PROCESS_SIMJOB {

    # iver = version index
    # jobid = index of random for multiple RANSEED keys in master
    
    my($iver,$jobid) = @_;

    my($mod, $NGENTOT_SUM, $NGEN_OFF );

    &set_version_names($iver,$jobid);

    &reset_counters();

    &set_jobFileNames($iver,$jobid);

    &create_version($iver);

    if ( $RESET_CIDOFF  ) {
	$CIDOFF    = 0 ;   # reset with each new job 
	$LAST_NGEN = 0 ; 
    }

    # Jun 5 2013:
    # if there are multiple RANSEED keys, then add CID offset so 
    # that each RANDOM sample has unique CIDs.
    $NGENTOT_SUM  = 0; 
    for ( $mod=0; $mod < $NGENMODEL[$iver]; $mod++ ) 
    { $NGENTOT_SUM += $GENMODEL_NGENTOT[$iver][$mod] ;   }
    $NGEN_OFF     = ($jobid-1) * $NGENTOT_SUM ;

    # July 17 2017: require SAMEFLAG_RANSEED to increment CIDOFF
    if ( $SAMEFLAG_RANSEED )  { $CIDOFF += $NGEN_OFF ; }

    # generate all models $m 
    for ( $mod=0; $mod < $NGENMODEL[$iver]; $mod++ ) { 
	if ( $IVERSION_DUPLICATE[$iver][$mod] >= 0 ) { 
	    next ;
	}
	else {
	    &simgen( $iver,$jobid,$mod);   # generate sim; wait to finish
	    &simcopy($iver,$jobid,$mod);   # copy files to VERSION_FINAL
	}
    }    # end of model loop

    &make_AUXFILE_IGNORE($iver,$jobid);
    &make_AUXFILE_README($iver,$jobid);
    &make_AUXFILE_LIST($iver,$jobid);
    &printSummary($iver);
    &cleanup($iver,$jobid);
    &make_DONEFILE($iver,$jobid);

    # append stats in table in final version 
    &write_GENSTAT_DRIVER($iver,$jobid); 

    # compress now if jobids are NOT combined later
    if ( $SAMEFLAG_RANSEED==0 ) { &compress_output($iver,$jobid); }

    return ;

} # end of PROCESS_SIMJOB



# =======================================
sub SUBMIT_NODES {

    # if SSH_NODELIST is specified, submit to all nodes.
    # if DO_SSH > 0, then this is already a launched job.
    # Launch sim_SNmix script, not snlc_sim job.
    # Apr 3 2019: don't allow NCPU > NJOBTOT

    my ($q,$cdir,$base, $iver,$j,$indx, $TEXT, $NJOBTOT );
    my ($cmd,$cmdNode, $keyNODEINDX, $str_indx, $cmdSetup);  

    if ( $NCPU       == 0 ) { return ; }
    if ( $DO_SSH      > 0 ) { return ; }
    if ( $DO_BATCH    > 0 ) { return ; }

    $q = '\'' ;
    $cdir  = "cd $currentDir";
    $base  = "$0 @ARGV[0 .. $NARG-1]" ;  # original user command

    $NJOBTOT = $NGENVERSION * $NRANJOB_SIMGEN ;
    if ( $NCPU > $NJOBTOT ) { 
	print " $NCPU CPUs > $NJOBTOT jobs; NCPU -> $NJOBTOT \n";
	$NCPU = $NJOBTOT ;  
    }

    &set_timeStamp();  # set Nsec5 = number of seconds since midnight

    # determine key based on SSH or BATCH system
    if ( $BATCH_NCORE > 0 ) 
    { $keyNODEINDX = "-NODEINDX_BATCH" ; $TEXT = "BATCH" ; }
    else
    { $keyNODEINDX = "-NODEINDX_SSH"   ; $TEXT = "SSH"  ;  }

    if ( length($SNANA_LOGIN_SETUP) > 0 ) 
    { $cmdSetup = "$SNANA_LOGIN_SETUP" ; }
    else
    { $cmdSetup = "set doNothing=1" ; }  # command to do nothing

    if ( length($SNANA_MODELPATH) > 0 ) 
	{ $cmdSetup = "$cmdSetup ; $ENVDEF_MODELPATH" ; }
   

    if ( $BATCH_NCORE && !(-e $BATCH_TEMPLATE)  ) {
	$MSGERR[0] = "Cannot find BATCH template file";
	$MSGERR[1] = "$BATCH_TEMPLATE";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);  
    }

    print "\n";
    print " Split SIMGEN jobs by nodes: \n";

    for( $j=0; $j < $NCPU;  $j++ ) {
	$indx   = $CPULIST_INDX[$j] ;
	$str_indx = &get_jobidString($indx );

	# append nodeindx and suffix to user command
	$cmdNode = "$base  $keyNODEINDX $indx  -SUFFIX ${Nsec5}" ; 

	if ( $BATCH_NCORE == 0 ) {   # SSH option
	    my $node     = $SSH_NODELIST[$j] ;
	    my $nodeLog  = "$LOGDIR/${GENPREFIX}-${str_indx}_${node}.stdout";
	    my $inQuotes = "$cdir ; $cmdSetup ; $cmdNode >& $nodeLog" ;
	    $cmd         = "ssh -x $node ${q} $inQuotes ${q}" ;
	    print "\t -> $node ($indx) \n" ;
	    system("$cmd &") ;
	}
	else {
	    # use batch system
	    my $batchName = "${GENPREFIX}_${str_indx}" ;
	    my $batchFile = "${GENPREFIX}_${str_indx}.BATCH" ;
	    my $batchLog  = "${GENPREFIX}_${str_indx}.LOG" ;
	    my $batchMem  = "$BATCH_MEM" ;
	    my $JOB       = "$cdir ; $cmdSetup ; $cmdNode" ;

	    print "\t prepare $batchFile  for $BATCH_COMMAND \n";
	   
	    sntools::make_batchFile($BATCH_TEMPLATE, $LOGDIR,
				    $batchName, $batchFile, $batchLog, 
				    $batchMem, $JOB);

	    qx(cd $LOGDIR ; $BATCH_COMMAND $batchFile);
	    usleep(100000);  # Feb 9 2018: 0.1 sec delays
	}
    }

    print "\n Done submiting $TEXT jobs. \n" ;
    exit(0);

} # end of SUBMIT_NODES

# ===============================
sub check_node {

    # set global DO_SIMJOB= 0 or 1 based on if this $INODE_GLOBAL
    # is a match
    # Note that $INODE_GLOBAL is globally incremented.

    my ($cpu_indx, $ctmp1, $ctmp2 );

    $DO_SIMJOB = 1;  

    if ( $DO_SSH || $DO_BATCH ) {
	$cpu_indx = $CPULIST_INDX[$INODE_GLOBAL];
	if ( $cpu_indx != $NODEINDX ) { $DO_SIMJOB = 0 ; }

	$INODE_GLOBAL++ ;  
	if ( $INODE_GLOBAL >= $NCPU ) { $INODE_GLOBAL = 0 ; }
    }


} # end of check_node

# ===============================
sub debug_abort {
    die "\n xxxxxx DEBUG ABORT xxxxxxxx \n";
}

# ========================================================
sub init_SIMGEN() {

    my ($m, $m0, $LEN_USER, $NCHAR_MOI ) ;

    $NGENMODEL_GLOBAL = 0 ;
    $mIa                = -9 ;
    $NGENMODEL_MAX      = 10 ;   

    $USER        = $ENV{'USER'};
    $HOST        = $ENV{'HOST'};
    $SNDATA_ROOT = $ENV{'SNDATA_ROOT'};
    $SNANA_DIR   = $ENV{'SNANA_DIR'};

# - - - - 
    $LEN_USER = length($USER);
    $NCHAR_MOI = 4 ;
    $MOI4      = substr($USER,0,$NCHAR_MOI);
    if ( $LEN_USER < $NCHAR_MOI ) 
    {  $MOI4 .= ('0' x ($NCHAR_MOI-$LEN_USER));  }

    $PREFIX_TEMP  = "TMP_${MOI4}" ;
    $SUFFIX_DUMP_TEMP = "DUMP_TEMP" ;
# - - - - -

    $PATH_SNDATA_SIM = "$SNDATA_ROOT/SIM" ;
    $USER_PATH_SNDATA_SIM = 0 ;

    $MISC_SDIR   = "misc";

    if ( $SNANA_DIR eq "" ) {
	$MSGERR[0] = "SNANA_DIR is not defined for HOST=$HOST ";
	sntools::FATAL_ERROR(@MSGERR);
    }


    $SSH_NNODE = $BATCH_NCORE = 0 ;
    $BATCH_TEMPLATE = "";

    $DOGEN_SNIa  = 0 ;
    $DOGEN_NONIa = 0 ;

    $OPT_ABORT = 0 ;
    $OPT_WARN  = 1 ;
    $OPT_QUIET = 2 ;

    $currentDir = `pwd`; 
    $currentDir =~ s/\s+$// ; 

    $LOGDIR         = "" ;
    $TOPDIR_SIMLOGS = "$currentDir"; 

    $INODE_GLOBAL = 0 ;

    $KILLJOBS_FLAG = 0 ;
    $DEBUG_GENSTAT = 0 ;
    $PROMPT_FLAG   = 1 ;

    @GENMODEL_NGENTOT  = ( ) ;
    
    @NVAR_SIMGEN_DUMP       = () ; 
    $CONVERT_SIMGEN_DUMP    = "" ;
    $NCLASS_SIMGEN_DUMP      = 0 ;
    $NCLASS_SIMGEN_DUMPALL   = 0 ;

# get max CID from source code (enabled Mar 9, 2012)
    my ($tmpFile, $key, @tmp);
    $tmpFile = "${SNANA_DIR}/src/snlc_sim.h";
    $key = "MXCID_SIM";
    @tmp = sntools::parse_line($tmpFile, 1, $key, $OPT_ABORT) ;
    $MXCID_SIM  = $tmp[0];
    $MXCID_SIM  =~ s/\s+$// ;  # remove trailing blanks

    $LAST_NGEN  = 0 ;
    $CIDOFF     = 0 ;
    $CIDRAN_MAX = 0 ;
    $CIDRAN_OFF = 0 ;

    @TOTAL_STRING = ();

    $DOSKIP_DUPLICATE_SIMJOBS = 1; 
    $NMODEL_DUPL_TOT = 0 ;

    $NABORT = 0 ;

    return ;

}   # end of init_SIMGEN

# ========================================================
sub parse_arg() {

    # parse command-line arguments.

    my ($arg, $argNext, $i, @USEARG );

    
    $NARG = scalar(@ARGV);
    if ( $NARG < 1 ) {
	$MSGERR[0] = "Must give master input filename as argument";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $SNMIX_INFILE_MASTER = $ARGV[0];
    $FASTFAC        = 1 ;
    $NODEINDX       = 0 ;
    $SUFFIX_TEMP    = "TEMP" ;
    $DO_SSH         = 0 ;
    $DO_BATCH       = 0 ;
    $SUBMIT_FLAG    = 1;
    $PROMPT_FLAG    = 1;

    for ( $i = 1; $i < $NARG ; $i++ ) { $USEARG[$i]=0; }

    for ( $i = 1; $i < $NARG ; $i++ ) {
	$arg     = $ARGV[$i];
	$argNext = $ARGV[$i+1] ; 

	if ( $arg eq "FAST10"    ) { $FASTFAC = 10  ;  $USEARG[$i]=1; }
	if ( $arg eq "FAST100"   ) { $FASTFAC = 100 ;  $USEARG[$i]=1; }
	if ( $arg eq "KILL"      ) { $KILLJOBS_FLAG=1; $USEARG[$i]=1; }	
	if ( $arg eq "GENSTAT"   ) { $DEBUG_GENSTAT=1; $USEARG[$i]=1; }	
	if ( $arg eq "NOSUBMIT"  ) { $SUBMIT_FLAG = 0; $USEARG[$i]=1; }
	if ( $arg eq "NOPROMPT"  ) { $PROMPT_FLAG = 0; $USEARG[$i]=1; }

	if ( "$arg" eq "-NODEINDX_SSH"   )  { 
	    $NODEINDX   = $argNext;  $DO_SSH = 1; 
	    $USEARG[$i] = $USEARG[$i+1] = 1 ; 
	}

	if ( "$arg" eq "-NODEINDX_BATCH" ) {
	    $NODEINDX  = $argNext; $DO_BATCH = 1;  
	    $USEARG[$i] = $USEARG[$i+1] = 1 ; 
	}

	# read time-stamp for suffix
	if ( "$arg" eq "-SUFFIX"  ) { 
	    $SUFFIX_TEMP  = $argNext; 
	    $USEARG[$i] = $USEARG[$i+1] = 1 ; 
	}

	if ( $arg eq "KICP" || $arg eq "kicp" ) 
	{ $BATCH_TEMPLATE = $BATCH_TEMPLATE_KICP ; }

    }

    # check that all args were used (Feb 2019)
    my $NERR=0;
    for ( $i = 1; $i < $NARG ; $i++ ) { 
	if ( $USEARG[$i] == 0 ) {
	    print " ERROR: unknown command-line arg: '$ARGV[$i]' \n";
	    $NERR++ ;
	}
    }
    if ( $NERR > 0 ) {
	$MSGERR[0] = "$NERR unknown command-line args.";
	$MSGERR[1] = "Scroll up to see list." ;
	sntools::FATAL_ERROR(@MSGERR);
    }
    
    # - - - - - 
    if ( $FASTFAC > 1 ) {
	print " Generate x$FASTFAC fewer SN to quick test. \n";
    }
    if ( $NODEINDX > 0  ) {
	print " Process SIMJOB(s) for nodeindx $NODEINDX \n";
    }


} # end of parse_arg

# ========================================================
sub parse_inFile_master() {

# check that file exists before parsing
# Jan 8 2019: echo INCLUDE files to allow for ENV

    my $inFile = "$SNMIX_INFILE_MASTER";
    my $cerr = "ERROR: Cannot open input file '$inFile'" ;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my ($key, @wdlist, $model,$imodel,$jtmp, $m, $iver, $ijob, $NTMP);
    my ($tmpFile1, $tmpFile2, @TMPNODES, $tmpLine, @tmp, $tmp);
    my (@CONTENTS_EXCLUDE_GENVERSION );

    $NSIMTYPE = 0 ;

    print "\n Parse master input file: \n";
    print "\t $inFile \n\n" ;

    sntools::loadArray_fromFile($inFile, \@CONTENTS_INFILE_MASTER);

    # read contents again, excluding GENVERSION info;
    # used below to read global SIMGEN_INFILE_NONIa[Ia] keys
    sntools::loadArray_excludeLines($inFile, 
				    "GENVERSION:", "ENDLIST_GENVERSION:",
				    \@CONTENTS_EXCLUDE_GENVERSION );


    $GENPREFIX                = "MIX" ;
    @SIMGEN_INFILE_GLOBAL     = () ;
    $INPUT_FILE_INCLUDE_Ia    = "" ;
    $INPUT_FILE_INCLUDE_NONIa = "" ;
    @GENOPT_FILE_INCLUDE      = () ;
    $INPUT_FILE_ZVAR_Ia    = "" ;
    $INPUT_FILE_ZVAR_NON1A = "" ;
    @INPUT_FILE_NON1AGRID  = () ;

    $SEARCHEFF_SPEC_FILE   = "" ;
    $HOSTNOISE_FILE        = "" ;
    $ZPHOTEFF_FILE         = "" ;
    $GENFILTERS = '' ;
    $NFILT = 0 ;
    $CLEANUP_FLAG    = 1  ; # remove TEMP directories
    $GENOPT_GLOBAL = '' ;
    $SNANA_LOGIN_SETUP = "" ;
    $RESET_CIDOFF = 0 ;  # change default flag to 0 (was 1) - July 2016
    $CIDRAN_MIN = 0 ;

    $WRFLAG_FITS = 0 ;
    $WRFLAG_TEXT = 0 ;
    $FORMAT_MASK = 0 ;

    $WRFLAG_CIDRAN = 0;

    # ---- start parsing

    $DONE_STAMP      = '' ;
    $DONE_STAMP_FLAG = 0;
    $key = "DONE_STAMP:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) {
	$DONE_STAMP      = "$tmp[0]" ;
	$DONE_STAMP_FLAG = 1 ;
	print " DONE_STAMP  = $DONE_STAMP \n"; 
#	if ( -e $DONE_STAMP )  { qx(rm $DONE_STAMP) ; }
    }

    $key = "NODELIST:";
# xxx    @TMPNODES = sntools::parse_line($inFile, 99, $key, $OPT_QUIET ) ;
    @TMPNODES = sntools::parse_array($key,99,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    @SSH_NODELIST = ();
    foreach $tmpLine ( @TMPNODES ) {
        @tmp = split(/\s+/,$tmpLine) ;
        @SSH_NODELIST =  ( @SSH_NODELIST , @tmp );
	$SSH_NNODE = scalar(@SSH_NODELIST);
	$NCPU      = $SSH_NNODE ;
    }

    $key = "JOBNAME_SIM:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) { 
	$JOBNAME_SIM = "$tmp[0]" ;  
	print " JOBNAME_SIM -> $JOBNAME_SIM \n";
    }
    $JOBNAME_SIM_FULLPATH = `which $JOBNAME_SIM` ;
    $JOBNAME_SIM_FULLPATH =~ s/\s+$// ;   # trim trailing whitespace

    $key = "PATH_SNDATA_SIM:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) { 
	$PATH_SNDATA_SIM = "$tmp[0]" ;  
        $PATH_SNDATA_SIM = qx(echo $PATH_SNDATA_SIM);
        $PATH_SNDATA_SIM =~ s/\s+$// ;   # trim trailing whitespace
	$USER_PATH_SNDATA_SIM = 1 ;
    }

    # check for batch command
    $key   = "BATCH_INFO:" ;
    @tmp   = sntools::parse_array($key,3,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) {
	@wdlist         = split(/\s+/,$tmp[0]) ;
	$BATCH_COMMAND  = $wdlist[0] ;
	if (length($BATCH_TEMPLATE)==0 ) { $BATCH_TEMPLATE = $wdlist[1] ; }
	$BATCH_NCORE    = $wdlist[2] ;

	# allow for ENV (May 22 2017) 
        $BATCH_TEMPLATE = qx(echo $BATCH_TEMPLATE);
        $BATCH_TEMPLATE  =~ s/\s+$// ;   # trim trailing whitespace

	$NCPU = $BATCH_NCORE ;  # use NCPU for both SSH and BATCH
	@SSH_NODELIST = ( "Chosen by batch" );

	my $BCMD = "$BATCH_COMMAND $BATCH_TEMPLATE";
	my $B3   = "$BATCH_NCORE cores";
	print " BATCH_COMMAND: $BCMD ($B3) \n";	
    }

    $key   = "BATCH_MEM:" ;  # MB
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) { $BATCH_MEM = $tmp[0] ; }

    if ( $SSH_NNODE > 0  && $BATCH_NCORE > 0 ) {
	$MSGERR[0] = "Cannot specify both NODELIST and BATCH_INFO keys.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }


    # now set CPULIST_INDX that is unique even if the node-name is repeated
    for ( $jtmp=1; $jtmp <= $NCPU ; $jtmp++ ) {
	@CPULIST_INDX = ( @CPULIST_INDX , $jtmp );
    }
    print " SSH_NODELIST(${NCPU}) = @SSH_NODELIST \n" ;
    print " CPULIST_INDX(${NCPU}) = @CPULIST_INDX \n" ;


    $key = "SNANA_LOGIN_SETUP:" ;
    @tmp = sntools::parse_array($key,99,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) {
	$SNANA_LOGIN_SETUP = "$tmp[0]" ;
	print " SNANA_LOGIN_SETUP = '${SNANA_LOGIN_SETUP}' \n" ;
    }


    $key = "TOPDIR_SIMLOGS:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) {
	$TOPDIR_SIMLOGS = $tmp[0] ;
	$TOPDIR_SIMLOGS = qx(echo $TOPDIR_SIMLOGS) ; # unpack ENV
	$TOPDIR_SIMLOGS =~ s/\s+$// ;   # trim trailing whitespace  
	print " TOPDIR_SIMLOGS = '${TOPDIR_SIMLOGS}' \n" ;
    }

    $key = "LOGDIR:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) {
	$LOGDIR = $tmp[0] ;
	$LOGDIR = qx(echo $LOGDIR) ; # unpack ENV
	$LOGDIR =~ s/\s+$// ;        # trim trailing whitespace  
	print " LOGDIR = '${LOGDIR}' \n" ;
    }

    $SNANA_MODELPATH  = "" ;
    $ENVDEF_MODELPATH = "" ;
    $key = "SNANA_MODELPATH:" ;
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) {
	$SNANA_MODELPATH  = "$tmp[0]" ;
	$SHELL_NAME       = sntools::shellName();
	$ENVDEF_MODELPATH = 
	    sntools::setEnvString($SHELL_NAME, 
				  "SNANA_MODELPATH", $SNANA_MODELPATH);
    }

    $key   = "CONVERT_SIMGEN_DUMP:" ;  # arg is ROOT or HBOOK
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) { $CONVERT_SIMGEN_DUMP = $tmp[0] ; }

    $key = "DOSKIP_DUPLICATE_SIMJOBS:";
    @tmp   = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) { $DOSKIP_DUPLICATE_SIMJOBS=$tmp[0] ; }

    $key = "SIMGEN_INFILE_Ia:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,
				@CONTENTS_EXCLUDE_GENVERSION );
    if ( scalar(@tmp) == 0 ) {
	# try alternate key
	$key = "SIMGEN_INFILE_SNIa:";
	@tmp = sntools::parse_array($key,1,$OPT_QUIET,
				    @CONTENTS_EXCLUDE_GENVERSION );
    }
    if ( scalar(@tmp) > 0 ) { 
	$tmpFile1 = $tmp[0];
	$tmpFile1 = qx(echo $tmpFile1); # unpack ENV
	$tmpFile1 =~ s/\s+$// ;   # trim trailing whitespace  

	$mIa=$NGENMODEL_GLOBAL ;  
	$DOGEN_SNIa=1 ; $NSIMTYPE++ ;  $NGENMODEL_GLOBAL++ ;
	$SIMGEN_INFILE_GLOBAL[$mIa]   = "$tmpFile1" ;
	$GENMODEL_NAME_GLOBAL[$mIa]   = "$tmpFile1" ;
	$GENMODEL_CLASS_GLOBAL[$mIa]  = "SNIa" ;   	
	
	print " SIMGEN_INFILE_Ia     = $tmpFile1  \n" ;
	unless (-e $tmpFile1 ) {
	    $MSGERR[0] = "'$tmpFile1' does not exist";
	    $MSGERR[1] = "Check argument of SIMGEN_INFILE_Ia:" ;
	    $MSGERR[1] = " or   argument of SIMGEN_INFILE_SNIa:" ;
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);  
	}
    }

    # note that multiple files can be listed after this key
    $key = "SIMGEN_INFILE_NONIa:"; $m = $NGENMODEL_GLOBAL;
    @tmp = sntools::parse_array($key,99,$OPT_QUIET,
				@CONTENTS_EXCLUDE_GENVERSION );
    $NTMP = scalar(@tmp); 
    if ( $NTMP > 0 ) { $NSIMTYPE++;  $DOGEN_NONIa = 1; }
    foreach $tmpLine ( @tmp ) {
	@wdlist = split(/\s+/,$tmpLine) ;
	foreach $tmpFile1 ( @wdlist) {
	    $tmpFile1 = qx(echo $tmpFile1); 
	    $tmpFile1 =~ s/\s+$// ;   # trim trailing whitespace  
	    print " SIMGEN_INFILE_NONIa  = $tmpFile1  \n" ;
	    unless (-e $tmpFile1 ) {
		$MSGERR[0] = "'$tmpFile1' does not exist";
		$MSGERR[1] = "Check argument of SIMGEN_INFILE_NONIa:" ;
		sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	    }
	    $GENMODEL_NAME_GLOBAL[$m]   = "$tmpFile1" ;
	    $GENMODEL_CLASS_GLOBAL[$m]  = "NONIaMODEL$m" ;
	    $SIMGEN_INFILE_GLOBAL[$m] = $tmpFile1 ;
	    $m++ ;  $NGENMODEL_GLOBAL++ ;
	}
    }



    $key = "RESET_CIDOFF:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) { $RESET_CIDOFF = $tmp[0]; }

    $key = "CIDRAN_MIN:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) { 
	$CIDRAN_MIN = $tmp[0]; 
	if ( $CIDRAN_MIN > 0 && $RESET_CIDOFF!=2 ) {
	    $MSGERR[0] = "If CIDRAN_MIN>0, must set RESET_CIDOFF: 2" ;
	    $MSGERR[1] = "Current CIDRAN_MIN: $CIDRAN_MIN ";
	    $MSGERR[2] = "Current RESET_CIDOFF: $RESET_CIDOFF" ;
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}

    }

    # ---------------------------------
    # parse keys that are required if simulating both Ia and CC,
    # but optional if simulating only 1 type.

    my ($OPT_PARSE, $SIMGEN_INCLUDE); 
    
    if ( $NSIMTYPE == 1 ) {
	$OPT_PARSE     = $OPT_QUIET ;
	sntools::loadArray_fromFile($SIMGEN_INFILE_GLOBAL[0],
				    \@CONTENTS_INFILE_SIMGEN);

        # parse include file into temporary variable for required keys
	$key = "INPUT_FILE_INCLUDE:";    
	@tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_SIMGEN);
	if ( scalar(@tmp) > 0 ) {
	    $SIMGEN_INCLUDE = $tmp[0];
	    $SIMGEN_INCLUDE = qx(echo $SIMGEN_INCLUDE); # unpack ENV
	    $SIMGEN_INCLUDE =~ s/\s+$// ;   # trim trailing whitespace  
	    sntools::loadArray_fromFile($SIMGEN_INCLUDE,
					\@CONTENTS_INFILE_INCLUDE);
	}
    }
    else 
    { $OPT_PARSE = $OPT_ABORT ; } 

    $SIMARG_FORMAT  = "" ;
    $SIMARG_RANSEED[0] = "" ;

    $key = "ZRANGE:" ; 
    @tmp = sntools::parse_array($key,2,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) {
	$MSGERR[0] = "$key  key is obsolete." ;
	$MSGERR[1] = "Remove 'ZRANGE' key from simgen-master file";
	$MSGERR[2] = "and instead define global redshift range with";
	$MSGERR[3] = "  GENOPT_GLOBAL: GENRANGE_REDSHIFT $tmp[0] $tmp[1]";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);	
    }


    $key = "H0:" ; 
    @tmp = sntools::parse_array($key,2,$OPT_QUIET,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) {
	$MSGERR[0] = "$key  key is obsolete." ;
	$MSGERR[1] = "Remove this key from simgen-master file";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);	
    }

    &parse_GENOPT_GLOBAL();
    &parse_GENVERSION();
    &parse_RANSEED();
    &parse_FORMAT_MASK();  
    &check_FORMAT_MASK();

    # --------------

    $key = "NGEN_UNIT:";
    $NGEN_UNIT_VAL = -9.0 ;
    @tmp = sntools::parse_array($key,2,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) {
	$NGEN_UNIT_VAL  = $tmp[0] ;
	$NGEN_UNIT_NAME = $tmp[1] ;
	print " NGEN_UNIT  = $NGEN_UNIT_VAL  $NGEN_UNIT_NAME \n" ;
    }

    $key = "GENPREFIX:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) {
	if ( "$GENPREFIX" ne "MIX" ) {
	    $MSGERR[0] = "GENPREFIX defined twice: $GENPREFIX and $tmp[0]" ;
	    $MSGERR[1] = "Only one GENPREFIX declaration allowed in " . 
		"master-input.";
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);	    
	}
	$GENPREFIX = $tmp[0]; 
    }
    print " GENPREFIX = $GENPREFIX  \n" ;

    $key = "CLEANUP_FLAG:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) {
	$CLEANUP_FLAG = $tmp[0];
	print " CLEANUP_FLAG     = $CLEANUP_FLAG  \n" ;
    }

    # get other optional info from sim-input file 
    &parse_INFILE_Ia();

    for($iver=1; $iver<=$NGENVERSION; $iver++ ) {
	for($m=0; $m < $NGENMODEL[$iver]; $m++ ) 
	{ &parse_INFILE_NONIa($iver,$m); }
    }

    if ( $NSIMTYPE == 0 ) {
	$MSGERR[0] = " Must specify sim-input file for Ia and/or NON1A." ;
	$MSGERR[1] = " Add 1 or more of the following keys in $SNMIX_INFILE_MASTER :" ;
	$MSGERR[2] = "  SIMGEN_INFILE_Ia: <inFile> " ;
	$MSGERR[3] = "  SIMGEN_INFILE_SNIa: <inFile> " ;
	$MSGERR[4] = "  SIMGEN_INFILE_NONIa: <inFile> " ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    print " DOGEN[SNIa,NONIa] = $DOGEN_SNIa,$DOGEN_NONIa \n";
    &check_SIMGEN_DUMP();
    
    if ( $USER_PATH_SNDATA_SIM ) {
	print " PATH_SNDATA_SIM -> $PATH_SNDATA_SIM \n";
	if ( !(-d $PATH_SNDATA_SIM) ) {
	    $MSGERR[0] = "PATH_SNDATA_SIM does not exist:" ;
	    $MSGERR[1] = "  $PATH_SNDATA_SIM" ;
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

    return ;

} # end of parse_inFile_master


# ===========================
sub parse_INFILE_Ia {

    my $inFile = $SIMGEN_INFILE_GLOBAL[$mIa] ;
    my( $key, @tmp, @wdlist );

    if ( $inFile eq "" ) { return ; } 

    print "\t Getting optional info from sim-input file ... \n";

    sntools::loadArray_fromFile($inFile,\@CONTENTS_INFILE_Ia);
	
    $key = "INPUT_FILE_INCLUDE:";    
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 ) {
	$INPUT_FILE_INCLUDE_Ia = "$tmp[0]" ;
	$INPUT_FILE_INCLUDE_Ia = qx(echo $INPUT_FILE_INCLUDE_Ia);
#	print "  INPUT_FILE_INCLUDE_Ia = $INPUT_FILE_INCLUDE_Ia \n";
	my (@INCLUDE);
	sntools::loadArray_fromFile($INPUT_FILE_INCLUDE_Ia,\@INCLUDE);
	@CONTENTS_INFILE_Ia = (@CONTENTS_INFILE_Ia, @INCLUDE);
    }

    # Aug 15 2020: check GENOPT_GLOBAL for include file
    if ( length($GENOPT_GLOBAL) > 0  ) {
	@wdlist  = split(/\s+/,$GENOPT_GLOBAL) ;
	$key = "INPUT_FILE_INCLUDE";
	@tmp = sntools::parse_array($key,1, $OPT_QUIET, @wdlist);
	if ( scalar(@tmp) > 0 ) {
	    $INPUT_FILE_INCLUDE_Ia = "$tmp[0]" ;
	    $INPUT_FILE_INCLUDE_Ia = qx(echo $INPUT_FILE_INCLUDE_Ia);
	    #print "  xxx INPUT_FILE_INCLUDE_Ia = $INPUT_FILE_INCLUDE_Ia \n";  
	    my (@INCLUDE);
	    sntools::loadArray_fromFile($INPUT_FILE_INCLUDE_Ia,\@INCLUDE);
	    @CONTENTS_INFILE_Ia = (@CONTENTS_INFILE_Ia, @INCLUDE);
	}
    }


    $key = "GENFILTERS:";
    @tmp = sntools::parse_array($key,1,$OPT_ABORT,@CONTENTS_INFILE_Ia);
    $GENFILTERS = $tmp[0];
#    print " GENFILTERS       = $GENFILTERS  \n" ;
    $NFILT = length($GENFILTERS);
	
    $key = "ZVARIATION_FILE:";    
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 ) {
	$INPUT_FILE_ZVAR_Ia = "$tmp[0]" ;
	print "  INPUT_FILE_ZVAR_Ia = $INPUT_FILE_ZVAR_Ia \n";
    }
    
    $key = "SEARCHEFF_SPEC_FILE:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 ) {
	$SEARCHEFF_SPEC_FILE = "$tmp[0]" ;
	print "  SEARCHEFF_SPEC_FILE = $SEARCHEFF_SPEC_FILE \n";
    }

    $key = "HOSTNOISE_FILE:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 ) {
	$HOSTNOISE_FILE = "$tmp[0]";
	print "  HOSTNOISE_FILE = $HOSTNOISE_FILE \n";
    }

    $key = "HOSTLIB_ZPHOTEFF_FILE:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 ) {
	$ZPHOTEFF_FILE = "$tmp[0]";
	print "  ZPHOTEFF_FILE = $ZPHOTEFF_FILE \n";
    }

    # check user CIDOFF unless random CIDs are requested.
    $key = "CIDOFF:" ;
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 && $WRFLAG_CIDRAN == 0 ) {
	$CIDOFF = $tmp[0];
	print "  CIDOFF = $CIDOFF   (RESET_CIDOFF = $RESET_CIDOFF) \n";
    }

    # check number of SIMGEN_DUMP variables.
    $key = "SIMGEN_DUMP:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 ) {
	$NVAR_SIMGEN_DUMP[$mIa] = $tmp[0];
	$NCLASS_SIMGEN_DUMP++ ; # write events that pass trigger+cuts
    }
    $key = "SIMGEN_DUMPALL:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 ) {
	$NVAR_SIMGEN_DUMP[$mIa] = $tmp[0];
	$NCLASS_SIMGEN_DUMPALL++ ; # write ALL generated events
    }

    $key = "GENMODEL:" ;	
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    $GENMODEL_NAME_GLOBAL[0] = $tmp[0] ;

    $key = "PATH_SNDATA_SIM:" ;	
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_Ia);
    if ( scalar(@tmp) > 0 ) { 
	$PATH_SNDATA_SIM = $tmp[0] ; 
	$PATH_SNDATA_SIM = qx(echo $PATH_SNDATA_SIM);
	$PATH_SNDATA_SIM =~ s/\s+$// ;   # trim trailing whitespace
	$USER_PATH_SNDATA_SIM = 1 ;
    }
    
    return ; 

} # end of parse_INFILE_Ia 


# ===========================
sub parse_INFILE_NONIa {
    my($iver,$m) = @_; 

    my $inFile = $SIMGEN_INFILE[$iver][$m];
    my( $key, @tmp, @wdlist );

    if ( $m == $mIa ) { return ; }
    @CONTENTS_INFILE_NONIa = () ;

    sntools::loadArray_fromFile($inFile,\@CONTENTS_INFILE_NONIa);

    $key = "INPUT_FILE_INCLUDE:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_NONIa);
    if ( scalar(@tmp) > 0 ) {
	$INPUT_FILE_INCLUDE_NONIa = $tmp[0];
	$INPUT_FILE_INCLUDE_NONIa = qx(echo $INPUT_FILE_INCLUDE_NONIa);
	my (@INCLUDE) ;
	sntools::loadArray_fromFile($INPUT_FILE_INCLUDE_NONIa,\@INCLUDE);
	@CONTENTS_INFILE_NONIa = (@CONTENTS_INFILE_NONIa, @INCLUDE);
    }

    # Aug 15 2020: check GENOPT_GLOBAL for include file
    if ( length($GENOPT_GLOBAL) > 0  ) {
	@wdlist         = split(/\s+/,$GENOPT_GLOBAL) ;
	$key = "INPUT_FILE_INCLUDE";
	@tmp = sntools::parse_array($key,1, $OPT_QUIET, @wdlist);
	if ( scalar(@tmp) > 0 ) {
	    $INPUT_FILE_INCLUDE_NONIa = $tmp[0];
	    $INPUT_FILE_INCLUDE_NONIa = qx(echo $INPUT_FILE_INCLUDE_NONIa);
	    my (@INCLUDE) ;
	    sntools::loadArray_fromFile($INPUT_FILE_INCLUDE_NONIa,\@INCLUDE);
	    @CONTENTS_INFILE_NONIa = (@CONTENTS_INFILE_NONIa, @INCLUDE);
	}
    }

    $key = "GENFILTERS:";
    @tmp = sntools::parse_array($key,1,$OPT_ABORT,@CONTENTS_INFILE_NONIa);
    if ( $NFILT > 0  &&  "$tmp[0]" ne "$GENFILTERS" ) {
	$MSGERR[0] = "GENFILTERS = '${tmp[0]}' in $inFile" ;
	$MSGERR[1] = "but GENFILTERS = '${GENFILTERS}' in Ia file.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }
    if ( $NFILT == 0 ) {
	$GENFILTERS = $tmp[0];
	print " GENFILTERS       = $GENFILTERS  \n" ;
	$NFILT = length($GENFILTERS);
    }

    $key = "ZVARIATION_FILE:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_NONIa);
    if ( scalar(@tmp) > 0 ) {
	$INPUT_FILE_ZVAR_NON1A = $tmp[0];
	print "  INPUT_FILE_ZVAR_NON1A = $INPUT_FILE_ZVAR_NON1A \n";
    }

    # check number of SIMGEN_DUMP variables.
    $key = "SIMGEN_DUMP:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_NONIa);
    if ( scalar(@tmp) > 0 ) {
	$NVAR_SIMGEN_DUMP[$m] = $tmp[0];
	$NCLASS_SIMGEN_DUMP++ ; # write events that pass trigger+cuts
    }
    $key = "SIMGEN_DUMPALL:";
    @tmp = sntools::parse_array($key,1,$OPT_QUIET,@CONTENTS_INFILE_NONIa);
    if ( scalar(@tmp) > 0 ) {
	$NVAR_SIMGEN_DUMP[$m] = $tmp[0];
	$NCLASS_SIMGEN_DUMPALL++ ; # write ALL generated events.
    }


    # store name of nonIa model (NONIa, FIXMAG, etc ...)
    $key = "GENMODEL:" ;	
    @tmp = sntools::parse_array($key,2,$OPT_QUIET,@CONTENTS_INFILE_NONIa);
    @wdlist         = split(/\s+/,$tmp[0]) ;
    my $tmpModel = $wdlist[0];
    my $tmpArg   = $wdlist[1];
    $GENMODEL_NAME_GLOBAL[$m]       = $tmpModel ;
    
    my $N1  = $NVAR_SIMGEN_DUMP[$mIa] ;
    my $N2  = $NVAR_SIMGEN_DUMP[$m] ;
    if ( $N1 > 0 && $N2 > 0 &&  $N1 != $N2 ) {
	$MSGERR[0] = "NVAR_SIMGEN_DUMP(Ia,NONIa) = $N1, $N2" ;
	$MSGERR[1] = "but they must be equal." ;
	$MSGERR[2] = "Check SIMGEN_DUMP key in" ;
	$MSGERR[3] = "  $SIMGEN_INFILE_GLOBAL[$mIa]  and" ;
	$MSGERR[4] = "  $SIMGEN_INFILE_GLOBAL[$m]  ." ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }
    if ( $GENMODEL_NAME_GLOBAL[$m] == "NON1AGRID" ) 
    { $INPUT_FILE_NON1AGRID[$m] = $tmpArg ; }  # maybe get rid of this?

    &syncVerify_files($m);

} # end parse_INFILE_NONIa


# ============================================
sub check_SIMGEN_DUMP {

    print " NCLASS_SIMGEN_DUMP,ALL = " . 
	"$NCLASS_SIMGEN_DUMP, $NCLASS_SIMGEN_DUMPALL \n" ;

    if ( $NCLASS_SIMGEN_DUMP || $NCLASS_SIMGEN_DUMPALL ) 
    { $DO_SIMGEN_DUMP = 1; }
    else
    { $DO_SIMGEN_DUMP = 0; }
   
} # end check_SIMGEN_DUMP

# ============================================
sub syncVerify_files {

    my ($m) = @_ ;

    # - - - - - - - - - - - - - - - - - - - - - - - - -
    #  for Ia+CC mix, make make sure certain gen-params
    #  are the same in each file: z-range, cosmo params, 
    #  peak-MJD range, etc ...
    #  Do check only for NONIA model; not for FIXMAG
    # - - - - - - - - - - - - - - - - - - - - - - - - -
    
    my ($ISCC,$G5);
    $G5   = substr($GENMODEL_NAME_GLOBAL[$m],0,5);  # e.g., allow NON1ASED
    $ISCC = ( $G5 ne "FIXMAG" );

    if ( $NSIMTYPE == 2 && $ISCC ) {
	my ($found);  
	$found =  &syncVerify_GENPAR("GENRANGE_PEAKMJD:",  1, $OPT_ABORT);
	$found =  &syncVerify_GENPAR("GENRANGE_REDSHIFT:", 2, $OPT_ABORT);

	$found =  &syncVerify_GENPAR("SOLID_ANGLE:",       1, $OPT_QUIET);
	if ( $found == 0 ) {
	    $found =  &syncVerify_GENPAR("GENRANGE_RA:",   2, $OPT_ABORT);
	    $found =  &syncVerify_GENPAR("GENRANGE_DECL:", 2, $OPT_ABORT);
	}
    }

} # end syncVerify

# ============================================
sub syncVerify_GENPAR {
    my ($KEY, $NVAL, $OPT_MISSING) = @_ ;

    # ABORT if Ia and NONIa values are differnt.
    # If OPT_MISSING = OPT_QUIET, then allow missing keys,
    # but always abort if values are different or if only
    # only one file is missing the key.
    # Return 1 of KEY is found; 0 otherwise.

    my (@VAL_Ia, @VAL_NONIa, $NVAL_Ia, $NVAL_NONIa, $ival, $VAL1, $VAL2);

    @VAL_Ia = sntools::parse_array($KEY,$NVAL,$OPT_MISSING,
				   @CONTENTS_INFILE_Ia);

    @VAL_NONIa = sntools::parse_array($KEY,$NVAL,$OPT_MISSING,
				      @CONTENTS_INFILE_NONIa);

    $NVAL_Ia    = scalar(@VAL_Ia);
    $NVAL_NONIa = scalar(@VAL_NONIa);

    if ( $OPT_MISSING == $OPT_QUIET ) {
	if ( $NVAL_Ia == 0 && $NVAL_NONIa == 0 ) { return 0 ; }
    }

    for($ival=0; $ival < $NVAL; $ival++ ) { 
	$VAL1 = $VAL_Ia[$ival] ;
	$VAL2 = $VAL_NONIa[$ival] ;

	if ( "$VAL1" ne "$VAL2" ) {	    
	    $MSGERR[0] = "Mis-matched '$KEY' keys not allowed: ";
	    $MSGERR[1] = "Ia-input file has " ;
	    $MSGERR[2] = "   $KEY @VAL_Ia \n";
	    $MSGERR[3] = "NONIa-input file has " ;
	    $MSGERR[4] = "   $KEY @VAL_NONIa \n" ;
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}
    }

    return 1 ;

} # end of syncVerify


# =====================================
sub parse_GENOPT_GLOBAL {

    my ($key, @tmp, $tmpLine, @values, $NVAL);

    $key = "GENOPT_GLOBAL:" ;
    @tmp = sntools::parse_array($key,99,$OPT_QUIET,@CONTENTS_INFILE_MASTER);
    if ( scalar(@tmp) > 0 ) {
	foreach $tmpLine ( @tmp ) {
	    if ( index($tmpLine,"#") > 0 ) {
		die "\n ERROR: comment not allowed after GENOPT_GLOBAL key\n" .
		    "   Invalid line: '$tmpLine' \n";
	    }
	    $GENOPT_GLOBAL = "$GENOPT_GLOBAL  $tmpLine" ;
	}
    }

    # check for special characters that need backslash (9/28 2017) 
    $GENOPT_GLOBAL =~ s/\(/\\(/g ;  # ( --> \(
    $GENOPT_GLOBAL =~ s/\)/\\)/g ;  # ) --> \)

    # check for FORMAT_MASK here, then later check master-input
    # file for <FORMAT_MASK: MASK>

    @values = sntools::parse_value_after_key("FORMAT_MASK",$GENOPT_GLOBAL);
    $NVAL = scalar(@values);
    if ( $NVAL == 1 ) {	$FORMAT_MASK = $values[0]; }

    # same for GENPREFIX
    @values = sntools::parse_value_after_key("GENPREFIX",$GENOPT_GLOBAL);
    $NVAL = scalar(@values);
    if ( $NVAL == 1 ) {	
	$GENPREFIX = $values[0]; 
	# remove GENPREFIX from GENOPT_GLOBAL because a modified
	# GENPREFIX will be automatically added to account for each core.
	print "\t (remove GENPREFIX $GENPREFIX from GENOPT_GLOBAL) \n";
	$GENOPT_GLOBAL =~ s/GENPREFIX//g ;
	$GENOPT_GLOBAL =~ s/$GENPREFIX//g ;

    }

    # print at end of function in case some parts are removed.
    print " GENOPT_GLOBAL  = '$GENOPT_GLOBAL' \n" ;    

} # end parse_GENOPT_GLOBAL

# =========================================
sub parse_GENVERSION {

    # parse "GENVERSION:"  and "GENOPT:" keys.
    # Allow comments after hash
    # GENOPT can be applied to specific models with
    #  GENOPT(SNIa):  bla bla
    #  GENOPT(NONIa):  bla bla
    #  GENOPT([inFile-subtring]):  bla bla
    #
    # Note that pad spaces don't matter, and thus "BLA bla" 
    # and "BLA   bla" arguments are considered duplicates.
    #
    # Feb 11 2019: check for SIMGEN_INFILE_NONIa[Ia]
    # Apr 18 2019: set DOGEN_NONIa[SNIa] flags (bug fix)
    # May 21 2020: increment NSIMTYPE if first time with Ia or NONIa


    my (@INFILE_LINES, @USE_GENOPT, @GENOPT_LINES, @NGENOPT, $NOPT );
    my ($VAL,$KEY, $key, $NWD, $GENOPT, $CLASS, $NAME );
    my ($iver, $iv, $jhash, $m, @mlist, $iwd, $narg, $opt, $iflag);
    my (@wdlist, $tmpLine, @argList, @INDX, @WDLIST, $indx );
    my ($keyArg, $tmpArg, $INFILE );
    my (@INFILE_NONIA_OVERRIDE);

    # need to read entire line ,not just sequence of words ->
    # re-read contents of master file.

    open  PTR_TEMP , "< ${SNMIX_INFILE_MASTER}" ;
    @INFILE_LINES = <PTR_TEMP> ;
    close PTR_TEMP ;

    $NGENVERSION        = 0;
    @GENVERSION_GENOPT  = ( ) ;
    @GENOPT_LINES       = ( ) ;

    foreach $tmpLine (@INFILE_LINES) {	

	$jhash = index($tmpLine,"#") ; 
	if ( $jhash > 0 ) { $tmpLine = substr($tmpLine,0,$jhash-1) ; }

	@wdlist   = split(/\s+/,$tmpLine) ;
	$KEY = $wdlist[0] ;
	$VAL = $wdlist[1] ;

	# store arguments of key, excluding key.
	@argList=(); $iwd=$narg=0; 
	foreach $tmpArg (@wdlist)  { 
	    if ( $iwd>0 ) { $argList[$narg]=$tmpArg; $narg++ ; }  
	    $iwd++; 
	}

	if ( $KEY eq "ENDLIST_GENVERSION:" )   { goto PROCESS_GENOPT ; }

	if ( $KEY eq "GENVERSION:" ) {	
	    
	    # abort if this GENVERSION name is already used
	    @INDX = 
		grep { $GENVERSION_NAME[$_] eq "$VAL" } 0 .. $#GENVERSION_NAME;
	    if ( scalar(@INDX) > 0 ) {
		$MSGERR[0] = "Control file GENVERSION specifier";
		$MSGERR[1] = "   GENVERSION: $VAL" ;
		$MSGERR[2] = "appears multiple times in";
		$MSGERR[3] = "$SNMIX_INFILE_MASTER ." ;
		$MSGERR[4] = "Each GENVERSION name must be unique.";
		sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	    }
	    $NGENVERSION++ ;  
	    $iver = $NGENVERSION ;
	    &store_SIMGEN_INFILE($iver, 0, \@argList);
	    $GENVERSION_NAME[$iver] = "$VAL" ;
	    $INFILE_NONIA_OVERRIDE[$iver] = 0 ;
	    $NGENOPT[$iver] = 0 ;

	    # do not allow dot in GENVERION (Nov 19 2015).
	    # sim runs fine with dot, but dot causes problems in the 
	    # merging after split_and_fit.
	    if ( index($VAL,".") > 0 ) {
		$MSGERR[0] = "Dot (.) not allowed in GENVERSION name";
		$MSGERR[1] = "GENVERSION: $VAL";
		sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	    }
	}

	# store GENOPT lines and process later after 
	# reading all GENOPT and SIMGEN-INPUT line for each GENVERSION.
	if ( substr($KEY,0,6) eq "GENOPT" ) {
	    $NOPT    = $NGENOPT[$iver];
	    $GENOPT_LINES[$iver][$NOPT] = "$tmpLine" ;
	    $NGENOPT[$iver]++ ;
	}
   
	# check for GENVERSION-dependent sim-input files
	if ( $KEY eq "SIMGEN_INFILE_Ia:" || $KEY eq "SIMGEN_INFILE_SNIa:" ) { 
            if ( !$DOGEN_SNIa ) { $NSIMTYPE++ ; }
	    $DOGEN_SNIa = 1;
	    &store_SIMGEN_INFILE($iver, 1, \@argList);  
	}
	if ( $KEY eq "SIMGEN_INFILE_NONIa:" )  { 
            if ( !$DOGEN_NONIa ) { $NSIMTYPE++ ; }
	    $DOGEN_NONIa = 1;
	    $iflag = 2 + $INFILE_NONIA_OVERRIDE[$iver] ;
	    &store_SIMGEN_INFILE($iver, $iflag, \@argList);  
	    $INFILE_NONIA_OVERRIDE[$iver]++ ;
	}
	
    }  # end loop over tmpLines
   

# ---------------------------------------------------------------
# process GENOPT now that we have all of the SIMGEN-input files.
  PROCESS_GENOPT:

    for ( $iv=1; $iv <= $NGENVERSION; $iv++ ) {
	for($opt=0; $opt< $NGENOPT[$iv]; $opt++ ) {  
	    $tmpLine = "$GENOPT_LINES[$iv][$opt]" ;
	    $tmpLine =~ s/\s+$// ;   # trim trailing whitespace
	    &PROCESS_GENOPT_LINE($iv,$tmpLine);	    
	}  

	# append GENOPT_GLOBAL
	for($m=0; $m < $NGENMODEL[$iver]; $m++ ) {
	    my $GENOPT_TMP  = "$GENVERSION_GENOPT[$iv][$m] $GENOPT_GLOBAL";
	    $GENVERSION_GENOPT[$iv][$m] = "$GENOPT_TMP";
	}
	
	# read few things 
	&read_GENOPT_STRING($iv);
    }  


# ---------------------------------------------------------------
  GENVERSION_SUMMARY:

    if ( $NGENVERSION == 0 ) {
	$MSGERR[0] = "Found no GENVERSION keys in ";
	$MSGERR[1] = "$SNMIX_INFILE_MASTER";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    print "\n" ;
    print "\t ***** GENVERSION SUMMARY ***** \n" ;
    for ( $iv=1; $iv <= $NGENVERSION; $iv++ ) {

	print " GENVERSION      = '$GENVERSION_NAME[$iv]' \n" ;

	for($m=0; $m < $NGENMODEL[$iv]; $m++ ) {

	    # remove leading spaces
	    $GENVERSION_GENOPT[$iv][$m]   =~ s/^\s+// ;  
	
	    # remove trailing spaces
	    $GENVERSION_GENOPT[$iv][$m]   =~ s/\s+$// ;

	    $INFILE = $SIMGEN_INFILE[$iv][$m] ;
	    $CLASS  = $GENMODEL_CLASS[$iv][$m];  
	    $GENOPT = $GENVERSION_GENOPT[$iv][$m] ;

	    print "   GENOPT($INFILE) = $GENOPT  ($CLASS) \n" ;
	}
    }
  
    print "\n" ;
    $| = 1;        # flush stdout  


    return ;
   
} # end of parse_GENVERSION


# ======================================
sub PROCESS_GENOPT_LINE {
    my ($iver,$tmpLine) = @_ ;
    
    # catenate GENOPT argument for each sub-model.

#    print " xxx iver=$iver  tmpLine='$tmpLine' \n";
    my (@wdlist, $KEY, $VAL, $key, $keyArg, $tmpArg, $NWD, $m );
    my ($FOUND, $GENOPT, $MATCH );
    my (@USE_GENOPT);
    my @keyArgList_Ia    = ( "SN1a", "SNIa", "SN1A", "SNIA", 
			     "1a",   "Ia",   "1A",   "IA" );
    my @keyArgList_NONIa = ( "NONIa", "NON1a", "NONIA", "NON1A" );
    

    @wdlist  = split(/\s+/,$tmpLine) ;
    $KEY = $wdlist[0] ;
    $VAL = $wdlist[1] ;
    $NWD = scalar(@wdlist);

    # for GENOPT key, allow argument in () to indicate Ia or NONIa,
    # or substring in any SIMGEN_INFILE name. See top of function.
    # $key is "GENOPT:" and keyArg is whatever is in ().
    
    $key = "" ;
    if ( substr($KEY,0,6) eq "GENOPT" ) {
	($key,$keyArg) = sntools::extractStringOpt($KEY);
    }

    if ( $key ne "GENOPT:" ) { return; }

    $FOUND   = 0 ;	   
    $tmpArg  = "@wdlist[1 .. $NWD]" ;  # snlc_sim.exe args

    # check for special characters that need backslash (9/28 2017) 
    $tmpArg =~ s/\(/\\(/g ;  # ( --> \(
    $tmpArg =~ s/\)/\\)/g ;  # ) --> \)
    
    # default is tmpArg applies to all sub-models	   
    for($m=0; $m<$NGENMODEL[$iver]; $m++ ) { $USE_GENOPT[$m]=1; }

    # check optional argument inside ()
    if ( length($keyArg) > 1 ) { 
	if ( grep( /^$keyArg$/, @keyArgList_Ia ) ) {
	    $FOUND = 1 ;
	    # for GENOPT(Ia), turn off options for NON1A models
	    for($m=0; $m<$NGENMODEL[$iver] ; $m++ ) 
	    { if ( $m != $mIa ) { $USE_GENOPT[$m] = 0; } }
	}	       
	elsif ( grep( /^$keyArg$/, @keyArgList_NONIa ) ) {
	    $FOUND = 1;
	    # for GENOPT(NONIa), turn off options for 1A model
	    if ( $mIa >= 0 ) { $USE_GENOPT[$mIa] = 0; }
	}
	else {
	    # keep only models where INFILE name contains keyArg string
	    for($m=0; $m<$NGENMODEL[$iver] ; $m++ ) {
		$MATCH=(index($SIMGEN_INFILE[$iver][$m],$keyArg)>=0);
		if ( $MATCH )  { $USE_GENOPT[$m] = 1 ; $FOUND=1; }
		else           { $USE_GENOPT[$m] = 0 ;     }
	    }
	}
    } # end of keyArg check
    
    # - - - - - append each GENOPT based on USE_GENOPT - - - - - 
    for($m=0; $m < $NGENMODEL[$iver]; $m++ ) {
	if ( $USE_GENOPT[$m] ) {
	    $GENOPT  = $GENVERSION_GENOPT[$iver][$m] ;
	    $GENVERSION_GENOPT[$iver][$m] = "$GENOPT $tmpArg" ;
	}
    }
    
    return ;

} # end PROCESS_GENOPT_LINE


# ==================================
sub read_GENOPT_STRING{

    # Created Oct 2019
    # check for a few optional items in GENOPT strings.

    my($iver) = @_ ;

    my ($m, $GENOPT, $val, @values, $INC);
    
    for($m=0; $m < $NGENMODEL[$iver]; $m++ ) {

	$GENOPT = "$GENVERSION_GENOPT[$iver][$m]" ;

	# parse for INCLUDE files (to copy)
	@values=sntools::parse_value_after_key("INPUT_FILE_INCLUDE",$GENOPT);
	# only add a new INCLUDE file to list
	foreach $val (@values) {
	    $INC  = "$GENOPT_FILE_INCLUDE[$iver]" ;
	    if(index($INC,$val)<0) {$GENOPT_FILE_INCLUDE[$iver]="$INC $val";}
	}
 	
    }
   

    return;

}   # end read_GENOPT_STRING

# ======================================
sub store_SIMGEN_INFILE(@) {

    my($iver,$iflag,$inFileList) = @_ ;

    # Feb 12 2019
    # Store simgen-input file(s) for this VERSION with index $iver.
    # iflag=0 -> load global SIMGEN files 
    # iflag=1 -> load version-dependent SNIa infile
    # iflag=2 -> load version-dependent NONIa infile
    # iflag>2 -> 2nd,3rd ... NON1a infile to override
    #
    # Mar 24 2020: if mIa<0, set mIa=0
    # Mar 31 2020: allow ENV for inFile.

    my ($m, $inFile, $inFile_Ia, $LDMP, $tmpFile );

    # ---------- BEGIN ----------
    
    $LDMP = 0 ;
    if ( $LDMP ) {
	print " xxx ------------------------------ \n";
	print " xxx iver=$iver  iflag=$iflag  NGENMODEL=$NGENMODEL[$iver]\n";
	if ( $iflag > 0 ) {
	    print " xxx inFileList = " ;
	    foreach $inFile (@$inFileList) { print "'$inFile' "; }
	    print "\n";
	}
    }

    # ----------------------------------------------------
    # make sure each inFile exists
    if ( $iflag > 0 ) {
	foreach $inFile (@$inFileList) { 
	    $inFile = qx(echo $inFile); 
	    $inFile =~ s/\s+$// ;   # trim trailing whitespace  
	    unless (-e $inFile ) {
		$MSGERR[0] = "Sim-input file '$inFile' does not exist;";
		$MSGERR[1] = "Check arguments of 'SIMGEN_INFILE_NONIa:'" ;
		$MSGERR[2] = "and 'SIMGEN_INFILE_Ia:' keys" ;
		sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	    }
	}
    }


    # ----------------------------------------------------
    if ( $iflag == 0 ) {
	$NGENMODEL[$iver] = $NGENMODEL_GLOBAL;
	for($m=0; $m < $NGENMODEL_GLOBAL; $m++ )  { 
	    $SIMGEN_INFILE[$iver][$m]  = $SIMGEN_INFILE_GLOBAL[$m]; 
	    $GENMODEL_CLASS[$iver][$m] =  $GENMODEL_CLASS_GLOBAL[$m] ;
	    $GENMODEL_NAME[$iver][$m]  =  $GENMODEL_NAME_GLOBAL[$m] ;
	}	
    }
    elsif ( $iflag == 1 ) {

	# make sure there is an SNIa model to override
	if ( $DOGEN_SNIa == 0 ) {
	    $MSGERR[0] = "Must define global SIMGEN_INFILE_Ia" ; 
	    $MSGERR[1] = "to override it under GENVERSION." ; 
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR); 
	}
	# overwrite or add SNIa model
	if ( $mIa < 0 ) {  # first SNIa input file
	    $mIa = 0; 
	    $NGENMODEL[$iver]++; 
	    $NGENMODEL_GLOBAL++ ; 
	}
	$SIMGEN_INFILE[$iver][$mIa]  =  @$inFileList[0] ; 
	$GENMODEL_CLASS[$iver][$mIa] =  "SNIa" ;
	$GENMODEL_NAME[$iver][$mIa]  =  $GENMODEL_NAME_GLOBAL[$mIa] ;
    }
    elsif ( $iflag == 2 ) {    
	# tricky: overwrite NON1A inFiles, but don't clobber SNIa
	$m=0;  
	if ( $mIa >= 0 ) { 
	    $inFile_Ia = $SIMGEN_INFILE[$iver][$mIa]; 
	    $SIMGEN_INFILE[$iver][$m]   = $inFile_Ia ; 
	    $GENMODEL_CLASS[$iver][$m]  = "SNIa" ;
	    $GENMODEL_NAME[$iver][$m]   = "$GENMODEL_NAME_GLOBAL[$m]" ;
	    $NGENMODEL[$iver]++ ;   
	    $m++ ;
	}
	foreach $inFile (@$inFileList)  { 
	    $SIMGEN_INFILE[$iver][$m]  = $inFile ; 
	    $GENMODEL_CLASS[$iver][$m] = "NONIaMODEL$m" ;
	    $GENMODEL_NAME[$iver][$m]  = "NotChecked" ;
	    $m++ ;     
	}
	$NGENMODEL[$iver] = $m ;
    }
    elsif ( $iflag > 2 ) {    
	# another SIMGEN_INFILE_NON1A key, so just add to list
	$m = $NGENMODEL[$iver];
	foreach $inFile (@$inFileList) {
	    $SIMGEN_INFILE[$iver][$m]  = $inFile ; 
	    $GENMODEL_CLASS[$iver][$m] = "NONIaMODEL$m" ;
	    $GENMODEL_NAME[$iver][$m]  = "NotChecked" ;
	    $m++ ;
	}
	$NGENMODEL[$iver] = $m;
    }

    return ;

} # end store_SIMGEN_INFILE 

# =========================================
sub parse_RANSEED {

    # Created July 2017
    # Parse the following:
    #  RANSEED  <SEED>
    #  RANSEED_REPEAT  <Nseed>  <Seed>   ! repeat same seed Nseed times
    #  RANSEED_CHANGE  <Nseed   <seed>   ! create Nseed different seeds
    #

    my ($key, @tmp, @wdlist,  $ijob);
    my ($RANSEED, $RANSEED_LAST, $RANSEED_FIRST );
    my $OPT_PARSE  = $OPT_QUIET ;

    $NRANJOB_SIMGEN = 1 ;
    $SAMEFLAG_RANSEED = 1; # default is all ranseeds the same

    $key = "RANSEED:";
    @tmp = sntools::parse_array($key,1,$OPT_PARSE,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) {
	$NRANJOB_SIMGEN = scalar(@tmp);
	for ( $ijob = 1; $ijob <= $NRANJOB_SIMGEN ; $ijob++ ) {
	    $RANSEED = $tmp[$ijob-1];
	    $SIMARG_RANSEED[$ijob] = "RANSEED $RANSEED" ;
#	    print " RANSEED[job $ijob]   = '$RANSEED'  \n" ;

	    # Turn off SAMEFLAG if ranseed changes.
	    # if SAMEFLAG_RANSEED =T, then all versions are combined;
	    # with different RANSEEDs, the versions stay separate.
	    if ( $ijob > 1 ) {
		$RANSEED_LAST = $tmp[$ijob-2] ;
		if ( $RANSEED != $RANSEED_LAST ) 
		{ $SAMEFLAG_RANSEED = 0; }
	    }
	}  # end ijob loop
    }  # end scalar if-block



    # ----------------------------------------------
    # check for command to releat the same RANSEED
    $key  = "RANSEED_REPEAT:" ;
    @tmp = sntools::parse_array($key,2,$OPT_PARSE,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) {
	@wdlist = split(/\s+/,$tmp[0]) ;
	$NRANJOB_SIMGEN = $wdlist[0] ;
	$RANSEED        = $wdlist[1] ;
	for ( $ijob = 1; $ijob <= $NRANJOB_SIMGEN ; $ijob++ ) 
	{ $SIMARG_RANSEED[$ijob] = "RANSEED $RANSEED" ; }
    }


    # ----------------------------------------------
    # check for command to releat the same RANSEED
    $key  = "RANSEED_CHANGE:" ;
    @tmp = sntools::parse_array($key,2,$OPT_PARSE,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0 ) {
	$SAMEFLAG_RANSEED = 0;        # each RANSEED is different
	@wdlist = split(/\s+/,$tmp[0]) ;
	$NRANJOB_SIMGEN = $wdlist[0] ;
	$RANSEED_FIRST  = $wdlist[1] ;  # first RANSEED given by user
	$DOSKIP_DUPLICATE_SIMJOBS = 0;  # Oct 29, 2019
	if ( $RESET_CIDOFF == 2 ) { $RESET_CIDOFF = 1; }
	for ( $ijob = 1; $ijob <= $NRANJOB_SIMGEN ; $ijob++ ) {
	    $RANSEED = $RANSEED_FIRST + 10000*$ijob + $ijob*$ijob + 13;;
	    $SIMARG_RANSEED[$ijob] = "RANSEED $RANSEED" ; 
	}
    }

    # ---------------------------------
    # print RANSEED summary

    for ( $ijob = 1; $ijob <= $NRANJOB_SIMGEN ; $ijob++ ) {
	print "\t  $SIMARG_RANSEED[$ijob]    ($ijob) \n" ;
    }

#    die "\n xxx DEBUG DIE xxxx \n";

    return ;


} # end parse_RANSEED

# =========================================
sub strip_GENOPT {

    my($FLAG,$GENOPT_IN) = @_ ;

    # FLAG=1 for Ia, FLAG=2 for NON1a

    my ($OPTFLAG_Ia, $OPTFLAG_NONIa, $GENOPT_OUT, @tmp);    
    $OPTFLAG_Ia = $OPTFLAG_NONIa = 0;
    $GENOPT_OUT = "$GENOPT_IN";

#    if ($string =~ m/foo/) { 

    # - - - - - - - 
    if ( $GENOPT_IN =~ m/NON1A/ ) { $OPTFLAG_NONIa = 1; }
    if ( $GENOPT_IN =~ m/NONIA/ ) { $OPTFLAG_NONIa = 1; }

    if ( $GENOPT_IN =~ m/SALT/  ) { $OPTFLAG_Ia = 1; }
    if ( $GENOPT_IN =~ m/GENMAG_SMEAR_MODELNAME/  ) { $OPTFLAG_Ia = 1; }

    # abort of both Ia and NON1A key are on same GENOPT line
    if ( $OPTFLAG_Ia && $OPTFLAG_NONIa ) {
	$MSGERR[0] = "Cannot mix Ia and NON1a keys on same GENOPT line.";
	$MSGERR[1] = "Check $SNMIX_INFILE_MASTER for ";
	$MSGERR[2] = "GENOPT:  $GENOPT_IN ";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    if ( $FLAG==1 && $OPTFLAG_NONIa ) { $GENOPT_OUT = ""; } # simIa with NON1A key
    if ( $FLAG==2 && $OPTFLAG_Ia    ) { $GENOPT_OUT = ""; } # simCC with Ia key

    return($GENOPT_OUT);

}  # end strip_GENOPT


# ========================================
sub parse_FORMAT_MASK {
    
    # Nov 22 2017: code moved from parse_inFile_master
    # Oct 24 2019: no longer check include file, and return if already set

    my ($key, @tmp) ;
    my $OPT_PARSE  = $OPT_QUIET ;

    $key = "FORMAT_MASK:" ;
    @tmp = sntools::parse_array($key,1,$OPT_PARSE,@CONTENTS_INFILE_MASTER );
    if ( scalar(@tmp) > 0  ) {
	
	if ( $FORMAT_MASK > 0 ) {
	    $MSGERR[0] = "FORMAT_MASK set twice ($FORMAT_MASK and $tmp[0])";
	    $MSGERR[1] = "FORMAT_MASK must be set once and only once";
	    $MSGERR[2] = "in simgen-master file.";
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);	    
	}
	else {
	    $FORMAT_MASK = $tmp[0] ;
	    $SIMARG_FORMAT = "FORMAT_MASK $FORMAT_MASK" ;
	}	
    }
    elsif ( $FORMAT_MASK > 0 ) {
	# do nothing if already read from GENOPT_GLOBAL 
    }
    else {
	# read format mask from sim-input file
	@tmp = sntools::parse_array($key,1,$OPT_ABORT,@CONTENTS_INFILE_SIMGEN);
	$FORMAT_MASK = $tmp[0] ;	
    }

} # end parse_FORMAT_MASK

# ========================================
sub check_FORMAT_MASK {

    # check that global $FORMAT_MASK is defined and if
    # it is text of FITS.

    my ($simCodeFile, $key, @tmp, $WRMASK_FITS, $WRMASK_CIDRAN, $cFORMAT );

    if ( $FORMAT_MASK <= 0 ) {
	$MSGERR[0] = "Invalid FORMAT_MASK = $FORMAT_MASK" ;
	$MSGERR[1] = "Must set 'FORMAT_MASK: <mask>' in $SNMIX_INFILE_MASTER";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }


    # read the sim-include file to check which bit is
    # for FITS format.
    $simCodeFile = "$SNANA_DIR/src/snlc_sim.h" ;

    $key = "WRMASK_FITS" ;
    @tmp = sntools::parse_line($simCodeFile, 1, $key, $OPT_ABORT) ;
    $WRMASK_FITS = $tmp[0] ;
    if ( $FORMAT_MASK & $WRMASK_FITS ) 
    { $WRFLAG_FITS = 1; $cFORMAT = "FITS"; }
    else
    { $WRFLAG_TEXT = 1; $cFORMAT = "TEXT";  }


    print " FORMAT_MASK = $FORMAT_MASK   ($cFORMAT)  \n";

    # ---------------------
    # check C-code include file if FORMAT_MASK includes RANDOM CID option.
    $key = "WRMASK_CIDRAN" ;
    @tmp = sntools::parse_line($simCodeFile, 1, $key, $OPT_ABORT) ;
    $WRMASK_CIDRAN = $tmp[0] ;
    if ( $FORMAT_MASK & $WRMASK_CIDRAN )  {
	$WRFLAG_CIDRAN = 1; 
	if ( $RESET_CIDOFF==0 ) { $RESET_CIDOFF = 1; }
    }

    # Dec 2 2106: abort if all random seeds are the same,
    #             but CIDRAN is not turned on
    if ( $SAMEFLAG_RANSEED && $WRFLAG_CIDRAN == 0 && $NRANJOB_SIMGEN>1 ) {
	$MSGERR[0] = "All RANSEED values are the same,";
	$MSGERR[1] = "but CIDRAN bit is off.";
	$MSGERR[2] = "Either add FORMAT_MASK += 16 ";
	$MSGERR[3] = "or set each RANSEED to different value.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;
    }	

} # end of check_FORMAT_MASK


# ===========================
sub set_timeStamp {

    # determine number of seconds since midnight to get unique
    # 'temp-versions' (avoids conflict if this script runs on
    # multuple nodes)

    my (@T_now, $Nsec);
    @T_now = localtime();
    $Nsec  = ($T_now[2] * 3600) + ($T_now[1] * 60) + $T_now[0] ;
    $Nsec5 = sprintf("%5.5d", $Nsec);
} 

# ============================
sub set_version_names {

    my($iver,$ijob) = @_;


    my ($VERSION, $verBase, $cjob, $ctmp, $version_temp, $version_final  );

    $VERSION = $GENVERSION_NAME[$iver];
    $cjob   = &get_jobidString($ijob );

    if ( $NRANJOB_SIMGEN == 1 ) {
	$verBase   = "${VERSION}" ;
    }
    else {
	$verBase   = "${VERSION}-${cjob}" ;
    }

    $version_final  = $verBase ;
    # multiple RANSEEDs
    if ( $NRANJOB_SIMGEN>1 && $SAMEFLAG_RANSEED ) 
    { $version_final = "${PREFIX_TEMP}_${verBase}" ; }
    
    $version_temp    = "${PREFIX_TEMP}_${verBase}_${SUFFIX_TEMP}" ;

    @VERSION_JOBLIST_FINAL[$ijob] = $version_final ;
    @VERSION_JOBLIST_TEMP[$ijob]  = $version_temp  ;

    $ctmp = "$version_temp -> $version_final";
    print "   VERSION[TEMP->FINAL] = $ctmp \n";

} # end of set_version_names 

# ====================================
sub get_jobidString {
    my($jobid) = @_ ;
    # Created July 2017
    # Return string for this jobid    
    my $cjobid = sprintf("%4.4d", $jobid);
    return($cjobid);
}

# ====================================
sub set_jobFileNames {

    my($iver,$jobid) = @_ ;
    my ($LDMP);

    # set global file names for this jobid

    $VERSION_TEMP     = "$VERSION_JOBLIST_TEMP[$jobid]" ;
    $SIMGEN_TEMPDIR   = "$PATH_SNDATA_SIM/${VERSION_TEMP}" ;

    $VERSION_FINAL    = "$VERSION_JOBLIST_FINAL[$jobid]" ;
    $SIMGEN_FINALDIR  = "$PATH_SNDATA_SIM/${VERSION_FINAL}" ;

    $SIMGEN_MISCDIR   = "${SIMGEN_FINALDIR}/$MISC_SDIR" ;

    $READMEFILE_FINAL = "${SIMGEN_FINALDIR}/${VERSION_FINAL}.README";
    $LISTFILE_FINAL   = "${SIMGEN_FINALDIR}/${VERSION_FINAL}.LIST";
    $IGNOREFILE_FINAL = "${SIMGEN_FINALDIR}/${VERSION_FINAL}.IGNORE";
    $DUMPFILE_FINAL   = "${SIMGEN_FINALDIR}/${VERSION_FINAL}.DUMP";

    $READMEFILE_TEMP = "${SIMGEN_TEMPDIR}/${VERSION_TEMP}.README";
    $DUMPFILE_TEMP   = "${SIMGEN_TEMPDIR}/${VERSION_TEMP}.DUMP";


    $LDMP = 0 ;
    if ( $LDMP ) {
	print "\n   xxxxxxxxxxxx JOBID = $jobid xxxxxxxxxxxxx \n";
	print "\t xxx SIMGEN_TEMPDIR   = $SIMGEN_TEMPDIR  \n" ;
	print "\t xxx SIMGEN_FINALDIR  = $SIMGEN_FINALDIR \n" ;
	print "\t xxx READMEFILE_FINAL = $READMEFILE_FINAL \n";
	print "\t xxx IGNOREFILE_FINAL = $IGNOREFILE_FINAL \n";
	print "\t xxx LISTFILE_FINAL   = $LISTFILE_FINAL \n";
	print "\t xxx DUMPFILE_FINAL   = $DUMPFILE_FINAL \n";
    }
} 


# ==============================
sub create_version {

    # Create $VERSION_FINAL where all SNe (Ia & NONIa)
    # will be moved to after generation.
    #
    # Apr 2012: copy SEARCHEFF_SPEC_FILE and replace system with qx.
    # Dec 2015: copy SIMLIB_FILE
    # Jan 2018: pass $iver argument to check for GENOPT_FILE_INCLUDE[$iver]
    # Jan 2019: remove <CR> from find to fix copy bug.

    my ($iver) = @_ ;

    my ($m, $cmdcp, $cmdclear, $cmd2, $cmd3, @FINC_LIST, @FINC_GENOPT, $finc);

    $cmdclear = "rm -r $SIMGEN_FINALDIR" ;
    $cmd2 = "mkdir $SIMGEN_FINALDIR" ;
    $cmd3 = "mkdir $SIMGEN_MISCDIR" ;

    if ( -d $SIMGEN_FINALDIR ) { qx($cmdclear); }
    qx($cmd2 ; $cmd3 );


    # copy input files to $SIMGEN_MISCDIR
    print "   Copy input  files to $SIMGEN_MISCDIR \n";

    $cmdcp = "cp $SNMIX_INFILE_MASTER $SIMGEN_MISCDIR/";
    qx($cmdcp);

    # create final DUMP file now so that new contents can
    # always be catenated to this existing file.
    if ( $DO_SIMGEN_DUMP ) {
	open  DUMPFILE , "> $DUMPFILE_FINAL" ;
	close DUMPFILE ;
    }

    for($m=0; $m < $NGENMODEL[$iver]; $m++ ) {
	$cmdcp = "cp $SIMGEN_INFILE[$iver][$m] $SIMGEN_MISCDIR/";
	qx($cmdcp);
    }


    # check all source of include file
    @FINC_GENOPT = split(/\s+/,$GENOPT_FILE_INCLUDE[$iver] ) ;
    @FINC_LIST = ( $INPUT_FILE_INCLUDE_Ia,  $INPUT_FILE_INCLUDE_NONIa,
		   @FINC_GENOPT );

    foreach $finc (@FINC_LIST ) {
	$finc =~ s/\n/ /;  # remove <CR>
	if ( $finc ne '' ) {
	    $cmdcp = "cp $finc $SIMGEN_MISCDIR/";
#	    print "\n\n 2. xxx copy-include-file command;\n $cmdcp \n\n";
	    qx($cmdcp);
	}
    }


    if ( $SEARCHEFF_SPEC_FILE ne '' ) {
	# copy file if it's there; otherwise will use public file
	if ( -e $SEARCHEFF_SPEC_FILE ) {
	    $cmdcp = "cp $SEARCHEFF_SPEC_FILE $SIMGEN_MISCDIR/";
	    qx($cmdcp);
	}
    }

    if ( $HOSTNOISE_FILE ne '' ) {
	$cmdcp = "cp $HOSTNOISE_FILE $SIMGEN_MISCDIR/";
	qx($cmdcp);
    }
    if ( $ZPHOTEFF_FILE ne '' ) {
	$cmdcp = "cp $ZPHOTEFF_FILE $SIMGEN_MISCDIR/";
	qx($cmdcp);
    }

    # copy optional ZVARIATION file
    if ( $INPUT_FILE_ZVAR_Ia ne '' ) {
	$cmdcp = "cp $INPUT_FILE_ZVAR_Ia $SIMGEN_MISCDIR/";
	qx($cmdcp);
    }
    if ( $INPUT_FILE_ZVAR_NON1A ne '' ) {
	$cmdcp = "cp $INPUT_FILE_ZVAR_NON1A $SIMGEN_MISCDIR/";
	qx($cmdcp);
    }

    # copy optional NON1AGRD (Mar 2016) if in current directory
    # (leave it  if it;s under $SNDATA_ROOT/models/NON1AGRID)
    if ( $INPUT_FILE_NON1AGRID[0] ne ''  && (-e $INPUT_FILE_NON1AGRID[0]) ) {
	$cmdcp = "cp $INPUT_FILE_NON1AGRID[0] $SIMGEN_MISCDIR/";
	qx($cmdcp);
    }

}  # end of create_version

# =======================================
sub init_duplicate {

    # Created Jan 2017
    # flag version-models with duplicate GENOPT.
    # Used later to create symbolic links to avoid 
    # re-generating a duplicate version-model.
    # This feature is most likely used by higher level scripts
    # such as split_and_fit & NEARNBR_pipeline.pl . 

    my ($iv0, $iv1, $m0, $m1 );

    print "\n\n Check for duplicate sim jobs: \n";

    # init array storing duplicate indices
    for ( $iv0=0 ; $iv0 <= $NGENVERSION; $iv0++ ) {
	$NMODEL_DUPL[$iv0] = 0;
	for($m0=0; $m0 < $NGENMODEL[$iv0]; $m0++ ) {
	    $IVERSION_DUPLICATE[$iv0][$m0]          = -9 ; 
	    $IMODEL_DUPLICATE[$iv0][$m0]            = -9 ; 
	    $IVERSION_DUPLICATE_INVERSE[$iv0][$m0]  = -9 ; 
	    $IMODEL_DUPLICATE_INVERSE[$iv0][$m0]    = -9 ; 
	}
    }

    if ( $DOSKIP_DUPLICATE_SIMJOBS == 0 ) { return; }

    # note that iv1 < iv0 --> iv0 is the duplicate
    for ( $iv0=2 ; $iv0 <= $NGENVERSION; $iv0++ ) {
	for($iv1=1; $iv1 < $iv0; $iv1++ ) {
	    for($m0=0; $m0 < $NGENMODEL[$iv0]; $m0++ ) { 
		for($m1=0; $m1 < $NGENMODEL[$iv1]; $m1++ ) { 
		    &match_duplicate($iv0,$iv1,$m0,$m1);		    
		} # end m0
	    } # end m0
	}     # end iv1
    }         # end iv0    

    if ( $NMODEL_DUPL_TOT > 0 && $RESET_CIDOFF != 2 ) {
	$MSGERR[0] = "To skip duplicate model jobs," ;
	$MSGERR[1] = "must set 'RESET_CIDOFF: 2' in $SNMIX_INFILE_MASTER";
	$MSGERR[2] = "so that every CID is unique among all GENVERSIONS.";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

#    die "\n\n xxx DIE init_duplicate_GENVERSION xxx \n";

} # end init_duplicate


# ============================================
sub match_duplicate {
    my ($iv0,$iv1,$m0,$m1) = @_ ;
    
    # Called from init_duplicate,
    # check if [version iv0, model m0] matches [version iv1, model m1].
    # "match" means same SIMGEN_INFILE and same GENOPT.
    # If there is a match, store info in DUPLICATE arrays.

    my ($INFILE0, $INFILE1, $GENOPTV0, $GENOPTV1, $IV0, $IV1, $string, $c2);

    $INFILE0 = $SIMGEN_INFILE[$iv0][$m0];
    $INFILE1 = $SIMGEN_INFILE[$iv1][$m1];

    # bail if sim-input files are different.
    if ( "$INFILE0" ne "$INFILE1" ) { return; }

    # next, check GENOPT
    $GENOPTV0 = "$GENVERSION_GENOPT[$iv0][$m0]" ;
    $GENOPTV1 = "$GENVERSION_GENOPT[$iv1][$m1]" ; 
    $IV0      = $IVERSION_DUPLICATE[$iv0][$m0] ;
    $IV1      = $IVERSION_DUPLICATE[$iv1][$m1] ;
    if ( ("$GENOPTV0" eq "$GENOPTV1") && $IV0 < 0 ) { 
	$NMODEL_DUPL_TOT++ ;
	$NMODEL_DUPL[$iv0]++ ;
	$c2 = sprintf("%2.2d", $NMODEL_DUPL_TOT);
	$IVERSION_DUPLICATE[$iv0][$m0]  = $iv1 ; 
	$IMODEL_DUPLICATE[$iv0][$m0]    = $m1 ; 
	$IVERSION_DUPLICATE_INVERSE[$iv1][$m1]  = $iv0 ; 
	$IMODEL_DUPLICATE_INVERSE[$iv1][$m1]    = $m0 ; 
	$string =
	    "DUPLICATE-$c2: " .
	    "$GENVERSION_NAME[$iv0]=$GENVERSION_NAME[$iv1] " .
	    "for $SIMGEN_INFILE[$iv0][$m0]";
	# store comment-string for README
	$STRING_COMMENT_DUPLICATE[$NMODEL_DUPL_TOT-1] = 
	    " $string";
	print "$string\n";
    }

    return ;

} # end match_duplicate {

# ============================================
sub final_SIMGEN_DUMP_duplicate {

    # Jan 2019
    # called after all sim jobs have finished;
    # make final duplicate update to SIMGEN_DUMP files.
    # Must wait for all jobs to finish to be sure that
    # all SIMGEN-DUMP files are ready to append.

    if ( $NMODEL_DUPL_TOT == 0 ) { return ; }
    if ( $DO_SIMGEN_DUMP  == 0 ) { return ; }

    print "\n Append DUMP files to account for skipped duplicates. \n";

    my($ivLink, $ivReal, $mLink, $mReal, $VERS_REAL, $VERS_LINK);
    my ($DUMPFILE_REAL, $DUMPFILE_LINK, $MM, $cat, $sedCmd, $gzFlag );

    for ( $ivLink=0 ; $ivLink <= $NGENVERSION; $ivLink++ ) {
	for($mLink=0; $mLink < $NGENMODEL[$ivLink]; $mLink++ ) {
	    $ivReal = $IVERSION_DUPLICATE[$ivLink][$mLink] ;
	    $mReal  = $IMODEL_DUPLICATE[$ivLink][$mLink] ;
	    if ( $ivReal < 0 ) { next ; }
	    if ( $mReal  < 0 ) { next ; }
	    $VERS_REAL = $GENVERSION_NAME[$ivReal] ;
	    $VERS_LINK = $GENVERSION_NAME[$ivLink] ;
	    $MM        = "MODEL$mReal" ;
	    $DUMPFILE_REAL  =  
		"$PATH_SNDATA_SIM/${VERS_REAL}/" . 
		"${VERS_REAL}.DUMP_${MM}" ;
	    $DUMPFILE_LINK  = 
		"$PATH_SNDATA_SIM/${VERS_LINK}/" . 
		"${VERS_LINK}.DUMP" ;
	    
	    $gzFlag = 0 ;	    	    
	    if ( -e "${DUMPFILE_LINK}.gz" ) 
	    { qx(gunzip ${DUMPFILE_LINK}.gz ); $gzFlag=1; }

# remove header keys and comments from real dump file
	    $sedCmd = 
		"sed -e '/NVAR/d' -e '/VARNAMES/d' -e '/#/d' " .
		" $DUMPFILE_REAL >> $DUMPFILE_LINK " ;
	    print " sed: $sedCmd \n";
	    system($sedCmd);
	    if($gzFlag) { qx(gzip $DUMPFILE_LINK); }
	}
    }

    # loop again and delete MM-dependent DUMP files
    for ( $ivLink=0 ; $ivLink <= $NGENVERSION; $ivLink++ ) {
	for($mLink=0; $mLink < $NGENMODEL[$ivLink]; $mLink++ ) { 
	    $ivReal = $IVERSION_DUPLICATE[$ivLink][$mLink] ;
	    $mReal  = $IMODEL_DUPLICATE[$ivLink][$mLink] ;
	    if ( $ivReal < 0 ) { next; }
	    if ( $mReal  < 0 ) { next; }
	    $VERS_REAL = $GENVERSION_NAME[$ivReal] ;
	    $MM        = "MODEL$mReal" ;
	    $DUMPFILE_REAL  = 
		"$PATH_SNDATA_SIM/${VERS_REAL}/" . 
		"${VERS_REAL}.DUMP_${MM}" ;
	    if ( $CLEANUP_FLAG ) { qx(rm $DUMPFILE_REAL); }
	}
    } 

    return ;

} # end final_SIMGEN_DUMP_duplicate

# =======================================
sub make_logDir {

    my ($cmd, $response);

    # always set name of log dir
    if ( $LOGDIR eq "" ) 
    { $LOGDIR       = "$TOPDIR_SIMLOGS/SIMLOGS_$GENPREFIX" ; }
    else {
	# if there are no slashes, then glue current pwd
	my $jslash  = rindex($LOGDIR,"/");  # location of last slash
	$LOGDIR = "$currentDir/$LOGDIR" ;
    }

    # if no done-stamp is specified, define a generic default stamp
    if ( $DONE_STAMP_FLAG == 0 ) {
	$DONE_STAMP      = "$LOGDIR/SIMJOB_ALL.DONE" ;
	$DONE_STAMP_FLAG = 1 ;
	print " set DONE_STAMP -> $DONE_STAMP \n";
    }

    # create log-directory for all batch command files and log files.

    if ( $DO_SSH || $DO_BATCH || $DEBUG_GENSTAT ) { return ; }

    # Aug 2017: if LOGDIR exists without DONE stamp, give SEVERE WARNING
    my $EXIST_LOGDIR = (-d $LOGDIR);
    my $EXIST_DONE   = (-e $DONE_STAMP );
    if ($PROMPT_FLAG &&  $EXIST_LOGDIR && $EXIST_DONE == 0 ) {
	print "\n\n ***** SEVERE WARNING ******** \n" ;
	print "LOGDIR = \n  '$LOGDIR' \n";
	print "exists without a done stamp. \n";
	print "Continuing could be very bad ... \n";
        print "Continue anyway y/[n] => " ;
        $response = <STDIN> ;
        $response =~ s/\s+$// ;
        unless ( $response eq "y" ) {
            print " Bye bye. \n";
            die " ***** ABORT ***** \n" ;
        }
    }


    if ( $EXIST_DONE   ) { qx(rm $DONE_STAMP) ; }
    if ( $EXIST_LOGDIR ) { qx(rm -r $LOGDIR)  ; }

    $cmd = "mkdir -p $LOGDIR ";
    print " Create log directory : $LOGDIR \n";
    qx($cmd) ;

    # copy master file to LOGDIR
    qx(cp $SNMIX_INFILE_MASTER $LOGDIR);

} # end of make_logDir


# ===================================
sub get_normalization {

    my($iver) = @_;

    # run SIM job with zero SNe to get calculated 
    # Number of SN per season. 
    # Note that the sim-job aborts, so grep the log-file
    # instead of the simulated README file.
    #
    # -------- BEGIN ----------
 
    my ($m, $NTMP);

    # make sure to get normalization only once per version
    if ( $STATUS_NORMALIZATION[$iver] >= 1 ) { return ; } 
    $STATUS_NORMALIZATION[$iver]++ ;

    # --------------------
    # if NGEN_UNIT is not given, just grep out NGENTOT_LC

    if ( $NGEN_UNIT_VAL < 0.0 ) {  
	for($m=0; $m < $NGENMODEL[$iver]; $m++ ) {
	    $NTMP = &get_NGEN($iver,$m); 
	    $GENMODEL_NGENTOT[$iver][$m] = $NTMP ;
	}
	goto DOSUM ;
    }


    for($m=0; $m < $NGENMODEL[$iver]; $m++ )  
    { &get_normalization_model($iver,$m); }

  DOSUM:

    &get_normalization_SUM($iver);

    return ;

} # end of &get_normalization


sub get_normalization_model {

    # Nov 12 2019: abort if normalization job fails
    # Apr 24 2020: abort if $NPER_SEASON is null string

    # iver = GENVERSION index, $m is model index
    my($iver,$m) = @_;

    my(@reqLine, $cmdNorm, $APPEND_NORM, @line, @wdlist, $NPER_SEASON );
    my($NGEN, $NGEN6, $NTMP);
    
    my $MODEL_CLASS  = "$GENMODEL_CLASS[$iver][$m]" ;
    my $MODEL_NAME   = "$GENMODEL_NAME[$iver][$m]" ; 
    my $GENVERSION   = "$GENVERSION_NAME[$iver]" ;
    my $VERSION      = "${GENVERSION}_norm" ;

    my $ARG_MODEL    = "GENMODEL $MODEL_NAME" ;
    my $ARG_VERS     = "GENVERSION $VERSION" ;
    my $ARG_NGEN     = "INIT_ONLY 1" ;

    my $SIMARG_GENOPT  = "$GENVERSION_GENOPT[$iver][$m] " ;
    my $SIMARG0   = "$SIMGEN_INFILE[$iver][$m]" ;
    my $SIMARG1   = "$ARG_VERS $ARG_NGEN " ;
    my $SIMARGS   = "$SIMARG0 $SIMARG1 " . 
	"$SIMARG_GENOPT " .   # xxx mark delete $GENOPT_GLOBAL " .
	"SIMLIB_MAXRANSTART 0" ;

    my $normLog  = "SIMnorm_${GENVERSION}_${MODEL_CLASS}.LOG" ;
    my $NORMLOG  = "$LOGDIR/$normLog" ;

    print "  Norm-STATUS($GENVERSION) = $STATUS_NORMALIZATION[$iver] \n";
    $| = 1;  # flush stdout

    # norm-log file may already exist if jobs are split among nodes.
    # Be careful that log file may exist but not be finished,
    # so sleep 2 extra seconds to make sure that log file is finished.
    if ( -e $NORMLOG ) {
	print "\t Norm-log file already exists:\n\t $normLog . \n";
	$NTMP = 0;
	# normLog exists, but make sure the QUIT line is there to
	# make sure it's really done
      CHECKDONE_NORM:	
	@reqLine = qx(grep QUIT $NORMLOG);
	if ( scalar(@reqLine) == 0 && $NTMP < 3)  
	{ sleep 1; $NTMP++; goto CHECKDONE_NORM;  }

	$cmdNorm = "sleep 1" ;
	$APPEND_NORM = 0 ;
    }
    else {     
	$cmdNorm  = "$JOBNAME_SIM $SIMARGS > $NORMLOG" ;
	$APPEND_NORM = 1 ;
    }
 
    print "   Get $MODEL_CLASS normalization ... " ;
    qx($cmdNorm) ;
    @line = qx { grep "per season ="  $NORMLOG }  ;
    if ( scalar(@line) == 0 ) {       
	$MSGERR[0] = "!!! SIMnorm job failed to get normalization !!! ";
	$MSGERR[1] = "See $normLog";
	$MSGERR[2] = "APPEND_NORM = $APPEND_NORM" ;
	$MSGERR[3] = "cmdNorm = '$cmdNorm' ";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;
    }

    @wdlist = split(/\s+/,$line[0]) ;
    $NPER_SEASON = $wdlist[7] ;
    if ( length($NPER_SEASON) == 0 ) {
	$MSGERR[0] = "Unable to extract NPER_SEASON from 7th word of";
	$MSGERR[1] = "grep line : $line[0]" ;
	$MSGERR[2] = "in NORMLOG file : $NORMLOG" ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;
    }
    print "NGENTOT_LC($MODEL_CLASS)/season = $NPER_SEASON \n";
    $| = 1;  # flush stdout
    if ( $APPEND_NORM ) { &append_normLog($NORMLOG); }
 

    # determine number to generate for SNIa model

    $NGEN = $NPER_SEASON * $NGEN_UNIT_VAL / $FASTFAC ;
    $NGEN = int($NGEN);
    $NGEN6 = sprintf("%6d", $NGEN);
    $GENMODEL_NGENTOT[$iver][$m] = $NGEN ;
    print "\t NGENTOT_LC = $NGEN6 for $MODEL_CLASS \n";
    $| = 1;  # flush stdout

    return ;

} # end of get_normalization_model

sub get_normalization_SUM {

    my($iver) = @_ ;

    my ($m, $NGENTOT_SUM, $NGENTOT, $NTMP, $CNN, $XTMP, $CIDRAN_ADD );

    $CNN='' ;
    for($m=0; $m < $NGENMODEL[$iver]; $m++ )  { 
	$NTMP     = $GENMODEL_NGENTOT[$iver][$m] ; 
	$NGENTOT += $NTMP;
	if ( $m == 0 ) 
	{ $CNN = "$NTMP"; }
	else
	{ $CNN = "$CNN,$NTMP"  } # increment error string.
    }

    $NGENTOT_SUM = int(($NGENTOT * $NRANJOB_SIMGEN)*1.01) + 10 ;
    if ( $NGENTOT_SUM < 100 ) { $NGENTOT_SUM = 100; }
    
    # Jan 10 2018: add 10% overhead to CIDRAN_ADD for more 
    #    efficient random-CID generation
    # BEWARE that random CIDs are very slow to generate 
    #  if CIDOFF is close to CIDRAN_MAX
    $CIDRAN_ADD  = $NGENTOT_SUM ;
    $CIDRAN_ADD *= 1.1 ;  $CIDRAN_ADD = int($CIDRAN_ADD);
    
    if ( $RESET_CIDOFF == 2 ) {
	if ( $iver==1 ) { $CIDRAN_MAX = $CIDRAN_MIN; } # jul 22 2018
	$CIDRAN_OFFSET[$iver] = $CIDRAN_OFF ;  # previous total
	$CIDRAN_MAX          += $CIDRAN_ADD ;
	$CIDRAN_OFF          += $NGENTOT_SUM ; # without 10% safety margin
    }
    else  { 
	$CIDRAN_MAX   = $CIDRAN_ADD ;
	$CIDRAN_OFFSET[$iver] = 0 ;
    }

	
    if ( $CIDRAN_MAX > $MXCID_SIM ) {
	$MSGERR[0] = "CIDRAN_MAX = $CIDRAN_MAX exceeds bound of";
	$MSGERR[1] = "MXCID_SIM = $MXCID_SIM (in snlc_sim.h). ";
	$MSGERR[2] = "NGENTOT_LC(Ia, NONIa) = $CNN ";
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR) ;
    }
    print "\n" ;

    return ;

} # end of get_normalization_SUM


# =========================
sub append_normLog { 

    # Jun 12 2013
    # to avoid confusion with ABORT message in the normLog,
    # add message that ABORT is normal. Masao has already
    # asked twice about these normal ABORTs.

    my ($normLogFile) = @_ ; 

    # open in append mode.
    open  PTR_NORMLOG , ">> $normLogFile" ;

    print PTR_NORMLOG "!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=! \n";
    print PTR_NORMLOG "\n";
    print PTR_NORMLOG "   Don't panic: abort was deliberate to read \n";
    print PTR_NORMLOG "   'Number of SN per season' above. \n";
    print PTR_NORMLOG "\n" ;
    print PTR_NORMLOG "!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=! \n";

    close PTR_NORMLOG

} # end of append_normLog

# ================
sub get_NGEN {

    my ($iver,$m) = @_ ;

    my $inFile = $SIMGEN_INFILE[$iver][$m];

    # May 2013
    # read NGEN_LC or NGENTOT_LC from sim-input files, and also
    # from GENOPT list.
    #
    # Jan 2019: refactor to use input model index $m

    my (@tmp1, @tmp2, $N1, $N2, $N, $key1, $key2, $key);

    $N1 = $N2 = $N = 0 ;

    $key1 = "NGEN_LC:";
    @tmp1 = sntools::parse_line($inFile, 1, $key1, $OPT_QUIET) ;
    if ( scalar(@tmp1) > 0 ) {
	if ( $WRFLAG_CIDRAN ) {
	    $MSGERR[0] = "NGEN_LC key is not compatible with random CIDs." ;
	    $MSGERR[1] = "Check sim-input file.";
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}
	$N1   = $tmp1[0] ;
    }

    $key2 = "NGENTOT_LC:";
    @tmp2 = sntools::parse_line($inFile, 1, $key2, $OPT_QUIET) ;
    if ( scalar(@tmp2) ) {  $N2   = $tmp2[0];  }

    if ( $N1 > 0 ) { $key = $key1; $N = $N1 ; }
    if ( $N2 > 0 ) { $key = $key2; $N = $N2 ; }

    if ( $N > 0 ) {
	print "   Will use $key $N  from $inFile\n";
    }
    else {
	$MSGERR[0] = "Could not find NGEN_LC or NGENTOT_LC in ";
	$MSGERR[1] = "$inFile" ;
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    
    # Jul 25 2017: check for NGENTOT_LC override in GENOPT list
    my ( $NVAL, @values, $GENOPT);
    $key2   = "NGENTOT_LC";
    $GENOPT = "$GENVERSION_GENOPT[$iver][$m]" ;
    @values = sntools::parse_value_after_key($key2,$GENOPT);
    $NVAL = scalar(@values) ;
    if ( $NVAL == 1 ) {
	$N = $values[0] ;
	print "   Will use $key2 $N  from GENOPT override.\n";
    }
    elsif ( $NVAL > 1 ) {
	$MSGERR[0] = 
	    "Found two NGENTOT_LC override values: $values[0] and $values[1]";
	$MSGERR[1] = "but only one allowed";
	$MSGERR[2] = "iver=$iver, $m=$m"; 
	sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
    }

    return $N ;

} # end of get_NGEN

# ================================
sub reset_counters {

    my ($m);

    $NSIMTOT_GEN  = 0 ;
    $NSIMTOT_WR   = 0 ;
    $NSIMTOT_SPEC = 0;

    for ( $m=0; $m < $NGENMODEL_MAX; $m++ ) {
	$NSIM_GEN[$m]   = 0 ;
	$NSIM_WR[$m]    = 0 ;
	$NSIM_SPEC[$m]  = 0 ;
    }

} # end  of reset_counters

# ===========================================
sub simgen {

    my($iver,$jobid,$m) = @_;

    # Construct arguments for snlc_sim.exe and  run the job
    # for this JOBID and  model (m). 
    # Generate SN for $iver(version), $jobid(RANSEED) and m=model 
    #   m=0 for SNIa
    #   m>0 for NONIa
    # Generate in the SIMGEN_TEMP area, and then copy
    # the files to the SIMGEN_FINAL area.
    #
    # -------- BEGIN --------
    
    my ($TMP_MODEL, $TMP_NGEN, $TMP_FILE, $TMP_TYPE, $TMP_GENOPT);
    my ($TMP_SIMARG_MODEL, $TMP_CIDOFF, $TMP_CIDRAN_MAX, $TMP_CLASS );
    my ($iarg, $c_NGEN, $c_CIDOFF, $c_JOBID, $c_MODEL);
    my ($logFile, $simcmd, $key, @tmp, @wdlist, $NWR, $NGEN, $NSPEC);
    my ($GENPREFIX_LOCAL, $cmdFile, $doneFile );
    my $qq = '"' ;

    $TMP_MODEL    = "$GENMODEL_NAME[$iver][$m]" ;
    $TMP_CLASS    = "$GENMODEL_CLASS[$iver][$m]" ;
    $TMP_NGEN     = "$GENMODEL_NGENTOT[$iver][$m]" ;
    $TMP_FILE     = "$SIMGEN_INFILE[$iver][$m]" ;
    $TMP_GENOPT   = "$GENVERSION_GENOPT[$iver][$m]" ;
    if ( $TMP_MODEL eq "FIXMAG" ) { $TMP_CLASS = "FIXMAG" ; }

    # set filenname prefix based on output type
    if ($WRFLAG_TEXT ) { 
	# Ia and CC text files have same prefix -> cannot tell them apart
	$GENPREFIX_LOCAL = $GENPREFIX ; 
    } 
    else { 
	# Ia and CC FITS file-names must be different
	$GENPREFIX_LOCAL = "${GENPREFIX}_${TMP_CLASS}" ; 
    }  

    # 4/30/2014 - for multiple RANSEEDs that are the same, 
    #             tack on RANSEED index to GENPREFIX
    if ( $SAMEFLAG_RANSEED && $NRANJOB_SIMGEN > 1 ) {
	my $cjob = &get_jobidString($jobid);
	$GENPREFIX_LOCAL = "${GENPREFIX_LOCAL}-${cjob}" ;
    }

    $NSIMARG        = 0;
    $CIDOFF        += $LAST_NGEN ;

    $NSIMARG++ ;
    $SIMARG_LIST[$NSIMARG]  = "$TMP_FILE";

    if ( $USER_PATH_SNDATA_SIM ) {
	$NSIMARG++ ;
	$SIMARG_LIST[$NSIMARG]  = "PATH_SNDATA_SIM $PATH_SNDATA_SIM" ;
    }

    # May 5 2017: new logic for adding CIDRAN_MAX so that it works
    #    with both NGEN_UNIT and NGENTOT_LC
    if ( $SAMEFLAG_RANSEED && $NRANJOB_SIMGEN>1 ) {
        $NSIMARG++ ;
        $SIMARG_LIST[$NSIMARG]  = 
	    "CIDRAN_MIN $CIDRAN_MIN  CIDRAN_MAX $CIDRAN_MAX" ;
    }

    $NSIMARG++ ;  $TMP_CIDOFF = $CIDOFF + $CIDRAN_OFFSET[$iver];
    $SIMARG_LIST[$NSIMARG]  = "CIDOFF $TMP_CIDOFF" ;

    $NSIMARG++ ;
    $SIMARG_LIST[$NSIMARG]  = "JOBID $jobid NJOBTOT $NRANJOB_SIMGEN";

    $NSIMARG++ ;
    $SIMARG_LIST[$NSIMARG]  = "$TMP_SIMARG_MODEL";

    $NSIMARG++ ;
    $SIMARG_LIST[$NSIMARG]  = "GENVERSION  ${VERSION_TEMP}" ;
    $NSIMARG++ ;
    $SIMARG_LIST[$NSIMARG]  = "GENPREFIX  ${GENPREFIX_LOCAL}" ;

    $NSIMARG++ ;
    $SIMARG_LIST[$NSIMARG]  = "$SIMARG_RANSEED[$jobid]" ;

    $NSIMARG++ ;
    $SIMARG_LIST[$NSIMARG]  = "$SIMARG_FORMAT" ;

    $NSIMARG++ ;
    $SIMARG_LIST[$NSIMARG]  = "$TMP_GENOPT" ;

    if ( $NGENMODEL[$iver] > 1 ) {
	# suppress model-dependent params from SIMGEN-DUMP and data files
	$NSIMARG++ ;
	$SIMARG_LIST[$NSIMARG]  = "WRFLAG_MODELPAR 0" ;
    }

    if ( $NGEN_UNIT_VAL > 0 ) {
	$NSIMARG++ ;
	$SIMARG_LIST[$NSIMARG]  = "NGENTOT_LC $TMP_NGEN  NGEN_LC 0";
    }

    $SIMARG_ALL = "" ;
    for ( $iarg=1; $iarg <= $NSIMARG; $iarg++ ) {
	$SIMARG_ALL = "$SIMARG_ALL " . " $SIMARG_LIST[$iarg]" ;
    }

    $c_NGEN   = sprintf("NGEN=%5d", $TMP_NGEN);
    $c_CIDOFF = sprintf("CIDOFF(MAX)=%5d(%5d)", $CIDOFF, $CIDRAN_MAX  );
    $c_JOBID  = &get_jobidString($jobid);

    $c_MODEL  = "MODEL=$TMP_MODEL";
    print "    Sim-$c_JOBID  $c_NGEN  $c_CIDOFF  $c_MODEL \n";

    $logFile  = "${LOGDIR}/${VERSION_TEMP}_${TMP_CLASS}.LOG" ;
    $cmdFile  = "${LOGDIR}/${VERSION_TEMP}_${TMP_CLASS}.README" ;
    $simcmd   = "$JOBNAME_SIM $SIMARG_ALL > $logFile " ;
    
    # write ENV variables and list of commands to README file 
    open  PTR_CMDFILE , "> $cmdFile" ;
    print PTR_CMDFILE " SNANA_DIR       = $SNANA_DIR \n";
    print PTR_CMDFILE " SNANA_MODELPATH = $SNANA_MODELPATH \n\n";
    print PTR_CMDFILE " Executed command: \n\t $simcmd \n";
    close PTR_CMDFILE ;

    # run job and wait for it to finish ...
    system("$simcmd");
    sleep(1); 

    # May 2020: make sure log file exists
    if ( !(-e $logFile) ) {
	print "\n xxx cannot find logFile = $logFile ??? xxx \n";
	$JOBABORT_FLAG = 1 ;   $NABORT++ ;
        goto LAST_NGEN ;
    }

    # Jun 7 2013: if sim aborted, exit here to make sure DONE stamp is made.
    $key = " ABORT " ;
    @tmp = sntools::parse_line($logFile, 0, $key, $OPT_QUIET) ;
    if ( scalar(@tmp) > 0 ) { 
	$JOBABORT_FLAG = 1 ;	$NABORT++ ;
#	print "\n xxx found ABORT in $logFile \n\n";
	goto LAST_NGEN ; 
    }

    # no abort if we get here.
    $JOBABORT_FLAG = 0 ; 

    # increment number of SN generated and written out
    $key = "Wrote" ;
    @tmp = sntools::parse_line($READMEFILE_TEMP, 1, $key, $OPT_ABORT) ;
    @wdlist   = split(/\s+/,$tmp[0]) ;
    $NWR      = $wdlist[0] ;

    $key = "Generated" ;   # Apr 24 2013
    @tmp = sntools::parse_line($READMEFILE_TEMP, 1, $key, $OPT_ABORT) ;
    @wdlist   = split(/\s+/,$tmp[0]) ;
    $NGEN     = $wdlist[0] ;

    $key = "Spectroscopic-type:" ;
    @tmp = sntools::parse_line($READMEFILE_TEMP, 3, $key, $OPT_ABORT) ;
    @wdlist   = split(/\s+/,$tmp[0]) ;
    $NSPEC    = $wdlist[2] ;

    $NSIMTOT_GEN   += $NGEN ; 
    $NSIMTOT_WR    += $NWR ;
    $NSIMTOT_SPEC  += $NSPEC ;
    $NSIM_GEN[$m]  += $TMP_NGEN ;
    $NSIM_WR[$m]   += $NWR ;
    $NSIM_SPEC[$m] += $NSPEC ;


    # - - - - - - - - - - - - - - - - - - - - - -
    # Apr 2019: store extra sums as follows:
    # SUM-SNIa:   $m = $Nm
    # SUM-NONIa:  $m = $Nm+1
    # SUM-ALL:    $m = $Nm+2  (same as NSIMTOT_XXX)
    my $Nm = $NGENMODEL[$iver];
    my ($mm);
    if($m == $mIa) {$mm = $Nm;}   else{$mm=$Nm+1;}
    $NSIM_GEN[$mm]  += $NGEN ; 
    $NSIM_WR[$mm]   += $NWR  ;  
    $NSIM_SPEC[$mm] += $NSPEC; 
    
    $NSIM_GEN[$Nm+2]  = $NSIMTOT_GEN; 
    $NSIM_WR[$Nm+2]   = $NSIMTOT_WR; 
    $NSIM_SPEC[$Nm+2] = $NSIMTOT_SPEC; 

  LAST_NGEN:
    $LAST_NGEN = $TMP_NGEN ; 

    return ;

} # end of simgen


# =================================
sub simcopy {

    # copy/move sim files (README and data) into final version-directory

    my($iver,$jobid,$m) = @_;

    my ($TMP_MODEL, $TMP_CLASS, $VERSION, $readme_archive);
    my ($cdd, $cmd, $cmdmv, $tempDump, $cmdcat, $dataFiles, $PREFIX ) ;
    my ($PREFIX_FITS_NEW, $PREFIX_FITS_OLD );

    $VERSION    = $VERSION_FINAL ;
    $TMP_MODEL  = "$GENMODEL_NAME[$iver][$m]" ; 
    $TMP_CLASS  = "$GENMODEL_CLASS[$iver][$m]" ;
    
    # save original README file
    $PREFIX = "${GENPREFIX}_${TMP_CLASS}";
    if ($NRANJOB_SIMGEN > 1) { 
	my $cjob = &get_jobidString($jobid); 
	$PREFIX = "$PREFIX-$cjob" ;
    }

    $readme_archive = "${SIMGEN_MISCDIR}/${PREFIX}.README" ;
    $cmd = "cp $READMEFILE_TEMP $readme_archive";
    system("$cmd");

    $cdd = "cd $SIMGEN_TEMPDIR";

    # move data files into FINAL version-directory;
    # note distinction between TEXT and FITS format.
    if ( $WRFLAG_TEXT ) {
	$dataFiles = "${GENPREFIX}*.DAT" ;
    }
    else {
	$dataFiles = "${GENPREFIX}*.FITS" ;
    }

    $cmdmv     = "mv ${dataFiles}  ${SIMGEN_FINALDIR}/" ;

    system("$cdd ; $cmdmv");

    if ( $DO_SIMGEN_DUMP == 0 ) { return ; }

    # catenate contents of DUMP file

    $tempDump = "${VERSION}_DUMP.TEMP" ;
    $cmdcat   = "cat $DUMPFILE_TEMP $DUMPFILE_FINAL > $tempDump" ;
    $cmdmv    = "mv $tempDump $DUMPFILE_FINAL" ;	
    system("$cmdcat ; $cmdmv ");

    # if this version+model is going to be duplicated later,
    # need to preserve a separate DUMP-FILE with MODEL$m subset
    # to add later into future DUMP versions.
    if ( $IVERSION_DUPLICATE_INVERSE[$iver][$m] >= 0 ) {
	my ($dumpFile_model_out, $dumpFile1, $dumpFile2);
	$dumpFile_model_out = "${DUMPFILE_FINAL}_MODEL${m}" ;
	$dumpFile1 = $DUMPFILE_TEMP;
	$dumpFile2 = "" ;
	if ( -e $dumpFile_model_out ) { $dumpFile2 = $dumpFile_model_out; }

	$cmdcat   = "cat $dumpFile1 $dumpFile2 > $tempDump" ;
	$cmdmv    = "mv $tempDump $dumpFile_model_out" ;
	system("$cmdcat ; $cmdmv ");	
    }

    
} # end of simcopy

# ============================
sub printSummary {

    my ($IVER) = @_  ;

    # Created Jan 5 2017
    # print stdout summary for current job.
    # [moved from end of make_AUXFILE_LIST, and added IVER_DUPL logic]

    my $IVER_STRING  = &get_jobidString($IVER);
    my $txt0 = " ${IVER_STRING}  SUMMARY for ${VERSION_FINAL}" ;

    my ($wtxt, $wtxt2);  $wtxt = $wtxt2 = "" ;
    $wtxt = "wrote $NSIMTOT_WR ($NSIMTOT_SPEC) of $NSIMTOT_GEN " ; 
    if ( $NSIMTOT_GEN == 0 ) { $wtxt2 = "(probably aborted)"; }
    print "${txt0} -> $wtxt  $wtxt2 \n" ;

    return ;

}  # end printSummary

# ===============================================
sub make_AUXFILE_LIST {

    # Sep 30 2019:
    #  list *SNIa* before *NONIaMODEL* ... matters for RANSEED_CHANGE.

    my($iver,$jobid) = @_;

    my ($NLISTFILE, $i, @tmp, $prefix, $cmd, $L, $headFiles);

    print "   Create AUXILIARY LIST   file for $VERSION_FINAL \n";

    open  PTR_LIST , "> $LISTFILE_FINAL" ;
    close PTR_LIST ;

    # create and fill LIST file.
    # ls in 10 steps (0-9) in case there are too many files
    # to ls them all in one step.

    if ( $WRFLAG_TEXT ) {
	for ( $i=0; $i <= 9 ; $i++ ) {	
	    $prefix = "${GENPREFIX}_SN${i}" ;
	    $cmd    = "ls ${prefix}* >> $LISTFILE_FINAL";
	    system("cd $SIMGEN_FINALDIR ; $cmd 2>/dev/null");
	}

	# count entries in LIST file to make sure that it matches
	# local counter (text only)

	@tmp = `cat $LISTFILE_FINAL` ;
	$NLISTFILE = scalar(@tmp);

	print "\t Finished Generating $NSIMTOT_GEN SNe. \n";
	$L = "${VERSION_FINAL}.LIST";
	print "\t Found $NLISTFILE entries in ${L}; expected $NSIMTOT_WR \n";
	if ( $NLISTFILE != $NSIMTOT_WR ) {
	    print "\t ==> ERROR in sim-file count !!!! \n";
	}

    }
    else {
	# FITS format SNIa, then NONIaMODEL 

	if ( $DOGEN_SNIa ) {
	    $headFiles = "${GENPREFIX}*SNIa*HEAD.FITS" ;
	    $cmd = "ls ${headFiles} >> $LISTFILE_FINAL";
	    system("cd $SIMGEN_FINALDIR ; $cmd 2>/dev/null");
	}

	if ( $DOGEN_NONIa ) {
	    $headFiles = "${GENPREFIX}*NONIaMODEL*HEAD.FITS" ;
	    $cmd = "ls ${headFiles} >> $LISTFILE_FINAL";
	    system("cd $SIMGEN_FINALDIR ; $cmd 2>/dev/null");
	}

	print "\t Finished Generating $NSIMTOT_GEN SNe. \n";
	print "\t Wrote $NSIMTOT_WR passing trigger and cuts.\n";
    }

    return ;

}  #  end of make_AUXFILE_LIST

# ===============================================
sub make_AUXFILE_IGNORE {

    my($iver,$jobid) = @_;
    print "   Create AUXILIARY IGNORE file for $VERSION_FINAL \n";
    open  IGNOREFILE , "> $IGNOREFILE_FINAL" ;
    close IGNOREFILE ;


} #  end of make_AUXFILE_IGNORE

# ===============================================
sub make_AUXFILE_README {

    my($iver,$jobid) = @_;

    my $Nm = $NGENMODEL[$iver];
    my ($tnow, $m, $NAME, $CLASS, $NGEN, $NWR, $NSPEC);
    my (@wdlist, $itag, $tag, $ctag, $ctmp, $GENOPT);

    print "   Create AUXILIARY README file for $VERSION_FINAL \n";

    open  READMEFILE , "> $READMEFILE_FINAL" ;

    $tnow    = localtime time ;

    print READMEFILE " Created $tnow \n";
    print READMEFILE " by script: $0 \n" ;
    print READMEFILE " FLUXCAL = 10^[-0.4*mag + 11] \n" ;

    print READMEFILE " \n";
    print READMEFILE " This directory has $NSIMTOT_WR total SN candidates \n";
    print READMEFILE " corresponding to $NGEN_UNIT_VAL $NGEN_UNIT_NAME \n" ;
    print READMEFILE " \n";

    print READMEFILE " Generation Breakdown: \n";
    print READMEFILE
	"#   model  CLASS          N(GEN)    N(CUTS)   N(spec) \n" ;
    print READMEFILE 
	"# ------------------------------------------------------ \n" ;


    for ( $m=0; $m < $Nm+3; $m++ ) {  
	$CLASS = $GENMODEL_CLASS[$iver][$m] ;
	$NGEN  = $NSIM_GEN[$m] ;
        $NWR   = $NSIM_WR[$m] ;
	$NSPEC = $NSIM_SPEC[$m] ;
	
	if ( $m == $Nm   ) { $CLASS = "SUM-SNIa";  }
	if ( $m == $Nm+1 ) { $CLASS = "SUM-NONIa"; }
	if ( $m == $Nm+2 ) { $CLASS = "SUM-ALL";   }

	$ctmp  = sprintf("     %2d   %-12s  %7d     %6d    %6d", 
			 $m, $CLASS,$NGEN,$NWR,$NSPEC);
	print READMEFILE " $ctmp \n";
    }

    &printLegend_INFILE($iver,\*READMEFILE);
    &printLegend_GENOPT($iver,\*READMEFILE);

    close READMEFILE ;

} #  end of make_AUXFILE_README


# ============================
sub write_GENSTAT_DRIVER {
    my ($iver,$jobid) = @_ ;

    # write job stats to table that can be analyzed later
    # This file contributes later to README and TOTAL_SUMMARY.LOG
    #
    # --> write_GENSTAT_open  ! one-time info for each GENVERSION
    # --> write_GENSTAT_job   ! update for each job
    # --> write_GENSTAT_sum   ! sum jobids & links
    #          -> write_GENSTAT_table  ! update sum-tables with sed
    #          -> move_GENSTAT_FILE
    
    my ($iv, $GV, $README_FILE, $GENVERSION, $STATFILE, $BUSYFILE);

    if ( $iver<0 && $jobid < 0 ) {
	# called at end to sum everything

	# sum everything over jobids & links
	&write_GENSTAT_sum();  

	if ( $SAMEFLAG_RANSEED ) {
	    print " Move GENSTAT files ... \n";
	    for($iv=1; $iv<=$NGENVERSION; $iv++ ) { &move_GENSTAT_FILE($iv); }
	}
	return ;
    }


    # - - - - - - - -
    $GENVERSION = "$GENVERSION_NAME[$iver]" ;
    $STATFILE   = &get_FILENAME_GENSTAT($iver);
    $BUSYFILE   = "$LOGDIR/BUSY_GENSTAT_${GENVERSION}" ;

    # avoid different batch jobs trying to open STATFILE
    # at the same time. Use BUSYFILE to control traffic.
    # Increase max BUSY wait from 10s to 60s (for NERSC).
    my $NSLEEP = 0 ;
    while ( -e $BUSYFILE ) { 
	sleep(1);  
	$NSLEEP++ ;
	if ($NSLEEP > 60 ) {
	    $MSGERR[0] = "BUSYFILE still there after $NSLEEP seconds";
	    $MSGERR[1] = "BUSYFILE = $BUSYFILE" ;
	    sntools::FATAL_ERROR_STAMP($DONE_STAMP,@MSGERR);
	}
    }



    if ( -e $STATFILE ) {
	# open in append mode
	open  PTR_TABLE , ">> $STATFILE" ;
    }
    else {
	# open new file
	open  PTR_TABLE , "> $STATFILE" ;
	&write_GENSTAT_open($iver,\*PTR_TABLE);
    }


    # update table for each job.
    open  PTR_BUSY  , "> $BUSYFILE" ;
    &write_GENSTAT_job($iver, $jobid, \*PTR_TABLE);
    close PTR_TABLE ;
    close PTR_BUSY ;   qx(rm $BUSYFILE);

    return ;

} # end write_GENSTAT_DRIVER


# ===============================
sub write_GENSTAT_open {
    my($iver,$FH) = @_ ;

    # write one-time info to newly opened GENSTAT file
    my ($m, $OPT, $inFile );
    my $GENVERSION = "$GENVERSION_NAME[$iver]" ;

    if ( $SAMEFLAG_RANSEED  ) {
	# all random versios combined into one
	print $FH "GENVERSION: $GENVERSION \n\n";
    }
    else {
	# each randomn version remains with jobid tag
	my $maxjobid  = &get_jobidString($NRANJOB_SIMGEN);
	print $FH "GENVERSION: ${GENVERSION}-0001 \n";
	print $FH "... \n";
	print $FH "GENVERSION: ${GENVERSION}-${maxjobid} \n";
    }

    # create table mapping legend between model (m) and INFILE
    &printLegend_INFILE($iver,\*$FH);
    &printLegend_GENOPT($iver,\*$FH);


    # leave space to insert summed table later
    print $FH "\n";
    $OPT = $OPT_GENSTAT_TABLE_plusLINKS ;
    print $FH "$TABLE_INSERT_MARK${OPT}:\n\n"; 
    
    print $FH "\n   *** TABLES BELOW ARE FOR DEBUG *** \n\n";

    if ( $NMODEL_DUPL_TOT ) {
	$OPT = $OPT_GENSTAT_TABLE_noLINKS ;
	print $FH "$TABLE_INSERT_MARK${OPT}:\n\n"; 
    }

    # - - - - - - - - - 
    # prepare table header for each job,versio,model
    print $FH 
	"$GENVERSION stats for each individiaul sim-job\n";
    print $FH
	"#   jobid iver  m   iver_link  m_link     " .
	"NGEN    NWRITE    NSPEC\n";
    print $FH
	"# ------------------------------------------" . 
	"----------------------\n";	

    $FH -> autoflush(1);

}  # end write_GENSTAT_open


# ===================
sub write_GENSTAT_job {

    # for input iver and jobid, write stats for each sub-model.

    my($iver,$jobid,$FH) = @_ ;

    my($m, $iver_link, $m_link, $CLASS, $NGEN, $NWR, $NSPEC, $LINE);

    for ( $m=0; $m < $NGENMODEL[$iver]; $m++ ) {  
	$iver_link = $IVERSION_DUPLICATE[$iver][$m] ;
	$m_link    = $IMODEL_DUPLICATE[$iver][$m]  ;
	$CLASS     = $GENMODEL_CLASS[$iver][$m] ;
	$NGEN      = $NSIM_GEN[$m] ;
        $NWR       = $NSIM_WR[$m] ;
	$NSPEC     = $NSIM_SPEC[$m] ;
	$LINE = sprintf("%2d   %3d  %2d     %3d    %5d   %9d %8d %8d ",
			$jobid, $iver, $m, $iver_link, $m_link,  
			$NGEN, $NWR, $NSPEC);
	print $FH "JOB: $LINE\n";
    }

    $FH -> autoflush(1);

} # end write_GENSTAT_job


# =============================================
sub write_GENSTAT_sum {

    # Apr 1 2019
    # Called after all jobs have completed, sum stats over jobids
    # and update GENSTAT* file for each version. Note that ALL
    # GENVERSIONs must be finished to account for links that avoid
    # duplicate generation.
    #
    # Each GENSTAT file is updated with two tables:
    #   * all stats including links (for analysis use)
    #   * gen-only stats without links (for CPU monitor)
    #
    # Writing is via 'sed', not print.
    #
    # For combining random sim-jobs, move each GENSTAT file to REAMDE
    #
    # JSTAT = 0,1,2 for NGEN,NWRITE,NSPEC

    my($iver, $GENSTAT_FILE, $GENVERSION, $KEY, $NKEY, $Nm);
    my (@tmp, @wdlist, $tmpLine, $j, $OPT, $README_FILE );
    my ($jobid, $iv, $m, $ivLink, $mLink, $mSum, $NGEN, $NWRITE, $NSPEC);
    my $JSTAT_NGEN   = 0;
    my $JSTAT_NWRITE = 1;
    my $JSTAT_NSPEC  = 2;

    print "\n";
    print "# -------------------------------------------------- \n";
    print "  Analyze & Update GENSTAT Tables \n";

    # - - - - -
    # first loop over each versoin and sum over jobid; note that
    # link have zero events.

    $KEY = "JOB:" ;
    for($iver=1; $iver <= $NGENVERSION; $iver++ ) {
	$GENVERSION     = $GENVERSION_NAME[$iver];
        $GENSTAT_FILE   = &get_FILENAME_GENSTAT($iver);
	@tmp     = sntools::parse_line($GENSTAT_FILE, 8, $KEY, $OPT_QUIET) ;
	$NKEY    = scalar(@tmp);
	foreach $tmpLine (@tmp) {
	    @wdlist  = split(/\s+/,$tmpLine ) ; 
	    $jobid      = $wdlist[0];
	    $iv         = $wdlist[1];
	    $m          = $wdlist[2];
	    $ivLink     = $wdlist[3];
	    $mLink      = $wdlist[4];
	    $NGEN       = $wdlist[5];
	    $NWRITE     = $wdlist[6];
	    $NSPEC      = $wdlist[7];
	    $GENSTAT_SUMJOB[$iv][$m][$JSTAT_NGEN]   += $NGEN ;  # sum over jobid
	    $GENSTAT_SUMJOB[$iv][$m][$JSTAT_NWRITE] += $NWRITE ; 
	    $GENSTAT_SUMJOB[$iv][$m][$JSTAT_NSPEC]  += $NSPEC ; 
	    $GENSTAT_ivLink[$iv][$m]         = $ivLink ;
	    $GENSTAT_mLink[$iv][$m]          = $mLink ;	   
	}	
    } # end iver loop over GENVERSIONs

    
# ------------------------------------------
# loop again and add in sums from links to get SUMTOT

    for($iv=1; $iv <= $NGENVERSION; $iv++ ) {
	$Nm = $NGENMODEL[$iv]; 
	for ( $m=0; $m < $Nm; $m++ ) {
	    for ( $j=0; $j <=2; $j++ ) {
		$GENSTAT_SUMTOT[$iv][$m][$j]   = $GENSTAT_SUMJOB[$iv][$m][$j];
		
		$ivLink = $GENSTAT_ivLink[$iv][$m] ;
		$mLink  = $GENSTAT_mLink[$iv][$m] ;
		if ( $ivLink >= 0 && $mLink >= 0 ) {
		    $GENSTAT_SUMTOT[$iv][$m][$j] += 
			$GENSTAT_SUMJOB[$ivLink][$mLink][$j];
		}
		
		# keep track of Ia and NONIa sums
		if ( $m == $mIa ) {$mSum=$Nm;} else {$mSum=$Nm+1;}
		$GENSTAT_SUMJOB[$iv][$mSum][$j]+=$GENSTAT_SUMJOB[$iv][$m][$j];
		$GENSTAT_SUMTOT[$iv][$mSum][$j]+=$GENSTAT_SUMTOT[$iv][$m][$j];
		
	    } # end j loop over NGEN, NWRITE, NSPEC
	}
    } # end iver loop

# - - - - - -
# Finally, update stats with & without links

    for($iv=1; $iv<=$NGENVERSION; $iv++ ) { 

	$OPT = $OPT_GENSTAT_TABLE_plusLINKS;
	&write_GENSTAT_table($iv,$OPT); 

	if ( $NMODEL_DUPL_TOT ) {
	    $OPT = $OPT_GENSTAT_TABLE_noLINKS;
	    &write_GENSTAT_table($iv,$OPT); 
	}
    }
    

    return ;

} # end write_GENSTAT_sum; 

# =======================
sub move_GENSTAT_FILE {
    # move GENSTAT file to README for combine-ranjob
    my($iver) = @_ ;
    my $GV = $GENVERSION_NAME[$iver];	
    my $GENSTAT_FILE  = &get_FILENAME_GENSTAT($iver);
    my $README_FILE   = "$PATH_SNDATA_SIM/${GV}/${GV}.README" ;
    qx(mv $GENSTAT_FILE $README_FILE); 
} # end  move_GENSTAT_FILE

# =========================
sub write_GENSTAT_table {

    my ($IVER,$OPT) = @_;

    # Apr 1 2019: 
    # write GENSTAT table using 'sed' to insert at mark in file.
    # OPT=1 -> total sums with links
    # OPT=2 -> actual generated sums without links

    my $GENVERSION = $GENVERSION_NAME[$IVER] ;
    my $STATFILE   = &get_FILENAME_GENSTAT($IVER) ;
    my $sedCmd = "sed " ;
    my @OUTLINES = ();
    my $NL = 0 ;

    my ($m, $Nm, $COMMENT, $cmodel, $ROWKEY, $INSERT_MARK );
    my ($NGEN, $NWR, $NSPEC, $LINE, $iline );

    if ( $OPT == $OPT_GENSTAT_TABLE_plusLINKS ) { 
	$ROWKEY="JOBSUM:"; 
	$COMMENT = "JOB-SUMS";
	if ( $NMODEL_DUPL_TOT ) 
	{ $COMMENT = "$COMMENT " . " with LINKS (for analysis)" ; }
    }
    else { 
	$ROWKEY="jobsum:"; 
	$COMMENT = "job-sums without LINKS (for CPU monitor)"; 
    }
    $INSERT_MARK = "$TABLE_INSERT_MARK${OPT}:" ;

    # print table header
    $NL++; $OUTLINES[$NL] = "";
    $NL++; $OUTLINES[$NL] = "Updat from HOST=$HOST";
    $NL++; $OUTLINES[$NL] = "$GENVERSION $COMMENT";
    $NL++; $OUTLINES[$NL] = 
	"#       model           NGEN    NWRITE    NSPEC " ;
    $NL++; $OUTLINES[$NL] = 
	"# -----------------------------------------------" ;

    $Nm = $NGENMODEL[$IVER]; 
    for ( $m=0; $m < $Nm+2; $m++ ) {

	if ( $OPT == $OPT_GENSTAT_TABLE_plusLINKS ) {
	    $NGEN  = $GENSTAT_SUMTOT[$IVER][$m][0];
	    $NWR   = $GENSTAT_SUMTOT[$IVER][$m][1];
	    $NSPEC = $GENSTAT_SUMTOT[$IVER][$m][2];
	}
	else {
	    $NGEN  = $GENSTAT_SUMJOB[$IVER][$m][0];
	    $NWR   = $GENSTAT_SUMJOB[$IVER][$m][1];
	    $NSPEC = $GENSTAT_SUMJOB[$IVER][$m][2];
	}

	if ( $m < $Nm ) 
	{ $cmodel = sprintf("%-10d", $m); }
	elsif ( $m == $Nm ) 
	{ $cmodel = sprintf("%-10s", "SUM-SNIa" ); }
	elsif ( $m == $Nm+1 ) 
	{ $cmodel = sprintf("%-10s", "SUM-NONIa" ); }
	  
	$LINE = sprintf("%s %s %9d %8d %8d", 
			$ROWKEY, $cmodel, $NGEN, $NWR, $NSPEC);
	$NL++ ; $OUTLINES[$NL] = "$LINE" ;
    }

    $NL++ ; $OUTLINES[$NL] = "";

    # prepare sed command to insert summed table in alread-existing
    # GENSTAT file
    my $TEMPFILE = "TEMP_${GENVERSION}.GENSTAT";
    for($iline=1; $iline <= $NL; $iline++ ) {
	$LINE = "$OUTLINES[$iline]";
	$sedCmd = "$sedCmd " . "-e '/$INSERT_MARK/a\\\b$LINE' " ;
    }    

    $sedCmd = "$sedCmd " . "-e '/$INSERT_MARK/d' " ; 

    print "     Update GENSTAT for $GENVERSION  (OPT=$OPT) \n";
    qx($sedCmd $STATFILE > $TEMPFILE; mv $TEMPFILE $STATFILE);

    return ;

} # end write_GENSTAT_table

# ========
sub get_FILENAME_GENSTAT {
    my($iver) = @_ ;
    my $GENVERSION = "$GENVERSION_NAME[$iver]" ;
    my $STATFILE   = "$LOGDIR/GENSTAT_${GENVERSION}.README" ;
    return($STATFILE);
}

# =============================================
sub printLegend_INFILE {
    my($IVER,$FH) = @_ ;  
    # print $m  $CLASS  $INFILE
    my ($m, $NAME, $INFILE, $CLASS, $LINE);
    print $FH "\n\n#     model  CLASS          INFILE \n";
    print $FH "# ------------------------------------------------ \n";
    for ( $m=0; $m < $NGENMODEL[$IVER]; $m++ ) {  
	$CLASS  = $GENMODEL_CLASS[$IVER][$m] ;
	$INFILE = $SIMGEN_INFILE[$IVER][$m] ;
	$LINE   = sprintf("INFILE: %2d   %-12s  %s", $m, $CLASS, $INFILE);
	print $FH "$LINE\n" ;
    }
    $FH -> autoflush(1);
} # end printLegend_INFILE

# =============================================
sub printLegend_GENOPT {
    my($IVER,$FH) = @_ ;  
    my ($m, $CLASS, $GENOPT, $LINE);
    # print $m  $CLASS  $GENOPT
    print $FH "\n\n#     model   CLASS         GENOPT    \n";
    print $FH "# ----------------------------------------------------- \n";
    for ( $m=0; $m < $NGENMODEL[$IVER]; $m++ ) { 
	$CLASS  = $GENMODEL_CLASS[$IVER][$m] ;    
	$GENOPT = "$GENVERSION_GENOPT[$IVER][$m]" ;
	if ( "$GENOPT" eq '' ) { $GENOPT = "none"; }
	$LINE   = sprintf("GENOPT: %2d   %-12s   %s", $m, $CLASS, $GENOPT );
	print $FH "$LINE\n" ;
    }
    $FH -> autoflush(1);
} # end printLegend_GENOPT

# ===============================================
sub make_AUXFILE_FILTERS {

    my($iver,$jobid) = @_;
    print "   Create AUXILIARY FILTER sdir for $VERSION_FINAL \n";

    my $fDir_temp  = "${SIMGEN_TEMPDIR}/${VERSION_TEMP}.FILTERS" ;
    my $fDir_final = "${SIMGEN_FINALDIR}/${VERSION_FINAL}.FILTERS" ;
    my $cmd = "cp -r $fDir_temp $fDir_final";
    qx($cmd);

} # end of make_AUXFILE_FILTERS


# ==============================
sub make_DONEFILE {
    my ($iver,$jobid) = @_ ;

    # create done file for this job(s).
    # One done file per Ia+CC mix.
    # If job aborted then apend _ABORT as part of done-file name.

    my ($GENVERSION, $str_ver, $str_job, $DONEFILE, $SUFFIX );

    $str_ver  = &get_jobidString($iver);
    $str_job  = &get_jobidString($jobid);

    $GENVERSION = "$GENVERSION_NAME[$iver]" ;

    $SUFFIX = "" ;
    if ( $NRANJOB_SIMGEN > 1 ) 
    { $SUFFIX = "-$str_job" ; }

    if ( $JOBABORT_FLAG      ) 
    { $SUFFIX = "${SUFFIX}_ABORT" ; }    

    $DONEFILE   = "$LOGDIR/SIMJOB${str_ver}_${GENVERSION}${SUFFIX}.DONE" ;

    qx(touch $DONEFILE);

} # end of make_DONEFILE

	    


# ========================
sub check_VersionDone{

    # Jan 24 2018
    # Returns 1 if this IVER is all done.
    # Use LOCK file to avoid conflicts when
    # two batch jobs both think they are done.
    
    my ($IVER) = @_ ;

    my $doneFlag = 0;  # init output to not done.
    my (@tmp, $str_ver, $doneList, $NDONE_FOUND, $NDONE_EXPECT );
    my ($LOCKFILE_DONE, $VERSION);
    
    $str_ver  = &get_jobidString($IVER);
    $doneList = "$LOGDIR/SIMJOB${str_ver}*.DONE" ;
    
    @tmp = qx(ls $doneList  2>/dev/null );
    $NDONE_FOUND  = scalar(@tmp);
    $NDONE_EXPECT = $NRANJOB_SIMGEN ;

    if ( $NDONE_FOUND == $NDONE_EXPECT ) {
	$VERSION = $GENVERSION_NAME[$IVER] ;
	$LOCKFILE_DONE = "$LOGDIR/LOCK_DONE_${VERSION}" ;

# May 9 2018: version-dependent sleep to help avoid conflict
	my $num_usec = $NODEINDX*(1000000/10); # .1 sec per nodeindx
	usleep($num_usec);
#	sleep($NODEINDX);

	# set done flag only if the LOCK file does not already exist.
	if ( ! ( -e $LOCKFILE_DONE ) ) {
	    open  PTR_LOCKFILE , "> $LOCKFILE_DONE" ;
	    close PTR_LOCKFILE ;
	    $doneFlag = 1;
	}
    }
    
    return($doneFlag);

} # end check_VersionDone{

# ==============================
sub check_ALLDONE {

    # Created Jun 7 2013 by RK
    # when all of the DONE files have been created then
    # - create SUMMARY.LOG file
    # - create global DONE stamp if DONE_STAMP: key is specified.
    # - make tarFile of logs, readme, .batch, etc ...
    #
    # Apr 3 2019: add BUSY FILE logic to avoid conflict.

    my ($doneList, $BUSUYFILE);
    my (@tmp, $NDONE_FOUND, $NDONE_EXPECT, $IVER, $JOBID );
    my $DONE_STATUS = "FAILED";

    $doneList = "$LOGDIR/SIMJOB*.DONE" ;
    
    @tmp = qx(ls $doneList  2>/dev/null );
    $NDONE_FOUND  = scalar(@tmp);
    $NDONE_EXPECT = $NGENVERSION  * $NRANJOB_SIMGEN ;


    if ( $NDONE_FOUND == $NDONE_EXPECT ) {

	sleep(2);

	# use BUSY file to avoid conflict from multiple batch scripts 
	my $BUSYFILE = "$LOGDIR/BUSY_ALLDONE_TASKS" ;
	if ( -e $BUSYFILE ) { return ; }
	open  PTR_BUSY  , "> $BUSYFILE" ;

	&final_SIMGEN_DUMP_duplicate();

	for($IVER=1; $IVER<= $NGENVERSION; $IVER++ )  { 
	    &convert_SIMGEN_DUMP($IVER);  # optional DUMP -> HBOOK or ROOT
	} 

	&write_GENSTAT_DRIVER(-1,-1); # make sums for each version
	&make_TOTAL_SUMMARY_FILE();

	# remove BUSY file before tarring things up
	close PTR_BUSY ;   qx(rm $BUSYFILE);

	# add DONE stamp, and make tar file in SIMLOGS
	if ( $DONE_STAMP_FLAG )	 {
	    qx(touch $DONE_STAMP) ;  
	    if($NABORT==0) {&makeTar_SIMLOGS(); $DONE_STATUS="SUCCESS"; }
	    qx(echo $DONE_STATUS >> $DONE_STAMP); 
	}
    }

} # end of check_ALLDONE


# ===============================
sub make_TOTAL_SUMMARY_FILE {

    # Created Apr 1 2019
    # For each version, write total stats for Ia and NONIa

    sleep 4 ;  # needed for very short jobs
    my $TOTAL_FILE   = "$LOGDIR/TOTAL_SUMMARY.LOG" ;
    my @TYPELIST     = ( "SNIa ", "NONIa" );
    my ($JJ, $SUMSTRING, $FLIST, @tmp, $tmpLine, @wdlist);
    my ($jslash, $vlen, $k0,$k1 );
    my($iv, $m, $GV, $gv, $LINE, $TYPE, $NGEN, $NWR, $NSPEC);
    my $format = "%-40s  %8d %7d %6d" ;

    print "\n Created $TOTAL_FILE \n";

    open  PTR_TMPFILE , "> $TOTAL_FILE" ;

    print PTR_TMPFILE "PATH_SNDATA_SIM: $PATH_SNDATA_SIM \n\n";

    for($JJ=0 ; $JJ < 2; $JJ++ ) {
	# JJ=0 for SNIa, JJ=1 for NONIa
	$TYPE = $TYPELIST[$JJ] ;     $SUMSTRING = "SUM-$TYPE" ;

	print PTR_TMPFILE
	    "                                                $SUMSTRING : \n" ;
	print PTR_TMPFILE 
	    "  GENVERSION                                  " . 
	    "NGEN     NWRITE  NSPEC   \n";
	print PTR_TMPFILE 
	    "# -----------------------------------------" . 
	    "------------------------ \n";

	if ( $SAMEFLAG_RANSEED  ) { 

	    # random subsets are combined; use GENSTAT sums
	    for($iv=1; $iv <= $NGENVERSION; $iv++ ) {
		$m     = $NGENMODEL[$iv] + $JJ ;  # SUM-TYPE index
		$GV    = $GENVERSION_NAME[$iv] ;
		$NGEN  = $GENSTAT_SUMTOT[$iv][$m][0] ;
		$NWR   = $GENSTAT_SUMTOT[$iv][$m][1] ;
		$NSPEC = $GENSTAT_SUMTOT[$iv][$m][2] ;
		$LINE = sprintf("$format", $GV,$NGEN,$NWR,$NSPEC); 
		print PTR_TMPFILE "$LINE\n";
	    }
	}
	else {
	    # random subsets remain separate; grep SUMs from README files.
	    for($iv=1; $iv <= $NGENVERSION; $iv++ ) {
		$GV    = $GENVERSION_NAME[$iv] ;
		$FLIST = "$PATH_SNDATA_SIM/${GV}-*/${GV}-*.README" ;
		@tmp   = qx(grep '$SUMSTRING' $FLIST);
		foreach $tmpLine (@tmp) {
		    @wdlist  = split(/\s+/,$tmpLine ) ; 
		    $k0      = rindex($wdlist[0],"/") + 1 ;
		    $k1      = index($wdlist[0],".README");
#		    $k0 = $jslash+1; $k1=$vlen-$k0-1 ;
		    $gv      = substr($wdlist[0],$k0,$k1-$k0); # genversion
		    $NGEN    = $wdlist[3] ;
		    $NWR     = $wdlist[4] ;
		    $NSPEC   = $wdlist[5] ;
		    $LINE = sprintf("$format", $gv,$NGEN,$NWR,$NSPEC); 
		    print PTR_TMPFILE "$LINE\n";
		}
	    }
	}

	print PTR_TMPFILE "\n\n" ;
	PTR_TMPFILE -> autoflush(1);
    }
	
    close PTR_TMPFILE ;
    return ;

} # end make_TOTAL_SUMMARY_FILE   



# ================================
sub makeTar_SIMLOGS {

    # create tarFile of logs and readme's under /SIMLIBS_[PREFIX].
    # then delete files.
    # Leave ALL.DONE and SIMJOB_SUMMARY.LOG so that they
    # always remain visible.

    my ($FLIST, $IVER, $VERSION, $TARFILE );

    $FLIST = "" ;
    if ( $NCPU > 0 ) { $FLIST = "${GENPREFIX}* " ; }  # batch files
    for ( $IVER=1 ; $IVER <= $NGENVERSION; $IVER++ ) {
	$VERSION = $GENVERSION_NAME[$IVER];
	$FLIST = "$FLIST *${VERSION}* " ;  # tack on log for each version
    }

    $TARFILE = "SIMLOGS.tar";
    print "\n Make tarFile from files = \n $FLIST \n\n";
    
    qx(cd $LOGDIR; tar -cf $TARFILE $FLIST );
    qx(cd $LOGDIR; gzip    $TARFILE  );
    sleep(1);
    qx(cd $LOGDIR; rm $FLIST  );    

} # end makeTar_SIMLOGS

# ================================
sub cleanup {

    # remove temp directories

    my($iver,$jobid) = @_;

    if ( $CLEANUP_FLAG == 0 ) { return ; }

    if ( length($SIMGEN_TEMPDIR) > 2 ) {
	my $cmd = "rm -r $SIMGEN_TEMPDIR" ;
	qx($cmd);
    }

} # end of cleanup


# ======================================
sub combine_random_simjobs {

    # Apr 30 2014 R.Kessler
    #
    # If multiple RANSEEDs are the same, then automatically
    # combine the jobs into a single version ... and delete
    # the individual RANSEED directories.

    my ($IVER) = @_  ;

    # continue only if we have multiple RANSEED jobs,
    # AND all RANSEEDs are the same !
    # If the ranseeds are different then leave the separate jobs alone.
    #
    # Aug 11 2017: tar the misc/ subdir.
    # Dec 01 2017: remove duplicate NVAR and VARNAME keys
    # May 10 2018: include NSPEC in summary
    # Jan 21 2019: write NVAR line only if it exists in original dump file
    #
    # --------- BEGIN ----------

    if ( $NRANJOB_SIMGEN   <= 1 ) { return ; }
    if ( $SAMEFLAG_RANSEED == 0 ) { return ; }

    my ($VERSION_ALL, $SIMDIR_ALL, $FOUND_NVAR_KEY);
    my (@VERSION_RAN, $SIMDIR_RAN, $VRAN, $cjob, $jobid, @FLIST);
    my ($mvList, $ls, @lsList, $iList, @summLines );
    my (@tmp, $iver, $m, $idup, $star,  $txtLine1, $txtLine2);

    $VERSION_ALL = "$GENVERSION_NAME[$IVER]" ;

    print "\n !!! COMBINE RANDOM SIMJOBS for $VERSION_ALL !!! \n" ;

    $SIMDIR_ALL  = "$PATH_SNDATA_SIM/$VERSION_ALL" ;
    if ( -d $SIMDIR_ALL ) { qx(rm -r $SIMDIR_ALL); }
    qx(mkdir $SIMDIR_ALL);

    my $CD = "cd $SIMDIR_ALL" ;

    qx($CD ; mkdir $MISC_SDIR );

    for($jobid=1; $jobid <= $NRANJOB_SIMGEN; $jobid++ ) {
	$cjob    = &get_jobidString($jobid);
	$VRAN    = "${PREFIX_TEMP}_${VERSION_ALL}-$cjob" ;
	$VERSION_RAN[$jobid] = "$VRAN" ;
    }

    # -------------------------------------------
    @lsList = (); $iList=0;

    if ( $WRFLAG_FITS ) {  
	$mvList = "${GENPREFIX}*.FITS" ;  

	if ( $DOGEN_SNIa  ) 
	{ $lsList[$iList] = "${GENPREFIX}*SNIa*HEAD.FITS*" ;  $iList++ ;} 

	if ( $DOGEN_NONIa ) 
	{ $lsList[$iList] = "${GENPREFIX}*NONIaMODEL*HEAD.FITS*" ; $iList++; }

    }
    if ( $WRFLAG_TEXT ) { 
	$mvList    = "${GENPREFIX}*.DAT"  ;  
	$lsList[0] = "${GENPREFIX}*.DAT" ; 
    }

    # -------------------------------------------------------------
    # loop over each random version and move contents to $VERSION

    for($jobid=1; $jobid <= $NRANJOB_SIMGEN; $jobid++ ) {
	$VRAN          = $VERSION_RAN[$jobid] ;
	$SIMDIR_RAN    = "$PATH_SNDATA_SIM/$VRAN" ;

	if ( !(-d $SIMDIR_RAN) ) {
	    print "\n\n ******** SEVERE WARNING ********** \n";
	    print "  Missing SIMDIR_RAN=\n  $SIMDIR_RAN \n\n";
	    next ;
	}

	# move the data files
	qx(cd $SIMDIR_RAN ; mv $mvList $SIMDIR_ALL/ );

	# create symbolic links for duplicate jobs ... 
	# which may still be running
	if ( $NMODEL_DUPL[$IVER] > 0 ) {
	    for($m=0; $m<$NGENMODEL[$IVER]; $m++ ) 
	    { &make_symLinks_duplicate($IVER,$m,$jobid); }
	}

	# move the misc/ files too (note some are redundant)
	qx(cd $SIMDIR_RAN/$MISC_SDIR ; mv *.* $SIMDIR_ALL/$MISC_SDIR) ;

	if ( $DO_SIMGEN_DUMP  ) { 
	    &move_random_SIMGEN_DUMP($IVER,$SIMDIR_RAN,$VRAN);
	}
    }  # end jobid


    # -------------------------------------------------
    # ------------------ AUX FILES --------------------
    # -------------------------------------------------
    

    my ($LIST_FILE, $IGNORE_FILE, $README_FILE, $DUMP_FILE);
    $LIST_FILE      = "${VERSION_ALL}.LIST" ;
    $IGNORE_FILE    = "${VERSION_ALL}.IGNORE" ;
    $DUMP_FILE      = "${VERSION_ALL}.DUMP" ;
    # note that README is created later after all jobs finish

    qx($CD; touch $IGNORE_FILE; touch $LIST_FILE );	

    # create combined LIST file; use ordered ls to make sure that
    # SNIa appear first in LIST file so that SIM_XXX params in the
    # output correspond to SNIa instead of NONIa.
    foreach $ls (@lsList)  { qx($CD ;  ls $ls >> $LIST_FILE); }

    # if there are duplicates, remove .gz extensions from list file
    if ( $NMODEL_DUPL_TOT > 0 ) {
	my $sedcmd = "sed -i 's/\.gz//g' $LIST_FILE" ;
	qx($CD ; $sedcmd);
    }


    # - - - - - - 
    if ( $DO_SIMGEN_DUMP )  { &cat_random_SIMGEN_DUMP($IVER); }


    # ------------
    # gzip stuff ... 
    &compress_output($IVER,0);

    # ------------------
    # remove RANSEED versions since we now have the combined version

  REMOVE_SIMDIRS:
    if ( $CLEANUP_FLAG ) {
	for($jobid=1; $jobid <= $NRANJOB_SIMGEN; $jobid++ ) {
	    $VRAN         = $VERSION_RAN[$jobid] ;
	    if ( $VRAN eq "" ) { next; }
	    $SIMDIR_RAN   = "$PATH_SNDATA_SIM/$VRAN" ;
	    qx(rm -r $SIMDIR_RAN); 
	}
    }

    return ;

} # end of combine_random_simjobs 

# ===========================================
sub move_random_SIMGEN_DUMP {

    # move DUMP file from input SIMDIR & VERSION to
    # final combined SIMDIR_ALL.

    my ($IVER, $SIMDIR, $VERSION) = @_ ;

    my $VERSION_ALL = "$GENVERSION_NAME[$IVER]" ;
    my $SIMDIR_ALL  = "$PATH_SNDATA_SIM/$VERSION_ALL" ;
    my $tmpDumpFile = "$SIMDIR_ALL/${VERSION}.${SUFFIX_DUMP_TEMP}" ;
    my ($m, $MM );

    # move DUMP file to .DUMP_TEMP 
    qx(cd $SIMDIR ; mv *.DUMP $tmpDumpFile);
	
    # check for MODEL-dependent DUMP files needed later
    # for duplicate versions
    for($m=0; $m<$NGENMODEL[$IVER]; $m++ ) {
	if ( $IVERSION_DUPLICATE_INVERSE[$IVER][$m] > 0 ) {
	    $MM = "MODEL$m" ;
	    qx(cd $SIMDIR ; mv *.DUMP_$MM ${tmpDumpFile}_${MM} );
	}
    }       

    return ;

} # end combine_random_SIMGEN_DUMP


# ===========================
sub cat_random_SIMGEN_DUMP {

    # catenate SIMGEN_DUMP files from random simjobs.

    my($IVER) = @_;

    
    my $VERSION_ALL   = "$GENVERSION_NAME[$IVER]" ;
    my $SIMDIR_ALL    = "$PATH_SNDATA_SIM/$VERSION_ALL" ;
    my $DUMP_FILE     = "${VERSION_ALL}.DUMP" ;
    my $CD            = "cd $SIMDIR_ALL" ;
    my ($m, $FOUND_NVAR_KEY );

    # catenate the DUMP files
    qx($CD ; cat *.${SUFFIX_DUMP_TEMP} > $DUMP_FILE);
    if($CLEANUP_FLAG) { qx($CD ; rm  *.${SUFFIX_DUMP_TEMP} ); }

    # catenate model-subsets to use later for duplicates
    for($m=0; $m<$NGENMODEL[$IVER]; $m++ ) {
	if ( $IVERSION_DUPLICATE_INVERSE[$IVER][$m] > 0 ) {
	    my ( $MM, $tmp_SUFFIX, $catCmd );
	    $MM         = "MODEL$m" ;
	    $tmp_SUFFIX = "${SUFFIX_DUMP_TEMP}_${MM}" ;
	    $catCmd     = "cat *.${tmp_SUFFIX} > ${DUMP_FILE}_${MM}" ;
	    qx($CD ; $catCmd ) ;
	    if($CLEANUP_FLAG) { qx($CD ; rm  *.${tmp_SUFFIX} ); }
	}
    }


    # Dec 1 2017:
    # Cleanup cat-DUMP file so that NVAR and VARNAMES appear only once.
    $FOUND_NVAR_KEY=0 ;
    my($line_NVAR, $line_VARNAMES, $sedcmd, @tmp);
    @tmp           = qx($CD; grep "NVAR: "     $DUMP_FILE);
    if ( scalar(@tmp) > 0 ) { $line_NVAR = "$tmp[0]"; $FOUND_NVAR_KEY=1; }

    @tmp           = qx($CD; grep "VARNAMES: " $DUMP_FILE);
    $line_VARNAMES = "$tmp[0]" ;    
#    $line_VARNAMES =~ s/\s+$// ;  # remove trailing blanks

    $sedcmd = "sed " ;
    if ( $FOUND_NVAR_KEY ) { 
	$sedcmd = "$sedcmd " . "-e '1i $line_NVAR'" ;
	$sedcmd = "$sedcmd " . "-e '/NVAR/d'"; 
    }
    $sedcmd = "$sedcmd " . "-e '1i $line_VARNAMES'" ;
    $sedcmd = "$sedcmd " . "-e '/VARNAMES/d'";
    qx($CD ; $sedcmd $DUMP_FILE > BLA.DUMP; mv BLA.DUMP $DUMP_FILE);

    return ;

} # end cat_random_SIMGEN_DUMP

# ===========================================
sub make_symLinks_duplicate {
    my($iverLink,$mLink,$jobid) = @_;

    # Created Jan 2019
    # Create FITS file symbolic link for duplicates.
    # Determine HEAD and PHOT filenames because if iverReal
    # jobs are still running then ls won't work.

    my $iverReal     = $IVERSION_DUPLICATE[$iverLink][$mLink]; 
    my $mReal        = $IMODEL_DUPLICATE[$iverLink][$mLink]; 
    if ( $iverReal         <  0 ) { return ; }
    if ( $mReal            <  0 ) { return ; }
    if ( $SAMEFLAG_RANSEED == 0 ) { return ; }

    my $VNAME_LINK   = "$GENVERSION_NAME[$iverLink]" ;
    my $VNAME_REAL   = "$GENVERSION_NAME[$iverReal]" ;
    my $SIMDIR_LINK  = "$PATH_SNDATA_SIM/$VNAME_LINK" ;
    my $SIMDIR_REAL  = "$PATH_SNDATA_SIM/$VNAME_REAL" ;
    my $CLASS_LINK   = $GENMODEL_CLASS[$iverLink][$mLink] ;
    my $CLASS_REAL   = $GENMODEL_CLASS[$iverReal][$mReal] ;

    # make symLink for each FITS file in this model class.

    my $cjob       = &get_jobidString($jobid) ; 
    my $fReal_HEAD = "${GENPREFIX}_${CLASS_REAL}-${cjob}" . "_HEAD.FITS.gz" ; 
    my $fReal_PHOT = "${GENPREFIX}_${CLASS_REAL}-${cjob}" . "_PHOT.FITS.gz" ; 
    my $fLink_HEAD = "${GENPREFIX}_${CLASS_LINK}-${cjob}" . "_HEAD.FITS.gz" ; 
    my $fLink_PHOT = "${GENPREFIX}_${CLASS_LINK}-${cjob}" . "_PHOT.FITS.gz" ; 

    qx(cd $SIMDIR_LINK; ln -s ../$VNAME_REAL/$fReal_HEAD $fLink_HEAD ) ;
    qx(cd $SIMDIR_LINK; ln -s ../$VNAME_REAL/$fReal_PHOT $fLink_PHOT ) ;

    # -------------------------------------------
    # -------------------------------------------
    # Feb 4 2019:
    # Legacy system: create new auxilliary file $VERSION.SYMLINK 
    # containing version name where symLinks are pointing.
    # This file is checked by split_and_fit so that fitting
    # can be skipped as well and FITRES outputs made with symLinks.
    #
    # This legacy system is planned to be replaced with a newer 
    # split_and_fit system that works for mixed Ia+CC sims.
    #
    if ( $DOGEN_SNIa>0 && $DOGEN_NONIa==0 ) {
	my $tmpFile = "${SIMDIR_LINK}/${VNAME_LINK}.SYMLINK" ;
	qx(echo $VNAME_REAL > $tmpFile);
    }
 

} # end make_symLinks_duplicate


# ==============================
sub compress_output {

    # Created Jan 23 2018 [code moved from combine_random_simjobs ]
    #  * /misc -> misc.tar
    #  * gzip FITS file
    #  * gzip DUMP file (if dump-all option)

    my ($IVER,$JOBID) = @_ ;

    my($CD, $SIMDIR_OUT, $VERSION, $DUMP_FILE, $CJOBID);
    my($tarFile, $tarCmd, $fullCmd, $gzipCmd );

    # -------- BEGIN ---------

    if ( $NRANJOB_SIMGEN   < 1 ) { return ; }

    $VERSION     = "$GENVERSION_NAME[$IVER]" ;

    if ( $SAMEFLAG_RANSEED == 0 ) {
	# each job is preserved with -[jobid] extension to versions
	$CJOBID      = &get_jobidString($JOBID);
	if ( $NRANJOB_SIMGEN > 1 && $SAMEFLAG_RANSEED == 0 ) 
	{ $VERSION = "$VERSION" . "-" . "$CJOBID" ;  }
    }

    $SIMDIR_OUT  = "$PATH_SNDATA_SIM/$VERSION" ;
    $DUMP_FILE   = "${VERSION}.DUMP" ;
    $CD          = "cd $SIMDIR_OUT" ;

    # Aug 11 2017: tar up the misc dir to save on file count quotas
    $tarFile = "${MISC_SDIR}" . ".tar" ;
    $tarCmd  = "tar -cf $tarFile $MISC_SDIR; gzip $tarFile" ;
    $fullCmd = "$CD ; $tarCmd ; rm -r $MISC_SDIR " ;
    print " tar /misc with command\n $fullCmd \n";
    qx($fullCmd) ;

    # for SIMGEN_DUMPALL, this file could be HUUUUGE, so gzip it.
    if ( $NCLASS_SIMGEN_DUMPALL > 0 ) { qx($CD ; gzip $DUMP_FILE ); }

    # Jan 23 2018: gzip FITS file
    if ( $WRFLAG_FITS ) { 
	$gzipCmd = "gzip *.FITS" ;
	$fullCmd = "$CD ; $gzipCmd " ;
	print " gzip *.FITS with command\n $fullCmd \n";
	qx($fullCmd) ;
    }


    return ;

}  # end compress_output


# ==================================
sub convert_SIMGEN_DUMP {
    my ($IVER) = @_  ;

    # Aug 2017: check option to convert SIMGEN_DUMP file into HBOOK or ROOT

    my $VERSION_ALL = "$GENVERSION_NAME[$IVER]" ;
    my $SIMDIR_ALL  = "$PATH_SNDATA_SIM/$VERSION_ALL" ;
    my $DUMP_FILE   = "${VERSION_ALL}.DUMP" ;
    my $CD          = "cd $SIMDIR_ALL" ;

    my($FMTARG, $CMD, $ARGS);

    if ( $DO_SIMGEN_DUMP == 0 ) { return ; }

    $FMTARG="" ;  
    if ( $CONVERT_SIMGEN_DUMP eq 'HBOOK' ) { $FMTARG='H'; }
    if ( $CONVERT_SIMGEN_DUMP eq 'ROOT'  ) { $FMTARG='R'; }

    if ( $FMTARG ne '' ) {
	$ARGS = "$DUMP_FILE $FMTARG -outprefix $VERSION_ALL" ;
	$CMD  = "combine_fitres.exe $ARGS" ;
	$CMD  = "$CMD" . " ; rm ${VERSION_ALL}.text" ;
	print " Convert $DUMP_FILE file to $CONVERT_SIMGEN_DUMP \n" ;
	qx($CD ; $CMD );
    }

} # end convert_SIMGEN_DUMP
