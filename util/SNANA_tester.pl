#!/usr/bin/perl
#
# Created June, 2011 by R.Kessler
#
# Usage:
#   SNANA_tester.pl
#   SNANA_tester.pl --TESTDIR <full path of private dir>
#   SNANA_tester.pl KILL     ! kill job on each node)
#   SNANA_tester.pl help     ! print info on adding new jobs)
#   SNANA_tester.pl LCPLOT   ! process LCPLOT*.nml files, make pdf plots
#                            ! --> runs interactively, no ssh
# Note that there are no required argumements;
# options are made via command line or by modifying
# the input files below.
#
# Output in $TESTDIR/logs/${cdate} :
#  *  MASTER.LOG  with summary
#  *  RESULTS_[OLD,NEW].log with parsed results
#  *  .LOG  for each test job and each snana version
#
#
# This script compares application results for two different
# versions of snana, mainly to make sure that a new snana version
# reproduces results from the previous version, and to identify
# cases where the new snana version gives different results.
#
# $SNANA_TESTS/SNANA_TESTS.LIST is the top-level file
# containing
#
#   OLD_SNANA_SETUP: setup  snana v6_16
#   NEW_SNANA_SETUP: setup  snana
#
#   TEST: TEST1
#   TEST: TEST2
#   etc ...
#
# where TEST[1,2 ..] are test-files that contain the
# following:  Note that simulation jobs should begin
# with "SIMGEN" because the SIMGEN jobs are run first
# in case there are fit-jobs that use the simulation
# output.
#
#  TESTJOB: <jobname>   [such as snlc_fit.exe or snlc_sim.exe]
#  TESTINPUT: <required input filename>  <optional 2ndary input>
#  TESTJOB_ARGS:  <optional commandLine args>
#  TESTNAME: <name of quantity to test>
#  TESTRESULT:  <grep command> to extract result-line for comparisons
#  WORDNUM:   <n1-n2>    specify range of words from  TESTRESULT above
#
#
#   HISTORY
#  ~~~~~~~~~~~
#
# Aug 02, 2011: use strict
#
# Dec 10, 2011: allow new argument in test-file,
#               TESTJOB_ARGS: a b c d TESTINPUT e f g ...
#               to give more complex instructions.
#               TESTINPUT will be replaced with NEW###_$TESTINPUT
#               and OLD###_$TESTINPUT. This fix is needed to test
#               SIMSED_extractSpec.exe that uses command-line args
#               instead of an input file.
#
# Aug 16, 2012: new command-line option 'help' to print info on
#                how to add more tests.
#
# Jan 1, 2013: 
#    add NEWFITS option to change kcor*.his files into kcor*.fits
#    files (to test new fits format for kcor files).
#
#    Allow multiple TESTINPUT: lines.
#
# Sep 21, 2103: fix run_TEST() so that blank results are flagged
#               and OLD and NEW blank-tests show as a failure
#               instead of a success.
#
# Mar 2014: print 3-digit ID in summary files
#
# May 14, 2014: new LCPLOT option ; see LCPLOT_DRIVER
#
# Jun 4 2014: fix to work with same NODE repeated in NODELISTs.
#             --> replace --NODE arg with --INODE argument.
#
#
# Jun 28 2016: copy input file only if it doesn't exist.
#              --> avoids error messages.
#
# July 20 2016: remove kcor_his2fits().
#
# May  07 2019: 
#   + read task files in $TESTDIR/tasks. 
#   + rename SNANA_TESTS.LIST -> SNANA_CodeTests_FNAL.LIST since this scheme
#     works only at FNAL
#   + this minor refactor is to prepare for a full python re-write
#     using same directories, but different LIST file that specifies
#     batch or SSH nodes.
#
# May 13 2019:
#   + remove leading blanks from grepped result string (to match python)
#   + -29s -> -40s format for RESULTS
#
# =============================================

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# ---------------------------
# declare subs

sub set_defaults;
sub parse_args;
sub init_stuff;
sub parse_LISTFILE ;
sub sort_TESTLIST;
sub assign_NODE ;
sub make_LOGDIR ;
sub runScript ;
sub parse_TEST ;
sub prep_TEST;
sub wait_for_TEST;
sub run_TEST;
sub analyze_RESULTS ;
sub print_help ;

sub LCPLOT_DRIVER ;

# ----------------------------
# globals

my ($SNANA_DIR, $SNDATA_ROOT, @SNANA_SETUP, $CDATE, $scriptName );
my (@OLD_SNANA_SETUP, @OLD_SNANA_NODELIST, $NNODE_OLD, $INDEX_OLD ) ;
my (@NEW_SNANA_SETUP, @NEW_SNANA_NODELIST, $NNODE_NEW, $INDEX_NEW ) ;
my (@NNODE, @NODELIST, $INODE_SSH, $SSHNODE_LOGFILE); 
#my $SSHNODE
my (@TESTLIST_NODE, @TESTLIST_INODE, @TESTLIST_IDCHAR);
my ($OLDFITS, $NEWFITS);
my (@TESTLIST_NAME, $KILLJOBS_FLAG, $PROC_VERSION, @OLDNEW_KEY );
my (@VERS_SUFFIX, $IVERSION, $SIMGEN_PREFIX, $SNTABLE_DUMP_PREFIX, $LISTFILE );
my ($NTEST_ALL, $NTEST_SIMGEN, $NTEST_NOTSIM);
my ($NTEST_SNTABLE, @SNTABLE_FLAG);
my ($TESTDIR, $TASKDIR, $TESTLOGDIR, $TESTINPUTDIR );
my ($TESTJOB, @TESTINPUT,  $TESTJOB_ARGS, $TESTNAME, $TESTRESULT, $WORDNUM );
my ($MASTER_LOGFILE, $q, $qq, $OPT_ABORT, $OPT_WARN, $OPT_QUIET );
my ($TEST_INFILE, $TEST_LOGFILE, $TEST_DONEFILE, $TEST_ARGLIST );
my ($TEST_WORDNUM_FIRST, $TEST_WORDNUM_LAST ) ;
my ($OPT_LCPLOT);

# --------------- BEGIN MAIN ----------------
my ($iver,$inode,$itest,$node);


&set_defaults();

# parse command-line args
&parse_args();

if ( $OPT_LCPLOT ) { &LCPLOT_DRIVER(); }

&init_stuff();

# read list of tests from LISTFILE
&parse_LISTFILE();

if ( $KILLJOBS_FLAG ) { 
    sntools::Killjobs(@OLD_SNANA_NODELIST) ; 
    sntools::Killjobs(@NEW_SNANA_NODELIST) ; 
    die "\n" ;
}

# sort tests with SIMGEN first, then the others
&sort_TESTLIST();

# assign node for each TEST
&assign_NODE();

if ( $PROC_VERSION eq ' ' ) {

    &make_LOGDIR();

    print "\n";
    for ( $iver = 0 ; $iver <= 1 ; $iver++ ) {
	for ( $inode=0; $inode < $NNODE[$iver]; $inode++ ) {
	    &runScript($iver,$inode);
	}
    }

    print "\n" ;
    print " Done submitting scripts for $NTEST_ALL test jobs.\n";
    print " SNANA Test Overview will be piped to: \n" ;
    print "    $MASTER_LOGFILE  \n" ;

    if ( $NEWFITS ) 
    { print " OPTION:  NEW  kcor*.his -> kcor*.fits \n"; }

    if ( $OLDFITS ) 
    { print " OPTION:  OLD  kcor*.his -> kcor*.fits \n"; }

    die   "\n" ;
   
}


# if we get here then process jobs for $SSHNODE

for ( $itest=1; $itest <= $NTEST_ALL ; $itest++ ) {
    #$node  = $TESTLIST_NODE[$IVERSION][$itest] ;
    $inode = $TESTLIST_INODE[$IVERSION][$itest] ;

    if ( $inode == $INODE_SSH )  { 

	&parse_TEST($itest);  # parse instructions 

	&prep_TEST($itest);   # prepare for test: copy files, etc ...

	&wait_for_TEST($itest);

	&run_TEST($itest);    # run the test

    }  # node
}  # itest  

# compare OLD vs. NEW for each TEST
&analyze_RESULTS();


# ===========================
#    END OF MAIN
# ===========================


# ===========================
sub set_defaults {

    $OPT_ABORT  = 0 ; # parsing option
    $OPT_WARN   = 1 ; # 1=> leave warning message
    $OPT_QUIET  = 2 ; # 1=> quiet

    $scriptName    = "$0" ;  # use current script for OLD and NEW versions
    $CDATE         = `date +%Y%m%d_%H%M` ;
    $CDATE         =~ s/\s+$// ;   # trim trailing whitespace

    $SNDATA_ROOT   = $ENV{'SNDATA_ROOT'};
    $SNANA_DIR     = $ENV{'SNANA_DIR'};
    $TESTDIR       = $ENV{'SNANA_TESTS'} ;
#    $TESTDIR       = "/data/des41.b/data/SNDATA_ROOT/SNANA_TESTS" ;

    $TASKDIR       = "${TESTDIR}/tasks" ;
    $TESTLOGDIR    = "${TESTDIR}/logs/${CDATE}" ;
    $TESTINPUTDIR  = "${TESTDIR}/inputs" ;

    @VERS_SUFFIX   = ( "OLD" ,  "NEW" ) ;

    $PROC_VERSION = " " ;
#    $SSHNODE      = " " ;
    $INODE_SSH    = -1 ;

    $SIMGEN_PREFIX       = "SIMGEN";
    $SNTABLE_DUMP_PREFIX = "SNTABLE_DUMP";

    $INDEX_OLD = 0 ;
    $INDEX_NEW = 1 ;
    @OLDNEW_KEY = ( "OLD" , "NEW" );

    $IVERSION     = -1 ;

    $q        = "\'" ;
    $qq       = '"' ;

    $SSHNODE_LOGFILE ;

    $KILLJOBS_FLAG = 0;

    $OLDFITS = 0 ;  # swap .his for .fits in OLD
    $NEWFITS = 0 ;  # swap .his for .fits in NEW

    $OPT_LCPLOT = 0 ;

} # end of set_defaults


# ===========================
sub parse_args {

    my ($NARG, $i);

    $NARG = scalar(@ARGV) ;

    for ( $i = 0; $i < $NARG ; $i++ ) {

#	if ( $ARGV[$i] eq "HELP" ) {  print_help();  }
#	if ( $ARGV[$i] eq "help" ) {  print_help();  }

	if ( $ARGV[$i] eq "--TESTDIR" ) {
	    $TESTDIR    = "$ARGV[$i+1]" ;
	    $TESTLOGDIR = "${TESTDIR}/logs/${CDATE}" ;
	}

	if ( $ARGV[$i] eq "--LOGDIR" ) {
	    $TESTLOGDIR = "$ARGV[$i+1]" ;
	}

#	if ( $ARGV[$i] eq "--NODE" ) { $SSHNODE = "$ARGV[$i+1]" ;  }

	if ( $ARGV[$i] eq "--INODE" ) { 
	    $INODE_SSH = "$ARGV[$i+1]" ;  
	}

	if ( $ARGV[$i] eq "KILL" ) {
	    $KILLJOBS_FLAG = 1;
	}

	if ( $ARGV[$i] eq "OLDFITS" ) {
	    $OLDFITS = 1;
	}

	if ( $ARGV[$i] eq "NEWFITS" ) {
	    $NEWFITS = 1;
	}

	if ( "$ARGV[$i]" eq "OLD" ) {
	    $PROC_VERSION = "OLD" ;
	    $IVERSION     = $INDEX_OLD ;
	}

	if ( "$ARGV[$i]" eq "NEW" ) {
	    $PROC_VERSION = "NEW" ;
	    $IVERSION     = $INDEX_NEW ;
	}

	if ( "$ARGV[$i]" eq "LCPLOT" ) {
	    $OPT_LCPLOT = 1;
	}
    }  # i

} # end of parse_args



# ========================
sub init_stuff {


    $MASTER_LOGFILE = "$TESTLOGDIR/MASTER.LOG" ;
    
    # construct full name of LISTFILE and make sure it's there !

    $LISTFILE = "${TESTDIR}/SNANA_CodeTests_FNAL.LIST" ;
    if ( !( -e $TESTDIR ) ) {
	print " ERROR: Cannot find LISTFILE: $LISTFILE \n" ;
	die   "  ***** ABORT *****  \n " ;
    }


} # end of init_stuff


# ========================
sub parse_LISTFILE {

    my ($key, $val, @tmp);
    my ($iline, $NLINE, @inFileContents, $tmpLine, @wdlist ) ;

    $key = "OLD_SNANA_SETUP:" ;
    @OLD_SNANA_SETUP = sntools::parse_line($LISTFILE, 99, $key, $OPT_ABORT ) ; 

    $key = "NEW_SNANA_SETUP:" ;
    @NEW_SNANA_SETUP = sntools::parse_line($LISTFILE, 99, $key, $OPT_ABORT ) ; 

    $key = "OLD_SNANA_NODES:" ;
    @tmp = sntools::parse_line($LISTFILE, 99, $key, $OPT_ABORT ) ; 
    @OLD_SNANA_NODELIST    = split(/\s+/,$tmp[0]) ;

    $key = "NEW_SNANA_NODES:" ;
    @tmp = sntools::parse_line($LISTFILE, 99, $key, $OPT_ABORT ) ; 
    @NEW_SNANA_NODELIST    = split(/\s+/,$tmp[0]) ;


    # define a few more things that depend on the above inputs

    $NNODE_OLD = scalar(@OLD_SNANA_NODELIST);
    $NNODE_NEW = scalar(@NEW_SNANA_NODELIST);
    $NNODE[$INDEX_OLD] = $NNODE_OLD ;
    $NNODE[$INDEX_NEW] = $NNODE_NEW ;

    $SNANA_SETUP[$INDEX_OLD] = "@OLD_SNANA_SETUP" ;
    $SNANA_SETUP[$INDEX_NEW] = "@NEW_SNANA_SETUP" ;

    for ( $inode=0; $inode < $NNODE_OLD ; $inode++ ) 
    { $NODELIST[$INDEX_OLD][$inode] = $OLD_SNANA_NODELIST[$inode]; }

    for ( $inode=0; $inode < $NNODE_NEW ; $inode++ ) 
    { $NODELIST[$INDEX_NEW][$inode] = $NEW_SNANA_NODELIST[$inode]; }

    # now read TEST: keys sequentially to look for END: key.

    @inFileContents = qx(cat ${LISTFILE}) ;
    $NLINE   = scalar(@inFileContents);
    $NTEST_ALL   = 0 ;

    $iline = 0;
    while ( $iline < $NLINE ) {
        $tmpLine = "$inFileContents[$iline]" ;
        @wdlist   = split(/\s+/,$tmpLine) ;
        $key = $wdlist[0] ;
        $val = $wdlist[1] ;

        if ( $key eq "END:" ) { return ; }

        if ( $key eq "TEST:" ) {
            $NTEST_ALL++ ;  
            $TESTLIST_NAME[$NTEST_ALL] = "$val" ;
        }
        $iline++ ;
    }

} # end of parse_LISTFILE


# ===================
sub sort_TESTLIST {

    # sort TESTLIST_NAME so that the SIMGEN tests are first.

    my ($TEST, $prefix6, $prefix7, $i, $j);
    my (@TESTLIST_NOTSIM, @TESTLIST_SNTABLE );
    my @TMPLIST = @TESTLIST_NAME ;


    $NTEST_SIMGEN = 0;
    $NTEST_NOTSIM = 0;
    $NTEST_SNTABLE = 0 ;

    for ( $i=1; $i <= $NTEST_ALL ; $i++ ) {
	$TEST    = $TMPLIST[$i] ;
	$prefix6 = substr($TEST,0,6);  # for SIMGEN
	$prefix7 = substr($TEST,0,7);  # for SNTABLE
	$SNTABLE_FLAG[$i] = 0 ;

	if ( $prefix6 eq $SIMGEN_PREFIX ) {
	    $NTEST_SIMGEN++ ;
	    $TESTLIST_NAME[$NTEST_SIMGEN] = $TEST ;
	}
	elsif ( $prefix7 eq $SNTABLE_DUMP_PREFIX ) {
	    $NTEST_SNTABLE++ ;
	    $TESTLIST_SNTABLE[$NTEST_SIMGEN] = $TEST ;
	    $SNTABLE_FLAG[$i] = 1 ;
	}
	else {
	    $NTEST_NOTSIM++ ;
	    $TESTLIST_NOTSIM[$NTEST_NOTSIM] = $TEST ;
	}
    }

    # now tack on the non-SIMGEN tests at the end of the official list

    for ( $i=1; $i <= $NTEST_NOTSIM ; $i++ ) {
	$j = $i + $NTEST_SIMGEN ;
	$TESTLIST_NAME[$j] = $TESTLIST_NOTSIM[$i] ;
    }

    for ( $i=1; $i <= $NTEST_SNTABLE ; $i++ ) {
	$j = $i + $NTEST_SIMGEN + $NTEST_NOTSIM ;
	$TESTLIST_NAME[$j] = $TESTLIST_SNTABLE[$i] ;
    }

    # July 2016: and tack on the SNTABLE tests (since LCFIT must run first)
    

} # end of sort_TESTLIST

# ======================
sub assign_NODE {

    # Assign node to each TEST in round-robin order.
    # Note that node depends on OLD/NEW snana version.
    # Fill 2-dim array TESTLIST_NODE[$iver][$itest].
    # Jun 2014: also fill TESTLIST_INODE

    my ($iver, $cver, $inode, $LDMP, $c3, $TEST, $i) ;

    if ( $INODE_SSH >= 0 ) {
#	print " xxx $PROC_VERSION NODELIST = @NODELIST \n";
	my $SSHNODE = "$NODELIST[$IVERSION][$INODE_SSH]" ;
	$SSHNODE_LOGFILE = "$TESTLOGDIR/${PROC_VERSION}_${SSHNODE}.LOG";
	open  PTR_NODELOG , "> $SSHNODE_LOGFILE" ; 
	print PTR_NODELOG " $PROC_VERSION SNANA_DIR = $SNANA_DIR \n\n";
    }

    $LDMP = 0 ;
    if  ( $PROC_VERSION eq " " ) { $LDMP = 1 ; }

    for ( $iver = 0 ; $iver <= 1 ; $iver++ ) {  # loop over OLD/NEW

	$inode = 0;

	for ( $i=1; $i <= $NTEST_ALL ; $i++ ) {

	    $c3   = sprintf("%3.3d", $i);

	    $node = $NODELIST[$iver][$inode];
	    $inode++ ;
	    if ( $inode >= $NNODE[$iver] ) { $inode = 0 ; }

	    $TESTLIST_NODE[$iver][$i]  = $node;
	    $TESTLIST_INODE[$iver][$i] = $inode;
	    $TESTLIST_IDCHAR[$i]       = $c3 ;  # formatted ID num

	    if ( $LDMP )  {
		$TEST = $TESTLIST_NAME[$i] ;
		$cver = $OLDNEW_KEY[$iver];
		print "\t SORTED TEST[$c3-$cver] = $TEST -> $node \n";
	    }
	}
    }

} # end of assign_NODE



# ========================
sub make_LOGDIR {

    if ( -d $TESTLOGDIR )  { qx(rm -r $TESTLOGDIR); }
    qx(mkdir  $TESTLOGDIR);

    # open logFile

    open  PTR_LOGFILE , ">> $MASTER_LOGFILE" ;  
 
    print PTR_LOGFILE " Overview log file for SNANA tests. \n";
    print PTR_LOGFILE " LISTFILE:  $LISTFILE \n";

    close PTR_LOGFILE ;

} # end of make_LOGDIR

# =======================
sub runScript {

    # run this same script via 'ssh -$node' with arguments 
    # OLD/NEW and --NODE $node.

    my ($iver,$inode) = @_;

    my ($VERS, $node, $cmd, $argList, $setup, $cdd, $OPTIONS );

    $VERS    = $OLDNEW_KEY[$iver] ;
    $node    = $NODELIST[$iver][$inode] ;
    $cdd     = "cd $TESTLOGDIR" ;

    $OPTIONS = "" ;
    if ( $NEWFITS ) { $OPTIONS = "NEWFITS" ; }
    if ( $OLDFITS ) { $OPTIONS = "OLDFITS" ; }

#    $argList = "$VERS --NODE $node --LOGDIR $TESTLOGDIR $OPTIONS" ;   
    $argList = "$VERS --INODE $inode --LOGDIR $TESTLOGDIR $OPTIONS" ;   
    $setup   = "$SNANA_SETUP[$iver]" ;
    $cmd     = "ssh -x $node ${qq}$setup ; $cdd ; $scriptName $argList${qq}" ;

    print "  Launch $VERS SNANA test jobs on node = $node \n" ;

    system("$cmd &"); # DDDDDDDDD

} # end of submit_jobs


# =======================
sub parse_TEST {

    # parse instructions from test file $itest.
    # Jan 1 2013: allow multiple TESTINPUT: lines.

    my ($itest) = @_;

    my ($key, @tmp, $tmpArg );
    my $TEST     = $TESTLIST_NAME[$itest] ;
    my $TESTFILE = "$TASKDIR/$TEST" ;
# xxx    my $TESTFILE = "$TESTDIR/$TEST" ; 

    $key = "TESTJOB:" ;
    @tmp = sntools::parse_line($TESTFILE, 1, $key, $OPT_ABORT ) ; 
    $TESTJOB = $tmp[0] ;

    $key = "TESTINPUT:" ;
    @TESTINPUT = () ;
    @tmp = sntools::parse_line($TESTFILE, 99, $key, $OPT_ABORT ) ; 
    foreach $tmpArg (@tmp) {
	@TESTINPUT = ( @TESTINPUT , split(/\s+/,$tmpArg) )  ;
    }


    # read option arguments if argument is not just an input filename
    $TESTJOB_ARGS = "" ;
    $key = "TESTJOB_ARGS:" ;
    @tmp = sntools::parse_line($TESTFILE, 99, $key, $OPT_QUIET ) ; 
    if( scalar(@tmp) > 0 ) { $TESTJOB_ARGS = "$tmp[0]" ; }

    $key = "TESTNAME:" ;
    @tmp = sntools::parse_line($TESTFILE, 1, $key, $OPT_ABORT ) ; 
    $TESTNAME = $tmp[0] ;


    $key = "TESTRESULT:" ;
    @tmp = sntools::parse_line($TESTFILE, 99, $key, $OPT_ABORT ) ; 
    $TESTRESULT = "$tmp[0]" ;

    my $key = "WORDNUM:" ;
    my @tmp = sntools::parse_line($TESTFILE, 1, $key, $OPT_ABORT ) ; 
    $WORDNUM = $tmp[0] ;

} # end of parse_TEST


# =======================
sub prep_TEST {

    # copy & prepare files.
    # Fill globals:  TEST_INFILE  TEST_LOGFILE

    my ($itest) = @_;

    my ($TEST, $ID, $LDMP, $VERS, $cmd, $sedcmd );
    my (@wdnum, @wdlist, $wd, $wd2) ;
    my ($inFile, $cpFile,  $PREFIX, $comment, $ifile, $F1tmp, $F2tmp);

    # ------------- BEGIN --------------

    $LDMP = ( $itest < 0 ) ;
    $TEST     = $TESTLIST_NAME[$itest] ;    # name of test
    $ID	      = $TESTLIST_IDCHAR[$itest] ;  # formatted test-ID number
    $VERS     = $OLDNEW_KEY[$IVERSION] ;    # OLD or NEW
    $PREFIX   = "${VERS}${ID}" ;    

    # on first test, update SNANA_DIR in the master log
    if ( $itest == 1 ) {
	$comment = " $VERS  SNANA_DIR = $SNANA_DIR";
	qx(echo "$comment" >> $MASTER_LOGFILE);
    }

    # copy input file(s) to logDir, except for SNTABLE_DUMP jobs
    if ( $SNTABLE_FLAG[$itest] ) { goto AFTER_COPY ; }

    $ifile = 0;
    foreach $inFile ( @TESTINPUT ) {

	if ( $ifile == 0 ) 
	{  $cpFile = "${PREFIX}_${inFile}_TEMP" ; }
	else
	{  $cpFile = "$inFile" ; }

	if ( !(-e "$TESTLOGDIR/$cpFile" ) )  { 
	    $cmd = "cp $TESTINPUTDIR/$inFile $cpFile" ;
	    qx(cd $TESTLOGDIR ; $cmd); 
	}

	$ifile++ ;
    }

    # Apr 2019: if XXX appears in TESTINPUT, then it's an output
    # of another test job, so XXX -> OLD/NEW; 

    # Only the first input file get _OLD[NEW] appended and XXX -> OLD/NEW    
    $inFile       = "$TESTINPUT[0]" ;

    # Apr 2019: HBOOK/ROOT input is output from previous job
    # PROBLEM: need to wait for HBOOK file to be created ???
#    if ( index($inFile,".ROOT")>0 || index($inFile,".HBOOK")>0 )  { 
#	$TEST_INFILE  = "${VERS}_${inFile}" ; 
#    }
#    else {
    $TEST_INFILE  = "${PREFIX}_${inFile}" ;
    $F1tmp        = "${TEST_INFILE}_TEMP" ;
    $F2tmp        = "${TEST_INFILE}" ;
    $sedcmd = "sed -e ${qq}s/XXX/${VERS}/g${qq}  $F1tmp > $F2tmp" ;    
    qx(cd $TESTLOGDIR ; $sedcmd ; rm $F1tmp );
#    }

  AFTER_COPY:

    # construct name of log-file and substitute logFile in TESTRESULT
    $TEST_LOGFILE = "${PREFIX}_${TEST}.LOG" ;
    $TESTRESULT   =~ s/logFile/$TEST_LOGFILE/;  

    # construct name of DONE  file
    $TEST_DONEFILE = "${PREFIX}_${TEST}.DONE" ;

    # extract word numbers
    $WORDNUM   =~ s/-/  / ;  
    @wdnum = split(/\s+/,$WORDNUM) ;
    $TEST_WORDNUM_FIRST = $wdnum[0] ;
    $TEST_WORDNUM_LAST  = $wdnum[1] ;

    # if TESTINPUT is specified in TESTJOB_ARGS,
    # replace it with TEST_INFILE
    $TEST_ARGLIST = "";
    if ( $TESTJOB_ARGS ne "" ) {
	@wdlist = split(/\s+/,$TESTJOB_ARGS) ;
	
	foreach $wd ( @wdlist ) {	   
	    $wd2 = $wd ;
	    if ( $wd eq "TESTINPUT" ) { $wd2 = $TEST_INFILE ; }		
	    $TEST_ARGLIST = "$TEST_ARGLIST $wd2 " ;
	}
    }

    if ( $LDMP ) {
	print "  Process TEST[$ID] = '${TEST}' \n";
	print "\t TESTJOB:      $TESTJOB \n";
	print "\t TESTJOB_ARGS: $TESTJOB_ARGS \n";
	print "\t TESTINPUT:    @TESTINPUT \n";
	print "\t TESTNAME:     $TESTNAME \n";
	print "\t TESTRESULT:   $TESTRESULT \n";
	print "\t WORDNUM:      $WORDNUM \n";
    }

} # prep_TEST


# =======================
sub run_TEST {

    # run TESTLIST_NAME[$itest]
    # Sep 21, 2013: if test-result value is blank, write BLANK-$VERS
    #               so that two blank results are not flagged as a match.

    my ($itest) = @_;

    my ($STRING_IDTEST);
    my ($TEST, $TEST_RESULT_STRING, $TEST_RESULT_VAL);
    my ($w1, $w2, $Nwd);
    my ($VERS, $ID, $NTMP, $cdd, $cmd_job, $ctmp1, $ctmp2);
    my ($cmd_DONE, $cmd_ABORT, $cmd_RESULT, @tmp, @tmp2);
    # -------------- BEGIN ------------

    $TEST     = $TESTLIST_NAME[$itest] ;
    $ID       = $TESTLIST_IDCHAR[$itest] ;  # formatted ID number
    $VERS     = $OLDNEW_KEY[$IVERSION] ;    # OLD or NEW

    $STRING_IDTEST = "${ID}-${TEST}:" ;

    $cdd      = "cd $TESTLOGDIR" ;

    if ( $TESTJOB_ARGS eq "" ) {
	$cmd_job  = "$TESTJOB $TEST_INFILE >& $TEST_LOGFILE" ;
    }
    else {
	$cmd_job  = "$TESTJOB $TEST_ARGLIST >& $TEST_LOGFILE" ;
    }
    $cmd_DONE = "touch $TEST_DONEFILE" ;

    print PTR_NODELOG " Process TEST[$ID] = $TEST \n";
    qx($cdd ; $cmd_job ; $cmd_DONE );

    # first check for ABORT; otherwise extract one-line RESULT from
    # log-file and paste into the DONE-file.

    $cmd_ABORT = "grep $qq ABORT $qq $TEST_LOGFILE" ;
    @tmp = qx($cdd ; $cmd_ABORT );
    $NTMP = scalar(@tmp);

    if ( $NTMP > 0 ) { 
	$TEST_RESULT_STRING = "${STRING_IDTEST}   ABORT-$VERS" ; 
    }
    else {
	# use the TESTRESULT command to extract the result for testing;
	# store this one-line result in the DONE file.
	# The 'tail -1' ensures that the last instance is used.
	$cmd_RESULT = "$TESTRESULT | tail -1" ;
	@tmp    = qx($cdd ; $cmd_RESULT );
	$tmp[0] =~ s/^\s+// ;   # remove leading spaces  (5.13.2019)
	$tmp[0] =~ s/\s+$// ;   # trim trailing whitespace
	@tmp2   = split(/\s+/,$tmp[0]) ;
	$w1     = $TEST_WORDNUM_FIRST ;
	$w2     = $TEST_WORDNUM_LAST  ;
	$Nwd    = scalar(@tmp2);

	if ( $Nwd > 0 ) 
	{ $TEST_RESULT_VAL = "@tmp2[$w1 .. $w2]" ; }
	else
	{ $TEST_RESULT_VAL = "BLANK-$VERS" ; }

	$ctmp1 = "${STRING_IDTEST}" ;
	$ctmp2 = "$TESTNAME = $TEST_RESULT_VAL" ;
	$TEST_RESULT_STRING = sprintf("%-40s  %s  ", $ctmp1 , $ctmp2 ) ;
    }

    # update RESULT in DONE file.
    qx($cdd ; echo "$TEST_RESULT_STRING" >> $TEST_DONEFILE );

} # end of run_TEST


# =======================
sub wait_for_TEST {

    # For SIMGEN test there is no wait.
    # After SIMGEN then ALL SIMGEN tests must be done to continue.

    my ($itest) = @_;

    my ($NDONE, $search, $VERS, $tmp, $cdate, $cdd) ;
    my (@DONELIST) ;

    # -------------- BEGIN ----------

    if ( $NTEST_SIMGEN == 0 ) { return ; }

    # No waiting for SIMGEN test ...
    if ( $itest <= $NTEST_SIMGEN ) { return ; }

    $cdd      = "cd $TESTLOGDIR";
    $VERS     = $OLDNEW_KEY[$IVERSION] ;    # OLD or NEW
    $search   = "${VERS}*${SIMGEN_PREFIX}*DONE" ;

  CHECK_DONE:

    # do ls on DONE files, but suppress STDOUT
    @DONELIST  = `$cdd ; ls ${search} 2>/dev/null` ;
    $NDONE     = scalar(@DONELIST);

    if ( $NDONE < $NTEST_SIMGEN ) {
	$cdate  = `date` ;  $cdate  =~ s/\s+$// ;   # trim  whitespace
	print PTR_NODELOG 
	    "   Wait for $NTEST_SIMGEN ${search} files " .
	    "(found $NDONE at $cdate)\n" ;
	sleep(10);
	goto CHECK_DONE ;
    }

    # for SNTABLE_DUMP, make sure all other jobs are done
    #   xyz
    if ( $SNTABLE_FLAG[$itest] ) { 
	my $NDONE_EXPECT = $NTEST_ALL - $NTEST_SNTABLE ;
	if ( $NDONE < $NDONE_EXPECT ) {
	    $cdate  = `date` ;  $cdate  =~ s/\s+$// ;   # trim  whitespace
	    print PTR_NODELOG 
		"   Wait for $NDONE_EXPECT ${VERS}*DONE files " .
		"(found $NDONE at $cdate)\n";
	    sleep(5);
	    goto CHECK_DONE ;
	}
    }

} # end of wait_for_TEST


# ============================
sub analyze_RESULTS {

    # if all jobs have finished then compare 
    # OLD vs. NEW and report discrepancies.
    # Otherwise just return;

    my ($cdd, $cmd, $iver, $VERS, $iline) ;
    my (@DONELIST, $OLD, $NEW, $TEST, $CTEST) ;
    my ($NDONE, $NDONE_EXPECT, @RESULTS_FILE );
    my ($NTEST_GOOD, $NTEST_BAD, $NLINE_RESULTS);
    my (@Contents_OLD, @tmp_OLD, $TEST_OLD);
    my (@Contents_NEW, @tmp_NEW, $TEST_NEW);
    
    # check number of DONE files .. should be $NTEST_ALL 
    # if all have finished

    $NDONE_EXPECT = 2 * $NTEST_ALL ;

    $cdd      = "cd $TESTLOGDIR";
    $cmd      = "ls *.DONE 2>/dev/null" ;
    @DONELIST = qx($cdd ; $cmd );
    $NDONE    = scalar(@DONELIST);

    if ( $NDONE < $NDONE_EXPECT ) { return ; }    

    # --------------------------------------------
    # all jobs have finished ... check  results.
    print " Analyze SNANA test results. \n";

    print PTR_NODELOG "\n Found all DONE files; analyze results. \n";

    open  PTR_LOGFILE , ">> $MASTER_LOGFILE" ;  
    print PTR_LOGFILE "\n";
    print PTR_LOGFILE " ----------------------------------- \n";
    print PTR_LOGFILE " Compare OLD-vs-NEW test results. \n\n";

    # catenate the results (OLD/NEW) into one RESULTS
    for ( $iver=0; $iver <=1; $iver++ ) {
	$VERS     = $OLDNEW_KEY[$iver] ;
	$RESULTS_FILE[$iver] = "RESULTS_${VERS}.LOG" ;
	$cmd = "cat $VERS*.DONE > $RESULTS_FILE[$iver]" ;
	qx($cdd ; $cmd);
    }

    @Contents_OLD  = qx($cdd ; cat $RESULTS_FILE[$INDEX_OLD]) ;
    @Contents_NEW  = qx($cdd ; cat $RESULTS_FILE[$INDEX_NEW]) ;
    $NLINE_RESULTS = scalar(@Contents_OLD);

    $iline = 0;
    $NTEST_GOOD = 0;
    $NTEST_BAD  = 0;

    while ( $iline < $NLINE_RESULTS ) {
        $OLD = "$Contents_OLD[$iline]" ;
        $NEW = "$Contents_NEW[$iline]" ;
	
	@tmp_OLD = split(/\s+/,$OLD) ;    $TEST_OLD = $tmp_OLD[0];
	@tmp_NEW = split(/\s+/,$NEW) ;    $TEST_NEW = $tmp_NEW[0];
	$TEST = $TEST_NEW ;
	$CTEST = sprintf("%-40s", $TEST );

	if ( $TEST_OLD ne $TEST_NEW ) {
	    print "\nERROR: line mis-match: \n";
	    print " TEST_OLD = '${TEST_OLD}' at line=$iline \n";
	    print " TEST_NEW = '${TEST_NEW}' at line=$iline \n";
	    die   " ***** ABORT ***** \n";
	}


	if ( "@tmp_OLD" eq "@tmp_NEW" ) {
	    print PTR_LOGFILE "$CTEST  Perfect match (OLD = NEW) \n";
	    $NTEST_GOOD++ ;
	}
	else {
	    print PTR_LOGFILE "$CTEST  mis-match (OLD != NEW) \n";
	    print PTR_LOGFILE "    ==> OLD: @tmp_OLD[1 .. 9] \n";
	    print PTR_LOGFILE "    ==> NEW: @tmp_NEW[1 .. 9] \n";
	    $NTEST_BAD++ ;
	}

	$iline++ ;	
    }

    print PTR_LOGFILE "\nTest Summary: \n";
    print PTR_LOGFILE "  $NTEST_GOOD tests have perfect match. \n";
    print PTR_LOGFILE "  $NTEST_BAD tests have mis-match. \n";

    close PTR_LOGFILE ;

    system("cat $MASTER_LOGFILE") ;

} # end  of  analyze_RESULTS 



# ==================================
sub LCPLOT_DRIVER {

    # Created May 2014
    # driver to interactively run the LCPLOT*.nml jobs 
    # in $TESTINPUTDIR, and then run mkfitplots.pl to 
    # create pdf files of lightcurve plots. User must
    # them visually scan the pdf files.
    # Each nml file must contain the key
    #  JOBNAME:  <jobname>
    # so that we know which program to run;
    # a missing key results in an ABORT.

    my ($cdinp, $NJOB, $ijob, @NMLFILE_LIST, $NMLFILE, $nmlFile) ;
    my (@JOBNAME_LIST, $jobName, $logFile, $CMD);
    my ($key, @tmp, @wdlist ) ;

    $cdinp = "cd $TESTINPUTDIR" ;

    # first get list of nml files and program names.
    
    @NMLFILE_LIST = qx($cdinp ; ls LCPLOT*.nml) ;
    $NJOB = 0;

    print " Preview LCPLOT Jobs: \n";
    foreach $nmlFile (@NMLFILE_LIST) {
	$nmlFile  =~ s/\s+$// ;   # trim trailing whitespace
	$NMLFILE  = "$TESTINPUTDIR/$nmlFile";
	$key      = "JOBNAME:" ;
	@tmp      = sntools::parse_line($NMLFILE, 1, $key, $OPT_ABORT) ;
	$jobName  = $tmp[0] ;
	$JOBNAME_LIST[$NJOB] = $jobName ;
	$NJOB++ ;
	my $line = sprintf("%-12.12s  %s", $jobName, $nmlFile);
	print "\t $line \n" ;
    }
    print " Finished preview. \n";
    
    # check Makefile for HBOOK and ROOT.
    # look for 'USE_HBOOK =  1'  and/or 'USE_ROOT = 1'.
    my ($OUTFILE_HBOOK, $OUTFILE_ROOT, $ARG_OUTFILE );
    my $USE_HBOOK = 0 ;
    my $USE_ROOT  = 0 ;
    my $Makefile  = "$SNANA_DIR/src/Makefile" ;

    @tmp  = sntools::parse_line($Makefile, 2, "USE_HBOOK", $OPT_QUIET) ;  
    @wdlist = split(/\s+/,$tmp[0]) ;
    if ( $wdlist[1] == 1 ) { $USE_HBOOK = 1; }

    @tmp  = sntools::parse_line($Makefile, 2, "USE_ROOT", $OPT_QUIET) ;  
    @wdlist = split(/\s+/,$tmp[0]) ;
    if ( $wdlist[1] == 1 ) { $USE_ROOT = 1; }

    print "\t USE(HBOOK,ROOT) = $USE_HBOOK, $USE_ROOT \n";

    my $LCPLOTDIR = "${TESTDIR}/logs/${CDATE}_LCPLOT" ;
    my $cdout     = "cd $LCPLOTDIR" ;
    qx(mkdir $LCPLOTDIR);
    print "\n Created output directory:\n $LCPLOTDIR \n\n";
    $| = 1 ;

    print " Run LCPLOT jobs: \n";
    for($ijob=0; $ijob < $NJOB; $ijob++ ) {
	$jobName = $JOBNAME_LIST[$ijob];
	$nmlFile = $NMLFILE_LIST[$ijob];
	qx(cp $TESTINPUTDIR/$nmlFile $LCPLOTDIR);

	# construct name of log file
	my $jdot   = index($nmlFile,".");
	my $prefix = substr($nmlFile,0,$jdot) ;
	$logFile   = "$prefix" . ".LOG" ;
	
	$OUTFILE_HBOOK = "${prefix}.HBOOK" ;
	$OUTFILE_ROOT  = "${prefix}.ROOT" ;

	$ARG_OUTFILE = "" ;
	if ( $USE_HBOOK ) 
	{ $ARG_OUTFILE = "$ARG_OUTFILE  HFILE_OUT $OUTFILE_HBOOK" ; }

	if ( $USE_ROOT ) 
	{ $ARG_OUTFILE = "$ARG_OUTFILE  ROOTFILE_OUT $OUTFILE_ROOT" ; }

	$CMD = "$jobName $nmlFile  $ARG_OUTFILE >& $logFile" ;
	print "     $jobName $nmlFile ... \n";
	qx($cdout; $CMD) ;

	# make pdf file(s)
	if($USE_HBOOK) {
	    $CMD = "mkfitplots.pl -h $OUTFILE_HBOOK";
	    print "\t $CMD ... \n";
	    qx($cdout ; $CMD); 
	}
	if($USE_ROOT ) { 
	    $CMD = "mkfitplots.pl -r $OUTFILE_ROOT";
	    print "\t $CMD ... \n";
	    qx($cdout ; $CMD); 
	}

	$| = 1 ;
    }

    my @PDFLIST = qx($cdout ; ls *.pdf );
    print "\n Created the following pdf files: \n @PDFLIST \n";

    # make tarball
    my $tarFile = "LCPLOT_pdfFiles.tar" ;
    qx($cdout ; tar -cf $tarFile *.pdf ; gzip $tarFile );
    print " LCPLOT-pdf tarball: \n $LCPLOTDIR/${tarFile}.gz \n";

    die "\n Finished all LCPLOT tests. \n";

} # end of LCPLOT_DRIVER
