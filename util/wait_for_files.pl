#!/usr/bin/perl
#
#
# Created Sep 2011 by R.Kessler
#
# wait for $NFILE files to exist of the form passed as argument.
# When all $NFILE files exist, then create the STAMP-file.
#
# Syntax
#
#  wait_for_files.pl  <NFILE>  <FILESPEC> <STAMPFILE>
#  wait_for_files.pl  <NFILE>  <FILESPEC> NOSTAMP  
#
# Example:
#  wait_for_files.pl  5  '*.DONE'  ALL.DONE
#  ==> when 5 files ending in ".DONE" exist, then ALL.DONE is created
#      and function exits.
#      This function sleeps until all 5 DONE files exist.
#
#  wait_for_files.pl  1  'MYFILE.TXT'  NOSTAMP
#   ==> when MYFILE.TXT exists, exit without creating stamp.
#
#  wait_for_files.pl  -5  '*.BUSY'  NOSTAMP
#  ==> exit when there are fewer than 5 *.BUSY files.
#
#
# HISTORY
#  Apr 6 2015: new NOSTAMP option
#  Jul 17 2017: allow N>NWANT
#  Dec 12 2017: 
#    + new logic to exit when fewer than <NBUSY> files.
#    + sleep time -> 5 seconds (was 10)
#
# ------------

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

sub parse_args ;
sub Nfile_exist ;

my ($NFILE_WANT, $FILESPEC, $STAMPFILE );
my ($LOGICTYPE_DONE, $LOGICTYPE_BUSY);
my $SLEEPTIME = 5 ;  # sleep time between each check

# ======== BEGIN MAIN ==========

&parse_args();

my ($N, @output, $errout, $FOUND) ;

CHECKFILES:
$N = &Nfile_exist();
$FOUND=0;

if ( $LOGICTYPE_DONE &&  ($N >= $NFILE_WANT) ) { $FOUND=1; }
if ( $LOGICTYPE_BUSY &&  ($N  < $NFILE_WANT) ) { $FOUND=1; }

if ( $FOUND ) {

    # of NOSTAMP is specified, then do nothing when all files are found.
    if ( $STAMPFILE ne "NOSTAMP" ) {
	print "wait_for_files.pl found $N files with '$FILESPEC' \n" ;
	@output = qx(touch $STAMPFILE);    
	$errout = "wait_for_files.pl (stamp-file creation)";
	print "   => create stamp-file: $STAMPFILE \n";
	sntools::checkFileExist($STAMPFILE, $errout);
    }
}
else {
    # print "Found only $N of  $NFILE_WANT files => snooze. \n";
    sleep($SLEEPTIME);
    goto CHECKFILES ;
}

exit($N);

# ======= END MAIN ===========

sub parse_args {


    my ($NARG, @MSGERR, $arg0);

    $LOGICTYPE_DONE=1;  # default
    $LOGICTYPE_BUSY=0; 

    $NARG = scalar(@ARGV);
    $arg0 = $ARGV[0];

    if ( $NARG < 3 ) {
	$MSGERR[0] = "3 arguments required (but NARG=$NARG).";
	$MSGERR[1] = "$0 <NFILE> <FILESPEC> <STAMPFILE> ";
	$MSGERR[2] = "  ARGV[0] = $ARGV[0] \n" ;
	$MSGERR[3] = "  ARGV[1] = $ARGV[1] \n" ;
	$MSGERR[3] = "  ARGV[2] = $ARGV[2] \n" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    $NFILE_WANT = $ARGV[0];
    $FILESPEC   = $ARGV[1];
    $STAMPFILE  = $ARGV[2];

    if ( substr($arg0,0,1) eq '-' ) {
	$NFILE_WANT = substr($arg0,1,4);	    
	$LOGICTYPE_DONE=0 ;
	$LOGICTYPE_BUSY=1 ;
    }

}  # end of parse_args


# ================
sub Nfile_exist {

    my (@tmp, $NTMP);

    @tmp = qx(ls $FILESPEC  2>/dev/null );
    $NTMP = scalar(@tmp);
    return $NTMP ;

} # end of Nfile_exist
