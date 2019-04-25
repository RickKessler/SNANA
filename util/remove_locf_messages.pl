#!/usr/bin/perl
#
# Remove annoying locf messages from [log] file.
# These messages come from 64-bit CERNLIB calls .
#
# Usage:
#   remove_locf_messages.pl  <file>
#      or
#   remove_locf_messages.pl  <file>  QUIET
#
# Note that <file> is modified, so if you want to
# keep the original then copy it first .
#
# The sed command is
#  sed -e '/locf/d' -e '/crash/d' -e '/!!!!!!/d' $INFILE 
#
# --------------------------------

use FindBin qw($Bin);
use lib "$Bin";
use sntools ;
use strict ;

my ($INFILE, $TMPFILE, $CMD, @MSGERR );

my $OPT_MSG = 1;

# ======== MAIN ============

$INFILE  = $ARGV[0] ;

if ( scalar(@ARGV) > 1 ) {
    if ( $ARGV[1] eq "QUIET" ) { $OPT_MSG = 0 ; }
    if ( $ARGV[1] eq "quiet" ) { $OPT_MSG = 0 ; }
}

if ( length($INFILE) == 0 ) {
    $MSGERR[0] = "Must give file name as argument." ;
    sntools::FATAL_ERROR(@MSGERR);
}

$TMPFILE = "${INFILE}.TEMP" ;

$CMD = "sed " ;
$CMD = "$CMD -e '/locf/d'" ;
$CMD = "$CMD  -e '/crash/d'" ; 
$CMD = "$CMD  -e '/!!!!!!/d'" ;

# print " CMD = $CMD \n";

qx($CMD $INFILE > $TMPFILE ; mv $TMPFILE $INFILE );

if ( $OPT_MSG ) {
    print "locf messages removed from $INFILE \n" ;
}

# ============== END ============
