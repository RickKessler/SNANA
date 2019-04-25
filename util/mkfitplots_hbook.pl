#!/usr/bin/perl
#
# Created Aug 2013 by R.Kessler
# make plots from hbook file using paw macro snana#fitres
#    (analog of mkfitplots_root.c )
# This script is called from mkfitplots.pl.
#
# Usage:
#   mkfitplots_hbook.pl $INFILE '$PLOTARGS'
#   mkfitplots_hbook.pl $INFILE '$PLOTARGS' PRIVATE  ! use private ~/kumacs
#
# Output is a pdf file name is ${INFILE}.pdf
#
#
# =========================================================

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

my ($INFILE, $PLOTARGS, $USE_PRIVATE_KUMACS);

sub parseArgs ;
sub runPaw ;
sub ps2pdf ;

# ======== START MAIN ======

&parseArgs();
&runPAW();
&ps2pdf();

# ======== END MAIN ========


sub parseArgs {

    $USE_PRIVATE_KUMACS = 0 ;

    # start with required args
    $INFILE   = $ARGV[0] ;  # 1st arg is name of hbook file
    $PLOTARGS = $ARGV[1] ;  # 2nd arg is list of paw-macro options

    # check for options
    my ($arg) ;
    foreach $arg ( @ARGV ) {
	if ( $arg eq "PRIVATE" ) {  $USE_PRIVATE_KUMACS = 1 ; }
    }

    print "\n";
    print " HBOOK filename: $INFILE \n";
    print " HBOOK/PAW arguments: '$PLOTARGS' \n";

}   # end of parseArgs

# ============
sub runPAW {

    my ($kumFile, $logFile, $suffix, @CMD, $NCMD, $cmd, $harg, $argtype ) ;
    my ($CMD);

    $harg = "hisfile=${INFILE}" ;
    $CMD  = "exec snana#fitres_mkplots ${harg} ${PLOTARGS}";

    $suffix = qx(date +%Y%m%d);
    $suffix =~ s/\s+$// ;   # trim trailing whitespace   

    $kumFile = "${INFILE}_${suffix}.kumac" ;
    $logFile = "${INFILE}_${suffix}.log" ;

    if ( -d $kumFile ) {  qx(rm $kumFile) ; }
    if ( -d $logFile ) {  qx(rm $logFile) ; }
    
    open KFILE , "> ${kumFile}" ;
    
    if ( $USE_PRIVATE_KUMACS == 0 ) 
    { sntools::write_loginMacroLines(\*KFILE,1); }
    else
    { print "   Use private macros in ~/kumacs \n"; }
    

    print KFILE "$CMD \n";
    
    close KFILE ;

    my @bla = qx(which paw);
    print "  paw command: @bla \n";

    print "\t extracting plots ... please be patient ... \n";

    qx(paw -b ${kumFile} >& $logFile ) ;
    qx(rm  -f ${kumFile} ${logFile} );
}


# ================================
sub ps2pdf {

    my ($tmp_psFile, $tmp_pdfFile ) ;

    # pdf filename is based on $INFILE

    $tmp_psFile  = "${INFILE}.ps";
    $tmp_pdfFile = "${INFILE}.pdf";


    # check for ps file ... won't be there if eps-per-SN option is used.
    if ( -e $tmp_psFile ) {
	print "\t convert ps -> pdf ... \n";
	qx(ps2pdf $tmp_psFile; rm $tmp_psFile );
    }

#    print " See plots in : $tmp_pdfFile \n" ;
}

# =========== END ============
