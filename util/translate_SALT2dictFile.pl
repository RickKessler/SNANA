#!/usr/bin/perl
#
# Created April 2013 by R.Kessler
# Translate the original salt2 dict-file output
# into snana "fitres" format.  The key "SN" is required.
#
# USAGE:
#    translate_SALT2dictFile.pl  $INFILE
#
# and the output SNANA-formatted file is SNANA_${INFILE}.FITRES
#
# Jun 21, 2013: fix a few problems reading header; allow/ignore @XXX
#
# ---------------

use FindBin qw($Bin);
use lib "$Bin";
use sntools;
use strict ;

sub get_INFILE ;
sub get_VARLIST ;
sub open_OUTFILE ;
sub update_OUTFILE ;

my $SNKEY = "SN" ;  # required key in dict file

my ($INFILE, $OUTFILE, @MSGERR, @CONTENTS );
my ($NVAR, @VARNAMES_DICT, @VARNAMES_SNANA, $IVAR_SNKEY );

# =========== BEGIN MAIN ==========

&get_INFILE() ;

my @CONTENTS = `cat $INFILE` ;

&get_VARLIST();

&open_OUTFILE();

my $line ;
foreach $line (@CONTENTS) { &update_OUTFILE($line) ; }

close PTR_OUTFILE ;
print " Done. \n";

# ========== END MAIN ============

sub get_INFILE {
    $INFILE = $ARGV[0] ;
    if ( !(-e $INFILE) ) {
	$MSGERR[0] = "Cannot find INFILE = '$INFILE'" ;
	sntools::FATAL_ERROR(@MSGERR);    
    }

    print " Input dict file: $INFILE \n";

} # end of get_INFILE


sub get_VARLIST {

    # read header to get list of variables.
    # load $NVAR and @VARNAMES_DICT
    # Each variable name has a semicolon ":" after it.
    # FOr header the first word must be either '#',
    # or @XXXX'

    my ($line, @wdlist, $varName, $VARNAME, $NLINE_READ, $wd0, $wd1, $wd2 );

    $IVAR_SNKEY = -9 ;
    $NVAR = 0 ;
    $NLINE_READ = 0 ;

    foreach $line (@CONTENTS) {
	
	$NLINE_READ++ ;
	@wdlist = split(/\s+/,$line);
	$wd0 =  $wdlist[0] ;
	$wd1 =  $wdlist[1] ;
	$wd2 =  $wdlist[2] ;
	
	if ( substr($wd0,0,1) eq "@" ) { next ; }

	if ( $wd0 ne '#' ) { goto DONE_HEADER ; }
	if ( $wd2 ne ':' ) { next   ; }

	$varName        = $wd1;
	$VARNAME        = uc($varName) ;
#	if ( $varName eq "$SNKEY" ) { $IVAR_SNKEY = $NVAR; }
	if ( $VARNAME eq "$SNKEY" ) { $IVAR_SNKEY = $NVAR; }

	$VARNAMES_DICT[$NVAR] = $varName ;
	$NVAR++ ;
    }
    
  DONE_HEADER:
    if ( $IVAR_SNKEY < 0 ) {
	$MSGERR[0] = "Could not find required '$SNKEY' key after " ;
	$MSGERR[1] = "reading $NLINE_READ lines from " ;
	$MSGERR[2] = "$INFILE" ;
	sntools::FATAL_ERROR(@MSGERR);    
    }

    print " Found $NVAR variables  (IVAR_SNKEY = $IVAR_SNKEY). \n" ;

} # end of get_VARLIST


sub open_OUTFILE {

    my ($varName);
    my $OUTFILE = "SNANA_${INFILE}.FITRES" ;
    
    print " Open $OUTFILE \n";

    open PTR_OUTFILE , "> $OUTFILE" ;

    print PTR_OUTFILE "# Original dict file: $INFILE \n" ;
    print PTR_OUTFILE "# Translated into SNANA/fitres format using \n";
    print PTR_OUTFILE "# $0 \n";
    print PTR_OUTFILE "#  \n";

    # move SN key to start of list, and change it to CCID
    @VARNAMES_SNANA = ("CID");
    foreach $varName (@VARNAMES_DICT) {
	if ( uc($varName) eq "$SNKEY" ) { next ; }
	@VARNAMES_SNANA = ( @VARNAMES_SNANA , $varName );
    }

    print PTR_OUTFILE "NVAR: $NVAR \n";
    print PTR_OUTFILE "VARNAMES: @VARNAMES_SNANA \n\n";

} # end of open_OUTFILE ;


sub  update_OUTFILE {

    # break $line into $NVAR terms, put CID at the beggining
    # and write SNANA/fitres lines starting SN: $CCID etc ...
    # If first char is '#' of @ then return without update.

    ($line) = @_ ;

    my (@wdlist, $SNID, $ivar, $wd0, $NWD );

    @wdlist = split(/\s+/,$line);
    $wd0    = $wdlist[0] ;

    if ( $wd0             eq '#' ) { return ; }
    if ( substr($wd0,0,1) eq "@" ) { return ; }

    $NWD = scalar(@wdlist);
    if ( $NWD != $NVAR ) {
	$MSGERR[0] = "Found NWD = $NWD but expeted NVAR=$NVAR words";
	$MSGERR[1] = "Check line = " ;
	$MSGERR[2] = "  $line" ;
	sntools::FATAL_ERROR(@MSGERR);    
    }



    $SNID = $wdlist[$IVAR_SNKEY];
#    $SNID+= 0; # xxxxxxxxx  for integer id

    print PTR_OUTFILE "SN:  $SNID  ";

    # write every value EXCEPT for SNID
    for ( $ivar=0 ; $ivar < $NVAR; $ivar++ )  {
	if ( $ivar == $IVAR_SNKEY ) { next ; }
	print PTR_OUTFILE "$wdlist[$ivar] " ;
    }

    print PTR_OUTFILE "\n";

} # end of update_OUTFILE

