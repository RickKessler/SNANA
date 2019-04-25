#!/usr/bin/perl
#
# Created Jul 12, 2011 by R.Kessler
#
#
# Extract value of any variable for one SN in fitres file.
# Print fitres value to screen with format
#  $varName(CID=$cid) = $VALUE  (IVAR=$IVAR)
#
# This script is useful for visual inspection of values,
# and can also be called by other scripts to extract
# values for additional analysis; the $VALUE is the
# 3rd word of the output.
#
#
# Usage:
#  get_fitresValue.pl <fitresFile> <CID>  <varName> 
#
#  where 
#  <fitresFile> is the name of the fitres file
#  <CID>        is the name of the SN/candidate
#  <varName>    is the name of the variable to evaluate.
#
#
#  HISTORY
# -----------------------------------
#
#
# ----------------------------------

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

my ($NARG, @MSGERR,$fitresFile,$CID, $varName);
my ($IVAR, $VALUE);

$NARG = scalar(@ARGV); 
if ( $NARG < 3 ) {
    $MSGERR[0] = "Invalid command-line arguments.";
    $MSGERR[1] = "Must give fitresFile, CID and varName as arguments." ;
    sntools::FATAL_ERROR(@MSGERR);
}

$fitresFile = $ARGV[0] ;
$CID        = $ARGV[1] ;
$varName    = $ARGV[2] ;


# make sure that fitresFile exists
if ( !(-e $fitresFile ) ) {
    $MSGERR[0] = "Could not find fitres file";
    $MSGERR[1] = "$fitresFile";
    sntools::FATAL_ERROR(@MSGERR);
}

$IVAR  = &IVAR_FITRES($varName);
$VALUE = &VALUE_FITRES($CID,$IVAR);

print " $varName(CID=$CID)  = $VALUE    (IVAR=$IVAR)\n";


# =====================================
# ============= END OF MAIN ===========
# =====================================


# =====================
sub IVAR_FITRES {
    my ($varName) = @_ ;

    my ($tmp, $ivar, @words, @VARLIST) ;

    # return IVAR index of varName in VARNAMES list.
    # first variable (CID) has IVAR=1.

    @VARLIST = qx(grep VARNAMES $fitresFile);
    @words   = split(/\s+/,$VARLIST[0]) ;


    $ivar = -9 ;
    foreach $tmp ( @words ) {
        $tmp       =~ s/\s+$// ;   # trim trailing whitespace
	if ( $tmp eq "VARNAMES:" ) { $ivar = 0 ; next ; }

	$ivar++ ;
#	print " xxxx ivar=$ivar at $tmp \n";
	if ( $tmp eq $varName ) { return $ivar ; }
    }

    # if we get here, then abort
    $MSGERR[0] = "Could not find VARNAME = '${varName}'";
    sntools::FATAL_ERROR(@MSGERR);

} # end of IVAR_FITRES


# ==========================
sub VALUE_FITRES {
    my ($cid,$ivar) = @_ ;

    my ($string, @words, $NWORDS);

    # grep out the line with CID
    $string  = qx(grep " $cid " $fitresFile);
    @words   = split(/\s+/,$string) ;
    $NWORDS  = scalar(@words);
    if ( $NWORDS == 0 ) {
	$MSGERR[0] = "Could not find CID = '${cid}' in  $fitresFile";
	sntools::FATAL_ERROR(@MSGERR);
    }

    return $words[$ivar];
#    print " line(CID=$cid, NWORDS=$NWORDS) = @words \n";

} # end of VALUE_FITRES
