#!/usr/bin/perl
#--------------
# script designed to extend wavelength range of
# SALT2 training SALT2color_dispersion.dat file 
# in order to use it for light curve sims. 
# 
# required inputs:
# 1) Name of file to be extended
# 2) Name of output file 
#
# JMosher 03/23/13
#
#########################

use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use sntools;

use File::Basename;
use File::Copy;
use POSIX;
use Cwd;
use strict;

sub get_command_args(@);
sub load_file(@);
sub naive_fix(@);

##### DECLARATIONS #####

my ($INFILE, $OUTFILE, $headstring);
my (@MSGERR);
my (@ARR_DAT, @ARR_HEAD);
my ($item);

##### MAIN ######

($INFILE, $OUTFILE) = get_command_args(@ARGV);

#define SALT2 file header delimiter
$headstring = "#";

load_file($INFILE,\@ARR_DAT, \@ARR_HEAD, $headstring);

naive_fix(\@ARR_DAT);

#write output
open (MYFILE, ">$OUTFILE");
foreach $item (@ARR_HEAD){
    print MYFILE "$item";
}
foreach $item (@ARR_DAT){
    print MYFILE "$item";
}
close (MYFILE);

print "Finished successfully. \n";

##### END OF MAIN #####

# ==========================================
sub naive_fix(@)
{
    my ($p_datarray) = (@_);

    my ($firstval, $lastval, $val, @wdList);
    my ($step);
    

    print "Extending file wavelength range ... \n";

    $step = 0.05;

    (@wdList) = split(' ', $$p_datarray[0]);
    if ( $wdList[0] == 2500 ) {
	$firstval = $wdList[1];
	$val = $firstval + $step;
	unshift(@$p_datarray, "2450 $val \n");
	$val = $val + $step;
	unshift(@$p_datarray, "2200 $val \n");
	$val = $val + $step;
	unshift(@$p_datarray, "2000 $val \n");
    } else {
	$MSGERR[0] = "SALT2extend_colordisp.pl infile ";
	$MSGERR[1] = "has unexpected value $wdList[0]. ";
	$MSGERR[2] = "Expected 2500. This code hasn't ";
	$MSGERR[3] = "been tested for other starting ";
	$MSGERR[4] = "wavelengths. ABORT. ";
	sntools::FATAL_ERROR(@MSGERR);
    }
	    

    $step = 0.005;

    (@wdList) = split(' ', $$p_datarray[-1]);
     if ( $wdList[0] == 8000 ) {
	$lastval = $wdList[1];
	$val = $lastval + $step;
	push(@$p_datarray, "9000 $val \n");
	$val = $val + $step;
	push(@$p_datarray, "10000 $val \n");
	$val = $val + $step;
	push(@$p_datarray, "11000 $val \n");
    } else {
	$MSGERR[0] = "SALT2extend_colordisp.pl infile ";
	$MSGERR[1] = "has unexpected value $wdList[0]. ";
	$MSGERR[2] = "Expected 8000. This code hasn't ";
	$MSGERR[3] = "been tested for other ending ";
	$MSGERR[4] = "wavelengths. ABORT. ";
	sntools::FATAL_ERROR(@MSGERR);
    }
    
    
} # end naive_fix(@)

# ==========================================

sub load_file(@)
{
    
    my ($infile, $p_arrdat, $p_arrhead, $headstring) = (@_);
    my ($ncol, $nline);
    my (@temparray);

    #scoop up the entire file for linear parsing
    #--a la SALT2train_pipeline.pl, RK-- 

    print "\nLoading file $infile ... \n";

    @temparray = qx(cat ${infile});
    
    ($ncol, $nline) = parse_filedata(\@temparray, $p_arrdat, $p_arrhead, $headstring);
    
} #end of load_file(@)

# ==========================================

sub get_command_args(@)
{

    my $NARG = scalar(@_);
    my ($j, $comment, $infile, $outfile);

    #check that user gave required inputs
    print "\nStarting SALT2extend_colordisp.pl ... \n";

    if ($NARG < 2) {
        $MSGERR[0] = "\nNot enough command line inputs \n";
	$MSGERR[1] = "Must enter name of file to be extended";
	$MSGERR[2] = "and desired output file name.";
	$MSGERR[3] = "Check inputs and try again.";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $infile = $ARGV[0];
    $outfile = $ARGV[1];

    ### Make sure file exists
    $comment = "SALT2extend_colordisp.pl input file";
    sntools::check_fileexist($infile, $comment);

    print "\t Extending wavelength range of $infile ... \n";
    
    return ($infile, $outfile);
    
} #end get_command_args(@)

# ==========================================

sub parse_filedata(@)
{
    my ($p_fitresdump, $p_arrtofill, $p_arrhead, $headstring) = (@_);
    my ($iline, $tmpLine, @wdList, $wd0, $itrain);
    my ($nline, $N_LINES, $N_COLS);

    # figure out number of lines in file
    $nline = scalar(@$p_fitresdump);

    print "\tNumber of lines to parse: $nline ... \n";
    
    $N_LINES = 0;
    $iline = 0;
    while ($iline < $nline) {
	$tmpLine = "$$p_fitresdump[$iline]";
	@wdList = split(/\s+/, $tmpLine);
	if ($wdList[0] =~ /$headstring/ ) {
	    push(@$p_arrhead, $tmpLine);
	} else {
	    $N_COLS = @wdList;
	    push(@$p_arrtofill, $tmpLine); 
	    $N_LINES++;
	}

	$iline++;
    }
    
    return ($N_LINES, $N_COLS);

} #end parse_filedata(@)
