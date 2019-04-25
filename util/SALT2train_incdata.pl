#!/usr/bin/perl
#
# --------------
# SALT2train_incdata.pl
#
# Written to copy real SN data into SIM/SALT2-training/ for SIM+DATA trains,
# enabling SALT2train_run.pl to treat real SN data as a SIM for 
# training purposes.
# 
# Usage:
#
# ./SALT2train_incdata.pl <input file>
#
# INPUT FILE KEYWORDS:
# 
# INCDAT_INPUT: <path/filename>
# INCDAT_ARGUMENTS: GENVERSION <GENVERSION(optional)> DATAPATH <path(optional)>
# END:
#
# Given the full address & name of a SALT2 traininglist file, <path/filename>,
#  a GENVERSION name, and an optional DATAPATH, 
#  this script will make directory  $SNDATA_ROOT/SIM/SALT2-training/<GENVERSION> ("TRAINDIR"),
#  copy all files from <path> to TRAINDIR, 
#  and copy SALT2 traininglist file <path/filename> to TRAINDIR/trainingsample.list .
#
##### HISTORY ##### 
#
# Created Sep 6 2012 to enable hybrid DAT-SIM SALT2 trainings. 
# 
# SEP 11 2012: altered check_fileexist syntax to ensure log file error message. 
#
###################
#------------------

use IO::Handle;
use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use sntools ;

use File::Copy;
use POSIX;
use Cwd;
use strict;

sub make_banner($);
sub make_subbanner($);
sub check_fileexist(@);
sub check_goodexit(@);

sub get_inFile(@);
sub parse_inFile($);
sub init_UNPACK(@);
sub prep_INCDATLOC($);

### DECLARATIONS ###
my ($NINCDATA, @INCDATA_DATAPATH);
my (@INCDATA_GENVERSION, @INCDATA_FILENAME);
my ($local_dir, $USER_inFile, $i);
my ($copydir, $fxcopy);
my @MSGERR;

### MAIN ###

$local_dir = getcwd;
$USER_inFile = get_inFile(@ARGV);
parse_inFile($USER_inFile);

for ($i=1; $i<=$NINCDATA; $i++) {
    $copydir = prep_INCDATLOC($i); ## prepare SIM/SALT2-training DIR
    qx(cp $INCDATA_DATAPATH[$i]/* $copydir); ## copy all files to new DIR
    qx(cp $INCDATA_FILENAME[$i] $copydir/trainingsample.list); ## copy trainlist file
}

print "\n FINISHED CREATING INCDATA DIRECTORIES \n\n";

### END OF MAIN ###


sub prep_INCDATLOC($)
{
    my $iincdat = shift;
    my ($UNPACK_dir, $genversion, $nofile_die);
    
    print "\n PREPPING DIRECTORY LOCATION FOR SET $iincdat \n";
    
    $genversion = $INCDATA_GENVERSION[$iincdat];
    $UNPACK_dir = $ENV{'SNDATA_ROOT'} . "/SIM/SALT2-training/" . $genversion;
    
    ### CLEAN UP OR CREATE SALT2 TRAINING DIR
    if ( -d $UNPACK_dir ) {
	print "\tCleaning up $UNPACK_dir \n";
        unlink(<$UNPACK_dir/*>);
    } else {
        print "\tCreating SALT2 TRAINING DIR \n";
        $nofile_die = "ABORT: Unable to create SALT2 TRAINING DIR. \n";
        mkdir($UNPACK_dir, 0777) || die $nofile_die;
    }
    return ($UNPACK_dir);

} #end prep_INCDATLOC()

sub parse_inFile($)
{
    use Switch;

    my $inFile_SALT2train = shift;
    my ( @line ) ; ## use to hold lines in file
    my ($trash_value, @VAR_ARRAY);
    my ($narg, $iarg, $key, @tmp, @tmp2);
    my ($iword, $nword, $islash, $temp_path);
    my ($iincdata);

    ### ENSURE inFILE existence HAS ALREADY BEEN CHECKED                                                                                                                                                                                               
    open (INFILE_SALT2TRAIN, $inFile_SALT2train) or  die "ERROR: Cannot open input file $inFile_SALT2train \n";
    make_banner("Parsing input file: $inFile_SALT2train");

    $NINCDATA = 0; ### initialize INCDAT counter

    while ( <INFILE_SALT2TRAIN> ){
        chomp;
        @line = split(' ', $_ );
        switch ($line[0]){
	    case("INCDAT_INPUT:") {
		$NINCDATA++;
		$INCDATA_FILENAME[$NINCDATA] = $line[1];
		# strip off name after last slash to obtain datapath
		$islash = rindex($INCDATA_FILENAME[$NINCDATA], "/");
		$temp_path = substr($INCDATA_FILENAME[$NINCDATA], 0, $islash);
		$INCDATA_DATAPATH[$NINCDATA] = $temp_path;
	    } # end case INCDAT_INPUT:
	    case("INCDAT_ARGUMENTS:") {
		$nword = @line;
		for ( $iword=1; $iword <= $nword; $iword++ ){
		    if ( $line[$iword] eq "GENVERSION" )    {
			$iword++;
			$INCDATA_GENVERSION[$NINCDATA] = $line[$iword];
		    } elsif ( $line[$iword] eq "DATAPATH" ) {
			$iword++;
			### override INCDATA_DATAPATH with opt input
			$INCDATA_DATAPATH[$NINCDATA] = $line[$iword];
		    } 
		}
	    } # end case INCDAT_ARGUMENTS:
	} # end of switch
    } # end of while

	### check input information
	for ( $iincdata=1; $iincdata<=$NINCDATA; $iincdata++){
	    if ($INCDATA_GENVERSION[$NINCDATA] eq "") {
		die "\n TRAINFILE: $INCDATA_FILENAME[$iincdata] missing GENVERSION keyword\n";
	    }
	    print "\t\tChecking that trainlist file $INCDATA_FILENAME[$iincdata] exists.\n";
	    check_fileexist("trainlist file", $INCDATA_FILENAME[$iincdata]);
	    print "\t\tChecking that DATAPATH $INCDATA_DATAPATH[$iincdata] exists.\n\n";
	    sntools::checkDirExist($INCDATA_DATAPATH[$iincdata],"INCDATA data path");
	} # end of for loop

	### summarize incdats
	print "\nCREATING $NINCDATA HYBRID SIM-DATA TRAINING SETS\n\n";
	for ( $iincdata=1; $iincdata<=$NINCDATA; $iincdata++){
	    print "\t\t SET $iincdata: \n";
	    print "\t\t\t GENVERSION: $INCDATA_GENVERSION[$iincdata]\n";
	    print "\t\t\t DATAPATH:   $INCDATA_DATAPATH[$iincdata]\n";
	    print "\t\t\t TRAINFILE:  $INCDATA_FILENAME[$iincdata]\n";
	} # end of for loop
    print "\n";

} # end sub parse_inFile

sub get_inFile(@)
{

    my $NARG = scalar(@_);
    my $tempFile;

    ## check that user gave input filename                                                                                                                                                                                                             
    if ( $NARG < 1 )  {
	print " Must give input filename as argument. \n" ;
	die " ***** ABORT ***** \n" ;
    }

    $tempFile = $ARGV[0];

    check_fileexist("input", $tempFile);

    return($tempFile);

} #end sub get_inFile(@)

sub make_banner($)
{

    my $banner_string = shift;
    my $SYMLINE = "\n _____________________________ \n\n" ;

    print "$SYMLINE";
    print "$banner_string";
    print "$SYMLINE";

} # End of make_banner

sub make_subbanner($)
{
    my $banner_string = shift;
    my $SYMLINE = "\n\t\t _____________________________ \n\n" ;

    print "$SYMLINE";
    print "\t\t $banner_string";
    print "$SYMLINE";

}

sub jobinfo_dump(@)
{ 
    my ($job, $input, $log) = @_;
    
    # if no job, log, or input is required, call with "none" in the relevant place

    print "JOB CALLED AS: \n";
    print "\t <JOB NAME> <JOB INPUT>  >& <JOB LOG> \n\n";
    print "\t USING \n";
    unless ($job eq 'none') {print "\t JOB NAME:  $job \n";}
    unless ($job eq 'none') {print "\t JOB INPUT: $input \n";}
    unless ($job eq 'none') {print "\t JOB LOG:   $log \n";}
}


sub check_fileexist(@)
{
    my ($key, $file) = @_;
    my $die_flag = 0;

    if ( -e $file) {$die_flag = 1;}

    unless ($die_flag) {
	$MSGERR[0] = "Cannot find $key file ";
	$MSGERR[1] = "     '$file'";
	$MSGERR[2] = "Check file name and path,";
	$MSGERR[3] = "and try again.";
	sntools::FATAL_ERROR(@MSGERR);
    }

} # End of check_fileexist

sub check_goodexit(@)
{
    ## ideally here, "key" is the name of the job whose log is being checked
    my $JOBname;
    my $logfile;
    my $ABORTkey;
    my @tmp; # for parse_line
    my $OPTQUIET = 2; 
    my $ABORT; 
    my $die_string;

    ($JOBname, $logfile, $ABORTkey) = @_;
    if ( $ABORTkey eq '' ) { $ABORTkey = "ABORT"; }
    
    @tmp = sntools::parse_line($logfile, 99, $ABORTkey, $OPTQUIET) ;
    $ABORT = @tmp; ### use number of entries in @tmp to test ABORT
    $die_string = "$JOBname job failed. ABORT simprep. \n\n";
    if ( $ABORT > 0 ) { 
	print "$ABORTkey string found in log $logfile .\n"; 
	die $die_string ;}

    
} # End of check_goodexit
