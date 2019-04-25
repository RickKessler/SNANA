#!/usr/bin/perl
#
# --------------
# SALT2train_translate.pl
#
# Written to translate SN simulations into SALT2 training format. 
#
# Usage:
#
# ./SALT2train_prep.pl  <input file>
#
# Given a light curve simulation input file and GENVERSION name, 
# script will translate that simulation into SALT2
# format in the directory $SNDATA_ROOT/SIM/SALT2-training/<GENVERSION> . 
# 
#  
###################
### For easy file navigation: sections are as follows:
#
### DECLARATIONS
### OUTPUT FILE FORMATS 
### MAIN CODE
### subroutines (listed individually)
#
##### HISTORY #####
#
# Sep 22 2011: (JLM) Altered run_DUMPFILES to correctly create
#                    trainingsample.list file independent of OS. 
#                    Used print statement in lieu of FORMAT for 
#                    improved readability. 
#
# OCT 17 2011 (JLM) Current version creates spectrum .DUMP file in 
#                    directory from which script is run.
#        
# OCT 18 2011 (JLM) Fixed small bug to make .DUMP format compatible with
#                    snana. 
#
# OCT 31 2011 (JLM) Spectrum .DUMP file now includes first wavelength, c, and x1. 
#
# NOV 7 2011 (JLM)  Parses correctly for either GRIDPARAM or PARAM x1 and c. 
#
# NOV 22 2011 (JLM) Fixed dumb bug in GRIDPARAM parse
#
# JUN 2 2012 (JLM) Make script quit if LC files are missing
#
# MAY 21 2013 (JLM) Convert nml file output from perl 'format' call
#                   to a series of 'print <FILE> "<text>";' 
# 
###################
# --------------

use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use sntools;

use File::Basename;
use File::Copy;
use POSIX;
use Cwd;
use strict;

sub make_banner($);
sub make_subbanner($);
sub check_fileexist(@);
sub check_goodexit(@);
#
sub get_inFile(@);
sub parse_inFile_SALT2train($);
sub parse_LCSIMLINE($);
#
sub init_UNPACK(@);
sub run_UNPACKLC();
sub run_DUMPFILES();
sub parse_LC($);
sub parse_LC_SALT2($);
sub parse_SPEC($);	     
sub prep_UNPACKLOC();
#
sub specfile_read($);
sub specfile_extract(@);
### DECLARATIONS ###

my $local_dir; ### current working directory
my $USER_inFile; ### User input file name
my $ABORT; ### Use for successful run check
my @MSGERR; ### Use for fatal error messages
###
###
my %LCSIM_HASH; ###
###
my $LCSIM_INPUT; ### User sim input file name
my $LCSIM_CALL; ### Full set of command line vars for lc_sim.exe
my $GENVERSION; ### sim input keyword value
my $GENMODEL; 
###
my $LCDUMP_FILE; ### lc sim file from which to grab SNe cids
my $SPECDUMP_FILE;
my $SPECSIM_dir; ### dir in which to put spec sims
my $UNPACK_dir; ### dir in which to unpack SALT2 formatted sims
###
my $SALT2_INSTRUMENT; ### required for unpack nml
my $SALT2_MAGSYS; ### required for unpack nml
my $SALT2_PREFIX; ### required for unpack nml
my $SALT2_FULLPREFIX; ### required for unpack nml ( includes SALT2 instrument )
my $OPT_REFORMAT_SALT2; ### required for unpack nml, default is 2
my $JOBNAME_SNANA; ### path to snana.exe
my $TRAINLIST_file; ### subroutine generates rfphase list of SEDS to make
###
my (@SPECLINES, @TEMP, $N_SPECFILE_DATA, $R_SPECFILE_DATA); ### globals for sub specfile_read 
###		
my $printprefix; ### all print names are for perl FORMAT calls
my $printphase;
my $printdate;
my $printz;
my $printCID;
my $printtype;
my $printfile;
my $printspecfile;


###  MAIN CODE ###

$local_dir = getcwd; 

$USER_inFile = get_inFile(@ARGV);

parse_inFile_SALT2train($USER_inFile);

while ( ($GENVERSION, $LCSIM_INPUT) = each(%LCSIM_HASH)) {

    init_UNPACK($USER_inFile, $LCSIM_INPUT);

    prep_UNPACKLOC();

    run_UNPACKLC();

    run_DUMPFILES();

}

print "FINISHED SUCCESSFULLY \n ";

###  END  ###
# ----------

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

}

    
sub init_UNPACK(@)
{
    my ($inFile, $SimInputFile) = @_;
    my $missing_die;
    my $nofile_die;
    my $LCSIMDIR;
    my $SIMREADME;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my $NARG7 = 7;
    my $NARGALL = 99;
    my @outlist; ### use for outputting EPOCH DIST variables
    my @tmp2; ### use for parsing @tmp in NARG2 cases
    my @tmp; ### consistent with rest of snana sntools.pl usage
    my $key; ### consistent with rest of snana sntools.pl usage
    my $GEN_SALT2INSTRUMENT;

    make_banner("PREPARING FOR UNPACKING");

    print "Grabbing required keys from input files. \n";

    ## parse_inFile_SALT2train has already run,
    ## GENVERSIONs have already been determined

    print "\nGrabbing required keys from $SimInputFile . \n\n";

    my $key = "SALT2_INSTRUMENT:" ;
    @tmp = sntools::parse_line($SimInputFile, $NARG1, $key, $OPTABORT) ;
    $SALT2_INSTRUMENT = @tmp[0];

    my $key = "SALT2_MAGSYS:" ;
    @tmp = sntools::parse_line($SimInputFile, $NARG1, $key, $OPTABORT) ;
    $SALT2_MAGSYS = @tmp[0];    

    print "\nGrabbing optional  keys. \n\n";

    my $key = "OPT_REFORMAT_SALT2:" ;
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $OPT_REFORMAT_SALT2 = @tmp[0];    
    if ( $OPT_REFORMAT_SALT2 eq '' ) { 
	$OPT_REFORMAT_SALT2 = 2; 
    }

    my $key = "JOBNAME_SNANA:" ; ## OPTIONAL
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $JOBNAME_SNANA = @tmp[0];
    if ( $JOBNAME_SNANA eq '' ) {
	$JOBNAME_SNANA = $ENV{'SNANA_DIR'} . "/bin" . "/snana.exe" ;
    }

    my $key = "SALT2_PREFIX:"; ### Must have this key to unpack
    @tmp = sntools::parse_line($inFile, $NARG1, $key, $OPTQUIET) ;
    $SALT2_PREFIX = @tmp[0];
    if ( $SALT2_PREFIX eq '' ) { 
	$SALT2_PREFIX = 'PHOTLC'; 
	print "\t USING DEFAULT SALT2 LC PREFIX \t $SALT2_PREFIX \n"; 
    } else {
	print "\t SALT2_PREFIX \t $SALT2_PREFIX \n"; 
    }

} ### end init_UNPACK($)

sub prep_UNPACKLOC()
{
    my $missing_die;
    my $nofile_die;
    my $LCSIMDIR;
    my $SIMREADME;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my $NARG7 = 7;
    my $NARGALL = 99;
    my @tmp2; ### use for parsing @tmp in NARG2 cases
    my @tmp; ### consistent with rest of snana sntools.pl usage
    my $key; ### consistent with rest of snana sntools.pl usage

    print "\nMaking required directory locations, checking for required files. \n";

    $UNPACK_dir = $ENV{'SNDATA_ROOT'} . "/SIM/SALT2-training/" . $GENVERSION;
    $LCSIMDIR= $ENV{'SNDATA_ROOT'} . "/SIM/" . $GENVERSION . "/";
    $LCDUMP_FILE = $LCSIMDIR . $GENVERSION . ".DUMP";
    $SIMREADME = $LCSIMDIR . $GENVERSION . ".README"; 

    print "\t LCSIM GENVERSION \t $GENVERSION \n";
    print "\t SALT2 TRAINING DIR $UNPACK_dir \n";

    ### MUST have access to dump file, readme file

    check_fileexist("LCSIM dump", $LCDUMP_FILE);
    print "\t LCSIM dump file located. \n";

    check_fileexist("LCSIM readme", $SIMREADME);
    print "\t LCSIM readme file located. \n";

    print "\nChecking for more required keys. \n";

    print "\nUSING required keys: \n\n";

    my $key = "BRIEF_DESCRIPTION:"; ### Must have this key to unpack
    @tmp = sntools::parse_line($SIMREADME, $NARG7, $key, $OPTABORT) ;
    @tmp2 = split(' ', @tmp[0]);
    $GENMODEL = @tmp2[6];
    print "\t GENMODEL \t $GENMODEL \n";

    print "\nDone finding required files and keys. \n"; 

    print "\nPreparing SALT2 TRAINING Directory. \n";

    ### CLEAN UP OR CREATE SALT2 TRAINING DIR
    if ( -d $UNPACK_dir ) { 
	print "\tCleaning up SALT2 TRAINING DIR \n"; 
	unlink(<$UNPACK_dir/*>);
    } else {
	print "\tCreating SALT2 TRAINING DIR \n"; 
	$nofile_die = "ABORT: Unable to create SALT2 TRAINING DIR. \n";
	mkdir($UNPACK_dir, 0777) || die $nofile_die;
    }

    $TRAINLIST_file = $UNPACK_dir . '/trainingsample.list'; 

} ### end prep_UNPACKLOC()

sub run_DUMPFILES()
{

    my @DUMPline;
    my $CID;
    my $LCSIM_DIR;
    my $filename;
    my $SEDSIM_lcfile; ## check for training lcfile
    my $LCDAT_FILE;
    my $LCSIM_datfile; ## use lc sim dat file to get params
    my $trainlist_file; ##
    my $readspec_file; 
    my $readlc_file;
    my @temp_LC_list; ## use to hold LC files initially
    my @temp_SPEC_list; ## use to hold SPEC files initially
    my $die_missing;
    my @LC_INFO;
    my @LC_PARAM;
    my $SPECDUMP_file;
    my $SPECDUMP_head1;
    my $SPECDUMP_head2;
    my ($nspec, $ispec, $spec_phase, $spec_date, $spec_firstwave);
    my @temp;
    my $name;
    my ($spectype, $specmjd, $firstwave, $lastwave);
    my ($bin1, $bin2, $firstdata_line, $wavecol);

    make_banner("CREATING < $GENVERSION > TRAININGSAMPLE.LIST FILE") ;

    print "\nMaking needed directory and file names\n";
    $LCSIM_DIR = $ENV{'SNDATA_ROOT'} . '/SIM/' . $GENVERSION . '/';
    $SPECDUMP_file = $local_dir . "/" . $GENVERSION . "_SPECINFO.DUMP"; 
    print "specdump_file is $SPECDUMP_file\n"; 
    $SPECSIM_dir = $ENV{'SNDATA_ROOT'} . '/SIM/SALT2-training/' . $GENVERSION . '/';
    $TRAINLIST_file = $SPECSIM_dir . '/trainingsample.list';


    print "\nReading SN IDs from LCSIM DUMP file. \n";
    print "\t Dump file: $LCDUMP_FILE \n";
    
    open(DUMPFILE, $LCDUMP_FILE) || die "ABORT: LCSIM DUMP FILE not found. \n ";
    open (DUMPOUT, ">$SPECDUMP_file"); 
    print DUMPOUT "NVAR: 10\n";
    print DUMPOUT "\n";
    print DUMPOUT "\#NOBS, PHMIN, PHMAX, SNRMAX are lc variables\n";
    print DUMPOUT "\#NSPEC is the number of spectra for that SN\n";
    print DUMPOUT "\#PHSPEC is the rest frame phase for the ISPEC spectrum\n"; 
    print DUMPOUT "\n";
    print DUMPOUT "VARNAMES: CID Z PEAKMJD NOBS PHMIN PHMAX SNRMAX NSPEC ISPEC PHSPEC FIRSTWAVE S2c S2x1 \n";
    open (TRAINLIST, ">$TRAINLIST_file");
    while(<DUMPFILE>){
        chomp;
        (@DUMPline) = split(' ', $_ );
        if ( $DUMPline[0] =~ /^SN:/ ) {
            ## for each SN, make file names based on CID
	    $CID = $DUMPline[1];
	    #initialize list vectors
	    @temp_LC_list = (); 
	    @temp_SPEC_list = (); 
	    opendir(DATA, $SPECSIM_dir) || die "Cannot open dir $SPECSIM_dir \n";
	    while ( $name = readdir(DATA)) {
		if ( $name =~ m/$CID/ && $name !~ m/spec/ ) {
		    push(@temp_LC_list, $name);
		} elsif ( $name =~ m/$CID/ && $name =~ m/spec/ ) {
		    push(@temp_SPEC_list, $name);
		}
	    }
	    closedir(DATA);
	    $printCID = $CID;
            $LCDAT_FILE = $LCSIM_DIR . $GENVERSION . '_SN0' . $CID . '.DAT';
	    @LC_PARAM = parse_LC_SALT2($LCDAT_FILE);
	    $printtype = "LC";		    
	    foreach (@temp_LC_list) { 
		$printfile = $_;
		$readlc_file = $SPECSIM_dir . $printfile;
		print TRAINLIST "$printCID $printtype $printfile \n";
		@LC_INFO = parse_LC($readlc_file);
		$nspec = @temp_SPEC_list;
		$ispec = 0;
		if ($nspec == 0 ) {
		    print DUMPOUT "SN: $CID @LC_INFO $nspec $ispec -99.9 -99.9 @LC_PARAM \n";
		}
	    }
	    $printtype = "SPEC";
	    foreach (@temp_SPEC_list) {
		$printfile = $_;
		print TRAINLIST "$printCID $printtype $printfile \n";
		$readspec_file = $SPECSIM_dir . $printfile;
		$ispec = $ispec + 1;
		($spec_date, $spec_firstwave) = parse_SPEC($readspec_file);
		$spec_phase = ($spec_date - $LC_INFO[1])/(1+$LC_INFO[0]);
		print DUMPOUT "SN: $CID @LC_INFO $nspec $ispec $spec_phase $spec_firstwave @LC_PARAM\n";
		### ($N_SPECFILE_DATA, $R_SPECFILE_DATA) = specfile_read($readspec_file);
		### @temp = specfile_extract($N_SPECFILE_DATA, $R_SPECFILE_DATA); 
		### ($spectype, $specmjd, $firstwave, $lastwave, $bin1, $bin2, $firstdata_line, $wavecol) = @temp;
	    }
	    if (@temp_LC_list < 1) {
		$MSGERR[0] = "LC file for $GENVERSION, SN $CID is missing.";
		$MSGERR[1] = "Check $LCSIM_DIR and $SPECSIM_dir for problems.";
		$MSGERR[2] = "ABORT";
		sntools::FATAL_ERROR(@MSGERR);
	    }

	}
    }
    
    close(TRAINLIST);
	close(DUMPFILE);
	close(DUMPOUT);
    
}

sub run_UNPACKLC()
{
    my $unpack_nmlfile; ### trash file
    my $unpack_logfile; ### trash file
    my $unpack_job;
    my @output;
    my $ABORT;
    my $ABORT_die;
    my $temp_hfile; ### trash file 
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG2 = 2;
    my @tmp; ### consistent with rest of snana sntools.pl usage
    my $key; ### consistent with rest of snana sntools.pl usage

make_banner("UNPACK < $GENVERSION > LIGHT CURVES AND SPECTRA INTO SPECSIM DIR");

print "\nSetting up locations of output files for writing and cleaning.\n";

    $temp_hfile = $UNPACK_dir . '/' . 'temp.his';
    $unpack_nmlfile = $UNPACK_dir . '/' . 'temp.nml';
    $unpack_logfile = $UNPACK_dir . '/' . 'unpack.log';

print "\nGenerating nml file. \n";
    
    $SALT2_FULLPREFIX = $SALT2_PREFIX . '_' . $SALT2_INSTRUMENT;

    $ABORT_die = "ABORT: Can't open nml file in SIMSED DIR\n";
    open(NMLFILE,">$unpack_nmlfile" ) or die $ABORT_die ;

    print NMLFILE  "&SNLCINP \n";
    print NMLFILE  "   VERSION_PHOTOMETRY = \'$GENVERSION\' \n";
    print NMLFILE  "   OPT_REFORMAT_SALT2 = $OPT_REFORMAT_SALT2 \n";
    print NMLFILE  "   REFORMAT_KEYS = \'\@INSTRUMENT $SALT2_INSTRUMENT \@MAGSYS $SALT2_MAGSYS \@PREFIX $SALT2_FULLPREFIX \' \n";
    print NMLFILE  "   HFILE_OUT = \'$temp_hfile\' \n";
    print NMLFILE  "   cutwin_cid = 0,100000 \n";
    print NMLFILE  "&END \n";

    close(NMLFILE);
    check_fileexist("Unpack nml", $unpack_nmlfile);

print "\nRUNNING LC UNPACKER IN SIMSED DIR\n";
    jobinfo_dump($JOBNAME_SNANA, "temp.nml", "unpack.log");

    @output = `cd $UNPACK_dir; $JOBNAME_SNANA temp.nml >& unpack.log ; cd $local_dir;` ;
    if ($? == -1 ) { 
       $ABORT_die = "\nABORT: Unpack routine failed. \n";
       die $ABORT_die; 
   }

    my $key = "ABORT" ; ## OPTIONAL
    check_goodexit("UNPACK", $unpack_nmlfile, $key); 

print "\nCleaning up unpack trash files\n";

   unlink($temp_hfile);
   unlink($unpack_nmlfile);
   unlink($unpack_logfile);

make_banner("UNPACK < $GENVERSION > LIGHT CURVES SUCCESSFUL");

}

sub parse_inFile_SALT2train($)
{
    use Switch;

    my $inFile_SALT2train = shift;
    my ( @line ) ; ## use to hold lines in file
    my $trash_value; 
    my $narg;
    my $iarg;
    my $key;
    my @tmp;
    my @tmp2;
    my $NARG1 = 1;
    my $NARGALL = 99;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $GEN_SALT2INSTRUMENT;
    my @VAR_ARRAY;


    ### ENSURE inFILE existence HAS ALREADY BEEN CHECKED

    open (INFILE_SALT2TRAIN, $inFile_SALT2train) or  die "ERROR: Cannot open input file $inFile_SALT2train \n";
    make_banner("Parsing input file: $inFile_SALT2train");
    

    while ( <INFILE_SALT2TRAIN> ){
        chomp;
        @line = split(' ', $_ );
        switch ($line[0]){
            case("LCSIM_INPUT:")    {
		parse_LCSIMLINE(\@line);
                if ( $GENVERSION eq '' ) {die "Must have GENVERSION as part of LCSIM_INPUT \n";}
                ### Make LCSIM_HASH, SPECSIM_HASH Entries Now ###
                ### Use Genversion as key - these should be unique ###
                $LCSIM_HASH{$GENVERSION} = $LCSIM_INPUT;
            }
            else {
		if ( $line[1] eq "LCSIM_INPUT:" ) {
		    $trash_value = shift(@line); 
		    parse_LCSIMLINE(\@line);
		    if ( $GENVERSION eq '' ) {die "Must have GENVERSION as part of LCSIM_INPUT \n";}
		    ### Make LCSIM_HASH, SPECSIM_HASH Entries Now ###
		    ### Use Genversion as key - these should be unique ###
		    $LCSIM_HASH{$GENVERSION} = $LCSIM_INPUT;
		}
	    }

	}
    }

    close(INFILE_SALT2TRAIN);

    make_banner("SIMULATIONS TO BE UNPACKED ARE:\n\n");

	while ( ($GENVERSION, $LCSIM_CALL) = each(%LCSIM_HASH)) {
	    print "\t\t$GENVERSION\n";
	}
		
}

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
    my $die_string; 

    $die_string = "\n\nABORT: Cannot find $key file $file \n\n";
    die $die_string unless -e $file; 

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

sub parse_LCSIMLINE($) 
{

    my $narg;
    my $iarg;
    my $key;
    my @tmp;
    my $NARG1 = 1;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
 
    my ( @lcsimline ) = @{$_[0]};

    $narg = @lcsimline; ## pull out number of elements in line vector

    if ( $narg == 2 ) {
	$LCSIM_CALL = @lcsimline[1];
	$LCSIM_INPUT = @lcsimline[1];

        ### parse LCSIM_INPUT file for GENVERSION ###
	check_fileexist("LCSIM_INPUT", $LCSIM_INPUT);
	my $key = "GENVERSION:" ;
	@tmp = sntools::parse_line($LCSIM_INPUT, $NARG1, $key, $OPTABORT) ;
	$GENVERSION = @tmp[0];
        ### done parsing ###
	if ( $GENVERSION eq '' ) {die "GENVERSION key not found in file $LCSIM_INPUT \n";}

    } elsif ($narg > 2 ) {
	$LCSIM_INPUT = @lcsimline[1];
	$LCSIM_CALL = join(' ', @lcsimline[1...$narg-1]);
	$iarg = 0;

        ### get GENVERSION from LCSIM_INPUT vector ###
	while ($iarg < $narg ) {
	    if ( @lcsimline[$iarg] eq 'GENVERSION' ) { $GENVERSION = @lcsimline[$iarg+1] }
	    $iarg = $iarg + 1;
	}

    }
}

sub specfile_read($)
{
    my  ($specfilename) = shift;
    my ($rspeclines);
    my ($nspeclines);
    print "\nOpening $specfilename. \n";
    open(SPECTRAFILE, $specfilename) || die "Couldn't open SPECTRAFILE $specfilename \n";
    @SPECLINES = <SPECTRAFILE>;
    close(SPECTRAFILE);
    $rspeclines = \@SPECLINES;
    $nspeclines = @SPECLINES;
    return ( $nspeclines, $rspeclines );
}

sub specfile_extract(@) 
{
    my ($nlines, $rlines) = @_ ;
    my ( $ilines ) = 0;
    my ( $count ) = -1;
    my ( $wavecol ) = 0;
    my ( $waveflag ) = 0;
    my ( @checkhead ) ;
    my ( $ncheckhead );
    my ( $icheckhead );
    my ( $SPECTYPE );
    my ( $SPECMJD );
    my ( $FIRSTWAVE );
    my ( $PREVWAVE );
    my ( $LASTWAVE );
    my ( $BIN1 );
    my ( $BIN2 );
    my ( $firstdata_line );
    my ( @data );

    while ( $ilines < $nlines) {
        chomp($rlines->[$ilines]);
        (@checkhead) = split(' ', $rlines->[$ilines]);
        ### START BY READING SPECTRUM FILE HEADER
        if ($checkhead[0] ne "" && $checkhead[0] =~ /^\@/){
                 ###
                 ### FIRST LINES STARTING WITH AMPERSANDS
                 ### EXTRACT SPECTYPE and Date.
                 ###
	    print ("amp line is @checkhead \n");
	    $ncheckhead = @checkhead;
	    $icheckhead = 0;
	    while ( $icheckhead < $ncheckhead){
		if ( $checkhead[$icheckhead]  =~ /FORMAT$/ ){
		    $SPECTYPE = cleanhead("@", $checkhead[$icheckhead + 1]);
		    print "spectype is $SPECTYPE\n";
		} elsif ($checkhead[$icheckhead] =~ /Date$/ ){
		    $SPECMJD = $checkhead[$icheckhead + 1];
		    print "date is $SPECMJD\n";
		}
		++$icheckhead;
	    }
	}  elsif ($checkhead[0] ne "" && $checkhead[0] =~ /^\#/){
                 ###
                 ### SECOND, LINES STARTING WITH POUND SIGNS
                 ### THESE LINES TELL FILE DATA COLUMN FORMAT
                 ### USE TO DETERMINE WAVELENGTH COLUMN
                 ### ASSUME THAT FLUX, FLUX_ERR ARE ALWAYS ADJACENT TO WAVELENGTH
                 ###
	    print ("pound line is @checkhead \n");
	    $ncheckhead = @checkhead;
	    $icheckhead = 0;
	    ++$count;
	    while ( $icheckhead < $ncheckhead) {
		if ( $checkhead[0] eq "#" && $ncheckhead == 1 ){
		    print "$count\n";
		    $count = $count -1 ;
		    print "$count\n";
		}
		if ( $checkhead[$icheckhead]  =~ /^WAVE/ || $checkhead[$icheckhead] =~ /WAVE$/ ){

		    $wavecol = $count;
                         print "wavecol is $wavecol\n"
			 }
		++$icheckhead;
	    }

	} else  {
                           ### READING THROUGH DATA
                           ### GET firstdata_line, FIRSTWAVE, LASTWAVE, BIN
                           ### EVENTUALLY NEED TO GET SNR AT VARIOUS INTERVALS
	    if ( $wavecol > -1 ){
		(@data) = split(' ', $rlines->[$ilines]);
		if ( $waveflag == 0 ) {
		    print "wave: $data[$wavecol], flux: $data[$wavecol + 1] ($data[$wavecol+2])\n";
		    $FIRSTWAVE = $data[$wavecol];
		    $firstdata_line = $ilines;
		    ++$waveflag;
		} elsif ($waveflag == 1){
		    $BIN1 = $data[$wavecol] - $FIRSTWAVE;
		    $LASTWAVE = $data[$wavecol];
		    ++$waveflag;
		    $ilines = $nlines - 10;
		} elsif ($ilines >= $nlines - 10 && $data[$wavecol] ne "" ){
		    $PREVWAVE = $LASTWAVE;
		    $LASTWAVE = $data[$wavecol];
		    $BIN2 = $LASTWAVE - $PREVWAVE;
		} else {
		    print "problem with spec_extract elseif loop \n";
		}

	    }
	}

	++$ilines;

    }


    return ($SPECTYPE, $SPECMJD, $FIRSTWAVE, $LASTWAVE, $BIN1, $BIN2, $firstdata_line, $wavecol);
}


sub parse_SPEC($)
{
    my $SPECtrainfile = shift;
    my @SPECLINE;
    my $date;
    my $firstwave;
    my $key;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my $NARG7 = 7;
    my @tmp;
    my $nline;

    my $key = "\@Date";
    @tmp = sntools::parse_line($SPECtrainfile, $NARG1, $key, $OPTABORT) ;
    $date = @tmp[0];

    open(SPECFILE, $SPECtrainfile);
    $nline = 0;
    while(<SPECFILE>){
	chomp;
	(@SPECLINE) = split(' ', $_ );
    ### skip spectrum file header lines
	unless ( $SPECLINE[0] =~ /^\@/ or $SPECLINE[0] =~ /^\#/ ) {
	     ### print "$nline, @SPECLINE \n";
            if ($nline == 0) {
		$firstwave = $SPECLINE[0];
		### print "file is $SPECtrainfile, firstwave is $firstwave\n";
		$nline = 1;
	    }
	}
    }
    close(SPECFILE);

    return ($date, $firstwave); 

}

sub parse_LC_SALT2($)
{
    my $LCdatfile = shift;
    my ( $s2c, $s2x1);
    my $key;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my @tmp;

    my $key = "SIMSED_GRIDONLY(S2x1):" ;
    @tmp = sntools::parse_line($LCdatfile, $NARG1, $key, $OPTQUIET) ;
    $s2x1 = @tmp[0];
    if ( $s2x1 eq '') {
	my $key = "SIMSED_PARAM(S2x1):" ;
	@tmp = sntools::parse_line($LCdatfile, $NARG1, $key, $OPTQUIET) ;
	$s2x1 = @tmp[0];	
    }

    my $key = "SIMSED_GRIDONLY(S2c):" ;
    @tmp = sntools::parse_line($LCdatfile, $NARG1, $key, $OPTQUIET) ;
    $s2c = @tmp[0];
    if ( $s2c eq '') {
	my $key = "SIMSED_PARAM(S2c):" ;
	@tmp = sntools::parse_line($LCdatfile, $NARG1, $key, $OPTQUIET) ;
	$s2c = @tmp[0];	
    }

    
    return ($s2c, $s2x1);
}

sub parse_LC($)
{
    my $LCtrainfile = shift;
    my ( $z, $peakmjd, $nobs, $tmin, $tmax, $snrmax ); 
    my $key;
    my $OPTABORT = 0;
    my $OPTWARN = 1;
    my $OPTQUIET = 2;
    my $NARG1 = 1;
    my @tmp;
    
    my $key = "\@Z_HELIO" ;
    @tmp = sntools::parse_line($LCtrainfile, $NARG1, $key, $OPTABORT) ;
    $z = @tmp[0];

    my $key = "\@DayMax" ;
    @tmp = sntools::parse_line($LCtrainfile, $NARG1, $key, $OPTABORT) ;
    $peakmjd = @tmp[0];

    my $key = "\@NOBS" ;
    @tmp = sntools::parse_line($LCtrainfile, $NARG1, $key, $OPTABORT) ;
    $nobs = @tmp[0];

    my $key = "\@TRESTMIN" ;
    @tmp = sntools::parse_line($LCtrainfile, $NARG1, $key, $OPTABORT) ;
    $tmin = @tmp[0];

    my $key = "\@TRESTMAX" ;
    @tmp = sntools::parse_line($LCtrainfile, $NARG1, $key, $OPTABORT) ;
    $tmax = @tmp[0];

    my $key = "\@SNRMAX" ;
    @tmp = sntools::parse_line($LCtrainfile, $NARG1, $key, $OPTABORT) ;
    $snrmax = @tmp[0];

    return($z, $peakmjd, $nobs, $tmin, $tmax, $snrmax);

}
