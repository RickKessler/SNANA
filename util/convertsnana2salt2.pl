#!/usr/bin/perl
#
# Nov 2011, R.Kessler
# Use snana.exe job to translate SNANA-formatted data files
# into SALT2 format. Also copy SNANA-auxilliary files into 
# SALT2 directory. Finally make tarball.
# By default, SALT2 files are put into a new directory named 
# $VERSION-SALT2format. The -salt2dir key can be used to specify
# a specific directory name.
#
# If the numbers of files exceeds MXSNLC then use 'split_version.pl'
# to run snana jobs that don't exceed the bound.
#
# Usage
#  convertsnana2salt2.pl  <nmlFile>
#
# Optional usage:
#
#  convertsnana2salt2.pl  <nmlFile> -salt2dir <salt2dir> 
#      (specify your directory name for salt2 files)
#
#  convertsnana2salt2.pl  <nmlFile> NOTAR  
#     (do not create tarball)
#
#
# The nmlFile must contain
#
# &SNLCINP
#     VERSION_PHOTOMETRY = '<VERSION>' 
#     OPT_REFORMAT_SALT2 = 2
#     REFORMAT_KEYS  = 
#       '@SURVEY <SURVEY>  @INSTRUMENT <INSTR>  @MAGSYS <MAGSYS>'
#
#     OPT_SETPKMJD = 4  ! optional: do NOT abort if PKMJD can't be found
# &END
#
#
#    HISTORY
#   ~~~~~~~~~~~~~
#
# Feb 01, 2012:  use new utility get_namelistValue()
#
# Oct 10, 2013:  
#    - new optional command-line input -salt2dir <salt2dir>
#    - new optional command-line input NOTAR
#    - check nmlFile for PRIVATE_DATA_PATH
#
# Dec 2 2015: fix bug found by DanS. Allow $SNDATA_ROOT/lcmerge/$VERSION/
#             directory
#    
# ---------------

use FindBin qw($Bin);
use lib "$Bin";
use sntools ;
use strict ;


sub parse_args ;
sub get_VERSION ;
sub get_DATADIR ;
sub get_NSPLIT ;
sub make_SALT2DIR ;
sub make_nmlFile ;
sub translate ;
sub make_tarBall ;
sub cleanup ;

# define globals
my ($nmlFile, $VERSION, @VERSION_LIST, $SIMFLAG, $TARFLAG);
my ($SALT2DIR, @MSGERR, $tarName, $DATADIR, $TOPDIR, $NFILE, $NSPLIT );

my $SNDATA_ROOT  = $ENV{'SNDATA_ROOT'};
my $SNANA_DIR    = $ENV{'SNANA_DIR'};

# =============== BEGIN MAIN ===============

&parse_args();

&get_VERSION();  # get name of photom version from nmlfile

&get_DATADIR();  # sim or data directory

&get_NSPLIT();

&make_SALT2DIR();

&translate();

&make_tarBall();

&cleanup();

# =========== END OF MAIN ==============


sub parse_args {

    my ($NARG, $i, $ii);

    $SALT2DIR = "" ;
    $nmlFile  = "" ;
    $TARFLAG  = 1  ; # make tarball by default

    $NARG = scalar(@ARGV);
    if ( $NARG < 1 ) {
	$MSGERR[0] = "Must give name of namelist file as argument.";
	sntools::FATAL_ERROR(@MSGERR);
    }
    $nmlFile    = $ARGV[0] ;

    for($i=1; $i < $NARG; $i++ ) {
	$ii = $i + 1 ;

	if ( $ARGV[$i] eq  "-salt2dir" ) { $SALT2DIR = $ARGV[$ii]; }
	if ( $ARGV[$i] eq "--salt2dir" ) { $SALT2DIR = $ARGV[$ii]; }

	if ( $ARGV[$i] eq "NOTAR" ) { $TARFLAG = 0 ; }
    }

} # end of parse_args



sub get_VERSION {

    # extract default version from &SNLCINP  namelist
    my $key     = "VERSION_PHOTOMETRY" ;
    $VERSION = sntools::get_namelistValue($nmlFile,$key);

    print " Prepare to translate $VERSION into SALT2 format.\n";

} # end of get_VERSION


sub get_DATADIR {

    my ($SIMDIR, $TMPDIR, $TMPFILE, $key, $PRIVATE_DATA_PATH );

    # check for private data path
    $key       = "PRIVATE_DATA_PATH" ;
    $PRIVATE_DATA_PATH = sntools::get_namelistValue($nmlFile,$key);

    $SIMDIR   = "$SNDATA_ROOT/SIM/$VERSION" ;
    $TMPDIR   = "$SNDATA_ROOT/lcmerge/$VERSION";
    $TMPFILE  = "$SNDATA_ROOT/lcmerge/$VERSION.README";

    if ( -d $SIMDIR ) {
	$DATADIR = "$SIMDIR" ;
	$TOPDIR  = "$SNDATA_ROOT/SIM" ;
    }
    elsif ( length($PRIVATE_DATA_PATH) > 0 ) {
	$DATADIR = "$PRIVATE_DATA_PATH" ;
	$TOPDIR  = "$DATADIR" ;
    }
    elsif ( -d $TMPDIR ) {
	$DATADIR = "$TMPDIR" ;
	$TOPDIR  = "$DATADIR" ;
    }
    elsif ( -e $TMPFILE ) {
	$DATADIR = "$SNDATA_ROOT/lcmerge" ;
	$TOPDIR  = "$DATADIR" ;
    }

   	
    print " DATADIR : $DATADIR \n";

} # end of get_DATADIR


sub get_NSPLIT {

    my ($listFile, @bla, $tmp, $MXSNLC, $snanaCode );
    # Now get number of data files from list file
    $listFile = "$DATADIR/${VERSION}.LIST";
    @bla = `cat $listFile`;
    $NFILE = scalar(@bla);

    # get MXSNLC from snana code
    $snanaCode = "$SNANA_DIR/src/snana.car" ;
    $tmp  = qx(grep ",MXSNLC " $snanaCode) ;
    @bla = split(/\s+/,$tmp) ;
    $MXSNLC = $bla[4] ;

##    $MXSNLC = 100 ; # DDDDDDDDDDDDDDD

    $NSPLIT = int( $NFILE / $MXSNLC) + 1 ;
    print " Number of data files to translate: $NFILE ==> $NSPLIT split jobs.\n";
    
    if ( $NSPLIT > 1 ) {
	qx(split_version.pl $VERSION $NSPLIT) ;
	@VERSION_LIST = qx(cd $TOPDIR; ls -d SPLIT*$VERSION );
    }
    else {
	@VERSION_LIST = "${VERSION}" ;
    }


} # end of get_NSPLIT


sub make_SALT2DIR {

    # set default SALT2DIR only if it is NOT already set; otherwise
    # leave it.
    if ( $SALT2DIR eq "" ) 
    {  $SALT2DIR = "${VERSION}-SALT2format"; }

    print " Create SALT2 directory: $SALT2DIR \n";

    if ( -d $SALT2DIR) { qx(rm -r $SALT2DIR); }
    qx(mkdir $SALT2DIR);

    qx(cp $nmlFile $SALT2DIR/ );
    qx(cp -r $DATADIR/$VERSION.*  $SALT2DIR);

    # remove SNANA's list file
    my ($LISTFILE);
    $LISTFILE = "$SALT2DIR/${VERSION}.LIST";
    if ( -e $LISTFILE ) { qx(rm $LISTFILE); }

}  # end of make_SALT2DIR


sub translate {

    my ($VERS, $arg, $rm, $logFile );

    foreach $VERS ( @VERSION_LIST ) {
        $VERS    =~ s/\s+$// ;   # trim trailing whitespace
	$arg     = "VERSION_PHOTOMETRY $VERS" ;
	$logFile = "${VERS}.log";
	print "\t Translating $VERS ... \n";
	qx(cd $SALT2DIR; snana.exe $nmlFile $arg > $logFile );
    }

    if ( $NSPLIT > 1 ) {
	$rm = "rm -r $TOPDIR/SPLIT*$VERSION* " ;
	qx($rm);
#	print " xxxx $rm \n";
    }

} # end of translate

sub make_tarBall {
    # create tarball
    
    if ( $TARFLAG == 0 ) { return ; }

    $tarName      = "${SALT2DIR}.tar";

    my $tmpFile = "${tarName}.gz";
    if ( -e $tmpFile ) { qx(rm $tmpFile); }

    qx(tar -cf $tarName $SALT2DIR ; gzip $tarName);
    print " Translated SALT2 files are in ${tmpFile} \n";

} # end of make_tarBall


sub cleanup {
    # Oct 2013: remove junk file(s)
    qx(rm $SALT2DIR/MNFIT*.LOG);
} 
