#!/usr/bin/perl
#
# Created July 13 2010 by R.Kessler
#
# Split original "VERSION" into smaller sub-versions
# to avoid the MXSNLC limit, and to simplify breaking
# a long job into multiple jobs.  If the split version
# size exceeds MXSNLC a warning is given ... it does
# not abort because non-SNANA codes may allow for
# more SNe per version.
#
# The output versions are
#  SPLIT01_[VERSION]
#  SPLIT02_[VERSION]
#  ...
#  SPLITNN_{VERSION]   (NN = NSPLIT argument)
#
# The data files are NOT copied; each split version
# contains only the auxilary files, and the list file
# names are ../[VERSION]/[name].
# Written to work for both data & SIM, 
# but initial version tested only on simulations.
#
#
# Usage:
#  split_version.pl version  Nsplit
#
#  split_version.pl version  Nsplit  PRIVATE_DATA_PATH <path>
#       (optional private data path)
#
# April 12, 2011: fix dumb bug that started at 2nd file instead of 1st.
#
# June 03, 2011: make symbolic link to original IGNORE file.
#                Previously the IGNORE-epochs were losts.
#
# Aug 03, 2011: use strict
#
# Aug 26, 2011:  organize into 'subs' and ABORT if split already exists.
#
# Apr 2013: add optional argumement PRIVATE_DATA_PATH
#
# Nov 11, 2014: fix to work with files in both lcmerge/ or lcmerge/VERSION.
#
# ==========================

use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use sntools ;
use strict ;

my  $SNDATA_ROOT  = $ENV{'SNDATA_ROOT'};
my  $SNANA_DIR    = $ENV{'SNANA_DIR'};
my  $TOPDIR_SPLIT = "$SNDATA_ROOT/lcmerge" ;

my ($NARG, $VERSION, $NSPLIT );
my ($SIMFLAG, $TOPDIR_VERSION, $sDir_version, $PRIVATE_DATA_PATH  );
my (@SNfiles, $Nfiles,$NperSplit, $MXSNLC );


sub initStuff ;
sub parse_args ;
sub Data_or_SIM ;
sub parseListFile ;
sub check_MXSNLC ;
sub make_split_version ;

# ------------------- BEGIN MAIN -------------------

&initStuff();

# parse command-line args
&parse_args();

# check if this VERSION is in the lcmerge(data)  or SIM directory
&Data_or_SIM();


# =======================================
# parse List file to get total number of files
&parseListFile();


# check with MXSNLC in snana.car to make sure  that snana jobs will run
&check_MXSNLC();

# =======
# create new version for each split
my ($isplit) ;
for ( $isplit = 1; $isplit <= $NSPLIT; $isplit++ ) {
    &make_split_version($isplit);
} 


# ========= END OF MAIN ==========


sub initStuff {

    $PRIVATE_DATA_PATH = "" ;
    
} # end of initStuff

# ==================

sub parse_args {

    $NARG = scalar(@ARGV); 
    if ( $NARG < 2 ) {
	print " Must give VERSION and Nsplit as arguments. \n" ;
	die " ***** BAD ARGS for split_version.pl ****** \n" ;
    }

    $VERSION = $ARGV[0] ;
    $NSPLIT  = $ARGV[1] ;

    my ($i) ;
    for($i=2; $i <= $NARG; $i++ ) {
	if ( $ARGV[$i] eq "PRIVATE_DATA_PATH" )
	{ $PRIVATE_DATA_PATH = $ARGV[$i+1]; }
    }

    print " Split $VERSION into $NSPLIT sub-versions \n" ;

    if ( length($PRIVATE_DATA_PATH) > 0 ) {
	print " PRIVATE_DATA_PATH = $PRIVATE_DATA_PATH \n";
    }
    
} # end of parse_args


# ==================
sub Data_or_SIM {

    # Nov 11 2014: fix to allow data in lcmerge or in lcmerge/$VERSION

    my ($tmpData1, $tmpData2, $tmpSIM);


    if ( length($PRIVATE_DATA_PATH) > 0 ) {
	$TOPDIR_VERSION = $PRIVATE_DATA_PATH ;
	$sDir_version  = "" ;  
	$SIMFLAG        = 0 ; 
	return ; 
    }

    $tmpData1 = "$SNDATA_ROOT/lcmerge/${VERSION}.README" ;
    $tmpData2 = "$SNDATA_ROOT/lcmerge/${VERSION}/${VERSION}.README" ;
    $tmpSIM   = "$SNDATA_ROOT/SIM/${VERSION}/${VERSION}.README" ;

    print " xxx tmpData2 = $tmpData2 \n";

    $SIMFLAG = 0;  # 0=data, 1=SIM
    if ( -e $tmpData1 ) {
	print " Found $VERSION in lcmerge (data) area \n";
	$TOPDIR_VERSION   = "$SNDATA_ROOT/lcmerge";
	$sDir_version     = "../" ;  # relative path for list file
    }
    elsif ( -e $tmpData2 ) {
	print " Found $VERSION in lcmerge (data) area \n";
	$TOPDIR_VERSION  = "$SNDATA_ROOT/lcmerge/$VERSION" ;
	$sDir_version    = "../${VERSION}/" ;
    }
    elsif ( -e $tmpSIM ) {
	print " Found $VERSION in SIM area \n";
	$TOPDIR_VERSION   = "$SNDATA_ROOT/SIM/${VERSION}" ;
	$sDir_version     = "../${VERSION}/" ;
	$SIMFLAG          = 1;
    }
    else {
	die " Could not find version $VERSION \n";
    }

} # end of Data_or_SIM


# ======================
sub parseListFile {

    my ($listFile);

    $listFile = "${TOPDIR_VERSION}/${VERSION}.LIST" ;
    @SNfiles  = `cat $listFile` ;
    $Nfiles   = scalar(@SNfiles) ;
    print " Found $Nfiles SN files for version $VERSION \n";
    
    $NperSplit = int($Nfiles / $NSPLIT) ;
    print " Each sub-version will have about $NperSplit SNe \n" ;

} # end of parseListFile


# =====================
sub check_MXSNLC {

    my ($snanaCode, $tmp, @line, $hashLine, $MM );

    $snanaCode = "$SNANA_DIR/src/snana.car" ;
    $tmp  = qx(grep ",MXSNLC " $snanaCode) ;
    @line = split(/\s+/,$tmp) ;
    $MXSNLC = $line[4] ;
    
    if ( $NperSplit > $MXSNLC ) {
	$hashLine = '#######################################################' ;
	$MM       = "MXSNLC=$MXSNLC" ;
	print "\n $hashLine \n" ;
	print "\n    WARNING: each split version exceeds $MM \n";
	print "\n $hashLine \n" ;
	print "\n" ;
    }
    
} # end of check_MXSNLC


# ===========================
sub make_split_version {

    my ($isplit) = @_ ;

    my ($csplit, $ifilemin, $ifilemax, $ifile, $file, $ctmp ) ;
    my ($cmd1, $cmd2, $cmd3, $cmd4 );
    my ($splitListFile, $VERSPLIT, @MSGERR, $splitDir );

    $csplit     = sprintf("%2.2d", $isplit);
    $VERSPLIT   = "SPLIT${csplit}_${VERSION}" ;

# determine splitDir based on data or SIM
    if ( $SIMFLAG )  { $splitDir = "$SNDATA_ROOT/SIM/${VERSPLIT}" ; }
    else             { $splitDir = "$SNDATA_ROOT/lcmerge/${VERSPLIT}" ; }

    # abort if any SPLIT LIST-FILE already exists.
    $splitListFile   = "${splitDir}/${VERSPLIT}.LIST" ;
    if ( -e $splitListFile ) {
	print "\n" ;
	$MSGERR[0] = "Version $VERSION is already split." ;
	$MSGERR[1] = "Must remove SPLIT*_${VERSION}* to re-make split." ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    # create new directory for DATA or SIM
    qx(mkdir $splitDir);

    $ifilemax = $isplit * $NperSplit ;
    $ifilemin = $ifilemax + 1 - $NperSplit ;
    if ( $isplit == $NSPLIT ) {
	$ifilemax = $Nfiles ;
    }

    $ctmp = "ifile=$ifilemin to $ifilemax" ;
    print "  Prepare split version '$VERSPLIT'  ($ctmp) \n";


# create auxilary files
    $cmd1 = "cd $splitDir" ;
    $cmd2 = "ln -s ${sDir_version}${VERSION}.IGNORE  ${VERSPLIT}.IGNORE" ;
    $cmd3 = "ln -s ${sDir_version}${VERSION}.README  ${VERSPLIT}.README" ;
    $cmd4 = "touch ${VERSPLIT}.LIST" ;

    system("$cmd1 ; $cmd2 ; $cmd3 ; $cmd4 ");
  
# now fill up list file

    open LISTFILE , ">> $splitListFile" ;
    $ifile = $ifilemin - 1 ;
    while ( $ifile < $ifilemax ) {
	$file = $SNfiles[$ifile] ;
	print LISTFILE "${sDir_version}$file" ;
	$ifile++ ;
    } # ifile
    close LISTFILE ;

} # end of make_split_version
