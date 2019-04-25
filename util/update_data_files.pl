#!/usr/bin/perl
#
# Created Jan 8 2018
#
# Usage:
#   update_data_files.pl <VERSION>  <updFile>  '<varList>'
#   update_data_files.pl <VERSION>  <updFile>  '<varList>' CLOBBER
# where
#   <VERSION> is the data file version; must be TEXT format.
#             Give full path if not in SNDATA_ROOT/lcmerge.
#
#   <updFile>  fitres-formatted update-file with variables  & values.
#
#   <varList>  space-separated list of variables to update (in quotes).
#              e.g.,  'VPEC VPEC_ERR'
#
#   CLOBBER  is optional argument to clobber original files
#            instead of making new directory. Be careful to 
#            make backup before doing this.
#
# If variable exists in data file, then replace it;
# If variable does not exist, then add it into header
#

my ($INPUT_VERSION, $INPUT_UPDATE_FILE, $INPUT_VARLIST );
my $INPUT_CLOBBER = 0 ;

my $SNDATA_ROOT       = $ENV{'SNDATA_ROOT'};
my $OPT_ABORT = 0 ; # parsing option  
my $OPT_WARN  = 1 ; # 1=> leave warning message; 2=>quiet  
my $OPT_QUIET = 2 ; # no warnong if key not found 

my $NUPD = 0 ;

my ($DATA_DIR, $VERSION, @CONTENTS_UPDATE_FILE, $UPDFILE_VARNAMES);
my (@VARLIST, @INDX_VARLIST, $NVARLIST );
my (@DATAFILE_LIST, $OUTDIR);

# ========================
sub parse_args ;
sub prep_INPUTS ;
sub updateDataFile_driver ;
sub makeDataFile ;

use IO::Handle ;
use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# ========== BEGIN MAIN ==============

&parse_args();
&prep_INPUTS();

print "\n" ;
my ($file);
foreach $file (@DATAFILE_LIST) {
    $file      =~ s/\s+$// ;   # trim trailing whitespace 
    &updateDataFile_driver($file);
}

if ( $INPUT_CLOBBER ) { qx(rm -r $OUTDIR); }

print "\n Done updating $NUPD data files. \n";

# ========== END MAIN ============

sub parse_args {

    # use echo to allow for ENV
    $INPUT_VERSION      = qx(echo $ARGV[0]);
    $INPUT_UPDATE_FILE  = qx(echo $ARGV[1]);
    $INPUT_VARLIST      = $ARGV[2] ; 

    if ( scalar(@ARGV) > 3 ) {
	if ( $ARGV[3] eq "CLOBBER" )  { $INPUT_CLOBBER = 1 ; }
    }


    $INPUT_VERSION      =~ s/\s+$// ;   # trim trailing whitespace 
    $INPUT_UPDATE_FILE  =~ s/\s+$// ;   # trim trailing whitespace             

    print " VERSION     = $INPUT_VERSION \n";
    print " UPDATE_FILE = $INPUT_UPDATE_FILE \n";
    print " VARLIST     = '$INPUT_VARLIST' \n";
    print " CLOBBER     = $INPUT_CLOBBER  \n";
    return ;

} # end parse_args


sub prep_INPUTS {
    
    # - - - - - - - - - - - - - - - - - - - - - - - - 
    # parse INPUT_VERSOIN to get DATA_DIR and VERSION

    my $jslash = rindex($INPUT_VERSION,"/");  # location of last slash
    if ( $jslash <= 0 ) {
	# VERSOIN is in nominal location for data files
	$VERSION  = $INPUT_VERSION ;
	$DATA_DIR = "$SNDATA_ROOT/lcmerge/$VERSION" ;
    }
    else {
	# private location for version
	$DATA_DIR = $INPUT_VERSION ;
	$VERSION  = substr($DATA_DIR,$jslash+1,99);
    }

    print "\n";
    print " DATA_DIR = '$DATA_DIR' \n" ;
    print " VERSION  = '$VERSION' \n" ;

    if ( $INPUT_CLOBBER ) 
    { $OUTDIR = "CLOBBER_${VERSION}" ; }
    else
    { $OUTDIR = "${VERSION}" ; }

    if ( -d $OUTDIR ) { qx(rm -r $OUTDIR); }
    qx(mkdir $OUTDIR) ;

    # - - - - - - - - - - - - - - - - - 
    # break VARLIST into separate strings
    @VARLIST  = split(/\s+/,$INPUT_VARLIST) ;

    # - - - - - - - - - - - - - - - - - - - - - - - - 
    # read update file
    @CONTENTS_UPDATE_FILE = ();
    sntools::loadArray_fromFile($INPUT_UPDATE_FILE, \@CONTENTS_UPDATE_FILE);

    my ($key, @tmp, $varName, @INDX, $indx, @VARNAMES );

    $key     = "VARNAMES:" ;
    @tmp     = sntools::parse_array($key, 99, $OPT_QUIET, 
				    @CONTENTS_UPDATE_FILE);
    $UPDFILE_VARNAMES = "$tmp[0]" ;
    print " --> found VARNAMES = '$UPDFILE_VARNAMES' \n";
    @VARNAMES  = split(/\s+/,$UPDFILE_VARNAMES) ;

    # locate each element of VARLIST in UPDFILE_VARNAMES
    $NVARLIST = 0 ;
    foreach $varName ( @VARLIST ) {
	@INDX =
	    grep { $VARNAMES[$_] eq "$varName" } 0 .. $#VARNAMES;
	$indx = $INDX[0] ;
	print "\t '$varName' index is $indx \n";

	$INDX_VARLIST[$NVARLIST] = $indx;
	$NVARLIST++ ;


    }


    # - - - - - - - - - -
    # read VERSION-LIST file 

    my $LISTFILE = "${DATA_DIR}/${VERSION}.LIST" ;
    open  PTR_TEMP , "< $LISTFILE" ;
    @DATAFILE_LIST = <PTR_TEMP> ;
    close PTR_TEMP ;

    return ;  

} # end prep_INPUTS


sub updateDataFile_driver {
    my ($file) = @_ ;

    my ($DATAFILE, $OUTFILE, $key, @tmp, $SNID, $VERB );
    my (@wdlist, $ivar, $indx, $varName, $VAL, $sedCmd, $updLine );

    $DATAFILE = "$DATA_DIR/$file" ;
    $OUTFILE  = "$OUTDIR/$file" ;

    $key  = "SNID:" ;
    @tmp  = sntools::parse_line($DATAFILE, 1, $key, $OPT_QUIET ) ;
    $SNID = $tmp[0] ;

    $key     = "$SNID" ;
    @tmp     = sntools::parse_array($key, 99, $OPT_QUIET, 
				    @CONTENTS_UPDATE_FILE);
    if ( scalar(@tmp) > 0 ) 
    { $VERB = "Updated"; @wdlist  = split(/\s+/,$tmp[0] ) ; }
    else
    { $VERB = "Skipped"; goto UPDATE_DONE ; }


    $sedCmd  = "sed  " . 
	"-e '/FILTERS:/a\#' " .
	"-e '/FILTERS:/a\# Variables copied from' " .
	"-e '/FILTERS:/a\#   $INPUT_UPDATE_FILE' " ;
    
    for ( $ivar=0; $ivar < $NVARLIST; $ivar++ ) {
	$varName = $VARLIST[$ivar] ;
	$indx    = $INDX_VARLIST[$ivar] - 1 ;
	$VAL     = $wdlist[$indx] ;
	$updLine = "${varName}:  $VAL";

	# check if this varName exists already in the data file
	$key     = "${varName}:" ;
	@tmp     = sntools::parse_array($key, 99, $OPT_QUIET, 
					@CONTENTS_UPDATE_FILE);
	if ( scalar(@tmp) > 0 ) {
	    $sedCmd = "$sedCmd " . "-e 's/$tmp[0]/$updLine/g' "
	}
	else {
	    $sedCmd = "$sedCmd " . "-e '/FILTERS:/a\ $updLine' " ;
	}
    }

    $sedCmd  = "$sedCmd  -e '/FILTERS:/a\#' " ;

    qx($sedCmd $DATAFILE > $OUTFILE);

    # check option to clobber original file with udpated file.
    if ( $INPUT_CLOBBER ) {
	qx(mv $OUTFILE $DATAFILE);
    }

#    die "\n\n xxx DEBUG DIE $sedCmd \n";

  UPDATE_DONE:
    $NUPD++ ;

    my $c3N = sprintf("%3.3d", $NUPD);
    print "   $c3N : $VERB $file  ($SNID) \n";

    return ;

} # end  updateDataFile


sub  makeDataFile {

    my ($DATAFILE, $varName, $VAL) = @_ ;

    # Created new data file.
    # If varName exists, replace it.
    # Otherwise, add it
    #

} # end makeDataFile
