#!/usr/bin/perl
#
# Created Aug 2011 by R.Kessler
#
# make tarball of /lcmerge/${VERSION}* to
# lcmerge/archive/${VERSION}_${cdate}.tar.gz
# where cdate = yyyy-nn-dd
#
# Usage
#   backup_SNDATA_VERSION.pl <VERSION>
#
#  where <VERSION> refers to a version of SNDATA files.
#
# Note that the archive backup is on the same disk
# as the original version, so this backup will NOT
# protect against disk failure ... only protects
# against stupid mistakes in the /lcmerge area.
#
# For disk-failure protection, see backup_SNDATA_ROOT.pl
# script and SNDATA_ROOT backups on $SNDATA_BACKUP.
#
#
# ------------------------------------------

# =========== Declarations =================
use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;


# declare subs
sub parse_args ;
sub miscInit ;
sub makeTarball ;

# globals

my ($SNDATA_ROOT, $DATA_DIR, $CDATE, $tarFile );
my (@SUFFIX_LIST, @MSGERR, $VERSION ) ;

# ========= BEGIN MAIN  =========

&parse_args();

&miscInit();

&makeTarball();

# ========= END OF MAIN ==========



# ======================
sub parse_args {

    my ($NARG);

    $NARG = scalar(@ARGV);
    if ( $NARG < 1 ) {
        $MSGERR[0]= " Must give data version as argument. " ;
        sntools::FATAL_ERROR(@MSGERR);
    }

    $VERSION = $ARGV[0];
    print " Backup data version '$VERSION' \n";

} # end of parse_args


# ===============================
sub miscInit {

    $SNDATA_ROOT = $ENV{'SNDATA_ROOT'};
    $DATA_DIR    = "$SNDATA_ROOT/lcmerge" ;
    $CDATE       = `date +%F` ;
    $CDATE       =~ s/\s+$// ;   # trim trailing whitespace
    $tarFile     = "${VERSION}_${CDATE}.tar";

    @SUFFIX_LIST = ( "README" , "LIST" , "IGNORE" );

} # end of miscInit

# ============================
sub makeTarball {

    my ($driverFileList, @dataFiles, $cmd1, $cmd2, $cmd3, $cdd );
    my ($suffix, $dataFileList, $file, $NDATA );

    print " Creating tarball \n";
    print "   \$SNDATA_ROOT/lcmerge/archive/${tarFile}.gz \n" ;

    $cdd = "cd $DATA_DIR";

    foreach $suffix ( @SUFFIX_LIST ) {
        $suffix  =~ s/\s+$// ;   # trim trailing whitespace
	$driverFileList = "$driverFileList ${VERSION}.${suffix}" ;
    }

    # convert list of @dataFiles into linear list that does not
    # have any <CR>
    @dataFiles = qx($cdd ; cat ${VERSION}.LIST);
    $NDATA = scalar(@dataFiles);
    $dataFileList = "" ;
    foreach $file (@dataFiles) {
        $file  =~ s/\s+$// ;   # trim trailing whitespace	
	$dataFileList = "$dataFileList $file";
    }

    $cmd1   = "tar -cf $tarFile ${driverFileList} ${dataFileList}" ;
    $cmd2   = "gzip $tarFile" ;
    $cmd3   = "mv ${tarFile}.gz archive" ;

#    print " cmd1 = $cmd1 \n" ; 
#    print " cmd2 = $cmd2 \n" ;
#    print " cmd3 = $cmd3 \n" ;

    qx($cdd ; $cmd1 ; $cmd2 ; $cmd3 );
    
    print " Done archiving $NDATA data files. \n";

}  # end of makeTarball
