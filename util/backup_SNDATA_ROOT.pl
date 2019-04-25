#!/usr/bin/perl
#
# Created Aug 2011 by R.Kessler
# [re-make backup_SNSATA_ROOT.cmd -> perl]
#
# Create a tarball of $SNDATA_ROOT
# and then copy the gzipped tarball to 
#    $SNDATA_BACKUP  
#    $SNANA_PUBLIC_LOCAL/downloads  ! staged public area, but still local
#    $SNANA_PUBLIC_URL/downloads    ! true public area via scp
#
# Usage:
#   backup_SNDATA_ROOT.pl            
#      (PUBLIC backup and copy to public URL)
#
#   backup_SNDATA_ROOT.pl INTERNAL   
#    (everything including $SNDATA_ROOT/INTERNAL and $SNDATA_ROOT/lcmerge;
#     backup only to local disks on dp62; nothing goes to URL ! )
#
#
#
#       HISTORY
# ~~~~~~~~~~~~~~~~~~~~~~
# Aug 17, 2012
#  Replace SURVEY argument with INTERNAL to backup 
#  $SNDATA_ROOT/INTERNAL*  
#
# Sep 2013: replace /analysis with /sample_input_files
#
# Jul 25 2014: 
#   fix to backup the *PHOT.FITS file for fits-formatted versions.
#
# Jul 28 2014: scp to new URL host on sdsslnx34 .
#       See $SNANA_PUBLIC_URL
#
# Nov 10 2014: fix checkDataVersions() to check version-dir as well
#              as list-file.
#
# Jun 8 2015:
#    include models listed under SIMSED/AAA_PUBLIC_VERSIONS.README
#    [based on interest at Aspen 2015 transient workshop]
#
# jan 2 2018: allow saving tar or tar.gz file in /lcmerge (see -e $VDIR)
#
# ---------------------
#
# =========== Declarations =================
use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;


# declare subs
sub miscInit ;
sub parse_args ;
sub get_tarList ;
sub makeTarball ;
sub addData2Tarball ;
sub copyTarball ;
sub publicRelease ;
sub query_continue ;
sub checkDataVersions ;

# globals
my $SDIR_INTERNAL = "INTERNAL" ;
my $DEVEL_HOST    = "fnal";

my $q = "\'" ;
my ($SNDATA_ROOT, $SNANA_PUBLIC_LOCAL, $SNANA_PUBLIC_URL);
my ($SNDATA_BACKUP, $CDATE );
my ($PUBLIC_FLAG );
my (@tarList, $tarFile, $TARFILE, $NOBJTAR, $SUFFIX);
my (@DATA_VERSIONS, @SIMSED_VERSIONS );
my ($EXCLUDE_FILE, @EXCLUDE_CONTENTS, @MSGERR );
my ($downloadDir, $cdroot, $DATA_DIR, $SIMSED_DIR, $PUBLIC_DATA_REAMDE );
my $TESTMODE = 0 ;

# =============== BEGIN MAIN =================

# hard-wired initializations
&miscInit();

# parse command-line args
&parse_args();

# get list of files and directories to archive
&get_tarList();

# check that all data versions exist before making tarball.
&checkDataVersions();

# ask to continue
&query_continue();

# make the tarball
&makeTarball();

# copy tarball to backup directory (env arg is defined)
&copyTarball("SNDATA_BACKUP"); 

# public release (from FNAL only)
&publicRelease();


# =========== END MAIN ============



# ====================
sub miscInit {

    $SNANA_PUBLIC_LOCAL  = $ENV{'SNANA_PUBLIC_LOCAL'} ;
    $SNANA_PUBLIC_URL    = $ENV{'SNANA_PUBLIC_URL'} ;
    $SNDATA_ROOT   = $ENV{'SNDATA_ROOT'}; 
    $cdroot        = "cd $SNDATA_ROOT" ;
    $DATA_DIR      = "$SNDATA_ROOT/lcmerge";
    $SIMSED_DIR    = "$SNDATA_ROOT/models/SIMSED";

    $EXCLUDE_FILE       = "EXCLUDE_FROM_TAR" ;

    $downloadDir        = "downloads" ;
    $PUBLIC_DATA_REAMDE = "AAA_PUBLIC_VERSIONS.README" ;
 
    $PUBLIC_FLAG   = 1 ; # default is for public release

    $CDATE  = `date +%F` ;
    $CDATE  =~ s/\s+$// ;   # trim trailing whitespace

    if ( $TESTMODE ) {
	$downloadDir = "downloads_test" ;
	if ( !(-d "$SNANA_PUBLIC_LOCAL/$downloadDir" ) ) 
	{ qx( mkdir -p $SNANA_PUBLIC_LOCAL/$downloadDir/archive ); }
    }

} # end of miscInit


# ====================
sub parse_args {

    my ($NARG, $tmpDir);

    $NARG = scalar(@ARGV);

    if ( $NARG <= 0 ) { return ; }


    if ( $ARGV[0] eq "INTERNAL" ) 	{  
	$PUBLIC_FLAG = 0 ; 
	# make sure that $SNDATA_ROOT/INTERNAL exists
	$tmpDir = "$SNDATA_ROOT/$SDIR_INTERNAL" ;
	if ( !(-d $tmpDir ) ) {
	    $MSGERR[0] = "\$SNDATA_ROOT/$SDIR_INTERNAL does not exist ?!?!?";
	    $MSGERR[1] = "Either create above directory or change argument.";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

} # end of parse_args


# ====================
sub get_tarList {

    my ($tmpFile);

    # get list of files and subdirs; then exlude
    # anything beginning with SIM

    @EXCLUDE_CONTENTS = `cat $SNDATA_ROOT/$EXCLUDE_FILE` ;

    $NOBJTAR = 0 ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "-X $EXCLUDE_FILE" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "SIM/README" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "AAA_README" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "standards" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "filters" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "simlib" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "snsed" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "sample_input_files" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "SURVEY.DEF" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "SNANA_TESTS" ;
    $NOBJTAR++ ;  $tarList[$NOBJTAR] = "INTERNAL" ;  # Aug 24 2012

    # include the big directories below unless in TESTMODE
    if ( $TESTMODE == 0 ) {
	$NOBJTAR++ ;  $tarList[$NOBJTAR] = "kcor" ;
	$NOBJTAR++ ;  $tarList[$NOBJTAR] = "models" ;
	$NOBJTAR++ ;  $tarList[$NOBJTAR] = "MWDUST" ;
    }

# load up the public data & SIMSED versions
    $tmpFile = "$DATA_DIR/$PUBLIC_DATA_REAMDE";
    @DATA_VERSIONS = qx(cat $tmpFile) ;
    
    $tmpFile = "$SIMSED_DIR/$PUBLIC_DATA_REAMDE";
    @SIMSED_VERSIONS = qx(cat $tmpFile) ;


    if ( $PUBLIC_FLAG ) {
	print " Prepare PUBLIC backup. \n" ;
	$SUFFIX = "PUBLIC" ;

	# exclude proprietary contents of INTERNAL (Aug 24, 2012)
	$NOBJTAR++ ;  $tarList[$NOBJTAR] = "--exclude=${q}INTERNAL/*${q}" ; 
    }
    else {
	$SUFFIX = "$SDIR_INTERNAL" ;
	print " Prepare PUBLIC + $SUFFIX  backup. \n" ;	
    }


    $tarFile = "SNDATA_ROOT-${SUFFIX}-${CDATE}.tar";
    $TARFILE = "$SNDATA_ROOT/${tarFile}.gz" ;

    print "\n";
    print " Tarball list =  \n   " ;
    sntools::printArray(@tarList);

    print "\n\n";
    print " Include Data Versions:  \n   " ;   
    sntools::printArray(@DATA_VERSIONS);

    print "\n\n";
    print " Include SIMSED Versions:  \n   " ;   
    sntools::printArray(@SIMSED_VERSIONS);

    print "\n\n";
    print " Files excluded from tarball: \n   ";
    sntools::printArray(@EXCLUDE_CONTENTS);

    print "\n\n";
    print " output tarFile: $tarFile \n" ;
    print "\n";
} # end of get_tarList


# ==============================
sub query_continue {

    my ($response);

    if ( $TESTMODE ) {
	print "\t !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n" ;
	print "\t !!! BEWARE: YOU ARE IN TEST-MODE !!! \n" ;
	print "\t !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n" ;
	print "\n";
    }

    print " Proceed to make $SUFFIX tarball y/[n] ? => " ;
    $response = <STDIN> ;
    $response =~ s/\s+$// ;
    unless ( $response eq "y" ) {
	print " Bye bye. \n";
	die " ***** ABORT ***** \n" ;
    }

} # end of query_continue

# ==========================
sub  checkDataVersions {

    # for public data versions, make sure that each data version exists.
    # Nov 10 2014: check both list file and version-dir.

    my ($VERSION, $LIST_FILE, $VDIR, $ISgz );
    my ($EXIST_VERSION );

    if ( $PUBLIC_FLAG == 0 ) { return ; }

    foreach $VERSION ( @DATA_VERSIONS ) {
	$VERSION   =~ s/\s+$// ;   # trim trailing whitespace
	$LIST_FILE = "$DATA_DIR/${VERSION}.LIST" ;
	$VDIR      = "$DATA_DIR/${VERSION}" ;

	$EXIST_VERSION = 0 ;
	if ( (-e $LIST_FILE) || ( -d $VDIR ) || ( -e $VDIR) ) 
	{ $EXIST_VERSION=1; } 

	if ( $EXIST_VERSION == 0 ) {
	    $MSGERR[0] = "Data version '$VERSION' does not exist !?!?!";
	    $MSGERR[1] = "Check \$SNDATA_ROOT/lcmerge/$PUBLIC_DATA_REAMDE";
	    sntools::FATAL_ERROR(@MSGERR);
	}	
    }

}   # end of  checkDataVersions 


# ===========================
sub makeTarball {

    my ($VERSION) ;
    my $DEBUG = 0;

    if ( $DEBUG ) { @tarList = ( "-X $EXCLUDE_FILE" , "models" , "filters");}

    print "\n ------------------------------------------- \n";
    print " Begin creating $tarFile ... \n" ;
    qx($cdroot ; tar -cf $tarFile @tarList);

    # For public option, include specific set of public data versions.
    if ( $PUBLIC_FLAG  ) {

	print "\n";
	foreach $VERSION ( @DATA_VERSIONS ) {
	    $VERSION  =~ s/\s+$// ;   # trim trailing whitespace
	    if ( $DEBUG == 0 ) { &addData2Tarball($VERSION); }
	}

	print "\n";
	foreach $VERSION ( @SIMSED_VERSIONS ) {
	    $VERSION  =~ s/\s+$// ;   # trim trailing whitespace
	    &addSIMSED2Tarball($VERSION);
	}
    }

    print " gzipping $tarFile ... \n" ;
    qx($cdroot ; gzip $tarFile);

    if ( $DEBUG ) { die "\n xxxx DIE DEBUG xxxx \n"; }

} # end of makeTarball


# =========================
sub addData2Tarball {

    # add data $VERSION to tarball.
    # read LIST file to get list of data files.
    # 
    # Jul 25 2014: if FITS format, add PHOT.FITS file (see $file2)
    # Nov 10 2014: check for version folder or list file.

    my ($VERSION) = @_ ;

    my (@dataFiles, $dataFileList, $file, $file2, $cdd, $NDATA);
    my ($LIST_FILE, $VDIR, $lcmergeTar );

    $cdd = "cd $DATA_DIR";
    
    # read the LIST file and tar the data files.

    $VDIR      = "$DATA_DIR/${VERSION}";
    if ( -d $VDIR ) 
    {  $LIST_FILE = "$VDIR/${VERSION}.LIST";    }
    else 
    {  $LIST_FILE = "$DATA_DIR/${VERSION}.LIST";  }


    # --------------------------

    if ( -d $VDIR ) {
	# just get the version folder
	$lcmergeTar = "lcmerge/${VERSION}";
    }
    elsif ( -e $VDIR ) {
	# probably a tar or tar.gz file
	$lcmergeTar = "lcmerge/${VERSION}";
    }
    else {
	# get explicit files under /lcmerge
	@dataFiles = qx(cat $LIST_FILE);
	$NDATA     = scalar(@dataFiles);

	$dataFileList = "" ;
	foreach $file (@dataFiles) {
	    $file  =~ s/\s+$// ;   # trim trailing whitespace
	    $dataFileList = "$dataFileList lcmerge/$file";
	    
	    # if this is a HEAD.FITS file, include PHOT.FITS file too.
	    my $j = index($file,"HEAD.FITS") ;
	    if( $j > 0 ) {
		$NDATA++ ;
		$file2 = substr($file,0,$j) . "PHOT.FITS" ;
		$dataFileList = "$dataFileList lcmerge/$file2";
	    }
	}
	$lcmergeTar = "lcmerge/${VERSION}.*  $dataFileList";
    }

    print "\t  Add public data version $VERSION  \n";
    qx($cdroot ; tar -rf $tarFile  $lcmergeTar );     

} # end of addData2Tarball {



# =========================
sub addSIMSED2Tarball {

    # add SIMSED $VERSION to tarball.

    my ($VERSION) = @_ ;

    my ($TMPDIR, $tmpSDIR );
    # --------------------------

    $TMPDIR  = "$SIMSED_DIR/${VERSION}";
    $tmpSDIR = "models/SIMSED/${VERSION}";
    unless ( -d $TMPDIR ) {
	$MSGERR[0] = "Cannot find SIMSED version" ;
	$MSGERR[1] = "$TMPDIR";
	sntools::FATAL_ERROR(@MSGERR);	
    }

    print "\t  Add SIMSED version $VERSION  \n";
    qx($cdroot ; tar -rf $tarFile  $tmpSDIR );     

} # end of addSIMSEDTarball {


# ===========================================
sub copyTarball {

    # copy tarball to $envName
    my ($envName) = @_ ;
    my ($cpDir);

    # decode name of environment variable into cpDir.
    $cpDir = $ENV{$envName} ; 

    if ( $cpDir eq ' ' ) {
	print " \$$envName does not exist => skip copy.\n";
    }
    else {
	print "\n Copy SNDATA_ROOT tarball to \$$envName \n" ;
	qx(cp $TARFILE $cpDir);
    }

} # end of copyTarball 


# ========================
sub publicRelease {

    # copy tarball to directory that is accessible to the public. 
    # This part runs only if $SNANA_PUBLIC_URL is defined and if the 
    # host machine is the development node.
    # 
    # July 27 2014: release in two stages.
    #    1) to SNANA_PUBLIC_LOCAL   (stage for public release)
    #    2) to SNANA_PUBLIC_URL     (via new scp )

    my ($HOST, $DEVEL_FLAG, $cdd, $cpDir, $versionFile, $VERSION_FILE );
    my ($archiveFile, $publicFile, $cmd1, $cmd2, $OLD_CDATE);
    my ($URLDIR, $scp );

    $HOST = $ENV{'HOST'} ; 
    $DEVEL_FLAG = index($HOST,$DEVEL_HOST) ;

    if ( $PUBLIC_FLAG        == 0   ) { return ; }
    if ( $SNANA_PUBLIC_LOCAL eq ' ' ) { return ; }
    if ( $DEVEL_FLAG         <  0   ) { return ; }

    print "\n";
    print " Copy SNDATA_ROOT tarball to \$SNANA_PUBLIC_LOCAL = \n";
    print "    $SNANA_PUBLIC_LOCAL \n";

    $cpDir        = "$SNANA_PUBLIC_LOCAL/$downloadDir" ;
    $cdd          = "cd $cpDir";
    $publicFile   = "SNDATA_ROOT.tar.gz" ;
    $versionFile  = "VERSION.SNDATA_ROOT";
    $VERSION_FILE = "$cpDir/$versionFile" ;

    # archive current/old SNDATA_ROOT before updating new SNDATA_ROOT.
    if ( -e $VERSION_FILE ) {
	$OLD_CDATE   = `cat $VERSION_FILE` ;
	$OLD_CDATE   =~ s/\s+$// ;   # trim trailing whitespace
	$archiveFile = "archive/SNDATA_ROOT-${OLD_CDATE}.tar.gz" ;
	print "\t Archive previous backup to $archiveFile \n";
	qx($cdd ; mv $publicFile $archiveFile );
	qx(rm $VERSION_FILE);
    }
    # update date-stamp in version file (for php code)
    print "\t Update SNDATA_ROOT time-stamp to '$CDATE' \n";
    qx(echo $CDATE > $VERSION_FILE );

    #  copy the tarball into local directory  that mirrors
    #  the URL directory structure.
    print "\t Copy tarball to $publicFile \n";
    qx($cdd ; cp $TARFILE  $publicFile);

    # -------------------------------------------------------
    # Jul 28 2014: copy to the new URL on a separate linux box

    print "\n" ;
    print " ======= Update public SNDATA_ROOT ======= \n";
    print " Copy tarball to URL at\n   $SNANA_PUBLIC_URL \n";

    $URLDIR = "$SNANA_PUBLIC_URL/$downloadDir" ;
    $scp    = "scp $VERSION_FILE  $cpDir/$publicFile  $URLDIR";
    qx($scp) ;

    # after debugging
    # remove kcor/* and MWDUST/ from EXCLUDE_FROM_TAR
    # put back  lcmerge/AAA_PUBLIC_VERSIONS.READ

} # end of publicRelease
