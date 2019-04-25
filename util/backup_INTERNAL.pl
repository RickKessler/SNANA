#!/usr/bin/perl
#
# Created Aug 26, 2012 by R.Kessler
#
# Backup $SDNATA_ROOT/INTERNAL/$SURVEY directories and store
# it in  $SNDATA_BACKUP
# 
# This script is independent of backup_SNDATA_ROOT.pl that
# backs up the PUBLIC $SNDATA_ROOT and puts it in $SNANA_PUBLIC.
#
# These backups are for internal-collaborator access and
# are NOT meant as robust backups for disk failures !!!
#
# Usage:
# backup_INTERNAL.pl  ALL    (do them all)
# backup_INTERNAL.pl  SDSS
# backup_INTERNAL.pl  DES
# backup_INTERNAL.pl  DES -encp <enstorePath>
# backup_INTERNAL.pl  LSST
#
#
#       HISTORY
# ~~~~~~~~~~~~~~~~~~~~~~
# Nov 4, 2012: copy tarball to $SNDATA_BACKUP
#              print backup locations to screen (for disk and via URL)
#
# Aug 18m 2013: exclude contents of EXCLUDE_FROM_TAR
#
# Sep 25, 2013: optional copy to enstore with -encp arg ... NOT WORKING
#
# Feb 22 2017: 
#   + remove copy to URL; only copy to $SNDATA_BACKUP
#   + abort if EXCLUDE_FROM_TAR file is missing
#
# ------------------------------------------

# =========== Declarations =================

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;


# declare subs

sub miscInit ;
sub parse_args ;
sub do_backup ;

# globals

my $URL_TOP = "http://das.sdss2.org/ge/sample/sdsssn/SNANA-PUBLIC" ;
#my $URL_TOP = "http://sdssdp62.fnal.gov/sdsssn/SNANA-PUBLIC" ;
my $SDIR_INTERNAL = "INTERNAL" ;
my $q = "\'" ;

my ($SNDATA_ROOT);
my ($CDATE, $SNDATA_BACKUP, $URL ) ;
my ($ENCP_DIR, $ENCP_PATH, $ENCP_FLAG);
my ($SDIR, $INPUT_SDIR, @SDIR_LIST );
my (@MSGERR, $PATH_INTERNAL, $cdpath, $DOALL_FLAG, @BACKUP_SDIR_LIST );
my ($exclude_file, $EXCLUDE_FILE, @EXCLUDE_CONTENTS );

# =============== BEGIN MAIN =================

# hard-wired initializations
&miscInit();

# parse command-line args
&parse_args();


foreach $SDIR ( @BACKUP_SDIR_LIST ) { do_backup($SDIR); }

# =========== END MAIN ============



# ====================
sub miscInit {

    my (@tmp, $i, $sdir, $LL);

    $SNDATA_ROOT   = $ENV{'SNDATA_ROOT'}; 
    $SNDATA_BACKUP = $ENV{'SNDATA_BACKUP'}; 
    $ENCP_DIR      = $ENV{'ENCP_DIR'}; 
    $PATH_INTERNAL = "$SNDATA_ROOT/$SDIR_INTERNAL";
    $cdpath        = "cd $PATH_INTERNAL" ;

    $CDATE         = `date +%y%m%d_%H%M` ;
    $CDATE  =~ s/\s+$// ;   # trim trailing whitespace


    $DOALL_FLAG = 0;

    # get list of valid SDIRs for ALL option or for help.
    @tmp = qx($cdpath ; ls -d */); 

    # store in global @SDIR_LIST after trimming trailing whitespace 
    # and removing / at the end.
    $i = 0;
    foreach $sdir (@tmp) {
	$sdir  =~ s/\s+$// ;   
	$LL = length($sdir);
	$SDIR_LIST[$i] = substr($sdir,0,$LL-1);
	$i++ ;
    }

    $exclude_file       = "EXCLUDE_FROM_TAR" ;
    
    $ENCP_PATH = "" ;
    $ENCP_FLAG = 0 ;

} # end of miscInit


# ====================
sub parse_args {

    my ($NARG, $tmpDir, $sdir, $i );

    $NARG = scalar(@ARGV);

    if ( $NARG <= 0 ) { 
	$MSGERR[0] = "Must give name of subdir to backup.";
	$MSGERR[1] = "Check \$SNDATA_ROOT/INTERNAL ";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $INPUT_SDIR = $ARGV[0] ;
    if ( $INPUT_SDIR eq "ALL" ) { $DOALL_FLAG = 1; }
    if ( $INPUT_SDIR eq "all" ) { $DOALL_FLAG = 1; }
    
    for($i=1; $i < $NARG; $i++ ) {

	if ( $ARGV[$i] eq "-encp" ) { 
	    $ENCP_PATH = $ARGV[$i+1] ; 
	    $ENCP_FLAG = 1;
	    if ( $ENCP_DIR eq "" ) {
		$MSGERR[0] = "Env var ENCP_DIR not defined.";
		$MSGERR[1] = "Must do `setup encp' ";
		sntools::FATAL_ERROR(@MSGERR);
	    }
	}
    }  # end loop over NARG

    if ( $DOALL_FLAG == 0 ) {
    # make sure that $SNDATA_ROOT/INTERNAL exists
	$tmpDir = "$SNDATA_ROOT/$SDIR_INTERNAL/$INPUT_SDIR" ;
	if ( !(-d $tmpDir ) ) {
	    $MSGERR[0] = "Cannot find directory $tmpDir" ;
	    $MSGERR[1] = "  $tmpDir ";
	    $MSGERR[2] = "Available subdir options: ";

	    $i=0;
	    foreach $sdir (@SDIR_LIST ) {
		$i++ ;
		$MSGERR[2+$i] = "\t $sdir";
	    }	    
	    sntools::FATAL_ERROR(@MSGERR);
	}
	$BACKUP_SDIR_LIST[0] = $INPUT_SDIR
    }
    else {
	@BACKUP_SDIR_LIST = @SDIR_LIST ;
    }

} # end of parse_args


# ====================
sub do_backup {

    my ($sdir_in) = @_ ;

    my ($sdir_out,  $tarFile, $tarFile_gz, $tarCmd, $mv, $OUTDIR );

    
    $sdir_out = "${sdir_in}_${CDATE}" ;
    print " Backup '$PATH_INTERNAL/$sdir_in'  \n\t -> $sdir_out \n";

    $EXCLUDE_FILE = "$SNDATA_ROOT/INTERNAL/$sdir_in/$exclude_file" ;
    unless ( -e $EXCLUDE_FILE ) {
	die "\n FATAL: missing exclude file:\n  $EXCLUDE_FILE \n";
    }
    @EXCLUDE_CONTENTS = `cat $EXCLUDE_FILE` ;
    print " Files excluded from tarball: \n   ";
    sntools::printArray(@EXCLUDE_CONTENTS);
    print "\n\n";

    $tarFile    = "${sdir_out}.tar" ;
    $tarFile_gz = "${tarFile}.gz" ;

    $tarCmd = "tar -cf ${tarFile} $sdir_in " .
	"-X $sdir_in/$exclude_file ; gzip $tarFile " ;
    $mv     = "mv ${tarFile_gz} $SNDATA_BACKUP" ;

    print " Creating ${tarFile_gz}  ... \n";
    qx($cdpath ; $tarCmd ; $mv );

    # make sure everyone has read+write priv
    qx(cd $SNDATA_BACKUP ; chmod g+rw ${tarFile_gz});

    # ----------
    print " Backups located at \n";
    print "  \$SNDATA_BACKUP/${tarFile_gz} \n" ;

    #    print "mv = $mv \n";

} # end of do_backup

