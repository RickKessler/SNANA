#!/usr/bin/perl
#
# Created Mar 25, 2012
# For development node only, do a few things to prepare
# the start of the next snana version;
#
# 1. start new section in README_UPDATES
# 2. update SNANA_DIR_NAME in snana.car
# 3. update 'old' version in SNANA_TESTS.LIST
#
# Usage 
#   startSNANAversion.pl v9_94
#
#
# Dec 10 2017: replace snana.car with sntools.h to update snana version.
# May 04 2019: little cleanup to prepare for git.
#
# =======================================

use IO::Handle;
use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use sntools ;
use strict ;


my $SNDATA_ROOT   = $ENV{'SNDATA_ROOT'};
my $SNANA_DIR     = $ENV{'SNANA_DIR'};

# xxx mark delete ymy $SNANA_PUBLIC  = $ENV{'SNANA_PUBLIC'};
# xxx mark delete my $SNANA_DIR     = "/home/s1/rkessler/snana" ;

my $README_UPDATES_FILE = "$SNANA_DIR/doc/README_UPDATES" ;
my $snana_file          = "$SNANA_DIR/src/sntools.h" ;
my $SNANA_TESTS_FILE    = "$SNDATA_ROOT/SNANA_TESTS/SNANA_TESTS.LIST" ;

my $CDATE       = `date +%F` ; $CDATE   =~ s/\s+$// ;
my $qq = '"' ;

my ($SNANA_VERSION, @MSGERR);

my $DEBUG_MODE = 0 ;

sub checkVersion ;
sub update_README_UPDATES ;
sub update_snana ;


# =========== BEGIN MAIN ================

my $NARG = scalar(@ARGV) ;
if ( $NARG < 1 ) {
    $MSGERR[0] = "Must give new snana version as argument." ;
    sntools::FATAL_ERROR(@MSGERR);
}

$SNANA_VERSION = $ARGV[0] ;

print "\n =========== Prepare to start snana $SNANA_VERSION =========== \n";

if ( $DEBUG_MODE ) {
    print "\t !!!! DEBUG MODE => print but do nothing !!!! \n";
}


# make sure that this version has not already been used somewhere
&checkVersion();

&update_README_UPDATES();

&update_snana();

# &update_SNANA_TESTS() ;  # removed Aug 27 2014

# ============ END MAIN =============

# -------------------------------
sub checkVersion {
    
    # abort if this version is already used somewhere

    my (@bla, $V, $file);
    my @CHECK_FILES = ( $README_UPDATES_FILE , $snana_file  );

    $V    = "$SNANA_VERSION" ;

    foreach $file ( @CHECK_FILES ) {
	print "\t check '$file' for '$V'\n";
	@bla = qx(grep $V $file) ;
	if ( scalar(@bla) > 0 ) {
	    $MSGERR[0] = "SNANA version '${V}' is already listed in";
	    $MSGERR[1] = "$file .";
	    $MSGERR[2] = "Pick new version or remove $V from above file." ;
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

} # end of checkVersion


sub update_README_UPDATES {

    my $file = $README_UPDATES_FILE ;
    my $V    = "$SNANA_VERSION";   
    my $tab  = "        " ;  # avoid read tabs to avoid parsing headahces

    print " ------------------------------------------- \n";
    print " Add $SNANA_VERSION to \n\t $file \n";

    if ( $DEBUG_MODE ) { return ; }

    # -------------------------------------------------------
    # Aug 27 2014: remove the END statement so that it can 
    #              be re-written at the new end of file.
    my $sedcmd  = "sed -e '/END:/d'" ;
    my $tmpFile = "${file}_noEND" ;
    qx($sedcmd $file > $tmpFile ; mv $tmpFile $file);
    # --------------------------------------------------------

    open PTR_UPD , ">> $file" ;

    print PTR_UPD "\n";
    print PTR_UPD 
	"# ======================================================================= \n";

    print PTR_UPD "$tab $V (  ) \n\n";
    print PTR_UPD "$tab ***** IMPORTANT($V) ***** \n\n";
    print PTR_UPD "$tab ***** USEFUL($V) ***** \n\n";
    print PTR_UPD "$tab ***** MISCELLANEOUS($V) ***** \n\n";
    print PTR_UPD "END: \n";

    close PTR_UPD ;

}  # update_README_UPDATES {


sub update_snana {

    my (@bla, $oldLine, $newLine, $V, $sedcmd);
    my $STARTKEY = "BEGIN SNANA PARAMS";
    my $TMPFILE = "TEMP_sntools.h";

    print " ------------------------------------------- \n";


    @bla = qx(grep $qq SNANA_VERSION_CURRENT $qq $snana_file | grep define );
    $oldLine   = "$bla[0]" ;
    $oldLine   =~ s/\s+$// ;

#    $newLine = "&   SNANA_VERSION = \"$V\"  ! Begin $CDATE";

    $V   = $SNANA_VERSION ;
    $newLine = "#define  SNANA_VERSION_CURRENT  \"$V\" ";
    print " Replace \n $oldLine\n with \n $newLine \n";

    $sedcmd  = "sed " ;
    $sedcmd  = "$sedcmd" . "-e 's/$oldLine/$newLine/' " ;

#    $sedcmd  = "$sedcmd" . "-e '/ SNANA_VERSION_CURRENT /d' ";
#    $sedcmd  = "$sedcmd" . "-e '/$STARTKEY/a\\     $newLine' " ;

    $sedcmd  = "$sedcmd" . " $snana_file > $TMPFILE";

    if ( $DEBUG_MODE == 0 ) {
	qx($sedcmd ; mv $TMPFILE $snana_file);
    }
    else {
	print "\t sedcmd = $sedcmd \n";
    }

}

