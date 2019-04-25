#!/usr/bin/perl
#
# Created Dec 17 2014 by R.Kessler
#
# Set C pre-processor command for HBOOK,ROOT, etc ..
# based on required ENV. This is called by Makefile,
# but can also be ran interactively.
# The file $SNANA_DIR/src/sntools_output.h is modified
# if the current pre-processor flag(s) are not consistent
# with the required ENV variables.
#
# Nominal pre-processor flags in sntools_output.h are
#   #define USE_HBOOK  (if $CERN_DIR exists)
#   #define USE_ROOT   (if $ROOT_DIR exists)
#
# If the associated ENV does not exist then the sed utility
# is used to substitute
#  USE_HBOOK --> USE_HBOOKxxx  and/or USE_ROOT --> USE_ROOTxxx
#
# --------------------

use strict ;

##my $SNANA_DIR  = $ENV{'SNANA_DIR'};
my $incFile    = "sntools_output.h" ;
##my $INCFILE    = "$SNANA_DIR/src/$incFile";
my $INCFILE    = "../src/$incFile";

my @ENVLIST    = ( "CERN_DIR",    "ROOT_DIR" );
my @cFLAGLIST  = ( "USE_HBOOK" ,  "USE_ROOT" );
my $NFLAG      = scalar(@cFLAGLIST);

sub checkFlag ;
sub setFlag ;

# ============== BEGIN MAIN ============


my($iflag);
for ( $iflag=0; $iflag < $NFLAG; $iflag++ ) {
    &checkFlag($iflag);
}

# =========== END MAIN ============

sub checkFlag {
    my ($iflag) = @_ ;

    my $ENV_NAME   = $ENVLIST[$iflag] ;
    my $cFLAG      = $cFLAGLIST[$iflag] ;
    my $ENV_VAL    = $ENV{$ENV_NAME};

    my $SETIT = length($ENV_VAL);

    # print " $ENV_NAME  = '$ENV_VAL'  ($SETIT)\n";

    &setFlag($SETIT, $cFLAG);

}  # end of checkFlag


sub setFlag {
    my ($SETIT, $cFLAG) = @_ ;

    # SETIT = 1 --> insert '#define $cFLAG' in incFile
    # SETIT = 0 --> undefine cFLAG

    my (@bla, @wdlist, $flag, $MATCH, $oldLine, $newLine);

    @bla = qx(grep define $incFile | grep $cFLAG );
    $oldLine  = "$bla[0]" ;
    $oldLine  =~ s/\s+$// ;   # trim trailing whitespace
    @wdlist   = split(/\s+/,$oldLine) ;
    $flag     = $wdlist[1] ;
    $MATCH    = 0 ;
    if ($flag eq $cFLAG) {$MATCH=1;}

    # bail if flag is already set correctly
    if ( $SETIT    && $MATCH    ) { return ; }
    if ( $SETIT==0 && $MATCH==0 ) { return ; }

    if ( $SETIT ) { 
	print "\t Set $cFLAG pre-processor flag in $incFile \n";
	$newLine = "#define ${cFLAG}" ;
    }
    else { 
	print "\t UNset $cFLAG pre-processor flag in $incFile \n";
	$newLine = "#define ${cFLAG}xxx" ;
    }

    # use 'sed' to make the substitution
    my $incFile_tmp = "${incFile}_TEMP" ;
    my ($CMD);
    
    $CMD = "sed -e 's/$oldLine/$newLine/1' " ;
        
    qx($CMD $INCFILE > $incFile_tmp ; mv $incFile_tmp $INCFILE);

} # end of setFlag
