#!/usr/bin/perl
#
# Created April 2013 by R.Kessler
# Utility to examine and dump analysis info from any SNTABLE
# in hbook (ntuple) or root (tree) format.
# The basic functions here are
#  - print list of all available variables in SNTABLE
#  - dump subset of variables for each event into text file.
#
# This utility is aimed for those who prefer NOT to use 
# hbook or root, and prefer to extract analysis variables
# into text files.
#
#
# USAGE:
#   sntable_dump.pl <inFile>
#      (list tables in <inFile>: hbook or root file)
#
#   sntable_dump.pl <inFile>  <tableID>
#      (list variable names in table, then exit)
#
#   sntable_dump.pl <inFile>  <tableID>  --v <list of variables>
#      (dump values for each table entry)
#
#   sntable_dump.pl <inFile>  <tableID>  --v <list of variables> --o <outFile>
#     (overwrite default outfile in ASCII/TEXT format)
#
#   sntable_dump.pl <inFile>  <tableID>  --v <list of variables> NOHEADER
#     (out file has no header info and no 'SN:' key)
#
#   sntable_dump.pl <inFile>  <tableID>  --v <vlist>  --a <textFile>
#     (append fitres text file with vlist variables)
#
#   sntable_dump.pl <inFile>  <tableID>  --outlier <Nsig1> <Nsig2>
#    (dump epochs for data-fit outlier between Nsig1 and Nsig2 sigma)
#
#   sntable_dump.pl <inFile>  <tableID>  --outlier_sim <Nsig1> <Nsig2>
#    (dump epochs for data-sim outlier between Nsig1 and Nsig2 sigma)
#
#   sntable_dump.pl <inFile>  <tableID>  --outlier <Nsig1> <Nsig2> \
#          --format IGNORE
#    (dump outliers in IGNORE format for SNANA data version)
#
#   sntable_dump.pl <inFile>  <tableID>  OBS
#     (dump each observation, mainly for fluxerrmap)
#
#  Note that --[v,o,a] and -[V,O,A] are accpeted (both lower,upper case)
#
#
#      HISTORY
#
# Apr 26 2014: allow HBOOK/hbook extenstion along with HIS/.his
#
# Aug 2014: new argument  '--outlier <n1> <n2>'
#
# Oct 28, 2014: .txt -> .text for output.
#
# Jan 10 2015: use H or R argument for combine_fitres so that output
#              table format (root or hbook) matches input table format.
#
# May 10 2016: 
#   for HBOOK files, allow root-string table names SNANA & FITRES.
#
# Feb 22 2018: new  OBS argument to dump observations for fluxerrmap
#
# ------------------------------------------

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# ----------- global declarations ------------

my $JOBNAME = "sntable_dump.exe" ;
my $IFILETYPE_HBOOK    = 1 ; 
my $IFILETYPE_ROOT     = 2 ; 
my @COMINE_FITRES_ARG  = ( "?", "H", "R" );
my @VALID_SUFFIX_HBOOK = ( ".HIS", ".HBOOK", ".TUP"  ) ;
my @VALID_SUFFIX_ROOT  = ( ".ROOT" ) ;

# inputs
my ($TABLE_INFILE, $TABLE_ID, @TABLE_VARLIST, $NVARLIST);
my ($OUTFILE_TEXT, $OPT_HEAD, $TEXTFILE_APPEND);
my (@OUTLIER_NSIGMA, $KEY_OUTLIER, $FORMAT_OUTFILE, $OPTOBS );


# misc. globals
my (@MSGERR, $CMD_DUMP, $IFILETYPE );

# -----------  function declarations -----------
sub parse_args ;
sub insert_ccid_varlist ;
sub get_fileType ;
sub check_TABLE_ID_ROOT ;
sub check_TABLE_ID_HBOOK ;
sub make_dump_command ;
sub run_dump_command ;
sub textfile_append ;

# ================ BEGIN MAIN ===============

# parse command-line args
&parse_args() ;


# make sure CCID is 1st element of VARLIST
&insert_ccid_varlist();  


# determine root or hbook file type
&get_fileType()  

# check table ID
&check_TABLE_ID_ROOT();
&check_TABLE_ID_HBOOK();

# construct C-program command 
&make_dump_command();  


# run C program
if ( &run_dump_command() != 0 ) { die "" ; }


# check option to append fitres file with table-extracted variables
&textfile_append();


# ================ END MAIN =================

sub parse_args {

    my ($NARG,  $i, $VFLAG, $USE, $NBAD, $KEY );
    $NARG = scalar(@ARGV); 

    if ( $NARG < 1 ) {
	$MSGERR[0] = "Must give table-File as arguments";
	sntools::FATAL_ERROR(@MSGERR); 
    }

    # init defaults
    $TABLE_INFILE    = "NULL" ;
    $TABLE_ID        = "" ;
    @TABLE_VARLIST   = ( );
    $OUTFILE_TEXT    = "" ;
    $TEXTFILE_APPEND = "" ;
    $OPT_HEAD        = "" ; # default (blank) is to add header
    $VFLAG = 0 ;
    
    @OUTLIER_NSIGMA = ( -9, -9 );
    $KEY_OUTLIER    = "xxx" ;
    $OPTOBS         = 0 ;

    $FORMAT_OUTFILE = "DEFAULT" ; # default is key format

    $NBAD = 0;

    # ------------------------------
    $TABLE_INFILE = $ARGV[0] ;

    if ( $NARG >= 2 ) { $TABLE_ID  = $ARGV[1] ; }
    
    $i = 2;
    while ($i < $NARG) {

	$USE = 0 ;

	if ( $ARGV[$i] eq "NOHEADER" ) 
	{ $USE=1; $OPT_HEAD = "NOHEADER" ; $VFLAG = 0 ; }

	if ( lc($ARGV[$i]) eq '--o' ) 
	{ $USE=1; $i++ ; $OUTFILE_TEXT = $ARGV[$i] ; $VFLAG=0; }

	if ( lc($ARGV[$i]) eq '--a' ) 
	{ $USE=1; $i++; $TEXTFILE_APPEND = $ARGV[$i] ; $VFLAG=0; }

	if ( $VFLAG ) 
	{ $USE=1; @TABLE_VARLIST = (@TABLE_VARLIST , $ARGV[$i]) ; }
	if ( lc($ARGV[$i]) eq '--v'  ) { $USE=1; $VFLAG = 1; }
	
	$KEY = "outlier" ;
	if ( lc($ARGV[$i]) eq "--$KEY" )  { 
	    $USE=1; 
	    $KEY_OUTLIER = "$KEY";
	    $i++ ; $OUTLIER_NSIGMA[0] = $ARGV[$i] ; 
	    $i++ ; $OUTLIER_NSIGMA[1] = $ARGV[$i] ; 
	}

	$KEY = "outlier_sim" ;
	if ( lc($ARGV[$i]) eq "--$KEY" )  { 
	    $USE=1; 
	    $KEY_OUTLIER = "$KEY";
	    $i++ ; $OUTLIER_NSIGMA[0] = $ARGV[$i] ; 
	    $i++ ; $OUTLIER_NSIGMA[1] = $ARGV[$i] ; 
	}


	$KEY = "format" ;
	if ( lc($ARGV[$i]) eq "--$KEY" )  { 
	    $USE = 1 ;
	    $i++ ; $FORMAT_OUTFILE = $ARGV[$i] ;
	}

	if ( lc($ARGV[$i]) eq "obs" )  { 
	    $USE=1; 
	    $OPTOBS = 1;
	}

	if ( $USE == 0 ) {
	    $NBAD++ ;
	    print " ERROR: Unknown argument '$ARGV[$i]' \n";
	}

	$i++ ;
    }

    if ( $NBAD > 0 ) 
    { die "\n FATAL ERROR: $NBAD Unknown args \n"; }



    # if OUTFILE_TEXT is not set, then set to default name using TABLE_ID
    if ( length($OUTFILE_TEXT) == 0 )  { 
	$OUTFILE_TEXT = "sntable_dump_${TABLE_ID}.fitres"; 
	if ( $OUTLIER_NSIGMA[0] >= 0 ) 
	{ $OUTFILE_TEXT = "sntable_dump_${KEY_OUTLIER}.fitres"; }
	if ( $OPTOBS ) 
	{ $OUTFILE_TEXT = "sntable_dump_obs.fitres"; }
    }
    

    $NVARLIST = scalar(@TABLE_VARLIST);

    print " Input table file : $TABLE_INFILE \n";
    print " Input table ID   : $TABLE_ID \n";

    if ( $NVARLIST > 0 ) { 
	print " List of variables to dump: @TABLE_VARLIST \n" ; 
	print " Output text file : $OUTFILE_TEXT \n" ;
    }

    if ( length($TEXTFILE_APPEND) > 0 ) {
	print " Append variables to : $TEXTFILE_APPEND \n";
	$OPT_HEAD = "" ;  # must have default header
    }

    print "\n";

} # end of parse_args


# ===================
sub insert_ccid_varlist {
    
    #
    # Make sure that CCID is first element of variable list.
    # -> first remove CCID or CID from wherever it is.
    # -> add CCID as first element.
    #

    if ( $NVARLIST == 0        ) { return ; }
    if ( $TABLE_ID eq "GLOBAL" ) { return ; }

    # don't add CID for outlier option since CID will be added 
    # internally by sntable_dump.exe (Aug 4 2014)
    if ( $OUTLIER_NSIGMA[0] >= 0 ) { return ; }
    if ( $OPTOBS                 ) { return ; }

    my ($VAR);
    my @TMPLIST = @TABLE_VARLIST ;

    @TABLE_VARLIST = ( "CCID" );

    foreach $VAR (@TMPLIST) {
	if ( $VAR eq "CCID" ) { next; }
	if ( $VAR eq "CID"  ) { next; }
	@TABLE_VARLIST  = ( @TABLE_VARLIST , $VAR ) ;
    }

} # end of insert_ccid_varlist

# ==============================
sub make_dump_command {

    $CMD_DUMP = "$JOBNAME  $TABLE_INFILE  $TABLE_ID  $OPT_HEAD" ;

    if ( $NVARLIST > 0 ) { 
	$CMD_DUMP = "$CMD_DUMP --v @TABLE_VARLIST" ; 
	$CMD_DUMP = "$CMD_DUMP --o $OUTFILE_TEXT" ; 
    }

    if ( $OUTLIER_NSIGMA[0] >= 0 ) {
	$CMD_DUMP = "$CMD_DUMP --$KEY_OUTLIER @OUTLIER_NSIGMA" ;
	$CMD_DUMP = "$CMD_DUMP --o $OUTFILE_TEXT" ; 
    }

    if ( $OPTOBS ) {
	$CMD_DUMP = "$CMD_DUMP  OBS  --o $OUTFILE_TEXT" ; 
    }

    if ( $FORMAT_OUTFILE ne "DEFAULT" ) {
	$CMD_DUMP = "$CMD_DUMP --format $FORMAT_OUTFILE" ; 
    }

    my $SNANA_DIR = $ENV{'SNANA_DIR'};
    print " SNANA_DIR:  $SNANA_DIR \n";
    print " DUMP COMMAND: \n    $CMD_DUMP \n";

} # end of make_dump_command



# ===========================
sub get_fileType {

    # set global IFILETYPE (hbook or root)

    my ($TF, $SUF, $suf);

    $TF = $TABLE_INFILE ;
    $IFILETYPE = -9;

    foreach $SUF ( @VALID_SUFFIX_HBOOK) {
	$suf = lc($SUF); # allow lower case
	if ( index($TF,$SUF)  > 0 ) { $IFILETYPE = $IFILETYPE_HBOOK ; }
	if ( index($TF,$suf)  > 0 ) { $IFILETYPE = $IFILETYPE_HBOOK ; }
    }

    foreach $SUF ( @VALID_SUFFIX_ROOT) {
	$suf = lc($SUF); # allow lower case
	if ( index($TF,$SUF)  > 0 ) { $IFILETYPE = $IFILETYPE_ROOT ; }
	if ( index($TF,$suf)  > 0 ) { $IFILETYPE = $IFILETYPE_ROOT ; }
    }

    if ( $IFILETYPE < 0 ) {
	$MSGERR[0] = "Cannot determine file type (hbook or root)" ;
	$MSGERR[1] = "for $TF ." ;
	$MSGERR[2] = "Allowed table-file extensions are: ";
	$MSGERR[3] = " @VALID_SUFFIX_HBOOK  @VALID_SUFFIX_ROOT \n";
	sntools::FATAL_ERROR(@MSGERR); 
    }


} # end of get_fileType


# ==============================
sub check_TABLE_ID_ROOT {

    # for root file, allow old NTUPLE IDs 7100 or 7788 by
    # translating them into their root names: SNANA or FITRES.

    my ($TID_ORIG, $REPLACE) ;

    if ( $IFILETYPE != $IFILETYPE_ROOT ) { return ; }

    $TID_ORIG = $TABLE_ID ;
    $REPLACE = 0 ;

    if ( $TABLE_ID == 7100 ) 
    { $TABLE_ID = "SNANA" ;  $REPLACE = 1;  }


    if ( $TABLE_ID == 7788 ) 
    { $TABLE_ID = "FITRES" ;  $REPLACE = 1;  }

    if ( $REPLACE ) {
	print "\n";
	print "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
	print "    WARNING: \n" ;
	print "\t Replacing hbook NTID=$TID_ORIG  \n";
	print "\t with root table name: $TABLE_ID\n";
	print "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
	print "\n";
    }


} # end of check_TABLE_ID_ROOT

# ==============================
sub check_TABLE_ID_HBOOK {

    # for HBOOK file, allow root-string table names FITRES & SNANA.

    my ($TID_ORIG, $REPLACE) ;

    if ( $IFILETYPE != $IFILETYPE_HBOOK ) { return ; }

    $TID_ORIG = $TABLE_ID ;
    $REPLACE = 0 ;

    if ( $TABLE_ID eq "SNANA" ) 
    { $TABLE_ID = 7100 ;  $REPLACE = 1;  }


    if ( $TABLE_ID eq "FITRES" ) 
    { $TABLE_ID = 7788 ;  $REPLACE = 1;  }

    if ( $REPLACE ) {
	print "\n";
	print "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
	print "    WARNING: \n" ;
	print "\t Replacing root table=$TID_ORIG  \n";
	print "\t with HBOOK ntuple ID: $TABLE_ID\n";
	print "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
	print "\n";
    }


} # end of check_TABLE_ID_HBOOK

# ===========================
sub run_dump_command {

    my ($istat);

    if ( -e $OUTFILE_TEXT ) { qx(rm $OUTFILE_TEXT) ; }

    $istat = system("$CMD_DUMP");

    return $istat ;

} # end of run_dump_command


# ============================
sub textfile_append {

    if ( length($TEXTFILE_APPEND) == 0 ) { return ; }

    my ($CMD, $ARGF, $ARGO, $OUTFILE_PREFIX, $TEXTFILE_PREFIX, $j);

    # use text file prefix to construct output prefix to make
    # sure that it is unique -> avoid clobbering with multiple jobs.

    $j = index($TEXTFILE_APPEND,".");
    $TEXTFILE_PREFIX = substr($TEXTFILE_APPEND,0,$j);
    $OUTFILE_PREFIX  = "sntable_append_${TEXTFILE_PREFIX}" ;

    print "\n Append $TEXTFILE_APPEND with $OUTFILE_TEXT ... \n";
    
    $ARGF = "$TEXTFILE_APPEND $OUTFILE_TEXT $COMINE_FITRES_ARG[$IFILETYPE]" ;
    $ARGO = "--outprefix $OUTFILE_PREFIX" ;
    $CMD  = "combine_fitres.exe $ARGF $ARGO" ;
    
    print "CMD = $CMD \n";
    system("$CMD");

    # ------------------------------
    # insert few lines of commentary at top of appended fitres file.

    my $qq = '"' ;
    my @VARPRINT = @TABLE_VARLIST[1 .. $#TABLE_VARLIST] ;

    my $comment1 = "# Appended variables '@VARPRINT'" ;
    my $comment2 = "# from $TABLE_INFILE  (table $TABLE_ID)";
    my $comment3 = "# ";

    my ($sedcmd, $outFile, $tmpFile);
    $outFile = "${OUTFILE_PREFIX}.text" ;
    $tmpFile = "TEMP_$outFile" ;

    $sedcmd = "sed ";
    $sedcmd = "$sedcmd -e 1i${qq}${comment1}${qq}";
    $sedcmd = "$sedcmd -e 1i${qq}${comment2}${qq}";
    $sedcmd = "$sedcmd -e 1i${qq}${comment3}${qq}";
    qx($sedcmd $outFile > $tmpFile ; mv $tmpFile $outFile);

} # end of  textfile_append
