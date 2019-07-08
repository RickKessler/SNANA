#!/usr/bin/perl
#
# Jan 2014 R.Kessler
#
# Translate csv formatted file into a file with
# the snana/fitres format:
#  VARNAMES: ROW <VARNAME1> <VARNAME2> <VARNAME3> ...
#  ROW: 1  <value1> <value2> <value3> ...
#  ROW  2  <value1> <value2> <value3> ...
#  etc ...
#
# If the first variable name is NOT CID or ROW or STARID, 
# then a ROW column is inserted.  Either comma or blank-space 
# delimeters are allowed.
#
#
# Usage
#  convertcsv2snana.pl <inFile>
#
#  convertcsv2snana.pl <inFile> CAPS
#       [all varnames -> uppercase]
#
#  convertcsv2snana.pl <inFile> CAPS   --VARNAMES <list of varnames>
#     [replace varnames list]
#
#
# May 12 2016: fix to allow '#' in first line.
# Jun 27 2019: remove NVAR key that is obsolete.
# ----------------------

use strict ;

my ($INFILE, $OPT_CAPS, @VARLIST_USER) ;
my ($OUTFILE, $TMPFILE, $SEDCMD);
my ($NVAR, $VARNAMES_ORIG, @VARLIST_ORIG, $VARNAMES_OUT);

sub parse_args ;
sub get_VARNAMES  ;
sub make_outFile ;

# ============== START MAIN ===============

&parse_args();


# remove commas and comment signs.
# ubstitute up to two ,, with ,-9999' in case of blanks; 
# this does not protect blank in last column.
$SEDCMD = "sed " .
    "-e 's/,,/,-9999,/g' " .
    "-e 's/,,/,-9999,/g' " . 
    "-e 's/,/ /g' " .
    "-e 's/#//g'" ; 
qx($SEDCMD $INFILE > $TMPFILE);

#  list of varnames
&get_VARNAMES();

&make_outFile();

print "\n Done. See $OUTFILE \n";

# =========== END MAIN =============


sub parse_args {
    
    my ($i, $NARG, $USERFLAG );

    $OPT_CAPS = 0 ;
    @VARLIST_USER = ()  ;
    $USERFLAG = 0 ;

    $NARG = scalar(@ARGV);

    if ( $NARG < 1 ) {
	die "\n ERROR: must specify input file \n";
    }

    
    $INFILE  = $ARGV[0];
    $OUTFILE = "${INFILE}_converted" ;
    $TMPFILE = "TMP_${INFILE}" ;

    for($i=1; $i < $NARG; $i++ ) {
	if ( $ARGV[$i] eq "CAPS"       ) { $OPT_CAPS = 1; }

	if ( $USERFLAG ) { @VARLIST_USER = (@VARLIST_USER, $ARGV[$i]); }
	if ( $ARGV[$i] eq "--VARNAMES" ) { $USERFLAG = 1; }
	if ( $ARGV[$i] eq "-VARNAMES"  ) { $USERFLAG = 1; }
    }

    print "\n";
    print " Translate $INFILE  -> $OUTFILE \n";


}  # end of parse_args


# =====================================
sub get_VARNAMES {

    open PTR_INFILE , "< $INFILE" ;
    $VARNAMES_ORIG = <PTR_INFILE> ;
    close PTR_INFILE ;

    if ( $OPT_CAPS ) { $VARNAMES_ORIG = uc($VARNAMES_ORIG); }

    # replace commas in VARNAMES with blank spaces so that
    # either comma or space delimeter is allowed
    my ($VARNAMES_noComma);
    $VARNAMES_noComma = $VARNAMES_ORIG ;  
    $VARNAMES_noComma =~ s/,/ /g ;
    $VARNAMES_noComma =~ s/#//g ;  # remove comment sign (May 2016)

    @VARLIST_ORIG =  split(" ", $VARNAMES_noComma);

    $NVAR = scalar(@VARLIST_ORIG) ;

    print " VARNAMES_ORIG($NVAR) = @VARLIST_ORIG \n";

    $VARNAMES_OUT = "@VARLIST_ORIG" ;

    # check for user varnames
    my $NVAR_USER = scalar(@VARLIST_USER);
    if ( $NVAR_USER > 0 ) {
	
	if ( $NVAR_USER != $NVAR ) 
	{ die "\n ERROR: NVAR=$NVAR but NVAR_USER=$NVAR_USER \n"; }

	$VARNAMES_OUT = "@VARLIST_USER" ;	
    }

} # end of get_VARNAMES 


# ========================
sub make_outFile {

    # transfer contents to OUTFILE, and add ROW variable
    # in first column so that 'combine_fitres.exe' works.

    my (@CONTENTS, $NROW, $NLINE, $LINE );

    # read everything from input file
    open PTR_TMPFILE , "< $TMPFILE" ;
    @CONTENTS = <PTR_TMPFILE> ;
    close PTR_TMPFILE ;
    qx(rm $TMPFILE);
    
    my ($FIRSTVAR, $VARNAME_ROW, $ADDROW, $cNROW);
    $ADDROW      = 1 ;
    $FIRSTVAR    = $VARLIST_ORIG[0];

    if ( $FIRSTVAR eq "CID"    ) { $ADDROW = 0 ; }
    if ( $FIRSTVAR eq "cid"    ) { $ADDROW = 0 ; }
    if ( $FIRSTVAR eq "ROW"    ) { $ADDROW = 0 ; }
    if ( $FIRSTVAR eq "GALID"  ) { $ADDROW = 0 ; }
    if ( $FIRSTVAR eq "STARID" ) { $ADDROW = 0 ; }

    if ( $ADDROW ) 
    { $VARNAME_ROW = "ROW" ; $NVAR = $NVAR + 1 ; }
    else
    { $VARNAME_ROW = "" ; }
	

    # open  OUTFILE
    open PTR_OUTFILE , "> $OUTFILE" ;
    
# xxx mark delete Jun 2019  print PTR_OUTFILE "NVAR: $NVAR \n";
    print PTR_OUTFILE "VARNAMES: $VARNAME_ROW $VARNAMES_OUT \n";

    $NROW = $NLINE = 0 ;
    foreach $LINE (@CONTENTS) {
	$NLINE++ ;
	if ( $NLINE == 1 ) { next ; } # skip variable definitions
	
	$LINE =~ s/\s+$// ; # remove trailing blanks
	$NROW++ ;	
	$cNROW = sprintf("%3d", $NROW); 

	if ( $ADDROW ) 
	{  print PTR_OUTFILE "ROW: $cNROW  $LINE \n"; }
	else
	{  print PTR_OUTFILE "ROW: $LINE \n"; }
    }

    print " Finished writing $NROW rows to: $OUTFILE \n";
    close PTR_OUTFILE ;

} # end of make_outFile
