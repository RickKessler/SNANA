#!/usr/bin/perl
#
# April 2014, R.Kessler
# Replace CERNLIB's ypatchy to convert .cra file into
# a fortran-compilable file. This script was written
# because of the difficulty in getting the original 
# ypatchy source code to work on 64 bit machines.
#
# Beware that this script replaces only part of the original 
# ypatchy functionality corresponding to what is used in snana. 
# As a precaution, this script will abort if any unrecogonized 
# patchy flags (e.g., +BLA) are found.
#
# Usage:
#  ypatchy.pl - <outFile>  <inFile>  .go
#
# + No options are currently supported, so 1st argument must be '-'
#    to indicate "no options."
# + If <inFile> has no extension, '.cra' is assumed.
#     Thus  <inFile> = 'snana' -> snana.cra
# + if .go is not given as 4th arg, script aborts since there
#     is no interactive mode.
#
# ================

use strict ;
use IO::Handle ;
# use List::MoreUtils 'any' ;

my $logFile = "ypatchy.log" ;
my $MXLINE_OUT  = 1000000 ;  # abort if more than this many lines written
my $MXCHAR_NAME = 24 ;       # max char-len of DECK or KEEP name

my ($USER_OPTIONS, $INPFILE, $OUTFILE, $GOFLAG);
my (@FILE_CONTENTS, @KEEPLINE);

my @USEFLAG_LIST   = () ;  # store args of USE,XXX.
my @PATCH_COM_LIST = () ;  # store list of common-block patches
my @PATCH_EXE_LIST = () ;  # store list of exe-code patches

my @DECKLIST_NAME     = () ;
my @DECKLIST_MINLINE  = () ;
my @DECKLIST_MAXLINE  = () ;

my @KEEPLIST_NAME     = () ;
my @KEEPLIST_MINLINE  = () ;
my @KEEPLIST_MAXLINE  = () ;

my $BLOCK_TYPE_CURRENT ;
my $BLOCK_TYPE_NONE = 0 ;
my $BLOCK_TYPE_KEEP = 1 ;
my $BLOCK_TYPE_DECK = 2 ;

my $PATCH_FLAG_CURRENT ;
my $PATCH_FLAG_NONE   =  0;
my $PATCH_FLAG_VALID  = +1;
my $PATCH_FLAG_IGNORE = -1;

my ($NFLAG, $NPATCH_COM, $NPATCH_EXE,  $NDECK, $NKEEP);
my ($NSELF_USE, $NSELF_IGNORE);
my ($NLINE_STORE ,$NLINE_OUT);
my ($SELF_IGNORE_FLAG );


sub parse_args ;
sub parseHead_INPFILE ;
sub initStuff ;

sub parseFile ;
sub parseLine ;
sub parseLine_SELF ;
sub parseLine_PATCH ;
sub parseLine_PAM ;
sub parseLine_USEFLAGS ;

sub parseFLAGLIST ;
sub storeLine_Fortran ;

sub cleanString ;
sub checkNameLen ;

sub write_forFile ;
sub writeDECK ;
sub insertKEEP ;

sub printSummary ;

# ============== BEGIN MAIN ===============

if ( -e $logFile ) { qx(rm $logFile); }
open PTR_LOG , "> $logFile" ;

&parse_args();

&initStuff();

&parseFile($INPFILE);

&write_forFile($OUTFILE);

&printSummary();


# ============ END MAIN ====================

sub parse_args {

    my ($NARG, $tmp, $jdot);

    $NARG = scalar(@ARGV);
    $INPFILE = $OUTFILE = $USER_OPTIONS = "" ;
    $GOFLAG  = 0 ;

    if ( $NARG < 4 ) {
	die "\n FATAL ERROR: execting 4 arguments.\n";
    }

    # ---------------------------------
    # parse options
    $USER_OPTIONS = substr($ARGV[0],1,99);
    if ( $ARGV[0] ne '-' ) {
	die "\n FATAL ERROR: first arg must be '-' since no options are supported.\n";
    }

    $OUTFILE = $ARGV[1];
    $INPFILE = $ARGV[2];    
    if ( $ARGV[3] eq '.go' ) { $GOFLAG = 1; }

    $jdot = index($INPFILE,".") ;

    if ( $jdot < 0 ) {
	# no extension -> assume .cra
	$INPFILE = "${INPFILE}.cra" ;
    }

    print PTR_LOG " Process input file $INPFILE \n";
    print PTR_LOG " Will create output file: $OUTFILE \n";
    print PTR_LOG "\n" ;

    unless ( $GOFLAG ) {
	die "\n FATAL ERROR: no interactive mode; must use '.go' option\n";
    }

} # end of parse_args


# =======================================
sub initStuff {
    
    $NFLAG = $NPATCH_COM = $NPATCH_EXE = 0;
    $NDECK = $NKEEP = 0 ;

    $NSELF_USE = $NSELF_IGNORE = 0 ;

    $PATCH_FLAG_CURRENT  = $PATCH_FLAG_NONE ;
    $SELF_IGNORE_FLAG = 0 ;

    $NLINE_STORE = $NLINE_OUT = 0 ;

}  # end of initStuff


# ========================
sub parseFile {
    my ($inFile) = @_ ;

    my ($NLINE, @contents, $line) ;

    open PTR_INFILE , "< $inFile" ;
    @contents = <PTR_INFILE> ;  $NLINE = scalar(@contents);
    close PTR_INFILE ;


    print PTR_LOG "   Parse file: $inFile  ($NLINE lines) \n" ;

    $PATCH_FLAG_CURRENT = $PATCH_FLAG_NONE ;  # reset global flag.
    $BLOCK_TYPE_CURRENT = $BLOCK_TYPE_NONE ;

    foreach $line (@contents) { 
	$line =~ s/\s+$// ;  # trim trailing whitespace
	&parseLine($inFile,$line); 
    }

    $PATCH_FLAG_CURRENT = $PATCH_FLAG_NONE ;  # reset global flag.

} # end of parseFile


# =======================
sub parseLine {

    # parse $line.
    # Note that input inFile is used only for printing comments.

    my ($inFile,$line) = @_ ;

    my (@PARTS, $ARG0, $ARG1, $ARGLAST, $NTMP);
    my ($L_USE, $L_PAM, $L_CDE, $L_SELF);
    my ($L_KEEP, $L_DECK, $L_PATCH, $L_NAMES);
    my ($jdot, $jIF, $KEEP, $DECK, $newFile );
    my ($FLAGLIST, $USEFLAG ) ;
    my $FF = "in file=$inFile";

    # bail if 1st char is not a '+'
    if ( substr($line,0,1) ne "+" ) { 
	&storeLine_Fortran($line,0) ;
	return ; 
    }

    # abort if last char is not a dot
    $jdot    = index($line,".");
    if ( $jdot < 0 ) {
	die "\n FATAL ERROR: Missing period at end of \n '$line' \n";
    }

    # strip of optional list of flags after IF= and set USEFLAG
    $jIF  = index($line,"IF=");
    if ( $jIF > 0 ) {
	my $j1       = $jIF + 3 ;
	my $jlen     = $jdot - $j1  ;
	$FLAGLIST = substr($line, $j1, $jlen) ; 
	$USEFLAG  = &parseFLAGLIST($FLAGLIST); 

	if ( index($line,"SNLCPLOT") > 0 ) {
#	    print " xxx FLAGLIST='$FLAGLIST for '$line' (USE=$USEFLAG) \n";
	}
    }
    else
    { $USEFLAG = 1 ; }


    # break line into comma-separated strings.
    @PARTS  = split(/,/,$line);
    $ARG0   = &cleanString($PARTS[0]);  # remove spaces, dot, *
    $ARG1   = &cleanString($PARTS[1]);

    # ----------------------------------------
    # define all valid patchy flags here;
    # anything else will cause an abort.
    $L_USE    = ( $ARG0 eq "+USE"   ) ;
    $L_PAM    = ( $ARG0 eq "+PAM"   ) ;
    $L_CDE    = ( $ARG0 eq "+CDE"   ) ;
    $L_SELF   = ( $ARG0 eq "+SELF"  ) ;
    $L_KEEP   = ( $ARG0 eq "+KEEP"  ) ;
    $L_DECK   = ( $ARG0 eq "+DECK"  ) ;
    $L_PATCH  = ( $ARG0 eq "+PATCH" ) ;
    $L_NAMES  = ( $ARG0 eq "+NAMES" ) ;   


    # turn off line-storage, except for CDE and SELF.
    unless ( $L_CDE || $L_SELF )  
    { $BLOCK_TYPE_CURRENT = $BLOCK_TYPE_NONE ; }

    my $VALID_PATCH = ($PATCH_FLAG_CURRENT == $PATCH_FLAG_VALID);

    if ( $L_USE ) {
	my $jP      = index($line,"P=");
	if ( $jP < 0 ) { 
	    # not a patch, hence it's a patchy flag
	    &parseLine_USEFLAGS($line);
	}
	else {
	    # it's a patch
	    unless ( $PATCH_FLAG_CURRENT == $PATCH_FLAG_IGNORE ) 
	    {  &parseLine_PATCH($inFile,$line); }
	}
    }  

    elsif ( $L_PAM ) {
	# check for file. 
	my $newFile = &parseLine_PAM($line);
	&parseFile($newFile);
    }

    elsif ( $L_PATCH ) {
	# found a patch of code or common blocks
	my ($PATCH, $ISCOM, $ISEXE);
	$PATCH  = $ARG1;

	# check if this patch is defined in the header
	$ISCOM = ( grep /$PATCH/, @PATCH_COM_LIST );
	$ISEXE = ( grep /$PATCH/, @PATCH_EXE_LIST );

	if ( $ISCOM || $ISEXE ) {
	    print PTR_LOG "\t Found defined PATCH = '$PATCH'  $FF \n";
	    $PATCH_FLAG_CURRENT = $PATCH_FLAG_VALID ;
	}
	else {
	    print PTR_LOG "\t Found undefined PATCH  = '$PATCH'  $FF \n";
	    $PATCH_FLAG_CURRENT = $PATCH_FLAG_IGNORE ;
	}
    }

    elsif ( $L_DECK && $VALID_PATCH  && $USEFLAG ) {
	# found DECK -> code
	&checkNameLen("DECK",$ARG1);
	$BLOCK_TYPE_CURRENT = $BLOCK_TYPE_DECK ;
	$NDECK++ ;
	$DECKLIST_NAME[$NDECK]    = $ARG1 ; 
	$DECKLIST_MINLINE[$NDECK] = -1 ;
	$DECKLIST_MAXLINE[$NDECK] = -1 ;
    }

    elsif ( $L_KEEP && $VALID_PATCH && $USEFLAG ) {
	# found common
	&checkNameLen("KEEP",$ARG1);
	$BLOCK_TYPE_CURRENT = $BLOCK_TYPE_KEEP ;
	$NKEEP++ ;
	$KEEPLIST_NAME[$NKEEP]    = $ARG1 ; 
	$KEEPLIST_MINLINE[$NKEEP] = -1 ;
	$KEEPLIST_MAXLINE[$NKEEP] = -1 ;
    }
    elsif ( $L_SELF ) {

	&parseLine_SELF($line);

    }
    elsif ( $L_CDE && $USEFLAG ) {
	&storeLine_Fortran($line,1) ; 
    }
    elsif ( $L_NAMES || $L_DECK || $L_KEEP || $L_CDE  ) {
	# allow, but ignore
    }
    else {
	my $CERR = "FATAL ERROR" ;
	die "\n $CERR: undefined patchy line\n\t $line \n $FF\n";
    }

} # end of parseLine


# =================================
sub storeLine_Fortran {

    # store this line of fortran.
    # isKEEP = 1 if line is +CDE to include a KEEP; 0 otherwise.

    my ($line, $isKEEP) = @_ ;

    my ($MINLINE, $VALID_PATCH);

    $VALID_PATCH = ($PATCH_FLAG_CURRENT == $PATCH_FLAG_VALID);

    if ( !$VALID_PATCH      ) { return ; }  # ignore patch
    if (  $SELF_IGNORE_FLAG ) { return ; }  # ignore SELF code

    
    $NLINE_STORE++ ;
    @FILE_CONTENTS[$NLINE_STORE] = "$line" ;
    @KEEPLINE[$NLINE_STORE]      = $isKEEP ;

    if ( $BLOCK_TYPE_CURRENT == $BLOCK_TYPE_KEEP ) {
	$MINLINE = $KEEPLIST_MINLINE[$NKEEP] ;
	if ( $MINLINE < 0 ) { $KEEPLIST_MINLINE[$NKEEP] = $NLINE_STORE; }
	$KEEPLIST_MAXLINE[$NKEEP] = $NLINE_STORE ;

    }
    elsif ( $BLOCK_TYPE_CURRENT == $BLOCK_TYPE_DECK ) {
	$MINLINE = $DECKLIST_MINLINE[$NDECK] ;
	if ( $MINLINE < 0 ) { $DECKLIST_MINLINE[$NDECK] = $NLINE_STORE; }
	$DECKLIST_MAXLINE[$NDECK] = $NLINE_STORE ;
    }

} # end of  storeLine_Fortran 

# =================================
sub parseLine_PATCH {

    # this line has already been identified as a PATCH-use
    # definition with +USE,P=patchName.
    # Don't confuse this with +PATCH.

    my ($inFile,$line) = @_ ;

    my (@PARTS, $PART, $part, $len, $PATCH, $IGNORE, $UU ) ;
    my $USEFLAG = "" ;
    my $TYPE    = "COM" ;  # COM=common block, EXE-> exec code

    @PARTS = split(/,/,$line);
    foreach $PART ( @PARTS ) {	       
	$part  = &cleanString($PART) ; # remove spaces and dot
	$len = length($part);
	if ( substr($part,0,2) eq "P=" )  { 
	    $PATCH  = substr($part,2,$len-2); 
	    $PATCH  =~ tr/\*//d ;  # remove asterisk
	}

	if ( $part eq "T=EXE" ) { $TYPE = "EXE" ;  }

	if ( substr($part,0,3) eq "IF=" ) 
	{ $USEFLAG = substr($part,3,$len-3); }
    }

    $IGNORE = 0 ;
    $UU = "" ;

    # check IF-conditional
    if ( length($USEFLAG) > 0 ) {
	unless ( grep /$USEFLAG/, @USEFLAG_LIST ) 
	{ $IGNORE = 1; }
	$UU = "(USEFLAG=$USEFLAG)" ;
    }

    # bail if patch flag is not defined.
    if ( $IGNORE ) { return ; }
    
    print PTR_LOG "\t Define $TYPE PATCH = $PATCH  in file=$inFile $UU\n" ;

    if ( $TYPE eq "COM" ) {
	$NPATCH_COM++ ;
	@PATCH_COM_LIST = ( @PATCH_COM_LIST, $PATCH );
    }
    else {
	$NPATCH_EXE++ ; 
	@PATCH_EXE_LIST = ( @PATCH_EXE_LIST, $PATCH );
    }   

    $BLOCK_TYPE_CURRENT = $BLOCK_TYPE_NONE ;

} # end of parseLine_PATCH


# ===============================
sub parseLine_PAM {

    # parse line containing 
    # +PAM ... T=A,C. <fileName>
    # and return the file name.
    
    my ($line) = @_ ;

    my (@PARTS, $endLine, $newFile, $TKEY, $JT);

    $TKEY = "T=A,C." ;
    $JT   = index($line, $TKEY );
    if ( $JT < 0 ) {
	my $c1err = "FATAL ERROR on line = '$line'";
	my $c2err = "+PAM line must contain $TKEY" ;
	die "\n $c1err \n $c2err \n";
    }

    # strip word after TKEY; allow for leading spaces and trailing comments.
    $endLine  = substr($line, $JT+length($TKEY), 99 );
    $endLine  =~ s/^\s+// ; # remove leading spaces

    @PARTS    = split(/ /,$endLine);
    $newFile  = $PARTS[0] ;  # first word after $TKEY

    return $newFile ;

} # end of parseLine_PAM 

# ===============================
sub parseLine_SELF {

    # the passed $line argument starts with '+SELF'.
    # parse the entire line and implement IF=SELF logic.
    # e.g., +SELF,IF=SNANA,SNFIT.
    # will check both the SNANA and SNFIT flags.
    # '+SELF.' ends the if-block.
    # Global flag $SELF_IGNORE_FLAG controls if the
    # lines insdie the IF=SELF block are included in the out-file.

    my ($line) = @_ ;

    my ($jdot, $jeq, $FLAG, $flag, $USEME);
    my ($FLAGLIST, @FLAGARRAY, $NEG );

    $jeq  = index($line,'=') ;
    $jdot = index($line,".");

    # extract comma-separate list of flags
    if ( $jeq > 0 ) { 
	# scoop up everything after '=' sign
	$FLAGLIST = substr($line,$jeq+1,$jdot-$jeq-1); 
    }
    else { 
	# just  plain +SELF. to end if-block
	$SELF_IGNORE_FLAG = 0 ; 
	return ; 
    }

    $USEME = &parseFLAGLIST($FLAGLIST);

    # ------------------------------------------------
    # apply SELF logic based on the local USEME flag.

    if ( $USEME ) {
	$NSELF_USE++ ;           # for monitor only
	$SELF_IGNORE_FLAG = 0 ;  # include code inside +SELF block
    }
    else {
	$NSELF_IGNORE++ ;        # for monitor only
	$SELF_IGNORE_FLAG = 1 ;  # ignore code inside +SELF block
    }


} # end of parseLine_SELF


# =======================================
sub parseLine_USEFLAGS {

    my ($line) = @_ ;

    # parse statement of the form
    #    +USE,FLAG1,FLAG2,FLAG3,ETC.
    # Add all flags to @USEFLAG_LIST.
    # Abort if any flag is repeated.

    my ($flagList, @flagList, $FLAG, $c1err, $c2err);

    if ( substr($line,0,4) ne "+USE" ) {
	$c1err = "FATAL ERROR: line must start with '+USE'" ;
	$c2err = "but line='$line'" ;
	die "\n$c1err \n\t $c2err\n";
    }

    $flagList = substr($line,5,99);
    $flagList = &cleanString($flagList);
    @flagList = split(/,/,$flagList);

    foreach $FLAG (@flagList) {
	$NFLAG++ ;
	$FLAG = &cleanString($FLAG);

	if ( grep /$FLAG/, @USEFLAG_LIST ) {
	    $c1err = "FATAL ERROR:  patchy-USE flag = '$FLAG'";
	    $c2err = "is defined multiple times.";
	    die "\n$c1err \n\t $c2err\n";
	}

	@USEFLAG_LIST = (@USEFLAG_LIST , $FLAG) ;
	print PTR_LOG "\t Define patchy flag  '$FLAG'  \n";
    }

} # end of parseUSEFLAGS

# =======================================================
sub parseFLAGLIST {

    # parse comma-separated list of patchy flags and return true
    # if any patchy is flag. Also account for negative logic.
    #
    # Examples:
    #  "SNANA,SNFIT"   returns true if either SNANA or SNFIT is set.
    #  "-HBOOK,SNANA"  returns true if HBOOK is not set or SNANA is set.

    my ($FLAGLIST) = @_ ;

    # local args
    my (@FLAGARRAY, $FLAG, $flag, $NEG, $USEFLAG);

    $USEFLAG = 0 ;  # init output logical

    # break up FLAGLIST into array
    @FLAGARRAY = split(/,/,$FLAGLIST);

    # loop over flags and check them all;
    # if any flag is set then -> USEFLAG=1

    foreach $FLAG (@FLAGARRAY) {
	$FLAG  = &cleanString($FLAG);
	&checkNameLen("FLAG",$FLAG);  # check string length
	$NEG = ( substr($FLAG,0,1) eq '-' ) ;
	if ( $NEG ) {
	    # inverted logic with -[FLAG]
	    $flag = substr($FLAG,1,30);  # remove '-' to get flag name
	    if ( !( grep /$flag/, @USEFLAG_LIST ) ) { $USEFLAG = 1 ; }
	}
	else {
	    # positive/normal logic
	    if ( grep /$FLAG/, @USEFLAG_LIST ) { $USEFLAG = 1; }
	}
    }
    
    return $USEFLAG ;

} # end of parseFLAGLIST


# ========================================
sub cleanString {
    my ($string) = @_ ;

    # remove 
    #  - extra spaces (leading or trailing), 
    #  - dot  
    #  - asterisk 
    # from string.
    
    my $outString = $string ;
    $outString   =~ tr/.\*//d ;  # remove dot and optional asterisk.
    $outString   =~ s/\s+$// ;   # trim trailing whitespace
    $outString   =~ s/^\s+// ;   # remove leading spaces

    return $outString ;

} # end of clean string

# ========================================
sub checkNameLen {

    # anort if length of NAME is too long.

    my ($WHAT,$NAME) = @_ ;

    if ( length($NAME) > $MXCHAR_NAME ) {
	my $MM = "MXCHAR_NAME=$MXCHAR_NAME" ;
	die "\n FATAL ERROR: $WHAT-name '$NAME' is too long ($MM)\n";
    }

} # end of checkNameLen

# ========================================
sub write_forFile {

    # write output fortran file
    my ($outFile) = @_ ;

    my ($ideck) ;

    print PTR_LOG "\n Open output file: $outFile \n";

    if ( -e $outFile ) { qx(rm $outFile); }
    open PTR_OUTFILE , "> $outFile" ;

    for($ideck=1; $ideck <= $NDECK; $ideck++ ) 
    {  &writeDECK($ideck); }
    
    close PTR_OUTFILE ;

} # end of write_forFile


# ========================================
sub writeDECK {
    my ($ideck) = @_ ;

    my ($NAME, $MINLINE, $MAXLINE, $LINE, $iline);

    $NAME    = $DECKLIST_NAME[$ideck] ;
    $MINLINE = $DECKLIST_MINLINE[$ideck] ;
    $MAXLINE = $DECKLIST_MAXLINE[$ideck] ;

    print PTR_OUTFILE "CDECK  ID>,  $NAME. \n";

    for($iline=$MINLINE ; $iline <= $MAXLINE; $iline++ ) {

	$LINE = "$FILE_CONTENTS[$iline]";
	if ( $KEEPLINE[$iline] ) {
	    &insertKEEP($LINE);
	}
	else {
	    print PTR_OUTFILE "$LINE\n";
	    $NLINE_OUT++ ;
	}

    }  # end of iline loop


#    PTR_OUTFILE -> autoflush(1);  # way slower !!!

} # end of writeDECK

# ========================================
sub insertKEEP {
    # LINE stars with +CDE,XXX ;
    # insert KEEP lines in PTR_OUTFILE.
    # Beware of IF-logic such as
    #  +CDE,XYZ,IF=ABC.
    # In this case, strip of the XYZ keep-name.
    
    my ($LINE) = @_ ;

    # local var
    my ($j1, $jdot, $jlen, $NAME, $indx, $keepLine, $tmpString, @tmpList);
    
    # --------------------------------
    # in case there is IF= logic, break string into comma-separated
    # strings and use first string.
    $j1   = index($LINE,",") + 1 ; # skip +CDE,
    $jdot = index($LINE,".");
    $jlen = $jdot - $j1 ;
    $tmpString = substr($LINE, $j1, $jlen);  # scoop everything after CDE,
    @tmpList   = split(/,/,$tmpString);      # comma-separate list
    $NAME      = $tmpList[0] ;               # first one on list

    # find keep-index
    my( @index )= grep { $KEEPLIST_NAME[$_] eq $NAME } 0..$#KEEPLIST_NAME;
    $indx = $index[0];

    # make sure we found keep-index
    if ( $indx < 0 || $indx eq "" ) {
	die "\n FATAL ERROR: cannot find index for $LINE \n";
    }

    # double check. If keep name is BLA, note that BLABLA will find index.
    # Thus require EXACT match
    my $tmpName = "$KEEPLIST_NAME[$indx]" ;
    if ( $NAME ne $tmpName ) {
	my $c1err = "FATAL ERROR: could not find KEEP $NAME ;" ;
	my $c2err = "Instead found KEEP $tmpName (indx=$indx)" ;
	die "\n$c1err \n\t $c2err \n";
    }

    my ($MINLINE, $MAXLINE, $iline);
    $MINLINE = $KEEPLIST_MINLINE[$indx] ;
    $MAXLINE = $KEEPLIST_MAXLINE[$indx] ;

    print PTR_OUTFILE "\nC $LINE inserted below by ypatchy. \n";

    for($iline=$MINLINE; $iline <= $MAXLINE; $iline++ ) {

	$keepLine = "$FILE_CONTENTS[$iline]" ;

	# check if this line is another CDE statement
	if ( substr($keepLine,0,4) eq "+CDE" ) {
	    &insertKEEP($keepLine);
	}
	else {
	    print PTR_OUTFILE "$keepLine\n" ;
	    $NLINE_OUT++ ;
	}
    }

    if ( $NLINE_OUT > $MXLINE_OUT ) {
	die "\n FATAL ERROR: $NLINE_OUT fortran lines exceeds bound.\n";
    }

} # end of insertKEEP

# ========================================
sub printSummary {

    print PTR_LOG "\n" ;
    print PTR_LOG " Done writing $NLINE_OUT fortran lines with \n";
    print PTR_LOG "\t $NDECK  DECKs \n";
    print PTR_LOG "\t $NKEEP  KEEPs \n";
    print PTR_LOG "\t $NSELF_USE    defined +SELFs --> included.\n";
    print PTR_LOG "\t $NSELF_IGNORE  undefined +SELFs --> ignored. \n";
    close PTR_LOG ;

    print " Done: see summary in $logFile \n";

} # end of &printSummary
