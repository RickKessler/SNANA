#!/usr/bin/env perl
#
# Created Jun 20, 2011 by R.Kessler
# Useful perl utilities.
#
# Aug 04, 2011: use strict
# Sep 10, 2011: add sed_strReplace
# Feb 01, 2012: add get_namelistValue
# Jun 28, 2012: add checkDirExist
# Feb 16, 2013: add make_batchFile
# Mar 25, 2013: add rdfitres_getcol
# Apr 30, 2013: add check_multipleJobs
# May 26, 2013: remove doneFile and XMCD args from make_batchFile
# Jun 04, 2013: add checkFileExist
# Sep 20, 2013: add uniq_entries, lininterp, weighted_avg
#                   arith_avg
# May 02, 2014: add parse_array
# Jan 05, 2017: add extractStringOpt
# Dec 03, 2017: add clean_fileList
# Apr 17, 2019: add FATAL_ERROR_STAMP function
# Jan 31, 2020: 
#   + replace a few "ERROR:" messages with FATAL_ERROR util.
#   + update FATAL_ERROR to look similar to C code error util.
#
# Aug 24 2020:  loadArray_fromFile aborts on CONFIG: key
#
package sntools ;
use List::Util qw(min max);
use strict ;

# declare utilities
sub parse_line(@) ;
sub parse_array(@);
sub parse_value_after_key(@);
sub loadArray_fromFile(@) ;
sub loadArray_excludeLines(@) ;
sub Killjobs(@) ;
sub FATAL_ERROR(@);
sub FATAL_ERROR_STAMP(@);
sub printArray(@);
sub sed_nmlInsert(@);
sub sed_strReplace(@) ;
sub shellName ;
sub setEnvString(@);
sub unsetEnvString(@);
sub userPrompt(@) ;
sub get_namelistValue(@);
sub checkDirExist(@);
sub checkFileExist(@);
sub write_loginMacroLines(@);
sub make_batchFile(@) ;
sub rdfitres_getcol(@); 
sub rdfitres_util_checkhead(@); 
sub check_multipleJobs(@) ;
sub uniq_entries(@) ;
sub weighted_avg(@) ;
sub lininterp(@) ;
sub arith_avg(@) ;
sub extractStringOpt(@) ;
sub clean_fileList(@) ;

sub request_submit_batch_jobs(@);

# ===============================================
sub parse_line(@) {

  my($inFile,$nwd_parse,$string,$opt) = @_;

# grep inFile_master for $string, and return argument after $string
# nwd_parse = number of words to return after $string
# If nwd_parse = large number (99), then just return all
# words after $string.
#
# $opt = 0 => ABORT if string not found
# $opt = 1 => warning if string not found, but do not abort
# $opt = 2 => no warning if string not found
#
# Mar 30, 2011: add opt=2 by J.Mosher
# Sep 07, 2011: if nwd_parse=0, just check that key exists.
# Dec 02, 2017: use zgrep on gzipped files

  my $fnam = "parse_line";

  my (@tmpall, @line, $Nline, $iline, @fields, $tmp ); 
  my ($iwd, $wd, $IWD_MATCH );
  my $qq = '"' ;
  my $q  = "'" ;
  # ----------- BEGIN -------------

  @tmpall ;
  if ( $inFile eq '' ) { return @tmpall ; }

  if ( index($inFile,".gz") > 0 )
  {  @line = qx { zgrep  ${qq}${string}${qq}  $inFile }  ; }
  else
  {  @line = qx { grep  ${qq}${string}${qq}  $inFile }  ; }      

  
  $Nline = scalar(@line) ;
  @tmpall = ( ) ;

  if ( $nwd_parse == 0 ) {
      @tmpall = @line ;
      $Nline  = 0;
  }


  for ($iline = 0; $iline < $Nline ; $iline++ ) {
      @fields = split(/\s+/,$line[$iline]) ;

      # make sure that first element exactly matches $string
      $iwd = 0; $IWD_MATCH = -9 ;
      foreach $wd ( @fields ) {
	  if ( "$wd" eq "$string" ) { $IWD_MATCH = $iwd ; }
	  $iwd++ ;
      }
      if ( $IWD_MATCH < 0 ) { next ; }

      my $Nfield  = scalar(@fields) - 1 - $IWD_MATCH ;
      my $nwd     = min ( $nwd_parse , $Nfield );
      my $firstwd = $IWD_MATCH + 1 ; 
      my $lastwd  = $IWD_MATCH + $nwd ;

      # start with word after string
      $tmp = "@fields[$firstwd .. $lastwd]" ;
      @tmpall = ( @tmpall , "$tmp" ) ;

  }


# abort if keywd not found

  if ( "$tmpall[0]" eq '' ) {      

      if ( $opt < 2 )  {
	  print "\n WARNING($fnam): could not find keyword '${string}' \n"; 
	  print "\t in file  $inFile \n"; 
	  $| = 1;  # flush stdout
      }

#      print "\t grep-line = '@{line}' \n";
      if ( $opt == 0 ) { die " ***** ABORT ***** \n" ; }
  }

  return @tmpall ;

}  # end of parse_line 


# ==================================================
sub parse_array(@) {

    my($KEY, $NWD_PARSE, $OPT, @ARRAY ) = @_ ;

    # May 2014
    # analogous to parse_line, but search input $KEY in input @ARRAY
    # instead of grepping file.  Return NWD_PARSE values after $KEY.
    # Input OPT controls error handling for missing key:
    #   OPT = 0 -> abort
    #   OPT = 1 -> give warning
    #   OPT = 2 -> do nothing
    #
    # One subtle difference is that you cannot pass NWD_PARSE=999 
    # to grab everything to the end-of-line;
    # to get everything on one line, use parse_line.
    # In short, you MUST know NWD_PARSE.
    #
    # Note that
    #   $ARRAY[0] = "everything after 1st $KEY"
    #   $ARRAY[1] = "everything after 2nd $KEY"
    # etc ...
    # The calling function must break up each ARRAY line
    # to get separate strings.
    #
    # Sep 12, 2014: return NWD_PARSE words, or words up to ENDLINE 
    #
    # Mar 3 2019: stop reading past comment indicator '#'

    my $fnam = "parse_array";

    my( @tmp_all,  $NMATCH, @MSGERR );
    my (@INDX, $indx, $indx0, $indx1, $j );
    
    @tmp_all = (); # init output

    if ( scalar(@ARRAY) == 0 ) {
	$MSGERR[0] = "Cannot parse empty \@FILE_CONTENTS array" ;
	$MSGERR[1] = "for KEY=$KEY";
	sntools::FATAL_ERROR(@MSGERR);
    }

    @INDX  = grep { $ARRAY[$_] eq "$KEY" } 0 .. $#ARRAY ;
    
    if ( $NWD_PARSE == 0 ) {
	# just check for string match
	if ( scalar(@INDX) > 0 ) { return "MATCH" ; }
    }

    $NMATCH = 0 ;
    foreach $indx (@INDX) {
	# make sure we have an exact match and not an approximate match
	if ( $ARRAY[$indx] eq $KEY ) {
	    $indx0 = $indx + 1;  
	    $indx1 = $indx + $NWD_PARSE;

	    # don't go past ENDLINE
	    # Mar 3 2019: don't go past comment symbo, '#'
	    $j = $indx0;

	    while ( $ARRAY[$j] ne "ENDLINE"  && $j <= $indx1 ) {
		if ( substr($ARRAY[$j],0,1) eq '#' ) { goto DONE_LINEPARSE ; }
		{ $j++ ; }
	    }

	  DONE_LINEPARSE:
	    $indx1 = $j - 1;
	    $tmp_all[$NMATCH] = "@ARRAY[$indx0 .. $indx1]" ; 
	    $NMATCH++ ;
	}
    }

    if ( $NMATCH == 0 ) {
	if ( $OPT < 2 ) {
	    print "\n WARNING($fnam): could not find keyword '${KEY}' \n"; 
	    print "\t in the input \@ARRAY \n"; 
	    print "\t ARRAY = @ARRAY \n";
	    $| = 1;  # flush stdout
	}
	if ( $OPT == 0 )  { die " ***** ABORT ***** \n" ; }
    }

    return @tmp_all ;    

}  # end of parse_array

sub parse_value_after_key(@) {

    # Created Oct 25 2019
    # for inString = "KEY1 VAL1 KEY2 VAL2 KEY3 VAL3"
    # and key = KEY2, function returns VAL2.
    # Note that multiple values are returned if 
    # key appears multiple times.

    my($KEY,$INPUT_STRING) = @_ ;

    my @valueList = ();
    my ($indx, @INDX, @WDLIST, $value );
    @WDLIST   = split(/\s+/,$INPUT_STRING) ;
    @INDX = grep{$WDLIST[$_] eq $KEY} 0 .. $#WDLIST ;	

    foreach $indx (@INDX) {
	$value = $WDLIST[$indx+1] ;
	@valueList = (@valueList, $value);
    }

    return(@valueList);

} # end sub parse_value_after_key

# ===============================================
sub loadArray_fromFile(@) {
    my($inFile, $contents ) = @_ ;

    # May 2 2014
    # for input @inFile, read file contents and
    # return these contents in @contents.
    # Each line of infile is broken into separate strings,
    # and thus @contents is an array of strings, not lines.
    #
    # Sep 12, 2014: skip comment lines and insert ENDLINE strings

    my ($tmp, @MSGERR, $LINE, $str0) ;

    print " Read/store contents of $inFile \n" ;

    $inFile     =~ s/\s+$// ;  # remove trailing spaces
    if ( !(-e $inFile ) ) {
	$MSGERR[0] = "Cannot open input file" ;
	$MSGERR[1] = "$inFile ";
	$MSGERR[2] = "from loadArray_fromFile";
	sntools::FATAL_ERROR(@MSGERR);
    }
    open (PTR_INFILE_TEMP, $inFile);

    foreach $LINE ( <PTR_INFILE_TEMP> ) {

        $LINE  =~ s/\s+$// ;   # trim trailing whitespace
	if ( length($LINE) == 0 ) { next ; }
	# skip comment line
	$str0 = substr($LINE,0,1);	
	if ( $str0 eq '#' ) { next ; }
	if ( $str0 eq '!' ) { next ; }
	if ( $str0 eq '%' ) { next ; } 

	@$contents = ( @$contents, split(/\s+/, $LINE) ) ;
	@$contents = ( @$contents, "ENDLINE" ) ;
#	print " xxx LINE = '$LINE' \n";

    }

    # Aug 2020: abort on CONFIG key   
    my $OPT_QUIET = 2 ;
    my @FOUND_CONFIG = sntools::parse_array("CONFIG:", 1, $OPT_QUIET, 
					    @$contents);
    if ( scalar(@FOUND_CONFIG) > 0 ) {
	$MSGERR[0] = "$0" ;
	$MSGERR[1] = "cannot process refactored input file." ;
	$MSGERR[2] = "Use submit_batch_jobs.sh" ;
	sntools::FATAL_ERROR(@MSGERR);	
    }

    return ;
}  


# ===============================================
sub loadArray_excludeLines(@) {
    my($inFile, $excludeKey_start, $excludeKey_end, $contents ) = @_ ;

    # Feb 2019
    # Same as loadArray_fromFile, but exclude lines between
    # start end keys.

    my ($tmp, @MSGERR, $LINE, $str0, @wdlist) ;
    my $SKIPLINE=0;

    print " Read/store contents of $inFile \n" ;
#    $msgerr="ERROR: Cannot open '$inFile' \n from loadArray_excludeLines\n";
#    open PTR_INFILE_TEMP, $inFile or die "\n $msgerr\n ";

    $inFile     =~ s/\s+$// ;  # remove trailing spaces
    if ( !(-e $inFile) ) {
	$MSGERR[0] = "Cannot open input file" ;
	$MSGERR[1] = "$inFile ";
	$MSGERR[2] = "from loadArray_excludeLines";
	sntools::FATAL_ERROR(@MSGERR);
    }

    open (PTR_INFILE_TEMP, $inFile ) ; 

    foreach $LINE ( <PTR_INFILE_TEMP> ) {

        $LINE  =~ s/\s+$// ;   # trim trailing whitespace
	if ( length($LINE) == 0 ) { next ; }
	# skip comment line
	$str0 = substr($LINE,0,1);	
	if ( $str0 eq '#' ) { next ; }
	if ( $str0 eq '!' ) { next ; }
	if ( $str0 eq '%' ) { next ; } 

	@wdlist    = split(/\s+/, $LINE);
	if ( $wdlist[0] eq "$excludeKey_start" ) { $SKIPLINE=1; }
	if ( $wdlist[0] eq "$excludeKey_end"   ) { $SKIPLINE=0; }
	if ( $SKIPLINE ) { next; }

	@$contents = ( @$contents, @wdlist ) ;
	@$contents = ( @$contents, "ENDLINE" ) ;
    }

    return ;
}  # end loadArray_excludeLines


# ===============================================
sub Killjobs(@) {

# Created Jun 19, 2011 by R.Kessler
# Generic function to kill all jobs on @NODELIST

  my(@NODELIST) = @_;

  my ($NTMP, $node, $q, $cmdkill, $i, @NODELIST_SPLIT ) ;
# ssh into each node and do "kill -KILL -1"

  $q = "\"";
  $NTMP = scalar(@NODELIST);
    
  for ( $i=0; $i < $NTMP ; $i++ ) {
      @NODELIST_SPLIT    = split(/\s+/,$NODELIST[$i]) ;
      foreach $node ( @NODELIST_SPLIT ) {         
	  $cmdkill = "ssh -x $node ${q}kill -KILL -1${q}" ;
	  print "\t $cmdkill \n";
	  system("$cmdkill &");
      }
  }
  
  printf " Done killing @NODELIST. \n" ;


}  # end of KillJobs


# ================================================
sub FATAL_ERROR_STAMP(@) {

    # April 17 2019
    # write FAILURE message to DONE_STAMP file if file is specified;
    # then call FATAL ERROR.

    my($DONE_STAMP_FILE,@MSGERR) = @_ ; 

    if ( length($DONE_STAMP_FILE) > 2 ) 
    { qx(echo FAILURE > $DONE_STAMP_FILE); }
    
    sntools::FATAL_ERROR(@MSGERR);
} 

# ================================================
sub FATAL_ERROR(@) {
    
    my (@ERRMSG) = @_ ;

    my ($NERRMSG, $i);

    $NERRMSG = scalar(@ERRMSG) ;
    
    print "\n" ;
    print "\n" ;
    print "\n   `|```````|`    ";
    print "\n   <| o\\ /o |>    ";
    print "\n    | ' ; ' |     ";
    print "\n    |  ___  |     ABORT program on Fatal Error. " ;
    print "\n    | |' '| |     " ;
    print "\n    | `---' |     " ;
    print "\n    \\_______/    " ;
    print "\n" ;    

# xxx mark delete     print "\n ERROR MESSAGES : \n" ;
    print "\n FATAL ERROR ABORT : \n" ;
    for ( $i = 0; $i < $NERRMSG; $i++ ) {
	print "\t $ERRMSG[$i] \n";
    }
    die "\n\n" ;

}  # ERRMSG

# ======================
sub printArray(@) {

    my (@ARRAY) = @_ ;
    # print array values on one line

    my ($elem) ;
    foreach $elem ( @ARRAY ) {
        $elem  =~ s/\s+$// ;   # trim trailing whitespace
	print "$elem ";
    }
 
} # end of printArray

# ====================
sub sed_nmlInsert(@) {

    my ($nmlFile, $nml, $nmlVar, $nmlValue) = @_ ;

    # return sed txt-command to insert namelist variable
    # If variable already exists then return ''.
    # - nmlFile  = name of file
    # - nml      = SNLCINP or FITINP
    # - nmlVar   = name of variable (MXLC_PLOT, etc ...)
    # - nmlValue = value

    my (@tmp, $NTMP, $sedcmd, $q1 );

    $sedcmd = "" ;
    @tmp = qx(grep $nmlVar $nmlFile);
    $NTMP = scalar(@tmp);
    if ( $NTMP == 0 ) {
        $sedcmd = "-e '/\&${nml}/a\ $nmlVar = $nmlValue' " ;
    }

    return $sedcmd ;

} # end of sed_nmlInsert


# ====================
sub sed_strReplace(@) {

    my ($inFile, $outFile, $NREPLACE, @stringList ) = @_ ;

    # Create Sep 10, 2011 by R.Kessler
    #
    # Generic utility to replace a set of KEY-STRING with another
    # set of strings. This is useful to modify template files that
    # have specific KEY-strings to replace.
    #
    # (I) $inFile          : input file with key-strings to replace
    # (I) $outFile         : output file with key-strings replaced by ...
    # (I) $NREPLACE        : Number of keys to replace
    # (I) @stringList      : first $NREPLACE values are keys to replace
    #                      : Next  $NREPLACE values are substitute-strings
    # 
    # Comments
    # - The lengths of $stringList must by 2*$NREPLACE.
    # - use semicolon to split string to multiple lines. 
    #   i.e.,   "job1 arg1 ; job2 arg2 ; touch jobs.DONE"
    #   gets written out as
    #     job1 arg1
    #     job2 arg2
    #     touch jobs.DONE
    #
    # Feb 17, 2013: replace ' with $qq (double quote) in sed commands
    #               so that we can use things like "MAGSHIFT 'g .02' "
    #
    # ------------------------
    
    my ($sedcmd, $sedAdd, $N, $N2, @MSGERR, $i, $i2, $ierr, $key, $string );
    my ($isplit, @splitString, $tmpSplit, @tmp );
    my $qq = '"' ;
    my $q  = "\'" ;

    $N2 = scalar(@stringList);
    $N  = $N2/2;

    if ( $N != $NREPLACE || $N == 0 ) { 
	$ierr = -1 ; 
	$ierr++ ; @MSGERR[$ierr] = "Expected $NREPLACE strings to replace" ;
	$ierr++ ; @MSGERR[$ierr] = "but found only $N. ";
	for ( $i = 0; $i < $N2 ; $i++ ) {
	    $ierr++ ; @MSGERR[$ierr] = 
		" - stringList[$i] = $stringList[$i] ";
	}
	sntools::FATAL_ERROR(@MSGERR);
    }

    $sedcmd = "sed" ;
    for ( $i=0 ; $i < $N ; $i++ ) {
	$i2     = $i + $N ;
	$key    = "$stringList[$i]" ;
	$string = "$stringList[$i2]" ;

	# abort if this key does not exist
##	@tmp  = sntools::parse_line($inFile, 0, $key, 0 ) ;

	$sedAdd = "" ;

	# if string has a semicolon => split string on separate lines
	# for human readability
	$isplit  = index($string,";") ;
	if ( $isplit >= 0 ) {
	    @splitString = split(';', $string);
	    foreach $tmpSplit ( @splitString ) {
		if ( $tmpSplit ne "" ) {
		    $sedAdd = "$sedAdd " . 
			"-e ${qq}/${key}/a\ $tmpSplit\\n${qq}" ;
		}
	    }
	    # delete the original key since it was not physically replaced
	    $sedAdd = "$sedAdd " . "-e ${qq}/$key/d${qq}";
	}
	else {
	    $sedAdd = "-e ${qq}s/$key/$string/g${qq}";  # plain replace/substitute
	}

	$sedcmd = "$sedcmd " . "$sedAdd";
    }

    # print "xxx sedcmd: $sedcmd \n";

    # replace keys and make new batch file.
#    print "xxx in sntools::sed_strReplace \n";
#    print "xxx will sed $inFile into $outFile \n";
    qx($sedcmd $inFile > $outFile );


} # end of sed_strReplace



sub shellName {

    my ($tmpName, $itmp, $shell, @MSGERR);
    my @shellNameList = ( "bash", "csh" );

    $tmpName = $ENV{SHELL} ; 
    
    foreach $shell ( @shellNameList ) {
	$itmp   = index($tmpName,$shell) ;
	if ( $itmp >= 0 ) { return $shell ; }
    }

    #  if we get here then abort
    $MSGERR[0] = "ERROR in shellName";
    $MSGERR[1] = "Unrecognized shell: '$tmpName'";
    $MSGERR[2] = "Valid shells : @shellNameList";
    sntools::FATAL_ERROR(@MSGERR);

} #  end of shellName


sub setVarString(@)  {

    my ($SHELL, $NAME, $VALUE) = @_ ;
    
    my (@MSGERR);

    if ( $SHELL eq "csh" ) {
	return "set $NAME = `$VALUE` ";
    } 
    elsif ( $SHELL eq "bash" ) {
	return "$NAME=`$VALUE` ";
    }
    else {
	$MSGERR[0] = "SETENV: Unrecognized shell '$SHELL'";
	sntools::FATAL_ERROR(@MSGERR);
    }

} # end of setVarString


sub setEnvString(@)  {

    my ($SHELL,$ENV_NAME,$ENV_VALUE) = @_ ;

    my (@MSGERR);

    if ( $SHELL eq "csh" ) {
	return "setenv $ENV_NAME $ENV_VALUE" ;
    }
    elsif ( $SHELL eq "bash" ) {
	return "export $ENV_NAME=$ENV_VALUE" ;
    }
    else {
	$MSGERR[0] = "SETENV: Unrecognized shell '$SHELL'";
	sntools::FATAL_ERROR(@MSGERR);
    }

} # end of setEnvString


sub unsetEnvString(@)  {

    my ($SHELL,$ENV_NAME) = @_ ;

    my (@MSGERR);

    if ( $SHELL eq "csh" ) {
	return "unsetenv $ENV_NAME " ;
    }
    elsif ( $SHELL eq "bash" ) {
	return "unset $ENV_NAME " ;
    }
    else {
	$MSGERR[0] = "UNSETENV: Unrecognized shell '$SHELL'";
	sntools::FATAL_ERROR(@MSGERR);
    }

} # end of unsetEnvString


sub userPrompt(@) {
    my ($comment) = @_ ;
    my ($response) ;

    print " $comment \n" ;
    print " Enter y to continue ... ";
    $response = <STDIN> ;
    $response =~ s/\s+$// ;
    unless ( $response eq "y" ) {
	die " Bye bye.\n" ;
    }

} # end of userPrompt



sub get_namelistValue(@) {

    my ($nmlFile,$key) = @_ ;

    # Created Feb 01, 2012 by R.Kessler
    #    return value between quotes in
    #    key = 'value'
    #    Input $nmlFile is the name of the namelist file.
    #
    # Dec 6, 2012: if $VALUE = '', then check next line since
    #              fortrans allows this.
    # 
    # May 1 2013: 
    #  Fix but to allow blank spaces inside quotes. Bug interpreted ' '
    #  as a blank line and then read next line. Bug fix is to leave
    #  quotes on for initial read, then remove quotes at the end.
    # 
    # Sep 12, 2014: also grep '=' to avoid conflict with comments.
    #
    # Nov 07, 2014: check mutiple lines with nml key and remove 
    #               lines with comment in front. ABORT if >1 nml key.
    #

    my ($catcmd, $sedcmd, $VALUE, $jeq, $string, $LENV );
    my (@lineList, $line, $LINE, $c1, $NLINE) ;
    my $qq = '"' ;
    my $q  = "\'" ;

    $catcmd   = "cat $nmlFile | grep $key | grep '=' " ;    
    @lineList = qx { $catcmd }  ;

    $NLINE = 0 ;
    foreach $line ( @lineList ) {	
	my $c1 = substr($line,0,1);
	if ( $c1 eq '#' ) { next ; }
	if ( $c1 eq '%' ) { next ; }
	if ( $c1 eq '!' ) { next ; }

	$LINE = $line ;	$NLINE++ ;
    }

    if ( $NLINE > 1 ) {
	my @MSGERR ;
      	$MSGERR[0] = "Found $NLINE valid namelist keys for:";
	$MSGERR[1] = "    key = $key";
	$MSGERR[2] = "    nmlFile = $nmlFile ";
	$MSGERR[3] = "Allow only 0 or 1 valid nml keys";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $LINE     =~ s/\s+$// ;  # remove trailing spaces

    # VALUE is the string after the '=' sign
    $jeq    = index($LINE,"=");

    # make sure that the string before = exactly matches key
    $string  = substr($LINE,0,$jeq);
    $string  =~ s/\s+$// ;  # remove trailing spaces
    $string  =~ s/^\s+// ;  # remove leading spaces
    if ( $string ne $key ) { return "" ; }

    $VALUE  = substr($LINE, $jeq+1, 99);
    $VALUE  =~ s/\s+$// ;  # remove trailing spaces
    $VALUE  =~ s/^\s+// ;  # remove leading spaces

    $LENV = length($VALUE) ;
    if ( $LENV == 0 ) {
	$sedcmd = "sed -n '/$LINE/{n;p;}' $nmlFile" ;
	$VALUE  = qx($sedcmd);
	$VALUE  =~ s/\s+$// ;  # remove trailing spaces
	$VALUE  =~ s/^\s+// ;  # remove leading spaces
    }

    $VALUE  =~ s/$q//g ;   # remove quotes

#    print " xxx sedcmd = $sedcmd \n";
#    print " xxx line   = '$line' \n";
#    print " xxx $key   = |${VALUE}| \n";
#    die "\n xxx DIE OK xxx\n";

    return $VALUE ;

}  # end of get_namelistValue(@) {


# ===========
sub get_nmlValue_fromArray(@) {
    my ($KEY,$ARRAY) = @_ ;
    # analogous to get_namelistValue, but here we parse
    # the input array instead of the file.
    # ??? maybe later ???
} # end of get_nmlValue_fromArray

# =======================================
sub checkDirExist(@) {

    # check that $dir exists from the top directory;
    # this is to ensure that $dir can be found from 
    # a remote login. ABORT if $dir cannot be found.
    # Inputs:
    #  $dir      = name of directory to check.
    #  $comment  = comment to associate what this directory is;
    #              used only for abort comment.
    #

    my ($dir,$comment) = @_ ;

    my ($pwd_dir, @MSGERR );

    $pwd_dir = qx( cd ; cd $dir 2>/dev/null ; pwd ) ;
    $pwd_dir =~ s/\s+$// ;   # trim trailing whitespace

    if ( $dir ne $pwd_dir  ) {
	$MSGERR[0] = "Could not find directory =";
	$MSGERR[1] = "     '$dir'";
	$MSGERR[2] = "after cd ~/ ";
	$MSGERR[3] = "(instead found pwd = '$pwd_dir')";
	$MSGERR[4] = "Make sure that full path is specified for $comment ";
	sntools::FATAL_ERROR(@MSGERR);
    }
} # done checkDirExist

# ==============================
sub checkFileExist(@) {

    # check that $file exists
    # ABORT if $file cannot be found.
    # Inputs:
    #  $file        = name of directory to check.
    #  $comment     = comment to associate what this file is;
    #                 used only for abort comment.
    #

    my ($file,$comment) = @_ ;

    my (@MSGERR );

    unless ( -e $file  ) {
	$MSGERR[0] = "Could not find file =";
	$MSGERR[1] = "     '$file'";
	$MSGERR[2] = "Check performance of subroutine $comment ";
	sntools::FATAL_ERROR(@MSGERR);
    }
} # done checkFileExist


# ==============================
sub write_loginMacroLines(@) {

    my($FP, $OPT_JournalFont) = @_ ;

    # Oct 30, 2012
    # write macro lines that would normally go into .pawlogin,
    # and set default macro directory to $SNANA_DIR/kumacs
    # instead of ~/kumacs.

    my $fnam          = "write_MacroLines" ;
    my $SNANA_DIR     = $ENV{'SNANA_DIR'};
    my $pawMacroDir   = "${SNANA_DIR}/kumacs" ;

    print "   ${fnam}: use macros in  $pawMacroDir \n";

    print $FP  "macro/default '.,$pawMacroDir' -auto  \n";
    print $FP  "filecase keep \n";
    print $FP  "opt zfl \n";
    print $FP  "alias/create hpl 'exec hpl' c\n";
    print $FP  "alias/create hcopy 'exec hpl#hcopy' c \n";
    print $FP  "alias/create hfit 'exec hpl#hfit' c \n";
    print $FP  "alias/translation on \n";
    if ( $OPT_JournalFont ) { 
	print "   ${fnam}: use journal-quality fonts. \n";
	print $FP  "exec journal\n"; 
    }


    print $FP  "* ------------------------------------------- \n";
    print $FP  "\n";

} # end of write_loginMacroLines




#=========================================
sub make_batchFile(@) {
    
    # May 2013:
    # Create batch file for batch system such as torque (qsub) or sbatch.
    # Same as previous make_batchFile, but without doneFile or XCMD args.
    #
    # Aug 6 2018: add "cd $batchDir" before $batchJob command
    # Sep 12 2019: add batchName argument

    my ($BATCH_TEMPLATE, $batchDir, $batchName, $batchFile, $batchLog, 
	$batchMEM, $batchJOB ) = @_ ;

    my (@tmp, $inF, $outF, $KEY, @MSGERR );
    my ($NKEY, @REPLACE_KEY, @REPLACE_STRING );

    # convert BATCH_TEMPLATE into $batchFile by substituting
    # the REPLACE keys in the TEMPLATE

    # prepare strings to replace keys in batch-template file.

    $REPLACE_KEY[0] = "REPLACE_NAME" ;
    $REPLACE_KEY[1] = "REPLACE_LOGFILE" ;
    $REPLACE_KEY[2] = "REPLACE_MEM" ;
    $REPLACE_KEY[3] = "REPLACE_JOB" ;    

    $REPLACE_STRING[0] = "$batchName";    # xxx mark delete "$batchFile" ;
    $REPLACE_STRING[1] = "$batchLog" ;
    $REPLACE_STRING[2] = "$batchMEM";
    $REPLACE_STRING[3] = "cd $batchDir ; $batchJOB" ;

    $inF   = "$BATCH_TEMPLATE" ;
    $outF  = "${batchDir}/$batchFile" ;

    # make sure that each key exists
    foreach $KEY ( @REPLACE_KEY ) {
	@tmp = qx(grep $KEY $inF);
	if ( scalar(@tmp) <= 0 ) {
	    $MSGERR[0] = "Could not find required key ";
	    $MSGERR[1] = "'$KEY'" ;
	    $MSGERR[2] = "in $inF";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

    $NKEY = scalar(@REPLACE_KEY);
    sntools::sed_strReplace($inF,$outF, $NKEY,@REPLACE_KEY,@REPLACE_STRING);

} # end of make_batchFile


# ==========================================
sub rdfitres_getcol(@)
{
    # This code extracts the column of name "$namestring" from an input file.
    # Requires sntools.pm utility function rdfitres_util_checkhead(@)
    #
    # INPUTS:
    #
    # infile      : snana-formatted file to parse
    # p_arrtofill : pointer to data vector to be filled
    # namestring  : name of column to fetch from file
    # headstring  : string identifying header line (commonly "VARNAMES:")
    # datstring   : string identifying data line (commonly "SN:")
    # errcheck    : integer. If == 1, will trap column values == 0.0.
    #               this is most useful for reading in uncertainties. 
    # verbose     : (optional) integer turning comments OFF if == 0
    #
    # RETURNS:
    # Number of rows of data found, or -1 if $namestring not in header. 
    #
    # JLM Mar 2013
    #
    #
    # May 23, 2013 : altered rdfitres_util_checkhead subroutine. 
    #                subroutine header-variable name comparison is now
    #                case-insensitive. (JLM)
    #
    # Oct 02, 2013 : added errcheck flag to enable val==0 error trapping (JLM)
    #

    my ($infile, $p_arrtofill, $namestring, $headstring, 
	$datstring, $errcheck, $verbose) = (@_);

    my ($ngrid, @wdLine, @header, $line);
    my ($headcol, @MSGERR);

    if ( !defined($verbose)) { $verbose = 1; }
    
    if ( $verbose == 1 ) { 
	print "\n Getting $namestring data from file $infile ...\n\n";

	if ( $errcheck == 1 ) {
	    print "\n Will abort if column entries == 0 are found ... \n\n";
	}
    }

    #grab the header line
    @header = qx(grep $headstring $infile);

    #determine correct column
    $headcol = sntools::rdfitres_util_checkhead(\@header, $namestring);

    #fetch data and ngrid. return -1 if column is not found. 
    $ngrid = 0;
    if ( $headcol == -1 ) {
	$ngrid = $headcol;
    } else {
	open (INFILE, $infile);
	foreach $line (<INFILE>) {
	    chomp($line);
	    @wdLine = split(/\s+/,$line) ;
	    if ( $wdLine[0] eq $datstring ) {
		if ( $wdLine[$headcol] == 0 && $errcheck == 1 ) {
		    $MSGERR[0] = "ENTRY==0 FOUND IN COLUMN $headcol ($namestring)";
		    $MSGERR[1] = "OF FILE $infile";
		    sntools::FATAL_ERROR(@MSGERR);
		} else {
		    push(@$p_arrtofill, $wdLine[$headcol]);
		}
		$ngrid++;
	    }
	}
	close(INFILE);
    } 
    
    return ($ngrid);

} #end rdfitres_getcol(@)

# ==========================================

sub rdfitres_util_checkhead(@){

    # Utility function for sntools::rdfitres_getcol(@). 
    # NOT STAND ALONE FUNCTION!
    #
    # Given a header line and a variable name, 
    # this code determines and returns the column location. 
    #
    # INPUTS:
    #
    # headstring   : pointer to header string to parse
    # varname      : name of variable
    #
    # RETURNS:
    # 
    # key_column   : column location of varname
    #
    # May 23, 2013 : change header-variable name comparison to be
    #                case-insensitive (changing all to lower case). 
    #
    # ------------------------

    my ($p_headstring, $varname) = @_;

    my $return_var_column;
    my (@readline, $line_flag, $iline);
    my ($var_column, $key_column, $ncol);
    my (@MSGERR);
    my ($lowcase_readline, $lowcase_varname);

    (@readline) = split(/\s+/, @$p_headstring[0]) ;
    $ncol = @readline;
    $line_flag = 0;
    $iline = 0;
    $lowcase_varname = lc $varname;
    while ( $line_flag < 1 and $iline <= $ncol ) {
        $lowcase_readline = lc $readline[$iline];
##        print "xxx  In VAR loop: $iline, $readline[$iline], $varname\n";
##        print "xxx  In VAR loop: $lowcase_readline, $lowcase_varname \n";
        if ($lowcase_readline eq $lowcase_varname ) {
            $line_flag = $line_flag+1;
            $key_column = $iline;
##            print "xxx MATCH between $readline[$iline] and $varname \n";
        }
        $iline = $iline + 1;
    }

    unless ($line_flag == 1){
	$key_column = -1;
	###$MSGERR[0] = "COLUMN HEADER $varname NOT FOUND ";
	###sntools::FATAL_ERROR(@MSGERR);	
    }    
    return($key_column);

} #end rdfitres_util_checkhead(@)


# ==========================================

sub check_fileexist(@)
{
    my ($file, $comment) = @_;
    my (@MSGERR);

    unless (-e $file) {
	$MSGERR[0] = "Could not find file =";
	$MSGERR[1] = "          '$file'";
	$MSGERR[2] = "Check that file name for $comment ";
	$MSGERR[3] = " is correct";
	sntools::FATAL_ERROR(@MSGERR);
    }

} #check_fileexist (JLM)

# ==========================================
sub check_multipleJobs(@) {

    my ($JOBNAME) = @_ ;

    # April 2013
    # if running more than 2 $JOBNAME jobs simultaneously,
    # then give warning and ask user to continue. This is to 
    # catch left-over  background jobs.

    my $MOI  = `whoami` ;   $MOI =~ s/\s+$// ;   # trim trailing spaces

    my (@BLA, $NTOT, $N, $job, $response, @wdlist, $pid, $pidKill );

    # make list of all split_and_fit jobs and count them
    @BLA = qx(ps -aef | grep $MOI | grep $JOBNAME | grep perl );
    $NTOT   = scalar(@BLA) ;

    # expect two jobs (current + grep), so give warning
    # if there are more than 2 

    if ( $NTOT > 2 ) {
	print "\n" ;
	print " WARNING: detected multiple $JOBNAME jobs: \n\n";
	$N = 0 ;
	foreach $job (@BLA) {
	    my $j    = index($job, "grep" ) ;
	    @wdlist  = split(/\s+/,$job) ;
	    $pid     = $wdlist[1];
	    if ( $j >= 0 ) { next ; }  # skip internal grep
	    print "\t $job \n"; 
	    $N++ ;	    
	    if ( $N == 1 ) { $pidKill = $pid ; }
	}

#	print " Suggested action: kill $pidKill \n";
	print " Continue anyway y/[n] => " ;
	$response = <STDIN> ;
	$response =~ s/\s+$// ;
	unless ( $response eq "y" ) {
	    print " Bye bye. \n";
	    die " ***** ABORT ***** \n" ;
	}
    }

}  # end of  check_multipleJobs

# ==========================================

sub uniq_entries(@) {

#
# Created Sep 16, 2013 by J.Mosher
# Given pointer $p_all to vector @all and 
# pointer $p_unique to vector @unique, 
# fills @unique with unique entries of all.
# Use in lieu of List::MoreUtils qw/ uniq /
#

    my ($p_all, $p_unique) = (@_);
    
    my ($entry);
    my (%temp_hash) = () ;
    
    foreach my $entry (@ {$p_all} ) {
	$temp_hash{$entry} = 1;
    }

    @ {$p_unique} = keys(%temp_hash);

}

# ==========================================

sub weighted_avg(@) 
{

    my ($p_datcol, $p_errcol, $subname, $varname) = (@_);

    # Given pointers to columns containing
    # data and data uncertainties, 
    # this code calculates and returns 
    # the data weighted average and
    # average on the mean. 

    my ($func, @MSGERR, $idat);
    my ($ndat, $nerr, $tempdat, $temperr, $tempw);
    my ($sumdat, $sumw);
    my ($mean, $meane);

    $func = "sntools::weighted_avg";
    $ndat = scalar @{ $p_datcol };
    $nerr = scalar @{ $p_errcol };

    unless ( $ndat == $nerr ) {
	$MSGERR[0] = "In subroutine $subname, ";
	$MSGERR[1] = " function $func ";
	$MSGERR[2] = " datcol and errcol have different ";
	$MSGERR[3] = " lengths. Check inputs. ";
        sntools::FATAL_ERROR(@MSGERR);
    }

    if ( $ndat == 0 ) {
	$MSGERR[0] = "In subroutine $subname, ";
	$MSGERR[1] = " function $func ";
	$MSGERR[2] = " datcol $varname is empty.  ";
	$MSGERR[3] = " Check inputs. ";
        sntools::FATAL_ERROR(@MSGERR);
    }

    $sumdat = $sumw = 0;
    for ( $idat = 0; $idat < $ndat ; $idat++ ) {
	$tempdat = $$p_datcol[$idat];
	$temperr = $$p_errcol[$idat];
	$tempw = 1.0/($temperr*$temperr);
	$sumdat = $sumdat + $tempw*$tempdat;
	$sumw = $sumw + $tempw;
    }
    $mean = $sumdat/$sumw;
    $meane = 1.0/sqrt($sumw);

    #return ($mean, $meane, $ndat, $sumw);
    return ($mean, $meane, $ndat);

} # end weighted_avg(@) ;

# ==========================================

sub lininterp(@)
{
    my ($x0, $y0, $x1, $y1, $x) = (@_);

    my ($percent);
    my ($val);

    $percent = ($x-$x0)/($x1-$x0);
    $val = $percent*($y1-$y0) + $y0;

    return($val);
    
} # end lininterp(@)

# ==========================================

sub arith_avg(@) 
{

    my ($p_datcol, $subname, $datname) = (@_);

    # This code calculates the arithmetic average
    # and rms for the given data column. 

    my ($func, @MSGERR, $idat);
    my ($ndat, $tempdat, $temprms);
    my ($sumdat, $denom);
    my ($mean, $rms);

    $func = "sntools::arith_avg";
    $ndat = scalar @{ $p_datcol };

    unless ( $ndat > 1 ) {
	$MSGERR[0] = "In function $func ";
	$MSGERR[1] = "called from $subname ";
	$MSGERR[2] = "datcol $datname has <= 1 entry ";
	$MSGERR[3] = "Check inputs. ";
        sntools::FATAL_ERROR(@MSGERR);
    }

    $sumdat = $denom = 0;
    for ( $idat = 0; $idat < $ndat ; $idat++ ) {
	$tempdat = $$p_datcol[$idat];
	$sumdat = $sumdat + $tempdat;
	$denom++;
    }
    $mean = $sumdat/$denom;

    $sumdat = $denom = 0;
    for ( $idat = 0; $idat < $ndat ; $idat++ ) {
	$tempdat = $$p_datcol[$idat];
	$temprms = ($tempdat - $mean);
	$sumdat = $sumdat + $temprms*$temprms;
	$denom++;
    }
    $rms = sqrt($sumdat/($denom-1));


    return ($mean, $rms, $denom);

} # end arith_avg(@) 

# ==================================
sub extractStringOpt(@) {

    my ($string) = (@_) ;

    # for input string = "BLABLA(XYZ):"
    # returns  key=BLABLA: and arg=XYZ

    # set default output args
    my $key = $string;  
    my $arg = "";  # default with no ()

    my $j1  = index($string,"(") ;
    my $j2  = index($string,")") ;
    my $narg = ($j2 - $j1)-1 ;

    if ( $j1 <= 0 || $j2 <= 0  ) { return($key,$arg) ; }

    $key = substr($string,0,$j1) . substr($string,$j2+1,99);
    $arg = substr($string,$j1+1,$narg);	

    return($key,$arg) ;

} # end extractStringOpt

# ===============================================

sub clean_fileList(@) {

    my ($PROMPT_FLAG,@FLIST) = @_ ;

    my ($SIZE_TOT, $cmd, $cmd_find, $ftmp, $nline, @tmp);
    my( @wdlist, $size, $response );

    print "\n# ===================================== \n";
    $SIZE_TOT = 0.0 ;
    foreach $ftmp ( @FLIST ) {
	$cmd_find = "find . -name $ftmp -exec du -mc {} +";
	@tmp = qx($cmd_find);
	$nline    = scalar(@tmp);
	if ( $nline > 0 ) { print "@tmp\n" ; }
	@wdlist   = split(/\s+/,$tmp[$nline-1]) ;
	$size     = $wdlist[0];
	$SIZE_TOT += $size;
    }

    if ( $PROMPT_FLAG == 0 ) {
	$response = 'y' ;
    }
    else {
	print "\n Enter 'y' to delete $SIZE_TOT MB " . 
	    "of files listed above => \n  ";
	$response = <STDIN> ;  $response =~ s/\s+$// ;
    }

    if ( $response eq "y" ) {
	foreach $ftmp ( @FLIST ) {
	    print " Delete $ftmp files ... \n";
	    $cmd_find = "find . -name $ftmp -exec rm -rf {} +";
	    @tmp = qx($cmd_find);
	}
    }
    else {
        print " Files NOT deleted.\n" ;
    }    

} # end clean_fileList


# ===============================================

sub request_submit_batch_jobs(@) {
    my ($perl_codeName, $date_disable ) = @_ ;

    my $new_codeName = "submit_batch_jobs.sh";
    my(@msgerr);
    $msgerr[0] = "$perl_codeName" ;
    $msgerr[1] = "has been disabled since $date_disable ;" ;
    $msgerr[2] = "--> Use $new_codeName instead." ;
    $msgerr[3] = "For help on translating master input file," ;
    $msgerr[4] = "   $new_codeName -H TRANSLATE" ;
    sntools::FATAL_ERROR(@msgerr);
} 

1;
