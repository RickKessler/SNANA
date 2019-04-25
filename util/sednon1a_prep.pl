#!/usr/bin/perl
#
# TODO: fix "export" to use export or setenv based on $SHELL
#
# Created Dec 2011 by R.Kessler
# Prepare non1a spectral templates for the SNANA simulations.
#  
#
# Usage/instructions:
#   Select a directory outside of $SNDATA_ROOT
#
#   sednon1a_prep.pl <input directory> 
#      (where input directory contains SED.LIST with list of SED files)
#
#  options:
#    sednon1a_prep.pl <input dir> NOPREP     (skip prep and verify)
#    sednon1a_prep.pl <input dir> NOVERIFY   (do prep, skip verify)
#    sednon1a_prep.pl <input dir> PRIVATE    (use ~/kumacs if needed)
#    sednon1a_prep.pl <input dir> NOCLEAN    (leave SIM/$USER4* dirs)
#
#
# Three output directories are created along with a symbolic link,
#
#   $inDir_SNANA
#   $inDir_VERIFY
#   $inDir_PEAKMAG
#   non1a -> $inDir_SNANA  (sym link needed for snana sim)
# 
# and the NON1A-part of a sim-input file is created,
#   $inDir_SNANA/$inDir_NON1A.INPUT
#   $inDir_SNANA/$inDir_NON1A_RAW.INPUT (all MAGOFF,MAGSMEAR = 0)
#
#
# The former ($inDir_SNANA) is needed by snana to run the simulation.
# You can run it as an unofficial version by setting
# environment variable
#   setenv $SNANA_MODELPATH <directory containing $inDir_SNANA>
# To install these SEDs into snana, 
#    mv $inDir_SNANA  $SNDATA_ROOT/snsed/non1a
#
# Each time you run this script, the symbolic link (non1a)
# is updated; however, you can manually change this link
# back to other directories.
#
# The $inDir_VERIFY directory contains  "lcplots_VERIFY.pdf"
# which shows the original data, the smoothed light curve,
# and the simulated light curve based on the warped SED.
#
#
#       HISTORY
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Jan 06, 2012: read new format for SED.LIST to fix index for each SED
#               and to specify MAGOFF & MAGSMEAR for each SNTYPE
#
# Feb 07, 2012: 
#   simulate PEAKMAG for both RAW & nominal (SNANA-tuned) where RAW 
#   has all of the MAGOFF=0 and MAGSMEAR=0. Then overlay both on the 
#   same plot to see how much shift/smear is added to match observations.
#
# Feb 14, 2012: 
#    - in wrSEDfile_magoff(), fix bug to handle lines
#      with or without leading blanks
#    - PEAKMAG plot has B,V & R instead of only B band.
#    - fix to work with "NOPREP NOVERIFY" => remake sim-input include
#      file and PEAKMAG plots
#
# Mar 23, 2012: parse COMMENT: fields in SED.LIST and transfer these
#               user-comments to SIMGEN*.INPUT files.
#
# Oct 30, 2012:  call sntools::write_loginMacroLines(\*PTR_KUMAC,1);
#                 (beware: update is untested)
#
# Nov 14, 2012:  fix macro bug in make_PEAKMAG_plots
#
# Nov 15, 2012: add NOCLEAN flag to leave SIM/$user4* dirs,
#               otherwise clean up by default. See new logical $DO_CLEAN.
#
# Feb 16, 2013: for kcor, replace OUTPUT_GRID_HIS with OUTFILE,
#               and use fits format instead of hbook.
#               WARNING: not tested.
#
# ------------------

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# -------------------------------
# hard-wired global args

my $USER  = `whoami` ;    $USER  =~ s/\s+$// ; 
my $USER4 = substr($USER,0,4);

my $SNANA_SYMLINK     = "non1a" ;  # link to OUTPUT_SED_DIR
my $OUTPUT_SUFFIX     = "SNANA";
my $VERIFY_SUFFIX     = "VERIFY" ;
my $PEAKMAG_SUFFIX    = "PEAKMAG" ;
my $GENV_PREFIX       = "${USER4}_NON1APREP" ;
my $SEDLIST_FILENAME  = "SED.LIST";
my $INPUTS_LINK       = "INPUTS"  ;
my $SED_SNIa          = "Hsiao07.dat" ;
my $SNDATA_ROOT       = $ENV{'SNDATA_ROOT'};
my $HOST              = $ENV{'HOST'};
my $CADENCE_SIMLIB_VERIFY  =  2.0 ;   # cadence (days) for simlib
my $CADENCE_SIMLIB_PEAKMAG = 20.0 ;   # cadence (days) for simlib
my $SIM_PEAKMJD       = 55300.0 ;

my $OPTABORT = 0 ;
my $OPTWARN  = 1 ;
my $OPTQUIET = 2 ;
my $MAGOFF_DEFAULT =   0.0 ;
my $TREF_EXPLODE   = -19.0 ; # apply offset if first TREST is 0
my $DO_VERIFY  = 1 ;
my $DO_SEDPREP = 1 ;
my $DO_PEAKMAG = 1 ;
my $NSURVEY    = 0 ;
my $SIM_CID    = 777777 ; # for easy string replacement

my $USE_PRIVATE_KUMACS   = 0;

my $SNANA_MODELPATH      = `pwd` ;   $SNANA_MODELPATH  =~ s/\s+$// ;  
my $set_SNANA_MODELPATH  = "export SNANA_MODELPATH=$SNANA_MODELPATH";

my $FILTLIST_PEAKMAG     = "UBVRI";
my $FILTPATH_PEAKMAG     = "Bessell90" ;
my $NFILT_PEAKMAG        = length($FILTLIST_PEAKMAG);
my @SIMINP_FILE_PEAKMAG  = 
    ( "SIMGEN_PEAKMAG_RAW.INPUT" , "SIMGEN_PEAKMAG.INPUT" );
my @SIMLOG_FILE_PEAKMAG  = 
    ( "SIMGEN_PEAKMAG_RAW.LOG", "SIMGEN_PEAKMAG.LOG" );
my $SIMLIB_FILE_PEAKMAG  = "NON1A_PEAKMAG.SIMLIB" ;
my $KCOR_FILE_PEAKMAG    = "KCOR_PEAKMAG.fits" ;
my @GENVERSION_PEAKMAG   = 
    ( "${GENV_PREFIX}_PEAKMAG_RAW", "${GENV_PREFIX}_PEAKMAG" ) ;

my $GENPREFIX_PEAKMAG    = "PEAKMAG";
my @NGEN_PEAKMAG         = ( 200, 10000 ) ;
my @NTUP_FILE_PEAKMAG    = ("NON1A_PEAKMAG_RAW.TUP", "NON1A_PEAKMAG.TUP" );
my @PEAKMAG_COMMENT      = 
    (
     "Simulate RAW SEDs with all MAGOFF = MAGSMEAR = 0",
     "Simulate SNANA SEDs with nominal MAGOFF and MAGSMEAR" 
     );


my $DO_CLEAN = 1; # default is to remove SIM subdirs

# misc global args


my (@MSGERR, $INPUT_SED_DIR, $OUTPUT_SED_DIR, $VERIFY_DIR );
my (@SEDFILE_LIST_INPUT, @SEDFILE_LIST_OUTPUT, @SEDINDX_LIST );
my (@SEDNAME_LIST, @SEDMAGOFF_LIST, @SNTYPE_LIST, @IDSURVEY_LIST  );
my (@SURVEY_NAMES, @SURVEY_PRIMARY, @SURVEY_SIMLIB, @SURVEY_KCORFILE );
my (@SURVEY_FILTERFILES, @SURVEY_FILTERLIST );
my (@SURVEY_FILTERZP, @SURVEY_FILTERPATH);
my (@SURVEY_VERSION, @SURVEY_HISFILE );
my ($NSEDFILE,  @NFILT_SURVEY, @simFile_include, @SIMFILE_INCLUDE );

my ($NSEDVER, @SEDVER_GENVERSION, @SEDVER_GENPREFIX );
my (@SEDVER_SIMFILE_INPUT, @SEDVER_SIMFILE_LOG );
my (@SEDVER_IDSURVEY, @SEDVER_INDEX, @SEDVER_SNID, @SEDVER_REDSHIFT);


my ($NTYPE_USER, @ITYPE_USER, @WTYPE_USER, @NPERTYPE_USER, @CTYPE_USER );
my (@MAGOFF_USER, @MAGSMEAR_USER, @COMMENT_USER );
my (@SEDITYPE_LIST, @SEDWTYPE_LIST );

my ($PEAKMAG_DIR);

# ------------------------
sub  parse_args ;
sub  parse_SEDLIST ;
sub  parse2_SEDLIST ;
sub  parse_SEDheaders ;
sub  make_SNANA_outDir ;
sub  make_SEDfile ;
sub  make_photFile ;
sub  wrSEDfile_magoff ;
sub  make_SNANA_listFile ;
sub  make_SNANA_simInputFile ;
sub  clean_PREP ;

sub  VERIFY ;
sub  get_IDSURVEY ;
sub  parse_NON1Afile ;
sub  parse_WARPfile ;
sub  make_VERIFY_kcorFile ;
sub  make_VERIFY_SIMLIB ; 
sub  make_VERIFY_simInputFile ;
sub  run_VERIFY_sim ;
sub  merge_VERIFY_sims ;
sub  run_VERIFY_snana ;
sub  make_VERIFY_lcplots ;


sub  PEAKMAG ;
sub  make_PEAKMAG_kcorFile ;
sub  make_PEAKMAG_SIMLIB ;
sub  make_PEAKMAG_simInputFile ;
sub  run_PEAKMAG_sim ;
sub  make_PEAKMAG_plots ;

sub clean ;

# ========== BEGIN MAIN ===========

&parse_args();

&parse_SEDLIST();
&parse_SEDheaders();

if ( $DO_SEDPREP ) {

    # make output directory
    &make_SNANA_outDir();

    # create rest-frame photometry file for each SED
    print "\n Create SED & photometry files for \n" ;
    my $ised;
    for ( $ised=1; $ised <= $NSEDFILE; $ised++ ) {
	print "$SEDNAME_LIST[$ised] ";
	&make_SEDfile($ised);  
	&make_photFile($ised);
    }
    printf("\n");
    
    # create final list file for SNANA
    &make_SNANA_listFile();

    &clean_PREP();

}  # end of DO_PREP if block


if ( $DO_SEDPREP  || $DO_PEAKMAG ) {
    # create sample sim-input with NON1A keys
    &parse2_SEDLIST();
    &make_SNANA_simInputFile(0);
    &make_SNANA_simInputFile(1);
}


if ( $DO_VERIFY  )  { &VERIFY(); }

if ( $DO_PEAKMAG )  { &PEAKMAG(); }

&clean();

# ========= END MAIN ==============

# -----------------------------------------------
sub parse_args {
    my ($NARG, $i);
    $NARG          = scalar(@ARGV);

    if ( $NARG < 1 ) {
	$MSGERR[0] = "Must give input SED directory as argument.";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $INPUT_SED_DIR  = $ARGV[0];

    # make sure that there are no slashes in SED_DIR
    my $islash = index($INPUT_SED_DIR,"/");
    if ( $islash > 0 ) {
	$MSGERR[0] = "Invalid '/' found in directory name" ;
	$MSGERR[1] = "'$INPUT_SED_DIR'" ;
	$MSGERR[2] = "Remove slashes.";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $OUTPUT_SED_DIR = "${INPUT_SED_DIR}_${OUTPUT_SUFFIX}";
    $VERIFY_DIR     = "${INPUT_SED_DIR}_${VERIFY_SUFFIX}";
    $PEAKMAG_DIR    = "${INPUT_SED_DIR}_${PEAKMAG_SUFFIX}";

    $simFile_include[0] = "SIMGEN_INCLUDE_NON1A_RAW.INPUT" ;
    $simFile_include[1] = "SIMGEN_INCLUDE_NON1A.INPUT" ;
    my $tmpDir = "$SNANA_MODELPATH/$OUTPUT_SED_DIR";
    $SIMFILE_INCLUDE[0] = "$tmpDir/$simFile_include[0]";
    $SIMFILE_INCLUDE[1] = "$tmpDir/$simFile_include[1]";


    for ( $i=1; $i <= $NARG; $i++ ) {
	if ( $ARGV[$i] eq "NOPREP"    ) { $DO_SEDPREP = 0 ; }
	if ( $ARGV[$i] eq "NOVERIFY"  ) { $DO_VERIFY  = 0 ; }
	if ( $ARGV[$i] eq "NOPEAKMAG" ) { $DO_PEAKMAG = 0 ; }
	if ( $ARGV[$i] eq "PRIVATE"   ) { $USE_PRIVATE_KUMACS = 1 ; }
	if ( $ARGV[$i] eq "NOCLEAN"   ) { $DO_CLEAN   = 0 ; }
    }


    if ( !(-e $INPUT_SED_DIR) ) {
	$MSGERR[0] = "Could not find input SED-directory";
	$MSGERR[1] = "'$INPUT_SED_DIR'";
	sntools::FATAL_ERROR(@MSGERR);
    }

    print " DO_SEDPREP:  $DO_SEDPREP  (translate SEDs into SNANA format) \n";
    print " DO_VERIFY:   $DO_VERIFY   (run sim to verify each input LC) \n";
    print " DO_PEAKMAG:  $DO_PEAKMAG  (mean and RMS of peak B band) \n";
    print " DO_CLEAN:    $DO_CLEAN    (rm -r SIM subdirs) \n";

} # end of parse_args


# ============================
sub  make_SNANA_outDir {

    if ( -d $OUTPUT_SED_DIR) { qx(rm -r $OUTPUT_SED_DIR); }
    qx(mkdir $OUTPUT_SED_DIR);

    # now make symbolic link to non1a ; needed for snana simulation.
    if ( -e $SNANA_SYMLINK ) { qx(rm  $SNANA_SYMLINK); }
    qx(ln -s $OUTPUT_SED_DIR $SNANA_SYMLINK);

} # end of     &make_SNANA_outDir


# ------------------------------
sub parse_SEDLIST {

    my (@tmp, $LISTFILE, $key, $ised, $sedIndx, $sedFile, @wdlist ) ;

    $LISTFILE = "$INPUT_SED_DIR/$SEDLIST_FILENAME" ;
    if ( !(-e $LISTFILE) ) {
	$MSGERR[0] = "Could not find SED-list file";
	$MSGERR[1] = "$LISTFILE";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $key   =  "SED:" ;
    @tmp   =  sntools::parse_line($LISTFILE, 2, $key, $OPTABORT );
    $NSEDFILE = scalar(@tmp);


    for ( $ised=1; $ised <= $NSEDFILE; $ised++ ) {
	@wdlist   = split(/\s+/,$tmp[$ised-1] ) ;
	$sedIndx     = $wdlist[0] ;
	$sedFile     = $wdlist[1] ;

	$SEDFILE_LIST_INPUT[$ised] = $sedFile ;
	$SEDINDX_LIST[$ised]       = $sedIndx ;
    }

    print " Found $NSEDFILE SED files to process. \n";


} # end of parse_SEDLIST


# ---------------------------
sub parse2_SEDLIST {

    # parse SNTYPE key(s) from SEDLIST ;
    # needed to make sample sim-input file with NON1A keys.
    # Each sntype key is as follows:
    #  SNTYPE:  <ITYPE>  <WGTSUM>  <MAGOFF>  <MAGSMEAR>  <CTYPE-LIST>
    #
    # where ITYPE is an integer type, and CTYPE-LIST is a space-separated
    # list of types to include such as 'II IIP IIp'  or 'IIl IIL'
    # (but without single quotes).
    # Note that the resulting file is NOT verification, 
    # but for users to simulate NON1A.

    my ($LISTFILE, $key, @tmp, $i, $ii, $NTMP, @wdlist, $ised );
    my ($ITYPE, $WTYPE, $CTYPE1, $CTYPE2, @iList, $name, $N );
    my ($MOFF, $MSMEAR);

    # init user ITYPE and WGT to 0 for each SED
    for ( $ised=0; $ised <= $NSEDFILE; $ised++ ) {
	$SEDITYPE_LIST[$ised] = 0 ;
	$SEDWTYPE_LIST[$ised] = 0.0 ;
    }
    $ITYPE_USER[0] = 0   ;
    $WTYPE_USER[0] = 0.0 ;
    $MAGOFF_USER[0]   = 0.0 ;
    $MAGSMEAR_USER[0] = 0.0 ;

    $LISTFILE = "$INPUT_SED_DIR/$SEDLIST_FILENAME" ;
    $key   =  "SNTYPE:" ;
    @tmp   =  sntools::parse_line($LISTFILE, 99, $key, $OPTQUIET );

    $key      =  "COMMENT:" ;
    @COMMENT_USER  =  sntools::parse_line($LISTFILE, 99, $key, $OPTQUIET );

    $NTYPE_USER = scalar(@tmp);
    if ( $NTYPE_USER == 0 ) { return ; }

    for ( $i=1; $i <= $NTYPE_USER; $i++ ) {

	@wdlist  = split(/\s+/,$tmp[$i-1] ) ;

	$ITYPE   = $wdlist[0];
	$WTYPE   = $wdlist[1];
	$MOFF    = $wdlist[2];
	$MSMEAR  = $wdlist[3];

	$NPERTYPE_USER[$i]   = 0 ;
	$ITYPE_USER[$i]      = $ITYPE ;
	$WTYPE_USER[$i]      = $WTYPE ;
	$MAGOFF_USER[$i]     = $MOFF ;
	$MAGSMEAR_USER[$i]   = $MSMEAR ;

	$NTMP            = scalar(@wdlist) - 4 ;
	$CTYPE_USER[$i]  = "$wdlist[4]";
	
	for ( $ii = 0; $ii < $NTMP; $ii++ ) {
	    $CTYPE1  = $wdlist[4+$ii] ;
	    $CTYPE1  =~ s/\s+$// ;  # remove trailing blanks

	    if ( $ii > 0 ) 
	    {  $CTYPE_USER[$i]   =  "$CTYPE_USER[$i]"  . "+$CTYPE1"; }
	    
	    for ( $ised=1; $ised <= $NSEDFILE; $ised++ ) {
		$CTYPE2 = $SNTYPE_LIST[$ised] ;
		if ( $CTYPE1 eq $CTYPE2 ) {
		    $iList[$ised]         = $i ;
		    $SEDITYPE_LIST[$ised] = $ITYPE ;
		    $NPERTYPE_USER[$i]++ ;
		}
	    }
	}
    }


    # set the weights now that we have NPERTYPE

    for ( $ised=1; $ised <= $NSEDFILE; $ised++ ) {
	$i = $iList[$ised];
	$N = $NPERTYPE_USER[$i] ;
	if ( $N == 0 ) { next ; }
	$WTYPE = $WTYPE_USER[$i]/$N ;
	$SEDWTYPE_LIST[$ised] = $WTYPE ;
    }


} # end of parse2_SEDLIST

# -------------------------------
sub parse_SEDheaders {

    # read list in list file and then parse each 
    # SED headers and store info.

    my ($inFile, $outFile, $ised, $c1, $key, $name, $magoff, $sntype );
    my ($inFile_full, $outFile_full);
    my (@tmp, $itmp);

    print "\n Get header info for: \n";

    for ( $ised=1; $ised <= $NSEDFILE; $ised++ ) {

	$inFile       =  $SEDFILE_LIST_INPUT[$ised] ;       
	$inFile_full  = "$INPUT_SED_DIR/$inFile" ;

	$key   =  "NAME:" ;
	@tmp   =  sntools::parse_line($inFile_full, 1, $key, $OPTQUIET );
	if ( scalar(@tmp) > 0 ) {
	    $name = $tmp[0];
	}
	else {
	    $itmp = index($inFile,".sed");
	    $name = substr($inFile,0,$itmp);
	}

	# get required SN TYPE
	$key   =  "SNTYPE:" ;
	@tmp   =  sntools::parse_line($inFile_full, 1, $key, $OPTABORT );
	$sntype = $tmp[0] ;

	# get optional mag offset
	$key   =  "MAGOFF:" ;
	@tmp   =  sntools::parse_line($inFile_full, 1, $key, $OPTQUIET );
	if ( scalar(@tmp) > 0 ) {
	    $magoff = $tmp[0];
	}
	else {
	    $magoff = $MAGOFF_DEFAULT ;
	}

	$outFile = "${name}.SED";
#	$outFile_full = "$OUTPUT_SED_DIR/$outFile" ;
	
	$SEDFILE_LIST_INPUT[$ised]  = $inFile ;
	$SEDFILE_LIST_OUTPUT[$ised] = $outFile ;
	$SEDNAME_LIST[$ised]        = $name ;
	$SEDMAGOFF_LIST[$ised]      = $magoff ;
	$SNTYPE_LIST[$ised]         = $sntype ;

	print "$name ";
    }

        
    print "\n Done parsing $NSEDFILE SED headers for SNANA simulation. \n";

} # end  of  parse_SEDheaders 


# ----------------------
sub make_SEDfile {

    # copy SED file if MAGOFF = 0;
    # otherwise re-write SED file with MAGOFF applied.

    my ($ised) = @_ ;

    my ($name, $magoff);
    my ($inFile, $inFile_full, $outFile, $outFile_full);

    $name   = $SEDNAME_LIST[$ised];
    $inFile = $SEDFILE_LIST_INPUT[$ised];
    $magoff = $SEDMAGOFF_LIST[$ised] ;
 
    $outFile      = "${name}.SED";
    $outFile_full = "$OUTPUT_SED_DIR/$outFile" ;
    $inFile_full  = "$INPUT_SED_DIR/$inFile";

    
    if ( $magoff == 0.0 ) {
	qx(cp $inFile_full $outFile_full); 
    }
    else {
	&wrSEDfile_magoff( $inFile_full, $magoff, $outFile_full );
    }

} # end of make_SEDfile


# --------------------------------------
sub wrSEDfile_magoff {

    # write SED file with MAGOFF applied to each flux.
    # Keep header info, but get rid of MAGOFF key.
    #
    # Feb 14, 2012: fix bug to handle lines with or without leading blanks
    #

    my ($inFile, $magoff, $outFile, $tmpFile ) = @_ ;

    # local variables

    my (@CONTENTS, $line, $fluxScale, $arg, @wdlist, $c0 );
    my ($Trest, $lam, $flux, $FLUX, $FIRST, $Trest_off );
    my ($cTrest, $cFLUX);

    @CONTENTS = `cat $inFile` ;
    $arg       = -0.4 * $magoff ;
    $fluxScale = 10.0**$arg;
    $tmpFile   = "${outFile}_TEMP";

#    print "\t\t fluxScale =  $fluxScale ... \n";

    open  PTR_SED , "> $tmpFile" ;

    $FIRST = 1;
    $Trest_off = 0.0 ;

    foreach $line ( @CONTENTS ) {
	$line  =~ s/\s+$// ;  # remove trailing blanks
	$line  =~ s/^\s+// ;  # remove leading spaces

	@wdlist     = split(/\s+/,$line) ;

	$c0  = substr($wdlist[0],0,1) ;
	if ( $c0 eq '#' ) {
	    print PTR_SED "$line \n" ;
	    next ;
	}
	
	$Trest = $wdlist[0];
	$lam   = $wdlist[1];
	$flux  = $wdlist[2];

	if  ( $FIRST && $Trest == 0.0 ) {
	    $Trest_off = $TREF_EXPLODE ;
	}

	$cTrest = sprintf("%9.4f",   $Trest + $Trest_off);
	$cFLUX  = sprintf("%10.4le", $flux * $fluxScale);

	print PTR_SED "$cTrest  $lam  $cFLUX \n" ;

	$FIRST = 0 ;
    }

    close PTR_SED ;

    # now strip out MAGOFF key and dump output to $outFile
    my $sedcmd  = "sed -e '/MAGOFF/d' ";
    qx($sedcmd $tmpFile > $outFile );

} # end of wrSEDfile_magoff {

# -------------------------------------
sub make_photFile {

    my ($ised) = @_ ;

    # Create photometry file from SED using kcor.exe

    my ($kcorFile, $photFile, $sedFile, $name, $sntype, $magoff );

    $name     = "$SEDNAME_LIST[$ised]" ;
    $kcorFile = "kcor_${name}.input";
    $photFile = "${name}.DAT";
    $sedFile  = "$SEDFILE_LIST_INPUT[$ised]";
    $sntype   = "$SNTYPE_LIST[$ised]";
    $magoff   = $SEDMAGOFF_LIST[$ised];

    open  PTR_KCOR , "> $OUTPUT_SED_DIR/$kcorFile";
    
    print PTR_KCOR "DUMP_SNMAG:   2 \n";
    print PTR_KCOR "SN_SED:       ../$INPUT_SED_DIR/$sedFile \n";
    print PTR_KCOR "SN_TYPE:      $sntype \n";
    print PTR_KCOR "SN_MAGOFF:    $magoff \n";
    print PTR_KCOR "TREST_RANGE:  -40  100 \n";
    print PTR_KCOR "\n";

    # use SDSS filters to define mags 
    print PTR_KCOR "MAGSYSTEM:   AB \n";
    print PTR_KCOR "FILTSYSTEM:  COUNT \n" ;
    print PTR_KCOR "FILTPATH:    SDSS \n"  ;
    print PTR_KCOR "FILTER:  a  u.dat  0.0 \n";
    print PTR_KCOR "FILTER:  b  g.dat  0.0 \n";
    print PTR_KCOR "FILTER:  c  r.dat  0.0 \n";
    print PTR_KCOR "FILTER:  d  i.dat  0.0 \n";
    print PTR_KCOR "FILTER:  e  z.dat  0.0 \n";
    print PTR_KCOR "\n";

    close PTR_KCOR ;

    my $kk = "kcor.exe $kcorFile"; 
    my $mv = "mv SNMAG_SDSS.TXT  ${name}.DAT" ;
    qx(cd $OUTPUT_SED_DIR ; $kk ; $mv );

} # end of make_photFile



sub make_SNANA_listFile {

    # create NON1A.LIST file for SNANA.

    my ($listFile, $ised, $name, $i3indx, $sntype );

    $listFile = "$OUTPUT_SED_DIR/NON1A.LIST" ;
    open  PTR_LIST , "> $listFile";
 
    print "\n Create final list file : $listFile \n";

    for ( $ised=1; $ised <= $NSEDFILE; $ised++ ) {
	$i3indx    = sprintf("%3.3d", $SEDINDX_LIST[$ised] );
	$name      = "$SEDNAME_LIST[$ised]" ;
	$sntype    = sprintf("%-4s" ,$SNTYPE_LIST[$ised] );
	print PTR_LIST "NON1A:  $i3indx  $sntype   $name \n";
    }

    close PTR_LIST ;

} # end of &make_SNANA_listFile


# ------------------------------
sub  make_SNANA_simInputFile {

    my ($OPT) = @_ ;

    # create sample sim-input file with NON1A keys.
    # This file should be included in a sim-input file.
    # OPT = 0 => all MAGOFF=  and MAGSMEAR=0 to get RAW distributions
    # OPT = 1 => use requested MAGOFF and MAGSMEAR values.

    my ($ised, $c3indx, $sntype, $name, $cdum, $comment );
    my ($ITYPE, $WTYPE, $CTYPE, $CWGT, $i, $NPRINT );
    my ($MOFF, $MSMEAR, $comment);

    
    if ( -e $SIMFILE_INCLUDE[$OPT] ) { qx(rm $SIMFILE_INCLUDE[$OPT]); }

    print "\n Create sample sim-input file : $simFile_include[$OPT] \n";

    open  PTR_SIMFILE , "> $SIMFILE_INCLUDE[$OPT]" ;

    # comments at top
    my $date   = `date +%F`;
    $date    =~ s/\s+$// ;  # remove trailing blanks
    print PTR_SIMFILE "# ----------------------------------------- \n";
    print PTR_SIMFILE "# Auto-generated comments: \n" ;
    print PTR_SIMFILE "#   NON1A keys for sim-input file. \n" ;
    print PTR_SIMFILE "#   Created $date   by  $USER    HOST=$HOST   \n";
    print PTR_SIMFILE "#   SCRIPT     : $0  $INPUT_SED_DIR \n";
    print PTR_SIMFILE "#   DIRECTORY  : $SNANA_MODELPATH  \n";
    print PTR_SIMFILE "#   NSED(NON1A): $NSEDFILE \n";
    print PTR_SIMFILE "# \n";


    print PTR_SIMFILE "# User-generated comments: \n" ;
    foreach $comment (@COMMENT_USER) {
	$comment  =~ s/\s+$// ;  # remove trailing blanks
	print PTR_SIMFILE "#   $comment \n";
    }

    print PTR_SIMFILE "\n";

    # ----------------

    # write user-SNTYPE keys
    for ( $i=1; $i <= $NTYPE_USER; $i++ ) {
	$ITYPE   = sprintf("%2.2d", $ITYPE_USER[$i] ) ;
	$CTYPE   = sprintf("%-10s", $CTYPE_USER[$i] ) ;
	$CWGT    = sprintf("(WGTSUM=%6.3f)", $WTYPE_USER[$i] ) ;
	print PTR_SIMFILE "#SNTYPE  $ITYPE  => $CTYPE   $CWGT \n" ;
    }
    print PTR_SIMFILE "\n\n";

    print PTR_SIMFILE "NON1A_KEYS: 5 \n";
    print PTR_SIMFILE "         INDEX   WGT    MAGOFF   MAGSMEAR  SNTYPE \n";

    $NPRINT = 0;

    # write templates sorted by the user-type

    for ( $i=0; $i <= $NTYPE_USER; $i++ ) {

	if ( $OPT == 0 ) {
	    $MOFF    = sprintf("%5.3f" , 0.0   ) ;
	    $MSMEAR  = sprintf("%5.3f" , 0.0   ) ;
	}
	else {
	    $MOFF    = sprintf("%5.3f" , $MAGOFF_USER[$i]   ) ;
	    $MSMEAR  = sprintf("%5.3f" , $MAGSMEAR_USER[$i]   ) ;
	}

	for ( $ised=1; $ised <= $NSEDFILE; $ised++ ) {

	    $ITYPE   = sprintf("%2.2d", $SEDITYPE_LIST[$ised] ) ;

	    if ( $ITYPE != $ITYPE_USER[$i] ) { next ; }
	    $WTYPE   = sprintf("%6.4f", $SEDWTYPE_LIST[$ised] ) ;
	    $CTYPE   = sprintf("%-4s" , $SNTYPE_LIST[$ised]   ) ;
	    $c3indx  = sprintf("%3.3d", $SEDINDX_LIST[$ised] );
	    $name    = "$SEDNAME_LIST[$ised]" ;
	    
	    $cdum    = "$WTYPE   $MOFF     $MSMEAR     $ITYPE" ; 
	    $comment = "# $CTYPE ($name)";
	    print PTR_SIMFILE "NON1A:    $c3indx   $cdum    $comment\n" ;
	    $NPRINT++ ;
	}
	print PTR_SIMFILE "\n";
    }

    close PTR_SIMFILE ;
    
    if ( $NPRINT != $NSEDFILE ) {
	$MSGERR[0] = "Wrote $NPRINT  NON1A keys";
	$MSGERR[1] = "but expected to write NSEDFILE = $NSEDFILE";
	sntools::FATAL_ERROR(@MSGERR);
    }


    $| = 1;  # flush stdout

} # end of     &make_SNANA_simInputFile

# ----------------------------
sub clean_PREP {
    # remove temp kcor*input files
    qx(cd $OUTPUT_SED_DIR ; rm kcor*.input *TEMP );
} 


# ----------------------------
sub VERIFY {

    my ($ID, @INFILES_WARP, $warpFile, $SURVEY, $ised, $dirList );
    # ----------------------------------------------

    print "\n" ;
    print "# ============================================= \n";
    print "  Prepare to VERIFY NON1A SEDs in simulations. \n" ;
    $| = 1;  # flush stdout

    # make local VERFIY_DIR  for kcor/simlib/sim-input files
    if ( -d $VERIFY_DIR ) { qx(rm -r $VERIFY_DIR); }
    qx(mkdir $VERIFY_DIR);

    # make symbolic link to input directory
    qx(cd $VERIFY_DIR; ln -s ../$INPUT_SED_DIR  $INPUTS_LINK );

    # remove old VERIFY versions to avoid confusion
    $dirList = "$SNDATA_ROOT/SIM/${GENV_PREFIX}_VERIFY*" ;
    qx(rm -r $dirList 2>/dev/null );

    &parse_NON1Afile() ;

    # ------------------------------------
    # now read/parse the WARP*input files.
    print "\n" ;
    @INFILES_WARP = `ls $INPUT_SED_DIR/WARP*` ;    
    foreach $warpFile ( @INFILES_WARP ) {
	$warpFile  =~ s/\s+$// ;  # remove trailing blanks
	&parse_WARPfile($warpFile);
    }

    # ------------------------------------
    # create kcor & simlib for each survey.
    print "\n";
    for ($ID = 1; $ID <= $NSURVEY; $ID++ ) { &make_VERIFY_kcorFile($ID); }
    print "\n";
    for ($ID = 1; $ID <= $NSURVEY; $ID++ ) { &make_VERIFY_SIMLIB($ID)  ; }
    print "\n";


    # make sim-input file for each SED template and run simulation
    print " Create/run simulation for: \n";
    for ($ised = 1; $ised <= $NSEDVER; $ised++ ) {
	&make_VERIFY_simInputFile($ised) ; 
	&run_VERIFY_sim($ised);
    }
    qx(chmod o-w $SNDATA_ROOT/SIM/${GENV_PREFIX}_VERIFY*);
    print "\n";

    print "\n";
    # merge sims separately for each survey,
    # then run snana job to get light curve plots
    for ( $ID=1; $ID <= $NSURVEY; $ID++ ) { 
	&merge_VERIFY_sims($ID); 
	&run_VERIFY_snana($ID);
    }

    # finally make light curve plots with overlays.
    &make_VERIFY_lcplots();


} # end of VERIFY


# ---------------------------------
sub get_IDSURVEY {

    # return integer id for this $SURVEY.
    # OPT = +1 => add new survey to list if not already on the list.
    # OPT = -1 => return -9 if survey is not on the list.


    my ($SURVEY, $OPT) =  @_ ;

    my ($ID, $survey_name);
    # check current surveys for match

    for ( $ID=1; $ID <= $NSURVEY; $ID++ ) {
	$survey_name = $SURVEY_NAMES[$ID];
	if ( $survey_name eq $SURVEY ) {
	    return $ID ;
	}
    }

    # if we get here then add a new survey to the list
    # and increment NSURVEY

    if ( $OPT < 0 ) { return -9 ; }


    $NSURVEY++ ;
    $SURVEY_NAMES[$NSURVEY] = $SURVEY ;
    return $NSURVEY ;


} # end of IDSURVEY

# --------------------
sub parse_NON1Afile {

    my ($ised, $listFile, @NON1A_LIST, $key, $tmp, @wdlist );
    my ($indx, $name, $c6, $itmp, @tmp, $TMP );
    my ($SURVEY, $SNID, $sedFile, $ZSN, $ID );

    $NSEDVER = 0;
    $listFile = "$OUTPUT_SED_DIR/NON1A.LIST" ;

    $key       =  "NON1A:" ;
    @NON1A_LIST =  sntools::parse_line($listFile, 3, $key, $OPTABORT );

    foreach $TMP (@NON1A_LIST) {
	@wdlist  = split(/\s+/,$TMP ) ;
	$indx    = $wdlist[0];
	$name    = $wdlist[2];	
	$c6      = substr($name,0,6) ;

	# skip the  Nugent ones.
	if ( $c6 eq "Nugent" ) { next ; }
	if ( $c6 eq "NUGENT" ) { next ; }

	# break the internal name into SURVEY and SNID
	$itmp   = index($name,"-");
	$SURVEY = substr($name,0,$itmp);
	$SNID   = substr($name,$itmp+1,50);

	# grep out the redshift from the SED file
	$key       =  "REDSHIFT:" ;
	$sedFile   = "$OUTPUT_SED_DIR/${name}.SED" ;
	@tmp       =  sntools::parse_line($sedFile, 1, $key, $OPTABORT );
	$ZSN       = $tmp[0] ;

	# get the survey ID 
	$ID = &get_IDSURVEY($SURVEY,1);

	print "\t Reading NON1A($indx) :  $SURVEY-$SNID  z=$ZSN \n";

	# store info
	$NSEDVER++ ;
	$SEDVER_IDSURVEY[$NSEDVER] = $ID;
	$SEDVER_SNID[$NSEDVER]     = $SNID;
	$SEDVER_REDSHIFT[$NSEDVER] = $ZSN ;
	$SEDVER_INDEX[$NSEDVER]    = $indx ;

    }

    if ( $NSURVEY <= 0 ) {
	$MSGERR[0] = "NSURVEY ${NSURVEY}" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    $| = 1;  # flush stdout

} # end of parse_NON1Afile 

# -----------------------------
sub parse_WARPfile {

    # parse WARP-input file for filters and primary ref that
    # are needed to make kcor-input, simlib, and sim-input files
    # for verification.

    my ($warpFile) = @_ ;

    my ($key, @tmp, $SURVEY, $IDSURVEY, @wdlist, $line);
    my ($f, $ffile, $zp, $last, $itmp, $NFILT );

    # start with the survey
    $key   = "SURVEY:" ;
    @tmp      =  sntools::parse_line($warpFile, 1, $key, $OPTABORT );
    $SURVEY   =  $tmp[0] ;
    $IDSURVEY = &get_IDSURVEY($SURVEY,-1);
    if ( $IDSURVEY < 0 ) { return ; }

    # Now we have a WARP-input file with a survey that is used.

    $NFILT_SURVEY[$IDSURVEY]       = 0 ;
    $SURVEY_FILTERLIST[$IDSURVEY]  = "" ;
    $SURVEY_FILTERFILES[$IDSURVEY] = "" ;
    $SURVEY_FILTERZP[$IDSURVEY]    = "" ;

    # get filter path
    $key    = "FILTPATH:" ;
    @tmp    =  sntools::parse_line($warpFile, 1, $key, $OPTABORT );
    $itmp   = index($tmp[0],"filters") + 8 ;
    $SURVEY_FILTERPATH[$IDSURVEY] = substr($tmp[0],$itmp,60);

    # get the filters.
    $key      = "FILTDEF:";
    @tmp      =  sntools::parse_line($warpFile, 3, $key, $OPTABORT );
    foreach $line ( @tmp ) {

	$NFILT_SURVEY[$IDSURVEY]++ ;
	$NFILT = $NFILT_SURVEY[$IDSURVEY] ;

	@wdlist  = split(/\s+/,$line ) ;
	$f     = $wdlist[0];
	$ffile = $wdlist[1];
	$zp    = $wdlist[2];

	$last = "$SURVEY_FILTERLIST[$IDSURVEY]";
	$SURVEY_FILTERLIST[$IDSURVEY] = "$last$f";

	$last = "$SURVEY_FILTERFILES[$IDSURVEY]";
	$SURVEY_FILTERFILES[$IDSURVEY] = "$last$ffile ";

	$last = "$SURVEY_FILTERZP[$IDSURVEY]";
	$SURVEY_FILTERZP[$IDSURVEY] = "$last$zp " ;
    }

    # get primary
    $key      = "PRIMARY:";
    @tmp      =  sntools::parse_line($warpFile, 1, $key, $OPTABORT );
    $itmp     = index($tmp[0],"standards") + 10 ;
    $SURVEY_PRIMARY[$IDSURVEY] = substr($tmp[0],$itmp,60);

    print " Found $SURVEY filters: $SURVEY_FILTERLIST[$IDSURVEY] \n";
    $| = 1;  # flush stdout

} # end of parse_WARPfile



# -----------------------------------
sub make_VERIFY_kcorFile {

    # construct kcor-input file and run kcor.exe for this survey.
    my ($IDSURVEY) = @_ ;

    my ($SURVEY, $PRIMARY);
    my ($inFile_kcor, $INFILE_KCOR, $file_kcor, $logFile_kcor, $prefix );
    my ($MAGSYS, $FILTPATH, $FILTLIST, @FILTFILES, @FILTZP );
    my ($NFILT, $ifilt, $cfilt, $file, $zp );

    $SURVEY      = $SURVEY_NAMES[$IDSURVEY];
    $PRIMARY     = $SURVEY_PRIMARY[$IDSURVEY];
    $FILTPATH    = $SURVEY_FILTERPATH[$IDSURVEY];
    $FILTLIST    = $SURVEY_FILTERLIST[$IDSURVEY];
    $NFILT       = $NFILT_SURVEY[$IDSURVEY];
    @FILTFILES   = split(/\s+/,$SURVEY_FILTERFILES[$IDSURVEY]) ;
    @FILTZP      = split(/\s+/,$SURVEY_FILTERZP[$IDSURVEY]) ;


    $prefix       = "kcor_${SURVEY}" ;
    $inFile_kcor  = "${prefix}.input" ;
    $INFILE_KCOR  = "$VERIFY_DIR/$inFile_kcor";
    $file_kcor    = "${prefix}.fits" ;
    $logFile_kcor = "${prefix}.log" ;

    $SURVEY_KCORFILE[$IDSURVEY] = "$file_kcor";

    $MAGSYS = "${SURVEY}-REF" ;

    print " Create $VERIFY_DIR/$file_kcor  \n";

    open  PTR_KCOR , "> $INFILE_KCOR" ;

    print PTR_KCOR "READ_ZPOFF:  0 \n";
    print PTR_KCOR "SN_SED:      $SED_SNIa \n"; 
    print PTR_KCOR "PRIMARY_SED: $MAGSYS  $PRIMARY \n"; 

    print PTR_KCOR "\n";
    print PTR_KCOR "MAGSYSTEM:  $MAGSYS \n"; 
    print PTR_KCOR "FILTSYSTEM: COUNT  \n";
    print PTR_KCOR "FILTPATH:   $FILTPATH \n";

    for( $ifilt=0; $ifilt < $NFILT; $ifilt++ ) {
	$cfilt = substr($FILTLIST,$ifilt,1);
	$file  = $FILTFILES[$ifilt] ;
	$zp    = $FILTZP[$ifilt] ;
	print PTR_KCOR "FILTER: $SURVEY-$cfilt  $file  $zp \n";
    }

    print PTR_KCOR "\n";
    print PTR_KCOR "LAMBDA_RANGE:  2100  12000 \n";
    print PTR_KCOR "OUTFILE:       $file_kcor \n";

    close PTR_KCOR ;

    # finally run the kcor job
    qx(cd $VERIFY_DIR; kcor.exe $inFile_kcor > $logFile_kcor);


} # end of make_VERIFY_kcorFile


# -----------------------------------
sub make_VERIFY_SIMLIB {

    # construct SIMLIB to simulate this survey.
    my ($IDSURVEY) = @_ ;

    my ($SIMLIB_FILE, $SIMLIB_FILE);
    my ($SURVEY, $FILTLIST, $NFILT );

    $SURVEY      = $SURVEY_NAMES[$IDSURVEY];
    $FILTLIST    = $SURVEY_FILTERLIST[$IDSURVEY];
    $NFILT       = $NFILT_SURVEY[$IDSURVEY];

    $SIMLIB_FILE = "VERIFY_${SURVEY}.SIMLIB" ;
    $SURVEY_SIMLIB[$IDSURVEY] = "$SIMLIB_FILE";

    print " Create $VERIFY_DIR/$SIMLIB_FILE \n" ;
    open  PTR_SIMLIB , "> $VERIFY_DIR/$SIMLIB_FILE" ;
    
    print PTR_SIMLIB "SURVEY:  $SURVEY      TELESCOPE: $SURVEY\n";
    print PTR_SIMLIB "FILTERS: $FILTLIST \n";

    print PTR_SIMLIB "USER: $USER    HOST: $HOST\n";
    print PTR_SIMLIB "COMMENT: $CADENCE_SIMLIB_VERIFY day cadence to verify non1a \n";
    print PTR_SIMLIB "BEGIN \n\n";

    my ($LIBID, $NOBS, $NMJD, $imjd, $RADEC, $MJD, $MJDOFF, $ifilt, $IDEXPT);
    my ($cMJD, $cPSF, $cZPT, $f);

    my $GAIN   = 1.0 ;
    my $RON    = 0.1 ; # readout noise
    my $SKYSIG = 100.0 ;
    my $PSFSIG = 1.5 ;

    $LIBID = 1;
    $NMJD  = int(120/$CADENCE_SIMLIB_VERIFY);
    $NOBS  = $NFILT * $NMJD ;
    $RADEC = "RA: 0.001  DECL: 0.001" ;
    $MJDOFF  = $SIM_PEAKMJD - 40.0 ;

    $cPSF = "$PSFSIG 0 0" ;
    $cZPT = "35 0.001   -99."  ;

    print PTR_SIMLIB "LIBID: $LIBID \n";
    print PTR_SIMLIB "$RADEC   NOBS: $NOBS    MWEBV: 0.00   PIXSIZE: 0.5 \n\n";

    for  ( $imjd = 1; $imjd <= $NMJD; $imjd++ ) {

	$MJDOFF += $CADENCE_SIMLIB_VERIFY ;
	$IDEXPT  = sprintf("%3.3d", $imjd);

	for ($ifilt=0; $ifilt < $NFILT; $ifilt++ ) {
	    $MJD  = $MJDOFF + $ifilt/100. ;
	    $cMJD = sprintf("%9.3f", $MJD);
	    $f    = substr($FILTLIST, $ifilt, 1);
	    print PTR_SIMLIB 
		"S: $cMJD  $IDEXPT  $f   $GAIN $RON  $SKYSIG  $cPSF  $cZPT\n";
	}	
    }
    print PTR_SIMLIB "END_LIBID: $LIBID\n\n";
    print PTR_SIMLIB "END_OF_SIMLIB: \n";
 

    close PTR_SIMLIB ;

} # end of make_VERIFY_SIMLIB


# ----------------------------------
sub make_VERIFY_simInputFile {

    my ($ISED) = @_ ;

    # create sim-input file with redshift, RA, DEC of original 
    # nonIa SN used to make SED template.


    my ($inputFile, $logFile, $SURVEY, $FILTLIST, $NFILT);
    my ($SIMLIB_FILE, $KCOR_FILE, $ZZ, $PP, $I3SED, $STRING);
    my ($IDSURVEY, $SNID, $ZSN, $CIDOFF, $I3NDX, $GENVER, $GENPRE );

    # strip off SED info
    $IDSURVEY    = $SEDVER_IDSURVEY[$ISED] ;
    $SNID        = $SEDVER_SNID[$ISED];
    $ZSN         = $SEDVER_REDSHIFT[$ISED];
    $I3NDX       = $SEDVER_INDEX[$ISED];

    # strip off SURVEY info
    $SURVEY      = $SURVEY_NAMES[$IDSURVEY];
    $FILTLIST    = $SURVEY_FILTERLIST[$IDSURVEY];
    $NFILT       = $NFILT_SURVEY[$IDSURVEY];
    $KCOR_FILE   = $SURVEY_KCORFILE[$IDSURVEY] ;
    $SIMLIB_FILE = $SURVEY_SIMLIB[$IDSURVEY] ;


    # construct strings needed for input file.
#    $CIDOFF = $I3NDX - 1;
    $CIDOFF = $SIM_CID - 1;


    $STRING = "${I3NDX}_${SURVEY}-${SNID}";

    $GENPRE    = "VERIFY${STRING}";
    $GENVER    = "${GENV_PREFIX}_${GENPRE}";
    $SEDVER_GENVERSION[$ISED] = $GENVER ;
    $SEDVER_GENPREFIX[$ISED]  = $GENPRE ;

    $ZZ = "$ZSN $ZSN";
    $PP = "$SIM_PEAKMJD $SIM_PEAKMJD";


    # start making sim-input file
    $inputFile = "SIMGEN${STRING}.INPUT";
    $logFile   = "SIMGEN${STRING}.LOG";
    $SEDVER_SIMFILE_INPUT[$ISED] = $inputFile ;
    $SEDVER_SIMFILE_LOG[$ISED]   = $logFile ;

    print "${SURVEY}-${SNID} ";
#    print " Create/run $VERIFY_DIR/$inputFile \n" ;

    open  PTR_INFILE , "> $VERIFY_DIR/$inputFile" ;

    print PTR_INFILE "NGENTOT_LC:   1 \n";
    print PTR_INFILE "SIMLIB_FILE:  $SIMLIB_FILE \n";
    print PTR_INFILE "GENVERSION:   $GENVER \n" ;
    print PTR_INFILE "GENPREFIX:    $GENPRE \n" ;
    print PTR_INFILE "GENSOURCE:    RANDOM \n";
    print PTR_INFILE "GENFILTERS:   $FILTLIST \n";
    print PTR_INFILE "KCOR_FILE:    $KCOR_FILE \n";

    print PTR_INFILE "\n";
    print PTR_INFILE "GENMODEL: NON1A \n";
    print PTR_INFILE "NON1A_KEYS: 5 \n" ;
    print PTR_INFILE "         INDEX   WGT   MAGOFF  MAGSMEAR  SNTYPE \n";
    print PTR_INFILE "NON1A:    $I3NDX    1.0    0.0      0.0       2   \n";

    print PTR_INFILE "\n" ;
    print PTR_INFILE "GENRANGE_REDSHIFT:  $ZZ \n";
    print PTR_INFILE "GENRANGE_PEAKMJD:   $PP \n";
    print PTR_INFILE "GENRANGE_TREST:     -18  60 \n";

    print PTR_INFILE "\n";
    print PTR_INFILE "RANSEED:         12345 \n";
    print PTR_INFILE "CLEARPROMPT:     0     \n";
    print PTR_INFILE "EXTINC_MILKYWAY: 0     \n";
    print PTR_INFILE "GENSIGMA_MWEBV_RATIO: 0.0 \n";
    print PTR_INFILE "GENTAU_AV:       0     \n";
    print PTR_INFILE "GENRANGE_AV:     0  0  \n";
    print PTR_INFILE "CIDOFF: $CIDOFF \n";
    print PTR_INFILE "FORMAT_MASK:     2   # text/terse \n";

    print PTR_INFILE "\n";

    print PTR_INFILE "OMEGA_MATTER:  0.3  \n"  ;
    print PTR_INFILE "OMEGA_LAMBDA:  0.7  \n"  ;
    print PTR_INFILE "W0_LAMBDA:    -1.00 \n"  ;
    print PTR_INFILE "H0:           70.0  \n"  ;

    close PTR_INFILE ;

} #  end of &make_VERIFY_simInputFile


# -----------------------
sub run_VERIFY_sim {

    # run the simulation for this ISED.
    # Then replace the integer CID with the string name.

    my ($ISED ) = @_ ;

    # local var
    my ($inputFile, $logFile, $cmd, $cdVER ) ;

    $inputFile = $SEDVER_SIMFILE_INPUT[$ISED] ;
    $logFile   = $SEDVER_SIMFILE_LOG[$ISED] ;
    $cdVER     = "cd $VERIFY_DIR";   
    $cmd       = "snlc_sim.exe $inputFile > $logFile";

    system("$cdVER ; $set_SNANA_MODELPATH ; $cmd");

    # now replace integer $SIM_CID with string name (SNID)

    my ($SNID, $GENV, $listFile, $simFile, $tmpFile, $sedcmd, $SIMFILE );
    $SNID     = $SEDVER_SNID[$ISED] ;
    $GENV     = $SEDVER_GENVERSION[$ISED];
    $listFile = "$SNDATA_ROOT/SIM/$GENV/$GENV.LIST";
    $simFile  = `cat $listFile` ;
    $simFile  =~ s/\s+$// ;  # remove trailing blanks
    $SIMFILE  = "$SNDATA_ROOT/SIM/$GENV/$simFile" ;

    $tmpFile = "TEMP_$simFile" ;
    $sedcmd  = "sed -e 's/$SIM_CID/$SNID/g'";
    qx($sedcmd $SIMFILE > $tmpFile ; mv $tmpFile $SIMFILE );


} # end of run_VERIFY_sim


sub merge_VERIFY_sims {

    my ($IDSURVEY) = @_ ;

    # merge simulated files for this survey

    my ($SIMDIR, $auxFile, $GENV, $SURVEY, $tmpList, $VV);

    $SURVEY = $SURVEY_NAMES[$IDSURVEY] ;
    $SURVEY_VERSION[$IDSURVEY] = "${GENV_PREFIX}_VERIFY_${SURVEY}" ;

    $GENV   = $SURVEY_VERSION[$IDSURVEY] ;
    $SIMDIR = "$SNDATA_ROOT/SIM/$GENV";

    print " Merge $SURVEY sims -> $GENV \n";

    # create merged-sim directory
    qx(mkdir $SIMDIR);

    # copy data files into merged area
    $VV      = "${GENV_PREFIX}_VERIFY" ;
    $tmpList = "$SNDATA_ROOT/SIM/${VV}*${SURVEY}*/VER*.DAT";
    qx(cp $tmpList $SIMDIR/);

    # create aux files
    $auxFile = "$SIMDIR/${GENV}.README";
    qx(touch $auxFile);

    $auxFile = "$SIMDIR/${GENV}.IGNORE";
    qx(touch $auxFile);

    $auxFile = "$SIMDIR/${GENV}.LIST";
    qx(cd $SIMDIR ; ls VER*.DAT >  $auxFile);

    
} # end of merge_VERIFY_sims


# ------------------------------
sub run_VERIFY_snana {

    my ($IDSURVEY) = @_ ;

    # run snana job on version containing SN from this IDSURVEY

    my ($nmlFile, $hisFile, $logFile, $GENV, $SURVEY, $prefix );

    $SURVEY = $SURVEY_NAMES[$IDSURVEY] ;
    $GENV   = $SURVEY_VERSION[$IDSURVEY] ;

    print "\t Run snana on $GENV to make light curve plots. \n";
    $| = 1;  # flush stdout

    # start by create input namelist file for snana.exe
    $prefix  = "SIMANA_VERIFY_${SURVEY}" ;
    $nmlFile = "${prefix}.nml";
    $hisFile = "${prefix}.his";
    $logFile = "${prefix}.log";

    @SURVEY_HISFILE[$IDSURVEY] = $hisFile ;

    open  PTR_NMLFILE , "> $VERIFY_DIR/$nmlFile" ;

    print PTR_NMLFILE "# Read $GENV and fill light curve plots.\n";
    print PTR_NMLFILE "# No fitting. \n";
#    print PTR_NMLFILE "# Do NOT abort when survey changes. \n";

    print PTR_NMLFILE "\n";
    print PTR_NMLFILE " &SNLCINP \n";
    print PTR_NMLFILE "   VERSION_PHOTOMETRY = '$GENV' \n";
    print PTR_NMLFILE "   HFILE_OUT          = '$hisFile' \n";
    print PTR_NMLFILE "   OPT_LCPLOT         = 1 \n";
#    print PTR_NMLFILE "   ABORT_ON_BADSURVEY = F \n";
    print PTR_NMLFILE "   NFIT_ITERATION     = 0 \n";
    print PTR_NMLFILE " &END \n";

    close PTR_NMLFILE ;

    # -----------------------------
    # run the SNANA job to get light curve plots
    qx(cd $VERIFY_DIR ; snana.exe $nmlFile > $logFile );

} # end of run_VERIFY_snana


# --------------------------------
sub make_VERIFY_lcplots {

    # Make and execute kumac file that plots light curves
    # for each SURVEY. Each light curve plot includes 
    # simulated points (every 2 days), original data,
    # and smoothed curve used to make the templates.

    my ($kFile, $logFile, $VDIR);
    my ($ID, $hisFile, $args, $plotPrefix, $psFile, $pdfFile );

    print "\n Make light curve plots ... \n";

    $plotPrefix = "NON1A_VERIFY" ;
    $kFile      = "$plotPrefix.kumac" ;
    $logFile    = "$plotPrefix.log" ;
    $psFile     = "$plotPrefix.ps" ;
    $pdfFile    = "$plotPrefix.pdf" ;
    $VDIR       = $VERIFY_DIR ;

    open  PTR_KUMAC , "> $VDIR/$kFile" ;

    if ( $USE_PRIVATE_KUMACS == 0 ) 
    { sntools::write_loginMacroLines(\*PTR_KUMAC,0); }

    $| = 1;  # flush stdout

    print PTR_KUMAC "exec journal \n";
    print PTR_KUMAC "fortran/file 66 $psFile \n";
    print PTR_KUMAC "metafile 66 -111 \n";
    print PTR_KUMAC "\n";
    $args = "ovnon1a=1 tmin=-25 tmax=60 prompt=0" ;

    for ( $ID=1; $ID <= $NSURVEY; $ID++ ) {
	$hisFile = $SURVEY_HISFILE[$ID] ;
	print PTR_KUMAC "h/file 1 $hisFile \n";
	print PTR_KUMAC "exec snana#fitres $args \n" ;
	print PTR_KUMAC "close 1 \n" ;
	print PTR_KUMAC "\n";
    }
    print PTR_KUMAC "close 66 \n";

    qx(cd $VDIR ; paw -b $kFile > $logFile);
    qx(cd $VDIR ; ps2pdf $psFile );

    print " See $VDIR/$pdfFile \n";

    close PTR_KUMAC ;

} # end of make_VERIFY_lcplots


# =====================================================
#
# Peak B band Mags in Rest frame
#
# =====================================================

sub PEAKMAG {

    # Jan 18, 2012
    # driver to make plot of peak B band mag in rest-frame.
    # Use simulation to generate B band at very low redshift,
    # and using the sim-input file created by this script.
    # Then plot MB_obs - MU
    # PeakMags are plotted for both RAW (MAGOFF=MAGSMEAR=0)
    # and for nominal to see how the raw distributions are
    # modified.

    my ($simFile, $incFile, $dirList, $OPT );

    print "\n" ;
    print "# ============================================= \n";
    print " Compute Peak B-band mags (mean and RMS) \n";
    print " in $PEAKMAG_DIR \n";
    $| = 1;  # flush stdout

    # make local PEAKMAG_DIR  for kcor/simlib/sim-input files
    if ( -d $PEAKMAG_DIR ) { qx(rm -r $PEAKMAG_DIR); }
    qx(mkdir $PEAKMAG_DIR);

    # remove old VERIFY versions to avoid confusion
    $dirList = "$SNDATA_ROOT/SIM/${GENV_PREFIX}_PEAKMAG*" ;
    qx(rm -r $dirList 2>/dev/null );
    


    &parse2_SEDLIST();

    &make_PEAKMAG_kcorFile();

    &make_PEAKMAG_SIMLIB();
    
    # OPT=0 for RAW (MAGOFF+MAGSMEAR=0)
    # OPT=1 or nominal

    for ($OPT=0; $OPT <=1; $OPT++ ) {
	print "   ----------------------------------------- \n";
	print "   $PEAKMAG_COMMENT[$OPT] \n";
	print "\t NON1A include File = $simFile_include[$OPT] \n";
	my $lnk = "ln -s $SIMFILE_INCLUDE[$OPT] $simFile_include[$OPT]";
	qx(cd $PEAKMAG_DIR ; $lnk );
	
	&make_PEAKMAG_simInputFile($OPT);

	&run_PEAKMAG_sim($OPT);
    }

    &make_PEAKMAG_plots();

} # end of PEAKMAG


# ------------------------------
sub make_PEAKMAG_kcorFile {


    my ($inFile_kcor, $INFILE_KCOR, $logFile_kcor );
    my ($ifilt, $cfilt, $file );

    my $MAGSYS   = "VEGA";
    my $FILTLIST = "$FILTLIST_PEAKMAG" ;    
    my $FILTPATH = "$FILTPATH_PEAKMAG" ;    

    print "\t Create dummy Kcor   file with  $FILTLIST ($NFILT_PEAKMAG)\n";
    $| = 1;  # flush stdout

    $logFile_kcor = "KCOR_PEAKMAG.log" ;
    $inFile_kcor  = "KCOR_PEAKMAG.input";
    $INFILE_KCOR  = "$PEAKMAG_DIR/$inFile_kcor";


    open  PTR_KCOR , "> $INFILE_KCOR" ;

    print PTR_KCOR "READ_ZPOFF:  0 \n";
    print PTR_KCOR "SN_SED:      $SED_SNIa \n"; 
    print PTR_KCOR "VEGA_SED:  alpha_lyr_stis_003.dat \n";

    print PTR_KCOR "\n";
    print PTR_KCOR "MAGSYSTEM:  $MAGSYS \n"; 
    print PTR_KCOR "FILTSYSTEM: COUNT  \n";
    print PTR_KCOR "FILTPATH:   $FILTPATH \n";

    for( $ifilt=0; $ifilt < $NFILT_PEAKMAG; $ifilt++ ) {
	$cfilt = substr($FILTLIST,$ifilt,1);
	$file  = "${FILTPATH}_${cfilt}.dat" ;
	print PTR_KCOR "FILTER:  $FILTPATH-$cfilt  $file  0.02 \n";
    }

    print PTR_KCOR "\n";
    print PTR_KCOR "LAMBDA_RANGE:  2100  12000 \n";
    print PTR_KCOR "OUTFILE:       $KCOR_FILE_PEAKMAG \n";

    close PTR_KCOR ;

    # finally run the kcor job
    qx(cd $PEAKMAG_DIR; kcor.exe $inFile_kcor > $logFile_kcor);


} # end of make_PEAKMAG_kcorFile


sub make_PEAKMAG_SIMLIB {

    print "\t Create dummy SIMLIB with  $FILTLIST_PEAKMAG \n";
    $| = 1;  # flush stdout

    open  PTR_SIMLIB , "> $PEAKMAG_DIR/$SIMLIB_FILE_PEAKMAG" ;
    
    print PTR_SIMLIB "SURVEY: LOWZ      TELESCOPE:  UNKNOWN \n";
    print PTR_SIMLIB "FILTERS: $FILTLIST_PEAKMAG \n";

    print PTR_SIMLIB "USER: $USER    HOST: $HOST\n";
    print PTR_SIMLIB "COMMENT: To get peakmags in rest-frame \n";
    print PTR_SIMLIB "BEGIN \n\n";

    my ($LIBID, $NOBS, $NMJD, $imjd, $RADEC, $MJD, $MJDOFF, $ifilt, $IDEXPT);
    my ($cMJD, $cPSF, $cZPT, $f);

    my $GAIN    = 1.0 ;
    my $RON     = 0.1 ; # readout noise
    my $SKYSIG  = 100.0 ;
    my $PSFSIG  = 1.5 ;
    my $CADENCE = $CADENCE_SIMLIB_PEAKMAG ;

    $LIBID = 1;
    $NMJD  = int(100./$CADENCE);
    $NOBS  = $NFILT_PEAKMAG * $NMJD ;
    $RADEC = "RA: 0.001  DECL: 0.001" ;
    $MJDOFF  = $SIM_PEAKMJD - 2.*$CADENCE ;

    $cPSF = "$PSFSIG 0 0" ;
    $cZPT = "35 0.001   -99."  ;

    print PTR_SIMLIB "LIBID: $LIBID \n";
    print PTR_SIMLIB "$RADEC   NOBS: $NOBS    MWEBV: 0.00   PIXSIZE: 0.5 \n\n";

    for  ( $imjd = 1; $imjd <= $NMJD; $imjd++ ) {

	$MJDOFF += $CADENCE ;
	$IDEXPT  = sprintf("%3.3d", $imjd);

	for ($ifilt=0; $ifilt < $NFILT_PEAKMAG; $ifilt++ ) {
	    $MJD  = $MJDOFF + $ifilt/100. ;
	    $cMJD = sprintf("%9.3f", $MJD);
	    $f    = substr($FILTLIST_PEAKMAG, $ifilt, 1);
	    print PTR_SIMLIB 
		"S: $cMJD  $IDEXPT  $f   $GAIN $RON  $SKYSIG  $cPSF  $cZPT\n";
	}	
    }
    print PTR_SIMLIB "END_LIBID: $LIBID\n\n";
    print PTR_SIMLIB "END_OF_SIMLIB: \n";

    close PTR_SIMLIB ;


} # end of make_PEAKMAG_SIMLIB 


sub make_PEAKMAG_simInputFile {

    my ($OPT) = @_ ;

    my ($GENVER, $GENPRE, $ZZ, $PP, $CIDOFF, $finp, $flog );

    print "\t Create dummy SIM-INPUT file with  $FILTLIST_PEAKMAG \n";
    $| = 1;  # flush stdout
    $ZZ = "0.0002 0.0002" ;
    $PP = "$SIM_PEAKMJD $SIM_PEAKMJD" ;

    $GENVER = $GENVERSION_PEAKMAG[$OPT] ;
    $GENPRE = $GENPREFIX_PEAKMAG ;
    $CIDOFF = 0 ;

    $finp = $SIMINP_FILE_PEAKMAG[$OPT] ;

    open  PTR_INFILE , "> $PEAKMAG_DIR/$finp" ;

    print PTR_INFILE "NGENTOT_LC:   $NGEN_PEAKMAG[$OPT] \n";
    print PTR_INFILE "SIMLIB_FILE:  $SIMLIB_FILE_PEAKMAG \n";
    print PTR_INFILE "GENVERSION:   $GENVER \n" ;
    print PTR_INFILE "GENPREFIX:    $GENPRE \n" ;
    print PTR_INFILE "GENSOURCE:    RANDOM \n";
    print PTR_INFILE "GENFILTERS:   $FILTLIST_PEAKMAG \n";
    print PTR_INFILE "KCOR_FILE:    $KCOR_FILE_PEAKMAG \n";

    print PTR_INFILE "\n";
    print PTR_INFILE "GENMODEL: NON1A \n";
    print PTR_INFILE "INPUT_FILE_INCLUDE: $simFile_include[$OPT] \n" ;
    print PTR_INFILE "\n" ;
    print PTR_INFILE "GENRANGE_REDSHIFT:  $ZZ \n";
    print PTR_INFILE "GENRANGE_PEAKMJD:   $PP \n";
    print PTR_INFILE "GENRANGE_TREST:     -18  40 \n";

    print PTR_INFILE "\n";
    print PTR_INFILE "RANSEED:         12345 \n";
    print PTR_INFILE "CLEARPROMPT:     0     \n";
    print PTR_INFILE "EXTINC_MILKYWAY: 0     \n";
    print PTR_INFILE "GENSIGMA_MWEBV_RATIO: 0.0 \n";
    print PTR_INFILE "GENTAU_AV:       0     \n";
    print PTR_INFILE "GENRANGE_AV:     0  0  \n";
    print PTR_INFILE "CIDOFF: $CIDOFF \n";
    print PTR_INFILE "FORMAT_MASK:     2   # text/terse \n";

    print PTR_INFILE "\n";

    print PTR_INFILE "OMEGA_MATTER:  0.3  \n"  ;
    print PTR_INFILE "OMEGA_LAMBDA:  0.7  \n"  ;
    print PTR_INFILE "W0_LAMBDA:    -1.00 \n"  ;
    print PTR_INFILE "H0:           70.0  \n"  ;

    # dump key
    my ($ifilt, $NVAR, $VARLIST, $f); 
    $NVAR    = 5 + $NFILT_PEAKMAG ;
    $VARLIST = "CID Z MU GENTYPE NON1A_INDEX ";
    for ($ifilt=0; $ifilt < $NFILT_PEAKMAG; $ifilt++ ) {
	$f    = substr($FILTLIST_PEAKMAG, $ifilt, 1);
	$VARLIST = "$VARLIST MAGT0_$f  ";
    }

    print PTR_INFILE "\n";    
    print PTR_INFILE "SIMGEN_DUMP: $NVAR  $VARLIST \n" ;

    close PTR_INFILE ;

} # end of make_PEAKMAG_simInputFile


sub run_PEAKMAG_sim {
    
    my ($OPT) = @_ ;

    my ($cdd, $cmd, $chmod, $dmpFile, $GENV, $mv   );

    my $finp = $SIMINP_FILE_PEAKMAG[$OPT] ;
    my $flog = $SIMLOG_FILE_PEAKMAG[$OPT] ;
    my $NGEN = $NGEN_PEAKMAG[$OPT] ;

    print "\t Generate $NGEN random nonIa using simulation ... \n";
    $GENV  = "$GENVERSION_PEAKMAG[$OPT]" ;
    $cdd   = "cd $PEAKMAG_DIR" ;

    $cmd   = "snlc_sim.exe $finp > $flog" ;
    $chmod = "chmod o-w $SNDATA_ROOT/SIM/$GENV" ;
    qx($cdd ; $set_SNANA_MODELPATH ; $cmd ; $chmod );
    
    # create ntuple from DUMP file
    $dmpFile = "$SNDATA_ROOT/SIM/${GENV}/${GENV}.DUMP";
    $cmd     = "combine_fitres.exe $dmpFile";
    $mv      = "mv combine_fitres.tup  $NTUP_FILE_PEAKMAG[$OPT]";
    qx($cdd ; $cmd; $mv );

} # end of run_PEAKMAG_sim



sub make_PEAKMAG_plots {

    # Plot rest-frame m_X = MAGT0_X - MU (where X = B,V,R)
    # Make plot separately for each GENTYPE,
    # and also for nominal (black) and raw (red).
    # RAW => MAGOFF = MAGSMEAR = 0.
    #
    # Feb 14, 2012: do B,V,R instead of just B
    # Nov 14, 2012: fix bug related to PRIVATE_KUMACS if-block

    my ($kFile, $logFile, $PKDIR, $i, $ITYPE, $CTYPE );
    my ($hisFile, $args, $plotPrefix, $psFile, $pdfFile );

    my $idraw    = 20 ;    # raw (MAGOFF=MAGSMEAR=0)
    my $idnom    = 21 ;    # nominal
    my $hid_raw  = $idraw ;
    my $hid_nom  = $idnom ;
    my $hid_tmp  = 30 ;
    my $magBins  = "56 -22 -8";
    my $q        = "\'";
    my $ntplot   = "7788.magt0_[FILT]-mu" ;

    print "   --------------------------------------- \n";
    print "   Make PEAKMAG plot for each SN Type. \n";

    $plotPrefix = "NON1A_PEAKMAG" ;
    $kFile      = "$plotPrefix.kumac" ;
    $logFile    = "$plotPrefix.log" ;
    $psFile     = "$plotPrefix.ps" ;
    $pdfFile    = "$plotPrefix.pdf" ;
    $PKDIR      = "$PEAKMAG_DIR" ;

    open  PTR_KUMAC , "> $PKDIR/$kFile" ;
    
    if ( $USE_PRIVATE_KUMACS == 0 ) { 
	sntools::write_loginMacroLines(\*PTR_KUMAC,0) ; 
	print PTR_KUMAC "exec NON1A_PEAKMAG#all \n" ;
	print PTR_KUMAC "return \n";
	print PTR_KUMAC "\n";
      }
    
    print PTR_KUMAC "macro all \n";
    print PTR_KUMAC "  exec journal \n";
    print PTR_KUMAC "  fortran/file 66 $psFile \n";
    print PTR_KUMAC "  metafile 66 -111 \n";
    print PTR_KUMAC "\n";
    
    print PTR_KUMAC "  h/file $idraw $NTUP_FILE_PEAKMAG[0] \n";
    print PTR_KUMAC "  h/file $idnom $NTUP_FILE_PEAKMAG[1] \n";


    print PTR_KUMAC "  opt stat ; set stat 001110 \n" ;
    print PTR_KUMAC "  set csiz .45 \n";
    print PTR_KUMAC "  set hwid 5 ; set lwid 5 \n" ;
    print PTR_KUMAC "  set yhti 1.1 ; set tsiz .4 \n" ;
    print PTR_KUMAC "\n" ;
    print PTR_KUMAC "  exec ${plotPrefix}#plotFilt B  \n";
    print PTR_KUMAC "  exec ${plotPrefix}#plotFilt V  \n";
    print PTR_KUMAC "  exec ${plotPrefix}#plotFilt R  \n";
    print PTR_KUMAC "\n";

    print PTR_KUMAC "  close 66 \n" ;
    print PTR_KUMAC "return \n";

    print PTR_KUMAC "\n\n";

    print PTR_KUMAC "* ============================  \n";
    print PTR_KUMAC "macro plotFilt \n" ;
    print PTR_KUMAC "  FILT = [1] \n" ;	
    print PTR_KUMAC "  zone 2 3 \n\n" ;

    for ( $i=1; $i <= $NTYPE_USER; $i++ ) {
	$CTYPE = $CTYPE_USER[$i]  ;
	$ITYPE = $ITYPE_USER[$i]  ;

	print PTR_KUMAC "* ------------------------------------ \n";
	print PTR_KUMAC "  cut 1 gentype=$ITYPE \n";
	print PTR_KUMAC "  chis = 'peak m?' // [FILT] // '! for $CTYPE ' \n";
	print PTR_KUMAC "  1dhis $hid_raw [chis] $magBins \n";
	print PTR_KUMAC "  1dhis $hid_nom [chis] $magBins \n";
	print PTR_KUMAC "  nt/proj $hid_raw //lun${idraw}/${ntplot} \$1  \n";
	print PTR_KUMAC "  nt/proj $hid_nom //lun${idnom}/${ntplot} \$1  \n";

	# plot nomina peak-mB distribution
	print PTR_KUMAC "  hpl $hid_nom \n";

	# get max bin contents into Bmax, and scale hid_raw to Bmax/2
	print PTR_KUMAC "  Bmax  = \$HINFO($hid_nom,'MAX') \n";
	print PTR_KUMAC "  Bmax2 = [Bmax]/2.001 \n" ;

	# overlay RAW distribution in red
	print PTR_KUMAC "  if ( [Bmax2] > 0 ) then \n";
	print PTR_KUMAC "    set hcol 1002 \n";
	print PTR_KUMAC "    hcopy [Bmax2]*$hid_raw/$hid_raw $hid_tmp \n";
	print PTR_KUMAC "    hpl $hid_tmp hist,s \n";
	print PTR_KUMAC "    set hcol 1 \n";
	print PTR_KUMAC "  endif \n";

	# re-plot nominal so that it's not covered by the red(raw)
	print PTR_KUMAC "  hpl $hid_nom s \n";

	print PTR_KUMAC "  h/del $hid_raw,$hid_nom,$hid_tmp ; cut 1 - \n" ;

	if ( $i == 1 ) {
	    print PTR_KUMAC "\n";
	    print PTR_KUMAC "* ========= LEGEND ============== \n";
	    print PTR_KUMAC "  null 0 1 0 1 -as \n";
	    print PTR_KUMAC "  key 0.40 1.25  1 'NOMINAL' 0.8 L \n"; 
	    print PTR_KUMAC "  set txci 2; set pmci 2 \n"; 
	    print PTR_KUMAC "  key 0.40 1.10 21 'RAW (MAGOFF=0, MAGSMEAR=0)' \n"; 
	    print PTR_KUMAC "  set txci 1; set pmci 1 \n"; 
	    print PTR_KUMAC "* ======================= \n";
	}
    }

    print PTR_KUMAC "return \n" ;

    close PTR_KUMAC ;

    qx(cd $PKDIR ; paw -b $kFile > $logFile) ;
    qx(cd $PKDIR ; ps2pdf $psFile ) ;

    print " See $PKDIR/$pdfFile \n";
    $| = 1;  # flush stdout


} # end of make_PEAKMAG_plots 


sub clean {
    if ( $DO_CLEAN == 0 ) { return ; }  
    my $cmdClean = "rm -r $SNDATA_ROOT/SIM/${GENV_PREFIX}*" ;
    print "\n Remove \$SNDATA_ROOT/SIM/${GENV_PREFIX}* \n";
    qx($cmdClean);

}  # end of clean

