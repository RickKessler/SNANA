#!/usr/bin/perl
#
# Created Mar 2013 by R.Kessler
# Run simulation (snlc_sim.exe) to generate SN-photometry templates 
# for a multi-dimensional grid of parameters:
#   redshift
#   shape   (stretch, x1, Delta, dm15 ...)
#   color     (AV or SALT2c)
#   colorLaw  (RV or Beta)
#   Trest   
# This grid is needed by the psnid.exe program, and can also be used
# by any program that needs a compact model description without
# model-dependent codes.
#
#
# Syntax:
#   sim_SNgrid.pl  <list of sim-input files>
#   sim_SNgrid.pl  <list of sim-input files>  PARALLEL
#   sim_SNgrid.pl  <list of sim-input files>  TEXT
#
# Comments on inputs:
# - list of sim-input files can include wildcards. For example,
#   the following are the same
#      sim_SNgrid.pl  sim_DES_SALT2.input  sim_DES_NON1A.input
#      sim_SNgrid.pl  sim_DES*
#   Files with a tilde (sim_DES_SALT2.input~) are ignored.
#
# - PARALLEL option launches all jobs in parallel on the SAME node.
#   Make sure there are enough cores and memory for the number
#   of jobs. 
#
# - TEXT is a debug-option to generate grids in TEXT format.
#   Note that psnid.exe CANNOT read TEXT format.
#
# Comments on output:
# - all output is generated in $OUTDIR = SNgrid/  .
#   $OUTDIR is created if not there, but otherwise this script writes
#   into the existing directory. However, if you re-run the same
#   sim_SNgrid job, then existing files in $OUTDIR are clobbered.
#   The idea is that you can keep adding more grids to the same
#   $OUTDIR.
#
# - Default grid file-name is  "GRID_$SURVEY_$GENMODEL.FITS" or
#   .FITS -> .TEXT for the TEXT option. SURVEY is taken from the SIMLIB 
#   and GENMODEL is the base of the GENMODEL key in the sim-input file. 
#   Thus for SALT2.Guy10, only the 'SALT2' part is used for the file name. 
#   You can over-write the default file-name prefix with the following 
#   sim-input key,
#       PREFIX_GRIDFILE:  BLABLA
#   This key is used only by sim_SNgrid.pl; the key is ignored by 
#   snlc_sim.exe.
#
# Sep 15 2017: allow ENV
# May 05 2022: fix a few bugs after a few years of rust.
#
# --------------------------------------

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;
use File::Basename;


my $OPT_ABORT = 0 ;  # flag to abort if key is missing in file
my $OPT_QUIET = 2 ;  # flag to staty quiet of key is missing
my $SNDATA_ROOT   = $ENV{'SNDATA_ROOT'};
my $HOST          = $ENV{'HOST'};
my $OUTDIR        = "SNgrid";
my $SIMJOB        = "snlc_sim.exe" ;
my $SCRIPTNAME    = "$0" ; 

my (@INFILE_LIST, @SURVEY_LIST, @GENMODEL_LIST, @PREFIX_LIST);
my (@GENVERSION_LIST, @GRIDFILE_LIST, @GRIDFILE_BASE_LIST, @INCFILE_LIST );
my ($FORMAT, $OPT_PARALLEL, @MSGERR);

my ($NJOB, $ijob, $inFile);

sub initStuff ;
sub parse_args ;
sub prep_simInput_Files ;
sub prep_outdir ;
sub genGrid ;

# ============ BEGIN MAIN =============


&initStuff();

&parse_args();

$NJOB = scalar(@INFILE_LIST);

&prep_outdir();

foreach $inFile ( @INFILE_LIST ) { &prep_simInput_Files($inFile); }

my $Twait = 3 ;
print "\n Will launch jobs in $Twait seconds ... \n";
sleep($Twait);

print "\n";
print "# =========================================== \n";
print "#       BEGIN PROCESSING GRIDs  \n" ;
print "# =========================================== \n";
print "\n";


for ( $ijob=0; $ijob < $NJOB; $ijob++ ) { genGrid($ijob) ; }


print "\n" ;
if ( $OPT_PARALLEL == 0 ) {
    print " Done. \n" ; 
}
else {
    print " Launched $NJOB parallel jobs on $HOST \n";
    print " Wait for all jobs to finished. \n" ;

  CHECK:
    my ($N, @bla) ;    
    sleep(10);
    @bla = qx(ps | grep $SIMJOB );
    $N = scalar(@bla) ;
    if ( $N > 0 ) { 
	print "\t $N  $SIMJOB remaining ... please wait ...\n";
	goto CHECK ; 
    }
}


my (@bla,$NOK);
# xxx mark delete @bla = qx(grep Wrote  $OUTDIR/*.LOG);
@bla = qx(grep NGENLC_WRITE  $OUTDIR/*.LOG); # May 5, 2022
$NOK = scalar(@bla);

print " Output is in  ${OUTDIR}/ \n";
print " $NOK of $NJOB jobs finished OK. \n";

if ( $NOK < $NJOB ) {
    print " WARNING: at least one sim job failed !!! Check $OUTDIR/*.LOG \n";
}

# ========== END MAIN ============

sub initStuff {
    $FORMAT        = "FITS" ;
    $OPT_PARALLEL  = 0 ;   # serial processing by default
    @INFILE_LIST   = () ;
    @SURVEY_LIST   = () ;
    @GENMODEL_LIST = () ;
    @PREFIX_LIST   = () ;
}

# ========================
sub parse_args {

    my ($NARG, $i, $arg);

    $NARG = scalar(@ARGV);
    if ( $NARG < 1 ) {
	$MSGERR[0] = "Must give at least one sim-input file as arg,";
	sntools::FATAL_ERROR(@MSGERR);
    }

 
    for ( $i = 0; $i < $NARG ; $i++ ) {
        $arg = $ARGV[$i];
	
	if ( $arg  eq "TEXT" ) {
	    $FORMAT = "TEXT" ;
	}
	elsif ( $arg  eq "PARALLEL" ) {
	    $OPT_PARALLEL = 1 ;  # submit jobs in paralle on same node
	}
	else {
	    # add to list unless there is a tilde at the end
	    my $j = index($arg,'~');
	    if ( $j < 0 ) {  @INFILE_LIST = (@INFILE_LIST , "$arg" ); }
	}
    }

} # end of parse_args


# ===================================
sub prep_simInput_Files {

    my ($inFile) = @_ ;

    # prepare sim-input file:
    # - Extract SURVEY, GENMODEL and optional PREFIX
    # - make sure at least one file exists
    # 

    my (@tmp, $KEY, $SIMLIB_FILE, $KCOR_FILE, $KCOR_SYMLINK, $ORIG );
    my ($SURVEY, $GENMODEL, $genmodel, $PREFIX, $GENVERSION );
    my ($GRIDFILE, $GRIDFILE_BASE, $INCFILE ) ;

    print " -------------------------------------- \n";
    print " Prepare sim-input file: $inFile \n" ;

    if ( !(-e $inFile) ) {
	$MSGERR[0] = "sim-input file '$inFile' does not exist." ;
	sntools::FATAL_ERROR(@MSGERR);	
    }

    # ------------------------------------
    # first check for required 'NGRID' keys in the in-input file
    my $DOCHECK = 1;
    if ( $DOCHECK ) {
	@tmp=sntools::parse_line($inFile,0, "NGRID_LOGZ:",     $OPT_ABORT);
	@tmp=sntools::parse_line($inFile,0, "NGRID_COLORPAR:", $OPT_ABORT);
	@tmp=sntools::parse_line($inFile,0, "NGRID_COLORLAW:", $OPT_ABORT);
	@tmp=sntools::parse_line($inFile,0, "NGRID_TREST:",    $OPT_ABORT);
    }

    # ---------------------------------------------
    # start by getting the survey from the simlib;
    # First grep the simlib, then grep the SURVEY
    $KEY  = "SIMLIB_FILE:" ;
    @tmp  = sntools::parse_line($inFile, 1, $KEY, $OPT_ABORT ) ;
    $SIMLIB_FILE = "$tmp[0]" ;
    $SIMLIB_FILE = qx(echo $SIMLIB_FILE); # unpack ENV
    $SIMLIB_FILE =~ s/\s+$// ;   # trim trailing whitespace  

    if ( -e $SIMLIB_FILE ) {
	qx(cp $SIMLIB_FILE  $OUTDIR/);
    }
    else {
	$ORIG        = $SIMLIB_FILE ;
	$SIMLIB_FILE = "$SNDATA_ROOT/simlib/$SIMLIB_FILE" ;
    }    
    if ( ! (-e $SIMLIB_FILE ) ) {
	$MSGERR[0] = "Could not find SIMLIB_FILE '$ORIG' ";	
	$MSGERR[1] = "needed by $inFile" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    # ---------------
    # xxxxxx mark delete May 5 2022 xxxxxxxxxx
    #
    # if KCOR_FILE is local, then make symbolic link in OUTDIR
    #$KEY  = "KCOR_FILE:" ;
    #@tmp  = sntools::parse_line($inFile, 1, $KEY, $OPT_ABORT ) ;
    #$KCOR_FILE = "$tmp[0]" ;
    #$KCOR_FILE = qx(echo $KCOR_FILE); # unpack ENV
    #$KCOR_FILE =~ s/\s+$// ;   # trim trailing whitespace  

    #$KCOR_SYMLINK = "$OUTDIR/$KCOR_FILE" ;
    #if ( (-e $KCOR_FILE) &&  !(-e $KCOR_SYMLINK ) )
    #{ qx(cd $OUTDIR ; ln -s ../$KCOR_FILE $KCOR_FILE); } 
    #
    # xxxxxxxxx end mark xxxxxxxx


    # -------------------
    $KEY  = "SURVEY:" ;
    @tmp  = sntools::parse_line($SIMLIB_FILE, 1, $KEY, $OPT_ABORT ) ;
    $SURVEY = "$tmp[0]" ;
    
    # -----------------------------
    # Now get GENMODEL
    $KEY  = "GENMODEL:" ;
    @tmp  = sntools::parse_line($inFile, 1, $KEY, $OPT_ABORT ) ;
    $GENMODEL = "$tmp[0]" ;
    $genmodel = $GENMODEL ;

    # extract sub-string before the dot
    my $idot = index($GENMODEL, '.') ;
    if ( $idot > 0 ) { $genmodel = substr($GENMODEL, 0, $idot ); }

    # -------------------------------
    # get GENVERSION so we know where the sim-output is
    $KEY  = "GENVERSION:" ;
    @tmp  = sntools::parse_line($inFile, 1, $KEY, $OPT_ABORT ) ;
    $GENVERSION = "$tmp[0]" ;

    # -----------------------
    # check for optional prefix; if not there then construct it.
    
    $KEY  = "PREFIX_GRIDFILE:" ;
    @tmp  = sntools::parse_line($inFile, 1, $KEY, $OPT_QUIET ) ;
    if ( scalar(@tmp) > 0 )  {
	$PREFIX = "$tmp[0]" ;
    }
    else {
	my $genmodel_base = basename($genmodel) ;
	$PREFIX = "GRID_${SURVEY}_${genmodel_base}" ;
    }



    # ----------
    # check for optional INCLUDE file
    $KEY     = "INPUT_FILE_INCLUDE:" ;
    $INCFILE = "" ; 
    @tmp  = sntools::parse_line($inFile, 1, $KEY, $OPT_QUIET ) ;    
    if ( scalar(@tmp) > 0 ) {  $INCFILE = "$tmp[0]" ; }
    $INCFILE = qx(echo $INCFILE); # unpack ENV
    $INCFILE =~ s/\s+$// ;   # trim trailing whitespace  

    # --------
    $GRIDFILE  = "${PREFIX}.$FORMAT" ;
    $GRIDFILE_BASE = basename($GRIDFILE) ;
	
    print "\t SURVEY   = $SURVEY \n";
    print "\t GENMODEL = $GENMODEL \n";
    print "\t GRIDFILE = $GRIDFILE \n" ;
    if ( length($INCFILE) > 0 )  { print "\t INCFILE  = $INCFILE \n" ; } 

    # global storage
    
    @SURVEY_LIST     = (@SURVEY_LIST ,     "$SURVEY"     );
    @GENMODEL_LIST   = (@GENMODEL_LIST ,   "$genmodel"   );
    @GENVERSION_LIST = (@GENVERSION_LIST , "$GENVERSION" );
    @PREFIX_LIST     = (@PREFIX_LIST ,     "$PREFIX"     );
    @GRIDFILE_LIST   = (@GRIDFILE_LIST ,   "$GRIDFILE"   );
    @GRIDFILE_BASE_LIST  = (@GRIDFILE_BASE_LIST , "$GRIDFILE_BASE"   );
    @INCFILE_LIST    = (@INCFILE_LIST  ,   "$INCFILE"    );

} # end of  prep_simInput_Files 


# =============================================
sub prep_outdir {

    print "\n ------------------------------------- \n";
    if ( -d $OUTDIR ) {
	print " OUTDIR $OUTDIR already exists.\n";
    }
    else {
	print " Created OUTDIR $OUTDIR .\n";
	qx(mkdir $OUTDIR);
    }

} # end of prep_outdir


# =========================================
sub genGrid {
    my ($ijob) = @_ ;

    my ($PREFIX, $INFILE, $GRIDFILE, $LOGFILE, $VERS );
    my ($CMD_SIM, $SIMARGS, $CMD_COPY, $GENGRID_FILE );
    my ($logFile, $gridFile, $gengridFile, $incFile, $gridFile_base );
    my (@tmp, $KEY);


    $INFILE   = $INFILE_LIST[$ijob] ;
    print " Process $INFILE ... \n";

    $PREFIX   = $PREFIX_LIST[$ijob] ;
    $VERS     = $GENVERSION_LIST[$ijob] ;

    $gridFile      = "$GRIDFILE_LIST[$ijob]" ;
    $gridFile_base = "$GRIDFILE_BASE_LIST[$ijob]" ;
    $GRIDFILE = "$OUTDIR/$gridFile_base" ;

    $logFile  = "${INFILE}.LOG" ;
    $LOGFILE  = "$OUTDIR/${logFile}" ;

    $gengridFile  = "${VERS}/${VERS}.GRID" ;
    $GENGRID_FILE = "$SNDATA_ROOT/SIM/$gengridFile" ;

    $incFile  = $INCFILE_LIST[$ijob] ;

    # clobber existing files.
    if ( -e $LOGFILE  ) { qx(rm $LOGFILE);  }
    if ( -e $GRIDFILE ) { qx(rm $GRIDFILE); }

    # copy sim-input file to outdir
    qx(cp  $INFILE $OUTDIR/ );

    # if there is an include file, copy that too.
    if ( length($incFile) > 0 ) { 
	print "\t copy INCLUDE file: '${incFile}' \n";
	qx(cp $incFile $OUTDIR/) ; 
    }

    $SIMARGS  = "$INFILE  GRID_FORMAT $FORMAT" ;
    $CMD_SIM  = "$SIMJOB  $SIMARGS" ;
    $CMD_COPY = "cp  $GENGRID_FILE $gridFile_base" ;

    print "\t Sim-job:  $CMD_SIM  \n" ;
    print "\t copy \$SNDATA_ROOT/SIM/$gengridFile \n" ;
    print "\t  --> $gridFile_base \n";

    # open log file and put global info up top

    my $whichJob = `which $SIMJOB`;
    $whichJob =~ s/\s+$// ;  # remove trailing blanks
    
    open  PTR_LOG , " > $LOGFILE" ;
    print PTR_LOG "####################################### \n";
    print PTR_LOG "# Simulate GRID with \n";
    print PTR_LOG "#   script:  $SCRIPTNAME \n";
    print PTR_LOG "#   sim-job: $whichJob \n" ;
    print PTR_LOG "#                            \n";
    print PTR_LOG "####################################### \n\n";
    close PTR_LOG ;

    # Create RUNJOB script
    my $cjob = sprintf("%2.2d", $ijob);
    my $RUNSCRIPT = "RUN${cjob}_${PREFIX}";
    open  PTR_RUNJOB , " > $OUTDIR/$RUNSCRIPT" ;
    print PTR_RUNJOB  "$CMD_SIM >> $logFile \n";
    print PTR_RUNJOB  "$CMD_COPY \n" ;
    close PTR_RUNJOB ;
    qx(chmod +x $OUTDIR/$RUNSCRIPT);


    if ( $OPT_PARALLEL == 0 ) 
    { qx(cd $OUTDIR ; ./$RUNSCRIPT ); }
    else
    { system("cd $OUTDIR ; ./$RUNSCRIPT &") ; }


} # end of genGrid
