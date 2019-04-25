#!/usr/bin/perl
#
# Created Jun 2011 by R.Kessler
#
# Generate simulated efficiency map as an arbitrary 
# function of SN model parameters.
# [replaces older SIMGEN_MAPGEN.cmd that was hard-wired for MLCS]
#
# Usage:
#   SIMEFF_MAPGEN.pl <inFile>
#   SIMEFF_MAPGEN.pl <inFile> --BINS <minbin> <maxbin>
#   SIMEFF_MAPGEN.pl <inFile> KILL
#
# Note that the --BINS option is called internally
# from this script.
#
# Aug 02, 2011: use strict
# Aug 29, 2011: call check_binsize to make sure that 
#               (NBIN-1)*Binsize = full range.
#
# Sep 16, 2011: set $BIN=0 if $VMAX=$VMIN 
#               (avoid divide-by-zero bug caught by Scolnic)
#
# Sep 19, 2011: allow interactive mode if there are no NODELIST keys
#
# Sep 23, 2011: allow the eff table to depend on sim-input keys that 
#               contain "OMEGA" or "LAMBDA". Such tables are NOT for
#               fitting, but for studying dependence on cosmo params.
#
# Mar 23, 2012: 
#     replace GEN_SNDATA_SIM 0  with  FORMAT_MASK 0
#     If sim-abort is detected, save log file to ABORT_XXX.LOG
#     and print usual error-message to screen.
#
# ------------------

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# -----------------------------
# declare subs

sub init_MAPGEN ;
sub parse_inFile ;
sub RUN_SIMEFF_MAPGEN ;
sub open_effmapFile;
sub assign_nodes ;
sub ssh_nodes ;
sub interactive_mapgen ;
sub check_binsize ;

# ----------------------------
# declare globals

my $q   = '"' ;

my ($MODE, $KILLJOBS_FLAG, @NODELIST, $NNODE, @MSGERR );
my ($MINBIN, $MAXBIN, $NBINTOT, @NODE_MINBIN, @NODE_MAXBIN );
my (@GENVAR_SCALE, @GENVAR_SIMNAME, @GENVAR_MAPNAME);
my (@GENVAR_NBIN, @GENVAR_MIN, @GENVAR_MAX, @GENVAR_BIN, @GENVAR_NPROD);
my ($INPUT_FILE, $SIMGEN_INFILE, $SIMGEN_VERSION_PREFIX, $NGENVAR );
my ($MAPGEN_OUTFILE,$SIMGEN_LOGFILE, $ABORT_LOGFILE, $SIMGEN_VERSION );
my ($EFFMAP_OUTFILE, $SIMGEN_EFFERR, $SNANA_LOGIN_SETUP);

# ----------- BEGIN MAIN ------------
&init_MAPGEN();

# parse command-line arguments
&parse_arg();

# parse input file
&parse_inFile();

if ( $KILLJOBS_FLAG ) { sntools::Killjobs(@NODELIST) ; die "\n"; }

if ( $MODE == 2 ) { &RUN_SIMEFF_MAPGEN(); exit ; }

# open output eff-map file
&open_effmapFile();

# assign bin-range for each node
&assign_nodes();

# do the ssh for each node
if ( $NNODE > 0 ) 
{ &ssh_nodes(); }
else
{ &interactive_mapgen(); }


# ==================================
#
#        END OF MAIN
#
# ===================================

# ========================
sub init_MAPGEN {

    $MODE   = 1 ;  # 1 = user-interactive;  2=ssh 
    $KILLJOBS_FLAG = 0 ;
    $MINBIN = -9 ;
    $MAXBIN = -9 ;
    $NBINTOT = 1;
} 

# ========================
sub parse_arg {

    # parse command-line arguments

    my ($NARG, $arg, $i) ;

    $NARG = scalar(@ARGV);
    if ( $NARG < 1 ) {
        die " Must give input filename as argument => ABORT \n";
    }

    $INPUT_FILE = $ARGV[0] ;    

    for ( $i = 1; $i < $NARG ; $i++ ) {
        $arg = $ARGV[$i];

        if ( $arg eq "KILL"   ) { $KILLJOBS_FLAG = 1 ; }

        if ( $arg eq "--BINS"   ) { 
	    $MODE      = 2 ;
	    $MINBIN    = $ARGV[$i+1] ;
	    $MAXBIN    = $ARGV[$i+2] ;
	}

    }
} # end of parse_arg



# ========================
sub parse_inFile {

    my ($key, $ivar, $tmpLine, @tmp, @words, $NBIN, @TMPNODES) ;
    my ($VMIN, $VMAX, $BIN, $SCALE, $V2V, $cvar, $cnbin );

    my $OPT_ABORT = 0 ; # parsing option
    my $OPT_WARN  = 2 ; # 1=> leave warning message; 2=>quiet
    my $inFile    = $INPUT_FILE ;

    # ------------

    @NODELIST ;
    $key = "NODELIST:" ;
    @TMPNODES    = sntools::parse_line($INPUT_FILE, 99, $key, $OPT_WARN ) ;
    foreach $tmpLine ( @TMPNODES ) {
        @tmp = split(/\s+/,$tmpLine) ;
        @NODELIST =  ( @NODELIST , @tmp );
    }
    $NNODE    = scalar(@NODELIST);


    $SNANA_LOGIN_SETUP = "" ;
    $key = "SNANA_LOGIN_SETUP:" ;
    @tmp = sntools::parse_line($inFile, 99, $key, $OPT_WARN) ;
    $SNANA_LOGIN_SETUP = "@tmp" ;

    # ------
    $key   = "SIMGEN_INFILE:" ;
    @tmp   = sntools::parse_line($inFile, 1, $key, $OPT_ABORT ) ; 
    $SIMGEN_INFILE = "$tmp[0]" ;

    # grep the GENVERSION  from the  simgen-input file
    $key   = "GENVERSION:" ;
    @tmp   = sntools::parse_line($SIMGEN_INFILE, 1, $key, $OPT_ABORT ) ; 
    $SIMGEN_VERSION_PREFIX = "$tmp[0]" ;
    
    $key   = "EFFMAP_OUTFILE:" ;
    @tmp   = sntools::parse_line($inFile, 1, $key, $OPT_ABORT ) ; 
    $EFFMAP_OUTFILE = "$tmp[0]" ;

    $key   = "SIMGEN_EFFERR:" ;
    @tmp   = sntools::parse_line($inFile, 1, $key, $OPT_ABORT ) ; 
    $SIMGEN_EFFERR = "$tmp[0]" ;

    $key   = "GENVAR:" ;
    @tmp   = sntools::parse_line($inFile, 6, $key, $OPT_ABORT ) ; 
    $NGENVAR = scalar(@tmp);
    for ( $ivar=0; $ivar < $NGENVAR; $ivar++  ) {
	$tmpLine = $tmp[$ivar];
	@words   = split(/\s+/,$tmpLine) ;

	$NBIN    = $words[3] ;

	@GENVAR_SCALE    = ( @GENVAR_SCALE    , $words[0] );
	@GENVAR_SIMNAME  = ( @GENVAR_SIMNAME  , $words[1] );
	@GENVAR_MAPNAME  = ( @GENVAR_MAPNAME  , $words[2] );
	@GENVAR_NBIN     = ( @GENVAR_NBIN     , $words[3] );
	@GENVAR_MIN      = ( @GENVAR_MIN      , $words[4] );
	@GENVAR_MAX      = ( @GENVAR_MAX      , $words[5] );


	$NBINTOT = $NBINTOT * $NBIN ; # total number of eff-bins

	# get bin size
	$VMIN  = $GENVAR_MIN[$ivar] ;
	$VMAX  = $GENVAR_MAX[$ivar] ;

	if ( $VMAX > $VMIN ) 
	{ $BIN   = ($VMAX-$VMIN)/($NBIN - 1.0) ; }
	else
	{ $BIN   = 0.0 ; }

	@GENVAR_BIN   = ( @GENVAR_BIN  , $BIN );

	# make sure that bin size x integer = total range
	&check_binsize($ivar);
    }


    # compute GENVAR_PROD
    $NBINTOT = 1;
    for ( $ivar = $NGENVAR-1; $ivar >= 0; $ivar--  ) {
	@GENVAR_NPROD    = ( @GENVAR_NPROD    , $NBINTOT );
	$NBINTOT *= $GENVAR_NBIN[$ivar];
    }

    # summarize inputs    
    if ( $MODE == 1 ) {
	print " SIMGEN_INFILE:   $SIMGEN_INFILE \n";
	print " SIMGEN_EFFERR:   $SIMGEN_EFFERR \n";
	print " SIMGEN_VERSION:  ${SIMGEN_VERSION_PREFIX}_XXX \n";
	print " NODELIST:        @NODELIST \n";
	print " GENVAR_NPROD     @GENVAR_NPROD \n";

	for ( $ivar=0; $ivar < $NGENVAR; $ivar++  ) {
	    $SCALE = $GENVAR_SCALE[$ivar] ;
	    $NBIN  = $GENVAR_NBIN[$ivar] ;
	    $VMIN  = $GENVAR_MIN[$ivar] ;
	    $VMAX  = $GENVAR_MAX[$ivar] ;

	    $V2V   = "$VMIN to $VMAX" ;
	    $cvar  = sprintf("%-6s", $GENVAR_MAPNAME[$ivar] );
	    $cnbin = sprintf("%3d", $NBIN );
	    print "  Generate $SCALE $cvar with $cnbin  bins from  $V2V\n";
	}
	print "\n";

    } # MODE=1
    


} # end of parse_inFile


# ======================
sub  check_binsize {
    my ($ivar) = @_ ;

    # check that when BIN size is trunctated to 3 decimal places,
    # NBIN * BIN = full range. The truncation comes in writing
    # out the values to the text (SIMEFF) file.

    my ($NBIN, $NN, $VMIN, $VMAX, $BIN, $BIN3, $NAME );
    my ($RANGE_EXACT, $RANGE_TMP, $RDIF, $cNxB );

    $NAME  = $GENVAR_MAPNAME[$ivar] ;
    $NBIN  = $GENVAR_NBIN[$ivar] ;
    $VMIN  = $GENVAR_MIN[$ivar] ;
    $VMAX  = $GENVAR_MAX[$ivar] ;
    $BIN   = $GENVAR_BIN[$ivar] ;
    $BIN3  = int($BIN*1000.0 + 0.5)/1000. ;

    $RANGE_EXACT = $VMAX - $VMIN ;
    $RANGE_TMP   = ($NBIN-1) * $BIN3 ;
    $RDIF        = abs($RANGE_EXACT - $RANGE_TMP) ;

#    print " xxx $NAME : BIN=$BIN  RANGE-dif = $RDIF \n";

    if ( $RDIF > 1.0E-10 ) {

	$cNxB      = "(NBIN-1)*BIN" ;
	$MSGERR[0] = "Input $NAME range is '$RANGE_EXACT'" ;
	$MSGERR[1] = "but $cNxB = ($NBIN-1) * $BIN3 = '$RANGE_TMP'" ;
	$MSGERR[2] = "Adjust $NAME binning so that $cNxB = $RANGE_EXACT" ;
	sntools::FATAL_ERROR(@MSGERR);
    }

} # end of check_binsize

# =========================
sub open_effmapFile {

    my ($ivar, $cdate, $who) ;
    my ($SCALE, $VMIN, $VMAX, $V2V, $NBIN, $NAME, $ctmp, $cvar, $cnbin );

    # open effmap file and write header info

    # remove current file if it exists
    if ( -e $EFFMAP_OUTFILE ) { qx(rm $EFFMAP_OUTFILE) ; }

    $cdate  = `date +%F` ;
    $who    = `whoami` ;

    $cdate         =~ s/\s+$// ;   # trim trailing whitespace
    $who           =~ s/\s+$// ;   # trim trailing whitespace

    print " Open output file for eff-map: $EFFMAP_OUTFILE \n";
    open PTR_EFFMAP , "> $EFFMAP_OUTFILE" ;

    print PTR_EFFMAP "Efficiency Map created by $who on $cdate \n";
    print PTR_EFFMAP "Command: $0 \n";
    print PTR_EFFMAP "\n";

    # setup table header

    print PTR_EFFMAP "NGENVAR: $NGENVAR \n";
    for ( $ivar=0; $ivar < $NGENVAR; $ivar++  ) {
	$SCALE = $GENVAR_SCALE[$ivar] ;
	$NBIN  = $GENVAR_NBIN[$ivar] ;
	$VMIN  = $GENVAR_MIN[$ivar] ;
	$VMAX  = $GENVAR_MAX[$ivar] ;

	$V2V   = sprintf("%8.3f  %8.3f", $VMIN, $VMAX) ;
	$cvar  = sprintf("%-8s", $GENVAR_MAPNAME[$ivar] );
	$cnbin = sprintf("%3d", $NBIN );
	print PTR_EFFMAP "GENVAR:  $SCALE   $cvar $cnbin  $V2V \n";
    }

    # create table header for efficiencies
    print PTR_EFFMAP "\n\n";
    print PTR_EFFMAP "#       ";

    for ( $ivar=0; $ivar < $NGENVAR; $ivar++  ) {
	$SCALE = $GENVAR_SCALE[$ivar] ;
	$NAME  = $GENVAR_MAPNAME[$ivar] ;

	if ( $SCALE eq "INV" ) 
	{  $ctmp  = "1/$NAME" ; }
	elsif ( $SCALE eq "LOG" ) 
	{  $ctmp  = "LOG($NAME)" ; }
	else
	{  $ctmp  = $NAME ; }

	$cvar = sprintf("%-8s", $ctmp);
	print PTR_EFFMAP "$cvar " ; 

    }

    $cvar = sprintf("%-6s", "EFF" );
    print PTR_EFFMAP "$cvar ";

    $cvar = sprintf("%-6s", "ERR(EFF)" );
    print PTR_EFFMAP "$cvar ";

    print PTR_EFFMAP "\n";

    print PTR_EFFMAP 
    "# ------------------------------------------------------------------ \n";

    close PTR_EFFMAP ;

} # end of open_effmapFile


# =================================
sub assign_nodes {

    # assigne range of GENVAR bins to each node.

    my ($inode, $NTMP, $imin, $imax, $cnode, $cmin, $cmax) ;

    if ( $NNODE == 0 ) { return ; }

    $NTMP = int(.5+($NBINTOT)/$NNODE) ;

    print "\n Split $NBINTOT GENVAR bins among $NNODE nodes ($NTMP per node). \n";

    for ( $inode=0; $inode < $NNODE; $inode++  ) {
	$imin =  $inode * $NTMP + 1;
	$imax =  $imin  + $NTMP - 1;
	if ( $inode == $NNODE - 1 ) {
	    $imax = $NBINTOT ;
	}
	@NODE_MINBIN = ( @NODE_MINBIN , $imin ) ;
	@NODE_MAXBIN = ( @NODE_MAXBIN , $imax ) ;

	$cnode = sprintf("%-10s", $NODELIST[$inode] ) ;
	$cmin  = sprintf("%6d", $imin );
	$cmax  = sprintf("%6d", $imax );
	print "    Absolute 1D bin-range($cnode) : $cmin to $cmax \n"; 
    }


} # end of assign_nodes


# =================================
sub ssh_nodes {

    my ($inode, $node, $imin, $imax, $workDir, $cdd, $cmd, $args) ;

    $workDir = `pwd` ;
    $cdd = "cd $workDir" ;

    print "\n Launch jobs ... \n";
    for ( $inode=0; $inode < $NNODE; $inode++  ) {

	$node =	$NODELIST[$inode] ;
	$imin = $NODE_MINBIN[$inode] ;
	$imax = $NODE_MAXBIN[$inode] ;

	$args = "$INPUT_FILE  --BINS $imin $imax";
	$cmd = "ssh -x $node ${q} $cdd ; $SNANA_LOGIN_SETUP ; $0 $args ${q}" ;
	print " $cmd \n";
	system("$cmd &");
	sleep(1);
    }

} # end of ssh_nodes


# ==============================
sub interactive_mapgen {

    my ( $imin, $imax, $args, $cmd );

    $imin = 1;
    $imax = $NBINTOT ;
    $args = "$INPUT_FILE  --BINS $imin $imax";
    $cmd = "$0 $args" ;
    print " $cmd \n";
    system("$cmd &");

} # end of interactive_mapgen

# =============================
sub RUN_SIMEFF_MAPGEN {

    # Driver routine to simulate the effic. map for MINBIN to MAXBIN.

    my ($IBIN1D, $prefix, $cmin, $cmax);

    # construct temp eff-storage file using bin-range and  GENVERSION
    $cmin = sprintf("%6.6d", $MINBIN );
    $cmax = sprintf("%6.6d", $MAXBIN );

    $prefix = "TEMP_${SIMGEN_VERSION_PREFIX}_${cmin}-${cmax}";
    $MAPGEN_OUTFILE = "${prefix}.DAT";
    $SIMGEN_LOGFILE = "${prefix}.LOG";
    $SIMGEN_VERSION = "${SIMGEN_VERSION_PREFIX}_${cmin}" ;

    # define ABORT log-file if it's needed
    $prefix = "ABORT_${SIMGEN_VERSION_PREFIX}_${cmin}-${cmax}";
    $ABORT_LOGFILE  = "${prefix}.LOG";

    if ( -e $MAPGEN_OUTFILE ) { qx(rm $MAPGEN_OUTFILE); }
    if ( -e $SIMGEN_LOGFILE ) { qx(rm $SIMGEN_LOGFILE); }
    if ( -e $ABORT_LOGFILE  ) { qx(rm $ABORT_LOGFILE);  }

    open PTR_MAPGEN , "> $MAPGEN_OUTFILE" ;

    for ( $IBIN1D=$MINBIN; $IBIN1D <= $MAXBIN ; $IBIN1D++ ) {
	&run_simeff($IBIN1D);
    }

    close PTR_MAPGEN ;

    &COMBINE_SIMEFF_MAPGEN();

} # end of RUN_SIMEFF_MAPGEN


# =======================
sub run_simeff {

    # translate 'IBIN1D' into bin-index for each GENVER,
    # and the run the simulation to get the efficiency.
    # Finally, update PTR_MAPGEN with the GENVAR values
    # and the simulated efficiency.
    #
    # IBIN1D-1 = i0 + N0*i1 + N0*N1*i2 + N0*N1*N2*i3 ...
    # and need to invert IBIN1D to get i0, i1 ..
    # Here IBIN1D=1 ... , and i1, i2 ... start at 0
    # 

    my ($IBIN1D) = @_;

    my ($ibin, $ivar, $LASTOFF, $value, $SIM_VALUE, $VV, $SCALE ) ;
    my ($simArgs, $cmd, $NPROD, $VARNAME );
    my (@ibinList , @VALUE_LIST, @cvalueList) ;
    my ($NARG, $VARG, $NOWRITE, $ERRARG) ;
    my ($jdum1, $jdum2, @words, @tmp) ;

    $LASTOFF = 0;

    # specify GENVERSION and turn off data-file output
    $VARG    = "GENVERSION $SIMGEN_VERSION";
#    $NOWRITE = "GEN_SNDATA_SIM 0  SIMGEN_DUMP -1" ;
    $NOWRITE = "FORMAT_MASK 0  SIMGEN_DUMP -1" ;
    $ERRARG  = "EFFERR_STOPGEN $SIMGEN_EFFERR EFFRANGE_STOPGEN 0 1";
    $NARG    = "NGEN_LC 100000  NGENTOT_LC 0" ;
    $simArgs = "$VARG $NOWRITE $ERRARG $NARG";
    
    # get ibin and value for each variable.
    for ( $ivar=0; $ivar < $NGENVAR; $ivar++  ) {

	# last ivar changes 1st; 1st ivar changes last.
	$NPROD   = $GENVAR_NPROD[$NGENVAR-1-$ivar]  ;
	$VARNAME = $GENVAR_SIMNAME[$ivar] ;

	$ibin    = int(( $IBIN1D - 1 - $LASTOFF) / $NPROD ) ;

	$value   = $GENVAR_MIN[$ivar] + $ibin * $GENVAR_BIN[$ivar] ;
	$LASTOFF += ( $ibin * $NPROD ) ;

	# $value might be the log or inverse, so get the actual SIM_VALUE.
	$SCALE = $GENVAR_SCALE[$ivar] ;
	if ( $SCALE eq "INV" ) 
	{  $SIM_VALUE  = 1./$value ; }
	elsif ( $SCALE eq "LOG" ) 
	{  $SIM_VALUE  = 10.0**$value ; }
	else
	{  $SIM_VALUE  = $value ; }

	$ibinList[$ivar]    = $ibin ;
	$VALUE_LIST[$ivar]  = $SIM_VALUE ;
	$cvalueList[$ivar]  = sprintf("%7.3f ", $value) ;

	
	# if VARNAME containes "OMEGA" or "LAMBDA", then this
	# is a fixed cosmological parameter and only one value
	# is given (since it's not a RANGE).

	$jdum1 = index($VARNAME,"OMEGA") ;
	$jdum2 = index($VARNAME,"LAMBDA") ;

	if ( $jdum1 < 0 && $jdum2 < 0 ) {
	    $VV = "$SIM_VALUE $SIM_VALUE" ;
	}
	else {
	    $VV = "$SIM_VALUE" ;
	}
	$simArgs = "$simArgs " . "$VARNAME $VV " ;
    }


    # run simulations
    $cmd = "snlc_sim.exe $SIMGEN_INFILE $simArgs > $SIMGEN_LOGFILE";
    qx($cmd);  

    # if we have a  sim-abort then mv TEMP logFile to ABORT logfile
    # and abort this script.   
    @tmp    = qx(grep " ABORT " $SIMGEN_LOGFILE );
    if ( scalar(@tmp) > 0 ) {
	qx(mv $SIMGEN_LOGFILE $ABORT_LOGFILE);
	$MSGERR[0] = "The following sim-job aborted:";
	$MSGERR[1] = "$cmd";
	$MSGERR[2] = "";
	$MSGERR[3] = "For details see '$ABORT_LOGFILE'";
	sntools::FATAL_ERROR(@MSGERR);
    }

    # parse log file to get efficiency
    my ($cEFF, $cERR) ;;
    @tmp    = qx(grep "SEARCH+CUTS" $SIMGEN_LOGFILE );
    @words  = split(/\s+/,$tmp[0]) ;
    $cEFF   = sprintf("%6.4f", $words[3]);
    $cERR   = sprintf("%6.4f", $words[5]);

    # update efficiency map
    print PTR_MAPGEN "EFF:  @cvalueList   $cEFF  $cERR \n" ;

} # end of run_simeff


# ==============================
sub COMBINE_SIMEFF_MAPGEN {
    
    # if/when all jobs have finsished, used 'cat' to combine
    # the all of the eff-map files into the one EFFMAP_OUTFILE
    # Also clean up junk.

    my ($NEFF, $tmpLine, @words) ;

    # all jobs are done when the number of EFF lines = NBINTOT

    # grep the number of 'EFF:' lines 
    $tmpLine = qx(grep "EFF:" TEMP*.DAT | wc );
    @words   = split(/\s+/,$tmpLine) ;
    $NEFF    = $words[1] ;

    # $HOST = `hostname` ;
    # print "\t NEFF = $NEFF on host = $HOST \n";

    if ( $NEFF == $NBINTOT ) {
	print "\n\n" ;
	print " All SIMEFF_MAPGEN jobs are done => \n";
	print " combine everything into $EFFMAP_OUTFILE \n";
	qx(cat TEMP*.DAT >> $EFFMAP_OUTFILE);
	# qx(rm TEMP*.DAT); 
	print "\n Done. Check EFF-MAP file : $EFFMAP_OUTFILE \n";

	&SIMEFF_CLEANUP();
    }


} # end of COMBINE_SIMEFF_MAPGEN


# ==========================
sub  SIMEFF_CLEANUP {
    # remove TEMP files
    qx(rm TEMP_${SIMGEN_VERSION_PREFIX}*) ;
} # end of SIMEFF_CLEANUP {
