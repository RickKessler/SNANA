#!/usr/bin/perl
#
#
# Created May 2007 by R.Kessler
#  [May 29, 2011: converted from C-shell to perl]
#
# Generic after-burner script to create pdf file with plot(s)
# of lightcurve fits,
# Input file can be HBOOK (--h $inFile) or ROOT (--r $inFile).
# Default output is a single pdf file with 1 set of SN light curves 
# (i.e., all filters) displayed per page. The 'eps' option generates
# a separate [prefix]_[cid].eps file per SN and these eps files
# can be embedded into a latex document. Note that you get either
# the global pdf file or one eps file per SN; to get both you must
# run this script twice: with and without the 'eps' option.
#
#
# Usage:
#   mkfitplots.pl -h $inFile   (for hbook file)
#   mkfitplots.pl -h $inFile   PRIVATE  (use ~/kumacs for debugging)
#
#   mkfitplots.pl -r  $inFile   (for root  file, Aug 2013)
#   mkfitplots.pl --r $inFile   (one or two dashed allowed for all keys)
#
#
# Usage Options:
#
# Set obs-epoch range relative to peak
#   mkfitplots.pl -h $inFile -tmin $tmin -tmax $tmax
#
# Plot just one CID,
#   mkfitplots.pl -h $inFile -cid <cid> 
#
# Make separate eps file for each SN (to embed in latex)
# With 'eps' option there is NO default ps file with all SNe plots.
#   mkfitplots.pl -h $inFile eps
#   mkfitplots.pl -h $inFile eps2gif  (convert each eps -> gif)
#
#
# separate output file for each SN (Nov 2014); eps for hbook, pdf for root
#   mkfitplots.pl -h $inFile multigif
#   mkfitplots.pl -h $inFile multipdf
#   mkfitplots.pl -h $inFile multieps
#
# Suppress fit-model overlay with 
#   mkfitplots.pl -h $inFile OVFITOFF
#
# Indicate peakFlux with
#   mkfitplots.pl -h $inFile PEAKFLUX
#
# Adjust fluxscale of plots with
#   mkfitplots.pl -h $inFile -fluxscale $fluxscale
#
# Plot mag instead of flux
#   mkfitplots.pl -h $inFile MAG  
#       or
#   mkfitplots.pl -h $inFile -fluxscale -1
#
# Suppress chi2 vs. epoch plots
#   mkfitplots.pl -h $inFile NOCHI2
#
# Suppress GRID
#   mkfitplots.pl -h $inFile NOGRID
#
# Adjust vertical page size (0 < vscale <= 1) with 
#   mkfitplots.pl -h $inFile -vscale $vscale
#
# Switch from vertical to horizontal layout (1 row)
#   mkfitplots.pl -h $inFile H1
#   mkfitplots.pl -h $inFile V1  (use vertical layout; default)
#
# Specify number of filters per page (default = min(7,Nfilt))
#   mkfitplots.pl -h $inFile -nrowpage $nrowpage
#
# Specify SN TYPE
#   mkfitplots.pl -h $inFile -type $type
#
#
# Plot residuals instead of light curves
#   mkfitplots.pl -h $inFile RESIDS
#
#
#
#
#  HISTORY
#
# May 29, 2011: extracting marginalized pdfs are now an option 
#               (see MARGPDF flag) instead of default. 
#
# Jun 22, 2011: add new options
#             --vscale <vscale> to scale vertical size of page
#             vert 0            switch from  vertical to horizontal layout
#             NOGRID            turn off grid (useful for papers)
#
#
# Jun 23, 2011: abort on any unrecognized command-line argument
#
# Jul 20, 2011: new argument --nrowpage
#
# Aug 02, 2011: use strict
#
# Aug 29, 2011: abort in non-existant his file  and update to use 
#               sntools::FATAL_ERROR utility
#
# Nov 28, 2011: add RESIDS option to plot residual intead of light curves.
#               Execute ps2pdf and remove ps file.
#
#               Add optional argument '--type $type'
#
# Oct 8 2012:  new eps2gif option (for web pages)
#
# Oct 30, 2012:  call sntools::write_loginMacroLines(\*PTR_KUMAC,1);
#
# Aug 22, 2013: major overhaul to work with either hbook or root.
#               Note that this script has been split into 3:
#    mkfitplots.pl (main),  mkfitplots_hbook.pl, mkfitplots_root.c
#    BEWARE: not yet tested with root.
#
#
# Aug 26, 2013: 
#   * new "MAG" option now does the same thing as "-fluxscale -1"
#   * Allow one or two dashes for each key,
#      e.g., '-r $rootFile$ and  '--r $rootFile' are both accepted
#
# Aug 7, 2014: new '-cid <cid>' option to plot just one LC.
#
# Nov 3, 2014: allow -R or -r, allow -H or -h
#              For multigif or multipdf, suppress 'snana' option for root.
#
# ---------------------------------------------

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;

use strict ;

# ---------------------------------
# define subs
sub set_defaults ;
sub set_defaults_HBOOK ;
sub set_defaults_ROOT ;
sub parse_args ;
sub isArgKey ;
sub set_plotArgs_HBOOK ;
sub set_plotArgs_ROOT ;
sub eps2gif ;
sub move_plotFiles ;
sub rename_pdf ;

# ------------------------------
# define globals

my $SNANA_DIR = $ENV{'SNANA_DIR'};

my $IFILE_HBOOK = 1 ;
my $IFILE_ROOT  = 2 ;

my $SCRIPT_HBOOK = "mkfitplots_hbook.pl" ;
my $SCRIPT_ROOT  = "$SNANA_DIR/util/mkfitplots_root.c"   ;

my $RESIDS_FLAG = 0 ;
my $q  = "\'" ;
my $qq = '"' ;

# user values (USRVAL) and options (USROPT)
my ($USRVAL_tmin, $USRVAL_tmax, $USRVAL_fscale, $USRVAL_vscale);
my ($USRVAL_type, $USRVAL_nrow, $USRVAL_cid );

my ($USROPT_ovfit, $USROPT_pchi2, $USROPT_grid, $USROPT_pkf );
my ($USROPT_vert,  $USROPT_eps2gif, $USROPT_margpdf );
my ($USROPT_multieps, $USROPT_multipdf, $USROPT_multigif );
my ($USROPT_redline0);

my ($USR_INFILE, $USR_INFILE_TYPE) ;


my ($plotArgs, @MSGERR, $CMDPLOT, @TMPOUT, $MULTIFLAG );

my $PRIVATE_FLAG = "";

# ---------------- BEGIN MAIN ----------------

print "Command: \n  $0 @ARGV \n\n";

&set_defaults();

&parse_args();

# set arguments for snana#fitres
if ( $USR_INFILE_TYPE == $IFILE_HBOOK )  {  
    &set_plotArgs_HBOOK(); 
    $CMDPLOT = "$SCRIPT_HBOOK $USR_INFILE $q$plotArgs$q $PRIVATE_FLAG" ;
}

if ( $USR_INFILE_TYPE == $IFILE_ROOT  )  {  
    &set_plotArgs_ROOT() ; 
    my $cid     = $USRVAL_cid ;
    my $rootCmd = "(${qq}${USR_INFILE}${qq},$cid,${qq}${plotArgs}${qq})";    
    $CMDPLOT = "root -l -b -q ${SCRIPT_ROOT}${q}${rootCmd}${q} 2>/dev/null" ;
}  

# execute plot command
print "\n CMDPLOT = \n$CMDPLOT \n" ;

@TMPOUT = qx($CMDPLOT) ;

if ( $USR_INFILE_TYPE == $IFILE_HBOOK )  {  
    print "@TMPOUT" ;  # print stdout from CMDPLOT 
}


# --------------------------------

    # first check eps -> gif option (for hbook only)
if (  $USROPT_eps2gif ) { &eps2gif(); }
 

if ( length($MULTIFLAG) > 0 ) {  
    &move_plotFiles();   # move plot-files into separate subdir
}
else {
    &rename_pdf();
}

# --------------------------------

print "\n Done.\n";

# ===========================
#  END OF MAIN
# ===========================


# =========================
sub set_defaults {

    $USR_INFILE      = "NULL" ;
    $USR_INFILE_TYPE = 0 ;

    $USRVAL_cid      =  -999 ;  # 0--> all
    $USRVAL_tmin     =  -30 ; 
    $USRVAL_tmax     =   80 ;
    $USRVAL_fscale   =  1.0 ;
    $USRVAL_vscale   =  1.0 ;
    $USRVAL_nrow     =   7 ;
    $USRVAL_type     =  -9 ;  # flag to include all types

    $USROPT_ovfit    =   1 ;  # include best-fit model
    $USROPT_grid     =   1 ;  # dashed grid
    $USROPT_pkf      =   0 ;  # show peak flux on plots
    $USROPT_margpdf  =   0 ;  # show marginalized pdfs (not working)

    $USROPT_eps2gif  =   0 ;  # separate gif file for each SN
    $USROPT_multieps =   0 ;  # separate eps file for each SN
    $USROPT_multipdf =   0 ;  # separate pdf-file for each SN
    $USROPT_multigif =   0 ;  # separate gif-file for each SN
    $USROPT_redline0 =   0 ;  # red horizontal line through 0

    $MULTIFLAG = "" ;

} # end of set_defaults

# =========================
sub set_defaults_HBOOK {

    $USROPT_pchi2    =  1 ;
    $USROPT_vert     =  1 ;  # vertical stacking
}

sub set_defaults_ROOT {
    $USROPT_pchi2    =  0 ;
    $USROPT_vert     =  0 ; # horizontal layout
}

# ======================
sub parse_args {

    my ($NARG, $ARG, $ARG2, $ii, $i, $iLast, $nrd, $ird, $NERR) ;
    my (@USE_ARG, @STACK) ;

    $NARG = scalar(@ARGV); 
    
    for ( $i=0; $i < $NARG; $i++ ) { $USE_ARG[$i] = 0 ; }
    $iLast = -1;

    for ( $i=0; $i < $NARG; $i++ ) {

	$ii = $i + 1 ;
	$nrd = 0;

	$ARG  = $ARGV[$i] ;
	$ARG2 = $ARGV[$ii] ;


	if (&isArgKey("h",$ARG)  || &isArgKey("H",$ARG) )  { 
	    $USR_INFILE  = $ARGV[$ii] ; $nrd=2; 
	    $USR_INFILE_TYPE = $IFILE_HBOOK ;
	    &set_defaults_HBOOK();
	}

	if (&isArgKey("r",$ARG) || &isArgKey("R",$ARG) )  {  
	    $USR_INFILE  = $ARGV[$ii] ; $nrd=2; 
	    $USR_INFILE_TYPE = $IFILE_ROOT ;
	    &set_defaults_ROOT();
	}

	if (&isArgKey("cid",      $ARG))  { $USRVAL_cid     = $ARG2; $nrd=2;}
	if (&isArgKey("tmin",     $ARG))  { $USRVAL_tmin    = $ARG2; $nrd=2;}
	if (&isArgKey("tmax",     $ARG))  { $USRVAL_tmax    = $ARG2; $nrd=2;}
	if (&isArgKey("fluxscale",$ARG))  { $USRVAL_fscale  = $ARG2; $nrd=2;}
	if (&isArgKey("vscale",   $ARG))  { $USRVAL_vscale  = $ARG2; $nrd=2;}
	if (&isArgKey("nrowpage", $ARG))  { $USRVAL_nrow    = $ARG2; $nrd=2;}
	if (&isArgKey("type",     $ARG))  { $USRVAL_type    = $ARG2; $nrd=2;}

	if ( $ARGV[$i] eq "OVFITOFF" )   
	{ $USROPT_ovfit   = 0; $USROPT_pchi2=0 ; $nrd=1;}

	if ( $ARGV[$i] eq "NOCHI2"   )   { $USROPT_pchi2   = 0; $nrd=1;}
	if ( $ARGV[$i] eq "NOGRID"   )   { $USROPT_grid    = 0; $nrd=1;}
	if ( $ARGV[$i] eq "PEAKFLUX" )   { $USROPT_pkf     = 1; $nrd=1;}
	if ( $ARGV[$i] eq "MARGPDF"  )   { $USROPT_margpdf = 1; $nrd=1;}
	if ( $ARGV[$i] eq "H1"       )   { $USROPT_vert    = 0; $nrd=1;}
	if ( $ARGV[$i] eq "V1"       )   { $USROPT_vert    = 1; $nrd=1;}

	if ( $ARGV[$i] eq "multipdf"  )   
	{ $USROPT_multipdf  = 1; $nrd=1; $MULTIFLAG= "pdf" ; }
	if ( $ARGV[$i] eq "multigif"  )   
	{ $USROPT_multigif  = 1; $nrd=1; $MULTIFLAG= "gif" ; }

	if ( $ARGV[$i] eq "redline0"  )   { $USROPT_redline0  = 1; $nrd=1;}


	# still support legacy eps option, but really should use 'multiFile' opt
	if ( $ARGV[$i] eq "eps"      )   
	{ $USROPT_multieps  = 1; $nrd=1; $MULTIFLAG="eps" ; }
	if ( $ARGV[$i] eq "multieps" )   
	{ $USROPT_multieps  = 1; $nrd=1; $MULTIFLAG="eps" ; }

	if ( $ARGV[$i] eq "eps2gif"  )   
	{ $USROPT_multieps = $USROPT_eps2gif = 1; $nrd=1; $MULTIFLAG="gif"; }

	if ( $ARGV[$i] eq "MAG"      )   { $USRVAL_fscale  = -1; $nrd=1;}

	if ( $ARGV[$i] eq "PRIVATE"  )   { $PRIVATE_FLAG = "PRIVATE"; $nrd=1; }

	for ( $ird=1; $ird <= $nrd; $ird++ ) {
	    $USE_ARG[$iLast+$ird] = 1 ;
	}
	$iLast = $i ;

    } # loop over $NARG


    # make sure that all arguments are recognized;
    # abort on any unknown argument

    $NERR = 0;
    for ( $i=0; $i < $NARG; $i++ ) {
	if ( $USE_ARG[$i] == 0 ) {
	    $NERR++ ;
	    print " WARNING: unknown argument '$ARGV[$i]' \n";
	}
    }
    if ( $NERR > 0 ) { 
	$MSGERR[0] = "unknown arguments above. ";
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ( $USR_INFILE eq "NULL" ) {
	$MSGERR[0] = "must give hbook/root filename as argument." ;
	$MSGERR[1] = "  mkfitplots.pl --h <USR_INFILEName> " ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ( ! (-e $USR_INFILE) ) {
	$MSGERR[0] = "'$USR_INFILE' does not exist." ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ( $USR_INFILE_TYPE == 0 ) {
	$MSGERR[0] = "Undefined USR_INFILE type";
	$MSGERR[1] = "Must be hbook (--h) or root (--r)";
	sntools::FATAL_ERROR(@MSGERR) ;
    }

    @STACK = ( "HORIZONTAL" , "VERTICAL" ) ;

    print "\n" ;
    print "  Make plots for         : $USR_INFILE \n" ;

    if ( $RESIDS_FLAG == 0 ) {
	my $USR_tt = "$USRVAL_tmin to $USRVAL_tmax";
	print "  CID to plot            : $USRVAL_cid  (<=0 -->all)\n";
	print "  Show Tobs range        : $USR_tt  days \n" ;
	print "  Overlay best-fit model : $USROPT_ovfit \n"   ;
	print "  Plot chi2 map          : $USROPT_pchi2 \n"   ;
	print "  Grid option            : $USROPT_grid \n"   ;
	print "  Show peak flux         : $USROPT_pkf    \n"  ;

	if ( $USRVAL_fscale < 0.0 ) 
	{ print "  Fluxscale              : $USRVAL_fscale -> mag \n"  ; }
	else 
	{ print "  Fluxscale              : $USRVAL_fscale \n"  ; }

	print "  Vert page size scale   : $USRVAL_vscale \n"  ;
	print "  Stacking option        : $STACK[$USROPT_vert]   \n"  ;
	print "  Plot marg. pdf         : $USROPT_margpdf \n" ;
	print "  Nrows per page         : $USRVAL_nrow \n" ;
	print "  plotFile per SN(gif,pdf,eps)  : $USROPT_multigif,$USROPT_multipdf,$USROPT_multieps \n";
	if ( $USROPT_multieps ) 
	{  print "  eps2gif                : $USROPT_eps2gif  \n" ; }
    }

} # end of parse_args


# ======================================
sub isArgKey {

    my($rawkey,$ARG) = @_ ;

    # For rawkey = xyz, return true if
    # $ARG = -xyz or --xyz
    # i.e., allow either 1 or 2 dashes.
    #

    my ($key1, $key2);
    $key1 = "-$rawkey";
    $key2 = "--$rawkey";

    if ( $key1 eq $ARG ) { return 1; }
    if ( $key2 eq $ARG ) { return 1; }

    return 0;

} # end of isArgKey

# ======================================
sub set_plotArgs_HBOOK {

    # set global plotArgs

    my ($argov, $argchi, $argt, $argf, $argv, $argrow, $argg );
    my ($argtype, $argcid);

    $argov  = "ovfit=${USROPT_ovfit}";
    $argchi = "pchi2=${USROPT_pchi2}";
    $argt   = "tmin=${USRVAL_tmin} tmax=${USRVAL_tmax}" ;
    $argf   = "fluxscale=$USRVAL_fscale" ;
    $argv   = "vscale=$USRVAL_vscale" ;
    $argrow = "nrowpage=$USRVAL_nrow" ;
    $argcid = "cid=${USRVAL_cid}";
    
    if ( $USRVAL_type != -9 )  
    { $argtype = "type=$USRVAL_type"; }
    else
    { $argtype = "" ; }

    if ( $USROPT_vert == 0 ) { $argv = "$argv" . " vert=0" ; }

    if ( $USROPT_grid == 0 ) 
    { $argg = "grid=0" ; }
    else
    { $argg = "" ; }

    $plotArgs = 
	"$argcid $argt $argov $argchi pkf=${USROPT_pkf} $argf $argv $argg $argrow $argtype" ;

    if ( $USROPT_multieps ) 
    { $plotArgs = "$plotArgs eps=1" ; }


} # end of set_plotArgs_HBOOK

# ======================================
sub set_plotArgs_ROOT {

    my ($NARG, $arg, @listArgs) ;

    $plotArgs = "" ;
    @listArgs = ( "snana" ) ; # special option to root-script

    # for multi-plotFile option, overwrite (clobber) snana option
    if ( $USROPT_multipdf   ) { @listArgs = ( "multipdf"); }
    if ( $USROPT_multigif   ) { @listArgs = ( "multigif"); }

    # options below are appended.
    if ( $USRVAL_fscale < 0 ) { @listArgs = ( @listArgs , "mag"    ); }
    if ( $USROPT_vert  == 1 ) { @listArgs = ( @listArgs , "column" ); }
    if ( $USROPT_ovfit == 0 ) { @listArgs = ( @listArgs , "data"   ); }
    if ( $USROPT_redline0   ) { @listArgs = ( @listArgs , "redline" ); }

    # convert list into comma-separated list needed for root script
    $NARG = 0;
    foreach $arg ( @listArgs) {
	
	if ( $NARG == 0 ) 
	{ $plotArgs = $listArgs[$NARG]; }
	else
	{ $plotArgs = "${plotArgs},$listArgs[$NARG]"; }

	$NARG++ ;
    }

#    die "\n xxx root plotArgs = '$plotArgs' xxx \n" ;

}  # end of set_plotArgs_ROOT


# ============================
sub eps2gif {

    my (@epsList, $epsFile, $gifFile, $idot, $prefix, $cmd);
    print " Convert eps -> gif ... \n";
    @epsList = qx(ls *.eps);
    foreach $epsFile (@epsList) {	   
	$epsFile =~ s/\s+$// ;   # trim trailing whitespace
	$idot    = index($epsFile,".") ;
	$prefix  = substr($epsFile,0,$idot);
	$gifFile = "${prefix}.gif" ;
	$cmd     = "convert $epsFile $gifFile" ;
	qx($cmd);
    }
    qx(rm *.eps);

} # end of eps2gif

# ============================
sub move_plotFiles {

    # move plot files into subdir
    my $sdir = "lcPlots" ;
    my ($Nfile, $ext, @bla );

    $ext = "${MULTIFLAG}" ;

    if ( !(-d $sdir) ) { qx(mkdir $sdir); }
    qx(mv *.${ext} $sdir/);
        
    @bla   = qx(ls $sdir) ;  
    $Nfile = scalar(@bla);
    print " Created  $Nfile  $ext files in $sdir \n";

}  # end of move_plotFiles


# ======================
sub  rename_pdf {

    # created pdf filename has the .his or .root in it;
    # here give a name based on the prefix so that we get
    # the same pdf file name for hbook or root.

    my ($idot, $prefix, $pdf_in, $pdf_out);

    
    $idot = index($USR_INFILE,".");
    $prefix = substr($USR_INFILE,0,$idot);

    $pdf_in   = "${USR_INFILE}.pdf";

    if ( $USRVAL_cid <= 0 ) 
    {  $pdf_out  = "${prefix}_fits.pdf" ; }
    else 
    {  $pdf_out  = "${prefix}_fit${USRVAL_cid}.pdf" ; }

    qx(mv $pdf_in $pdf_out);
    print " rename $pdf_in -> $pdf_out \n" ;

} # end of rename_pdf
