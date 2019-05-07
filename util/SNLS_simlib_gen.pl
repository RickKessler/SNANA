#!/usr/bin/env perl
#
# Created Nov 2011 by R.Kessler
#
# Translate SNLS list of observations  into SNANA-simlib format.
#
# Usage:
#  SNLS_simlib_gen.pl <obsFile>
#
#
# Feb 28, 2012: IDEXPT is now the field number
#
# -----------------------------------------

use IO::Handle;
use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# hard wire model-uncertainty parameters from
# PAW > errcalc#noisefit  survey=SNLS
# PAW > errcalc#zpterrfit survey=SNLS
my @ZPTERR_SIMLIB         = ( -9.0, -9.0, 0.001, 0.001, 0.001, 0.001 ) ;
my $FLUXERR_COR_FILE      = "FLUXERR_COR_SNLS.MAP" ;

# define NEA/(4*PI*PSF^2) from J.Guy e-mail (Feb 29, 2012)
# See GFSEEING in the [private] SNLS headers.
# For g band it is about 1.43*.96 = 1.37 
my @NEA_over_4piPSF2 = ( -9., -9., 1.37, 1.39, 1.40, 1.44 ); 
my $PSFRATIO_SIMLIB  = 0.1 ; # hard wire, then compute PSFSIG2

#hard-wire names of keys in the input file
my $KEY_IDEXPO  = "expo";
my $KEY_FIELD   = "field";
my $KEY_FILT    = "band";
my $KEY_MJD     = "mjd";
my $KEY_PSFSIG  = "seeing";
my $KEY_SKYSIG  = "skysig" ;
my $KEY_ZPTAVG  = "zp" ;
my $KEY_ZPTERR  = "zpe" ;

#misc. hard-wire
my @FILTERLIST = ( "?" , "u", "g", "r" , "i", "z" );
my $SNANA_SIMLIB_FILENAME = "SNLS_3year.SIMLIB";
my $SURVEY  = "SNLS" ;
my $TELE    = "CFHT" ;
my $NFIELD  = 4 ;
my @RA  = ( 0.0, 36.7, 150.1, 215.0, 334.0 );  # Deg (D1 ... D4)
my @DEC = ( 0.0, -4.6, 1.80 , 52.3, -17.8 );   # idem
my $GAIN     = 1.6    ;  # photoelectrons / ADU
my $CCDNOISE = 3.2    ;  # photoelectrons
my $PIXSIZE  = 0.187  ;  # arcsec
my $MJD_MIN  = 52700. ;
my $MJD_MAX  = 54000. ;  # SNLS 3year

# misc. globals
my ($NARG, @MSGERR, $SNLS_OBS_FILENAME ) ;
my ($IVAR_IDEXPO, $IVAR_FIELD, $IVAR_FILT, $IVAR_MJD);
my ($IVAR_PSFSIG, $IVAR_SKYSIG, $IVAR_ZPTAVG, $IVAR_ZPTERR );
my (@CONTENTS, @KEYLINES, $NKEYTOT, $NKEY_ERR );

# -------------------
sub init_snmake ;
sub parse_args ;
sub parse_KEYS ;
sub check_KEY ;
sub open_SIMLIB ;
sub update_SIMLIB ;
sub solve_2ndPSF ;

# ==================== BEGIN MAIN ====================

&init_simlib();

# parse user args
&parse_args();

# parse key in the input SNLS-obs file
&parse_KEYS();

# open output simlib and write header info
&open_SIMLIB();

@CONTENTS = `cat $SNLS_OBS_FILENAME` ;

my ($IDFLD, $line);
for ( $IDFLD=1; $IDFLD <= $NFIELD; $IDFLD++ ) {
    &update_SIMLIB($IDFLD)
}

# finally run the coadd program to combine exposures within each night
# into a single effective exposure.

print "\n";
print " Coadd exposures within each night:\n";
qx(simlib_coadd.exe $SNANA_SIMLIB_FILENAME);
print " Final SNLS simlib: ${SNANA_SIMLIB_FILENAME}.COADD \n";

# ===========================
#
#  END OF MAIN
#
# ===========================


# ===========================
sub init_simlib {

    $NKEY_ERR = 0;
    $IVAR_IDEXPO = -9 ;
    $IVAR_FIELD  = -9 ;
    $IVAR_FILT   = -9 ;
    $IVAR_MJD    = -9 ;
    $IVAR_PSFSIG = -9 ;
    $IVAR_SKYSIG = -9 ;
    $IVAR_ZPTAVG   = -9 ;
    $IVAR_ZPTERR   = -9 ;

} # end of init_snmake


# ===========================
sub parse_args {

    $NARG = scalar(@ARGV) ;

    if ( $NARG == 0 ) {
	$MSGERR[0] = "Must give name of SNLS obs-file as argument." ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    $SNLS_OBS_FILENAME = $ARGV[0];


} # end of parse_args

# =======================
sub parse_KEYS {

    # parse header keys and fill $IVAR_XXX integers
    # with the relative location of variables XXX.

    my ($keyLine, $key, @wdlist, $wd0, $IVAR );

    @KEYLINES = qx(grep '#' $SNLS_OBS_FILENAME);

    $NKEYTOT = scalar(@KEYLINES);
    $IVAR = 0 ;

    foreach $keyLine (@KEYLINES) {
	$IVAR++ ;
	@wdlist     = split(/\s+/,$keyLine) ;
	$wd0        = $wdlist[0] ;
	$key        = substr($wd0, 1, 99);

	if ( $key eq $KEY_IDEXPO ) { $IVAR_IDEXPO = $IVAR ; }
	if ( $key eq $KEY_FIELD  ) { $IVAR_FIELD  = $IVAR ; }
	if ( $key eq $KEY_FILT   ) { $IVAR_FILT   = $IVAR ; }
	if ( $key eq $KEY_MJD    ) { $IVAR_MJD    = $IVAR ; }
	if ( $key eq $KEY_PSFSIG ) { $IVAR_PSFSIG = $IVAR ; }
	if ( $key eq $KEY_SKYSIG ) { $IVAR_SKYSIG = $IVAR ; }
	if ( $key eq $KEY_ZPTAVG ) { $IVAR_ZPTAVG = $IVAR ; }
	if ( $key eq $KEY_ZPTERR ) { $IVAR_ZPTERR = $IVAR ; }
    }

    # make sure that each key is defined; else abort.

    &check_KEY($KEY_IDEXPO, $IVAR_IDEXPO  );
    &check_KEY($KEY_FIELD,  $IVAR_FIELD  );
    &check_KEY($KEY_FILT,   $IVAR_FILT   );
    &check_KEY($KEY_MJD,    $IVAR_MJD    );
    &check_KEY($KEY_PSFSIG, $IVAR_PSFSIG );
    &check_KEY($KEY_SKYSIG, $IVAR_SKYSIG );
    &check_KEY($KEY_ZPTAVG, $IVAR_ZPTAVG );
    &check_KEY($KEY_ZPTERR, $IVAR_ZPTERR );

    if ( $NKEY_ERR > 0 ) {
	$MSGERR[0] = "Missing $NKEY_ERR  header keys in $SNLS_OBS_FILENAME";
	$MSGERR[1] = "See ERROR messages above.";
	sntools::FATAL_ERROR(@MSGERR);
    }
    print " Found all required header keys. \n\n";


} # end of parse_KEYS

sub check_KEY {
    my ($KEY, $IVAR) = @_;

    if ( $IVAR <= 0 ) {
	print " ERROR: could not find header key '$KEY' \n";
	$NKEY_ERR++ ;
    }
    else {
	print "\t Found header key '$KEY' at position $IVAR \n";
    }

} # end of check_KEY


# =================
sub open_SIMLIB {

    my ($ifilt, $i, $line );

    if ( -e $SNANA_SIMLIB_FILENAME ) { qx(rm $SNANA_SIMLIB_FILENAME); }
    open  PTR_SIMLIB , " > $SNANA_SIMLIB_FILENAME";

    print PTR_SIMLIB "SURVEY: $SURVEY   FILTERS: griz   TELESCOPE: $TELE\n";
    print PTR_SIMLIB "COMMENT: created by $0 \n";
    
    print PTR_SIMLIB "\n\n";
   
    my @BLA = `cat $FLUXERR_COR_FILE` ;
    my $NBLA = scalar(@BLA);
    for ($i=0; $i < $NBLA; $i++ ) {
	$line = $BLA[$i];
	$line =~ s/\s+$// ;  # remove trailing blanks
	print PTR_SIMLIB "$line\n";
    }


    print PTR_SIMLIB "\nBEGIN LIBGEN \n";
    print PTR_SIMLIB "\n";
    

} # end of open_SIMLIB

# ===============
sub update_SIMLIB {

    my ($IDFLD) = @_ ;

    # ------

    my (@wdlist, $wd0, $line );
    my ($IDEXPO, $IFIELD, $FIELD, $FILT, $IFILT, $MJD, $SKYSIG);
    my ($PSFSIG1, $PSFSIG2 );
    my ($ZPTERR, $ZPTAVG);
    my ($NOBS, @OBSLINES, $obs, $LIBID, $RADEC );

    $FIELD = "D$IDFLD";
    print " Prepare LIBID for Field $FIELD => ";
   
    $NOBS = 0;
    foreach $line (@CONTENTS) {

	# bail if this  is a header key
	@wdlist     = split(/\s+/,$line) ;
	$wd0        = $wdlist[0] ;
	if ( substr($wd0,0,1) eq '#' ) { next ; }

	$IFIELD   = $wdlist[$IVAR_FIELD-1];
	if ( $IFIELD != $IDFLD ) { next ; }

	$IFILT   = $wdlist[$IVAR_FILT-1];
	if ( $IFILT <= 1 ) { next ; }

	$MJD     = $wdlist[$IVAR_MJD-1];
	if ( $MJD < $MJD_MIN ) { next ; }
	if ( $MJD > $MJD_MAX ) { next ; }

	$IDEXPO   = $wdlist[$IVAR_IDEXPO-1];
	$SKYSIG   = $wdlist[$IVAR_SKYSIG-1];
	$PSFSIG1  = $wdlist[$IVAR_PSFSIG-1];
	$PSFSIG2  = &solve_2ndPSF($IFILT, $PSFSIG1,$PSFRATIO_SIMLIB);
	$ZPTAVG   = $wdlist[$IVAR_ZPTAVG-1];

#	$ZPTERR   = $wdlist[$IVAR_ZPTERR-1]; 
	$ZPTERR   = $ZPTERR_SIMLIB[$IFILT];  # from fit


	$FILT = $FILTERLIST[$IFILT];

	$NOBS++ ;

	my $cMJD = sprintf("%9.3f", $MJD);
	my $cID  = sprintf("%d%6.6d", $IDFLD, $IDFLD );
	my $cG   = sprintf("%4.2f", $GAIN) ;
	my $cRON = sprintf("%4.2f", $CCDNOISE);
	my $cSKY = sprintf("%5.2f", $SKYSIG);
	my $cPSF = sprintf("%5.3f %5.3f %3.1f", 
			   $PSFSIG1, $PSFSIG2, $PSFRATIO_SIMLIB );
	my $cZPT = sprintf("%6.3f %5.3f", $ZPTAVG, $ZPTERR);

	$OBSLINES[$NOBS] = 
	    "S: $cMJD  $cID  $FILT   $cG  $cRON  $cSKY $cPSF   $cZPT  -99.";
    }

    if ( $NOBS == 0 ) {
	print " WARNING: NOBS=$NOBS for FIELD $FIELD \n";
	return ;
    }
    else {
	print " found $NOBS exposures. \n";
    }
    $RADEC = "RA: $RA[$IDFLD]   DECL: $DEC[$IDFLD]";
    $LIBID = $IDFLD;

    print PTR_SIMLIB "# -------------------------------------------- \n";
    print PTR_SIMLIB "LIBID: $LIBID      NOBS: $NOBS\n";
    print PTR_SIMLIB "$RADEC   MWEBV: 0.0  PIXSIZE: $PIXSIZE \n";
    print PTR_SIMLIB "\n";
    print PTR_SIMLIB 
	"#                           CCD  CCD         PSF1 PSF2 PSF2/1       \n";
    print PTR_SIMLIB 
	"#     MJD      IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG \n";


    for ( $obs=1; $obs <= $NOBS ; $obs++ ) {
	print PTR_SIMLIB "$OBSLINES[$obs] \n";
    }
    print PTR_SIMLIB "END_LIBID:  $LIBID \n\n";


} # end of update_SIMLIB


sub solve_2ndPSF {
    my ($ifilt, $PSFSIG1, $PSFRATIO ) = @_ ;

    my ($PSFSIG2, $NEA_RATIO_FIX);
    my ($R, $RR, $RRr,  $Rsig, $Rmin, $Rmax, $Rstep );
    my ($NEA_RATIO, $DIF, $DIFMIN );

    $NEA_RATIO_FIX = $NEA_over_4piPSF2[$ifilt];
    if ( $NEA_RATIO_FIX < 1 || $NEA_RATIO_FIX > 2.0 ) {
	$MSGERR[0] = "Invalid NEA_RATIO = $NEA_RATIO_FIX for ifilt=$ifilt";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $PSFSIG2 = 2. * $PSFSIG1 ;

    $Rstep = .01 ;  $Rmin = 1.8 ; $Rmax=2.5 ;
    $DIFMIN = 99999. ; $Rsig = -9.0 ;
    for ( $R = $Rmin ; $R < $Rmax; $R += $Rstep) {
	$RR  = $R*$R ;
	$RRr = $RR * $PSFRATIO ;
	my $dum1 = 1.0 + $RR ;
	my $dum2 = (1.0 + $RRr)**2 ;
	my $dum3 = $RR*((1.0 + $PSFRATIO)**2) ;
	$NEA_RATIO = ($dum1*$dum2)/($dum3+$dum2);
	$DIF = abs($NEA_RATIO-$NEA_RATIO_FIX);
	if ( $DIF < $DIFMIN ) {
	    $DIFMIN = $DIF ; $Rsig = $R ;
	}

    }

    if ( $Rsig > $Rmax - $Rstep || $Rsig < $Rmin+$Rstep ) {
	$MSGERR[0] = "Rsig = $Rsig is at edge of grid.";
	$MSGERR[1] = "R-grid: $Rmin to $Rmax with Rstep=$Rstep" ;
	$MSGERR[2] = "ifilt=$ifilt  PSFSIG1=$PSFSIG1  PSFRATIO=$PSFRATIO";
	sntools::FATAL_ERROR(@MSGERR);	
    }

    $PSFSIG2 = $Rsig * $PSFSIG1;
#    die "\n xxx PSFSIG[1,2] = $PSFSIG1 , $PSFSIG2 xxxxx \n\n";
    return $PSFSIG2 ;

} # end of solve_2ndPSF
