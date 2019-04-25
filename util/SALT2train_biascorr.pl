#!/usr/bin/perl
#--------------
# script designed to do bias corrections
# JMosher 03/22/13
#### 
# required inputs:
# 1) fitres file whose distances are to be corrected ["TEST" set]
# 2) bias-correction fitres file, simulated to determine biases ["FIX" set]
# 3) output file name ["CORR"]
# 4) correction type -- current options are R13 and P11
#    ** if using P11, must also enter sim alpha and sim beta **
# optional inputs:
# 1) binsize -- width of z bin over which to avg bias corrections
# 2) survey flag and binsize -- if turned on the script will repeat bias corrections
#    N+1 times, first with the full sample, then once for each of the N unique survey 
#    subsamples. 
#    
####
# procedure:
# -- TEST distances are averaged as a function of redshift
# -- bias corrections are calculated for the same redshift bins as TEST data
# -- CORRECTED average distances <mu>_{CORR} and average distance errors <mue>_CORR 
#    are calculated as follows:
#     --> <mu>_{CORR}  =  <mu>_{TEST} + <corr>_{FIX}
#     --> <mue>_{CORR} =  sqrt(<mue>_{TEST}^2 + <corre>_{FIX}^2) 
# -- CORRECTED average distances are output in redshift bins
####
# comments on averaging DM in redshift bins:
# -- average redshifts are weighted by DM uncertainty
# -- average redshift error is sqrt(1.0/sum(1/ze^2))
# -- average DM is < DM_i - SIMMU_i > + REF_MU(<z>)
# -- REF_MU is obtained by linearly interpolation from
#       unique (SIM_Z, SIM_MU) values grabbed from fitres file(s) 
####
# Correction types:
#  R13 - as described in the Kessler et al. 2013 mag smear paper
#  P11 - inspired by SNLS3 malmquist bias analysis (Perrett et al. 2011)
####
# Default binsize is 0.05
####
# 
# ****** HISTORY ********
#
# May 26 2013: For P11 bias corrections, added checks. If detected SN sample
#              is empty, abort. If BINSIZE is too small (less than 
#              $Nlimit SNe in bin) increase BINSIZE by BINSIZEstep. 
#              $Nlimit, $BINSIZE, $BINSIZEstep are all (optional) input variables. 
#              Checks may be extended to R13 corrections, but haven't been
#              at this time.(JLM)
#
# May 31 2013: To aid user debugging, added fail checks on 
#              sntools::rdfitres_getcol calls. (JLM)
#
# Jun 25 2013: Added recalculation of ZMIN, ZMAX to P11 bias correction binsize check.
#              Avoids binsize check bug that happens with small numbers of FIX SNe
#              (i.e. during QUICKTEST runs) (JLM). 
#
# Sep 20 2013: Extensive overhaul to fix uncertainty calculations by BINNING all DM data. 
#              Changes to many subroutines. Introduction to code has been extensively 
#              rewritten. (JLM)
#
# Oct 08 2013: Added -survey <survey binsize> option to produce survey-by-survey results
#              in addition to full TEST sample results. Cleaned up (make_R13, make_P11
#              cosmo array consistency, ) (JLM)
#
# Mar 19 2014: Updated sub build_fitres_array(@) to take into account new fitres keywords
#              SIMSED_S2x1 and SIMSED_S2c. (JLM)
#
#########################

use List::Util qw(min max);
use FindBin qw($Bin);
use lib "$Bin";
use sntools;

use File::Basename;
use File::Copy;
use POSIX;
use Cwd;
use strict;

# initialization subroutines
sub get_command_args(@);
sub survey_setup(@);
# key analysis subroutines
sub make_P11_corrarray(@);
sub make_R13_corrarray(@);
sub make_outfile(@);
# most important utilities
sub build_cosmo_array(@);
sub get_cosmomu(@);
sub build_fitres_array(@);
sub rebin_data(@);
sub grab_corr(@);
sub lininterp(@);
# bookkeeping utilities
sub get_zrange(@);
sub get_binzrange(@);
sub check_z(@);
sub getcol_failcheck(@);
sub check_binsize(@);
sub check_P11binsize(@);
sub add_to_filename(@);

#### DECL ####

my (@MSGERR, $CODEname, $VERBOSE);
my ($SURVFLAG, $SURVBINSIZE);
my (@SURVEYBINSIZES, @UNIQSURVEYS);
my (@SURVOUTFILES, @SURVDATFILES);
my ($isurv, $NSURV, $survval);
my ($OUTTYPE, $P11TYPE, $ZERRFLAG);
my ($datfile, $bcorrfile, $goutfile);
my ($gdatfilebin, $datfilebin, $outfile);
my (@ZDAT, @ZBCORR);
my ($NCOL_DAT, $NLINE_DAT);
my ($NCOL_BCORR, $NLINE_BCORR);
my ($CORRTYPE, $SIMA, $SIMB);
my (@ARR_CORRDATA, @ARR_TESTDATA);
my ($BINSIZEstep, $Nlimit, $BINSIZE);
my ($BZMIN, $BZMAX, $DZMIN, $DZMAX);
my ($GZMIN, $GZMAX, $GBINSIZE); #MASTER LOOP VARS
my ($ZMIN, $ZMAX); #SURVEY LOOP Z 
my ($DUMPFLAG);

my (@wdLine, $temp_mean, $temp_errmean);
my ($item1, $item2);

#### MAIN ####

$CODEname = "SALT2train_biascorr.pl";

########
######## SETTING UP VARIABLES
########

## get file names, optional args
## perform dummy checks
 get_command_args(@ARGV); 

## get Z columns, 
## use full TEST BINSIZE to set full redshift range for corrections

 sntools::rdfitres_getcol($datfile, \@ZDAT, "Z", "VARNAMES:", "SN:", 0, $VERBOSE);
 sntools::rdfitres_getcol($bcorrfile, \@ZBCORR, "Z", "VARNAMES:", "SN:", 0, $VERBOSE);

 $BINSIZE = $GBINSIZE;
 ($DZMIN, $DZMAX) = get_zrange(\@ZDAT); #check TEST data z limits
 ($BZMIN, $BZMAX) = get_zrange(\@ZBCORR); #check FIX data z limits
 ($GZMIN, $GZMAX) = check_z($DZMIN, $DZMAX, $BZMIN, $BZMAX, -1); #set global ZMIN, ZMAX (survey = -1 => ALL data)

## set up survey-related variables 
## if surveyflag == 0, $NSURV = 1, only FULL TEST SET is analyzed
 $NSURV = survey_setup($datfile, \@SURVEYBINSIZES, \@UNIQSURVEYS, 
    \@SURVOUTFILES, \@SURVDATFILES);

########
######## RUNNING BIAS CORRECIONS
########

## loop over surveys
 for ( $isurv=0; $isurv<$NSURV; $isurv++ ) {

    # initialize survey-dependent variables
    $survval = $UNIQSURVEYS[$isurv];
    $BINSIZE = $SURVEYBINSIZES[$isurv];
    $outfile = $SURVOUTFILES[$isurv]; 
    $datfilebin = $SURVDATFILES[$isurv];
    $ZMIN = $GZMIN;
    $ZMAX = $GZMAX;
    @ARR_TESTDATA = ();

    print "\n\tStarting survey $survval ... \n";
    print "\t\t Survey binsize is $BINSIZE \n";
    print "\t\t Survey outfile is $outfile, BINNED dat file is $datfilebin \n";
    print "\t\t Survey initial ZMIN, ZMAX are $ZMIN, $ZMAX \n\n";

    ### 
    ### GRAB NEEDED DATA FROM DATFILE, PACKAGE IN ARRAY
    ### format: <z><ze><mu><mue><sim_mu><effmask><survey id><cid>
    ($DZMIN, $DZMAX) = build_fitres_array($datfile, "DAT", \@ARR_TESTDATA, $survval);

    ## make the correction array
    ## steps: build ref DM array (build_cosmo_array), 
    ##        load fitres data (build_fitres_array),
    ##        rebin data (rebin_data).
    if ($CORRTYPE eq "R13") { 
       make_R13_corrarray(\@ARR_CORRDATA, $survval);
    } elsif ($CORRTYPE eq "P11") {
       make_P11_corrarray(\@ARR_CORRDATA, $survval);
    } 

    # write corrected data to output file
    make_outfile(\@ARR_TESTDATA, $survval);
    
} # end loop over surveys

print "\n $CODEname finished successfully \n\n";

#### END MAIN ####

sub get_command_args(@)
{

    my $NARG = scalar(@_);

    my ($j, $comment, $SUBname);

    $SUBname = "get_command_args";

    $VERBOSE = 0; # set default verbosity type
    $P11TYPE = "APPX"; # set default P11 correction type
    $GBINSIZE = 0.05; # SET DEFAULT BINSIZE
    $BINSIZEstep = 0.01; # SET DEFAULT BINSIZEstep
    $Nlimit = 3; # SET DEFAULT Nlimit

    #check that user gave required inputs
    print "\nStarting SALT2train_biascorr.pl ... \n";

    if ($NARG < 8) {
	# 8 = minimum entries x 2
        $MSGERR[0] = "\nNot enough command line inputs. \n";
	$MSGERR[1] = "Must enter -datfile <name> -bcorrfile <name> -outfile <name> \n";
	$MSGERR[2] = "           -corrtype <name> \n";
	$MSGERR[3] = "Valid corrtypes are R13, P11 \n";
	$MSGERR[4] = "P11 REQUIRES -sima <val> and -simb <val> \n";
	$MSGERR[5] = "Optional flags: -binsize <val> -binsizestep <val> \n";
	$MSGERR[6] = "                -Nlimit <val>  -P11type <val> \n";
	$MSGERR[7] = "                -survey <binsize>  -verbose \n";
	$MSGERR[8] = "Default redshift binsize=0.05 \n";
	$MSGERR[9] = "Valid P11types are TRUE, APPX(default) \n";
	$MSGERR[10] = "Using -survey will produce extra SURV<#>.fitres \n";
	$MSGERR[11] = "  files to be created \n";
	$MSGERR[12] = "Using -verbose will dump extra info to log \n";
	$MSGERR[13] = "Check inputs and try again. \n";
	sntools::FATAL_ERROR(@MSGERR);
    }

    for ($j=0; $j<$NARG; $j++) {
	###print "xxx ARGV[$j] is $ARGV[$j] \n";
	if ("$ARGV[$j]" eq "-datfile"){
	    $datfile = $ARGV[$j+1];
	    ###print "xxx found -datfile. val is $datfile \n";
	    $j++;
	} elsif ("$ARGV[$j]" eq "-bcorrfile"){
	    $bcorrfile = $ARGV[$j+1];
	    ###print "xxx found -bcorrfile. val is $bcorrfile \n";
	    $j++;
	} elsif ("$ARGV[$j]" eq "-outfile"){
	    $goutfile = $ARGV[$j+1];
	    ###print "xxx found -outfile. val is $outfile \n";
	    $j++;
	} elsif ("$ARGV[$j]" eq "-P11type"){
	    $P11TYPE = $ARGV[$j+1];
	    $j++;
	} elsif ("$ARGV[$j]" eq "-sima"){
	    $SIMA = $ARGV[$j+1];
	    $j++;
	} elsif ("$ARGV[$j]" eq "-simb"){
	    $SIMB = $ARGV[$j+1];
	    $j++;
        } elsif ("$ARGV[$j]" eq "-corrtype"){
	    $CORRTYPE = $ARGV[$j+1];
	    $j++;
	} elsif ("$ARGV[$j]" eq "-binsize"){
	    $GBINSIZE = $ARGV[$j+1];
	    $j++;
	} elsif ("$ARGV[$j]" eq "-binsizestep"){
	    $BINSIZEstep = $ARGV[$j+1];
	    $j++;
	} elsif ("$ARGV[$j]" eq "-Nlimit"){
	    $Nlimit = $ARGV[$j+1];
	    $j++;
	} elsif ("$ARGV[$j]" eq "-survey"){
	    $SURVFLAG = 1;
	    $SURVBINSIZE = $ARGV[$j+1];
	    $j++;
	} elsif ("$ARGV[$j]" eq "-verbose"){
	    $VERBOSE = 1;
	} else{
	    $MSGERR[0] = "In subroutine $SUBname";
	    $MSGERR[1] = "Unused command line flag $ARGV[$j]";
	    $MSGERR[2] = "Check input syntax";
	    sntools::FATAL_ERROR(@MSGERR);
	}
    }

    ### DUMMY CHECKS
    if ($CORRTYPE ne "P11" && $CORRTYPE ne "R13" ) { 
	$MSGERR[0] = "In subroutine $SUBname";
	$MSGERR[1] = "INVALID -corrtype: $CORRTYPE";
	$MSGERR[2] = "VALID choices are: P11, R13";
	$MSGERR[3] = "Check inputs and try again";
	$MSGERR[4] = "ABORT";
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ($P11TYPE ne "TRUE" && $P11TYPE ne "APPX" ) { 
	$MSGERR[0] = "In subroutine $SUBname";
	$MSGERR[1] = "INVALID -P11type: $P11TYPE";
	$MSGERR[2] = "VALID choices are: TRUE, APPX(default)";
	$MSGERR[3] = "Check inputs and try again";
	$MSGERR[4] = "ABORT";
	sntools::FATAL_ERROR(@MSGERR);
    }

    if ($CORRTYPE eq "P11" && ($SIMA eq "" || $SIMB eq "" )){
	$MSGERR[0] = "In subroutine $SUBname";
	$MSGERR[1] = "-corrtype P11 REQUIRES sima, simb";
	$MSGERR[2] = "Change -corrtype or use -sima, -simb flags";
	$MSGERR[3] = "ABORT";
	sntools::FATAL_ERROR(@MSGERR);
    }

    $comment = "$CODEname -datfile ";
    sntools::check_fileexist($datfile, $comment);
    $comment = "$CODEname -bcorrfile ";
    sntools::check_fileexist($bcorrfile, $comment);

    ## MAKE BINNED DATFILE NAME, KEEPING ORIGINAL EXTENSION
    $gdatfilebin = add_to_filename($datfile, "BINNED");

    ## OUTPUT SUMMARY FOR LOG FILE
    print "\n*************************************************\n";
    print "\nBIAS CORRECTION SUMMARY \n\n";
    print "\n*************************************************\n";
    print "\n\tBias correction type is $CORRTYPE \n";
    if ( $CORRTYPE eq "P11" ) {
	print "\t\tP11 correction type is $P11TYPE \n";
    }
    print "\n\tRedshift binsize is $GBINSIZE \n";
    print "\n\tIf less than $Nlimit SNe per bin, \n";
    print "\t\tBINSIZE will be increased by $BINSIZEstep \n\n";
    print "\n\tObtaining bias corrections from file $bcorrfile \n";
    print "\n\tApplying bias corrections to file $datfile \n";
    print "\n\tBinned corrected distances will be written to $outfile \n";
    print "\n\tBinned uncorrected distances will be written to $datfilebin \n";
    print "\n\tSurvey setting is $SURVFLAG (1=on, 0=off) \n\n";
    if ( $SURVFLAG == 1 ) {
	print "\n\tSurvey binsize will be $SURVBINSIZE \n\n";
    }
    print "\n\tVerbose setting is $VERBOSE (1=on, 0=off) \n\n";
    print "\n*************************************************\n";
    
} #end get_command_args(@)

sub make_R13_corrarray(@)
{
    my ($p_R13corrarray, $survid) = (@_);
    my (@cosmo, @corr, @sortcorrALL, @sortcorr);
    my ($iarr, $narr, $tempcorr, $tempcorrerr);
    my ($iz, $NZ, $binz, $zlow, $zhigh, $refmu);
    my ($mZ, $mZe, $tcorr, $tcorre, $mSMU, $mSMUe, $mID, $nBIN);
    my ($survey_bzmin, $survey_bzmax);

    print "\n\tGenerating R13 corrections ... \n";

    #########
    #########
    #### R13 BIAS CORRECTION: 
    ####
    ####   CORRECTIONS: mu_{corr,TEST}(zbin) = mu_{TEST}(zbin) + Dmu_{FIX}(zbin)
    ####   
    ####   FIX data: simulated with TEST alpha, beta
    ####             includes only detected ("D") SNe
    ####
    ####   CALCULATING CORRECTIONS Dmu_{FIX}(zbin):
    ####
    ####     dmu_i = simmu_i - mu_i
    ####
    ####     Dmu(zbin)  = < dmu >(zbin)
    ####     
    ####     Dmue(zbin) = sig_dmu = sqrt(1.0/C), where
    ####     
    ####         C is defined as sigma(i=1 to N_FIX) 1.0/(dmu_i^2)
    ####     
    #########
    #########    

    ### BUILD COSMOLOGY ARRAY
    build_cosmo_array($bcorrfile, \@cosmo); 

    ### GRAB NEEDED DATA FROM BCORRFILE, DETERMINE R13 MU' VALUES, PACKAGE IN ARRAY
    ($survey_bzmin, $survey_bzmax) = build_fitres_array($bcorrfile, "R13", \@corr, $survid);

    ### REDETERMINE ZMIN, ZMAX VALUES FOR THIS SURVEY
    ($ZMIN, $ZMAX) = check_z($DZMIN, $DZMAX, $survey_bzmin, $survey_bzmax, $survid);

    ### SORT ARRAY, SELECT FOR APPROPRIATE DATA TYPES
    @sortcorrALL = sort { $a->[0] <=> $b->[0] } @corr; #arrange in order of z
    @sortcorr = grep { $_->[5] == 3 } @sortcorrALL; #select effmask==3

    ### REBIN DATA USING SURVEY BINSIZE
    print "\n\tRebinning R13 bias corr data ... \n";
    $NZ = 0;
    for ( $iz = $ZMIN; $iz <= $ZMAX ; $iz=$iz+$BINSIZE ) {
	$binz = $iz;
        ($mZ, $mZe, $tcorr, $tcorre, $mSMU, $mSMUe, $mID, $nBIN, $zlow, $zhigh, $refmu) 
	    = rebin_data(\@sortcorr, \@cosmo, $iz, $BINSIZE, "R13");
	$$p_R13corrarray[$NZ] = [$binz, $mZ, $tcorr, $tcorre, $nBIN, $nBIN, $nBIN];
	$NZ++;
    }
    

} #end make_R13_corrarray(@)

sub check_binsize(@)
{

    #### UTILITY to determine and return the number of SNe in a given 
    #### redshift bin. 

    my ($p_corrarr, $zval, $zrange) = (@_);
    my (@subarray, $zhigh, $zlow, $ntemp);

    $zhigh = $zval + $BINSIZE/2.0;
    $zlow = $zval - $BINSIZE/2.0;

    @subarray = grep { $_->[0] > $zlow && $_->[0] <= $zhigh } @$p_corrarr;
    $ntemp = @subarray;

    return ($ntemp);

} #end check_binsize(@)

sub rebin_data(@)
{
    ####
    #### UTILITY to calculate redshift-binned averages of fitres array data 
    ####

    my ($p_datarr, $p_cosmoarr, $zval, $zrange, $type) = (@_);

    my (@subarray, $zhigh, $zlow);
    my (@zcol, @zecol, @mu1col, @muecol, @mu2col, @mu3col); 
    my (@S2x1col, @S2ccol, @NDOFcol, @SNRMAX1col, @x1col, @x1ecol);
    my (@ccol, @cecol, @mBcol, @mBecol, @FITPROBcol, @idcol);
    my (@altmu1col, @altmu2col, $iarr, $narr);
    # GENERIC OUTPUT NAMES
    my ($zout, $zeout, $mu1out, $mue1out);
    my ($mu2out, $mue2out);
    # MU TYPE NAMES
    my ($meanzMU, $sigzMU);
    my ($meanMU1, $sigMU1, $meanMU2, $sigMU2);
    my ($meanALTMU1, $sigALTMU1, $meanALTMU2, $sigALTMU2);
    # others
    my ($meanzze, $sigzze);
    my ($meanx1, $sigx1, $meanc, $sigc);
    my ($meanmB, $sigmB);
    my ($meanid, $rmsid, $NBIN, $Nindex);
    my ($meanS2x1, $rmsS2x1, $meanS2c, $rmsS2c);
    my ($meanNDOF, $rmsNDOF, $meanSNRMAX1, $rmsSNRMAX1);
    my ($meanFITPROB, $rmsFITPROB);
    my ($binmu);
    my ($arithz);
    
    my $SUBname = "rebin_data";

    ########
    ######## SUB rebin_data: Rebins SN data packaged by sub build_fitres_array
    ########                    (special version to test binning)           
    ########
    ######## build_fitres_array data is [FOR BIAS CORR][FOR TEST INFO]:
    ########   z, ze, mu1, mue, mu2, type[D/ND], survey id, cid, mu3,  
    ########          S2x1, S2c, NDOF, SNRMAX1, x1, x1e, c, ce, mB, mBe, FITPROB                  
    ########
    ######## mu1, mu2 values depend on correction type (see sub build_fitres_array for more),
    ######## mu3 is always simmu
    ######## 

    # determine boundaries of redshift bin
    ($zlow, $zhigh) = get_binzrange($zval);

    # grab subset of data falling within redshift bin
    @subarray = grep { $_->[0] > $zlow && $_->[0] <= $zhigh } @$p_datarr;
    $Nindex = @subarray;
    if ( $VERBOSE ) {
	print "\t\tNumber of elements in rebin_data subarray (>$zlow, <=$zhigh): $Nindex \n";
    }
    if ( $Nindex == 0 ) {
	$MSGERR[0] = "$SUBname failure:";
	$MSGERR[1] = "No subarray elements between zlow:$zlow";
	$MSGERR[2] = "and zhigh:$zhigh.";
	$MSGERR[3] = "Check input and try again";
	sntools::FATAL_ERROR(@MSGERR);
    }
	
    # grab needed columns from data array
    @zcol = map $_->[ 0 ], @subarray;
    @zecol = map $_->[ 1 ], @subarray;
    @mu1col = map $_->[ 2 ], @subarray;
    @muecol = map $_->[ 3 ], @subarray;
    @mu2col = map $_->[ 4 ], @subarray;
    @idcol = map $_->[ 6 ], @subarray;
    @mu3col = map $_->[ 8 ], @subarray;
    @S2x1col = map $_->[ 9 ], @subarray;
    @S2ccol = map $_->[ 10 ], @subarray;
    @NDOFcol = map $_->[ 11 ], @subarray;
    @SNRMAX1col = map $_->[ 12 ], @subarray;
    @x1col = map $_->[ 13 ], @subarray;
    @x1ecol = map $_->[ 14 ], @subarray;
    @ccol = map $_->[ 15 ], @subarray;
    @cecol = map $_->[ 16 ], @subarray;
    @mBcol = map $_->[ 17 ], @subarray;
    @mBecol = map $_->[ 18 ], @subarray;
    @FITPROBcol = map $_->[ 19 ], @subarray;
    
    # subtract SIMMU (MU3) from MU1 and MU2 to make ALTMU1 and ALTMU2
    # current types are :
    #
    #  R13: MU1 = $simmu_dat - $mu_dat                    MU2 = $mu_dat
    #  P11: MU1 = $mu_dat + ($da*$x1_dat) + ($db*$c_dat)  MU2 = $mu_dat
    #  DAT: MU1 = $mu_dat                                 MU2 = $simmu_dat
    #
    $narr = @mu3col;
    for ($iarr = 0; $iarr < $narr; $iarr++) {
	$altmu1col[$iarr] = $mu1col[$iarr] - $mu3col[$iarr];
	$altmu2col[$iarr] = $mu2col[$iarr] - $mu3col[$iarr];
    }

    # get averages
    ($meanzze, $sigzze, $NBIN) = sntools::weighted_avg(\@zcol, \@zecol, $SUBname, "zz");
    ($meanzMU, $sigzMU, $NBIN) = sntools::weighted_avg(\@zcol, \@muecol, $SUBname, "zMU");

    ($meanMU1, $sigMU1, $NBIN) = sntools::weighted_avg(\@mu1col, \@muecol, $SUBname, "MU1");
    ($meanMU2, $sigMU2, $NBIN) = sntools::weighted_avg(\@mu2col, \@muecol, $SUBname, "MU2");

    ($meanALTMU1, $sigALTMU1, $NBIN) = sntools::weighted_avg(\@altmu1col, \@muecol, $SUBname, "ALTMU1");
    ($meanALTMU2, $sigALTMU2, $NBIN) = sntools::weighted_avg(\@altmu2col, \@muecol, $SUBname, "ALTMU2");

    ($meanx1, $sigx1, $NBIN) = sntools::weighted_avg(\@x1col, \@x1ecol, $SUBname, "x1");
    ($meanc, $sigc, $NBIN) = sntools::weighted_avg(\@ccol, \@cecol, $SUBname, "c");
    ($meanmB, $sigmB, $NBIN) = sntools::weighted_avg(\@mBcol, \@mBecol, $SUBname, "mB");

    ($meanid, $rmsid, $NBIN) = sntools::arith_avg(\@idcol, $SUBname, "id");
    ($arithz, $NBIN, $NBIN) = sntools::arith_avg(\@zcol, $SUBname, "z");

    ($meanS2x1, $rmsS2x1, $NBIN) = sntools::arith_avg(\@S2x1col, $SUBname, "S2x1");
    ($meanS2c, $rmsS2c, $NBIN) = sntools::arith_avg(\@S2ccol, $SUBname, "S2c");
    ($meanNDOF, $rmsNDOF, $NBIN) = sntools::arith_avg(\@NDOFcol, $SUBname, "NDOF");
    ($meanSNRMAX1, $rmsSNRMAX1, $NBIN) = sntools::arith_avg(\@SNRMAX1col, $SUBname, "SNRMAX1");
    ($meanFITPROB, $rmsFITPROB, $NBIN) = sntools::arith_avg(\@FITPROBcol, $SUBname, "FITPROB");


    # MAKE OUTPUT VARIABLES
    #
    #  R13: MU1    = $simmu_dat - $mu_dat                    MU2 = $mu_dat
    #       MUALT1 = MU1 - $simmu_dat = -$mu_dat             MUALT2 = MU2 - $simmu_dat
    #
    #  P11: MU1    = $mu_dat + ($da*$x1_dat) + ($db*$c_dat)  MU2 = $mu_dat
    #       MUALT1 = MU1 - $simmu_dat                        MUALT2 = MU2 - $simmu_dat
    #
    #  DAT: MU1    = $mu_dat                                 MU2 = $simmu_dat
    #       MUALT1 = MU1 - $simmu_dat                        MU2 = MU2 - $simmu_dat
    #

    $zout = $meanzMU; 
    $zeout = $sigzze;
    $binmu = get_cosmomu($p_cosmoarr, $zout, $zlow, $zhigh);
    if ( $type eq "R13" ) {
	$mu1out = $meanMU1;
	$mue1out = $sigMU1;
    } else {
	$mu1out = $meanALTMU1 + $binmu;
	$mue1out = $sigALTMU1;
    }
    $mu2out = $meanALTMU2 + $binmu;
    $mue2out = $sigALTMU2;

    if ( $type eq "DAT" ) {
	#### FULL LIST OF OUTPUT VARIABLES REQUIRED FOR BINNED TEST DATA FILES
	return ($zout, $zeout, $mu1out, $mue1out, $mu2out, $mue2out, $meanid , $NBIN, $zlow, $zhigh, $binmu, 
	        $meanx1, $sigx1, $meanc, $sigc, $meanS2x1, $meanS2c, $meanNDOF, $meanSNRMAX1, $meanFITPROB);
    } else {
	#### ONLY SUBSET OF OUTPUT VARIABLES REQUIRED FOR CORR DATA FILES
	return ($zout, $zeout, $mu1out, $mue1out, $mu2out, $mue2out, $meanid , $NBIN, $zlow, $zhigh, $binmu); #SEP16
    }

} #end rebin_data(@)

sub make_outfile(@)
{
    
    my ($p_TESTDAT, $survey) = (@_);

    my ($zval, $corr, $dcorr, $ncorr, $binz);
    my ($header1, $header2, $header3, $header4, $header);
    my ($newmu, $newdmu, $nval, $ival);
    my (@sortTESTDAT, @BINTESTDAT, $ndat, $idat, $iz, $snid);
    my ($zerr, $mu, $muerr, $simmu, $simmuerr, $meanid, $nbindat);
    my ($x1, $x1e, $c, $ce, $S2x1, $S2c, $NDOF, $SNRMAX1, $FITPROB);
    my ($mB, $mBe, @cosmo);
    my ($dummyid, $zhigh, $zlow, $refmu, $corrz);
    # SURVEY RELATED VARIABLES
    my ($temp_zlow, $temp_zhigh, $temp_narray);

    print "\n\tCorrecting $datfile MU \n";

    ### BUILD COSMOLOGY ARRAY
    build_cosmo_array($datfile, \@cosmo);

    ### SORT DATA BY ASCENDING Z
    @sortTESTDAT = sort { $a->[0] <=> $b->[0] } @{$p_TESTDAT};
    $temp_narray = @sortTESTDAT; 

    #########
    #########
    #### CALCULATING BIAS CORRECTIONS:
    ####  outfile: 
    ####      mu_i,corr(z)  = mu_i(z) + Dmu(zbin) 
    ####      mue_i,corr(z) = sqrt(mue_i(z)^2 + Ncorr*Dmue(zbin)^2) 
    ####
    ####      here, Ncorr = nD(zbin) for P11
    ####            Ncorr =  N(zbin) for R13
    #### 
    ####      AUG 30 2013: this technique gives WRONG mue_i,corr(z)
    ####                   IGNORES correlations between mu_i(z) and Dmu(zbin)
    ####                   fix in progress...
    ####
    ####  outfile.binned:
    ####      <mu,corr>(zbin) = <mu,test>(zbin) + Dmu(zbin)
    ####       mue,corr(zbin) = sqrt(sigma(<mu,corr>,zbin)^2 + sigma(Dmu,zbin)^2)
    ####       
    ####       here, both sigmas are ERRORS ON THE MEAN VALUE (i.e. prop to 1/sqrt(N))
    ####
    ####       AUG 30 2013: using binned file requires minor changes to SALT2mu.exe
    ########
    ########

    ########
    ######## Create BINNED outfile
    ########

    print "\tWriting binned corrected distances to $outfile \n";
    print "\tAlso writing binned uncorrected distances to \n";
    print "\t\t $datfilebin \n";

    #### OPEN OUTPUT FILES, WRITE HEADER TO FILES

    ### CORR file header
    $header1 = "NVAR: 17 \n";
    $header2 = "VARNAMES: CID Z ZERR MU MUERR IDAVG IDSURVEY MUSHIFT ";
    $header3 = "MUSHIFTERR SIM_MU SIMMUERR ORIGMU ORIGMUERR NBIN NBIAS SIMA SIMB \n\n";
    $header = $header1 . $header2 . $header3;

    open (MYFILE, ">$outfile");
    print MYFILE $header;

    ### DAT file header
    $header1 = "NVAR: 28 \n";
    $header2 = "VARNAMES: CID Z ZERR MU MUERR IDAVG IDSURVEY MUSHIFT ";
    $header3 = "MUSHIFTERR SIM_MU SIMMUERR ORIGMU ORIGMUERR NBIN NBIAS SIMA SIMB ";
    $header4 = "S2x1 S2c x1 x1ERR c cERR mB mBERR NDOF SNRMAX1 FITPROB \n\n";
    $header = $header1 . $header2 . $header3 . $header4;

    open (CONTROLFILE, ">$datfilebin");
    print CONTROLFILE $header;

    #### GENERATE BINNED TEST SN DATA WITH BINS SAME AS BIASCORR BINS
    ### TESTDAT format: <z><ze><mu><mue><sim_mu><effmask><survey id><cid>
    $dummyid = -9;
    $snid = ($survey+1)*100;
    for ( $iz =$ZMIN; $iz <= $ZMAX ; $iz=$iz+$BINSIZE ) {
	$binz = $iz;
        ($zval, $zerr, $mu, $muerr, $simmu, $simmuerr, $meanid, $nbindat, $zlow, $zhigh, $refmu, 
                $x1, $x1e, $c, $ce, $S2x1, $S2c, $NDOF, $SNRMAX1, $FITPROB) = 
	        rebin_data(\@sortTESTDAT, \@cosmo, $iz, $BINSIZE, "DAT");
         #printf "xxxTESTDATINFO <meanz> <meanzerr> <mu_ND> <mue_ND> <nND> <zlow> <zhigh> <muref>\n";
         #printf "xxxTESTDAT %0.4f %0.4g %0.4f %0.4f %d %0.4f %0.4f %0.4f \n", 
	         #$zval, $zerr, $mu, $muerr, $nbindat, $zlow, $zhigh, $refmu;
	($corr, $dcorr, $ncorr, $corrz) = grab_corr(\@ARR_CORRDATA, $binz);
	$newmu = $mu+$corr;
	$newdmu = sqrt($muerr*$muerr+$dcorr*$dcorr);
	#printf "xxxCORRDAT(meanz, meanze, mucorr, mucorre) %0.4f %0.4g %0.4f %0.4f %0.4f %0.4f %0.4f \n\n", 
	#     $zval, $zerr, $newmu, $newdmu, $corrz, $corr, $dcorr;
	
	printf MYFILE "SN: %d  %0.4f  %0.4g  ", $snid, $zval, $zerr;
	printf MYFILE "%0.4f  %0.4g  %0.4f  %d    ", $newmu, $newdmu, $meanid, $dummyid; 
	printf MYFILE "%0.4f  %0.4f  %0.4f  %0.4f  ", $corr, $dcorr, $simmu, $simmuerr;
	printf MYFILE "%0.4f  %0.4f  %d    %d    ", $mu, $muerr, $nbindat, $ncorr;
	printf MYFILE "%0.4f  %0.4f \n", $SIMA,  $SIMB;

	$newmu = $mu;
	$newdmu = sqrt($muerr*$muerr);
	
	printf CONTROLFILE "SN: %d  %0.4f  %0.4g  ", $snid, $zval, $zerr;
	printf CONTROLFILE "%0.4f  %0.4g  %0.4f  %d    ", $newmu, $newdmu, $meanid, $dummyid; 
	printf CONTROLFILE "%0.4f  %0.4f  %0.4f  %0.4f  ", $corr, $dcorr, $simmu, $simmuerr;
	printf CONTROLFILE "%0.4f  %0.4f  %d    %d    ", $mu, $muerr, $nbindat, $ncorr;
	printf CONTROLFILE "%0.4f  %0.4f  ", $SIMA,  $SIMB;
	printf CONTROLFILE "%0.4f  %0.4f  %0.4f  %0.4f  ", $S2x1, $S2c, $x1, $x1e; 
        printf CONTROLFILE "%0.4f  %0.4f  %0.4f  %0.4f  ", $c, $ce, $mB, $mBe; 
	printf CONTROLFILE "%0.4f  %0.4f  %0.4f  \n", $NDOF, $SNRMAX1, $FITPROB;
    
	$snid++;
    }

    close (MYFILE);
    close (CONTROLFILE);
    
}

sub grab_corr(@)
{
    ##### UTILITY: for a given central redshift bin value, 
    ##### grabs appropriate bias correction from binned
    ##### bias correction data array. For this to work, 
    ##### TEST data and BIAS CORRECTION data must be binned
    ##### identically in redshift. (JLM OCT 2013)

    my ($p_corrarray, $zval) = (@_);
    my (@subarray, $ntemp, $tempcorrZ);
    my ($tempcorrNEW, $tempcorrerrNEW, $ncorrD, $ncorrND, $ncorrALL);

    # chooses row of corrarray with appropriate redshift
    @subarray = grep { $_->[0] == $zval } @$p_corrarray;
    $ntemp = @subarray;
    
    if ( $ntemp == 1 ) {
	$tempcorrZ = $subarray[0][1];
	$tempcorrNEW = $subarray[0][2];
	$tempcorrerrNEW = $subarray[0][3];
	$ncorrD = $subarray[0][4];
	$ncorrND = $subarray[0][5];
	$ncorrALL = $subarray[0][6];
    } else {
	    $MSGERR[0] = "problem with grab_corr extraction routine";
	    $MSGERR[1] = "number of elements is $ntemp ";
	    $MSGERR[2] = "expected number to be 1 ";
	    $MSGERR[3] = "check code ";
	    sntools::FATAL_ERROR(@MSGERR);    	
    }

    return($tempcorrNEW, $tempcorrerrNEW, $ncorrD, $tempcorrZ);

} # sub grab_corr(@)

sub get_cosmomu(@)
{

    #### UTILITY to calculate cosmological distance modulus from 
    #### cosmo array (see build_cosmo_array) for a given
    #### redshift range.

    my ($p_cosmoarray, $zval) = (@_);

    my (@subarray, $nval, $ival, $minaddress, $tempmin, $min);
    my (@ztemp, $ntemp);
    my ($muval, $zlow, $zhigh);

    # determine bin redshift range
    ($zlow, $zhigh) = get_binzrange($zval);

    # pick out elements of cosmo array falling at or near
    # center of bin redshift. 
    @subarray = grep { $_->[0] == $zval } @$p_cosmoarray; #at bin redshift
    $ntemp = @subarray;
    
    if ( $ntemp == 1 ) {
	$muval = $subarray[0][1];
    } elsif ( $ntemp == 0 ) {
	# determine array values within bin redshift range
	@subarray = grep { $_->[0] >= $zlow && $_->[0] <= $zhigh } @$p_cosmoarray;
	@ztemp = map $_->[0], @subarray;
	$nval = @ztemp;
	$min = 10000.;
	$minaddress = -9;
	# select values immediately above and below desired bin redshift value
	for ( $ival = 0; $ival < $nval; $ival++ ) {
	    $tempmin = abs($ztemp[$ival]-$zval);
	    if ( $tempmin < $min && $ztemp[$ival] < $zval )  {
		$min = $tempmin;
		$minaddress = $ival;
	    }
	}

	# abort if can't find cosmo values within bin redshift range
	if ( $minaddress < 0 ) {
	    $MSGERR[0] = "problem with get_cosmomu extraction routine";
	    $MSGERR[1] = "invalid minaddress $minaddress ";
	    $MSGERR[2] = "check code ";
	    sntools::FATAL_ERROR(@MSGERR);    	
	} else {
	    # otherwise, use linear interpolation to get distance modulus at redshift val
	    $muval = lininterp($subarray[$minaddress][0], $subarray[$minaddress][1],
		   $subarray[$minaddress+1][0], $subarray[$minaddress+1][1],
		   $zval);   
	}
    }else {
	    $MSGERR[0] = "problem with get_cosmomu extraction routine";
	    $MSGERR[1] = "number of elements is $ntemp ";
	    $MSGERR[2] = "expected number to be 0 or 1 ";
	    $MSGERR[3] = "check code ";
	    sntools::FATAL_ERROR(@MSGERR);    	
    }
    
    return($muval);
    
} # sub get_cosmomu(@)

sub lininterp(@)
{
    #### UTILITY to perform
    #### basic linear interpolation between 
    #### two data points

    my ($x0, $y0, $x1, $y1, $x) = (@_);

    my ($percent);
    my ($val);

    $percent = ($x-$x0)/($x1-$x0);
    $val = $percent*($y1-$y0) + $y0;

    return($val);
    
} # end lininterp(@)

sub survey_setup(@)
{

    my ($infile, $p_survbins, $p_survvec, $p_outfiles, $p_datfiles) = (@_);
    
    my ($temp_ngrid, @temp_vec);
    my ($addstring, $iadd, $sval); 

    if ( $SURVFLAG == 0 ) {
	$temp_ngrid = 1;
	$$p_survvec[0] = -1;
	$$p_survbins[0] = $GBINSIZE; 
	$$p_outfiles[0] = $goutfile;
	$$p_datfiles[0] = $gdatfilebin;
    } else {
	#### GRAB IDSURVEY COLUMN FROM FITRES FILE
	$temp_ngrid = sntools::rdfitres_getcol($infile, \@temp_vec, "IDSURVEY", "VARNAMES:", "SN:", 0, $VERBOSE);
	
	#### GRAB UNIQUE ELEMENTS OF idsurvey
	sntools::uniq_entries(\@temp_vec, $p_survvec);
	##print "xxxDEBUG survey_setup @{$p_survvec} \n";

	#### MAKE NEW NAMES FOR OUTPUT FILES
	$iadd = 0;
	foreach $sval (@{$p_survvec}) {
	    $addstring = "S" . $sval . "S";
	    $$p_outfiles[$iadd] = add_to_filename($goutfile, $addstring);
	    $$p_datfiles[$iadd] = add_to_filename($gdatfilebin, $addstring);
	    $iadd++;
	}
	$temp_ngrid = unshift(@{$p_outfiles}, $goutfile);
	$temp_ngrid = unshift(@{$p_datfiles}, $gdatfilebin);

	@{$p_survbins} = (($SURVBINSIZE) x @{$p_survvec} ); #initialize surveybins vector of length @p_survvec
        #### ADD FULL DATA SET $BINSIZE to beginning of survey BINSIZE vector
	$temp_ngrid = unshift(@{$p_survbins}, $GBINSIZE); 

	#### ADD -1(=all entries) to beginning of UNIQUE survey vector
	$temp_ngrid = unshift(@{$p_survvec}, "-1");
	

    }

    return ($temp_ngrid);

} # end survey_setup(@)

sub get_zrange(@)
{
    my ($p_zarr) = (@_);
    my ($col, $narr, $iarr, @wdLine, $min, $max);
    my ($altwd);

    $max = -1.0;
    $min = 10.0;

    $narr = @$p_zarr; 
    $iarr = 0;
    while ($iarr < $narr) {
	$altwd = $$p_zarr[$iarr];
	if ( $altwd < $min ) {$min = $altwd;}
	if ( $altwd > $max ) {$max = $altwd;}
	$iarr++;
    }

    return ($min, $max);

} #end get_zrange(@)

sub check_P11binsize(@)
{
    #### UTILITY to determine and return the number of SNe in a given 
    #### redshift bin. For P11 bias corrections, checks all three 
    #### data arrays (ALL, DETECTED, and NON-DETECTED).

    my ($p_corr, $p_corrD, $p_corrND, $survid) = (@_);
    
    my ($BINSIZEflag, $iz, $NZ, $ALLflag, $Dflag);
    my ($ncorrALL, $ncorrD, $ncorrND);
    my ($SUBname);

    $SUBname = "check_P11binsize";

    #### CHECK ADEQUACY OF INPUT BINSIZE
    print "\nChecking P11 bias corr BINSIZE ... \n";
    print "\t Requiring nALL and nD to be at least Nlimit:$Nlimit SNe\n";
    print "\t in each bin.\n\n";
    $BINSIZEflag = -1;
    while ( $BINSIZEflag < 0 ) {
	$NZ = $ALLflag = $Dflag = 0;
        ## LOOP TO COMPARE nALL, nD IN EACH BIN WITH $Nlimit
	for ( $iz =$ZMIN; $iz <= $ZMAX ; $iz=$iz+$BINSIZE ) {
	    ($ncorrALL) = check_binsize($p_corr, $iz, $BINSIZE);
	    ($ncorrD) = check_binsize($p_corrD, $iz, $BINSIZE);
	    ($ncorrND) = check_binsize($p_corrND, $iz, $BINSIZE);
	    if ($ncorrALL < $Nlimit) {$ALLflag = $ALLflag + 1;}
	    if ($ncorrD < $Nlimit) {$Dflag = $Dflag + 1;}
	    $NZ++;
	}
	## CHECK flag results
	if ($ALLflag == 0 & $Dflag == 0 ) {
	    $BINSIZEflag = 1;
	    print "\t For BINSIZE = $BINSIZE, all bins have N >= $Nlimit.\n";
	    print "\t Continuing with bias correction ... \n";
	} elsif ( $BINSIZE > ($ZMAX-$ZMIN) ) {
	    $MSGERR[0] = "In subroutine $SUBname ";
	    $MSGERR[1] = "invalid BINSIZE  ";
	    $MSGERR[2] = "Check input fitres files";
	    sntools::FATAL_ERROR(@MSGERR);
	} else {
	    print "\t For BINSIZE = $BINSIZE, some bins fail $Nlimit cut.\n";
	    print "\t\t  Of ALL      SNe: $ALLflag of $NZ bins fail cut.\n";
	    print "\t\t  Of DETECTED SNe: $Dflag of $NZ bins fail cut.\n\n";
	    $BINSIZE = $BINSIZE + $BINSIZEstep;
	    #($ZMIN, $ZMAX) = check_z(); OBSOLETE?
	    print "\t Checking with new BINSIZE: $BINSIZE (step: $BINSIZEstep) ... \n";
	    #print "\t Checking with new ZMIN, ZMAX: $ZMIN, $ZMAX ... \n\n";
	} #end of flag check if statement
    } #end of while( $BINSIZEflag < 0 ) loop

} #end check_P11binsize(@)

sub check_z(@)
{

    #### UTILITY to determine largest possible redshift
    #### range based on TEST data zmin, zmax ($dzmin, $dzmax) 
    #### and BIAS CORR data zmin, zmax ($bzmin, $bzmax)
    #### UTILITY uses survey loop level binsize ($BINSIZE) variable

    my ($loc_dzmin, $loc_dzmax, $loc_bzmin, $loc_bzmax, $survid) = (@_);
    
    my($dzmin, $dzmax, $bzmin, $bzmax);
    my($zmin, $zmax);

    print "\n\t\tChecking redshift ranges for SURVEY ID $survid ... \n";
    #get DZMIN and DZMAX to nearest 0.01 
    $dzmin = floor($loc_dzmin*100)/100;
    $dzmax = ceil($loc_dzmax*100)/100;
    $bzmin = floor($loc_bzmin*100)/100;
    $bzmax = ceil($loc_bzmax*100)/100;
    
    #make sure BZMIN and BZMAX cover range
    if ( $dzmin < $bzmin ) {
	print "\t\tWARNING: survey bdata doesn't extend to data min $dzmin!!!! \n";
	print "\t\t         survey bdata min is $loc_bzmin, data min is $loc_dzmin \n";
	$zmin = $bzmin;
    } else {
	$zmin = $dzmin; 
    }

    if ( $dzmax > $bzmax ) {
	print "\t\tWARNING: survey bdata doesn't extend to data max $dzmax!!! \n";
	print "\t\t         survey bdata max is $loc_bzmax, data max is $loc_dzmax \n";
	$zmax = $bzmax;
    } else {
	$zmax = $dzmax;
    }

    # ALLOW zmax to cover full range of data zmax, even if bias corrections don't fully
    # overlap. This allows bin of width BINSIZE to exist even if top end isn't totally 
    # covered by bias correction data. 
    if ($survid == -1 ) {$zmax = $dzmax;}

    print "\t\t Running survey $survid corrections from ZMIN: $zmin to ZMAX: $zmax \n";
    $zmin = $zmin+($BINSIZE/2.); #zmin is the center of the bin...
    print "\t\t\t Center of first bin will be $zmin \n";

    return($zmin, $zmax);

} # end check_z(@)

sub getcol_failcheck(@)
{
    #### UTILITY to provide additional debugging info in case of 
    #### sntools::rdfitres_getcol failure

    my ($checkflag, $VARNAME, $SUBNAME, $verbose) = (@_);

    if ($checkflag == -1 ) {
	$MSGERR[0] = "In subroutine $SUBNAME ";
	$MSGERR[1] = "for column varname $VARNAME ";
	$MSGERR[2] = "rdfitres_getcol failed ";
	$MSGERR[3] = "Check that snana header names match expectations";
	sntools::FATAL_ERROR(@MSGERR);
    } else {
	if ( $verbose ) { print "\t successfully found $VARNAME data ...\n\n";}
    }

} # end getcol_failcheck(@)

sub make_P11_corrarray(@)
{
    my ($p_P11corrarray, $survid) = (@_); 

    ## utility variables
    my ($SUBname);
    my (@cosmo, @corr, @sortcorr, @sortcorrD, @sortcorrND);
    my ($ngrid, $iarr, $narr, $iz, $NZ);
    ## generic variables for rebin_data output
    my ($mSMU, $mSMUe, $mid);
    ## specific variables for rebin_data output
    my ($mZALL, $mZeALL, $mZND, $mZeND, $mZD, $mZeD);
    my ($zlowD, $zhighD, $zlowND, $zhighND, $zlowALL, $zhighALL);
    my ($murefALL, $murefD, $murefND);
    ## variables for correction calculation
    my ($tcorrALL, $tcorreALL, $ncorrALL);
    my ($tcorrD, $tcorreD, $ncorrD);
    my ($tcorrND, $tcorreND, $ncorrND);
    my ($coeff, $termALL, $termND, $termD);
    my ($tcorr, $tcorre);
    ## index variable for correction array
    my ($binz);
    ## survey related variables
    my ($survey_bzmin, $survey_bzmax);

    $SUBname = "make_P11_corrarray";
	
    #########
    #########
    #### P11 BIAS CORRECTION: 
    ####
    ####   CORRECTIONS: mu_{corr,TEST}(zbin) = mu_{TEST}(zbin) + Dmu_{FIX}(zbin)
    ####
    ####   FIX data: simulated with TEST alpha, beta
    ####             includes ALL SNe, detected ("D") and not detected ("ND")
    ####
    ####   SNANA EFFMASK VARIABLES USED TO SORT "D" and "ND":
    ####           D: effmask == 3
    ####          ND: effmask <= 2
    ####
    ####   CALCULATING CORRECTIONS Dmu_{FIX}(zbin):
    ####
    ####     All mu are first corrected to mu' by replacing fitted FIX alpha/beta
    ####     with simulated (i.e. fitted TEST) alpha/beta.
    ####
    ####            mu' = mu - (a_sim - a_fit)*x1_fit - (b_fit - b_sim)*c_fit
    ####
    ####     Then binned corrections Dmu(zbin) and uncertainties Dmue(zbin)
    ####     are calculated.
    ####
    ####     Dmu(zbin)  =  < mu' ALL >(zbin) - < mu' DETECTED >(zbin)
    ####                               -OR- 
    ####                =  C'[< mu' ND >(zbin) - < mu' D >(zbin)], where
    ####
    ####             C' = C_ND/C_ALL with C defined as the sum of weights.
    ####                  [i.e. C_ALL = sigma(i=1 to i=N_ALL) 1.0/(dmu_i^2)] 
    ####   
    ####     The second definition of Dmu(zbin) facilitates Dmue calculation.
    ####
    ####     Dmue(zbin) =  C'[(sig_muND)^2 + (sig_muD)^2]^(1/2), where
    ####
    ####             C' is defined as above and sig_<type> are
    ####             sig_muND = sqrt(1.0/C_ND) and sig_muD  = sqrt(1.0/C_D)
    ####             
    ####   TWO TYPES OF CORRECTIONS, "TRUE" and "APPX", ARE AVAILABLE:
    ####       
    ####       Both use Dmu(zbin) = < mu' ALL >(zbin) - < mu' DETECTED >(zbin) 
    ####
    ####       TRUE: Dmue(zbin) =  [C'*((sig_muND)^2 + (sig_muD)^2)]^(1/2)
    ####
    ####       APPX: Dmue(zbin) = sqrt((N_ALL-N_D)/N_ALL)*$tcorreD;
    ####
    #########
    #########

    print "\n\tGenerating P11 corrections ... \n";

    ### GRAB NEEDED DATA FROM BCORRFILE, MAKE SORTED COSMO (z, simmu) ARRAY
    build_cosmo_array($bcorrfile, \@cosmo);

    ### GRAB NEEDED DATA FROM BCORRFILE, DETERMINE P11 MU' VALUES, PACKAGE IN ARRAY
    ($survey_bzmin, $survey_bzmax) = build_fitres_array($bcorrfile, "P11", \@corr, $survid);

    ### REDETERMINE ZMIN, ZMAX VALUES FOR THIS SURVEY
    ($ZMIN, $ZMAX) = check_z($DZMIN, $DZMAX, $survey_bzmin, $survey_bzmax, $survid);

    #### SORT ARRAY IN ASCENDING ORDER OF Z, MAKE D, ND SUBARRAYS
    @sortcorr = sort { $a->[0] <=> $b->[0] } @corr;
    @sortcorrD = grep { $_->[5] == 3 } @sortcorr;
    @sortcorrND = grep {$_->[5] < 3 } @sortcorr;

    #### CHECK THAT SURVEY BINSIZE GIVES ADEQUATE NUMBERS OF N_D,N_ND PER BIN
    print "\n\tChecking number of detected SNe in FIX fitres file ...";
    $ngrid = @sortcorrD;
    if ( $ngrid == 0 ) {
	$MSGERR[0] = "In subroutine $SUBname ";
	$MSGERR[1] = "NO DETECTED SNE (EFF_MASK == 3)";
	$MSGERR[2] = "in FIX fitres file. ";
	$MSGERR[3] = "Increase NSNe in DMBIASCORR SIMS";
	sntools::FATAL_ERROR(@MSGERR);    	
    } 
    check_P11binsize(\@sortcorr, \@sortcorrD, \@sortcorrND);
    print "\t\tUSING BINSIZE $BINSIZE \n\n";

    ########
    ########  CALCULATE P11 BIAS CORR Dmu(zbin), Dmue(zbin)
    ########

    $NZ = 0;
    for ( $iz =$ZMIN; $iz <= $ZMAX ; $iz=$iz+$BINSIZE ) {
	$binz = $iz; 
	if ( $VERBOSE ) { print "\t\tRebinning ALL FIX data ... \n"; }
	($mZALL, $mZeALL, $tcorrALL, $tcorreALL, $mSMU, $mSMUe, 
	         $mid, $ncorrALL, $zlowALL, $zhighALL, $murefALL) 
	         = rebin_data(\@sortcorr, \@cosmo, $iz, $BINSIZE, "P11");
	if ( $VERBOSE ) { print "\t\tRebinning ALL DETECTED FIX data ... \n"; }
	($mZD, $mZeD, $tcorrD, $tcorreD, $mSMU, $mSMUe, 
                 $mid, $ncorrD, $zlowD, $zhighD, $murefD) 
	         = rebin_data(\@sortcorrD, \@cosmo, $iz, $BINSIZE, "P11");
	if ( $VERBOSE ) { 
	    printf "\t\txxxP11 BCORRINFO <meanz> <meanzerr> <mu_ND> <mue_ND> <nND> <zlow> <zhigh> <muref>\n";
	    printf "\t\txxxP11 BCORRDETECTED %0.4f %0.4g %0.4f %0.4f %d %0.4f %0.4f %0.4f \n", 
	      $mZALL, $mZeALL, $tcorrALL, $tcorreALL, $ncorrALL, $zlowALL, $zhighALL, $murefALL;
	    printf "\t\txxxP11 BCORRDETECTED %0.4f %0.4g %0.4f %0.4f %d %0.4f %0.4f %0.4f \n", 
	      $mZD, $mZeD, $tcorrD, $tcorreD, $ncorrD, $zlowD, $zhighD, $murefD;
	}

	# CORRECTION Dmu(zbin)
	$tcorr = $tcorrALL-$tcorrD;

	# ERROR CALCULATION Dmue(zbin)
	if ( $P11TYPE eq "TRUE" ) {
	    $termALL = $tcorreALL*$tcorreALL;
	    $termD = $tcorreD*$tcorreD;
	    #TRUE REQUIRES ND FOR ERROR CALCULATION
	    ($ncorrND) = check_binsize(\@sortcorrND, $iz, $BINSIZE);
	     if ( $ncorrND > 2 ) {
		 if ( $VERBOSE ) { print "\t\tRebinning ALL UNDETECTED FIX data ... \n"; }
		 ($mZND, $mZeND, $tcorrND, $tcorreND, $mSMU, $mSMUe, $mid, $ncorrND, $zlowND, $zhighND, $murefND) 
		       = rebin_data(\@sortcorrND, \@cosmo, $iz, $BINSIZE, "P11");
		  if ( $VERBOSE ) {printf "\t\txxxP11 BCORRNOTDETECTED %0.4f %0.4g %0.4f %0.4f %d %0.4f %0.4f %0.4f \n", 
				   $mZND, $mZeND, $tcorrND, $tcorreND, $ncorrND, $zlowND, $zhighND, $murefND; }
		 $termND = $tcorreND*$tcorreND;
		 $coeff = $termALL/$termND;
		 $tcorre = sqrt($coeff*($termND + $termD));
	     } else {
		 if ( $VERBOSE ) {printf "\t\txxxP11 BCORRNOTDETECTED %0.4f %0.4g %0.4f %0.4f %d %0.4f %0.4f %0.4f \n", 
				  0.0, 0.0, 0.0, 0.0, $ncorrND, $zlowND, $zhighND, 0.0; }
		 $tcorrND = 0.0;
		 $tcorre = 0.0;
	     } 
        } elsif ( $P11TYPE eq "APPX" ) {
	    #ONLY D, ALL USED FOR ERROR CALCULATION
	    $tcorre = sqrt(($ncorrALL-$ncorrD)/$ncorrALL)*$tcorreD;
	} else {
	    $MSGERR[0] = "INVALID P11TYPE. ";
	    $MSGERR[1] = "Valid options are TRUE, APPX ";
	    $MSGERR[2] = "Check inputs and try again ";
	    sntools::FATAL_ERROR(@MSGERR);
	}
	if ( $VERBOSE ) {
	    printf "\t\txxxP11 BCORRFINAL (binz, meanz, corr, corre)  %0.4f %0.4f %0.4f %0.4f \n\n", $binz, $mZD, $tcorr, $tcorre;
	}
	$$p_P11corrarray[$NZ] = [$binz, $mZD, $tcorr, $tcorre, $ncorrD, $ncorrND, $ncorrALL];
	$NZ++;
    }
    
} #end make_P11_corrarray(@)

sub build_fitres_array(@)
{

    #### UTILITY to extract, package needed columns from 
    #### input fitres file $fitresfile. 
    #### 
    #### This subroutine is set up to filter data based on 
    #### survey and to screen for zeros in mBERR values 
    ####
    #### JLM, Oct 08 2013 -- certain pathological SNe fits pass 
    #### fit cuts and end up in fitres file. These pathological
    #### fits seem to be more likely with C11 color smearing, 
    #### and feature (among other things) maxed out PEAKMJD 
    #### errors and mBERR == 0. Allowing these SNe into 
    #### bias correction calculations crashes weighted average
    #### routines used in binning. 
    #### 
    
    my ($fitresfile, $type, $p_fitresarray, $survid) = (@_);

    my ($temp_ngrid, $SUBname);
    my (@cid_dat, @z_dat, @zerr_dat, @mu_dat, @muerr_dat);
    my (@id_dat, @simmu_dat, @effmask_dat, @x1_dat, @c_dat); 
    my (@S2x1_dat, @S2c_dat, @NDOF_dat, @SNRMAX1_dat, @mB_dat);
    my (@FITPROB_dat, @x1err_dat, @cerr_dat, @mBerr_dat);
    my (@checkmB);
    my ($fita, $fitb, $da, $db);
    my ($narr, $iarr, $ioutput);
    my ($mu1, $mu2, $mu3);
    my ($temp_outvec, @z_survey);
    my ($survey_zmin, $survey_zmax);

    print "\n\tGathering raw data from $fitresfile \n";

    ########
    ######## SUB build_fitres_array: Collects and bins FITRES SNe data
    ########
    ######## output array [FOR BIAS CORR][FOR TEST INFO]:
    ########   z, ze, mu1, mue, mu2, type[D/ND], survey id, cid, mu3,  
    ########          S2x1, S2c, NDOF, SNRMAX1, x1, x1e, c, ce, 
    ########          mB, mBe, FITPROB                  
    ########
    ########  where mu1: 
    ########
    ########      $mu_dat + ($da*$x1_dat) + ($db*$c_dat) (P11)
    ########
    ########      $simmu_dat - $mu_dat (R13)
    ########
    ########      $mu_dat (all other types)
    ########
    ########  where mu2: 
    ########
    ########      $mu_dat (P11, R13)
    ########
    ########      $simmu_dat (all other types)
    ########
    ########  where mu3:
    ########      
    ########      $simmu_dat (all types)
    ########

    #### GRAB DATA FROM FITRES FILE: Z, ZERR, MU, MUERR, IDSURVEY, SIM_MU, EFFMASK, CID, x1, c 
    ####                                x1ERR, cERR, S2x1, S2c, NDOF, SNRMAX1, mB, FITPROB                  
    $SUBname = "build_fitres_array";
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@z_dat, "Z", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "Z", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@zerr_dat, "ZERR", "VARNAMES:", "SN:", 1, $VERBOSE);
     getcol_failcheck($temp_ngrid, "ZERR", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@mu_dat, "MU", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "MU", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@muerr_dat, "MUERR", "VARNAMES:", "SN:", 1, $VERBOSE);
     getcol_failcheck($temp_ngrid, "MUERR", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@id_dat, "IDSURVEY", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "IDSURVEY", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@simmu_dat, "SIM_MU", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "SIM_MU", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@effmask_dat, "EFFMASK", "VARNAMES:", "SN:", 0, $VERBOSE);
     if ($temp_ngrid == -1) {	
	$temp_ngrid = 
	    sntools::rdfitres_getcol($fitresfile, \@effmask_dat, "SIM_EFFMASK", "VARNAMES:", "SN:", 0, $VERBOSE);
	getcol_failcheck($temp_ngrid, "SIM_EFFMASK and EFFMASK", $SUBname, $VERBOSE);
     } else {
	getcol_failcheck($temp_ngrid, "EFFMASK", $SUBname, $VERBOSE);
     }
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@cid_dat, "CID", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "CID", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@x1_dat, "x1", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "x1", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@c_dat, "c", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "c", $SUBname, $VERBOSE);
    ### VARIABLES FOR TEST DATA
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@x1err_dat, "x1ERR", "VARNAMES:", "SN:", 1, $VERBOSE);
     getcol_failcheck($temp_ngrid, "x1ERR", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@cerr_dat, "cERR", "VARNAMES:", "SN:", 1, $VERBOSE);
     getcol_failcheck($temp_ngrid, "cERR", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@mBerr_dat, "mBERR", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "mBERR", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@S2x1_dat, "SIMx1", "VARNAMES:", "SN:", 0, $VERBOSE);    
    if ( $temp_ngrid == -1 ) {
	    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@S2x1_dat, "S2x1", "VARNAMES:", "SN:", 0, $VERBOSE);    
	    if ( $temp_ngrid == -1 ) {
		$temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@S2x1_dat, "SIMSED_S2x1", "VARNAMES:", "SN:", 0, $VERBOSE);    
	    }
	    getcol_failcheck($temp_ngrid, "SIMx1, S2x1, and SIMSED_S2x1", $SUBname, $VERBOSE);
    } 
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@S2c_dat, "SIMc", "VARNAMES:", "SN:", 0, $VERBOSE);    
    if ( $temp_ngrid == -1 ) {
	    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@S2c_dat, "S2c", "VARNAMES:", "SN:", 0, $VERBOSE);    
	    if ( $temp_ngrid == -1 ) {
		$temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@S2x1_dat, "SIMSED_S2c", "VARNAMES:", "SN:", 0, $VERBOSE);    
	    }
	    getcol_failcheck($temp_ngrid, "SIMc, S2c, and SIMSED_S2c", $SUBname, $VERBOSE);
    }
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@NDOF_dat, "NDOF", "VARNAMES:", "SN:", 0, $VERBOSE);    
     getcol_failcheck($temp_ngrid, "NDOF", $SUBname, $VERBOSE);
    $temp_ngrid = 
	sntools::rdfitres_getcol($fitresfile, \@SNRMAX1_dat, "SNRMAX1", "VARNAMES:", "SN:", 0, $VERBOSE);    
     getcol_failcheck($temp_ngrid, "SNRMAX1", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@mB_dat, "mB", "VARNAMES:", "SN:", 0, $VERBOSE);    
     getcol_failcheck($temp_ngrid, "mB", $SUBname, $VERBOSE);
    $temp_ngrid = 
	sntools::rdfitres_getcol($fitresfile, \@FITPROB_dat, "FITPROB", "VARNAMES:", "SN:", 0, $VERBOSE);    
     getcol_failcheck($temp_ngrid, "FIPROB", $SUBname, $VERBOSE);


    #### GRAB FITTED ALPHA AND BETA FROM TEST SN INPUT FILE
    #### ONLY USED FOR P11 SIMCORR
    ####  (AUG 30 2013: FIXED BUG USING TEST FILE FOR BOTH SIM, fit VALUES)
     $fita = qx(grep alpha0 $fitresfile | awk '{print \$4}');
     $fitb = qx(grep beta0 $fitresfile | awk '{print \$4}');
     $da = -1.0*($SIMA-$fita); ##check OCT7 2013 JLM
     $db = -1.0*($fitb-$SIMB); ##check OCT7 2013 JLM

    print "\t\tUSING sima $SIMA, fita $fita \n";
    print "\t\tUSING simb $SIMB, fitb $fitb \n";

    #### MAKE ARRAY OF NEEDED VARIABLES: z, ze, mu1, mue, mu2, type[D/ND], survey id, cid  
    $narr = @z_dat;
    $ioutput = 0;
    for ($iarr = 0; $iarr<$narr; $iarr++) {
	if ( $type eq "P11" ) {
	    # use simulated alpha, beta rather than fitted alpha beta
	    # $mu1 = $mu_dat[$iarr]+ ($da*$x1_dat[$iarr]) + ($db*$c_dat[$iarr]);
	    $mu1 = $mu_dat[$iarr] - ($da*$x1_dat[$iarr]) - ($db*$c_dat[$iarr]);
	    $mu2 = $mu_dat[$iarr];
	} elsif ( $type eq "R13" ) {
	    $mu1 = $simmu_dat[$iarr] - $mu_dat[$iarr];
	    $mu2 = $mu_dat[$iarr];
	} else {
	    $mu1 = $mu_dat[$iarr];
	    $mu2 = $simmu_dat[$iarr];
	}
	$mu3 = $simmu_dat[$iarr];
	$temp_outvec =  [$z_dat[$iarr], $zerr_dat[$iarr], $mu1, $muerr_dat[$iarr], $mu2, 
			$effmask_dat[$iarr], $id_dat[$iarr], $cid_dat[$iarr], $mu3, 
	                $S2x1_dat[$iarr], $S2c_dat[$iarr], $NDOF_dat[$iarr], $SNRMAX1_dat[$iarr], 
			$x1_dat[$iarr], $x1err_dat[$iarr], $c_dat[$iarr], $cerr_dat[$iarr], 
			$mB_dat[$iarr], $mBerr_dat[$iarr], $FITPROB_dat[$iarr]];

	# ANALYZE FULL SET OF DATA ( $survid == -1 )
	if ( $survid < 0 && $mBerr_dat[$iarr] > 0 ) {
	    $$p_fitresarray[$ioutput] = $temp_outvec;
	    push(@z_survey, $z_dat[$iarr]);
	    push(@checkmB, $mBerr_dat[$iarr]);
	    $ioutput++;
	} else {
	# ANALYZE SPECIFIC SURVEY $survid DATA
        # REQUIRE mBerr > 0! 
	    if ( $id_dat[$iarr] == $survid && $mBerr_dat[$iarr] > 0 ) {
		push(@z_survey, $z_dat[$iarr]);
         	push(@checkmB, $mBerr_dat[$iarr]);
		$$p_fitresarray[$ioutput] = $temp_outvec;
		$ioutput++;
	    }
	}
    } ## LOOP OVER ALL SN in FITRES FILE -> BUILDING OUTPUT FITRES ARRAY
    
    ($survey_zmin, $survey_zmax) = get_zrange(\@checkmB); #check output mBerr
    if ( $VERBOSE ) {
	print "\t\t\t xxxbuild_fitres : survey min, max mBerr are $survey_zmin, $survey_zmax \n";
    }

    ($survey_zmin, $survey_zmax) = get_zrange(\@z_survey); #check new survey z limits
    if ( $VERBOSE ) {
	print "\t\t\t xxxbuild_fitres : survey min, max are $survey_zmin, $survey_zmax \n";
    }
    
    return($survey_zmin, $survey_zmax);

} #end build_fitres_array(@)

sub build_cosmo_array(@)
{

    my ($fitresfile, $p_cosmoarray) = (@_);

    my ($temp_ngrid, $SUBname, @temp_array, @avg_array, @uniq_entries);
    my (@z_cosmo, @zerr_cosmo, @simmu_cosmo, @temp);
    my (@z_unique, $uniq_zval);
    my (@sub_z, @sub_zerr, @sub_simmu, $nsub);
    my ($zavg, $zsig, $simmuavg, $tempavg, $NBIN);
    my ($iarr, $narr);

    print "\n\t\tMaking cosmo array from $fitresfile \n"; 

    ########
    ######## Collect and bin FITRES SNe data
    ########

    #### GRAB DATA FROM FITRES FILE: SIMZ, SIM_MU
    $SUBname = "build_cosmo_array";
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@z_cosmo, "SIMZ", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "SIMZ", $SUBname, $VERBOSE);
    $temp_ngrid = sntools::rdfitres_getcol($fitresfile, \@simmu_cosmo, "SIM_MU", "VARNAMES:", "SN:", 0, $VERBOSE);
     getcol_failcheck($temp_ngrid, "SIM_MU", $SUBname, $VERBOSE);

    #### GRAB UNIQUE ELEMENTS OF z
    sntools::uniq_entries(\@z_cosmo, \@z_unique);

    #### MAKE ARRAY OF NEEDED VARIABLES: z, zerr, simmu
    $narr = $nsub = @z_cosmo;
    for ($iarr = 0; $iarr < $narr; $iarr++) {
	$temp_array[$iarr] = [$z_cosmo[$iarr], $simmu_cosmo[$iarr]];
    }

    #### MAKE ARRAY OF AVERAGED QUANTITIES zavg, simmuavg
    $narr = @z_unique;
    for ($iarr = 0; $iarr < $narr; $iarr++) {
	$uniq_zval = $z_unique[$iarr]; #get uniq_zval from array
	@uniq_entries = grep { $_->[0] == $uniq_zval } @temp_array; #pick out entries with z==uniq_zval
	@sub_z = map $_->[ 0 ],  @uniq_entries; #pick out columns
	@sub_simmu = map $_->[ 1 ],  @uniq_entries;	    
	$avg_array[$iarr] = [$sub_z[0], $sub_simmu[0]];
    }

    #### SORT ARRAY BY REDSHIFT, LOWEST TO HIGHEST
    @{ $p_cosmoarray } = sort { $a->[0] <=> $b->[0] } @avg_array;

    print "\t\tDone building cosmo array \n";

} #end build_cosmo_array(@)

sub get_binzrange(@)
{
    ### UTILITY CALCULATING bin zmin, zmax
    ### for a given central $zval. 
    ### UTILITY makes use of global (survey loop level)
    ### $BINSIZE variable. 

    my ($zval) = (@_);

    my ($zlow, $zhigh);

    $zlow = $zval - $BINSIZE/2.0;
    $zhigh = $zval + $BINSIZE/2.0;

    return ($zlow, $zhigh);

} #end get_binzrange(@)

sub add_to_filename(@)
{
    ### UTILITY ALLOWING SURVEY or DATATYPE specific
    ### suffix "$addstring" to be inserted into base 
    ### "$infilename" filename

    my ($infilename, $addstring) = (@_);
    
    my (@infiletoken, $nfiletok, $filebin);
    my ($ifiletok, $outfilename);

    @infiletoken = split(/\./, $infilename);
    $nfiletok = @infiletoken; 
    $filebin = $infiletoken[0];
    for ($ifiletok=1; $ifiletok<($nfiletok-1); $ifiletok++) {
	$filebin = join(".", $filebin , $infiletoken[$ifiletok]);
    }

    $outfilename = join(".", $filebin, $addstring, $infiletoken[-1]);

    return($outfilename);

} #end add_to_filename
