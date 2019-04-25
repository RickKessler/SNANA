#!/usr/bin/env perl
#
# R. Biswas, Nov 2010
# used to convert SALT2 format light curves to SNANA format. The usage can
#be see by running the code with no arguments.
#
# R. Biswas, Feb 2011
# Changing the use of filters printed out. From the filters in the file, 
# it is now a specified global filter for the entire version
#
# Changed to correctly handle empty lines in lightcurves.
# R. Biswas, ,Sat Jan 25 15:13:05 CST 2014
# Changed the subroutine readlcdata so that the lineno variable is only 
# incremented if the line is non-empty. This is used to set NOBS 
# R. Biswas Tue Jan 28 11:50:55 CST 2014
#
# Apr 18 2014 RK 
#   - nicer format to README file (easier to read)
#   - do NOT require DayMax (of no DayMax, PEAKMJD -> 0)
#
# =======================================
use strict;
use warnings;
use POSIX;
use Cwd;
use File::Basename;
use File::Path;
use Text::Wrap;
use subs qw(locatefile calcmag calcmagerr sortfilter calc_snanaflux printusage qsalt2file recordcandidate getargs setdirs checkinitialization getlcheader recordnonsaltfiles getcomments parselightcurveline readlcdata formatheader printsnanaformatfile getsnid);

sub printusage{
	print "$0 converts SALT2 format light curve files  (version 2.3.5 onwards)\n";
	print "which have a single datafile, rather than multiple datafiles  \n";	
	print "as in the previous versions \n"; 
	print "\n";

	my $numargs= $#ARGV +1;
	print "You have $numargs arguments \n" ;
	print "usage of $0 \n";
	print "perl $0 -survey <SURVEY> -version <VERSION>  -filters <FILTERLIST> \
[-filelist <path to filelist.txt>] [-datadir <DATADIR>] \
[-make] [-outdir OUTDIR] [-filearray FILEARRAY ]\n";
	print "square brackets denote optional arguments, in which case \n";
	print "datadir = current directory \n";
	print "outdir = datadir/SNANA_format \n";
	print "Eg. usage : \n";
	print "perl $0  -survey SNLS -version SNLS3year -make >& ../out \n";
	#if (@filearray) {print "filearray given \n";}

	print "$0 is meant to work with the newer SALT2format files, \n";
	print "but will not work with older formats which have multiple \n";
	print "files for each light curve \n";

	print "\n\n $0 attempts to distinguish files that are not in the SALT2 \n lightcurve format in the data directory (eg. it can distinguish the SALT2 spectra files) \n and parse only the light curve files. \n This is done using some simple tests \n and files which do not look like salt2 light curve files \n are listed in .NOTSALTFILE along with a reason. It is still easy for non salt files to pass this test, in which case the parsing operations of $0 would fail. Users should be careful not to leave extra files in the datadirectory which may fool the code.\n" 
}

	###global variables
###set globalsurveyname
my $cmdfilearray =0;
my $globalsurvey='';
my $globalversion='';
my $specifiedglobalfilters='';
my $prefix = '';
my $datadir='';#='/users/astro/rbiswas/doc/Salt2Training/docs/SNLS_Data/salt2-SNLS3-sample-1b';
my $mkdir='';
my $outdir='';
my $filelist='';
my $surveyset='';
my $versionset='';
my $fname='' ; 
my $dry='';
my $header_filt;
my @filearray='';
my @filenames;
	#Auxilliary files
my $ignorefile ;
my $listfile ;
my $readmefile ;
my $notsaltfilelist ;
################################################################################################
	###Main
	###Start
#print "start by getting args\n";
&getargs(@ARGV);
&setdirs();
&checkinitialization();

open LF,">",$listfile or die "Cannot find lc file:$listfile $_";
#print LF '';
open NS,">",$notsaltfilelist or die "Cannot find lc file:$notsaltfilelist $_";
open IL,">",$ignorefile or die "Cannot find lc file:$ignorefile $_";
print IL '';
my $locatecalls=0;

if ($cmdfilearray){
	my $numfilearray = @filearray;
	print "$0 obtained array of $numfilearray candidate files from the commandline array \n";
	for my $file( @filearray){
		#print $file,"\n" ;
		push (@filenames,$file);
		}
	}
else	{
	@filenames = &locatefile($datadir,""); ###list of light curve files in dir 
	print "$0 obtained array of filenames by reading directory \n";
	}
for my $file (@filenames){
	#print "File: $file \n";
	}
#exit ();
my $globalfilter='';
my @allfilters=('U','B','V','g','r','i','z','R','I');
print "Set filter ordering \n";
my $lightcurve;
my $lcfield='NULL';
my $header;
##############################LC loops###########################
		####loop over the light curve files
my $lcindex=1;
my $comments='';
print "comments initialized \n";

my $numfile =0;

print "Entering loop to deal with every file \n";
for my $file  (@filenames){ 
	next if ($file =~ /^\.\.?$/);
	next if ($file =~ /^\.?$/);
	next if ($file eq $0 );
	if ($file =~/.*\.pl/) {print "getting perl scripts \n";}
	$numfile++ ;
	#next if ('_'.$file=~/\_\./);
	print "numfile  = $numfile, file = $file, lcindex = $lcindex \n";
	&recordcandidate($numfile, $file );
	$file=$datadir.'/'.$file;
		###For each light curve reset those variables that will change
	$lightcurve='';
	my $magref;
	$header='';
		###variable $comments valid across all files
		### reset for individual light curve unless one lc is done
	#if ($lcindex<2) {$comments='';}

		###read the light curve data file
	print "starting file loops with", $file,"\n";
	open LC,$file or die "Cannot find lc file:$file $_";
	my @lcdata =<LC>;
	close LC;

		###declare header variables for the lc	
	my ($survey,$snid,$iauc,$photometry_version,$sntype);
	my ($ra,$dec,$magtyper,$magsys,$mwebv,$peakmjd,$field);
	my ($redshift_helio,$redshift_cmb,$redshift_final,$redshift_spec,$redshift_status);
	my  ($foundsn,$foundZ, $founddaymax);
	$header_filt='';
	
		###go through each line and set the info
		###reset the lineno at the beginning of the lc
	my $lineno=0;
	my $nobs;
	#my $commentsover='';
	my $issalt2file=1;

	my @listvars=($survey,$snid,$iauc,$photometry_version,$sntype,$ra,$dec,$magtyper,$magsys,$mwebv,$peakmjd,$field,$redshift_helio,$redshift_final,$redshift_spec,$redshift_status);
		
	my @headerinfo = &getlcheader(@lcdata);
	my $redshift_err='';
	($snid, $ra, $dec , $mwebv, $redshift_helio, $redshift_cmb, $peakmjd,$field,$foundZ,$foundsn,$founddaymax,$redshift_err) = @headerinfo;
	#if (!$redshift_err){$redshift_err='NULL';}
	print "header read off \n";
	print $foundsn;
	print $foundZ;
	print $founddaymax;
	print $foundsn, $foundZ, $founddaymax, "expect correct vals \n";
		##Check if it is salt2file, accordingly set $issalt2file
	print "will leave filename  $file if not legit \n";
	$issalt2file=&qsalt2file($foundsn, $foundZ, $founddaymax); 
	#print "returned, issalt2file = $issalt2file \n";
	&recordnonsaltfiles ($issalt2file,$foundsn,$foundZ,$founddaymax, $file);
	next if (0==$issalt2file);	
	if (!$field){
		$field ='NULL'; #Some valid lc don't have fields recorded
		}
#	if (!$xfoc){
#		$xfoc = 'NULL';
#	}
	
	print "$snid, $ra, $dec, $mwebv, $redshift_helio, $peakmjd,last =$field \n";
	$lcfield = $field;
	#print "Got lcfield \n";
	#print "issalt2file line 212", $issalt2file, "\n" ;
	if (1==$lcindex) {&getcomments(@lcdata);}

	print "got comments  \n";

	if (1==$issalt2file) {
		($nobs,$lightcurve) = &readlcdata(@lcdata);
		}
	$header = &formatheader($nobs ,@headerinfo);
	if (not($dry)){
		print $file ;
		print "CALL LIGHT CURVE BEGIN \n";	
		print $file ,"\n";
		print $lightcurve ;
		print "CALL LIGHT CURVE END \n";	
		&printsnanaformatfile ($snid,$header,$lightcurve);
	}
	$lcindex++;
}
sub getargs{
	my @input =@_;
		###Exit if no arguments;
	if (0==@input) {
		&printusage;
		exit();
		}
		###CHECK INPUT (comment out when done)
	#####for my $argind (0..$#input){
	#####	print "ARGIND: $argind $input[$argind] \n" 
	#####	}
	#####exit();
		###set global variables from command line arguments
		###remember a later (on command line) variable overrides 
		###earlier var, (so in perl -make -dry  dry overrides make)
	for my $argind (0..$#input){
		if  ($input[$argind]=~/-survey/){
			$globalsurvey = $input[$argind+1];
			$surveyset='True';
			print "set Survey = $globalsurvey \n";
			}		
		if  ($input[$argind]=~/-version/){
			$globalversion = $input[$argind+1];
			$prefix = $globalversion;
			$versionset='TRUE';
			print "set Version = $globalversion \n";
			}		
		if  ($input[$argind]=~/-filelist/){
			$filelist = $input[$argind+1];
			}		
		if  ($input[$argind]=~/-datadir/){
			$datadir = $input[$argind+1];
			print "command line argument datadir set to $datadir \n";
			}		
		if  ($input[$argind]=~/-outdir/){
			$outdir = $input[$argind+1];
			}		
		if  ($input[$argind]=~/-make/){
			$mkdir ='TRUE';
			print "made mkdir true \n";
			}		
		if  ($input[$argind]=~/-filearray/){
			print "getargs getting array of possible SN files \n";
			for my $ind ($argind+1..$#input){ 
				push (@filearray,$input[$ind]);
				}
			my $numfilearray =@filearray ;
			print "There are $numfilearray such files found by getarg \n";				
			$cmdfilearray = 1;
			}		
		if  ($input[$argind]=~/-filename/){
			$fname = $input[$argind+1];
			}		
		if  ($input[$argind]=~/-filters/){
			$specifiedglobalfilters = $input[$argind+1];
			print "+++++++++\n" ;
			print $specifiedglobalfilters , "\n";
			}		
		if  ($input[$argind]=~/-dry/){
			$dry = 'true';
			$mkdir ='';
			print "now dry = $dry and mkdir = $mkdir \n";
			}		
		}
	print "mkdir = $mkdir \n";
	print "got command line arguments \n"
}
sub setdirs{
	if (($versionset) && ($surveyset)){
		;}
	else{
		&printusage();
		print "\n\n\n";
		die "Cannot proceed without survey name or version";
	}
	my $datalength = length $datadir;
	print "the datadir is set to $datadir with len $datalength	\n";
	if (0 == length $datadir){
		$datadir=getcwd;
		print "using datadir = current working directory $datadir \n";
		}
	else {
		print "datadir taken from command line to be $datadir \n";
		}
	
	
	if (0 == length $outdir){
		$outdir=$datadir.'/SNANA_format';
		print "using outdir = $outdir \n";
	}
	
	if (-d $outdir) {
		print "The output directory exists \n";
		}
	elsif ($mkdir){
		print "getting mkdir = $mkdir \n";
		mkpath $outdir or die "output directory $outdir does not exist and cannot create it: $_";
		print "created $outdir \n" ;
		}
	else {
		print "Output directory does not exist and you have to create $outdir manually \n";
	die "Cannot proceed without output directory \n";
	}
	$ignorefile = $outdir.'/'.$globalversion.'.IGNORE';
	$listfile = $outdir.'/'.$globalversion.'.LIST';
	$readmefile = $outdir.'/'.$globalversion.'.README';
	$notsaltfilelist =$outdir.'/'.$globalversion.'.NOTSALTFILE';
}
sub checkinitialization(){
	print "datadir = $datadir \n";
	print "The output directory is $outdir \n";
	
		###Name auxilliary files
	
	print "Set auxilliary file names \n";
	print "ignorefile =$ignorefile \n";
	print "listfile = $listfile \n";
	print "readmefile = $readmefile \n";
	print "notsaltfilelist = $notsaltfilelist \n";
}

sub getlcheader{
	my @lc = @_;
	my ($survey,$snid,$iauc,$photometry_version,$sntype);
	my ($ra,$dec,$magtyper,$magsys,$mwebv,$peakmjd,$field);
	my ($redshift_helio,$redshift_final,$redshift_spec,$redshift_status,$redshift_err);
		###Filter parameters to decide what is a salt2file
	my $foundsn=0;
	my $founddaymax= 1 ; # 0;
	my $foundZ=0;
	my $redshift_cmb ='';
	$redshift_err ='';
	$peakmjd = 0.0 ; # RK Apr 16 2014

	for my $line (@lc){
			###get header info
		#print $line ;
		#print "entered loop of LC \n";
		if ($line=~'@SURVEY') {
			$survey=(split(/\s/,$line))[1];
			print $survey,"\n";
			}
		if ($line=~'@Z_HELIO') {#use to set salt2file filter
			$redshift_helio=(split(/\s/,$line))[1];
			print "redshift obtained = $redshift_helio \n";
			$foundZ =1;
			}
		if ($line=~'@Z_CMB') {#use to set salt2file filter
			$redshift_cmb=(split(/\s/,$line))[1];
			}
		if ($line=~'@SN ') {#use to set salt2file filter
			my $snidtmp=(split(/\s/,$line))[1];
			$snid = &getsnid($snidtmp);
			print "obtained SNID = $snid \n";
			$foundsn =1;
			}
		if ($line=~'@MWEBV') {
			$mwebv=(split(/\s/,$line))[1];
			print "mwebv is $mwebv \n";
			}
		if ($line=~'@DayMax') {#use to set salt2file filter
			$peakmjd=(split(/\s/,$line))[1];
			print "pkmjd is $peakmjd \n";
			$founddaymax=1;
			print "daymax found  $founddaymax \n";
			}
		if ($line=~'@RA') {
			$ra=(split(/\s/,$line))[1];
			print  "RA =$ra \n";
			}
#NEW ADDITION, RB, 
#Fri Dec  7 18:00:47 CST 2012
#		my $xfoc;
#		if ($line=~'@X_FOCAL_PLANE '){
#			$xfoc = (split(/\s/,$line))[1];
#			print "X_FOCAL_PLANE = $xfoc \n" ;
#			}
			
#END NEW ADDITION, RB
		if ($line=~'@z_source') {
			$redshift_status=(split(/\s/,$line))[1];
			print $redshift_status,"\n";
			}
		if ($line=~'@Redshift_err') {
			$redshift_err=(split(/\s/,$line))[1];
			print $redshift_err,"\n";
			}
		if ($line=~'@DEC') {
			$dec=(split(/\s/,$line))[1];
			print $dec,"\n";
			}
		if ($line=~'@FIELD') {
			$field=(split(/\s/,$line))[1];
			print "field found = $field \n";
			}
#Change now
	}
	#set variable to signify that header is read
		#if ($line=~'#end') {
		#	#$lineno++;
		#	}
	my @headerinfo;
	if (!$redshift_err){$redshift_err='NULL';}
	if ($redshift_cmb){
		@headerinfo = ($snid, $ra, $dec , $mwebv, $redshift_helio, $redshift_cmb,$peakmjd,$field, $foundZ,$foundsn,$founddaymax,$redshift_err);
		}
	else {
		@headerinfo = ($snid, $ra, $dec , $mwebv, $redshift_helio, 'NULL', $peakmjd,$field,$foundZ,$foundsn,$founddaymax,$redshift_err);
	}
	print $foundsn , $founddaymax, $foundZ, "saltq\n" ;
		return @headerinfo ;
}

sub recordnonsaltfiles{
	my ($issalt2file,$foundsnfile ,$foundZfile,$founddaymaxfile, $file)  = @_;
	if (0==$issalt2file){
		print NS $file ,"\n"; 
		print NS $foundsnfile, $foundZfile, $founddaymaxfile ,"\t";
		if (0==$foundsnfile){
			print  NS "keyword  SN not found \t"; 
			}
		if (0==$foundZfile){
			print  NS "keyword  Z not found \t";
			}
		if (0==$founddaymaxfile){
			print  NS "keyword  Daymax not found \t"; 
			}
		print NS "\n";
		}

	#elsif (1==$issalt2file) {
	#		print LF $file ,"\n";
	#	}
	
}	

sub getcomments{
	my @LC = @_;
	my $comments=''; 
	my $commentsover=0;

	for my $line (@LC){
		print $line ;
		if ($line=~'@'){
			$commentsover=1;
			
			}
		if (0==$commentsover){
			if ($line=~'#'){
				$comments=$comments.$line."\n";
				print "jelp \n";
				}
			}
		}

	###Print Comments
	open CF,">",$readmefile or die "Cannot find lc file:$readmefile $_";

	print CF "SNANA version created with command:\n" ;
	print CF "  $0 @ARGV\n" ;
	print CF " by user = " . `whoami` . " \n" ;
	print CF " from SALT2 files in Directory: \n";
	print CF "    $datadir \n" ;
	print CF " on " . `date`. "\n";

#	 print CF "\n \#printing comments read in by $0 from original files\n";
	print CF $comments;
	close CF;
}

sub parselightcurveline{
	my ($line) = @_;

	#print "Will format the light curve variables \n";
		### format light curve variables
	#print "$line, \n";
	#print "LINE\n";
	#print $line ;
	if ($line =~ /^\s*$/){
		return '';}
	my ($mjd,$flux,$fluxerr,$zp,$filt,$magsys)= split(/\s+/,$line);

	#Lines to test fix the code 
	#	print $mjd , "\n";
	#	print $flux , "\n"; 
	#	print $fluxerr, "\n"; 
	#	print $zp , "\n";
	#	print $filt, "\n" ;
	#	print $magsys, "\n" ;
	#	print "\n END\n";
	#print "mjd =$mjd, flux =$flux, fluxerr =$fluxerr, zp =$zp, filt=$filt , magsys =$magsys \n";
	my $snanaflux = &calc_snanaflux($flux,$zp);
	#$magref = $magsys;
	my $snanafluxerr = &calc_snanaflux($fluxerr,$zp);
	my $instrument = (split(/::/,$filt))[0];
	my $filter = (split(/::/,$filt))[1];
	my $mag = &calcmag($flux, $zp);
	my $magerr = &calcmagerr($flux, $fluxerr, $zp);
	#print FH "$mjd \t $filter \t $field \t $flux \t $fluxerr \t $mag \t $magerr", "\n";
	#print $field,"\n";
	#my $line_reformatted = "$mjd \t $filter \t $field\t";# $flux \t $fluxerr\t";
	my $snr=0;
	#print "flux = $flux, fluxerr = $fluxerr, snr =$snr \n";
		#calculate snr
	if (($fluxerr <0.0000000001) or ($flux <0)) {
		$snr = -9.0;
		}
	else {
		$snr= $flux/$fluxerr;
		}

	#print "here we have \n";
		###format the numbers we will print
	my $line_reformatted = sprintf("%4s %-5.3f  %1s %3s  %006.3f    %03.3f     %006.3f  %006.3f   %06.3f     %04.3f", 'OBS:', $mjd  ,$filter,$lcfield ,$snanaflux,$snanafluxerr, $snr ,$mag,$magerr,$zp); 

		###check which filters appear in the light curve
	if ($header_filt!~/$filter/){
		$header_filt=$header_filt.$filter;
		}
	#print "in parselightcurve $header_filt , $filter\n";
	if ($globalfilter!~/$filter/){
		$globalfilter=$globalfilter.$filter;
		}
	#print "filt \t",$filter,"\t",$header_filt, "\n";
	#$lightcurve=$lightcurve.$line_reformatted."\n";
	return $line_reformatted;
}

sub readlcdata{
	###reads the lightcurve data after end 
	my @lc= @_;
	my $startparsing=0;
	#print "reading lightcurve data on $#lc lines \n";
	my $lcurve='';
	my $lineno =0;
	for my $line (@lc){
		#print "on lc lineno = $lineno \n";
		#print $line ;
		if (1==$startparsing) {
			print "Will be parsing lc lines \n";
			#print "NOW\n" ;
			#print $line ;
			my $parsedlc = &parselightcurveline($line);
			#print "done \n";
			if ($parsedlc) {
				$lineno++;
			}
			$lcurve =$lcurve.$parsedlc."\n";
			}
		if ($line=~'#end'){
			$startparsing=1;
			$lineno=0;
			}
		}
	my $nobs = $lineno;
	my @retvar = ($lineno, $lcurve);
	#print "RETVAR" ;
	#print @retvar ;
	return @retvar ;
}

sub formatheader{
	my ($nobs, $snid, $ra, $dec , $mwebv, $redshift_helio, $redshift_cmb,$peakmjd,$field,$foundZ,$foundsn,$founddaymax,$redshift_err) = @_;
	my $rstatus = 'OK';
	if ($redshift_err eq 'NULL') {
		$redshift_err =0.0010;
		$rstatus= 'warning: redshift err fudged';
	}
	my $versionphotometry;
	my $survey = $globalsurvey; #for public 
	my $header ='SURVEY:  '.$survey."\n";
	$versionphotometry = $globalversion ; #for public
	$header =$header.'SNID: '.$snid."\n"; 
	$header =$header.'IAUC: UNKNOWN'."\n"; 
	$header =$header.'PHOTOMETRY_VERSION: '."$versionphotometry \n"; 
		#$header =$header.'SNTYPE: '."\n"; 
	$header =$header.'RA: '.$ra." deg"."\n"; 
	$header =$header.'DECL: '.$dec." deg"."\n"; 
	$header =$header.'MWEBV: '.$mwebv." MW E(B-V)"."\n"; 
	$header =$header.'REDSHIFT_HELIO: '.$redshift_helio." +- $redshift_err    (HELIO)\n"; 
	if ($redshift_cmb =~/NULL/){
		$header =$header.'REDSHIFT_FINAL: '.$redshift_helio." +- $redshift_err    (HELIO)\n";
		$header =$header.'REDSHIFT_Status:'.$rstatus."\n"; }
	else
		{
	$header =$header.'REDSHIFT_CMB: '.$redshift_cmb." +- $redshift_err    (CMB)\n"; 
	$header =$header.'REDSHIFT_FINAL: '.$redshift_cmb." +- $redshift_err    (CMB)\n"; 
	$header =$header.'REDSHIFT_Status:'.$rstatus."\n"; 
	}
	$header =$header.'SEARCH_PEAKMJD: '.$peakmjd."\n"; 
	my $sorted_headerfilter = &sortfilter($header_filt);
	print "File: $snid, $sorted_headerfilter, $header_filt \n";
	#$header =$header.'FILTERS: '.$sorted_headerfilter."\n"; 
		#Changed RB, 02.28.11, Rick wants the same filters
		#in all files in a version
	$header =$header.'FILTERS: '.$specifiedglobalfilters."\n";
	$header =$header.'# ============================================'."\n";
	$header =$header.'# TERSE LIGHT CURVE OUTPUT:'."\n";
	$header =$header.'#'."\n";
	$header =$header.'NOBS: '.$nobs."\n";
	$header =$header.'NVAR: 9'."\n";
	my $lclabels = sprintf("%6s %3s %3s %3s %1s %1s   %1s     %1s     %1s       %1s", 'VARLIST:','MJD','FLT','FIELD','FLUXCAL','FLUXCALERR','SNR','MAG','MAGERR','Zpt');
	$lclabels=$lclabels."\n";
	$header =$header.$lclabels;
	return $header;
}

			#print OL $lightcurve;
sub printsnanaformatfile{
	my ($snid,$header,$lightcurve)  = @_;
	print "printing snanafile for $snid \n";
	my $outfilename =$prefix.'_'.$snid.'.dat';
	my $outfile=$outdir.'/'.$outfilename;
	#print "issaltfile =$issalt2file \t, writing $outfile \n";
	open OL,">",$outfile or die "Cannot write file:$outfile $_";
	print OL $header; 
	print "LIGHT CURV" ; "\n";
	#print $lightcurve ;
	print "END LIGHT CURV" ; "\n";
	print OL $lightcurve;
	#print OL "DONE \n";
	close OL;
	print LF $outfilename ,"\n";
}
#close FH;
#close LZ;
#close HZ;

close LF;
close IL;
close NS;
print "the list of all filters used is $globalfilter \n";
#print "Ending convertsalt2snana.pl gracefully for survey = $survey and version = $version \n";


sub recordcandidate {
	my ($num,$fl) = @_;
	if (1==$num){
		open FH, ">", 'candidatesaltfile.txt'or die;
		}
	else{
		open FH, ">>", 'candidatesaltfile.txt'or die;
		}
	print FH "filename $fl \n";
	close FH;
}

sub calc_snanaflux {
	my ($flux,$zp)= @_;
	#print "FLUX , zp", $flux , $zp, "\n";
	my $snanaflux = $flux*(10**(0.4*(27.5 -$zp)));
		###From Rick's post on the wiki 
		###http://supernovae.in2p3.fr/jla/doku.php?id=posts:salt2_vs._snana_format
	return $snanaflux;
}

sub calcmag {
	my ($flux,$zp) = @_;
	my $mag;
	if ($flux>0){
		$mag = -2.5*(log($flux)/log(10)) +  $zp;
	}
	else {$mag =666;}
return $mag;
}

sub qsalt2file{
	my ($var1,$var2,$var3)=@_;
	if (($var1 && $var2 && $var3) ==1) {
		print "This is a salt2file \n";
		return 1;
		}
	else {
		print "This is not a salt2file \n";
                return 0;
	}
}

sub calcmagerr{
	my ($flux, $fluxerr, $zp) =@_;
	my $magerr;
	if ($flux >$fluxerr) {
		my $magpluserr= &calcmag($flux+$fluxerr, $zp);
		my $magminuserr= &calcmag($flux-$fluxerr, $zp);
		$magerr=abs($magpluserr-$magminuserr)/2.;
	}
	else {$magerr=-9.0;}
return $magerr;
}

sub getsnid{
	my ($tmp)= @_;
	my ($snid,$casesnid);
	if ($tmp=~/sn.*/){
		my @tmparr=split(/sn/,$tmp);
		$casesnid = $tmparr[1];
		}
	else
		{
		$casesnid=$tmp;
		}
	if ($casesnid=~/\s*\d\d\d\d[a-z]\s*$/){
		$snid = uc($casesnid);
		}
	else {
		$snid = $casesnid;
		}
		
	return $snid;
}
	
		
sub sortfilter {
	my ($filter) = @_;
	my $outfilt='';
	for my $f (@allfilters){
		if ($filter =~$f){
			$outfilt=$outfilt.$f;
			}
		}	
	return $outfilt;
}

sub locatefile{
	#call syntax &locatefile (directory, pattern) returns a list of 
	#filenames found in directory matching pattern. 
	
	my ($dir,$patt) = @_;
	my $file;
	my @flist;
	$locatecalls=$locatecalls+1;

	my $found=0;
	opendir (DIR, $dir) or die "Error in opening directory $dir:_$ \n";
		#print $dir,"\n";
	while ($file = readdir(DIR)){
		#print $file,$patt,"\n";
		if ($file =~ m/^$patt/) {
			if ($file=~/^\./){next;}
			#print $file,"\n";
			#print GH $patt, "\t", $file ,"\n";
			#print GH $file,"\n";
			$found =$found+1;
			@flist=(@flist,$file);

		}
	}   
	closedir(DIR);
	if ($found==0) {
		#	print NF "on call $locatecalls Not found","\t", $patt,"\n";
	}
#close (NF);
	#print @filelist;
	return @flist;
}
sub calc_zcmb{
	my ($ra,$dec,$zhel)=@_;
	my $bapex = 1;
	my $vapex = 1;
	my $lapex = 1;
	my $rads=1;
	my ($vdotn,$lightc)=0;
	my $ss = sin($rads*$dec)*sin($rads*$bapex);
	my $ccc= cos($rads*$dec)*cos($rads*$bapex)*cos($rads*($ra-$lapex)); 

}
