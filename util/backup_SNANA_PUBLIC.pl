#!/usr/bin/perl
#
# Created Aug 2011 by R.Kessler
# re-make back_SNANA_PUBLIC.cmd -> perl
#
# Usage
#   backup_SNANA_PUBLIC.pl <SNANA version>
#
#  where <version> refers to a cut version of snana.
#  The output is two files:
#    $SNANA_PUBLIC_LOCAL/downloads/SNANA.tar.gz
#    $SNANA_PUBLIC_LOCAL/VERSION.SNANA   (contains version number) ??
#
# and the previous version is moved to  downloads/archive
# with the version appended to the tarball name.
# => ensures that there is only one tarball in the /downloads area/
#
# Note that files are copied, rather than tar'ed,
# so that empty /obj & /bin dirs can be made without
# the hastle of EXCLUDE files.
#
# Oct 28, 2011: copy manual to web-page 
#               (left out by acident when changing from csh to perl)
#
# Oct 23, 2012: no more STAMP SDSS -> PUBLIC since SDSS programs are removed
#
# Nov 07, 2012: copy snana_install.pdf in addition to snana_manual.pdf
#               The public install manual was quite stale (Nov 2012)
#
# Jul 27 2014:
#   Backup now involves another stop: scp to another machine.
#   See $SNANA_PUBLIC_LOCAL and $SNANA_PUBLIC_URL
#   See new function archiveURL().
#
# Feb 21 2015: add /eups to cpDirList
#
# =========== Declarations =================

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;


sub parse_args ;
sub miscInit ;
sub archiveOld_local ;
sub archiveNew_local ;
sub archiveURL ;
sub cpList ;

my ($SNDATA_ROOT, $SNANA_DIR, $SNANA_PUBLIC_LOCAL, $SNANA_PUBLIC_URL);
my ($NARG, $VERSION, @MSGERR, $versionFile );
my ($UPS_snanaDir, $TARNAME, $TARFILE, $tarDir, $downloadsDir, $docDir);
my (@mkdirList, @cpdirList, @srcList, @utiList, @docList) ;
my ($sdirSNANA, $sdirSNANA_full);

my (@URL_FILELIST_downloads ) ;
my (@URL_FILELIST_doc) ;

# ========= BEGIN MAIN  =========

&parse_args();

&miscInit();

&archiveOld_local();

&archiveNew_local();

&archiveURL();

# ========= END OF MAIN ==========


# ======================
sub parse_args {

    $NARG = scalar(@ARGV);

    if ( $NARG < 1 ) {
	$MSGERR[0]= " Must give snana version as argument. " ;
	sntools::FATAL_ERROR(@MSGERR);
    }

    $VERSION = $ARGV[0];
    print " Backup snana version '$VERSION' \n";

} # end of parse_args


# ===============================
sub miscInit {

    $SNANA_PUBLIC_LOCAL  = $ENV{'SNANA_PUBLIC_LOCAL'};
    $SNANA_PUBLIC_URL    = $ENV{'SNANA_PUBLIC_URL'};
    print " \$SNANA_PUBLIC_LOCAL = $SNANA_PUBLIC_LOCAL \n";
    print " \$SNANA_PUBLIC_URL   = $SNANA_PUBLIC_URL   \n";

    # define snana dir in SDSS products area.
    $UPS_snanaDir = "/sdss/ups/prd/snana/${VERSION}/Linux" ;
    $TARNAME = "SNANA" ;
    $TARFILE = "${TARNAME}.tar" ;

    # file containing name of version
    $versionFile = "VERSION.SNANA" ; 

    @URL_FILELIST_downloads = 	( "${TARFILE}.gz" , "$versionFile" );

# -----------------------
# cat PUBLIC.LIST files for list of files to copy

    @mkdirList = ( "src" ,  "obj" ,  "bin",  "lib" , "util",  "doc" ) ;
    @cpdirList = ( "ups" ,  "eups",  "kumacs" );
    @srcList   = `cat ${UPS_snanaDir}/src/PUBLIC.LIST`  ;
    @utiList   = `cat ${UPS_snanaDir}/util/PUBLIC.LIST` ;
    @docList   = `cat ${UPS_snanaDir}/doc/PUBLIC.LIST`  ;


    # define tarDir to be directory for SNANA tarball
    $downloadsDir = "downloads" ;
    $docDir       = "doc" ;
    $tarDir       = "$SNANA_PUBLIC_LOCAL/$downloadsDir";

} # end of miscInit


# =====================
sub archiveNew_local {

    # make archive of snana $VERSION

    my ($cdd, $sdir, $cdsrc );
    my ($cmd1, $cmd2, $cmd3, $cpManual, $manual, $inF, $DOCDIR);

# define snana directory
    $sdirSNANA       = "${TARNAME}_${VERSION}" ;
    $sdirSNANA_full  = "${tarDir}/${sdirSNANA}" ;
    $cdd             = "cd $sdirSNANA_full";
    $cdsrc           = "cd $sdirSNANA_full/src";

# create new VERSION file
    $cmd1 = "rm $versionFile" ;
    $cmd2 = "echo $VERSION > $versionFile" ;
    qx(cd $tarDir ; $cmd1 ; $cmd2 ) ;


    # create temp snana topdir, and needed subdirs
    print " Create archive snana dir in : \n" ;
    print " \t \$SNANA_PUBLIC_LOCAL/$downloadsDir/$sdirSNANA \n";
#    print "  $sdirSNANA_full \n" ;
    print "\n" ;

    qx(mkdir $sdirSNANA_full);
    
    foreach $sdir ( @mkdirList ) {
	print " Create subdir '${sdir}' \n" ;
	qx($cdd; mkdir $sdir;  chmod -R +w $sdir/ );
    }

    foreach $sdir ( @cpdirList ) {
	print " Copy subdir '${sdir}' \n" ;
	qx($cdd ; cp -r $UPS_snanaDir/$sdir $sdir;  chmod -R +w $sdir/ );
	qx($cdd ; cd $sdir ;  rm -rf CVS */CVS Makefile  */Makefile);
    }


    # now copy public source codes
    cpList("src",  @srcList) ;
    cpList("util", @utiList) ;
    cpList("doc",  @docList) ;

    # replace SDSS stamp with PUBLIC stamp
    $cmd1 = "sed -e 's/STAMP := SDSS/STAMP := PUBLIC/g' Makefile > TMPMAKE";
    $cmd2 = "mv TMPMAKE Makefile";
##    qx($cdsrc ; $cmd1 ; $cmd2 );

    # ----------------------------------------
    # make tarball

    print "\n Create tarball : $sdirSNANA -> ${TARFILE}.gz\n" ;
    $cdd  = "cd $tarDir" ;
    $cmd1 = "tar -cf ${TARFILE} $sdirSNANA" ; 
    $cmd2 = "gzip ${TARFILE}" ;
    $cmd3 = "rm -rf $sdirSNANA" ;
    qx($cdd ; $cmd1 ; $cmd2 ; $cmd3 );

    # xxxxxxxxxxxxxx
    # copy Joe's web stuff to cvs since SNANA_PUBLIC is not backed up
    # print " Backup html files ... \n";
    # qx(cp $SNANA_PUBLIC_LOCAL/*.*  ~/snana/web/ ) ;
    # xxxxxxxxxxxxxx

    
    # update the manuals 
    @URL_FILELIST_doc =  ( "snana_manual.pdf", "snana_install.pdf" );

    $DOCDIR = "$SNANA_PUBLIC_LOCAL/$docDir" ;
    foreach $manual ( @URL_FILELIST_doc ) {
	$inF       = "$UPS_snanaDir/$docDir/$manual" ;
	$cpManual  = "cp $inF  $DOCDIR" ;
	qx($cpManual);
    }


} # end of archiveNew


# =================
sub cpList {

    my ($sdir, @fileList) = @_ ;
    # copy each file from $sdir into archive-sdir

    my ($file, $cpCmd, $NCOPY, $cN, $privCmd);

    $NCOPY = scalar(@fileList);
    $cN    = sprintf("%2d" , $NCOPY );
    print "\t Copy  $cN  '$sdir' files ...\n";

    foreach $file ( @fileList ) {
	$file  =~ s/\s+$// ;   # trim trailing whitespace
	$cpCmd   = "cp $UPS_snanaDir/$sdir/${file}  $sdirSNANA_full/$sdir/" ;
	qx($cpCmd);	
    }

    # set write-priv on all copied files.
    $privCmd = "chmod -R +w $sdirSNANA_full/$sdir/*" ;
    qx($privCmd);

} # end of cpList



# =====================
sub archiveOld_local {

    # backup current archive in /archive.

    my ($OLD_VERSION, $archiveFile );

    $OLD_VERSION = qx(cd $tarDir ; cat $versionFile);
    $OLD_VERSION =~ s/\s+$// ;   # trim trailing whitespace
    $archiveFile = "archive/SNANA_${OLD_VERSION}.tar.gz";

    print " Move previous SNANA tarball to  $archiveFile \n" ;
    qx(cd $tarDir ; mv ${TARNAME}.tar.gz  ${archiveFile} );


} # end of archiveOld_local


# =========================
sub archiveURL {

    # Created July 28 2014 
    # copy to separate linux box defined by
    # ENV $SNANA_PUBLIC_URL
    # [need to get away from dp62 machine]

    my ($file, $scp, $inFile, $outDir);

    print "\n scp files to public archive:\n   $SNANA_PUBLIC_URL \n";

    
    foreach $file ( @URL_FILELIST_downloads ) {
	$inFile  = "$SNANA_PUBLIC_LOCAL/$downloadsDir/$file";
	$outDir  = "$SNANA_PUBLIC_URL/$downloadsDir" ;
	$scp     = "scp $inFile $outDir" ;
	print "\t scp $file \n";
	qx($scp);
    }

    foreach $file ( @URL_FILELIST_doc ) {
	$inFile  = "$SNANA_PUBLIC_LOCAL/$docDir/$file";
	$outDir  = "$SNANA_PUBLIC_URL/$docDir" ;
	$scp     = "scp $inFile $outDir" ;
	print "\t scp $file \n";
	qx($scp);
    }

} # end of archiveURL
