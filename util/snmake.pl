#!/usr/bin/perl
#
# Created Aug 2011 by R.Kessler
# [translated from C-shell snmake and updated with recent libs ]
#
# Utility to compile and link private ${jobname}.cra file 
# to create ${jobname}.exe 
#
# Usage:
#  snmake.pl snana_private
#  snmake.pl snana_private  --lib  MY_STRANGELIB.a
#  snmake.pl snana_private  --lib  MY_STRANGEOBJ.o --lib SECONDLIB.o
#  snmake.pl snlc_fit_private
#  snmake.pl snlc_fit_private  NOCLEAN  ! leave .lnk, ypatch.log, etc ...
#  snmake.pl psnid_private
#
#
# NOTES:
#
#  *  the --lib argument can be a library or an object file
#
#  * if the jobname contains the string 'fit', then libsnfit.a is used
#    to link; otherwise libsnana.a is used. Also, the 'fit' key uses
#    the genmag_${GENMAG}.o objects.
#
#
#
#   HISTORY
#
# Apr 21, 2014: replace ypatchy with local ypatchy.pl;
#               no need for ':' since file-name case is preserved.
#
# Nov 28 2014:  add ROOTLIB to linking.
#
# Feb 19 2015: add miniut.o to list of object files. Remove MNLIB.
#
# Mar 2 2017: add sntools_gridread.o to list of .o files
#             Needed to make snlc_fit_private.
#
# Feb 15 2018: few fixes in get_SNobj (after fluxerrmap refactor)
#
# Nov 7 2018: added -m64 and -O1 flags for compilation
#
# Jan 3 2019: 
#  + fix link command in sub link().
#  + if not MAKEFIT, remove some .o files from OBJTOOLS (sub get_SNOBJ)
#
# --------------------------------------------------------

use FindBin qw($Bin);
use lib "$Bin" ;
use sntools ;
use strict ;

# ====================
sub init_snmake ;
sub parse_args ;
sub get_SNOBJ ;
sub get_SNANALIB ;
sub get_EXTLIBDIR ;
sub run_patchy ;
sub makeObjFile ;
sub link;
sub checkTimeStamp ;

my ($NARG, $MAKESNANA, $MAKEFIT, $MAKEPSNID, $HELP, $CLEANFLAG ) ;
my ($jobname, $JOBNAME );
my ($MYLIBS, $MYBIN, $UNAME) ;
my ($FC, $LMFLAG, $FFLAGS) ;
my ($SNANA_DIR, $CERN_DIR, $ROOT_DIR, $MINUIT_DIR, $CFITSIO_DIR );
my ($OBJDIR, $LIBDIR, $SNOBJ, $SNANALIB, $EXTLIBDIR, $LG2C ) ;

# ================ BEGIN MAIN ====================

&init_snmake();

# parse user args
&parse_args();

$JOBNAME = "$MYBIN/${jobname}.exe" ;
if ( -e $JOBNAME) { qx(rm $JOBNAME) ; }

# get list of object files
&get_SNOBJ();

# get list of libraries to use
&get_SNANALIB();

# get libdir for external libraries
&get_EXTLIBDIR();

# make .f file with patchy
&run_patchy();

# make .o file
&makeObjFile();

# link (make .exe binary)
&link();

# make  sure that exe has recent time stamp.
&checkTimeStamp();



# ===========================
#
#  END OF MAIN
#
# ===========================


# ===========================
sub init_snmake {

    $HELP   = 0 ;
    $CLEANFLAG = 1; 
    $MYLIBS = ""   ;
    $MYBIN  = "." ;
    $FC     = "gfortran" ;
    $UNAME  = `uname -s` ;

    $FFLAGS    = "-c -static -fno-automatic -m64 -fsecond-underscore -O1" ;


    $SNANA_DIR    = $ENV{'SNANA_DIR'} ;
    $CERN_DIR     = $ENV{'CERN_DIR'} ;
    $ROOT_DIR     = $ENV{'ROOT_DIR'} ;
    $MINUIT_DIR   = $ENV{'MINUIT_DIR'} ;
    $CFITSIO_DIR  = $ENV{'CFITSIO_DIR'} ;

    $OBJDIR = "$SNANA_DIR/obj" ;
    $LIBDIR = "$SNANA_DIR/lib" ;    

} # end of init_snmake


# ===========================
sub parse_args {

    $NARG = scalar(@ARGV) ;
    
    if ( $NARG == 0 ) {
	$HELP = 1 ;
    }
    else {
	if ( "$ARGV[1]" eq "help"    ) { $HELP      = 1 ; }
	if ( "$ARGV[1]" eq "NOCLEAN" ) { $CLEANFLAG = 0 ; }
	if ( "$ARGV[1]" eq "noclean" ) { $CLEANFLAG = 0 ; }
    }


    if ( $HELP ) {
	my $ctmp;
	$ctmp = "snmake.pl mysn" ;
	print " $ctmp                (uses mysn.cra to create mysn.exe) \n" ;
	print " $ctmp --bin MYBIN    (idem, but creates MYBIN/mysn.exe) \n" ;
	print " $ctmp --lib MYLIBS   (idem, but link to MYLIBS) \n" ;
	die "\n";
    }

    $jobname = $ARGV[0] ;

    # parse remaining items
    my $i;
    $i = 0 ;
    while ( $i < $NARG ) {
	$i++ ;
	if ( "$ARGV[$i]" eq "--bin" ) {
	    $MYBIN  = "$ARGV[$i+1]" ;
	}
	if ( "$ARGV[$i]" eq "--lib" ) {
	    $MYLIBS  =  "$MYLIBS " . "$ARGV[$i+1]" ;
	}
    }


    # determine if this job is a fitter job or snana job
    my $lcJob = lc($jobname);
    $MAKESNANA   = $MAKEFIT = $MAKEPSNID = 0 ;
    if ( index($lcJob,"fit")   >= 0 ) { $MAKEFIT   = 1 ; }
    if ( index($lcJob,"psnid") >= 0 ) { $MAKEPSNID = 1 ; } 
    if ( ( $MAKEFIT==0 ) && ($MAKEPSNID==0) ) { $MAKESNANA = 1; }

    print " MakeFlag(SNANA,FIT,PSNID) = $MAKESNANA, $MAKEFIT, $MAKEPSNID \n";

} # end of parse_args



# ============================
sub get_SNOBJ {

    # fill global SNOBJ with .o files
    # May 22, 2012: add genSmear_models.o for MAKEFIT option
    # Jan 03, 2019: check $MAKEFIT

    my (@objList, $OBJTOOLS, $genmag_OBJLIST,  $obj) ;

    # be careful because snana_private links to libsnana.a,
    # but snlc_fit does not.

    if ( $MAKEFIT ) {
	$OBJTOOLS  = 
	    "${OBJDIR}/sntools.o " . 
	    "${OBJDIR}/sntools_fluxErrModels.o " .
	    "${OBJDIR}/sntools_fluxErrModels_legacy.o " .
	    "${OBJDIR}/sntools_nonlinearity.o " .
	    "${OBJDIR}/sntools_fitsio.o " ;
    }

    # tools for both snana and snlc_jobs
    $OBJTOOLS = "$OBJTOOLS " .
	"${OBJDIR}/sntools_gridread.o " .
	"${OBJDIR}/sntools_output.o " .
	"${OBJDIR}/sntools_nearnbr.o " .
	"${OBJDIR}/MWgaldust.o " .
	"${OBJDIR}/multiseason.o " .
	"${OBJDIR}/minuit.o " 
	;


    $SNOBJ = "$OBJTOOLS" ;

    if ( $MAKEPSNID ) {
	$SNOBJ = 
	    "$SNOBJ " . 
	    "${OBJDIR}/psnid_tools.o " .   
	    "${OBJDIR}/psnid_BEST.o "
	    ;
    }

    if ( $MAKEFIT ) {
	$genmag_OBJLIST = 
	    qx(cd $OBJDIR ; ls genmag*.o gloes.o sntools_genSmear.o) ;
	@objList   = split(/\s+/,$genmag_OBJLIST) ;
	foreach $obj ( @objList ) {
	    $obj   =~ s/\s+$// ;   # trim trailing whitespace
	    $SNOBJ = "$SNOBJ " . "${OBJDIR}/$obj" ;
	}
    }

} # end of get_SNOBJ


# ============================
sub get_SNANALIB {

    if ( $MAKESNANA  || $MAKEPSNID ) 
    { 	$SNANALIB = "${LIBDIR}/libsnana.a"; }
    elsif ( $MAKEFIT   ) 
    {  	$SNANALIB = "${LIBDIR}/libsnfit.a"; }
    else 
    {	$SNANALIB = "" ;  }

} # end of get_SNOBJ


sub get_EXTLIBDIR {

    my ($BITNESS, $TESTDIR, $J64);

    $EXTLIBDIR = "lib";  # default

    # check for bitness-dependent libs
    $BITNESS = qx(uname -p);
    $BITNESS   =~ s/\s+$// ;   # trim trailing whitespace
    $TESTDIR = "${CERN_DIR}/lib$BITNESS" ;
    if ( -d $TESTDIR )  { $EXTLIBDIR = "lib${BITNESS}" ; }

    $J64 = index($BITNESS,"64");
    $LMFLAG = "-m32" ; 
    if ( $J64 > 0 ) { $LMFLAG = "-m64" ; }

    # check if libg2c exists
    $LG2C = "" ;

#  Remove Aug 8 2014
#    my @BLA = qx(ls /usr/lib/libg2c.*  2>/dev/null ) ;   
#    if ( scalar(@BLA) > 0 )  { $LG2C = "-lg2c"; }
# 

} # end of get_EXTLIBDIR

# ===================
sub run_patchy {

    # create .f file with patchy

    my ($forFile, $craFile, $carFile);

    $forFile = "${jobname}.f" ;
    $craFile = "${jobname}.cra" ;

    print " Use ypatchy.pl to create ${jobname}.f ... \n";
    
    $carFile = "$SNANA_DIR/src/snana.car" ;
    qx(ln -s $carFile     snana.lnk   ); 

    $carFile = "$SNANA_DIR/src/snlc_fit.car" ;
    qx(ln -s $carFile  snlc_fit.lnk);

    if ( $MAKEPSNID ) {
	$carFile = "$SNANA_DIR/src/psnid.car" ;
	qx(ln -s $carFile  psnid.lnk);
    }

    qx(ypatchy.pl - ${forFile}  ${craFile} .go);

    # clean up
    if ( $CLEANFLAG ) 
    { qx(rm *.lnk ; rm ypatchy.log ); }   

    # abort if zero-length file is created
    my @bla = `cat $forFile` ;
    if ( scalar(@bla) == 0 ) {
	die "\n FATAL ERROR: $forFile has zero-length  ?!?!? \n";
    }

} # end of run_patchy

# =====================
sub makeObjFile {

    my (@errFC, $cmdFC, @MSGERR, $err, $ierr );
    
    print "\n";
    print " Begin compilation with: $FFLAGS \n " ;    

    $cmdFC  = "${FC} -c -static $FFLAGS ${jobname}.f";
    @errFC  = qx($cmdFC 2>&1 );

    if ( scalar(@errFC) > 0 ) {
	$MSGERR[0] = "Error compiling $jobname" ;

	$ierr = 0;
	foreach $err ( @errFC ) {
	    $ierr++ ; $MSGERR[$ierr] = "$err" ;
	}
	sntools::FATAL_ERROR(@MSGERR);
    }

}    # end of makeOjbFile

# ===============================
sub link {

    print "\n" ;
    print " Begin linking $MYBIN/${jobname}.exe ... \n" ;
    print " SNANA_DIR    = $SNANA_DIR \n" ;
    print " CERN_DIR     = $CERN_DIR \n";
    print " ROOT_DIR     = $ROOT_DIR \n";
    print " MINUIT_DIR   = $MINUIT_DIR \n";
    print " CFITSIO_DIR  = $CFITSIO_DIR \n";
    print " SNANA LIB    = $SNANALIB \n" ;
    print " PRIVATE LIB  = $MYLIBS \n" ;
    print " \n" ;

    my ($LGSL, $CERNLIB, $ROOTLIB, $SNLIBS, $ALL_OBJ, $ALL_LIBS);
    my ($FITSIOLIB, $LINK, $GSL_DIR);

    # require gsl
    $GSL_DIR  = $ENV{'GSL_DIR'} ;
    $LGSL  =  
	"$GSL_DIR/$EXTLIBDIR/libgsl.a  ". 
	"$GSL_DIR/$EXTLIBDIR/libgslcblas.a" ;

    $FITSIOLIB  = "-L${CFITSIO_DIR}/$EXTLIBDIR -lcfitsio";

    $CERNLIB = $ROOTLIB = "" ;
    if ( length($CERN_DIR) > 0 ) {
	$CERNLIB  = "-L${CERN_DIR}/$EXTLIBDIR -lpacklib -lmathlib -lkernlib -lncurses" ;
    }
    if ( length($ROOT_DIR) > 0 ) {
	$ROOTLIB = "-L${ROOT_DIR}/$EXTLIBDIR -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -ldl" ;
    }

    $SNLIBS     = "${SNANALIB} ${MYLIBS}" ;
    $ALL_OBJ    = "${jobname}.o ${SNOBJ}" ;
    $ALL_LIBS   = "$SNLIBS $FITSIOLIB  $CERNLIB $ROOTLIB $LGSL" ;
    $JOBNAME    = "${MYBIN}/${jobname}.exe";
    

#    $LINK = "${FC} $LMFLAG -lpthread -o  $JOBNAME $LG2C $ALL_OBJ $ALL_LIBS";

# Nov 2018: try to re-arrainge order to fix problem found by Cinabro
    $LINK = "${FC} $LMFLAG -lpthread -o  $JOBNAME ${jobname}.o " .
	"$SNLIBS $SNOBJ " .
	"$SNANALIB " .  # added Jan 3 2019 (for ge2dex & in2dex)
	"$FITSIOLIB  $CERNLIB $ROOTLIB $LGSL" ;

    print " LINK command: \n\t $LINK \n\n";
    
    qx($LINK);

}  # end of link


# ======================
sub checkTimeStamp {

    my $cdate ;
    my $bla ;
    my @wdlist ;

    $cdate = `date`;
    $cdate   =~ s/\s+$// ;   # trim trailing whitespace

    print "\n" ;
    print "******************************* \n" ;
    print " Done linking $MYBIN/${jobname}.exe  \n" ;
    print " Current Date --> $cdate \n" ;
    
    $bla      = qx(ls -l ${MYBIN}/${jobname}.exe) ;
    @wdlist   = split(/\s+/,$bla) ;

    print " $JOBNAME TIME STAMP --> @wdlist[5 ... 7] \n" ;
} 
