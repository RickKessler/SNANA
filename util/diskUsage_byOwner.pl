#!/usr/bin/perl
#
# Created Dec 2016
# Make table of disk usage by own. Default directory
# to check is $SNDATA_ROOT/SIM. Ouput file has following
# format:
#  [user]  [Nfile]   [space(Mb)]
#
#
# Usage:
#   diskUsage_byOwner.pl 
#     [examine default directory $SNDATA_ROOT/SIM]
#
#   diskUsage_byOwner.pl <dirName>
#     [examine dirName instead of $SNDATA_ROOT/SIM]
#

use strict ;
use IO::Handle ;

# declarations
my $outFile     = "DISK_USAGE.LOG";
my $SNDATA_ROOT = $ENV{'SNDATA_ROOT'};
my $TOPDIR      = "$SNDATA_ROOT/SIM";
my $CDATE    = `date +%Y-%m-%d_%H%M` ;  $CDATE =~ s/\s+$// ;
my $OUTFILE  = "$TOPDIR/$outFile" ;

my (@SORTED_USERLIST, $USER, @FILELIST, $NFILE_TOT, $SIZEMB_TOT, $SIZEGB_TOT );

sub parseArgs ;
sub open_OUTFILE ;
sub close_OUTFILE ;
sub get_USERLIST ;
sub get_USERDISK ;

sub uniq ;

# ============ START MAIN ============

&parseArgs();

&open_OUTFILE();

print " Get disk usage by owner in \n\t $TOPDIR\n";

&get_USERLIST();

@FILELIST = qx(cd $TOPDIR; ls -Ro );
$NFILE_TOT=0 ; $SIZEMB_TOT=0.0;

foreach $USER ( @SORTED_USERLIST) {
    &get_USERDISK($USER);
}

&close_OUTFILE();

print "\n See summary in\n\t $OUTFILE\n";

# ========= END MAIN =============

sub parseArgs {

    my $NARG = scalar(@ARGV);
    if ( $NARG == 1 ) {
	$TOPDIR  = $ARGV[0];
	$OUTFILE = "$TOPDIR/$outFile"; 
    }
}


# ===============================
sub open_OUTFILE {

    print " Open $OUTFILE \n";

    open  PTR_OUTFILE , "> $OUTFILE" ;
    print PTR_OUTFILE " $CDATE\n Get disk usage by owner in \n\t $TOPDIR\n\n ";

    print PTR_OUTFILE "    User                NFile  space(Mb) \n";
    print PTR_OUTFILE " -------------------------------------------- \n";
}


sub close_OUTFILE  {

    $SIZEGB_TOT = int($SIZEMB_TOT/1000.) ;

    print PTR_OUTFILE " -------------------------------------------- \n";
    print PTR_OUTFILE "\n";
    print PTR_OUTFILE " Total number of files: $NFILE_TOT \n";
    print PTR_OUTFILE " Total size   of files: $SIZEGB_TOT  GigaBytes\n";

close PTR_OUTFILE ;
} 


# =======================================
sub get_USERLIST{

    my (@TMPLIST_LINES, @TMPLIST_USER, $tmp, $N, @wdlist);

    @TMPLIST_LINES = qx(cd $TOPDIR; ls -o);
    $N=0;
    foreach $tmp (@TMPLIST_LINES) {
	@wdlist = split(/\s+/,$tmp ) ;	
	if ( scalar(@wdlist) < 3 ) { next; }
	$TMPLIST_USER[$N] = $wdlist[2];
	$N++ ;
    }
    
    my @USERLIST = uniq(@TMPLIST_USER);

    @SORTED_USERLIST = sort(@USERLIST);

    print "\n sorted list of owners: \n  @SORTED_USERLIST \n" ;

}     # end get_USERLIST


# =======================================
sub get_USERDISK {
    my ($user) = @_ ;

    my (@tmp, @wdlist, $line, $NFILE, $size, @matches);

    print " ===> Examine disk usage for $user \n";

    my @matches = grep { / $user / } @FILELIST;

###    @tmp = qx(cd $TOPDIR; find . -user $user -type f -exec du -smc {} +) ;

    $NFILE = scalar(@matches);
    $size = 0.0 ;

    # sum disk sizes
    foreach $line ( @matches) {
	@wdlist = split(/\s+/, $line ) ;
	$size += $wdlist[3];
    }
    $size /= 1.0E6 ;  # convert to MB


    # update outFile
    my $txtUser = sprintf(" %-20.20s   %6d   %6d", 
			  $user, $NFILE, $size);

    print PTR_OUTFILE "$txtUser\n";

    $NFILE_TOT   += $NFILE ;
    $SIZEMB_TOT  += int($size); 

}   # end get_USERDISK

# =======================================
sub uniq { my %seen;  grep !$seen{$_}++, @_ ; }
