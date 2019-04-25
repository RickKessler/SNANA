#!/usr/bin/perl -w
#
# Radomly re-order the GAL entries in a SNANA "HOSTLIB".
# Motivated by the LSST catalog which comes strangely ordered
# David Cinabro 28 June 2011
#
# USAGE: >galaxy_randomizer.pl input.dat output.dat
# NOTES: It needs both input and output specificed.
#        If both not given it exits with usage notes.
#        output can equal input in which case it will overwrite
#        the old file.
#        It expects that all the header material appears before the data
#        and that the data all begins with the keyword GAL:
#        To be sure that the file is SNANA galaxy catalog it checks
#        for the keyword NVAR: in the header.  If not found it
#        aborts with an error message.
#        It reads in the input dividing it into a header and data.
#        It then randomizes the order of the data using a well known,
#        highly efficient algorithm.
#        It then writes the header followed by the ransomized data
#        to the output.
#        If one would like to debug with a fixed randomization, uncomment
#        the srand line to give the random number generator a fixed random
#        seed.

if ($#ARGV != 1) {
   print ("Usage:  galaxy_randomizer.pl INPUT.DAT OUTPUT.DAT\n");
   print ("Note: if INPUT.DAT and OUTPUT.DAT are the same\n");
   print ("      it will overwrite with randomly ordered data.  Take care.\n");
   exit;
} 

print ("Randomly re-order the data in a SNANA HOSTLIB galaxy catalog. \n");

$INPUT_FILE = $ARGV[0];
$OUTPUT_FILE = $ARGV[1];

# read input file

open(INPUT_FILE);
@array = <INPUT_FILE>;
close(INPUT_FILE);

# divide input file into header and data

$dcount = -1;
$hcount = -1;
foreach $line (@array) {
  $pos = index($line,"GAL:");
  if ($pos != 0) {
    $hcount = $hcount + 1;
    $head[$hcount] = $line;
  }
  else {
    $dcount = $dcount + 1;
    $data[$dcount] = $line;
  }
}

# check header for NVAR:

$snanafile = 0;
foreach $line (@head){
  $pos = index($line,"NVAR:");
  if ($pos != 0) {
     next;
  }
  elsif ($pos == 0) {
    $snanafile = 1;
  }
}
if ($snanafile != 1) {
    print("galaxy_randomizer.pl cannot find the NVAR: keyword\n");
    print("in the header of ",$INPUT_FILE,".  This file does not\n");
    print("look like a SNANA galaxy host file.  Stopping.\n");
    exit;
}

# randomize the data
print ("Randomizing ",$dcount+1," galaxy entries in ",$INPUT_FILE,"\n");

# fisher_yates_shuffle( \@array ) : 
# generate a random permutation of @array in place
# stolen from the PERL FAQ
sub fisher_yates_shuffle {
  my $array = shift;
  my $i;
  for ($i = @$array; --$i; ) {
    my $j = int rand ($i+1);
      next if $i == $j;
      @$array[$i,$j] = @$array[$j,$i];
  }
}

#uncomment next line to get a fixed random seed for debugging
#srand(123456789);
fisher_yates_shuffle( \@data);

# write output

open(OUT,">".$OUTPUT_FILE);
print OUT @head;
print OUT @data;
close(OUT);

print ("Done Randomizing galaxy entries.  Written to ",$OUTPUT_FILE,"\n");

exit;
