#!/usr/bin/perl
use strict;
use Cwd qw();

# usage  ch_val.pl pdbs/Pucci2018dirN.txt
# ch_val.pl pdbs/Pucci2018invN.txt  > Pucci2018invN.txt
# ch_val.pl pdbs/Pucci2018dirN.txt  > Pucci2018dirN.txt
# ch_val.pl pdbs/Pucci2018N.txt  > Pucci2018N.txt


my @myAA;
my @myAA0;
my $line;

open(INPUT, "<", $ARGV[0]) or die "\nFailed to open $ARGV[0]\n";
while(<INPUT>)
{
 next if ($_ =~ /^#/); # Skip "#" begining lines

 # Parse input file to get relevant data
 $line = $_;
 chomp($line);

 # Get data from Mutations Database
 my @dbdata = split /\s+/, $line; # Load mutation data array
 
 # Get mutation data stuff
 push (@myAA0, \@dbdata);

 }
close(INPUT);


open(INPUT, "<", "changes.txt") or die "\nFailed to open changes.txt\n";

my $cont=0;
while(<INPUT>)
{
 next if ($_ =~ /^#/); # Skip "#" begining lines

 # Parse input file to get relevant data
 $line = $_;
 chomp($line);

 # Get data from Mutations Database
 my @dbdata = split /\s+/, $line; # Load mutation data array
 
 # Get mutation data stuff
 push (@myAA, \@dbdata);
 
 my $tag= $dbdata[2];
 
 
 
 my $wt = substr( $dbdata[1], 0, 1 ); # WT aminoacid

 my $chain = substr( $dbdata[1], 1, 1 ); # chain aminoacid
 
 my $mut = substr( $dbdata[1], -1, 1 ); # Mut aminoacid
 
 my $pos = substr($dbdata[1], 2); 
 chop($pos);

 my $rever=$mut.$chain.$pos.$wt;

 #printf("%4s %-7s\n",  $dbdata[1], $rever); 
   

for (my $k=0; $k <= $#myAA0; $k++) {

 if ($dbdata[1] eq $myAA0[$k][1]) {   
 $myAA0[$k][2]=$dbdata[3];
# printf("%4s %-7s %6.3f  %6.3f  %6.3f %-7s \n", $dbdata[0], $dbdata[1], $dbdata[2],$myAA0[$k][2],$dbdata[3], $myAA0[$k][1]); 
 } 
 if ($rever eq $myAA0[$k][1]) {   
 $myAA0[$k][2]=$dbdata[3]*-1.0;
# printf("%4s %-7s %6.3f  %6.3f  %6.3f %-7s \n", $dbdata[0], $rever, $dbdata[2],$myAA0[$k][2],$dbdata[3], $myAA0[$k][1]); 
 }  
 
}


 
 }
 
 
 for (my $k=0; $k <= $#myAA0; $k++) {
 # for ORIGINAL Pucci_files
 # printf("%4s %-7s %6.3f %3d %6.2f\n", $myAA0[$k][0], $myAA0[$k][1], $myAA0[$k][2],$myAA0[$k][3], $myAA0[$k][4]);
 # printf("%4s %-7s %6.3f\n", $myAA0[$k][0], $myAA0[$k][1], $myAA0[$k][2]);
  
  printf("%4s %-7s %7.3f %7.3f\n", $myAA0[$k][0], $myAA0[$k][1], $myAA0[$k][2],$myAA0[$k][3]);
 
 } 
 
 
 
