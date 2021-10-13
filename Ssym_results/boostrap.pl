#!/usr/bin/perl
use strict;
use warnings;

use List::Util qw/shuffle/;



# usage  boostrap.pl KORPM_Pucci2018_sym.txt

my @myP;
my @myN;

my $line;

open(INPUT, "<", $ARGV[0]) or die "\nFailed to open $ARGV[0]\n";

 my $lineo;
my $cont=0;
while(<INPUT>)
{
 next if ($_ =~ /^#/); # Skip "#" begining lines

 # Parse input file to get relevant data
 $line = $_;
 chomp($line);

 # Get data from Mutations Database
 my @dbdata = split /\s+/, $line; # Load mutation data array
  
 if (($dbdata[4] =~ /-/) ) {
  $lineo = "$dbdata[0] $dbdata[2] $dbdata[4] $dbdata[5]"; 
  push (@myN, $lineo);
  $dbdata[4] *=-1.0;
  $lineo = "$dbdata[1] $dbdata[3] $dbdata[4] $dbdata[6]"; 
  push (@myP, $lineo);
 }
 else {
  $lineo = "$dbdata[0] $dbdata[2] $dbdata[4] $dbdata[5]"; 
  push (@myP, $lineo);
  $dbdata[4] *=-1.0;
  $lineo = "$dbdata[1] $dbdata[3] $dbdata[4] $dbdata[6]"; 
  push (@myN, $lineo);
 }

 }
 close(INPUT);
 
 
 my @arr;
 
for(my $i = 0; $i <= $#myN; $i++)
{
@arr[$i]=$i;
} 



for(my $b = 0; $b <= 10; $b++) 
{

@arr = shuffle(@arr);
# print @arr;



 
 #  print "deteted $#myN  $#myP\n";
 open(OUTPUT, ">", "kk$b.txt") or die "\nFailed to open $ARGV[0]\n";
 for(my $i = 0; $i <= $#myN; $i++)
{
  if ($i<86) {
  print OUTPUT "$myP[$arr[$i]]\n";
  }
  else {
  print OUTPUT "$myN[$arr[$i]]\n";
  }
}
 close(OUTPUT);
 
 system ("~/korpm/Mstat.pl kk$b.txt 3 4 2 | grep \"X   342\"\n");
 system ("~/korpm/confusion.pl kk$b.txt 3 4 1000 \> kkC$b\_prc.txt \n");
 
}
 

