#!/usr/bin/perl
use strict;

my @list=("KORPM","Cartddg","FoldX","Evo","Dynamut2","PopMs","DDGun3D","TherNet","ACDCNN");

system ("rm out.txt temp temp2; touch out.txt");

my $cmd=" sort -k 1 -k 2 KORPM_Pucci2018_sym.txt > temp2;    awk \'{printf \"%4s %-7s %7.3f %4s %-7s %7.3f\\n\",\$1, \$3, \$5, \$2, \$4, -\$5}\' temp2 > out.txt"; 
print $cmd;

system ($cmd);

foreach (@list) {
#my $cmd=" sort -k 1 -k 2 $_\_Pucci2018_sym.txt > temp2;    awk \'{printf \"%4s %-7s %7.3f %7.3f %7.3f\\n\",\$1, \$3, \$5, \$6, \$7}\' temp2 > temp"; 
my $cmd=" sort -k 1 -k 2 $_\_Pucci2018_sym.txt > temp2;    awk \'{printf \" %7.3f %7.3f\\n\", \$6, \$7}\' temp2 > temp"; 
system ($cmd);
system ("paste out.txt temp > temp2; mv temp2 out.txt");
 
print $cmd;

}
