#!/usr/bin/perl
use strict;

my @list=("CartDD","FoldX","EvoFF","Dynamut2","PopMs","DDgun3D","ThermoNet","ACDCNN");

system ("rm out.txt temp temp2; touch out.txt");

my $cmd=" sort -k 1 -k 2 S461_KORPM_DIR.txt > temp2;    awk \'{printf \"%4s %-7s %7.3f %7.3f\\n\",\$1, \$2,\$10,\$11}\' temp2 > out.txt"; 
print $cmd;

system ($cmd);

foreach (@list) {
#my $cmd=" sort -k 1 -k 2 S461_$_\_DIR.txt  > temp2;    awk \'{printf \" %7.3f %7.3f\\n\", \$3, \$4}\' temp2 > temp"; 
my $cmd=" sort -k 1 -k 2 S461_$_\_DIR.txt  > temp2;    awk \'{printf \"%7.3f\\n\", \$4}\' temp2 > temp"; 
system ($cmd);
system ("paste out.txt temp > temp2; mv temp2 out.txt");
 
print $cmd;

}
