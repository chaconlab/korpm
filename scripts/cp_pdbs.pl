#!/usr/bin/perl
use strict;



my %pdbsN;
my $cmd;
# read korp mutation file

open(MYFILE, $ARGV[0]);
system("mkdir $ARGV[1]");
while(<MYFILE>) {
  my @f = split(/\s+/, $_);
  $pdbsN{$f[0]}++;
  $cmd ="cp /home/pablo/korpm/thermomutdb/pdbs/$f[0].pdb $ARGV[1]/\n";
  $cmd ="cp /home/pablo/PDBs/$f[0].pdb $ARGV[1]/\n";
  $cmd ="cp /home/pablo/PDBs/$f[0].pdb $ARGV[1]/\n";
  $cmd ="cp Db/$f[0].pdb $ARGV[1]/\n";
  $cmd ="cp ../R1/Db/$f[0].pdb $ARGV[1]/\n";
  system($cmd);
}
close(MYFILE);

