#!/usr/bin/perl
use strict;
my %nclass=();
my @nmutclass=();
my %clustPDB=();
my @f;
my %pdbsclass=();
my @pdbs;
my @tag;
my @dde;


my %pdbsN;
my $cmd;
# read korp mutation file

open(MYFILE, $ARGV[0]);

while(<MYFILE>) {
  my @f = split(/\s+/, $_);
  $pdbsN{$f[0]}++;
}
close(MYFILE);

my $npdbs= scalar keys %pdbsN;

print " Number of distint pdbs $npdbs\n";
$cmd = "rm temp; touch temp; rm -r tmp";
system($cmd);
    
for(keys %pdbsN) {
  my $fasta_file="$ARGV[1]/".$_.".fasta";
  $cmd = "cat $fasta_file >> temp\n";
  system($cmd);
}

my $ID=$ARGV[2];
my $C=$ARGV[3];
$cmd ="sed s/_prot/\"\"/g temp | sed s/\'> \'\/\'>tr\|\'\/g > all.fasta";
system($cmd);
  
print " mmseqs clustering ID cutoff $ID $C\n";  

$cmd ="mmseqs easy-cluster all.fasta all tmp --min-seq-id $ID -c $C --cov-mode 0 > log";
system($cmd);
print("$cmd\n");
  



# read mmseqs tsv clastering
open(MYFILE, "all_cluster.tsv");
while(<MYFILE>) {
  chop;
  @f = split(/\t/, $_);
  $nclass{$f[0]}++; 
  push  @{$clustPDB{ $f[0] }}, $f[1];
}
close(MYFILE);

my $npdbs= scalar keys %nclass;
print " Number of clases $npdbs\n";

my @keyC = sort ( { ($nclass{$b} <=> $nclass{$a}) or ( $a cmp $b) } keys %nclass );

my $i=0;



for (my $k=0; $k <= $#keyC; $k++) {         
#print "class $k $keyC[$k] $nclass{$keyC[$k]} ";
for (my $j=0; $j < $nclass{$keyC[$k]}; $j++) {
     $pdbsclass{$clustPDB{$keyC[$k]}[$j]}=$k; 
#     print " $clustPDB{$keyC[$k]}[$j] $k";
     $i++;
  }
#print "\n";  
}

print " Number of pdbs in clases $i\n";




 open(MYFILE2, ">", "cluster.txt")
  or die "cluster.txt";

# read korp mutation file
open(MYFILE, $ARGV[0]);
$i=0; 
while(<MYFILE>) {
  @f = split(/\s+/, $_);
  $pdbs[$i]=$f[0];
  $tag[$i]=$f[1];
  $dde[$i]=$f[2];
  $nmutclass[$pdbsclass{$f[0]}]++; 
  $f[3]= $pdbsclass{$f[0]};
  #print MYFILE2 join( "\t" => @f ), $/;
  
   my $mtdat = sprintf "%4s %-7s %7.3f %-3s %2d %2d %2d %4s %4.2f %2d %-24s %2d %-4s %1d %-15s %1d %-7s %s %4s %5s %6s %6s %5s %10s %-8s %s %s", 
                       $f[0], $f[1], $f[2], $f[3], $f[4], $f[5], $f[6], $f[7], $f[8], $f[9], 
                       $f[10], $f[11], $f[12], $f[13], $f[14], $f[15], 
  			$f[16], $f[17], $f[18], $f[19], $f[20], $f[21], $f[22], $f[23], $f[24], $f[25], $f[26]; # Dump mutation data
 
  
  
  print MYFILE2 "$mtdat\n"; 
  

  $i++;
}
close(MYFILE);
close(MYFILE2);



$i=0; my $nmuts;
$nmuts=0;
for (my $j=0; $j <= $#keyC; $j++) { 
#for (my $j=0; $j < $#nmutclass; $j++) {
  if  ($nmutclass[$j] > 0)  {
  print "class $keyC[$j] $j $nmutclass[$j] $nclass{$keyC[$j]}";
  for (my $jj=0; $jj < $nclass{$keyC[$j]}; $jj++) {
     print " $clustPDB{$keyC[$j]}[$jj]($pdbsN{$clustPDB{$keyC[$j]}[$jj]})";
  }
  print "\n";
  
  $i++;
  $nmuts+=$nmutclass[$j];
  } 
}

print "Total $i $nmuts\n"; 



