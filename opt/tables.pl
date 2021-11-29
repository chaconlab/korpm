#!/usr/bin/perl
use strict;

# usage  tables.pl



my @Kfold=("5","10");
my @kweight=("5","3");
my $cmd ="";
printf "K-fold cluster by homology  Kw=5->20 weights   Kw=3->REDUCED ALPHABET DE+ILV+NQST+FY+KH\n"; 
printf "   Kw  Kf      A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y  |  RMSE   MAE   PPV   Sen   Spe   PCC   ACC   NPV   MCC   ROC   PRC   AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH      Sen     Spe\n";



foreach my $kw (@kweight) {
foreach my $k (@Kfold) {
my $out="r".$kw."_k".$k.".log";

# optimize
$cmd="korpm_opt Id50c08_1merNCLB.opt -r $kw --rsa --class -k $k -n 100 > $out;";
system ($cmd); 

# get prc
my $out2="r".$kw."_k".$k.".prc";
$cmd="../../confusion.pl  out.txt 2 3 100 >  $out2\n"; 
system ($cmd);
#print $cmd;
# print
$cmd="printf \"%4d %4d \" $kw $k;  grep AVG $out \| sed  \'s\/AVG\/\/\' \| tr -d \'\\n\';  tail -1 $out2 \| sed  \'s\/\#\/\/\' \n";
#print $cmd;
system ($cmd);
$cmd="printf \"%4d %4d \" $kw $k;  grep SIG $out \| sed  \'s\/SIG\/\/\';";  
system ($cmd);


#print $cmd;

}
}
printf "\nK-fold pure random Kw=5->20 weights   Kw=3->REDUCED ALPHABET DE+ILV+NQST+FY+KH\n\n"; 
printf "   Kw  Kf      A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y  |  RMSE   MAE   PPV   Sen   Spe   PCC   ACC   NPV   MCC   ROC   PRC   AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH      Sen     Spe\n";
foreach my $kw (@kweight) {
foreach my $k (@Kfold) {
my $out="r".$kw."_k".$k.".log";

# optimize
$cmd="\nkorpm_opt Id50c08_1merNCLB.opt -r $kw --rsa --class -k $k -n 100 --krand > $out;\n\n";
system ($cmd); 

# get prc
my $out2="r".$kw."_k".$k.".prc";
$cmd="../../confusion.pl  out.txt 2 3 100 >  $out2\n"; 
system ($cmd);
#print $cmd;
# print
$cmd="printf \"%4d %4d \" $kw $k;  grep AVG $out \| sed  \'s\/AVG\/\/\' \| tr -d \'\\n\';  tail -1 $out2 \| sed  \'s\/\#\/\/\' \n";
#print $cmd;
system ($cmd);
$cmd="printf \"%4d %4d \" $kw $k;  grep SIG $out \| sed  \'s\/SIG\/\/\';";  
system ($cmd);

}
}









 

  
 
 
 
 
