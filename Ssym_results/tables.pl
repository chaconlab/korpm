#!/usr/bin/perl
use strict;

# usage  tables.pl

my @list_th=("Cartddg","FoldX");

my @list=("Evo","PopMs","Dynamut","DDGun3D","TherNet");

my $cmd ="";

system ("../../Mstat.pl KORPM_Pucci2018N.txt 10 11 2 > NKORPM.Mlog");
system ("../../confusion.pl  KORPM_Pucci2018N.txt 10 11 1000 >  KORPM_Pucci2018N_prc.txt");


system ("../../Mstat.pl KORPM_Pucci2018revN.txt 10 11 2 > RKORPM.Mlog");
system ("../../confusion.pl  KORPM_Pucci2018revN.txt 10 11 1000 >  KORPM_Pucci2018revN_prc.txt");


system ("../../Mstat.pl KORPM_Pucci2018dirN.txt 10 11 2 > DKORPM.Mlog");
system ("../../confusion.pl  KORPM_Pucci2018dirN.txt 10 11 1000 >  KORPM_Pucci2018dirN_prc.txt");

   

foreach (@list_th) {
 $cmd="../../Mstat.pl  $_\_Pucci2018N_th8.txt 3 4 2 >  N$_.Mlog\n"; system ($cmd);
 $cmd="../../confusion.pl  $_\_Pucci2018N_th8.txt 3 4 1000 >  $_\_Pucci2018N_th8_prc.txt"; system ($cmd); 
 $cmd="../../Mstat.pl  $_\_Pucci2018dirN_th8.txt 3 4 2 > D$_.Mlog\n"; system ($cmd);
 $cmd="../../confusion.pl  $_\_Pucci2018dirN_th8.txt 3 4 1000 >  $_\_Pucci2018dirN_th8_prc.txt"; system ($cmd); 
 $cmd="../../Mstat.pl  $_\_Pucci2018revN_th8.txt 3 4 2 > R$_.Mlog\n"; system ($cmd);
 $cmd="../../confusion.pl  $_\_Pucci2018revN_th8.txt 3 4 1000 >  $_\_Pucci2018revN_th8_prc.txt"; system ($cmd); 
 $cmd="paste $_\_Pucci2018dirN_th8.txt  $_\_Pucci2018revN_th8.txt  > temp"; system ($cmd); 
 $cmd="awk 'function abs(x){return (x < 0) ? -x : x;} {printf \"%4s %4s %-7s %-7s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\\n\",  \$1,\$5, \$2, \$6, \$3, \$4, \$8, (\$4-\$8), abs((\$4-\$8)), (\$3+\$7)  }' temp > $_\_Pucci2018_th8_sym.txt"; system ($cmd); 

 
 
# print "$cmd\n"; 
  # system ($cmd);  
}

foreach (@list) {
 $cmd="../../Mstat.pl  $_\_Pucci2018N.txt 3 4 2 > N$_.Mlog\n"; system ($cmd);
 $cmd="../../confusion.pl  $_\_Pucci2018N.txt 3 4 1000 >  $_\_Pucci2018N_prc.txt"; system ($cmd);
 $cmd="../../Mstat.pl  $_\_Pucci2018dirN.txt 3 4 2 > D$_.Mlog\n"; system ($cmd);
 $cmd="../../confusion.pl  $_\_Pucci2018dirN.txt 3 4 1000 >  $_\_Pucci2018dirN_prc.txt"; system ($cmd); 
 $cmd="../../Mstat.pl  $_\_Pucci2018revN.txt 3 4 2 > R$_.Mlog\n"; system ($cmd);
 $cmd="../../confusion.pl  $_\_Pucci2018revN.txt 3 4 1000 >  $_\_Pucci2018revN_prc.txt"; system ($cmd); 
 $cmd="paste $_\_Pucci2018dirN.txt  $_\_Pucci2018revN.txt  > temp"; system ($cmd); 
 $cmd="awk 'function abs(x){return (x < 0) ? -x : x;} {printf \"%4s %4s %-7s %-7s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\\n\",  \$1,\$5, \$2, \$6, \$3, \$4, \$8, (\$4-\$8), abs((\$4-\$8)), (\$3+\$7)  }' temp > $_\_Pucci2018_sym.txt"; system ($cmd); 

# print "$cmd\n"; 
  
  
}

push(@list_th, @list);

system (' echo "SYM" ');
system ('echo "              #     S     D     T   TP  avg  err   FP   TN  avg  err   FN   NC    P     N    TPR   FPR   SPE   PPV   NPV   ACC   ERR  accn  RMSE   MAE    PCC    Sc   Ob1   Ob2   MCC #   AUC_R   AUC_P   THROC    TPR     FPR     BMCC    TH      TPR     FPR"' );

my $K="KORPM";
 $cmd="printf \"%-7s\" $K; grep \"\^ X \" N$K.Mlog | tr -d \'\\n\'; echo -n \" \"; tail -1 $K\_Pucci2018N_prc.txt"; system ($cmd);
foreach (@list_th) {
 $cmd="printf \"%-7s\" $_; grep \"\^ X \" N$_.Mlog | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018N*prc.txt"; system ($cmd);
}
system ('echo "REVERSE"');
system ('echo "              #     S     D     T   TP  avg  err   FP   TN  avg  err   FN   NC    P     N    TPR   FPR   SPE   PPV   NPV   ACC   ERR  accn  RMSE   MAE    PCC    Sc   Ob1   Ob2   MCC #   AUC_R   AUC_P   THROC    TPR     FPR     BMCC    TH      TPR     FPR"' );

$cmd="printf \"%-7s\" $K; grep \"\^ X \" R$K.Mlog | tr -d \'\\n\'; echo -n \" \"; tail -1 $K\_Pucci2018revN_prc.txt"; system ($cmd);
foreach (@list_th) {
 $cmd="printf \"%-7s\" $_; grep \"\^ X \" R$_.Mlog | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018revN*prc.txt"; system ($cmd);
}
system ('echo "DIRECT " ');
system ('echo "              #     S     D     T   TP  avg  err   FP   TN  avg  err   FN   NC    P     N    TPR   FPR   SPE   PPV   NPV   ACC   ERR  accn  RMSE   MAE    PCC    Sc   Ob1   Ob2   MCC #   AUC_R   AUC_P   THROC    TPR     FPR     BMCC    TH      TPR     FPR"' );

$cmd="printf \"%-7s\" $K; grep \"\^ X \" D$K.Mlog | tr -d \'\\n\'; echo -n \" \"; tail -1 $K\_Pucci2018dirN_prc.txt"; system ($cmd);
foreach (@list_th) {
 $cmd="printf \"%-7s\" $_; grep \"\^ X \" D$_.Mlog | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018dirN*prc.txt"; system ($cmd);
}





 

  
 
 
 
 
