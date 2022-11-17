#!/usr/bin/perl
use strict;

# usage  tables.pl > table

my @list_th=("Cartddg","FoldX");

my @list=("Evo","PopMs","Dynamut2","DDGun3D","TherNet","ACDCNN");

my @listK=("KORPM");

my $cmd ="";


my $bin ="100";

my $ws="1.875  0.718  0.876  0.876  2.440  0.973  1.030  2.512  1.030  2.512  1.643  1.113  1.411  1.113  1.608  1.113  1.113  2.512  1.411  2.440"; # B 3    5
my $ws="1.531  0.659  0.730  0.730  2.202  0.735  0.883  2.232  0.883  2.232  1.188  0.961  1.348  0.961  1.281  0.961  0.961  2.232  1.348  2.202"; # exluded Ssym  Kd 3    5




system ("../sbg/bin/korpm  ../Ssym.txt -r 5 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPM_Pucci2018N.txt --dexp --rsa --class --weights \"$ws\" > log");
system ("../scripts/Mstat.pl KORPM_Pucci2018N.txt 10 11 2 > NKORPM.Mlog");
system ("../scripts/confusion.pl  KORPM_Pucci2018N.txt 10 11 $bin >  KORPM_Pucci2018N_prc.txt");
system ("../sbg/bin/korpm  ../SsymR.txt -r 5 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPM_Pucci2018revN.txt --dexp --rsa --class --weights \"$ws\" > log");
system ("../scripts/Mstat.pl KORPM_Pucci2018revN.txt 10 11 2 > RKORPM.Mlog");
system ("../scripts/confusion.pl  KORPM_Pucci2018revN.txt 10 11 $bin >  KORPM_Pucci2018revN_prc.txt");

system ("../sbg/bin/korpm  ../SsymD.txt -r 5 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPM_Pucci2018dirN.txt --dexp --rsa --class --weights \"$ws\" > log");
system ("../scripts/Mstat.pl KORPM_Pucci2018dirN.txt 10 11 2 > DKORPM.Mlog");
system ("../scripts/confusion.pl  KORPM_Pucci2018dirN.txt 10 11 $bin >  KORPM_Pucci2018dirN_prc.txt");
system ("paste KORPM_Pucci2018dirN.txt  KORPM_Pucci2018revN.txt > temp"); 
system ("awk \'function abs(x){return (x < 0) ? -x : x;} {printf \"\%4s \%4s \%-7s \%-7s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \\n\",\$1,\$19, \$2, \$20, \$10, \$11,\$29, (\$11+\$29), abs((\$11+\$29)), \$3+\$4  }\' temp > KORPM_Pucci2018_sym.txt");



foreach (@list_th) {
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018N_th8.txt 3 4 2 >  N$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018N_th8.txt 3 4 $bin >  $_\_Pucci2018N_th8_prc.txt"; system ($cmd); 
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018dirN_th8.txt 3 4 2 > D$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018dirN_th8.txt 3 4 $bin >  $_\_Pucci2018dirN_th8_prc.txt"; system ($cmd); 
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018revN_th8.txt 3 4 2 > R$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018revN_th8.txt 3 4 $bin >  $_\_Pucci2018revN_th8_prc.txt"; system ($cmd); 
 $cmd="paste $_\_Pucci2018dirN_th8.txt  $_\_Pucci2018revN_th8.txt  > temp"; system ($cmd); 
 $cmd="awk 'function abs(x){return (x < 0) ? -x : x;} {printf \"%4s %4s %-7s %-7s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\\n\",  \$1,\$5, \$2, \$6, \$3, \$4, \$8, (\$4+\$8), abs((\$4+\$8)), (\$3+\$7)  }' temp > $_\_Pucci2018_sym.txt"; system ($cmd); 

 
 
# print "$cmd\n"; 
# system ($cmd);  
}

foreach (@list) {
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018N.txt 3 4 2 > N$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018N.txt 3 4 $bin >  $_\_Pucci2018N_prc.txt"; system ($cmd);
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018dirN.txt 3 4 2 > D$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018dirN.txt 3 4 $bin >  $_\_Pucci2018dirN_prc.txt"; system ($cmd); 
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018revN.txt 3 4 2 > R$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018revN.txt 3 4 $bin >  $_\_Pucci2018revN_prc.txt"; system ($cmd); 
 $cmd="paste $_\_Pucci2018dirN.txt  $_\_Pucci2018revN.txt  > temp"; system ($cmd); 
 $cmd="awk 'function abs(x){return (x < 0) ? -x : x;} {printf \"%4s %4s %-7s %-7s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\\n\",  \$1,\$5, \$2, \$6, \$3, \$4, \$8, (\$4+\$8), abs((\$4+\$8)), (\$3+\$7)  }' temp > $_\_Pucci2018_sym.txt"; system ($cmd); 

# print "$cmd\n"; 
  
  
}

push(@list_th, @list);

system (' echo "\n SYM" ');
system ('echo "              S     D     T   TP  avg  err   FP   TN  avg  err   FN    P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE   MAE   PCC    Sc    Ob1   Ob2  MCC     AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe"' );



foreach (@listK) {
 $cmd="printf \"%-8s\" $_; grep \"\^ X \" N$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018N_prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
 #print $cmd;
} 
 
foreach (@list_th) {
 $cmd="printf \"%-8s\" $_; grep \"\^ X \" N$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018N*prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}
system ('echo "\n REVERSE"');
system ('echo "              S     D     T   TP  avg  err   FP   TN  avg  err   FN    P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE   MAE   PCC    Sc    Ob1   Ob2  MCC     AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe"' );

foreach (@listK) {
 $cmd="printf \"%-8s\" $_; grep \"\^ X \" R$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018revN_prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}

foreach (@list_th) {
 $cmd="printf \"%-8s\" $_; grep \"\^ X \" R$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018revN*prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}

system ('echo "\n DIRECT " ');
system ('echo "              S     D     T   TP  avg  err   FP   TN  avg  err   FN    P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE   MAE   PCC    Sc    Ob1   Ob2  MCC     AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe"' );

foreach (@listK) {
$cmd="printf \"%-8s\" $_; grep \"\^ X \" D$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018dirN_prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}
foreach (@list_th) {
 $cmd="printf \"%-8s\" $_; grep \"\^ X \" D$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018dirN*prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}





 

  
 
 
 
 
