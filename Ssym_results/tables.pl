#!/usr/bin/perl
use strict;

# usage  tables.pl > table

my @list_th=("Cartddg","FoldX");

my @list=("Evo","PopMs","Dynamut","DDGun3D","TherNet","ACDCNN");

my @listK=("KORPM");

my $cmd ="";

#system ("korpm  ../Pucci2018N.txt -r 5 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPM_Pucci2018N.txt --dexp --rsa --class  > log");
#system ("../scripts/Mstat.pl KORPM_Pucci2018N.txt 10 11 2 > NKORPM.Mlog");
#system ("../scripts/confusion.pl  KORPM_Pucci2018N.txt 10 11 1000 >  KORPM_Pucci2018N_prc.txt");
#system ("korpm  ../Pucci2018revN.txt -r 5 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPM_Pucci2018revN.txt --dexp --rsa --class > log");
#system ("../scripts/Mstat.pl KORPM_Pucci2018revN.txt 10 11 2 > RKORPM.Mlog");
#system ("../scripts/confusion.pl  KORPM_Pucci2018revN.txt 10 11 1000 >  KORPM_Pucci2018revN_prc.txt");
#system ("korpm  ../Pucci2018dirN.txt -r 5 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPM_Pucci2018dirN.txt --dexp --rsa --class > log");
#system ("../scripts/Mstat.pl KORPM_Pucci2018dirN.txt 10 11 2 > DKORPM.Mlog");
#system ("../scripts/confusion.pl  KORPM_Pucci2018dirN.txt 10 11 1000 >  KORPM_Pucci2018dirN_prc.txt");
#system ("paste KORPM_Pucci2018dirN.txt  KORPM_Pucci2018revN.txt > temp"); 
#system ("awk \'function abs(x){return (x < 0) ? -x : x;} {printf \"\%4s \%4s \%-7s \%-7s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \\n\",\$1,\$19, \$2, \$20, \$10, \$11,\$29, (\$11+\$29), abs((\$11+\$29)), \$3+\$4  }\' temp > KORPM_Pucci2018_sym.txt");


#system ("korpm  ../Pucci2018N.txt -r 21 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPMt_Pucci2018N.txt --dexp --rsa --class > log");
#system ("../scripts/Mstat.pl KORPMt_Pucci2018N.txt 10 11 2 > NKORPMt.Mlog");
#system ("../scripts/confusion.pl  KORPMt_Pucci2018N.txt 10 11 1000 >  KORPMt_Pucci2018N_prc.txt");
#system ("korpm  ../Pucci2018revN.txt -r 21 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPMt_Pucci2018revN.txt --dexp --rsa --class > log");
#system ("../scripts/Mstat.pl KORPMt_Pucci2018revN.txt 10 11 2 > RKORPMt.Mlog");
#system ("../scripts/confusion.pl  KORPMt_Pucci2018revN.txt 10 11 1000 >  KORPMt_Pucci2018revN_prc.txt");
#system ("korpm  ../Pucci2018dirN.txt -r 21 --dir ../Ssym --score_file ../pot/korp6Dv1.bin -o KORPMt_Pucci2018dirN.txt --dexp --rsa --class > log");
#system ("../scripts/Mstat.pl KORPMt_Pucci2018dirN.txt 10 11 2 > DKORPMt.Mlog");
#system ("../scripts/confusion.pl  KORPMt_Pucci2018dirN.txt 10 11 1000 >  KORPMt_Pucci2018dirN_prc.txt");
#system ("paste KORPMt_Pucci2018dirN.txt  KORPMt_Pucci2018revN.txt > temp"); 
#system ("awk \'function abs(x){return (x < 0) ? -x : x;} {printf \"\%4s \%4s \%-7s \%-7s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \\n\",\$1,\$19, \$2, \$20, \$10, \$11,\$29, (\$11+\$29), abs((\$11+\$29)), \$3+\$4  }\' temp > KORPMt_Pucci2018_sym.txt");
   

foreach (@list_th) {
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018N_th8.txt 3 4 2 >  N$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018N_th8.txt 3 4 1000 >  $_\_Pucci2018N_th8_prc.txt"; system ($cmd); 
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018dirN_th8.txt 3 4 2 > D$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018dirN_th8.txt 3 4 1000 >  $_\_Pucci2018dirN_th8_prc.txt"; system ($cmd); 
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018revN_th8.txt 3 4 2 > R$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018revN_th8.txt 3 4 1000 >  $_\_Pucci2018revN_th8_prc.txt"; system ($cmd); 
 $cmd="paste $_\_Pucci2018dirN_th8.txt  $_\_Pucci2018revN_th8.txt  > temp"; system ($cmd); 
 $cmd="awk 'function abs(x){return (x < 0) ? -x : x;} {printf \"%4s %4s %-7s %-7s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\\n\",  \$1,\$5, \$2, \$6, \$3, \$4, \$8, (\$4+\$8), abs((\$4+\$8)), (\$3+\$7)  }' temp > $_\_Pucci2018_th8_sym.txt"; system ($cmd); 

 
 
# print "$cmd\n"; 
# system ($cmd);  
}

foreach (@list) {
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018N.txt 3 4 2 > N$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018N.txt 3 4 1000 >  $_\_Pucci2018N_prc.txt"; system ($cmd);
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018dirN.txt 3 4 2 > D$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018dirN.txt 3 4 1000 >  $_\_Pucci2018dirN_prc.txt"; system ($cmd); 
 $cmd="../scripts/Mstat.pl  $_\_Pucci2018revN.txt 3 4 2 > R$_.Mlog\n"; system ($cmd);
 $cmd="../scripts/confusion.pl  $_\_Pucci2018revN.txt 3 4 1000 >  $_\_Pucci2018revN_prc.txt"; system ($cmd); 
 $cmd="paste $_\_Pucci2018dirN.txt  $_\_Pucci2018revN.txt  > temp"; system ($cmd); 
 $cmd="awk 'function abs(x){return (x < 0) ? -x : x;} {printf \"%4s %4s %-7s %-7s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\\n\",  \$1,\$5, \$2, \$6, \$3, \$4, \$8, (\$4+\$8), abs((\$4+\$8)), (\$3+\$7)  }' temp > $_\_Pucci2018_sym.txt"; system ($cmd); 

# print "$cmd\n"; 
  
  
}

push(@list_th, @list);

system (' echo "SYM" ');
system ('echo "             S     D     T   TP  avg  err   FP   TN  avg  err   FN    P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE   MAE   PCC    Sc    Ob1   Ob2  MCC     AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe"' );





foreach (@listK) {
 $cmd="printf \"%-7s\" $_; grep \"\^ X \" N$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018N_prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
 #print $cmd;
} 
 
foreach (@list_th) {
 $cmd="printf \"%-7s\" $_; grep \"\^ X \" N$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018N*prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}
system ('echo "REVERSE"');
system ('echo "             S     D     T   TP  avg  err   FP   TN  avg  err   FN    P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE   MAE   PCC    Sc    Ob1   Ob2  MCC     AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe"' );

foreach (@listK) {
 $cmd="printf \"%-7s\" $_; grep \"\^ X \" R$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018revN_prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}

foreach (@list_th) {
 $cmd="printf \"%-7s\" $_; grep \"\^ X \" R$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018revN*prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}

system ('echo "DIRECT " ');
system ('echo "DIRECT " ');
system ('echo "             S     D     T   TP  avg  err   FP   TN  avg  err   FN    P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE   MAE   PCC    Sc    Ob1   Ob2  MCC     AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe"' );

foreach (@listK) {
$cmd="printf \"%-7s\" $_; grep \"\^ X \" D$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018dirN_prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}
foreach (@list_th) {
 $cmd="printf \"%-7s\" $_; grep \"\^ X \" D$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 $_\_Pucci2018dirN*prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}





 

  
 
 
 
 
