#!/usr/bin/perl
use strict;

# usage  tables.pl > table


#awk -F "," '{ printf "%s %9s %8.3f %8.3f\n",substr($2,0,4),substr($3,0,1)substr($2,length($2))substr($3,2,length($3)),$7,$12}' Data_s669_with_predictions.csv >  S669_FOLDX_DIR.txt
#awk -F "," '{ printf "%s %9s %8.3f %8.3f\n",substr($2,0,4),substr($3,0,1)substr($2,length($2))substr($3,2,length($3)),$7,$14}' Data_s669_with_predictions.csv >  S669_Dynamut2_DIR.txt
#awk -F "," '{ printf "%s %9s %8.3f %8.3f\n",substr($2,0,4),substr($3,0,1)substr($2,length($2))substr($3,2,length($3)),$7,$22}' Data_s669_with_predictions.csv >  S669_DDgun3D_DIR.txt
#awk -F "," '{ printf "%s %9s %8.3f %8.3f\n",substr($2,0,4),substr($3,0,1)substr($2,length($2))substr($3,2,length($3)),$7,$25}' Data_s669_with_predictions.csv >  S669_ACDCNN_DIR.txt
#awk -F "," '{ printf "%s %9s %8.3f %8.3f\n",substr($2,0,4),substr($3,0,1)substr($2,length($2))substr($3,2,length($3)),$7,$48}' Data_s669_with_predictions.csv >  S669_ThermoNet_DIR.txt
#awk -F "," '{ printf "%s %9s %8.3f %8.3f\n",substr($2,0,4),substr($3,0,1)substr($2,length($2))substr($3,2,length($3)),$7,$42}' Data_s669_with_predictions.csv >  S669_PopMx_DIR.txt

#awk 'function abs(x){return (x < 0) ? -x : x;}   function thd(x){return (x < -8.0) ? -8 : x;}  { printf "%4s %-7s %7.2f %7.2f\n",  $1 , $2, $3, thd($4)  }' Cartesianddg_S669_new.txt > Cartesianddg_S669_new_th8.txt
#awk 'function abs(x){return (x < 0) ? -x : x;}   function thd(x){return (x < -8.0) ? -8 : x;}  { printf "%4s %-7s %7.2f %7.2f\n",  $1 , $2, $3, thd($4)  }' FoldX_S669_new.txt > FoldX_S669_new_th8.txt
#awk 'function abs(x){return (x < 0) ? -x : x;}   function thd(x){return (x < -8.0) ? -8 : x;}  { printf "%4s %-7s %7.2f %7.2f\n",  $1 , $2, $3, thd($4)  }' Evo_S669_new.txt > Evo_S669_new_th8.txt

my @list0=("Dynamut2","PopMs","DDgun3D","ThermoNet","ACDCNN");
my @list=("CartDD","FoldX","EvoFF","Dynamut2","PopMs","DDgun3D","ThermoNet","ACDCNN");
my @listF=("KORPM","CartDD","FoldX","EvoFF","Dynamut2","PopMs","DDgun3D","ThermoNet","ACDCNN");


my $cmd ="";

my $ws="1.875  0.718  0.876  0.876  2.440  0.973  1.030  2.512  1.030  2.512  1.643  1.113  1.411  1.113  1.608  1.113  1.113  2.512  1.411  2.440"; # B standad weights
my $ws="1.889  0.710  0.894  0.894  2.333  0.941  1.047  2.526  1.047  2.526  1.638  1.114  1.400  1.114  1.638  1.114  1.114  2.526  1.400  2.333"; # excluded S461 (from S669)



print "$ws\n"; 

# removed cases
my $dump="2JIE\\|1XZO\\|1N18\\|2VY0\\|1O1U\\|4YEE\\|4YEF\\|3BCI\\|3O39\\|2MPC\\|3FIS\\|2KJ3\\|1FH5\\|3D2A\\|3K82\\|1XWS\\|1IR3\\|5VP3\\|1PFL\\|3S92\\|4BJX\\|3BCI\\|2MPC\\|3S92\\|4BJX\\|2CLR\\|1X0J\\|1R6R\\|1A7V\\|4N6V\\|1PRG\\|1MN1\\|2PR5\\|1GWY\\|2KS4\\|1GLU\\|1HCQ\\|3ECU\\|1PRE\\|1SPD\\|2JUC\\|2DVV\\|2OUO\\|1OSI\\|1F8I\\|1D5G\\|4WAA\\|1DXX\\|2RPN\\|3G1G\\|2BJD";

 $cmd="grep -v \"I129\\|I130\\|$dump\"  Evo_S669_new_th8.txt > S461_EvoFF_DIR.txt"; system ($cmd);
 #print "$cmd\n"; 
 $cmd="grep -v  \"I129\\|I130\\|$dump\" FoldX_S669_new_th8.txt > S461_FoldX_DIR.txt"; system ($cmd);
 #print "$cmd\n"; 
 $cmd="grep -v  \"I129\\|I130\\|$dump\" Cartesianddg_S669_new_th8.txt > S461_CartDD_DIR.txt"; system ($cmd);

foreach (@list0) {
 $cmd="./ch_val.pl  S669_$_\_DIR.txt >  temp"; system ($cmd);
 $cmd="grep -v \"I129\\|I130\\|$dump\" temp  > S461_$_\_DIR.txt";
 #print "$cmd\n";  
 system ($cmd);
}


### KORMP
# change values --> changes.txt 
$cmd="./ch_val.pl  S669.txt >  temp"; system ($cmd);
$cmd="grep -v  \"I129\\|I130\\|$dump\" temp > S461.txt"; system ($cmd);
$cmd="korpm  S461.txt -r 5 --dir ../S461 --score_file ../pot/korp6Dv1.bin -o S461_KORPM_DIR.txt --dexp --rsa --class --weights \"$ws\" > log";system ($cmd);
system ("../scripts/Mstat.pl S461_KORPM_DIR.txt  10 11 2 > DKORPM.Mlog");
system ("../scripts/confusion.pl   S461_KORPM_DIR.txt 10 11 200 >  S461_KORPM_DIR_prc.txt");

#print $cmd;
 


### REST

foreach (@list) {


 $cmd="../scripts/Mstat.pl   S461_$_\_DIR.txt 3 4 2 > D$_.Mlog\n"; system ($cmd);
 # print "$cmd\n";
 $cmd="../scripts/confusion.pl  S461_$_\_DIR.txt 3 4 200 >  S461_$_\_DIR_prc.txt"; system ($cmd);
}




system ('echo "\n               S     D     T   TP  avg  err   FP   TN  avg  err   FN    P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE   MAE   PCC    Sc    Ob1   Ob2  MCC     AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe"' );

foreach (@listF) {
$cmd="printf \"%-9s\" $_; grep \"\^ X \" D$_.Mlog | sed  \'s\/\X\/\/\' | tr -d \'\\n\'; echo -n \" \"; tail -1 S461_$_\_DIR_prc.txt | sed  \'s\/\#\/\/\' "; system ($cmd);
}





 

  
 
 
 
 
