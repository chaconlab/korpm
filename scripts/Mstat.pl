#!/usr/bin/perl
use strict;
use Cwd qw();

if ( not ($#ARGV == 3 or $#ARGV == 4 or $#ARGV == 5 or $#ARGV == 6 or $#ARGV == 7) ) 
{ 
  print "USAGE:\n\t$0 <data.txt> [col_exp] [col_pred] [col_mut] [thr](optional)\n";
  print "\nSINTAX:\n";
  print "\t<data.txt> --> Input text file with data\n";
  print "\t[col_exp]  --> Column index for Experimental truth (>0 positives and <=0 negatives).\n";
  print "\t[col_pred] --> Column index for Predicted results (>0 positives and <=0 negatives).\n";
  print "\t[col_mut]  --> Column index for mutation code (e.g. CA30S).\n";
  print "\t[thr]      --> (optional) Threshold (ddG >= thr are not-destabilizing [positives] and ddG < thr are destabilizing [negatives]). By default it is 0. Absolute value of the predicted data will be used instead of raw data if the \"abs\" string is appended to the thr value (e.g. 0.89abs) (WARNING: the Positive Experimental values must be greater than this threshold).\n";
  print "\t[thr3stat] --> (optional) Threshold for 3-states stuff (ddG >= thr3stat is Stabilizing, ddG <= thr3stat is Destabilizing, and ddG between +thr3stat and -thr3stat are Neutral.\n";
  print "\t[rsaL] -->threshold for remove rsa < rsaL\n";
  print "\t[rsaG] -->threshold for remove rsa < rsaG\n";
 
  print "\nDESCRIPTION:\n";
  print "\tCompute several Binary clasifier statistics (confusion matrix stuff).\n";
  print "\nEXAMPLE:\n";
  print "\tKORPM format --> Mstat.pl out_all.txt 10 11 2\n";
  print "\t3-cols format --> Mstat.pl resultN.txt 3 4 2\n";
  print "\tUsing custom threshold --> Mstat.pl out_all.txt 10 11 2 -1.0\n";
  print "\tUsing absolute value threshold --> Mstat.pl out_all.txt 10 11 2 0.9abs\n";
  
  print "\n";
  exit
}

# Get parser input
my $file = $ARGV[0]; # Input text file with data
my $col_exp = int($ARGV[1]) - 1; # Column index (0,1,...) for Experimental truth (>0 positives and <=0 negatives)
my $col_pred = int($ARGV[2]) - 1; # Column index (0,1,...) for Predicted results (>0 positives and <=0 negatives)
my $col_mut = int($ARGV[3]) - 1; # Column index (0,1,...) for Experimental truth (>0 positives and <=0 negatives)
my $thr = 0.0; # (optional) Threshold (ddG >= thr are not-destabilizing [positives] and ddG < thr are destabilizing [negatives]).\n";
my $thrE = 0.0; # Current threshold for Experimental Positives or Negatives
my $thr3stat = 1.0; # 3-states threshold (default= 1.0, as in Rosetta paper)
my $rsaL = -10000.0; # threshold for remove rsa < $rsaL
my $rsaG = +10000.0;   # threshold for remove rsa > $rsaL
 $thr = $ARGV[4] if($ARGV[4]);


my $doabs = 0;
if($thr =~ /abs/)
{
  $doabs = 1;
  $thr =~ s/abs//ig; # Remove "abs" tag
  $thr *= 1.0; # Force float number
  print "# Current Positives|Negatives threshold (thr) is $thr (ddG >= +$thr and ddG <= -$thr are Positives and ddG between +$thr and -$thr are Negatives.\n";
}
else
{
  print "# Current Positives|Negatives threshold (thr) is $thr (ddG >= $thr are not-destabilizing [positives] and ddG < $thr are destabilizing [negatives]).\n";
}
$thr3stat = $ARGV[5] if($ARGV[5]);

$rsaL = $ARGV[6] if($ARGV[6]);
$rsaG = $ARGV[7] if($ARGV[7]);

my $thr_stab = $thr3stat; # Lower threshold to consider a mutation as Stabilizing (stabilizing >= 1), see: Rosetta's 768 mutants paper
my $thr_des = -$thr3stat; # Upper threshold to consider a mutation as Destabilizing (destabilizing <= 1)

# Some variables
my @f;
my @seq;
my @ch1;
my @ch2;
my @dde;
my @ddc;
my @classExp; # Experimental classification: 0= Stabilizing, 1= Neutral, 2= Destabilizing
my @classCalc; # Calculated (prediction) classification: 0= Stabilizing, 1= Neutral, 2= Destabilizing
my $i;
my $rmse;
# Performance statistics
my $TPddG = 0; # True Positives for stabilizing ddG prediction
my $FPddG = 0; # False Positives for stabilizing ddG prediction
my $TNddG = 0; # True Negatives for stabilizing ddG prediction
my $FNddG = 0; # False Negatives for stabilizing ddG prediction
my $NCddG = 0; # Non-Conclusive
my $TPavg = 0.0; # True Positives ddG average
my $TPerr = 0.0; # True Positives ddG error average
my $TNavg = 0.0; # True Negatives ddG average
my $TNerr = 0.0; # True Negatives ddG error average
my $nSc = 0; # Number of "Same class"
my $nOb1 = 0; # Number of "Off by one" class
my $nOb2 = 0; # Number of "Off by two" class
my $cont = 0;
  
$i=0; 
$rmse=0;

open(MYFILE, $file) or die;
while(<MYFILE>) 
{
  next if($_ =~ /^#/);
  
  @f = split(/\s+/, $_);
  
  $dde[$i] = $f[$col_exp] * 1.0;
#  $dde[$i] += 5 if($dde[$i] > 0);
#  $dde[$i] += 5 if($dde[$i] > 0);
#  $dde[$i] *= -1;

   
  if($doabs == 1)
  {
    $ddc[$i] = abs( $f[$col_pred] * 1.0 );
  }
  else
  {
    $ddc[$i] = $f[$col_pred] * 1.0;
  }

  
# rsa 
    if ( ($f[17]*1.0) < $rsaL) {next};
    if ( ($f[17]*1.0) > $rsaG) {next};

#  if (($f[17] > 0.7) or ($f[17] < 0.3)) {next};
#  print "$f[0] $f[$col_mut] $f[$col_exp]  $f[17]\n";
#  $ddc[$i] += 0.3;
 # if ( $ddc[$i] < 0) {$ddc[$i] += 0.1}
 # if ( $ddc[$i] > 0) {$ddc[$i] -= 0.1}

  $seq[$i] = $f[$col_mut];
  $ch1[$i] = substr( $seq[$i], 0, 1 ); # WT aminoacid
  $ch2[$i] = substr( $seq[$i], -1, 1 ); # Mut aminoacid
 
 if ( $ddc[$i] > 0.0 )
 {  
  # $ddc[$i] -= 0.3;
 } else 
 {
  # $ddc[$i] += 0.3;
 }

       
 if ((( $ch1[$i] eq 'A' ) or ( $ch1[$i] eq 'G' ) or ( $ch1[$i] eq 'S' ) or ( $ch1[$i] eq 'C' )) and (( $ch2[$i] eq 'I' ) or ( $ch2[$i] eq 'M' ) or ( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'K' ) or ( $ch2[$i] eq 'L' ) or ( $ch2[$i] eq 'F' ) or ( $ch2[$i] eq 'W' ) or ( $ch2[$i] eq 'Y' )) )
        {
 #        $ddc[$i] += 0.2;
          
        }
# bs
 if ((( $ch2[$i] eq 'A' ) or ( $ch2[$i] eq 'G' ) or ( $ch2[$i] eq 'S' ) or ( $ch2[$i] eq 'C' )) and (( $ch1[$i] eq 'I' ) or ( $ch1[$i] eq 'M' ) or ( $ch1[$i] eq 'R' ) or ( $ch1[$i] eq 'K' ) or ( $ch1[$i] eq 'L' ) or ( $ch1[$i] eq 'F' ) or ( $ch1[$i] eq 'W' ) or ( $ch1[$i] eq 'Y' )) )
        {
     #     $ddc[$i] -= 0.2;
       }

# //hw
      
 if ((( $ch1[$i] eq 'A' ) or ( $ch1[$i] eq 'C' ) or ( $ch1[$i] eq 'I' ) or ( $ch1[$i] eq 'L' ) or ( $ch1[$i] eq 'M' ) or ( $ch1[$i] eq 'F' ) or ( $ch1[$i] eq 'W' ) or ( $ch1[$i] eq 'V' ))
        and  (( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'N' ) or ( $ch2[$i] eq 'D' ) or ( $ch2[$i] eq 'E' ) or ( $ch2[$i] eq 'K' ) ) )
       {
#          $ddc[$i] += 0.2;
       }     


# //hn
             
 if ((( $ch1[$i] eq 'A' ) or ( $ch1[$i] eq 'C' ) or ( $ch1[$i] eq 'I' ) or ( $ch1[$i] eq 'L' ) or ( $ch1[$i] eq 'M' ) or ( $ch1[$i] eq 'F' ) or ( $ch1[$i] eq 'W' ) or ( $ch1[$i] eq 'V' ))
        and (( $ch2[$i] eq 'G' ) or ( $ch2[$i] eq 'H' ) or ( $ch2[$i] eq 'P' ) or ( $ch2[$i] eq 'S' ) or ( $ch2[$i] eq 'T' )  or ( $ch2[$i] eq 'Y' )))
        {
#          $ddc[$i] += 0.2;
        }
          
      if ((( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'H' ) or ( $ch2[$i] eq 'K' ) ) and (( $ch1[$i] eq 'D' ) or ( $ch1[$i] eq 'E' ) ) )
        {
      # $ddc[$i] += 0.2;
        } 
        
        if ((( $ch1[$i] eq 'R' ) or ( $ch1[$i] eq 'H' ) or ( $ch1[$i] eq 'K' ) ) and (( $ch2[$i] eq 'D' ) or ( $ch2[$i] eq 'E' ) ) )
        {
        
      # $ddc[$i] += 0.2;

        } 
        
# if($ch1[$i] eq 'R')  {$ddc[$i] += .5;}
# if($ch1[$i] eq 'E')  {$ddc[$i] -= .2}
# if($ch2[$i] eq 'H')  {$ddc[$i] += .2}
#  if($ch2[$i] eq 'R')  {$ddc[$i] -= .4}
#  if($ch2[$i] eq 'K')  {$ddc[$i] -= .2}

 # if (($ch1[$i] eq 'A') and ($ch2[$i] eq 'D') or ($ch2[$i] eq 'E')) {$ddc[$i] *= -1;}
 # if (($ch2[$i] eq 'A') and ($ch1[$i] eq 'D') or ($ch1[$i] eq 'E')) {$ddc[$i] *= 0.5;} 

#     next if($ch1[$i] eq 'W'); 
#     next if($ch2[$i] eq 'W'); # quitar 
#     next if($ch1[$i] eq 'F'); # quitar 
#      next if($ch2[$i] eq 'F');  
#      next if($ch1[$i] eq 'Q'); #  
#     next if($ch2[$i] eq 'Q'); #  quitar 
#     next if($ch1[$i] eq 'Y'); # quitar
#     next if($ch2[$i] eq 'Y'); # 
#     next if($ch1[$i] eq 'P'); # quitar
#   next if($ch2[$i] eq 'P'); # quitar
#    next if($ch1[$i] eq 'R');
#    next if($ch2[$i] eq 'R');
#    next if($ch1[$i] eq 'E'); # *
#    next if($ch2[$i] eq 'E');
#    next if($ch1[$i] eq 'D'); # *
#    next if($ch2[$i] eq 'D');
#    next if($ch1[$i] eq 'H'); # quitar
#    next if($ch2[$i] eq 'H');
#    next if($ch1[$i] eq 'C'); # quitar
#    next if($ch2[$i] eq 'C'); # quitar 


  # Classification of Experiemental data
  if($dde[$i] >= $thr_stab) { $classExp[$i] = 0; } # Stabilizing mutation
  elsif($dde[$i] <= $thr_des) { $classExp[$i] = 2; } # Destabilizing mutation
  else { $classExp[$i] = 1; } # Neutral mutation

  # Classification of Calculated (Predicted) data
  if($ddc[$i] >= $thr_stab) { $classCalc[$i] = 0; } # Stabilizing mutation
  elsif($ddc[$i] <= $thr_des) { $classCalc[$i] = 2; } # Destabilizing mutation
  else { $classCalc[$i] = 1; } # Neutral mutation

  $i++;
 # print "$f[0] $f[$col_mut] $f[$col_exp] $f[$col_pred] $f[16]\n";
 # print "$f[0] $f[$col_mut] $f[$col_exp]  $f[16]\n";
}
close(MYFILE);

my $ti = $i;
my $cont = 0;
# my @aa = ("X","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","p");
my @aa = ("X","A","V","I","L","M","F","W","Y","R","H","K","D","E","S","T","N","Q","C","G","P","p", "1", "2", "3","4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18");
print "\naa     #     S     D     T   TP  avg  err   FP   TN  avg  err   FN   NC    P     N    TPR   FPR   SPE   PPV   NPV   ACC   ERR  accn  RMSE   MAE    PCC    Sc   Ob1   Ob2   MCC\n";
foreach  my $key( @aa ) 
{

  $TPddG = 0; # True Positives for stabilizing ddG prediction
  $FPddG = 0; # False Positives for stabilizing ddG prediction
  $TNddG = 0; # True Negatives for stabilizing ddG prediction
  $FNddG = 0; # False Negatives for stabilizing ddG prediction
  $NCddG = 0; # Non-Conclusive
  $TPavg = 0.0; # True Positives ddG average
  $TPerr = 0.0; # True Positives ddG error average
  $TNavg = 0.0; # True Negatives ddG average
  $TNerr = 0.0; # True Negatives ddG error average
  
  $nSc = 0; # Number of "Same class"
  $nOb1 = 0; # Number of "Off by one" class
  $nOb2 = 0; # Number of "Off by two" class
  my $cont = 0;
  my $rmse = 0; # Root Mean Squared Error
  my $MAE = 0.0; # Mean Absolute Error
  my @exp = (); # Experimental energies
  my @calc = (); # Calculated energies
  
  my $contzero = 0; 
  for (my $i=0; $i < $ti ; $i++) 
  {
  # if( ($key eq 'X') or ($key eq $ch1[$i]) ) 
  #  if( ($key eq 'X') or ($key eq $ch2[$i]) ) 
    if (($key eq 'X') or ( $key > 0) or ($key eq $ch1[$i]) or ($key eq $ch2[$i])) 
    {



      if ( $key eq '1')   {   # small to big   very small A, G, S Small [108-117]: N, D, C, P, T. Medium [138-154]: Q, E, H, V. Large [162-174]: R, I, L, K, M. Very large [189-228]: 

        #  A, G, S vs  F, W, Y. R, I, L, K, M.
        if ((( $ch1[$i] eq 'A' ) or ( $ch1[$i] eq 'G' ) or ( $ch1[$i] eq 'S' ) or ( $ch1[$i] eq 'C' )) and (( $ch2[$i] eq 'I' ) or ( $ch2[$i] eq 'M' ) or ( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'K' ) or ( $ch2[$i] eq 'L' ) or ( $ch2[$i] eq 'F' ) or ( $ch2[$i] eq 'W' ) or ( $ch2[$i] eq 'Y' )) )
        {
        #  $ddc[$i] += 0.2 ;
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}       
       
      } elsif ( $key eq '2')  # big to small
       {
       
        if ((( $ch2[$i] eq 'A' ) or ( $ch2[$i] eq 'G' ) or ( $ch2[$i] eq 'S' ) or ( $ch2[$i] eq 'C' )) and (( $ch1[$i] eq 'I' ) or ( $ch1[$i] eq 'M' ) or ( $ch1[$i] eq 'R' ) or ( $ch1[$i] eq 'K' ) or ( $ch1[$i] eq 'L' ) or ( $ch1[$i] eq 'F' ) or ( $ch1[$i] eq 'W' ) or ( $ch1[$i] eq 'Y' )) )
        {
         # $ddc[$i] += 0.3;
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {
                    #$ddc[$i] -= 0.3;
        next;
        }      
      
      } elsif ( $key eq '3')   # positve nega
       {
       
        if ((( $ch1[$i] eq 'R' ) or ( $ch1[$i] eq 'H' ) or ( $ch1[$i] eq 'K' ) ) and (( $ch2[$i] eq 'D' ) or ( $ch2[$i] eq 'E' ) ) )
        {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;} 
     } elsif ( $key eq '4')  # nega positve
       {
       
        if ((( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'H' ) or ( $ch2[$i] eq 'K' ) ) and (( $ch1[$i] eq 'D' ) or ( $ch1[$i] eq 'E' ) ) )
        {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;} 
     } elsif ( $key eq '5')    # w to h  Hydrophobic    (A, C, I, L, M, F, W, V) ro  Hydrophilic (R, N, D, Q, E, K).  Neutral (G, H, P, S, T, Y)
     
       {
       
        if ((( $ch2[$i] eq 'A' ) or ( $ch2[$i] eq 'C' ) or ( $ch2[$i] eq 'I' ) or ( $ch2[$i] eq 'L' ) or ( $ch2[$i] eq 'M' ) or ( $ch2[$i] eq 'F' ) or ( $ch2[$i] eq 'W' ) or ( $ch2[$i] eq 'V' ))
        and (( $ch1[$i] eq 'R' ) or ( $ch1[$i] eq 'N' ) or ( $ch1[$i] eq 'D' ) or ( $ch1[$i] eq 'E' ) or ( $ch1[$i] eq 'K' ) ) )
        {
          #$ddc[$i] -= 0.0;
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      

      } elsif ( $key eq '6')  # h 2 w 
       {
       
         if ((( $ch1[$i] eq 'A' ) or ( $ch1[$i] eq 'C' ) or ( $ch1[$i] eq 'I' ) or ( $ch1[$i] eq 'L' ) or ( $ch1[$i] eq 'M' ) or ( $ch1[$i] eq 'F' ) or ( $ch1[$i] eq 'W' ) or ( $ch1[$i] eq 'V' ))
        and  (( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'N' ) or ( $ch2[$i] eq 'D' ) or ( $ch2[$i] eq 'E' ) or ( $ch2[$i] eq 'K' ) ) )
       {
          #$ddc[$i] += 0.3;
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
    
       } elsif ( $key eq '7')    # n to h Neutral (G, H, P, S, T, Y) to Hydrophobic (A, C, I, L, M, F, W, V) 
     
       {
       
        if ((( $ch2[$i] eq 'A' ) or ( $ch2[$i] eq 'C' ) or ( $ch2[$i] eq 'I' ) or ( $ch2[$i] eq 'L' ) or ( $ch2[$i] eq 'M' ) or ( $ch2[$i] eq 'F' ) or ( $ch2[$i] eq 'W' ) or ( $ch2[$i] eq 'V' ))
        and (( $ch1[$i] eq 'G' ) or ( $ch1[$i] eq 'H' ) or ( $ch1[$i] eq 'P' ) or ( $ch1[$i] eq 'S' ) or ( $ch1[$i] eq 'T' )  or ( $ch1[$i] eq 'Y' )))
        {

          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}    
          
       } elsif ( $key eq '8')    #  h to n  Hydrophobic (A, C, I, L, M, F, W, V) Neutral (G, H, P, S, T, Y) to
     
       {

       
        if ((( $ch1[$i] eq 'A' ) or ( $ch1[$i] eq 'C' ) or ( $ch1[$i] eq 'I' ) or ( $ch1[$i] eq 'L' ) or ( $ch1[$i] eq 'M' ) or ( $ch1[$i] eq 'F' ) or ( $ch1[$i] eq 'W' ) or ( $ch1[$i] eq 'V' ))
        and (( $ch2[$i] eq 'G' ) or ( $ch2[$i] eq 'H' ) or ( $ch2[$i] eq 'P' ) or ( $ch2[$i] eq 'S' ) or ( $ch2[$i] eq 'T' )  or ( $ch2[$i] eq 'Y' )))
        {
          #$ddc[$i] += 0.2;
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
      
       } elsif ( $key eq '9')     #   n to w
       {
       
        if ((( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'N' ) or ( $ch2[$i] eq 'D' ) or ( $ch2[$i] eq 'E' ) or ( $ch2[$i] eq 'K' ) )
        and (( $ch1[$i] eq 'G' ) or ( $ch1[$i] eq 'H' ) or ( $ch1[$i] eq 'P' ) or ( $ch1[$i] eq 'S' ) or ( $ch1[$i] eq 'T' ) or ( $ch1[$i] eq 'Y' )))
       {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
    
    } elsif ( $key eq '10')  #  w to n
       {
       
         if ( (( $ch1[$i] eq 'R' ) or ( $ch1[$i] eq 'N' ) or ( $ch1[$i] eq 'D' ) or ( $ch1[$i] eq 'E' ) or ( $ch1[$i] eq 'K' ) )
        and  (( $ch2[$i] eq 'G' ) or ( $ch2[$i] eq 'H' ) or ( $ch2[$i] eq 'P' ) or ( $ch2[$i] eq 'S' ) or ( $ch2[$i] eq 'T' )  or ( $ch2[$i] eq 'Y' )))
       {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
       } else {next;}
       
     } elsif ( $key eq '11')   # h
       {
       
        if ((( $ch1[$i] eq 'A' ) or ( $ch1[$i] eq 'C' ) or ( $ch1[$i] eq 'I' ) or ( $ch1[$i] eq 'L' ) or ( $ch1[$i] eq 'M' ) or ( $ch1[$i] eq 'F' ) or ( $ch1[$i] eq 'W' ) or ( $ch1[$i] eq 'V' )))
       {

          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
           
          
    } elsif ( $key eq '12')   # w
       {
       
         if ((( $ch1[$i] eq 'R' ) or ( $ch1[$i] eq 'N' ) or ( $ch1[$i] eq 'D' ) or ( $ch1[$i] eq 'E' ) or ( $ch1[$i] eq 'K' ) ) )
       {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
        
    } elsif ( $key eq '13')  # n
       {

         if ((( $ch1[$i] eq 'G' ) or ( $ch1[$i] eq 'H' ) or ( $ch1[$i] eq 'P' ) or ( $ch1[$i] eq 'S' ) or ( $ch1[$i] eq 'T' ) or ( $ch1[$i] eq 'Y' )  ) )
       {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
          #  Polar (R, N, D, Q, E, H, K, S, T, Y)
          #  Nonpolar (A, C, G, I, L, M, F, P, W, V)     
             
    }  elsif ( $key eq '14')   # 0h
       {
       
        if ((( $ch2[$i] eq 'A' ) or ( $ch2[$i] eq 'C' ) or ( $ch2[$i] eq 'I' ) or ( $ch2[$i] eq 'L' ) or ( $ch2[$i] eq 'M' ) or ( $ch2[$i] eq 'F' ) or ( $ch2[$i] eq 'W' ) or ( $ch2[$i] eq 'V' )))
       {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
           
          
    } elsif ( $key eq '15')   # 0w
       {
       
         if ((( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'N' ) or ( $ch2[$i] eq 'D' ) or ( $ch2[$i] eq 'E' ) or ( $ch2[$i] eq 'K' ) ) )
       {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
        
    } elsif ( $key eq '16')  # 0n
       {

         if ((( $ch2[$i] eq 'G' ) or ( $ch2[$i] eq 'H' ) or ( $ch2[$i] eq 'P' ) or ( $ch2[$i] eq 'S' ) or ( $ch2[$i] eq 'T' ) or ( $ch2[$i] eq 'Y' )  ) )
       {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;}      
          #  Polar (R, N, D, Q, E, H, K, S, T, Y)
          #  Nonpolar (A, C, G, I, L, M, F, P, W, V)        
     } elsif ( $key eq '17')   # pos pos       
        {
       
        if ((( $ch1[$i] eq 'R' ) or ( $ch1[$i] eq 'H' ) or ( $ch1[$i] eq 'K' ) ) and (( $ch2[$i] eq 'R' ) or ( $ch2[$i] eq 'H' ) or ( $ch2[$i] eq 'K' ) ) )
        {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;} 
        
     } elsif ( $key eq '18')  # nega nega
       {
       
        if ((( $ch1[$i] eq 'D' ) or ( $ch1[$i] eq 'E' )  ) and (( $ch2[$i] eq 'D' ) or ( $ch2[$i] eq 'E' ) ) )
        {
          push(@exp,$dde[$i]); 
          push(@calc,$ddc[$i]);
        } else {next;} 
         
     } else {
      
      push(@exp,$dde[$i]); 
      push(@calc,$ddc[$i]);
      
      }

      #print "Chains $f[6] exp $f[7]\n";
      my $dE = $dde[$i]-$ddc[$i]; # delta Energy (energy error)
      $rmse += ($dE)**2;
      $MAE += abs($dE);
      
#     if (( $dde[$i] == 0) and (abs($ddc[$i]) > 0.01 ))
#      {
#        $NCddG++; # Count Non-Conclusive
#      }
#      else 
#      {

#        if( $dde[$i] >= $thr) # Stabilizing mutation, ddG >0 (Positive)

        if( $dde[$i] > $thrE or ($dde[$i] == 0 and $contzero % 2 == 0) ) # Stabilizing mutation, ddG >0 (Positive)
        {
          if($ddc[$i] > $thr )
          {
            $TPddG++; # Count True Positives
            $TPavg += $dde[$i]; # True Positives ddG average
            $TPerr += abs($dde[$i]-$ddc[$i]); # True Positives ddG error average |Eexp-Ecalc|
          }
          else
          {
            $FNddG++; # Count False Negatives
          }
        }
        else # Destabilizing mutation, ddG < 0 (Negative)
        {
          if($ddc[$i] <= $thr)
          {
            $TNddG++; # Count True Negatives
            $TNavg += $dde[$i]; # True Negatives ddG average
            $TNerr += abs($dde[$i]-$ddc[$i]); # True Negatives ddG error average |Eexp-Ecalc|
          }
          else
          {
            $FPddG++; # Count False Positives
          }
        }
#      }
      
      $contzero++ if( $dde[$i] == 0 ); # "Solving" the Experimental Zeros in ddG issue
      
      # Classification stuff (according to Rosetta ddG_monomer paper)
      my $diff = abs($classExp[$i] - $classCalc[$i]); # Compute diference
      if( $diff == 0 ) { $nSc++; } # Count number of "Same class"
      elsif( $diff == 1 ) { $nOb1++; } # Count number of "Off by one"
      elsif( $diff == 2 ) { $nOb2++; } # Count number of "Off by two"

      $cont++;
    }
  }

  # Compute Pearson's Cross-Correlation
  my $PCC = 0; # Pearson's Cross-Correlation
  if($cont > 1)
  {
    my $eea = average(\@exp); # Experimental Energy Average
    my $cea = average(\@calc); # Experimental Energy Average
    my $ees = sigma(\@exp,$eea); # Experimental Energy Sigma (standard deviation)
    my $ces = sigma(\@calc,$cea); # Experimental Energy Sigma (standard deviation)
    my $dot = dot(\@exp, \@calc, $eea, $cea); # Dot product sum
    if($ees * $ces != 0)
    {
      $PCC = $dot / ($cont * $ees * $ces); # Pearson's Cross-Correlation
    }
  }

  if($TPddG == 0)
  {
    $TPavg = 0.0; # True Positives ddG average
    $TPerr = 0.0; # True Positives ddG error average |Eexp-Ecalc|
  }
  else
  {
    $TPavg /= $TPddG; # True Positives ddG average
    $TPerr /= $TPddG; # True Positives ddG error average |Eexp-Ecalc|
  }

  if($TNddG == 0)
  {
    $TNavg = 0.0; # True Negatives ddG average
    $TNerr = 0.0; # True Negatives ddG error average |Eexp-Ecalc|
  }
  else
  {
    $TNavg /= $TNddG; # True Negatives ddG average
    $TNerr /= $TNddG; # True Negatives ddG error average |Eexp-Ecalc|
  }

  my $TPR; # RECall = Sensitivity (SN) = True Positive Rate (TPR)
  if (($TPddG+$FNddG)==0) { $TPR=0; }
  else {
    $TPR = $TPddG / ( $TPddG + $FNddG ); 
  }

  my $FPR; # False Positive Rate (FPR)
  if (($TNddG + $FPddG)==0) { $TPR=0; }
  else {
    $FPR = $FPddG / ( $TNddG + $FPddG ); 
  }

  my $SPE;
  if( ($TNddG+$FPddG) == 0) { $SPE=0; }
  else {
    $SPE = $TNddG / ($TNddG+$FPddG); # SPecificity = 1 - False Positive Rate (FPR)
  }

  my $pre;
  if (($TPddG+$FPddG)==0) {$pre=0;}
  else {
  $pre = $TPddG/($TPddG+$FPddG);
  }

  my $npv;
  if (($TNddG+$FNddG)==0) {$npv=0;}
  else {
  $npv = $TNddG/($TNddG+$FNddG);
  }

  if ($cont==0)
  {
    $rmse = 0;
    $MAE = 0;
  }
  else 
  {
    $rmse = sqrt($rmse/$cont);
    $MAE /= $cont;
  }

  my $MCC;
  if(($TPddG + $FPddG)*($TPddG + $FNddG)*($TNddG + $FPddG)*($TNddG + $FNddG) != 0)
  {
    $MCC = ($TPddG * $TNddG - $FPddG * $FNddG) / ((($TPddG + $FPddG)*($TPddG + $FNddG)*($TNddG + $FPddG)*($TNddG + $FNddG))**0.5); # Matthews Correlation Coefficient (MCC)
  }
  
  my $tot = $TPddG + $FPddG + $TNddG + $FNddG + $NCddG;
  #printf " Total= %d  TP= %d  FP= %d  TN= %d  FN= %d  NC= %d\n", $tot, $TPddG, $FPddG, $TNddG, $FNddG, $NCddG;

  my $fl=" ".$key;
  
  if ($key=="1") { $fl="sb"; }
  if ($key=="2") { $fl="bs"; }
  if ($key=="3") { $fl="pn"; }
  if ($key=="4") { $fl="np"; }
  if ($key=="5") { $fl="wh";}
  if ($key=="6") { $fl="hw";}
  if ($key=="7") { $fl="nh"; }
  if ($key=="8") { $fl="hn"; }
  if ($key=="9") { $fl="nw"; }  
  if ($key=="10") { $fl="wn"; }  
  if ($key=="11") { $fl="h0"; }
  if ($key=="12") { $fl="n0"; }
  if ($key=="13") { $fl="w0"; }  
  if ($key=="14") { $fl="0h"; }
  if ($key=="15") { $fl="0w"; }
  if ($key=="16") { $fl="0n"; }  
  if ($key=="17") { $fl="pp"; }
  if ($key=="18") { $fl="nn"; } 
  if ($tot != 0)
  {
    my $ERR = ($FPddG + $FNddG) / ($TPddG + $TNddG + $FNddG + $FPddG); # ERror Rate
    my $ACC = ($TPddG + $TNddG) / ($TPddG + $TNddG + $FNddG + $FPddG); # ACCuracy

    printf "%2s %5d %5d %5d %5d %4d %4.1f %4.1f %4d %4d %4.1f %4.1f %4d %4d %5d %5d %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %6.3f %5.1f %5.1f %5.1f %5.2f\n",$fl, $cont, ($TPddG + $FNddG), ($FPddG + $TNddG), $tot, $TPddG, $TPavg, $TPerr, $FPddG, $TNddG, $TNavg, $TNerr, $FNddG, $NCddG, $TPddG+$FPddG, $TNddG+$FNddG, $TPR, $FPR, $SPE, $pre,$npv, $ACC, $ERR,($TPR + $SPE)/2.0, $rmse, $MAE, $PCC, 100*$nSc/$tot, 100*$nOb1/$tot, 100*$nOb2/$tot, $MCC;
  }
  else
  {
    printf "%2s %5d %5d %5d %5d %4d %4.1f %4.1f %4d %4d %4.1f %4.1f %4d %4d %5d %5d %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %6.3f %5.1f %5.1f %5.1f %5.2f\n",$fl, $cont, ($TPddG + $FNddG), ($FPddG + $TNddG), $tot, $TPddG, $TPavg, $TPerr, $FPddG, $TNddG, $TNavg, $TNerr, $FNddG, $NCddG, $TPddG+$FPddG, 0, 0, 0, $pre, 0, 0, 0, 0, $rmse, 0, 0, 0, 0, 0, 0, 0;
  }
}


if(1)
{

pop(@aa); # remove last element, i.e. 'p'
printf "\n ";
 #foreach my $aa( @aa )

for my $uno (0..20) 
{
  my $aa=$aa[$uno];
  if($aa eq 'X')
  {
    printf " %4s %-4s",$aa,'#';
  }
  else
  {
    printf " %4s %-3s",$aa,'#';
  }
}
printf "\n";

my $view = 'PCC';
my $disp;

#foreach my $aa1( @aa ) 

for my $uno (0..20) 

{
  my $aa1=$aa[$uno];

  printf "%s",$aa1;
 
#  foreach my $aa2( @aa ) 
 for my $dos (0..20)  
  {
      my $aa2=$aa[$dos];
#    my ($TP,$TN,$FP,$FN,$EP,$EN,$PP,$PN) = binarystat2(\@exp,\@pred,$thrE,$thrP,\@ch1,\@ch2,$aa1,$aa2);
    my ($num, $MAE, $PCC, $nSc, $nOb1, $nOb2 ) = stats(\@dde,\@ddc,\@ch1,\@ch2,$aa1,$aa2);
    
    $disp = $MAE if($view =~ 'MAE');    
    $disp = $PCC if($view =~ 'PCC');    
#    $disp = $nSc;   

    
    # printf "num= $num  MAE= $MAE  PPC= $PCC\n";
    if( $num >= 3 )
    {
      if($aa2 eq 'X')
      {
        printf " %4.1f %-4d",$disp,$num;
      }
      else
      {
        printf " %4.1f %-3d",$disp,$num;
      }
    }
    else
    {
      if( $num > 0 )
      {
        if($aa2 eq 'X')
        {
          printf " %4s %-4d",'*',$num;
        }
        else
        {
          printf " %4s %-3d",'*',$num;
        }
      }
      else
      {
        if($aa2 eq 'X')
        {
          printf " %4s %-4d",'0',$num;
        }
        else
        {
          printf " %4s %-3d",'0',$num;
        }
      }
    }
  }
  printf "\n";
}

}




# Compute average
sub average
{
  my $data = shift;
  my $avg = 0.0;
  foreach (@{$data})
  {
    $avg += $_;
  }
  return $avg / ($#{$data} + 1);
}

# Compute standard deviation given Average (avg)
sub sigma
{
  my $data = shift;
  my $avg = shift;
  my $sig = 0.0;
  foreach (@{$data})
  {
#    printf "data= $_\n";
    $sig += ($_-$avg)**2;
  }
#  print "avg= $avg  sig= $sig\n";
  return ( $sig / ($#{$data} + 1) )**0.5;
}

# Compute the raw dot product given two arrays and their averages
sub dot
{
  my $a = shift; # data "a"
  my $b = shift; # data "b"
  my $aa = shift; # "a" average
  my $ba = shift; # "b" average
  my $dot = 0.0;
  for( my $i = 0; $i <= $#{$a}; $i++ )
  {
    $dot += ($a->[$i] - $aa) * ($b->[$i] - $ba);
  }
  return $dot;
}


# Compute basic binary statistics (TP, TN, FP and FN) for given thresholds "thrE" and "thrP",
# required to compute the confusion matrix of a binary classifier. 
sub binarystat2
{
  # Input data
  my $exp = shift;
  my $pred = shift;
  my $thrE = shift; # positive/negative threshold for experimental values
  my $thrP = shift; # positive/negative threshold for predicted values
  my $id1 = shift; # 1st aminoacid id
  my $id2 = shift; # 2nd aminoacid id
  my $aa1 = shift; # 1st aminoacid to match
  my $aa2 = shift; # 2nd aminoacid to match
  
  # Counters
  my $EP = 0; # Experimental Positives
  my $EN = 0; # Experimental Negatives
  my $PP = 0; # Predicted Positives
  my $PN = 0; # Predicted Negatives
  my $TP = 0; # True Positives
  my $TN = 0; # True Negatives
  my $FP = 0; # False Positives
  my $FN = 0; # False Negatives

  # Simple number-of-elements check  
  my $nexp = $#{$exp}+1; # Number of experimental values
  my $npred = $#{$pred}+1; # Number of predicted values
  if($nexp != $npred)
  {
    print "Error, the number of experimental ($nexp) and predicted ($npred) values is different! Forcing exit!\n";
    exit;
  }
  # print "Number of experimental ($nexp) and predicted ($npred) values read\n";

  for(my $i = 0; $i < $nexp; $i++)
  {
    if( ($id1->[$i] eq $aa1 and $id2->[$i] eq $aa2) or ($id1->[$i] eq 'X' or $id2->[$i] eq $aa2) or ($id1->[$i] eq $aa1 or $id2->[$i] eq 'X') )
    {
      if($pred->[$i] > $thrP) # Prediction is Positive?
      { # Positive
        $PP++; # count Predicted Positive
        if($exp->[$i] > $thrE) # True? (comparison with experiment)
        {
          $TP++; # count True Positive
          $EP++; # count Experimental Positive
        }
        else
        {
          $FP++; # count False Positive
          $EN++; # count Experimental Negative
        }
      }
      else 
      { # Negative
        $PN++; # count Predicted Negative
        if($exp->[$i] <= $thrE) # True? (comparison with experiment)
        {
          $TN++; # count True Negative
          $EN++; # count Experimental Negative
        }
        else
        {
          $FN++; # count False Negative
          $EP++; # count Experimental Positives
        }
      }
    }
  }
 # print "TP= $TP  TN= $TN  FP= $FP  FN= $FN  EP= $EP  EN= $EN  PP= $PP  PN= $PN\n";

 return ($TP,$TN,$FP,$FN,$EP,$EN,$PP,$PN);
}


# Compute basic statistics (e.g. MAE, PCC)
sub stats
{
  # Input data
  my $exp = shift;
  my $pred = shift;
  my $id1 = shift; # 1st aminoacid id (read from mutations file)
  my $id2 = shift; # 2nd aminoacid id (read from mutations file)
  my $aa1 = shift; # 1st aminoacid to match
  my $aa2 = shift; # 2nd aminoacid to match
  
  # Simple number-of-elements check  
  my $nexp = $#{$exp}+1; # Number of experimental values
  my $npred = $#{$pred}+1; # Number of predicted values
  if($nexp != $npred)
  {
    print "Error, the number of experimental ($nexp) and predicted ($npred) values is different! Forcing exit!\n";
    exit;
  }
  # print "Number of experimental ($nexp) and predicted ($npred) values read\n";

  my @exp2 = ();
  my @pred2 = ();
  my $MAE = 0.0;
  my $num = 0; 
  my $nSc = 0; # Number of "Same class"
  my $nOb1 = 0; # Number of "Off by one" class
  my $nOb2 = 0; # Number of "Off by two" class
  
  for(my $i = 0; $i < $nexp; $i++)
  {
#    print "($id1->[$i] eq $aa1 and $id2->[$i] eq $aa2) or ($id1->[$i] eq X and $id2->[$i] eq $aa2) or ($id1->[$i] eq $aa1 and $id2->[$i] eq X)\n";
    
#    if( ($id1->[$i] eq $aa1 and $id2->[$i] eq $aa2) or ($aa1 eq 'X' and $id2->[$i] eq $aa2) or ($id1->[$i] eq $aa1 and $aa2 eq 'X') )
    # print "($id1->[$i] eq $aa1 and $id2->[$i] eq $aa2) )";
    
    if( ($id1->[$i] eq $aa1 and $id2->[$i] eq $aa2) or ($aa1 eq 'X' and $id2->[$i] eq $aa2) or ($id1->[$i] eq $aa1 and $aa2 eq 'X')  or ($aa1 eq 'X' and $aa2 eq 'X')  )
    {
      # print "TRUE\n";
      push(@exp2, $exp->[$i]);
      push(@pred2, $pred->[$i]);
      $MAE += abs($exp->[$i] - $pred->[$i]);
      # print "$MAE += abs($exp->[$i] - $pred->[$i])\n";
      $num++;
      my $diff = abs($classExp[$i] - $classCalc[$i]); # Compute diference
      if( $diff == 0 ) { $nSc++; } # Count number of "Same class"
      elsif( $diff == 1 ) { $nOb1++; } # Count number of "Off by one"
      elsif( $diff == 2 ) { $nOb2++; } # Count number of "Off by two"
    }
    else
    {
      # print "FALSE\n";
    }
  }

   if($num > 0)
  {
  $nSc /= $num;
   $nOb1 /= $num;
   $nOb2 /= $num;
  }    
  # print "->$nSc $nOb1 $nOb2\n";
  
  my $PCC = 0.0; # Pearson's Cross-Correlation
  if($num > 0)
  {
    $MAE /= $num;
    # print "num= $num  -->  MAE= $MAE\n";
  
    # Compute Pearson's Cross-Correlation
    if($num > 1)
    {
      my $eea = average(\@exp2); # Experimental Energy Average
      my $cea = average(\@pred2); # Experimental Energy Average
      my $ees = sigma(\@exp2,$eea); # Experimental Energy Sigma (standard deviation)
      my $ces = sigma(\@pred2,$cea); # Experimental Energy Sigma (standard deviation)
      my $dot = dot(\@exp2, \@pred2, $eea, $cea); # Dot product sum
      if($ees * $ces != 0)
      {
        $PCC = $dot / ($num * $ees * $ces); # Pearson's Cross-Correlation
      }
    }
  }

  return ($num, $MAE, $PCC, $nSc, $nOb1, $nOb2);
}

