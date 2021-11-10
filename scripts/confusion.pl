#!/usr/bin/perl

#
# Mon (6/11/2020)
# Confusion Matrix stuff
# Compute Receiver-Operating (ROC) and Precision-Recall (PRC) curves
#

use strict;

if ( $#ARGV != 3 ) 
{ 
  print "USAGE:\n\t$0 <data.txt> [col_exp] [col_pred] [num]\n";
  print "\nSINTAX:\n";
  print "\t<data.txt> --> Input text file with data\n";
  print "\t[col_exp]  --> Column index for Experimental truth (>0 positives and <=0 negatives).\n";
  print "\t[col_pred] --> Column index for Predicted results (>0 positives and <=0 negatives).\n";
  print "\t[num]      --> Number of points (Pos/Neg prediction thresholds) scanned.\n";
  print "\nDESCRIPTION:\n";
  print "\tCompute Receiver-Operating (ROC) and Precision-Recall (PRC) curves from Confusion Matrix. To get:\n";
  print "\t\tReceiver-Operating Curve (ROC) --> Plot REC|TPR (y-axis) vs. FPR (x-axis).\n";
  print "\t\tPrecision-Recall Curve (PRC) --> Plot PREC|PPV (y-axis) vs. REC|TPR (x-axis).\n";
  print "\tFor example, in Gnuplot try the following line to plot PRC and ROC curves, respectively: \n\t\tplot \"prc.txt\" u 2:3 w l, \"\" u 4:2 w l, x\n";
  print "\n";
  exit 
}

# Loading stuff        
my $debug = 0; # 0= false, 1= true
my $file = $ARGV[0]; # Input text file with data
my $col_exp = int($ARGV[1]); # Column index (1,2,...) for Experimental truth (>0 positives and <=0 negatives)
my $col_pred = int($ARGV[2]); # Column index (1,2,...) for Predicted results (>0 positives and <=0 negatives)
my $num = $ARGV[3]; # Number of points (Pos/Neg prediction thresholds) scanned

# Dump log
print "# COMMAND> @ARGV\n";

my $doabs = 0;
if($num =~ /abs/)
{
  $doabs = 1;
  $num =~ s/abs//ig; # Remove "abs" tag
  print "# Using the Absolute value for the Predicted data!\n";
}
$num *= 1.0; # Force float number


# Variables
my @exp; # Experimental values
my @pred; # Predicted values
my $nentry = 0; # Number of data entries
my $thrP = 0.0; # Current threshold for Predicted Positives or Negatives
my $thrE = 0.0; # Current threshold for Experimental Positives or Negatives

# Parse input file to load relevant data
open(INPUT,"$file") or die "\nFailed to open $file\n";
while(<INPUT>)
{
  next if ($_ =~ /^#/); # Skip "#" begining lines
 
  my @line = split /\s+/, $_; # Split input line
#  push(@pred,$line[$col_pred-1]); # Load predicted value
  if($doabs == 1)
  {
    push(@pred, abs($line[$col_pred-1]) ); # Load predicted data as Absolute value 
  }
  else
  {
    push(@pred, $line[$col_pred-1] ); # Load predicted Raw value
  }
  push(@exp,$line[$col_exp-1]); # Load experimental value
  $nentry++; # Count the number of data entries
}
close INPUT;
print "# Total data entries loaded: $nentry\n";

# Check input read
if($debug)
{
  for(my $i = 0; $i < $nentry; $i++)
  {
    printf "exp= %7.3f  pred= %7.3f\n", $exp[$i], $pred[$i];
  }
}

# Get some general statistics
my ($TP,$TN,$FP,$FN,$EP,$EN,$PP,$PN) = binarystat(\@exp,\@pred,$thrE,$thrP);
printf "# Total Experimental Positives $EP and Negatives $EN using $thrE threshold.\n";
printf "# Experimental Positives ratio EP/(EN+EP)= %5.3f\n", $EP/($EN+$EP);
printf "# Total Predicted Positives $PP and Negatives $PN using $thrP threshold.\n";
printf "# Predicted Positives ratio PP/(PN+PP)= %5.3f\n", $PP/($PN+$PP);

my $minP = min(\@pred); # Minimum predicted value
my $maxP = max(\@pred); # Maximum predicted value

# prevent exact comparations 
$minP-=0.000001;
$maxP+=0.000001;
my $rangeP = $maxP - $minP; # Range of predicted values
my $delta = $rangeP / $num; # Threshold Increment
print "# Predicted data goes from $minP to $maxP (range $rangeP) --> $num samples taken in $delta steps.\n";




# open(PRC, '>', "prc.txt") or die $!;
# printf PRC "#%9s %8s %8s %8s %8s %8s\n", 'thrP', 'REC|TPR', 'PREC|PPV', 'FPR', 'SP', 'MCC';
printf "#\n#%9s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n", 'thrP', 'REC|TPR', 'PREC|PPV', 'FPR', 'SP', 'MCC', 'ACC', 'ERR', 'F05', 'F1', 'F2', 'dist2';

# Stuff to compute AUCs
my @rec; # RECall array (TPR)
my @ppv; # PREcision array (PPV)
my @fpr; # False Positives Rate array

# Stuff to compute optimal threshold from ROC
# "Using the ROC point distance to the top-left corner, we establish the best disease classification ΔΔG value for each predictor when assessing general perturbation."
# See: Greiner, M., Pfeiffer, D. & Smith, R. D. Principles and practical application of the receiver-operating characteristic analysis for diagnostic tests. Prev. Vet. Med. 45, 23–41 (2000).
my $mindist2 = 1000; # minimum distance to corner in ROC
my $thrbestROC = 0; # Best thresshold corresponding to minimum distance to corner in ROC
my $TPRbestROC = 0; # Best True Positive Rate wrt Min ROC distance
my $FPRbestROC = 0; # Best False Positive Rate wrt Min ROC distance
my $bestMCC = 0; # Best (highest) MCC
my $thrbestMCC = 0; # Best thresshold corresponding to minimum MCC
my $TPRbestMCC = 0; # Best True Positive Rate for minimum MCC
my $FPRbestMCC = 0; # Best False Positive Rate for minimum MCC

# MAIN LOOP
for(my $i = 0; $i < $num; $i++)
{
  # Update predicted threshold for current point of the curve
  $thrP = $minP + $i * $delta;
  
  # Compute Confusion Matrix basic data
  my ($TP,$TN,$FP,$FN,$EP,$EN,$PP,$PN) = binarystat(\@exp,\@pred,$thrE,$thrP);
  # print "TP= $TP  TN= $TN  FP= $FP  FN= $FN  EP= $EP  EN= $EN  PP= $PP  PN= $PN --> ";

  # Compute Basic Evaluation Measures from the confusion matrix
  # See Table 1 in: "The Precision-Recall Plot is more informative than the ROC plot when evauating binary classifiers on imbalanced datasets".
  #                 Saito and Rehmsmeier (Plos One 2015)
  my $ACC = ($TP + $TN) / ($TP + $TN + $FN + $FP); # ACCuracy
  my $ERR = ($FP + $FN) / ($TP + $TN + $FN + $FP); # ERror Rate

  my $REC = 0.0;
  $REC = $TP / ($TP + $FN) if( $TP + $FN != 0); # RECall = Sensitivity (SN) = True Positive Rate (TPR)
  my $SP = 0.0;
  $SP = $TN / ($TN + $FP) if( $TN + $FP != 0); # SPecificity = 1 - False Positive Rate (FPR)
  my $FPR = 0.0;
  $FPR = $FP / ($TN + $FP) if( $TN + $FP != 0); # False Positive Rate (FPR)
  my $PPV = 0.0;
  $PPV = $TP / ($TP + $FP) if( $TP + $FP != 0); # Positive Predictive Value = PRECision
  my $MCC = 0.0;
  my $MCC = ($TP * $TN - $FP * $FN) / ((($TP + $FP)*($TP + $FN)*($TN + $FP)*($TN + $FN))**0.5) if( ($TP + $FP)*($TP + $FN)*($TN + $FP)*($TN + $FN) != 0); # Matthews correlation coefficient
  my $F05; # F-score(0.5)
  my $F1; # F-score(1)
  my $F2; # F-score(2)
  if($PPV == 0 or $REC == 0)
  {
    $F05 = 0.0; # F-score(0.5)
    $F1 = 0.0; # F-score(1)
    $F2 = 0.0; # F-score(2)
  }
  else
  {
    $F05 = 1.5 * ($PPV * $REC) / (0.25 * $PPV + $REC); # F-score(0.5)
    $F1 = 2 * ($PPV * $REC) / ($PPV + $REC); # F-score(1)
    $F2 = 5 * ($PPV * $REC) / (4 * $PPV + $REC); # F-score(2)
  }
  
  push(@rec, $REC);
  push(@ppv, $PPV);
  push(@fpr, $FPR);

  # Best threshold stuff
  my $dist2 = ( 0 - $FPR )**2 + ( 1 - $REC )**2;
  if($dist2 < $mindist2)
  {
    $mindist2 = $dist2;
    $thrbestROC = $thrP;
    $TPRbestROC = $REC; # Best True Positive Rate (Recall)
    $FPRbestROC = $FPR; # Best False Positive Rate
  }
  
  if($MCC > $bestMCC) # Best (highest) MCC
  {
    $bestMCC = $MCC;
    $thrbestMCC = $thrP; # Best thresshold corresponding to minimum MCC
    $TPRbestMCC = $REC; # Best True Positive Rate for minimum MCC
    $FPRbestMCC = $FPR; # Best False Positive Rate for minimum MCC
  }
    
  # Dump results
  # printf PRC "%10.3f %8.6f %8.6f %8.6f %8.6f %8.6f\n", $thrP, $REC, $PPV, $FPR, $SP, $MCC;
  printf "%10.3f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f\n", $thrP, $REC, $PPV, $FPR, $SP, $MCC, $ACC, $ERR, $F05, $F1, $F2, $dist2;

# MON COMMENT
#  ($TP,$TN,$FP,$FN,$EP,$EN,$PP,$PN) = binarystat(\@exp,\@pred,0,0);
#  print "TP= $TP  TN= $TN  FP= $FP  FN= $FN  EP= $EP  EN= $EN  PP= $PP  PN= $PN --> ";

}

# Compute AUCs
my $AUC_ROC = 0.0; # Area Under Curve for ROC
my $AUC_PRC = 0.0; # Area Under Curve for PRC
for(my $i = 0; $i < scalar(@rec) - 1; $i++)
{
  $AUC_ROC += ($fpr[$i] - $fpr[$i+1]) * ($rec[$i] + $rec[$i+1]) / 2; 
  $AUC_PRC += ($rec[$i] - $rec[$i+1]) * ($ppv[$i] + $ppv[$i+1]) / 2; 
}
printf "# AUC_ROC= %-7.3f  AUC_PRC= %-7.3f  Best_ROC thr= %-7.3f TPR= %-7.3f FPR= %-7.3f  BestMCC= %-7.3f thr= %-7.3f TPR= %-7.3f FPR= %-7.3f\n",$AUC_ROC,$AUC_PRC,$thrbestROC,$TPRbestROC,$FPRbestROC,$bestMCC,$thrbestMCC,$TPRbestMCC,$FPRbestMCC;
printf "#   AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe\n";
printf "# %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",$AUC_ROC,$AUC_PRC,$thrbestROC,$TPRbestROC,(1-$FPRbestROC),$bestMCC,$thrbestMCC,$TPRbestMCC,(1-$FPRbestMCC);
exit;

#################################################################

# Compute basic binary statistics (TP, TN, FP and FN) for given thresholds "thrE" and "thrP",
# required to compute the confusion matrix of a binary classifier. 
sub binarystat
{
  # Input data
  my $exp = shift;
  my $pred = shift;
  my $thrE = shift; # positive/negative threshold for experimental values
  my $thrP = shift; # positive/negative threshold for predicted values
  
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

  my $contzero = 0;
  
  for(my $i = 0; $i < $nexp; $i++)
  {
    if( $exp[$i] > $thrE or ($exp[$i] == 0 and $contzero % 2 == 0) ) # Stabilizing mutation, ddG >0 (Positive)
    { # Stabilizing
      $EP++; # count Experimental Stabilizing (Positive)
      if($pred[$i] > $thrP) # Prediction is True?
      {
        $TP++; # count True Positive
        $PP++; # count Predicted Positive
      }
      else
      {
        $FN++; # count False Negative
        $PN++; # count Predicted Destabilizing (Negative)
      }
    }
    else
    { # Destabilizing
      $EN++; # count Experimental Destabilizing (Negative)
      if($pred[$i] <= $thrP) # Prediction is True?
      {
        $TN++; # count True Negative
        $PN++; # count Predicted Destabilizing (Negative)
      }
      else
      {
        $FP++; # count False Positive
        $PP++; # count Predicted Positive
      }
    }
    
    $contzero++ if( $exp[$i] == 0 ); # "Solving" the Experimental Zeros in ddG issue
  }
  
  
 # print "TP= $TP  TN= $TN  FP= $FP  FN= $FN  EP= $EP  EN= $EN  PP= $PP  PN= $PN\n";

 return ($TP,$TN,$FP,$FN,$EP,$EN,$PP,$PN);
}

# Get maximum value from array
sub max
{
  my $data = shift; # get data
  
  my $max = ${$data}[0]; # Load first value into max
  foreach (@{$data})
  {
    $max = $_ if($_ > $max);
  }
  
  return $max;
}

# Get minimum value from array
sub min
{
  my $data = shift; # get data
  
  my $min = ${$data}[0]; # Load first value into min
  foreach (@{$data})
  {
    $min = $_ if($_ < $min);
  }
  
  return $min;
}
