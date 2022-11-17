
#Results from S461 dataset

```sh
S461_KORPM_DIR.txt   Direct subset
         *_prc.txt   Corresponding PRC Curves      
```

Equivalent result files for Cartddg, FoldX, Evo, PopMs, Dynamut2, DDGun3D, ACDCNN methods

#Comparative TABLES

## Direct case 
```sh
               S     D     T   TP  avg  err   FP   TN  avg  err   FN    P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE   MAE   PCC    Sc    Ob1   Ob2  MCC     AUC_R   AUC_P   THROC    Sen     Spe     BMCC    TH     Sen     Spe
KORPM         74   387   461   43  0.6  0.4   81  306 -1.6  0.9   31   124   337 0.581 0.791 0.347 0.908 0.757 0.686 1.208 0.906 0.570  67.0  32.5   0.4  0.31    0.771   0.362  -0.198   0.730   0.703   0.330  -0.198   0.730   0.703
KORPM*        74   387   461   43  0.6  0.4   79  308 -1.6  0.9   31   122   339 0.581 0.796 0.352 0.909 0.761 0.688 1.204 0.904 0.575  66.6  33.0   0.4  0.31    0.775   0.381  -0.201   0.730   0.698   0.325  -0.201   0.730   0.698
CartDD        74   386   460   26  0.6  0.9   27  359 -1.5  3.2   48    53   407 0.351 0.930 0.491 0.882 0.837 0.641 3.590 2.928 0.605  58.0  39.6   2.4  0.32    0.783   0.407  -2.396   0.784   0.689   0.363  -0.154   0.405   0.925
FoldX         74   387   461   30  0.6  0.6   82  305 -1.5  1.1   44   112   349 0.405 0.788 0.268 0.874 0.727 0.597 1.909 1.263 0.301  65.9  29.3   4.8  0.17    0.674   0.232  -0.481   0.703   0.618   0.237  -0.481   0.703   0.618
EvoFF         74   387   461   33  0.5  0.5   86  301 -1.6  0.8   41   119   342 0.446 0.778 0.277 0.880 0.725 0.612 1.275 0.971 0.463  66.4  31.5   2.2  0.19    0.679   0.258  -0.605   0.770   0.561   0.243  -0.605   0.770   0.561
Dynamut2      74   387   461   27  0.6  0.4   99  288 -1.6  0.9   47   126   335 0.365 0.744 0.214 0.860 0.683 0.555 1.274 0.961 0.500  64.6  34.1   1.3  0.09    0.649   0.222  -0.352   0.730   0.527   0.205  -0.673   0.865   0.403
PopMs         74   387   461   23  0.4  0.2   31  356 -1.5  0.7   51    54   407 0.311 0.920 0.426 0.875 0.822 0.615 1.022 0.763 0.611  72.0  27.5   0.4  0.26    0.768   0.349  -0.471   0.689   0.711   0.318  -0.373   0.635   0.762
DDgun3D       74   387   461   33  0.6  0.5   48  339 -1.5  0.7   41    81   380 0.446 0.876 0.407 0.892 0.807 0.661 1.110 0.806 0.631  73.1  25.8   1.1  0.31    0.765   0.361  -0.287   0.676   0.682   0.311   0.010   0.446   0.876
ThermoNet     74   387   461   37  0.6  0.5  106  281 -1.6  0.9   37   143   318 0.500 0.726 0.259 0.884 0.690 0.613 1.236 0.931 0.552  68.3  30.8   0.9  0.18    0.697   0.274  -0.227   0.784   0.540   0.238  -0.227   0.784   0.540
ACDCNN        74   387   461   22  0.6  0.4   43  344 -1.5  0.7   52    65   396 0.297 0.889 0.338 0.869 0.794 0.593 1.065 0.775 0.605  72.7  27.1   0.2  0.20    0.736   0.340  -0.215   0.595   0.749   0.279  -0.166   0.541   0.793
```

-KORPM  results trained with balanced dataset [Id25c03_1merNCLB_S461.txt](../Id25c03_1merNCLB_S461.txt) exluding all proteins with a sequence identity <0.25 respect S461 proteins 
-KORPM* results trained with balanced dataset [Id25c03_1merNCLB.txt](../Id25c03_1merNCLB.txt) in where just the Ssym mutations where excluded (default korpm weights)  


