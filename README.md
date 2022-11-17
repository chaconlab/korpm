# korpm

Fast method for predicting the stability change upon mutation from 3D structure. Predicting protein stability changes upon mutation using a simple orientational potential. I. Martín-Hernández, Y. Dehouck, U. Bastolla, J.R. López-Blanco and P. Chacón (submitted).

## Usage 

The usage is very simple:  

```sh
sbg/bin/korpm input.txt --dir Ssym --score_file pot/korp6Dv1.bin -o out.txt
```
it only requires: 1) two column input file specifying both PDB file and mutation, i.e. 1BNI IA76A  and 2) the path where the PDB files are located (--dir) and the KORP potential file (--score_file).  The results are also stored in the out.txt (-o option) file.

```sh
more out.txt
1BNI IA76A       -1.396
1EY0 TA44V        0.108
1IHB FA82Q       -0.727
```
The mutation columns stands for: 1st letter is the wild type amino acid, 2nd is the chain ID, digits corresponds to PDB residue position, and the last letter is the mutated amino acid. We follow the standard convention ΔΔG >= 0 (positives) are stabilizing and ΔΔG < 0 (negatives) are destabilizing.

## ΔΔG Curated Databases

We extracted from [Thermomut](http://biosig.unimelb.edu.au/thermomutdb/) and [ProThermDB](https://web.iitm.ac.in/bioinfo2/prothermdb/index.html) unique mutations trying to avoid entries that potentially interact with ligands or belong to a protein-protein interfaces, and removing entries measured at extreme temperature or pH conditions. The initial curated database data comprise 3766 mutations from 141 proteins families (seq. identity <30%) with an average of ΔΔG -1.0 Kcal/mol and a standard deviation of 1.6 Kcal/mol. In total, 73% are destabilizing (ΔΔG>0) and 27% are stabilizing (ΔΔG<0). By removing mainly alanines' destabilizing mutations, we obtain a more balanced subset that includes 2344 mutations from 129 proteins families, 58% destabilizing and 42% stabilizing with an average of ΔΔG -0.7 Kcal/mol and a standard deviation of 1.6 Kcal/mol. This subset, named [Id30c08_1merNCLB.txt](Id30c08_1merNCLB.txt), was used for extract training and validation datasets for k-fold cross-validation experiments. Note that this subset is far from being perfectly balanced, e.g., the most frequent amino acid involved in the mutation still is alanine and cysteines, tryptophans, and, prolines still are underpopulated. 

<table border="0">

 <tr>
    <td>
     <img src="images/unbalanced.jpg">  </td>
    <td> 
      <img src="images/balanced.jpg">  </td>
 </tr>
  <tr>
    <td align="center" ><b style="font-size:30px"><a href="Id25c03_1merNCL.txt">Id25c03_1merNCL.txt</a> </b></td>
    <td align="center" ><b style="font-size:30px"><a href="Id25c03_1merNCLB.txt">Id25c03_1merNCLB.txt</a> </b></td>
 </tr></table>

In the directory [Db](Db) you can find all the correspond PDB files. 

## Results with Ssym

Ssym is a data set with equal number of stabilizing and destabilizing mutations compiled by Pucci et al. (https://doi.org/10.1093/bioinformatics/bty348) for which the structure of both the wild-type and mutant protein are available.  

```sh
sbg/bin/korpm Ssym.txt --dexp --dir Ssym --score_file pot/korp6Dv1.bin -o Ssym_all.txt
```
Where [Ssym.txt](Ssym.txt) is the mutations input file and the [Ssym](Ssym) directory in where the input PDB files are store. Since this input contains the experimental ΔΔG (see Appendix for small corrections) you can cross-check the predictions by: 
```sh
scripts/Mstat.pl Ssym_all.txt 10 11 2

# Current Positives|Negatives threshold (thr) is 0 (ddG >= 0 are not-destabilizing [positives] and ddG < 0 are destabilizing [negatives]).

aa     S     D     T   TP  avg  err   FP   TN  avg  err   FN     P     N   SEN   SPE   PPV   NPV   ACC  accn  RMSE  MAE   PCC    Sc    Ob1   Ob2  MCC
 X   342   342   684  264  1.5  0.9   71  271 -1.5  0.8   78   335   349 0.772 0.792 0.788 0.777 0.782 0.782 1.331 0.969 0.695  64.6  34.6   0.7  0.56
 A    97    97   194   80  1.9  1.0   15   82 -1.9  0.9   17    95    99 0.825 0.845 0.842 0.828 0.835 0.835 1.525 1.089 0.744  67.0  32.0   1.0  0.67
 V   106   106   212   82  1.2  0.8   24   82 -1.2  0.8   24   106   106 0.774 0.774 0.774 0.774 0.774 0.774 1.187 0.871 0.691  64.6  35.4   0.0  0.55
 I    68    68   136   57  1.6  0.8   11   57 -1.6  0.7   11    68    68 0.838 0.838 0.838 0.838 0.838 0.838 1.128 0.890 0.810  66.9  31.6   1.5  0.68
 L    41    41    82   30  1.7  1.0   12   29 -1.7  1.0   11    42    40 0.732 0.707 0.714 0.725 0.720 0.720 1.507 1.123 0.638  62.2  37.8   0.0  0.44
 .....etc
 
```

 
### Check ΔΔG Anti-symmetry in Ssym

```sh
sbg/bin/korpm SsymD.txt --dexp --dir Ssym --score_file pot/korp6Dv1.bin -o Ssym_dir.txt
sbg/bin/korpm SsymR.txt --dexp --dir Ssym --score_file pot/korp6Dv1.bin -o Ssym_rev.txt
paste Ssym_dir.txt  Ssym_rev.txt  > temp
awk 'function abs(x){return (x < 0) ? -x : x;} {printf "%s %s %s %s %s %s %s %f  %f %s %s\n",$1,$19, $2, $20, $10, $11,$29, ($11+$29), abs(($11+$29)), $3, $4  }' temp > KORPM_Ssym.txt
```

you can see the results in your favourite plot, for example in gnuplot:


<table border="0">

 <tr>
    <td>
<pre>
plot  "KORPM_Ssym.txt" u 6:7
stat "KORPM_Ssym.txt" u 6:7
...
  Linear Model:       y = -0.8207 x + 0.03672
  Slope:              -0.8207 +- 0.02468
  Intercept:          0.03672 +- 0.03807
  Correlation:        r = -0.8745
...
</pre>
  </td>
    <td> 
      <img src="images/gnuplot.jpg" alt="Italian Trulli">  </td>
 </tr>
</table>

### Comparative results Ssym

Here you can find some compartive results with state of the art stability prediction programs:
<font size="8" face="Courier New" >
<table border="1">
<tr><td>METHOD</td><td>RMSE</td><td>MAE</td><td>PCC</td><td>Sc</td><td>Ob1</td><td>Ob2</td><td>TPR</td><td>TNR</td><td> PPV</td><td>NPV</td><td>ACC</td><td>MCC</td><td>AROC</td><td>APRC</td></tr>
<tr><td>KORPM</td><td>1.33</td><td>0.97</td><td>0.69</td><td>64.6</td><td>34.6</td><td>0.7</td><td>0.77</td><td>0.79</td><td>0.79</td><td>0.78</td><td>0.78</td><td>0.56</td><td>0.86</td><td>0.86</td></tr>
<tr><td>Cartddg</td><td>3.44</td><td>2.63</td><td>0.63</td><td>52.3</td><td>41.1</td><td>6.6</td><td>0.58</td><td>0.87</td><td>0.82</td><td>0.67</td><td>0.73</td><td>0.47</td><td>0.81</td><td>0.82</td></tr>
<tr><td>ACDCNN</td><td>1.38</td><td>1.01</td><td>0.69</td><td>61.5</td><td>38.1</td><td>0.0</td><td>0.70</td><td>0.70</td><td>0.70</td><td>0.70</td><td>0.70</td><td>0.40</td><td>0.80</td><td>0.80</td></tr>
<tr><td>FoldX</td><td>1.86</td><td>1.29</td><td>0.54</td><td>60.1</td><td>34.5</td><td>5.4</td><td>0.55</td><td>0.78</td><td>0.71</td><td>0.63</td><td>0.66</td><td>0.33</td><td>0.74</td><td>0.75</td></tr>
<tr><td>EvoFF</td><td>1.56</td><td>1.12</td><td>0.54</td><td>61.7</td><td>34.9</td><td>3.4</td><td>0.61</td><td>0.70</td><td>0.67</td><td>0.64</td><td>0.66</td><td>0.31</td><td>0.74</td><td>0.75</td></tr>
<tr><td>PopMusic-S</td><td>1.58</td><td>1.15</td><td>0.52</td><td>56.6</td><td>42.4</td><td>1.0</td><td>0.67</td><td>0.71</td><td>0.70</td><td>0.68</td><td>0.69</td><td>0.38</td><td>0.76</td><td>0.74</td></tr>
<tr><td>Dynamut</td><td>1.88</td><td>1.37</td><td>0.38</td><td>54.4</td><td>38.2</td><td>7.5</td><td>0.21</td><td>0.88</td><td>0.64</td><td>0.53</td><td>0.55</td><td>0.13</td><td>0.62</td><td>0.62</td></tr>
<tr><td>DDGun3D</td><td>1.43</td><td>1.04</td><td>0.63</td><td>61.8</td><td>37.4</td><td>0.7</td><td>0.68</td><td>0.69</td><td>0.69</td><td>0.69</td><td>0.69</td><td>0.37</td><td>0.75</td><td>0.76</td></tr>
<tr><td>ThermoNet</td><td>1.53</td><td>1.09</td><td>0.55</td><td>58.2</td><td>40.9</td><td>0.9</td><td>0.65</td><td>0.70</td><td>0.69</td><td>0.67</td><td>0.68</td><td>0.35</td><td>0.75</td><td>0.74</td></tr>
</table>
</font>

you can find complete results in [Ssym_results](Ssym_results)  

## Results with S461

you can find complete results in [S461_results](S461_results)

### Appendix.

Corrections of original Ssym dataset based on ThermoMutDB data.

<font size="8" face="Courier New" >
<table border="1">
<tr><td>PDB</td><td>Mutation</td><td>Original </td><td>Corrected </td><td> Medline References from ThermoMutDB</td></tr>
<tr><td>1BNI</td><td>IA96V</td><td>-3.1</td><td>-0.9</td><td>2669964 (-0.90); 1569557 (0.95);  9551101 (-0.80)</td></tr>
<tr><td>1BNI</td><td>SA91A</td><td>-2.4</td><td>-1.8</td><td>14516751</td></tr>
<tr><td>1L63</td><td>SA44T</td><td>0.0</td><td>0.01</td><td>8289284</td> </tr>
<tr><td>1L63</td><td>SA38N</td><td>0.0</td><td>-0.01</td><td>1911773</td></tr>
<tr><td>1L63</td><td>LA91A</td><td>-3.9</td><td>-2.6</td><td>10545167</td></tr>
<tr><td>1L63</td><td>AA130S</td><td>1.0</td><td>-1.0</td><td>8218201</td></tr>
<tr><td>1LZ1</td><td>VA2G</td><td>-1.3</td><td>-2.29</td><td>11087397</td></tr>
<tr><td>1LZ1</td><td>VA2L</td><td>0.3</td><td>-0.05</td><td>11087397</td></tr>
<tr><td>1LZ1</td><td>IA56T</td><td>-4.3</td><td>-3.6</td><td>9010773; 10556244</td></tr>
<tr><td>1LZ1</td><td>VA74I</td><td>-1.9</td><td>0.45</td><td>11087397</td></tr>
<tr><td>1LZ1</td><td>VA74L</td><td>-0.4</td><td>0.19</td><td>11087397</td></tr>
<tr><td>1LZ1</td><td>VA74M</td><td>-0.4</td><td>0.65</td><td>11927576; 11087397</td></tr>
<tr><td>1LZ1</td><td>VA110G</td><td>-2.2</td><td>0.48</td><td>11927576; 11087397</td></tr>
<tr><td>1LZ1</td><td>VA110I</td><td>-0.8</td><td>0.86</td><td>11927576; 11087397</td></tr>
<tr><td>1LZ1</td><td>VA110F</td><td>-1.9</td><td>-0.05</td><td>11927576; 11087397</td></tr>
<tr><td>2LZM</td><td>IA3C</td><td>0.0</td><td>1.2</td><td>3405287 </td></tr>
<tr><td>2LZM</td><td>RA119E</td><td>0.0</td><td>-0.04</td><td>1942034 </td></tr> 
<tr><td>4LYZ</td><td>GA49A</td><td>-0.7</td><td>-1.9</td><td>11112507; 8771183</td></tr>
<tr><td>4LYZ</td><td>GA71A</td><td>-2.1</td><td>-0.38</td><td>11112507</td></tr>
<tr><td>4LYZ</td><td>GA102A</td><td>-1.2</td><td>0.02</td><td>11112507</td></tr>
<tr><td>4LYZ</td><td>GA117A</td><td>-0.8</td><td>-1.46</td><td>11112507</td></tr>
<tr><td>1RN1</td><td>QC25K</td><td>1.4</td><td>0.93</td><td>2663837</td></tr>
<tr><td>2LZM</td><td>RA96K</td><td>0.0</td><td>-0.001</td><td>Avoiding zero for the binary classification</td></tr>
<tr><td>1L63</td><td>SA44E</td><td>0.0</td><td>0.001</td><td>Avoiding zero for the binary classification</td></tr>
<tr><td>2LZM</td><td>KA60P</td><td>0.0</td><td>-0.001</td><td>Avoiding zero for the binary classification</td></tr>
</table>
</font>
