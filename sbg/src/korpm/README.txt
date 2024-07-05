
###################################################
#  --- KORP v1.18 release - August 1st, 2018 ---  #
#  --- http://chaconlab.org/modeling/korp    ---  #
###################################################


----------------------------------------------
KORP (Knowledge based ORientational Potential)
----------------------------------------------

Pre-compiled 64-bit LINUX binaries (v1.18) are available in the bin/ directory. 
Please, choose the appropriate release for your system:

 [Program]        Compiler              Libraries     Linkage     Version
 ------------------------------------------------------------------------
  korpe*          Intel icpx (v16.0.1)     -           static        1.18 
  korpe_gcc       GNU gcc     (v4.6.3)     -           static        1.18   

* Intel compiled binary is the faster alternative (10-20%).



##################
# KORPE TUTORIAL #
##################

The basic instructions to compute the KORP energies for both loops and proteins
are described here. 

More info about KORP is available at http://chaconlab.org/modeling/korp


---------------
 Loop modeling
---------------

The KORP energy computation for loop ensembles is fast and easy with "korpe".


> Single run
------------

For a single case run, just provide a PDB file (filaname.pdb) of the protein 
and a Multi-PDB file with all the loops to be assessed (multi_loops.pdb).

#> korpe filename.pdb --loops multi_loops.pdb --score_file path_to_map/korp_v101_loop.bin -o basename

In the output basename.txt file you will find a table with the KORP energies.

For example in the directory rcd12/ you can run:

#> korpe 1h3n12A559.pdb --loops 1h3n12A559_closed.pdb --score_file ../korp_loop_v101.bin -o output

and check the output.txt file. 


> Multiple run 
--------------

For a multiple cases run, use a plain text file with the base names of the 
protein PDBs, i.e. a list of PDB file names without the .pdb extension 
(filename_list.txt), and the common suffix of the loops Multi-PDB files 
(_loops.pdb) instead.

#> korpe filename_list.txt --loops _loops.pdb --score_file ../korp_loop_v101.bin -o suffix 

For each XXX entry in the filename_list.txt you will find a file name 
XXX_suffix.txt with the corresponding a table with the KORP energies.

For example in the directory rcd12/ you can run:

#> korpe rcd12_ids.txt --loops _closed.pdb --score_file ../korp_loop_v101.bin -o output 


*IMPORTANT
----------

- The presence of one N-terminal and one C-terminal residue (anchors)
in the Multi-PDB file is mandatory.

- Both the protein PDB and the loops Multi-PDB must be numerated 
consistently.

- Note that the loop modeling energy map korp_loop_v101.bin is used here.

- KORP energy is side-chain independent, it only requires the coordinates
of the N, CA, and C mainchain atoms.

- If benchmarking, add the --rmsd flag and provide a complete PDB (with native
loop) to obtain complete KORP energy statistics in the output_score.txt file.
For example, use these comands for single and multiple runs, respectively:

#> korpe 1h3n12A559.pdb --loops 1h3n12A559_closed.pdb --score_file ../korp6Dv101loop.bin -o output --rmsd

#> korpe rcd12_ids.txt --loops _closed.pdb --score_file ../korp6Dv101loop.bin -o output --rmsd 


------------------
 Protein modeling
------------------

The KORP energy computation for protein model is fast and easy with "korpe", 


> Single case run
-----------------
For a single case run, just provide a PDB file of the protein and the
corresponding energy map. For example, in directory CASP12DCsel20/ run:

#> korpe T0902D1_s432m2.pdb --score_file ../korp_prot_v101.bin

The computed energy will be prompted to screen.


> Multiple run
--------------
For a multiple cases run, just input a plain text file with the base
names of the PDBs, i.e. a list of PDB file names without the .pdb
extension. In directory CASP12DCsel20/ run:

#> korpe casp12dcsel20.txt --score_file ../korp_prot_v101.bin -o korp6D
                                           

*IMPORTANT
----------

- Note that the protein modeling energy map korp_prot_v101.bin is used here.

- Use a correct residue numeration in the PDB, it may be important for
bonding/non-bonding contacts discrimination.

- KORP energy is side-chain independent, it only requires the coordinates 
of the N, CA, and C mainchain atoms.


##################
# KORPE timmings #
##################

Benchmark     Targets  Cases   Time*
                [#]     [#]     [s]
-------------------------------------------
rcd6            502   502000    197
rcd8            214   214000    104
rcd10           100   100000     59
rcd12            42    42000     29
casp10sel20      99     1834     16  (15)
casp11sel20      81     1610     14  (13)
casp12sel20      44      889      9   (8)
casp10best150    99    11729     88  (76)
casp11best150    81    11423     84  (71)
casp12best150    44     5849     45  (39)
rosetta41        41     4599     23  (22)
rosetta3DR41     41     4141     17  (16)
itasser56        56    24650    101  (99)
itasser3DR56     56    22456     73  (65)
3DRobot         200    60200    335 (293)

*Data obtained using the Intel compiled version of korpe.


----------
REFERENCES
----------

Please, cite our work if our tools or benchmarks become useful for your
research.

> KORP energy and benchmarks
----------------------------
Lopez-Blanco JR and Chacon P (2018). KORP: Knowledge-based 6D
potential for protein structure prediction. Bioinformatics (to be
published)

> Improved RCD method
---------------------
Lopez-Blanco JR, Canosa-Valls AJ, Li Y, and Chacon P (2016). RCD+: Fast loop
modeling server. NAR (DOI: 10.1093/nar/gkw395).

> Original RCD method
---------------------
Chys P and Chacon P (2013). Random coordinate descent with spinor-matrices and
geometric filters for efficient loop closure. J. Chem. Theory Comput.
9:1821-1829.


-------
CONTACT
-------

Please, feel free to contact with us!
(Any suggestion or bug report is welcome)

Jose Ramon Lopez-Blanco (PhD.)
jrlopez@iqfr.csic.es

Pablo Chacon (PhD.)
pablo@chaconlab.org








######################################
#  ---  BENCHMARKS DESCRIPTION  ---  #
# http://chaconlab.org/modeling/korp #
######################################


---------------------------
Protein modeling benchmarks
---------------------------


> CASP benchmarks
-----------------
The server predicted decoys for CASP 10, 11, and 12 protein modeling
challenges were obtained from, for example,
http://predictioncenter.org/download_area/CASP12/server_predictions/.

The structures for CASP stages 1 (sel20) and 2 (best150) were trimmed
to single domains according to the CASP's domain definition found at
their website (e.g.
http://predictioncenter.org/casp12/domains_summary.cgi). This way, and
for a fair comparison, both native structures and the corresponding
decoys have the same number of atoms.

Benchmark     Targets    Cases
id              [#]      [#]
-----------------------------
casp10dcsel20    99      1834
casp11dcsel20    81      1610
casp12dcsel20    44       889
casp10dcbest150  99     11729
casp11dcbest150  81     11423
casp12dcbest150  44      5849

The input files with all the decoy names required for korpe multi-run
are provided. They are named <benchmark_id>.txt, i.e. casp12dcsel20.txt.

The CASPxxDCyy/ directories contain the corresponding CASP decoys
trimmed to domains for the sel20 and best150 sets, where the "xx"
indicates the CASP edition, i.e. 10, 11, or 12, and "yy" the benchmark
sets, i.e. sel20 or best150.

Inside you will find all the native structures (e.g. T0676D1.pdb) and
all the decoys (e.g. T0676D1_s114m3.pdb). File name indicates,
respectively, the CASP target number (i.e. T0676), the domain (i.e.
D1), the server number (i.e. s114), and the model number (i.e. m3).
The input lists of cases for korpe are named caspXXdcYY.txt.

Note that a plain-text table ("TMall" file, e.g. casp12dcsel20_TMall.txt) with
the C-alpha RMSDs, scores (TM, MS and GDT_TS), and the energies
computed using others potentials (ICOSA, RW+, GOAP, ORDER_AVE, VoroQ,
DOSP, and CSF) can be found inside each benchmark sub-directory.

*Commands used to obtain paper's results
----------------------------------------
time korpe casp12dcsel20.txt --score_file ../korp6Dv101prot.bin -o korp6D
time korpe casp12dcbest150.txt --score_file ../korp6Dv101prot.bin -o korp6D
time korpe casp11dcsel20.txt --score_file ../korp6Dv101prot.bin -o korp6D
time korpe casp11dcbest150.txt --score_file ../korp6Dv101prot.bin -o korp6D
time korpe casp10dcsel20.txt --score_file ../korp6Dv101prot.bin -o korp6D
time korpe casp10dcbest150.txt --score_file ../korp6Dv101prot.bin -o korp6D


-------
CONTACT
-------

Please, feel free to contact with us!
(Any suggestion or bug report is welcome)

Jose Ramon Lopez-Blanco (PhD.)
jrlopez@iqfr.csic.es

Pablo Chacon (PhD.)
pablo@chaconlab.org




######################################
#  ---  BENCHMARKS DESCRIPTION  ---  #
# http://chaconlab.org/modeling/korp #
######################################


> Classic benchmarks
--------------------

The classic Rosetta and I-Tasser protein modeling benchmarks were
downloaded from http://fons.bakerlab.org/rosetta_decoys_62proteins.tgz
and https://zhanglab.ccmb.med.umich.edu/decoys, respectively.

Benchmark     Targets  Cases
id              [#]     [#]
----------------------------
rosetta41        41     4599
itasser56        56    24650

The rosetta41/ directory contains well-scoring Rosetta protein decoys
(e.g. 1scj_073.pdb) and their native crystal structures (e.g.
1scj.pdb). Natives were refined in Rosetta to produce 20 relaxed
decoys (e.g. 1scj_rx_15.pdb) relieving clashes and, in some cases,
repacking rotamers. See the README in the downloaded file for more
info.

The itasser56/ directory contains the decoys for the 56 non-homologous
small proteins of the Decoy Set II. This decoy set is a non-rudundant
subset of the I-TASSER decoys (Set-I). The raw decoys were first
generated by I-TASSER ab initio simulation, which were then
structurally refined by GROMACS4.0 using OPLS-AA force field with the
purpose of removing steric clashes and refining torsion angles. Each
set includes 300-500 decoys, which was used to test the atomic
potentials, including RW and RW+ potentials, see: J Zhang and Y Zhang,
A Distance-Dependent Atomic Potential Derived from Random-Walk Ideal
Chain Reference State for Protein Fold Selection and Structure
Prediction. PLoS One, vol 5, e15386 (2010).

Note that a plain-text table ("TMall" file, e.g. rosetta41_TMall.txt) with
the C-alpha RMSDs, scores (TM, MS and GDT_TS), and the energies
computed using others potentials (ICOSA, RW+, GOAP, ORDER_AVE, VoroQ,
DOSP, and CSF) can be found inside each benchmark sub-directory.

*Commands used to obtain paper's results
----------------------------------------
time korpe rosetta41.txt --score_file ../korp6Dv101prot.bin -o korp6D
time korpe itasser56.txt --score_file ../korp6Dv101prot.bin -o korp6D


-------
CONTACT
-------

Please, feel free to contact with us!
(Any suggestion or bug report is welcome)

Jose Ramon Lopez-Blanco (PhD.)
jrlopez@iqfr.csic.es

Pablo Chacon (PhD.)
pablo@chaconlab.org




######################################
#  ---  BENCHMARKS DESCRIPTION  ---  #
# http://chaconlab.org/modeling/korp #
######################################


> 3DRobot benchmarks
--------------------

The original 3DRobot and the 3DRobot versions of the classic
benchmarks were downloaded from:
https://zhanglab.ccmb.med.umich.edu/3DRobot/decoys/

Benchmark     Targets  Cases
id              [#]     [#]
----------------------------
3DRobot         200    60200
rosetta3DR41     41     4141
itasser3DR56     56    22456

Note that a plain-text table ("TMall" file, e.g. 3DR200_TMall.txt) with
the C-alpha RMSDs, scores (TM, MS and GDT_TS), and the energies
computed using others potentials (ICOSA, RW+, GOAP, ORDER_AVE, VoroQ,
DOSP, and CSF) can be found inside each benchmark sub-directory.

*Commands used to obtain paper's results
----------------------------------------
time korpe 3DR200.txt --score_file ../korp6Dv101prot.bin -o korp6D
time korpe rosetta3DR41.txt --score_file ../korp6Dv101prot.bin -o korp6D
time korpe itasser3DR56.txt --score_file ../korp6Dv101prot.bin -o korp6D


-------
CONTACT
-------

Please, feel free to contact with us!
(Any suggestion or bug report is welcome)

Jose Ramon Lopez-Blanco (PhD.)
jrlopez@iqfr.csic.es

Pablo Chacon (PhD.)
pablo@chaconlab.org




######################################
#  ---  BENCHMARKS DESCRIPTION  ---  #
# http://chaconlab.org/modeling/korp #
######################################


------------------------
Loop modeling benchmarks
------------------------

The 1000 loops per target were generated by RCD for the 6, 8, 10, and
12 residues long cases from the SOAP-loop training set.

The directories named rcdXX/ contain our RCD loop modeling benchmark
(based on the SOAP-loop training set cases), where "XX" indicates loop
length.

Benchmark     Targets  Cases
dir             [#]     [#]
----------------------------
rcd6            502   502000
rcd8            214   214000
rcd10           100   100000
rcd12            42    42000

The following RCD command generates more loop modeling benchmarks
easily. Note that the generated loops are randomly generated by RCD.

#> cd rcdXX
#> rcd rcdXX.txt -r --bench -n 1000 -x ~/PROCESS/rcd/dunbrack.bin -d 0.1 -t 0.99 -o <output_directory>

Replace any "XX" in this command by 6, 8, 10, or 12 as required. Feel free
to generate more than 1000 loops changing the -n parameter.

Detailed info about RCD usage can be found here: 
http://chaconlab.org/modeling/rcd

*Commands used to obtain paper's results
----------------------------------------
time korpe rcd6_ids.txt --loops _closed.pdb --rmsd --score_file ../korp6Dv101loop.bin -o korp6D
time korpe rcd8_ids.txt --loops _closed.pdb --rmsd --score_file ../korp6Dv101loop.bin -o korp6D
time korpe rcd10_ids.txt --loops _closed.pdb --rmsd --score_file ../korp6Dv101loop.bin -o korp6D
time korpe rcd12_ids.txt --loops _closed.pdb --rmsd --score_file ../korp6Dv101loop.bin -o korp6D


----------
REFERENCES
----------

Please, cite our work if our tools or benchmarks become useful for your
research.

> KORP energy and benchmarks
----------------------------
Lopez-Blanco JR and Chacon P (2018). KORP: Knowledge-based 6D
potential for protein structure prediction. Bioinformatics (to be
published)


-------
CONTACT
-------

Please, feel free to contact with us!
(Any suggestion or bug report is welcome)

Jose Ramon Lopez-Blanco (PhD.)
jrlopez@iqfr.csic.es

Pablo Chacon (PhD.)
pablo@chaconlab.org
