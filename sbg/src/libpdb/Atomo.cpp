
/*Implementation of the different methods of the Atom Class*/

#include <stdio.h>


#include "Atomo.h"

using namespace std;
/*Table with the principal elements*/
Element Table_Elements::table_Elements[NUM_ELEMENTS]=
{
//                        Symbol   Symbol group period number      weight      atomic   cov    vdw     en
  Element((char*)"Carbon",    C , (char *)"C ",      14,    2,   6.0,       12.011,     0.77,  0.77,  1.85,  2.55), // 0
  Element((char*)"Hydrogen",  H,  (char *)"H ",       1,    1,   1.0,       1.00797,    0.78,  0.3,   1.2,   2.2),  // 1
  Element((char*)"Nitrogen",  N,  (char *)"N ",      15,    2,   7.0,       14.00674,   0.71,  0.7,   1.54,  3.04), // 2
  Element((char*)"Oxygen",    O,  (char *)"O ",      16,    2,   8.0,       15.9994,    0.6,   0.66,  1.4,   3.44), // 3
  Element((char*)"Phosphorus",P,  (char *)"P ",      15,    3,   15.0,      30.973762,  1.15,  1.10,  1.9,   2.19), // 4
  Element((char*)"Sulphur",   S,  (char *)"S ",      16,    3,   16.0,      32.066,     1.04,  1.04,  1.85,  2.58), // 5
  Element((char*)"Calcium",   CA, (char *)"CA",       2,    4,   20.0,      40.078,     1.97,  1.74, 1.367,   1.0), // 6 Mon: vdw radius taken from Rosetta
  Element((char*)"Iron",      FE, (char *)"FE",       8,    4,   26.0,      55.845,     1.24,  1.16, 0.650,   1.83),// 7 Mon: vdw radius taken from Rosetta
  Element((char*)"Magnesium", MG, (char *)"MG",       2,    3,   12.0,      24.30506,   1.6,   1.36, 1.185,  1.31), // 8 Mon: vdw radius taken from Rosetta
  Element((char*)"Manganese", MN, (char *)"MN",       7,    4,   25.0,      54.93805,   1.24,  1.77,  1.0,   1.55), // 9 Mon: vdw radius manually set to 100pm
  Element((char*)"Sodium",    NA, (char *)"NA",       1,    3,   11.0,      22.989768,  1.54,  0.0,  1.364,  0.93), // 10 Mon: 2.31A is a huge vdw radius for an ion...
  Element((char*)"Zinc",      ZN, (char *)"ZN",      12,    4,   30.0,      65.39,      1.33,  1.25, 1.090,   1.65),// 11 Mon: vdw radius taken from Rosetta
  Element((char*)"Nickel",    NI, (char *)"NI",      10,    4,   28.0,      58.6934,    1.25,  1.15,  0.0,   1.91), // 12
  Element((char*)"Copper",    CU, (char *)"CU",      11,    4,   29.0,      63.546,     1.28,  1.17,  0.70,   1.9), // 13 Mon: vdw radius manually set to 70pm
  Element((char*)"Potassium", K,  (char *)"K ",       1,    4,   19.0,      39.0983,    2.27,  2.03, 1.764,  0.82), // 14 Mon: 2.31A is a huge vdw radius for an ion...
  Element((char*)"Cobalt",    CO, (char *)"CO",       9,    4,   27.0,      58.9332,    1.25,  1.16,  0.8 ,  1.88), // 15 Mon: vdw radius manually set to 80pm
  Element((char*)"Aluminum",  AL, (char *)"AL",      13,    3,   13.0,      26.981539,  1.43,  1.25,  2.05,  1.61), // 16
  Element((char*)"Bromine",   BR, (char *)"BR",      17,    4,   35.0,      79.904,     0.0,   1.14,  1.95,  2.96), // 17
  Element((char*)"Chlorine",  CL, (char *)"CL",      17,    3,   17.0,      35.4527,    0.0,   0.99,  1.81,  3.16), // 18
  Element((char*)"Chromium",  CR, (char *)"CR",       6,    4,   24.0,      51.9961,    1.25,  0.0,   0.0,   1.66), // 19
  Element((char*)"Silicon",   SI, (char *)"SI",      14,    3,   14.0,      28.0855,    1.17,  1.17,  2.0,   1.9),  // 20
  Element((char*)"Cadmium",   CD, (char *)"CD",      12,    5,   48.0,     112.411,     1.49,  1.41,  0.0,   1.69), // 21
  Element((char*)"Gold",      AU, (char *)"AU",      11,    6,   79.0,     196.96654,   1.44,  1.34,  0.0,   2.0),  // 22
  Element((char*)"Silver",    AG, (char *)"AG",      11,    5,   47.0,     107.8682,    1.44,  1.34,  0.0,   1.93), // 23
  Element((char*)"Platinum",  PT, (char *)"PT",      10,    6,   78.0,     195.08,      1.38,  1.29,  0.0,   2.54), // 24
  Element((char*)"Mercury",   HG, (char *)"HG",      12,    6,   80.0,     200.59,      1.60,  1.44,  0.0,   1.8),  // 25
  Element((char*)"Iodine",     I, (char *)"I ",      17,    5,   53.0,     126.904,     1.40,  1.39,  1.98,  2.96), // 26
  Element((char*)"Fluorine",   F, (char *)"F ",      17,    2,    9.0,      18.998,     0.42,  0.71,  1.47,  3.98), // 27
  Element((char*)"Deuterium",  D, (char *)"D ",       1,    1,    1.0,       2.01410,   0.78,  0.3,   1.2,   2.2)   // 28
};
// PyRosetta's metallic ions with available .parms file: (must have vdw radius, otherwise they do not generate a density map, e.g. in rcd)
// PyRosetta.namespace.ubuntu.release-72/database/chemical/residue_type_sets/fa_standard/residue_types/metal_ions
//residue_types/metal_ions/CA.params
//residue_types/metal_ions/CO.params
//residue_types/metal_ions/CU.params
//residue_types/metal_ions/FE.params
//residue_types/metal_ions/FE2.params
//residue_types/metal_ions/K.params
//residue_types/metal_ions/MG.params
//residue_types/metal_ions/MN.params
//residue_types/metal_ions/NA.params
//residue_types/metal_ions/ZN.params


Atom_type *atom_types;
int num_atom_type;

Atom_type atom_types_Rosseta[54]=
{
// at   name   hybridation polar    chrge   vdw     soft    deep   acept  donor    Avol    Asolpar
  {0 , "CNH2", SP2_HYBRID, APOLAR,  0.550, 2.0000, 2.0000, 0.1200, false, false, 0.01918, -0.0010},  //  1   CNH2                    // carbonyl C in Asn and Gln and guanidyl C in Arg
  {0 , "COO ", SP2_HYBRID, APOLAR,  0.620, 2.0000, 2.0000, 0.1200, false, false, 0.01918, -0.0010},  //  2   COO                     // carboxyl C in Asp and Glu
  {0 , "CH1 ", SP3_HYBRID, APOLAR, -0.090, 2.0000, 2.1400, 0.0486, false, false, 0.01918, -0.0010},  //  3   CH1                     // aliphatic C with one H (Val, Ile, Thr)
  {0 , "CH2 ", SP3_HYBRID, APOLAR, -0.180, 2.0000, 2.1400, 0.1142, false, false, 0.01918, -0.0010},  //  4   CH2                     // aliphatic C with two H (other residues)
  {0 , "CH3 ", SP3_HYBRID, APOLAR, -0.270, 2.0000, 2.1400, 0.1811, false, false, 0.01918, -0.0010},  //  5   CH3                     // aliphatic C with three H (Ala)
  {0 , "aroC", SP2_HYBRID, APOLAR, -0.115, 2.0000, 2.1400, 0.1200, false, false, 0.1108, -0.0005},  //  6   aroC                    // aromatic ring C (His, Phe, Tyr, Trp)
  {2 , "Ntrp", SP2_HYBRID, POLAR,  -0.610, 1.7500, 1.7500, 0.2384, false,  true, -0.03910, -0.0016},  //  7   Ntrp                    // N in Trp side-chain
  {2 , "Nhis", RING_HYBRID,POLAR,  -0.530, 1.7500, 1.7500, 0.2384, true,  false, -0.03910, -0.0016},  //  8   Nhis                    // N in His side-chain
  {2 , "NH2O", SP2_HYBRID, POLAR,  -0.470, 1.7500, 1.7500, 0.2384, false,  true, -0.03910, -0.0016},  //  9   NH2O                    // N in Asn and Gln side-chain
  {2 , "Nlys", SP3_HYBRID, POLAR,  -0.620, 1.7500, 1.7500, 0.2384, false,  true, -0.12604, -0.0016},  // 10   NLYS                    // N in Lys side-chain, N-terminus?
  {2 , "Narg", SP2_HYBRID, POLAR,  -0.750, 1.7500, 1.7500, 0.2384, false,  true, -0.06256, -0.0016},  // 11   Narg                    // N in Arg side-chain   **** -7.0, 07/08/01 ... too many buried Arg
  {2 , "Npro", SP2_HYBRID, APOLAR, -0.370, 1.7500, 1.8725, 0.2384, false,  true, -0.03910, -0.0016},  // 12   Npro                    // N in Pro backbone
  {3 , "OH  ", SP3_HYBRID, POLAR,  -0.660, 1.5500, 1.6585, 0.1591, true,   true, -0.04255, -0.0025},  // 13   OH                      // hydroxyl O in Ser, Thr and Tyr
  {3 , "ONH2", SP2_HYBRID, POLAR,  -0.550, 1.5500, 1.5500, 0.1591, true,  false, -0.03128, -0.0025},  // 14   ONH2                    // carbonyl O in Asn and Gln  **** -5.85, 07/08/01 ... too many buried Asn,Arg
  {3 , "OOC ", SP2_HYBRID, POLAR,  -0.760, 1.5500, 1.5500, 0.2100, true,  false, -0.06877, -0.0025},  // 15   OOC                     // carboyxl O in Asp and Glu
  {5 , "S   ", SP3_HYBRID, POLAR,  -0.160, 1.9000, 2.0330, 0.1600, false, false, 0.02576,  -0.0021},  // 16   S                       // sulfur in Cys and Met
  {2 , "Nbb ", SP2_HYBRID, POLAR,  -0.470, 1.7500, 1.8725, 0.2384, false,  true, -0.03910, -0.0016},  // 17   Nbb                     // backbone N'
  {0 , "CAbb", SP3_HYBRID, APOLAR,  0.070, 2.0000, 2.1400, 0.0486, false, false, 0.01918, -0.0010},  // 18   CAbb                    // backbone CA
  {0 , "CObb", SP2_HYBRID, APOLAR,  0.510, 2.0000, 2.1400, 0.1400, false, false, 0.01918, -0.0010},  // 19   CObb                    // backbone C'
  {0 , "OCbb", SP2_HYBRID, POLAR,  -0.510, 1.5500, 1.6585, 0.1591, true,  false, -0.03128, -0.0025},  // 20   OCbb                    // backbone O'
  {4 , "Phos", SP3_HYBRID, APOLAR, -0.160, 1.9000, 2.0330, 0.3182, false, false,  0.000, -0.0011},  // 21   Phos                    // nucleic acid P (from S)
  {1 , "Hpol", H_HYBRID,   APOLAR,  0.430, 1.0000, 1.0700, 0.0500, false, false,  0.000,  0.0005},  // 22   Hpol                    // polar H
  {1 , "Hapo", H_HYBRID,   APOLAR,  0.095, 1.2000, 1.2840, 0.0500, false, false,  0.000,  0.0005},  // 23   Hapo                    // nonpolar H
  {1 , "Haro", H_HYBRID,   APOLAR,  0.115, 1.2000, 1.2840, 0.0500, false, false,  0.000,  0.0005},  // 24   Haro                    // aromatic H
  {1 , "HNbb", H_HYBRID,   APOLAR,  0.310, 1.0000, 1.0700, 0.0500, false, false,  0.000,  0.0005},  // 25   HNbb                    // backbone HN
  {3 , "HOH ", SP3_HYBRID, POLAR,   0.000, 1.4000, 1.4000, 0.0500, true,   true,  0.000,  0.0005},  // 26   H2O                     // H2O
  {99, "F   ", SP3_HYBRID, APOLAR, -0.250, 1.7100, 1.7100, 0.0750, false, false,  0.000, -0.0011},  // 27   F                       // F    wild guess
  {18, "Cl  ", SP3_HYBRID, APOLAR, -0.130, 2.0700, 2.0700, 0.2400, false, false,  0.000, -0.0011},  // 28   Cl                      // Cl   wild guess
  {17, "Br  ", SP3_HYBRID, APOLAR, -0.100, 2.2200, 2.2200, 0.3200, false, false,  0.000, -0.0011},  // 29   Br                      // Br   wild guess
  {99, "I   ", SP3_HYBRID, APOLAR, -0.090, 2.3600, 2.3600, 0.4240, false, false,  0.000, -0.0011},  // 30   I                       // I    wild guess
  {11, "Zn2p", SP3_HYBRID, POLAR,   2.000, 1.0900, 1.0900, 0.2500, false, false,  0.000, -0.0011},  // 31   Zn2p                    // Zn2p wild guess
  {7 , "Fe2p", SP3_HYBRID, POLAR,   2.000, 0.7800, 0.7800, 0.0000, false, false,  0.000, -0.0011},  // 32   Fe2p                    // Fe2p wild guess
  {7 , "Fe3p", SP3_HYBRID, POLAR,   3.000, 0.6500, 0.6500, 0.0000, false, false,  0.000, -0.0011},  // 33   Fe3p                    // Fe3p wild guess
  {8 , "Mg2p", SP3_HYBRID, POLAR,   2.000, 1.1850, 1.1850, 0.0150, false, false,  0.000, -0.0011},  // 34   Mg2p                    // Mg2p wild guess
  {6 , "Ca2p", SP3_HYBRID, POLAR,   2.000, 1.3670, 1.3670, 0.1200, false, false,  0.000, -0.0011},  // 35   Ca2p                    // Ca2p wild guess
  {10, "Na1p", SP3_HYBRID, POLAR,   1.000, 1.3638, 1.3638, 0.0469, false, false,  0.000,  0.000},  // 36   Na1p                    // Na1p wild guess
  {14, "K1p ", SP3_HYBRID, POLAR,   1.000, 1.7638, 1.7638, 0.0870, false, false,  0.000,  0.000},  // 37   K1p                     // K1p  wild guess
  {99, "VOOC", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 38    1 ASP/GLU              // V01
  {99, "VCOO", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 39    2 ASP/GLU              // V02
  {99, "VOCN", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 40    3 ASN/GLN or BB        // V03
  {99, "VNOC", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 41    4 ASN/GLN or BB        // V04
  {99, "VCON", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 42    5 ASN/GLN or BB        // V05
  {99, "VSOG", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 43    6 SER OG               // V06
  {99, "VSCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 44    7 SER CB               // V07
  {99, "VCSG", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 45    8 CYS SG               // V08
  {99, "VCCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 46    9 CYS CB               // V09
  {99, "VRNH", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 47   10 ARG NH               // V10
  {99, "VRNE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 48   11 ARG NE               // V11
  {99, "VKNZ", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 49   12 LYS NZ               // V12
  {99, "VKCE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 50   13 LYS CE               // V13
  {99, "VHND", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 51   14 HIS ND               // V14
  {99, "VHNE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 52   15 HIS NE               // V15
  {99, "VHCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 53   16 HIS CB               // V16
  {99, "VHPO", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000}   // 54   17 HPOL                 // V17
};




//NUEVOS RADIOS DE VDW
Atom_type atom_types_ICM[54]=
{
// at   name   hybridation polar    chrge   vdw     soft    deep   acept  donor    Avol    Asolpar
  {0 , "CNH2", SP2_HYBRID, APOLAR,  0.550, 1.8100, 2.0000, 0.1200, false, false, 0.01918, -0.0010},  //  1   CNH2		      // carbonyl C in Asn and Gln and guanidyl C in Arg
  {0 , "COO ", SP2_HYBRID, APOLAR,  0.620, 1.7600, 2.0000, 0.1200, false, false, 0.01918, -0.0010},  //  2   COO		      // carboxyl C in Asp and Glu
  {0 , "CH1 ", SP3_HYBRID, APOLAR, -0.090, 2.0100, 2.1400, 0.0486, false, false, 0.01918, -0.0010},  //  3   CH1		      // aliphatic C with one H (Val, Ile, Thr)
  {0 , "CH2 ", SP3_HYBRID, APOLAR, -0.180, 1.9200, 2.1400, 0.1142, false, false, 0.01918, -0.0010},  //  4   CH2		      // aliphatic C with two H (other residues)
  {0 , "CH3 ", SP3_HYBRID, APOLAR, -0.270, 1.9200, 2.1400, 0.1811, false, false, 0.01918, -0.0010},  //  5   CH3		      // aliphatic C with three H (Ala)
  {0 , "aroC", SP2_HYBRID, APOLAR, -0.115, 1.7400, 2.1400, 0.1200, false, false, 0.1108, -0.0005},  //  6   aroC		      // aromatic ring C (His, Phe, Tyr, Trp)
  {2 , "Ntrp", SP2_HYBRID, POLAR,  -0.610, 1.6600, 1.7500, 0.2384, false,  true, -0.03910, -0.0016},  //  7   Ntrp		      // N in Trp side-chain
  {2 , "Nhis", RING_HYBRID,POLAR,  -0.530, 1.6500, 1.7500, 0.2384, true,  false, -0.03910, -0.0016},  //  8   Nhis		      // N in His side-chain
  {2 , "NH2O", SP2_HYBRID, POLAR,  -0.470, 1.6200, 1.7500, 0.2384, false,  true, -0.03910, -0.0016},  //  9   NH2O		      // N in Asn and Gln side-chain
  {2 , "Nlys", SP3_HYBRID, POLAR,  -0.620, 1.6700, 1.7500, 0.2384, false,  true, -0.12604, -0.0016},  // 10   NLYS		      // N in Lys side-chain, N-terminus?
  {2 , "Narg", SP2_HYBRID, POLAR,  -0.750, 1.6700, 1.7500, 0.2384, false,  true, -0.06256, -0.0016},  // 11   Narg		      // N in Arg side-chain   **** -7.0, 07/08/01 ... too many buried Arg
  {2 , "Npro", SP2_HYBRID, APOLAR, -0.370, 1.6700, 1.8725, 0.2384, false,  true, -0.03910, -0.0016},  // 12   Npro		      // N in Pro backbone
  {3 , "OH  ", SP3_HYBRID, POLAR,  -0.660, 1.5400, 1.6585, 0.1591, true,   true, -0.04255, -0.0025},  // 13   OH		      // hydroxyl O in Ser, Thr and Tyr
  {3 , "ONH2", SP2_HYBRID, POLAR,  -0.550, 1.5200, 1.5500, 0.1591, true,  false, -0.03128, -0.0025},  // 14   ONH2		      // carbonyl O in Asn and Gln  **** -5.85, 07/08/01 ... too many buried Asn,Arg
  {3 , "OOC ", SP2_HYBRID, POLAR,  -0.760, 1.4900, 1.5500, 0.2100, true,  false, -0.06877, -0.0025},  // 15   OOC		      // carboyxl O in Asp and Glu
  {4 , "S   ", SP3_HYBRID, POLAR,  -0.160, 1.9400, 2.0330, 0.1600, false, false, 0.02576,  -0.0021},  // 16   S 		      // sulfur in Cys and Met
  {2 , "Nbb ", SP2_HYBRID, POLAR,  -0.470, 1.7000, 1.8725, 0.2384, false,  true, -0.03910, -0.0016},  // 17   Nbb		      // backbone N'
  {0 , "CAbb", SP3_HYBRID, APOLAR,  0.070, 1.9000, 2.1400, 0.0486, false, false, 0.01918, -0.0010},  // 18   CAbb		      // backbone CA
  {0 , "CObb", SP2_HYBRID, APOLAR,  0.510, 1.7500, 2.1400, 0.1400, false, false, 0.01918, -0.0010},  // 19   CObb		      // backbone C'
  {0 , "OCbb", SP2_HYBRID, POLAR,  -0.510, 1.4900, 1.6585, 0.1591, true,  false, -0.03128, -0.0025},  // 20   OCbb		      // backbone O'
  {4 , "Phos", SP3_HYBRID, APOLAR, -0.160, 1.9000, 2.0330, 0.3182, false, false,  0.000, -0.0011},  // 21   Phos		      // nucleic acid P (from S)
  {1 , "Hpol", H_HYBRID,   APOLAR,  0.430, 1.0000, 1.0700, 0.0500, false, false,  0.000,  0.0005},  // 22   Hpol		      // polar H
  {1 , "Hapo", H_HYBRID,   APOLAR,  0.095, 1.2000, 1.2840, 0.0500, false, false,  0.000,  0.0005},  // 23   Hapo		      // nonpolar H
  {1 , "Haro", H_HYBRID,   APOLAR,  0.115, 1.2000, 1.2840, 0.0500, false, false,  0.000,  0.0005},  // 24   Haro		      // aromatic H
  {1 , "HNbb", H_HYBRID,   APOLAR,  0.310, 1.0000, 1.0700, 0.0500, false, false,  0.000,  0.0005},  // 25   HNbb		      // backbone HN
  {3 , "HOH ", SP3_HYBRID, POLAR,   0.000, 1.4000, 1.4000, 0.0500, true,   true,  0.000,  0.0005},  // 26   H2O		      // H2O
  {27, "F   ", SP3_HYBRID, APOLAR, -0.250, 1.7100, 1.7100, 0.0750, false, false,  0.000, -0.0011},  // 27   F 		      // F    wild guess
  {18, "Cl  ", SP3_HYBRID, APOLAR, -0.130, 2.0700, 2.0700, 0.2400, false, false,  0.000, -0.0011},  // 28   Cl		      // Cl   wild guess
  {17, "Br  ", SP3_HYBRID, APOLAR, -0.100, 2.2200, 2.2200, 0.3200, false, false,  0.000, -0.0011},  // 29   Br		      // Br   wild guess
  {26, "I   ", SP3_HYBRID, APOLAR, -0.090, 2.3600, 2.3600, 0.4240, false, false,  0.000, -0.0011},  // 30   I 		      // I    wild guess
  {11, "Zn2p", SP3_HYBRID, POLAR,   2.000, 1.0900, 1.0900, 0.2500, false, false,  0.000, -0.0011},  // 31   Zn2p		      // Zn2p wild guess
  {7 , "Fe2p", SP3_HYBRID, POLAR,   2.000, 0.7800, 0.7800, 0.0000, false, false,  0.000, -0.0011},  // 32   Fe2p		      // Fe2p wild guess
  {7 , "Fe3p", SP3_HYBRID, POLAR,   3.000, 0.6500, 0.6500, 0.0000, false, false,  0.000, -0.0011},  // 33   Fe3p		      // Fe3p wild guess
  {8 , "Mg2p", SP3_HYBRID, POLAR,   2.000, 1.1850, 1.1850, 0.0150, false, false,  0.000, -0.0011},  // 34   Mg2p		      // Mg2p wild guess
  {6 , "Ca2p", SP3_HYBRID, POLAR,   2.000, 1.3670, 1.3670, 0.1200, false, false,  0.000, -0.0011},  // 35   Ca2p		      // Ca2p wild guess
  {10, "Na1p", SP3_HYBRID, POLAR,   1.000, 1.3638, 1.3638, 0.0469, false, false,  0.000,  0.000},  // 36   Na1p		      // Na1p wild guess
  {14, "K1p ", SP3_HYBRID, POLAR,   1.000, 1.7638, 1.7638, 0.0870, false, false,  0.000,  0.000},  // 37   K1p		      // K1p  wild guess
  {99, "VOOC", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 38    1 ASP/GLU	      // V01
  {99, "VCOO", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 39    2 ASP/GLU	      // V02
  {99, "VOCN", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 40    3 ASN/GLN or BB        // V03
  {99, "VNOC", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 41    4 ASN/GLN or BB        // V04
  {99, "VCON", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 42    5 ASN/GLN or BB        // V05
  {99, "VSOG", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 43    6 SER OG 	      // V06
  {99, "VSCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 44    7 SER CB 	      // V07
  {99, "VCSG", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 45    8 CYS SG 	      // V08
  {99, "VCCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 46    9 CYS CB 	      // V09
  {99, "VRNH", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 47   10 ARG NH 	      // V10
  {99, "VRNE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 48   11 ARG NE 	      // V11
  {99, "VKNZ", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 49   12 LYS NZ 	      // V12
  {99, "VKCE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 50   13 LYS CE 	      // V13
  {99, "VHND", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 51   14 HIS ND 	      // V14
  {99, "VHNE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 52   15 HIS NE 	      // V15
  {99, "VHCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  // 53   16 HIS CB 	      // V16
  {99, "VHPO", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000}   // 54   17 HPOL		      // V17
};

Atom_type atom_types_EEF1[31]=
{

  // Solo he rellenado con EEF1 at, vdw, soft, deep
// at   name   hybridation polar    chrge   vdw     soft    deep   acept  donor
//											Rmin/2	Rmin/2	emin
  {1 , "H",    SP2_HYBRID, APOLAR,  0.550, 0.8000, 0.8000, 0.04980, false, false},  //  1
  {1 , "HC",   SP2_HYBRID, APOLAR,  0.620, 0.6000, 0.6000, 0.04980, false, false},  //  2
  {1 , "HA",   SP3_HYBRID, APOLAR, -0.090, 1.4680, 1.4680, 0.04500, false, false},  //  3
  {0 , "CT",   SP3_HYBRID, APOLAR, -0.180, 2.4900, 2.4900, 0.02620, false, false},  //  4
  {0 , "C",    SP3_HYBRID, APOLAR, -0.270, 2.1000, 2.1000, 0.12000, false, false},  //  5
  {0 , "CH1E", SP2_HYBRID, APOLAR, -0.115, 2.3650, 2.3650, 0.04860, false, false},  //  6
  {0 , "CH2E", SP2_HYBRID, POLAR,  -0.610, 2.2350, 2.2350, 0.11420, false,  true},  //  7
  {0 , "CH3E", RING_HYBRID,POLAR,  -0.530, 2.1650, 2.1650, 0.18110, true,  false},  //  8
  {0 , "CR1E", SP2_HYBRID, POLAR,  -0.470, 2.1000, 2.1000, 0.12000, false,  true},  //  9
  {2 , "N",    SP3_HYBRID, POLAR,  -0.620, 1.6000, 1.6000, 0.23840, false,  true},  // 10
  {2 , "NR",  SP2_HYBRID, POLAR,  -0.750, 1.6000, 1.6000, 0.23840, false,  true},  // 11
  {2 , "NP",  SP2_HYBRID, APOLAR, -0.370, 1.6000, 1.6000, 0.23840, false,  true},  // 12
  {2 , "NH1", SP3_HYBRID, POLAR,  -0.660, 1.6000, 1.6000, 0.23840, true,   true},  // 13
  {2 , "NH2", SP2_HYBRID, POLAR,  -0.550, 1.6000, 1.6000, 0.23840, true,  false},  // 14
  {2 , "NH3", SP2_HYBRID, POLAR,  -0.760, 1.6000, 1.6000, 0.23840, true,  false},  // 15
  {2 , "NC2", SP3_HYBRID, POLAR,  -0.160, 1.6000, 1.6000, 0.23840, false, false},  // 16
  {3 , "O",   SP2_HYBRID, POLAR,  -0.470, 1.6000, 1.6000, 0.15910, false,  true},  // 17
  {3 , "OC",  SP3_HYBRID, APOLAR,  0.070, 1.6000, 1.6000, 0.64690, false, false},  // 18
  {3 , "OH1", SP2_HYBRID, APOLAR,  0.510, 1.6000, 1.6000, 0.15910, false, false},  // 19
  {3 , "OH2", SP2_HYBRID, POLAR,  -0.510, 1.7398, 1.7398, 0.07580, true,  false},  // 20
  {5 , "S",   SP3_HYBRID, APOLAR, -0.160, 1.8900, 1.8900, 0.04300, false, false},  // 21
  {5 , "SH1E",H_HYBRID,   APOLAR,  0.430, 1.8900, 1.8900, 0.04300, false, false},  // 22
  {7 , "FE",  H_HYBRID,   APOLAR,  0.095, 0.6500, 0.6500, 0.00000, false, false},  // 23
  {99 , "OS",  H_HYBRID,   APOLAR,  0.115, 1.6000, 1.6000, 0.15910, false, false},  // 24
  {99 , "CR",  H_HYBRID,   APOLAR,  0.310, 2.1000, 2.1000, 0.12000, false, false},  // 25
  {99 , "CM",  SP3_HYBRID, POLAR,   0.000, 2.4900, 2.4900, 0.02620, true,   true},  // 26
  {99, "OM",   SP3_HYBRID, APOLAR, -0.250, 1.6000, 1.6000, 0.15910, false, false},  // 27
  {99, "LP",   SP3_HYBRID, APOLAR, -0.130, 0.2245, 0.2245, 0.04598, false, false},  // 28
  {99, "HT",   SP3_HYBRID, APOLAR, -0.100, 0.8000, 0.8000, 0.04980, false, false},  // 29
  {99, "OT",   SP3_HYBRID, APOLAR, -0.090, 1.6000, 1.6000, 0.15910, false, false},  // 30
  {99, "CAL",  SP3_HYBRID, POLAR,   2.000, 1.7100, 1.7100, 0.12000, false, false}  // 31
  };

Atom_type atom_types_Sybil[45]=
{
// at   name   hybridation polar    chrge   vdw     soft    deep   acept  donor    Avol    Asolpar
  {0 , "C2  ", SP2_HYBRID, APOLAR, 0.510, 2.0000, 2.1400, 0.1400, false, false, 0.01918, -0.0010},   // 1
  {0 , "C3  ", SP3_HYBRID, APOLAR,-0.270, 2.0000, 2.1400, 0.1811, false, false, 0.01918, -0.0010},   // 2
  {0 , "Car ", SP2_HYBRID, APOLAR,-0.115, 2.0000, 2.1400, 0.1200, false, false, 0.1108, -0.0005},    // 3
  {0 , "Ccat", SP2_HYBRID, APOLAR,-0.115, 2.0000, 2.1400, 0.1200, false, false, 0.1108, -0.0005},    // 4
  {2 , "N3  ", SP3_HYBRID, POLAR, -0.620, 1.7500, 1.7500, 0.2384, true,  true, -0.12604, -0.0016},   // 5
  {2 , "Nam ", SP2_HYBRID, POLAR, -0.470, 1.7500, 1.8725, 0.2384, false,  true,-0.03910, -0.0016},   // 6
  {2 , "Npl3", SP2_HYBRID, POLAR, -0.370, 1.7500, 1.8725, 0.2384, false,  true,-0.03910, -0.0016},   // 7
  {3 , "O2  ", SP2_HYBRID, POLAR, -0.510, 1.5500, 1.6585, 0.1591, true,  false,-0.03128, -0.0025},   // 8
  {3 , "O3  ", SP3_HYBRID, POLAR, -0.660, 1.5500, 1.6585, 0.1591, true,   true,-0.04255, -0.0025},   // 9
  {3 , "Oco2", SP2_HYBRID, POLAR, -0.760, 1.5500, 1.5500, 0.2100, true,  false,-0.06877, -0.0025},   //10
  {5 , "S3  ", SP3_HYBRID, POLAR, -0.160, 1.9000, 2.0330, 0.1600, false, false, 0.02576, -0.0021},   //11  //from here it is the same as Rosseta
  {4 , "Phos", SP3_HYBRID, APOLAR, -0.160, 1.9000, 2.0330, 0.3182, false, false,  0.000, -0.0011},  // 12
  {1 , "Hpol", H_HYBRID,   APOLAR,  0.430, 1.0000, 1.0700, 0.0500, false, false,  0.000,  0.0005},  // 13
  {1 , "Hapo", H_HYBRID,   APOLAR,  0.095, 1.2000, 1.2840, 0.0500, false, false,  0.000,  0.0005},  // 14
  {1 , "Haro", H_HYBRID,   APOLAR,  0.115, 1.2000, 1.2840, 0.0500, false, false,  0.000,  0.0005},  // 15
  {1 , "HNbb", H_HYBRID,   APOLAR,  0.310, 1.0000, 1.0700, 0.0500, false, false,  0.000,  0.0005},  // 16
  {3 , "HOH ", SP3_HYBRID, POLAR,   0.000, 1.4000, 1.4000, 0.0500, true,   true,  0.000,  0.0005},  // 17
  {27, "F   ", SP3_HYBRID, APOLAR, -0.250, 1.7100, 1.7100, 0.0750, false, false,  0.000, -0.0011},  // 18
  {18, "Cl  ", SP3_HYBRID, APOLAR, -0.130, 2.0700, 2.0700, 0.2400, false, false,  0.000, -0.0011},  // 19
  {17, "Br  ", SP3_HYBRID, APOLAR, -0.100, 2.2200, 2.2200, 0.3200, false, false,  0.000, -0.0011},  // 20
  {26, "I   ", SP3_HYBRID, APOLAR, -0.090, 2.3600, 2.3600, 0.4240, false, false,  0.000, -0.0011},  // 21
  {11, "Zn2p", SP3_HYBRID, POLAR,   2.000, 1.0900, 1.0900, 0.2500, false, false,  0.000, -0.0011},  // 22
  {7 , "Fe2p", SP3_HYBRID, POLAR,   2.000, 0.7800, 0.7800, 0.0000, false, false,  0.000, -0.0011},  // 23
  {7 , "Fe3p", SP3_HYBRID, POLAR,   3.000, 0.6500, 0.6500, 0.0000, false, false,  0.000, -0.0011},  // 24
  {8 , "Mg2p", SP3_HYBRID, POLAR,   2.000, 1.1850, 1.1850, 0.0150, false, false,  0.000, -0.0011},  // 25
  {6 , "Ca2p", SP3_HYBRID, POLAR,   2.000, 1.3670, 1.3670, 0.1200, false, false,  0.000, -0.0011},  // 26
  {10, "Na1p", SP3_HYBRID, POLAR,   1.000, 1.3638, 1.3638, 0.0469, false, false,  0.000,  0.000},  //  27
  {14, "K1p ", SP3_HYBRID, POLAR,   1.000, 1.7638, 1.7638, 0.0870, false, false,  0.000,  0.000},  //  28
  {99, "VOOC", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  29
  {99, "VCOO", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  30
  {99, "VOCN", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  31
  {99, "VNOC", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  32
  {99, "VCON", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  33
  {99, "VSOG", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  34
  {99, "VSCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  35
  {99, "VCSG", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  36
  {99, "VCCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  37
  {99, "VRNH", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  38
  {99, "VRNE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  39
  {99, "VKNZ", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  40
  {99, "VKCE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  41
  {99, "VHND", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  42
  {99, "VHNE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  43
  {99, "VHCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000},  //  44
  {99, "VHPO", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false,  0.000,  0.000}   //  45
};


/*
//NUEVOS RADIOS DE VDW a partir de tablas de Julio Kovacs
Atom_type atom_types[54]=
{
// at   name   hybridation polar    chrge   vdw     soft    deep   acept  donor
  {0 , "CNH2", SP2_HYBRID, APOLAR,  0.550, 1.7000, 2.0000, 0.1200, false, false},  //  1   CNH2                    // carbonyl C in Asn and Gln and guanidyl C in Arg
  {0 , "COO ", SP2_HYBRID, APOLAR,  0.620, 1.7000, 2.0000, 0.1200, false, false},  //  2   COO                     // carboxyl C in Asp and Glu
  {0 , "CH1 ", SP3_HYBRID, APOLAR, -0.090, 1.7000, 2.1400, 0.0486, false, false},  //  3   CH1                     // aliphatic C with one H (Val, Ile, Thr)
  {0 , "CH2 ", SP3_HYBRID, APOLAR, -0.180, 1.7000, 2.1400, 0.1142, false, false},  //  4   CH2                     // aliphatic C with two H (other residues)
  {0 , "CH3 ", SP3_HYBRID, APOLAR, -0.270, 1.7000, 2.1400, 0.1811, false, false},  //  5   CH3                     // aliphatic C with three H (Ala)
  {0 , "aroC", SP2_HYBRID, APOLAR, -0.115, 1.7400, 2.1400, 0.1200, false, false},  //  6   aroC                    // aromatic ring C (His, Phe, Tyr, Trp)
  {1 , "Ntrp", SP2_HYBRID, POLAR,  -0.610, 1.5000, 1.7500, 0.2384, false,  true},  //  7   Ntrp                    // N in Trp side-chain
  {1 , "Nhis", RING_HYBRID,POLAR,  -0.530, 1.5000, 1.7500, 0.2384, true,  false},  //  8   Nhis                    // N in His side-chain
  {1 , "NH2O", SP2_HYBRID, POLAR,  -0.470, 1.5000, 1.7500, 0.2384, false,  true},  //  9   NH2O                    // N in Asn and Gln side-chain
  {1 , "Nlys", SP3_HYBRID, POLAR,  -0.620, 1.5000, 1.7500, 0.2384, false,  true},  // 10   NLYS                    // N in Lys side-chain, N-terminus?
  {1 , "Narg", SP2_HYBRID, POLAR,  -0.750, 1.5000, 1.7500, 0.2384, false,  true},  // 11   Narg                    // N in Arg side-chain   **** -7.0, 07/08/01 ... too many buried Arg
  {1 , "Npro", SP2_HYBRID, APOLAR, -0.370, 1.5000, 1.8725, 0.2384, false,  true},  // 12   Npro                    // N in Pro backbone
  {2 , "OH  ", SP3_HYBRID, POLAR,  -0.660, 1.4000, 1.6585, 0.1591, true,   true},  // 13   OH                      // hydroxyl O in Ser, Thr and Tyr
  {2 , "ONH2", SP2_HYBRID, POLAR,  -0.550, 1.4000, 1.5500, 0.1591, true,  false},  // 14   ONH2                    // carbonyl O in Asn and Gln  **** -5.85, 07/08/01 ... too many buried Asn,Arg
  {2 , "OOC ", SP2_HYBRID, POLAR,  -0.760, 1.4000, 1.5500, 0.2100, true,  false},  // 15   OOC                     // carboyxl O in Asp and Glu
  {3 , "S   ", SP3_HYBRID, POLAR,  -0.160, 1.8500, 2.0330, 0.1600, false, false},  // 16   S                       // sulfur in Cys and Met
  {2 , "Nbb ", SP2_HYBRID, POLAR,  -0.470, 1.5000, 1.8725, 0.2384, false,  true},  // 17   Nbb                     // backbone N'
  {0 , "CAbb", SP3_HYBRID, APOLAR,  0.070, 1.7000, 2.1400, 0.0486, false, false},  // 18   CAbb                    // backbone CA
  {0 , "CObb", SP2_HYBRID, APOLAR,  0.510, 1.7000, 2.1400, 0.1400, false, false},  // 19   CObb                    // backbone C'
  {0 , "OCbb", SP2_HYBRID, POLAR,  -0.510, 1.5000, 1.6585, 0.1591, true,  false},  // 20   OCbb                    // backbone O'
  {4 , "Phos", SP3_HYBRID, APOLAR, -0.160, 1.9000, 2.0330, 0.3182, false, false},  // 21   Phos                    // nucleic acid P (from S)
  {1 , "Hpol", H_HYBRID,   APOLAR,  0.430, 1.0000, 1.0700, 0.0500, false, false},  // 22   Hpol                    // polar H
  {1 , "Hapo", H_HYBRID,   APOLAR,  0.095, 1.0000, 1.2840, 0.0500, false, false},  // 23   Hapo                    // nonpolar H
  {1 , "Haro", H_HYBRID,   APOLAR,  0.115, 1.0000, 1.2840, 0.0500, false, false},  // 24   Haro                    // aromatic H
  {1 , "HNbb", H_HYBRID,   APOLAR,  0.310, 1.0000, 1.0700, 0.0500, false, false},  // 25   HNbb                    // backbone HN
  {3 , "HOH ", SP3_HYBRID, POLAR,   0.000, 1.0000, 1.4000, 0.0500, true,   true},  // 26   H2O                     // H2O
  {99, "F   ", SP3_HYBRID, APOLAR, -0.250, 1.7100, 1.7100, 0.0750, false, false},  // 27   F                       // F    wild guess
  {18, "Cl  ", SP3_HYBRID, APOLAR, -0.130, 2.0700, 2.0700, 0.2400, false, false},  // 28   Cl                      // Cl   wild guess
  {17, "Br  ", SP3_HYBRID, APOLAR, -0.100, 2.2200, 2.2200, 0.3200, false, false},  // 29   Br                      // Br   wild guess
  {99, "I   ", SP3_HYBRID, APOLAR, -0.090, 2.3600, 2.3600, 0.4240, false, false},  // 30   I                       // I    wild guess
  {11, "Zn2p", SP3_HYBRID, POLAR,   2.000, 1.0900, 1.0900, 0.2500, false, false},  // 31   Zn2p                    // Zn2p wild guess
  {7 , "Fe2p", SP3_HYBRID, POLAR,   2.000, 0.7800, 0.7800, 0.0000, false, false},  // 32   Fe2p                    // Fe2p wild guess
  {7 , "Fe3p", SP3_HYBRID, POLAR,   3.000, 0.6500, 0.6500, 0.0000, false, false},  // 33   Fe3p                    // Fe3p wild guess
  {8 , "Mg2p", SP3_HYBRID, POLAR,   2.000, 1.1850, 1.1850, 0.0150, false, false},  // 34   Mg2p                    // Mg2p wild guess
  {6 , "Ca2p", SP3_HYBRID, POLAR,   2.000, 1.3670, 1.3670, 0.1200, false, false},  // 35   Ca2p                    // Ca2p wild guess
  {10, "Na1p", SP3_HYBRID, POLAR,   1.000, 1.3638, 1.3638, 0.0469, false, false},  // 36   Na1p                    // Na1p wild guess
  {14, "K1p ", SP3_HYBRID, POLAR,   1.000, 1.7638, 1.7638, 0.0870, false, false},  // 37   K1p                     // K1p  wild guess
  {99, "VOOC", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 38    1 ASP/GLU              // V01
  {99, "VCOO", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 39    2 ASP/GLU              // V02
  {99, "VOCN", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 40    3 ASN/GLN or BB        // V03
  {99, "VNOC", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 41    4 ASN/GLN or BB        // V04
  {99, "VCON", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 42    5 ASN/GLN or BB        // V05
  {99, "VSOG", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 43    6 SER OG               // V06
  {99, "VSCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 44    7 SER CB               // V07
  {99, "VCSG", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 45    8 CYS SG               // V08
  {99, "VCCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 46    9 CYS CB               // V09
  {99, "VRNH", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 47   10 ARG NH               // V10
  {99, "VRNE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 48   11 ARG NE               // V11
  {99, "VKNZ", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 49   12 LYS NZ               // V12
  {99, "VKCE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 50   13 LYS CE               // V13
  {99, "VHND", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 51   14 HIS ND               // V14
  {99, "VHNE", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 52   15 HIS NE               // V15
  {99, "VHCB", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false},  // 53   16 HIS CB               // V16
  {99, "VHPO", H_HYBRID,   APOLAR,  0.000, 0.0000, 0.0000, 0.0000, false, false}   // 54   17 HPOL                 // V17
};
*/

/*Default constructor*/
Atom::Atom()
{
  /*The atom is created with the firts element of the table*/
  element= Table_Elements::getElement(0);
  pdbSerialNumber=0;
  memcpy(pdbName,(const char *)" C  ",5);
  //strcpy(pdbName,(const char *)" C  ");
  pdbocc=0.0;
  pdbfact=0.0;
  position[0]=0.0;
  position[1]=0.0;
  position[2]=0.0;
  charge=0;
  father=NULL;
  bonds=NULL;
  num_bonds=0;
}

/*Constructor.
 e: Element Object of the atom. The Object element is copied in another object, not assignated
 p: three-dimension position of the atom
 ch: charge of the atom*/
Atom::Atom(Element *e, Tcoor p,float ch,char name[5],int serial, float occ, float fact)
{
  element=e;
  memcpy(pdbName,(const char *)name,5);
//  strcpy(pdbName,(const char *)name);
  position[0]=p[0];
  position[1]=p[1];
  position[2]=p[2];
  charge=ch;
  pdbSerialNumber=serial;
  pdbocc=occ;
  pdbfact=fact;
  father=NULL;
  bonds=NULL;
  num_bonds=0;
}

/*Constructor from template */
Atom::Atom(Atom_type atom_type, char name[5], float *p, int serial)
{
  element= Table_Elements::getElement(atom_type.at);
  memcpy(pdbName,(const char *)name,5);

//  strcpy(pdbName,(const char *)name);

  position[0]=p[0];
  position[1]=p[1];
  position[2]=p[2];
  //printf("%s %f %f %f\n",name, position[0],position[1],position[2]);

  charge=atom_type.chrge;
  pdbSerialNumber=serial;
  pdbocc=1.0;
  pdbfact=1.0;
  father=NULL;
  bonds=NULL;
  num_bonds=0;

}



/*Constructor.
  Copy of another atom. The only difference between both atoms is the name
  n: name of the new atom.
  a: atom to copy
*/
Atom::Atom(Atom *a)
{
  element=a->getElement();
  // strcpy(pdbName,a->getPdbName());
  memcpy(pdbName,(const char *)a->getPdbName(),5);

  Tcoor pos;
  a->getPosition(pos);
  position[0]=pos[0];
  position[1]=pos[1];
  position[2]=pos[2];

  charge=a->getCharge();
  pdbSerialNumber=a->getPdbSerial();
  pdbocc=a->getPdbocc();
  pdbfact=a->getPdbfact();
  father=NULL;
  bonds=NULL;
  num_bonds=0;

}

/*Destructor*/
Atom::~Atom()
{
	//std::cout<<"Atomo eliminado"<<std::endl;
	if(num_bonds>0)
		free(bonds);
}


/*Get the position of the atom*/
void Atom::getPosition(Tcoor coor)
{
  coor[0]=position[0];
  coor[1]=position[1];
  coor[2]=position[2];
}


/*Get the charge of the atom*/
float Atom::getCharge()
{
  return charge;
}

/*Get the object Element of the atom*/
Element* Atom::getElement()
{
  return element;
}


int Atom::getPdbSerial()
{
  return pdbSerialNumber;
}

char* Atom::getPdbName()
{
  return pdbName;
}

float Atom::getPdbocc()
{
  return pdbocc;
}

float Atom::getPdbfact()
{
  return pdbfact;
}

/*Set new position for the atom*/
void Atom::setPosition(Tcoor pos)
{
  position[0]=pos[0];
  position[1]=pos[1];
  position[2]=pos[2];

}

/*Set new Charge for the atom*/
void Atom::setCharge(float ch)
{
  charge=ch;
}

void Atom::setPdbSerial(int serial)
{
  pdbSerialNumber=serial;
}

void Atom::setPdbName(char n[5])
{
  strcpy(pdbName,n);
}

void Atom::setPdbocc(float occ)
{
  pdbocc=occ;
}


void Atom::setPdbfact(float fact)
{
  pdbfact=fact;
}

/*move an offset the position of the atom */
bool Atom::move(Tcoor offset)
{
  position[0]+=offset[0];
  position[1]+=offset[1];
  position[2]+=offset[2];
  return true;
}
/*move an offset the position of the atom */
bool Atom::moven(Tcoor offset)
{
  position[0]-=offset[0];
  position[1]-=offset[1];
  position[2]-=offset[2];
  return true;
}



char *Atom::getName()
{
  return pdbName;
}

bool Atom::initAll()
{
  return true;
}

bool Atom::moveAll(Tcoor offset)
{
  return move(offset);
}

bool Atom::moveAlln(Tcoor offset)
{
  return moven(offset);
}

TElement Atom::getClass()
{
  return pdb_atom;
}


TMOL Atom::getMolType()
{
	PDB_Container *f;
	f=(PDB_Container*)getFather();
	if(f!=NULL)
		return f->getMolType();
	else
		return tmol_null;
}

int Atom::get_numBonds()
{
	return num_bonds;
}
Bond* Atom::getBond(int i)
{
	if( i>num_bonds )
		return NULL;
	else
	{
		return bonds[i];
	}
}
void Atom::insertBond(Bond *b)
{
	num_bonds++;

	bonds=(Bond**)realloc(bonds,sizeof(Bond*)*num_bonds);
	bonds[num_bonds-1]=b;
}

bool Atom::removeBond(Bond *b)
{
	int i,j;

	for(i=0;i<num_bonds;i++)
	{
		if(b==bonds[i])
		{
			for(j=i+1;j<num_bonds;j++)
				bonds[j-1]=bonds[j];
			num_bonds--;
			bonds=(Bond**)realloc(bonds,sizeof(Bond*)*num_bonds);
			return true;
		}
	}
	return false;
}

Bond::Bond(Atom* i, Atom *f, int l)
{
	final=f;
	init=i;
	link=l;
}

Bond::~Bond()
{

}

int Bond::getLink()
{
	return link;
}

Atom * Bond::getInit()
{
	return init;
}

Atom *Bond::getFinal()
{
	return final;
}


