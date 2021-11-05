/*
 * korpm.h
 *
 * KORP's Mapping (otherwise duplicate definition issues occur)
 *
 *  Created on: Feb 16, 2018
 *      Author: mon
 */

#ifndef INCLUDE_KORPM_H_
#define INCLUDE_KORPM_H_

char aasAA[]={ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG, SER, THR, VAL, TRP, TYR};
char frasAA[]={ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG, SER, THR, VAL, TRP, TYR};

// OUR RESIDUE TYPES (from ResIni.h)
//	ALA =  0 ;
//	CYS =  1 ;
//	ASP =  2 ;
//	GLU =  3 ;
//	PHE =  4 ;
//	GLY =  5 ;
//	HIS =  6 ;
//	ILE =  7 ;
//	LYS =  8 ;
//	LEU =  9 ;
//	MET =  10 ;
//	ASN =  11 ;
//	PRO =  12 ;
//	GLN =  13 ;
//	ARG =  14 ;
//	SER =  15 ;
//	THR =  16 ;
//	VAL =  17 ;
//	TRP =  18 ;
//	TYR =  19 ;
//	ASH =  20 ;       //Neutral ASP
//	CYX =  21 ;       //SS-bonded CYS
//	CYM =  22 ;       //Negative CYS
//	GLH =  23 ;       //Neutral GLU
//	HIP =  24 ;       //Positive HIS
//	HID =  25 ;       //Neutral HIS, proton HD1 present
//	HIE =  26 ;       //Neutral HIS, proton HE2 present
//	LYN =  27 ;       //Neutral LYS
//	TYM =  28 ;       //Negative TYR
//	MSE =  29 ; 		// Seleno-Methionine
//	NtE =  30 ;
//	CtE =  31 ;
//	DGUA =  32 ;
//	DADE =  33 ;
//	DCYT =  34 ;
//	DTHY =  35 ;
//	GUA =  36 ;
//	ADE =  37 ;
//	CYT =  38 ;
//	URA =  39 ;

// Interaction frames mapping:
//   mapping[][0]--> Number of interaction frames for current residue (to screen all interactions)
//   mapping[][1]--> Interaction-specific potential map index (to map different residue interactions into the required potential map)
//   mapping[][2]--> "N" atom index of the interaction frame (for orthogonal framework definition)
//   mapping[][3]--> "CA" atom index of the interaction frame (for orthogonal framework definition and pairwise distance evaluation)
//   mapping[][4]--> "C" atom index of the interaction frame (for orthogonal framework definition)

// INTERACTION FRAMES: 20 standard aminoacids model (mappingAA) using CA for distance evaluation and N,C for interaction frame definition
char a01[] =  {1,  0, 0, 1, 2}; // ALA
char a02[] =  {1,  1, 0, 1, 2}; // CYS
char a03[] =  {1,  2, 0, 1, 2}; // ASP
char a04[] =  {1,  3, 0, 1, 2}; // GLU
char a05[] =  {1,  4, 0, 1, 2}; // PHE
char a06[] =  {1,  5, 0, 1, 2}; // GLY
char a07[] =  {1,  6, 0, 1, 2}; // HIS
char a08[] =  {1,  7, 0, 1, 2}; // ILE
char a09[] =  {1,  8, 0, 1, 2}; // LYS
char a10[] =  {1,  9, 0, 1, 2}; // LEU
char a11[] =  {1, 10, 0, 1, 2}; // MET
char a12[] =  {1, 11, 0, 1, 2}; // ASN
char a13[] =  {1, 12, 0, 1, 2}; // PRO
char a14[] =  {1, 13, 0, 1, 2}; // GLN
char a15[] =  {1, 14, 0, 1, 2}; // ARG
char a16[] =  {1, 15, 0, 1, 2}; // SER
char a17[] =  {1, 16, 0, 1, 2}; // THR
char a18[] =  {1, 17, 0, 1, 2}; // VAL
char a19[] =  {1, 18, 0, 1, 2}; // TRP
char a20[] =  {1, 19, 0, 1, 2}; // TYR
// Warning! This "chapa" is strictly required since in C you can't directly initialize an array from another. E.g.:
// int *myarray[] = { {3,4,5}, {32,54,6} }; <-- NOT ALLOWED IN C
char *mappingAA[] = { a01, a02, a03, a04, a05, a06, a07, a08, a09, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20 };

// INTERACTION FRAMES: 20 standard aminoacids model (mappingAA) using CB for distance evaluation and N,C for interaction frame definition
char c01[] =  {1,  0, 0, 4, 2}; // ALA
char c02[] =  {1,  1, 0, 4, 2}; // CYS
char c03[] =  {1,  2, 0, 4, 2}; // ASP
char c04[] =  {1,  3, 0, 4, 2}; // GLU
char c05[] =  {1,  4, 0, 4, 2}; // PHE
char c06[] =  {1,  5, 0, 1, 2}; // GLY (using CA instead of CB)
char c07[] =  {1,  6, 0, 4, 2}; // HIS
char c08[] =  {1,  7, 0, 4, 2}; // ILE
char c09[] =  {1,  8, 0, 4, 2}; // LYS
char c10[] =  {1,  9, 0, 4, 2}; // LEU
char c11[] =  {1, 10, 0, 4, 2}; // MET
char c12[] =  {1, 11, 0, 4, 2}; // ASN
char c13[] =  {1, 12, 0, 4, 2}; // PRO
char c14[] =  {1, 13, 0, 4, 2}; // GLN
char c15[] =  {1, 14, 0, 4, 2}; // ARG
char c16[] =  {1, 15, 0, 4, 2}; // SER
char c17[] =  {1, 16, 0, 4, 2}; // THR
char c18[] =  {1, 17, 0, 4, 2}; // VAL
char c19[] =  {1, 18, 0, 4, 2}; // TRP
char c20[] =  {1, 19, 0, 4, 2}; // TYR
char *mappingAACB[] = { c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20 };

// INTERACTION FRAMES: 20 standard aminoacids model (mappingAA) using O for distance evaluation and N,C for interaction frame definition
char d01[] =  {1,  0, 0, 3, 2}; // ALA
char d02[] =  {1,  1, 0, 3, 2}; // CYS
char d03[] =  {1,  2, 0, 3, 2}; // ASP
char d04[] =  {1,  3, 0, 3, 2}; // GLU
char d05[] =  {1,  4, 0, 3, 2}; // PHE
char d06[] =  {1,  5, 0, 3, 2}; // GLY
char d07[] =  {1,  6, 0, 3, 2}; // HIS
char d08[] =  {1,  7, 0, 3, 2}; // ILE
char d09[] =  {1,  8, 0, 3, 2}; // LYS
char d10[] =  {1,  9, 0, 3, 2}; // LEU
char d11[] =  {1, 10, 0, 3, 2}; // MET
char d12[] =  {1, 11, 0, 3, 2}; // ASN
char d13[] =  {1, 12, 0, 3, 2}; // PRO
char d14[] =  {1, 13, 0, 3, 2}; // GLN
char d15[] =  {1, 14, 0, 3, 2}; // ARG
char d16[] =  {1, 15, 0, 3, 2}; // SER
char d17[] =  {1, 16, 0, 3, 2}; // THR
char d18[] =  {1, 17, 0, 3, 2}; // VAL
char d19[] =  {1, 18, 0, 3, 2}; // TRP
char d20[] =  {1, 19, 0, 3, 2}; // TYR
char *mappingAAO[] = { d01, d02, d03, d04, d05, d06, d07, d08, d09, d10, d11, d12, d13, d14, d15, d16, d17, d18, d19, d20 };

// INTERACTION FRAMES: 20 standard aminoacids model (mappingAAC) using C for distance evaluation and N,CA for interaction frame definition
char f01[] =  {1,  0, 0, 2, 1}; // ALA
char f02[] =  {1,  1, 0, 2, 1}; // CYS
char f03[] =  {1,  2, 0, 2, 1}; // ASP
char f04[] =  {1,  3, 0, 2, 1}; // GLU
char f05[] =  {1,  4, 0, 2, 1}; // PHE
char f06[] =  {1,  5, 0, 2, 1}; // GLY
char f07[] =  {1,  6, 0, 2, 1}; // HIS
char f08[] =  {1,  7, 0, 2, 1}; // ILE
char f09[] =  {1,  8, 0, 2, 1}; // LYS
char f10[] =  {1,  9, 0, 2, 1}; // LEU
char f11[] =  {1, 10, 0, 2, 1}; // MET
char f12[] =  {1, 11, 0, 2, 1}; // ASN
char f13[] =  {1, 12, 0, 2, 1}; // PRO
char f14[] =  {1, 13, 0, 2, 1}; // GLN
char f15[] =  {1, 14, 0, 2, 1}; // ARG
char f16[] =  {1, 15, 0, 2, 1}; // SER
char f17[] =  {1, 16, 0, 2, 1}; // THR
char f18[] =  {1, 17, 0, 2, 1}; // VAL
char f19[] =  {1, 18, 0, 2, 1}; // TRP
char f20[] =  {1, 19, 0, 2, 1}; // TYR
char *mappingAAC[] = {f01, f02, f03, f04, f05, f06, f07, f08, f09, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20 };

// INTERACTION FRAMES: 20 standard aminoacids model (mappingAAC) using N for distance evaluation and CA,C for interaction frame definition
char b01[] =  {1,  0, 1, 0, 2}; // ALA
char b02[] =  {1,  1, 1, 0, 2}; // CYS
char b03[] =  {1,  2, 1, 0, 2}; // ASP
char b04[] =  {1,  3, 1, 0, 2}; // GLU
char b05[] =  {1,  4, 1, 0, 2}; // PHE
char b06[] =  {1,  5, 1, 0, 2}; // GLY
char b07[] =  {1,  6, 1, 0, 2}; // HIS
char b08[] =  {1,  7, 1, 0, 2}; // ILE
char b09[] =  {1,  8, 1, 0, 2}; // LYS
char b10[] =  {1,  9, 1, 0, 2}; // LEU
char b11[] =  {1, 10, 1, 0, 2}; // MET
char b12[] =  {1, 11, 1, 0, 2}; // ASN
char b13[] =  {1, 12, 1, 0, 2}; // PRO
char b14[] =  {1, 13, 1, 0, 2}; // GLN
char b15[] =  {1, 14, 1, 0, 2}; // ARG
char b16[] =  {1, 15, 1, 0, 2}; // SER
char b17[] =  {1, 16, 1, 0, 2}; // THR
char b18[] =  {1, 17, 1, 0, 2}; // VAL
char b19[] =  {1, 18, 1, 0, 2}; // TRP
char b20[] =  {1, 19, 1, 0, 2}; // TYR
char *mappingAAN[] = { b01, b02, b03, b04, b05, b06, b07, b08, b09, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20 };

// INTERACTION FRAMES:
char nftAP = 2; // Number of frame types
// 20 cores (CA for distances and N,C for frames) + 20 peptide groups (C for distances and CA,O for frames)
char e01[] =  {2,  0, 0, 1, 2,  0, 1, 2, 3}; // ALA
char e02[] =  {2,  1, 0, 1, 2,  1, 1, 2, 3}; // CYS
char e03[] =  {2,  2, 0, 1, 2,  2, 1, 2, 3}; // ASP
char e04[] =  {2,  3, 0, 1, 2,  3, 1, 2, 3}; // GLU
char e05[] =  {2,  4, 0, 1, 2,  4, 1, 2, 3}; // PHE
char e06[] =  {2,  5, 0, 1, 2,  5, 1, 2, 3}; // GLY
char e07[] =  {2,  6, 0, 1, 2,  6, 1, 2, 3}; // HIS
char e08[] =  {2,  7, 0, 1, 2,  7, 1, 2, 3}; // ILE
char e09[] =  {2,  8, 0, 1, 2,  8, 1, 2, 3}; // LYS
char e10[] =  {2,  9, 0, 1, 2,  9, 1, 2, 3}; // LEU
char e11[] =  {2, 10, 0, 1, 2, 10, 1, 2, 3}; // MET
char e12[] =  {2, 11, 0, 1, 2, 11, 1, 2, 3}; // ASN
char e13[] =  {2, 12, 0, 1, 2, 12, 1, 2, 3}; // PRO
char e14[] =  {2, 13, 0, 1, 2, 13, 1, 2, 3}; // GLN
char e15[] =  {2, 14, 0, 1, 2, 14, 1, 2, 3}; // ARG
char e16[] =  {2, 15, 0, 1, 2, 15, 1, 2, 3}; // SER
char e17[] =  {2, 16, 0, 1, 2, 16, 1, 2, 3}; // THR
char e18[] =  {2, 17, 0, 1, 2, 17, 1, 2, 3}; // VAL
char e19[] =  {2, 18, 0, 1, 2, 18, 1, 2, 3}; // TRP
char e20[] =  {2, 19, 0, 1, 2, 19, 1, 2, 3}; // TYR
char *mappingAP_CAC[] = { e01, e02, e03, e04, e05, e06, e07, e08, e09, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19, e20 };

// 20 cores (CA for distances and N,C for frames) + 20 peptide groups (O for distances and CA,C for frames)
char g01[] =  {2,  0, 0, 1, 2,  0, 1, 3, 2}; // ALA
char g02[] =  {2,  1, 0, 1, 2,  1, 1, 3, 2}; // CYS
char g03[] =  {2,  2, 0, 1, 2,  2, 1, 3, 2}; // ASP
char g04[] =  {2,  3, 0, 1, 2,  3, 1, 3, 2}; // GLU
char g05[] =  {2,  4, 0, 1, 2,  4, 1, 3, 2}; // PHE
char g06[] =  {2,  5, 0, 1, 2,  5, 1, 3, 2}; // GLY
char g07[] =  {2,  6, 0, 1, 2,  6, 1, 3, 2}; // HIS
char g08[] =  {2,  7, 0, 1, 2,  7, 1, 3, 2}; // ILE
char g09[] =  {2,  8, 0, 1, 2,  8, 1, 3, 2}; // LYS
char g10[] =  {2,  9, 0, 1, 2,  9, 1, 3, 2}; // LEU
char g11[] =  {2, 10, 0, 1, 2, 10, 1, 3, 2}; // MET
char g12[] =  {2, 11, 0, 1, 2, 11, 1, 3, 2}; // ASN
char g13[] =  {2, 12, 0, 1, 2, 12, 1, 3, 2}; // PRO
char g14[] =  {2, 13, 0, 1, 2, 13, 1, 3, 2}; // GLN
char g15[] =  {2, 14, 0, 1, 2, 14, 1, 3, 2}; // ARG
char g16[] =  {2, 15, 0, 1, 2, 15, 1, 3, 2}; // SER
char g17[] =  {2, 16, 0, 1, 2, 16, 1, 3, 2}; // THR
char g18[] =  {2, 17, 0, 1, 2, 17, 1, 3, 2}; // VAL
char g19[] =  {2, 18, 0, 1, 2, 18, 1, 3, 2}; // TRP
char g20[] =  {2, 19, 0, 1, 2, 19, 1, 3, 2}; // TYR
char *mappingAP_CAO[] = { g01, g02, g03, g04, g05, g06, g07, g08, g09, g10, g11, g12, g13, g14, g15, g16, g17, g18, g19, g20 };

// 20 cores (CB for distances and N,C for frames) + 20 peptide groups (O for distances and CA,C for frames)
char h01[] =  {2,  0, 0, 4, 2,  0, 1, 3, 2}; // ALA
char h02[] =  {2,  1, 0, 4, 2,  1, 1, 3, 2}; // CYS
char h03[] =  {2,  2, 0, 4, 2,  2, 1, 3, 2}; // ASP
char h04[] =  {2,  3, 0, 4, 2,  3, 1, 3, 2}; // GLU
char h05[] =  {2,  4, 0, 4, 2,  4, 1, 3, 2}; // PHE
char h06[] =  {2,  5, 0, 1, 2,  5, 1, 3, 2}; // GLY (CA used instead CB)
char h07[] =  {2,  6, 0, 4, 2,  6, 1, 3, 2}; // HIS
char h08[] =  {2,  7, 0, 4, 2,  7, 1, 3, 2}; // ILE
char h09[] =  {2,  8, 0, 4, 2,  8, 1, 3, 2}; // LYS
char h10[] =  {2,  9, 0, 4, 2,  9, 1, 3, 2}; // LEU
char h11[] =  {2, 10, 0, 4, 2, 10, 1, 3, 2}; // MET
char h12[] =  {2, 11, 0, 4, 2, 11, 1, 3, 2}; // ASN
char h13[] =  {2, 12, 0, 4, 2, 12, 1, 3, 2}; // PRO
char h14[] =  {2, 13, 0, 4, 2, 13, 1, 3, 2}; // GLN
char h15[] =  {2, 14, 0, 4, 2, 14, 1, 3, 2}; // ARG
char h16[] =  {2, 15, 0, 4, 2, 15, 1, 3, 2}; // SER
char h17[] =  {2, 16, 0, 4, 2, 16, 1, 3, 2}; // THR
char h18[] =  {2, 17, 0, 4, 2, 17, 1, 3, 2}; // VAL
char h19[] =  {2, 18, 0, 4, 2, 18, 1, 3, 2}; // TRP
char h20[] =  {2, 19, 0, 4, 2, 19, 1, 3, 2}; // TYR
char *mappingAP_CBO[] = { h01, h02, h03, h04, h05, h06, h07, h08, h09, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20 };

// 20 Peptide-frames (O for distances and CA,C for frames) +
// 20 Nitrogen-frames (N for distances and CA,CB for frames) +
char p01[] =  {2,  0, 1, 3, 2,  0, 1, 0, 4}; // ALA
char p02[] =  {2,  1, 1, 3, 2,  1, 1, 0, 4}; // CYS
char p03[] =  {2,  2, 1, 3, 2,  2, 1, 0, 4}; // ASP
char p04[] =  {2,  3, 1, 3, 2,  3, 1, 0, 4}; // GLU
char p05[] =  {2,  4, 1, 3, 2,  4, 1, 0, 4}; // PHE
char p06[] =  {1,  5, 1, 3, 2};                           // GLY (no CB)
char p07[] =  {2,  6, 1, 3, 2,  6, 1, 0, 4}; // HIS
char p08[] =  {2,  7, 1, 3, 2,  7, 1, 0, 4}; // ILE
char p09[] =  {2,  8, 1, 3, 2,  8, 1, 0, 4}; // LYS
char p10[] =  {2,  9, 1, 3, 2,  9, 1, 0, 4}; // LEU
char p11[] =  {2, 10, 1, 3, 2, 10, 1, 0, 4}; // MET
char p12[] =  {2, 11, 1, 3, 2, 11, 1, 0, 4}; // ASN
char p13[] =  {2, 12, 1, 3, 2, 12, 1, 0, 4}; // PRO
char p14[] =  {2, 13, 1, 3, 2, 13, 1, 0, 4}; // GLN
char p15[] =  {2, 14, 1, 3, 2, 14, 1, 0, 4}; // ARG
char p16[] =  {2, 15, 1, 3, 2, 15, 1, 0, 4}; // SER
char p17[] =  {2, 16, 1, 3, 2, 16, 1, 0, 4}; // THR
char p18[] =  {2, 17, 1, 3, 2, 17, 1, 0, 4}; // VAL
char p19[] =  {2, 18, 1, 3, 2, 18, 1, 0, 4}; // TRP
char p20[] =  {2, 19, 1, 3, 2, 19, 1, 0, 4}; // TYR
char *mappingAP_ON[] = { p01, p02, p03, p04, p05, p06, p07, p08, p09, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20 };

// 20 Peptide-frames (O for distances and CA,C for frames) +
// 20 Nitrogen-frames (N for distances and CA,CB for frames) +
// 20 Residue-frames (CB for distances and N,CA for frames)
char q01[] =  {3,  0, 1, 3, 2,  0, 1, 0, 4,  0, 0, 4, 1}; // ALA
char q02[] =  {3,  1, 1, 3, 2,  1, 1, 0, 4,  1, 0, 4, 1}; // CYS
char q03[] =  {3,  2, 1, 3, 2,  2, 1, 0, 4,  2, 0, 4, 1}; // ASP
char q04[] =  {3,  3, 1, 3, 2,  3, 1, 0, 4,  3, 0, 4, 1}; // GLU
char q05[] =  {3,  4, 1, 3, 2,  4, 1, 0, 4,  4, 0, 4, 1}; // PHE
char q06[] =  {1,  5, 1, 3, 2};                           // GLY (no CB)
char q07[] =  {3,  6, 1, 3, 2,  6, 1, 0, 4,  6, 0, 4, 1}; // HIS
char q08[] =  {3,  7, 1, 3, 2,  7, 1, 0, 4,  7, 0, 4, 1}; // ILE
char q09[] =  {3,  8, 1, 3, 2,  8, 1, 0, 4,  8, 0, 4, 1}; // LYS
char q10[] =  {3,  9, 1, 3, 2,  9, 1, 0, 4,  9, 0, 4, 1}; // LEU
char q11[] =  {3, 10, 1, 3, 2, 10, 1, 0, 4, 10, 0, 4, 1}; // MET
char q12[] =  {3, 11, 1, 3, 2, 11, 1, 0, 4, 11, 0, 4, 1}; // ASN
char q13[] =  {3, 12, 1, 3, 2, 12, 1, 0, 4, 12, 0, 4, 1}; // PRO
char q14[] =  {3, 13, 1, 3, 2, 13, 1, 0, 4, 13, 0, 4, 1}; // GLN
char q15[] =  {3, 14, 1, 3, 2, 14, 1, 0, 4, 14, 0, 4, 1}; // ARG
char q16[] =  {3, 15, 1, 3, 2, 15, 1, 0, 4, 15, 0, 4, 1}; // SER
char q17[] =  {3, 16, 1, 3, 2, 16, 1, 0, 4, 16, 0, 4, 1}; // THR
char q18[] =  {3, 17, 1, 3, 2, 17, 1, 0, 4, 17, 0, 4, 1}; // VAL
char q19[] =  {3, 18, 1, 3, 2, 18, 1, 0, 4, 18, 0, 4, 1}; // TRP
char q20[] =  {3, 19, 1, 3, 2, 19, 1, 0, 4, 19, 0, 4, 1}; // TYR
char *mappingAP_ONCB[] = { q01, q02, q03, q04, q05, q06, q07, q08, q09, q10, q11, q12, q13, q14, q15, q16, q17, q18, q19, q20 };

// INTERACTION FRAMES:
char nft5F = 5; // Number of frame types
// Distance:             CA           C            O            N            CB
// Frame:             N     C      CA    O      CA    C      CA    CB     N     C
// 1:CA/N^C, 2:C/CA^O, 4:O/CA^C 4:N/CA^C, and 5:CB/N^C
char i01[] =  {5,  0, 0, 1, 2,  0, 1, 2, 3,  0, 1, 3, 2,  0, 1, 0, 4,  0, 0, 4, 2}; // ALA
char i02[] =  {5,  1, 0, 1, 2,  1, 1, 2, 3,  1, 1, 3, 2,  1, 1, 0, 4,  1, 0, 4, 2}; // CYS
char i03[] =  {5,  2, 0, 1, 2,  2, 1, 2, 3,  2, 1, 3, 2,  2, 1, 0, 4,  2, 0, 4, 2}; // ASP
char i04[] =  {5,  3, 0, 1, 2,  3, 1, 2, 3,  3, 1, 3, 2,  3, 1, 0, 4,  3, 0, 4, 2}; // GLU
char i05[] =  {5,  4, 0, 1, 2,  4, 1, 2, 3,  4, 1, 3, 2,  4, 1, 0, 4,  4, 0, 4, 2}; // PHE
char i06[] =  {3,  5, 0, 1, 2,  5, 1, 2, 3,  5, 1, 3, 2};                           // GLY (no CB)
char i07[] =  {5,  6, 0, 1, 2,  6, 1, 2, 3,  6, 1, 3, 2,  6, 1, 0, 4,  6, 0, 4, 2}; // HIS
char i08[] =  {5,  7, 0, 1, 2,  7, 1, 2, 3,  7, 1, 3, 2,  7, 1, 0, 4,  7, 0, 4, 2}; // ILE
char i09[] =  {5,  8, 0, 1, 2,  8, 1, 2, 3,  8, 1, 3, 2,  8, 1, 0, 4,  8, 0, 4, 2}; // LYS
char i10[] =  {5,  9, 0, 1, 2,  9, 1, 2, 3,  9, 1, 3, 2,  9, 1, 0, 4,  9, 0, 4, 2}; // LEU
char i11[] =  {5, 10, 0, 1, 2, 10, 1, 2, 3, 10, 1, 3, 2, 10, 1, 0, 4, 10, 0, 4, 2}; // MET
char i12[] =  {5, 11, 0, 1, 2, 11, 1, 2, 3, 11, 1, 3, 2, 11, 1, 0, 4, 11, 0, 4, 2}; // ASN
char i13[] =  {5, 12, 0, 1, 2, 12, 1, 2, 3, 12, 1, 3, 2, 12, 1, 0, 4, 12, 0, 4, 2}; // PRO
char i14[] =  {5, 13, 0, 1, 2, 13, 1, 2, 3, 13, 1, 3, 2, 13, 1, 0, 4, 13, 0, 4, 2}; // GLN
char i15[] =  {5, 14, 0, 1, 2, 14, 1, 2, 3, 14, 1, 3, 2, 14, 1, 0, 4, 14, 0, 4, 2}; // ARG
char i16[] =  {5, 15, 0, 1, 2, 15, 1, 2, 3, 15, 1, 3, 2, 15, 1, 0, 4, 15, 0, 4, 2}; // SER
char i17[] =  {5, 16, 0, 1, 2, 16, 1, 2, 3, 16, 1, 3, 2, 16, 1, 0, 4, 16, 0, 4, 2}; // THR
char i18[] =  {5, 17, 0, 1, 2, 17, 1, 2, 3, 17, 1, 3, 2, 17, 1, 0, 4, 17, 0, 4, 2}; // VAL
char i19[] =  {5, 18, 0, 1, 2, 18, 1, 2, 3, 18, 1, 3, 2, 18, 1, 0, 4, 18, 0, 4, 2}; // TRP
char i20[] =  {5, 19, 0, 1, 2, 19, 1, 2, 3, 19, 1, 3, 2, 19, 1, 0, 4, 19, 0, 4, 2}; // TYR
char *mapping5F[] = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19, i20 };

// INTERACTION FRAMES for 4D:
char frasOA[]={ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG, SER, THR, VAL, TRP, TYR};
// See Table S1 form Order_AVE's paper for further info. (3rd atom is neglected in Order_AVE, only 1st and 2nd are used)
// Distance:             N            CA           C            O            CB
// Frame:             CA    C      N     C      CA    O      CA    C      N     C
char j01[] =  {5,  0,-2, 0, 2,  0, 0, 1, 2,  0, 1, 2, 3,  0, 2, 3, 2,  0, 1, 4, 2}; // ALA
char j02[] =  {5,  1,-2, 0, 2,  1, 0, 1, 2,  1, 1, 2, 3,  1, 2, 3, 2,  1, 1, 4, 2}; // CYS
char j03[] =  {5,  2,-2, 0, 2,  2, 0, 1, 2,  2, 1, 2, 3,  2, 2, 3, 2,  2, 1, 4, 2}; // ASP
char j04[] =  {5,  3,-2, 0, 2,  3, 0, 1, 2,  3, 1, 2, 3,  3, 2, 3, 2,  3, 1, 4, 2}; // GLU
char j05[] =  {5,  4,-2, 0, 2,  4, 0, 1, 2,  4, 1, 2, 3,  4, 2, 3, 2,  4, 1, 4, 2}; // PHE
char j06[] =  {4,  5,-2, 0, 2,  5, 0, 1, 2,  5, 1, 2, 3,  5, 2, 3, 2};              // GLY (CB term neglected)
char j07[] =  {5,  6,-2, 0, 2,  6, 0, 1, 2,  6, 1, 2, 3,  6, 2, 3, 2,  6, 1, 4, 2}; // HIS
char j08[] =  {5,  7,-2, 0, 2,  7, 0, 1, 2,  7, 1, 2, 3,  7, 2, 3, 2,  7, 1, 4, 2}; // ILE
char j09[] =  {5,  8,-2, 0, 2,  8, 0, 1, 2,  8, 1, 2, 3,  8, 2, 3, 2,  8, 1, 4, 2}; // LYS
char j10[] =  {5,  9,-2, 0, 2,  9, 0, 1, 2,  9, 1, 2, 3,  9, 2, 3, 2,  9, 1, 4, 2}; // LEU
char j11[] =  {5, 10,-2, 0, 2, 10, 0, 1, 2, 10, 1, 2, 3, 10, 2, 3, 2, 10, 1, 4, 2}; // MET
char j12[] =  {5, 11,-2, 0, 2, 11, 0, 1, 2, 11, 1, 2, 3, 11, 2, 3, 2, 11, 1, 4, 2}; // ASN
char j13[] =  {5, 12,-2, 0, 2, 12, 0, 1, 2, 12, 1, 2, 3, 12, 2, 3, 2, 12, 1, 4, 2}; // PRO
char j14[] =  {5, 13,-2, 0, 2, 13, 0, 1, 2, 13, 1, 2, 3, 13, 2, 3, 2, 13, 1, 4, 2}; // GLN
char j15[] =  {5, 14,-2, 0, 2, 14, 0, 1, 2, 14, 1, 2, 3, 14, 2, 3, 2, 14, 1, 4, 2}; // ARG
char j16[] =  {5, 15,-2, 0, 2, 15, 0, 1, 2, 15, 1, 2, 3, 15, 2, 3, 2, 15, 1, 4, 2}; // SER
char j17[] =  {5, 16,-2, 0, 2, 16, 0, 1, 2, 16, 1, 2, 3, 16, 2, 3, 2, 16, 1, 4, 2}; // THR
char j18[] =  {5, 17,-2, 0, 2, 17, 0, 1, 2, 17, 1, 2, 3, 17, 2, 3, 2, 17, 1, 4, 2}; // VAL
char j19[] =  {5, 18,-2, 0, 2, 18, 0, 1, 2, 18, 1, 2, 3, 18, 2, 3, 2, 18, 1, 4, 2}; // TRP
char j20[] =  {5, 19,-2, 0, 2, 19, 0, 1, 2, 19, 1, 2, 3, 19, 2, 3, 2, 19, 1, 4, 2}; // TYR
char *mappingOA[] = { j01, j02, j03, j04, j05, j06, j07, j08, j09, j10, j11, j12, j13, j14, j15, j16, j17, j18, j19, j20 };

// Frame type: 54 ---> Order_AVE's N,CA,C,O frames (All but CB frame)
char nft4Df54 = 4; // Number of frame types for Order_AVE's N,CA,C,O frames
// See Table S1 form Order_AVE's paper for further info. (3rd atom is neglected in Order_AVE, only 1st and 2nd are used)
// Distance:             N            CA           C            O
// Frame:      #      CA    C      N     C      CA    O      CA    C
char k01[] =  {4,  0,-2, 0, 2,  0, 0, 1, 2,  0, 1, 2, 3,  0, 2, 3, 2}; // ALA
char k02[] =  {4,  1,-2, 0, 2,  1, 0, 1, 2,  1, 1, 2, 3,  1, 2, 3, 2}; // CYS
char k03[] =  {4,  2,-2, 0, 2,  2, 0, 1, 2,  2, 1, 2, 3,  2, 2, 3, 2}; // ASP
char k04[] =  {4,  3,-2, 0, 2,  3, 0, 1, 2,  3, 1, 2, 3,  3, 2, 3, 2}; // GLU
char k05[] =  {4,  4,-2, 0, 2,  4, 0, 1, 2,  4, 1, 2, 3,  4, 2, 3, 2}; // PHE
char k06[] =  {4,  5,-2, 0, 2,  5, 0, 1, 2,  5, 1, 2, 3,  5, 2, 3, 2}; // GLY
char k07[] =  {4,  6,-2, 0, 2,  6, 0, 1, 2,  6, 1, 2, 3,  6, 2, 3, 2}; // HIS
char k08[] =  {4,  7,-2, 0, 2,  7, 0, 1, 2,  7, 1, 2, 3,  7, 2, 3, 2}; // ILE
char k09[] =  {4,  8,-2, 0, 2,  8, 0, 1, 2,  8, 1, 2, 3,  8, 2, 3, 2}; // LYS
char k10[] =  {4,  9,-2, 0, 2,  9, 0, 1, 2,  9, 1, 2, 3,  9, 2, 3, 2}; // LEU
char k11[] =  {4, 10,-2, 0, 2, 10, 0, 1, 2, 10, 1, 2, 3, 10, 2, 3, 2}; // MET
char k12[] =  {4, 11,-2, 0, 2, 11, 0, 1, 2, 11, 1, 2, 3, 11, 2, 3, 2}; // ASN
char k13[] =  {4, 12,-2, 0, 2, 12, 0, 1, 2, 12, 1, 2, 3, 12, 2, 3, 2}; // PRO
char k14[] =  {4, 13,-2, 0, 2, 13, 0, 1, 2, 13, 1, 2, 3, 13, 2, 3, 2}; // GLN
char k15[] =  {4, 14,-2, 0, 2, 14, 0, 1, 2, 14, 1, 2, 3, 14, 2, 3, 2}; // ARG
char k16[] =  {4, 15,-2, 0, 2, 15, 0, 1, 2, 15, 1, 2, 3, 15, 2, 3, 2}; // SER
char k17[] =  {4, 16,-2, 0, 2, 16, 0, 1, 2, 16, 1, 2, 3, 16, 2, 3, 2}; // THR
char k18[] =  {4, 17,-2, 0, 2, 17, 0, 1, 2, 17, 1, 2, 3, 17, 2, 3, 2}; // VAL
char k19[] =  {4, 18,-2, 0, 2, 18, 0, 1, 2, 18, 1, 2, 3, 18, 2, 3, 2}; // TRP
char k20[] =  {4, 19,-2, 0, 2, 19, 0, 1, 2, 19, 1, 2, 3, 19, 2, 3, 2}; // TYR
char *mapping4Df54[] = { k01, k02, k03, k04, k05, k06, k07, k08, k09, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20 };

// Frame type: 52 ---> Order_AVE's N,O frames (only H-bonding donor/acceptor frames)
char nft4Df52 = 2; // Number of frame types for Order_AVE's N,O frames
// See Table S1 form Order_AVE's paper for further info. (3rd atom is neglected in Order_AVE, only 1st and 2nd are used)
// Distance:             N            O
// Frame:      #      CA    C      CA    C
char l01[] =  {2,  0,-2, 0, 2,  0, 2, 3, 2}; // ALA
char l02[] =  {2,  1,-2, 0, 2,  1, 2, 3, 2}; // CYS
char l03[] =  {2,  2,-2, 0, 2,  2, 2, 3, 2}; // ASP
char l04[] =  {2,  3,-2, 0, 2,  3, 2, 3, 2}; // GLU
char l05[] =  {2,  4,-2, 0, 2,  4, 2, 3, 2}; // PHE
char l06[] =  {2,  5,-2, 0, 2,  5, 2, 3, 2}; // GLY
char l07[] =  {2,  6,-2, 0, 2,  6, 2, 3, 2}; // HIS
char l08[] =  {2,  7,-2, 0, 2,  7, 2, 3, 2}; // ILE
char l09[] =  {2,  8,-2, 0, 2,  8, 2, 3, 2}; // LYS
char l10[] =  {2,  9,-2, 0, 2,  9, 2, 3, 2}; // LEU
char l11[] =  {2, 10,-2, 0, 2, 10, 2, 3, 2}; // MET
char l12[] =  {2, 11,-2, 0, 2, 11, 2, 3, 2}; // ASN
char l13[] =  {2, 12,-2, 0, 2, 12, 2, 3, 2}; // PRO
char l14[] =  {2, 13,-2, 0, 2, 13, 2, 3, 2}; // GLN
char l15[] =  {2, 14,-2, 0, 2, 14, 2, 3, 2}; // ARG
char l16[] =  {2, 15,-2, 0, 2, 15, 2, 3, 2}; // SER
char l17[] =  {2, 16,-2, 0, 2, 16, 2, 3, 2}; // THR
char l18[] =  {2, 17,-2, 0, 2, 17, 2, 3, 2}; // VAL
char l19[] =  {2, 18,-2, 0, 2, 18, 2, 3, 2}; // TRP
char l20[] =  {2, 19,-2, 0, 2, 19, 2, 3, 2}; // TYR
char *mapping4Df52[] = { l01, l02, l03, l04, l05, l06, l07, l08, l09, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 };

// Frame type: 50 ---> Order_AVE's CA frame (minimalist 4D single frame)
char nft4Df50 = 1; // Number of frame types for Order_AVE's CA frames
// See Table S1 form Order_AVE's paper for further info. (3rd atom is neglected in Order_AVE, only 1st and 2nd are used)
// Distance:             CA
// Frame:      #      N     C
char m01[] =  {1,  0, 0, 1, 2}; // ALA
char m02[] =  {1,  1, 0, 1, 2}; // CYS
char m03[] =  {1,  2, 0, 1, 2}; // ASP
char m04[] =  {1,  3, 0, 1, 2}; // GLU
char m05[] =  {1,  4, 0, 1, 2}; // PHE
char m06[] =  {1,  5, 0, 1, 2}; // GLY
char m07[] =  {1,  6, 0, 1, 2}; // HIS
char m08[] =  {1,  7, 0, 1, 2}; // ILE
char m09[] =  {1,  8, 0, 1, 2}; // LYS
char m10[] =  {1,  9, 0, 1, 2}; // LEU
char m11[] =  {1, 10, 0, 1, 2}; // MET
char m12[] =  {1, 11, 0, 1, 2}; // ASN
char m13[] =  {1, 12, 0, 1, 2}; // PRO
char m14[] =  {1, 13, 0, 1, 2}; // GLN
char m15[] =  {1, 14, 0, 1, 2}; // ARG
char m16[] =  {1, 15, 0, 1, 2}; // SER
char m17[] =  {1, 16, 0, 1, 2}; // THR
char m18[] =  {1, 17, 0, 1, 2}; // VAL
char m19[] =  {1, 18, 0, 1, 2}; // TRP
char m20[] =  {1, 19, 0, 1, 2}; // TYR
char *mapping4Df50[] = { m01, m02, m03, m04, m05, m06, m07, m08, m09, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20 };

// Frame type: 53 ---> Order_AVE's N,O,CB frames (the most important frames?)
char nft4Df53 = 3; // Number of frame types for Order_AVE's N,O,CB frames
// See Table S1 form Order_AVE's paper for further info. (3rd atom is neglected in Order_AVE, only 1st and 2nd are used)
// Distance:             N            O            CB
// Frame:             CA    C      CA    C      N     C
char n01[] =  {3,  0,-2, 0, 2,  0, 2, 3, 2,  0, 1, 4, 2}; // ALA
char n02[] =  {3,  1,-2, 0, 2,  1, 2, 3, 2,  0, 1, 4, 2}; // CYS
char n03[] =  {3,  2,-2, 0, 2,  2, 2, 3, 2,  2, 1, 4, 2}; // ASP
char n04[] =  {3,  3,-2, 0, 2,  3, 2, 3, 2,  3, 1, 4, 2}; // GLU
char n05[] =  {3,  4,-2, 0, 2,  4, 2, 3, 2,  4, 1, 4, 2}; // PHE
char n06[] =  {2,  5,-2, 0, 2,  5, 2, 3, 2};              // GLY (CB term neglected)
char n07[] =  {3,  6,-2, 0, 2,  6, 2, 3, 2,  6, 1, 4, 2}; // HIS
char n08[] =  {3,  7,-2, 0, 2,  7, 2, 3, 2,  7, 1, 4, 2}; // ILE
char n09[] =  {3,  8,-2, 0, 2,  8, 2, 3, 2,  8, 1, 4, 2}; // LYS
char n10[] =  {3,  9,-2, 0, 2,  9, 2, 3, 2,  9, 1, 4, 2}; // LEU
char n11[] =  {3, 10,-2, 0, 2, 10, 2, 3, 2, 10, 1, 4, 2}; // MET
char n12[] =  {3, 11,-2, 0, 2, 11, 2, 3, 2, 11, 1, 4, 2}; // ASN
char n13[] =  {3, 12,-2, 0, 2, 12, 2, 3, 2, 12, 1, 4, 2}; // PRO
char n14[] =  {3, 13,-2, 0, 2, 13, 2, 3, 2, 13, 1, 4, 2}; // GLN
char n15[] =  {3, 14,-2, 0, 2, 14, 2, 3, 2, 14, 1, 4, 2}; // ARG
char n16[] =  {3, 15,-2, 0, 2, 15, 2, 3, 2, 15, 1, 4, 2}; // SER
char n17[] =  {3, 16,-2, 0, 2, 16, 2, 3, 2, 16, 1, 4, 2}; // THR
char n18[] =  {3, 17,-2, 0, 2, 17, 2, 3, 2, 17, 1, 4, 2}; // VAL
char n19[] =  {3, 18,-2, 0, 2, 18, 2, 3, 2, 18, 1, 4, 2}; // TRP
char n20[] =  {3, 19,-2, 0, 2, 19, 2, 3, 2, 19, 1, 4, 2}; // TYR
char *mapping4Df53[] = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20 };

// Frame type: 56 ---> Order_AVE's N*,O,CB frames (the most important frames with alternative N* frame definition)
char nft4Df56 = 3; // Number of frame types for Order_AVE's N*,O,CB frames
// See Table S1 form Order_AVE's paper for further info. (3rd atom is neglected in Order_AVE, only 1st and 2nd are used)
// Distance:             N            O            CB
// Frame:             CA    C      CA    C      N     C
char o01[] =  {3,  0, 1, 0, 2,  0, 2, 3, 2,  0, 1, 4, 2}; // ALA
char o02[] =  {3,  1, 1, 0, 2,  1, 2, 3, 2,  0, 1, 4, 2}; // CYS
char o03[] =  {3,  2, 1, 0, 2,  2, 2, 3, 2,  2, 1, 4, 2}; // ASP
char o04[] =  {3,  3, 1, 0, 2,  3, 2, 3, 2,  3, 1, 4, 2}; // GLU
char o05[] =  {3,  4, 1, 0, 2,  4, 2, 3, 2,  4, 1, 4, 2}; // PHE
char o06[] =  {2,  5, 1, 0, 2,  5, 2, 3, 2};              // GLY (CB term neglected)
char o07[] =  {3,  6, 1, 0, 2,  6, 2, 3, 2,  6, 1, 4, 2}; // HIS
char o08[] =  {3,  7, 1, 0, 2,  7, 2, 3, 2,  7, 1, 4, 2}; // ILE
char o09[] =  {3,  8, 1, 0, 2,  8, 2, 3, 2,  8, 1, 4, 2}; // LYS
char o10[] =  {3,  9, 1, 0, 2,  9, 2, 3, 2,  9, 1, 4, 2}; // LEU
char o11[] =  {3, 10, 1, 0, 2, 10, 2, 3, 2, 10, 1, 4, 2}; // MET
char o12[] =  {3, 11, 1, 0, 2, 11, 2, 3, 2, 11, 1, 4, 2}; // ASN
char o13[] =  {3, 12, 1, 0, 2, 12, 2, 3, 2, 12, 1, 4, 2}; // PRO
char o14[] =  {3, 13, 1, 0, 2, 13, 2, 3, 2, 13, 1, 4, 2}; // GLN
char o15[] =  {3, 14, 1, 0, 2, 14, 2, 3, 2, 14, 1, 4, 2}; // ARG
char o16[] =  {3, 15, 1, 0, 2, 15, 2, 3, 2, 15, 1, 4, 2}; // SER
char o17[] =  {3, 16, 1, 0, 2, 16, 2, 3, 2, 16, 1, 4, 2}; // THR
char o18[] =  {3, 17, 1, 0, 2, 17, 2, 3, 2, 17, 1, 4, 2}; // VAL
char o19[] =  {3, 18, 1, 0, 2, 18, 2, 3, 2, 18, 1, 4, 2}; // TRP
char o20[] =  {3, 19, 1, 0, 2, 19, 2, 3, 2, 19, 1, 4, 2}; // TYR
char *mapping4Df56[] = { o01, o02, o03, o04, o05, o06, o07, o08, o09, o10, o11, o12, o13, o14, o15, o16, o17, o18, o19, o20 };




#endif /* INCLUDE_KORPM_H_ */
