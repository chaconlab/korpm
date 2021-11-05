/*
 * pd2.h
 * SBG's implementation of the PD2 bump filter (MacDonald et al. 2013)
 *  Created on: Apr 29, 2015
 *      Author: mon
 */

#ifndef PD2_H_
#define PD2_H_

#include <libpdb/include/Macromolecule.h>
#include <libpdb/include/pdbIter.h>

// PD2 ENERGY STUFF...
#define myN  0
#define myCA 1
#define myC  2
#define myO  3
#define myCB 4

// "extern" required to prevent multiple definition because a new instance of the data structure is created each time this is included.
extern double pd2[5][5];

// #########################################################################
// FUNCTIONS DECLARATION ZONE
// #########################################################################

// The "PD2" bump energy calculation stuff...
//  coord --> Coordinates (1D array)
//  type --> Atom type array (size= number of atoms)
//  res --> Residue index array (size= number of atoms)
//  natoms --> Number of atoms
double pd2_bump(float *coord, char *type, int *res, int natoms);
// The "PD2" bump energy calculation stuff... for the Minimal Backbone Model (N, CA, and C atoms only)
//  coord --> Coordinates (1D array or the N, CA, and C atom coordinates)
//  type --> Atom type array (size= number of atoms)
//  res --> Residue index array (size= number of atoms)
//  natoms --> Number of atoms
double pd2_bumpNCAC(float *coord, char *type, int *res, int natoms);
// The "PD2" bump energy calculation stuff: "Pocket vs. Loop"
//  coord --> Pocket coordinates (1D array)
//  type --> Atom type array for pocket (size= pocket number of atoms)
//  res --> Residue index array for pocket (size= pocket number of atoms)
//  npocket --> Pocket number of atoms
//  co --> Loop coordinates for N, CA, C atoms
//  pco --> Loop coordinates for O and CB atoms
//  nloop --> Loop number of atoms (N, CA, C only!)
//  loopNt --> Loop first residue index
//  cg --> Coarse grainin level
//  residuemarker --> Array with the residue types. GLY=1 (size= loop number of atoms, N, CA, C only!)
double pd2_extrabump(float *coord, char *type, int *res, int npocket, double **co, double **pco, int nloop, int loopNt, int cg, int *residuemarker);
// The "PD2" bump energy calculation stuff: "Loop vs. Loop"
//  co --> Loop coordinates for N, CA, C atoms
//  pco --> Loop coordinates for O and CB atoms
//  nloop --> Loop number of atoms (N, CA, C only!)
//  cg --> Coarse graining level
//  residuemarker --> Array with the residue types. GLY=1 (size= loop number of atoms, N, CA, C only!)
double pd2_intrabump(double **co, double **pco, int nloop, int cg, int *residuemarker);
// Creates an array with the atomic type (N,CA,etc..) for the whole Macromolecule ("PD2" energy calculation stuff)
void pd2_type(Macromolecule *mol, char **p_type);
// Creates an array with the residue index of the whole Macromolecule ("PD2" energy calculation stuff)
void pd2_res(Macromolecule *mol, int **p_res);
// Returns the internal residue index (pos_fragment) corresponding to the input index (resnum PDB) ("PD2" energy calculation stuff)
//  resnum --> Residue index in PDB
int resindex(Macromolecule *mol, int resnum);


#endif /* PD2_H_ */
