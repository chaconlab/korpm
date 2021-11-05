/*
 * rama.h
 *
 *  Created on: Jan 4, 2019
 *      Author: mon
 */

#ifndef INCLUDE_RAMA_H_
#define INCLUDE_RAMA_H_

#include <stdio.h> // needed by some compilers
#include <stdlib.h>
#include <cmath>
//#include <libpdb/include/ResIni.h>
#include "libpdb/include/Macromolecule.h"
//#include "libga/include/chain_manip.h"

#define NAAS 21 // Number of aminoacids (20 + Cis-PRO) // TEST

// read H3 maps
float ***read_h3(char *file, int *h3size, float *h3step); // Read a Dunbrack map (generated with our "rang" program)


// Dunbrack stuff
float *****read_dunbrack(char *file, int *p_size, float *p_step); // Read a Dunbrack map (generated with our "rang" program)
// Returns the sequence indices array from the 1-letter code sequence (seq) for our Rama potential (automatic memory allocation)
//  seq --> Sequence in 1-letter code
//  nresl --> Number of residues
int *seq2iaa(char *seq, int nres);

// Allocate and zero initializes a size x size matrix
float **MatrixInit(int size);

// Element-wise accumulation of M2 (weighted or not) into M1 (size x size squared matrices)
// Optionally, M2 can be weighted by some factor "w"
void MatrixAccum(float **M1, float **M2, int size, float w = 1.0);

float **MatrixProduct(float **M1, float **M2, int size); // required by "map_gen"
float **MatrixDivision(float **M1, float **M2, int size); // required by "map_gen"
void norm_map(float **map, int size);
void show_map(float **map,int size, char *file);

float **map_gen(int C, int L, int R, float *****pdfs, int size, float **ref=NULL); // generate our own Dunbrack maps

// Generates a "Total" reference map from Dunbrack's tabulated PDFs (allocating output map memory).
//  pdfs --> Dunbrack's PDFs
//  size --> Number of side bins in Dunbrack's PDFs (size x size 2D maps)
//  weights --> Array of molar fractions per residue (if NULL all weights 1.0)
float **ref_gen(float *****pdfs, int size, float *weights = NULL);

// Avoid quantization problems when "angle" is -180 or +180 (-M_PI or +M_PI)
// 	angle = ]-2*M_PI,+2*M_PI[
// 	returned index is bounded to [0,size-1]
int angle2index(double angle, int size);

// Ramachandran energy derived from neighbor-dependant PDFs from Dunbrack's paper.
float rama_energy(double *dihedral_angle, int nr_atoms_loop,float ***maps, int size);

// Compute the dot product between two 3D vectors (float)
float dotprod(float *vector1, float *vector2);

// Normalize one 3D vector (float)
float norm(float *bond);

// Get the dihedrals array (dihedral) for supplied coordinates (co)
//  co --> floats array of atomic coordinates (the first non-zero dihedral defines the 4th atom position)
//  natoms --> Number of atoms of the kinematic chain
//  dihedral --> floats array of dihedral angles (the first 3 are zero by convention)
void finddihedral(float *co, int natoms, float *dihedral);

// Ramachandran energy derived from neighbor-dependent PDFs from Dunbrack's paper.
float rama_energy2(float *dihedral_angle, int nres, float ***maps, int size);
float rama_energy2(double *dihedral_angle, int nres, float ***maps, int size);

// Generate the Complete Ramachandran potential (20x20x20 = 8000 72x72 maps)
//  pdfs       --> Dunbrack's Neighbor dependent PDFs, as read by read_dunbrack()
//  rama_model --> Ramachandran energy model
//  size       --> Sampling of Dunbrack's PDFs (Optiona, 72 by default)
float *****gen_allrama(float *****pdfs, int rama_model, int size = 72);

// Set the cis-Pro identifier (p) into a "nres" residues sequence "seq" from its dihedral angles "dihedrals"
int setCisPro(float *dihedrals, char *seq, int nres);

#endif /* INCLUDE_RAMA_H_ */
