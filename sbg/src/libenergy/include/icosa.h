/*
 * locoicos.h
 *
 *  Created on: Apr 29, 2015
 *      Author: mon
 */

#ifndef ENERGY_H_
#define ENERGY_H_

#define PI 3.1415926535897932384626433832795

// Definitions for Loco_icos energy model
#define NAASloco 20 // Number of aminoacids
#define NTRI 20 // Number of triangular faces for an icosahedron
#define NVERT 12 // Number of vertices for an icosahedron
#define NDIS 12 // Number of distance bins
#define CAdist 2.8 // Minimal CA-CA distance
#define GAP 6 // Gap in sequence to consider interactions

// Distance
// #define DIST(I,J)	(sqrtf(powf(I[0]-J[0],2)+powf(I[1]-J[1],2)+powf(I[2]-J[2],2)))
#define DIST(I,J)	(sqrt(pow(I[0]-J[0],2)+pow(I[1]-J[1],2)+pow(I[2]-J[2],2)))
// Distance^2
//#define DIST2(I,J)	(powf(I[0]-J[0],2)+powf(I[1]-J[1],2)+powf(I[2]-J[2],2))
#define DIST2(I,J)	(pow(I[0]-J[0],2)+pow(I[1]-J[1],2)+pow(I[2]-J[2],2))

// Cross product
#define CROSS(X,Y,V)	{									\
							V[0] = X[1]*Y[2]-X[2]*Y[1];		\
							V[1] = X[2]*Y[0]-X[0]*Y[2];		\
							V[2] = X[0]*Y[1]-X[1]*Y[0];		\
						}

// Create a vector V
#define VECTOR(X,Y,V)	{						\
							V[0] = X[0]-Y[0];	\
							V[1] = X[1]-Y[1];	\
							V[2] = X[2]-Y[2];	\
						}

// Dot product X * Y
#define DOT(X,Y)	(X[0]*Y[0] + X[1]*Y[1] + X[2]*Y[2])

// Reverse sign of X into V
#define REVERSE(X,V)	{ V[0] = -X[0]; V[1] = -X[1]; V[2] = -X[2]; }

// "extern" required to prevent multiple definition because a new instance of the data structure is created each time this is included.
extern char aa[];
extern int icosTri[][3];

// Read loco energy from text file
void read_loco(double ****loco, char *name);
// Return aminoacid index given its 1-letter code character.
int aa2index(char aa, char *list);
// Initialize icosahedron data
void init_icos(double icosVertices[][3]);
// The "loco" energy calculation for single macromolecule
double loco_energy(float *xyz, int *iseq, int nres, double ****loco, double icosVertices[][3]);
// The "loco" energy calculation for protein-protein docking
double loco_energy(float *xyz, int *iseq, int nres, float *xyzL, int *iseqL, int nresL, double ****loco, double icosVertices[][3]);
// The "loco" energy calculation stuff...
// lindex --> residue index where the loop has been extracted from (in "coord" array), i.e. index of the first right-side residue.
double loco_energy(float *coordl, int *iseql, int nresl, float *coord, int *iseq, int nres, int lindex, double ****loco, double icosVertices[][3]);
// The "loco" energy calculation stuff...
double loco_energy(float *xyz, int *iseq, int nres, int nresl, int lindex, double ****loco, double icosVertices[][3]);
// Compute CA position in local coordinate frame
void getcoor(float *Nref, float *CAref, float *Cref, float *contactCAref, double *cen);
// Determine if a ray going through a triangle
bool rayIntersectsTriangle(double *dref, double *v0, double *v1, double *v2);
// Convert "co->el" 2D-matrix (Pieter's) into linear array (coordMatrix-like)
//  el --> 2D-matrix pointer
//  coord --> linear array
//  natoms --> number of atoms to copy coordinates
void co2coord(double **el, float *coord, int natoms);
// Convert a 1-letter code sequence into integer LOCO code sequence
//  seq --> 1-letter code sequence array
//  p_iseq --> Pointer to the array of integers with LOCO code sequence (if NULL, it allocates memory)
//  nres --> Number of residues
void seq2iseq(char *seq, int **p_iseq, int nres);



#endif /* ENERGY_H_ */
