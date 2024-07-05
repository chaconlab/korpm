/*
 * korp.h
 *
 *  Created on: Feb 15, 2018
 *      Author: mon
 */

#ifndef INCLUDE_KORPE_H_
#define INCLUDE_KORPE_H_

#include <stdio.h> // needed by some compilers
#include "libpdb/include/pdbIter.h"
//#include "libpdb/include/Atomo.h"

//#include "libpdb/include/Macromolecule.h"


#define DIST2SMALL 0.5           // Too-small inter-frame-center distance
#define CONTBLOCK (size_t)50000  // Number of contacts that are re-allocated at a time.
#define BONDING_THR 9            // Number of residues limit to consider bonding interactions (greater values will be non-bonding)
#define MAXSHELLS 100  // Maximum number of radial shells (for parser input)

typedef struct
{
   int ncells;      // Number of constant areas (cells)
   int nring;       // Output total number of rings (including polar caps, "nring" size)
   int *ncellsring; // Number of cells per ring (including polar caps, "nring" size)
   int *icell;      // Index of the first cell of each ring from north pole (including polar caps, "nring" size)
   float *theta;    // Ring boundaries array, i.e. Latitudes (Theta) [Radians] (including polar caps, "nring" size)
   float *dpsi;     // Psi increment per ring [Radians] (including polar caps, "nring" size)
   float dchi;      // Chi increment [Radians] (constant for all rings)
} mesh;

// Intraction Frame structure
typedef struct
{
	int i;       // Residue index
	char t;      // Frame type (depending on mapping)
	char ch;       // Chain id (required for bonding/non-bonding distinction)
	float p[3];  // Position
	float vx[3];
	float vy[3];
	float vz[3];
} frame;

// Contact structure
typedef struct
{
//	char pdbid[5]; // Pointer to PDB-IDs list (4-letter pdb-id + '\0' character)
	char ssA;      // Secondary structure of A
	char ssB;      // Secondary structure of B
	char seqA;     // Sequence of A (using our unique residue identifier, according to ResIni.h)
	char seqB;     // Sequence of B (using our unique residue identifier, according to ResIni.h)
	char froitA;    // FRame Of Interest Type (sequential index of the A interaction frame) for A residue
	char froitB;    // FRame Of Interest Type (sequential index of the B interaction frame) for B Residue
	char sd;       // Sequential Distance between A and B (bonded)
	float d;       // Distance A-B (module of Rab or Rba)
	float thetaA;  // Angle between Rab vector and vA z-axis or Order_AVE's Alpha angle
	float thetaB;  // Angle between Rba vector and vB z-axis or Order_AVE's Beta angle
	float psiA;    // Angle between the projection of Rab vector into the xy-plane and vA x-axis
	float psiB;    // Angle between the projection of Rba vector into the xy-plane and vB x-axis
	float chi;     // Dihedral angle between vA z-axis and vB z-axis or Order_AVE's Gamma angle
	bool intra;    // true --> Intra-chain contact, false --> Inter-chain contact (A belongs to the requested chain but B not)
    int  fai;      //residue index of A
    int  fbi;      //residue index of B
    char fach;     // Chain A
    char fbch;     // Chain B
	//	char sasA;     // Solvent Accessible Surface of A
//	char sasB;     // Solvent Accessible Surface of B
//	float vA[3];   // Cartesian projection of Rab into local frame A
//	float vB[3];   // Cartesian projection of Rba into local frame B
} contact;

// KORP energy map
typedef struct
{
	int model; // CG-model
	int dimensions; // Number of dimensions
	float cutoff; // Maximum distance-to-be-considered cutoff
	char frame_model; // Frame model index
	int nonbonding;
	int nonbonding2;
	float bonding_factor; // Bonding energy factor (typically 0.3)
	int ngauss;
	bool fullgauss;
	bool use_ji; // Using j-i information in 3D KORP
	bool each_bonding;
	bool use_bonding;

	char nft; // Number of frame types
	int nintres; // Number of interacting residues (to convert "aas" array into "iaa" array)
	int nintfra; // Number of interacting frames (some residues would be mapped together)
	char *aas; // Mapping identifiers to indices (WARNING, the sequential order must match "mapping" array)
	char *fras; // Frame IDs (just to dump some info)
	char **mapping; // Interaction frames mapping
	char *smapping; // smapping[10] = {0,0,0,0,0,0,0,0,0,0}; // Sequential distance Mapping. It returns the "s" index for the corresponding Bonding (or Non-Bonding) interaction

	int nr; // Number of radial bins (number of shells)
	int nchis; // Number of bins in Chi dimension
	float *br; // br[MAXSHELLS+1]; // Radial shells boundaries
	// int *ncells; // ncells[MAXSHELLS]; // Array with the number of cells per shell
	int ncells; // Number of cells in the "cutoff" shell
	int *scell; // Number of cells in each shell array
	float minr; // Minimum radial distance
	mesh **meshes; // Mesh data structures array
	char *iaa; // Given a residue unique identifier (an integer number 1 byte long) it returns the corresponding residue index in mapping
	int nslices; // Total number of interaction maps/data. Some day it will be easy adding i+1, i+2, etc. independent bonding maps.
	int nsmaps; // Number of independent maps (0= Non-bonding, 1= Bonding), when many bonding are considered it must be increased accordingly
	float *fmapping; // fmapping[10] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}; // Sequential distance bonding factor. It returns the "bonding factor" for the corresponding Bonding (or Non-Bonding) interaction

	// KORP potential map(s)
	float *******maps3D; // KORP xy-ij-3D energy Map (the full potential map)
	float *********maps; // Non-bonding KORP xy-ij-6D/4D energy Map (the full potential map)

} korp;


// Pablo's progress indicator bar
void indicator(char *text, int index, int size_loop);

float norm3D(float *a);
float *normalize3D(float *a);
float *sum3D(float *a, float *b, float *c);
float *diff3D(float *a, float *b, float *c);
void mult3D(float *a, float m); // Multiply a 3D vector "a" by one scalar "m"
float *cross3D(float *a, float *b, float *c);
// "Negative" cross-product (a = -a) (Useful for KORP's frames2ic function)
float *cross3DN(float *a, float *b, float *c);
float ang3D(float *a, float *b);
float *proj3D(float *a, float *b, float *c); // Project one 3D vector "a" into another "b"
float dihedral3Dunit(float *ua, float *ub, float *uc); // Computes the dihedral angle (in radians) formed by three consecutive 3D unit vectors.
// Computes the dihedral angle (in radians) formed by three consecutive 3D unit vectors.
// Using "Negative" cross-product (a = -a) (Useful for KORP's frames2ic function)
float dihedral3DunitN(float *ua, float *ub, float *uc);


// Read header data into map file (common to Binned maps and GMM)
void readMapHeader(FILE *f_map, int *dimensions, float *cutoff, char *frame_model, int *nonbonding, int *nonbonding2, float *bonding_factor, int *ngauss, bool *fullgauss, bool *use_ji, bool *each_bonding);

// Initialize frames mapping for a given frame model
// INPUT:
//  frame_model  Frame model
// OUTPUT:
//  p_mapping    Frame mapping
//  model        Coarse Graining model
//  dimensions   Number of dimensions (only in 4D frame models dimensions will be updated)
//  aas          Residue indices array (for residues)
//  fras         Residue indices array (for frames)
//  nft          Number of frame types
//  nintres      Number of interacting residues (to convert "aas" array into "iaa" array)
//  nintfra      Number of interacting frames (some residues would be mapped together)
void initFrames(char frame_model, char ***p_mapping, int *model, int *dimensions, char **aas, char **fras, char *nft, int *nintres, int *nintfra);


// Read meshes and related variables
//  nr     --> Number of radial bins
//  br     --> Boundaries for radial bins
//  scell  --> Number of cells per shell array
mesh **readMeshes(FILE *f_map, int *nr, float **p_br, int **scell );

// Read Order_AVE's radial boundaries and number or radial, alpha, beta, or gamma bins
//  nr  --> Number of radial bins
//  br  --> Boundaries for radial bins
//  nab --> Number of Alpha or Beta bins
//  ng  --> Number of Gamma bins
void readMeshOA(FILE *f_map, int *nr, float **p_br, int *nab, int *ng);

// Get raw meshes from binary FILE (automatic memory allocation)
// INPUT: f_map --> File handle (already open)
//        nr    --> Number of shells
// OUTPUT: Meshes array read
mesh **getMeshes(FILE *f_map, int nr);

// Reads a complete KORP energy map (valid for 3/4/6D)
//  INPUT:  file --> KORP map file name
//          bonding_factor --> (optional) if >= 0.0 --> Using input bonding factor, otherwise using map's bonding factor.
//  OUTPUT: KORP map pointer (memory automatically allocated within this function)
korp *readKORP(char *file, float bonding_factor = -1.0);


// 7D Integer/Float Map initialization: xy-frame-type + 20x20 (2D) + r(1D) + thetaA,psiA(2D)
void initMap3D(int *******p_map, int side);
void initMap3D(float *******p_map, int side);
// 5D Integer/Float Map initialization: 20x20 (2D) + r(1D) + thetaA,psiA(2D)
void initMap3D(int *****p_map, int side);
void initMap3D(float *****p_map, int side);
// 8D Integer/Float Map initialization: 20x20 (2D) + r(1D) + thetaA,psiA(2D) + thetaA,psiA(2D) + chi(1D)
void initMap6D(int *******p_map, int side);
void initMap6D(float *******p_map, int side);
// 10D Integer/Float Map initialization: xy-frame-type + 20x20 (2D) + r(1D) + thetaA,psiA(2D) + thetaA,psiA(2D) + chi(1D)
void initMap8D(int *********p_map, int side);
void initMap8D(float *********p_map, int side);

// Get raw Map from binary FILE handle (automatic memory allocation)
// INPUT: f_map --> File handle (already open)
//        meshes --> Meshes array
//        nr    --> Number of shells
// OUTPUT: p_smap --> Size of map [bytes]
//         return --> Map read
float **getMapICOSA(FILE *f_map, mesh **meshes, int nr, int *p_smap);
// Get raw Map from binary FILE handle (automatic memory allocation)
// INPUT: f_map --> File handle (already open)
//        meshes --> Meshes array
//        nr    --> Number of shells
// OUTPUT: p_smap --> Size of map [bytes]
//         return --> Map read
float ****getMap(FILE *f_map, mesh **meshes, int nr, int *p_smap);
// Get Order_AVE 4D raw Map from binary FILE handle (automatic memory allocation)
// INPUT:  f_map  --> File handle (already open)
//         nr     --> Number of radial bins
//         nab    --> Number of Alpha or Beta bins
//         ng     --> Number of Gamma bins
// OUTPUT: p_smap --> Size of map [bytes]
//         return --> Map read
float ****getMap(FILE *f_map, int nr, int nab, int ng, int *p_smap);


// Compute contacts for energy evaluation (automatic memory allocation)
//  pdb         --> Macromolecule
//  p_contacts  --> Pointer to the contacts array (if NULL, it allocates memory, otherwise memory pre-allocation expected)
//  cutoff      --> Maximum distance cutoff
//  iaa         --> Residues mapping into the mapping table
//  mapping     --> Mapping table with the atom indices to define the orthogonal framework and obtain the pairwise distances
//  nintres     --> Number of interacting residues
int contactPDB(Macromolecule *pdb, contact **p_contacts, float cutoff, char *iaa, char **mapping, int nintres);

// Compute frames for FAST energy evaluation (automatic memory allocation) [ONLY for single frame N-CA-C model (-f 10)]
// (Mapping not used for maximum speed. Only N-CA-C frame model allowed!)
//  coord       --> Coordinates (naive array of floats) of the 1st macromolecule
//  num_res     --> Number of residues (or frames) in the coordinates array
//  resnums     --> PDB's per residue numeration (required for proper bonding/non-bonding discrimination)
//  reschain    --> PDB's per residue chain-id (required for proper bonding/non-bonding discrimination)
frame *frameCoord(float *coord, int num_res, int *resnums, char *reschains);

// Compute contacts from frames for FAST energy evaluation (automatic memory allocation) [ONLY for single frame N-CA-C model (-f 10)]
// (Mapping not used for maximum speed. Only N-CA-C frame model allowed!)
//  p_contacts  --> Pointer to the contacts array (if NULL, it allocates memory, otherwise memory pre-allocation expected)
//  cutoff      --> Maximum distance cutoff
//  -- All vs. all (also valid for Loop vs. itself and Environment vs. itself modes) ---------------------------------------
//  frames      --> Frames (array of frames) of the 1st macromolecule
//  nres        --> Number of residues (or frames) of the 1st macromolecule
//  seq         --> Residues sequence specified by the unique identifier (as defined in ResIni.h or korpm.h)
//  -- Loop vs. Environment mode -------------------------------------------------------------------------------------------
//  anchorNt    --> (OPTIONAL) Absolute index of the Nt-anchor residue (the last non-mobile from Nt-end)
//  frames2     --> (OPTIONAL) Frames (array of frames) of the 2nd macromolecule (the Loop, only mobile residues)
//  nres2       --> (OPTIONAL) Number of mobile residues in the loop
//  seq2        --> (OPTIONAL) Residues sequence specified by the unique identifier (as defined in ResIni.h or korpm.h)
int contactCoord(contact **p_contacts, float cutoff, frame *frames, int nres, int *seq,
		                               int anchorNt=999999999, frame *frames2=NULL, int nres2=0, int *seq2=NULL);

// Compute contacts from frames for FAST energy evaluation (automatic memory allocation) [ONLY for single frame N-CA-C model (-f 10)]
// (Mapping not used for maximum speed. Only N-CA-C frame model allowed!)
//  p_contacts  --> Pointer to the contacts array (if NULL, it allocates memory, otherwise memory pre-allocation expected)
//  cutoff      --> Maximum distance cutoff
//  -- All vs. all (also valid for Loop vs. itself and Environment vs. itself modes) ---------------------------------------
//  frames      --> Frames (array of frames) of the 1st macromolecule
//  nres        --> Number of residues (or frames) of the 1st macromolecule
//  seq         --> Residues sequence specified by the unique identifier (as defined in ResIni.h or korpm.h)
//  -- --- -------------------------------------------------------------------------------------------
//  contact posM mUt
//  posM index of Mutation
int contactCoordM(contact **p_contactsM, float cutoff, frame *frames, int nres, int *seq, int posM, char chainM);

// the same with list
int contactCoordML(contact **p_contactsM, float cutoff, frame *frames, int nres, int *seq, int *posM, char *chainM, int nmut);

// Compute contacts ONLY for single frame N-CA-C model (-f 10) for FAST energy evaluation (automatic memory allocation)
// (Mapping not used for maximum speed. Only N-CA-C frame model allowed!)
//  coord       --> Coordinates (naive array of floats) of the 1st macromolecule
//  num_res     --> Number of residues of the 1st macromolecule
//  seq         --> Residues sequence specified by the unique identifier (as defined in ResIni.h)
//  p_contacts  --> Pointer to the contacts array (if NULL, it allocates memory, otherwise memory pre-allocation expected)
//  cutoff      --> Maximum distance cutoff
//  resnums     --> PDB's per residue numeration (required for proper bonding/non-bonding discrimination)
//  reschain    --> PDB's per residue chain-id (required for proper bonding/non-bonding discrimination)
//  -- Environment vs. Environment mode (Loop-less --------------------
//  anchorNt    --> (OPTIONAL) Absolute index of the Nt-anchor residue (the last non-mobile from Nt-end)
//  nresloop    --> (OPTIONAL) Number of mobile residues in the loop
//  -- Loop vs. Environment mode --------------------------------------
//  coord2      --> (OPTIONAL) Coordinates (naive array of floats) of the 2nd macromolecule (the Loop, only mobile residues)
//  seq2        --> (OPTIONAL) Residues sequence specified by the unique identifier (as defined in ResIni.h)
//  resnumsloop --> (OPTIONAL) Loop's PDB per residue numeration (required for proper bonding/non-bonding discrimination)
//  reschainloop--> (OPTIONAL) Loop's PDB per residue chain-id (required for proper bonding/non-bonding discrimination)
int contactCoord(float *coord, int num_res, int *seq, contact **p_contacts, float cutoff, int *resnums, char *reschains,
		int anchorNt=999999999, int nresloop=0, float *coord2=NULL, int *seq2=NULL, int *resnumsloop=NULL, char *reschainsloop=NULL);

// Compute Order_AVE (4D) contacts for energy evaluation (automatic memory allocation)
//  pdb         --> Macromolecule
//  p_contacts  --> Pointer to the contacts array (if NULL, it allocates memory, otherwise memory pre-allocation expected)
//  cutoff      --> Maximum distance cutoff
//  iaa         --> Residues mapping into the mapping table
//  mapping     --> Mapping table with the atom indices to define the orthogonal framework and obtain the pairwise distances
//  nintres     --> Number of interacting residues
int contactPDB_OA(Macromolecule *pdb, contact **p_contacts, float cutoff, char *iaa, char **mapping, int nintres);

// Computes the 5 internal coordinates (all 6 except distance) between A and B interacting frames
//  fa,fb --> interacting frame pointers
//  ta,tb --> Theta angles
//  pa,pb --> Psi angles
//  chi   --> Chi angle
void frames2ic(frame *fa, frame *fb, float *ta, float *tb, float *pa, float *pb, float *chi);

// Computes the 3 internal coordinates for the Order_AVE like potentials between (i,ir) and (j,jr) interacting frames
//  i,j   --> interacting atom coordinates (previously employed in distance evaluation)
//  ir,jr --> reference atoms to define the reference frame vector
//  a,b,g --> Alpha, Beta and Gamma angles
void frames2ic_OA(Tcoor i, Tcoor j, Tcoor ir, Tcoor jr, float *a, float *b, float *g);


// Return the 6D bin indices from one contact (INT version)
//   INPUT:   c      --> Contact
//            ms     --> Meshes
//            br     --> Boundaries of shell Radius
//   OUTPUT:  p_ir   --> Index of Radius
//            p_ita  --> Index of Theta-A
//            p_itb  --> Index of Theta-B
//            p_ipa  --> Index of Psi-A
//            p_ipb  --> Index of Psi-B
//            p_ic   --> Index of Chi
void contact2bins(contact *c, mesh **ms, float *br, int *p_ir, int *p_ita, int *p_itb, int *p_ipa, int *p_ipb, int *p_ic);

// Return the 4D bin indices from one Order_AVE contact (INT version)
//   INPUT:   c      --> Contact
//            br     --> Boundaries of shell Radius
//            nab    --> Number of Alpha or Beta bins
//            ng     --> Number of Gamma bins
//   OUTPUT:  p_ir   --> Index of Radius
//            p_ia   --> Index of Alpha
//            p_ib   --> Index of Beta
//            p_ig   --> Index of Gamma
void contact2bins_OA(contact *c, float *br, int nab, int ng, int *p_ir, int *p_ia, int *p_ib, int *p_ig);;

// KORP 3D energy function
float korp3D(contact *contacts, int icont, korp *map);
float korp3DM(contact *contacts, int icont, korp *map);
float korp3DM(contact *contacts, int icont, korp *map, int pos, char ch, int aa);
float korp3DM(contact *contacts, int icont, korp *map, int pos, char ch, int aa, float **faa);
float korp3DMW(contact *contacts, int icont, korp *map, int pos, char ch, int aa, double *W);


// KORP 4D energy function
float korp4D(contact *contacts, int icont, korp *map);
// KORP 6D energy function
double korp6D(contact *contacts, int icont, korp *map);
double korp6DW(contact *contacts, int icont, korp *map, double *W);
double korp6DM(contact *contacts, int icont, korp *map);
double korp6DM(contact *contacts, int icont, korp *map, int pos, char ch, int aa);
double korp6DM(contact *contacts, int icont, korp *map, int pos, char ch, int aa, float *faa);
double korp6DM_BIND(contact *contacts, int icont, korp *map, int pos, char ch, int aa, float *faa);
double korp6DMW(contact *contacts, int icont, korp *map, int pos, char ch, int aa, double *W);
double korp6DMW_MPM(contact *contacts, int icont, korp *map, int *pos, char *ch, int *aa, int nmut, double *W);
// MUT with 21 WEIGHTs KORP 6D energy function (21st is an additive term)
double korp6DMW21(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa, double *W);
double korp6DMWRSA(contact *contacts, int icont, korp *map, int pos, char ch, int aa, double *W, float RSA);

double korp6DMW_BIND(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa, double *W);

// KORP 6D energy function including weights (W)
double korp6DW2(contact *contacts, int icont, korp *map, double *W);

// Update Multiple Point Mutations (MPMs) in contacts list
int contact_MPM(contact *contacts, int icont, int *posM, char *chM,  int *Maa, int nmut);

// Get a per residue array of PDB's residue numbers (automatic memory allocation)
int *getResNums(Macromolecule *mol);
// Get a per residue array of PDB's chain-ids (automatic memory allocation)
char *getResChainIds(Macromolecule *mol);

// Prints into "f" file-handle all contact information in the "ncont" contacts
//  header -> Set true to output header as well
void print_contacts(FILE *f, contact *contacts, long ncont, bool header=false);


#endif /* INCLUDE_KORPE_H_ */
