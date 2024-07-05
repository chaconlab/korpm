/*
 * korpe.cpp
 *
 *  Created on: Feb 15, 2018
 *      Author: mon
 */

#include "korpe.h"
#include "korpm.h" // KORP mapping

// Pablo's progress indicator bar
void indicator(char *text, int index, int size_loop)
/*======================================================================
 *INDICATOR shows at the terminal screen the progression of a number
 *of cycles
 *----------------------------------------------------------------------
 *	INPUT:  - int SIZE_LOOP indicating the integer value that corresponds
 *            with a 100% bar
 * 			- int INDEX which is the current integer index and will be
 * 		      converted to a percentage from the 100% bar
 *          - char *pointer TEXT with text which will be in front of bar
 *----------------------------------------------------------------------
 *  OUTPUT: - explicitly VOID
 * 			- output to terminal screen upon execution
 *----------------------------------------------------------------------
 *  REMARK: - in the computational cycles nothing should be calculated
 *            since otherwise the bar representation is interrupted
 *          - uses C++ strings
 *          - SIZE_BAR could be made variable and hence is possible source
 *            of BUGs
 *=====================================================================*/
{
	int bar;
	int size_bar = 40;
	float ratio;
	ratio = (float) index/ (float) size_loop;
	bar = (int) floor( ratio * (float) size_bar );

	fprintf(stdout,"\r%s [",text);

	for(int i=0; i<bar; i++)
		fprintf(stdout,"=");
	fprintf(stdout,">");

	for(int i=bar+1; i<size_bar; i++)
		fprintf(stdout,"-");

	if(bar == size_bar)
		fprintf(stdout,"] %3.0f%% Finished!\n",ratio*100);
	else
		fprintf(stdout,"] %3.0f%%",ratio*100);
	fflush(stdout);
}


float norm3D(float *a)
{
	return sqrtf( powf(a[0],2) + powf(a[1],2) + powf(a[2],2) );
}

float *normalize3D(float *a)
{
	float d = sqrtf( powf(a[0],2) + powf(a[1],2) + powf(a[2],2) );
	a[0] /= d;
	a[1] /= d;
	a[2] /= d;
	return a;
}

float *sum3D(float *a, float *b, float *c)
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
	return c;
}

float *diff3D(float *a, float *b, float *c)
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
	return c;
}

// Multiply a 3D vector "a" by one scalar "m"
void mult3D(float *a, float m)
{
	a[0] *= m;
	a[1] *= m;
	a[2] *= m;
}

float *cross3D(float *a, float *b, float *c)
{
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
	return c;
}

// "Negative" cross-product (a = -a) (Useful for KORP's frames2ic function)
float *cross3DN(float *a, float *b, float *c)
{
	c[0] = -a[1]*b[2]+a[2]*b[1];
	c[1] = -a[2]*b[0]+a[0]*b[2];
	c[2] = -a[0]*b[1]+a[1]*b[0];
	return c;
}

float dot3D(float *a, float *b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Compute the angle between two 3D vectors (in radians)
float ang3D(float *a, float *b)
{
	//	return acosf ( dot3D(a,b)/(norm3D(a)*norm3D(b)) );
	float c = dot3D(a,b)/(norm3D(a)*norm3D(b));
	// Required to prevent Bug due to numerical instability
	if(c > 1.0)
		c = 1.0;
	if(c < -1.0)
		c = -1.0;
	return acosf( c );
}
// Project one 3D vector "a" into another "b"
float *proj3D(float *a, float *b, float *c)
{
	c[0] = b[0];
	c[1] = b[1];
	c[2] = b[2];
	normalize3D(c);
	mult3D(c,dot3D(a,c));
	return c;
}

// Computes the dihedral angle (in radians) formed by three consecutive 3D unit vectors.
float dihedral3Dunit(float *ua, float *ub, float *uc)
{
	//	Tcoor v1,v2,v3; // dummy vectors
	float v1[3],v2[3],v3[3]; // dummy vectors

	// CrossProduct(b_a, b_c, n1);
	cross3D(ua,ub,v1);
	// CrossProduct(b_c, c_d, n2);
	cross3D(ub,uc,v2);
	// CrossProduct(n1, b_c, m);
	cross3D(v1,ub,v3);

	// x = DotProduct(n1, n2);
	// y = DotProduct(m, n2);
	// atan2(y, x);
	return atan2( dot3D(v3,v2), dot3D(v1,v2) );
}

// Computes the dihedral angle (in radians) formed by three consecutive 3D unit vectors.
// Using "Negative" cross-product (a = -a) (Useful for KORP's frames2ic function)
float dihedral3DunitN(float *ua, float *ub, float *uc)
{
	//	Tcoor v1,v2,v3; // dummy vectors
	float v1[3],v2[3],v3[3]; // dummy vectors

	// CrossProduct(b_a, b_c, n1);
	cross3DN(ua,ub,v1); // "Negative" version because three consecutive 3D unit vectors
	// CrossProduct(b_c, c_d, n2);
	cross3D(ub,uc,v2);
	// CrossProduct(n1, b_c, m);
	cross3D(v1,ub,v3);

	// x = DotProduct(n1, n2);
	// y = DotProduct(m, n2);
	// atan2(y, x);
	return atan2( dot3D(v3,v2), dot3D(v1,v2) );
}

// Read header data into map file (common to Binned maps and GMM)
void readMapHeader(FILE *f_map, int *dimensions, float *cutoff, char *frame_model, int *nonbonding, int *nonbonding2, float *bonding_factor, int *ngauss, bool *fullgauss, bool *use_ji, bool *each_bonding)
{
	fread(dimensions, sizeof(int), 1, f_map); // Number of dimensions (3-, 4-, or 6-D)
	fread(cutoff, sizeof(float), 1, f_map); // Distance Cutoff used in map
	fread(frame_model, sizeof(char), 1, f_map); // Frame model used to obtain contacts
	fread(nonbonding, sizeof(int), 1, f_map); // Lower cutoff (nb < |i-j|)
	fread(nonbonding2, sizeof(int), 1, f_map); // Upper cutoff (|i-j| < nb2)
	fread(bonding_factor, sizeof(float), 1, f_map); // Bonding map factor

	// Mon: Add here --symNB and/or --useJI ????
	fread(ngauss, sizeof(int), 1, f_map); // Number of GMM gaussians (<=0 for Binned maps)

	char dummychar;
	fread(&dummychar, sizeof(char), 1, f_map); // Full (1) or Diagonal (0) gaussians
	if(dummychar > 0)
		*fullgauss = true; // Set either true or false for using Full (asymmetric) or Diagonal (spherical) Gaussian functions in the GMM
	else
		*fullgauss = false; // Set either true or false for using Full (asymmetric) or Diagonal (spherical) Gaussian functions in the GMM

	fread(&dummychar, sizeof(char), 1, f_map); // Use ji (1) or Not-use ji (0) interactions (only for 3D)
	if(dummychar > 0)
		*use_ji= true; // Set either true or false for using ji (1) or Not-using ji (0) interactions (only for 3D)
	else
		*use_ji = false; // Set either true or false for using ji (1) or Not-using ji (0) interactions (only for 3D)

	fread(&dummychar, sizeof(char), 1, f_map); // "Each bonding", to separate bonding contacts (i+1, i+2, etc...)
	if(dummychar > 0)
		*each_bonding = true; // Set either true or false for using ji (1) or Not-using ji (0) interactions (only for 3D)
	else
		*each_bonding = false; // Set either true or false for using ji (1) or Not-using ji (0) interactions (only for 3D)

	fread(&dummychar, sizeof(char), 1, f_map); // Some extra unused char... for future use

	int dummyint;
	fread(&dummyint, sizeof(int), 1, f_map); // Some extra unused integer... for future use
	fread(&dummyint, sizeof(int), 1, f_map); // Some extra unused integer... for future use
}

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
void initFrames(char frame_model, char ***p_mapping, int *model, int *dimensions, char **aas, char **fras, char *nft, int *nintres, int *nintfra)
{
	// Common parameters
	*aas = aasAA;
	*fras = frasAA;
	*nintres = 20;
	*nintfra = 20;

	// Specific parameters
	switch(frame_model)
	{
	case 10:  // 1-frame/aa using CA for distance-checking
		printf("Frame model %d => 1-frame/aa using CA for distances and N,C for frames\n",frame_model);
		*model = 0;
		*nft = 1; // number of frame types
		*p_mapping = mappingAA;
		break;
	case 11: // 1-frame/aa using CB for distance-checking
		printf("Frame model %d => 1-frame/aa using CB for distance and N,C for frames\n",frame_model);
		*model = 1;
		*nft = 1; // number of frame types
		*p_mapping = mappingAACB;
		break;
	case 12: // 1-frame/aa using O for distance-checking
		printf("Frame model %d => 1-frame/aa using O for distances and N,C for frames\n",frame_model);
		*model = 1;
		*nft = 1; // number of frame types
		*p_mapping = mappingAAO;
		break;
	case 13: // 1-frame/aa using C for distance-checking
		printf("Frame model %d => 1-frame/aa using C for distances and N,CA for frames\n",frame_model);
		*model = 0;
		*nft = 1; // number of frame types
		*p_mapping = mappingAAC;
		break;
	case 14: // 1-frame/aa using N for distance-checking
		printf("Frame model %d => 1-frame/aa using N for distances and CA,C for frames\n",frame_model);
		*model = 0;
		*nft = 1; // number of frame types
		*p_mapping = mappingAAN;
		break;
	case 20: // 2-frames/aa using CA (core-frame) and C (peptide-frame) for distance-checking
		printf("Frame model %d => 2-frames/aa: Core (CA for distances and N,C for frames) and Peptide (C for dist. and CA,O for frames)\n",frame_model);
		*model = 1; // MON: this should be "0", shouldn't it?
		*nft = 2; // number of frame types
		*p_mapping = mappingAP_CAC;
		break;
	case 21: // 2-frames/aa using CB (core-frame) and O (peptide-frame) for distance-checking
		printf("Frame model %d => 2-frames/aa: Core (CB for distances and N,C for frames) and Peptide (O for dist. and CA,C for frames)\n",frame_model);
		*model = 1;
		*nft = 2; // number of frame types
		*p_mapping = mappingAP_CBO;
		break;
	case 22: // 2-frames/aa using CA (core-frame) and O (peptide-frame) for distance-checking
		printf("Frame model %d => 2-frames/aa: Core (CA for distances and N,C for frames) and Peptide (O for dist. and CA,C for frames)\n",frame_model);
		*model = 1;
		*nft = 2; // number of frame types
		*p_mapping = mappingAP_CAO;
		break;
	case 26: // 2-frames/aa: O:CA/C (O for distances and CA,C for frames), N:CA/CB (N for dist. and CA,CB for fram.)
		printf("Frame model %d => 2-frames/aa: O:CA/C (O for distances and CA,C for frames), N:CA/CB (N for dist. and CA,CB for fram.)\n",frame_model);
		*model = 1;
		*nft = 2; // number of frame types
		*p_mapping = mappingAP_ON;
		break;
	case 36: // 3-frames/aa: O:CA/C (O for distances and CA,C for frames), N:CA/CB (N for dist. and CA,CB for fram.), CB:N/CA (CB for dist. and N/CA for fram.)
		printf("Frame model %d => 3-frames/aa: O:CA/C (O for distances and CA,C for frames), N:CA/CB (N for dist. and CA,CB for fram.), CB:N/CA (CB for dist. and N/CA for fram.)\n",frame_model);
		*model = 1;
		*nft = 3; // number of frame types
		*p_mapping = mappingAP_ONCB;
		break;
	case 55: // 5 Frames/aa:  1:CA/N^C, 2:C/CA^O, 4:O/CA^C 4:N/CA^C, and 5:CB/N^C (Distance/Frame)
		printf("Frame model %d => 5-frames/aa: 1:CA/N^C, 2:C/CA^O, 4:O/CA^C 4:N/CA^C, and 5:CB/N^C (Distance/Frame)\n",frame_model);
		*model = 1;
		*nft = 5; // number of frame types
		*p_mapping = mapping5F;
		break;
	case 51: // 5-frames/aa: 1:N/C-1, 2:CA/N, 3:C/CA, 4:O/C, and 5:CB/CA (Distance/Frame) (Only valid for Order_AVE's 4D method)"
		printf("Frame model %d => 5-frames/aa: 1:N/C-1, 2:CA/N, 3:C/CA, 4:O/C, and 5:CB/CA (Distance/Frame) (Only valid for 4D)\n",frame_model);
		*model = 1;
		*dimensions = 4; // This frame model is exclusive for 4D
		*nft = 5; // number of frame types
		*p_mapping = mappingOA;
		break;
	case 54: // 4-frames/aa: 1:N/C-1, 2:CA/N, 3:C/CA, and 4:O/C (Distance/Frame) (Only valid for Order_AVE's 4D method)"
		printf("Frame model %d => 4-frames/aa: 1:N/C-1, 2:CA/N, 3:C/CA, and 4:O/C (Distance/Frame) (Only valid for 4D)\n",frame_model);
		*model = 1;
		*dimensions = 4; // This frame model is exclusive for 4D
		*nft = 4; // number of frame types
		*p_mapping = mapping4Df54;
		break;
	case 52: // 2-frames/aa: 1:N/C-1 and 4:O/C (Distance/Frame) (Only valid for Order_AVE's 4D method)"
		printf("Frame model %d => 2-frames/aa: 1:N/C-1 and 4:O/C (Distance/Frame) (Only valid for 4D)\n",frame_model);
		*model = 1;
		*dimensions = 4; // This frame model is exclusive for 4D
		*nft = 2; // number of frame types
		*p_mapping = mapping4Df52;
		break;
	case 50: // 1-frames/aa: 1:CA/N (Distance/Frame) (Only valid for Order_AVE's 4D method)"
		printf("Frame model %d => 1-frame/aa: 1:CA/N (Distance/Frame) (Only valid for 4D)\n",frame_model);
		*model = 1;
		*dimensions = 4; // This frame model is exclusive for 4D
		*nft = 1; // number of frame types
		*p_mapping = mapping4Df50;
		break;
	case 53: // 3-frames/aa: 1:N/C-1, 2:O/C, and 3:CB/CA (Distance/Frame) (Only valid for Order_AVE's 4D method)"
		printf("Frame model %d => 3-frames/aa: 1:N/C-1, 2:O/C, and 3:CB/CA (Distance/Frame) (Only valid for 4D)\n",frame_model);
		*model = 1;
		*dimensions = 4; // This frame model is exclusive for 4D
		*nft = 3; // number of frame types
		*p_mapping = mapping4Df53;
		break;
	case 56: // 3-frames/aa: 1:N*/CA, 2:O/C, and 3:CB/CA (Distance/Frame) (Only valid for Order_AVE's 4D method)"
		printf("Frame model %d => 3-frames/aa: 1:N*/CA, 2:O/C, and 3:CB/CA (Distance/Frame) (Only valid for 4D)\n",frame_model);
		*model = 1;
		*dimensions = 4; // This frame model is exclusive for 4D
		*nft = 3; // number of frame types
		*p_mapping = mapping4Df56;
		break;

	default:
		printf("Please, select a valid frame model! (you've entered %d)\n Forcing exit!\n",frame_model);
		exit(2);
		break;
	}
}

// Read meshes and related variables
//  nr     --> Number of radial bins
//  br     --> Boundaries for radial bins
//  scell  --> Number of cells per shell array
mesh **readMeshes(FILE *f_map, int *nr, float **p_br, int **scell )
{
	// Getting "nr" from Map file
	fread(nr, sizeof(int), 1, f_map); // 3rd data is an integer with the Number of radial bins expected in file
	// Allocate memory for "br"
	*p_br = (float *) malloc( sizeof(float) * ((*nr) + 1) );
	// Getting "br" from Map file
	fread(*p_br, sizeof(float), *nr+1, f_map); // 4th, "nr+1" data is an array of floats with the Boundaries of radial bins
	// Computing Number of cells in each shell
	*scell = (int *) malloc( sizeof(int) * (*nr) );
	fread(*scell, sizeof(int), *nr, f_map); // 5th, "nr" data is an array of integers with the Number of cells within each shell

	// Get meshes
	mesh **meshes;
	meshes = getMeshes(f_map, *nr);
	return meshes;
}

// Read Order_AVE's radial boundaries and number or radial, alpha, beta, or gamma bins
//  nr  --> Number of radial bins
//  br  --> Boundaries for radial bins
//  nab --> Number of Alpha or Beta bins
//  ng  --> Number of Gamma bins
void readMeshOA(FILE *f_map, int *nr, float **p_br, int *nab, int *ng)
{
	// Getting "nr" from Map file
	fread(nr, sizeof(int), 1, f_map); // integer with the Number of radial bins expected in file
	// Allocate memory for "br"
	*p_br = (float *) malloc( sizeof(float) * ((*nr) + 1) );
	// Getting "br" from Map file
	fread(*p_br, sizeof(float), *nr+1, f_map); // array of floats with the Boundaries of radial bins
	// Getting "nab" from Map file
	fread(nab, sizeof(int), 1, f_map); // integer with the Number of Alpha or Beta bins expected in file
	// Getting "ng" from Map file
	fread(ng, sizeof(int), 1, f_map); // integer with the Number of Gamma bins expected in file
}

// Get raw meshes from binary FILE (automatic memory allocation)
// INPUT: f_map --> File handle (already open)
//        nr    --> Number of shells
// OUTPUT: Meshes array read
mesh **getMeshes(FILE *f_map, int nr)
{
	mesh **meshes;
	if( !(meshes = (mesh **) malloc( sizeof(mesh *) * nr ) ) )
	{
		fprintf(stderr,"getMeshes> Memory Error in %d shells allocation, forcing exit!\n",nr);
		exit(1);
	}

	for(int i=0; i<nr; i++) // Screen shells
	{
		if( !(meshes[i] = (mesh *) malloc( sizeof(mesh) ) ) )
		{
			fprintf(stderr,"getMeshes> Memory Error in mesh structure allocation, forcing exit!\n");
			exit(1);
		}

		fread(&meshes[i]->ncells, sizeof(int), 1, f_map);                   // 1 INT --> Number of cells (total)
		fread(&meshes[i]->nring, sizeof(int), 1, f_map);                    // 2 INT --> Number of rings
		meshes[i]->ncellsring = (int *) malloc( sizeof(int) * meshes[i]->nring); // Memory allocation
		fread(meshes[i]->ncellsring, sizeof(int), meshes[i]->nring, f_map); // 3 Array of INTs --> Number of cells per ring
		meshes[i]->icell = (int *) malloc( sizeof(int) * meshes[i]->nring); // Memory allocation
		fread(meshes[i]->icell, sizeof(int), meshes[i]->nring, f_map);      // 4 Array of INTs --> Index of the first cell per ring
		meshes[i]->theta = (float *) malloc( sizeof(float) * meshes[i]->nring); // Memory allocation
		fread(meshes[i]->theta, sizeof(float), meshes[i]->nring, f_map);    // 5 Array of FLOATs --> Theta angles per ring
		meshes[i]->dpsi = (float *) malloc( sizeof(float) * meshes[i]->nring); // Memory allocation
		fread(meshes[i]->dpsi, sizeof(float), meshes[i]->nring, f_map);     // 6 Array of FLOATs -->  Psi increment per ring
		fread(&meshes[i]->dchi, sizeof(float), 1, f_map);     				// 7 FLOAT -->  Chi increment per shell
	}

	return meshes;
}

// Reads a complete KORP energy map (valid for 3/4/6D)
//  INPUT:  file --> KORP map file name
//          bonding_factor --> (optional) if >= 0.0 --> Using input bonding factor, otherwise using map's bonding factor.
//  OUTPUT: KORP map pointer (memory automatically allocated within this function)
korp *readKORP(char *file, float bonding_factor)
{
	bool debug = false;
	char prog[] = "readKORP";
	FILE *f_map; // Binary potential map file (it includes Mesh info too)
	korp *map; // KORP's map data structure

	// Allocate one instance of the KORP's map data structure
	map = (korp *) malloc( sizeof(korp) * 1);

	// Opening Map file (if necessary)
	if(debug)
		printf("%s> Reading data from input potential map file: %s\n",prog,file);

	// Open for reading a binary map
	if( !(f_map = fopen(file,"rb") ) )
	{
		printf( "%s> Error, I can't read %s binary file! Forcing exit!\n", prog, file );
		exit(1);
	}

	// Read header data into map file (common to Binned maps and GMM)
	readMapHeader(f_map, &map->dimensions, &map->cutoff, &map->frame_model, &map->nonbonding, &map->nonbonding2, &map->bonding_factor, &map->ngauss, &map->fullgauss, &map->use_ji, &map->each_bonding);

	if(bonding_factor >= 0.0)
		map->bonding_factor = bonding_factor; // Overriding bonding_factor by user input

	// maxr = cutoff; // Maximum value of Rab

	if(map->nonbonding2 >= 0) // Trick to trigger Bonding/Non-bonding energy maps
		map->use_bonding = true;
	else
		map->use_bonding = false;

	//	if(!bf_specified) // use file's bonding_factor if parser's bonding_factor was not provided
	//		map->bonding_factor = dummybf; // map's bonding factor
	//	else
	//		map->bonding_factor = bonding_factor; // MON: do this?

	// Show read parameters
	fprintf(stdout,"%s> dimensions= %d  cutoff= %.3f  Frame model %d  nb= %d  nb2= %d  use_bonding= %d  bonding_factor= %.3f  ngauss= %d  fullgauss= %d  useJI= %d  each_bonding= %d found in %s map.\n",
			prog, map->dimensions, map->cutoff, map->frame_model, map->nonbonding, map->nonbonding2, map->use_bonding, map->bonding_factor, map->ngauss, map->fullgauss, map->use_ji, map->each_bonding, file);

	initFrames(map->frame_model, &map->mapping, &map->model, &map->dimensions, &map->aas, &map->fras, &map->nft, &map->nintres, &map->nintfra); // Frame model initialization

	// Get raw meshes from binary FILE (automatic memory allocation)
	if(map->dimensions != 4) // If non Order_AVE's (non 4D)
	{
		printf("%s> Reading raw meshes from %s\n",prog,file);
		map->meshes = readMeshes(f_map, &map->nr, &map->br, &map->scell );
		map->ncells = map->scell[map->nr-1]; // Number of cells in the "cutoff" shell
	}
	else
	{
		// Read Order_AVE's radial boundaries and number or radial, alpha, beta, or gamma bins
		printf("%s> Reading 4D mesh from %s\n",prog,file);
		readMeshOA(f_map, &map->nr, &map->br, &map->ncells, &map->nchis);
		printf("%s> nr= %d  nab= %d  ng= %d\n",prog,map->nr,map->ncells,map->nchis);
	}

	// the other parameters have been read above...
	map->minr = map->br[0]; // Set "minr"

	// Getting "nchis" form Map file
	printf("%s> Current potential map is %dD and frame model %d (parser input overridden by map)\n",prog,map->dimensions,map->frame_model);

	// Initialize 1-letter code into ASCII-code indices (fast trick ;-)
	map->iaa = (char *) malloc( sizeof(char) * 256);
	for(int i=0; i<256; i++)
		map->iaa[i] = -1; // Initialization (if iaa[i] == -1, then the residue is not mapped)
	for(int i=0; i<map->nintres; i++)
		map->iaa[(int)map->aas[i] ] = i;

	if(debug)
	{
		for(int i=0; i<256; i++)
			fprintf(stderr,"%d ",map->iaa[i]);
		fprintf(stderr,"\n");
	}

	// MON: nslices can be determined first and then allocate memory... consider review: Warning in 3D...further check...
	map->nslices = 10;
	map->fmapping = (float *) malloc( sizeof(float) * map->nslices ); // = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}; // Sequential distance bonding factor. It returns the "bonding factor" for the corresponding Bonding (or Non-Bonding) interaction
	for(int x=0; x < map->nslices; x++)
		map->fmapping[x] = 0.0; // Initialization... to avoid some valgrind Warnings...

	map->fmapping[0] = 1.0; // Non-Bonding always 1.0 (as reference point)
	if(map->each_bonding) // Consider each i+nb, i+nb+1, ..., i+nb+bnb, bonding interaction individually.
	{
		map->nslices = 1 + map->nonbonding2 - map->nonbonding; // |i-j| > nb to be considered, and |i-j| > bnb for non-bonding
		for(int x=1; x<map->nslices; x++)
			map->fmapping[x] = map->bonding_factor;
	}
	else // All i+nb, i+nb+1, ..., i+nb+bnb, bonding interactions are considered at the same time.
	{
		if(map->use_bonding)
		{
			map->nslices = 2; // 0= non-bonding and 1= bonding
			map->fmapping[1] = map->bonding_factor; // Bonding factor
		}
		else
			map->nslices = 1; // 0= non-bonding
	}

	map->nsmaps = map->nslices; // nsmaps is required for bonding/nonbonding handling...
	if(map->dimensions == 3 && map->use_ji) // in this case, the number of maps is doubled
		map->nslices *= 2;

	// Sequential Mapping (for bonding/non-bonding stuff)
	map->smapping = (char *) malloc( sizeof(char) * 10 );
	for(int x=0; x<map->nslices; x++)
		map->smapping[x] = 0; // Initialization... to avoid some valgrind Warnings...

	map->smapping[0] = 0; // Non-Bonding always go at "0" position
	fprintf(stdout,"%s> S-mapping: %d",prog,map->smapping[0]);
	int bi=1; // bonding index
	for(int i=1; i<10; i++)
	{
		if(i<=map->nonbonding)
			map->smapping[i] = -1; // Not-considered
		else
			if(i>map->nonbonding2 || !map->use_bonding)
				map->smapping[i] = 0; // Non-bonding
			else // Bonding
			{
				if(map->each_bonding)
				{
					map->smapping[i] = bi; // Bonding i+1, i+2, etc...
					bi++; // increase bonding index for "slice" mapping
				}
				else
					map->smapping[i] = 1; // Bonding (one map for all bonding)
			}
		fprintf(stdout," %d",map->smapping[i]);
	}
	fprintf(stdout,"\n");
	// This should change for considering i+1, i+2, etc... independent maps

	fprintf(stdout,"%s> F-mapping: ",prog);
	for(int i=0; i<10; i++)
		fprintf(stdout," %5.2f",map->fmapping[i]);
	fprintf(stdout,"\n");

	map->maps3D = (float *******) malloc( sizeof(float ******) * map->nslices);
	map->maps = (float *********) malloc( sizeof(float ********) * map->nslices);
	int smap = 0; // Size of Map [bytes]

	int nmaps; // Total number of individual potential maps
	int imaps = 0; // Current map index
	nmaps = map->nslices * map->nft * map->nft * map->nintfra * map->nintfra;

	printf("%s> Reading Bin-wise %dD maps:\n",prog,map->dimensions);
	for(int s=0; s<map->nslices; s++)
	{
		if(map->dimensions==3) // 3D
		{
			map->maps3D[s] = NULL; // required for automatic memory allocation
			initMap3D(&map->maps3D[s],map->nft); // allocate and initialize the first 2D (xy)

			smap=0;
			// printf("%s> Reading Bin-wise 3D maps:\n",prog);
			for(int x=0; x<map->nft; x++)
				for(int y=0; y<map->nft; y++) // AP and PA
				{
					initMap3D(&map->maps3D[s][x][y],map->nintfra); // allocate and initialize the first 2D (ij)

					for(int i=0; i<map->nintfra; i++)
					{
						// Skip missing intereacting frames, e.g. GLY does not have CB...
						if(x >= map->mapping[i][0]) // number of interacting frames for i-th residue
							continue;

						// for(int j=i; j<nintfra; j++)
						for(int j=0; j<map->nintfra; j++)
						{
							// Skip missing intereacting frames, e.g. GLY does not have CB...
							if(y >= map->mapping[j][0]) // number of interacting frames for j-th residue
								continue;

							// printf("%s> Reading Bin-wise 3D map: %1d %2d %2d %2d %2d\n",prog,s,x,y,i,j);

							map->maps3D[s][x][y][i][j] = getMapICOSA(f_map, map->meshes, map->nr, &smap); // Get raw Map from binary FILE (automatic memory allocation)
							char *msg   = (char *)"energy> Reading potential map";
							if(imaps % 10 == 0)
								indicator(msg,imaps,nmaps);
							imaps++; // Count maps read so far
						}
					}
				}
		}
		else // 4D or 6D
		{
			map->maps[s] = NULL;
			initMap8D(&map->maps[s],map->nft); // allocate and initialize the first 2D (xy)

			// printf("%s> Reading Bin-wise %dD maps:\n",prog,dimensions);
			smap=0;
			for(int x=0; x<map->nft; x++)
				for(int y=0; y<map->nft; y++) // AP and PA
				{
					initMap6D(&map->maps[s][x][y],map->nintfra); // allocate and initialize the first 2D (ij)

					// char iname[5],jname[5];

					for(int i=0; i<map->nintfra; i++)
					{
						// Skip missing intereacting frames, e.g. GLY does not have CB...
						if(x >= map->mapping[i][0]) // number of interacting frames for i-th residue
							continue;

						// resname_from_resnum(aas[i],iname);

						for(int j=0; j<map->nintfra; j++)
						{
							// Skip missing intereacting frames, e.g. GLY does not have CB...
							if(y >= map->mapping[j][0]) // number of interacting frames for j-th residue
								continue;

							// resname_from_resnum(aas[j],jname);
							// printf("%s> Reading Bin-wise %dD map (%s-%s): %1d %2d %2d %2d %2d\n",prog,dimensions,iname,jname,s,x,y,i,j);

							// fprintf(stdout,"%s> Reading 6D-map %c-%c (i-j and j-i)\n",prog,aas[i],aas[j]);
							if(map->dimensions == 4) // 4D
								map->maps[s][x][y][i][j] = getMap(f_map, map->nr, map->ncells, map->nchis, &smap); // Get raw 6D Map from binary FILE (automatic memory allocation)
							else // 6D
								map->maps[s][x][y][i][j] = getMap(f_map, map->meshes, map->nr, &smap); // Get raw 6D Map from binary FILE (automatic memory allocation)
							char *msg   = (char *)"energy> Reading potential map";
							if(imaps % 10 == 0)
								indicator(msg,imaps,nmaps);
							imaps++; // Count maps read so far
						}
					}
				}
		}
		// printf("%s> All %dD-maps read! (%d bytes, %f MB)\n",prog,dimensions,smap,smap/pow(2,20));
	}
	char *msg   = (char *)"energy> Reading potential map";
	indicator(msg,nmaps,nmaps);
	printf("%s> Map reading stuff finished!\n",prog);

	return map;
}


// 7D Integer/Float Map initialization: xy-frame-type + 20x20 (2D) + r(1D) + thetaA,psiA(2D)
void initMap3D(int *******p_map, int side)
{
	int ******aamaps; // 3D Maps for the 210 non-redundant interactions (each aamaps[i][j] will contain a complete map)
	aamaps = (int ******) malloc( sizeof(int *****) * side );
	for(int i=0; i<side; i++)
	{
		aamaps[i] = (int *****) malloc( sizeof(int ****) * side );
		for(int j=0; j<side; j++)
			aamaps[i][j] = NULL; // NULL required for automatic memory allocation later...
	}
	*p_map = aamaps; // output
}
void initMap3D(float  *******p_map, int side)
{
	float ******aamaps; // 3D Maps for the 210 non-redundant interactions (each aamaps[i][j] will contain a complete map)
	aamaps = (float ******) malloc( sizeof(float *****) * side );
	for(int i=0; i<side; i++)
	{
		aamaps[i] = (float *****) malloc( sizeof(float ****) * side );
		for(int j=0; j<side; j++)
			aamaps[i][j] = NULL;
	}
	*p_map = aamaps; // output
}

// 6D Integer/Float Map initialization: 20x20 (2D) + r(1D) + thetaA,psiA(2D)
void initMap3D(int *****p_map, int side)
{
	int ****aamaps; // 3D Maps for the 210 non-redundant interactions (each aamaps[i][j] will contain a complete map)
	aamaps = (int ****) malloc( sizeof(int ***) * side );
	for(int i=0; i<side; i++)
	{
		aamaps[i] = (int ***) malloc( sizeof(int **) * side );
		for(int j=0; j<side; j++)
			aamaps[i][j] = NULL; // NULL required for automatic memory allocation later...
	}
	*p_map = aamaps; // output
}
void initMap3D(float  *****p_map, int side)
{
	float ****aamaps; // 3D Maps for the 210 non-redundant interactions (each aamaps[i][j] will contain a complete map)
	aamaps = (float ****) malloc( sizeof(float ***) * side );
	for(int i=0; i<side; i++)
	{
		aamaps[i] = (float ***) malloc( sizeof(float **) * side );
		for(int j=0; j<side; j++)
			aamaps[i][j] = NULL;
	}
	*p_map = aamaps; // output
}

// 8D Integer Map initialization: 20x20 (2D) + r(1D) + thetaA,psiA(2D) + thetaA,psiA(2D) + chi(1D)
void initMap6D(int *******p_map, int side)
{
	int ******aamaps; // 6D Maps for the 210 non-redundant interactions (each aamaps[i][j] will contain a complete map)
	aamaps = (int ******) malloc( sizeof(int *****) * side );
	for(int i=0; i<side; i++)
	{
		aamaps[i] = (int *****) malloc( sizeof(int ****) * side );
		for(int j=0; j<side; j++)
			aamaps[i][j] = NULL; // NULL required for automatic memory allocation later...
	}
	*p_map = aamaps; // output
}

// 8D Float Map initialization: 20x20 (2D) + r(1D) + thetaA,psiA(2D) + thetaA,psiA(2D) + chi(1D)
void initMap6D(float  *******p_map, int side)
{
	float ******aamaps; // 6D Maps for the 210 non-redundant interactions (each aamaps[i][j] will contain a complete map)
	aamaps = (float ******) malloc( sizeof(float *****) * side );
	for(int i=0; i<side; i++)
	{
		aamaps[i] = (float *****) malloc( sizeof(float ****) * side );
		for(int j=0; j<side; j++)
			aamaps[i][j] = NULL;
	}
	*p_map = aamaps; // output
}


// 10D Integer Map initialization: xy-frame-type + 20x20 (2D) + r(1D) + thetaA,psiA(2D) + thetaA,psiA(2D) + chi(1D)
void initMap8D(int *********p_map, int side)
{
	int ********aamaps; // 6D Maps for the 210 non-redundant interactions (each aamaps[i][j] will contain a complete map)
	aamaps = (int ********) malloc( sizeof(int *******) * side );
	for(int i=0; i<side; i++)
	{
		aamaps[i] = (int *******) malloc( sizeof(int ******) * side );
		for(int j=0; j<side; j++)
			aamaps[i][j] = NULL; // NULL required for automatic memory allocation later...
	}
	*p_map = aamaps; // output
}

// 10D Float Map initialization: xy-frame-type + 20x20 (2D) + r(1D) + thetaA,psiA(2D) + thetaA,psiA(2D) + chi(1D)
void initMap8D(float  *********p_map, int side)
{
	float ********aamaps; // 6D Maps for the 210 non-redundant interactions (each aamaps[i][j] will contain a complete map)
	aamaps = (float ********) malloc( sizeof(float *******) * side );
	for(int i=0; i<side; i++)
	{
		aamaps[i] = (float *******) malloc( sizeof(float ******) * side );
		for(int j=0; j<side; j++)
			aamaps[i][j] = NULL;
	}
	*p_map = aamaps; // output
}



// Get raw Map from binary FILE handle (automatic memory allocation)
// INPUT: f_map --> File handle (already open)
//        meshes --> Meshes array
//        nr    --> Number of shells
// OUTPUT: p_smap --> Size of map [bytes]
//         return --> Map read
float **getMapICOSA(FILE *f_map, mesh **meshes, int nr, int *p_smap)
{
	float **map;
	int smap; // size of map [bytes]

	map = (float **) malloc( sizeof(float *) * nr); // allocate radial dimension (R)
	smap = sizeof(float *) * nr;

	for(int r=0; r<nr; r++) // screen shells (radial bins)
	{
		map[r] = (float *) malloc( sizeof(float) * meshes[r]->ncells); // allocate shell for A residue (ThetaA and PsiA dimensions)
		smap += sizeof(float) * meshes[r]->ncells;
		fread(map[r], sizeof(float), meshes[r]->ncells, f_map); // Read all orientation bins for A residue
	}

	*p_smap = smap;
	return map;
}

// Get raw Map from binary FILE handle (automatic memory allocation)
// INPUT: f_map --> File handle (already open)
//        meshes --> Meshes array
//        nr    --> Number of shells
// OUTPUT: p_smap --> Size of map [bytes]
//         return --> Map read
float ****getMap(FILE *f_map, mesh **meshes, int nr, int *p_smap)
{
	float ****map;
	int smap=0; // size of map [bytes]
	int nc;

	map = (float ****) malloc( sizeof(float ***) * nr); // allocate radial dimension (R)
	smap += sizeof(float ***) * nr;

	for(int r=0; r<nr; r++) // Screen meshes (shells)
	{
		nc = (int) roundf( 2*M_PI / meshes[r]->dchi );

		map[r] = (float ***) malloc( sizeof(float **) * meshes[r]->ncells ); // Allocate dimensions (ThetaA,PsiA) for A residue
		smap += sizeof(float **) * meshes[r]->ncells;

		for(int a=0; a<meshes[r]->ncells; a++)
		{
			map[r][a] = (float **) malloc( sizeof(float *) * meshes[r]->ncells ); // Allocate dimensions (ThetaB,PsiB) for B residue
			smap += sizeof(float *) * meshes[r]->ncells;
			for(int b=0; b<meshes[r]->ncells; b++)
			{
				map[r][a][b] = (float *) malloc( sizeof(float) * nc ); // Allocate cells for Chis
				smap += sizeof(float) * nc;
				for(int c=0; c<nc; c++)
				{
					fread(&map[r][a][b][c], sizeof(float), 1, f_map); // Read all orientation bins for A residue
					// fprintf(stderr,"r= %d  a= %d  b= %d  c= %d\n",r,a,b,c);
				}
			}
		}
	}

	*p_smap += smap;
	return map;
}

// Get Order_AVE 4D raw Map from binary FILE handle (automatic memory allocation)
// INPUT:  f_map  --> File handle (already open)
//         nr     --> Number of radial bins
//         nab    --> Number of Alpha or Beta bins
//         ng     --> Number of Gamma bins
// OUTPUT: p_smap --> Size of map [bytes]
//         return --> Map read
float ****getMap(FILE *f_map, int nr, int nab, int ng, int *p_smap)
{
	float ****map;
	int smap = 0; // size of map [bytes]

	map = (float ****) malloc( sizeof(float ***) * nr); // allocate radial dimension (R)
	smap += sizeof(float ***) * nr;

	for(int r=0; r<nr; r++) // Screen meshes (shells)
	{
		map[r] = (float ***) malloc( sizeof(float **) * nab ); // Allocate Alpha bins
		smap += sizeof(float **) * nab;

		for(int a=0; a<nab; a++)
		{
			map[r][a] = (float **) malloc( sizeof(float *) * nab ); // Allocate Beta bins
			smap += sizeof(float *) * nab;
			for(int b=0; b<nab; b++)
			{
				map[r][a][b] = (float *) malloc( sizeof(float) * ng ); // Allocate Gamma bins
				smap += sizeof(float) * ng;
				for(int g=0; g<ng; g++)
				{
					fread(&map[r][a][b][g], sizeof(float), 1, f_map); // Read all orientation bins
					// fprintf(stderr,"r= %d  a= %d  b= %d  g= %d\n",r,a,b,g);
				}
			}
		}
	}

	*p_smap += smap;
	return map;
}


// Compute contacts for energy evaluation (automatic memory allocation)
//  pdb         --> Macromolecule
//  p_contacts  --> Pointer to the contacts array (if NULL, it allocates memory, otherwise memory pre-allocation expected)
//  cutoff      --> Maximum distance cutoff
//  iaa         --> Residues mapping into the mapping table
//  mapping     --> Mapping table with the atom indices to define the orthogonal framework and obtain the pairwise distances
//  nintres     --> Number of interacting residues
int contactPDB(Macromolecule *pdb, contact **p_contacts, float cutoff, char *iaa, char **mapping, int nintres)
{
	//	fprintf(stderr,"Hi, I'm contactPDB\n");

	char prog[] = "contactPDB";
	//int num_atoms = pdb->get_num_atoms();
	int num_res = pdb->get_num_fragments();

	//	fprintf(stderr,"nintres= %d  num_atoms= %d  num_res= %d\n",nintres,num_atoms,num_res);
	//	pdb->writePDB("contactPDB.pdb");
	//	fprintf(stderr,"Written contactPDB.pdb\n");
	//	exit(0);

	// Maximum number of interaction frames per residue (to allocate memory for the maximum number of FRames of Interest)
	int maxint = 0;
	for(int i=0; i<nintres; i++)
	{
		if(mapping[i][0] > maxint)
			maxint = mapping[i][0];
		//		fprintf(stderr,"mapping[%d][0] = %d\n",i,mapping[i][0]);
	}
	//	for(int i=0; i<30;i++)
	//		fprintf(stderr,"iaa[%d] = %d\n",i,iaa[i]);
	//	fprintf(stderr,"nintres= %d  maxint= %d  num_atoms= %d  num_res= %d\n",nintres,maxint,num_atoms,num_res);

	// Allocate memory for residue frames
	frame *frames;
	frames = (frame *) malloc(sizeof(frame) * num_res * maxint); // allocate the maximum possible

	// The residues sequence must be specified by unique identifier (as defined in ResIni.h)
	char *seq;
	seq = (char *) malloc( sizeof(char) * num_res * maxint);

	// Some variables...
	pdbIter *iter_chain,*iter_res,*iter_atom;
	Residue *res;
	//	Segment *seg;
	Chain *ch;
	Tcoor N,CA,C,r12,r13;

	//	fprintf(stderr,"\n");

	// COMPUTING RESIDUE FRAMES
	//	iter_seg = new pdbIter( pdb ); // iter to screen segments
	iter_chain = new pdbIter( pdb ); // iter to screen chains
	//	int ires = 0; // Residue index
	int ifram = 0; // Frame index
	char *imap; // current interactions Mapping
	char ichain; // current chain
	frame *fram; // current frame pointer
	//	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // Iter segments to compute sequential distances fine
	for ( iter_chain->pos_chain = 0; !iter_chain->gend_chain(); iter_chain->next_chain() ) // Iter chains to compute sequential distances fine
	{
		//		seg = (Segment *) iter_seg->get_segment();
		//		iter_res = new pdbIter( seg );
		ch = (Chain *) iter_chain->get_chain();
		iter_res = new pdbIter( ch );

		// Getting the chain id of current segment
		//		iter_seg->pos_segment = 0;
		//		fprintf(stderr,"Before\n");

		//		ch = (Chain *) seg->getFather(); // ->get_chain();
		ichain = (ch->getName())[0]; // Current chain id (1-letter)
		//		fprintf(stderr,"ichain= %c\n",ichain);

		for ( iter_res->pos_fragment = 0; !iter_res->gend_fragment(); iter_res->next_fragment() ) // iter frags
		{
			//			if(ires < 83 || ires > 90) // LOOP IS: 83, 84, 85, 86, 87, 88, 89, 90
			//			if(ires >= 83 && ires <= 90) // LOOP IS: 83, 84, 85, 86, 87, 88, 89, 90
			// LEFT SIDE + LOOP
			//			if(ires <= 90) // LOOP IS: 83, 84, 85, 86, 87, 88, 89, 90
			//			{

			res = ( Residue * ) iter_res->get_fragment();

			//			if(res->getIdNumber() >= 42 && res->getIdNumber() <= 53) // LOOP IS:
			//			if(res->getIdNumber() < 42 || res->getIdNumber() > 53) // LOOP IS:
			//			{

			iter_atom = new pdbIter( res ); // screens residue atoms

			// fprintf(stderr,"res->id= %s\n",res->getName());
			// fprintf(stderr,"res->IdNumber= %d\n",res->getIdNumber());
			// fprintf(stderr,"res->fragid= %d\n",res->fragid);

			imap = mapping[(int)iaa[ (int)res->fragid ] ];

			// fprintf(stderr,"res->id= %s res->fragid= %d  imap[0]= %d\n",res->getName(), (int) res->fragid, imap[0]);

			for(int t=0; t<imap[0]; t++)
			{
				fram = &frames[ifram];

				seq[ifram] = res->fragid; // Get residue identifier
				// fprintf(stderr,"%d ",seq[ifram]);

				// MON: use PDB residue number instead of sequential "ires" id
				// fram->i = ires; // store residue index (required later)
				fram->i = res->getIdNumber(); // Get PDB numeration (required for correct bonding/non-bonding discrimination)
				fram->t = t; // store frame type (required later)
				fram->ch = ichain; // store chain id (required later)

				// fprintf(stderr,"fram->i= %d  fram->t= %d\n",fram->i,fram->t);

				iter_atom->pos_atom = imap[1+4*t+1]; // N's position
				(iter_atom->get_atom())->getPosition(N);
				iter_atom->pos_atom = imap[1+4*t+2]; // CA's position
				(iter_atom->get_atom())->getPosition(CA);
				iter_atom->pos_atom = imap[1+4*t+3]; // C's position
				(iter_atom->get_atom())->getPosition(C);

				// Compute residue frames (Eq.1 from GOAP's paper Zhou and Skolnik 2011. A1=N, A=CA, A2=C)
				diff3D(N,CA,r12);
				diff3D(C,CA,r13);
				sum3D(r12, r13, fram->vz);
				normalize3D(fram->vz); // Vz computed
				cross3D(fram->vz,r13,fram->vy);
				normalize3D(fram->vy); // Vy computed
				cross3D(fram->vy,fram->vz,fram->vx); // Vx computed

				// Store interaction distances (CA-positions) into frames array
				fram->p[0] = CA[0];
				fram->p[1] = CA[1];
				fram->p[2] = CA[2];

				// fprintf(stderr,"i=%4d   fram->p= %f %f %f\n",ifram,fram->p[0],fram->p[1],fram->p[2]);

				ifram++; // update frame index
			}

			delete iter_atom;

			//			}
			//			ires++; // update current residue index
			//			}
		}
		delete iter_res;
	}
	delete iter_chain;
	// frames = (frame *) realloc(frames, sizeof(frame) * ifram ); // Reallocating to exactly store "ifram" interaction frames.
	// seq = (char *) realloc(seq, sizeof(char) * ifram ); // Reallocating to exactly store "ifram" sequence chars.
	// exit(0);
	//	fprintf(stderr,"\n");

	// Allocating initial memory for contacts
	contact *contacts; // Contacts array
	if(*p_contacts == NULL) // Contact list is empty (first PDB typically)
		*p_contacts = contacts = (contact *) malloc( sizeof(contact) * CONTBLOCK); // Allocating first memory block to store contacts
	else
		contacts = *p_contacts; // already "mallocated"

	// COMPUTING CONTACT DATA
	// Screen all inter-residue contacts
	float dab,ta,tb,pa,pb,chi; // inter-frame internal coordinates
	int bond;
	int icont=0; // contact index for current protein
	frame *fa,*fb;
	float cutoff2;
	cutoff2 = cutoff * cutoff; // cutoff squared for efficient distance evaluation

	for(int a = 0; a < ifram; a++) // Screen all interaction frames
		for(int b = a+1; b < ifram; b++) // Non-redundant screening
		{
			// dab = dist3D(frames[b].p,frames[a].p);
			// if(dab < cutoff)
			fa = &frames[a];
			fb = &frames[b];

			// This one seems the fastest alternative... but all are very close...
			if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
					if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
					{
						dab = sqrtf(dab);

						if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
						{
							fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
							// exit(2);
						}

						frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);

						// Bonding/Non-bonding discrimination
						if( fa->ch != fb->ch ) // If residue indices belong to different chains
							bond = 0; // always non-bonding interaction between frames from different chains
						else
						{
							// bond = (int)fabsf(fa->i - fb->i); // Mon: abs is not necessary, isn't it?
							// if( bond > BONDING_THR)
							//	bond = 0; // non-bonding

							if( fb->i >= fa->i ) // Sanity check...
							{
								bond = fb->i - fa->i; // (bonding or non-bonding)
								if( bond > BONDING_THR)
									bond = 0; // non-bonding
							}
							else
							{
								fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. (fb->ch= %c  fa->ch= %c) Forcing exit!\n",fb->i,fa->i,fb->ch,fa->ch);
								exit(2);
							}
						}

						// When contact memory is over this allocates a new contacts block
						if(icont % CONTBLOCK == 0)
						{
							// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
							fflush(stdout);
							if( !(*p_contacts = contacts = (contact *) realloc(contacts, sizeof(contact) * (icont+CONTBLOCK))) )  // Reallocating to store new contacts
							{
								fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icont/powf(2,20));
								exit(1);
							}
						}

						// Store in RAM the computed data
						contacts[icont].ssA = '-'; // Not computed
						contacts[icont].ssB = '-';
						contacts[icont].seqA = seq[a];
						contacts[icont].seqB = seq[b];
						// fprintf(stderr,"seq[a]= %d  seq[b]= %d  seqA= %d  seqB= %d\n",seq[a],seq[b],contacts[icont].seqA,contacts[icont].seqB);
						contacts[icont].froitA = fa->t;
						contacts[icont].froitB = fb->t;
						contacts[icont].sd = bond;
						contacts[icont].thetaA = ta;
						contacts[icont].thetaB = tb;
						contacts[icont].psiA = pa;
						contacts[icont].psiB = pb;
						contacts[icont].chi = chi;
						contacts[icont].d = dab;

						icont++; // Update contact index
					}
		}
	//	printf("%s> %d contacts found (%.1f/residue at %.0f A cutoff), %.1f%% from %d available.\n",
	//			prog,icont,(float)icont/num_res,cutoff,100*(float)icont/(num_res*(num_res-1)/2),num_res*(num_res-1)/2);
	//	printf("%s> Memory used for contacts: %ld bytes (%.1f MB)\n\n", prog, sizeof(contact)*icont, (float)sizeof(contact)*icont/powf(2,20));

	// Free memory
	free(seq);
	free(frames);

	return icont; // Successful exit
}

// Compute frames for FAST energy evaluation (automatic memory allocation) [ONLY for single frame N-CA-C model (-f 10)]
// (Mapping not used for maximum speed. Only N-CA-C frame model allowed!)
//  coord       --> Coordinates (naive array of floats) of the 1st macromolecule
//  num_res     --> Number of residues (or frames) in the coordinates array
//  resnums     --> PDB's per residue numeration (required for proper bonding/non-bonding discrimination)
//  reschain    --> PDB's per residue chain-id (required for proper bonding/non-bonding discrimination)
frame *frameCoord(float *coord, int num_res, int *resnums, char *reschains)
{
	// Allocate memory for residue frames
	frame *frames;
	frames = (frame *) malloc(sizeof(frame) * num_res); // allocate the maximum possible

	// The residues sequence must be specified by unique identifier (as defined in ResIni.h)
	//	char *seq;
	//	seq = (char *) malloc( sizeof(char) * num_res);

	// Some variables...
	Tcoor r12,r13;

	//	fprintf(stderr,"seq=\n");
	//	for(int i=0; i<num_res; i++)
	//		fprintf(stderr,"%d ",seq[i]);
	//	fprintf(stderr,"\n");
	//	fprintf(stderr,"anchorNt %d %d", anchorNt, nresloop);

	// COMPUTING MACROMOLECULE FRAMES (FULL OR ENVIRONMENT)
	frame *fram; // current frame pointer
	for(int i=0; i<num_res; i++) // Screen coordinates of 1st macromolecule
	{
		fram = &frames[i];
		//		if(i>anchorNt)
		//			fram->i = i + nresloop; // store residue index (required later)
		//		else
		//			fram->i = i; // store residue index (required later)
		fram->i = resnums[i]; // using PDB's residue numeration
		fram->ch = reschains[i]; // using PDB's chain-ids
		fram->t = 0; // store frame type (required later)
		// fprintf(stderr,"fram->i= %d  fram->t= %d\n",fram->i,fram->t);

		// Compute residue frames (Eq.1 from GOAP's paper Zhou and Skolnik 2011. A1=N, A=CA, A2=C)
		diff3D(coord+9*i,coord+9*i+3,r12);   // diff3D(N,CA,r12);
		diff3D(coord+9*i+6,coord+9*i+3,r13); // diff3D(C,CA,r13);
		sum3D(r12, r13, fram->vz);
		normalize3D(fram->vz); // Vz computed
		cross3D(fram->vz,r13,fram->vy);
		normalize3D(fram->vy); // Vy computed
		cross3D(fram->vy,fram->vz,fram->vx); // Vx computed

		// Store interaction distances (CA-positions) into frames array
		fram->p[0] = coord[9*i+3];
		fram->p[1] = coord[9*i+4];
		fram->p[2] = coord[9*i+5];
	}
	return frames;
}

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
		int anchorNt, frame *frames2, int nres2, int *seq2)
{
	char prog[] = "contactCoord";

	// Allocating initial memory for contacts
	contact *contacts; // Contacts array
	if(*p_contacts == NULL) // Contact list is empty (first PDB typically)
		*p_contacts = contacts = (contact *) malloc( sizeof(contact) * CONTBLOCK); // Allocating first memory block to store contacts
	else
		contacts = *p_contacts; // already "mallocated"

	// COMPUTING CONTACT DATA
	// Screen all inter-residue contacts
	float dab,ta,tb,pa,pb,chi; // inter-frame internal coordinates
	int bond;
	int icont=0; // contact index for current protein
	frame *fa,*fb;
	float cutoff2;
	cutoff2 = cutoff * cutoff; // cutoff squared for efficient distance evaluation

	if(frames2 == NULL) // Standard mode (all vs. all mode)
	{
		for(int a = 0; a < nres; a++) // Screen all interaction frames
			for(int b = a+1; b < nres; b++) // Non-redundant screening
			{
				fa = &frames[a];
				fb = &frames[b];

				// This one seems the fastest alternative... but all are very close...
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
					if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
						if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
						{
							dab = sqrtf(dab);

							if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
							{
								fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
								// exit(2);
							}

							frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);


							//							bond = (int)fabsf(fa->i - fb->i); // Mon: abs is not necessary, isn't it?
							//							if( bond > BONDING_THR)
							//								bond = 0;

							// Bonding/Non-bonding discrimination
							if( fa->ch != fb->ch ) // If residue indices belong to different chains
								bond = 0; // always non-bonding interaction between frames from different chains
							else
							{
								if( fb->i >= fa->i ) // Sanity check...
								{
									bond = fb->i - fa->i; // (bonding or non-bonding)
									if( bond > BONDING_THR)
										bond = 0; // non-bonding
								}
								else
								{
									//									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. Forcing exit!\n",fb->i,fa->i);
									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. (fb->ch= %c  fa->ch= %c) Forcing exit!\n",fb->i,fa->i,fb->ch,fa->ch);
									exit(2);
								}
							}

							// When contact memory is over this allocates a new contacts block
							if(icont % CONTBLOCK == 0)
							{
								// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
								fflush(stdout);
								if( !(*p_contacts = contacts = (contact *) realloc(contacts, sizeof(contact) * (icont+CONTBLOCK))) )  // Reallocating to store new contacts
								{
									fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icont/powf(2,20));
									exit(1);
								}
							}

							// Store in RAM the computed data
							contacts[icont].ssA = 'N'; // Not computed
							contacts[icont].ssB = 'N';
							contacts[icont].seqA = seq[a];
							contacts[icont].seqB = seq[b];



							contacts[icont].froitA = fa->t;
							contacts[icont].froitB = fb->t;
							//							contacts[icont].froitA = (unsigned int) fa->i;
							//							contacts[icont].froitB = (unsigned int) fb->i;
							contacts[icont].sd = bond;
							contacts[icont].thetaA = ta;
							contacts[icont].thetaB = tb;
							contacts[icont].psiA = pa;
							contacts[icont].psiB = pb;
							contacts[icont].chi = chi;
							contacts[icont].d = dab;


							icont++; // Update contact index
						}
			}
	}
	else // Loop vs. Environment mode
	{
		//		fprintf(stderr,"\nseq2=\n");
		//		for(int i=0; i<nresloop; i++)
		//			fprintf(stderr,"%d %8.3f\n",seq2[i],coord2[9*i+3]);
		//		fprintf(stderr,"\n");

		// Environment left side vs. Loop
		for(int a = 0; a <= anchorNt; a++) // Screen all interaction frames
			for(int b = 0; b < nres2; b++) // Non-redundant screening
			{
				fa = &frames[a]; // environment
				fb = &frames2[b]; // loop

				// This one seems the fastest alternative... but all are very close...
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
					if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
						if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
						{
							dab = sqrtf(dab);

							if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
							{
								fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
								// exit(2);
							}


							frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);


							// Bonding/Non-bonding discrimination
							if( fa->ch != fb->ch ) // If residue indices belong to different chains
								bond = 0; // always non-bonding interaction between frames from different chains
							else
							{
								if( fb->i >= fa->i ) // Sanity check...
								{
									bond = fb->i - fa->i; // (bonding or non-bonding)
									if( bond > BONDING_THR)
										bond = 0; // non-bonding
								}
								else
								{
									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. Forcing exit!\n",fb->i,fa->i);
									exit(2);
								}
							}

							// When contact memory is over this allocates a new contacts block
							if(icont % CONTBLOCK == 0)
							{
								// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
								fflush(stdout);
								if( !(*p_contacts = contacts = (contact *) realloc(contacts, sizeof(contact) * (icont+CONTBLOCK))) )  // Reallocating to store new contacts
								{
									fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icont/powf(2,20));
									exit(1);
								}
							}

							// Store in RAM the computed data
							contacts[icont].ssA = 'N'; // Not computed
							contacts[icont].ssB = 'N';
							contacts[icont].seqA = seq[a]; // environment
							contacts[icont].seqB = seq2[b]; // loop
							contacts[icont].froitA = fa->t;
							contacts[icont].froitB = fb->t;
							//							contacts[icont].froitA = (unsigned int) fa->i;
							//							contacts[icont].froitB = (unsigned int) fb->i;
							contacts[icont].sd = bond;
							contacts[icont].thetaA = ta;
							contacts[icont].thetaB = tb;
							contacts[icont].psiA = pa;
							contacts[icont].psiB = pb;
							contacts[icont].chi = chi;
							contacts[icont].d = dab;

							//							if(a==1)
							//							{
							//								fprintf(stderr,"%4d %4d %4d %4d %4d %8.3f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
							//										contacts[icont].sd,  contacts[icont].froitA,  contacts[icont].froitB,  a, b,  contacts[icont].d,
							//												  contacts[icont].thetaA, contacts[icont].thetaB, contacts[icont].psiA, contacts[icont].psiB,contacts[icont].chi);
							//							}

							//							fprintf(stderr,"Env_vs_Loop %4d  A= %8.3f %8.3f %8.3f  B= %8.3f %8.3f %8.3f  a= %3d  b= %3d  i= %3d %3d  seq= %2d %2d  sd= %d  d= %8.3f\n",
							//									icont, fa->p[0], fa->p[1], fa->p[2], fb->p[0], fb->p[1], fb->p[2], a,b, fa->i, fb->i,
							//									contacts[icont].seqA, contacts[icont].seqB, contacts[icont].sd, contacts[icont].d);

							icont++; // Update contact index
						}
			}

		// Loop vs. Environment right side
		for(int a = 0; a < nres2; a++) // Non-redundant screening
			for(int b = anchorNt+1; b < nres ; b++) // Screen all interaction frames
			{
				fa = &frames2[a]; // loop
				fb = &frames[b]; // environment

				// This one seems the fastest alternative... but all are very close...
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
					if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
						if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
						{
							dab = sqrtf(dab);

							if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
							{
								fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
								// exit(2);
							}

							frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);

							//							bond = (int)fabsf(fa->i - fb->i);
							//							if( bond > BONDING_THR)
							//								bond = 0;

							// Bonding/Non-bonding discrimination
							if( fa->ch != fb->ch ) // If residue indices belong to different chains
								bond = 0; // always non-bonding interaction between frames from different chains
							else
							{
								if( fb->i >= fa->i ) // Sanity check...
								{
									bond = fb->i - fa->i; // (bonding or non-bonding)
									if( bond > BONDING_THR)
										bond = 0; // non-bonding
								}
								else
								{
									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. Forcing exit!\n",fb->i,fa->i);
									exit(2);
								}
							}

							// When contact memory is over this allocates a new contacts block
							if(icont % CONTBLOCK == 0)
							{
								// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
								fflush(stdout);
								if( !(*p_contacts = contacts = (contact *) realloc(contacts, sizeof(contact) * (icont+CONTBLOCK))) )  // Reallocating to store new contacts
								{
									fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icont/powf(2,20));
									exit(1);
								}
							}

							// Store in RAM the computed data
							contacts[icont].ssA = 'N'; // Not computed
							contacts[icont].ssB = 'N';
							contacts[icont].seqA = seq2[a]; // loop
							contacts[icont].seqB = seq[b]; // environment
							contacts[icont].froitA = fa->t;
							contacts[icont].froitB = fb->t;
							//							contacts[icont].froitA = (unsigned int) fa->i;
							//							contacts[icont].froitB = (unsigned int) fb->i;
							contacts[icont].sd = bond;
							contacts[icont].thetaA = ta;
							contacts[icont].thetaB = tb;
							contacts[icont].psiA = pa;
							contacts[icont].psiB = pb;
							contacts[icont].chi = chi;
							contacts[icont].d = dab;

							//							fprintf(stderr,"Loop_vs_Env %4d  A= %8.3f %8.3f %8.3f  B= %8.3f %8.3f %8.3f  i= %3d %3d  seq= %2d %2d  sd= %d\n",
							//									icont, fa->p[0], fa->p[1], fa->p[2], fb->p[0], fb->p[1], fb->p[2], fa->i, fb->i,
							//									contacts[icont].seqA, contacts[icont].seqB, contacts[icont].sd);

							icont++; // Update contact index
						}
			}

		//		free(frames2);
	}

	//	printf("%s> %d contacts found (%.1f/residue at %.0f A cutoff), %.1f%% from %d available.\n",
	//			prog,icont,(float)icont/num_res,cutoff,100*(float)icont/(num_res*(num_res-1)/2),num_res*(num_res-1)/2);
	//	printf("%s> Memory used for contacts: %ld bytes (%.1f MB)\n\n", prog, sizeof(contact)*icont, (float)sizeof(contact)*icont/powf(2,20));

	// Free memory
	//	free(frames);

	return icont; // Return the total number of contacts upon successful exit
}

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
int contactCoordM(contact **p_contactsM, float cutoff, frame *frames, int nres, int *seq, int posM, char chainM)
{
	char prog[] = "contactCoord";

	// Allocating initial memory for contacts
	contact *contactsM; // Contacts array
	if(*p_contactsM == NULL) // Contact list is empty (first PDB typically)
		*p_contactsM = contactsM = (contact *) malloc( sizeof(contact) * CONTBLOCK); // Allocating first memory block to store contacts
	else
		contactsM = *p_contactsM; // already "mallocated"

	// COMPUTING CONTACT DATA
	// Screen all inter-residue contacts
	float dab,ta,tb,pa,pb,chi; // inter-frame internal coordinates
	int bond;
	int icontM = 0; // contact index for current protein

	frame *fa,*fb;
	float cutoff2;
	cutoff2 = cutoff * cutoff; // cutoff squared for efficient distance evaluation

	for(int a = 0; a < nres; a++) // Screen all interaction frames
		for(int b = a+1; b < nres; b++) // Non-redundant screening
		{
			fa = &frames[a];
			fb = &frames[b];

			// Chapa to select just the mutant contacts
			if ( ((fa->ch == chainM) and (fa->i == posM )) or ((fb->ch == chainM) and (fb->i == posM )) )
			{
				// This one seems the fastest alternative... but all are very close...
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
					if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
						if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
						{
							dab = sqrtf(dab);

							if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
							{
								fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
								// exit(2);
							}

							frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);

							// Bonding/Non-bonding discrimination
							if( fa->ch != fb->ch ) // If residue indices belong to different chains
								bond = 0; // always non-bonding interaction between frames from different chains
							else
							{
								if( fb->i >= fa->i ) // Sanity check...
								{
									bond = fb->i - fa->i; // (bonding or non-bonding)
									if( bond > BONDING_THR)
										bond = 0; // non-bonding
								}
								else
								{
									//									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. Forcing exit!\n",fb->i,fa->i);
									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. (fb->ch= %c  fa->ch= %c) Forcing exit!\n",fb->i,fa->i,fb->ch,fa->ch);
									exit(2);
								}
							}


							// printf("%d %d %d %d\n", icontM, fb->i, fa->i, posM);
							// When contact memory is over this allocates a new contacts block
							if(icontM % CONTBLOCK == 0)
							{
								// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
								fflush(stdout);
								if( !(*p_contactsM = contactsM = (contact *) realloc(contactsM, sizeof(contact) * (icontM+CONTBLOCK))) )  // Reallocating to store new contacts
								{
									fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icontM/powf(2,20));
									exit(1);
								}
							}

							// Store in RAM the computed data
							// contactsM[icontM].ssA = 'N'; // Not computed
							// contactsM[icontM].ssB = 'N';
							contactsM[icontM].fai = fa->i;  // added pablo 2020
							contactsM[icontM].fbi = fb->i;  // added pablo 2020
							contactsM[icontM].fach = fa->ch;  // added pablo 2020
							contactsM[icontM].fbch = fb->ch;  // added pablo 2020
							contactsM[icontM].seqA = seq[a];
							contactsM[icontM].seqB = seq[b];
							contactsM[icontM].froitA = fa->t;
							contactsM[icontM].froitB = fb->t;
							// printf ("%c %c\n",contactsM[icontM].froitA,contactsM[icontM].froitB);
							//							contacts[icont].froitA = (unsigned int) fa->i;
							//							contacts[icont].froitB = (unsigned int) fb->i;
							contactsM[icontM].sd = bond;
							contactsM[icontM].thetaA = ta;
							contactsM[icontM].thetaB = tb;
							contactsM[icontM].psiA = pa;
							contactsM[icontM].psiB = pb;
							contactsM[icontM].chi = chi;
							contactsM[icontM].d = dab;
							icontM++; // Update contact index

						}
			}
		}

	return icontM; // Return the total number of contacts upon successful exit
}

// MON: Multi Point Mutations? Yes.
int contactCoordML(contact **p_contactsM, float cutoff, frame *frames, int nres, int *seq, int *posM, char *chainM, int nmut)
{
	char prog[] = "contactCoord";

	// Allocating initial memory for contacts
	contact *contactsM; // Contacts array
	if(*p_contactsM == NULL) // Contact list is empty (first PDB typically)
		*p_contactsM = contactsM = (contact *) malloc( sizeof(contact) * CONTBLOCK); // Allocating first memory block to store contacts
	else
		contactsM = *p_contactsM; // already "mallocated"

	// COMPUTING CONTACT DATA
	// Screen all inter-residue contacts
	float dab,ta,tb,pa,pb,chi; // inter-frame internal coordinates
	int bond;
	int icontM = 0; // contact index for current protein

	frame *fa,*fb;
	float cutoff2;
	cutoff2 = cutoff * cutoff; // cutoff squared for efficient distance evaluation

	for(int a = 0; a < nres; a++) // Screen all interaction frames
		for(int b = a+1; b < nres; b++) // Non-redundant screening
		{
			fa = &frames[a];
			fb = &frames[b];

			// Chapa to select just the mutant contacts
			for(int m = 0; m < nmut; m++) // MON: "<" instead of "<=" required
				if ( ((fa->ch == chainM[m]) and (fa->i == posM[m] )) or ((fb->ch == chainM[m]) and (fb->i == posM[m] )) )
				{
					// This one seems the fastest alternative... but all are very close...
					if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
						if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
							if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
							{
								dab = sqrtf(dab);

								if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
								{
									fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
									// exit(2);
								}

								frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);

								// Bonding/Non-bonding discrimination
								if( fa->ch != fb->ch ) // If residue indices belong to different chains
									bond = 0; // always non-bonding interaction between frames from different chains
								else
								{
									if( fb->i >= fa->i ) // Sanity check...
									{
										bond = fb->i - fa->i; // (bonding or non-bonding)
										if( bond > BONDING_THR)
											bond = 0; // non-bonding
									}
									else
									{
										//									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. Forcing exit!\n",fb->i,fa->i);
										fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. (fb->ch= %c  fa->ch= %c) Forcing exit!\n",fb->i,fa->i,fb->ch,fa->ch);
										exit(2);
									}
								}


								// printf("%d %d %d %d\n", icontM, fb->i, fa->i, posM);
								// When contact memory is over this allocates a new contacts block
								if(icontM % CONTBLOCK == 0)
								{
									// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
									fflush(stdout);
									if( !(*p_contactsM = contactsM = (contact *) realloc(contactsM, sizeof(contact) * (icontM+CONTBLOCK))) )  // Reallocating to store new contacts
									{
										fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icontM/powf(2,20));
										exit(1);
									}
								}

								// Store in RAM the computed data
								// contactsM[icontM].ssA = 'N'; // Not computed
								// contactsM[icontM].ssB = 'N';
								contactsM[icontM].fai = fa->i;  // added pablo 2020
								contactsM[icontM].fbi = fb->i;  // added pablo 2020
								contactsM[icontM].fach = fa->ch;  // added pablo 2020
								contactsM[icontM].fbch = fb->ch;  // added pablo 2020
								contactsM[icontM].seqA = seq[a];
								contactsM[icontM].seqB = seq[b];
								contactsM[icontM].froitA = fa->t;
								contactsM[icontM].froitB = fb->t;
								// printf ("%c %c\n",contactsM[icontM].froitA,contactsM[icontM].froitB);
								//							contacts[icont].froitA = (unsigned int) fa->i;
								//							contacts[icont].froitB = (unsigned int) fb->i;
								contactsM[icontM].sd = bond;
								contactsM[icontM].thetaA = ta;
								contactsM[icontM].thetaB = tb;
								contactsM[icontM].psiA = pa;
								contactsM[icontM].psiB = pb;
								contactsM[icontM].chi = chi;
								contactsM[icontM].d = dab;
								icontM++; // Update contact index

							}
					break;
				}
		}

	return icontM; // Return the total number of contacts upon successful exit
}



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
		int anchorNt, int nresloop, float *coord2, int *seq2, int *resnumsloop, char *reschainsloop)
{
	//	fprintf(stderr,"Hi, I'm contactPDB\n");

	char prog[] = "contactCoord";

	//	fprintf(stderr,"nintres= %d  num_atoms= %d  num_res= %d\n",nintres,num_atoms,num_res);
	//	pdb->writePDB("contactPDB.pdb");
	//	fprintf(stderr,"Written contactPDB.pdb\n");
	// exit(0);

	// Maximum number of interaction frames per residue (to allocate memory for the maximum number of FRames of Interest)
	//	for(int i=0; i<30;i++)
	//		fprintf(stderr,"iaa[%d] = %d\n",i,iaa[i]);
	//	fprintf(stderr,"nintres= %d  maxint= %d  num_atoms= %d  num_res= %d\n",nintres,maxint,num_atoms,num_res);

	// Allocate memory for residue frames
	frame *frames,*frames2;
	frames = (frame *) malloc(sizeof(frame) * num_res); // allocate the maximum possible

	// The residues sequence must be specified by unique identifier (as defined in ResIni.h)
	//	char *seq;
	//	seq = (char *) malloc( sizeof(char) * num_res);

	// Some variables...
	Tcoor r12,r13;

	//	fprintf(stderr,"seq=\n");
	//	for(int i=0; i<num_res; i++)
	//		fprintf(stderr,"%d ",seq[i]);
	//	fprintf(stderr,"\n");
	//	fprintf(stderr,"anchorNt %d %d", anchorNt, nresloop);

	// COMPUTING MACROMOLECULE FRAMES (FULL OR ENVIRONMENT)
	frame *fram; // current frame pointer
	for(int i=0; i<num_res; i++) // Screen coordinates of 1st macromolecule
	{
		fram = &frames[i];

		//		if(i>anchorNt)
		//			fram->i = i + nresloop; // store residue index (required later)
		//		else
		//			fram->i = i; // store residue index (required later)
		fram->i = resnums[i]; // using PDB's residue numeration
		fram->ch = reschains[i]; // using PDB's chain-ids

		fram->t = 0; // store frame type (required later)
		// fprintf(stderr,"fram->i= %d  fram->t= %d\n",fram->i,fram->t);

		// Compute residue frames (Eq.1 from GOAP's paper Zhou and Skolnik 2011. A1=N, A=CA, A2=C)

		//		N[0]=*(coord+9*i);
		//		N[1]=*(coord+9*i+1);
		//		N[2]=*(coord+9*i+2);
		//
		//		CA[0]=*(coord+9*i+3);
		//		CA[1]=*(coord+9*i+4);
		//		CA[2]=*(coord+9*i+5);
		//
		//		C[0]=*(coord+9*i+6);
		//		C[1]=*(coord+9*i+7);
		//	    C[2]=*(coord+9*i+8);
		//
		//	    diff3D(N,CA,r12);
		//      diff3D(C,CA,r13);

		diff3D(coord+9*i,coord+9*i+3,r12);   // diff3D(N,CA,r12);
		diff3D(coord+9*i+6,coord+9*i+3,r13); // diff3D(C,CA,r13);
		sum3D(r12, r13, fram->vz);
		normalize3D(fram->vz); // Vz computed
		cross3D(fram->vz,r13,fram->vy);
		normalize3D(fram->vy); // Vy computed
		cross3D(fram->vy,fram->vz,fram->vx); // Vx computed


		// Store interaction distances (CA-positions) into frames array
		fram->p[0] = coord[9*i+3];
		fram->p[1] = coord[9*i+4];
		fram->p[2] = coord[9*i+5];

		//		fprintf(stderr,"i=%4d   fram->p= %f %f %f   seq= %d\n",i,fram->p[0],fram->p[1],fram->p[2], seq[i]);

		//		if(i==1)
		//		{
		//			fprintf(stderr,"i=%4d   fram->p= %8.3f %8.3f %8.3f  seq= %d  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f \n",i,
		//					fram->p[0],fram->p[1],fram->p[2],seq[i],
		//					fram->vx[0],fram->vx[1],fram->vx[2],
		//					fram->vy[0],fram->vy[1],fram->vy[2],
		//					fram->vz[0],fram->vz[1],fram->vz[2] );
		//		}
	}

	// Allocating initial memory for contacts
	contact *contacts; // Contacts array
	if(*p_contacts == NULL) // Contact list is empty (first PDB typically)
		*p_contacts = contacts = (contact *) malloc( sizeof(contact) * CONTBLOCK); // Allocating first memory block to store contacts
	else
		contacts = *p_contacts; // already "mallocated"

	// COMPUTING CONTACT DATA
	// Screen all inter-residue contacts
	float dab,ta,tb,pa,pb,chi; // inter-frame internal coordinates
	int bond;
	int icont=0; // contact index for current protein
	frame *fa,*fb;
	float cutoff2;
	cutoff2 = cutoff * cutoff; // cutoff squared for efficient distance evaluation

	if(coord2 == NULL) // Standard mode (all vs. all mode)
	{
		for(int a = 0; a < num_res; a++) // Screen all interaction frames
			for(int b = a+1; b < num_res; b++) // Non-redundant screening
			{
				fa = &frames[a];
				fb = &frames[b];

				// This one seems the fastest alternative... but all are very close...
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
					if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
						if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
						{
							dab = sqrtf(dab);

							if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
							{
								fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
								// exit(2);
							}

							frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);

							//							bond = (int)fabsf(fa->i - fb->i); // Mon: abs is not necessary, isn't it?
							//							if( bond > BONDING_THR)
							//								bond = 0;

							// Bonding/Non-bonding discrimination
							if( fa->ch != fb->ch ) // If residue indices belong to different chains
								bond = 0; // always non-bonding interaction between frames from different chains
							else
							{
								if( fb->i >= fa->i ) // Sanity check...
								{
									bond = fb->i - fa->i; // (bonding or non-bonding)
									if( bond > BONDING_THR)
										bond = 0; // non-bonding
								}
								else
								{
									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. (fb->ch= %c  fa->ch= %c) Forcing exit!\n",fb->i,fa->i,fb->ch,fa->ch);
									exit(2);
								}
							}

							// When contact memory is over this allocates a new contacts block
							if(icont % CONTBLOCK == 0)
							{
								// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
								fflush(stdout);
								if( !(*p_contacts = contacts = (contact *) realloc(contacts, sizeof(contact) * (icont+CONTBLOCK))) )  // Reallocating to store new contacts
								{
									fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icont/powf(2,20));
									exit(1);
								}
							}

							// Store in RAM the computed data
							contacts[icont].ssA = 'N'; // Not computed
							contacts[icont].ssB = 'N';
							contacts[icont].seqA = seq[a];
							contacts[icont].seqB = seq[b];
							contacts[icont].froitA = fa->t;
							contacts[icont].froitB = fb->t;
							//							contacts[icont].froitA = (unsigned int) fa->i;
							//							contacts[icont].froitB = (unsigned int) fb->i;
							contacts[icont].sd = bond;
							contacts[icont].thetaA = ta;
							contacts[icont].thetaB = tb;
							contacts[icont].psiA = pa;
							contacts[icont].psiB = pb;
							contacts[icont].chi = chi;
							contacts[icont].d = dab;

							//							if(a==1 && b==83)
							//							{
							//								fprintf(stderr,"%4d %4d %4d %4d %4d %8.3f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
							//										contacts[icont].sd,  contacts[icont].froitA,  contacts[icont].froitB,  a, b,  contacts[icont].d,
							//												  contacts[icont].thetaA, contacts[icont].thetaB, contacts[icont].psiA, contacts[icont].psiB,contacts[icont].chi);
							//
							//								fprintf(stderr,"a=%4d  OK     fram->p= %8.3f %8.3f %8.3f  seq= %d  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f \n",
							//										a,
							//										fa->p[0], fa->p[1], fa->p[2],seq[a],
							//										fa->vx[0],fa->vx[1],fa->vx[2],
							//										fa->vy[0],fa->vy[1],fa->vy[2],
							//										fa->vz[0],fa->vz[1],fa->vz[2] );
							//								fprintf(stderr,"b=%4d  OK     fram->p= %8.3f %8.3f %8.3f  seq= %d  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f \n",
							//										b,
							//										fb->p[0],fb->p[1],fb->p[2],seq[b],
							//										fb->vx[0],fb->vx[1],fb->vx[2],
							//										fb->vy[0],fb->vy[1],fb->vy[2],
							//										fb->vz[0],fb->vz[1],fb->vz[2] );
							//							}

							//							fprintf(stderr,"Env_vs_Env  %4d  A= %8.3f %8.3f %8.3f  B= %8.3f %8.3f %8.3f  i= %3d %3d  seq= %2d %2d  sd= %d\n",
							//									icont, fa->p[0], fa->p[1], fa->p[2], fb->p[0], fb->p[1], fb->p[2], fa->i, fb->i,
							//									contacts[icont].seqA, contacts[icont].seqB, contacts[icont].sd);

							icont++; // Update contact index
						}
			}
	}
	else // Loop vs. Environment mode
	{
		//		fprintf(stderr,"\nseq2=\n");
		//		for(int i=0; i<nresloop; i++)
		//			fprintf(stderr,"%d %8.3f\n",seq2[i],coord2[9*i+3]);
		//		fprintf(stderr,"\n");

		frames2 = (frame *) malloc(sizeof(frame) * nresloop); // allocate loop frames

		// COMPUTING LOOP FRAMES
		for(int i=0; i<nresloop; i++) // Screen coordinates of 1st macromolecule
		{
			fram = &frames2[i];

			//			fram->i = i + anchorNt + 1; // store residue index (required later)

			fram->i = resnumsloop[i]; // using PDB's residue numeration
			fram->ch = reschainsloop[i]; // using PDB's chain-ids

			fram->t = 0; // store frame type (required later)
			//			fprintf(stderr,"fram->i= %d  fram->t= %d\n",fram->i,fram->t);

			// Compute residue frames (Eq.1 from GOAP's paper Zhou and Skolnik 2011. A1=N, A=CA, A2=C)
			diff3D(coord2+9*i,coord2+9*i+3,r12);   // diff3D(N,CA,r12);
			diff3D(coord2+9*i+6,coord2+9*i+3,r13); // diff3D(C,CA,r13);
			sum3D(r12, r13, fram->vz);
			normalize3D(fram->vz); // Vz computed
			cross3D(fram->vz,r13,fram->vy);
			normalize3D(fram->vy); // Vy computed
			cross3D(fram->vy,fram->vz,fram->vx); // Vx computed

			// Store interaction distances (CA-positions) into frames array
			fram->p[0] = coord2[9*i+3];
			fram->p[1] = coord2[9*i+4];
			fram->p[2] = coord2[9*i+5];
			// fprintf(stderr,"i=%4d   fram->p= %f %f %f\n",i,fram->p[0],fram->p[1],fram->p[2]);
		}

		// Environment left side vs. Loop
		for(int a = 0; a <= anchorNt; a++) // Screen all interaction frames
			for(int b = 0; b < nresloop; b++) // Non-redundant screening
			{
				fa = &frames[a]; // environment
				fb = &frames2[b]; // loop

				// This one seems the fastest alternative... but all are very close...
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
					if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
						if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
						{
							dab = sqrtf(dab);

							if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
							{
								fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
								// exit(2);
							}

							//							if(a==1)
							//							{
							//								fprintf(stderr,"a=%4d  Before  fram->p= %8.3f %8.3f %8.3f  seq= %d  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f \n",
							//										a,
							//										fa->p[0],fa->p[1],fa->p[2],seq[a],
							//										fa->vx[0],fa->vx[1],fa->vx[2],
							//										fa->vy[0],fa->vy[1],fa->vy[2],
							//										fa->vz[0],fa->vz[1],fa->vz[2] );
							//								fprintf(stderr,"b=%4d  Before  fram->p= %8.3f %8.3f %8.3f  seq= %d  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f \n",
							//										b,
							//										fb->p[0],fb->p[1],fb->p[2],seq2[b],
							//										fb->vx[0],fb->vx[1],fb->vx[2],
							//										fb->vy[0],fb->vy[1],fb->vy[2],
							//										fb->vz[0],fb->vz[1],fb->vz[2] );
							//							}

							frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);

							//							if(a==1)
							//							{
							//								fprintf(stderr,"a=%4d  After  fram->p= %8.3f %8.3f %8.3f  seq= %d  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f \n",
							//										a,
							//										fa->p[0],fa->p[1],fa->p[2],seq[a],
							//										fa->vx[0],fa->vx[1],fa->vx[2],
							//										fa->vy[0],fa->vy[1],fa->vy[2],
							//										fa->vz[0],fa->vz[1],fa->vz[2] );
							//
							//								fprintf(stderr,"b=%4d  After  fram->p= %8.3f %8.3f %8.3f  seq= %d  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f \n",
							//										b,
							//										fb->p[0],fb->p[1],fb->p[2],seq2[b],
							//										fb->vx[0],fb->vx[1],fb->vx[2],
							//										fb->vy[0],fb->vy[1],fb->vy[2],
							//										fb->vz[0],fb->vz[1],fb->vz[2] );
							//							}

							//							bond = (int)fabsf(fa->i - fb->i);
							//							if( bond > BONDING_THR)
							//								bond = 0;

							// Bonding/Non-bonding discrimination
							if( fa->ch != fb->ch ) // If residue indices belong to different chains
								bond = 0; // always non-bonding interaction between frames from different chains
							else
							{
								if( fb->i >= fa->i ) // Sanity check...
								{
									bond = fb->i - fa->i; // (bonding or non-bonding)
									if( bond > BONDING_THR)
										bond = 0; // non-bonding
								}
								else
								{
									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. Forcing exit!\n",fb->i,fa->i);
									exit(2);
								}
							}

							// When contact memory is over this allocates a new contacts block
							if(icont % CONTBLOCK == 0)
							{
								// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
								fflush(stdout);
								if( !(*p_contacts = contacts = (contact *) realloc(contacts, sizeof(contact) * (icont+CONTBLOCK))) )  // Reallocating to store new contacts
								{
									fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icont/powf(2,20));
									exit(1);
								}
							}

							// Store in RAM the computed data
							contacts[icont].ssA = 'N'; // Not computed
							contacts[icont].ssB = 'N';
							contacts[icont].seqA = seq[a]; // environment
							contacts[icont].seqB = seq2[b]; // loop
							contacts[icont].froitA = fa->t;
							contacts[icont].froitB = fb->t;
							//							contacts[icont].froitA = (unsigned int) fa->i;
							//							contacts[icont].froitB = (unsigned int) fb->i;
							contacts[icont].sd = bond;
							contacts[icont].thetaA = ta;
							contacts[icont].thetaB = tb;
							contacts[icont].psiA = pa;
							contacts[icont].psiB = pb;
							contacts[icont].chi = chi;
							contacts[icont].d = dab;

							//							if(a==1)
							//							{
							//								fprintf(stderr,"%4d %4d %4d %4d %4d %8.3f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
							//										contacts[icont].sd,  contacts[icont].froitA,  contacts[icont].froitB,  a, b,  contacts[icont].d,
							//												  contacts[icont].thetaA, contacts[icont].thetaB, contacts[icont].psiA, contacts[icont].psiB,contacts[icont].chi);
							//							}

							//							fprintf(stderr,"Env_vs_Loop %4d  A= %8.3f %8.3f %8.3f  B= %8.3f %8.3f %8.3f  a= %3d  b= %3d  i= %3d %3d  seq= %2d %2d  sd= %d  d= %8.3f\n",
							//									icont, fa->p[0], fa->p[1], fa->p[2], fb->p[0], fb->p[1], fb->p[2], a,b, fa->i, fb->i,
							//									contacts[icont].seqA, contacts[icont].seqB, contacts[icont].sd, contacts[icont].d);

							icont++; // Update contact index
						}
			}

		// Loop vs. Environment right side
		for(int a = 0; a < nresloop; a++) // Non-redundant screening
			for(int b = anchorNt+1; b < num_res ; b++) // Screen all interaction frames
			{
				fa = &frames2[a]; // loop
				fb = &frames[b]; // environment

				// This one seems the fastest alternative... but all are very close...
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
					if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
						if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
						{
							dab = sqrtf(dab);

							if(dab < DIST2SMALL) // To see crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
							{
								fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
								// exit(2);
							}

							frames2ic(fa,fb,&ta,&tb,&pa,&pb,&chi);

							//							bond = (int)fabsf(fa->i - fb->i);
							//							if( bond > BONDING_THR)
							//								bond = 0;

							// Bonding/Non-bonding discrimination
							if( fa->ch != fb->ch ) // If residue indices belong to different chains
								bond = 0; // always non-bonding interaction between frames from different chains
							else
							{
								if( fb->i >= fa->i ) // Sanity check...
								{
									bond = fb->i - fa->i; // (bonding or non-bonding)
									if( bond > BONDING_THR)
										bond = 0; // non-bonding
								}
								else
								{
									fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. Forcing exit!\n",fb->i,fa->i);
									exit(2);
								}
							}

							// When contact memory is over this allocates a new contacts block
							if(icont % CONTBLOCK == 0)
							{
								// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
								fflush(stdout);
								if( !(*p_contacts = contacts = (contact *) realloc(contacts, sizeof(contact) * (icont+CONTBLOCK))) )  // Reallocating to store new contacts
								{
									fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icont/powf(2,20));
									exit(1);
								}
							}

							// Store in RAM the computed data
							contacts[icont].ssA = 'N'; // Not computed
							contacts[icont].ssB = 'N';
							contacts[icont].seqA = seq2[a]; // loop
							contacts[icont].seqB = seq[b]; // environment
							contacts[icont].froitA = fa->t;
							contacts[icont].froitB = fb->t;
							//							contacts[icont].froitA = (unsigned int) fa->i;
							//							contacts[icont].froitB = (unsigned int) fb->i;
							contacts[icont].sd = bond;
							contacts[icont].thetaA = ta;
							contacts[icont].thetaB = tb;
							contacts[icont].psiA = pa;
							contacts[icont].psiB = pb;
							contacts[icont].chi = chi;
							contacts[icont].d = dab;

							//							fprintf(stderr,"Loop_vs_Env %4d  A= %8.3f %8.3f %8.3f  B= %8.3f %8.3f %8.3f  i= %3d %3d  seq= %2d %2d  sd= %d\n",
							//									icont, fa->p[0], fa->p[1], fa->p[2], fb->p[0], fb->p[1], fb->p[2], fa->i, fb->i,
							//									contacts[icont].seqA, contacts[icont].seqB, contacts[icont].sd);

							icont++; // Update contact index
						}
			}

		free(frames2);
	}

	//	printf("%s> %d contacts found (%.1f/residue at %.0f A cutoff), %.1f%% from %d available.\n",
	//			prog,icont,(float)icont/num_res,cutoff,100*(float)icont/(num_res*(num_res-1)/2),num_res*(num_res-1)/2);
	//	printf("%s> Memory used for contacts: %ld bytes (%.1f MB)\n\n", prog, sizeof(contact)*icont, (float)sizeof(contact)*icont/powf(2,20));

	// Free memory
	free(frames);

	return icont; // Return the total number of contacts upon successful exit
}


// Compute Order_AVE (4D) contacts for energy evaluation (automatic memory allocation)
//  pdb         --> Macromolecule
//  p_contacts  --> Pointer to the contacts array (if NULL, it allocates memory, otherwise memory pre-allocation expected)
//  cutoff      --> Maximum distance cutoff
//  iaa         --> Residues mapping into the mapping table
//  mapping     --> Mapping table with the atom indices to define the orthogonal framework and obtain the pairwise distances
//  nintres     --> Number of interacting residues
int contactPDB_OA(Macromolecule *pdb, contact **p_contacts, float cutoff, char *iaa, char **mapping, int nintres)
{
	char prog[] = "contactPDB_OA";
	//int num_atoms = pdb->get_num_atoms();
	int num_res = pdb->get_num_fragments();

	// Maximum number of interaction frames per residue (to allocate memory for the maximum number of FRames of Interest)
	int maxint = 0;
	for(int i=0; i<nintres; i++)
	{
		if(mapping[i][0] > maxint)
			maxint = mapping[i][0];
	}

	// Allocate memory for residue frames
	frame *frames;
	frames = (frame *) malloc(sizeof(frame) * num_res * maxint); // allocate the maximum possible

	// The residues sequence must be specified by unique identifier (as defined in ResIni.h)
	char *seq;
	seq = (char *) malloc( sizeof(char) * num_res * maxint);

	// Some variables...
	pdbIter *iter_ch,*iter_res,*iter_atom;
	pdbIter *iter_prev; // Previous residue iterator
	Residue *res;
	Chain *ch;
	Tcoor N,CA;

	// COMPUTING RESIDUE FRAMES
	iter_ch = new pdbIter( pdb ); // iter to screen segments
	//	int ires = 0; // Residue index
	int ifram = 0; // Frame index
	char *imap; // current interactions Mapping
	char ichain; // current chain
	frame *fram; // current frame pointer
	for ( iter_ch->pos_chain = 0; !iter_ch->gend_chain(); iter_ch->next_chain() ) // Iter chains to compute sequential distances fine
	{
		ch = (Chain *) iter_ch->get_chain();
		iter_res = new pdbIter( ch );

		ichain = (ch->getName())[0]; // Current chain id (1-letter)

		for ( iter_res->pos_fragment = 0; !iter_res->gend_fragment(); iter_res->next_fragment() ) // iter frags
		{
			res = ( Residue * ) iter_res->get_fragment();
			iter_atom = new pdbIter( res ); // screens residue atoms

			imap = mapping[ (int)iaa[ (int)res->fragid ] ];

			for(int t=0; t<imap[0]; t++) // screen frame types according to mapping
			{
				if(imap[1+4*t+1] < 0 ) // If reference atom (N) comes from previous residue
				{
					if(iter_res->pos_fragment > 0) // If there exists previous residue
					{
						iter_prev->pos_atom = -imap[1+4*t+1]; // atom index from previous residue must be sign reversed
						(iter_prev->get_atom())->getPosition(N); // N's position (from previous residue)
					}
					else // Skip frame type
						continue;
				}
				else
				{
					iter_atom->pos_atom = imap[1+4*t+1];
					(iter_atom->get_atom())->getPosition(N); // N's position (from current residue)
				}

				iter_atom->pos_atom = imap[1+4*t+2]; // CA's position
				(iter_atom->get_atom())->getPosition(CA);


				seq[ifram] = res->fragid; // Get residue identifier

				fram = &frames[ifram];

				//				fram->i = ires; // store residue index (required later)
				fram->i = res->getIdNumber(); // Get PDB numeration (required for correct bonding/non-bonding discrimination)
				fram->t = t; // store frame type (required later)
				fram->ch = ichain; // store chain id (required later)

				// Store reference atom position (N)
				fram->vx[0] = N[0];
				fram->vx[1] = N[1];
				fram->vx[2] = N[2];

				// Store interaction distances (CA-positions) into frames array
				fram->p[0] = CA[0];
				fram->p[1] = CA[1];
				fram->p[2] = CA[2];

				// The other "frame" fields are ignored for Order_AVE 4D

				ifram++; // update frame index
			}

			delete iter_atom;
			//			ires++; // update current residue index

			if(iter_res->pos_fragment > 0) // Non first residue
				delete iter_prev; // delete previous residue iterator
			iter_prev = new pdbIter( res ); // screens residue atoms
		}
		delete iter_res;
	}
	delete iter_ch;
	// frames = (frame *) realloc(frames, sizeof(frame) * ifram ); // Reallocating to exactly store "ifram" interaction frames.
	// seq = (char *) realloc(seq, sizeof(char) * ifram ); // Reallocating to exactly store "ifram" sequence chars.

	// Allocating initial memory for contacts
	contact *contacts; // Contacts array
	if(*p_contacts == NULL) // Contact list is empty (first PDB typically)
		*p_contacts = contacts = (contact *) malloc( sizeof(contact) * CONTBLOCK); // Allocating first memory block to store contacts
	else
		contacts = *p_contacts; // already "mallocated"

	// COMPUTING CONTACT DATA
	// Screen all inter-residue contacts
	float dab,ta,tb,chi; // inter-frame internal coordinates
	int bond;
	int icont=0; // contact index for current protein
	frame *fa,*fb;
	float cutoff2;
	cutoff2 = cutoff * cutoff; // cutoff squared for efficient distance evaluation

	for(int a = 0; a < ifram; a++) // Screen all interaction frames
		for(int b = a+1; b < ifram; b++) // Non-redundant screening
		{
			fa = &frames[a];
			fb = &frames[b];

			// This one seems the fastest alternative... but all are very close...
			if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) ) < cutoff2 )
				if( ( (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) ) < cutoff2 )
					if( ( dab = (fb->p[0]-fa->p[0])*(fb->p[0]-fa->p[0]) + (fb->p[1]-fa->p[1])*(fb->p[1]-fa->p[1]) + (fb->p[2]-fa->p[2])*(fb->p[2]-fa->p[2]) ) < cutoff2 )
					{
						dab = sqrtf(dab);

						if(dab < DIST2SMALL) // To detect crap PDBs (ITASSER benchmark set II has some 1.5 A CA-CA distance...not too bad)
						{
							fprintf(stderr,"contactPDB> Warning, very small inter-CA distance %f (< %.2f A) in PDB between residues %d and %d.\n",dab,DIST2SMALL*1.0, a,b);
							// exit(2);
						}

						// Computes the 3 internal coordinates for the Order_AVE like potentials between (i,ir) and (j,jr) interacting frames
						frames2ic_OA(fa->p, fb->p, fa->vx, fb->vx, &ta, &tb, &chi);

						// Bonding/Non-bonding discrimination
						if( fa->ch != fb->ch ) // If residue indices belong to different chains
							bond = 0; // always non-bonding interaction between frames from different chains
						else
						{
							if( fb->i >= fa->i ) // Sanity check...
							{
								bond = fb->i - fa->i; // (bonding or non-bonding)
								if( bond > BONDING_THR)
									bond = 0; // non-bonding
							}
							else
							{
								fprintf(stderr,"Some PDB numeration error was found... fb->i= %d < fa->i= %d. (fb->ch= %c  fa->ch= %c) Forcing exit!\n",fb->i,fa->i,fb->ch,fa->ch);
								exit(2);
							}
						}

						// When contact memory is over this allocates a new contacts block
						if(icont % CONTBLOCK == 0)
						{
							// printf("\r%s> Contacts added so far: %6d k (%.0f MB)\n", prog, icont/1000, (float)sizeof(contact)*icont/powf(2,20));
							fflush(stdout);
							if( !(*p_contacts = contacts = (contact *) realloc(contacts, sizeof(contact) * (icont+CONTBLOCK))) )  // Reallocating to store new contacts
							{
								fprintf(stderr,"%s> Sorry, unable to allocate memory (%.1f MB) for contacts!!! Forcing exit!", prog, (float)sizeof(contact)*icont/powf(2,20));
								exit(1);
							}
						}

						// Store in RAM the computed data
						contacts[icont].ssA = 'N'; // Not computed
						contacts[icont].ssB = 'N';
						contacts[icont].seqA = seq[a];
						contacts[icont].seqB = seq[b];
						contacts[icont].froitA = fa->t;
						contacts[icont].froitB = fb->t;
						contacts[icont].sd = bond;
						contacts[icont].thetaA = ta;
						contacts[icont].thetaB = tb;
						// contacts[icont].psiA = 0.0;
						// contacts[icont].psiB = 0.0;
						contacts[icont].chi = chi;
						contacts[icont].d = dab;

						icont++; // Update contact index
					}
		}
	//	printf("%s> %d contacts found (%.1f/residue at %.0f A cutoff), %.1f%% from %d available.\n",
	//			prog,icont,(float)icont/num_res,cutoff,100*(float)icont/(num_res*(num_res-1)/2),num_res*(num_res-1)/2);
	//	printf("%s> Memory used for contacts: %ld bytes (%.1f MB)\n\n", prog, sizeof(contact)*icont, (float)sizeof(contact)*icont/powf(2,20));

	// Free memory
	free(seq);
	free(frames);

	return icont; // Successful exit
}




// Computes the 5 internal coordinates (all 6 except distance) between A and B interacting frames
//  fa,fb --> interacting frame pointers
//  ta,tb --> Theta angles
//  pa,pb --> Psi angles
//  chi   --> Chi angle
void frames2ic(frame *fa, frame *fb, float *ta, float *tb, float *pa, float *pb, float *chi)
{
	Tcoor rab,rba,v1,v2;

	// fprintf(stderr,"fa %f %f %f  fb %f %f %f\n",fa->p[0],fa->p[1],fa->p[2],fb->p[0],fb->p[1],fb->p[2]);

	// Local coordinates "a"
	diff3D(fb->p,fa->p,rab); // b-a = rab

	// Local coordinates "b"
	diff3D(fa->p,fb->p,rba); // a-b = rba

	// Angles between vectors Rab and Vz for a and b residues (Theta_a and Theta_b)
	//						ta = 90 - ang3D(frames[a].vz,rab) * 180.0 / M_PI; // "90" is required for Gnuplot's spherical mapping (set map spherical)
	//						tb = 90 - ang3D(frames[b].vz,rab) * 180.0 / M_PI; // "90" is required for Gnuplot's spherical mapping (set map spherical)
	*ta = ang3D(fa->vz,rab); // "90" is required for Gnuplot's spherical mapping (set map spherical)
	*tb = ang3D(fb->vz,rba); // "90" is required for Gnuplot's spherical mapping (set map spherical)
	// fprintf(stderr,"tb= %f\n",*tb);

	// Angles between the projection of Rab into the plane Vx-Vy and vector Vx (Psi_a and Psi_b)
	proj3D(rab,fa->vz,v1); // V1 = projection of Rab into Vz
	diff3D(rab,v1,v2); // V2 = projection of Rab into Vx-Vy plane
	*pa = ang3D(fa->vx,v2);
	if( dot3D(fa->vy,v2) < 0.0 )
		*pa *= -1.0;
	*pa += M_PI;

	proj3D(rba,fb->vz,v1); // V1 = projection of Rab into Vz
	diff3D(rba,v1,v2); // V2 = projection of Rab into Vx-Vy plane
	*pb = ang3D(fb->vx,v2);
	if( dot3D(fb->vy,v2) < 0.0 )
		*pb *= -1.0;
	*pb += M_PI;

	// Dihedral angle Chi between Vz vectors wrt Rab
	//	mult3D(fa->vz, -1); // Reverse fa->vz to obtain three "consecutive" vectors (fa->vz, urab, fb->vz)
	//	float *kk;
	//	kk=(float *) malloc(sizeof(float)*3);
	//	kk[0] = -fa->vz[0];
	//	kk[1] = -fa->vz[1];
	//	kk[2] = -fa->vz[2];

	normalize3D(rab);
	//	*chi = M_PI + dihedral3Dunit(fa->vz, rab, fb->vz);
	//	*chi = M_PI + dihedral3Dunit(kk, rab, fb->vz);

	// Using "Negative" version to reverse fa->vz to obtain three "consecutive" vectors (fa->vz, urab, fb->vz)
	// BUG: 26/4/2018
	*chi = M_PI + dihedral3DunitN(fa->vz, rab, fb->vz);

	// fprintf(stderr,"chi= %f\n",*chi);
}

// Computes the 3 internal coordinates for the Order_AVE like potentials between (i,ir) and (j,jr) interacting frames
//  i,j   --> interacting atom coordinates (previously employed in distance evaluation)
//  ir,jr --> reference atoms to define the reference frame vector
//  a,b,g --> Alpha, Beta and Gamma angles
void frames2ic_OA(Tcoor i, Tcoor j, Tcoor ir, Tcoor jr, float *a, float *b, float *g)
{
	Tcoor rab,rba,iri,jrj,iir;

	// Local coordinates "i"
	diff3D(i,j,rab); // i-j = rab
	diff3D(ir,i,iri); // ir-i = iri
	diff3D(i,ir,iir); // i-ir = iir

	// Local coordinates "j"
	diff3D(j,i,rba); // j-i = rba
	diff3D(jr,j,jrj); // jr-j = jrj

	// Angles between iri and jrj vectors wrt Rab
	//						ta = 90 - ang3D(frames[a].vz,rab) * 180.0 / M_PI; // "90" is required for Gnuplot's spherical mapping (set map spherical)
	//						tb = 90 - ang3D(frames[b].vz,rab) * 180.0 / M_PI; // "90" is required for Gnuplot's spherical mapping (set map spherical)
	*a = ang3D(iri,rba); // "90" is required for Gnuplot's spherical mapping (set map spherical)
	*b = ang3D(jrj,rab); // "90" is required for Gnuplot's spherical mapping (set map spherical)
	// fprintf(stdout,"a= %f  b= %f\n",*a,*b);

	// Dihedral angle Gamma between iri and jrj vectors wrt Rab
	normalize3D(iir);
	normalize3D(jrj);
	normalize3D(rba);
	// *g = M_PI + dihedral3Dunit(iri, rab, jrj);
	*g = M_PI - dihedral3Dunit(iir, rba, jrj); // ir-->i, i-->j, j-->jr
	// fprintf(stderr,"g= %f\n",*g);
}

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
void contact2bins(contact *c, mesh **ms, float *br, int *p_ir, int *p_ita, int *p_itb, int *p_ipa, int *p_ipb, int *p_ic)
{
	int ir,ita,itb,ipa,ipb,ic;
	mesh *m;

	ir = 0; // Radial bin
	while(c->d > br[ir+1] ) // Mon: watch out
		ir++;

	m = ms[ir];

	ita = 0; // Theta-A bin
	while(c->thetaA > m->theta[ita] && ita < m->nring-1)
		ita++;

	ipa = (int) (c->psiA / m->dpsi[ita]); // Psi-A bin
	if( ipa >= m->ncellsring[ita] ) // this only happens rarely...but it does...
		ipa = m->ncellsring[ita]-1;

	itb = 0; // Theta-B bin
	while(c->thetaB > m->theta[itb] && itb < m->nring-1)
		itb++;

	ipb = (int) (c->psiB / m->dpsi[itb]); // Psi-B bin
	if( ipb >= m->ncellsring[itb] ) // this only happens rarely...but it does...
		ipb = m->ncellsring[itb]-1;

	ic = (int) (c->chi / m->dchi);
	if( ic >= (int) roundf(2*M_PI/m->dchi) ) // this only happens rarely...but it does...
		ic = (int) roundf(2*M_PI/m->dchi)-1;

	*p_ir = ir;
	*p_ita = ita;
	*p_itb = itb;
	*p_ipa = ipa;
	*p_ipb = ipb;
	*p_ic = ic;
}

// Return the 4D bin indices from one Order_AVE contact (INT version)
//   INPUT:   c      --> Contact
//            br     --> Boundaries of shell Radius
//            nab    --> Number of Alpha or Beta bins
//            ng     --> Number of Gamma bins
//   OUTPUT:  p_ir   --> Index of Radius
//            p_ia   --> Index of Alpha
//            p_ib   --> Index of Beta
//            p_ig   --> Index of Gamma
void contact2bins_OA(contact *c, float *br, int nab, int ng, int *p_ir, int *p_ia, int *p_ib, int *p_ig)
{
	float dab = (M_PI/(float)nab); // Alpha or Beta increment

	*p_ir = 0; // Radial bin
	while(c->d > br[(*p_ir)+1] ) // Mon: watch out
		(*p_ir)++;

	*p_ia = c->thetaA / dab; // Alpha bin
	if(*p_ia >= nab)
		*p_ia = nab-1; // this should prevent overflow

	*p_ib = c->thetaB / dab; // Beta bin
	if(*p_ib >= nab)
		*p_ib = nab-1; // this should prevent overflow

	*p_ig = c->chi / (2*M_PI/(float)ng); // Gamma bin
	if(*p_ig >= ng)
		*p_ig = ng-1; // this should prevent overflow
}

// KORP 3D energy function
float korp3D(contact *contacts, int icont, korp *map)
{
	float *******maps = map->maps3D;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;
	int nsmaps = map->nsmaps;
	bool use_ji = map->use_ji;

	int ir,ita,ipa,itb,ipb,ic;
	float energy = 0.0;
	float f;

	// iaa[ aas[i] ] = i;
	int i,j,x,y,s;
	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{
			f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts
			x = contacts[c].froitA;
			y = contacts[c].froitB;
			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1+4*contacts[c].froitA]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1+4*contacts[c].froitB]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			//			fprintf(stderr,"c %d  i %d  j %d  seqA %d  seqB %d  froitA %d  froitB %d  ir %d  ita %d  itb %d  ipa %d  ipb %d  ic %d\n",
			//					c,i,j,contacts[c].seqA,contacts[c].seqB,contacts[c].froitA,contacts[c].froitB,ir,ita,itb,ipa,ipb,ic);

			energy += f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa]; // Take into account i-th (A) residue

			if(use_ji)
				energy += f * maps[s+nsmaps][x][y][i][j][ir][meshes[ir]->icell[itb]+ipb]; // Take into account j-th (B) residue

			//fprintf(stderr,"korp3D:  Eij= %8f  Eji= %8f\n",maps[i][j][ir][meshes[ir]->icell[ita]+ipa],maps[j][i][ir][meshes[ir]->icell[itb]+ipb]);

			// fprintf(stdout,"i= %c  j= %c  Eij= %10f  Eji= %10f  E= %10f\n", contacts[c].seqA, contacts[c].seqB, myicosa[i][j][ir][meshes[ir]->icell[ita]+ipa], myicosa[j][i][ir][meshes[ir]->icell[itb]+ipb], energy);
			//					fprintf(stderr,"c= %5d  E= %8.3f  E= %8.3f  d= %5.2f  ir= %d  ta= %5.1f  ita= %2d   pa= %5.1f  ipa= %2d  tb= %5.1f  itb= %2d  pb= %5.1f  ipb= %2d  ic= %d  ia= %d  ib= %d\n",
			//							c,myicosa[i][j][ir][meshes[ir]->icell[ita]+ipa],energy,contacts[c].d,ir,contacts[c].thetaA, ita, contacts[c].psiA, ipa, contacts[c].thetaB, itb, contacts[c].psiB,ipb,ic,
			//							meshes[ir]->icell[ita]+ipa, meshes[ir]->icell[itb]+ipb);
			//
			//					fprintf(stdout,"%s> ICOSA-map %c-%c  %c-%c     %f\n",prog,aas[i],aas[j],contacts[c].seqA,contacts[c].seqB,myicosa[i][j][ir][meshes[ir]->icell[ita]+ipa]);
			//					printMap(stdout, myicosa[i][j], meshes, nr, nchis, norm);
			//				exit(0);
		}
	}

	return energy;
}

// KORP 3D energy function
float korp3DM(contact *contacts, int icont, korp *map)
{
	float *******maps = map->maps3D;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;
	int nsmaps = map->nsmaps;
	bool use_ji = map->use_ji;

	int ir,ita,ipa,itb,ipb,ic;
	float energy = 0.0;
	float f;

	// iaa[ aas[i] ] = i;
	int i,j,s;
	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{
			f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts
			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			energy += f * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa]; // Take into account i-th (A) residue

			if(use_ji)
				energy += f * maps[s+nsmaps][0][0][i][j][ir][meshes[ir]->icell[itb]+ipb]; // Take into account j-th (B) residue

		}
	}

	return energy;
}


// KORP 3D energy function
float korp3DM(contact *contacts, int icont, korp *map, int posM, char chM, int Maa )
{
	float *******maps = map->maps3D;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;
	int nsmaps = map->nsmaps;
	bool use_ji = map->use_ji;

	int ir,ita,ipa,itb,ipb,ic;
	float energy = 0.0;
	float f;

	// iaa[ aas[i] ] = i;
	int i,j,s;
	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			if ( (contacts[c].fai==posM) and (contacts[c].fach==chM) )   {
				contacts[c].seqA= Maa;
			}
			if ( (contacts[c].fbi==posM) and (contacts[c].fbch==chM) )  {
				contacts[c].seqB= Maa;
			}

			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts

			energy += f * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa]; // Take into account i-th (A) residue

			if(use_ji)
				energy += f * maps[s+nsmaps][0][0][i][j][ir][meshes[ir]->icell[itb]+ipb]; // Take into account j-th (B) residue

		}
	}

	return energy;
}

// KORP 3D energy function
float korp3DM(contact *contacts, int icont, korp *map, int posM, char chM, int Maa, float **faa )
{
	float *******maps = map->maps3D;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;
	int nsmaps = map->nsmaps;
	bool use_ji = map->use_ji;

	int ir,ita,ipa,itb,ipb,ic;
	float energy = 0.0;
	float f;
	float *FAA;

	FAA = (float *) malloc( sizeof(float) * 400 );

	for(int i=0;i<20;i++) {
		for(int j=0;j<20;j++) {
			FAA[i+20*j]=0;
		}
	}
	// iaa[ aas[i] ] = i;
	int i,j,s;
	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			if ( (contacts[c].fai==posM) and (contacts[c].fach==chM) )   {
				contacts[c].seqA= Maa;
			}
			if ( (contacts[c].fbi==posM) and (contacts[c].fbch==chM) )  {
				contacts[c].seqB= Maa;
			}

			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts

			int index;
			index= (int) contacts[c].seqA+20*( (int) contacts[c].seqB);

			energy += f * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa]; // Take into account i-th (A) residue
			FAA[index] += f * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa]; // Take into account i-th (A) residue

			if(use_ji) {
				energy += f * maps[s+nsmaps][0][0][i][j][ir][meshes[ir]->icell[itb]+ipb]; // Take into account j-th (B) residue
				FAA[index] += f * maps[s+nsmaps][0][0][i][j][ir][meshes[ir]->icell[itb]+ipb]; // Take into account j-th (B) residue

			}


		}
	}
	*faa=FAA;
	return energy;
}


// KORP 3D energy function for mutations
float korp3DMW(contact *contacts, int icont, korp *map, int posM, char chM, int Maa, double *W )
{

	float *******maps = map->maps3D;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;
	int nsmaps = map->nsmaps;
	bool use_ji = map->use_ji;

	int ir,ita,ipa,itb,ipb,ic;
	float energy = 0.0;
	float f;

	// iaa[ aas[i] ] = i;
	int i,j,s;
	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			if ( (contacts[c].fai==posM) and (contacts[c].fach==chM) )   {
				contacts[c].seqA= Maa;
			}
			if ( (contacts[c].fbi==posM) and (contacts[c].fbch==chM) )  {
				contacts[c].seqB= Maa;
			}

			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			f = fmapping[s]* W[(int)contacts[c].seqA]*W[(int)contacts[c].seqB]; // Custom weighting factors for Bonding and Non-bonding contacts

			energy += f * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa]; // Take into account i-th (A) residue

			if(use_ji)
				energy += f * maps[s+nsmaps][0][0][i][j][ir][meshes[ir]->icell[itb]+ipb]; // Take into account j-th (B) residue

		}
	}

	return energy/100.0;
}



// KORP 4D energy function
float korp4D(contact *contacts, int icont, korp *map)
{
	float *********maps = map->maps;
	int nab = map->ncells;
	int ng = map->nchis;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	float energy = 0.0;
	int ir,ita,itb,ic,x,y,i,j,s;
	float f;

	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{
			f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts

			x = contacts[c].froitA;
			y = contacts[c].froitB;

			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1+4*x]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1+4*y]; // j-th interaction frame map

			contact2bins_OA(&contacts[c],br,nab,ng,&ir,&ita,&itb,&ic); // Return the 4D bin indices from one contact

			energy += f * maps[s][x][y][i][j][ir][ita][itb][ic];
		}
	}
	return energy;
}

// KORP 6D energy function
double korp6D(contact *contacts, int icont, korp *map)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,x,y,i,j,s;

	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{
			x = contacts[c].froitA;
			y = contacts[c].froitB;

			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1+4*contacts[c].froitA]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1+4*contacts[c].froitB]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			//			fprintf(stdout,"contact2bins> ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",ir,ita,itb,ipa,ipb,ic,
			//					maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

			// f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts
			energy += fmapping[s] * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic];

			//			fprintf(stderr,"%4d %4d %4d %4d %4d %2d %2d %8.3f  ita %4d ipa %4d itb %4d ir %4d  E %8.5f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
			//					          s,  contacts[c].froitA,  contacts[c].froitB,  i, j, contacts[c].seqA, contacts[c].seqB, contacts[c].d,      ita,    ipa,    itb, ir,
			//							  f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic],
			//							  contacts[c].thetaA, contacts[c].thetaB, contacts[c].psiA, contacts[c].psiB,contacts[c].chi);
		}
	}
	return energy;
}

// KORP 6D energy function
double korp6DW(contact *contacts, int icont, korp *map, double  *W)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,x,y,i,j,s;
	float f;

	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{
			x = contacts[c].froitA;
			y = contacts[c].froitB;

			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1+4*contacts[c].froitA]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1+4*contacts[c].froitB]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			//			fprintf(stdout,"contact2bins> ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",ir,ita,itb,ipa,ipb,ic,
			//					maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

			f = fmapping[s]* W[(int)contacts[c].seqA]*W[(int)contacts[c].seqB]; // Custom weighting factors for Bonding and Non-bonding contacts

			energy += f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic];

			//			fprintf(stderr,"%4d %4d %4d %4d %4d %2d %2d %8.3f  ita %4d ipa %4d itb %4d ir %4d  E %8.5f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
			//					          s,  contacts[c].froitA,  contacts[c].froitB,  i, j, contacts[c].seqA, contacts[c].seqB, contacts[c].d,      ita,    ipa,    itb, ir,
			//							  f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic],
			//							  contacts[c].thetaA, contacts[c].thetaB, contacts[c].psiA, contacts[c].psiB,contacts[c].chi);
		}
	}
	return energy;
}


// KORP 6D energy function removed froitAB
double korp6DM(contact *contacts, int icont, korp *map)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	float f;

	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map


			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			//			fprintf(stdout,"contact2bins> ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",ir,ita,itb,ipa,ipb,ic,
			//					maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

			f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts
			energy += fmapping[s] * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic];

			//			fprintf(stderr,"%4d %4d %4d %4d %4d %2d %2d %8.3f  ita %4d ipa %4d itb %4d ir %4d  E %8.5f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
			//					          s,  contacts[c].froitA,  contacts[c].froitB,  i, j, contacts[c].seqA, contacts[c].seqB, contacts[c].d,      ita,    ipa,    itb, ir,
			//							  f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic],
			//							  contacts[c].thetaA, contacts[c].thetaB, contacts[c].psiA, contacts[c].psiB,contacts[c].chi);

		}
	}
	return energy;
}

// MUT KORP 6D energy function
double korp6DM(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;


	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	float f;

	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			// printf("%d fai %d %d pos %d Maa %c cont %c\n", c, contacts[c].fai,contacts[c].fbi, posIN, Maa, contacts[c].seqA);
			if ( (contacts[c].fai==posM) and (contacts[c].fach==chM) )   {
				contacts[c].seqA= Maa;
			}
			if ( (contacts[c].fbi==posM) and (contacts[c].fbch==chM) )  {
				contacts[c].seqB= Maa;
			}
			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map



			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			//			fprintf(stdout,"contact2bins> ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",ir,ita,itb,ipa,ipb,ic,
			//					maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

			f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts
			energy += fmapping[s] * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic];

			//			fprintf(stderr,"%4d %4d %4d %4d %4d %2d %2d %8.3f  ita %4d ipa %4d itb %4d ir %4d  E %8.5f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
			//					          s,  contacts[c].froitA,  contacts[c].froitB,  i, j, contacts[c].seqA, contacts[c].seqB, contacts[c].d,      ita,    ipa,    itb, ir,
			//							  f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic],
			//							  contacts[c].thetaA, contacts[c].thetaB, contacts[c].psiA, contacts[c].psiB,contacts[c].chi);

		}
	}
	return energy;
}


// MUT KORP 6D energy function with partial contributions faa....
double korp6DM(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa, float *faa)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	int c_seqA, c_seqB;

	//	float *FAA;
	//	FAA = (float *) malloc( sizeof(float) * 400 );

	for(int i=0;i<400;i++)
		//    	FAA[i] = 0.0;
		faa[i] = 0.0;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	float f;

	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
	    if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			// printf("%d fai %d %d pos %d Maa %c cont %c\n", c, contacts[c].fai,contacts[c].fbi, posIN, Maa, contacts[c].seqA);


			if ( (contacts[c].fai == posM) and (contacts[c].fach == chM) )
			{
				c_seqA = Maa;
			} else c_seqA = contacts[c].seqA;

			if ( (contacts[c].fbi == posM) and (contacts[c].fbch == chM) )
			{
				c_seqB = Maa;
			} else c_seqB = contacts[c].seqB;

			int index;

			index= (int) c_seqA+20*( (int) c_seqB);

			/*
			if (c_seqA==19) c_seqA=18;
		    if (c_seqB==19) c_seqB=18;

		    if (c_seqA==4) c_seqA=18;
		    if (c_seqB==4) c_seqB=18;
            */

			i = mapping[ (int)iaa[c_seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[c_seqB] ] [1]; // j-th interaction frame map



		    contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			//			fprintf(stdout,"contact2bins> ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",ir,ita,itb,ipa,ipb,ic,
			//					maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

			f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts
			energy += fmapping[s] * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic];

			//			fprintf(stderr,"%4d %4d %4d %4d %4d %2d %2d %8.3f  ita %4d ipa %4d itb %4d ir %4d  E %8.5f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
			//					          s,  contacts[c].froitA,  contacts[c].froitB,  i, j, contacts[c].seqA, contacts[c].seqB, contacts[c].d,      ita,    ipa,    itb, ir,
			//							  f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic],
			//							  contacts[c].thetaA, contacts[c].thetaB, contacts[c].psiA, contacts[c].psiB,contacts[c].chi);
			//printf (" %d ", index);

			//			 FAA[index] += fmapping[s] * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic];
			faa[index] += fmapping[s] * maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic];
		}
	}

	//	*faa = FAA;
	return energy;
}

// MUT KORP 6D energy function with partial contributions faa....
double korp6DM_BIND(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa, float *faa)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	//	float *FAA;
	//	FAA = (float *) malloc( sizeof(float) * 400 );

	for(int i=0;i<400;i++)
		//    	FAA[i] = 0.0;
		faa[i] = 0.0;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	float f;

	//	int kk = 0;

	for(long c=0; c<icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			// printf("%d fai %d %d pos %d Maa %c cont %c\n", c, contacts[c].fai,contacts[c].fbi, posIN, Maa, contacts[c].seqA);
			if( contacts[c].fai == posM and contacts[c].fach == chM )
				contacts[c].seqA = Maa;

			if( contacts[c].fbi == posM and contacts[c].fbch == chM )
				contacts[c].seqB = Maa;

			if( !(contacts[c].fach==chM and contacts[c].fbch==chM) )
			{
				i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
				j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map

				int index;
				index = (int) contacts[c].seqA + 20 * ( (int) contacts[c].seqB );

				contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

				//			fprintf(stdout,"contact2bins> ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",ir,ita,itb,ipa,ipb,ic,
				//					maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

				f = fmapping[s]; // Custom weighting factors for Bonding and Non-bonding contacts
				float value=maps[s][0][0][i][j][ir][ meshes[ir]->icell[ita] + ipa ][ meshes[ir]->icell[itb] + ipb ][ic];
				// if (value > 3 ) value=3;

				energy += fmapping[s] * value;
				//kk++;

				//			fprintf(stderr,"%4d %4d %4d %4d %4d %2d %2d %8.3f  ita %4d ipa %4d itb %4d ir %4d  E %8.5f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
				//					          s,  contacts[c].froitA,  contacts[c].froitB,  i, j, contacts[c].seqA, contacts[c].seqB, contacts[c].d,      ita,    ipa,    itb, ir,
				//							  f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic],
				//							  contacts[c].thetaA, contacts[c].thetaB, contacts[c].psiA, contacts[c].psiB,contacts[c].chi);
				//printf (" %d ", index);

				faa[index] += fmapping[s] * value;
			}
		}
	}
	//fprintf(stderr,"icont= %d  --> kk= %d\n",icont, kk);
	//	*faa = FAA;
	return energy;
}

// MUT WEIGHTED  KORP 6D energy function
double korp6DMW(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa, double *W)
{
	// double W[20]={ 0.023754, 0.022738, 0.016045, 0.021554, 0.038330, 0.022697, 0.001958, 0.039904, 0.020681, 0.049354, 0.027223, 0.026858, 0.024015, 0.012121, 0.034633, 0.022557, 0.025499, 0.038887, 0.022436, 0.039417,};
	// double W[20]={ 2.301813, 2.262442, 1.575543, 2.141481, 3.790123, 2.194824, 0.205695, 3.913127, 2.050392, 4.761978, 2.681753, 2.698784, 2.370083, 1.206800, 3.442918, 2.212829, 2.519377, 3.794067, 2.234677, 3.926169};
	// double W[20]={ 1.592,  1.613,  1.320,  1.750,  1.997,  1.552,  1.323,  2.453,  1.790,  2.678,  1.621,  1.846,  1.578,  1.201,  2.084,  1.748,  2.028,  2.295,  2.085,  2.212};
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	int c_seqA, c_seqB;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	double f;

	for(long c = 0; c < icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			//printf("%d fai %d %d Maa %c cont %c  %f %f\n", c, contacts[c].fai,contacts[c].fbi, posM, Maa, contacts[c].seqA, contacts[c].d, br[0]);
			if ( (contacts[c].fai == posM) and (contacts[c].fach == chM) )
			{
				c_seqA = Maa;
			} else c_seqA = contacts[c].seqA;

			if ( (contacts[c].fbi == posM) and (contacts[c].fbch == chM) )
			{
				c_seqB = Maa;
			} else c_seqB = contacts[c].seqB;

			i = mapping[ (int)iaa[c_seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[c_seqB] ] [1]; // j-th interaction frame map


			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

		//	fprintf(stdout,"contact2bins> %2d-%2d ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",c_seqA, c_seqB, ir,ita,itb,ipa,ipb,ic,
		//						maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

			f = fmapping[s] * W[c_seqA] * W[c_seqB]; // Custom weighting factors for Bonding and Non-bonding contacts
			float value = maps[s][0][0][i][j][ir][ meshes[ir]->icell[ita] + ipa ][ meshes[ir]->icell[itb] + ipb ][ic];
			// if (value > 3 ) value=3;

			energy += f * value;



		}
	}


	return energy/100; // MON: What is this? It may be not said at Methods...
}


// Update contacts for Multiple Point Mutations (MPMs)
int contact_MPM(contact *contacts, int icont, int *posM, char *chM,  int *Maa, int nmut)
{
	int nup = 0; // Number of Updated contacts

	// Update Mutant Residues Id in contact list
	for(long c = 0; c < icont; c++) // screen all contacts
		for(int m = 0; m < nmut; m++)
		{
			if ( (contacts[c].fai == posM[m]) and (contacts[c].fach == chM[m]) )
			{
				contacts[c].seqA = Maa[m];
				nup++;
			}

			if ( (contacts[c].fbi == posM[m]) and (contacts[c].fbch == chM[m]) )
			{
				contacts[c].seqB = Maa[m];
				nup++;
			}
		}

	return(nup);
}

// KORP 6D energy function for Multiple Point Mutations (MPMs) including weights (W)
// (Set "nmut" to 0 for WT)
double korp6DMW_MPM(contact *contacts, int icont, korp *map,  int *posM, char *chM,  int *Maa, int nmut, double *W)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	int c_seqA, c_seqB;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	double f;

	if(nmut<1)
	{
		fprintf(stderr,"Error, invalid number of mutations: nmut= %d\n", nmut);
		exit(1);
	}

	for(long c = 0; c < icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{
			//printf("%d fai %d %d Maa %c cont %c  %f %f\n", c, contacts[c].fai,contacts[c].fbi, posM, Maa, contacts[c].seqA, contacts[c].d, br[0]);
			c_seqA = contacts[c].seqA;
			c_seqB = contacts[c].seqB;

			for(int m = 0; m < nmut; m++)
			{
				if ( (contacts[c].fai == posM[m]) and (contacts[c].fach == chM[m]) )
					c_seqA = Maa[m];
				else if ( (contacts[c].fbi == posM[m]) and (contacts[c].fbch == chM[m]) )
					c_seqB = Maa[m];
			}

			i = mapping[ (int)iaa[c_seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[c_seqB] ] [1]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

		//	fprintf(stdout,"contact2bins> %2d-%2d ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",c_seqA, c_seqB, ir,ita,itb,ipa,ipb,ic,
		//		maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

			f = fmapping[s] * W[c_seqA] * W[c_seqB]; // Custom weighting factors for Bonding and Non-bonding contacts
			float value = maps[s][0][0][i][j][ir][ meshes[ir]->icell[ita] + ipa ][ meshes[ir]->icell[itb] + ipb ][ic];

			energy += f * value;

		//	fprintf(stderr, "korp6DMW_MPM>  i= %3d  j= %3d  c_seqA= %d  c_seqB= %d  f= %10.6f  energy= %10.6f\n", i, j, c_seqA, c_seqB, f, energy);
		}
	}

	return energy/100; // MON: What is this? It may be not said at Methods...
}


// KORP 6D energy function including weights (W)
double korp6DW2(contact *contacts, int icont, korp *map, double *W)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	int c_seqA, c_seqB;

	double energy = 0.0;
	int ir, ita, ipa, itb, ipb, ic, i, j, s;
	double f;

	for(long c = 0; c < icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{
			c_seqA = contacts[c].seqA;
			c_seqB = contacts[c].seqB;

			i = mapping[ (int)iaa[c_seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[c_seqB] ] [1]; // j-th interaction frame map

			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

		//	fprintf(stdout,"contact2bins> %2d-%2d ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",c_seqA, c_seqB, ir,ita,itb,ipa,ipb,ic,
		//		maps[s][0][0][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

			f = fmapping[s] * W[c_seqA] * W[c_seqB]; // Custom weighting factors for Bonding and Non-bonding contacts
			float value = maps[s][0][0][i][j][ir][ meshes[ir]->icell[ita] + ipa ][ meshes[ir]->icell[itb] + ipb ][ic];
			energy += f * value;

			// fprintf(stderr, "korp6DMW_MPM>  i= %3d  j= %3d  c_seqA= %d  c_seqB= %d  f= %10.6f  energy= %10.6f\n", i, j, c_seqA, c_seqB, f, energy);
		}
	}

	return energy/100; // MON: What is this? It may be not said at Methods...
}

// MUT with 21 WEIGHTs KORP 6D energy function (21st is an additive term)
double korp6DMW21(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa, double *W)
{
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	double f;
	int c_seqA, c_seqB;

	for(long c = 0; c < icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like

		{


			if ( (contacts[c].fai == posM) and (contacts[c].fach == chM) )
			{
				c_seqA = Maa;
			} else c_seqA = contacts[c].seqA;

			if ( (contacts[c].fbi == posM) and (contacts[c].fbch == chM) )
			{
				c_seqB = Maa;
			} else c_seqB = contacts[c].seqB;

			i = mapping[ (int)iaa[c_seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[c_seqB] ] [1]; // j-th interaction frame map


			// fprintf(stdout,"\n dist A %2d  B %2d  D %f %f\n", c_seqA, c_seqB, contacts[c].d, br[0] );

//			if (c_seqA>=18 || c_seqB>=18 || c_seqA==6 || c_seqB ==6 && contacts[c].d < 3.5  ) {
//				energy += 0;
//			}
//			else {
			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			f = fmapping[s] * W[c_seqA] * W[c_seqB] * maps[s][0][0][i][j][ir][ meshes[ir]->icell[ita] + ipa ][ meshes[ir]->icell[itb] + ipb ][ic]; // Custom weighting factors for Bonding and Non-bonding contacts
	        // energy +=  f + (W[20] + exp(-1.0*f)/(1.0+pow(exp(-1.0*f),2)));
			energy +=  f + W[20];
//			}

		//	f = fmapping[s] * W[contacts[c].seqA] * W[contacts[c].seqB]; // Custom weighting factors for Bonding and Non-bonding contacts
		//	energy += W[20] + f * maps[s][0][0][i][j][ir][ meshes[ir]->icell[ita] + ipa ][ meshes[ir]->icell[itb] + ipb ][ic];
		}
	}

	return energy/100;
}


// MUT WEIGHTED  KORP 6D energy function
double korp6DMWRSA(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa, double *W, float RSA)
{
	// double W[20]={ 0.023754, 0.022738, 0.016045, 0.021554, 0.038330, 0.022697, 0.001958, 0.039904, 0.020681, 0.049354, 0.027223, 0.026858, 0.024015, 0.012121, 0.034633, 0.022557, 0.025499, 0.038887, 0.022436, 0.039417,};
	// double W[20]={ 2.301813, 2.262442, 1.575543, 2.141481, 3.790123, 2.194824, 0.205695, 3.913127, 2.050392, 4.761978, 2.681753, 2.698784, 2.370083, 1.206800, 3.442918, 2.212829, 2.519377, 3.794067, 2.234677, 3.926169};
	// double W[20]={ 1.592,  1.613,  1.320,  1.750,  1.997,  1.552,  1.323,  2.453,  1.790,  2.678,  1.621,  1.846,  1.578,  1.201,  2.084,  1.748,  2.028,  2.295,  2.085,  2.212};
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	double f;

	for(long c = 0; c < icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			// printf("%d fai %d %d pos %d Maa %c cont %c\n", c, contacts[c].fai,contacts[c].fbi, posIN, Maa, contacts[c].seqA);
			if ( (contacts[c].fai == posM) and (contacts[c].fach == chM) )
			{
				//fprintf(stderr,"seqA= %2d  Maa= %2d  posM= %4d  chM= %c  chB= %c\n", contacts[c].seqA, Maa, posM, chM, contacts[c].fbch);
				contacts[c].seqA = Maa;
			}
			if ( (contacts[c].fbi == posM) and (contacts[c].fbch == chM) )
			{
				//fprintf(stderr,"seqB= %2d  Maa= %2d  posM= %4d  chM= %c  chA= %c\n", contacts[c].seqB, Maa, posM, chM, contacts[c].fach);
				contacts[c].seqB= Maa;
			}
			i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
			j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map


			contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

			//			fprintf(stdout,"contact2bins> ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",ir,ita,itb,ipa,ipb,ic,
			//					maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);
			if (RSA > 20.0)
				f = fmapping[s] * W[(int)contacts[c].seqA] * W[(int)contacts[c].seqB]; // Custom weighting factors for Bonding and Non-bonding contacts
			else
				f = fmapping[s] * W[(int)contacts[c].seqA+20] * W[(int)contacts[c].seqB+20];

			energy += f * maps[s][0][0][i][j][ir][ meshes[ir]->icell[ita] + ipa ][ meshes[ir]->icell[itb] + ipb ][ic];

			//			fprintf(stderr,"%4d %4d %4d %4d %4d %2d %2d %8.3f  ita %4d ipa %4d itb %4d ir %4d  E %8.5f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
			//					          s,  contacts[c].froitA,  contacts[c].froitB,  i, j, contacts[c].seqA, contacts[c].seqB, contacts[c].d,      ita,    ipa,    itb, ir,
			//							  f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic],
			//							  contacts[c].thetaA, contacts[c].thetaB, contacts[c].psiA, contacts[c].psiB,contacts[c].chi);
			//printf (" %d ", index);

		}
	}


	return energy/100;
}



double korp6DMW_BIND(contact *contacts, int icont, korp *map,  int posM, char chM,  int Maa, double *W)
{
	// double W[20]={ 0.023754, 0.022738, 0.016045, 0.021554, 0.038330, 0.022697, 0.001958, 0.039904, 0.020681, 0.049354, 0.027223, 0.026858, 0.024015, 0.012121, 0.034633, 0.022557, 0.025499, 0.038887, 0.022436, 0.039417,};
	// double W[20]={ 2.301813, 2.262442, 1.575543, 2.141481, 3.790123, 2.194824, 0.205695, 3.913127, 2.050392, 4.761978, 2.681753, 2.698784, 2.370083, 1.206800, 3.442918, 2.212829, 2.519377, 3.794067, 2.234677, 3.926169};
	// double W[20]={ 1.592,  1.613,  1.320,  1.750,  1.997,  1.552,  1.323,  2.453,  1.790,  2.678,  1.621,  1.846,  1.578,  1.201,  2.084,  1.748,  2.028,  2.295,  2.085,  2.212};
	float *********maps = map->maps;
	mesh **meshes = map->meshes;
	char **mapping = map->mapping;
	char *smapping = map->smapping;
	float *fmapping = map->fmapping;
	float *br = map->br;
	char *iaa = map->iaa;

	double energy = 0.0;
	int ir,ita,ipa,itb,ipb,ic,i,j,s;
	double f;

	for(long c = 0; c < icont; c++) // screen all contacts
	{
		s = smapping[(int)contacts[c].sd]; // 0= Non-bonding, >0= Bonding: 1= i+x, 2= i+x+1, 3= i+x+2, etc...
		if( s >= 0 && contacts[c].d > br[0]  ) // ICOSA-like
		{

			// printf("%d fai %d %d pos %d Maa %c cont %c\n", c, contacts[c].fai,contacts[c].fbi, posIN, Maa, contacts[c].seqA);
			if ( contacts[c].fai == posM and contacts[c].fach == chM )
			{
				//fprintf(stderr,"seqA= %2d  Maa= %2d  posM= %4d  chM= %c  chB= %c\n", contacts[c].seqA, Maa, posM, chM, contacts[c].fbch);
				contacts[c].seqA = Maa;
			}

			if ( contacts[c].fbi == posM and contacts[c].fbch == chM )
			{
				//fprintf(stderr,"seqB= %2d  Maa= %2d  posM= %4d  chM= %c  chA= %c\n", contacts[c].seqB, Maa, posM, chM, contacts[c].fach);
				contacts[c].seqB= Maa;
			}

			// Only inter-chain contacts is the same as "Binding" experiment
			if ( !(contacts[c].fach==chM and contacts[c].fbch==chM) )
			{
				i = mapping[ (int)iaa[(int)contacts[c].seqA] ] [1]; // i-th interaction frame map
				j = mapping[ (int)iaa[(int)contacts[c].seqB] ] [1]; // j-th interaction frame map

				contact2bins(&contacts[c],meshes,br,&ir,&ita,&itb,&ipa,&ipb,&ic); // Return the 6D bin indices from one contact

				//			fprintf(stdout,"contact2bins> ir= %2d  ita= %2d  itb= %2d  ipa= %2d  ipb= %2d  ic= %2d  energy= %f\n",ir,ita,itb,ipa,ipb,ic,
				//					maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic]);

				f = fmapping[s] * W[(int)contacts[c].seqA] * W[(int)contacts[c].seqB]; // Custom weighting factors for Bonding and Non-bonding contacts

				energy += f * maps[s][0][0][i][j][ir][ meshes[ir]->icell[ita] + ipa ][ meshes[ir]->icell[itb] + ipb ][ic];

				//			fprintf(stderr,"%4d %4d %4d %4d %4d %2d %2d %8.3f  ita %4d ipa %4d itb %4d ir %4d  E %8.5f  ta= %8.6f  tb= %8.6f  pa= %8.6f  pb= %8.6f  chi= %8.6f\n",
				//					          s,  contacts[c].froitA,  contacts[c].froitB,  i, j, contacts[c].seqA, contacts[c].seqB, contacts[c].d,      ita,    ipa,    itb, ir,
				//							  f * maps[s][x][y][i][j][ir][meshes[ir]->icell[ita]+ipa][meshes[ir]->icell[itb]+ipb][ic],
				//							  contacts[c].thetaA, contacts[c].thetaB, contacts[c].psiA, contacts[c].psiB,contacts[c].chi);
				//printf (" %d ", index);
			}
		}
	}


	return energy/100;
}


// Get a per residue array of PDB's residue numbers (automatic memory allocation)
int *getResNums(Macromolecule *mol)
{
	int *resnums = (int *) malloc( sizeof(int) * mol->get_num_fragments() );
	pdbIter *iter_res = new pdbIter( mol );

	for ( iter_res->pos_fragment = 0; !iter_res->gend_fragment(); iter_res->next_fragment() ) // iter frags
	{resnums[iter_res->pos_fragment] = (iter_res->get_fragment())->getIdNumber(); // Get PDB numeration (required for correct bonding/non-bonding discrimination)
	// printf ("%d %d\n", iter_res->pos_fragment, resnums[iter_res->pos_fragment]);
	}
	delete iter_res;
	return resnums;
}

// Get a per residue array of PDB's chain-ids (automatic memory allocation)
char *getResChainIds(Macromolecule *mol)
{
	char *reschainids = (char *) malloc( sizeof(char) * mol->get_num_fragments() );
	pdbIter *iter_chain = new pdbIter( mol );
	pdbIter *iter_res;
	Chain *ch;
	char ichain; // current chain
	int ires = 0; // global residue index
	for ( iter_chain->pos_chain = 0; !iter_chain->gend_chain(); iter_chain->next_chain() ) // Iter chains to compute sequential distances fine
	{
		ch = (Chain *) iter_chain->get_chain();
		iter_res = new pdbIter( ch );
		ichain = (ch->getName())[0]; // Get current chain id (1-letter)

		for ( iter_res->pos_fragment = 0; !iter_res->gend_fragment(); iter_res->next_fragment() ) // iter frags
		{
			reschainids[ires] = ichain; // Store chain id
			ires++; // update residue index
		}
		delete iter_res;
	}
	delete iter_chain;

	return reschainids;
}



// Prints into "f" file-handle all contact information in the "ncont" contacts
void print_contacts(FILE *f, contact *contacts, long ncont, bool header)
{
	fprintf(f,"contacts_inside= %p (Hex)  ncont= %ld\n",contacts,ncont);
	// PABLO REMOVE WARNING 2023
	// fprintf(f,"contacts_inside= %x (Hex)  ncont= %ld\n",contacts,ncont);


	if(header)
		fprintf(f,"#%2s %5s %3s %2s %3s  %6s %5s %5s %5s %5s %5s\n",
				"SS","ResID","FrT","In","A-B","Dab[A]","ThA","ThB","PsiA","PsiB","Chi");

	contact *c;
	for(long i=0; i<ncont; i++)
	{
		c = &(contacts[i]);
		fprintf(f," %c%c %2d %2d %1d %1d %2d %3d  %6.3f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
				c->ssA,c->ssB,c->seqA,c->seqB,c->froitA,c->froitB,c->intra,c->sd,
				c->d,c->thetaA,c->thetaB,c->psiA,c->psiB,c->chi);
	}
}
