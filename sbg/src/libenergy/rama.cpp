/*
 * rama.cpp
 *
 *  Created on: Jan 4, 2019
 *      Author: mon
 */

#include "libenergy/include/rama.h" // Ramachandran energies

// Reads all Dunbrack's PDFs from file (created with "rang" program)
float *****read_dunbrack(char *file, int *p_size, float *p_step)
{
	FILE *f_map; // File handler
	if( !(f_map = fopen(file,"rb")) ) // Read binary file
	{
		printf("Reading Dunbrack's PDFs error file \"%s\" missing\n",file);
		exit(1);
	}

	// Write header with "size" and "step"
	fread(p_size,sizeof(int),1,f_map);
	fread(p_step,sizeof(float),1,f_map);

	// Memory allocation
	float *****pdfs; // 5*pointer containing all pdfs
	int mem=0;
	pdfs = (float *****) malloc( sizeof(float ****) * NAAS); // Central aminoacid level
	mem += sizeof(float ****) * NAAS;
	for( int c=0; c<NAAS; c++)
	{
		pdfs[c] = (float ****) malloc( sizeof(float ***) * 2); // Whether it is a Left-side (0) or Right-side (1) aminoacid
		mem += sizeof(float ***) * 2;
		for( int lr=0; lr<2; lr++)
		{
			pdfs[c][lr] = (float ***) malloc( sizeof(float **) * (NAAS+1));  // Left/Right-side aminoacid (20aas + 1 cis-PRO + 1 "ALL")
			mem += sizeof(float **) * NAAS+1;
			for( int s=0; s<NAAS+1; s++)
			{
				pdfs[c][lr][s] = (float **) malloc( sizeof(float *) * (*p_size)); // Pointer array
				mem += sizeof(float *) * (*p_size);
				for( int i=0; i<*p_size; i++)
				{
					pdfs[c][lr][s][i] = (float *) malloc( sizeof(float) * (*p_size)); // data
					mem += sizeof(float) * (*p_size);
					for( int j=0; j<*p_size; j++)
						pdfs[c][lr][s][i][j] = 0.0; // initialize pdf
				}
			}
		}
	}
	// fprintf(stdout,"Msg(read_dunbrack): Allocated %d bytes!\n",mem);

	for( int c=0; c<NAAS; c++) // 20=cis-Pro, 21=ALL
		for( int lr=0; lr<2; lr++)
			for( int s=0; s<NAAS+1; s++)
				for( int i=0; i<*p_size; i++)
					fread(pdfs[c][lr][s][i], sizeof(float), *p_size, f_map); // read data  ("i" is the Phi angle index)

	fclose(f_map);
	return pdfs;
}


//  Reads all Dunbrack's H3 PDFs from file
float ***read_h3(char *file, int *p_size, float *p_step)
{

	int ntot, nnom, nnt, nct, bins;
	float paso;

	float ***mapas;

	FILE *f; // File handler
	if( !(f = fopen(file,"rb")) ) // Read binary file
	{
		printf("Reading H3 map input error! file %s\n",file);
		exit(1);
	}



	fread(&ntot, sizeof(int), 1, f);
	fread(&nnom, sizeof(int), 1, f);
	fread(&nnt, sizeof(int), 1, f);
	fread(&nct, sizeof(int), 1, f);

	printf("rcd> reading %d maps aa %d nt %d ct %d\n", ntot, nnom, nnt, nct);

	mapas = (float ***)malloc(sizeof(float **)*ntot);

	char mapaname[ntot][20];
	for (int i = 0; i < ntot; ++i)
	{
		fread(&bins, sizeof(int), 1, f);  // must be a int....
		fread(&paso, sizeof(float), 1, f);
		fread(&mapaname[i], sizeof(char), 20, f);
		mapas[i] = (float **)malloc(sizeof(float *)*bins);
		// printf("# maps %s %d %f\n", mapaname[i], bins, paso );

		for(int j=0; j < bins; ++j)
		{
			mapas[i][j] = (float *)malloc(sizeof(float)*bins);
			for(int k=0; k < bins; ++k)
			{
				fread(&mapas[i][j][k], sizeof(float), 1, f);

			}
		}
	}

	*p_size=bins;
	*p_step=paso;
	fclose(f);

	printf("rcd> ");
	float contador[ntot];
	for (int i = 0; i < ntot; ++i )
	{
		contador[i] = 0;
		for (int j = 0; j < bins; ++j )
		{
			for (int k = 0; k < bins; ++k )
			{
				contador[i] = contador[i] + mapas[i][j][k];
			}
		}
		printf("%2d-%s(%5.3f) ", i, mapaname[i], contador[i]);
	}

	printf("\n");



	return mapas;
}




// Returns the sequence indices array from the 1-letter code sequence (seq) for our Rama potential (automatic memory allocation)
//  seq --> Sequence in 1-letter code
//  nresl --> Number of residues
int *seq2iaa(char *seq, int nres)
{
	int *iaa;

	//    Dunbrack based PDFs
	if( !(iaa = (int *) malloc( sizeof(int) * nres ) ) )
	{
		printf("seq2iaa> Error: Sorry, unable to allocate iaa memory!!!\n");
		exit(1);
	}

	// Loop AA ID (including anchors)
	for (int i = 0; i < nres; i++)
	{
		if(seq[i]=='p') // If cis-Proline
		{
			printf("seq2iaa> Remark: cis-Proline requested by user in residue %d\n",i+1);
			iaa[i] = 20; // cis-Proline index
		}
		else // other aminoacids
			for (int j = 0; j < 20;j++)
				if(seq[i] == AA[j].aa_name1)
				{
					iaa[i]= j; // Array with ID number of AA except C-anchor +1 or N-anchor -1 (added previously in line 952)
					break;
				}
	}

	return iaa;
}

// Allocate and zero initializes a size x size matrix
float **MatrixInit(int size)
{
	float **prod;
	prod = (float **) malloc( sizeof(float *) * size);
	for( int i=0; i<size; i++)
	{
		prod[i] = (float *) malloc( sizeof(float) * size);
		for( int j=0; j<size; j++)
			prod[i][j] = 0.0; // initialize
	}
	return prod;
}

// Element-wise accumulation of M2 (weighted or not) into M1 (size x size squared matrices)
// Optionally, M2 can be weighted by some factor "w"
void MatrixAccum(float **M1, float **M2, int size, float w)
{
	for( int i=0; i<size; i++)
		for( int j=0; j<size; j++)
			M1[i][j] += w * M2[i][j];
}

// Element-wise multiplication of two squared matrices (size x size)
float **MatrixProduct(float **M1, float **M2, int size)
{
	float **prod;
	if( !(prod = (float **) malloc( sizeof(float *) * size)) )
	{
		printf("Sorry, unable to allocate prod memory!!!\n");
		exit(1);
	}

	for( int i=0; i<size; i++)
	{
		if( !(prod[i] = (float *) malloc( sizeof(float) * size)) )
		{
			printf("Sorry, unable to allocate prod[i] memory!!!\n");
			exit(1);
		}

		for( int j=0; j<size; j++)
			prod[i][j] = 0.0; // initialize
	}

	for( int i=0; i<size; i++)
		for( int j=0; j<size; j++)
			prod[i][j] = M1[i][j]*M2[i][j];

	return prod;
}

// Element-wise division of two squared matrices (size x size)
float **MatrixDivision(float **M1, float **M2, int size)
{
	float **divgc;
	if( !(divgc = (float **) malloc( sizeof(float *) * size)) )
	{
		printf("Sorry, unable to allocate divgc memory!!!\n");
		exit(1);
	}

	for( int i=0; i<size; i++)
	{
		if( !(divgc[i] = (float *) malloc( sizeof(float) * size)) )
		{
			printf("Sorry, unable to allocate divgc[i] memory!!!\n");
			exit(1);
		}
		for( int j=0; j<size; j++)
			divgc[i][j] = 0.0; // initialize
	}

	for( int i=0; i<size; i++)
		for( int j=0; j<size; j++)
			divgc[i][j] = M1[i][j]/M2[i][j];

	return divgc;
}

// Normalizes "in situ" a map (size x size) to sum=1.0. It does not allocate memory.
void norm_map(float **map, int size)
{
	double acum=0.0;

	// Get sum
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			acum += map[i][j];

	// Normalize
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			map[i][j] /= acum;
}

// Writes a pdf map in some file (size x size)
void show_map(float **map,int size, char *file)
{
	FILE *f_file;
	if( !(f_file = fopen(file,"w")) )
	{
		printf("Reading input error! file %s\n",file);
		exit(1);
	}

	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			fprintf(f_file, "%5d %5d %12e \n", i, j, map[i][j]);

	fclose(f_file);
}


// Generate the conditional PDF that a Central (C) residue has Left (L) and Right (R) neighbors ("Plcr") from the Dunbrack's PDFs.
// (it allocates output map memory)
//  C --> central residue index (residue id)
//  L --> left-side residue index (residue id)
//  R --> right-side residue index (residue id)
//  pdfs --> Dunbrack's PDFs
//  size --> Number of side bins in Dunbrack's PDFs (size x size 2D maps)
//  ref --> Reference state map (size x size 2D map)
//          if ref=NULL --> Plcr/1
//          if ref=Pc   --> Plcr/Pc (Pc == independent probabilities for Central "c" residue)
//          if ref=Pt   --> Plcr/Pt (Pt == total Ramachandran plot generated with "ref_gen()")
float **map_gen(int C, int L, int R, float *****pdfs, int size, float **ref)
{
	float **prod;
	float **map;

	prod = MatrixProduct(pdfs[C][0][L],pdfs[C][1][R],size);
	// Bug (4/1/2019)
	// map = MatrixDivision(prod,pdfs[C][1][8],size); // ALL-R [21][1]
	if(ref == NULL)
		map = prod; // For Plcr/1 version
	//map = MatrixDivision(prod,pdfs[C][1][21],size); // L or R = ALL --> [21] (both seem equivalent..., see "rang" program)
	else
		map = MatrixDivision(prod,ref,size); // Using user-supplied reference map

	norm_map(map,size); // Really required? It seems not needed for Paxb/Px version nor Px/Ptot

	//	float dummy=0.0;
	//	for(int i=0; i<size; i++)
	//		for(int j=0; j<size; j++)
	////			dummy += pdfs[C][0][L][i][j];
	//			dummy += pdfs[C][1][21][i][j];
	//	fprintf(stderr,"debug> SUM(pdfs[C][1][21])= %f\n", dummy);

	if(ref != NULL)
	{
		for(int i=0; i<size; i++)
			free(prod[i]);
		free(prod);
	}

	return map;
}

// Generates a "Total" reference map from Dunbrack's tabulated PDFs (allocating output map memory).
//  pdfs --> Dunbrack's PDFs
//  size --> Number of side bins in Dunbrack's PDFs (size x size 2D maps)
//  weights --> Array of molar fractions per residue
float **ref_gen(float *****pdfs, int size, float *weights)
{
	float **map = MatrixInit(size); // Allocate and zero initializes a size x size matrix

	if(weights == NULL)
		for(int i=0; i<20; i++)
			MatrixAccum(map, pdfs[i][1][21], size, 1.0); // R_ALL_map --> pdfs[i][1][21],  L_ALL_map --> pdfs[i][0][21]
	else
		for(int i=0; i<20; i++)
			MatrixAccum(map, pdfs[i][1][21], size, weights[i]); // R_ALL_map --> pdfs[i][1][21],  L_ALL_map --> pdfs[i][0][21]
	norm_map(map,size); // Normalize map to Sum=1.0

	return map;
}

int angle2index(double angle, int size)
// Avoid quantization problems when "angle" is -180 or +180 (-M_PI or +M_PI)
// 	angle = ]-M_PI,+M_PI[
// 	returned index is bounded to [0,size-1]
{
	// double angle0 = angle; // Print luego
	if (angle  >=  M_PI)
		angle += -2*M_PI + 1.0E-6; // Avoid (angle = -M_PI)
	else
		if (angle <= -M_PI)
			angle += 2*M_PI - 1.0E-6; // Avoid (angle = +M_PI)

	int index = angle * size/(2*M_PI) + size/2; // Convert from continuous angle to discrete bin

	return index; // Index is bounded to [0,size-1]
}

// Ramachandran energy derived from neighbor-dependent PDFs from Dunbrack's paper.
float rama_energy(double *dihedral_angle,int nr_atoms_loop,float ***maps,int size)
{
	//float c = 180/M_PI;
	float energy = 0.0;

	for(int i=0; i < (nr_atoms_loop/3)-1; i++) // screen residues
	{
		//		fprintf(stderr,"Phi%d= %6.2f  Psi%d= %6.2f  p= %f\n",i,dihedral_angle[i*3+2]*c,i,dihedral_angle[i*3+3]*c,
		//				maps[ i ][ angle2index(dihedral_angle[i*3+2], size) ][ angle2index(dihedral_angle[i*3+3], size) ] );
		energy -= log10(maps[ i ][ angle2index(dihedral_angle[i*3+2], size) ][ angle2index(dihedral_angle[i*3+3], size) ]);
	}
	//	fprintf(stderr,"Phi%d= %6.2f, Psi%d= %6.2f  p= %f\n",(nr_atoms_loop/3)-1, dihedral_angle[nr_atoms_loop-1]*c, (nr_atoms_loop/3)-1, dihedral_angle[0]*c,
	//			maps[ (nr_atoms_loop/3)-1 ][ angle2index(dihedral_angle[nr_atoms_loop-1], size) ][ angle2index(dihedral_angle[0], size) ] );
	energy -= log10(maps[ (nr_atoms_loop/3)-1 ][ angle2index(dihedral_angle[nr_atoms_loop-1], size) ][ angle2index(dihedral_angle[0], size) ]);

	return energy;
}

// Compute the dot product between two 3D vectors (float)
float dotprod(float *vector1, float *vector2)
{
	float dot = 0.0;
	for(int i = 0; i < 3; i++)
		dot += vector1[i]*vector2[i];
	return dot;
}

// Normalize one 3D vector (float)
float norm(float *bond)
{
	return sqrtf( dotprod(bond,bond) );
}

// Get the dihedrals array (dihedral) for supplied coordinates (co)
//  co --> floats array of atomic coordinates (the first non-zero dihedral defines the 4th atom position)
//  natoms --> Number of atoms of the kinematic chain
//  dihedral --> floats array of dihedral angles (the first 3 are zero by convention)
void finddihedral(float *co, int natoms, float *dihedral)
{
	float angle_dihedral;
	float s;
	float ct;
	float bond1[3];
	float bond2[3];
	float bond3[3];
	float normal1[3];
	float normal2[3];

	//	First three atoms have formally no dihedral angle
	dihedral[0] = 0;
	dihedral[1] = 0;
	dihedral[2] = 0;

	for(int j = 0; j < natoms - 3; j++ )
	{
		/* Calculate vector C,B,C */
		for (int i = 0; i < 3; i++)
		{
			//		  bond1[i] = co[0][i] - co[1][i];
			//		  bond2[i] = co[2][i] - co[1][i]; // Bond that contains current dihedral angle
			//		  bond3[i] = co[2][i] - co[3][i];
			bond1[i] = co[j*3     + i] - co[(j+1)*3 + i];
			bond2[i] = co[(j+2)*3 + i] - co[(j+1)*3 + i]; // Bond that contains current dihedral angle
			bond3[i] = co[(j+2)*3 + i] - co[(j+3)*3 + i];

			//	  if (i==0) printf("->b1 %8.3f %8.3f b2 %8.3f %8.3f b3 %8.3f %8.3f\n",co[0][i],co[1][i],co[2][i],co[1][i],co[2][i],co[3][i]);
		}

		/* Normal to plane 1 */
		normal1[0] = bond1[1] * bond2[2] - bond1[2] * bond2[1];
		normal1[1] = bond1[2] * bond2[0] - bond1[0] * bond2[2];
		normal1[2] = bond1[0] * bond2[1] - bond1[1] * bond2[0];
		/* Normal to plane 2 */
		normal2[0] = bond2[2] * bond3[1] - bond2[1] * bond3[2];
		normal2[1] = bond2[0] * bond3[2] - bond2[2] * bond3[0];
		normal2[2] = bond2[1] * bond3[0] - bond2[0] * bond3[1];

		ct = dotprod(normal1,normal2)/(norm(normal1)*norm(normal2));

		//	In case normalization fails
		if ( ct > 1.0 )
		{
			ct = 1.0;
		}
		else if ( ct < ( -1.0 ) )
		{
			ct = -1.0;
		}

		angle_dihedral = acosf( ct );

		s = bond2[0] * ( normal1[2] * normal2[1] - normal1[1] * normal2[2] ) +
				bond2[1] * ( normal1[0] * normal2[2] - normal1[2] * normal2[0] ) +
				bond2[2] * ( normal1[1] * normal2[0] - normal1[0] * normal2[1] );

		if ( s < 0.0 )
		{
			angle_dihedral = -angle_dihedral;
		}

		angle_dihedral = ( angle_dihedral > 0.0 ) ? M_PI - angle_dihedral : -( M_PI + angle_dihedral );

		//	Correction algorithm to get value in [-Pi,+Pi] range --- BUT IS THIS REALLY NEEDED?
		if (angle_dihedral > M_PI)
		{
			angle_dihedral = angle_dihedral - 2*M_PI;
		}
		else if (angle_dihedral < -M_PI)
		{
			angle_dihedral = angle_dihedral + 2*M_PI;
		}

		dihedral[j + 3] = angle_dihedral;

		//	  co++;
	}

	//delete co;			//	WRONG DON'T DELETE ORIGINAL OBJECT
	//return angle_dihedral ;
}

// Ramachandran energy derived from neighbor-dependent PDFs from Dunbrack's paper.
// It computes Plcr/Pc, Plcr/Pt, or Pc/Pt potentials depending on the reference used in map_gen().
float rama_energy2(float *dihedral_angle, int nres, float ***maps, int size)
{
	//float c = 180/M_PI; // [rad] to [deg]
	float energy = 0.0;

	for(int i=1; i < nres-1; i++) // screen residues
	{
		//		fprintf(stderr,"Phi%d= %6.2f  Psi%d= %6.2f  p= %f\n",i,dihedral_angle[i*3+2]*c,i,dihedral_angle[i*3+3]*c,
		//				maps[ i ][ angle2index(dihedral_angle[i*3+2], size) ][ angle2index(dihedral_angle[i*3+3], size) ] );
		energy -= log10(maps[ i ][ angle2index(dihedral_angle[i*3+2], size) ][ angle2index(dihedral_angle[i*3+3], size) ]);
	}
	//	fprintf(stderr,"Phi%d= %6.2f, Psi%d= %6.2f  p= %f\n",(nr_atoms_loop/3)-1, dihedral_angle[nr_atoms_loop-1]*c, (nr_atoms_loop/3)-1, dihedral_angle[0]*c,
	//			maps[ (nr_atoms_loop/3)-1 ][ angle2index(dihedral_angle[nr_atoms_loop-1], size) ][ angle2index(dihedral_angle[0], size) ] );
	//	energy -= log10(maps[ nres-1 ][ angle2index(dihedral_angle[nres*3-1], size) ][ angle2index(dihedral_angle[0], size) ]);

	return energy;
}

// Ramachandran energy derived from neighbor-dependent PDFs from Dunbrack's paper.
// It computes Plcr/Pc, Plcr/Pt, or Pc/Pt potentials depending on the reference used in map_gen().
float rama_energy2(double *dihedral_angle, int nres, float ***maps, int size)
{
	// float c = 180/M_PI; // [rad] to [deg]
	float energy = 0.0;

	for(int i=1; i < nres-1; i++) // screen residues
	{
		//		fprintf(stderr,"Phi%d= %6.2f  Psi%d= %6.2f  p= %f\n",i,dihedral_angle[i*3+2]*c,i,dihedral_angle[i*3+3]*c,
		//				maps[ i ][ angle2index(dihedral_angle[i*3+2], size) ][ angle2index(dihedral_angle[i*3+3], size) ] );
		energy -= log10(maps[ i ][ angle2index(dihedral_angle[i*3+2], size) ][ angle2index(dihedral_angle[i*3+3], size) ]);
	}
	//	fprintf(stderr,"Phi%d= %6.2f, Psi%d= %6.2f  p= %f\n",(nr_atoms_loop/3)-1, dihedral_angle[nr_atoms_loop-1]*c, (nr_atoms_loop/3)-1, dihedral_angle[0]*c,
	//			maps[ (nr_atoms_loop/3)-1 ][ angle2index(dihedral_angle[nr_atoms_loop-1], size) ][ angle2index(dihedral_angle[0], size) ] );
	//	energy -= log10(maps[ nres-1 ][ angle2index(dihedral_angle[nres*3-1], size) ][ angle2index(dihedral_angle[0], size) ]);

	return energy;
}

// Generate the Complete Ramachandran potential (20x20x20 = 8000 72x72 maps)
//  pdfs       --> Dunbrack's Neighbor dependent PDFs, as read by read_dunbrack()
//  rama_model --> Ramachandran energy model
//  size       --> Sampling of Dunbrack's PDFs (Optiona, 72 by default)
float *****gen_allrama(float *****pdfs, int rama_model, int size)
{
	// AA's molar fraction (from 37k "Pisces3Aid90" proteins)
	float aax[] = {0.0830,0.0134,0.0586,0.0664,0.0407,0.0734,0.0233,0.0579,0.0571,0.0939,0.0180,0.0428,0.0459,0.0371,0.0507,0.0611,0.0554,0.0715,0.0141,0.0357};
	// Data obtained by running "korp v1.63" in mon@yolem:~/korp/Pisces3Aid90NCAC the following command:
	// #> time korp PISCES_id90_r3A.txt -c 16 --nocont --nosas --notest --disk -f 10 --stride .stride -o deleteme
	//
	//	ALA =  0 ; 	//	CYS =  1 ;	//	ASP =  2 ;	//	GLU =  3 ;	//	PHE =  4 ;	//	GLY =  5 ;	//	HIS =  6 ;	//	ILE =  7 ;	//	LYS =  8 ;
	//	LEU =  9 ;	//	MET =  10 ;	//	ASN =  11 ;	//	PRO =  12 ;	//	GLN =  13 ;	//	ARG =  14 ;	//	SER =  15 ;	//	THR =  16 ;	//	VAL =  17 ;
	//	TRP =  18 ;	//	TYR =  19 ;
	//	korp> Total number of interacting frames = 9092777
	//	korp> Frame IDs     =       0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19
	//	korp> Residue count =  754937 122021 532685 604106 370097 667096 211702 526486 519258 853976 163964 389053 416917 337443 460686 555425 503561 650261 128222 324881
	//	korp> Molar fraction=  0.0830 0.0134 0.0586 0.0664 0.0407 0.0734 0.0233 0.0579 0.0571 0.0939 0.0180 0.0428 0.0459 0.0371 0.0507 0.0611 0.0554 0.0715 0.0141 0.0357
	//	korp> TOTAL_________________________________________________________________________
	//	korp> The 9092791 residues from 36851 PDBs were processed.
	//	korp> Total SS [%]:  34.2 Helix,  22.4 Sheet,  18.0 Coil,  20.3 Turn,   1.1 Bridge,   4.0 Helix3/10, and   0.1 Rest
	//	korp> 266284346 contacts found among 9092777 frames (29.3/residue at 16 A cutoff).
	//	korp> Total memory used for contacts: 9586236456 bytes (9142.1 MB)
	//	korp> 266284346 contacts with a 16.00 A cutoff written into deleteme_cont.bin

	float **refmap; // "Total" reference map from Dunbrack's tabulated PDFs (allocating output map memory).
	if(rama_model == 1 || rama_model == 2) // Weighted total probability map required
		refmap = ref_gen(pdfs, size, aax);
	else if(rama_model == 5 || rama_model == 6) // Non-weighted total probability map required
		refmap = ref_gen(pdfs, size, NULL);

	// show_map(refmap,size,"refmap.txt"); // Gnuplot command:  plot "refmap.txt" using 1:2:3 with image

	float *****allmaps; // Complete Ramachandran potential (20x20x20 = 8000 72x72 maps)
	if( !(allmaps = (float *****) malloc( sizeof(float ****) * 20)) )
	{
		printf("Sorry, unable to allocate allmaps memory!!!\n");
		exit(1);
	}
	for(int c = 0; c < NAAS; c++) // Screen Central AAs (NAAS=21 --> 20 std.aas + cis-Pro)
	{
		if( !(allmaps[c] = (float ****) malloc( sizeof(float ***) * 20)) )  // 20 std.aas
		{
			printf("Sorry, unable to allocate allmaps memory!!!\n");
			exit(1);
		}
		for(int l = 0; l < 20; l++) // Screen Left AAs
		{
			if( !(allmaps[c][l] = (float ***) malloc( sizeof(float **) * 20)) )  // 20 std.aas
			{
				printf("Sorry, unable to allocate allmaps memory!!!\n");
				exit(1);
			}
			for(int r = 0; r < 20; r++) // Screen Right AAs
			{
				// fprintf(stderr,"Generating  c=%2d  l=%2d  r=%2d  maps\n",c,l,r);
				switch(rama_model)
				{
				case 0:
					allmaps[c][l][r] = pdfs[c][1][21]; // Px/1 ("21" means the ALL probability for each aminoacid "x")
					break;
				case 1:
				case 5:
					allmaps[c][l][r] = MatrixDivision(pdfs[c][1][21],refmap,size); // Px/Pt or Px/Pt1
					break;
				case 2:
				case 6:
					allmaps[c][l][r] = map_gen(c, l, r, pdfs, size, refmap); // Plcr/Pt or Plcr/Pt1 (Plcr/Pt is slightly worse than Plcr/Pc)
					break;
				case 3:
					allmaps[c][l][r] = map_gen(c, l, r, pdfs, size, pdfs[c][1][21]); // Plcr/Px
					break;
				case 4:
					allmaps[c][l][r] = map_gen(c, l, r, pdfs, size, NULL); // Plcr/1
					break;
				default:
					fprintf(stderr,"Please, insert a valid Ramachandran energy model (%d is not valid). Forcing exit!\n",rama_model);
					exit(1);
				}

				// Plot all 8000 maps
				// char kk[100];
				// sprintf(kk,"allmaps_c%02dl%02dr%02d.txt",c,l,r);
				// show_map(allmaps[c][l][r],size,kk); // Gnuplot command:  plot "mierda_01.txt" using 1:2:3 with image
			}
		}
	}

	if(rama_model == 1 || rama_model == 2) // Total probability map required
		free(refmap);

	return allmaps;
}


//	for(int c = 0; c < 20; c++) // Screen Central AAs
//		for(int l = 0; l < 20; l++) // Screen Left AAs
//			for(int r = 0; r < 20; r++) // Screen Right AAs
//			{
//				fprintf(stderr,"Dumping  c=%2d  l=%2d  r=%2d  maps\n",c,l,r);
//
//				// Plot all 8000 maps
//				 char kk[100];
//				 sprintf(kk,"mierda_c%02dl%02dr%02d.txt",c,l,r);
//				 show_map(allmaps[c][l][r],size,kk); // Gnuplot command:  plot "mierda_01.txt" using 1:2:3 with image
//			}


// Set the cis-Pro identifier (p) into a "nres" residues sequence "seq" from its dihedral angles "dihedrals"
int setCisPro(float *dihedrals, char *seq, int nres)
{
	int ncispro = 0; // Number of cis-Prolines detected
	for(int x=0; x<nres; x++)
		if(seq[x] == 'P' && fabsf(dihedrals[x*3+1]) <= 0.5*M_PI) // if Omega is +- 90 deg (0.5 rad)
		{
			seq[x] = 'p';
			ncispro++; // count number of cis-Prolines
		}
	return ncispro;
}
