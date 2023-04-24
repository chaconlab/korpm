/*
 * pd2.cpp
 * SBG's implementation of the PD2 bump filter (MacDonald et al. 2013)
 *  Created on: Apr 29, 2015
 *      Author: mon
 */

#include "pd2.h"

double pd2[5][5] = {
		// N vs. something...
		{  3.21750,   // N-N
		   3.49664,   // N-CA
		   3.15000,   // N-C
		   2.81000,   // N-O
		   3.52911 }, // N-CB
		// CA vs. something...
		{  3.49664,   // CA-N
		   3.80000,   // CA-CA
		   3.72500,   // CA-C
		   3.11000,   // CA-O
		   3.83529 }, // CA-CB
		// C vs. something...
		{  3.15000,   // C-N
		   3.72500,   // C-CA
		   3.71250,   // C-C
		   2.75000,   // C-O
		   3.79087 }, // C-CB
		// O vs. something...
		{  2.81000,   // O-N
		   3.11000,   // O-CA
		   2.75000,   // O-C
		   2.93040,   // O-O
		   3.17000 }, // O-CB
		// CB vs. something...
		{  3.52911,   // CB-N
		   3.83529,   // CB-CA
		   3.79087,   // CB-C
		   3.17000,   // CB-O
		   3.87090 }  // CB-CB
		};

// The "PD2" bump energy calculation stuff...
//  coord --> Coordinates (1D array)
//  type --> Atom type array (size= number of atoms)
//  res --> Residue index array (size= number of atoms)
//  natoms --> Number of atoms
double pd2_bump(float *coord, char *type, int *res, int natoms)
{

	double radius,radius2;
	double energy = 0.0;

	// FORWARD (i<j) interactions for atom pairs more than FIVE covalent bonds apart
	// -----------------------------------------------------------------------------
	//      Oi  B1        O2  B3
	//      |   |         |   |
	// Ni   Ci  A1   N2   C2  A3
	//  \  / \ /  \ / \  / \ /  \
	//   Ai   N1   C1  A2   N3   C3
	//   |         |   |         |
	//   Bi        O1  B2        O3
	//
	// *Note: A=CA and B=CB
	//
	//  Ni/Bi vs. all BB atoms from residues j >= i+2
	//  Ni/Bi vs. O atoms from residues j >= i+1
	//  Ni/Bi vs. CB atoms from residues j >= i+2
	//  Ai/Oi vs. all BB atoms away from CA (included) of residue j >= i+2
	//  Ai/Oi vs. O atoms from residues j >= i+2
	//  Ai/Oi vs. CB atoms from residues j >= i+2
	//  Ci vs. all BB atoms away from C (included) of residue j >= i+2
	//  Ci vs. O atoms from residues j >= i+2
	//  Ci vs. CB atoms from residues j >= i+2
	//
	int typei,resi;
	float *coordj;
	float coordix,coordiy,coordiz;
	double disa2;
	for(int i=0; i<natoms-4; i++)
	{
		typei = type[i];
		resi = res[i];
		coordix = coord[3*i];
		coordiy = coord[3*i + 1];
		coordiz = coord[3*i + 2];
		if( typei == myN || typei == myCB)
		{
			for(int j=i+4; j<natoms; j++)
			{
				// if( res[j]-res[i] >= gap)
				if( res[j]-resi >= 2 || (res[j]-resi == 1 && type[j] == myO) )
				{
					coordj = coord+3*j;
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ (int) typei ][ (int) type[j] ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						// fprintf(stderr,"i=%4d (%3d)  j=%4d (%3d)  d2= %10.4f  d= %10.4f  radius2= %10.4f\n",i,res[i],j,res[j],dist2,sqrt(dist2),radius2);
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"N/CB i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,res[i],type[i],j,res[j],type[j],disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						// fprintf(stderr,"i %10f %10f %10f  j %10f %10f %10f\n",coordix,coordiy,coordiz,coordj[0],coordj[1],coordj[2]);
					}
				}
			}
		}
		else if( typei == myCA || typei == myO )
		{
			for(int j=i+4; j<natoms; j++)
			{
				if( res[j]-resi >= 3 || (res[j]-resi == 2 && type[j] != myN) )
				// if( res[j]-res[i] >= gap)
				{
					coordj = coord+3*j;
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ (int)type[j] ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						// fprintf(stderr,"i=%4d (%3d)  j=%4d (%3d)  d2= %10.4f  d= %10.4f  radius2= %10.4f\n",i,res[i],j,res[j],dist2,sqrt(dist2),radius2);
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"CA/O i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,res[i],type[i],j,res[j],type[j],disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						// fprintf(stderr,"i %10f %10f %10f  j %10f %10f %10f\n",coordix,coordiy,coordiz,coordj[0],coordj[1],coordj[2]);
					}
				}
			}
		}
		else
		{
			for(int j=i+4; j<natoms; j++)
			{
				if( res[j]-resi >= 3 || (res[j]-resi == 2 && !(type[j] == myN || type[j] == myCA)) )
				// if( res[j]-res[i] >= gap)
				{
					coordj = coord+3*j;
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ (int)type[j] ];
					radius2 = radius*radius;

					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						// fprintf(stderr,"i=%4d (%3d)  j=%4d (%3d)  d2= %10.4f  d= %10.4f  radius2= %10.4f\n",i,res[i],j,res[j],dist2,sqrt(dist2),radius2);
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"C    i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,res[i],type[i],j,res[j],type[j],disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}

				}
			}
		}
	}
	return energy;
}


// The "PD2" bump energy calculation stuff... for the Minimal Backbone Model (N, CA, and C atoms only)
//  coord --> Coordinates (1D array or the N, CA, and C atom coordinates)
//  type --> Atom type array (size= number of atoms)
//  res --> Residue index array (size= number of atoms)
//  natoms --> Number of atoms
double pd2_bumpNCAC(float *coord, char *type, int *res, int natoms)
{
	double radius,radius2;
	double energy = 0.0;

	// FORWARD (i<j) interactions for atom pairs more than FIVE covalent bonds apart
	// -----------------------------------------------------------------------------
	//
	// Ni   Ci  A1   N2   C2  A3
	//  \  / \ /  \ / \  / \ /  \
	//   Ai   N1   C1  A2   N3   C3
	//
	// *Note: A=CA
	//
	//  Ni vs. all BB atoms from residues j >= i+2
	//  Ai vs. all BB atoms away from CA (included) of residue j >= i+2
	//  Ci vs. all BB atoms away from C (included) of residue j >= i+2
	//
	int typei,resi;
	float *coordj;
	float coordix,coordiy,coordiz;
	double disa2;
	for(int i=0; i<natoms-4; i++)
	{
		typei = type[i];
		resi = res[i];
		coordix = coord[3*i];
		coordiy = coord[3*i + 1];
		coordiz = coord[3*i + 2];
		if( typei == myN )
		{
			for(int j=i+4; j<natoms; j++)
			{
				// if( res[j]-res[i] >= gap)
				if( res[j]-resi >= 2 )
				{
					coordj = coord+3*j;
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ (int)type[j] ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						// fprintf(stderr,"i=%4d (%3d)  j=%4d (%3d)  d2= %10.4f  d= %10.4f  radius2= %10.4f\n",i,res[i],j,res[j],dist2,sqrt(dist2),radius2);
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"N/CB i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,res[i],type[i],j,res[j],type[j],disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						// fprintf(stderr,"i %10f %10f %10f  j %10f %10f %10f\n",coordix,coordiy,coordiz,coordj[0],coordj[1],coordj[2]);
					}
				}
			}
		}
		else if( typei == myCA )
		{
			for(int j=i+4; j<natoms; j++)
			{
				if( res[j]-resi >= 3 || (res[j]-resi == 2 && type[j] != myN) )
				// if( res[j]-res[i] >= gap)
				{
					coordj = coord+3*j;
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ (int)type[j] ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						// fprintf(stderr,"i=%4d (%3d)  j=%4d (%3d)  d2= %10.4f  d= %10.4f  radius2= %10.4f\n",i,res[i],j,res[j],dist2,sqrt(dist2),radius2);
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"CA/O i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,res[i],type[i],j,res[j],type[j],disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						// fprintf(stderr,"i %10f %10f %10f  j %10f %10f %10f\n",coordix,coordiy,coordiz,coordj[0],coordj[1],coordj[2]);
					}
				}
			}
		}
		else
		{
			for(int j=i+4; j<natoms; j++)
			{
				if( res[j]-resi >= 3 || (res[j]-resi == 2 && !(type[j] == myN || type[j] == myCA)) )
				// if( res[j]-res[i] >= gap)
				{
					coordj = coord+3*j;
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ (int)type[j] ];
					radius2 = radius*radius;

					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						// fprintf(stderr,"i=%4d (%3d)  j=%4d (%3d)  d2= %10.4f  d= %10.4f  radius2= %10.4f\n",i,res[i],j,res[j],dist2,sqrt(dist2),radius2);
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"C    i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,res[i],type[i],j,res[j],type[j],disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}

				}
			}
		}
	}
	return energy;
}


// The "PD2" bump energy calculation stuff: "Loop vs. Loop"
//  co --> Loop coordinates for N, CA, C atoms
//  pco --> Loop coordinates for O and CB atoms
//  nloop --> Loop number of atoms (N, CA, C only!)
//  cg --> Coarse graining level
//  residuemarker --> Array with the residue types. GLY=1 (size= loop number of atoms, N, CA, C only!)
double pd2_intrabump(double **co, double **pco, int nloop, int cg, int *residuemarker)
{
	// int gap = 6; // Sequential cutoff
	double radius,radius2;
	double energy = 0.0;

	// FORWARD (i<j) interactions for atom pairs more than FIVE covalent bonds apart
	// -----------------------------------------------------------------------------
	//      Oi  B1        O2  B3
	//      |   |         |   |
	// Ni   Ci  A1   N2   C2  A3
	//  \  / \ /  \ / \  / \ /  \
	//   Ai   N1   C1  A2   N3   C3
	//   |         |   |         |
	//   Bi        O1  B2        O3
	//
	// *Note: A=CA and B=CB
	//
	//  Ni/Bi vs. all BB atoms from residues j >= i+2
	//  Ni/Bi vs. O atoms from residues j >= i+1
	//  Ni/Bi vs. CB atoms from residues j >= i+2
	//  Ai/Oi vs. all BB atoms away from CA (included) of residue j >= i+2
	//  Ai/Oi vs. O atoms from residues j >= i+2
	//  Ai/Oi vs. CB atoms from residues j >= i+2
	//  Ci vs. all BB atoms away from C (included) of residue j >= i+2
	//  Ci vs. O atoms from residues j >= i+2
	//  Ci vs. CB atoms from residues j >= i+2
	//
	int typei,typej,resi,resj;
	double *coordj;
	float coordix,coordiy,coordiz;
	double disa2;

	for(int i=3; i<nloop-3; i++) // Screen loop BB, i.e. "co"
	{
		// N,CA,C vs. N,CA,C (and O,CB)
		typei = i%3; // 0=N, 1=CA, 2=C
		resi = i/3; // sequential residue index
		coordix = co[i][0];
		coordiy = co[i][1];
		coordiz = co[i][2];
		if( typei == myN )
		{
			// N vs. N,CA,C
			for(int j=i+6; j<nloop-3; j++) // Screen loop
			{
				coordj = co[j];
				typej = j%3; // 0=N, 1=CA, 2=C
				resj = j/3; // sequential residue index
				if( resj-resi >= 2 )
				{
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"N    vs. CO  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}
				}
			}

			// N vs. "pco"'s O
			if(cg == 1 || cg == 3) //  1: N,CA,C,O    3: N,CA,C,O,CB
			{
				typej = 3; // 3=O
				for(int j=i+5; j<nloop-3; j+=3) // Screen all O's 1 residue away
				{
					coordj = pco[j];
					// resj = j/3; // sequential residue index
					// if( resj-resi >= 1 )
					{
						// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
						radius = pd2[ typei ][ typej ];
						radius2 = radius*radius;
						if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
								(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
								(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
						{
							energy += radius - ( disa2 / radius );
							// fprintf(stderr,"N    vs. O   i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						}
					}
				}
			}

			// N vs. "pco"'s CB
			if(cg >= 2) //  2: N,CA,C,CB    3: N,CA,C,O,CB
			{
				typej = 4; // 4=CB
				for(int j=i+7; j<nloop-3; j+=3) // Screen all CB's 2 residues away
				{
					if(residuemarker[j/3] != 1) // If non-GLY, it has CB
					{
						// resj = j/3; // sequential residue index
						coordj = pco[j];
						// if( resj-resi >= 2 )
						{
							// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
							radius = pd2[ typei ][ typej ];
							radius2 = radius*radius;
							if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
									(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
									(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
							{
								energy += radius - ( disa2 / radius );
								// fprintf(stderr,"N    vs. CB  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
							}
						}
					}
				}
			}
		}
		else if( typei == myCA )
		{
			// CA vs. N,CA,C
			for(int j=i+6; j<nloop-3; j++) // Screen loop BB atoms 2 residues away (non N)
			{
				coordj = co[j];
				typej = j%3; // 0=N, 1=CA, 2=C
				// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
				radius = pd2[ typei ][ typej ];
				radius2 = radius*radius;
				if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
					(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
					(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
				{
					energy += radius - ( disa2 / radius );
					// fprintf(stderr,"CA   vs. CO  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
				}
			}

			// CA vs. "pco"'s O
			if(cg == 1 || cg == 3) //  1: N,CA,C,O    3: N,CA,C,O,CB
			{
				typej = 3; // 3=O
				for(int j=i+7; j<nloop-3; j+=3) // Screen all O's 2 residues away
				{
					coordj = pco[j];
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"CA   vs. O   i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}
				}
			}

			// CA vs. "pco"'s CB
			if(cg >= 2) //  2: N,CA,C,CB    3: N,CA,C,O,CB
			{
				typej = 4; // 4=CB
				for(int j=i+6; j<nloop-3; j+=3) // Screen all CB's 2 residues away
				{
					if(residuemarker[j/3] != 1) // If non-GLY, it has CB
					{
						coordj = pco[j];
						// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
						radius = pd2[ typei ][ typej ];
						radius2 = radius*radius;
						if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
							(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
							(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
						{
							energy += radius - ( disa2 / radius );
							// fprintf(stderr,"CA   vs. CB  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						}
					}
				}
			}
		}
		else // for "C"
		{
			// C vs. N,CA,C
			for(int j=i+6; j<nloop-3; j++) // Screen loop BB atoms 2 residues away (non-N and non-CA)
			{
				coordj = co[j];
				typej = j%3; // 0=N, 1=CA, 2=C
				// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
				// fprintf(stderr,"typei= %d  typej= %d  pd2= %f\n",typei,typej,pd2[ typei ][ typej ]);
				radius = pd2[ typei ][ typej ];
				radius2 = radius*radius;
				if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
					(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
					(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
				{
					energy += radius - ( disa2 / radius );
					// fprintf(stderr,"C    vs. CO  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
				}
			}

			// C vs. "pco"'s O
			if(cg == 1 || cg == 3) //  1: N,CA,C,O    3: N,CA,C,O,CB
			{
				typej = 3; // 3=O
				for(int j=i+6; j<nloop-3; j+=3) // Screen all O's 2 residues away
				{
					coordj = pco[j];
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"C    vs. O   i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}
				}
			}

			// C vs. "pco"'s CB
			if(cg >= 2) //  2: N,CA,C,CB    3: N,CA,C,O,CB
			{
				typej = 4; // 4=CB
				for(int j=i+5; j<nloop-3; j+=3) // Screen all CB's 2 residues away
				{
					if(residuemarker[j/3] != 1) // If non-GLY, it has CB
					{
						coordj = pco[j];
						// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
						radius = pd2[ typei ][ typej ];
						radius2 = radius*radius;
						if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
							(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
							(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
						{
							energy += radius - ( disa2 / radius );
							// fprintf(stderr,"C    vs. CB  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						}
					}
				}
			}
		}

		// O vs. O and CB
		//  Ai/Oi vs. all BB atoms away from CA (included) of residue j >= i+2
		//  Ai/Oi vs. O atoms from residues j >= i+2
		//  Ai/Oi vs. CB atoms from residues j >= i+2
		if( typei == myC && (cg == 1 || cg == 3) ) // C has O    //  cg= 1: N,CA,C,O  3: N,CA,C,O,CB
		{
			typei = 3; // 3=O
			coordix = pco[i][0];
			coordiy = pco[i][1];
			coordiz = pco[i][2];

			// O vs. N,CA,C
			for(int j=i+5; j<nloop-3; j++) // Screen loop BB atoms 2 residues away (non-N)
			{
				coordj = co[j];
				typej = j%3; // 0=N, 1=CA, 2=C
				// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
				// fprintf(stderr,"typei= %d  typej= %d  pd2= %f\n",typei,typej,pd2[ typei ][ typej ]);
				radius = pd2[ typei ][ typej ];
				radius2 = radius*radius;
				if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
					(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
					(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
				{
					energy += radius - ( disa2 / radius );
					// fprintf(stderr,"C    vs. CO  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
				}
			}

			// O vs. "pco"'s O
			typej = 3; // 3=O
			for(int j=i+6; j<nloop-3; j+=3) // Screen O's 2 residues away
			{
				coordj = pco[j];
				// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
				radius = pd2[ typei ][ typej ];
				radius2 = radius*radius;
				if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
					(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
					(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
				{
					energy += radius - ( disa2 / radius );
					// fprintf(stderr,"O    vs. O   i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
				}
			}
			// O vs. "pco"'s CB
			typej = 4; // 4=CB
			for(int j=i+5; j<nloop-3; j+=3) // Screen CB's 2 residues away
			{
				if(residuemarker[j/3] != 1) // If non-GLY, it has CB
				{
					coordj = pco[j];
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"O    vs. CB  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}
				}
			}
		}

		// CB vs. O and CB
		//  Ni/Bi vs. all BB atoms from residues j >= i+2
		//  Ni/Bi vs. O atoms from residues j >= i+1
		//  Ni/Bi vs. CB atoms from residues j >= i+2
		if(typei == myCA && residuemarker[i/3] != 1 && cg >= 2)  // CA has CB if non-GLY  //  cg= 2: N,CA,C,CB  3: N,CA,C,O,CB
		{
			typei = 4; // 4=CB
			coordix = pco[i][0];
			coordiy = pco[i][1];
			coordiz = pco[i][2];

			// CB vs. N,CA,C
			for(int j=i+5; j<nloop-3; j++) // Screen loop BB atoms 2 residues away
			{
				coordj = co[j];
				typej = j%3; // 0=N, 1=CA, 2=C
				// resj = j/3; // sequential residue index
				// if( resj-resi >= 2 )
				{
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"CB   vs. CO  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}
				}
			}

			// CB vs. "pco"'s O
			typej = 3; // 3=O
			for(int j=i+4; j<nloop-3; j+=3) // Screen O's 1 residues away
			{
				coordj = pco[j];
				// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
				radius = pd2[ typei ][ typej ];
				radius2 = radius*radius;
				if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
					(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
					(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
				{
					energy += radius - ( disa2 / radius );
					// fprintf(stderr,"CB   vs. O   i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
				}
			}
			// CB vs. "pco"'s CB
			typej = 4; // 4=CB
			for(int j=i+6; j<nloop-3; j+=3) // Screen CB's 2 residues away
			{
				if(residuemarker[j/3] != 1) // If non-GLY, it has CB
				{
					coordj = pco[j];
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"CB   vs. CB  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,j/3,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}
				}
			}
		}
	}
	return 100*energy;
}


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
double pd2_extrabump(float *coord, char *type, int *res, int npocket, double **co, double **pco, int nloop, int loopNt, int cg, int *residuemarker)
{
	// int gap = 6; // Sequential cutoff
	double radius,radius2;
	double energy = 0.0;

	// FORWARD (i<j) interactions for atom pairs more than FIVE covalent bonds apart
	// -----------------------------------------------------------------------------
	//      Oi  B1        O2  B3
	//      |   |         |   |
	// Ni   Ci  A1   N2   C2  A3
	//  \  / \ /  \ / \  / \ /  \
	//   Ai   N1   C1  A2   N3   C3
	//   |         |   |         |
	//   Bi        O1  B2        O3
	//
	// *Note: A=CA and B=CB
	//
	//  Ni/Bi vs. all BB atoms from residues j >= i+2
	//  Ni/Bi vs. O atoms from residues j >= i+1
	//  Ni/Bi vs. CB atoms from residues j >= i+2
	//  Ai/Oi vs. all BB atoms away from CA (included) of residue j >= i+2
	//  Ai/Oi vs. O atoms from residues j >= i+2
	//  Ai/Oi vs. CB atoms from residues j >= i+2
	//  Ci vs. all BB atoms away from C (included) of residue j >= i+2
	//  Ci vs. O atoms from residues j >= i+2
	//  Ci vs. CB atoms from residues j >= i+2
	//
	// BACKWARD (i>j) interactions for atom pairs more than FIVE covalent bonds apart
	// ------------------------------------------------------------------------------
	//      O3  B2        O1  Bi
	//      |   |         |   |
	// N3   C3  A2   N1   C1  Ai
	//  \  / \ /  \ / \  / \ /  \
	//   A3   N2   C2  A1   Ni   Ci
	//   |         |   |         |
	//   B3        O2  B1        Oi
	//
	// *Note: A=CA and B=CB
	//
	//  Ni vs. all BB atoms away from N (included) of residue j <= i-2
	//  Ni vs. O atoms from residues j <= i-3
	//  Ni vs. CB atoms from residues j <= i-2
	//  Ai vs. all BB atoms away from CA (included) of residue j <= i-2
	//  Ai vs. O atoms from residues j <= i-2
	//  Ai vs. CB atoms from residues j <= i-2
	//  Ci/Bi vs. all atoms from residue j <= i-2
	//  Oi vs. all BB atoms away from N (included) of residue j <= i-1
	//  Oi vs. O atoms from residues j <= i-1
	//  Oi vs. CB atoms from residues j <= i-1
	//
	int typei,typej,resi,resj;
	double *coordj;
	float coordix,coordiy,coordiz;
	double disa2;
	for(int i=0; i<npocket; i++) // Screen pocket
	{
		typei = type[i];
		resi = res[i];
		if(resi >= loopNt)
			resi += nloop/3-2; // adds the number of loop residues
		coordix = coord[3*i];
		coordiy = coord[3*i + 1];
		coordiz = coord[3*i + 2];
		if( typei == myN || typei == myCB )
		{
			// Pocket vs. "co"
			for(int j=3; j<nloop-3; j++) // Screen loop
			{
				resj = loopNt + j/3 - 1;
				// if( res[j]-res[i] >= gap)
				// if( abs(resj-resi) >= 2 )
				if( resj-resi >= 2 || (resi-resj >= 3 || (resi-resj == 2 && typej != myCA && typej != myC) ) )
				{
					coordj = co[j];
					typej = j%3; // 0=N, 1=CA, 2=C
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"N/CB i= %4d  %3d  %d  j= %4d  %3d %d  d2= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,radius2, radius - ( disa2 / radius ), energy);
						// fprintf(stderr,"P:N/CB vs. CO  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						// fprintf(stderr,"i %10f %10f %10f  j %10f %10f %10f\n",coordix,coordiy,coordiz,coordj[0],coordj[1],coordj[2]);
					}
				}
			}

			// Pocket vs. "pco"'s O
			if(cg == 1 || cg == 3) //  1: N,CA,C,O    3: N,CA,C,O,CB
			{
				typej = 3; // 3=O
				for(int j=5; j<nloop-3; j+=3) // Screen O's
				{
					resj = loopNt + j/3 - 1;
					// if( res[j]-res[i] >= gap)
					// if( resj-resi >= 1 || (resi-resj >= 3 || (resi-resj == 2 && typej != myCA && typej != myC) ))
					if( resj-resi >= 1 || (resi-resj >= 3 && typei == myN) || (resi-resj >= 2 && typei == myCB) )
					{
						coordj = pco[j];
						// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
						radius = pd2[ typei ][ typej ];
						radius2 = radius*radius;
						if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
								(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
								(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
						{
							energy += radius - ( disa2 / radius );
							// fprintf(stderr,"P:N/CB vs. O   i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
							// fprintf(stderr,"i %10f %10f %10f  j %10f %10f %10f\n",coordix,coordiy,coordiz,coordj[0],coordj[1],coordj[2]);
						}
					}
				}
			}

			// Pocket vs. "pco"'s CB
			if(cg >= 2) //  2: N,CA,C,CB    3: N,CA,C,O,CB
			{
				typej = 4; // 4=CB
				for(int j=4; j<nloop-3; j+=3) // Screen CB's
				{
					if(residuemarker[j/3] != 1) // If non-GLY, it has CB
					{
						resj = loopNt + j/3 - 1;
						// if( res[j]-res[i] >= gap)
						//if( abs(resj-resi) >= 2 )
						if( resj-resi >= 2 || resi-resj >= 2)
						{
							coordj = pco[j];
							// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
							radius = pd2[ typei ][ typej ];
							radius2 = radius*radius;
							if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
									(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
									(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
							{
								energy += radius - ( disa2 / radius );
								// fprintf(stderr,"P:N/CB vs. CB  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
							}
						}
					}
				}
			}
		}
		else if( typei == myCA || typei == myO )
		{
			// Pocket vs. "co"
			for(int j=3; j<nloop-3; j++) // Screen loop
			{
				typej = j%3; // 0=N, 1=CA, 2=C
				resj = loopNt + j/3 - 1;
				// if( res[j]-res[i] >= gap)
				if( (resj-resi >= 3 || (resj-resi == 2 && typej != myN))  ||
					(typei == myCA && (resi-resj >= 3 || (resi-resj == 2 && typej != myC))) ||
					(typei == myO && (resi-resj >= 2 || (resi-resj == 1 && typej != myCA && typej != myC))) )
				{
					coordj = co[j];
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"P:CA/O vs. CO  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}
				}
			}

			// Pocket vs. "pco"'s O
			if(cg == 1 || cg == 3) //  1: N,CA,C,O    3: N,CA,C,O,CB
			{
				typej = 3; // 3=O
				for(int j=5; j<nloop-3; j+=3) // Screen O's
				{
					resj = loopNt + j/3 - 1;
					// if( res[j]-res[i] >= gap)
					if( resj-resi >= 2 || (typei == myCA && resi-resj >= 2) || (typei == myO && resi-resj >= 1) )
					{
						coordj = pco[j];
						// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
						radius = pd2[ typei ][ typej ];
						radius2 = radius*radius;
						if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
							(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
							(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
						{
							energy += radius - ( disa2 / radius );
							// fprintf(stderr,"P:CA/O vs. O   i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						}
					}
				}
			}

			// Pocket vs. "pco"'s CB
			if(cg >= 2) //  2: N,CA,C,CB    3: N,CA,C,O,CB
			{
				typej = 4; // 4=CB
				for(int j=4; j<nloop-3; j+=3) // Screen CB's
				{
					if(residuemarker[j/3] != 1) // If non-GLY, it has CB
					{
						resj = loopNt + j/3 - 1;
						// if( res[j]-res[i] >= gap)
						if( resj-resi >= 2 || (typei == myCA && resi-resj >= 2) || (typei == myO && resi-resj >= 1) )
						{
							coordj = pco[j];
							// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
							radius = pd2[ typei ][ typej ];
							radius2 = radius*radius;
							if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
								(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
								(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
							{
								energy += radius - ( disa2 / radius );
								// fprintf(stderr,"P:CA/O vs. CB  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
							}
						}
					}
				}
			}
		}
		else // for "C"
		{
			// Pocket vs. "co"
			for(int j=3; j<nloop-3; j++) // Screen loop
			{
				typej = j%3; // 0=N, 1=CA, 2=C
				resj = loopNt + j/3 - 1;
				// if( res[j]-res[i] >= gap)
				// if( resj-resi >= 3 || (resj-resi == 2 && !(typej == myN || typej == myCA)) )
				if( (resj-resi >= 3 || (resj-resi == 2 && typej == myC)) || (resi-resj >= 2))
				{
					coordj = co[j];
					// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
					// fprintf(stderr,"typei= %d  typej= %d  pd2= %f\n",typei,typej,pd2[ typei ][ typej ]);
					radius = pd2[ typei ][ typej ];
					radius2 = radius*radius;
					if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
						(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
						(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
					{
						energy += radius - ( disa2 / radius );
						// fprintf(stderr,"P:C    vs. CO  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
					}
				}
			}

			// Pocket vs. "pco"'s O
			if(cg == 1 || cg == 3) //  1: N,CA,C,O    3: N,CA,C,O,CB
			{
				typej = 3; // 3=O
				for(int j=5; j<nloop-3; j+=3) // Screen O's
				{
					resj = loopNt + j/3 - 1;
					// if( res[j]-res[i] >= gap)
					if( resj-resi >= 2 || resi-resj >= 2 )
					{
						coordj = pco[j];
						// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
						radius = pd2[ typei ][ typej ];
						radius2 = radius*radius;
						if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
							(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
							(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
						{
							energy += radius - ( disa2 / radius );
							// fprintf(stderr,"P:C    vs. O   i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
						}
					}
				}
			}

			// Pocket vs. "pco"'s CB
			if(cg >= 2) //  2: N,CA,C,CB    3: N,CA,C,O,CB
			{
				typej = 4; // 4=CB
				for(int j=4; j<nloop-3; j+=3) // Screen CB's
				{
					if(residuemarker[j/3] != 1) // If non-GLY, it has CB
					{
						resj = loopNt + j/3 - 1;
						// if( res[j]-res[i] >= gap)
						if( resj-resi >= 2 || resi-resj >= 2)
						{
							coordj = pco[j];
							// PD2 "soft bump" energy term formula: (radius_sq - (dist_sq)) / radius     // Equivalent to:  radius - (dist_sq/radius)
							radius = pd2[ typei ][ typej ];
							radius2 = radius*radius;
							if( (disa2 =  pow(coordix-coordj[0],2)) < radius2 &&
								(disa2 += pow(coordiy-coordj[1],2)) < radius2 &&
								(disa2 += pow(coordiz-coordj[2],2)) < radius2 )
							{
								energy += radius - ( disa2 / radius );
								// fprintf(stderr,"P:C    vs. CB  i= %4d  %3d  %d  j= %4d  %3d  %d  d2= %10.4f  d= %10.4f  radius2= %10.4f PD2= %10.4f Accum= %10.4f\n",i,resi,typei,j,resj,typej,disa2,sqrt(disa2),radius2,(radius2 - disa2)/sqrt(radius2), energy);
							}
						}
					}
				}
			}
		}
	}
	return 100*energy;
}


// Creates an array with the atomic type (N,CA,etc..) for the whole Macromolecule ("PD2" energy calculation stuff)
// p_type --> Pointer to the "type array" (*p_type == NULL for automatic memory allocation)
void pd2_type(Macromolecule *mol, char **p_type)
{
	char *type;
	int natoms;
	natoms = mol->get_num_atoms();
	if(*p_type == NULL) // Allocate memory (if necessary)
	{
		type = (char *) malloc(sizeof(char) * natoms);
		*p_type = type;
	}

	pdbIter *iter_atom = new pdbIter(mol);
	Atom *atom;
	for( iter_atom->pos_atom = 0; !iter_atom->gend_atom(); iter_atom->next_atom() )
	{
		atom = iter_atom->get_atom();
		// fprintf(stderr,"pos_atom= %4d  name= %s\n",iter_atom->pos_atom,atom->getName());
		if(strncmp(atom->getName()+1,"N ",2) == 0)
			type[iter_atom->pos_atom] = myN;
		else if(strncmp(atom->getName()+1,"CA",2) == 0)
			type[iter_atom->pos_atom] = myCA;
		else if(strncmp(atom->getName()+1,"C ",2) == 0)
			type[iter_atom->pos_atom] = myC;
		else if(strncmp(atom->getName()+1,"O ",2) == 0)
			type[iter_atom->pos_atom] = myO;
		else if(strncmp(atom->getName()+1,"CB",2) == 0)
			type[iter_atom->pos_atom] = myCB;
	}
	delete iter_atom;
}

// Creates an array with the residue index of the whole Macromolecule ("PD2" energy calculation stuff)
// p_res --> Pointer to the "residue index array" (*p_res == NULL for automatic memory allocation)
void pd2_res(Macromolecule *mol, int **p_res)
{
	int *res;
	int natoms;
	natoms = mol->get_num_atoms();
	if(*p_res == NULL) // Allocate memory (if necessary)
	{
		res = (int *) malloc(sizeof(int) * natoms);
		*p_res = res;
	}

	pdbIter *iter_frag = new pdbIter(mol);
	pdbIter *iter_atom;
	Fragment *frag;

	int cont=0; // index counter
	for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
	{
		frag = iter_frag->get_fragment(); // Get residue
		iter_atom = new pdbIter(frag);
		for( iter_atom->pos_atom = 0; !iter_atom->gend_atom(); iter_atom->next_atom() )
		{
			res[cont] = iter_frag->pos_fragment;
			cont++; // atom index
		}
		delete iter_atom;
	}
	delete iter_frag;
}

// Returns the internal residue index (pos_fragment) corresponding to the input index (resnum PDB) ("PD2" energy calculation stuff)
//  resnum --> Residue index in PDB
int resindex(Macromolecule *mol, int resnum)
{
	int buf = -1;

	pdbIter *iter_frag = new pdbIter(mol);
	Fragment *frag;

	for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
	{
		frag = iter_frag->get_fragment(); // Get residue
		if(frag->getIdNumber() == resnum)
		{
			buf = iter_frag->pos_fragment;
			delete iter_frag;
			return buf;
		}
	}
	delete iter_frag;
	return buf;
}

