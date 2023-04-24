/*
 * tobi.cpp
 *
 * Implementation of the pairwise potential described in:
 * Designing coarse grained-and atom based-potentials for protein-protein docking.
 * Dror Tobi. BMC Structural Biology (2010).
 * DOI: 10.1186/1472-6807-10-40. Published: 15 November 2010
 *
 * This code was adapted from Pablo's implementation of TOBI energy found in Frodock.
 *
 *  Created on: Feb 8, 2017
 *      Author: mon
 */

#include <libpdb/include/Macromolecule.h>
#include "tobi.h"

// WARNING: TOBI energy requires that in MacroInfo.cpp there exist the following method:
// int Macromolecule::pdbmatrices(float **coord, int **nres, int **pfirst, int **cas, int **tpatom)

// TOBI energy calculation "all in one" wrapper.
//  *molR and *molL  --> Receptor and Ligand Macromolecules
//  tobi             --> 1 or 2 for Tobi-1 (ADP-I) or Tobi-2 (ADP-II) energy models (2 by default)
float tobi_energy(Macromolecule *molR, Macromolecule *molL, char tobi)
{
	float energy;
	int Rnumres, *Rnres, *Rpfirst, *Rcas, *Rtpatom, Lnumres, *Lnres, *Lpfirst, *Lcas, *Ltpatom;
	float *Rcoord, *Lcoord;

	// Getting coordinates (for all atoms and CAs), number of atoms per residue, and atom types for TOBI energy.
	Rnumres = molR->pdbmatrices(&Rcoord, &Rnres, &Rpfirst, &Rcas, &Rtpatom);
	Lnumres = molL->pdbmatrices(&Lcoord, &Lnres, &Lpfirst, &Lcas, &Ltpatom);

	// The TOBI energy calculation.
	energy = tobi_energy(Rcoord, Rnumres, Rnres, Rpfirst, Rcas, Rtpatom, Lcoord, Lnumres, Lnres, Lpfirst, Lcas, Ltpatom);

	// Free all working arrays
	free(Rnres);
	free(Rpfirst);
	free(Rcas);
	free(Rtpatom);
	free(Rcoord);
	free(Lnres);
	free(Lpfirst);
	free(Lcas);
	free(Ltpatom);
	free(Lcoord);

	return energy;
}


// TOBI energy calculation.
//  *R/Lcoord    --> Atomic coordinates (R or L stand for either Receptor or Ligand)
//  R/Lnumres    --> Number of residues
//  *R/Lnres     --> Number of atoms for each residue
//  *R/Lpfirst   --> Index of the first atom for each residue
//  *R/Lcas      --> Index of the CA atom for each residue
//  *R/Ltpatom   --> Tobi's atom type
//  tobi         --> 1 or 2 for Tobi-1 (ADP-I) or Tobi-2 (ADP-II) energy models (2 by default)
float tobi_energy(float *Rcoord, int Rnumres, int *Rnres, int *Rpfirst, int *Rcas, int *Rtpatom,
		float *Lcoord, int Lnumres, int *Lnres, int *Lpfirst, int *Lcas, int *Ltpatom, char tobi)
{
	double pot = 0.0;
	//	float dist, distLx, distLy, distLpx, distLpy, dlimit=6*6+2, dlimitCA=24*24+1;
	float dist, distLx, distLy, distLpx, distLpy;
	float dlimit = 6*6+2;
	float dlimitCA = 24*24+1;
	float px,py,pz, rx,ry,rz, Rpx, Rpy, Rpz, px2,py2,pz2;

	for (int i1 = 0; i1 < Rnumres; i1++ )
	{
		// pos CA receptor
		rx = *(Rcoord+Rcas[i1]);
		ry = *(Rcoord+Rcas[i1]+1);
		rz = *(Rcoord+Rcas[i1]+2);

		for (int i2 = 0; i2 < Lnumres; i2++ )
		{
			// pos CA Ligand
			px = *(Lcoord+Lcas[i2]);
			py = *(Lcoord+Lcas[i2]+1);
			pz = *(Lcoord+Lcas[i2]+2);

//			Lx= px*matrix[0][0] + py*matrix[0][1] + pz*matrix[0][2];
//			Lx+=point[0]+0.5+origin_receptor.x();

			//  if ((i1==0)&&(i2==0))
			//  fprintf(stderr," coord L %f %f %f  %f %f %f    %f %f %f\n",  px, py, pz, Lx, Ly, Lz, rx, ry, rz);

//			if ( ((distLx=(Lx-rx)*(Lx-rx))) < dlimitCA)
			if ( ((distLx = (px-rx)*(px-rx))) < dlimitCA)
			{

//				Ly= px*matrix[1][0] + py*matrix[1][1] + pz*matrix[1][2];
//				Ly+=point[1]+0.5+origin_receptor.y();

//				if ( ((distLy=distLx+ (Ly-ry)*(Ly-ry)))  < dlimitCA)
				if ( ((distLy = distLx+ (py-ry)*(py-ry)))  < dlimitCA)
				{

//					Lz= px*matrix[2][0] + py*matrix[2][1] + pz*matrix[2][2];
//					Lz+=point[2]+0.5+origin_receptor.z();

					// dist= (Lx-px)*(Lx-px) +(Ly-py)*(Ly-py) +(Lz-pz)*(Lz-pz) ;

//					if ( (distLy + (Lz-rz)*(Lz-rz) ) < dlimitCA)
					if ( (distLy + (pz-rz)*(pz-rz) ) < dlimitCA)
					{

						// fprintf(stderr,"-%7.3f res %d %d  %7.3f %7.3f %7.3f  %7.3f %7.3f %7.3f\n",
						//		distLy + (Lz-rz)*(Lz-rz), i1, i2, Lx, Ly, Lz, rx, ry, rz);

						// Amino acids at less than dlimit

						for (int i3 = 0; i3 < Rnres[i1]; i3++ )
						{

							// receptor atoms
							Rpx = *(Rcoord+Rpfirst[i1]+3*i3);
							Rpy = *(Rcoord+Rpfirst[i1]+1+3*i3);
							Rpz = *(Rcoord+Rpfirst[i1]+2+3*i3);

							for (int i4 = 0; i4 < Lnres[i2]; i4++ )
							{

								// Ligand atoms
								px2 = *(Lcoord+Lpfirst[i2]+3*i4);
								py2 = *(Lcoord+Lpfirst[i2]+1+3*i4);
								pz2 = *(Lcoord+Lpfirst[i2]+2+3*i4);

//								Lpx= px2*matrix[0][0] + py2*matrix[0][1] + pz2*matrix[0][2];
//								Lpx+=point[0]+0.5+origin_receptor.x();

//								if ( ((distLpx=(Lpx-Rpx)*(Lpx-Rpx))) <= dlimit)
								if ( ((distLpx = (px2-Rpx)*(px2-Rpx))) <= dlimit)
								{
//									Lpy= px2*matrix[1][0] + py2*matrix[1][1] + pz2*matrix[1][2];
//									Lpy+=point[1]+0.5+origin_receptor.y();

//									if ( ((distLpy=distLpx+ (Lpy-Rpy)*(Lpy-Rpy)))  <= dlimit)
									if ( ((distLpy = distLpx+ (py2-Rpy)*(py2-Rpy))) <= dlimit)
									{
//										Lpz= px2*matrix[2][0] + py2*matrix[2][1] + pz2*matrix[2][2];
//										Lpz+=point[2]+0.5+origin_receptor.z();

//										if ( ((dist = distLpy + (Lpz-Rpz)*(Lpz-Rpz)) ) <= dlimit)
										if ( ((dist = distLpy + (pz2-Rpz)*(pz2-Rpz)) ) <= dlimit)
										{
											if(tobi == 1) // Tobi-1
											{
												if (dist <= 16.0) // 16=4*4
													pot += tobi1_one[Rtpatom[Rpfirst[i1]/3+i3]][Ltpatom[Lpfirst[i2]/3+i4]];
												else if (dist <= 36.0) // 36=6*6
													pot += tobi1_two[Rtpatom[Rpfirst[i1]/3+i3]][Ltpatom[Lpfirst[i2]/3+i4]];
											}
											else // Tobi-2
											{
												if (dist <= 16.0) // 16=4*4
													pot += tobi2_one[Rtpatom[Rpfirst[i1]/3+i3]][Ltpatom[Lpfirst[i2]/3+i4]];
												else if (dist <= 36.0) // 36=6*6
													pot += tobi2_two[Rtpatom[Rpfirst[i1]/3+i3]][Ltpatom[Lpfirst[i2]/3+i4]];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return (float)pot;
}
