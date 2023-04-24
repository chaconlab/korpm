/*
 * icosa.cpp


 *
 *  Created on: Apr 29, 2015
 *      Author: mon
 */
#include <icosa.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// AA list
char aa[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};

// Icosaedral faces (triangles) list defined from icosVertices
int icosTri[][3] = {
		{0, 2, 1},
		{0, 3, 2},
		{0, 4, 3},
		{0, 5, 4},
		{0, 1, 5},
		{1, 2, 7},
		{2, 3, 8},
		{3, 4, 9},
		{4, 5, 10},
		{5, 1, 6},
		{1, 7, 6},
		{2, 8, 7},
		{3, 9, 8},
		{4, 10, 9},
		{5, 6, 10},
		{6, 7, 11},
		{7, 8, 11},
		{8, 9, 11},
		{9, 10, 11},
		{10, 6, 11}
};

// Read loco energy from text file
void read_loco(double ****loco, char *name)
{
	FILE *f_file;
	char line[150];
	char aai,aaj;
	int tri,dist,nlines=0;
	double ener;

	if( !(f_file = fopen(name, "r") ) )
	{
		printf("Msg(read_loco): I'm sorry, unable to open %s file\nForcing exit\n",name);
		exit(1);
	}

	while( fgets(line, 150, f_file) != NULL )
	{
		// fprintf(stderr,"%s\n",line);
//		sscanf(line,"%c%c %d %d %f",&aai,&aaj,&dist,&tri,&ener);
		sscanf(line,"%c%c %d %d %lf",&aai,&aaj,&dist,&tri,&ener);
//		fprintf(stderr,"aai= %c  aaj= %c  dist= %d  tri= %d  ener= %f\n",aai,aaj,dist,tri,ener);
		loco[aa2index(aai,aa)][aa2index(aaj,aa)][dist][tri] = ener;
		nlines++;
	}

	if(nlines != NAASloco*NAASloco*NDIS*NTRI)
	{
		fprintf(stderr,"Some mismatch occurred in %s (%d found but %d expected)\n",name,nlines,NAASloco*NAASloco*NDIS*NTRI);
		exit(1);
	}

	fclose(f_file);
}

// Return aminoacid index given its 1-letter code character.
int aa2index(char aa, char *list)
{
	for(int i=0; i<NAASloco; i++)
//		if(aa==list[i])
		if(aa == list[i] || aa-32 == list[i]) // this handles lowercase characters...
			return i;
	fprintf(stderr,"\naa2index> Sorry, %c aminoacid not found!!! Forcing exit!\n\n",aa);
	exit(1);
	return -1;
}

// Initialize icosahedron data
void init_icos(double icosVertices[][3])
{
  double theta = 26.56505117707799 * PI / 180.0;
//  double stheta = sinf(theta);
//  double ctheta = cosf(theta);
  double stheta = sin(theta);
  double ctheta = cos(theta);
  double phi = PI/5.0;

  // Lower vertex (1 vertex)
  icosVertices[0][0] = 0.0;
  icosVertices[0][1] = 0.0;
  icosVertices[0][2] = -1.0;

  // Lower pentagon (5 vertices)
  for (int i=1; i<6; i++)
  {
	  // push @icosVertices, [cos($phi)*$ctheta, sin($phi)*$ctheta, -$stheta];
//	  icosVertices[i][0] = cosf(phi)*ctheta;
//	  icosVertices[i][1] = sinf(phi)*ctheta;
	  icosVertices[i][0] = cos(phi)*ctheta;
	  icosVertices[i][1] = sin(phi)*ctheta;
	  icosVertices[i][2] = -stheta;

	  phi += 2.0*PI/5.0;
  }

  // Upper pentagon (5 vertices)
  phi=0.0;
  for (int i=6; i<11; i++)
  {
	  // push @icosVertices, [cos($phi)*$ctheta, sin($phi)*$ctheta, $stheta];
//	  icosVertices[i][0] = cosf(phi)*ctheta;
//	  icosVertices[i][1] = sinf(phi)*ctheta;
	  icosVertices[i][0] = cos(phi)*ctheta;
	  icosVertices[i][1] = sin(phi)*ctheta;
	  icosVertices[i][2] = stheta;
	  phi += 2.0*PI/5.0;
  }

  // Upper vertex
  icosVertices[NVERT-1][0] = 0.0;
  icosVertices[NVERT-1][1] = 0.0;
  icosVertices[NVERT-1][2] = 1.0;
}

// The "loco" energy calculation stuff...
double loco_energy(float *xyz, int *iseq, int nres, double ****loco, double icosVertices[][3])
{
	double dist;
	double cutoff2=pow(12+CAdist,2);
	double remo[3]; // atomic position in local coordinate frame
	int di,tri;
	double contactscore = 0.0;

	for (int i=0; i<nres; i++)
	{
		for (int j=i+GAP; j<nres; j++)
		{
			dist = DIST2((xyz+i*9+3),(xyz+j*9+3));

//			// If it involves loop...
//			if( (i >= 288 && i <= 299) || (j >= 288 && j <= 299) )

			if (dist <= cutoff2)
			{
				dist=sqrt(dist);
				// Handle "i" residue
//				@remo = getcoor(@{$peptide{$faresnum[$i]}{"N"}}, @{$peptide{$faresnum[$i]}{"CA"}}, @{$peptide{$faresnum[$i]}{"C"}}, @{$peptide{$faresnum[$j]}{"CA"}});
				getcoor(xyz+(i*9), xyz+(i*9)+3, xyz+(i*9)+6, xyz+(j*9)+3, remo); // Compute CA position in local coordinate frame

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Get icoshedron number
				//	my ($d, $tri) = geticoscoor($remo[0], $remo[1], $remo[2], $dist);
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
//				fprintf(stderr,"tri-1= %d\n",tri);
//				fprintf(stderr,"i= %d  ",iseq[i]);
//				fprintf(stderr,"j= %d  ",iseq[j]);
//				fprintf(stderr,"di= %d  ",di);
//				fprintf(stderr,"tri= %d  ",tri);
//				fprintf(stderr,"loco= %f\n", loco[ iseq[i] ][ iseq[j] ][ di ][ tri ]);

//				$contactscore+=$locoscore{substr($fasta, $i, 1).substr($fasta, $j, 1)}{$d}{$tri};
//				fprintf(stderr,"i= %d  j= %d  is= %d  js= %d  di= %d  tri= %d  loco= %f\n", i, j, iseq[i], iseq[j], di, tri, loco[ iseq[i] ][ iseq[j] ][ di ][ tri ]);
				contactscore += loco[ iseq[i] ][ iseq[j] ][ di ][ tri ];
//				fprintf(stderr,"ii= %d  jj= %d  loco= %f\n",i,j,loco[ iseq[i] ][ iseq[j] ][ di ][ tri ]);
//				fprintf(stderr,"ii= %d  jj= %d  loco= %f  coord-i= %f %f %f  coord-j= %f %f %f\n",i,j,loco[ iseq[i] ][ iseq[j] ][ di ][ tri ],xyz[i*9+3],xyz[i*9+4],xyz[i*9+5],xyz[j*9+3],xyz[j*9+4],xyz[j*9+5]);


				// Handle "j" residue
//				@remo = getcoor(@{$peptide{$faresnum[$j]}{"N"}}, @{$peptide{$faresnum[$j]}{"CA"}}, @{$peptide{$faresnum[$j]}{"C"}}, @{$peptide{$faresnum[$i]}{"CA"}});
				getcoor(xyz+(j*9), xyz+(j*9)+3, xyz+(j*9)+6, xyz+(i*9)+3, remo); // Compute CA position in local coordinate frame

				// Get icoshedron number
//				($d, $tri) = geticoscoor($remo[0], $remo[1], $remo[2], $dist);
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
//				fprintf(stderr,"tri-2= %d\n",tri);

//				$contactscore+=$locoscore{substr($fasta, $j, 1).substr($fasta, $i, 1)}{$d}{$tri};
//				fprintf(stderr,"i= %d  j= %d  is= %d  js= %d  di= %d  tri= %d  loco= %f\n", j, i, iseq[j], iseq[i], di, tri, loco[ iseq[j] ][ iseq[i] ][ di ][ tri ]);
				contactscore += loco[ iseq[j] ][ iseq[i] ][ di ][ tri ];
//				double myener =loco[ iseq[j] ][ iseq[i] ][ di ][ tri ];
//				fprintf(stderr,"ii= %d  jj= %d  loco= %f\n",i,j,loco[ iseq[j] ][ iseq[i] ][ di ][ tri ]);
//				fprintf(stderr,"ii= %d  jj= %d  loco= %f  di= %d tri=%d i=%d j=%d seqi=%d seqj=%d d=%f\n",i,j,myener,di,tri,i,j,iseq[i],iseq[j],dist);

			}
		}
	}

	return ( (double)contactscore );
}

// The "loco" energy calculation for protein-protein docking
double loco_energy(float *xyz, int *iseq, int nres, float *xyzL, int *iseqL, int nresL, double ****loco, double icosVertices[][3])
{
	double dist;
	double cutoff2;
	double remo[3]; // atomic position in local coordinate frame
	int di,tri;
	double contactscore = 0.0;
	cutoff2 = pow(12+CAdist,2);

	for (int i=0; i<nres; i++) // screen receptor
	{
		for (int j=0; j<nresL; j++) // screen ligand
		{
			dist = DIST2((xyz+i*9+3),(xyzL+j*9+3));

			if (dist <= cutoff2)
			{
				dist=sqrt(dist);
				// Handle "i" residue
				getcoor(xyz+(i*9), xyz+(i*9)+3, xyz+(i*9)+6, xyzL+(j*9)+3, remo); // Compute CA position in local coordinate frame

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseq[i] ][ iseqL[j] ][ di ][ tri ];


				// Handle "j" residue
				getcoor(xyzL+(j*9), xyzL+(j*9)+3, xyzL+(j*9)+6, xyz+(i*9)+3, remo); // Compute CA position in local coordinate frame

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseqL[j] ][ iseq[i] ][ di ][ tri ];
			}
		}
	}

	return ( (double)contactscore );
}

// The "loco" energy calculation stuff when loop is not present
double loco_energy(float *xyz, int *iseq, int nres, int nresl, int lindex, double ****loco, double icosVertices[][3])
{
	double dist;
	double cutoff2=pow(12+CAdist,2);
	double remo[3]; // atomic position in local coordinate frame
	int di,tri;
	double contactscore = 0.0;

	int myj,myi;
	for (int i=0; i<nres; i++)
	{
		for (int j=i+1; j<nres; j++)
		{

			if(j<lindex)
				myj = j;
			else
				myj = j+nresl;

			if(i<lindex) // MON: this may go above, into between i-for and j-for...
				myi = i;
			else
				myi = i+nresl;

			if(myj-myi >= GAP)
			// if( !( (myj>lindex) && (myj<lindex+nresl) ) || !( (myi>lindex) && (myi<lindex+nresl) ) )

			// if(i > anchorCt-1 || i < anchorCt)
//			if( (i < anchorCt) && (j < anchorCt) &&   )

			// If it does not involve the loop...
			// if( !((i >= lindex && i < lindex+nresl) || (j >= lindex && j < lindex+nresl)) )
//			if( j-i >= GAP || (nresl >= GAP-1 && i < lindex && j >= lindex) || (nresl < GAP-1 && i < lindex && j >= lindex && j-i+nres >=GAP ) )
			{
				dist = DIST2((xyz+i*9+3),(xyz+j*9+3));
				if (dist <= cutoff2)
				{
					dist=sqrt(dist);
					// Handle "i" residue
					//				@remo = getcoor(@{$peptide{$faresnum[$i]}{"N"}}, @{$peptide{$faresnum[$i]}{"CA"}}, @{$peptide{$faresnum[$i]}{"C"}}, @{$peptide{$faresnum[$j]}{"CA"}});
					getcoor(xyz+(i*9), xyz+(i*9)+3, xyz+(i*9)+6, xyz+(j*9)+3, remo); // Compute CA position in local coordinate frame

					// Get CA-CA distance bin
					di = (int) floor(dist-CAdist);
					if(di < 0)
						di = 0;

					// Get icoshedron number
					//	my ($d, $tri) = geticoscoor($remo[0], $remo[1], $remo[2], $dist);
					tri = -1;
					for (int k = 0; k < NTRI; k++)
						if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
						{
							tri = k;
							break;
						}
					//				fprintf(stderr,"tri-1= %d\n",tri);
					//				fprintf(stderr,"i= %d  ",iseq[i]);
					//				fprintf(stderr,"j= %d  ",iseq[j]);
					//				fprintf(stderr,"di= %d  ",di);
					//				fprintf(stderr,"tri= %d  ",tri);
					//				fprintf(stderr,"loco= %f\n", loco[ iseq[i] ][ iseq[j] ][ di ][ tri ]);

					//				$contactscore+=$locoscore{substr($fasta, $i, 1).substr($fasta, $j, 1)}{$d}{$tri};
					//				fprintf(stderr,"i= %d  j= %d  is= %d  js= %d  di= %d  tri= %d  loco= %f\n", i, j, iseq[i], iseq[j], di, tri, loco[ iseq[i] ][ iseq[j] ][ di ][ tri ]);
					contactscore += loco[ iseq[i] ][ iseq[j] ][ di ][ tri ];

					// Handle "j" residue
					//				@remo = getcoor(@{$peptide{$faresnum[$j]}{"N"}}, @{$peptide{$faresnum[$j]}{"CA"}}, @{$peptide{$faresnum[$j]}{"C"}}, @{$peptide{$faresnum[$i]}{"CA"}});
					getcoor(xyz+(j*9), xyz+(j*9)+3, xyz+(j*9)+6, xyz+(i*9)+3, remo); // Compute CA position in local coordinate frame

					// Get icoshedron number
					//				($d, $tri) = geticoscoor($remo[0], $remo[1], $remo[2], $dist);
					tri = -1;
					for (int k = 0; k < NTRI; k++)
						if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
						{
							tri = k;
							break;
						}
					//				fprintf(stderr,"tri-2= %d\n",tri);

					//				$contactscore+=$locoscore{substr($fasta, $j, 1).substr($fasta, $i, 1)}{$d}{$tri};
					//				fprintf(stderr,"i= %d  j= %d  is= %d  js= %d  di= %d  tri= %d  loco= %f\n", j, i, iseq[j], iseq[i], di, tri, loco[ iseq[j] ][ iseq[i] ][ di ][ tri ]);
					contactscore += loco[ iseq[j] ][ iseq[i] ][ di ][ tri ];
				}
			}
		}
	}

	return ( (double)contactscore );
}



// The "loco" energy calculation stuff...
// lindex --> residue index where the loop has been extracted from (in "coord" array), i.e. index of the first right-side residue.
double loco_energy2(float *coordl, int *iseql, int nresl, float *coord, int *iseq, int nres, int lindex, double ****loco, double icosVertices[][3])
{
	double dist;
	double cutoff2=pow(12+CAdist,2);
	double remo[3]; // atomic position in local coordinate frame
	int di,tri,limit;
	double contactscore = 0.0;

	// Loop vs. environment
	for (int i=0; i<nresl; i++) // screen all loop residues
	{
		if( i < GAP )
			limit = lindex-GAP;
		else
			limit = lindex;

		// Left side
		for (int j = 0; j < limit ; j++) // left side environment residues
		{
			dist = DIST2((coordl+i*9+3),(coord+j*9+3));
			if (dist <= cutoff2)
			{
				dist=sqrt(dist);

				// Handle "i" residue
				getcoor(coordl+(i*9), coordl+(i*9)+3, coordl+(i*9)+6, coord+(j*9)+3, remo); // Compute CA position in local coordinate frame

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseql[i] ][ iseq[j] ][ di ][ tri ];

				// Handle "j" residue
				getcoor(coord+(j*9), coord+(j*9)+3, coord+(j*9)+6, coord+(i*9)+3, remo); // Compute CA position in local coordinate frame

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Get icosahedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseq[j] ][ iseql[i] ][ di ][ tri ];
			}
		}

		if( i > nresl-GAP )
			limit = lindex + (nresl-i);
		else
			limit = lindex;

		// Right side
		for (int j = limit; j < nres; j++) // right side environment residues
		{
			dist = DIST2((coordl+i*9+3),(coord+j*9+3));
			if (dist <= cutoff2)
			{
				dist=sqrt(dist);

				// Handle "i" residue (loop vs. environment)
				getcoor(coordl+(i*9), coordl+(i*9)+3, coordl+(i*9)+6, coord+(j*9)+3, remo); // Compute CA position in local coordinate frame

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseql[i] ][ iseq[j] ][ di ][ tri ];

				// Handle "j" residue (environment vs. loop)
				getcoor(coord+(j*9), coord+(j*9)+3, coord+(j*9)+6, coordl+(i*9)+3, remo); // Compute CA position in local coordinate frame

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseq[j] ][ iseql[i] ][ di ][ tri ];

			}
		}
	}


	// Loop vs. loop
	for (int i=0; i<nresl; i++) // screen all loop residues
	{
		for (int j = i+GAP; j < nresl; j++) // screen all loop residues
		{
			dist = DIST2((coordl+i*9+3),(coordl+j*9+3));
			if (dist <= cutoff2)
			{
				dist=sqrt(dist);

				// Handle "i" residue
				getcoor(coordl+(i*9), coordl+(i*9)+3, coordl+(i*9)+6, coordl+(j*9)+3, remo); // Compute CA position in local coordinate frame

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseql[i] ][ iseql[j] ][ di ][ tri ];

				// Handle "j" residue
				getcoor(coordl+(j*9), coordl+(j*9)+3, coordl+(j*9)+6, coordl+(i*9)+3, remo); // Compute CA position in local coordinate frame

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseq[j] ][ iseql[i] ][ di ][ tri ];
			}
		}
	}

	return ( (double)contactscore );
}

// The "loco" energy calculation stuff...
// lindex --> residue index where the loop has been extracted from (in "coord" array), i.e. index of the first right-side residue.
double loco_energy(float *coordl, int *iseql, int nresl, float *coord, int *iseq, int nres, int lindex, double ****loco, double icosVertices[][3])
{
	double dist;
	double cutoff2=pow(12+CAdist,2);
	double remo[3]; // atomic position in local coordinate frame
	int di,tri;
	double contactscore = 0.0;

	int myj,myi,swap;

	// Left side
	for (int j = 0; j < nres; j++) // environment residues
	{
		// Environment vs. loop
		for (int i = 0; i < nresl; i++) // screen all loop residues
		{
			myi = lindex+i;

			if(j<lindex)
				myj = j;
			else
				myj = j+nresl;

			if(myi > myj)
			{
				swap = myi;
				myi = myj;
				myj = swap;
			}

			if (myj-myi >= GAP)

//			if( (lindex+i) - j >= GAP || (j+nresl)-(lindex+i) >= GAP )
			{
				dist = DIST2((coordl+i*9+3),(coord+j*9+3)); // CA(i)-CA(j) distance
				if (dist < cutoff2)
				{
					dist=sqrt(dist);

					// Get CA-CA distance bin
					di = (int) floor(dist-CAdist);
					if(di < 0)
						di = 0;

					// Handle "i" residue
					getcoor(coordl+(i*9), coordl+(i*9)+3, coordl+(i*9)+6, coord+(j*9)+3, remo); // Compute CA position in local coordinate frame

					// Get icoshedron number
					tri = -1;
					for (int k = 0; k < NTRI; k++)
						if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
						{
							tri = k;
							break;
						}
					contactscore += loco[ iseql[i] ][ iseq[j] ][ di ][ tri ];

//					myener = loco[ iseql[i] ][ iseq[j] ][ di ][ tri ];
//					fprintf(stderr,"ii= %d  jj= %d  loco= %f  coordl= %f %f %f  coord= %f %f %f\n",myi,myj,myener,coordl[i*9+3],coordl[i*9+4],coordl[i*9+5],coord[j*9+3],coord[j*9+4],coord[j*9+5]);

					// Handle "j" residue
					getcoor(coord+(j*9), coord+(j*9)+3, coord+(j*9)+6, coordl+(i*9)+3, remo); // Compute CA position in local coordinate frame

					// Get icoshedron number
					tri = -1;
					for (int k = 0; k < NTRI; k++)
						if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
						{
							tri = k;
							break;
						}
					contactscore += loco[ iseq[j] ][ iseql[i] ][ di ][ tri ];

//					myener = loco[ iseq[j] ][ iseql[i] ][ di ][ tri ];
//					fprintf(stderr,"ii= %d  jj= %d  loco= %f\n",myi,myj,myener);
//					fprintf(stderr,"ii= %d  jj= %d  loco= %f  di= %d tri=%d i=%d j=%d seqi=%d seqj=%d d=%f\n",myi,myj,myener,di,tri,i,j,iseql[i],iseq[j],dist);
				}
			}
		}
	}

	// Loop vs. loop
	for (int i=0; i<nresl; i++) // screen all loop residues
	{
		for (int j = i+GAP; j < nresl; j++) // screen all loop residues
		{
			dist = DIST2((coordl+i*9+3),(coordl+j*9+3));
			if (dist <= cutoff2)
			{
				dist=sqrt(dist);

				// Get CA-CA distance bin
				di = (int) floor(dist-CAdist);
				if(di < 0)
					di = 0;

				// Handle "i" residue
				getcoor(coordl+(i*9), coordl+(i*9)+3, coordl+(i*9)+6, coordl+(j*9)+3, remo); // Compute CA position in local coordinate frame

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseql[i] ][ iseql[j] ][ di ][ tri ];

//				myener = loco[ iseql[i] ][ iseql[j] ][ di ][ tri ];
//				fprintf(stderr,"ii= %d  jj= %d  loco= %f\n",lindex+i,lindex+j,myener);

				// Handle "j" residue
				getcoor(coordl+(j*9), coordl+(j*9)+3, coordl+(j*9)+6, coordl+(i*9)+3, remo); // Compute CA position in local coordinate frame

				// Get icoshedron number
				tri = -1;
				for (int k = 0; k < NTRI; k++)
					if( rayIntersectsTriangle(remo, icosVertices[icosTri[k][0]], icosVertices[icosTri[k][1]], icosVertices[icosTri[k][2]]) )
					{
						tri = k;
						break;
					}
				contactscore += loco[ iseql[j] ][ iseql[i] ][ di ][ tri ];

//				myener = loco[ iseql[j] ][ iseql[i] ][ di ][ tri ];
//				fprintf(stderr,"ii= %d  jj= %d  loco= %f\n",lindex+i,lindex+j,myener);

			}
		}
	}
// exit(0);
	return ( (double)contactscore );
}




// Compute CA position in local coordinate frame
void getcoor(float *Nref, float *CAref, float *Cref, float *contactCAref, double *cen)
{
	double CANdist = DIST(CAref,Nref);
	double CACdist = DIST(CAref,Cref);
	double N[3],CA[3],C[3];
	CA[0] = 0.0;
	CA[1] = 0.0;
	CA[2] = 0.0;
	N[0] = CANdist;
	N[1] = 0.0;
	N[2] = 0.0;
	double CACvector[3] = { CAref[0]-Cref[0], CAref[1]-Cref[1], CAref[2]-Cref[2] };
	double CANvector[3] = { CAref[0]-Nref[0], CAref[1]-Nref[1], CAref[2]-Nref[2] };
	double NCACangle = acosf( (CACvector[0]*CANvector[0]+CACvector[1]*CANvector[1]+CACvector[2]*CANvector[2])/(CANdist*CACdist) );
	C[0] = CACdist*cosf(NCACangle);
	C[1] = CACdist*sinf(NCACangle);
//	double NCACangle = acos( (CACvector[0]*CANvector[0]+CACvector[1]*CANvector[1]+CACvector[2]*CANvector[2])/(CANdist*CACdist) );
//	C[0] = CACdist*cos(NCACangle);
//	C[1] = CACdist*sin(NCACangle);
	C[2] = 0.0;

	double d[3];
	d[0] = DIST2(contactCAref,Nref);
	d[1] = DIST2(contactCAref,CAref);
	d[2] = DIST2(contactCAref,Cref);

//fprintf(stderr,"getcoor> d= %f %f %f\n",d[0],d[1],d[2]);

//  # set up a, b, c arrays
//  my (@a, @b, @c);
//  push(@a, $N[2]);
//  push(@b, $N[1]);
//  push(@c, $N[0]);
//  push(@a, $CA[2]);
//  push(@b, $CA[1]);
//  push(@c, $CA[0]);
//  push(@a, $C[2]);
//  push(@b, $C[1]);
//  push(@c, $C[0]);
//

//	my $D2=(($d[0]-$d[1])-($b[0]*$b[0]-$b[1]*$b[1])-($a[0]*$a[0]-$a[1]*$a[1])-($c[0]*$c[0]-$c[1]*$c[1]))/2.0;
	double D2 = ( (d[0]-d[1]) - (N[1]*N[1]-CA[1]*CA[1]) - (N[2]*N[2]-CA[2]*CA[2]) - (N[0]*N[0]-CA[0]*CA[0]) ) / 2.0;

//fprintf(stderr,"getcoor> D2= %f\n",D2);

//	my $D3=(($d[0]-$d[2])-($b[0]*$b[0]-$b[2]*$b[2])-($a[0]*$a[0]-$a[2]*$a[2])-($c[0]*$c[0]-$c[2]*$c[2]))/2.0;
	double D3 = ( (d[0]-d[2]) - (N[1]*N[1]-C[1]*C[1]) - (N[2]*N[2]-C[2]*C[2]) - (N[0]*N[0]-C[0]*C[0]) ) / 2.0;

//fprintf(stderr,"getcoor> D3= %f\n",D3);

//  my $A=(($a[2]-$a[0])*($c[1]-$c[0])-($a[1]-$a[0])*($c[2]-$c[0])) / (($b[1]-$b[0])*($c[2]-$c[0])-($b[2]-$b[0])*($c[1]-$c[0]));
	double A = ((C[2]-N[2])*(CA[0]-N[0])-(CA[2]-N[2])*(C[0]-N[0])) / ((CA[1]-N[1])*(C[0]-N[0])-(C[1]-N[1])*(CA[0]-N[0]));
//  my $B=($D2*($c[2]-$c[0])-$D3*($c[1]-$c[0])) / (($b[1]-$b[0])*($c[2]-$c[0])-($b[2]-$b[0])*($c[1]-$c[0]));
	double B = (D2*(C[0]-N[0]) - D3*(CA[0]-N[0]) ) / ( (CA[1]-N[1])*(C[0]-N[0]) - (C[1]-N[1])*(CA[0]-N[0]) );

//  my $AA=($a[0]-$a[1]+$A*($b[0]-$b[1]))/($c[1]-$c[0]);
	double AA = (N[2]-CA[2]+A*(N[1]-CA[1])) / (CA[0]-N[0]);

//  my $BB=($D2-$B*($b[1]-$b[0]))/($c[1]-$c[0]);
	double BB = (D2-B*(CA[1]-N[1])) / (CA[0]-N[0]);

//  my $AAA=1.0+$A*$A+$AA*$AA;
	double AAA = 1.0 + A*A + AA*AA;

//  my $BBB=2*$A*($B-$b[0])+2*$AA*($BB-$c[0])-2*$a[0];
	double BBB = 2*A*(B*N[1]) + 2*AA*(BB-N[0]) - 2*N[2];

//  my $CCC=$a[0]*$a[0]+($B-$b[0])*($B-$b[0])+($BB-$c[0])*($BB-$c[0])-$d[0];
	double CCC = N[2]*N[2]+(B-N[1])*(B-N[1])+(BB-N[0])*(BB-N[0])-d[0];

//  my @x = ([0, 0, 0], [0, 0, 0]);
	double x[2][3];
//
//  if ($BBB*$BBB-4*$AAA*$CCC>=0)
	if(BBB*BBB-4*AAA*CCC >= 0)
		x[0][2] = (-BBB+sqrtf(BBB*BBB-4*AAA*CCC))/2.0/AAA; // Mon: 2nd order equation (+ solution): Ax^2 + Bx + C = 0
//		x[0][2] = (-BBB+sqrt(BBB*BBB-4*AAA*CCC))/2.0/AAA; // Mon: 2nd order equation (+ solution): Ax^2 + Bx + C = 0
	else
		x[0][2] = (-BBB)/2.0/AAA;
	x[0][1]=A*x[0][2]+B;
	x[0][0]=AA*x[0][2]+BB;

	double NCA_org[3];    // vector N->CA
	double CCA_org[3];    // vector C->CA
	double contactCA_org[3];  // vector contactCA->CA
	double normal_org[3];	// normal of the plane N-CA-C

	for (int i=0; i<3; i++)
	{
		NCA_org[i]=Nref[i]-CAref[i];
		CCA_org[i]=Cref[i]-CAref[i];
		contactCA_org[i]=contactCAref[i]-CAref[i];
	}

	normal_org[0]=NCA_org[1]*CCA_org[2]-NCA_org[2]*CCA_org[1];
	normal_org[1]=NCA_org[2]*CCA_org[0]-NCA_org[0]*CCA_org[2];
	normal_org[2]=NCA_org[0]*CCA_org[1]-NCA_org[1]*CCA_org[0];

	// the trick here is we don't need to normalize it
	double angle_org = normal_org[0]*contactCA_org[0]+normal_org[1]*contactCA_org[1]+normal_org[2]*contactCA_org[2];

	double NCA_new[3]; // vector N->CA
	double CCA_new[3]; //vector C->CA
	double contactCA_new[3]; // vector contactCA->CA
	double normal_new[3]; // normal of the plane N-CA-C

	for (int i=0; i<3; i++)
	{
		NCA_new[i] = N[i]-CA[i];
		CCA_new[i] = C[i]-CA[i];
		contactCA_new[i] = x[0][i]-CA[i];
	}

	normal_new[0] = NCA_new[1]*CCA_new[2] - NCA_new[2]*CCA_new[1];
	normal_new[1]=NCA_new[2]*CCA_new[0] - NCA_new[0]*CCA_new[2];
	normal_new[2]=NCA_new[0]*CCA_new[1] - NCA_new[1]*CCA_new[0];

	double angle_new = normal_new[0]*contactCA_new[0]+normal_new[1]*contactCA_new[1]+normal_new[2]*contactCA_new[2];

	if (angle_new*angle_org>0)
	{
		for (int i=0;i<3;i++)
			cen[i]=x[0][i];
	}
	else
	{
		if (BBB*BBB-4*AAA*CCC >= 0)
			x[1][2]=(-BBB-sqrtf(BBB*BBB-4*AAA*CCC))/2.0/AAA;
//			x[1][2]=(-BBB-sqrt(BBB*BBB-4*AAA*CCC))/2.0/AAA;
		else
			x[1][2]=(-BBB)/2.0/AAA;

		x[1][1] = A*x[1][2] + B;
		x[1][0] = AA*x[1][2] + BB;

		for (int i=0; i<3; i++)
			cen[i] = x[1][i];
	}

//	fprintf(stderr,"getcoor> cen= %f %f %f\n",cen[0],cen[1],cen[2]);
}

// Determine if a ray going through a triangle.
// 	( from: http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/ )
bool rayIntersectsTriangle(double *dref, double *v0, double *v1, double *v2)
{
	double e1[3], e2[3], h[3], s[3], q[3];
	double a, f, u, v, t;

	// @e1 = vector(@{$v1ref}, @{$v0ref});
	//  e1[0] = v1[0]-v0[0];
	//  e1[1] = v1[1]-v0[1];
	//  e1[2] = v1[2]-v0[2];
	VECTOR(v1,v0,e1); // vector v1-v0

//fprintf(stderr,"rayIntersectsTriangle> e1= %f %f %f\n",e1[0],e1[1],e1[2]);

	// @e2 = vector(@{$v2ref}, @{$v0ref});
	//  e2[0] = v2[0]-v0[0];
	//  e2[1] = v2[1]-v0[1];
	//  e2[2] = v2[2]-v0[2];
	VECTOR(v2,v0,e2); // vector v2-v0

//fprintf(stderr,"rayIntersectsTriangle> e2= %f %f %f\n",e2[0],e2[1],e2[2]);

	// @h = cross(@{$dref}, @e2);
	CROSS(dref,e2,h); // cross product

//fprintf(stderr,"rayIntersectsTriangle> dref= %f %f %f\n",dref[0],dref[1],dref[2]);
//fprintf(stderr,"rayIntersectsTriangle> h= %f %f %f\n",h[0],h[1],h[2]);

	// $a = dot(@e1, @h);
	a = DOT(e1,h);

//fprintf(stderr,"rayIntersectsTriangle> a= %f\n",a);

	if (a > -0.00001 && a < 0.00001)
		return(false);

	// @s = vector(@{$pref}, @{$v0ref}); // pref={0,0,0}
	REVERSE(v0,s);

	f = 1/a;
	u = f * DOT(s,h);

//fprintf(stderr,"rayIntersectsTriangle> u= %f\n",u);

	if (u < 0.0 || u > 1.0)
		return(false);

	// @q=cross(@s, @e1);
	CROSS(s,e1,q);

	// $v = $f * dot(@{$dref}, @q);
	v = f * DOT(dref, q);

//fprintf(stderr,"rayIntersectsTriangle> v= %f\n",v);

	if (v < 0.0 || u+v > 1.0)
		return(false);

	// at this stage we can compute t to find out where the intersection point is on the line
	// $t = $f * dot(@e2, @q);
	t = f * DOT(e2, q);

//fprintf(stderr,"rayIntersectsTriangle> t= %f\n",t);

	if (t > 0)
		return(true);
	else
		return(false);
}

// Convert "co->el" 2D-matrix (Pieter's) into linear array (coordMatrix-like)
//  el --> 2D-matrix pointer
//  coord --> linear array
//  natoms --> number of atoms to copy coordinates
void co2coord(double **el, float *coord, int natoms)
{
	for(int i=0; i<natoms; i++)
	{
		coord[3*i]   = el[i][0];
		coord[3*i+1] = el[i][1];
		coord[3*i+2] = el[i][2];
	}
}

// Convert a 1-letter code sequence into integer LOCO code sequence
//  seq --> 1-letter code sequence array
//  p_iseq --> Pointer to the array of integers with LOCO code sequence (if NULL, it allocates memory)
//  nres --> Number of residues
void seq2iseq(char *seq, int **p_iseq, int nres)
{
	// Initialize indices sequence (array of indices)
	if(*p_iseq == NULL)
	{
		*p_iseq = (int *) malloc(sizeof(int)*nres);
//		if( !(*p_iseq = (int *) malloc(sizeof(int)*nres)) )
//		{
//			fprintf(stderr,"seq2iseq> Sorry, memory allocation error (nres= %d)\n",nres);
//			exit(1);
//		}
	}
	for(int i=0; i<nres; i++)
	{
		// fprintf(stderr,"seq[%d] = %c\n", i,seq[i]);
		(*p_iseq)[i] = aa2index(seq[i],aa); // Index for a given AA
	}

}
