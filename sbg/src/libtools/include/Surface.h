// FastSqrt.h
//

// taken from
// T. Can, C.-I. Chen, Y.-F. Wang, "Efficient molecular surface generation using level-set methods," Journal of Molecular Graphics and Modelling (JMGM)

#ifndef _SURFACE_H
#define _SURFACE_H


typedef struct S_prot
{
        float *xpoints;
        float *ypoints;
        float *zpoints;
        float *rpoints;
        int npoints;
} *S_Molecule;


typedef struct p
{
        short x;
        short y;
        short z;
} S_PPoint;

typedef struct gp
{
        S_PPoint point;
        char phi;
        int from;
        float dist;
} S_GridPoint;

typedef struct g
{
        short N;
        S_GridPoint*** matrix;
        short stepSize;
} *S_Grid;


typedef struct nb
{
        short x;
        short y;
        short z;
        struct nb* next;
} nbNode;

typedef struct
{
	float fX;
	float fY;
	float fZ;
} GLvector;




S_Grid createGrid(short i);
S_Grid copyGrid(S_Grid g2);
void destroyGrid(S_Grid g);
/**
* PR: Number of voxels
*/
int findProbesMol(S_Grid g, float PR);
void shrink(S_Grid g, float PR);
int fastMarching(S_Grid g, bool inner);
void expand(S_Grid g, float PR);
void marchingCube(S_Grid grid);

#endif

