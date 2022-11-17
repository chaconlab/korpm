//============================================================================
// Name        : korpm_opt.cpp
// Author      : Pablo and Mon
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

//  KORPM OPTIMIZATION
//  icpc  -L/usr/local/lib test.cpp -o opt -lnlopt

#define VERSION "v1.16"

#include <math.h>
#include <nlopt.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmdl/CmdLine.h>
#include <random>

// #define NVAR 41
#define NVAR 401
int nvar = 0;

//#define KN 10   // Kn fold
//#define NRUNS 10
//#define KNNRUNS KN*NRUNS
#define NFACTOR 0.8
#define MAXCASES 10000
int BINSAUC=1000;

// Required by "quicksort"
#define swap(a,b,t) ((t)=(a),(a)=(b),(b)=(t))

// AA list    0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20
char aa[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','p'};
//char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
char map[]= { 0,  1,  2,  2,  3,  4,  5,  6,  5,  6,  7,  8,  9,  8, 10, 8, 8, 6, 9, 3, 11};  // HK + FY + VIL + NQST + DE 0.525


int hyb []= { 0,  0,  1,  1,  0,  2,  2,  0,  1,  0,  0,   1,  2,  1,  1,  2,  2,  0,  0,  2, 20};   // A, C, I, L, M, F, W, V  Hydrophilic R, N, D, Q, E, K  neutral  G, H, P, S, T, Y
int pnn []= { 2,  2,  0,  0,  2,  2,  1,  2,  1,  2,  2,   2,  2,  2,  1,  2,  2,  2,  2,  2, 20};   //	Positive charged (R, H, K) negative charged (D, E).
float vol []= {88.6, 108.5, 111.1, 138.4, 189.9, 193.6, 60.1, 153.2, 166.7, 168.6, 166.7, 162.9, 114.1, 112.7, 143.8, 173.4, 89.0, 116.1, 140.0, 227.8, 193.6};
float dol []= {1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9	, 3.8	, 1.9, 	-3.5, -1.6, -3.5, -4.5, -0.8,-0.7	, 4.2, -0.9	, -1.3	, -1.6};

float *arec; // RECall array (TPR)
float *appv; // PREcision array (PPV)
float *afpr; // False Positives Rate array

float **Arec; // RECall array (TPR)
float **Appv; // PREcision array (PPV)
float **Afpr; // False Positives Rate array

float minP=10000;
float maxP=-1000;




//char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  7,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19};   // I L 0.527
// char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  11, 12, 13, 14, 15, 15, 16, 17, 18, 19}; // S T   0.527
// char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  11, 12, 11, 13, 14, 15, 16, 17, 18, 19}; // N Q  OK
//   char map[]= { 0,  1,  2,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19};  // D E 0.528
//   char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  11, 12, 13,  8, 14, 15, 16, 17, 18, 19}; // K R  0.523 *
//    char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  11, 12, 13, 14, 15, 16, 17, 18,  4, 19}; // F Y 0.527
//    char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  11, 12, 13, 14, 15, 16, 17, 4,  18, 19}; // F W  0.526
//     char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  11, 12, 13, 14, 15, 16, 17, 18, 18, 19}; //  W Y  0.519
// char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  7,  9, 10, 11, 10, 12, 13, 14, 15, 16, 17, 18};  // IL + NQ 0.528
// char map[]= { 0,  1,  2,  3,  4,  5,  6,  7,  8,  7,  9, 10, 11, 10, 12, 13, 13, 14, 15, 16, 17};  // IL + NQ + ST 0.528
// char map[]= { 0,  1,  2,  2,  3,  4,  5,  6,  7,  6,  8,  9, 10,  9, 11, 12, 12, 13, 14, 15, 16};  // IL + NQ + ST + DE 0.527
// char map[]= { 0,  1,  2,  2,  3,  4,  5,  6,  7,  6,  8,  9, 10,  9, 11, 12, 12, 13, 14, 3, 15};  // FY + IL + NQ + ST + DE 0.525




int Score[20][20] = {
		{  9, -2, -2, -2, -1, -1, -1, -1, -3, -3,  0, -1, -1, -3, -3, -3, -4, -3, -3, -3 },
		{ -2,  6,  3,  1,  0,  0,  0, -1, -3, -4, -2, -2, -2, -3, -1, -3, -3, -3, -3, -3 },
		{ -2,  3,  7,  2, -1, -1, -1, -1, -3, -3, -2, -2, -2, -2,  2, -1, -2, -3, -2, -2 },
		{ -2,  1,  2, 11, -1, -2, -3, -3, -2, -4, -3, -2, -3, -4, -2, -2, -3, -4, -3, -3 },
		{ -1,  0, -1, -1,  5,  2,  1,  1, -3, -2, -1, -1, -1, -2, -2,  0, -2, -3, -1, -1 },
		{ -1,  0, -1, -2,  2,  4,  2,  1, -4, -3, -1, -1, -2, -3, -3, -2, -3, -4, -2, -2 },
		{ -1,  0, -1, -3,  1,  2,  4,  3, -4, -3, -1, -1, -2, -3, -3, -3, -3, -3, -3, -3 },
		{ -1, -1, -1, -3,  1,  1,  3,  4, -3, -2,  0,  0, -2, -3, -3, -2, -2, -3, -3, -2 },
		{ -3, -3, -3, -2, -3, -4, -4, -3,  6, -2,  0, -2,  0,  0, -2, -2, -2, -1, -2, -2 },
		{ -3, -4, -3, -4, -2, -3, -3, -2, -2,  7, -1, -1, -1, -2, -2, -1, -1, -1, -2, -1 },
		{  0, -2, -2, -3, -1, -1, -1,  0,  0, -1,  4,  0,  1, -2, -2, -1, -1, -2, -1, -1 },
		{ -1, -2, -2, -2, -1, -1, -1,  0, -2, -1,  0,  5,  1,  0, -2, -1, -1, -1, -1, -1 },
		{ -1, -2, -2, -3, -1, -2, -2, -2,  0, -1,  1,  1,  4,  1, -1,  0,  0,  0, -1,  0 },
		{ -3, -3, -2, -4, -2, -3, -3, -3,  0, -2, -2,  0,  1,  6,  1,  0,  0,  1,  0,  0 },
		{ -3, -1,  2, -2, -2, -3, -3, -3, -2, -2, -2, -2, -1,  1,  8,  0,  0, -1,  0, -1 },
		{ -3, -3, -1, -2,  0, -2, -3, -2, -2, -1, -1, -1,  0,  0,  0,  5,  2,  0,  1,  1 },
		{ -4, -3, -2, -3, -2, -3, -3, -2, -2, -1, -1, -1,  0,  0,  0,  2,  5,  2,  0,  1 },
		{ -3, -3, -3, -4, -3, -4, -3, -3, -1, -1, -2, -1,  0,  1, -1,  0,  2,  6, -2, -1 },
		{ -3, -3, -2, -3, -1, -2, -3, -3, -2, -2, -1, -1, -1,  0,  0,  1,  0, -2,  5,  2 },
		{ -3, -3, -2, -3, -1, -2, -3, -2, -2, -1, -1, -1,  0,  0, -1,  1,  1, -1,  2,  5 }
};




typedef struct {
	int nexp;
	double aveE;
	double aveC;
	double sumay;
	float fact;
	float dexp[MAXCASES][40]; // Wild-type sums
	float dcom[MAXCASES][40]; // Mutant sums
	float exp[MAXCASES]; // Experimental ddG
	float com[MAXCASES]; // Predicted ddG
	int aa1[MAXCASES]; // Wild-type aminoacid index
	int aa2[MAXCASES]; // Mutant aminoacid index
	int index[MAXCASES]; // Index for shuffle
	int index2[MAXCASES]; // Index for shuffle
	int Hclass[MAXCASES]; // homology class type;
	float rsa[MAXCASES]; // Mutation RSA
	int nrange; // 0-nrange training nrange-nexp testing
	float EramaWT[MAXCASES][3]; // Rama Wild-Type energy partial sums
	float EramaMut[MAXCASES][3]; // Rama Mutant energy partial sums
	int iaaWT[MAXCASES][3]; // Rama Wild-Type residue type index
	int iaaMut[MAXCASES][3]; // Rama Mutant residue type index
} fdata;


double PCC(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;
	double dumpx, ave=0, dumpy, sumax=0, sumay=0, suma=0;
	// int nrange;

	// nrange=d->nexp*d->fact;

	ave=0;
	for(int ii=0; ii < d->nrange; ii++)
	{
		int i=d->index[ii];

		double sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			sumaa+=x[d->aa1[i]]*x[j]*(d->dexp[i][j*2]+d->dexp[i][2*j+1]);
			sumab+=x[d->aa2[i]]*x[j]*(d->dcom[i][j*2]+d->dcom[i][2*j+1]);
		}
		dumpx=sumab-sumaa;

		dumpx/=100;
		ave+=dumpx;
		d->com[i]=dumpx;
	}

	ave/=(d->nrange);

	for(int ii=0;ii<d->nrange;ii++)
	{
		int i=d->index[ii];
		dumpx= d->com[i]-ave;
		dumpy= d->exp[i]-d->aveE;
		suma+=dumpx*dumpy;
		sumax+=dumpx*dumpx;
		sumay+=dumpy*dumpy;
	}

	ave = suma/(sqrt(sumax)*sqrt(sumay));


	//     for(int j=0;j<20;j++)
	//       fprintf(stdout," %c %.2f ", aa[j], x[j]);
	//      fprintf(stdout," suma %f %f\n", ave, d->fact);

	return ave;
}

double PPV(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0, ppv=0, sen=0, mcc=0;
	//int nrange;

	//nrange = d->nexp * d->fact;

	ave=0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];

		double sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			sumaa += x[ d->aa1[i] ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type
			sumab += x[ d->aa2[i] ] * x[j] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated
		}
		dumpx = sumab - sumaa; // ddG (check sign)

		d->com[i]=dumpx/100;
	}

	float TPddG=0;
	float FNddG=0;
	float TNddG=0;
	float FPddG=0;
	float thr = 0.0;

	suma = 0.0; float dump=0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];

		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)

		// suma += dumpx * dumpx;
		suma += fabsf(dumpx);
		// Huber loss
		//float delta=2.0;
		//dump=fabsf(dumpx);
		//if (dump>delta)  suma+=0.5*dumpx* dumpx;
		//else suma += delta*dump-0.5*delta*delta;

		if(d->exp[i] >= thr) // Stabilizing mutation, ddG<0 (Positive)
		{
			if(d->com[i] >= thr )
				TPddG++; // Count True Positives

			if(d->com[i] < thr)
				FNddG++; // Count False Negatives
		}
		else // Destabilizing mutation, ddG > 0 (Negative)
		{
			if(d->com[i] < thr)
				TNddG++; // Count True Negatives

			if(d->com[i] >= thr)
				FPddG++; // Count False Positives
		}


		if(d->exp[i] >= 0) // Stabilizing mutation, ddG<0 (Positive)
		{
			if(d->com[i] >= 0 )
			{
				TPddG++; // Count True Positives
			}

			if(d->com[i] < 0)
			{
				FNddG++; // Count False Negatives
			}
		}
		else // Destabilizing mutation, ddG > 0 (Negative)
		{
			if(d->com[i] < 0)
			{
				TNddG++; // Count True Negatives
			}

			if(d->com[i] >= 0)
			{
				FPddG++; // Count False Positives
			}
		}



	}


	// rmse = sqrt(suma/nrange);
	rmse = (suma/d->nrange);


	if( TPddG + FPddG == 0 )
		ppv=0;
	else
		ppv = TPddG/(TPddG+FPddG);

	if( TPddG + FNddG == 0 )
		sen=0;
	else
		sen = TPddG/(TPddG+FNddG);
	//  fprintf(stdout," %f %f\n", TPddG, ppv);
	/*       for(int j=0;j<20;j++)
  fprintf(stdout," %c %.2f ", aa[j], x[j]);
  fprintf(stdout," suma %f\n", rmse);
	 */
	//	if( (TPddG + FPddG)*(TPddG + FNddG)*(TNddG + FPddG)*(TNddG + FNddG) != 0 )
	//					mcc = (float) (TPddG * TNddG - FPddG * FNddG) / ( sqrt( (float) (TPddG + FPddG)*(TPddG + FNddG)*(TNddG + FPddG)*(TNddG + FNddG)) ); // Matthews correlation coefficient
	//				else
	//					mcc = 0.0;
	//fprintf(stdout," suma %f %f\n", rmse, mcc);

	// return (rmse+(1-ppv) );
	return (rmse);
}

double TEST0(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx,  sumax=0, suma=0;
	double sumaa=0, sumab=0;
	int i;




	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];

		sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			sumaa +=  x[j] ; // Wild-Type
			//sumab += x[ d->aa2[i] ] * x[j] ; // Mutated
		}
		//dumpx = sumab - sumaa; // ddG (check sign)
		d->com[i] =sumaa*(x[ d->aa2[i] ] - x[ d->aa1[i] ]) ; // ddG (check sign)

		//d->com[i]=dumpx; ///100;


		suma += fabsf(( d->com[i] - d->exp[i] ));
	}


	return (suma/d->nrange);
}



double TEST0_dic(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx,  sumax=0, suma=0;
	double sumaa=0, sumab=0;
	int i;
	int wt, mut, aa;
	//int nrange;

	//nrange = d->nexp * d->fact;

	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];

		if (nvar != 1 )
		{

			wt = map [d->aa1[i]];
			mut = map [d->aa2[i]];
			d->com[i] = x[ mut ] - x[ wt ];

			//			sumaa=0, sumab=0;
			//			for(int j=0;j<20;j++)
			//			{
			//				aa = map[j];
			//				sumaa += x[ wt ] * x[aa] ; // Wild-Type
			//				sumab += x[ mut ] * x[aa] ; // Mutated
			//
			//			}
			//			dumpx = sumab - sumaa; // ddG (check sign)
			//			d->com[i]=dumpx;


		}
		else
		{
			d->com[i]=x[0];
		}
	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		suma += fabsf(( d->com[i] - d->exp[i] ));
	}

	// fprintf(stdout," %10.6f %10.6f\n", suma/d->nrange, suma);

	return (suma/d->nrange);
}

double TEST0H_dic(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx,  sumax=0, suma=0;
	double sumaa=0, sumab=0;
	int i;
	int wt, mut, aa;
	//int nrange;

	//nrange = d->nexp * d->fact;

	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];

		if (nvar != 1 )
		{

			wt = map [d->aa1[i]];
			mut = map [d->aa2[i]];
			d->com[i] = x[ mut ] - x[ wt ];

			//			sumaa=0, sumab=0;
			//			for(int j=0;j<20;j++)
			//			{
			//				aa = map[j];
			//				sumaa += x[ wt ] * x[aa] ; // Wild-Type
			//				sumab += x[ mut ] * x[aa] ; // Mutated
			//
			//			}
			//			dumpx = sumab - sumaa; // ddG (check sign)
			//			d->com[i]=dumpx;


		}
		else
		{
			d->com[i]=x[0];
		}
	}



	suma = 0.0; float dump=0; float delta=3.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		dump=fabsf(dumpx);
		if (dump<=delta)  suma+=0.5*dumpx* dumpx;
		else suma += delta*dump-0.5*delta*delta;
	}


	// fprintf(stdout," %10.6f %10.6f\n", suma/d->nrange, suma);

	return (suma/d->nrange);
}



double MAE_dic(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx,  sumax=0, suma=0;
	double sumaa=0, sumab=0;
	int i;
	int wt, mut, aa;
	//int nrange;

	//nrange = d->nexp * d->fact;

	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];

		wt = map [d->aa1[i]];
		mut = map [d->aa2[i]];


		sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{

			aa = map[j];
			sumaa += x[ wt ] * x[aa] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type
			sumab += x[ mut ] * x[aa]  * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated

		}
		dumpx = sumab - sumaa; // ddG (check sign)
		d->com[i]=dumpx/100;
	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		suma += fabsf(( d->com[i] - d->exp[i] ));
		//if (d->exp[i] > 0 )  suma += 3.0*fabsf(( d->com[i] - d->exp[i] ));
		//    else 		suma += fabsf(( d->com[i] - d->exp[i] ));
	}


	return (suma/d->nrange);
}


double AUC_dic(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx,  sumax=0, suma=0;
	double sumaa=0, sumab=0;
	int i;
	int wt, mut, aa;
	//int nrange;

	//nrange = d->nexp * d->fact;

	float minP=10000;
	float maxP=-1000;



	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];

		wt = map [d->aa1[i]];
		mut = map [d->aa2[i]];

		sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			aa = map[j];
			sumaa += x[ wt ] * x[aa] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type
			sumab += x[ mut ] * x[aa]  * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated

		}
		dumpx = sumab - sumaa; // ddG (check sign)
		d->com[i]=dumpx/100;

		//		if(d->exp[i] < minP) minP=d->exp[i];
		//		if(d->exp[i] > maxP) maxP=d->exp[i];



	}




	// AUC
	float  TPddG = 0;
	float  FPddG = 0;
	float  FNddG = 0;
	float  TNddG = 0;
	float deltau=  (maxP - minP)/BINSAUC;
	for(int au=0; au<BINSAUC; au++)
	{
		TPddG = 0;
		FPddG = 0;
		FNddG = 0;
		TNddG = 0;
		float thr = minP + au * deltau;
		for(int ii=0;ii<d->nrange;ii++)
		{
			i=d->index[ii];

			// Contingency table (confusion matrix)
			if(d->exp[i] >= 0.0) // Stabilizing mutation, ddG >=0 (Positive)
			{
				if(d->com[i] >= thr )
					TPddG++; // Count True Positives
				else
					FNddG++; // Count False Negatives

			}
			else // Destabilizing mutation, ddG < 0 (Negative)
			{
				if(d->com[i] < thr)
					TNddG++; // Count True Negatives
				else
					FPddG++; // Count False Positives
			}
		}

		if( TPddG + FNddG > 0 )
			arec[au]=(float) TPddG / ((float) TPddG + (float) FNddG);  // RECall = Sensitivity (SN) = True Positive Rate (TPR)
		else
			arec[au]=0;
		if( TPddG + FPddG > 0 )
			appv[au]=(float) TPddG / ((float) TPddG + (float) FPddG); // Positive Predictive Value = PRECision
		else
			appv[au]=0;
		if( TNddG + FPddG > 0 )
			afpr[au]=(float) FPddG / (float) (TNddG + FPddG); // False Positive Rate (FPR)
		else
			afpr[au]=0;
	}
	float AUC_ROC = 0.0; // Area Under Curve for ROC
	float AUC_PRC = 0.0; // Area Under Curve for PRC

	// for(int au= -minP/deltau-200; au<-minP/deltau+200; au++)

	for(int au=0; au<(BINSAUC-1); au++)

	{
		AUC_ROC += (afpr[au] - afpr[au+1]) * (arec[au] + arec[au+1]) / 2;
		AUC_PRC += (arec[au] - arec[au+1]) * (appv[au] + appv[au+1]) / 2;
	}
	/*
	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		suma += fabsf(( d->com[i] - d->exp[i] ));
	}
	 */
	suma = 0.0; float dump=0; float delta=2.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		dump=fabsf(dumpx);
		if (dump<=delta)  suma+=0.5*dumpx* dumpx;
		else suma += delta*dump-0.5*delta*delta;
	}

	// return (1-suma/d->nrange);

	return (AUC_PRC+0.001*(2-suma/d->nrange));
}



double HUBERL_dic(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx,  sumax=0, suma=0;
	double sumaa=0, sumab=0;
	int i;
	int wt, mut, aa;


	//int nrange;

	//nrange = d->nexp * d->fact;

	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];

		wt = map [d->aa1[i]];
		mut = map [d->aa2[i]];


		sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{

			aa = map[j];
			sumaa += x[ wt ] * x[aa] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type
			sumab += x[ mut ] * x[aa]  * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated



		}
		dumpx = sumab - sumaa; // ddG (check sign)
		d->com[i]=dumpx/100;
	}

	suma = 0.0; float dump=0; float delta=2.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		dump=fabsf(dumpx);
		if (dump<=delta)  suma+=0.5*dumpx* dumpx;
		else suma += delta*dump-0.5*delta*delta;
	}


	return (suma/d->nrange);
}




double HUBERL(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx,  sumax=0, suma=0;
	double sumaa=0, sumab=0;
	int i,wt,mut;

	//int nrange;

	//nrange = d->nexp * d->fact;

	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];
		wt =  d->aa1[i];
		mut = d->aa2[i];

		sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			sumaa += x[ wt ] * x[j]  * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type
			sumab += x[ mut ] * x[j]  * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated
		}
		dumpx = sumab - sumaa;
		d->com[i]=dumpx/100;
	}



	suma = 0.0; float dump=0; float delta=2.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		dump=fabsf(dumpx);
		if (dump<=delta)  suma+=0.5*dumpx* dumpx;
		else suma += delta*dump-0.5*delta*delta;

	}



	return (suma/d->nrange);
}




double MAE(unsigned n, const double *x, double *grad, void *data)
{

	fdata *d = (fdata *) data;

	double dumpx,  sumax=0, suma=0;
	double sumaa=0, sumab=0;
	int i;
	int wt, mut;


	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];

		wt=d->aa1[i];
		mut=d->aa2[i];


		sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{

			sumaa += x[ wt ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] )  ; // Wild-Type
			sumab += x[ mut ] * x[j]  * ( d->dcom[i][j*2] + d->dcom[i][2*j+1]) ; // Mutated

		}
		dumpx = sumab - sumaa; // ddG (check sign)
		d->com[i]=dumpx/100;
	}



	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		suma += fabsf(( d->com[i] - d->exp[i] ));
		//if (d->exp[i] > 0 )  suma += 1.0*fabsf(( d->com[i] - d->exp[i] ));
		//else 		suma += fabsf(( d->com[i] - d->exp[i] ));
	}


	return (suma/d->nrange);
}






// Mean Absolute Error --> 20 factors + 1 additive scale value
double MAE21(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;
	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0, ppv, sen;
	int		wt, mut;

	ave=0;
	for(int ii=0;ii<d->nrange;ii++) // Screen training mutation cases
	{
		int i = d->index[ii];
		wt=d->aa1[i];
		mut=d->aa2[i];

		double sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			sumaa += x[ wt ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type Sums
			sumab += x[ mut ] * x[j] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated Sums

		}
		dumpx = (sumab - sumaa)/100.0;
		d->com[i] = x[n-1] + dumpx; // Add some constant scale value
		//d->com[i] = x[20]* fabs(dumpx) + dumpx; // Add some constant scale value
		// d->com[i] = dumpx +  (x[20] + exp(-1.0*dumpx)/(1.0+pow(exp(-1.0*dumpx),2)));

	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		suma += fabsf(dumpx);
		//if (d->exp[i] > 0 )  suma += 1.0*fabsf(dumpx);
		//else 		suma += fabsf(dumpx);

	}

	// rmse = sqrt(suma/nrange);

	return (suma/d->nrange);
}

// Mean Absolute Error --> 20 factors + 1 additive scale value
double MAE22(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;
	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0, ppv, sen;
	int		wt, mut;

	ave=0;
	for(int ii=0;ii<d->nrange;ii++) // Screen training mutation cases
	{
		int i = d->index[ii];
		wt=d->aa1[i];
		mut=d->aa2[i];

		double sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			sumaa += x[ wt ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type Sums
			sumab += x[ mut ] * x[j] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated Sums



		}
		dumpx = (sumab - sumaa)/100.0; // Mutation ddG (check sign)


		d->com[i] = dumpx  + x[21];

	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		suma += fabsf(dumpx);


	}

	// rmse = sqrt(suma/nrange);

	return (suma/d->nrange);
}

// Mean Absolute Error --> 20 factors + 1 additive scale value
double MAE21_dic(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;
	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0, ppv, sen;
	int wt, i, mut, aa;
	ave=0;
	for(int ii=0;ii<d->nrange;ii++) // Screen training mutation cases
	{
		i = d->index[ii];
		wt = map [d->aa1[i]];
		mut = map [d->aa2[i]];

		double sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			aa = map[j];
			sumaa += x[ wt ] * x[aa] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type Sums
			sumab += x[ mut ] * x[aa] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated Sums
		}
		dumpx = (sumab - sumaa)/100; // Mutation ddG (check sign)


		d->com[i] = x[n-1] + dumpx; // Add some constant scale value





	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		suma += fabsf(dumpx);
		//if (d->exp[i] > 0 )  suma += 1.0*fabsf(dumpx);
		//else 		suma += fabsf(dumpx);

	}



	// rmse = sqrt(suma/nrange);

	//fprintf(stdout," suma %f %f %d\n", suma/d->nrange, x[n-1], n);
	//getchar();
	return (suma/d->nrange);
}


double MAE21rama(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;
	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0;
	int wt, mut, i;

	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];

		wt=d->aa1[i];
		mut=d->aa2[i];


		double sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			sumaa += x[ wt  ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type Sums
			sumab += x[ mut ] * x[j] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated Sums
		}

		for(int j=0; j<3; j++)
		{
			sumaa += 100*x[ 20 ] * ( d->EramaWT[i][j] ); // Wild-Type
			sumab += 100*x[ 20 ] * ( d->EramaMut[i][j] ); // Mutant
		}
		dumpx = sumab - sumaa; // ddG (check sign)
		d->com[i]=dumpx/100;
	}




	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		suma += fabsf(dumpx);
	}

	// rmse = sqrt(suma/nrange);

	return (suma/d->nrange);
}



double MAERSA(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;
	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0, ppv, sen;

	ave=0;
	for(int ii=0;ii<d->nrange;ii++) // Screen training mutation cases
	{
		int i = d->index[ii];

		double sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			sumaa += x[ d->aa1[i] ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type Sums
			sumab += x[ d->aa2[i] ] * x[j] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated Sums
		}
		dumpx = sumab - sumaa; // Mutation ddG (check sign)

		if (d->rsa[i] > 20.0)
			d->com[i] = x[20] * dumpx / 100; // Add some constant scale value
		else
			d->com[i] = x[21] * dumpx / 100; // Add some constant scale value
	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		suma += fabsf(dumpx);
	}

	// rmse = sqrt(suma/nrange);
	rmse = (suma/d->nrange);

	return (rmse);
}


double MAE40(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;
	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0, v;
	int wt, mut;

	ave=0;
	for(int ii=0;ii<d->nrange;ii++) // Screen training mutation cases
	{
		int i = d->index[ii];
		wt=d->aa1[i];
		mut=d->aa2[i];

		double sumaa=0, sumab=0;
		ave=0;
		for(int j=0;j<20;j++)
		{
			v=( d->dexp[i][j*2] + d->dexp[i][2*j+1] );
			sumaa += x[ wt ]* x[j] /(1+ exp(-1.0*(v-x[j+20])))/(1+ exp(-1.0*(v-x[wt+20])));
			//sumaa += x[ wt ]* x[j] /(1+ exp(-v*fabs(x[j+20])))/(1+ exp(-v*fabs(x[wt+20])));

			v=( d->dcom[i][j*2] + d->dcom[i][2*j+1] );
			sumab += x[ mut ]*x[j] /(1+ exp(-1.0*(v-x[j+20])))/(1+ exp(-1.0*(v-x[mut+20])));
			//sumab += x[ mut ]*x[j] /(1+ exp(-v*fabs(x[j+20]))) /(1+ exp(-v*fabs(x[mut+20])));

		}
		dumpx = sumab - sumaa; // ddG (check sign)
		//dumpx = (sumab - 100*fabs(x[ mut+ 20])-ave) - (sumaa - 100*fabs(x[ wt + 20]-ave)); // Mutation ddG (check sign)
		// dumpx = (sumab - 100*(x[ mut+ 20])) - (sumaa - 100*(x[ wt + 20])); // Mutation ddG (check sign)
		// dumpx = (sumab + 100*fabs(x[ mut+ 20])) - (sumaa + 100*fabs(x[ wt + 20])); // Mutation ddG (check sign)

		d->com[i] = dumpx / 100; // Add some constant scale value
	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		suma += fabsf(dumpx);
	}

	// rmse = sqrt(suma/nrange);
	rmse = (suma/d->nrange);

	return (rmse);
}


double MAE40_SUM(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;
	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0;
	int wt, mut;

	ave=0;
	for(int ii=0;ii<d->nrange;ii++) // Screen training mutation cases
	{
		int i = d->index[ii];
		wt=d->aa1[i];
		mut=d->aa2[i];

		double sumaa=0, sumab=0;
		ave=0;
		for(int j=0;j<20;j++)
		{
			sumaa += x[ wt ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ) ; // Wild-Type Sums
			sumab += x[ mut ] * x[j] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ) ; // Mutated Sums
			ave += (x[ j+ 20]);
		}
		ave/=20;
		dumpx = (sumab - 100*(x[ mut+ 20])-ave) - (sumaa - 100*(x[ wt + 20]-ave)); // Mutation ddG (check sign)
		//dumpx = (sumab - 100*fabs(x[ mut+ 20])-ave) - (sumaa - 100*fabs(x[ wt + 20]-ave)); // Mutation ddG (check sign)
		// dumpx = (sumab - 100*(x[ mut+ 20])) - (sumaa - 100*(x[ wt + 20])); // Mutation ddG (check sign)
		// dumpx = (sumab + 100*fabs(x[ mut+ 20])) - (sumaa + 100*fabs(x[ wt + 20])); // Mutation ddG (check sign)

		d->com[i] = dumpx / 100; // Add some constant scale value
	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		suma += fabsf(dumpx);
	}

	// rmse = sqrt(suma/nrange);
	rmse = (suma/d->nrange);

	return (rmse);
}






double MAERSA40(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;
	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0, ppv, sen;

	ave=0;
	for(int ii=0;ii<d->nrange;ii++) // Screen training mutation cases
	{
		int i = d->index[ii];

		double sumaa=0, sumab=0;
		if (d->rsa[i] > 20.0)
			for(int j=0;j<20;j++)
			{
				sumaa += x[ d->aa1[i] ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type Sums
				sumab += x[ d->aa2[i] ] * x[j] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated Sums
			}
		else
			for(int j=0;j<20;j++)
			{
				sumaa += x[ d->aa1[i] + 20] * x[j + 20] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type Sums
				sumab += x[ d->aa2[i] + 20] * x[j + 20] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated Sums
			}
		dumpx = sumab - sumaa; // Mutation ddG (check sign)

		d->com[i] = dumpx / 100; // Add some constant scale value
	}

	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		int i = d->index[ii];
		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)
		suma += fabsf(dumpx);
	}

	// rmse = sqrt(suma/nrange);
	rmse = (suma/d->nrange);

	return (rmse);
}



double MAE400(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;

	double dumpx, rmse, ave=0, dumpy, sumax=0, sumay=0, suma=0, ppv, sen;
	int nrange;

	nrange = d->nexp * d->fact; // Bootstap range

	ave=0;
	for(int ii=0;ii<nrange;ii++)
	{
		int i = d->index[ii]; // Bootstrap index from original input data

		double sumaa=0, sumab=0;
		for(int j=0;j<20;j++)
		{
			//			sumaa += x[ d->aa1[i] ] * x[j] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type
			//			sumab += x[ d->aa2[i] ] * x[j] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated
			sumaa += x[ d->aa1[i] + 20 * j ] * ( d->dexp[i][j*2] + d->dexp[i][2*j+1] ); // Wild-Type
			sumab += x[ d->aa2[i] + 20 * j ] * ( d->dcom[i][j*2] + d->dcom[i][2*j+1] ); // Mutated
		}
		dumpx = sumab - sumaa; // ddG (check sign)

		d->com[i]=dumpx/100;
	}


	suma=0;
	for(int ii=0;ii<nrange;ii++)
	{
		int i = d->index[ii];

		dumpx = ( d->com[i] - d->exp[i] ); // Computing Score (e.g. RMSE, PPV, ...)

		// suma += dumpx * dumpx;
		suma += fabsf(dumpx);

	}


	// rmse = sqrt(suma/nrange);
	rmse = (suma/nrange); // MAE



	//	return (rmse+(1-ppv));
	//	fprintf(stdout," %10.6f", rmse);
	return (rmse);

}


double MAErama(unsigned n, const double *x, double *grad, void *data)
{
	fdata *d = (fdata *) data;

	double dumpx, sumax=0, sumay=0, suma=0;
	int i;


	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];

		double sumaa = 0.0, sumab = 0.0;
		for(int j=0; j<3; j++)
		{
			sumaa += x[ 0 ] * ( d->EramaWT[i][j] ); // Wild-Type
			sumab += x[ 0 ] * ( d->EramaMut[i][j] ); // Mutant
			// sumaa += x[ d->iaaWT[i][j] ] * ( d->EramaWT[i][j] ); // Wild-Type
			// sumab += x[ d->iaaMut[i][j] ] * ( d->EramaMut[i][j] ); // Mutant
		}
		dumpx = sumab - sumaa; // ddG (check sign)

		d->com[i] = dumpx;
	}


	suma = 0.0;
	for(int ii=0;ii<d->nrange;ii++)
	{
		i = d->index[ii];
		suma += fabsf(( d->com[i] - d->exp[i] ));
	}





	return (suma/d->nrange);
	//	return (rmse+(1-ppv));

}



// Compute the average of a float array "data" of "ndata" elements
float average(float *data, int ndata)
{
	double dummy = 0.0;
	for(int i=0; i<ndata; i++)
		dummy += data[i];
	dummy /= ndata;
	return (float)dummy;
}
float average(int *data, int ndata)
{
	double dummy = 0.0;
	for(int i=0; i<ndata; i++)
		dummy += data[i];
	dummy /= ndata;
	return (float)dummy;
}

// Compute the sigma of a float array "data" of "ndata" elements. Optionally, the average can be provided if already known.
float sigma(float *data, int ndata, float avg=0.0)
{
	double dummy = 0.0;
	if(avg == 0.0)
		avg = average(data,ndata);
	for(int i=0; i<ndata; i++)
		dummy += pow(data[i]-avg,2);
	dummy /= ndata;
	return (float)sqrt(dummy);
}
float sigma(int *data, int ndata, float avg=0.0)
{
	double dummy = 0.0;
	if(avg == 0.0)
		avg = average(data,ndata);
	for(int i=0; i<ndata; i++)
		dummy += pow(data[i]-avg,2);
	dummy /= ndata;
	return (float)sqrt(dummy);
}

// Quick-sort recursive routine to sort everything in between `low' <-> `high'
//  arr   --> Input array of values (float) [v0,v1,v2,...vN-1]
//  index --> Input array of indices (int) [0,1,2,...,N-1]
//  low   --> Low pivot point (for recursion)
//  high  --> High pivot point (for recursion)
void quicksort(float *arr, int *index, int low, int high)
{
	int y; // dummy variable for swap ints
	int i = low; // left
	int j = high; // right
	float z = arr[ index[(low + high) / 2] ]; // pivot

	do {
		while(arr[index[i]] < z) // find member above
			i++;

		while(arr[index[j]] > z) // // find member below
			j--;

		if(i <= j) // Swap two elements
		{
			swap(index[i],index[j],y);
			i++;
			j--;
		}
	} while(i <= j);

	// recurse
	if(low < j)
		quicksort(arr, index, low, j);

	if(i < high)
		quicksort(arr, index, i, high);
}

// Compute the median from "data" array of "ndata" elements
float median(float *data, int ndata)
{
	// Sort array of RMSDs
	int *sorted = (int *) malloc(sizeof(int) * ndata); // Allocate memory for the SORTed indices array for RMSDs (sortr)
	for(int u = 0; u < ndata; u++) // Initialize sorted array of RMSDs
		sorted[u] = u; // Load with the sequential indices
	quicksort(data, sorted, 0, ndata-1); // Quick sort algorithm

	if(ndata == 1)
		return data[0]; // The unique value
	else
		if(ndata % 2 == 0) // even
			return (data[ sorted[ndata/2] ] + data[ sorted[(ndata/2)-1] ]) / 2; // Average of the two median values
		else // odd
			return data[ sorted[ndata/2] ]; // The median value
}

// Compute the Pearson's correlation coefficient applied to a sample ("r" factor)
float corrPearson(float *x, float *y, int n)
{
	float ax = average(x,n);
	float ay = average(y,n);
	float sx = sigma(x,n,ax);
	float sy = sigma(y,n,ay);

	double r = 0.0;
	for(int i=0; i<n; i++)
		r += (x[i]-ax)*(y[i]-ay);
	r /= n;

	return (float) (r/(sx*sy));
}

// Compute the z-score of "n" data elements ("data") wrt some reference value ("ref")
float zscore(float *data, float ref, int n)
{
	float avg = average(data, n);
	float sig = sigma(data, n, avg);
	return (avg - ref) / sig;
}

// Get the maximum value
float getmax(float *a, int n)
{
	float max = a[0];
	for(int i=1; i<n; i++)
		if(a[i] > max)
			max = a[i];
	return max;
}

// Get the minimum value
float getmin(float *a, int n)
{
	float min = a[0];
	for(int i=1; i<n; i++)
		if(a[i] < min)
			min = a[i];
	return min;
}


/**************************************************/

int main( int argc, char * argv[] )
{
	bool debug = false;
	bool sym_switch = false; // Enable symmetric mutations matrix
	bool pclase = false;
	bool pkrand = false;
	bool rsa_switch = false; //Read Relative Surface Area (RSA) from input file. RSA value is expected at 4th or 5th (if --class) columns
	bool avg_switch = true;  // average switch if true average KD
	bool save_switch = true;  //  save results
	int verbose=1;

	std::string temp;
	char input[256]; // Input text file
	int run_mode = -1; // run mode
	float wrange = 40.0; // Weights range for optimization

	float classratio = 0.2; // Class overflow factor
	float minfactor = 0.8; // Factor to determine the minimum number of cases per k-fold bin
	int minmut = 1; // minum numeber of aa cases per k-fold
	int KN = 10;   // K in K-fold test
	int NRUNS = 10; // Number of runs
	int KNNRUNS = KN * NRUNS; // Number of items for statistics
	int ENRUNS=NRUNS; // effective average stats

	FILE *fout,*fout2 ;

	if ( (fout=fopen("kd_stats.txt", "w"))==NULL)
	{
		fprintf(stderr, "\n  Error->Cannot open file\n");
	}

	if (save_switch)
		if ( (fout2=fopen("out.txt", "w"))==NULL)
		{
			fprintf(stderr, "\n  Error->Cannot out.txt open file\n");
		}



	// COMMAND-LINE PARSER:
	using namespace TCLAP;
	CmdLine cmd(argv[0],"   ", VERSION );
	try {

		UnlabeledValueArg<std::string> Input("list","List of PDB/MUT\n format \"1BNI IA88V -0.931\" ","default","mutations");
		cmd.add( Input );

		ValueArg<int> Runmode("r","runmode","Energy/Run mode: 5= KORPweight, 6= RAMAweight, 56= KORPweight+RAMAweight, 21= 20weights+1(additive factor)", false, 5,"int");
		cmd.add( Runmode );

		ValueArg<float> Wrange("","wrange","Upper/lower optimization range for weights", false, 5.0,"float");
		cmd.add( Wrange );

		SwitchArg Sym("","sym", "Enable symmetric mutations matrix", false);
		cmd.add( Sym );

		SwitchArg Pclase("","class", "Classified mutations input file. Each line from input mutations list should look like this: \"1BNI IA88V -0.931 27\", where \"27\" (the 4th column) is the class index.", false);
		cmd.add( Pclase );

		ValueArg<int> Kfold("k","kfold","K in a K-fold test (default=5)", false, 5,"int");
		cmd.add( Kfold );

		ValueArg<int> Nruns("n","nruns","Number of runs (default=10)", false, 10,"int");
		cmd.add( Nruns );

		ValueArg<int> Ver("","verbose","verbose mode 0/1/2", false, 1,"int");
		cmd.add( Ver );

		ValueArg<float> ClassRatio("","classratio","Class overflow ratio (classratio), i.e. aproxkbins = classratio * N/k. Values around 0.2 should be Ok. (default=0.2)", false, 0.2,"float");
		cmd.add( ClassRatio );

		ValueArg<float> MinFactor("","minfactor","Factor to determine the minimum number of cases per k-fold bin (mincases), i.e. mincases = minfactor * aproxkbins. (default=0.9)", false, 0.9,"float");
		cmd.add( MinFactor );

		ValueArg<int> MinMut("","minmut","Minimum number of cases of mutation per aa type in the k-fold , e.g 1 at leat one mut per aa (default=1)", false, 0.9,"float");
		cmd.add( MinMut );

		SwitchArg RSA("","rsa", "Read Relative Surface Area (RSA) from input file. RSA value is expected at 4th or 5th (if --class) columns.", false);
		cmd.add( RSA);

		SwitchArg Krand("","krand", "force random kfold", false);
		cmd.add( Krand);

		SwitchArg Avg("","kdavg", "average inside KD", false);
		cmd.add( Avg );


		// Parse the command line.
		cmd.parse(argc,argv);

		strcpy(input,((temp=Input.getValue()).c_str()));

		run_mode = Runmode.getValue();
		wrange = Wrange.getValue();
		verbose = Ver.getValue();;

		sym_switch = Sym.isSet(); // Enable symmetric mutations matrix
		avg_switch = Avg.isSet(); // unabable average kd

		if (avg_switch)
			fprintf(stderr, " Average all K-folds\n");

		if(Pclase.isSet())
			pclase = true;

		if(Krand.isSet())
			pkrand = true;

		if(RSA.isSet())
			rsa_switch = true;

		KN = Kfold.getValue();   // K in K-fold test

		NRUNS = Nruns.getValue(); // Number of runs

		KNNRUNS = KN*NRUNS; // Number of items for statistics

		classratio = ClassRatio.getValue(); // Class overflow ratio
		minmut = MinMut.getValue(); // Factor to determine the minimum number of cases per k-fold bin
		minfactor = MinFactor.getValue(); // Factor to determine the minimum number of cases per k-fold bin
	}
	catch ( ArgException& e )
	{
		std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl;
	}
	using namespace std;

	// INPUT DATA
	fdata d; // Main data structure with all stuff
	time_t t;
	srand((unsigned) time(&t));
	char dumpc;
	char dumps[20];
	int stats[20][20];

	for(int i=0;i<20;i++)
		for(int j=0;j<20;j++)
			stats[j][i]=0;

	// Classes statistics
	int statsH[500][1000]; // type / list of clasess in [0] store the #
	for(int j=0;j<500;j++)
		statsH[j][0]=0;

	int aa1, aa2;
	float dumpF;
	float exp1, com, exp2, com2;
	FILE *f;

	// Store stuff to compute general statistics
	float *Rmse; // Root Mean Squared Error array
	float *Mae; // Mean Absolute Error array
	float *Ppv; // Positive Predictive Power array
	float *Sen; // Sensitivity array
	float *Pcor; // Pearson's cross-correlation
	float *Acc; // ACCuracy
	float *Err; // ERror Rate
	float *Mcc; // Matthews correlation coefficient
	float *Fpr; // False Positive Rate (FPR)
	float *auc_roc; // auc roc
	float *auc_prc; // auc prc


	arec= (float *) malloc(sizeof(float) * BINSAUC); if (arec==NULL) exit (1);
	appv= (float *) malloc(sizeof(float) * BINSAUC); if (appv==NULL) exit (1);
	afpr= (float *) malloc(sizeof(float) * BINSAUC); if (afpr==NULL) exit (1);


	if (avg_switch)
		ENRUNS=NRUNS;
	else
		ENRUNS=KNNRUNS;


	float **Weights = (float**)malloc(NVAR * sizeof(float*));
	for (int index=0;index<NVAR;++index)
	{
		Weights[index] = (float*)malloc(KNNRUNS * sizeof(float));
	}

	Arec= (float**)malloc(sizeof(float*) * BINSAUC);
	Appv= (float**)malloc(sizeof(float*) * BINSAUC);
	Afpr= (float**)malloc(sizeof(float*) * BINSAUC);

	for (int index=0;index<BINSAUC;++index)
	{
		Arec[index] = (float*)malloc(KNNRUNS * sizeof(float));
		Appv[index] = (float*)malloc(KNNRUNS * sizeof(float));
		Afpr[index] = (float*)malloc(KNNRUNS * sizeof(float));
	}

	for (int i=0;i<BINSAUC;i++)
		for (int j=0;j<KNNRUNS;j++){
			Arec[i][j]=0;
			Appv[i][j]=0;
			Afpr[i][j]=0;
		}




	Rmse= (float *) malloc(sizeof(float) * ENRUNS); if (Rmse==NULL) exit (1);
	Mae= (float *) malloc(sizeof(float) * ENRUNS); if (Mae==NULL) exit (1);
	Ppv= (float *) malloc(sizeof(float) * ENRUNS); if (Ppv==NULL) exit (1);
	Sen= (float *) malloc(sizeof(float) * ENRUNS); if (Sen==NULL) exit (1);
	Pcor= (float *) malloc(sizeof(float) * ENRUNS); if (Pcor==NULL) exit (1);
	Acc= (float *) malloc(sizeof(float) * ENRUNS); if (Acc==NULL) exit (1);
	Err= (float *) malloc(sizeof(float) * ENRUNS); if (Err==NULL) exit (1);
	Mcc= (float *) malloc(sizeof(float) * ENRUNS); if (Mcc==NULL) exit (1);
	Fpr= (float *) malloc(sizeof(float) * ENRUNS); if (Fpr==NULL) exit (1);
	auc_prc= (float *) malloc(sizeof(float) * ENRUNS); if (auc_prc==NULL) exit (1);
	auc_roc= (float *) malloc(sizeof(float) * ENRUNS); if (auc_roc==NULL) exit (1);





	// Open input file
	strcpy(input,argv[1]);
	if ( (f=fopen(input, "r"))==NULL)
	{
		fprintf(stderr, "\n  Error->Cannot open file\n");
	}

	// Main program loop
	int i = 0; // case index
	while( !feof(f) )
	{
		int check;
		if (pclase)
			if(rsa_switch)
				check  = fscanf(f, "%*s %s %*c %*c %d %d %d %d %d %f %f %*f %*f %*f %*f %*d %d %f", &dumps, &aa1, &aa2, &d.iaaWT[i][0], &d.iaaWT[i][1], &d.iaaWT[i][2], &exp1, &com, &d.Hclass[i], &d.rsa[i]);
			else
				check  = fscanf(f, "%*s %s %*c %*c %d %d %d %d %d %f %f %*f %*f %*f %*f %*d %d", &dumps, &aa1, &aa2, &d.iaaWT[i][0], &d.iaaWT[i][1], &d.iaaWT[i][2], &exp1, &com, &d.Hclass[i]);
		else
			if(rsa_switch)
				check  = fscanf(f, "%*s %s %*c %*c %d %d %d %d %d %f %f %*f %*f %*f %*f %*d %f", &dumps, &aa1, &aa2, &d.iaaWT[i][0], &d.iaaWT[i][1], &d.iaaWT[i][2], &exp1, &com, &d.rsa[i]);
			else
				check  = fscanf(f, "%*s %s %*c %*c %d %d %d %d %d %f %f %*f %*f %*f %*f %*d", &dumps, &aa1, &aa2, &d.iaaWT[i][0], &d.iaaWT[i][1], &d.iaaWT[i][2], &exp1, &com);

		if( check < 0 )
			break;


		// if (exp1<0.0) continue;

		//		fprintf(stderr,"i= %4d  rsa= %7.2f\n",i,d.rsa[i]);

		// Store mutant type indices for Rama
		d.iaaMut[i][0] = d.iaaWT[i][0];
		d.iaaMut[i][1] = aa2;
		d.iaaMut[i][2] = d.iaaWT[i][2];
		d.aa1[i]=aa1;
		d.aa2[i]=aa2;
		d.exp[i]=exp1;
		d.com[i]=com;
		d.aveE += exp1;
		d.aveC += com;
		stats[aa1][aa2]++;

		if( pclase ) // Classes stuff
		{
			statsH[ d.Hclass[i] ][ 0 ]++;
			statsH[ d.Hclass[i] ][ statsH[d.Hclass[i]][0] ] = i; // storing the position of the class
		}

		// Inverse stuff
		if(sym_switch) // Enable symmetric mutations matrix
		{
			d.iaaMut[i+1][0] = d.iaaWT[i][0];
			d.iaaMut[i+1][1] = aa2;
			d.iaaMut[i+1][2] = d.iaaWT[i][2];
			d.aa1[i+1]=aa2;
			d.aa2[i+1]=aa1;
			d.exp[i+1] = -exp1;
			d.com[i+1] = -com;
			d.aveE += -exp1;
			d.aveC += -com;
			stats[aa2][aa1]++;
		}

		// Get wild-type energies sums
		for(int j=0;j<20;j++)
		{
			fscanf(f," %f %f", &exp1,&com); // NOW exp and com are dummies
			d.dexp[i][j*2] = exp1; // dummy
			d.dexp[i][j*2+1] = com; // dummy

			// Inverse stuff
			if(sym_switch) // Enable symmetric mutations matrix
			{
				d.dcom[i+1][j*2] = exp1;
				d.dcom[i+1][j*2+1] = com;
			}
		}

		// Get mutant energies sums
		for(int j=0;j<20;j++)
		{
			fscanf(f," %f %f", &exp1,&com); // NOW exp and com are dummies
			d.dcom[i][j*2] = exp1;
			d.dcom[i][j*2+1] = com;

			// Inverse stuff
			if(sym_switch) // Enable symmetric mutations matrix
			{
				d.dexp[i+1][j*2] = exp1; // dummy
				d.dexp[i+1][j*2+1] = com; // dummy
			}
		}

		// Get Rama sums
		for(int j = 0; j < 3; j++)
		{
			fscanf(f," %f %f", &d.EramaWT[i][j], &d.EramaMut[i][j] );
			// fprintf(stderr," %f %f", d.EramaWT[i][j], d.EramaMut[i][j] );

			if(sym_switch) // Enable symmetric mutations matrix
			{
				d.EramaWT[i+1][j] = d.EramaMut[i][j];
				d.EramaMut[i+1][j] = d.EramaWT[i][j];
			}
		}

		// if ((d.exp[i]<-0.1)&&((aa1==2)||(aa2==2)||(aa1==3)||(aa2==3))) {
		//    if ((d.exp[i]>-0.1)&&((aa1==2)||(aa2==2)||(aa1==3)||(aa2==3)||(aa1==14)||(aa2==14)||(aa1==8)||(aa2==8)||(aa1==6)||(aa2==6))) {

		// if (d.exp[i]> -1000) {
		//	 if (	 d.rsa[i] > 20 && (d.exp[i]>0.0001)) {
		// pending add input flag
		if (	 d.rsa[i] < 9990.7) {

			//		if ((d.exp[i]>0.0001)&&(aa1!=14)&&(aa2!=14)&&(aa1!=8)&&(aa2!=8)&&(aa1!=16)&&(aa2!=16)&&(aa1!=15)&&(aa2!=15)&&(aa1!=1)&&(aa2!=1)&&(aa1!=2)&&(aa2!=2)&&(aa1!=3)&&(aa2!=3)&&(aa1!=5)&&(aa2!=5))   {
			//	fprintf(stderr,"%d %d %d %f %f %d\n",i,aa1,aa2,d.exp[i], d.rsa[i], d.Hclass[i]);

			if(sym_switch) // Enable symmetric mutations matrix
				i += 2; // Update main index
			else
				i++; // Update main index
		}
	}
	fclose(f);

	d.nexp = i;
	d.aveE /= (float) i;
	d.aveC /= (float) i;

	if ( (f=fopen("table1.csv", "w"))==NULL)
		fprintf(stderr, "\n  Error->Cannot table file\n");


	// Dump aminoacids count statistics
	int c1[20],c2[20];
	if (verbose>0) fprintf(stdout,"\n    ");
	for(int j=0;j<20;j++)
	{
		if (verbose>0) fprintf(stdout,"  %c  ", aa[j]);
		fprintf(f,",%c", aa[j]);

	}
	if (verbose>0) fprintf(stdout,",X\n");
	fprintf(f,"\n");


	for(int i=0;i<20;i++)
	{
		c2[i]=0; c2[i]=0;

	}
	for(int i=0;i<20;i++)
	{
		if (verbose>0) fprintf(stdout,"%c ", aa[i]);
		fprintf(f,"%c,", aa[i]);

		c1[i]=0;
		for(int j=0;j<20;j++)
		{
			if (verbose>0) fprintf(stdout,"%5d", stats[i][j]);
			fprintf(f,"%d,", stats[i][j]);

			c1[i]+=stats[i][j];
			c2[i]+=stats[j][i];
		}
		if (verbose>0) fprintf(stdout,"%5d\n", c1[i]);
		fprintf(f,"%d\n", c1[i]);

	}
	if (verbose>0) fprintf(stdout,"");
	fprintf(f,"X");

	for(int j=0;j<20;j++) {
		if (verbose>0) fprintf(stdout,"%5d", c2[j]);
		if (i==19)
			fprintf(f,"%d", c2[j]);
		else
			fprintf(f,",%d", c2[j]);

	}
	if (verbose>0) fprintf(stdout,"\n");
	fprintf(f,",\n");
	fclose(f);

	if (verbose>0) {
		fprintf(stdout,"\n    ");
		for(int j=0;j<20;j++)
			fprintf(stdout,"  %c  ", aa[j]);
		fprintf(stdout,"\n");
	}

	for(int i=0;i<20;i++)
	{
		if (verbose>0)  fprintf(stdout,"%c ", aa[i]);
		for(int j=0;j<20;j++)
		{
			int dump=0;
			dump = stats[i][j] + stats[j][i];
			if(i >= j)
				if (verbose>0) fprintf(stdout,"%5d", dump);
		}
		if (verbose>0) fprintf(stdout,"\n");
	}

	if (verbose>0) fprintf(stdout,"\n    ");

	if ( (f=fopen("table2.csv", "w"))==NULL)
		fprintf(stderr, "\n  Error->Cannot table file\n");

	for(int i=0;i<20;i++)
	{
		if (verbose>0) fprintf(stdout,"%c    ", aa[i]);
		if (i==0) 		fprintf(f,"%c", aa[i]);
		else fprintf(f,",%c", aa[i]);
	}
	if (verbose>0) fprintf(stdout,"\n");
	fprintf(f,"\n");


	for(int i=0;i<20;i++)
	{
		if (verbose>0) fprintf(stdout,"%5d", (c1[i]+c2[i]));
		if (i==0) fprintf(f,"%d", (c1[i]+c2[i]));
		fprintf(f,",%d", (c1[i]+c2[i]));
	}
	if (verbose>0) fprintf(stdout," \t aa <-> aa\n");
	fprintf(f,"\n");


	for(int i=0;i<20;i++)
	{
		if (verbose>0) fprintf(stdout,"%5d", c1[i]);
		if (i==0) fprintf(f,"%d", c1[i]);
		fprintf(f,",%d", c1[i]);
	}
	if (verbose>0) fprintf(stdout," \t aa <-> xx\n");
	fprintf(f,"\n");

	for(int i=0;i<20;i++)
	{
		if (verbose>0) fprintf(stdout,"%5d", c2[i]);
		if (i==0) fprintf(f,"%d", c2[i]);
		fprintf(f,",%d", c2[i]);

	}
	if (verbose>0) fprintf(stdout," \t xx <-> aa\n");
	fprintf(f,"\n");
	fclose(f);


	int nclasses = 0;
	int imaxclass = 0;
	int nmutstat = 0; // number of mutations for cross-checking

	if( pclase ) // Classes stuff
	{
		if (verbose>0) fprintf(stdout,"\n");
		if ( (f=fopen("table3.csv", "w"))==NULL)
			fprintf(stderr, "\n  Error->Cannot table file\n");
		fprintf(f,"index, class, N\n");

		for(int j=0;j<500;j++)
		{
			if( statsH[j][0] != 0 )
			{
				if (verbose>0) fprintf(stdout,"%3d(%3d)  ", j, statsH[j][0]);
				fprintf(f,"%3d,%3d,%3d\n", nclasses+1, j, statsH[j][0]);

				nclasses++;
				if( nclasses%15==0 && nclasses!=0 )
					if (verbose>0) fprintf(stdout,"\n");

				if(imaxclass < j)
					imaxclass = j;

				nmutstat += statsH[j][0]; // count number of mutations for cross-checking
			}
		}
		if (verbose>0) fprintf(stdout,"\n\n#classes %3d  imaxclass %3d  nmutstat %4d %4d\n", nclasses, imaxclass, nmutstat,d.nexp);
		fclose(f);

	}


	int posi=0;
	int nega = 0;
	int nzero =0;

	for(int ii=0;ii<d.nexp;ii++)
	{
		if (d.exp[ii]>0) posi++;
		else if (d.exp[ii]<0) nega++;
		else nzero++;
		if(d.exp[ii] < minP) minP=d.exp[ii];
		if(d.exp[ii] > maxP) maxP=d.exp[ii];
	}
	if (verbose>0) fprintf(stdout,"#S %5d(%3.1f)  D %5d(%3.1f)  ZERO %5d TOTAL %5d  DDexp MAX %.2f MIN %.2f \n", posi, (100.0*posi/d.nexp), nega, (100.0*nega/d.nexp), nzero, d.nexp, maxP,minP );

	// getchar();

	// MON: WATCH OUT!!!
	// grep " T T" out.txt
	// pdbs/1PGA.pdb    TA53T T T 16 16  4 16 17

	//
	//  OPTIMIZATION
	//
	double f0;  // fitness function
	double f0min;
	double f0avg;
	double f0corre;
	double f0correavg;
	double sen; // Sensitivity (SN) = RECall = True Positive Rate (TPR)
	double acc; // ACCuracy
	double err; // ERror Rate
	double mcc; // Matthews correlation coefficient
	double fpr; // False Positive Rate (FPR)
	double ppv; // Positive Predictive Value = PRECision
	double ppvavg;
	double ppvmax;
	double x[NVAR]; // opt vars
	double lb[NVAR]; // lower bounds
	double ub[NVAR];   // upper bounds
	double xmax[NVAR];
	double xavg[NVAR];
	int NEFF=0;

	// NEW ALL DATA
	// double xi[NVAR]={ 1.62906,1.25222,1.16860,1.38073,2.36187,1.41764,1.08243,2.26649,1.41488,2.85751,1.49585,1.61682,1.43767,0.99600,1.93309,1.22990,1.59651,2.31484,1.44050,2.58258};
	double xi[NVAR];

	// INITIALIZATION
	int index=0;
	int id=0;
	d.fact=NFACTOR;

	int nrange = d.nexp*d.fact;


	switch(run_mode)
	{
	case 3:                // map korp
	case 4:               // map single weight

		nvar = -1;
		for(int i=0; i<21; i++)
		{
			if (map[i]>nvar) nvar=map[i];
		}
		break;


	case 5:
		nvar = 20;
		break;

	case 6:
		nvar = 1;
		break;

	case 56:
		nvar = 21; // MAE
		break;

	case 400:
		nvar = 400;
		break;

	case 21:
		nvar = -1;
		for(int i=0; i<21; i++)
		{
			if (map[i]>nvar) nvar=map[i];
		}
		nvar++;
		//		nvar = 21;


		break;


	case 22:
		nvar = 22;
		break;

	case 40:
		nvar = 40;
		break;

	default:
		fprintf(stderr,"Error, unknown run_mode %d, forcing exit!\n", run_mode);
		exit(1);
		break;
	}
	if (verbose>0) fprintf(stdout, "-> nvar %d\n\n", nvar);


	char fmt[10];
	if(run_mode == 400)
		strcpy(fmt, " %5.1f");
	else
		strcpy(fmt, " %7.3f");

	if (verbose>0) fprintf(stdout,"init         ");
	for(int j=0;j<nvar;j++)
	{
		if (verbose>0) fprintf(stdout, fmt, xi[j]);
		xavg[j] = 0.0;
	}
	if (verbose>0) fprintf(stdout,"\n");


	f0min = 10000;
	f0avg = 0.0;
	ppvmax = -1;
	ppvavg = 0;
	f0correavg = 0;


	// Random number
	std::random_device rnd;
	std::mt19937 gen( rnd() );
	// std::mt19937 gen(0); // static seed
	// std::uniform_real_distribution<> dis(-1.,1.);

	// Initializing of uniform_int_distribution class
	uniform_int_distribution<int> distribution(0, d.nexp - 1);
	uniform_real_distribution<double> distrFloat(0, 1);


	// Weights initialization
	for(int i=0;i<nvar;i++)
		xi[i] = 1.0;


	int kn_start[KN];
	int knbins[KN];

	// Sort randomly by class
	int *classindex; // Aux list the class
	int *classcountN;
	int *classcountT;
	float *classcountF; // Aux list the class

	classcountF= (float *) malloc(sizeof(float) * (imaxclass+1)); if (classcountF==NULL) exit (1);
	classcountT= (int *) malloc(sizeof(int) * (imaxclass+1)); if (classcountT==NULL) exit (1);
	classcountN= (int *) malloc(sizeof(int) * (imaxclass+1)); if (classcountN==NULL) exit (1);
	classindex= (int *) malloc(sizeof(int) * (imaxclass+1)); if (classindex==NULL) exit (1);



	//	// Statistics per class
	//	float statclass[1000]; // Statistics per class
	//	for (int i = 0; i < 1000; i++)
	//		statclass[i] = 0.0;

	int reset_class=0;
	int reset_mut=0;

	if (pkrand) pclase=false;

	// Main runs loop
	for (int r = 0; r < NRUNS; r++)
	{
		if (avg_switch)
		{
			Rmse[r] = 0;
			Mae[r]  = 0; // Mean Absolute Error array
			Ppv[r]  = 0; // Positive Predictive Power array
			Sen[r]  = 0; // Sensitivity array
			Pcor[r] = 0; // Pearson's cross-correlation
			Acc[r]  = 0; // ACCuracy
			Err[r]  = 0; // ERror Rate
			Mcc[r]  = 0; // Matthews correlation coefficient
			Fpr[r]  = 0; // False Positive Rate (FPR)
			auc_prc[r] = 0;
			auc_roc[r] = 0;
		}

		int cont_mut[20][KN];
		bool checkZero=false;



		if ((!pclase))  // Standard k-fold (no classes)
			//		if(true)
		{
			int delta_kn = d.nexp / KN;

			if ((r==0) &&(verbose>0))
			{
				fprintf(stdout,"\n");
				fprintf(stdout,"Number of cases %d K-%-d fold training subset %d\n",d.nexp, KN, delta_kn );
			}
			// Reset auxiliary indices array
			for (int i = 0; i < d.nexp; i++)
				d.index2[i] = i;

			// Split Kn training and test sets
			for (int kn = 0; kn < KN; kn++)
			{
				kn_start[kn] = 0;
				knbins[kn] = delta_kn;
			}
			knbins[KN-1] = d.nexp - (KN-1) * delta_kn;

			for (int kn = 1; kn < KN; kn++)
				kn_start[kn] = kn_start[kn-1] + knbins[kn];

			// Shuffle (randomize indices)
			for(int k = 0; k < 5; k++) // repeat 5 times for "good" randomness
			{
				for (int i = 0; i < d.nexp; i++)
				{
					// int myrand = rand() / (RAND_MAX / (d.nexp - i) + 1);
					int j = distribution(gen); // get a random number between 0 and d.nexp-1

					// Swapping a pair of values
					int id =  d.index2[j];
					d.index2[j] =  d.index2[i];
					d.index2[i] = id;

					// fprintf(stderr,"myrand= %6d  i= %4d  j= %4d\n",myrand,i,j);
					// fprintf(stderr,"i= %4d  j= %4d\n",i,j);
				}


				//fprintf(stdout,"\n");
				// CHECK FOR MINUM AA  counts....
				if (k==4) {
					checkZero=false;

					for (int kn = 0; kn < KN; kn++)
					{
						for (int i = 0; i < 20; i++)
							cont_mut[i][kn]=0;
						// fprintf(stdout,"\n Kd %4d ", kn);

						for (int j = kn_start[kn]; j < kn_start[kn] + knbins[kn]; j++)
						{
							cont_mut[d.aa1[d.index2[j]]][kn]++;
							cont_mut[d.aa2[d.index2[j]]][kn]++;
						}
						for (int i = 0; i < 20; i++)
						{

							// fprintf(stdout,"%4d ", cont_mut[i][kn]);
							if (cont_mut[i][kn] < minmut)
							{
								checkZero=true;
								break;
							}
						}
						// fprintf(stdout,"\n");
					}


					if (checkZero)
					{
						k=0;
						// fprintf(stdout," Zero mutation in classes\n");
						reset_mut++;

						if (reset_mut>5000) {
							fprintf(stderr,"\n Unbalanced data set impossible to generate K-fold, check --minmun\n");
							exit(1);
						}

						// getchar();

					}
				}

			}



		}
		else // Select by class
		{



			for (int i = 0; i <= imaxclass; i++)
			{
				classindex[i] = 0;
				classcountN[i] = 0;
				classcountT[i] = 0;
				classcountF[i] = 0.0;

			}

			// reset clases....maybe out of the loop
			int realnclass = 0;
			int nmutstat2 = 0;

			for (int i = 0; i <= imaxclass; i++)
			{
				if (statsH[i][0] != 0)
				{
					classindex[realnclass] = i;
					realnclass++;

					nmutstat2 += statsH[i][0];
				}
			}
			//fprintf(stderr, "imaxclass %3d  nmutstat2 %4d  realnclass %3d\n", imaxclass, nmutstat2, realnclass);


			// Shuffle (randomize indices of classes)
			uniform_int_distribution<int> distribution2(0, realnclass-1);


			for(int k = 0; k < 5; k++) // repeat 5 times for "good" randomness
			{
				for (int i = 0; i < realnclass; i++)
				{
					// int myrand = rand() / (RAND_MAX / (d.nexp - i) + 1);
					int j = distribution2(gen); // get a random number between 0 and realnclass

					// Swapping a pair of values
					int id =  classindex[j]; // backup old value
					classindex[j] =  classindex[i]; // set new value
					classindex[i] = id; // set old value swapped
				}
			}


			int aproxkbins = (1-classratio) * d.nexp / KN;
			int aproxkbinsAvg =aproxkbins;

			int countkn = 0;
			int kn = 0;
			int countknT = 0;
			int classcount = 0;
			int countkn_old = 0;
			for (int kn = 0; kn < KN; kn++)
			{
				knbins[kn]=0;
				kn_start[kn]=0;

			}

			for (int i = 0; i < realnclass; i++) // Screen "real" classes
			{
				classcount++;
				countkn_old = countkn;
				countkn += statsH[classindex[i]][0]; // Accumulating number of mutations for each class

				if (countkn > aproxkbins ) // if its over "aproxkbins"...
				{
					if( countkn > aproxkbins * (1+classratio) and classcount != 1 and countkn_old > aproxkbins*(1-classratio) )
					{   // exceeds too much go one back
						countkn = countkn_old;
						classcount--;
						i--;
					}

					knbins[kn] = countkn;

					if(kn == 0)
						kn_start[kn] = 0;
					else
						kn_start[kn] = countknT;


					countknT += countkn;
					countkn = 0; // reset
					classcountN[kn] = classcount;
					classcount = 0; // reset
					aproxkbins= (1-classratio)*(d.nexp-kn_start[kn] - knbins[kn])/(KN-kn-1);
					aproxkbinsAvg += aproxkbins;

					//fprintf(stderr,"%d %4d %4d %4d %f\n", KN-kn, kn_start[kn] + knbins[kn], aproxkbinsAvg,aproxkbins,classratio-(kn)*classratio/(KN*1.0));
					kn++;

				}



				if (kn == KN-1)
					break;
			}

			classcount = 0;
			for (int kn = 0; kn < KN-1; kn++)
				classcount += classcountN[kn];
			classcountN[KN-1] = realnclass - classcount;

			kn_start[KN-1] = kn_start[KN-2] + knbins[KN-2];
			knbins[KN-1] = d.nexp - kn_start[KN-1];
			checkZero=false;
			if (knbins[KN-1] == 0) 	{
				//fprintf(stdout,"Kbins zero %4d ", knbins[KN-1]);
				checkZero=true;
			}

			// Minimum number of cases per bin test
			int mincases = minfactor * aproxkbinsAvg/KN;
			if (mincases<d.nexp /KN*0.5)
				checkZero=true;


			for(int kn = 0; kn < KN; kn++) {
				if( knbins[kn] < mincases) {
					//fprintf(stdout,"Knbins kd %4d->%d %d %d\n", kn, kn_start[kn] + knbins[kn], knbins[kn], mincases);
					checkZero=true;

				}
				//getchar();
			}


			if (checkZero)
			{
				r--;
				//				fprintf(stdout," Zero classes restart\n");
				reset_class++;
				if (reset_class>5000) {
					fprintf(stderr,"\n Unbalanced data set impossible to generate K-fold, check minfactor\n");
					exit(1);
				}

				// getchar();
				continue;
			}



			// sort cases by randomized cluster list
			int indexP = 0;
			int countclass; // number of mutations per class


			for (int i = 0; i < realnclass; i++)
			{
				countclass = 0; // initialize the number of mutations per class

				for (int j = 0; j < d.nexp; j++)
				{
					if( d.Hclass[j] == classindex[i] )
					{
						d.index2[indexP] = j;
						indexP++;
						countclass++; // count the number of mutations per class
					}

					if( countclass == statsH[ classindex[i] ][0] )
						break;
				}
			}

			int cont_mut[20][KN];

			//fprintf(stdout,"\n");
			for (int kn = 0; kn < KN; kn++)
			{
				for (int i = 0; i < 20; i++)
					cont_mut[i][kn]=0;
				// fprintf(stdout,"Kd %4d ", kn);

				for (int j = kn_start[kn]; j < kn_start[kn] + knbins[kn]; j++)
				{
					cont_mut[d.aa1[j]][kn]++;
					cont_mut[d.aa2[j]][kn]++;
				}
				for (int i = 0; i < 20; i++)
				{

					if (cont_mut[i][kn] < minmut)
					{
						// fprintf(stdout,"min %4d\n", cont_mut[i][kn]);
						checkZero=true;
						break;
					}
				}
				// fprintf(stdout,"\n");
			}


			if (checkZero)
			{
				r--;
				// fprintf(stdout," Zero mutation in classes\n");
				reset_mut++;

				if (reset_mut>5000) {
					fprintf(stderr,"\n Unbalanced datase set imposible to generarate K-fold, check --minmun\n");
					exit(1);
				}

				// getchar();
				continue;
			}
			if (verbose>0)
				fprintf(stdout,"\nTraining subset from %d\n",d.nexp );


			int nclasscheck = 0;
			if (verbose>0) fprintf(stdout,"\n\nK-fold exp %d/%d # class limit  %4d resets %d %d\n", r, NRUNS, aproxkbinsAvg/(KN-1), reset_class, reset_mut);
			for (int kn = 0; kn < KN; kn++)
			{
				if (verbose>0) fprintf(stdout,"kd %4d %4d %4d %4d %4d->", kn, kn_start[kn], knbins[kn], kn_start[kn] + knbins[kn], classcountN[kn]);
				for (int i = nclasscheck; i < nclasscheck+classcountN[kn]; i++)
					if (verbose>0) fprintf(stdout,"%4d ", classindex[i]);
				if (verbose>0) fprintf(stdout,"\n");

				nclasscheck += classcountN[kn];
			}

			if (verbose>0)  fprintf(stdout,"\nTotal classes %4d\n\n", nclasscheck);
			reset_class =0;
			reset_mut=0;

			// exit(0);
		}  // else selected by class

		// BOOSTRAP
		if (r==0) {

			if (run_mode==21 or run_mode==3)
			{
				if (verbose>1) fprintf(stdout,"                 ");
				int laa;
				if (run_mode==21) laa=1; else laa=0;
				for(int i=0;i<nvar-laa;i++)
					if (verbose>1) fprintf(stdout,"%c      ", aa[map[i]]);
				for(int i=0;i<laa;i++)
					if (verbose>1)  fprintf(stdout,"       " );
				if (verbose>1)  fprintf(stdout,"rmse   mae   ppv   sen   spe   pcc   mcc   ROC   PRC\n");


			} else

				if (verbose>1) 	fprintf(stdout,"                   A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y     rmse   mae   ppv   sen   spe  pcc   mcc   ROC   PRC\n");

		}




		// for (int kn = 0; kn < KN; kn++)
		//	fprintf(stderr,"kd %4d %4d %4d\n",kn,kn_start[kn],knbins[kn]);


		// K-fold loop
		int nclasscheck = 0;
		for (int kn = 0; kn < KN; kn++)
		{
			// reorder index by kn_start
			// copy randomized initial values
			for (int i = 0; i < d.nexp; i++)
			{
				d.index[i] =  d.index2[i];
			}

			// mov testing at the end
			for (int i = kn_start[kn]; i < kn_start[kn] + knbins[kn]; i++)
			{
				// Swapping a pair of values
				int j = d.nexp - knbins[kn] + i - kn_start[kn];
				int id =  d.index[j];
				d.index[j] =  d.index[i];
				d.index[i] = id;
				//fprintf(stderr,"i= %4d  j= %4d  %d\n",i,j,knbins[kn]);
			}
			//getchar();
			d.nrange = d.nexp - knbins[kn]; // Number of elements in TRAINING set

			// initialization

			for(int i=0;i<nvar;i++)
			{
				x[i] = xi[i]+distrFloat(gen);
				// x[i] = 0.3+distrFloat(gen);

				//fprintf(stdout,"%f ",x[i]);
				lb[i]=  xi[i] - wrange;

				ub[i]=  xi[i] + wrange;
				// ub[i]=  wrange/wrange*3.0 +0.00001;

			}

			for(int i=20;i<nvar;i++)
			{
				x[i] = 0.001 + distrFloat(gen)*0.5;

				// fprintf(stdout,"%f ",x[i]);
				lb[i]=  - wrange*0.5 -0.01;
				// if (lb[i]>-0.5) lb[i]=-0.6;

				ub[i]=   + wrange*0.5 + 0.01;

			}

			//
			//			x[1]=0.3;
			//			lb[1]=x[1]-0.01;
			//			ub[1]=x[1]+0.01;
			//
			//			x[18]=0.3;
			//			lb[18]=x[18]-0.01;
			//			ub[18]=x[18]+0.01;

			//			x[19]=2.3;
			//			lb[19]=x[19]-0.01;
			//			ub[19]=x[19]+0.01;

			//			x[4]=1.0;
			//			lb[19]=x[4]-0.01;
			//			ub[19]=x[4]+0.01;



			//			x[20]=0.2;
			//			lb[20]=x[20]-0.01;
			//			ub[20]=x[20]+0.01;


			//			x[5]=0.6;
			//			lb[5]=x[5]-0.01;
			//			ub[5]=x[5]+0.01;
			//
			//			x[2]=0.6;
			//			lb[2]=x[2]-0.01;
			//			ub[2]=x[2]+0.01;
			//
			//			x[3]=0.6;
			//		    lb[3]=x[3]-0.01;
			//			ub[3]=x[3]+0.01;


			//			x[11]=0.900;
			//			lb[11]=x[11]-0.01;
			//			ub[11]=x[11]+0.01;
			//
			//			x[16]=1.30;
			//			lb[16]=x[16]-0.01;
			//			ub[16]=x[16]+0.01;

			//			x[10]=0.900;
			//			lb[10]=x[10]-0.01;
			//			ub[10]=x[10]+0.01;
			//
			//			x[14]=1.30;
			//			lb[14]=x[14]-0.01;
			//			ub[14]=x[14]+0.01;
			//

			//			x[12]=0.900;
			//			lb[12]=x[12]-0.01;
			//			ub[12]=x[12]+0.01;

			//			x[18]=1.0;
			//			lb[18]=x[18]-0.01;
			//			ub[18]=x[18]+0.01;


			if (verbose>1) fprintf(stdout,"%-4d",kn);
			fprintf(fout,"%-5d %-4d",r, kn);


			// Set optimiztion function stuff
			nlopt_opt func;
			func = nlopt_create(NLOPT_LN_BOBYQA, nvar);
			// func = nlopt_create(NLOPT_LN_AUGLAG, nvar); // set number of optimizaiton variables
			//func = nlopt_create(NLOPT_LN_SBPLX, nvar);
			// func = nlopt_create(NLOPT_LN_NEWUOA, nvar);

			nlopt_set_lower_bounds(func,lb);
			nlopt_set_upper_bounds(func,ub);
			// nlopt_set_ftol_rel(func, 1e-7); // End tolerance
			// nlopt_set_ftol_rel(func, 1e-1); // End tolerance?
			nlopt_set_ftol_abs(func, 1e-8); // End tolerance?
			nlopt_set_maxeval(func, 50000);

			double dumpx, dumpy, sumax=0, sumay=0, suma=0;
			//		for(int i=0;i<d.nexp;i++)
			//		{
			//			dumpx = d.com[i]-d.aveC;
			//			dumpy = d.exp[i]-d.aveE;
			//			suma += dumpx * dumpy;
			//			sumax += dumpx * dumpx;
			//			sumay += dumpy * dumpy;
			//			// fprintf(stdout, "%d %f %f %f %f %f\n", i, dumpx, dumpy, sumax, sumay, suma);
			//		}

			// double corre;
			//		corre = suma / ( sqrt(sumax) * sqrt(sumay) );
			// fprintf(stdout," corre %f\n", corre);

			float rmse; // Root Mean Squared Error
			float rmse_abs; // Mean Absolute Error


			switch(run_mode)
			{
			case 3:
				nlopt_set_min_objective(func, MAE_dic, &d);
				//nlopt_set_min_objective(func, HUBERL_dic, &d);

				// nlopt_set_max_objective(func, AUC_dic, &d);


				break;

			case 4:
				// nlopt_set_min_objective(func, TEST0H_dic, &d);
				nlopt_set_min_objective(func, TEST0_dic, &d);

				break;

			case 5:
				// nlopt_set_max_objective(func, PCC, &d);
				// nlopt_set_min_objective(func, PPV, &d);
				// nlopt_set_min_objective(func, RMSE_1, &d); // set nvar=1
				// nlopt_set_min_objective(func, RMSE_2, &d); // set nvar=1
				// nlopt_set_min_objective(func, RMSE, &d);
				// nlopt_set_min_objective(func, PPV, &d);
				// nlopt_set_min_objective(func, HUBERL, &d);
				nlopt_set_min_objective(func, MAE, &d);

				break;

			case 6:
				nlopt_set_min_objective(func, MAErama, &d);
				break;

			case 56:
				// nlopt_set_min_objective(func, PPV56, &d); // nvar =40 pending
				nlopt_set_min_objective(func, MAE21rama, &d);
				break;

			case 400:
				nlopt_set_min_objective(func, MAE400, &d);
				break;

			case 21:

				//nlopt_set_min_objective(func, MAE21, &d);
				nlopt_set_min_objective(func, MAE21_dic, &d);

				break;

			case 22:
				// nlopt_set_min_objective(func, MAERSA, &d);
				nlopt_set_min_objective(func, MAE22, &d);

				break;

			case 40:
				// nlopt_set_min_objective(func, MAERSA40, &d);
				nlopt_set_min_objective(func, MAE40, &d);

				break;

			default:
				fprintf(stderr,"Error, unknown run_mode %d, forcing exit!\n", run_mode);
				exit(1);
				break;

			}


			// Performing the optimization
			if( nlopt_optimize(func, x, &f0) < 0 )
			{
				printf("%d nlopt failed! %f\n", r, f0);
				getchar();
			}
			else
			{
				// display res
				// Compute KORPM energies in TEST set
				double ave = 0; int aa1, wt, mut;
				double dumpx = 0, sumaa=0, sumab=0;

				for(int ii=d.nrange; ii<d.nexp; ii++) // Screen test set mutations
				{
					int i = d.index[ii];

					dumpx = 0, sumaa=0, sumab=0;

					switch(run_mode)
					{
					case 3:


						wt = map [d.aa1[i]];
						mut = map [d.aa2[i]];
						for(int j=0;j<20;j++)
						{
							aa1 = map[j];
							sumaa += x[wt]  * x[aa1] * ( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
							sumab += x[mut] * x[aa1] * ( d.dcom[i][j*2] + d.dcom[i][2*j+1] );
						}
						dumpx=sumab-sumaa;
						d.com[i]=dumpx/100;
						break;



					case 4:  // simple weight

						if (nvar != 1 )
						{
							wt = map [d.aa1[i]];
							mut = map [d.aa2[i]];
							d.com[i] = x[ mut] - x[ wt ];
						}
						else
							d.com[i]=x[0]; //100;						break;
						break;

					case 5:
						wt =d.aa1[i];
						mut = d.aa2[i];

						double cte, cte2;
						for(int j=0;j<20;j++)
						{
							cte=0; cte2=0;

							//							if (((wt == 2)||(wt == 3)) && ((j == 8)||(j == 14)||(j == 6)))  cte =-.1;
							//						    if (((mut == 2)||(mut == 3)) && ((j == 8)||(j == 14)||(j == 6)))  cte2 =-.1;

							//										if (((j == 2)||(j == 3)) && ((wt == 8)||(wt == 14)|(wt == 6)))      cte =-0.1;
							//									    if (((j == 2)||(j == 3)) && ((mut == 8)||(mut == 14)||(mut == 6)))  cte2 =-0.1;

							sumaa += x[wt] * x[j] * ( d.dexp[i][j*2] + d.dexp[i][2*j+1]  );
							sumab += x[mut] * x[j] * ( d.dcom[i][j*2] + d.dcom[i][2*j+1] );

							//sumaa += 2/(1+exp(-0.5*x[wt] * x[j])) * ( d.dexp[i][j*3] + d.dexp[i][2*j+1]  );
							//sumab += 2/(1+exp(-0.5*x[mut] * x[j])) * ( d.dcom[i][j*3] + d.dcom[i][2*j+1] );

							//sumaa+=x[0]*(d.dexp[i][j*2]+d.dexp[i][2*j+1])+x[1]; // RMSE_2
							//sumab+=x[0]*(d.dcom[i][j*2]+d.dcom[i][2*j+1])+x[1]; // RMSE_2
						}
						dumpx=sumab-sumaa;
						d.com[i]=dumpx/100;
						// d.com[i] *= .5*(Score[bi[wt]][bi[wt]] - Score[bi[wt]][bi[mut]]);
						//						if ((( mut== 2)||( mut == 3)) && ((wt == 8)||(wt == 14)|(wt == 6)))   d.com[i] -= 0.5;
						//					    if ((( wt== 2)||( wt == 3)) && ((mut == 8)||(mut == 14)|(mut == 6)))  d.com[i] -= 0.5;


						break;

					case 6:

						for(int j=1; j<2; j++)
						{
							sumaa += x[ 0 ] *  ( d.EramaWT[i][j] ); // Wild-Type
							sumab += x[ 0 ]  * ( d.EramaMut[i][j] ); // Mutant
							// sumaa += x[ d.iaaWT[i][j] ] * ( d.EramaWT[i][j] ); // Wild-Type
							// sumab += x[ d.iaaMut[i][j] ] * ( d.EramaMut[i][j] ); // Mutant

						}
						dumpx=sumab-sumaa;
						d.com[i] = dumpx;
						break;



					case 56:
						wt = map [d.aa1[i]];
						mut = map [d.aa2[i]];
						for(int j=0;j<20;j++)
						{
							sumaa += x[wt] * x[j] * ( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
							sumab += x[mut] * x[j] * ( d.dcom[i][j*2] + d.dcom[i][2*j+1] );
						}


						for(int j=0; j<3; j++)
						{
							// sumaa += x[ d.iaaWT[i][j] + 20 ] * ( d.EramaWT[i][j] ); // Wild-Type
							// sumab += x[ d.iaaMut[i][j] +20 ] * ( d.EramaMut[i][j] ); // Mutant
							sumaa += 100*x[  20 ] * ( d.EramaWT[i][j] ); // Wild-Type
							sumab += 100*x[  20 ] * ( d.EramaMut[i][j] ); // Mutant

						}
						dumpx=sumab-sumaa;
						d.com[i] = dumpx/100;
						break;

					case 400:
						for(int j=0;j<20;j++)
						{
							sumaa += x[d.aa1[i] + 20 * j] *( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
							sumab += x[d.aa2[i] + 20 * j] *( d.dcom[i][j*2] + d.dcom[i][2*j+1] );
						}
						dumpx = sumab - sumaa;
						d.com[i] = dumpx / 100;
						break;

					case 21:


						wt = map [d.aa1[i]];
						mut = map [d.aa2[i]];
						for(int j=0;j<20;j++)
						{
							aa1 = map[j];
							sumaa += x[wt]  * x[aa1] * ( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
							sumab += x[mut] * x[aa1] * ( d.dcom[i][j*2] + d.dcom[i][2*j+1] );
						}
						dumpx=sumab-sumaa;
						d.com[i]= dumpx/100 +  x[nvar-1];
						break;


					case 22:
						wt =d.aa1[i];
						mut = d.aa2[i];

						for(int j=0;j<20;j++)
						{
							sumaa += x[wt] * x[j] * ( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
							sumab += x[mut] * x[j] * ( d.dcom[i][j*2] + d.dcom[i][2*j+1] );



						}
						dumpx = (sumab - sumaa)/100;

						//						 if ((hyb[wt] == 0 ) and (hyb[mut] == 1 ))
						//							d.com[i] = x[nvar-2] + dumpx ; // Add some constant scale value
						//						else if ((hyb[wt] == 1 ) and (hyb[wt] == 0 ))
						//							d.com[i] = -x[nvar-2] + dumpx;
						//						else  d.com[i] = x[nvar-1] + dumpx ;
						//
						//						 if ((pnn[wt] == 0 ) and (pnn[mut] == 1 ))
						//						//if (d.rsa[i] > 0.2)
						//							d.com[i] = x[nvar-2] + dumpx ; // Add some constant scale value
						//						else if ((pnn[wt] == 1 ) and (pnn[wt] == 0 ))
						//							d.com[i] = -x[nvar-2] + dumpx;
						//						else  d.com[i] = x[nvar-1] + dumpx ;

						//						if (d.rsa[i]>0.3)
						//						d.com[i] = (1-x[nvar-2]*d.rsa[i])*dumpx;
						//						else
						//						d.com[i] = (1-x[nvar-1]*d.rsa[i])*dumpx;

						//						if ( pnn[wt] == 0 )
						//							d.com[i] = (1.0-x[20]*d.rsa[i])*dumpx;
						//						else if ( pnn[wt] == 1 )
						//							 d.com[i] = (1.0+x[20]*d.rsa[i])*dumpx;
						//						else 	 d.com[i] = (1.0-x[21]*d.rsa[i])*dumpx;

						//d.com[i] = 0.9/(1+exp(-x[20]*(d.rsa[i]-x[21]) ) )*dumpx;

						if ( hyb[wt] <= 1 )
							d.com[i] = dumpx + (dol[ wt ] - dol[ mut ])/60.1*x[20] + x[21] ;
						else
							d.com[i] = dumpx  + x[21];

						// d.com[i] = dumpx / 100; // Add some constant scale value
						break;

					case 39:
						if (d.rsa[i] > 20.0)
							for(int j=0;j<20;j++)
							{
								sumaa += x[d.aa1[i]] * x[j] * ( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
								sumab += x[d.aa2[i]] * x[j] * ( d.dcom[i][j*2] + d.dcom[i][2*j+1] );
							}
						else
							for(int j=0;j<20;j++)
							{
								sumaa += x[d.aa1[i] + 20] * x[j + 20] * ( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
								sumab += x[d.aa2[i] + 20] * x[j + 20] * ( d.dcom[i][j*2] + d.dcom[i][2*j+1] );
							}
						dumpx = sumab - sumaa;

						d.com[i] = dumpx / 100; // Add some constant scale value
						break;

					case 40:
						wt = d.aa1[i];
						mut = d.aa2[i];
						double ave=0, v, expf;
						for(int j=0;j<20;j++)
						{
							//	sumaa += x[wt] * x[j] * ( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
							//	sumab += x[mut] * x[j] * ( d.dcom[i][j*2] + d.dcom[i][2*j+1] );

							v=( d.dexp[i][j*2] + d.dexp[i][2*j+1] );
							sumaa += x[ wt ]* x[j] /(1+ exp(-1.0*(v-x[j+20])))/(1+ exp(-1.0*(v-x[wt+20])));
							//	sumaa += x[ wt ]* x[j] /(1+ exp(-v*fabs(x[j+20])))/(1+ exp(-v*fabs(x[wt+20])));

							v=( d.dcom[i][j*2] + d.dcom[i][2*j+1] );
							sumab += x[ mut ]*x[j] /(1+ exp(-1.0*(v-x[j+20])))/(1+ exp(-1.0*(v-x[mut+20])));
							//	sumab += x[ mut ]* x[j] /(1+ exp(-v*fabs(x[j+20])))/(1+ exp(-v*fabs(x[mut+20])));

							//	ave += (x[ j+ 20]);
						}
						// ave/=20.0;
						dumpx = sumab - sumaa;

						//dumpx = (sumab - 100*fabs(x[ mut+ 20])-ave) - (sumaa - 100*fabs(x[ wt + 20]-ave)); // Mutation ddG (check sign)
						//dumpx = (sumab - 100*(x[ mut+ 20])-ave) - (sumaa - 100*(x[ wt + 20]-ave)); // Mutation ddG (check sign)

						// dumpx = (sumab - 100*(x[ mut+ 20])) - (sumaa - 100*(x[ wt + 20])); // Mutation ddG (check sign)
						// dumpx = (sumab + 100*fabs(x[ mut+ 20])) - (sumaa + 100*fabs(x[ wt + 20])); // Mutation ddG (check sign)

						d.com[i] = dumpx / 100; // Add some constant scale value
						break;

					default:
						fprintf(stderr,"Error, unknown run_mode %d, forcing exit!\n", run_mode);
						exit(1);
						break;
					}

					ave += d.com[i];
				}

				int TPddG = 0;
				int FPddG = 0;
				int FNddG = 0;
				int TNddG = 0;
				float thr = 0.0;

				// thr =0;
				suma = 0.0;
				double suma_sqrt = 0.0;

				// minP=1000;
				// maxP=-1000;
				for(int ii=d.nrange; ii<d.nexp; ii++)
				{
					int i=d.index[ii];

					dumpx= (d.com[i]-d.exp[i]); // E_mut - E_wt



					if (save_switch)
						fprintf(fout2,"%c99%c %7.3f %7.3f\n", aa[d.aa1[i]], aa[d.aa2[i]], d.exp[i], d.com[i] );

					suma_sqrt += dumpx*dumpx;
					suma += fabsf(dumpx);

					// if(d.com[i] < minP) minP=d.com[i];
					// if(d.com[i] > maxP) maxP=d.com[i];


					// Contingency table (confusion matrix)
					if(d.exp[i] >= 0.0) // Stabilizing mutation, ddG >=0 (Positive)
					{
						if(d.com[i] >= thr )
							TPddG++; // Count True Positives
						else
							FNddG++; // Count False Negatives

					}
					else // Destabilizing mutation, ddG < 0 (Negative)
					{
						if(d.com[i] < thr)
							TNddG++; // Count True Negatives
						else
							FPddG++; // Count False Positives
					}
				}

				// float rmse = sqrt(suma/nrange);
				rmse = sqrt(suma_sqrt/(d.nexp-d.nrange)); // Root Mean Squared Error
				rmse_abs = (suma/(d.nexp-d.nrange)); // Mean Absolute Error
				// fprintf(stderr, "suma_sqrt= %f  suma= %f  d.nexp= %d  nrange= %d\n", suma_sqrt, suma, d.nexp, nrange);

				//f0avg += f0;
				f0avg += rmse;
				ave /= (d.nexp-d.nrange);

				NEFF += (d.nexp-d.nrange);

				if(TPddG > 0)
				{
					ppv = (float) TPddG / ((float) TPddG + (float) FPddG); // Positive Predictive Value = PRECision
					sen = (float) TPddG / ((float) TPddG + (float) FNddG);  // Sensitivity (SN) = RECall = True Positive Rate (TPR)
				}
				else
				{
					ppv = 0.0;
					sen = 0.0;
				}

				if( FPddG != 0 )
					fpr = (float) FPddG / (float) (TNddG + FPddG); // False Positive Rate (FPR)
				else
					fpr = 0.0;

				ppvavg += ppv;

				if( TPddG + TNddG > 0 )
					acc = (float) (TPddG + TNddG) / (float) (TPddG + TNddG + FNddG + FPddG); // ACCuracy
				else
					acc = 0.0;

				/*
				if( FPddG + FNddG > 0 )
					err = (float) (FPddG + FNddG) / (float) (TPddG + TNddG + FNddG + FPddG); // ERror Rate
				else
					err = 0.0;
				 */

				if( TNddG + FNddG > 0 )
					err = (float) (TNddG) / (float) ( FNddG + TNddG); // now NPV
				else
					err = 0.0;

				if( (TPddG + FPddG)*(TPddG + FNddG)*(TNddG + FPddG)*(TNddG + FNddG) != 0 )
					mcc = (float) (TPddG * TNddG - FPddG * FNddG) / ( sqrt( (float) (TPddG + FPddG)*(TPddG + FNddG)*(TNddG + FPddG)*(TNddG + FNddG)) ); // Matthews correlation coefficient
				else
					mcc = 0.0;

				// AUC

				float deltau=  (maxP - minP)/BINSAUC;
				for(int au=0; au<BINSAUC; au++)
				{
					TPddG = 0;
					FPddG = 0;
					FNddG = 0;
					TNddG = 0;
					thr = minP + au * deltau;
					for(int ii=d.nrange; ii<d.nexp; ii++)
					{
						int i=d.index[ii];


						// Contingency table (confusion matrix)
						if(d.exp[i] >= 0.0) // Stabilizing mutation, ddG >=0 (Positive)
						{
							if(d.com[i] >= thr )
								TPddG++; // Count True Positives
							else
								FNddG++; // Count False Negatives

						}
						else // Destabilizing mutation, ddG < 0 (Negative)
						{
							if(d.com[i] < thr)
								TNddG++; // Count True Negatives
							else
								FPddG++; // Count False Positives
						}
					}

					if( TPddG + FNddG > 0 )
						arec[au]=(float) TPddG / ((float) TPddG + (float) FNddG);  // RECall = Sensitivity (SN) = True Positive Rate (TPR)
					else
						arec[au]=0;

					if( TPddG + FPddG > 0 )
						appv[au]=(float) TPddG / ((float) TPddG + (float) FPddG); // Positive Predictive Value = PRECision
					else
						appv[au]=0;

					if( TNddG + FPddG > 0 )
						afpr[au]=(float) FPddG / (float) (TNddG + FPddG); // False Positive Rate (FPR)
					else
						afpr[au]=0;

					Afpr[au][KN*r+kn]=	afpr[au];
					Arec[au][KN*r+kn]=	arec[au];
					Appv[au][KN*r+kn]=	appv[au];


				}
				float AUC_ROC = 0.0; // Area Under Curve for ROC
				float AUC_PRC = 0.0; // Area Under Curve for PRC
				for(int au=0; au<BINSAUC-1; au++)
				{
					//fprintf(stderr,"%7.3f %7.3f %7.3f %7.3f\n", afpr[au+1], arec[au+1], appv[au+1],  (minP + (au+1) * deltau) );

					AUC_ROC += (afpr[au] - afpr[au+1]) * (arec[au] + arec[au+1]) / 2;
					AUC_PRC += (arec[au] - arec[au+1]) * (appv[au] + appv[au+1]) / 2;
				}


				//  cross-correlation
				suma = 0.0;
				sumax = 0.0;
				sumay = 0.0;
				for(int ii=d.nrange; ii<d.nexp; ii++)
				{
					int i=d.index[ii];
					dumpx = d.com[i] - ave;
					dumpy = d.exp[i] - d.aveE;
					suma += dumpx*dumpy;
					sumax += dumpx*dumpx;
					sumay += dumpy*dumpy;
				}
				double corre;

				if( sumax*sumay != 0.0 )
					corre = suma/(sqrt(sumax)*sqrt(sumay));
				else
					corre = 0.0;

				f0correavg += corre;


				// Keep maximum
				if( ppv > ppvmax )
				{
					ppvmax = ppv;
					//f0min = f0;
					f0min = rmse;
					f0corre = corre;
					for(int i=0;i<nvar;i++)
						xmax[i] = x[i];
				}

				// Weights stuff
				if (verbose>1) fprintf(stdout," %7.5f ", f0);
				fprintf(fout," %7.5f ", f0);

				if (nvar<=20)
				{
					for(int i=0;i<nvar;i++)
					{
						if (verbose>1) fprintf(stdout,"%7.3f", x[i]);
						fprintf(fout,"%7.3f", x[i]);

						xavg[i] += x[i];
						if (avg_switch)
							Weights[i][r] = x[i]; // Store Weights into the corresponding array
						else
							Weights[i][KN*r+kn] = x[i]; // Store Weights into the corresponding array

					}
				} else
				{

					for(int i=0;i<20;i++)
					{
						if (verbose>1) fprintf(stdout,"%7.3f", x[i]);
						fprintf(fout,"%7.3f", x[i]);

						xavg[i] += x[i];
						if (avg_switch)
							Weights[i][r] = x[i]; // Store Weights into the corresponding array
						else
							Weights[i][KN*r+kn] = x[i]; // Store Weights into the corresponding array

					}

					if(nvar != 21 and nvar!=22)
					{
						if (verbose>1) fprintf(stdout,"\n    ");
						if (verbose>1) fprintf(stdout," %7s ", "");
					}
					// fprintf(stdout," %7.5f ", f0);
					for(int i=20;i<nvar;i++)
					{
						if (verbose>1) fprintf(stdout,"%7.3f", x[i]);
						fprintf(fout,"%7.3f", x[i]);

						xavg[i] += x[i];
						if (avg_switch)
							Weights[i][r] = x[i]; // Store Weights into the corresponding array
						else
							Weights[i][KN*r+kn] = x[i]; // Store Weights into the corresponding array

					}
				}

				// Dump some info to screen
				if (verbose>1) fprintf(stdout," -> %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", rmse, rmse_abs, ppv, sen, 1-fpr, corre, mcc, AUC_ROC, AUC_PRC);
				fprintf(fout," %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", rmse, rmse_abs, ppv, sen, 1-fpr,corre, mcc, AUC_ROC, AUC_PRC);

				//fprintf(stderr,"kd %4d %4d %4d %4d %4d->", kn, kn_start[kn], knbins[kn], kn_start[kn] + knbins[kn], classcountN[kn]);

				if (pclase)
				{
					for (int i = nclasscheck; i < nclasscheck+classcountN[kn]; i++)
					{
						classcountT[classindex[i]]++;
						classcountF[classindex[i]] +=AUC_PRC;
						//	fprintf(stderr,"%4d ", classindex[i]);

					}
					//fprintf(stderr,"\n");
					nclasscheck += classcountN[kn];
				}


				// Store stuff to compute general statistics

				if (avg_switch) {
					Rmse[r] += rmse;
					Mae[r]  += rmse_abs; // Mean Absolute Error array
					Ppv[r]  += ppv; // Positive Predictive Power array
					Sen[r]  += sen; // Sensitivity array
					Pcor[r] += corre; // Pearson's cross-correlation
					Acc[r]  += acc; // ACCuracy
					Err[r]  += err; // ERror Rate
					Mcc[r]  += mcc; // Matthews correlation coefficient
					Fpr[r]  += fpr; // False Positive Rate (FPR)
					auc_roc[r]  += AUC_ROC;
					auc_prc[r]  += AUC_PRC;
				} else
				{
					Rmse[KN*r+kn] = rmse; // Root Mean Squared Error array
					Mae[KN*r+kn] = rmse_abs; // Mean Absolute Error array
					Ppv[KN*r+kn] = ppv; // Positive Predictive Power array
					Sen[KN*r+kn] = sen; // Sensitivity array
					Pcor[KN*r+kn] = corre; // Pearson's cross-correlation
					Acc[KN*r+kn] = acc; // ACCuracy
					Err[KN*r+kn] = err; // ERror Rate
					Mcc[KN*r+kn] = mcc; // Matthews correlation coefficient
					Fpr[KN*r+kn] = fpr; // False Positive Rate (FPR)
					Rmse[KN*r+kn] = rmse; // Root Mean Squared Error array
					Rmse[KN*r+kn] = rmse; // Root Mean Squared Error array
					auc_roc[KN*r+kn]  += AUC_ROC;
					auc_prc[KN*r+kn]  += AUC_PRC;
				}

			}


		} // kn
		if (avg_switch) {
			Rmse[r] /= KN;
			Mae[r]  /= KN;
			Ppv[r]  /= KN;
			Pcor[r] /= KN;
			Sen[r]  /= KN; // Sensitivity array
			Acc[r]  /= KN;
			Err[r]  /= KN;
			Mcc[r]  /= KN;
			Fpr[r]  /= KN;
			auc_roc[r]  /= KN;
			auc_prc[r]  /= KN;
		}

	} // Repeat

	fclose(fout);
	if (save_switch)
		fclose(fout2);

	if (verbose>1) fprintf(stdout," %7.5f %7.5f %7.5f->", f0min, ppvmax, f0corre );
	int laa=20;

	if (run_mode<=4) nvar=20;
	else {

		if (run_mode==21) {
			map[20]=nvar-1;
			nvar=21;
			laa=21;
		}
		else {		// no mapping or 21

			for(int i=0;i<20;i++)
				map[i]=i;
		}
	}
	for(int i=0;i<laa;i++)
	{
		if (verbose>1) fprintf(stdout,"%7.5f ", xmax[map[i]]);
	}
	for(int i=laa;i<nvar;i++)
	{
		if (verbose>1) fprintf(stdout,"%7.5f ", xmax[i]);
	}
	if (verbose>1) fprintf(stdout,"\n");

	float fact=(KNNRUNS*1.0);
	if (verbose>1) fprintf(stdout," %7.5f %7.5f %7.5f->", f0avg /fact, ppvavg/fact, f0correavg/fact);
	for(int i=0;i<laa;i++)
	{
		if (verbose>1) fprintf(stdout,"%7.5f ", xavg[map[i]]/fact);
	}
	for(int i=laa;i<nvar;i++)
	{
		if (verbose>1) fprintf(stdout,"%7.5f ", xavg[i]/fact);
	}

	if (verbose>0) {
		fprintf(stdout,"\n");

		if (nvar!=40) {

			// Statistics output
			fprintf(stdout,"\n%-3s","  ");
			for(int i=0;i<nvar;i++)
				fprintf(stdout," %6c", aa[i] );
			fprintf(stdout," | %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s\n","RMSE","MAE","PPV","Sen","Spe", "PCC","ACC","NPV","MCC","ROC","PRC" );

			fprintf(stdout,"%-3s", "AVG");
			for(int i=0;i<laa;i++)
				fprintf(stdout," %6.3f", average(Weights[map[i]],ENRUNS) );
			for(int i=laa;i<nvar;i++)
				fprintf(stdout," %6.3f", average(Weights[i],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", average(Rmse,ENRUNS), average(Mae,ENRUNS), average(Ppv,ENRUNS), average(Sen,ENRUNS), (1-average(Fpr,ENRUNS)), average(Pcor,ENRUNS), average(Acc,ENRUNS), average(Err,ENRUNS), average(Mcc,ENRUNS), average(auc_roc,ENRUNS), average(auc_prc,ENRUNS) );

			fprintf(stdout,"%-3s", "SIG");
			for(int i=0;i<laa;i++)
				fprintf(stdout," %6.3f", sigma(Weights[map[i]],ENRUNS) );
			for(int i=laa;i<nvar;i++)
				fprintf(stdout," %6.3f", sigma(Weights[i],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", sigma(Rmse,ENRUNS), sigma(Mae,ENRUNS), sigma(Ppv,ENRUNS), sigma(Sen,ENRUNS), sigma(Fpr,ENRUNS), sigma(Pcor,ENRUNS), sigma(Acc,ENRUNS), sigma(Err,ENRUNS), sigma(Mcc,ENRUNS), sigma(auc_roc,ENRUNS), sigma(auc_prc,ENRUNS) );

			fprintf(stdout,"%-3s", "MED");
			for(int i=0;i<laa;i++)
				fprintf(stdout," %6.3f", median(Weights[map[i]],ENRUNS) );
			for(int i=laa;i<nvar;i++)
				fprintf(stdout," %6.3f", median(Weights[i],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", median(Rmse,ENRUNS), median(Mae,ENRUNS), median(Ppv,ENRUNS), median(Sen,ENRUNS), (1-median(Fpr,ENRUNS)), median(Pcor,ENRUNS), median(Acc,ENRUNS), median(Err,ENRUNS), median(Mcc,ENRUNS), median(auc_roc,ENRUNS), median(auc_prc,ENRUNS) );

			fprintf(stdout,"%-3s", "MAX");
			for(int i=0;i<laa;i++)
				fprintf(stdout," %6.3f", getmax(Weights[map[i]],ENRUNS) );
			for(int i=laa;i<nvar;i++)
				fprintf(stdout," %6.3f", getmax(Weights[i],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", getmax(Rmse,ENRUNS), getmax(Mae,ENRUNS), getmax(Ppv,ENRUNS), getmax(Sen,ENRUNS), (1-getmax(Fpr,ENRUNS)), getmax(Pcor,ENRUNS), getmax(Acc,ENRUNS), getmax(Err,ENRUNS), getmax(Mcc,ENRUNS), getmax(auc_roc,ENRUNS), getmax(auc_prc,ENRUNS));

			fprintf(stdout,"%-3s", "MIN");
			for(int i=0;i<laa;i++)
				fprintf(stdout," %6.3f", getmin(Weights[map[i]],ENRUNS) );
			for(int i=laa;i<nvar;i++)
				fprintf(stdout," %6.3f", getmin(Weights[i],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", getmin(Rmse,ENRUNS), getmin(Mae,ENRUNS), getmin(Ppv,ENRUNS), getmin(Sen,ENRUNS), (1-getmin(Fpr,ENRUNS)) , getmin(Pcor,ENRUNS), getmin(Acc,ENRUNS), getmin(Err,ENRUNS), getmin(Mcc,ENRUNS),  getmin(auc_roc,ENRUNS), getmin(auc_prc,ENRUNS));

		} else
		{
			// Statistics output
			if (nvar < 22) fprintf(stdout,"\n%-3s","      ");
			else fprintf(stdout,"%-3s", "  ");
			for(int i=0;i<20;i++)
				fprintf(stdout," %6c", aa[i] );
			fprintf(stdout," | %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s\n","RMSE","MAE","PPV","Sen","FPR", "PCC","ACC","NPV","MCC","ROC","PRC" );

			if (nvar < 22)  fprintf(stdout,"%-3s", "  ");
			fprintf(stdout,"%-3s", "AVG");
			for(int i=0;i<20;i++)
				fprintf(stdout," %6.3f", average(Weights[map[i]],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", average(Rmse,ENRUNS), average(Mae,ENRUNS), average(Ppv,ENRUNS), average(Sen,ENRUNS), (1-average(Fpr,ENRUNS)), average(Pcor,ENRUNS), average(Acc,ENRUNS), average(Err,ENRUNS), average(Mcc,ENRUNS), average(auc_roc,ENRUNS), average(auc_prc,ENRUNS) );

			fprintf(stdout,"%-3s", "  ");
			for(int i=20;i<nvar;i++)
				fprintf(stdout," %6.3f", average(Weights[i],ENRUNS) );
			if (nvar > 22)
				fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", average(Rmse,ENRUNS), average(Mae,ENRUNS), average(Ppv,ENRUNS), average(Sen,ENRUNS), (1-average(Fpr,ENRUNS)), average(Pcor,ENRUNS), average(Acc,ENRUNS), average(Err,ENRUNS), average(Mcc,ENRUNS), average(auc_roc,ENRUNS), average(auc_prc,ENRUNS) );



			fprintf(stdout,"%-3s", "SIG");
			for(int i=0;i<20;i++)
				fprintf(stdout," %6.3f", sigma(Weights[map[i]],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", sigma(Rmse,ENRUNS), sigma(Mae,ENRUNS), sigma(Ppv,ENRUNS), sigma(Sen,ENRUNS), sigma(Fpr,ENRUNS), sigma(Pcor,ENRUNS), sigma(Acc,ENRUNS), sigma(Err,ENRUNS), sigma(Mcc,ENRUNS), sigma(auc_roc,ENRUNS), sigma(auc_prc,ENRUNS) );
			fprintf(stdout,"%-3s", "  ");
			for(int i=20;i<nvar;i++)
				fprintf(stdout," %6.3f", sigma(Weights[i],ENRUNS) );
			if (nvar > 22)
				fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", sigma(Rmse,ENRUNS), sigma(Mae,ENRUNS), sigma(Ppv,ENRUNS), sigma(Sen,ENRUNS), sigma(Fpr,ENRUNS), sigma(Pcor,ENRUNS), sigma(Acc,ENRUNS), sigma(Err,ENRUNS), sigma(Mcc,ENRUNS), sigma(auc_roc,ENRUNS), sigma(auc_prc,ENRUNS) );


			fprintf(stdout,"%-3s", "MED");
			for(int i=0;i<20;i++)
				fprintf(stdout," %6.3f", median(Weights[map[i]],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", median(Rmse,ENRUNS), median(Mae,ENRUNS), median(Ppv,ENRUNS), median(Sen,ENRUNS), (1-median(Fpr,ENRUNS)), median(Pcor,ENRUNS), median(Acc,ENRUNS), median(Err,ENRUNS), median(Mcc,ENRUNS), median(auc_roc,ENRUNS), median(auc_prc,ENRUNS) );
			fprintf(stdout,"%-3s", "  ");
			for(int i=20;i<nvar;i++)
				fprintf(stdout," %6.3f", median(Weights[i],ENRUNS) );
			if (nvar > 22)
				fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", median(Rmse,ENRUNS), median(Mae,ENRUNS), median(Ppv,ENRUNS), median(Sen,ENRUNS), (1-median(Fpr,ENRUNS)), median(Pcor,ENRUNS), median(Acc,ENRUNS), median(Err,ENRUNS), median(Mcc,ENRUNS), median(auc_roc,ENRUNS), median(auc_prc,ENRUNS) );


			fprintf(stdout,"%-3s", "MAX");
			for(int i=0;i<20;i++)
				fprintf(stdout," %6.3f", getmax(Weights[map[i]],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", getmax(Rmse,ENRUNS), getmax(Mae,ENRUNS), getmax(Ppv,ENRUNS), getmax(Sen,ENRUNS), (1-getmax(Fpr,ENRUNS)), getmax(Pcor,ENRUNS), getmax(Acc,ENRUNS), getmax(Err,ENRUNS), getmax(Mcc,ENRUNS), getmax(auc_roc,ENRUNS), getmax(auc_prc,ENRUNS));
			fprintf(stdout,"%-3s", "  ");
			for(int i=20;i<nvar;i++)
				fprintf(stdout," %6.3f", getmax(Weights[i],ENRUNS) );
			if (nvar > 22)
				fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", getmax(Rmse,ENRUNS), getmax(Mae,ENRUNS), getmax(Ppv,ENRUNS), getmax(Sen,ENRUNS), (1-getmax(Fpr,ENRUNS)), getmax(Pcor,ENRUNS), getmax(Acc,ENRUNS), getmax(Err,ENRUNS), getmax(Mcc,ENRUNS), getmax(auc_roc,ENRUNS), getmax(auc_prc,ENRUNS));

			fprintf(stdout,"%-3s", "MIN");
			for(int i=0;i<20;i++)
				fprintf(stdout," %6.3f", getmin(Weights[map[i]],ENRUNS) );
			fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", getmin(Rmse,ENRUNS), getmin(Mae,ENRUNS), getmin(Ppv,ENRUNS), getmin(Sen,ENRUNS), (1-getmin(Fpr,ENRUNS)) , getmin(Pcor,ENRUNS), getmin(Acc,ENRUNS), getmin(Err,ENRUNS), getmin(Mcc,ENRUNS),  getmin(auc_roc,ENRUNS), getmin(auc_prc,ENRUNS));
			fprintf(stdout,"%-3s", "  ");
			for(int i=20;i<nvar;i++)
				fprintf(stdout," %6.3f", getmin(Weights[i],ENRUNS) );
			if (nvar > 22)
				fprintf(stdout," | %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", getmin(Rmse,ENRUNS), getmin(Mae,ENRUNS), getmin(Ppv,ENRUNS), getmin(Sen,ENRUNS), (1-getmin(Fpr,ENRUNS)) , getmin(Pcor,ENRUNS), getmin(Acc,ENRUNS), getmin(Err,ENRUNS), getmin(Mcc,ENRUNS),  getmin(auc_roc,ENRUNS), getmin(auc_prc,ENRUNS));



		}

	}


	{

		if (verbose>0) fprintf(stdout," save rocs.txt %d\n", NRUNS);

		if ( (f=fopen("rocs.txt", "w"))==NULL)
		{
			fprintf(stderr, "\n  Error->Cannot open file\n");
		}
		fprintf(f,"# FPR REC PPV  1 vs 2 and 2 vs 3\n");

		float avgroc[BINSAUC][3];
		float sigroc[BINSAUC][3];

		for(int i=0; i<BINSAUC; i++)
		{

			// save all plots
//			for(int j=0; j<KNNRUNS; j++)
//			{
//				fprintf(f,"%7.3f %7.3f %7.3f ", Afpr[i][j], Arec[i][j], Appv[i][j]);
//			}

			// avg
			avgroc[i][0]=0; avgroc[i][1]=0; avgroc[i][2]=0;
			for(int j=0; j<KNNRUNS; j++)
			{
				avgroc[i][0]+=Afpr[i][j];
				avgroc[i][1]+=Arec[i][j];
				avgroc[i][2]+=Appv[i][j];
			}
			avgroc[i][0]/=(KNNRUNS*1.0);
			avgroc[i][1]/=(KNNRUNS*1.0);
			avgroc[i][2]/=(KNNRUNS*1.0);



			// sigma
			sigroc[i][0]=0; sigroc[i][1]=0; sigroc[i][2]=0;
			for(int j=0; j<KNNRUNS; j++)
			{
				sigroc[i][0]+=pow(Afpr[i][j]-avgroc[i][0],2);
				sigroc[i][1]+=pow(Arec[i][j]-avgroc[i][1],2);
				sigroc[i][2]+=pow(Appv[i][j]-avgroc[i][2],2);
			}
			sigroc[i][0]=sqrt(sigroc[i][0]/(KNNRUNS*1.0));
			sigroc[i][1]=sqrt(sigroc[i][1]/(KNNRUNS*1.0));
			sigroc[i][2]=sqrt(sigroc[i][2]/(KNNRUNS*1.0));

			fprintf(f,"%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f ",
					avgroc[i][0],
					avgroc[i][0]+1.96*sigroc[i][0]/sqrt(KNNRUNS*1.0),
					avgroc[i][0]-1.96*sigroc[i][0]/sqrt(KNNRUNS*1.0),
					avgroc[i][1],
					avgroc[i][1]+1.96*sigroc[i][1]/sqrt(KNNRUNS*1.0),
					avgroc[i][1]-1.96*sigroc[i][1]/sqrt(KNNRUNS*1.0),
					avgroc[i][2],
					avgroc[i][2]+1.96*sigroc[i][2]/sqrt(KNNRUNS*1.0),
					avgroc[i][2]-1.96*sigroc[i][2]/sqrt(KNNRUNS*1.0));





			fprintf(f,"\n");

		}
		fclose(f);

	}





	if (pclase) {

		if (verbose>0) fprintf(stdout," save classes.txt");

		if ( (f=fopen("classes.txt", "w"))==NULL)
		{
			fprintf(stderr, "\n  Error->Cannot open file\n");
		}

		for (int i = 0 ; i < nclasses; i++)
		{
			fprintf(f,"%5d %5.3f %5d\n", classindex[i], classcountF[classindex[i]]/classcountT[classindex[i]], statsH[ classindex[i] ][0] );

		}
		fclose(f);

	}


	free(Rmse);
	free(Mae);
	free(Ppv);
	free(Sen);
	free(Pcor);
	free(Acc);
	free(Err);
	free(Mcc);
	free(classcountF);
	free(classcountT);
	free(classcountN);
	free(classindex);
	for (int index=0;index<NVAR;++index)
		free(Weights[index]);
	free(Weights);


}
