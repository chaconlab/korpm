#include <stdio.h>


#include <Macromolecule.h>
#include <pdbIter.h>
#include <libtools/include/Memstuff.h>
// PDB atoms to prevent char warnings
char *Nstr   = (char *)" N  ";
char *CAstr  = (char *)" CA ";
char *Cstr   = (char *)" C  ";
char *Ostr   = (char *)" O  ";
char *CBstr  = (char *)" CB ";
char *CGstr  = (char *)" CG ";
char *CDstr  = (char *)" CD ";
char *CEstr  = (char *)" CE ";
char *CZstr  = (char *)" CZ ";
char *HAstr  = (char *)" HA ";
char *Hstr   = (char *)" H  ";

using namespace std;

//Auxiliar functions
#define f_a(A,B,C,D,E)    ((((A-2*B+C)/(D-2*A+B))+2.881525456)/(((E-2*D+A)/(D-2*A+B))-2.532088887))
#define f_pi(A,B,C,D,E)   ((((A-2*B+C)/(D-2*A+B))+1.0)/((E-2*D+A)/(D-2*A+B)))
#define f_310(A,B,C,D,E)  ((((A-2*B+C)/(D-2*A+B))-3.513337098)/(((E-2*D+A)/(D-2*A+B))+3.228707461))
#define f_b(A,B,C,D)      ((A-B)/(C-D))
void aux_min_w (double **d, int p, int i, int start, int end,  float *min1, float *min2);

// --------------------------------------------------------------------------------------
// minRmsd's Auxiliar functions and defines declaration (needed by: minRmsd and minWRmsd)
// --------------------------------------------------------------------------------------
#define PI 3.14159265358979323846
#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
		h = a[k][l]; \
		a[i][j] = g-s*(h+g*tau); \
		a[k][l] = h+s*(g-h*tau); }
/* vector functions using c arrays */
void normalize(double a[3]);
double dot(double a[3], double b[3]);
void cross(double a[3], double b[3], double c[3]);
/*
 * setup_rotation()
 *
 *      given two lists of x,y,z coordinates, constructs
 * the correlation R matrix and the E value needed to calculate the
 * least-squares rotation matrix.
 */
void setup_rotation(double **ref_xlist, double **mov_xlist, int n_list, double mov_com[3],
		double mov_to_ref[3], double R[3][3], double* E0,
		double *weights=NULL);
/*
 * jacobi3
 *
 *    computes eigenval and eigen_vec of a real 3x3
 * symmetric matrix. On output, elements of a that are above
 * the diagonal are destroyed. d[1..3] returns the
 * eigenval of a. v[1..3][1..3] is a matrix whose
 * columns contain, on output, the normalized eigen_vec of
 * a. n_rot returns the number of Jacobi rotations that were required.
 */
int jacobi3(double a[3][3], double d[3], double v[3][3], int* n_rot);
/*
 * diagonalize_symmetric
 *
 *    Diagonalize a 3x3 matrix & sort eigenval by size
 */
int diagonalize_symmetric(double matrix[3][3], double eigen_vec[3][3], double eigenval[3]);
/*
 * calculate_rotation_matrix()
 *
 *   calculates the rotation matrix U and the
 * rmsd from the R matrix and E0:
 */
int calculate_rotation_matrix(double R[3][3], double U[3][3], double E0, double* residual);
// --------------------------------------------------------
// END minRmsd's Auxiliar functions and defines declaration
// --------------------------------------------------------

void Macromolecule::info( FILE *file )
{
	bool contP, contCh, contF, contR, contA;
	int cont,cont_p,cont_ch,cont_seg,tot_mol;
	int tot_ch=0,tot_seg=0,tot_frag=0,tot_at=0;
	int csmol=0;
	Segment *seg;
	Chain *cad;
	Molecule *p;
	Fragment *frag;
	// Atom *a;
	char lineM[100];
	char lineMv[100];
	char lineC[100];
	char lineCv[100];
	char lineS[100];

	fprintf( file, "molinf>\n");
	fprintf( file, "molinf> Mol Chains   Segments    Atoms\n");
	fprintf( file, "molinf> --- ------ -----------  -------\n");

	initAll();

	tot_mol = this->getLimit();
	contP = true;
	cont_p = 0;
	while ( contP != false )
	{
		cont_p++;
		p = ( Molecule * ) this->getCurrent();

		sprintf( lineM, "molinf> %d",cont_p);

		if(p->getMolType()!=tmol_smol)
		{
			//    	if(cont_p == 1)
			//    	{
			if(p->getMolType()==tmol_protein)
				sprintf( lineM, "molinf>%3d",cont_p);
			if(p->getMolType()==tmol_nacid)
				sprintf( lineM, "molinf>%3d",cont_p);
			//    	}
			//    	else
			//   			sprintf( lineM,     "molinf>         %3d",cont_p);
			sprintf( lineMv,        "molinf>          ");

			cont_ch = 0;
			contCh = true;
			while ( contCh != false )
			{
				cont_ch++;
				tot_ch++;
				cad = ( Chain * ) p->getCurrent();

				if(cont_ch == 1)
					sprintf(lineC, "%s %3d->%c %d",lineM,cont_ch,cad->getName()[0], cad->getLimit());
				else
					sprintf(lineC, "%s %3d->%c %d",lineMv,cont_ch,cad->getName()[0], cad->getLimit());
				sprintf(lineCv, "%s",lineMv);
				cont_seg = 0;
				contF = true;
				while ( contF != false )
				{
					cont_seg++;
					tot_seg++;
					seg = ( Segment * ) cad->getCurrent();
					if(cont_seg == 1)
					{
						if(seg->getMolType()==tmol_protein)
							sprintf(lineS, "%s\n%s %3d %5d aa",lineC,lineCv,cont_seg,seg->getLimit());
						if(seg->getMolType()==tmol_nacid)
							sprintf(lineS, "%s\n%s %3d %5d na",lineC,lineCv,cont_seg,seg->getLimit());
					}
					else
					{
						if(seg->getMolType()==tmol_protein)
							sprintf(lineS, "%s %3d %5d aa",lineCv,cont_seg,seg->getLimit());
						if(seg->getMolType()==tmol_nacid)
							sprintf(lineS, "%s %3d %5d na",lineCv,cont_seg,seg->getLimit());
					}
					cont = 0;
					contR = true;
					while ( contR != false )
					{
						tot_frag++;
						frag = ( Fragment * ) seg->getCurrent();
						contA = true;
						while ( contA != false )
						{
							// a = ( Atom * ) frag->getCurrent();
							cont++;
							contA = frag->next();
						}
						contR = seg->next();
					}
					fprintf( file, "%s %6d\n", lineS,cont );
					tot_at += cont;
					contF = cad->next();
				}
				contCh = p->next();
			}
		}
		else
		{
			csmol++;
			// fprintf( file, "\tSmall Molecule:\n\tNumber of atoms: %d\n", p->getLimit() );
		}
		contP = this->next();
	}
	fprintf( file, "molinf>\n");
	fprintf( file, "molinf> --- ------ -----------  -------\n");
	fprintf( file, "molinf> SUMMARY:\n");
	fprintf( file, "molinf> Number of Molecules ..... %d\n",tot_mol );
	fprintf( file, "molinf> Number of Chain ......... %d\n",tot_ch);
	fprintf( file, "molinf> Number of Segments ...... %d\n",tot_seg);
	fprintf( file, "molinf> Number of Groups ........ %d\n",tot_frag);
	fprintf( file, "molinf> Number of Atoms ......... %d\n",tot_at);
	fprintf( file, "molinf> Number of Hetero Mol .... %d\n",csmol);
	fprintf( file, "molinf> Number of Hetero Atoms... %d\n",csmol);

	fprintf( file, "molinf>\n");
}




bool Macromolecule::geoCenter( Tcoor position )
{
	Atom * a;

	Tcoor pos;
	int cont = 0;

	position[0] = 0;
	position[1] = 0;
	position[2] = 0;

	this->initAll();

	do
	{

		a = this->getCurrentAtom();
		a->getPosition(pos);
		cont++;

		position[0] += pos[0];
		position[1] += pos[1];
		position[2] += pos[2];

	}
	while ( this->nextAtom() != false );





	if ( cont > 0 )
	{
		position[0] /= ( float )cont;
		position[1] /= ( float )cont;
		position[2] /= ( float )cont;
		return true;
	}
	else
	{
		fprintf(stdout,"Error_geoCenter: Division by zero\n");
		return false;
	}


}
/**
 * Computes the dimension of the system that are stored in the macromoleculeDimension structure
 * The center stored is computed as the middle point between the farthest atoms of the system (included hydrogens)
 */
void Macromolecule::geoBox()
{
	bool fin;
	Atom * a;
	Tcoor pos;
	initAll();

	macromoleculeDimension.Xmax = -1000000;
	macromoleculeDimension.Ymax = -1000000;
	macromoleculeDimension.Zmax = -1000000;

	macromoleculeDimension.Xmin = 1000000;
	macromoleculeDimension.Ymin = 1000000;
	macromoleculeDimension.Zmin = 1000000;
	fin = true;
	while ( fin != false )
	{
		a = getCurrentAtom();
		a->getPosition(pos);

		if ( pos[0] > macromoleculeDimension.Xmax )
			macromoleculeDimension.Xmax = pos[0];
		if ( pos[0] < macromoleculeDimension.Xmin )
			macromoleculeDimension.Xmin = pos[0];

		if ( pos[1] > macromoleculeDimension.Ymax )
			macromoleculeDimension.Ymax = pos[1];
		if ( pos[1] < macromoleculeDimension.Ymin )
			macromoleculeDimension.Ymin = pos[1];

		if ( pos[2] > macromoleculeDimension.Zmax )
			macromoleculeDimension.Zmax = pos[2];
		if ( pos[2] < macromoleculeDimension.Zmin )
			macromoleculeDimension.Zmin = pos[2];

		fin = nextAtom();
	}
	macromoleculeDimension.geometricCenter[0] =
			( macromoleculeDimension.Xmax - macromoleculeDimension.Xmin ) / 2 + macromoleculeDimension.Xmin;
	macromoleculeDimension.geometricCenter[1] =
			( macromoleculeDimension.Ymax - macromoleculeDimension.Ymin ) / 2 + macromoleculeDimension.Ymin;
	macromoleculeDimension.geometricCenter[2] =
			( macromoleculeDimension.Zmax - macromoleculeDimension.Zmin ) / 2 + macromoleculeDimension.Zmin;
}


float Macromolecule::maxLength()
{
	float maxL = 0, currentL;
	Tcoor pos;
	Atom * a;

	initAll();
	do
	{
		a = getCurrentAtom();
		a->getPosition(pos);
		currentL = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );
		if ( currentL > maxL )
			maxL = currentL;

	}
	while ( nextAtom() != false );
	return ( maxL * 2 );
}

float Macromolecule::maxLength_backbone()
{
	float maxL = 0, currentL;
	Tcoor pos;
	Atom * a;



	initAll();
	do
	{
		a = getCurrentAtom();
		if(   strcmp(a->getName()," C  ")==0
				||strcmp(a->getName()," CA ")==0
				||strcmp(a->getName()," N  ")==0
				||strcmp(a->getName()," O  ")==0
				//||strcmp(a->getName()," CB ")==0
		)
		{
			a->getPosition(pos);
			currentL = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );
			if ( currentL > maxL )
				maxL = currentL;
		}
	}
	while ( nextAtom() != false );
	return ( maxL * 2 );
}

// Mon modified (26/3/2010) and (8/4/2010)
// (Apply each table entry to the corresponding atom b-factor)
void Macromolecule::exchange_Pdbfact(double *table, bool fragment)
{
	int i;
	Segment * seg;
	Atom * atom;
	pdbIter *iter1, *iter2;
	Fragment * res;
	//char * at_name;


	//Iterador para recorrer segmentos
	iter1 = new pdbIter( this );

	//Bucle para recorrer segmentos
	i=0;
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );

		//Bucle para recorrer residuos
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			res = ( Fragment * ) iter2->get_fragment();

			//i++;
			//Iterador sobre los atomos del residuo
			pdbIter * iter3 = new pdbIter( res );

			for ( iter3->pos_atom = 0; !iter3->gend_atom(); iter3->next_atom() )
			{
				atom = ( Atom * ) iter3->get_atom();
				// at_name = atom->getName();
				atom->setPdbfact(table[i]);
				if(!fragment) // Mon added (8/4/2010)
					i++; // Mon added (26/3/2010)
			}
			iter3->~pdbIter();
			// Mon modified (26/3/2010)
			//      i++;
			if(fragment) // Mon added (8/4/2010)
				i++;
		}
		iter2->~pdbIter();
	}
	iter1->~pdbIter();





}


void Macromolecule::get_Pdbfact(double *table)
{
	int i;
	Atom * a;

	initAll();
	i=0;
	do
	{
		a = getCurrentAtom();
		table[i]=a->getPdbfact();
		i++;
	}
	while ( nextAtom() != false );

}



float Macromolecule::maxDist()
{
	pdbIter *iter1,*iter2;
	Tcoor pos1,pos2;

	iter1=new pdbIter(this);
	iter2=new pdbIter(this);

	float maxL = 0, currentL;
	float sum1,sum2,sum3;


	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		(iter1->get_atom())->getPosition(pos1);

		for(iter2->pos_atom=iter1->pos_atom+1;!iter2->gend_atom();iter2->next_atom())
		{
			(iter2->get_atom())->getPosition(pos2);
			sum1=pos1[0] - pos2[0];
			sum2=pos1[1] - pos2[1];
			sum3=pos1[2] - pos2[2];

			currentL = (float)sqrt( sum1*sum1 + sum2*sum2 + sum3*sum3);
			if ( currentL > maxL )
				maxL = currentL;
		}
	}

	delete iter1;
	delete iter2;

	return maxL;
}



bool Macromolecule::geoProperty( float * position, char propierty )
{
	Atom * a;
	Element * e;
	Tcoor pos;
	float value, total = 0;
	position[0] = 0;
	position[1] = 0;
	position[2] = 0;


	initAll();
	do
	{
		a = getCurrentAtom();
		a->getPosition(pos);
		e = a->getElement();

		switch ( propierty )
		{
		case( 'W' ):
        								value = e->weight;
		break;
		case( 'A' ):
        								value = e->atomic;
		break;
		default:
			fprintf(stdout,"Error_geoPropierty: Incorrect propierty identifier\n");
			return false;
			break;
		}

		total += value;
		position[0] += pos[0] * value;
		position[1] += pos[1] * value;
		position[2] += pos[2] * value;

	}
	while ( nextAtom() != false );

	if ( total > 0 )
	{
		position[0] /= total;
		position[1] /= total;
		position[2] /= total;
		return true;
	}
	else
	{
		fprintf(stdout,"Error_geoPropierty: Division by zero\n");
		return false;
	}
}

float Macromolecule::eDensity()
{
	bool end = true;
	Atom * at;
	float suma = 0;

	initAll();
	while ( end != false )
	{
		/* Obtengo un atomo */
		at = getCurrentAtom();
		/* Obtencion de la posicion real del atomo en el PDb en amstrongs */
		suma += ( at->getElement() )->number;
		end = nextAtom();
	}
	return ( suma );
}



Macromolecule * Macromolecule::select( Conditions * cond, bool positive, bool include_hetero)
{

	bool contM, contCh, contS, contR, contA;
	bool insert;
	int cM,cCh,cS, cR, cA ;
	cM=cCh=cS= cR= cA=-1;

	Segment * seg;
	Chain * cad;
	Molecule * mol;
	Fragment * res;
	Atom * a;

	Segment * seg2;
	Chain * cad2;
	Molecule * mol2;
	Fragment * res2;


	Macromolecule * m = new Macromolecule();

	initAll();
	contM = true;
	while ( contM != false )
	{
		mol = ( Molecule * ) this->getCurrent();

		cM++;
		switch(mol->getClass())
		{
		case pdb_protein:
			mol2 = new Protein();
			break;
		case pdb_nacid:
			mol2 = new NAcid();
			break;
		case pdb_smol:
			mol2 = new SMol( mol->getName(),mol->getIdNumber() );
			break;
		default:
			break;
		}

		if(mol->getClass()==pdb_smol)
		{
			contA = true;
			while ( contA != false )
			{
				cA++;
				a = ( Atom * ) mol->getCurrent();

				if(!include_hetero) // Chapa to include all Hetero-atoms in selection...
				{
					insert = false;

					for ( int i = 0; i < cond->limit && !insert; i++ )
					{
						if ( ( cond->getCondition( i ) )->fitCondition( cM, -1, -1, -1, cA, a->getPdbName() ) == true )
							insert = true;
					}

					if(!positive && insert)
						insert=false;
					else
						if(!positive && !insert)
							insert=true;

					if ( insert == true )
						mol2->add( a );
				}
				else
					mol2->add( a ); // always adds the atom

				contA = mol->next();
			}
			if(mol2->getLimit()>0)
				m->add(mol2);
			else
				delete mol2;

		}
		else
		{
			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) mol->getCurrent();

				cCh++;
				cad2 = new Chain( cad->getName() );

				contS = true;
				while ( contS != false )
				{
					seg = ( Segment * ) cad->getCurrent();

					cS++;
					seg2 = new Segment();

					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();

						//cR++;
						cR=res->getIdNumber();
						if(res->getClass()==pdb_residue)
							res2 = new Residue( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());
						if(res->getClass()==pdb_nucleotide)
							res2 = new Nucleotide( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());

						contA = true;
						while ( contA != false )
						{
							cA++;
							a = ( Atom * ) res->getCurrent();
							insert = false;

							for ( int i = 0; i < cond->limit && !insert; i++ )
							{
								if ( ( cond->getCondition( i ) )->fitCondition( cM, cS, cCh, cR, cA, a->getPdbName() ) == true )
									insert = true;
							}

							if(!positive && insert)
								insert=false;
							else
								if(!positive && !insert)
									insert=true;


							if ( insert == true )
							{
								// cout<<"Insertando atomo: "<<a->getPdbName()<<endl;
								res2->add( a );
							}

							contA = res->next();
						}

						if ( res2->getLimit() > 0 )
							seg2->add( res2 );
						else
							delete res2;

						contR = seg->next();
					}

					if ( seg2->getLimit() > 0 )
						cad2->add( seg2 );
					else
						delete seg2;

					contS = cad->next();
				}

				if ( cad2->getLimit() > 0 )
					mol2->add( cad2 );
				else
					delete cad2;

				contCh = mol->next();
			}

			if ( mol2->getLimit() > 0 )
				m->add( mol2 );
			else
				delete mol2;
		}

		contM = this->next();

	}
	return m;

}

Macromolecule * Macromolecule::select_cpy( Conditions * cond, bool positive )
{

	bool contM, contCh, contS, contR, contA;
	bool insert;
	int cM,cCh,cS, cR, cA ;
	cM=cCh=cS= cR= cA=-1;

	Segment * seg;
	Chain * cad;
	Molecule * mol;
	Fragment * res;
	Atom * a;

	Segment * seg2;
	Chain * cad2;
	Molecule * mol2;
	Fragment * res2;


	Macromolecule * m = new Macromolecule();

	initAll();
	contM = true;
	while ( contM != false )
	{
		mol = ( Molecule * ) this->getCurrent();
		cM++;

		switch(mol->getClass())
		{
		case pdb_protein:
			mol2 = new Protein();
			break;
		case pdb_nacid:
			mol2 = new NAcid();
			break;
		case pdb_smol:
			mol2 = new SMol( mol->getName(),mol->getIdNumber() );
			break;
		default:
			break;
		}

		if(mol->getClass()==pdb_smol)
		{

			contA = true;
			while ( contA != false )
			{
				cA++;
				a = ( Atom * ) mol->getCurrent();
				insert = false;

				for ( int i = 0; i < cond->limit && !insert; i++ )
				{
					if ( ( cond->getCondition( i ) )->fitCondition( cM, -1, -1, -1, cA, a->getPdbName() ) == true )
						insert = true;
				}

				if(!positive && insert)
					insert=false;
				else
					if(!positive && !insert)
						insert=true;


				if ( insert == true )
				{
					mol2->add( new Atom(a) );
				}

				contA = mol->next();
			}


		}
		else
		{
			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) mol->getCurrent();

				cCh++;
				cad2 = new Chain( cad->getName() );

				contS = true;
				while ( contS != false )
				{
					seg = ( Segment * ) cad->getCurrent();

					cS++;
					seg2 = new Segment();

					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();

						//cR++;
						cR=res->getIdNumber();
						if(res->getClass()==pdb_residue)
							res2 = new Residue( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());
						if(res->getClass()==pdb_nucleotide)
							res2 = new Nucleotide( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());

						contA = true;
						while ( contA != false )
						{
							cA++;
							a = ( Atom * ) res->getCurrent();
							insert = false;

							for ( int i = 0; i < cond->limit && !insert; i++ )
							{
								if ( ( cond->getCondition( i ) )->fitCondition( cM, cS, cCh, cR, cA, a->getPdbName() ) == true )
									insert = true;
							}

							if(!positive && insert)
								insert=false;
							else
								if(!positive && !insert)
									insert=true;


							if ( insert == true )
							{
								// cout<<"Insertando atomo: "<<a->getPdbName()<<endl;
								res2->add( new Atom(a) );
							}


							contA = res->next();
						}

						if ( res2->getLimit() > 0 )
							seg2->add( res2 );
						else
							delete res2;

						contR = seg->next();
					}

					if ( seg2->getLimit() > 0 )
						cad2->add( seg2 );
					else
						delete seg2;

					contS = cad->next();
				}

				if ( cad2->getLimit() > 0 )
					mol2->add( cad2 );
				else
					delete cad2;

				contCh = mol->next();
			}
		}
		if ( mol2->getLimit() > 0 )
			m->add( mol2 );
		else
			delete mol2;
		contM = this->next();
	}
	return m;

}


// Mon created (27/04/2007)
// Selects atoms if name matches.
Macromolecule * Macromolecule::select_aaName( char **aaName )
{
	bool contM, contCh, contS, contR, contA;
	bool insert;


	Segment * seg;
	Chain * cad;
	Molecule * mol;
	Fragment * res;
	Atom * a;

	Segment * seg2;
	Chain * cad2;
	Molecule * mol2;
	Fragment * res2;


	Macromolecule * m = new Macromolecule();

	initAll();
	contM = true;
	while ( contM != false )
	{
		mol = ( Molecule * ) this->getCurrent();

		switch(mol->getClass())
		{
		case pdb_protein:
			mol2 = new Protein();
			break;
		case pdb_nacid:
			mol2 = new NAcid();
			break;
		case pdb_smol:
			mol2 = new SMol( mol->getName(),mol->getIdNumber() );
			break;
		default:
			break;
		}

		if(mol->getClass()==pdb_smol)
		{

			contA = true;
			while ( contA != false )
			{
				a = ( Atom * ) mol->getCurrent();
				insert = false;

				for( int j = 0; (aaName[j] != NULL) && !insert; j++ ) // screens whether some name matches
					if( strncmp(mol->getName(),aaName[j],3) == 0 )
					{
						insert = true;
						mol2->add( a );
					}

				contA = mol->next();
			}

		}
		else
		{

			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) mol->getCurrent();
				cad2 = new Chain( cad->getName() );

				contS = true;
				while ( contS != false )
				{
					seg = ( Segment * ) cad->getCurrent();

					seg2 = new Segment();

					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();

						if(res->getClass()==pdb_residue)
							res2 = new Residue( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());
						if(res->getClass()==pdb_nucleotide)
							res2 = new Nucleotide( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());

						contA = true;
						while ( contA != false )
						{
							a = ( Atom * ) res->getCurrent();
							insert = false; // to avoid unnecessary comparations
							for( int j = 0; (aaName[j] != NULL) && !insert; j++ ) // screens whether some name matches
								if( strncmp(res->getName(),aaName[j],3) == 0 )
								{
									insert = true;
									res2->add( a );
								}

							contA = res->next();
						}

						if ( res2->getLimit() > 0 )
							seg2->add( res2 );

						contR = seg->next();
					}

					if ( seg2->getLimit() > 0 )
						cad2->add( seg2 );

					contS = cad->next();
				}

				if ( cad2->getLimit() > 0 )
					mol2->add( cad2 );

				contCh = mol->next();
			}
		}

		if ( mol2->getLimit() > 0 )
			m->add( mol2 );

		contM = this->next();
	}
	return m;

}

// me parece una patata...esto hay que cambiarlo
int Macromolecule::get_num_atoms()
{
	/* int value;
  pdbIter *iter1;
  iter1=new pdbIter(this);
  value=(int) iter1->num_atom();
  delete iter1;
  return value;*/

	return num_atoms();

}

int Macromolecule::get_num_fragments()
{
	int value;
	pdbIter *iter1;
	iter1=new pdbIter(this);
	value=(int) iter1->num_fragment();
	delete iter1;
	return value;
}

// New "minRmsd" implementation (Mon 1/6/2010) from:
// http://www.koders.com/cpp/fid5705387AFD679C2E0A9F25C24204E2B0C3256EB8.aspx?s=%22Marcelo+Silveira%22
// (It fixes some problems related to highly symmetric structures)
float Macromolecule::minRmsd( Macromolecule * mol2, float mdest[4][4] )
{
	double **ref_xlist;
	double **mov_xlist;
	bool end;
	bool end2;
	Tcoor pos1,pos2;
	Atom *at,*at2;
	int num_atoms=0;
	int n,i,j;
	num_atoms = get_num_atoms();

	end = true;
	end2 = true;

	//	fprintf(stderr,"USING new minRmsd()\n");
	ref_xlist = (double **) malloc(sizeof(double *) * num_atoms );
	mov_xlist = (double **) malloc(sizeof(double *) * num_atoms );
	for(n = 0; n < num_atoms; n++ )
	{
		ref_xlist[n] = (double *) malloc(sizeof(double)*3 );
		mov_xlist[n] = (double *) malloc(sizeof(double)*3 );
	}

	initAll();
	mol2->initAll();
	int natoms = 0;
	while ( end != false )
	{
		/* get atoms */
		at = getCurrentAtom();
		at2 = mol2->getCurrentAtom();

		/* get positions */
		at->getPosition(pos1);
		at2->getPosition(pos2);

		for (int i = 0; i < 3; i++ )
		{
			ref_xlist[natoms][i] = pos1[i];
			mov_xlist[natoms][i] = pos2[i];
		}

		end = nextAtom();
		end2 = (bool)mol2->nextAtom();
		natoms++;
	}

	int n_list = num_atoms;
	double mov_com[3];
	double mov_to_ref[3];
	double U[3][3];
	double Eo, residual;
	double R[3][3];

	setup_rotation(ref_xlist, mov_xlist, n_list, mov_com, mov_to_ref, R, &Eo);
	calculate_rotation_matrix(R, U, Eo, &residual);

	residual = fabs(residual); /* avoids the awkward case of -0.0 */
	double rmsd = sqrt( fabs((double) (residual)*2.0/((double)n_list)) );

	// Output rotation matrix
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			mdest[i][j] = (float) U[i][j];

	// get ref_com[3] vector (I don't know why it's not output by setup_rotation)
	mov_to_ref[0] += mov_com[0];
	mov_to_ref[1] += mov_com[1];
	mov_to_ref[2] += mov_com[2];

	// Translational
	mdest[0][3] = mov_to_ref[0]  - (mdest[0][0]*mov_com[0] + mdest[0][1]*mov_com[1] + mdest[0][2]*mov_com[2]);
	mdest[1][3] = mov_to_ref[1]  - (mdest[1][0]*mov_com[0] + mdest[1][1]*mov_com[1] + mdest[1][2]*mov_com[2]);
	mdest[2][3] = mov_to_ref[2]  - (mdest[2][0]*mov_com[0] + mdest[2][1]*mov_com[1] + mdest[2][2]*mov_com[2]);
	mdest[3][3] = 1.0;

	// Free memory
	for(n = 0; n < num_atoms; n++ )
	{
		free( ref_xlist[n] );
		free( mov_xlist[n] );
	}
	free( ref_xlist );
	free( mov_xlist );

	return rmsd;
}

// Overloading input to take into account alignment masks at atomic level (17/7/2012)
// New "minRmsd" implementation (Mon 1/6/2010) from:
// http://www.koders.com/cpp/fid5705387AFD679C2E0A9F25C24204E2B0C3256EB8.aspx?s=%22Marcelo+Silveira%22
// (It fixes some problems related to highly symmetric structures)
float Macromolecule::minRmsd( Macromolecule * mol2, float mdest[4][4], bool *mask1, bool *mask2 )
{
	bool debug = false;
	double **ref_xlist;
	double **mov_xlist;
	Tcoor pos1,pos2;
	Atom *at1,*at2;
	int num_atoms=0;
	int n,i,j,ncommon=0;
	num_atoms = get_num_atoms();

	ref_xlist = (double **) malloc(sizeof(double *) * num_atoms );
	mov_xlist = (double **) malloc(sizeof(double *) * num_atoms );

	pdbIter *iter1,*iter2;
	iter1 = new pdbIter(this,false);
	iter2 = new pdbIter(mol2,false);

	iter2->pos_atom=0;
	for(iter1->pos_atom=0; !iter1->gend_atom(); iter1->next_atom())
	{
		if(mask1[iter1->pos_atom]) // Mol 1 matching atom
		{
			while(!iter2->gend_atom() && !mask2[iter2->pos_atom]) // this places the index into the corresponding atom
				iter2->next_atom();

			if(!iter2->gend_atom())
			{
				// get atoms
				at1 = iter1->get_atom();
				at2 = iter2->get_atom();

				// get positions
				at1->getPosition(pos1);
				at2->getPosition(pos2);

				if(debug)
					printf("Getting atoms %d (%s) and %d (%s)\n",iter1->pos_atom,at1->getName(),iter2->pos_atom,at2->getName());

				// allocating memory
				ref_xlist[ncommon] = (double *) malloc(sizeof(double)*3 );
				mov_xlist[ncommon] = (double *) malloc(sizeof(double)*3 );

				for (int i = 0; i < 3; i++ )
				{
					ref_xlist[ncommon][i] = pos1[i];
					mov_xlist[ncommon][i] = pos2[i];
				}

				ncommon++;

				iter2->next_atom();
			}
		}
	}

	delete(iter1);
	delete(iter2);

	int n_list = ncommon;
	double mov_com[3];
	double mov_to_ref[3];
	double U[3][3];
	double Eo, residual;
	double R[3][3];

	setup_rotation(ref_xlist, mov_xlist, n_list, mov_com, mov_to_ref, R, &Eo);
	calculate_rotation_matrix(R, U, Eo, &residual);

	residual = fabs(residual); /* avoids the awkward case of -0.0 */
	double rmsd = sqrt( fabs((double) (residual)*2.0/((double)n_list)) );

	// Output rotation matrix
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			mdest[i][j] = (float) U[i][j];

	// get ref_com[3] vector (I don't know why it's not output by setup_rotation)
	mov_to_ref[0] += mov_com[0];
	mov_to_ref[1] += mov_com[1];
	mov_to_ref[2] += mov_com[2];

	// Translational
	mdest[0][3] = mov_to_ref[0]  - (mdest[0][0]*mov_com[0] + mdest[0][1]*mov_com[1] + mdest[0][2]*mov_com[2]);
	mdest[1][3] = mov_to_ref[1]  - (mdest[1][0]*mov_com[0] + mdest[1][1]*mov_com[1] + mdest[1][2]*mov_com[2]);
	mdest[2][3] = mov_to_ref[2]  - (mdest[2][0]*mov_com[0] + mdest[2][1]*mov_com[1] + mdest[2][2]*mov_com[2]);
	mdest[3][3] = 1.0;

	// Free memory
	for(n = 0; n < ncommon; n++ )
	{
		free( ref_xlist[n] );
		free( mov_xlist[n] );
	}
	free( ref_xlist );
	free( mov_xlist );

	return rmsd;
}


// Mon made (1/6/2010)
// Weighted RMSD minimization routine
// "profile" --> residue level weights (one per residue)
// See: Damm & Carlson (2005)
// http://www.koders.com/cpp/fid5705387AFD679C2E0A9F25C24204E2B0C3256EB8.aspx?s=%22Marcelo+Silveira%22
// (It fixes some problems related to highly symmetric structures)
float Macromolecule::minWRmsd( Macromolecule * mol2, float mdest[4][4], double *profile )
{
	double **ref_xlist;
	double **mov_xlist;
	bool end = true;
	bool end2 = true;
	Tcoor pos1,pos2;
	Atom *at,*at2;
	int num_atoms=0;
	int n,i,j;
	num_atoms = get_num_atoms();

	//	fprintf(stderr,"USING new minRmsd()\n");
	ref_xlist = (double **) malloc(sizeof(double *) * num_atoms );
	mov_xlist = (double **) malloc(sizeof(double *) * num_atoms );
	for(n = 0; n < num_atoms; n++ )
	{
		ref_xlist[n] = (double *) malloc(sizeof(double)*3 );
		mov_xlist[n] = (double *) malloc(sizeof(double)*3 );
	}

	initAll();
	mol2->initAll();
	int natoms = 0;
	while ( end != false )
	{
		/* get atoms */
		at = getCurrentAtom();
		at2 = mol2->getCurrentAtom();

		/* get positions */
		at->getPosition(pos1);
		at2->getPosition(pos2);

		for (i = 0; i < 3; i++ )
		{
			ref_xlist[natoms][i] = pos1[i];
			mov_xlist[natoms][i] = pos2[i];
		}

		end = nextAtom();
		end2 = mol2->nextAtom();
		natoms++;
	}

	int n_list = num_atoms;
	double mov_com[3];
	double mov_to_ref[3];
	double U[3][3];
	double Eo, residual;
	double R[3][3];

	setup_rotation(ref_xlist, mov_xlist, n_list, mov_com, mov_to_ref, R, &Eo, profile);
	calculate_rotation_matrix(R, U, Eo, &residual);

	residual = fabs(residual); /* avoids the awkward case of -0.0 */
	double rmsd = sqrt( fabs((double) (residual)*2.0/((double)n_list)) );

	// Output rotation matrix
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			mdest[i][j] = (float) U[i][j];

	// get ref_com[3] vector (I don't know why it's not output by setup_rotation)
	mov_to_ref[0] += mov_com[0];
	mov_to_ref[1] += mov_com[1];
	mov_to_ref[2] += mov_com[2];

	// Translational
	mdest[0][3] = mov_to_ref[0]  - (mdest[0][0]*mov_com[0] + mdest[0][1]*mov_com[1] + mdest[0][2]*mov_com[2]);
	mdest[1][3] = mov_to_ref[1]  - (mdest[1][0]*mov_com[0] + mdest[1][1]*mov_com[1] + mdest[1][2]*mov_com[2]);
	mdest[2][3] = mov_to_ref[2]  - (mdest[2][0]*mov_com[0] + mdest[2][1]*mov_com[1] + mdest[2][2]*mov_com[2]);

	// Free memory
	for(n = 0; n < num_atoms; n++ )
	{
		free( ref_xlist[n] );
		free( mov_xlist[n] );
	}
	free( ref_xlist );
	free( mov_xlist );

	return rmsd;
}


// Overloading input to take into account alignment masks at atomic level (17/7/2012)
// Weighted RMSD minimization routine (Mon made: 1/6/2010)
// "profile" --> residue level weights (one per residue)
// See: Damm & Carlson (2005)
// http://www.koders.com/cpp/fid5705387AFD679C2E0A9F25C24204E2B0C3256EB8.aspx?s=%22Marcelo+Silveira%22
// (It fixes some problems related to highly symmetric structures)
float Macromolecule::minWRmsd( Macromolecule * mol2, float mdest[4][4], double *profile, bool *mask1, bool *mask2 )
{
	//	bool debug = true;
	double **ref_xlist;
	double **mov_xlist;
	double *profile2;
	Tcoor pos1,pos2;
	Atom *at1,*at2;
	int num_atoms=0, ncommon=0;
	int n,i,j;
	num_atoms = get_num_atoms();

	//	if(debug)
	//		printf("num_atoms= %d\n",num_atoms);

	//	fprintf(stderr,"USING new minRmsd()\n");
	ref_xlist = (double **) malloc(sizeof(double *) * num_atoms );
	mov_xlist = (double **) malloc(sizeof(double *) * num_atoms );
	profile2 = (double *) malloc(sizeof(double) * num_atoms );

	if(!ref_xlist || !mov_xlist || !profile2)
	{
		printf("Msg(minWRmsd): Memory allocation failed! Forcing exit!\n");
		exit(1);
	}

	pdbIter *iter1,*iter2;
	iter1 = new pdbIter(this,false);
	iter2 = new pdbIter(mol2,false);

	iter2->pos_atom=0;
	for(iter1->pos_atom=0; !iter1->gend_atom(); iter1->next_atom())
	{
		if(mask1[iter1->pos_atom]) // Mol 1 matching atom
		{
			while(!iter2->gend_atom() && !mask2[iter2->pos_atom]) // this places the index into the corresponding atom
				iter2->next_atom();

			if(!iter2->gend_atom())
			{
				// get atoms
				at1 = iter1->get_atom();
				at2 = iter2->get_atom();

				// get positions
				at1->getPosition(pos1);
				at2->getPosition(pos2);

				// allocating memory
				ref_xlist[ncommon] = (double *) malloc(sizeof(double)*3 );
				mov_xlist[ncommon] = (double *) malloc(sizeof(double)*3 );

				for (i = 0; i < 3; i++ )
				{
					ref_xlist[ncommon][i] = pos1[i];
					mov_xlist[ncommon][i] = pos2[i];
				}
				profile2[ncommon] = profile[iter1->pos_atom]; // loading weights into a "compacted array"

				//				if(debug)
				//				{
				//					printf("Getting atoms %d (%s) and %d (%s)\n",iter1->pos_atom,at1->getName(),iter2->pos_atom,at2->getName());
				//					printf("\tref_xlist= %f %f %f    mov_xlist= %f %f %f\n",ref_xlist[ncommon][0],ref_xlist[ncommon][1],ref_xlist[ncommon][2]
				//					                                                       ,mov_xlist[ncommon][0],mov_xlist[ncommon][1],mov_xlist[ncommon][2]);
				//				}

				ncommon++;

				iter2->next_atom();
			}
		}
	}

	delete(iter1);
	delete(iter2);

	int n_list = ncommon;
	double mov_com[3];
	double mov_to_ref[3];
	double U[3][3];
	double Eo, residual;
	double R[3][3];

	setup_rotation(ref_xlist, mov_xlist, n_list, mov_com, mov_to_ref, R, &Eo, profile2);
	free(profile2);
	calculate_rotation_matrix(R, U, Eo, &residual);

	//	if(debug)
	//		printf("AFTER: n_list= %d  mov_com= %f %f %f  mov_to_ref= %f %f %f\n",n_list,mov_com[0],mov_com[1],mov_com[2],mov_to_ref[0],mov_to_ref[1],mov_to_ref[2]);

	residual = fabs(residual); /* avoids the awkward case of -0.0 */
	double rmsd = sqrt( fabs((double) (residual)*2.0/((double)n_list)) );

	// Output rotation matrix
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			mdest[i][j] = (float) U[i][j];

	// get ref_com[3] vector (I don't know why it's not output by setup_rotation)
	mov_to_ref[0] += mov_com[0];
	mov_to_ref[1] += mov_com[1];
	mov_to_ref[2] += mov_com[2];

	// Translational
	mdest[0][3] = mov_to_ref[0]  - (mdest[0][0]*mov_com[0] + mdest[0][1]*mov_com[1] + mdest[0][2]*mov_com[2]);
	mdest[1][3] = mov_to_ref[1]  - (mdest[1][0]*mov_com[0] + mdest[1][1]*mov_com[1] + mdest[1][2]*mov_com[2]);
	mdest[2][3] = mov_to_ref[2]  - (mdest[2][0]*mov_com[0] + mdest[2][1]*mov_com[1] + mdest[2][2]*mov_com[2]);

	//	if(debug)
	//		printf("%f %f %f  %f %f %f    %f %f %f\n", mdest[0][0], mdest[0][1], mdest[0][2], mdest[1][0], mdest[1][1], mdest[1][2], mdest[2][0], mdest[2][1], mdest[2][2]);

	// Free memory
	for(n = 0; n < ncommon; n++ )
	{
		free( ref_xlist[n] );
		free( mov_xlist[n] );
	}
	free( ref_xlist );
	free( mov_xlist );

	return rmsd;
}


float Macromolecule::rmsd( Macromolecule * mol2 )
{
	bool end = true, end2 = true;
	Atom * at, * at2;
	Tcoor pos1, pos2;
	double suma = 0, tmp;
	int a, natoms;

	initAll();
	mol2->initAll();
	natoms = 0;
	while ( end != false && end2!=false )
	{

		/* get atoms */
		at = getCurrentAtom();
		at2 = mol2->getCurrentAtom();

		/* get positions */
		at->getPosition(pos1);
		at2->getPosition(pos2);

		tmp=0;
		for ( a = 0; a < 3; a++ )
		{
			tmp += (pos1[a] - pos2[a])*(pos1[a] - pos2[a]);
		}
		suma += tmp;


		end = nextAtom(); end2 = mol2->nextAtom();
		natoms++;
	}

	//printf("natoms %d %f\n", natoms, suma);

	if ( suma != 0 ) return sqrt( ( float )suma / ( float )natoms );
	else
		return ( ( float )suma );
}

float Macromolecule::rmsd( Macromolecule * mol2, int **list, int total )
{
	bool end = true, end2 = true;
	Fragment *frag1,*frag2;
	Atom * at, * at2;
	Tcoor pos1, pos2;
	double suma = 0, tmp;
	int i,a, natoms;
	pdbIter *iter1,*iter2;


	iter1=new pdbIter(this,false);
	iter2=new pdbIter(mol2,false);

	natoms=0;
	for(i=0;i<total;i++)
	{
		if(list[0][i]>=0 && list[1][i]>=0)
		{
			iter1->pos_fragment=list[0][i];
			iter2->pos_fragment=list[1][i];
			frag1=iter1->get_fragment();
			frag2=iter2->get_fragment();
			frag1->initAll();
			frag2->initAll();

			//fprintf(stderr,"comparando %d %s con %d %s\n",list[0][i],frag1->getName(),list[1][i],frag2->getName());
			end=true;
			while(end != false)
			{
				/* get atoms */
				at = (Atom*)(frag1->getCurrent());
				at2 = (Atom*)(frag2->getCurrent());

				//fprintf(stderr,"Atomos> comparando  %s con  %s\n",at->getName(),at2->getName());

				/* get positions */
				at->getPosition(pos1);
				at2->getPosition(pos2);

				tmp=0;
				for ( a = 0; a < 3; a++ )
				{
					tmp += (pos1[a] - pos2[a])*(pos1[a] - pos2[a]);
				}
				suma += tmp;

				end = frag1->next(); end2 = frag2->next();
				natoms++;

			}
		}
	}

	delete(iter1);
	delete(iter2);

	if ( suma != 0 ) return sqrt(( float )suma / ( float )natoms );
	else
		return ( ( float )suma );

}

float Macromolecule::rmsd( Macromolecule * mol2, int **list, int total, Macromolecule *other1, Macromolecule *other2,int **list2, int total2)
{
	bool end = true, end2 = true;
	Fragment *frag1,*frag2;
	Atom * at, * at2;
	Tcoor pos1, pos2;
	double suma = 0,suma2, tmp;
	int cont2;
	int i,a, natoms;
	pdbIter *iter1,*iter2;


	natoms=0;
	iter1=new pdbIter(this,false);
	iter2=new pdbIter(mol2,false);

	for(i=0;i<total;i++)
	{
		if(list[0][i]>=0 && list[1][i]>=0)
		{
			iter1->pos_fragment=list[0][i];
			iter2->pos_fragment=list[1][i];
			frag1=iter1->get_fragment();
			frag2=iter2->get_fragment();
			frag1->initAll();
			frag2->initAll();

			//fprintf(stderr,"comparando %d %s con %d %s\n",list[0][i],frag1->getName(),list[1][i],frag2->getName());
			end=true;
			suma2=0;
			cont2=0;
			while(end != false)
			{
				/* get atoms */
				at = (Atom*)(frag1->getCurrent());
				at2 = (Atom*)(frag2->getCurrent());

				//fprintf(stderr,"Atomos> comparando  %s con  %s\n",at->getName(),at2->getName());

				/* get positions */
				at->getPosition(pos1);
				at2->getPosition(pos2);

				tmp=0;
				for ( a = 0; a < 3; a++ )
				{
					tmp += (pos1[a] - pos2[a])*(pos1[a] - pos2[a]);
				}
				suma += tmp;
				suma2 += tmp;
				cont2 ++;
				end = frag1->next(); end2 = frag2->next();
				natoms++;

			}
			//fprintf(stdout,"lig dist=%f pos1=%d pos2=%d natoms=%d\n",sqrt(( float )suma2 / ( float )cont2 ),list[0][i],list[1][i],cont2);

		}
	}

	delete(iter1);
	delete(iter2);

	iter1=new pdbIter(other1,false);
	iter2=new pdbIter(other2,false);
	//fprintf(stderr,"Other\n");
	for(i=0;i<total2;i++)
	{
		if(list2[0][i]>=0 && list2[1][i]>=0)
		{
			iter1->pos_fragment=list2[0][i];
			iter2->pos_fragment=list2[1][i];
			frag1=iter1->get_fragment();
			frag2=iter2->get_fragment();
			frag1->initAll();
			frag2->initAll();

			//fprintf(stderr,"comparando %d %s con %d %s\n",list2[0][i],frag1->getName(),list2[1][i],frag2->getName());
			end=true;
			suma2=0;
			cont2=0;
			while(end != false)
			{
				// get atoms
				at = (Atom*)(frag1->getCurrent());
				at2 = (Atom*)(frag2->getCurrent());

				//fprintf(stderr,"Atomos> comparando  %s con  %s\n",at->getName(),at2->getName());

				// get positions
				at->getPosition(pos1);
				at2->getPosition(pos2);

				tmp=0;
				for ( a = 0; a < 3; a++ )
				{
					tmp += (pos1[a] - pos2[a])*(pos1[a] - pos2[a]);
				}
				suma += tmp;
				suma2 += tmp;
				cont2 ++;
				end = frag1->next(); end2 = frag2->next();
				natoms++;
			}
			//fprintf(stdout,"res dist=%f pos1=%d pos2=%d natom=%d\n",sqrt(( float )suma2 / ( float )cont2 ),list2[0][i],list2[1][i],cont2);

		}
	}

	delete(iter1);
	delete(iter2);

	if ( suma != 0 ) return sqrt(( float )suma / ( float )natoms );
	else
		return ( ( float )suma );

}

float Macromolecule::rmsd( Macromolecule * mol2, char **list )
{
	bool end = false, end_atom,end_atom2;
	Fragment *frag1,*frag2;
	Atom * at, * at2;
	Tcoor pos1, pos2;
	double suma = 0, tmp;
	int i,a, natoms=0;
	pdbIter *iter1,*iter2;


	int cont=0;

	iter1=new pdbIter(this,false);
	iter2=new pdbIter(mol2,false);

	iter1->pos_fragment=0;
	iter2->pos_fragment=0;

	i=0;
	do
	{
		while((list[0][i]=='0' || list[1][i]=='0')&& list[0][i]!='\0' && list[1][i]!='\0')
		{
			if(list[0][i]=='1')
				iter1->next_fragment();
			if(list[1][i]=='1')
				iter2->next_fragment();
			i++;
		}

		if(list[0][i]=='\0' || list[1][i]=='\0')
			end=true;

		if(!end)
		{
			frag1=iter1->get_fragment();
			frag2=iter2->get_fragment();

			//fprintf(stderr,"%d - %d\n",iter1->pos_fragment,iter2->pos_fragment);
			end_atom=end_atom2=true;
			//if(strcmp(frag1->getName(),frag2->getName())==0)
			while(end_atom && end_atom2)
			{

				// get atoms
				at = (Atom*)(frag1->getCurrent());
				at2 = (Atom*)(frag2->getCurrent());

				//fprintf(stderr,"Atomos> comparando  %s con  %s\n",at->getName(),at2->getName());
				// get positions
				at->getPosition(pos1);
				at2->getPosition(pos2);

				tmp=0;
				for ( a = 0; a < 3; a++ )
				{
					tmp += (pos1[a] - pos2[a])*(pos1[a] - pos2[a]);
				}

				natoms++;
				suma += tmp;
				cont++;

				end_atom = frag1->next(); end_atom2 = frag2->next();


			}

			i++;
			iter1->next_fragment();
			iter2->next_fragment();

			if(list[0][i]=='\0' || list[1][i]=='\0')
				end=true;
		}
	}while(!end);

	fprintf(stderr,"cont=%d\n",cont);
	delete(iter1);
	delete(iter2);

	if ( suma != 0 ) return sqrt(( float )suma / ( float )natoms) ;
	else
		return ( ( float )suma );

}

// Mon made (17/7/2012)
// RMSD between to molecules with different number of residues
// @param mol2: Reference to the second molecule
// @param mask1: Boolean mask containing matching atoms information for molecule 1
// @param mask2: Boolean mask containing matching atoms information for molecule 2
float Macromolecule::rmsd( Macromolecule *mol2, bool *mask1, bool *mask2 )
{
	bool debug = false;
	Atom *at1, *at2;
	Tcoor pos1, pos2;
	double suma = 0.0;
	int ncommon = 0;

	pdbIter *iter1,*iter2;
	iter1 = new pdbIter(this,false);
	iter2 = new pdbIter(mol2,false);

	iter2->pos_atom=0;
	for(iter1->pos_atom=0; !iter1->gend_atom(); iter1->next_atom())
	{
		if(mask1[iter1->pos_atom]) // Mol 1 matching atom
		{
			//			while(!mask2[iter2->pos_atom] && !iter2->gend_atom()) // this places the index into the corresponding atom
			//			fprintf(stderr,"iter2->pos_atom= %d\n",iter2->pos_atom);
			while(!iter2->gend_atom() && !mask2[iter2->pos_atom]) // this places the index into the corresponding atom
				iter2->next_atom();

			if(!iter2->gend_atom())
			{
				// get atoms
				at1 = iter1->get_atom();
				at2 = iter2->get_atom();

				// get positions
				at1->getPosition(pos1);
				at2->getPosition(pos2);

				if(debug)
					printf("Getting atoms %d (%s) and %d (%s)\n",iter1->pos_atom,at1->getName(),iter2->pos_atom,at2->getName());

				// Computing RMSD
				suma += (pos1[0]-pos2[0])*(pos1[0]-pos2[0]) + (pos1[1]-pos2[1])*(pos1[1]-pos2[1]) + (pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
				ncommon++;

				iter2->next_atom();
			}
		}
	}

	delete(iter1);
	delete(iter2);

	if(debug)
		printf("ncommon %d    suma %f\n", ncommon, suma);

	if ( suma != 0.0 )
		return (float) sqrt( suma / ( double ) ncommon );
	else
		return ( 0.0 );
}

float Macromolecule::rmsd_file( Macromolecule *mol2, bool *mask1, bool *mask2, char *name )
{

	Atom *at1, *at2;
	Tcoor pos1, pos2;
	double suma = 0.0;

    FILE *file;

	pdbIter *iter1,*iter2;
	iter1 = new pdbIter(this,false);
	iter2 = new pdbIter(mol2,false);

	file=fopen(name,"wt");
	if(file==NULL) {
	     fprintf(stderr, " Can not open  %s\n",name);
		return false;
	}

	iter2->pos_atom=0;
	for(iter1->pos_atom=0; !iter1->gend_atom(); iter1->next_atom())
	{
		if(mask1[iter1->pos_atom]) // Mol 1 matching atom
		{
			//			while(!mask2[iter2->pos_atom] && !iter2->gend_atom()) // this places the index into the corresponding atom
			//			fprintf(stderr,"iter2->pos_atom= %d\n",iter2->pos_atom);
			while(!iter2->gend_atom() && !mask2[iter2->pos_atom]) // this places the index into the corresponding atom
				iter2->next_atom();

			if(!iter2->gend_atom())
			{
				// get atoms
				at1 = iter1->get_atom();
				at2 = iter2->get_atom();

				// get positions
				at1->getPosition(pos1);
				at2->getPosition(pos2);
				// Computing RMSD
				suma = (pos1[0]-pos2[0])*(pos1[0]-pos2[0]) + (pos1[1]-pos2[1])*(pos1[1]-pos2[1]) + (pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
				// 	printf("Getting atoms %d (%s) and %d (%s)\n",iter1->pos_atom,at1->getName(),iter2->pos_atom,at2->getName());

//				fprintf(file, "%d %.3f\n",iter1->pos_atom, sqrt(suma) );
				fprintf(file, "%d %.3f %d %s %d %s\n",iter1->pos_atom, sqrt(suma),
						((Residue*)(at1->getFather()))->getIdNumber(),
						((Residue*)(at1->getFather()))->getName(),
						((Residue*)(at2->getFather()))->getIdNumber(),
						((Residue*)(at2->getFather()))->getName()
				);
				iter2->next_atom();
			}
		}
	}

	delete(iter1);
	delete(iter2);
    fclose(file);


   return ( 0.0 );
}

// Mon made (8/5/2008)
// The individual squared deviations will be stored into Bfactor PDB register
// (Inside first Macromolecule, not "mol2")
float Macromolecule::rmsd_bf( Macromolecule * mol2 )
{
	bool end = true, end2 = true;
	Atom * at, * at2;
	Tcoor pos1, pos2;
	double suma = 0, tmp, max=0, tmp_sqrt;
	int a, natoms;

	initAll();
	mol2->initAll();
	natoms = 0;
	while ( end != false && end2!=false )
	{

		/* get atoms */
		at = getCurrentAtom();
		at2 = mol2->getCurrentAtom();

		/* get positions */
		at->getPosition(pos1);
		at2->getPosition(pos2);

		tmp=0;
		for ( a = 0; a < 3; a++ )
		{
			tmp += (pos1[a] - pos2[a])*(pos1[a] - pos2[a]);
		}
		tmp_sqrt = sqrt(tmp);
		at2->setPdbfact( tmp_sqrt ); // Storing Root Square Deviation (RSD) inside B-factor register
		if(tmp_sqrt > max) max = tmp_sqrt; // Stores the maximum RSD
		suma += tmp;

		end = nextAtom(); end2 = mol2->nextAtom();
		natoms++;
	}

	mol2->initAll();
	end2 = true;

	/* (RASMOL NORMALIZES!)
  // RSD (B-factor) Normalization (0-99.99)
//  max = sqrt(max);
  while ( end2 != false )
  {
    at2 = mol2->getCurrentAtom();
  	at2->setPdbfact( 99.99 * ( at2->getPdbfact() / max ) );
//  	printf("Setting B-fact to: %f\n",at2->getPdbfact() );
    end2 = mol2->nextAtom();
  }
	 */

	//printf("natoms %d %f\n", natoms, suma);

	if ( suma != 0 ) return sqrt( ( float )suma / ( float )natoms );
	else
		return ( ( float )suma );
}

// MON: Weighted RMSD computation
float Macromolecule::wrmsd2(Macromolecule *mol2, double *weights)
{
	pdbIter *it,*it2;
	int num,num2;
	Tcoor pos,pos2;
	double dist2;
	double diff[3];
	double sum=0.0;
	double sum_w=0.0;

	it = new pdbIter(this);
	it2 = new pdbIter(mol2);

	// some checking
	num = it->num_atom();
	num2 = it2->num_atom();
	fprintf(stderr,"Msg(wrmsd2): mol= %d mol2= %d\n",num,num2);
	if(num != num2)
	{
		fprintf(stderr,"Msg(wrmsd2): Sorry, number of atoms mismatch: mol= %d mol2= %d\n",num,num2);
		exit(1);
	}

	for(it->pos_atom=0,it2->pos_atom=0; !it->gend_atom(); it->next_atom(),it2->next_atom())
	{
		(it->get_atom())->getPosition(pos);
		(it2->get_atom())->getPosition(pos2);
		diff[0] = pos[0] - pos2[0];
		diff[1] = pos[1] - pos2[1];
		diff[2] = pos[2] - pos2[2];
		dist2 = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
		sum += dist2 * weights[it->pos_atom];
		sum_w += weights[it->pos_atom];
	}
	delete it;
	delete it2;

	return sqrt(sum/sum_w);
}


// Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
// allocating profile memory if (*p_profile==NULL).
// Weights are computed according to ec.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
// Damm & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
// Comparison for Flexible Proteins and Predicted Protein Structures".
// Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
void Macromolecule::gaussian_weight(Macromolecule *mol2, double **p_profile, double c)
{
	int num_atoms;
	num_atoms = this->get_num_atoms();

	// some initial checking
	if(num_atoms != mol2->get_num_atoms())
	{
		fprintf(stderr,"Msg(gaussian_weight): Sorry, different number of atoms %d and %d\n",num_atoms,mol2->get_num_atoms());
		exit(1);
	}

	double *profile;
	if(*p_profile == NULL)
	{
		profile = (double *) malloc(sizeof(double) * num_atoms);
		char *msg  = (char *)"gaussian_weight weights profile allocation";
		check_pointer(profile,msg);
		*p_profile = profile;
	}
	else
		profile = *p_profile;

	pdbIter *iter,*iter2;
	Tcoor pos,pos2;
	double d2;
	iter = new pdbIter(this);
	iter2 = new pdbIter(mol2);
	//	fprintf(stderr,"profile:\n");
	for( iter->pos_atom=0,iter2->pos_atom=0; !iter->gend_atom(); iter->next_atom(),iter2->next_atom() )
	{
		(iter->get_atom())->getPosition(pos);
		(iter2->get_atom())->getPosition(pos2);
		d2 = pow(pos[0]-pos2[0],2) + pow(pos[1]-pos2[1],2) + pow(pos[2]-pos2[2],2);
		profile[iter->pos_atom] = exp( -d2/c );
		//		fprintf(stderr,"%f ",profile[iter->pos_atom]);
	}
	//	fprintf(stderr,"\n");
	delete iter;
	delete iter2;
	//	fprintf(stderr,"Msg(gaussian_weight): Success!\n");
}

// Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
// allocating profile memory if (*p_profile==NULL).
// (17/7/2012) --> Taking into account sequence alignment (at atomic level) for different size molecules.
//                 *The Weights profile its relative to the first molecule (mol)
// Weights are computed according to eq.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
// Damm & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
//                   Comparison for Flexible Proteins and Predicted Protein Structures".
//                   Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
void Macromolecule::gaussian_weight(Macromolecule *mol2, double **p_profile, bool *mask1, bool *mask2, double c)
{
	bool debug = false;
	int num_atoms;
	num_atoms = this->get_num_atoms();

	//	// some initial checking
	//	if(num_atoms != mol2->get_num_atoms())
	//	{
	//		fprintf(stderr,"Msg(gaussian_weight): Sorry, different number of atoms %d and %d\n",num_atoms,mol2->get_num_atoms());
	//		exit(1);
	//	}

	double *profile;
	if(*p_profile == NULL)
	{
		profile = (double *) malloc(sizeof(double) * num_atoms);
		char *msg  = (char *)"gaussian_weight weights profile allocation";
		check_pointer(profile,msg);
		*p_profile = profile;
	}
	else
		profile = *p_profile;

	// Initialization (mandatory if masks are provided!!!)
	for(int i=0; i<num_atoms; i++)
		profile[i] = 0.0; // some negative number for marking atoms that should not be considered

	pdbIter *iter1,*iter2;
	Atom *at1,*at2;
	Tcoor pos1,pos2;
	double d2;
	iter1 = new pdbIter(this);
	iter2 = new pdbIter(mol2);
	int ncommon = 0;
	//	fprintf(stderr,"profile:\n");


	//	pdbIter *iter1,*iter2;
	//	iter1 = new pdbIter(this,false);
	//	iter2 = new pdbIter(mol2,false);

	iter2->pos_atom=0;
	for(iter1->pos_atom=0; !iter1->gend_atom(); iter1->next_atom())
	{
		if(mask1[iter1->pos_atom]) // Mol 1 matching atom
		{
			while(!iter2->gend_atom() && !mask2[iter2->pos_atom]) // this places the index into the corresponding atom
				iter2->next_atom();

			if(!iter2->gend_atom())
			{
				// get atoms
				at1 = iter1->get_atom();
				at2 = iter2->get_atom();

				// get positions
				at1->getPosition(pos1);
				at2->getPosition(pos2);

				if(debug)
					printf("Getting atoms %d (%s) and %d (%s)\n",iter1->pos_atom,at1->getName(),iter2->pos_atom,at2->getName());

				// Computing Weighted RMSD profile
				d2 = pow(pos1[0]-pos2[0],2) + pow(pos1[1]-pos2[1],2) + pow(pos1[2]-pos2[2],2);
				profile[iter1->pos_atom] = exp( -d2/c );
				//				profile[ncommon] = exp( -d2/c );

				if(debug)
					printf("profile[%d]= %f\n",iter1->pos_atom,profile[iter1->pos_atom]);
				//					printf("profile[%d]= %f\n",iter1->pos_atom,profile[ncommon]);

				ncommon++;
				iter2->next_atom();
			}
		}
	}

	//	// Some checking
	//	if(num_atoms != iter1->pos_atom)


	delete iter1;
	delete iter2;
	//	fprintf(stderr,"Msg(gaussian_weight): Success!\n");
}

// MON: Converts an alignment mask at residue level (made with "read_aln" function) into one at atomic level (CA)
bool *Macromolecule::maskres2maskatom(bool *mask)
{
	bool debug = false;
	pdbIter *iter_frag,*iter_atom;
	Atom *atom;
	int i,natoms;
	bool *out_mask;

	// Creating mask
	natoms = this->get_num_atoms();

	if( !(out_mask = (bool *) malloc(sizeof(bool) * natoms) ) )
	{
		printf("Msg(maskres2maskatom): Memory allocation failed!\nForcing exit!\n");
		exit(1);
	}

	// Iterator to iterate through fragments
	iter_frag = new pdbIter(this);

	i = 0; // output mask index
	for(iter_frag->pos_fragment=0; !iter_frag->gend_fragment(); iter_frag->next_fragment())
	{
		// Iterator to iterate through fragment's atoms
		iter_atom = new pdbIter( iter_frag->get_fragment() );

		for(iter_atom->pos_atom=0; !iter_atom->gend_atom(); iter_atom->next_atom())
		{
			atom = iter_atom->get_atom();
			//			printf("%4d  %s\n",i,atom->getName());
			if(strncmp(atom->getName()," CA ",4) == 0 || strncmp(atom->getName()," P  ",4) == 0)
				if(mask[iter_frag->pos_fragment])
					out_mask[i] = true;
				else
					out_mask[i] = false;
			else
				out_mask[i] = false;
			i++;
		}
		delete(iter_atom);
	}
	delete(iter_frag);

	if(natoms != i)
	{
		printf("Msg(maskres2maskatom): Number of atoms mismatch between alignment mask (%d) and macromolecule (%d)!\nForcing exit!\n",natoms,i);
		exit(2);
	}

	// Some checking...
	if(debug)
	{
		printf("\nAtom level mask:\n");
		for(i=0; i<natoms; i++)
		{
			if(out_mask[i])
				printf("1 ");
			else
				printf("0 ");
		}
		printf("\n");
	}

	return(out_mask);
}


bool Macromolecule::collision(Macromolecule *mol,float limit,float acc)
{
	pdbIter *iter1,*iter2;
	Tcoor pos1,pos2;
	float sum1,sum2,sum3;
	float limit2=limit*limit;

	iter1=new pdbIter(this);
	iter2=new pdbIter(mol);

	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		if(iter1->get_atom()->getPdbocc()>acc)
		{
			(iter1->get_atom())->getPosition(pos1);
			for(iter2->pos_atom=0;!iter2->gend_atom();iter2->next_atom())
			{
				if(iter2->get_atom()->getPdbocc()>acc)
				{
					(iter2->get_atom())->getPosition(pos2);
					sum1=pos1[0] - pos2[0];
					sum2=pos1[1] - pos2[1];
					sum3=pos1[2] - pos2[2];

					if((sum1*sum1+sum2*sum2+sum3*sum3)<limit2)
						return true;
				}
			}
		}
	}
	delete iter1;
	delete iter2;

	return false;
}

int Macromolecule::count_clashes(Macromolecule *mol,float limit)
{
	pdbIter *iter1,*iter2;
	Tcoor pos1,pos2;
	float sum1,sum2,sum3;
	float limit2=limit*limit;
	int cont_clashes=0;

	iter1=new pdbIter(this);
	iter2=new pdbIter(mol);

	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		if(iter1->get_atom()->getElement()->symbol != H)
		{
			(iter1->get_atom())->getPosition(pos1);

			for(iter2->pos_atom=0;!iter2->gend_atom();iter2->next_atom())
			{
				if(iter2->get_atom()->getElement()->symbol != H)
				{
					(iter2->get_atom())->getPosition(pos2);
					sum1=pos1[0] - pos2[0];
					sum2=pos1[1] - pos2[1];
					sum3=pos1[2] - pos2[2];

					if((sum1*sum1+sum2*sum2+sum3*sum3)<limit2)
						cont_clashes++;
				}
			}
		}
	}
	delete iter1;
	delete iter2;

	return cont_clashes;
}

float Macromolecule::minDistance(Macromolecule *mol,int *index1, int *index2)
{
	pdbIter *iter1,*iter2;
	Tcoor pos1,pos2;
	float sum1,sum2,sum3;
	float limit=10000000;
	float dist;


	iter1=new pdbIter(this);
	iter2=new pdbIter(mol);

	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		(iter1->get_atom())->getPosition(pos1);
		for(iter2->pos_atom=0;!iter2->gend_atom();iter2->next_atom())
		{
			(iter2->get_atom())->getPosition(pos2);
			sum1=pos1[0] - pos2[0];
			sum2=pos1[1] - pos2[1];
			sum3=pos1[2] - pos2[2];
			dist=(sum1*sum1+sum2*sum2+sum3*sum3);
			if(dist<limit)
			{
				limit=dist;
				//*index1=iter1->pos_atom;
				*index1=(iter1->get_atom())->getPdbSerial();
				//*index2=iter2->pos_atom;
				*index2=(iter2->get_atom())->getPdbSerial();

			}
		}
	}

	delete iter1;
	delete iter2;

	limit=sqrt(limit);
	return(limit);

}


// Computes de pair-pair matrix distance
void Macromolecule::distanceMatrix( double **dmatrix)
{
	pdbIter *iter1,*iter2;
	Tcoor pos1,pos2;
	int index=0,natoms;

	iter1=new pdbIter(this);
	iter2=new pdbIter(this);
	natoms=(int) iter1->num_atom();

	// std::cout<<"atomos "<<index<<std::endl;

	double *matriz;

	// Changed by MON (23/4/2007)
	/*  matriz =(double *) malloc(sizeof(double)*(natoms*(natoms+1)/2));
	memset(matriz,0,sizeof(double)*(natoms*(natoms+1)/2)); */
	matriz =(double *) malloc(sizeof(double)*(natoms*(natoms-1)/2));
	memset(matriz,0,sizeof(double)*(natoms*(natoms-1)/2));

	double sum1,sum2,sum3;
	index=0;
	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		(iter1->get_atom())->getPosition(pos1);
		// i=j
		index=natoms*iter1->pos_atom-1-iter1->pos_atom*(iter1->pos_atom+3)/2;
		for(iter2->pos_atom=iter1->pos_atom+1;!iter2->gend_atom();iter2->next_atom())
		{
			(iter2->get_atom())->getPosition(pos2);
			sum1=pos1[0] - pos2[0];
			sum2=pos1[1] - pos2[1];
			sum3=pos1[2] - pos2[2];
			matriz[index+iter2->pos_atom] = sqrt( sum1*sum1 + sum2*sum2 + sum3*sum3);
			//index++;
		}


	}

	delete iter1;
	delete iter2;

	*dmatrix=matriz;
	//std::cout<<"pair-pair"<<index<<std::endl;

	// return maxL;
}

// Save coordinates in a single row pointer
void Macromolecule::coordMatrix(float **coord)
{
	pdbIter *iter1;
	Tcoor pos1;
	int index=0;

	iter1=new pdbIter(this);
	index=(int) iter1->num_atom();


	float *cxyz;
	cxyz =(float *) malloc(index*3*sizeof(float));



	index=0;
	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		(iter1->get_atom())->getPosition(pos1);
		cxyz[index++]=pos1[0];
		cxyz[index++]=pos1[1];
		cxyz[index++]=pos1[2];
	}

	delete iter1;

	*coord=cxyz;

	// return maxL;
}

// Set macromolecule coordinates from single row coordinates
void Macromolecule::coordMatrixSet(float *coord)
{
	pdbIter *iter1;
	Tcoor pos1;

	iter1 = new pdbIter(this);

	for(iter1->pos_atom=0; !iter1->gend_atom(); iter1->next_atom())
	{
		pos1[0] = coord[iter1->pos_atom*3 ];
		pos1[1] = coord[iter1->pos_atom*3 + 1];
		pos1[2] = coord[iter1->pos_atom*3 + 2];
		(iter1->get_atom())->setPosition(pos1);
	}

	delete iter1;
}

// get Atom property in a single row pointer
void Macromolecule::getPtrAtomProperty(float **coord)
{
	pdbIter *iter1;
	int index=0;

	iter1=new pdbIter(this);
	index=(int) iter1->num_atom();


	float *P;
	P =(float *) malloc(index*sizeof(float));



	index=0;
	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		P[index++]=(iter1->get_atom())->getElement()->number;


	}

	delete iter1;

	*coord=P;

	// return maxL;
}



void Macromolecule::getPtrAtomPropertyBFS(float **coord)
{
	int index=0;
	float BFS[N_AMINO][14];
	pdbIter *iter_frag; // Iter to screen all fragments
	pdbIter *iter_res; // Iter to screen all segments
	Residue * res;
	Atom * atom;
	char * at_name;
	int resn;
	int cont;
	cont=0;
	float numberE;
	float BFS_L, BFS_VL, BFS_M, BFS_H;

	iter_frag = new pdbIter( this );
	index=(int) iter_frag->num_atom();

	float *P;
	P =(float *) malloc(index*sizeof(float));

	for (int k = 0; k < 30; k++ )
		for (int l = 0; l < 14; l++ )
			BFS[k][l]=0;

	// lateral chain factors
	BFS_VL=0.02;
	BFS_L= 0.20;
	BFS_M= 0.30;
	BFS_H= 0.40;

	for (int k = 0; k < 30; k++ )
	{

		BFS[k][0]=1.2; // N
		BFS[k][1]=1.2; // CA
		BFS[k][2]=1.2; // C
		BFS[k][3]=0.3; // 0
		if (k!=5)  // not GLY
		{
			BFS[k][4]=0.5; // CB

			switch (k)
			{
			// less
			case 14: case 8: case 2: case 3: case 10:
				// ARG (14) LYS(8)  ASP(2) GLU (3) MET (10)
				for(int n=5;n<AA[k].natoms;n++)
					BFS[k][n]=BFS_L;
				break;
				// medium
			case 1: case 21: case 13:  case 11: case 16: case 15:
				// 1 CYS, 21 CYX, 13 GLN, 6 HIS, 25 HID, 26 HIE, 11 ASN, 16 THR, 15 SER
				for(int n=5;n<AA[k].natoms;n++)
					BFS[k][n]=BFS_M;
				break;
			default:
				// 19 TYR, 18 TRP, 0 ALA, 4 PHE, 12 PRO, 7 ILE, 9  LEU, 17 VAL
				for(int n=5;n<AA[k].natoms;n++)
					BFS[k][n]=BFS_H;
				break;
			}
		}

	}
	// oxigens.....
	BFS[2][6]=BFS_VL; // ASP OD1
	BFS[2][7]=BFS_VL; // ASP OD2
	BFS[3][7]=BFS_VL; // GLU OE1
	BFS[3][8]=BFS_VL; // ASP OE2
	BFS[13][7]=BFS_VL; // GLN OE1
	BFS[11][6]=BFS_VL; // ASN OD1
	BFS[16][5]=BFS_VL;  // THR OE1
	BFS[15][5]=BFS_VL;  // SER OD1
	BFS[19][11]=BFS_VL; // TYR OD1

	cont=0;
	for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
	{
		res = ( Residue * ) iter_frag->get_fragment();
		iter_res = new pdbIter(res);

		if(res->getClass() != pdb_residue) {

			if(res->getClass() == pdb_nucleotide) {
				/// NUCLEOTIDS
				for(iter_res->pos_atom=0; !iter_res->gend_atom(); iter_res->next_atom() ) {
					atom = ( Atom * ) iter_res->get_atom();
					at_name = atom->getName();
					numberE=atom->getElement()->number;
					if (numberE == 8 ) // O
						P[cont]=0.3;
					else if (numberE == 15 ) // P
						P[cont]=2.0; //
					else if (numberE == 7 ) // N
						P[cont]=1.4;
					else if (numberE == 6 ) // C
						P[cont]=1.4;
					else
						P[cont]=0.00;
					//fprintf(stdout, "+--> %s %s %f\n",res->getName(), at_name, atom->getElement()->number);
					//fprintf(stdout, "+--> %d %s %s %f\n",cont, res->getName(), at_name, P[cont]);
					cont++;
				}
			}
			else 	{

				/// HETERO
				for(iter_res->pos_atom=0; !iter_res->gend_atom(); iter_res->next_atom() ) {
					atom = ( Atom * ) iter_res->get_atom();
					at_name = atom->getName();
					numberE=atom->getElement()->number;
					if (numberE == 8 ) // O
						P[cont]=0.1;
					else if (numberE == 15 ) // P
						P[cont]=2.0; //
					else if (numberE == 7 ) // N
						P[cont]=1.2;
					else if (numberE == 6 ) // C
						P[cont]=1.2;
					else if (numberE == 16 ) // S
						P[cont]=2.0;
					else if (numberE == 20 ) // CA
						P[cont]=0.1;
					else if (numberE == 26 ) // FE
						P[cont]=0.1;
					else if (numberE == 12 ) // MG
						P[cont]=0.1;
					else if (numberE == 25 ) // MN
						P[cont]=0.1;
					else if (numberE == 11 ) // NA
						P[cont]=0.1;
					else if (numberE == 30 ) // Zn
						P[cont]=0.1;
					else if (numberE == 28)  // Ni
						P[cont]=0.1;
					else if (numberE == 29 ) // Cu
						P[cont]=0.1;
					else if (numberE == 19.0 ) // K
						P[cont]=0.1;
					else
						P[cont]=0.35;
					//fprintf(stdout, "---> %d %s %s %f\n",cont, res->getName(), at_name, P[cont]);
					cont++;
				}
			}
		}
		else  // PROTEIN
		{
			resn = resnum_from_resname( res->getName() );
			for(iter_res->pos_atom=0; !iter_res->gend_atom(); iter_res->next_atom() ) {
				atom = ( Atom * ) iter_res->get_atom();
				at_name = atom->getName();
				int detect=0;
				for(int n=0;n<AA[resn].natoms;n++)
				{
					if(strcmp(at_name,AA[resn].atom[n].atom_name)==0)
					{
						P[cont]=BFS[resn][n];
						cont++;detect=1;
						break;
					}
				}
				if (detect==0) {
					P[cont]=0;  // if not detected (??)
					if(strcmp(at_name," OXT")==0) {
						P[cont]=BFS_VL;
						//fprintf(stdout, "#--> %d %s %s %f\n",cont, res->getName(), at_name, P[cont]);
					}
					cont++;
				}
			}
		}
		delete iter_res;
	}
	delete iter_frag;

	*coord=P;

	// return maxL;
}




// Save coordinates in a single row pointer
int Macromolecule::pdbmatrices(float **coord, int **nres, int **pfirst, int **cas, int **tpatom)
{
	pdbIter *iter1, *iter2, *iter0, *iter3;
	Tcoor pos1;
	int natom=0, numres=0, index=0;
	Fragment * resi;
	char * at_name, * resn;
	Segment * seg;
	Atom * atomi;


	iter0=new pdbIter(this);
	natom=(int) iter0->num_atom();
	delete iter0;

	float *cxyz;
	cxyz =(float *) malloc(natom*3*sizeof(float));

	int *Tpatom;
	Tpatom =(int *) malloc(natom*sizeof(int));


	//Iterador para recorrer segmentos
	iter1 = new pdbIter( this );

	numres=0; index=0;
	//Bucle para recorrer segmentos
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );
		fprintf( stderr, "Processing aa %s  %d\n", resi->getName(), iter2->pos_fragment );

		//Bucle para recorrer residuos
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			resi = ( Fragment * ) iter2->get_fragment();
			//resn = resnum_from_resname( resi->getName() );
			// fprintf( stderr, "Processing aa %s resn= %d  %d\n", resi->getName(), resn, iter2->pos_fragment );
			iter3 = new pdbIter( resi );
			numres++;

			for ( iter3->pos_atom = 0; !iter3->gend_atom(); iter3->next_atom() )
			{
				atomi = ( Atom * ) iter3->get_atom();
				//at_name = atomi->getName();
				atomi->getPosition(pos1);
				cxyz[index++]=pos1[0];
				cxyz[index++]=pos1[1];
				cxyz[index++]=pos1[2];

			}
			delete iter3;
		}
		delete iter2;
	}

	int *Nres;
	Nres =(int *) malloc(numres*sizeof(int));

	int *Pfirst;
	Pfirst   =(int *) malloc(numres*sizeof(int));

	int *Cas;
	Cas =(int *) malloc(numres*sizeof(int));


	numres=-1; index=0; int numatom=-1;
	//Bucle para recorrer segmentos
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );
		//Bucle para recorrer residuos
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			resi = ( Fragment * ) iter2->get_fragment();
			//resn = resnum_from_resname( resi->getName() );
			resn = resi->getName() ;
			// fprintf( stderr, "Processing aa %s resn= %d  %d\n", resi->getName(), resn, iter2->pos_fragment );
			iter3 = new pdbIter( resi );

			// first atom
			numres++;
			Pfirst[numres]=index;

			int natom_res=0;
			for ( iter3->pos_atom = 0; !iter3->gend_atom(); iter3->next_atom() )
			{
				atomi = ( Atom * ) iter3->get_atom();
				at_name = atomi->getName();
				numatom++;
				natom_res++;
				Tpatom[numatom]=18;
				if ( strcmp(at_name," N  ") == 0  )
				{
					Tpatom[numatom]=0;
				}
				else if ( strcmp(at_name," CA ") == 0  )
				{
					Cas[numres]=index;
					if(strcmp(resn,"GLY")==0)
						Tpatom[numatom]=4;
					else
						Tpatom[numatom]=1;
				}
				else if ( strcmp(at_name," C  ") == 0  )
				{
					Tpatom[numatom]=2;
				}
				else if ( strcmp(at_name," O  ") == 0  )
				{
					Tpatom[numatom]=3;
				}
				else if ( strcmp(at_name," CB ") == 0  )
				{
					if(strcmp(resn,"SER")==0)
						Tpatom[numatom]=12;
					else
						Tpatom[numatom]=5;
				}
				else if ( strcmp(at_name," CE ") == 0  )
				{
					if(strcmp(resn,"LYS")==0)
						Tpatom[numatom]=6;
					else if(strcmp(resn,"MET")==0)
						Tpatom[numatom]=16;
				}
				else if ( strcmp(at_name," NZ ") == 0  )
				{
					if(strcmp(resn,"LYS")==0)
						Tpatom[numatom]=6;
				}
				else if ( strcmp(at_name," CD ") == 0  )
				{
					if(strcmp(resn,"PRO")==0)
						Tpatom[numatom]=5;
					else if(strcmp(resn,"LYS")==0)
						Tpatom[numatom]=7;
					else if ((strcmp(resn,"GLU")==0)||(strcmp(resn,"GLH")==0))
						Tpatom[numatom]=8;
					else if(strcmp(resn,"GLN")==0)
						Tpatom[numatom]=10;
					else if(strcmp(resn,"ARG")==0)
						Tpatom[numatom]=11;
					else if(strcmp(resn,"ILE")==0)
						Tpatom[numatom]=16;
				}
				else if ( strcmp(at_name," CG ") == 0  )
				{
					if(strcmp(resn,"PRO")==0)
						Tpatom[numatom]=5;
					else if((strcmp(resn,"ASP")==0)||(strcmp(resn,"ASH")==0))
						Tpatom[numatom]=8;
					else if(strcmp(resn,"ASN")==0)
						Tpatom[numatom]=10;
					else if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
						Tpatom[numatom]=13;
					else if(strcmp(resn,"ARG")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"GLN")==0)
						Tpatom[numatom]=15;
					else if ((strcmp(resn,"GLU")==0)||(strcmp(resn,"GLH")==0))
						Tpatom[numatom]=15;
					else if(strcmp(resn,"LEU")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"LYS")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"MET")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"PHE")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"TYR")==0)
						Tpatom[numatom]=15;
				}
				else if(strcmp(at_name," OD1") == 0 )
				{
					if((strcmp(resn,"ASP")==0)||(strcmp(resn,"ASH")==0))
						Tpatom[numatom]=8;
					else if(strcmp(resn,"ASN")==0)
						Tpatom[numatom]=10;
				}
				else if(strcmp(at_name," OD2") == 0  )
				{
					if((strcmp(resn,"ASP")==0)||(strcmp(resn,"ASH")==0))
						Tpatom[numatom]=8;
				}
				else if(strcmp(at_name," OE1") == 0  )
				{
					if ((strcmp(resn,"GLU")==0)||(strcmp(resn,"GLH")==0))
						Tpatom[numatom]=8;
					else if(strcmp(resn,"GLN")==0)
						Tpatom[numatom]=10;
				}
				else if(strcmp(at_name," OE2") == 0  )
				{
					if ((strcmp(resn,"GLU")==0)||(strcmp(resn,"GLH")==0))
						Tpatom[numatom]=8;
				}
				else if(strcmp(at_name," CZ ") == 0  )
				{
					if(strcmp(resn,"ARG")==0)
						Tpatom[numatom]=9;
					else if(strcmp(resn,"TYR")==0)
						Tpatom[numatom]=14;
					else if(strcmp(resn,"PHE")==0)
						Tpatom[numatom]=15;
				}
				else if( (strcmp(at_name," NH1") == 0 )||( strcmp(at_name," NH2") == 0)  )
				{
					if(strcmp(resn,"ARG")==0)
						Tpatom[numatom]=9;
				}
				else if(strcmp(at_name," ND2") == 0  )
				{
					if(strcmp(resn,"ASN")==0)
						Tpatom[numatom]=10;
				}
				else if(strcmp(at_name," NE ") == 0  )
				{
					if(strcmp(resn,"ARG")==0)
						// Tpatom[numatom]=10;
						Tpatom[numatom]=11;
				}
				else if(strcmp(at_name," OG ") == 0  )
				{
					if(strcmp(resn,"SER")==0)
						Tpatom[numatom]=12;
				}
				else if(strcmp(at_name," OG1") == 0  )
				{
					if(strcmp(resn,"THR")==0)
						Tpatom[numatom]=12;
				}
				else if(strcmp(at_name," OH ") == 0  )
				{
					if(strcmp(resn,"TYR")==0)
						Tpatom[numatom]=12;
				}
				else if(strcmp(at_name," ND1") == 0  )
				{
					if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
						Tpatom[numatom]=13;
				}
				else if(strcmp(at_name," CD1") == 0  )
				{
					if(strcmp(resn,"PHE")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"TYR")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"LEU")==0)
						Tpatom[numatom]=16;
					// This CD1 from ILE (type 16) is not in Table 1, but seems required...
					else if(strcmp(resn,"ILE")==0)
						Tpatom[numatom]=16;
				}
				else if(strcmp(at_name," CD2") == 0  )
				{
					if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
						Tpatom[numatom]=13;
					else if(strcmp(resn,"PHE")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"TYR")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"LEU")==0)
						Tpatom[numatom]=16;
				}
				else if(strcmp(at_name," NE2") == 0  )
				{
					if(strcmp(resn,"GLN")==0)
						Tpatom[numatom]=10;
					else if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
						Tpatom[numatom]=13;
				}
				else if(strcmp(at_name," NE1") == 0  )
				{
					if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=13;
				}
				else if(strcmp(at_name," CE1") == 0  )
				{
					if(strcmp(resn,"TYR")==0)
						Tpatom[numatom]=14;
					else if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
						Tpatom[numatom]=13;
					else if(strcmp(resn,"PHE")==0)
						Tpatom[numatom]=15;
				}
				else if(strcmp(at_name," CE2") == 0  )
				{
					if(strcmp(resn,"TYR")==0)
						Tpatom[numatom]=14;
					else if(strcmp(resn,"PHE")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=15;
				}
				else if(strcmp(at_name," CE3") == 0  )
				{
					if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=15;
				}
				else if(strcmp(at_name," CZ2") == 0  )
				{
					if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=15;
				}
				else if(strcmp(at_name," CZ3") == 0  )
				{
					if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=15;
				}
				else if(strcmp(at_name," CG1") == 0  )
				{
					if(strcmp(resn,"ILE")==0)
						Tpatom[numatom]=15;
					if(strcmp(resn,"VAL")==0)
						Tpatom[numatom]=16;
				}
				else if(strcmp(at_name," SD ") == 0  )
				{
					if(strcmp(resn,"MET")==0)
						Tpatom[numatom]=15;
				}
				else if(strcmp(at_name," CG2") == 0  )
				{
					if(strcmp(resn,"THR")==0)
						Tpatom[numatom]=15;
					else if(strcmp(resn,"ILE")==0)
						Tpatom[numatom]=16;
					else if(strcmp(resn,"VAL")==0)
						Tpatom[numatom]=16;
				}
				else if(strcmp(at_name," CH2") == 0  )
				{
					if(strcmp(resn,"TRP")==0)
						Tpatom[numatom]=15;
				}
				else if(strcmp(at_name," SG ") == 0  )
				{
					if((strcmp(resn,"CYS")==0)||(strcmp(resn,"CYX")==0)||(strcmp(resn,"CYM")==0))
						Tpatom[numatom]=17;
				}
				else if(strcmp(at_name," OXT") == 0  )
				{
					Tpatom[numatom]=3;
				}
				if (Tpatom[numatom]==18)
					fprintf(stderr," RES %s atomo %s typ %d    %d\n",  resn,  at_name, Tpatom[numatom], numatom);
				index+=3;
			}

			// Number
			Nres[numres]= natom_res;

			delete iter3;
		}
		delete iter2;
	}

	// fprintf(stderr," coord 0 %f %f %f\n",  cxyz[0],  cxyz[1], cxyz[2]);
	// for (int i1 = 0; i1 < numatom; i1++ ) {
	// 	fprintf(stderr," %d %d\n", Tpatom[i1], i1);
	// }

	// for (int i1 = 0; i1 < numres; i1++ ) {
	// 	fprintf(stderr,"-%d %d %d\n", Pfirst[i1], i1, Nres[i1]);
	// }
	// pos CA Ligand
	// px = *(cxyz+Cas[i1]);
	// py = *(cxyz+Cas[i1]+1);
	// pz = *(cxyz+Cas[i1]+2);

	// fprintf(stderr," res %d %d %d %d\n",  i1, Cas[i1],  Nres[i1], Pfirst[i1]);
	// getchar();
	//  }

	delete iter1;

	*coord=cxyz;
	*nres = Nres;
	*pfirst = Pfirst;
	*cas = Cas;
	*tpatom = Tpatom;

	return numres;
}


int Soap[30][14];


void soap_potential()
{
	int p;

	p=0; //ALA
	Soap[p][0] =   0; //    N
	Soap[p][1] =   1; //    CA
	Soap[p][2] =   2; //    C
	Soap[p][3] =   3; //    O
	Soap[p][4] =   4; //    CB

	p=14; //ARG
	Soap[p][0] =   5; //    N
	Soap[p][1] =   6; //    CA
	Soap[p][2] =   7; //    C
	Soap[p][3] =   8; //    O
	Soap[p][4] =   9; //    CB
	Soap[p][5] =  10; //    CG
	Soap[p][6] =  11; //    CD
	Soap[p][7] =  12; //    NE
	Soap[p][8] =  13; //    CZ
	Soap[p][9] =  14; //    NH1
	Soap[p][10] = 14; //    NH2

	p=11; //ASN
	Soap[p][0] =  15; //    N
	Soap[p][1] =  16; //    CA
	Soap[p][2] =  17; //    C
	Soap[p][3] =  18; //    O
	Soap[p][4] =  19; //    CB
	Soap[p][5] =  20; //    CG
	Soap[p][6] =  21; //    OD1
	Soap[p][7] =  22; //    ND2

	p=2;  //ASP
	Soap[p][0] =  23; //    N
	Soap[p][1] =  24; //    CA
	Soap[p][2] =  25; //    C
	Soap[p][3] =  26; //    O
	Soap[p][4] =  27; //    CB
	Soap[p][5] =  28; //    CG
	Soap[p][6] =  29; //    OD1
	Soap[p][7] =  29; //    OD2


	p=1; //CYS
	Soap[p][0] =  30; //    N
	Soap[p][1] =  31; //    CA
	Soap[p][2] =  32; //    C
	Soap[p][3] =  33; //    O
	Soap[p][4] =  34; //    CB
	Soap[p][5] =  35; //    SG

	p=21; //CYX
	Soap[p][0] =  30; //    N
	Soap[p][1] =  31; //    CA
	Soap[p][2] =  32; //    C
	Soap[p][3] =  33; //    O
	Soap[p][4] =  34; //    CB
	Soap[p][5] =  35; //    SG

	p=22; //CYM
	Soap[p][0] =  30; //    N
	Soap[p][1] =  31; //    CA
	Soap[p][2] =  32; //    C
	Soap[p][3] =  33; //    O
	Soap[p][4] =  34; //    CB
	Soap[p][5] =  35; //    SG


	p=13;  //GLN
	Soap[p][0] =  36; //    N
	Soap[p][1] =  37; //    CA
	Soap[p][2] =  38; //    C
	Soap[p][3] =  39; //    O
	Soap[p][4] =  40; //    CB
	Soap[p][5] =  41; //    CG
	Soap[p][6] =  42; //    CD
	Soap[p][7] =  43; //    OE1
	Soap[p][8] =  44; //    NE2

	p=3;  //GLU
	Soap[p][0] =  45; //    N
	Soap[p][1] =  46; //    CA
	Soap[p][2] =  47; //    C
	Soap[p][3] =  48; //    O
	Soap[p][4] =  49; //    CB
	Soap[p][5] =  50; //    CG
	Soap[p][6] =  51; //    CD
	Soap[p][7] =  52; //    OE1
	Soap[p][8] =  52; //    OE2


	p=23;  //GLN
	Soap[p][0] =  45; //    N
	Soap[p][1] =  46; //    CA
	Soap[p][2] =  47; //    C
	Soap[p][3] =  48; //    O
	Soap[p][4] =  49; //    CB
	Soap[p][5] =  50; //    CG
	Soap[p][6] =  51; //    CD
	Soap[p][7] =  52; //    OE1
	Soap[p][8] =  52; //    OE2


	p=5;  //GLY
	Soap[p][0] =  53; //    N
	Soap[p][1] =  54; //    CA
	Soap[p][2] =  55; //    C
	Soap[p][3] =  56; //    O

	p=6;  //HIS
	Soap[p][0] =  57; //    N
	Soap[p][1] =  58; //    CA
	Soap[p][2] =  59; //    C
	Soap[p][3] =  60; //    O
	Soap[p][4] =  61; //    CB
	Soap[p][5] =  62; //    CG
	Soap[p][6] =  63; //    ND1
	Soap[p][7] =  64; //    CD2
	Soap[p][8] =  65; //    CE1
	Soap[p][9] =  66; //    NE2

	p=24;  //HIP
	Soap[p][0] =  57; //    N
	Soap[p][1] =  58; //    CA
	Soap[p][2] =  59; //    C
	Soap[p][3] =  60; //    O
	Soap[p][4] =  61; //    CB
	Soap[p][5] =  62; //    CG
	Soap[p][6] =  63; //    ND1
	Soap[p][7] =  64; //    CD2
	Soap[p][8] =  65; //    CE1
	Soap[p][9] =  66; //    NE2


	p=25;  //HID
	Soap[p][0] =  57; //    N
	Soap[p][1] =  58; //    CA
	Soap[p][2] =  59; //    C
	Soap[p][3] =  60; //    O
	Soap[p][4] =  61; //    CB
	Soap[p][5] =  62; //    CG
	Soap[p][6] =  63; //    ND1
	Soap[p][7] =  64; //    CD2
	Soap[p][8] =  65; //    CE1
	Soap[p][9] =  66; //    NE2


	p=26;  //HIE
	Soap[p][0] =  57; //    N
	Soap[p][1] =  58; //    CA
	Soap[p][2] =  59; //    C
	Soap[p][3] =  60; //    O
	Soap[p][4] =  61; //    CB
	Soap[p][5] =  62; //    CG
	Soap[p][6] =  63; //    ND1
	Soap[p][7] =  64; //    CD2
	Soap[p][8] =  65; //    CE1
	Soap[p][9] =  66; //    NE2


	p=7;  //ILE
	Soap[p][0] =  67; //    N
	Soap[p][1] =  68; //    CA
	Soap[p][2] =  69; //    C
	Soap[p][3] =  70; //    O
	Soap[p][4] =  71; //    CB
	Soap[p][5] =  72; //    CG1
	Soap[p][6] =  73; //    CG2
	Soap[p][7] =  74; //    CD1

	p=9;  //LEU
	Soap[p][0] =  75; //    N
	Soap[p][1] =  76; //    CA
	Soap[p][2] =  77; //    C
	Soap[p][3] =  78; //    O
	Soap[p][4] =  79; //    CB
	Soap[p][5] =  80; //    CG
	Soap[p][6] =  81; //    CD1
	Soap[p][7] =  81; //    CD2

	p=8;  //LYS
	Soap[p][0] =  82; //    N
	Soap[p][1] =  83; //    CA
	Soap[p][2] =  84; //    C
	Soap[p][3] =  85; //    O
	Soap[p][4] =  86; //    CB
	Soap[p][5] =  87; //    CG
	Soap[p][6] =  88; //    CD
	Soap[p][7] =  89; //    CE
	Soap[p][8] =  90; //    NZ


	p=27;  //LYN
	Soap[p][0] =  82; //    N
	Soap[p][1] =  83; //    CA
	Soap[p][2] =  84; //    C
	Soap[p][3] =  85; //    O
	Soap[p][4] =  86; //    CB
	Soap[p][5] =  87; //    CG
	Soap[p][6] =  88; //    CD
	Soap[p][7] =  89; //    CE
	Soap[p][8] =  90; //    NZ


	p=10;  //MET
	Soap[p][0] =  91; //    N
	Soap[p][1] =  92; //    CA
	Soap[p][2] =  93; //    C
	Soap[p][3] =  94; //    O
	Soap[p][4] =  95; //    CB
	Soap[p][5] =  96; //    CG
	Soap[p][6] =  97; //    SD
	Soap[p][7] =  98; //    CE


	p=29;  //MSE
	Soap[p][0] =  91; //    N
	Soap[p][1] =  92; //    CA
	Soap[p][2] =  93; //    C
	Soap[p][3] =  94; //    O
	Soap[p][4] =  95; //    CB
	Soap[p][5] =  96; //    CG
	Soap[p][6] =  97; //    SD
	Soap[p][7] =  98; //    CE


	p=4;  //PHE
	Soap[p][0] =  99; //    N
	Soap[p][1] = 100; //    CA
	Soap[p][2] = 101; //    C
	Soap[p][3] = 102; //    O
	Soap[p][4] = 103; //    CB
	Soap[p][5] = 104; //    CG
	Soap[p][6] = 105; //    CD1
	Soap[p][7] = 105; //    CD2
	Soap[p][8] = 106; //    CE1
	Soap[p][9] = 106; //    CE2
	Soap[p][10] =107; //    CZ

	p=12;  //PRO
	Soap[p][0] = 108; //    N
	Soap[p][1] = 109; //    CA
	Soap[p][2] = 110; //    C
	Soap[p][3] = 111; //    O
	Soap[p][4] = 112; //    CB
	Soap[p][5] = 113; //    CG
	Soap[p][6] = 114; //    CD

	p=15;  //SER
	Soap[p][0] = 115; //    N
	Soap[p][1] = 116; //    CA
	Soap[p][2] = 117; //    C
	Soap[p][3] = 118; //    O
	Soap[p][4] = 119; //    CB
	Soap[p][5] = 120; //    OG

	p=16;  //THR
	Soap[p][0] = 121; //    N
	Soap[p][1] = 122; //    CA
	Soap[p][2] = 123; //    C
	Soap[p][3] = 124; //    O
	Soap[p][4] = 125; //    CB
	Soap[p][5] = 126; //    OG1
	Soap[p][6] = 127; //    CG2

	p=18;  //TRP
	Soap[p][0] = 128; //    N
	Soap[p][1] = 129; //    CA
	Soap[p][2] = 130; //    C
	Soap[p][3] = 131; //    O
	Soap[p][4] = 132; //    CB
	Soap[p][5] = 133; //    CG
	Soap[p][6] = 134; //    CD1
	Soap[p][7] = 135; //    CD2
	Soap[p][8] = 136; //    NE1
	Soap[p][9] = 137; //    CE2
	Soap[p][10]= 138; //    CE3
	Soap[p][11]= 139; //    CZ2
	Soap[p][12]= 140; //    CZ3
	Soap[p][13]= 141; //    CH2

	p=19;  //TYR
	Soap[p][0] = 142; //    N
	Soap[p][1] = 143; //    CA
	Soap[p][2] = 144; //    C
	Soap[p][3] = 145; //    O
	Soap[p][4] = 146; //    CB
	Soap[p][5] = 147; //    CG
	Soap[p][6] = 148; //    CD1
	Soap[p][7] = 148; //    CD2
	Soap[p][8] = 149; //    CE1
	Soap[p][9] = 149; //    CE2
	Soap[p][10]= 150; //    CZ
	Soap[p][11]= 151; //    OH


	p=28;  //TYR
	Soap[p][0] = 142; //    N
	Soap[p][1] = 143; //    CA
	Soap[p][2] = 144; //    C
	Soap[p][3] = 145; //    O
	Soap[p][4] = 146; //    CB
	Soap[p][5] = 147; //    CG
	Soap[p][6] = 148; //    CD1
	Soap[p][7] = 148; //    CD2
	Soap[p][8] = 149; //    CE1
	Soap[p][9] = 149; //    CE2
	Soap[p][10]= 150; //    CZ
	Soap[p][11]= 151; //    OH

	p=17;  //VAL
	Soap[p][0] = 152; //    N
	Soap[p][1] = 153; //    CA
	Soap[p][2] = 154; //    C
	Soap[p][3] = 155; //    O
	Soap[p][4] = 156; //    CB
	Soap[p][5] = 157; //    CG1
	Soap[p][6] = 157; //    CG2

}


int tobi_atom(char *at_name, char* resn ) {

	//fprintf(stderr," RES %s atomo \"%s\"\n",  resn,  at_name);
	//if(strcmp(at_name," OXT") == 0  )
	//	fprintf(stderr,"--> RES %s atomo \"%s\"\n",  resn,  at_name);

	if ( strcmp(at_name," N  ") == 0  )
	{
		return 0;
	}
	else if ( strcmp(at_name," CA ") == 0  )
	{
		// Cas[numres]=index;
		if(strcmp(resn,"GLY")==0)
			return 4;
		else
			return 1;
	}
	else if ( strcmp(at_name," C  ") == 0  )
	{
		return 2;
	}
	else if ( strcmp(at_name," O  ") == 0  )
	{
		return 3;
	}
	else if ( strcmp(at_name," CB ") == 0  )
	{
		if(strcmp(resn,"SER")==0)
			return 12;
		else
			return 5;
	}
	else if ( strcmp(at_name," CE ") == 0  )
	{
		if(strcmp(resn,"LYS")==0)
			return 6;
		else if(strcmp(resn,"MET")==0)
			return 16;
	}
	else if ( strcmp(at_name," NZ ") == 0  )
	{
		if(strcmp(resn,"LYS")==0)
			return 6;
	}
	else if ( strcmp(at_name," CD ") == 0  )
	{
		if(strcmp(resn,"PRO")==0)
			return 5;
		else if(strcmp(resn,"LYS")==0)
			return 7;
		else if ((strcmp(resn,"GLU")==0)||(strcmp(resn,"GLH")==0))
			return 8;
		else if(strcmp(resn,"GLN")==0)
			return 10;
		else if(strcmp(resn,"ARG")==0)
			return 11;
		else if(strcmp(resn,"ILE")==0)
			return 16;
	}
	else if ( strcmp(at_name," CG ") == 0  )
	{
		if(strcmp(resn,"PRO")==0)
			return 5;
		else if((strcmp(resn,"ASP")==0)||(strcmp(resn,"ASH")==0))
			return 8;
		else if(strcmp(resn,"ASN")==0)
			return 10;
		else if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
			return 13;
		else if(strcmp(resn,"ARG")==0)
			return 15;
		else if(strcmp(resn,"GLN")==0)
			return 15;
		else if ((strcmp(resn,"GLU")==0)||(strcmp(resn,"GLH")==0))
			return 15;
		else if(strcmp(resn,"LEU")==0)
			return 15;
		else if(strcmp(resn,"LYS")==0)
			return 15;
		else if(strcmp(resn,"MET")==0)
			return 15;
		else if(strcmp(resn,"PHE")==0)
			return 15;
		else if(strcmp(resn,"TRP")==0)
			return 15;
		else if(strcmp(resn,"TYR")==0)
			return 15;
	}
	else if(strcmp(at_name," OD1") == 0 )
	{
		if((strcmp(resn,"ASP")==0)||(strcmp(resn,"ASH")==0))
			return 8;
		else if(strcmp(resn,"ASN")==0)
			return 10;
	}
	else if(strcmp(at_name," OD2") == 0  )
	{
		if((strcmp(resn,"ASP")==0)||(strcmp(resn,"ASH")==0))
			return 8;
	}
	else if(strcmp(at_name," OE1") == 0  )
	{
		if ((strcmp(resn,"GLU")==0)||(strcmp(resn,"GLH")==0))
			return 8;
		else if(strcmp(resn,"GLN")==0)
			return 10;
	}
	else if(strcmp(at_name," OE2") == 0  )
	{
		if ((strcmp(resn,"GLU")==0)||(strcmp(resn,"GLH")==0))
			return 8;
	}
	else if(strcmp(at_name," CZ ") == 0  )
	{
		if(strcmp(resn,"ARG")==0)
			return 9;
		else if(strcmp(resn,"TYR")==0)
			return 14;
		else if(strcmp(resn,"PHE")==0)
			return 15;
	}
	else if( (strcmp(at_name," NH1") == 0 )||( strcmp(at_name," NH2") == 0)  )
	{
		if(strcmp(resn,"ARG")==0)
			return 9;
	}
	else if(strcmp(at_name," ND2") == 0  )
	{
		if(strcmp(resn,"ASN")==0)
			return 10;
	}
	else if(strcmp(at_name," NE ") == 0  )
	{
		if(strcmp(resn,"ARG")==0)
			// return 10;
			return 11;
	}
	else if(strcmp(at_name," OG ") == 0  )
	{
		if(strcmp(resn,"SER")==0)
			return 12;
	}
	else if(strcmp(at_name," OG1") == 0  )
	{
		if(strcmp(resn,"THR")==0)
			return 12;
	}
	else if(strcmp(at_name," OH ") == 0  )
	{
		if(strcmp(resn,"TYR")==0)
			return 12;
	}
	else if(strcmp(at_name," ND1") == 0  )
	{
		if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
			return 13;
	}
	else if(strcmp(at_name," CD1") == 0  )
	{
		if(strcmp(resn,"PHE")==0)
			return 15;
		else if(strcmp(resn,"TRP")==0)
			return 15;
		else if(strcmp(resn,"TYR")==0)
			return 15;
		else if(strcmp(resn,"LEU")==0)
			return 16;
		// This CD1 from ILE (type 16) is not in Table 1, but seems required...
		else if(strcmp(resn,"ILE")==0)
			return 16;
	}
	else if(strcmp(at_name," CD2") == 0  )
	{
		if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
			return 13;
		else if(strcmp(resn,"PHE")==0)
			return 15;
		else if(strcmp(resn,"TRP")==0)
			return 15;
		else if(strcmp(resn,"TYR")==0)
			return 15;
		else if(strcmp(resn,"LEU")==0)
			return 16;
	}
	else if(strcmp(at_name," NE2") == 0  )
	{
		if(strcmp(resn,"GLN")==0)
			return 10;
		else if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
			return 13;
	}
	else if(strcmp(at_name," NE1") == 0  )
	{
		if(strcmp(resn,"TRP")==0)
			return 13;
	}
	else if(strcmp(at_name," CE1") == 0  )
	{
		if(strcmp(resn,"TYR")==0)
			return 14;
		else if ((strcmp(resn,"HIS")==0)||(strcmp(resn,"HID")==0)||(strcmp(resn,"HIP")==0)||(strcmp(resn,"HIE")==0))
			return 13;
		else if(strcmp(resn,"PHE")==0)
			return 15;
	}
	else if(strcmp(at_name," CE2") == 0  )
	{
		if(strcmp(resn,"TYR")==0)
			return 14;
		else if(strcmp(resn,"PHE")==0)
			return 15;
		else if(strcmp(resn,"TRP")==0)
			return 15;
	}
	else if(strcmp(at_name," CE3") == 0  )
	{
		if(strcmp(resn,"TRP")==0)
			return 15;
	}
	else if(strcmp(at_name," CZ2") == 0  )
	{
		if(strcmp(resn,"TRP")==0)
			return 15;
	}
	else if(strcmp(at_name," CZ3") == 0  )
	{
		if(strcmp(resn,"TRP")==0)
			return 15;
	}
	else if(strcmp(at_name," CG1") == 0  )
	{
		if(strcmp(resn,"ILE")==0)
			return 15;
		if(strcmp(resn,"VAL")==0)
			return 16;
	}
	else if(strcmp(at_name," SD ") == 0  )
	{
		if(strcmp(resn,"MET")==0)
			return 15;
	}
	else if(strcmp(at_name," CG2") == 0  )
	{
		if(strcmp(resn,"THR")==0)
			return 15;
		else if(strcmp(resn,"ILE")==0)
			return 16;
		else if(strcmp(resn,"VAL")==0)
			return 16;
	}
	else if(strcmp(at_name," CH2") == 0  )
	{
		if(strcmp(resn,"TRP")==0)
			return 15;
	}
	else if(strcmp(at_name," SG ") == 0  )
	{
		if((strcmp(resn,"CYS")==0)||(strcmp(resn,"CYX")==0)||(strcmp(resn,"CYM")==0))
			return 17;
	}
	else if(strcmp(at_name," OXT") == 0  )
	{
		return 3;
	}

	//fprintf(stderr,"++RES %s atomo %s\n",  resn,  at_name);

	return  18;

}

// Save coordinates in a single row pointer
int Macromolecule::ppimatrices(float **cxyz, int **tpatom, int **spatom,  int **pfirst, int **cas, int **nres)
{
	pdbIter *iter1, *iter2, *iter3;
	pdbIter *iter_prot,*iterAP;
	int *fatomP, *fatomP_tobi, numres, natom, natomres, cont;
	float *coor_atomP;
	Segment * seg;
	Fragment *fragP;
	Fragment * resi;
	Atom *atomP;
	float posP[3];
	char * at_name;
	bool detect_at=false;
	bool first_detected=false;

	soap_potential();


	//Iterador para recorrer segmentos
	iter1 = new pdbIter( this );

	numres=0;  natom=0;
	//Bucle para recorrer segmentos
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );

		//Bucle para recorrer residuos
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			resi = ( Fragment * ) iter2->get_fragment();
			//resn = resnum_from_resname( resi->getName() );
			// fprintf( stderr, "Processing aa %s resn= %d  %d\n", resi->getName(), resn, iter2->pos_fragment );
			iter3 = new pdbIter( resi );
			numres++;
			natom+=iter3->num_atom();
			delete iter3;
		}
		delete iter2;
	}

 	int *Nres;
	Nres =(int *) malloc(numres*sizeof(int));

	int *Pfirst;
	Pfirst   =(int *) malloc(numres*sizeof(int));

	int *Cas;
	Cas =(int *) malloc(numres*sizeof(int));

	fatomP =(int*)malloc(sizeof(int)*natom);  // Soap
	fatomP_tobi =(int*)malloc(sizeof(int)*natom); /// tobi
	coor_atomP=(float*)malloc(sizeof(float*)*natom*3);

    // fprintf( stderr, "Processing  %d  %d\n", numres, natom );

	natom=0; numres=0; cont=0;
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter_prot = new pdbIter( seg );

		for(iter_prot->pos_fragment=0; !iter_prot->gend_fragment(); iter_prot->next_fragment() )
		{
			fragP = iter_prot->get_fragment();
			fragP->getName();
			//fprintf(stderr, " AA=%s\n", fragP->getName());
			iterAP = new pdbIter(fragP); // iterate  atoms
			natomres=iterAP->num_atom();
			

			first_detected=true;
			for (int m=0; m<N_AMINO; m++)
			{
				if (strcmp (fragP->getName(),AA[m].aa_name3)==0)
				{
					for(iterAP->pos_atom=0; !iterAP->gend_atom(); iterAP->next_atom() )
					{
						atomP = iterAP->get_atom();
						atomP->getPosition(posP);
						at_name = atomP->getName();

						detect_at=false;
						for(int n=0;n<AA[m].natoms;n++)
						{
							if(strcmp(at_name,AA[m].atom[n].atom_name)==0)
							{
								fatomP[cont]=Soap[m][n];
								detect_at=true;
								fatomP_tobi[cont]=tobi_atom(atomP->getPdbName(), fragP->getName() );
								break;
							}
						}

						if (detect_at) {
							if (first_detected) {
								Pfirst[numres]=cont;
								Nres[numres]=natomres;
								first_detected=false;
								Cas[numres]=cont;
                                numres++;
							}
							if ( strcmp(at_name," CA ") == 0  ) {
								Cas[numres-1]=cont;
							}
							coor_atomP[cont*3+0]=posP[0];
							coor_atomP[cont*3+1]=posP[1];
							coor_atomP[cont*3+2]=posP[2];
							cont++;

						} else
							if(strcmp(at_name," OXT") == 0 )
							{
								coor_atomP[cont*3+0]=posP[0];
								coor_atomP[cont*3+1]=posP[1];
								coor_atomP[cont*3+2]=posP[2];
								fatomP[cont]=29;
								fatomP_tobi[cont]=3;
								cont++;

							}
							else fprintf(stdout, "frodock_opt> Atom Receptor %s not recognized \n",atomP->getPdbName());

					}
					m=N_AMINO+1;
				}
			}
			delete iterAP;
		}
           delete iter_prot;
	}
	
    //fprintf( stderr, "Processing2  %d  %d\n", numres, cont );

	*cxyz=coor_atomP;
	*tpatom = fatomP_tobi;
	*spatom = fatomP;
	*pfirst = Pfirst;
	*cas = Cas;
	*nres = Nres;
       // fprintf( stderr, " pasado...." );
	return numres;
//	return cont;
}






float Macromolecule::minRmsd( Macromolecule * mol2 )
{
	bool end = true, end2 = true;
	Atom * at, * at2;
	Tcoor pos1, pos2;
	float center1[3], center2[3];
	double cent[3] [3], det, tmp, tmp1, handedness;
	int natoms;
	double m[3][3], aa[3][3],err;
	int i, j;

	// initialize
	for ( i = 0; i < 3; i++ )
	{
		for ( j = 0; j < 3; j++ )
		{
			m[i] [j] = 0.0F;
			aa[i] [j] = 0.0F;
		}
		m[i] [i] = 1.0F; center1[i]=0.0F;center2[i]=0.0F;
	}


	/* if ((natoms=getLimit())!=mol2->getLimit())  { printf("Error different atom numbers %d %d\n", natoms,
  mol2->getLimit()); return -1; } */

	// center pdbs

	geoCenter( center1 );
	mol2->geoCenter( center2 );

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			cent[i] [j] = (double) center1[i] * center2[j];

	/* Calculate correlation matrix */
	initAll();
	mol2->initAll();
	natoms = 0;
	while ( end != false )
	{
		/* get atoms */
		at = getCurrentAtom();
		at2 = mol2->getCurrentAtom();

		/* get positions */
		at->getPosition(pos1);
		at2->getPosition(pos2);

		/* sum (pos1-c1)*(pos2-c2) = sum (pos1*pos2-c1*c2)
    for ( i = 0; i < 3; i++ )
      for ( j = 0; j < 3; j++ )
        aa[i] [j] += (double) pos1[j] * pos2[i] - cent[i] [j];
		 */

		double x1[3], x2[3];

		for ( j = 0; j < 3; j++ )
		{
			x1[j] = (double) pos1[j]-center1[j];
			x2[j] = (double) pos2[j]-center2[j];
		}

		for ( i = 0; i < 3; i++ )
			for ( j = 0; j < 3; j++ )
				aa[i] [j] += x1[i]*x2[j];



		end = nextAtom(); end2 = mol2->nextAtom();
		natoms++;
	}


	/* rescale cross moments, helps keeps numbers real */
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			aa[i] [j] /= ( double )natoms;

	/* determinat */
	det = aa[0] [2] * ( aa[1] [0] * aa[2] [1] - aa[1] [1] * aa[2] [0] ) - aa[1] [2]
																				 * ( aa[0] [0] * aa[2] [1] - aa[0] [1] * aa[2] [0] ) + aa[2] [2] * ( aa[0] [0] * aa[1] [1] - aa[0] [1] * aa[1] [0] );

	if ( fabs( det ) <= 1.0E-24 )
	{
		return 20000;

	}
	handedness = ( det < 0 ) ? -1.0 : 1.0;

	//  multiply cross moments by itself
	for ( i = 0; i < 3; i++ )
	{
		for ( j = i; j < 3; j++ )
		{
			m[j] [i] = m[i] [j] = // well it is symmetric afterall
					aa[0] [i] * aa[0] [j] + aa[1] [i] * aa[1] [j] + aa[2] [i] * aa[2] [j];
		}
	}



	double ev[3];
	double a, b, c, r, q, sq, theta;
	double xx, yy, zz, xy, xz, yz;

	xx = m[0] [0];
	yy = m[1] [1];
	zz = m[2] [2];
	xy = m[0] [1];
	xz = m[0] [2];
	yz = m[1] [2];

	// coefficients of characterisitic polynomial
	a = xx + yy + zz;
	b = -xx * zz - xx * yy - yy * zz + xy * xy + xz * xz + yz * yz;
	c = xx * yy * zz - xz * xz * yy - xy * xy * zz - yz * yz * xx + 2 * xy * xz * yz;

	/* eigenvals are the roots of the characteristic polymonial  0 = c+b*e +a*e^2  - e^3 */
	/* solving for the three roots now: */
	/* dont try to follow this in detail: its just a tricky */
	/* factorization of the formulas for cubic equation roots. */

	/* real roots of cubic (sort of butt ugly) */
	b = -b;
	q = ( a * a - 3.0 * b ) / 9.0;
	r = ( 2.0 * a * a * a - 9.0 * a * b + 27.0 * c ) / 54.0;
	if ( q * q * q < r * r )
	{
		return 1;
		/** catch this later on * */
	}
	sq = sqrt( q );
	theta = acos( r / ( sq * sq * sq ) );
	sq *= 2.0;
	a /= -3.0;

	ev[0] = sq * cos( theta / 3.0 ) - a;
	ev[1] = sq * cos( ( theta + 2.0 * PDB_PI ) / 3.0 ) - a;
	ev[2] = sq * cos( ( theta - 2.0 * PDB_PI ) / 3.0 ) - a;

	// printf( "%f %f %f\n", ev[0], ev[1], ev[2] );

	if ( handedness < 0 )
	{

		// reorder eigen values  so that ev(3) is the smallest eigenvalue

		if ( ev[1] > ev[2] )
		{
			if ( ev[2] > ev[0] )
			{
				tmp = ev[2];
				ev[2] = ev[0];
				ev[0] = tmp;
			}
		}
		else
		{
			if ( ev[1] > ev[0] )
			{
				tmp = ev[2];
				ev[2] = ev[0];
				ev[0] = tmp;
			}
			else
			{
				tmp = ev[2];
				ev[2] = ev[1];
				ev[1] = tmp;
			}
		}
	}
	//                 // ev(3) is now the smallest eigen value.  the other two are not
	//                 //  sorted.

	//            // now we must catch the special case of the rotation with inversion.
	//            // we cannot allow inversion rotations.
	//            // fortunatley, and curiously, the optimal non-inverted rotation
	//            // matrix
	//            // will have the similar eigen values.


	double rms_ctx, rms_sum;
	rms_ctx = sqrt( fabs( ev[0] ) ) + sqrt( fabs( ev[1] ) ) + handedness * sqrt( fabs( ev[2] ) );


	end = true; end2 = true;
	initAll();
	mol2->initAll();
	natoms = 0;
	rms_sum = 0.0;
	while ( end != false )
	{
		/* get atoms */
		at = getCurrentAtom();
		at2 = mol2->getCurrentAtom();

		/* get positions */
		at->getPosition(pos1);
		at2->getPosition(pos2);

		for ( i = 0; i < 3; i++ )
		{
			tmp = ( pos2[i] - center2[i] );
			tmp1 = ( pos1[i] - center1[i] );
			rms_sum += (tmp * tmp + tmp1 * tmp1);
		}

		end = nextAtom(); end2 = mol2->nextAtom();
		natoms++;
	}

	rms_sum /= ( double )natoms;

	/* and combine the outer and cross terms into the final calculation. */
	/* (the abs() just saves us a headache when the roundoff error accidantally makes the sum negative) */

	//   printf( "%f %f %f\n", rms_sum, rms_ctx * 2, fabs( rms_sum - rms_ctx * 2.0) );
	err = sqrt( ( tmp = rms_sum - rms_ctx * 2., fabs( tmp ) ) );
	return float( err );

}



bool closeResidues(Fragment * res1,Fragment *res2, float square_distance,float square_limit)
{
	bool contA,contA2;
	Atom *a,*a2;
	Tcoor pos_atom1,pos_atom2;
	float dist;
	bool first=true;
	contA = true;
	res1->initAll();
	while ( contA != false )
	{
		a = ( Atom * ) res1->getCurrent();
		a->getPosition(pos_atom1);

		res2->initAll();
		contA2=true;
		while( contA2 != false )
		{
			a2 = (Atom*) res2->getCurrent();
			a2->getPosition(pos_atom2);
			dist=(pos_atom2[0]-pos_atom1[0])*(pos_atom2[0]-pos_atom1[0])+
					(pos_atom2[1]-pos_atom1[1])*(pos_atom2[1]-pos_atom1[1])+
					(pos_atom2[2]-pos_atom1[2])*(pos_atom2[2]-pos_atom1[2]);

			if(dist<square_distance)
				return true;

			if(first)
			{
				if (dist > square_limit)
					return false;
			}
			first=false;

			contA2=res2->next();
		}

		contA=res1->next();
	}

	return false;

}

bool sameResidues(Fragment * res1,Fragment *res2, float square_distance, float *distance, bool name)
{
	bool contA;
	Atom *a;
	Tcoor pos_atom1,pos_atom2;
	float dist;
	char *name1,*name2;
	char *atom_name1, *atom_name2;
	bool foundCA1=false, foundCA2=false;
	char main_atom[10];

	if(res1->getClass()!=res2->getClass())
		return false;
	else
		switch (res1->getClass())
		{
		case pdb_residue:
			strcpy(main_atom," CA ");
			break;
		case pdb_nucleotide:
			strcpy(main_atom," P  ");
			break;
		default :
			return false;
			break;
		}

	name1=res1->getName();
	name2=res2->getName();
	if(strcmp(name1,name2)==0 || !name)
	{
		contA = true;
		res1->initAll();
		while ( contA != false )
		{
			a = ( Atom * ) res1->getCurrent();
			a->getPosition(pos_atom1);
			atom_name1=a->getName();
			if(strcmp(atom_name1,main_atom)==0)
			{
				contA=false;
				foundCA1=true;
			}
			else
				contA=res1->next();
		}

		contA = true;
		res2->initAll();
		while ( contA != false )
		{
			a = ( Atom * ) res2->getCurrent();
			a->getPosition(pos_atom2);
			atom_name2=a->getName();
			if(strcmp(atom_name2,main_atom)==0)
			{
				contA=false;
				foundCA2=true;
			}
			else
				contA=res2->next();
		}

		if(foundCA1 && foundCA2)
		{

			dist=(pos_atom2[0]-pos_atom1[0])*(pos_atom2[0]-pos_atom1[0])+
					(pos_atom2[1]-pos_atom1[1])*(pos_atom2[1]-pos_atom1[1])+
					(pos_atom2[2]-pos_atom1[2])*(pos_atom2[2]-pos_atom1[2]);

			if(dist<square_distance)
			{
				*distance=dist;
				return true;
			}
		}
		else
			if(!foundCA1)
			{
				*distance=-1;
				return false;
			}

	}
	*distance=0;
	return false;

}


int** Macromolecule::getInterface( Macromolecule * mol2, float distance, int *numContacts,
		char ***interfaceC, bool withHET )
{
	pdbIter *iterS1, *iterS2, *iter1, *iter2;
	Fragment * res1,*res2;
	Chain * seg, * seg2;

	int cont_contacts=0;
	int **interface;
	float square_distance=distance*distance;
	float square_limit=((30.0-distance)*(30.0-distance));
	bool closeEnough;

	interface=(int**)malloc(sizeof(int*)*2);
	(*interfaceC)=(char**)malloc(sizeof(char*)*2);

	interface[0]=(int*)malloc(0);
	interface[1]=(int*)malloc(0);
	(*interfaceC)[0]=(char*)malloc(0);
	(*interfaceC)[1]=(char*)malloc(0);


	//Iterador para recorrer segmentos
	iterS1 = new pdbIter( this );
	iterS2=  new pdbIter( mol2);

	char ichain1, ichain2; // current chain

	for ( iterS1->pos_chain = 0; !iterS1->gend_chain(); iterS1->next_chain() )
		{
			seg = (Chain *) iterS1->get_chain();
			ichain1 = (seg->getName())[0]; // Get current chain id (1-letter)
			iter1 = new pdbIter( seg );

			//Bucle para recorrer residuos
			for ( iter1->pos_fragment = 0; !iter1->gend_fragment(); iter1->next_fragment() )
			{
				res1 = ( Fragment * ) iter1->get_fragment();

				for ( iterS2->pos_chain = 0; !iterS2->gend_chain(); iterS2->next_chain() )
						{

							seg2 = ( Chain * ) iterS2->get_chain();
							ichain2 = (seg2->getName())[0]; // Get current chain id (1-letter)
							iter2 = new pdbIter( seg2 );


							//Bucle para recorrer residuos
							for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
							{
								res2 = ( Fragment * ) iter2->get_fragment();

								closeEnough=closeResidues(res1,res2, square_distance,square_limit);

								if(closeEnough)
								{
									interface[0]=(int*)realloc(interface[0],(cont_contacts+1)*sizeof(int));
									interface[1]=(int*)realloc(interface[1],(cont_contacts+1)*sizeof(int));
									(*interfaceC)[0]=(char*)realloc((*interfaceC)[0],(cont_contacts+1)*sizeof(char));
									(*interfaceC)[1]=(char*)realloc((*interfaceC)[1],(cont_contacts+1)*sizeof(char));

									interface[0][cont_contacts]=res1->getIdNumber();
									interface[1][cont_contacts]=res2->getIdNumber();
									(*interfaceC)[0][cont_contacts]=ichain1;
									(*interfaceC)[1][cont_contacts]=ichain2;
									cont_contacts++;
								}

							}
							delete iter2;

						}

			}
			delete iter1;

		}
	delete iterS1;



//		iter2=new pdbIter(mol2,withHET);
//    iter1->pos_fragment=0;
//	do
//	{
//		res1=iter1->get_fragment();
//
//
//		iter2->pos_fragment=0;
//		closeEnough=false;
//		cont2=0;
//		do{
//			res2=iter2->get_fragment();
//
//			closeEnough=closeResidues(res1,res2, square_distance,square_limit);
//
//			if(closeEnough)
//			{
//				interface[0]=(int*)realloc(interface[0],(cont_contacts+1)*sizeof(int));
//				interface[1]=(int*)realloc(interface[1],(cont_contacts+1)*sizeof(int));
//
//				interface[0][cont_contacts]=res1->getIdNumber();
//				interface[1][cont_contacts]=res2->getIdNumber();
//				cont_contacts++;
//			}
//
//			cont2++;
//			iter2->next_fragment();
//		}while(!iter2->gend_fragment());
//
//		cont++;
//		iter1->next_fragment();
//	}while(!iter1->gend_fragment());
//
//	delete(iter1);
//	delete(iter2);
	*numContacts=cont_contacts;
	return interface;
}

int* Macromolecule::relative(Macromolecule *mol2, int *list,float distance, int total)
{
	pdbIter *iter1, *iter2,*iter_aux1,*iter_aux2;
	Fragment * res1,*res2,*resAux1,*resAux2;
	int cont=0, cont2=0;
	float Distance,dist;
	float square_distance=distance*distance;
	bool closeEnough;
	int *relatives;
	int i,j,step;
	char *name1,*name2;
	bool wrong;

	iter1=new pdbIter(this,false);
	iter2=new pdbIter(mol2,false);

	relatives=(int*)malloc(sizeof(int)*total);

	for(i=0;i<total;i++)
	{
		cont=list[i];
		iter1->pos_fragment=list[i];
		res1=iter1->get_fragment();
		relatives[i]=-1;

		Distance=square_distance;
		iter2->pos_fragment=0;
		closeEnough=false;
		cont2=0;
		do{
			res2=iter2->get_fragment();
			closeEnough=sameResidues(res1,res2, Distance,&dist,true);
			if(closeEnough)
			{
				Distance=dist;
				relatives[i]=cont2;
			}
			else
				if(dist==-1)
				{
					relatives[i]=-2;
				}
			cont2++;
			iter2->next_fragment();
		}while(!iter2->gend_fragment());
	}

	iter_aux1=new pdbIter(this,false);
	iter_aux2=new pdbIter(mol2,false);

	for(step=1;step<=3;step++)
	{
		square_distance=distance+(3*step);
		square_distance*=square_distance;
		for(i=0;i<total;i++)
		{
			iter1->pos_fragment=list[i];
			res1=iter1->get_fragment();

			if (relatives[i]==-1)
			{
				Distance=square_distance;
				iter2->pos_fragment=0;
				do{
					res2=iter2->get_fragment();
					closeEnough=sameResidues(res1,res2, Distance,&dist,false);
					if(closeEnough)
					{
						wrong=false;
						iter_aux1->pos_fragment=list[i];
						iter_aux2->pos_fragment=iter2->pos_fragment;
						for(j=-2;j<0 && !wrong;j++)
						{
							iter_aux1->back_fragment();
							iter_aux2->back_fragment();
							resAux1=iter_aux1->get_fragment();
							resAux2=iter_aux2->get_fragment();

							if(resAux1==NULL || resAux2==NULL)
								j=0;
							else
							{
								name1=resAux1->getName();
								name2=resAux2->getName();
								if(strcmp(name1,name2)!=0)
									wrong=true;
							}
						}
						iter_aux1->pos_fragment=list[i];
						iter_aux2->pos_fragment=iter2->pos_fragment;
						for(j=2;j>0 && !wrong;j--)
						{
							iter_aux1->next_fragment();
							iter_aux2->next_fragment();
							resAux1=iter_aux1->get_fragment();
							resAux2=iter_aux2->get_fragment();

							if(resAux1==NULL || resAux2==NULL)
								j=0;
							else
							{
								name1=resAux1->getName();
								name2=resAux2->getName();
								if(strcmp(name1,name2)!=0)
									wrong=true;
							}

						}

						if(!wrong)
						{
							Distance=dist;
							relatives[i]=iter2->pos_fragment;
						}
					}
					iter2->next_fragment();
				}while(!iter2->gend_fragment());
			}

		}
	}
	delete(iter1);
	delete(iter2);
	delete(iter_aux1);
	delete(iter_aux2);
	return relatives;

}


char * Macromolecule::getResInfo()
{
	pdbIter *iter1;
	int max;
	char *lista,*name;
	Fragment *res;
	iter1=new pdbIter(this,false);
	max=iter1->num_fragment();
	lista=(char*)malloc(sizeof(char)*(3*max+1));
	lista[0]='\0';

	iter1->pos_fragment=0;
	do
	{
		res=iter1->get_fragment();
		name=res->getName();
		strcat(lista,name);
		iter1->next_fragment();
	}while(!iter1->gend_fragment());

	delete iter1;
	return lista;
}


// Save coordinates of N, CA, C, 0 in a vector
void Macromolecule::coordBackbone(vector < vector<float> > & backbone)
{
	pdbIter *iter1;
	Tcoor pos1;

	iter1=new pdbIter(this);

	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		(iter1->get_atom())->getPosition(pos1);

		if ( strcmp((iter1->get_atom())->getName(), " N  " ) == 0 ||
				strcmp((iter1->get_atom())->getName(), " CA " ) == 0 ||
				strcmp((iter1->get_atom())->getName(), " C  " ) == 0) {

			vector < float > a( 3 );
			a[0]=pos1[0];
			a[1]=pos1[1];
			a[2]=pos1[2];
			backbone.push_back(a);
		}
	}

	delete iter1;

}

// Save coordinates of CA atoms
void Macromolecule::loadCAs(vector < vector<float> > & CAs)
{
	pdbIter *iter1;
	Tcoor pos1;

	iter1=new pdbIter(this);

	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		(iter1->get_atom())->getPosition(pos1);

		if ( strcmp((iter1->get_atom())->getName(), " CA " ) == 0 ) {

			vector < float > a( 3 );
			a[0]=pos1[0];
			a[1]=pos1[1];
			a[2]=pos1[2];
			CAs.push_back(a);
		}
	}

	delete iter1;

}


// Load the coordinates N, CA and C of the first aa int fname PDB
void Macromolecule::loadNCACFirstAA( vector <vector <float> > & CAs ) {

	pdbIter *iter1;
	Tcoor pos1;

	iter1=new pdbIter(this);

	for(iter1->pos_atom=0;!iter1->gend_atom();iter1->next_atom())
	{
		(iter1->get_atom())->getPosition(pos1);

		if ( strcmp((iter1->get_atom())->getName(), " N  " ) == 0 ||
				strcmp((iter1->get_atom())->getName(), " CA " ) == 0 ||
				strcmp((iter1->get_atom())->getName(), " C  " ) == 0 ) {

			vector < float > a( 3 );
			a[0]=pos1[0];
			a[1]=pos1[1];
			a[2]=pos1[2];
			CAs.push_back(a);

			if ( strcmp((iter1->get_atom())->getName(), " C  " ) == 0 )
				break;

		}
	}

	delete iter1;


}

char* Macromolecule::get_sequence()
{
	int n_fragments,i;
	char *seq,*seq2;

	Conditions *conds = new Conditions();
	Condition *cond = new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
	// CA selection
	cond->add(CAstr);
	conds->add(cond);

	Macromolecule *mol_ca = this->select_cpy(conds);
	seq2 = mol_ca->getResInfo();
	n_fragments = mol_ca->get_num_fragments();


	seq=(char *)malloc(sizeof(char)*(n_fragments+1));
	seq[0] = '\0';

	for(int n=0; n<n_fragments; n++)
	{
		for(i=0;i<N_AMINO;i++)
		{
			if(strncmp( seq2 + 3*n, AA[i].aa_name3,3)==0)
			{
				*(seq+n) = AA[i].aa_name1;
				break;
			}
		}
	}
	free(seq2);
	*(seq+n_fragments)='\0';
	delete mol_ca;
	delete conds;
	return seq;
}

bool Macromolecule::check_asa()
{
	int cont_atom=0;
	float cont_asa=0;
	Atom *a;

	initAll();
	do
	{
		cont_atom++;
		a = getCurrentAtom();
		cont_asa+=a->getPdbocc();
	}
	while ( nextAtom() != false );



	if(cont_asa>(float)(cont_atom+10))
		return true;
	else
		return false;

}


char* Macromolecule::secondary_structure(bool output_char)
{
	Macromolecule *prot_ca;
	Conditions *conds;
	Condition *cond;
	pdbIter *iter;
	Atom *atom[5];
	char *ss;
	float atom_coor[5][3];

	int i=0, j, i_d;
	double f[3][3];
	float  min1, min2;
	double f_beta[3];
	double a_helix[3];
	double a_beta;
	char helix0,helix1,helix2,beta,coil;

	helix0='H';
	helix1='G';
	helix2='I';
	beta='E';
	coil='C';

	//Select only Alpha Carbons
	conds= new Conditions();
	cond= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
	cond->add(CAstr);
	conds->add(cond);
	prot_ca = this->select_cpy(conds);
	//delete(conds);
	//Atom iterator
	iter= new pdbIter(prot_ca);

	//____________________allocate memory_____________
	ss = (char *) malloc( sizeof(char ) * ( prot_ca->get_num_atoms() + 1 ) ); // Mon 2016/09/20: Minor bug revomed (+1 required)
	ss[prot_ca->get_num_atoms()]='\0'; // Termination character is mandatory

	double **d;
	d = (double **) malloc( sizeof(double *) * prot_ca->get_num_atoms() );
	for (int p=0; p<prot_ca->get_num_atoms(); p++)
	{
		d[p] = (double *) malloc( sizeof(double) * 4 );
		d[p][0]=0.0;
		d[p][1]=0.0;
		d[p][2]=0.0;
		d[p][3]=0.0;
	}

	//__________________For all residues except the last 4__________________
	for(iter->pos_atom=0,i_d=0; iter->pos_atom < (iter->num_atom()-4) ; iter->next_atom(),i_d++)
	{
		atom[0]=iter->get_atom();
		atom[0]->getPosition(atom_coor[0]);

		for(i=1;i<5;i++)
		{
			iter->pos_atom++;
			atom[i]=iter->get_atom();
			atom[i]->getPosition(atom_coor[i]);
		}

		//__________________Ecuation for the helix__________________

		/* Ecuation (8) */

		/*X coordenate*/
		f[0][0] = f_a(atom_coor[2][0], atom_coor[1][0], atom_coor[0][0], atom_coor[3][0], atom_coor[4][0]);
		f[0][1] = f_pi(atom_coor[2][0], atom_coor[1][0], atom_coor[0][0], atom_coor[3][0], atom_coor[4][0]);
		f[0][2] = f_310(atom_coor[2][0], atom_coor[1][0], atom_coor[0][0], atom_coor[3][0], atom_coor[4][0]);

		/*Y coordenate*/
		f[1][0] = f_a(atom_coor[2][1], atom_coor[1][1], atom_coor[0][1], atom_coor[3][1], atom_coor[4][1]);
		f[1][1] = f_pi(atom_coor[2][1], atom_coor[1][1], atom_coor[0][1], atom_coor[3][1], atom_coor[4][1]);
		f[1][2] = f_310(atom_coor[2][1], atom_coor[1][1], atom_coor[0][1], atom_coor[3][1], atom_coor[4][1]);

		/*Z coordenate*/
		f[2][0] = f_a(atom_coor[2][2], atom_coor[1][2], atom_coor[0][2], atom_coor[3][2], atom_coor[4][2]);
		f[2][1] = f_pi(atom_coor[2][2], atom_coor[1][2], atom_coor[0][2], atom_coor[3][2], atom_coor[4][2]);
		f[2][2] = f_310(atom_coor[2][2], atom_coor[1][2], atom_coor[0][2], atom_coor[3][2], atom_coor[4][2]);

		/*Ecuation (13)*/
		d[i_d][0] = sqrt( ((f[0][0] + 1)*(f[0][0] + 1))+((f[1][0] + 1)*(f[1][0] + 1))+((f[2][0] + 1)*(f[2][0] + 1)) );
		d[i_d][1] = sqrt( ((f[0][1] + 1)*(f[0][1] + 1))+((f[1][1] + 1)*(f[1][1] + 1))+((f[2][1] + 1)*(f[2][1] + 1)) );
		d[i_d][2] = sqrt( ((f[0][2] + 1)*(f[0][2] + 1))+((f[1][2] + 1)*(f[1][2] + 1))+((f[2][2] + 1)*(f[2][2] + 1)) );

		//__________________Ecuation for the β-strand__________________

		/*Ecuation (11)*/

		/*X coordenate*/
		f_beta[0] = f_b(atom_coor[3][0], atom_coor[1][0], atom_coor[2][0], atom_coor[0][0]);

		/*Y coordenate*/
		f_beta[1] = f_b(atom_coor[3][1], atom_coor[1][1], atom_coor[2][1], atom_coor[0][1]);

		/*Z coordenate*/
		f_beta[2] = f_b(atom_coor[3][2], atom_coor[1][2], atom_coor[2][2], atom_coor[0][2]);


		/*Ecuation (14)*/
		d[i_d][3] = sqrt( ((f_beta[0] - 1)*(f_beta[0] - 1))+((f_beta[1] - 1)*(f_beta[1] - 1))+((f_beta[2] - 1)*(f_beta[2] - 1)) );

		iter->pos_atom-=4;

	}
	//_____________________For the last 4 residues____________________
	// MON: WARNING WITH "-4" 's ...
	for(iter->pos_atom=iter->num_atom()-4,i_d=iter->num_atom()-4; iter->pos_atom < (iter->num_atom()) ; iter->next_atom(),i_d++)
	{
		atom[4]=iter->get_atom();
		atom[4]->getPosition(atom_coor[0]);

		int atom_coor_pos=0;

		for(i=3;i>=0;i--)
		{
			iter->pos_atom--;
			atom_coor_pos++;
			atom[i]=iter->get_atom();
			atom[i]->getPosition(atom_coor[atom_coor_pos]);
		}

		//__________________Ecuation for the helix__________________

		/* Ecuation (8) */

		/*X coordenate*/
		f[0][0] = f_a(atom_coor[2][0], atom_coor[1][0], atom_coor[0][0], atom_coor[3][0], atom_coor[4][0]);
		f[0][1] = f_pi(atom_coor[2][0], atom_coor[1][0], atom_coor[0][0], atom_coor[3][0], atom_coor[4][0]);
		f[0][2] = f_310(atom_coor[2][0], atom_coor[1][0], atom_coor[0][0], atom_coor[3][0], atom_coor[4][0]);

		/*Y coordenate*/
		f[1][0] = f_a(atom_coor[2][1], atom_coor[1][1], atom_coor[0][1], atom_coor[3][1], atom_coor[4][1]);
		f[1][1] = f_pi(atom_coor[2][1], atom_coor[1][1], atom_coor[0][1], atom_coor[3][1], atom_coor[4][1]);
		f[1][2] = f_310(atom_coor[2][1], atom_coor[1][1], atom_coor[0][1], atom_coor[3][1], atom_coor[4][1]);

		/*Z coordenate*/
		f[2][0] = f_a(atom_coor[2][2], atom_coor[1][2], atom_coor[0][2], atom_coor[3][2], atom_coor[4][2]);
		f[2][1] = f_pi(atom_coor[2][2], atom_coor[1][2], atom_coor[0][2], atom_coor[3][2], atom_coor[4][2]);
		f[2][2] = f_310(atom_coor[2][2], atom_coor[1][2], atom_coor[0][2], atom_coor[3][2], atom_coor[4][2]);

		/*Ecuation (13)*/
		d[i_d][0] = sqrt( ((f[0][0] + 1)*(f[0][0] + 1))+((f[1][0] + 1)*(f[1][0] + 1))+((f[2][0] + 1)*(f[2][0] + 1)) );
		d[i_d][1] = sqrt( ((f[0][1] + 1)*(f[0][1] + 1))+((f[1][1] + 1)*(f[1][1] + 1))+((f[2][1] + 1)*(f[2][1] + 1)) );
		d[i_d][2] = sqrt( ((f[0][2] + 1)*(f[0][2] + 1))+((f[1][2] + 1)*(f[1][2] + 1))+((f[2][2] + 1)*(f[2][2] + 1)) );

		//__________________Ecuation for the β-strand__________________

		/*Ecuation(11)*/

		/*X coordenate*/
		f_beta[0] = f_b(atom_coor[3][0], atom_coor[1][0], atom_coor[2][0], atom_coor[0][0]);

		/*Y coordenate*/
		f_beta[1] = f_b(atom_coor[3][1], atom_coor[1][1], atom_coor[2][1], atom_coor[0][1]);

		/*Z coordenate*/
		f_beta[2] = f_b(atom_coor[3][2], atom_coor[1][2], atom_coor[2][2], atom_coor[0][2]);

		/*Ecuation (14)*/
		d[i_d][3] = sqrt( ((f_beta[0] - 1)*(f_beta[0] - 1))+((f_beta[1] - 1)*(f_beta[1] - 1))+((f_beta[2] - 1)*(f_beta[2] - 1)) );
		iter->pos_atom+=4;
	}
	//__________________Search minimum(first 4 residues)__________________

	for(i=0;i<4;i++)
	{
		/*helix*/
		aux_min_w(d,0,i,i,0,&min1,&min2);
		a_helix[0] = (min1 + min2)/2.0;

		aux_min_w(d,1,i,i,0,&min1,&min2);
		a_helix[1] = (min1 + min2)/2.0;

		aux_min_w(d,2,i,i,0,&min1,&min2);
		a_helix[2] = (min1 + min2)/2.0;
		/*β-strand*/
		aux_min_w(d,3,i,i,0,&min1,&min2);
		a_beta = (min1 + min2)/2.0;

		/* Print results */
		if (a_helix[0]<0.30)
			ss[i]=helix0;
		else if (a_helix[1]<0.30)
			ss[i]=helix0;
		else if (a_helix[2]<0.30)
			ss[i]=helix0;
		else if (a_beta<0.75)
			ss[i]=beta;
		else
			ss[i]=coil;

	}
	//____________________Search minimum_____________________

	for(i=4;i<iter->num_atom()-4;i++)
	{
		/*helix*/
		aux_min_w(d,0,i,i,i-4,&min1,&min2);
		a_helix[0] = (min1 + min2)/2.0;

		aux_min_w(d,1,i,i,i-4,&min1,&min2);
		a_helix[1] = (min1 + min2)/2.0;

		aux_min_w(d,2,i,i,i-4,&min1,&min2);
		a_helix[2] = (min1 + min2)/2.0;

		/*β-strand*/
		aux_min_w(d,3,i,i,i-3,&min1,&min2);
		a_beta = (min1 + min2)/2.0;

		/* Print results */
		if (a_helix[0]<0.30)
			ss[i]=helix0;
		else if (a_helix[1]<0.30)
			ss[i]=helix0;
		else if (a_helix[2]<0.30)
			ss[i]=helix0;
		else if (a_beta<0.75)
			ss[i]=beta;
		else
			ss[i]=coil;
	}
	//__________________Search minimum(last 4 residues)__________________

	for(i=iter->num_atom()-4;i<iter->num_atom();i++)
	{
		/*helix*/
		aux_min_w(d,0,i,iter->num_atom()-1,i,&min1,&min2);
		a_helix[0] = (min1 + min2)/2.0;

		aux_min_w(d,1,i,iter->num_atom()-1,i,&min1,&min2);
		a_helix[1] = (min1 + min2)/2.0;

		aux_min_w(d,2,i,iter->num_atom()-1,i,&min1,&min2);
		a_helix[2] = (min1 + min2)/2.0;

		/*β-strand*/
		aux_min_w(d,3,i,iter->num_atom()-1,i,&min1,&min2);
		a_beta = (min1 + min2)/2.0;

		/* Print results */
		if (a_helix[0]<0.30)
			ss[i]=helix0;
		else if (a_helix[1]<0.30)
			ss[i]=helix0;
		else if (a_helix[2]<0.30)
			ss[i]=helix0;
		else if (a_beta<0.75)
			ss[i]=beta;
		else
			ss[i]=coil;
	}
	char lon;
	int cl;
	int clold=1;
	lon=ss[4];
	cl=1;
	for(i=4;i<iter->num_atom()-4;i++)
	{
		if (ss[i]==lon) cl++;
		else
		{
			if ((lon==beta)&&(ss[i]==helix0))
			{
				ss[i-1]=coil;
				if (ss[i-2]==beta)
					ss[i-2]=coil;
			}

			if ((lon==helix0)&&(ss[i]==beta))
			{
				ss[i]=coil;
				if (ss[i+1]==beta)
					ss[i+1]=coil;
			}
			lon=ss[i];
			clold=cl;
			cl=1;
		}
	}

	lon=ss[0];
	cl=1;
	for(i=1;i<iter->num_atom();i++)
	{
		if (ss[i]==lon) cl++;
		else
		{
			if ((cl<=1)&&(lon!=coil)&&(lon!=helix0)) {
				for(j=0;j<cl;j++)
					ss[i-j-1]=coil;
			}

			if ((cl<=3)&&(lon!=coil)&&(lon!=beta)) {
				for(j=0;j<cl;j++)
					ss[i-j-1]=coil;
			}
			lon=ss[i];
			cl=1;
		}
	}

	for(i=0,iter->pos_atom=0; iter->pos_atom < iter->num_atom() ; i++,iter->next_atom())
	{
		((Residue*)((iter->get_atom())->getFather()))->set_ss(ss[i]);
	}


	delete iter;
	for (int p=0; p<prot_ca->get_num_atoms(); p++)
	{
		free(d[p]);
	}
	free(d);

	//	prot_ca->erase_level(pdb_atom);
	delete prot_ca;
	delete conds;

	if( output_char==false)
	{
		free(ss);
		return NULL;
	}

	return ss;
}

char *Macromolecule::get_ss()
{
	Conditions *conds;
	Condition *cond;
	pdbIter *iter;
	Macromolecule *prot_ca;
	int i,num;
	char *cad;

	//Select only Alpha Carbons
	conds= new Conditions();
	cond= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
	cond->add(CAstr);
	conds->add(cond);
	prot_ca=this->select_cpy(conds);
	delete(conds);
	//Atom iterator
	iter= new pdbIter(prot_ca);

	num=prot_ca->get_num_atoms();
	cad=(char*)malloc(sizeof(char)*num+1);

	for(i=0,iter->pos_atom=0; iter->pos_atom < iter->num_atom() ; i++,iter->next_atom())
	{
		cad[i]=((Residue*)((iter->get_atom())->getFather()))->get_ss();
		//fprintf(stderr,"cad[%d]=%c\n",i,cad[i]);
	}
	cad[num]='\0';

	//	prot_ca->erase_level(pdb_atom);
	delete(prot_ca);
	delete(iter);
	delete(conds);

	return cad;
}

// It will make an inter-atomic distance profile between two macromolecules
// (It allocates memory! When *profile = NULL) (It also returns the RMSD)
// "mask" and "mask2" are the atom-wise masks to compute inter atomic distances (if != NULL, only CA atoms will be considered)
double Macromolecule::dist_profile( Macromolecule *mol2, double **profile, bool *mask1, bool *mask2 )
{
	bool debug = false;
	Atom *at1, *at2;
	Tcoor pos1, pos2;
	double suma = 0.0, tmp, dist;
	int a, natoms, totatoms=0, res_index=0, nres=0;
	double *p_profile = *profile;

	pdbIter *iter_res,*iter_atom,*iter2_res,*iter2_atom;
	Residue *res, *res2;

	iter_res = new pdbIter(this,false);
	//	nres = this->get_num_fragments();
	nres = iter_res->num_fragment();

	if(p_profile == NULL) // If NULL, allocates memory
	{
		p_profile = (double *) malloc( sizeof(double) * nres);
		if(!p_profile)
		{
			printf("Msg(dist_profile): Unable to allocate memory!\n");
			exit(1);
		}
	}
	else
		p_profile = *profile;

	// Initializing profile array
	for( int i=0; i<nres; i++)
		p_profile[i] = 0.0;

	pdbIter *iter1,*iter2;
	int atom_index = 0;
	if(mask1 != NULL)
	{
		//		iter_res = new pdbIter(this,false);
		iter2 = new pdbIter(mol2,false);
		iter2->pos_atom=0;

		for(iter_res->pos_fragment=0; !iter_res->gend_fragment(); iter_res->next_fragment())
		{
			iter1 = new pdbIter(iter_res->get_fragment(),false);
			for(iter1->pos_atom=0; !iter1->gend_atom(); iter1->next_atom())
			{
				if(mask1[atom_index]) // Mol 1 matching atom
				{
					while(!iter2->gend_atom() && !mask2[iter2->pos_atom]) // this places the index into the corresponding Mol 2 atom
						iter2->next_atom();

					if(!iter2->gend_atom())
					{
						// get atoms
						at1 = iter1->get_atom();
						at2 = iter2->get_atom();

						// get positions
						at1->getPosition(pos1);
						at2->getPosition(pos2);

						tmp=0.0;
						for ( a = 0; a < 3; a++ )
							tmp += (pos1[a] - pos2[a])*(pos1[a] - pos2[a]);
						p_profile[res_index] = sqrt(tmp); // Distance

						if(debug)
							printf("Getting atoms %d (%s) and %d (%s): dist= %f  index= %d  iter2_pos_atom= %d\n",iter1->pos_atom,at1->getName(),iter2->pos_atom,at2->getName(),p_profile[res_index],res_index,iter2->pos_atom);

						suma += tmp;
						iter2->next_atom(); // next atom Mol 2
					}
				}
				atom_index++; // Mol 1 atom index
			}
			delete(iter1);
			res_index++; // next residue in Mol 1
		}
		delete(iter2);
		if(res_index != nres)
		{
			printf("Sorry, number of aligned atoms mismatch. Aligned CA/P atoms (%d) and Residues (%d) should be equal.\n",res_index,nres);
			exit(2);
		}
	}
	else
	{
		iter_res = new pdbIter( this );
		iter2_res = new pdbIter( mol2 );
		for( iter_res->pos_fragment=0, iter2_res->pos_fragment=0;
				!iter_res->gend_fragment() && !iter2_res->gend_fragment();
				iter_res->next_fragment(), iter2_res->next_fragment() )
		{
			res_index = iter_res->pos_fragment; // shorter...
			res = (Residue *) iter_res->get_fragment();
			res2 = (Residue *) iter2_res->get_fragment();
			iter_atom = new pdbIter( res );
			iter2_atom = new pdbIter( res2 );
			natoms=0;
			for( iter_atom->pos_atom=0, iter2_atom->pos_atom=0;
					!iter_atom->gend_atom() && !iter2_atom->gend_atom();
					iter_atom->next_atom(), iter2_atom->next_atom() )
			{
				// get atoms
				at1 = iter_atom->get_atom();
				at2 = iter2_atom->get_atom();

				// get positions
				at1->getPosition(pos1);
				at2->getPosition(pos2);

				tmp=0;
				for ( a = 0; a < 3; a++ )
					tmp += (pos1[a] - pos2[a])*(pos1[a] - pos2[a]);
				dist = sqrt(tmp); // Distance
				p_profile[res_index] += dist;
				suma += tmp;
				natoms++; // counting atoms per residue
			}
			delete iter_atom;
			delete iter2_atom;
			totatoms += natoms; // counts all atoms
			p_profile[res_index] /= natoms; // Distance average between residue atom pairs.
			// Maximum distance??? Would it be better???

			//		printf("profile %d: %f\n",index,p_profile[index]);
		}
		delete iter_res;
		delete iter2_res;
	}
	*profile = p_profile; // Outputs profile...

	//printf("natoms %d %f\n", natoms, suma);

	if ( suma != 0.0 ) return sqrt( ( double )suma / ( double )totatoms );
	else
		return ( ( double )suma );
}

//***********************************final**************************************


void aux_min_w (double **d, int p, int i, int start, int end,  float *min1, float *min2)
{
	int  minpos=-9999;


	*min1=10000;
	*min2=1000;

	for(int j=start;j>=end;j--)
	{
		if (*min1 > d[j][p])
		{
			*min1=d[j][p];
			minpos=j;
		}
	}
	for(int j=start;j>=end;j--)
	{
		if ( (*min2 > d[j][p])  && (j!=minpos) )
			*min2 = d[j][p];
	}

}



// --------------------------------------------------------------
// minRmsd's Auxiliar functions (needed by: minRmsd and minWRmsd)
// --------------------------------------------------------------
void normalize(double a[3])
{
	double  b;

	b = sqrt((double)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
	a[0] /= b;
	a[1] /= b;
	a[2] /= b;
}
double dot(double a[3], double b[3])
{
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}
void cross(double a[3], double b[3], double c[3])
{
	a[0] = b[1]*c[2] - b[2]*c[1];
	a[1] = b[2]*c[0] - b[0]*c[2];
	a[2] = b[0]*c[1] - b[1]*c[0];
}
/*
 * setup_rotation()
 *
 *      given two lists of x,y,z coordinates, constructs
 * the correlation R matrix and the E value needed to calculate the
 * least-squares rotation matrix.
 */
void setup_rotation(double **ref_xlist,
		double **mov_xlist,
		int n_list,
		double mov_com[3],
		double mov_to_ref[3],
		double R[3][3],
		double *E0,
		double *weights)
{
	int i, j, n;
	double ref_com[3];
	double sum_w=0.0; // weights summation

	/* calculate the centre of mass */
	for (i=0; i<3; i++)
	{
		mov_com[i] = 0.0;
		ref_com[i] = 0.0;
	}

	if(weights == NULL) // un-weighted coords.
	{
		for (n=0; n<n_list; n++) // atoms
			for (i=0; i<3; i++) // coords.
			{
				mov_com[i] += mov_xlist[n][i];
				ref_com[i] += ref_xlist[n][i];
			}
		sum_w = (double) n_list;
	}
	else // weighted coords.
		for (n=0; n<n_list; n++) // atoms
		{
			for (i=0; i<3; i++) // coords.
			{
				mov_com[i] += mov_xlist[n][i] * weights[n];
				ref_com[i] += ref_xlist[n][i] * weights[n];
			}
			sum_w += weights[n];
		}

	for (i=0; i<3; i++)
	{
		//	  mov_com[i] /= n_list;
		//	  ref_com[i] /= n_list;
		mov_com[i] /= sum_w;
		ref_com[i] /= sum_w;
		mov_to_ref[i] = ref_com[i] - mov_com[i];
	}

	/* shift mov_xlist and ref_xlist to centre of mass */
	for (n=0; n<n_list; n++)
		for (i=0; i<3; i++)
		{
			mov_xlist[n][i] -= mov_com[i];
			ref_xlist[n][i] -= ref_com[i];
		}

	/* initialize */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			R[i][j] = 0.0;
	*E0 = 0.0;

	// Un-weighted
	if(weights == NULL)
		for (n=0; n<n_list; n++)
		{
			// E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n)
			for (i=0; i<3; i++)
				*E0 +=  mov_xlist[n][i] * mov_xlist[n][i]
													   + ref_xlist[n][i] * ref_xlist[n][i];

			// correlation matrix R:
			// R[i,j) = sum(over n): y(n,i) * x(n,j)
			// where x(n) and y(n) are two vector sets
			for (i=0; i<3; i++)
				for (j=0; j<3; j++)
					R[i][j] += mov_xlist[n][i] * ref_xlist[n][j];
		}
	else // Weighted
		for (n=0; n<n_list; n++)
		{
			// Warning; I don't know if "E0" expression is OK:
			// E0 = 1/2 * sum(over n): y(n)*y(n)*w(n) + x(n)*x(n)*w(n)
			for (i=0; i<3; i++)
				*E0 +=  mov_xlist[n][i] * mov_xlist[n][i] * weights[n]
																	+ ref_xlist[n][i] * ref_xlist[n][i] * weights[n];

			// correlation matrix R:
			// R[i,j) = sum(over n): y(n,i) * x(n,j) * w(n)
			// where x(n) and y(n) are two vector sets
			for (i=0; i<3; i++)
				for (j=0; j<3; j++)
					R[i][j] += mov_xlist[n][i] * ref_xlist[n][j] * weights[n];
		}

	*E0 *= 0.5;
}

/*
 * jacobi3
 *
 *    computes eigenval and eigen_vec of a real 3x3
 * symmetric matrix. On output, elements of a that are above
 * the diagonal are destroyed. d[1..3] returns the
 * eigenval of a. v[1..3][1..3] is a matrix whose
 * columns contain, on output, the normalized eigen_vec of
 * a. n_rot returns the number of Jacobi rotations that were required.
 */
int jacobi3(double a[3][3], double d[3], double v[3][3], int* n_rot)
{
	int count, k, i, j;
	double tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];

	/*Initialize v to the identity matrix.*/
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
			v[i][j] = 0.0;
		v[i][i] = 1.0;
	}

	/* Initialize b and d to the diagonal of a */
	for (i=0; i<3; i++)
		b[i] = d[i] = a[i][i];

	/* z will accumulate terms */
	for (i=0; i<3; i++)
		z[i] = 0.0;

	*n_rot = 0;

	/* 50 tries */
	for (count=0; count<50; count++)
	{

		/* sum off-diagonal elements */
		sum = 0.0;
		for (i=0; i<2; i++)
		{
			for (j=i+1; j<3; j++)
				sum += fabs(a[i][j]);
		}

		/* if converged to machine underflow */
		if (sum == 0.0)
			return(1);

		/* on 1st three sweeps... */
		if (count < 3)
			tresh = sum * 0.2 / 9.0;
		else
			tresh = 0.0;

		for (i=0; i<2; i++)
		{
			for (j=i+1; j<3; j++)
			{
				g = 100.0 * fabs(a[i][j]);

				/*  after four sweeps, skip the rotation if
				 *   the off-diagonal element is small
				 */
				if ( count > 3  &&  fabs(d[i])+g == fabs(d[i])
						&&  fabs(d[j])+g == fabs(d[j]) )
				{
					a[i][j] = 0.0;
				}
				else if (fabs(a[i][j]) > tresh)
				{
					h = d[j] - d[i];

					if (fabs(h)+g == fabs(h))
					{
						t = a[i][j] / h;
					}
					else
					{
						theta = 0.5 * h / (a[i][j]);
						t = 1.0 / ( fabs(theta) +
								(double)sqrt(1.0 + theta*theta) );
						if (theta < 0.0)
							t = -t;
					}

					c = 1.0 / (double) sqrt(1 + t*t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[i][j];

					z[i] -= h;
					z[j] += h;
					d[i] -= h;
					d[j] += h;

					a[i][j] = 0.0;

					for (k=0; k<=i-1; k++)
						ROTATE(a, k, i, k, j)

						for (k=i+1; k<=j-1; k++)
							ROTATE(a, i, k, k, j)

							for (k=j+1; k<3; k++)
								ROTATE(a, i, k, j, k)

								for (k=0; k<3; k++)
									ROTATE(v, k, i, k, j)

									++(*n_rot);
				}
			}
		}

		for (i=0; i<3; i++)
		{
			b[i] += z[i];
			d[i] = b[i];
			z[i] = 0.0;
		}
	}

	//printf("Too many iterations in jacobi3\n");
	return (0);
}
/*
 * diagonalize_symmetric
 *
 *    Diagonalize a 3x3 matrix & sort eigenval by size
 */
int diagonalize_symmetric(double matrix[3][3],
		double eigen_vec[3][3],
		double eigenval[3])
{
	int n_rot, i, j, k;
	double vec[3][3];
	double val;

	if (!jacobi3(matrix, eigenval, vec, &n_rot))
	{
		//printf("convergence failed\n");
		return (0);
	}

	/* sort solutions by eigenval */
	for (i=0; i<3; i++)
	{
		k = i;
		val = eigenval[i];

		for (j=i+1; j<3; j++)
			if (eigenval[j] >= val)
			{
				k = j;
				val = eigenval[k];
			}

		if (k != i)
		{
			eigenval[k] = eigenval[i];
			eigenval[i] = val;
			for (j=0; j<3; j++)
			{
				val = vec[j][i];
				vec[j][i] = vec[j][k];
				vec[j][k] = val;
			}
		}
	}

	/* transpose such that first index refers to solution index */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			eigen_vec[i][j] = vec[j][i];

	return (1);
}
/*
 * calculate_rotation_matrix()
 *
 *   calculates the rotation matrix U and the
 * rmsd from the R matrix and E0:
 */
int calculate_rotation_matrix(double R[3][3],
		double U[3][3],
		double E0,
		double* residual)
{
	int i, j, k;
	double Rt[3][3], RtR[3][3];
	double left_eigenvec[3][3], right_eigenvec[3][3], eigenval[3];
	double v[3];
	double sigma;

	/* build Rt, transpose of R  */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			Rt[i][j] = R[j][i];

	/* make symmetric RtR = Rt X R */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
		{
			RtR[i][j] = 0.0;
			for (k = 0; k<3; k++)
				RtR[i][j] += Rt[k][i] * R[j][k];
		}

	if (!diagonalize_symmetric(RtR, right_eigenvec, eigenval))
		return(0);

	/* right_eigenvec's should be an orthogonal system but could be left
	 * or right-handed. Let's force into right-handed system.
	 */
	cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);

	/* From the Kabsch algorithm, the eigenvec's of RtR
	 * are identical to the right_eigenvec's of R.
	 * This means that left_eigenvec = R x right_eigenvec
	 */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			left_eigenvec[i][j] = dot(&right_eigenvec[i][0], &Rt[j][0]);

	for (i=0; i<3; i++)
		normalize(&left_eigenvec[i][0]);

	/*
	 * Force left_eigenvec[2] to be orthogonal to the other vectors.
	 * First check if the rotational matrices generated from the
	 * orthogonal eigenvectors are in a right-handed or left-handed
	 * co-ordinate system - given by sigma. Sigma is needed to
	 * resolve this ambiguity in calculating the RMSD.
	 */
	cross(v, &left_eigenvec[0][0], &left_eigenvec[1][0]);
	if (dot(v, &left_eigenvec[2][0]) < 0.0)
		sigma = -1.0;
	else
		sigma = 1.0;
	for (i=0; i<3; i++)
		left_eigenvec[2][i] = v[i];

	/* calc optimal rotation matrix U that minimises residual */
	for (i=0;i<3; i++)
		for (j=0; j<3; j++)
		{
			U[i][j] = 0.0;
			for (k=0; k<3; k++)
				U[i][j] += left_eigenvec[k][i] * right_eigenvec[k][j];
		}

	*residual = E0 - (double) sqrt(fabs(eigenval[0]))
                								 - (double) sqrt(fabs(eigenval[1]))
												 - sigma * (double) sqrt(fabs(eigenval[2]));

	return (1);
}
// ------------------------------------------------------------------
// END minRmsd's Auxiliar functions (needed by: minRmsd and minWRmsd)
// ------------------------------------------------------------------

