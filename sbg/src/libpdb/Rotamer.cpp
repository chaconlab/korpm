#include <stdio.h>

#include <fstream>
//#include "ResIni.h"
#include "Macromolecule.h"
#include "pdbIter.h"

using namespace std;

float dihedral_bk( float * a1, float * a2, float * a3, float * a4 )
{
	float xij, yij, zij, xkj, ykj, zkj, xkl, ykl, zkl, dxi, dyi, dzi, gxi, gyi, gzi, bi, bk, ct,
	z1, z2, ap, s;
	/* fprintf( stderr, "\n%f %f %f\n", a1[0],a1[1],a1[2]); fprintf( stderr, "%f %f %f\n", a2[0],a2[1],a2[2]);
  fprintf( stderr, "%f %f %f\n", a3[0],a3[1],a3[2]); fprintf( stderr, "%f %f %f\n", a4[0],a4[1],a4[2]); */


	/* Calculate the vectors C,B,C */
	xij = a1[0] - a2[0];
	yij = a1[1] - a2[1];
	zij = a1[2] - a2[2];

	xkj = a3[0] - a2[0];
	ykj = a3[1] - a2[1];
	zkj = a3[2] - a2[2];

	xkl = a3[0] - a4[0];
	ykl = a3[1] - a4[1];
	zkl = a3[2] - a4[2];

	//fprintf( stderr, "%f %f %f\n", xij,yij,zij);
	//fprintf( stderr, "%f %f %f\n", xkj,ykj,zkj);
	//fprintf( stderr, "%f %f %f\n", xkl,ykl,zkl);


	/* Calculate the normals to the two planes n1 and n2 this is given as the cross products: AB x BC --------- = n1 |AB x BC|
  BC x CD --------- = n2 |BC x CD| */
	dxi = yij * zkj - zij * ykj;
	/* Normal to plane 1 */
	dyi = zij * xkj - xij * zkj;
	dzi = xij * ykj - yij * xkj;
	gxi = zkj * ykl - ykj * zkl;
	/* Mormal to plane 2 */
	gyi = xkj * zkl - zkj * xkl;
	gzi = ykj * xkl - xkj * ykl;



	/* Calculate the length of the two normals */
	bi = dxi * dxi + dyi * dyi + dzi * dzi;
	bk = gxi * gxi + gyi * gyi + gzi * gzi;
	ct = dxi * gxi + dyi * gyi + dzi * gzi;


	bi = ( float )sqrt( ( double )bi );
	bk = ( float )sqrt( ( double )bk );

	z1 = 1. / bi;
	z2 = 1. / bk;

	ct = ct * z1 * z2;
	if ( ct > 1.0 ) ct = 1.0;
	if ( ct < ( -1.0 ) ) ct = -1.0;
	ap = acos( ct );
	//fprintf( stderr, "%f %f %f\n", ap,z1,z2);

	s = xkj * ( dzi * gyi - dyi * gzi ) + ykj * ( dxi * gzi - dzi * gxi ) + zkj * ( dyi * gxi - dxi * gyi );

	if ( s < 0.0 ) ap = -ap;

	ap = ( ap > 0.0 ) ? PDB_PI - ap : -( PDB_PI + ap );

	//fprintf( stderr, "d %f %f\n", ap, ap*180/M_PI);
	//getchar();

	return ( ap * 180.0 / PDB_PI );
}


void rotamer_from_chi( float * chi, int aan, int * rot )
{

	int MAX_CHI = 4;
	int nchis;

	for ( int i = 0; i < MAX_CHI; i++ )
	{
		rot[i] = 0;
	}

	nchis = AA[aan].nchi;


	//------------------------------------------------------------------------------
	// if not protein  rot = 0


	if ( aan == PRO )
	{ // chi1 of pro
		if ( chi[0] >= 0.0 )
		{
			rot[0] = 1; rot[1] = 1; rot[2] = 0; rot[3] = 0;
		}
		if ( chi[0] <= 0.0 )
		{
			rot[0] = 2; rot[1] = 1; rot[2] = 0; rot[3] = 0;
		}
	}

	else if ( aan == GLY || aan == ALA )
	{
		rot[0] = 0; rot[1] = 0; rot[2] = 0; rot[3] = 0;

	}

	else
	{
		//// default to 0 for chi's greater than nchi, etc
		for ( int i = 0; i < nchis; i++ )
		{
			rot[i] = 0;

			if ( i == 0 ) // chi 1
			{
				// chi1 of all residues except P,G,A,RNA
				if ( ( chi[i] >= 0.0 ) && ( chi[i] <= 120.0 ) )
					rot[i] = 1;
				if ( fabs( chi[i] ) >= 120.0 )
					rot[i] = 2;
				if ( ( chi[i] >= -120.0 ) && ( chi[i] <= 0.0 ) )
					rot[i] = 3;
			}

			else if ( i == 1 )
			{ // chi 2
				if ( aan == ARG )
				{
					if ( ( chi[i] >= 0.0 ) && ( chi[i] <= 120.0 ) )
						rot[i] = 1;
					if ( fabs( chi[i] ) >= 120.0 )
						rot[i] = 2;
					if ( ( chi[i] >= -120.0 ) && ( chi[i] <= 0.0 ) )
						rot[i] = 3;
				}

				else if ( aan == GLU || aan == HIS || aan == ILE || aan == LYS || aan == LEU || aan == MET || aan == GLN )
				{
					if ( ( chi[i] >= 0.0 ) && ( chi[i] <= 120.0 ) )
						rot[i] = 1;
					if ( fabs( chi[i] ) > 120.0 )
						rot[i] = 2;
					if ( ( chi[i] >= -120.0 ) && ( chi[i] <= 0.0 ) )
						rot[i] = 3;
				}

				else if ( aan == ASP )
				{
					if ( ( ( chi[i] >= 30.0 ) && ( chi[i] <= 90.0 ) ) || ( ( chi[i] <= -90.0 ) && ( chi[i] >= -150.0 ) ) )
						rot[i] = 1;
					if ( ( ( chi[i] >= -30.0 ) && ( chi[i] <= 30.0 ) ) || ( fabs( chi[i] ) >= 150.0 ) )
						rot[i] = 2;
					if ( ( ( chi[i] >= -90.0 ) && ( chi[i] <= -30.0 ) ) || ( ( chi[i] >= 90.0 ) && ( chi[i] <= 150.0 ) ) )
						rot[i] = 3;
				}

				else if ( ( aan == PHE ) || ( aan == TYR ) )
				{
					if ( ( ( chi[i] >= 30.0 ) && ( chi[i] <= 150.0 ) ) || ( ( chi[i] <= -30.0 ) && ( chi[i] >= -150.0 ) ) )
						rot[i] = 1;
					if ( ( ( chi[i] >= -30.0 ) && ( chi[i] <= 30.0 ) ) || ( fabs( chi[i] ) >= 150.0 ) )
						rot[i] = 2;
				}

				else if ( aan == TRP )
				{
					if ( ( chi[i] >= -180.0 ) && ( chi[i] <= -60.0 ) )
						rot[i] = 1;
					if ( ( chi[i] >= -60.0 ) && ( chi[i] <= 60.0 ) )
						rot[i] = 2;
					if ( ( chi[i] >= 60.0 ) && ( chi[i] <= 180.0 ) )
						rot[i] = 3;
				}
				else if ( aan == ASN )
				{ // chi2 of asn
					//// ctsa - note this is a special case of the new dunbrack rotamer set
					////
					if ( rot[0] == 1 )
					{
						if ( chi[i] >= -150.0 && chi[i] < -90.0 )
							rot[i] = 1;
						if ( chi[i] >= -90.0 && chi[i] < -30.0 )
							rot[i] = 2;
						if ( chi[i] >= -30.0 && chi[i] < 30.0 )
							rot[i] = 3;
						if ( chi[i] >= 30.0 && chi[i] < 90.0 )
							rot[i] = 4;
						if ( chi[i] >= 90.0 && chi[i] < 150.0 )
							rot[i] = 5;
						if ( fabs( chi[i] ) >= 150.0 )
							rot[i] = 6;
					}
					else if ( rot[0] == 2 )
					{
						if ( chi[i] >= -180.0 && chi[i] < -90.0 )
							rot[i] = 1;
						if ( chi[i] >= -90.0 && chi[i] < -45.0 )
							rot[i] = 2;
						if ( chi[i] >= -45.0 && chi[i] < 0.0 )
							rot[i] = 3;
						if ( chi[i] >= 0.0 && chi[i] < 45.0 )
							rot[i] = 4;
						if ( chi[i] >= 45.0 && chi[i] < 90.0 )
							rot[i] = 5;

						if ( chi[i] >= 90.0 && chi[i] <= 180.0 )
							rot[i] = 6;
					}
					else
					{
						if ( chi[i] >= -180.0 && chi[i] < -105.0 )
							rot[i] = 1;
						if ( chi[i] >= -105.0 && chi[i] < -45.0 )
							rot[i] = 2;
						if ( chi[i] >= -45.0 && chi[i] < 15.0 )
							rot[i] = 3;
						if ( chi[i] >= 15.0 && chi[i] < 60.0 )
							rot[i] = 4;
						if ( chi[i] >= 60.0 && chi[i] < 120.0 )
							rot[i] = 5;
						if ( chi[i] >= 120.0 && chi[i] <= 180.0 )
							rot[i] = 6;
					}
				}
			}
			else if ( i == 2 )
			{ // chi 3
				if ( aan == GLU )
				{
					if ( ( ( chi[i] >= 30.0 ) && ( chi[i] <= 90.0 ) ) || ( ( chi[i] <= -90.0 ) && ( chi[i] >= -150.0 ) ) )
						rot[i] = 1;
					if ( ( ( chi[i] >= -30.0 ) && ( chi[i] <= 30.0 ) ) || ( fabs( chi[i] ) >= 150.0 ) )
						rot[i] = 2;
					if ( ( ( chi[i] >= -90.0 ) && ( chi[i] <= -30.0 ) ) || ( ( chi[i] >= 90.0 ) && ( chi[i] <= 150.0 ) ) )
						rot[i] = 3;
				}
				else if ( aan == ARG || aan == LYS || aan == MET )
				{
					if ( ( chi[i] >= 0.0 ) && ( chi[i] <= 120.0 ) )
						rot[i] = 1;
					if ( fabs( chi[i] ) > 120.0 )
						rot[i] = 2;
					if ( ( chi[i] >= -120.0 ) && ( chi[i] <= 0.0 ) )
						rot[i] = 3;
				}
				else if ( aan == GLN )
				{ // chi3 of gln
					//// ctsa - note this is a special case of the new dunbrack rotamer set
					////
					if ( rot[1] == 2 )
					{
						if ( chi[i] >= 135.0 || chi[i] < -135.0 )
							rot[i] = 1;
						if ( chi[i] >= -135.0 && chi[i] < -45.0 )
							rot[i] = 2;
						if ( chi[i] >= -45.0 && chi[i] < 45.0 )
							rot[i] = 3;
						if ( chi[i] >= 45.0 && chi[i] < 135.0 )
							rot[i] = 4;
					}
					else
					{
						if ( chi[i] >= -180.0 && chi[i] < -90.0 )
							rot[i] = 1;
						if ( chi[i] >= -90.0 && chi[i] < 0.0 )
							rot[i] = 2;
						if ( chi[i] >= 0.0 && chi[i] < 90.0 )
							rot[i] = 3;
						if ( chi[i] >= 90.0 && chi[i] <= 180.0 )
							rot[i] = 4;
					}
				}
			}
			else if ( i == 3 )
			{ // chi 4
				if ( aan == ARG || aan == LYS )
				{
					if ( ( chi[i] >= 0.0 ) && ( chi[i] <= 120.0 ) )
						rot[i] = 1;
					if ( fabs( chi[i] ) > 120.0 )
						rot[i] = 2;
					if ( ( chi[i] >= -120.0 ) && ( chi[i] <= 0.0 ) )
						rot[i] = 3;
				}
			}
		}

	}


}


////Calculate all dihedral angles
void Macromolecule::all_dihedrals( float * * dihedral )
{
	Residue * res;
	int resn, resnp=0, nres;
	Segment * seg;
	pdbIter * iter1, * iter2;
	Tcoor atN, atNn, atC, atCp, atCA, atCAn;
	double phi, psi, omega;
	int chino, i;
	double chis[4];
	pdbIter * iter3;
	//Iterador
	 iter1 = new pdbIter( this, false, false, true, false ); // only prot...


	nres = ( int )iter1->num_fragment();

	float * dang;
	dang = ( float * ) malloc( sizeof( float ) * nres * ( 3 + 4 ) );
	int dangi = 0;



	//Bucle para recorrer segmentos
	i = 1;
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{

		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );

		if (  seg->getMolType() == tmol_protein )  // if not hetero
		{

			//
			// get first residue
			//
			iter2->pos_fragment = 0;
			res = ( Residue * ) iter2->get_fragment();
			resn = resnum_from_resname( res->getName() );
			// get atoms
			iter3 = new pdbIter( res );
			iter3->pos_atom = 0;
			// N
			( iter3->get_atom() )->getPosition( atN );
			// CA
			iter3->next_atom();
			( iter3->get_atom() )->getPosition( atCA );
			// C
			iter3->next_atom();
			( iter3->get_atom() )->getPosition( atC );
			iter3->~pdbIter();
			// end get atoms

			// compute chis
			for ( chino = 0; chino < 4; chino++ )
				chis[chino]=9999;
			res->get_achis( resn, chis );
			resnp = resn;

			if (iter2->num_fragment()==1) {  // only one residue case
				dang[dangi + 0] = 9999;
				dang[dangi + 1] = 9999;
				dang[dangi + 2] = 9999;
				res->get_achis( resn, chis );
				for ( chino = 0; chino < AA[resnp].nchi; chino++ )
					dang[dangi + 3 + chino] = chis[chino];
				dangi += 7;

			} else {
				//
				// get second residue
				//
				iter2->next_fragment();
				res = ( Residue * ) iter2->get_fragment();
				resn = resnum_from_resname( res->getName() );
				// get atoms
				iter3 = new pdbIter( res );
				iter3->pos_atom = 0;
				// N_next
				( iter3->get_atom() )->getPosition( atNn );
				// CA_next
				iter3->next_atom();
				( iter3->get_atom() )->getPosition( atCAn );

				// compute psi first   N, CA, C, N_next
				// phi    Cn-1, N, CA, C
				psi = dihedral_bk( atN, atCA, atC, atNn );
				phi = 9999.0;
				omega = dihedral_bk( atCA, atC, atNn, atCAn );

				dang[dangi + 0] = phi;
				dang[dangi + 1] = psi;
				dang[dangi + 2] = omega;
				for ( chino = 0; chino < AA[resnp].nchi; chino++ )
					dang[dangi + 3 + chino] = chis[chino];
				dangi += 7;


				// N
				for ( int i = 0; i < 3; i++ )
					atN[i] = atNn[i];
				// CA
				for ( int i = 0; i < 3; i++ )
					atCA[i] = atCAn[i];
				// C previous
				for ( int i = 0; i < 3; i++ )
					atCp[i] = atC[i];
				// C
				iter3->next_atom();
				( iter3->get_atom() )->getPosition( atC );
				iter3->~pdbIter();
				// end get atoms

				// compute chis
				res->get_achis( resn, chis );
				resnp = resn;


				for ( iter2->pos_fragment=2; !iter2->gend_fragment(); iter2->next_fragment() )
				{


					phi = dihedral_bk( atCp, atN, atCA, atC );
					//
					// get current residue
					//
					res = ( Residue * ) iter2->get_fragment();
					resn = resnum_from_resname( res->getName() );

					// get atoms
					iter3 = new pdbIter( res );
					iter3->pos_atom = 0;

					// N_next
					( iter3->get_atom() )->getPosition( atNn );
					// CA_next
					iter3->next_atom();
					( iter3->get_atom() )->getPosition( atCAn );

					// printf ("hola pos_fragment =%d %d\n",iter2->pos_fragment,iter2->num_fragment());

					psi = dihedral_bk( atN, atCA, atC, atNn );
					omega = dihedral_bk( atCA, atC, atNn, atCAn );

					// compute psi first   N, CA, C, N_next
					// N
					for ( int i = 0; i < 3; i++ )
						atN[i] = atNn[i];
					// CA
					for ( int i = 0; i < 3; i++ )
						atCA[i] = atCAn[i];
					// C previous
					for ( int i = 0; i < 3; i++ )
						atCp[i] = atC[i];
					// C
					iter3->next_atom();
					( iter3->get_atom() )->getPosition( atC );
					iter3->~pdbIter();
					// end get atoms

					dang[dangi + 0] = phi;
					dang[dangi + 1] = psi;
					dang[dangi + 2] = omega;
					for ( chino = 0; chino < AA[resnp].nchi; chino++ )
						dang[dangi + 3 + chino] = chis[chino];
					dangi += 7;

					// compute chis
					res->get_achis( resn, chis );
					resnp = resn;

				}


				// last residue
				dang[dangi + 0] = dihedral_bk( atCp, atN, atCA, atC );;
				dang[dangi + 1] = 9999.0;
				dang[dangi + 2] = 9999.0;

				for ( chino = 0; chino < AA[resnp].nchi; chino++ )
					dang[dangi + 3 + chino] = chis[chino];
				dangi += 7;

				iter2->~pdbIter();

			}
		} // end not t_smol
	}



	iter1->~pdbIter();

	* dihedral = dang;

}


///Store all dihedral angles calculated in a file with the format
//AA PHI PSI OMEGA ROT_INDEX
void Macromolecule::show_all_dihedrals( float * dang, char * filename )
{

	Residue * res;
	int resn, resnp, nres;
	Segment * seg;
	pdbIter * iter1, * iter2;
	int chino, i;
	float * chis;
	chis = new float[4];
	int r1234[4];
	int rot_index;
	char buffer[60];
	Tcoor  atCp;

	ofstream fout( filename );
	if ( !fout.is_open() )
	{
		fprintf(stdout,"Could not output dihedrals to: %s\n",filename);
		return;
	}


	//Iterador
	iter1 = new pdbIter( this, false, false );
	nres = ( int )iter1->num_fragment();


	int dangi = 0;


	//fout << "\n------------------------------------------------------------------------------------" << endl;
	//fout << "             PHI      PSI     OMEGA     ROT_INDEX" << endl;
	//fout << "------------------------------------------------------------------------------------" << endl;



	//Bucle para recorrer segmentos
	i = 1;
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );


		// Si hay mas de un segmento entonces guardar el C del segmento anterior para
		// poder calcular la rotacion y la translacion al primer aminoacido del nuevo
		// segmento
		if (iter1->pos_segment > 0) {
			pdbIter *iter3 = new pdbIter( res );
			iter3->pos_atom = 3;
			// C prev
			( iter3->get_atom() )->getPosition(atCp);
			iter3->~pdbIter();
		}


		//
		// get first residue
		//
		iter2->pos_fragment = 0;
		res = ( Residue * ) iter2->get_fragment();
		resn = resnum_from_resname( res->getName() );

		// Guardar las coordenadas del primer N del nuevo segmento para poder calcular
		// la rotacion y la translacion a el.
		if (iter1->pos_segment > 0) {

			sprintf(buffer, "#\n");
			fout << buffer;

		}



		sprintf( buffer, "%3s ", res->getName() );
		fout << buffer;
		resnp = resn;
		//
		// get second residue
		//
		iter2->next_fragment();
		res = ( Residue * ) iter2->get_fragment();
		//resnp = resn;
		resn = resnum_from_resname( res->getName() );

		for ( chino = 0; chino < 3; chino++ )
		{
			sprintf( buffer, "%8.3f ", dang[dangi + chino] );
			fout << buffer;
		}

		for ( chino = 0; chino < 4; chino++ )
		{
			//sprintf( buffer, "%8.3f ", dang[dangi + 3 + chino] );
			//fout << buffer;
			chis[chino] = dang[dangi + 3 + chino];
		}
		rotamer_from_chi( chis, resnp, r1234 );
		//fprintf( stderr, "-> %1d %1d %1d %1d\n", r1234[0], r1234[1], r1234[2], r1234[3] );

		rot_index = r1234[0] * 1000 + r1234[1] * 100 + r1234[2] * 10 + r1234[3];

		sprintf( buffer, "%4d\n", rot_index );
		fout << buffer;

		//fprintf( stderr, "\n%5d %3s ", i++, res->getName() );
		sprintf( buffer, "%3s ", res->getName() );
		fout << buffer;
		resnp = resn;

		for ( iter2->pos_fragment = 2; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			//
			// get nd residue
			//
			res = ( Residue * ) iter2->get_fragment();
			resn = resnum_from_resname( res->getName() );
			dangi += 7;

			for ( chino = 0; chino < 3; chino++ )
			{
				sprintf( buffer, "%8.3f ", dang[dangi + chino] );
				fout << buffer;
			}

			for ( chino = 0; chino < 4; chino++ )
			{
				//sprintf( buffer, "%8.3f ", dang[dangi + 3 + chino] );
				//fout << buffer;
				chis[chino] = dang[dangi + 3 + chino];
			}
			rotamer_from_chi( chis, resnp, r1234 );

			//fprintf( stderr, "-> %1d %1d %1d %1d", r1234[0], r1234[1], r1234[2], r1234[3] );
			rot_index = r1234[0] * 1000 + r1234[1] * 100 + r1234[2] * 10 + r1234[3];

			sprintf( buffer, "%4d\n", rot_index );
			fout << buffer;

			//fprintf( stderr, "\n%5d %3s ", i++, res->getName() );
			sprintf( buffer, "%3s ", res->getName() );
			fout << buffer;
			resnp = resn;

		}
		iter2->~pdbIter();

		//
		// get last residue
		//
		dangi += 7;

		for ( chino = 0; chino < 3; chino++ )
		{
			sprintf( buffer, "%8.3f ", dang[dangi + chino] );
			fout << buffer;
		}

		for ( chino = 0; chino < 4; chino++ )
		{
			// sprintf( buffer, "%8.3f ", dang[dangi + 3 + chino] );
			// fout << buffer;
			chis[chino] = dang[dangi + 3 + chino];
		}
		rotamer_from_chi( chis, resnp, r1234 );


		//fprintf( stderr, "-> %1d %1d %1d %1d\n", r1234[0], r1234[1], r1234[2], r1234[3] );
		if ( iter1->pos_segment == iter1->num_segment() - 1 ) {
			rot_index = r1234[0] * 1000 + r1234[1] * 100 + r1234[2] * 10 + r1234[3];
			sprintf( buffer, "%4d\n", rot_index );
		}
		else {
			rot_index = r1234[0] * 1000 + r1234[1] * 100 + r1234[2] * 10 + r1234[3];
			sprintf( buffer, "%4d ", rot_index );
		}

		fout << buffer;
		dangi += 7;

	}
	iter1->~pdbIter();

	//fprintf( stderr, "------------------------------------------------------------------------------------\n" );
	// fout << "------------------------------------------------------------------------------------" << endl;

	fout.close();
}

void Macromolecule::CalculateFirstNCACofSegments( vector < vector<float> > & NCACsPDB) {

	pdbIter * iter1, * iter2;
	Residue * res;
	Segment * seg;
	Tcoor at;

	//Iterador
	iter1 = new pdbIter( this,false,false );

	//Bucle para recorrer segmentos
	for ( iter1->pos_segment = 1; !iter1->gend_segment(); iter1->next_segment() ) {

		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );

		// get first residue
		iter2->pos_fragment = 0;
		res = ( Residue * ) iter2->get_fragment();

		// Guardar las coordenadas del primer N del nuevo segmento
		pdbIter *iter3 = new pdbIter( res );

		// Guardar el N, CA y C
		for (iter3->pos_atom = 0; iter3->pos_atom < 3; iter3->pos_atom++) {
			( iter3->get_atom() )->getPosition(at);
			vector < float > a( 3 );
			a[0] = at[0];
			a[1] = at[1];
			a[2] = at[2];
			NCACsPDB.push_back( a );
		}

		iter3->~pdbIter();
		iter2->~pdbIter();


	}

	iter1->~pdbIter();

	//fprintf( stderr, "------------------------------------------------------------------------------------\n" );
	// fout << "------------------------------------------------------------------------------------" << endl;
}

// Mon added (25/4/2007)
// Some refinement of "::show_alldihedrals"...
//
///Store all dihedral angles calculated in a file with the format
//AA PHI PSI OMEGA ROT_INDEX
void Macromolecule::show_all_dihedrals_mon( float * dang, char * filename )
{
	Residue * res;
	int resn;
	Segment * seg;
	pdbIter * iter1, * iter2;
	int chino;
	float * chis;
	chis = new float[4];
	//int *r1234;
	//r1234 = ( int * ) malloc( 4 * sizeof( int ) );
	int r1234[4];
	int rot_index;
	char buffer[60];

	ofstream fout( filename );
	if ( !fout.is_open() )
	{
		fprintf(stdout,"Could not output dihedrals to: %s\n",filename);
		return;
	}

	//Iterador
	iter1 = new pdbIter( this,false,false );
	int dangi = 0;

	//fout << "\n------------------------------------------------------------------------------------" << endl;
	//fout << "             PHI      PSI     OMEGA     ROT_INDEX" << endl;
	//fout << "------------------------------------------------------------------------------------" << endl;



	//Bucle para recorrer segmentos
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{

		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );

		// Main loop to iterate fragments inside a segment
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			// get nd residue
			res = ( Residue * ) iter2->get_fragment();
			resn = resnum_from_resname( res->getName() );

			sprintf( buffer, "%3s ", res->getName() );
			fout << buffer;

			for ( chino = 0; chino < 3; chino++ )
			{
				sprintf( buffer, "%8.3f ", dang[dangi + chino] );
				fout << buffer;
			}
			for ( chino = 0; chino < 4; chino++ ) {
				//        sprintf( buffer, "%8.3f ", dang[dangi + 3 + chino] );
				//        fout << buffer;
				chis[chino] = dang[dangi + 3 + chino];
			}
			rotamer_from_chi( chis, resn, r1234 ); // rotamer_from_chi( chis, resnp, r1234 );
			//fprintf( stderr, "-> %1d %1d %1d %1d", r1234[0], r1234[1], r1234[2], r1234[3] );
			rot_index = r1234[0] * 1000 + r1234[1] * 100 + r1234[2] * 10 + r1234[3];
			sprintf( buffer, "%4d\n", rot_index );
			fout << buffer;

			dangi += 7;
		}
		iter2->~pdbIter();
	}
	iter1->~pdbIter();

	//fprintf( stderr, "------------------------------------------------------------------------------------\n" );
	// fout << "------------------------------------------------------------------------------------" << endl;

	fout.close();
}
