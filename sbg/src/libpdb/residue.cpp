#include <stdio.h>

#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "Macromolecule.h"
#include "Residue.h"
//#include "ResIni.h"
#include <math.h>
#include <vector>

using namespace std;

///Auxiliar function
float periodic_range( float a, float x )
{
  float const halfx = 0.5f * x;
  return ( ( a >= halfx || a < -halfx ) ? fmod( fmod( a, x ) + ( x + halfx ), x ) - halfx : a );
}


///Auxiliar function
float set_chi_to_periodic_range( double chi, int aa, int chino )
{


  //------------------------------------------------------------------------------
  if ( chino == 2 )
  { // handle symmetry cases for chi2
    //  -- not for HIS in new dunbrack version
    if ( aa == 4 || aa == 19 )
    {
      return periodic_range( chi - 60., 180. ) + 60.;
    }
    else if ( aa == 2 )
    {
      return periodic_range( chi, 180. );
    }
  }
  else if ( chino == 3 )
  { // handle symmetry cases for chi3
    if ( aa == 3 )
    {
      return periodic_range( chi, 180. );
    }
  }
  return periodic_range( chi, 360. );
}


///Auxiliar function
float dihedral( float * a1, float * a2, float * a3, float * a4 )
{
  float xij, yij, zij, xkj, ykj, zkj, xkl, ykl, zkl, dxi, dyi, dzi, gxi, gyi, gzi, bi, bk, ct,z1, z2, ap, s;
  //float boi2, boj2, bioj, bjoi;

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

 // boi2 = 1. / bi;
 // boj2 = 1. / bk;
  bi = ( float )sqrt( ( double )bi );
  bk = ( float )sqrt( ( double )bk );

  z1 = 1. / bi;
  z2 = 1. / bk;
 // bioj = bi * z2;
  // bjoi = bk * z1;
  ct = ct * z1 * z2;
  if ( ct > 1.0 ) ct = 1.0;
  if ( ct < ( -1.0 ) ) ct = -1.0;
  ap = acos( ct );
  //fprintf( stderr, "%f %f %f\n", ap,z1,z2);

  s = xkj * ( dzi * gyi - dyi * gzi ) + ykj * ( dxi * gzi - dzi * gxi ) + zkj * ( dyi * gxi - dxi * gyi );

  if ( s < 0.0 ) ap = -ap;

  ap = ( ap > 0.0 ) ? PDB_PI - ap : -( PDB_PI + ap );

  //fprintf( stderr, "d %f %f\n", ap, ap*180/PDB_PI);
  //getchar();

  return ( ap * 180.0 / PDB_PI );
}

///Auxiliar function
float dihedral( double * a1, double * a2, double * a3, double * a4 )
{
  double xij, yij, zij, xkj, ykj, zkj, xkl, ykl, zkl, dxi, dyi, dzi, gxi, gyi, gzi, bi, bk, ct,
       z1, z2, ap, s;
 // double  boi2, boj2, bioj, bjoi;

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

//  boi2 = 1. / bi;
//  boj2 = 1. / bk;
  bi = sqrt( bi );
  bk = sqrt( bk );

  z1 = 1. / bi;
  z2 = 1. / bk;
//  bioj = bi * z2;
//  bjoi = bk * z1;
  ct = ct * z1 * z2;
  if ( ct > 1.0 ) ct = 1.0;
  if ( ct < ( -1.0 ) ) ct = -1.0;
  ap = acos( ct );
  //fprintf( stderr, "%f %f %f\n", ap,z1,z2);

  s = xkj * ( dzi * gyi - dyi * gzi ) + ykj * ( dxi * gzi - dzi * gxi ) + zkj * ( dyi * gxi - dxi * gyi );

  if ( s < 0.0 ) ap = -ap;

  ap = ( ap > 0.0 ) ? PDB_PI - ap : -( PDB_PI + ap );

  //fprintf( stderr, "d %f %f\n", ap, ap*180/PDB_PI);
  //getchar();

  return ( ap * 180.0 / PDB_PI );
}


inline float sin_cos_range( float const x, float const tol =.001f )
{
  assert( ( tol >= 0.0f ) );
  if ( ( x >= -1.0f ) && ( x <= 1.0f ) )
  { // In valid [-1,+1] range
    return x;
  }
  else if ( ( x < -1.0f ) && ( x >= -( 1.0f + tol ) ) )
  { // Within tolerance: Adjust
    return -1.0f;
  }
  else if ( ( x > 1.0f ) && ( x <= 1.0f + tol ) )
  { // Within tolerance: Adjust
    return 1.0f;
  }
  else
  { // Out of range
	  fprintf(stdout,"sin_cos_range ERROR: %f is outside of [-1,+1] range\n",x);
    exit( EXIT_FAILURE );
  }
}

Residue::Residue( Tname name, int in_nid, int in_npos, char i_letter):Fragment(name, in_nid, in_npos, i_letter)
	  {
	  };

Residue::Residue( Residue * old, bool with_elements ):Fragment( old, with_elements )
	  {
	  };


// build from template
Residue::Residue( int aa,  int in_nid, int in_npos,char in_letter,  TMOL t)
{
  Atom * at;
  // strcpy( id, ( char * ) AA[aa].aa_name3 );
  memcpy(id,AA[aa].aa_name3,4);

  limit = AA[aa].natoms;
  nid=in_nid;
  npos=in_npos;
  letter=in_letter;
  fragid = resnum_from_resname(id); // Mon added...

  elements = ( PDB_Contained * * ) malloc( sizeof( PDB_Contained * ) * limit );
  if ( elements == NULL )
	  fprintf(stdout,"Memory allocation error\n");;

  for ( int i = 0; i < limit; i++ )
  {
    at = new Atom( atom_types[AA[aa].atom[i].fullatom_type-1], AA[aa].atom[i].atom_name, AA[aa].atom[i].icoor, i );
    elements[i] = at;
  }

  if ( limit != 0 )
    currentE = 0;
  else
    currentE = -1;
}



void Residue::get_achis( int aan, double * chis )
{
  int chino;
  Tcoor  p1,  p2,  p3,  p4;
  int MAXCHI = 4;

  for ( int i = 0; i < MAXCHI; i++ )
    chis[i] = 0.0;


  for ( chino = 0; chino < AA[aan].nchi; chino++ )
  {
    if ( AA[aan].chi_atoms[chino] [0] >= limit )
    {
      chis[chino] = 9999.0; continue;
    }
    ( ( Atom * ) elements[(int) AA[aan].chi_atoms[chino] [0]] )->getPosition(p1);
    if ( AA[aan].chi_atoms[chino] [1] >= limit )
    {
      chis[chino] = 9999.0; continue;
    }
    ( ( Atom * ) elements[(int)AA[aan].chi_atoms[chino] [1]] )->getPosition(p2);
    if ( AA[aan].chi_atoms[chino] [2] >= limit )
    {
      chis[chino] = 9999.0; continue;
    }
    ( ( Atom * ) elements[(int)AA[aan].chi_atoms[chino] [2]] )->getPosition(p3);
    if ( AA[aan].chi_atoms[chino] [3] >= limit )
    {
      chis[chino] = 9999.0; continue;
    }
    ( ( Atom * ) elements[(int)AA[aan].chi_atoms[chino] [3]] )->getPosition(p4);
    chis[chino] = set_chi_to_periodic_range( dihedral( p1, p2, p3, p4 ), aan, chino );


  }

}

void Residue::rotamerize2( float * newchi )
{
  int chino, i;
  Tcoor c1, c2, c3, c4, p1;
  float cpos[3];
  float ichi, rot;
  float * * mat;
  float * vec;
  int aan;

  //cout << "ROTAMERIZE2 " << this->getName() << " " << newchi[0] << " " << newchi[1] << " " << newchi[2] << " " << newchi[3] << endl;
  //getchar();

  mat = new float * [3];
  for ( int i = 0; i < 3; i++ )
    mat[i] = new float[3];

  vec = new float[3];

  aan = resnum_from_resname( id );

  for ( chino = 0; chino < AA[aan].nchi; chino++ )
  {

    ( ( Atom * ) elements[(int)AA[aan].chi_atoms[chino] [0]] )->getPosition(c1);
    ( ( Atom * ) elements[(int)AA[aan].chi_atoms[chino] [1]] )->getPosition(c2);
    ( ( Atom * ) elements[(int)AA[aan].chi_atoms[chino] [2]] )->getPosition(c3);
    ( ( Atom * ) elements[(int)AA[aan].chi_atoms[chino] [3]] )->getPosition(c4);


    // detemine chi angle
    ichi = dihedral( c1, c2, c3, c4 );

    //** rotate angle by new-ial degrees
    rot = ( newchi[chino] - ichi ) * ( PDB_PI / 180.0 );


    //** generate rotation vector and matrix
    getrotS_bk( c2, c3, rot, mat, vec );


    for ( i = 0; i < AA[aan].natoms; i++ )
    {

      ( ( Atom * ) elements[i] )->getPosition(p1);

      if ( AA[aan].chi_required[chino] [i] )
      {

        // rot a1 coord c output
        cpos[0] = mat[0] [0] * p1[0] + mat[0] [1] * p1[1] + mat[0] [2] * p1[2];
        cpos[1] = mat[1] [0] * p1[0] + mat[1] [1] * p1[1] + mat[1] [2] * p1[2];
        cpos[2] = mat[2] [0] * p1[0] + mat[2] [1] * p1[1] + mat[2] [2] * p1[2];

        // + move
        cpos[0] += vec[0];
        cpos[1] += vec[1];
        cpos[2] += vec[2];

        ( ( Atom * ) elements[i] )->setPosition( cpos );
      }

    }

  }


  /*cout << "COMPROBACION " << endl;
  for ( chino = 0; chino < AA[aan].nchi; chino++ )
 {

   ( ( Atom * ) elements[AA[aan].chi_atoms[chino] [0]] )->getPosition(c1);
   ( ( Atom * ) elements[AA[aan].chi_atoms[chino] [1]] )->getPosition(c2);
   ( ( Atom * ) elements[AA[aan].chi_atoms[chino] [2]] )->getPosition(c3);
   ( ( Atom * ) elements[AA[aan].chi_atoms[chino] [3]] )->getPosition(c4);


   // detemine chi angle
   ichi = dihedral( c1, c2, c3, c4 );

  cout << "chi[" << chino << "] = " << ichi << endl;
 }*/

}


// inline ??
char reschar_from_resname( char * resname )
{

  int k;

  for ( k = 0; k < 19; k++ )
  {
    if ( strcmp( resname, AA[k].aa_name3 ) == 0 )
      return AA[k].aa_name1;
  }

  if ( k >= 20 ) {
    fprintf( stderr, "  Warning %s residue not identified\n", resname );
    exit(-1);
  }
  //
  return  0;

}



//PARA ALINEAR AMINOACIDOS

inline void subvec( float * v1, float * v2, float * v3 )
{
  v3[0] = v1[0] - v2[0];
  v3[1] = v1[1] - v2[1];
  v3[2] = v1[2] - v2[2];
}


inline void unitvec( float * v1, float * v2 )
{

  //Objexx: Casts to double needed to match Fortran precision
  double cc = sqrt( ( static_cast < double > ( v1[0] ) * v1[0] ) + ( static_cast < double > ( v1[1] ) * v1[1] )
       + ( static_cast < double > ( v1[2] ) * v1[2] ) );

  if ( cc == 0.0 ) cc = 0.0001;
  // cems changed division to inverse multiply
  double const invX = 1.0 / cc;
  v2[0] = v1[0] * invX;
  v2[1] = v1[1] * invX;
  v2[2] = v1[2] * invX;
}


inline void cros( float * v1, float * v2, float * v3 )
{

  //Objexx: Casts to double needed to match Fortran precision
  v3[0] = static_cast < double > ( v1[1] ) * v2[2] - static_cast < double > ( v1[2] ) * v2[1];
  v3[1] = static_cast < double > ( v1[2] ) * v2[0] - static_cast < double > ( v1[0] ) * v2[2];
  v3[2] = static_cast < double > ( v1[0] ) * v2[1] - static_cast < double > ( v1[1] ) * v2[0];
}


void Stranspose( float * * mat )
{

  float temp;

  for ( int i = 0; i <= 1; ++i )
  {
    for ( int j = i + 1; j <= 2; ++j )
    {
      temp = mat[i] [j];
      mat[i] [j] = mat[j] [i];
      mat[j] [i] = temp;
    }
  }

}


void angles_coord_sys( float * p1, float * p2, float * p3, float * * mat )
{


  subvec( p1, p2, & mat[0] [0] );

  unitvec( & mat[0] [0], & mat[0] [0] );

  subvec( p3, p2, & mat[1] [0] );

  cros( & mat[0] [0], & mat[1] [0], & mat[2] [0] );

  unitvec( & mat[2] [0], & mat[2] [0] );

  cros( & mat[2] [0], & mat[0] [0], & mat[1] [0] );

  Stranspose( mat );


}


void Smat_multiply( float * * a, float * * b, float * * c_out )
{

  for ( int i = 0; i <= 2; ++i )
  {
    for ( int j = 0; j <= 2; ++j )
    {
      c_out[i] [j] = static_cast < double > ( b[i] [0] ) * a[0] [j] + static_cast < double
           > ( b[i] [1] ) * a[1] [j] + static_cast < double > ( b[i] [2] ) * a[2] [j];
    }
  }
}



inline void Srotate( float * * mat, float * v1, float * v2 )
{

  // it is okay if v1 and v2 are the same array, or use rotate_in_place to make this explicit
  double const a = v1[0]; // v1(1)
  double const b = v1[1]; // v1(2)
  double const c = v1[2]; // v1(3)
  v2[0] = mat[0] [0] * a + mat[0] [1] * b + mat[0] [2] * c; // v2(1)
  v2[1] = mat[1] [0] * a + mat[1] [1] * b + mat[1] [2] * c; // v2(2)
  v2[2] = mat[2] [0] * a + mat[2] [1] * b + mat[2] [2] * c; // v2(3)
}


void Smat_multiply_transpose( float * * a, float * * b, float * * c_out )
{

  // C_out = A^transpose * B
  c_out[0] [0] = static_cast < double > ( b[0] [0] ) * a[0] [0] + static_cast < double > ( b[0] [1] ) * a[0] [1] + static_cast
       < double > ( b[0] [2] ) * a[0] [2];
  c_out[0] [1] = static_cast < double > ( b[0] [0] ) * a[1] [0] + static_cast < double > ( b[0] [1] ) * a[1] [1] + static_cast
       < double > ( b[0] [2] ) * a[1] [2];
  c_out[0] [2] = static_cast < double > ( b[0] [0] ) * a[2] [0] + static_cast < double > ( b[0] [1] ) * a[2] [1] + static_cast
       < double > ( b[0] [2] ) * a[2] [2];
  c_out[1] [0] = static_cast < double > ( b[1] [0] ) * a[0] [0] + static_cast < double > ( b[1] [1] ) * a[0] [1] + static_cast
       < double > ( b[1] [2] ) * a[0] [2];
  c_out[1] [1] = static_cast < double > ( b[1] [0] ) * a[1] [0] + static_cast < double > ( b[1] [1] ) * a[1] [1] + static_cast
       < double > ( b[1] [2] ) * a[1] [2];
  c_out[1] [2] = static_cast < double > ( b[1] [0] ) * a[2] [0] + static_cast < double > ( b[1] [1] ) * a[2] [1] + static_cast
       < double > ( b[1] [2] ) * a[2] [2];
  c_out[2] [0] = static_cast < double > ( b[2] [0] ) * a[0] [0] + static_cast < double > ( b[2] [1] ) * a[0] [1] + static_cast
       < double > ( b[2] [2] ) * a[0] [2];
  c_out[2] [1] = static_cast < double > ( b[2] [0] ) * a[1] [0] + static_cast < double > ( b[2] [1] ) * a[1] [1] + static_cast
       < double > ( b[2] [2] ) * a[1] [2];
  c_out[2] [2] = static_cast < double > ( b[2] [0] ) * a[2] [0] + static_cast < double > ( b[2] [1] ) * a[2] [1] + static_cast
       < double > ( b[2] [2] ) * a[2] [2];
}



void angles_align_transform( float * p1, float * p2, float * p3, float * q1, float * q2, float * q3,
     float * * Mat, float * vector )
     {

       float * * MatQ;
       float * * MatP;
       float * off;

       MatQ = new float * [3];
       for ( int i = 0; i < 3; i++ )
         MatQ[i] = new float[3];

       MatP = new float * [3];
       for ( int i = 0; i < 3; i++ )
         MatP[i] = new float[3];


       off = new float[3];


       angles_coord_sys( p1, p2, p3, MatP );
       angles_coord_sys( q1, q2, q3, MatQ );

       Smat_multiply( MatQ, MatP, Mat );

       //debug_mat("QT *",Mat);
       Stranspose( MatP );

       Smat_multiply( MatQ, MatP, Mat );

       //debug_mat("QPT *",Mat);
       Stranspose( Mat );

       Srotate( Mat, q1, off );
       //for ( int i = 1; i <= 3; ++i ) {
       //        cout << SS( off[i]-p1[i] );
       //} cout << endl;

       Srotate( Mat, q2, off );
       //for ( int i = 1; i <= 3; ++i ) {
       //        cout << SS( off[i]-p2[i] );
       //} cout << endl;

       Srotate( Mat, q3, off );
       //for ( int i = 1; i <= 3; ++i ) {
       //        cout << SS( off[i]-p3[i] );
       //} cout << endl;

       Stranspose( MatQ );

       Smat_multiply( MatQ, MatP, Mat );

       //debug_mat("QTPT *",Mat);
       Stranspose( MatP );

       Smat_multiply( MatQ, MatP, Mat );

       //debug_mat("QTP *",Mat);
       Stranspose( MatQ );

       Smat_multiply_transpose( MatQ, MatP, Mat );

       //debug_mat("QP TP*",Mat);
       Stranspose( Mat );

       Srotate( Mat, q1, off );
       //for ( int i = 1; i <= 3; ++i ) {
       //        cout << SS( off[i]-p1[i] );
       //} cout << endl;

       Srotate( Mat, q2, off );
       //for ( int i = 1; i <= 3; ++i ) {
       //        cout << SS( off[i]-p2[i] );
       //} cout << endl;

       Srotate( Mat, q3, off );
       //for ( int i = 1; i <= 3; ++i ) {
       //        cout << SS( off[i]-p3[i] );
       //} cout << endl;


       vector[0] = q2[0] - ( Mat[0] [0] * p2[0] + Mat[0] [1] * p2[1] + Mat[0] [2] * p2[2] );
       vector[1] = q2[1] - ( Mat[1] [0] * p2[0] + Mat[1] [1] * p2[1] + Mat[1] [2] * p2[2] );
       vector[2] = q2[2] - ( Mat[2] [0] * p2[0] + Mat[2] [1] * p2[1] + Mat[2] [2] * p2[2] );

}

//FIN PARA ALINEAR AMINOACIDOS


//PARA ALINEAR AMINOACIDOS

//==================================================================================================
// Adjust a sine or cosine Value to the Valid [-1,1] Range if Within Specified Tolerance
//==================================================================================================
inline double sin_cos_range( double const x, double const tol =.001 )
{
  assert( ( tol >= 0.0 ) );
  if ( ( x >= -1.0 ) && ( x <= 1.0 ) )
  { // In valid [-1,+1] range
    return x;
  }
  else if ( ( x < -1.0 ) && ( x >= -( 1.0 + tol ) ) )
  { // Within tolerance: Adjust
    return -1.0;
  }
  else if ( ( x > 1.0 ) && ( x <= 1.0 + tol ) )
  { // Within tolerance: Adjust
    return 1.0;
  }
  else
  { // Out of range
	  fprintf(stdout,"sin_cos_range ERROR: %f is outside of [-1,+1] range\n",x);
    exit( EXIT_FAILURE );
  }
}


inline void Cross_bk( float * a, // dimension( 3 )
     float * b, // dimension( 3 )
     float * c, // dimension( 3 )
     float & d )
     {

       //* a X b = c, length of c  = d
       //* cross product A X B
       c[0] = static_cast < double > ( a[1] ) * b[2] - static_cast < double > ( a[2] ) * b[1];
       c[1] = static_cast < double > ( a[2] ) * b[0] - static_cast < double > ( a[0] ) * b[2];
       c[2] = static_cast < double > ( a[0] ) * b[1] - static_cast < double > ( a[1] ) * b[0];
       d = sqrt( c[0] * c[0] + c[1] * c[1] + c[2] * c[2] );
}

inline void getmat_S_bk( float phi, float psi, float kappa, float * * aa )
{

  //***cb Try to avoid underflows...
  if ( fabs( phi ) < 1.0E-4 ) phi = 0.0;
  if ( fabs( psi ) < 1.0E-4 ) psi = 0.0;
  if ( fabs( kappa ) < 1.0E-4 )
  { // Return identity
    aa[0] [1] = aa[0] [2] = aa[1] [0] = aa[1] [2] = aa[2] [0] = aa[2] [1] = 0.0;
    aa[0] [0] = aa[1] [1] = aa[2] [2] = 1.0; // aa(1,1) = aa(2,2) = aa(3,3) = 1.0;
    return;
  }
  float const sf = sin( phi );
  float const cf = cos( phi );
  float const ss = sin( psi );
  float const cs = cos( psi );
  float const sk = sin( -kappa );
  float const ck = cos( -kappa );

  //* now calculate the product matrix of the five rotation matrices
  aa[0] [0] = ( cf * cf + sf * sf * cs * cs ) * ck + sf * sf * ss * ss; // aa(1,1)
  aa[1] [0] = ss * ss * sf * cf * ( ck - 1.0f ) - cs * sk; // aa(1,2)
  aa[2] [0] = sf * cs * ss * ( ck - 1.0f ) + cf * sk * ss; // aa(1,3)
  aa[0] [1] = ss * ss * cf * sf * ( ck - 1.0f ) + cs * sk; // aa(2,1)
  aa[1] [1] = ( sf * sf + cf * cf * cs * cs ) * ck + cf * cf * ss * ss; // aa(2,2)
  aa[2] [1] = cf * ss * cs * ( 1.0f - ck ) + sf * ss * sk; // aa(2,3)
  aa[0] [2] = cs * ss * sf * ( ck - 1.0f ) - cf * ss * sk; // aa(3,1)
  aa[1] [2] = cf * ss * cs * ( 1.0f - ck ) - sf * ss * sk; // aa(3,2)
  aa[2] [2] = ss * ss * ck + cs * cs; // aa(3,3)
  //* done.
}


inline void rotate_S_bk( float * * mat, float * v1, float * v2 )
{
  v2[0] = static_cast < double > ( mat[0] [0] ) * v1[0] + static_cast < double > ( mat[0] [1] ) * v1[1] + static_cast < double
       > ( mat[0] [2] ) * v1[2];
  v2[1] = static_cast < double > ( mat[1] [0] ) * v1[0] + static_cast < double > ( mat[1] [1] ) * v1[1] + static_cast < double
       > ( mat[1] [2] ) * v1[2];
  v2[2] = static_cast < double > ( mat[2] [0] ) * v1[0] + static_cast < double > ( mat[2] [1] ) * v1[1] + static_cast < double
       > ( mat[2] [2] ) * v1[2];
}


void lineup_bk( float * a1, float * a2, float * b1, float * b2, float * * mat, float * vec )
{

  static float a[3];
  static float b[3];
  static float c[3];

  //* get difference vectors
  a[0] = a2[0] - a1[0];
  a[1] = a2[1] - a1[1];
  a[2] = a2[2] - a1[2];
  b[0] = b2[0] - b1[0];
  b[1] = b2[1] - b1[1];
  b[2] = b2[2] - b1[2];


  //* get vector lengths
  //* get a.dot.b
  float const da = sqrt( a[0] * a[0] + a[1] * a[1] + a[2] * a[2] );
  float const db = sqrt( b[0] * b[0] + b[1] * b[1] + b[2] * b[2] );
  float dab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];


  //* get chi
  //cb The following is to protect against underflows on the Mac:
  float const chi = ( ( fabs( dab - ( da * db ) ) >= 0.000001 ) ? acos( sin_cos_range( dab / ( da * db ) ) ) : 0.0 );

  float psi;
  float phi;
  if ( chi >= 1.0E-4 )
  {
    //* cross product A X B = C
    Cross_bk( a, b, c, dab );

    //* get phi, psi
    psi = acos( sin_cos_range( c[2] / dab ) );
    phi = ( ( c[0] != 0.0 || c[1] != 0.0 ) ? atan2( -c[0], c[1] ) : 0.0 );
  }
  else
  {
    psi = 0.0;
    phi = 0.0;
  }

  //* get matrix
  getmat_S_bk( phi, psi, chi, mat );

  //* get vector, = -Mat*a1 + b1
  rotate_S_bk( mat, a1, c );
  vec[0] = -c[0] + b1[0];
  vec[1] = -c[1] + b1[1];
  vec[2] = -c[2] + b1[2];
}



void align_bk( float * a1, float * a2, float * a3, float * b1, float * * mat, float * vec )
{

  float a[3];
  float b[3];
  float c[3];
  float ab[3];
  float ac[3];
  float aba[3];
  float dab, dac, daba, ab_ac, aba_ac, chi;

  //* get difference vectors
  for ( int i = 0; i < 3; ++i )
  {
    a[i] = a1[i] - a2[i];
    b[i] = a3[i] - a2[i];
    c[i] = b1[i] - a2[i];
  }

  Cross_bk( b, a, ab, dab );
  Cross_bk( c, a, ac, dac );
  Cross_bk( ab, a, aba, daba );

  //* get vector lengths
  //* get ab.dot.ac
  ab_ac = ab[0] * ac[0] + ab[1] * ac[1] + ab[2] * ac[2];

  //*    aba.dot.ac
  aba_ac = aba[0] * ac[0] + aba[1] * ac[1] + aba[2] * ac[2];

  //* cosine ab_ac
  ab_ac /= dab;

  //* cosine aba_ac = sine ab_ac
  aba_ac /= daba;

  //* get chi
  chi = atan2( aba_ac, ab_ac );

  //* get matrix for chi rotation
  getrotS_bk( a1, a2, chi, mat, vec );
}

//FIN PARA ALINEAR AMINOACIDOS

TElement Residue::getClass()
{
  return pdb_residue;
}

void Residue::Align( Residue * resi, COORDS positions )
{

  Tcoor pos1_NH, pos2_NH, pos1_CA, pos2_CA, pos1_C, pos2_C, position;
  float * * matrix;
  float cpositions[3];
  float vector[3];
  float *position2;

  matrix = new float * [3];
  for ( int i = 0; i < 3; i++ )
    matrix[i] = new float[3];


  ( ( Atom * ) elements[0] )->getPosition(pos1_NH);
  ( ( Atom * ) resi->getE( 0 ) )->getPosition(pos2_NH);

  ( ( Atom * ) elements[1] )->getPosition(pos1_CA);
  ( ( Atom * ) resi->getE( 1 ) )->getPosition(pos2_CA);

  //angles_align_transform(pos1_NH, pos1_CA, pos1_C, pos2_NH, pos2_CA, pos2_C, matrix, vector);

  lineup_bk( pos1_CA, pos1_NH, pos2_CA, pos2_NH, matrix, vector );

  for ( int i = 0; i < limit; i++ )
  {

    ( ( Atom * ) elements[i] )->getPosition(position);

    //Rotation
    cpositions[0] = matrix[0] [0] * position[0] + matrix[0] [1] * position[1] + matrix[0] [2] * position[2];
    cpositions[1] = matrix[1] [0] * position[0] + matrix[1] [1] * position[1] + matrix[1] [2] * position[2];
    cpositions[2] = matrix[2] [0] * position[0] + matrix[2] [1] * position[1] + matrix[2] [2] * position[2];

    //Trasnlation
    positions[i] [0] = cpositions[0] + vector[0];
    positions[i] [1] = cpositions[1] + vector[1];
    positions[i] [2] = cpositions[2] + vector[2];

  }


  ( ( Atom * ) elements[0] )->getPosition(pos1_NH);
  ( ( Atom * ) elements[1] )->getPosition(pos1_CA);
  ( ( Atom * ) elements[2] )->getPosition(pos1_C);

  ( ( Atom * ) resi->getE( 2 ) )->getPosition(pos2_C);


  //Si alineamos con la funcion angles_align_transform no hace falta la llamada a la siguiente funcion
  //ni el siguiente bucle for
  align_bk( pos1_NH, pos1_CA, pos1_C, pos2_C, matrix, vector );

  for ( int i = 2; i < limit; i++ )
  {
    position2 = positions[i];

    //Rotation
    cpositions[0] = matrix[0] [0] * position2[0] + matrix[0] [1] * position2[1] + matrix[0] [2] * position2[2];
    cpositions[1] = matrix[1] [0] * position2[0] + matrix[1] [1] * position2[1] + matrix[1] [2] * position2[2];
    cpositions[2] = matrix[2] [0] * position2[0] + matrix[2] [1] * position2[1] + matrix[2] [2] * position2[2];

    //Trasnlation
    positions[i] [0] = cpositions[0] + vector[0];
    positions[i] [1] = cpositions[1] + vector[1];
    positions[i] [2] = cpositions[2] + vector[2];
  }


}


float periodic_range2( float a, float x )
{
  float const halfx = 0.5f * x;
  return ( ( a >= halfx || a < -halfx ) ? fmod( fmod( a, x ) + ( x + halfx ), x ) - halfx : a );
}

void Residue::Align_backbone( float NH[3], float CA[3], float C[3], char atom[2] )
{

  Tcoor pos_NH, pos_CA, pos_C, position;
  float * * matrix;
  float * cpositions;
  float vector[3];

  matrix = new float * [3];
  for ( int i = 0; i < 3; i++ )
    matrix[i] = new float[3];

  cpositions = new float[3];

  ( ( Atom * ) elements[0] )->getPosition(pos_NH);
  ( ( Atom * ) elements[1] )->getPosition(pos_CA);
  ( ( Atom * ) elements[2] )->getPosition(pos_C);

  if (!strcmp(atom, "NH"))
    lineup_bk( pos_NH, pos_CA, NH, CA, matrix, vector );
  else if (!strcmp(atom, "CA"))
    lineup_bk( pos_CA, pos_NH, CA, NH, matrix, vector );
  else if (!strcmp(atom, " C"))
    lineup_bk( pos_C, pos_CA, C, CA, matrix, vector );

  for ( int i = 0; i < limit; i++ )
  {
    ( ( Atom * ) elements[i] )->getPosition(position);

    //Rotation
    cpositions[0] = matrix[0] [0] * position[0] + matrix[0] [1] * position[1] + matrix[0] [2] * position[2];
    cpositions[1] = matrix[1] [0] * position[0] + matrix[1] [1] * position[1] + matrix[1] [2] * position[2];
    cpositions[2] = matrix[2] [0] * position[0] + matrix[2] [1] * position[1] + matrix[2] [2] * position[2];

    //Trasnlation
    cpositions[0] += vector[0];
    cpositions[1] += vector[1];
    cpositions[2] += vector[2];

    ( ( Atom * ) elements[i] )->setPosition( cpositions );

  }


  ( ( Atom * ) elements[0] )->getPosition(pos_NH);
  ( ( Atom * ) elements[1] )->getPosition(pos_CA);
  ( ( Atom * ) elements[2] )->getPosition(pos_C);


  if (!strcmp(atom, "NH"))
    align_bk( pos_NH, pos_CA, pos_C, C, matrix, vector );
  else if (!strcmp(atom, "CA"))
    align_bk( pos_CA, pos_NH, pos_C, C, matrix, vector );
  else if (!strcmp(atom, " C"))
    align_bk( pos_C, pos_CA, pos_NH, NH, matrix, vector );


  for ( int i = 0; i < limit; i++ )
  {

    ( ( Atom * ) elements[i] )->getPosition(position);

    //Rotation
    cpositions[0] = matrix[0] [0] * position[0] + matrix[0] [1] * position[1] + matrix[0] [2] * position[2];
    cpositions[1] = matrix[1] [0] * position[0] + matrix[1] [1] * position[1] + matrix[1] [2] * position[2];
    cpositions[2] = matrix[2] [0] * position[0] + matrix[2] [1] * position[1] + matrix[2] [2] * position[2];

    //Trasnlation
    cpositions[0] += vector[0];
    cpositions[1] += vector[1];
    cpositions[2] += vector[2];

    ( ( Atom * ) elements[i] )->setPosition( cpositions );

  }
}

// Created by Mon (25/03/2008)
void Residue::Align_backbone( float NH[3], float CA[3], float C[3], char atom[2], float ***matrix_out, float **vector_out,float ***matrix_out2, float **vector_out2 )
{

  Tcoor pos_NH, pos_CA, pos_C, position;
  float **matrix,**matrix2;
//  float **matrix_tot;
  float *cpositions;
  float vector[3],vector2[3];
//  float *vector_tot;

//  vector_tot = new float [3];

  matrix = new float * [3];
  for ( int i = 0; i < 3; i++ )
    matrix[i] = new float[3];

  matrix2 = new float * [3];
  for ( int i = 0; i < 3; i++ )
    matrix2[i] = new float[3];

/*  matrix_tot = new float * [3];
  for ( int i = 0; i < 3; i++ )
    matrix_tot[i] = new float[3]; */

  cpositions = new float[3];

  ( ( Atom * ) elements[0] )->getPosition(pos_NH);
  ( ( Atom * ) elements[1] )->getPosition(pos_CA);
  ( ( Atom * ) elements[2] )->getPosition(pos_C);

  if (!strcmp(atom, "NH"))
    lineup_bk( pos_NH, pos_CA, NH, CA, matrix, vector );
  else if (!strcmp(atom, "CA"))
    lineup_bk( pos_CA, pos_NH, CA, NH, matrix, vector );
  else if (!strcmp(atom, " C"))
    lineup_bk( pos_C, pos_CA, C, CA, matrix, vector );

printf("Vector1: %f %f %f\n",vector[0],vector[1],vector[2]);

  for ( int i = 0; i < limit; i++ )
  {
    ( ( Atom * ) elements[i] )->getPosition(position);

    //Rotation
    cpositions[0] = matrix[0] [0] * position[0] + matrix[0] [1] * position[1] + matrix[0] [2] * position[2];
    cpositions[1] = matrix[1] [0] * position[0] + matrix[1] [1] * position[1] + matrix[1] [2] * position[2];
    cpositions[2] = matrix[2] [0] * position[0] + matrix[2] [1] * position[1] + matrix[2] [2] * position[2];

    //Trasnlation
    cpositions[0] += vector[0];
    cpositions[1] += vector[1];
    cpositions[2] += vector[2];

    ( ( Atom * ) elements[i] )->setPosition( cpositions );
  }

  ( ( Atom * ) elements[0] )->getPosition(pos_NH);
  ( ( Atom * ) elements[1] )->getPosition(pos_CA);
  ( ( Atom * ) elements[2] )->getPosition(pos_C);

  if (!strcmp(atom, "NH"))
    align_bk( pos_NH, pos_CA, pos_C, C, matrix2, vector2 );
  else if (!strcmp(atom, "CA"))
    align_bk( pos_CA, pos_NH, pos_C, C, matrix2, vector2 );
  else if (!strcmp(atom, " C"))
    align_bk( pos_C, pos_CA, pos_NH, NH, matrix2, vector2 );

printf("Vector2: %f %f %f\n",vector2[0],vector2[1],vector2[2]);
  for ( int i = 0; i < limit; i++ )
  {

    ( ( Atom * ) elements[i] )->getPosition(position);

    //Rotation
    cpositions[0] = matrix2[0] [0] * position[0] + matrix2[0] [1] * position[1] + matrix2[0] [2] * position[2];
    cpositions[1] = matrix2[1] [0] * position[0] + matrix2[1] [1] * position[1] + matrix2[1] [2] * position[2];
    cpositions[2] = matrix2[2] [0] * position[0] + matrix2[2] [1] * position[1] + matrix2[2] [2] * position[2];

    //Trasnlation
    cpositions[0] += vector2[0];
    cpositions[1] += vector2[1];
    cpositions[2] += vector2[2];

    ( ( Atom * ) elements[i] )->setPosition( cpositions );
  }

/*  // Single rotation computation
  for(int i=0;i<3;i++)
  	for(int j=0;j<3;j++)
  		matrix_tot[i][j] =  matrix[i][j] * matrix2[i][j];

  // Single traslation computation
  for(int i=0;i<3;i++)
	vector_tot[i] = vector[i] + vector2[i];
*/

  *matrix_out=matrix; // Outputs Rotation matrix
  *vector_out=vector; // Outputs Traslation vector
  *matrix_out2=matrix2; // Outputs Rotation matrix
  *vector_out2=vector2; // Outputs Traslation vector
}

// Simple Residue Rotation (no traslation)
// Created by Mon (25/03/2008)
void Residue::Rotate( float **matrix )
{

  Tcoor position;
  float *cpositions;
  cpositions = new float[3];

  for ( int i = 0; i < limit; i++ )
  {
    ( ( Atom * ) elements[i] )->getPosition(position);

    //Rotation
    cpositions[0] = matrix[0] [0] * position[0] + matrix[0] [1] * position[1] + matrix[0] [2] * position[2];
    cpositions[1] = matrix[1] [0] * position[0] + matrix[1] [1] * position[1] + matrix[1] [2] * position[2];
    cpositions[2] = matrix[2] [0] * position[0] + matrix[2] [1] * position[1] + matrix[2] [2] * position[2];

    ( ( Atom * ) elements[i] )->setPosition( cpositions );
  }
}

// Simple Residue Rotation and Translation
// Created by Mon (25/03/2008)
void Residue::Rotrans( float **matrix, float *vector )
{

  Tcoor position;
  float *cpositions;
  cpositions = new float[3];

  for ( int i = 0; i < limit; i++ )
  {
    ( ( Atom * ) elements[i] )->getPosition(position);

    //Rotation
    cpositions[0] = matrix[0] [0] * position[0] + matrix[0] [1] * position[1] + matrix[0] [2] * position[2];
    cpositions[1] = matrix[1] [0] * position[0] + matrix[1] [1] * position[1] + matrix[1] [2] * position[2];
    cpositions[2] = matrix[2] [0] * position[0] + matrix[2] [1] * position[1] + matrix[2] [2] * position[2];

    //Trasnlation
    cpositions[0] += vector[0];
    cpositions[1] += vector[1];
    cpositions[2] += vector[2];

    ( ( Atom * ) elements[i] )->setPosition( cpositions );
  }
}

void Residue::OriginCA()
{

  float translation[3]; // Is the position of the CA atom
  Tcoor position, new_position;

  ( ( Atom * ) elements[1] )->getPosition(position);
  for ( int i = 0; i < 3; i++ )
    translation[i] = position[i];

  for ( int i = 0; i < limit; i++ )
  {
    ( ( Atom * ) elements[i] )->getPosition(position);

    new_position[0] = position[0] - translation[0];
    new_position[1] = position[1] - translation[1];
    new_position[2] = position[2] - translation[2];

   ( ( Atom * ) elements[i] )->setPosition(new_position);


  }
}



void Residue::OriginCA( COORDS positions )
{

  Tcoor translation; // Is the position of the CA atom
  Tcoor position;

  ( ( Atom * ) elements[1] )->getPosition(translation);

  for ( int i = 0; i < limit; i++ )
  {
    ( ( Atom * ) elements[i] )->getPosition(position);

    positions[i] [0] = position[0] - translation[0];
    positions[i] [1] = position[1] - translation[1];
    positions[i] [2] = position[2] - translation[2];
  }
}


void Residue::OriginCA2( COORDS positions )
{

  float translation[3]; // Is the position of the CA atom

  for ( int i = 0; i < 3; i++ )
  {
    translation[i] = positions[1] [i];
  }

  for ( int i = 0; i < limit; i++ )
  {
    positions[i] [0] -= translation[0];
    positions[i] [1] -= translation[1];
    positions[i] [2] -= translation[2];
  }
}

void Residue::OriginC()
{

  float translation[3]; // Is the position of the C atom
  Tcoor position, new_position;

  ( ( Atom * ) elements[2] )->getPosition(position);
  for ( int i = 0; i < 3; i++ )
    translation[i] = position[i];

  for ( int i = 0; i < limit; i++ )
  {
    ( ( Atom * ) elements[i] )->getPosition(position);

    new_position[0] = position[0] - translation[0];
    new_position[1] = position[1] - translation[1];
    new_position[2] = position[2] - translation[2];

   ( ( Atom * ) elements[i] )->setPosition(new_position);

  }
}


void Residue::place_atom( int atom )
{

  float * * matrix;
  float * vector;
  float * * template_coords;
  int template_atoms[3];
  int id;
  float cpositions[3];
  Tcoor p1, pos1, pos2;

  id = resnum_from_resname( getName() );

  matrix = new float * [3];
  for ( int i = 0; i < 3; i++ )
    matrix[i] = new float[3];

  vector = new float[3];

  template_coords = new float * [4];
  for ( int i = 0; i < 4; i++ )
    template_coords[i] = new float[3];

  // Almacena en template_atoms el indice de los tres atomos a partir
  // de los cuales se puede reconstruir el atomo
  for ( int i = 0; i < 3; i++ )
  {
    template_atoms[i] = AA[id].atom[atom].ta[i];
  }

  // Lee en las coordenadas del aminoacido template las coordenadas
  // de esos tres atomos
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
    {
      template_coords[i] [j] = AA[id].atom[template_atoms[i]].icoor[j];
    }

  for ( int i = 0; i < 3; i++ )
    template_coords[3] [i] = AA[id].atom[atom].icoor[i];


  // ALinea el template
  (( Atom * ) elements[template_atoms[0]])->getPosition(pos1);
  (( Atom * ) elements[template_atoms[1]])->getPosition(pos2);
  lineup_bk( template_coords[0], template_coords[1], pos1, pos2, matrix, vector );

  for ( int i = 0; i < 4; i++ )
  {

    //Rotation
    cpositions[0] = matrix[0] [0] * template_coords[i] [0] + matrix[0] [1] * template_coords[i] [1] + matrix[0] [2]
         * template_coords[i] [2];
    cpositions[1] = matrix[1] [0] * template_coords[i] [0] + matrix[1] [1] * template_coords[i] [1] + matrix[1] [2]
         * template_coords[i] [2];
    cpositions[2] = matrix[2] [0] * template_coords[i] [0] + matrix[2] [1] * template_coords[i] [1] + matrix[2] [2]
         * template_coords[i] [2];

    //Trasnlation
    cpositions[0] += vector[0];
    cpositions[1] += vector[1];
    cpositions[2] += vector[2];

    for ( int j = 0; j < 3; j++ )
      template_coords[i] [j] = cpositions[j];

  }


  (( Atom * ) elements[template_atoms[2]])->getPosition(pos1);
  align_bk( template_coords[0], template_coords[1], template_coords[2], pos1, matrix, vector );

  for ( int i = 2; i < 4; i++ )
  {

    //Rotation
    cpositions[0] = matrix[0] [0] * template_coords[i] [0] + matrix[0] [1] * template_coords[i] [1] + matrix[0] [2]
         * template_coords[i] [2];
    cpositions[1] = matrix[1] [0] * template_coords[i] [0] + matrix[1] [1] * template_coords[i] [1] + matrix[1] [2]
         * template_coords[i] [2];
    cpositions[2] = matrix[2] [0] * template_coords[i] [0] + matrix[2] [1] * template_coords[i] [1] + matrix[2] [2]
         * template_coords[i] [2];

    //Trasnlation
    cpositions[0] += vector[0];
    cpositions[1] += vector[1];
    cpositions[2] += vector[2];

    for ( int j = 0; j < 3; j++ )
      template_coords[i] [j] = cpositions[j];

  }


  // Almacena la coordenada final
  ( ( Atom * ) elements[atom] )->setPosition( template_coords[3] );
  ( ( Atom * ) elements[atom] )->getPosition(p1);
}

TMOL Residue::getMolType()
{
	return tmol_protein;
}


/*
int Create_Exclusions_Pairs( int aa, int atom, int atom_bonded, int nbonds, vector < int > & atoms )
{
  int Nbonded;
  int neighbor;
  bool bonded_atoms = false;
  bool is;

  ofstream fout;

  fout.open( "EXCLUSIONS", ofstream::app );
  if ( !fout.is_open() )
  {
    printf( "Error opening the file EXCLUSIONS\n" );
    exit( -1 );
  }


  Nbonded = AA[aa].atom[atom_bonded].nbonded_neighbors;


  for ( int i = 0; i < Nbonded; i++ )
  {
    neighbor = AA[aa].atom[atom_bonded].bonded_neighbor[i];
    is = false;

    if ( neighbor != atom )
    {
      for ( int i = 0; i < atoms.size(); i++ )
      {
        if ( atoms[i] == neighbor )
        {
          is = true;
          break;
        }
      }
      if ( !is )
      {
        atoms.push_back( neighbor );
        bonded_atoms = true;
        if ( atom < neighbor )
        {
          fout << " 0 " << AA[aa].aa_name3 << " " << AA[aa].atom[atom].atom_name << " " << AA[aa].aa_name3 << " "
               << AA[aa].atom[neighbor].atom_name << " " << nbonds << endl;
        }
      }

    }
  }

  fout.close();
  return bonded_atoms;

}*/


//find the matirx (mat) and vector (vec) that rotate chi deg about an axis
//define by a1-->a2 rotation is right-handed. That is, it is a clockwise
//rotation when looking from a1 to a2 chi in radians
void getrotS_bk( float * a1, float * a2, float chi, float * * mat, float * vec )
{
  float a[3];
  float c[3];
  int i, j;

  if ( fabs( chi ) < 1.0E-4 )
  {
    for ( i = 0; i < 3; i++ )
      for ( j = 0; j < 3; j++ )
        mat[i] [j] = 0.0;
    for ( i = 0; i < 3; i++ )
      vec[i] = 0.0;
    mat[0] [0] = mat[1] [1] = mat[2] [2] = 1.0;

  }
  else
  {

    //* get difference vector: a
    for ( i = 0; i < 3; i++ )
      a[i] = a2[i] - a1[i];

    //* get length of vector A
    float const da = sqrt( a[0] * a[0] + a[1] * a[1] + a[2] * a[2] );
    //* get phi, psi
    float psi = acos( sin_cos_range( a[2] / da ) );
    float phi = 0.0;
    if ( a[0] != 0.0 || a[1] != 0.0 ) phi = atan2( -a[0], a[1] );
    //* get matrix

    //Try to avoid underflows...
    if ( fabs( phi ) < 1.0E-4 ) phi = 0.0;
    if ( fabs( psi ) < 1.0E-4 ) psi = 0.0;
    if ( fabs( chi ) < 1.0E-4 )
    {

      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          if ( i != j || i == 1 )
            mat[i] [j] = 1.0;

    }


    float const sf = sin( phi );
    float const cf = cos( phi );
    float const ss = sin( psi );
    float const cs = cos( psi );
    float const sk = sin( -chi );
    float const ck = cos( -chi );
    //* now calculate the product matrix of the five rotation matrices
    mat[0] [0] = ( cf * cf + sf * sf * cs * cs ) * ck + sf * sf * ss * ss; // aa(1,1)
    mat[1] [0] = ss * ss * sf * cf * ( ck - 1.0 ) - cs * sk; // aa(1,2)
    mat[2] [0] = sf * cs * ss * ( ck - 1.0 ) + cf * sk * ss; // aa(1,3)
    mat[0] [1] = ss * ss * cf * sf * ( ck - 1.0 ) + cs * sk; // aa(2,1)
    mat[1] [1] = ( sf * sf + cf * cf * cs * cs ) * ck + cf * cf * ss * ss; // aa(2,2)
    mat[2] [1] = cf * ss * cs * ( 1.0 - ck ) + sf * ss * sk; // aa(2,3)
    mat[0] [2] = cs * ss * sf * ( ck - 1.0 ) - cf * ss * sk; // aa(3,1)
    mat[1] [2] = cf * ss * cs * ( 1.0 - ck ) - sf * ss * sk; // aa(3,2)
    mat[2] [2] = ss * ss * ck + cs * cs; // aa(3,3)


    //* get vector, = -Mat*a1 + a1
    c[0] = mat[0] [0] * a1[0] + mat[0] [1] * a1[1] + mat[0] [2] * a1[2];
    c[1] = mat[1] [0] * a1[0] + mat[1] [1] * a1[1] + mat[1] [2] * a1[2];
    c[2] = mat[2] [0] * a1[0] + mat[2] [1] * a1[1] + mat[2] [2] * a1[2];


    for ( i = 0; i < 3; i++ )
      vec[i] = a1[i] - c[i];
  }

}


void angle_bk(float *a1, float *a2, float *b1, float *b2, float & ang){

  float a[3], b[3];
  float da,db,dab, angd;

 // determine the angle between the two vectors a1->a2,b1->b2
    ang = 0.0;

 // * get difference vectors
    for ( int i = 0; i < 3; ++i ) {
       a[i] = a2[i] - a1[i];
       b[i] = b2[i] - b1[i];
    }

  // * get vector lengths
  // * get a.dot.b
    da = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    db = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
    dab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    da = sqrt(da);
    db = sqrt(db);


  // * get angle
  // The following is to protect against underflows on the Mac:
   if ( fabs(dab-(da*db)) < 0.000001 ) {
      angd = 0.0;
   }
   else {
      angd = acos(sin_cos_range(dab/(da*db)));
   }

  ang = (float)angd;

}

void angle_bk(double *a1, double *a2, double *b1, double *b2, double & ang){

  double a[3], b[3];
  double da,db,dab, angd;

 // determine the angle between the two vectors a1->a2,b1->b2
    ang = 0.0;

 // * get difference vectors
    for ( int i = 0; i < 3; ++i ) {
       a[i] = a2[i] - a1[i];
       b[i] = b2[i] - b1[i];
    }

  // * get vector lengths
  // * get a.dot.b
    da = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    db = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
    dab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    da = sqrt(da);
    db = sqrt(db);

  // * get angle
  // The following is to protect against underflows on the Mac:
   if ( fabs(dab-(da*db)) < 0.000001 ) {
      angd = 0.0;
   }
   else {
      angd = acos(sin_cos_range(dab/(da*db)));
   }

  ang = angd;

}


